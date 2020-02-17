"""
**************************************************************************
**  rG4-seeker
**  A pipeline for processing and analyzing rG4-seq data
**
**  Version 1.0 -- August 10, 2019
**
**  Copyright (C) 2019 by Eugene Yui-Ching Chow, Ting-Fung Chan, All rights reserved.
**  Contact:  eugenechow823@gmail.com, tf.chan@cuhk.edu.hk
**  Organization:  School of Life Sciences, The Chinese University of Hong Kong,
**                 Shatin, NT, Hong Kong SAR
**
**  This file is part of rG4-seeker.
**
**  rG4-seeker is free software; you can redistribute it and/or
**  modify it under the terms of the GNU General Public License
**  as published by the Free Software Foundation; either version
**  3 of the License, or (at your option) any later version.
**
**  rG4-seeker is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public
**  License along with rG4-seeker; if not, see
**  <http://www.gnu.org/licenses/>.
**************************************************************************/
"""

import HTSeq
from xopen import xopen
from collections import namedtuple, Counter, OrderedDict
from itertools import combinations
import re
import pyfaidx
import csv
import time
import logging
# rg4seeker imports
from rg4seeker import gff_parser, fsr_precompute

FSRTSV_entry = namedtuple('fsrtsv_entry', ['Interval', 'treatment_coverage', 'treatment_readstart', 'control_coverage',
                                           'control_readstart', 'first_stage_alpha', 'second_stage_alpha_list'])
RTS_csv_header = ['chromosome',
                  'start',
                  'end',
                  'rG4_structural_class',
                  'avg_sequencing_coverage',
                  'strand',
                  'rG4_intervals',
                  'splicing',
                  'sequence_diagram', ]
FIRST_STAGE_ALPHA_DICT = {'0.2': 0, '0.15874': 1, '0.12599': 2, '0.1': 3, '0.07937': 4, '0.062996': 5, '0.05': 6}
COVERAGE_MIN_FILTER = 6
FIRST_STAGE_ALPHA_ALLOWANCE = 0.2
FIRST_STAGE_MIN_FSR_DIFF = 0.05
FIRST_STAGE_ALPHA_VALUES = tuple([0.2, 0.15874, 0.12599, 0.1, 0.07937, 0.062996, 0.05])
SECOND_STAGE_EXPECTED_READSTOP_PROB = 0.025
RE_g6 = re.compile('G{6,}')
RE_digest_iv_str = re.compile('^(.+):([0-9]+)-([0-9]+):([+-.])$')
RTS_3PRIME_ALLOWANCE = 2
FIVEPRIME_BLACKLIST = 5


class FaidxGenome:
    def __init__(self, in_fasta_path):
        self.fasta = pyfaidx.Fasta(in_fasta_path)

    def get_fasta(self, _interval):
        if _interval.strand == '+':
            return str(self.fasta[_interval.chrom][_interval.start:_interval.end].seq).upper()
        elif _interval.strand == '-':
            return str(self.fasta[_interval.chrom][_interval.start:_interval.end].reverse.complement.seq).upper()

    def get_g_percent(self, _interval):
        flanking_fasta = self.get_fasta(_interval)
        best_g_fraction = 0.0
        best_i = -1
        for i in range(0, 40):
            curr_g_fraction = max(flanking_fasta[i:].count('G') / (50 - i), flanking_fasta[i:-1].count('G') / (50 - i))
            if best_g_fraction < curr_g_fraction:
                best_g_fraction = curr_g_fraction
                best_i = i
        return best_g_fraction * 100, flanking_fasta[best_i:]


def parse_fsrtsv(line, reverse_treatment_control=False):
    fields = line.strip().split('\t')
    if len(fields) < 8:
        raise Exception('Incorrect fsr tsv format')
    if len(fields) < 9:
        return None
    assert int(fields[2]) - int(fields[1]) == 1
    _interval = gff_parser.Interval(fields[0], int(fields[1]), int(fields[2]), fields[3])
    if not reverse_treatment_control:
        _treatment_coverage = int(fields[4])
        _treatment_readstart = int(fields[5])
        _control_coverage = int(fields[6])
        _control_readstart = int(fields[7])
    else:
        _treatment_coverage = int(fields[6])
        _treatment_readstart = int(fields[7])
        _control_coverage = int(fields[4])
        _control_readstart = int(fields[5])
    if len(fields) > 8:
        _first_stage_alpha = float(fields[8])
    else:
        _first_stage_alpha = 0.0
    if len(fields) > 9:
        _second_stage_alpha = [float(x) for x in fields[9:]]
    else:
        _second_stage_alpha = None
    return FSRTSV_entry(_interval, _treatment_coverage, _treatment_readstart, _control_coverage, _control_readstart,
                        _first_stage_alpha, _second_stage_alpha)


class FSRMergeBuffer:
    def __init__(self, _entry):
        self.buffer = [_entry]
        self.Interval = gff_parser.Interval(_entry.Interval.chrom, _entry.Interval.start, _entry.Interval.end,
                                            _entry.Interval.strand)

    def is_adjacent_to(self, _new_entry):
        if not (
                self.Interval.strand == _new_entry.Interval.strand and
                self.Interval.chrom == _new_entry.Interval.chrom and
                self.Interval.end == _new_entry.Interval.start
        ):
            return False
        # Check if the coverage values are correct (no >10% change in coverage not caused by readstarts)
        if self.Interval.strand == '+':
            exp_curr_treatment_coverage = _new_entry.treatment_coverage - _new_entry.treatment_readstart
            if abs(exp_curr_treatment_coverage - self.buffer[-1].treatment_coverage) / float(
                    exp_curr_treatment_coverage) > 0.1:
                return False
            exp_curr_control_coverage = _new_entry.control_coverage - _new_entry.control_readstart
            if abs(exp_curr_control_coverage - self.buffer[-1].control_coverage) / float(
                    exp_curr_control_coverage) > 0.1:
                return False
        else:
            exp_next_treatment_coverage = self.buffer[-1].treatment_coverage - self.buffer[-1].treatment_readstart
            if abs(exp_next_treatment_coverage - _new_entry.treatment_coverage) / float(
                    exp_next_treatment_coverage) > 0.1:
                return False
            exp_next_control_coverage = self.buffer[-1].control_coverage - self.buffer[-1].control_readstart
            if abs(exp_next_control_coverage - _new_entry.control_coverage) / float(exp_next_control_coverage) > 0.1:
                return False

        return True

    def append_entry(self, _new_entry):
        self.buffer.append(_new_entry)
        self.Interval.end += 1
        assert self.Interval.end == _new_entry.Interval.end

    def flush(self, _fsr_store, _second_alpha_field, _second_alpha):
        if len(self.buffer) <= 1:
            return 0
        if self.Interval.strand == '+':
            # Rightmost has max coverage
            treatment_coverage = self.buffer[-1].treatment_coverage
            control_coverage = self.buffer[-1].control_coverage
        else:
            # Leftmost has max coverage
            treatment_coverage = self.buffer[0].treatment_coverage
            control_coverage = self.buffer[0].control_coverage

        treatment_readstart = min(sum([x.treatment_readstart for x in self.buffer]), treatment_coverage)
        control_readstart = min(sum([x.control_readstart for x in self.buffer]), control_coverage)

        # Now, run the fsr compute stages
        first_alpha = fsr_precompute.first_stage(tmnt_rss=treatment_readstart,
                                                 tmnt_cov=treatment_coverage,
                                                 ctrl_rss=control_readstart,
                                                 ctrl_cov=control_coverage)

        if first_alpha <= fsr_precompute.FIRST_STAGE_ALPHA_ALLOWANCE:
            second_alpha_list = fsr_precompute.second_stage(tmnt_rss=treatment_readstart,
                                                            tmnt_cov=treatment_coverage,
                                                            ctrl_rss=control_readstart,
                                                            ctrl_cov=control_coverage,
                                                            _first_alpha=first_alpha)
        else:
            second_alpha_list = None

        combined_entry = FSRTSV_entry(self.Interval, treatment_coverage, treatment_readstart, control_coverage,
                                      control_readstart,
                                      first_alpha, second_alpha_list)
        if combined_entry.first_stage_alpha <= first_alpha:
            if (combined_entry.second_stage_alpha_list and len(
                    combined_entry.second_stage_alpha_list) > _second_alpha_field):
                if combined_entry.second_stage_alpha_list[_second_alpha_field] <= _second_alpha:
                    _fsr_store.append((combined_entry.Interval.htseq_iv(),
                                       treatment_coverage,
                                       combined_entry.first_stage_alpha,
                                       combined_entry.second_stage_alpha_list))
                    return 1
        return 0


def load_fsrtsv(_fsr_store, _fsrtsv_path, _cutoff_config_path, _whitelist_array, _faidx_genome, disable_merging=False,
                disable_whitelisting=False, load_decoy_instead=False, call_lithium=False):
    def is_adjacent_to_guanine(_entry):
        if _entry.Interval.strand == '+':
            if _faidx_genome.get_fasta(_entry.Interval) == 'G':
                return True
            else:
                return False
        elif _entry.Interval.strand == '-':
            next_position = gff_parser.Interval(_entry.Interval.chrom, _entry.Interval.start + 1,
                                                _entry.Interval.end + 1,
                                                _entry.Interval.strand)
            if _faidx_genome.get_fasta(next_position) == 'G':
                return True
            else:
                return False

    def pass_whitelist(_htseqiv, _whitelist_array):
        for iv, boolean in _whitelist_array[_htseqiv].steps():
            if boolean:
                return True
        return False

    logging.debug('[{0}] Loading {1}'.format(time.ctime(), _fsrtsv_path))
    if _cutoff_config_path == 'None':
        sample_first_alpha_cutoff = FIRST_STAGE_ALPHA_ALLOWANCE
        sample_second_alpha_cutoff = FIRST_STAGE_ALPHA_ALLOWANCE
        second_alpha_field =  FIRST_STAGE_ALPHA_DICT[str(FIRST_STAGE_ALPHA_ALLOWANCE)]
    else:
        with open(_cutoff_config_path) as f:
            lines = f.readlines()
            str_sample_first_alpha_cutoff = lines[0].strip().rstrip("0")
            sample_first_alpha_cutoff = float(str_sample_first_alpha_cutoff)
            str_sample_second_alpha_cutoff = lines[0].strip().rstrip("0")
            sample_second_alpha_cutoff = float(str_sample_second_alpha_cutoff)
            second_alpha_field = FIRST_STAGE_ALPHA_DICT[str_sample_first_alpha_cutoff]

    i = 0
    with xopen(_fsrtsv_path) as f:
        fsr_merge_buffer = None
        for line in f:
            entry = parse_fsrtsv(line, call_lithium)
            if not entry:
                continue
            treatment_fsr = entry.treatment_readstart / entry.treatment_coverage
            control_fsr = entry.control_readstart / entry.control_coverage

            if load_decoy_instead:
                if treatment_fsr >= control_fsr:
                    continue
                report_coverage = entry.control_coverage
            else:
                if treatment_fsr <= control_fsr:
                    continue
                report_coverage = entry.treatment_coverage

            # New function: merge adjacent G peaks and try-again
            if entry.first_stage_alpha <= sample_first_alpha_cutoff:
                if entry.second_stage_alpha_list and len(entry.second_stage_alpha_list) > second_alpha_field:
                    if entry.second_stage_alpha_list[second_alpha_field] <= sample_second_alpha_cutoff:
                        if disable_whitelisting or pass_whitelist(entry.Interval.htseq_iv(), _whitelist_array):
                            if fsr_merge_buffer:
                                fsr_merge_buffer.flush(_fsr_store, second_alpha_field, sample_second_alpha_cutoff)
                                fsr_merge_buffer = None
                            _fsr_store.append((entry.Interval.htseq_iv(),
                                               report_coverage,
                                               entry.first_stage_alpha,
                                               entry.second_stage_alpha_list))
                            i += 1
                            continue
            if not disable_merging and not load_decoy_instead:  # OMG dont merge decoys to make even huger background
                if entry.first_stage_alpha <= FIRST_STAGE_ALPHA_ALLOWANCE and (
                        treatment_fsr - control_fsr) > FIRST_STAGE_MIN_FSR_DIFF * 2:
                    # Pass first stage but fail second -> rescue by merging
                    can_merge = False
                    if fsr_merge_buffer:
                        if (fsr_merge_buffer.is_adjacent_to(entry) and
                                (disable_whitelisting or pass_whitelist(entry.Interval.htseq_iv(), _whitelist_array))):
                            fsr_merge_buffer.append_entry(entry)
                            can_merge = True
                        else:
                            i += fsr_merge_buffer.flush(_fsr_store, second_alpha_field, sample_second_alpha_cutoff)
                            fsr_merge_buffer = None
                            can_merge = False
                    if not can_merge:
                        # is_adjacent_to_guanine only applicable to 1st nucleotide
                        if (is_adjacent_to_guanine(entry) and (
                                disable_whitelisting or pass_whitelist(entry.Interval.htseq_iv(), _whitelist_array))):
                            fsr_merge_buffer = FSRMergeBuffer(entry)
                        else:
                            pass

    logging.debug('[{0}] Loaded {1} from {2}'.format(time.ctime(), i, _fsrtsv_path))
    return


def load_all_tsv(prefix, analysis_conditions, no_of_replicates, bootstrap_no, _whitelist_array, _faidx_genome,
                 call_lithium=False):
    fsr_store = {}
    fsrtsv_manifest = "-{0}-rep{1}.bootstrap_{2:02d}.fsrtsv.{3}.gz"
    cutoff_config_manifest = "-{0}-rep{1}.bootstrap_{2:02d}.cutoff.conf"
    for sample_condition in analysis_conditions:
        fsr_store[sample_condition] = {}
        for replicate_id in range(1, no_of_replicates + 1):
            forward_in_fsrtsv = prefix + fsrtsv_manifest.format(sample_condition, replicate_id, bootstrap_no, 'forward')
            reverse_in_fsrtsv = prefix + fsrtsv_manifest.format(sample_condition, replicate_id, bootstrap_no, 'reverse')
            cutoff_config = prefix + cutoff_config_manifest.format(sample_condition, replicate_id, bootstrap_no)
            fsr_store[sample_condition][replicate_id] = {}
            fsr_store[sample_condition][replicate_id]['+'] = []
            fsr_store[sample_condition][replicate_id]['-'] = []
            load_fsrtsv(fsr_store[sample_condition][replicate_id]['+'], forward_in_fsrtsv, cutoff_config,
                        _whitelist_array, _faidx_genome, call_lithium=call_lithium)
            load_fsrtsv(fsr_store[sample_condition][replicate_id]['-'], reverse_in_fsrtsv, cutoff_config,
                        _whitelist_array, _faidx_genome, call_lithium=call_lithium)

    return fsr_store


def load_whitelist(_whitelist_array, _in_whitelist_bed):
    f = HTSeq.BED_Reader(_in_whitelist_bed)
    logging.debug('[{0}] Loading whitelist bed {1}'.format(time.ctime(), _in_whitelist_bed))
    for line in f:
        _whitelist_array[line.iv] = True


def do_combine_rts_arrays(_fsr_store, _single_cond_rts_arrays):
    combined_rts_array = HTSeq.GenomicArrayOfSets(chroms='auto', stranded=True, storage='step')
    if 'K' in _single_cond_rts_arrays:
        for replicate_id in _single_cond_rts_arrays['K'].keys():
            for strand_symbol in ('+', '-'):
                for items in _fsr_store['K'][replicate_id][strand_symbol]:
                    iv = items[0]
                    combined_rts_array[iv] += replicate_id

    if 'KPDS' in _single_cond_rts_arrays:
        for replicate_id in _single_cond_rts_arrays['KPDS'].keys():
            for strand_symbol in ('+', '-'):
                for items in _fsr_store['KPDS'][replicate_id][strand_symbol]:
                    iv = items[0]
                    combined_rts_array[iv] += replicate_id * 10
    return combined_rts_array


def construct_rts_arrays(_fsr_store, no_of_replicates, analysis_conditions):
    single_cond_rts_arrays = {}
    assert no_of_replicates <= 8
    for condition in analysis_conditions:
        if condition in _fsr_store:
            single_cond_rts_arrays[condition] = {}
            for replicate_id in range(1, no_of_replicates + 1):
                single_cond_rts_arrays[condition][replicate_id] = HTSeq.GenomicArray(chroms='auto', stranded=True,
                                                                                     typecode='i')
    if 'K' in single_cond_rts_arrays:
        for replicate_id in single_cond_rts_arrays['K'].keys():
            for strand_symbol in ('+', '-'):
                for items in _fsr_store['K'][replicate_id][strand_symbol]:
                    iv = items[0]
                    coverage = items[1]
                    single_cond_rts_arrays['K'][replicate_id][iv] = coverage

    if 'KPDS' in single_cond_rts_arrays:
        for replicate_id in single_cond_rts_arrays['KPDS'].keys():
            for strand_symbol in ('+', '-'):
                for items in _fsr_store['KPDS'][replicate_id][strand_symbol]:
                    iv = items[0]
                    coverage = items[1]
                    single_cond_rts_arrays['KPDS'][replicate_id][iv] = coverage

    #    k_rts_array, kpds_rts_array, combined_rts_array = combine_rts_arrays(single_cond_rts_arrays)
    combined_rts_array = do_combine_rts_arrays(_fsr_store, single_cond_rts_arrays)

    return combined_rts_array, single_cond_rts_arrays


def load_rg4_motif(_motif_bed_list, _all_rts_array=None):
    # As of 17-4-2019, the speedhack using RTS array input is broken. Pay the RAM tax.
    rg4_motif_information = []
    rg4_motif_array = HTSeq.GenomicArrayOfSets(chroms='auto', stranded=True, storage='step')
    rg4_threeprime_array = HTSeq.GenomicArrayOfSets(chroms='auto', stranded=True, storage='step')
    rg4_motif_serial = 0
    for in_motif_bed in _motif_bed_list:
        logging.debug('[{0}] Loading {1}'.format(time.ctime(), in_motif_bed))
        with xopen(in_motif_bed) as f:
            for line in f:
                fields = line.strip().split('\t')
                motif_chrom = fields[0]
                motif_start = int(fields[1])
                motif_end = int(fields[2])
                motif_strand = fields[4]
                motif_iv = HTSeq.GenomicInterval(motif_chrom, motif_start, motif_end, motif_strand)
                if _all_rts_array:
                    have_motif = False
                    for iv, rts_set in _all_rts_array[motif_iv].steps():
                        if rts_set:
                            have_motif = True
                        continue
                    if not have_motif:
                        continue
                rg4_intervals = [gff_parser.Interval(*RE_digest_iv_str.match(x).groups()) for x in fields[7].split(',')]

                rg4_motif_information.append(fields)
                if len(fields) >= 12:
                    rg4_3flank_intervals = [gff_parser.Interval(*RE_digest_iv_str.match(x).groups()) for x in
                                            fields[11].split(',')]
                    if motif_strand == '+':
                        iv3 = HTSeq.GenomicInterval(rg4_3flank_intervals[0].chrom, rg4_3flank_intervals[0].start,
                                                    rg4_3flank_intervals[0].start + 1 + RTS_3PRIME_ALLOWANCE, '+')
                    elif motif_strand == '-':
                        iv3 = HTSeq.GenomicInterval(rg4_3flank_intervals[0].chrom,
                                                    rg4_3flank_intervals[0].end - 1 - RTS_3PRIME_ALLOWANCE,
                                                    rg4_3flank_intervals[0].end, '-')
                    else:
                        continue
                    rg4_threeprime_array[iv3] += rg4_motif_serial
                    rg4_motif_array[iv3] += rg4_motif_serial
                if motif_strand == '+':
                    rg4iv3 = HTSeq.GenomicInterval(rg4_intervals[-1].chrom,
                                                   rg4_intervals[-1].end - 1 - RTS_3PRIME_ALLOWANCE,
                                                   rg4_intervals[-1].end,
                                                   '+')
                elif motif_strand == '-':
                    rg4iv3 = HTSeq.GenomicInterval(rg4_intervals[0].chrom, rg4_intervals[0].start,
                                                   rg4_intervals[0].start + 1 + RTS_3PRIME_ALLOWANCE, '-')
                else:
                    continue
                rg4_threeprime_array[rg4iv3] += rg4_motif_serial
                for iv in rg4_intervals:
                    rg4_motif_array[iv.htseq_iv()] += rg4_motif_serial
                rg4_motif_serial += 1
    return rg4_motif_information, rg4_motif_array, rg4_threeprime_array


def find_valid_rts_sites(combined_rts_array, analysis_conditions):
    valid_rts_sites = []
    for chrom in combined_rts_array.chrom_vectors.keys():
        for strand_symbol in ('+', '-'):
            curr_iv = None
            curr_set = None
            for iv, replicate_set in combined_rts_array[chrom][strand_symbol].steps():
                if not replicate_set:
                    continue
                if not curr_iv:
                    curr_iv = iv
                    curr_set = replicate_set
                    continue
                elif iv.start == curr_iv.end and iv.chrom == curr_iv.chrom:  # i.e. adjacent
                    assert iv.strand == curr_iv.strand
                    curr_iv = HTSeq.GenomicInterval(curr_iv.chrom, curr_iv.start, iv.end, curr_iv.strand)
                    curr_set |= replicate_set
                    continue
                elif iv.start == curr_iv.end + 1 and iv.chrom == curr_iv.chrom:  # i.e. adjacent-jump1bp
                    assert iv.strand == curr_iv.strand
                    curr_iv = HTSeq.GenomicInterval(curr_iv.chrom, curr_iv.start, iv.end, curr_iv.strand)
                    curr_set |= replicate_set
                    continue
                else:
                    if 'K' in analysis_conditions and 'KPDS' in analysis_conditions:
                        if (
                                len(curr_set & {1, 2, 3, 4}) >= 1
                                and len(curr_set & {10, 20, 30, 40}) >= 1
                                and len(curr_set) >= 3
                        ):
                            valid_rts_sites.append(tuple([curr_iv, curr_set]))

                    elif len(curr_set) >= 2:
                        valid_rts_sites.append(tuple([curr_iv, curr_set]))
                    curr_iv = iv
                    curr_set = replicate_set
            if curr_iv and curr_set:
                if 'K' in analysis_conditions and 'KPDS' in analysis_conditions:
                    if (
                            len(curr_set & {1, 2, 3, 4}) >= 1
                            and len(curr_set & {10, 20, 30, 40}) >= 1
                            and len(curr_set) >= 3
                    ):
                        valid_rts_sites.append(tuple([curr_iv, curr_set]))
                elif len(curr_set) >= 2:
                    valid_rts_sites.append(tuple([curr_iv, curr_set]))
    return valid_rts_sites


def single_cond_find_valid_rts_sites(_single_cond_rts_array):
    valid_rts_sites = []
    for chrom in _single_cond_rts_array.chrom_vectors.keys():
        for strand_symbol in ('+', '-'):
            curr_iv = None
            curr_set = None
            for iv, replicate_set in _single_cond_rts_array[chrom][strand_symbol].steps():
                if not replicate_set:
                    continue
                if not curr_iv:
                    curr_iv = iv
                    curr_set = replicate_set
                    continue
                elif iv.start == curr_iv.end and iv.chrom == curr_iv.chrom:  # i.e. adjacent
                    assert iv.strand == curr_iv.strand
                    curr_iv = HTSeq.GenomicInterval(curr_iv.chrom, curr_iv.start, iv.end, curr_iv.strand)
                    curr_set |= replicate_set
                    continue
                elif iv.start == curr_iv.end + 1 and iv.chrom == curr_iv.chrom:  # i.e. adjacent-jump1bp
                    assert iv.strand == curr_iv.strand
                    curr_iv = HTSeq.GenomicInterval(curr_iv.chrom, curr_iv.start, iv.end, curr_iv.strand)
                    curr_set |= replicate_set
                    continue
                else:
                    valid_rts_sites.append(tuple([curr_iv, curr_set]))
                    curr_iv = iv
                    curr_set = replicate_set
            if curr_iv and curr_set:
                valid_rts_sites.append(tuple([curr_iv, curr_set]))
    return valid_rts_sites


def find_g_extensions(seq):
    assert len(seq) >= 1
    accu = ['']
    accu_index = 0
    gext_accu = [False]
    gext = False
    next_char = ''
    for i in range(len(seq) - 1):
        char = seq[i]
        next_char = seq[i + 1]
        if (char == 'G' and next_char == 'G' and not gext) or (char != 'G' and gext):
            accu.append('')
            gext_accu.append(not gext)
            accu_index += 1
            gext = not gext
        accu[accu_index] += char
    if (gext and next_char == 'G') or (not gext):
        accu[accu_index] += next_char
    else:
        gext_accu.append(False)
        accu.append(next_char)
        accu_index += 1

    result = []
    gext_result = []
    for i, x in enumerate(accu):
        if not x:
            continue
        if gext_accu[i] and RE_g6.match(x):
            result += ['GG', 'G', x[4:]]
            gext_result += [True, False, True]
            continue
        result.append(x)
        gext_result.append(gext_accu[i])

    if gext_result.count(True) >= 4:
        if not gext_result[-1]:
            if len(result[-1]) > 2:  # The RTS need to be adjacent to 3'
                return None
            return result[-8:]
        else:
            return result[-7:]
    elif gext_result.count(True) >= 3:
        if not gext_result[-1]:
            if len(result[-1]) > 2:  # The RTS need to be adjacent to 3'
                return None
            return result[-6:]
        else:
            return result[-5:]
    else:
        return None


def check_novel_rg4(check_iv, _faidx_genome):
    if check_iv.strand == '+':
        g40_check_iv = HTSeq.GenomicInterval(check_iv.chrom, check_iv.end - 60, check_iv.end, '+')
        flank_iv = gff_parser.Interval(check_iv.chrom, check_iv.end, check_iv.end + 10, '+')
    elif check_iv.strand == '-':
        g40_check_iv = HTSeq.GenomicInterval(check_iv.chrom, check_iv.start, check_iv.start + 60, '-')
        flank_iv = gff_parser.Interval(check_iv.chrom, check_iv.start - 10, check_iv.start, '-')
    else:
        return

    gratio, seq = _faidx_genome.get_g_percent(g40_check_iv)
    result = find_g_extensions(_faidx_genome.get_fasta(g40_check_iv))

    verdict = 'Unknown'
    if gratio >= 40:
        if result:
            if len(result) in (5, 6):
                verdict = 'potential G-triplex & G>=40%'
            if len(result) in (7, 8):
                verdict = 'potential G-quadruplex & G>=40%'
        else:
            verdict = 'G>=40%'

    return_iv = gff_parser.Interval(g40_check_iv.chrom, g40_check_iv.start, g40_check_iv.end, g40_check_iv.strand)
    return verdict, [return_iv], [flank_iv]


def ivlist_fetch_annotation(_ivlist, _transcriptome_dict):
    output_accu = []
    for annotation_name, transcriptome in _transcriptome_dict.items():
        if 'REFSEQ' in annotation_name.upper():
            gene_names, transcript_ids, cdsutr = transcriptome.annotation_query(_ivlist)
            transcript_categories = set()
            final_transcript_ids = []
            if transcript_ids:
                for transcript_id in transcript_ids:
                    transcript_categories.add(transcript_id[0:2])
                for category in ['NM', 'NR', 'XM', 'XR']:  # Prioritize refseq annotation types:
                    if category in transcript_categories:
                        for tid in transcript_ids:
                            if tid[0:2] == category:
                                final_transcript_ids.append(tid)
                        break
            transcript_ids = final_transcript_ids
        else:
            gene_names, transcript_ids, cdsutr = transcriptome.annotation_query(_ivlist)
        if gene_names:
            output_accu.append(';'.join(gene_names))
        else:
            output_accu.append('N/A')
        if transcript_ids:
            output_accu.append(';'.join(transcript_ids))
        else:
            output_accu.append('N/A')
        if cdsutr:
            output_accu.append(';'.join(cdsutr))
        else:
            output_accu.append('N/A')
    return output_accu


def rts_analyze(_valid_rts_sites, _rg4_motif_information, _rg4_threeprime_array, _faidx_genome,
                _single_cond_rts_arrays, _transcriptome_dict, replicate_breakdown=False):
    rts_entry_store = []
    # Need to accept 2 types of input. Detection is based on available items in _single_cond_rts_arrays
    # 1. All replicates, all conditions
    # 2. One replicate, one condition

    available_conditions = _single_cond_rts_arrays.keys()

    # Process valid rts sites one by one in a large for-loop
    verdict_counter = Counter()
    for rts_iv, replicate_set in _valid_rts_sites:
        rts_query_interval = gff_parser.Interval(rts_iv.chrom, rts_iv.start, rts_iv.end, rts_iv.strand)
        rg4_motif_set = set()
        for __, rg4id_set in _rg4_threeprime_array[rts_iv].steps():
            rg4_motif_set |= rg4id_set

        # Logic: To define "One Motif"
        # 1. Find the best rG4 structural class first to simplify problem
        # 2. One RTS site can span across >1 bp, and on a transcript region with >1 valid isoforms
        # 3. To simplify 2, we first merge "identical" rG4 motifs with identical nucleotide sequence
        # 3.

        # Begin select best rG4 structural class
        final_rg4_motifs = []
        rg4_type_verdict = None
        seq_set = set()
        rg4_types = set([_rg4_motif_information[item][5] for item in sorted(rg4_motif_set)])
        for tcode in ('canonical/G3L1-7', 'longloop', 'bulges', 'two-quartet'):
            if tcode in rg4_types:
                for motif in [_rg4_motif_information[item] for item in sorted(rg4_motif_set)]:
                    if motif[5] == tcode:
                        # Check: if splicing of rG4 match rG4 (i.e. full contained)
                        motif_rg4_iv = [gff_parser.Interval(*RE_digest_iv_str.match(x).groups()) for x in
                                        motif[7].split(',')]
                        if len(motif) >= 12:
                            motif_threeprime_iv = [gff_parser.Interval(*RE_digest_iv_str.match(x).groups()) for x in
                                                   motif[11].split(',')]
                        else:
                            motif_threeprime_iv = []
                        fail_splicing_match = False
                        for i in rts_query_interval.return_position_list():
                            found = False
                            for seq_interval in motif_rg4_iv + motif_threeprime_iv:
                                if seq_interval.start <= i < seq_interval.end:
                                    found = True
                            if not found:
                                fail_splicing_match = True
                        if not fail_splicing_match:
                            seq_set.add(motif[7])
                            final_rg4_motifs.append(motif)
            if final_rg4_motifs:
                rg4_type_verdict = tcode
                break

        # End select best rG4 structural class

        if not final_rg4_motifs and not rg4_type_verdict:
            # No known motifs matching RTS site
            rg4_type_verdict, rg4_iv, best_flank3seq_iv = check_novel_rg4(rts_iv, _faidx_genome)
            rg4seq = _faidx_genome.get_fasta(rg4_iv[0])
            best_flank3seq = _faidx_genome.get_fasta(best_flank3seq_iv[0])
            best_flank3seq_length = best_flank3seq_iv[0].len
            curr_rg4_annotation_info = ivlist_fetch_annotation(rg4_iv, _transcriptome_dict)

        else:

            # Begin select best rG4 motif from best structural class

            selected_rg4_motifs = []
            if len(seq_set) == 1:
                # One sequence is left after merging -> select all motifs
                selected_rg4_motifs = final_rg4_motifs

            elif len(seq_set) > 1:
                # More than one sequence is left after merging -> select 3' most motifs
                rg4_iv_str = list(sorted(seq_set))
                first_iv = gff_parser.Interval(*RE_digest_iv_str.match(rg4_iv_str[0].split(',')[-1]).groups())
                # As of newest version, ivs are arrange by 5' to 3'
                max_end = first_iv.end
                min_start = first_iv.start
                best_index = 0
                if first_iv.strand == '+':
                    for index, iv_str in enumerate(rg4_iv_str):
                        curr_iv = gff_parser.Interval(*RE_digest_iv_str.match(iv_str).groups())
                        if curr_iv.end > max_end or (curr_iv.start < min_start and curr_iv.end == max_end):
                            best_index = index
                            max_end = curr_iv.end
                            min_start = curr_iv.start
                elif first_iv.strand == '-':
                    for index, iv_str in enumerate(rg4_iv_str):
                        curr_iv = gff_parser.Interval(*RE_digest_iv_str.match(iv_str).groups())
                        if curr_iv.start < min_start or (curr_iv.end > max_end and curr_iv.start == min_start):
                            best_index = index
                            max_end = curr_iv.end
                            min_start = curr_iv.start
                for motif in final_rg4_motifs:
                    if motif[7] == rg4_iv_str[best_index]:
                        selected_rg4_motifs.append(motif)

            assert selected_rg4_motifs
            # End select best rG4 motif from best structural class

            # Begin select most representative rG4 motif from best motifs
            best_selected_motif = selected_rg4_motifs[0]
            best_flank3seq_length = 0
            best_flank3seq_iv = None
            best_flank3seq = ''
            best_flank3seq_match_query = False
            for motif in selected_rg4_motifs:
                flank5seq, rg4seq, flank3seq = motif[8].split('|')

                if not flank3seq:
                    continue

                curr_flank3_intervals = [gff_parser.Interval(*RE_digest_iv_str.match(x).groups()) for x in
                                         motif[11].split(',')]
                curr_flank3seq_match_query = False

                for flank3_iv in curr_flank3_intervals:
                    if rts_query_interval.overlap_with(flank3_iv):
                        curr_flank3seq_match_query = True

                # Replacement logic:
                # 1. Unconditionally prefer flank3seq that overlaps with query RTS site over non-overlapping ones
                # 2. Prefer longer flank3seq
                # 3. If no rG4 motif has flank3seq (which is unlikely) -> just choose the first one
                accept_replace = False
                if best_flank3seq_match_query:
                    if curr_flank3seq_match_query:
                        if len(flank3seq) > best_flank3seq_length:
                            accept_replace = True
                else:
                    if curr_flank3seq_match_query:
                        accept_replace = True
                    elif len(flank3seq) > best_flank3seq_length:
                        accept_replace = True

                if accept_replace:
                    best_selected_motif = motif
                    best_flank3seq_iv = curr_flank3_intervals
                    best_flank3seq_length = len(flank3seq)
                    best_flank3seq = flank3seq

            rg4_iv = [gff_parser.Interval(*RE_digest_iv_str.match(x).groups()) for x in
                      best_selected_motif[7].split(',')]
            rg4seq = best_selected_motif[6]
            curr_rg4_annotation_info = ivlist_fetch_annotation(rg4_iv, _transcriptome_dict)
            # End select most representative rG4 motif

        # Begin layout of sequence/RTS site diagram
        rg4_iv_section_offset = [0] + [x.len for x in rg4_iv[:-1]]
        rts_row = OrderedDict()
        rts_row_label = OrderedDict()
        all_sample_coverage = []

        if not replicate_breakdown:
            if 'K' in available_conditions:
                rts_row['K'] = [['-'] * (len(rg4seq)), ['-'] * best_flank3seq_length]
                rts_row_label['K'] = 'rG4seq-K '
            if 'KPDS' in available_conditions:
                rts_row['KPDS'] = [['-'] * (len(rg4seq)), ['-'] * best_flank3seq_length]
                rts_row_label['KPDS'] = 'rG4seq-KPDS '
        else:
            for condition in available_conditions:
                rts_row[condition] = OrderedDict()
                rts_row_label[condition] = OrderedDict()
                for replicate in sorted(_single_cond_rts_arrays[condition].keys()):
                    rts_row[condition][replicate] = [['-'] * (len(rg4seq)), ['-'] * best_flank3seq_length]
                    rts_row_label[condition][replicate] = 'rG4seq-' + condition + '-rep{0}'.format(replicate)
        for condition in ['K', 'KPDS']:
            if condition not in available_conditions:
                continue
            for replicate in sorted(_single_cond_rts_arrays[condition].keys()):
                if not replicate_breakdown:
                    rts_store = rts_row[condition]
                else:
                    rts_store = rts_row[condition][replicate]
                for i in rts_query_interval.return_position_list():
                    coverage = _single_cond_rts_arrays[condition][replicate][
                        HTSeq.GenomicPosition(rts_iv.chrom, i, rts_iv.strand)]
                    assert coverage == 0 or coverage >= 6
                    if coverage != 0:
                        all_sample_coverage.append(coverage)
                        within_rg4 = False
                        for index, seq_interval in enumerate(rg4_iv):
                            if seq_interval.start <= i < seq_interval.end:
                                if rts_iv.strand == '+':
                                    rts_store[0][i - seq_interval.start + rg4_iv_section_offset[index]] = '*'
                                else:
                                    rts_store[0][seq_interval.end - i - 1 + rg4_iv_section_offset[index]] = '*'
                                within_rg4 = True
                                break

                        if not best_flank3seq_iv:
                            continue
                        best_flank3seq_iv_section_offset = [0] + [x.len for x in best_flank3seq_iv[:-1]]
                        for index, seq_interval in enumerate(best_flank3seq_iv):
                            if seq_interval.start <= i < seq_interval.end:
                                assert not within_rg4
                                if rts_iv.strand == '+':
                                    rts_store[1][i - seq_interval.start + best_flank3seq_iv_section_offset[index]] = '*'
                                else:
                                    rts_store[1][
                                        seq_interval.end - i - 1 + best_flank3seq_iv_section_offset[index]] = '*'
                                break

        flank3_displaylength = min(len(best_flank3seq), 10)
        verdict_counter[rg4_type_verdict] += 1

        # First, print the sequence and stuff which are fixed anyways
        seq_diag_store = []

        if rts_iv.strand == '+':
            seq_diag_store.append('{:>16}'.format('5\'- ') + rg4seq + best_flank3seq[0:flank3_displaylength] + ' -3\'')
        elif rts_iv.strand == '-':
            seq_diag_store.append('{:>16}'.format('5\'- ') + rg4seq + best_flank3seq[0:flank3_displaylength] + ' -3\'')

        # Then, plot the RTS site asterisks
        for condition in ['K', 'KPDS']:
            if condition not in available_conditions:
                continue
            for replicate in sorted(_single_cond_rts_arrays[condition].keys()):
                if replicate_breakdown:
                    rts_store = rts_row[condition][replicate]
                    rts_label = rts_row_label[condition][replicate]
                else:
                    rts_store = rts_row[condition]
                    rts_label = rts_row_label[condition]
                if rts_iv.strand == '+':
                    seq_diag_store.append('{:>16}'.format(rts_label) + ''.join(rts_store[0]) + ''.join(
                        rts_store[1][0:flank3_displaylength]))
                elif rts_iv.strand == '-':
                    seq_diag_store.append('{:>16}'.format(rts_label) + ''.join(rts_store[0]) + ''.join(
                        rts_store[1][0:flank3_displaylength]))
                if not replicate_breakdown:
                    break
        seq_diag = '\n'.join(seq_diag_store)
        rg4_interval_str = ','.join([x.return_coordinates() for x in rg4_iv])
        rg4_has_splicing = True if len(rg4_iv) > 1 else False
        # Finally, gather the pieces and report
        this_rts_entry = [rts_iv.chrom,
                          rts_iv.start,
                          rts_iv.end,
                          rg4_type_verdict,
                          int(round(sum(all_sample_coverage) / len(all_sample_coverage))),
                          rts_iv.strand,
                          rg4_interval_str,
                          rg4_has_splicing,
                          seq_diag, ] + curr_rg4_annotation_info
        rts_entry_store.append(this_rts_entry)
    return rts_entry_store


def export_csv(_rts_analyze_result, output_path):
    with open(output_path, 'w', newline='') as csvfw:
        csv_writer = csv.writer(csvfw, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csv_writer.writerow(RTS_csv_header)
        for line in _rts_analyze_result:
            csv_writer.writerow(line)


def main(prefix, analysis_conditions, no_of_replicates, bootstrap_no, genome_fasta, whitelist_bed_list, motif_bed_list,
         annotation_names, annotation_gffs, call_lithium=False):
    gff_dict = OrderedDict()
    for i in range(len(annotation_names)):
        gff_dict[annotation_names[i]] = annotation_gffs[i]

    logging.info('[{0}] Loading whitelists...'.format(time.ctime()))
    faidx_genome = FaidxGenome(genome_fasta)
    whitelist_array = HTSeq.GenomicArray(chroms='auto', stranded=True, storage='step', typecode='b')
    for in_whitelist_bed in whitelist_bed_list:
        load_whitelist(whitelist_array, in_whitelist_bed)

    logging.info('[{0}] Loading annotations...'.format(time.ctime()))
    transcriptome_dict = OrderedDict()
    for annotation_name, gff_path in sorted(gff_dict.items()):
        transcriptome_dict[annotation_name] = gff_parser.Transcriptome()
        transcriptome_dict[annotation_name].load_annotation(gff_path)
        transcriptome_dict[annotation_name].finalize()
        # Update CSV header
        RTS_csv_header.append('{0}_gene_names'.format(annotation_name))
        RTS_csv_header.append('{0}_coding_regions'.format(annotation_name))
        RTS_csv_header.append('{0}_transcript_ids'.format(annotation_name))

    logging.info('[{0}] Loading RTS site data...'.format(time.ctime()))
    fsr_store = load_all_tsv(prefix, analysis_conditions, no_of_replicates, bootstrap_no, whitelist_array, faidx_genome,
                             call_lithium=call_lithium)

    # Full combined
    full_combined_rts_array, single_cond_rts_arrays = construct_rts_arrays(fsr_store, no_of_replicates,
                                                                           analysis_conditions)
    rg4_motif_information, rg4_motif_array, rg4_threeprime_array = load_rg4_motif(motif_bed_list,
                                                                                  _all_rts_array=None)
    if call_lithium:  # Reverse the treatment and control group, call Lithiums-specific RTS instead
        prefix += '.reverse_lithium'

    if no_of_replicates > 1:
        logging.info('[{0}] Analyzing RTS combining all replicates...'.format(time.ctime()))
        full_valid_rts_sites = find_valid_rts_sites(full_combined_rts_array, analysis_conditions)
        full_result = rts_analyze(full_valid_rts_sites, rg4_motif_information, rg4_threeprime_array, faidx_genome,
                                  single_cond_rts_arrays, transcriptome_dict)
        full_csv = prefix + '.rG4_list.combined.csv'
        export_csv(full_result, full_csv)
        del full_result
        full_breakdown_result = rts_analyze(full_valid_rts_sites, rg4_motif_information, rg4_threeprime_array,
                                            faidx_genome,
                                            single_cond_rts_arrays, transcriptome_dict, replicate_breakdown=True)
        full_breakdown_csv = prefix + '.rG4_list.combined_breakdown.csv'
        export_csv(full_breakdown_result, full_breakdown_csv)
        del full_breakdown_result

    else:
        logging.info('[{0}] Analyzing RTS independently for each replicate...'.format(time.ctime()))
        for condition in fsr_store.keys():
            for replicate in fsr_store[condition].keys():
                current_single_cond_rts_array = single_cond_rts_arrays[condition][replicate]
                current_valid_rts_sites = single_cond_find_valid_rts_sites(current_single_cond_rts_array)
                current_result = rts_analyze(current_valid_rts_sites, rg4_motif_information, rg4_threeprime_array,
                                             faidx_genome, {condition: {replicate: current_single_cond_rts_array}},
                                             transcriptome_dict)
                current_csv = prefix + '.rG4_list.{0}-rep{1}.csv'.format(condition, replicate)
                export_csv(current_result, current_csv)

    logging.info('[{0}] Finished analyzing RTS'.format(time.ctime()))
    return


if __name__ == "__main__":
    exit()
