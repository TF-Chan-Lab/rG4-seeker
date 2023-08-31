"""
**************************************************************************
**  rG4-seeker
**  A pipeline for processing and analyzing rG4-seq data
**
**  Version 1.0 -- August 10, 2019
**
**  Copyright (C) 2019 by Eugene Yui-Ching Chow, Ting-Fung Chan, All rights reserved.
**  Contact:  eugene.chow@link.cuhk.edu.hk, tf.chan@cuhk.edu.hk
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

from collections import OrderedDict, namedtuple
import multiprocessing
import HTSeq
from xopen import xopen
import pyfaidx
import array
import re
import time
import os
import logging
# rg4seeker imports
from rg4seeker import gff_parser

# Global Variables
SEQ_FILTER_SEARCH_RANGE = 10
EXON_END_BLACKLIST_RANGE = 10
FLANKING_OFFSET = 50
gff_parser.SEQ_FILTER_SEARCH_RANGE = 10
gff_parser.EXON_END_BLACKLIST_RANGE = 10

CHROMOSOME_ESCAPE_LIST = frozenset()


def merge_queue_to_gz(out_file_path, _file_output_queues, _chromosome_order):
    bed_line = '{0}\t{1}\t{2}\t.\t0\t{3}\n'

    def write_line(_fw, _chrom, _range, _strand, _bed_line):
        fw.write(_bed_line.format(_chrom, _range[0], _range[1], _strand).encode('UTF-8'))

    line_count = 0
    with xopen(out_file_path, mode='wb') as fw:
        for _chrom in _chromosome_order:
            curr_forward_queue = _file_output_queues[_chrom + '+']
            curr_reverse_queue = _file_output_queues[_chrom + '-']
            forward_range = curr_forward_queue.get()
            reverse_range = curr_reverse_queue.get()

            while True:
                line_count += 1
                if not forward_range and not reverse_range:
                    break
                if not reverse_range:
                    while True:
                        forward_range = curr_forward_queue.get()
                        if not forward_range:
                            break
                        write_line(fw, _chrom, forward_range, '+', bed_line)
                    break

                if not forward_range:
                    while True:
                        reverse_range = curr_reverse_queue.get()
                        if not reverse_range:
                            break
                        write_line(fw, _chrom, reverse_range, '-', bed_line)
                    break

                write_forward = True
                if forward_range[0] > reverse_range[0] or (
                        forward_range[0] == reverse_range[0] and forward_range[1] > reverse_range[1]):
                    write_forward = False
                if write_forward:
                    write_line(fw, _chrom, forward_range, '+', bed_line)
                    forward_range = curr_forward_queue.get()
                    continue
                else:
                    write_line(fw, _chrom, reverse_range, '-', bed_line)
                    reverse_range = curr_reverse_queue.get()
                    continue
    return


class FastaGenome:
    def __init__(self, _input_genome_fasta_path):
        self.pyfaidx_obj = pyfaidx.Fasta(_input_genome_fasta_path)

    def retrieve_sequence(self, _interval):
        if _interval.strand == '+':
            return str(self.pyfaidx_obj[_interval.chrom][_interval.start:_interval.end].seq)
        elif _interval.strand == '-':
            return str(self.pyfaidx_obj[_interval.chrom][_interval.start:_interval.end].reverse.complement.seq)

    def chromosome_order(self):
        return self.pyfaidx_obj.keys()


class Whitelister:
    def __init__(self, _chromosome_order, _output_whitelist_bed_prefix, _input_gff, _fasta_path, overwrite=True):
        self.chromosome_order = _chromosome_order
        self.prefix = _output_whitelist_bed_prefix
        self.output_whitelist_bed_path = self.prefix + '.rG4_whitelist.bed.gz'
        self.input_gff = _input_gff
        self.fasta_path = _fasta_path
        self.overwrite = overwrite
        self.wl_thread_list = []
        for chrom in self.chromosome_order:
            new_wl_thread = WLThread(chrom, self.prefix, self.input_gff, self.fasta_path)
            self.wl_thread_list.append(new_wl_thread)

    def run(self, processes=None):
        if not self.overwrite and os.path.exists(self.output_whitelist_bed_path) and (
                os.path.isfile(self.output_whitelist_bed_path) or os.path.islink(self.output_whitelist_bed_path)):
            logging.debug(
                '[{0}] whitelist already exists for annotation {1} at {2}'.format(time.ctime(), self.input_gff,
                                                                                  self.output_whitelist_bed_path))
            return self.output_whitelist_bed_path
        logging.info('[{0}] Begin building whitelist for annotation {1}'.format(time.ctime(), self.input_gff, ))
        with multiprocessing.Pool(processes) as pool:
            multiple_results = [pool.apply_async(wl_thread.run, args=()) for wl_thread in self.wl_thread_list]
            temp_beds = [res.get() for res in multiple_results]

        with xopen(self.output_whitelist_bed_path, 'wb', compresslevel=9) as fw:
            for temp_bed_file in temp_beds:
                if temp_bed_file:
                    with xopen(temp_bed_file, 'rb') as f:
                        for line in f:
                            fw.write(line)
                    os.remove(temp_bed_file)
        logging.info('[{0}] Completed building whitelist for annotation {1}'.format(time.ctime(), self.input_gff))
        return self.output_whitelist_bed_path


class WLThread:
    def __init__(self, chrom, prefix, _input_gff, _fasta_path):
        self.chrom = chrom
        self.output_temp_bed = prefix + '.{0}.rG4_whitelist.bed.gz'.format(self.chrom)
        self.results = {'+': [], '-': []}
        self.bed_output = {'+': [], '-': []}
        self.input_gff = _input_gff
        self.fasta_path = _fasta_path

    @staticmethod
    def _check_for_whitelist(seq):
        if not len(seq) >= 3:
            logging.warning(
                'Warning: whitelist checking sequencing length not equal to {0}'.format(SEQ_FILTER_SEARCH_RANGE))
            return False
        if not ('G' in seq[-3:]):
            return False
        if not ('GG' in seq or seq.count('G') >= 3):
            return False
        return True

    def _exon_whitelist(self, interval, fasta_sequence):
        if not len(fasta_sequence) == interval.len:
            logging.critical('Fasta subsequence length does not match interval {0}'.format(interval))
            raise Exception('Fasta subsequence length does not match interval'.format(interval))

        for pos in range(EXON_END_BLACKLIST_RANGE, interval.len):
            if self._check_for_whitelist(fasta_sequence[max(0, pos - SEQ_FILTER_SEARCH_RANGE):pos]):
                if interval.strand == '+':
                    self.results['+'].append(interval.start + pos)
                elif interval.strand == '-':
                    self.results['-'].append(interval.end - 1 - pos)

        return

    def _junction_whitelist(self, _seq_filter_search_range, _junction_genomic_position_array, _junction_fasta_sequence):
        for pos in range(_seq_filter_search_range, len(_junction_fasta_sequence)):
            if self._check_for_whitelist(_junction_fasta_sequence[max(0, pos - _seq_filter_search_range):pos]):
                curr_position = _junction_genomic_position_array[pos]
                strand = curr_position[1]
                self.results[strand].append(curr_position[0])
        return

    def _sort_results(self):
        for strand in ('+', '-'):
            if not self.results[strand]:
                continue
            self.results[strand] = sorted(list(frozenset(self.results[strand])))
            current_strand_result = self.results[strand]

            bed_span_start = 0
            bed_span_end = 0
            bed_span_endval = current_strand_result[bed_span_end] + 1
            for i in range(1, len(current_strand_result)):
                if current_strand_result[i] == bed_span_endval:
                    bed_span_end += 1
                    bed_span_endval += 1
                else:
                    self.bed_output[strand].append((current_strand_result[bed_span_start], bed_span_endval))
                    bed_span_start = i
                    bed_span_end = i
                    bed_span_endval = current_strand_result[bed_span_end] + 1
            self.bed_output[strand].append((current_strand_result[bed_span_start], bed_span_endval))

    def _export_bedfile(self):
        with xopen(self.output_temp_bed, 'wb') as fw:
            for strand in ('+', '-'):
                for start, end in self.bed_output[strand]:
                    fw.write('{0}\t{1}\t{2}\t.\t0\t{3}\n'.format(self.chrom, start, end, strand).encode('UTF-8'))

    def run(self):
        logging.debug('Start whitelisting chromosome {0}'.format(self.chrom))
        fasta_genome = FastaGenome(self.fasta_path)
        input_transcriptome = gff_parser.Transcriptome()
        input_transcriptome.load_annotation(self.input_gff, chromosome_restriction=self.chrom)
        input_transcriptome.finalize()
        start = time.time()
        if self.chrom in input_transcriptome.chromosome_dict and len(
                input_transcriptome.chromosome_dict[self.chrom]) > 0:
            logging.debug(
                'Chromosome {0} has {1} genes'.format(self.chrom, len(input_transcriptome.chromosome_dict[self.chrom])))
        else:
            logging.debug('Chromosome {0} has no annotated genes'.format(self.chrom))
            return None

        exon_array = HTSeq.GenomicArray(chroms=[self.chrom], typecode='b', storage='step', stranded=True)
        for gene_obj in input_transcriptome.chromosome_dict[self.chrom]:
            for transcript_id, transcript_obj in gene_obj.transcript_dict.items():
                for exon_obj in transcript_obj.exon_list:
                    _htseq_iv = exon_obj.get_iv()
                    if _htseq_iv:
                        exon_array[_htseq_iv] = True

                for junction_obj_list in transcript_obj.dump_junction_groups():
                    junction_fasta_sequence = ''
                    junction_genomic_position_array = []
                    # Stores (pos,strand) tuples, since junctions can link opposite stranded exons
                    for junction_order, junction in enumerate(junction_obj_list):
                        iv_left, iv_right = junction.get_whitelist_iv_pair()
                        if junction_order == 0:
                            # Only first junction in group can add the LHS components, otherwise will duplicate:
                            junction_fasta_sequence += fasta_genome.retrieve_sequence(iv_left)
                            for item in iv_left.parse_generate_pos_strand_tuple():
                                junction_genomic_position_array.append(item)
                        junction_fasta_sequence += fasta_genome.retrieve_sequence(iv_right)
                        for item in iv_right.parse_generate_pos_strand_tuple():
                            junction_genomic_position_array.append(item)
                    if not len(junction_fasta_sequence) == len(junction_genomic_position_array):
                        raise Exception('Fasta sub-sequence length does not match interval')
                    self._junction_whitelist(SEQ_FILTER_SEARCH_RANGE,
                                             junction_genomic_position_array,
                                             junction_fasta_sequence, )

        for strand in ('+', '-'):
            for _htseq_iv, _bool in exon_array.chrom_vectors[self.chrom][strand].steps():
                if _bool:
                    curr_interval = gff_parser.Interval(_htseq_iv.chrom, _htseq_iv.start, _htseq_iv.end,
                                                        _htseq_iv.strand)
                    curr_fasta = fasta_genome.retrieve_sequence(curr_interval)
                    self._exon_whitelist(curr_interval, curr_fasta)

        self._sort_results()
        self._export_bedfile()
        logging.debug('Completed whitelisting chromosome {0} in {1}s'.format(self.chrom, time.time() - start))
        return self.output_temp_bed


def whitelisting_driver(input_annotation_gff_path, input_genome_fasta_path, output_whitelist_bed_prefix, processes,
                        overwrite):
    # Whitelisting
    chromosome_set = set()
    chromosome_order = []

    pyfaidx_obj = pyfaidx.Fasta(input_genome_fasta_path)
    with xopen(input_annotation_gff_path)as f:
        for line in f:
            if line[0] == '#':
                continue
            chrom = line.strip().split('\t')[0]
            if chrom not in pyfaidx_obj:
                continue
            if chrom not in chromosome_set:
                chromosome_set.add(chrom)
                chromosome_order.append(chrom)
    whitelisting = Whitelister(chromosome_order, output_whitelist_bed_prefix, input_annotation_gff_path,
                               input_genome_fasta_path, overwrite)
    return whitelisting.run(processes)


def rg4_finding_driver(input_annotation_gff_path, input_genome_fasta_path, output_motif_bed_prefix, processes,
                       overwrite):
    chromosome_set = set()
    chromosome_order = []
    pyfaidx_obj = pyfaidx.Fasta(input_genome_fasta_path)
    with xopen(input_annotation_gff_path)as f:
        for line in f:
            if line[0] == '#':
                continue
            chrom = line.strip().split('\t')[0]
            if chrom not in pyfaidx_obj:
                continue
            if chrom not in chromosome_set:
                chromosome_set.add(chrom)
                chromosome_order.append(chrom)
    rg4finder = RG4Finder(chromosome_order, output_motif_bed_prefix, input_annotation_gff_path, input_genome_fasta_path,
                          overwrite)
    return rg4finder.run(processes)


class RG4Finder:
    def __init__(self, _chromosome_order, _output_motif_bed_prefix, _input_gff, _fasta_path, overwrite=True):
        self.chromosome_order = _chromosome_order
        self.prefix = _output_motif_bed_prefix
        self.output_motif_bed_path = self.prefix + '.rG4_motif.bed.gz'
        self.input_gff = _input_gff
        self.fasta_path = _fasta_path

        self.rf_thread_list = []
        logging.debug(
            "Begin rG4 motif finder for chromosomes: {0}".format(','.join([str(x) for x in self.chromosome_order])))
        for chrom in self.chromosome_order:
            new_rf_thread = RFThread(chrom, self.prefix, self.input_gff, self.fasta_path)
            self.rf_thread_list.append(new_rf_thread)
        self.overwrite = overwrite

    def run(self, processes=None):
        if not self.overwrite and os.path.exists(self.output_motif_bed_path) and (
                os.path.isfile(self.output_motif_bed_path) or os.path.islink(self.output_motif_bed_path)):
            logging.warning(
                '[{0}] rG4 list already exists for annotation {1} at {2}'.format(time.ctime(), self.input_gff,
                                                                                 self.output_motif_bed_path))
            return self.output_motif_bed_path
        logging.info('[{0}] Begin rG4 finding for annotation {1}'.format(time.ctime(), self.input_gff, ))
        with multiprocessing.Pool(processes) as pool:
            multiple_results = [pool.apply_async(wl_thread.run, args=()) for wl_thread in self.rf_thread_list]
            temp_beds = [res.get() for res in multiple_results]

        with xopen(self.output_motif_bed_path, 'wb', compresslevel=9) as fw:
            for temp_bed_file in temp_beds:
                if temp_bed_file:
                    with xopen(temp_bed_file, 'rb') as f:
                        for line in f:
                            fw.write(line)
                    os.remove(temp_bed_file)
        logging.info('[{0}] Completed rG4 finding for annotation {1}'.format(time.ctime(), self.input_gff, ))
        return self.output_motif_bed_path


class RFThread:
    def __init__(self, chrom, prefix, _input_gff, _fasta_path):
        self.chrom = chrom
        self.output_temp_motif = prefix + '.{0}.rG4_motif.bed.gz'.format(self.chrom)
        self.motif_output = {'+': [], '-': []}
        self.cluster_output = {'+': [], '-': []}
        self.input_gff = _input_gff
        self.fasta_path = _fasta_path

    @staticmethod
    def rg4_regular_expression_generate():

        gtract3 = '(GGG)'
        gtract2 = '(GG)'
        loop1to7 = '([ATCGN]{1,7})'
        loop1to9 = '([ATCGN]{1,9})'
        loop8to12 = '([ATCGN]{8,12})'
        loop8to21 = '([ATCGN]{8,21})'
        gtract_b1 = '(G[ATCN]GG|GG[ATCN]G)'
        gtract_b7 = '(GG[ATCN]{1,7}G|G[ATCN]{1,7}GG)'

        motifs = OrderedDict()
        motifs['canonical'] = (gtract3 + loop1to7 + gtract3 + loop1to7 + gtract3 + loop1to7 + gtract3 + '$', 1)
        motifs['longloop_leftlateral'] = (
            gtract3 + loop8to12 + gtract3 + loop1to7 + gtract3 + loop1to7 + gtract3 + '$', 2)
        motifs['longloop_rightlateral'] = (
            gtract3 + loop1to7 + gtract3 + loop1to7 + gtract3 + loop8to12 + gtract3 + '$', 2)
        motifs['longloop_central'] = (gtract3 + loop1to7 + gtract3 + loop8to21 + gtract3 + loop1to7 + gtract3 + '$', 2)
        motifs['1nt_bulges'] = (gtract_b1 + loop1to9 + gtract_b1 + loop1to9 + gtract_b1 + loop1to9 + gtract_b1 + '$', 3)
        motifs['multi_nt_bulge_tract1'] = (
            gtract_b7 + loop1to9 + gtract3 + loop1to9 + gtract3 + loop1to9 + gtract3 + '$', 3)
        motifs['multi_nt_bulge_tract2'] = (
            gtract3 + loop1to9 + gtract_b7 + loop1to9 + gtract3 + loop1to9 + gtract3 + '$', 3)
        motifs['multi_nt_bulge_tract3'] = (
            gtract3 + loop1to9 + gtract3 + loop1to9 + gtract_b7 + loop1to9 + gtract3 + '$', 3)
        motifs['multi_nt_bulge_tract4'] = (
            gtract3 + loop1to9 + gtract3 + loop1to9 + gtract3 + loop1to9 + gtract_b7 + '$', 3)
        motifs['two-quartet'] = (gtract2 + loop1to9 + gtract2 + loop1to9 + gtract2 + loop1to9 + gtract2 + '$', 4)

        """
        order = 0
        for all_gtract in sorted(frozenset(permutations([gtract3] * 3 + [gtract2to3], 4))):
            for all_loops in sorted(frozenset(permutations([loop1to7] * 2 + [loop8to28]))):
                order += 1
                motifs['ex_longloop_motif_' + str(order)] = (
                    ''.join([x for pair in zip_longest(all_gtract, all_loops) for x in pair if x]) + '$', 5)
        """
        return tuple([tuple([re.compile(re_pattern), rg4_priority]) for re_pattern, rg4_priority in motifs.values()])

    @staticmethod
    def merge_positions(_position_array, _strand_array, array_start_pos, array_end_pos, _chrom):
        output_interval_list = []
        curr_start_pos = array_start_pos
        curr_end_pos = array_start_pos
        curr_strand = _strand_array[array_start_pos]
        curr_interval_length = 1
        for pos in range(array_start_pos + 1, array_end_pos):
            if _strand_array[pos] == '+' and curr_strand == '+':
                if _position_array[pos] == _position_array[curr_end_pos] + 1:
                    curr_end_pos = pos
                    curr_interval_length += 1
                    continue
            if _strand_array[pos] == '-' and curr_strand == '-':
                if _position_array[pos] == _position_array[curr_end_pos] - 1:
                    curr_end_pos = pos
                    curr_interval_length += 1
                    continue

            # 2 positions are not adjacent, flush buffer
            if curr_strand == '+':
                output_interval_list.append(
                    gff_parser.Interval(_chrom, _position_array[curr_start_pos], _position_array[curr_end_pos] + 1,
                                        curr_strand))
            elif curr_strand == '-':
                output_interval_list.append(
                    gff_parser.Interval(_chrom, _position_array[curr_end_pos], _position_array[curr_start_pos] + 1,
                                        curr_strand))
            assert output_interval_list[-1].len == curr_interval_length
            curr_start_pos = pos
            curr_end_pos = pos
            curr_strand = _strand_array[pos]
            curr_interval_length = 1

        # Flush Last Buffer
        if curr_strand == '+':
            output_interval_list.append(
                gff_parser.Interval(_chrom, _position_array[curr_start_pos], _position_array[curr_end_pos] + 1,
                                    curr_strand))
        elif curr_strand == '-':
            output_interval_list.append(
                gff_parser.Interval(_chrom, _position_array[curr_end_pos], _position_array[curr_start_pos] + 1,
                                    curr_strand))
        return output_interval_list

    @staticmethod
    def get_interval_list_boundary(_interval_list):
        _chrom = _interval_list[0].chrom
        _leftmost_pos = _interval_list[0].start
        _rightmost_pos = _interval_list[0].end - 1
        _overall_strand = _interval_list[0].strand
        for next_interval in _interval_list[1:]:
            assert next_interval.chrom == _chrom
            if next_interval.start < _leftmost_pos:
                _leftmost_pos = next_interval.start
            if next_interval.end - 1 > _rightmost_pos:
                _rightmost_pos = next_interval.end - 1
            if next_interval.strand != _overall_strand:
                _overall_strand = '.'
        return _chrom, _leftmost_pos, _rightmost_pos, _overall_strand

    @staticmethod
    def get_rg4_type_name(_rg4_type_id):
        rg4_type_id_to_name = {1: 'canonical/G3L1-7',
                               2: 'longloop',
                               3: 'bulges',
                               4: 'two-quartet',
                               5: 'non-canonical', }
        if 1 <= _rg4_type_id <= 5:
            return rg4_type_id_to_name[_rg4_type_id]
        else:
            return 'unknown'

    def gene_find_rg4(self, _gene_obj, _transcript_ids, _transcript_items, _rg4_patterns, _chrom):
        # Finds all rG4 motifs and cluster per transcript and retain parent cluster information
        # Then attempt to merge clusters regardless of splicing (only by genomic interval)
        # But retain the exact per-transcript cluster coordinate

        Local_rg4_motif_object = namedtuple('local_rg4_motif_object',
                                            ['motif_id', 'start', 'end', 'priority', 'flank_start', 'flank_end'])
        Global_rg4_cluster_object = namedtuple('global_rg4_cluster_object',
                                               ['parent_tid', 'cluster_id', 'interval_list'])
        Global_rg4_motif_object = namedtuple('global_rg4_motif_object',
                                             ['interval_list', 'rg4_type_id', 'sequence', 'cluster_id_tuple_list',
                                              'left_flank_interval_list', 'right_flank_interval_list',
                                              'flank_sequence'])

        _RG4_MIN_LEN = 11
        _RG4_MAX_LEN = 54
        _RTS_3PRIME_ALLOWANCE = 1
        _NON_CANONICAL_3PRIME_ALLOWANCE = 2

        per_transcript_rg4_clusters = {}
        per_cluster_rg4_motifs = {}

        for tid in _transcript_ids:
            next_local_motif_id = 1
            next_local_cluster_id = 1
            for sequence, exon_list in _transcript_items[tid]:
                assert len(sequence) == sum([x.interval.len for x in exon_list])
                position_array = array.array('l',
                                             [pos for sublist in [x.interval.return_position_list() for x in exon_list]
                                              for pos in sublist])
                strand_array = array.array('u', [_strand for sublist in
                                                 [[x.interval.strand] * x.interval.len for x in exon_list]
                                                 for _strand in sublist])
                rg4_motif_array = [None] * len(sequence)
                rg4_cluster_array = HTSeq.GenomicArray(typecode='i', stranded=False, storage='ndarray',
                                                       chroms={tid: len(sequence)})
                rg4_cluster_array[HTSeq.GenomicInterval(tid, 0, len(sequence), '.')] = 0
                non_canonical_list = []

                # rG4 min length = 11, max length = 54
                for three_prime_end_pos in range(len(sequence) - 1, _RG4_MIN_LEN - 1,
                                                 -1):  # Equivalent to list(range(RG4_MIN_LEN,len(sequence))[::-1]
                    if sequence[three_prime_end_pos - 1] != 'G':
                        continue
                    sub_sequence = sequence[max(0, three_prime_end_pos - _RG4_MAX_LEN):three_prime_end_pos]
                    if sub_sequence.count('G') < 8:
                        continue

                    best_start = three_prime_end_pos
                    best_priority = 9999
                    most_distant_start = three_prime_end_pos
                    have_rg4 = False
                    for re_pattern, rg4_priority in _rg4_patterns:  # Sliding window
                        result = re_pattern.search(sub_sequence)
                        if result:
                            have_rg4 = True
                            curr_start = three_prime_end_pos - len(result.group())
                            if rg4_priority <= best_priority:
                                best_priority = rg4_priority
                                if best_start == three_prime_end_pos or curr_start > best_start:
                                    best_start = curr_start
                            if curr_start < most_distant_start:
                                most_distant_start = curr_start
                    if have_rg4:
                        # Logic: merge all canonicals into clusters. for non-canonicals never merge into clusters,
                        # retain only if their 3' end does not overlap with any canonical motifs' 3' end
                        # Note: we record the BEST motif but cluster the LONGEST motif.
                        if best_priority <= 4:  # Canonical rG4s
                            curr_local_rg4_motif_obj = Local_rg4_motif_object(
                                motif_id=next_local_motif_id,
                                start=best_start,
                                end=three_prime_end_pos,
                                priority=best_priority,
                                flank_start=max(0, best_start - FLANKING_OFFSET),
                                flank_end=min(three_prime_end_pos + FLANKING_OFFSET, len(sequence) - 1),
                            )
                            # The idea is, we try to allow only best rG4 motif for each G gtract.
                            # But of course there are bulged motifs that need some extra "care"

                            pos = three_prime_end_pos + _RTS_3PRIME_ALLOWANCE
                            while True:
                                if pos < len(sequence):
                                    if not rg4_motif_array[pos]:
                                        rg4_motif_array[pos] = curr_local_rg4_motif_obj
                                    elif rg4_motif_array[pos].priority > best_priority:
                                        j = pos
                                        while True:
                                            rg4_motif_array[j] = curr_local_rg4_motif_obj
                                            j += 1
                                            if not j < len(sequence):
                                                break
                                            if sequence[j] != 'G':
                                                break
                                    else:
                                        break
                                pos -= 1
                                if (pos < three_prime_end_pos - _RTS_3PRIME_ALLOWANCE) and sequence[pos] != 'G':
                                    break

                            rg4_cluster_array[
                                HTSeq.GenomicInterval(tid, most_distant_start, three_prime_end_pos, '.')] = 1
                        else:  # Non-canonical rG4s, delay until complete search to find overlaps with clusters
                            non_canonical_list.append(
                                Local_rg4_motif_object(motif_id=next_local_motif_id,
                                                       start=best_start,
                                                       end=three_prime_end_pos,
                                                       priority=best_priority,
                                                       flank_start=max(0, best_start - FLANKING_OFFSET),
                                                       flank_end=min(three_prime_end_pos + FLANKING_OFFSET,
                                                                     len(sequence) - 1),
                                                       )
                            )
                        next_local_motif_id += 1

                # Check non_canonical motifs, add back good ones back
                for non_canonical_rg4_motif_object in non_canonical_list:
                    overlap_with_canonical_motifs = False
                    for pos in range(max(0, non_canonical_rg4_motif_object.start - _NON_CANONICAL_3PRIME_ALLOWANCE - 1),
                                     non_canonical_rg4_motif_object.end + _NON_CANONICAL_3PRIME_ALLOWANCE):
                        if not pos < len(sequence):
                            continue
                        if rg4_motif_array[pos]:
                            overlap_with_canonical_motifs = True
                            break
                    if not overlap_with_canonical_motifs:
                        for pos in range(non_canonical_rg4_motif_object.end - _RTS_3PRIME_ALLOWANCE,
                                         non_canonical_rg4_motif_object.end + _RTS_3PRIME_ALLOWANCE + 1):
                            if not pos < len(sequence):
                                continue
                            if not rg4_motif_array[pos]:
                                rg4_motif_array[pos] = non_canonical_rg4_motif_object
                                rg4_cluster_array[
                                    HTSeq.GenomicInterval(tid, non_canonical_rg4_motif_object.start,
                                                          non_canonical_rg4_motif_object.end, '.')] = 1
                            else:
                                raise Exception('Assigning non canonical rg4 motif object on occupied positions')

                # Transform clusters into interval lists and add ID
                per_transcript_rg4_clusters[tid] = []
                for current_htseq_iv, val in list(rg4_cluster_array.chrom_vectors[tid]['.'].steps()):
                    if val == 1:
                        rg4_cluster_array[current_htseq_iv] = next_local_cluster_id
                        per_transcript_rg4_clusters[tid].append(Global_rg4_cluster_object(tid,
                                                                                          next_local_cluster_id,
                                                                                          self.merge_positions(
                                                                                              position_array,
                                                                                              strand_array,
                                                                                              current_htseq_iv.start,
                                                                                              current_htseq_iv.end,
                                                                                              _chrom)))
                        next_local_cluster_id += 1


                # Transform motifs into intervals lists, assign parent and add ID
                processed_local_motifs = set()
                for curr_local_rg4_motif_obj in rg4_motif_array:
                    if not curr_local_rg4_motif_obj:
                        continue
                    if curr_local_rg4_motif_obj.motif_id not in processed_local_motifs:
                        cluster_id = rg4_cluster_array[
                            HTSeq.GenomicPosition(tid, curr_local_rg4_motif_obj.end - 1, '.')]
                        assert cluster_id != 0 or curr_local_rg4_motif_obj.priority > 4
                        id_tuple = (tid, cluster_id)
                        if id_tuple not in per_cluster_rg4_motifs:
                            per_cluster_rg4_motifs[id_tuple] = []

                        curr_left_flank_interval_list = []
                        curr_right_flank_interval_list = []
                        if curr_local_rg4_motif_obj.flank_start != curr_local_rg4_motif_obj.start:
                            curr_left_flank_interval_list = self.merge_positions(position_array,
                                                                                 strand_array,
                                                                                 curr_local_rg4_motif_obj.flank_start,
                                                                                 curr_local_rg4_motif_obj.start, _chrom)
                        if curr_local_rg4_motif_obj.flank_end != curr_local_rg4_motif_obj.end:
                            curr_right_flank_interval_list = self.merge_positions(position_array,
                                                                                  strand_array,
                                                                                  curr_local_rg4_motif_obj.end,
                                                                                  curr_local_rg4_motif_obj.flank_end,
                                                                                  _chrom)

                        per_cluster_rg4_motifs[id_tuple].append(
                            Global_rg4_motif_object(
                                interval_list=self.merge_positions(position_array,
                                                                   strand_array,
                                                                   curr_local_rg4_motif_obj.start,
                                                                   curr_local_rg4_motif_obj.end,
                                                                   _chrom),
                                rg4_type_id=curr_local_rg4_motif_obj.priority,
                                sequence=''.join(sequence[
                                                 curr_local_rg4_motif_obj.start:curr_local_rg4_motif_obj.end]),
                                cluster_id_tuple_list=[id_tuple],
                                left_flank_interval_list=curr_left_flank_interval_list,
                                right_flank_interval_list=curr_right_flank_interval_list,

                                flank_sequence=(
                                        sequence[curr_local_rg4_motif_obj.flank_start:curr_local_rg4_motif_obj.start] +
                                        '|' +
                                        sequence[curr_local_rg4_motif_obj.start:curr_local_rg4_motif_obj.end] +
                                        '|' +
                                        sequence[curr_local_rg4_motif_obj.end:curr_local_rg4_motif_obj.flank_end]
                                )

                            )
                        )
                        processed_local_motifs.add(curr_local_rg4_motif_obj.motif_id)

        final_cluster_groups = []
        rg4_cluster_by_id_tuple = {}
        child_rg4_motif_by_cluster_group = {}
        if len(_transcript_ids) > 1:  # Merge clusters
            pairwise_overlaps = {}
            completed_tids = set()
            # First, acquire pairwise overlapping relationships
            for self_tid in _transcript_ids:
                for self_rg4_cluster_obj in per_transcript_rg4_clusters[self_tid]:
                    rg4_cluster_by_id_tuple[(self_tid, self_rg4_cluster_obj.cluster_id)] = self_rg4_cluster_obj
                    for query_tid in _transcript_ids:
                        if query_tid in completed_tids:
                            continue
                        for query_self_rg4_cluster_obj in per_transcript_rg4_clusters[query_tid]:
                            have_overlap = False
                            for s_interval in self_rg4_cluster_obj.interval_list:
                                if have_overlap:
                                    break
                                for q_interval in query_self_rg4_cluster_obj.interval_list:
                                    if s_interval.overlap_with(q_interval):
                                        have_overlap = True
                                        break
                            if have_overlap:
                                self_id_tuple = (self_rg4_cluster_obj.parent_tid, self_rg4_cluster_obj.cluster_id)
                                query_id_tuple = (
                                    query_self_rg4_cluster_obj.parent_tid, query_self_rg4_cluster_obj.cluster_id)
                                if self_id_tuple not in pairwise_overlaps:
                                    pairwise_overlaps[self_id_tuple] = set()
                                if query_id_tuple not in pairwise_overlaps:
                                    pairwise_overlaps[query_id_tuple] = set()
                                pairwise_overlaps[self_id_tuple].add(query_id_tuple)
                                pairwise_overlaps[query_id_tuple].add(self_id_tuple)
                completed_tids.add(self_tid)

            # Expand pairwise to overall relationships
            while True:
                further_merging_occured = False
                for self_id_tuple in pairwise_overlaps.keys():
                    for x_member in pairwise_overlaps[self_id_tuple]:
                        for y_member in pairwise_overlaps[self_id_tuple]:
                            if x_member != y_member:
                                if x_member not in pairwise_overlaps[y_member]:
                                    further_merging_occured = True
                                    pairwise_overlaps[x_member].add(y_member)
                                    pairwise_overlaps[y_member].add(x_member)
                if not further_merging_occured:
                    break

            # Construct groups of rg4 clusters
            processed_clusters = set()
            for self_tid in _transcript_ids:
                for self_rg4_cluster_obj in per_transcript_rg4_clusters[self_tid]:
                    self_id_tuple = (self_rg4_cluster_obj.parent_tid, self_rg4_cluster_obj.cluster_id)
                    if self_id_tuple in processed_clusters:
                        continue
                    if self_id_tuple not in pairwise_overlaps:
                        final_cluster_groups.append(self_id_tuple)
                        processed_clusters.add(self_id_tuple)
                    else:
                        final_cluster_groups.append(tuple(sorted(pairwise_overlaps[self_id_tuple])))
                        processed_clusters |= pairwise_overlaps[self_id_tuple]
            del processed_clusters

            # For each rG4 cluster groups, merge the child rg4 motifs whenever possible
            for cluster_group_id, cluster_list in enumerate(final_cluster_groups):
                child_rg4_motif_by_cluster_group[cluster_group_id] = []
                child_rg4_motif_list = child_rg4_motif_by_cluster_group[cluster_group_id]

                child_rg4_motif_list_size = 0
                for cluster_id_tuple in cluster_list:
                    for self_rg4_motif_obj in per_cluster_rg4_motifs[cluster_id_tuple]:
                        if child_rg4_motif_list_size == 0:
                            child_rg4_motif_list.append(self_rg4_motif_obj)
                        else:
                            same_motif = False
                            for query_rg4_motif_obj in child_rg4_motif_list[0:child_rg4_motif_list_size]:
                                if (len(query_rg4_motif_obj.interval_list) ==
                                        len(self_rg4_motif_obj.interval_list) and
                                        len(query_rg4_motif_obj.left_flank_interval_list) ==
                                        len(self_rg4_motif_obj.left_flank_interval_list) and
                                        len(query_rg4_motif_obj.right_flank_interval_list) ==
                                        len(self_rg4_motif_obj.right_flank_interval_list)):
                                    same_motif = True
                                    if same_motif:
                                        for pos, query_interval in enumerate(query_rg4_motif_obj.interval_list):
                                            if not query_interval.is_same(self_rg4_motif_obj.interval_list[pos]):
                                                same_motif = False
                                                break
                                    if same_motif:
                                        for pos, query_interval in enumerate(
                                                query_rg4_motif_obj.left_flank_interval_list):
                                            if not query_interval.is_same(
                                                    self_rg4_motif_obj.left_flank_interval_list[pos]):
                                                same_motif = False
                                                break
                                    if same_motif:
                                        for pos, query_interval in enumerate(
                                                query_rg4_motif_obj.right_flank_interval_list):
                                            if not query_interval.is_same(
                                                    self_rg4_motif_obj.right_flank_interval_list[pos]):
                                                same_motif = False
                                                break
                                    if same_motif:
                                        assert query_rg4_motif_obj.sequence == self_rg4_motif_obj.sequence
                                        query_rg4_motif_obj.cluster_id_tuple_list.append(cluster_id_tuple)
                                        break
                                if same_motif:
                                    break
                            if not same_motif:
                                child_rg4_motif_list.append(self_rg4_motif_obj)
                        child_rg4_motif_list_size = len(child_rg4_motif_list)

        else:  # 1 transcript then doesnt need to merge
            for self_tid in _transcript_ids:
                for self_rg4_cluster_obj in per_transcript_rg4_clusters[self_tid]:
                    self_id_tuple = (self_rg4_cluster_obj.parent_tid, self_rg4_cluster_obj.cluster_id)
                    rg4_cluster_by_id_tuple[self_id_tuple] = self_rg4_cluster_obj
                    final_cluster_groups.append(tuple([self_id_tuple]))
                    child_rg4_motif_by_cluster_group[len(final_cluster_groups) - 1] = []
                    child_rg4_motif_list = child_rg4_motif_by_cluster_group[len(final_cluster_groups) - 1]

                    for self_rg4_motif_obj in per_cluster_rg4_motifs[self_id_tuple]:
                        child_rg4_motif_list.append(self_rg4_motif_obj)

        # "Forcefully" assign overall genomic interval boundaries for each cluster group
        # (which ensure the interval bounds everything regardless of splicing and can sacrifice stranded-ness)

        for cluster_group_id, cluster_list in enumerate(final_cluster_groups):
            cluster_group_left_pos = -1
            cluster_group_right_pos = -1
            cluster_group_strand = ''
            cluster_transcript_information_accumulator = []

            for cluster_id_tuple in cluster_list:
                # Determine cluster group boundaries
                curr_rg4_cluster_obj = rg4_cluster_by_id_tuple[cluster_id_tuple]
                _chrom, leftmost_pos, rightmost_pos, overall_strand = self.get_interval_list_boundary(
                    curr_rg4_cluster_obj.interval_list)

                if cluster_group_left_pos == -1 and cluster_group_right_pos == -1:
                    cluster_group_left_pos = leftmost_pos
                    cluster_group_right_pos = rightmost_pos
                    cluster_group_strand = overall_strand
                else:
                    assert cluster_group_left_pos != -1 and cluster_group_right_pos != -1
                    if leftmost_pos < cluster_group_left_pos:
                        cluster_group_left_pos = leftmost_pos
                    if rightmost_pos > cluster_group_right_pos:
                        cluster_group_right_pos = rightmost_pos
                    if cluster_group_strand != overall_strand:
                        cluster_group_strand = '.'

                # Add per-transcript information in form of transcript_id|cluster_id|UTRCDS|interval,interval;
                cluster_transcript_information_accumulator.append(
                    '|'.join(
                        [
                            str(curr_rg4_cluster_obj.parent_tid),
                            # str(curr_rg4_cluster_obj.cluster_id),
                            str(','.join(
                                _gene_obj.transcript_dict[curr_rg4_cluster_obj.parent_tid].get_utrcds_annotation(
                                    curr_rg4_cluster_obj.interval_list)
                            )),
                            ','.join([str(x.return_coordinates()) for x in curr_rg4_cluster_obj.interval_list]),
                        ]
                    )
                )

            cluster_bedline_accumulator = [_chrom,
                                           cluster_group_left_pos,
                                           cluster_group_right_pos + 1,
                                           _gene_obj.gene_name + '|' + _gene_obj.gene_id,
                                           cluster_group_id,
                                           cluster_group_strand]
            self.cluster_output[cluster_group_strand].append(
                cluster_bedline_accumulator + [';'.join(cluster_transcript_information_accumulator)])

            # In comparison the motif are handled with slightly more courtesy......
            for child_rg4_motif_obj in child_rg4_motif_by_cluster_group[cluster_group_id]:
                _chrom, _left, _right, _strand = self.get_interval_list_boundary(
                    child_rg4_motif_obj.interval_list)
                motif_bedline_accumulator = [_chrom,
                                             _left,
                                             _right + 1,
                                             _gene_obj.gene_name + '|' + _gene_obj.gene_id,
                                             cluster_group_strand,
                                             self.get_rg4_type_name(child_rg4_motif_obj.rg4_type_id),
                                             child_rg4_motif_obj.sequence,
                                             ','.join([str(x.return_coordinates()) for x in
                                                       child_rg4_motif_obj.interval_list]),
                                             child_rg4_motif_obj.flank_sequence,
                                             ','.join([str(x.return_coordinates()) for x in
                                                       child_rg4_motif_obj.left_flank_interval_list]),
                                             ','.join([str(x.return_coordinates()) for x in
                                                       child_rg4_motif_obj.right_flank_interval_list]),
                                             ]
                motif_transcript_information_accumulator = []
                for cluster_id_tuple in child_rg4_motif_obj.cluster_id_tuple_list:
                    parent_transcript_id = cluster_id_tuple[0]
                    motif_transcript_information_accumulator.append(
                        '|'.join(
                            [
                                str(parent_transcript_id),
                                str(','.join(
                                    _gene_obj.transcript_dict[parent_transcript_id].get_utrcds_annotation(
                                        child_rg4_motif_obj.interval_list))),
                            ]
                        )
                    )
                self.motif_output[cluster_group_strand].append(motif_bedline_accumulator[:-2] + [
                    ';'.join(motif_transcript_information_accumulator)] + motif_bedline_accumulator[-2:])

    def _export_sorted_motif_bed(self):
        with xopen(self.output_temp_motif, 'wb') as fw:
            for strand in ('+', '-'):
                for next_line in sorted(self.motif_output[strand]):
                    fw.write(('\t'.join([str(x) for x in next_line]) + '\n').encode('utf-8'))

    def run(self):
        logging.debug('Start rG4 finding in chromosome {0}'.format(self.chrom))
        fasta_genome = FastaGenome(self.fasta_path)
        input_transcriptome = gff_parser.Transcriptome()
        input_transcriptome.load_annotation(self.input_gff, chromosome_restriction=self.chrom)
        input_transcriptome.finalize()
        start = time.time()
        if self.chrom in input_transcriptome.chromosome_dict and len(
                input_transcriptome.chromosome_dict[self.chrom]) > 0:
            logging.debug(
                'Chromosome {0} has {1} genes'.format(self.chrom, len(input_transcriptome.chromosome_dict[self.chrom])))
        else:
            logging.debug('Chromosome {0} has no annotated genes'.format(self.chrom))
            return None

        rg4_patterns = self.rg4_regular_expression_generate()
        for i, gene_obj in enumerate(input_transcriptome.chromosome_dict[self.chrom]):
            all_transcript_ids = list(gene_obj.transcript_dict.keys())
            all_transcript_items = {}
            for transcript_id, transcript_obj in gene_obj.transcript_dict.items():
                if transcript_obj.junction_list:  # Have junctions, i.e. exons properly connected
                    full_transcript_sequence = ''
                    full_transcript_exons = []
                    for exon_obj in transcript_obj.exon_list:
                        full_transcript_sequence += fasta_genome.retrieve_sequence(exon_obj.interval)
                        full_transcript_exons.append(exon_obj)
                    all_transcript_items[transcript_id] = tuple(
                        [tuple([full_transcript_sequence, tuple(full_transcript_exons)])])
                else:
                    temp_accumulator = []
                    for exon_obj in transcript_obj.exon_list:
                        temp_accumulator.append(
                            tuple([fasta_genome.retrieve_sequence(exon_obj.interval), tuple([exon_obj])]))
                    all_transcript_items[transcript_id] = tuple(temp_accumulator)
            self.gene_find_rg4(gene_obj, all_transcript_ids, all_transcript_items, rg4_patterns, self.chrom)
        self._export_sorted_motif_bed()
        logging.debug('Completed rG4 finding chromosome {0} in {1}s'.format(self.chrom, time.time() - start))
        return self.output_temp_motif


def execute(input_annotation_gff_path, input_genome_fasta_path, prefix, processes, overwrite):
    whitelist_bed = whitelisting_driver(input_annotation_gff_path, input_genome_fasta_path, prefix, processes,
                                        overwrite)
    rg4_motif_bed = rg4_finding_driver(input_annotation_gff_path, input_genome_fasta_path, prefix, processes, overwrite)
    return whitelist_bed, rg4_motif_bed


if __name__ == "__main__":
    exit(0)
