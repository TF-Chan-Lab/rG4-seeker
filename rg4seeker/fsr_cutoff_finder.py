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

from collections import namedtuple
import HTSeq
import numpy as np
import warnings
import logging
import time
import os
# rg4seeker imports
from rg4seeker import apriori_rts_analysis

warnings.simplefilter('ignore', np.RankWarning)  # Turn off polyfit Rankwarning

FSRTSV_entry = namedtuple('fsrtsv_entry', ['Interval', 'treatment_coverage', 'treatment_readstart', 'control_coverage',
                                           'control_readstart', 'first_stage_alpha', 'second_stage_alpha_list'])
FIRST_STAGE_ALPHA_DICT = {'0.2': 0, '0.15874': 1, '0.12599': 2, '0.1': 3, '0.07937': 4, '0.062996': 5, '0.05': 6}
COVERAGE_MIN_FILTER = 6
FIRST_STAGE_ALPHA_ALLOWANCE = 0.2
FIRST_STAGE_MIN_FSR_DIFF = 0.05
FIRST_STAGE_ALPHA_VALUES = tuple([0.2, 0.15874, 0.12599, 0.1, 0.07937, 0.062996, 0.05])
SECOND_STAGE_EXPECTED_READSTOP_PROB = 0.025

Cutoff_point = namedtuple('fpr_entry',
                          ['first_alpha_cutoff', 'second_alpha_cutoff', 'not_unknown_count', 'unknown_count', 'FDR'])


def single_sample_load_fsrtsv(forward_in_fsrtsv, reverse_in_fsrtsv, _whitelist_array, _faidx_genome):
    unfiltered_treatment_fsr_store = {'+': [], '-': []}
    treatment_fsr_store = {'+': [], '-': []}
    unfiltered_decoy_fsr_store = {'+': [], '-': []}
    decoy_fsr_store = {'+': [], '-': []}

    apriori_rts_analysis.load_fsrtsv(treatment_fsr_store['+'], forward_in_fsrtsv, "None", _whitelist_array,
                                     _faidx_genome, disable_merging=True, disable_whitelisting=False,
                                     load_decoy_instead=False)
    apriori_rts_analysis.load_fsrtsv(treatment_fsr_store['-'], reverse_in_fsrtsv, "None", _whitelist_array,
                                     _faidx_genome, disable_merging=True, disable_whitelisting=False,
                                     load_decoy_instead=False)
    return unfiltered_treatment_fsr_store, treatment_fsr_store, unfiltered_decoy_fsr_store, decoy_fsr_store


def unknown_fdr_cut(_treatment_fsr, _rg4_motif_information, _rg4_threeprime_array, _criteria_structural_classes,
                    _faidx_genome, _target_fdr):
    def calculate_cutoff_points(not_unknown_list, unknown_list, first_alpha_cutoff):
        if not not_unknown_list:
            return None

        _result_cutoff_points = []
        not_unknown_pointer = 0
        unknown_pointer = 0
        _possible_cutoffs = sorted(list(set(unknown_list)), reverse=True) # Smaller cutoff value = more strict
        for second_alpha_cutoff in _possible_cutoffs:
            while True:
                if not_unknown_list[not_unknown_pointer] > second_alpha_cutoff:
                    not_unknown_pointer += 1
                else:
                    break
            not_unknown_count = len(not_unknown_list) - not_unknown_pointer

            if unknown_list:
                while True:
                    if unknown_list[unknown_pointer] > second_alpha_cutoff:
                        unknown_pointer += 1
                    else:
                        break
                unknown_count = len(unknown_list) - unknown_pointer
            else:
                unknown_count = 0

            fdr = unknown_count / (not_unknown_count + unknown_count)
            _result_cutoff_points.append(
                Cutoff_point(first_alpha_cutoff, second_alpha_cutoff, not_unknown_count, unknown_count, fdr))
        return _result_cutoff_points

    def find_best_cutoff(_result_cutoff_points, _target_fdr):
        if not _result_cutoff_points:
            return None
        best_cutoff_point = None
        best_not_unknown_count = 0
        if result_cutoff_points[0].FDR <= _target_fdr:
            return [result_cutoff_points[0]]
        for cutoff_point in _result_cutoff_points:
            if cutoff_point.FDR <= _target_fdr and cutoff_point.not_unknown_count > best_not_unknown_count:
                best_cutoff_point = cutoff_point
                best_not_unknown_count = cutoff_point.not_unknown_count
        if not best_cutoff_point:
            return [sorted(_result_cutoff_points, key=lambda x: x.FDR)[0]]
        else:
            return [best_cutoff_point]

    def find_global_cutoff(_result_cutoff_points, _target_fdr):
        if not _result_cutoff_points:
            return None
        best_cutoff_point = None
        best_not_unknown_count = 0
        for cutoff_point in _result_cutoff_points:
            if cutoff_point.FDR <= _target_fdr and cutoff_point.not_unknown_count > best_not_unknown_count:
                best_cutoff_point = cutoff_point
                best_not_unknown_count = cutoff_point.not_unknown_count
        if not best_cutoff_point:
            return sorted(_result_cutoff_points, key=lambda x: x.unknown_count)[0]
        else:
            return best_cutoff_point

    def load_alpha_values(_fsr_store, _positive_alpha_accu, _negative_alpha_accu, _rg4_motif_information,
                          _rg4_threeprime_array, _criteria_structural_classes, _faidx_genome):
        for strand_symbol in ['+', '-']:
            for item in _fsr_store[strand_symbol]:
                iv = item[0]
                curr_first_stage_alpha = item[2]
                curr_second_stage_alpha_list = item[3]
                if curr_first_stage_alpha and curr_second_stage_alpha_list:
                    status = 0 # 2: struct motifs, 1: G>=40%, 0: Others
                    for __, valset in _rg4_threeprime_array[iv].steps():
                        if valset:
                            for motif_serial in valset:
                                if _rg4_motif_information[motif_serial][5] in _criteria_structural_classes:
                                    status = 2

                    # Treat G>=40% as not_unknown as well
                    if status != 2:
                        verdict, __, __ = apriori_rts_analysis.check_novel_rg4(iv, _faidx_genome)
                        if verdict != 'Unknown':
                        #if 'potential' in verdict:
                            status = 1

                    if status == 0:
                        for k, val in enumerate(curr_second_stage_alpha_list):
                            _negative_alpha_accu[k].append(val)
                    if status == 2:
                        for k, val in enumerate(curr_second_stage_alpha_list):
                            _positive_alpha_accu[k].append(val)
        for k in _positive_alpha_accu.keys():
            _positive_alpha_accu[k] = sorted(_positive_alpha_accu[k], reverse=True) # Smaller value = more strict
        for k in _negative_alpha_accu.keys():
            _negative_alpha_accu[k] = sorted(_negative_alpha_accu[k], reverse=True) # Smaller value = more strict

    positive_alpha_accu = {k: [] for k in range(len(FIRST_STAGE_ALPHA_VALUES))}
    negative_alpha_accu = {k: [] for k in range(len(FIRST_STAGE_ALPHA_VALUES))}
    load_alpha_values(_treatment_fsr, positive_alpha_accu, negative_alpha_accu, _rg4_motif_information,
                      _rg4_threeprime_array, _criteria_structural_classes, _faidx_genome)

    all_possible_cutoffs = []
    all_optimal_cutoffs = []
    error_dump = []
    for k in range(len(FIRST_STAGE_ALPHA_VALUES)):
        result_cutoff_points = calculate_cutoff_points(positive_alpha_accu[k], negative_alpha_accu[k],
                                                       FIRST_STAGE_ALPHA_VALUES[k])
        if not result_cutoff_points:
            all_possible_cutoffs.append('No cutoff points available at {0}'.format(FIRST_STAGE_ALPHA_VALUES[k]))
        else:
            all_possible_cutoffs += [str(x) for x in result_cutoff_points]
        current_optimal_cutoffs = find_best_cutoff(result_cutoff_points, _target_fdr)
        if current_optimal_cutoffs:
            all_optimal_cutoffs += current_optimal_cutoffs
        else:
            error_dump.append('No optimal cutoff points available at {0}'.format(FIRST_STAGE_ALPHA_VALUES[k]))
    logging.debug(all_optimal_cutoffs)
    global_best_cutoff = find_global_cutoff(all_optimal_cutoffs, _target_fdr)

    return global_best_cutoff, error_dump + all_optimal_cutoffs + [''] + all_possible_cutoffs


class FSRCutoffFinder:
    def __init__(self, prefix_list, genome_fasta, whitelist_bed_list, motif_bed_list, bootstrap_no=0, debug_dump=False,
                 overwrite=True, target_fdr=0.015):
        self.prefix_list = prefix_list
        self.debug_dump = debug_dump
        self.io_tuples = []
        for prefix in prefix_list:
            forward_in_fsrtsv = prefix + '.bootstrap_{0:02d}.fsrtsv.{1}.gz'.format(bootstrap_no, 'forward')
            reverse_in_fsrtsv = prefix + '.bootstrap_{0:02d}.fsrtsv.{1}.gz'.format(bootstrap_no, 'reverse')
            output_cutoff_config = prefix + '.bootstrap_{0:02d}.cutoff.conf'.format(bootstrap_no)
            if self.debug_dump:
                output_cutoff_dump = prefix + '.bootstrap_{0:02d}.cutoff.dump'.format(bootstrap_no)
            else:
                output_cutoff_dump = None
            self.io_tuples.append([forward_in_fsrtsv, reverse_in_fsrtsv, output_cutoff_config, output_cutoff_dump])
        self.genome_fasta = genome_fasta
        self.motif_bed_list = motif_bed_list
        self.whitelist_bed_list = whitelist_bed_list
        self.overwrite = overwrite
        self.target_fdr = target_fdr

    def run(self):
        logging.info('[{0}] Target cutoff at FDR={1}'.format(time.ctime(), self.target_fdr))
        faidx_genome = apriori_rts_analysis.FaidxGenome(self.genome_fasta)
        whitelist_array = HTSeq.GenomicArray(chroms='auto', stranded=True, storage='step', typecode='b')
        for in_whitelist_bed in self.whitelist_bed_list:
            apriori_rts_analysis.load_whitelist(whitelist_array, in_whitelist_bed)

        rg4_motif_information, rg4_motif_array, rg4_threeprime_array = apriori_rts_analysis.load_rg4_motif(
            self.motif_bed_list, _all_rts_array=None)

        structural_classes = ['canonical/G3L1-7', 'longloop', 'bulges', 'two-quartet']
        for forward_in_fsrtsv, reverse_in_fsrtsv, output_cutoff_config, output_cutoff_dump in self.io_tuples:
            if not self.overwrite and os.path.exists(output_cutoff_config) and os.path.isfile(output_cutoff_config):
                continue
            u_treatment_fsr, treatment_fsr, u_decoy_fsr, decoy_fsr = single_sample_load_fsrtsv(forward_in_fsrtsv,
                                                                                               reverse_in_fsrtsv,
                                                                                               whitelist_array,
                                                                                               faidx_genome)
            rts_array = HTSeq.GenomicArray(chroms='auto', stranded=True, typecode='b')
            for strand_symbol in ('+', '-'):
                for fsr_store in [u_treatment_fsr]:
                    for items in fsr_store[strand_symbol]:
                        iv = items[0]
                        rts_array[iv] = True
            result, result_dump = unknown_fdr_cut(treatment_fsr, rg4_motif_information, rg4_threeprime_array,
                                                  frozenset(structural_classes), faidx_genome, self.target_fdr)
            assert result
            with open(output_cutoff_config, 'w') as fw:
                fw.write(str(result.first_alpha_cutoff) + '\n')
                fw.write(str(result.second_alpha_cutoff) + '\n')

            if self.debug_dump:
                with open(output_cutoff_dump, 'w') as fw:
                    for line in result_dump:
                        fw.write(str(line) + '\n')

            del u_treatment_fsr, treatment_fsr, u_decoy_fsr, decoy_fsr
            del rts_array
            del result, result_dump


if __name__ == "__main__":
    exit(0)
