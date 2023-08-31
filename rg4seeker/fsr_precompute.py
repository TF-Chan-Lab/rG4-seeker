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
import multiprocessing
from xopen import xopen
import os
import logging
import time
from statsmodels.stats.proportion import proportion_confint
COVRSS_position = namedtuple('genomic_position_with_value',
                             ['chrom', 'pos', 'strand', 'coverage_value', 'readstart_value'])
Chrom_info = namedtuple('chrom_info', ['chrom_set', 'chrom_list'])

COVERAGE_MIN_FILTER = 6
FIRST_STAGE_ALPHA_ALLOWANCE = 0.2
FIRST_STAGE_MIN_FSR_DIFF = 0.05
FIRST_STAGE_ALPHA_VALUES = tuple([0.2, 0.15874, 0.12599, 0.1, 0.07937, 0.062996, 0.05])
SECOND_STAGE_EXPECTED_READSTOP_PROB = 0.025
SECOND_STAGE_ITERATION_LIMIT = 16
SECOND_STAGE_ERROR_ALLOWANCE = 0.001

def first_stage(tmnt_rss, tmnt_cov, ctrl_rss, ctrl_cov):
    tmnt_fsr = tmnt_rss / tmnt_cov
    if tmnt_fsr < FIRST_STAGE_MIN_FSR_DIFF:
        return 0.0
    pass_alpha = 0.0
    for test_alpha in FIRST_STAGE_ALPHA_VALUES:
        t_ci_min, t_ci_max = proportion_confint(tmnt_rss, tmnt_cov, test_alpha, method='wilson')
        c_ci_min, c_ci_max = proportion_confint(ctrl_rss, ctrl_cov, test_alpha, method='wilson')
        delta_fsr = t_ci_min - c_ci_max
        if delta_fsr >= 0:
            pass_alpha = test_alpha
        else:
            break
    return pass_alpha


def second_stage(tmnt_rss, tmnt_cov, ctrl_rss, ctrl_cov, _first_alpha):
    def get_first_stage_offset(_ctrl_rss, _ctrl_cov, _tmnt_rss, _tmnt_cov, _alpha):
        t_ci_min, t_ci_max = proportion_confint(_tmnt_rss, _tmnt_cov, _alpha, method='wilson')
        c_ci_min, c_ci_max = proportion_confint(_ctrl_rss, _ctrl_cov, _alpha, method='wilson')
        return (t_ci_min - c_ci_max)

    def second_stage_function(_alpha, _remainder_fsr, _n):
        ci_min, ci_max = proportion_confint(n * SECOND_STAGE_EXPECTED_READSTOP_PROB, _n, _alpha, method='wilson')
        second_stage_offset = ci_max
        return _remainder_fsr - second_stage_offset

    def second_stage_optimize(_remainder_fsr, _tmnt_cov):
        if second_stage_function(FIRST_STAGE_ALPHA_VALUES[0], _remainder_fsr, _tmnt_cov) < 0:
            return None
        if second_stage_function(FIRST_STAGE_ALPHA_VALUES[-1], _remainder_fsr, _tmnt_cov) >= 0:
            return FIRST_STAGE_ALPHA_VALUES[-1]
        upper_limit = FIRST_STAGE_ALPHA_VALUES[0]
        lower_limit = FIRST_STAGE_ALPHA_VALUES[-1]
        iteration_count = 0
        while True:
            second_stage_test_alpha = (upper_limit + lower_limit) / 2
            result = second_stage_function(second_stage_test_alpha, _remainder_fsr, _tmnt_cov)
            iteration_count += 1
            if abs(result) < SECOND_STAGE_ERROR_ALLOWANCE or iteration_count > SECOND_STAGE_ITERATION_LIMIT:
                return second_stage_test_alpha
            if result > 0:
                upper_limit = second_stage_test_alpha
            else:
                lower_limit = second_stage_test_alpha

    n = tmnt_cov  # sample size is determined by treatment coverage
    output_second_alpha_list = []

    for first_stage_test_alpha in FIRST_STAGE_ALPHA_VALUES:
        remainder_fsr = get_first_stage_offset(ctrl_rss, ctrl_cov, tmnt_rss, tmnt_cov, first_stage_test_alpha)
        output_second_alpha_list.append(second_stage_optimize(remainder_fsr, tmnt_cov))
        if not output_second_alpha_list[-1]:
            output_second_alpha_list.pop()
            break
    if not output_second_alpha_list:
        return None
    return output_second_alpha_list


class FSRPrecompute:
    def __init__(self, treatment_prefix, control_prefix, chromosome_list, bootstrap_no=0):
        self.precompute_thread_list = []
        self.concat_dict = {}
        for strand_symbol in ('+', '-'):
            if strand_symbol == '+':
                strand = 'forward'
            elif strand_symbol == '-':
                strand = 'reverse'
            else:
                raise Exception('Stand = {0} ???'.format(strand_symbol))
            if 'Li-rep' in treatment_prefix:
                concat_fsrtsv = treatment_prefix + '-rep' + control_prefix[
                    -1] + '.bootstrap_{0:02d}.fsrtsv.{1}.gz'.format(bootstrap_no, strand)
            else:
                concat_fsrtsv = treatment_prefix + '.bootstrap_{0:02d}.fsrtsv.{1}.gz'.format(bootstrap_no, strand)
            self.concat_dict[concat_fsrtsv] = []
            for chrom in chromosome_list:
                self.precompute_thread_list.append(
                    PrecomputeThread(treatment_prefix, control_prefix, chrom, bootstrap_no, strand_symbol))
                self.concat_dict[concat_fsrtsv].append(self.precompute_thread_list[-1].out_fsrtsv)

    def run(self, processes=8):
        with multiprocessing.Pool(processes) as pool:
            multiple_results = [pool.apply_async(i.run, args=()) for i in self.precompute_thread_list]
            [res.get() for res in multiple_results]
        for concat_fsrtsv, split_fsrtsv_list in self.concat_dict.items():
            with xopen(concat_fsrtsv, 'wb', compresslevel=9) as fw:
                for split_fsrtsv in split_fsrtsv_list:
                    with xopen(split_fsrtsv, 'rb') as f:
                        for line in f:
                            fw.write(line)
                    os.remove(split_fsrtsv)


class PrecomputeThread:
    def __init__(self, treatment_prefix, control_prefix, chrom, bootstrap_no, strand_symbol):

        if strand_symbol == '+':
            strand = 'forward'
        elif strand_symbol == '-':
            strand = 'reverse'
        else:
            raise Exception('Stand = {0} ???'.format(strand_symbol))
        self.strand_symbol = strand_symbol
        self.in_treatment_covrss = treatment_prefix + '.{0}.bootstrap_{1:02d}.covrss.{2}.gz'.format(chrom, bootstrap_no,
                                                                                                    strand)
        self.in_control_covrss = control_prefix + '.{0}.bootstrap_{1:02d}.covrss.{2}.gz'.format(chrom, bootstrap_no,
                                                                                                strand)
        if 'Li-rep' in treatment_prefix:
            self.out_fsrtsv = treatment_prefix + '-rep' + control_prefix[
                -1] + '.{0}.bootstrap_{1:02d}.fsrtsv.{2}.gz'.format(chrom, bootstrap_no, strand)
        else:
            self.out_fsrtsv = treatment_prefix + '.{0}.bootstrap_{1:02d}.fsrtsv.{2}.gz'.format(chrom, bootstrap_no,
                                                                                               strand)

    def run(self):
        logging.debug(
            '{0} Executing FSR Precompute Thread, {1} strand: using {2},{3}'.format(time.ctime(), self.strand_symbol,
                                                                                    self.in_treatment_covrss,
                                                                                    self.in_control_covrss))
        with xopen(self.out_fsrtsv, 'wb', compresslevel=9) as fw:
            for line in self.alpha_compute():
                if line:
                    # print(line)
                    fw.write(('\t'.join(line) + '\n').encode('UTF-8'))
        return

    @staticmethod
    def covrss_generator(in_covrss):
        with xopen(in_covrss, 'rb') as f:
            for line in f:
                fields = line.decode('UTF-8').strip().split('\t')
                if len(fields) != 6:
                    raise Exception('Incorrect covrss tsv format')
                assert int(fields[2]) - int(fields[1]) == 1
                yield COVRSS_position(fields[0], int(fields[1]), fields[5], int(fields[3]), int(fields[4]))

    def covrss_muxer(self):
        ftreatment = self.covrss_generator(self.in_treatment_covrss)
        fcontrol = self.covrss_generator(self.in_control_covrss)
        while True:
            try:
                treatment_covrss_entry = ftreatment.__next__()
                control_covrss_entry = fcontrol.__next__()
                while True:
                    assert treatment_covrss_entry.chrom == control_covrss_entry.chrom

                    if treatment_covrss_entry.pos == control_covrss_entry.pos:
                        # Speed hack: check if there is something to test with
                        if (treatment_covrss_entry.readstart_value >= 1) or (control_covrss_entry.readstart_value >= 1):
                            yield (treatment_covrss_entry, control_covrss_entry)
                        treatment_covrss_entry = ftreatment.__next__()
                        control_covrss_entry = fcontrol.__next__()
                        continue
                    elif treatment_covrss_entry.pos < control_covrss_entry.pos:
                        treatment_covrss_entry = ftreatment.__next__()
                        continue
                    elif treatment_covrss_entry.pos > control_covrss_entry.pos:
                        control_covrss_entry = fcontrol.__next__()
                        continue
                    else:
                        raise Exception("Unable to compare positions of entries")
            except StopIteration:
                return

    def alpha_compute(self):
        for treatment_covrss, control_covrss in self.covrss_muxer():
            if (treatment_covrss.coverage_value < COVERAGE_MIN_FILTER or
                    control_covrss.coverage_value < COVERAGE_MIN_FILTER):
                """
                yield ([str(x) for x in [treatment_covrss.chrom,
                                         treatment_covrss.pos,
                                         treatment_covrss.pos + 1,
                                         treatment_covrss.strand,
                                         treatment_covrss.coverage_value,
                                         treatment_covrss.readstart_value,
                                         control_covrss.coverage_value,
                                         control_covrss.readstart_value, ]])
                """

                continue
            # print(treatment_covrss, control_covrss)
            treatment_fsr = treatment_covrss.readstart_value / treatment_covrss.coverage_value
            control_fsr = control_covrss.readstart_value / control_covrss.coverage_value

            first_alpha = None
            second_alpha_list = []
            if treatment_fsr > control_fsr:
                first_alpha = first_stage(tmnt_rss=treatment_covrss.readstart_value,
                                          tmnt_cov=treatment_covrss.coverage_value,
                                          ctrl_rss=control_covrss.readstart_value,
                                          ctrl_cov=control_covrss.coverage_value)

                if first_alpha <= FIRST_STAGE_ALPHA_ALLOWANCE:
                    second_alpha_list = second_stage(
                        tmnt_rss=treatment_covrss.readstart_value,
                        tmnt_cov=treatment_covrss.coverage_value,
                        ctrl_rss=control_covrss.readstart_value,
                        ctrl_cov=control_covrss.coverage_value,
                        _first_alpha=first_alpha)
            elif treatment_fsr < control_fsr:
                first_alpha = first_stage(ctrl_rss=treatment_covrss.readstart_value,
                                          ctrl_cov=treatment_covrss.coverage_value,
                                          tmnt_rss=control_covrss.readstart_value,
                                          tmnt_cov=control_covrss.coverage_value)

                if first_alpha <= FIRST_STAGE_ALPHA_ALLOWANCE:
                    second_alpha_list = second_stage(
                        ctrl_rss=treatment_covrss.readstart_value,
                        ctrl_cov=treatment_covrss.coverage_value,
                        tmnt_rss=control_covrss.readstart_value,
                        tmnt_cov=control_covrss.coverage_value,
                        _first_alpha=first_alpha)
            if first_alpha:
                yield ([str(x) for x in [treatment_covrss.chrom,
                                         treatment_covrss.pos,
                                         treatment_covrss.pos + 1,
                                         treatment_covrss.strand,
                                         treatment_covrss.coverage_value,
                                         treatment_covrss.readstart_value,
                                         control_covrss.coverage_value,
                                         control_covrss.readstart_value,
                                         first_alpha]] +
                       ([str(x) for x in second_alpha_list] if second_alpha_list else []))
        return


if __name__ == "__main__":
    exit(0)
