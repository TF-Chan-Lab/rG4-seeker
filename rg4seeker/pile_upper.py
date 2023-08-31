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
**************************************************************************
"""

import multiprocessing
from xopen import xopen
import array
from collections import namedtuple
import re
import logging
import gzip
import io

BOOTSTRAPPER_RINGBUFFER_SIZE = 200
Bamtsv_split = re.compile('^([^:]+):([0-9]+)-([0-9]+):([+-]):([0-9]+)\n')
Bamtsv_iv = namedtuple('bamtsv_iv', ['chrom', 'start', 'end', 'strand', 'subread_order'])


class PileUpper:
    def __init__(self, prefix, chromosome_list, bootstrap_no=0):
        self.rb_thread_list = []
        self.concat_dict = {}
        for strand_symbol in ('+', '-'):
            if strand_symbol == '+':
                strand = 'forward'
            elif strand_symbol == '-':
                strand = 'reverse'
            else:
                raise Exception('Stand = {0} ???'.format(strand_symbol))
            concat_covrss = prefix + '.bootstrap_{0:02d}.covrss.{1}.gz'.format(bootstrap_no, strand)
            self.concat_dict[concat_covrss] = []
            for chrom in chromosome_list:
                self.rb_thread_list.append(RBThread(prefix, chrom, bootstrap_no, strand_symbol))
                self.concat_dict[concat_covrss].append(self.rb_thread_list[-1].out_covrss)

    def run(self, processes=1):
        with multiprocessing.Pool(processes) as pool:
            for i in self.rb_thread_list:
                pool.apply_async(i.run, args=())
            pool.close()
            pool.join()

    def concat(self):
        for concat_covrss, split_covrss_list in self.concat_dict.items():
            with xopen(concat_covrss, 'wb', compresslevel=9) as fw:
                for split_covrss in split_covrss_list:
                    with xopen(split_covrss, 'rb') as f:
                        for line in f:
                            fw.write(line)


class RBThread:
    def __init__(self, prefix, chrom, bootstrap_no, strand_symbol):
        if strand_symbol == '+':
            strand = 'forward'
        elif strand_symbol == '-':
            strand = 'reverse'
        else:
            raise Exception('Stand = {0} ???'.format(strand_symbol))
        self.strand_symbol = strand_symbol
        self.in_tsv = prefix + '.{0}.bamtsv.gz'.format(chrom)
        self.in_randtable = prefix + '.{0}.bootstrap_{1:02d}.randtable.gz'.format(chrom, bootstrap_no)
        self.out_covrss = prefix + '.{0}.bootstrap_{1:02d}.covrss.{2}.gz'.format(chrom, bootstrap_no, strand)

    @staticmethod
    def parse_bamtsv(line):
        digest = Bamtsv_split.match(line)
        if not digest:
            raise Exception('Incorrect BAMTSV format, {0}'.format(line))
        item = digest.groups()
        return Bamtsv_iv(item[0], int(item[1]), int(item[2]), item[3], int(item[4]))

    @staticmethod
    def parse_random_table(line):
        fields = line.strip().split('\t')
        if len(fields) != 2:
            raise Exception('Incorrect random table format')
        bamtsv_line_id, bamtsv_line_repeats = [int(x) for x in fields]
        return bamtsv_line_id, bamtsv_line_repeats

    def get_next_read(self):
        with gzip.open(self.in_tsv, 'rb') as gztsv, gzip.open(self.in_randtable, 'rb') as gzrand:
            ftsv = io.BufferedReader(gztsv)
            frand = io.BufferedReader(gzrand)
            try:
                curr_bamtsv = self.parse_bamtsv(ftsv.readline().decode('UTF-8'))
            except StopIteration:
                raise Exception('Empty BAMTSV files {0}'.format(self.in_tsv))
            bamtsv_index = 0
            for next_random_table_line in frand.readlines():
                query_bamtsv_line_id, query_bamtsv_line_repeats = self.parse_random_table(
                    next_random_table_line.decode('UTF-8'))
                while bamtsv_index < query_bamtsv_line_id:
                    curr_bamtsv = self.parse_bamtsv(ftsv.readline().decode('UTF-8'))
                    bamtsv_index += 1
                if curr_bamtsv.strand == self.strand_symbol:
                    yield (curr_bamtsv, query_bamtsv_line_repeats)

    def run(self):
        logging.debug('Executing FSR Precompute Thread, {0} strand: using {1} '.format(self.strand_symbol,
                                                                                       self.in_tsv, ))
        with xopen(self.out_covrss, 'wb') as fcovrssw:
            coverage_ring_buffer = array.array('l', [0] * BOOTSTRAPPER_RINGBUFFER_SIZE)
            readstart_ring_buffer = array.array('l', [0] * BOOTSTRAPPER_RINGBUFFER_SIZE)
            rb_start_index = 0
            rb_end_index = 0
            chr_name = ''
            chr_strand = ''
            chr_pos_start = 0
            chr_pos_end = -1
            processed_chromosomes = set()
            for next_bamtsv_iv, count in self.get_next_read():

                if (rb_end_index - rb_start_index) == 0 or (
                        next_bamtsv_iv.start >= chr_pos_end and chr_name == next_bamtsv_iv.chrom) or (
                        chr_name != next_bamtsv_iv.chrom and next_bamtsv_iv.chrom not in processed_chromosomes):
                    if (rb_end_index - rb_start_index) != 0:
                        self.flush_rb(coverage_ring_buffer,
                                      readstart_ring_buffer,
                                      rb_start_index,
                                      rb_end_index,
                                      chr_pos_start,
                                      chr_pos_end,
                                      chr_name,
                                      chr_strand,
                                      fcovrssw)
                    rb_start_index, rb_end_index, chr_pos_start, chr_pos_end = self.initialize_rb(next_bamtsv_iv,
                                                                                                  count,
                                                                                                  coverage_ring_buffer,
                                                                                                  readstart_ring_buffer)
                    chr_name = next_bamtsv_iv.chrom
                    processed_chromosomes.add(next_bamtsv_iv.chrom)
                    chr_strand = next_bamtsv_iv.strand
                    continue

                if next_bamtsv_iv.start >= chr_pos_start and chr_name == next_bamtsv_iv.chrom:
                    rb_start_index, rb_end_index, chr_pos_start, chr_pos_end = self.update_rb(next_bamtsv_iv,
                                                                                              count,
                                                                                              coverage_ring_buffer,
                                                                                              readstart_ring_buffer,
                                                                                              rb_start_index,
                                                                                              rb_end_index,
                                                                                              chr_pos_start,
                                                                                              chr_pos_end,
                                                                                              chr_name,
                                                                                              chr_strand,
                                                                                              fcovrssw)
                    continue
                raise Exception('Unsorted BAMTSV input')
            if chr_pos_end != -1:  # Somehow there will be no reads coming through - nothing to flush
                self.flush_rb(coverage_ring_buffer,
                              readstart_ring_buffer,
                              rb_start_index,
                              rb_end_index,
                              chr_pos_start,
                              chr_pos_end,
                              chr_name,
                              chr_strand, fcovrssw)
            return

    @staticmethod
    def initialize_rb(_next_bamtsv_iv, _count, _coverage_ring_buffer, _readstart_ring_buffer):
        _chr_pos_start = _next_bamtsv_iv.start
        _chr_pos_end = _next_bamtsv_iv.end
        _rb_end_index = 0
        for i in range(BOOTSTRAPPER_RINGBUFFER_SIZE):
            _coverage_ring_buffer[i] = 0

        for i in range(_next_bamtsv_iv.start, _next_bamtsv_iv.end):
            _coverage_ring_buffer[_rb_end_index] = _count
            _rb_end_index += 1

        for i in range(BOOTSTRAPPER_RINGBUFFER_SIZE):
            _readstart_ring_buffer[i] = 0

        if _next_bamtsv_iv.strand == '+':
            if _next_bamtsv_iv.subread_order == 1:
                _readstart_ring_buffer[0] = _count
        elif _next_bamtsv_iv.strand == '-':
            if _next_bamtsv_iv.subread_order == 1:
                assert _readstart_ring_buffer[_rb_end_index - 1] == 0
                _readstart_ring_buffer[_rb_end_index - 1] = _count
                if _readstart_ring_buffer[_rb_end_index - 1] > _coverage_ring_buffer[_rb_end_index - 1]:
                    raise Exception('RSS>COV at initialize_rb')
                    pass
        else:
            raise Exception('Unstranded BAMTSV entry')
        return 0, _rb_end_index, _chr_pos_start, _chr_pos_end

    @staticmethod
    def flush_first_position_to_file(_coverage_ring_buffer, _readstart_ring_buffer,
                                     _rb_start_index, _chr_pos_start, _fcovrssw,
                                     _chr_name, _chr_strand):
        coverage_value = _coverage_ring_buffer[_rb_start_index % BOOTSTRAPPER_RINGBUFFER_SIZE]
        _coverage_ring_buffer[_rb_start_index % BOOTSTRAPPER_RINGBUFFER_SIZE] = 0
        readstart_value = _readstart_ring_buffer[_rb_start_index % BOOTSTRAPPER_RINGBUFFER_SIZE]
        _readstart_ring_buffer[_rb_start_index % BOOTSTRAPPER_RINGBUFFER_SIZE] = 0
        if not coverage_value:
            raise Exception('Zero value in coverage ring buffer')
        _fcovrssw.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(_chr_name,
                                                                _chr_pos_start,
                                                                _chr_pos_start + 1,
                                                                coverage_value,
                                                                readstart_value,
                                                                _chr_strand).encode('UTF-8'))
        return _rb_start_index + 1, _chr_pos_start + 1

    def flush_rb(self, _coverage_ring_buffer, _readstart_ring_buffer, _rb_start_index, _rb_end_index, _chr_pos_start,
                 _chr_pos_end, _chr_name, _chr_strand, _fcovrssw):
        # Flush everything in ring buffer
        assert ((_rb_end_index - _rb_start_index) == (_chr_pos_end - _chr_pos_start))
        while (_rb_end_index - _rb_start_index) > 0:
            _rb_start_index, _chr_pos_start = self.flush_first_position_to_file(_coverage_ring_buffer,
                                                                                _readstart_ring_buffer,
                                                                                _rb_start_index,
                                                                                _chr_pos_start,
                                                                                _fcovrssw,
                                                                                _chr_name,
                                                                                _chr_strand)
            assert _rb_start_index <= _rb_end_index

        return

    def update_rb(self, _next_bamtsv_iv, _count, _coverage_ring_buffer, _readstart_ring_buffer, _rb_start_index,
                  _rb_end_index, _chr_pos_start, _chr_pos_end, _chr_name, _chr_strand, _fcovrssw, ):
        # Update the rb according to the next bamtsv iv
        roll_shift = _next_bamtsv_iv.start - _chr_pos_start
        # Flush editable leftmost values

        for i in range(roll_shift):
            _rb_start_index, _chr_pos_start = self.flush_first_position_to_file(_coverage_ring_buffer,
                                                                                _readstart_ring_buffer,
                                                                                _rb_start_index, _chr_pos_start,
                                                                                _fcovrssw, _chr_name,
                                                                                _chr_strand)

        assert _chr_pos_start == _next_bamtsv_iv.start and (_rb_end_index - _rb_start_index) == (
                _chr_pos_end - _chr_pos_start)

        # Update overlapping coverage values
        overlapping_size = min(_next_bamtsv_iv.end, _chr_pos_end) - _chr_pos_start
        for i in range(overlapping_size):
            _coverage_ring_buffer[(_rb_start_index + i) % BOOTSTRAPPER_RINGBUFFER_SIZE] += _count

        # Append new coverage values
        right_extension_size = _next_bamtsv_iv.end - _chr_pos_end
        if right_extension_size > 0:
            for i in range(right_extension_size):
                _coverage_ring_buffer[_rb_end_index % BOOTSTRAPPER_RINGBUFFER_SIZE] += _count
                _rb_end_index += 1
            _chr_pos_end = _next_bamtsv_iv.end

        # Add new readstart value
        if _next_bamtsv_iv.strand == '+':
            if _next_bamtsv_iv.subread_order == 1:
                _readstart_ring_buffer[_rb_start_index % BOOTSTRAPPER_RINGBUFFER_SIZE] += _count
                if _readstart_ring_buffer[_rb_start_index % BOOTSTRAPPER_RINGBUFFER_SIZE] > _coverage_ring_buffer[
                    _rb_start_index % BOOTSTRAPPER_RINGBUFFER_SIZE]:
                    raise Exception('RSS>COV at update_rb, +ve strand')

        elif _next_bamtsv_iv.strand == '-':
            if _next_bamtsv_iv.subread_order == 1:
                _readstart_ring_buffer[(_rb_start_index + (
                        _next_bamtsv_iv.end - _next_bamtsv_iv.start) - 1) % BOOTSTRAPPER_RINGBUFFER_SIZE] += _count
                if _readstart_ring_buffer[(_rb_start_index + (
                        _next_bamtsv_iv.end - _next_bamtsv_iv.start) - 1) % BOOTSTRAPPER_RINGBUFFER_SIZE] > \
                        _coverage_ring_buffer[(_rb_start_index + (
                                _next_bamtsv_iv.end - _next_bamtsv_iv.start) - 1) % BOOTSTRAPPER_RINGBUFFER_SIZE]:
                    raise Exception('RSS>COV at update_rb, -ve strand')

        else:
            raise Exception('Unstranded BAMTSV entry')

        assert (_rb_end_index - _rb_start_index) == (_chr_pos_end - _chr_pos_start)
        return _rb_start_index, _rb_end_index, _chr_pos_start, _chr_pos_end


if __name__ == "__main__":
    exit(0)
