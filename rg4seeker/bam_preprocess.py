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
**************************************************************************
"""

import pysam
import sys
import array
import random
from xopen import xopen
import gc
import re
import logging
import time

RANDOM_SEED = 2018
MD_TAG_DIGEST = re.compile('([0-9]+|[^]?[A-Za-z])')


class BamPreprocess:
    def __init__(self, in_bam, out_prefix, aligner, reads_type, no_of_bootstraps=0, ):
        self.in_bam = in_bam
        self.out_prefix = out_prefix
        self.no_of_bootstraps = no_of_bootstraps
        self.total_no_of_processed_reads = 0
        self.no_of_reads_by_chromosome = []
        self.all_line_id_to_read_id = {}
        self.aligner = aligner
        self.reads_type = reads_type
        return

    def run(self):
        logging.info('[{0}] Begin preprocessing BAM file {1}'.format(time.ctime(), self.in_bam))
        self._read_processor()
        self._random_generate()
        with open(self.out_prefix + '.chromosomes.txt', 'wb') as fw:
            for chrom in self.get_chromosomes():
                fw.write((chrom + '\n').encode('UTF-8'))
        logging.info('[{0}] Completed preprocessing BAM file {1}'.format(time.ctime(), self.in_bam))

    def get_chromosomes(self):
        return [x[0] for x in self.no_of_reads_by_chromosome]

    def _flusher(self, _chrom_index, _read_info):
        """
        Flush reads from 1 chromosome in buffer into output file.
        Returns the number of flushed reads.

        keyword arguments:
        fw : writable gzip file object, of the output file
        _chrom_index : dict, key=int, chromosome index ; value=string, chromosome name
        _read_info : list, of reads to flush out
        total_no_of_processed_reads : int, no of flushed reads before, used as an offset for read id

        """

        assert type(_chrom_index) is dict
        assert type(_read_info) is list

        curr_chrom_n = _read_info[0][0][0]  # Getting chromosome serial from first read
        curr_chrom = _chrom_index[curr_chrom_n]
        line_id_to_read_id = []

        def strand_convert(strand):
            if strand:
                return '+'
            else:
                return '-'

        logging.debug("[{0}] Flushing {1} reads in Chromosome {2}...".format(time.ctime(), len(_read_info), curr_chrom))
        sys.stdout.flush()

        # Transform _read_info into a list with read_id index for sorting
        iv_list = []
        for i in range(len(_read_info)):
            for subread_n, item in enumerate(_read_info[i]):
                assert item[0] == curr_chrom_n  # Verify chromosome index of all reads are same
                if item[2]:  # Forward Strand
                    subread_order = subread_n + 1
                else:
                    subread_order = len(_read_info[i]) - subread_n
                iv_list.append((item[1][0], item[1][1], item[2], i + self.total_no_of_processed_reads, subread_order))

        # Sort the list by start,end position regardless of strand
        iv_list = sorted(iv_list, key=lambda x: (x[0], x[1]))

        # Flush the sorted list to file
        with xopen(self.out_prefix + '.{0}.bamtsv.gz'.format(curr_chrom), mode='wb') as fw:
            for i in range(len(iv_list)):
                line_id_to_read_id.append(iv_list[i][3])
                fw.write(
                    (("{0}:{1}-{2}:{3}:{4}".format(curr_chrom,
                                                   iv_list[i][0],
                                                   iv_list[i][1],
                                                   strand_convert(iv_list[i][2]), iv_list[i][4])) + '\n').encode(
                        'utf-8')
                )

        self.no_of_reads_by_chromosome.append((curr_chrom, len(_read_info)))
        self.total_no_of_processed_reads += len(_read_info)
        self.all_line_id_to_read_id[curr_chrom] = array.array('l', line_id_to_read_id)
        logging.debug("[{0}] Chromosome {1} OK".format(time.ctime(), curr_chrom))
        sys.stdout.flush()

        return

    def _read_processor(self):
        """
        Takes input of a coordinate-sorted input BAM file (except single-end reads only)

        Key arguments:
        in_bam : string, full path to input BAM file
        out_tsv : string, full path to output BAMTSV.gz file

        """

        # Cigar Operations
        # 0M 1I 2D 3N 4S 5H 6P 7= 8X 9B
        # cigar_refoffset = {0: 1, 1: 0, 2: 1, 3: 1, 4: 0, 5: 0}

        readname_index = {}  # Dictionary of {read_name_index:read_name}
        read_info = []  # List of reads
        chrom_index = {}  # Dictionary of (chrom_index:chrom_name)
        processed_chroms = set()  # Set of processed chrom_name

        read_n = 0  # Counts no. of reads processed
        chrom_n = 0  # Counts no. of chromosomes processed

        if self.aligner == 'STAR':
            unique_quality_const = 255
        else:
            raise Exception('Aligner not supported')
        if self.reads_type == 'SE':
            require_all_flags = 0
            exclude_any_flags = 3844
        elif self.reads_type == 'PE':
            require_all_flags = 67
            exclude_any_flags = 3980
        else:
            raise Exception('Reads type not supported')

        with pysam.AlignmentFile(self.in_bam, 'rb') as bamfile:
            try:
                if bamfile.header['HD']['SO'] != 'coordinate':
                    raise KeyError
            except KeyError:
                raise Exception('BAM file {0} does not appear to be coordinate-sorted.'.format(self.in_bam))

            for read in bamfile:

                # Require proper alignment flag
                if (read.flag & require_all_flags) != require_all_flags:
                    continue
                if (read.flag & exclude_any_flags) > 0:
                    continue

                # Uniquely mapping reads only
                if self.aligner == 'STAR':
                    if not read.mapping_quality == unique_quality_const:
                        if read.get_tag('NH') != 1:
                            continue

                # NO short clips or errors allowed at 5'!
                if not read.is_reverse:
                    cigar_operation = read.cigartuples[0][0]
                else:
                    cigar_operation = read.cigartuples[-1][0]
                if not cigar_operation == 0 or cigar_operation == 7:
                    continue
                md_string = MD_TAG_DIGEST.findall(read.get_tag('MD'))

                if not read.is_reverse:
                    fiveprime_md = int(md_string[0])
                else:
                    fiveprime_md = int(md_string[-1])
                if fiveprime_md == 0:
                    continue

                # Get Read Information
                readid = read.qname
                chrom = read.reference_name
                if not read.is_reverse:
                    strand = True
                else:
                    strand = False

                # Check for change in processing chromosome
                if chrom not in processed_chroms:
                    chrom_index[chrom_n] = chrom
                    processed_chroms.add(chrom)
                    chrom_n += 1
                    # Flush information to ram when chromosome changed
                    if read_info:
                        self._flusher(chrom_index, read_info)
                        read_info = []  # Clear read buffer
                        readname_index = {}  # Clear readname index
                        read_n = 0  # Reset read index
                        gc.collect()

                    logging.debug("[{0}] Processing Chromosome {1}".format(time.ctime(), chrom))
                    sys.stdout.flush()

                # Index reads by their input order, reset once read buffer is flushed
                if readid not in readname_index:
                    readname_index[readid] = read_n
                    read_info.append([])
                    read_n += 1
                curr_index = readname_index[readid]
                aligned_regions = read.get_blocks()

                # Store information
                for iv in aligned_regions:
                    read_info[curr_index].append([chrom_n - 1, iv, strand])

            self._flusher(chrom_index, read_info)
        return

    def _random_generate(self):
        """
        Generate sorted random number tables
        """
        for chrom, no_of_reads in self.no_of_reads_by_chromosome:
            with xopen(self.out_prefix + '.{0}.bootstrap_{1:02d}.randtable.gz'.format(chrom, 0), 'wb') as fw:
                for line_id, read_id in enumerate(self.all_line_id_to_read_id[chrom]):
                    fw.write("{0}\t{1}\n".format(line_id, 1).encode('utf-8'))

        random.seed(a=RANDOM_SEED, version=2)
        bootstrap_seeds = []
        for i in range(self.no_of_bootstraps):
            bootstrap_seeds.append(random.randrange(self.total_no_of_processed_reads))

        for i in range(self.no_of_bootstraps):
            logging.debug(
                "[{0}] Generating random number tables for bootstrap {0}, BAM file {1}".format(time.ctime(), i + 1,
                                                                                               self.in_bam))
            random.seed(a=bootstrap_seeds[i], version=2)
            count_table = array.array('l', [0] * self.total_no_of_processed_reads)
            for n in range(self.total_no_of_processed_reads):
                count_table[random.randrange(self.total_no_of_processed_reads)] += 1
            for chrom, no_of_reads in self.no_of_reads_by_chromosome:
                with xopen(self.out_prefix + '.{0}.bootstrap_{1:02d}.randtable.gz'.format(chrom, i + 1),
                           mode='wb') as fw:
                    for line_id, read_id in enumerate(self.all_line_id_to_read_id[chrom]):
                        if count_table[read_id]:
                            fw.write("{0}\t{1}\n".format(line_id, count_table[read_id]).encode('utf-8'))


if __name__ == "__main__":
    exit(0)
