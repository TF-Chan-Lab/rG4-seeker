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

import HTSeq
import time
from xopen import xopen
import logging
from collections import OrderedDict, deque

EXON_END_BLACKLIST_RANGE = 10
SEQ_FILTER_SEARCH_RANGE = 10
REFSEQ_TRANSCRIPT_TYPE = {"antisense_RNA", "guide_RNA", "lnc_RNA", "mRNA", "primary_transcript", "RNase_MRP_RNA",
                          "RNase_P_RNA", "rRNA", "scRNA", "snoRNA", "snRNA", "telomerase_RNA", "transcript", "tRNA",
                          "vault_RNA", "Y_RNA", 'C_gene_segment', 'D_gene_segment', 'J_gene_segment', 'V_gene_segment'}


class Transcriptome:
    def __init__(self):
        self.chromosome_dict = OrderedDict()
        self.gene_dict = OrderedDict()
        self.transcript_parent_dict = {}
        self.editable = True
        self.transcript_obj_array = None

    def __getitem__(self, key):
        return self.gene_dict[key]

    def _add_gene(self, _gene_id, _gene_interval, _gene_type, _gene_name):
        if _gene_id in self.gene_dict:
            if self.gene_dict[_gene_id].interval.is_same(_gene_interval):
                logging.warning('Warning: Duplicated input gene entry with gene_id {0}'.format(_gene_id))
            else:
                logging.critical(
                    'Error: Duplicated input gene entry with gene_id {0} and conflicting interval'.format(_gene_id))
                raise Exception(
                    'Error: Duplicated input gene entry with gene_id {0} and conflicting interval'.format(_gene_id))
        if _gene_interval.chrom not in self.chromosome_dict:
            self.chromosome_dict[_gene_interval.chrom] = []

        self.gene_dict[_gene_id] = Gene(_gene_id, _gene_interval, _gene_type, _gene_name)
        self.chromosome_dict[_gene_interval.chrom].append(self.gene_dict[_gene_id])

    def _add_transcript(self, _gene_id, _transcript_id, _transcript_interval, _transcript_type):
        if _gene_id not in self.gene_dict:
            logging.debug('Gene ID ' + _gene_id + ' NOT FOUND')
            raise KeyError
        self.transcript_parent_dict[_transcript_id] = self.gene_dict[_gene_id]
        self.gene_dict[_gene_id].add_transcript(_transcript_id, _transcript_interval, _transcript_type)

    def _add_exon(self, _transcript_id, _exon_interval, _exon_number):
        if _transcript_id not in self.transcript_parent_dict:
            logging.debug('Transcript ID ' + _transcript_id + ' NOT FOUND')
            raise KeyError
        self.transcript_parent_dict[_transcript_id][_transcript_id].add_exon(_exon_interval, _exon_number)

    def _add_cds(self, _transcript_id, _feature_interval):
        if _transcript_id not in self.transcript_parent_dict:
            logging.debug('Transcript ID ' + _transcript_id + ' NOT FOUND')
            raise KeyError
        self.transcript_parent_dict[_transcript_id][_transcript_id].add_cds(_feature_interval)

    def _add_utr5(self, _transcript_id, _feature_interval):
        if _transcript_id not in self.transcript_parent_dict:
            logging.debug('Transcript ID ' + _transcript_id + ' NOT FOUND')
            raise KeyError
        self.transcript_parent_dict[_transcript_id][_transcript_id].add_utr5(_feature_interval)

    def _add_utr3(self, _transcript_id, _feature_interval):
        if _transcript_id not in self.transcript_parent_dict:
            logging.debug('Transcript ID ' + _transcript_id + ' NOT FOUND')
            raise KeyError
        self.transcript_parent_dict[_transcript_id][_transcript_id].add_utr3(_feature_interval)

    def _record_gene(self, _annotation_line):
        if 'ID' in _annotation_line.attr:
            gene_id = _annotation_line.attr['ID']
        elif 'gene_id' in _annotation_line.attr:
            gene_id = _annotation_line.attr['gene_id']
        else:
            raise Exception(
                'Unable to find corresponding gene_id attribute in gene annotation, {0}'.format(_annotation_line.attr))
        if 'gene_type' in _annotation_line.attr:
            gene_type = _annotation_line.attr['gene_type']
        else:
            gene_type = _annotation_line.type
        if 'gene_name' in _annotation_line.attr:
            gene_name = _annotation_line.attr['gene_name']
        elif 'gene' in _annotation_line.attr:
            gene_name = _annotation_line.attr['gene']
        elif 'Name' in _annotation_line.attr:
            gene_name = _annotation_line.attr['Name']
        elif 'description' in _annotation_line.attr:
            gene_name = _annotation_line.attr['description']
        else:
            gene_name = ''

        self._add_gene(gene_id, Interval(_annotation_line.iv.chrom, _annotation_line.iv.start,
                                         _annotation_line.iv.end,
                                         _annotation_line.iv.strand), gene_type, gene_name)
        if (_annotation_line.type == 'transposable_element_gene' or
                _annotation_line.type == 'blocked_reading_frame' or
                _annotation_line.type == 'pseudogene'):
            self._add_transcript(gene_id, gene_id,
                                 Interval(_annotation_line.iv.chrom, _annotation_line.iv.start,
                                          _annotation_line.iv.end,
                                          _annotation_line.iv.strand), gene_type)

    def _record_transcript(self, _annotation_line):
        transcript_interval = Interval(_annotation_line.iv.chrom, _annotation_line.iv.start,
                                       _annotation_line.iv.end,
                                       _annotation_line.iv.strand)
        if 'Parent' in _annotation_line.attr:
            gene_id = _annotation_line.attr['Parent']
        elif 'gene_id' in _annotation_line.attr:
            gene_id = _annotation_line.attr['gene_id']
        else:
            logging.critical('Unable to find corresponding gene_id attribute in transcript annotation, {0}'.format(
                _annotation_line.attr))
            raise Exception('Unable to find corresponding gene_id attribute in transcript annotation, {0}'.format(
                _annotation_line.attr))

        if 'ID' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['ID']
        elif 'transcript_id' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['transcript_id']
        else:
            _transcript_id = gene_id + ':transcript:' + transcript_interval.return_coordinates()
            logging.warning('Warning: Unable to find corresponding transcript_id attribute in transcript {0}'.format(
                _transcript_id))

        if 'transcript_type' in _annotation_line.attr:
            transcript_type = _annotation_line.attr['transcript_type']
        else:
            transcript_type = _annotation_line.type
        self._add_transcript(gene_id, _transcript_id,
                             transcript_interval,
                             transcript_type)

    def _record_exon(self, _annotation_line):
        if 'Parent' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['Parent']
        elif 'transcript_id' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['transcript_id']
        else:
            logging.critical(
                'Unable to find corresponding transcript_id attribute in exon annotation for exon ID {0}'.format(
                    _annotation_line.attr['ID']))
            raise Exception(
                'Unable to find corresponding transcript_id attribute in exon annotation for exon ID {0}'.format(
                    _annotation_line.attr['ID']))
        if 'exon_number' in _annotation_line.attr:
            exon_number = int(_annotation_line.attr['exon_number'])
        else:
            exon_number = None
        self._add_exon(_transcript_id, Interval(_annotation_line.iv.chrom, _annotation_line.iv.start,
                                                _annotation_line.iv.end,
                                                _annotation_line.iv.strand), exon_number)

    def _record_cds(self, _annotation_line):
        if 'Parent' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['Parent']
        elif 'transcript_id' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['transcript_id']
        else:
            raise Exception('Unable to find corresponding transcript_id attribute in CDS annotation')
        if ',' in _transcript_id and False:
            for sub_transcript_id in _transcript_id.split(','):
                self._add_cds(sub_transcript_id, Interval(_annotation_line.iv.chrom,
                                                          _annotation_line.iv.start,
                                                          _annotation_line.iv.end,
                                                          _annotation_line.iv.strand))
        else:
            self._add_cds(_transcript_id, Interval(_annotation_line.iv.chrom,
                                                   _annotation_line.iv.start,
                                                   _annotation_line.iv.end,
                                                   _annotation_line.iv.strand))
        return

    def _record_utr5(self, _annotation_line):
        if 'transcript_id' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['transcript_id']
        elif 'Parent' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['Parent']
        else:
            raise Exception('Unable to find corresponding transcript_id attribute in UTR5 annotation')
        self._add_utr5(_transcript_id, Interval(_annotation_line.iv.chrom,
                                                _annotation_line.iv.start,
                                                _annotation_line.iv.end,
                                                _annotation_line.iv.strand))
        return

    def _record_utr3(self, _annotation_line):
        if 'transcript_id' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['transcript_id']
        elif 'Parent' in _annotation_line.attr:
            _transcript_id = _annotation_line.attr['Parent']
        else:
            raise Exception('Unable to find corresponding transcript_id attribute in UTR3 annotation')
        self._add_utr3(_transcript_id, Interval(_annotation_line.iv.chrom,
                                                _annotation_line.iv.start,
                                                _annotation_line.iv.end,
                                                _annotation_line.iv.strand))
        return

    def _process_annotation_line(self, annotation_line, is_retry=False):
        if 'ID' in annotation_line.attr and 'PAR_Y' in annotation_line.attr['ID']:
            return
        if (
                annotation_line.type == 'gene' or
                annotation_line.type == 'pseudogene'
        ):
            self._record_gene(annotation_line)
        elif (
                annotation_line.type == 'transcript' or
                annotation_line.type in REFSEQ_TRANSCRIPT_TYPE
        ):
            self._record_transcript(annotation_line)
        elif annotation_line.type == 'exon':
            if is_retry and 'ID' in annotation_line.attr and 'MIR' in annotation_line.attr['ID']:
                return
            self._record_exon(annotation_line)
        elif annotation_line.type == 'CDS':
            if is_retry and 'Parent' in annotation_line.attr and 'gene' in annotation_line.attr['Parent']:
                return
            self._record_cds(annotation_line)
        return

    def load_annotation(self, _in_gff_path, chromosome_restriction=None):
        if not self.editable:
            logging.error("Error: Transcriptome already finalized.")
            return

        i = 0
        input_buffer = deque()
        start = time.time()
        with xopen(_in_gff_path) as f:
            rawlines = []
            if not chromosome_restriction:
                f = HTSeq.GFF_Reader(f)

            else:
                chrom_name = chromosome_restriction + '\t'
                chrom_name_length = len(chromosome_restriction) + 1

                for rawline in f:
                    if rawline[0:chrom_name_length] == chrom_name:
                        rawlines.append(rawline)
                f = HTSeq.GFF_Reader(rawlines)

            for line in f:
                i += 1
                if i % 100000 == 0:
                    if not chromosome_restriction:
                        logging.debug("Processed {0} GFF lines in {1:2f}s".format(i, time.time() - start))
                    else:
                        pass

                try:
                    self._process_annotation_line(line)
                except KeyError:
                    input_buffer.append(line)

            del rawlines
        initial_input_buffer_size = len(input_buffer)

        i = 0
        while True:
            if len(input_buffer) == 0:
                break
            line = input_buffer.pop()
            try:
                self._process_annotation_line(line, is_retry=True)
            except KeyError:
                print(line)
                input_buffer.append(line)
            i += 1
            if i % initial_input_buffer_size == 0:
                if len(input_buffer) == initial_input_buffer_size or i >= initial_input_buffer_size * 10:
                    for item in input_buffer:
                        print(item)
                    logging.critical('Cannot digest remaining {0} lines in annotation'.format(len(input_buffer)))
                    raise Exception('Cannot digest remaining {0} lines in annotation'.format(len(input_buffer)))
                else:
                    initial_input_buffer_size = len(input_buffer)
                    i = 0

        return

    def finalize(self):
        if not self.editable:
            logging.error("Error: Transcriptome already finalized.")
            return
        for gene_object in self.gene_dict.values():
            gene_object.finalize()
        self.editable = False

        return

    def _generate_transcript_object_array(self):
        self.transcript_obj_array = HTSeq.GenomicArrayOfSets(chroms='auto', stranded=True, storage='step')
        for gene_id, gene_obj in self.gene_dict.items():
            for transcript_id, transcript_obj in gene_obj.transcript_dict.items():
                self.transcript_obj_array[transcript_obj.interval.htseq_iv()] += transcript_obj

    def annotation_query(self, _query_interval_list):
        if not self.transcript_obj_array:
            self._generate_transcript_object_array()
        result_tobj = set()
        result_gene_dict = {}
        for query_iv in _query_interval_list:
            query_htseqiv = HTSeq.GenomicInterval(query_iv.chrom, query_iv.start, query_iv.end, query_iv.strand)
            for iv, tobj_set in self.transcript_obj_array[query_htseqiv].steps():
                result_tobj |= tobj_set

        for tobj in result_tobj:
            gene_obj = self.transcript_parent_dict[tobj.transcript_id]
            gene_id = gene_obj.gene_id
            gene_name = gene_obj.gene_name
            if 'gene' not in gene_name:
                full_gene_name = '{0}|{1}'.format(gene_name, gene_id)
            else:
                full_gene_name = gene_name
            if full_gene_name not in result_gene_dict:
                result_gene_dict[full_gene_name] = [[], set()]
            result_gene_dict[full_gene_name][0].append(tobj.transcript_id)
            result_gene_dict[full_gene_name][1] |= set(tobj.get_utrcds_annotation(_query_interval_list))

        if not result_gene_dict:
            return [None, None, None]
        else:
            all_gene_names = []
            transcript_ids = []
            coding_status = set()
            for full_gene_name in sorted(result_gene_dict.keys()):
                all_gene_names.append(full_gene_name)
                transcript_ids += sorted(list(result_gene_dict[full_gene_name][0]))
                cdsutr = result_gene_dict[full_gene_name][1]
                cdsutr &= {'5UTR', 'CDS', '3UTR'}
                if len(cdsutr) == 0:
                    cdsutr = {'non_coding'}
                coding_status |= cdsutr
            return [sorted(list(x)) if x else [] for x in [all_gene_names, transcript_ids, coding_status]]


class Gene:

    def __init__(self, _gene_id, _gene_interval, _gene_type, _gene_name):
        self.gene_id = _gene_id
        self.gene_type = _gene_type
        self.gene_name = _gene_name
        self.interval = _gene_interval
        self.transcript_dict = OrderedDict()

    def __getitem__(self, key):
        return self.transcript_dict[key]

    def add_transcript(self, _transcript_id, _transcript_interval, _transcript_type):
        if _transcript_id in self.transcript_dict:
            raise Exception('Duplicated input transcript')
        if _transcript_interval.chrom != self.interval.chrom:
            raise Exception('Gene/Transcript chromosome mismatch in {0}'.format(self.gene_id))
        self.transcript_dict[_transcript_id] = Transcript(_transcript_id, _transcript_interval, _transcript_type)

    def finalize(self):
        transcript_type_list = []
        for transcript_object in self.transcript_dict.values():
            transcript_type_list.append(transcript_object.finalize())
        if not self.gene_type:
            transcript_type_list = list(frozenset(transcript_type_list))
            if len(transcript_type_list) == 1:
                self.gene_type = transcript_type_list[0]
            elif 'protein_coding' in transcript_type_list:
                self.gene_type = 'protein_coding'
            elif 'mRNA' in transcript_type_list:
                self.gene_type = 'mRNA'
            else:
                self.gene_type = ','.join(sorted(transcript_type_list))


class Transcript:

    def __init__(self, _transcript_id, _transcript_interval, _transcript_type):
        self.transcript_id = _transcript_id
        self.transcript_type = _transcript_type
        self.interval = _transcript_interval
        self.utr5_list = []
        self.cds_list = []
        self.utr3_list = []
        self.exon_list = []
        self.junction_list = []
        self.editable = True

    def add_exon(self, _exon_interval, _exon_number):
        if not self.editable:
            raise Exception('Exon addition on editable transcript {0}'.format(self.transcript_id))
        if _exon_interval.chrom != self.interval.chrom:
            raise Exception('Transcript/Exon chromosome mismatch in {0}'.format(self.transcript_id))
        self.exon_list.append(Exon(_exon_interval, _exon_number))

    def add_utr5(self, _feature_interval):
        if not self.editable:
            raise Exception('5UTR addition on editable transcript {0}'.format(self.transcript_id))
        self.utr5_list.append(_feature_interval)

    def add_utr3(self, _feature_interval):
        if not self.editable:
            raise Exception('3UTR addition on editable transcript {0}'.format(self.transcript_id))
        self.utr3_list.append(_feature_interval)

    def add_cds(self, _feature_interval):
        if not self.editable:
            raise Exception('CDS addition on editable transcript {0}'.format(self.transcript_id))
        self.cds_list.append(_feature_interval)

    def finalize(self):
        if not self.editable:
            raise Exception('Finalization on editable transcript {0}'.format(self.transcript_id))
        # 0. fix transcripts without exons like pre-miRNA
        if not self.exon_list:
            self.add_exon(self.interval, 1)

        # 1. confirm exon in correct sequence
        have_exon_numbering = True
        exon_sortable = True
        for exon in self.exon_list:
            if not exon.exon_number:
                have_exon_numbering = False
                break
        if have_exon_numbering:
            self.exon_list = sorted(self.exon_list, key=lambda v: v.exon_number)
            for i, exon in enumerate(self.exon_list):
                if not exon.exon_number == i + 1:
                    logging.critical([x.exon_number for x in self.exon_list])
                    raise Exception('Unexpected exon numbering in {0}'.format(self.transcript_id))
        elif self.interval.strand == '+':
            self.exon_list = sorted(self.exon_list, key=lambda v: v.interval.start)
            for i in range(len(self.exon_list)):
                self.exon_list[i].exon_number = i + 1
        elif self.interval.strand == '-':
            self.exon_list = sorted(self.exon_list, key=lambda v: v.interval.end, reverse=True)
            for i in range(len(self.exon_list)):
                self.exon_list[i].exon_number = i + 1
        else:  # This must be some trans-spliced transcript
            exon_sortable = False
            logging.warning('Unable to infer splicing junction in {0}'.format(self.transcript_id))

        if exon_sortable:
            for i in range(len(self.exon_list) - 1):
                left_exon = self.exon_list[i]
                right_exon = self.exon_list[i + 1]
                junction_chrom = left_exon.interval.chrom

                junction_start_strand = left_exon.interval.strand
                if left_exon.interval.strand == '+':
                    junction_start = left_exon.interval.end - 1  # 0-based
                    junction_exon_leftmost = left_exon.interval.start
                elif left_exon.interval.strand == '-':
                    junction_start = left_exon.interval.start
                    junction_exon_leftmost = left_exon.interval.end - 1
                else:
                    raise Exception(
                        'No strand for exon {0} in transcript {1}'.format(left_exon.exon_number, self.transcript_id))

                junction_end_strand = right_exon.interval.strand
                if right_exon.interval.strand == '+':
                    junction_end = right_exon.interval.start
                    junction_exon_rightmost = right_exon.interval.end - 1
                elif right_exon.interval.strand == '-':
                    junction_end = right_exon.interval.end - 1  # 0-based
                    junction_exon_rightmost = right_exon.interval.start
                else:
                    raise Exception(
                        'No strand for exon {0} in transcript {1}'.format(right_exon.exon_number, self.transcript_id))
                self.junction_list.append(
                    Junction(junction_chrom, junction_start, junction_end, junction_start_strand, junction_end_strand,
                             junction_exon_leftmost, junction_exon_rightmost))
            if not len(self.junction_list) == len(self.exon_list) - 1:
                raise Exception(
                    'Number of junctions and exons does not match in transcript {0}'.format(self.transcript_id))
        if not self.transcript_type:
            if self.cds_list:
                self.transcript_type = 'protein_coding'
            elif not self.utr3_list and not self.utr5_list and not self.cds_list:
                self.transcript_type = 'unknown_non_coding'

        # 2. Infer 5'UTR and/or 3'UTR when only CDS annotation is available
        if exon_sortable and self.cds_list:
            self.utr5_list = []
            self.utr3_list = []
            at_5utr = True
            at_3utr = False
            for exon in self.exon_list:
                overlap_cds = False
                for cds_interval in self.cds_list:
                    if cds_interval.overlap_with(exon.interval):
                        overlap_cds = True
                        assert at_5utr != at_3utr
                        if exon.interval.strand == '+':
                            if at_5utr:
                                if exon.interval.start < cds_interval.start:
                                    self.utr5_list.append(Interval(exon.interval.chrom,
                                                                   exon.interval.start,
                                                                   cds_interval.start,
                                                                   exon.interval.strand))
                                at_5utr = False
                                at_3utr = True
                            if at_3utr:  # 1 exon can have all 5'UTR, CDS and 3'UTR
                                if exon.interval.end > cds_interval.end:
                                    self.utr3_list.append(Interval(exon.interval.chrom,
                                                                   cds_interval.end,
                                                                   exon.interval.end,
                                                                   exon.interval.strand))
                        elif exon.interval.strand == '-':
                            if at_5utr:
                                if exon.interval.end > cds_interval.end:
                                    self.utr5_list.append(Interval(exon.interval.chrom,
                                                                   cds_interval.end,
                                                                   exon.interval.end,
                                                                   exon.interval.strand))
                                at_5utr = False
                                at_3utr = True
                            if at_3utr:  # 1 exon can have all 5'UTR, CDS and 3'UTR
                                if exon.interval.start < cds_interval.start:
                                    self.utr3_list.append(Interval(exon.interval.chrom,
                                                                   exon.interval.start,
                                                                   cds_interval.start,
                                                                   exon.interval.strand))
                        break
                if not overlap_cds:
                    assert at_5utr != at_3utr
                    if at_5utr:
                        self.utr5_list.append(exon.interval)
                    elif at_3utr:
                        self.utr3_list.append(exon.interval)
        self.editable = False
        return self.transcript_type  # In case no gene_type is given in annotation, infer from transcript types

    def get_utrcds_annotation(self, _query_interval_list):
        if self.editable:
            raise Exception('UTR/CDS query on non-editable transcript')
        if not self.cds_list and not self.utr3_list and not self.utr5_list:
            return ['non_coding']
        output_accumulator = []
        for feature_interval in self.utr5_list:
            if '5UTR' in output_accumulator:
                break
            for query_interval in _query_interval_list:
                if feature_interval.overlap_with(query_interval):
                    output_accumulator.append('5UTR')
                    break

        for feature_interval in self.cds_list:
            if 'CDS' in output_accumulator:
                break
            for query_interval in _query_interval_list:
                if feature_interval.overlap_with(query_interval):
                    output_accumulator.append('CDS')
                    break

        for feature_interval in self.utr3_list:
            if '3UTR' in output_accumulator:
                break
            for query_interval in _query_interval_list:
                if feature_interval.overlap_with(query_interval):
                    output_accumulator.append('3UTR')
                    break

        return output_accumulator

    def dump_junction_groups(self):  # Generator Function
        # Junctions need different treatment from exons for whitelisting,
        # since exon length may <10 (limit of seq filter)
        # For this, we fallback to original implementation of array of genomic positions + sequence in string
        if self.editable:
            raise Exception('Dumping Junctions from non-editable transcript annotation.')

        # For transcripts that does not have valid junctions
        if not self.junction_list:
            return

        all_junction_groups = []
        current_junction_group = []
        for exon_index, exon_obj in enumerate(self.exon_list):
            if exon_index == 0:  # Length of first exon doesnt matter in any case
                continue
            junction_index = exon_index - 1
            current_junction_group.append(junction_index)
            if exon_obj.interval.len >= SEQ_FILTER_SEARCH_RANGE:
                all_junction_groups.append(current_junction_group)
                current_junction_group = []
        if current_junction_group:
            all_junction_groups.append(current_junction_group)

        for junction_group in all_junction_groups:
            yield (self.junction_list[i] for i in junction_group)


class Exon:
    def __init__(self, _exon_interval, _exon_number):
        self.interval = _exon_interval
        self.exon_number = _exon_number

    def get_whitelist_iv(self):
        if self.interval.strand == '+':
            new_start = min(self.interval.start + EXON_END_BLACKLIST_RANGE, self.interval.end)
            if new_start == self.interval.end:
                return None
            return HTSeq.GenomicInterval(self.interval.chrom,
                                         max(self.interval.start + EXON_END_BLACKLIST_RANGE, self.interval.end - 1),
                                         self.interval.end, self.interval.strand)
        elif self.interval.strand == '-':
            new_end = max(self.interval.end - EXON_END_BLACKLIST_RANGE, self.interval.start)
            if new_end == self.interval.start:
                return None
            return HTSeq.GenomicInterval(self.interval.chrom, self.interval.start,
                                         self.interval.end - EXON_END_BLACKLIST_RANGE, self.interval.strand)

    def get_iv(self):
        return HTSeq.GenomicInterval(self.interval.chrom, self.interval.start, self.interval.end, self.interval.strand)


class Junction:
    def __init__(self, _chrom, _start, _end, _start_strand, _end_strand, _exon_leftmost, _exon_rightmost):
        self.chrom = _chrom
        self.start = _start
        self.end = _end
        self.start_strand = _start_strand
        self.end_strand = _end_strand
        self.exon_leftmost = _exon_leftmost
        self.exon_rightmost = _exon_rightmost

    def get_whitelist_iv_pair(self):
        if self.start_strand == '+':
            start_iv = Interval(self.chrom,
                                max(self.start - SEQ_FILTER_SEARCH_RANGE + 1, self.exon_leftmost),
                                self.start + 1, self.start_strand)
        elif self.start_strand == '-':
            start_iv = Interval(self.chrom, self.start,
                                min(self.start + SEQ_FILTER_SEARCH_RANGE, self.exon_leftmost + 1),
                                self.start_strand)
        else:
            raise Exception('Junction start does not have strand')

        if self.end_strand == '+':
            end_iv = Interval(self.chrom, self.end,
                              min(self.end + EXON_END_BLACKLIST_RANGE, self.exon_rightmost + 1),
                              self.end_strand)
        elif self.end_strand == '-':
            end_iv = Interval(self.chrom,
                              max(self.end - EXON_END_BLACKLIST_RANGE + 1, self.exon_rightmost),
                              self.end + 1, self.end_strand)
        else:
            raise Exception('Junction end does not have strand')

        return start_iv, end_iv


class Interval:
    def __init__(self, _chrom, _start, _end, _strand):
        self.chrom = _chrom
        self.start = int(_start)
        self.end = int(_end)
        self.strand = _strand
        if (self.end - self.start) < 1:
            raise ValueError('Interval end is smaller than start.')
        self.len = self.end - self.start

    def __str__(self):
        return self.return_coordinates()

    def is_same(self, _interval):
        if (self.chrom == _interval.chrom and
                self.start == _interval.start and
                self.end == _interval.end and
                self.strand == _interval.strand):
            return True
        return False

    def overlap_with(self, _interval):
        if (self.chrom == _interval.chrom and
                self.strand == _interval.strand and
                self.start < _interval.end and
                self.end > _interval.start):
            return True
        return False

    def parse_generate_pos_strand_tuple(self):
        if self.strand == '+':
            for i in range(self.start, self.end):
                yield (i, self.strand)
        elif self.strand == '-':
            for i in range(self.end - 1, self.start - 1, -1):
                yield (i, self.strand)

    def return_position_list(self):
        if self.strand == '+':
            return list(range(self.start, self.end))
        elif self.strand == '-':
            return list(range(self.end - 1, self.start - 1, -1))

    def return_coordinates(self):
        return '{0}:{1}-{2}:{3}'.format(self.chrom, self.start, self.end, self.strand)

    def htseq_iv(self):
        return HTSeq.GenomicInterval(self.chrom, self.start, self.end, self.strand)


if __name__ == "__main__":
    exit(0)
