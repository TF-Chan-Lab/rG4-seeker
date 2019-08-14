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

import configparser
import sys
import multiprocessing

import logging
import time
# rg4_seeker imports
from rg4seeker import apriori_rts_analysis, bam_preprocess, pile_upper, whitelister, fsr_precompute, fsr_cutoff_finder


def main_driver(ini_file):
    logging.basicConfig(level='INFO')
    config = configparser.ConfigParser()
    config.read(ini_file)

    # global section
    working_dir = config.get('global', 'WORKING_DIR')
    sample_name = config.get('global', 'SAMPLE_NAME')
    try:
        threads = config.getint('global', 'THREADS')
    except configparser.NoOptionError:
        threads = 1
        logging.warning(
            'Caution: number of threads to use not specified, assuming running rG4-seeker single-threaded (slow)')
    threads = max(threads, 1)
    annotation_number = config.getint('global', 'NO_OF_ANNOTATIONS')
    replicate_number = config.getint('global', 'NO_OF_REPLICATES')
    have_kpds_condition = config.getboolean('global', 'HAVE_KPDS_CONDITION')
    aligner = config.get('global', 'ALIGNER')
    if aligner not in ['STAR']:
        logging.critical('Sorry, currently only BAM files from STAR aligner output is supported')
        raise Exception('Sorry, currently only BAM files from STAR aligner output is supported')
    reads_type = config.get('global', 'READS_TYPE')
    if reads_type not in ['SE', 'PE']:
        logging.critical('Please specify reads type as either "SE" or "PE".')
        raise Exception('Please specify reads type as either "SE" or "PE".')
    if 'TARGET_FDR' in config['global']:
        target_fdr = float(config.get('global', 'TARGET_FDR'))
    else:
        target_fdr = 0.015

    # genome setting
    genome_fasta = config.get('genome', 'GENOME_FASTA')

    # annotation
    annotation_names = []
    annotation_gffs = []
    for i in range(1, annotation_number + 1):
        annotation_names.append(config.get('annotation_{0}'.format(i), 'ANNOTATION_NAME'))
        annotation_gffs.append(config.get('annotation_{0}'.format(i), 'ANNOTATION_GFF'))

    # debug section
    no_of_bootstraps = 0
    parallel_preprocess = 1
    skip_steps = set()
    start_bootstrap_id = 0
    dump_fsrcutofffinder = False
    overwrite_reusables = False
    starting_replicate = 1
    if 'debug' in config:
        if 'PARALLEL_PREPROCESS' in config['debug']:
            parallel_preprocess = int(config.get('debug', 'PARALLEL_PREPROCESS'))
        if 'OMIT' in config['debug']:
            skip_steps = [x.strip() for x in config.get('debug', 'OMIT').split(',')]
        if 'DUMP_FSRCUTOFFFINDER' in config['debug']:
            dump_fsrcutofffinder = config.getboolean('debug', 'DUMP_FSRCUTOFFFINDER')
        if 'OVERWRITE_REUSABLES' in config['debug']:
            overwrite_reusables = config.getboolean('debug', 'OVERWRITE_REUSABLES')

    # replicate bam file section
    replicate_bam_files = {}
    for i in range(starting_replicate, replicate_number + 1):
        rep = 'replicate_{0}'.format(i)
        replicate_bam_files[i] = {'Li': config.get(rep, 'LI_BAM_FILE'),
                                  'K': config.get(rep, 'K_BAM_FILE'), }
        if have_kpds_condition:
            replicate_bam_files[i]['KPDS'] = config.get(rep, 'KPDS_BAM_FILE')
    if skip_steps:
        logging.info('[Debug] omit steps: {0}'.format(','.join(skip_steps)))

    # Execution
    logging.info(
        '[{0}] rG4-seeker: Pre-computing putative quadruplex sequences from reference genome and gene annotations'.format(
            time.ctime()))
    whitelist_beds = []
    rg4_motif_beds = []
    for i in range(annotation_number):
        annotation_name = annotation_names[i]
        annotation_gff = annotation_gffs[i]
        prefix = working_dir + '/' + sample_name + '.' + annotation_name
        if 'whitelister' not in skip_steps:
            whitelist_bed, rg4_motif_bed = whitelister.execute(annotation_gff, genome_fasta, prefix, threads,
                                                               overwrite_reusables)
            whitelist_beds.append(whitelist_bed)
            rg4_motif_beds.append(rg4_motif_bed)
        else:
            whitelist_beds.append(prefix + '.rG4_whitelist.bed.gz')
            rg4_motif_beds.append(prefix + '.rG4_motif.bed.gz')

    if have_kpds_condition:
        preprocess_conditions = ['Li', 'K', 'KPDS']
        analysis_conditions = ['K', 'KPDS']
    else:
        preprocess_conditions = ['Li', 'K']
        analysis_conditions = ['K']

    logging.info('[{0}] rG4-seeker: Pre-processing input aligned reads'.format(time.ctime()))
    if 'bam_preprocess' not in skip_steps:
        preprocess_thread_list = []
        for replicate_id in range(starting_replicate, replicate_number + 1):
            for condition in preprocess_conditions:
                prefix = working_dir + '/{0}-{1}-rep{2}'.format(sample_name, condition, replicate_id)
                preprocess_thread_list.append(bam_preprocess.BamPreprocess(replicate_bam_files[replicate_id][condition],
                                                                           prefix,
                                                                           aligner,
                                                                           reads_type,
                                                                           no_of_bootstraps))
        with multiprocessing.Pool(parallel_preprocess) as pool:
            multiple_results = [pool.apply_async(i.run, args=()) for i in preprocess_thread_list]
            [res.get() for res in multiple_results]

    logging.info('[{0}] rG4-seeker: Piling read coverage and 5\' start site information'.format(time.ctime()))
    if 'pile_upper' not in skip_steps:
        for bootstrap_id in range(start_bootstrap_id, no_of_bootstraps + 1):
            pu_obj_list = []
            for replicate_id in range(starting_replicate, replicate_number + 1):
                for condition in preprocess_conditions:
                    logging.info(
                        '[{0}] Piling read coverage and 5\' start site information for sample {1}-rep{2} in {3} CPU threads'.format(
                            time.ctime(), condition, replicate_id, threads))
                    prefix = working_dir + '/{0}-{1}-rep{2}'.format(sample_name, condition, replicate_id)
                    chroms_list = []
                    with open(prefix + '.chromosomes.txt') as f:
                        for line in f:
                            chroms_list.append(line.strip())
                    curr_pu_obj = pile_upper.PileUpper(prefix, chroms_list, bootstrap_id)
                    curr_pu_obj.run(threads)
                    pu_obj_list.append(curr_pu_obj)
            with multiprocessing.Pool(threads) as pool:
                multiple_results = [pool.apply_async(i.concat, args=()) for i in pu_obj_list]
                [res.get() for res in multiple_results]

    logging.info('[{0}] rG4-seeker: Computing RSR peaks'.format(time.ctime()))
    if 'fsr_precompute' not in skip_steps:
        for bootstrap_id in range(start_bootstrap_id, no_of_bootstraps + 1):
            for replicate_id in range(starting_replicate, replicate_number + 1):
                for condition in analysis_conditions:
                    logging.info('[{0}] Computing RSR peaks for sample {1}-rep{2}'.format(time.ctime(), condition,
                                                                                          replicate_id))
                    treatment_prefix = working_dir + '/{0}-{1}-rep{2}'.format(sample_name, condition, replicate_id)
                    control_prefix = working_dir + '/{0}-{1}-rep{2}'.format(sample_name, 'Li', replicate_id)
                    treatment_chroms_list = []
                    with open(treatment_prefix + '.chromosomes.txt') as f:
                        for line in f:
                            treatment_chroms_list.append(line.strip())
                    control_chroms_list = []
                    with open(control_prefix + '.chromosomes.txt') as f:
                        for line in f:
                            control_chroms_list.append(line.strip())
                    control_chroms_list = set(control_chroms_list)
                    chroms_list = []
                    for chrom in treatment_chroms_list:
                        if chrom in control_chroms_list:
                            chroms_list.append(chrom)
                    fsrp_obj = fsr_precompute.FSRPrecompute(treatment_prefix, control_prefix, chroms_list, bootstrap_id)
                    fsrp_obj.run(threads)

    logging.info('[{0}] rG4-seeker: Determining per-replicate detection cutoffs'.format(time.ctime()))
    if 'fsr_cutoff_finder' not in skip_steps:
        prefix_list = []
        for bootstrap_id in range(start_bootstrap_id, no_of_bootstraps + 1):
            for replicate_id in range(starting_replicate, replicate_number + 1):
                for condition in analysis_conditions:
                    prefix = working_dir + '/{0}-{1}-rep{2}'.format(sample_name, condition, replicate_id)
                    prefix_list.append(prefix)
            ffc_obj = fsr_cutoff_finder.FSRCutoffFinder(prefix_list, genome_fasta, whitelist_beds, rg4_motif_beds,
                                                        bootstrap_id, dump_fsrcutofffinder, overwrite_reusables,
                                                        target_fdr)
            ffc_obj.run()

    logging.info('[{0}] rG4-seeker: Analyzing RTS sites'.format(time.ctime()))
    if 'apriori_rts_analysis' not in skip_steps:
        for bootstrap_id in range(start_bootstrap_id, no_of_bootstraps + 1):
            prefix = working_dir + '/{0}'.format(sample_name)
            apriori_rts_analysis.main(prefix, analysis_conditions, replicate_number, bootstrap_id, genome_fasta,
                                      whitelist_beds, rg4_motif_beds, annotation_names, annotation_gffs)

    logging.info('[{0}] rG4-seeker: All done!'.format(time.ctime()))


def main():
    if len(sys.argv) == 2:
        main_driver(sys.argv[1])
    else:
        print('rG4-seeker : A pipeline for processing and analyzing rG4-seq data')
        print('Usage:')
        print('( From manual installation ) rG4-seeker CONFIGURATION_FILE')
        print('( From docker container )    docker -v WORKING_DIR:WORKING_DIR run rg4_seeker CONFIGURATION_FILE')
        print('                             *Please put configuration file in the working directory*')
    return 0


if __name__ == "__main__":
    main()
