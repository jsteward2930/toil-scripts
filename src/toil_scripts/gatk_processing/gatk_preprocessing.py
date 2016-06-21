#!/usr/bin/env python2.7
"""
@author Jacob Pfeil
@data 01/13/2016

Toil pipeline for processing bam files for GATK halpotype calling

1 = download shared data
2 = reference preprocessing
3 = download sample
4 = remove secondary alignments
5 = index sample
6 = sort sample
7 = mark duplicates
8 = index
9 = realigner target
10 = indel realignment
11 = index
12 = base recalibration
13 = output bqsr fil6
"""
import argparse
import collections
import multiprocessing
import os
import sys
import logging
import textwrap
import yaml

from copy import deepcopy
from urlparse import urlparse

from toil.job import Job
from toil_scripts.lib import require
from toil_scripts.lib.files import copy_to
from toil_scripts.lib.urls import download_url_job, s3am_upload
from toil_scripts.lib.programs import docker_call

from bd2k.util.processes import which

_log = logging.getLogger(__name__)

# Convenience functions used in the pipeline


def upload_or_move(filename, work_dir=None, output_dir=None, s3_dir=None, s3_key_path=None):
    if output_dir:
        copy_to(filename, output_dir, work_dir)

    elif s3_dir:
        s3am_upload(fpath=os.path.join(work_dir, filename),
                    s3_dir=s3_dir,
                    s3_key_path=s3_key_path)

    else:
        raise ValueError('No output_directory or s3_dir defined. Cannot determine where to store %s' % filename)


def get_files_from_filestore(job, workDir, inputDict):
    """
    Given one or more strings representing file_names, return the paths to those files. Each item must be unpacked!

    work_dir: str       Current working directory
    ids: dict           Dictionary of fileStore IDs
    *args: str(s)       for every file in *args, place file in work_dir via FileStore
    """
    paths = collections.OrderedDict()
    for name, fileStoreID in inputDict.iteritems():
        if not os.path.exists(os.path.join(workDir, name)):
            file_path = job.fileStore.readGlobalFile(fileStoreID, os.path.join(workDir, name))
        else:
            file_path = name
        inputDict[name] = file_path
    return inputDict


def samtools_faidx(job, fileStoreID, mock=False):
    """
    Uses Samtools to create reference index file (.fasta.fai)

    ref_id: str     The fileStore ID of the reference
    """
    work_dir = job.fileStore.getLocalTempDir()

    inputs = {'genome.fa': fileStoreID}
    get_files_from_filestore(job, work_dir, inputs)

    # Call: Samtools
    command = ['faidx', 'genome.fa']
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/samtools:0.1.19--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=inputs.keys(),
                outputs={'genome.fa.fai': None},
                mock=mock)
    output = os.path.join(work_dir, 'genome.fa.fai')
    assert os.path.exists(output)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(output)


def picard_CreateSequenceDictionary(job, fileStoreID, memory=8, mock=False):
    """
    Uses Picardtools to create reference dictionary (.dict) for the sample

    ref_id: str     The fileStore ID of the reference
    input_args: dict        Dictionary of input arguments (from main())
    """
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'genome.fa': fileStoreID}
    get_files_from_filestore(job, work_dir, inputs)
    # Call: picardtools
    command = ['CreateSequenceDictionary', 'R=genome.fa', 'O=genome.dict']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=inputs.keys(),
                outputs={'genome.dict': None},
                mock=mock)
    # Write to fileStore
    return job.fileStore.writeGlobalFile(os.path.join(work_dir, 'genome.dict'))


def download_gatk_files(job, config):
    """
    Downloads files shared by all samples in the pipeline

    input_args: dict        Dictionary of input arguments (from main())
    """
    reference_names = ['genome.fa', 'phase.vcf', 'mills.vcf', 'dbsnp.vcf']
    for name in reference_names:
        key, _ = name.split('.')
        config[name] = job.addChildJobFn(download_url_job, url=config[key], name=name,
                                             s3_key_path=config['ssec']).rv()
    return config


def save_file(job, fileStoreID, fname, config):
    work_dir = job.fileStore.getLocalTempDir()
    work_path = os.path.join(work_dir, fname)
    job.fileStore.readGlobalFile(fileStoreID, work_path)
    upload_or_move(job, work_dir, config, fname)


def reference_preprocessing(job, config, mock=False):
    """
    Create index and dict file for reference

    input_args: dict        Dictionary of input argumnets
    shared_ids: dict        Dictionary of fileStore IDs
    """
    fileStoreID = config['genome.fa']
    if 'genome.fa.fai' not in config:
        config['genome.fa.fai'] = job.addChildJobFn(samtools_faidx, fileStoreID, mock=mock).rv()
    if 'genome.dict' not in config:
        config['genome.dict'] = job.addChildJobFn(picard_CreateSequenceDictionary, fileStoreID, mock=mock).rv()
    return config


def remove_supplementary_alignments(job, fileStoreID, mock=False):
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'sample.bam': fileStoreID}
    get_files_from_filestore(job, work_dir, inputs)
    outpath = os.path.join(work_dir, 'sample.nosuppl.bam')

    command = ['view',
               '-b', '-o', '/data/sample.nosuppl.bam',
               '-F', '0x800',
               '/data/sample.bam']
               
    docker_call(work_dir=work_dir, parameters=command,
                tool='quay.io/ucsc_cgl/samtools:1.3--256539928ea162949d8a65ca5c79a72ef557ce7c',
                inputs=inputs.keys(),
                outputs={'sample.nosuppl.bam': None},
                mock=mock)
    return job.fileStore.writeGlobalFile(outpath)


def _move_bai(outpath):

    # bam index gets written without bam extension
    name = outpath.replace('bam', 'bai')
    newname = outpath + '.bai'
    os.rename(name, newname)
    return newname


def picard_SortSam(job, fileStoreID, config):
    """
    Uses picardtools SortSam to sort a sample bam file

    Assumes machine has at least 4G of RAM
    """
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'sample.nosuppl.bam': fileStoreID}
    memory = config['memory'] if 'memory' in config else '4'
    get_files_from_filestore(job, work_dir, inputs)
    outpath = os.path.join(work_dir, 'sample.sorted.bam')

    #Call: picardtools
    command = ['SortSam',
               'INPUT=sample.nosuppl.bam',
               'OUTPUT=sample.sorted.bam',
               'SORT_ORDER=coordinate',
               'CREATE_INDEX=true']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=inputs.keys(),
                outputs={'sample.sorted.bam': None, 'sample.sorted.bai': None},
                mock=config['mock'])
    bam_id = job.fileStore.writeGlobalFile(outpath)
    bai_id = job.fileStore.writeGlobalFile(_move_bai(outpath))
    return bam_id, bai_id


def picard_MarkDuplicates(job, bamFileStoreID, baiFileStoreID, config):
    """
    Uses picardtools MarkDuplicates
    """
    work_dir = job.fileStore.getLocalTempDir()
    inputs = {'sample.sorted.bam': bamFileStoreID,
              'sample.sorted.bai': baiFileStoreID}

    memory = config['memory'] if config['memory'] is not None else '4'

    # Retrieve file path
    get_files_from_filestore(job, work_dir, inputs)
    outpath = os.path.join(work_dir, 'sample.mkdups.bam')
    # Call: picardtools
    command = ['MarkDuplicates',
               'INPUT=sample.sorted.bam',
               'OUTPUT=sample.mkdups.bam',
               'METRICS_FILE=metrics.txt',
               'ASSUME_SORTED=true',
               'CREATE_INDEX=true']
    docker_call(work_dir=work_dir, parameters=command,
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                tool='quay.io/ucsc_cgl/picardtools:1.95--dd5ac549b95eb3e5d166a5e310417ef13651994e',
                inputs=inputs.keys(), mock=config['mock'],
                outputs={'sample.mkdups.bam': None, 'sample.mkdups.bai': None})

    # picard writes the index for file.bam at file.bai, not file.bam.bai
    _move_bai(outpath)
    bam_id = job.fileStore.writeGlobalFile(outpath)
    bai_id = job.fileStore.writeGlobalFile(outpath + ".bai")
    return bam_id, bai_id
    # job.addChildJobFn(realigner_target_creator, shared_ids, input_args, cores = multiprocessing.cpu_count())


def gatk_RealignerTargetCreator(job, bamFileStoreID, baiFileStoreID, config):
    """
    Creates <type>.intervals file needed for indel realignment

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()
    memory = config['memory'] if config['memory'] is not None else '4'
    inputs = {'sample.sorted.bam': bamFileStoreID,
              'sample.sorted.bai': baiFileStoreID}
    references = {key: config[key] for key in ['genome.fa','genome.fa.fai', 'genome.dict', 'phase.vcf', 'mills.vcf']}
    inputs.update(references)
    get_files_from_filestore(job, work_dir, inputs)
    # Output file path
    output = os.path.join(work_dir, 'sample.intervals')
    # Call: GATK -- RealignerTargetCreator
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'RealignerTargetCreator',
                  '-nt', str(cores),
                  '-R', 'ref.fa',
                  '-I', 'sample.mkdups.bam',
                  '-known', 'phase.vcf',
                  '-known', 'mills.vcf',
                  '--downsampling_type', 'NONE',
                  '-o', 'sample.intervals']

    docker_call(work_dir=work_dir, parameters=parameters,
                tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                inputs=inputs.keys(),
                outputs={'sample.intervals': None},
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                mock=config['mock'])
    return job.fileStore.writeGlobalFile(output)
    # job.addChildJobFn(indel_realignment, shared_ids, input_args)


def gatk_IndelRealigner(job, bamFileStoreID, baiFileStoreID, intervalsFileStoreID, config):
    """
    Creates realigned bams using <sample>.intervals file from previous step

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    # Unpack convenience variables for job
    work_dir = job.fileStore.getLocalTempDir()
    memory = '10'
    # Retrieve input file paths
    inputs = {'sample.mkdups.bam': bamFileStoreID, 'sample.mkdups.bam.bai': baiFileStoreID,
              'sample.intervals': intervalsFileStoreID}
    inputs.update({key: config[key] for key in
                   ['genome.fa','genome.fa.fai', 'genome.dict', 'phase.vcf', 'mills.vcf']})
    get_files_from_filestore(job, work_dir, inputs)
    # Output file path
    outpath = os.path.join(work_dir, 'sample.indel.bam')
    # Call: GATK -- IndelRealigner
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'IndelRealigner',
                  '-R', 'ref.fa',
                  '-I', 'sample.mkdups.bam',
                  '-known', 'phase.vcf',
                  '-known', 'mills.vcf',
                  '-targetIntervals', 'sample.intervals',
                  '--downsampling_type', 'NONE',
                  '-maxReads', str(720000),
                  '-maxInMemory', str(5400000),
                  '-o', 'sample.indel.bam']

    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                work_dir=work_dir, parameters=parameters,
                inputs=inputs,
                outputs={'sample.indel.bam': None, 'sample.indel.bai': None},
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                mock=config['mock'])

    # Write to fileStore
    bam_id = job.fileStore.writeGlobalFile(outpath)
    bai_id = job.fileStore.writeGlobalFile(_move_bai(outpath))
    return bam_id, bai_id
    # job.addChildJobFn(base_recalibration, shared_ids, input_args, cores = multiprocessing.cpu_count())


def gatk_BaseRecalibrator(job, bamFileStoreID, baiFileStoreID, config):
    """
    Creates recal table to perform Base Quality Score Recalibration

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """

    job.fileStore.logToMaster('Running GATK BaseRecalibrator: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()
    # TODO Do someting better with allocating memory
    memory = '10'
    # Retrieve input file paths
    inputs = {'sample.indel.bam': bamFileStoreID, 'sample.indel.bam.bai': baiFileStoreID}
    inputs.update({key: config[key] for key in
                   ['genome.fa','genome.fa.fai', 'genome.dict', 'dbsnp.vcf']})
    get_files_from_filestore(job, work_dir, inputs)
    # Output file path
    output = os.path.join(work_dir, 'sample.recal.table')

    # Call: GATK -- IndelRealigner
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'BaseRecalibrator',
                  '-nct', str(cores),
                  '-R', 'genome.fa',
                  '-I', 'sample.indel.bam',
                  '-knownSites', 'dbsnp.vcf',
                  '-o', 'sample.recal.table']
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                work_dir=work_dir, parameters=parameters,
                inputs=inputs.keys(),
                outputs={'sample.recal.table': None},
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                mock=config['mock'])
    # Write to fileStore
    return job.fileStore.writeGlobalFile(output)
    # job.addChildJobFn(print_reads, shared_ids, input_args, cores = cores)


def gatk_PrintReads(job, bamFileStoreID, baiFileStoreID, recalFileStoreID, config):
    """
    Create bam that has undergone Base Quality Score Recalibration (BQSR)

    job_vars: tuple     Contains the input_args and ids dictionaries
    sample: str         Either "normal" or "tumor" to track which one is which
    """
    job.fileStore.logToMaster('Running GATK PrintReads: {}'.format(config['uuid']))
    work_dir = job.fileStore.getLocalTempDir()
    cores = multiprocessing.cpu_count()
    # TODO Do someting better with allocating memory
    memory = '10'
    # Retrieve input file paths
    inputs = {'sample.indel.bam': bamFileStoreID, 'sample.indel.bam.bai': baiFileStoreID,
              'recal.table': recalFileStoreID}
    inputs.update({key: config[key] for key in
                   ['genome.fa','genome.fa.fai', 'genome.dict']})
    get_files_from_filestore(job, work_dir, inputs)
    # Output file
    # Unpack convenience variables for job
    uuid = config['uuid']
    suffix = config['suffix'] if config['suffix'] is not None else ''
    outfile = '{}{}.bam'.format(uuid, suffix)
    gatk_outfile_idx = '{}{}.bai'.format(uuid, suffix)
    outpath = os.path.join(work_dir, outfile)
    # Call: GATK -- PrintReads
    parameters = ['-U', 'ALLOW_SEQ_DICT_INCOMPATIBILITY', # RISKY! (?) See #189
                  '-T', 'PrintReads',
                  '-nct', str(cores),
                  '-R', 'ref.fa',
                  '--emit_original_quals',
                  '-I', 'sample.indel.bam',
                  '-BQSR', 'sample.recal.table',
                  '-o', outfile]
    docker_call(tool='quay.io/ucsc_cgl/gatk:3.5--dba6dae49156168a909c43330350c6161dc7ecc2',
                work_dir=work_dir, parameters=parameters,
                inputs=inputs,
                outputs={outfile: None, gatk_outfile_idx: None},
                env={'JAVA_OPTS':'-Xmx%sg' % memory},
                mock=config['mock'])
    
    upload_or_move(outpath, output_dir=config['output_dir'], s3_dir=config['s3_dir'], s3_key_path=config['ssec'])
    upload_or_move(_move_bai(outpath), output_dir=config['output_dir'],
                   s3_dir=config['s3_dir'], s3_key_path=config['ssec'])

def generate_file(file_path, generate_func):
    require(not os.path.exists(file_path), file_path + ' already exists!')
    with open(file_path, 'w') as f:
        f.write(generate_func())
    print('\t{} has been generated in the current working directory.'.format(os.path.basename(file_path)))

def generate_config():
    return textwrap.dedent("""
    # CGL Germline Variant Pipeline configuration file
    # This configuration file is formatted in YAML. Simply write the value (at least one space) after the colon.
    # Edit the values in this configuration file and then rerun the pipeline: "toil-variant run"
    # URLs can take the form: http://, file://, s3://, gnos://.
    # Comments (beginning with #) do not need to be removed. Optional parameters may be left blank
    ####################################################################################################################
    # Required: URL to reference genome
    genome: s3://cgl-pipeline-inputs/variant_hg19/hg19.fa

    # Required: URL to phase indels VCF
    phase: s3://cgl-pipeline-inputs/variant_hg19/1000G_phase1.indels.hg19.sites.vcf

    # Required: URL to Mills indel VCF
    mills: s3://cgl-pipeline-inputs/variant_hg19/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf

    # Required: URL to dbsnp VCF
    dbsnp: s3://cgl-pipeline-inputs/variant_hg19/dbsnp_138.hg19.vcf

    # Required: URL to cosmic VCF
    cosmic: s3://cgl-pipeline-inputs/variant_hg19/cosmic.hg19.vcf

    # Optional: Provide a full path to where results will appear
    output-dir:

    # Optional: Provide an s3 path (s3://bucket/dir) where results will appear
    s3-dir:

    # Optional: Provide a full path to a 32-byte key used for SSE-C Encryption in Amazon
    ssec:

    # Optional: Provide a full path to a CGHub Key used to access GNOS hosted data
    gtkey:

    # Optional: Provide a suffix for naming output files
    suffix:

    # Optional: Provide a RAM specification
    memory:
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        #   Edit this manifest to include information pertaining to each sample pair to be run.
        #   There are 2tab-separated columns: UUID, Sample BAM URL
        #
        #   UUID            This should be a unique identifier for the sample to be processed
        #   Sample URL      A URL (http://, file://, s3://, gnos://) pointing to the normal bam
        #
        #   Examples of several combinations are provided below. Lines beginning with # are ignored.
        #
        #   UUID_1  file:///path/to/sample.bam
        #   UUID_2  http://sample-depot.com/sample.bam
        #   UUID_3  s3://my-bucket-name/directory/sample.bam
        #
        #   Place your samples below, one per line.
        """[1:])

def parse_manifest(path_to_manifest):
    """
    Parses samples, specified in either a manifest or listed with --samples

    :param str path_to_manifest: Path to configuration file
    :return: Samples and their attributes as defined in the manifest
    :rtype: list[list]
    """
    samples = []
    with open(path_to_manifest, 'r') as f:
        for line in f.readlines():
            if not line.isspace() and not line.startswith('#'):
                sample = line.strip().split('\t')
                require(len(sample) == 2, 'Bad manifest format! '
                                          'Expected 3 tab separated columns, got: {}'.format(sample))
                uuid, url = sample
                require(urlparse(url).scheme and urlparse(url), 'Invalid URL passed for {}'.format(url))
                samples.append(sample)
    return samples

def gatk_preprocessing_pipeline(job, uuid, url, config):
    config = deepcopy(config)
    config['uuid'] = uuid
    download = job.wrapJobFn(download_url_job, url, name='sample.bam', s3_key_path=config['ssec'])
    rm_secondary = job.wrapJobFn(remove_supplementary_alignments, download.rv(), mock=config['mock'])
    picard_sort = job.wrapJobFn(picard_SortSam, rm_secondary.rv(), config)
    mdups = job.wrapJobFn(picard_MarkDuplicates, picard_sort.rv(0), picard_sort.rv(1), config)
    target = job.wrapJobFn(gatk_RealignerTargetCreator, mdups.rv(0), mdups.rv(1), config)
    indel = job.wrapJobFn(gatk_IndelRealigner, mdups.rv(0), mdups.rv(1), target.rv(), config)
    base = job.wrapJobFn(gatk_BaseRecalibrator, indel.rv(0), indel.rv(1), config)
    print_reads = job.wrapJobFn(gatk_PrintReads, indel.rv(0), indel.rv(1), base.rv(), config)


    job.addChild(download)
    download.addChild(rm_secondary)
    rm_secondary.addChild(picard_sort)
    picard_sort.addChild(mdups)
    mdups.addChild(target)
    target.addChild(indel)
    indel.addChild(base)
    base.addChild(print_reads)


def main():
    """
    GATK Pre-processing Script
    """
    parser = argparse.ArgumentParser(description=main.__doc__, formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(dest='command')
    subparsers.add_parser('generate-config', help='Generates an editable config in the current working directory.')
    subparsers.add_parser('generate-manifest', help='Generates an editable manifest in the current working directory.')
    subparsers.add_parser('generate', help='Generates a config and manifest in the current working directory.')


    # Run subparser
    parser_run = subparsers.add_parser('run', help='Runs the CGL germline pipeline')
    parser_run.add_argument('--config', default='config-toil-germline.yaml', type=str,
                            help='Path to the (filled in) config file, generated with "generate-config". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--manifest', default='manifest-toil-germline.tsv', type=str,
                            help='Path to the (filled in) manifest file, generated with "generate-manifest". '
                                 '\nDefault value: "%(default)s"')
    parser_run.add_argument('--bam', default=None, type=str,
                            help='URL for the sample BAM. URLs can take the form: http://, file://, s3://, '
                                 'and gnos://. The UUID for the sample must be given with the "--uuid" flag.')
    parser_run.add_argument('--uuid', default=None, type=str, help='Provide the UUID of a sample when using the'
                                                                   '"--bam" option')

    # If no arguments provided, print full help menu
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    Job.Runner.addToilOptions(parser)
    args = parser.parse_args()

    cwd = os.getcwd()
    if args.command == 'generate-config' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'config-toil-germline.yaml'), generate_config)
    if args.command == 'generate-manifest' or args.command == 'generate':
        generate_file(os.path.join(cwd, 'manifest-toil-germline.tsv'), generate_manifest)
    if 'generate' in args.command:
        sys.exit()
    if args.command == 'run':
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).iteritems()}

        require(set(config) > {'genome', 'mills', 'dbsnp', 'phase'},
                'Missing inputs for preprocessing, check config file')

        # Program checks
        for program in ['curl', 'docker']:
            require(next(which(program)), program + ' must be installed on every node.'.format(program))

        if args.bam or args.uuid:
            require(args.bam and args.uuid, '"--bam" and "--uuid" must all be supplied')
            samples = [[args.uuid, args.normal, args.tumor]]
        else:
            samples = parse_manifest(args.manifest)

        root = Job.wrapJobFn(download_gatk_files, config)
        ref = Job.wrapJobFn(reference_preprocessing, root.rv(), mock=config['mock'])

        root.addFollowOn(ref)

        for uuid, url in samples:
            ref.addFollowOnJobFn(gatk_preprocessing_pipeline, uuid, url, ref.rv())

        Job.Runner.startToil(root, args)

if __name__ == '__main__':
    main()
