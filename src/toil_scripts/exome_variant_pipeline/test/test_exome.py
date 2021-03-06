import os
import subprocess
import shutil
import textwrap

from boto.s3.connection import S3Connection, Bucket, Key


def test_exome(tmpdir):
    work_dir = str(tmpdir)
    create_config_and_manifest(work_dir)
    subdir = '/mnt/ephemeral/toil-scripts/exome/'
    os.makedirs(os.path.join(subdir, 'workDir'))
    # normal_bam = 's3://cgl-pipeline-inputs/exome/ci/chr6.normal.bam'
    # tumor_bam = 's3://cgl-pipeline-inputs/exome/ci/chr6.tumor.bam'
    # Call Pipeline
    try:
        base_command = ['toil-exome', 'run',
                        '--config', os.path.join(work_dir, 'toil-exome.config'),
                        os.path.join(subdir, 'jstore'),
                        '--workDir', os.path.join(subdir, 'workDir')]
        # Run with --samples
        # subprocess.check_call(base_command + ['--samples', sample])
        # Run with manifest
        subprocess.check_call(base_command + ['--manifest', os.path.join(work_dir, 'toil-exome-manifest.tsv')])
    finally:
        shutil.rmtree(subdir)
        conn = S3Connection()
        b = Bucket(conn, 'cgl-driver-projects')
        k = Key(b)
        k.key = 'test/ci/exome-ci-test.tar.gz'
        k.delete()


def create_config_and_manifest(path):
    """Creates config file for test at path"""
    config_path = os.path.join(path, 'toil-exome.config')
    manifest_path = os.path.join(path, 'toil-exome-manifest.tsv')
    with open(config_path, 'w') as f:
        f.write(generate_config())
    with open(manifest_path, 'w') as f:
        f.write(generate_manifest())


def generate_config():
    return textwrap.dedent("""
    reference: s3://cgl-pipeline-inputs/exome/ci/chr6.hg19.fa
    phase: s3://cgl-pipeline-inputs/exome/ci/chr6.phase.hg19.vcf
    mills: s3://cgl-pipeline-inputs/exome/ci/chr6.Mills.indels.hg19.sites.vcf
    dbsnp: s3://cgl-pipeline-inputs/exome/ci/chr6.dbsnp.hg19.vcf
    cosmic: s3://cgl-pipeline-inputs/exome/ci/chr6.cosmic.hg19.vcf
    run-mutect: true
    run-pindel: true
    run-muse: true
    preprocessing: true
    output-dir:
    s3-dir: s3://cgl-driver-projects/test/ci
    ssec:
    gtkey:
    ci-test: true
    """[1:])


def generate_manifest():
    return textwrap.dedent("""
        exome-ci-test\ts3://cgl-pipeline-inputs/exome/ci/chr6.normal.bam\ts3://cgl-pipeline-inputs/exome/ci/chr6.tumor.bam
        """[1:])
