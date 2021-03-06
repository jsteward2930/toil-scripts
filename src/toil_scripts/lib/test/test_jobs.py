import os
from toil.job import Job


def test_sample_batcher(tmpdir):
    from toil_scripts.lib.jobs import map_job
    options = Job.Runner.getDefaultOptions(os.path.join(str(tmpdir), 'test_store'))
    samples = [x for x in xrange(200)]
    j = Job.wrapJobFn(map_job, _test_batch, samples, 'a', 'b', 'c')
    Job.Runner.startToil(j, options)


def _test_batch(job, sample, a, b, c):
    assert str(sample).isdigit()
    assert a == 'a'
    assert b == 'b'
    assert c == 'c'
