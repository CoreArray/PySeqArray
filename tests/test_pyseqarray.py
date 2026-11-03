"""Tests for PySeqArray.

Run under pytest (``pytest tests/``), with the example file
1KG_phase1_release_v3_chr22.gds (1092 samples, 19773 variants). Covers:
opening the file, GetData() with sample and variant filters, the filter
stack, Apply(), RunParallel() with forked processes and with a Pool which
opens the file again, and AlleleFreq().
"""

import multiprocessing as mp

import numpy as np
import pytest

import PySeqArray as ps


N_SAMPLE = 1092
N_VARIANT = 19773

has_fork = 'fork' in mp.get_all_start_methods()


@pytest.fixture
def f():
	f = ps.SeqArrayFile()
	f.open(ps.seqExample('1KG_phase1_release_v3_chr22.gds'))
	yield f
	f.close()


# defined at module level, so picklable for a Pool which opens the file again
def get_position(file, param):
	return file.GetData('position')

def get_sample_num(file, param):
	return [ len(file.GetData('sample.id')) ]


# ---------------------------------------------------------------------------

def test_open_example(f):
	sid = f.GetData('sample.id')
	assert sid.shape == (N_SAMPLE,)
	assert list(sid[:3]) == ['HG00096', 'HG00097', 'HG00099']
	pos = f.GetData('position')
	assert pos.dtype == np.int32 and pos.shape == (N_VARIANT,)
	assert list(pos[:2]) == [16051497, 16059752]
	assert list(f.GetData('$chrom_pos')[:2]) == ['22_16051497', '22_16059752']


def test_get_data(f):
	f.FilterSet2(sample=np.arange(10), variant=np.arange(5), verbose=False)
	g = f.GetData('genotype')
	assert g.dtype == np.uint8 and g.shape == (5, 10, 2)
	# dosage of the reference allele
	d = f.GetData('$dosage')
	assert list(d[0]) == [2, 1, 1, 1, 1, 1, 2, 2, 1, 0]
	assert np.array_equal(d, (g == 0).sum(axis=2))
	allele = [ s.split(',') for s in f.GetData('allele') ]
	assert list(f.GetData('$ref')) == [ s[0] for s in allele ]
	assert list(f.GetData('$alt')) == [ s[1] for s in allele ]
	assert list(f.GetData('$num_allele')) == [ len(s) for s in allele ]


def test_filter_set(f):
	sid = f.GetData('sample.id')
	vid = f.GetData('variant.id')
	f.FilterSet(sample_id=sid[[2, 5]], variant_id=vid[10:20], verbose=False)
	assert list(f.GetData('sample.id')) == list(sid[[2, 5]])
	assert np.array_equal(f.GetData('variant.id'), vid[10:20])
	assert f.FilterGet(True).sum() == 2
	assert f.FilterGet(False).sum() == 10
	f.FilterReset(verbose=False)
	assert len(f.GetData('sample.id')) == N_SAMPLE
	assert len(f.GetData('variant.id')) == N_VARIANT


def test_filter_set2(f):
	sid = f.GetData('sample.id')
	pos = f.GetData('position')
	# a range object and a list
	f.FilterSet2(sample=range(5), variant=[3, 1, 7], verbose=False)
	assert list(f.GetData('sample.id')) == list(sid[:5])
	assert list(f.GetData('position')) == list(pos[[1, 3, 7]])
	f.FilterReset(verbose=False)
	# a bool vector, not contiguous in memory
	b = np.zeros(2*N_SAMPLE, dtype=bool)
	b[::4] = True
	f.FilterSet2(sample=b[::2], verbose=False)
	assert list(f.GetData('sample.id')) == list(sid[::2])
	# the selected samples of the selected ones
	f.FilterSet2(sample=np.arange(3), intersect=True, verbose=False)
	assert list(f.GetData('sample.id')) == list(sid[0:6:2])


def test_filter_push_pop(f):
	f.FilterSet2(sample=np.arange(10), verbose=False)
	f.FilterPush()
	assert len(f.GetData('sample.id')) == N_SAMPLE
	f.FilterPop()
	assert len(f.GetData('sample.id')) == 10
	f.FilterPush(False)
	assert len(f.GetData('sample.id')) == 10
	f.FilterPop()


def test_apply(f):
	f.FilterSet2(variant=np.arange(100), verbose=False)
	pos = f.GetData('position')
	v = f.Apply('position', lambda x: x, asis='unlist', bsize=30)
	assert np.array_equal(v, pos)
	# more than one variable, with a parameter
	v = f.Apply(['position', 'genotype'], lambda p, x, g: (len(x) + p, g.shape),
		param=1, asis='list', bsize=30)
	assert v == [ (31, (30, N_SAMPLE, 2)) ] * 3 + [ (11, (10, N_SAMPLE, 2)) ]
	assert f.Apply('position', lambda x: None) is None


@pytest.mark.skipif(not has_fork, reason='no fork')
def test_run_parallel_fork(f):
	f.FilterSet2(variant=np.arange(100, 1100), verbose=False)
	pos = f.GetData('position')
	assert np.array_equal(f.RunParallel(get_position, ncpu=3), pos)
	assert np.array_equal(f.RunParallel(get_position, ncpu=1), pos)
	# a lambda function, which is not pickled
	v = f.RunParallel(lambda file, param: file.GetData('position') + param,
		param=1, ncpu=2)
	assert np.array_equal(v, pos + 1)
	v = f.RunParallel(get_sample_num, ncpu=4, split='by.sample', combine='list')
	assert len(v) == 4 and sum(sum(v, [])) == N_SAMPLE
	v = f.RunParallel(get_position, ncpu=2, combine=(lambda x, y: np.hstack([y, x])))
	assert np.array_equal(np.sort(v), np.sort(pos))
	assert f.RunParallel(get_position, ncpu=2, combine='none') is None
	# the filter is kept
	assert np.array_equal(f.GetData('position'), pos)


def test_run_parallel_pool(f):
	f.FilterSet2(variant=np.arange(100, 1100), verbose=False)
	pos = f.GetData('position')
	with mp.get_context('spawn').Pool(2) as pool:
		assert np.array_equal(f.RunParallel(get_position, ncpu=pool), pos)
		v = f.RunParallel(get_sample_num, ncpu=pool, split='by.sample', combine='list')
		assert sum(sum(v, [])) == N_SAMPLE


def test_fork_jobs_in_one_process(f):
	# a forked process may run more than one job
	f.FilterSet2(variant=np.arange(1000), verbose=False)
	pos = f.GetData('position')
	ps._fork_param = (f, get_position, None, 'by.variant')
	try:
		v = [ ps._proc_fork_func([ i, 3 ]) for i in range(3) ]
	finally:
		ps._fork_param = None
	assert np.array_equal(np.hstack(v), pos)
	assert np.array_equal(f.GetData('position'), pos)


@pytest.mark.parametrize('ncpu', [ 1, 3 ])
def test_allele_freq(f, ncpu):
	f.FilterSet2(sample=np.arange(0, N_SAMPLE, 3), variant=np.arange(500, 3500),
		verbose=False)
	g = f.GetData('genotype').reshape(3000, -1)
	af = np.column_stack([ (g == 0).mean(axis=1), (g == 1).mean(axis=1) ])
	v = f.AlleleFreq(ncpu=ncpu)
	assert v.shape == (3000,) and np.allclose(v, af[:, 0])
	assert np.allclose(f.AlleleFreq(1, ncpu), af[:, 1])
	v = f.AlleleFreq(None, ncpu)
	assert len(v) == 3000 and np.allclose(np.vstack(v), af)
	with pytest.raises(ValueError):
		f.AlleleFreq(-1)


def test_allele_freq_missing():
	g = np.array([ [[0, 1], [255, 1]], [[255, 255], [255, 255]],
		[[2, 0], [0, 1]] ], dtype=np.uint8)
	assert np.allclose(ps._calc_af_ref(0, g), [ 1/3, np.nan, 1/2 ], equal_nan=True)
	v = ps._calc_af_all(g, np.array([ 2, 2, 3 ]))
	assert np.allclose(v[0], [ 1/3, 2/3 ])
	assert np.all(np.isnan(v[1]))
	assert np.allclose(v[2], [ 1/2, 1/4, 1/4 ])
