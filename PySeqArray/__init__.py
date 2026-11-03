# import numpy
import numpy as np
# import os
import os
# import multiprocessing
import multiprocessing as mp
import multiprocessing.pool as pl
# import pycorearray
import pycorearray
# import c library
import PySeqArray.ccall as cc
# other ...
from functools import reduce


## export version number
__version__ = '0.2.0'



# ===========================================================================

def seqExample(filename=None):
	"""Example files

	Return a file name in the folder of example data.

	Parameters
	----------
	filename : str
		a file name in the folder of example data, or None for returning the path of example folder

	Returns
	-------
	string

	Examples
	--------
	>>> seqExample('1KG_phase1_release_v3_chr22.gds')
	"""
	import PySeqArray
	s = os.path.dirname(PySeqArray.__file__)
	if filename is None:
		return os.path.join(s, 'data')
	else:
		return os.path.join(s, 'data', filename)



# ===========================================================================

# the file, user-defined function, parameter and split passed to the forked
# processes, which inherit them (and the opened file) from the parent process
_fork_param = None

# define internal function using forking
def _proc_fork_func(x):
	i = x[0]; ncpu = x[1]
	file, fun, param, split = _fork_param
	# a process may run more than one job, so restore the filter afterward
	cc.flt_push(file.fileid, False)
	try:
		cc.flt_split(file.fileid, i, ncpu, split)
		return fun(file, param)
	finally:
		cc.flt_pop(file.fileid)

# define a process function
def _proc_func(x):
	i = x[0]; ncpu = x[1]
	fn = x[2]; fun = x[3]; param = x[4]; sel = x[5]; split = x[6]
	import PySeqArray
	import PySeqArray.ccall as cc
	file = PySeqArray.SeqArrayFile()
	file.open(fn, allow_dup=True)
	file.FilterSet2(sel[0], sel[1], verbose=False)
	cc.flt_split(file.fileid, i, ncpu, split)
	return fun(file, param)



# ===========================================================================

# frequencies of the allele 'ref' per variant in a block of genotypes
# (variant x sample x ploidy, 0xFF for missing)
def _calc_af_ref(ref, geno):
	g = geno.reshape(geno.shape[0], -1)
	n = np.sum(g != 0xFF, axis=1)
	with np.errstate(invalid='ignore'):
		return np.sum(g == ref, axis=1) / n

# frequencies of all alleles per variant in a block of genotypes
def _calc_af_all(geno, num_allele):
	rv = []
	with np.errstate(invalid='ignore'):
		for g, n in zip(geno.reshape(geno.shape[0], -1), num_allele):
			g = g[g != 0xFF]
			rv.append(np.bincount(g, minlength=n)[:n] / g.size)
	return rv

# allele frequencies in a process of RunParallel()
def _proc_af(file, param):
	ref = param[0]; verbose = param[1]
	if ref is None:
		v = file.Apply(['genotype', '$num_allele'], _calc_af_all, asis='list',
			verbose=verbose)
		return [ x for blk in v for x in blk ]
	else:
		return file.Apply('genotype', _calc_af_ref, param=ref, asis='unlist',
			verbose=verbose)



# ===========================================================================

class SeqArrayFile(pycorearray.gdsfile):
	"""
	Class for SeqArray GDS files
	"""

	def __init__(self):
		pycorearray.gdsfile.__init__(self)

	def __del__(self):
		cc.file_done(self.fileid)
		pycorearray.gdsfile.__del__(self)


	def create(self, filename, allow_dup=False):
		raise Exception('not supported!')


	def open(self, filename, readonly=True, allow_dup=False):
		"""Open an SeqArray file

		Open an existing file of SeqArray GDS for reading or writing.

		Parameters
		----------
		filename : str
			the file name of a new GDS file to be created
		readonly : bool
			if True, the file is opened read-only; otherwise, it is allowed to write data to the file
		allow_dup : bool
			if True, it is allowed to open a GDS file with read-only mode when it has been opened in the same session

		Returns
		-------
		None

		See Also
		--------
		close: close a SeqArray file
		"""
		pycorearray.gdsfile.open(self, filename, readonly, allow_dup)
		cc.file_init(self.fileid)
		# TODO: file checking


	def close(self):
		"""Close a SeqArray file

		Close a SeqArray GDS file.

		Returns
		-------
		None

		See Also
		--------
		open : open an existing SeqArray file
		"""
		cc.file_done(self.fileid)
		pycorearray.gdsfile.close(self)


	def FilterSet(self, sample_id=None, variant_id=None, intersect=False, verbose=True):
		"""Set a filter

		Set a filter to sample and/or variant with IDs.

		Parameters
		----------
		sample_id : str
			sample id to be selected
		variant_id : bool
			variant id to be selected
		intersect : bool
			if False, the candidate variants for selection are all possible variants (by default);
			if True, the candidate variants are from the selected variants defined via the previous call
		verbose : bool
			if True, show information

		Returns
		-------
		None

		See Also
		--------
		FilterReset : reset the filter
		"""
		s = not sample_id is None
		v = not variant_id is None
		if s or v:
			if s:
				cc.set_sample(self.fileid, sample_id, intersect, verbose)
			if v:
				cc.set_variant(self.fileid, variant_id, intersect, verbose)


	def FilterSet2(self, sample=None, variant=None, intersect=False, verbose=True):
		"""Set a filter

		Set a filter to sample and/or variant with a bool vector or an index vector.

		Parameters
		----------
		sample : vector, range
			a bool vector, an indexing vecot for selecting samples or a range object
		variant : vector, range
			a bool vector, an indexing vecot for selecting variants or a range object
		intersect : bool
			if False, the candidate variants for selection are all possible variants (by default);
			if True, the candidate variants are from the selected variants defined via the previous call
		verbose : bool
			if True, show information

		Returns
		-------
		None

		See Also
		--------
		FilterSet : set a filter
		FilterReset : reset the filter
		"""
		# the C code requires a contiguous numpy array
		if not sample is None:
			sample = np.ascontiguousarray(sample)
			cc.set_sample2(self.fileid, sample, intersect, verbose)
		if not variant is None:
			variant = np.ascontiguousarray(variant)
			cc.set_variant2(self.fileid, variant, intersect, verbose)


	def FilterReset(self, sample=True, variant=True, verbose=True):
		"""Reset the filter

		Clear the existing filter on sample and/or variant.

		Parameters
		----------
		sample : str
			if True, reset the filter of sample
		variant : bool
			if True, reset the filter of variant
		verbose : bool
			if True, show information

		Returns
		-------
		None

		See Also
		--------
		FilterSet : set a filter
		"""
		if sample:
			cc.set_sample(self.fileid, None, False, verbose)
		if variant:
			cc.set_variant(self.fileid, None, False, verbose)


	def FilterPush(self, reset=True):
		"""Push a filter

		Push the current filter to the stack.

		Parameters
		----------
		reset : bool
			if True, reset the filter of sample and variant after pushing

		Returns
		-------
		None

		See Also
		--------
		FilterPop : recover the last filter
		"""
		cc.flt_push(self.fileid, reset)


	def FilterPop(self):
		"""Pop a filter

		Pop or recover the last filter in the stack.

		Returns
		-------
		None

		See Also
		--------
		FilterPush : push the current filter to the stack
		"""
		cc.flt_pop(self.fileid)


	def FilterGet(self, sample=True):
		"""Get a sample/variant filter

		Get a sample or variant filter.

		Parameters
		----------
		sample : bool
			If True, return the sample filter; otherwise, return the variant filter

		Returns
		-------
		A numpy object (a bool vector)

		See Also
		--------
		FilterSet : set a filter
		"""
		return(cc.get_filter(self.fileid, sample))


	def GetData(self, name):
		"""Get data

		Get data from a SeqArray file with a given variable name and a sample/variant filter

		Parameters
		----------
		name : str
			the variable name

		Returns
		-------
		a numpy array object

		See Also
		--------
		FilterSet : set a filter
		"""
		return cc.get_data(self.fileid, name)


	def Apply(self, name, fun, param=None, asis='none', bsize=1024, verbose=False):
		"""Apply function over array margins

		Apply a user-defined function to margins of genotypes and annotations via blocking

		Parameters
		----------
		name : str, list
			the variable name, or a list of variable names
		fun : function
			the user-defined function
		param: object
			the parameter passed to the user-defined function if it is not None
		asis : str
			'none', no return; 'list', a list of the returned values from the user-defined function;
			'unlist', flatten the returned values from the user-defined function
		bsize : int
			block size
		verbose : bool
			show progress information if True

		Returns
		-------
		None, a list or a numpy array object

		See Also
		--------
		FilterSet : set a filter
		"""
		v = cc.apply(self.fileid, name, fun, param, asis, bsize, verbose)
		if asis == 'unlist':
			v = np.hstack(v)
		return(v)


	def RunParallel(self, fun, param=None, ncpu=0, split='by.variant', combine='unlist'):
		"""Apply Functions in Parallel

		Apply a user-defined function in parallel over array margins

		Parameters
		----------
		fun : function
			the user-defined function
		param : object
			the parameter passed to the user-defined function if it is not None
		ncpu : int
			the number of cores or an instance of 'multiprocessing.pool.Pool';
			0 to use the number of cores minus 1
		split : str
			'by.variant', 'by.sample', 'none': split the dataset by variant or sample according to multiple processes, or "none" for no split
		combine : str, function
			'none', no return; 'list', a list of the returned values from the user-defined function;
			'unlist', flatten the returned values from the user-defined function

		Returns
		-------
		None, a list or a numpy array object

		Notes
		-----
		If `ncpu` is a number, the processes are forked on Linux and macOS,
		and share the opened file and the filter. Otherwise (a Pool, or on
		Windows), each process opens the file again, so `fun` and `param`
		should be picklable (e.g., `fun` is not a lambda function).
		"""
		# check
		if not isinstance(ncpu, (int, float, pl.Pool)):
			raise ValueError('`ncpu` should be a numeric value or `multiprocessing.pool.Pool`.')
		if not (combine is None or isinstance(combine, str) or callable(combine)):
			raise ValueError('`combine` should be None, a string or a function.')
		# the number of processes
		if isinstance(ncpu, pl.Pool):
			pool = ncpu
			ncpu = pool._processes
		else:
			pool = None
			ncpu = int(ncpu)
			if ncpu <= 0:
				ncpu = max(mp.cpu_count() - 1, 1)
		# run
		if ncpu > 1:
			# fork explicitly, since it is not the default start method on
			# macOS (Python >= 3.8) or Linux (Python >= 3.14)
			if pool is None and 'fork' in mp.get_all_start_methods():
				global _fork_param
				_fork_param = (self, fun, param, split)
				try:
					with mp.get_context('fork').Pool(processes=ncpu) as pa:
						v = pa.map(_proc_fork_func, [ [ i,ncpu ] for i in range(ncpu) ])
				finally:
					_fork_param = None
			else:
				sel = [ self.FilterGet(True), self.FilterGet(False) ]
				pm = [ [ i,ncpu,self.filename,fun,param,sel,split ] for i in range(ncpu) ]
				if pool is None:
					with pl.Pool(processes=ncpu) as pa:
						v = pa.map(_proc_func, pm)
				else:
					v = pool.map(_proc_func, pm)
			# output
			if combine is None or combine == 'none':
				v = None
			elif combine == 'unlist':
				v = np.hstack(v)
			elif callable(combine):
				v = reduce(combine, v)
			elif combine != 'list':
				raise ValueError('`combine` is invalid.')
			return v
		else:
			return fun(self, param)


	####  Methods  ####

	def AlleleFreq(self, ref=0, ncpu=1, verbose=False):
		"""Allele frequencies

		Calculate the allele frequencies of the selected samples per selected variant

		Parameters
		----------
		ref : int, None
			the index of allele (0 for the reference allele, 1 for the first alternative allele, ...),
			or None for the frequencies of all alleles
		ncpu : int
			the number of cores or an instance of 'multiprocessing.pool.Pool', see RunParallel()
		verbose : bool
			show progress information if True

		Returns
		-------
		a numpy array if `ref` is an integer, or a list of numpy arrays if `ref` is None;
		NaN if all genotypes of a variant are missing

		See Also
		--------
		RunParallel : apply functions in parallel
		"""
		if not (ref is None or (isinstance(ref, (int, np.integer)) and ref >= 0)):
			raise ValueError('`ref` should be a non-negative integer or None.')
		if ref is None:
			return self.RunParallel(_proc_af, param=[ ref, verbose ], ncpu=ncpu,
				combine=(lambda x, y: x + y))
		else:
			return self.RunParallel(_proc_af, param=[ ref, verbose ], ncpu=ncpu)

