# import numpy
import numpy as np
# import os
import os
# import multiprocessing
import multiprocessing as mp
import multiprocessing.pool as pl
# import pygds
import pygds
# import c library
import PySeqArray.ccall as cc
# other ...
from sys import platform
from functools import reduce


## export version number
__version__ = '0.1'



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

# define internal function using forking
def _proc_fork_func(x):
	i = x[0]; ncpu = x[1]
	file = x[2]; fun = x[3]; param = x[4]; split = x[5]
	cc.flt_split(file.fileid, i, ncpu, split)
	return fun(file, param)

# define a process function
def _proc_func(x):
	i = x[0]; ncpu = x[1]
	fn = x[2]; fun = x[3]; param = x[4]; sel = x[5]; split = x[6]
	import PySeqArray
	file = PySeqArray.SeqArrayFile()
	file.open(fn, allow_dup=True)
	file.FilterSet2(sel[0], sel[1], verbose=False)
	cc.flt_split(file.fileid, i, ncpu, split)
	return fun(file, param)



# ===========================================================================

class SeqArrayFile(pygds.gdsfile):
	"""
	Class for SeqArray GDS files
	"""

	def __init__(self):
		pygds.gdsfile.__init__(self)

	def __del__(self):
		cc.file_done(self.fileid)
		pygds.gdsfile.__del__(self)


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
		pygds.gdsfile.open(self, filename, readonly, allow_dup)
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
		pygds.gdsfile.close(self)


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
		if not sample is None:
			cc.set_sample2(self.fileid, sample, intersect, verbose)
		if not variant is None:
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


	def FilterPush(reset=True):
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


	def FilterPop():
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


	def Apply(self, name, fun, param=None, as_is='none', bsize=1024, verbose=False):
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
		as_is : str
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
		v = cc.apply(self.fileid, name, fun, param, as_is, bsize, verbose)
		if as_is == 'unlist':
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
		"""
		# check
		if not isinstance(ncpu, (int, float, pl.Pool)):
			raise ValueError('`ncpu` should be a numeric value or `multiprocessing.pool.Pool`.')
		if not (combine is None or isinstance(combine, str) or callable(combine)):
			raise ValueError('`combine` should be None, a string or a function.')
		# run
		if isinstance(ncpu, (int, float)):
			if ncpu <= 0:
				ncpu = mp.cpu_count() - 1
				if ncpu <= 0:
					ncpu = 1
			if ncpu >= 1:
				pa = pl.Pool(processes=ncpu)
		else:
			pa = ncpu
			ncpu = pa._processes
		# run
		if ncpu >= 1:
			# direct forking or not
			is_fork = False
			if isinstance(ncpu, (int, float)):
				is_fork = (platform=="linux" or platform=="linux2" or
					platform=="unix" or platform=="darwin")
			if is_fork:
				pm = [ [ i,ncpu,self,fun,param,split ] for i in range(ncpu) ]
				v = pa.map(_proc_fork_func, pm)
			else:
				sel = [ self.FilterGet(True), self.FilterGet(False) ]
				pm = [ [ i,ncpu,self.filename,fun,param,sel,split ] for i in range(ncpu) ]
				v = pa.map(_proc_func, pm)
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

