# import numpy
import numpy as np
# import os
import os
# import pygds
import pygds
# import c library
import PySeqArray.ccall as cc




# ===========================================================================

def get_example_path(filename=None):
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
	>>> get_example_path('1KG_phase1_release_v3_chr22.gds')
	"""
	import PySeqArray
	s = os.path.dirname(PySeqArray.__file__)
	if filename == None:
		return os.path.join(s, 'data')
	else:
		return os.path.join(s, 'data', filename)


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
		if sample_id!=None or variant_id!=None:
			if sample_id!=None:
				cc.set_sample(self.fileid, sample_id, intersect, verbose)
			if variant_id!=None:
				cc.set_variant(self.fileid, variant_id, intersect, verbose)

	def FilterReset(self, sample=True, variant=True, verbose=True):
		if sample:
			cc.set_sample(self.fileid, None, False, verbose)
		if variant:
			cc.set_variant(self.fileid, None, False, verbose)

	def FilterPush(reset=True):
		"""
		"""
		cc.flt_push(self.fileid, reset)

	def FilterPop():
		cc.flt_pop(self.fileid)

	def GetData(self, name):
		return cc.get_data(self.fileid, name)



