# import numpy
import numpy as np
# import pygds
import pygds
# import c library
import PySeqArray.ccall as cc



# ===========================================================================

class SeqArray(pygds.gdsfile):
	'''
	Class for SeqArray GDS files
	'''

	def __init__(self):
		pygds.gdsfile.__init__(self)

	def __del__(self):
		cc.file_done(self.fileid)
		pygds.gdsfile.__del__(self)

	def create(self, filename, allow_dup=False):
		raise Exception('not supported!')

	def open(self, filename, readonly=True, allow_dup=False):
		pygds.gdsfile.open(self, filename, readonly, allow_dup)
		cc.file_init(self.fileid)
		# TODO: file checking

	def close(self):
		cc.file_done(self.fileid)
		pygds.gdsfile.close(self)

	def GetData(self, name):
		return cc.get_data(self.fileid, name)



