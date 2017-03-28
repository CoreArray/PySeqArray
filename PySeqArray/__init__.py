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
		id = self.fileid
		pygds.gdsfile.close(self)
		cc.file_done(id)

	def sync(self):
		cc.sync_gds(self.fileid)

	def filesize(self):
		return cc.filesize(self.fileid)

	def root(self):
		v = gdsnode()
		(v.idx, v.pid) = cc.root_gds(self.fileid)
		return v

	def index(self, path, silent=False):
		v = gdsnode()
		v.idx, v.pid = cc.index_gds(self.fileid, path, silent)
		if v.idx >= 0:
			return v
		else:
			return None

	def show(self):
		print('File:', self.filename)
		self.root().show()




