# import numpy
import numpy as np
# import pygds
import pygds



# ===========================================================================

class SeqArray(pygds.gdsfile):
	'''
	Common base class for SeqArray GDS files
	'''

	def __init__(self):
		self.filename = ''
		self.fileid = -1

	def __del__(self):
		self.close()

	def create(self, filename, allow_dup=False):
		self.fileid = cc.create_gds(filename, allow_dup)
		self.filename = filename

	def open(self, filename, readonly=True, allow_dup=False):
		self.fileid = cc.open_gds(filename, readonly, allow_dup)
		self.filename = filename

	def close(self):
		cc.close_gds(self.fileid)
		self.fileid = -1
		self.filename = ''

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




