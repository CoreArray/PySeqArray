from distutils.core import setup, Extension
import os
import pygds


src_fnlst = [ os.path.join('src', fn) for fn in [
	'LinkGDS.c' ] ]


setup(name='PySeqArray',
	version = '0.1',
	description = 'Python Interface to GDS Files for Data Management of Whole-Genome Sequence Variant Calls',
	url = 'http://github.com/CoreArray/PySeqArray',
	author = 'Xiuwen Zheng',
	author_email = 'zhengxwen@gmail.com',
	license = 'GPLv3',
	packages = [ 'PySeqArray' ],
	install_requires = [ 'numpy', 'pandas', 'pygds' ],
	ext_modules = [ Extension('PySeqArray.ccall',
		src_fnlst,
		include_dirs = [ pygds.get_include() ],
		define_macros = [ ('USING_PYTHON', None) ],
	) ]
)
