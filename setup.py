from distutils.core import setup, Extension
import numpy
import os
import pycorearray


src_fnlst = [ os.path.join('src', fn) for fn in [
	'GetData.cpp', 'Index.cpp', 'Methods.cpp', 'ReadByVariant.cpp',
	'PySeqArray.cpp', 'LinkGDS.c', 'vectorization.c' ] ]


setup(name='PySeqArray',
	version = '0.1.0',
	description = 'Python Interface to SeqArray Files for Data Management of Whole-Genome Sequence Variant Calls',
	url = 'http://github.com/CoreArray/PySeqArray',
	author = 'Xiuwen Zheng',
	author_email = 'zhengxwen@gmail.com',
	license = 'GPLv3',
	packages = [ 'PySeqArray' ],
	install_requires = [ 'numpy', 'pycorearray' ], # 'multiprocessing' ],
	ext_modules = [ Extension('PySeqArray.ccall',
		src_fnlst,
		include_dirs = [ pycorearray.get_include(), numpy.get_include() ],
		define_macros = [
			('USING_PYTHON', None),
			('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
			# oldest NumPy the binary runs on
			('NPY_TARGET_VERSION', 'NPY_1_25_API_VERSION'),
		],
	) ],
	package_data = {
		'PySeqArray': [ 'data/*.gds' ]
	}
)
