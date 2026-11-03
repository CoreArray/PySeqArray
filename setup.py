import os

from setuptools import setup, Extension
import numpy
import pycorearray


src_fnlst = [ os.path.join('src', fn) for fn in [
	'GetData.cpp', 'Index.cpp', 'Methods.cpp', 'ReadByVariant.cpp',
	'PySeqArray.cpp', 'LinkGDS.c', 'vectorization.c' ] ]


# the project metadata is in pyproject.toml
ext = Extension(
	'PySeqArray.ccall',
	src_fnlst,
	include_dirs = [ pycorearray.get_include(), numpy.get_include() ],
	define_macros = [
		('USING_PYTHON', None),
		('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
		# oldest NumPy the binary runs on; keep in sync with pyproject.toml
		('NPY_TARGET_VERSION', 'NPY_1_25_API_VERSION'),
	],
)

setup(ext_modules=[ext])
