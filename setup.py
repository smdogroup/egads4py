import os
from glob import glob
from subprocess import check_output

# Import numpy
import numpy

# Import distutils
from setuptools import setup
from distutils.core import Extension as Ext
from Cython.Build import cythonize

# Convert from local to absolute directories
def get_global_dir(files):
    tmr_root = os.path.abspath(os.path.dirname(__file__))
    new = []
    for f in files:
        new.append(os.path.join(tmr_root, f))
    return new

# Relative paths for the include/library directories
rel_inc_dirs = ['include', 'src', 'src/Surreal']
rel_lib_dirs = ['lib']
libs = ['egads']
sources = ['egads4py/egads.pyx']

# Convert from relative to absolute directories
inc_dirs = get_global_dir(rel_inc_dirs)

lib_dirs = [os.path.join(os.environ['CASROOT'],
                         os.environ['CASARCH'], 'lib')]
runtime_lib_dirs = [os.path.join(os.environ['CASROOT'],
                                 os.environ['CASARCH'], 'lib')]

# Add the numpy directories
inc_dirs.extend([numpy.get_include()])
lib_dirs.extend(get_global_dir(rel_lib_dirs))
runtime_lib_dirs.extend(get_global_dir(rel_lib_dirs))

# Add the include directories from OpenCascade
for sufix in ['include/oce', 'inc', 'include']:
    cas_inc = os.path.join(os.environ['CASROOT'], sufix)
    if os.path.isdir(cas_inc):
        inc_dirs.append(cas_inc)
        break

# Add the libraries from OpenCascade
libs.extend(['TKBool', 'TKernel', 'TKFeat', 'TKBO', 'TKGeomAlgo',
             'TKMath', 'TKOffset', 'TKPrim', 'TKPShape', 'TKTopAlgo',
             'TKBRep', 'TKG2d', 'TKG3d', 'TKGeomBase', 'TKShHealing',
             'TKSTEP', 'TKSTEP209', 'TKSTEPBase', 'TKSTEPAttr',
             'TKXSBase', 'TKIGES', 'TKFillet', 'PTKernel', 'dl' ])

exts = []

exts.append(Ext('egads4py.egads', sources=sources,
                include_dirs=inc_dirs, libraries=libs,
                library_dirs=lib_dirs,
                runtime_library_dirs=runtime_lib_dirs))
setup(name='egads4py',
      version=0.1,
      description='python interface to egads',
      author='Graeme J. Kennedy',
      author_email='graeme.kennedy@ae.gatech.edu',
      ext_modules=cythonize(exts, language='c',
                            include_path=inc_dirs))
