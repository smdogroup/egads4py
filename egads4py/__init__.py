'''
egads4py import
'''

import os

def get_cython_include():
    '''
    Get the include directory for the Cython .pxd files in TACS
    '''
    return [os.path.abspath(os.path.dirname(__file__))]

def get_include():
    '''
    Get the include directory for the Cython .pxd files in TACS
    '''
    root_path, tail = os.path.split(os.path.abspath(os.path.dirname(__file__)))

    rel_inc_dirs = ['include']

    inc_dirs = []
    for path in rel_inc_dirs:
    	inc_dirs.append(os.path.join(root_path, path))

    return inc_dirs

def get_libraries():
    '''
    Get the library directories
    '''
    root_path, tail = os.path.split(os.path.abspath(os.path.dirname(__file__)))

    rel_lib_dirs = ['lib']
    libs = ['egads']
    lib_dirs = []
    for path in rel_lib_dirs:
    	lib_dirs.append(os.path.join(root_path, path))

    return lib_dirs, libs
