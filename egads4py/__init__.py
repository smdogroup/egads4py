'''
egads4py import
'''

import os

def get_cython_include():
    '''
    Get the include directory for the Cython .pxd files in TACS
    '''
    return [os.path.abspath(os.path.dirname(__file__))]

