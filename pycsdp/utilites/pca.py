"""Uses the PCA reduction method from C-Library (pca/)."""

import ctypes
import os
import numpy as np
from pycsdp.utilites import defaults

__VERSION__ = 'Release'


def nl_pca_reduction(arr, ver=__VERSION__):

    # TODO add steps for make nl_pca
    shape = arr.shape
    dtype = arr.dtype

    nl_path = defaults.get_nlpca(ver)
    nl_lib = ctypes.CDLL(nl_path)

    nl_fn = nl_lib.my_c_function
    nl_fn.argtypes = [np.ctypeslib.ndpointer(dtype=dtype, ndim=len(shape), shape=shape, flags='C_CONTIGUOUS')]
    nl_fn.restype = None

    nl_fn(arr)
    # TODO add return outputs of nl_pca
    return



def l_pca_reduction(arr, ver=__VERSION__):

    # TODO add steps for make nl_pca
    shape = arr.shape
    dtype = arr.dtype

    l_path = defaults.get_lpca(ver)
    l_lib = ctypes.CDLL(l_path)

    l_fn = l_lib.my_c_function
    l_fn.argtypes = [np.ctypeslib.ndpointer(dtype=dtype, ndim=len(shape), shape=shape, flags='C_CONTIGUOUS')]
    l_fn.restype = None

    l_fn(arr)
    # TODO add return outputs of nl_pca
    return
