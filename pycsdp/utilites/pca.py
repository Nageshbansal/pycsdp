"""Uses the PCA reduction method from C-Library (pca/)."""

import ctypes
import os
import numpy as np
from pycsdp.utilites import defaults, make

__VERSION__ = 'Release'


def pca_reduction(model, arr, ver=__VERSION__):

    if model == 'nl_pca':
        nl_path = defaults.get_nlpca(ver)
        nl_make = make.make_obj_red(nl_path)

        # add checks for process reliablity

        f_arr = nl_pca_reduction(arr, nl_path)
        return f_arr

    elif model == 'l_pca':
        l_path = defaults.get_lpca(ver)
        l_make = make.make_obj_red(l_path)

        # add checks for process reliablity
        print(l_path)
        f_arr = l_pca_reduction(arr, l_path)
        return f_arr


def nl_pca_reduction(arr, path, ver=__VERSION__):

    # TODO add steps for make nl_pca
    shape = arr.shape
    dtype = arr.dtype
    nl_lib = ctypes.CDLL(os.path.join(path, 'nl_mvu_pca_dmoss_kemoss.so'))

    nl_fn = nl_lib.my_c_function
    nl_fn.argtypes = [np.ctypeslib.ndpointer(dtype=dtype, ndim=len(shape), shape=shape, flags='C_CONTIGUOUS')]
    nl_fn.restype = None

    nl_fn(arr)
    # TODO add return outputs of nl_pca
    return


def l_pca_reduction(arr, path, ver=__VERSION__):

    # TODO add steps for make nl_pca
    shape = arr.shape
    dtype = arr.dtype
    print(path)
    l_lib = ctypes.CDLL('/home/nageshbansal/soc/lpb/pycsdp/pycsdp/pca/l_pca_dmoss_kemoss/Release/l_pca_dmoss_kemoss.so')

    l_fn = l_lib.my_c_function
    l_fn.argtypes = [np.ctypeslib.ndpointer(dtype=dtype, ndim=len(shape), shape=shape, flags='C_CONTIGUOUS')]
    l_fn.restype = None

    l_fn(arr)
    # TODO add return outputs of nl_pca
    return
