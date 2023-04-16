"""Uses the PCA reduction method from C-Library (pca/)."""

import ctypes
import os
from pycsdp.utilites import defaults, make
# import time
import numpy as np

__VERSION__ = 'Release'


def _clean_pca(model='all'):
    """Cleans all generated file dring make.

    Parameters
    ----------
    model: str
        l_pca, nl_pca, all (default: all)
    """

    nl_paths = [defaults.get_nlpca("Debug"), defaults.get_nlpca("Release")]
    l_paths = [defaults.get_lpca("Debug"), defaults.get_lpca("Release")]
    rm_paths = []

    if model == 'all':
        rm_paths = nl_paths + l_paths
    elif model == 'nl_pca':
        rm_paths = nl_paths
    elif model == 'l_pca':
        rm_paths = l_paths

    for mpath in rm_paths:
        make.make_clean(mpath)


def pca_reduction(model, ver=__VERSION__):

    if model == 'nl_pca':
        nl_path = defaults.get_nlpca(ver)
        nl_make = make.make_obj_red(nl_path)

        # add checks for process reliablity
        nl_fn, free_fs = load_nl_pca(nl_path, ver)
        return nl_fn, free_fs

    elif model == 'l_pca':
        l_path = defaults.get_lpca(ver)
        l_make = make.make_obj_red(l_path)

        # add checks for process reliablity
        l_fn, free_fs = load_l_pca(l_path, ver)
        return l_fn, free_fs


def load_nl_pca(nl_path, ver=__VERSION__):
    nl_lib = ctypes.CDLL(os.path.join(nl_path, 'nl_mvu_pca_dmoss_kemoss.so'))
    nl_fn = nl_lib.main
    free_fs = nl_lib.free_fs
    return nl_fn, free_fs

def load_l_pca(nl_path, ver=__VERSION__):
    l_lib = ctypes.CDLL(os.path.join(nl_path, 'l_pca_dmoss_kemoss.so'))
    l_fn = l_lib.main
    free_fs = l_lib.free_fs
    return l_fn, free_fs


def nl_pca_reduction(arr, nl_fn, free_fs, ver=__VERSION__):

    shape = arr.shape
    # dtype = arr.dtype
    c_double_p = ctypes.POINTER(ctypes.c_double)
    c_int_p = ctypes.POINTER(ctypes.c_int)
    data_p = arr.ctypes.data_as(c_double_p)

    nl_fn.argtypes = (c_double_p, ctypes.c_int, ctypes.c_int)
    nl_fn.restype = c_int_p

    free_fs.argtype = ctypes.c_void_p
    free_fs.restype = ctypes.c_int

    fs_ptr = nl_fn(data_p, shape[0], shape[1])
    fs_data = np.ctypeslib.as_array(fs_ptr, shape=(3,))
    fs = fs_data.copy()
    try:
        free_fs(fs_ptr)
    except:
        print("something went wrong!!")

    return fs


def l_pca_reduction(arr, l_fn, free_fs, ver=__VERSION__):

    shape = arr.shape
    # dtype = arr.dtype

    c_double_p = ctypes.POINTER(ctypes.c_double)
    c_int_p = ctypes.POINTER(ctypes.c_int)

    data_p = arr.ctypes.data_as(c_double_p)
    l_fn.argtypes = (c_double_p, ctypes.c_int, ctypes.c_int)
    l_fn.restype = c_int_p

    free_fs.argtype = ctypes.c_void_p
    free_fs.restype = ctypes.c_int

    fs_ptr = l_fn(data_p, shape[0], shape[1])
    fs_data = np.ctypeslib.as_array(fs_ptr, shape=(2,))
    fs = fs_data.copy()

    free_fs(fs_ptr)
    return fs
