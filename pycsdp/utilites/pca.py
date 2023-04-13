"""Uses the PCA reduction method from C-Library (pca/)."""

import ctypes
import os
from pycsdp.utilites import defaults, make
import time

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

    shape = arr.shape
    # dtype = arr.dtype
    nl_lib = ctypes.CDLL(os.path.join(path, 'nl_mvu_pca_dmoss_kemoss.so'))

    nl_fn = nl_lib.main

    c_double_p = ctypes.POINTER(ctypes.c_double)
    c_int_p = ctypes.POINTER(ctypes.c_int)
    # data = arr.astype(np.float32)
    data_p = arr.ctypes.data_as(c_double_p)
    nl_fn.argtypes = (c_double_p, ctypes.c_int, ctypes.c_int)
    nl_fn.restypes = c_int_p

    fs_ptr = nl_fn(data_p, shape[0], shape[1])
    fs_data = np.ctypeslib.as_array(fs_ptr, shape=(2,))
    return fs_data



def l_pca_reduction(arr, path, ver=__VERSION__):

    shape = arr.shape
    # dtype = arr.dtype
    l_lib = ctypes.CDLL(os.path.join(path, 'l_pca_dmoss_kemoss.so'))

    l_fn = l_lib.main
    free_fs = l_lib.free_fs

    c_double_p = ctypes.POINTER(ctypes.c_double)
    c_int_p = ctypes.POINTER(ctypes.c_int)

    data_p = arr.ctypes.data_as(c_double_p)
    l_fn.argtypes = (c_double_p, ctypes.c_int, ctypes.c_int)
    l_fn.restype = c_int_p

    free_fs.argtype = ctypes.c_void_p
    free_fs.restype = ctypes.c_int

    fs_ptr = l_fn(data_p, shape[0], shape[1])
    fs_data = np.ctypeslib.as_array(fs_ptr, shape=(2,))

    # Free fs pointer (causes seg fault)
    # time.sleep(2)
    free_fs(fs_ptr)

    return fs_data
