import numpy as np
import sys
import subprocess
from os.path import dirname, abspath, join

# target:  is to execute both makefiles in nl_mvu_pca_dmoss and
# l_pca_dmoss_kemoss. get the output...that is a string. We need
# to convert the string output to list/ndarray in python

# Step 1: Run the makefile
# Step 2: pass .txt to executable


def make_obj_red(path):
    """To execute the makefile using subprocess

    Returns
    -------
    outfile
    """

    process = subprocess.Popen(
        f"make clean -C {path} && make -C {path}", stderr=subprocess.STDOUT,
        shell=True)

    # subprocess.Popen(
    #     "make -C {path}", stderr=subprocess.STDOUT, shell=True
    # )

    if process.wait() != 0:
        print('Something went wrong, idk what!')


def execute_model_deprecated(model, data):
    """This method shows how the older vesion works

    model: str
        Model Name (l_pca or nl_pca)
    data: str, ndarray
        Path to the data (.txt) file.

    """
    if isinstance(data, np.ndarray):
        with open('pf_data_nl.txt', 'w') as f:
            f.write(data)

    if model == 'nl_pca':
        path = 'src/pca/nl_mvu_pca_dmoss_kemoss/Debug/nl_mvu_pca_dmoss_kemoss'
    elif model == 'l_pca':
        path = 'src/pca/l_pca_dmoss_kemoss/Debug/l_pca_dmoss_kemoss'

    process = subprocess.Popen(
        f"./{path} {data}", stderr=subprocess.STDOUT,
        shell=True)
    if process.wait() != 0:
        print('Model run terminated!')
    # pass


def execute_model(model, data):
    """
    mode: str
        Model Name (l_pca or nl_pca)
    data: ndarray
        Numpy array non dominated solution.
    """
    # process = subprocess.Popen()
    pass


l_pca = make_obj_red(
    'src/pca/l_pca_dmoss_kemoss/Debug')
# nl_pca = make_obj_red(
#     'src/pca/nl_mvu_pca_dmoss_kemoss/Debug')

execute_model_deprecated("l_pca", 'src/dtlz5_2_5_moea.out')




# NOT NEEDED ----- TRY RUN
# with open('src/dtlz5_2_5_moea.out', 'r') as f:
#     input_data = f.read()

# print(input_data)
# # nd_data = input_data.split(' ')
# # print(nd_data)
# nd_array = np.array(list(input_data))

# print(nd_array)

# print(nd_array)
# TRY RUN ---- NUMPY TO CTYPE DATA CONversion
