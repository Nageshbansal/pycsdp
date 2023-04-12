import numpy as np
import subprocess


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

    stdout, stderr = process.communicate()

    if process.wait() != 0:
        print('Something went wrong, idk what!')
        return stdout, stderr

    return


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
