=======================================
pycsdp
=======================================


A Python Implementation  for integrating the CSDP and LHFiD

Installing the Current Release
------------------------------

If you have Python installed you can install the current release using either pip: ::

   pip install git+https://github.com/Nageshbansal/pycsdp.git


Installing the package from source
----------------------------------

1. Get the latest source by cloning this repo: ::

      git clone https://github.com/Nageshbansal/pycsdp.git

2. Install the dependencies: ::

      pip install -r requirements.txt

3. Install neonwranglerpy: ::

      pip install .


Quick Start Guide
----------------------------------
1. Import the pca module: ::

      from pycsdp.utilites import pca

2. Use the pca_reduction function to build the C libraries and import pca function: ::

      nl_fn, free_fs = pca.pca_reduction('nl_pca', arr)

3. Use the nl_pca or l_pca_reduction function to use the NL_PCA_MVU or L_PCA reduction: ::

      out = pca.nl_pca_reduction(arr, nl_fn, free_fs)
