import numpy as np
from pycsdp.utilites import pca



arr = np.loadtxt('/home/nageshbansal/soc/lpb/objective_reduction/dtlz5_2_5_moea.out')


print(arr.shape)
nl_fn, free_fs = pca.pca_reduction('nl_pca', arr)
out = pca.nl_pca_reduction(arr, nl_fn, free_fs)
pca._clean_pca('nl_pca')
