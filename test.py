import numpy as np
from pycsdp.utilites import pca

arr = np.loadtxt('/home/nageshbansal/soc/lpb/objective_reduction/dtlz5_2_5_moea.out')

print(arr.shape)

nl_pca = pca.nl_pca_reduction(arr)
