import numpy as np
from pycsdp.utilites import pca


arr = np.loadtxt('/home/nageshbansal/soc/lpb/objective_reduction/dtlz5_2_5_moea.out')

print(arr.shape)

nl_pca = pca.pca_reduction('nl_pca', arr)
# nl_pca_debug = pca.pca_reduction('nl_pca', arr, 'Debug')
# l_pca = pca.pca_reduction('l_pca', arr)
# l_pca_debug = pca.pca_reduction('nl_pca', arr, 'Debug')
