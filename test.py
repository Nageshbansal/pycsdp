import numpy as np
from pycsdp.utilites import pca



arr = np.loadtxt('/home/shivam/Programming/skip-the-beats/src/dtlz5_2_5_moea.out')


print(arr.shape)

# nl_pca = pca.pca_reduction('nl_pca', arr)
# nl_pca_debug = pca.pca_reduction('nl_pca', arr, 'Debug')
# l_pca = pca.pca_reduction('l_pca', arr)
# print(l_pca)
# l_pca_debug = pca.pca_reduction('nl_pca', arr, 'Debug')

pca._clean_pca()
