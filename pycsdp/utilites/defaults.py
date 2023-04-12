import os
from pycsdp import get_pca


def get_nlpca(ver="Release"):
    nl_path = get_pca("nl_mvu_pca_dmoss_kemoss")
    return os.path.join(nl_path, ver, "nl_mvu_pca_dmoss_kemoss.so")

def get_lpca(ver="Release"):
    l_path = get_pca("l_pca_dmoss_kemoos")
    return os.path.join(nl_path, ver, "l_pca_dmoss_kemoss.so")
