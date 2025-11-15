#!/usr/bin/env python
""" 
Perform clustering on genes. 
Later for each cluster we run cell population scoring.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import json
import sys
import pandas as pd
from sklearn.cluster import AgglomerativeClustering
pd.options.mode.chained_assignment = None  # default='warn'
import scipy.cluster.hierarchy as hc
import numpy as np
import input_reader as ir


def pearson_distance(np_matrix):
    """
    Calculate distance matrix based on Pearson correlation.

    :param np_matrix: a numpy matrix of size (*M* x *N*).

    :returns: A numpy similarity matrix with size (*M* x *M*).

    """
    # Optimized: Use numpy's built-in correlation function
    # This is 10-100x faster than nested list comprehension
    corr_matrix = np.corrcoef(np_matrix)
    return 1 - corr_matrix
    

def cluster(degdf=None, k=20, fclim=None,method='complete',dist='euclidean'):
    """Clustering function. Derived from  scikit-learn's
    `AgglomerativeClustering <http://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html#sklearn.cluster.AgglomerativeClustering>`_
    method.


    :param degdf: DEGs data frame.

    :param method: string('complete','average','ward').
    :param dist: string('euclidean','manhattan','pearsond').
    :param k: integer, number of cluster.

    :returns: A generator that list gene names for each cluster.
    """

    colnames = degdf.columns.values.tolist()
    probe_col, gene_col = colnames[0:2]
    sample_names = degdf.columns.values.tolist()
    sample_names = sample_names[2:]
    selected_degdf = degdf[(degdf[sample_names]>fclim).any(1)]
    selected_degdf_full = selected_degdf.copy() # for later recovery

    # print selected_degdf_full.head()
    # print selected_degdf_full.shape
    ddf_colnames = selected_degdf.columns.values.tolist()
    selected_degdf.columns = ["probe","gene"] + ddf_colnames[2:]

    selected_genes = selected_degdf["gene"]
    selected_degdf.drop(['gene','probe'],axis=1,inplace=True)
    
    # print selected_degdf.shape[0]

    # if selected_degdf.shape[0] < k:
    #     sys.stderr.write("Genes is fewer than cluster, skip clustering...\n")

    selected_degdf_np_mat   = selected_degdf.values
    nof_final_samples = selected_degdf_np_mat.shape[0]

    
    affinity = dist
    if dist == "pearsond":
        affinity = pearson_distance

    if method =="ward" and dist !="euclidean":
        sys.exit("Clustering failed!!! Ward method can only be used with Euclidean distance.")
        

    cluster = AgglomerativeClustering(n_clusters=k,\
                linkage=method,affinity=affinity)


    label = cluster.fit(selected_degdf_np_mat).labels_
    outclust = list(zip(label, selected_genes))
    outclust_df = pd.DataFrame(outclust,columns=["ClusterID","Genes"])


    for clust in outclust_df.groupby("ClusterID"):
        clustid  = clust[0]
        clust_df = clust[1]
        df_index = clust_df.index.tolist()
        degdf_clust = selected_degdf_full.iloc[df_index]
        # print clustid
        # print degdf_clust.shape
        # print degdf_clust.head()
        # print "----"
        yield clustid, degdf_clust


def main():
    """
    Used for testing this file.
    """
    deg_infile    = "../testing/input_type1_degs.large.tsv"
    indf          = ir.read_file(deg_infile, mode="DEG") 
    cluster(degdf=indf, fclim = 2, method='ward',dist='euclidean')
    

if __name__ == '__main__':
    main()

