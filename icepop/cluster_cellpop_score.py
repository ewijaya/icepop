#!/usr/bin/env python
""" 
Here we take clustered genes (``degdf``), then 
for each cluster compute cell population score.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import pandas as pd
from sklearn.preprocessing import scale
import input_reader as ir
import cellpop_score as cp
import species_cell_proportion as scp
import clustering as cl

def cluster_cellpop_score(cellpopdf=None, degdf=None, fclim=None,\
        gene_count=False,logscale=False, k=None,\
        method="ward", dist="euclidean"):
    """
    :param cellpopdf: Cell population data frame.
    :param degdf: DEGs data frame.

    :param fclim: float, foldchange lower limit thresold.
    :param logscale: boolean, transform fold change to logscale.
                  This does not affect the result substansially.
    :param gene_count: boolean, normalization method. 
                  If *True*, divide by sum of product. The total weight
                  will have to result to 1.00.
    :param method: string('complete','average','ward')
    :param dist: string('euclidean','manhattan','pearsond')
    :param k: integer, number of cluster

    :returns: Pandas data frame that contains cell population scores 
            for each cluster and clustered genes with ClusterID added.

    Usage:

    >>> from icepop import cluster_cellpop_score as ccp
    >>> fold_change_lim = 1.5
    >>> nof_clust       = 15
    >>> full_clust_cpopdf, full_clust_degdf =
    >>> ccp.cluster_cellpop_score(cellpopdf=cellpop_df,
            >>> degdf=indf,fclim=fold_change_lim,
            >>> gene_count=False,logscale=False, k=nof_clust, method="ward",
            >>> dist="euclidean")
    >>> full_clust_cpopdf.to_csv("full_clust.tsv",sep="\t",index=False,float_format='%.3f')
    
    """
    allclusters_degs_df = []
    allclusters_cellpop_df = []
    for clustid, degdf_clust in cl.cluster(degdf=degdf, k=k, fclim=fclim,\
            method=method,dist=dist):

        clustid += 1
        tmp_degdf_clust = degdf_clust.copy()
        tmp_degdf_clust.insert(0,"ClusterID",clustid)
        allclusters_degs_df.append(tmp_degdf_clust)

        _, cpop_score_df = cp.deg_cellpopscore_df(cellpopdf, degdf_clust, \
                fclim=fclim, gene_count=False, logscale=False)

        # Here we scale across cell type, using Z-score
        cpop_score_scale_df = pd.DataFrame(scale(cpop_score_df,axis=0),index=cpop_score_df.index, \
                columns=cpop_score_df.columns)
        # Here we sum the Z-score across samples
        cpop_score_scale_series = cpop_score_scale_df.sum(axis=1)
        cpop_score_scale_sumdf = pd.DataFrame(cpop_score_scale_series,columns=["zscore"]).T
        cpop_score_scale_sumdf.insert(0,"Cluster",clustid)
        allclusters_cellpop_df.append(cpop_score_scale_sumdf)
        
    full_cluster_cellpop_df = pd.concat(allclusters_cellpop_df)
    full_cluster_deg_df = pd.concat(allclusters_degs_df)
    return full_cluster_cellpop_df, full_cluster_deg_df 

def main():
    """
    Used for testing this file.
    """
    deg_infile    = "../testing/input_type1_degs.tsv"
    cellpop_df    = scp.get_prop(species="mouse",mode="pandas_df")
    indf          = ir.read_file(deg_infile, mode="DEG") 
    fold_change_lim = 1.5   
    nof_clust       = 10
    full_clust_cpopdf, full_clust_degdf = \
    cluster_cellpop_score(cellpopdf=cellpop_df, degdf=indf,\
            fclim=fold_change_lim,\
        gene_count=False,logscale=False, k=nof_clust,\
        method="ward", dist="euclidean")
   
    print full_clust_cpopdf
    print full_clust_degdf
     

if __name__ == '__main__':
    main()
