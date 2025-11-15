#!/usr/bin/env python
""" 
Various function to calculate specificity
score of every gene in proportion data (ImmGen/IRIS)

"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import json
import sys
import os
import numpy as np
import operator
from collections import defaultdict
from numpy import linalg as LA
import pandas as pd
import input_reader as ir


def bind_sparseness_score(indf=None):
    """
    Given a Pandas data frame with genes
    as row and column of cell types or samples
    append sparseness score at the end
    
    :indf: Pandas data frame.
    :returns: Pandas data frame with appended sparseness score.
    """

    nindf = indf.copy()
    genelist   = nindf.iloc[:,0]
    val_df  = nindf.iloc[:,1:]
    
    sparse_list = sparseness(val_df.values)
    nindf["sparseness_score"] = sparse_list

    return nindf




def jsd_vectorized(xs=None):
    """
    Using Jensen-Shannon distance calculate 
    the specificity score of every gene with respect to the  cell type.

    :param xs: Numpy matrix, which contain list of expression from all cell types.
    :returns: Numpy array (1D),  value with sparsity scores for every genes.
    """

    xs = normalize(xs)
    nof_celltype = xs.shape[1]
    I = np.eye(nof_celltype)
    
    # Perform "(x+y)" in a vectorized manner
    xsI = xs+I[:,None]

    # Calculate d1 and d2 using xsI again in vectorized manner
    d1 = xs*np.log2(2*xs/xsI)
    d2 = I[:,None,:]*np.log2((2*I[:,None,:])/xsI)

    # Use np.nansum to ignore NaNs & sum along rows to
    # get all distances
    dists = np.nansum(d1,2) + np.nansum(d2,2)

    # Pack the argmax IDs and the
    # corresponding scores as final output   
    ID = dists.argmax(0)

    # First column is JS distance,
    # Second column is cell type index
    out = np.vstack((0.5*dists[ID,np.arange(dists.shape[1])],ID)).T
    return out[:,0]



def normalize(dfv=None): 
    """
    Normalize every row.

    :param dfv: Numpy matrix.
    :returns: Numpy matrix, with normalize score.
    """
    dfvlog2 = np.log2(dfv + 1)
    final = dfvlog2 / dfvlog2.sum(axis=1)[:,None]
    return final
    

def sparseness(xs=None):
    """
    A function to calculate sparseness score
    of a given gene.
    The value range from 0 to 1. Score 1 means
    that the gene is perfectly expressed in one
    cell type. Calcuated the following way (Hoyer 2004):

        .. math::
            
            sparseness(x) = \\frac{\sqrt{n}-\\frac{\\sum{||x_{i}||}}{\\sqrt{\\sum{x_{i}^2}}}}{\\sqrt{n}-1}

    Here :math:`n` refers to the number of cell types and :math:`x_i`
    expression of a gene in cell type :math:`i`.

    :param xs: Numpy matrix, which contain list of expression from all cell types.
    :returns: Numpy array (1D),  value with sparsity scores for every genes.
    """

   
    nr = np.sqrt(xs.shape[1])
    a  = np.sum(np.abs(xs),axis=1)
    b  = np.sqrt(np.sum(np.square(xs),axis=1))
    sparseness = (nr - a/b) / (nr -1)
    return sparseness

def assign_specificity_score(df, method='sparseness'):
    """
    Assign specificity score to the ImmGen/IRIS data. 
    :param df: Pandas data frame, generally proportion file.
    :param method: str("sparseness","js")

    :returns: a Panda data frame with where every genes will have its specificity score.
    """
    genelist = df.iloc[:,0]
    express_df   = df.iloc[:,1:]

    # This is faster implementation than using 'apply'
    # The function is implemented on 2-D numpy array.
    marker_genes_df  = None

    if method=="sparseness":
        # more suited for detecting marker genes
        marker_genes_df  = sparseness(express_df.values)
    else:
        marker_genes_df  = jsd_vectorized(express_df.values)

    mg_df = df.copy()
    mg_df[method] = marker_genes_df
    return mg_df
    

def find_topk_marker_genes(df, method='sparseness',to_exclude=None,lim=0.8,top_k=1):
    """
    Find marker genes from expression cell population by
    some specificity score. Then select to 10 genes, 
    one gene which is most specific in one cell type. 
    

    :param df: Pandas data frame, generally proportion file.
    :param to_exclude: list of cell type to exclude. If None,
                       average between abTcell and gdTcell will be performed.
    :param method: str("sparseness","js")

    :returns: a Panda data frame with selected marker genes.
    """
    outerdict = ir.read_specificity_pickle()
     
     
    mg_df = df.copy()
    # get top_k genes for every cell types
    all_sel_genes = []
    for ct,topgenes in outerdict.items():
        sel_genes = topgenes[0:top_k] 
        sel_genes = [ x.split()[0] for x in sel_genes] 
        # print ct, ":",  ",".join(sel_genes) 
        all_sel_genes += sel_genes
        
    mg_df = mg_df[mg_df['Genes'].isin(all_sel_genes)]
    tmp_marker_mat = mg_df.iloc[:,1:-1].values
    condn          = LA.cond(tmp_marker_mat)
    mg_df.drop(to_exclude[0],axis=1,inplace=True)
    # print mg_df
    # print "CN:", condn
    return mg_df
    

def find_marker_genes(df, method="sparseness", lim=0.8):
    """
    Find marker genes from expression cell population by
    some specificity score.

    :param df: Pandas data frame, generally proportion file.
    :param method: str("sparseness","js")
    :param lim:  float, threshold to select 
    
    """
    mg_df = assign_specificity_score(df,method=method)    
    mg_df = mg_df.loc[mg_df[method] >= lim]
    mg_df.drop(method,axis=1,inplace=True)
   
    return mg_df

def condn_thres_mat(df,method='sparseness',verbose=False):
    """
    Function to enumerate condition number from series of thresholds.

    :param df: Pandas data frame, generally proportion file.
    :returns: numpy (3 x 2) matrix, which stores list of threshold, condition \
        number and number of markers.
    """
    conds = np.empty((10,3))
    for i, splim in enumerate(np.arange(0,1,0.1)):
        marker_df_test = find_marker_genes(df,method=method, lim=splim)
        tmp_marker_mat = marker_df_test.iloc[:,1:].values
        condn          = LA.cond(tmp_marker_mat)
        nof_markers    = tmp_marker_mat.shape[0]
        conds[i,0] = splim
        conds[i,1] = condn
        conds[i,2] = nof_markers
        if verbose:
            condout = str(splim) + " " + str(tmp_marker_mat.shape[0]) + " " +str(condn)
            sys.stderr.write(condout + "\n") 

    return conds  
    
def find_best_marker_genes(df, method="sparseness",verbose=False):
    """
    We use `condition number` to choose the best
    specificity threshold.

    :param df: Pandas data frame, generally proportion file.
    :param method: str("sparseness","js")
    :param verbose: boolean
    :returns: a Panda data frame with selected marker genes.
    
    """
    # print "FBMG"
    # print df.head()
    conds = condn_thres_mat(df,method=method,verbose=verbose)
    # print repr(conds)
    # Index where the condition number is minimum 
    minid =  conds.argmin(0)[1]
    min_lim, min_connum, nof_marker_genes = conds[minid]

    if verbose:
        sys.stderr.write("Best threshold    : " +  str(min_lim) +  "\n")
        sys.stderr.write("Condition number  : " +  str(min_connum) +   "\n")
        sys.stderr.write("Marker genes count: " +  str(nof_marker_genes) +   "\n")
    # print min_lim, minconid
    return find_marker_genes(df,method=method, lim=min_lim)
    

if __name__ == '__main__':
    main()
