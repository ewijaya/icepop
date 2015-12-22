#!/usr/bin/env python
""" 
Functions for calculating cell population.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"

import numpy as np
import sys
import pandas as pd
import species_cell_proportion as scp
import input_reader as ir

def deg_cellpopscore_df(cellpopdf=None, degdf=None, fclim=None, gene_count=False, logscale=False):
    """
    Calculate the cell population score per sample.
    All done in Pandas data frame. 
    
    :param cellpopdf: Pandas data frame, cell population.
    :param degdf: Pandas data frame, differentially expressed genes (DEGs).
    :param fclim: float, foldchange lower limit thresold.
    :param logscale: boolean, transform fold change to logscale.
                  This does not affect the result substansially.
    :param gene_count: boolean, normalization method. 
                  If *True*, divide by sum of product. The total weight
                  will have to result to 1.00.

    :returns: Score per cell population as data frame.


    Usage:

    >>> from icepop import cellpop_score as cp
    >>> cpop_score_df = cp.deg_cellpopscore_df(cellpop_df, indf, fclim=0, gene_count=False, logscale=False)

    """
    
    # print cellpopdf.head()
    sample_names = degdf.columns.values.tolist()
    sample_names = sample_names[2:]

    degdf.columns = ["probe","Genes"] + sample_names 

    

    all_dfs = []
    for i, sample in enumerate(sample_names):
        unwanted_samples = list(set(sample_names) - set([sample]))
        sample_df = degdf[degdf.columns.difference(unwanted_samples)]
        sample_df = sample_df[["probe","Genes", sample]]
        # Standardized the colum name so that
        # it can be fitted with cellpopdf
        # sample_df.columns = ["probe","Genes", sample]
        sample_df = sample_df[["Genes", sample]]
        
        sample_df = sample_df[ sample_df[sample] >=  fclim ]
        sample_df["Genes"] = sample_df["Genes"].str.strip()
        #print sample_df.shape

        # prune out sample with FC < 0 
        if logscale:
            sample_df = sample_df[ sample_df[sample] > 0 ]
            sample_df[sample] = np.log(sample_df[sample])

        # print cellpopdf.head()
        cellpop_score_series = get_population_values(cellpop_df = cellpopdf,\
                                             sample_df  = sample_df,\
                                             input_gene_count = gene_count)

        cellpop_prop_df = pd.DataFrame(cellpop_score_series, columns=[sample])
        all_dfs.append(cellpop_prop_df)
    

    final_df = pd.concat(all_dfs,axis=1).fillna(0)
    #print final_df
    return final_df
    

def get_population_values(cellpop_df=None, sample_df=None, input_gene_count=None):
    """
    Here is where the actual calculation take place.
    
    :param cellpop_df: Pandas data frame, cell population.
    :param sample_df: Pandas data frame, sample fold change. 

    :param gene_count: boolean, normalization method. 
                If *True*,divide by sum of product. The total weight
                will have to result to 1.00.

    :returns: Pandas series, with score per cell population.
        If the output contains all zeros, it means the genes defined
        in ``sample_df`` is not included in ImmGen/IRIS database.
    """
    
    cellpop_colnames = cellpop_df.columns.values.tolist()
    celltypes = cellpop_colnames[1:]

    # print sample_df
    sample_name = sample_df.columns.values.tolist()
    sample_name = sample_name[1]

    # Case sensetive merge
    # merged_df = pd.merge(sample_df, cellpop_df, on="Genes")
    # print merged_df.head()


    # Case insensetive merge is more robust
    tmp_cellpop_df = cellpop_df.copy()
    tmp_sample_df  = sample_df.copy()
    tmp_cellpop_df['genes_lower'] = tmp_cellpop_df['Genes'].str.strip().str.lower()
    tmp_sample_df['genes_lower'] = tmp_sample_df['Genes'].str.strip().str.lower()
    merged_df = pd.merge(tmp_sample_df,tmp_cellpop_df,left_on="genes_lower",right_on="genes_lower",how="left")
    merged_df.drop(["Genes_y","genes_lower"],axis=1,inplace=True)
    merged_df.rename(columns={"Genes_x":"Genes"}, inplace=True)


    # Multiple each cell type with fold change
    multp_df = merged_df[celltypes].multiply(merged_df[sample_name],axis="index")

    # Get total values
    tot_prod_by_adjvtype = multp_df.sum(numeric_only=True).sum()
    tot_prod_by_adjv_by_celltype = multp_df.sum(numeric_only=True)
    adjvtype_gene_count = sample_df[sample_df[sample_name] > 0].shape[0]


    cellpop_perc = tot_prod_by_adjv_by_celltype/tot_prod_by_adjvtype
    if input_gene_count:
        cellpop_perc = tot_prod_by_adjv_by_celltype/adjvtype_gene_count

    # print cellpop_perc.head()
        
    # print cellpop_perc.sum()
    return cellpop_perc

def main():
    """Used for testing this file. """ 
    cellpop_df = scp.get_prop(species="mouse",mode="pandas_df") 
    deg_inputfile="../testing/input_type1_degs.tsv"
    indf = ir.read_file(deg_inputfile, mode="DEG")
    print indf.head()
    cpop_score_df = deg_cellpopscore_df(cellpopdf=cellpop_df, degdf=indf, fclim=0, gene_count=False, \
            logscale=False)
    print cpop_score_df.head()


if __name__ == '__main__':
    main()
