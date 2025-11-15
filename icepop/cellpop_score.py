#!/usr/bin/env python
""" 
Functions for calculating cell population.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import json
import math
import numpy as np
import warnings
import sys
from collections import defaultdict
import pandas as pd
import species_cell_proportion as scp
import input_reader as ir
from . import constants as const


def _process_sample(cellpopdf, sample_df, sample, fclim, logscale, gene_count):
    """
    Helper function to process a single sample.
    Extracted to reduce code duplication between deg_cellpopscore_generank and deg_cellpopscore_df.

    :param cellpopdf: Cell population DataFrame
    :param sample_df: Sample DataFrame with fold change data
    :param sample: Sample name
    :param fclim: Fold change limit
    :param logscale: Whether to use log scale
    :param gene_count: Whether to use gene count normalization
    :returns: Filtered sample DataFrame and number of genes
    """
    sample_df = sample_df[[sample_df.columns[0], sample]]
    sample_df = sample_df[sample_df[sample] >= fclim]
    sample_df[sample_df.columns[0]] = sample_df[sample_df.columns[0]].str.strip()
    nof_genes = sample_df.shape[0]

    # Prune out sample with FC < 0 if logscale
    if logscale:
        sample_df = sample_df[sample_df[sample] > 0]
        sample_df[sample] = np.log(sample_df[sample])

    return sample_df, nof_genes


def deg_cellpopscore_generank(cellpopdf=None, degdf=None, fclim=None, gene_count=False, logscale=False,
        cpop_thres='median'):
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
    :param cpop_thres: string('uniform','median','q1'), threshod to choose which cell population
                    is induced. This threshold is for cell population score.
                    It affects the cell type list to be reported. If uniform,
                    we apply 1/(number of cell types).

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
    all_gene_x_ctypeweight_df = []
    sample_ct_generank_dict = defaultdict(lambda: defaultdict(list))
    threshold_dict = {}
    
    for i, sample in enumerate(sample_names):
        # Optimized: Use list comprehension instead of set operations
        unwanted_samples = [s for s in sample_names if s != sample]
        sample_df = degdf[degdf.columns.difference(unwanted_samples)]
        sample_df = sample_df[["probe","Genes", sample]]
        # Standardized the column name so that it can be fitted with cellpopdf
        sample_df = sample_df[["Genes", sample]]

        # Use helper function to process sample
        sample_df, nof_genes = _process_sample(cellpopdf, sample_df, sample, fclim, logscale, gene_count)

        cellpop_score_series = get_population_values(cellpop_df = cellpopdf,\
                                             sample_df  = sample_df,\
                                             input_gene_count = gene_count)
        cellpop_score_series_round = cellpop_score_series.round(decimals=const.CELLPOP_SCORE_DECIMALS)

        # Calculate median here.
        # The threshold is Median + Stdev(Median[idx-1],Median[idx],Median[idx+1])
        cpop_score_list     = cellpop_score_series.tolist()
        cpop_score_list.sort()
        cpop_score_thres = find_threshold(cpop_score_list,method=cpop_thres)
            
        threshold_dict[sample] = cpop_score_thres

        filt_cellpop_score_series = cellpop_score_series_round[cellpop_score_series_round > cpop_score_thres]
        wanted_cellpop = filt_cellpop_score_series.index.values.tolist()
        ct_generank_dict, genes_x_ctypeweight_df = get_genebygene_fc_score_matrix_product(cellpop_df = cellpopdf,\
                                             sample_df  = sample_df,\
                                             input_gene_count =
                                             gene_count,filt_cellpop=wanted_cellpop)

        genes_x_ctypeweight_df["Sample"] = sample
        sample_ct_generank_dict[sample] = ct_generank_dict
        cellpop_prop_df = pd.DataFrame(cellpop_score_series, columns=[sample])


        all_gene_x_ctypeweight_df.append(genes_x_ctypeweight_df)
        all_dfs.append(cellpop_prop_df)
    
    final_df = pd.concat(all_dfs,axis=1).fillna(0)
    final_gene_x_ctypeweight_df = pd.concat(all_gene_x_ctypeweight_df,axis=0)
    return sample_ct_generank_dict, final_df, threshold_dict, final_gene_x_ctypeweight_df


def sample_response_score(score_list,method=None):
    """

    Given a list of values, calculate the *sample response score* (SR).
    It calls :func:`find_threshold` function that calculates the *cell type response
    score* (CRT).


        .. math::
    
            SR = -log \\Bigg(\\frac{CRT}{(1/nof\_cell\_type)} \\Bigg) \\\\
                = -log (CRT \\cdot  nof\_cell\_type)

    """
    
    cpop_thres = find_threshold(score_list, method=method)
    nof_cell_type = len(score_list)

    cpop_thres_ct_prod = cpop_thres * nof_cell_type 

    sample_response = 0
    # Remove warnings
    if cpop_thres_ct_prod != 0:
        sample_response = -1 * np.log(cpop_thres_ct_prod)
    
    return sample_response

def find_threshold(score_list,method=None):
    """
    Given a list of values, calculate
    the median, median's lower and upper neighbors.
    This is the *cell type response threshold* (CRT).
    Then return the standard deviation of the three values above.
    """
    score_list              = [0 if math.isnan(x) else x for x in score_list]
    score_list              = sorted(score_list)
    score_median            = np.median(score_list)
    idx_low                 = (np.abs(np.array(score_list)-score_median)).argmin()
    idx_high                = idx_low 
    medians                 = score_list[idx_low], score_median, score_list[idx_high]
    score_median_std        = np.std(medians)
    score_thres             = score_median + score_median_std

    # Quartile based
    score_1quartile         = np.percentile(score_list, const.FIRST_QUARTILE_PERCENTILE, axis=0)
    idx_low_q               = (np.abs(np.array(score_list)-score_1quartile)).argmin()
    # idx_high_q              = idx_low_q  + 1
    # print score_list
    # print score_1quartile
    first_quartiles         = score_list[idx_low_q], \
                              score_1quartile,\
                              score_list[idx_low_q + 1]
    score_1quartile_std     = np.std(first_quartiles)
    score_thres_q           = score_1quartile + score_1quartile_std
        

    # print "ME:", score_thres
    # print "1Q:", score_thres_q


    if method == 'median':
        return np.round(score_thres, decimals=const.CELLPOP_SCORE_DECIMALS)
    elif method == 'q1':
        return np.round(score_thres_q, decimals=const.CELLPOP_SCORE_DECIMALS)
    elif method == 'q1std':
        return np.round(score_thres_q - score_1quartile_std, decimals=const.CELLPOP_SCORE_DECIMALS)
    elif method == 'uniform':
        return (1/(len(score_list)+0.00)) 




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


    all_dfs                 = []
    nof_genes_dict          = {}
    sample_response_dict    = {}
    celltype_response_dict  = {}
    for i, sample in enumerate(sample_names):
        # Optimized: Use list comprehension instead of set operations
        unwanted_samples = [s for s in sample_names if s != sample]
        sample_df = degdf[degdf.columns.difference(unwanted_samples)]
        sample_df = sample_df[["probe","Genes", sample]]
        # Standardized the column name so that it can be fitted with cellpopdf
        sample_df = sample_df[["Genes", sample]]

        # Use helper function to process sample
        sample_df, nof_genes = _process_sample(cellpopdf, sample_df, sample, fclim, logscale, gene_count)
        nof_genes_dict[sample] = nof_genes

        # print cellpopdf.head()
        cellpop_score_series = get_population_values(cellpop_df = cellpopdf,\
                                             sample_df  = sample_df,\
                                             input_gene_count = gene_count)

        cellpop_score_list = cellpop_score_series.tolist()
        
        celltype_response_thres = find_threshold(cellpop_score_list, method='q1std') 
        sample_response_val   = sample_response_score(cellpop_score_list, method='q1std') 
        sample_response_val   = "{:.3f}".format(sample_response_val)
        celltype_response_thres = "{:.3f}".format(celltype_response_thres)

        if sample_response_val == '-inf' or sample_response_val == 'inf':
            sample_response_val = const.INFINITY_REPLACEMENT

        sample_response_dict[sample] = sample_response_val
        celltype_response_dict[sample] = celltype_response_thres

            

        cellpop_prop_df = pd.DataFrame(cellpop_score_series, columns=[sample])
        all_dfs.append(cellpop_prop_df)
    

    final_df = pd.concat(all_dfs,axis=1).fillna(0)
    return nof_genes_dict, sample_response_dict, celltype_response_dict, final_df
    

def merge_scorematrix_samplefc(cellpop_df=None, sample_df=None, input_gene_count=None):
    """
    Merge Immgen/IRIS matrix with sample fold change
    data. It is possible that the gene list in sample
    does not exist at all in IRIS/ImmGen. In that case 
    we return the data frames with 0 values.

    :param cellpop_df: Pandas data frame, cell population.
    :param sample_df: Pandas data frame, sample fold change. 

    :param gene_count: boolean, normalization method. 
                If *True*,divide by sum of product. The total weight
                will have to result to 1.00.

    :returns: Pandas series, merged data.
    
    """


    # Case insensetive merge is more robust
    # We do that by lower casing both left and right keys
    tmp_cellpop_df = cellpop_df.copy()
    tmp_sample_df  = sample_df.copy()
    tmp_cellpop_df['genes_lower'] = tmp_cellpop_df['Genes'].str.strip().str.lower()
    tmp_sample_df['genes_lower'] = tmp_sample_df['Genes'].str.strip().str.lower()
    merged_df = pd.merge(tmp_sample_df,tmp_cellpop_df,left_on="genes_lower",right_on="genes_lower",how="left")
    merged_df.drop(["Genes_y","genes_lower"],axis=1,inplace=True)
    merged_df.rename(columns={"Genes_x":"Genes"}, inplace=True)
    merged_df = merged_df.fillna(0)

    return merged_df
    
def get_genebygene_fc_score_matrix_product(cellpop_df=None, sample_df=None, input_gene_count=None,filt_cellpop=None):
    """
    Calculate product of FC with ImmGen/IRIS weight
    gene by gene.
    
    :param cellpop_df: Pandas data frame, cell population.
    :param sample_df: Pandas data frame, sample fold change. 

    :param gene_count: boolean, normalization method. 
                If *True*,divide by sum of product. The total weight
                will have to result to 1.00.

    :param filt_cellpop: list of cell type that passed threshold.
    :returns: 
    """
    
    merged_df = merge_scorematrix_samplefc(cellpop_df, sample_df,input_gene_count)

    cellpop_colnames = cellpop_df.columns.values.tolist()
    celltypes = cellpop_colnames[1:]

    # print sample_df
    sample_name = sample_df.columns.values.tolist()
    sample_name = sample_name[1]

    # Multiple each cell type with fold change
    multp_df = merged_df[celltypes].multiply(merged_df[sample_name],axis="index")
    multp_df.dropna(axis=0, how='all',inplace=True)
    genes_df = merged_df[["Genes"]]
    genes_multp_df = pd.merge(genes_df,multp_df,left_index=True,right_index=True)


    # print genes_multp_df.head()
    # tmp_path = "/Users/ewijaya-macmini-ifrec/Desktop/foldchange_x_celltype_weight/"
    # tmp_path = "/Users/ewijaya-macmini-ifrec/Various_Projects/CellPopulation_Framework_in_making/testing/degs_based_analysis/foldchange_x_celltype_weight/"
    # tmp_outfile = tmp_path + "/" + sample_name + "_fc_x_ctypeweight.xlsx"
    # tmp_genes_multp_df = genes_multp_df.copy()
    # tmp_genes_multp_df.to_excel(tmp_outfile,index=False,float_format='%.3f')
    # print sample_name
    # print tmp_outfile
    # print tmp_genes_multp_df


    outerdict    = defaultdict(list)
    for fct in filt_cellpop:
        # print fct
        genes_ct_df = genes_multp_df[["Genes", fct]]
        genes_ct_df  = genes_ct_df[np.isfinite(genes_ct_df[fct])]
        genes_ct_df.sort_values(by=fct,ascending=False,inplace=True)
        # print genes_ct_df
        ct_gene_list = genes_ct_df["Genes"].tolist()
        outerdict[fct] = ct_gene_list
        
    # print json.dumps(outerdict, indent=4)
    return outerdict, genes_multp_df 

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
    
    merged_df = merge_scorematrix_samplefc(cellpop_df, sample_df,input_gene_count)

    cellpop_colnames = cellpop_df.columns.values.tolist()
    celltypes = cellpop_colnames[1:]

    # print sample_df
    sample_name = sample_df.columns.values.tolist()
    sample_name = sample_name[1]

    # Multiple each cell type with fold change
    # Replace NaN and Inf with 0 (must check correctness)
    # Optimized: Combine replace operations
    multp_df = merged_df[celltypes].multiply(merged_df[sample_name],axis="index").replace([np.nan, np.inf], 0)
    # print multp_df
    

    # Get total values
    tot_prod_by_adjvtype = multp_df.sum(numeric_only=True).sum()
    tot_prod_by_adjv_by_celltype = multp_df.sum(numeric_only=True)
    adjvtype_gene_count = sample_df[sample_df[sample_name] > 0].shape[0]


    cellpop_perc = tot_prod_by_adjv_by_celltype/tot_prod_by_adjvtype
    # print "===="
    # print tot_prod_by_adjvtype
    # print tot_prod_by_adjv_by_celltype
    # print cellpop_perc
    # print "===="
    if input_gene_count:
        cellpop_perc = tot_prod_by_adjv_by_celltype/adjvtype_gene_count

    # print cellpop_perc.head()
        
    # print cellpop_perc.sum()
    return cellpop_perc

# Test code moved to tests/ directory
