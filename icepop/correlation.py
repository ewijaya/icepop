#!/usr/bin/env python
""" 
This module is used for cell purity type of 
analysis. Correlation of input user expression
profile with all ImmGen/IRIS cell type at phenotype
level is needed, before performing linear programming
later.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
from collections import defaultdict
import json
import pandas as pd

def average_expression(phendict=None,cpop_phen_exp_df=None, ctype_phen_df=None,by=None):
    """
    After obtaining good phenotype from correlation calculation,
    we would like to create matrix which average expression for 
    every cell type.

    :param phhendict: Dictionary of list that contain sample and its selected phenotypes.
    :param cpop_phen_exp_df: Pandas dataframe of cell type expression at phenotype level.
    :param ctype_phen_df: Pandas dataframe of two colums (cell type and phenotype).

    :returns: A generator.
    """
    
    # Select only 
    nctype_phen_df = ctype_phen_df[["celltype","phenotype"]]
    # print cpop_phen_exp_df.head()
    # print json.dumps(phendict, indent=4)

    
    for sample,phens in phendict.iteritems():
        new_df = cpop_phen_exp_df[["Genes"] + phens]
        new_df.set_index("Genes",inplace=True)
        col_dict = nctype_phen_df.set_index('phenotype').squeeze().to_dict()
        avg_df = new_df.groupby(col_dict, axis=1).mean()
        avg_df.reset_index(inplace=True)
        yield sample, avg_df

        

def correlation(userinput_df=None, cpop_phen_exp_df=None, ctype_phen_df=None, method='pearson'):
    """
    Perform correlation.
    
    :param userinput_df: Pandas dataframe of gene expression from user input.
    :param cpop_phen_exp_df: Pandas dataframe of cell type expression at phenotype level.
    :param ctype_phen_df: Pandas dataframe of two colums (cell type and phenotype).
    :param method: {'pearson','kendall','spearman'}
    """
    
    # Clean input data
    input_colnames = userinput_df.columns.values.tolist()
    input_samples = input_colnames[2:]
    userinput_df = userinput_df.drop(userinput_df.columns[0],axis=1,inplace=False)
    userinput_df.columns = ["Genes"] + input_samples
    
    outerdict = defaultdict(list)
    tmp_alldf = []
    for sample in input_samples:
        sample_df = userinput_df[["Genes",sample]]
        merged_df = pd.merge(cpop_phen_exp_df, sample_df, on=["Genes"],
                left_index=True, right_index=True)

        
        # Normalize across rows for every columns
        # merged_df.set_index(keys="Genes",inplace=True)
        # merged_df = (merged_df/merged_df.sum()) 
        # merged_df.reset_index(inplace=True)
        # print merged_df.head()


        # Calculate correlation of input with every cell type
        corr_series = merged_df.corr(method=method).iloc[:-1,-1] 
        corr_df = pd.DataFrame(corr_series)
        corr_df.sort_values(by=sample, axis=0,ascending=False, inplace=True)

        # To show Aoshi
        tmp_corr_df = corr_df.copy()
        tmp_corr_df.index.name = "phenotype"
        tmp_corr_df.reset_index(inplace=True)
        tmp_merge_df = pd.merge(tmp_corr_df, ctype_phen_df, on=["phenotype"])
        tmp_merge_df = tmp_merge_df[['phenotype','celltype',sample]]
        tmp_alldf.append(tmp_merge_df)

        # Check if all correlation are negative
        corr_vals = corr_df[sample].tolist()
        all_neg = all(x <0 for x in corr_vals)

        if all_neg:
            # Pick top 10
            corr_df = corr_df.head(n=10L)
        else:
            # Choose positive correlation
            # corr_df = corr_df.head(n=100L) # version 3
            corr_df = corr_df[corr_df[sample] > 0] # version 2

        corr_df.index.name = "phenotype"
        corr_df.reset_index(inplace=True)

        # Merge with cell type dict
        ncorr_df = pd.merge(corr_df, ctype_phen_df, on=["phenotype"])
        phenlist = ncorr_df["phenotype"].tolist()
        outerdict[sample] = phenlist

    merged_df = reduce(lambda ldf, rdf: pd.merge(ldf,rdf, on=["phenotype","celltype"]), tmp_alldf) 
    # print merged_df.head()
    return outerdict, merged_df
    

if __name__ == '__main__':
    main()
