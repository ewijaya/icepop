#!/usr/bin/env python
""" 
This code is only to be used once. 
It creates Pickle or `HDF5 <https://en.wikipedia.org/wiki/Hierarchical_Data_Format>`_ format files.
So that faster loading of data structure can be made.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import csv
import json
import os
import argparse
import re
import pandas as pd
import pickle
import h5py
from collections import defaultdict
from clock import clockit
import cell_proportion as cellprop
import input_reader as ir
import specificity as sp


@clockit
def make_pickle(indict):
    """
    Function to make pickle file. 
    This is considered not safe and slow.
    """
    pickle_output = "proportion_data/immgen_mouse.pkl"
    output = open(pickle_output,"wb")
    pickle.dump(indict,output)
    output.close()
    return    

def make_immgen_expr_cpop_phendict(out_type="HDF5"):
    """
    From original file, we create
    a HDF5 that list all the cell population and phenotype as
    data frame
    """
    cell_prop_file = './proportion_data/ImmgenCons_all_celltypes_MicroarrayExp.csv'
    df = pd.read_csv(cell_prop_file,header=None,index_col=[1,2], nrows=3, skiprows=[4])
    df = df.ix[:, 1:].reset_index(drop=True).T
    df.columns = ['celltype','organs', 'phenotype']
    df.reset_index(drop=True,inplace=True)
    hdf_output = "./proportion_data/immgen_mouse_celltype_phenotype.h5"
    df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
    return
    
def make_immgen_expr_phenotype(out_type="HDF5"):
    """
    Obtain gene expression without summarizing, 
    i.e. include all phenotypes.
    """
    cell_prop_file = './proportion_data/ImmgenCons_all_celltypes_MicroarrayExp.csv'
    df = pd.io.parsers.read_table(cell_prop_file,sep=",",skiprows=[0,1])
    df.drop(['ProbeSetID','Description'],axis=1,inplace=True)
    colnames = df.columns.values.tolist()
    phenlist = colnames[1:]
    df.columns = ["Genes"] + phenlist
    hdf_output = "./proportion_data/immgen_mouse_expr_by_phenotype.h5"
    df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
    return 

    

def make_immgen_expr(out_type="HDF5", group_by="organs"):
    """
    Make the persistent files for ImmGen.
    We implement only for HDF5 as Pandas table
    But here we just use expression without weighting them.
    """

    cell_prop_file = '../proportion_data/ImmgenCons_all_celltypes_MicroarrayExp.csv'
    ig             = cellprop.ImmGen(cell_prop_file)

    if out_type =="HDF5":
        # Only IMMGEN contain organ
        if group_by == "organ":
            df            = ig.ComputeAverageByOrgans()

            hdf_output_organ = "./proportion_data/immgen_mouse_expr_organ.h5"
            df.to_hdf(hdf_output_organ, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
        else:
            ig_expr        = ig.ComputeAverage()
            ig_expr_dict   = ig_expr[0]
            df = pd.DataFrame.from_dict(ig_expr_dict).fillna(0).T

            df = df.convert_objects(convert_numeric=True)
            df.index.name = "Genes"
            df.reset_index(inplace=True)
            hdf_output = "./proportion_data/immgen_mouse_expr.h5"
            df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
    else:
        pass
    

def make_immgen_weight(out_type="HDF5",group_by=None):
    """
    Make the persistent files for ImmGen.
    Now we implement for HDF5 as Pandas table
    and Pickle as nested dictionary.
    """
    cell_prop_file = './proportion_data/ImmgenCons_all_celltypes_MicroarrayExp.csv'
    ig = cellprop.ImmGen(cell_prop_file)
    igweight_dict, celltype_list  = ig.ComputePercentageWeight()

    if out_type =="HDF5":
        if group_by == 'organ':
            df            = ig.ComputeAverageByOrgans(weighted=True)
            df = (df/df.sum()) * 100
            df = df.fillna(0)
            hdf_output = "./proportion_data/immgen_mouse_organ.h5"
            df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
        else:
            df = pd.DataFrame.from_dict(igweight_dict).fillna(0).T

            # Normalize column so that it sums to 1.
            # This makes the interpretation easier. 
            # So that when the fold change is constant for all genes 
            # the histogram will also be flat.
            df = (df/df.sum()) * 100

            df.index.name = "Genes"
            df.reset_index(inplace=True)
            hdf_output = "./proportion_data/immgen_mouse.h5"
            df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
    elif out_type == "pickle":
        pickle_output = "./proportion_data/immgen_mouse.pkl"
        output = open(pickle_output,"wb")
        pickle.dump(igweight_dict,output)
        output.close()
        
    return
    
def make_iris(out_type="HDF5"):
    """
    Make the persistent files for IRIS.
    Now we implement for HDF5 as Pandas table
    and Pickle as nested dictionary.
    """
    cell_prop_file = './proportion_data/IRIS.csv'
    iris = cellprop.Iris(cell_prop_file)
    irisweight_df, celltype_list  = iris.ComputePercentageWeight()

    irisweight_df.set_index("Genes",inplace=True)
    irisweight_df = (irisweight_df/irisweight_df.sum()) * 100
    irisweight_df.index.name="Genes"
    irisweight_df.reset_index(inplace=True)
    irisweight_df["Genes"] = irisweight_df['Genes'].str.rstrip()
    irisweight_df.drop_duplicates('Genes', keep='first',inplace=True)
    # print irisweight_df.head(n=30L)

    if out_type=="HDF5":
        hdf_output = "./proportion_data/iris_human.h5"
        irisweight_df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)

        
def make_iris_expr(out_type="HDF5"):
    """
    Make the persistent files for IRIS.
    Now we implement for HDF5 as Pandas table
    and Pickle as nested dictionary.
    """
    cell_prop_file = './proportion_data/IRIS.csv'
    iris = cellprop.Iris(cell_prop_file)
    iris_expr_df  = iris.ComputeAverage()
    iris_expr_df  = iris_expr_df.reset_index()
    iris_expr_df.drop('Probes',axis=1,inplace=True)
    if out_type=="HDF5":
        hdf_output = "./proportion_data/iris_human_expr.h5"
        iris_expr_df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
    
def make_immgen_specificity(species='mouse',method='sparseness',to_exclude=['gdTCells']):
    """
    Make specificity file in HDF5.
    The data structure is mouse.
    """
    cellpop_df = ir.read_hdf_expr(species=species)
    mg_df = cellpop_df
    if to_exclude == None:
        tcells_subset   = ["abTcells","gdTCells"]
        mg_df['Tcells'] = mg_df[tcells_subset].mean(axis=1)
        mg_df.drop(tcells_subset,axis=1,inplace=True)
        col_list = list(mg_df)
        # swap last two columns
        col_list[-1], col_list[-2] = col_list[-2], col_list[-1]
        mg_df.columns = col_list
        # pass
    else:
        mg_df.drop(to_exclude[0],axis=1,inplace=True)
    

    # Using circular package calling
    # possibly harmful
    mg_df = sp.assign_specificity_score(mg_df,method=method)    
    # Sort genes by specificity desending
    mg_df.sort_values(method,axis=0,ascending=False, inplace=True)
    colnames = mg_df.columns.values.tolist()
    celltypes = colnames[1:-1]

    outerdict = defaultdict(list)
    for i, row in mg_df.iterrows():
        rdf =  pd.DataFrame(row)
        # capture gene row
        genedf = rdf.iloc[[0]].to_dict(orient='records')
        genename = genedf[0].values()[0]
        # remove gene row
        rdf = rdf.ix[1:]

        # capture specificity score row
        score_df = rdf.iloc[[-1]].to_dict(orient='records')
        score    = score_df[0].values()[0]
        # remove specificity score row
        rdf = rdf.ix[:-1]
        rdf_dict = rdf.to_dict(orient='dict')
        rdf_dict = rdf_dict.values()[0]

        # pick cell type with highest expression
        sorted_rdf_dict = sorted(rdf_dict.items(), key=lambda x: (-x[1], x[0])) 
        highest_ctexpr = sorted_rdf_dict[0]
        highct = highest_ctexpr[0]
        outerdict[highct].append(genename + " " + "({:1.3f})".format(score)  )

    to_include = {'gdTCells':'abTcells', 'abTcells':'gdTCells'}
    
    pickle_output = "./proportion_data/immgen_mouse_celltype_specificity_" + to_include[to_exclude[0]] + ".pkl"
    output = open(pickle_output,"wb")
    pickle.dump(outerdict,output)
    output.close()
    

        
    
    

def main():
    """Make the actual HDF files for both ImmGen and IRIS."""
    # make_iris_expr(out_type="HDF5")
    # make_immgen_expr(out_type="HDF5",group_by="organ")
    # make_immgen_expr_organ(out_type="HDF5",group_by="organ")
    # make_immgen_expr_cpop_phendict(out_type="HDF5")
    # make_immgen_expr_phenotype(out_type="HDF5")
    # make_immgen_specificity(species='mouse')

    # make_immgen_weight(out_type="HDF5")
    # make_immgen_weight(out_type="HDF5",group_by="organ")

    make_iris(out_type="HDF5")
    # make_immgen(out_type="pickle")
    
    


if __name__ == '__main__':
    main()
