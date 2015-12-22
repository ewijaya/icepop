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
from clock import clockit
import cell_proportion as cellprop


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
            hdf_output_organ = "../proportion_data/immgen_mouse_expr_organ.h5"
            df.to_hdf(hdf_output_organ, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
        else:
            ig_expr        = ig.ComputeAverage()
            ig_expr_dict   = ig_expr[0]
            df = pd.DataFrame.from_dict(ig_expr_dict).fillna(0).T
            df = df.convert_objects(convert_numeric=True)
            df.index.name = "Genes"
            df.reset_index(inplace=True)
            hdf_output = "../proportion_data/immgen_mouse_expr.h5"
            df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
    else:
        pass
    

def make_immgen_weight(out_type="HDF5",group_by=None):
    """
    Make the persistent files for ImmGen.
    Now we implement for HDF5 as Pandas table
    and Pickle as nested dictionary.
    """
    cell_prop_file = '../proportion_data/ImmgenCons_all_celltypes_MicroarrayExp.csv'
    ig = cellprop.ImmGen(cell_prop_file)
    igweight_dict, celltype_list  = ig.ComputePercentageWeight()

    if out_type =="HDF5":
        if group_by == 'organ':
            df            = ig.ComputeAverageByOrgans(weighted=True)
            hdf_output = "../proportion_data/immgen_mouse_organ.h5"
            df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
        else:
            df = pd.DataFrame.from_dict(igweight_dict).fillna(0).T
            df.index.name = "Genes"
            df.reset_index(inplace=True)
            hdf_output = "../proportion_data/immgen_mouse.h5"
            df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
    elif out_type == "pickle":
        pickle_output = "../proportion_data/immgen_mouse.pkl"
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
    cell_prop_file = '../proportion_data/IRIS.csv'
    iris = cellprop.Iris(cell_prop_file)
    irisweight_df, celltype_list  = iris.ComputePercentageWeight()
    if out_type=="HDF5":
        hdf_output = "../proportion_data/iris_human.h5"
        irisweight_df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)

        
def make_iris_expr(out_type="HDF5"):
    """
    Make the persistent files for IRIS.
    Now we implement for HDF5 as Pandas table
    and Pickle as nested dictionary.
    """
    cell_prop_file = '../proportion_data/IRIS.csv'
    iris = cellprop.Iris(cell_prop_file)
    iris_expr_df  = iris.ComputeAverage()
    iris_expr_df  = iris_expr_df.reset_index()
    iris_expr_df.drop('Probes',axis=1,inplace=True)
    if out_type=="HDF5":
        hdf_output = "../proportion_data/iris_human_expr.h5"
        iris_expr_df.to_hdf(hdf_output, 'fixed', mode='w', complib='zlib',complevel=2, append=False)
    

def main():
    """Make the actual HDF files for both ImmGen and IRIS."""
    # make_iris_expr(out_type="HDF5")
    # make_immgen_expr(out_type="HDF5",group_by="organ")
    # make_immgen_weight(out_type="HDF5")
    make_immgen_weight(out_type="HDF5",group_by="organ")
    # make_iris(out_type="HDF5")
    # make_immgen(out_type="pickle")
    
    


if __name__ == '__main__':
    main()
