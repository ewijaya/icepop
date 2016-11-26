#!/usr/bin/env python
""" 
Functions to read input with various format:

    - DEG (fold change)
    - RAW expression data from GEO format.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import pandas as pd
import pickle
import sys
import os
import os.path

def read_hdf_from_file(infile=None):
    """
    HDF reading from input file.

    :param infile: A HDF file, usually from ImmGen/IRIS proportion.
                Made by serialize_cell_proportion.py.

    """
    return pd.read_hdf(infile,'fixed')
   
def read_specificity_pickle():
    """
    Description of read from pickle file.
    The output is dictionary of list which 
    contain list of genes sorted by the most specific one.
    The file is created by serialized_cell_proportion.py
    """
    tmp = os.path.abspath(__file__)
    package_path = os.path.dirname(tmp)
    pkl_filevar = os.path.join(package_path, "proportion_data/immgen_mouse_celltype_specificity_abTcells.pkl")

    outerdict = pickle.load(open(pkl_filevar,"rb"))
    return outerdict
    
    
    

def read_hdf(species=None,organs=False):
    """ 
    HDF reading is faster than pickle.
    Here given the species, return the weighted proportion
    as dataframe.

    :param species: species, used to indicate the h5 file.
    :returns: Pandas data frame.

    """
    # With __file__ we set the relative path where we can
    # access the proportion data, relative to this file.


    tmp = os.path.abspath(__file__)
    package_path = os.path.dirname(tmp)

    if species=='mouse' and organs:
        infile = os.path.join(package_path, "proportion_data/immgen_mouse_organ.h5")
        out_df = pd.read_hdf(infile,'fixed')
        return out_df
    elif species=='human' and organs:
        sys.exit("H. Sapiens data has no organs option")
    else:
        infile = os.path.join(package_path, "proportion_data/immgen_mouse.h5")
        if species=="human":
            infile = os.path.join(package_path,"proportion_data/iris_human.h5")
        out_df = pd.read_hdf(infile,'fixed')
        return out_df

def read_hdf_expr(species=None,organs=False):
    """ 
    Here given the species, return the weighted expression
    for each celltypes as dataframe.

    :param species: species, used to indicate the h5 file.
    :returns: Pandas data frame.

    """
    # With __file__ we set the relative path where we can
    # access the proportion data, relative to this file.
    tmp = os.path.abspath(__file__)
    package_path = os.path.dirname(tmp)


    if species=='mouse' and organs:
        infile = os.path.join(package_path, "proportion_data/immgen_mouse_expr_organ.h5")
        out_df = read_hdf_from_file(infile=infile)
        return out_df
    else:
        infile = os.path.join(package_path, "proportion_data/immgen_mouse_expr.h5")
        if species=="human":
            infile = os.path.join(package_path,"proportion_data/iris_human_expr.h5")
        out_df = pd.read_hdf(infile,'fixed')
        return out_df


def read_hdf_pheno():
    """
    Used for reading only phenotype granularity
    type of ImmGen data. It returns two
    data frames: 1) scoring matrix; 2) phenotype celltype correspondence.
    """
    # With __file__ we set the relative path where we can
    # access the proportion data, relative to this file.
    tmp = os.path.abspath(__file__)
    package_path = os.path.dirname(tmp)
    exp_infile    = os.path.join(package_path, "proportion_data/immgen_mouse_expr_by_phenotype.h5")
    cpphen_infile = os.path.join(package_path, "proportion_data/immgen_mouse_celltype_phenotype.h5")

    out_exp_df = pd.read_hdf(exp_infile)
    out_cpphen_df = pd.read_hdf(cpphen_infile)
    return out_exp_df, out_cpphen_df
    

def load_hdf_cvfilter():
    """
    Reading HDF file that
    contain gene list and the corresponding CV filter,
    for the given organ. Calculated using `BARCODE 3.0 <http://barcode.luhs.org/>`_ .
    """
    tmp = os.path.abspath(__file__)
    package_path = os.path.dirname(tmp)
    cvfilt_infile    = os.path.join(package_path, "proportion_data/mouse_organ_genes.cvfilter.h5")
    out_cvfilt_df = pd.read_hdf(cvfilt_infile)
    return out_cvfilt_df
    

def read_file(infile, mode="DEG"):
    """Reading various input file.
    Supported format are tab delimited (TSV),
    comma delimited (CSV), Excel (.xlsx or .xls).

    :param mode: string('DEG','RAW').
        DEG (default) is in fold change format,
        RAW in GEO format of raw expression values. 
    
    :returns: Pandas data frame.

    Usage:

    >>> from icepop import input_reader as ir
    >>> deg_infile = "input_type1_degs.tsv"
    >>> indf       = ir.read_file(deg_infile, mode="DEG")

    """
    if mode=="DEG" or mode == "RAW":
        outdf = pd.DataFrame.empty
        if infile.endswith((".xlsx",".xls")):
            outdf = pd.read_excel(infile)
        elif infile.endswith(".csv"):
            outdf = pd.io.parsers.read_csv(infile)
        elif infile.endswith(".tsv"):
            outdf = pd.io.parsers.read_table(infile)
        return outdf
    else:
        # To be added
        pass
        

def main():
    """Description of main
    
    """
    outdf = load_hdf_cvfilter()
    print outdf
    

if __name__ == '__main__':
    main()    
