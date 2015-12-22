#!/usr/bin/env python
""" 
This modules has interface to let user
choose which species they want. And return
the cell population weight accordingly.
Additional option include choosing between
Pandas dataframe or standard nested dictionary.
It behaves as the extension of ``cell_proportion`` module.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import csv
import json
import os
import pandas as pd
import pickle
from clock import clockit
from input_reader import read_hdf




def choose_organ(cp_organ_df, organ_name):
    """
    Description of choose_organ
    """
    cp_organ_df = cp_organ_df.loc[:, organ_name]
    cp_organ_df.reset_index(level=0,inplace=True)
    cp_organ_df =  cp_organ_df.rename(columns = {'cell':"Genes"})
    return cp_organ_df
    

def get_prop(species="mouse",mode=""):
    """
    Get cell type proportion. 
    
    :param species: string('mouse','human')
        mouse (default) obtained from ImmGen database and human from IRIS.
    :param mode: string('pandas_df','dict')
    :returns: function to read persistence file.

    Usage:

    >>> from icepop import species_cell_proportion as scp
    >>> cellpop_df = scp.get_prop(species="mouse",mode="pandas_df")

    """
    if (species=="mouse" or species=="human") and mode=="pandas_df":
        return read_hdf(species=species)
    elif species=="mouse" and mode=="dict":
        return read_pickle()
    else:
        pass
        # IRIS to be added
        


def main():
    """Used to execute this file."""
    # outdf   = get_prop(species="mouse",mode="pandas_df") 
    outdf   = get_prop(species="human",mode="pandas_df") 
    print outdf.head()

    # outdict = get_prop(species="mouse",mode="dict") 
    # print json.dumps(outdict, indent=4)
    
    
if __name__ == '__main__':
    main()
