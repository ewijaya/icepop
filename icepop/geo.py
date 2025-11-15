#!/usr/bin/env python
""" 
Stores method to parse the `GEO <http://www.ncbi.nlm.nih.gov/geo/>`_ file
in various formats. These functions are inherited
from `GEOparse <https://geoparse.readthedocs.org>`_ package.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import json
import GEOparse
import pandas as pd

def load(geoid=None, destdir=None, filepath=None):
    """
    Function to download GEO file in SOFT format.

    :param geo: string, GEO id (e.g. "GSE69886")
    :param destdir: string, destination to store  the files.

    Usage:

    >>> from icepop import geo
    >>> gse = geo.load(geoid="GSE74306", destdir="./")
    
    or if the data is already downloaded

    >>> gse = geo.load(filepath="./GSE74306.soft.gz")

    """
    gse = None
    if filepath:
       gse = GEOparse.get_GEO(filepath=filepath) 
    else:
       gse = GEOparse.get_GEO(geo=geoid, destdir=destdir)
    return gse

def get_gpl(handle=None,id=None):
    """
    A function to return meta data from the given
    GSE object.

    :param id:int, List id of which you want to return the GPL.

    """
    if id:
        return handle.gpls.values()[id]
    else:
        return handle.gpls.values()[0]



def iterate(handle=None, type="gsm",anncol="Gene Symbol"):
    """
    Functions to iterate content of GEO object. 

    :param handle: a GSE handler.
    :param type: string, iterator type. 
    :param anncol: string, annotation column. Derived from GPL class.
    :returns: Iterator that holds GSM name and expression values data frame.
              Probe names and gene symbols are included.

    Usage:

    >>> for gsm_name, gsm_df in geo.iterate_gsm(handle=gse):
    >>>     print gsm_name
    >>>     print gsm_df.head()
    
    """
    if type=="gsm":
        gpl = get_gpl(handle=handle)
        for gsm_name, gsm in handle.gsms.items():
            # yield gsm_name, gsm
            annotation_column = anncol
            outdf = gsm.annotate(gpl, annotation_column)
            outdf = outdf[["ID_REF",annotation_column, "VALUE"]]
            yield gsm_name, outdf
    else:
        pass

def accumulate(handle=None, type="gsm",anncol="Gene Symbol",gpl_id=None):
    """
    Instead of iterating content of GEO object,
    we return one single Data frame with combined
    GSMs.

    :param handle: a GSE handler.
    :param anncol: string, annotation column. Derived from GPL class.
    :param type: string, iterator type. 
    :param gpl_id: int, index of GPL id you want to return.
    :returns: Data frame from multiple
              GSMs name and expression values data frame.
              Probe names and gene symbols are included.

    Usage:

    >>> full_df = geo.accumulate(handle=gse,type="gsm",anncol="Gene Symbol")
    
    """
    if type=="gsm":
        gpl = get_gpl(handle=handle,id=gpl_id)
        
        # print dir(gpl)
        # print gpl.table[["ID",anncol]].head()
        # print gpl.show_columns
        all_dfs = []
        annotation_column = anncol
        for gsm_name, gsm in handle.gsms.items():
            # yield gsm_name, gsm
            outdf = gsm.annotate(gpl, annotation_column)
            outdf = outdf[["ID_REF",annotation_column, "VALUE"]]
            outdf.columns = ["ID_REF",annotation_column, gsm_name]
            if outdf.empty:
                continue
            # print gsm_name
            # print outdf.head()
            all_dfs.append(outdf)
        
        merged_df = reduce(lambda ldf, rdf: pd.merge(ldf,rdf, \
                on=["ID_REF", annotation_column]), all_dfs).fillna("NoSymbol")

        # Keep only rows that has gene symbol
        unwanted = ["---","NoSymbol"]
        merged_df = merged_df[~merged_df[annotation_column].isin(unwanted)]
        return merged_df 
    else:
        pass
        

if __name__ == '__main__':
    main()

