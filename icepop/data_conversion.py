#!/usr/bin/env python
""" 
Various data conversion functions.

"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import csv
import json
import os
import collections

def df_to_dol(df):
    """
    Function to convert dataframe into 
    JSON format. Later to be readable for Plottable.js

    Used especially for deconvolutin pie chart.
    """
    s = df.to_json(orient='columns')
    sdict = json.loads(s)
    alldict =[]
    for k,v in sdict.iteritems():
        sampledict = {}
        sampledict['sample'] = k
        sampledict['pies_pct'] = []


        # sort by cell types
        od = collections.OrderedDict(sorted(v.items()))
        for ct, pct in od.iteritems():
            tmpdict = {}
            tmpdict['score'] = pct
            tmpdict['celltype'] = ct

            sampledict['pies_pct'].append(tmpdict)
        alldict.append(sampledict)

    return alldict

