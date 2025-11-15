#!/usr/bin/env python
""" 
Given a fold change threshold, 
we enumerate all the answers, per sample.
The output later will be used for 
JavaScript processing at the fron end.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015, The Cool Project"
import os
import json
import sys
import collections
import pandas as pd
import targetmine_query_urllib as tgm
import species_cell_proportion as scp
import cellpop_score as cp
import cluster_cellpop_score as ccp
import input_reader as ir
from collections import defaultdict

def enumerate_geneclust_go_output(cellpopdf, degdf, gene_count=False, \
        outfilename=None, k=None, logscale=None, fc_range=None,\
        method="ward", dist="euclidean",\
        species=None,pvalim=1,cormeth=None,immune=True,verbose=False,\
        tsv=False):

    """
    Enumerate all the gene cluster cell population score given 
    the range of fold change threshold. 
    This file output is later to be used to make JavaScript rendering, modeled according to:

    
    :param cellpopdf: Cell population data frame.
    :param degdf: DEGs data frame.
    :param fc_range: list of floats, range of fold change 
    :param logscale: boolean, transform fold change to logscale.
                      This does not affect the result substansially.
    :param gene_count: boolean, normalization method. 
                      If 'True', divide by sum of product. 
                      The total weight will have to result to 1.00.
    :param outfilename: output filename. The output type depends on the extension.
                    They are (".json", ".tsv")
    :param method: string('complete','average','ward')
    :param dist: string('euclidean','manhattan','pearsond')
    :param k: integer, number of cluster

    :param gene_ontology: boolean, whether or not to run GO.

    The following parameters only take effect when gene_ontology is set to true.

    :param species: string('mouse','human','rat') 
    :param useSymbol: string('true','false'), whether query is done using gene
                      symbol, otherwise ProbeID
    :param pvalim: string, lowerbound of P-value for a GO term to be
                      displayed.
    :param cormeth: string('HOLM_BONFERRONI','BENJAMINI_HOCHBERG','BONFERRONI')
    :param immune: boolean, show only immune related GO terms.
    :param verbose: boolean, display processing stage.

    :returns: JSON format output for both GO as list of list. 


    Usage:

    >>> from icepop import enumerate_output as eo
    >>> out_go_json = "input_cluster_type1_degs.large.go.json" 
    >>> eo.enumerate_geneclust_go_output(cellpop_df, indf,
            >>> fc_range=foldchange_range, outfilename = out_go_json,
            >>> gene_count=True, logscale=False,k=nof_clust, method="ward",
            >>> dist="euclidean",
            >>> species='mouse',pvalim=1.00,cormeth="HOLM_BONFERRONI",immune=True)

    """
    outerdict = defaultdict(list)
    go_lol = []
    for fc_lim in fc_range:

        sample_names = degdf.columns.values.tolist()
        sample_names = sample_names[2:]
        selected_degdf = degdf[(degdf[sample_names]>fc_lim).any(1)]
        nof_final_samples = selected_degdf.shape[0]


        # This condition is important, to
        # prevent crash at clustering.py on AgglomerativeClustering
        # because final genes count is less than number of cluster.
        if nof_final_samples < k:
            if verbose:
                sys.stderr.write("Skip clustering for fold change:" + str(fc_lim) +  "\n")
            continue
            
        if verbose:
            sys.stderr.write("Clustering for foldchange: " + str(fc_lim) + "\n")
            sys.stderr.write("Number of genes after fold change filter:" +  selected_degdf.shape[0] + "\n")

        full_clust_cpopdf, full_clust_degdf = \
        ccp.cluster_cellpop_score(cellpopdf=cellpopdf, degdf=degdf,\
                fclim=fc_lim,gene_count=gene_count,logscale=False, k=k,\
                method=method, dist=dist)

        #------------------------- 
        # GO for cluster 
        #------------------------- 
        # sys.stderr.write("Process GO analysis for GO\n")
        tmp_fcddf = full_clust_degdf.copy()
        tmpcn = tmp_fcddf.columns.values.tolist()

        # rename column for consistency
        tmp_fcddf.columns._data[1] = 'probe'
        tmp_fcddf.columns._data[2] = 'gene'

        tmp_fcddf = tmp_fcddf[["ClusterID","gene"]]
        for clusterid, tdf in tmp_fcddf.groupby(['ClusterID']):
            genelist = tdf['gene'].tolist()
            # Make sure the list only contain string
            genelist = list(filter(lambda x: isinstance(x, str), genelist))
            genelist_str =  ",".join(genelist)
            # print clusterid, genelist_str
            for vals in tgm.get_tgm_data(genelist_str,immune=immune, species=species, \
                    useSymbol="true",pvalim=pvalim,cormeth=cormeth):

                theme, term, pval, genelist = vals
                nvals = [clusterid, fc_lim, term, pval, genelist]
                go_lol.append(nvals)
                # print "\t", nvals

    headers = ["Cluster ID","Fold Change","GO Term","P-value","Gene List"]
    if outfilename == None:
        outdf = pd.DataFrame(go_lol,columns=headers)
        return outdf
    else:
        _, extension = os.path.splitext(outfilename)
        if extension==".tsv":
            if verbose:
                 sys.stderr.write("Creating TSV file\n")
            outdf = pd.DataFrame(go_lol,columns=headers)
            # print outdf.head()
            outdf.to_csv(outfilename,sep="\t",index=False,float_format='%.3f')
        else:
            if verbose:
                 sys.stderr.write("Creating JSON file\n")
            # Write cell population cluster to JSON
            with open(outfilename,'w') as jsonout:
                json.dump(go_lol, jsonout, indent=4)
        

def enumerate_geneclust_output(cellpopdf, degdf, gene_count=False, \
        outfilename=None, k=None, logscale=None, fc_range=None,\
        method="ward", dist="euclidean"):

    """
    Enumerate all the gene cluster cell population score given 
    the range of fold change threshold. 
    
    :param cellpopdf: Cell population data frame.
    :param degdf: DEGs data frame.
    :param fc_range: list of floats, range of fold change 
    :param logscale: boolean, transform fold change to logscale.
                      This does not affect the result substansially.
    :param gene_count: boolean, normalization method. 
                      If 'True', divide by sum of product. 
                      The total weight will have to result to 1.00.
    :param outfilename: output filename.
    :param method: string('complete','average','ward')
    :param dist: string('euclidean','manhattan','pearsond')
    :param k: integer, number of cluster


    :returns: JSON format output for both cell population.

    Usage:

    >>> from icepop import enumerate_output as eo
    >>> out_json = "input_cluster_type1_degs.large.json" 
    >>> nof_clust        = 15
    >>> foldchange_range = [1.5,2,2.5,3,3.5,4,4.5,5]
    >>> eo.enumerate_geneclust_output(cellpop_df, indf,
            >>> fc_range=foldchange_range, outfilename = out_json,
            >>> gene_count=True,  logscale=False,k=nof_clust, method="ward",
            >>> dist="euclidean")

    """

    lod = []
    fc_range.sort()
    for fc_lim in fc_range:
        # print fc_lim
        outerdict = defaultdict(list)
        sample_names = degdf.columns.values.tolist()
        sample_names = sample_names[2:]
        selected_degdf = degdf[(degdf[sample_names]>fc_lim).any(1)]
        nof_final_samples = selected_degdf.shape[0]


        # This condition is important, to
        # prevent crash at clustering.py on AgglomerativeClustering
        if nof_final_samples < k:
            continue
            
        full_clust_cpopdf, full_clust_degdf = \
        ccp.cluster_cellpop_score(cellpopdf=cellpopdf, degdf=degdf,\
                fclim=fc_lim,gene_count=gene_count,logscale=False, k=k,\
                method=method, dist=dist)

        #------------------------- 
        # Cell population  score for cluster
        #------------------------- 
        tmp_df = full_clust_cpopdf.copy()
        tmp_df.set_index("Cluster", drop=True, inplace=True)
        mdict = tmp_df.to_dict()
        tmp_list  = []
        for celltype, clustval in mdict.items():
            for clustid, val in clustval.items():
                xy_name = celltype + ";" + str(clustid)
                tmp_list.append({'celltype;clustid':xy_name, 'value':{'zscore':val}})
        lod.append({'threshold':fc_lim, 'heatdata':tmp_list})

    # Write cell population cluster to JSON
    with open(outfilename,'w') as jsonout:
        json.dump(lod, jsonout, indent=4)
    
    return


def enumerate_output(cellpopdf, degdf, gene_count=False, \
        outfilename=None, logscale=None, fc_range=None):
    """
    Enumerate all the cell population score given 
    the range of fold change threshold.
    
    :param cellpopdf: Cell population data frame.
    :param degdf: DEGs data frame.

    :param fc_range: list of floats, range of fold change 
    :param logscale: boolean, transform fold change to logscale.
                  This does not affect the result substansially.
    :param gene_count: boolean, normalization method. 
                  If 'True', divide by sum of product. 
                  The total weight will have to result to 1.00.
    :param outfilename: output filename.

    :returns: JSON format output.

    Usage:

    >>> from icepop import enumerate_output as eo
    >>> foldchange_range =[1.5,2,2.5,3,3.5,4,4.5,5]
    >>> out_json = "input_type1_degs.large.json"
    >>> eo.enumerate_output(cellpop_df, indf, fc_range = foldchange_range, outfilename = out_json,  gene_count=True, logscale=False)
    
    """

    final = []
    for fc_lim in fc_range:
        # sys.stderr.write("EO:" +  str(gene_count) + "\n")
        nof_genes_dict, sample_response_dict, celltype_response_dict, cpop_score_df = cp.deg_cellpopscore_df(cellpopdf, degdf, fclim=fc_lim, \
                            gene_count=gene_count, logscale=logscale)
        mdict = cpop_score_df.to_dict()
        outlist = []
        for sample, cellpopvals in mdict.items():
            celltypelist            = []
            nof_genes               = nof_genes_dict[sample]
            celltype_response_thres = float(celltype_response_dict[sample])
            sample_response_val   = float(sample_response_dict[sample])


            # print sample
            # print "CRT : ", celltype_response_thres
            # print "SR  : ", sample_response_val

            # sort by cell types
            od = collections.OrderedDict(sorted(cellpopvals.items()))
            for celltype, val in od.items():
                celltypelist.append({"celltype":celltype,"score":val})
            outlist.append({
                            "values":celltypelist, 
                            "sample": sample,
                            "nof_genes": nof_genes, 
                            "celltype_response_thres":celltype_response_thres,
                            "sample_response_score":sample_response_val
                            })
        
        
        final.append({"threshold":fc_lim,"histograms":outlist}) 
    

    if outfilename == None:
        pass
    else:
        with open(outfilename,'w') as jsonout:
            json.dump(final, jsonout, indent=4)
    return

# Test code moved to tests/ directory
