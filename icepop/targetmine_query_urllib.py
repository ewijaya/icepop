#!/usr/bin/env python
""" 

Programmatic way to do GO analysis.
It uses `TargetMine <http://targetmine.mizuguchilab.org/>`_.  
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2014, The Cool Project"
import sys
import os
import re
import requests
import decimal


def get_tgm_data(genelist_string,species="mouse",\
        immune=True,\
        useSymbol="true",pvalim="1.00",cormeth="HOLM_BONFERRONI"):
    """
    Get GO enrichment given list of genes.

    :param genelist_string: string, comma separated (e.g. "Ssb,Eny2,Hsp90b1,Ube3a,Cry1,Prkaa1").
    :param species: string('mouse','human','rat'). 
    :param useSymbol: string('true','false'), whether query is done using gene symbol, otherwise ProbeID.
    :param pvalim: string('1.00','0.01'), lowerbound of P-value for a GO term to be displayed.
    :param cormeth: string('HOLM_BONFERRONI','BENJAMINI_HOCHBERG','BONFERRONI')
    :param immune: boolean, show only immune related GO terms.

    :returns: Generator for GO with enrichment Pvalues

    """

    theme_dict = { #"Pathway": 'All data set',  
                   "Gene Ontology":"Biological_Process",
                   #"GO Slim":"Biological_Process",
                   #'KEGG Pathway',
                   # "Integrated Pathway Cluster": ""
                 }
    
    immune_go_list = ['stress','cell death','wounding','inflammatory',
               'migration', 'cytokine', 'defense response','chemotaxis', 
               'immune', 'interferon', 'interleukin']



    taxonid_dict = { 'human': '9606',
                     'mouse': '10090',
                     'rat'  : '10116'}

    # ========= Variables, you can modify them accordingly ==========
    enrichment_url = "targetmine.mizuguchilab.org"
    cut_off = pvalim
    taxon_id = taxonid_dict[species]
    correctionMethod = cormeth
    showBg = "false"
    useSymbol = useSymbol
    # username = "tmext"
    # password = "saito27x"

    for theme,data_set in theme_dict.iteritems():
        # mainurl = "http://"+ enrichment_url + "/enrich/analysis/enrichment"  
        mainurl = "http://"+ enrichment_url + "/pipeline/analysis/enrichment"  

        payload = {'theme' : theme, 
                   'dataSet' : data_set,
                   'taxonId' : taxon_id,
                   'correctionMethod' : correctionMethod,
                   'cutOff'  : cut_off,
                   'useSymbol' : useSymbol,
                   'list' : genelist_string,
                   'showBg' : showBg}



        # r = requests.post(mainurl,auth=(username,password),data=payload)
        r = requests.post(mainurl,data=payload)
        for res in r.text.split("\n"):
            tmp =  res.split("\t")
            if len(tmp) != 4: continue
            term = tmp[1]
            pval = tmp[2]
            genelist = tmp[3]


            # Format p-value into scientific style
            pval = float(pval)
            pval = "{:.2E}".format(decimal.Decimal(pval))
            pval = str(pval)

            m = re.search('\[(.*)\]',genelist)
            genelist = m.group(1)

            if immune and inlist(immune_go_list, term):
                yield theme, term, pval, genelist
            elif immune == False:
                yield theme, term,pval,genelist


    


def inlist(substring_list, query):
    """
    Check if the sustring in the list is contained
    in the string.
    """
    return any(substring in query for substring in substring_list)
    

def main():
    """
    Use for testing this file.
    """
    correction = "BONFERRONI"
    genelist_str = "Ssb,Eny2,Hsp90b1,Ube3a,Cry1,Prkaa1"
    pvalim = 1
    species = 'mouse'
    for vals in get_tgm_data(genelist_str,immune=True, species=species, useSymbol="true",pvalim=pvalim,cormeth=correction):
        print "\t".join(vals)
        

if __name__ == '__main__':
    main()
