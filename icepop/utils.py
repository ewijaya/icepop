#!/usr/bin/env python
""" 
Packages for various utilities.

"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
from collections import defaultdict
import re
import shutil
import os
import input_reader as ir

def copytree(src, dst, symlinks=False, ignore=None):
    """
    Copy source directory to target directory
    """
    for item in os.listdir(src):
        s = os.path.join(src, item)
        d = os.path.join(dst, item)
        if os.path.isdir(s):
            shutil.copytree(s, d, symlinks, ignore)
        else:
            shutil.copy2(s, d)

def cv_filter(degdf, organ='sp'):
    """
    Given list input DEG data and organ, 
    we throw genes that has FC lim < 1.

    :param degdf: Pandas data frame, differentially expressed genes (DEGs).
    :param organ: string('ln','lv','sp').

    :returns: new DEG data with pruned genes.
    """
    organ = organ.upper()
    cvfilt_df = ir.load_hdf_cvfilter()
    cvfilt_df = cvfilt_df[ cvfilt_df['Organ'] == organ ]
    cvfilt_df = cvfilt_df[ cvfilt_df['cv_filter'] < 1 ]
    wanted_genes = cvfilt_df["Genes"].tolist()
    degdf = degdf[ degdf['Genes'].isin(wanted_genes) ]
    return degdf

def mix_color(ctype_color_dict):
    """
    Given the RGB color. Mix
    them by taking average of R, G, B value
    """
    outdict = {}
    for ct1, rgb1 in ctype_color_dict.items():
        m = re.search('([\d+\,]+)',rgb1)
        rgb_1 = [int (x) for x in m.group(1).split(",")]
        for ct2, rgb2 in ctype_color_dict.items():
                
            m2 = re.search('([\d+\,]+)',rgb2)
            rgb_2 = [int (y) for y in m2.group(1).split(",")]
            mixed_rgb = [(rgb_1[i] + rgb_2[i])/2 for i,_ in enumerate(rgb_1)]
            mixed_rgb_str = "(" + ",".join([str(x) for x in mixed_rgb]) + ")"
            # print ct1, rgb1, ct2,rgb2, mixed_rgb_str
            outdict[ct1 + " " + ct2] = mixed_rgb_str
        
    return outdict
    

def get_sample_gene_list( genectype_df= None):
    """
    Get list of genes for each sample after taking top 100
    values from all of the samples.
    Input is a data frame that contain fc x cell type weight
    before summing up

    """
    colnames = genectype_df.columns.values.tolist()
    celltypes = colnames[1:-1]
    val_df = genectype_df[celltypes]
    allval_lst = val_df.values.flatten().tolist()
    allval_lst = list(set(allval_lst))
    allval_lst.sort(reverse=True)
    topk  = allval_lst[0:100]
    topk_vals = [ round(x,4)  for x in list(set(topk))]
    topk_vals.sort(reverse=True)
    # print json.dumps(topk_vals, indent=4)
    samples = list(set(genectype_df["Sample"].tolist()))

    outerdict = defaultdict(lambda: defaultdict(list))
    
    for sample in samples:
        idf = genectype_df[ genectype_df['Sample'] == sample ]
        idf = idf.copy()
        idf.drop('Sample',axis=1,inplace=True)
        genes_df = idf["Genes"]
        full_genelist =  genes_df.tolist()


        indone = {}
        # Optimized: Use itertuples() instead of iterrows() - 10-100x faster
        for row in idf.itertuples(index=False):
            cellvaldict = row._asdict()
            gene = cellvaldict["Genes"]
            for cell,val in  cellvaldict.items():
                if cell == "Genes": continue

                if round(val,4) in topk_vals:

                    if gene + " " + cell in indone: continue

                    ngene = gene[0] +  gene[1:].lower()
                    # print gene, ngene
                    outerdict[sample][cell].append(ngene)
                    indone[ngene + " " + cell] = 1
    
    return outerdict

if __name__ == '__main__':
    main()
