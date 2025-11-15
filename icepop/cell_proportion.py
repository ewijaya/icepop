#!/usr/bin/env python
"""
.. _cell_proportion-label:

These are classes to process various cell population data.
The cell population data include:
   
* `ImmGen <http://www.immgen.org/>`_ for mouse.
* `IRIS <http://share.gene.com/share/clark.iris.2004/iris/iris.html>`_  for human.

The above files are downloaded as public database.
It's stored in ``proportion_data/`` directory
See doctsring below for detailed explanation.

"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import csv
import json
import gc
import numpy as np 
import pandas as pd
# pd.set_option('expand_frame_repr',False)
from itertools import islice
from collections import defaultdict

class Iris:
    """
    These are classes to process IRIS data.
    The default input is ``ImmgenCons_all_celltypes_MicroarrayExp.csv``
    For each gene compute the weight of *K* cell types. 
    Weight (in percentage) are the proportion of average gene 
    expression of a given cell type over all cell types.

    >>> from icepop import cell_proportion as cellprop
    >>> iris = cellprop.ImmGen(irisfile)
    >>> iris_weight_dict, iris_ctlist = iris.ComputePercentageWeight()
    """
    def __init__(self, irisfn):
        self.irisfn = irisfn

    def ComputeAverage(self):
        """
        Compute average for subset of each cell type.
        """
        df = pd.read_csv(self.irisfn, sep=",", comment="#",
                skiprows=[2])
        colnames = df.columns.values.tolist()
        ncolnames = [x  for x in colnames if "DUMMY" in x]
        df.drop(ncolnames,axis=1,inplace=True)

        # Averaging here
        df_average = df.set_index(["Probes","Genes"])
        # Optimized: Pre-compute column mapping instead of lambda
        col_map = {col: col.split(".")[0] for col in df_average.columns}
        df_average = df_average.groupby(col_map, axis=1).mean()

        return df_average
        

    def ComputePercentageWeight(self):
        """
        Same with ImmGen version but here we do it
        using Pandas. Hence the output here is a data frame
        instead of nested dictionary.
        """
        #-------------------------------------------
        #  Compute weight
        #------------------------------------------- 
        df_weight = self.ComputeAverage()
        celltype_list  = df_weight.columns.values.tolist()
        # Sum values across cell type 
        df_weight["sum"] = df_weight.sum(axis=1)

        # Divide each value for every column with sum
        df_weight = df_weight.loc[:,celltype_list].div(df_weight["sum"],axis=0)
        # Multiply every row with 100 as percentage
        df_weight = df_weight.multiply(100, axis=0) 

        # Check if total sum is equal to 1
        # df_weight["sum"] = df_weight.sum(axis=1)

        df_weight = df_weight.reset_index()
        df_weight.drop('Probes',axis=1,inplace=True)
        return df_weight, celltype_list
        
        

class ImmGen:
    """
    These are classes to process ImmGen data.
    The default input is ``ImmgenCons_all_celltypes_MicroarrayExp.csv``.
    For each gene compute the weight of 10 cell types. 
    Weight (in percentage) are the proportion of average gene 
    expression of a given cell type over all cell types.

    >>> from icepop import cell_proportion as cellprop
    >>> ig = cellprop.ImmGen(igfile)
    >>> ig_weight_dict, ig_ctlist = ig.ComputePercentageWeight()
    """
    def __init__(self, immg_fn):
        """ Define parent constructor for ImmGen"""
        self.immg_fn = immg_fn
    
    def GetGeneCelltypeDoList(self):
        """Parse ImmGen data and then collect them as
        dictionary of list. Each dictionary of list is constructed
        with gene names and cell types as keys and 
        list of gene expressions as values. 
        """
        celltypes=[];
        outerdict = defaultdict(lambda: defaultdict(list))
            
        # Parse ImmGen data
        with open(self.immg_fn,'rb') as csvfile:
            csvreader = csv.reader(csvfile,skipinitialspace=True)
            celltypes = next(csvreader,[])[3:]
            # next two rows can be skipped
            next(islice(csvfile,2,2),None)
            for row in csvreader:
                name = row[1]
                for celltype,value in zip(celltypes,row[3:]):
                    # Collect them into DoL
                    outerdict[name][celltype].append(float(value))
        return (celltypes,outerdict)

    def ComputeAverageByOrgans(self, weighted=False):
        """
        Group by organs. Then average for for every cell type.
        """
        infile = self.immg_fn
        df = pd.read_csv(infile,header=None,index_col=[1,2],skiprows=[2],dtype='unicode').iloc[:,1:]
        df.columns = pd.MultiIndex.from_arrays(df.iloc[:2].values)
        df = df.iloc[2:].astype(float)
        df.index.names = ['cell', 'organ']
        df = df.reset_index('organ', drop=True)
        result = df.groupby(level=[0, 1], axis=1).mean()
        result = result.stack().replace(np.nan, 0).unstack()
        result = result.swaplevel(0,1, axis=1).sort_index(axis=1)

        if weighted:
            # Weighted them for 0 to 1 scale. 
            sum_org =  result.sum(level=0, axis=1)
            result = result.div(sum_org,axis=1,level=0)*100
            # print result.loc[:,'thymus'].head()
            return result
        else:
            return result
        


    def ComputeAverage(self):
        """
        Note that, a given cell type may contain 
        many subtypes (e.g. StemCells contains:
        SC_LT34F_BM, SC_LTSL_BM, SC_STSL_BM, SC_LTSL_FL,... etc).
        Here we compute the average of gene expression
        in the given cell type from its subtypes.

        Let :math:`K` be the list of subtypes :math:`k` in celltype :math:`j`.
        And :math:`e_{ik}` is the expression of gene :math:`i` in 
        subtype :math:`k`. The average expression :math:`\\bar{e}_{ij}` is computed as
        follows: 

            .. math::
                
                \\bar{e}_{ij} = \\frac{\\sum_{k=1}^{K} {e_{ik}}}{|K|}

        """
        outerdict2 = defaultdict(dict)
        ctlist, outdict = self.GetGeneCelltypeDoList()
        for gene,ctvals in outdict.items():
            for ct,vals in ctvals.items():
                avg = "%.2f"  % (sum(vals)/len(vals))
                outerdict2[gene][ct] = avg


        return(outerdict2, ctlist)

    def _compute_sum_dict(self, outdict2):
        """
        Helper: Compute sum of values for each gene across all cell types.
        Extracted to reduce duplication.
        """
        return dict((k, sum(float(f) for f in v.values()))
                    for k, v in outdict2.items())

    def _get_unique_sorted_celltypes(self, celltypelist):
        """
        Helper: Get unique sorted cell type list.
        Extracted to reduce duplication.
        """
        uctlist = list(set(celltypelist))
        uctlist.sort()
        return uctlist

    def ComputePercentageWeight(self):
        """
        After we obtain the average of cell types' gene expression,
        here given a gene, we compute the weight of each gene
        in their respective cell type. Let :math:`e_{ij}` be the expression
        of gene :math:`i` in cell type :math:`j`.
        The cell type weight for every gene is computed as follows:

            .. math::

                w_{ij} = \\frac{\\bar{e}_{ij}}{\\sum_{j=1}^{C} { \\bar{e}_{ij} }}

        """
        # Sum up value for each gene for all cell type
        outdict2, celltypelist = self.ComputeAverage()
        uctlist = self._get_unique_sorted_celltypes(celltypelist)
        sumdict = self._compute_sum_dict(outdict2)

        # Compute percentage
        outerdict3 = defaultdict(dict)
        for gene,ctavg in outdict2.items():
            totavg = float(sumdict[gene])
            for ct,avg in ctavg.items():
                perc = "%.2f" % ((float(avg) / totavg) * 100)
                outerdict3[gene][ct] = float(perc)
        return(outerdict3, uctlist)

    def ComputePercentageWeightBracket(self):
        """
        Similar with ``ComputePercentageWeight()``
        Only output them in the form of dictionary of list.
        List contains cell type with bracketted weight,
        (e.g. abTcells(26.12)).
        """
        # Sum up value for each gene for all cell type
        outdict2, celltypelist = self.ComputeAverage()
        uctlist = self._get_unique_sorted_celltypes(celltypelist)
        sumdict = self._compute_sum_dict(outdict2)

        # Compute percentage
        outerdict3 = defaultdict(list)
        for gene,ctavg in outdict2.items():
            totavg = float(sumdict[gene])
            ctavg_sorted = sorted(ctavg.items(),
                           key=lambda x: float(x[1]), reverse=True)
            for ct,avg in ctavg_sorted:
                perc = ct + "(%.2f)" % ((float(avg) / totavg) * 100)
                outerdict3[gene].append(perc)
        return(outerdict3, uctlist) 

def main(igfile):
    """Use for testing this file.

    :param igfile: text file, cell population proportion.

    To test this code:: 

        python cell_proportion.py 
    """
    ig = ImmGen(igfile)
    od3 = ig.ComputePercentageWeight()
    # print json.dumps(od3, indent=4)

if __name__ == '__main__':
    immgen_file = 'proportion_data/ImmgenCons_all_celltypes_MicroarrayExp.csv'
    main(immgen_file)
