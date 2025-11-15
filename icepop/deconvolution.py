#!/usr/bin/env python
""" 
Various functions to perform deconvolution of 
raw gene expression.

"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"
import pandas as pd
import numpy as np
from scipy.optimize import nnls
from scipy.optimize import minimize


def normalize_data(A, B, method='v3'):
    """
    Normalize matrices A and B using various methods.
    Refactored from inline code to improve maintainability.

    :param A: Matrix A (cell type proportion from ImmGen/IRIS)
    :param B: Vector B (mixed sample expression)
    :param method: Normalization method ('v1', 'v2', 'v3', 'v4')
    :returns: Tuple of (normalized_A, normalized_B)
    """
    if method == 'v1':
        # Normalize by maximum value across both A and B
        max_val = max(A.max(axis=0).max(), B.max())
        return A / max_val, B / max_val
    elif method == 'v2':
        # Normalize by sum
        return A / A.sum(axis=0), B / B.sum()
    elif method == 'v3':
        # Normalize by maximum (Best performing method)
        return A / A.max(axis=0), B / B.max()
    elif method == 'v4':
        # Normalize rows then columns (Better than V3 in some cases)
        new_A = A / A.sum(axis=1)[:, None]  # normalize row
        new_A = new_A / new_A.sum(axis=0)   # normalize column
        return new_A, B / B.sum()
    else:
        raise ValueError(f"Unknown normalization method: {method}")


def deconvolve_expr(mixedexp_df=None, cellpop_df=None,method=None,norm_method='v3'):
    """
    Perform deconvolution of raw expression data based on non-negative least-squares.

        .. math::
            
            Ax = B

    Here :math:`A` is the ImmGen/Iris matrix,
    :math:`B` is the user input sample expression where the concentration
    of each population is mixed. And :math:`x_{i}` is the value we are trying to
    estimate, that is the proportion of the cell type :math:`i` 
    in the sample user input.

    We obtain the solution by minimizing Euclidean 2-norm

        .. math::

            min_{x}(|| Ax - B ||^{2}), s.t. \\left\\{
                \\begin{array}{c l}      
                    \\sum_{i} x_{i} = 1\\\
                     x_{i} \\ge0, \\forall i
                \\end{array}\\right.

       

    which is essentially a `quadratic programming <https://en.wikipedia.org/wiki/Non-negative_least_squares>`_.

    We used the optimization method implemented in
    scipy `NNLS <http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.nnls.html>`_
    and scipy `SQLSP <http://docs.scipy.org/doc/scipy-0.16.0/reference/generated/scipy.optimize.minimize.html>`_.

    :param mixedexp_df: Pandas data frame of gene expression of single sample (:math:`B`).
    :param cellpop_df:  Pandas dataframe from ImmGen/IRIS cellpopulation (:math:`A`).
    :param method: string ('nnls','sqlsp')
    :param norm_method: string ('v1','v2','v3','v4') normalization method. Default 'v3'.

    :returns: a dictionary with deconvolved weight for every
              cell types (:math:`x`).
    """
    mixedexp_df = mixedexp_df.copy()
    cellpop_df = cellpop_df.copy()
    colnames = mixedexp_df.columns.values.tolist()
    sample_names = colnames[1:]
    mixedexp_df.columns = ["Genes"] + sample_names 
    mixedexp_df["Genes"] = mixedexp_df["Genes"].str.strip().str.upper()
    cellpop_df["Genes"] = cellpop_df["Genes"].str.strip().str.upper()
    nof_samples = mixedexp_df.shape[1] - 1

    # print mixedexp_df.head()
    # print cellpop_df.head()
   
    nof_celltypes = cellpop_df.shape[1] - 1
    cellpop_df_colnames = cellpop_df.columns.values.tolist()
    celltype_names = cellpop_df_colnames[1:]
    # print celltype_names

    # Ensure the input and cell proportion 
    # are aligned based on common genes.
    merged_df = reduce(lambda ldf, rdf: pd.merge(ldf,rdf, on="Genes"),
            [mixedexp_df,cellpop_df])


    B_series = merged_df[sample_names][sample_names[0]]
    A_df = merged_df[celltype_names]

    B = B_series
    A = A_df.values

    # print A
    # print B

    # Normalize the data so that they are comparable
    # Uses configurable normalization method (default: v3 - best performing)
    A, B = normalize_data(A, B, method=norm_method)


    if method == "sqlsp":
        x = by_sqlsp(A, B)
    else:
        x = by_nnls(A, B)

    # combine them
    x = x.tolist()
    final = dict(zip(celltype_names, x))

    # convert to pandas data frame
    final_df = pd.DataFrame.from_dict(final, orient="index")
    final_df.columns = [sample_names]
    final_df = np.round(final_df, decimals=10) 
    return final_df

def by_nnls(A=None, B=None):
    """
    Solution with Scipy NNLS. This method
    doesn't cater equalities explicitly. 
    Hence the solution will not sums up to 1 inherently. 
    We need to do normalization at the end.
   
    :param A: Pandas data frame (cell type proportion Immgen/IRIS)
    :param B: Pandas data frame (mixed samples)

    :return x: solution
    
    """
    # Add ones so that \sum x_{i} = 1
    # This is the last condition and must 
    # only be added at the end.
    number_of_celltypes = A.shape[-1]
    A = np.nan_to_num(A)
    B = np.nan_to_num(B)
    # print A.shape
    # print B.shape
    # print "B"
    # print repr(B)
    # print "A"
    # print repr(A)

    # A = np.vstack([A,np.ones(number_of_celltypes)])
    # B = np.hstack([B,1.0])

    # Tikhonov regularization
    lamb = 0.5
    n_variables = A.shape[1]
    A = np.concatenate([A, np.sqrt(lamb)*np.eye(n_variables)])
    B = np.concatenate([B, np.zeros(n_variables)])


    # NNLS
    x, rnorm = nnls(A,B)
    
    # Further ensure that the solution sums to 1
    # because earlier nnls doesn't give sum to 1
    x = x / x.sum()
    # print x.sum()
    return x
    

def by_sqlsp(A=None, B=None):
    """
    Method using  Sequential Least SQuares Programming (SQLSP).
    It allows the inclusion of equality constraint.
    Hence more exact than NNLS. The solution is guaranteed
    to sum to 1.

    :param A: Pandas data frame (cell type proportion Immgen/IRIS)
    :param B: Pandas data frame (mixed samples)

    :return x: solution
    
    """
    number_of_celltypes = A.shape[-1]
    A = np.nan_to_num(A)
    B = np.nan_to_num(B)
    def f(x):
        """ 

        This means:

        .. math::
            || Ax - B ||


        Produces better solution than

        .. math::
            || Ax - B ||^2

        """
        # return np.linalg.norm((A.dot(x) - B)^2)    
        return np.linalg.norm(A.dot(x) - B)    

    cons ={'type': 'eq',
           'fun': lambda x: sum(x) - 1}

    # initial guess
    x0 = [1] + [0] * (number_of_celltypes - 1)

    solution = minimize(f, x0, method='SLSQP', \
            bounds=((0, np.inf),)* number_of_celltypes, constraints=cons)
    x = solution.x
    return x
    
    

