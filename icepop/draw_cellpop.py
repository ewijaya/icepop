#!/usr/bin/env python
from __future__ import division
""" 
Given a data frame that looks like this:

                    FoldChange_Sample2  FoldChange_Sample3
    Bcells                    0.041977            0.059211
    DendriticCells            0.071067            0.076072
    ...
    StemCells                 0.048593            0.065881
    StromalCells              0.108675            0.189875
    abTcells                  0.059072            0.067225
    gdTCells                  0.062106            0.066489


We perform  various plotting.
"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015."
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
matplotlib.style.use('ggplot')
import seaborn as sns
sns.set(style="white")

brewer_Dark2_ten = ["#1B9E77", "#AE6D1C", "#A16864", "#9B58A5", "#D8367D", 
                    "#749829", "#BBA90B", "#C9930D", "#97722D", "#666666"]

brewer_Accent_ten = ["#7FC97F", "#B0B4C1", "#E1B8A8", "#FDD58C", "#E8EE9B",
                     "#4E7CAD", "#B2258F", "#DA2950", "#AB5D28", "#666666"]

brewer_Paired_ten = ["#A6CEE3", "#3F8EAA", "#79C360", "#B89B74", "#E52829", 
                     "#FDB762", "#ED8F47", "#9471B4", "#DDD399", "#B15928"]

tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]   


# Scale the RGB values to the [0, 1] range, which is the format matplotlib
# accepts.    
for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.)   


def bar_plot(indf, title=None,outfile=None):
    """
    Make barplot. All in one figure. Not very pretty.
    Based on `Pandas method
    <http://pandas.pydata.org/pandas-docs/version/0.15.0/visualization.html#histograms>`_.
    
    :param indf:    cell population data frame.
    :param title:   string.
    :param outfile: string (e.g. 'plot.png').

    :returns: Histogram plot.
        Image file type is defined automatically from output file extension.
    """

    toplot_df = indf
    barplot   = toplot_df.plot(kind="bar", figsize=(15,15), \
                          color = brewer_Paired_ten, \
                          width=0.7,\
                          rot=45, fontsize = 20,\
                          legend=False,
                          subplots=False)

    plt.xlabel("")
    plt.ylabel("Score", fontsize=25, fontweight="bold")
    fig = barplot.get_figure()
    fig.suptitle(title,fontsize=30,fontweight="bold")
    fig.savefig(outfile)
    plt.close()
    return

def condn_plot(condn_mat,outfile=None):
    """
    A function to plot condition number
    with respect to series of specificity threshold.

    :param condn_mat: a numpy matrix generated by :func:`condn_thres_mat() <specificity.condn_thres_mat>`.
    :returns: line plot.
    """
   
    sns.set(style="darkgrid")
    spec_lim = condn_mat[:,0]
    condn_vec = condn_mat[:,1]
    minid = condn_mat.argmin(0)[1]
    min_lim, min_connum, nof_marker_genes = condn_mat[minid]


    # Actual plotting
    # Color from Tableau 20 scheme
    maincol = "#D12B31"
    axcol   = "#42B1AF"
    plt.plot(spec_lim, condn_vec, lw=3.5, color=maincol)
    plt.xlabel('Threshold')
    plt.ylabel('Condition Number')
    plt.axhline(y=min_connum,linewidth=1.2,linestyle='--',color=axcol)
    plt.axvline(x=min_lim,linewidth=1.2,linestyle='--',color=axcol)
    plt.plot([min_lim],[min_connum],'o')
    in_text = str(nof_marker_genes.astype(int)) + ' marker genes' 
    plt.text(min_lim - 0.05, min_connum + 5, in_text.encode('string-escape') )
    plt.suptitle('Best Specificity Threshold', fontweight='bold')

    # Save to file
    plt.savefig(outfile)
    plt.close()

    return 


def bar_plot_facet(indf, title=None, outfile=None):
    """
    Make barplot in trellis (facet), using `Seaborn <http://stanford.edu/~mwaskom/software/seaborn/>`_.
    Each sample one figure.

    :param indf:    cell population data frame.
    :param title:   string.
    :param outfile: string (e.g. 'plot.png').

    :returns: Histogram plot facet.
        Image file type is defined automatically from output file extension.

    Usage:

    >>> from icepop import draw_cellpop as dc
    >>> outfile = "barplot.png"
    >>> dc.bar_plot_facet(cpop_score_df, outfile=outfile)

    """

    sample_names = indf.columns.values.tolist()

    tmp_df = indf
    tmp_df["Celltype"] = tmp_df.index
    ndf = pd.melt(tmp_df,id_vars=["Celltype"], var_name="Samples", \
                  value_name="Score", value_vars=sample_names )

    # Initialize a grid of plots with an Axes for each walk
    # Use factorplot
    grid = sns.factorplot(x="Celltype", y="Score", col="Samples", size=10, 
            palette="Set3", kind="bar", saturation=.8, data=ndf, col_wrap=5,
            legend=True)
    (grid.set_titles("{col_name}", name="Helvetica", va="top", weight="bold",size=30)
        .set_ylabels(size=25)
        .set_axis_labels("","Score")
        .set_xticklabels(size=20) #[] we hide the x-ticks label
        .set_yticklabels(size=20)
        .despine(bottom=True))

    #print dir(grid._legend_out)
    for ax in grid.axes:
        plt.setp(ax.get_xticklabels(), rotation=45)
    grid.savefig(outfile)
    plt.close()
    return

def cellpop_cluster_heatmap(cellpop_clust_df, fig_width=5.25, \
        fig_height=10.25, title=None, outfile=None):
    """
    Make heatmap for cell population clusters, using `Seaborn <http://stanford.edu/~mwaskom/software/seaborn/>`_.

    :param cellpop_clust_df: cell population cluster data frame.
    :param title: string
    :param outfile: string (e.g. 'plot.png').

    :returns: Heatmap, with cluster number as *y*-axis and cell population as *x*-axis.
        Image file type is defined automatically from output file extension.

    Usage:

    >>> from icepop import draw_cellpop as dc
    >>> dc.cellpop_cluster_heatmap(full_clust_cpopdf, outfile="heatmap.png")

    """
    df = cellpop_clust_df.copy()
    df.set_index("Cluster",drop=True,inplace=True)

    nof_row = len(df)
    fig_height_ratio = nof_row / 50 
     
    cmap = sns.diverging_palette(220, 20, sep=5, as_cmap=True)
    grid_kws = {"height_ratios": (1.5, .05), "hspace": .4}
    fig, (ax, cbar_ax) = plt.subplots(2, gridspec_kw=grid_kws)
    ax = sns.heatmap(df, square=False, cmap=cmap,
                     ax = ax, cbar_ax=cbar_ax,
                     cbar_kws={"orientation": "horizontal"}
            )
    plt.setp(ax.get_xticklabels(), rotation=90)
    cbar_ax.set_xlabel('Z-score')


    # Because it gets ugly when row is greater than 50
    # Then we need to resize it accordingly
    if nof_row > 50:
        fig_height = fig_height * fig_height_ratio
        
    fig.set_size_inches(fig_width,fig_height)
    fig.savefig(outfile)
    plt.close()
    return

def pieplot_deconv(indf,nrows=3, ncols=4, outfile=None):
    """
    Plot pie chart for deconvolution of
    gene expression as facet

    :param indf:    cell population data frame.
    :param nrows:   integer
    :param ncols:   integer
    :param outfile: string (e.g. 'plot.png').
    """
    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(12,12))
    for ax in axes.flat:
        ax.axis('off')

    label_incl_lim = 0.3
    for ax, col in zip(axes.flat, indf.columns):
        data = indf[col]
        labels = [n if v > data.sum() * label_incl_lim else '' for n, v in
                zip(indf.index,data)]
        
        ax.pie(data, labels=labels, autopct=choose_autopct(label_incl_lim*100),  colors=tableau20)
        ax.set(ylabel='',  aspect='equal' )
        ax.set_title(col,fontweight='bold')

    colnames = indf.columns.values.tolist()
    nof_samples = len(colnames)

    legend_row_id = (nrows * ncols) % nof_samples

    legend_col_id = None
    if nof_samples < ncols:
        legend_col_id = nof_samples - 1
    else:
        # Not perfect need to be tested
        legend_col_id = ncols - ((nrows * ncols) % nof_samples ) -1 



    # for 10 samples
    # axes[2, 1].legend(indf.index, bbox_to_anchor=(3.0,1.5),fontsize=15)
    axes[legend_row_id, legend_col_id].legend(indf.index, bbox_to_anchor=(2.5,1.5),fontsize=15)

    fig.savefig(outfile,dpi=300)
    plt.close()
    return 


def choose_autopct(limit):
    """
    Function to select autopct based on threshold. 
    """
    def inner_autopct(pct):
        return ('%.2f' % pct) if pct > limit else ''
    return inner_autopct
    
    
