ICEPOP (Immune CEll POPulation) is the method for estimating immune cell
population in the expressed genes. It enable analysis of differentially
expressed genes (DEGs) or raw expression data. These APIs and scripts will let
you have a fine grained control on the data analysis lacking in the web
version. On top of the core immune cell population deconvolution, it allows
you to download the raw data from NCBI GEO gene expression database,
normalize, and plot the data using command line interface.


Installation
============
ICEPOP Is best installed via `pip <https://pip.pypa.io/en/stable/>`_ (highly recommended) or 
`easy_install <https://wiki.python.org/moin/CheeseShopTutorial>`_ (older, but still works fine)::

    $ pip install git+https://github.com/ewijaya/icepop.git 

or:: 

    $ pip install git+git://github.com/ewijaya/icepop.git

 
Examples
========

Calculating immune response score from DEGs (in CSV, TSV or Excel file).
Make the output as a plot.


Alternative access 
==================
* `Web Application <https://sysimg.ifrec.osaka-u.ac.jp/icepop/>`_
* `API documentation <https://sysimg.ifrec.osaka-u.ac.jp/icepop/static//apidoc/html/index.html>`_

