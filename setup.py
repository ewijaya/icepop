#!/usr/bin/env python
""" 


"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"

from distutils.core import setup
import warnings

setup(name='icepop',
      #version.release.build
      version='1.4.2',
      author='Edward Wijaya',
      description='Python package for analysing immune cell population in expressed genes',
      license='MIT License',
      platforms=['linux','macosx'],
      author_email= 'ewijaya@gmail.com',
      url = 'https://github.com/ewijaya/icepop',
      keywords = ['bioinformatics','immunology','gene expression', 'microarray', 'reads'],
      install_requires=['pandas',
                        'sklearn',
                        'matplotlib',
                        'scipy',
                        'GEOparse',
                        'seaborn' ],
      packages=['icepop'],
      package_data={'icepop':[
                       'proportion_data/immgen_mouse_celltype_phenotype.h5',
                       'proportion_data/immgen_mouse_expr_by_phenotype.h5',
                       'proportion_data/immgen_mouse_organ.h5',
                       'proportion_data/immgen_mouse.h5',
                       'proportion_data/iris_human.h5',
                       'proportion_data/mouse_organ_genes.cvfilter.h5',
                       'proportion_data/ImmgenCons_all_celltypes_MicroarrayExp.csv',
                       'proportion_data/IRIS.csv'
                       'circos_conf/etc/*.*'
                       ]},
     scripts=[ 'scripts/icepop_degs',
               'scripts/go_degs_cluster',
               'scripts/icepop_degs_cluster',
               'scripts/icepop_pure' 
               ]

)
