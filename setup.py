#!/usr/bin/env python
""" 


"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"

from distutils.core import setup

setup(name='icepop',
      version='0.1',
      author='Edward Wijaya',
      maintainer_email= 'ewijaya@gmail.com',
      install_requires=['pandas',
                        'sklearn',
                        'matplotlib',
                        'scipy',
                        'GEOparse',
                        'seaborn'
          
          ],
      package_data = {'proportion':[
                       'proportion_data/immgen_mouse.h5',
                       'proportion_data/iris_human.h5',
                       'proportion_data/ImmgenCons_all_celltypes_MicroarrayExp.csv',
                       'proportion_data/IRIS.csv'
                       ]},

      scripts=[ 'scripts/icepop_degs',
                'scripts/go_degs_cluster',
                'scripts/icepop_degs_cluster' ]
)
