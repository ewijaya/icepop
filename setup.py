#!/usr/bin/env python
""" 


"""
__author__    = "Edward Wijaya <ewijaya@gmail.com>"
__copyright__ = "Copyright 2015"

from distutils.core import setup

setup(name='icepop',
      #major.minor.release.build
      version='1.0.0.0',
      author='Edward Wijaya',
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

)
