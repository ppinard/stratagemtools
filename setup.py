#!/usr/bin/env python

# Standard library modules.
import codecs

# Third party modules.
from setuptools import setup, find_packages

# Local modules.
import versioneer

# Globals and constants variables.


# Get the long description from the relevant file
with codecs.open('README.rst', encoding='utf-8') as f:
    long_description = f.read()

setup(name="stratagemtools",
      version=versioneer.get_version(),
      description="Python interface to SAMx STRATAGem",
      long_description=long_description,

      author="Philippe T. Pinard",
      author_email="philippe.pinard@gmail.com",

      url='https://github.com/ppinard/stratagemtools',
      license="MIT",
      keywords='microscopy microanalysis thin-film stratagem samx',

      classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: Microsoft :: Windows',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Physics',
        ],

      packages=find_packages(),

      install_requires=['pyparsing'],

      test_suite='nose.collector',

      cmdclass=versioneer.get_cmdclass(),
)

