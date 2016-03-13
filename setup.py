#!/usr/bin/env python

# Script information for the file.
__author__ = "Philippe T. Pinard"
__email__ = "philippe.pinard@gmail.com"
__version__ = "0.1"
__copyright__ = "Copyright (c) 2013 Philippe T. Pinard"
__license__ = "GPL v3"

# Standard library modules.
import os
import re
import codecs

# Third party modules.
from setuptools import setup, find_packages

# Local modules.

# Globals and constants variables.
BASEDIR = os.path.abspath(os.path.dirname(__file__))

def find_version(*file_paths):
    """
    Read the version number from a source file.

    .. note::

       Why read it, and not import?
       see https://groups.google.com/d/topic/pypa-dev/0PkjVpcxTzQ/discussion
    """
    # Open in Latin-1 so that we avoid encoding errors.
    # Use codecs.open for Python 2 compatibility
    with codecs.open(os.path.join(BASEDIR, *file_paths), 'r', 'latin1') as f:
        version_file = f.read()

    # The version line must have the form
    # __version__ = 'ver'
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

# Get the long description from the relevant file
with codecs.open('README.rst', encoding='utf-8') as f:
    long_description = f.read()

setup(name="stratagemtools",
      version=find_version('stratagemtools', 'stratagem.py'),
      description="Interface to STRATAGem",
      long_description=long_description,

      author="Philippe T. Pinard",
      author_email="philippe.pinard@gmail.com",

      url='http://stratagemtools.bitbucket.org',
      license="GPL v3",
      keywords='microscopy microanalysis thin-film stratagem samx',

      classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: Microsoft :: Windows',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Topic :: Scientific/Engineering :: Physics',
        ],

      packages=find_packages(),
      zip_safe=False,
      
      install_requires=['pyparsing'],

      test_suite='nose.collector',
)

