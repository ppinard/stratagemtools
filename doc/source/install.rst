Installation guide
==================

Requirements
------------

*stratagemtools* requires:

  - Microsoft Windows
  - STRATAGem >= 4.1.0, stratadll.dll >= 4.7.0
  - Python 3.x (32-bit, does not work with 64-bit)

and the following Python libraires, which should be automatically installed by
pip.

  - pyparsing
  - nose (for testing only)

Installation of release version
-------------------------------

Run in a command prompt::

  py -3 -m pip install stratagemtools

Installation of developer version
---------------------------------

Run in a command prompt::

  git clone https://github.com/ppinard/stratagemtools.git
  cd stratagemtools
  py -3 -m pip install -e .
 
