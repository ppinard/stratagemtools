#!/usr/bin/env python
"""
================================================================================
:mod:`experiment` -- Experiment
================================================================================

.. module:: experiment
   :synopsis: Experiment

.. inheritance-diagram:: experiment

"""

# Script information for the file.
__author__ = "Philippe T. Pinard"
__email__ = "philippe.pinard@gmail.com"
__version__ = "0.1"
__copyright__ = "Copyright (c) 2011 Philippe T. Pinard"
__license__ = "GPL v3"

# Standard library modules.

# Third party modules.

# Local modules.

# Globals and constants variables.
LINE_KA = 0
LINE_KB = 1
LINE_LA = 2
LINE_LB = 3
LINE_MA = 4
LINE_MB = 5

class Experiment:
    def __init__(self, z, line, energy_eV, kratio=0.0, standard='', analyzed=True):
        """
        Creates a new experiment. 
        An experiment indicates the measurement conditions and/or results.
        
        :arg z: atomic number
        :arg line: x-ray characteristic line (use constant)
        :type line: :class:`int`
        :arg energy_eV: beam energy (in eV)
        :arg kratio: measured k-ratio
        :arg analyzed: whether to use this experiment in the calculations
        """
        self._z = z
        self._line = line
        self._energy_eV = energy_eV
        self._kratio = kratio
        self._standard = standard
        self._analyzed = analyzed

    def is_analyzed(self):
        """
        Whether to use this experiment in the calculations.
        """
        return self._analyzed

    @property
    def z(self):
        """
        Returns the atomic number.
        """
        return self._z

    @property
    def line(self):
        """
        Returns the x-ray characteristic line.
        """
        return self._line

    @property
    def energy_eV(self):
        """
        Return the beam energy (in eV).
        """
        return self._energy_eV

    @property
    def kratio(self):
        """
        Returns the measured k-ratio or ``0.0`` if no k-ratio was measured.
        """
        return self._kratio

    @property
    def standard(self):
        return self._standard
