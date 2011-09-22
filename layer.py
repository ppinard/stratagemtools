#!/usr/bin/env python
"""
================================================================================
:mod:`layer` -- A layer
================================================================================

.. module:: layer
   :synopsis: A layer

.. inheritance-diagram:: layer

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

class Layer:
    def __init__(self, elements={}, thickness=0.0, mass_thickness=0.0, density=0.0):
        """
        Creates a new layer or substrate.
        
        :arg elements: :class:`dict` of atomic number / weight fraction. 
            If the concentration is negative, it is unknown.
        :arg thickness: thickness (in nm).
            If the thickness is negative, it is unknown.
        :arg mass_thickness: mass thickness (in ug/cm2)
            If the mass thickness is negative, it is unknown.
        :arg density: density (in g/cm3)
            If the density is less or equal to zero, the weighted density based
            on the concentration is used.
        """
        self._elements = elements
        self._thickness = thickness
        self._mass_thickness = mass_thickness
        self._density = density

    def iter_elements(self):
        """
        Iterator over the atomic number and weight fraction.
        
        Example::
        
          for z, wf in layer.iter_elements():
            print z, wf
        """
        return self._elements.iteritems()

    def is_thickness_known(self):
        return self._thickness >= 0 and self._mass_thickness >= 0

    @property
    def thickness(self):
        """
        Returns thickness in nm. 
        A negative thickness can be set, it indicates that the thickness is
        unknown.
        """
        if self._thickness > 0.0:
            return self._thickness
        elif self._mass_thickness > 0.0 and self._density > 0.0:
            return self._mass_thickness / self._density * 10
        else:
            return 0.0

    @property
    def mass_thickness(self):
        """
        Returns mass thickness in ug/cm2.
        A negative mass thickness can be set, it indicates that the thickness is
        unknown.
        """
        if self._mass_thickness > 0.0:
            return self._mass_thickness
        elif self._thickness > 0.0 and self._density > 0.0:
            return self._thickness * self._density / 10
        else:
            return 0.0

    @property
    def density(self):
        """
        Returns density in g/cm3.
        """
        if self._density > 0.0:
            return self._density
        elif self._thickness > 0.0 and self._mass_thickness > 0.0:
            return 10 * self._mass_thickness / self._thickness
        else:
            return 0.0

