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

DENSITIES = [
    0.0899, 0.1787, 0.5300, 1.8500, 2.3400, 2.6200, 1.2510, 1.4290,
    1.6960, 0.9010, 0.9700, 1.7400, 2.7000, 2.3300, 1.8200, 2.0700,
    3.1700, 1.7840, 0.8600, 1.5500, 3.0000, 4.5000, 5.8000, 7.1900,
    7.4300, 7.8600, 8.9000, 8.9000, 8.9600, 7.1400, 5.9100, 5.3200,
    5.7200, 4.8000, 3.1200, 3.7400, 1.5300, 2.6000, 4.5000, 6.4900,
    8.5500, 10.200, 11.500, 12.200, 12.400, 12.000, 10.500, 8.6500,
    7.3100, 7.3000, 6.6800, 6.2400, 4.9200, 5.8900, 1.8700, 3.5000,
    6.7000, 6.7800, 6.7700, 7.0000, 6.4750, 7.5400, 5.2600, 7.8900,
    8.2700, 8.5400, 8.8000, 9.0500, 9.3300, 6.9800, 9.8400, 13.100,
    16.600, 19.300, 21.000, 22.400, 22.500, 21.400, 19.300, 13.530,
    11.850, 11.400, 9.8000, 9.4000, 1.0000, 9.9100, 1.0000, 5.0000,
    10.070, 11.700, 15.400, 18.900, 20.400, 19.800, 13.600, 13.511
]

class Layer:
    def __init__(self, elements={}, thickness=0.0, mass_thickness=0.0, density=0.0):
        """
        Creates a new layer or substrate.
        
        It is recommended to only specify either the thickness or the mass
        thickness. If both are specified, their value must be correct.
        
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

        if density <= 0.0:
            density = 0.0
            for z, wf in elements.iteritems():
                density += DENSITIES[z - 1] * wf
        self._density = density

        if thickness <= 0.0 and mass_thickness <= 0.0:
            thickness = -1.0
            mass_thickness = -1.0
        elif thickness > 0.0 and mass_thickness <= 0.0:
            mass_thickness = thickness * density / 10
        elif thickness <= 0.0 and mass_thickness > 0.0:
            thickness = mass_thickness / density * 10
        else:
            assert thickness == mass_thickness / density * 10

        self._thickness = thickness
        self._mass_thickness = mass_thickness

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
        return self._thickness

    @property
    def mass_thickness(self):
        """
        Returns mass thickness in ug/cm2.
        A negative mass thickness can be set, it indicates that the thickness is
        unknown.
        """
        return self._mass_thickness

    @property
    def density(self):
        """
        Returns density in g/cm3.
        """
        return self._density

