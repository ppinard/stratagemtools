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
import stratagemtools.element_properties as ep

# Globals and constants variables.

class Layer:

    def __init__(self, composition={}, thickness_m=None,
                 mass_thickness_kg_m2=None, density_kg_m3=None):
        """
        Creates a new layer or substrate.
        
        It is recommended to only specify either the thickness or the mass
        thickness. If both are specified, their value must be correct.
        
        :arg composition: :class:`dict` of atomic number / weight fraction. 
            If the concentration is ``None``, it is unknown.
        :arg thickness: thickness (in nm).
            If the thickness is ``None``, it is unknown.
        :arg mass_thickness: mass thickness (in ug/cm2)
            If the mass thickness is ``None``, it is unknown.
        :arg density: density (in g/cm3)
            If the density is less or equal to zero, the weighted density based
            on the concentration is used.
        """
        is_composition_known = True
        for z, wf in composition.items():
            if wf is None:
                is_composition_known = False
                break
        self._composition = composition.copy()

        if density_kg_m3 is None and is_composition_known:
            density_kg_m3 = 0.0
            for z, wf in composition.items():
                density_kg_m3 += wf / ep.mass_density_kg_m3(z)
            density_kg_m3 = 1.0 / density_kg_m3
        self._density_kg_m3 = density_kg_m3
        self._is_density_known = density_kg_m3 is not None

        if thickness_m is not None and \
                mass_thickness_kg_m2 is None and \
                density_kg_m3 is not None:
            mass_thickness_kg_m2 = thickness_m * density_kg_m3
        elif thickness_m is None and \
                mass_thickness_kg_m2 is not None and \
                density_kg_m3 is not None:
            thickness_m = mass_thickness_kg_m2 / density_kg_m3

        self._thickness_m = thickness_m
        self._mass_thickness_kg_m2 = mass_thickness_kg_m2
        self._is_thickness_known = thickness_m is not None and \
            mass_thickness_kg_m2 is not None

    def __repr__(self):
        comp_str = ', '.join('%s: %s' % (ep.symbol(z), wf) \
                             for z, wf in self.composition.items())
        thickness_str = '%s nm' % (self.thickness_m * 1e9,) \
            if self.is_thickness_known() else 'unknown'
        return '<Layer(composition={%s}, thickness=%s)>' % (comp_str, thickness_str)

    def is_thickness_known(self):
        return self._is_thickness_known

    def is_density_known(self):
        return self._is_density_known

    @property
    def composition(self):
        return self._composition.copy()

    @property
    def thickness_m(self):
        """
        Returns thickness in meters. 
        """
        return self._thickness_m

    @property
    def mass_thickness_kg_m2(self):
        """
        Returns mass thickness in kg/m2.
        """
        return self._mass_thickness_kg_m2

    @property
    def density_kg_m3(self):
        """
        Returns density in kg/m3.
        """
        return self._density_kg_m3

