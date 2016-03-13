#!/usr/bin/env python

# Script information for the file.
__author__ = "Philippe T. Pinard"
__email__ = "philippe.pinard@gmail.com"
__copyright__ = "Copyright (c) 2015 Philippe T. Pinard"
__license__ = "GPL v3"

# Standard library modules.
import string
from collections import defaultdict

# Third party modules.
from pyparsing import Word, Group, Optional, OneOrMore

# Local modules.
import stratagemtools.element_properties as ep

# Globals and constants variables.
CONC_UNKNOWN = None
CONC_DIFF = '?'

_symbol = Word(string.ascii_uppercase, string.ascii_lowercase)
_digit = Word(string.digits + ".")
_elementRef = Group(_symbol + Optional(_digit, default="1"))
CHEMICAL_FORMULA_PARSER = OneOrMore(_elementRef)

def composition_from_formula(formula):
    # Parse chemical formula
    formulaData = CHEMICAL_FORMULA_PARSER.parseString(formula)

    zs = []
    atomicfractions = []
    for symbol, atomicfraction in formulaData:
        zs.append(ep.atomic_number(symbol=symbol))
        atomicfractions.append(float(atomicfraction))

    # Calculate total atomic mass
    totalatomicmass = 0.0
    for z, atomicfraction in zip(zs, atomicfractions):
        atomicmass = ep.atomic_mass_kg_mol(z)
        totalatomicmass += atomicfraction * atomicmass

    # Create composition
    composition = defaultdict(float)

    for z, atomicfraction in zip(zs, atomicfractions):
        atomicmass = ep.atomic_mass_kg_mol(z)
        weightfraction = atomicfraction * atomicmass / totalatomicmass
        composition[z] += weightfraction

    return composition

class _Layer:

    def __init__(self, composition, thickness_m=None,
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
        self._is_composition_known = \
            all(isinstance(wf, float) for wf in composition.values())
        self._composition = composition.copy()

        if density_kg_m3 is None and self._is_composition_known:
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
        comp_str = ', '.join('%s: %.4f' % (ep.symbol(z), wf) \
                             for z, wf in self.composition.items())
        thickness_str = '%s nm' % (self.thickness_m * 1e9,) \
            if self.is_thickness_known() else 'unknown'
        return '<Layer(composition={%s}, thickness=%s)>' % (comp_str, thickness_str)

    def is_composition_known(self):
        return self._is_composition_known

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

class Sample:

    def __init__(self, composition):
        """
        Creates a sample.
        
        :arg composition: composition of the substrate as :class:`dict` where
            the keys are atomic numbers and the values, weight fractions. 
            If the weight fraction is not known, set is to ``None``.
        """
        self._substrate = _Layer(composition)
        self._layers = []

    def __repr__(self):
        comp_str = ', '.join('%s: %.4f' % (ep.symbol(z), wf) \
                             for z, wf in self.composition.items())
        return '<%s(%s, %i layers)>' % (self.__class__.__name__, comp_str,
                                        len(self._layers))

    def add_layer(self, composition, thickness_m=None,
                  mass_thickness_kg_m2=None, density_kg_m3=None):
        layer = _Layer(composition, thickness_m, mass_thickness_kg_m2, density_kg_m3)
        self._layers.append(layer)
        return layer

    def pop_layer(self, index):
        self._layers.pop(index)

    def get_layer(self, index):
        return self._layers[index]

    @property
    def composition(self):
        return self._substrate.composition

    @property
    def substrate(self):
        return self._substrate

    @property
    def layers(self):
        return tuple(self._layers)
