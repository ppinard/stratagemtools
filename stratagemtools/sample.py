"""
Definition of a sample, corresponding to a multilayer.
A multilayer a series of horizontal layers deposited on a substrate.
A multilayer may also have no layer and simply be considered as a substrate.
The substrate composition, and the composition, thickness and mass thickness 
of each layer are defined by the :class:`Sample` class.
"""

__all__ = ['CONC_UNKNOWN', 'CONC_DIFF', 'composition_from_formula',
           'Layer', 'Sample']

# Standard library modules.
import string
from collections import defaultdict

# Third party modules.
from pyparsing import Word, Group, Optional, OneOrMore

# Local modules.
import stratagemtools.element_properties as ep

# Globals and constants variables.
CONC_UNKNOWN = None
"""Flag when the composition of an element is unknown."""

CONC_DIFF = '?'
"""Flag when the composition of an element should be calculated by difference."""

_symbol = Word(string.ascii_uppercase, string.ascii_lowercase)
_digit = Word(string.digits + ".")
_elementRef = Group(_symbol + Optional(_digit, default="1"))
CHEMICAL_FORMULA_PARSER = OneOrMore(_elementRef)

def composition_from_formula(formula):
    """
    Calculates the composition (expressed in weight fractions) of a chemical 
    formula.
    
    Example::
    
        >>> composition_from_formula('Al2O3')
        ... {8: 0.4707492883573059, 13: 0.5292507116426941}
    
    :arg formula: a valid chemical formula
    :type formula: :class:`str`
    
    :return: composition (expressed in weight fractions). The keys of the 
        :class:`dict` are atomic numbers and the values, weight fractions.
    :rtype: :class:`dict`
    """

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

class Layer:
    """
    Object to store a layer definition.
    Once created a layer object is immutable.
    
    .. note:: Should not be used directly, but via the :class:`Sample`'s methods
       :meth:`add_layer <.Sample.add_layer>` and 
       :meth:`get_layer <.Sample.get_layer>`.
    """

    def __init__(self, composition, thickness_m=None,
                 mass_thickness_kg_m2=None, density_kg_m3=None):
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
        """
        Returns whether the composition is known, i.e. contains no 
        :data:`CONC_UNKNOWN` or :data:`CONC_DIFF` flag.
        
        :rtype: :class:`bool`
        """
        return self._is_composition_known

    def is_thickness_known(self):
        """
        Returns whether the thickness is known.
        
        :rtype: :class:`bool`
        """
        return self._is_thickness_known

    def is_density_known(self):
        """
        Returns whether the density is known.
        
        :rtype: :class:`bool`
        """
        return self._is_density_known

    @property
    def composition(self):
        """
        Returns a copy of the composition :class:`dict`. 
        The composition, and any other parameters, cannot be modified.
        
        :return: composition (expressed in weight fractions). The keys of the 
            :class:`dict` are atomic numbers and the values, weight fractions.
        :rtype: :class:`dict`
        """
        return self._composition.copy()

    @property
    def thickness_m(self):
        """
        Returns thickness in meters.
        It can be ``None`` if not defined.
        
        :rtype: :class:`float`
        """
        return self._thickness_m

    @property
    def mass_thickness_kg_m2(self):
        """
        Returns mass thickness in kg/m2.
        It can be ``None`` if not defined.
        
        :rtype: :class:`float`
        """
        return self._mass_thickness_kg_m2

    @property
    def density_kg_m3(self):
        """
        Returns density in kg/m3.
        It can be ``None`` if not defined.
        
        :rtype: :class:`float`
        """
        return self._density_kg_m3

class Sample:
    """
    Object to store a multilayer sample definition.
    """

    def __init__(self, composition, density_kg_m3=None):
        """
        :arg composition: composition of the substrate as :class:`dict` where
            the keys are atomic numbers and the values, weight fractions. 
            If the weight fraction is not known, set it to :data:`CONC_UNKNOWN`,
            if the weight fraction should be calculated by difference, set it
            to :data:`CONC_DIFF`.
        :type composition: :class:`dict`
        
        :arg density_kg_m3: mass density in kilograms per cubic meter (optional).
            If the composition is known, the density will be automatically
            calculated based on:
            
            .. math:: 
                
                \\frac{1}{\\rho} = \\sum \\frac{w_i}{\\rho_i}
            
            where :math:`w_i` and :math:`\\rho_i` are respectively the weight 
            fraction and mass density of element :math:`i`.
        :type density_kg_m3: :class:`float`
        """
        self._substrate = Layer(composition, density_kg_m3=density_kg_m3)
        self._layers = []

    def __repr__(self):
        comp_str = ', '.join('%s: %.4f' % (ep.symbol(z), wf) \
                             for z, wf in self.composition.items())
        return '<%s(%s, %i layers)>' % (self.__class__.__name__, comp_str,
                                        len(self._layers))

    def add_layer(self, composition, thickness_m=None,
                  mass_thickness_kg_m2=None, density_kg_m3=None):
        """
        Adds a layer below the previously added layer, or if no layer was added
        on top of the substrate.
        
        :arg composition: composition of the layer as :class:`dict` where
            the keys are atomic numbers and the values, weight fractions. 
            If the weight fraction is not known, set it to :data:`CONC_UNKNOWN`,
            if the weight fraction should be calculated by difference, set it
            to :data:`CONC_DIFF`.
        :type composition: :class:`dict`
        
        :arg thickness_m: thickness of the layer in meters (optional). If the
            *mass_thickness_kg_m2* and *density_kg_m3* are known, the thickness
            will be automatically calculated.
        :type thickness_m: :class:`float`
        
        :arg mass_thickness_kg_m2: mass thickness of the layer in kilograms per
            square meter (optional). The mass thickness is defined as the
            thickness times the density. If *thickness_m* and *density_kg_m3*
            are known the mass thickness will be automatically calculated.
        :type mass_thickness_kg_m2: :class:`float`
        
        :arg density_kg_m3: mass density in kilograms per cubic meter (optional).
            If the composition is known, the density will be automatically
            calculated based on:
            
            .. math:: 
                
                \\frac{1}{\\rho} = \\sum \\frac{w_i}{\\rho_i}
            
            where :math:`w_i` and :math:`\\rho_i` are respectively the weight 
            fraction and mass density of element :math:`i`.
        :type density_kg_m3: :class:`float`
        
        :return: a layer
        :rtype: :class:`Layer`
        """
        layer = Layer(composition, thickness_m, mass_thickness_kg_m2, density_kg_m3)
        self._layers.append(layer)
        return layer

    def pop_layer(self, index):
        """
        Removes the layer at *index*.
        
        :arg index: index of the layer to be removed
        :type index: :class:`int`
        """
        self._layers.pop(index)

    def get_layer(self, index):
        """
        Returns the layer at *index*.
        
        :arg index: index of the layer
        :type index: :class:`int`
        
        :return: a layer
        :rtype: :class:`Layer`
        """
        return self._layers[index]

    @property
    def composition(self):
        """
        Returns the a copy of the composition of the substrate.
        The composition, and any other parameters, cannot be modified.
        
        :return: composition (expressed in weight fractions). The keys of the 
            :class:`dict` are atomic numbers and the values, weight fractions.
        :rtype: :class:`dict`
        """
        return self._substrate.composition

    @property
    def substrate(self):
        """
        Returns the "layer" of the substrate, a :class:`Layer` object 
        corresponding to the composition and density of the substrate.
        
        :rtype: :class:`Layer`
        """
        return self._substrate

    @property
    def layers(self):
        """
        Returns a copy of layers of this sample.
        It cannot be modified.
        
        :rtype: :class:`tuple`
        """
        return tuple(self._layers)
