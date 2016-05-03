"""
Definition of an experiment, structure to setup experimental parameters and
measurements in STRATAGem.
"""

__all__ = ['LINE_KA', 'LINE_KB', 'LINE_LA', 'LINE_LB', 'LINE_MA', 'LINE_MB',
           'Experiment']

# Standard library modules.

# Third party modules.

# Local modules.
from stratagemtools.element_properties import symbol

# Globals and constants variables.

LINE_KA = 0
"""X-ray line :math:`\\text{K}\\alpha`"""

LINE_KB = 1
"""X-ray line :math:`\\text{K}\\beta`"""

LINE_LA = 2
"""X-ray line :math:`\\text{L}\\alpha`"""

LINE_LB = 3
"""X-ray line :math:`\\text{L}\\beta`"""

LINE_MA = 4
"""X-ray line :math:`\\text{M}\\alpha`"""

LINE_MB = 5
"""X-ray line :math:`\\text{M}\\beta`"""

class Experiment:
    """
    Object to store experimental parameters and measurements.
    Once created an experiment object is immutable.
    """

    def __init__(self, z, line, energy_eV, kratio=0.0, standard='', analyzed=True):
        """
        :arg z: atomic number
        :type z: :class:`int`
        
        :arg line: X-ray characteristic line, either
            
            * :data:`LINE_KA`
            * :data:`LINE_KB`
            * :data:`LINE_LA`
            * :data:`LINE_LB`
            * :data:`LINE_MA`
            * :data:`LINE_MB`
            
            Note that no other X-ray lines are supported.
        :type line: :class:`int`
        
        :arg energy_eV: beam energy (in eV)
        :type energy_eV: :class:`float`
        
        :arg kratio: measured k-ratio (optional)
        :type kratio: :class:`float`
        
        :arg standard: three options
        
            * empty string for pure standard (e.g. ``Fe``)
            * standard name which correspond to the filename of the standard
              saved in the standard directory 
              (see :attr:`Stratagem.standard_directory`)
            * a :class:`Sample`
            
        :arg analyzed: whether to use this experiment in the calculations
        :type analyzed: :class:`bool`
        """
        self._z = z
        self._line = line
        self._energy_eV = energy_eV
        self._kratio = kratio
        self._standard = standard
        self._analyzed = analyzed

    def __repr__(self):
        line = {LINE_KA: 'Ka', LINE_KB: 'Kb',
                LINE_LA: 'La', LINE_LB: 'Lb',
                LINE_MA: 'Ma', LINE_MB: 'Mb'}[self.line]
        energy_keV = self.energy_eV / 1e3
        standard = self.standard if self.standard else 'pure'
        extra = 'analyzed' if self.is_analyzed() else 'not analyzed'
        return '<Experiment(%s %s, %s kV, kratio=%.4f, standard=%s, %s)>' % \
            (symbol(self.z), line, energy_keV, self.kratio, standard, extra)

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
        """
        REturns the standard.
        """
        return self._standard
