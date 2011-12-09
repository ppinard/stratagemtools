#!/usr/bin/env python
"""
================================================================================
:mod:`stratagem` -- Stratagem runner
================================================================================

.. module:: stratagem
   :synopsis: Stratagem runner

.. inheritance-diagram:: stratagemtools.stratagem

"""

# Script information for the file.
__author__ = "Philippe T. Pinard"
__email__ = "philippe.pinard@gmail.com"
__version__ = "0.1"
__copyright__ = "Copyright (c) 2011 Philippe T. Pinard"
__license__ = "GPL v3"

# Standard library modules.
import ctypes as c
import logging as l
from ConfigParser import SafeConfigParser

# Third party modules.

# Local modules.

# Globals and constants variables.
_SECTION_STRATAGEM = "Stratagem"
_OPTION_DLLPATH = "dllPath"

PRZMODE_XPP = 0
PRZMODE_PAP = 1
PRZMODE_GAU = 2

FLUORESCENCE_NONE = 0
FLUORESCENCE_LINE = 1
FLUORESCENCE_LINE_CONT = 2

class StratagemError(Exception):
    pass

class Stratagem:
    def __init__(self, configfile=None, dll_path=None):
        """
        Initializes the connection to the Stratagem DLL.
        One of the following argument must be specified:
        
        :arg configfile: file-object to the configuration file
        :arg dll_path: complete path to the location of the ``stratadll.dll``
        """
        if configfile:
            config = SafeConfigParser()
            config.readfp(configfile)

            if config.has_section(_SECTION_STRATAGEM):
                if config.has_option(_SECTION_STRATAGEM, _OPTION_DLLPATH):
                    dll_path = config.get(_SECTION_STRATAGEM, _OPTION_DLLPATH)

        l.debug("dll=%s", dll_path)
        self._lib = c.WinDLL(dll_path)

        self._key = c.create_string_buffer('DemoSample')
        self._stobjectnew(self._key)
        self._stenableerrordisplay(False)

        self.reset()

    def _stobjectnew(self, key):
        bNormal_ = c.c_bool(True)
        iniFlags_ = c.c_int(0)

        l.debug("StObjectNew(key, %r, %i)", True, 0)
        if not self._lib.StObjectNew(self._key, bNormal_, iniFlags_):
            raise StratagemError, "Cannot create object"

    def _stenableerrordisplay(self, enable):
        enable_ = c.c_bool(enable)

        l.debug("StEnableErrorDisplay(%r)", enable)
        self._lib.StEnableErrorDisplay(enable_)

    def close(self):
        """
        Closes the connection to the Statagem DLL.
        """
        self._lib.StObjectDelete(self._key)
        del self._lib

    def reset(self):
        """
        Resets all parameters to the defaults, remove all layers and experiments.
        """
        self._lib.StObjectReset(self._key)
        self._layers = {} # layer: index
        self._substrate = None
        self._experiments = {} # analyzed experiments

    def add_layer(self, layer, substrate=False):
        """
        Adds a layer from top to bottom. 
        The last layer added is considered as the substrate.
        
        :arg layer: layer
        :type layer: :class:`Layer`
        
        :return: index of the layer
        """
        l.debug("StSdAddLayer(key)")
        iLayer_ = self._lib.StSdGetNbLayers(self._key)

        l.debug("StSdAddLayer(key, %i)", iLayer_)
        if not self._lib.StSdAddLayer(self._key, iLayer_):
            raise StratagemError, "Cannot add layer"

        for i, element in enumerate(layer.iter_elements()):
            iElt_ = c.c_int(i)
            l.debug("StSdAddElt(key, %i, %i)", iLayer_, i)
            if not self._lib.StSdAddElt(self._key, iLayer_, iElt_):
                raise StratagemError, "Cannot add element"

            z, conc = element
            nra_ = c.c_int(z)
            l.debug("StSdSetNrAtom(key, %i, %i, %i)", iLayer_, i, z)
            if not self._lib.StSdSetNrAtom(self._key, iLayer_, iElt_, nra_):
                raise StratagemError, "Cannot set atomic number"

            if conc >= 0:
                flag = 0

                wf_ = c.c_double(conc)
                l.debug("StSdSetConc(key, %i, %i, %f)", iLayer_, i, conc)
                if not self._lib.StSdSetConc(self._key, iLayer_, iElt_, wf_):
                    raise StratagemError, "Cannot set concentration"
            else:
                flag = 1

            l.debug("StSdSetConcFlag(key, %i, %i, %i)", iLayer_, i, flag)
            if not self._lib.StSdSetConcFlag(self._key, iLayer_, iElt_, c.c_int(flag)):
                raise StratagemError, "Cannot set concentration flag"

        if not substrate:
            thickKnown = layer.is_thickness_known()
            thickKnown_ = c.c_bool(thickKnown)
            thickness = layer.thickness * 10 # Angstroms
            thickness_ = c.c_double(thickness)
            mass_thickness = layer.mass_thickness / 1e6 # g/cm2
            mass_thickness_ = c.c_double(mass_thickness)
            density = layer.density
            density_ = c.c_double(density)
            if density <= 0.0: # calculate density
                l.debug('StSdDefaultDensity(key, %i)', iLayer_)
                self._lib.StSdDefaultDensity(self._key, iLayer_, c.byref(density_))

            l.debug("StSdSetThick(key, %i, %r, %d, %d, %d)", iLayer_, thickKnown,
                    mass_thickness, thickness, density)
            if not self._lib.StSdSetThick(self._key, iLayer_, thickKnown_,
                                          mass_thickness_, thickness_, density_):
                raise StratagemError, "Cannot set thickness"

            self._layers.setdefault(layer, int(iLayer_))

    def add_substrate(self, layer):
        """
        Adds a layer as the substrate.
        """
        if layer is None:
            raise ValueError, "The layer cannot be None"
        if self._substrate is not None:
            raise ValueError, "A substrate was already defined."

        self.add_layer(layer, substrate=True)
        self._substrate = layer

    def add_experiment(self, experiment):
        """
        Add an experiment, measurements of k-ratio at different energies.
        
        :arg experiment: experiment
        :type experiment: :class:`Experiment`
        """
        nra_ = c.c_int(experiment.z)
        klm_ = c.c_int(experiment.line)
        hv_ = c.c_double(experiment.hv)
        iElt_ = c.c_int()
        iLine_ = c.c_int()
        iExpK_ = c.c_int()
        l.debug('StEdAddNrAtomLineHV(key, %i, %i)', experiment.z, experiment.line)
        if not self._lib.StEdAddNrAtomLineHV(self._key, nra_, klm_, hv_,
                                             c.byref(iElt_), c.byref(iLine_), c.byref(iExpK_)):
            raise StratagemError, "Cannot add atomic number and line"

        analyzed = experiment.is_analyzed()
        analyzed_ = c.c_bool(analyzed)
        l.debug("StEdSetAnalyzedFlag(key, %i, %r)", iElt_.value, analyzed)
        if not self._lib.StEdSetAnalyzedFlag(self._key, iElt_, analyzed_):
            raise StratagemError, "Cannot add experiment analyzed flag"

        hv_ = c.c_double(experiment.hv)
        kratio_ = c.c_double(experiment.kratio)
        l.debug("StEdSetExpK(key, %i, %i, %i, %f, %f, %f, 0.0, 2)",
                iElt_.value, iLine_.value, iExpK_.value, experiment.hv,
                experiment.hv, experiment.kratio)
        if not self._lib.StEdSetExpK(self._key, iElt_, iLine_, iExpK_,
                                     hv_, hv_, kratio_, c.c_double(0.0),
                                     c.c_int(2)):
            raise StratagemError, "Cannot set experiment k-ratio"

        if analyzed:
            indexes = (iElt_.value, iLine_.value, iExpK_.value)
            self._experiments.setdefault(experiment, indexes)

    def set_geometry(self, toa, tilt, azimuth):
        """
        Sets the geometry.
        
        :arg toa: take off angle (in radians)
        :arg tilt: tilt angle (in radians)
        :arg azimuth: azimuthal angle (in radians)
        """
        toa_ = c.c_double(toa)
        tilt_ = c.c_double(tilt)
        azimuth_ = c.c_double(azimuth)
        l.debug('StSetGeomParams(key, %f, %f, %f)', toa, tilt, azimuth)
        if not self._lib.StSetGeomParams(self._key, toa_, tilt_, azimuth_):
            raise StratagemError, "Cannot set geometry parameters"

    def set_prz_mode(self, mode):
        """
        Sets the type of model to use for the phi-rho-z.
        """
        mode_ = c.c_int(mode)
        l.debug('StSetPrzMode(%i)', mode)
        self._lib.StSetPrzMode(mode_)

    def set_fluorescence(self, flag):
        """
        Sets whether to consider characteristic fluorescence, characteristic
        and continuum fluorescence or no fluoresence.
        """
        flag_ = c.c_int(flag)
        l.debug('StSetFluorFlg(%i)', flag)
        self._lib.StSetFluorFlg(flag_)

    def compute_kratio_vs_thickness(self, layer, thickness_low, thickness_high, step):
        """
        Computes the variation of the k-ratio as a function of the mass 
        thickness for a layer.
        
        :arg layer: layer (must have been previously added)
        :arg thickness_low: lower limit of the thickness (in nm)
        :arg thickness_high: upper limit of the thickness (in nm)
        :arg step: number of steps
        
        :return: :class:`list` of thicknesses, :class:`dict` of experiment-kratios
        """
        l.debug('StSetKvsThicknessUnit(2)')
        self._lib.StSetKvsThicknessUnit(2) # unit in nm

        if layer not in self._layers:
            raise ValueError, "Unknown layer"
        iLayer = self._layers[layer]
        iLayer_ = c.c_int(iLayer)

        step_ = c.c_int(step)
        l.debug('StSetNbComputedHV(%i)', step)
        self._lib.StSetNbComputedHV(step_)

        # Compute
        low_ = c.c_double(thickness_low)
        high_ = c.c_double(thickness_high)
        l.debug('StComputeKvsThickness(key, %i, %f, %f)',
                iLayer, thickness_low, thickness_high)
        if not self._lib.StComputeKvsThickness(self._key, iLayer_, low_, high_):
            raise StratagemError, "Cannot compute k-ratio vs thickness"

        # Fetch results
        thicknesses = []
        kratios = {}

        thick_ = c.c_double()
        k_ = c.c_double()
        for i in range(step + 1):
            i_ = c.c_int(i)

            if not self._lib.StGetKvsT_Thick(self._key, i_, c.byref(thick_)):
                raise StratagemError, "Cannot get thickness"
            thicknesses.append(thick_.value)

            for experiment, indexes in self._experiments.iteritems():
                iElt_ = c.c_int(indexes[0])
                iLine_ = c.c_int(indexes[1])
                iHv_ = c.c_int(indexes[2])

                if not self._lib.StGetKvsT_K(self._key, i_, iElt_, iLine_,
                                             iHv_, c.byref(k_)):
                    raise StratagemError, "Cannot get k-ratio"
                kratios.setdefault(experiment, []).append(k_.value)

        return thicknesses, kratios

    def compute_kratios(self):
        """
        Computes the kratios of the different experiments.
        
        :return: :class:`dict` of experiment-kratios
        """
        for i, layer in enumerate(self._layers.keys()):
            if not layer.is_thickness_known():
                raise ValueError, "Thickness of layer %i is unknown" % i

        # Compute
        layer = self._layers.keys()[0]
        thickness_low = layer.thickness
        thickness_high = layer.thickness * 10
        step = 1

        _thicknesses, kratios = \
            self.compute_kratio_vs_thickness(layer, thickness_low, thickness_high, step)

        # Reorganize results
        output = {}
        for experiment, kratio in kratios.iteritems():
            output.setdefault(experiment, kratio[0])

        return output

    def compute_thicknesses(self, iteration_max=50):
        """
        Computes the thicknesses of each layer.
        
        :return: :class:`dict` of layer-thicknesses (in nm)
        """
        iteration_max_ = c.c_int(iteration_max)
        l.debug('StSetMaxNbIter(%i)', iteration_max)
        self._lib.StSetMaxNbIter(iteration_max_)

        # Compute
        l.debug('StComputeIterpStart(key)')
        if not self._lib.StComputeIterpStart(self._key):
            raise StratagemError, "Cannot start iteration"

        continue_ = c.c_bool(True)
        iteration = 0

        while(True):
            iteration += 1

            l.debug('StComputeIterpNext(key, %r)' % continue_.value)
            if not self._lib.StComputeIterpNext(self._key, c.byref(continue_)):
                break

            if not continue_.value:
                break

        # Fetch results
        thicknesses = {}

        thickKnown = c.c_bool()
        massThickness = c.c_double()
        thickness = c.c_double()
        density = c.c_double()

        for layer, iLayer in self._layers.iteritems():
            iLayer_ = c.c_int(iLayer)
            l.debug("StSdGetThick(key, %i)", iLayer)

            if not self._lib.StSdGetThick(self._key, iLayer_, c.byref(thickKnown),
                                          c.byref(massThickness), c.byref(thickness),
                                          c.byref(density)):
                raise StratagemError, "Cannot get thickness"

            thicknesses.setdefault(layer, thickness.value / 10)

        return thicknesses

