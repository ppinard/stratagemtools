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
from configparser import SafeConfigParser
from operator import attrgetter

# Third party modules.

# Local modules.

# Globals and constants variables.
from stratagemtools.layer import Layer, DENSITIES, ATOMIC_MASSES

_SECTION_STRATAGEM = "Stratagem"
_OPTION_DLLPATH = "dllPath"

PRZMODE_XPP = 0
PRZMODE_PAP = 1
PRZMODE_GAU = 2

FLUORESCENCE_NONE = 0
FLUORESCENCE_LINE = 1
FLUORESCENCE_LINE_CONT = 2

CONCENTRATION_FLAG_KNOWN = 0
CONCENTRATION_FLAG_UNKNOWN = 1
CONCENTRATION_FLAG_STOICHIOMETRIC = 2
CONCENTRATION_FLAG_TRACE = 3
CONCENTRATION_FLAG_DIFFERENCE = 4

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

        self._key = c.create_string_buffer(b'DemoSample')
        self._stobjectnew(self._key)
        self._stenableerrordisplay(False)

        self.reset()

    def _stobjectnew(self, key):
        bNormal_ = c.c_bool(True)
        iniFlags_ = c.c_int(0)

        l.debug("StObjectNew(key, %r, %i)", True, 0)
        if not self._lib.StObjectNew(self._key, bNormal_, iniFlags_):
            raise StratagemError("Cannot create object")

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
            raise StratagemError("Cannot add layer")

        for i, value in enumerate(layer.composition.items()):
            iElt_ = c.c_int(i)
            l.debug("StSdAddElt(key, %i, %i)", iLayer_, i)
            if not self._lib.StSdAddElt(self._key, iLayer_, iElt_):
                raise StratagemError("Cannot add element")

            z, wf = value
            nra_ = c.c_int(z)
            l.debug("StSdSetNrAtom(key, %i, %i, %i)", iLayer_, i, z)
            if not self._lib.StSdSetNrAtom(self._key, iLayer_, iElt_, nra_):
                raise StratagemError("Cannot set atomic number")

            if wf is not None:
                flag = CONCENTRATION_FLAG_KNOWN

                wf_ = c.c_double(wf)
                l.debug("StSdSetConc(key, %i, %i, %f)", iLayer_, i, wf)
                if not self._lib.StSdSetConc(self._key, iLayer_, iElt_, wf_):
                    raise StratagemError("Cannot set concentration")
            else:
                flag = CONCENTRATION_FLAG_UNKNOWN

            l.debug("StSdSetConcFlag(key, %i, %i, %i)", iLayer_, i, flag)
            if not self._lib.StSdSetConcFlag(self._key, iLayer_, iElt_, c.c_int(flag)):
                raise StratagemError("Cannot set concentration flag")

        if not substrate:
            thick_known = layer.is_thickness_known()
            thick_known_ = c.c_bool(thick_known)

            if layer.is_density_known():
                density = layer.density_kg_m3 / 1e3 # g/cm3
            else:
                density = 10.0
            density_ = c.c_double(density)

            if thick_known:
                thickness = layer.thickness_m * 1e10  # Angstroms
                mass_thickness = layer.mass_thickness_kg_m2 * 0.1 # g/cm2
            else:
                thickness = 0.0
                mass_thickness = 0.0
            thickness_ = c.c_double(thickness)
            mass_thickness_ = c.c_double(mass_thickness)

            l.debug("StSdSetThick(key, %i, %r, %d, %d, %d)", iLayer_,
                    thick_known, mass_thickness, thickness, density)
            if not self._lib.StSdSetThick(self._key, iLayer_, thick_known_,
                                          mass_thickness_, thickness_, density_):
                raise StratagemError("Cannot set thickness")

            self._layers.setdefault(layer, int(iLayer_))

    def add_substrate(self, layer):
        """
        Adds a layer as the substrate.
        """
        if layer is None:
            raise ValueError("The layer cannot be None")
        if self._substrate is not None:
            raise ValueError("A substrate was already defined.")

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
        hv_ = c.c_double(experiment.energy_eV / 1e3)
        iElt_ = c.c_int()
        iLine_ = c.c_int()
        iExpK_ = c.c_int()
        l.debug('StEdAddNrAtomLineHV(key, %i, %i)', experiment.z, experiment.line)
        if not self._lib.StEdAddNrAtomLineHV(self._key, nra_, klm_, hv_,
                                             c.byref(iElt_), c.byref(iLine_), c.byref(iExpK_)):
            raise StratagemError("Cannot add atomic number and line")

        standard = experiment.standard
        standard_ = c.create_string_buffer(standard.encode('ascii'))
        l.debug('StEdSetLine(key, %i, %i, %i, %s)', iElt_.value, iLine_.value, klm_.value, standard)
        if not self._lib.StEdSetLine(self._key, iElt_, iLine_, klm_, standard_):
            raise StratagemError("Cannot set standard")

        analyzed = experiment.is_analyzed()
        analyzed_ = c.c_bool(analyzed)
        l.debug("StEdSetAnalyzedFlag(key, %i, %r)", iElt_.value, analyzed)
        if not self._lib.StEdSetAnalyzedFlag(self._key, iElt_, analyzed_):
            raise StratagemError("Cannot add experiment analyzed flag")

        kratio_ = c.c_double(experiment.kratio)
        l.debug("StEdSetExpK(key, %i, %i, %i, %f, %f, %f, 0.0, 2)",
                iElt_.value, iLine_.value, iExpK_.value,
                experiment.energy_eV / 1e3, experiment.energy_eV / 1e3,
                experiment.kratio)
        if not self._lib.StEdSetExpK(self._key, iElt_, iLine_, iExpK_,
                                     hv_, hv_, kratio_, c.c_double(0.0),
                                     c.c_int(2)):
            raise StratagemError("Cannot set experiment k-ratio")

        if experiment.is_analyzed():
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
            raise StratagemError("Cannot set geometry parameters")

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
        and continuum fluorescence or no fluorescence.
        """
        flag_ = c.c_int(flag)
        l.debug('StSetFluorFlg(%i)', flag)
        self._lib.StSetFluorFlg(flag_)

    def set_standard_directory(self, dirpath):
        dirpath_ = c.create_string_buffer(dirpath.encode('ascii'))
        self._lib.StSetDirectory(c.c_int(1), dirpath_)

    def get_standard_directory(self):
        dirpath = (c.c_char * 256)()
        self._lib.StGetDirectory(c.c_int(1), c.byref(dirpath), 256)
        return dirpath.value.decode('ascii')

    def compute_kratio_vs_thickness(self, layer,
                                    thickness_low_m, thickness_high_m, step):
        """
        Computes the variation of the k-ratio as a function of the thickness 
        for a layer.
        
        :arg layer: layer (must have been previously added)
        :arg thickness_low_m: lower limit of the thickness (in nm)
        :arg thickness_high: upper limit of the thickness (in nm)
        :arg step: number of steps
        
        :return: :class:`list` of thicknesses, :class:`dict` of experiment-kratios
        """
        l.debug('StSetKvsThicknessUnit(2)')
        self._lib.StSetKvsThicknessUnit(2) # unit in nm

        if layer not in self._layers:
            raise ValueError("Unknown layer")
        iLayer = self._layers[layer]
        iLayer_ = c.c_int(iLayer)

        step_ = c.c_int(step)
        l.debug('StSetNbComputedHV(%i)', step)
        self._lib.StSetNbComputedHV(step_)

        # Compute
        low_ = c.c_double(thickness_low_m * 1e9)
        high_ = c.c_double(thickness_high_m * 1e9)
        l.debug('StComputeKvsThickness(key, %i, %f, %f)',
                iLayer, thickness_low_m * 1e9, thickness_high_m * 1e9)
        if not self._lib.StComputeKvsThickness(self._key, iLayer_, low_, high_):
            raise StratagemError("Cannot compute k-ratio vs thickness")

        # Fetch results
        thicknesses = []
        kratios = {}

        thick_ = c.c_double()
        k_ = c.c_double()
        for i in range(step + 1):
            i_ = c.c_int(i)

            if not self._lib.StGetKvsT_Thick(self._key, i_, c.byref(thick_)):
                raise StratagemError("Cannot get thickness")
            thicknesses.append(thick_.value)

            for experiment, indexes in self._experiments.items():
                iElt_ = c.c_int(indexes[0])
                iLine_ = c.c_int(indexes[1])
                iHv_ = c.c_int(indexes[2])

                if not self._lib.StGetKvsT_K(self._key, i_, iElt_, iLine_,
                                             iHv_, c.byref(k_)):
                    raise StratagemError("Cannot get k-ratio")
                kratios.setdefault(experiment, []).append(k_.value)

        return thicknesses, kratios

    def compute_kratio_vs_energy(self, energy_high_eV, step):
        """
        Computes the variation of the k-ratio as a function of the incident
        energy. 
        Note that the computation also starts at 0 keV up to the specified energy.
        
        :arg energy_high: upper limit of the thickness (in eV)
        :arg step: number of steps
        
        :return: :class:`list` of energies, :class:`dict` of experiment-kratios
        """
        step_ = c.c_int(step)
        l.debug('StSetNbComputedHV(%i)', step)
        self._lib.StSetNbComputedHV(step_)

        energy_ = c.c_double(energy_high_eV / 1e3)
        l.debug('StSetMaxHV(%f)' % (energy_high_eV / 1e3,))
        self._lib.StSetMaxHV(energy_)

        # Compute
        l.debug('StComputeKvsHV(key)')
        if not self._lib.StComputeKvsHV(self._key):
            raise StratagemError("Cannot compute k-ratio vs energy")

        # Fetch results
        energies = []
        kratios = {}

        k_ = c.c_double()
        bHV_ = c.c_bool(True)
        increment = float(energy_high_eV / 1e3) / step

        for i in range(step + 1):
            hv = i * increment
            hv_ = c.c_double(hv)

            for experiment, indexes in self._experiments.items():
                iElt_ = c.c_int(indexes[0])
                iLine_ = c.c_int(indexes[1])

                if not self._lib.StKvsHvOrRx(self._key, iElt_, iLine_, hv_, bHV_, c.byref(k_)):
                    raise StratagemError("Cannot get k-ratio")

                kratios.setdefault(experiment, []).append(k_.value)

            energies.append(hv)

        return energies, kratios

    def compute_kratios(self):
        """
        Computes the kratios of the different experiments.
        
        :return: :class:`dict` of experiment-kratios
        """
        if len(self._layers) == 0:
            return self._compute_kratios_substrate()
        else:
            return self._compute_kratios_multilayers()

    def _compute_kratios_multilayers(self):
        """
        Computes the kratios using the :meth:`compute_kratio_vs_thickness`.
        """
        for i, layer in enumerate(self._layers.keys()):
            if not layer.is_thickness_known():
                raise ValueError("Thickness of layer %i is unknown" % i)

        # Compute
        layer = list(self._layers.keys())[0]
        thickness_low_m = layer.thickness_m
        thickness_high_m = layer.thickness_m * 10
        step = 1

        _thicknesses, kratios = \
            self.compute_kratio_vs_thickness(layer, thickness_low_m,
                                             thickness_high_m, step)

        # Reorganize results
        output = {}
        for experiment, kratio in kratios.items():
            output.setdefault(experiment, kratio[0])

        return output

    def _compute_kratios_substrate(self):
        """
        Computes the kratios using the :meth:`compute_kratio_vs_energy`.
        """
        output = {}

        step = 2
        for experiment in self._experiments:
            energy_high_eV = experiment.energy_eV

            _energies, kratios = \
                self.compute_kratio_vs_energy(energy_high_eV, step)

            kratio = kratios[experiment][-1]
            if (kratio < 0): # Bug in strategem that some energy don't work
                l.warn("STRATAGem returns a negative k-ratio, re-try with energy + 1 eV")
                _energies, kratios = \
                    self.compute_kratio_vs_energy(energy_high_eV + 1.0, step)
                kratio = kratios[experiment][-1]

            output.setdefault(experiment, kratio)

        return output

    def compute(self, iteration_max=50):
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
            raise StratagemError("Cannot start iteration")

        continue_ = c.c_bool(True)
        iteration = 0

        l.debug('Start iteration')
        while True:
            iteration += 1
            l.debug('Iteration #%i' % iteration)

            l.debug('StComputeIterpNext(key, %r)' % continue_.value)
            if not self._lib.StComputeIterpNext(self._key, c.byref(continue_)):
                break

            if not continue_.value:
                break

        l.debug('Iteration completed')

        # Fetch results
        layers = {}

        thick_known = c.c_bool()
        mass_thickness = c.c_double()
        thickness = c.c_double()
        density = c.c_double()

        for layer, iLayer in self._layers.items():
            iLayer_ = c.c_int(iLayer)

            l.debug('StSdGetNbElts(key, %i)' % iLayer)
            nbelt = self._lib.StSdGetNbElts(self._key, iLayer_)
            if nbelt == -1:
                raise StratagemError("Cannot get number of elements")

            flag_ = (c.c_int * nbelt)()
            wfs_ = (c.c_double * nbelt)()
            l.debug('StSdGetLayRawConcs(key, %i, flag, wfs)' % iLayer)
            if not self._lib.StSdGetLayRawConcs(self._key, iLayer_,
                                                flag_, wfs_):
                raise StratagemError("Cannot get layer concentration")

            composition = {}
            for z in layer.composition.keys():
                nra_ = c.c_int(z)
                l.debug('StSdGetEltIdx(key, %i, %i)' % (iLayer, z))
                zindex = self._lib.StSdGetEltIdx(self._key, iLayer_, nra_)
                composition[z] = wfs_[zindex]

            l.debug("StSdGetThick(key, %i)", iLayer)
            if not self._lib.StSdGetThick(self._key, iLayer_, c.byref(thick_known),
                                          c.byref(mass_thickness), c.byref(thickness),
                                          c.byref(density)):
                raise StratagemError("Cannot get thickness")

            newlayer = Layer(composition,
                             thickness.value / 1e10,
                             mass_thickness.value * 10.0,
                             density.value * 1e3)
            layers[layer] = newlayer

        return layers

    def compute_prz(self, maxdepth_m=None, bins=100):
        """
        Compute :math:`\\phi(\\rho z)` of all experiments.
        
        .. warning::
        
           Only available for substrate (no layers).
        
        :arg maxdepth_m: maximum depth of the :math:`\\phi(\\rho z)` 
          distribution in meters. If ``None``, Kanaya-Okayama electron range
          is used with a safety factor of 1.5.
        :type maxdepth_m: :class:`float`
        
        :arg bins: number of bins in the :math:`\\phi(\\rho z)` distribution
        :type bins: :class:`int`
        
        :return: a :class:`dict` where the keys are the experiments and the 
            values are a tuple containing three lists:
            
                * :math:`\rho z` coordinates (in g/cm2)
                * generated intensities of :math:`\\phi(\\rho z)` (no absorption)
                * emitted intensites of :math:`\\phi(\\rho z)`
        """
        if len(self._layers) > 0:
            raise RuntimeError('PRZ can only be computed for substrate')

        # Set scaling
        hvs_eV = map(attrgetter('energy_eV'), self._experiments.keys())
        maxhv_eV = max(hvs_eV)
        maxhv_ = c.c_double(maxhv_eV / 1e3)
        l.debug('StSetScaleHV(%s)', maxhv_eV / 1e3)
        self._lib.StSetScaleHV(maxhv_)

        # Compute
        l.debug('StComputePrz(key)')
        if not self._lib.StComputePrz(self._key):
            raise StratagemError('Cannot compute prz')

        # Get values
        przs = {}

        for experiment, indexes in self._experiments.items():
            # Size of each bin
            if maxdepth_m is None:
                # Calculate max depth using Kanaya-Okayama
                maxdepth_m = 0.0
                energy_keV = experiment.energy_eV / 1e3

                for z, fraction in self._substrate.composition.items():
                    dr = (0.0276 * ATOMIC_MASSES[z + 1] * energy_keV ** 1.67) / \
                          (z ** 0.89 * DENSITIES[z + 1])
                    maxdepth_m += fraction / (dr * 1e-6)

                maxdepth_m = 1.0 / maxdepth_m
                maxdepth_m *= 1.5 # safety factor

            increment_kg_m2 = (maxdepth_m * self._substrate.density_kg_m3) / bins

            # Indexes
            iElt_ = c.c_int(indexes[0])
            iLine_ = c.c_int(indexes[1])
            iHV_ = c.c_int(0)

            rzs = []
            ys_generated = []
            ys_emitted = []

            for i in range(bins):
                rz_ = c.c_double(i * increment_kg_m2 * 0.1)
                rzs.append(i * increment_kg_m2)

                y_ = c.c_double()
                bUseExp_ = c.c_bool(True)
                self._lib.StPhiRhoZ(self._key, iElt_, iLine_, iHV_, rz_,
                                    bUseExp_, c.byref(y_))
                ys_emitted.append(y_.value)

                y_ = c.c_double()
                bUseExp_ = c.c_bool(False)
                self._lib.StPhiRhoZ(self._key, iElt_, iLine_, iHV_, rz_,
                                    bUseExp_, c.byref(y_))
                ys_generated.append(y_.value)

            przs.setdefault(experiment, (rzs, ys_generated, ys_emitted))

        return przs

