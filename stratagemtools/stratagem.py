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
__email__ = "philippe.pinard@gmailogger.com"
__version__ = "0.1"
__copyright__ = "Copyright (c) 2011 Philippe T. Pinard"
__license__ = "GPL v3"

# Standard library modules.
import os
import ctypes as c
import logging
logger = logging.getLogger(__name__)
from operator import attrgetter
try:
    import winreg
except ImportError:
    try:
        import _winreg as winreg
    except ImportError:
        class winreg:

            HKEY_CURRENT_USER = None

            class _PyHKEY(object):

                def __enter__(self):
                    return self

                def __exit__(self, exc_type, exc_value, traceback):
                    pass

            def OpenKey(self, key, sub_key, res, sam):
                return self._PyHKEY()

            def QueryValueEx(self, key, value_name):
                return None

# Third party modules.

# Local modules.
from stratagemtools.sample import Sample
from stratagemtools.element_properties import \
    atomic_mass_kg_mol, mass_density_kg_m3

# Globals and constants variables.
_REGISTRY_KEY = "Software\SAMx\Stratagem\Configuration"
_REGISTRY_VALUENAME = 'InstallOEMDirectory'

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
    def __init__(self, dll_path=None):
        """
        Initializes the connection to the Stratagem DLlogger.
        One of the following argument must be specified:
        
        :arg configfile: file-object to the configuration file
        :arg dll_path: complete path to the location of the ``stratadllogger.dll``
        """
        if dll_path is None:
            with winreg.OpenKey(winreg.HKEY_CURRENT_USER, _REGISTRY_KEY) as key: #@UndefinedVariable
                basedir = winreg.QueryValueEx(key, _REGISTRY_VALUENAME)[0] #@UndefinedVariable
            dll_path = os.path.join(basedir, 'bin', 'stratadll.dll')

        logger.debug("dll=%s", dll_path)
        self._lib = c.WinDLL(dll_path)

        self._key = c.create_string_buffer(b'DemoSample')
        self._stobjectnew(self._key)
        self._stenableerrordisplay(False)

        self.reset()

    def _stobjectnew(self, key):
        bNormal_ = c.c_bool(True)
        iniFlags_ = c.c_int(0)

        logger.debug("StObjectNew(key, %r, %i)", True, 0)
        if not self._lib.StObjectNew(self._key, bNormal_, iniFlags_):
            self._raise_error("Cannot create object")

    def _stenableerrordisplay(self, enable):
        enable_ = c.c_bool(enable)

        logger.debug("StEnableErrorDisplay(%r)", enable)
        self._lib.StEnableErrorDisplay(enable_)

    def _raise_error(self, alternate=''):
        errnum_ = c.c_ulong()
        errtype_ = c.c_int()

        self._lib.StGetLastError(c.byref(errnum_), c.byref(errtype_))

        if errnum_.value != 0:
            if errtype_.value == 0:
                buf_ = c.create_string_buffer(256)
                self._lib.StGetMsg(errnum_, buf_, 256)
                raise StratagemError(buf_.value.decode('ascii'))
            elif errtype_.value == 1:
                raise c.WinError(errtype_.value)
            else:
                raise StratagemError('Error %i' % errnum_.value)
        else:
            raise StratagemError(alternate)

    def close(self):
        """
        Closes the connection to the Statagem DLlogger.
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

    def set_sample(self, sample):
        self.reset()

        for layer in sample.layers:
            self._add_layer(layer, substrate=False)

        self._add_layer(sample.substrate, substrate=True)

    def get_sample(self):
        sample = Sample(self._substrate.composition)

        for layer in self._layers:
            sample.add_layer(layer.composition, layer.thickness_m,
                             layer.mass_thickness_kg_m2, layer.density_kg_m3)

        return sample

    def _add_layer(self, layer, substrate=False):
        """
        Adds a layer from top to bottom. 
        The last layer added is considered as the substrate.
        
        :arg layer: layer
        :type layer: :class:`Layer`
        
        :return: index of the layer
        """
        logger.debug("StSdAddLayer(key)")
        iLayer_ = self._lib.StSdGetNbLayers(self._key)

        logger.debug("StSdAddLayer(key, %i)", iLayer_)
        if not self._lib.StSdAddLayer(self._key, iLayer_):
            self._raise_error("Cannot add layer")

        for i, value in enumerate(layer.composition.items()):
            iElt_ = c.c_int(i)
            logger.debug("StSdAddElt(key, %i, %i)", iLayer_, i)
            if not self._lib.StSdAddElt(self._key, iLayer_, iElt_):
                self._raise_error("Cannot add element")

            z, wf = value
            nra_ = c.c_int(z)
            logger.debug("StSdSetNrAtom(key, %i, %i, %i)", iLayer_, i, z)
            if not self._lib.StSdSetNrAtom(self._key, iLayer_, iElt_, nra_):
                self._raise_error("Cannot set atomic number")

            if wf is not None:
                flag = CONCENTRATION_FLAG_KNOWN

                wf_ = c.c_double(wf)
                logger.debug("StSdSetConc(key, %i, %i, %f)", iLayer_, i, wf)
                if not self._lib.StSdSetConc(self._key, iLayer_, iElt_, wf_):
                    self._raise_error("Cannot set concentration")
            else:
                flag = CONCENTRATION_FLAG_UNKNOWN

            logger.debug("StSdSetConcFlag(key, %i, %i, %i)", iLayer_, i, flag)
            if not self._lib.StSdSetConcFlag(self._key, iLayer_, iElt_, c.c_int(flag)):
                self._raise_error("Cannot set concentration flag")

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

            logger.debug("StSdSetThick(key, %i, %r, %d, %d, %d)", iLayer_,
                    thick_known, mass_thickness, thickness, density)
            if not self._lib.StSdSetThick(self._key, iLayer_, thick_known_,
                                          mass_thickness_, thickness_, density_):
                self._raise_error("Cannot set thickness")

            self._layers.setdefault(layer, int(iLayer_))
        else:
            self._substrate = layer

    def _create_standard(self, composition):
        pass

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
        logger.debug('StEdAddNrAtomLineHV(key, %i, %i)', experiment.z, experiment.line)
        if not self._lib.StEdAddNrAtomLineHV(self._key, nra_, klm_, hv_,
                                             c.byref(iElt_), c.byref(iLine_), c.byref(iExpK_)):
            self._raise_error("Cannot add atomic number and line")

        standard = experiment.standard
        if isinstance(standard, dict):
            standard = self._create_standard(standard)
        standard_ = c.create_string_buffer(standard.encode('ascii'))
        logger.debug('StEdSetLine(key, %i, %i, %i, %s)', iElt_.value, iLine_.value, klm_.value, standard)
        if not self._lib.StEdSetLine(self._key, iElt_, iLine_, klm_, standard_):
            self._raise_error("Cannot set standard")

        analyzed = experiment.is_analyzed()
        analyzed_ = c.c_bool(analyzed)
        logger.debug("StEdSetAnalyzedFlag(key, %i, %r)", iElt_.value, analyzed)
        if not self._lib.StEdSetAnalyzedFlag(self._key, iElt_, analyzed_):
            self._raise_error("Cannot add experiment analyzed flag")

        kratio_ = c.c_double(experiment.kratio)
        logger.debug("StEdSetExpK(key, %i, %i, %i, %f, %f, %f, 0.0, 2)",
                iElt_.value, iLine_.value, iExpK_.value,
                experiment.energy_eV / 1e3, experiment.energy_eV / 1e3,
                experiment.kratio)
        if not self._lib.StEdSetExpK(self._key, iElt_, iLine_, iExpK_,
                                     hv_, hv_, kratio_, c.c_double(0.0),
                                     c.c_int(2)):
            self._raise_error("Cannot set experiment k-ratio")

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
        logger.debug('StSetGeomParams(key, %f, %f, %f)', toa, tilt, azimuth)
        if not self._lib.StSetGeomParams(self._key, toa_, tilt_, azimuth_):
            self._raise_error("Cannot set geometry parameters")

    def set_prz_mode(self, mode):
        """
        Sets the type of model to use for the phi-rho-z.
        """
        mode_ = c.c_int(mode)
        logger.debug('StSetPrzMode(%i)', mode)
        self._lib.StSetPrzMode(mode_)

    def set_fluorescence(self, flag):
        """
        Sets whether to consider characteristic fluorescence, characteristic
        and continuum fluorescence or no fluorescence.
        """
        flag_ = c.c_int(flag)
        logger.debug('StSetFluorFlg(%i)', flag)
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
        logger.debug('StSetKvsThicknessUnit(2)')
        self._lib.StSetKvsThicknessUnit(2) # unit in nm

        if layer not in self._layers:
            raise ValueError("Unknown layer")
        iLayer = self._layers[layer]
        iLayer_ = c.c_int(iLayer)

        step_ = c.c_int(step)
        logger.debug('StSetNbComputedHV(%i)', step)
        self._lib.StSetNbComputedHV(step_)

        # Compute
        low_ = c.c_double(thickness_low_m * 1e9)
        high_ = c.c_double(thickness_high_m * 1e9)
        logger.debug('StComputeKvsThickness(key, %i, %f, %f)',
                iLayer, thickness_low_m * 1e9, thickness_high_m * 1e9)
        if not self._lib.StComputeKvsThickness(self._key, iLayer_, low_, high_):
            self._raise_error("Cannot compute k-ratio vs thickness")

        # Fetch results
        thicknesses = []
        kratios = {}

        thick_ = c.c_double()
        k_ = c.c_double()
        for i in range(step + 1):
            i_ = c.c_int(i)

            if not self._lib.StGetKvsT_Thick(self._key, i_, c.byref(thick_)):
                self._raise_error("Cannot get thickness")
            thicknesses.append(thick_.value)

            for experiment, indexes in self._experiments.items():
                iElt_ = c.c_int(indexes[0])
                iLine_ = c.c_int(indexes[1])
                iHv_ = c.c_int(indexes[2])

                if not self._lib.StGetKvsT_K(self._key, i_, iElt_, iLine_,
                                             iHv_, c.byref(k_)):
                    self._raise_error("Cannot get k-ratio")
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
        logger.debug('StSetNbComputedHV(%i)', step)
        self._lib.StSetNbComputedHV(step_)

        energy_ = c.c_double(energy_high_eV / 1e3)
        logger.debug('StSetMaxHV(%f)' % (energy_high_eV / 1e3,))
        self._lib.StSetMaxHV(energy_)

        # Compute
        logger.debug('StComputeKvsHV(key)')
        if not self._lib.StComputeKvsHV(self._key):
            self._raise_error("Cannot compute k-ratio vs energy")

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
                    self._raise_error("Cannot get k-ratio")

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
                logger.warn("STRATAGem returns a negative k-ratio, re-try with energy + 1 eV")
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
        logger.debug('StSetMaxNbIter(%i)', iteration_max)
        self._lib.StSetMaxNbIter(iteration_max_)

        # Compute
        logger.debug('StComputeIterpStart(key)')
        if not self._lib.StComputeIterpStart(self._key):
            self._raise_error()
            self._raise_error("Cannot start iteration")

        continue_ = c.c_bool(True)
        iteration = 0

        logger.debug('Start iteration')
        while True:
            iteration += 1
            logger.debug('Iteration #%i' % iteration)

            logger.debug('StComputeIterpNext(key, %r)' % continue_.value)
            if not self._lib.StComputeIterpNext(self._key, c.byref(continue_)):
                break

            if not continue_.value:
                break

        logger.debug('Iteration completed')

        # Fetch results
        thick_known = c.c_bool()
        mass_thickness = c.c_double()
        thickness = c.c_double()
        density = c.c_double()

        def get_layer(ilayer, layer):
            ilayer_ = c.c_int(ilayer)

            logger.debug('StSdGetNbElts(key, %i)' % ilayer)
            nbelt = self._lib.StSdGetNbElts(self._key, ilayer_)
            if nbelt == -1:
                self._raise_error("Cannot get number of elements")

            flag_ = (c.c_int * nbelt)()
            wfs_ = (c.c_double * nbelt)()
            logger.debug('StSdGetLayRawConcs(key, %i, flag, wfs)' % ilayer)
            if not self._lib.StSdGetLayRawConcs(self._key, ilayer_,
                                                flag_, wfs_):
                self._raise_error("Cannot get layer concentration")

            composition = {}
            for z in layer.composition.keys():
                nra_ = c.c_int(z)
                logger.debug('StSdGetEltIdx(key, %i, %i)' % (ilayer, z))
                zindex = self._lib.StSdGetEltIdx(self._key, ilayer_, nra_)
                composition[z] = wfs_[zindex]

            logger.debug("StSdGetThick(key, %i)", ilayer)
            if not self._lib.StSdGetThick(self._key, ilayer_, c.byref(thick_known),
                                          c.byref(mass_thickness), c.byref(thickness),
                                          c.byref(density)):
                self._raise_error("Cannot get thickness")

            return (composition, thickness.value / 1e10,
                    mass_thickness.value * 10.0, density.value * 1e3)

        sample = Sample(get_layer(len(self._layers), self._substrate)[0])

        for layer, ilayer in self._layers.items():
            sample.add_layer(*get_layer(ilayer, layer))

        return sample

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
        logger.debug('StSetScaleHV(%s)', maxhv_eV / 1e3)
        self._lib.StSetScaleHV(maxhv_)

        # Compute
        logger.debug('StComputePrz(key)')
        if not self._lib.StComputePrz(self._key):
            self._raise_error('Cannot compute prz')

        # Get values
        przs = {}

        for experiment, indexes in self._experiments.items():
            # Size of each bin
            if maxdepth_m is None:
                # Calculate max depth using Kanaya-Okayama
                maxdepth_m = 0.0
                energy_keV = experiment.energy_eV / 1e3

                for z, fraction in self._substrate.composition.items():
                    dr = (0.0276 * atomic_mass_kg_mol(z) * 1e3 * energy_keV ** 1.67) / \
                          (z ** 0.89 * mass_density_kg_m3(z) / 1e3)
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

