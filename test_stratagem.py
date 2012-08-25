#!/usr/bin/env python
""" """

# Script information for the file.
__author__ = "Philippe T. Pinard"
__email__ = "philippe.pinard@gmail.com"
__version__ = "0.1"
__copyright__ = "Copyright (c) 2011 Philippe T. Pinard"
__license__ = "GPL v3"

# Standard library modules.
import unittest
import logging
import os
import sys
import ctypes as c
import math

# Third party modules.

# Local modules.
import DrixUtilities.Files as Files

from stratagem import Stratagem
from layer import Layer
from experiment import Experiment

# Globals and constants variables.
from experiment import LINE_KA, LINE_LA, LINE_MA
from stratagem import PRZMODE_PAP, FLUORESCENCE_NONE

path = Files.getCurrentModulePath(__file__, 'stratagem.cfg')
path = os.path.abspath(path)

class TestStratagem(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)

        configfile = open(path, 'r')
        self.s = Stratagem(configfile)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

        self.s.close()

    def testskeleton(self):
        self.assertTrue(True)

    def _setup_known_thickness(self):
        film = Layer({13:1.0}, thickness=20)
        self.s.add_layer(film)

        subs = Layer({79:1.0})
        self.s.add_substrate(subs)

        exp0 = Experiment(13, LINE_KA, 25.0, 0.25)
        self.s.add_experiment(exp0)

        exp1 = Experiment(79, LINE_MA, 25.0)
        self.s.add_experiment(exp1)

        self.s.set_geometry(40.0 / 180 * math.pi, 0.0, 0.0)

        return film, subs, exp0, exp1

    def _setup_unknown_thickness(self):
        film = Layer({13:1.0}, thickness= -1.0)
        self.s.add_layer(film)

        subs = Layer({79:1.0})
        self.s.add_substrate(subs)

        exp0 = Experiment(13, LINE_KA, 25.0, 0.008347815)
        self.s.add_experiment(exp0)

        exp1 = Experiment(79, LINE_MA, 25.0, 0.981343858)
        self.s.add_experiment(exp1)

        self.s.set_geometry(40.0 / 180 * math.pi, 0.0, 0.0)

        return film, subs, exp0, exp1

    def _setup_substrate(self):
        subs = Layer({13: 0.5, 79: 0.5})
        self.s.add_substrate(subs)

        exp0 = Experiment(13, LINE_KA, 10.0)
        self.s.add_experiment(exp0)

        exp1 = Experiment(79, LINE_MA, 10.0)
        self.s.add_experiment(exp1)

        self.s.set_geometry(40.0 / 180 * math.pi, 0.0, 0.0)

        return subs, exp0, exp1

    def testadd_layer(self):
        layer = Layer({13: 0.2, 29: 0.8}, mass_thickness=3.55, density=2.7)
        self.s.add_layer(layer)

        layer = Layer({6:-1.0, 74:-1.0}, thickness= -1.0)
        self.s.add_layer(layer)

        layer = Layer({79:1.0})
        self.s.add_layer(layer)

        # Tests
        self.assertEqual(3, self.s._lib.StSdGetNbLayers(self.s._key))

        self.assertEqual(2, self.s._lib.StSdGetNbElts(self.s._key, c.c_int(0)))
        self.assertEqual(2, self.s._lib.StSdGetNbElts(self.s._key, c.c_int(1)))
        self.assertEqual(1, self.s._lib.StSdGetNbElts(self.s._key, c.c_int(2)))

        self.assertEqual(13, self.s._lib.StSdGetNrAtom(self.s._key, c.c_int(0), c.c_int(0)))
        self.assertEqual(29, self.s._lib.StSdGetNrAtom(self.s._key, c.c_int(0), c.c_int(1)))
        self.assertEqual(74, self.s._lib.StSdGetNrAtom(self.s._key, c.c_int(1), c.c_int(0)))
        self.assertEqual(6, self.s._lib.StSdGetNrAtom(self.s._key, c.c_int(1), c.c_int(1)))
        self.assertEqual(79, self.s._lib.StSdGetNrAtom(self.s._key, c.c_int(2), c.c_int(0)))

        flags = (c.c_int * 2)()
        wfs = (c.c_double * 2)()
        afs = (c.c_double * 2)()
        self.s._lib.StSdGetLayConcs(self.s._key, c.c_int(0), flags, wfs, afs)
        self.assertEqual(0, flags[0])
        self.assertEqual(0, flags[1])
        self.assertEqual(0.2, wfs[0])
        self.assertEqual(0.8, wfs[1])

        self.s._lib.StSdGetLayConcs(self.s._key, c.c_int(1), flags, wfs, afs)
        self.assertEqual(1, flags[0])
        self.assertEqual(1, flags[1])

        flags = (c.c_int * 1)()
        wfs = (c.c_double * 1)()
        afs = (c.c_double * 1)()
        self.s._lib.StSdGetLayConcs(self.s._key, c.c_int(2), flags, wfs, afs)
        self.assertEqual(0, flags[0])
        self.assertEqual(1.0, wfs[0])

        thickKnown = c.c_bool()
        massThickness = c.c_double()
        thickness = c.c_double()
        density = c.c_double()
        self.s._lib.StSdGetThick(self.s._key, c.c_int(0), c.byref(thickKnown),
                                 c.byref(massThickness), c.byref(thickness),
                                 c.byref(density))
        self.assertTrue(thickKnown.value)
        self.assertAlmostEqual(13.14815, thickness.value / 10, 3)
        self.assertAlmostEqual(3.55, massThickness.value * 1e6, 3)
        self.assertAlmostEqual(2.70, density.value, 3)

        self.s._lib.StSdGetThick(self.s._key, c.c_int(1), c.byref(thickKnown),
                                 c.byref(massThickness), c.byref(thickness),
                                 c.byref(density))
        self.assertFalse(thickKnown.value)
        self.assertAlmostEqual(-1.0, thickness.value / 10, 3)
        self.assertAlmostEqual(-1.0, massThickness.value * 1e6, 3)
        self.assertAlmostEqual(0.0, density.value, 3)

        self.s._lib.StSdGetThick(self.s._key, c.c_int(2), c.byref(thickKnown),
                                 c.byref(massThickness), c.byref(thickness),
                                 c.byref(density))
        self.assertFalse(thickKnown.value)
        self.assertAlmostEqual(-1.0, thickness.value / 10, 3)
        self.assertAlmostEqual(-1.0, massThickness.value * 1e6, 3)
        self.assertAlmostEqual(19.30, density.value, 3)

    def testadd_experiment(self):
        self.s.add_experiment(Experiment(29, LINE_KA, 25.0, 0.5))
        self.s.add_experiment(Experiment(29, LINE_KA, 20.0, 0.35))

        self.s.add_experiment(Experiment(79, LINE_LA, 20.0, analyzed=False))

        self.s.add_experiment(Experiment(79, LINE_MA, 5.0, analyzed=False))
        self.s.add_experiment(Experiment(79, LINE_MA, 20.0, analyzed=False))

        # Tests
        self.assertEqual(2, self.s._lib.StEdGetNbElts(self.s._key))

        self.assertEqual(29, self.s._lib.StEdGetNrAtom(self.s._key, c.c_int(0)))
        self.assertEqual(79, self.s._lib.StEdGetNrAtom(self.s._key, c.c_int(1)))

        self.assertEqual(1, self.s._lib.StEdGetNbLines(self.s._key, c.c_int(0)))
        self.assertEqual(2, self.s._lib.StEdGetNbLines(self.s._key, c.c_int(1)))

        line_ = c.c_int()
        stdName_ = (c.c_char * 1)()
        self.s._lib.StEdGetLine(self.s._key, c.c_int(0), c.c_int(0),
                                c.byref(line_), c.byref(stdName_), c.c_int())
        self.assertEqual(0, line_.value)

        self.s._lib.StEdGetLine(self.s._key, c.c_int(1), c.c_int(0),
                                c.byref(line_), c.byref(stdName_), c.c_int())
        self.assertEqual(2, line_.value)

        self.s._lib.StEdGetLine(self.s._key, c.c_int(1), c.c_int(1),
                                c.byref(line_), c.byref(stdName_), c.c_int())
        self.assertEqual(4, line_.value)

        analyzed_ = c.c_bool()
        self.s._lib.StEdGetAnalyzedFlag(self.s._key, c.c_int(0), c.byref(analyzed_))
        self.assertTrue(analyzed_.value)

        self.s._lib.StEdGetAnalyzedFlag(self.s._key, c.c_int(1), c.byref(analyzed_))
        self.assertFalse(analyzed_.value)

        self.assertEqual(2, self.s._lib.StEdGetNbExpKs(self.s._key, c.c_int(0), c.c_int(0)))
        self.assertEqual(1, self.s._lib.StEdGetNbExpKs(self.s._key, c.c_int(1), c.c_int(0)))
        self.assertEqual(2, self.s._lib.StEdGetNbExpKs(self.s._key, c.c_int(1), c.c_int(1)))

        hv_ = c.c_double()
        stdHv_ = c.c_double()
        k_ = c.c_double()
        kr_ = c.c_double()
        self.s._lib.StEdGetExpK(self.s._key, c.c_int(0), c.c_int(0), c.c_int(1),
                                c.byref(hv_), c.byref(stdHv_), c.byref(k_),
                                c.byref(kr_), c.byref(c.c_int()))
        self.assertAlmostEqual(25.0, hv_.value, 3)
        self.assertAlmostEqual(0.5, k_.value, 3)

        self.s._lib.StEdGetExpK(self.s._key, c.c_int(0), c.c_int(0), c.c_int(0),
                                c.byref(hv_), c.byref(stdHv_), c.byref(k_),
                                c.byref(kr_), c.byref(c.c_int()))
        self.assertAlmostEqual(20.0, hv_.value, 3)
        self.assertAlmostEqual(0.35, k_.value, 3)

        self.s._lib.StEdGetExpK(self.s._key, c.c_int(1), c.c_int(0), c.c_int(0),
                                c.byref(hv_), c.byref(stdHv_), c.byref(k_),
                                c.byref(kr_), c.byref(c.c_int()))
        self.assertAlmostEqual(20.0, hv_.value, 3)
        self.assertAlmostEqual(0.0, k_.value, 3)

        self.s._lib.StEdGetExpK(self.s._key, c.c_int(1), c.c_int(1), c.c_int(0),
                                c.byref(hv_), c.byref(stdHv_), c.byref(k_),
                                c.byref(kr_), c.byref(c.c_int()))
        self.assertAlmostEqual(5.0, hv_.value, 3)
        self.assertAlmostEqual(0.0, k_.value, 3)

        self.s._lib.StEdGetExpK(self.s._key, c.c_int(1), c.c_int(1), c.c_int(1),
                                c.byref(hv_), c.byref(stdHv_), c.byref(k_),
                                c.byref(kr_), c.byref(c.c_int()))
        self.assertAlmostEqual(20.0, hv_.value, 3)
        self.assertAlmostEqual(0.0, k_.value, 3)

    def testset_geometry(self):
        self.s.set_geometry(0.1, 0.2, 0.3)

        toa_ = c.c_double()
        tilt_ = c.c_double()
        azimuth_ = c.c_double()
        self.s._lib.StGetGeomParams(self.s._key, c.byref(toa_), c.byref(tilt_),
                                    c.byref(azimuth_))
        self.assertAlmostEqual(0.1, toa_.value, 3)
        self.assertAlmostEqual(0.2, tilt_.value, 3)
        self.assertAlmostEqual(0.3, azimuth_.value, 3)

    def testset_prz_mode(self):
        self.s.set_prz_mode(PRZMODE_PAP)
        self.assertEqual(1, self.s._lib.StGetPrzMode())

    def testset_fluorescence(self):
        self.s.set_fluorescence(FLUORESCENCE_NONE)
        self.assertEqual(0, self.s._lib.StGetFluorFlg())

    def testcompute_kratio_vs_thickness(self):
        film, _subs, exp0, exp1 = self._setup_known_thickness()

        thicknesses, kratios = \
            self.s.compute_kratio_vs_thickness(film, 0.0, 100.0, 10)

        self.assertEqual(2, self.s._lib.StGetKvsThicknessUnit())

        self.assertEqual(11, len(thicknesses))
        self.assertEqual(2, len(kratios))

        self.assertEqual(11, len(kratios[exp0]))
        self.assertEqual(11, len(kratios[exp1]))

    def testcompute_kratio_vs_energy(self):
        _film, _subs, exp0, exp1 = self._setup_known_thickness()

        energies, kratios = self.s.compute_kratio_vs_energy(30.0, 10)

        exp0_kratios = kratios[exp0]
        exp1_kratios = kratios[exp1]

        self.assertEqual(11, len(energies))
        self.assertEqual(2, len(kratios))
        self.assertEqual(11, len(exp0_kratios))
        self.assertEqual(11, len(exp1_kratios))

        self.assertAlmostEqual(-1.0, exp0_kratios[0], 3)
        self.assertAlmostEqual(-1.0, exp1_kratios[0], 3)

        self.assertAlmostEqual(0.0076135, exp0_kratios[-2], 5)
        self.assertAlmostEqual(0.9819639, exp1_kratios[-2], 5)

    def testcompute_kratios(self):
        _film, _subs, exp0, exp1 = self._setup_known_thickness()

        kratios = self.s.compute_kratios()

        self.assertEqual(2, len(kratios))
        self.assertAlmostEqual(0.008347815, kratios[exp0], 3)
        self.assertAlmostEqual(0.981343858, kratios[exp1], 3)

    def testcompute_kratios_substrate(self):
        _subs, exp0, exp1 = self._setup_substrate()

        kratios = self.s.compute_kratios()

        self.assertEqual(2, len(kratios))
        self.assertAlmostEqual(0.4947, kratios[exp0], 3)
        self.assertAlmostEqual(0.36845, kratios[exp1], 3)

    def testcompute_thicknesses(self):
        film, _subs, _exp0, _exp1 = self._setup_unknown_thickness()

        thicknesses = self.s.compute_thicknesses()

        self.assertEqual(1, len(thicknesses))
        self.assertAlmostEqual(20.0, thicknesses[film], 1)

    def testcompute_prz(self):
        _subs, exp0, _exp1 = self._setup_substrate()

        przs = self.s.compute_prz(None, 100)

        rzs, generated, emitted = przs[exp0]
        self.assertEqual(100, len(rzs))
        self.assertEqual(100, len(generated))
        self.assertEqual(100, len(emitted))

        self.assertAlmostEqual(0.0, rzs[0], 4)
        self.assertAlmostEqual(1.969945, generated[0], 4)
        self.assertAlmostEqual(1.969945, emitted[0], 4)

        self.assertAlmostEqual(1.24978e-5, rzs[1], 4)
        self.assertAlmostEqual(2.23834, generated[1], 4)
        self.assertAlmostEqual(2.18035, emitted[1], 4)

if __name__ == '__main__': #pragma: no cover
    if os.name == 'nt':
        logging.getLogger().setLevel(logging.DEBUG)
        unittest.main()
    else:
        print >> sys.stderr, 'Tests can only be ran under Windows platform'

