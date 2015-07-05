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
import ctypes as c
import math
import tempfile
import shutil

# Third party modules.

# Local modules.
from stratagemtools.stratagem import Stratagem, StratagemError
from stratagemtools.sample import Sample
from stratagemtools.experiment import Experiment

# Globals and constants variables.
from stratagemtools.experiment import LINE_KA, LINE_LA, LINE_MA
from stratagemtools.stratagem import PRZMODE_PAP, FLUORESCENCE_NONE, FLUORESCENCE_LINE_CONT

class TestStratagem(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)

        self.s = Stratagem()

        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        unittest.TestCase.tearDown(self)

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _create_sample_known_thickness(self):
        sample = Sample({79:1.0})
        sample.add_layer({13:1.0}, thickness_m=20e-9)

        exp0 = Experiment(13, LINE_KA, 25.0e3, 0.25)
        exp1 = Experiment(79, LINE_MA, 25.0e3)

        return sample, [exp0, exp1]

    def _create_sample_unknown_thickness(self):
        sample = Sample({79:1.0})
        sample.add_layer({13:1.0}, thickness_m=None)

        exp0 = Experiment(13, LINE_KA, 25.0e3, 0.008347815)
        exp1 = Experiment(79, LINE_MA, 25.0e3, 0.981343858)

        return sample, [exp0, exp1]

    def _create_sample_substrate(self):
        sample = Sample({13: 0.5, 79: 0.5})

        exp0 = Experiment(13, LINE_KA, 10.0e3)
        exp1 = Experiment(79, LINE_MA, 10.0e3)

        return sample, [exp0, exp1]

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def test_raise_error(self):
        self.assertRaises(StratagemError, self.s.compute)

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testset_sample1(self):
        sample = Sample({79:1.0})

        with self.s:
            self.s.set_sample(sample)
            self.assertEqual(sample.composition, self.s.get_sample().composition)

        self.assertRaises(StratagemError, self.s.set_sample, sample) # Not initialized

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testset_sample2(self):
        sample = Sample({1: 1.0})
        sample.add_layer({13: 0.2, 29: 0.8}, mass_thickness_kg_m2=3.55e-5,
                         density_kg_m3=2.7e3)
        sample.add_layer({6: None, 74: None}, thickness_m=None)
        sample.add_layer({79:1.0})

        with self.s:
            self.s.set_sample(sample)

            # Tests
            self.assertEqual(4, self.s._lib.StSdGetNbLayers(self.s._key))

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
            self.assertAlmostEqual(0.0, thickness.value, 3)
            self.assertAlmostEqual(0.0, massThickness.value, 3)
            self.assertAlmostEqual(10.0, density.value, 3)

            self.s._lib.StSdGetThick(self.s._key, c.c_int(2), c.byref(thickKnown),
                                     c.byref(massThickness), c.byref(thickness),
                                     c.byref(density))
            self.assertFalse(thickKnown.value)
            self.assertAlmostEqual(0.0, thickness.value, 3)
            self.assertAlmostEqual(0.0, massThickness.value, 3)
            self.assertAlmostEqual(19.30, density.value, 3)

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testadd_experiment(self):
        with self.s:
            self.s.add_experiment(Experiment(29, LINE_KA, 25.0e3, 0.5))
            self.s.add_experiment(Experiment(29, LINE_KA, 20.0e3, 0.35))

            self.s.add_experiment(Experiment(79, LINE_LA, 20.0e3, analyzed=False))

            self.s.add_experiment(Experiment(79, LINE_MA, 5.0e3, analyzed=False))
            self.s.add_experiment(Experiment(79, LINE_MA, 20.0e3, analyzed=False))

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

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testgeometry(self):
        with self.s:
            self.s.set_geometry(0.1, 0.2, 0.3)
            toa, tilt, azimuth = self.s.get_geometry()

        self.assertAlmostEqual(0.1, toa, 3)
        self.assertAlmostEqual(0.2, tilt, 3)
        self.assertAlmostEqual(0.3, azimuth, 3)

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testprz_mode(self):
        with self.s:
            self.s.set_prz_mode(PRZMODE_PAP)
            self.assertEqual(1, self.s.get_prz_mode())

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testfluorescence(self):
        with self.s:
            self.s.set_fluorescence(FLUORESCENCE_NONE)
            self.assertEqual(0, self.s.get_fluorescence())

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testcompute_kratio_vs_thickness(self):
        sample, (exp0, exp1) = self._create_sample_known_thickness()

        with self.s:
            self.s.set_sample(sample)
            self.s.add_experiment(exp0)
            self.s.add_experiment(exp1)

            self.s.set_geometry(math.radians(40.0), 0.0, 0.0)

            thicknesses, kratios = \
                self.s.compute_kratio_vs_thickness(sample.get_layer(0), 0.0, 100.0e-9, 10)

        self.assertEqual(2, self.s._lib.StGetKvsThicknessUnit())

        self.assertEqual(11, len(thicknesses))
        self.assertEqual(2, len(kratios))

        self.assertEqual(11, len(kratios[exp0]))
        self.assertEqual(11, len(kratios[exp1]))

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testcompute_kratio_vs_energy(self):
        sample, (exp0, exp1) = self._create_sample_known_thickness()

        with self.s:
            self.s.set_sample(sample)
            self.s.add_experiment(exp0)
            self.s.add_experiment(exp1)

            self.s.set_geometry(math.radians(40.0), 0.0, 0.0)

            energies, kratios = self.s.compute_kratio_vs_energy(30.0e3, 10)

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

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testcompute_kratios(self):
        sample, (exp0, exp1) = self._create_sample_known_thickness()

        with self.s:
            self.s.set_sample(sample)
            self.s.add_experiment(exp0)
            self.s.add_experiment(exp1)

            self.s.set_geometry(math.radians(40.0), 0.0, 0.0)

            kratios = self.s.compute_kratios()

        self.assertEqual(2, len(kratios))
        self.assertAlmostEqual(0.008347815, kratios[exp0], 3)
        self.assertAlmostEqual(0.981343858, kratios[exp1], 3)

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testcompute_kratios_substrate(self):
        sample, (exp0, exp1) = self._create_sample_substrate()

        with self.s:
            self.s.set_sample(sample)
            self.s.add_experiment(exp0)
            self.s.add_experiment(exp1)

            self.s.set_geometry(math.radians(40.0), 0.0, 0.0)

            kratios = self.s.compute_kratios()

        self.assertEqual(2, len(kratios))
        self.assertAlmostEqual(0.4947, kratios[exp0], 3)
        self.assertAlmostEqual(0.36845, kratios[exp1], 3)

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testcompute(self):
        sample, (exp0, exp1) = self._create_sample_unknown_thickness()

        with self.s:
            self.s.set_sample(sample)
            self.s.add_experiment(exp0)
            self.s.add_experiment(exp1)

            self.s.set_geometry(math.radians(40.0), 0.0, 0.0)

            newsample = self.s.compute()

        self.assertEqual(1, len(newsample.layers))
        self.assertAlmostEqual(20.0e-9, newsample.get_layer(0).thickness_m, 9)
        self.assertAlmostEqual(1.0, newsample.get_layer(0).composition[13], 4)

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def testcompute_prz(self):
        sample, (exp0, exp1) = self._create_sample_substrate()

        with self.s:
            self.s.set_sample(sample)
            self.s.add_experiment(exp0)
            self.s.add_experiment(exp1)

            self.s.set_fluorescence(FLUORESCENCE_LINE_CONT)
            self.s.set_prz_mode(PRZMODE_PAP)
            self.s.set_geometry(math.radians(40.0), 0.0, 0.0)

            przs = self.s.compute_prz(None, 100)

        rzs, generated, emitted = przs[exp0]
        self.assertEqual(100, len(rzs))
        self.assertEqual(100, len(generated))
        self.assertEqual(100, len(emitted))

        self.assertAlmostEqual(0.0, rzs[0], 4)
        self.assertAlmostEqual(1.969945, generated[0], 4)
        self.assertAlmostEqual(1.969945, emitted[0], 4)

        self.assertAlmostEqual(1.24978e-5, rzs[1], 4)
        self.assertAlmostEqual(2.0575, generated[1], 4)
        self.assertAlmostEqual(2.0438, emitted[1], 4)

    @unittest.skipUnless(os.name == 'nt', 'Test can only be ran under Windows platform')
    def teststandard_sample(self):
        sample = Sample({79:1.0})
        sample.add_layer({13:1.0}, thickness_m=None)

        std = Sample({13: 0.5293, 8: 0.4707})
        exp0 = Experiment(13, LINE_KA, 25.0e3, 0.008347815, standard=std)
        exp1 = Experiment(79, LINE_MA, 25.0e3, 0.981343858)

        with self.s:
            self.s.set_standard_directory(self.tmpdir)

            self.s.set_sample(sample)
            self.s.add_experiment(exp0)
            self.s.add_experiment(exp1)

            self.s.set_geometry(math.radians(40.0), 0.0, 0.0)

            newsample = self.s.compute()

        self.assertEqual(1, len(newsample.layers))
        self.assertAlmostEqual(13.59e-9, newsample.get_layer(0).thickness_m, 9)
        self.assertAlmostEqual(1.0, newsample.get_layer(0).composition[13], 4)

if __name__ == '__main__': #pragma: no cover
    logging.getLogger().setLevel(logging.DEBUG)
    unittest.main()
