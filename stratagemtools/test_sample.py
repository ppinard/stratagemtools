#!/usr/bin/env python
""" """

# Standard library modules.
import unittest
import logging

# Third party modules.

# Local modules.
from stratagemtools.sample import Sample, composition_from_formula

# Globals and constants variables.

class TestSample(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)

        self.sample = Sample({1: 1.0})

        composition = {29: 0.5, 30: 0.4}
        density_kg_m3 = 7.336 * 1000.0
        self.l1 = self.sample.add_layer(composition, thickness_m=5.0e-4, density_kg_m3=density_kg_m3)
        self.l2 = self.sample.add_layer(composition, mass_thickness_kg_m2=3.668, density_kg_m3=density_kg_m3)
        self.l3 = self.sample.add_layer(composition, density_kg_m3=density_kg_m3)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def testskeleton(self):
        self.assertEqual(2, len(list(self.l1.composition)))
        self.assertAlmostEqual(5.0e-4, self.l1.thickness_m, 3)
        self.assertAlmostEqual(3.668, self.l1.mass_thickness_kg_m2, 3)
        self.assertAlmostEqual(7336, self.l1.density_kg_m3, 3)
        self.assertTrue(self.l1.is_thickness_known())
        self.assertTrue(self.l1.is_composition_known())

        self.assertEqual(2, len(list(self.l2.composition)))
        self.assertAlmostEqual(5.0e-4, self.l2.thickness_m, 3)
        self.assertAlmostEqual(3.668, self.l2.mass_thickness_kg_m2, 3)
        self.assertAlmostEqual(7336, self.l2.density_kg_m3, 3)
        self.assertTrue(self.l2.is_thickness_known())
        self.assertTrue(self.l2.is_composition_known())

        self.assertEqual(2, len(list(self.l3.composition)))
        self.assertIsNone(self.l3.thickness_m)
        self.assertIsNone(self.l3.mass_thickness_kg_m2)
        self.assertAlmostEqual(7336, self.l3.density_kg_m3, 3)
        self.assertFalse(self.l3.is_thickness_known())
        self.assertTrue(self.l3.is_composition_known())

    def testis_composition_known(self):
        s = Sample({8: None, 10: 0.5, 12: '?'})
        self.assertFalse(s.substrate.is_composition_known())

    def testcomposition_from_formula(self):
        weightFractionAl = 0.21358626371988801
        weightFractionNa = 0.27298103136883051
        weightFractionB = 0.51343270491128157

        comp = composition_from_formula('Al2Na3B12')
        self.assertAlmostEqual(weightFractionAl, comp[13], 4)
        self.assertAlmostEqual(weightFractionNa, comp[11], 4)
        self.assertAlmostEqual(weightFractionB, comp[5], 4)

        comp = composition_from_formula('Al 2 Na 3 B 12')
        self.assertAlmostEqual(weightFractionAl, comp[13], 4)
        self.assertAlmostEqual(weightFractionNa, comp[11], 4)
        self.assertAlmostEqual(weightFractionB, comp[5], 4)

        comp = composition_from_formula('Al2 Na3 B12')
        self.assertAlmostEqual(weightFractionAl, comp[13], 4)
        self.assertAlmostEqual(weightFractionNa, comp[11], 4)
        self.assertAlmostEqual(weightFractionB, comp[5], 4)

        self.assertRaises(ValueError, composition_from_formula, 'Aq2 Na3 B12')

        comp = composition_from_formula('Al2')
        self.assertAlmostEqual(1.0, comp[13], 4)

if __name__ == '__main__': #pragma: no cover
    logging.getLogger().setLevel(logging.DEBUG)
    unittest.main()
