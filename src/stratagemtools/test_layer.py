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

# Third party modules.

# Local modules.
from layer import Layer

# Globals and constants variables.

class TestLayer(unittest.TestCase):

    def setUp(self):
        unittest.TestCase.setUp(self)

        elements = {29: 0.5, 30: 0.4}
        self.l1 = Layer(elements, thickness=5.0)
        self.l2 = Layer(elements, mass_thickness=3.668)
        self.l3 = Layer(elements)

    def tearDown(self):
        unittest.TestCase.tearDown(self)

    def testskeleton(self):
        self.assertEqual(2, len(list(self.l1.iter_elements())))
        self.assertAlmostEqual(5.0, self.l1.thickness, 3)
        self.assertAlmostEqual(3.668, self.l1.mass_thickness, 3)
        self.assertAlmostEqual(7.336, self.l1.density, 3)
        self.assertTrue(self.l1.is_thickness_known())

        self.assertEqual(2, len(list(self.l2.iter_elements())))
        self.assertAlmostEqual(5.0, self.l2.thickness, 3)
        self.assertAlmostEqual(3.668, self.l2.mass_thickness, 3)
        self.assertAlmostEqual(7.336, self.l2.density, 3)
        self.assertTrue(self.l2.is_thickness_known())

        self.assertEqual(2, len(list(self.l3.iter_elements())))
        self.assertAlmostEqual(-1.0, self.l3.thickness, 3)
        self.assertAlmostEqual(-1.0, self.l3.mass_thickness, 3)
        self.assertAlmostEqual(7.336, self.l3.density, 3)
        self.assertFalse(self.l3.is_thickness_known())

if __name__ == '__main__': #pragma: no cover
    logging.getLogger().setLevel(logging.DEBUG)
    unittest.main()
