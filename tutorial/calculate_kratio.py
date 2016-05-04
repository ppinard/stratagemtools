"""
This is a tutorial to calculate the Al, O and Si k-ratios from a multilayer
sample consisting of a Si substrate and 30-nm Al2O3 layer.
"""

from stratagemtools.sample import Sample
from stratagemtools.experiment import Experiment
from stratagemtools.stratagem import Stratagem

sample = Sample({14: 1.0})

from stratagemtools.sample import composition_from_formula
comp = composition_from_formula('Al2O3')
sample.add_layer(comp, 30e-9, density_kg_m3=3950.0)

from stratagemtools.experiment import LINE_KA
energy_eV = 15e3
exp_si = Experiment(14, LINE_KA, energy_eV)
exp_al = Experiment(13, LINE_KA, energy_eV)
exp_o = Experiment(8, LINE_KA, energy_eV)

with Stratagem() as strata:
    strata.set_sample(sample)
    strata.add_experiments(exp_si, exp_al, exp_o)
    kratios = strata.compute_kratios()

import stratagemtools.element_properties as ep
for exp, kratio  in kratios.items():
    print('{0}: {1:.3f}'.format(ep.symbol(exp.z), kratio))

import math
from stratagemtools.stratagem import PRZMODE_XPP, FLUORESCENCE_LINE_CONT
with Stratagem() as strata:
    strata.set_geometry(math.radians(40), 0.0, 0.0)
    strata.set_prz_mode(PRZMODE_XPP)
    strata.set_fluorescence(FLUORESCENCE_LINE_CONT)
