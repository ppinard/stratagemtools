"""
This tutorial shows how to to compute the mass thickness of a layer from 
experimentally measured k-ratios.
"""

import math
from stratagemtools.sample import Sample, composition_from_formula
from stratagemtools.experiment import Experiment, LINE_KA
from stratagemtools.stratagem import Stratagem, PRZMODE_XPP, FLUORESCENCE_LINE_CONT

unknown = Sample({14: 1.0})
unknown.add_layer(composition_from_formula('Al2O3'), density_kg_m3=3950.0)

std_al2o3 = Sample(composition_from_formula('Al2O3'), density_kg_m3=3950.0)

energy_eV = 15e3
exp_si = Experiment(14, LINE_KA, energy_eV, analyzed=False)
exp_al = Experiment(13, LINE_KA, energy_eV, kratio=0.034, standard=std_al2o3)
exp_o = Experiment(8, LINE_KA, energy_eV, kratio=0.058, standard=std_al2o3)

with Stratagem() as strata:
    strata.set_geometry(math.radians(40), 0.0, 0.0)
    strata.set_prz_mode(PRZMODE_XPP)
    strata.set_fluorescence(FLUORESCENCE_LINE_CONT)

    strata.set_sample(unknown)
    strata.add_experiments(exp_si, exp_al, exp_o)
    newsample = strata.compute()

thickness_m = newsample.get_layer(0).thickness_m
print('Thickness of first layer: {0:.2f} nm'.format(thickness_m * 1e9))
