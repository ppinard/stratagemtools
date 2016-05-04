"""
This tutorial shows how to use different standards in the definition of the
:class:`Experiment`.
"""

import math
from stratagemtools.sample import Sample, composition_from_formula
from stratagemtools.experiment import Experiment, LINE_KA
from stratagemtools.stratagem import Stratagem, PRZMODE_XPP, FLUORESCENCE_LINE_CONT
import stratagemtools.element_properties as ep

unknown = Sample({14: 1.0})
comp = composition_from_formula('Al2O3')
unknown.add_layer(comp, 30e-9, density_kg_m3=3950.0)

std_fe2o3 = Sample(composition_from_formula('Fe2O3'), density_kg_m3=5240.0)
std_sio2 = Sample(composition_from_formula('SiO2'), density_kg_m3=2650.0)
std_al2o3 = Sample(composition_from_formula('Al2O3'), density_kg_m3=3950.0)

energy_eV = 15e3
exp_si = Experiment(14, LINE_KA, energy_eV, standard=std_sio2)
exp_al = Experiment(13, LINE_KA, energy_eV, standard=std_al2o3)

with Stratagem() as strata:
    strata.set_geometry(math.radians(40), 0.0, 0.0)
    strata.set_prz_mode(PRZMODE_XPP)
    strata.set_fluorescence(FLUORESCENCE_LINE_CONT)

    for std_name, std_o in [('Fe2O3', std_fe2o3),
                            ('SiO2', std_sio2),
                            ('Al2O3', std_al2o3)]:
        exp_o = Experiment(8, LINE_KA, energy_eV, standard=std_o)

        strata.reset()
        strata.set_sample(unknown)
        strata.add_experiments(exp_si, exp_al, exp_o)
        kratios = strata.compute_kratios()

        print(std_name)
        for exp, kratio  in kratios.items():
            print('{0}: {1:.3f}'.format(ep.symbol(exp.z), kratio))
        print('-' * 80)
