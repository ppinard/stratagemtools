Calculate k-ratios
==================

This is a tutorial to calculate the Al and O :math:`\text{K}\alpha` k-ratios 
from a multilayer sample consisting of a Si substrate and 30-nm Al2O3 layer at 
an accelerating voltage of 15 kV.

Let's start by importing the three classes that will later need: 
:class:`.Sample`, :class:`.Experiment` and :class:`.Stratagem` classes.

.. literalinclude:: /../../tutorial/calculate_kratio.py
   :lines: 6-9

The :class:`.Sample` class is used to define the composition of the substrate 
and the composition, thickness and mass thickness of each layer.
In this example, we first define the substrate composition as follows:

.. literalinclude:: /../../tutorial/calculate_kratio.py
   :lines: 10
   
The argument ``{14: 1.0}`` defines the substrate composition, consisting of
silicon (atomic number 14) and with a weight fraction of 1.0 (pure silicon).
To help us define the layer composition, we can use the utility function
:func:`.composition_from_formula` in the :mod:`.sample` module.
The function returns the composition expressed in weight fraction for a given 
chemical formula.

.. literalinclude:: /../../tutorial/calculate_kratio.py
   :lines: 12-13
   
The calculated composition can then be used to define a layer, using the method
:meth:`add_layer <.Sample.add_layer>`.
The thickness and density (taken from Wikipedia) are also specified.
Note that the thickness is expressed in meters and the density in kilograms per
cubic meter.

.. literalinclude:: /../../tutorial/calculate_kratio.py
   :lines: 14
   
The next step is to define :class:`.Experiment`'s.
:class:`.Experiment` specifies the experimental parameters used to analyze
or to use to calculate the k-ratio of each element.
As such, one experiment must be created for each element in the sample, 
whether or not it is analyzed or of interest to be calculated.
For this example, the experiments are:

.. literalinclude:: /../../tutorial/calculate_kratio.py
   :lines: 16-20

The constant :data:`.LINE_KA` is imported from the :mod:`.experiment` module to
specify the :math:`\text{K}\alpha` X-ray line.
For the moment, the standards (the denominator of the k-ratio) are all assumed 
to be pure sample, i.e. pure silicon, pure aluminum and (yes!) pure oxygen.
The use of custom standards is addressed in another tutorial, 
:ref:`custom_standard`.

The :class:`Sample` and :class:`Experiment` objects should then be added to the
:class:`.Stratagem` interface.
The interface works as a context manager (``with`` statement) in order to 
establish and properly close the connection to the STRATAGem's DLL.
All operations on a sample and experiments should be performed inside the 
``with`` statement.
The following lines of code set the sample, add the experiment and compute the
k-ratios.

.. literalinclude:: /../../tutorial/calculate_kratio.py
   :lines: 22-25

The :meth:`compute_kratios <.Stratagem.compute_kratios>` method returns a 
:class:`dict` where the keys are the experiments and the values, k-ratios.
To help printing the results, the utility function 
:meth:`symbol <.element_properties.symbol>` can be used to convert atomic number
into element symbol.

.. literalinclude:: /../../tutorial/calculate_kratio.py
   :lines: 27-29
   
Going back to the ``with`` statement, it is advisable to always get the 
geometry, type of :math:`\phi(\rho z)` and fluorescence flag, as the default
values may be changed.
*stratagemtools* relies on the default values from STRATAGem.

.. literalinclude:: /../../tutorial/calculate_kratio.py
   :lines: 31-36



