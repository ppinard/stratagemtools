.. _compute_thickness:

Compute thickness
=================

This tutorial shows how to use the main functionality of STRATAGem, to compute
the mass thickness of a layer from experimentally measured k-ratios.
For this example, we used the calculated k-ratios from the previous tutorial
(:ref:`custom_standard`) as "experimental k-ratios".
We should therefore compute a Al2O3 layer with a thickness of 30 nm.

As usual, let's start by importing the important classes and constants.

.. literalinclude:: /../../tutorial/compute_thickness.py
   :lines: 6-9
   
The next step is to define our :class:`.Sample`. 
Note here that we do not specify any thickness for the Al2O3 layer.
No argument is given, which is equivalent to setting the thickness to ``None``.

.. literalinclude:: /../../tutorial/compute_thickness.py
   :lines: 11-12
   
From the previous tutorial (:ref:`custom_standard`), the following Al and O 
:math:`\text{K}\alpha` k-ratios were respectively calculated 0.034 and 0.058
using the Al2O3 standard.
We now use these values to define the experiments.
We do not need to know the k-ratio for the Si :math:`\text{K}\alpha`, but
we need to specify this element as not analyzed.

.. literalinclude:: /../../tutorial/compute_thickness.py
   :lines: 16-19

We then execute the method :meth:`compute <.Stratagem.compute>` from the
:class:`.Stratagem` interface to compute the unknown thickness.
The method returns a new :class:`.Sample` object with the calculated thickness.
Note that the new :class:`.Sample` object could also be retrieved from the 
method :meth:`get_sample <.Stratagem.get_sample>` or the property 
:attr:`sample <.Stratagem.sample>` after executing the 
:meth:`compute <.Stratagem.compute>` method.

.. literalinclude:: /../../tutorial/compute_thickness.py
   :lines: 21-28
   
Finally, we can print the calculated thickness.

.. literalinclude:: /../../tutorial/compute_thickness.py
   :lines: 30-31
   
Getting back on our feet, we get 30.02 nm!
