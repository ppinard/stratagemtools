.. _custom_standard:

Custom standard
===============

This tutorial shows how to use different standards in the definition of the
:class:`Experiment`.
The possibly to easily define standards is certainly an interesting feature of 
*stratagemtools*, in comparison to the graphical interface of STRATAGem where
the standards must be manually defined and saved in separate files.

In *stratagemtools*, a standard is a :class:`.Sample`.
Any sample definition can be used as a standard, as long as the composition
and thickness of every layer are known.
This example shows how to use three different standards to calculate the Al, O 
and Si :math:`\text{K}\alpha` k-ratios from the previous tutorial, 
:ref:`calculate_kratio`.

We import the same packages and constants as the last tutorial

.. literalinclude:: /../../tutorial/custom_standard.py
   :lines: 6-10
   
and create the unknown sample, a 30-nm Al2O3 layer over a Si substrate

 .. literalinclude:: /../../tutorial/custom_standard.py
   :lines: 12-14
   
Now we define three standards, Fe2O3, SiO2 and Al2O3:

 .. literalinclude:: /../../tutorial/custom_standard.py
   :lines: 16-18
   
We then create experiments for Si and Al using the new standards. 
Note that in the previous tutorial pure standards were assumed.

 .. literalinclude:: /../../tutorial/custom_standard.py
   :lines: 20-22
   
For the O :math:`\text{K}\alpha`, all three standards can be used.
To illustrate how easy it is to change the standard, we will loop over the
standards and define a new experiment using each standard.
This gives

 .. literalinclude:: /../../tutorial/custom_standard.py
   :lines: 24-42


