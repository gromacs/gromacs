Embedding proteins into the membranes
-------------------------------------

|Gromacs| is capable of inserting the protein into pre-equilibrated lipid
bilayers with minimal perturbation of the lipids using the method, which
was initially described as a ProtSqueeze technique, \ :ref:`157 <refYesylevskyy2007>`
and later implemented as g_membed tool \ :ref:`158 <refWolf2010>`. Currently the
functionality of g_membed is available in mdrun as described in the
user guide.

This method works by first artificially shrinking the protein in the
:math:`xy`-plane, then it removes lipids that overlap with that much
smaller core. Then the protein atoms are gradually resized back to their
initial configuration, using normal dynamics for the rest of the system,
so the lipids adapt to the protein. Further lipids are removed as
required.

.. raw:: latex

    \clearpage


