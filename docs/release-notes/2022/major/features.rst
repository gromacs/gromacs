New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Hybrid Quantum-Classical simulations (QM/MM) with CP2K interface
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Simulations of chemical reactions pathways can provide an atomistic insight into many 
biological and chemical processes. To perform such kind of modelling in complex systems, 
that includes solvent and/or proteins Multi-scale Quantum Mechanics / Molecular Mechanics 
(QM/MM) approaches are often used. Here we introduce a whole new interface to perform QM/MM 
simulations in fully periodic systems using MDModule that couples |Gromacs| with CP2K 
quantum chemistry package. This enables hybrid simulations of systems in systems 
where chemical reactions occurs. The interface supports most of the simulations techniques 
available in |Gromacs| including energy minimization, classical MD and enhanced sampling methods
such as umbrella sampling and accelerated weight histogram method.

Transformation pull coordinate for mathematical transformations of pull coordinates
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

A new pull coordinate type named transformation has been added. This enables mathematical
transformation of previously defined pull coordinates using a user supplied formula
in a string. This allows for example non-linear transformation of a distance, e.g.
a contact coordinate or (non-)linear combinations of multiple pull coordinates.
This is a powerful tool for defining complex reaction coordinates and it can be combined
with the Accelerated Weight Histogram Method to enhance sampling.

Replica-exchange molecular dynamics simulations with GPU update
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Replica-exchange molecular dynamics now works with GPU update.

A new formulation of soft-core interactions for free energy calculations
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

With this addition Gromacs allows to choose from two schemes to soften
non-bonded interactions during alchemical perturbations:
Beutler *et al.*\ :ref:`100 <refBeutler94>` and Gapsys *et al.*\ :ref:`185 <refGapsys2012>` soft-core functions.

More flexible sharing of biases in AWH
""""""""""""""""""""""""""""""""""""""

With the accelerated weight histogram method, biases can now be shared between
subsets of all simulations, without restrictions. The allows for more flexible
ensemble simulation setups, as well as simpler launches of sets of simulations.

More features implemented in modular simulator
""""""""""""""""""""""""""""""""""""""""""""""

Several features were added to the modular simulator, including all temperature
and pressure coupling algorithms available in the legacy simulator, expanded
ensemble and pull.

Free energy calculations now support all non-perturbed bonded interactions
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Previously |Gromacs| did not permit any usage of a few more special bonded
interactions (restricted angles/dihedrals or combined bending-torsion potentials)
in free energy calculations. These are now allowed, as long as the interaction
itself is not perturbed.

:issue:`3691`
