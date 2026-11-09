New and improved features
^^^^^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Multiple molecule types with identical SETTLE parameters are now supported
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

This can be useful when one, for instance, wants to apply position restraints
to a subset of water molecules.

Periodic SCC-DFTB, GFN1-xTB and GFN2-xTB are available for CP2K QM/MM
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

The CP2K QM/MM interface can now generate periodic SCC-DFTB, GFN1-xTB and
GFN2-xTB inputs using short-range regularized point-charge coupling with SPME
electrostatics. For GFN1-xTB and GFN2-xTB CP2K should be compiled with support 
for the tblite interface. 
The user can select Gaussian Expansion of Electrostatic Potential (GEEP) 
for MM charges in tight-binding methods through the
``qmmm-cp2k-dftb-electrostatic-coupling = gauss`` MDP option, 
which requires a CP2K 2027.1 or higher, whereas the default is 
short-range regularized point-charge coupling 
``qmmm-cp2k-dftb-electrostatic-coupling = point-charge``.
Regular DFT methods are still using GEEP (Gauss) for electrostatics.

CMAP interactions now support free energy perturbation
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Energy correction map (CMAP) torsion interactions can now be
perturbed between an A state and a B state using the
``bonded-lambdas`` free-energy lambda component. The A and B
state CMAP grids are linearly interpolated at each lambda
value, and the corresponding :math:`\partial H/\partial\lambda`
contribution is accumulated for BAR/TI analysis.

To use CMAP FEP, specify a B-state CMAP type in the ``[ cmap ]``
section of the topology. The type can be given as a 1-based index
or as an explicit name token (e.g. ``GLY``, ``PZQ``) that matches
the optional name field on the corresponding ``[ cmaptypes ]``
header line. When the A and B state CMAP types differ, the
interaction is treated as perturbed; when they are identical (or
no B-state is specified) the standard non-perturbed code path is
used.

Neural network potentials can run on Intel GPUs
"""""""""""""""""""""""""""""""""""""""""""""""

The NNPot interface now runs model inference on Intel GPUs when a LibTorch
with XPU support (PyTorch 2.5 or later) is used. The model runs on the first
Intel GPU LibTorch enumerates, by default and with ``GMX_NN_DEVICE=gpu``.
Configuring warns when a SYCL build for Intel GPUs is paired with a LibTorch
without XPU support.

:issue:`5754`

Reading and writing of mmCIF/PDBx structure and trajectory files
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

|Gromacs| tools can now read and write structures in the mmCIF/PDBx format
(``.cif``), the primary distribution format of the wwPDB. This covers entries
that the legacy PDB format cannot express, including extended entry identifiers
(the ``pdb_########`` form, up to 12 characters) and chemical-component
identifiers of up to 5 characters. Multi-model ``.cif`` files can be read and
written as trajectories, for example with :ref:`gmx trjconv`.

With the ``-conect`` option of :ref:`gmx editconf` and :ref:`gmx trjconv`, atom
connectivity is written as ``struct_conn`` records, the mmCIF counterpart of the
PDB ``LINK`` and ``SSBOND`` records, and read back again.

mmCIF allows chain identifiers of more than one character, while |Gromacs| stores
a single character. Such identifiers are truncated to their first character, and
a warning naming them is printed.

:issue:`4753`
New I/O format for neural network potentials
""""""""""""""""""""""""""""""""""""""""""""

The I/O format for neural network potentials in the NNPot interface has been changed to a 
dictionary-based format. This allows for more flexibility in the input and output of the NNP model,
and will make it easier to add new features in the future. In particular, the ``model-inputX``
:ref:`mdp` options have been removed, and instead all available inputs are always passed to
the model as a dictionary. Similarly, the model is now expected to return a dictionary containing
the energy and forces, instead of a tuple. Note that this change is not backward-compatible, and 
models exported for use with past GROMACS versions will no longer work. For more details on how to
use the new I/O format, please see the :ref:`NNPot section in the reference manual<nnpot>`.

Stochastic dynamics simulations with GPU update
"""""""""""""""""""""""""""""""""""""""""""""""

Simulations using the stochastic dynamics (langevin) integrator
can now be run with the update step done on the GPU for the CUDA,
HIP and SYCL backends. This feature currently needs to be enabled
by explicitly setting ``-update gpu`` on the command line. 
