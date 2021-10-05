.. _qmmm:

Hybrid Quantum-Classical simulations (QM/MM) with CP2K interface
----------------------------------------------------------------

In a molecular mechanics (MM) force field, the influence of electrons is
expressed by empirical parameters that are assigned on the basis of
experimental data, or on the basis of results from high-level quantum
chemistry calculations. These are valid for the ground state of a given
covalent structure, and the MM approximation is usually sufficiently
accurate for ground-state processes in which the overall connectivity
between the atoms in the system remains unchanged. However, for
processes in which the connectivity does change, such as chemical
reactions, or processes that involve multiple electronic states, such as
photochemical conversions, electrons can no longer be ignored, and a
quantum mechanical description is required for at least those parts of
the system in which the reaction takes place.

One approach to the simulation of chemical reactions in solution, or in
enzymes, is to use a combination of quantum mechanics (QM) and molecular
mechanics (MM). The reacting parts of the system are treated quantum
mechanically, with the remainder being modeled using the force field.
The current version of |Gromacs| provides an interface to the popular
Quantum Chemistry package CP2K :ref:`188 <refcp2k2020>`.

Overview
^^^^^^^^

|Gromacs| interactions between the QM and the MM subsystems are handled using
the GEEP approach as described by Laino et al. :ref:`189 <refLaino2005>`. 
This method of evaluating interactions between the QM and MM subsystems 
is a variant of the "electrostatic embedding" scheme. The electrostatic 
interactions between the electrons of the QM region and the MM atoms 
and between the QM nuclei and the MM atoms are explicitly included into
the Hamiltonian for the QM subsystem:

   .. math::

      H^{QM/MM} =
      H^{QM}_e-\sum_i^n\sum_J^M\frac{e^2Q_J}{4\pi\epsilon_0r_{iJ}}+\sum_A^N\sum_J^M\frac{e^2Z_AQ_J}{e\pi\epsilon_0R_{AJ}},

where :math:`n` and :math:`N` are the number of electrons and nuclei
in the QM region, respectively, and :math:`M` is the number of
charged MM atoms. The first term on the right hand side is the
original electronic Hamiltonian of an isolated QM system. The first
of the double sums is the total electrostatic interaction between the
QM electrons and the MM atoms. The total electrostatic interaction of
the QM nuclei with the MM atoms is given by the second double sum. 
An important advantage of using the CP2K/GEEP combination is that it allows
evaluation of forces for both QM-QM and QM-MM interactions,
in the case of systems with periodic boundary conditions (PBC).
To avoid double accounting for electrostatic interactions and LJ,
classical MM charge on the QM atoms are zeroed out as well as LJ
interactions between QM-QM atoms are excluded. It should be noted that 
LJ interactions between QM-MM atoms are kept and still calculated by |Gromacs|.
Bonded interactions between QM and MM atoms are described at the MM
level by the appropriate force-field terms. All bonds,
consisting of 2 QM atoms, angles and settles containing 2 or 3 QM atoms, 
dihedrals containing 3 or 4 QM atoms are excluded from the forcefield
evaluation. Broken chemical bonds between QM and MM subsystems needs to be capped
in the QM calculation. This is done within CP2K by adding a hydrogen atom to 
complete the valence of the QM region. The force on this atom, which is present 
in the QM region only, is distributed over the two atoms of the bond. 
The cap atom is usually referred to as a link atom. Within the interface 
all described topology modifications are performed automatically during :ref:`gmx grompp` pre-processing.

Software prerequisites
^^^^^^^^^^^^^^^^^^^^^^

CP2K version 8.1 (or later) should be linked into |Gromacs| as libcp2k.
For a specific installation instructions please follow the :ref:`installing with CP2K` guide.

Limitations in simulation techniques
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The QM/MM interface limits simulations in two ways.
First, no topology modifications are possible during the simulations in the QM region.
Second, interface completely ignores "B" state parameters in the topology, making
double topology setups impossible, e.g. free-energy perturbation simulations (:ref:`dgimplement`).

In addition it should be noted that the contribution of forces from QM/MM to the system 
virial are not accounted for. The size of the effect on the pressure-coupling algorithm 
grows with the total summed force due to QM-MM interactions and might produce artifacts 
in simulations with the NPT ensemble.

Usage
^^^^^

QM/MM simulations with CP2K interface are controlled by setting :ref:`mdp` file options and,
in some cases, providing an additional input file for :ref:`gmx grompp` with the ``-qmi``
command-line option. All options that are related to QM/MM simulations with CP2K 
are prefixed with ``qmmm-cp2k``.

Setting :mdp-value:`qmmm-cp2k-active=true` will trigger a QM/MM simulation using the whole
system as QM part and default parameters for all other options.

Choosing atoms for QM calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The QM part of your system is chosen with a name that corresponds to an atom group
in the index file of |Gromacs| to the :mdp:`qmmm-cp2k-qmgroup` option in :ref:`mdp` file.
The typical QM part should consist of atoms that are interesting from the chemical point of view,
i.e. part of the system where reaction happens. To make computation of the
QM part feasible, it should be small and
as compact as possible in a space. DFT simulations often scale as 3rd order of 
the number of atoms in the QM part. This means increasing number of atoms in the QM part 
by a factor of 2 will slow down the simulation by a factor of 8.

In addition user should provide total charge of your QM subsystem with
:mdp:`qmmm-cp2k-qmcharge` option and spin-state (multiplicity) with :mdp:`qmmm-cp2k-qmmultiplicity`
option.

Supported QM methods
^^^^^^^^^^^^^^^^^^^^

The QM method is chosen with :mdp:`qmmm-cp2k-qmmethod` in the :ref:`mdp` file.
Currently the following QM methods are supported:

#. :mdp-value:`qmmm-cp2k-qmmethod=PBE` - DFT using PBE functional and DZVP-MOLOPT basis set.
#. :mdp-value:`qmmm-cp2k-qmmethod=BLYP` - DFT using BLYP functional and DZVP-MOLOPT basis set.

That list will be updated with a new methods once they are tested and included into the
interface.

Providing your own CP2K input file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In addition it is possible to use custom external CP2K input file with 
:mdp-value:`qmmm-cp2k-qmmethod=INPUT` and providing file with 
:ref:`gmx grompp` with ``-qmi`` option. The external file will be incorporated into the
:ref:`tpr` file of the simulation and are subject to the following restrictions:

#. ``RUN_TYPE`` option in the CP2K input should be equal to ``ENERGY_FORCE``.
#. ``CHARGE`` option should be present.
#. ``MULTIPILICTY`` option should be present.
#. ``COORD_FILE_NAME`` option should be present pointing towards :ref:`pdb` file. 
#. Both ``CHARGE_EXTENDED TRUE`` and ``COORD_FILE_FORMAT PDB`` options should be present.
#. Incremental includes (``@INCLUDE`` directive) are not allowed in the CP2K input file . 

Changing names of CP2K files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

During :ref:`gmx mdrun` simulation additional files will be produced with ``.inp``, ``.out`` and
``.pdb``. They contain CP2K input, CP2K output and :ref:`pdb` file with point charges of MM atoms
in the extended beta field. By default all CP2K related files names will be deduced from :ref:`tpr` 
simulation file name by adding ``_cp2k`` suffix. In order to change it manually 
:mdp:`qmmm-cp2k-qmfilenames` option should be used.

Output
^^^^^^

The energy output file will contain an additional "Quantum En." term.
This is the energy that is added to the system from the QM/MM interactions.
In addition, a file containing CP2K output will appear in the simulation directory 
with the ``.out`` extension.

Future developments
^^^^^^^^^^^^^^^^^^^

support of additional DFT methods will be added in the future, as well as semi-empirical and 
DFTB description of the QM subsystem will be allowed. Support of the multiple 
time-stepping approach to speed-up simulation will be added. Excited state simulations
will be implemented with TD-DFT description of the wavefunction.
