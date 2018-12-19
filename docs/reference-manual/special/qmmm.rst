.. _qmmm:

Mixed Quantum-Classical simulation techniques
---------------------------------------------

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
The current version of |Gromacs| provides interfaces to several popular
Quantum Chemistry packages (MOPAC :ref:`150 <refmopac>`,
GAMESS-UK \ :ref:`151 <refgamess-uk>`, Gaussian \ :ref:`152 <refg03>` and
CPMD \ :ref:`153 <refCar85a>`).

|Gromacs| interactions between the two subsystems are either handled as
described by Field et al. :ref:`154 <refField90a>` or within
the ONIOM approach by Morokuma and coworkers \ :ref:`155 <refMaseras96a>`,
:ref:`156 <refSvensson96a>`.

Overview
^^^^^^^^

Two approaches for describing the interactions between the QM and MM
subsystems are supported in this version:

#. **Electronic Embedding** The electrostatic interactions between the
   electrons of the QM region and the MM atoms and between the QM nuclei
   and the MM atoms are included in the Hamiltonian for the QM
   subsystem:

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
   Bonded interactions between QM and MM atoms are described at the MM
   level by the appropriate force-field terms. Chemical bonds that
   connect the two subsystems are capped by a hydrogen atom to complete
   the valence of the QM region. The force on this atom, which is
   present in the QM region only, is distributed over the two atoms of
   the bond. The cap atom is usually referred to as a link atom.

#. **ONIOM** In the ONIOM approach, the energy and gradients are first
   evaluated for the isolated QM subsystem at the desired level of *ab
   initio* theory. Subsequently, the energy and gradients of the total
   system, including the QM region, are computed using the molecular
   mechanics force field and added to the energy and gradients
   calculated for the isolated QM subsystem. Finally, in order to
   correct for counting the interactions inside the QM region twice, a
   molecular mechanics calculation is performed on the isolated QM
   subsystem and the energy and gradients are subtracted. This leads to
   the following expression for the total QM/MM energy (and gradients
   likewise):

   .. math::

      E_{tot} = E_{I}^{QM}
      +E_{I+II}^{MM}-E_{I}^{MM},

   where the subscripts I and II refer to the QM and MM subsystems,
   respectively. The superscripts indicate at what level of theory the
   energies are computed. The ONIOM scheme has the advantage that it is
   not restricted to a two-layer QM/MM description, but can easily
   handle more than two layers, with each layer described at a different
   level of theory.

Usage
^^^^^

To make use of the QM/MM functionality in |Gromacs|, one needs to:

#. introduce link atoms at the QM/MM boundary, if needed;

#. specify which atoms are to be treated at a QM level;

#. specify the QM level, basis set, type of QM/MM interface and so on.

Adding link atoms
^^^^^^^^^^^^^^^^^

At the bond that connects the QM and MM subsystems, a link atoms is
introduced. In |Gromacs| the link atom has special atomtype, called LA.
This atomtype is treated as a hydrogen atom in the QM calculation, and
as a virtual site in the force-field calculation. The link atoms, if
any, are part of the system, but have no interaction with any other
atom, except that the QM force working on it is distributed over the two
atoms of the bond. In the topology, the link atom (LA), therefore, is
defined as a virtual site atom:

::

    [ virtual_sites2 ]
    LA QMatom MMatom 1 0.65

See sec. :ref:`vsitetop` for more details on how virtual sites are
treated. The link atom is replaced at every step of the simulation.

In addition, the bond itself is replaced by a constraint:

::

    [ constraints ]
    QMatom MMatom 2 0.153

**Note** that, because in our system the QM/MM bond is a carbon-carbon
bond (0.153 nm), we use a constraint length of 0.153 nm, and dummy
position of 0.65. The latter is the ratio between the ideal C-H bond
length and the ideal C-C bond length. With this ratio, the link atom is
always 0.1 nm away from the ``QMatom``, consistent with the carbon-hydrogen
bond length. If the QM and MM subsystems are connected by a different
kind of bond, a different constraint and a different dummy position,
appropriate for that bond type, are required.

Specifying the QM atoms
^^^^^^^^^^^^^^^^^^^^^^^

Atoms that should be treated at a QM level of theory, including the link
atoms, are added to the index file. In addition, the chemical bonds
between the atoms in the QM region are to be defined as connect bonds
(bond type 5) in the topology file:

::

    [ bonds ]
    QMatom1 QMatom2 5
    QMatom2 QMatom3 5

Specifying the QM/MM simulation parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the :ref:`mdp` file, the following parameters control a
QM/MM simulation.

``QMMM = no``
    | If this is set to ``yes``, a QM/MM simulation is
      requested. Several groups of atoms can be described at different
      QM levels separately. These are specified in the QMMM-grps field
      separated by spaces. The level of *ab initio* theory at which the
      groups are described is specified by ``QMmethod`` and
      ``QMbasis`` Fields. Describing the groups at different
      levels of theory is only possible with the ONIOM QM/MM scheme,
      specified by ``QMMMscheme``.

``QMMM-grps =``
    | groups to be described at the QM level

``QMMMscheme = normal``
    | Options are ``normal`` and ``ONIOM``. This
      selects the QM/MM interface. ``normal`` implies that
      the QM subsystem is electronically embedded in the MM subsystem.
      There can only be one ``QMMM-grps`` that is modeled at
      the ``QMmethod`` and ``QMbasis`` level of
      * ab initio* theory. The rest of the system is described at the MM
      level. The QM and MM subsystems interact as follows: MM point
      charges are included in the QM one-electron Hamiltonian and all
      Lennard-Jones interactions are described at the MM level. If
      ``ONIOM`` is selected, the interaction between the
      subsystem is described using the ONIOM method by Morokuma and
      co-workers. There can be more than one QMMM-grps each modeled at a
      different level of QM theory (QMmethod and QMbasis).

``QMmethod =``
    | Method used to compute the energy and gradients on the QM atoms.
      Available methods are AM1, PM3, RHF, UHF, DFT, B3LYP, MP2, CASSCF,
      MMVB and CPMD. For CASSCF, the number of electrons and orbitals
      included in the active space is specified by
      ``CASelectrons`` and ``CASorbitals``. For
      CPMD, the plane-wave cut-off is specified by the
      ``planewavecutoff`` keyword.

``QMbasis =``
    | Gaussian basis set used to expand the electronic wave-function.
      Only Gaussian basis sets are currently available, i.e. STO-3G,
      3-21G, 3-21G\*, 3-21+G\*, 6-21G, 6-31G, 6-31G\*, 6-31+G\*, and
      6-311G. For CPMD, which uses plane wave expansion rather than
      atom-centered basis functions, the ``planewavecutoff``
      keyword controls the plane wave expansion.

``QMcharge =``
    | The total charge in *e* of the ``QMMM-grps``. In case
      there are more than one ``QMMM-grps``, the total
      charge of each ONIOM layer needs to be specified separately.

``QMmult =``
    | The multiplicity of the ``QMMM-grps``. In case there
      are more than one ``QMMM-grps``, the multiplicity of
      each ONIOM layer needs to be specified separately.

``CASorbitals =``
    | The number of orbitals to be included in the active space when
      doing a CASSCF computation.

``CASelectrons =``
    | The number of electrons to be included in the active space when
      doing a CASSCF computation.

``SH = no``
    | If this is set to yes, a QM/MM MD simulation on the excited
      state-potential energy surface and enforce a diabatic hop to the
      ground-state when the system hits the conical intersection
      hyperline in the course the simulation. This option only works in
      combination with the CASSCF method.

Output
^^^^^^

The energies and gradients computed in the QM calculation are added to
those computed by |Gromacs|. In the :ref:`edr` file there is a
section for the total QM energy.

Future developments
^^^^^^^^^^^^^^^^^^^

Several features are currently under development to increase the
accuracy of the QM/MM interface. One useful feature is the use of
delocalized MM charges in the QM computations. The most important
benefit of using such smeared-out charges is that the Coulombic
potential has a finite value at interatomic distances. In the point
charge representation, the partially-charged MM atoms close to the QM
region tend to “over-polarize” the QM system, which leads to artifacts
in the calculation.

What is needed as well is a transition state optimizer.
