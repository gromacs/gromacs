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

QMMM is currently not supported in GROMACS. 
