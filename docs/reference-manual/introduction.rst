Introduction
============

.. _compchem:

Computational Chemistry and Molecular Modeling
----------------------------------------------

|Gromacs| is an engine to perform molecular dynamics simulations and
energy minimization. These are two of the many techniques that belong to
the realm of computational chemistry and molecular modeling.
*Computational chemistry* is just a name to indicate the use of
computational techniques in chemistry, ranging from quantum mechanics of
molecules to dynamics of large complex molecular aggregates. *Molecular
modeling* indicates the general process of describing complex chemical
systems in terms of a realistic atomic model, with the goal being to
understand and predict macroscopic properties based on detailed
knowledge on an atomic scale. Often, molecular modeling is used to
design new materials, for which the accurate prediction of physical
properties of realistic systems is required.

Macroscopic physical properties can be distinguished by

#. *static equilibrium properties*, such as the binding constant of an
   inhibitor to an enzyme, the average potential energy of a system, or the
   radial distribution function of a liquid, and 

#. *dynamic or non-equilibrium properties*, such as the viscosity of a liquid,
   diffusion processes in membranes, the dynamics of phase changes,
   reaction kinetics, or the dynamics of defects in crystals. 

The choice of
technique depends on the question asked and on the feasibility of the
method to yield reliable results at the present state of the art.
Ideally, the (relativistic) time-dependent Schrödinger equation
describes the properties of molecular systems with high accuracy, but
anything more complex than the equilibrium state of a few atoms cannot
be handled at this *ab initio* level. Thus, approximations are
necessary; the higher the complexity of a system and the longer the time
span of the processes of interest is, the more severe the required
approximations are. At a certain point (reached very much earlier than
one would wish), the *ab initio* approach must be augmented or replaced
by *empirical* parameterization of the model used. Where simulations
based on physical principles of atomic interactions still fail due to
the complexity of the system, molecular modeling is based entirely on a
similarity analysis of known structural and chemical data. The QSAR
methods (Quantitative Structure-Activity Relations) and many
homology-based protein structure predictions belong to the latter
category.

Macroscopic properties are always ensemble averages over a
representative statistical ensemble (either equilibrium or
non-equilibrium) of molecular systems. For molecular modeling, this has
two important consequences:

-  The knowledge of a single structure, even if it is the structure of
   the global energy minimum, is not sufficient. It is necessary to
   generate a representative ensemble at a given temperature, in order
   to compute macroscopic properties. But this is not enough to compute
   thermodynamic equilibrium properties that are based on free energies,
   such as phase equilibria, binding constants, solubilities, relative
   stability of molecular conformations, etc. The computation of free
   energies and thermodynamic potentials requires special extensions of
   molecular simulation techniques.

-  While molecular simulations, in principle, provide atomic details of
   the structures and motions, such details are often not relevant for
   the macroscopic properties of interest. This opens the way to
   simplify the description of interactions and average over irrelevant
   details. The science of statistical mechanics provides the
   theoretical framework for such simplifications. There is a hierarchy
   of methods ranging from considering groups of atoms as one unit,
   describing motion in a reduced number of collective coordinates,
   averaging over solvent molecules with potentials of mean force
   combined with stochastic dynamics :ref:`9 <refGunsteren90>`, to
   *mesoscopic dynamics* describing densities rather than atoms and
   fluxes as response to thermodynamic gradients rather than velocities
   or accelerations as response to forces \ :ref:`10 <refFraaije93>`.

For the generation of a representative equilibrium ensemble two methods
are available: 

#. *Monte Carlo simulations* and

#. *Molecular Dynamics simulations*. 
  
For the generation of non-equilibrium
ensembles and for the analysis of dynamic events, only the second method
is appropriate. While Monte Carlo simulations are more simple than MD
(they do not require the computation of forces), they do not yield
significantly better statistics than MD in a given amount of computer
time. Therefore, MD is the more universal technique. If a starting
configuration is very far from equilibrium, the forces may be
excessively large and the MD simulation may fail. In those cases, a
robust *energy minimization* is required. Another reason to perform an
energy minimization is the removal of all kinetic energy from the
system: if several “snapshots” from dynamic simulations must be
compared, energy minimization reduces the thermal noise in the
structures and potential energies so that they can be compared better.

Molecular Dynamics Simulations
------------------------------

MD simulations solve Newton’s equations of motion for a system of
:math:`N` interacting atoms:

.. math:: m_i \frac{\partial^2 \mathbf{r}_i}{\partial t^2}  = \mathbf{F}_i, \;i=1 \ldots N.
          :label: eqnnewtonslaws

The forces are the negative derivatives of a potential function
:math:`V(\mathbf{r}_1, \mathbf{r}_2, \ldots, \mathbf{r}_N)`:

.. math:: \mathbf{F}_i = - \frac{\partial V}{\partial \mathbf{r}_i}
          :label: eqnmdforces

The equations are solved simultaneously in small time steps. The system
is followed for some time, taking care that the temperature and pressure
remain at the required values, and the coordinates are written to an
output file at regular intervals. The coordinates as a function of time
represent a *trajectory* of the system. After initial changes, the
system will usually reach an *equilibrium state*. By averaging over an
equilibrium trajectory, many macroscopic properties can be extracted
from the output file.

It is useful at this point to consider the limitations of MD
simulations. The user should be aware of those limitations and always
perform checks on known experimental properties to assess the accuracy
of the simulation. We list the approximations below.

**The simulations are classical**

-     Using Newton’s equation of motion automatically implies the use of
      *classical mechanics* to describe the motion of atoms. This is all
      right for most atoms at normal temperatures, but there are
      exceptions. Hydrogen atoms are quite light and the motion of
      protons is sometimes of essential quantum mechanical character.
      For example, a proton may *tunnel* through a potential barrier in
      the course of a transfer over a hydrogen bond. Such processes
      cannot be properly treated by classical dynamics! Helium liquid at
      low temperature is another example where classical mechanics
      breaks down. While helium may not deeply concern us, the high
      frequency vibrations of covalent bonds should make us worry! The
      statistical mechanics of a classical harmonic oscillator differs
      appreciably from that of a real quantum oscillator when the
      resonance frequency :math:`\nu` approximates or exceeds
      :math:`k_BT/h`. Now at room temperature the wavenumber
      :math:`\sigma = 1/\lambda = \nu/c` at which :math:`h
      \nu = k_BT` is approximately 200 cm\ :math:`^{-1}`. Thus, all
      frequencies higher than, say, 100 cm\ :math:`^{-1}` may misbehave
      in classical simulations. This means that practically all bond and
      bond-angle vibrations are suspect, and even hydrogen-bonded
      motions as translational or librational H-bond vibrations are
      beyond the classical limit (see :numref:`Table %s <tab-vibrations>`)
      What can we do? 

.. |H2CX| replace:: H\ :math:`_2`\ CX
.. |OHO1| replace:: O-H\ :math:`\cdots`\ O
.. |INCM| replace:: :math:`\mathrm{cm}~^{-1}`

.. _tab-vibrations:

.. table::
        Typical vibrational frequencies (wavenumbers) in molecules and hydrogen-bonded
        liquids. Compare :math:`kT/h = 200~\mathrm{cm}^{-1}` at 300 K.
        :widths: auto
        :align: center

        +---------------+-------------+------------+
        |               | type of     | wavenumber |
        | type of bond  | vibration   | |INCM|     |
        +===============+=============+============+
        | C-H, O-H, N-H | stretch     | 3000--3500 |
        +---------------+-------------+------------+
        | C=C, C=O      | stretch     | 1700--2000 |
        +---------------+-------------+------------+
        | HOH           | bending     | 1600       |
        +---------------+-------------+------------+
        | C-C           | stretch     | 1400--1600 |
        +---------------+-------------+------------+
        | |H2CX|        | sciss, rock | 1000--1500 |
        +---------------+-------------+------------+
        | CCC           | bending     |  800--1000 |
        +---------------+-------------+------------+
        | |OHO1|        | libration   |  400--700  |
        +---------------+-------------+------------+
        | |OHO1|        | stretch     |   50--200  |
        +---------------+-------------+------------+



-     Well, apart from real quantum-dynamical simulations, we can do one
      of two things:

      (a)   If we perform MD simulations using harmonic oscillators for
            bonds, we should make corrections to the total internal energy
            :math:`U = E_{kin} + E_{pot}` and specific heat :math:`C_V` (and
            to entropy :math:`S` and free energy :math:`A` or :math:`G` if
            those are calculated). The corrections to the energy and specific
            heat of a one-dimensional oscillator with frequency :math:`\nu`
            are: \ :ref:`11 <refMcQuarrie76>`

            .. math:: U^{QM} = U^{cl} +kT \left( {\frac{1}{2}}x - 1 + \frac{x}{e^x-1} \right)
                      :label: eqnmdqmcorr

            .. math:: C_V^{QM} = C_V^{cl} + k \left( \frac{x^2e^x}{(e^x-1)^2} - 1 \right)
                      :label: eqnmdqmcorr2

            where :math:`x=h\nu /kT`. The classical oscillator absorbs too
            much energy (:math:`kT`), while the high-frequency quantum
            oscillator is in its ground state at the zero-point energy level
            of :math:`\frac{1}{2} h\nu`.

      (b)   We can treat the bonds (and bond angles) as
            *constraints* in the equations of
            motion. The rationale behind this is that a quantum oscillator in
            its ground state resembles a constrained bond more closely than a
            classical oscillator. A good practical reason for this choice is
            that the algorithm can use larger time steps when the highest
            frequencies are removed. In practice the time step can be made
            four times as large when bonds are constrained than when they are
            oscillators \ :ref:`12 <refGunsteren77>`. |Gromacs| has this
            option for the bonds and bond angles. The flexibility of the
            latter is rather essential to allow for the realistic motion and
            coverage of configurational space \ :ref:`13 <refGunsteren82>`.

**Electrons are in the ground state**
      In MD we use a *conservative* force field that is a function of
      the positions of atoms only. This means that the electronic
      motions are not considered: the electrons are supposed to adjust
      their dynamics instantly when the atomic positions change (the
      *Born-Oppenheimer*
      approximation), and remain in their ground state. This is really
      all right, almost always. But of course, electron transfer
      processes and electronically excited states can not be treated.
      Neither can chemical reactions be treated properly, but there are
      other reasons to shy away from reactions for the time being.

**Force fields are approximate**
      Force fields
      provide the forces.
      They are not really a part of the simulation method and their
      parameters can be modified by the user as the need arises or
      knowledge improves. But the form of the forces that can be used in
      a particular program is subject to limitations. The force field
      that is incorporated in |Gromacs| is described in Chapter 4. In the
      present version the force field is pair-additive (apart from
      long-range Coulomb forces), it cannot incorporate
      polarizabilities, and it does not contain fine-tuning of bonded
      interactions. This urges the inclusion of some limitations in this
      list below. For the rest it is quite useful and fairly reliable
      for biologically-relevant macromolecules in aqueous solution!

**The force field is pair-additive**
      This means that all *non-bonded* forces result from the sum of
      non-bonded pair interactions. Non pair-additive interactions, the
      most important example of which is interaction through atomic
      polarizability, are represented by *effective pair potentials*.
      Only average non pair-additive contributions are incorporated.
      This also means that the pair interactions are not pure, *i.e.*,
      they are not valid for isolated pairs or for situations that
      differ appreciably from the test systems on which the models were
      parameterized. In fact, the effective pair potentials are not that
      bad in practice. But the omission of polarizability also means
      that electrons in atoms do not provide a dielectric constant as
      they should. For example, real liquid alkanes have a dielectric
      constant of slightly more than 2, which reduce the long-range
      electrostatic interaction between (partial) charges. Thus, the
      simulations will exaggerate the long-range Coulomb terms. Luckily,
      the next item compensates this effect a bit.

**Long-range interactions are cut off**
      In this version, |Gromacs| always uses a
      cut-off
      radius for the Lennard-Jones
      interactions and sometimes for the Coulomb interactions as well.
      The “minimum-image convention” used by |Gromacs| requires that only
      one image of each particle in the periodic boundary conditions is
      considered for a pair interaction, so the cut-off radius cannot
      exceed half the box size. That is still pretty big for large
      systems, and trouble is only expected for systems containing
      charged particles. But then truly bad things can happen, like
      accumulation of charges at the cut-off boundary or very wrong
      energies! For such systems, you should consider using one of the
      implemented long-range electrostatic algorithms, such as
      particle-mesh Ewald \ :ref:`14 <refDarden93>`,
      :ref:`15 <refEssmann95>`.

**Boundary conditions are unnatural**
      Since system size is small (even 10,000 particles is small), a
      cluster of particles will have a lot of unwanted boundary with its
      environment (vacuum). We must avoid this condition if we wish to
      simulate a bulk system. As such, we use periodic boundary
      conditions to avoid real phase boundaries. Since liquids are not
      crystals, something unnatural remains. This item is mentioned last
      because it is the least of the evils. For large systems, the
      errors are small, but for small systems with a lot of internal
      spatial correlation, the periodic boundaries may enhance internal
      correlation. In that case, beware of, and test, the influence of
      system size. This is especially important when using lattice sums
      for long-range electrostatics, since these are known to sometimes
      introduce extra ordering.

Energy Minimization and Search Methods
--------------------------------------

As mentioned in sec. :ref:`Compchem`, in many cases energy minimization
is required. |Gromacs| provides a number of methods for local energy
minimization, as detailed in sec. :ref:`EM`.

The potential energy function of a (macro)molecular system is a very
complex landscape (or *hypersurface*) in a large number of dimensions.
It has one deepest point, the *global minimum* and a very large number
of *local minima*, where all derivatives of the potential energy
function with respect to the coordinates are zero and all second
derivatives are non-negative. The matrix of second derivatives, which is
called the *Hessian matrix*, has non-negative eigenvalues; only the
collective coordinates that correspond to translation and rotation (for
an isolated molecule) have zero eigenvalues. In between the local minima
there are *saddle points*, where the Hessian matrix has only one
negative eigenvalue. These points are the mountain passes through which
the system can migrate from one local minimum to another.

Knowledge of all local minima, including the global one, and of all
saddle points would enable us to describe the relevant structures and
conformations and their free energies, as well as the dynamics of
structural transitions. Unfortunately, the dimensionality of the
configurational space and the number of local minima is so high that it
is impossible to sample the space at a sufficient number of points to
obtain a complete survey. In particular, no minimization method exists
that guarantees the determination of the global minimum in any practical
amount of time. Impractical methods exist, some much faster than
others \ :ref:`16 <refGeman84>`. However, given a starting configuration,
it is possible to find the *nearest local minimum*. “Nearest” in this
context does not always imply “nearest” in a geometrical sense (*i.e.*,
the least sum of square coordinate differences), but means the minimum
that can be reached by systematically moving down the steepest local
gradient. Finding this nearest local minimum is all that |Gromacs| can do
for you, sorry! If you want to find other minima and hope to discover
the global minimum in the process, the best advice is to experiment with
temperature-coupled MD: run your system at a high temperature for a
while and then quench it slowly down to the required temperature; do
this repeatedly! If something as a melting or glass transition
temperature exists, it is wise to stay for some time slightly below that
temperature and cool down slowly according to some clever scheme, a
process called *simulated annealing*. Since no physical truth is
required, you can use your imagination to speed up this process. One
trick that often works is to make hydrogen atoms heavier (mass 10 or
so): although that will slow down the otherwise very rapid motions of
hydrogen atoms, it will hardly influence the slower motions in the
system, while enabling you to increase the time step by a factor of 3 or
4. You can also modify the potential energy function during the search
procedure, *e.g.* by removing barriers (remove dihedral angle functions
or replace repulsive potentials by *soft-core*
potentials \ :ref:`17 <refNilges88>`), but always take care to restore the correct
functions slowly. The best search method that allows rather drastic
structural changes is to allow excursions into four-dimensional
space \ :ref:`18 <refSchaik93>`, but this requires some extra programming
beyond the standard capabilities of |Gromacs|.

Three possible energy minimization methods are:

-  Those that require only function evaluations. Examples are the
   simplex method and its variants. A step is made on the basis of the
   results of previous evaluations. If derivative information is
   available, such methods are inferior to those that use this
   information.

-  Those that use derivative information. Since the partial derivatives
   of the potential energy with respect to all coordinates are known in
   MD programs (these are equal to minus the forces) this class of
   methods is very suitable as modification of MD programs.

-  Those that use second derivative information as well. These methods
   are superior in their convergence properties near the minimum: a
   quadratic potential function is minimized in one step! The problem is
   that for :math:`N` particles a :math:`3N\times 3N` matrix must be
   computed, stored, and inverted. Apart from the extra programming to
   obtain second derivatives, for most systems of interest this is
   beyond the available capacity. There are intermediate methods that
   build up the Hessian matrix on the fly, but they also suffer from
   excessive storage requirements. So |Gromacs| will shy away from this
   class of methods.

The *steepest descent* method, available in |Gromacs|, is of the second
class. It simply takes a step in the direction of the negative gradient
(hence in the direction of the force), without any consideration of the
history built up in previous steps. The step size is adjusted such that
the search is fast, but the motion is always downhill. This is a simple
and sturdy, but somewhat stupid, method: its convergence can be quite
slow, especially in the vicinity of the local minimum! The
faster-converging *conjugate gradient method* (see *e.g.*
:ref:`19 <refZimmerman91>`) uses gradient information from previous steps. In general,
steepest descents will bring you close to the nearest local minimum very
quickly, while conjugate gradients brings you *very* close to the local
minimum, but performs worse far away from the minimum. |Gromacs| also
supports the L-BFGS minimizer, which is mostly comparable to *conjugate
gradient method*, but in some cases converges faster.

.. raw:: latex

    \clearpage


