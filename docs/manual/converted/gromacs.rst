|image|

| 

| Contributions from

 Emile Apol, Rossen Apostolov, Herman J.C. Berendsen,
Aldert van Buuren, Pär Bjelkmar, Rudi van Drunen,
Anton Feenstra, Sebastian Fritsch, Gerrit Groenhof,
Christoph Junghans, Jochen Hub, Peter Kasson,
Carsten Kutzner, Brad Lambeth, Per Larsson, Justin A. Lemkul,
Viveca Lindahl, Magnus Lundborg, Erik Marklund, Pieter Meulenhoff,
Teemu Murtola, Szilárd Páll, Sander Pronk,
Roland Schulz, Michael Shirts, Alfons Sijbers,
Peter Tieleman, Christian Wennberg and Maarten Wolf.

 Mark Abraham, Berk Hess, David van der Spoel, and Erik Lindahl.

| © 1991–2000: Department of Biophysical Chemistry, University of
  Groningen.
| Nijenborgh 4, 9747 AG Groningen, The Netherlands.

| © 2001–: The GROMACS development teams at the Royal Institute of
  Technology and
| Uppsala University, Sweden.

More information can be found on our website:
`www.gromacs.org <http://www.gromacs.org>`__.

Preface & Disclaimer
--------------------

This manual is not complete and has no pretention to be so due to lack
of time of the contributors – our first priority is to improve the
software. It is worked on continuously, which in some cases might mean
the information is not entirely correct.

Comments on form and content are welcome, please send them to one of the
mailing lists (see `www.gromacs.org <http://www.gromacs.org>`__), or
open an issue at `redmine.gromacs.org <http://redmine.gromacs.org>`__.
Corrections can also be made in the GROMACS git source repository and
uploaded to `gerrit.gromacs.org <http://gerrit.gromacs.org>`__.

We release an updated version of the manual whenever we release a new
version of the software, so in general it is a good idea to use a manual
with the same major and minor release number as your GROMACS
installation.

On-line Resources
-----------------

You can find more documentation and other material at our homepage
`www.gromacs.org <http://www.gromacs.org>`__. Among other things there
is an on-line reference, several GROMACS mailing lists with archives and
contributed topologies/force fields.

Citation information
--------------------

When citing this document in any scientific publication please refer to
it as:

    M.J. Abraham, D. van der Spoel, E. Lindahl, B. Hess, and the GROMACS
    development team,

    ,

    ()

However, we prefer that you cite (some of) the GROMACS
papers \ `1 <#ref-Bekker93a>`__\ `8 <#ref-Abraham2015>`__ when you
publish your results. Any future development depends on academic
research grants, since the package is distributed as free software!

GROMACS is *Free Software*
--------------------------

The entire GROMACS package is available under the GNU Lesser General
Public License (LGPL), version 2.1. This means it’s free as in free
speech, not just that you can use it without paying us money. You can
redistribute GROMACS and/or modify it under the terms of the LGPL as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version. For details, check the
COPYING file in the source code or consult
http://www.gnu.org/licenses/old-licenses/lgpl-2.1.html.

The GROMACS source code and and selected set of binary packages are
available on our homepage, `www.gromacs.org <http://www.gromacs.org>`__.
Have fun.

Introduction
============

Computational Chemistry and Molecular Modeling
----------------------------------------------

GROMACS is an engine to perform molecular dynamics simulations and
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

Macroscopic physical properties can be distinguished by (:math:`a`)
*static equilibrium properties*, such as the binding constant of an
inhibitor to an enzyme, the average potential energy of a system, or the
radial distribution function of a liquid, and (:math:`b`) *dynamic or
non-equilibrium properties*, such as the viscosity of a liquid,
diffusion processes in membranes, the dynamics of phase changes,
reaction kinetics, or the dynamics of defects in crystals. The choice of
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
   combined with stochastic dynamics `9 <#ref-Gunsteren90>`__, to
   *mesoscopic dynamics* describing densities rather than atoms and
   fluxes as response to thermodynamic gradients rather than velocities
   or accelerations as response to forces \ `10 <#ref-Fraaije93>`__.

For the generation of a representative equilibrium ensemble two methods
are available: (:math:`a`) *Monte Carlo simulations* and (:math:`b`)
*Molecular Dynamics simulations*. For the generation of non-equilibrium
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

.. math:: m_i \frac{\partial^2 {\mbox{\boldmath ${r}$}}_i}{\partial t^2}  = {\mbox{\boldmath ${F}$}}_i, \;i=1 \ldots N.

 The forces are the negative derivatives of a potential function
:math:`V({\mbox{\boldmath ${r}$}}_1, 
{\mbox{\boldmath ${r}$}}_2, \ldots, {\mbox{\boldmath ${r}$}}_N)`:

.. math:: {\mbox{\boldmath ${F}$}}_i = - \frac{\partial V}{\partial {\mbox{\boldmath ${r}$}}_i}

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
    | 
    | Using Newton’s equation of motion automatically implies the use of
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
      :math:`\sigma = 1/\lambda =
      \nu/c` at which :math:`h
      \nu = k_BT` is approximately 200 cm\ :math:`^{-1}`. Thus, all
      frequencies higher than, say, 100 cm\ :math:`^{-1}` may misbehave
      in classical simulations. This means that practically all bond and
      bond-angle vibrations are suspect, and even hydrogen-bonded
      motions as translational or librational H-bond vibrations are
      beyond the classical limit (see Table [tab:vibrations]). What can
      we do?

    | Well, apart from real quantum-dynamical simulations, we can do one
      of two things:
    | (a) If we perform MD simulations using harmonic oscillators for
      bonds, we should make corrections to the total internal energy
      :math:`U = E_{kin} + E_{pot}` and specific heat :math:`C_V` (and
      to entropy :math:`S` and free energy :math:`A` or :math:`G` if
      those are calculated). The corrections to the energy and specific
      heat of a one-dimensional oscillator with frequency :math:`\nu`
      are: \ `11 <#ref-McQuarrie76>`__

      .. math:: U^{QM} = U^{cl} +kT \left( {\frac{1}{2}}x - 1 + \frac{x}{e^x-1} \right)

      .. math:: C_V^{QM} = C_V^{cl} + k \left( \frac{x^2e^x}{(e^x-1)^2} - 1 \right),

       where :math:`x=h\nu /kT`. The classical oscillator absorbs too
      much energy (:math:`kT`), while the high-frequency quantum
      oscillator is in its ground state at the zero-point energy level
      of :math:`\frac{1}{2} h\nu`.
    | (b) We can treat the bonds (and bond angles) as *constraints* in
      the equations of motion. The rationale behind this is that a
      quantum oscillator in its ground state resembles a constrained
      bond more closely than a classical oscillator. A good practical
      reason for this choice is that the algorithm can use larger time
      steps when the highest frequencies are removed. In practice the
      time step can be made four times as large when bonds are
      constrained than when they are
      oscillators \ `12 <#ref-Gunsteren77>`__. GROMACS has this option
      for the bonds and bond angles. The flexibility of the latter is
      rather essential to allow for the realistic motion and coverage of
      configurational space \ `13 <#ref-Gunsteren82>`__.

**Electrons are in the ground state**
    | 
    | In MD we use a *conservative* force field that is a function of
      the positions of atoms only. This means that the electronic
      motions are not considered: the electrons are supposed to adjust
      their dynamics instantly when the atomic positions change (the
      *Born-Oppenheimer* approximation), and remain in their ground
      state. This is really all right, almost always. But of course,
      electron transfer processes and electronically excited states can
      not be treated. Neither can chemical reactions be treated
      properly, but there are other reasons to shy away from reactions
      for the time being.

**Force fields are approximate**
    | 
    | Force fields provide the forces. They are not really a part of the
      simulation method and their parameters can be modified by the user
      as the need arises or knowledge improves. But the form of the
      forces that can be used in a particular program is subject to
      limitations. The force field that is incorporated in GROMACS is
      described in Chapter 4. In the present version the force field is
      pair-additive (apart from long-range Coulomb forces), it cannot
      incorporate polarizabilities, and it does not contain fine-tuning
      of bonded interactions. This urges the inclusion of some
      limitations in this list below. For the rest it is quite useful
      and fairly reliable for biologically-relevant macromolecules in
      aqueous solution!

**The force field is pair-additive**
    | 
    | This means that all *non-bonded* forces result from the sum of
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
    | 
    | In this version, GROMACS always uses a cut-off radius for the
      Lennard-Jones interactions and sometimes for the Coulomb
      interactions as well. The “minimum-image convention” used by
      GROMACS requires that only one image of each particle in the
      periodic boundary conditions is considered for a pair interaction,
      so the cut-off radius cannot exceed half the box size. That is
      still pretty big for large systems, and trouble is only expected
      for systems containing charged particles. But then truly bad
      things can happen, like accumulation of charges at the cut-off
      boundary or very wrong energies! For such systems, you should
      consider using one of the implemented long-range electrostatic
      algorithms, such as particle-mesh Ewald \ `14 <#ref-Darden93>`__,
      `15 <#ref-Essmann95>`__.

**Boundary conditions are unnatural**
    | 
    | Since system size is small (even 10,000 particles is small), a
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

As mentioned in sec. [sec:Compchem], in many cases energy minimization
is required. GROMACS provides a number of methods for local energy
minimization, as detailed in sec. [sec:EM].

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
others \ `16 <#ref-Geman84>`__. However, given a starting configuration,
it is possible to find the *nearest local minimum*. “Nearest” in this
context does not always imply “nearest” in a geometrical sense (*i.e.*,
the least sum of square coordinate differences), but means the minimum
that can be reached by systematically moving down the steepest local
gradient. Finding this nearest local minimum is all that GROMACS can do
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
potentials \ `17 <#ref-Nilges88>`__), but always take care to restore
the correct functions slowly. The best search method that allows rather
drastic structural changes is to allow excursions into four-dimensional
space \ `18 <#ref-Schaik93>`__, but this requires some extra programming
beyond the standard capabilities of GROMACS.

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
   excessive storage requirements. So GROMACS will shy away from this
   class of methods.

The *steepest descent* method, available in GROMACS, is of the second
class. It simply takes a step in the direction of the negative gradient
(hence in the direction of the force), without any consideration of the
history built up in previous steps. The step size is adjusted such that
the search is fast, but the motion is always downhill. This is a simple
and sturdy, but somewhat stupid, method: its convergence can be quite
slow, especially in the vicinity of the local minimum! The
faster-converging *conjugate gradient method* (see *e.g.*
`19 <#ref-Zimmerman91>`__) uses gradient information from previous
steps. In general, steepest descents will bring you close to the nearest
local minimum very quickly, while conjugate gradients brings you *very*
close to the local minimum, but performs worse far away from the
minimum. GROMACS also supports the L-BFGS minimizer, which is mostly
comparable to *conjugate gradient method*, but in some cases converges
faster.

Definitions and Units
=====================

Notation
--------

The following conventions for mathematical typesetting are used
throughout this document:

We define the *lowercase* subscripts :math:`i`, :math:`j`, :math:`k` and
:math:`l` to denote particles: :math:`{{\mbox{\boldmath ${r}$}}_i}` is
the *position vector* of particle :math:`i`, and using this notation:

.. math::

   \begin{aligned}
   {{\mbox{\boldmath ${r}$}}_{ij}}=	{{\mbox{\boldmath ${r}$}}_j}-{{\mbox{\boldmath ${r}$}}_i}\\
   {r_{ij}}=	| {{\mbox{\boldmath ${r}$}}_{ij}}|\end{aligned}

 The force on particle :math:`i` is denoted by
:math:`{\mbox{\boldmath ${F}$}}_i` and

.. math:: {\mbox{\boldmath ${F}$}}_{ij} = \mbox{force on $i$ exerted by $j$}

 Please note that we changed notation as of version 2.0 to
:math:`{{\mbox{\boldmath ${r}$}}_{ij}}={{\mbox{\boldmath ${r}$}}_j}-{{\mbox{\boldmath ${r}$}}_i}`
since this is the notation commonly used. If you encounter an error, let
us know.

MD units
--------

GROMACS uses a consistent set of units that produce values in the
vicinity of unity for most relevant molecular quantities. Let us call
them *MD units*. The basic units in this system are nm, ps, K, electron
charge (e) and atomic mass unit (u), see Table [tab:basicunits]. The
values used in GROMACS are taken from the CODATA Internationally
recommended 2010 values of fundamental physical constants (see
``http://nist.gov``).

Consistent with these units are a set of derived units, given in
Table [tab:derivedunits].

The **electric conversion factor** :math:`f=\frac{1}{4 \pi
\varepsilon_o}={138.935\,458}` kJ mol\ :math:`^{-1}` nm e:math:`^{-2}`.
It relates the mechanical quantities to the electrical quantities as in

.. math:: V = f \frac{q^2}{r} \mbox{\ \ or\ \ } F = f \frac{q^2}{r^2}

Electric potentials :math:`\Phi` and electric fields
:math:`{\mbox{\boldmath ${E}$}}` are intermediate quantities in the
calculation of energies and forces. They do not occur inside GROMACS. If
they are used in evaluations, there is a choice of equations and related
units. We strongly recommend following the usual practice of including
the factor :math:`f` in expressions that evaluate :math:`\Phi` and
:math:`{\mbox{\boldmath ${E}$}}`:

.. math::

   \begin{aligned}
   \Phi({\mbox{\boldmath ${r}$}}) = f \sum_j \frac{q_j}{|{\mbox{\boldmath ${r}$}}-{\mbox{\boldmath ${r}$}}_j|} 	\\
   {\mbox{\boldmath ${E}$}}({\mbox{\boldmath ${r}$}}) = f \sum_j q_j \frac{({\mbox{\boldmath ${r}$}}-{\mbox{\boldmath ${r}$}}_j)}{|{\mbox{\boldmath ${r}$}}-{\mbox{\boldmath ${r}$}}_j|^3}\end{aligned}

 With these definitions, :math:`q\Phi` is an energy and
:math:`q{\mbox{\boldmath ${E}$}}` is a force. The units are those given
in Table [tab:derivedunits]: about 10 mV for potential. Thus, the
potential of an electronic charge at a distance of 1 nm equals
:math:`f \approx 140` units :math:`\approx
1.4` V. (exact value: :math:`1.439\,964\,5` V)

**Note** that these units are mutually consistent; changing any of the
units is likely to produce inconsistencies and is therefore *strongly
discouraged*! In particular: if Å are used instead of nm, the unit of
time changes to 0.1 ps. If kcal mol\ :math:`^{-1}` (= 4.184 kJ
mol\ :math:`^{-1}`) is used instead of kJ mol\ :math:`^{-1}` for energy,
the unit of time becomes 0.488882 ps and the unit of temperature changes
to 4.184 K. But in both cases all electrical energies go wrong, because
they will still be computed in kJ mol\ :math:`^{-1}`, expecting nm as
the unit of length. Although careful rescaling of charges may still
yield consistency, it is clear that such confusions must be rigidly
avoided.

In terms of the MD units, the usual physical constants take on different
values (see Table [tab:consts]). All quantities are per mol rather than
per molecule. There is no distinction between Boltzmann’s constant
:math:`k` and the gas constant :math:`R`: their value is
:math:`0.008\,314\,462\,1` kJ mol:math:`^{-1}` K:math:`^{-1}`.

Reduced units
-------------

When simulating Lennard-Jones (LJ) systems, it might be advantageous to
use reduced units (*i.e.*, setting
:math:`\epsilon_{ii}=\sigma_{ii}=m_i=k_B=1` for one type of atoms). This
is possible. When specifying the input in reduced units, the output will
also be in reduced units. The one exception is the *temperature*, which
is expressed in :math:`0.008\,314\,462\,1` reduced units. This is a
consequence of using Boltzmann’s constant in the evaluation of
temperature in the code. Thus not :math:`T`, but :math:`k_BT`, is the
reduced temperature. A GROMACS temperature :math:`T=1` means a reduced
temperature of :math:`0.008\ldots` units; if a reduced temperature of 1
is required, the GROMACS temperature should be :math:`120.272\,36`.

In Table [tab:reduced] quantities are given for LJ potentials:

.. math:: V_{LJ} = 4\epsilon \left[ \left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6} \right]

Mixed or Double precision
-------------------------

GROMACS can be compiled in either mixed or double precision.
Documentation of previous GROMACS versions referred to “single
precision”, but the implementation has made selective use of double
precision for many years. Using single precision for all variables would
lead to a significant reduction in accuracy. Although in “mixed
precision” all state vectors, i.e. particle coordinates, velocities and
forces, are stored in single precision, critical variables are double
precision. A typical example of the latter is the virial, which is a sum
over all forces in the system, which have varying signs. In addition, in
many parts of the code we managed to avoid double precision for
arithmetic, by paying attention to summation order or reorganization of
mathematical expressions. The default configuration uses mixed
precision, but it is easy to turn on double precision by adding the
option -DGMX\_DOUBLE=on to cmake. Double precision will be 20 to 100%
slower than mixed precision depending on the architecture you are
running on. Double precision will use somewhat more memory and run
input, energy and full-precision trajectory files will be almost twice
as large.

The energies in mixed precision are accurate up to the last decimal, the
last one or two decimals of the forces are non-significant. The virial
is less accurate than the forces, since the virial is only one order of
magnitude larger than the size of each element in the sum over all atoms
(sec. [sec:virial]). In most cases this is not really a problem, since
the fluctuations in the virial can be two orders of magnitude larger
than the average. Using cut-offs for the Coulomb interactions cause
large errors in the energies, forces, and virial. Even when using a
reaction-field or lattice sum method, the errors are larger than, or
comparable to, the errors due to the partial use of single precision.
Since MD is chaotic, trajectories with very similar starting conditions
will diverge rapidly, the divergence is faster in mixed precision than
in double precision.

For most simulations, mixed precision is accurate enough. In some cases
double precision is required to get reasonable results:

-  normal mode analysis, for the conjugate gradient or l-bfgs
   minimization and the calculation and diagonalization of the Hessian

-  long-term energy conservation, especially for large systems

Algorithms
==========

Introduction
------------

In this chapter we first give describe some general concepts used in
GROMACS: *periodic boundary conditions* (sec. [sec:pbc]) and the *group
concept* (sec. [sec:groupconcept]). The MD algorithm is described in
sec. [sec:MD]: first a global form of the algorithm is given, which is
refined in subsequent subsections. The (simple) EM (Energy Minimization)
algorithm is described in sec. [sec:EM]. Some other algorithms for
special purpose dynamics are described after this.

A few issues are of general interest. In all cases the *system* must be
defined, consisting of molecules. Molecules again consist of particles
with defined interaction functions. The detailed description of the
*topology* of the molecules and of the *force field* and the calculation
of forces is given in chapter [ch:ff]. In the present chapter we
describe other aspects of the algorithm, such as pair list generation,
update of velocities and positions, coupling to external temperature and
pressure, conservation of constraints. The *analysis* of the data
generated by an MD simulation is treated in chapter [ch:analysis].

Periodic boundary conditions
----------------------------

.. figure:: plots/pbctric
   :alt: Periodic boundary conditions in two dimensions.
   :width: 9.00000cm

   Periodic boundary conditions in two dimensions.

The classical way to minimize edge effects in a finite system is to
apply *periodic boundary conditions*. The atoms of the system to be
simulated are put into a space-filling box, which is surrounded by
translated copies of itself (Fig. [fig:pbc]). Thus there are no
boundaries of the system; the artifact caused by unwanted boundaries in
an isolated cluster is now replaced by the artifact of periodic
conditions. If the system is crystalline, such boundary conditions are
desired (although motions are naturally restricted to periodic motions
with wavelengths fitting into the box). If one wishes to simulate
non-periodic systems, such as liquids or solutions, the periodicity by
itself causes errors. The errors can be evaluated by comparing various
system sizes; they are expected to be less severe than the errors
resulting from an unnatural boundary with vacuum.

There are several possible shapes for space-filling unit cells. Some,
like the *rhombic dodecahedron* and the *truncated
octahedron* `20 <#ref-Adams79>`__ are closer to being a sphere than a
cube is, and are therefore better suited to the study of an
approximately spherical macromolecule in solution, since fewer solvent
molecules are required to fill the box given a minimum distance between
macromolecular images. At the same time, rhombic dodecahedra and
truncated octahedra are special cases of *triclinic* unit cells; the
most general space-filling unit cells that comprise all possible
space-filling shapes \ `21 <#ref-Bekker95>`__. For this reason, GROMACS
is based on the triclinic unit cell.

GROMACS uses periodic boundary conditions, combined with the *minimum
image convention*: only one – the nearest – image of each particle is
considered for short-range non-bonded interaction terms. For long-range
electrostatic interactions this is not always accurate enough, and
GROMACS therefore also incorporates lattice sum methods such as Ewald
Sum, PME and PPPM.

GROMACS supports triclinic boxes of any shape. The simulation box (unit
cell) is defined by the 3 box vectors :math:`{\bf a}`,\ :math:`{\bf b}`
and :math:`{\bf c}`. The box vectors must satisfy the following
conditions:

.. math::

   \label{eqn:box_rot}
   a_y = a_z = b_z = 0

.. math::

   \label{eqn:box_shift1}
   a_x>0,~~~~b_y>0,~~~~c_z>0

.. math::

   \label{eqn:box_shift2}
   |b_x| \leq \frac{1}{2} \, a_x,~~~~
   |c_x| \leq \frac{1}{2} \, a_x,~~~~
   |c_y| \leq \frac{1}{2} \, b_y

 Equations [eqn:box\_rot] can always be satisfied by rotating the box.
Inequalities ([eqn:box\_shift1]) and ([eqn:box\_shift2]) can always be
satisfied by adding and subtracting box vectors.

Even when simulating using a triclinic box, GROMACS always keeps the
particles in a brick-shaped volume for efficiency, as illustrated in
Fig. [fig:pbc] for a 2-dimensional system. Therefore, from the output
trajectory it might seem that the simulation was done in a rectangular
box. The program trjconv can be used to convert the trajectory to a
different unit-cell representation.

It is also possible to simulate without periodic boundary conditions,
but it is usually more efficient to simulate an isolated cluster of
molecules in a large periodic box, since fast grid searching can only be
used in a periodic system.

|A rhombic dodecahedron and truncated octahedron (arbitrary
orientations).|     |A rhombic dodecahedron and truncated octahedron
(arbitrary orientations).|

Some useful box types
~~~~~~~~~~~~~~~~~~~~~

The three most useful box types for simulations of solvated systems are
described in Table [tab:boxtypes]. The rhombic dodecahedron
(Fig. [fig:boxshapes]) is the smallest and most regular space-filling
unit cell. Each of the 12 image cells is at the same distance. The
volume is 71% of the volume of a cube having the same image distance.
This saves about 29% of CPU-time when simulating a spherical or flexible
molecule in solvent. There are two different orientations of a rhombic
dodecahedron that satisfy equations [eqn:box\_rot], [eqn:box\_shift1]
and [eqn:box\_shift2]. The program editconf produces the orientation
which has a square intersection with the xy-plane. This orientation was
chosen because the first two box vectors coincide with the x and y-axis,
which is easier to comprehend. The other orientation can be useful for
simulations of membrane proteins. In this case the cross-section with
the xy-plane is a hexagon, which has an area which is 14% smaller than
the area of a square with the same image distance. The height of the box
(:math:`c_z`) should be changed to obtain an optimal spacing. This box
shape not only saves CPU time, it also results in a more uniform
arrangement of the proteins.

Cut-off restrictions
~~~~~~~~~~~~~~~~~~~~

The minimum image convention implies that the cut-off radius used to
truncate non-bonded interactions may not exceed half the shortest box
vector:

.. math::

   \label{eqn:physicalrc}
     R_c < {\frac{1}{2}}\min(\|{\bf a}\|,\|{\bf b}\|,\|{\bf c}\|),

 because otherwise more than one image would be within the cut-off
distance of the force. When a macromolecule, such as a protein, is
studied in solution, this restriction alone is not sufficient: in
principle, a single solvent molecule should not be able to ‘see’ both
sides of the macromolecule. This means that the length of each box
vector must exceed the length of the macromolecule in the direction of
that edge *plus* two times the cut-off radius :math:`R_c`. It is,
however, common to compromise in this respect, and make the solvent
layer somewhat smaller in order to reduce the computational cost. For
efficiency reasons the cut-off with triclinic boxes is more restricted.
For grid search the extra restriction is weak:

.. math::

   \label{eqn:gridrc}
   R_c < \min(a_x,b_y,c_z)

 For simple search the extra restriction is stronger:

.. math::

   \label{eqn:simplerc}
   R_c < {\frac{1}{2}}\min(a_x,b_y,c_z)

Each unit cell (cubic, rectangular or triclinic) is surrounded by 26
translated images. A particular image can therefore always be identified
by an index pointing to one of 27 *translation vectors* and constructed
by applying a translation with the indexed vector (see [subsec:forces]).
Restriction ([eqn:gridrc]) ensures that only 26 images need to be
considered.

The group concept
-----------------

The GROMACS MD and analysis programs use user-defined *groups* of atoms
to perform certain actions on. The maximum number of groups is 256, but
each atom can only belong to six different groups, one each of the
following:

temperature-coupling group
    The temperature coupling parameters (reference temperature, time
    constant, number of degrees of freedom, see [subsec:update]) can be
    defined for each T-coupling group separately. For example, in a
    solvated macromolecule the solvent (that tends to generate more
    heating by force and integration errors) can be coupled with a
    shorter time constant to a bath than is a macromolecule, or a
    surface can be kept cooler than an adsorbing molecule. Many
    different T-coupling groups may be defined. See also center of mass
    groups below.

freeze group
    Atoms that belong to a freeze group are kept stationary in the
    dynamics. This is useful during equilibration, *e.g.* to avoid badly
    placed solvent molecules giving unreasonable kicks to protein atoms,
    although the same effect can also be obtained by putting a
    restraining potential on the atoms that must be protected. The
    freeze option can be used, if desired, on just one or two
    coordinates of an atom, thereby freezing the atoms in a plane or on
    a line. When an atom is partially frozen, constraints will still be
    able to move it, even in a frozen direction. A fully frozen atom can
    not be moved by constraints. Many freeze groups can be defined.
    Frozen coordinates are unaffected by pressure scaling; in some cases
    this can produce unwanted results, particularly when constraints are
    also used (in this case you will get very large pressures).
    Accordingly, it is recommended to avoid combining freeze groups with
    constraints and pressure coupling. For the sake of equilibration it
    could suffice to start with freezing in a constant volume
    simulation, and afterward use position restraints in conjunction
    with constant pressure.

accelerate group
    On each atom in an “accelerate group” an acceleration
    :math:`{\mbox{\boldmath ${a}$}}^g` is imposed. This is equivalent to
    an external force. This feature makes it possible to drive the
    system into a non-equilibrium state and enables the performance of
    non-equilibrium MD and hence to obtain transport properties.

energy-monitor group
    Mutual interactions between all energy-monitor groups are compiled
    during the simulation. This is done separately for Lennard-Jones and
    Coulomb terms. In principle up to 256 groups could be defined, but
    that would lead to 256\ :math:`\times`\ 256 items! Better use this
    concept sparingly.

    All non-bonded interactions between pairs of energy-monitor groups
    can be excluded (see details in the User Guide). Pairs of particles
    from excluded pairs of energy-monitor groups are not put into the
    pair list. This can result in a significant speedup for simulations
    where interactions within or between parts of the system are not
    required.

center of mass group
    In GROMACS the center of mass (COM) motion can be removed, for
    either the complete system or for groups of atoms. The latter is
    useful, *e.g.* for systems where there is limited friction (*e.g.*
    gas systems) to prevent center of mass motion to occur. It makes
    sense to use the same groups for temperature coupling and center of
    mass motion removal.

Compressed position output group
    In order to further reduce the size of the compressed trajectory
    file (.xtc or .tng), it is possible to store only a subset of all
    particles. All x-compression groups that are specified are saved,
    the rest are not. If no such groups are specified, than all atoms
    are saved to the compressed trajectory file.

The use of groups in GROMACS tools is described in
sec. [sec:usinggroups].

Molecular Dynamics
------------------

**THE GLOBAL MD ALGORITHM**

--------------

| 
| **1. Input initial conditions**
| Potential interaction :math:`V` as a function of atom positions
| Positions :math:`{\mbox{\boldmath ${r}$}}` of all atoms in the system
| Velocities :math:`{\mbox{\boldmath ${v}$}}` of all atoms in the system
| :math:`\Downarrow`

--------------

| 
| **repeat 2,3,4** for the required number of steps:

--------------

| 
| **2. Compute forces**
| The force on any atom
| :math:`{\mbox{\boldmath ${F}$}}_i = - \displaystyle\frac{\partial V}{\partial {\mbox{\boldmath ${r}$}}_i}`
| is computed by calculating the force between non-bonded atom pairs:
| :math:`{\mbox{\boldmath ${F}$}}_i = \sum_j {\mbox{\boldmath ${F}$}}_{ij}`
| plus the forces due to bonded interactions (which may depend on 1, 2,
  3, or 4 atoms), plus restraining and/or external forces.
| The potential and kinetic energies and the pressure tensor may be
  computed.
| :math:`\Downarrow`
| **3. Update configuration**
| The movement of the atoms is simulated by numerically solving Newton’s
  equations of motion
| :math:`\displaystyle
  \frac {{\mbox{d}}^2{\mbox{\boldmath ${r}$}}_i}{{\mbox{d}}t^2} = \frac{{\mbox{\boldmath ${F}$}}_i}{m_i} `
| or
| :math:`\displaystyle
  \frac{{\mbox{d}}{\mbox{\boldmath ${r}$}}_i}{{\mbox{d}}t} = {\mbox{\boldmath ${v}$}}_i ; \;\;
  \frac{{\mbox{d}}{\mbox{\boldmath ${v}$}}_i}{{\mbox{d}}t} = \frac{{\mbox{\boldmath ${F}$}}_i}{m_i} `
| :math:`\Downarrow`
| **4.** if required: **Output step**
| write positions, velocities, energies, temperature, pressure, etc.

A global flow scheme for MD is given in Fig. [fig:global]. Each MD or EM
run requires as input a set of initial coordinates and – optionally –
initial velocities of all particles involved. This chapter does not
describe how these are obtained; for the setup of an actual MD run check
the online manual at `www.gromacs.org <http://www.gromacs.org>`__.

Initial conditions
~~~~~~~~~~~~~~~~~~

Topology and force field
^^^^^^^^^^^^^^^^^^^^^^^^

The system topology, including a description of the force field, must be
read in. Force fields and topologies are described in chapter [ch:ff]
and [ch:top], respectively. All this information is static; it is never
modified during the run.

Coordinates and velocities
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: plots/maxwell
   :alt: A Maxwell-Boltzmann velocity distribution, generated from
   random numbers.
   :width: 8.00000cm

   A Maxwell-Boltzmann velocity distribution, generated from random
   numbers.

Then, before a run starts, the box size and the coordinates and
velocities of all particles are required. The box size and shape is
determined by three vectors (nine numbers)
:math:`{\mbox{\boldmath ${b}$}}_1, {\mbox{\boldmath ${b}$}}_2, {\mbox{\boldmath ${b}$}}_3`,
which represent the three basis vectors of the periodic box.

If the run starts at :math:`t=t_0`, the coordinates at :math:`t=t_0`
must be known. The *leap-frog algorithm*, the default algorithm used to
update the time step with :math:`{{\Delta t}}` (see [subsec:update]),
also requires that the velocities at
:math:`t=t_0 - {{\frac{1}{2}}{{\Delta t}}}` are known. If velocities are
not available, the program can generate initial atomic velocities
:math:`v_i, i=1\ldots 3N` with a (Fig. [fig:maxwell]) at a given
absolute temperature :math:`T`:

.. math:: p(v_i) = \sqrt{\frac{m_i}{2 \pi kT}}\exp\left(-\frac{m_i v_i^2}{2kT}\right)

 where :math:`k` is Boltzmann’s constant (see chapter [ch:defunits]). To
accomplish this, normally distributed random numbers are generated by
adding twelve random numbers :math:`R_k` in the range
:math:`0 \le R_k < 1` and subtracting 6.0 from their sum. The result is
then multiplied by the standard deviation of the velocity distribution
:math:`\sqrt{kT/m_i}`. Since the resulting total energy will not
correspond exactly to the required temperature :math:`T`, a correction
is made: first the center-of-mass motion is removed and then all
velocities are scaled so that the total energy corresponds exactly to
:math:`T` (see eqn. [eqn:E-T]).

Center-of-mass motion
^^^^^^^^^^^^^^^^^^^^^

The center-of-mass velocity is normally set to zero at every step; there
is (usually) no net external force acting on the system and the
center-of-mass velocity should remain constant. In practice, however,
the update algorithm introduces a very slow change in the center-of-mass
velocity, and therefore in the total kinetic energy of the system –
especially when temperature coupling is used. If such changes are not
quenched, an appreciable center-of-mass motion can develop in long runs,
and the temperature will be significantly misinterpreted. Something
similar may happen due to overall rotational motion, but only when an
isolated cluster is simulated. In periodic systems with filled boxes,
the overall rotational motion is coupled to other degrees of freedom and
does not cause such problems.

Neighbor searching
~~~~~~~~~~~~~~~~~~

As mentioned in chapter [ch:ff], internal forces are either generated
from fixed (static) lists, or from dynamic lists. The latter consist of
non-bonded interactions between any pair of particles. When calculating
the non-bonded forces, it is convenient to have all particles in a
rectangular box. As shown in Fig. [fig:pbc], it is possible to transform
a triclinic box into a rectangular box. The output coordinates are
always in a rectangular box, even when a dodecahedron or triclinic box
was used for the simulation. Equation [eqn:box\_rot] ensures that we can
reset particles in a rectangular box by first shifting them with box
vector :math:`{\bf c}`, then with :math:`{\bf b}` and finally with
:math:`{\bf a}`. Equations [eqn:box\_shift2], [eqn:physicalrc] and
[eqn:gridrc] ensure that we can find the 14 nearest triclinic images
within a linear combination that does not involve multiples of box
vectors.

Pair lists generation
^^^^^^^^^^^^^^^^^^^^^

The non-bonded pair forces need to be calculated only for those pairs
:math:`i,j` for which the distance :math:`r_{ij}` between :math:`i` and
the nearest image of :math:`j` is less than a given cut-off radius
:math:`R_c`. Some of the particle pairs that fulfill this criterion are
excluded, when their interaction is already fully accounted for by
bonded interactions. GROMACS employs a *pair list* that contains those
particle pairs for which non-bonded forces must be calculated. The pair
list contains particles :math:`i`, a displacement vector for particle
:math:`i`, and all particles :math:`j` that are within ``rlist`` of this
particular image of particle :math:`i`. The list is updated every
``nstlist`` steps.

To make the neighbor list, all particles that are close (*i.e.* within
the neighbor list cut-off) to a given particle must be found. This
searching, usually called neighbor search (NS) or pair search, involves
periodic boundary conditions and determining the *image* (see
sec. [sec:pbc]). The search algorithm is :math:`O(N)`, although a
simpler :math:`O(N^2)` algorithm is still available under some
conditions.

Cut-off schemes: group versus Verlet
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

From version 4.6, GROMACS supports two different cut-off scheme setups:
the original one based on particle groups and one using a Verlet buffer.
There are some important differences that affect results, performance
and feature support. The group scheme can be made to work (almost) like
the Verlet scheme, but this will lead to a decrease in performance. The
group scheme is especially fast for water molecules, which are abundant
in many simulations, but on the most recent x86 processors, this
advantage is negated by the better instruction-level parallelism
available in the Verlet-scheme implementation. The group scheme is
deprecated in version 5.0, and will be removed in a future version. For
practical details of choosing and setting up cut-off schemes, please see
the User Guide.

In the group scheme, a neighbor list is generated consisting of pairs of
groups of at least one particle. These groups were originally charge
groups (see sec. [sec:chargegroup]), but with a proper treatment of
long-range electrostatics, performance in unbuffered simulations is
their only advantage. A pair of groups is put into the neighbor list
when their center of geometry is within the cut-off distance.
Interactions between all particle pairs (one from each charge group) are
calculated for a certain number of MD steps, until the neighbor list is
updated. This setup is efficient, as the neighbor search only checks
distance between charge-group pair, not particle pairs (saves a factor
of :math:`3 \times 3 = 9` with a three-particle water model) and the
non-bonded force kernels can be optimized for, say, a water molecule
“group”. Without explicit buffering, this setup leads to energy drift as
some particle pairs which are within the cut-off don’t interact and some
outside the cut-off do interact. This can be caused by

-  particles moving across the cut-off between neighbor search steps,
   and/or

-  for charge groups consisting of more than one particle, particle
   pairs moving in/out of the cut-off when their charge group center of
   geometry distance is outside/inside of the cut-off.

Explicitly adding a buffer to the neighbor list will remove such
artifacts, but this comes at a high computational cost. How severe the
artifacts are depends on the system, the properties in which you are
interested, and the cut-off setup.

The Verlet cut-off scheme uses a buffered pair list by default. It also
uses clusters of particles, but these are not static as in the group
scheme. Rather, the clusters are defined spatially and consist of 4 or 8
particles, which is convenient for stream computing, using e.g. SSE, AVX
or CUDA on GPUs. At neighbor search steps, a pair list is created with a
Verlet buffer, ie. the pair-list cut-off is larger than the interaction
cut-off. In the non-bonded kernels, interactions are only computed when
a particle pair is within the cut-off distance at that particular time
step. This ensures that as particles move between pair search steps,
forces between nearly all particles within the cut-off distance are
calculated. We say *nearly* all particles, because GROMACS uses a fixed
pair list update frequency for efficiency. A particle-pair, whose
distance was outside the cut-off, could possibly move enough during this
fixed number of steps that its distance is now within the cut-off. This
small chance results in a small energy drift, and the size of the chance
depends on the temperature. When temperature coupling is used, the
buffer size can be determined automatically, given a certain tolerance
on the energy drift.

The Verlet cut-off scheme is implemented in a very efficient fashion
based on clusters of particles. The simplest example is a cluster size
of 4 particles. The pair list is then constructed based on cluster
pairs. The cluster-pair search is much faster searching based on
particle pairs, because :math:`4 \times 4 = 16` particle pairs are put
in the list at once. The non-bonded force calculation kernel can then
calculate many particle-pair interactions at once, which maps nicely to
SIMD or SIMT units on modern hardware, which can perform multiple
floating operations at once. These non-bonded kernels are much faster
than the kernels used in the group scheme for most types of systems,
particularly on newer hardware.

Additionally, when the list buffer is determined automatically as
described below, we also apply dynamic pair list pruning. The pair list
can be constructed infrequently, but that can lead to a lot of pairs in
the list that are outside the cut-off range for all or most of the life
time of this pair list. Such pairs can be pruned out by applying a
cluster-pair kernel that only determines which clusters are in range.
Because of the way the non-bonded data is regularized in GROMACS, this
kernel is an order of magnitude faster than the search and the
interaction kernel. On the GPU this pruning is overlapped with the
integration on the CPU, so it is free in most cases. Therefore we can
prune every 4-10 integration steps with little overhead and
significantly reduce the number of cluster pairs in the interaction
kernel. This procedure is applied automatically, unless the user set the
pair-list buffer size manually.

Energy drift and pair-list buffering
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a canonical (NVT) ensemble, the average energy error caused by
diffusion of :math:`j` particles from outside the pair-list cut-off
:math:`r_\ell` to inside the interaction cut-off :math:`r_c` over the
lifetime of the list can be determined from the atomic displacements and
the shape of the potential at the cut-off. The displacement distribution
along one dimension for a freely moving particle with mass :math:`m`
over time :math:`t` at temperature :math:`T` is a Gaussian :math:`G(x)`
of zero mean and variance :math:`\sigma^2 = t^2 k_B T/m`. For the
distance between two particles, the variance changes to
:math:`\sigma^2 = \sigma_{12}^2 =
t^2 k_B T(1/m_1+1/m_2)`. Note that in practice particles usually
interact with (bump into) other particles over time :math:`t` and
therefore the real displacement distribution is much narrower. Given a
non-bonded interaction cut-off distance of :math:`r_c` and a pair-list
cut-off :math:`r_\ell=r_c+r_b` for :math:`r_b` the Verlet buffer size,
we can then write the average energy error after time :math:`t` for all
missing pair interactions between a single :math:`i` particle of type 1
surrounded by all :math:`j` particles that are of type 2 with number
density :math:`\rho_2`, when the inter-particle distance changes from
:math:`r_0` to :math:`r_t`, as:

.. math::

   \langle \Delta V \rangle =
   \int_{0}^{r_c} \int_{r_\ell}^\infty 4 \pi r_0^2 \rho_2 V(r_t) G\!\left(\frac{r_t-r_0}{\sigma}\right) d r_0\, d r_t

 To evaluate this analytically, we need to make some approximations.
First we replace :math:`V(r_t)` by a Taylor expansion around
:math:`r_c`, then we can move the lower bound of the integral over
:math:`r_0` to :math:`-\infty` which will simplify the result:

.. math::

   \begin{aligned}
   \langle \Delta V \rangle &\approx&
   \int_{-\infty}^{r_c} \int_{r_\ell}^\infty 4 \pi r_0^2 \rho_2 \Big[ V'(r_c) (r_t - r_c) +
   \nonumber\\
   & &
   \phantom{\int_{-\infty}^{r_c} \int_{r_\ell}^\infty 4 \pi r_0^2 \rho_2 \Big[}
   V''(r_c)\frac{1}{2}(r_t - r_c)^2 +
   \nonumber\\
   & &
   \phantom{\int_{-\infty}^{r_c} \int_{r_\ell}^\infty 4 \pi r_0^2 \rho_2 \Big[}
     V'''(r_c)\frac{1}{6}(r_t - r_c)^3 +
     \nonumber\\
   & &
   \phantom{\int_{-\infty}^{r_c} \int_{r_\ell}^\infty 4 \pi r_0^2 \rho_2 \Big[}
     O \! \left((r_t - r_c)^4 \right)\Big] G\!\left(\frac{r_t-r_0}{\sigma}\right) d r_0 \, d r_t\end{aligned}

 Replacing the factor :math:`r_0^2` by :math:`(r_\ell + \sigma)^2`,
which results in a slight overestimate, allows us to calculate the
integrals analytically:

.. math::

   \begin{aligned}
   \langle \Delta V \rangle \!
   &\approx&
   4 \pi (r_\ell+\sigma)^2 \rho_2
   \int_{-\infty}^{r_c} \int_{r_\ell}^\infty \Big[ V'(r_c) (r_t - r_c) +
   \nonumber\\
   & &
   \phantom{4 \pi (r_\ell+\sigma)^2 \rho_2 \int_{-\infty}^{r_c} \int_{r_\ell}^\infty \Big[}
   V''(r_c)\frac{1}{2}(r_t - r_c)^2 +
   \nonumber\\
   & &
   \phantom{4 \pi (r_\ell+\sigma)^2 \rho_2 \int_{-\infty}^{r_c} \int_{r_\ell}^\infty \Big[}
   V'''(r_c)\frac{1}{6}(r_t - r_c)^3 \Big] G\!\left(\frac{r_t-r_0}{\sigma}\right)
   d r_0 \, d r_t\\
   &=&
   4 \pi (r_\ell+\sigma)^2 \rho_2 \bigg\{
   \frac{1}{2}V'(r_c)\left[r_b \sigma G\!\left(\frac{r_b}{\sigma}\right) - (r_b^2+\sigma^2)E\!\left(\frac{r_b}{\sigma}\right) \right] +
   \nonumber\\
   & &
   \phantom{4 \pi (r_\ell+\sigma)^2 \rho_2 \bigg\{ }
   \frac{1}{6}V''(r_c)\left[ \sigma(r_b^2+2\sigma^2) G\!\left(\frac{r_b}{\sigma}\right) - r_b(r_b^2+3\sigma^2 ) E\!\left(\frac{r_b}{\sigma}\right) \right] +
   \nonumber\\
   & &
   \phantom{4 \pi (r_\ell+\sigma)^2 \rho_2 \bigg\{ }
   \frac{1}{24}V'''(r_c)\bigg[ r_b\sigma(r_b^2+5\sigma^2) G\!\left(\frac{r_b}{\sigma}\right)
   \nonumber\\
   & &
   \phantom{4 \pi (r_\ell+\sigma)^2 \rho_2 \bigg\{ \frac{1}{24}V'''(r_c)\bigg[ }
    - (r_b^4+6r_b^2\sigma^2+3\sigma^4 ) E\!\left(\frac{r_b}{\sigma}\right) \bigg]
   \bigg\}\end{aligned}

where :math:`G(x)` is a Gaussian distribution with 0 mean and unit
variance and :math:`E(x)=\frac{1}{2}\mathrm{erfc}(x/\sqrt{2})`. We
always want to achieve small energy error, so :math:`\sigma` will be
small compared to both :math:`r_c` and :math:`r_\ell`, thus the
approximations in the equations above are good, since the Gaussian
distribution decays rapidly. The energy error needs to be averaged over
all particle pair types and weighted with the particle counts. In
GROMACS we don’t allow cancellation of error between pair types, so we
average the absolute values. To obtain the average energy error per unit
time, it needs to be divided by the neighbor-list life time
:math:`t = ({\tt nstlist} - 1)\times{\tt dt}`. The function can not be
inverted analytically, so we use bisection to obtain the buffer size
:math:`r_b` for a target drift. Again we note that in practice the error
we usually be much smaller than this estimate, as in the condensed phase
particle displacements will be much smaller than for freely moving
particles, which is the assumption used here.

When (bond) constraints are present, some particles will have fewer
degrees of freedom. This will reduce the energy errors. For simplicity,
we only consider one constraint per particle, the heaviest particle in
case a particle is involved in multiple constraints. This simplification
overestimates the displacement. The motion of a constrained particle is
a superposition of the 3D motion of the center of mass of both particles
and a 2D rotation around the center of mass. The displacement in an
arbitrary direction of a particle with 2 degrees of freedom is not
Gaussian, but rather follows the complementary error function:

.. math::

   \frac{\sqrt{\pi}}{2\sqrt{2}\sigma}\,\mathrm{erfc}\left(\frac{|r|}{\sqrt{2}\,\sigma}\right)
   \label{eqn:2D_disp}

 where :math:`\sigma^2` is again :math:`t^2 k_B T/m`. This distribution
can no longer be integrated analytically to obtain the energy error. But
we can generate a tight upper bound using a scaled and shifted Gaussian
distribution (not shown). This Gaussian distribution can then be used to
calculate the energy error as described above. The rotation displacement
around the center of mass can not be more than the length of the arm. To
take this into account, we scale :math:`\sigma` in eqn. [eqn:2D\_disp]
(details not presented here) to obtain an overestimate of the real
displacement. This latter effect significantly reduces the buffer size
for longer neighborlist lifetimes in e.g. water, as constrained
hydrogens are by far the fastest particles, but they can not move
further than 0.1 nm from the heavy atom they are connected to.

There is one important implementation detail that reduces the energy
errors caused by the finite Verlet buffer list size. The derivation
above assumes a particle pair-list. However, the GROMACS implementation
uses a cluster pair-list for efficiency. The pair list consists of pairs
of clusters of 4 particles in most cases, also called a
:math:`4 \times 4` list, but the list can also be :math:`4 \times 8`
(GPU CUDA kernels and AVX 256-bit single precision kernels) or
:math:`4 \times 2` (SSE double-precision kernels). This means that the
pair-list is effectively much larger than the corresponding
:math:`1 \times 1` list. Thus slightly beyond the pair-list cut-off
there will still be a large fraction of particle pairs present in the
list. This fraction can be determined in a simulation and accurately
estimated under some reasonable assumptions. The fraction decreases with
increasing pair-list range, meaning that a smaller buffer can be used.
For typical all-atom simulations with a cut-off of 0.9 nm this fraction
is around 0.9, which gives a reduction in the energy errors of a factor
of 10. This reduction is taken into account during the automatic Verlet
buffer calculation and results in a smaller buffer size.

.. figure:: plots/verlet-drift
   :alt: Energy drift per atom for an SPC/E water system at 300K with a
   time step of 2 fs and a pair-list update period of 10 steps
   (pair-list life time: 18 fs). PME was used with ewald-rtol set to
   10\ :math:`^{-5}`; this parameter affects the shape of the potential
   at the cut-off. Error estimates due to finite Verlet buffer size are
   shown for a :math:`1 \times 1` atom pair list and :math:`4 \times 4`
   atom pair list without and with (dashed line) cancellation of
   positive and negative errors. Real energy drift is shown for
   simulations using double- and mixed-precision settings. Rounding
   errors in the SETTLE constraint algorithm from the use of single
   precision causes the drift to become negative at large buffer size.
   Note that at zero buffer size, the real drift is small because
   positive (H-H) and negative (O-H) energy errors cancel.
   :width: 9.00000cm

   Energy drift per atom for an SPC/E water system at 300K with a time
   step of 2 fs and a pair-list update period of 10 steps (pair-list
   life time: 18 fs). PME was used with ewald-rtol set to
   10\ :math:`^{-5}`; this parameter affects the shape of the potential
   at the cut-off. Error estimates due to finite Verlet buffer size are
   shown for a :math:`1 \times 1` atom pair list and :math:`4 \times 4`
   atom pair list without and with (dashed line) cancellation of
   positive and negative errors. Real energy drift is shown for
   simulations using double- and mixed-precision settings. Rounding
   errors in the SETTLE constraint algorithm from the use of single
   precision causes the drift to become negative at large buffer size.
   Note that at zero buffer size, the real drift is small because
   positive (H-H) and negative (O-H) energy errors cancel.

In Fig. [fig:verletdrift] one can see that for small buffer sizes the
drift of the total energy is much smaller than the pair energy error
tolerance, due to cancellation of errors. For larger buffer size, the
error estimate is a factor of 6 higher than drift of the total energy,
or alternatively the buffer estimate is 0.024 nm too large. This is
because the protons don’t move freely over 18 fs, but rather vibrate.

Cut-off artifacts and switched interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With the Verlet scheme, the pair potentials are shifted to be zero at
the cut-off, which makes the potential the integral of the force. This
is only possible in the group scheme if the shape of the potential is
such that its value is zero at the cut-off distance. However, there can
still be energy drift when the forces are non-zero at the cut-off. This
effect is extremely small and often not noticeable, as other integration
errors (e.g. from constraints) may dominate. To completely avoid cut-off
artifacts, the non-bonded forces can be switched exactly to zero at some
distance smaller than the neighbor list cut-off (there are several ways
to do this in GROMACS, see sec. [sec:mod\_nb\_int]). One then has a
buffer with the size equal to the neighbor list cut-off less the longest
interaction cut-off.

Simple search
^^^^^^^^^^^^^

Due to eqns. [eqn:box\_rot] and [eqn:simplerc], the vector
:math:`{{\mbox{\boldmath ${r}$}}_{ij}}` connecting images within the
cut-off :math:`R_c` can be found by constructing:

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${r}$}}'''   & = & {\mbox{\boldmath ${r}$}}_j-{\mbox{\boldmath ${r}$}}_i \\
   {\mbox{\boldmath ${r}$}}''    & = & {\mbox{\boldmath ${r}$}}''' - {\bf c}*\verb'round'(r'''_z/c_z) \\
   {\mbox{\boldmath ${r}$}}'     & = & {\mbox{\boldmath ${r}$}}'' - {\bf b}*\verb'round'(r''_y/b_y) \\
   {\mbox{\boldmath ${r}$}}_{ij} & = & {\mbox{\boldmath ${r}$}}' - {\bf a}*\verb'round'(r'_x/a_x)\end{aligned}

 When distances between two particles in a triclinic box are needed that
do not obey eqn. [eqn:box\_rot], many shifts of combinations of box
vectors need to be considered to find the nearest image.

.. figure:: plots/nstric
   :alt: Grid search in two dimensions. The arrows are the box vectors.
   :width: 8.00000cm

   Grid search in two dimensions. The arrows are the box vectors.

Grid search
^^^^^^^^^^^

The grid search is schematically depicted in Fig. [fig:grid]. All
particles are put on the NS grid, with the smallest spacing :math:`\ge`
:math:`R_c/2` in each of the directions. In the direction of each box
vector, a particle :math:`i` has three images. For each direction the
image may be -1,0 or 1, corresponding to a translation over -1, 0 or +1
box vector. We do not search the surrounding NS grid cells for neighbors
of :math:`i` and then calculate the image, but rather construct the
images first and then search neighbors corresponding to that image of
:math:`i`. As Fig. [fig:grid] shows, some grid cells may be searched
more than once for different images of :math:`i`. This is not a problem,
since, due to the minimum image convention, at most one image will “see”
the :math:`j`-particle. For every particle, fewer than 125 (5:math:`^3`)
neighboring cells are searched. Therefore, the algorithm scales linearly
with the number of particles. Although the prefactor is large, the
scaling behavior makes the algorithm far superior over the standard
:math:`O(N^2)` algorithm when there are more than a few hundred
particles. The grid search is equally fast for rectangular and triclinic
boxes. Thus for most protein and peptide simulations the rhombic
dodecahedron will be the preferred box shape.

Charge groups
^^^^^^^^^^^^^

Charge groups were originally introduced to reduce cut-off artifacts of
Coulomb interactions. When a plain cut-off is used, significant jumps in
the potential and forces arise when atoms with (partial) charges move in
and out of the cut-off radius. When all chemical moieties have a net
charge of zero, these jumps can be reduced by moving groups of atoms
with net charge zero, called charge groups, in and out of the neighbor
list. This reduces the cut-off effects from the charge-charge level to
the dipole-dipole level, which decay much faster. With the advent of
full range electrostatics methods, such as particle-mesh Ewald
(sec. [sec:pme]), the use of charge groups is no longer required for
accuracy. It might even have a slight negative effect on the accuracy or
efficiency, depending on how the neighbor list is made and the
interactions are calculated.

But there is still an important reason for using “charge groups”:
efficiency with the group cut-off scheme. Where applicable, neighbor
searching is carried out on the basis of charge groups which are defined
in the molecular topology. If the nearest image distance between the
*geometrical centers* of the atoms of two charge groups is less than the
cut-off radius, all atom pairs between the charge groups are included in
the pair list. The neighbor searching for a water system, for instance,
is :math:`3^2=9` times faster when each molecule is treated as a charge
group. Also the highly optimized water force loops (see
sec. [sec:waterloops]) only work when all atoms in a water molecule form
a single charge group. Currently the name *neighbor-search group* would
be more appropriate, but the name charge group is retained for
historical reasons. When developing a new force field, the advice is to
use charge groups of 3 to 4 atoms for optimal performance. For all-atom
force fields this is relatively easy, as one can simply put hydrogen
atoms, and in some case oxygen atoms, in the same charge group as the
heavy atom they are connected to; for example: CH\ :math:`_3`,
CH\ :math:`_2`, CH, NH\ :math:`_2`, NH, OH, CO\ :math:`_2`, CO.

With the Verlet cut-off scheme, charge groups are ignored.

Compute forces
~~~~~~~~~~~~~~

Potential energy
^^^^^^^^^^^^^^^^

When forces are computed, the potential energy of each interaction term
is computed as well. The total potential energy is summed for various
contributions, such as Lennard-Jones, Coulomb, and bonded terms. It is
also possible to compute these contributions for *energy-monitor groups*
of atoms that are separately defined (see sec. [sec:groupconcept]).

Kinetic energy and temperature
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The temperature is given by the total kinetic energy of the
:math:`N`-particle system:

.. math:: E_{kin} = {\frac{1}{2}}\sum_{i=1}^N m_i v_i^2

 From this the absolute temperature :math:`T` can be computed using:

.. math::

   {\frac{1}{2}}N_{\mathrm{df}} kT = E_{\mathrm{kin}}
   \label{eqn:E-T}

 where :math:`k` is Boltzmann’s constant and :math:`N_{df}` is the
number of degrees of freedom which can be computed from:

.. math:: N_{\mathrm{df}}  ~=~     3 N - N_c - N_{\mathrm{com}}

 Here :math:`N_c` is the number of *constraints* imposed on the system.
When performing molecular dynamics :math:`N_{\mathrm{com}}=3` additional
degrees of freedom must be removed, because the three center-of-mass
velocities are constants of the motion, which are usually set to zero.
When simulating in vacuo, the rotation around the center of mass can
also be removed, in this case :math:`N_{\mathrm{com}}=6`. When more than
one temperature-coupling group is used, the number of degrees of freedom
for group :math:`i` is:

.. math:: N^i_{\mathrm{df}}  ~=~  (3 N^i - N^i_c) \frac{3 N - N_c - N_{\mathrm{com}}}{3 N - N_c}

The kinetic energy can also be written as a tensor, which is necessary
for pressure calculation in a triclinic system, or systems where shear
forces are imposed:

.. math:: {\bf E}_{\mathrm{kin}} = {\frac{1}{2}}\sum_i^N m_i {{\mbox{\boldmath ${v}$}}_i}\otimes {{\mbox{\boldmath ${v}$}}_i}

Pressure and virial
^^^^^^^^^^^^^^^^^^^

The pressure tensor **P** is calculated from the difference between
kinetic energy :math:`E_{\mathrm{kin}}` and the virial
:math:`{\bf \Xi}`:

.. math::

   {\bf P} = \frac{2}{V} ({\bf E}_{\mathrm{kin}}-{\bf \Xi})
   \label{eqn:P}

 where :math:`V` is the volume of the computational box. The scalar
pressure :math:`P`, which can be used for pressure coupling in the case
of isotropic systems, is computed as:

.. math:: P       = {\rm trace}({\bf P})/3

The virial :math:`{\bf \Xi}` tensor is defined as:

.. math:: {\bf \Xi} = -{\frac{1}{2}}\sum_{i<j} {{\mbox{\boldmath ${r}$}}_{ij}}\otimes {{\mbox{\boldmath ${F}$}}_{ij}}\label{eqn:Xi}

The GROMACS implementation of the virial computation is described in
sec. [sec:virial].

The leap-frog integrator
~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: plots/leapfrog
   :alt: The Leap-Frog integration method. The algorithm is called
   Leap-Frog because :math:`{\mbox{\boldmath ${r}$}}` and
   :math:`{\mbox{\boldmath ${v}$}}` are leaping like frogs over each
   other’s backs.
   :width: 8.00000cm

   The Leap-Frog integration method. The algorithm is called Leap-Frog
   because :math:`{\mbox{\boldmath ${r}$}}` and
   :math:`{\mbox{\boldmath ${v}$}}` are leaping like frogs over each
   other’s backs.

The default MD integrator in GROMACS is the so-called *leap-frog*
algorithm \ `22 <#ref-Hockney74>`__ for the integration of the equations
of motion. When extremely accurate integration with temperature and/or
pressure coupling is required, the velocity Verlet integrators are also
present and may be preferable (see [subsec:vverlet]). The leap-frog
algorithm uses positions :math:`{\mbox{\boldmath ${r}$}}` at time
:math:`t` and velocities :math:`{\mbox{\boldmath ${v}$}}` at time
:math:`t-{{\frac{1}{2}}{{\Delta t}}}`; it updates positions and
velocities using the forces :math:`{\mbox{\boldmath ${F}$}}(t)`
determined by the positions at time :math:`t` using these relations:

.. math::

   \begin{aligned}
   \label{eqn:leapfrogv}
   {\mbox{\boldmath ${v}$}}(t+{{\frac{1}{2}}{{\Delta t}}})  &~=~&   {\mbox{\boldmath ${v}$}}(t-{{\frac{1}{2}}{{\Delta t}}})+\frac{{{\Delta t}}}{m}{\mbox{\boldmath ${F}$}}(t)   \\
   {\mbox{\boldmath ${r}$}}(t+{{\Delta t}})   &~=~&   {\mbox{\boldmath ${r}$}}(t)+{{\Delta t}}{\mbox{\boldmath ${v}$}}(t+{{\frac{1}{2}}{{\Delta t}}})\end{aligned}

 The algorithm is visualized in Fig. [fig:leapfrog]. It produces
trajectories that are identical to the Verlet \ `23 <#ref-Verlet67>`__
algorithm, whose position-update relation is

.. math:: {\mbox{\boldmath ${r}$}}(t+{{\Delta t}})~=~2{\mbox{\boldmath ${r}$}}(t) - {\mbox{\boldmath ${r}$}}(t-{{\Delta t}}) + \frac{1}{m}{\mbox{\boldmath ${F}$}}(t){{\Delta t}}^2+O({{\Delta t}}^4)

 The algorithm is of third order in :math:`{\mbox{\boldmath ${r}$}}` and
is time-reversible. See ref. \ `24 <#ref-Berendsen86b>`__ for the merits
of this algorithm and comparison with other time integration algorithms.

The equations of motion are modified for temperature coupling and
pressure coupling, and extended to include the conservation of
constraints, all of which are described below.

The velocity Verlet integrator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The velocity Verlet algorithm \ `25 <#ref-Swope82>`__ is also
implemented in GROMACS, though it is not yet fully integrated with all
sets of options. In velocity Verlet, positions
:math:`{\mbox{\boldmath ${r}$}}` and velocities
:math:`{\mbox{\boldmath ${v}$}}` at time :math:`t` are used to integrate
the equations of motion; velocities at the previous half step are not
required.

.. math::

   \begin{aligned}
   \label{eqn:velocityverlet1}
   {\mbox{\boldmath ${v}$}}(t+{{\frac{1}{2}}{{\Delta t}}})  &~=~&   {\mbox{\boldmath ${v}$}}(t)+\frac{{{\Delta t}}}{2m}{\mbox{\boldmath ${F}$}}(t)   \\
   {\mbox{\boldmath ${r}$}}(t+{{\Delta t}})   &~=~&   {\mbox{\boldmath ${r}$}}(t)+{{\Delta t}}\,{\mbox{\boldmath ${v}$}}(t+{{\frac{1}{2}}{{\Delta t}}}) \\
   {\mbox{\boldmath ${v}$}}(t+{{\Delta t}})   &~=~&   {\mbox{\boldmath ${v}$}}(t+{{\frac{1}{2}}{{\Delta t}}})+\frac{{{\Delta t}}}{2m}{\mbox{\boldmath ${F}$}}(t+{{\Delta t}})\end{aligned}

 or, equivalently,

.. math::

   \begin{aligned}
   \label{eqn:velocityverlet2}
   {\mbox{\boldmath ${r}$}}(t+{{\Delta t}})   &~=~&   {\mbox{\boldmath ${r}$}}(t)+ {{\Delta t}}\,{\mbox{\boldmath ${v}$}} + \frac{{{\Delta t}}^2}{2m}{\mbox{\boldmath ${F}$}}(t) \\
   {\mbox{\boldmath ${v}$}}(t+{{\Delta t}})   &~=~&   {\mbox{\boldmath ${v}$}}(t)+ \frac{{{\Delta t}}}{2m}\left[{\mbox{\boldmath ${F}$}}(t) + {\mbox{\boldmath ${F}$}}(t+{{\Delta t}})\right]\end{aligned}

 With no temperature or pressure coupling, and with *corresponding*
starting points, leap-frog and velocity Verlet will generate identical
trajectories, as can easily be verified by hand from the equations
above. Given a single starting file with the *same* starting point
:math:`{\mbox{\boldmath ${x}$}}(0)` and
:math:`{\mbox{\boldmath ${v}$}}(0)`, leap-frog and velocity Verlet will
*not* give identical trajectories, as leap-frog will interpret the
velocities as corresponding to :math:`t=-{{\frac{1}{2}}{{\Delta t}}}`,
while velocity Verlet will interpret them as corresponding to the
timepoint :math:`t=0`.

Understanding reversible integrators: The Trotter decomposition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To further understand the relationship between velocity Verlet and
leap-frog integration, we introduce the reversible Trotter formulation
of dynamics, which is also useful to understanding implementations of
thermostats and barostats in GROMACS.

A system of coupled, first-order differential equations can be evolved
from time :math:`t = 0` to time :math:`t` by applying the evolution
operator

.. math::

   \begin{aligned}
   \Gamma(t) &=& \exp(iLt) \Gamma(0) \nonumber \\
   iL &=& \dot{\Gamma}\cdot \nabla_{\Gamma},\end{aligned}

 where :math:`L` is the Liouville operator, and :math:`\Gamma` is the
multidimensional vector of independent variables (positions and
velocities). A short-time approximation to the true operator, accurate
at time :math:`{{\Delta t}}= t/P`, is applied :math:`P` times in
succession to evolve the system as

.. math:: \Gamma(t) = \prod_{i=1}^P \exp(iL{{\Delta t}}) \Gamma(0)

 For NVE dynamics, the Liouville operator is

.. math::

   \begin{aligned}
   iL = \sum_{i=1}^{N} {{{\mbox{\boldmath{$v$}}}}}_i \cdot \nabla_{{{{\mbox{\boldmath{$r$}}}}}_i} + \sum_{i=1}^N \frac{1}{m_i}{{{\mbox{\boldmath{$F$}}}}}(r_i) \cdot \nabla_{{{{\mbox{\boldmath{$v$}}}}}_i}.\end{aligned}

 This can be split into two additive operators

.. math::

   \begin{aligned}
   iL_1 &=& \sum_{i=1}^N \frac{1}{m_i}{{{\mbox{\boldmath{$F$}}}}}(r_i) \cdot \nabla_{{{{\mbox{\boldmath{$v$}}}}}_i} \nonumber \\
   iL_2 &=& \sum_{i=1}^{N} {{{\mbox{\boldmath{$v$}}}}}_i \cdot \nabla_{{{{\mbox{\boldmath{$r$}}}}}_i} \end{aligned}

 Then a short-time, symmetric, and thus reversible approximation of the
true dynamics will be

.. math::

   \begin{aligned}
   \exp(iL{{\Delta t}}) = \exp(iL_2{{\frac{1}{2}}{{\Delta t}}}) \exp(iL_1{{\Delta t}}) \exp(iL_2{{\frac{1}{2}}{{\Delta t}}}) + \mathcal{O}({{\Delta t}}^3).
   \label{eq:NVE_Trotter}\end{aligned}

 This corresponds to velocity Verlet integration. The first exponential
term over :math:`{{\frac{1}{2}}{{\Delta t}}}` corresponds to a velocity
half-step, the second exponential term over :math:`{{\Delta t}}`
corresponds to a full velocity step, and the last exponential term over
:math:`{{\frac{1}{2}}{{\Delta t}}}` is the final velocity half step. For
future times :math:`t = n{{\Delta t}}`, this becomes

.. math::

   \begin{aligned}
   \exp(iLn{{\Delta t}}) &\approx&  \left(\exp(iL_2{{\frac{1}{2}}{{\Delta t}}}) \exp(iL_1{{\Delta t}}) \exp(iL_2{{\frac{1}{2}}{{\Delta t}}})\right)^n \nonumber \\
                &\approx&  \exp(iL_2{{\frac{1}{2}}{{\Delta t}}}) \bigg(\exp(iL_1{{\Delta t}}) \exp(iL_2{{\Delta t}})\bigg)^{n-1} \nonumber \\
                &       &  \;\;\;\; \exp(iL_1{{\Delta t}}) \exp(iL_2{{\frac{1}{2}}{{\Delta t}}}) \end{aligned}

 This formalism allows us to easily see the difference between the
different flavors of Verlet integrators. The leap-frog integrator can be
seen as starting with Eq. [eq:NVE\_Trotter] with the
:math:`\exp\left(iL_1 {\Delta t}\right)` term, instead of the half-step
velocity term, yielding

.. math::

   \begin{aligned}
   \exp(iLn{\Delta t}) &=& \exp\left(iL_1 {\Delta t}\right) \exp\left(iL_2 {{\Delta t}}\right) + \mathcal{O}({{\Delta t}}^3).\end{aligned}

 Here, the full step in velocity is between
:math:`t-{{\frac{1}{2}}{{\Delta t}}}` and
:math:`t+{{\frac{1}{2}}{{\Delta t}}}`, since it is a combination of the
velocity half steps in velocity Verlet. For future times
:math:`t = n{{\Delta t}}`, this becomes

.. math::

   \begin{aligned}
   \exp(iLn{\Delta t}) &\approx& \bigg(\exp\left(iL_1 {\Delta t}\right) \exp\left(iL_2 {{\Delta t}}\right)  \bigg)^{n}.\end{aligned}

 Although at first this does not appear symmetric, as long as the full
velocity step is between :math:`t-{{\frac{1}{2}}{{\Delta t}}}` and
:math:`t+{{\frac{1}{2}}{{\Delta t}}}`, then this is simply a way of
starting velocity Verlet at a different place in the cycle.

Even though the trajectory and thus potential energies are identical
between leap-frog and velocity Verlet, the kinetic energy and
temperature will not necessarily be the same. Standard velocity Verlet
uses the velocities at the :math:`t` to calculate the kinetic energy and
thus the temperature only at time :math:`t`; the kinetic energy is then
a sum over all particles

.. math::

   \begin{aligned}
   KE_{\mathrm{full}}(t) &=& \sum_i \left(\frac{1}{2m_i}{\mbox{\boldmath ${v}$}}_i(t)\right)^2 \nonumber\\ 
         &=& \sum_i \frac{1}{2m_i}\left(\frac{1}{2}{\mbox{\boldmath ${v}$}}_i(t-{{\frac{1}{2}}{{\Delta t}}})+\frac{1}{2}{\mbox{\boldmath ${v}$}}_i(t+{{\frac{1}{2}}{{\Delta t}}})\right)^2,\end{aligned}

 with the square on the *outside* of the average. Standard leap-frog
calculates the kinetic energy at time :math:`t` based on the average
kinetic energies at the timesteps :math:`t+{{\frac{1}{2}}{{\Delta t}}}`
and :math:`t-{{\frac{1}{2}}{{\Delta t}}}`, or the sum over all particles

.. math::

   \begin{aligned}
   KE_{\mathrm{average}}(t) &=& \sum_i \frac{1}{2m_i}\left(\frac{1}{2}{\mbox{\boldmath ${v}$}}_i(t-{{\frac{1}{2}}{{\Delta t}}})^2+\frac{1}{2}{\mbox{\boldmath ${v}$}}_i(t+{{\frac{1}{2}}{{\Delta t}}})^2\right),\end{aligned}

 where the square is *inside* the average.

A non-standard variant of velocity Verlet which averages the kinetic
energies :math:`KE(t+{{\frac{1}{2}}{{\Delta t}}})` and
:math:`KE(t-{{\frac{1}{2}}{{\Delta t}}})`, exactly like leap-frog, is
also now implemented in GROMACS (as .mdp file option md-vv-avek).
Without temperature and pressure coupling, velocity Verlet with
half-step-averaged kinetic energies and leap-frog will be identical up
to numerical precision. For temperature- and pressure-control schemes,
however, velocity Verlet with half-step-averaged kinetic energies and
leap-frog will be different, as will be discussed in the section in
thermostats and barostats.

The half-step-averaged kinetic energy and temperature are slightly more
accurate for a given step size; the difference in average kinetic
energies using the half-step-averaged kinetic energies (*md* and
*md-vv-avek*) will be closer to the kinetic energy obtained in the limit
of small step size than will the full-step kinetic energy (using
*md-vv*). For NVE simulations, this difference is usually not
significant, since the positions and velocities of the particles are
still identical; it makes a difference in the way the the temperature of
the simulations are *interpreted*, but *not* in the trajectories that
are produced. Although the kinetic energy is more accurate with the
half-step-averaged method, meaning that it changes less as the timestep
gets large, it is also more noisy. The RMS deviation of the total energy
of the system (sum of kinetic plus potential) in the half-step-averaged
kinetic energy case will be higher (about twice as high in most cases)
than the full-step kinetic energy. The drift will still be the same,
however, as again, the trajectories are identical.

For NVT simulations, however, there *will* be a difference, as discussed
in the section on temperature control, since the velocities of the
particles are adjusted such that kinetic energies of the simulations,
which can be calculated either way, reach the distribution corresponding
to the set temperature. In this case, the three methods will not give
identical results.

Because the velocity and position are both defined at the same time
:math:`t` the velocity Verlet integrator can be used for some methods,
especially rigorously correct pressure control methods, that are not
actually possible with leap-frog. The integration itself takes
negligibly more time than leap-frog, but twice as many communication
calls are currently required. In most cases, and especially for large
systems where communication speed is important for parallelization and
differences between thermodynamic ensembles vanish in the :math:`1/N`
limit, and when only NVT ensembles are required, leap-frog will likely
be the preferred integrator. For pressure control simulations where the
fine details of the thermodynamics are important, only velocity Verlet
allows the true ensemble to be calculated. In either case, simulation
with double precision may be required to get fine details of
thermodynamics correct.

Multiple time stepping
~~~~~~~~~~~~~~~~~~~~~~

Several other simulation packages uses multiple time stepping for bonds
and/or the PME mesh forces. In GROMACS we have not implemented this
(yet), since we use a different philosophy. Bonds can be constrained
(which is also a more sound approximation of a physical quantum
oscillator), which allows the smallest time step to be increased to the
larger one. This not only halves the number of force calculations, but
also the update calculations. For even larger time steps, angle
vibrations involving hydrogen atoms can be removed using virtual
interaction sites (see sec. [sec:rmfast]), which brings the shortest
time step up to PME mesh update frequency of a multiple time stepping
scheme.

Temperature coupling
~~~~~~~~~~~~~~~~~~~~

While direct use of molecular dynamics gives rise to the NVE (constant
number, constant volume, constant energy ensemble), most quantities that
we wish to calculate are actually from a constant temperature (NVT)
ensemble, also called the canonical ensemble. GROMACS can use the
*weak-coupling* scheme of Berendsen \ `26 <#ref-Berendsen84>`__,
stochastic randomization through the Andersen
thermostat \ `27 <#ref-Andersen80>`__, the extended ensemble Nosé-Hoover
scheme \ `28 <#ref-Nose84>`__, `29 <#ref-Hoover85>`__, or a
velocity-rescaling scheme \ `30 <#ref-Bussi2007a>`__ to simulate
constant temperature, with advantages of each of the schemes laid out
below.

There are several other reasons why it might be necessary to control the
temperature of the system (drift during equilibration, drift as a result
of force truncation and integration errors, heating due to external or
frictional forces), but this is not entirely correct to do from a
thermodynamic standpoint, and in some cases only masks the symptoms
(increase in temperature of the system) rather than the underlying
problem (deviations from correct physics in the dynamics). For larger
systems, errors in ensemble averages and structural properties incurred
by using temperature control to remove slow drifts in temperature appear
to be negligible, but no completely comprehensive comparisons have been
carried out, and some caution must be taking in interpreting the
results.

When using temperature and/or pressure coupling the total energy is no
longer conserved. Instead there is a conserved energy quantity the
formula of which will depend on the combination or temperature and
pressure coupling algorithm used. For all coupling algorithms, except
for Andersen temperature coupling and Parrinello-Rahman pressure
coupling combined with shear stress, the conserved energy quantity is
computed and stored in the energy and log file. Note that this quantity
will not be conserved when external forces are applied to the system,
such as pulling on group with a changing distance or an electric field.
Furthermore, how well the energy is conserved depends on the accuracy of
all algorithms involved in the simulation. Usually the algorithms that
cause most drift are constraints and the pair-list buffer, depending on
the parameters used.

Berendsen temperature coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Berendsen algorithm mimics weak coupling with first-order kinetics
to an external heat bath with given temperature :math:`T_0`. See
ref. \ `31 <#ref-Berendsen91>`__ for a comparison with the Nosé-Hoover
scheme. The effect of this algorithm is that a deviation of the system
temperature from :math:`T_0` is slowly corrected according to:

.. math::

   \frac{{\mbox{d}}T}{{\mbox{d}}t} = \frac{T_0-T}{\tau}
   \label{eqn:Tcoupling}

 which means that a temperature deviation decays exponentially with a
time constant :math:`\tau`. This method of coupling has the advantage
that the strength of the coupling can be varied and adapted to the user
requirement: for equilibration purposes the coupling time can be taken
quite short (*e.g.* 0.01 ps), but for reliable equilibrium runs it can
be taken much longer (*e.g.* 0.5 ps) in which case it hardly influences
the conservative dynamics.

The Berendsen thermostat suppresses the fluctuations of the kinetic
energy. This means that one does not generate a proper canonical
ensemble, so rigorously, the sampling will be incorrect. This error
scales with :math:`1/N`, so for very large systems most ensemble
averages will not be affected significantly, except for the distribution
of the kinetic energy itself. However, fluctuation properties, such as
the heat capacity, will be affected. A similar thermostat which does
produce a correct ensemble is the velocity rescaling
thermostat \ `30 <#ref-Bussi2007a>`__ described below.

The heat flow into or out of the system is affected by scaling the
velocities of each particle every step, or every :math:`n_\mathrm{TC}`
steps, with a time-dependent factor :math:`\lambda`, given by:

.. math::

   \lambda = \left[ 1 + \frac{n_\mathrm{TC} \Delta t}{\tau_T}
   \left\{\frac{T_0}{T(t -  {{\frac{1}{2}}{{\Delta t}}})} - 1 \right\} \right]^{1/2}
   \label{eqn:lambda}

 The parameter :math:`\tau_T` is close, but not exactly equal, to the
time constant :math:`\tau` of the temperature coupling
(eqn. [eqn:Tcoupling]):

.. math:: \tau = 2 C_V \tau_T / N_{df} k

 where :math:`C_V` is the total heat capacity of the system, :math:`k`
is Boltzmann’s constant, and :math:`N_{df}` is the total number of
degrees of freedom. The reason that :math:`\tau \neq \tau_T` is that the
kinetic energy change caused by scaling the velocities is partly
redistributed between kinetic and potential energy and hence the change
in temperature is less than the scaling energy. In practice, the ratio
:math:`\tau / \tau_T` ranges from 1 (gas) to 2 (harmonic solid) to 3
(water). When we use the term “temperature coupling time constant,” we
mean the parameter :math:`\tau_T`. **Note** that in practice the scaling
factor :math:`\lambda` is limited to the range of 0.8
:math:`<= \lambda <=` 1.25, to avoid scaling by very large numbers which
may crash the simulation. In normal use, :math:`\lambda` will always be
much closer to 1.0.

The thermostat modifies the kinetic energy at each scaling step by:

.. math:: \Delta E_k = (\lambda - 1)^2 E_k

 The sum of these changes over the run needs to subtracted from the
total energy to obtain the conserved energy quantity.

Velocity-rescaling temperature coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The velocity-rescaling thermostat \ `30 <#ref-Bussi2007a>`__ is
essentially a Berendsen thermostat (see above) with an additional
stochastic term that ensures a correct kinetic energy distribution by
modifying it according to

.. math::

   {\mbox{d}}K = (K_0 - K) \frac{{\mbox{d}}t}{\tau_T} + 2 \sqrt{\frac{K K_0}{N_f}} \frac{{\mbox{d}}W}{\sqrt{\tau_T}},
   \label{eqn:vrescale}

 where :math:`K` is the kinetic energy, :math:`N_f` the number of
degrees of freedom and :math:`{\mbox{d}}W` a Wiener process. There are
no additional parameters, except for a random seed. This thermostat
produces a correct canonical ensemble and still has the advantage of the
Berendsen thermostat: first order decay of temperature deviations and no
oscillations.

Andersen thermostat
^^^^^^^^^^^^^^^^^^^

One simple way to maintain a thermostatted ensemble is to take an
:math:`NVE` integrator and periodically re-select the velocities of the
particles from a Maxwell-Boltzmann
distribution. \ `27 <#ref-Andersen80>`__ This can either be done by
randomizing all the velocities simultaneously (massive collision) every
:math:`\tau_T/{{\Delta t}}` steps (andersen-massive), or by randomizing
every particle with some small probability every timestep (andersen),
equal to :math:`{{\Delta t}}/\tau`, where in both cases
:math:`{{\Delta t}}` is the timestep and :math:`\tau_T` is a
characteristic coupling time scale. Because of the way constraints
operate, all particles in the same constraint group must be randomized
simultaneously. Because of parallelization issues, the andersen version
cannot currently (5.0) be used in systems with constraints.
andersen-massive can be used regardless of constraints. This thermostat
is also currently only possible with velocity Verlet algorithms, because
it operates directly on the velocities at each timestep.

This algorithm completely avoids some of the ergodicity issues of other
thermostatting algorithms, as energy cannot flow back and forth between
energetically decoupled components of the system as in velocity scaling
motions. However, it can slow down the kinetics of system by randomizing
correlated motions of the system, including slowing sampling when
:math:`\tau_T` is at moderate levels (less than 10 ps). This algorithm
should therefore generally not be used when examining kinetics or
transport properties of the system. \ `32 <#ref-Basconi2013>`__

Nosé-Hoover temperature coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Berendsen weak-coupling algorithm is extremely efficient for
relaxing a system to the target temperature, but once the system has
reached equilibrium it might be more important to probe a correct
canonical ensemble. This is unfortunately not the case for the
weak-coupling scheme.

To enable canonical ensemble simulations, GROMACS also supports the
extended-ensemble approach first proposed by Nosé `28 <#ref-Nose84>`__
and later modified by Hoover \ `29 <#ref-Hoover85>`__. The system
Hamiltonian is extended by introducing a thermal reservoir and a
friction term in the equations of motion. The friction force is
proportional to the product of each particle’s velocity and a friction
parameter, :math:`\xi`. This friction parameter (or “heat bath”
variable) is a fully dynamic quantity with its own momentum
(:math:`p_{\xi}`) and equation of motion; the time derivative is
calculated from the difference between the current kinetic energy and
the reference temperature.

In this formulation, the particles’ equations of motion in
Fig. [fig:global] are replaced by:

.. math::

   \frac {{\mbox{d}}^2{\mbox{\boldmath ${r}$}}_i}{{\mbox{d}}t^2} = \frac{{\mbox{\boldmath ${F}$}}_i}{m_i} - 
   \frac{p_{\xi}}{Q}\frac{{\mbox{d}}{\mbox{\boldmath ${r}$}}_i}{{\mbox{d}}t} ,
   \label{eqn:NH-eqn-of-motion}

where the equation of motion for the heat bath parameter :math:`\xi` is:

.. math:: \frac {{\mbox{d}}p_{\xi}}{{\mbox{d}}t} = \left( T - T_0 \right).

 The reference temperature is denoted :math:`T_0`, while :math:`T` is
the current instantaneous temperature of the system. The strength of the
coupling is determined by the constant :math:`Q` (usually called the
“mass parameter” of the reservoir) in combination with the reference
temperature.  [1]_

The conserved quantity for the Nosé-Hoover equations of motion is not
the total energy, but rather

.. math::

   \begin{aligned}
   H = \sum_{i=1}^{N} \frac{{{{\mbox{\boldmath{$p$}}}}}_i}{2m_i} + U\left({{{\mbox{\boldmath{$r$}}}}}_1,{{{\mbox{\boldmath{$r$}}}}}_2,\ldots,{{{\mbox{\boldmath{$r$}}}}}_N\right) +\frac{p_{\xi}^2}{2Q} + N_fkT\xi,\end{aligned}

 where :math:`N_f` is the total number of degrees of freedom.

In our opinion, the mass parameter is a somewhat awkward way of
describing coupling strength, especially due to its dependence on
reference temperature (and some implementations even include the number
of degrees of freedom in your system when defining :math:`Q`). To
maintain the coupling strength, one would have to change :math:`Q` in
proportion to the change in reference temperature. For this reason, we
prefer to let the GROMACS user work instead with the period
:math:`\tau_T` of the oscillations of kinetic energy between the system
and the reservoir instead. It is directly related to :math:`Q` and
:math:`T_0` via:

.. math:: Q = \frac {\tau_T^2 T_0}{4 \pi^2}.

 This provides a much more intuitive way of selecting the Nosé-Hoover
coupling strength (similar to the weak-coupling relaxation), and in
addition :math:`\tau_T` is independent of system size and reference
temperature.

It is however important to keep the difference between the weak-coupling
scheme and the Nosé-Hoover algorithm in mind: Using weak coupling you
get a strongly damped *exponential relaxation*, while the Nosé-Hoover
approach produces an *oscillatory relaxation*. The actual time it takes
to relax with Nosé-Hoover coupling is several times larger than the
period of the oscillations that you select. These oscillations (in
contrast to exponential relaxation) also means that the time constant
normally should be 4–5 times larger than the relaxation time used with
weak coupling, but your mileage may vary.

Nosé-Hoover dynamics in simple systems such as collections of harmonic
oscillators, can be *nonergodic*, meaning that only a subsection of
phase space is ever sampled, even if the simulations were to run for
infinitely long. For this reason, the Nosé-Hoover chain approach was
developed, where each of the Nosé-Hoover thermostats has its own
Nosé-Hoover thermostat controlling its temperature. In the limit of an
infinite chain of thermostats, the dynamics are guaranteed to be
ergodic. Using just a few chains can greatly improve the ergodicity, but
recent research has shown that the system will still be nonergodic, and
it is still not entirely clear what the practical effect of
this \ `33 <#ref-Cooke2008>`__. Currently, the default number of chains
is 10, but this can be controlled by the user. In the case of chains,
the equations are modified in the following way to include a chain of
thermostatting particles \ `34 <#ref-Martyna1992>`__:

.. math::

   \begin{aligned}
   \frac {{\mbox{d}}^2{\mbox{\boldmath ${r}$}}_i}{{\mbox{d}}t^2} &~=~& \frac{{\mbox{\boldmath ${F}$}}_i}{m_i} - \frac{p_{{\xi}_1}}{Q_1} \frac{{\mbox{d}}{\mbox{\boldmath ${r}$}}_i}{{\mbox{d}}t} \nonumber \\
   \frac {{\mbox{d}}p_{{\xi}_1}}{{\mbox{d}}t} &~=~& \left( T - T_0 \right) - p_{{\xi}_1} \frac{p_{{\xi}_2}}{Q_2} \nonumber \\
   \frac {{\mbox{d}}p_{{\xi}_{i=2\ldots N}}}{{\mbox{d}}t} &~=~& \left(\frac{p_{\xi_{i-1}}^2}{Q_{i-1}} -kT\right) - p_{\xi_i} \frac{p_{\xi_{i+1}}}{Q_{i+1}} \nonumber \\
   \frac {{\mbox{d}}p_{\xi_N}}{{\mbox{d}}t} &~=~& \left(\frac{p_{\xi_{N-1}}^2}{Q_{N-1}}-kT\right)
   \label{eqn:NH-chain-eqn-of-motion}\end{aligned}

The conserved quantity for Nosé-Hoover chains is

.. math::

   \begin{aligned}
   H = \sum_{i=1}^{N} \frac{{{{\mbox{\boldmath{$p$}}}}}_i}{2m_i} + U\left({{{\mbox{\boldmath{$r$}}}}}_1,{{{\mbox{\boldmath{$r$}}}}}_2,\ldots,{{{\mbox{\boldmath{$r$}}}}}_N\right) +\sum_{k=1}^M\frac{p^2_{\xi_k}}{2Q^{\prime}_k} + N_fkT\xi_1 + kT\sum_{k=2}^M \xi_k \end{aligned}

 The values and velocities of the Nosé-Hoover thermostat variables are
generally not included in the output, as they take up a fair amount of
space and are generally not important for analysis of simulations, but
by setting an mdp option the values of all the positions and velocities
of all Nosé-Hoover particles in the chain are written to the .edr file.
Leap-frog simulations currently can only have Nosé-Hoover chain lengths
of 1, but this will likely be updated in later version.

As described in the integrator section, for temperature coupling, the
temperature that the algorithm attempts to match to the reference
temperature is calculated differently in velocity Verlet and leap-frog
dynamics. Velocity Verlet (*md-vv*) uses the full-step kinetic energy,
while leap-frog and *md-vv-avek* use the half-step-averaged kinetic
energy.

We can examine the Trotter decomposition again to better understand the
differences between these constant-temperature integrators. In the case
of Nosé-Hoover dynamics (for simplicity, using a chain with :math:`N=1`,
with more details in Ref. \ `35 <#ref-Martyna1996>`__), we split the
Liouville operator as

.. math:: iL = iL_1 + iL_2 + iL_{\mathrm{NHC}},

 where

.. math::

   \begin{aligned}
   iL_1 &=& \sum_{i=1}^N \left[\frac{{{{\mbox{\boldmath{$p$}}}}}_i}{m_i}\right]\cdot \frac{\partial}{\partial {{{\mbox{\boldmath{$r$}}}}}_i} \nonumber \\
   iL_2 &=& \sum_{i=1}^N {{{\mbox{\boldmath{$F$}}}}}_i\cdot \frac{\partial}{\partial {{{\mbox{\boldmath{$p$}}}}}_i} \nonumber \\
   iL_{\mathrm{NHC}} &=& \sum_{i=1}^N-\frac{p_{\xi}}{Q}{{{\mbox{\boldmath{$v$}}}}}_i\cdot \nabla_{{{{\mbox{\boldmath{$v$}}}}}_i} +\frac{p_{\xi}}{Q}\frac{\partial }{\partial \xi} + \left( T - T_0 \right)\frac{\partial }{\partial p_{\xi}}\end{aligned}

 For standard velocity Verlet with Nosé-Hoover temperature control, this
becomes

.. math::

   \begin{aligned}
   \exp(iL{\Delta t}) &=& \exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) \exp\left(iL_2 {\Delta t}/2\right) \nonumber \\
   &&\exp\left(iL_1 {\Delta t}\right) \exp\left(iL_2 {\Delta t}/2\right) \exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) + \mathcal{O}({{\Delta t}}^3).\end{aligned}

 For half-step-averaged temperature control using *md-vv-avek*, this
decomposition will not work, since we do not have the full step
temperature until after the second velocity step. However, we can
construct an alternate decomposition that is still reversible, by
switching the place of the NHC and velocity portions of the
decomposition:

.. math::

   \begin{aligned}
   \exp(iL{\Delta t}) &=& \exp\left(iL_2 {\Delta t}/2\right) \exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right)\exp\left(iL_1 {\Delta t}\right)\nonumber \\
   &&\exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) \exp\left(iL_2 {\Delta t}/2\right)+ \mathcal{O}({{\Delta t}}^3)
   \label{eq:half_step_NHC_integrator}\end{aligned}

 This formalism allows us to easily see the difference between the
different flavors of velocity Verlet integrator. The leap-frog
integrator can be seen as starting with
Eq. [eq:half\_step\_NHC\_integrator] just before the
:math:`\exp\left(iL_1
{\Delta t}\right)` term, yielding:

.. math::

   \begin{aligned}
   \exp(iL{\Delta t}) &=&  \exp\left(iL_1 {\Delta t}\right) \exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) \nonumber \\
   &&\exp\left(iL_2 {\Delta t}\right) \exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) + \mathcal{O}({{\Delta t}}^3)\end{aligned}

 and then using some algebra tricks to solve for some quantities are
required before they are actually calculated \ `36 <#ref-Holian95>`__.

Group temperature coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^

In GROMACS temperature coupling can be performed on groups of atoms,
typically a protein and solvent. The reason such algorithms were
introduced is that energy exchange between different components is not
perfect, due to different effects including cut-offs etc. If now the
whole system is coupled to one heat bath, water (which experiences the
largest cut-off noise) will tend to heat up and the protein will cool
down. Typically 100 K differences can be obtained. With the use of
proper electrostatic methods (PME) these difference are much smaller but
still not negligible. The parameters for temperature coupling in groups
are given in the mdp file. Recent investigation has shown that small
temperature differences between protein and water may actually be an
artifact of the way temperature is calculated when there are finite
timesteps, and very large differences in temperature are likely a sign
of something else seriously going wrong with the system, and should be
investigated carefully \ `37 <#ref-Eastwood2010>`__.

One special case should be mentioned: it is possible to
temperature-couple only part of the system, leaving other parts without
temperature coupling. This is done by specifying :math:`{-1}` for the
time constant :math:`\tau_T` for the group that should not be
thermostatted. If only part of the system is thermostatted, the system
will still eventually converge to an NVT system. In fact, one suggestion
for minimizing errors in the temperature caused by discretized timesteps
is that if constraints on the water are used, then only the water
degrees of freedom should be thermostatted, not protein degrees of
freedom, as the higher frequency modes in the protein can cause larger
deviations from the “true” temperature, the temperature obtained with
small timesteps \ `37 <#ref-Eastwood2010>`__.

Pressure coupling
~~~~~~~~~~~~~~~~~

In the same spirit as the temperature coupling, the system can also be
coupled to a “pressure bath.” GROMACS supports both the Berendsen
algorithm \ `26 <#ref-Berendsen84>`__ that scales coordinates and box
vectors every step, the extended-ensemble Parrinello-Rahman
approach \ `38 <#ref-Parrinello81>`__, `39 <#ref-Nose83>`__, and for the
velocity Verlet variants, the Martyna-Tuckerman-Tobias-Klein (MTTK)
implementation of pressure control \ `35 <#ref-Martyna1996>`__.
Parrinello-Rahman and Berendsen can be combined with any of the
temperature coupling methods above. MTTK can only be used with
Nosé-Hoover temperature control. From 5.1 afterwards, it can only used
when the system does not have constraints.

Berendsen pressure coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Berendsen algorithm rescales the coordinates and box vectors every
step, or every :math:`n_\mathrm{PC}` steps, with a matrix :math:`\mu`,
which has the effect of a first-order kinetic relaxation of the pressure
towards a given reference pressure :math:`{\bf P}_0` according to

.. math:: \frac{{\mbox{d}}{\bf P}}{{\mbox{d}}t} = \frac{{\bf P}_0-{\bf P}}{\tau_p}.

 The scaling matrix :math:`\mu` is given by

.. math::

   \mu_{ij}
   = \delta_{ij} - \frac{n_\mathrm{PC}\Delta t}{3\, \tau_p} \beta_{ij} \{P_{0ij} - P_{ij}(t) \}.
   \label{eqn:mu}

 Here, :math:`\beta` is the isothermal compressibility of the system. In
most cases this will be a diagonal matrix, with equal elements on the
diagonal, the value of which is generally not known. It suffices to take
a rough estimate because the value of :math:`\beta` only influences the
non-critical time constant of the pressure relaxation without affecting
the average pressure itself. For water at 1 atm and 300 K
:math:`\beta = 4.6 \times 10^{-10}`
Pa\ :math:`^{-1} = 4.6 \times 10^{-5}` bar\ :math:`^{-1}`, which is
:math:`7.6 \times 10^{-4}` MD units (see chapter [ch:defunits]). Most
other liquids have similar values. When scaling completely
anisotropically, the system has to be rotated in order to obey
eqn. [eqn:box\_rot]. This rotation is approximated in first order in the
scaling, which is usually less than :math:`10^{-4}`. The actual scaling
matrix :math:`\mu'` is

.. math::

   \mbox{\boldmath $\mu'$} = 
   \left(\begin{array}{ccc}
   \mu_{xx} & \mu_{xy} + \mu_{yx} & \mu_{xz} + \mu_{zx} \\
   0        & \mu_{yy}            & \mu_{yz} + \mu_{zy} \\
   0        & 0                   & \mu_{zz}
   \end{array}\right).

 The velocities are neither scaled nor rotated. Since the equations of
motion are modified by pressure coupling, the conserved energy quantity
also needs to be modified. For first order pressure coupling, the work
the barostat applies to the system every step needs to be subtracted
from the total energy to obtain the conserved energy quantity:

.. math::

   - \sum_{i,j} (\mu_{ij} -\delta_{ij}) P_{ij} V =
   \sum_{i,j} 2(\mu_{ij} -\delta_{ij}) \Xi_{ij}

 where :math:`\delta_{ij}` is the Kronecker delta and :math:`{\bf \Xi}`
is the virial. Note that the factor 2 originates from the factor
:math:`\frac{1}{2}` in the virial definition (eqn. [eqn:Xi]).

In GROMACS, the Berendsen scaling can also be done isotropically, which
means that instead of :math:`{\mbox{\boldmath ${P}$}}` a diagonal matrix
with elements of size trace\ :math:`({\mbox{\boldmath ${P}$}})/3` is
used. For systems with interfaces, semi-isotropic scaling can be useful.
In this case, the :math:`x/y`-directions are scaled isotropically and
the :math:`z` direction is scaled independently. The compressibility in
the :math:`x/y` or :math:`z`-direction can be set to zero, to scale only
in the other direction(s).

If you allow full anisotropic deformations and use constraints you might
have to scale more slowly or decrease your timestep to avoid errors from
the constraint algorithms. It is important to note that although the
Berendsen pressure control algorithm yields a simulation with the
correct average pressure, it does not yield the exact NPT ensemble, and
it is not yet clear exactly what errors this approximation may yield.

Parrinello-Rahman pressure coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In cases where the fluctuations in pressure or volume are important *per
se* (*e.g.* to calculate thermodynamic properties), especially for small
systems, it may be a problem that the exact ensemble is not well defined
for the weak-coupling scheme, and that it does not simulate the true NPT
ensemble.

GROMACS also supports constant-pressure simulations using the
Parrinello-Rahman approach \ `38 <#ref-Parrinello81>`__,
`39 <#ref-Nose83>`__, which is similar to the Nosé-Hoover temperature
coupling, and in theory gives the true NPT ensemble. With the
Parrinello-Rahman barostat, the box vectors as represented by the matrix
obey the matrix equation of motion [2]_

.. math:: \frac{{\mbox{d}}{\mbox{\boldmath ${b}$}}^2}{{\mbox{d}}t^2}= V {\mbox{\boldmath ${W}$}}^{-1} {\mbox{\boldmath ${b}$}}'^{-1} \left( {\mbox{\boldmath ${P}$}} - {\mbox{\boldmath ${P}$}}_{ref}\right).

The volume of the box is denoted :math:`V`, and
:math:`{\mbox{\boldmath ${W}$}}` is a matrix parameter that determines
the strength of the coupling. The matrices and :math:`_{ref}` are the
current and reference pressures, respectively.

The equations of motion for the particles are also changed, just as for
the Nosé-Hoover coupling. In most cases you would combine the
Parrinello-Rahman barostat with the Nosé-Hoover thermostat, but to keep
it simple we only show the Parrinello-Rahman modification here. The
modified Hamiltonian, which will be conserved, is:

.. math::

   E_\mathrm{pot} + E_\mathrm{kin} +  \sum_i P_{ii} V +
   \sum_{i,j} \frac{1}{2} W_{ij}  \left( \frac{{\mbox{d}}b_{ij}}{{\mbox{d}}t} \right)^2

 The equations of motion for the atoms, obtained from the Hamiltonian
are:

.. math::

   \begin{aligned}
    \frac {{\mbox{d}}^2{\mbox{\boldmath ${r}$}}_i}{{\mbox{d}}t^2} & = & \frac{{\mbox{\boldmath ${F}$}}_i}{m_i} -
   {\mbox{\boldmath ${M}$}} \frac{{\mbox{d}}{\mbox{\boldmath ${r}$}}_i}{{\mbox{d}}t} , \\ {\mbox{\boldmath ${M}$}} & = & {\mbox{\boldmath ${b}$}}^{-1} \left[
     {\mbox{\boldmath ${b}$}} \frac{{\mbox{d}}{\mbox{\boldmath ${b}$}}'}{{\mbox{d}}t} + \frac{{\mbox{d}}{\mbox{\boldmath ${b}$}}}{{\mbox{d}}t} {\mbox{\boldmath ${b}$}}'
     \right] {\mbox{\boldmath ${b}$}}'^{-1}.
     \end{aligned}

 This extra term has the appearance of a friction, but it should be
noted that it is ficticious, and rather an effect of the
Parrinello-Rahman equations of motion being defined with all particle
coordinates represented relative to the box vectors, while GROMACS uses
normal Cartesian coordinates for positions, velocities and forces. It is
worth noting that the kinetic energy too should formally be calculated
based on velocities relative to the box vectors. This can have an effect
e.g. for external constant stress, but for now we only support coupling
to constant external pressures, and for any normal simulation the
velocities of box vectors should be extremely small compared to particle
velocities. Gang Liu has done some work on deriving this for Cartesian
coordinates\ `40 <#ref-Liu2015>`__ that we will try to implement at some
point in the future together with support for external stress.

The (inverse) mass parameter matrix
:math:`{\mbox{\boldmath ${W}$}}^{-1}` determines the strength of the
coupling, and how the box can be deformed. The box restriction
([eqn:box\_rot]) will be fulfilled automatically if the corresponding
elements of :math:`{\mbox{\boldmath ${W}$}}^{-1}` are zero. Since the
coupling strength also depends on the size of your box, we prefer to
calculate it automatically in GROMACS. You only have to provide the
approximate isothermal compressibilities :math:`\beta` and the pressure
time constant :math:`\tau_p` in the input file (:math:`L` is the largest
box matrix element):

.. math::

   \left(
   {\mbox{\boldmath ${W}$}}^{-1} \right)_{ij} = \frac{4 \pi^2 \beta_{ij}}{3 \tau_p^2 L}.

Just as for the Nosé-Hoover thermostat, you should realize that the
Parrinello-Rahman time constant is *not* equivalent to the relaxation
time used in the Berendsen pressure coupling algorithm. In most cases
you will need to use a 4–5 times larger time constant with
Parrinello-Rahman coupling. If your pressure is very far from
equilibrium, the Parrinello-Rahman coupling may result in very large box
oscillations that could even crash your run. In that case you would have
to increase the time constant, or (better) use the weak-coupling scheme
to reach the target pressure, and then switch to Parrinello-Rahman
coupling once the system is in equilibrium. Additionally, using the
leap-frog algorithm, the pressure at time :math:`t` is not available
until after the time step has completed, and so the pressure from the
previous step must be used, which makes the algorithm not directly
reversible, and may not be appropriate for high precision thermodynamic
calculations.

Surface-tension coupling
^^^^^^^^^^^^^^^^^^^^^^^^

When a periodic system consists of more than one phase, separated by
surfaces which are parallel to the :math:`xy`-plane, the surface tension
and the :math:`z`-component of the pressure can be coupled to a pressure
bath. Presently, this only works with the Berendsen pressure coupling
algorithm in GROMACS. The average surface tension :math:`\gamma(t)` can
be calculated from the difference between the normal and the lateral
pressure

.. math::

   \begin{aligned}
   \gamma(t) & = & 
   \frac{1}{n} \int_0^{L_z}
   \left\{ P_{zz}(z,t) - \frac{P_{xx}(z,t) + P_{yy}(z,t)}{2} \right\} \mbox{d}z \\
   & = &
   \frac{L_z}{n} \left\{ P_{zz}(t) - \frac{P_{xx}(t) + P_{yy}(t)}{2} \right\},\end{aligned}

 where :math:`L_z` is the height of the box and :math:`n` is the number
of surfaces. The pressure in the z-direction is corrected by scaling the
height of the box with :math:`\mu_{zz}`

.. math:: \Delta P_{zz} = \frac{\Delta t}{\tau_p} \{ P_{0zz} - P_{zz}(t) \}

.. math:: \mu_{zz} = 1 + \beta_{zz} \Delta P_{zz}

 This is similar to normal pressure coupling, except that the factor of
:math:`1/3` is missing. The pressure correction in the
:math:`z`-direction is then used to get the correct convergence for the
surface tension to the reference value :math:`\gamma_0`. The correction
factor for the box length in the :math:`x`/:math:`y`-direction is

.. math::

   \mu_{x/y} = 1 + \frac{\Delta t}{2\,\tau_p} \beta_{x/y}
           \left( \frac{n \gamma_0}{\mu_{zz} L_z}
           - \left\{ P_{zz}(t)+\Delta P_{zz} - \frac{P_{xx}(t) + P_{yy}(t)}{2} \right\} 
           \right)

 The value of :math:`\beta_{zz}` is more critical than with normal
pressure coupling. Normally an incorrect compressibility will just scale
:math:`\tau_p`, but with surface tension coupling it affects the
convergence of the surface tension. When :math:`\beta_{zz}` is set to
zero (constant box height), :math:`\Delta P_{zz}` is also set to zero,
which is necessary for obtaining the correct surface tension.

MTTK pressure control algorithms
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As mentioned in the previous section, one weakness of leap-frog
integration is in constant pressure simulations, since the pressure
requires a calculation of both the virial and the kinetic energy at the
full time step; for leap-frog, this information is not available until
*after* the full timestep. Velocity Verlet does allow the calculation,
at the cost of an extra round of global communication, and can compute,
mod any integration errors, the true NPT ensemble.

The full equations, combining both pressure coupling and temperature
coupling, are taken from Martyna *et al.* `35 <#ref-Martyna1996>`__ and
Tuckerman \ `41 <#ref-Tuckerman2006>`__ and are referred to here as MTTK
equations (Martyna-Tuckerman-Tobias-Klein). We introduce for convenience
:math:`\epsilon = (1/3)\ln (V/V_0)`, where :math:`V_0` is a reference
volume. The momentum of :math:`\epsilon` is
:math:`{v_{\epsilon}}= p_{\epsilon}/W =
\dot{\epsilon} = \dot{V}/3V`, and define :math:`\alpha = 1 + 3/N_{dof}`
(see Ref \ `41 <#ref-Tuckerman2006>`__)

The isobaric equations are

.. math::

   \begin{aligned}
   \dot{{{{\mbox{\boldmath{$r$}}}}}}_i &=& \frac{{{{\mbox{\boldmath{$p$}}}}}_i}{m_i} + \frac{{p_{\epsilon}}}{W} {{{\mbox{\boldmath{$r$}}}}}_i \nonumber \\
   \frac{\dot{{{{\mbox{\boldmath{$p$}}}}}}_i}{m_i} &=& \frac{1}{m_i}{{{\mbox{\boldmath{$F$}}}}}_i - \alpha\frac{{p_{\epsilon}}}{W} \frac{{{{\mbox{\boldmath{$p$}}}}}_i}{m_i} \nonumber \\
   \dot{\epsilon} &=& \frac{{p_{\epsilon}}}{W} \nonumber \\
   \frac{\dot{{p_{\epsilon}}}}{W} &=& \frac{3V}{W}(P_{\mathrm{int}} - P) + (\alpha-1)\left(\sum_{n=1}^N\frac{{{{\mbox{\boldmath{$p$}}}}}_i^2}{m_i}\right),\\\end{aligned}

 where

.. math::

   \begin{aligned}
   P_{\mathrm{int}} &=& P_{\mathrm{kin}} -P_{\mathrm{vir}} = \frac{1}{3V}\left[\sum_{i=1}^N \left(\frac{{{{\mbox{\boldmath{$p$}}}}}_i^2}{2m_i} - {{{\mbox{\boldmath{$r$}}}}}_i \cdot {{{\mbox{\boldmath{$F$}}}}}_i\
   \right)\right].\end{aligned}

 The terms including :math:`\alpha` are required to make phase space
incompressible \ `41 <#ref-Tuckerman2006>`__. The :math:`\epsilon`
acceleration term can be rewritten as

.. math::

   \begin{aligned}
   \frac{\dot{{p_{\epsilon}}}}{W} &=& \frac{3V}{W}\left(\alpha P_{\mathrm{kin}} - P_{\mathrm{vir}} - P\right)\end{aligned}

 In terms of velocities, these equations become

.. math::

   \begin{aligned}
   \dot{{{{\mbox{\boldmath{$r$}}}}}}_i &=& {{{\mbox{\boldmath{$v$}}}}}_i + {v_{\epsilon}}{{{\mbox{\boldmath{$r$}}}}}_i \nonumber \\
   \dot{{{{\mbox{\boldmath{$v$}}}}}}_i &=& \frac{1}{m_i}{{{\mbox{\boldmath{$F$}}}}}_i - \alpha{v_{\epsilon}}{{{\mbox{\boldmath{$v$}}}}}_i \nonumber \\
   \dot{\epsilon} &=& {v_{\epsilon}}\nonumber \\
   \dot{{v_{\epsilon}}} &=& \frac{3V}{W}(P_{\mathrm{int}} - P) + (\alpha-1)\left( \sum_{n=1}^N \frac{1}{2} m_i {{{\mbox{\boldmath{$v$}}}}}_i^2\right)\nonumber \\
   P_{\mathrm{int}} &=& P_{\mathrm{kin}} - P_{\mathrm{vir}} = \frac{1}{3V}\left[\sum_{i=1}^N \left(\frac{1}{2} m_i{{{\mbox{\boldmath{$v$}}}}}_i^2 - {{{\mbox{\boldmath{$r$}}}}}_i \cdot {{{\mbox{\boldmath{$F$}}}}}_i\right)\right]\end{aligned}

 For these equations, the conserved quantity is

.. math::

   \begin{aligned}
   H = \sum_{i=1}^{N} \frac{{{{\mbox{\boldmath{$p$}}}}}_i^2}{2m_i} + U\left({{{\mbox{\boldmath{$r$}}}}}_1,{{{\mbox{\boldmath{$r$}}}}}_2,\ldots,{{{\mbox{\boldmath{$r$}}}}}_N\right) + \frac{p_\epsilon}{2W} + PV\end{aligned}

 The next step is to add temperature control. Adding Nosé-Hoover chains,
including to the barostat degree of freedom, where we use :math:`\eta`
for the barostat Nosé-Hoover variables, and :math:`Q^{\prime}` for the
coupling constants of the thermostats of the barostats, we get

.. math::

   \begin{aligned}
   \dot{{{{\mbox{\boldmath{$r$}}}}}}_i &=& \frac{{{{\mbox{\boldmath{$p$}}}}}_i}{m_i} + \frac{{p_{\epsilon}}}{W} {{{\mbox{\boldmath{$r$}}}}}_i \nonumber \\
   \frac{\dot{{{{\mbox{\boldmath{$p$}}}}}}_i}{m_i} &=& \frac{1}{m_i}{{{\mbox{\boldmath{$F$}}}}}_i - \alpha\frac{{p_{\epsilon}}}{W} \frac{{{{\mbox{\boldmath{$p$}}}}}_i}{m_i} - \frac{p_{\xi_1}}{Q_1}\frac{{{{\mbox{\boldmath{$p$}}}}}_i}{m_i}\nonumber \\
   \dot{\epsilon} &=& \frac{{p_{\epsilon}}}{W} \nonumber \\
   \frac{\dot{{p_{\epsilon}}}}{W} &=& \frac{3V}{W}(\alpha P_{\mathrm{kin}} - P_{\mathrm{vir}} - P) -\frac{p_{\eta_1}}{Q^{\prime}_1}{p_{\epsilon}}\nonumber \\
   \dot{\xi}_k &=& \frac{p_{\xi_k}}{Q_k} \nonumber \\ 
   \dot{\eta}_k &=& \frac{p_{\eta_k}}{Q^{\prime}_k} \nonumber \\
   \dot{p}_{\xi_k} &=& G_k - \frac{p_{\xi_{k+1}}}{Q_{k+1}} \;\;\;\; k=1,\ldots, M-1 \nonumber \\ 
   \dot{p}_{\eta_k} &=& G^\prime_k - \frac{p_{\eta_{k+1}}}{Q^\prime_{k+1}} \;\;\;\; k=1,\ldots, M-1 \nonumber \\
   \dot{p}_{\xi_M} &=& G_M \nonumber \\
   \dot{p}_{\eta_M} &=& G^\prime_M, \nonumber \\\end{aligned}

 where

.. math::

   \begin{aligned}
   P_{\mathrm{int}} &=& P_{\mathrm{kin}} - P_{\mathrm{vir}} = \frac{1}{3V}\left[\sum_{i=1}^N \left(\frac{{{{\mbox{\boldmath{$p$}}}}}_i^2}{2m_i} - {{{\mbox{\boldmath{$r$}}}}}_i \cdot {{{\mbox{\boldmath{$F$}}}}}_i\right)\right] \nonumber \\
   G_1  &=& \sum_{i=1}^N \frac{{{{\mbox{\boldmath{$p$}}}}}^2_i}{m_i} - N_f kT \nonumber \\
   G_k  &=&  \frac{p^2_{\xi_{k-1}}}{2Q_{k-1}} - kT \;\; k = 2,\ldots,M \nonumber \\
   G^\prime_1 &=& \frac{{p_{\epsilon}}^2}{2W} - kT \nonumber \\
   G^\prime_k &=& \frac{p^2_{\eta_{k-1}}}{2Q^\prime_{k-1}} - kT \;\; k = 2,\ldots,M\end{aligned}

 The conserved quantity is now

.. math::

   \begin{aligned}
   H = \sum_{i=1}^{N} \frac{{{{\mbox{\boldmath{$p$}}}}}_i}{2m_i} + U\left({{{\mbox{\boldmath{$r$}}}}}_1,{{{\mbox{\boldmath{$r$}}}}}_2,\ldots,{{{\mbox{\boldmath{$r$}}}}}_N\right) + \frac{p^2_\epsilon}{2W} + PV + \nonumber \\
   \sum_{k=1}^M\frac{p^2_{\xi_k}}{2Q_k} +\sum_{k=1}^M\frac{p^2_{\eta_k}}{2Q^{\prime}_k} + N_fkT\xi_1 +  kT\sum_{i=2}^M \xi_k + kT\sum_{k=1}^M \eta_k\end{aligned}

 Returning to the Trotter decomposition formalism, for pressure control
and temperature control \ `35 <#ref-Martyna1996>`__ we get:

.. math::

   \begin{aligned}
   iL = iL_1 + iL_2 + iL_{\epsilon,1} + iL_{\epsilon,2} + iL_{\mathrm{NHC-baro}} + iL_{\mathrm{NHC}}\end{aligned}

 where “NHC-baro” corresponds to the Nosè-Hoover chain of the barostat,
and NHC corresponds to the NHC of the particles,

.. math::

   \begin{aligned}
   iL_1 &=& \sum_{i=1}^N \left[\frac{{{{\mbox{\boldmath{$p$}}}}}_i}{m_i} + \frac{{p_{\epsilon}}}{W}{{{\mbox{\boldmath{$r$}}}}}_i\right]\cdot \frac{\partial}{\partial {{{\mbox{\boldmath{$r$}}}}}_i} \\
   iL_2 &=& \sum_{i=1}^N {{{\mbox{\boldmath{$F$}}}}}_i - \alpha \frac{{p_{\epsilon}}}{W}{{{\mbox{\boldmath{$p$}}}}}_i \cdot \frac{\partial}{\partial {{{\mbox{\boldmath{$p$}}}}}_i} \\
   iL_{\epsilon,1} &=& \frac{p_\epsilon}{W} \frac{\partial}{\partial \epsilon}\\
   iL_{\epsilon,2} &=& G_{\epsilon} \frac{\partial}{\partial p_\epsilon}\end{aligned}

 and where

.. math::

   \begin{aligned}
   G_{\epsilon} = 3V\left(\alpha P_{\mathrm{kin}} - P_{\mathrm{vir}} - P\right)\end{aligned}

 Using the Trotter decomposition, we get

.. math::

   \begin{aligned}
   \exp(iL{\Delta t}) &=& \exp\left(iL_{\mathrm{NHC-baro}}{\Delta t}/2\right)\exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) \nonumber \nonumber \\
   &&\exp\left(iL_{\epsilon,2}{\Delta t}/2\right) \exp\left(iL_2 {\Delta t}/2\right) \nonumber \nonumber \\
   &&\exp\left(iL_{\epsilon,1}{\Delta t}\right) \exp\left(iL_1 {\Delta t}\right) \nonumber \nonumber \\
   &&\exp\left(iL_2 {\Delta t}/2\right) \exp\left(iL_{\epsilon,2}{\Delta t}/2\right) \nonumber \nonumber \\
   &&\exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right)\exp\left(iL_{\mathrm{NHC-baro}}{\Delta t}/2\right) + \mathcal{O}({\Delta t}^3)\end{aligned}

 The action of :math:`\exp\left(iL_1 {\Delta t}\right)` comes from the
solution of the the differential equation
:math:`\dot{{{{\mbox{\boldmath{$r$}}}}}}_i = {{{\mbox{\boldmath{$v$}}}}}_i + {v_{\epsilon}}{{{\mbox{\boldmath{$r$}}}}}_i`
with
:math:`{{{\mbox{\boldmath{$v$}}}}}_i = {{{\mbox{\boldmath{$p$}}}}}_i/m_i`
and :math:`{v_{\epsilon}}` constant with initial condition
:math:`{{{\mbox{\boldmath{$r$}}}}}_i(0)`, evaluate at
:math:`t=\Delta t`. This yields the evolution

.. math:: {{{\mbox{\boldmath{$r$}}}}}_i({\Delta t}) = {{{\mbox{\boldmath{$r$}}}}}_i(0)e^{{v_{\epsilon}}{\Delta t}} + \Delta t {{{\mbox{\boldmath{$v$}}}}}_i(0) e^{{v_{\epsilon}}{\Delta t}/2} {\frac{\sinh{\left( {v_{\epsilon}}{\Delta t}/2\right)}}{{v_{\epsilon}}{\Delta t}/2}}.

 The action of :math:`\exp\left(iL_2 {\Delta t}/2\right)` comes from the
solution of the differential equation
:math:`\dot{{{{\mbox{\boldmath{$v$}}}}}}_i = \frac{{{{\mbox{\boldmath{$F$}}}}}_i}{m_i} -
\alpha{v_{\epsilon}}{{{\mbox{\boldmath{$v$}}}}}_i`, yielding

.. math:: {{{\mbox{\boldmath{$v$}}}}}_i({\Delta t}/2) = {{{\mbox{\boldmath{$v$}}}}}_i(0)e^{-\alpha{v_{\epsilon}}{\Delta t}/2} + \frac{\Delta t}{2m_i}{{{\mbox{\boldmath{$F$}}}}}_i(0) e^{-\alpha{v_{\epsilon}}{\Delta t}/4}{\frac{\sinh{\left( \alpha{v_{\epsilon}}{\Delta t}/4\right)}}{\alpha{v_{\epsilon}}{\Delta t}/4}}.

 *md-vv-avek* uses the full step kinetic energies for determining the
pressure with the pressure control, but the half-step-averaged kinetic
energy for the temperatures, which can be written as a Trotter
decomposition as

.. math::

   \begin{aligned}
   \exp(iL{\Delta t}) &=& \exp\left(iL_{\mathrm{NHC-baro}}{\Delta t}/2\right)\nonumber \exp\left(iL_{\epsilon,2}{\Delta t}/2\right) \exp\left(iL_2 {\Delta t}/2\right) \nonumber \\
   &&\exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) \exp\left(iL_{\epsilon,1}{\Delta t}\right) \exp\left(iL_1 {\Delta t}\right) \exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) \nonumber \\
   &&\exp\left(iL_2 {\Delta t}/2\right) \exp\left(iL_{\epsilon,2}{\Delta t}/2\right) \exp\left(iL_{\mathrm{NHC-baro}}{\Delta t}/2\right) + \mathcal{O}({\Delta t}^3)\end{aligned}

With constraints, the equations become significantly more complicated,
in that each of these equations need to be solved iteratively for the
constraint forces. Before GROMACS 5.1, these iterative constraints were
solved as described in \ `42 <#ref-Yu2010>`__. From GROMACS 5.1 onward,
MTTK with constraints has been removed because of numerical stability
issues with the iterations.

Infrequent evaluation of temperature and pressure coupling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Temperature and pressure control require global communication to compute
the kinetic energy and virial, which can become costly if performed
every step for large systems. We can rearrange the Trotter decomposition
to give alternate symplectic, reversible integrator with the coupling
steps every :math:`n` steps instead of every steps. These new
integrators will diverge if the coupling time step is too large, as the
auxiliary variable integrations will not converge. However, in most
cases, long coupling times are more appropriate, as they disturb the
dynamics less \ `35 <#ref-Martyna1996>`__.

Standard velocity Verlet with Nosé-Hoover temperature control has a
Trotter expansion

.. math::

   \begin{aligned}
   \exp(iL{\Delta t}) &\approx& \exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right) \exp\left(iL_2 {\Delta t}/2\right) \nonumber \\
   &&\exp\left(iL_1 {\Delta t}\right) \exp\left(iL_2 {\Delta t}/2\right) \exp\left(iL_{\mathrm{NHC}}{\Delta t}/2\right).\end{aligned}

 If the Nosé-Hoover chain is sufficiently slow with respect to the
motions of the system, we can write an alternate integrator over
:math:`n` steps for velocity Verlet as

.. math::

   \begin{aligned}
   \exp(iL{\Delta t}) &\approx& (\exp\left(iL_{\mathrm{NHC}}(n{\Delta t}/2)\right)\left[\exp\left(iL_2 {\Delta t}/2\right)\right. \nonumber \\
   &&\left.\exp\left(iL_1 {\Delta t}\right) \exp\left(iL_2 {\Delta t}/2\right)\right]^n \exp\left(iL_{\mathrm{NHC}}(n{\Delta t}/2)\right).\end{aligned}

 For pressure control, this becomes

.. math::

   \begin{aligned}
   \exp(iL{\Delta t}) &\approx& \exp\left(iL_{\mathrm{NHC-baro}}(n{\Delta t}/2)\right)\exp\left(iL_{\mathrm{NHC}}(n{\Delta t}/2)\right) \nonumber \nonumber \\
   &&\exp\left(iL_{\epsilon,2}(n{\Delta t}/2)\right) \left[\exp\left(iL_2 {\Delta t}/2\right)\right. \nonumber \nonumber \\
   &&\exp\left(iL_{\epsilon,1}{\Delta t}\right) \exp\left(iL_1 {\Delta t}\right) \nonumber \nonumber \\
   &&\left.\exp\left(iL_2 {\Delta t}/2\right)\right]^n \exp\left(iL_{\epsilon,2}(n{\Delta t}/2)\right) \nonumber \nonumber \\
   &&\exp\left(iL_{\mathrm{NHC}}(n{\Delta t}/2)\right)\exp\left(iL_{\mathrm{NHC-baro}}(n{\Delta t}/2)\right),\end{aligned}

 where the box volume integration occurs every step, but the auxiliary
variable integrations happen every :math:`n` steps.

The complete update algorithm
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**THE UPDATE ALGORITHM**

--------------

| 
| Given:
| Positions :math:`{\mbox{\boldmath ${r}$}}` of all atoms at time
  :math:`t`
| Velocities :math:`{\mbox{\boldmath ${v}$}}` of all atoms at time
  :math:`t-{{\frac{1}{2}}{{\Delta t}}}`
| Accelerations :math:`{\mbox{\boldmath ${F}$}}/m` on all atoms at time
  :math:`t`.
| (Forces are computed disregarding any constraints)
| Total kinetic energy and virial at :math:`t-{{\Delta t}}`
| :math:`\Downarrow`
| **1.** Compute the scaling factors :math:`\lambda` and :math:`\mu`
| according to eqns. [eqn:lambda] and [eqn:mu]
| :math:`\Downarrow`
| **2.** Update and scale velocities:
  :math:`{\mbox{\boldmath ${v}$}}' =  \lambda ({\mbox{\boldmath ${v}$}} +
  {\mbox{\boldmath ${a}$}} \Delta t)`
| :math:`\Downarrow`
| **3.** Compute new unconstrained coordinates:
  :math:`{\mbox{\boldmath ${r}$}}' = {\mbox{\boldmath ${r}$}} + {\mbox{\boldmath ${v}$}}'
  \Delta t`
| :math:`\Downarrow`
| **4.** Apply constraint algorithm to coordinates:
  constrain(\ :math:`{\mbox{\boldmath ${r}$}}^{'} \rightarrow  {\mbox{\boldmath ${r}$}}'';
  \,  {\mbox{\boldmath ${r}$}}`)
| :math:`\Downarrow`
| **5.** Correct velocities for constraints:
  :math:`{\mbox{\boldmath ${v}$}} = ({\mbox{\boldmath ${r}$}}'' -
  {\mbox{\boldmath ${r}$}}) / \Delta t`
| :math:`\Downarrow`
| **6.** Scale coordinates and box:
  :math:`{\mbox{\boldmath ${r}$}} = \mu {\mbox{\boldmath ${r}$}}''; {\mbox{\boldmath ${b}$}} =
  \mu  {\mbox{\boldmath ${b}$}}`

The complete algorithm for the update of velocities and coordinates is
given using leap-frog in Fig. [fig:complete-update]. The SHAKE algorithm
of step 4 is explained below.

GROMACS has a provision to “freeze” (prevent motion of) selected
particles, which must be defined as a “freeze group.” This is
implemented using a *freeze factor :math:`{\mbox{\boldmath ${f}$}}_g`*,
which is a vector, and differs for each freeze group (see
sec. [sec:groupconcept]). This vector contains only zero (freeze) or one
(don’t freeze). When we take this freeze factor and the external
acceleration :math:`{\mbox{\boldmath ${a}$}}_h` into account the update
algorithm for the velocities becomes

.. math:: {\mbox{\boldmath ${v}$}}(t+{\frac{\Delta t}{2}})~=~{\mbox{\boldmath ${f}$}}_g * \lambda * \left[ {\mbox{\boldmath ${v}$}}(t-{\frac{\Delta t}{2}}) +\frac{{\mbox{\boldmath ${F}$}}(t)}{m}\Delta t + {\mbox{\boldmath ${a}$}}_h \Delta t \right],

 where :math:`g` and :math:`h` are group indices which differ per atom.

Output step
~~~~~~~~~~~

The most important output of the MD run is the *trajectory file*, which
contains particle coordinates and (optionally) velocities at regular
intervals. The trajectory file contains frames that could include
positions, velocities and/or forces, as well as information about the
dimensions of the simulation volume, integration step, integration time,
etc. The interpretation of the time varies with the integrator chosen,
as described above. For Velocity Verlet integrators, velocities labeled
at time :math:`t` are for that time. For other integrators (e.g.
leap-frog, stochastic dynamics), the velocities labeled at time
:math:`t` are for time :math:`t - {{\frac{1}{2}}{{\Delta t}}}`.

Since the trajectory files are lengthy, one should not save every step!
To retain all information it suffices to write a frame every 15 steps,
since at least 30 steps are made per period of the highest frequency in
the system, and Shannon’s sampling theorem states that two samples per
period of the highest frequency in a band-limited signal contain all
available information. But that still gives very long files! So, if the
highest frequencies are not of interest, 10 or 20 samples per ps may
suffice. Be aware of the distortion of high-frequency motions by the
*stroboscopic effect*, called *aliasing*: higher frequencies are
mirrored with respect to the sampling frequency and appear as lower
frequencies.

GROMACS can also write reduced-precision coordinates for a subset of the
simulation system to a special compressed trajectory file format. All
the other tools can read and write this format. See the User Guide for
details on how to set up your .mdp file to have mdrun use this feature.

Shell molecular dynamics
------------------------

GROMACS can simulate polarizability using the shell model of Dick and
Overhauser \ `43 <#ref-Dick58>`__. In such models a shell particle
representing the electronic degrees of freedom is attached to a nucleus
by a spring. The potential energy is minimized with respect to the shell
position at every step of the simulation (see below). Successful
applications of shell models in GROMACS have been published for
:math:`N_2` `44 <#ref-Jordan95>`__ and
water \ `45 <#ref-Maaren2001a>`__.

Optimization of the shell positions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The force :math:`_S` on a shell particle :math:`S` can be decomposed
into two components

.. math:: {\mbox{\boldmath ${F}$}}_S ~=~ {\mbox{\boldmath ${F}$}}_{bond} + {\mbox{\boldmath ${F}$}}_{nb}

 where :math:`_{bond}` denotes the component representing the
polarization energy, usually represented by a harmonic potential and
:math:`_{nb}` is the sum of Coulomb and van der Waals interactions. If
we assume that :math:`_{nb}` is almost constant we can analytically
derive the optimal position of the shell, i.e. where :math:`_S` = 0. If
we have the shell S connected to atom A we have

.. math:: {\mbox{\boldmath ${F}$}}_{bond} ~=~ k_b \left( {\mbox{\boldmath ${x}$}}_S - {\mbox{\boldmath ${x}$}}_A\right).

 In an iterative solver, we have positions :math:`_S(n)` where :math:`n`
is the iteration count. We now have at iteration :math:`n`

.. math:: {\mbox{\boldmath ${F}$}}_{nb} ~=~ {\mbox{\boldmath ${F}$}}_S - k_b \left( {\mbox{\boldmath ${x}$}}_S(n) - {\mbox{\boldmath ${x}$}}_A\right)

 and the optimal position for the shells :math:`x_S(n+1)` thus follows
from

.. math:: {\mbox{\boldmath ${F}$}}_S - k_b \left( {\mbox{\boldmath ${x}$}}_S(n) - {\mbox{\boldmath ${x}$}}_A\right) + k_b \left( {\mbox{\boldmath ${x}$}}_S(n+1) - {\mbox{\boldmath ${x}$}}_A\right) = 0

 if we write

.. math:: \Delta {\mbox{\boldmath ${x}$}}_S = {\mbox{\boldmath ${x}$}}_S(n+1) - {\mbox{\boldmath ${x}$}}_S(n)

 we finally obtain

.. math:: \Delta {\mbox{\boldmath ${x}$}}_S = {\mbox{\boldmath ${F}$}}_S/k_b

 which then yields the algorithm to compute the next trial in the
optimization of shell positions

.. math:: {\mbox{\boldmath ${x}$}}_S(n+1) ~=~ {\mbox{\boldmath ${x}$}}_S(n) + {\mbox{\boldmath ${F}$}}_S/k_b.

Constraint algorithms
---------------------

Constraints can be imposed in GROMACS using LINCS (default) or the
traditional SHAKE method.

SHAKE
~~~~~

The SHAKE \ `46 <#ref-Ryckaert77>`__ algorithm changes a set of
unconstrained coordinates :math:`{\mbox{\boldmath ${r}$}}^{'}` to a set
of coordinates :math:`{\mbox{\boldmath ${r}$}}''` that fulfill a list of
distance constraints, using a set :math:`{\mbox{\boldmath ${r}$}}`
reference, as

.. math:: {\rm SHAKE}({\mbox{\boldmath ${r}$}}^{'} \rightarrow {\mbox{\boldmath ${r}$}}'';\, {\mbox{\boldmath ${r}$}})

 This action is consistent with solving a set of Lagrange multipliers in
the constrained equations of motion. SHAKE needs a *relative tolerance*;
it will continue until all constraints are satisfied within that
relative tolerance. An error message is given if SHAKE cannot reset the
coordinates because the deviation is too large, or if a given number of
iterations is surpassed.

Assume the equations of motion must fulfill :math:`K` holonomic
constraints, expressed as

.. math:: \sigma_k({\mbox{\boldmath ${r}$}}_1 \ldots {\mbox{\boldmath ${r}$}}_N) = 0; \;\; k=1 \ldots K.

 For example,
:math:`({\mbox{\boldmath ${r}$}}_1 - {\mbox{\boldmath ${r}$}}_2)^2 - b^2 = 0`.
Then the forces are defined as

.. math::

   - \frac{\partial}{\partial {\mbox{\boldmath ${r}$}}_i} \left( V + \sum_{k=1}^K \lambda_k
   \sigma_k \right),

 where :math:`\lambda_k` are Lagrange multipliers which must be solved
to fulfill the constraint equations. The second part of this sum
determines the *constraint forces* :math:`{\mbox{\boldmath ${G}$}}_i`,
defined by

.. math::

   {\mbox{\boldmath ${G}$}}_i = -\sum_{k=1}^K \lambda_k \frac{\partial \sigma_k}{\partial
   {\mbox{\boldmath ${r}$}}_i}

 The displacement due to the constraint forces in the leap-frog or
Verlet algorithm is equal to
:math:`({\mbox{\boldmath ${G}$}}_i/m_i)({{\Delta t}})^2`. Solving the
Lagrange multipliers (and hence the displacements) requires the solution
of a set of coupled equations of the second degree. These are solved
iteratively by SHAKE. [subsec:SETTLE] For the special case of rigid
water molecules, that often make up more than 80% of the simulation
system we have implemented the SETTLE
algorithm \ `47 <#ref-Miyamoto92>`__ (sec. [sec:constraints]).

For velocity Verlet, an additional round of constraining must be done,
to constrain the velocities of the second velocity half step, removing
any component of the velocity parallel to the bond vector. This step is
called RATTLE, and is covered in more detail in the original Andersen
paper \ `48 <#ref-Andersen1983a>`__.

LINCS
~~~~~

The LINCS algorithm
^^^^^^^^^^^^^^^^^^^

LINCS is an algorithm that resets bonds to their correct lengths after
an unconstrained update \ `49 <#ref-Hess97>`__. The method is
non-iterative, as it always uses two steps. Although LINCS is based on
matrices, no matrix-matrix multiplications are needed. The method is
more stable and faster than SHAKE, but it can only be used with bond
constraints and isolated angle constraints, such as the proton angle in
OH. Because of its stability, LINCS is especially useful for Brownian
dynamics. LINCS has two parameters, which are explained in the
subsection parameters. The parallel version of LINCS, P-LINCS, is
described in subsection [subsec:plincs].

The LINCS formulas
^^^^^^^^^^^^^^^^^^

We consider a system of :math:`N` particles, with positions given by a
:math:`3N` vector :math:`{\mbox{\boldmath ${r}$}}(t)`. For molecular
dynamics the equations of motion are given by Newton’s Law

.. math::

   \label{eqn:c1}
   {{\mbox{d}}^2 {\mbox{\boldmath ${r}$}} \over {\mbox{d}}t^2} = {{{\mbox{\boldmath ${M}$}}}^{-1}}{\mbox{\boldmath ${F}$}},

 where :math:`{\mbox{\boldmath ${F}$}}` is the :math:`3N` force vector
and :math:`{{\mbox{\boldmath ${M}$}}}` is a :math:`3N \times 3N`
diagonal matrix, containing the masses of the particles. The system is
constrained by :math:`K` time-independent constraint equations

.. math::

   \label{eqn:c2}
   g_i({\mbox{\boldmath ${r}$}}) = | {\mbox{\boldmath ${r}$}}_{i_1}-{\mbox{\boldmath ${r}$}}_{i_2} | - d_i = 0 ~~~~~~i=1,\ldots,K.

In a numerical integration scheme, LINCS is applied after an
unconstrained update, just like SHAKE. The algorithm works in two steps
(see figure Fig. [fig:lincs]). In the first step, the projections of the
new bonds on the old bonds are set to zero. In the second step, a
correction is applied for the lengthening of the bonds due to rotation.
The numerics for the first step and the second step are very similar. A
complete derivation of the algorithm can be found in
`49 <#ref-Hess97>`__. Only a short description of the first step is
given here.

.. figure:: plots/lincs
   :alt: The three position updates needed for one time step. The dashed
   line is the old bond of length :math:`d`, the solid lines are the new
   bonds. :math:`l=d
   \cos \theta` and :math:`p=(2 d^2 - l^2)^{1 \over 2}`.
   :height: 5.00000cm

   The three position updates needed for one time step. The dashed line
   is the old bond of length :math:`d`, the solid lines are the new
   bonds. :math:`l=d
   \cos \theta` and :math:`p=(2 d^2 - l^2)^{1 \over 2}`.

A new notation is introduced for the gradient matrix of the constraint
equations which appears on the right hand side of this equation:

.. math::

   \label{eqn:c3}
   B_{hi} = {{\partial}g_h \over {\partial}r_i}

 Notice that :math:`{{\mbox{\boldmath ${B}$}}}` is a :math:`K \times 3N`
matrix, it contains the directions of the constraints. The following
equation shows how the new constrained coordinates
:math:`{\mbox{\boldmath ${r}$}}_{n+1}` are related to the unconstrained
coordinates :math:`{\mbox{\boldmath ${r}$}}_{n+1}^{unc}` by

.. math::

   \label{eqn:m0}
   \begin{array}{c}
     {\mbox{\boldmath ${r}$}}_{n+1}=({\mbox{\boldmath ${I}$}}-{{\mbox{\boldmath ${T}$}}}_n {\mbox{\boldmath ${B}$}}_n) {\mbox{\boldmath ${r}$}}_{n+1}^{unc} + {{\mbox{\boldmath ${T}$}}}_n {{\mbox{\boldmath ${d}$}}}=  
     \\[2mm]
     {\mbox{\boldmath ${r}$}}_{n+1}^{unc} - 
   {{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}_n ({{\mbox{\boldmath ${B}$}}}_n {{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}_n^T)^{-1} ({{\mbox{\boldmath ${B}$}}}_n {\mbox{\boldmath ${r}$}}_{n+1}^{unc} - {{\mbox{\boldmath ${d}$}}}) 
   \end{array}

 where
:math:`{{\mbox{\boldmath ${T}$}}}= {{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}^T ({{\mbox{\boldmath ${B}$}}}{{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}^T)^{-1}`.
The derivation of this equation from eqns. [eqn:c1] and [eqn:c2] can be
found in `49 <#ref-Hess97>`__.

This first step does not set the real bond lengths to the prescribed
lengths, but the projection of the new bonds onto the old directions of
the bonds. To correct for the rotation of bond :math:`i`, the projection
of the bond, :math:`p_i`, on the old direction is set to

.. math::

   \label{eqn:m1a}
   p_i=\sqrt{2 d_i^2 - l_i^2},

 where :math:`l_i` is the bond length after the first projection. The
corrected positions are

.. math::

   \label{eqn:m1b}
   {\mbox{\boldmath ${r}$}}_{n+1}^*=({\mbox{\boldmath ${I}$}}-{{\mbox{\boldmath ${T}$}}}_n {{\mbox{\boldmath ${B}$}}}_n){\mbox{\boldmath ${r}$}}_{n+1} + {{\mbox{\boldmath ${T}$}}}_n {\mbox{\boldmath ${p}$}}.

 This correction for rotational effects is actually an iterative
process, but during MD only one iteration is applied. The relative
constraint deviation after this procedure will be less than 0.0001 for
every constraint. In energy minimization, this might not be accurate
enough, so the number of iterations is equal to the order of the
expansion (see below).

Half of the CPU time goes to inverting the constraint coupling matrix
:math:`{{\mbox{\boldmath ${B}$}}}_n {{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}_n^T`,
which has to be done every time step. This :math:`K \times K` matrix has
:math:`1/m_{i_1} + 1/m_{i_2}` on the diagonal. The off-diagonal elements
are only non-zero when two bonds are connected, then the element is
:math:`\cos \phi /m_c`, where :math:`m_c` is the mass of the atom
connecting the two bonds and :math:`\phi` is the angle between the
bonds.

The matrix :math:`{{\mbox{\boldmath ${T}$}}}` is inverted through a
power expansion. A :math:`K \times K` matrix
:math:`{\mbox{\boldmath ${S}$}}` is introduced which is the inverse
square root of the diagonal of
:math:`{{\mbox{\boldmath ${B}$}}}_n {{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}_n^T`.
This matrix is used to convert the diagonal elements of the coupling
matrix to one:

.. math::

   \label{eqn:m2}
   \begin{array}{c}
   ({{\mbox{\boldmath ${B}$}}}_n {{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}_n^T)^{-1}
   = {{\mbox{\boldmath ${S}$}}}{{\mbox{\boldmath ${S}$}}}^{-1} ({{\mbox{\boldmath ${B}$}}}_n {{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}_n^T)^{-1} {{\mbox{\boldmath ${S}$}}}^{-1} {{\mbox{\boldmath ${S}$}}}\\[2mm]
   = {{\mbox{\boldmath ${S}$}}}({{\mbox{\boldmath ${S}$}}}{{\mbox{\boldmath ${B}$}}}_n {{{\mbox{\boldmath ${M}$}}}^{-1}}{{\mbox{\boldmath ${B}$}}}_n^T {{\mbox{\boldmath ${S}$}}})^{-1} {{\mbox{\boldmath ${S}$}}}=
     {{\mbox{\boldmath ${S}$}}}({\mbox{\boldmath ${I}$}} - {\mbox{\boldmath ${A}$}}_n)^{-1} {{\mbox{\boldmath ${S}$}}}\end{array}

 The matrix :math:`{\mbox{\boldmath ${A}$}}_n` is symmetric and sparse
and has zeros on the diagonal. Thus a simple trick can be used to
calculate the inverse:

.. math::

   \label{eqn:m3}
   ({\mbox{\boldmath ${I}$}}-{\mbox{\boldmath ${A}$}}_n)^{-1}= 
           {\mbox{\boldmath ${I}$}} + {\mbox{\boldmath ${A}$}}_n + {\mbox{\boldmath ${A}$}}_n^2 + {\mbox{\boldmath ${A}$}}_n^3 + \ldots

This inversion method is only valid if the absolute values of all the
eigenvalues of :math:`{\mbox{\boldmath ${A}$}}_n` are smaller than one.
In molecules with only bond constraints, the connectivity is so low that
this will always be true, even if ring structures are present. Problems
can arise in angle-constrained molecules. By constraining angles with
additional distance constraints, multiple small ring structures are
introduced. This gives a high connectivity, leading to large
eigenvalues. Therefore LINCS should NOT be used with coupled
angle-constraints.

For molecules with all bonds constrained the eigenvalues of :math:`A`
are around 0.4. This means that with each additional order in the
expansion eqn. [eqn:m3] the deviations decrease by a factor 0.4. But for
relatively isolated triangles of constraints the largest eigenvalue is
around 0.7. Such triangles can occur when removing hydrogen angle
vibrations with an additional angle constraint in alcohol groups or when
constraining water molecules with LINCS, for instance with flexible
constraints. The constraints in such triangles converge twice as slow as
the other constraints. Therefore, starting with GROMACS 4, additional
terms are added to the expansion for such triangles

.. math::

   \label{eqn:m3_ang}
   ({\mbox{\boldmath ${I}$}}-{\mbox{\boldmath ${A}$}}_n)^{-1} \approx
           {\mbox{\boldmath ${I}$}} + {\mbox{\boldmath ${A}$}}_n + \ldots + {\mbox{\boldmath ${A}$}}_n^{N_i} +
           \left({\mbox{\boldmath ${A}$}}^*_n + \ldots + {{\mbox{\boldmath ${A}$}}_n^*}^{N_i} \right) {\mbox{\boldmath ${A}$}}_n^{N_i}

 where :math:`N_i` is the normal order of the expansion and
:math:`{\mbox{\boldmath ${A}$}}^*` only contains the elements of
:math:`{\mbox{\boldmath ${A}$}}` that couple constraints within rigid
triangles, all other elements are zero. In this manner, the accuracy of
angle constraints comes close to that of the other constraints, while
the series of matrix vector multiplications required for determining the
expansion only needs to be extended for a few constraint couplings. This
procedure is described in the P-LINCS paper\ `50 <#ref-Hess2008a>`__.

The LINCS Parameters
^^^^^^^^^^^^^^^^^^^^

The accuracy of LINCS depends on the number of matrices used in the
expansion eqn. [eqn:m3]. For MD calculations a fourth order expansion is
enough. For Brownian dynamics with large time steps an eighth order
expansion may be necessary. The order is a parameter in the .mdp file.
The implementation of LINCS is done in such a way that the algorithm
will never crash. Even when it is impossible to to reset the constraints
LINCS will generate a conformation which fulfills the constraints as
well as possible. However, LINCS will generate a warning when in one
step a bond rotates over more than a predefined angle. This angle is set
by the user in the .mdp file.

Simulated Annealing
-------------------

The well known simulated annealing (SA) protocol is supported in
GROMACS, and you can even couple multiple groups of atoms separately
with an arbitrary number of reference temperatures that change during
the simulation. The annealing is implemented by simply changing the
current reference temperature for each group in the temperature
coupling, so the actual relaxation and coupling properties depends on
the type of thermostat you use and how hard you are coupling it. Since
we are changing the reference temperature it is important to remember
that the system will NOT instantaneously reach this value - you need to
allow for the inherent relaxation time in the coupling algorithm too. If
you are changing the annealing reference temperature faster than the
temperature relaxation you will probably end up with a crash when the
difference becomes too large.

The annealing protocol is specified as a series of corresponding times
and reference temperatures for each group, and you can also choose
whether you only want a single sequence (after which the temperature
will be coupled to the last reference value), or if the annealing should
be periodic and restart at the first reference point once the sequence
is completed. You can mix and match both types of annealing and
non-annealed groups in your simulation.

Stochastic Dynamics
-------------------

Stochastic or velocity Langevin dynamics adds a friction and a noise
term to Newton’s equations of motion, as

.. math::

   \label{SDeq}
   m_i {{\mbox{d}}^2 {\mbox{\boldmath ${r}$}}_i \over {\mbox{d}}t^2} =
   - m_i \gamma_i {{\mbox{d}}{\mbox{\boldmath ${r}$}}_i \over {\mbox{d}}t} + {\mbox{\boldmath ${F}$}}_i({\mbox{\boldmath ${r}$}}) + {\stackrel{\circ}{{\mbox{\boldmath ${r}$}}}}_i,

 where :math:`\gamma_i` is the friction constant :math:`[1/\mbox{ps}]`
and :math:`{\stackrel{\circ}{{\mbox{\boldmath ${r}$}}}}_i\!\!(t)` is a
noise process with
:math:`\langle {\stackrel{\circ}{r}}_i\!\!(t) {\stackrel{\circ}{r}}_j\!\!(t+s) \rangle = 
    2 m_i \gamma_i k_B T \delta(s) \delta_{ij}`. When :math:`1/\gamma_i`
is large compared to the time scales present in the system, one could
see stochastic dynamics as molecular dynamics with stochastic
temperature-coupling. But any processes that take longer than
:math:`1/\gamma_i`, e.g. hydrodynamics, will be dampened. Since each
degree of freedom is coupled independently to a heat bath, equilibration
of fast modes occurs rapidly. For simulating a system in vacuum there is
the additional advantage that there is no accumulation of errors for the
overall translational and rotational degrees of freedom. When
:math:`1/\gamma_i` is small compared to the time scales present in the
system, the dynamics will be completely different from MD, but the
sampling is still correct.

In GROMACS there is one simple and efficient implementation. Its
accuracy is equivalent to the normal MD leap-frog and Velocity Verlet
integrator. It is nearly identical to the common way of discretizing the
Langevin equation, but the friction and velocity term are applied in an
impulse fashion \ `51 <#ref-Goga2012>`__. It can be described as:

.. math::

   \begin{aligned}
   \label{eqn:sd_int1}
   {\mbox{\boldmath ${v}$}}'  &~=~&   {\mbox{\boldmath ${v}$}}(t-{{\frac{1}{2}}{{\Delta t}}}) + \frac{1}{m}{\mbox{\boldmath ${F}$}}(t){{\Delta t}}\\
   \Delta{\mbox{\boldmath ${v}$}}     &~=~&   -\alpha \, {\mbox{\boldmath ${v}$}}'(t+{{\frac{1}{2}}{{\Delta t}}}) + \sqrt{\frac{k_B T}{m}(1 - \alpha^2)} \, {{\mbox{\boldmath ${r}$}}^G}_i \\
   {\mbox{\boldmath ${r}$}}(t+{{\Delta t}})   &~=~&   {\mbox{\boldmath ${r}$}}(t)+\left({\mbox{\boldmath ${v}$}}' +\frac{1}{2}\Delta {\mbox{\boldmath ${v}$}}\right){{\Delta t}}\label{eqn:sd1_x_upd}\\
   {\mbox{\boldmath ${v}$}}(t+{{\frac{1}{2}}{{\Delta t}}})  &~=~&   {\mbox{\boldmath ${v}$}}' + \Delta {\mbox{\boldmath ${v}$}} \\
   \alpha &~=~& 1 - e^{-\gamma {{\Delta t}}}\end{aligned}

 where :math:`{{\mbox{\boldmath ${r}$}}^G}_i` is Gaussian distributed
noise with :math:`\mu = 0`, :math:`\sigma = 1`. The velocity is first
updated a full time step without friction and noise to get
:math:`{\mbox{\boldmath ${v}$}}'`, identical to the normal update in
leap-frog. The friction and noise are then applied as an impulse at step
:math:`t+{{\Delta t}}`. The advantage of this scheme is that the
velocity-dependent terms act at the full time step, which makes the
correct integration of forces that depend on both coordinates and
velocities, such as constraints and dissipative particle dynamics (DPD,
not implented yet), straightforward. With constraints, the coordinate
update eqn. [eqn:sd1\_x\_upd] is split into a normal leap-frog update
and a :math:`\Delta {\mbox{\boldmath ${v}$}}`. After both of these
updates the constraints are applied to coordinates and velocities.

When using SD as a thermostat, an appropriate value for :math:`\gamma`
is e.g. 0.5 ps\ :math:`^{-1}`, since this results in a friction that is
lower than the internal friction of water, while it still provides
efficient thermostatting.

Brownian Dynamics
-----------------

In the limit of high friction, stochastic dynamics reduces to Brownian
dynamics, also called position Langevin dynamics. This applies to
over-damped systems, *i.e.* systems in which the inertia effects are
negligible. The equation is

.. math:: {{\mbox{d}}{\mbox{\boldmath ${r}$}}_i \over {\mbox{d}}t} = \frac{1}{\gamma_i} {\mbox{\boldmath ${F}$}}_i({\mbox{\boldmath ${r}$}}) + {\stackrel{\circ}{{\mbox{\boldmath ${r}$}}}}_i

 where :math:`\gamma_i` is the friction coefficient
:math:`[\mbox{amu/ps}]` and
:math:`{\stackrel{\circ}{{\mbox{\boldmath ${r}$}}}}_i\!\!(t)` is a noise
process with
:math:`\langle {\stackrel{\circ}{r}}_i\!\!(t) {\stackrel{\circ}{r}}_j\!\!(t+s) \rangle = 
    2 \delta(s) \delta_{ij} k_B T / \gamma_i`. In GROMACS the equations
are integrated with a simple, explicit scheme

.. math::

   {\mbox{\boldmath ${r}$}}_i(t+\Delta t) = {\mbox{\boldmath ${r}$}}_i(t) +
           {\Delta t \over \gamma_i} {\mbox{\boldmath ${F}$}}_i({\mbox{\boldmath ${r}$}}(t)) 
           + \sqrt{2 k_B T {\Delta t \over \gamma_i}}\, {{\mbox{\boldmath ${r}$}}^G}_i,

 where :math:`{{\mbox{\boldmath ${r}$}}^G}_i` is Gaussian distributed
noise with :math:`\mu = 0`, :math:`\sigma = 1`. The friction
coefficients :math:`\gamma_i` can be chosen the same for all particles
or as :math:`\gamma_i = m_i\,\gamma_i`, where the friction constants
:math:`\gamma_i` can be different for different groups of atoms. Because
the system is assumed to be over-damped, large timesteps can be used.
LINCS should be used for the constraints since SHAKE will not converge
for large atomic displacements. BD is an option of the mdrun program.

Energy Minimization
-------------------

Energy minimization in GROMACS can be done using steepest descent,
conjugate gradients, or l-bfgs (limited-memory
Broyden-Fletcher-Goldfarb-Shanno quasi-Newtonian minimizer...we prefer
the abbreviation). EM is just an option of the mdrun program.

Steepest Descent
~~~~~~~~~~~~~~~~

Although steepest descent is certainly not the most efficient algorithm
for searching, it is robust and easy to implement.

We define the vector :math:`{\mbox{\boldmath ${r}$}}` as the vector of
all :math:`3N` coordinates. Initially a maximum displacement :math:`h_0`
(*e.g.* 0.01 nm) must be given.

| First the forces :math:`{\mbox{\boldmath ${F}$}}` and potential energy
  are calculated. New positions are calculated by

  .. math:: {\mbox{\boldmath ${r}$}}_{n+1} =  {\mbox{\boldmath ${r}$}}_n + \frac{{\mbox{\boldmath ${F}$}}_n}{\max (|{\mbox{\boldmath ${F}$}}_n|)} h_n,

   where :math:`h_n` is the maximum displacement and
  :math:`{\mbox{\boldmath ${F}$}}_n` is the force, or the negative
  gradient of the potential :math:`V`. The notation :math:`\max
  (|{\mbox{\boldmath ${F}$}}_n|)` means the largest scalar force on any
  atom. The forces and energy are again computed for the new positions
| If (:math:`V_{n+1} < V_n`) the new positions are accepted and
  :math:`h_{n+1} = 1.2
  h_n`.
| If (:math:`V_{n+1} \geq V_n`) the new positions are rejected and
  :math:`h_n = 0.2 h_n`.

The algorithm stops when either a user-specified number of force
evaluations has been performed (*e.g.* 100), or when the maximum of the
absolute values of the force (gradient) components is smaller than a
specified value :math:`\epsilon`. Since force truncation produces some
noise in the energy evaluation, the stopping criterion should not be
made too tight to avoid endless iterations. A reasonable value for
:math:`\epsilon` can be estimated from the root mean square force
:math:`f` a harmonic oscillator would exhibit at a temperature
:math:`T`. This value is

.. math:: f = 2 \pi \nu \sqrt{ 2mkT},

 where :math:`\nu` is the oscillator frequency, :math:`m` the (reduced)
mass, and :math:`k` Boltzmann’s constant. For a weak oscillator with a
wave number of 100 cm\ :math:`^{-1}` and a mass of 10 atomic units, at a
temperature of 1 K, :math:`f=7.7` kJ mol\ :math:`^{-1}` nm:math:`^{-1}`.
A value for :math:`\epsilon` between 1 and 10 is acceptable.

Conjugate Gradient
~~~~~~~~~~~~~~~~~~

Conjugate gradient is slower than steepest descent in the early stages
of the minimization, but becomes more efficient closer to the energy
minimum. The parameters and stop criterion are the same as for steepest
descent. In GROMACS conjugate gradient can not be used with constraints,
including the SETTLE algorithm for water \ `47 <#ref-Miyamoto92>`__, as
this has not been implemented. If water is present it must be of a
flexible model, which can be specified in the .mdp file by define =
-DFLEXIBLE.

This is not really a restriction, since the accuracy of conjugate
gradient is only required for minimization prior to a normal-mode
analysis, which cannot be performed with constraints. For most other
purposes steepest descent is efficient enough.

L-BFGS
~~~~~~

The original BFGS algorithm works by successively creating better
approximations of the inverse Hessian matrix, and moving the system to
the currently estimated minimum. The memory requirements for this are
proportional to the square of the number of particles, so it is not
practical for large systems like biomolecules. Instead, we use the
L-BFGS algorithm of Nocedal \ `52 <#ref-Byrd95a>`__,
`53 <#ref-Zhu97a>`__, which approximates the inverse Hessian by a fixed
number of corrections from previous steps. This sliding-window technique
is almost as efficient as the original method, but the memory
requirements are much lower - proportional to the number of particles
multiplied with the correction steps. In practice we have found it to
converge faster than conjugate gradients, but due to the correction
steps it is not yet parallelized. It is also noteworthy that switched or
shifted interactions usually improve the convergence, since sharp
cut-offs mean the potential function at the current coordinates is
slightly different from the previous steps used to build the inverse
Hessian approximation.

Normal-Mode Analysis
--------------------

Normal-mode analysis \ `54 <#ref-Levitt83>`__\ `56 <#ref-BBrooks83b>`__
can be performed using GROMACS, by diagonalization of the mass-weighted
Hessian :math:`H`:

.. math::

   \begin{aligned}
   R^T M^{-1/2} H M^{-1/2} R   &=& \mbox{diag}(\lambda_1,\ldots,\lambda_{3N})
   \\
   \lambda_i &=& (2 \pi \omega_i)^2\end{aligned}

 where :math:`M` contains the atomic masses, :math:`R` is a matrix that
contains the eigenvectors as columns, :math:`\lambda_i` are the
eigenvalues and :math:`\omega_i` are the corresponding frequencies.

First the Hessian matrix, which is a :math:`3N \times 3N` matrix where
:math:`N` is the number of atoms, needs to be calculated:

.. math::

   \begin{aligned}
   H_{ij}  &=&     \frac{\partial^2 V}{\partial x_i \partial x_j}\end{aligned}

 where :math:`x_i` and :math:`x_j` denote the atomic x, y or z
coordinates. In practice, this equation is not used, but the Hessian is
calculated numerically from the force as:

.. math::

   \begin{aligned}
   H_{ij} &=& -
     \frac{f_i({\bf x}+h{\bf e}_j) - f_i({\bf x}-h{\bf e}_j)}{2h}
   \\
   f_i     &=& - \frac{\partial V}{\partial x_i}\end{aligned}

 where :math:`{\bf e}_j` is the unit vector in direction :math:`j`. It
should be noted that for a usual normal-mode calculation, it is
necessary to completely minimize the energy prior to computation of the
Hessian. The tolerance required depends on the type of system, but a
rough indication is 0.001 kJ mol\ :math:`^{-1}`. Minimization should be
done with conjugate gradients or L-BFGS in double precision.

A number of GROMACS programs are involved in these calculations. First,
the energy should be minimized using mdrun. Then, mdrun computes the
Hessian. **Note** that for generating the run input file, one should use
the minimized conformation from the full precision trajectory file, as
the structure file is not accurate enough. gmx nmeig does the
diagonalization and the sorting of the normal modes according to their
frequencies. Both mdrun and gmx nmeig should be run in double precision.
The normal modes can be analyzed with the program gmx anaeig. Ensembles
of structures at any temperature and for any subset of normal modes can
be generated with gmx nmens. An overview of normal-mode analysis and the
related principal component analysis (see sec. [sec:covanal]) can be
found in \ `57 <#ref-Hayward95b>`__.

Free energy calculations
------------------------

Slow-growth methods
~~~~~~~~~~~~~~~~~~~

Free energy calculations can be performed in GROMACS using a number of
methods, including “slow-growth.” An example problem might be
calculating the difference in free energy of binding of an inhibitor
**I** to an enzyme **E** and to a mutated enzyme
**E\ :math:`^{\prime}`**. It is not feasible with computer simulations
to perform a docking calculation for such a large complex, or even
releasing the inhibitor from the enzyme in a reasonable amount of
computer time with reasonable accuracy. However, if we consider the free
energy cycle in Fig. [fig:free]A we can write:

.. math::

   \Delta G_1 - \Delta G_2 =       \Delta G_3 - \Delta G_4
   \label{eqn:ddg}

 If we are interested in the left-hand term we can equally well compute
the right-hand term.

|Free energy cycles. **A:** to calculate :math:`\Delta G_{12}`, the free
energy difference between the binding of inhibitor **I** to enzymes
**E** respectively **E\ :math:`^{\prime}`**. **B:** to calculate
:math:`\Delta G_{12}`, the free energy difference for binding of
inhibitors **I** respectively **I\ :math:`^{\prime}`** to enzyme
**E**.| |Free energy cycles. **A:** to calculate :math:`\Delta G_{12}`,
the free energy difference between the binding of inhibitor **I** to
enzymes **E** respectively **E\ :math:`^{\prime}`**. **B:** to calculate
:math:`\Delta G_{12}`, the free energy difference for binding of
inhibitors **I** respectively **I\ :math:`^{\prime}`** to enzyme **E**.|

If we want to compute the difference in free energy of binding of two
inhibitors **I** and **I\ :math:`^{\prime}`** to an enzyme **E**
(Fig. [fig:free]B) we can again use eqn. [eqn:ddg] to compute the
desired property.

Free energy differences between two molecular species can be calculated
in GROMACS using the “slow-growth” method. Such free energy differences
between different molecular species are physically meaningless, but they
can be used to obtain meaningful quantities employing a thermodynamic
cycle. The method requires a simulation during which the Hamiltonian of
the system changes slowly from that describing one system (A) to that
describing the other system (B). The change must be so slow that the
system remains in equilibrium during the process; if that requirement is
fulfilled, the change is reversible and a slow-growth simulation from B
to A will yield the same results (but with a different sign) as a
slow-growth simulation from A to B. This is a useful check, but the user
should be aware of the danger that equality of forward and backward
growth results does not guarantee correctness of the results.

The required modification of the Hamiltonian :math:`H` is realized by
making :math:`H` a function of a *coupling parameter* :math:`\lambda:
H=H(p,q;\lambda)` in such a way that :math:`\lambda=0` describes system
A and :math:`\lambda=1` describes system B:

.. math:: H(p,q;0)=H{^{\mathrm{A}}}(p,q);~~~~ H(p,q;1)=H{^{\mathrm{B}}}(p,q).

 In GROMACS, the functional form of the :math:`\lambda`-dependence is
different for the various force-field contributions and is described in
section sec. [sec:feia].

The Helmholtz free energy :math:`A` is related to the partition function
:math:`Q` of an :math:`N,V,T` ensemble, which is assumed to be the
equilibrium ensemble generated by a MD simulation at constant volume and
temperature. The generally more useful Gibbs free energy :math:`G` is
related to the partition function :math:`\Delta` of an :math:`N,p,T`
ensemble, which is assumed to be the equilibrium ensemble generated by a
MD simulation at constant pressure and temperature:

.. math::

   \begin{aligned}
    A(\lambda) &=&  -k_BT \ln Q \\
    Q &=& c \int\!\!\int \exp[-\beta H(p,q;\lambda)]\,dp\,dq \\
    G(\lambda) &=&  -k_BT \ln \Delta \\
    \Delta &=& c \int\!\!\int\!\!\int \exp[-\beta H(p,q;\lambda) -\beta
   pV]\,dp\,dq\,dV \\
   G &=& A + pV, \end{aligned}

 where :math:`\beta = 1/(k_BT)` and :math:`c = (N! h^{3N})^{-1}`. These
integrals over phase space cannot be evaluated from a simulation, but it
is possible to evaluate the derivative with respect to :math:`\lambda`
as an ensemble average:

.. math::

   \frac{dA}{d\lambda} =  \frac{\int\!\!\int (\partial H/ \partial
   \lambda) \exp[-\beta H(p,q;\lambda)]\,dp\,dq}{\int\!\!\int \exp[-\beta
   H(p,q;\lambda)]\,dp\,dq} = 
   \left\langle \frac{\partial H}{\partial \lambda} \right\rangle_{NVT;\lambda},

 with a similar relation for :math:`dG/d\lambda` in the :math:`N,p,T`
ensemble. The difference in free energy between A and B can be found by
integrating the derivative over :math:`\lambda`:

.. math::

   \begin{aligned}
     A{^{\mathrm{B}}}(V,T)-A{^{\mathrm{A}}}(V,T) &=& \int_0^1 \left\langle \frac{\partial
   H}{\partial \lambda} \right\rangle_{NVT;\lambda} \,d\lambda 
   \label{eq:delA} \\
    G{^{\mathrm{B}}}(p,T)-G{^{\mathrm{A}}}(p,T) &=& \int_0^1 \left\langle \frac{\partial
   H}{\partial \lambda} \right\rangle_{NpT;\lambda} \,d\lambda.
   \label{eq:delG}\end{aligned}

 If one wishes to evaluate
:math:`G{^{\mathrm{B}}}(p,T)-G{^{\mathrm{A}}}(p,T)`, the natural choice
is a constant-pressure simulation. However, this quantity can also be
obtained from a slow-growth simulation at constant volume, starting with
system A at pressure :math:`p` and volume :math:`V` and ending with
system B at pressure :math:`p_B`, by applying the following small (but,
in principle, exact) correction:

.. math::

   G{^{\mathrm{B}}}(p)-G{^{\mathrm{A}}}(p) =
   A{^{\mathrm{B}}}(V)-A{^{\mathrm{A}}}(V) - \int_p^{p{^{\mathrm{B}}}}[V{^{\mathrm{B}}}(p')-V]\,dp'

 Here we omitted the constant :math:`T` from the notation. This
correction is roughly equal to
:math:`-\frac{1}{2} (p{^{\mathrm{B}}}-p)\Delta V=(\Delta V)^2/(2
\kappa V)`, where :math:`\Delta V` is the volume change at :math:`p` and
:math:`\kappa` is the isothermal compressibility. This is usually small;
for example, the growth of a water molecule from nothing in a bath of
1000 water molecules at constant volume would produce an additional
pressure of as much as 22 bar, but a correction to the Helmholtz free
energy of just -1 kJ mol\ :math:`^{-1}`. In Cartesian coordinates, the
kinetic energy term in the Hamiltonian depends only on the momenta, and
can be separately integrated and, in fact, removed from the equations.
When masses do not change, there is no contribution from the kinetic
energy at all; otherwise the integrated contribution to the free energy
is :math:`-\frac{3}{2} k_BT \ln
(m{^{\mathrm{B}}}/m{^{\mathrm{A}}})`. **Note** that this is only true in
the absence of constraints.

Thermodynamic integration
~~~~~~~~~~~~~~~~~~~~~~~~~

GROMACS offers the possibility to integrate eq. [eq:delA] or eq.
[eq:delG] in one simulation over the full range from A to B. However, if
the change is large and insufficient sampling can be expected, the user
may prefer to determine the value of :math:`\langle
dG/d\lambda \rangle` accurately at a number of well-chosen intermediate
values of :math:`\lambda`. This can easily be done by setting the
stepsize delta\_lambda to zero. Each simulation can be equilibrated
first, and a proper error estimate can be made for each value of
:math:`dG/d\lambda` from the fluctuation of :math:`\partial H/\partial
\lambda`. The total free energy change is then determined afterward by
an appropriate numerical integration procedure.

GROMACS now also supports the use of Bennett’s Acceptance
Ratio \ `58 <#ref-Bennett1976>`__ for calculating values of
:math:`\Delta`\ G for transformations from state A to state B using the
program gmx bar. The same data can also be used to calculate free
energies using MBAR \ `59 <#ref-Shirts2008>`__, though the analysis
currently requires external tools from the external pymbar package, at
https://SimTK.org/home/pymbar.

The :math:`\lambda`-dependence for the force-field contributions is
described in detail in section sec. [sec:feia].

Replica exchange
----------------

Replica exchange molecular dynamics (REMD) is a method that can be used
to speed up the sampling of any type of simulation, especially if
conformations are separated by relatively high energy barriers. It
involves simulating multiple replicas of the same system at different
temperatures and randomly exchanging the complete state of two replicas
at regular intervals with the probability:

.. math::

   P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
   \left(\frac{1}{k_B T_1} - \frac{1}{k_B T_2}\right)(U_1 - U_2)
    \right] \right)

 where :math:`T_1` and :math:`T_2` are the reference temperatures and
:math:`U_1` and :math:`U_2` are the instantaneous potential energies of
replicas 1 and 2 respectively. After exchange the velocities are scaled
by :math:`(T_1/T_2)^{\pm0.5}` and a neighbor search is performed the
next step. This combines the fast sampling and frequent barrier-crossing
of the highest temperature with correct Boltzmann sampling at all the
different temperatures \ `60 <#ref-Hukushima96a>`__,
`61 <#ref-Sugita99>`__. We only attempt exchanges for neighboring
temperatures as the probability decreases very rapidly with the
temperature difference. One should not attempt exchanges for all
possible pairs in one step. If, for instance, replicas 1 and 2 would
exchange, the chance of exchange for replicas 2 and 3 not only depends
on the energies of replicas 2 and 3, but also on the energy of replica
1. In GROMACS this is solved by attempting exchange for all “odd” pairs
on “odd” attempts and for all “even” pairs on “even” attempts. If we
have four replicas: 0, 1, 2 and 3, ordered in temperature and we attempt
exchange every 1000 steps, pairs 0-1 and 2-3 will be tried at steps
1000, 3000 etc. and pair 1-2 at steps 2000, 4000 etc.

How should one choose the temperatures? The energy difference can be
written as:

.. math:: U_1 - U_2 =  N_{df} \frac{c}{2} k_B (T_1 - T_2)

 where :math:`N_{df}` is the total number of degrees of freedom of one
replica and :math:`c` is 1 for harmonic potentials and around 2 for
protein/water systems. If :math:`T_2 = (1+\epsilon) T_1` the probability
becomes:

.. math::

   P(1 \leftrightarrow 2)
     = \exp\left( -\frac{\epsilon^2 c\,N_{df}}{2 (1+\epsilon)} \right)
   \approx \exp\left(-\epsilon^2 \frac{c}{2} N_{df} \right)

 Thus for a probability of :math:`e^{-2}\approx 0.135` one obtains
:math:`\epsilon \approx 2/\sqrt{c\,N_{df}}`. With all bonds constrained
one has :math:`N_{df} \approx 2\, N_{atoms}` and thus for :math:`c` = 2
one should choose :math:`\epsilon` as :math:`1/\sqrt{N_{atoms}}`.
However there is one problem when using pressure coupling. The density
at higher temperatures will decrease, leading to higher
energy \ `62 <#ref-Seibert2005a>`__, which should be taken into account.
The GROMACS website features a so-called “REMD calculator,” that lets
you type in the temperature range and the number of atoms, and based on
that proposes a set of temperatures.

An extension to the REMD for the isobaric-isothermal ensemble was
proposed by Okabe *et al.* `63 <#ref-Okabe2001a>`__. In this work the
exchange probability is modified to:

.. math::

   P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
   \left(\frac{1}{k_B T_1} - \frac{1}{k_B T_2}\right)(U_1 - U_2) +
   \left(\frac{P_1}{k_B T_1} - \frac{P_2}{k_B T_2}\right)\left(V_1-V_2\right)
    \right] \right)

 where :math:`P_1` and :math:`P_2` are the respective reference
pressures and :math:`V_1` and :math:`V_2` are the respective
instantaneous volumes in the simulations. In most cases the differences
in volume are so small that the second term is negligible. It only plays
a role when the difference between :math:`P_1` and :math:`P_2` is large
or in phase transitions.

Hamiltonian replica exchange is also supported in GROMACS. In
Hamiltonian replica exchange, each replica has a different Hamiltonian,
defined by the free energy pathway specified for the simulation. The
exchange probability to maintain the correct ensemble probabilities is:

.. math::

   P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
       \left(\frac{1}{k_B T} - \frac{1}{k_B T}\right)((U_1(x_2) - U_1(x_1)) + (U_2(x_1) - U_2(x_2)))
   \right]
   \right)

 The separate Hamiltonians are defined by the free energy functionality
of GROMACS, with swaps made between the different values of
:math:`\lambda` defined in the mdp file.

Hamiltonian and temperature replica exchange can also be performed
simultaneously, using the acceptance criteria:

.. math::

   P(1 \leftrightarrow 2)=\min\left(1,\exp\left[
   \left(\frac{1}{k_B T} - \right)(\frac{U_1(x_2) - U_1(x_1)}{k_B T_1} + \frac{U_2(x_1) - U_2(x_2)}{k_B T_2})
    \right] \right)

Gibbs sampling replica exchange has also been implemented in
GROMACS `64 <#ref-Chodera2011>`__. In Gibbs sampling replica exchange,
all possible pairs are tested for exchange, allowing swaps between
replicas that are not neighbors.

Gibbs sampling replica exchange requires no additional potential energy
calculations. However there is an additional communication cost in Gibbs
sampling replica exchange, as for some permutations, more than one round
of swaps must take place. In some cases, this extra communication cost
might affect the efficiency.

All replica exchange variants are options of the mdrun program. It will
only work when MPI is installed, due to the inherent parallelism in the
algorithm. For efficiency each replica can run on a separate rank. See
the manual page of mdrun on how to use these multinode features.

Essential Dynamics sampling
---------------------------

The results from Essential Dynamics (see sec. [sec:covanal]) of a
protein can be used to guide MD simulations. The idea is that from an
initial MD simulation (or from other sources) a definition of the
collective fluctuations with largest amplitude is obtained. The position
along one or more of these collective modes can be constrained in a
(second) MD simulation in a number of ways for several purposes. For
example, the position along a certain mode may be kept fixed to monitor
the average force (free-energy gradient) on that coordinate in that
position. Another application is to enhance sampling efficiency with
respect to usual MD `65 <#ref-Degroot96a>`__, `66 <#ref-Degroot96b>`__.
In this case, the system is encouraged to sample its available
configuration space more systematically than in a diffusion-like path
that proteins usually take.

Another possibility to enhance sampling is flooding. Here a flooding
potential is added to certain (collective) degrees of freedom to expel
the system out of a region of phase space `67 <#ref-Lange2006a>`__.

The procedure for essential dynamics sampling or flooding is as follows.
First, the eigenvectors and eigenvalues need to be determined using
covariance analysis (gmx covar) or normal-mode analysis (gmx nmeig).
Then, this information is fed into make\_edi, which has many options for
selecting vectors and setting parameters, see gmx make\_edi -h. The
generated edi input file is then passed to mdrun.

Expanded Ensemble
-----------------

In an expanded ensemble simulation \ `68 <#ref-Lyubartsev1992>`__, both
the coordinates and the thermodynamic ensemble are treated as
configuration variables that can be sampled over. The probability of any
given state can be written as:

.. math:: P(\vec{x},k) \propto \exp\left(-\beta_k U_k + g_k\right),

 where :math:`\beta_k = \frac{1}{k_B T_k}` is the :math:`\beta`
corresponding to the :math:`k`\ th thermodynamic state, and :math:`g_k`
is a user-specified weight factor corresponding to the :math:`k`\ th
state. This space is therefore a *mixed*, *generalized*, or *expanded*
ensemble which samples from multiple thermodynamic ensembles
simultaneously. :math:`g_k` is chosen to give a specific weighting of
each subensemble in the expanded ensemble, and can either be fixed, or
determined by an iterative procedure. The set of :math:`g_k` is
frequently chosen to give each thermodynamic ensemble equal probability,
in which case :math:`g_k` is equal to the free energy in non-dimensional
units, but they can be set to arbitrary values as desired. Several
different algorithms can be used to equilibrate these weights, described
in the mdp option listings.

In GROMACS, this space is sampled by alternating sampling in the
:math:`k` and :math:`\vec{x}` directions. Sampling in the
:math:`\vec{x}` direction is done by standard molecular dynamics
sampling; sampling between the different thermodynamics states is done
by Monte Carlo, with several different Monte Carlo moves supported. The
:math:`k` states can be defined by different temperatures, or choices of
the free energy :math:`\lambda` variable, or both. Expanded ensemble
simulations thus represent a serialization of the replica exchange
formalism, allowing a single simulation to explore many thermodynamic
states.

Parallelization
---------------

The CPU time required for a simulation can be reduced by running the
simulation in parallel over more than one core. Ideally, one would want
to have linear scaling: running on :math:`N` cores makes the simulation
:math:`N` times faster. In practice this can only be achieved for a
small number of cores. The scaling will depend a lot on the algorithms
used. Also, different algorithms can have different restrictions on the
interaction ranges between atoms.

Domain decomposition
--------------------

Since most interactions in molecular simulations are local, domain
decomposition is a natural way to decompose the system. In domain
decomposition, a spatial domain is assigned to each rank, which will
then integrate the equations of motion for the particles that currently
reside in its local domain. With domain decomposition, there are two
choices that have to be made: the division of the unit cell into domains
and the assignment of the forces to domains. Most molecular simulation
packages use the half-shell method for assigning the forces. But there
are two methods that always require less communication: the eighth
shell \ `69 <#ref-Liem1991>`__ and the midpoint \ `70 <#ref-Shaw2006>`__
method. GROMACS currently uses the eighth shell method, but for certain
systems or hardware architectures it might be advantageous to use the
midpoint method. Therefore, we might implement the midpoint method in
the future. Most of the details of the domain decomposition can be found
in the GROMACS 4 paper \ `5 <#ref-Hess2008b>`__.

Coordinate and force communication
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the most general case of a triclinic unit cell, the space in divided
with a 1-, 2-, or 3-D grid in parallelepipeds that we call domain
decomposition cells. Each cell is assigned to a particle-particle rank.
The system is partitioned over the ranks at the beginning of each MD
step in which neighbor searching is performed. Since the neighbor
searching is based on charge groups, charge groups are also the units
for the domain decomposition. Charge groups are assigned to the cell
where their center of geometry resides. Before the forces can be
calculated, the coordinates from some neighboring cells need to be
communicated, and after the forces are calculated, the forces need to be
communicated in the other direction. The communication and force
assignment is based on zones that can cover one or multiple cells. An
example of a zone setup is shown in Fig. [fig:ddcells].

.. figure:: plots/dd-cells
   :alt:  A non-staggered domain decomposition grid of
   3\ :math:`\times`\ 2\ :math:`\times`\ 2 cells. Coordinates in zones 1
   to 7 are communicated to the corner cell that has its home particles
   in zone 0. :math:`r_c` is the cut-off radius. [fig:ddcells]
   :width: 6.00000cm

    A non-staggered domain decomposition grid of
   3\ :math:`\times`\ 2\ :math:`\times`\ 2 cells. Coordinates in zones 1
   to 7 are communicated to the corner cell that has its home particles
   in zone 0. :math:`r_c` is the cut-off radius. [fig:ddcells] 

The coordinates are communicated by moving data along the “negative”
direction in :math:`x`, :math:`y` or :math:`z` to the next neighbor.
This can be done in one or multiple pulses. In Fig. [fig:ddcells] two
pulses in :math:`x` are required, then one in :math:`y` and then one in
:math:`z`. The forces are communicated by reversing this procedure. See
the GROMACS 4 paper \ `5 <#ref-Hess2008b>`__ for details on determining
which non-bonded and bonded forces should be calculated on which rank.

Dynamic load balancing
~~~~~~~~~~~~~~~~~~~~~~

When different ranks have a different computational load (load
imbalance), all ranks will have to wait for the one that takes the most
time. One would like to avoid such a situation. Load imbalance can occur
due to four reasons:

-  inhomogeneous particle distribution

-  inhomogeneous interaction cost distribution (charged/uncharged,
   water/non-water due to GROMACS water innerloops)

-  statistical fluctuation (only with small particle numbers)

-  differences in communication time, due to network topology and/or
   other jobs on the machine interfering with our communication

So we need a dynamic load balancing algorithm where the volume of each
domain decomposition cell can be adjusted *independently*. To achieve
this, the 2- or 3-D domain decomposition grids need to be staggered.
Fig. [fig:ddtric] shows the most general case in 2-D. Due to the
staggering, one might require two distance checks for deciding if a
charge group needs to be communicated: a non-bonded distance and a
bonded distance check.

.. figure:: plots/dd-tric
   :alt:  The zones to communicate to the rank of zone 0, see the text
   for details. :math:`r_c` and :math:`r_b` are the non-bonded and
   bonded cut-off radii respectively, :math:`d` is an example of a
   distance between following, staggered boundaries of cells.
   [fig:ddtric]
   :width: 7.00000cm

    The zones to communicate to the rank of zone 0, see the text for
   details. :math:`r_c` and :math:`r_b` are the non-bonded and bonded
   cut-off radii respectively, :math:`d` is an example of a distance
   between following, staggered boundaries of cells. [fig:ddtric] 

By default, mdrun automatically turns on the dynamic load balancing
during a simulation when the total performance loss due to the force
calculation imbalance is 2% or more. **Note** that the reported force
load imbalance numbers might be higher, since the force calculation is
only part of work that needs to be done during an integration step. The
load imbalance is reported in the log file at log output steps and when
the -v option is used also on screen. The average load imbalance and the
total performance loss due to load imbalance are reported at the end of
the log file.

There is one important parameter for the dynamic load balancing, which
is the minimum allowed scaling. By default, each dimension of the domain
decomposition cell can scale down by at least a factor of 0.8. For 3-D
domain decomposition this allows cells to change their volume by about a
factor of 0.5, which should allow for compensation of a load imbalance
of 100%. The minimum allowed scaling can be changed with the -dds option
of mdrun.

The load imbalance is measured by timing a single region of the MD step
on each MPI rank. This region can not include MPI communication, as
timing of MPI calls does not allow separating wait due to imbalance from
actual communication. The domain volumes are then scaled, with
under-relaxation, inversely proportional with the measured time. This
procedure will decrease the load imbalance when the change in load in
the measured region correlates with the change in domain volume and the
load outside the measured region does not depend strongly on the domain
volume. In CPU-only simulations, the load is measured between the
coordinate and the force communication. In simulations with non-bonded
work on GPUs, we overlap communication and work on the CPU with
calculation on the GPU. Therefore we measure from the last communication
before the force calculation to when the CPU or GPU is finished,
whichever is last. When not using PME ranks, we subtract the time in PME
from the CPU time, as this includes MPI calls and the PME load is
independent of domain size. This generally works well, unless the
non-bonded load is low and there is imbalance in the bonded
interactions. Then two issues can arise. Dynamic load balancing can
increase the imbalance in update and constraints and with PME the
coordinate and force redistribution time can go up significantly.
Although dynamic load balancing can significantly improve performance in
cases where there is imbalance in the bonded interactions on the CPU,
there are many situations in which some domains continue decreasing in
size and the load imbalance increases and/or PME coordinate and force
redistribution cost increases significantly. As of version 2016.1, mdrun
disables the dynamic load balancing when measurement indicates that it
deteriorates performance. This means that in most cases the user will
get good performance with the default, automated dynamic load balancing
setting.

Constraints in parallel
~~~~~~~~~~~~~~~~~~~~~~~

Since with domain decomposition parts of molecules can reside on
different ranks, bond constraints can cross cell boundaries. Therefore a
parallel constraint algorithm is required. GROMACS uses the P-LINCS
algorithm \ `50 <#ref-Hess2008a>`__, which is the parallel version of
the LINCS algorithm \ `49 <#ref-Hess97>`__ (see [subsec:lincs]). The
P-LINCS procedure is illustrated in Fig. [fig:plincs]. When molecules
cross the cell boundaries, atoms in such molecules up to (lincs\_order +
1) bonds away are communicated over the cell boundaries. Then, the
normal LINCS algorithm can be applied to the local bonds plus the
communicated ones. After this procedure, the local bonds are correctly
constrained, even though the extra communicated ones are not. One
coordinate communication step is required for the initial LINCS step and
one for each iteration. Forces do not need to be communicated.

.. figure:: plots/par-lincs2
   :alt:  Example of the parallel setup of P-LINCS with one molecule
   split over three domain decomposition cells, using a matrix expansion
   order of 3. The top part shows which atom coordinates need to be
   communicated to which cells. The bottom parts show the local
   constraints (solid) and the non-local constraints (dashed) for each
   of the three cells. [fig:plincs]
   :width: 6.00000cm

    Example of the parallel setup of P-LINCS with one molecule split
   over three domain decomposition cells, using a matrix expansion order
   of 3. The top part shows which atom coordinates need to be
   communicated to which cells. The bottom parts show the local
   constraints (solid) and the non-local constraints (dashed) for each
   of the three cells. [fig:plincs] 

Interaction ranges
~~~~~~~~~~~~~~~~~~

Domain decomposition takes advantage of the locality of interactions.
This means that there will be limitations on the range of interactions.
By default, mdrun tries to find the optimal balance between interaction
range and efficiency. But it can happen that a simulation stops with an
error message about missing interactions, or that a simulation might run
slightly faster with shorter interaction ranges. A list of interaction
ranges and their default values is given in Table [tab:dd\_ranges].

In most cases the defaults of mdrun should not cause the simulation to
stop with an error message of missing interactions. The range for the
bonded interactions is determined from the distance between bonded
charge-groups in the starting configuration, with 10% added for
headroom. For the constraints, the value of :math:`r_{\mathrm{con}}` is
determined by taking the maximum distance that (lincs\_order + 1) bonds
can cover when they all connect at angles of 120 degrees. The actual
constraint communication is not limited by :math:`r_{\mathrm{con}}`, but
by the minimum cell size :math:`L_C`, which has the following lower
limit:

.. math:: L_C \geq \max(r_{\mathrm{mb}},r_{\mathrm{con}})

 Without dynamic load balancing the system is actually allowed to scale
beyond this limit when pressure scaling is used. **Note** that for
triclinic boxes, :math:`L_C` is not simply the box diagonal component
divided by the number of cells in that direction, rather it is the
shortest distance between the triclinic cells borders. For rhombic
dodecahedra this is a factor of :math:`\sqrt{3/2}` shorter along
:math:`x` and :math:`y`.

When :math:`r_{\mathrm{mb}} > r_c`, mdrun employs a smart algorithm to
reduce the communication. Simply communicating all charge groups within
:math:`r_{\mathrm{mb}}` would increase the amount of communication
enormously. Therefore only charge-groups that are connected by bonded
interactions to charge groups which are not locally present are
communicated. This leads to little extra communication, but also to a
slightly increased cost for the domain decomposition setup. In some
cases, *e.g.* coarse-grained simulations with a very short cut-off, one
might want to set :math:`r_{\mathrm{mb}}` by hand to reduce this cost.

Multiple-Program, Multiple-Data PME parallelization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Electrostatics interactions are long-range, therefore special algorithms
are used to avoid summation over many atom pairs. In GROMACS this is
usually PME (sec. [sec:pme]). Since with PME all particles interact with
each other, global communication is required. This will usually be the
limiting factor for scaling with domain decomposition. To reduce the
effect of this problem, we have come up with a Multiple-Program,
Multiple-Data approach \ `5 <#ref-Hess2008b>`__. Here, some ranks are
selected to do only the PME mesh calculation, while the other ranks,
called particle-particle (PP) ranks, do all the rest of the work. For
rectangular boxes the optimal PP to PME rank ratio is usually 3:1, for
rhombic dodecahedra usually 2:1. When the number of PME ranks is reduced
by a factor of 4, the number of communication calls is reduced by about
a factor of 16. Or put differently, we can now scale to 4 times more
ranks. In addition, for modern 4 or 8 core machines in a network, the
effective network bandwidth for PME is quadrupled, since only a quarter
of the cores will be using the network connection on each machine during
the PME calculations.

.. figure:: plots/mpmd-pme
   :alt:  Example of 8 ranks without (left) and with (right) MPMD. The
   PME communication (red arrows) is much higher on the left than on the
   right. For MPMD additional PP - PME coordinate and force
   communication (blue arrows) is required, but the total communication
   complexity is lower. [fig:mpmd\_pme]
   :width: 12.00000cm

    Example of 8 ranks without (left) and with (right) MPMD. The PME
   communication (red arrows) is much higher on the left than on the
   right. For MPMD additional PP - PME coordinate and force
   communication (blue arrows) is required, but the total communication
   complexity is lower. [fig:mpmd\_pme] 

mdrun will by default interleave the PP and PME ranks. If the ranks are
not number consecutively inside the machines, one might want to use
mdrun -ddorder pp\_pme. For machines with a real 3-D torus and proper
communication software that assigns the ranks accordingly one should use
mdrun -ddorder cartesian.

To optimize the performance one should usually set up the cut-offs and
the PME grid such that the PME load is 25 to 33% of the total
calculation load. grompp will print an estimate for this load at the end
and also mdrun calculates the same estimate to determine the optimal
number of PME ranks to use. For high parallelization it might be
worthwhile to optimize the PME load with the mdp settings and/or the
number of PME ranks with the -npme option of mdrun. For changing the
electrostatics settings it is useful to know the accuracy of the
electrostatics remains nearly constant when the Coulomb cut-off and the
PME grid spacing are scaled by the same factor. **Note** that it is
usually better to overestimate than to underestimate the number of PME
ranks, since the number of PME ranks is smaller than the number of PP
ranks, which leads to less total waiting time.

The PME domain decomposition can be 1-D or 2-D along the :math:`x`
and/or :math:`y` axis. 2-D decomposition is also known as pencil
decomposition because of the shape of the domains at high
parallelization. 1-D decomposition along the :math:`y` axis can only be
used when the PP decomposition has only 1 domain along :math:`x`. 2-D
PME decomposition has to have the number of domains along :math:`x`
equal to the number of the PP decomposition. mdrun automatically chooses
1-D or 2-D PME decomposition (when possible with the total given number
of ranks), based on the minimum amount of communication for the
coordinate redistribution in PME plus the communication for the grid
overlap and transposes. To avoid superfluous communication of
coordinates and forces between the PP and PME ranks, the number of DD
cells in the :math:`x` direction should ideally be the same or a
multiple of the number of PME ranks. By default, mdrun takes care of
this issue.

Domain decomposition flow chart
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In Fig. [fig:dd\_flow] a flow chart is shown for domain decomposition
with all possible communication for different algorithms. For simpler
simulations, the same flow chart applies, without the algorithms and
communication for the algorithms that are not used.

.. figure:: plots/flowchart
   :alt:  Flow chart showing the algorithms and communication (arrows)
   for a standard MD simulation with virtual sites, constraints and
   separate PME-mesh ranks. [fig:dd\_flow]
   :width: 12.00000cm

    Flow chart showing the algorithms and communication (arrows) for a
   standard MD simulation with virtual sites, constraints and separate
   PME-mesh ranks. [fig:dd\_flow] 

Interaction function and force fields
=====================================

To accommodate the potential functions used in some popular force fields
(see [sec:ff]), GROMACS offers a choice of functions, both for
non-bonded interaction and for dihedral interactions. They are described
in the appropriate subsections.

The potential functions can be subdivided into three parts

#. *Non-bonded*: Lennard-Jones or Buckingham, and Coulomb or modified
   Coulomb. The non-bonded interactions are computed on the basis of a
   neighbor list (a list of non-bonded atoms within a certain radius),
   in which exclusions are already removed.

#. *Bonded*: covalent bond-stretching, angle-bending, improper
   dihedrals, and proper dihedrals. These are computed on the basis of
   fixed lists.

#. *Restraints*: position restraints, angle restraints, distance
   restraints, orientation restraints and dihedral restraints, all based
   on fixed lists.

#. *Applied Forces*: externally applied forces, see
   chapter [ch:special].

Non-bonded interactions
-----------------------

Non-bonded interactions in GROMACS are pair-additive:

.. math:: V({\mbox{\boldmath ${r}$}}_1,\ldots {\mbox{\boldmath ${r}$}}_N) = \sum_{i<j}V_{ij}({{\mbox{\boldmath ${r}$}}_{ij}});

.. math:: {\mbox{\boldmath ${F}$}}_i = -\sum_j \frac{dV_{ij}(r_{ij})}{dr_{ij}} \frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}

 Since the potential only depends on the scalar distance, interactions
will be centro-symmetric, i.e. the vectorial partial force on particle
:math:`i` from the pairwise interaction :math:`V_{ij}(r_{ij})` has the
opposite direction of the partial force on particle :math:`j`. For
efficiency reasons, interactions are calculated by loops over
interactions and updating both partial forces rather than summing one
complete nonbonded force at a time. The non-bonded interactions contain
a repulsion term, a dispersion term, and a Coulomb term. The repulsion
and dispersion term are combined in either the Lennard-Jones (or 6-12
interaction), or the Buckingham (or exp-6 potential). In addition,
(partially) charged atoms act through the Coulomb term.

The Lennard-Jones interaction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Lennard-Jones potential :math:`V_{LJ}` between two atoms equals:

.. math::

   V_{LJ}({r_{ij}}) =  \frac{C_{ij}^{(12)}}{{r_{ij}}^{12}} -
                           \frac{C_{ij}^{(6)}}{{r_{ij}}^6}

 See also Fig. [fig:lj] The parameters :math:`C^{(12)}_{ij}` and
:math:`C^{(6)}_{ij}` depend on pairs of *atom types*; consequently they
are taken from a matrix of LJ-parameters. In the Verlet cut-off scheme,
the potential is shifted by a constant such that it is zero at the
cut-off distance.

.. figure:: plots/f-lj
   :alt: The Lennard-Jones interaction.
   :width: 8.00000cm

   The Lennard-Jones interaction.

The force derived from this potential is:

.. math::

   {\mbox{\boldmath ${F}$}}_i({{\mbox{\boldmath ${r}$}}_{ij}}) = \left( 12~\frac{C_{ij}^{(12)}}{{r_{ij}}^{13}} -
                                    6~\frac{C_{ij}^{(6)}}{{r_{ij}}^7} \right) {\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{{r_{ij}}}}

The LJ potential may also be written in the following form:

.. math::

   V_{LJ}({{\mbox{\boldmath ${r}$}}_{ij}}) = 4\epsilon_{ij}\left(\left(\frac{\sigma_{ij}} {{r_{ij}}}\right)^{12}
                   - \left(\frac{\sigma_{ij}}{{r_{ij}}}\right)^{6} \right)
   \label{eqn:sigeps}

In constructing the parameter matrix for the non-bonded LJ-parameters,
two types of combination rules can be used within GROMACS, only
geometric averages (type 1 in the input section of the force-field
file):

.. math::

   \begin{array}{rcl}
   C_{ij}^{(6)}    &=& \left({C_{ii}^{(6)} \, C_{jj}^{(6)}}\right)^{1/2}    \\
   C_{ij}^{(12)}   &=& \left({C_{ii}^{(12)} \, C_{jj}^{(12)}}\right)^{1/2}
   \label{eqn:comb}
   \end{array}

 or, alternatively the Lorentz-Berthelot rules can be used. An
arithmetic average is used to calculate :math:`\sigma_{ij}`, while a
geometric average is used to calculate :math:`\epsilon_{ij}` (type 2):

.. math::

   \begin{array}{rcl}
    \sigma_{ij}   &=& \frac{1}{ 2}(\sigma_{ii} + \sigma_{jj})        \\
    \epsilon_{ij} &=& \left({\epsilon_{ii} \, \epsilon_{jj}}\right)^{1/2}
    \label{eqn:lorentzberthelot}
   \end{array}

 finally an geometric average for both parameters can be used (type 3):

.. math::

   \begin{array}{rcl}
    \sigma_{ij}   &=& \left({\sigma_{ii} \, \sigma_{jj}}\right)^{1/2}        \\
    \epsilon_{ij} &=& \left({\epsilon_{ii} \, \epsilon_{jj}}\right)^{1/2}
   \end{array}

 This last rule is used by the OPLS force field.

Buckingham potential
~~~~~~~~~~~~~~~~~~~~

The Buckingham potential has a more flexible and realistic repulsion
term than the Lennard-Jones interaction, but is also more expensive to
compute. The potential form is:

.. math::

   V_{bh}({r_{ij}}) = A_{ij} \exp(-B_{ij} {r_{ij}}) -
                           \frac{C_{ij}}{{r_{ij}}^6}

.. figure:: plots/f-bham
   :alt: The Buckingham interaction.
   :width: 8.00000cm

   The Buckingham interaction.

See also Fig. [fig:bham]. The force derived from this is:

.. math::

   {\mbox{\boldmath ${F}$}}_i({r_{ij}}) = \left[ A_{ij}B_{ij}\exp(-B_{ij} {r_{ij}}) -
                                    6\frac{C_{ij}}{{r_{ij}}^7} \right] {\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{{r_{ij}}}}

Coulomb interaction
~~~~~~~~~~~~~~~~~~~

The Coulomb interaction between two charge particles is given by:

.. math::

   V_c({r_{ij}}) = f \frac{q_i q_j}{{\varepsilon_r}{r_{ij}}}
   \label{eqn:vcoul}

 See also Fig. [fig:coul], where
:math:`f = \frac{1}{4\pi \varepsilon_0} =
{138.935\,458}` (see chapter [ch:defunits])

.. figure:: plots/vcrf
   :alt: The Coulomb interaction (for particles with equal signed
   charge) with and without reaction field. In the latter case
   :math:`{\varepsilon_r}` was 1, :math:`{\varepsilon_{rf}}` was 78, and
   :math:`r_c` was 0.9 nm. The dot-dashed line is the same as the dashed
   line, except for a constant.
   :width: 8.00000cm

   The Coulomb interaction (for particles with equal signed charge) with
   and without reaction field. In the latter case
   :math:`{\varepsilon_r}` was 1, :math:`{\varepsilon_{rf}}` was 78, and
   :math:`r_c` was 0.9 nm. The dot-dashed line is the same as the dashed
   line, except for a constant.

The force derived from this potential is:

.. math:: {\mbox{\boldmath ${F}$}}_i({{\mbox{\boldmath ${r}$}}_{ij}}) = f \frac{q_i q_j}{{\varepsilon_r}{r_{ij}}^2}{\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{{r_{ij}}}}

A plain Coulomb interaction should only be used without cut-off or when
all pairs fall within the cut-off, since there is an abrupt, large
change in the force at the cut-off. In case you do want to use a
cut-off, the potential can be shifted by a constant to make the
potential the integral of the force. With the group cut-off scheme, this
shift is only applied to non-excluded pairs. With the Verlet cut-off
scheme, the shift is also applied to excluded pairs and self
interactions, which makes the potential equivalent to a reaction field
with :math:`{\varepsilon_{rf}}=1` (see below).

In GROMACS the relative dielectric constant :math:`{\varepsilon_r}` may
be set in the in the input for grompp.

Coulomb interaction with reaction field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Coulomb interaction can be modified for homogeneous systems by
assuming a constant dielectric environment beyond the cut-off
:math:`r_c` with a dielectric constant of :math:`{\varepsilon_{rf}}`.
The interaction then reads:

.. math::

   V_{crf} ~=~
     f \frac{q_i q_j}{{\varepsilon_r}{r_{ij}}}\left[1+\frac{{\varepsilon_{rf}}-{\varepsilon_r}}{2{\varepsilon_{rf}}+{\varepsilon_r}}
     \,\frac{{r_{ij}}^3}{r_c^3}\right]
     - f\frac{q_i q_j}{{\varepsilon_r}r_c}\,\frac{3{\varepsilon_{rf}}}{2{\varepsilon_{rf}}+{\varepsilon_r}}
   \label{eqn:vcrf}

 in which the constant expression on the right makes the potential zero
at the cut-off :math:`r_c`. For charged cut-off spheres this corresponds
to neutralization with a homogeneous background charge. We can rewrite
eqn. [eqn:vcrf] for simplicity as

.. math:: V_{crf} ~=~     f \frac{q_i q_j}{{\varepsilon_r}}\left[\frac{1}{{r_{ij}}} + k_{rf}~ {r_{ij}}^2 -c_{rf}\right]

 with

.. math::

   \begin{aligned}
   k_{rf}  &=&     \frac{1}{r_c^3}\,\frac{{\varepsilon_{rf}}-{\varepsilon_r}}{(2{\varepsilon_{rf}}+{\varepsilon_r})}   \label{eqn:krf}\\
   c_{rf}  &=&     \frac{1}{r_c}+k_{rf}\,r_c^2 ~=~ \frac{1}{r_c}\,\frac{3{\varepsilon_{rf}}}{(2{\varepsilon_{rf}}+{\varepsilon_r})}
   \label{eqn:crf}\end{aligned}

 For large :math:`{\varepsilon_{rf}}` the :math:`k_{rf}` goes to
:math:`r_c^{-3}/2`, while for :math:`{\varepsilon_{rf}}` =
:math:`{\varepsilon_r}` the correction vanishes. In Fig. [fig:coul] the
modified interaction is plotted, and it is clear that the derivative
with respect to :math:`{r_{ij}}` (= -force) goes to zero at the cut-off
distance. The force derived from this potential reads:

.. math:: {\mbox{\boldmath ${F}$}}_i({{\mbox{\boldmath ${r}$}}_{ij}}) = f \frac{q_i q_j}{{\varepsilon_r}}\left[\frac{1}{{r_{ij}}^2} - 2 k_{rf}{r_{ij}}\right]{\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{{r_{ij}}}}\label{eqn:fcrf}

 The reaction-field correction should also be applied to all excluded
atoms pairs, including self pairs, in which case the normal Coulomb term
in eqns. [eqn:vcrf] and [eqn:fcrf] is absent.

Tironi *et al.* have introduced a generalized reaction field in which
the dielectric continuum beyond the cut-off :math:`r_c` also has an
ionic strength :math:`I` `71 <#ref-Tironi95>`__. In this case we can
rewrite the constants :math:`k_{rf}` and :math:`c_{rf}` using the
inverse Debye screening length :math:`\kappa`:

.. math::

   \begin{aligned}
   \kappa^2  &=&     
      \frac{2 I \,F^2}{\varepsilon_0 {\varepsilon_{rf}}RT}
      = \frac{F^2}{\varepsilon_0 {\varepsilon_{rf}}RT}\sum_{i=1}^{K} c_i z_i^2     \\
   k_{rf}  &=&     \frac{1}{r_c^3}\,
       \frac{({\varepsilon_{rf}}-{\varepsilon_r})(1 + \kappa r_c) + {\frac{1}{2}}{\varepsilon_{rf}}(\kappa r_c)^2}
            {(2{\varepsilon_{rf}}+ {\varepsilon_r})(1 + \kappa r_c) + {\varepsilon_{rf}}(\kappa r_c)^2}
       \label{eqn:kgrf}\\
   c_{rf}  &=&     \frac{1}{r_c}\,
       \frac{3{\varepsilon_{rf}}(1 + \kappa r_c + {\frac{1}{2}}(\kappa r_c)^2)}
            {(2{\varepsilon_{rf}}+{\varepsilon_r})(1 + \kappa r_c) + {\varepsilon_{rf}}(\kappa r_c)^2}
       \label{eqn:cgrf}\end{aligned}

 where :math:`F` is Faraday’s constant, :math:`R` is the ideal gas
constant, :math:`T` the absolute temperature, :math:`c_i` the molar
concentration for species :math:`i` and :math:`z_i` the charge number of
species :math:`i` where we have :math:`K` different species. In the
limit of zero ionic strength (:math:`\kappa=0`) eqns. [eqn:kgrf] and
[eqn:cgrf] reduce to the simple forms of eqns. [eqn:krf] and [eqn:crf]
respectively.

Modified non-bonded interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In GROMACS, the non-bonded potentials can be modified by a shift
function, also called a force-switch function, since it switches the
force to zero at the cut-off. The purpose of this is to replace the
truncated forces by forces that are continuous and have continuous
derivatives at the cut-off radius. With such forces the time integration
produces smaller errors. But note that for Lennard-Jones interactions
these errors are usually smaller than other errors, such as integration
errors at the repulsive part of the potential. For Coulomb interactions
we advise against using a shifted potential and for use of a reaction
field or a proper long-range method such as PME.

There is *no* fundamental difference between a switch function (which
multiplies the potential with a function) and a shift function (which
adds a function to the force or potential) \ `72 <#ref-Spoel2006a>`__.
The switch function is a special case of the shift function, which we
apply to the *force function* :math:`F(r)`, related to the electrostatic
or van der Waals force acting on particle :math:`i` by particle
:math:`j` as:

.. math:: {\mbox{\boldmath ${F}$}}_i = c \, F(r_{ij}) \frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}

 For pure Coulomb or Lennard-Jones interactions
:math:`F(r) = F_\alpha(r) = \alpha \, r^{-(\alpha+1)}`. The switched
force :math:`F_s(r)` can generally be written as:

.. math::

   \begin{array}{rcl}
   \vspace{2mm}
   F_s(r)~=&~F_\alpha(r)   & r < r_1               \\
   \vspace{2mm}
   F_s(r)~=&~F_\alpha(r)+S(r)      & r_1 \le r < r_c       \\
   F_s(r)~=&~0             & r_c \le r     
   \end{array}

 When :math:`r_1=0` this is a traditional shift function, otherwise it
acts as a switch function. The corresponding shifted potential function
then reads:

.. math:: V_s(r) =  \int^{\infty}_r~F_s(x)\, dx

The GROMACS **force switch** function :math:`S_F(r)` should be smooth at
the boundaries, therefore the following boundary conditions are imposed
on the switch function:

.. math::

   \begin{array}{rcl}
   S_F(r_1)          &=&0            \\
   S_F'(r_1)         &=&0            \\
   S_F(r_c)          &=&-F_\alpha(r_c)       \\
   S_F'(r_c)         &=&-F_\alpha'(r_c)
   \end{array}

 A 3\ :math:`^{rd}` degree polynomial of the form

.. math:: S_F(r) = A(r-r_1)^2 + B(r-r_1)^3

 fulfills these requirements. The constants A and B are given by the
boundary condition at :math:`r_c`:

.. math::

   \begin{array}{rcl}
   \vspace{2mm}
   A &~=~& -\alpha \, \displaystyle
           \frac{(\alpha+4)r_c~-~(\alpha+1)r_1} {r_c^{\alpha+2}~(r_c-r_1)^2} \\
   B &~=~& \alpha \, \displaystyle
           \frac{(\alpha+3)r_c~-~(\alpha+1)r_1}{r_c^{\alpha+2}~(r_c-r_1)^3}
   \end{array}

 Thus the total force function is:

.. math:: F_s(r) = \frac{\alpha}{r^{\alpha+1}} + A(r-r_1)^2 + B(r-r_1)^3

 and the potential function reads:

.. math:: V_s(r) = \frac{1}{r^\alpha} - \frac{A}{3} (r-r_1)^3 - \frac{B}{4} (r-r_1)^4 - C

 where

.. math:: C =  \frac{1}{r_c^\alpha} - \frac{A}{3} (r_c-r_1)^3 - \frac{B}{4} (r_c-r_1)^4

The GROMACS **potential-switch** function :math:`S_V(r)` scales the
potential between :math:`r_1` and :math:`r_c`, and has similar boundary
conditions, intended to produce smoothly-varying potential and forces:

.. math::

   \begin{array}{rcl}
   S_V(r_1)          &=&1 \\
   S_V'(r_1)         &=&0 \\
   S_V''(r_1)        &=&0 \\
   S_V(r_c)          &=&0 \\
   S_V'(r_c)         &=&0 \\
   S_V''(r_c)        &=&0
   \end{array}

The fifth-degree polynomial that has these properties is

.. math:: S_V(r; r_1, r_c) = \frac{1 - 10(r-r_1)^3(r_c-r_1)^2 + 15(r-r_1)^4(r_c-r_1) - 6(r-r_1)}{(r_c-r_1)^5}

This implementation is found in several other simulation
packages,\ `73 <#ref-Ohmine1988>`__\ `75 <#ref-Guenot1993>`__ but
differs from that in CHARMM.\ `76 <#ref-Steinbach1994>`__ Switching the
potential leads to artificially large forces in the switching region,
therefore it is not recommended to switch Coulomb interactions using
this function,\ `72 <#ref-Spoel2006a>`__ but switching Lennard-Jones
interactions using this function produces acceptable results.

Modified short-range interactions with Ewald summation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When Ewald summation or particle-mesh Ewald is used to calculate the
long-range interactions, the short-range Coulomb potential must also be
modified. Here the potential is switched to (nearly) zero at the
cut-off, instead of the force. In this case the short range potential is
given by:

.. math:: V(r) = f \frac{\mbox{erfc}(\beta r_{ij})}{r_{ij}} q_i q_j,

 where :math:`\beta` is a parameter that determines the relative weight
between the direct space sum and the reciprocal space sum and
erfc\ :math:`(x)` is the complementary error function. For further
details on long-range electrostatics, see sec. [sec:lr\_elstat].

Bonded interactions
-------------------

Bonded interactions are based on a fixed list of atoms. They are not
exclusively pair interactions, but include 3- and 4-body interactions as
well. There are *bond stretching* (2-body), *bond angle* (3-body), and
*dihedral angle* (4-body) interactions. A special type of dihedral
interaction (called *improper dihedral*) is used to force atoms to
remain in a plane or to prevent transition to a configuration of
opposite chirality (a mirror image).

Bond stretching
~~~~~~~~~~~~~~~

Harmonic potential
^^^^^^^^^^^^^^^^^^

The bond stretching between two covalently bonded atoms :math:`i` and
:math:`j` is represented by a harmonic potential:

.. figure:: plots/f-bond
   :alt: Principle of bond stretching (left), and the bond stretching
   potential (right).
   :width: 7.00000cm

   Principle of bond stretching (left), and the bond stretching
   potential (right).

.. math:: V_b~({r_{ij}}) = {\frac{1}{2}}k^b_{ij}({r_{ij}}-b_{ij})^2

See also Fig. [fig:bstretch1], with the force given by:

.. math:: {\mbox{\boldmath ${F}$}}_i({{\mbox{\boldmath ${r}$}}_{ij}}) = k^b_{ij}({r_{ij}}-b_{ij}) {\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{{r_{ij}}}}

Fourth power potential
^^^^^^^^^^^^^^^^^^^^^^

In the GROMOS-96 force field \ `77 <#ref-gromos96>`__, the covalent bond
potential is, for reasons of computational efficiency, written as:

.. math:: V_b~({r_{ij}}) = \frac{1}{4}k^b_{ij}\left({r_{ij}}^2-b_{ij}^2\right)^2

 The corresponding force is:

.. math:: {\mbox{\boldmath ${F}$}}_i({{\mbox{\boldmath ${r}$}}_{ij}}) = k^b_{ij}({r_{ij}}^2-b_{ij}^2)~{{\mbox{\boldmath ${r}$}}_{ij}}

 The force constants for this form of the potential are related to the
usual harmonic force constant :math:`k^{b,\mathrm{harm}}`
(sec. [sec:bondpot]) as

.. math:: 2 k^b b_{ij}^2 = k^{b,\mathrm{harm}}

 The force constants are mostly derived from the harmonic ones used in
GROMOS-87 `78 <#ref-biomos>`__. Although this form is computationally
more efficient (because no square root has to be evaluated), it is
conceptually more complex. One particular disadvantage is that since the
form is not harmonic, the average energy of a single bond is not equal
to :math:`{\frac{1}{2}}kT` as it is for the normal harmonic potential.

Morse potential bond stretching
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For some systems that require an anharmonic bond stretching potential,
the Morse potential \ `79 <#ref-Morse29>`__ between two atoms *i* and
*j* is available in GROMACS. This potential differs from the harmonic
potential in that it has an asymmetric potential well and a zero force
at infinite distance. The functional form is:

.. math:: \displaystyle V_{morse} (r_{ij}) = D_{ij} [1 - \exp(-\beta_{ij}(r_{ij}-b_{ij}))]^2,

 See also Fig. [fig:morse], and the corresponding force is:

.. math::

   \begin{array}{rcl}
   \displaystyle {\bf F}_{morse} ({\bf r}_{ij})&=&2 D_{ij} \beta_{ij} \exp(-\beta_{ij}(r_{ij}-b_{ij})) * \\
   \displaystyle \: & \: &[1 - \exp(-\beta_{ij}(r_{ij}-b_{ij}))] \frac{\displaystyle {\bf r}_{ij}}{\displaystyle r_{ij}},
   \end{array}

 where :math:` \displaystyle D_{ij} ` is the depth of the well in
kJ/mol, :math:` \displaystyle \beta_{ij} ` defines the steepness of the
well (in nm\ :math:`^{-1} `), and :math:` \displaystyle b_{ij} ` is the
equilibrium distance in nm. The steepness parameter
:math:` \displaystyle \beta_{ij}
` can be expressed in terms of the reduced mass of the atoms *i* and
*j*, the fundamental vibration frequency :math:` \displaystyle
\omega_{ij} ` and the well depth :math:` \displaystyle D_{ij} `:

.. math:: \displaystyle \beta_{ij}= \omega_{ij} \sqrt{\frac{\mu_{ij}}{2 D_{ij}}}

 and because :math:` \displaystyle \omega = \sqrt{k/\mu} `, one can
rewrite :math:` \displaystyle \beta_{ij} ` in terms of the harmonic
force constant :math:` \displaystyle k_{ij} `:

.. math::

   \displaystyle \beta_{ij}= \sqrt{\frac{k_{ij}}{2 D_{ij}}}
   \label{eqn:betaij}

 For small deviations :math:` \displaystyle (r_{ij}-b_{ij}) `, one can
approximate the :math:` \displaystyle \exp `-term to first-order using a
Taylor expansion:

.. math::

   \displaystyle \exp(-x) \approx 1-x
   \label{eqn:expminx}

 and substituting eqn. [eqn:betaij] and eqn. [eqn:expminx] in the
functional form:

.. math::

   \begin{array}{rcl}
   \displaystyle V_{morse} (r_{ij})&=&D_{ij} [1 - \exp(-\beta_{ij}(r_{ij}-b_{ij}))]^2\\
   \displaystyle \:&=&D_{ij} [1 - (1 -\sqrt{\frac{k_{ij}}{2 D_{ij}}}(r_{ij}-b_{ij}))]^2\\
   \displaystyle \:&=&\frac{1}{2} k_{ij} (r_{ij}-b_{ij}))^2
   \end{array}

 we recover the harmonic bond stretching potential.

.. figure:: plots/f-morse
   :alt: The Morse potential well, with bond length 0.15 nm.
   :width: 7.00000cm

   The Morse potential well, with bond length 0.15 nm.

Cubic bond stretching potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another anharmonic bond stretching potential that is slightly simpler
than the Morse potential adds a cubic term in the distance to the simple
harmonic form:

.. math:: V_b~({r_{ij}}) = k^b_{ij}({r_{ij}}-b_{ij})^2 + k^b_{ij}k^{cub}_{ij}({r_{ij}}-b_{ij})^3

 A flexible water model (based on the SPC water
model \ `80 <#ref-Berendsen81>`__) including a cubic bond stretching
potential for the O-H bond was developed by
Ferguson \ `81 <#ref-Ferguson95>`__. This model was found to yield a
reasonable infrared spectrum. The Ferguson water model is available in
the GROMACS library (flexwat-ferguson.itp). It should be noted that the
potential is asymmetric: overstretching leads to infinitely low
energies. The integration timestep is therefore limited to 1 fs.

The force corresponding to this potential is:

.. math:: {\mbox{\boldmath ${F}$}}_i({{\mbox{\boldmath ${r}$}}_{ij}}) = 2k^b_{ij}({r_{ij}}-b_{ij})~{\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{{r_{ij}}}}+ 3k^b_{ij}k^{cub}_{ij}({r_{ij}}-b_{ij})^2~{\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{{r_{ij}}}}

FENE bond stretching potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In coarse-grained polymer simulations the begineqnarrayds are often
connected by a FENE (finitely extensible nonlinear elastic)
potential \ `82 <#ref-Warner72>`__:

.. math::

   V_{\mbox{\small FENE}}({r_{ij}}) =
     -{\frac{1}{2}}k^b_{ij} b^2_{ij} \log\left(1 - \frac{{r_{ij}}^2}{b^2_{ij}}\right)

 The potential looks complicated, but the expression for the force is
simpler:

.. math::

   F_{\mbox{\small FENE}}({{\mbox{\boldmath ${r}$}}_{ij}}) =
     -k^b_{ij} \left(1 - \frac{{r_{ij}}^2}{b^2_{ij}}\right)^{-1} {{\mbox{\boldmath ${r}$}}_{ij}}

 At short distances the potential asymptotically goes to a harmonic
potential with force constant :math:`k^b`, while it diverges at distance
:math:`b`.

Harmonic angle potential
~~~~~~~~~~~~~~~~~~~~~~~~

The bond-angle vibration between a triplet of atoms :math:`i` -
:math:`j` - :math:`k` is also represented by a harmonic potential on the
angle :math:`{\theta_{ijk}}`

.. figure:: plots/f-angle
   :alt: Principle of angle vibration (left) and the bond angle
   potential (right).
   :width: 7.00000cm

   Principle of angle vibration (left) and the bond angle potential
   (right).

.. math:: V_a({\theta_{ijk}}) = {\frac{1}{2}}k^{\theta}_{ijk}({\theta_{ijk}}-{\theta_{ijk}}^0)^2

As the bond-angle vibration is represented by a harmonic potential, the
form is the same as the bond stretching (Fig. [fig:bstretch1]).

The force equations are given by the chain rule:

.. math::

   \begin{array}{l}
   {{\mbox{\boldmath ${F}$}}_i}~=~ -\displaystyle\frac{d V_a({\theta_{ijk}})}{d {{\mbox{\boldmath ${r}$}}_i}}   \\
   {{\mbox{\boldmath ${F}$}}_k}~=~ -\displaystyle\frac{d V_a({\theta_{ijk}})}{d {{\mbox{\boldmath ${r}$}}_k}}   \\
   {{\mbox{\boldmath ${F}$}}_j}~=~ -{{\mbox{\boldmath ${F}$}}_i}-{{\mbox{\boldmath ${F}$}}_k}\end{array}
   ~ \mbox{ ~ where ~ } ~
    {\theta_{ijk}}= \arccos \frac{({{\mbox{\boldmath ${r}$}}_{ij}}\cdot {\mbox{\boldmath ${r}$}}_{kj})}{r_{ij}r_{kj}}

 The numbering :math:`i,j,k` is in sequence of covalently bonded atoms.
Atom :math:`j` is in the middle; atoms :math:`i` and :math:`k` are at
the ends (see Fig. [fig:angle]). **Note** that in the input in topology
files, angles are given in degrees and force constants in
kJ/mol/rad\ :math:`^2`.

Cosine based angle potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the GROMOS-96 force field a simplified function is used to represent
angle vibrations:

.. math::

   V_a({\theta_{ijk}}) = {\frac{1}{2}}k^{\theta}_{ijk}\left(\cos({\theta_{ijk}}) - \cos({\theta_{ijk}}^0)\right)^2
   \label{eq:G96angle}

 where

.. math:: \cos({\theta_{ijk}}) = \frac{{{\mbox{\boldmath ${r}$}}_{ij}}\cdot{\mbox{\boldmath ${r}$}}_{kj}}{{r_{ij}}r_{kj}}

 The corresponding force can be derived by partial differentiation with
respect to the atomic positions. The force constants in this function
are related to the force constants in the harmonic form
:math:`k^{\theta,\mathrm{harm}}` ([subsec:harmonicangle]) by:

.. math:: k^{\theta} \sin^2({\theta_{ijk}}^0) = k^{\theta,\mathrm{harm}}

 In the GROMOS-96 manual there is a much more complicated conversion
formula which is temperature dependent. The formulas are equivalent at 0
K and the differences at 300 K are on the order of 0.1 to 0.2%. **Note**
that in the input in topology files, angles are given in degrees and
force constants in kJ/mol.

Restricted bending potential
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The restricted bending (ReB) potential \ `83 <#ref-MonicaGoga2013>`__
prevents the bending angle :math:`\theta` from reaching the
:math:`180^{\circ}` value. In this way, the numerical instabilities due
to the calculation of the torsion angle and potential are eliminated
when performing coarse-grained molecular dynamics simulations.

To systematically hinder the bending angles from reaching the
:math:`180^{\circ}` value, the bending potential [eq:G96angle] is
divided by a :math:`\sin^2\theta` factor:

.. math::

   V_{\rm ReB}(\theta_i) = \frac{1}{2} k_{\theta} \frac{(\cos\theta_i - \cos\theta_0)^2}{\sin^2\theta_i}.
   \label{eq:ReB}

 Figure  Fig. [fig:ReB] shows the comparison between the ReB potential,
[eq:ReB], and the standard one [eq:G96angle].

.. figure:: plots/fig-02
   :alt: Bending angle potentials: cosine harmonic (solid black line),
   angle harmonic (dashed black line) and restricted bending (red) with
   the same bending constant :math:`k_{\theta}=85` kJ mol\ :math:`^{-1}`
   and equilibrium angle :math:`\theta_0=130^{\circ}`. The orange line
   represents the sum of a cosine harmonic (:math:`k =50` kJ
   mol\ :math:`^{-1}`) with a restricted bending (:math:`k =25` kJ
   mol\ :math:`^{-1}`) potential, both with
   :math:`\theta_0=130^{\circ}`.
   :width: 10.00000cm

   Bending angle potentials: cosine harmonic (solid black line), angle
   harmonic (dashed black line) and restricted bending (red) with the
   same bending constant :math:`k_{\theta}=85` kJ mol\ :math:`^{-1}` and
   equilibrium angle :math:`\theta_0=130^{\circ}`. The orange line
   represents the sum of a cosine harmonic (:math:`k =50` kJ
   mol\ :math:`^{-1}`) with a restricted bending (:math:`k =25` kJ
   mol\ :math:`^{-1}`) potential, both with
   :math:`\theta_0=130^{\circ}`.

The wall of the ReB potential is very repulsive in the region close to
:math:`180^{\circ}` and, as a result, the bending angles are kept within
a safe interval, far from instabilities. The power :math:`2` of
:math:`\sin\theta_i` in the denominator has been chosen to guarantee
this behavior and allows an elegant differentiation:

.. math::

   F_{\rm ReB}(\theta_i) = \frac{2k_{\theta}}{\sin^4\theta_i}(\cos\theta_i - \cos\theta_0) (1 - \cos\theta_i\cos\theta_0) \frac{\partial \cos\theta_i}{\partial \vec r_{k}}.
   \label{eq:diff_ReB}

 Due to its construction, the restricted bending potential cannot be
used for equilibrium :math:`\theta_0` values too close to
:math:`0^{\circ}` or :math:`180^{\circ}` (from experience, at least
:math:`10^{\circ}` difference is recommended). It is very important
that, in the starting configuration, all the bending angles have to be
in the safe interval to avoid initial instabilities. This bending
potential can be used in combination with any form of torsion potential.
It will always prevent three consecutive particles from becoming
collinear and, as a result, any torsion potential will remain free of
singularities. It can be also added to a standard bending potential to
affect the angle around :math:`180^{\circ}`, but to keep its original
form around the minimum (see the orange curve in Fig. [fig:ReB]).

Urey-Bradley potential
~~~~~~~~~~~~~~~~~~~~~~

The Urey-Bradley bond-angle vibration between a triplet of atoms
:math:`i` - :math:`j` - :math:`k` is represented by a harmonic potential
on the angle :math:`{\theta_{ijk}}` and a harmonic correction term on
the distance between the atoms :math:`i` and :math:`k`. Although this
can be easily written as a simple sum of two terms, it is convenient to
have it as a single entry in the topology file and in the output as a
separate energy term. It is used mainly in the CHARMm force
field \ `84 <#ref-BBrooks83>`__. The energy is given by:

.. math:: V_a({\theta_{ijk}}) = {\frac{1}{2}}k^{\theta}_{ijk}({\theta_{ijk}}-{\theta_{ijk}}^0)^2 + {\frac{1}{2}}k^{UB}_{ijk}(r_{ik}-r_{ik}^0)^2

The force equations can be deduced from sections [subsec:harmonicbond]
and [subsec:harmonicangle].

Bond-Bond cross term
~~~~~~~~~~~~~~~~~~~~

The bond-bond cross term for three particles :math:`i, j, k` forming
bonds :math:`i-j` and :math:`k-j` is given
by \ `85 <#ref-Lawrence2003b>`__:

.. math::

   V_{rr'} ~=~ k_{rr'} \left(\left|{\mbox{\boldmath ${r}$}}_{i}-{\mbox{\boldmath ${r}$}}_j\right|-r_{1e}\right) \left(\left|{\mbox{\boldmath ${r}$}}_{k}-{\mbox{\boldmath ${r}$}}_j\right|-r_{2e}\right)
   \label{crossbb}

 where :math:`k_{rr'}` is the force constant, and :math:`r_{1e}` and
:math:`r_{2e}` are the equilibrium bond lengths of the :math:`i-j` and
:math:`k-j` bonds respectively. The force associated with this potential
on particle :math:`i` is:

.. math:: {\mbox{\boldmath ${F}$}}_{i} = -k_{rr'}\left(\left|{\mbox{\boldmath ${r}$}}_{k}-{\mbox{\boldmath ${r}$}}_j\right|-r_{2e}\right)\frac{{\mbox{\boldmath ${r}$}}_i-{\mbox{\boldmath ${r}$}}_j}{\left|{\mbox{\boldmath ${r}$}}_{i}-{\mbox{\boldmath ${r}$}}_j\right|}

 The force on atom :math:`k` can be obtained by swapping :math:`i` and
:math:`k` in the above equation. Finally, the force on atom :math:`j`
follows from the fact that the sum of internal forces should be zero:
:math:`{\mbox{\boldmath ${F}$}}_j = -{\mbox{\boldmath ${F}$}}_i-{\mbox{\boldmath ${F}$}}_k`.

Bond-Angle cross term
~~~~~~~~~~~~~~~~~~~~~

The bond-angle cross term for three particles :math:`i, j, k` forming
bonds :math:`i-j` and :math:`k-j` is given
by \ `85 <#ref-Lawrence2003b>`__:

.. math:: V_{r\theta} ~=~ k_{r\theta} \left(\left|{\mbox{\boldmath ${r}$}}_{i}-{\mbox{\boldmath ${r}$}}_k\right|-r_{3e} \right) \left(\left|{\mbox{\boldmath ${r}$}}_{i}-{\mbox{\boldmath ${r}$}}_j\right|-r_{1e} + \left|{\mbox{\boldmath ${r}$}}_{k}-{\mbox{\boldmath ${r}$}}_j\right|-r_{2e}\right)

 where :math:`k_{r\theta}` is the force constant, :math:`r_{3e}` is the
:math:`i-k` distance, and the other constants are the same as in
Equation [crossbb]. The force associated with the potential on atom
:math:`i` is:

.. math::

   {\mbox{\boldmath ${F}$}}_{i} ~=~ -k_{r\theta}\left[\left(\left|{\mbox{\boldmath ${r}$}}_{i}-{\mbox{\boldmath ${r}$}}_{k}\right|-r_{3e}\right)\frac{{\mbox{\boldmath ${r}$}}_i-{\mbox{\boldmath ${r}$}}_j}{\left|{\mbox{\boldmath ${r}$}}_{i}-{\mbox{\boldmath ${r}$}}_j\right|} \\
   + \left(\left|{\mbox{\boldmath ${r}$}}_{i}-{\mbox{\boldmath ${r}$}}_j\right|-r_{1e} + \left|{\mbox{\boldmath ${r}$}}_{k}-{\mbox{\boldmath ${r}$}}_j\right|-r_{2e}\right)\frac{{\mbox{\boldmath ${r}$}}_i-{\mbox{\boldmath ${r}$}}_k}{\left|{\mbox{\boldmath ${r}$}}_{i}-{\mbox{\boldmath ${r}$}}_k\right|}\right]

Quartic angle potential
~~~~~~~~~~~~~~~~~~~~~~~

For special purposes there is an angle potential that uses a fourth
order polynomial:

.. math:: V_q({\theta_{ijk}}) ~=~ \sum_{n=0}^5 C_n ({\theta_{ijk}}-{\theta_{ijk}}^0)^n

Improper dihedrals
~~~~~~~~~~~~~~~~~~

Improper dihedrals are meant to keep planar groups (*e.g.* aromatic
rings) planar, or to prevent molecules from flipping over to their
mirror images, see Fig. [fig:imp].

|Principle of improper dihedral angles. Out of plane bending for rings
(left), substituents of rings (middle), out of tetrahedral (right). The
improper dihedral angle :math:`\xi` is defined as the angle between
planes (i,j,k) and (j,k,l) in all cases.|  |Principle of improper
dihedral angles. Out of plane bending for rings (left), substituents of
rings (middle), out of tetrahedral (right). The improper dihedral angle
:math:`\xi` is defined as the angle between planes (i,j,k) and (j,k,l)
in all cases.| |Principle of improper dihedral angles. Out of plane
bending for rings (left), substituents of rings (middle), out of
tetrahedral (right). The improper dihedral angle :math:`\xi` is defined
as the angle between planes (i,j,k) and (j,k,l) in all cases.|

Improper dihedrals: harmonic type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The simplest improper dihedral potential is a harmonic potential; it is
plotted in Fig. [fig:imps].

.. math:: V_{id}(\xi_{ijkl}) = {\frac{1}{2}}k_{\xi}(\xi_{ijkl}-\xi_0)^2

 Since the potential is harmonic it is discontinuous, but since the
discontinuity is chosen at 180\ :math:`^\circ` distance from
:math:`\xi_0` this will never cause problems. **Note** that in the input
in topology files, angles are given in degrees and force constants in
kJ/mol/rad\ :math:`^2`.

.. figure:: plots/f-imps.pdf
   :alt: Improper dihedral potential.
   :width: 10.00000cm

   Improper dihedral potential.

Improper dihedrals: periodic type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This potential is identical to the periodic proper dihedral (see below).
There is a separate dihedral type for this (type 4) only to be able to
distinguish improper from proper dihedrals in the parameter section and
the output.

Proper dihedrals
~~~~~~~~~~~~~~~~

For the normal dihedral interaction there is a choice of either the
GROMOS periodic function or a function based on expansion in powers of
:math:`\cos \phi` (the so-called Ryckaert-Bellemans potential). This
choice has consequences for the inclusion of special interactions
between the first and the fourth atom of the dihedral quadruple. With
the periodic GROMOS potential a special 1-4 LJ-interaction must be
included; with the Ryckaert-Bellemans potential *for alkanes* the 1-4
interactions must be excluded from the non-bonded list. **Note:**
Ryckaert-Bellemans potentials are also used in *e.g.* the OPLS force
field in combination with 1-4 interactions. You should therefore not
modify topologies generated by pdb2gmx in this case.

Proper dihedrals: periodic type
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Proper dihedral angles are defined according to the IUPAC/IUB
convention, where :math:`\phi` is the angle between the :math:`ijk` and
the :math:`jkl` planes, with **zero** corresponding to the *cis*
configuration (:math:`i` and :math:`l` on the same side). There are two
dihedral function types in GROMACS topology files. There is the standard
type 1 which behaves like any other bonded interactions. For certain
force fields, type 9 is useful. Type 9 allows multiple potential
functions to be applied automatically to a single dihedral in the
section when multiple parameters are defined for the same atomtypes in
the section.

.. figure:: plots/f-dih
   :alt: Principle of proper dihedral angle (left, in *trans* form) and
   the dihedral angle potential (right).
   :width: 7.00000cm

   Principle of proper dihedral angle (left, in *trans* form) and the
   dihedral angle potential (right).

.. math:: V_d(\phi_{ijkl}) = k_{\phi}(1 + \cos(n \phi - \phi_s))

Proper dihedrals: Ryckaert-Bellemans function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| For alkanes, the following proper dihedral potential is often used
  (see Fig. [fig:rbdih]):

  .. math:: V_{rb}(\phi_{ijkl}) = \sum_{n=0}^5 C_n( \cos(\psi ))^n,

   where :math:`\psi = \phi - 180^\circ`.
| **Note:** A conversion from one convention to another can be achieved
  by multiplying every coefficient :math:` \displaystyle C_n ` by
  :math:` \displaystyle (-1)^n `.

An example of constants for :math:`C` is given in Table [tab:crb].

.. figure:: plots/f-rbs
   :alt: Ryckaert-Bellemans dihedral potential.
   :width: 8.00000cm

   Ryckaert-Bellemans dihedral potential.

(**Note:** The use of this potential implies exclusion of LJ
interactions between the first and the last atom of the dihedral, and
:math:`\psi` is defined according to the “polymer convention”
(:math:`\psi_{trans}=0`).)

| The RB dihedral function can also be used to include Fourier dihedrals
  (see below):

  .. math::

     V_{rb} (\phi_{ijkl}) ~=~ \frac{1}{2} \left[F_1(1+\cos(\phi)) + F_2(
     1-\cos(2\phi)) + F_3(1+\cos(3\phi)) + F_4(1-\cos(4\phi))\right]

   Because of the equalities :math:` \cos(2\phi) = 2\cos^2(\phi) - 1 `,
  :math:` \cos(3\phi) = 4\cos^3(\phi) - 3\cos(\phi) ` and
  :math:` \cos(4\phi) = 8\cos^4(\phi) - 8\cos^2(\phi) + 1 ` one can
  translate the OPLS parameters to Ryckaert-Bellemans parameters as
  follows:

  .. math::

     \displaystyle
     \begin{array}{rcl}
     \displaystyle C_0&=&F_2 + \frac{1}{2} (F_1 + F_3)\\
     \displaystyle C_1&=&\frac{1}{2} (- F_1 + 3 \, F_3)\\
     \displaystyle C_2&=& -F_2 + 4 \, F_4\\
     \displaystyle C_3&=&-2 \, F_3\\
     \displaystyle C_4&=&-4 \, F_4\\
     \displaystyle C_5&=&0
     \end{array}

   with OPLS parameters in protein convention and RB parameters in
  polymer convention (this yields a minus sign for the odd powers of
  cos\ :math:`(\phi)`).
| **Note:** Mind the conversion from **kcal mol\ :math:`^{-1}`** for
  literature OPLS and RB parameters to **kJ mol\ :math:`^{-1}`** in
  GROMACS.

Proper dihedrals: Fourier function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

| The OPLS potential function is given as the first three
   `86 <#ref-Jorgensen1996>`__ or four \ `87 <#ref-Robertson2015a>`__
  cosine terms of a Fourier series. In GROMACS the four term function is
  implemented:

  .. math::

     V_{F} (\phi_{ijkl}) ~=~ \frac{1}{2} \left[C_1(1+\cos(\phi)) + C_2(
     1-\cos(2\phi)) + C_3(1+\cos(3\phi)) + C_4(1-\cos(4\phi))\right],

   Internally, GROMACS uses the Ryckaert-Bellemans code to compute
  Fourier dihedrals (see above), because this is more efficient.
| **Note:** Mind the conversion from *k*\ cal mol\ :math:`^{-1}` for
  literature OPLS parameters to **kJ mol\ :math:`^{-1}`** in GROMACS.

Proper dihedrals: Restricted torsion potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In a manner very similar to the restricted bending potential (see
[subsec:ReB]), a restricted torsion/dihedral potential is introduced:

.. math::

   V_{\rm ReT}(\phi_i) = \frac{1}{2} k_{\phi} \frac{(\cos\phi_i - \cos\phi_0)^2}{\sin^2\phi_i}
   \label{eq:ReT}

 with the advantages of being a function of :math:`\cos\phi` (no
problems taking the derivative of :math:`\sin\phi`) and of keeping the
torsion angle at only one minimum value. In this case, the factor
:math:`\sin^2\phi` does not allow the dihedral angle to move from the
[:math:`-180^{\circ}`:0] to [0::math:`180^{\circ}`] interval, i.e. it
cannot have maxima both at :math:`-\phi_0` and :math:`+\phi_0` maxima,
but only one of them. For this reason, all the dihedral angles of the
starting configuration should have their values in the desired angles
interval and the the equilibrium :math:`\phi_0` value should not be too
close to the interval limits (as for the restricted bending potential,
described in [subsec:ReB], at least :math:`10^{\circ}` difference is
recommended).

Proper dihedrals: Combined bending-torsion potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When the four particles forming the dihedral angle become collinear
(this situation will never happen in atomistic simulations, but it can
occur in coarse-grained simulations) the calculation of the torsion
angle and potential leads to numerical instabilities. One way to avoid
this is to use the restricted bending potential (see [subsec:ReB]) that
prevents the dihedral from reaching the :math:`180^{\circ}` value.

Another way is to disregard any effects of the dihedral becoming
ill-defined, keeping the dihedral force and potential calculation
continuous in entire angle range by coupling the torsion potential (in a
cosine form) with the bending potentials of the adjacent bending angles
in a unique expression:

.. math::

   V_{\rm CBT}(\theta_{i-1}, \theta_i, \phi_i) = k_{\phi} \sin^3\theta_{i-1} \sin^3\theta_{i} \sum_{n=0}^4 { a_n \cos^n\phi_i}.
   \label{eq:CBT}

 This combined bending-torsion (CBT) potential has been proposed
by \ `88 <#ref-BulacuGiessen2005>`__ for polymer melt simulations and is
extensively described in \ `83 <#ref-MonicaGoga2013>`__.

This potential has two main advantages:

-  it does not only depend on the dihedral angle :math:`\phi_i` (between
   the :math:`i-2`, :math:`i-1`, :math:`i` and :math:`i+1`
   begineqnarrayds) but also on the bending angles :math:`\theta_{i-1}`
   and :math:`\theta_i` defined from three adjacent begineqnarrayds
   (:math:`i-2`, :math:`i-1` and :math:`i`, and :math:`i-1`, :math:`i`
   and :math:`i+1`, respectively). The two :math:`\sin^3\theta`
   pre-factors, tentatively suggested
   by \ `89 <#ref-ScottScheragator1966>`__ and theoretically discussed
   by \ `90 <#ref-PaulingBond>`__, cancel the torsion potential and
   force when either of the two bending angles approaches the value of
   :math:`180^\circ`.

-  its dependence on :math:`\phi_i` is expressed through a polynomial in
   :math:`\cos\phi_i` that avoids the singularities in
   :math:`\phi=0^\circ` or :math:`180^\circ` in calculating the
   torsional force.

These two properties make the CBT potential well-behaved for MD
simulations with weak constraints on the bending angles or even for
steered / non-equilibrium MD in which the bending and torsion angles
suffer major modifications. When using the CBT potential, the bending
potentials for the adjacent :math:`\theta_{i-1}` and :math:`\theta_i`
may have any form. It is also possible to leave out the two angle
bending terms (:math:`\theta_{i-1}` and :math:`\theta_{i}`) completely.
Fig. [fig:CBT] illustrates the difference between a torsion potential
with and without the :math:`\sin^{3}\theta` factors (blue and gray
curves, respectively).

.. figure:: plots/fig-04
   :alt: Blue: surface plot of the combined bending-torsion potential
   ([eq:CBT] with :math:`k = 10` kJ mol\ :math:`^{-1}`,
   :math:`a_0=2.41`, :math:`a_1=-2.95`, :math:`a_2=0.36`,
   :math:`a_3=1.33`) when, for simplicity, the bending angles behave the
   same (:math:`\theta_1=\theta_2=\theta`). Gray: the same torsion
   potential without the :math:`\sin^{3}\theta` terms
   (Ryckaert-Bellemans type). :math:`\phi` is the dihedral angle.
   :width: 10.00000cm

   Blue: surface plot of the combined bending-torsion potential
   ([eq:CBT] with :math:`k = 10` kJ mol\ :math:`^{-1}`,
   :math:`a_0=2.41`, :math:`a_1=-2.95`, :math:`a_2=0.36`,
   :math:`a_3=1.33`) when, for simplicity, the bending angles behave the
   same (:math:`\theta_1=\theta_2=\theta`). Gray: the same torsion
   potential without the :math:`\sin^{3}\theta` terms
   (Ryckaert-Bellemans type). :math:`\phi` is the dihedral angle.

Additionally, the derivative of :math:`V_{CBT}` with respect to the
Cartesian variables is straightforward:

.. math::

   \frac{\partial V_{\rm CBT}(\theta_{i-1},\theta_i,\phi_i)} {\partial \vec r_{l}} = \frac{\partial V_{\rm CBT}}{\partial \theta_{i-1}} \frac{\partial \theta_{i-1}}{\partial \vec r_{l}} +
                                                                                     \frac{\partial V_{\rm CBT}}{\partial \theta_{i  }} \frac{\partial \theta_{i  }}{\partial \vec r_{l}} +
                                                                                     \frac{\partial V_{\rm CBT}}{\partial \phi_{i    }} \frac{\partial \phi_{i    }}{\partial \vec r_{l}}
   \label{eq:force_cbt}

 The CBT is based on a cosine form without multiplicity, so it can only
be symmetrical around :math:`0^{\circ}`. To obtain an asymmetrical
dihedral angle distribution (e.g. only one maximum in
[:math:`-180^{\circ}`::math:`180^{\circ}`] interval), a standard torsion
potential such as harmonic angle or periodic cosine potentials should be
used instead of a CBT potential. However, these two forms have the
inconveniences of the force derivation (:math:`1/\sin\phi`) and of the
alignment of begineqnarrayds (:math:`\theta_i` or
:math:`\theta_{i-1} = 0^{\circ}, 180^{\circ}`). Coupling such
non-\ :math:`\cos\phi` potentials with :math:`\sin^3\theta` factors does
not improve simulation stability since there are cases in which
:math:`\theta` and :math:`\phi` are simultaneously :math:`180^{\circ}`.
The integration at this step would be possible (due to the cancelling of
the torsion potential) but the next step would be singular
(:math:`\theta` is not :math:`180^{\circ}` and :math:`\phi` is very
close to :math:`180^{\circ}`).

Tabulated bonded interaction functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

| For full flexibility, any functional shape can be used for bonds,
  angles and dihedrals through user-supplied tabulated functions. The
  functional shapes are:

  .. math::

     \begin{aligned}
     V_b(r_{ij})      &=& k \, f^b_n(r_{ij}) \\
     V_a({\theta_{ijk}})       &=& k \, f^a_n({\theta_{ijk}}) \\
     V_d(\phi_{ijkl}) &=& k \, f^d_n(\phi_{ijkl})\end{aligned}

   where :math:`k` is a force constant in units of energy and :math:`f`
  is a cubic spline function; for details see [subsec:cubicspline]. For
  each interaction, the force constant :math:`k` and the table number
  :math:`n` are specified in the topology. There are two different types
  of bonds, one that generates exclusions (type 8) and one that does not
  (type 9). For details see Table [tab:topfile2]. The table files are
  supplied to the mdrun program. After the table file name an
  underscore, the letter “b” for bonds, “a” for angles or “d” for
  dihedrals and the table number must be appended. For example, a
  tabulated bond with :math:`n=0` can be read from the file
  table\_b0.xvg. Multiple tables can be supplied simply by adding files
  with different values of :math:`n`, and are applied to the appropriate
  bonds, as specified in the topology (Table [tab:topfile2]). The format
  for the table files is three fixed-format columns of any suitable
  width. These columns must contain :math:`x`, :math:`f(x)`,
  :math:`-f'(x)`, and the values of :math:`x` should be uniformly
  spaced. Requirements for entries in the topology are given
  in Table [tab:topfile2]. The setup of the tables is as follows:
| **bonds**: :math:`x` is the distance in nm. For distances beyond the
  table length, mdrun will quit with an error message.
| **angles**: :math:`x` is the angle in degrees. The table should go
  from 0 up to and including 180 degrees; the derivative is taken in
  degrees.
| **dihedrals**: :math:`x` is the dihedral angle in degrees. The table
  should go from -180 up to and including 180 degrees; the IUPAC/IUB
  convention is used, *i.e.* zero is cis, the derivative is taken in
  degrees.

Restraints
----------

Special potentials are used for imposing restraints on the motion of the
system, either to avoid disastrous deviations, or to include knowledge
from experimental data. In either case they are not really part of the
force field and the reliability of the parameters is not important. The
potential forms, as implemented in GROMACS, are mentioned just for the
sake of completeness. Restraints and constraints refer to quite
different algorithms in GROMACS.

Position restraints
~~~~~~~~~~~~~~~~~~~

These are used to restrain particles to fixed reference positions
:math:`{\mbox{\boldmath ${R}$}}_i`. They can be used during
equilibration in order to avoid drastic rearrangements of critical parts
(*e.g.* to restrain motion in a protein that is subjected to large
solvent forces when the solvent is not yet equilibrated). Another
application is the restraining of particles in a shell around a region
that is simulated in detail, while the shell is only approximated
because it lacks proper interaction from missing particles outside the
shell. Restraining will then maintain the integrity of the inner part.
For spherical shells, it is a wise procedure to make the force constant
depend on the radius, increasing from zero at the inner boundary to a
large value at the outer boundary. This feature has not, however, been
implemented in GROMACS.

The following form is used:

.. math:: V_{pr}({\mbox{\boldmath ${r}$}}_i) = {\frac{1}{2}}k_{pr}|{{\mbox{\boldmath ${r}$}}_i}-{\mbox{\boldmath ${R}$}}_i|^2

 The potential is plotted in Fig. [fig:positionrestraint].

.. figure:: plots/f-pr
   :alt: Position restraint potential.
   :width: 8.00000cm

   Position restraint potential.

The potential form can be rewritten without loss of generality as:

.. math:: V_{pr}({\mbox{\boldmath ${r}$}}_i) = {\frac{1}{2}} \left[ k_{pr}^x (x_i-X_i)^2 ~{\hat{\bf x}} + k_{pr}^y (y_i-Y_i)^2 ~{\hat{\bf y}} + k_{pr}^z (z_i-Z_i)^2 ~{\hat{\bf z}}\right]

Now the forces are:

.. math::

   \begin{array}{rcl}
   F_i^x &=& -k_{pr}^x~(x_i - X_i) \\
   F_i^y &=& -k_{pr}^y~(y_i - Y_i) \\
   F_i^z &=& -k_{pr}^z~(z_i - Z_i)
   \end{array}

 Using three different force constants the position restraints can be
turned on or off in each spatial dimension; this means that atoms can be
harmonically restrained to a plane or a line. Position restraints are
applied to a special fixed list of atoms. Such a list is usually
generated by the pdb2gmx program.

Flat-bottomed position restraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Flat-bottomed position restraints can be used to restrain particles to
part of the simulation volume. No force acts on the restrained particle
within the flat-bottomed region of the potential, however a harmonic
force acts to move the particle to the flat-bottomed region if it is
outside it. It is possible to apply normal and flat-bottomed position
restraints on the same particle (however, only with the same reference
position :math:`{\mbox{\boldmath ${R}$}}_i`). The following general
potential is used (Figure [fig:fbposres]A):

.. math:: V_\mathrm{fb}({\mbox{\boldmath ${r}$}}_i) = \frac{1}{2}k_\mathrm{fb} [d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) - r_\mathrm{fb}]^2\,H[d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) - r_\mathrm{fb}],

 where :math:`{\mbox{\boldmath ${R}$}}_i` is the reference position,
:math:`r_\mathrm{fb}` is the distance from the center with a flat
potential, :math:`k_\mathrm{fb}` the force constant, and :math:`H` is
the Heaviside step function. The distance
:math:`d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i)` from
the reference position depends on the geometry :math:`g` of the
flat-bottomed potential.

.. figure:: plots/fbposres
   :alt: Flat-bottomed position restraint potential. (A) Not inverted,
   (B) inverted.
   :width: 10.00000cm

   Flat-bottomed position restraint potential. (A) Not inverted, (B)
   inverted.

| The following geometries for the flat-bottomed potential are
  supported: (:math:`g =1`): The particle is kept in a sphere of given
  radius. The force acts towards the center of the sphere. The following
  distance calculation is used:

  .. math:: d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) = |{\mbox{\boldmath ${r}$}}_i-{\mbox{\boldmath ${R}$}}_i|

   **Cylinder** (:math:`g=6,7,8`): The particle is kept in a cylinder of
  given radius parallel to the :math:`x` (:math:`g=6`), :math:`y`
  (:math:`g=7`), or :math:`z`-axis (:math:`g=8`). For backwards
  compatibility, setting :math:`g=2` is mapped to :math:`g=8` in the
  code so that old .tpr files and topologies work. The force from the
  flat-bottomed potential acts towards the axis of the cylinder. The
  component of the force parallel to the cylinder axis is zero. For a
  cylinder aligned along the :math:`z`-axis:

  .. math:: d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) = \sqrt{ (x_i-X_i)^2 + (y_i - Y_i)^2 }

   **Layer** (:math:`g=3,4,5`): The particle is kept in a layer defined
  by the thickness and the normal of the layer. The layer normal can be
  parallel to the :math:`x`, :math:`y`, or :math:`z`-axis. The force
  acts parallel to the layer normal.
| 

  .. math::

     d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) = |x_i-X_i|, \;\;\;\mbox{or}\;\;\; 
      d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) = |y_i-Y_i|, \;\;\;\mbox{or}\;\;\; 
     d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) = |z_i-Z_i|.

It is possible to apply multiple independent flat-bottomed position
restraints of different geometry on one particle. For example, applying
a cylinder and a layer in :math:`z` keeps a particle within a disk.
Applying three layers in :math:`x`, :math:`y`, and :math:`z` keeps the
particle within a cuboid.

In addition, it is possible to invert the restrained region with the
unrestrained region, leading to a potential that acts to keep the
particle *outside* of the volume defined by
:math:`{\mbox{\boldmath ${R}$}}_i`, :math:`g`, and
:math:`r_\mathrm{fb}`. That feature is switched on by defining a
negative :math:`r_\mathrm{fb}` in the topology. The following potential
is used (Figure [fig:fbposres]B):

.. math::

   V_\mathrm{fb}^{\mathrm{inv}}({\mbox{\boldmath ${r}$}}_i) = \frac{1}{2}k_\mathrm{fb}
     [d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) - |r_\mathrm{fb}|]^2\,
     H[ -(d_g({\mbox{\boldmath ${r}$}}_i;{\mbox{\boldmath ${R}$}}_i) - |r_\mathrm{fb}|)].

Angle restraints
~~~~~~~~~~~~~~~~

These are used to restrain the angle between two pairs of particles or
between one pair of particles and the :math:`z`-axis. The functional
form is similar to that of a proper dihedral. For two pairs of atoms:

.. math::

   V_{ar}({\mbox{\boldmath ${r}$}}_i,{\mbox{\boldmath ${r}$}}_j,{\mbox{\boldmath ${r}$}}_k,{\mbox{\boldmath ${r}$}}_l)
           = k_{ar}(1 - \cos(n (\theta - \theta_0))
           )
   ,~~~~\mbox{where}~~
   \theta = \arccos\left(\frac{{\mbox{\boldmath ${r}$}}_j -{\mbox{\boldmath ${r}$}}_i}{\|{\mbox{\boldmath ${r}$}}_j -{\mbox{\boldmath ${r}$}}_i\|}
    \cdot \frac{{\mbox{\boldmath ${r}$}}_l -{\mbox{\boldmath ${r}$}}_k}{\|{\mbox{\boldmath ${r}$}}_l -{\mbox{\boldmath ${r}$}}_k\|} \right)

 For one pair of atoms and the :math:`z`-axis:

.. math::

   V_{ar}({\mbox{\boldmath ${r}$}}_i,{\mbox{\boldmath ${r}$}}_j) = k_{ar}(1 - \cos(n (\theta - \theta_0))
           )
   ,~~~~\mbox{where}~~
   \theta = \arccos\left(\frac{{\mbox{\boldmath ${r}$}}_j -{\mbox{\boldmath ${r}$}}_i}{\|{\mbox{\boldmath ${r}$}}_j -{\mbox{\boldmath ${r}$}}_i\|}
    \cdot \left( \begin{array}{c} 0 \\ 0 \\ 1 \\ \end{array} \right) \right)

 A multiplicity (:math:`n`) of 2 is useful when you do not want to
distinguish between parallel and anti-parallel vectors. The equilibrium
angle :math:`\theta` should be between 0 and 180 degrees for
multiplicity 1 and between 0 and 90 degrees for multiplicity 2.

Dihedral restraints
~~~~~~~~~~~~~~~~~~~

These are used to restrain the dihedral angle :math:`\phi` defined by
four particles as in an improper dihedral (sec. [sec:imp]) but with a
slightly modified potential. Using:

.. math::

   \phi' = \left(\phi-\phi_0\right) ~{\rm MOD}~ 2\pi
   \label{eqn:dphi}

 where :math:`\phi_0` is the reference angle, the potential is defined
as:

.. math::

   V_{dihr}(\phi') ~=~ \left\{
   \begin{array}{lcllll}
   {\frac{1}{2}}k_{dihr}(\phi'-\phi_0-\Delta\phi)^2      
                   &\mbox{for}&     \phi' & >   & \Delta\phi       \\[1.5ex]
   0               &\mbox{for}&     \phi' & \le & \Delta\phi       \\[1.5ex]
   \end{array}\right.
   \label{eqn:dihre}

 where :math:`\Delta\phi` is a user defined angle and :math:`k_{dihr}`
is the force constant. **Note** that in the input in topology files,
angles are given in degrees and force constants in
kJ/mol/rad\ :math:`^2`.

Distance restraints
~~~~~~~~~~~~~~~~~~~

Distance restraints add a penalty to the potential when the distance
between specified pairs of atoms exceeds a threshold value. They are
normally used to impose experimental restraints from, for instance,
experiments in nuclear magnetic resonance (NMR), on the motion of the
system. Thus, MD can be used for structure refinement using NMR data. In
GROMACS there are three ways to impose restraints on pairs of atoms:

-  Simple harmonic restraints: use type 6 (see sec. [sec:excl]).

-  [subsec:harmonicrestraint]Piecewise linear/harmonic restraints: type
   10.

-  Complex NMR distance restraints, optionally with pair, time and/or
   ensemble averaging.

The last two options will be detailed now.

The potential form for distance restraints is quadratic below a
specified lower bound and between two specified upper bounds, and linear
beyond the largest bound (see Fig. [fig:dist]).

.. math::

   V_{dr}(r_{ij}) ~=~ \left\{
   \begin{array}{lcllllll}
   {\frac{1}{2}}k_{dr}(r_{ij}-r_0)^2      
                   &\mbox{for}&     &     & r_{ij} & < & r_0       \\[1.5ex]
   0               &\mbox{for}& r_0 & \le & r_{ij} & < & r_1       \\[1.5ex]
   {\frac{1}{2}}k_{dr}(r_{ij}-r_1)^2      
                   &\mbox{for}& r_1 & \le & r_{ij} & < & r_2       \\[1.5ex]
   {\frac{1}{2}}k_{dr}(r_2-r_1)(2r_{ij}-r_2-r_1)  
                   &\mbox{for}& r_2 & \le & r_{ij} &   &
   \end{array}\right.
   \label{eqn:disre}

.. figure:: plots/f-dr
   :alt: Distance Restraint potential.
   :width: 8.00000cm

   Distance Restraint potential.

The forces are

.. math::

   {\mbox{\boldmath ${F}$}}_i~=~ \left\{
   \begin{array}{lcllllll}
   -k_{dr}(r_{ij}-r_0)\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}} 
                   &\mbox{for}&     &     & r_{ij} & < & r_0       \\[1.5ex]
   0               &\mbox{for}& r_0 & \le & r_{ij} & < & r_1       \\[1.5ex]
   -k_{dr}(r_{ij}-r_1)\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}} 
                   &\mbox{for}& r_1 & \le & r_{ij} & < & r_2       \\[1.5ex]
   -k_{dr}(r_2-r_1)\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}    
                   &\mbox{for}& r_2 & \le & r_{ij} &   &
   \end{array} \right.

For restraints not derived from NMR data, this functionality will
usually suffice and a section of type 10 can be used to apply individual
restraints between pairs of atoms, see [subsec:topfile]. For applying
restraints derived from NMR measurements, more complex functionality
might be required, which is provided through the section and is
described below.

Time averaging
^^^^^^^^^^^^^^

Distance restraints based on instantaneous distances can potentially
reduce the fluctuations in a molecule significantly. This problem can be
overcome by restraining to a *time averaged*
distance \ `91 <#ref-Torda89>`__. The forces with time averaging are:

.. math::

   {\mbox{\boldmath ${F}$}}_i~=~ \left\{
   \begin{array}{lcllllll}
   -k^a_{dr}(\bar{r}_{ij}-r_0)\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}   
                   &\mbox{for}&     &     & \bar{r}_{ij} & < & r_0 \\[1.5ex]
   0               &\mbox{for}& r_0 & \le & \bar{r}_{ij} & < & r_1 \\[1.5ex]
   -k^a_{dr}(\bar{r}_{ij}-r_1)\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}   
                   &\mbox{for}& r_1 & \le & \bar{r}_{ij} & < & r_2 \\[1.5ex]
   -k^a_{dr}(r_2-r_1)\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}    
                   &\mbox{for}& r_2 & \le & \bar{r}_{ij} &   &
   \end{array} \right.

 where :math:`\bar{r}_{ij}` is given by an exponential running average
with decay time :math:`\tau`:

.. math::

   \bar{r}_{ij} ~=~ < r_{ij}^{-3} >^{-1/3}
   \label{eqn:rav}

 The force constant :math:`k^a_{dr}` is switched on slowly to compensate
for the lack of history at the beginning of the simulation:

.. math:: k^a_{dr} = k_{dr} \left(1-\exp\left(-\frac{t}{\tau}\right)\right)

 Because of the time averaging, we can no longer speak of a distance
restraint potential.

This way an atom can satisfy two incompatible distance restraints *on
average* by moving between two positions. An example would be an amino
acid side-chain that is rotating around its :math:`\chi` dihedral angle,
thereby coming close to various other groups. Such a mobile side chain
can give rise to multiple NOEs that can not be fulfilled by a single
structure.

The computation of the time averaged distance in the mdrun program is
done in the following fashion:

.. math::

   \begin{array}{rcl}
   \overline{r^{-3}}_{ij}(0)       &=& r_{ij}(0)^{-3}      \\
   \overline{r^{-3}}_{ij}(t)       &=& \overline{r^{-3}}_{ij}(t-\Delta t)~\exp{\left(-\frac{\Delta t}{\tau}\right)} + r_{ij}(t)^{-3}\left[1-\exp{\left(-\frac{\Delta t}{\tau}\right)}\right]
   \label{eqn:ravdisre}
   \end{array}

When a pair is within the bounds, it can still feel a force because the
time averaged distance can still be beyond a bound. To prevent the
protons from being pulled too close together, a mixed approach can be
used. In this approach, the penalty is zero when the instantaneous
distance is within the bounds, otherwise the violation is the square
root of the product of the instantaneous violation and the time averaged
violation:

.. math::

   {\mbox{\boldmath ${F}$}}_i~=~ \left\{
   \begin{array}{lclll}
   k^a_{dr}\sqrt{(r_{ij}-r_0)(\bar{r}_{ij}-r_0)}\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}   
       & \mbox{for} & r_{ij} < r_0 & \mbox{and} & \bar{r}_{ij} < r_0 \\[1.5ex]
   -k^a _{dr} \,
     \mbox{min}\left(\sqrt{(r_{ij}-r_1)(\bar{r}_{ij}-r_1)},r_2-r_1\right)
     \frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}   
       & \mbox{for} & r_{ij} > r_1 & \mbox{and} & \bar{r}_{ij} > r_1 \\[1.5ex]
   0               &\mbox{otherwise}
   \end{array} \right.

Averaging over multiple pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes it is unclear from experimental data which atom pair gives
rise to a single NOE, in other occasions it can be obvious that more
than one pair contributes due to the symmetry of the system, *e.g.* a
methyl group with three protons. For such a group, it is not possible to
distinguish between the protons, therefore they should all be taken into
account when calculating the distance between this methyl group and
another proton (or group of protons). Due to the physical nature of
magnetic resonance, the intensity of the NOE signal is inversely
proportional to the sixth power of the inter-atomic distance. Thus, when
combining atom pairs, a fixed list of :math:`N` restraints may be taken
together, where the apparent “distance” is given by:

.. math::

   r_N(t) = \left [\sum_{n=1}^{N} \bar{r}_{n}(t)^{-6} \right]^{-1/6}
   \label{eqn:rsix}

 where we use :math:`r_{ij}` or eqn. [eqn:rav] for the
:math:`\bar{r}_{n}`. The :math:`r_N` of the instantaneous and
time-averaged distances can be combined to do a mixed restraining, as
indicated above. As more pairs of protons contribute to the same NOE
signal, the intensity will increase, and the summed “distance” will be
shorter than any of its components due to the reciprocal summation.

There are two options for distributing the forces over the atom pairs.
In the conservative option, the force is defined as the derivative of
the restraint potential with respect to the coordinates. This results in
a conservative potential when time averaging is not used. The force
distribution over the pairs is proportional to :math:`r^{-6}`. This
means that a close pair feels a much larger force than a distant pair,
which might lead to a molecule that is “too rigid.” The other option is
an equal force distribution. In this case each pair feels :math:`1/N` of
the derivative of the restraint potential with respect to :math:`r_N`.
The advantage of this method is that more conformations might be
sampled, but the non-conservative nature of the forces can lead to local
heating of the protons.

It is also possible to use *ensemble averaging* using multiple (protein)
molecules. In this case the bounds should be lowered as in:

.. math::

   \begin{array}{rcl}
   r_1     &~=~&   r_1 * M^{-1/6}  \\
   r_2     &~=~&   r_2 * M^{-1/6}
   \end{array}

 where :math:`M` is the number of molecules. The GROMACS preprocessor
grompp can do this automatically when the appropriate option is given.
The resulting “distance” is then used to calculate the scalar force
according to:

.. math::

   {\mbox{\boldmath ${F}$}}_i~=~\left\{
   \begin{array}{rcl}
   ~& 0 \hspace{4cm}  & r_{N} < r_1         \\
    & k_{dr}(r_{N}-r_1)\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}} & r_1 \le r_{N} < r_2 \\
    & k_{dr}(r_2-r_1)\frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{r_{ij}}    & r_{N} \ge r_2 
   \end{array} \right.

 where :math:`i` and :math:`j` denote the atoms of all the pairs that
contribute to the NOE signal.

Using distance restraints
^^^^^^^^^^^^^^^^^^^^^^^^^

A list of distance restrains based on NOE data can be added to a
molecule definition in your topology file, like in the following
example:

::

    [ distance_restraints ]
    ; ai   aj   type   index   type'      low     up1     up2     fac
    10     16      1       0       1      0.0     0.3     0.4     1.0
    10     28      1       1       1      0.0     0.3     0.4     1.0
    10     46      1       1       1      0.0     0.3     0.4     1.0
    16     22      1       2       1      0.0     0.3     0.4     2.5
    16     34      1       3       1      0.0     0.5     0.6     1.0

In this example a number of features can be found. In columns ai and aj
you find the atom numbers of the particles to be restrained. The type
column should always be 1. As explained in  [subsec:distancerestraint],
multiple distances can contribute to a single NOE signal. In the
topology this can be set using the index column. In our example, the
restraints 10-28 and 10-46 both have index 1, therefore they are treated
simultaneously. An extra requirement for treating restraints together is
that the restraints must be on successive lines, without any other
intervening restraint. The type’ column will usually be 1, but can be
set to 2 to obtain a distance restraint that will never be time- and
ensemble-averaged; this can be useful for restraining hydrogen bonds.
The columns low, up1, and up2 hold the values of :math:`r_0`,
:math:`r_1`, and :math:`r_2` from  eqn. [eqn:disre]. In some cases it
can be useful to have different force constants for some restraints;
this is controlled by the column fac. The force constant in the
parameter file is multiplied by the value in the column fac for each
restraint. Information for each restraint is stored in the energy file
and can be processed and plotted with gmx nmr.

Orientation restraints
~~~~~~~~~~~~~~~~~~~~~~

This section describes how orientations between vectors, as measured in
certain NMR experiments, can be calculated and restrained in MD
simulations. The presented refinement methodology and a comparison of
results with and without time and ensemble averaging have been
published \ `92 <#ref-Hess2003>`__.

Theory
^^^^^^

In an NMR experiment, orientations of vectors can be measured when a
molecule does not tumble completely isotropically in the solvent. Two
examples of such orientation measurements are residual dipolar couplings
(between two nuclei) or chemical shift anisotropies. An observable for a
vector :math:`{\mbox{\boldmath ${r}$}}_i` can be written as follows:

.. math:: \delta_i = \frac{2}{3} \mbox{tr}({{\mathbf S}}{{\mathbf D}}_i)

 where :math:`{{\mathbf S}}` is the dimensionless order tensor of the
molecule. The tensor :math:`{{\mathbf D}}_i` is given by:

.. math::

   \label{orient_def}
   {{\mathbf D}}_i = \frac{c_i}{\|{\mbox{\boldmath ${r}$}}_i\|^\alpha} \left(
   \begin{array}{lll}
   3 x x - 1 & 3 x y     & 3 x z     \\
   3 x y     & 3 y y - 1 & 3 y z     \\
   3 x z     & 3 y z     & 3 z z - 1 \\
   \end{array} \right)

.. math::

   \mbox{with:} \quad 
   x=\frac{r_{i,x}}{\|{\mbox{\boldmath ${r}$}}_i\|}, \quad
   y=\frac{r_{i,y}}{\|{\mbox{\boldmath ${r}$}}_i\|}, \quad 
   z=\frac{r_{i,z}}{\|{\mbox{\boldmath ${r}$}}_i\|}

 For a dipolar coupling :math:`{\mbox{\boldmath ${r}$}}_i` is the vector
connecting the two nuclei, :math:`\alpha=3` and the constant :math:`c_i`
is given by:

.. math:: c_i = \frac{\mu_0}{4\pi} \gamma_1^i \gamma_2^i \frac{\hbar}{4\pi}

 where :math:`\gamma_1^i` and :math:`\gamma_2^i` are the gyromagnetic
ratios of the two nuclei.

The order tensor is symmetric and has trace zero. Using a rotation
matrix :math:`{\mathbf T}` it can be transformed into the following
form:

.. math::

   {\mathbf T}^T {{\mathbf S}}{\mathbf T} = s \left( \begin{array}{ccc}
   -\frac{1}{2}(1-\eta) & 0                    & 0 \\
   0                    & -\frac{1}{2}(1+\eta) & 0 \\
   0                    & 0                    & 1
   \end{array} \right)

 where :math:`-1 \leq s \leq 1` and :math:`0 \leq \eta \leq 1`.
:math:`s` is called the order parameter and :math:`\eta` the asymmetry
of the order tensor :math:`{{\mathbf S}}`. When the molecule tumbles
isotropically in the solvent, :math:`s` is zero, and no orientational
effects can be observed because all :math:`\delta_i` are zero.

Calculating orientations in a simulation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For reasons which are explained below, the :math:`{{\mathbf D}}`
matrices are calculated which respect to a reference orientation of the
molecule. The orientation is defined by a rotation matrix
:math:`{{\mathbf R}}`, which is needed to least-squares fit the current
coordinates of a selected set of atoms onto a reference conformation.
The reference conformation is the starting conformation of the
simulation. In case of ensemble averaging, which will be treated later,
the structure is taken from the first subsystem. The calculated
:math:`{{\mathbf D}}_i^c` matrix is given by:

.. math::

   \label{D_rot}
   {{\mathbf D}}_i^c(t) = {{\mathbf R}}(t) {{\mathbf D}}_i(t) {{\mathbf R}}^T(t)

 The calculated orientation for vector :math:`i` is given by:

.. math:: \delta^c_i(t) = \frac{2}{3} \mbox{tr}({{\mathbf S}}(t){{\mathbf D}}_i^c(t))

 The order tensor :math:`{{\mathbf S}}(t)` is usually unknown. A
reasonable choice for the order tensor is the tensor which minimizes the
(weighted) mean square difference between the calculated and the
observed orientations:

.. math::

   \label{S_msd}
   MSD(t) = \left(\sum_{i=1}^N w_i\right)^{-1} \sum_{i=1}^N w_i (\delta_i^c (t) -\delta_i^{exp})^2

 To properly combine different types of measurements, the unit of
:math:`w_i` should be such that all terms are dimensionless. This means
the unit of :math:`w_i` is the unit of :math:`\delta_i` to the power
:math:`-2`. **Note** that scaling all :math:`w_i` with a constant factor
does not influence the order tensor.

Time averaging
^^^^^^^^^^^^^^

Since the tensors :math:`{{\mathbf D}}_i` fluctuate rapidly in time,
much faster than can be observed in an experiment, they should be
averaged over time in the simulation. However, in a simulation the time
and the number of copies of a molecule are limited. Usually one can not
obtain a converged average of the :math:`{{\mathbf D}}_i` tensors over
all orientations of the molecule. If one assumes that the average
orientations of the :math:`{\mbox{\boldmath ${r}$}}_i` vectors within
the molecule converge much faster than the tumbling time of the
molecule, the tensor can be averaged in an axis system that rotates with
the molecule, as expressed by equation ([D\_rot]). The time-averaged
tensors are calculated using an exponentially decaying memory function:

.. math::

   {{\mathbf D}}^a_i(t) = \frac{\displaystyle
   \int_{u=t_0}^t {{\mathbf D}}^c_i(u) \exp\left(-\frac{t-u}{\tau}\right)\mbox{d} u
   }{\displaystyle
   \int_{u=t_0}^t \exp\left(-\frac{t-u}{\tau}\right)\mbox{d} u
   }

 Assuming that the order tensor :math:`{{\mathbf S}}` fluctuates slower
than the :math:`{{\mathbf D}}_i`, the time-averaged orientation can be
calculated as:

.. math:: \delta_i^a(t) = \frac{2}{3} \mbox{tr}({{\mathbf S}}(t) {{\mathbf D}}_i^a(t))

 where the order tensor :math:`{{\mathbf S}}(t)` is calculated using
expression ([S\_msd]) with :math:`\delta_i^c(t)` replaced by
:math:`\delta_i^a(t)`.

Restraining
^^^^^^^^^^^

The simulated structure can be restrained by applying a force
proportional to the difference between the calculated and the
experimental orientations. When no time averaging is applied, a proper
potential can be defined as:

.. math:: V = \frac{1}{2} k \sum_{i=1}^N w_i (\delta_i^c (t) -\delta_i^{exp})^2

 where the unit of :math:`k` is the unit of energy. Thus the effective
force constant for restraint :math:`i` is :math:`k w_i`. The forces are
given by minus the gradient of :math:`V`. The force
:math:`{\mbox{\boldmath ${F}$}}\!_i` working on vector
:math:`{\mbox{\boldmath ${r}$}}_i` is:

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${F}$}}\!_i(t) 
   & = & - \frac{\mbox{d} V}{\mbox{d}{\mbox{\boldmath ${r}$}}_i} \\
   & = & -k w_i (\delta_i^c (t) -\delta_i^{exp}) \frac{\mbox{d} \delta_i (t)}{\mbox{d}{\mbox{\boldmath ${r}$}}_i} \\
   & = & -k w_i (\delta_i^c (t) -\delta_i^{exp})
   \frac{2 c_i}{\|{\mbox{\boldmath ${r}$}}\|^{2+\alpha}} \left(2 {{\mathbf R}}^T {{\mathbf S}}{{\mathbf R}}{\mbox{\boldmath ${r}$}}_i - \frac{2+\alpha}{\|{\mbox{\boldmath ${r}$}}\|^2} \mbox{tr}({{\mathbf R}}^T {{\mathbf S}}{{\mathbf R}}{\mbox{\boldmath ${r}$}}_i {\mbox{\boldmath ${r}$}}_i^T) {\mbox{\boldmath ${r}$}}_i \right)\end{aligned}

Ensemble averaging
^^^^^^^^^^^^^^^^^^

Ensemble averaging can be applied by simulating a system of :math:`M`
subsystems that each contain an identical set of orientation restraints.
The systems only interact via the orientation restraint potential which
is defined as:

.. math::

   V = M \frac{1}{2} k \sum_{i=1}^N w_i 
   \langle \delta_i^c (t) -\delta_i^{exp} \rangle^2

 The force on vector :math:`{\mbox{\boldmath ${r}$}}_{i,m}` in subsystem
:math:`m` is given by:

.. math::

   {\mbox{\boldmath ${F}$}}\!_{i,m}(t) = - \frac{\mbox{d} V}{\mbox{d}{\mbox{\boldmath ${r}$}}_{i,m}} =
   -k w_i \langle \delta_i^c (t) -\delta_i^{exp} \rangle \frac{\mbox{d} \delta_{i,m}^c (t)}{\mbox{d}{\mbox{\boldmath ${r}$}}_{i,m}} \\

Time averaging
^^^^^^^^^^^^^^

When using time averaging it is not possible to define a potential. We
can still define a quantity that gives a rough idea of the energy stored
in the restraints:

.. math::

   V = M \frac{1}{2} k^a \sum_{i=1}^N w_i 
   \langle \delta_i^a (t) -\delta_i^{exp} \rangle^2

 The force constant :math:`k_a` is switched on slowly to compensate for
the lack of history at times close to :math:`t_0`. It is exactly
proportional to the amount of average that has been accumulated:

.. math::

   k^a =
    k \, \frac{1}{\tau}\int_{u=t_0}^t \exp\left(-\frac{t-u}{\tau}\right)\mbox{d} u

 What really matters is the definition of the force. It is chosen to be
proportional to the square root of the product of the time-averaged and
the instantaneous deviation. Using only the time-averaged deviation
induces large oscillations. The force is given by:

.. math::

   {\mbox{\boldmath ${F}$}}\!_{i,m}(t) =
   \left\{ \begin{array}{ll}
   0 & \quad \mbox{for} \quad a\, b \leq 0 \\
   \displaystyle
   k^a w_i \frac{a}{|a|} \sqrt{a\, b} \, \frac{\mbox{d} \delta_{i,m}^c (t)}{\mbox{d}{\mbox{\boldmath ${r}$}}_{i,m}}
   & \quad \mbox{for} \quad a\, b > 0 
   \end{array}
   \right.

.. math::

   \begin{aligned}
   a &=& \langle \delta_i^a (t) -\delta_i^{exp} \rangle \\
   b &=& \langle \delta_i^c (t) -\delta_i^{exp} \rangle\end{aligned}

Using orientation restraints
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Orientation restraints can be added to a molecule definition in the
topology file in the section . Here we give an example section
containing five N-H residual dipolar coupling restraints:

::

    [ orientation_restraints ]
    ; ai   aj  type  exp.  label  alpha    const.     obs.   weight
    ;                                Hz      nm^3       Hz    Hz^-2
      31   32     1     1      3      3     6.083    -6.73      1.0
      43   44     1     1      4      3     6.083    -7.87      1.0
      55   56     1     1      5      3     6.083    -7.13      1.0
      65   66     1     1      6      3     6.083    -2.57      1.0
      73   74     1     1      7      3     6.083    -2.10      1.0

The unit of the observable is Hz, but one can choose any other unit. In
columns ai and aj you find the atom numbers of the particles to be
restrained. The type column should always be 1. The exp. column denotes
the experiment number, starting at 1. For each experiment a separate
order tensor :math:`{{\mathbf S}}` is optimized. The label should be a
unique number larger than zero for each restraint. The alpha column
contains the power :math:`\alpha` that is used in
equation ([orient\_def]) to calculate the orientation. The const. column
contains the constant :math:`c_i` used in the same equation. The
constant should have the unit of the observable times
nm\ :math:`^\alpha`. The column obs. contains the observable, in any
unit you like. The last column contains the weights :math:`w_i`; the
unit should be the inverse of the square of the unit of the observable.

Some parameters for orientation restraints can be specified in the
grompp.mdp file, for a study of the effect of different force constants
and averaging times and ensemble averaging see \ `92 <#ref-Hess2003>`__.
Information for each restraint is stored in the energy file and can be
processed and plotted with gmx nmr.

Polarization
------------

Polarization can be treated by GROMACS by attaching shell (Drude)
particles to atoms and/or virtual sites. The energy of the shell
particle is then minimized at each time step in order to remain on the
Born-Oppenheimer surface.

Simple polarization
~~~~~~~~~~~~~~~~~~~

This is implemented as a harmonic potential with equilibrium distance 0.
The input given in the topology file is the polarizability
:math:`\alpha` (in GROMACS units) as follows:

::

    [ polarization ]
    ; Atom i  j  type  alpha
    1         2  1     0.001

in this case the polarizability volume is 0.001 nm\ :math:`^3` (or 1
Å\ :math:`^3`). In order to compute the harmonic force constant
:math:`k_{cs}` (where :math:`cs` stands for core-shell), the following
is used \ `45 <#ref-Maaren2001a>`__:

.. math:: k_{cs} ~=~ \frac{q_s^2}{\alpha}

 where :math:`q_s` is the charge on the shell particle.

Anharmonic polarization
~~~~~~~~~~~~~~~~~~~~~~~

For the development of the Drude force field by Roux and
McKerell \ `93 <#ref-Lopes2013a>`__ it was found that some particles can
overpolarize and this was fixed by introducing a higher order term in
the polarization energy:

.. math::

   \begin{aligned}
   V_{pol} ~=& \frac{k_{cs}}{2} r_{cs}^2 & r_{cs} \le \delta \\
               =& \frac{k_{cs}}{2} r_{cs}^2 + k_{hyp} (r_{cs}-\delta)^4 & r_{cs} > \delta\end{aligned}

 where :math:`\delta` is a user-defined constant that is set to 0.02 nm
for anions in the Drude force field \ `94 <#ref-HYu2010>`__. Since this
original introduction it has also been used in other atom
types \ `93 <#ref-Lopes2013a>`__.

::

    [ polarization ]
    ;Atom i j    type   alpha (nm^3)    delta  khyp
    1       2       2       0.001786     0.02  16.736e8

The above force constant :math:`k_{hyp}` corresponds to
4\ :math:`\cdot`\ 10\ :math:`^8` kcal/mol/nm\ :math:`^4`, hence the
strange number.

Water polarization
~~~~~~~~~~~~~~~~~~

A special potential for water that allows anisotropic polarization of a
single shell particle \ `45 <#ref-Maaren2001a>`__.

Thole polarization
~~~~~~~~~~~~~~~~~~

Based on early work by Thole `95 <#ref-Thole81>`__, Roux and coworkers
have implemented potentials for molecules like
ethanol \ `96 <#ref-Lamoureux2003a>`__\ `98 <#ref-Noskov2005a>`__.
Within such molecules, there are intra-molecular interactions between
shell particles, however these must be screened because full Coulomb
would be too strong. The potential between two shell particles :math:`i`
and :math:`j` is:

.. math:: V_{thole} ~=~ \frac{q_i q_j}{r_{ij}}\left[1-\left(1+\frac{{\bar{r}_{ij}}}{2}\right){\rm exp}^{-{\bar{r}_{ij}}}\right]

**Note** that there is a sign error in Equation 1 of Noskov *et
al.* `98 <#ref-Noskov2005a>`__:

.. math:: {\bar{r}_{ij}}~=~ a\frac{r_{ij}}{(\alpha_i \alpha_j)^{1/6}}

 where :math:`a` is a magic (dimensionless) constant, usually chosen to
be 2.6 \ `98 <#ref-Noskov2005a>`__; :math:`\alpha_i` and
:math:`\alpha_j` are the polarizabilities of the respective shell
particles.

Free energy interactions
------------------------

This section describes the :math:`\lambda`-dependence of the potentials
used for free energy calculations (see sec. [sec:fecalc]). All common
types of potentials and constraints can be interpolated smoothly from
state A (:math:`\lambda=0`) to state B (:math:`\lambda=1`) and vice
versa. All bonded interactions are interpolated by linear interpolation
of the interaction parameters. Non-bonded interactions can be
interpolated linearly or via soft-core interactions.

Starting in GROMACS 4.6, :math:`\lambda` is a vector, allowing different
components of the free energy transformation to be carried out at
different rates. Coulomb, Lennard-Jones, bonded, and restraint terms can
all be controlled independently, as described in the .mdp options.

Harmonic potentials
^^^^^^^^^^^^^^^^^^^

The example given here is for the bond potential, which is harmonic in
GROMACS. However, these equations apply to the angle potential and the
improper dihedral potential as well.

.. math::

   \begin{aligned}
   V_b     &=&{\frac{1}{2}}\left[{(1-{\lambda})}k_b^A + 
                   {\lambda}k_b^B\right] \left[b - {(1-{\lambda})}b_0^A - {\lambda}b_0^B\right]^2  \\
   {\frac{\partial V_b}{\partial {\lambda}}}&=&{\frac{1}{2}}(k_b^B-k_b^A)
                   \left[b - {(1-{\lambda})}b_0^A + {\lambda}b_0^B\right]^2 + 
   		\nonumber\\
           & & \phantom{{\frac{1}{2}}}(b_0^A-b_0^B) \left[b - {(1-{\lambda})}b_0^A -{\lambda}b_0^B\right]
   		\left[{(1-{\lambda})}k_b^A + {\lambda}k_b^B \right]\end{aligned}

GROMOS-96 bonds and angles
^^^^^^^^^^^^^^^^^^^^^^^^^^

Fourth-power bond stretching and cosine-based angle potentials are
interpolated by linear interpolation of the force constant and the
equilibrium position. Formulas are not given here.

Proper dihedrals
^^^^^^^^^^^^^^^^

For the proper dihedrals, the equations are somewhat more complicated:

.. math::

   \begin{aligned}
   V_d     &=&\left[{(1-{\lambda})}k_d^A + {\lambda}k_d^B \right]
           \left( 1+ \cos\left[n_{\phi} \phi - 
   		    {(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B
   		    \right]\right)\\
   {\frac{\partial V_d}{\partial {\lambda}}}&=&(k_d^B-k_d^A) 
            \left( 1+ \cos
   		 \left[
   		    n_{\phi} \phi- {(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B
   		 \right]
   	 \right) +
   	 \nonumber\\
           &&(\phi_s^B - \phi_s^A) \left[{(1-{\lambda})}k_d^A - {\lambda}k_d^B\right] 
           \sin\left[  n_{\phi}\phi - {(1-{\lambda})}\phi_s^A - {\lambda}\phi_s^B \right]\end{aligned}

 **Note:** that the multiplicity :math:`n_{\phi}` can not be
parameterized because the function should remain periodic on the
interval :math:`[0,2\pi]`.

Tabulated bonded interactions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For tabulated bonded interactions only the force constant can
interpolated:

.. math::

   \begin{aligned}
         V  &=& ({(1-{\lambda})}k^A + {\lambda}k^B) \, f \\
   {\frac{\partial V}{\partial {\lambda}}} &=& (k^B - k^A) \, f\end{aligned}

Coulomb interaction
^^^^^^^^^^^^^^^^^^^

The Coulomb interaction between two particles of which the charge varies
with :math:`{\lambda}` is:

.. math::

   \begin{aligned}
   V_c &=& \frac{f}{{\varepsilon_{rf}}{r_{ij}}}\left[{(1-{\lambda})}q_i^A q_j^A + {\lambda}\, q_i^B q_j^B\right] \\
   {\frac{\partial V_c}{\partial {\lambda}}}&=& \frac{f}{{\varepsilon_{rf}}{r_{ij}}}\left[- q_i^A q_j^A + q_i^B q_j^B\right]\end{aligned}

 where :math:`f = \frac{1}{4\pi \varepsilon_0} = {138.935\,458}` (see
chapter [ch:defunits]).

Coulomb interaction with reaction field
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The Coulomb interaction including a reaction field, between two
particles of which the charge varies with :math:`{\lambda}` is:

.. math::

   \begin{aligned}
   V_c     &=& f\left[\frac{1}{{r_{ij}}} + k_{rf}~ {r_{ij}}^2 -c_{rf}\right]
                \left[{(1-{\lambda})}q_i^A q_j^A + {\lambda}\, q_i^B q_j^B\right] \\
   {\frac{\partial V_c}{\partial {\lambda}}}&=& f\left[\frac{1}{{r_{ij}}} + k_{rf}~ {r_{ij}}^2 -c_{rf}\right]
                  \left[- q_i^A q_j^A + q_i^B q_j^B\right]
   	       \label{eq:dVcoulombdlambda}\end{aligned}

 **Note** that the constants :math:`k_{rf}` and :math:`c_{rf}` are
defined using the dielectric constant :math:`{\varepsilon_{rf}}` of the
medium (see sec. [sec:coulrf]).

Lennard-Jones interaction
^^^^^^^^^^^^^^^^^^^^^^^^^

For the Lennard-Jones interaction between two particles of which the
*atom type* varies with :math:`{\lambda}` we can write:

.. math::

   \begin{aligned}
   V_{LJ}  &=&     \frac{{(1-{\lambda})}C_{12}^A + {\lambda}\, C_{12}^B}{{r_{ij}}^{12}} -
                   \frac{{(1-{\lambda})}C_6^A + {\lambda}\, C_6^B}{{r_{ij}}^6}   \\
   {\frac{\partial V_{LJ}}{\partial {\lambda}}}&=&\frac{C_{12}^B - C_{12}^A}{{r_{ij}}^{12}} -
                   \frac{C_6^B - C_6^A}{{r_{ij}}^6}
   		\label{eq:dVljdlambda}\end{aligned}

 It should be noted that it is also possible to express a pathway from
state A to state B using :math:`\sigma` and :math:`\epsilon` (see
eqn. [eqn:sigeps]). It may seem to make sense physically to vary the
force field parameters :math:`\sigma` and :math:`\epsilon` rather than
the derived parameters :math:`C_{12}` and :math:`C_{6}`. However, the
difference between the pathways in parameter space is not large, and the
free energy itself does not depend on the pathway, so we use the simple
formulation presented above.

Kinetic Energy
^^^^^^^^^^^^^^

When the mass of a particle changes, there is also a contribution of the
kinetic energy to the free energy (note that we can not write the
momentum as m, since that would result in the sign of
:math:`{\frac{\partial E_k}{\partial {\lambda}}}` being
incorrect \ `99 <#ref-Gunsteren98a>`__):

.. math::

   \begin{aligned}
   E_k      &=&     {\frac{1}{2}}\frac{{\mbox{\boldmath ${p}$}}^2}{{(1-{\lambda})}m^A + {\lambda}m^B}        \\
   {\frac{\partial E_k}{\partial {\lambda}}}&=&    -{\frac{1}{2}}\frac{{\mbox{\boldmath ${p}$}}^2(m^B-m^A)}{({(1-{\lambda})}m^A + {\lambda}m^B)^2}\end{aligned}

after taking the derivative, we *can* insert = m, such that:

.. math:: {\frac{\partial E_k}{\partial {\lambda}}}~=~    -{\frac{1}{2}}{\mbox{\boldmath ${v}$}}^2(m^B-m^A)

Constraints
^^^^^^^^^^^

The constraints are formally part of the Hamiltonian, and therefore they
give a contribution to the free energy. In GROMACS this can be
calculated using the LINCS or the SHAKE algorithm. If we have
:math:`k = 1 \ldots K` constraint equations :math:`g_k` for LINCS, then

.. math:: g_k     =       |{\mbox{\boldmath ${r}$}}_{k}| - d_{k}

 where :math:`{\mbox{\boldmath ${r}$}}_k` is the displacement vector
between two particles and :math:`d_k` is the constraint distance between
the two particles. We can express the fact that the constraint distance
has a :math:`{\lambda}` dependency by

.. math:: d_k     =       {(1-{\lambda})}d_{k}^A + {\lambda}d_k^B

Thus the :math:`{\lambda}`-dependent constraint equation is

.. math:: g_k     =       |{\mbox{\boldmath ${r}$}}_{k}| - \left({(1-{\lambda})}d_{k}^A + {\lambda}d_k^B\right).

The (zero) contribution :math:`G` to the Hamiltonian from the
constraints (using Lagrange multipliers :math:`\lambda_k`, which are
logically distinct from the free-energy :math:`{\lambda}`) is

.. math::

   \begin{aligned}
   G           &=&     \sum^K_k \lambda_k g_k    \\
   {\frac{\partial G}{\partial {\lambda}}}    &=&     \frac{\partial G}{\partial d_k} {\frac{\partial d_k}{\partial {\lambda}}} \\
               &=&     - \sum^K_k \lambda_k \left(d_k^B-d_k^A\right)\end{aligned}

For SHAKE, the constraint equations are

.. math:: g_k     =       {\mbox{\boldmath ${r}$}}_{k}^2 - d_{k}^2

 with :math:`d_k` as before, so

.. math::

   \begin{aligned}
   {\frac{\partial G}{\partial {\lambda}}}    &=&     -2 \sum^K_k \lambda_k \left(d_k^B-d_k^A\right)\end{aligned}

Soft-core interactions
~~~~~~~~~~~~~~~~~~~~~~

.. figure:: plots/softcore
   :alt: Soft-core interactions at :math:`{\lambda}=0.5`, with
   :math:`p=2` and :math:`C_6^A=C_{12}^A=C_6^B=C_{12}^B=1`.
   :height: 6.00000cm

   Soft-core interactions at :math:`{\lambda}=0.5`, with :math:`p=2` and
   :math:`C_6^A=C_{12}^A=C_6^B=C_{12}^B=1`.

In a free-energy calculation where particles grow out of nothing, or
particles disappear, using the the simple linear interpolation of the
Lennard-Jones and Coulomb potentials as described in
Equations [eq:dVljdlambda] and [eq:dVcoulombdlambda] may lead to poor
convergence. When the particles have nearly disappeared, or are close to
appearing (at :math:`{\lambda}` close to 0 or 1), the interaction energy
will be weak enough for particles to get very close to each other,
leading to large fluctuations in the measured values of
:math:`\partial V/\partial {\lambda}` (which, because of the simple
linear interpolation, depends on the potentials at both the endpoints of
:math:`{\lambda}`).

To circumvent these problems, the singularities in the potentials need
to be removed. This can be done by modifying the regular Lennard-Jones
and Coulomb potentials with “soft-core” potentials that limit the
energies and forces involved at :math:`{\lambda}` values between 0 and
1, but not *at* :math:`{\lambda}=0` or 1.

In GROMACS the soft-core potentials :math:`V_{sc}` are shifted versions
of the regular potentials, so that the singularity in the potential and
its derivatives at :math:`r=0` is never reached:

.. math::

   \begin{aligned}
   V_{sc}(r) &=& {(1-{\lambda})}V^A(r_A) + {\lambda}V^B(r_B)
       \\
   r_A &=& \left(\alpha \sigma_A^6 {\lambda}^p + r^6 \right)^\frac{1}{6}
       \\
   r_B &=& \left(\alpha \sigma_B^6 {(1-{\lambda})}^p + r^6 \right)^\frac{1}{6}\end{aligned}

 where :math:`V^A` and :math:`V^B` are the normal “hard core” Van der
Waals or electrostatic potentials in state A (:math:`{\lambda}=0`) and
state B (:math:`{\lambda}=1`) respectively, :math:`\alpha` is the
soft-core parameter (set with sc\_alpha in the .mdp file), :math:`p` is
the soft-core :math:`{\lambda}` power (set with sc\_power),
:math:`\sigma` is the radius of the interaction, which is
:math:`(C_{12}/C_6)^{1/6}` or an input parameter (sc\_sigma) when
:math:`C_6` or :math:`C_{12}` is zero.

For intermediate :math:`{\lambda}`, :math:`r_A` and :math:`r_B` alter
the interactions very little for :math:`r > \alpha^{1/6} \sigma` and
quickly switch the soft-core interaction to an almost constant value for
smaller :math:`r` (Fig. [fig:softcore]). The force is:

.. math::

   F_{sc}(r) = -\frac{\partial V_{sc}(r)}{\partial r} =
    {(1-{\lambda})}F^A(r_A) \left(\frac{r}{r_A}\right)^5 +
   {\lambda}F^B(r_B) \left(\frac{r}{r_B}\right)^5

 where :math:`F^A` and :math:`F^B` are the “hard core” forces. The
contribution to the derivative of the free energy is:

.. math::

   \begin{aligned}
   {\frac{\partial V_{sc}(r)}{\partial {\lambda}}} & = &
    V^B(r_B) -V^A(r_A)  + 
   	{(1-{\lambda})}\frac{\partial V^A(r_A)}{\partial r_A}
   		   \frac{\partial r_A}{\partial {\lambda}} + 
   	{\lambda}\frac{\partial V^B(r_B)}{\partial r_B}
   		   \frac{\partial r_B}{\partial {\lambda}}
   \nonumber\\
   &=&
    V^B(r_B) -V^A(r_A)  + \nonumber \\
    & &
    \frac{p \alpha}{6}
          \left[ {\lambda}F^B(r_B) r^{-5}_B \sigma_B^6 {(1-{\lambda})}^{p-1} -
   	       {(1-{\lambda})}F^A(r_A) r^{-5}_A \sigma_A^6 {\lambda}^{p-1} \right]\end{aligned}

The original GROMOS Lennard-Jones soft-core
function \ `100 <#ref-Beutler94>`__ uses :math:`p=2`, but :math:`p=1`
gives a smoother :math:`\partial H/\partial{\lambda}` curve. Another
issue that should be considered is the soft-core effect of hydrogens
without Lennard-Jones interaction. Their soft-core :math:`\sigma` is set
with sc-sigma in the .mdp file. These hydrogens produce peaks in
:math:`\partial H/\partial{\lambda}` at :math:`{\lambda}` is 0 and/or 1
for :math:`p=1` and close to 0 and/or 1 with :math:`p=2`. Lowering will
decrease this effect, but it will also increase the interactions with
hydrogens relative to the other interactions in the soft-core state.

When soft-core potentials are selected (by setting sc-alpha >0), and the
Coulomb and Lennard-Jones potentials are turned on or off sequentially,
then the Coulombic interaction is turned off linearly, rather than using
soft-core interactions, which should be less statistically noisy in most
cases. This behavior can be overwritten by using the mdp option sc-coul
to yes. Note that the sc-coul is only taken into account when lambda
states are used, not with couple-lambda0 / couple-lambda1, and you can
still turn off soft-core interactions by setting sc-alpha=0.
Additionally, the soft-core interaction potential is only applied when
either the A or B state has zero interaction potential. If both A and B
states have nonzero interaction potential, default linear scaling
described above is used. When both Coulombic and Lennard-Jones
interactions are turned off simultaneously, a soft-core potential is
used, and a hydrogen is being introduced or deleted, the sigma is set to
sc-sigma-min, which itself defaults to sc-sigma-default.

Recently, a new formulation of the soft-core approach has been derived
that in most cases gives lower and more even statistical variance than
the standard soft-core path described above. \ `101 <#ref-Pham2011>`__,
`102 <#ref-Pham2012>`__ Specifically, we have:

.. math::

   \begin{aligned}
   V_{sc}(r) &=& {(1-{\lambda})}V^A(r_A) + {\lambda}V^B(r_B)
       \\
   r_A &=& \left(\alpha \sigma_A^{48} {\lambda}^p + r^{48} \right)^\frac{1}{48}
       \\
   r_B &=& \left(\alpha \sigma_B^{48} {(1-{\lambda})}^p + r^{48} \right)^\frac{1}{48}\end{aligned}

 This “1-1-48” path is also implemented in GROMACS. Note that for this
path the soft core :math:`\alpha` should satisfy
:math:`0.001 < \alpha < 0.003`, rather than :math:`\alpha \approx
0.5`.

Methods
-------

Exclusions and 1-4 Interactions.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Atoms within a molecule that are close by in the chain, *i.e.* atoms
that are covalently bonded, or linked by one or two atoms are called
*first neighbors, second neighbors* and *third neighbors*, respectively
(see Fig. [fig:chain]). Since the interactions of atom **i** with atoms
**i+1** and **i+2** are mainly quantum mechanical, they can not be
modeled by a Lennard-Jones potential. Instead it is assumed that these
interactions are adequately modeled by a harmonic bond term or
constraint (**i, i+1**) and a harmonic angle term (**i, i+2**). The
first and second neighbors (atoms **i+1** and **i+2**) are therefore
*excluded* from the Lennard-Jones interaction list of atom **i**; atoms
**i+1** and **i+2** are called *exclusions* of atom **i**.

.. figure:: plots/chain
   :alt: Atoms along an alkane chain.
   :width: 8.00000cm

   Atoms along an alkane chain.

For third neighbors, the normal Lennard-Jones repulsion is sometimes
still too strong, which means that when applied to a molecule, the
molecule would deform or break due to the internal strain. This is
especially the case for carbon-carbon interactions in a
*cis*-conformation (*e.g.* *cis*-butane). Therefore, for some of these
interactions, the Lennard-Jones repulsion has been reduced in the GROMOS
force field, which is implemented by keeping a separate list of 1-4 and
normal Lennard-Jones parameters. In other force fields, such as
OPLS \ `103 <#ref-Jorgensen88>`__, the standard Lennard-Jones parameters
are reduced by a factor of two, but in that case also the dispersion
(r:math:`^{-6}`) and the Coulomb interaction are scaled. GROMACS can use
either of these methods.

Charge Groups
~~~~~~~~~~~~~

In principle, the force calculation in MD is an :math:`O(N^2)` problem.
Therefore, we apply a cut-off for non-bonded force (NBF) calculations;
only the particles within a certain distance of each other are
interacting. This reduces the cost to :math:`O(N)` (typically
:math:`100N` to :math:`200N`) of the NBF. It also introduces an error,
which is, in most cases, acceptable, except when applying the cut-off
implies the creation of charges, in which case you should consider using
the lattice sum methods provided by GROMACS.

Consider a water molecule interacting with another atom. If we would
apply a plain cut-off on an atom-atom basis we might include the
atom-oxygen interaction (with a charge of :math:`-0.82`) without the
compensating charge of the protons, and as a result, induce a large
dipole moment over the system. Therefore, we have to keep groups of
atoms with total charge 0 together. These groups are called *charge
groups*. Note that with a proper treatment of long-range electrostatics
(e.g. particle-mesh Ewald (sec. [sec:pme]), keeping charge groups
together is not required.

Treatment of Cut-offs in the group scheme
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GROMACS is quite flexible in treating cut-offs, which implies there can
be quite a number of parameters to set. These parameters are set in the
input file for grompp. There are two sort of parameters that affect the
cut-off interactions; you can select which type of interaction to use in
each case, and which cut-offs should be used in the neighbor searching.

For both Coulomb and van der Waals interactions there are interaction
type selectors (termed vdwtype and coulombtype) and two parameters, for
a total of six non-bonded interaction parameters. See the User Guide for
a complete description of these parameters.

In the group cut-off scheme, all of the interaction functions in
Table [tab:funcparm] require that neighbor searching be done with a
radius at least as large as the :math:`r_c` specified for the functional
form, because of the use of charge groups. The extra radius is typically
of the order of 0.25 nm (roughly the largest distance between two atoms
in a charge group plus the distance a charge group can diffuse within
neighbor list updates).

Virtual interaction sites
-------------------------

Virtual interaction sites (called dummy atoms in GROMACS versions before
3.3) can be used in GROMACS in a number of ways. We write the position
of the virtual site :math:`{\mbox{\boldmath ${r}$}}_s` as a function of
the positions of other particles :math:`_i`:
:math:`{\mbox{\boldmath ${r}$}}_s =
f({\mbox{\boldmath ${r}$}}_1..{\mbox{\boldmath ${r}$}}_n)`. The virtual
site, which may carry charge or be involved in other interactions, can
now be used in the force calculation. The force acting on the virtual
site must be redistributed over the particles with mass in a consistent
way. A good way to do this can be found in
ref. \ `104 <#ref-Berendsen84b>`__. We can write the potential energy
as:

.. math:: V = V({{\mbox{\boldmath ${r}$}}_s},{\mbox{\boldmath ${r}$}}_1,\ldots,{\mbox{\boldmath ${r}$}}_n) = V^*({\mbox{\boldmath ${r}$}}_1,\ldots,{\mbox{\boldmath ${r}$}}_n)

 The force on the particle :math:`i` is then:

.. math::

   {\mbox{\boldmath ${F}$}}_i = -\frac{\partial V^*}{\partial {\mbox{\boldmath ${r}$}}_i} 
            = -\frac{\partial V}{\partial {\mbox{\boldmath ${r}$}}_i} - 
               \frac{\partial V}{\partial {{\mbox{\boldmath ${r}$}}_s}} 
               \frac{\partial {{\mbox{\boldmath ${r}$}}_s}}{\partial {\mbox{\boldmath ${r}$}}_i}
            = {\mbox{\boldmath ${F}$}}_i^{direct} + {{\mbox{\boldmath ${F}$}}_i'}

 The first term is the normal force. The second term is the force on
particle :math:`i` due to the virtual site, which can be written in
tensor notation:

.. math::

   {{\mbox{\boldmath ${F}$}}_i'}= \left[\begin{array}{ccc}
   {\displaystyle\frac{\partial x_s}{\partial x_i}} & {\displaystyle\frac{\partial y_s}{\partial x_i}} & {\displaystyle\frac{\partial z_s}{\partial x_i}}        \\[1ex]
   {\displaystyle\frac{\partial x_s}{\partial y_i}} & {\displaystyle\frac{\partial y_s}{\partial y_i}} & {\displaystyle\frac{\partial z_s}{\partial y_i}}        \\[1ex]
   {\displaystyle\frac{\partial x_s}{\partial z_i}} & {\displaystyle\frac{\partial y_s}{\partial z_i}} & {\displaystyle\frac{\partial z_s}{\partial z_i}}
   \end{array}\right]{{\mbox{\boldmath ${F}$}}_{s}}\label{eqn:fvsite}

where :math:`{{\mbox{\boldmath ${F}$}}_{s}}` is the force on the virtual
site and :math:`x_s`, :math:`y_s` and :math:`z_s` are the coordinates of
the virtual site. In this way, the total force and the total torque are
conserved \ `104 <#ref-Berendsen84b>`__.

The computation of the virial (eqn. [eqn:Xi]) is non-trivial when
virtual sites are used. Since the virial involves a summation over all
the atoms (rather than virtual sites), the forces must be redistributed
from the virtual sites to the atoms (using  eqn. [eqn:fvsite]) *before*
computation of the virial. In some special cases where the forces on the
atoms can be written as a linear combination of the forces on the
virtual sites (types 2 and 3 below) there is no difference between
computing the virial before and after the redistribution of forces.
However, in the general case redistribution should be done first.

.. figure:: plots/dummies
   :alt: The six different types of virtual site construction in . The
   constructing atoms are shown as black circles, the virtual sites in
   gray.
   :width: 15.00000cm

   The six different types of virtual site construction in . The
   constructing atoms are shown as black circles, the virtual sites in
   gray.

There are six ways to construct virtual sites from surrounding atoms in
GROMACS, which we classify by the number of constructing atoms. **Note**
that all site types mentioned can be constructed from types 3fd
(normalized, in-plane) and 3out (non-normalized, out of plane). However,
the amount of computation involved increases sharply along this list, so
we strongly recommended using the first adequate virtual site type that
will be sufficient for a certain purpose. Fig. [fig:vsites] depicts 6 of
the available virtual site constructions. The conceptually simplest
construction types are linear combinations:

.. math:: {{\mbox{\boldmath ${r}$}}_s}= \sum_{i=1}^N w_i \, {\mbox{\boldmath ${r}$}}_i

 The force is then redistributed using the same weights:

.. math:: {{\mbox{\boldmath ${F}$}}_i'}= w_i \, {{\mbox{\boldmath ${F}$}}_{s}}

The types of virtual sites supported in GROMACS are given in the list
below. Constructing atoms in virtual sites can be virtual sites
themselves, but only if they are higher in the list, i.e. virtual sites
can be constructed from “particles” that are simpler virtual sites.

-  [subsec:vsite2]As a linear combination of two atoms
   (Fig. [fig:vsites] 2):

   .. math:: w_i = 1 - a ~,~~ w_j = a

    In this case the virtual site is on the line through atoms :math:`i`
   and :math:`j`.

-  [subsec:vsite3]As a linear combination of three atoms
   (Fig. [fig:vsites] 3):

   .. math:: w_i = 1 - a - b ~,~~ w_j = a ~,~~ w_k = b

    In this case the virtual site is in the plane of the other three
   particles.

-  [subsec:vsite3fd]In the plane of three atoms, with a fixed distance
   (Fig. [fig:vsites] 3fd):

   .. math::

      {{\mbox{\boldmath ${r}$}}_s}~=~ {\mbox{\boldmath ${r}$}}_i + b \frac{  {{\mbox{\boldmath ${r}$}}_{ij}}+ a {{\mbox{\boldmath ${r}$}}_{jk}}}
                                          {| {{\mbox{\boldmath ${r}$}}_{ij}}+ a {{\mbox{\boldmath ${r}$}}_{jk}}|}

    In this case the virtual site is in the plane of the other three
   particles at a distance of :math:`|b|` from :math:`i`. The force on
   particles :math:`i`, :math:`j` and :math:`k` due to the force on the
   virtual site can be computed as:

   .. math::

      \begin{array}{lcr}
              {{\mbox{\boldmath ${F}$}}_i'}&=& \displaystyle {{\mbox{\boldmath ${F}$}}_{s}}- \gamma ( {{\mbox{\boldmath ${F}$}}_{s}}- {\mbox{\boldmath ${p}$}} ) \\[1ex]
              {{\mbox{\boldmath ${F}$}}_j'}&=& \displaystyle (1-a)\gamma ({{\mbox{\boldmath ${F}$}}_{s}}- {\mbox{\boldmath ${p}$}})      \\[1ex]
              {{\mbox{\boldmath ${F}$}}_k'}&=& \displaystyle a \gamma ({{\mbox{\boldmath ${F}$}}_{s}}- {\mbox{\boldmath ${p}$}})         \\
              \end{array}
              ~\mbox{~ where~ }~
              \begin{array}{c}
      \displaystyle \gamma = \frac{b}{| {{\mbox{\boldmath ${r}$}}_{ij}}+ a {{\mbox{\boldmath ${r}$}}_{jk}}|} \\[2ex]
      \displaystyle {\mbox{\boldmath ${p}$}} = \frac{ {{\mbox{\boldmath ${r}$}}_{is}}\cdot {{\mbox{\boldmath ${F}$}}_{s}}}
                            { {{\mbox{\boldmath ${r}$}}_{is}}\cdot {{\mbox{\boldmath ${r}$}}_{is}}} {{\mbox{\boldmath ${r}$}}_{is}}\end{array}

-  [subsec:vsite3fad]In the plane of three atoms, with a fixed angle and
   distance (Fig. [fig:vsites] 3fad):

   .. math::

      \label{eqn:vsite2fad-F}
               {{\mbox{\boldmath ${r}$}}_s}~=~ {\mbox{\boldmath ${r}$}}_i +
                          d \cos \theta \frac{{{\mbox{\boldmath ${r}$}}_{ij}}}{|{{\mbox{\boldmath ${r}$}}_{ij}}|} +
                          d \sin \theta \frac{{\mbox{\boldmath ${r}$}}_\perp}{|{\mbox{\boldmath ${r}$}}_\perp|}
              ~\mbox{~ where~ }~
              {\mbox{\boldmath ${r}$}}_\perp ~=~ {{\mbox{\boldmath ${r}$}}_{jk}}- 
                              \frac{ {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${r}$}}_{jk}}}
                                   { {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${r}$}}_{ij}}}
                               {{\mbox{\boldmath ${r}$}}_{ij}}

    In this case the virtual site is in the plane of the other three
   particles at a distance of :math:`|d|` from :math:`i` at an angle of
   :math:`\alpha` with :math:`{{\mbox{\boldmath ${r}$}}_{ij}}`. Atom
   :math:`k` defines the plane and the direction of the angle. **Note**
   that in this case :math:`b` and :math:`\alpha` must be specified,
   instead of :math:`a` and :math:`b` (see also sec. [sec:vsitetop]).
   The force on particles :math:`i`, :math:`j` and :math:`k` due to the
   force on the virtual site can be computed as (with
   :math:`{\mbox{\boldmath ${r}$}}_\perp` as defined in
   eqn. [eqn:vsite2fad-F]):

   .. math::

      \begin{array}{c}
              \begin{array}{lclllll}
              {{\mbox{\boldmath ${F}$}}_i'}&=& {{\mbox{\boldmath ${F}$}}_{s}}&-& 
                      {\displaystyle\frac}{d \cos \theta}{|{{\mbox{\boldmath ${r}$}}_{ij}}|} {\mbox{\boldmath ${F}$}}_1 &+&
                      {\displaystyle\frac}{d \sin \theta}{|{\mbox{\boldmath ${r}$}}_\perp|} \left( 
                      {\displaystyle\frac}{ {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${r}$}}_{jk}}}
                           { {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${r}$}}_{ij}}} {\mbox{\boldmath ${F}$}}_2     +
                      {\mbox{\boldmath ${F}$}}_3 \right)                                \\[3ex]
              {{\mbox{\boldmath ${F}$}}_j'}&=& &&
                      {\displaystyle\frac}{d \cos \theta}{|{{\mbox{\boldmath ${r}$}}_{ij}}|} {\mbox{\boldmath ${F}$}}_1 &-&
                      {\displaystyle\frac}{d \sin \theta}{|{\mbox{\boldmath ${r}$}}_\perp|} \left(
                       {\mbox{\boldmath ${F}$}}_2 + 
                       {\displaystyle\frac}{ {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${r}$}}_{jk}}}
                              { {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${r}$}}_{ij}}} {\mbox{\boldmath ${F}$}}_2 +
                      {\mbox{\boldmath ${F}$}}_3 \right)                                \\[3ex]
              {{\mbox{\boldmath ${F}$}}_k'}&=& && &&
                      {\displaystyle\frac}{d \sin \theta}{|{\mbox{\boldmath ${r}$}}_\perp|} {\mbox{\boldmath ${F}$}}_2  \\[3ex]
              \end{array}                                             \\[5ex]
              \mbox{where ~}
              {\mbox{\boldmath ${F}$}}_1 = {{\mbox{\boldmath ${F}$}}_{s}}-
                        {\displaystyle\frac}{ {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${F}$}}_{s}}}
                              { {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${r}$}}_{ij}}} {{\mbox{\boldmath ${r}$}}_{ij}}\mbox{\,, ~}
              {\mbox{\boldmath ${F}$}}_2 = {\mbox{\boldmath ${F}$}}_1 -
                        {\displaystyle\frac}{ {\mbox{\boldmath ${r}$}}_\perp \cdot {{\mbox{\boldmath ${F}$}}_{s}}}
                              { {\mbox{\boldmath ${r}$}}_\perp \cdot {\mbox{\boldmath ${r}$}}_\perp } {\mbox{\boldmath ${r}$}}_\perp
              \mbox{~and ~}
              {\mbox{\boldmath ${F}$}}_3 = {\displaystyle\frac}{ {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${F}$}}_{s}}}
                               { {{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${r}$}}_{ij}}} {\mbox{\boldmath ${r}$}}_\perp
      \end{array}

-  [subsec:vsite3out]As a non-linear combination of three atoms, out of
   plane (Fig. [fig:vsites] 3out):

   .. math::

      {{\mbox{\boldmath ${r}$}}_s}~=~ {\mbox{\boldmath ${r}$}}_i + a {{\mbox{\boldmath ${r}$}}_{ij}}+ b {{\mbox{\boldmath ${r}$}}_{ik}}+
                      c ({{\mbox{\boldmath ${r}$}}_{ij}}\times {{\mbox{\boldmath ${r}$}}_{ik}})

    This enables the construction of virtual sites out of the plane of
   the other atoms. The force on particles :math:`i,j` and :math:`k` due
   to the force on the virtual site can be computed as:

   .. math::

      \begin{array}{lcl}
      \vspace{4mm}
      {{\mbox{\boldmath ${F}$}}_j'}&=& \left[\begin{array}{ccc}
       a              &  -c\,z_{ik}   & c\,y_{ik}     \\[0.5ex]
       c\,z_{ik}      &   a           & -c\,x_{ik}    \\[0.5ex]
      -c\,y_{ik}      &   c\,x_{ik}   & a
      \end{array}\right]{{\mbox{\boldmath ${F}$}}_{s}}\\
      \vspace{4mm}
      {{\mbox{\boldmath ${F}$}}_k'}&=& \left[\begin{array}{ccc}
       b              &   c\,z_{ij}   & -c\,y_{ij}    \\[0.5ex]
      -c\,z_{ij}      &   b           & c\,x_{ij}     \\[0.5ex]
       c\,y_{ij}      &  -c\,x_{ij}   & b
      \end{array}\right]{{\mbox{\boldmath ${F}$}}_{s}}\\
      {{\mbox{\boldmath ${F}$}}_i'}&=& {{\mbox{\boldmath ${F}$}}_{s}}- {{\mbox{\boldmath ${F}$}}_j'}- {{\mbox{\boldmath ${F}$}}_k'}\end{array}

-  [subsec:vsite4fdn]From four atoms, with a fixed distance, see
   separate Fig. [fig:vsite-4fdn]. This construction is a bit complex,
   in particular since the previous type (4fd) could be unstable which
   forced us to introduce a more elaborate construction:

   .. figure:: plots/vsite-4fdn
      :alt: The new 4fdn virtual site construction, which is stable even
      when all constructing atoms are in the same plane.
      :width: 5.00000cm

      The new 4fdn virtual site construction, which is stable even when
      all constructing atoms are in the same plane.

   .. math::

      \begin{aligned}
      \mathbf{r}_{ja} &=& a\, \mathbf{r}_{ik} - \mathbf{r}_{ij} = a\, (\mathbf{x}_k - \mathbf{x}_i) - (\mathbf{x}_j - \mathbf{x}_i) \nonumber \\
      \mathbf{r}_{jb} &=& b\, \mathbf{r}_{il} - \mathbf{r}_{ij} = b\, (\mathbf{x}_l - \mathbf{x}_i) - (\mathbf{x}_j - \mathbf{x}_i) \nonumber \\
      \mathbf{r}_m &=& \mathbf{r}_{ja} \times \mathbf{r}_{jb} \nonumber \\
      \mathbf{x}_s &=& \mathbf{x}_i + c \frac{\mathbf{r}_m}{|\mathbf{r}_m|}
      \label{eq:vsite}\end{aligned}

   In this case the virtual site is at a distance of :math:`|c|` from
   :math:`i`, while :math:`a` and :math:`b` are parameters. **Note**
   that the vectors :math:`\mathbf{r}_{ik}` and :math:`\mathbf{r}_{ij}`
   are not normalized to save floating-point operations. The force on
   particles :math:`i`, :math:`j`, :math:`k` and :math:`l` due to the
   force on the virtual site are computed through chain rule derivatives
   of the construction expression. This is exact and conserves energy,
   but it does lead to relatively lengthy expressions that we do not
   include here (over 200 floating-point operations). The interested
   reader can look at the source code in ``vsite.c``. Fortunately, this
   vsite type is normally only used for chiral centers such as
   :math:`C_{\alpha}` atoms in proteins.

   The new 4fdn construct is identified with a ‘type’ value of 2 in the
   topology. The earlier 4fd type is still supported internally (‘type’
   value 1), but it should not be used for new simulations. All current
   GROMACS tools will automatically generate type 4fdn instead.

-  [subsec:vsiteN] A linear combination of :math:`N` atoms with relative
   weights :math:`a_i`. The weight for atom :math:`i` is:

   .. math:: w_i = a_i \left(\sum_{j=1}^N a_j \right)^{-1}

    There are three options for setting the weights:

   -  center of geometry: equal weights

   -  center of mass: :math:`a_i` is the mass of atom :math:`i`; when in
      free-energy simulations the mass of the atom is changed, only the
      mass of the A-state is used for the weight

   -  center of weights: :math:`a_i` is defined by the user

Long Range Electrostatics
-------------------------

Ewald summation
~~~~~~~~~~~~~~~

The total electrostatic energy of :math:`N` particles and their periodic
images is given by

.. math::

   V=\frac{f}{2}\sum_{n_x}\sum_{n_y}
   \sum_{n_{z}*} \sum_{i}^{N} \sum_{j}^{N}
   \frac{q_i q_j}{{\bf r}_{ij,{\bf n}}}.
   \label{eqn:totalcoulomb}

 :math:`(n_x,n_y,n_z)={\bf n}` is the box index vector, and the star
indicates that terms with :math:`i=j` should be omitted when
:math:`(n_x,n_y,n_z)=(0,0,0)`. The distance :math:`{\bf r}_{ij,{\bf n}}`
is the real distance between the charges and not the minimum-image. This
sum is conditionally convergent, but very slow.

Ewald summation was first introduced as a method to calculate long-range
interactions of the periodic images in
crystals \ `105 <#ref-Ewald21>`__. The idea is to convert the single
slowly-converging sum eqn. [eqn:totalcoulomb] into two
quickly-converging terms and a constant term:

.. math::

   \begin{aligned}
   V &=& V_{\mathrm{dir}} + V_{\mathrm{rec}} + V_{0} \\[0.5ex]
   V_{\mathrm{dir}} &=& \frac{f}{2} \sum_{i,j}^{N}
   \sum_{n_x}\sum_{n_y}
   \sum_{n_{z}*} q_i q_j \frac{\mbox{erfc}(\beta {r}_{ij,{\bf n}} )}{{r}_{ij,{\bf n}}} \\[0.5ex]
   V_{\mathrm{rec}} &=& \frac{f}{2 \pi V} \sum_{i,j}^{N} q_i q_j
   \sum_{m_x}\sum_{m_y}
   \sum_{m_{z}*} \frac{\exp{\left( -(\pi {\bf m}/\beta)^2 + 2 \pi i
         {\bf m} \cdot ({\bf r}_i - {\bf r}_j)\right)}}{{\bf m}^2} \\[0.5ex]
   V_{0} &=& -\frac{f \beta}{\sqrt{\pi}}\sum_{i}^{N} q_i^2,\end{aligned}

 where :math:`\beta` is a parameter that determines the relative weight
of the direct and reciprocal sums and :math:`{\bf m}=(m_x,m_y,m_z)`. In
this way we can use a short cut-off (of the order of :math:`1` nm) in
the direct space sum and a short cut-off in the reciprocal space sum
(*e.g.* 10 wave vectors in each direction). Unfortunately, the
computational cost of the reciprocal part of the sum increases as
:math:`N^2` (or :math:`N^{3/2}` with a slightly better algorithm) and it
is therefore not realistic for use in large systems.

Using Ewald
^^^^^^^^^^^

Don’t use Ewald unless you are absolutely sure this is what you want -
for almost all cases the PME method below will perform much better. If
you still want to employ classical Ewald summation enter this in your
.mdp file, if the side of your box is about :math:`3` nm:

::

    coulombtype     = Ewald
    rvdw            = 0.9
    rlist           = 0.9
    rcoulomb        = 0.9
    fourierspacing  = 0.6
    ewald-rtol      = 1e-5

The ratio of the box dimensions and the fourierspacing parameter
determines the highest magnitude of wave vectors :math:`m_x,m_y,m_z` to
use in each direction. With a 3-nm cubic box this example would use
:math:`11` wave vectors (from :math:`-5` to :math:`5`) in each
direction. The ewald-rtol parameter is the relative strength of the
electrostatic interaction at the cut-off. Decreasing this gives you a
more accurate direct sum, but a less accurate reciprocal sum.

PME
~~~

Particle-mesh Ewald is a method proposed by Tom
Darden \ `14 <#ref-Darden93>`__ to improve the performance of the
reciprocal sum. Instead of directly summing wave vectors, the charges
are assigned to a grid using interpolation. The implementation in
GROMACS uses cardinal B-spline interpolation \ `15 <#ref-Essmann95>`__,
which is referred to as smooth PME (SPME). The grid is then Fourier
transformed with a 3D FFT algorithm and the reciprocal energy term
obtained by a single sum over the grid in k-space.

The potential at the grid points is calculated by inverse
transformation, and by using the interpolation factors we get the forces
on each atom.

The PME algorithm scales as :math:`N \log(N)`, and is substantially
faster than ordinary Ewald summation on medium to large systems. On very
small systems it might still be better to use Ewald to avoid the
overhead in setting up grids and transforms. For the parallelization of
PME see the section on MPMD PME ([subsec:mpmd\_pme]).

With the Verlet cut-off scheme, the PME direct space potential is
shifted by a constant such that the potential is zero at the cut-off.
This shift is small and since the net system charge is close to zero,
the total shift is very small, unlike in the case of the Lennard-Jones
potential where all shifts add up. We apply the shift anyhow, such that
the potential is the exact integral of the force.

Using PME
^^^^^^^^^

As an example for using Particle-mesh Ewald summation in GROMACS,
specify the following lines in your .mdp file:

::

    coulombtype     = PME
    rvdw            = 0.9
    rlist           = 0.9
    rcoulomb        = 0.9
    fourierspacing  = 0.12
    pme-order       = 4
    ewald-rtol      = 1e-5

In this case the fourierspacing parameter determines the maximum spacing
for the FFT grid (i.e. minimum number of grid points), and pme-order
controls the interpolation order. Using fourth-order (cubic)
interpolation and this spacing should give electrostatic energies
accurate to about :math:`5\cdot10^{-3}`. Since the Lennard-Jones
energies are not this accurate it might even be possible to increase
this spacing slightly.

Pressure scaling works with PME, but be aware of the fact that
anisotropic scaling can introduce artificial ordering in some systems.

P3M-AD
~~~~~~

The Particle-Particle Particle-Mesh methods of Hockney & Eastwood can
also be applied in GROMACS for the treatment of long range electrostatic
interactions \ `106 <#ref-Hockney81>`__. Although the P3M method was the
first efficient long-range electrostatics method for molecular
simulation, the smooth PME (SPME) method has largely replaced P3M as the
method of choice in atomistic simulations. One performance disadvantage
of the original P3M method was that it required 3 3D-FFT back transforms
to obtain the forces on the particles. But this is not required for P3M
and the forces can be derived through analytical differentiation of the
potential, as done in PME. The resulting method is termed P3M-AD. The
only remaining difference between P3M-AD and PME is the optimization of
the lattice Green influence function for error minimization that P3M
uses. However, in 2012 it has been shown that the SPME influence
function can be modified to obtain P3M \ `107 <#ref-Ballenegger2012>`__.
This means that the advantage of error minimization in P3M-AD can be
used at the same computational cost and with the same code as PME, just
by adding a few lines to modify the influence function. However, at
optimal parameter setting the effect of error minimization in P3M-AD is
less than 10%. P3M-AD does show large accuracy gains with interlaced
(also known as staggered) grids, but that is not supported in GROMACS
(yet).

P3M is used in GROMACS with exactly the same options as used with PME by
selecting the electrostatics type:

::

    coulombtype     = P3M-AD

Optimizing Fourier transforms and PME calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is recommended to optimize the parameters for calculation of
electrostatic interaction such as PME grid dimensions and cut-off radii.
This is particularly relevant to do before launching long production
runs.

gmx mdrun will automatically do a lot of PME optimization, and GROMACS
also includes a special tool, gmx tune\_pme, which automates the process
of selecting the optimal number of PME-only ranks.

Long Range Van der Waals interactions
-------------------------------------

Dispersion correction
~~~~~~~~~~~~~~~~~~~~~

In this section, we derive long-range corrections due to the use of a
cut-off for Lennard-Jones or Buckingham interactions. We assume that the
cut-off is so long that the repulsion term can safely be neglected, and
therefore only the dispersion term is taken into account. Due to the
nature of the dispersion interaction (we are truncating a potential
proportional to :math:`-r^{-6}`), energy and pressure corrections are
both negative. While the energy correction is usually small, it may be
important for free energy calculations where differences between two
different Hamiltonians are considered. In contrast, the pressure
correction is very large and can not be neglected under any
circumstances where a correct pressure is required, especially for any
NPT simulations. Although it is, in principle, possible to parameterize
a force field such that the pressure is close to the desired
experimental value without correction, such a method makes the
parameterization dependent on the cut-off and is therefore undesirable.

Energy
^^^^^^

The long-range contribution of the dispersion interaction to the virial
can be derived analytically, if we assume a homogeneous system beyond
the cut-off distance :math:`r_c`. The dispersion energy between two
particles is written as:

.. math:: V({r_{ij}}) ~=~- C_6\,{r_{ij}}^{-6}

 and the corresponding force is:

.. math:: {{\mbox{\boldmath ${F}$}}_{ij}}~=~- 6\,C_6\,{r_{ij}}^{-8}{{\mbox{\boldmath ${r}$}}_{ij}}

 In a periodic system it is not easy to calculate the full potentials,
so usually a cut-off is applied, which can be abrupt or smooth. We will
call the potential and force with cut-off :math:`V_c` and
:math:`{\mbox{\boldmath ${F}$}}_c`. The long-range contribution to the
dispersion energy in a system with :math:`N` particles and particle
density :math:`\rho` = :math:`N/V` is:

.. math::

   \label{eqn:enercorr}
   V_{lr} ~=~ {\frac{1}{2}}N \rho\int_0^{\infty}   4\pi r^2 g(r) \left( V(r) -V_c(r) \right) {{{\rm d}r}}

 We will integrate this for the shift function, which is the most
general form of van der Waals interaction available in GROMACS. The
shift function has a constant difference :math:`S` from 0 to :math:`r_1`
and is 0 beyond the cut-off distance :math:`r_c`. We can integrate
eqn. [eqn:enercorr], assuming that the density in the sphere within
:math:`r_1` is equal to the global density and the radial distribution
function :math:`g(r)` is 1 beyond :math:`r_1`:

.. math::

   \begin{aligned}
   \nonumber
   V_{lr}  &=& {\frac{1}{2}}N \left(
     \rho\int_0^{r_1}  4\pi r^2 g(r) \, C_6 \,S\,{{{\rm d}r}}
   + \rho\int_{r_1}^{r_c}  4\pi r^2 \left( V(r) -V_c(r) \right) {{{\rm d}r}}
   + \rho\int_{r_c}^{\infty}  4\pi r^2 V(r) \, {{{\rm d}r}}
   \right) \\
   & = & {\frac{1}{2}}N \left(\left(\frac{4}{3}\pi \rho r_1^{3} - 1\right) C_6 \,S
   + \rho\int_{r_1}^{r_c} 4\pi r^2 \left( V(r) -V_c(r) \right) {{{\rm d}r}}
   -\frac{4}{3} \pi N \rho\, C_6\,r_c^{-3}
   \right)\end{aligned}

 where the term :math:`-1` corrects for the self-interaction. For a
plain cut-off we only need to assume that :math:`g(r)` is 1 beyond
:math:`r_c` and the correction reduces to \ `108 <#ref-Allen87>`__:

.. math::

   \begin{aligned}
   V_{lr} & = & -\frac{2}{3} \pi N \rho\, C_6\,r_c^{-3}\end{aligned}

 If we consider, for example, a box of pure water, simulated with a
cut-off of 0.9 nm and a density of 1 g cm\ :math:`^{-3}` this correction
is :math:`-0.75` kJ mol\ :math:`^{-1}` per molecule.

For a homogeneous mixture we need to define an *average dispersion
constant*:

.. math::

   \label{eqn:avcsix}
   {\left< C_6 \right>}= \frac{2}{N(N-1)}\sum_i^N\sum_{j>i}^N C_6(i,j)\\

 In GROMACS, excluded pairs of atoms do not contribute to the average.

In the case of inhomogeneous simulation systems, *e.g.* a system with a
lipid interface, the energy correction can be applied if
:math:`{\left< C_6 \right>}` for both components is comparable.

Virial and pressure
^^^^^^^^^^^^^^^^^^^

The scalar virial of the system due to the dispersion interaction
between two particles :math:`i` and :math:`j` is given by:

.. math:: \Xi~=~-{\frac{1}{2}}{{\mbox{\boldmath ${r}$}}_{ij}}\cdot {{\mbox{\boldmath ${F}$}}_{ij}}~=~ 3\,C_6\,{r_{ij}}^{-6}

 The pressure is given by:

.. math:: P~=~\frac{2}{3\,V}\left(E_{kin} - \Xi\right)

 The long-range correction to the virial is given by:

.. math:: \Xi_{lr} ~=~ {\frac{1}{2}}N \rho \int_0^{\infty} 4\pi r^2 g(r) (\Xi -\Xi_c) \,{{\rm d}r}

 We can again integrate the long-range contribution to the virial
assuming :math:`g(r)` is 1 beyond :math:`r_1`:

.. math::

   \begin{aligned}
   \Xi_{lr}&=&	{\frac{1}{2}}N \rho \left(
       \int_{r_1}^{r_c}  4 \pi r^2 (\Xi -\Xi_c)  \,{{\rm d}r}+ \int_{r_c}^{\infty} 4 \pi r^2 3\,C_6\,{r_{ij}}^{-6}\,  {{\rm d}r}\right)	\nonumber\\
           &=&     {\frac{1}{2}}N \rho \left(
       \int_{r_1}^{r_c} 4 \pi r^2 (\Xi -\Xi_c) \, {{\rm d}r}+ 4 \pi C_6 \, r_c^{-3} \right)\end{aligned}

 For a plain cut-off the correction to the pressure
is \ `108 <#ref-Allen87>`__:

.. math:: P_{lr}~=~-\frac{4}{3} \pi C_6\, \rho^2 r_c^{-3}

 Using the same example of a water box, the correction to the virial is
0.75 kJ mol\ :math:`^{-1}` per molecule, the corresponding correction to
the pressure for SPC water is approximately :math:`-280` bar.

For homogeneous mixtures, we can again use the average dispersion
constant :math:`{\left< C_6 \right>}` (eqn. [eqn:avcsix]):

.. math::

   P_{lr}~=~-\frac{4}{3} \pi {\left< C_6 \right>}\rho^2 r_c^{-3}
   \label{eqn:pcorr}

 For inhomogeneous systems, eqn. [eqn:pcorr] can be applied under the
same restriction as holds for the energy (see sec. [sec:ecorr]).

Lennard-Jones PME
~~~~~~~~~~~~~~~~~

In order to treat systems, using Lennard-Jones potentials, that are
non-homogeneous outside of the cut-off distance, we can instead use the
Particle-mesh Ewald method as discussed for electrostatics above. In
this case the modified Ewald equations become

.. math::

   \begin{aligned}
   V &=& V_{\mathrm{dir}} + V_{\mathrm{rec}} + V_{0} \\[0.5ex]
   V_{\mathrm{dir}} &=& -\frac{1}{2} \sum_{i,j}^{N}
   \sum_{n_x}\sum_{n_y}
   \sum_{n_{z}*} \frac{C^{ij}_6 g(\beta {r}_{ij,{\bf n}})}{{r_{ij,{\bf n}}}^6}
   \label{eqn:ljpmerealspace}\\[0.5ex]
   V_{\mathrm{rec}} &=& \frac{{\pi}^{\frac{3}{2}} \beta^{3}}{2V} \sum_{m_x}\sum_{m_y}\sum_{m_{z}*}
   f(\pi |{\mathbf m}|/\beta) \times \sum_{i,j}^{N} C^{ij}_6 {\mathrm{exp}}\left[-2\pi i {\bf m}\cdot({\bf r_i}-{\bf r_j})\right] \\[0.5ex]
   V_{0} &=& -\frac{\beta^{6}}{12}\sum_{i}^{N} C^{ii}_6\end{aligned}

where :math:`{\bf m}=(m_x,m_y,m_z)`, :math:`\beta` is the parameter
determining the weight between direct and reciprocal space, and
:math:`{C^{ij}_6}` is the combined dispersion parameter for particle
:math:`i` and :math:`j`. The star indicates that terms with
:math:`i = j` should be omitted when :math:`((n_x,n_y,n_z)=(0,0,0))`,
and :math:`{\bf r}_{ij,{\bf n}}` is the real distance between the
particles. Following the derivation by
Essmann \ `15 <#ref-Essmann95>`__, the functions :math:`f` and :math:`g`
introduced above are defined as

.. math::

   \begin{aligned}
   f(x)&=&1/3\left[(1-2x^2){\mathrm{exp}}(-x^2) + 2{x^3}\sqrt{\pi}\,{\mathrm{erfc}}(x) \right] \\
   g(x)&=&{\mathrm{exp}}(-x^2)(1+x^2+\frac{x^4}{2}).\end{aligned}

The above methodology works fine as long as the dispersion parameters
can be combined geometrically (eqn. [eqn:comb]) in the same way as the
charges for electrostatics

.. math:: C^{ij}_{6,\mathrm{geom}} = \left(C^{ii}_6 \, C^{jj}_6\right)^{1/2}

 For Lorentz-Berthelot combination rules (eqn. [eqn:lorentzberthelot]),
the reciprocal part of this sum has to be calculated seven times due to
the splitting of the dispersion parameter according to

.. math:: C^{ij}_{6,\mathrm{L-B}} = (\sigma_i+\sigma_j)^6=\sum_{n=0}^{6} P_{n}\sigma_{i}^{n}\sigma_{j}^{(6-n)},

 for :math:`P_{n}` the Pascal triangle coefficients. This introduces a
non-negligible cost to the reciprocal part, requiring seven separate
FFTs, and therefore this has been the limiting factor in previous
attempts to implement LJ-PME. A solution to this problem is to use
geometrical combination rules in order to calculate an approximate
interaction parameter for the reciprocal part of the potential, yielding
a total interaction of

.. math::

   \begin{aligned}
   V(r<r_c) & = & \underbrace{C^{\mathrm{dir}}_6 g(\beta r) r^{-6}}_{\mathrm{Direct \  space}} + \underbrace{C^\mathrm{recip}_{6,\mathrm{geom}} [1 - g(\beta r)] r^{-6}}_{\mathrm{Reciprocal \  space}} \nonumber \\
   &=& C^\mathrm{recip}_{6,\mathrm{geom}}r^{-6} + \left(C^{\mathrm{dir}}_6-C^\mathrm{recip}_{6,\mathrm{geom}}\right)g(\beta r)r^{-6} \\
   V(r>r_c) & = & \underbrace{C^\mathrm{recip}_{6,\mathrm{geom}} [1 - g(\beta r)] r^{-6}}_{\mathrm{Reciprocal \  space}}.\end{aligned}

 This will preserve a well-defined Hamiltonian and significantly
increase the performance of the simulations. The approximation does
introduce some errors, but since the difference is located in the
interactions calculated in reciprocal space, the effect will be very
small compared to the total interaction energy. In a simulation of a
lipid bilayer, using a cut-off of 1.0 nm, the relative error in total
dispersion energy was below 0.5%. A more thorough discussion of this can
be found in `109 <#ref-Wennberg13>`__.

In GROMACS we now perform the proper calculation of this interaction by
subtracting, from the direct-space interactions, the contribution made
by the approximate potential that is used in the reciprocal part

.. math::

   V_\mathrm{dir} = C^{\mathrm{dir}}_6 r^{-6} - C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}.
   \label{eqn:ljpmedirectspace}

 This potential will reduce to the expression in
eqn. [eqn:ljpmerealspace] when
:math:`C^{\mathrm{dir}}_6 = C^\mathrm{recip}_6`, and the total
interaction is given by

.. math::

   \begin{aligned}
   \nonumber V(r<r_c) &=& \underbrace{C^{\mathrm{dir}}_6 r^{-6} - C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}}_{\mathrm{Direct \  space}} + \underbrace{C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}}_{\mathrm{Reciprocal \  space}} \\ 
   &=&C^{\mathrm{dir}}_6 r^{-6}
   \label {eqn:ljpmecorr2} \\
   V(r>r_c) &=& C^\mathrm{recip}_6 [1 - g(\beta r)] r^{-6}.\end{aligned}

 For the case when :math:`C^{\mathrm{dir}}_6 \neq C^\mathrm{recip}_6`
this will retain an unmodified LJ force up to the cut-off, and the error
is an order of magnitude smaller than in simulations where the
direct-space interactions do not account for the approximation used in
reciprocal space. When using a VdW interaction modifier of
potential-shift, the constant

.. math:: \left(-C^{\mathrm{dir}}_6 + C^\mathrm{recip}_6 [1 - g(\beta r_c)]\right) r_c^{-6}

 is added to eqn. [eqn:ljpmecorr2] in order to ensure that the potential
is continuous at the cutoff. Note that, in the same way as
eqn. [eqn:ljpmedirectspace], this degenerates into the expected
:math:`-C_6g(\beta r_c)r^{-6}_c` when :math:`C^{\mathrm{dir}}_6 =
C^\mathrm{recip}_6`. In addition to this, a long-range dispersion
correction can be applied to correct for the approximation using a
combination rule in reciprocal space. This correction assumes, as for
the cut-off LJ potential, a uniform particle distribution. But since the
error of the combination rule approximation is very small this
long-range correction is not necessary in most cases. Also note that
this homogenous correction does not correct the surface tension, which
is an inhomogeneous property.

Using LJ-PME
^^^^^^^^^^^^

As an example for using Particle-mesh Ewald summation for Lennard-Jones
interactions in GROMACS, specify the following lines in your .mdp file:

::

    vdwtype          = PME
    rvdw             = 0.9
    vdw-modifier     = Potential-Shift
    rlist            = 0.9
    rcoulomb         = 0.9
    fourierspacing   = 0.12
    pme-order        = 4
    ewald-rtol-lj    = 0.001
    lj-pme-comb-rule = geometric

The same Fourier grid and interpolation order are used if both LJ-PME
and electrostatic PME are active, so the settings for fourierspacing and
pme-order are common to both. ewald-rtol-lj controls the splitting
between direct and reciprocal space in the same way as ewald-rtol. In
addition to this, the combination rule to be used in reciprocal space is
determined by lj-pme-comb-rule. If the current force field uses
Lorentz-Berthelot combination rules, it is possible to set
lj-pme-comb-rule = geometric in order to gain a significant increase in
performance for a small loss in accuracy. The details of this
approximation can be found in the section above.

Note that the use of a complete long-range dispersion correction means
that as with Coulomb PME, rvdw is now a free parameter in the method,
rather than being necessarily restricted by the force-field
parameterization scheme. Thus it is now possible to optimize the cutoff,
spacing, order and tolerance terms for accuracy and best performance.

Naturally, the use of LJ-PME rather than LJ cut-off adds computation and
communication done for the reciprocal-space part, so for best
performance in balancing the load of parallel simulations using PME-only
ranks, more such ranks should be used. It may be possible to improve
upon the automatic load-balancing used by mdrun.

Force field
-----------

A force field is built up from two distinct components:

-  The set of equations (called the *potential functions*) used to
   generate the potential energies and their derivatives, the forces.
   These are described in detail in the previous chapter.

-  The parameters used in this set of equations. These are not given in
   this manual, but in the data files corresponding to your GROMACS
   distribution.

Within one set of equations various sets of parameters can be used. Care
must be taken that the combination of equations and parameters form a
consistent set. It is in general dangerous to make *ad hoc* changes in a
subset of parameters, because the various contributions to the total
force are usually interdependent. This means in principle that every
change should be documented, verified by comparison to experimental data
and published in a peer-reviewed journal before it can be used.

GROMACS @GMX\_VERSION\_STRING@ includes several force fields, and
additional ones are available on the website. If you do not know which
one to select we recommend GROMOS-96 for united-atom setups and
OPLS-AA/L for all-atom parameters. That said, we describe the available
options in some detail.

All-hydrogen force field
^^^^^^^^^^^^^^^^^^^^^^^^

The GROMOS-87-based all-hydrogen force field is almost identical to the
normal GROMOS-87 force field, since the extra hydrogens have no
Lennard-Jones interaction and zero charge. The only differences are in
the bond angle and improper dihedral angle terms. This force field is
only useful when you need the exact hydrogen positions, for instance for
distance restraints derived from NMR measurements. When citing this
force field please read the previous paragraph.

GROMOS-96
~~~~~~~~~

GROMACS supports the GROMOS-96 force fields \ `77 <#ref-gromos96>`__.
All parameters for the 43A1, 43A2 (development, improved alkane
dihedrals), 45A3, 53A5, and 53A6 parameter sets are included. All
standard building blocks are included and topologies can be built
automatically by pdb2gmx.

The GROMOS-96 force field is a further development of the GROMOS-87
force field. It has improvements over the GROMOS-87 force field for
proteins and small molecules. **Note** that the sugar parameters present
in 53A6 do correspond to those published in
2004\ `110 <#ref-Oostenbrink2004>`__, which are different from those
present in 45A4, which is not included in GROMACS at this time. The 45A4
parameter set corresponds to a later revision of these parameters. The
GROMOS-96 force field is not, however, recommended for use with long
alkanes and lipids. The GROMOS-96 force field differs from the GROMOS-87
force field in a few respects:

-  the force field parameters

-  the parameters for the bonded interactions are not linked to atom
   types

-  a fourth power bond stretching potential ([subsec:G96bond])

-  an angle potential based on the cosine of the angle
   ([subsec:G96angle])

There are two differences in implementation between GROMACS and
GROMOS-96 which can lead to slightly different results when simulating
the same system with both packages:

-  in GROMOS-96 neighbor searching for solvents is performed on the
   first atom of the solvent molecule. This is not implemented in
   GROMACS, but the difference with searching by centers of charge
   groups is very small

-  the virial in GROMOS-96 is molecule-based. This is not implemented in
   GROMACS, which uses atomic virials

The GROMOS-96 force field was parameterized with a Lennard-Jones cut-off
of 1.4 nm, so be sure to use a Lennard-Jones cut-off (rvdw) of at least
1.4. A larger cut-off is possible because the Lennard-Jones potential
and forces are almost zero beyond 1.4 nm.

GROMOS-96 files
^^^^^^^^^^^^^^^

GROMACS can read and write GROMOS-96 coordinate and trajectory files.
These files should have the extension .g96. Such a file can be a
GROMOS-96 initial/final configuration file, a coordinate trajectory
file, or a combination of both. The file is fixed format; all floats are
written as 15.9, and as such, files can get huge. GROMACS supports the
following data blocks in the given order:

-  Header block:

   ::

       TITLE (mandatory)

-  Frame blocks:

   ::

       TIMESTEP (optional)
       POSITION/POSITIONRED (mandatory)
       VELOCITY/VELOCITYRED (optional)
       BOX (optional)

See the GROMOS-96 manual \ `77 <#ref-gromos96>`__ for a complete
description of the blocks. **Note** that all GROMACS programs can read
compressed (.Z) or gzipped (.gz) files.

OPLS/AA
~~~~~~~

AMBER
~~~~~

GROMACS provides native support for the following AMBER force fields:

-  AMBER94 \ `111 <#ref-Cornell1995>`__

-  AMBER96 \ `112 <#ref-Kollman1996>`__

-  AMBER99 \ `113 <#ref-Wang2000>`__

-  AMBER99SB \ `114 <#ref-Hornak2006>`__

-  AMBER99SB-ILDN \ `115 <#ref-Lindorff2010>`__

-  AMBER03 \ `116 <#ref-Duan2003>`__

-  AMBERGS \ `117 <#ref-Garcia2002>`__

CHARMM
~~~~~~

GROMACS supports the CHARMM force field for
proteins \ `118 <#ref-mackerell04>`__, `119 <#ref-mackerell98>`__,
lipids \ `120 <#ref-feller00>`__ and nucleic
acids \ `121 <#ref-foloppe00>`__, `122 <#ref-Mac2000>`__. The protein
parameters (and to some extent the lipid and nucleic acid parameters)
were thoroughly tested – both by comparing potential energies between
the port and the standard parameter set in the CHARMM molecular
simulation package, as well by how the protein force field behaves
together with GROMACS-specific techniques such as virtual sites
(enabling long time steps) recently
implemented \ `123 <#ref-Larsson10>`__ – and the details and results are
presented in the paper by Bjelkmar et al. \ `124 <#ref-Bjelkmar10>`__.
The nucleic acid parameters, as well as the ones for HEME, were
converted and tested by Michel Cuendet.

When selecting the CHARMM force field in pdb2gmx the default option is
to use CMAP (for torsional correction map). To exclude CMAP, use
-nocmap. The basic form of the CMAP term implemented in GROMACS is a
function of the :math:`\phi` and :math:`\psi` backbone torsion angles.
This term is defined in the .rtp file by a statement at the end of each
residue supporting CMAP. The following five atom names define the two
torsional angles. Atoms 1-4 define :math:`\phi`, and atoms 2-5 define
:math:`\psi`. The corresponding atom types are then matched to the
correct CMAP type in the cmap.itp file that contains the correction
maps.

A port of the CHARMM36 force field for use with GROMACS is also
available at http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs.

For branched polymers or other topologies not supported by pdb2gmx, it
is possible to use TopoTools \ `125 <#ref-kohlmeyer2016>`__ to generate
a GROMACS top file.

Coarse-grained force fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~

[sec:cg-forcefields] Coarse-graining is a systematic way of reducing the
number of degrees of freedom representing a system of interest. To
achieve this, typically whole groups of atoms are represented by single
begineqnarrayds and the coarse-grained force fields describes their
effective interactions. Depending on the choice of parameterization, the
functional form of such an interaction can be complicated and often
tabulated potentials are used.

Coarse-grained models are designed to reproduce certain properties of a
reference system. This can be either a full atomistic model or even
experimental data. Depending on the properties to reproduce there are
different methods to derive such force fields. An incomplete list of
methods is given below:

-  Conserving free energies

   -  Simplex method

   -  MARTINI force field (see next section)

-  Conserving distributions (like the radial distribution function),
   so-called structure-based coarse-graining

   -  (iterative) Boltzmann inversion

   -  Inverse Monte Carlo

-  Conversing forces

   -  Force matching

Note that coarse-grained potentials are state dependent (e.g.
temperature, density,...) and should be re-parametrized depending on the
system of interest and the simulation conditions. This can for example
be done using the Versatile Object-oriented Toolkit for Coarse-Graining
Applications (VOTCA) **???**. The package was designed to assists in
systematic coarse-graining, provides implementations for most of the
algorithms mentioned above and has a well tested interface to GROMACS.
It is available as open source and further information can be found at
`www.votca.org <http://www.votca.org>`__.

MARTINI
~~~~~~~

The MARTINI force field is a coarse-grain parameter set that allows for
the construction of many systems, including proteins and membranes.

PLUM
~~~~

The PLUM force field `126 <#ref-bereau12>`__ is an example of a
solvent-free protein-membrane model for which the membrane was derived
from structure-based coarse-graining \ `127 <#ref-wang_jpcb10>`__. A
GROMACS implementation can be found at
`code.google.com/p/plumx <http://code.google.com/p/plumx/>`__.

Topologies
==========

Introduction
------------

GROMACS must know on which atoms and combinations of atoms the various
contributions to the potential functions (see chapter [ch:ff]) must act.
It must also know what parameters must be applied to the various
functions. All this is described in the *topology* file .top, which
lists the *constant attributes* of each atom. There are many more atom
types than elements, but only atom types present in biological systems
are parameterized in the force field, plus some metals, ions and
silicon. The bonded and special interactions are determined by fixed
lists that are included in the topology file. Certain non-bonded
interactions must be excluded (first and second neighbors), as these are
already treated in bonded interactions. In addition, there are *dynamic
attributes* of atoms - their positions, velocities and forces. These do
not strictly belong to the molecular topology, and are stored in the
coordinate file .gro (positions and velocities), or trajectory file .trr
(positions, velocities, forces).

This chapter describes the setup of the topology file, the .top file and
the database files: what the parameters stand for and how/where to
change them if needed. First, all file formats are explained. Section
[subsec:fffiles] describes the organization of the files in each force
field.

**Note:** if you construct your own topologies, we encourage you to
upload them to our topology archive at
`www.gromacs.org <http://www.gromacs.org>`__! Just imagine how thankful
you’d have been if your topology had been available there before you
started. The same goes for new force fields or modified versions of the
standard force fields - contribute them to the force field archive!

Particle type
-------------

In GROMACS, there are three types of particles, see Table [tab:ptype].
Only regular atoms and virtual interaction sites are used in GROMACS;
shells are necessary for polarizable models like the Shell-Water
models \ `45 <#ref-Maaren2001a>`__.

Atom types
~~~~~~~~~~

Each force field defines a set of atom types, which have a
characteristic name or number, and mass (in a.m.u.). These listings are
found in the atomtypes.atp file (.atp = **a**\ tom **t**\ ype
**p**\ arameter file). Therefore, it is in this file that you can begin
to change and/or add an atom type. A sample from the gromos43a1.ff force
field is listed below.

::

        O  15.99940 ;     carbonyl oxygen (C=O)
       OM  15.99940 ;     carboxyl oxygen (CO-)
       OA  15.99940 ;     hydroxyl, sugar or ester oxygen
       OW  15.99940 ;     water oxygen
        N  14.00670 ;     peptide nitrogen (N or NH)
       NT  14.00670 ;     terminal nitrogen (NH2)
       NL  14.00670 ;     terminal nitrogen (NH3)
       NR  14.00670 ;     aromatic nitrogen
       NZ  14.00670 ;     Arg NH (NH2)
       NE  14.00670 ;     Arg NE (NH)
        C  12.01100 ;     bare carbon
      CH1  13.01900 ;     aliphatic or sugar CH-group
      CH2  14.02700 ;     aliphatic or sugar CH2-group
      CH3  15.03500 ;     aliphatic CH3-group

**Note:** GROMACS makes use of the atom types as a name, *not* as a
number (as *e.g.* in GROMOS).

Virtual sites
~~~~~~~~~~~~~

Some force fields use virtual interaction sites (interaction sites that
are constructed from other particle positions) on which certain
interactions are located (*e.g.* on benzene rings, to reproduce the
correct quadrupole). This is described in sec. [sec:virtual\_sites].

To make virtual sites in your system, you should include a section (for
backward compatibility the old name can also be used) in your topology
file, where the ‘?’ stands for the number constructing particles for the
virtual site. This will be ‘2’ for type 2, ‘3’ for types 3, 3fd, 3fad
and 3out and ‘4’ for type 4fdn. The last of these replace an older 4fd
type (with the ‘type’ value 1) that could occasionally be unstable;
while it is still supported internally in the code, the old 4fd type
should not be used in new input files. The different types are explained
in sec. [sec:virtual\_sites].

Parameters for type 2 should look like this:

::

    [ virtual_sites2 ]
    ; Site  from        funct  a
    5       1     2     1      0.7439756

for type 3 like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          b
    5       1     2     3      1       0.7439756  0.128012

for type 3fd like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          d
    5       1     2     3      2       0.5        -0.105

for type 3fad like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   theta      d
    5       1     2     3      3       120        0.5

for type 3out like this:

::

    [ virtual_sites3 ]
    ; Site  from               funct   a          b          c
    5       1     2     3      4       -0.4       -0.4       6.9281

for type 4fdn like this:

::

    [ virtual_sites4 ]
    ; Site  from                      funct   a          b          c
    5       1     2     3     4       2       1.0        0.9       0.105

This will result in the construction of a virtual site, number 5 (first
column ‘Site’), based on the positions of the atoms whose indices are 1
and 2 or 1, 2 and 3 or 1, 2, 3 and 4 (next two, three or four columns
‘from’) following the rules determined by the function number (next
column ‘funct’) with the parameters specified (last one, two or three
columns ‘a b . .’). Obviously, the atom numbers (including virtual site
number) depend on the molecule. It may be instructive to study the
topologies for TIP4P or TIP5P water models that are included with the
GROMACS distribution.

**Note** that if any constant bonded interactions are defined between
virtual sites and/or normal atoms, they will be removed by grompp
(unless the option tt -normvsbds is used). This removal of bonded
interactions is done after generating exclusions, as the generation of
exclusions is based on “chemically” bonded interactions.

Virtual sites can be constructed in a more generic way using basic
geometric parameters. The directive that can be used is . Required
parameters are listed in Table [tab:topfile2]. An example entry for
defining a virtual site at the center of geometry of a given set of
atoms might be:

::

    [ virtual_sitesn ]
    ; Site   funct    from
    5        1        1     2     3     4

Parameter files
---------------

Atoms
~~~~~

The *static* properties (see Table [tab:statprop] assigned to the atom
types are assigned based on data in several places. The mass is listed
in atomtypes.atp (see [subsec:atomtype]), whereas the charge is listed
in .rtp (.rtp = **r**\ esidue **t**\ opology **p**\ arameter file,
see [subsec:rtp]). This implies that the charges are only defined in the
building blocks of amino acids, nucleic acids or otherwise, as defined
by the user. When generating a topology (.top) using the pdb2gmx
program, the information from these files is combined.

Non-bonded parameters
~~~~~~~~~~~~~~~~~~~~~

The non-bonded parameters consist of the van der Waals parameters V (c6
or :math:`\sigma`, depending on the combination rule) and W (c12 or
:math:`\epsilon`), as listed in the file ffnonbonded.itp, where ptype is
the particle type (see Table [tab:ptype]). As with the bonded
parameters, entries in directives are applied to their counterparts in
the topology file. Missing parameters generate warnings, except as noted
below in section [subsec:pairinteractions].

::

    [ atomtypes ]
    ;name   at.num      mass      charge   ptype         V(c6)        W(c12)
        O        8  15.99940       0.000       A   0.22617E-02   0.74158E-06
       OM        8  15.99940       0.000       A   0.22617E-02   0.74158E-06
       .....

    [ nonbond_params ]
      ; i    j func       V(c6)        W(c12)
        O    O    1 0.22617E-02   0.74158E-06
        O   OA    1 0.22617E-02   0.13807E-05
        .....

**Note** that most of the included force fields also include the at.num.
column, but this same information is implied in the OPLS-AA bond\_type
column. The interpretation of the parameters V and W depends on the
combination rule that was chosen in the section of the topology file
(see [subsec:topfile]):

.. math::

   \begin{aligned}
   \mbox{for combination rule 1}: & &
   \begin{array}{llllll}
     \mbox{V}_{ii} & = & C^{(6)}_{i}  & = & 4\,\epsilon_i\sigma_i^{6} &
     \mbox{[ kJ mol$^{-1}$ nm$^{6}$ ]}\\
     \mbox{W}_{ii} & = & C^{(12)}_{i} & = & 4\,\epsilon_i\sigma_i^{12} &
     \mbox{[ kJ mol$^{-1}$ nm$^{12}$ ]}\\
   \end{array}
   \\
   \mbox{for combination rules 2 and 3}: & &
   \begin{array}{llll}
     \mbox{V}_{ii} & = & \sigma_i   & \mbox{[ nm ]} \\
     \mbox{W}_{ii} & = & \epsilon_i & \mbox{[ kJ mol$^{-1}$ ]}
   \end{array}\end{aligned}

 Some or all combinations for different atom types can be given in the
section, again with parameters V and W as defined above. Any combination
that is not given will be computed from the parameters for the
corresponding atom types, according to the combination rule:

.. math::

   \begin{aligned}
   \mbox{for combination rules 1 and 3}: & &
   \begin{array}{lll}
     C^{(6)}_{ij}  & = & \left(C^{(6)}_i\,C^{(6)}_j\right)^{\frac{1}{2}} \\
     C^{(12)}_{ij} & = & \left(C^{(12)}_i\,C^{(12)}_j\right)^{\frac{1}{2}}
   \end{array}
   \\
   \mbox{for combination rule 2}: & &
   \begin{array}{lll}
     \sigma_{ij}   & = & \frac{1}{2}(\sigma_i+\sigma_j) \\
     \epsilon_{ij} & = & \sqrt{\epsilon_i\,\epsilon_j}
   \end{array}\end{aligned}

 When :math:`\sigma` and :math:`\epsilon` need to be supplied (rules 2
and 3), it would seem it is impossible to have a non-zero :math:`C^{12}`
combined with a zero :math:`C^6` parameter. However, providing a
negative :math:`\sigma` will do exactly that, such that :math:`C^6` is
set to zero and :math:`C^{12}` is calculated normally. This situation
represents a special case in reading the value of :math:`\sigma`, and
nothing more.

There is only one set of combination rules for Buckingham potentials:

.. math::

   \begin{array}{rcl}
   A_{ij}   &=& \left(A_{ii} \, A_{jj}\right)^{1/2}    \\
   B_{ij}   &=& 2 / \left(\frac{1}{B_{ii}} + \frac{1}{B_{jj}}\right)        \\
   C_{ij}   &=& \left(C_{ii} \, C_{jj}\right)^{1/2}
   \end{array}

Bonded parameters
~~~~~~~~~~~~~~~~~

The bonded parameters (*i.e.* bonds, bond angles, improper and proper
dihedrals) are listed in ffbonded.itp.  The entries in this database
describe, respectively, the atom types in the interactions, the type of
the interaction, and the parameters associated with that interaction.
These parameters are then read by grompp when processing a topology and
applied to the relevant bonded parameters, *i.e.* bondtypes are applied
to entries in the directive, etc. Any bonded parameter that is missing
from the relevant directive generates a fatal error. The types of
interactions are listed in Table [tab:topfile2]. Example excerpts from
such files follow:

::

    [ bondtypes ]
      ; i    j func        b0          kb
        C    O    1   0.12300     502080.
        C   OM    1   0.12500     418400.
        ......

    [ angletypes ]
      ; i    j    k func       th0         cth
       HO   OA    C    1   109.500     397.480
       HO   OA  CH1    1   109.500     397.480
       ......

    [ dihedraltypes ]
      ; i    l func        q0          cq
     NR5*  NR5    2     0.000     167.360
     NR5* NR5*    2     0.000     167.360
     ......

    [ dihedraltypes ]
      ; j    k func      phi0          cp   mult
        C   OA    1   180.000      16.736      2
        C    N    1   180.000      33.472      2
        ......

    [ dihedraltypes ]
    ;
    ; Ryckaert-Bellemans Dihedrals
    ;
    ; aj    ak      funct
    CP2     CP2     3       9.2789  12.156  -13.120 -3.0597 26.240  -31.495

In the ffbonded.itp file, you can add bonded parameters. If you want to
include parameters for new atom types, make sure you define them in
atomtypes.atp as well.

For most interaction types, bonded parameters are searched and assigned
using an exact match for all type names and allowing only a single set
of parameters. The exception to this rule are dihedral parameters. For
wildcard atom type names can be specified with the letter X in one or
more of the four positions. Thus one can for example assign proper
dihedral parameters based on the types of the middle two atoms. The
parameters for the entry with the most exact matches, i.e. the least
wildcard matches, will be used. Note that GROMACS versions older than
5.1.3 used the first match, which means that a full match would be
ignored if it is preceded by an entry that matches on wildcards. Thus it
is suggested to put wildcard entries at the end, in case someone might
use a forcefield with older versions of GROMACS. In addition there is a
dihedral type 9 which adds the possibility of assigning multiple
dihedral potentials, useful for combining terms with different
multiplicities. The different dihedral potential parameter sets should
be on directly adjacent lines in the section.

Molecule definition
-------------------

Moleculetype entries
~~~~~~~~~~~~~~~~~~~~

An organizational structure that usually corresponds to molecules is the
entry. This entry serves two main purposes. One is to give structure to
the topology file(s), usually corresponding to real molecules. This
makes the topology easier to read and writing it less labor intensive. A
second purpose is computational efficiency. The system definition that
is kept in memory is proportional in size of the moleculetype
definitions. If a molecule is present in 100000 copies, this saves a
factor of 100000 in memory, which means the system usually fits in
cache, which can improve performance tremendously. Interactions that
correspond to chemical bonds, that generate exclusions, can only be
defined between atoms within a moleculetype. It is allowed to have
multiple molecules which are not covalently bonded in one moleculetype
definition. Molecules can be made infinitely long by connecting to
themselves over periodic boundaries. When such periodic molecules are
present, an option in the mdp file needs to be set to tell GROMACS not
to attempt to make molecules that are broken over periodic boundaries
whole again.

Intermolecular interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some cases, one would like atoms in different molecules to also
interact with other interactions than the usual non-bonded interactions.
This is often the case in binding studies. When the molecules are
covalently bound, e.g. a ligand binding covalently to a protein, they
are effectively one molecule and they should be defined in one entry.
Note that pdb2gmx has an option to put two or more molecules in one
entry. When molecules are not covalently bound, it is much more
convenient to use separate moleculetype definitions and specify the
intermolecular interactions in the section. In this section, which is
placed at the end of the topology (see Table [tab:topfile1]), normal
bonded interactions can be specified using global atom indices. The only
restrictions are that no interactions can be used that generates
exclusions and no constraints can be used.

Intramolecular pair interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Extra Lennard-Jones and electrostatic interactions between pairs of
atoms in a molecule can be added in the section of a molecule
definition. The parameters for these interactions can be set
independently from the non-bonded interaction parameters. In the GROMOS
force fields, pairs are only used to modify the 1-4 interactions
(interactions of atoms separated by three bonds). In these force fields
the 1-4 interactions are excluded from the non-bonded interactions (see
sec. [sec:excl]).

::


    [ pairtypes ]
      ; i    j func         cs6          cs12 ; THESE ARE 1-4 INTERACTIONS
        O    O    1 0.22617E-02   0.74158E-06
        O   OM    1 0.22617E-02   0.74158E-06
        .....

The pair interaction parameters for the atom types in ffnonbonded.itp
are listed in the section. The GROMOS force fields list all these
interaction parameters explicitly, but this section might be empty for
force fields like OPLS that calculate the 1-4 interactions by uniformly
scaling the parameters. Pair parameters that are not present in the
section are only generated when gen-pairs is set to “yes” in the
directive of forcefield.itp (see [subsec:topfile]). When gen-pairs is
set to “no,” grompp will give a warning for each pair type for which no
parameters are given.

The normal pair interactions, intended for 1-4 interactions, have
function type 1. Function type 2 and the are intended for free-energy
simulations. When determining hydration free energies, the solute needs
to be decoupled from the solvent. This can be done by adding a B-state
topology (see sec. [sec:fecalc]) that uses zero for all solute
non-bonded parameters, *i.e.* charges and LJ parameters. However, the
free energy difference between the A and B states is not the total
hydration free energy. One has to add the free energy for reintroducing
the internal Coulomb and LJ interactions in the solute when in vacuum.
This second step can be combined with the first step when the Coulomb
and LJ interactions within the solute are not modified. For this
purpose, there is a pairs function type 2, which is identical to
function type 1, except that the B-state parameters are always identical
to the A-state parameters. For searching the parameters in the section,
no distinction is made between function type 1 and 2. The pairs section
is intended to replace the non-bonded interaction. It uses the unscaled
charges and the non-bonded LJ parameters; it also only uses the A-state
parameters. **Note** that one should add exclusions for all atom pairs
listed in , otherwise such pairs will also end up in the normal neighbor
lists.

Alternatively, this same behavior can be achieved without ever touching
the topology, by using the couple-moltype, couple-lambda0,
couple-lambda1, and couple-intramol keywords. See sections
sec. [sec:fecalc] and sec. [sec:dgimplement] for more information.

All three pair types always use plain Coulomb interactions, even when
Reaction-field, PME, Ewald or shifted Coulomb interactions are selected
for the non-bonded interactions. Energies for types 1 and 2 are written
to the energy and log file in separate “LJ-14” and “Coulomb-14” entries
per energy group pair. Energies for are added to the “LJ-(SR)” and
“Coulomb-(SR)” terms.

Exclusions
~~~~~~~~~~

The exclusions for non-bonded interactions are generated by grompp for
neighboring atoms up to a certain number of bonds away, as defined in
the section in the topology file (see [subsec:topfile]). Particles are
considered bonded when they are connected by “chemical” bonds ( types 1
to 5, 7 or 8) or constraints ( type 1). Type 5 can be used to create a
connection between two atoms without creating an interaction. There is a
harmonic interaction ( type 6) that does not connect the atoms by a
chemical bond. There is also a second constraint type ( type 2) that
fixes the distance, but does not connect the atoms by a chemical bond.
For a complete list of all these interactions, see Table [tab:topfile2].

Extra exclusions within a molecule can be added manually in a section.
Each line should start with one atom index, followed by one or more atom
indices. All non-bonded interactions between the first atom and the
other atoms will be excluded.

When all non-bonded interactions within or between groups of atoms need
to be excluded, is it more convenient and much more efficient to use
energy monitor group exclusions (see sec. [sec:groupconcept]).

Constraint algorithms
---------------------

Constraints are defined in the section. The format is two atom numbers
followed by the function type, which can be 1 or 2, and the constraint
distance. The only difference between the two types is that type 1 is
used for generating exclusions and type 2 is not (see sec. [sec:excl]).
The distances are constrained using the LINCS or the SHAKE algorithm,
which can be selected in the .mdp file. Both types of constraints can be
perturbed in free-energy calculations by adding a second constraint
distance (see [subsec:constraintforce]). Several types of bonds and
angles (see Table [tab:topfile2]) can be converted automatically to
constraints by grompp. There are several options for this in the .mdp
file.

We have also implemented the SETTLE
algorithm \ `47 <#ref-Miyamoto92>`__, which is an analytical solution of
SHAKE, specifically for water. SETTLE can be selected in the topology
file. See, for instance, the SPC molecule definition:

::

    [ moleculetype ]
    ; molname       nrexcl
    SOL             1

    [ atoms ]
    ; nr    at type res nr  ren nm  at nm   cg nr   charge
    1       OW      1       SOL     OW1     1       -0.82
    2       HW      1       SOL     HW2     1        0.41
    3       HW      1       SOL     HW3     1        0.41

    [ settles ]
    ; OW    funct   doh     dhh
    1       1       0.1     0.16333

    [ exclusions ]
    1       2       3
    2       1       3
    3       1       2

The directive defines the first atom of the water molecule. The settle
funct is always 1, and the distance between O-H and H-H distances must
be given. **Note** that the algorithm can also be used for TIP3P and
TIP4P \ `128 <#ref-Jorgensen83>`__. TIP3P just has another geometry.
TIP4P has a virtual site, but since that is generated it does not need
to be shaken (nor stirred).

pdb2gmx input files
-------------------

The GROMACS program pdb2gmx generates a topology for the input
coordinate file. Several formats are supported for that coordinate file,
but .pdb is the most commonly-used format (hence the name pdb2gmx).
pdb2gmx searches for force fields in sub-directories of the GROMACS
share/top directory and your working directory. Force fields are
recognized from the file forcefield.itp in a directory with the
extension .ff. The file forcefield.doc may be present, and if so, its
first line will be used by pdb2gmx to present a short description to the
user to help in choosing a force field. Otherwise, the user can choose a
force field with the -ff xxx command-line argument to pdb2gmx, which
indicates that a force field in a xxx.ff directory is desired. pdb2gmx
will search first in the working directory, then in the GROMACS
share/top directory, and use the first matching xxx.ff directory found.

Two general files are read by pdb2gmx: an atom type file (extension
.atp, see [subsec:atomtype]) from the force-field directory, and a file
called residuetypes.dat from either the working directory, or the
GROMACS share/top directory. residuetypes.dat determines which residue
names are considered protein, DNA, RNA, water, and ions.

pdb2gmx can read one or multiple databases with topological information
for different types of molecules. A set of files belonging to one
database should have the same basename, preferably telling something
about the type of molecules (*e.g.* aminoacids, rna, dna). The possible
files are:

-  <basename>.rtp

-  <basename>.r2b (optional)

-  <basename>.arn (optional)

-  <basename>.hdb (optional)

-  <basename>.n.tdb (optional)

-  <basename>.c.tdb (optional)

Only the .rtp file, which contains the topologies of the building
blocks, is mandatory. Information from other files will only be used for
building blocks that come from an .rtp file with the same base name. The
user can add building blocks to a force field by having additional files
with the same base name in their working directory. By default, only
extra building blocks can be defined, but calling pdb2gmx with the -rtpo
option will allow building blocks in a local file to replace the default
ones in the force field.

Residue database
~~~~~~~~~~~~~~~~

The files holding the residue databases have the extension .rtp.
Originally this file contained building blocks (amino acids) for
proteins, and is the GROMACS interpretation of the rt37c4.dat file of
GROMOS. So the residue database file contains information (bonds,
charges, charge groups, and improper dihedrals) for a frequently-used
building block. It is better *not* to change this file because it is
standard input for pdb2gmx, but if changes are needed make them in the
.top file (see [subsec:topfile]), or in a .rtp file in the working
directory as explained in sec. [sec:pdb2gmxfiles]. Defining topologies
of new small molecules is probably easier by writing an include topology
file .itp directly. This will be discussed in section [subsec:molitp].
When adding a new protein residue to the database, don’t forget to add
the residue name to the residuetypes.dat file, so that grompp, make\_ndx
and analysis tools can recognize the residue as a protein residue (see
[subsec:defaultgroups]).

The .rtp files are only used by pdb2gmx. As mentioned before, the only
extra information this program needs from the .rtp database is bonds,
charges of atoms, charge groups, and improper dihedrals, because the
rest is read from the coordinate input file. Some proteins contain
residues that are not standard, but are listed in the coordinate file.
You have to construct a building block for this “strange” residue,
otherwise you will not obtain a .top file. This also holds for molecules
in the coordinate file such as ligands, polyatomic ions, crystallization
co-solvents, etc. The residue database is constructed in the following
way:

::

    [ bondedtypes ]  ; mandatory
    ; bonds  angles  dihedrals  impropers
         1       1          1          2  ; mandatory

    [ GLY ]  ; mandatory

     [ atoms ]  ; mandatory 
    ; name  type  charge  chargegroup 
         N     N  -0.280     0
         H     H   0.280     0
        CA   CH2   0.000     1
         C     C   0.380     2
         O     O  -0.380     2

     [ bonds ]  ; optional
    ;atom1 atom2      b0      kb
         N     H
         N    CA
        CA     C
         C     O
        -C     N

     [ exclusions ]  ; optional
    ;atom1 atom2

     [ angles ]  ; optional
    ;atom1 atom2 atom3    th0    cth

     [ dihedrals ]  ; optional
    ;atom1 atom2 atom3 atom4   phi0     cp   mult

     [ impropers ]  ; optional
    ;atom1 atom2 atom3 atom4     q0     cq
         N    -C    CA     H
        -C   -CA     N    -O

    [ ZN ]

     [ atoms ]
        ZN    ZN   2.000     0

The file is free format; the only restriction is that there can be at
most one entry on a line. The first field in the file is the field,
which is followed by four numbers, indicating the interaction type for
bonds, angles, dihedrals, and improper dihedrals. The file contains
residue entries, which consist of atoms and (optionally) bonds, angles,
dihedrals, and impropers. The charge group codes denote the charge group
numbers. Atoms in the same charge group should always be ordered
consecutively. When using the hydrogen database with pdb2gmx for adding
missing hydrogens (see [subsec:hdb]), the atom names defined in the .rtp
entry should correspond exactly to the naming convention used in the
hydrogen database. The atom names in the bonded interaction can be
preceded by a minus or a plus, indicating that the atom is in the
preceding or following residue respectively. Explicit parameters added
to bonds, angles, dihedrals, and impropers override the standard
parameters in the .itp files. This should only be used in special cases.
Instead of parameters, a string can be added for each bonded
interaction. This is used in GROMOS-96 .rtp files. These strings are
copied to the topology file and can be replaced by force-field
parameters by the C-preprocessor in grompp using #define statements.

pdb2gmx automatically generates all angles. This means that for most
force fields the field is only useful for overriding .itp parameters.
For the GROMOS-96 force field the interaction number of all angles needs
to be specified.

pdb2gmx automatically generates one proper dihedral for every rotatable
bond, preferably on heavy atoms. When the field is used, no other
dihedrals will be generated for the bonds corresponding to the specified
dihedrals. It is possible to put more than one dihedral function on a
rotatable bond. In the case of CHARMM27 FF pdb2gmx can add correction
maps to the dihedrals using the default -cmap option. Please refer to
[subsec:charmmff] for more information.

pdb2gmx sets the number of exclusions to 3, which means that
interactions between atoms connected by at most 3 bonds are excluded.
Pair interactions are generated for all pairs of atoms that are
separated by 3 bonds (except pairs of hydrogens). When more interactions
need to be excluded, or some pair interactions should not be generated,
an field can be added, followed by pairs of atom names on separate
lines. All non-bonded and pair interactions between these atoms will be
excluded.

Residue to building block database
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Each force field has its own naming convention for residues. Most
residues have consistent naming, but some, especially those with
different protonation states, can have many different names. The .r2b
files are used to convert standard residue names to the force-field
build block names. If no .r2b is present in the force-field directory or
a residue is not listed, the building block name is assumed to be
identical to the residue name. The .r2b can contain 2 or 5 columns. The
2-column format has the residue name in the first column and the
building block name in the second. The 5-column format has 3 additional
columns with the building block for the residue occurring in the
N-terminus, C-terminus and both termini at the same time (single residue
molecule). This is useful for, for instance, the AMBER force fields. If
one or more of the terminal versions are not present, a dash should be
entered in the corresponding column.

There is a GROMACS naming convention for residues which is only apparent
(except for the pdb2gmx code) through the .r2b file and specbond.dat
files. This convention is only of importance when you are adding residue
types to an .rtp file. The convention is listed in Table [tab:r2b]. For
special bonds with, for instance, a heme group, the GROMACS naming
convention is introduced through specbond.dat (see [subsec:specbond]),
which can subsequently be translated by the .r2b file, if required.

Atom renaming database
~~~~~~~~~~~~~~~~~~~~~~

Force fields often use atom names that do not follow IUPAC or PDB
convention. The .arn database is used to translate the atom names in the
coordinate file to the force-field names. Atoms that are not listed keep
their names. The file has three columns: the building block name, the
old atom name, and the new atom name, respectively. The residue name
supports question-mark wildcards that match a single character.

An additional general atom renaming file called xlateat.dat is present
in the share/top directory, which translates common non-standard atom
names in the coordinate file to IUPAC/PDB convention. Thus, when writing
force-field files, you can assume standard atom names and no further
atom name translation is required, except for translating from standard
atom names to the force-field ones.

Hydrogen database
~~~~~~~~~~~~~~~~~

The hydrogen database is stored in .hdb files. It contains information
for the pdb2gmx program on how to connect hydrogen atoms to existing
atoms. In versions of the database before GROMACS 3.3, hydrogen atoms
were named after the atom they are connected to: the first letter of the
atom name was replaced by an ‘H.’ In the versions from 3.3 onwards, the
H atom has to be listed explicitly, because the old behavior was
protein-specific and hence could not be generalized to other molecules.
If more than one hydrogen atom is connected to the same atom, a number
will be added to the end of the hydrogen atom name. For example, adding
two hydrogen atoms to ``ND2`` (in asparagine), the hydrogen atoms will
be named ``HD21`` and ``HD22``. This is important since atom naming in
the ``.rtp`` file (see [subsec:rtp]) must be the same. The format of the
hydrogen database is as follows:

::

    ; res   # additions
            # H add type    H       i       j       k
    ALA     1
            1       1       H       N       -C      CA
    ARG     4
            1       2       H       N       CA      C
            1       1       HE      NE      CD      CZ
            2       3       HH1     NH1     CZ      NE
            2       3       HH2     NH2     CZ      NE

On the first line we see the residue name (ALA or ARG) and the number of
kinds of hydrogen atoms that may be added to this residue by the
hydrogen database. After that follows one line for each addition, on
which we see:

-  The number of H atoms added

-  The method for adding H atoms, which can be any of:

   #. | *one planar hydrogen, *e.g.* rings or peptide bond*
      | One hydrogen atom (n) is generated, lying in the plane of atoms
        (i,j,k) on the plane bisecting angle (j-i-k) at a distance of
        0.1 nm from atom i, such that the angles (n-i-j) and (n-i-k) are
        :math:`>` 90\ :math:`^{\rm o}`.

   #. | *one single hydrogen, *e.g.* hydroxyl*
      | One hydrogen atom (n) is generated at a distance of 0.1 nm from
        atom i, such that angle (n-i-j)=109.5 degrees and dihedral
        (n-i-j-k)=trans.

   #. | *two planar hydrogens, *e.g.* ethylene -C=CH:math:`_2`, or amide
        -C(=O)NH:math:`_2`*
      | Two hydrogens (n1,n2) are generated at a distance of 0.1 nm from
        atom i, such that angle (n1-i-j)=(n2-i-j)=120 degrees and
        dihedral (n1-i-j-k)=cis and (n2-i-j-k)=trans, such that names
        are according to IUPAC standards \ `129 <#ref-iupac70>`__.

   #. | *two or three tetrahedral hydrogens, *e.g.* -CH:math:`_3`*
      | Three (n1,n2,n3) or two (n1,n2) hydrogens are generated at a
        distance of 0.1 nm from atom i, such that angle
        (n1-i-j)=(n2-i-j)=(n3-i-j)=109.47:math:`^{\rm o}`, dihedral
        (n1-i-j-k)=trans, (n2-i-j-k)=trans+120 and
        (n3-i-j-k)=trans+240:math:`^{\rm o}`.

   #. | *one tetrahedral hydrogen, *e.g.* C\ :math:`_3`\ CH*
      | One hydrogen atom (n:math:`^\prime`) is generated at a distance
        of 0.1 nm from atom i in tetrahedral conformation such that
        angle
        (n:math:`^\prime`-i-j)=(n:math:`^\prime`-i-k)=(n:math:`^\prime`-i-l)=109.47:math:`^{\rm o}`.

   #. | *two tetrahedral hydrogens, *e.g.* C-CH\ :math:`_2`-C*
      | Two hydrogen atoms (n1,n2) are generated at a distance of 0.1 nm
        from atom i in tetrahedral conformation on the plane bisecting
        angle j-i-k with angle
        (n1-i-n2)=(n1-i-j)=(n1-i-k)=109.47:math:`^{\rm o}`.

   #. | *two water hydrogens*
      | Two hydrogens are generated around atom i according to
        SPC \ `80 <#ref-Berendsen81>`__ water geometry. The symmetry
        axis will alternate between three coordinate axes in both
        directions.

   #. | *three water “hydrogens”*
      | Two hydrogens are generated around atom i according to
        SPC \ `80 <#ref-Berendsen81>`__ water geometry. The symmetry
        axis will alternate between three coordinate axes in both
        directions. In addition, an extra particle is generated on the
        position of the oxygen with the first letter of the name
        replaced by ‘M’. This is for use with four-atom water models
        such as TIP4P \ `128 <#ref-Jorgensen83>`__.

   #. | *four water “hydrogens”*
      | Same as above, except that two additional particles are
        generated on the position of the oxygen, with names ‘LP1’ and
        ‘LP2.’ This is for use with five-atom water models such as
        TIP5P \ `130 <#ref-Mahoney2000a>`__.

-  The name of the new H atom (or its prefix, *e.g.* HD2 for the
   asparagine example given earlier).

-  Three or four control atoms (i,j,k,l), where the first always is the
   atom to which the H atoms are connected. The other two or three
   depend on the code selected. For water, there is only one control
   atom.

Some more exotic cases can be approximately constructed from the above
tools, and with suitable use of energy minimization are good enough for
beginning MD simulations. For example secondary amine hydrogen, nitrenyl
hydrogen (C=NH) and even ethynyl hydrogen could be approximately
constructed using method 2 above for hydroxyl hydrogen.

Termini database
~~~~~~~~~~~~~~~~

The termini databases are stored in aminoacids.n.tdb and
aminoacids.c.tdb for the N- and C-termini respectively. They contain
information for the pdb2gmx program on how to connect new atoms to
existing ones, which atoms should be removed or changed, and which
bonded interactions should be added. Their format is as follows (from
gromos43a1.ff/aminoacids.c.tdb):

::

    [ None ]

    [ COO- ]
    [ replace ]
    C	C	C	12.011	0.27
    O 	O1	OM	15.9994	-0.635
    OXT	O2	OM	15.9994	-0.635
    [ add ]
    2	8	O	C	CA	N
    	OM	15.9994	-0.635
    [ bonds ]
    C	O1	gb_5
    C	O2	gb_5
    [ angles ]
    O1	C	O2	ga_37
    CA	C	O1	ga_21
    CA	C	O2	ga_21
    [ dihedrals ]
    N	CA	C	O2	gd_20
    [ impropers ]
    C	CA	O2	O1	gi_1

The file is organized in blocks, each with a header specifying the name
of the block. These blocks correspond to different types of termini that
can be added to a molecule. In this example is the first block,
corresponding to changing the terminal carbon atom into a deprotonated
carboxyl group. is the second terminus type, corresponding to a terminus
that leaves the molecule as it is. Block names cannot be any of the
following: replace, add, delete, bonds, angles, dihedrals, impropers.
Doing so would interfere with the parameters of the block, and would
probably also be very confusing to human readers.

For each block the following options are present:

-  | 
   | Replace an existing atom by one with a different atom type, atom
     name, charge, and/or mass. This entry can be used to replace an
     atom that is present both in the input coordinates and in the .rtp
     database, but also to only rename an atom in the input coordinates
     such that it matches the name in the force field. In the latter
     case, there should also be a corresponding section present that
     gives instructions to add the same atom, such that the position in
     the sequence and the bonding is known. Such an atom can be present
     in the input coordinates and kept, or not present and constructed
     by pdb2gmx. For each atom to be replaced on line should be entered
     with the following fields:

   -  name of the atom to be replaced

   -  new atom name (optional)

   -  new atom type

   -  new mass

   -  new charge

-  | 
   | Add new atoms. For each (group of) added atom(s), a two-line entry
     is necessary. The first line contains the same fields as an entry
     in the hydrogen database (name of the new atom, number of atoms,
     type of addition, control atoms, see [subsec:hdb]), but the
     possible types of addition are extended by two more, specifically
     for C-terminal additions:

   #. | *two carboxyl oxygens, -COO:math:`^-`*
      | Two oxygens (n1,n2) are generated according to rule 3, at a
        distance of 0.136 nm from atom i and an angle
        (n1-i-j)=(n2-i-j)=117 degrees

   #. | *carboxyl oxygens and hydrogen, -COOH*
      | Two oxygens (n1,n2) are generated according to rule 3, at
        distances of 0.123 nm and 0.125 nm from atom i for n1 and n2,
        respectively, and angles (n1-i-j)=121 and (n2-i-j)=115 degrees.
        One hydrogen (n:math:`^\prime`) is generated around n2 according
        to rule 2, where n-i-j and n-i-j-k should be read as
        n\ :math:`^\prime`-n2-i and n\ :math:`^\prime`-n2-i-j,
        respectively.

   After this line, another line follows that specifies the details of
   the added atom(s), in the same way as for replacing atoms, *i.e.*:

   -  atom type

   -  mass

   -  charge

   -  charge group (optional)

   Like in the hydrogen database (see [subsec:rtp]), when more than one
   atom is connected to an existing one, a number will be appended to
   the end of the atom name. **Note** that, like in the hydrogen
   database, the atom name is now on the same line as the control atoms,
   whereas it was at the beginning of the second line prior to GROMACS
   version 3.3. When the charge group field is left out, the added atom
   will have the same charge group number as the atom that it is bonded
   to.

-  | 
   | Delete existing atoms. One atom name per line.

-  | , , and
   | Add additional bonded parameters. The format is identical to that
     used in the .rtp file, see [subsec:rtp].

Virtual site database
~~~~~~~~~~~~~~~~~~~~~

Since we cannot rely on the positions of hydrogens in input files, we
need a special input file to decide the geometries and parameters with
which to add virtual site hydrogens. For more complex virtual site
constructs (*e.g.* when entire aromatic side chains are made rigid) we
also need information about the equilibrium bond lengths and angles for
all atoms in the side chain. This information is specified in the .vsd
file for each force field. Just as for the termini, there is one such
file for each class of residues in the .rtp file.

The virtual site database is not really a very simple list of
information. The first couple of sections specify which mass centers
(typically called MCH\ :math:`_3`/MNH:math:`_3`) to use for
CH\ :math:`_3`, NH\ :math:`_3`, and NH\ :math:`_2` groups. Depending on
the equilibrium bond lengths and angles between the hydrogens and heavy
atoms we need to apply slightly different constraint distances between
these mass centers. **Note** that we do *not* have to specify the actual
parameters (that is automatic), just the type of mass center to use. To
accomplish this, there are three sections names ``[ CH3 ]``,
``[ NH3 ]``, and ``[ NH2 ]``. For each of these we expect three columns.
The first column is the atom type bound to the 2/3 hydrogens, the second
column is the next heavy atom type which this is bound, and the third
column the type of mass center to use. As a special case, in the
``[ NH2 ]`` section it is also possible to specify ``planar`` in the
second column, which will use a different construction without mass
center. There are currently different opinions in some force fields
whether an NH\ :math:`_2` group should be planar or not, but we try hard
to stick to the default equilibrium parameters of the force field.

The second part of the virtual site database contains explicit
equilibrium bond lengths and angles for pairs/triplets of atoms in
aromatic side chains. These entries are currently read by specific
routines in the virtual site generation code, so if you would like to
extend it *e.g.* to nucleic acids you would also need to write new code
there. These sections are named after the short amino acid names
(``[ PHE ]``, ``[ TYR ]``, ``[ TRP ]``, ``[ HID ]``, ``[ HIE ]``,
``[ HIP ]``), and simply contain 2 or 3 columns with atom names,
followed by a number specifying the bond length (in nm) or angle (in
degrees). **Note** that these are approximations of the equilibrated
geometry for the entire molecule, which might not be identical to the
equilibrium value for a single bond/angle if the molecule is strained.

Special bonds
~~~~~~~~~~~~~

The primary mechanism used by pdb2gmx to generate inter-residue bonds
relies on head-to-tail linking of backbone atoms in different residues
to build a macromolecule. In some cases (*e.g.* disulfide bonds, a heme
group, branched polymers), it is necessary to create inter-residue bonds
that do not lie on the backbone. The file specbond.dat takes care of
this function. It is necessary that the residues belong to the same .
The -merge and -chainsep functions of pdb2gmx can be useful when
managing special inter-residue bonds between different chains.

The first line of specbond.dat indicates the number of entries that are
in the file. If you add a new entry, be sure to increment this number.
The remaining lines in the file provide the specifications for creating
bonds. The format of the lines is as follows:

resA atomA nbondsA resB atomB nbondsB length newresA newresB

The columns indicate:

#. resA The name of residue A that participates in the bond.

#. atomA The name of the atom in residue A that forms the bond.

#. nbondsA The total number of bonds atomA can form.

#. resB The name of residue B that participates in the bond.

#. atomB The name of the atom in residue B that forms the bond.

#. nbondsB The total number of bonds atomB can form.

#. length The reference length for the bond. If atomA and atomB are not
   within length :math:`\pm` 10% in the coordinate file supplied to
   pdb2gmx, no bond will be formed.

#. newresA The new name of residue A, if necessary. Some force fields
   use *e.g.* CYS2 for a cysteine in a disulfide or heme linkage.

#. newresB The new name of residue B, likewise.

File formats
------------

Topology file
~~~~~~~~~~~~~

The topology file is built following the GROMACS specification for a
molecular topology. A .top file can be generated by pdb2gmx. All
possible entries in the topology file are listed in Tables
[tab:topfile1] and [tab:topfile2]. Also tabulated are: all the units of
the parameters, which interactions can be perturbed for free energy
calculations, which bonded interactions are used by grompp for
generating exclusions, and which bonded interactions can be converted to
constraints by grompp.

Description of the file layout:

-  Semicolon (;) and newline characters surround comments

-  On a line ending with :math:`\backslash` the newline character is
   ignored.

-  Directives are surrounded by [ and ]

-  The topology hierarchy (which must be followed) consists of three
   levels:

   -  the parameter level, which defines certain force-field
      specifications (see Table [tab:topfile1])

   -  the molecule level, which should contain one or more molecule
      definitions (see Table [tab:topfile2])

   -  the system level, containing only system-specific information (
      and )

-  Items should be separated by spaces or tabs, not commas

-  Atoms in molecules should be numbered consecutively starting at 1

-  Atoms in the same charge group must be listed consecutively

-  The file is parsed only once, which implies that no forward
   references can be treated: items must be defined before they can be
   used

-  Exclusions can be generated from the bonds or overridden manually

-  The bonded force types can be generated from the atom types or
   overridden per bond

-  It is possible to apply multiple bonded interactions of the same type
   on the same atoms

-  Descriptive comment lines and empty lines are highly recommended

-  Starting with GROMACS version 3.1.3, all directives at the parameter
   level can be used multiple times and there are no restrictions on the
   order, except that an atom type needs to be defined before it can be
   used in other parameter definitions

-  If parameters for a certain interaction are defined multiple times
   for the same combination of atom types the last definition is used;
   starting with GROMACS version 3.1.3 grompp generates a warning for
   parameter redefinitions with different values

-  Using one of the , , , , etc. without having used before is
   meaningless and generates a warning

-  Using without having used before is meaningless and generates a
   warning.

-  After the only allowed directive is

-  Using an unknown string in causes all the data until the next
   directive to be ignored and generates a warning

Here is an example of a topology file, urea.top:

::

    ;
    ;       Example topology file
    ;
    ; The force-field files to be included
    #include "amber99.ff/forcefield.itp"

    [ moleculetype ]
    ; name  nrexcl
    Urea         3

    [ atoms ]
       1  C  1  URE      C      1     0.880229  12.01000   ; amber C  type
       2  O  1  URE      O      2    -0.613359  16.00000   ; amber O  type
       3  N  1  URE     N1      3    -0.923545  14.01000   ; amber N  type
       4  H  1  URE    H11      4     0.395055   1.00800   ; amber H  type
       5  H  1  URE    H12      5     0.395055   1.00800   ; amber H  type
       6  N  1  URE     N2      6    -0.923545  14.01000   ; amber N  type
       7  H  1  URE    H21      7     0.395055   1.00800   ; amber H  type
       8  H  1  URE    H22      8     0.395055   1.00800   ; amber H  type

    [ bonds ]
        1	2
        1	3	
        1   6
        3	4
        3	5
        6	7
        6	8

    [ dihedrals ] 
    ;   ai    aj    ak    al funct  definition
         2     1     3     4   9     
         2     1     3     5   9     
         2     1     6     7   9     
         2     1     6     8   9     
         3     1     6     7   9     
         3     1     6     8   9     
         6     1     3     4   9     
         6     1     3     5   9     

    [ dihedrals ] 
         3     6     1     2   4     
         1     4     3     5   4	 
         1     7     6     8   4

    [ position_restraints ]
    ; you wouldn't normally use this for a molecule like Urea,
    ; but we include it here for didactic purposes
    ; ai   funct    fc
       1     1     1000    1000    1000 ; Restrain to a point
       2     1     1000       0    1000 ; Restrain to a line (Y-axis)
       3     1     1000       0       0 ; Restrain to a plane (Y-Z-plane)

    [ dihedral_restraints ]
    ; ai   aj    ak    al  type  phi  dphi  fc
        3    6     1    2     1  180     0  10
        1    4     3    5     1  180     0  10

    ; Include TIP3P water topology
    #include "amber99/tip3p.itp"

    [ system ]
    Urea in Water

    [ molecules ]
    ;molecule name   nr.
    Urea             1
    SOL              1000

Here follows the explanatory text.

**#include “amber99.ff/forcefield.itp” :** this includes the information
for the force field you are using, including bonded and non-bonded
parameters. This example uses the AMBER99 force field, but your
simulation may use a different force field. grompp will automatically go
and find this file and copy-and-paste its content. That content can be
seen in , and it is

::

    #define _FF_AMBER
    #define _FF_AMBER99

    [ defaults ]
    ; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ
    1               2               yes             0.5     0.8333

    #include "ffnonbonded.itp"
    #include "ffbonded.itp"

The two #define statements set up the conditions so that future parts of
the topology can know that the AMBER 99 force field is in use.

 **:**

-  nbfunc is the non-bonded function type. Use 1 (Lennard-Jones) or 2
   (Buckingham)

-  comb-rule is the number of the combination rule (see [subsec:nbpar]).

-  gen-pairs is for pair generation. The default is ‘no’, *i.e.* get 1-4
   parameters from the pairtypes list. When parameters are not present
   in the list, stop with a fatal error. Setting ‘yes’ generates 1-4
   parameters that are not present in the pair list from normal
   Lennard-Jones parameters using fudgeLJ

-  fudgeLJ is the factor by which to multiply Lennard-Jones 1-4
   interactions, default 1

-  fudgeQQ is the factor by which to multiply electrostatic 1-4
   interactions, default 1

-  :math:`N` is the power for the repulsion term in a 6-\ :math:`N`
   potential (with nonbonded-type Lennard-Jones only), starting with
   GROMACS version 4.5, mdrun also reads and applies :math:`N`, for
   values not equal to 12 tabulated interaction functions are used (in
   older version you would have to use user tabulated interactions).

**Note** that gen-pairs, fudgeLJ, fudgeQQ, and :math:`N` are optional.
fudgeLJ is only used when generate pairs is set to ‘yes’, and fudgeQQ is
always used. However, if you want to specify :math:`N` you need to give
a value for the other parameters as well.

Then some other #include statements add in the large amount of data
needed to describe the rest of the force field. We will skip these and
return to urea.top. There we will see

 **:** defines the name of your molecule in this .top and nrexcl = 3
stands for excluding non-bonded interactions between atoms that are no
further than 3 bonds away.

 **:** defines the molecule, where nr and type are fixed, the rest is
user defined. So atom can be named as you like, cgnr made larger or
smaller (if possible, the total charge of a charge group should be
zero), and charges can be changed here too.

 **:** no comment.

 **:** LJ and Coulomb 1-4 interactions

 **:** no comment

 **:** in this case there are 9 proper dihedrals (funct = 1), 3 improper
(funct = 4) and no Ryckaert-Bellemans type dihedrals. If you want to
include Ryckaert-Bellemans type dihedrals in a topology, do the
following (in case of *e.g.* decane):

::

    [ dihedrals ]
    ;  ai    aj    ak    al funct       c0       c1       c2
        1    2     3     4     3 
        2    3     4     5     3

In the original implementation of the potential for
alkanes \ `131 <#ref-Ryckaert78>`__ no 1-4 interactions were used, which
means that in order to implement that particular force field you need to
remove the 1-4 interactions from the section of your topology. In most
modern force fields, like OPLS/AA or Amber the rules are different, and
the Ryckaert-Bellemans potential is used as a cosine series in
combination with 1-4 interactions.

 **:** harmonically restrain the selected particles to reference
positions ([subsec:positionrestraint]). The reference positions are read
from a separate coordinate file by grompp.

 **:** restrain selected dihedrals to a reference value. The
implementation of dihedral restraints is described in section
[subsec:dihedralrestraint] of the manual. The parameters specified in
the [dihedral\_restraints] directive are as follows:

-  type has only one possible value which is 1

-  phi is the value of :math:`\phi_0` in eqn. [eqn:dphi] and
   eqn. [eqn:dihre] of the manual.

-  dphi is the value of :math:`\Delta\phi` in eqn. [eqn:dihre] of the
   manual.

-  fc is the force constant :math:`k_{dihr}` in eqn. [eqn:dihre] of the
   manual.

**#include “tip3p.itp” :** includes a topology file that was already
constructed (see section [subsec:molitp]).

 **:** title of your system, user-defined

 **:** this defines the total number of (sub)molecules in your system
that are defined in this .top. In this example file, it stands for 1
urea molecule dissolved in 1000 water molecules. The molecule type SOL
is defined in the tip3p.itp file. Each name here must correspond to a
name given with earlier in the topology. The order of the blocks of
molecule types and the numbers of such molecules must match the
coordinate file that accompanies the topology when supplied to grompp.
The blocks of molecules do not need to be contiguous, but some tools
(e.g. genion) may act only on the first or last such block of a
particular molecule type. Also, these blocks have nothing to do with the
definition of groups (see sec. [sec:groupconcept] and
sec. [sec:usinggroups]).

Molecule.itp file
~~~~~~~~~~~~~~~~~

If you construct a topology file you will use frequently (like the water
molecule, tip3p.itp, which is already constructed for you) it is good to
make a molecule.itp file. This only lists the information of one
particular molecule and allows you to re-use the in multiple systems
without re-invoking pdb2gmx or manually copying and pasting. An example
urea.itp follows:

::

    [ moleculetype ]
    ; molname	nrexcl
    URE		3

    [ atoms ]
       1  C  1  URE      C      1     0.880229  12.01000   ; amber C  type
    ...
       8  H  1  URE    H22      8     0.395055   1.00800   ; amber H  type

    [ bonds ]
        1	2
    ...
        6	8
    [ dihedrals ] 
    ;   ai    aj    ak    al funct  definition
         2     1     3     4   9     
    ...
         6     1     3     5   9     
    [ dihedrals ] 
         3     6     1     2   4     
         1     4     3     5   4	 
         1     7     6     8   4

Using .itp files results in a very short .top file:

::

    ;
    ;       Example topology file
    ;
    ; The force field files to be included
    #include "amber99.ff/forcefield.itp"

    #include "urea.itp"

    ; Include TIP3P water topology
    #include "amber99/tip3p.itp"

    [ system ]
    Urea in Water

    [ molecules ]
    ;molecule name   nr.
    Urea             1
    SOL              1000

Ifdef statements
~~~~~~~~~~~~~~~~

A very powerful feature in GROMACS is the use of #ifdef statements in
your .top file. By making use of this statement, and associated #define
statements like were seen in earlier, different parameters for one
molecule can be used in the same .top file. An example is given for TFE,
where there is an option to use different charges on the atoms: charges
derived by De Loof *et al.* `132 <#ref-Loof92>`__ or by Van Buuren and
Berendsen \ `133 <#ref-Buuren93a>`__. In fact, you can use much of the
functionality of the C preprocessor, cpp, because grompp contains
similar pre-processing functions to scan the file. The way to make use
of the #ifdef option is as follows:

-  either use the option define = -DDeLoof in the .mdp file (containing
   grompp input parameters), or use the line #define DeLoof early in
   your .top or .itp file; and

-  put the #ifdef statements in your .top, as shown below:

::

    ...



    [ atoms ]
    ; nr     type     resnr    residu     atom      cgnr      charge        mass
    #ifdef DeLoof
    ; Use Charges from DeLoof
       1        C        1        TFE        C         1        0.74        
       2        F        1        TFE        F         1       -0.25        
       3        F        1        TFE        F         1       -0.25        
       4        F        1        TFE        F         1       -0.25        
       5      CH2        1        TFE      CH2         1        0.25        
       6       OA        1        TFE       OA         1       -0.65        
       7       HO        1        TFE       HO         1        0.41        
    #else
    ; Use Charges from VanBuuren
       1        C        1        TFE        C         1        0.59        
       2        F        1        TFE        F         1       -0.2         
       3        F        1        TFE        F         1       -0.2         
       4        F        1        TFE        F         1       -0.2         
       5      CH2        1        TFE      CH2         1        0.26        
       6       OA        1        TFE       OA         1       -0.55        
       7       HO        1        TFE       HO         1        0.3         
    #endif

    [ bonds ]
    ;  ai    aj funct           c0           c1
        6     7     1 1.000000e-01 3.138000e+05 
        1     2     1 1.360000e-01 4.184000e+05 
        1     3     1 1.360000e-01 4.184000e+05 
        1     4     1 1.360000e-01 4.184000e+05 
        1     5     1 1.530000e-01 3.347000e+05 
        5     6     1 1.430000e-01 3.347000e+05 
    ...

This mechanism is used by pdb2gmx to implement optional position
restraints ([subsec:positionrestraint]) by #include-ing an .itp file
whose contents will be meaningful only if a particular #define is set
(and spelled correctly!)

Topologies for free energy calculations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Free energy differences between two systems, A and B, can be calculated
as described in sec. [sec:fecalc]. Systems A and B are described by
topologies consisting of the same number of molecules with the same
number of atoms. Masses and non-bonded interactions can be perturbed by
adding B parameters under the directive. Bonded interactions can be
perturbed by adding B parameters to the bonded types or the bonded
interactions. The parameters that can be perturbed are listed in Tables
[tab:topfile1] and [tab:topfile2]. The :math:`\lambda`-dependence of the
interactions is described in section sec. [sec:feia]. The bonded
parameters that are used (on the line of the bonded interaction
definition, or the ones looked up on atom types in the bonded type
lists) is explained in Table [tab:topfe]. In most cases, things should
work intuitively. When the A and B atom types in a bonded interaction
are not all identical and parameters are not present for the B-state,
either on the line or in the bonded types, grompp uses the A-state
parameters and issues a warning. For free energy calculations, all or no
parameters for topology B (:math:`\lambda = 1`) should be added on the
same line, after the normal parameters, in the same order as the normal
parameters. From GROMACS 4.6 onward, if :math:`\lambda` is treated as a
vector, then the bonded-lambdas component controls all bonded terms that
are not explicitly labeled as restraints. Restrain terms are controlled
by the restraint-lambdas component.

| Below is an example of a topology which changes from 200 propanols to
  200 pentanes using the GROMOS-96 force field.

::

     
    ; Include force field parameters
    #include "gromos43a1.ff/forcefield.itp"

    [ moleculetype ]
    ; Name            nrexcl
    PropPent          3

    [ atoms ]
    ; nr type resnr residue atom cgnr  charge    mass  typeB chargeB  massB
      1    H    1     PROP    PH    1   0.398    1.008  CH3     0.0  15.035
      2   OA    1     PROP    PO    1  -0.548  15.9994  CH2     0.0  14.027
      3  CH2    1     PROP   PC1    1   0.150   14.027  CH2     0.0  14.027
      4  CH2    1     PROP   PC2    2   0.000   14.027
      5  CH3    1     PROP   PC3    2   0.000   15.035

    [ bonds ]
    ;  ai    aj funct    par_A  par_B 
        1     2     2    gb_1   gb_26
        2     3     2    gb_17  gb_26
        3     4     2    gb_26  gb_26
        4     5     2    gb_26

    [ pairs ]
    ;  ai    aj funct
        1     4     1
        2     5     1

    [ angles ]
    ;  ai    aj    ak funct    par_A   par_B
        1     2     3     2    ga_11   ga_14
        2     3     4     2    ga_14   ga_14
        3     4     5     2    ga_14   ga_14

    [ dihedrals ]
    ;  ai    aj    ak    al funct    par_A   par_B
        1     2     3     4     1    gd_12   gd_17
        2     3     4     5     1    gd_17   gd_17

    [ system ]
    ; Name
    Propanol to Pentane

    [ molecules ]
    ; Compound        #mols
    PropPent          200

Atoms that are not perturbed, PC2 and PC3, do not need B-state parameter
specifications, since the B parameters will be copied from the A
parameters. Bonded interactions between atoms that are not perturbed do
not need B parameter specifications, as is the case for the last bond in
the example topology. Topologies using the OPLS/AA force field need no
bonded parameters at all, since both the A and B parameters are
determined by the atom types. Non-bonded interactions involving one or
two perturbed atoms use the free-energy perturbation functional forms.
Non-bonded interactions between two non-perturbed atoms use the normal
functional forms. This means that when, for instance, only the charge of
a particle is perturbed, its Lennard-Jones interactions will also be
affected when lambda is not equal to zero or one.

**Note** that this topology uses the GROMOS-96 force field, in which the
bonded interactions are not determined by the atom types. The bonded
interaction strings are converted by the C-preprocessor. The force-field
parameter files contain lines like:

::

    #define gb_26       0.1530  7.1500e+06

    #define gd_17     0.000       5.86          3

Constraint forces
~~~~~~~~~~~~~~~~~

| The constraint force between two atoms in one molecule can be
  calculated with the free energy perturbation code by adding a
  constraint between the two atoms, with a different length in the A and
  B topology. When the B length is 1 nm longer than the A length and
  lambda is kept constant at zero, the derivative of the Hamiltonian
  with respect to lambda is the constraint force. For constraints
  between molecules, the pull code can be used, see sec. [sec:pull].
  Below is an example for calculating the constraint force at 0.7 nm
  between two methanes in water, by combining the two methanes into one
  “molecule.” **Note** that the definition of a “molecule” in GROMACS
  does not necessarily correspond to the chemical definition of a
  molecule. In GROMACS, a “molecule” can be defined as any group of
  atoms that one wishes to consider simultaneously. The added constraint
  is of function type 2, which means that it is not used for generating
  exclusions (see sec. [sec:excl]). Note that the constraint free energy
  term is included in the derivative term, and is specifically included
  in the bonded-lambdas component. However, the free energy for changing
  constraints is *not* included in the potential energy differences used
  for BAR and MBAR, as this requires reevaluating the energy at each of
  the constraint components. This functionality is planned for later
  versions.

::

    ; Include force-field parameters
    #include "gromos43a1.ff/forcefield.itp"

    [ moleculetype ]
    ; Name            nrexcl
    Methanes               1

    [ atoms ]
    ; nr   type   resnr  residu   atom    cgnr     charge    mass
       1    CH4     1     CH4      C1       1          0    16.043
       2    CH4     1     CH4      C2       2          0    16.043
    [ constraints ]
    ;  ai    aj funct   length_A  length_B
        1     2     2        0.7       1.7

    #include "gromos43a1.ff/spc.itp"

    [ system ]
    ; Name
    Methanes in Water

    [ molecules ]
    ; Compound        #mols
    Methanes              1
    SOL                2002

Coordinate file
~~~~~~~~~~~~~~~

Files with the .gro file extension contain a molecular structure in
GROMOS-87 format. A sample piece is included below:

::

    MD of 2 waters, reformat step, PA aug-91
        6
        1WATER  OW1    1   0.126   1.624   1.679  0.1227 -0.0580  0.0434
        1WATER  HW2    2   0.190   1.661   1.747  0.8085  0.3191 -0.7791
        1WATER  HW3    3   0.177   1.568   1.613 -0.9045 -2.6469  1.3180
        2WATER  OW1    4   1.275   0.053   0.622  0.2519  0.3140 -0.1734
        2WATER  HW2    5   1.337   0.002   0.680 -1.0641 -1.1349  0.0257
        2WATER  HW3    6   1.326   0.120   0.568  1.9427 -0.8216 -0.0244
       1.82060   1.82060   1.82060

This format is fixed, *i.e.* all columns are in a fixed position. If you
want to read such a file in your own program without using the GROMACS
libraries you can use the following formats:

**C-format:** “%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f”

Or to be more precise, with title *etc.* it looks like this:

::

      "%s\n", Title
      "%5d\n", natoms
      for (i=0; (i<natoms); i++) {
        "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n",
          residuenr,residuename,atomname,atomnr,x,y,z,vx,vy,vz
      }
      "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n",
        box[X][X],box[Y][Y],box[Z][Z],
        box[X][Y],box[X][Z],box[Y][X],box[Y][Z],box[Z][X],box[Z][Y]

**Fortran format:** (i5,2a5,i5,3f8.3,3f8.4)

So confin.gro is the GROMACS coordinate file and is almost the same as
the GROMOS-87 file (for GROMOS users: when used with ntx=7). The only
difference is the box for which GROMACS uses a tensor, not a vector.

Force field organization 
-------------------------

Force-field files
~~~~~~~~~~~~~~~~~

Many force fields are available by default. Force fields are detected by
the presence of <name>.ff directories in the $GMXLIB/share/gromacs/top
sub-directory and/or the working directory. The information regarding
the location of the force field files is printed by pdb2gmx so you can
easily keep track of which version of a force field is being called, in
case you have made modifications in one location or another. The force
fields included with GROMACS are:

-  AMBER03 protein, nucleic AMBER94 (Duan et al., J. Comp. Chem. 24,
   1999-2012, 2003)

-  AMBER94 force field (Cornell et al., JACS 117, 5179-5197, 1995)

-  AMBER96 protein, nucleic AMBER94 (Kollman et al., Acc. Chem. Res. 29,
   461-469, 1996)

-  AMBER99 protein, nucleic AMBER94 (Wang et al., J. Comp. Chem. 21,
   1049-1074, 2000)

-  AMBER99SB protein, nucleic AMBER94 (Hornak et al., Proteins 65,
   712-725, 2006)

-  AMBER99SB-ILDN protein, nucleic AMBER94 (Lindorff-Larsen et al.,
   Proteins 78, 1950-58, 2010)

-  AMBERGS force field (Garcia & Sanbonmatsu, PNAS 99, 2782-2787, 2002)

-  CHARMM27 all-atom force field (CHARM22 plus CMAP for proteins)

-  GROMOS96 43a1 force field

-  GROMOS96 43a2 force field (improved alkane dihedrals)

-  GROMOS96 45a3 force field (Schuler JCC 2001 22 1205)

-  GROMOS96 53a5 force field (JCC 2004 vol 25 pag 1656)

-  GROMOS96 53a6 force field (JCC 2004 vol 25 pag 1656)

-  GROMOS96 54a7 force field (Eur. Biophys. J. (2011), 40,, 843-856,
   DOI: 10.1007/s00249-011-0700-9)

-  OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)

A force field is included at the beginning of a topology file with an
#include statement followed by <name>.ff/forcefield.itp. This statement
includes the force-field file, which, in turn, may include other
force-field files. All the force fields are organized in the same way.
An example of the amber99.ff/forcefield.itp was shown in
[subsec:topfile].

For each force field, there several files which are only used by
pdb2gmx. These are: residue databases (.rtp, see [subsec:rtp]) the
hydrogen database (.hdb, see [subsec:hdb]), two termini databases
(.n.tdb and .c.tdb, see [subsec:tdb]) and the atom type database (.atp,
see [subsec:atomtype]), which contains only the masses. Other optional
files are described in sec. [sec:pdb2gmxfiles].

Changing force-field parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If one wants to change the parameters of few bonded interactions in a
molecule, this is most easily accomplished by typing the parameters
behind the definition of the bonded interaction directly in the .top
file under the section (see [subsec:topfile] for the format and units).
If one wants to change the parameters for all instances of a certain
interaction one can change them in the force-field file or add a new
section after including the force field. When parameters for a certain
interaction are defined multiple times, the last definition is used. As
of GROMACS version 3.1.3, a warning is generated when parameters are
redefined with a different value. Changing the Lennard-Jones parameters
of an atom type is not recommended, because in the GROMOS force fields
the Lennard-Jones parameters for several combinations of atom types are
not generated according to the standard combination rules. Such
combinations (and possibly others that do follow the combination rules)
are defined in the section, and changing the Lennard-Jones parameters of
an atom type has no effect on these combinations.

Adding atom types
~~~~~~~~~~~~~~~~~

As of GROMACS version 3.1.3, atom types can be added in an extra section
after the the inclusion of the normal force field. After the definition
of the new atom type(s), additional non-bonded and pair parameters can
be defined. In pre-3.1.3 versions of GROMACS, the new atom types needed
to be added in the section of the force-field files, because all
non-bonded parameters above the last section would be overwritten using
the standard combination rules.

Special Topics
==============

Free energy implementation
--------------------------

For free energy calculations, there are two things that must be
specified; the end states, and the pathway connecting the end states.
The end states can be specified in two ways. The most straightforward is
through the specification of end states in the topology file. Most
potential forms support both an :math:`A` state and a :math:`B` state.
Whenever both states are specified, then the :math:`A` state corresponds
to the initial free energy state, and the :math:`B` state corresponds to
the final state.

In some cases, the end state can also be defined in some cases without
altering the topology, solely through the .mdp file, through the use of
the couple-moltype,couple-lambda0, couple-lambda1, and couple-intramol
mdp keywords. Any molecule type selected in couple-moltype will
automatically have a :math:`B` state implicitly constructed (and the
:math:`A` state redefined) according to the couple-lambda keywords.
couple-lambda0 and couple-lambda1 define the non-bonded parameters that
are present in the :math:`A` state (couple-lambda0) and the :math:`B`
state (couple-lambda1). The choices are ’q’,’vdw’, and ’vdw-q’; these
indicate the Coulombic, van der Waals, or both parameters that are
turned on in the respective state.

Once the end states are defined, then the path between the end states
has to be defined. This path is defined solely in the .mdp file.
Starting in 4.6, :math:`\lambda` is a vector of components, with
Coulombic, van der Waals, bonded, restraint, and mass components all
able to be adjusted independently. This makes it possible to turn off
the Coulombic term linearly, and then the van der Waals using soft core,
all in the same simulation. This is especially useful for replica
exchange or expanded ensemble simulations, where it is important to
sample all the way from interacting to non-interacting states in the
same simulation to improve sampling.

fep-lambdas is the default array of :math:`\lambda` values ranging from
0 to 1. All of the other lambda arrays use the values in this array if
they are not specified. The previous behavior, where the pathway is
controlled by a single :math:`\lambda` variable, can be preserved by
using only fep-lambdas to define the pathway.

For example, if you wanted to first to change the Coulombic terms, then
the van der Waals terms, changing bonded at the same time rate as the
van der Waals, but changing the restraints throughout the first
two-thirds of the simulation, then you could use this :math:`\lambda`
vector:

::

    coul-lambdas           = 0.0 0.2 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    vdw-lambdas            = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
    bonded-lambdas         = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
    restraint-lambdas      = 0.0 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0

This is also equivalent to:

::

    fep-lambdas            = 0.0 0.0 0.0 0.0 0.4 0.5 0.6 0.7 0.8 1.0
    coul-lambdas           = 0.0 0.2 0.5 1.0 1.0 1.0 1.0 1.0 1.0 1.0
    restraint-lambdas      = 0.0 0.0 0.1 0.2 0.3 0.5 0.7 1.0 1.0 1.0

The fep-lambda array, in this case, is being used as the default to fill
in the bonded and van der Waals :math:`\lambda` arrays. Usually, it’s
best to fill in all arrays explicitly, just to make sure things are
properly assigned.

If you want to turn on only restraints going from :math:`A` to
:math:`B`, then it would be:

::

    restraint-lambdas      = 0.0 0.1 0.2 0.4 0.6 1.0

and all of the other components of the :math:`\lambda` vector would be
left in the :math:`A` state.

To compute free energies with a vector :math:`\lambda` using
thermodynamic integration, then the TI equation becomes vector equation:

.. math:: \Delta F = \int \langle \nabla H \rangle \cdot d\vec{\lambda}

 or for finite differences:

.. math:: \Delta F \approx \int \sum \langle \nabla H \rangle \cdot \Delta\lambda

The external pymbar script downloaded from https://SimTK.org/home/pymbar
can compute this integral automatically from the GROMACS dhdl.xvg
output.

Potential of mean force
-----------------------

A potential of mean force (PMF) is a potential that is obtained by
integrating the mean force from an ensemble of configurations. In
GROMACS, there are several different methods to calculate the mean
force. Each method has its limitations, which are listed below.

-  **pull code:** between the centers of mass of molecules or groups of
   molecules.

-  **AWH code:** currently acts on coordinates provided by the pull
   code.

-  **free-energy code with harmonic bonds or constraints:** between
   single atoms.

-  **free-energy code with position restraints:** changing the
   conformation of a relatively immobile group of atoms.

-  **pull code in limited cases:** between groups of atoms that are part
   of a larger molecule for which the bonds are constrained with SHAKE
   or LINCS. If the pull group if relatively large, the pull code can be
   used.

The pull and free-energy code a described in more detail in the
following two sections.

Entropic effects
^^^^^^^^^^^^^^^^

When a distance between two atoms or the centers of mass of two groups
is constrained or restrained, there will be a purely entropic
contribution to the PMF due to the rotation of the two
groups \ `134 <#ref-RMNeumann1980a>`__. For a system of two
non-interacting masses the potential of mean force is:

.. math:: V_{pmf}(r) = -(n_c - 1) k_B T \log(r)

 where :math:`n_c` is the number of dimensions in which the constraint
works (i.e. :math:`n_c=3` for a normal constraint and :math:`n_c=1` when
only the :math:`z`-direction is constrained). Whether one needs to
correct for this contribution depends on what the PMF should represent.
When one wants to pull a substrate into a protein, this entropic term
indeed contributes to the work to get the substrate into the protein.
But when calculating a PMF between two solutes in a solvent, for the
purpose of simulating without solvent, the entropic contribution should
be removed. **Note** that this term can be significant; when at 300K the
distance is halved, the contribution is 3.5 kJ mol\ :math:`^{-1}`.

Non-equilibrium pulling
-----------------------

When the distance between two groups is changed continuously, work is
applied to the system, which means that the system is no longer in
equilibrium. Although in the limit of very slow pulling the system is
again in equilibrium, for many systems this limit is not reachable
within reasonable computational time. However, one can use the Jarzynski
relation \ `135 <#ref-Jarzynski1997a>`__ to obtain the equilibrium
free-energy difference :math:`\Delta G` between two distances from many
non-equilibrium simulations:

.. math::

   \Delta G_{AB} = -k_BT \log \left\langle e^{-\beta W_{AB}} \right\rangle_A
      \label{eq:Jarz}

 where :math:`W_{AB}` is the work performed to force the system along
one path from state A to B, the angular bracket denotes averaging over a
canonical ensemble of the initial state A and :math:`\beta=1/k_B T`.

The pull code
-------------

[sec:pull] The pull code applies forces or constraints between the
centers of mass of one or more pairs of groups of atoms. Each pull
reaction coordinate is called a “coordinate” and it operates on usually
two, but sometimes more, pull groups. A pull group can be part of one or
more pull coordinates. Furthermore, a coordinate can also operate on a
single group and an absolute reference position in space. The distance
between a pair of groups can be determined in 1, 2 or 3 dimensions, or
can be along a user-defined vector. The reference distance can be
constant or can change linearly with time. Normally all atoms are
weighted by their mass, but an additional weighting factor can also be
used.

.. figure:: plots/pull
   :alt: Schematic picture of pulling a lipid out of a lipid bilayer
   with umbrella pulling. :math:`V_{rup}` is the velocity at which the
   spring is retracted, :math:`Z_{link}` is the atom to which the spring
   is attached and :math:`Z_{spring}` is the location of the spring.
   :width: 6.00000cm

   Schematic picture of pulling a lipid out of a lipid bilayer with
   umbrella pulling. :math:`V_{rup}` is the velocity at which the spring
   is retracted, :math:`Z_{link}` is the atom to which the spring is
   attached and :math:`Z_{spring}` is the location of the spring.

Several different pull types, i.e. ways to apply the pull force, are
supported, and in all cases the reference distance can be constant or
linearly changing with time.

#. **Umbrella pulling** A harmonic potential is applied between the
   centers of mass of two groups. Thus, the force is proportional to the
   displacement.

#. **Constraint pulling** The distance between the centers of mass of
   two groups is constrained. The constraint force can be written to a
   file. This method uses the SHAKE algorithm but only needs 1 iteration
   to be exact if only two groups are constrained.

#. **Constant force pulling** A constant force is applied between the
   centers of mass of two groups. Thus, the potential is linear. In this
   case there is no reference distance of pull rate.

#. **Flat bottom pulling** Like umbrella pulling, but the potential and
   force are zero for coordinate values below (pull-coord?-type =
   flat-bottom) or above (pull-coord?-type = flat-bottom-high) a
   reference value. This is useful for restraining e.g. the distance
   between two molecules to a certain region.

In addition, there are different types of reaction coordinates,
so-called pull geometries. These are set with the .mdp option
pull-coord?-geometry.

Definition of the center of mass
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In GROMACS, there are three ways to define the center of mass of a
group. The standard way is a “plain” center of mass, possibly with
additional weighting factors. With periodic boundary conditions it is no
longer possible to uniquely define the center of mass of a group of
atoms. Therefore, a reference atom is used. For determining the center
of mass, for all other atoms in the group, the closest periodic image to
the reference atom is used. This uniquely defines the center of mass. By
default, the middle (determined by the order in the topology) atom is
used as a reference atom, but the user can also select any other atom if
it would be closer to center of the group.

For a layered system, for instance a lipid bilayer, it may be of
interest to calculate the PMF of a lipid as function of its distance
from the whole bilayer. The whole bilayer can be taken as reference
group in that case, but it might also be of interest to define the
reaction coordinate for the PMF more locally. The .mdp option
pull-coord?-geometry = cylinder does not use all the atoms of the
reference group, but instead dynamically only those within a cylinder
with radius pull-cylinder-r around the pull vector going through the
pull group. This only works for distances defined in one dimension, and
the cylinder is oriented with its long axis along this one dimension. To
avoid jumps in the pull force, contributions of atoms are weighted as a
function of distance (in addition to the mass weighting):

.. math::

   \begin{aligned}
   w(r < r_\mathrm{cyl}) & = &
   1-2 \left(\frac{r}{r_\mathrm{cyl}}\right)^2 + \left(\frac{r}{r_\mathrm{cyl}}\right)^4 \\
   w(r \geq r_\mathrm{cyl}) & = & 0\end{aligned}

 Note that the radial dependence on the weight causes a radial force on
both cylinder group and the other pull group. This is an undesirable,
but unavoidable effect. To minimize this effect, the cylinder radius
should be chosen sufficiently large. The effective mass is 0.47 times
that of a cylinder with uniform weights and equal to the mass of uniform
cylinder of 0.79 times the radius.

.. figure:: plots/pullref
   :alt: Comparison of a plain center of mass reference group versus a
   cylinder reference group applied to interface systems. C is the
   reference group. The circles represent the center of mass of two
   groups plus the reference group, :math:`d_c` is the reference
   distance.
   :width: 6.00000cm

   Comparison of a plain center of mass reference group versus a
   cylinder reference group applied to interface systems. C is the
   reference group. The circles represent the center of mass of two
   groups plus the reference group, :math:`d_c` is the reference
   distance.

For a group of molecules in a periodic system, a plain reference group
might not be well-defined. An example is a water slab that is connected
periodically in :math:`x` and :math:`y`, but has two liquid-vapor
interfaces along :math:`z`. In such a setup, water molecules can
evaporate from the liquid and they will move through the vapor, through
the periodic boundary, to the other interface. Such a system is
inherently periodic and there is no proper way of defining a “plain”
center of mass along :math:`z`. A proper solution is to using a cosine
shaped weighting profile for all atoms in the reference group. The
profile is a cosine with a single period in the unit cell. Its phase is
optimized to give the maximum sum of weights, including mass weighting.
This provides a unique and continuous reference position that is nearly
identical to the plain center of mass position in case all atoms are all
within a half of the unit-cell length. See ref `136 <#ref-Engin2010a>`__
for details.

When relative weights :math:`w_i` are used during the calculations,
either by supplying weights in the input or due to cylinder geometry or
due to cosine weighting, the weights need to be scaled to conserve
momentum:

.. math::

   w'_i = w_i
   \left. \sum_{j=1}^N w_j \, m_j \right/ \sum_{j=1}^N w_j^2 \, m_j

 where :math:`m_j` is the mass of atom :math:`j` of the group. The mass
of the group, required for calculating the constraint force, is:

.. math:: M = \sum_{i=1}^N w'_i \, m_i

 The definition of the weighted center of mass is:

.. math:: {\mbox{\boldmath ${r}$}}_{com} = \left. \sum_{i=1}^N w'_i \, m_i \, {\mbox{\boldmath ${r}$}}_i \right/ M

 From the centers of mass the AFM, constraint, or umbrella force
:math:`{\mbox{\boldmath ${F}$}}_{\!com}` on each group can be
calculated. The force on the center of mass of a group is redistributed
to the atoms as follows:

.. math:: {\mbox{\boldmath ${F}$}}_{\!i} = \frac{w'_i \, m_i}{M} \, {\mbox{\boldmath ${F}$}}_{\!com}

Definition of the pull direction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most common setup is to pull along the direction of the vector
containing the two pull groups, this is selected with
pull-coord?-geometry = distance. You might want to pull along a certain
vector instead, which is selected with pull-coord?-geometry = direction.
But this can cause unwanted torque forces in the system, unless you pull
against a reference group with (nearly) fixed orientation, e.g. a
membrane protein embedded in a membrane along x/y while pulling along z.
If your reference group does not have a fixed orientation, you should
probably use pull-coord?-geometry = direction-relative, see
Fig. [fig:pulldirrel]. Since the potential now depends on the
coordinates of two additional groups defining the orientation, the
torque forces will work on these two groups.

.. figure:: plots/pulldirrel
   :alt: The pull setup for geometry direction-relative. The “normal”
   pull groups are 1 and 2. Groups 3 and 4 define the pull direction and
   thus the direction of the normal pull forces (red). This leads to
   reaction forces (blue) on groups 3 and 4, which are perpendicular to
   the pull direction. Their magnitude is given by the “normal” pull
   force times the ratio of :math:`d_p` and the distance between groups
   3 and 4.
   :width: 5.00000cm

   The pull setup for geometry direction-relative. The “normal” pull
   groups are 1 and 2. Groups 3 and 4 define the pull direction and thus
   the direction of the normal pull forces (red). This leads to reaction
   forces (blue) on groups 3 and 4, which are perpendicular to the pull
   direction. Their magnitude is given by the “normal” pull force times
   the ratio of :math:`d_p` and the distance between groups 3 and 4.

Definition of the angle and dihedral pull geometries
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Four pull groups are required for pull-coord?-geometry = angle. In the
same way as for geometries with two groups, each consecutive pair of
groups :math:`i` and :math:`i+1` define a vector connecting the COMs of
groups :math:`i` and :math:`i+1`. The angle is defined as the angle
between the two resulting vectors. E.g., the .mdp option
pull-coord?-groups = 1 2 2 4 defines the angle between the vector from
the COM of group 1 to the COM of group 2 and the vector from the COM of
group 2 to the COM of group 4. The angle takes values in the closed
interval [0, 180] deg. For pull-coord?-geometry = angle-axis the angle
is defined with respect to a reference axis given by pull-coord?-vec and
only two groups need to be given. The dihedral geometry requires six
pull groups. These pair up in the same way as described above and so
define three vectors. The dihedral angle is defined as the angle between
the two planes spanned by the two first and the two last vectors.
Equivalently, the dihedral angle can be seen as the angle between the
first and the third vector when these vectors are projected onto a plane
normal to the second vector (the axis vector). As an example, consider a
dihedral angle involving four groups: 1, 5, 8 and 9. Here, the .mdp
option pull-coord?-groups = 8 1 1 5 5 9 specifies the three vectors that
define the dihedral angle: the first vector is the COM distance vector
from group 8 to 1, the second vector is the COM distance vector from
group 1 to 5, and the third vector is the COM distance vector from group
5 to 9. The dihedral angle takes values in the interval (-180, 180] deg
and has periodic boundaries.

Limitations
^^^^^^^^^^^

There is one theoretical limitation: strictly speaking, constraint
forces can only be calculated between groups that are not connected by
constraints to the rest of the system. If a group contains part of a
molecule of which the bond lengths are constrained, the pull constraint
and LINCS or SHAKE bond constraint algorithms should be iterated
simultaneously. This is not done in GROMACS. This means that for
simulations with constraints = all-bonds in the .mdp file pulling is,
strictly speaking, limited to whole molecules or groups of molecules. In
some cases this limitation can be avoided by using the free energy code,
see sec. [sec:fepmf]. In practice, the errors caused by not iterating
the two constraint algorithms can be negligible when the pull group
consists of a large amount of atoms and/or the pull force is small. In
such cases, the constraint correction displacement of the pull group is
small compared to the bond lengths.

Adaptive biasing with AWH
-------------------------

[sec:awh] The accelerated weight histogram method
(AWH) `137 <#ref-lindahl2014accelerated>`__ calculates the PMF along a
reaction coordinate by adding an adaptively determined biasing
potential. AWH flattens free energy barriers along the reaction
coordinate by applying a history-dependent potential to the system that
“fills up” free energy minima. This is similar in spirit to other
adaptive biasing potential methods, e.g. the
Wang-Landau \ `138 <#ref-wang2001efficient>`__, local
elevation \ `139 <#ref-huber1994local>`__ and
metadynamics \ `140 <#ref-laio2002escaping>`__ methods. The initial
sampling stage of AWH makes the method robust against the choice of
input parameters. Furthermore, the target distribution along the
reaction coordinate may be chosen freely.

Basics of the method
~~~~~~~~~~~~~~~~~~~~

Rather than biasing the reaction coordinate :math:`\xi(x)` directly, AWH
acts on a *reference coordinate* :math:`\lambda`. The reaction
coordinate :math:`\xi(x)` is coupled to :math:`\lambda` with a harmonic
potential

.. math:: Q(\xi,\lambda) = \frac{1}{2} \beta k (\xi - \lambda)^2,

 so that for large force constants :math:`k`,
:math:`\xi \approx \lambda`. Note the use of dimensionless energies for
compatibility with previously published work. Units of energy are
obtained by multiplication with :math:`k_BT=1/\beta`. In the simulation,
:math:`\lambda` samples the user-defined sampling interval :math:`I`.
For a multidimensional reaction coordinate :math:`\xi`, the sampling
interval is the Cartesian product :math:`I=\Pi_\mu I_\mu` (a rectangular
domain). The connection between atom coordinates and :math:`\lambda` is
established through the extended
ensemble \ `68 <#ref-Lyubartsev1992>`__,

.. math::

   \label{eq:awh:pxlambda}
   P(x,\lambda) = \frac{1}{\mathcal{Z}}e^{g(\lambda) - Q(\xi(x),\lambda) - V(x)},

 where :math:`g(\lambda)` is a bias function (a free variable) and
:math:`V(x)` is the unbiased potential energy of the system. The
distribution along :math:`\lambda` can be tuned to be any predefined
*target distribution* :math:`\rho(\lambda)` (often chosen to be flat) by
choosing :math:`g(\lambda)` wisely. This is evident from

.. math::

   \label{eq:awh:plambda}
   P(\lambda) = \int P(x,\lambda)  dx = 
   \frac{1}{\mathcal{Z}}e^{g(\lambda)} \int e^{- Q(\xi(x),\lambda) - V(x)}  dx 
   \equiv \frac{1}{\mathcal{Z}}e^{g(\lambda) - F(\lambda)},

 where :math:`F(\lambda)` is the free energy

.. math::

   \label{eq:awh:flambda}
   F(\lambda) = -\ln \int e^{- Q(\xi(x),\lambda) - V(x)}  dx.

 Being the convolution of the PMF with the Gaussian defined by the
harmonic potential, :math:`F(\lambda)` is a smoothened version of the
PMF. Eq. [eq:awh:plambda] shows that in order to obtain
:math:`P(\lambda)=\rho(\lambda)`, :math:`F(\lambda)` needs to be
determined accurately. Thus, AWH adaptively calculates
:math:`F(\lambda)` and simultaneously converges :math:`P(\lambda)`
toward :math:`\rho(\lambda)`.

The free energy update
^^^^^^^^^^^^^^^^^^^^^^

AWH is initialized with an estimate of the free energy
:math:`F_0(\lambda)`. At regular time intervals this estimate is updated
using data collected in between the updates. At update :math:`n`, the
applied bias :math:`g_n(\lambda)` is a function of the current free
energy estimate :math:`F_n(\lambda)` and target distribution
:math:`\rho_n(\lambda)`,

.. math::

   \label{eq:awh:grhofrelation}
   g_n(\lambda) = \ln \rho_n(\lambda) +F_n(\lambda),

 which is consistent with Eq. [eq:awh:plambda]. Note that also the
target distribution may be updated during the simulation (see examples
in section [sec:awh:targets]). Substituting this choice of :math:`g=g_n`
back into Eq. [eq:awh:plambda] yields the simple free energy update

.. math::

   \label{eq:awh:dfnaive}
   \Delta F_n(\lambda) 
   = F(\lambda) - F_n(\lambda) 
   = -\ln\frac{P_n(\lambda)}{\rho_n(\lambda)},

 which would yield a better estimate :math:`F_{n+1} = F_n + \Delta F_n`,
assuming :math:`P_n(\lambda)` can be measured accurately. AWH estimates
:math:`P_n(\lambda)` by regularly calculating the conditional
distribution

.. math::

   \label{eq:awh:omega}
   \omega_n(\lambda|x) \equiv P_n(\lambda|x) = \frac{e^{g_n(\lambda) - Q(\xi(x), \lambda)}}{\sum_{\lambda'} e^{g_n(\lambda') - Q(\xi(x),\lambda')}}.

 Accumulating these probability weights yields
:math:`\sum_t \omega(\lambda|x(t)) \sim P_n(\lambda)`, where
:math:`\int P_n(\lambda|x) P_n(x) dx = P_n(\lambda)` has been used. The
:math:`\omega_n(\lambda|x)` weights are thus the samples of the AWH
method. With the limited amount of sampling one has in practice, update
scheme [eq:awh:dfnaive] yields very noisy results. AWH instead applies a
free energy update that has the same form but which can be applied
repeatedly with limited and localized sampling,

.. math:: \Delta F_n = -\ln \frac{W_n(\lambda) + \sum_t \omega_n(\lambda|x(t))}{W_n(\lambda) + \sum_t\rho_n(\lambda)) }.

 Here :math:`W_n(\lambda)` is the *reference weight histogram*
representing prior sampling. The update for :math:`W(\lambda)`,
disregarding the initial stage (see section [sec:awh:initial-stage]), is

.. math::

   \label{eq:awh:w-update}
   W_{n+1}(\lambda) = W_n(\lambda) + \sum_t\rho_n(\lambda).

 Thus, the weight histogram equals the targeted, “ideal” history of
samples. There are two important things to note about the free energy
update. First, sampling is driven away from oversampled, currently local
regions. For such :math:`\lambda` values,
:math:`\omega_n(\lambda) > \rho_n(\lambda)` and
:math:`\Delta F_n(\lambda) < 0`, which by Eq. [eq:awh:grhofrelation]
implies :math:`\Delta g_n(\lambda) < 0` (assuming
:math:`\Delta \rho_n \equiv 0`). Thus, the probability to sample
:math:`\lambda` decreases after the update (see Eq. [eq:awh:plambda]).
Secondly, the normalization of the histogram
:math:`N_n=\sum_\lambda W_n(\lambda)`, determines the update size
:math:`|\Delta F(\lambda)|`. For instance, for a single sample
:math:`\omega(\lambda|x)`, the shape of the update is approximately a
Gaussian function of width :math:`\sigma=1/\sqrt{\beta k}` and height
:math:`\propto 1/N_n` `137 <#ref-lindahl2014accelerated>`__,

.. math::

   \label{eq:awh:dfsize}
   |\Delta F_n(\lambda)| \propto \frac{1}{N_n} e^{-\frac{1}{2} \beta k (\xi(x) - \lambda)^2}.

 Therefore, as samples accumulate in :math:`W(\lambda)` and :math:`N_n`
grows, the updates get smaller, allowing for the free energy to
converge.

Note that quantity of interest to the user is not :math:`F(\lambda)` but
the PMF :math:`\Phi(\xi)`. :math:`\Phi(\xi)` is extracted by reweighting
samples :math:`\xi(t)` on the
fly \ `137 <#ref-lindahl2014accelerated>`__ (see also
section [sec:awh:reweight]) and will converge at the same rate as
:math:`F(\lambda)`, see Fig. [fig:awh:bias-evolution]. The PMF will be
written to output (see section [sec:awh:usage]).

Applying the bias to the system
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The bias potential can be applied to the system in two ways. Either by
applying a harmonic potential centered at :math:`\lambda(t)`, which is
sampled using (rejection-free) Monte-Carlo sampling from the conditional
distribution :math:`\omega_n(\lambda|x(t)) = P_n(\lambda|x(t))`, see
Eq. [eq:awh:omega]. This is also called Gibbs sampling or independence
sampling. Alternatively, and by default in the code, the following
*convolved bias potential* can be applied,

.. math::

   \label{eq:awh:biaspotential}
   U_n(\xi) = -\ln \int e^{ g_n(\lambda) -Q(\xi,\lambda)} d \lambda.

 These two approaches are equivalent in the sense that they give rise to
the same biased probabilities :math:`P_n(x)` (cf. [eq:awh:pxlambda])
while the dynamics are clearly different in the two cases. This choice
does not affect the internals of the AWH algorithm, only what force and
potential AWH returns to the MD engine.

|AWH evolution in time for a Brownian particle in a double-well
potential. The reaction coordinate :math:`\xi(t)` traverses the sampling
interval multiple times in the initial stage before exiting and entering
the final stage (top left). In the final stage, the dynamics of
:math:`\xi` becomes increasingly diffusive. The times of covering are
shown as :math:`\times`-markers of different colors. At these times the
free energy update size :math:`\sim 1/N`, where :math:`N` is the size of
the weight histogram, is decreased by scaling :math:`N` by a factor of
:math:`\gamma=3` (top right). In the final stage, :math:`N` grows at the
sampling rate and thus :math:`1/N\sim1/t`. The exit from the final stage
is determined on the fly by ensuring that the effective sample weight
:math:`s` of data collected in the final stage exceeds that of initial
stage data (bottom left; note that :math:`\ln s(t)` is plotted). An
estimate of the PMF is also extracted from the simulation (bottom
right), which after exiting the initial stage should estimate global
free energy differences fairly accurately. | |AWH evolution in time for
a Brownian particle in a double-well potential. The reaction coordinate
:math:`\xi(t)` traverses the sampling interval multiple times in the
initial stage before exiting and entering the final stage (top left). In
the final stage, the dynamics of :math:`\xi` becomes increasingly
diffusive. The times of covering are shown as :math:`\times`-markers of
different colors. At these times the free energy update size
:math:`\sim 1/N`, where :math:`N` is the size of the weight histogram,
is decreased by scaling :math:`N` by a factor of :math:`\gamma=3` (top
right). In the final stage, :math:`N` grows at the sampling rate and
thus :math:`1/N\sim1/t`. The exit from the final stage is determined on
the fly by ensuring that the effective sample weight :math:`s` of data
collected in the final stage exceeds that of initial stage data (bottom
left; note that :math:`\ln s(t)` is plotted). An estimate of the PMF is
also extracted from the simulation (bottom right), which after exiting
the initial stage should estimate global free energy differences fairly
accurately. |

|AWH evolution in time for a Brownian particle in a double-well
potential. The reaction coordinate :math:`\xi(t)` traverses the sampling
interval multiple times in the initial stage before exiting and entering
the final stage (top left). In the final stage, the dynamics of
:math:`\xi` becomes increasingly diffusive. The times of covering are
shown as :math:`\times`-markers of different colors. At these times the
free energy update size :math:`\sim 1/N`, where :math:`N` is the size of
the weight histogram, is decreased by scaling :math:`N` by a factor of
:math:`\gamma=3` (top right). In the final stage, :math:`N` grows at the
sampling rate and thus :math:`1/N\sim1/t`. The exit from the final stage
is determined on the fly by ensuring that the effective sample weight
:math:`s` of data collected in the final stage exceeds that of initial
stage data (bottom left; note that :math:`\ln s(t)` is plotted). An
estimate of the PMF is also extracted from the simulation (bottom
right), which after exiting the initial stage should estimate global
free energy differences fairly accurately. | |AWH evolution in time for
a Brownian particle in a double-well potential. The reaction coordinate
:math:`\xi(t)` traverses the sampling interval multiple times in the
initial stage before exiting and entering the final stage (top left). In
the final stage, the dynamics of :math:`\xi` becomes increasingly
diffusive. The times of covering are shown as :math:`\times`-markers of
different colors. At these times the free energy update size
:math:`\sim 1/N`, where :math:`N` is the size of the weight histogram,
is decreased by scaling :math:`N` by a factor of :math:`\gamma=3` (top
right). In the final stage, :math:`N` grows at the sampling rate and
thus :math:`1/N\sim1/t`. The exit from the final stage is determined on
the fly by ensuring that the effective sample weight :math:`s` of data
collected in the final stage exceeds that of initial stage data (bottom
left; note that :math:`\ln s(t)` is plotted). An estimate of the PMF is
also extracted from the simulation (bottom right), which after exiting
the initial stage should estimate global free energy differences fairly
accurately. |

The initial stage
~~~~~~~~~~~~~~~~~

Initially, when the bias potential is far from optimal, samples will be
highly correlated. In such cases, letting :math:`W(\lambda)` accumulate
samples as prescribed by Eq. [eq:awh:w-update], entails a too rapid
decay of the free energy update size. This motivates splitting the
simulation into an *initial stage* where the weight histogram grows
according to a more restrictive and robust protocol, and a *final stage*
where the the weight histogram grows linearly at the sampling rate
(Eq. [eq:awh:w-update]). The AWH initial stage takes inspiration from
the well-known Wang-Landau algorithm \ `138 <#ref-wang2001efficient>`__,
although there are differences in the details.

In the initial stage the update size is kept constant (by keeping
:math:`N_n` constant) until a transition across the sampling interval
has been detected, a “covering”. For the definition of a covering, see
Eq. [eq:awh:covering] below. After a covering has occurred, :math:`N_n`
is scaled up by a constant “growth factor” :math:`\gamma`, chosen
heuristically as :math:`\gamma=3`. Thus, in the initial stage
:math:`N_n` is set dynamically as :math:`N_{n} = \gamma^{m} N_0`, where
:math:`m` is the number of coverings. Since the update size scales as
:math:`1/N` ( Eq. [eq:awh:dfsize]) this leads to a close to exponential
decay of the update size in the initial stage, see
Fig. [fig:awh:bias-evolution].

The update size directly determines the rate of change of
:math:`F_n(\lambda)` and hence, from Eq. [eq:awh:grhofrelation], also
the rate of change of the bias funcion :math:`g_n(\lambda)` Thus
initially, when :math:`N_n` is kept small and updates large, the system
will be driven along the reaction coordinate by the constantly
fluctuating bias. If :math:`N_0` is set small enough, the first
transition will typically be fast because of the large update size and
will quickly give a first rough estimate of the free energy. The second
transition, using :math:`N_1=\gamma N_0` refines this estimate further.
Thus, rather than very carefully filling free energy minima using a
small initial update size, the sampling interval is sweeped
back-and-forth multiple times, using a wide range of update sizes, see
Fig. [fig:awh:bias-evolution]. This way, the initial stage also makes
AWH robust against the choice of :math:`N_0`.

The covering criterion
^^^^^^^^^^^^^^^^^^^^^^

In the general case of a multidimensional reaction coordinate
:math:`\lambda=(\lambda_\mu)`, the sampling interval :math:`I` is
considered covered when all dimensions have been covered. A dimension
:math:`d` is covered if all points :math:`\lambda_\mu` in the
one-dimensional sampling interval :math:`I_\mu` have been “visited”.
Finally, a point :math:`\lambda_\mu \in I_\mu` has been visited if there
is at least one point :math:`\lambda^*\in I` with
:math:`\lambda^*_\mu = \lambda_\mu` that since the last covering has
accumulated probability weight corresponding to the peak of a
multidimensional Gaussian distribution

.. math::

   \label{eq:awh:covering}
   \Delta W(\lambda^*)
   \ge w_{\mathrm{peak}}
    \equiv \prod_\mu \frac{\Delta \lambda_\mu}{\sqrt{2\pi}\sigma_k}.

 Here, :math:`\Delta \lambda_\mu` is the point spacing of the
discretized :math:`I_\mu` and :math:`\sigma_k=1/\sqrt{\beta k_\mu}`
(where :math:`k_\mu` is the force constant) is the Gaussian width.

Exit from the initial stage
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For longer times, when major free energy barriers have largely been
flattened by the converging bias potential, the histogram
:math:`W(\lambda)` should grow at the actual sampling rate and the
initial stage needs to be exited \ `141 <#ref-belardinelli2007fast>`__.
There are multiple reasonable (heuristic) ways of determining when this
transition should take place. One option is to postulate that the number
of samples in the weight histogram :math:`N_n` should never exceed the
actual number of collected samples, and exit the initial stage when this
condition breaks \ `137 <#ref-lindahl2014accelerated>`__. In the initial
stage, :math:`N` grows close to exponentially while the collected number
of samples grows linearly, so an exit will surely occur eventually. Here
we instead apply an exit criterion based on the observation that
“artifically” keeping :math:`N` constant while continuing to collect
samples corresponds to scaling down the relative weight of old samples
relative to new ones. Similarly, the subsequent scaling up of :math:`N`
by a factor :math:`\gamma` corresponds to scaling up the weight of old
data. Briefly, the exit criterion is devised such that the weight of a
sample collected *after* the initial stage is always larger or equal to
the weight of a sample collected *during* the initial stage, see
Fig. [fig:awh:bias-evolution]. This is consistent with scaling down
early, noisy data.

The initial stage exit criterion will now be described in detail. We
start out at the beginning of a covering stage, so that :math:`N` has
just been scaled by :math:`\gamma` and is now kept constant. Thus, the
first sample of this stage has the weight :math:`s= 1/\gamma` relative
to the last sample of the previous covering stage. We assume that
:math:`\Delta N` samples are collected and added to :math:`W` for each
update . To keep :math:`N` constant, :math:`W` needs to be scaled down
by a factor :math:`N/(N + \Delta N)` after every update. Equivalently,
this means that new data is scaled up relative to old data by the
inverse factor. Thus, after :math:`\Delta n` updates a new sample has
the relative weight
:math:`s=(1/\gamma) [(N_n + \Delta N)/N_n]^{\Delta n}`. Now assume
covering occurs at this time. To continue to the next covering stage,
:math:`N` should be scaled by :math:`\gamma`, which corresponds to again
multiplying :math:`s` by :math:`1/\gamma`. If at this point
:math:`s \ge \gamma`, then after rescaling :math:`s \ge 1`; i.e. overall
the relative weight of a new sample relative to an old sample is still
growing fast. If on the contrary :math:`s < \gamma`, and this defines
the exit from the initial stage, then the initial stage is over and from
now :math:`N` simply grows at the sampling rate (see
Eq. [eq:awh:w-update]). To really ensure that :math:`s\ge 1` holds
before exiting, so that samples after the exit have at least the sample
weight of older samples, the last covering stage is extended by a
sufficient number of updates.

Choice of target distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The target distribution :math:`\rho(\lambda)` is traditionally chosen to
be uniform

.. math:: \rho_{\mathrm{const}}(\lambda) = \mathrm{const.}

 This choice exactly flattens :math:`F(\lambda)` in user-defined
sampling interval :math:`I`. Generally,
:math:`\rho(\lambda)=0, \lambda\notin I`. In certain cases other choices
may be preferable. For instance, in the multidimensional case the
rectangular sampling interval is likely to contain regions of very high
free energy, e.g. where atoms are clashing. To exclude such regions,
:math:`\rho(\lambda)` can specified by the following function of the
free energy

.. math::

   \label{eq:awh:rhocut}
   \rho_{\mathrm{cut}}(\lambda) \propto \frac{1}{1+ e^{F(\lambda) - F_{\mathrm{cut}}}}, 

 where :math:`F_{\mathrm{cut}}` is a free energy cutoff (relative to
:math:`\min_\lambda F(\lambda)`). Thus, regions of the sampling interval
where :math:`F(\lambda) > F_{\mathrm{cut}}` will be exponentially
suppressed (in a smooth fashion). Alternatively, very high free energy
regions could be avoided while still flattening more moderate free
energy barriers by targeting a Boltzmann distribution corresponding to
scaling :math:`\beta=1/k_BT` by a factor :math:`0<s_\beta<1`,

.. math::

   \label{eq:awh:rhoboltz}
   \rho_{\mathrm{Boltz}}(\lambda) \propto e^{-s_\beta F(\lambda)}, 

 The parameter :math:`s_\beta` determines to what degree the free energy
landscape is flattened; the lower :math:`s_\beta`, the flatter. Note
that both :math:`\rho_{\mathrm{cut}}(\lambda)` and
:math:`\rho_{\mathrm{Boltz}}(\lambda)` depend on :math:`F(\lambda)`,
which needs to be substituted by the current best estimate
:math:`F_n(\lambda)`. Thus, the target distribution is also updated
(consistently with Eq. [eq:awh:grhofrelation]).

There is in fact an alternative approach to obtaining
:math:`\rho_{\mathrm{Boltz}}(\lambda)` as the limiting target
distribution in AWH, which is particular in the way the weight histogram
:math:`W(\lambda)` and the target distribution :math:`\rho` are updated
and coupled to each other. This yields an evolution of the bias
potential which is very similar to that of well-tempered
metadynamics \ `142 <#ref-barducci2008well>`__,
see \ `137 <#ref-lindahl2014accelerated>`__ for details. Because of the
popularity and success of well-tempered metadynamics, this is a special
case worth considering. In this case :math:`\rho` is a function of the
reference weight histogram

.. math:: \rho_{\mathrm{Boltz,loc}}(\lambda) \propto W(\lambda), 

 and the update of the weight histogram is modified (cf.
Eq. [eq:awh:w-update])

.. math:: W_{n+1}(\lambda) =  W_{n}(\lambda) + s_{\beta}\sum_t \omega(\lambda|x(t)).

 Thus, here the weight histogram equals the real history of samples, but
scaled by :math:`s_\beta`. This target distribution is called *local*
Boltzmann since :math:`W` is only modified locally, where sampling has
taken place. We see that when :math:`s_\beta \approx 0` the histogram
essentially does not grow and the size of the free energy update will
stay at a constant value (as in the original formulation of
metadynamics). Thus, the free energy estimate will not converge, but
continue to fluctuate around the correct value. This illustrates the
inherent coupling between the convergence and choice of target
distribution for this special choice of target. Furthermore note that
when using :math:`\rho=\rho_{\mathrm{Boltz,loc}}` there is no initial
stage (section [sec:awh:initial-stage]). The rescaling of the weight
histogram applied in the initial stage is a global operation, which is
incompatible :math:`\rho_{\mathrm{Boltz,loc}}` only depending locally on
the sampling history.

Lastly, the target distribution can be modulated by arbitrary
probability weights

.. math:: \rho(\lambda) = \rho_0(\lambda) w_{\mathrm{user}}(\lambda).

 where :math:`w_{\mathrm{user}}(\lambda)` is provided by user data and
in principle :math:`\rho_0(\lambda)` can be any of the target
distributions mentioned above.

Multiple independent or sharing biases
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Multiple independent bias potentials may be applied within one
simulation. This only makes sense if the biased coordinates
:math:`\xi^{(1)}`, :math:`\xi^{(2)}`, :math:`\ldots` evolve essentially
independently from one another. A typical example of this would be when
applying an independent bias to each monomer of a protein. Furthermore,
multiple AWH simulations can be launched in parallel, each with a (set
of) indepedendent biases.

If the defined sampling interval is large relative to the diffusion time
of the reaction coordinate, traversing the sampling interval multiple
times as is required by the initial stage
(section [sec:awh:initial-stage]) may take an infeasible mount of
simulation time. In these cases it could be advantageous to parallelize
the work and have a group of multiple “walkers” :math:`\xi^{(i)}(t)`
share a single bias potential. This can be achieved by collecting
samples from all :math:`\xi^{(i)}` of the same sharing group into a
single histogram and update a common free energy estimate. Samples can
be shared between walkers within the simulation and/or between multiple
simulations. However, currently only sharing between simulations is
supported in the code while all biases within a simulation are
independent.

Note that when attempting to shorten the simulation time by using
bias-sharing walkers, care must be taken to ensure the simulations are
still long enough to properly explore and equilibrate all regions of the
sampling interval. To begin, the walkers in a group should be
decorrelated and distributed approximately according to the target
distribution before starting to refine the free energy. This can be
achieved e.g. by “equilibrating” the shared weight histogram before
letting it grow; for instance, :math:`W(\lambda)/N\approx \rho(\lambda)`
with some tolerance.

Furthermore, the “covering” or transition criterion of the initial stage
should to be generalized to detect when the sampling interval has been
collectively traversed. One alternative is to just use the same
criterion as for a single walker (but now with more samples), see
Eq. [eq:awh:covering]. However, in contrast to the single walker case
this does not ensure that any real transitions across the sampling
interval has taken place; in principle all walkers could be sampling
only very locally and still cover the whole interval. Just as with a
standard umbrella sampling procedure, the free energy may appear to be
converged while in reality simulations sampling closeby :math:`\lambda`
values are sampling disconnected regions of phase space. A stricter
criterion, which helps avoid such issues, is to require that before a
simulation marks a point :math:`\lambda_\mu` along dimension :math:`\mu`
as visited, and shares this with the other walkers, also all points
within a certain diameter :math:`D_{\mathrm{cover}}` should have been
visited (i.e.fulfill Eq. [eq:awh:covering]). Increasing
:math:`D_{\mathrm{cover}}` increases robustness, but may slow down
convergence. For the maximum value of :math:`D_{\mathrm{cover}}`, equal
to the length of the sampling interval, the sampling interval is
considered covered when at least one walker has independently traversed
the sampling interval.

Reweighting and combining biased data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Often one may want to, post-simulation, calculate the unbiased PMF
:math:`\Phi(u)` of another variable :math:`u(x)`. :math:`\Phi(u)` can be
estimated using :math:`\xi`-biased data by reweighting (“unbiasing”) the
trajectory using the bias potential :math:`U_{n(t)}`, see
Eq. [eq:awh:biaspotential]. Essentially, one bins the biased data along
:math:`u` and removes the effect of :math:`U_{n(t)}` by dividing the
weight of samples :math:`u(t)` by :math:`e^{-U_{n(t)}(\xi(t))}`,

.. math::

   \label{eq:awh:unbias}
   \hat{\Phi}(u)  = -\ln 
   \sum_t 1_u(u(t))e^{U_{n(t)}(\xi(t)} \mathcal{Z}_{n(t)}.

 Here the indicator function :math:`1_u` denotes the binning procedure:
:math:`1_u(u') = 1` if :math:`u'` falls into the bin labeled by
:math:`u` and :math:`0` otherwise. The normalization factor
:math:`\mathcal{Z}_n = \int e^{-\Phi(\xi) - U_{n}(\xi)}d \xi` is the
partition function of the extended ensemble. As can be seen
:math:`\mathcal{Z}_n` depends on :math:`\Phi(\xi)`, the PMF of the
(biased) reaction coordinate :math:`\xi` (which is calculated and
written to file by the AWH simulation). It is advisable to use only
final stage data in the reweighting procedure due to the rapid change of
the bias potential during the initial stage. If one would include
initial stage data, one should use the sample weights that are inferred
by the repeated rescaling of the histogram in the initial stage, for the
sake of consistency. Initial stage samples would then in any case be
heavily scaled down relative to final stage samples. Note that
Eq. [eq:awh:unbias] can also be used to combine data from multiple
simulations (by adding another sum also over the trajectory set).
Furthermore, when multiple independent AWH biases have generated a set
of PMF estimates :math:`\{\hat{\Phi}^{(i)}(\xi)\}`, a combined best
estimate :math:`\hat{\Phi}(\xi)` can be obtained by applying
self-consistent exponential averaging. More details on this procedure
and a derivation of Eq. [eq:awh:unbias] (using slightly different
notation) can be found in `143 <#ref-lindahl2017sequence>`__.

The friction metric
~~~~~~~~~~~~~~~~~~~

During the AWH simulation, the following time-integrated force
correlation function is calculated,

.. math::

   \label{eq:awh:metric}
   \eta_{\mu\nu}(\lambda) =
   \beta
   \int_0^\infty
   \frac{
   {\left<{\delta \mathcal{F}_{\mu}(x(t),\lambda)
   \delta \mathcal{F}_\nu(x(0),\lambda)
   \omega(\lambda|x(t)) \omega(\lambda|x(0))}\right>}}
   {{\left<{\omega^2(\lambda|x)}\right>}}
   dt.

Here
:math:`\mathcal F_\mu(x,\lambda) = k_\mu (\xi_\mu(x) - \lambda_\mu)` is
the force along dimension :math:`\mu` from an harmonic potential
centered at :math:`\lambda` and
:math:`\delta \mathcal F_\mu(x,\lambda) = \mathcal F_\mu(x,\lambda) - {\left<{\mathcal F_\mu(x,\lambda)}\right>}`
is the deviation of the force. The factors :math:`\omega(\lambda|x(t))`,
see Eq. [eq:awh:omega], reweight the samples.
:math:`\eta_{\mu\nu}(\lambda)` is a friction
tensor \ `144 <#ref-sivak2012thermodynamic>`__. Its matrix elements are
inversely proportional to local diffusion coefficients. A measure of
sampling (in)efficiency at each :math:`\lambda` is given by

.. math::

   \label{eq:awh:sqrt-metric}
   \eta^{\frac{1}{2}}(\lambda) = \sqrt{\det\eta_{\mu\nu}(\lambda)}.

 A large value of :math:`\eta^{\frac{1}{2}}(\lambda)` indicates slow
dynamics and long correlation times, which may require more sampling.

Usage
~~~~~

AWH stores data in the energy file (.edr) with a frequency set by the
user. The data – the PMF, the convolved bias, distributions of the
:math:`\lambda` and :math:`\xi` coordinates, etc. – can be extracted
after the simulation using the gmx awh tool. Furthermore, the trajectory
of the reaction coordinate :math:`\xi(t)` is printed to the pull output
file :math:`{\tt pullx.xvg}`. The log file (.log) also contains
information; check for messages starting with “awh”, they will tell you
about covering and potential sampling issues.

Setting the initial update size
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The initial value of the weight histogram size :math:`N` sets the
initial update size (and the rate of change of the bias). When :math:`N`
is kept constant, like in the initial stage, the average variance of the
free energy scales as :math:`\varepsilon^2 \sim 1/(ND)`
`137 <#ref-lindahl2014accelerated>`__, for a simple model system with
constant diffusion :math:`D` along the reaction coordinate. This
provides a ballpark estimate used by AWH to initialize :math:`N` in
terms of more meaningful quantities

.. math::

   \label{eq:awh:n0}
   \frac{1}{N_0} = \frac{1}{N_0(\varepsilon_0, D)} \sim D\varepsilon_0^2.

 Essentially, this tells us that a slower system (small :math:`D`)
requires more samples (larger :math:`N^0`) to attain the same level of
accuracy (:math:`\varepsilon_0`) at a given sampling rate. Conversely,
for a system of given diffusion, how to choose the initial biasing rate
depends on how good the initial accuracy is. Both the initial error
:math:`\varepsilon_0` and the diffusion :math:`D` only need to be
roughly estimated or guessed. In the typical case, one would only tweak
the :math:`D` parameter, and use a default value for
:math:`\varepsilon_0`. For good convergence, :math:`D` should be chosen
as large as possible (while maintaining a stable system) giving large
initial bias updates and fast initial transitions. Choosing :math:`D`
too small can lead to slow initial convergence. It may be a good idea to
run a short trial simulation and after the first covering check the
maximum free energy difference of the PMF estimate. If this is much
larger than the expected magnitude of the free energy barriers that
should be crossed, then the system is probably being pulled too hard and
:math:`D` should be decreased. :math:`\varepsilon_0` on the other hand,
would only be tweaked when starting an AWH simulation using a fairly
accurate guess of the PMF as input.

Tips for efficient sampling
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The force constant :math:`k` should be larger than the curvature of the
PMF landscape. If this is not the case, the distributions of the
reaction coordinate :math:`\xi` and the reference coordinate
:math:`\lambda`, will differ significantly and warnings will be printed
in the log file. One can choose :math:`k` as large as the time step
supports. This will neccessarily increase the number of points of the
discretized sampling interval :math:`I`. In general however, it should
not affect the performance of the simulation noticeably because the AWH
update is implemented such that only sampled points are accessed at free
energy update time.

As with any method, the choice of reaction coordinate(s) is critical. If
a single reaction coordinate does not suffice, identifying a second
reaction coordinate and sampling the two-dimensional landscape may help.
In this case, using a target distribution with a free energy cutoff (see
Eq. [eq:awh:rhocut]) might be required to avoid sampling uninteresting
regions of very high free energy. Obtaining accurate free energies for
reaction coordinates of much higher dimensionality than 3 or possibly 4
is generally not feasible.

Monitoring the transition rate of :math:`\xi(t)`, across the sampling
interval is also advisable. For reliable statistics (e.g. when
reweighting the trajectory as described in section [sec:awh:reweight]),
one would generally want to observe at least a few transitions after
having exited the initial stage. Furthermore, if the dynamics of the
reaction coordinate suddenly changes, this may be a sign of e.g. a
reaction coordinate problem.

Difficult regions of sampling may also be detected by calculating the
friction tensor :math:`\eta_{\mu\nu}(\lambda)` in the sampling interval,
see section [sec:awh:friction]. :math:`\eta_{\mu\nu}(\lambda)` as well
as the sampling efficiency measure :math:`\eta^{\frac{1}{2}}(\lambda)`
(Eq. [eq:awh:sqrt-metric]) are written to the energy file and can be
extracted with gmx awh. A high peak in
:math:`\eta^{\frac{1}{2}}(\lambda)` indicates that this region requires
longer time to sample properly.

Enforced Rotation
-----------------

[sec:rotation]

="2D

This module can be used to enforce the rotation of a group of atoms, as
*e.g.* a protein subunit. There are a variety of rotation potentials,
among them complex ones that allow flexible adaptations of both the
rotated subunit as well as the local rotation axis during the
simulation. An example application can be found in ref.
`145 <#ref-Kutzner2011>`__.

.. figure:: plots/rotation.pdf
   :alt: Comparison of fixed and flexible axis rotation. A: Rotating the
   sketched shape inside the white tubular cavity can create artifacts
   when a fixed rotation axis (dashed) is used. More realistically, the
   shape would revolve like a flexible pipe-cleaner (dotted) inside the
   begineqnarrayring (gray). B: Fixed rotation around an axis with a
   pivot point specified by the vector . C: Subdividing the rotating
   fragment into slabs with separate rotation axes (:math:`\uparrow`)
   and pivot points (:math:`\bullet`) for each slab allows for
   flexibility. The distance between two slabs with indices :math:`n`
   and :math:`n+1` is :math:`\Delta x`.
   :width: 13.00000cm

   Comparison of fixed and flexible axis rotation. A: Rotating the
   sketched shape inside the white tubular cavity can create artifacts
   when a fixed rotation axis (dashed) is used. More realistically, the
   shape would revolve like a flexible pipe-cleaner (dotted) inside the
   begineqnarrayring (gray). B: Fixed rotation around an axis with a
   pivot point specified by the vector . C: Subdividing the rotating
   fragment into slabs with separate rotation axes (:math:`\uparrow`)
   and pivot points (:math:`\bullet`) for each slab allows for
   flexibility. The distance between two slabs with indices :math:`n`
   and :math:`n+1` is :math:`\Delta x`.

.. figure:: plots/equipotential.pdf
   :alt: Selection of different rotation potentials and definition of
   notation. All four potentials :math:`V` (color coded) are shown for a
   single atom at position :math:`{\mbox{\boldmath ${x}$}}_j(t)`. A:
   Isotropic potential :math:`V\rotiso`, B: radial motion potential
   :math:`V\rotrm` and flexible potential :math:`V\rotflex`, C–D: radial
   motion2 potential :math:`V\rotrmtwo` and flexible2 potential
   :math:`V\rotflextwo` for :math:`\epsilon' = 0`\ nm\ :math:`^2` (C)
   and :math:`\epsilon' = 0.01`\ nm\ :math:`^2` (D). The rotation axis
   is perpendicular to the plane and marked by :math:`\otimes`. The
   light gray contours indicate Boltzmann factors :math:`e^{-V/(k_B T)}`
   in the :math:`{\mbox{\boldmath ${x}$}}_j`-plane for :math:`T=300`\ K
   and :math:`k=200`\ kJ/(mol\ :math:`\cdot`\ nm\ :math:`^2`). The green
   arrow shows the direction of the force
   :math:`{\mbox{\boldmath ${F}$}}_{\!j}` acting on atom :math:`j`; the
   blue dashed line indicates the motion of the reference position.
   :width: 13.00000cm

   Selection of different rotation potentials and definition of
   notation. All four potentials :math:`V` (color coded) are shown for a
   single atom at position :math:`{\mbox{\boldmath ${x}$}}_j(t)`. A:
   Isotropic potential :math:`V\rotiso`, B: radial motion potential
   :math:`V\rotrm` and flexible potential :math:`V\rotflex`, C–D: radial
   motion2 potential :math:`V\rotrmtwo` and flexible2 potential
   :math:`V\rotflextwo` for :math:`\epsilon' = 0`\ nm\ :math:`^2` (C)
   and :math:`\epsilon' = 0.01`\ nm\ :math:`^2` (D). The rotation axis
   is perpendicular to the plane and marked by :math:`\otimes`. The
   light gray contours indicate Boltzmann factors :math:`e^{-V/(k_B T)}`
   in the :math:`{\mbox{\boldmath ${x}$}}_j`-plane for :math:`T=300`\ K
   and :math:`k=200`\ kJ/(mol\ :math:`\cdot`\ nm\ :math:`^2`). The green
   arrow shows the direction of the force
   :math:`{\mbox{\boldmath ${F}$}}_{\!j}` acting on atom :math:`j`; the
   blue dashed line indicates the motion of the reference position.

Fixed Axis Rotation
~~~~~~~~~~~~~~~~~~~

Stationary Axis with an Isotropic Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the fixed axis approach (see Fig. [fig:rotation]B), torque on a group
of :math:`N` atoms with positions :math:`{\mbox{\boldmath ${x}$}}_i`
(denoted “rotation group”) is applied by rotating a reference set of
atomic positions – usually their initial positions
:math:`{\mbox{\boldmath ${y}$}}_i^0` – at a constant angular velocity
:math:`\omega` around an axis defined by a direction vector
:math:`\hat{{\mbox{\boldmath ${v}$}}}` and a pivot point . To that aim,
each atom with position :math:`{\mbox{\boldmath ${x}$}}_i` is attracted
by a “virtual spring” potential to its moving reference position
:math:`{\mbox{\boldmath ${y}$}}_i = \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})`,
where :math:`\mathbf{\Omega}(t)` is a matrix that describes the rotation
around the axis. In the simplest case, the “springs” are described by a
harmonic potential,

.. math::

   V\rotiso = \frac{k}{2} \sum_{i=1}^{N} w_i \left[ \mathbf{\Omega}(t)
   ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}})  \right]^2 ,
   \label{eqn:potiso}

 with optional mass-weighted prefactors :math:`w_i = N \, m_i/M` with
total mass :math:`M = \sum_{i=1}^N m_i`. The rotation matrix
:math:`\mathbf{\Omega}(t)` is

.. math::

   \mathbf{\Omega}(t) =  
   \left(   
   \begin{array}{ccc}
   \cos\omega t + v_x^2{\,\xi\,}& v_x v_y{\,\xi\,}- v_z\sin\omega t  & v_x v_z{\,\xi\,}+ v_y\sin\omega t\\
   v_x v_y{\,\xi\,}+ v_z\sin\omega t  & \cos\omega t + v_y^2{\,\xi\,}& v_y v_z{\,\xi\,}- v_x\sin\omega t\\
   v_x v_z{\,\xi\,}- v_y\sin\omega t  & v_y v_z{\,\xi\,}+ v_x\sin\omega t  & \cos\omega t + v_z^2{\,\xi\,}\\
   \end{array}
   \right) ,

where :math:`v_x`, :math:`v_y`, and :math:`v_z` are the components of
the normalized rotation vector :math:`\hat{{\mbox{\boldmath ${v}$}}}`,
and :math:`{\,\xi\,}:= 1-\cos(\omega t)`. As illustrated in
Fig. [fig:equipotential]A for a single atom :math:`j`, the rotation
matrix :math:`\mathbf{\Omega}(t)` operates on the initial reference
positions
:math:`{\mbox{\boldmath ${y}$}}_j^0 = {\mbox{\boldmath ${x}$}}_j(t_0)`
of atom :math:`j` at :math:`t=t_0`. At a later time :math:`t`, the
reference position has rotated away from its initial place (along the
blue dashed line), resulting in the force

.. math::

   {\mbox{\boldmath ${F}$}}_{\!j}\rotiso 
   = -\nabla_{\!j} \, V\rotiso 
   = k \, w_j \left[
   \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_j^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_j - {\mbox{\boldmath ${u}$}} ) \right] ,
   \label{eqn:force_fixed}

 which is directed towards the reference position.

Pivot-Free Isotropic Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Instead of a fixed pivot vector this potential uses the center of mass
:math:`{\mbox{\boldmath ${x}$}}_c` of the rotation group as pivot for
the rotation axis,

.. math::

   {\mbox{\boldmath ${x}$}}_c   = \frac{1}{M} \sum_{i=1}^N m_i {\mbox{\boldmath ${x}$}}_i 
   \label{eqn:com}
   \mbox{\hspace{4ex}and\hspace{4ex}}
   {\mbox{\boldmath ${y}$}}_c^0 = \frac{1}{M} \sum_{i=1}^N m_i {\mbox{\boldmath ${y}$}}_i^0 \ ,

 which yields the “pivot-free” isotropic potential

.. math::

   V\rotisopf = \frac{k}{2} \sum_{i=1}^{N} w_i \left[ \mathbf{\Omega}(t)
   ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c) \right]^2 ,
   \label{eqn:potisopf}

 with forces

.. math::

   \mathbf{F}_{\!j}\rotisopf = k \, w_j 
   \left[ 
   \mathbf{\Omega}(t) ( {\mbox{\boldmath ${y}$}}_j^0 - {\mbox{\boldmath ${y}$}}_c^0) 
                    - ( {\mbox{\boldmath ${x}$}}_j   - {\mbox{\boldmath ${x}$}}_c )
   \right] .
   \label{eqn:force_isopf}

 Without mass-weighting, the pivot :math:`{\mbox{\boldmath ${x}$}}_c` is
the geometrical center of the group. [sec:fixed]

Parallel Motion Potential Variant
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The forces generated by the isotropic potentials (eqns. [eqn:potiso] and
[eqn:potisopf]) also contain components parallel to the rotation axis
and thereby restrain motions along the axis of either the whole rotation
group (in case of :math:`V\rotiso`) or within the rotation group (in
case of :math:`V\rotisopf`). For cases where unrestrained motion along
the axis is preferred, we have implemented a “parallel motion” variant
by eliminating all components parallel to the rotation axis for the
potential. This is achieved by projecting the distance vectors between
reference and actual positions

.. math:: {\mbox{\boldmath ${r}$}}_i = \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}})

 onto the plane perpendicular to the rotation vector,

.. math::

   \label{eqn:project}
   {\mbox{\boldmath ${r}$}}_i^\perp :=  {\mbox{\boldmath ${r}$}}_i - ({\mbox{\boldmath ${r}$}}_i \cdot \hat{{\mbox{\boldmath ${v}$}}})\hat{{\mbox{\boldmath ${v}$}}} \ ,

 yielding

.. math::

   \begin{aligned}
   \nonumber
   V\rotpm &=& \frac{k}{2} \sum_{i=1}^{N} w_i ( {\mbox{\boldmath ${r}$}}_i^\perp )^2 \\
           &=& \frac{k}{2} \sum_{i=1}^{N} w_i
    \left\lbrace
    \mathbf{\Omega}(t)
      ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}})  \right. \nonumber \\
   && \left. - \left\lbrace
   \left[ \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}}) \right] \cdot\hat{{\mbox{\boldmath ${v}$}}}
     \right\rbrace\hat{{\mbox{\boldmath ${v}$}}} \right\rbrace^2 ,
   \label{eqn:potpm}\end{aligned}

 and similarly

.. math::

   {\mbox{\boldmath ${F}$}}_{\!j}\rotpm = k \, w_j \, {\mbox{\boldmath ${r}$}}_j^\perp .
   \label{eqn:force_pm}

Pivot-Free Parallel Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Replacing in eqn. [eqn:potpm] the fixed pivot by the center of mass
:math:`{\mbox{\boldmath ${x_c}$}}` yields the pivot-free variant of the
parallel motion potential. With

.. math:: {\mbox{\boldmath ${s}$}}_i = \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0) - ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c)

 the respective potential and forces are

.. math::

   \begin{aligned}
   V\rotpmpf &=& \frac{k}{2} \sum_{i=1}^{N} w_i ( {\mbox{\boldmath ${s}$}}_i^\perp )^2 \ , \\
   \label{eqn:potpmpf}
   {\mbox{\boldmath ${F}$}}_{\!j}\rotpmpf &=& k \, w_j \, {\mbox{\boldmath ${s}$}}_j^\perp .
   \label{eqn:force_pmpf}\end{aligned}

Radial Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^

In the above variants, the minimum of the rotation potential is either a
single point at the reference position
:math:`{\mbox{\boldmath ${y}$}}_i` (for the isotropic potentials) or a
single line through :math:`{\mbox{\boldmath ${y}$}}_i` parallel to the
rotation axis (for the parallel motion potentials). As a result, radial
forces restrict radial motions of the atoms. The two subsequent types of
rotation potentials, :math:`V\rotrm` and :math:`V\rotrmtwo`, drastically
reduce or even eliminate this effect. The first variant, :math:`V\rotrm`
(Fig. [fig:equipotential]B), eliminates all force components parallel to
the vector connecting the reference atom and the rotation axis,

.. math::

   V\rotrm = \frac{k}{2} \sum_{i=1}^{N} w_i \left[
   {\mbox{\boldmath ${p}$}}_i
   \cdot({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}}) \right]^2 ,
   \label{eqn:potrm}

 with

.. math::

   {\mbox{\boldmath ${p}$}}_i := 
   \frac{\hat{{\mbox{\boldmath ${v}$}}}\times \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})} {\| \hat{{\mbox{\boldmath ${v}$}}}\times \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})\|} \ .

 This variant depends only on the distance
:math:`{\mbox{\boldmath ${p}$}}_i \cdot ({\mbox{\boldmath ${x}$}}_i -
{\mbox{\boldmath ${u}$}})` of atom :math:`i` from the plane spanned by
:math:`\hat{{\mbox{\boldmath ${v}$}}}` and
:math:`\mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})`.
The resulting force is

.. math::

   \mathbf{F}_{\!j}\rotrm =
    -k \, w_j \left[ {\mbox{\boldmath ${p}$}}_j\cdot({\mbox{\boldmath ${x}$}}_j - {\mbox{\boldmath ${u}$}}) \right] \,{\mbox{\boldmath ${p}$}}_j \,  .
   \label{eqn:potrm_force}

Pivot-Free Radial Motion Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Proceeding similar to the pivot-free isotropic potential yields a
pivot-free version of the above potential. With

.. math::

   {\mbox{\boldmath ${q}$}}_i := 
   \frac{\hat{{\mbox{\boldmath ${v}$}}}\times \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0)} {\| \hat{{\mbox{\boldmath ${v}$}}}\times \mathbf{\Omega}(t) ({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0)\|} \, ,

 the potential and force for the pivot-free variant of the radial motion
potential read

.. math::

   \begin{aligned}
   V\rotrmpf & = & \frac{k}{2} \sum_{i=1}^{N} w_i \left[
   {\mbox{\boldmath ${q}$}}_i
   \cdot({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c)
   \right]^2 \, , \\
   \label{eqn:potrmpf}
   \mathbf{F}_{\!j}\rotrmpf & = &
    -k \, w_j \left[ {\mbox{\boldmath ${q}$}}_j\cdot({\mbox{\boldmath ${x}$}}_j - {\mbox{\boldmath ${x}$}}_c) \right] \,{\mbox{\boldmath ${q}$}}_j 
    + k   \frac{m_j}{M} \sum_{i=1}^{N} w_i \left[
    {\mbox{\boldmath ${q}$}}_i\cdot({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c) \right]\,{\mbox{\boldmath ${q}$}}_i \, .
   \label{eqn:potrmpf_force}\end{aligned}

Radial Motion 2 Alternative Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

As seen in Fig. [fig:equipotential]B, the force resulting from
:math:`V\rotrm` still contains a small, second-order radial component.
In most cases, this perturbation is tolerable; if not, the following
alternative, :math:`V\rotrmtwo`, fully eliminates the radial
contribution to the force, as depicted in Fig. [fig:equipotential]C,

.. math::

   V\rotrmtwo = 
   \frac{k}{2} \sum_{i=1}^{N} w_i\, 
   \frac{\left[ (\hat{{\mbox{\boldmath ${v}$}}} \times ( {\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}} ))
   \cdot \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}}) \right]^2}
   {\| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}}) \|^2 +
   \epsilon'} \, ,
   \label{eqn:potrm2}

 where a small parameter :math:`\epsilon'` has been introduced to avoid
singularities. For :math:`\epsilon'=0`\ nm\ :math:`^2`, the
equipotential planes are spanned by :math:`{\mbox{\boldmath ${x}$}}_i -
{\mbox{\boldmath ${u}$}}` and :math:`\hat{{\mbox{\boldmath ${v}$}}}`,
yielding a force perpendicular to
:math:`{\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}}`, thus not
contracting or expanding structural parts that moved away from or toward
the rotation axis.

Choosing a small positive :math:`\epsilon'` (*e.g.*,
:math:`\epsilon'=0.01`\ nm\ :math:`^2`, Fig. [fig:equipotential]D) in
the denominator of eqn. [eqn:potrm2] yields a well-defined potential and
continuous forces also close to the rotation axis, which is not the case
for :math:`\epsilon'=0`\ nm\ :math:`^2` (Fig. [fig:equipotential]C).
With

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${r}$}}_i & := & \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${u}$}})\\
   {\mbox{\boldmath ${s}$}}_i & := & \frac{\hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i -
   {\mbox{\boldmath ${u}$}} ) }{ \| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${u}$}})
   \| } \equiv \; \Psi_{i} \;\; {\hat{{\mbox{\boldmath ${v}$}}} \times
   ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${u}$}} ) }\\
   \Psi_i^{*}   & := & \frac{1}{ \| \hat{{\mbox{\boldmath ${v}$}}} \times
   ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${u}$}}) \|^2 + \epsilon'}\end{aligned}

 the force on atom :math:`j` reads

.. math::

   {\mbox{\boldmath ${F}$}}_{\!j}\rotrmtwo  = 
   - k\; 
   \left\lbrace w_j\;
   ({\mbox{\boldmath ${s}$}}_j\cdot{\mbox{\boldmath ${r}$}}_{\!j})\;
   \left[ \frac{\Psi_{\!j}^*   }{\Psi_{\!j}  }  {\mbox{\boldmath ${r}$}}_{\!j} 
        - \frac{\Psi_{\!j}^{*2}}{\Psi_{\!j}^3}
        ({\mbox{\boldmath ${s}$}}_j\cdot{\mbox{\boldmath ${r}$}}_{\!j}){\mbox{\boldmath ${s}$}}_j \right]
   \right\rbrace \times \hat{{\mbox{\boldmath ${v}$}}} .
   \label{eqn:potrm2_force}

Pivot-Free Radial Motion 2 Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The pivot-free variant of the above potential is

.. math::

   V{^\mathrm{rm2\mhyphen pf}}= 
   \frac{k}{2} \sum_{i=1}^{N} w_i\, 
   \frac{\left[ (\hat{{\mbox{\boldmath ${v}$}}} \times ( {\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c ))
   \cdot \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c) \right]^2}
   {\| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c) \|^2 +
   \epsilon'} \, .
   \label{eqn:potrm2pf}

 With

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${r}$}}_i & := & \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c)\\
   {\mbox{\boldmath ${s}$}}_i & := & \frac{\hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i -
   {\mbox{\boldmath ${x}$}}_c ) }{ \| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c)
   \| } \equiv \; \Psi_{i} \;\; {\hat{{\mbox{\boldmath ${v}$}}} \times
   ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${x}$}}_c ) }\\ \Psi_i^{*}   & := & \frac{1}{ \| \hat{{\mbox{\boldmath ${v}$}}} \times
   ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${x}$}}_c) \|^2 + \epsilon'}\end{aligned}

 the force on atom :math:`j` reads

.. math::

   \begin{aligned}
   \nonumber
   {\mbox{\boldmath ${F}$}}_{\!j}{^\mathrm{rm2\mhyphen pf}}& = &
   - k\; 
   \left\lbrace w_j\;
   ({\mbox{\boldmath ${s}$}}_j\cdot{\mbox{\boldmath ${r}$}}_{\!j})\;
   \left[ \frac{\Psi_{\!j}^*   }{\Psi_{\!j}  } {\mbox{\boldmath ${r}$}}_{\!j} 
        - \frac{\Psi_{\!j}^{*2}}{\Psi_{\!j}^3}
        ({\mbox{\boldmath ${s}$}}_j\cdot{\mbox{\boldmath ${r}$}}_{\!j}){\mbox{\boldmath ${s}$}}_j \right]
   \right\rbrace \times \hat{{\mbox{\boldmath ${v}$}}}\\
        & &
   + k\;\frac{m_j}{M} \left\lbrace \sum_{i=1}^{N}
   w_i\;({\mbox{\boldmath ${s}$}}_i\cdot{\mbox{\boldmath ${r}$}}_i) \; 
   \left[ \frac{\Psi_i^*   }{\Psi_i  }  {\mbox{\boldmath ${r}$}}_i
        - \frac{\Psi_i^{*2}}{\Psi_i^3} ({\mbox{\boldmath ${s}$}}_i\cdot{\mbox{\boldmath ${r}$}}_i )\;
        {\mbox{\boldmath ${s}$}}_i \right] \right\rbrace \times \hat{{\mbox{\boldmath ${v}$}}} \, .
   \label{eqn:potrm2pf_force}\end{aligned}

Flexible Axis Rotation
~~~~~~~~~~~~~~~~~~~~~~

As sketched in Fig. [fig:rotation]A–B, the rigid body behavior of the
fixed axis rotation scheme is a drawback for many applications. In
particular, deformations of the rotation group are suppressed when the
equilibrium atom positions directly depend on the reference positions.
To avoid this limitation, eqns. [eqn:potrmpf] and [eqn:potrm2pf] will
now be generalized towards a “flexible axis” as sketched in
Fig. [fig:rotation]C. This will be achieved by subdividing the rotation
group into a set of equidistant slabs perpendicular to the rotation
vector, and by applying a separate rotation potential to each of these
slabs. Fig. [fig:rotation]C shows the midplanes of the slabs as dotted
straight lines and the centers as thick black dots.

To avoid discontinuities in the potential and in the forces, we define
“soft slabs” by weighing the contributions of each slab :math:`n` to the
total potential function :math:`V\rotflex` by a Gaussian function

.. math::

   \label{eqn:gaussian}
   g_n({\mbox{\boldmath ${x}$}}_i) = \Gamma \ \mbox{exp} \left(
   -\frac{\beta_n^2({\mbox{\boldmath ${x}$}}_i)}{2\sigma^2}  \right) ,

 centered at the midplane of the :math:`n`\ th slab. Here :math:`\sigma`
is the width of the Gaussian function, :math:`\Delta x` the distance
between adjacent slabs, and

.. math:: \beta_n({\mbox{\boldmath ${x}$}}_i) := {\mbox{\boldmath ${x}$}}_i \cdot \hat{{\mbox{\boldmath ${v}$}}} - n \, \Delta x \, .

.. figure:: plots/gaussians.pdf
   :alt: Gaussian functions :math:`g_n` centered at
   :math:`n \, \Delta x` for a slab distance :math:`\Delta x = 1.5` nm
   and :math:`n \geq -2`. Gaussian function :math:`g_0` is highlighted
   in bold; the dashed line depicts the sum of the shown Gaussian
   functions.
   :width: 6.50000cm

   Gaussian functions :math:`g_n` centered at :math:`n \, \Delta x` for
   a slab distance :math:`\Delta x = 1.5` nm and :math:`n \geq -2`.
   Gaussian function :math:`g_0` is highlighted in bold; the dashed line
   depicts the sum of the shown Gaussian functions.

A most convenient choice is :math:`\sigma = 0.7 \Delta x` and

.. math::

   1/\Gamma = \sum_{n \in Z}
   \mbox{exp}
   \left(-\frac{(n - \frac{1}{4})^2}{2\cdot 0.7^2}\right)
   \approx 1.75464 \, ,

 which yields a nearly constant sum, essentially independent of
:math:`{\mbox{\boldmath ${x}$}}_i` (dashed line in
Fig. [fig:gaussians]), *i.e.*,

.. math::

   \sum_{n \in Z} g_n({\mbox{\boldmath ${x}$}}_i) =  1 + \epsilon({\mbox{\boldmath ${x}$}}_i) \, ,
   \label{eqn:normal}

 with
:math:` | \epsilon({\mbox{\boldmath ${x}$}}_i) | < 1.3\cdot 10^{-4}`.
This choice also implies that the individual contributions to the force
from the slabs add up to unity such that no further normalization is
required.

To each slab center :math:`{\mbox{\boldmath ${x}$}}_c^n`, all atoms
contribute by their Gaussian-weighted (optionally also mass-weighted)
position vectors
:math:`g_n({\mbox{\boldmath ${x}$}}_i) \, {\mbox{\boldmath ${x}$}}_i`.
The instantaneous slab centers :math:`{\mbox{\boldmath ${x}$}}_c^n` are
calculated from the current positions
:math:`{\mbox{\boldmath ${x}$}}_i`,

.. math::

   \label{eqn:defx0} 
   {\mbox{\boldmath ${x}$}}_c^n =
   \frac{\sum_{i=1}^N g_n({\mbox{\boldmath ${x}$}}_i) \, m_i \, {\mbox{\boldmath ${x}$}}_i}
        {\sum_{i=1}^N g_n({\mbox{\boldmath ${x}$}}_i) \, m_i} \, ,\\

 while the reference centers :math:`{\mbox{\boldmath ${y}$}}_c^n` are
calculated from the reference positions
:math:`{\mbox{\boldmath ${y}$}}_i^0`,

.. math::

   \label{eqn:defy0}
   {\mbox{\boldmath ${y}$}}_c^n =
   \frac{\sum_{i=1}^N g_n({\mbox{\boldmath ${y}$}}_i^0) \, m_i \, {\mbox{\boldmath ${y}$}}_i^0}
        {\sum_{i=1}^N g_n({\mbox{\boldmath ${y}$}}_i^0) \, m_i} \, .

 Due to the rapid decay of :math:`g_n`, each slab will essentially
involve contributions from atoms located within :math:`\approx
3\Delta x` from the slab center only.

Flexible Axis Potential
^^^^^^^^^^^^^^^^^^^^^^^

We consider two flexible axis variants. For the first variant, the slab
segmentation procedure with Gaussian weighting is applied to the radial
motion potential (eqn. [eqn:potrmpf]/Fig. [fig:equipotential]B),
yielding as the contribution of slab :math:`n`

.. math::

   V^n = 
   \frac{k}{2} \sum_{i=1}^{N} w_i \, g_n({\mbox{\boldmath ${x}$}}_i) 
   \left[
   {\mbox{\boldmath ${q}$}}_i^n
   \cdot
    ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n) 
   \right]^2  ,
   \label{eqn:flexpot}

 and a total potential function

.. math::

   V\rotflex = \sum_n V^n \, .
   \label{eqn:potflex}

 Note that the global center of mass :math:`{\mbox{\boldmath ${x}$}}_c`
used in eqn. [eqn:potrmpf] is now replaced by
:math:`{\mbox{\boldmath ${x}$}}_c^n`, the center of mass of the slab.
With

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${q}$}}_i^n & := & \frac{\hat{{\mbox{\boldmath ${v}$}}} \times
   \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^n) }{ \| \hat{{\mbox{\boldmath ${v}$}}}
   \times \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^n) \| } \\
   b_i^n         & := & {\mbox{\boldmath ${q}$}}_i^n \cdot ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n) \, ,\end{aligned}

 the resulting force on atom :math:`j` reads

.. math::

   \begin{aligned}
   \nonumber\hspace{-15mm}
   {\mbox{\boldmath ${F}$}}_{\!j}\rotflex &=&
   - \, k \, w_j \sum_n g_n({\mbox{\boldmath ${x}$}}_j) \, b_j^n \left\lbrace  {\mbox{\boldmath ${q}$}}_j^n -
   b_j^n \frac{\beta_n({\mbox{\boldmath ${x}$}}_j)}{2\sigma^2} \hat{{\mbox{\boldmath ${v}$}}} \right\rbrace \\ & &
   + \, k \, m_j \sum_n \frac{g_n({\mbox{\boldmath ${x}$}}_j)}{\sum_h g_n({\mbox{\boldmath ${x}$}}_h)}
   \sum_{i=1}^{N} w_i \, g_n({\mbox{\boldmath ${x}$}}_i) \, b_i^n \left\lbrace 
   {\mbox{\boldmath ${q}$}}_i^n -\frac{\beta_n({\mbox{\boldmath ${x}$}}_j)}{\sigma^2}
   \left[ {\mbox{\boldmath ${q}$}}_i^n \cdot ({\mbox{\boldmath ${x}$}}_j - {\mbox{\boldmath ${x}$}}_c^n )\right]
   \hat{{\mbox{\boldmath ${v}$}}} \right\rbrace .
   \label{eqn:potflex_force}\end{aligned}

 Note that for :math:`V\rotflex`, as defined, the slabs are fixed in
space and so are the reference centers
:math:`{\mbox{\boldmath ${y}$}}_c^n`. If during the simulation the
rotation group moves too far in :math:`{\mbox{\boldmath ${v}$}}`
direction, it may enter a region where – due to the lack of nearby
reference positions – no reference slab centers are defined, rendering
the potential evaluation impossible. We therefore have included a
slightly modified version of this potential that avoids this problem by
attaching the midplane of slab :math:`n=0` to the center of mass of the
rotation group, yielding slabs that move with the rotation group. This
is achieved by subtracting the center of mass
:math:`{\mbox{\boldmath ${x}$}}_c` of the group from the positions,

.. math::

   \tilde{{\mbox{\boldmath ${x}$}}}_i = {\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c \, , \mbox{\ \ \ and \ \ } 
   \tilde{{\mbox{\boldmath ${y}$}}}_i^0 = {\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^0 \, ,
   \label{eqn:trafo}

 such that

.. math::

   \begin{aligned}
   V\rotflext 
     & = & \frac{k}{2} \sum_n \sum_{i=1}^{N} w_i \, g_n(\tilde{{\mbox{\boldmath ${x}$}}}_i)
     \left[ \frac{\hat{{\mbox{\boldmath ${v}$}}} \times \mathbf{\Omega}(t)(\tilde{{\mbox{\boldmath ${y}$}}}_i^0
     - \tilde{{\mbox{\boldmath ${y}$}}}_c^n) }{ \| \hat{{\mbox{\boldmath ${v}$}}} \times
   \mathbf{\Omega}(t)(\tilde{{\mbox{\boldmath ${y}$}}}_i^0 -
   \tilde{{\mbox{\boldmath ${y}$}}}_c^n) \| }
   \cdot
    (\tilde{{\mbox{\boldmath ${x}$}}}_i - \tilde{{\mbox{\boldmath ${x}$}}}_c^n) 
   \right]^2 .
   \label{eqn:potflext}\end{aligned}

 To simplify the force derivation, and for efficiency reasons, we here
assume :math:`{\mbox{\boldmath ${x}$}}_c` to be constant, and thus
:math:`\partial {\mbox{\boldmath ${x}$}}_c / \partial x =
\partial {\mbox{\boldmath ${x}$}}_c / \partial y = \partial {\mbox{\boldmath ${x}$}}_c / \partial z = 0`.
The resulting force error is small (of order :math:`O(1/N)` or
:math:`O(m_j/M)` if mass-weighting is applied) and can therefore be
tolerated. With this assumption, the forces
:math:`{\mbox{\boldmath ${F}$}}\rotflext` have the same form as
eqn. [eqn:potflex\_force].

Flexible Axis 2 Alternative Potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this second variant, slab segmentation is applied to
:math:`V\rotrmtwo` (eqn. [eqn:potrm2pf]), resulting in a flexible axis
potential without radial force contributions
(Fig. [fig:equipotential]C),

.. math::

   V{^\mathrm{flex2}}= 
   \frac{k}{2} \sum_{i=1}^{N} \sum_n w_i\,g_n({\mbox{\boldmath ${x}$}}_i) 
   \frac{\left[ (\hat{{\mbox{\boldmath ${v}$}}} \times ( {\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n ))
   \cdot \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^n) \right]^2}
   {\| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n) \|^2 +
   \epsilon'} \, .
   \label{eqn:potflex2}

 With

.. math::

   \begin{aligned}
   {\mbox{\boldmath ${r}$}}_i^n & := & \mathbf{\Omega}(t)({\mbox{\boldmath ${y}$}}_i^0 - {\mbox{\boldmath ${y}$}}_c^n)\\
   {\mbox{\boldmath ${s}$}}_i^n & := & \frac{\hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i -
   {\mbox{\boldmath ${x}$}}_c^n ) }{ \| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i - {\mbox{\boldmath ${x}$}}_c^n)
   \| } \equiv \; \psi_{i} \;\; {\hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${x}$}}_c^n ) }\\
   \psi_i^{*}     & := & \frac{1}{ \| \hat{{\mbox{\boldmath ${v}$}}} \times ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${x}$}}_c^n) \|^2 + \epsilon'}\\ 
   W_j^n          & := & \frac{g_n({\mbox{\boldmath ${x}$}}_j)\,m_j}{\sum_h g_n({\mbox{\boldmath ${x}$}}_h)\,m_h}\\
   {\mbox{\boldmath ${S}$}}^n   & := & 
   \sum_{i=1}^{N} w_i\;g_n({\mbox{\boldmath ${x}$}}_i)
   \; ({\mbox{\boldmath ${s}$}}_i^n\cdot{\mbox{\boldmath ${r}$}}_i^n)
   \left[ \frac{\psi_i^*   }{\psi_i  }  {\mbox{\boldmath ${r}$}}_i^n
        - \frac{\psi_i^{*2}}{\psi_i^3} ({\mbox{\boldmath ${s}$}}_i^n\cdot{\mbox{\boldmath ${r}$}}_i^n )\;
        {\mbox{\boldmath ${s}$}}_i^n \right] \label{eqn:Sn}\end{aligned}

 the force on atom :math:`j` reads

.. math::

   \begin{aligned}
   \nonumber
   {\mbox{\boldmath ${F}$}}_{\!j}{^\mathrm{flex2}}& = &
   - k\; 
   \left\lbrace \sum_n w_j\;g_n({\mbox{\boldmath ${x}$}}_j)\;
   ({\mbox{\boldmath ${s}$}}_j^n\cdot{\mbox{\boldmath ${r}$}}_{\!j}^n)\;
   \left[ \frac{\psi_j^*   }{\psi_j  }  {\mbox{\boldmath ${r}$}}_{\!j}^n 
        - \frac{\psi_j^{*2}}{\psi_j^3} ({\mbox{\boldmath ${s}$}}_j^n\cdot{\mbox{\boldmath ${r}$}}_{\!j}^n)\;
        {\mbox{\boldmath ${s}$}}_{\!j}^n \right] \right\rbrace \times \hat{{\mbox{\boldmath ${v}$}}} \\
   \nonumber
   & &
   + k \left\lbrace \sum_n W_{\!j}^n \, {\mbox{\boldmath ${S}$}}^n \right\rbrace \times
   \hat{{\mbox{\boldmath ${v}$}}}
   - k \left\lbrace \sum_n W_{\!j}^n \; \frac{\beta_n({\mbox{\boldmath ${x}$}}_j)}{\sigma^2} \frac{1}{\psi_j}\;\; 
   {\mbox{\boldmath ${s}$}}_j^n \cdot 
   {\mbox{\boldmath ${S}$}}^n \right\rbrace \hat{{\mbox{\boldmath ${v}$}}}\\ 
   & & 
   + \frac{k}{2} \left\lbrace \sum_n w_j\;g_n({\mbox{\boldmath ${x}$}}_j)
   \frac{\beta_n({\mbox{\boldmath ${x}$}}_j)}{\sigma^2} 
   \frac{\psi_j^*}{\psi_j^2}( {\mbox{\boldmath ${s}$}}_j^n \cdot {\mbox{\boldmath ${r}$}}_{\!j}^n )^2 \right\rbrace
   \hat{{\mbox{\boldmath ${v}$}}} .
   \label{eqn:potflex2_force}\end{aligned}

Applying transformation ([eqn:trafo]) yields a “translation-tolerant”
version of the flexible2 potential,
:math:`V{^\mathrm{flex2\mhyphen t}}`. Again, assuming that
:math:`\partial {\mbox{\boldmath ${x}$}}_c / \partial x`,
:math:`\partial {\mbox{\boldmath ${x}$}}_c /
\partial y`, :math:`\partial {\mbox{\boldmath ${x}$}}_c / \partial z`
are small, the resulting equations for
:math:`V{^\mathrm{flex2\mhyphen t}}` and
:math:`{\mbox{\boldmath ${F}$}}{^\mathrm{flex2\mhyphen t}}` are similar
to those of :math:`V\rotflextwo` and
:math:`{\mbox{\boldmath ${F}$}}\rotflextwo`.

Usage
~~~~~

To apply enforced rotation, the particles :math:`i` that are to be
subjected to one of the rotation potentials are defined via index groups
rot-group0, rot-group1, etc., in the .mdp input file. The reference
positions :math:`{\mbox{\boldmath ${y}$}}_i^0` are read from a special
.trr file provided to grompp. If no such file is found,
:math:`{\mbox{\boldmath ${x}$}}_i(t=0)` are used as reference positions
and written to .trr such that they can be used for subsequent setups.
All parameters of the potentials such as :math:`k`, :math:`\epsilon'`,
etc. (Table [tab:vars]) are provided as .mdp parameters; rot-type
selects the type of the potential. The option rot-massw allows to choose
whether or not to use mass-weighted averaging. For the flexible
potentials, a cutoff value :math:`g_n^\mathrm{min}` (typically
:math:`g_n^\mathrm{min}=0.001`) makes shure that only significant
contributions to :math:`V` and are evaluated, *i.e.* terms with
:math:`g_n({\mbox{\boldmath ${x}$}}) < g_n^\mathrm{min}` are omitted.
Table [tab:quantities] summarizes observables that are written to
additional output files and which are described below.

Angle of Rotation Groups: Fixed Axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For fixed axis rotation, the average angle :math:`\theta_\mathrm{av}(t)`
of the group relative to the reference group is determined via the
distance-weighted angular deviation of all rotation group atoms from
their reference positions,

.. math::

   \theta_\mathrm{av} = \left. \sum_{i=1}^{N} r_i \ \theta_i \right/ \sum_{i=1}^N r_i \ .
   \label{eqn:avangle}

 Here, :math:`r_i` is the distance of the reference position to the
rotation axis, and the difference angles :math:`\theta_i` are determined
from the atomic positions, projected onto a plane perpendicular to the
rotation axis through pivot point :math:`{\mbox{\boldmath ${u}$}}` (see
eqn. [eqn:project] for the definition of :math:`\perp`),

.. math::

   \cos \theta_i = 
   \frac{({\mbox{\boldmath ${y}$}}_i-{\mbox{\boldmath ${u}$}})^\perp \cdot ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${u}$}})^\perp}
        { \| ({\mbox{\boldmath ${y}$}}_i-{\mbox{\boldmath ${u}$}})^\perp \cdot ({\mbox{\boldmath ${x}$}}_i-{\mbox{\boldmath ${u}$}})^\perp
        \| } \ .

 The sign of :math:`\theta_\mathrm{av}` is chosen such that
:math:`\theta_\mathrm{av} > 0` if the actual structure rotates ahead of
the reference.

Angle of Rotation Groups: Flexible Axis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For flexible axis rotation, two outputs are provided, the angle of the
entire rotation group, and separate angles for the segments in the
slabs. The angle of the entire rotation group is determined by an RMSD
fit of :math:`{\mbox{\boldmath ${x}$}}_i` to the reference positions
:math:`{\mbox{\boldmath ${y}$}}_i^0` at :math:`t=0`, yielding
:math:`\theta_\mathrm{fit}` as the angle by which the reference has to
be rotated around :math:`\hat{{\mbox{\boldmath ${v}$}}}` for the optimal
fit,

.. math::

   \mathrm{RMSD} \big( {\mbox{\boldmath ${x}$}}_i,\ \mathbf{\Omega}(\theta_\mathrm{fit})
   {\mbox{\boldmath ${y}$}}_i^0 \big) \stackrel{!}{=} \mathrm{min} \, .
   \label{eqn:rmsdfit}

 To determine the local angle for each slab :math:`n`, both reference
and actual positions are weighted with the Gaussian function of slab
:math:`n`, and :math:`\theta_\mathrm{fit}(t,n)` is calculated as in
eqn. [eqn:rmsdfit]) from the Gaussian-weighted positions.

For all angles, the .mdp input option rot-fit-method controls whether a
normal RMSD fit is performed or whether for the fit each position
:math:`{\mbox{\boldmath ${x}$}}_i` is put at the same distance to the
rotation axis as its reference counterpart
:math:`{\mbox{\boldmath ${y}$}}_i^0`. In the latter case, the RMSD
measures only angular differences, not radial ones.

Angle Determination by Searching the Energy Minimum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, for rot-fit-method = potential, the angle of the rotation
group is determined as the angle for which the rotation potential energy
is minimal. Therefore, the used rotation potential is additionally
evaluated for a set of angles around the current reference angle. In
this case, the rotangles.log output file contains the values of the
rotation potential at the chosen set of angles, while rotation.xvg lists
the angle with minimal potential energy.

Torque
^^^^^^

The torque :math:`{\mbox{\boldmath ${\tau}$}}(t)` exerted by the
rotation potential is calculated for fixed axis rotation via

.. math::

   {\mbox{\boldmath ${\tau}$}}(t) = \sum_{i=1}^{N} {\mbox{\boldmath ${r}$}}_i(t) \times {\mbox{\boldmath ${f}$}}_{\!i}^\perp(t) ,
   \label{eqn:torque}

 where :math:`{\mbox{\boldmath ${r}$}}_i(t)` is the distance vector from
the rotation axis to :math:`{\mbox{\boldmath ${x}$}}_i(t)` and
:math:`{\mbox{\boldmath ${f}$}}_{\!i}^\perp(t)` is the force component
perpendicular to :math:`{\mbox{\boldmath ${r}$}}_i(t)` and
:math:`\hat{{\mbox{\boldmath ${v}$}}}`. For flexible axis rotation,
torques :math:`{\mbox{\boldmath ${\tau}$}}_{\!n}` are calculated for
each slab using the local rotation axis of the slab and the
Gaussian-weighted positions.

Electric fields
---------------

A pulsed and oscillating electric field can be applied according to:

.. math::

   E(t) = E_0 \exp\left[-\frac{(t-t_0)^2}{2\sigma^2}\right]\cos\left[\omega (t-t_0)\right]
   \label{eq_efield}

 where :math:`E_0` is the field strength, the angular frequency ,
:math:`t_0` is the time at of the peak in the field strength and
:math:`\sigma` is the with of the pulse. Special cases occur when
:math:`\sigma` = 0 (non-pulsed field) and for :math:`\omega` is 0
(static field).

This simulated laser-pulse was applied to simulations of melting
ice \ `146 <#ref-Caleman2008a>`__. A pulsed electric field may look ike
Fig. [fig:field]. In the supporting information of that paper the impact
of an applied electric field on a system under periodic boundary
conditions is analyzed. It is described that the effective electric
field under PBC is larger than the applied field, by a factor depending
on the size of the box and the dielectric properties of molecules in the
box. For a system with static dielectric properties this factor can be
corrected for. But for a system where the dielectric varies over time,
for example a membrane protein with a pore that opens and closes during
the simulatippn, this way of applying an electric field is not useful.
In such cases one can use the computational electrophysiology protocol
described in the next section (sec. [sec:compel]).

.. figure:: plots/field
   :alt: A simulated laser pulse in GROMACS.
   :width: 8.00000cm

   A simulated laser pulse in GROMACS.

Electric fields are applied when the following options are specified in
the grompp.mdp file. You specify, in order, :math:`E_0`, :math:`\omega`,
:math:`t_0` and :math:`\sigma`:

::

    electric-field-x = 0.04 0       0     0

yields a static field with :math:`E_0` = 0.04 V/nm in the X-direction.
In contrast,

::

    electric-field-x = 2.0  150     5     0

yields an oscillating electric field with :math:`E_0` = 2 V/nm,
:math:`\omega` = 150/ps and :math:`t_0` = 5 ps. Finally

::

    electric-field-x = 2.0  150     5     1

yields an pulsed-oscillating electric field with :math:`E_0` = 2 V/nm,
:math:`\omega` = 150/ps and :math:`t_0` = 5 ps and :math:`\sigma` = 1
ps. Read more in ref. \ `146 <#ref-Caleman2008a>`__. Note that the input
file format is changed from the undocumented older version. A figure
like Fig. [fig:field] may be produced by passing the -field option to
gmx mdrun.

Computational Electrophysiology
-------------------------------

The Computational Electrophysiology (CompEL) protocol
`147 <#ref-Kutzner2011b>`__ allows the simulation of ion flux through
membrane channels, driven by transmembrane potentials or ion
concentration gradients. Just as in real cells, CompEL establishes
transmembrane potentials by sustaining a small imbalance of charges
:math:`\Delta q` across the membrane, which gives rise to a potential
difference :math:`\Delta U` according to the membrane capacitance:

.. math:: \Delta U = \Delta q / C_{membrane}

 The transmembrane electric field and concentration gradients are
controlled by .mdp options, which allow the user to set reference counts
for the ions on either side of the membrane. If a difference between the
actual and the reference numbers persists over a certain time span,
specified by the user, a number of ion/water pairs are exchanged between
the compartments until the reference numbers are restored. Alongside the
calculation of channel conductance and ion selectivity, CompEL
simulations also enable determination of the channel reversal potential,
an important characteristic obtained in electrophysiology experiments.

In a CompEL setup, the simulation system is divided into two
compartments **A** and **B** with independent ion concentrations. This
is best achieved by using double bilayer systems with a copy (or copies)
of the channel/pore of interest in each bilayer (Fig. [fig:compelsetup]
A, B). If the channel axes point in the same direction, channel flux is
observed simultaneously at positive and negative potentials in this way,
which is for instance important for studying channel rectification.

.. figure:: plots/compelsetup.pdf
   :alt: Typical double-membrane setup for CompEL simulations (A, B).
   Ion/water molecule exchanges will be performed as needed between the
   two light blue volumes around the dotted black lines (A). Plot (C)
   shows the potential difference :math:`\Delta U` resulting from the
   selected charge imbalance :math:`\Delta q_{ref}` between the
   compartments.
   :width: 13.50000cm

   Typical double-membrane setup for CompEL simulations (A, B).
   Ion/water molecule exchanges will be performed as needed between the
   two light blue volumes around the dotted black lines (A). Plot (C)
   shows the potential difference :math:`\Delta U` resulting from the
   selected charge imbalance :math:`\Delta q_{ref}` between the
   compartments.

The potential difference :math:`\Delta U` across the membrane is easily
calculated with the gmx potential utility. By this, the potential drop
along :math:`z` or the pore axis is exactly known in each time interval
of the simulation (Fig. [fig:compelsetup] C). Type and number of ions
:math:`n_i` of charge :math:`q_i`, traversing the channel in the
simulation, are written to the swapions.xvg output file, from which the
average channel conductance :math:`G` in each interval :math:`\Delta t`
is determined by:

.. math:: G = \frac{\sum_{i} n_{i}q_{i}}{\Delta t \, \Delta U} \, .

 The ion selectivity is calculated as the number flux ratio of different
species. Best results are obtained by averaging these values over
several overlapping time intervals.

The calculation of reversal potentials is best achieved using a small
set of simulations in which a given transmembrane concentration gradient
is complemented with small ion imbalances of varying magnitude. For
example, if one compartment contains 1M salt and the other 0.1M, and
given charge neutrality otherwise, a set of simulations with
:math:`\Delta q = 0\,e`, :math:`\Delta q = 2\,e`,
:math:`\Delta q = 4\,e` could be used. Fitting a straight line through
the current-voltage relationship of all obtained :math:`I`-:math:`U`
pairs near zero current will then yield :math:`U_{rev}`.

Usage
~~~~~

The following .mdp options control the CompEL protocol:

::

    swapcoords     = Z        ; Swap positions: no, X, Y, Z
    swap-frequency = 100      ; Swap attempt frequency

Choose Z if your membrane is in the :math:`xy`-plane
(Fig. [fig:compelsetup]). Ions will be exchanged between compartments
depending on their :math:`z`-positions alone. swap-frequency determines
how often a swap attempt will be made. This step requires that the
positions of the split groups, the ions, and possibly the solvent
molecules are communicated between the parallel processes, so if chosen
too small it can decrease the simulation performance. The Position
swapping entry in the cycle and time accounting table at the end of the
md.log file summarizes the amount of runtime spent in the swap module.

::

    split-group0   = channel0 ; Defines compartment boundary
    split-group1   = channel1 ; Defines other compartment boundary
    massw-split0   = no       ; use mass-weighted center?
    massw-split1   = no

split-group0 and split-group1 are two index groups that define the
boundaries between the two compartments, which are usually the centers
of the channels. If massw-split0 or massw-split1 are set to yes, the
center of mass of each index group is used as boundary, here in
:math:`z`-direction. Otherwise, the geometrical centers will be used
(:math:`\times` in Fig. [fig:compelsetup] A). If, such as here, a
membrane channel is selected as split group, the center of the channel
will define the dividing plane between the compartments (dashed
horizontal lines). All index groups must be defined in the index file.

If, to restore the requested ion counts, an ion from one compartment has
to be exchanged with a water molecule from the other compartment, then
those molecules are swapped which have the largest distance to the
compartment-defining boundaries (dashed horizontal lines). Depending on
the ion concentration, this effectively results in exchanges of
molecules between the light blue volumes. If a channel is very
asymmetric in :math:`z`-direction and would extend into one of the swap
volumes, one can offset the swap exchange plane with the bulk-offset
parameter. A value of 0.0 means no offset :math:`b`, values
:math:`-1.0 < b < 0` move the swap exchange plane closer to the lower,
values :math:`0 < b < 1.0` closer to the upper membrane.
Fig. [fig:compelsetup] A (left) depicts that for the **A** compartment.

::

    solvent-group  = SOL      ; Group containing the solvent molecules
    iontypes       = 3        ; Number of different ion types to control
    iontype0-name  = NA       ; Group name of the ion type
    iontype0-in-A  = 51       ; Reference count of ions of type 0 in A
    iontype0-in-B  = 35       ; Reference count of ions of type 0 in B
    iontype1-name  = K
    iontype1-in-A  = 10
    iontype1-in-B  = 38
    iontype2-name  = CL
    iontype2-in-A  = -1
    iontype2-in-B  = -1

The group name of solvent molecules acting as exchange partners for the
ions has to be set with solvent-group. The number of different ionic
species under control of the CompEL protocol is given by the iontypes
parameter, while iontype0-name gives the name of the index group
containing the atoms of this ionic species. The reference number of ions
of this type can be set with the iontype0-in-A and iontype0-in-B options
for compartments **A** and **B**, respectively. Obviously, the sum of
iontype0-in-A and iontype0-in-B needs to equal the number of ions in the
group defined by iontype0-name. A reference number of -1 means: use the
number of ions as found at the beginning of the simulation as the
reference value.

::

    coupl-steps    = 10       ; Average over these many swap steps
    threshold      = 1        ; Do not swap if < threshold

If coupl-steps is set to 1, then the momentary ion distribution
determines whether ions are exchanged. coupl-steps > 1 will use the
time-average of ion distributions over the selected number of attempt
steps instead. This can be useful, for example, when ions diffuse near
compartment boundaries, which would lead to numerous unproductive ion
exchanges. A threshold of 1 means that a swap is performed if the
average ion count in a compartment differs by at least 1 from the
requested values. Higher thresholds will lead to toleration of larger
differences. Ions are exchanged until the requested number :math:`\pm`
the threshold is reached.

::

    cyl0-r         = 5.0      ; Split cylinder 0 radius (nm)
    cyl0-up        = 0.75     ; Split cylinder 0 upper extension (nm)
    cyl0-down      = 0.75     ; Split cylinder 0 lower extension (nm)
    cyl1-r         = 5.0      ; same for other channel
    cyl1-up        = 0.75
    cyl1-down      = 0.75

The cylinder options are used to define virtual geometric cylinders
around the channel’s pore to track how many ions of which type have
passed each channel. Ions will be counted as having traveled through a
channel according to the definition of the channel’s cylinder radius,
upper and lower extension, relative to the location of the respective
split group. This will not affect the actual flux or exchange, but will
provide you with the ion permeation numbers across each of the channels.
Note that an ion can only be counted as passing through a particular
channel if it is detected *within* the defined split cylinder in a swap
step. If swap-frequency is chosen too high, a particular ion may be
detected in compartment **A** in one swap step, and in compartment **B**
in the following swap step, so it will be unclear through which of the
channels it has passed.

A double-layered system for CompEL simulations can be easily prepared by
duplicating an existing membrane/channel MD system in the direction of
the membrane normal (typically :math:`z`) with gmx editconf -translate 0
0 <l\_z>, where l\_z is the box length in that direction. If you have
already defined index groups for the channel for the single-layered
system, gmx make\_ndx -n index.ndx -twin will provide you with the
groups for the double-layered system.

To suppress large fluctuations of the membranes along the swap
direction, it may be useful to apply a harmonic potential (acting only
in the swap dimension) between each of the two channel and/or bilayer
centers using umbrella pulling (see section [sec:pull]).

Multimeric channels
~~~~~~~~~~~~~~~~~~~

If a split group consists of more than one molecule, the correct PBC
image of all molecules with respect to each other has to be chosen such
that the channel center can be correctly determined. GROMACS assumes
that the starting structure in the .tpr file has the correct PBC
representation. Set the following environment variable to check whether
that is the case:

-  GMX\_COMPELDUMP: output the starting structure after it has been made
   whole to .pdb file.

Calculating a PMF using the free-energy code
--------------------------------------------

The free-energy coupling-parameter approach (see sec. [sec:fecalc])
provides several ways to calculate potentials of mean force. A potential
of mean force between two atoms can be calculated by connecting them
with a harmonic potential or a constraint. For this purpose there are
special potentials that avoid the generation of extra exclusions,
see sec. [sec:excl]. When the position of the minimum or the constraint
length is 1 nm more in state B than in state A, the restraint or
constraint force is given by :math:`\partial H/\partial \lambda`. The
distance between the atoms can be changed as a function of
:math:`\lambda` and time by setting delta-lambda in the .mdp file. The
results should be identical (although not numerically due to the
different implementations) to the results of the pull code with umbrella
sampling and constraint pulling. Unlike the pull code, the free energy
code can also handle atoms that are connected by constraints.

Potentials of mean force can also be calculated using position
restraints. With position restraints, atoms can be linked to a position
in space with a harmonic potential (see [subsec:positionrestraint]).
These positions can be made a function of the coupling parameter
:math:`\lambda`. The positions for the A and the B states are supplied
to grompp with the -r and -rb options, respectively. One could use this
approach to do targeted MD; note that we do not encourage the use of
targeted MD for proteins. A protein can be forced from one conformation
to another by using these conformations as position restraint
coordinates for state A and B. One can then slowly change
:math:`\lambda` from 0 to 1. The main drawback of this approach is that
the conformational freedom of the protein is severely limited by the
position restraints, independent of the change from state A to B. Also,
the protein is forced from state A to B in an almost straight line,
whereas the real pathway might be very different. An example of a more
fruitful application is a solid system or a liquid confined between
walls where one wants to measure the force required to change the
separation between the boundaries or walls. Because the boundaries (or
walls) already need to be fixed, the position restraints do not limit
the system in its sampling.

Removing fastest degrees of freedom
-----------------------------------

The maximum time step in MD simulations is limited by the smallest
oscillation period that can be found in the simulated system.
Bond-stretching vibrations are in their quantum-mechanical ground state
and are therefore better represented by a constraint instead of a
harmonic potential.

For the remaining degrees of freedom, the shortest oscillation period
(as measured from a simulation) is 13 fs for bond-angle vibrations
involving hydrogen atoms. Taking as a guideline that with a Verlet
(leap-frog) integration scheme a minimum of 5 numerical integration
steps should be performed per period of a harmonic oscillation in order
to integrate it with reasonable accuracy, the maximum time step will be
about 3 fs. Disregarding these very fast oscillations of period 13 fs,
the next shortest periods are around 20 fs, which will allow a maximum
time step of about 4 fs.

Removing the bond-angle degrees of freedom from hydrogen atoms can best
be done by defining them as virtual interaction sites instead of normal
atoms. Whereas a normal atom is connected to the molecule with bonds,
angles and dihedrals, a virtual site’s position is calculated from the
position of three nearby heavy atoms in a predefined manner (see also
sec. [sec:virtual\_sites]). For the hydrogens in water and in hydroxyl,
sulfhydryl, or amine groups, no degrees of freedom can be removed,
because rotational freedom should be preserved. The only other option
available to slow down these motions is to increase the mass of the
hydrogen atoms at the expense of the mass of the connected heavy atom.
This will increase the moment of inertia of the water molecules and the
hydroxyl, sulfhydryl, or amine groups, without affecting the equilibrium
properties of the system and without affecting the dynamical properties
too much. These constructions will shortly be described in
sec. [sec:vsitehydro] and have previously been described in full
detail \ `148 <#ref-feenstra99>`__.

Using both virtual sites and modified masses, the next bottleneck is
likely to be formed by the improper dihedrals (which are used to
preserve planarity or chirality of molecular groups) and the peptide
dihedrals. The peptide dihedral cannot be changed without affecting the
physical behavior of the protein. The improper dihedrals that preserve
planarity mostly deal with aromatic residues. Bonds, angles, and
dihedrals in these residues can also be replaced with somewhat elaborate
virtual site constructions.

All modifications described in this section can be performed using the
GROMACS topology building tool pdb2gmx. Separate options exist to
increase hydrogen masses, virtualize all hydrogen atoms, or also
virtualize all aromatic residues. **Note** that when all hydrogen atoms
are virtualized, those inside the aromatic residues will be virtualized
as well, *i.e.* hydrogens in the aromatic residues are treated
differently depending on the treatment of the aromatic residues.

Parameters for the virtual site constructions for the hydrogen atoms are
inferred from the force-field parameters (*vis*. bond lengths and
angles) directly by grompp while processing the topology file. The
constructions for the aromatic residues are based on the bond lengths
and angles for the geometry as described in the force fields, but these
parameters are hard-coded into pdb2gmx due to the complex nature of the
construction needed for a whole aromatic group.

Hydrogen bond-angle vibrations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Construction of virtual sites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. figure:: plots/dumtypes
   :alt: The different types of virtual site constructions used for
   hydrogen atoms. The atoms used in the construction of the virtual
   site(s) are depicted as black circles, virtual sites as gray ones.
   Hydrogens are smaller than heavy atoms. A: fixed bond angle, note
   that here the hydrogen is not a virtual site; B: in the plane of
   three atoms, with fixed distance; C: in the plane of three atoms,
   with fixed angle and distance; D: construction for amine groups
   (-NH:math:`_2` or -NH:math:`_3^+`), see text for details.
   :width: 11.00000cm

   The different types of virtual site constructions used for hydrogen
   atoms. The atoms used in the construction of the virtual site(s) are
   depicted as black circles, virtual sites as gray ones. Hydrogens are
   smaller than heavy atoms. A: fixed bond angle, note that here the
   hydrogen is not a virtual site; B: in the plane of three atoms, with
   fixed distance; C: in the plane of three atoms, with fixed angle and
   distance; D: construction for amine groups (-NH:math:`_2` or
   -NH:math:`_3^+`), see text for details.

The goal of defining hydrogen atoms as virtual sites is to remove all
high-frequency degrees of freedom from them. In some cases, not all
degrees of freedom of a hydrogen atom should be removed, *e.g.* in the
case of hydroxyl or amine groups the rotational freedom of the hydrogen
atom(s) should be preserved. Care should be taken that no unwanted
correlations are introduced by the construction of virtual sites, *e.g.*
bond-angle vibration between the constructing atoms could translate into
hydrogen bond-length vibration. Additionally, since virtual sites are by
definition massless, in order to preserve total system mass, the mass of
each hydrogen atom that is treated as virtual site should be added to
the bonded heavy atom.

Taking into account these considerations, the hydrogen atoms in a
protein naturally fall into several categories, each requiring a
different approach (see also Fig. [fig:vsitehydro]).

-  *hydroxyl (-OH) or sulfhydryl (-SH) hydrogen:* The only internal
   degree of freedom in a hydroxyl group that can be constrained is the
   bending of the C-O-H angle. This angle is fixed by defining an
   additional bond of appropriate length, see Fig. [fig:vsitehydro]A.
   Doing so removes the high-frequency angle bending, but leaves the
   dihedral rotational freedom. The same goes for a sulfhydryl group.
   **Note** that in these cases the hydrogen is not treated as a virtual
   site.

-  *single amine or amide (-NH-) and aromatic hydrogens (-CH-):* The
   position of these hydrogens cannot be constructed from a linear
   combination of bond vectors, because of the flexibility of the angle
   between the heavy atoms. Instead, the hydrogen atom is positioned at
   a fixed distance from the bonded heavy atom on a line going through
   the bonded heavy atom and a point on the line through both second
   bonded atoms, see Fig. [fig:vsitehydro]B.

-  *planar amine (-NH:math:`_2`) hydrogens:* The method used for the
   single amide hydrogen is not well suited for planar amine groups,
   because no suitable two heavy atoms can be found to define the
   direction of the hydrogen atoms. Instead, the hydrogen is constructed
   at a fixed distance from the nitrogen atom, with a fixed angle to the
   carbon atom, in the plane defined by one of the other heavy atoms,
   see Fig. [fig:vsitehydro]C.

-  *amine group (umbrella -NH:math:`_2` or -NH:math:`_3^+`) hydrogens:*
   Amine hydrogens with rotational freedom cannot be constructed as
   virtual sites from the heavy atoms they are connected to, since this
   would result in loss of the rotational freedom of the amine group. To
   preserve the rotational freedom while removing the hydrogen
   bond-angle degrees of freedom, two “dummy masses” are constructed
   with the same total mass, moment of inertia (for rotation around the
   C-N bond) and center of mass as the amine group. These dummy masses
   have no interaction with any other atom, except for the fact that
   they are connected to the carbon and to each other, resulting in a
   rigid triangle. From these three particles, the positions of the
   nitrogen and hydrogen atoms are constructed as linear combinations of
   the two carbon-mass vectors and their outer product, resulting in an
   amine group with rotational freedom intact, but without other
   internal degrees of freedom. See Fig. [fig:vsitehydro]D.

.. figure:: plots/dumaro
   :alt: The different types of virtual site constructions used for
   aromatic residues. The atoms used in the construction of the virtual
   site(s) are depicted as black circles, virtual sites as gray ones.
   Hydrogens are smaller than heavy atoms. A: phenylalanine; B: tyrosine
   (note that the hydroxyl hydrogen is *not* a virtual site); C:
   tryptophan; D: histidine.
   :width: 15.00000cm

   The different types of virtual site constructions used for aromatic
   residues. The atoms used in the construction of the virtual site(s)
   are depicted as black circles, virtual sites as gray ones. Hydrogens
   are smaller than heavy atoms. A: phenylalanine; B: tyrosine (note
   that the hydroxyl hydrogen is *not* a virtual site); C: tryptophan;
   D: histidine.

Out-of-plane vibrations in aromatic groups
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The planar arrangements in the side chains of the aromatic residues
lends itself perfectly to a virtual-site construction, giving a
perfectly planar group without the inherently unstable constraints that
are necessary to keep normal atoms in a plane. The basic approach is to
define three atoms or dummy masses with constraints between them to fix
the geometry and create the rest of the atoms as simple virtual sites
type (see sec. [sec:virtual\_sites]) from these three. Each of the
aromatic residues require a different approach:

-  *Phenylalanine:* C\ :math:`_\gamma`, C\ :math:`_{{\epsilon}1}`, and
   C\ :math:`_{{\epsilon}2}` are kept as normal atoms, but with each a
   mass of one third the total mass of the phenyl group. See
   Fig. [fig:vsitehydro]A.

-  *Tyrosine:* The ring is treated identically to the phenylalanine
   ring. Additionally, constraints are defined between
   C\ :math:`_{{\epsilon}1}`, C\ :math:`_{{\epsilon}2}`, and
   O\ :math:`_{\eta}`. The original improper dihedral angles will keep
   both triangles (one for the ring and one with O\ :math:`_{\eta}`) in
   a plane, but due to the larger moments of inertia this construction
   will be much more stable. The bond-angle in the hydroxyl group will
   be constrained by a constraint between C\ :math:`_\gamma` and
   H\ :math:`_{\eta}`. **Note** that the hydrogen is not treated as a
   virtual site. See Fig. [fig:vsitehydro]B.

-  *Tryptophan:* C\ :math:`_\beta` is kept as a normal atom and two
   dummy masses are created at the center of mass of each of the rings,
   each with a mass equal to the total mass of the respective ring
   (C\ :math:`_{{\delta}2}` and C\ :math:`_{{\epsilon}2}` are each
   counted half for each ring). This keeps the overall center of mass
   and the moment of inertia almost (but not quite) equal to what it
   was. See Fig. [fig:vsitehydro]C.

-  *Histidine:* C\ :math:`_\gamma`, C\ :math:`_{{\epsilon}1}` and
   N\ :math:`_{{\epsilon}2}` are kept as normal atoms, but with masses
   redistributed such that the center of mass of the ring is preserved.
   See Fig. [fig:vsitehydro]D.

Viscosity calculation
---------------------

The shear viscosity is a property of liquids that can be determined
easily by experiment. It is useful for parameterizing a force field
because it is a kinetic property, while most other properties which are
used for parameterization are thermodynamic. The viscosity is also an
important property, since it influences the rates of conformational
changes of molecules solvated in the liquid.

The viscosity can be calculated from an equilibrium simulation using an
Einstein relation:

.. math::

   \eta = \frac{1}{2}\frac{V}{k_B T} \lim_{t \rightarrow \infty}
   \frac{\mbox{d}}{\mbox{d} t} \left\langle 
   \left( \int_{t_0}^{{t_0}+t} P_{xz}(t') \mbox{d} t' \right)^2
   \right\rangle_{t_0}

 This can be done with gmx energy. This method converges very
slowly \ `149 <#ref-Hess2002a>`__, and as such a nanosecond simulation
might not be long enough for an accurate determination of the viscosity.
The result is very dependent on the treatment of the electrostatics.
Using a (short) cut-off results in large noise on the off-diagonal
pressure elements, which can increase the calculated viscosity by an
order of magnitude.

GROMACS also has a non-equilibrium method for determining the
viscosity \ `149 <#ref-Hess2002a>`__. This makes use of the fact that
energy, which is fed into system by external forces, is dissipated
through viscous friction. The generated heat is removed by coupling to a
heat bath. For a Newtonian liquid adding a small force will result in a
velocity gradient according to the following equation:

.. math:: a_x(z) + \frac{\eta}{\rho} \frac{\partial^2 v_x(z)}{\partial z^2} = 0

 Here we have applied an acceleration :math:`a_x(z)` in the
:math:`x`-direction, which is a function of the :math:`z`-coordinate. In
GROMACS the acceleration profile is:

.. math:: a_x(z) = A \cos\left(\frac{2\pi z}{l_z}\right)

 where :math:`l_z` is the height of the box. The generated velocity
profile is:

.. math:: v_x(z) = V \cos\left(\frac{2\pi z}{l_z}\right)

.. math:: V = A \frac{\rho}{\eta}\left(\frac{l_z}{2\pi}\right)^2

 The viscosity can be calculated from :math:`A` and :math:`V`:

.. math::

   \label{visc}
   \eta = \frac{A}{V}\rho \left(\frac{l_z}{2\pi}\right)^2

In the simulation :math:`V` is defined as:

.. math::

   V = \frac{\displaystyle \sum_{i=1}^N m_i v_{i,x} 2 \cos\left(\frac{2\pi z}{l_z}\right)}
            {\displaystyle \sum_{i=1}^N m_i}

 The generated velocity profile is not coupled to the heat bath.
Moreover, the velocity profile is excluded from the kinetic energy. One
would like :math:`V` to be as large as possible to get good statistics.
However, the shear rate should not be so high that the system gets too
far from equilibrium. The maximum shear rate occurs where the cosine is
zero, the rate being:

.. math::

   \mbox{sh}_{\max} =  \max_z \left| \frac{\partial v_x(z)}{\partial z} \right|
   = A \frac{\rho}{\eta} \frac{l_z}{2\pi}

 For a simulation with: :math:`\eta=10^{-3}`
[kgm:math:`^{-1}`\ s\ :math:`^{-1}`],
:math:`\rho=10^3`\ [kgm:math:`^{-3}`] and :math:`l_z=2\pi`\ [nm],
:math:`\mbox{sh}_{\max}=1`\ [psnm:math:`^{-1}`] :math:`A`. This shear
rate should be smaller than one over the longest correlation time in the
system. For most liquids, this will be the rotation correlation time,
which is around 10 ps. In this case, :math:`A` should be smaller than
0.1[nmps\ :math:`^{-2}`]. When the shear rate is too high, the observed
viscosity will be too low. Because :math:`V` is proportional to the
square of the box height, the optimal box is elongated in the
:math:`z`-direction. In general, a simulation length of 100 ps is enough
to obtain an accurate value for the viscosity.

The heat generated by the viscous friction is removed by coupling to a
heat bath. Because this coupling is not instantaneous the real
temperature of the liquid will be slightly lower than the observed
temperature. Berendsen derived this temperature
shift \ `31 <#ref-Berendsen91>`__, which can be written in terms of the
shear rate as:

.. math:: T_s = \frac{\eta\,\tau}{2 \rho\,C_v} \mbox{sh}_{\max}^2

 where :math:`\tau` is the coupling time for the Berendsen thermostat
and :math:`C_v` is the heat capacity. Using the values of the example
above, :math:`\tau=10^{-13}` [s] and :math:`C_v=2 \cdot 10^3`\ [J
kg\ :math:`^{-1}`\ K\ :math:`^{-1}`], we get:
:math:`T_s=25`\ [Kps:math:`^{-2}`]sh\ :math:`_{\max}^2`. When we want
the shear rate to be smaller than :math:`1/10`\ [ps:math:`^{-1}`],
:math:`T_s` is smaller than 0.25[K], which is negligible.

**Note** that the system has to build up the velocity profile when
starting from an equilibrium state. This build-up time is of the order
of the correlation time of the liquid.

Two quantities are written to the energy file, along with their averages
and fluctuations: :math:`V` and :math:`1/\eta`, as obtained from
([visc]).

Tabulated interaction functions
-------------------------------

Cubic splines for potentials
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In some of the inner loops of GROMACS, look-up tables are used for
computation of potential and forces. The tables are interpolated using a
cubic spline algorithm. There are separate tables for electrostatic,
dispersion, and repulsion interactions, but for the sake of caching
performance these have been combined into a single array. The cubic
spline interpolation for :math:`x_i \leq x < x_{i+1}` looks like this:

.. math::

   V_s(x) = A_0 + A_1 \,\epsilon + A_2 \,\epsilon^2 + A_3 \,\epsilon^3
   \label{eqn:spline}

 where the table spacing :math:`h` and fraction :math:`\epsilon` are
given by:

.. math::

   \begin{aligned}
   h	&=&	x_{i+1} - x_i	\\
   \epsilon&=&	(x - x_i)/h\end{aligned}

 so that :math:`0 \le \epsilon < 1`. From this, we can calculate the
derivative in order to determine the forces:

.. math::

   -V_s'(x) ~=~ 
   -\frac{{\rm d}V_s(x)}{{\rm d}\epsilon}\frac{{\rm d}\epsilon}{{\rm d}x} ~=~
   -(A_1 + 2 A_2 \,\epsilon + 3 A_3 \,\epsilon^2)/h

 The four coefficients are determined from the four conditions that
:math:`V_s` and :math:`-V_s'` at both ends of each interval should match
the exact potential :math:`V` and force :math:`-V'`. This results in the
following errors for each interval:

.. math::

   \begin{aligned}
   |V_s  - V  |_{max} &=& V'''' \frac{h^4}{384} + O(h^5) \\
   |V_s' - V' |_{max} &=& V'''' \frac{h^3}{72\sqrt{3}} + O(h^4) \\
   |V_s''- V''|_{max} &=& V'''' \frac{h^2}{12}  + O(h^3)\end{aligned}

 V and V’ are continuous, while V” is the first discontinuous
derivative. The number of points per nanometer is 500 and 2000 for
mixed- and double-precision versions of GROMACS, respectively. This
means that the errors in the potential and force will usually be smaller
than the mixed precision accuracy.

GROMACS stores :math:`A_0`, :math:`A_1`, :math:`A_2` and :math:`A_3`.
The force routines get a table with these four parameters and a scaling
factor :math:`s` that is equal to the number of points per nm. (**Note**
that :math:`h` is :math:`s^{-1}`). The algorithm goes a little something
like this:

#. Calculate distance vector (:math:`_{ij}`) and distance
   r\ :math:`_{ij}`

#. Multiply r\ :math:`_{ij}` by :math:`s` and truncate to an integer
   value :math:`n_0` to get a table index

#. Calculate fractional component (:math:`\epsilon` =
   :math:`s`\ r\ :math:`_{ij} - n_0`) and :math:`\epsilon^2`

#. Do the interpolation to calculate the potential :math:`V` and the
   scalar force :math:`f`

#. Calculate the vector force by multiplying :math:`f` with
   :math:`_{ij}`

**Note** that table look-up is significantly *slower* than computation
of the most simple Lennard-Jones and Coulomb interaction. However, it is
much faster than the shifted Coulomb function used in conjunction with
the PPPM method. Finally, it is much easier to modify a table for the
potential (and get a graphical representation of it) than to modify the
inner loops of the MD program.

User-specified potential functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can also use your own potential functions without editing the
GROMACS code. The potential function should be according to the
following equation

.. math:: V(r_{ij}) ~=~ \frac{q_i q_j}{4 \pi\epsilon_0} f(r_{ij}) + C_6 \,g(r_{ij}) + C_{12} \,h(r_{ij})

 where :math:`f`, :math:`g`, and :math:`h` are user defined functions.
**Note** that if :math:`g(r)` represents a normal dispersion
interaction, :math:`g(r)` should be :math:`<` 0. C\ :math:`_6`,
C\ :math:`_{12}` and the charges are read from the topology. Also note
that combination rules are only supported for Lennard-Jones and
Buckingham, and that your tables should match the parameters in the
binary topology.

When you add the following lines in your .mdp file:

::

    rlist           = 1.0
    coulombtype     = User
    rcoulomb        = 1.0
    vdwtype         = User
    rvdw            = 1.0

mdrun will read a single non-bonded table file, or multiple when
energygrp-table is set (see below). The name of the file(s) can be set
with the mdrun option -table. The table file should contain seven
columns of table look-up data in the order: :math:`x`, :math:`f(x)`,
:math:`-f'(x)`, :math:`g(x)`, :math:`-g'(x)`, :math:`h(x)`,
:math:`-h'(x)`. The :math:`x` should run from 0 to :math:`r_c+1` (the
value of table\_extension can be changed in the .mdp file). You can
choose the spacing you like; for the standard tables GROMACS uses a
spacing of 0.002 and 0.0005 nm when you run in mixed and double
precision, respectively. In this context, :math:`r_c` denotes the
maximum of the two cut-offs rvdw and rcoulomb (see above). These
variables need not be the same (and need not be 1.0 either). Some
functions used for potentials contain a singularity at :math:`x = 0`,
but since atoms are normally not closer to each other than 0.1 nm, the
function value at :math:`x = 0` is not important. Finally, it is also
possible to combine a standard Coulomb with a modified LJ potential (or
vice versa). One then specifies *e.g.* coulombtype = Cut-off or
coulombtype = PME, combined with vdwtype = User. The table file must
always contain the 7 columns however, and meaningful data (i.e. not
zeroes) must be entered in all columns. A number of pre-built table
files can be found in the GMXLIB directory for 6-8, 6-9, 6-10, 6-11, and
6-12 Lennard-Jones potentials combined with a normal Coulomb.

If you want to have different functional forms between different groups
of atoms, this can be set through energy groups. Different tables can be
used for non-bonded interactions between different energy groups pairs
through the .mdp option energygrp-table (see details in the User Guide).
Atoms that should interact with a different potential should be put into
different energy groups. Between group pairs which are not listed in
energygrp-table, the normal user tables will be used. This makes it easy
to use a different functional form between a few types of atoms.

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
The current version of GROMACS provides interfaces to several popular
Quantum Chemistry packages (MOPAC `150 <#ref-mopac>`__,
GAMESS-UK \ `151 <#ref-gamess-uk>`__, Gaussian \ `152 <#ref-g03>`__ and
CPMD \ `153 <#ref-Car85a>`__).

GROMACS interactions between the two subsystems are either handled as
described by Field *et al.* `154 <#ref-Field90a>`__ or within the ONIOM
approach by Morokuma and coworkers \ `155 <#ref-Maseras96a>`__,
`156 <#ref-Svensson96a>`__.

Overview
~~~~~~~~

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
~~~~~

To make use of the QM/MM functionality in GROMACS, one needs to:

#. introduce link atoms at the QM/MM boundary, if needed;

#. specify which atoms are to be treated at a QM level;

#. specify the QM level, basis set, type of QM/MM interface and so on.

Adding link atoms
^^^^^^^^^^^^^^^^^

At the bond that connects the QM and MM subsystems, a link atoms is
introduced. In GROMACS the link atom has special atomtype, called LA.
This atomtype is treated as a hydrogen atom in the QM calculation, and
as a virtual site in the force-field calculation. The link atoms, if
any, are part of the system, but have no interaction with any other
atom, except that the QM force working on it is distributed over the two
atoms of the bond. In the topology, the link atom (LA), therefore, is
defined as a virtual site atom:

::

    [ virtual_sites2 ]
    LA QMatom MMatom 1 0.65

See sec. [sec:vsitetop] for more details on how virtual sites are
treated. The link atom is replaced at every step of the simulation.

In addition, the bond itself is replaced by a constraint:

::

    [ constraints ]
    QMatom MMatom 2 0.153

**Note** that, because in our system the QM/MM bond is a carbon-carbon
bond (0.153 nm), we use a constraint length of 0.153 nm, and dummy
position of 0.65. The latter is the ratio between the ideal C-H bond
length and the ideal C-C bond length. With this ratio, the link atom is
always 0.1 nm away from the QMatom, consistent with the carbon-hydrogen
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

In the .mdp file, the following parameters control a QM/MM simulation.

QMMM = no
    | 
    | If this is set to yes, a QM/MM simulation is requested. Several
      groups of atoms can be described at different QM levels
      separately. These are specified in the QMMM-grps field separated
      by spaces. The level of *ab initio* theory at which the groups are
      described is specified by QMmethod and QMbasis Fields. Describing
      the groups at different levels of theory is only possible with the
      ONIOM QM/MM scheme, specified by QMMMscheme.

QMMM-grps =
    | 
    | groups to be described at the QM level

QMMMscheme = normal
    | 
    | Options are normal and ONIOM. This selects the QM/MM interface.
      normal implies that the QM subsystem is electronically embedded in
      the MM subsystem. There can only be one QMMM-grps that is modeled
      at the QMmethod and QMbasis level of * ab initio* theory. The rest
      of the system is described at the MM level. The QM and MM
      subsystems interact as follows: MM point charges are included in
      the QM one-electron Hamiltonian and all Lennard-Jones interactions
      are described at the MM level. If ONIOM is selected, the
      interaction between the subsystem is described using the ONIOM
      method by Morokuma and co-workers. There can be more than one
      QMMM-grps each modeled at a different level of QM theory (QMmethod
      and QMbasis).

QMmethod =
    | 
    | Method used to compute the energy and gradients on the QM atoms.
      Available methods are AM1, PM3, RHF, UHF, DFT, B3LYP, MP2, CASSCF,
      MMVB and CPMD. For CASSCF, the number of electrons and orbitals
      included in the active space is specified by CASelectrons and
      CASorbitals. For CPMD, the plane-wave cut-off is specified by the
      planewavecutoff keyword.

QMbasis =
    | 
    | Gaussian basis set used to expand the electronic wave-function.
      Only Gaussian basis sets are currently available, i.e. STO-3G,
      3-21G, 3-21G\*, 3-21+G\*, 6-21G, 6-31G, 6-31G\*, 6-31+G\*, and
      6-311G. For CPMD, which uses plane wave expansion rather than
      atom-centered basis functions, the planewavecutoff keyword
      controls the plane wave expansion.

QMcharge =
    | 
    | The total charge in *e* of the QMMM-grps. In case there are more
      than one QMMM-grps, the total charge of each ONIOM layer needs to
      be specified separately.

QMmult =
    | 
    | The multiplicity of the QMMM-grps. In case there are more than one
      QMMM-grps, the multiplicity of each ONIOM layer needs to be
      specified separately.

CASorbitals =
    | 
    | The number of orbitals to be included in the active space when
      doing a CASSCF computation.

CASelectrons =
    | 
    | The number of electrons to be included in the active space when
      doing a CASSCF computation.

SH = no
    | 
    | If this is set to yes, a QM/MM MD simulation on the excited
      state-potential energy surface and enforce a diabatic hop to the
      ground-state when the system hits the conical intersection
      hyperline in the course the simulation. This option only works in
      combination with the CASSCF method.

Output
~~~~~~

The energies and gradients computed in the QM calculation are added to
those computed by GROMACS. In the .edr file there is a section for the
total QM energy.

Future developments
~~~~~~~~~~~~~~~~~~~

Several features are currently under development to increase the
accuracy of the QM/MM interface. One useful feature is the use of
delocalized MM charges in the QM computations. The most important
benefit of using such smeared-out charges is that the Coulombic
potential has a finite value at interatomic distances. In the point
charge representation, the partially-charged MM atoms close to the QM
region tend to “over-polarize” the QM system, which leads to artifacts
in the calculation.

What is needed as well is a transition state optimizer.

Using VMD plug-ins for trajectory file I/O
------------------------------------------

tools are able to use the plug-ins found in an existing installation of
`VMD <http://www.ks.uiuc.edu/Research/vmd>`__ in order to read and write
trajectory files in formats that are not native to GROMACS. You will be
able to supply an AMBER DCD-format trajectory filename directly to
GROMACS tools, for example.

This requires a VMD installation not older than version 1.8, that your
system provides the dlopen function so that programs can determine at
run time what plug-ins exist, and that you build shared libraries when
building GROMACS. CMake will find the vmd executable in your path, and
from it, or the environment variable VMDDIR at configuration or run
time, locate the plug-ins. Alternatively, the VMD\_PLUGIN\_PATH can be
used at run time to specify a path where these plug-ins can be found.
Note that these plug-ins are in a binary format, and that format must
match the architecture of the machine attempting to use them.

Interactive Molecular Dynamics
------------------------------

GROMACS supports the interactive molecular dynamics (IMD) protocol as
implemented by `VMD <http://www.ks.uiuc.edu/Research/vmd>`__ to control
a running simulation in NAMD. IMD allows to monitor a running GROMACS
simulation from a VMD client. In addition, the user can interact with
the simulation by pulling on atoms, residues or fragments with a mouse
or a force-feedback device. Additional information about the GROMACS
implementation and an exemplary GROMACS IMD system can be found `on this
homepage <http://www.mpibpc.mpg.de/grubmueller/interactivemd>`__.

Simulation input preparation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The GROMACS implementation allows transmission and interaction with a
part of the running simulation only, e.g. in cases where no water
molecules should be transmitted or pulled. The group is specified via
the .mdp option IMD-group. When IMD-group is empty, the IMD protocol is
disabled and cannot be enabled via the switches in mdrun. To interact
with the entire system, IMD-group can be set to System. When using
grompp, a .gro file to be used as VMD input is written out (-imd switch
of grompp).

Starting the simulation
~~~~~~~~~~~~~~~~~~~~~~~

Communication between VMD and GROMACS is achieved via TCP sockets and
thus enables controlling an mdrun running locally or on a remote
cluster. The port for the connection can be specified with the -imdport
switch of mdrun, 8888 is the default. If a port number of 0 or smaller
is provided, GROMACS automatically assigns a free port to use with IMD.

Every :math:`N` steps, the mdrun client receives the applied forces from
VMD and sends the new positions to the client. VMD permits increasing or
decreasing the communication frequency interactively. By default, the
simulation starts and runs even if no IMD client is connected. This
behavior is changed by the -imdwait switch of mdrun. After startup and
whenever the client has disconnected, the integration stops until
reconnection of the client. When the -imdterm switch is used, the
simulation can be terminated by pressing the stop button in VMD. This is
disabled by default. Finally, to allow interacting with the simulation
(i.e. pulling from VMD) the -imdpull switch has to be used. Therefore, a
simulation can only be monitored but not influenced from the VMD client
when none of -imdwait, -imdterm or -imdpull are set. However, since the
IMD protocol requires no authentication, it is not advisable to run
simulations on a host directly reachable from an insecure environment.
Secure shell forwarding of TCP can be used to connect to running
simulations not directly reachable from the interacting host. Note that
the IMD command line switches of mdrun are hidden by default and show up
in the help text only with gmx mdrun -h -hidden.

Connecting from VMD
~~~~~~~~~~~~~~~~~~~

In VMD, first the structure corresponding to the IMD group has to be
loaded (*File :math:`\rightarrow` New Molecule*). Then the IMD
connection window has to be used (*Extensions :math:`\rightarrow`
Simulation :math:`\rightarrow` IMD Connect (NAMD)*). In the IMD
connection window, hostname and port have to be specified and followed
by pressing *Connect*. *Detach Sim* allows disconnecting without
terminating the simulation, while *Stop Sim* ends the simulation on the
next neighbor searching step (if allowed by -imdterm).

The timestep transfer rate allows adjusting the communication frequency
between simulation and IMD client. Setting the keep rate loads every
:math:`N^\mathrm{th}` frame into VMD instead of discarding them when a
new one is received. The displayed energies are in SI units in contrast
to energies displayed from NAMD simulations.

Embedding proteins into the membranes
-------------------------------------

GROMACS is capable of inserting the protein into pre-equilibrated lipid
bilayers with minimal perturbation of the lipids using the method, which
was initially described as a ProtSqueeze
technique,\ `157 <#ref-Yesylevskyy2007>`__ and later implemented as
g\_membed tool.\ `158 <#ref-Wolf2010>`__ Currently the functionality of
g\_membed is available in mdrun as described in the user guide.

This method works by first artificially shrinking the protein in the
:math:`xy`-plane, then it removes lipids that overlap with that much
smaller core. Then the protein atoms are gradually resized back to their
initial configuration, using normal dynamics for the rest of the system,
so the lipids adapt to the protein. Further lipids are removed as
required.

Run parameters and Programs
===========================

Online documentation
--------------------

More documentation is available online from the GROMACS web site,
http://manual.gromacs.org/documentation.

In addition, we install standard UNIX man pages for all the programs. If
you have sourced the GMXRC script in the GROMACS binary directory for
your host they should already be present in your MANPATH environment
variable, and you should be able to type *e.g.* man gmx-grompp. You can
also use the -h flag on the command line (e.g. gmx grompp -h) to see the
same information, as well as gmx help grompp. The list of all programs
are available from gmx help.

File types
----------

Table [tab:form] lists the file types used by GROMACS along with a short
description, and you can find a more detail description for each file in
your HTML reference, or in our online version.

GROMACS files written in XDR format can be read on any architecture with
GROMACS version 1.6 or later if the configuration script found the XDR
libraries on your system. They should always be present on UNIX since
they are necessary for NFS support.

Run Parameters
--------------

The descriptions of .mdp parameters can be found at
http://manual.gromacs.org/current/mdp-options.html or in your
installation at share/gromacs/html/mdp-options.html

Analysis
========

In this chapter different ways of analyzing your trajectory are
described. The names of the corresponding analysis programs are given.
Specific information on the in- and output of these programs can be
found in the online manual at
`www.gromacs.org <http://www.gromacs.org>`__. The output files are often
produced as finished Grace/Xmgr graphs.

First, in sec. [sec:usinggroups], the group concept in analysis is
explained. [subsec:selections] explains a newer concept of dynamic
selections, which is currently supported by a few tools. Then, the
different analysis tools are presented.

Using Groups
------------

| 
| In chapter [ch:algorithms], it was explained how *groups of atoms* can
  be used in mdrun (see sec. [sec:groupconcept]). In most analysis
  programs, groups of atoms must also be chosen. Most programs can
  generate several default index groups, but groups can always be read
  from an index file. Let’s consider the example of a simulation of a
  binary mixture of components A and B. When we want to calculate the
  radial distribution function (RDF) :math:`g_{AB}(r)` of A with respect
  to B, we have to calculate:

  .. math:: 4\pi r^2 g_{AB}(r)      ~=~     V~\sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} P(r)

   where :math:`V` is the volume and :math:`P(r)` is the probability of
  finding a B atom at distance :math:`r` from an A atom.

By having the user define the *atom numbers* for groups A and B in a
simple file, we can calculate this :math:`g_{AB}` in the most general
way, without having to make any assumptions in the RDF program about the
type of particles.

Groups can therefore consist of a series of *atom numbers*, but in some
cases also of *molecule numbers*. It is also possible to specify a
series of angles by *triples* of *atom numbers*, dihedrals by
*quadruples* of *atom numbers* and bonds or vectors (in a molecule) by
*pairs* of *atom numbers*. When appropriate the type of index file will
be specified for the following analysis programs. To help creating such
index files (index.ndx), there are a couple of programs to generate
them, using either your input configuration or the topology. To generate
an index file consisting of a series of *atom numbers* (as in the
example of :math:`g_{AB}`), use gmx make\_ndx or gmx select. To generate
an index file with angles or dihedrals, use gmx mk\_angndx. Of course
you can also make them by hand. The general format is presented here:

::

    [ Oxygen ]
       1       4       7

    [ Hydrogen ]
       2       3       5       6
       8       9

First, the group name is written between square brackets. The following
atom numbers may be spread out over as many lines as you like. The atom
numbering starts at 1.

Each tool that can use groups will offer the available alternatives for
the user to choose. That choice can be made with the number of the
group, or its name. In fact, the first few letters of the group name
will suffice if that will distinguish the group from all others. There
are ways to use Unix shell features to choose group names on the command
line, rather than interactively. Consult
`www.gromacs.org <http://www.gromacs.org>`__ for suggestions.

Default Groups
~~~~~~~~~~~~~~

When no index file is supplied to analysis tools or grompp, a number of
default groups are generated to choose from:

System
    | 
    | all atoms in the system

Protein
    | 
    | all protein atoms

Protein-H
    | 
    | protein atoms excluding hydrogens

C-alpha
    | 
    | C\ :math:`_{\alpha}` atoms

Backbone
    | 
    | protein backbone atoms; N, C\ :math:`_{\alpha}` and C

MainChain
    | 
    | protein main chain atoms: N, C\ :math:`_{\alpha}`, C and O,
      including oxygens in C-terminus

MainChain+Cb
    | 
    | protein main chain atoms including C\ :math:`_{\beta}`

MainChain+H
    | 
    | protein main chain atoms including backbone amide hydrogens and
      hydrogens on the N-terminus

SideChain
    | 
    | protein side chain atoms; that is all atoms except N,
      C\ :math:`_{\alpha}`, C, O, backbone amide hydrogens, oxygens in
      C-terminus and hydrogens on the N-terminus

SideChain-H
    | 
    | protein side chain atoms excluding all hydrogens

Prot-Masses
    | 
    | protein atoms excluding dummy masses (as used in virtual site
      constructions of NH\ :math:`_3` groups and tryptophan
      side-chains), see also sec. [sec:vsitetop]; this group is only
      included when it differs from the “Protein” group

Non-Protein
    | 
    | all non-protein atoms

DNA
    | 
    | all DNA atoms

RNA
    | 
    | all RNA atoms

Water
    | 
    | water molecules (names like SOL, WAT, HOH, etc.) See
      residuetypes.dat for a full listing

non-Water
    | 
    | anything not covered by the Water group

Ion
    | 
    | any name matching an Ion entry in residuetypes.dat

Water\_and\_Ions
    | 
    | combination of the Water and Ions groups

molecule\_name
    | 
    | for all residues/molecules which are not recognized as protein,
      DNA, or RNA; one group per residue/molecule name is generated

Other
    | 
    | all atoms which are neither protein, DNA, nor RNA.

Empty groups will not be generated. Most of the groups only contain
protein atoms. An atom is considered a protein atom if its residue name
is listed in the residuetypes.dat file and is listed as a “Protein”
entry. The process for determinding DNA, RNA, etc. is analogous. If you
need to modify these classifications, then you can copy the file from
the library directory into your working directory and edit the local
copy.

Selections
~~~~~~~~~~

| gmx select
| Currently, a few analysis tools support an extended concept of
  *(dynamic) selections*. There are three main differences to
  traditional index groups:

-  The selections are specified as text instead of reading fixed atom
   indices from a file, using a syntax similar to VMD. The text can be
   entered interactively, provided on the command line, or from a file.

-  The selections are not restricted to atoms, but can also specify that
   the analysis is to be performed on, e.g., center-of-mass positions of
   a group of atoms. Some tools may not support selections that do not
   evaluate to single atoms, e.g., if they require information that is
   available only for single atoms, like atom names or types.

-  The selections can be dynamic, i.e., evaluate to different atoms for
   different trajectory frames. This allows analyzing only a subset of
   the system that satisfies some geometric criteria.

As an example of a simple selection, resname ABC and within 2 of resname
DEF selects all atoms in residues named ABC that are within 2nm of any
atom in a residue named DEF.

Tools that accept selections can also use traditional index files
similarly to older tools: it is possible to give an .ndx file to the
tool, and directly select a group from the index file as a selection,
either by group number or by group name. The index groups can also be
used as a part of a more complicated selection.

To get started, you can run gmx select with a single structure, and use
the interactive prompt to try out different selections. The tool
provides, among others, output options -on and -ofpdb to write out the
selected atoms to an index file and to a .pdb file, respectively. This
does not allow testing selections that evaluate to center-of-mass
positions, but other selections can be tested and the result examined.

The detailed syntax and the individual keywords that can be used in
selections can be accessed by typing help in the interactive prompt of
any selection-enabled tool, as well as with gmx help selections. The
help is divided into subtopics that can be accessed with, e.g., help
syntax / gmx help selections syntax. Some individual selection keywords
have extended help as well, which can be accessed with, e.g., help
keywords within.

The interactive prompt does not currently provide much editing
capabilities. If you need them, you can run the program under rlwrap.

For tools that do not yet support the selection syntax, you can use gmx
select -on to generate static index groups to pass to the tool. However,
this only allows for a small subset (only the first bullet from the
above list) of the flexibility that fully selection-aware tools offer.

It is also possible to write your own analysis tools to take advantage
of the flexibility of these selections: see the template.cpp file in the
share/gromacs/template directory of your installation for an example.

Looking at your trajectory
--------------------------

|The window of gmx view showing a box of water.|

| gmx view
| Before analyzing your trajectory it is often informative to look at
  your trajectory first. GROMACS comes with a simple trajectory viewer
  gmx view; the advantage with this one is that it does not require
  OpenGL, which usually isn’t present on *e.g.* supercomputers. It is
  also possible to generate a hard-copy in Encapsulated Postscript
  format (see Fig. [fig:ngmxdump]). If you want a faster and more fancy
  viewer there are several programs that can read the GROMACS trajectory
  formats – have a look at our homepage
  (`www.gromacs.org <http://www.gromacs.org>`__) for updated links.

General properties
------------------

| gmx energy, gmx traj
| To analyze some or all *energies* and other properties, such as *total
  pressure*, *pressure tensor*, *density*, *box-volume* and *box-sizes*,
  use the program gmx energy. A choice can be made from a list a set of
  energies, like potential, kinetic or total energy, or individual
  contributions, like Lennard-Jones or dihedral energies.

The *center-of-mass velocity*, defined as

.. math:: {\bf v}_{com} = {1 \over M} \sum_{i=1}^N m_i {\bf v}_i

 with :math:`M = \sum_{i=1}^N m_i` the total mass of the system, can be
monitored in time by the program gmx traj -com -ov. It is however
recommended to remove the center-of-mass velocity every step (see
chapter [ch:algorithms])!

Radial distribution functions
-----------------------------

| gmx rdf
| The *radial distribution function* (RDF) or pair correlation function
  :math:`g_{AB}(r)` between particles of type :math:`A` and :math:`B` is
  defined in the following way:

.. math::

   \begin{array}{rcl}
   g_{AB}(r)&=&    {\displaystyle \frac{\langle \rho_B(r) \rangle}{\langle\rho_B\rangle_{local}}}         \\
            &=&    {\displaystyle \frac{1}{\langle\rho_B\rangle_{local}}}{\displaystyle \frac{1}{N_A}}
                   \sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} 
                   {\displaystyle \frac{\delta( r_{ij} - r )}{4 \pi r^2}}         \\
   \end{array}

with :math:`\langle\rho_B(r)\rangle` the particle density of type
:math:`B` at a distance :math:`r` around particles :math:`A`, and
:math:`\langle\rho_B\rangle_{local}` the particle density of type
:math:`B` averaged over all spheres around particles :math:`A` with
radius :math:`r_{max}` (see Fig. [fig:rdfex]C).

|Definition of slices in gmx rdf: A. :math:`g_{AB}(r)`. B.
:math:`g_{AB}(r,\theta)`. The slices are colored gray. C. Normalization
:math:`\langle\rho_B\rangle_{local}`. D. Normalization
:math:`\langle\rho_B\rangle_{local,\:\theta }`. Normalization volumes
are colored gray.|

Usually the value of :math:`r_{max}` is half of the box length. The
averaging is also performed in time. In practice the analysis program
gmx rdf divides the system into spherical slices (from :math:`r` to
:math:`r+dr`, see Fig. [fig:rdfex]A) and makes a histogram in stead of
the :math:`\delta`-function. An example of the RDF of oxygen-oxygen in
SPC water \ `80 <#ref-Berendsen81>`__ is given in Fig. [fig:rdf].

|:math:`g_{OO}(r)` for Oxygen-Oxygen of SPC-water.|

With gmx rdf it is also possible to calculate an angle dependent rdf
:math:`g_{AB}(r,\theta)`, where the angle :math:`\theta` is defined with
respect to a certain laboratory axis :math:`{\bf e}`, see
Fig. [fig:rdfex]B.

.. math::

   \begin{aligned}
   g_{AB}(r,\theta) &=& {1 \over \langle\rho_B\rangle_{local,\:\theta }} {1 \over N_A} \sum_{i \in A}^{N_A} \sum_{j \in B}^{N_B} {\delta( r_{ij} - r ) \delta(\theta_{ij} -\theta) \over 2 \pi r^2 sin(\theta)}\\
   cos(\theta_{ij}) &=& {{\bf r}_{ij} \cdot {\bf e} \over \|r_{ij}\| \;\| e\| }\end{aligned}

 This :math:`g_{AB}(r,\theta)` is useful for analyzing anisotropic
systems. **Note** that in this case the normalization
:math:`\langle\rho_B\rangle_{local,\:\theta}` is the average density in
all angle slices from :math:`\theta` to :math:`\theta + d\theta` up to
:math:`r_{max}`, so angle dependent, see Fig. [fig:rdfex]D.

Correlation functions
---------------------

Theory of correlation functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The theory of correlation functions is well
established \ `108 <#ref-Allen87>`__. We describe here the
implementation of the various correlation function flavors in the
GROMACS code. The definition of the autocorrelation function (ACF)
:math:`C_f(t)` for a property :math:`f(t)` is:

.. math::

   C_f(t)  ~=~     \left\langle f(\xi) f(\xi+t)\right\rangle_{\xi}
   \label{eqn:corr}

 where the notation on the right hand side indicates averaging over
:math:`\xi`, *i.e.* over time origins. It is also possible to compute
cross-correlation function from two properties :math:`f(t)` and
:math:`g(t)`:

.. math:: C_{fg}(t) ~=~   \left\langle f(\xi) g(\xi+t)\right\rangle_{\xi}

 however, in GROMACS there is no standard mechanism to do this
(**note:** you can use the xmgr program to compute cross correlations).
The integral of the correlation function over time is the correlation
time :math:`\tau_f`:

.. math::

   \tau_f  ~=~     \int_0^{\infty} C_f(t) {\rm d} t
   \label{eqn:corrtime}

In practice, correlation functions are calculated based on data points
with discrete time intervals :math:`\Delta`\ t, so that the ACF from an
MD simulation is:

.. math::

   C_f(j\Delta t)  ~=~     \frac{1}{N-j}\sum_{i=0}^{N-1-j} f(i\Delta t) f((i+j)\Delta t)
   \label{eqn:corrmd}

 where :math:`N` is the number of available time frames for the
calculation. The resulting ACF is obviously only available at time
points with the same interval :math:`\Delta`\ t. Since, for many
applications, it is necessary to know the short time behavior of the ACF
(*e.g.* the first 10 ps) this often means that we have to save the data
with intervals much shorter than the time scale of interest. Another
implication of eqn. [eqn:corrmd] is that in principle we can not compute
all points of the ACF with the same accuracy, since we have :math:`N-1`
data points for :math:`C_f(\Delta t)` but only 1 for
:math:`C_f((N-1)\Delta t)`. However, if we decide to compute only an ACF
of length :math:`M\Delta t`, where :math:`M \leq N/2` we can compute all
points with the same statistical accuracy:

.. math:: C_f(j\Delta t)  ~=~ \frac{1}{M}\sum_{i=0}^{N-1-M} f(i\Delta t)f((i+j)\Delta t)

 Here of course :math:`j < M`. :math:`M` is sometimes referred to as the
time lag of the correlation function. When we decide to do this, we
intentionally do not use all the available points for very short time
intervals (:math:`j << M`), but it makes it easier to interpret the
results. Another aspect that may not be neglected when computing ACFs
from simulation is that usually the time origins :math:`\xi`
(eqn. [eqn:corr]) are not statistically independent, which may introduce
a bias in the results. This can be tested using a block-averaging
procedure, where only time origins with a spacing at least the length of
the time lag are included, *e.g.* using :math:`k` time origins with
spacing of :math:`M\Delta t` (where :math:`kM \leq N`):

.. math:: C_f(j\Delta t)  ~=~ \frac{1}{k}\sum_{i=0}^{k-1} f(iM\Delta t)f((iM+j)\Delta t)

 However, one needs very long simulations to get good accuracy this way,
because there are many fewer points that contribute to the ACF.

Using FFT for computation of the ACF
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The computational cost for calculating an ACF according to
eqn. [eqn:corrmd] is proportional to :math:`N^2`, which is considerable.
However, this can be improved by using fast Fourier transforms to do the
convolution \ `108 <#ref-Allen87>`__.

Special forms of the ACF
~~~~~~~~~~~~~~~~~~~~~~~~

There are some important varieties on the ACF, *e.g.* the ACF of a
vector :

.. math::

   C_{{\mbox{\boldmath ${p}$}}}(t) ~=~       \int_0^{\infty} P_n(\cos\angle\left({\mbox{\boldmath ${p}$}}(\xi),{\mbox{\boldmath ${p}$}}(\xi+t)\right) {\rm d} \xi
   \label{eqn:corrleg}

 where :math:`P_n(x)` is the :math:`n^{th}` order Legendre
polynomial. [3]_ Such correlation times can actually be obtained
experimentally using *e.g.* NMR or other relaxation experiments. GROMACS
can compute correlations using the 1\ :math:`^{st}` and 2\ :math:`^{nd}`
order Legendre polynomial (eqn. [eqn:corrleg]). This can also be used
for rotational autocorrelation (gmx rotacf) and dipole autocorrelation
(gmx dipoles).

In order to study torsion angle dynamics, we define a dihedral
autocorrelation function as \ `159 <#ref-Spoel97a>`__:

.. math::

   C(t)    ~=~     \left\langle \cos(\theta(\tau)-\theta(\tau+t))\right\rangle_{\tau}
   \label{eqn:coenk}

 **Note** that this is not a product of two functions as is generally
used for correlation functions, but it may be rewritten as the sum of
two products:

.. math::

   C(t)    ~=~     \left\langle\cos(\theta(\tau))\cos(\theta(\tau+t))\,+\,\sin(\theta(\tau))\sin(\theta(\tau+t))\right\rangle_{\tau}
   \label{eqn:cot}

Some Applications
~~~~~~~~~~~~~~~~~

The program gmx velacc calculates the *velocity autocorrelation
function*.

.. math:: C_{{\mbox{\boldmath ${v}$}}} (\tau) ~=~ \langle {{\mbox{\boldmath ${v}$}}}_i(\tau) \cdot {{\mbox{\boldmath ${v}$}}}_i(0) \rangle_{i \in A}

 The self diffusion coefficient can be calculated using the Green-Kubo
relation \ `108 <#ref-Allen87>`__:

.. math:: D_A ~=~ {1\over 3} \int_0^{\infty} \langle {\bf v}_i(t) \cdot {\bf v}_i(0) \rangle_{i \in A} \; dt

 which is just the integral of the velocity autocorrelation function.
There is a widely-held belief that the velocity ACF converges faster
than the mean square displacement (sec. [sec:msd]), which can also be
used for the computation of diffusion constants. However, Allen &
Tildesley \ `108 <#ref-Allen87>`__ warn us that the long-time
contribution to the velocity ACF can not be ignored, so care must be
taken.

Another important quantity is the dipole correlation time. The *dipole
correlation function* for particles of type :math:`A` is calculated as
follows by gmx dipoles:

.. math::

   C_{\mu} (\tau) ~=~
   \langle {\bf \mu}_i(\tau) \cdot {\bf \mu}_i(0) \rangle_{i \in A}

 with :math:`{\bf \mu}_i = \sum_{j \in i} {\bf r}_j q_j`. The dipole
correlation time can be computed using eqn. [eqn:corrtime]. For some
applications see \ **???**.

The viscosity of a liquid can be related to the correlation time of the
Pressure tensor
:math:`{\mbox{\boldmath ${P}$}}` `160 <#ref-PSmith93c>`__,
`161 <#ref-Balasubramanian96>`__. gmx energy can compute the viscosity,
but this is not very accurate \ `149 <#ref-Hess2002a>`__, and actually
the values do not converge.

Curve fitting in GROMACS
------------------------

Sum of exponential functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it is useful to fit a curve to an analytical function, for
example in the case of autocorrelation functions with noisy tails.
GROMACS is not a general purpose curve-fitting tool however and
therefore GROMACS only supports a limited number of functions.
Table [tab:fitfn] lists the available options with the corresponding
command-line options. The underlying routines for fitting use the
Levenberg-Marquardt algorithm as implemented in the lmfit
package \ `162 <#ref-lmfit>`__ (a bare-bones version of which is
included in GROMACS in which an option for error-weighted fitting was
implemented).

Error estimation
~~~~~~~~~~~~~~~~

Under the hood GROMACS implements some more fitting functions, namely a
function to estimate the error in time-correlated data due to
Hess \ `149 <#ref-Hess2002a>`__:

.. math::

   \varepsilon^2(t) =
   \alpha\tau_1\left(1+\frac{\tau_1}{t}\left(e^{-t/\tau_1}-1\right)\right)
         + (1-\alpha)\tau_2\left(1+\frac{\tau_2}{t}\left(e^{-t/\tau_2}-1\right)\right)

 where :math:`\tau_1` and :math:`\tau_2` are time constants (with
:math:`\tau_2 \ge \tau_1`) and :math:`\alpha` usually is close to 1 (in
the fitting procedure it is enforced that :math:`0\leq\alpha\leq 1`).
This is used in gmx analyze for error estimation using

.. math:: \lim_{t\rightarrow\infty}\varepsilon(t) = \sigma\sqrt{\frac{2(\alpha\tau_1+(1-\alpha)\tau_2)}{T}}

 where :math:`\sigma` is the standard deviation of the data set and
:math:`T` is the total simulation time \ `149 <#ref-Hess2002a>`__.

Interphase boundary demarcation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to determine the position and width of an interface,
Steen-Sæthre *et al.* fitted a density profile to the following function

.. math::

   f(x) ~=~ \frac{a_0+a_1}{2} - \frac{a_0-a_1}{2}{\rm
     erf}\left(\frac{x-a_2}{a_3^2}\right)

 where :math:`a_0` and :math:`a_1` are densities of different phases,
:math:`x` is the coordinate normal to the interface, :math:`a_2` is the
position of the interface and :math:`a_3` is the width of the
interface \ `163 <#ref-Steen-Saethre2014a>`__. This is implemented in
gmx densorder.

Transverse current autocorrelation function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In order to establish the transverse current autocorrelation function
(useful for computing viscosity \ `164 <#ref-Palmer1994a>`__) the
following function is fitted:

.. math::

   f(x) ~=~ e^{-\nu}\left({\rm cosh}(\omega\nu)+\frac{{\rm
       sinh}(\omega\nu)}{\omega}\right)

 with :math:`\nu = x/(2a_0)` and :math:`\omega = \sqrt{1-a_1}`. This is
implemented in gmx tcaf.

Viscosity estimation from pressure autocorrelation function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The viscosity is a notoriously difficult property to extract from
simulations \ `149 <#ref-Hess2002a>`__, `165 <#ref-Wensink2003a>`__. It
is *in principle* possible to determine it by integrating the pressure
autocorrelation function \ `160 <#ref-PSmith93c>`__, however this is
often hampered by the noisy tail of the ACF. A workaround to this is
fitting the ACF to the following function \ `166 <#ref-Guo2002b>`__:

.. math::

   f(t)/f(0) = (1-C) {\rm cos}(\omega t) e^{-(t/\tau_f)^{\beta_f}} + C
   e^{-(t/\tau_s)^{\beta_s}}

 where :math:`\omega` is the frequency of rapid pressure oscillations
(mainly due to bonded forces in molecular simulations), :math:`\tau_f`
and :math:`\beta_f` are the time constant and exponent of fast
relaxation in a stretched-exponential approximation, :math:`\tau_s` and
:math:`\beta_s` are constants for slow relaxation and :math:`C` is the
pre-factor that determines the weight between fast and slow relaxation.
After a fit, the integral of the function :math:`f(t)` is used to
compute the viscosity:

.. math:: \eta = \frac{V}{k_B T}\int_0^{\infty} f(t) dt

 This equation has been applied to computing the bulk and shear
viscosity using different elements from the pressure
tensor \ `167 <#ref-Fanourgakis2012a>`__. This is implemented in gmx
viscosity.

Mean Square Displacement
------------------------

| gmx msd
| To determine the self diffusion coefficient :math:`D_A` of particles
  of type :math:`A`, one can use the Einstein
  relation `108 <#ref-Allen87>`__:

  .. math::

     \lim_{t \rightarrow \infty} \langle
     \|{\bf r}_i(t) - {\bf r}_i(0)\|^2 \rangle_{i \in A} ~=~ 6 D_A t

   This *mean square displacement* and :math:`D_A` are calculated by the
  program gmx msd. Normally an index file containing atom numbers is
  used and the MSD is averaged over these atoms. For molecules
  consisting of more than one atom, :math:`{\bf r}_i` can be taken as
  the center of mass positions of the molecules. In that case, you
  should use an index file with molecule numbers. The results will be
  nearly identical to averaging over atoms, however. The gmx msd program
  can also be used for calculating diffusion in one or two dimensions.
  This is useful for studying lateral diffusion on interfaces.

An example of the mean square displacement of SPC water is given in
Fig. [fig:msdwater].

|Mean Square Displacement of SPC-water.|

Bonds/distances, angles and dihedrals
-------------------------------------

| gmx distance, gmx angle, gmx gangle
| To monitor specific *bonds* in your modules, or more generally
  distances between points, the program gmx distance can calculate
  distances as a function of time, as well as the distribution of the
  distance. With a traditional index file, the groups should consist of
  pairs of atom numbers, for example:

::

    [ bonds_1 ]
     1     2
     3     4
     9    10

    [ bonds_2 ]
    12    13

Selections are also supported, with first two positions defining the
first distance, second pair of positions defining the second distance
and so on. You can calculate the distances between CA and CB atoms in
all your residues (assuming that every residue either has both atoms, or
neither) using a selection such as:

::

    name CA CB

The selections also allow more generic distances to be computed. For
example, to compute the distances between centers of mass of two
residues, you can use:

::

    com of resname AAA plus com of resname BBB

The program gmx angle calculates the distribution of *angles* and
*dihedrals* in time. It also gives the average angle or dihedral. The
index file consists of triplets or quadruples of atom numbers:

::

    [ angles ]
     1     2     3
     2     3     4
     3     4     5

    [ dihedrals ]
     1     2     3     4
     2     3     5     5

For the dihedral angles you can use either the “biochemical convention”
(:math:`\phi = 0 \equiv cis`) or “polymer convention”
(:math:`\phi = 0 \equiv trans`), see Fig. [fig:dih\_def].

|Dihedral conventions: A. “Biochemical convention”. B. “Polymer
convention”.|

The program gmx gangle provides a selection-enabled version to compute
angles. This tool can also compute angles and dihedrals, but does not
support all the options of gmx angle, such as autocorrelation or other
time series analyses. In addition, it supports angles between two
vectors, a vector and a plane, two planes (defined by 2 or 3 points,
respectively), a vector/plane and the :math:`z` axis, or a vector/plane
and the normal of a sphere (determined by a single position). Also the
angle between a vector/plane compared to its position in the first frame
is supported. For planes, gmx gangle uses the normal vector
perpendicular to the plane. See Fig. [fig:sgangle]A, B, C) for the
definitions.

|Angle options of gmx gangle: A. Angle between two vectors. B. Angle
between two planes. C. Angle between a vector and the :math:`z` axis. D.
Angle between a vector and the normal of a sphere. Also other
combinations are supported: planes and vectors can be used
interchangeably.|

Radius of gyration and distances
--------------------------------

| gmx gyrate, gmx distance, gmx mindist, gmx mdmat, gmx pairdist, gmx
  xpm2ps
| To have a rough measure for the compactness of a structure, you can
  calculate the *radius of gyration* with the program gmx gyrate as
  follows:

  .. math::

     R_g ~=~ \left({\frac{\sum_i \|{\bf r}_i\|^2 m_i}{\sum_i m_i}}\right)^{{\frac{1}{2}}}
     \label{eqn:rg}

   where :math:`m_i` is the mass of atom :math:`i` and :math:`{\bf r}_i`
  the position of atom :math:`i` with respect to the center of mass of
  the molecule. It is especially useful to characterize polymer
  solutions and proteins. The program will also provide the radius of
  gyration around the coordinate axis (or, optionally, principal axes)
  by only summing the radii components orthogonal to each axis, for
  instance

  .. math::

     R_{g,x} ~=~ \left({\frac{\sum_i \left( r_{i,y}^2 + r_{i,z}^2 \right) m_i}{\sum_i m_i}}\right)^{{\frac{1}{2}}}
     \label{eqn:rgaxis}

Sometimes it is interesting to plot the *distance* between two atoms, or
the *minimum* distance between two groups of atoms (*e.g.*: protein
side-chains in a salt bridge). To calculate these distances between
certain groups there are several possibilities:

:math:`\bullet`
    The *distance between the geometrical centers* of two groups can be
    calculated with the program gmx distance, as explained in
    sec. [sec:bad].

:math:`\bullet`
    The *minimum distance* between two groups of atoms during time can
    be calculated with the program gmx mindist. It also calculates the
    *number of contacts* between these groups within a certain radius
    :math:`r_{max}`.

:math:`\bullet`
    gmx pairdist is a selection-enabled version of gmx mindist.

:math:`\bullet`
    To monitor the *minimum distances between amino acid residues*
    within a (protein) molecule, you can use the program gmx mdmat. This
    minimum distance between two residues A\ :math:`_i` and
    A\ :math:`_j` is defined as the smallest distance between any pair
    of atoms (i :math:`\in` A\ :math:`_i`, j :math:`\in` A\ :math:`_j`).
    The output is a symmetrical matrix of smallest distances between all
    residues. To visualize this matrix, you can use a program such as
    xv. If you want to view the axes and legend or if you want to print
    the matrix, you can convert it with xpm2ps into a Postscript
    picture, see Fig. [fig:mdmat].

    .. figure:: plots/distm
       :alt: A minimum distance matrix for a
       peptide \ `168 <#ref-Spoel96b>`__.
       :width: 6.50000cm

       A minimum distance matrix for a
       peptide \ `168 <#ref-Spoel96b>`__.

    Plotting these matrices for different time-frames, one can analyze
    changes in the structure, and *e.g.* forming of salt bridges.

Root mean square deviations in structure
----------------------------------------

| gmx rms, gmx rmsdist
| The *root mean square deviation* (:math:`RMSD`) of certain atoms in a
  molecule with respect to a reference structure can be calculated with
  the program gmx rms by least-square fitting the structure to the
  reference structure (:math:`t_2 = 0`) and subsequently calculating the
  :math:`RMSD` (eqn. [eqn:rmsd]).

  .. math::

     RMSD(t_1,t_2) ~=~ \left[\frac{1}{M} \sum_{i=1}^N m_i \|{\bf r}_i(t_1)-{\bf r}_i(t_2)\|^2 \right]^{\frac{1}{2}}
     \label{eqn:rmsd}

   where :math:`M = \sum_{i=1}^N m_i` and :math:`{\bf r}_i(t)` is the
  position of atom :math:`i` at time :math:`t`. **Note** that fitting
  does not have to use the same atoms as the calculation of the
  :math:`RMSD`; *e.g.* a protein is usually fitted on the backbone atoms
  (N,C:math:`_{\alpha}`,C), but the :math:`RMSD` can be computed of the
  backbone or of the whole protein.

Instead of comparing the structures to the initial structure at time
:math:`t=0` (so for example a crystal structure), one can also calculate
eqn. [eqn:rmsd] with a structure at time :math:`t_2=t_1-\tau`. This
gives some insight in the mobility as a function of :math:`\tau`. A
matrix can also be made with the :math:`RMSD` as a function of
:math:`t_1` and :math:`t_2`, which gives a nice graphical interpretation
of a trajectory. If there are transitions in a trajectory, they will
clearly show up in such a matrix.

Alternatively the :math:`RMSD` can be computed using a fit-free method
with the program gmx rmsdist:

.. math::

   RMSD(t) ~=~     \left[\frac{1}{N^2}\sum_{i=1}^N \sum_{j=1}^N    \|{\bf r}_{ij}(t)-{\bf r}_{ij}(0)\|^2\right]^{\frac{1}{2}}
   \label{eqn:rmsdff}

 where the *distance* **r**\ :math:`_{ij}` between atoms at time
:math:`t` is compared with the distance between the same atoms at time
:math:`0`.

Covariance analysis
-------------------

Covariance analysis, also called principal component analysis or
essential dynamics `169 <#ref-Amadei93>`__\ , can find correlated
motions. It uses the covariance matrix :math:`C` of the atomic
coordinates:

.. math::

   C_{ij} = \left \langle 
   M_{ii}^{\frac{1}{2}} (x_i - \langle x_i \rangle)
   M_{jj}^{\frac{1}{2}}  (x_j - \langle x_j \rangle)
   \right \rangle

 where :math:`M` is a diagonal matrix containing the masses of the atoms
(mass-weighted analysis) or the unit matrix (non-mass weighted
analysis). :math:`C` is a symmetric :math:`3N \times 3N` matrix, which
can be diagonalized with an orthonormal transformation matrix :math:`R`:

.. math::

   R^T C R = \mbox{diag}(\lambda_1,\lambda_2,\ldots,\lambda_{3N})
   ~~~~\mbox{where}~~\lambda_1 \geq \lambda_2 \geq \ldots \geq \lambda_{3N}

 The columns of :math:`R` are the eigenvectors, also called principal or
essential modes. :math:`R` defines a transformation to a new coordinate
system. The trajectory can be projected on the principal modes to give
the principal components :math:`p_i(t)`:

.. math:: {\bf p}(t) = R^T M^{\frac{1}{2}} ({\bf x}(t) - \langle {\bf x} \rangle)

 The eigenvalue :math:`\lambda_i` is the mean square fluctuation of
principal component :math:`i`. The first few principal modes often
describe collective, global motions in the system. The trajectory can be
filtered along one (or more) principal modes. For one principal mode
:math:`i` this goes as follows:

.. math::

   {\bf x}^f(t) =
   \langle {\bf x} \rangle + M^{-\frac{1}{2}} R_{*i} \, p_i(t)

When the analysis is performed on a macromolecule, one often wants to
remove the overall rotation and translation to look at the internal
motion only. This can be achieved by least square fitting to a reference
structure. Care has to be taken that the reference structure is
representative for the ensemble, since the choice of reference structure
influences the covariance matrix.

One should always check if the principal modes are well defined. If the
first principal component resembles a half cosine and the second
resembles a full cosine, you might be filtering noise (see below). A
good way to check the relevance of the first few principal modes is to
calculate the overlap of the sampling between the first and second half
of the simulation. **Note** that this can only be done when the same
reference structure is used for the two halves.

A good measure for the overlap has been defined
in \ `170 <#ref-Hess2002b>`__. The elements of the covariance matrix are
proportional to the square of the displacement, so we need to take the
square root of the matrix to examine the extent of sampling. The square
root can be calculated from the eigenvalues :math:`\lambda_i` and the
eigenvectors, which are the columns of the rotation matrix :math:`R`.
For a symmetric and diagonally-dominant matrix :math:`A` of size
:math:`3N \times 3N` the square root can be calculated as:

.. math::

   A^\frac{1}{2} = 
   R \, \mbox{diag}(\lambda_1^\frac{1}{2},\lambda_2^\frac{1}{2},\ldots,\lambda_{3N}^\frac{1}{2}) \, R^T

 It can be verified easily that the product of this matrix with itself
gives :math:`A`. Now we can define a difference :math:`d` between
covariance matrices :math:`A` and :math:`B` as follows:

.. math::

   \begin{aligned}
   d(A,B) & = & \sqrt{\mbox{tr}\left(\left(A^\frac{1}{2} - B^\frac{1}{2}\right)^2\right)
   }
   \\ & = &
   \sqrt{\mbox{tr}\left(A + B - 2 A^\frac{1}{2} B^\frac{1}{2}\right)}
   \\ & = &
   \left( \sum_{i=1}^N \left( \lambda_i^A + \lambda_i^B \right)
   - 2 \sum_{i=1}^N \sum_{j=1}^N \sqrt{\lambda_i^A \lambda_j^B}
   \left(R_i^A \cdot R_j^B\right)^2 \right)^\frac{1}{2}\end{aligned}

 where tr is the trace of a matrix. We can now define the overlap
:math:`s` as:

.. math:: s(A,B) = 1 - \frac{d(A,B)}{\sqrt{\mbox{tr}A + \mbox{tr} B}}

 The overlap is 1 if and only if matrices :math:`A` and :math:`B` are
identical. It is 0 when the sampled subspaces are completely orthogonal.

A commonly-used measure is the subspace overlap of the first few
eigenvectors of covariance matrices. The overlap of the subspace spanned
by :math:`m` orthonormal vectors :math:`{\bf w}_1,\ldots,{\bf w}_m` with
a reference subspace spanned by :math:`n` orthonormal vectors
:math:`{\bf v}_1,\ldots,{\bf v}_n` can be quantified as follows:

.. math::

   \mbox{overlap}({\bf v},{\bf w}) =
   \frac{1}{n} \sum_{i=1}^n \sum_{j=1}^m ({\bf v}_i \cdot {\bf w}_j)^2

 The overlap will increase with increasing :math:`m` and will be 1 when
set :math:`{\bf v}` is a subspace of set :math:`{\bf w}`. The
disadvantage of this method is that it does not take the eigenvalues
into account. All eigenvectors are weighted equally, and when degenerate
subspaces are present (equal eigenvalues), the calculated overlap will
be too low.

Another useful check is the cosine content. It has been proven that the
the principal components of random diffusion are cosines with the number
of periods equal to half the principal component
index \ `170 <#ref-Hess2002b>`__, `171 <#ref-Hess2000>`__. The
eigenvalues are proportional to the index to the power :math:`-2`. The
cosine content is defined as:

.. math::

   \frac{2}{T}
   \left( \int_0^T \cos\left(\frac{i \pi t}{T}\right) \, p_i(t) \mbox{d} t \right)^2
   \left( \int_0^T p_i^2(t) \mbox{d} t \right)^{-1}

 When the cosine content of the first few principal components is close
to 1, the largest fluctuations are not connected with the potential, but
with random diffusion.

The covariance matrix is built and diagonalized by gmx covar. The
principal components and overlap (and many more things) can be plotted
and analyzed with gmx anaeig. The cosine content can be calculated with
gmx analyze.

Dihedral principal component analysis
-------------------------------------

| gmx angle, gmx covar, gmx anaeig
| Principal component analysis can be performed in dihedral
  space \ `172 <#ref-Mu2005a>`__ using GROMACS. You start by defining
  the dihedral angles of interest in an index file, either using gmx
  mk\_angndx or otherwise. Then you use the gmx angle program with the
  -or flag to produce a new .trr file containing the cosine and sine of
  each dihedral angle in two coordinates, respectively. That is, in the
  .trr file you will have a series of numbers corresponding to:
  cos(\ :math:`\phi_1`), sin(\ :math:`\phi_1`), cos(\ :math:`\phi_2`),
  sin(\ :math:`\phi_2`), ..., cos(\ :math:`\phi_n`),
  sin(\ :math:`\phi_n`), and the array is padded with zeros, if
  necessary. Then you can use this .trr file as input for the gmx covar
  program and perform principal component analysis as usual. For this to
  work you will need to generate a reference file (.tpr, .gro, .pdb
  etc.) containing the same number of “atoms” as the new .trr file, that
  is for :math:`n` dihedrals you need 2\ :math:`n`/3 atoms (rounded up
  if not an integer number). You should use the -nofit option for gmx
  covar since the coordinates in the dummy reference file do not
  correspond in any way to the information in the .trr file. Analysis of
  the results is done using gmx anaeig.

Hydrogen bonds
--------------

| gmx hbond
| The program gmx hbond analyzes the *hydrogen bonds* (H-bonds) between
  all possible donors D and acceptors A. To determine if an H-bond
  exists, a geometrical criterion is used, see also Fig. [fig:hbond]:

  .. math::

     \begin{array}{rclcl}
     r       & \leq  & r_{HB}        & = & 0.35~\mbox{nm}    \\
     \alpha  & \leq  & \alpha_{HB}   & = & 30^o              \\
     \end{array}

.. figure:: plots/hbond
   :alt: Geometrical Hydrogen bond criterion.
   :width: 2.50000cm

   Geometrical Hydrogen bond criterion.

The value of :math:`r_{HB} = 0.35` nm corresponds to the first minimum
of the RDF of SPC water (see also Fig. [fig:rdf]).

The program gmx hbond analyzes all hydrogen bonds existing between two
groups of atoms (which must be either identical or non-overlapping) or
in specified donor-hydrogen-acceptor triplets, in the following ways:

|Insertion of water into an H-bond. (1) Normal H-bond between two
residues. (2) H-bonding bridge via a water molecule.|

-  Donor-Acceptor distance (:math:`r`) distribution of all H-bonds

-  Hydrogen-Donor-Acceptor angle (:math:`\alpha`) distribution of all
   H-bonds

-  The total number of H-bonds in each time frame

-  The number of H-bonds in time between residues, divided into groups
   :math:`n`-:math:`n`\ +\ :math:`i` where :math:`n` and
   :math:`n`\ +\ :math:`i` stand for residue numbers and :math:`i` goes
   from 0 to 6. The group for :math:`i=6` also includes all H-bonds for
   :math:`i>6`. These groups include the
   :math:`n`-:math:`n`\ +\ :math:`3`, :math:`n`-:math:`n`\ +\ :math:`4`
   and :math:`n`-:math:`n`\ +\ :math:`5` H-bonds, which provide a
   measure for the formation of :math:`\alpha`-helices or
   :math:`\beta`-turns or strands.

-  The lifetime of the H-bonds is calculated from the average over all
   autocorrelation functions of the existence functions (either 0 or 1)
   of all H-bonds:

   .. math::

      C(\tau) ~=~ \langle s_i(t)~s_i (t + \tau) \rangle
      \label{eqn:hbcorr}

    with :math:`s_i(t) = \{0,1\}` for H-bond :math:`i` at time
   :math:`t`. The integral of :math:`C(\tau)` gives a rough estimate of
   the average H-bond lifetime :math:`\tau_{HB}`:

   .. math::

      \tau_{HB} ~=~ \int_{0}^{\infty} C(\tau) d\tau
      \label{eqn:hblife}

    Both the integral and the complete autocorrelation function
   :math:`C(\tau)` will be output, so that more sophisticated analysis
   (*e.g.* using multi-exponential fits) can be used to get better
   estimates for :math:`\tau_{HB}`. A more complete analysis is given in
   ref. \ `173 <#ref-Spoel2006b>`__; one of the more fancy option is the
   Luzar and Chandler analysis of hydrogen bond
   kinetics \ `174 <#ref-Luzar96b>`__, `175 <#ref-Luzar2000a>`__.

-  An H-bond existence map can be generated of dimensions
   *# H-bonds*\ :math:`\times`\ *# frames*. The ordering is identical to
   the index file (see below), but reversed, meaning that the last
   triplet in the index file corresponds to the first row of the
   existence map.

-  Index groups are output containing the analyzed groups, all
   donor-hydrogen atom pairs and acceptor atoms in these groups,
   donor-hydrogen-acceptor triplets involved in hydrogen bonds between
   the analyzed groups and all solvent atoms involved in insertion.

Protein-related items
---------------------

| gmx do\_dssp, gmx rama, gmx wheel
| To analyze structural changes of a protein, you can calculate the
  radius of gyration or the minimum residue distances over time (see
  sec. [sec:rg]), or calculate the RMSD (sec. [sec:rmsd]).

You can also look at the changing of *secondary structure elements*
during your run. For this, you can use the program gmx do\_dssp, which
is an interface for the commercial program DSSP `176 <#ref-Kabsch83>`__.
For further information, see the DSSP manual. A typical output plot of
gmx do\_dssp is given in Fig. [fig:dssp].

.. figure:: plots/dssp
   :alt: Analysis of the secondary structure elements of a peptide in
   time.
   :width: 12.00000cm

   Analysis of the secondary structure elements of a peptide in time.

One other important analysis of proteins is the so-called *Ramachandran
plot*. This is the projection of the structure on the two dihedral
angles :math:`\phi` and :math:`\psi` of the protein backbone, see
Fig. [fig:phipsi].

.. figure:: plots/phipsi
   :alt: Definition of the dihedral angles :math:`\phi` and :math:`\psi`
   of the protein backbone.
   :width: 5.00000cm

   Definition of the dihedral angles :math:`\phi` and :math:`\psi` of
   the protein backbone.

To evaluate this Ramachandran plot you can use the program gmx rama. A
typical output is given in Fig. [fig:rama].

|Ramachandran plot of a small protein.|

When studying :math:`\alpha`-helices it is useful to have a *helical
wheel* projection of your peptide, to see whether a peptide is
amphipathic. This can be done using the gmx wheel program. Two examples
are plotted in Fig. [fig:wheel].

.. figure:: plots/hpr-wheel
   :alt: Helical wheel projection of the N-terminal helix of HPr.

   Helical wheel projection of the N-terminal helix of HPr.

Interface-related items
-----------------------

| gmx order, gmx density, gmx potential, gmx traj
| When simulating molecules with long carbon tails, it can be
  interesting to calculate their average orientation. There are several
  flavors of order parameters, most of which are related. The program
  gmx order can calculate order parameters using the equation:

.. math::

   S_{z} = \frac{3}{2}\langle {\cos^2{\theta_z}} \rangle - \frac{1}{2}
   \label{eqn:Sgr}

where :math:`\theta_z` is the angle between the :math:`z`-axis of the
simulation box and the molecular axis under consideration. The latter is
defined as the vector from C\ :math:`_{n-1}` to C\ :math:`_{n+1}`. The
parameters :math:`S_x` and :math:`S_y` are defined in the same way. The
brackets imply averaging over time and molecules. Order parameters can
vary between 1 (full order along the interface normal) and :math:`-1/2`
(full order perpendicular to the normal), with a value of zero in the
case of isotropic orientation.

The program can do two things for you. It can calculate the order
parameter for each CH\ :math:`_2` segment separately, for any of three
axes, or it can divide the box in slices and calculate the average value
of the order parameter per segment in one slice. The first method gives
an idea of the ordering of a molecule from head to tail, the second
method gives an idea of the ordering as function of the box length.

The electrostatic potential (:math:`\psi`) across the interface can be
computed from a trajectory by evaluating the double integral of the
charge density (:math:`\rho(z)`):

.. math::

   \psi(z) - \psi(-\infty) = - \int_{-\infty}^z dz' \int_{-\infty}^{z'} \rho(z'')dz''/ \epsilon_0 
   \label{eqn:elpotgr}

 where the position :math:`z=-\infty` is far enough in the bulk phase
such that the field is zero. With this method, it is possible to “split”
the total potential into separate contributions from lipid and water
molecules. The program gmx potential divides the box in slices and sums
all charges of the atoms in each slice. It then integrates this charge
density to give the electric field, which is in turn integrated to give
the potential. Charge density, electric field, and potential are written
to xvgr input files.

The program gmx traj is a very simple analysis program. All it does is
print the coordinates, velocities, or forces of selected atoms. It can
also calculate the center of mass of one or more molecules and print the
coordinates of the center of mass to three files. By itself, this is
probably not a very useful analysis, but having the coordinates of
selected molecules or atoms can be very handy for further analysis, not
only in interfacial systems.

The program gmx density calculates the mass density of groups and gives
a plot of the density against a box axis. This is useful for looking at
the distribution of groups or atoms across the interface.

Some implementation details
===========================

In this chapter we will present some implementation details. This is far
from complete, but we deemed it necessary to clarify some things that
would otherwise be hard to understand.

Single Sum Virial in GROMACS
----------------------------

The virial :math:`\Xi` can be written in full tensor form as:

.. math:: \Xi~=~-{\frac{1}{2}}~\sum_{i < j}^N~{{\mbox{\boldmath ${r}$}}_{ij}}\otimes{{\mbox{\boldmath ${F}$}}_{ij}}

 where :math:`\otimes` denotes the *direct product* of two vectors. [4]_
When this is computed in the inner loop of an MD program 9
multiplications and 9 additions are needed. [5]_

Here it is shown how it is possible to extract the virial calculation
from the inner loop \ `177 <#ref-Bekker93b>`__.

Virial
~~~~~~

In a system with periodic boundary conditions, the periodicity must be
taken into account for the virial:

.. math:: \Xi~=~-{\frac{1}{2}}~\sum_{i < j}^{N}~{{\mbox{\boldmath ${r}$}}_{ij}^n}\otimes{{\mbox{\boldmath ${F}$}}_{ij}}

 where :math:`{{\mbox{\boldmath ${r}$}}_{ij}^n}` denotes the distance
vector of the *nearest image* of atom :math:`i` from atom :math:`j`. In
this definition we add a *shift vector* :math:`\delta_i` to the position
vector :math:`{{\mbox{\boldmath ${r}$}}_i}` of atom :math:`i`. The
difference vector :math:`{{\mbox{\boldmath ${r}$}}_{ij}^n}` is thus
equal to:

.. math:: {{\mbox{\boldmath ${r}$}}_{ij}^n}~=~{{\mbox{\boldmath ${r}$}}_i}+\delta_i-{{\mbox{\boldmath ${r}$}}_j}

 or in shorthand:

.. math:: {{\mbox{\boldmath ${r}$}}_{ij}^n}~=~{{\mbox{\boldmath ${r}$}}_i^n}-{{\mbox{\boldmath ${r}$}}_j}

 In a triclinic system, there are 27 possible images of :math:`i`; when
a truncated octahedron is used, there are 15 possible images.

Virial from non-bonded forces
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here the derivation for the single sum virial in the *non-bonded force*
routine is given. There are a couple of considerations that are special
to GROMACS that we take into account:

-  When calculating short-range interactions, we apply the *minimum
   image convention* and only consider the closest image of each
   neighbor - and in particular we never allow interactions between a
   particle and any of its periodic images. For all the equations below,
   this means :math:`i \neq j`.

-  In general, either the :math:`i` or :math:`j` particle might be
   shifted to a neighbor cell to get the closest interaction (shift
   :math:`\delta_{ij}`). However, with minimum image convention there
   can be at most 27 different shifts for particles in the central cell,
   and for typical (very short-ranged) biomolecular interactions there
   are typically only a few different shifts involved for each particle,
   not to mention that each interaction can only be present for one
   shift.

-  For the GROMACS nonbonded interactions we use this to split the
   neighborlist of each :math:`i` particle into multiple separate lists,
   where each list has a constant shift :math:`\delta_i` for the
   :math:`i` partlcle. We can represent this as a sum over shifts (for
   which we use index :math:`s`), with the constraint that each particle
   interaction can only contribute to one of the terms in this sum, and
   the shift is no longer dependent on the :math:`j` particles. For any
   sum that does not contain complex dependence on :math:`s`, this means
   the sum trivially reduces to just the sum over :math:`i` and/or
   :math:`j`.

-  To simplify some of the sums, we replace sums over :math:`j<i` with
   double sums over all particles (remember, :math:`i \neq j`) and
   divide by 2.

Starting from the above definition of the virial, we then get

.. math::

   \begin{aligned}
   \Xi
   &~=~&-{\frac{1}{2}}~\sum_{i < j}^{N}~{\mathbf r}^n_{ij} \otimes {\mathbf F}_{ij} \nonumber \\
   &~=~&-{\frac{1}{2}}~\sum_{i < j}^{N}~\left( {\mathbf r}_i + \delta_{ij} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij} \nonumber \\
   &~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{j=1}^{N}~\left( {\mathbf r}_i + \delta_{ij} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij} \nonumber \\
   &~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} - {\mathbf r}_j \right) \otimes {\mathbf F}_{ij,s} \nonumber \\
   &~=~&-{\frac{1}{4}}~\sum_{i=}^{N}~\sum_{s}~\sum_{j=1}^{N}~\left( \left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} -{\mathbf r}_j \otimes {\mathbf F}_{ij,s} \right) \nonumber \\
   &~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^N ~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^{N} {\mathbf r}_j \otimes {\mathbf F}_{ij,s} \nonumber \\
   &~=~&-{\frac{1}{4}}~\sum_{i=1}^{N}~\sum_{s}~\sum_{j=1}^N ~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{i=1}^{N}~\sum_{j=1}^{N} {\mathbf r}_j \otimes {\mathbf F}_{ij} \nonumber \\
   &~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes ~\sum_{j=1}^N {\mathbf F}_{ij,s} + {\frac{1}{4}}\sum_{j=1}^N {\mathbf r}_j \otimes \sum_{i=1}^{N} {\mathbf F}_{ij} \nonumber \\
   &~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes ~\sum_{j=1}^N {\mathbf F}_{ij,s} - {\frac{1}{4}}\sum_{j=1}^N {\mathbf r}_j \otimes \sum_{i=1}^{N} {\mathbf F}_{ji} \nonumber \\
   &~=~&-{\frac{1}{4}}~\sum_{s}~\sum_{i=1}^{N}~\left( {\mathbf r}_i + \delta_{i,s} \right) \otimes {\mathbf F}_{i,s} - {\frac{1}{4}}\sum_{j=1}^N~{\mathbf r}_j \otimes {\mathbf F}_{j}  \nonumber \\
   &~=~&-{\frac{1}{4}}~\left(\sum_{i=1}^{N}~{\mathbf r}_i  \otimes {\mathbf F}_{i} + \sum_{j=1}^N~{\mathbf r}_j \otimes {\mathbf F}_{j} \right) - {\frac{1}{4}}\sum_{s}~\sum_{i=1}^{N} \delta_{i,s} \otimes {\mathbf F}_{i,s}  \nonumber \\
   &~=~&-{\frac{1}{2}}\sum_{i=1}^{N}~{\mathbf r}_i \otimes {\mathbf F}_{i} -{\frac{1}{4}}\sum_{s}~\sum_{i=1}^{N}~\delta_{i,s} \otimes {\mathbf F}_{i,s} \nonumber \\
   &~=~&-{\frac{1}{2}}\sum_{i=1}^{N}~{\mathbf r}_i \otimes {\mathbf F}_{i} -{\frac{1}{4}}\sum_{s}~\delta_{s} \otimes {\mathbf F}_{s} \nonumber \\
   &~=~&\Xi_0 + \Xi_1\end{aligned}

In the second-last stage, we have used the property that each shift
vector itself does not depend on the coordinates of particle :math:`i`,
so it is possible to sum up all forces corresponding to each shift
vector (in the nonbonded kernels), and then just use a sum over the
different shift vectors outside the kernels. We have also used

.. math::

   \begin{aligned}
   {{\mbox{\boldmath ${F}$}}_i}&~=~&\sum_{j=1}^N~{{\mbox{\boldmath ${F}$}}_{ij}}\\
   {{\mbox{\boldmath ${F}$}}_j}&~=~&\sum_{i=1}^N~{{\mbox{\boldmath ${F}$}}_{ji}}\end{aligned}

 which is the total force on :math:`i` with respect to :math:`j`.
Because we use Newton’s Third Law:

.. math:: {{\mbox{\boldmath ${F}$}}_{ij}}~=~-{{\mbox{\boldmath ${F}$}}_{ji}}

 we must, in the implementation, double the term containing the shift
:math:`\delta_i`. Similarly, in a few places we have summed the
shift-dependent force over all shifts to come up with the total force
per interaction or particle.

This separates the total virial :math:`\Xi` into a component
:math:`\Xi_0` that is a single sum over particles, and a second
component :math:`\Xi_1` that describes the influence of the particle
shifts, and that is only a sum over the different shift vectors.

The intra-molecular shift (mol-shift)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the bonded forces and SHAKE it is possible to make a *mol-shift*
list, in which the periodicity is stored. We simple have an array mshift
in which for each atom an index in the shiftvec array is stored.

The algorithm to generate such a list can be derived from graph theory,
considering each particle in a molecule as a begineqnarrayd in a graph,
the bonds as edges.

#. Represent the bonds and atoms as bidirectional graph

#. Make all atoms white

#. Make one of the white atoms black (atom :math:`i`) and put it in the
   central box

#. Make all of the neighbors of :math:`i` that are currently white, gray

#. Pick one of the gray atoms (atom :math:`j`), give it the correct
   periodicity with respect to any of its black neighbors and make it
   black

#. Make all of the neighbors of :math:`j` that are currently white, gray

#. If any gray atom remains, go to [5]

#. If any white atom remains, go to [3]

Using this algorithm we can

-  optimize the bonded force calculation as well as SHAKE

-  calculate the virial from the bonded forces in the single sum method
   again

Find a representation of the bonds as a bidirectional graph.

Virial from Covalent Bonds
~~~~~~~~~~~~~~~~~~~~~~~~~~

Since the covalent bond force gives a contribution to the virial, we
have:

.. math::

   \begin{aligned}
   b	&~=~&	\|{{\mbox{\boldmath ${r}$}}_{ij}^n}\|					\\
   V_b	&~=~&	{\frac{1}{2}}k_b(b-b_0)^2				\\
   {{\mbox{\boldmath ${F}$}}_i}&~=~&	-\nabla V_b					\\
   	&~=~&	k_b(b-b_0)\frac{{{\mbox{\boldmath ${r}$}}_{ij}^n}}{b}			\\
   {{\mbox{\boldmath ${F}$}}_j}&~=~&	-{{\mbox{\boldmath ${F}$}}_i}\end{aligned}

 The virial contribution from the bonds then is:

.. math::

   \begin{aligned}
   \Xi_b	&~=~&	-{\frac{1}{2}}({{\mbox{\boldmath ${r}$}}_i^n}\otimes{{\mbox{\boldmath ${F}$}}_i}~+~{{\mbox{\boldmath ${r}$}}_j}\otimes{{\mbox{\boldmath ${F}$}}_j})	\\
   	&~=~&	-{\frac{1}{2}}{{\mbox{\boldmath ${r}$}}_{ij}^n}\otimes{{\mbox{\boldmath ${F}$}}_i}\end{aligned}

Virial from SHAKE
~~~~~~~~~~~~~~~~~

An important contribution to the virial comes from shake. Satisfying the
constraints a force **G** that is exerted on the particles “shaken.” If
this force does not come out of the algorithm (as in standard SHAKE) it
can be calculated afterward (when using *leap-frog*) by:

.. math::

   \begin{aligned}
   \Delta{{\mbox{\boldmath ${r}$}}_i}&~=~&{{\mbox{\boldmath ${r}$}}_i}(t+{{\Delta t}})-
   [{{\mbox{\boldmath ${r}$}}_i}(t)+{\bf v}_i(t-\frac{{{\Delta t}}}{2}){{\Delta t}}+\frac{{{\mbox{\boldmath ${F}$}}_i}}{m_i}{{\Delta t}}^2]	\\
   {\bf G}_i&~=~&\frac{m_i\Delta{{\mbox{\boldmath ${r}$}}_i}}{{{\Delta t}}^2}\end{aligned}

 This does not help us in the general case. Only when no periodicity is
needed (like in rigid water) this can be used, otherwise we must add the
virial calculation in the inner loop of SHAKE.

When it *is* applicable the virial can be calculated in the single sum
way:

.. math:: \Xi~=~-{\frac{1}{2}}\sum_i^{N_c}~{{\mbox{\boldmath ${r}$}}_i}\otimes{{\mbox{\boldmath ${F}$}}_i}

 where :math:`N_c` is the number of constrained atoms.

Optimizations
-------------

Here we describe some of the algorithmic optimizations used in GROMACS,
apart from parallelism.

Inner Loops for Water
~~~~~~~~~~~~~~~~~~~~~

GROMACS uses special inner loops to calculate non-bonded interactions
for water molecules with other atoms, and yet another set of loops for
interactions between pairs of water molecules. There highly optimized
loops for two types of water models. For three site models similar to
SPC \ `80 <#ref-Berendsen81>`__, *i.e.*:

#. There are three atoms in the molecule.

#. The whole molecule is a single charge group.

#. The first atom has Lennard-Jones (sec. [sec:lj]) and Coulomb
   (sec. [sec:coul]) interactions.

#. Atoms two and three have only Coulomb interactions, and equal
   charges.

These loops also works for the SPC/E \ `178 <#ref-Berendsen87>`__ and
TIP3P \ `128 <#ref-Jorgensen83>`__ water models. And for four site water
models similar to TIP4P \ `128 <#ref-Jorgensen83>`__:

#. There are four atoms in the molecule.

#. The whole molecule is a single charge group.

#. The first atom has only Lennard-Jones (sec. [sec:lj]) interactions.

#. Atoms two and three have only Coulomb (sec. [sec:coul]) interactions,
   and equal charges.

#. Atom four has only Coulomb interactions.

The benefit of these implementations is that there are more
floating-point operations in a single loop, which implies that some
compilers can schedule the code better. However, it turns out that even
some of the most advanced compilers have problems with scheduling,
implying that manual tweaking is necessary to get optimum performance.
This may include common-sub-expression elimination, or moving code
around.

Averages and fluctuations
=========================

Formulae for averaging
----------------------

**Note:** this section was taken from ref \ `179 <#ref-Gunsteren94a>`__.

When analyzing a MD trajectory averages :math:`\left<x\right>` and
fluctuations

.. math::

   \left<(\Delta x)^2\right>^{{\frac{1}{2}}} ~=~ \left<[x-\left<x\right>]^2\right>^{{\frac{1}{2}}}
   \label{eqn:var0}

 of a quantity :math:`x` are to be computed. The variance
:math:`\sigma_x` of a series of N\ :math:`_x` values, {x:math:`_i`}, can
be computed from

.. math::

   \sigma_x~=~ \sum_{i=1}^{N_x} x_i^2 ~-~  \frac{1}{N_x}\left(\sum_{i=1}^{N_x}x_i\right)^2
   \label{eqn:var1}

 Unfortunately this formula is numerically not very accurate, especially
when :math:`\sigma_x^{{\frac{1}{2}}}` is small compared to the values of
:math:`x_i`. The following (equivalent) expression is numerically more
accurate

.. math:: \sigma_x ~=~ \sum_{i=1}^{N_x} [x_i  - \left<x\right>]^2

 with

.. math::

   \left<x\right> ~=~ \frac{1}{N_x} \sum_{i=1}^{N_x} x_i
   \label{eqn:var2}

 Using  eqns. [eqn:var1] and [eqn:var2] one has to go through the series
of :math:`x_i` values twice, once to determine :math:`\left<x\right>`
and again to compute :math:`\sigma_x`, whereas eqn. [eqn:var0] requires
only one sequential scan of the series {x:math:`_i`}. However, one may
cast eqn. [eqn:var1] in another form, containing partial sums, which
allows for a sequential update algorithm. Define the partial sum

.. math:: X_{n,m} ~=~ \sum_{i=n}^{m} x_i

 and the partial variance

.. math::

   \sigma_{n,m} ~=~ \sum_{i=n}^{m}  \left[x_i - \frac{X_{n,m}}{m-n+1}\right]^2  
   \label{eqn:sigma}

 It can be shown that

.. math::

   X_{n,m+k} ~=~  X_{n,m} + X_{m+1,m+k}         
   \label{eqn:Xpartial}

 and

.. math::

   \begin{aligned}
   \sigma_{n,m+k} &=& \sigma_{n,m} + \sigma_{m+1,m+k} + \left[~\frac {X_{n,m}}{m-n+1} - \frac{X_{n,m+k}}{m+k-n+1}~\right]^2~* \nonumber\\
      && ~\frac{(m-n+1)(m+k-n+1)}{k}
   \label{eqn:varpartial}\end{aligned}

 For :math:`n=1` one finds

.. math::

   \sigma_{1,m+k} ~=~ \sigma_{1,m} + \sigma_{m+1,m+k}~+~
     \left[~\frac{X_{1,m}}{m} - \frac{X_{1,m+k}}{m+k}~\right]^2~ \frac{m(m+k)}{k}
   \label{eqn:sig1}

 and for :math:`n=1` and :math:`k=1`  (eqn. [eqn:varpartial]) becomes

.. math::

   \begin{aligned}
   \sigma_{1,m+1}  &=& \sigma_{1,m} + 
                           \left[\frac{X_{1,m}}{m} - \frac{X_{1,m+1}}{m+1}\right]^2 m(m+1)\\
                   &=& \sigma_{1,m} + 
                           \frac {[~X_{1,m} - m x_{m+1}~]^2}{m(m+1)}
   \label{eqn:simplevar0}\end{aligned}

 where we have used the relation

.. math::

   X_{1,m+1} ~=~  X_{1,m} + x_{m+1}                       
   \label{eqn:simplevar1}

 Using formulae (eqn. [eqn:simplevar0]) and  (eqn. [eqn:simplevar1]) the
average

.. math:: \left<x\right> ~=~ \frac{X_{1,N_x}}{N_x}

 and the fluctuation

.. math:: \left<(\Delta x)^2\right>^{{\frac{1}{2}}} = \left[\frac {\sigma_{1,N_x}}{N_x}\right]^{{\frac{1}{2}}}

 can be obtained by one sweep through the data.

Implementation
--------------

In GROMACS the instantaneous energies :math:`E(m)` are stored in the
energy file, along with the values of :math:`\sigma_{1,m}` and
:math:`X_{1,m}`. Although the steps are counted from 0, for the energy
and fluctuations steps are counted from 1. This means that the equations
presented here are the ones that are implemented. We give somewhat
lengthy derivations in this section to simplify checking of code and
equations later on.

Part of a Simulation
~~~~~~~~~~~~~~~~~~~~

It is not uncommon to perform a simulation where the first part, *e.g.*
100 ps, is taken as equilibration. However, the averages and
fluctuations as printed in the log file are computed over the whole
simulation. The equilibration time, which is now part of the simulation,
may in such a case invalidate the averages and fluctuations, because
these numbers are now dominated by the initial drift towards
equilibrium.

Using eqns. [eqn:Xpartial] and [eqn:varpartial] the average and standard
deviation over part of the trajectory can be computed as:

.. math::

   \begin{aligned}
   X_{m+1,m+k}     &=& X_{1,m+k} - X_{1,m}                 \\
   \sigma_{m+1,m+k} &=& \sigma_{1,m+k}-\sigma_{1,m} - \left[~\frac{X_{1,m}}{m} - \frac{X_{1,m+k}}{m+k}~\right]^{2}~ \frac{m(m+k)}{k}\end{aligned}

or, more generally (with :math:`p \geq 1` and :math:`q \geq p`):

.. math::

   \begin{aligned}
   X_{p,q}         &=&     X_{1,q} - X_{1,p-1}     \\
   \sigma_{p,q}    &=&     \sigma_{1,q}-\sigma_{1,p-1} - \left[~\frac{X_{1,p-1}}{p-1} - \frac{X_{1,q}}{q}~\right]^{2}~ \frac{(p-1)q}{q-p+1}\end{aligned}

 **Note** that implementation of this is not entirely trivial, since
energies are not stored every time step of the simulation. We therefore
have to construct :math:`X_{1,p-1}` and :math:`\sigma_{1,p-1}` from the
information at time :math:`p` using eqns. [eqn:simplevar0] and
[eqn:simplevar1]:

.. math::

   \begin{aligned}
   X_{1,p-1}       &=&     X_{1,p} - x_p   \\
   \sigma_{1,p-1}  &=&     \sigma_{1,p} -  \frac {[~X_{1,p-1} - (p-1) x_{p}~]^2}{(p-1)p}\end{aligned}

Combining two simulations
~~~~~~~~~~~~~~~~~~~~~~~~~

Another frequently occurring problem is, that the fluctuations of two
simulations must be combined. Consider the following example: we have
two simulations (A) of :math:`n` and (B) of :math:`m` steps, in which
the second simulation is a continuation of the first. However, the
second simulation starts numbering from 1 instead of from :math:`n+1`.
For the partial sum this is no problem, we have to add :math:`X_{1,n}^A`
from run A:

.. math::

   X_{1,n+m}^{AB} ~=~ X_{1,n}^A + X_{1,m}^B
   \label{eqn:pscomb}

 When we want to compute the partial variance from the two components we
have to make a correction :math:`\Delta\sigma`:

.. math:: \sigma_{1,n+m}^{AB} ~=~ \sigma_{1,n}^A + \sigma_{1,m}^B +\Delta\sigma

 if we define :math:`x_i^{AB}` as the combined and renumbered set of
data points we can write:

.. math:: \sigma_{1,n+m}^{AB} ~=~ \sum_{i=1}^{n+m}  \left[x_i^{AB} - \frac{X_{1,n+m}^{AB}}{n+m}\right]^2

 and thus

.. math::

   \sum_{i=1}^{n+m}  \left[x_i^{AB} - \frac{X_{1,n+m}^{AB}}{n+m}\right]^2  ~=~
   \sum_{i=1}^{n}  \left[x_i^{A} - \frac{X_{1,n}^{A}}{n}\right]^2  +
   \sum_{i=1}^{m}  \left[x_i^{B} - \frac{X_{1,m}^{B}}{m}\right]^2  +\Delta\sigma

 or

.. math::

   \begin{aligned}
   \sum_{i=1}^{n+m}  \left[(x_i^{AB})^2 - 2 x_i^{AB}\frac{X^{AB}_{1,n+m}}{n+m} + \left(\frac{X^{AB}_{1,n+m}}{n+m}\right)^2  \right] &-& \nonumber \\
   \sum_{i=1}^{n}  \left[(x_i^{A})^2 - 2 x_i^{A}\frac{X^A_{1,n}}{n} + \left(\frac{X^A_{1,n}}{n}\right)^2  \right] &-& \nonumber \\
   \sum_{i=1}^{m}  \left[(x_i^{B})^2 - 2 x_i^{B}\frac{X^B_{1,m}}{m} + \left(\frac{X^B_{1,m}}{m}\right)^2  \right] &=& \Delta\sigma\end{aligned}

 all the :math:`x_i^2` terms drop out, and the terms independent of the
summation counter :math:`i` can be simplified:

.. math::

   \begin{aligned}
   \frac{\left(X^{AB}_{1,n+m}\right)^2}{n+m} \,-\, 
   \frac{\left(X^A_{1,n}\right)^2}{n} \,-\, 
   \frac{\left(X^B_{1,m}\right)^2}{m} &-& \nonumber \\
   2\,\frac{X^{AB}_{1,n+m}}{n+m}\sum_{i=1}^{n+m}x_i^{AB} \,+\,
   2\,\frac{X^{A}_{1,n}}{n}\sum_{i=1}^{n}x_i^{A} \,+\,
   2\,\frac{X^{B}_{1,m}}{m}\sum_{i=1}^{m}x_i^{B} &=& \Delta\sigma\end{aligned}

 we recognize the three partial sums on the second line and use
eqn. [eqn:pscomb] to obtain:

.. math:: \Delta\sigma ~=~ \frac{\left(mX^A_{1,n} - nX^B_{1,m}\right)^2}{nm(n+m)}

 if we check this by inserting :math:`m=1` we get back
eqn. [eqn:simplevar0]

Summing energy terms
~~~~~~~~~~~~~~~~~~~~

The gmx energy program can also sum energy terms into one, *e.g.*
potential + kinetic = total. For the partial averages this is again easy
if we have :math:`S` energy components :math:`s`:

.. math::

   X_{m,n}^S ~=~ \sum_{i=m}^n \sum_{s=1}^S x_i^s ~=~ \sum_{s=1}^S \sum_{i=m}^n x_i^s ~=~ \sum_{s=1}^S X_{m,n}^s
   \label{eqn:sumterms}

 For the fluctuations it is less trivial again, considering for example
that the fluctuation in potential and kinetic energy should cancel.
Nevertheless we can try the same approach as before by writing:

.. math:: \sigma_{m,n}^S ~=~ \sum_{s=1}^S \sigma_{m,n}^s + \Delta\sigma

 if we fill in eqn. [eqn:sigma]:

.. math::

   \sum_{i=m}^n \left[\left(\sum_{s=1}^S x_i^s\right) - \frac{X_{m,n}^S}{m-n+1}\right]^2 ~=~
   \sum_{s=1}^S \sum_{i=m}^n \left[\left(x_i^s\right) - \frac{X_{m,n}^s}{m-n+1}\right]^2 + \Delta\sigma
   \label{eqn:sigmaterms}

 which we can expand to:

.. math::

   \begin{aligned}
   &~&\sum_{i=m}^n \left[\sum_{s=1}^S (x_i^s)^2 + \left(\frac{X_{m,n}^S}{m-n+1}\right)^2 -2\left(\frac{X_{m,n}^S}{m-n+1}\sum_{s=1}^S x_i^s + \sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'} \right)\right]    \nonumber \\
   &-&\sum_{s=1}^S \sum_{i=m}^n \left[(x_i^s)^2 - 2\,\frac{X_{m,n}^s}{m-n+1}\,x_i^s + \left(\frac{X_{m,n}^s}{m-n+1}\right)^2\right] ~=~\Delta\sigma \end{aligned}

 the terms with :math:`(x_i^s)^2` cancel, so that we can simplify to:

.. math::

   \begin{aligned}
   &~&\frac{\left(X_{m,n}^S\right)^2}{m-n+1} -2 \frac{X_{m,n}^S}{m-n+1}\sum_{i=m}^n\sum_{s=1}^S x_i^s -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, -        \nonumber \\
   &~&\sum_{s=1}^S \sum_{i=m}^n \left[- 2\,\frac{X_{m,n}^s}{m-n+1}\,x_i^s + \left(\frac{X_{m,n}^s}{m-n+1}\right)^2\right] ~=~\Delta\sigma \end{aligned}

 or

.. math:: -\frac{\left(X_{m,n}^S\right)^2}{m-n+1}  -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, +  \sum_{s=1}^S \frac{\left(X_{m,n}^s\right)^2}{m-n+1}  ~=~\Delta\sigma

 If we now expand the first term using eqn. [eqn:sumterms] we obtain:

.. math:: -\frac{\left(\sum_{s=1}^SX_{m,n}^s\right)^2}{m-n+1}  -2\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\, +      \sum_{s=1}^S \frac{\left(X_{m,n}^s\right)^2}{m-n+1}  ~=~\Delta\sigma

 which we can reformulate to:

.. math:: -2\left[\sum_{s=1}^S \sum_{s'=s+1}^S X_{m,n}^s X_{m,n}^{s'}\,+\sum_{i=m}^n\sum_{s=1}^S \sum_{s'=s+1}^S x_i^s x_i^{s'}\right] ~=~\Delta\sigma

 or

.. math:: -2\left[\sum_{s=1}^S X_{m,n}^s \sum_{s'=s+1}^S X_{m,n}^{s'}\,+\,\sum_{s=1}^S \sum_{i=m}^nx_i^s \sum_{s'=s+1}^S x_i^{s'}\right] ~=~\Delta\sigma

 which gives

.. math:: -2\sum_{s=1}^S \left[X_{m,n}^s \sum_{s'=s+1}^S \sum_{i=m}^n x_i^{s'}\,+\,\sum_{i=m}^n x_i^s \sum_{s'=s+1}^S x_i^{s'}\right] ~=~\Delta\sigma

 Since we need all data points :math:`i` to evaluate this, in general
this is not possible. We can then make an estimate of
:math:`\sigma_{m,n}^S` using only the data points that are available
using the left hand side of eqn. [eqn:sigmaterms]. While the average can
be computed using all time steps in the simulation, the accuracy of the
fluctuations is thus limited by the frequency with which energies are
saved. Since this can be easily done with a program such as xmgr this is
not built-in in GROMACS.

.. raw:: html

   <div id="refs" class="references">

.. raw:: html

   <div id="ref-Bekker93a">

:sup:`1` H. Bekker, H.J.C. Berendsen, E.J. Dijkstra, S. Achterop, R. van
Drunen, D. van der Spoel, A. Sij:raw-latex:`\-`bers, and H. Keegstra *et
al.*, “Gromacs: A parallel computer for molecular dynamics simulations”;
pp. 252–256 in *Physics computing 92*. Edited by R.A. de Groot and J.
Nadrchal. World Scientific, Singapore, 1993.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Berendsen95a">

:sup:`2` H.J.C. Berendsen, D. van der Spoel, and R. van Drunen,
“GROMACS: A message-passing parallel molecular dynamics implementation,”
*Comp. Phys. Comm.*, **91** 43–56 (1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lindahl2001a">

:sup:`3` E. Lindahl, B. Hess, and D. van der Spoel, “GROMACS 3.0: A
package for molecular simulation and trajectory analysis,” *J. Mol.
Mod.*, **7** 306–317 (2001).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Spoel2005a">

:sup:`4` D. van der Spoel, E. Lindahl, B. Hess, G. Groenhof, A.E. Mark,
and H.J.C. Berendsen, “GROMACS: Fast, Flexible and Free,” *J. Comp.
Chem.*, **26** 1701–1718 (2005).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hess2008b">

:sup:`5` B. Hess, C. Kutzner, D. van der Spoel, and E. Lindahl, “GROMACS
4: Algorithms for Highly Efficient, Load-Balanced, and Scalable
Molecular Simulation,” *J. Chem. Theory Comput.*, **4** [3] 435–447
(2008).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Pronk2013">

:sup:`6` S. Pronk, S. Páll, R. Schulz, P. Larsson, P. Bjelkmar, R.
Apostolov, M.R. Shirts, and J.C. Smith *et al.*, “GROMACS 4.5: A
high-throughput and highly parallel open source molecular simulation
toolkit,” *Bioinformatics*, **29** [7] 845–854 (2013).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Pall2015">

:sup:`7` S. Páll, M.J. Abraham, C. Kutzner, B. Hess, and E. Lindahl,
“Tackling exascale software challenges in molecular dynamics simulations
with GROMACS”; pp. 3–27 in *Solving software challenges for exascale*.
Edited by S. Markidis and E. Laure. Springer International Publishing
Switzerland, London, 2015.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Abraham2015">

:sup:`8` M.J. Abraham, T. Murtola, R. Schulz, S. Páll, J.C. Smith, B.
Hess, and E. Lindahl, “GROMACS: High performance molecular simulations
through multi-level parallelism from laptops to supercomputers,”
*SoftwareX*, **1–2** 19–25 (2015).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Gunsteren90">

:sup:`9` W.F. van Gunsteren and H.J.C. Berendsen, “Computer simulation
of molecular dynamics: Methodology, applications, and perspectives in
chemistry,” *Angew. Chem. Int. Ed. Engl.*, **29** 992–1023 (1990).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Fraaije93">

:sup:`10` J.G.E.M. Fraaije, “Dynamic density functional theory for
microphase separation kinetics of block copolymer melts,” *J. Chem.
Phys.*, **99** 9202–9212 (1993).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-McQuarrie76">

:sup:`11` D.A. McQuarrie, *Statistical mechanics*. Harper & Row, New
York, 1976.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Gunsteren77">

:sup:`12` W.F. van Gunsteren and H.J.C. Berendsen, “Algorithms for
macromolecular dynamics and constraint dynamics,” *Mol. Phys.*, **34**
1311–1327 (1977).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Gunsteren82">

:sup:`13` W.F. van Gunsteren and M. Karplus, “Effect of constraints on
the dynamics of macromolecules,” *Macromolecules*, **15** 1528–1544
(1982).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Darden93">

:sup:`14` T. Darden, D. York, and L. Pedersen, “Particle mesh Ewald: An
N\ :math:`\bullet`\ log(N) method for Ewald sums in large systems,” *J.
Chem. Phys.*, **98** 10089–10092 (1993).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Essmann95">

:sup:`15` U. Essmann, L. Perera, M.L. Berkowitz, T. Darden, H. Lee, and
L.G. Pedersen, “A smooth particle mesh ewald potential,” *J. Chem.
Phys.*, **103** 8577–8592 (1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Geman84">

:sup:`16` S. Geman and D. Geman, “Stochastic relaxation, Gibbs
distributions and the Bayesian restoration of images,” *IEEE Trans.
Patt. Anal. Mach. Int.*, **6** 721 (1984).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Nilges88">

:sup:`17` M. Nilges, G.M. Clore, and A.M. Gronenborn, “Determination of
three-dimensional structures of proteins from interproton distance data
by dynamical simulated annealing from a random array of atoms,” *FEBS
Lett.*, **239** 129–136 (1988).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Schaik93">

:sup:`18` R.C. van Schaik, H.J.C. Berendsen, A.E. Torda, and W.F. van
Gunsteren, “A structure refinement method based on molecular dynamics in
4 spatial dimensions,” *J. Mol. Biol.*, **234** 751–762 (1993).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Zimmerman91">

:sup:`19` K. Zimmerman, “All purpose molecular mechanics simulator and
energy minimizer,” *J. Comp. Chem.*, **12** 310–319 (1991).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Adams79">

:sup:`20` D.J. Adams, E.M. Adams, and G.J. Hills, “The computer
simulation of polar liquids,” *Mol. Phys.*, **38** 387–400 (1979).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Bekker95">

:sup:`21` H. Bekker, E.J. Dijkstra, M.K.R. Renardus, and H.J.C.
Berendsen, “An efficient, box shape independent non-bonded force and
virial algorithm for molecular dynamics,” *Mol. Sim.*, **14** 137–152
(1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hockney74">

:sup:`22` R.W. Hockney, S.P. Goel, and J. Eastwood, “Quiet High
Resolution Computer Models of a Plasma,” *J. Comp. Phys.*, **14**
148–158 (1974).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Verlet67">

:sup:`23` L. Verlet., “Computer experiments on classical fluids. I.
Thermodynamical properties of Lennard-Jones molecules,” *Phys. Rev.*,
**159** 98–103 (1967).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Berendsen86b">

:sup:`24` H.J.C. Berendsen and W.F. van Gunsteren, “Practical algorithms
for dynamics simulations”; in 1986.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Swope82">

:sup:`25` W.C. Swope, H.C. Andersen, P.H. Berens, and K.R. Wilson, “A
computer-simulation method for the calculation of equilibrium-constants
for the formation of physical clusters of molecules: Application to
small water clusters,” *J. Chem. Phys.*, **76** 637–649 (1982).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Berendsen84">

:sup:`26` H.J.C. Berendsen, J.P.M. Postma, A. DiNola, and J.R. Haak,
“Molecular dynamics with coupling to an external bath,” *J. Chem.
Phys.*, **81** 3684–3690 (1984).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Andersen80">

:sup:`27` H.C. Andersen, “Molecular dynamics simulations at constant
pressure and/or temperature,” *J. Chem. Phys.*, **72** 2384 (1980).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Nose84">

:sup:`28` S. Nosé, “A molecular dynamics method for simulations in the
canonical ensemble,” *Mol. Phys.*, **52** 255–268 (1984).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hoover85">

:sup:`29` W.G. Hoover, “Canonical dynamics: Equilibrium phase-space
distributions,” *Phys. Rev. **A***, **31** 1695–1697 (1985).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Bussi2007a">

:sup:`30` G. Bussi, D. Donadio, and M. Parrinello, “Canonical sampling
through velocity rescaling,” *J. Chem. Phys.*, **126** 014101 (2007).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Berendsen91">

:sup:`31` H.J.C. Berendsen, “Transport properties computed by linear
response through weak coupling to a bath”; pp. 139–155 in *Computer
simulations in material science*. Edited by M. Meyer and V. Pontikis.
Kluwer, 1991.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Basconi2013">

:sup:`32` J.E. Basconi and M.R. Shirts, “Effects of temperature control
algorithms on transport properties and kinetics in molecular dynamics
simulations,” *J. Chem. Theory Comput.*, **9** [7] 2887–2899 (2013).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Cooke2008">

:sup:`33` B. Cooke and S.J. Schmidler, “Preserving the Boltzmann
ensemble in replica-exchange molecular dynamics,” *J. Chem. Phys.*,
**129** 164112 (2008).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Martyna1992">

:sup:`34` G.J. Martyna, M.L. Klein, and M.E. Tuckerman, “Nosé-Hoover
chains: The canonical ensemble via continuous dynamics,” *J. Chem.
Phys.*, **97** 2635–2643 (1992).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Martyna1996">

:sup:`35` G.J. Martyna, M.E. Tuckerman, D.J. Tobias, and M.L. Klein,
“Explicit reversible integrators for extended systems dynamics,” *Mol.
Phys.*, **87** 1117–1157 (1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Holian95">

:sup:`36` B.L. Holian, A.F. Voter, and R. Ravelo, “Thermostatted
molecular dynamics: How to avoid the Toda demon hidden in Nosé-Hoover
dynamics,” *Phys. Rev. E*, **52** [3] 2338–2347 (1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Eastwood2010">

:sup:`37` M.P. Eastwood, K.A. Stafford, R.A. Lippert, M.Ø. Jensen, P.
Maragakis, C. Predescu, R.O. Dror, and D.E. Shaw, “Equipartition and the
calculation of temperature in biomolecular simulations,” *J. Chem.
Theory Comput.*, **ASAP** DOI: 10.1021/ct9002916 (2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Parrinello81">

:sup:`38` M. Parrinello and A. Rahman, “Polymorphic transitions in
single crystals: A new molecular dynamics method,” *J. Appl. Phys.*,
**52** 7182–7190 (1981).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Nose83">

:sup:`39` S. Nosé and M.L. Klein, “Constant pressure molecular dynamics
for molecular systems,” *Mol. Phys.*, **50** 1055–1076 (1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Liu2015">

:sup:`40` G. Liu, “Dynamical equations for the period vectors in a
periodic system under constant external stress,” *Can. J. Phys.*, **93**
974–978 (2015).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Tuckerman2006">

:sup:`41` M.E. Tuckerman, J. Alejandre, R. López-Rendón, A.L. Jochim,
and G.J. Martyna, “A Liouville-operator derived measure-preserving
integrator for molecular dynamics simulations in the isothermal-isobaric
ensemble,” *J. Phys. A.*, **59** 5629–5651 (2006).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Yu2010">

:sup:`42` T.-Q. Yu, J. Alejandre, R. Lopez-Rendon, G.J. Martyna, and
M.E. Tuckerman, “Measure-preserving integrators for molecular dynamics
in the isothermal-isobaric ensemble derived from the liouville
operator,” *Chem. Phys.*, **370** 294–305 (2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Dick58">

:sup:`43` B.G. Dick and A.W. Overhauser, “Theory of the dielectric
constants of alkali halide crystals,” *Phys. Rev.*, **112** 90–103
(1958).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Jordan95">

:sup:`44` P.C. Jordan, P.J. van Maaren, J. Mavri, D. van der Spoel, and
H.J.C. Berendsen, “Towards phase transferable potential functions:
Methodology and application to nitrogen,” *J. Chem. Phys.*, **103**
2272–2285 (1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Maaren2001a">

:sup:`45` P.J. van Maaren and D. van der Spoel, “Molecular dynamics
simulations of a water with a novel shell-model potential,” *J. Phys.
Chem. B.*, **105** 2618–2626 (2001).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Ryckaert77">

:sup:`46` J.P. Ryckaert, G. Ciccotti, and H.J.C. Berendsen, “Numerical
integration of the cartesian equations of motion of a system with
constraints; molecular dynamics of n-alkanes,” *J. Comp. Phys.*, **23**
327–341 (1977).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Miyamoto92">

:sup:`47` S. Miyamoto and P.A. Kollman, “SETTLE: An analytical version
of the SHAKE and RATTLE algorithms for rigid water models,” *J. Comp.
Chem.*, **13** 952–962 (1992).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Andersen1983a">

:sup:`48` H.C. Andersen, “RATTLE: A ‘Velocity’ version of the SHAKE
algorithm for molecular dynamics calculations,” *J. Comp. Phys.*, **52**
24–34 (1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hess97">

:sup:`49` B. Hess, H. Bekker, H.J.C. Berendsen, and J.G.E.M. Fraaije,
“LINCS: A linear constraint solver for molecular simulations,” *J. Comp.
Chem.*, **18** 1463–1472 (1997).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hess2008a">

:sup:`50` B. Hess, “P-LINCS: A parallel linear constraint solver for
molecular simulation,” *J. Chem. Theory Comput.*, **4** 116–122 (2007).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Goga2012">

:sup:`51` N. Goga, A.J. Rzepiela, A.H. de Vries, S.J. Marrink, and
H.J.C. Berendsen, “Efficient algorithms for Langevin and DPD dynamics,”
*J. Chem. Theory Comput.*, **8** 3637–3649 (2012).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Byrd95a">

:sup:`52` R.H. Byrd, P. Lu, and J. Nocedal, “A limited memory algorithm
for bound constrained optimization,” *SIAM J. Scientif. Statistic.
Comput.*, **16** 1190–1208 (1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Zhu97a">

:sup:`53` C. Zhu, R.H. Byrd, and J. Nocedal, “L-BFGS-B: Algorithm 778:
L-BFGS-B, FORTRAN routines for large scale bound constrained
optimization,” *ACM Trans. Math. Softw.*, **23** 550–560 (1997).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Levitt83">

:sup:`54` M. Levitt, C. Sander, and P.S. Stern, “The normal modes of a
protein: Native bovine pancreatic trypsin inhibitor,” *Int. J. Quant.
Chem: Quant. Biol. Symp.*, **10** 181–199 (1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Go83">

:sup:`55` N. G\ :math:`\bar{\rm o}`, T. Noguti, and T. Nishikawa,
“Dynamics of a small globular protein in terms of low-frequency
vibrational modes,” *Proc. Natl. Acad. Sci. USA*, **80** 3696–3700
(1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-BBrooks83b">

:sup:`56` B. Brooks and M. Karplus, “Harmonic dynamics of proteins:
Normal modes and fluctuations in bovine pancreatic trypsin inhibitor,”
*Proc. Natl. Acad. Sci. USA*, **80** 6571–6575 (1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hayward95b">

:sup:`57` S. Hayward and N. G\ :math:`\bar{\rm o}`, “Collective variable
description of native protein dynamics,” *Annu. Rev. Phys. Chem.*,
**46** 223–250 (1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Bennett1976">

:sup:`58` C.H. Bennett, “Efficient Estimation of Free Energy Differences
from Monte Carlo Data,” *J. Comp. Phys.*, **22** 245–268 (1976).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Shirts2008">

:sup:`59` M.R. Shirts and J.D. Chodera, “Statistically optimal analysis
of multiple equilibrium simulations,” *J. Chem. Phys.*, **129** 124105
(2008).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hukushima96a">

:sup:`60` K. Hukushima and K. Nemoto, “Exchange Monte Carlo Method and
Application to Spin Glass Simulations,” *J. Phys. Soc. Jpn.*, **65**
1604–1608 (1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Sugita99">

:sup:`61` Y. Sugita and Y. Okamoto, “Replica-exchange molecular dynamics
method for protein folding,” *Chem. Phys. Lett.*, **314** 141–151
(1999).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Seibert2005a">

:sup:`62` M. Seibert, A. Patriksson, B. Hess, and D. van der Spoel,
“Reproducible polypeptide folding and structure prediction using
molecular dynamics simulations,” *J. Mol. Biol.*, **354** 173–183
(2005).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Okabe2001a">

:sup:`63` T. Okabe, M. Kawata, Y. Okamoto, and M. Mikami,
“Replica-exchange Monte Carlo method for the isobaric-isothermal
ensemble,” *Chem. Phys. Lett.*, **335** 435–439 (2001).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Chodera2011">

:sup:`64` J.D. Chodera and M.R. Shirts, “Replica exchange and expanded
ensemble simulations as gibbs sampling: Simple improvements for enhanced
mixing,” *J. Chem. Phys.*, **135** 194110 (2011).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Degroot96a">

:sup:`65` B.L. de Groot, A. Amadei, D.M.F. van Aalten, and H.J.C.
Berendsen, “Towards an exhaustive sampling of the configurational spaces
of the two forms of the peptide hormone guanylin,” *J. Biomol. Str.
Dyn.*, **13** [5] 741–751 (1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Degroot96b">

:sup:`66` B.L. de Groot, A. Amadei, R.M. Scheek, N.A.J. van Nuland, and
H.J.C. Berendsen, “An extended sampling of the configurational space of
HPr from *E. coli*,” *PROTEINS: Struct. Funct. Gen.*, **26** 314–322
(1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lange2006a">

:sup:`67` O.E. Lange, L.V. Schafer, and H. Grubmuller, “Flooding in
GROMACS: Accelerated barrier crossings in molecular dynamics,” *J. Comp.
Chem.*, **27** 1693–1702 (2006).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lyubartsev1992">

:sup:`68` A.P. Lyubartsev, A.A. Martsinovski, S.V. Shevkunov, and P.N.
Vorontsov-Velyaminov, “New approach to Monte Carlo calculation of the
free energy: Method of expanded ensembles,” *J. Chem. Phys.*, **96**
1776–1783 (1992).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Liem1991">

:sup:`69` S.Y. Liem, D. Brown, and J.H.R. Clarke, “Molecular dynamics
simulations on distributed memory machines,” *Comput. Phys. Commun.*,
**67** [2] 261–267 (1991).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Shaw2006">

:sup:`70` K.J. Bowers, R.O. Dror, and D.E. Shaw, “The midpoint method
for parallelization of particle simulations,” *J. Chem. Phys.*, **124**
[18] 184109–184109 (2006).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Tironi95">

:sup:`71` I.G. Tironi, R. Sperb, P.E. Smith, and W.F. van Gunsteren, “A
generalized reaction field method for molecular dynamics simulations,”
*J. Chem. Phys.*, **102** 5451–5459 (1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Spoel2006a">

:sup:`72` D. van der Spoel and P.J. van Maaren, “The origin of layer
structure artifacts in simulations of liquid water,” *J. Chem. Theory
Comput.*, **2** 1–11 (2006).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Ohmine1988">

:sup:`73` I. Ohmine, H. Tanaka, and P.G. Wolynes, “Large local energy
fluctuations in water. II. Cooperative motions and fluctuations,” *J.
Chem. Phys.*, **89** 5852–5860 (1988).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Kitchen1990">

:sup:`74` D.B. Kitchen, F. Hirata, J.D. Westbrook, R. Levy, D. Kofke,
and M. Yarmush, “Conserving energy during molecular dynamics simulations
of water, proteins, and proteins in water,” *J. Comp. Chem.*, **11**
1169–1180 (1990).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Guenot1993">

:sup:`75` J. Guenot and P.A. Kollman, “Conformational and energetic
effects of truncating nonbonded interactions in an aqueous protein
dynamics simulation,” *J. Comp. Chem.*, **14** 295–311 (1993).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Steinbach1994">

:sup:`76` P.J. Steinbach and B.R. Brooks, “New spherical-cutoff methods
for long-range forces in macromolecular simulation,” *J. Comp. Chem.*,
**15** 667–683 (1994).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-gromos96">

:sup:`77` W.F. van Gunsteren, S.R. Billeter, A.A. Eising, P.H.
Hünenberger, P. Krüger, A.E. Mark, W.R.P. Scott, and I.G. Tironi,
*Biomolecular simulation: The GROMOS96 manual and user guide*.
Hochschulverlag AG an der ETH Zürich, Zürich, Switzerland, 1996.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-biomos">

:sup:`78` W.F. van Gunsteren and H.J.C. Berendsen, *Gromos-87 manual*.
Biomos BV, Nijenborgh 4, 9747 AG Groningen, The Netherlands, 1987.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Morse29">

:sup:`79` P.M. Morse, “Diatomic molecules according to the wave
mechanics. II. vibrational levels.” *Phys. Rev.*, **34** 57–64 (1929).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Berendsen81">

:sup:`80` H.J.C. Berendsen, J.P.M. Postma, W.F. van Gunsteren, and J.
Hermans, “Interaction models for water in relation to protein
hydration”; pp. 331–342 in *Intermolecular forces*. Edited by B.
Pullman. D. Reidel Publishing Company, Dordrecht, 1981.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Ferguson95">

:sup:`81` D.M. Ferguson, “Parametrization and evaluation of a flexible
water model,” *J. Comp. Chem.*, **16** 501–511 (1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Warner72">

:sup:`82` H.R. Warner Jr., “Kinetic theory and rheology of dilute
suspensions of finitely extendible dumbbells,” *Ind. Eng. Chem.
Fundam.*, **11** [3] 379–387 (1972).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-MonicaGoga2013">

:sup:`83` M. Bulacu, N. Goga, W. Zhao, G. Rossi, L. Monticelli, X.
Periole, D. Tieleman, and S. Marrink, “Improved angle potentials for
coarse-grained molecular dynamics simulations,” *J. Chem. Phys.*,
**123** [11] (2005).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-BBrooks83">

:sup:`84` B.R. Brooks, R.E. Bruccoleri, B.D. Olafson, D.J. States, S.
Swaminathan, and M. Karplus, “CHARMM: A program for macromolecular
energy, minimization, and dynamics calculation,” *J. Comp. Chem.*, **4**
187–217 (1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lawrence2003b">

:sup:`85` C.P. Lawrence and J.L. Skinner, “Flexible TIP4P model for
molecular dynamics simulation of liquid water,” *Chem. Phys. Lett.*,
**372** 842–847 (2003).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Jorgensen1996">

:sup:`86` W.L. Jorgensen, D.S. Maxwell, and J. Tirado-Rives,
“Development and testing of the oPLS all-atom force field on
conformational energetics and properties of organic liquids,” *J. Am.
Chem. Soc.*, **118** 11225–11236 (1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Robertson2015a">

:sup:`87` M.J. Robertson, J. Tirado-Rives, and W.L. Jorgensen, “Improved
peptide and protein torsional energetics with the oPLS-aA force field,”
*J. Chem. Theory Comput.*, **11** 3499–3509 (2015).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-BulacuGiessen2005">

:sup:`88` M. Bulacu and E. van der Giessen, “Effect of bending and
torsion rigidity on self-diffusion in polymer melts: A
molecular-dynamics study,” *JCTC*, **9** [8] 3282–3292 (2013).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-ScottScheragator1966">

:sup:`89` R.A. Scott and H. Scheraga, “Conformational analysis of
macromolecules,” *J. Chem. Phys.*, **44** 3054–3069 (1966).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-PaulingBond">

:sup:`90` L. Pauling, *The nature of chemical bond*. Cornell University
Press, Ithaca; New York, 1960.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Torda89">

:sup:`91` A.E. Torda, R.M. Scheek, and W.F. van Gunsteren,
“Time-dependent distance restraints in molecular dynamics simulations,”
*Chem. Phys. Lett.*, **157** 289–294 (1989).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hess2003">

:sup:`92` B. Hess and R.M. Scheek, “Orientation restraints in molecular
dynamics simulations using time and ensemble averaging,” *J. Magn.
Reson.*, **164** 19–27 (2003).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lopes2013a">

:sup:`93` P.E.M. Lopes, J. Huang, J. Shim, Y. Luo, H. Li, B. Roux, and
J. MacKerell Alexander D., “Polarizable force field for peptides and
proteins based on the classical drude oscillator,” *J. Chem. Theory
Comput*, **9** 5430–5449 (2013).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-HYu2010">

:sup:`94` H. Yu, T.W. Whitfield, E. Harder, G. Lamoureux, I. Vorobyov,
V.M. Anisimov, A.D. MacKerell, Jr., and B. Roux, “Simulating Monovalent
and Divalent Ions in Aqueous Solution Using a Drude Polarizable Force
Field,” *J. Chem. Theory Comput.*, **6** 774–786 (2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Thole81">

:sup:`95` B.T. Thole, “Molecular polarizabilities with a modified dipole
interaction,” *Chem. Phys.*, **59** 341–345 (1981).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lamoureux2003a">

:sup:`96` G. Lamoureux and B. Roux, “Modeling induced polarization with
classical drude oscillators: Theory and molecular dynamics simulation
algorithm,” *J. Chem. Phys.*, **119** 3025–3039 (2003).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lamoureux2003b">

:sup:`97` G. Lamoureux, A.D. MacKerell, and B. Roux, “A simple
polarizable model of water based on classical drude oscillators,” *J.
Chem. Phys.*, **119** 5185–5197 (2003).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Noskov2005a">

:sup:`98` S.Y. Noskov, G. Lamoureux, and B. Roux, “Molecular dynamics
study of hydration in ethanol-water mixtures using a polarizable force
field,” *J. Phys. Chem. B.*, **109** 6705–6713 (2005).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Gunsteren98a">

:sup:`99` W.F. van Gunsteren and A.E. Mark, “Validation of molecular
dynamics simulations,” *J. Chem. Phys.*, **108** 6109–6116 (1998).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Beutler94">

:sup:`100` T.C. Beutler, A.E. Mark, R.C. van Schaik, P.R. Greber, and
W.F. van Gunsteren, “Avoiding singularities and numerical instabilities
in free energy calculations based on molecular simulations,” *Chem.
Phys. Lett.*, **222** 529–539 (1994).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Pham2011">

:sup:`101` T.T. Pham and M.R. Shirts, “Identifying low variance pathways
for free energy calculations of molecular transformations in solution
phase,” *J. Chem. Phys.*, **135** 034114 (2011).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Pham2012">

:sup:`102` T.T. Pham and M.R. Shirts, “Optimal pairwise and non-pairwise
alchemical pathways for free energy calculations of molecular
transformation in solution phase,” *J. Chem. Phys.*, **136** 124120
(2012).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Jorgensen88">

:sup:`103` W.L. Jorgensen and J. Tirado-Rives, “The OPLS potential
functions for proteins. energy minimizations for crystals of cyclic
peptides and crambin,” *J. Am. Chem. Soc.*, **110** 1657–1666 (1988).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Berendsen84b">

:sup:`104` H.J.C. Berendsen and W.F. van Gunsteren, “Molecular dynamics
simulations: Techniques and approaches”; pp. 475–500 in *Molecular
liquids-dynamics and interactions*. Edited by A.J.B. et al. Reidel,
Dordrecht, The Netherlands, 1984.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Ewald21">

:sup:`105` P.P. Ewald, “Die Berechnung optischer und elektrostatischer
Gitterpotentiale,” *Ann. Phys.*, **64** 253–287 (1921).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hockney81">

:sup:`106` R.W. Hockney and J.W. Eastwood, *Computer simulation using
particles*. McGraw-Hill, New York, 1981.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Ballenegger2012">

:sup:`107` V. Ballenegger, J.J. Cerdà, and C. Holm, “How to convert SPME
to P3M: Influence functions and error estimates,” *J. Chem. Theory
Comput.*, **8** [3] 936–947 (2012).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Allen87">

:sup:`108` M.P. Allen and D.J. Tildesley, *Computer simulations of
liquids*. Oxford Science Publications, Oxford, 1987.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Wennberg13">

:sup:`109` C.L. Wennberg, T. Murtola, B. Hess, and E. Lindahl,
“Lennard-Jones Lattice Summation in Bilayer Simulations Has Critical
Effects on Surface Tension and Lipid Properties,” *J. Chem. Theory
Comput.*, **9** 3527–3537 (2013).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Oostenbrink2004">

:sup:`110` C. Oostenbrink, A. Villa, A.E. Mark, and W.F. Van Gunsteren,
“A biomolecular force field based on the free enthalpy of hydration and
solvation: The GROMOS force-field parameter sets 53A5 and 53A6,”
*Journal of Computational Chemistry*, **25** [13] 1656–1676 (2004).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Cornell1995">

:sup:`111` W.D. Cornell, P. Cieplak, C.I. Bayly, I.R. Gould, K.R. Merz
Jr., D.M. Ferguson, D.C. Spellmeyer, and T. Fox *et al.*, “A Second
Generation Force Field for the Simulation of Proteins, Nucleic Acids,
and Organic Molecules,” *J. Am. Chem. Soc.*, **117** [19] 5179–5197
(1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Kollman1996">

:sup:`112` P.A. Kollman, “Advances and Continuing Challenges in
Achieving Realistic and Predictive Simulations of the Properties of
Organic and Biological Molecules,” *Acc. Chem. Res.*, **29** [10]
461–469 (1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Wang2000">

:sup:`113` J. Wang, P. Cieplak, and P.A. Kollman, “How Well Does a
Restrained Electrostatic Potential (RESP) Model Perform in Calculating
Conformational Energies of Organic and Biological Molecules?” *J. Comp.
Chem.*, **21** [12] 1049–1074 (2000).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hornak2006">

:sup:`114` V. Hornak, R. Abel, A. Okur, B. Strockbine, A. Roitberg, and
C. Simmerling, “Comparison of Multiple Amber Force Fields and
Development of Improved Protein Backbone Parameters,” *PROTEINS: Struct.
Funct. Gen.*, **65** 712–725 (2006).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Lindorff2010">

:sup:`115` K. Lindorff-Larsen, S. Piana, K. Palmo, P. Maragakis, J.L.
Klepeis, R.O. Dorr, and D.E. Shaw, “Improved side-chain torsion
potentials for the AMBER ff99SB protein force field,” *PROTEINS: Struct.
Funct. Gen.*, **78** 1950–1958 (2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Duan2003">

:sup:`116` Y. Duan, C. Wu, S. Chowdhury, M.C. Lee, G. Xiong, W. Zhang,
R. Yang, and P. Cieplak *et al.*, “A Point-Charge Force Field for
Molecular Mechanics Simulations of Proteins Based on Condensed-Phase
Quantum Mechanical Calculations,” *J. Comp. Chem.*, **24** [16]
1999–2012 (2003).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Garcia2002">

:sup:`117` A.E. García and K.Y. Sanbonmatsu, “\ :math:`\alpha`-Helical
stabilization by side chain shielding of backbone hydrogen bonds,”
*Proc. Natl. Acad. Sci. USA*, **99** [5] 2782–2787 (2002).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-mackerell04">

:sup:`118` J. MacKerell A. D., M. Feig, and C.L. Brooks III, “Extending
the treatment of backbone energetics in protein force fields:
Limitations of gas-phase quantum mechanics in reproducing protein
conformational distributions in molecular dynamics simulations,” *J.
Comp. Chem.*, **25** [11] 1400–15 (2004).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-mackerell98">

:sup:`119` A.D. MacKerell, D. Bashford, Bellott, R.L. Dunbrack, J.D.
Evanseck, M.J. Field, S. Fischer, and J. Gao *et al.*, “All-atom
empirical potential for molecular modeling and dynamics studies of
proteins,” *J. Phys. Chem. B.*, **102** [18] 3586–3616 (1998).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-feller00">

:sup:`120` S.E. Feller and A.D. MacKerell, “An improved empirical
potential energy function for molecular simulations of phospholipids,”
*J. Phys. Chem. B.*, **104** [31] 7510–7515 (2000).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-foloppe00">

:sup:`121` N. Foloppe and A.D. MacKerell, “All-atom empirical force
field for nucleic acids: I. Parameter optimization based on small
molecule and condensed phase macromolecular target data,” *J. Comp.
Chem.*, **21** [2] 86–104 (2000).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Mac2000">

:sup:`122` A.D. MacKerell and N.K. Banavali, “All-atom empirical force
field for nucleic acids: II. application to molecular dynamics
simulations of DNA and RNA in solution,” *J. Comp. Chem.*, **21** [2]
105–120 (2000).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Larsson10">

:sup:`123` P. Larsson and E. Lindahl, “A High-Performance
Parallel-Generalized Born Implementation Enabled by Tabulated
Interaction Rescaling,” *J. Comp. Chem.*, **31** [14] 2593–2600 (2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Bjelkmar10">

:sup:`124` P. Bjelkmar, P. Larsson, M.A. Cuendet, B. Hess, and E.
Lindahl, “Implementation of the CHARMM force field in GROMACS: Analysis
of protein stability effects from correction maps, virtual interaction
sites, and water models,” *J. Chem. Theory Comput.*, **6** 459–466
(2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-kohlmeyer2016">

:sup:`125` A. Kohlmeyer and J. Vermaas, *TopoTools: Release 1.6 with
CHARMM export in topogromacs*, (2016).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-bereau12">

:sup:`126` T. Bereau, Z.-J. Wang, and M. Deserno, *Solvent-free
coarse-grained model for unbiased high-resolution protein-lipid
interactions*, (n.d.).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-wang_jpcb10">

:sup:`127` Z.-J. Wang and M. Deserno, “A systematically coarse-grained
solvent-free model for quantitative phospholipid bilayer simulations,”
*J. Phys. Chem. B.*, **114** [34] 11207–11220 (2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Jorgensen83">

:sup:`128` W.L. Jorgensen, J. Chandrasekhar, J.D. Madura, R.W. Impey,
and M.L. Klein, “Comparison of simple potential functions for simulating
liquid water,” *J. Chem. Phys.*, **79** 926–935 (1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-iupac70">

:sup:`129` IUPAC-IUB Commission on Biochemical Nomenclature,
“Abbreviations and Symbols for the Description of the Conformation of
Polypeptide Chains. Tentative Rules (1969),” *Biochemistry*, **9**
3471–3478 (1970).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Mahoney2000a">

:sup:`130` M.W. Mahoney and W.L. Jorgensen, “A five-site model for
liquid water and the reproduction of the density anomaly by rigid,
nonpolarizable potential functions,” *J. Chem. Phys.*, **112** 8910–8922
(2000).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Ryckaert78">

:sup:`131` J.P. Ryckaert and A. Bellemans, “Molecular dynamics of liquid
alkanes,” *Far. Disc. Chem. Soc.*, **66** 95–106 (1978).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Loof92">

:sup:`132` H. de Loof, L. Nilsson, and R. Rigler, “Molecular dynamics
simulations of galanin in aqueous and nonaqueous solution,” *J. Am.
Chem. Soc.*, **114** 4028–4035 (1992).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Buuren93a">

:sup:`133` A.R. van Buuren and H.J.C. Berendsen, “Molecular Dynamics
simulation of the stability of a 22 residue alpha-helix in water and 30%
trifluoroethanol,” *Biopolymers*, **33** 1159–1166 (1993).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-RMNeumann1980a">

:sup:`134` R.M. Neumann, “Entropic approach to Brownian Movement,” *Am.
J. Phys.*, **48** 354–357 (1980).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Jarzynski1997a">

:sup:`135` C. Jarzynski, “Nonequilibrium equality for free energy
differences,” *Phys. Rev. Lett.*, **78** [14] 2690–2693 ().

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Engin2010a">

:sup:`136` M.S. O. Engin A. Villa and B. Hess, “Driving forces for
adsorption of amphiphilic peptides to air-water interface,” *J. Phys.
Chem. B.*, (2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-lindahl2014accelerated">

:sup:`137` V. Lindahl, J. Lidmar, and B. Hess, “Accelerated weight
histogram method for exploring free energy landscapes,” *The Journal of
chemical physics*, **141** [4] 044110 (2014).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-wang2001efficient">

:sup:`138` F. Wang and D. Landau, “Efficient, multiple-range random walk
algorithm to calculate the density of states,” *Physical review
letters*, **86** [10] 2050 (2001).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-huber1994local">

:sup:`139` T. Huber, A.E. Torda, and W.F. van Gunsteren, “Local
elevation: A method for improving the searching properties of molecular
dynamics simulation,” *Journal of computer-aided molecular design*,
**8** [6] 695–708 (1994).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-laio2002escaping">

:sup:`140` A. Laio and M. Parrinello, “Escaping free-energy minima,”
*Proceedings of the National Academy of Sciences*, **99** [20]
12562–12566 (2002).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-belardinelli2007fast">

:sup:`141` R. Belardinelli and V. Pereyra, “Fast algorithm to calculate
density of states,” *Physical Review E*, **75** [4] 046701 (2007).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-barducci2008well">

:sup:`142` A. Barducci, G. Bussi, and M. Parrinello, “Well-tempered
metadynamics: A smoothly converging and tunable free-energy method,”
*Physical review letters*, **100** [2] 020603 (2008).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-lindahl2017sequence">

:sup:`143` V. Lindahl, A. Villa, and B. Hess, “Sequence dependency of
canonical base pair opening in the dNA double helix,” *PLoS
computational biology*, **13** [4] e1005463 (2017).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-sivak2012thermodynamic">

:sup:`144` D.A. Sivak and G.E. Crooks, “Thermodynamic metrics and
optimal paths,” *Physical review letters*, **108** [19] 190602 (2012).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Kutzner2011">

:sup:`145` C. Kutzner, J. Czub, and H. Grubmüller, “Keep it flexible:
Driving macromolecular rotary motions in atomistic simulations with
GROMACS,” *J. Chem. Theory Comput.*, **7** 1381–1393 (2011).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Caleman2008a">

:sup:`146` C. Caleman and D. van der Spoel, “Picosecond Melting of Ice
by an Infrared Laser Pulse - A simulation study,” *Angew. Chem., Int.
Ed. Engl.*, **47** 1417–1420 (2008).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Kutzner2011b">

:sup:`147` C. Kutzner, H. Grubmüller, B.L. de Groot, and U. Zachariae,
“Computational electrophysiology: The molecular dynamics of ion channel
permeation and selectivity in atomistic detail,” *Biophys. J.*, **101**
809–817 (2011).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-feenstra99">

:sup:`148` K.A. Feenstra, B. Hess, and H.J.C. Berendsen, “Improving
efficiency of large time-scale molecular dynamics simulations of
hydrogen-rich systems,” *J. Comp. Chem.*, **20** 786–798 (1999).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hess2002a">

:sup:`149` B. Hess, “Determining the shear viscosity of model liquids
from molecular dynamics,” *J. Chem. Phys.*, **116** 209–217 (2002).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-mopac">

:sup:`150` M.J.S. Dewar, “Development and status of MINDO/3 and MNDO,”
*J. Mol. Struct.*, **100** 41 (1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-gamess-uk">

:sup:`151` M.F. Guest, R.J. Harrison, J.H. van Lenthe, and L.C.H. van
Corler, “Computational chemistry on the FPS-X64 scientific computers -
Experience on single- and multi-processor systems,” *Theor. Chim. Act.*,
**71** 117 (1987).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-g03">

:sup:`152` M.J. Frisch, G.W. Trucks, H.B. Schlegel, G.E. Scuseria, M.A.
Robb, J.R. Cheeseman, J.A. Montgomery Jr., and T. Vreven *et al.*,
*Gaussian 03, Revision C.02*, (n.d.).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Car85a">

:sup:`153` R. Car and M. Parrinello, “Unified approach for molecular
dynamics and density-functional theory,” *Phys. Rev. Lett.*, **55**
2471–2474 (1985).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Field90a">

:sup:`154` M. Field, P.A. Bash, and M. Karplus, “A combined quantum
mechanical and molecular mechanical potential for molecular dynamics
simulation,” *J. Comp. Chem.*, **11** 700 (1990).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Maseras96a">

:sup:`155` F. Maseras and K. Morokuma, “IMOMM: A New Ab Initio +
Molecular Mechanics Geometry Optimization Scheme of Equilibrium
Structures and Transition States,” *J. Comp. Chem.*, **16** 1170–1179
(1995).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Svensson96a">

:sup:`156` M. Svensson, S. Humbel, R.D.J. Froes, T. Matsubara, S.
Sieber, and K. Morokuma, “ONIOM a multilayered integrated MO + MM method
for geometry optimizations and single point energy predictions. a test
for Diels-Alder reactions and Pt(P(t-Bu)3)2 + H2 oxidative addition,”
*J. Phys. Chem.*, **100** 19357 (1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Yesylevskyy2007">

:sup:`157` S. Yesylevskyy, “ProtSqueeze: Simple and effective automated
tool for setting up membrane protein simulations,” *J. Chem. Inf.
Model.*, **47** 1986–1994 (2007).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Wolf2010">

:sup:`158` M. Wolf, M. Hoefling, C. Aponte-Santamaría, H. Grubmüller,
and G. Groenhof, “g\_membed: Efficient insertion of a membrane protein
into an equilibrated lipid bilayer with minimal perturbation,” *J. Comp.
Chem.*, **31** 2169–2174 (2010).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Spoel97a">

:sup:`159` D. van der Spoel and H.J.C. Berendsen, “Molecular dynamics
simulations of Leu-enkephalin in water and DMSO,” *Biophys. J.*, **72**
2032–2041 (1997).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-PSmith93c">

:sup:`160` P.E. Smith and W.F. van Gunsteren, “The Viscosity of SPC and
SPC/E Water,” *Chem. Phys. Lett.*, **215** 315–318 (1993).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Balasubramanian96">

:sup:`161` S. Balasubramanian, C.J. Mundy, and M.L. Klein, “Shear
viscosity of polar fluids: Miolecular dynamics calculations of water,”
*J. Chem. Phys.*, **105** 11190–11195 (1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-lmfit">

:sup:`162` J. Wuttke, *Lmfit*, (2013).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Steen-Saethre2014a">

:sup:`163` B. Steen-Sæthre, A.C. Hoffmann, and D. van der Spoel, “Order
parameters and algorithmic approaches for detection and demarcation of
interfaces in hydrate-fluid and ice-fluid systems,” *J. Chem. Theor.
Comput.*, **10** 5606–5615 (2014).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Palmer1994a">

:sup:`164` B.J. Palmer, “Transverse-current autocorrelation-function
calculations of the shear viscosity for molecular liquids.” *Phys. Rev.
E*, **49** 359–366 (1994).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Wensink2003a">

:sup:`165` E.J.W. Wensink, A.C. Hoffmann, P.J. van Maaren, and D. van
der Spoel, “Dynamic properties of water/alcohol mixtures studied by
computer simulation,” *J. Chem. Phys.*, **119** 7308–7317 (2003).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Guo2002b">

:sup:`166` G.-J. Guo, Y.-G. Zhang, K. Refson, and Y.-J. Zhao, “Viscosity
and stress autocorrelation function in supercooled water: A molecular
dynamics study,” *Mol. Phys.*, **100** 2617–2627 (2002).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Fanourgakis2012a">

:sup:`167` G.S. Fanourgakis, J.S. Medina, and R. Prosmiti, “Determining
the bulk viscosity of rigid water models,” *J. Phys. Chem. A*, **116**
2564–2570 (2012).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Spoel96b">

:sup:`168` D. van der Spoel, H.J. Vogel, and H.J.C. Berendsen,
“Molecular dynamics simulations of N-terminal peptides from a nucleotide
binding protein,” *PROTEINS: Struct. Funct. Gen.*, **24** 450–466
(1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Amadei93">

:sup:`169` A. Amadei, A.B.M. Linssen, and H.J.C. Berendsen, “Essential
dynamics of proteins,” *PROTEINS: Struct. Funct. Gen.*, **17** 412–425
(1993).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hess2002b">

:sup:`170` B. Hess, “Convergence of sampling in protein simulations,”
*Phys. Rev. **E***, **65** 031910 (2002).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Hess2000">

:sup:`171` B. Hess, “Similarities between principal components of
protein dynamics and random diffusion,” *Phys. Rev. **E***, **62**
8438–8448 (2000).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Mu2005a">

:sup:`172` Y. Mu, P.H. Nguyen, and G. Stock, “Energy landscape of a
small peptide revelaed by dihedral angle principal component analysis,”
*PROTEINS: Struct. Funct. Gen.*, **58** 45–52 (2005).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Spoel2006b">

:sup:`173` D. van der Spoel, P.J. van Maaren, P. Larsson, and N.
Timneanu, “Thermodynamics of hydrogen bonding in hydrophilic and
hydrophobic media,” *J. Phys. Chem. B.*, **110** 4393–4398 (2006).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Luzar96b">

:sup:`174` A. Luzar and D. Chandler, “Hydrogen-bond kinetics in liquid
water,” *Nature*, **379** 55–57 (1996).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Luzar2000a">

:sup:`175` A. Luzar, “Resolving the hydrogen bond dynamics conundrum,”
*J. Chem. Phys.*, **113** 10663–10675 (2000).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Kabsch83">

:sup:`176` W. Kabsch and C. Sander, “Dictionary of protein secondary
structure: Pattern recognition of hydrogen-bonded and geometrical
features,” *Biopolymers*, **22** 2577–2637 (1983).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Bekker93b">

:sup:`177` H. Bekker, H.J.C. Berendsen, E.J. Dijkstra, S. Achterop, R.
v. Drunen, D. v. d. Spoel, A. Sij:raw-latex:`\-`bers, and H. Keegstra
*et al.*, “Gromacs Method of Virial Calculation Using a Single Sum”; pp.
257–261 in *Physics computing 92*. Edited by R.A. de Groot and J.
Nadrchal. World Scientific, Singapore, 1993.

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Berendsen87">

:sup:`178` H.J.C. Berendsen, J.R. Grigera, and T.P. Straatsma, “The
missing term in effective pair potentials,” *J. Phys. Chem.*, **91**
6269–6271 (1987).

.. raw:: html

   </div>

.. raw:: html

   <div id="ref-Gunsteren94a">

:sup:`179` W.F. van Gunsteren and H.J.C. Berendsen, *Molecular dynamics
of simple systems*, (1994).

.. raw:: html

   </div>

.. raw:: html

   </div>

.. [1]
   Note that some derivations, an alternative notation
   :math:`\xi_{\mathrm{alt}} = v_{\xi} = p_{\xi}/Q` is used.

.. [2]
   The box matrix representation in GROMACS corresponds to the transpose
   of the box matrix representation in the paper by Nosé and Klein.
   Because of this, some of our equations will look slightly different.

.. [3]
   :math:`P_0(x) = 1`, :math:`P_1(x) = x`, :math:`P_2(x) = (3x^2-1)/2`

.. [4]
   :math:`({\bf u}\otimes{\bf v})^{{\alpha\beta}}~=~{\bf u}_{{\alpha}}{\bf v}_{{\beta}}`

.. [5]
   The calculation of Lennard-Jones and Coulomb forces is about 50
   floating point operations.

.. |image| image:: plots/peregrine
   :height: 5.00000in
.. |A rhombic dodecahedron and truncated octahedron (arbitrary orientations).| image:: plots/rhododec
   :width: 5.00000cm
.. |A rhombic dodecahedron and truncated octahedron (arbitrary orientations).| image:: plots/truncoct
   :width: 5.00000cm
.. |Free energy cycles. **A:** to calculate :math:`\Delta G_{12}`, the free energy difference between the binding of inhibitor **I** to enzymes **E** respectively **E\ :math:`^{\prime}`**. **B:** to calculate :math:`\Delta G_{12}`, the free energy difference for binding of inhibitors **I** respectively **I\ :math:`^{\prime}`** to enzyme **E**.| image:: plots/free1
   :width: 6.00000cm
.. |Free energy cycles. **A:** to calculate :math:`\Delta G_{12}`, the free energy difference between the binding of inhibitor **I** to enzymes **E** respectively **E\ :math:`^{\prime}`**. **B:** to calculate :math:`\Delta G_{12}`, the free energy difference for binding of inhibitors **I** respectively **I\ :math:`^{\prime}`** to enzyme **E**.| image:: plots/free2
   :width: 6.00000cm
.. |Principle of improper dihedral angles. Out of plane bending for rings (left), substituents of rings (middle), out of tetrahedral (right). The improper dihedral angle :math:`\xi` is defined as the angle between planes (i,j,k) and (j,k,l) in all cases.| image:: plots/ring-imp
   :width: 4.00000cm
.. |Principle of improper dihedral angles. Out of plane bending for rings (left), substituents of rings (middle), out of tetrahedral (right). The improper dihedral angle :math:`\xi` is defined as the angle between planes (i,j,k) and (j,k,l) in all cases.| image:: plots/subst-im
   :width: 3.00000cm
.. |Principle of improper dihedral angles. Out of plane bending for rings (left), substituents of rings (middle), out of tetrahedral (right). The improper dihedral angle :math:`\xi` is defined as the angle between planes (i,j,k) and (j,k,l) in all cases.| image:: plots/tetra-im
   :width: 3.00000cm
.. |AWH evolution in time for a Brownian particle in a double-well potential. The reaction coordinate :math:`\xi(t)` traverses the sampling interval multiple times in the initial stage before exiting and entering the final stage (top left). In the final stage, the dynamics of :math:`\xi` becomes increasingly diffusive. The times of covering are shown as :math:`\times`-markers of different colors. At these times the free energy update size :math:`\sim 1/N`, where :math:`N` is the size of the weight histogram, is decreased by scaling :math:`N` by a factor of :math:`\gamma=3` (top right). In the final stage, :math:`N` grows at the sampling rate and thus :math:`1/N\sim1/t`. The exit from the final stage is determined on the fly by ensuring that the effective sample weight :math:`s` of data collected in the final stage exceeds that of initial stage data (bottom left; note that :math:`\ln s(t)` is plotted). An estimate of the PMF is also extracted from the simulation (bottom right), which after exiting the initial stage should estimate global free energy differences fairly accurately. | image:: plots/awh-traj
   :width: 6.00000cm
.. |AWH evolution in time for a Brownian particle in a double-well potential. The reaction coordinate :math:`\xi(t)` traverses the sampling interval multiple times in the initial stage before exiting and entering the final stage (top left). In the final stage, the dynamics of :math:`\xi` becomes increasingly diffusive. The times of covering are shown as :math:`\times`-markers of different colors. At these times the free energy update size :math:`\sim 1/N`, where :math:`N` is the size of the weight histogram, is decreased by scaling :math:`N` by a factor of :math:`\gamma=3` (top right). In the final stage, :math:`N` grows at the sampling rate and thus :math:`1/N\sim1/t`. The exit from the final stage is determined on the fly by ensuring that the effective sample weight :math:`s` of data collected in the final stage exceeds that of initial stage data (bottom left; note that :math:`\ln s(t)` is plotted). An estimate of the PMF is also extracted from the simulation (bottom right), which after exiting the initial stage should estimate global free energy differences fairly accurately. | image:: plots/awh-invN
   :width: 6.00000cm
.. |AWH evolution in time for a Brownian particle in a double-well potential. The reaction coordinate :math:`\xi(t)` traverses the sampling interval multiple times in the initial stage before exiting and entering the final stage (top left). In the final stage, the dynamics of :math:`\xi` becomes increasingly diffusive. The times of covering are shown as :math:`\times`-markers of different colors. At these times the free energy update size :math:`\sim 1/N`, where :math:`N` is the size of the weight histogram, is decreased by scaling :math:`N` by a factor of :math:`\gamma=3` (top right). In the final stage, :math:`N` grows at the sampling rate and thus :math:`1/N\sim1/t`. The exit from the final stage is determined on the fly by ensuring that the effective sample weight :math:`s` of data collected in the final stage exceeds that of initial stage data (bottom left; note that :math:`\ln s(t)` is plotted). An estimate of the PMF is also extracted from the simulation (bottom right), which after exiting the initial stage should estimate global free energy differences fairly accurately. | image:: plots/awh-sampleweights
   :width: 6.00000cm
.. |AWH evolution in time for a Brownian particle in a double-well potential. The reaction coordinate :math:`\xi(t)` traverses the sampling interval multiple times in the initial stage before exiting and entering the final stage (top left). In the final stage, the dynamics of :math:`\xi` becomes increasingly diffusive. The times of covering are shown as :math:`\times`-markers of different colors. At these times the free energy update size :math:`\sim 1/N`, where :math:`N` is the size of the weight histogram, is decreased by scaling :math:`N` by a factor of :math:`\gamma=3` (top right). In the final stage, :math:`N` grows at the sampling rate and thus :math:`1/N\sim1/t`. The exit from the final stage is determined on the fly by ensuring that the effective sample weight :math:`s` of data collected in the final stage exceeds that of initial stage data (bottom left; note that :math:`\ln s(t)` is plotted). An estimate of the PMF is also extracted from the simulation (bottom right), which after exiting the initial stage should estimate global free energy differences fairly accurately. | image:: plots/awh-pmfs
   :width: 6.00000cm
.. |The window of gmx view showing a box of water.| image:: plots/ngmxdump
   :width: 8.00000cm
.. |Definition of slices in gmx rdf: A. :math:`g_{AB}(r)`. B. :math:`g_{AB}(r,\theta)`. The slices are colored gray. C. Normalization :math:`\langle\rho_B\rangle_{local}`. D. Normalization :math:`\langle\rho_B\rangle_{local,\:\theta }`. Normalization volumes are colored gray.| image:: plots/rdf
   :width: 7.00000cm
.. |:math:`g_{OO}(r)` for Oxygen-Oxygen of SPC-water.| image:: plots/rdfO-O
   :width: 8.00000cm
.. |Mean Square Displacement of SPC-water.| image:: plots/msdwater
   :width: 8.00000cm
.. |Dihedral conventions: A. “Biochemical convention”. B. “Polymer convention”.| image:: plots/dih-def
   :width: 3.50000cm
.. |Angle options of gmx gangle: A. Angle between two vectors. B. Angle between two planes. C. Angle between a vector and the :math:`z` axis. D. Angle between a vector and the normal of a sphere. Also other combinations are supported: planes and vectors can be used interchangeably.| image:: plots/sgangle
.. |Insertion of water into an H-bond. (1) Normal H-bond between two residues. (2) H-bonding bridge via a water molecule.| image:: plots/hbond-insert
   :width: 5.00000cm
.. |Ramachandran plot of a small protein.| image:: plots/rama
   :width: 8.00000cm
