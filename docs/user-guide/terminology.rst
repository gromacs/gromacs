Terminology
===========

.. _gmx-pressure:

Pressure
--------

The pressure in molecular dynamics can be computed from the kinetic energy and
the virial. 

Fluctuation
^^^^^^^^^^^

Whether or not pressure coupling is used within a simulation, the pressure
value for the simulation box will oscillate significantly. Instantaneous
pressure is meaningless, and not well-defined. Over a picosecond time scale it
usually will not be a good indicator of the true pressure. This variation is
entirely normal due to the fact that pressure is a macroscopic property and can
only be measured properly as time average, while it is being measured and/or
adjusted with pressure coupling on the microscopic scale. How much it varies
and the speed at which it does depends on the number of atoms in the system,
the type of pressure coupling used and the value of the coupling constants.
Fluctuations of the order of hundreds of bar are typical. For a box of 216
waters, fluctuations of 500-600 bar are standard. Since the fluctuations go
down with the square root of the number of particles, a system of 21600 water
molecules (100 times larger) will still have pressure fluctuations of 50-60 bar.

.. _gmx-pbc:

Periodic boundary conditions
----------------------------

Periodic boundary conditions (PBC) are used in molecular dynamics simulations
to avoid problems with boundary effects caused by finite size, and make the
system more like an infinite one, at the cost of possible periodicity effects.

Beginners visualizing a trajectory sometimes think they are observing a problem
when

* the molecule(s) does not stay in the centre of the box, or
* it appears that (parts of) the molecule(s) diffuse out of the box, or
* holes are created, or
* broken molecules appear, or
* their unit cell was a rhombic dodecahedron or cubic octahedron but it looks
  like a slanted cube after the simulation, or
* crazy bonds all across the simulation cell appear.

This is not a problem or error that is occurring, it is what you should expect.

The existence of PBC means that any atom that leaves a simulation box by, say,
the right-hand face, then enters the simulation box by the left-hand face. In
the example of a large protein, if you look at the face of the simulation box
that is opposite to the one from which the protein is protruding, then a hole
in the solvent will be visible. The reason that the molecule(s) move from where
they were initially located within the box is (for the vast majority of
simulations) they are free to diffuse around. And so they do. They are not held
in a magic location of the box. The box is not centered around anything while
performing the simulation. Molecules are not made whole as a matter of course.
Moreover, any periodic cell shape can be expressed as a parallelepiped (a.k.a.
triclinic cell), and |Gromacs| does so internally regardless of the initial
shape of the box.

These visual issues can be fixed after the conclusion of the simulation by
judicious use of the optional inputs to :ref:`gmx trjconv` to process the
trajectory files. Similarly, analyses such as RMSD of atomic positions can be
flawed when a reference structure is compared with a structure that needs
adjusting for periodicity effects, and the solution with :ref:`gmx trjconv`
follows the same lines. Some complex cases needing more than one operation will
require more than one invocation of :ref:`gmx trjconv` in order to work.

For further information, see the corresponding section in the :ref:`Reference Manual <pbc>`.

Suggested workflow
^^^^^^^^^^^^^^^^^^

Fixing periodicity effects with :ref:`gmx trjconv` to suit visualization or
analysis can be tricky. Multiple invocations can be necessary. You may need to
create custom index groups (e.g. to keep your ligand with your protein)
Following the steps below in order (omitting those not required) should help
get a pleasant result. You will need to consult ``gmx trjconv -h`` to find out
the details for each step. That's deliberate -- there is no magic "do what I
want" recipe. You have to decide what you want, first. :-)

#. First make your molecules whole if you want them whole.
#. Cluster your molecules/particles if you want them clustered.
#. If you want jumps removed, extract the first frame from the trajectory to
   use as the reference, and then use ``-pbc nojump`` with that first
   frame as reference.
#. Center your system using some criterion. Doing so shifts the system, so
   don't use ``-pbc nojump`` after this step.
#. Perhaps put everything in some box with the other ``-pbc`` or ``-ur``
   options.
#. Fit the resulting trajectory to some (other) reference structure (if
   desired), and don't use any PBC related option afterwards.

With point three, the issue is that :ref:`gmx trjconv` removes the jumps from
the first frame using the reference structure provided with -s. If the reference
structure (run input file) is not clustered/whole, using ``-pbc nojump``
will undo steps 1 and 2.

.. _gmx-thermostats:

Thermostats
-----------

Thermostats are designed to help a simulation sample from the correct ensemble
(i.e. NVT or NPT) by modulating the temperature of the system in some fashion.
First, we need to establish what we mean by temperature. In simulations, the
"instantaneous (kinetic) temperature" is usually computed from the kinetic
energy of the system using the equipartition theorem. In other words, the
temperature is computed from the system's total kinetic energy.

So, what's the goal of a thermostat? Actually, it turns out the goal is not to
keep the temperature constant, as that would mean fixing the total kinetic
energy, which would be silly and not the aim of NVT or NPT. Rather, it's to
ensure that the average temperature of a system be correct.

To see why this is the case, imagine a glass of water sitting in a room.
Suppose you can look very closely at a few molecules in some small region of
the glass, and measure their kinetic energies. You would not expect the kinetic
energy of this small number of particles to remain precisely constant; rather,
you'd expect fluctuations in the kinetic energy due to the small number of
particles. As you average over larger and larger numbers of particles, the
fluctuations in the average get smaller and smaller, so finally by the time you
look at the whole glass, you say it has "constant temperature".

Molecular dynamics simulations are often fairly small compared to a glass of
water, so we have bigger fluctuations. So it's really more appropriate here to
think of the role of a thermostat as ensuring that we have

(a) the correct average temperature, and
(b) the fluctuations of the correct size.

See the relevant section in the :ref:`Reference Manual <temp-coupling>`
for details on how temperature coupling is applied and
the types currently available.

.. _gmx-thermostats-do:

What to do
^^^^^^^^^^

Some hints on practices that generally are a good idea:

* Preferably, use a thermostat that samples the correct distribution of
  temperatures (for examples, see the corresponding manual section), in addition
  to giving you the correct average temperature.
* At least: use a thermostat that gives you the correct average temperature,
  and apply it to components of your system for which they are justified (see
  the first bullet in `What not to do`_). In some cases, using
  ``tc-grps = System`` may lead to the "hot solvent/cold solute" problem
  described in the 3rd reference in `Further reading`_.

.. _gmx-thermostats-dont:

What not to do
^^^^^^^^^^^^^^

Some hints on practices that generally not a good idea to use:

* Do not use separate thermostats for every component of your system. Some
  molecular dynamics thermostats only work well in the thermodynamic limit. A
  group must be of sufficient size to justify its own thermostat. If you use one
  thermostat for, say, a small molecule, another for protein, and another for
  water, you are likely introducing errors and artifacts that are hard to
  predict. In particular, do not couple ions in aqueous solvent in a separate
  group from that solvent. For a protein simulation, using ``tc-grps = Protein
  Non-Protein`` is usually best.
* Do not use thermostats that work well only in the limit of a large number of
  degrees of freedom for systems with few degrees of freedom. For example, do
  not use Nosé-Hoover or Berendsen thermostats for types of free energy
  calculations where you will have a component of the system with very few
  degrees of freedom in an end state (i.e. a noninteracting small molecule).

Further reading
^^^^^^^^^^^^^^^

#. Cheng, A. & Merz, K. M. Application of the Nosé-Hoover chain algorithm to
   the study of protein dynamics. *J. Phys. Chem.* **100** (5), 1927–1937
   (`1996 <http://pubs.acs.org/doi/abs/10.1021/jp951968y>`_).
#. Mor, A., Ziv, G. & Levy, Y. Simulations of proteins with inhomogeneous
   degrees of freedom: the effect of thermostats. *J. Comput. Chem.* **29**
   (12), 1992–1998 (`2008 <https://doi.org/10.1002/jcc.20951>`_).
#. Lingenheil, M., Denschlag, R., Reichold, R. & Tavan, P. The
   "hot-solvent/cold-solute" problem revisited. *J. Chem. Theory Comput.* **4**
   (8), 1293–1306 (`2008 <http://pubs.acs.org/doi/abs/10.1021/ct8000365>`__).

Energy conservation
-------------------

In principle, a molecular dynamics simulation should conserve the total energy,
the total momentum and (in a non-periodic system) the total angular momentum. A
number of algorithmic and numerical issues make that this is not always the
case:

* Cut-off treatment and/or long-range electrostatics treatment (see Van Der
  Spoel, D. & van Maaren, P. J. The origin of layer structure artifacts in
  simulations of liquid water. *J. Chem. Theor. Comp.* **2**, 1–11
  (`2006 <https://doi.org/10.1021/ct0502256>`_).)
* Treatment of pair lists,
* Constraint algorithms (see e.g. Hess, B. P-LINCS: A parallel linear constraint
  solver for molecular simulation. *J. Chem. Theor. Comp.* **4**, 116–122
  (`2008 <https://doi.org/10.1021/ct700200b>`__).).
* The integration timestep.
* :ref:`Temperature coupling <gmx-thermostats>` and :ref:`pressure coupling <gmx-pressure>`.
* Round-off error (in particular in single precision), for example subtracting
  large numbers (Lippert, R. A. et al. A common, avoidable source of error in
  molecular dynamics integrators. *J. Chem. Phys.* **126**, 046101 (`2007 <https://doi.org/10.1063/1.2431176>`_).).
* The choice of the integration algorithm (in |Gromacs| this is normally
  leap-frog).
* Removal of center of mass motion: when doing this in more than one group the
  conservation of energy will be violated.

Average structure
-----------------

Various |Gromacs| utilities can compute average structures. Presumably the idea
for this comes from something like an ensemble-average NMR structure. In some
cases, it makes sense to calculate an average structure (as a step on the way
to calculating root-mean-squared fluctuations (RMSF), for example, one needs
the average position of all of the atoms).

However, it's important to remember that an average structure isn't necessarily
meaningful. By way of analogy, suppose I alternate holding a ball in my left
hand, then in my right hand. What's the average position of the ball? Halfway
in between -- even though I always have it either in my left hand or my right
hand. Similarly, for structures, averages will tend to be meaningless anytime
there are separate metastable conformational states. This can happen on a
sidechain level, or for some regions of backbone, or even whole helices or
components of the secondary structure.

Thus, if you derive an average structure from a molecular dynamics simulation,
and find artifacts like unphysical bond lengths, weird structures, etc., this
doesn't necessarily mean something is wrong. It just shows the above: an
average structure from a simulation is not necessarily a physically meaningful
structure.

.. _blowing-up:

Blowing up
----------

*Blowing up* is a highly technical term used to describe a common sort of
simulation failure. In brief, it describes a failure typically due to an
unacceptably large force that ends up resulting in a failure of the integrator.

To give a bit more background, it's important to remember that molecular
dynamics numerically integrates Newton's equations of motion by taking small,
discrete timesteps, and using these timesteps to determine new velocities and
positions from velocities, positions, and forces at the previous timestep. If
forces become too large at one timestep, this can result in extremely large
changes in velocity/position when going to the next timestep. Typically, this
will result in a cascade of errors: one atom experiences a very large force one
timestep, and thus goes shooting across the system in an uncontrolled way in
the next timestep, overshooting its preferred location or landing on top of
another atom or something similar. This then results in even larger forces the
next timestep, more uncontrolled motions, and so on. Ultimately, this will
cause the simulation package to crash in some way, since it can't cope with
such situations. In simulations with constraints, the first symptom of this
will usually be some LINCS or SHAKE warning or error -- not because the
constraints are the source of the problem, but just because they're the first
thing to crash. Similarly, in simulations with domain decomposition, you may
see messages about particles being more than a cell length out of the domain
decomposition cell of their charge group, which are symptomatic of your
underlying problem, and not the domain decomposition algorithm itself. Likewise
for warnings about tabulated or 1-4 interactions being outside the distance
supported by the table. This can happen on one computer system while another
resulted in a stable simulation because of the impossibility of numerical
reproducibility of these calculations on different computer systems.

Possible causes include:

* you didn't minimize well enough,
* you have a bad starting structure, perhaps with steric clashes,
* you are using too large a timestep (particularly given your choice of
  constraints),
* you are doing particle insertion in free energy calculations without using
  soft core,
* you are using inappropriate pressure coupling (e.g. when you are not in
  equilibrium, Berendsen can be best while relaxing the volume, but you will
  need to switch to a more accurate pressure-coupling algorithm later),
* you are using inappropriate temperature coupling, perhaps on inappropriate
  groups, or
* your position restraints are to coordinates too different from those present
  in the system, or
* you have a single water molecule somewhere within the system that is
  isolated from the other water molecules, or
* you are experiencing a bug in :ref:`gmx mdrun`.

Because blowing up is due, typically, to forces that are too large for a
particular timestep size, there are a couple of basic solutions:

* make sure the forces don't get that large, or
* use a smaller timestep.

Better system preparation is a way to make sure that forces don't get large, if
the problems are occurring near the beginning of a simulation.

.. _system-diagnosis:

Diagnosing an unstable system
-----------------------------

Troubleshooting a system that is blowing up can be challenging, especially for
an inexperienced user. Here are a few general tips that one may find useful
when addressing such a scenario:

#. If the crash is happening relatively early (within a few steps), set
   ``nstxout`` (or ``nstxout-compressed``) to 1, capturing all possible frames.
   Watch the resulting trajectory to see which atoms/residues/molecules become
   unstable first.
#. Simplify the problem to try to establish a cause:

   * If you have a new box of solvent, try minimizing and simulating a single
     molecule to see if the instability is due to some inherent problem with
     the molecule's topology or if instead there are clashes in your starting
     configuration.
   * If you have a protein-ligand system, try simulating the protein alone in
     the desired solvent. If it is stable, simulate the ligand in vacuo to see
     if its topology gives stable configurations, energies, etc.
   * Remove the use of fancy algorithms, particularly if you haven't
     equilibrated thoroughly first

#. Monitor various components of the system's energy using :ref:`gmx energy`.
   If an intramolecular term is spiking, that may indicate improper bonded
   parameters, for example.
#. Make sure you haven't been ignoring error messages (missing atoms when
   running :ref:`gmx pdb2gmx`, mismatching names when running :ref:`gmx grompp`,
   etc.) or using work-arounds (like using ``gmx grompp -maxwarn`` when you
   shouldn't be) to make sure your topology is intact and being interpreted
   correctly.
#. Make sure you are using appropriate settings in your :ref:`mdp` file for the
   force field you have chosen and the type of system you have. Particularly
   important settings are treatment of cutoffs, proper neighbor searching
   interval (``nstlist``), and temperature coupling. Improper settings can lead
   to a breakdown in the model physics, even if the starting configuration of
   the system is reasonable.

When using no explict solvent, starting your equilibration with a smaller time
step than your production run can help energy equipartition more stably.

There are several common situations in which instability frequently arises,
usually in the introduction of new species (ligands or other molecules) into
the system. To determine the source of the problem, simplify the system (e.g.
the case of a protein-ligand complex) in the following way.

#. Does the protein (in water) minimize adequately by itself? This is a test of
   the integrity of the coordinates and system preparation. If this fails,
   something probably went wrong when running :ref:`gmx pdb2gmx` (see below), or
   maybe :ref:`gmx genion` placed an ion very close to the protein (it is
   random, after all).
#. Does the ligand minimize in vacuo? This is a test of the topology. If it
   does not, check your parameterization of the ligand and any implementation of
   new parameters in force field files.
#. (If previous item is successful) Does the ligand minimize in water, and/or
   does a short simulation of the ligand in water succeed?

Other sources of possible problems are in the biomolecule topology itself.

#. Did you use ``-missing`` when running :ref:`gmx pdb2gmx`? If so, don't.
   Reconstruct missing coordinates rather than ignoring them.
#. Did you override long/short bond warnings by changing the lengths? If so,
   don't. You probably have missing atoms or some terrible input geometry.

.. _gmx-md:

Molecular dynamics
------------------

Molecular dynamics (MD) is computer simulation with atoms and/or molecules
interacting using some basic laws of physics.
The |Gromacs| :ref:`Reference Manual <md>` provides a good general introduction to this area,
as well as specific material for use with |Gromacs|. The first few chapters are mandatory reading
for anybody wishing to use |Gromacs| and not waste time.

* Introduction to molecular modeling (`slides`_, `video`_)] - theoretical framework, modeling levels,
  limitations and possibilities, systems and methods (Erik Lindahl).

Books
^^^^^

There are several text books around.

Good introductory books are:

* \A. Leach (2001) Molecular Modeling: Principles and Applications.
* \T. Schlick (2002) Molecular Modeling and Simulation

With programming background:

* \D. Rapaport (1996) The Art of Molecular Dynamics Simulation
* \D. Frenkel, B. Smith (2001) Understanding Molecular Simulation

More from the physicist's view:

* \M. Allen, D. Tildesley (1989) Computer simulation of liquids
* \H.J.C. Berendsen (2007) Simulating the Physical World: Hierarchical Modeling from Quantum Mechanics to Fluid Dynamics

Types / Ensembles
^^^^^^^^^^^^^^^^^
* NVE - number of particles (N), system volume (V) and energy (E) are constant / conserved.
* NVT - number of particles (N), system volume (V) and temperature (T) are
  constant / conserved. (See :ref:`thermostats <gmx-thermostats>` for more on *constant* temperature).
* NPT - number of particles (N), system pressure (P) and temperature (T) are constant / conserved.
  (See :ref:`pressure coupling <gmx-pressure>` for more on *constant* pressure).

.. _slides: https://extras.csc.fi/chem/courses/gmx2007/Erik_Talks/preworkshop_tutorial_introduction.pdf
.. _video:  https://video.csc.fi/playlist/dedicated/0_7z3nas0q/0_tccn9xof

.. _gmx-force-field:

Force field
-----------

Force fields are sets of potential functions and parametrized interactions that can be used to study
physical systems. A general introduction to their history, function and use is beyond the scope of this
guide, and the user is asked to consult either the relevant literature or 
try to start at the relevant `Wikipedia page`_.

.. _Wikipedia page: https://en.wikipedia.org/wiki/Force_field_(chemistry)
