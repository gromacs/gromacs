Non-bonded cut-off schemes
==========================

The default cut-off scheme in |Gromacs| |version| is based on classical
buffered Verlet lists. These are implemented extremely efficiently
on modern CPUs and accelerators, and support nearly all of the
algorithms used in |Gromacs|.

Before version 4.6, |Gromacs| always used pair-lists based on groups of
particles. These groups of particles were originally charge-groups, which were
necessary with plain cut-off electrostatics. With the use of PME (or
reaction-field with a buffer), charge groups are no longer necessary
(and are ignored in the Verlet scheme). In |Gromacs| 4.6 and later, the
group-based cut-off scheme is still available, but is **deprecated since
5.0**. It is still available mainly for backwards
compatibility, to support the algorithms that have not yet been
converted, and for the few cases where it may allow faster simulations
with bio-molecular systems dominated by water.

Without PME, the group cut-off scheme should generally be combined
with a buffered pair-list to help avoid artifacts. However, the
group-scheme kernels that can implement this are much slower than
either the unbuffered group-scheme kernels, or the buffered
Verlet-scheme kernels. Use of the Verlet scheme is strongly encouraged
for all kinds of simulations, because it is easier and faster to run
correctly. In particular, GPU acceleration is available only with the
Verlet scheme.

The Verlet scheme uses properly buffered lists with exact cut-offs.
The size of the buffer is chosen with :mdp:`verlet-buffer-tolerance`
to permit a certain level of drift.  Both the LJ and Coulomb potential
are shifted to zero by subtracting the value at the cut-off. This
ensures that the energy is the integral of the force. Still it is
advisable to have small forces at the cut-off, hence to use PME or
reaction-field with infinite epsilon.

Non-bonded scheme feature comparison
------------------------------------

All |Gromacs| |version| features not directly related to non-bonded
interactions are supported in both schemes. Eventually, all non-bonded
features will be supported in the Verlet scheme. A table describing
the compatibility of just non-bonded features with the two schemes is
given below.

Table: Support levels within the group and Verlet cut-off schemes
for features related to non-bonded interactions

====================================  ============ =======
Feature                               group        Verlet
====================================  ============ =======
unbuffered cut-off scheme             default      not by default
exact cut-off                         shift/switch always
potential-shift interactions          yes          yes
potential-switch interactions         yes          yes
force-switch interactions             yes          yes
switched potential                    yes          yes
switched forces                       yes          yes
non-periodic systems                  yes          Z + walls
implicit solvent                      yes          no
free energy perturbed non-bondeds     yes          yes
energy group contributions            yes          only on CPU
energy group exclusions               yes          no
OpenMP multi-threading                only PME     all
native GPU support                    no           yes
Coulomb PME                           yes          yes
Lennard-Jones PME                     yes          yes
virtual sites                         yes          yes
User-supplied tabulated interactions  yes          no
Buckingham VdW interactions           yes          no
rcoulomb != rvdw                      yes          yes
twin-range                            no           no
====================================  ============ =======

Performance
-----------

The performance of the group cut-off scheme depends very much on the
composition of the system and the use of buffering. There are
optimized kernels for interactions with water, so anything with a lot
of water runs very fast. But if you want properly buffered
interactions, you need to add a buffer that takes into account both
charge-group size and diffusion, and check each interaction against
the cut-off length each time step. This makes simulations much
slower. The performance of the Verlet scheme with the new non-bonded
kernels is independent of system composition and is intended to always
run with a buffered pair-list. Typically, buffer size is 0 to 10% of
the cut-off, so you could win a bit of performance by reducing or
removing the buffer, but this might not be a good trade-off of
simulation quality.

The table below shows a performance comparison of most of the relevant
setups. Any atomistic model will have performance comparable to tips3p
(which has LJ on the hydrogens), unless a united-atom force field is
used. The performance of a protein in water will be between the tip3p
and tips3p performance. The group scheme is optimized for water
interactions, which means a single charge group containing one particle
with LJ, and 2 or 3 particles without LJ. Such kernels for water are
roughly twice as fast as a comparable system with LJ and/or without
charge groups. The implementation of the Verlet cut-off scheme has no
interaction-specific optimizations, except for only calculating half
of the LJ interactions if less than half of the particles have LJ. For
molecules solvated in water the scaling of the Verlet scheme to higher
numbers of cores is better than that of the group scheme, because the
load is more balanced. On the most recent Intel CPUs, the absolute
performance of the Verlet scheme exceeds that of the group scheme,
even for water-only systems.

Table: Performance in ns/day of various water systems under different
non-bonded setups in |Gromacs| using either 8 thread-MPI ranks (group
scheme), or 8 OpenMP threads (Verlet scheme). 3000 particles, 1.0 nm
cut-off, PME with 0.11 nm grid, dt=2 fs, Intel Core i7 2600 (AVX), 3.4
GHz + Nvidia GTX660Ti

========================  =================  ===============  ================  =====================
system                    group, unbuffered  group, buffered  Verlet, buffered  Verlet, buffered, GPU
========================  =================  ===============  ================  =====================
tip3p, charge groups      208                116              170               450
tips3p, charge groups     129                63               162               450
tips3p, no charge groups  104                75               162               450
========================  =================  ===============  ================  =====================

How to use the Verlet scheme
----------------------------

The Verlet scheme is enabled by default with option :mdp:`cutoff-scheme`.
The value of [.mdp] option :mdp:`verlet-buffer-tolerance` will add a
pair-list buffer whose size is tuned for the given energy drift (in
kJ/mol/ns per particle). The effective drift is usually much lower, as
:ref:`gmx grompp` assumes constant particle velocities. (Note that in single
precision for normal atomistic simulations constraints cause a drift
somewhere around 0.0001 kJ/mol/ns per particle, so it doesn't make sense
to go much lower.) Details on how the buffer size is chosen can be
found in the reference below and in the `reference manual`_.

.. _reference manual: gmx-manual-parent-dir_

For constant-energy (NVE) simulations, the buffer size will be
inferred from the temperature that corresponds to the velocities
(either those generated, if applicable, or those found in the input
configuration). Alternatively, :mdp:`verlet-buffer-tolerance` can be set
to -1 and a buffer set manually by specifying :mdp:`rlist` greater than
the larger of :mdp:`rcoulomb` and :mdp:`rvdw`. The simplest way to get a
reasonable buffer size is to use an NVT mdp file with the target
temperature set to what you expect in your NVE simulation, and
transfer the buffer size printed by :ref:`gmx grompp` to your NVE [.mdp] file.

When a GPU is used, nstlist is automatically increased by :ref:`gmx mdrun`,
usually to 20 or more; rlist is increased along to stay below the
target energy drift. Further information on running :ref:`gmx mdrun` with
GPUs :ref:`is available<gmx-mdrun-on-gpu>`.

Further information
-------------------

For further information on algorithmic and implementation details of
the Verlet cut-off scheme and the MxN kernels, as well as detailed
performance analysis, please consult the following article:

`Páll, S. and Hess, B. A flexible algorithm for calculating pair
interactions on SIMD architectures. Comput. Phys. Commun. 184,
2641–2650 (2013). <http://dx.doi.org/10.1016/j.cpc.2013.06.003>`__
