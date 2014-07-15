# Non-bonded cut-off schemes

The default cut-off scheme in GROMACS 5.0 is based on classical
buffered Verlet lists. These are implemented extremely efficiently
on modern CPUs and accelerators, and support nearly all of the
algorithms used in GROMACS.

Before version 4.6, GROMACS always used pair-lists based on groups of
atoms. These groups of atoms were orginally charge-groups, which were
necessary with plain cut-off electrostatics. With the use of PME (or
reaction-field with a buffer), charge groups are no longer necessary
(and are ignored in the Verlet scheme). In GROMACS 4.6 and later, the
group-based cut-off scheme is still available, but is **deprecated in
5.0**. It is still available mainly for backwards compatibility, to
support the algorithms that have not yet been converted, and for the
few cases where it may allow faster simulations with bio-molecular
systems dominated by water.

The group cut-off scheme should generally be combined with a buffered
pair-list to help avoid artefacts. If a buffer is not used, then
energy drift can result from either atoms moving across the cut-off in
the period between neighbour-search steps, or, in the case of
multi-atom charge groups, from atoms moving across the cut-off
distance when the center of geometry has not. However, the
group-scheme kernels that can implement this are much slower than
either the unbuffered group-scheme kernels, or the buffered
Verlet-scheme kernels. Use of the Verlet scheme is strongly
encouraged, because it is easier and faster to do a correct
simulation. In particular, GPU acceleration is available only with the
Verlet scheme.

The Verlet scheme uses properly buffered lists with exact cut-offs.
Both the LJ and Coulomb potential are shifted to zero by subtracting
the value at the cut-off. This ensures that the energy is the integral
of the force. Still it is advisable to have small forces at the
cut-off, hence to use PME or reaction-field with infinite epsilon.

In the Verlet scheme, particle-pair forces (and energies when
necessary) are calculated in groups of `M x N` particles, where `M`
and `N` are typically 2, 4 or 8. This is convenient for modern
streaming processors but leads to the computation of a significant
number of interactions of zero strength, because their particles lie
outside the cut-off. This is still more efficient than the buffered
group scheme, however.

## Non-bonded scheme feature comparison

All GROMACS features not directly related to non-bonded interactions
are supported in both schemes. Eventually, most non-bonded features
will be supported in the Verlet scheme. A table describing the
compatibility of just non-bonded features with the two schemes is
given below.

Table: Support levels within the group and Verlet cut-off schemes
for features related to non-bonded interactions

Feature                               group        Verlet
------------------------------        -----        -------
unbuffered cut-off scheme             default      not by default
exact cut-off                         shift/switch always
shifted interactions                  force+energy energy
switched potential                    yes          yes
switched forces                       yes          yes
non-periodic systems                  yes          Z + walls
implicit solvent                      yes          no
free energy perturbed non-bondeds     yes          yes
energy group contributions            yes          only on CPU
energy group exclusions               yes          no
AdResS multi-scale                    yes          no
OpenMP multi-threading                only PME     all
native GPU support                    no           yes
Coulomb PME                           yes          yes
Lennard-Jones PME                     yes          yes
virtual sites                         yes          yes
User-supplied tabulated interactions  yes          no
Buckingham VdW interactions           yes (slow)   no
rcoulomb != rvdw                      yes          no
twin-range                            yes          no

## Performance

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
the cut-off, so you could win a bit of peformance by reducing or
removing the buffer, but this might not be a good trade-off of
simulation quality.

The table below shows a performance comparison of most of the relevant
setups. Any atomistic model will have performance comparable to tips3p
(which has LJ on the hydrogens), unless a united-atom force field is
used. The performance of a protein in water will be between the tip3p
and tips3p performance. The group scheme is optimized for water
interactions, which means a single charge group containing one atom
with LJ, and 2 or 3 atoms without LJ. Such kernels for water are
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
non-bonded setups in GROMACS using either 8 thread-MPI ranks (group
scheme), or 8 OpenMP threads (Verlet scheme). 3000 atoms, 1.0 nm
cut-off, PME with 0.11 nm grid, dt=2 fs, Intel Core i7 2600, 3.4 GHz +
Nvidia GTX660Ti

------            ----------  --------  --------  -------------
system            group,      group,    Verlet,   Verlet
                  unbuffered  buffered  buffered  buffered, GPU
------            ----------  --------  --------  -------------
tip3p,            208         116       170       450
charge groups

tips3p,           129         63        162       450
charge groups

tips3p,           104         75        162       450
no charge groups
------            ----------  --------  --------  -------------

## How to use the Verlet scheme
You can use the Verlet cut-off scheme simply by setting in your mdp file:
    cutoff-scheme           = Verlet
    verlet-buffer-tolerance = 0.005

The value of [.mdp] option [`verlet-buffer-tolerance`] will add a
pair-list buffer whose size is tuned for the given energy drift (in
kJ/mol/ns per atom). The effective drift is usually much lower, as
[grompp] assumes constant particle velocities. (Note that in single
precision for normal atomistic simulations constraints cause a drift
somewhere around 0.0001 kJ/mol/ns per atom, so it doesn't make sense
to go much lower.) Details on how the buffer size is chosen can be
found in the reference below and in the Reference Manual.

For constant-energy (NVE) simulations, the buffer size will be
inferred from the temperature that corresponds to the velocities
(either those generated, if applicable, or those found in the input
configuration). Alternatively, [`verlet-buffer-tolerance`] can be set
to -1 and a buffer set manually by specifying [`rlist`] greater than
the larger of [`rcoulomb`] and [`rvdw`]. The simplest way to get a
reasonable buffer size is to use an NVT mdp file with the target
temperature set to what you expect in your NVE simulation, and
transfer the buffer size printed by grompp to your NVE [.mdp] file.

When a GPU is used, nstlist is automatically increased by mdrun,
usually to 20 or more; rlist is increased along to stay below the
target energy drift. Further information on [running mdrun with
GPUs] is available.

## Further information

For further information on algorithmic and implementation details of
the Verlet cut-off scheme and the MxN kernels, as well as detailed
performance analysis, please consult the following article:

Páll, S. and Hess, B. A flexible algorithm for calculating pair
interactions on SIMD architectures. Comput. Phys. Commun. 184,
2641–2650 (2013). <http://dx.doi.org/10.1016/j.cpc.2013.06.003>
