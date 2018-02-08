Performance improvements
^^^^^^^^^^^^^^^^^^^^^^^^

Implemented support for PME long-ranged interactions on GPUs
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
A single GPU can now be used to accelerate the computation of the
long-ranged PME interactions. This feature provides excellent
performance improvements, in particular that only 2-4 CPU cores per
GPU will be about as fast as the 2016 version that needed many more
CPU cores to balance the GPU. Performance on hardware that had good
balance of GPU and CPU also shows minor improvements, and the capacity
for hardware with strong GPUs to run effective simulations is now
greatly improved.

Currently, the GPU used for PME must be either the same GPU as used
for the short-ranged interactions and in the same single rank of the
simulation, or any GPU used from a PME-only rank. mdrun -pme gpu now
requires that PME runs on a GPU, if supported. All CUDA versions and
hardware generations supported by |Gromacs| can run this code path,
including CUDA 9.0 and Volta GPUs. However, not all combinations
of features are supported with PME on GPUs - notably FEP calculations
are not yet available.

The user guide is updated to reflect the new capabilities, and more
documentation will be forthcoming.

Added more SIMD intrinsics support for PME spread and gather
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Achieved speedup on Intel KNL processors of around 11% for PME
spread/gather on typical simulation systems.

Added SIMD intrinsics version of simple update
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
In the simple case of leap-frog without pressure coupling and with at
most one temperature-coupling group, the update of velocities and
coordinates is now implemented with SIMD intrinsics for improved
simulation rate.

Add SIMD intrinsics version of Urey-Bradley angle kernel
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
For steps where energies and shift forces are not required, this kernel
improves performance, which can otherwise be rate limiting in GPU-accelerated
runs, particularly with CHARMM force fields.

Use OpenMP up to 16 threads with AMD Ryzen when automating run setup
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
AMD Ryzen appears to always perform slightly better with OpenMP
than MPI, up to using all 16 threads on the 8-core die.

128-bit AVX2 SIMD for AMD Ryzen
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
While Ryzen supports 256-bit AVX2, the internal units are organized
to execute either a single 256-bit instruction or two 128-bit SIMD
instruction per cycle. Since most of our kernels are slightly
less efficient for wider SIMD, this improves performance by roughly
10%.

Choose faster nbnxn SIMD kernels on AMD Zen
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
On AMD Zen, tabulated Ewald kernels are always faster than analytical.
And with AVX2_256 2xNN kernels are faster than 4xN.
These faster choices are now made based on CpuInfo at run time.

Refs :issue:`2328`

Enabled group-scheme SIMD with GMX_SIMD=AVX2_128
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The group-scheme kernels can use AVX instructions from either the
AVX_128_FMA and AVX_256 extensions. But hardware that supports the new
AVX2_128 extensions also supports AVX_256, so we enable such support
for the group-scheme kernels.

Detect AVX-512 FMA units to choose best SIMD
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Recent Intel x86 hardware can have multiple AVX-512 FMA units, and the
number of those units and the way their use interacts with the way the
CPU chooses its clock speed mean that it can be advantageous to avoid
using AVX-512 SIMD support in |Gromacs| if there is only one such
unit.  Because there is no way to query the hardware to count the
number of such units, we run code at CMake and mdrun time to compare
the performance from using such units, and recommend the version that
is best. This may mean that building |Gromacs| on the front-end node
of the cluster might not suit the compute nodes, even when they are
all from the same generation of Intel's hardware.

Speed up nbnxn buffer clearing
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Tweaked conditional in the nonbonded GPU kernels
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
GPU compilers miss an easy optimization of a loop invariant in the
inner-lop conditional. Precomputing part of the conditional together
with using bitwise instead of logical and/or improves performance with
most compilers by up to 5%.

