Highlights
^^^^^^^^^^

|Gromacs| 2016 was released on August 4, 2016. Patch releases
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

* As always, we've got several useful performance improvements, with or
  without GPUs. CPU-side SIMD and threading enhancements will
  make GPU-accelerated simulations faster even if we'd left the GPU
  code alone! Thanks to these and additional GPU kernel improvements,
  in GPU-accelerated runs expect around 15% improvement
  in throughput. (And not just for plain vanilla MD, either... the
  pull code now supports OpenMP threading throughout, and
  multi-simulations have less coupling between simulations.)
* We have a new C++11 portability layer permitting us to accelerate in
  SIMD on the CPU lots of minor routines. These will also often
  improve runs that use accelerators or many nodes through better load
  balancing. POWER8, ARM64, AVX512 (KNL), and more are fully SIMD accelerated now
  because they are supported in the new portability layer!
* We made further SIMD acceleration of bonded interactions which
  reduces their calculation time by about a factor of 2. This improves
  load balance at high parallelization by a factor of 2, and shows
  significantly better scaling.
* Similarly, SIMD acceleration of SETTLE reduces the time for
  constraints by a factor of 3 to 5 - which has a strong effect for GPU runs.
* OpenCL GPU support is now available with all combinations of MPI,
  thread-MPI and GPU sharing (ie. the same as CUDA). Kernel performance
  has improved by up to 60%. AMD GPUs benefit the most, OpenCL on NVIDIA is
  generally still slow.
* Tools in the new analysis framework can handle trajectories that
  are subsets of the simulation system.
* New pull coordinate geometries angle-axis, dihedral, and normal angle.
* Checkpoint restarts work only in the cases where the implementation
  can always do what the user wants.
* The version numbering has changed to be the year of the release,
  plus (in future) a patch number. |Gromacs| 2016 will be the initial
  release from this branch, then |Gromacs| 2016.1 will have the set of
  bugs that have been fixed in |Gromacs| 2016, etc.
