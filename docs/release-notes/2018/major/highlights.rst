Highlights
^^^^^^^^^^

|Gromacs| 2018 was released on January 10, 2018. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, and all enabled and automated by default. We are
extremely interested in your feedback on how well this worked on your
simulations and hardware. They are:

* PME long-ranged interactions can now run on a single GPU, which
  means many fewer CPU cores are needed for good performance.
* Optimized SIMD support for recent CPU architectures:
  AMD Zen, Intel Skylake-X and Skylake Xeon-SP.

There are some new features available also:

* The AWH (Accelerated Weight Histogram) method is now supported,
  which is an adaptive biasing method used for overcoming free energy
  barriers and calculating free energies (see
  https://doi.org/10.1063/1.4890371).
* A new dual-list dynamic-pruning algorithm for the short-ranged
  interactions, that uses an inner and outer list to permit a longer-lived
  outer list, while doing less work overall and making runs
  less sensitive to the choice of the "nslist" parameter.
* A physical validation suite is added, which runs a series of short
  simulations, to verify the expected statistical properties,
  e.g. of energy distributions between the simulations, as a sensitive
  test that the code correctly samples the expected ensemble.
* Conserved quantities are computed and reported for more integration
  schemes - now including all Berendsen and Parrinello-Rahman schemes.
