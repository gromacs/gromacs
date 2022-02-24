Highlights
^^^^^^^^^^

|Gromacs| 2022 was released on February 22st, 2022. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. In addition,
several new features are available for running simulations. We are extremely
interested in your feedback on how well the new release works on your
simulations and hardware. The new features are:

* Free-energy kernels are accelerated using SIMD, which make free-energy
  calculations up to three times as fast when using GPUs
* A new formulation of the soft-cored non-bonded interactions for free-energy calculations allows for a finer control of the alchemical transformation pathways
* New transformation pull coordinate allows arbitrary mathematical transformations of one of more other pull coordinates
* New interface for multi-scale Quantum Mechanics / Molecular Mechanics (QM/MM) simulations with the CP2K quantum 
  chemistry package, supporting periodic boundary conditions.
* grompp performance improvements
* `Cool quotes music playlist <https://open.spotify.com/playlist/4oj41X9tgIAJuLgfWPq6ZX>`_
* Additional features were ported to modular simulator
* Added AMD GPU support with SYCL via hipSYCL_
* More GPU offload features supported with SYCL (PME, GPU update). 
* Improved parallelization with GPU-accelerated runs using CUDA and extended GPU direct communication to support multi-node simulation using CUDA-aware MPI.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!
