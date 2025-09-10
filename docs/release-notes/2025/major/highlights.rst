Highlights
^^^^^^^^^^

|Gromacs| 2025.0 was released on February 11th, 2025. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. In addition,
several new features are available for running simulations. We are extremely
interested in your feedback on how well the new release works on your
simulations and hardware. The new features are:

* A feature-limited version of `PLUMED <https://www.plumed.org/>`_ can be used,
  on a non-Windows build of |Gromacs|, without needing to apply a patch.

* Basic support for running simulations with Neural Network Potential (NNP) models
  trained in `PyTorch <https://pytorch.org/>`_.

* Extended OpenMP parallelization of pair search and domain decomposition can improve the performance, especially relevant with fast GPUs.

* Added support for using AMD HIP as GPU backend. This is currently limited to the
  NBNxM (nonbonded interactions within cut-off) kernels.

* Support for continuing expanded ensemble equilibration across simulations by
  enabling initialization of `init-lambda-counts` and `init-wl-histogram-counts`
  through mdp options.

* GPU-direct communication is now used by default when the MPI library 
  supports it.

* Enhanced PP halo exchange using GPU kernel-initiated communication implemented using NVSHMEM, improving
  performance when scaling to multiple NVIDIA GPUs.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!
