Highlights
^^^^^^^^^^

|Gromacs| 2023 was released on February 6th, 2023. Patch releases may
have been made since then, please use the updated versions!  Here are
some highlights of what you can expect, along with more detail in the
links below!

As always, we've got several useful performance improvements, with or
without GPUs, all enabled and automated by default. In addition,
several new features are available for running simulations. We are extremely
interested in your feedback on how well the new release works on your
simulations and hardware. The new features are:

* The SYCL GPU implementation, which is the GPU portability layer that
  supports all major GPU platforms, has received major extensions
  in support for both platforms and features. To ensure portability
  in practice, the |Gromacs| GPU portability layer
  is actively developed with multiple SYCL implementations (hipSYCL,
  oneAPI DPC++, IntelLLVM) and regularly tested on multiple GPU backends.

  * SYCL supports more GPU offload features: bonded forces and
    direct GPU-GPU communication with GPU-aware MPI.
  * SYCL hardware support includes AMD (including RDNA support added here)
    and Intel for production as well as NVIDIA GPUs (not for production).
  * SYCL optimizations targeting important HPC platforms.

* PME decomposition has been optimized and extended to support offloading the entire
  PME calculation to multiple GPUs, including the FFT computation; when combined with
  cuFFTmp or heFFTe this enables much improved strong scaling (experimental feature).
* CUDA Graph support has been added to execute GPU-resident single-/multi-GPU
  simulations using thread-MPI entirely on the GPU to improve performance
  (experimental feature).
* Apple M1/M2 GPUs are now supported via the OpenCL GPU backend.
* New ensemble temperature mdp options allow setting the temperature of
  the ensemble for simulations without temperature coupling or with
  different reference temperatures.
* With :ref:`gmx dssp`, |Gromacs| now has a native implementation of the DSSP
  algorithm, which replaces ``gmx do_dssp``.

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!
