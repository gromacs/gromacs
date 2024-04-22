.. _install guide exotic:

********************************************
Installation guide for exotic configurations
********************************************

.. highlight:: bash

Special instructions for building |Gromacs| on less-common systems
------------------------------------------------------------------

These instructions pertain to building |Gromacs| |version|.
This document is complementary to the `up-to-date installation instructions`_ instructions.

The configurations listed here are expected to work, but **are not recommended for typical users**.

.. _install guide exotic sycl:

SYCL GPU acceleration for AMD and NVIDIA GPUs using Intel oneAPI DPC++
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

AMD and NVIDIA GPUs can also be used with Intel oneAPI BaseKit and Codeplay oneAPI plugins.

For most users, we recommend using :ref:`CUDA <CUDA GPU acceleration>` for NVIDIA GPUs and
:ref:`AdaptiveCpp <SYCL GPU acceleration AMD>` for AMD GPUs instead.

With some versions of oneAPI, you might receive "The compiler you are using does not support OpenMP parallelism"
error from CMake. In this case, please add the following options to your CMake command:

- For oneAPI 2024.x: ``-DCMAKE_C_FLAGS="-isystem /opt/intel/oneapi/compiler/latest/opt/compiler/include"  -DCMAKE_CXX_FLAGS="-isystem /opt/intel/oneapi/compiler/latest/opt/compiler/include"``
- For oneAPI 2023.x: ``-DCMAKE_C_FLAGS="-isystem /opt/intel/oneapi/compiler/latest/linux/compiler/include"  -DCMAKE_CXX_FLAGS="-isystem /opt/intel/oneapi/compiler/latest/linux/compiler/include"``

AMD GPUs
""""""""

After installing Intel oneAPI toolkit 2023.0 or newer, a compatible ROCm version,
and the `Codeplay plugin <https://developer.codeplay.com/products/oneapi/amd/home/>`_,
set up the environment by running ``source /opt/intel/oneapi/setvars.sh --include-intel-llvm``
or loading an appropriate :command:`module load` on an HPC system.

Then, configure |Gromacs| using the following command (replace ``gfxXYZ`` with the target architecture):

::

   cmake .. -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
            -DGMX_GPU=SYCL -DGMX_SYCL=DPCPP \
            -DGMX_GPU_NB_CLUSTER_SIZE=8 -DGMX_GPU_FFT_LIBRARY=vkfft \
            -DSYCL_CXX_FLAGS_EXTRA='-fsycl-targets=amdgcn-amd-amdhsa;-Xsycl-target-backend;--offload-arch=gfxXYZ'


NVIDIA GPUs
"""""""""""

After installing Intel oneAPI toolkit 2023.0 or newer, a compatible CUDA version,
and the `Codeplay plugin <https://developer.codeplay.com/products/oneapi/nvidia/home/>`__,
set up the environment by running ``source /opt/intel/oneapi/setvars.sh --include-intel-llvm``
or loading an appropriate :command:`module load` on an HPC system.

 Then, configure |Gromacs| using the following command:

::

   cmake .. -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
            -DGMX_GPU=SYCL -DGMX_SYCL=DPCPP \
            -DGMX_GPU_NB_CLUSTER_SIZE=8 -DGMX_GPU_FFT_LIBRARY=vkfft \
            -DSYCL_CXX_FLAGS_EXTRA=-fsycl-targets=nvptx64-nvidia-cuda

.. _install guide exotic adaptivecpp:

SYCL GPU acceleration for NVIDIA GPUs using AdaptiveCpp (hipSYCL)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


For most users, we recommend using :ref:`CUDA <CUDA GPU acceleration>` for NVIDIA GPUs.

Build and install AdaptiveCpp_ with CUDA backend (we recommend using the mainline Clang, not the ROCm-bundled one).

Then, use the following command to build |Gromacs| (make sure to use the same compiler and set target GPU architecture
instead of ``sm_XY``):

::

   cmake .. -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ \
            -DGMX_GPU=SYCL -DGMX_SYCL=ACPP -DHIPSYCL_TARGETS='cuda:sm_XY'

.. _install guide static linking:

Static linking
~~~~~~~~~~~~~~

Dynamic linking of the |Gromacs| executables will lead to a
smaller disk footprint when installed, and so is the default on
platforms where we believe it has been tested repeatedly and found to work.
In general, this includes Linux, Windows, Mac OS X and BSD systems.
Static binaries take more space, but on some hardware and/or under
some conditions they are recommended or even necessary, most commonly when you are running large parallel
simulation using MPI libraries (e.g. Cray).

* To link |Gromacs| binaries statically against the internal |Gromacs|
  libraries, set ``-DBUILD_SHARED_LIBS=OFF``.
* To link statically against external (non-system) libraries as well,
  set ``-DGMX_PREFER_STATIC_LIBS=ON``. Note, that in
  general ``cmake`` picks up whatever is available, so this option only
  instructs ``cmake`` to prefer static libraries when both static and
  shared are available. If no static version of an external library is
  available, even when the aforementioned option is ``ON``, the shared
  library will be used. Also note that the resulting binaries will
  still be dynamically linked against system libraries on platforms
  where that is the default. To use static system libraries,
  additional compiler/linker flags are necessary, e.g. ``-static-libgcc
  -static-libstdc++``.
* To attempt to link a fully static binary set
  ``-DGMX_BUILD_SHARED_EXE=OFF``. This will prevent CMake from explicitly
  setting any dynamic linking flags. This option also sets
  ``-DBUILD_SHARED_LIBS=OFF`` and ``-DGMX_PREFER_STATIC_LIBS=ON`` by
  default, but the above caveats apply. For compilers which don't
  default to static linking, the required flags have to be specified. On
  Linux, this is usually ``CFLAGS=-static CXXFLAGS=-static``.


Building on Solaris
~~~~~~~~~~~~~~~~~~~

The built-in |Gromacs| processor detection does not work on Solaris,
so it is strongly recommended that you build |Gromacs| with
``-DGMX_HWLOC=on`` and ensure that the ``CMAKE_PREFIX_PATH`` includes
the path where the hwloc headers and libraries can be found. At least
version 1.11.8 of hwloc is recommended.
