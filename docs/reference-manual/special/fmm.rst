.. _fmm:

Fast Multipole Method (FMM)
---------------------------

The **Fast Multipole Method (FMM)** is a numerical algorithm designed to accelerate 
the computation of long-range electrostatic interactions in particle systems. It 
achieves this by grouping distant particles into cells and approximating their collective 
influence through multipole and local expansions. This reduces the computational complexity 
from :math:`O(N^2)` to approximately :math:`O(N)`.

Within the hierarchical subdivision of the FMM, interactions are divided into long-range and 
short-range components. Cells are considered well separated when they do not share a boundary 
point, meaning that only their collective multipole contributions are used for interaction. These 
multipole-based computations are referred to as far-field interactions, evaluated using multipole 
kernels (also known as far-field kernels). Conversely, cells that do share a boundary are not well 
separated, and their particle-particle interactions must be computed explicitly. These explicit 
computations between particles in neighboring cells are known as **direct interactions**, and they 
are evaluated using direct kernels (also called near-field kernels) without any multipole approximation.

In practice, the number of neighboring cells included in the direct calculation depends on the 
specific implementation. In FMM implementations, neighborhood may include one or two layers of 
adjacent cells, depending on the chosen expansion order, desired accuracy, and the relative 
performance of direct and multipole kernels.

A detailed description of the method can be found in Ref. \ :ref:`196 <refGreengard1987>`.

Integration with |Gromacs|
~~~~~~~~~~~~~~~~~~~~~~~~~~

|Gromacs| provides an extensible **FMM module interface** that allows linking external FMM libraries. 
The interface supports flexible configurations, enabling forward and backward compatibility between 
|Gromacs| and external FMM libraries.

Specifically:

- Short-range electrostatics can be computed either by |Gromacs| kernels or by **FMM-provided kernels**.
- Long-range electrostatics can be offloaded entirely to an external FMM library.

Note that, while both configurations are supported by design, the current **|Gromacs| kernels** (PlainCPU, 
SIMD, and CUDA) compute only Lennard-Jones (LJ-only) interactions and LJ interactions are handled solely 
by |Gromacs|, not the FMM library. Support for electrostatic Coulomb interactions in |Gromacs| will be 
added in future versions.

In the current setup, FMM shares the same MPI ranks as |Gromacs|, meaning there are no dedicated FMM ranks.
Although both |Gromacs| and the FMM library may use GPU-GPU communication internally, the current 
interface requires particle coordinates to be transferred from the GPU to the CPU before being passed to 
the FMM library.

Supported Libraries
~~~~~~~~~~~~~~~~~~~

The following FMM implementations ("backends") are currently supported in the main |Gromacs| distribution:

- **ExaFMM** (`https://github.com/exafmm/ <https://github.com/exafmm/>`_)
- **FMSolvr** (`https://www.fz-juelich.de/en/jsc/projects/fmhub <https://www.fz-juelich.de/en/jsc/projects/fmhub>`_)

Usage
~~~~~

To use an external FMM library with |Gromacs|, set the mdp parameter ``fmm-backend`` to either ``exafmm`` or 
``fmsolvr``. The behaviour of the selected backend can be controlled by mdp options. 

See :ref:`this section of the documentation <mdp-fmm>` for detailed usage of these options.

Linking an FMM Library
~~~~~~~~~~~~~~~~~~~~~~

External FMM source code can be placed under :file:`src/external/` in the |Gromacs| source tree.
To integrate it, extend the ``if(GMX_USE_EXT_FMM)`` block in :file:`src/gromacs/fmm/CMakeLists.txt`,
define the correct source path, and link the built library target (for example, ``exafmm-lib`` or ``fmsolvr-lib``)
to the ``fmm`` target.

CUDA specific sources can be conditionally added within the same block using ``GMX_USE_CUDA``.

Multiple FMM backends (for example, ExaFMM or FMSolvr) can coexist under :file:`src/external/`
or in its subdirectories. Each backend provides its own implementation of :cpp:class:`FmmForceProvider`
(for example, ``FmmForceProvider::ExafmmImpl`` or ``FmmForceProvider::FmSolvrImpl``). While multiple backends 
can be compiled simultaneously, only one backend is active at runtime.

An example CMake snippet is provided below:

.. code-block:: cmake

    if(GMX_USE_EXT_FMM)
      # Specify the external FMM source path
      set(FMM_DIR ${CMAKE_SOURCE_DIR}/src/external/<fmm_backend>)

      # Collect all source files
      file(GLOB_RECURSE FMM_SOURCES ${FMM_DIR}/src/*.cpp)
      list(APPEND FMM_INCLUDE_DIRS ${FMM_DIR}/include)

      # Optionally include GPU sources
      if(GMX_GPU_CUDA)
        file(GLOB_RECURSE FMM_CUDA_SOURCES ${FMM_DIR}/src/*.cu)
        list(APPEND FMM_SOURCES ${FMM_CUDA_SOURCES})
        list(APPEND FMM_INCLUDE_DIRS ${FMM_DIR}/include/cuda)
      endif()

      # Define the external FMM library target
      add_library(fmm-lib STATIC ${FMM_SOURCES})
      target_include_directories(fmm-lib PUBLIC ${FMM_INCLUDE_DIRS})

      # Link it with the GROMACS FMM target
      target_link_libraries(fmm PRIVATE fmm-lib)
    endif()
