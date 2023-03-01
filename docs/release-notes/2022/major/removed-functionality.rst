Removed functionality
^^^^^^^^^^^^^^^^^^^^^

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Removed mdrun-only build configuration
""""""""""""""""""""""""""""""""""""""

The need for the mdrun-only build of |Gromacs| has expired, as it has
the same set of dependencies as regular |Gromacs|. It was deprecated
in |Gromacs| 2021. Removing it will simplify maintenance, testing,
documentation, installation, and teaching new users.

:issue:`3808`

Removed support for x86 MIC, ARMv7, Sparc64 HPC-ACE, and IBM VMX SIMD
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These platforms are dead in HPC and so are no longer supported. The
KNL platform is unaffected by this change.

:issue:`3891`

Removed deprecated environment variables
""""""""""""""""""""""""""""""""""""""""

The following environment variables were removed after being deprecated
in favor of better-named alternatives:

* ``GMX_CUDA_NB_ANA_EWALD`` and ``GMX_OCL_NB_ANA_EWALD`` (use ``GMX_GPU_NB_ANA_EWALD``)
* ``GMX_CUDA_NB_TAB_EWALD`` and ``GMX_OCL_NB_TAB_EWALD`` (use ``GMX_GPU_NB_TAB_EWALD``)
* ``GMX_CUDA_NB_EWALD_TWINCUT`` and ``GMX_OCL_NB_EWALD_TWINCUT`` (use ``GMX_GPU_NB_EWALD_TWINCUT``)

:issue:`3803`

Removed the ability for gmx wham to read .pdo files
"""""""""""""""""""""""""""""""""""""""""""""""""""

Files in .pdo format were written by |Gromacs| versions prior to 4.0.
That is so long ago that being able to read them is no longer
relevant, so this capability was deprecated in version 2021. If you do
need to read such files, please use an older version of |Gromacs|.

Removed 32bit support
"""""""""""""""""""""

We deprecated 32bit support in 2020 and have had no way to test it ourselves for a while
before that. Those architectures are no longer relevant in HPC, so we officially no longer
support building |Gromacs| on them.
