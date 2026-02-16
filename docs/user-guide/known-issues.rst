Known issues affecting users of |Gromacs|
=========================================

.. _gmx-users-known-issues:

Here is a non-exhaustive list of issues that are we are aware of that are
affecting regular users of |Gromacs|.

"Cannot find a working standard library" error with ROCm Clang
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Some Clang installations don't contain a compatible C++ standard library.
In such cases, you might have to install ``g++`` and help CMake find it
by setting ``-DGMX_GPLUSGPLUS_PATH=/path/to/bin/g++``.

On Ubuntu 22.04, installing GCC 12 standard library (with
``sudo apt install libstdc++-12-dev``) usually works well even without
setting ``-DGMX_GPLUSGPLUS_PATH``.

:issue:`4679`

Expanded ensemble does not checkpoint correctly
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the legacy simulator, because of shortcomings in the
implementation, successful expanded-ensemble MC steps that occurred on
checkpoint steps were not recorded in the checkpoint. If that
checkpoint was used for a restart, then it would not necessarily
behave correctly and reproducibly afterwards. So checkpointing of
expanded-ensemble simulations is disabled for the legacy simulator.

Checkpointing of expanded ensemble in the modular simulator works
correctly.

To work around the issue, either avoid ``-update gpu`` (so that it
uses the modular simulator path which does not have
the bug), or use an older version of |Gromacs|
(which does do the buggy checkpointing), or refrain from
restarting from checkpoints in the affected case.

:issue:`4629`

Compiling with GCC 12-14 on POWER9 architectures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There are multiple failing unit tests after compilation with some versions of
GCC 12-14 on POWER9 architectures.
The problem is fixed in GCC 12.5, GCC 13.4, and GCC 14.2, respectively.
GCC 11 and earlier, as well as GCC 15, are *not* affected.

:issue:`4823`

NbnxmTest crash with oneAPI 2024.1
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When building with oneAPI 2024.1, the ``NbnxmTest`` test can segfault in
some cases. Using oneAPI 2024.2 or newer should resolve the issue.

:issue:`5247`

Severe performance regression with SVE and LLVM 20
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is a major performance regression on AArch64 CPUs
when using SVE SIMD and building with LLVM 20.

The problem can be worked around by using LLVM 19
or setting ``-DGMX_SIMD=ARM_NEON_ASIMD``.

The issue is expected to be fixed in LLVM 21.

:issue:`5390`

