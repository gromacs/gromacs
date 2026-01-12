Deprecated functionality
------------------------

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without
   a space between the colon and number!

Changes anticipated to |Gromacs| 2026 functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functionality deprecated in |Gromacs| 2026
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generating angle constraints with grompp is deprecated
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Google Native Client (NaCl) build deprecated
""""""""""""""""""""""""""""""""""""""""""""

Last release of the sandbox has been 10 years ago
and the support for NaCl has been removed from LLVM 22.
Building option ``GMX_NACL`` is now deprecated in |Gromacs|.


Functionality deprecated in |Gromacs| 2025
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MTTK pressure coupling is deprecated
""""""""""""""""""""""""""""""""""""

:issue:`5072`

The TNG trajectory format is deprecated
"""""""""""""""""""""""""""""""""""""""

The TNG file format will be removed in a future release. It will be replaced by the more
widely used HDF5-based H5MD format. There will be at least one version of |Gromacs|
supporting both H5MD and TNG to allow conversions between the two formats.

:issue:`5225`

Functionality deprecated in |Gromacs| 2024
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The analysis tool ``gmx gyrate-legacy`` deprecated
""""""""""""""""""""""""""""""""""""""""""""""""""

``gmx gyrate`` has been partly re-implemented in the modern |Gromacs| analysis framework.
The old implementation is still available as ``gmx gyrate-legacy``. Please plan to update to the new version.
Please let the |Gromacs| developers know of desired functionality missing from, or broken in, the new implementation.

:issue:`4927`

The analysis tool ``gmx hbond-legacy`` deprecated
"""""""""""""""""""""""""""""""""""""""""""""""""

``gmx hbond`` has been partly re-implemented in the modern |Gromacs| analysis framework.
The old implementation is still available as ``gmx hbond-legacy``. Please plan to update to the new version.
Please let the |Gromacs| developers know of desired functionality missing from, or broken in, the new implementation.

:issue:`4927`

The analysis tools ``gmx sans`` and ``gmx saxs`` deprecated
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

``gmx sans`` and ``gmx saxs`` has been partly re-implemented under new name ``gmx scattering`` in the modern |Gromacs| analysis framework.
The old implementations are still available as ``gmx sans-legacy`` and ``gmx saxs-legacy``. Please plan to update to the new version.
Please let the |Gromacs| developers know of desired functionality missing from, or broken in, the new implementation.

:issue:`4927`

Functionality deprecated in |Gromacs| 2022
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GMX_OPENCL_NB_CLUSTER_SIZE CMake variable deprecated in favor of GMX_GPU_NB_CLUSTER_SIZE
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Both OpenCL and SYCL support different cluster sizes, so GMX_GPU_NB_CLUSTER_SIZE should
be used going forward.

Guessing masses and atomic radii from atom names is deprecated
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

When atom masses or van-der-Waals radii are needed, we suggest building
a proper |Gromacs| topology instead of using PDB files directly, even
if the tool supports it.

:issue:`3368`
:issue:`4288`

Functionality deprecated in |Gromacs| 2021
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``mdrun -deffnm`` to be removed
"""""""""""""""""""""""""""""""

This functionality is convenient when running very simple simulations,
because it permits grouping of a set of files that then differ only
their suffix. However, it does not work in the wider case of an
``mdrun`` module (or modules) writing multiple ``.xvg`` output
files. The resulting filenames collide. That, and its interaction with
checkpointing and appending, have led to quite a few bug reports.

Because users can use a folder to group files (a standard mechanism
that they understand from experience outside of |Gromacs|), we can
build and test better software for them if we remove the erstwhile
convenience of ``mdrun -deffnm``. Please update your workflows
accordingly.

:issue:`3818`

OpenCL to be removed as a GPU framework
"""""""""""""""""""""""""""""""""""""""

Since there are no prospects of an emerging GPU vendor in HPC needing
OpenCL support, we intend to focus our efforts on CUDA, SYCL, and HIP
for GPU acceleration in |Gromacs|. OpenCL support has been deprecated
for a while, and will be removed in a future release.

:issue:`3818`

Support for version 1 of the hardware locality library ``hwloc``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Version 2 has been supported in |Gromacs| for several
years. The capabilities of newer hardware and hardware-support APIs
are of most interest for |Gromacs| moving forward, so we should
minimize our testing work and encourage clusters to upgrade older
``hwloc`` installations.

:issue:`3818`

Legacy API
""""""""""

The legacy installed headers have been deprecated for a
while, however we wish to state more broadly that all headers found
within the ``src`` directory tree of |Gromacs| are intended for
internal consumption only, and are thus subject to change without
notice. Further, the form and contents of the ``libgromacs`` library
and related CMake targets may change as we move towards building APIs
and supporting machinery that can be stable and supported in the long
term.

:issue:`3818`

Functionality deprecated in |Gromacs| 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generation of virtual sites to replace aromatic rings in standard residues
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

These are thought to produce artefacts under some circumstances
(unpublished results), were never well tested, are not widely used,
and we need to simplify pdb2gmx.

:issue:`3254`

Benchmarking options only available with ``gmx benchmark``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

Options such as ``-confout``, ``-resethway``, ``-resetstep`` are not
intended for use by regular mdrun users, so making them only available
with a dedicated tool is more clear. Also, this permits us to customize
defaults for e.g. writing files at the end of a simulation part in ways
that suit the respective mdrun and benchmark use cases, so ``-confout``
will no longer be required.

:issue:`3255`

``gmx mdrun -nsteps``
"""""""""""""""""""""

The number of simulation steps described by the .tpr file can be
changed with ``gmx convert-tpr``, or altered in .mdp file before the
call to ``gmx grompp``. The convenience of this mdrun option was
outweighed by the doubtful quality of its implementation, no clear
record in the log file, and lack of maintenance.

:issue:`3256`
