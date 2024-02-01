Deprecated functionality
------------------------

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

Changes anticipated to |Gromacs| 2022 functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Functionality deprecated in |Gromacs| 2022
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

GMX_OPENCL_NB_CLUSTER_SIZE CMake variable deprecated in favor of GMX_GPU_NB_CLUSTER_SIZE
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Both OpenCL and SYCL support different cluster sizes, so GMX_GPU_NB_CLUSTER_SIZE should
be used going forward.

The built-in viewer ``gmx view`` will be removed
""""""""""""""""""""""""""""""""""""""""""""""""

There is little use and no tests of this functionality, so it is not worth attempting to
maintain moving forward.

:issue:`4296`

The analysis tool ``gmx chi`` will be removed
"""""""""""""""""""""""""""""""""""""""""""""

This tool has not been functional for a few years.
Please comment at the linked issue if you have any interest in it.

:issue:`4108`

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
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

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
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3818` Work is underway for ports to AMD and Intel GPUs, and it
is likely that those ports will not be based on the current |Gromacs|
OpenCL port. Nvidia GPUs are targeted by the CUDA port, and no changes
are expectd there. The core team can't maintain, test, and extend up
to 4 ports with current resource levels. Since there are no prospects
of an emerging GPU vendor in HPC needing OpenCL support, we will
remove the OpenCL port once AMD and Intel support is established in
other ways.

Support for version 1 of the hardware locality library ``hwloc``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3818` Version 2 has been supported in |Gromacs| for several
years. The capabilities of newer hardware and hardware-support APIs
are of most interest for |Gromacs| moving forward, so we should
minimize our testing work and encourage clusters to upgrade older
``hwloc`` installations.

Legacy API
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3818` The legacy installed headers have been deprecated for a
while, however we wish to state more broadly that all headers found
within the ``src`` directory tree of |Gromacs| are intended for
internal consumption only, and are thus subject to change without
notice. Further, the form and contents of the ``libgromacs`` library
and related CMake targets may change as we move towards building APIs
and supporting machinery that can be stable and supported in the long
term.

Functionality deprecated in |Gromacs| 2019
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Generation of virtual sites to replace aromatic rings in standard residues
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3254` These are thought to produce artefacts under some circumstances
(unpublished results), were never well tested, are not widely used,
and we need to simplify pdb2gmx.

Benchmarking options only available with ``gmx benchmark``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3255` Options such as ``-confout``, ``-resethway``, ``-resetstep`` are not
intended for use by regular mdrun users, so making them only available
with a dedicated tool is more clear. Also, this permits us to customize
defaults for e.g. writing files at the end of a simulation part in ways
that suit the respective mdrun and benchmark use cases, so ``-confout``
will no longer be required.

``gmx mdrun -nsteps``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3256` The number of simulation steps described by the .tpr file can be
changed with ``gmx convert-tpr``, or altered in .mdp file before the
call to ``gmx grompp``. The convenience of this mdrun option was
outweighted by the doubtful quality of its implementation, no clear
record in the log file, and lack of maintenance.

