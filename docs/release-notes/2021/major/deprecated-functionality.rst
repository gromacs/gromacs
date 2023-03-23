Deprecated functionality
------------------------

.. Note to developers!
   Please use """"""" to underline the individual entries for fixed issues in the subfolders,
   otherwise the formatting on the webpage is messed up.
   Also, please use the syntax :issue:`number` to reference issues on GitLab, without the
   a space between the colon and number!

The core |Gromacs| team wants to let users and downstream developers
know about impending changes so that disruption is minimized. Do get
in touch if you feel something inappropriate is planned!

Deprecated functionality often remains in |Gromacs| for a year or
more, but this should not be relied upon.

.. Note to maintainers!
   The sections below should general copy the contents from the previous major release,
   except where appropriate when code or planning changes have happened.

Changes anticipated to |Gromacs| 2021 functionality
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``gmx mdrun -membed``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The feature for embedding a protein in a membrane will be retained,
but probably in a different form, such as ``gmx membed``.

``gmx mdrun -rerun``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The feature for computing potential energy quantities from a
trajectory will be retained, but probably in a different form, such as
``gmx rerun`` and ``gmx test-particle-insertion``.

Integrator .mdp options will only contain dynamical integrators
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Energy minimization will be accessed in a differt form, perhaps with
``gmx minimize`` and interpret an .mdp field for which minimizer to
use. Normal-mode analysis may be accessed with e.g. ``gmx
normal-modes``. The command-line help for these tools will then
be better able to document which functionality is supported when.

Much functionality in ``trjconv``, ``editconf``, ``eneconv`` and ``trjcat``
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The functionality in such tools is being separated to make it
available in composable modules, that we plan to make available as
simpler tools, and eventually via the |Gromacs| API that is under
development.

``gmx do_dssp`` to be replaced
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This tool is deprecated, because it is problematic for some users to
obtain and install a separate DSSP binary, so we plan to replace the
implementation at some point with a native implementation, likely
based upon xssp, and make it available under a new gmx tool name.

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

Intel KNC (MIC) support
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3818` This architecture is nearly extinct in HPC. Note that
KNL support will continue and is not affected by this deprecation.

Sparc64 HPC ACE
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This architecture is nearly extinct in HPC.

Legacy SIMD architecture support
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3818` We occasionally need to extend the |Gromacs| SIMD
framework, and so should slowly remove older architectures that are
difficult or impossible to test. The following implementations are
deprecated and will not support new functionality in future.

* Power 7
* ARMv7 (this platform was deprecated in |Gromacs| 2020)
* x86 MIC (this platform was deprecated in |Gromacs| 2021)
* Sparc64 HPC ACE  (this platform was deprecated in |Gromacs| 2021)

The mdrun-only build of |Gromacs|
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`3808` Before |Gromacs| had the ``gmx`` wrapper binary, the
``mdrun`` binary could be built independently of the many other binary
tools that were built by default. That was useful for installing on
compute clusters because dependencies for ``mdrun`` were
minimized. However, we now manage such dependencies better with CMake,
and an mdrun-only build is no longer easier to build. The mdrun-only
build is also harder to test, and introduces complexity into
documenting |Gromacs| and teaching users to use it. So it is time to
remove that build.

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

Constant-acceleration MD
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
:issue:`1354` This has been broken for many years, and will be removed
as nobody has been found with interest to fix it.

Reading .pdo files in ``gmx wham``
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The pull code in |Gromacs| before version 4.0 wrote files in ``.pdo``
format. Analyses of such files are likely no longer relevant, and if
they are, using any older |Gromacs| version will work. ``gmx wham`` will be
simpler to maintain and extend if we no longer support reading
``.pdo`` files.

Functionality deprecated in |Gromacs| 2020
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Support for 32bit architectures
"""""""""""""""""""""""""""""""
:issue:`3252` There are no current or planned large scale resources using 32bit architectures,
and we have no ability to properly test and evaluate them.

Free-energy soft-core power 48
""""""""""""""""""""""""""""""
:issue:`3253` Free-energy soft-core power 48 is almost never used and is therefore deprecated.

Support for Armv7
"""""""""""""""""
:issue:`2990` There are several issues with current code for the architecture, and we don't
have the resources for support and fix issues related to it. As the architecture has no
large HPC impact it is thus deprecated.

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

