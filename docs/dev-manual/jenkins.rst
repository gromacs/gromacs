Understanding Jenkins builds
============================

This page documents what different Jenkins builds actually run from the
|Gromacs| source tree.  The purpose is two-fold:

* Provide information on how to interpret Jenkins failures and how to run the
  same tasks locally to diagnose issues (in most cases, referring to the
  special targets described in :doc:`build-system`).
* Provide information on what changes in the build system (or other parts of
  the repository) need special care to not break Jenkins builds.

Separate page documents how to interact with the Jenkins UI for these builds:
:doc:`releng/jenkins-ui`.
:doc:`releng/jenkins-howto` has information on how to do common things with
Jenkins builds.

.. TODO: Add a link to a wiki page about general Jenkins documentation, once
   there is more of that.

Pre-submit verification
-----------------------

The following builds are triggered for each patch set uploaded to Gerrit.

Compilation and tests
^^^^^^^^^^^^^^^^^^^^^

The main build compiles |Gromacs| with different configurations and runs the
tests.  The configurations used for Jenkins verification are specified in
:file:`admin/builds/pre-submit-matrix.txt`.

The exact build sequence can be found in :file:`admin/builds/gromacs.py`,
including the logic that translates the build options in the matrix file to
CMake options.

Documentation
^^^^^^^^^^^^^

This build builds various types of documentation:

* PDF reference manual using LaTeX
* Doxygen documentation extracted from the source code
* Set of HTML pages containing an installation guide, a user guide, and a
  developer guide, as well as links to the above.  This set of HTML pages can
  be browsed from Jenkins.
* Man pages
* INSTALL text file

The last three require building the :file:`gmx` binary and running it, so
compilation failures will also show in this build.
All log files that contain warnings are archived as artifacts in the build, and
presence of any warnings marks the build unstable.  Brief description of which
part failed is reported back to Gerrit.

Additionally, the build runs some source code checks that rely on the Doxygen
documentation.  See the description of the ``check-source`` target in
:doc:`gmxtree`.

:doc:`doxygen` provides general guidelines for Doxygen usage, which can be
helpful in understanding and solving Doxygen warnings and some of the
``check-source`` issues.
:doc:`includestyle` provides guidelines for #include order and style, which is
another part of ``check-source`` checks.

The exact build sequence is in :file:`admin/builds/documentation.py`.
See that file for details of what it exactly builds and how.  Most changes in the
documentation build system will require changes in this script, but Jenkins
configuration should be more static.

clang static analysis
^^^^^^^^^^^^^^^^^^^^^

The file :file:`admin/builds/clang-analyzer.py` specifies the exact build
sequence and the CMake cache variables used for clang static analysis.  This
file also specifies the clang version used for the analysis, as well as the C++
compiler used (``clang-static-analyzer-<version>``).

To run the analysis outside Jenkins, you should run both ``cmake`` and ``make``
under ``scan-build`` command using the same CMake cache variables as in the
build script. When you do the initial CMake configuration with ``scan-build``,
it sets the C++ compiler to the analyzer. Note that using ``scan-build`` like
this will also analyze C code, but Jenkins ignores C code for analysis. This
can result in extra warnings, which can be suppressed by manually setting
CMAKE_C_COMPILER to a value other than Clang static analyzer.

cppcheck
^^^^^^^^

This build runs the :command:`cppcheck` static analysis tool.  Any issues found
mark the build unstable, and can be browsed in Jenkins.

It runs :command:`cmake` to generate the build system, and then builds the
``cppcheck`` target.  Nothing is compiled by this target, it only runs
:command:`cppcheck` for the designated source files.  The CMake configuration
options do not affect the set of files checked, but they do affect the checked
code through :file:`config.h` and such.

The exact build sequence and the CMake configuration used is in
:file:`admin/builds/cppcheck.py`.

uncrustify
^^^^^^^^^^

This build checks the source code for formatting such as consistent indentation
and use of braces, as well as for copyright headers.  See :doc:`formatting` for
the guidelines that are enforced.

The exact build sequence is in :file:`admin/builds/uncrustify.py`, which
essentially just runs ::

  admin/uncrustify.sh check --rev=HEAD^

If the any changes are required, the build is marked unstable.
If the script completely fails (should be rare), the build fails.
A file with issues found by the script is archived as an artifact in the build,
and a summary is reported back to Gerrit (or the actual issues if there are
only a few).
See :doc:`uncrustify` for more details on uncrustify and on scripts to run it.

On-demand builds
----------------

These builds can be triggered on request for certain changes in Gerrit, or
manually from Jenkins.  See :ref:`releng-triggering-builds` for details on
how to trigger these.

Coverage
^^^^^^^^

This build compiles one configuration of |Gromacs| with instrumentation for
coverage, runs the tests, and produces a coverage report using gcovr.
The report can be browsed on Jenkins.

The exact build sequence is in :file:`admin/builds/coverage.py`, including
specification of the configuration tested.

Source tarball
^^^^^^^^^^^^^^

This build creates the source tarball for distribution.  Some of the content
that is put into the tarball is generated by executing the :command:`gmx`
binary, so this build also compiles the source code (with a minimal set of
options).

The build compiles the code and those targets that generate content necessary
for the tarball, followed by building the ``package_source`` target.
After that, it just generates a file that is used by other builds.

The exact build sequence is in :file:`admin/builds/source-package.py`.

Release workflow
^^^^^^^^^^^^^^^^

This build creates source and regressiontest tarballs, builds, installs, and
tests a few configuration using those, and builds documentation to be placed on
the documentation web site for a new release.  The set of configurations tested
is specified in :file:`admin/builds/release-matrix.txt`.

The exact build sequence is desribed in :ref:`releng-workflow-release`.
The build uses the source tarball build as a subbuild, and parts of the build
are executed using :file:`admin/builds/gromacs.py` and
:file:`admin/builds/documentation.py`.

:file:`admin/builds/get-version-info.py` is used for getting the version
information from the source tree as part of this workflow.

:file:`admin/builds/update-regtest-hash.py` has logic to update the
regressiontests tarball MD5 sum for the released tarball automatically.
