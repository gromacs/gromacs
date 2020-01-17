Understanding Jenkins builds
============================

This page documents what different Jenkins builds actually run from the
|Gromacs| source tree.  The purpose is two-fold:

* Provide information on how to interpret Jenkins failures and how to run the
  same tasks locally to diagnose issues (in most cases, referring to the
  special targets described in :doc:`build-system`).
* Provide information on what changes in the build system (or other parts of
  the repository) need special care to not break Jenkins builds.

.. todo:: Add a link to a wiki page about general Jenkins documentation, once
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

uncrustify
^^^^^^^^^^

This build checks for source code formatting issues with uncrustify, and enforces
the copyright style.  See :doc:`formatting` for the guidelines that are enforced.

The exact build sequence is in :file:`admin/builds/uncrustify.py`, which
essentially just runs ::

  admin/uncrustify.sh check --rev=HEAD^

If the any changes are required, the build is marked unstable.
If the script completely fails (should be rare), the build fails.
A file with issues found by the script is archived as an artifact in the build,
and a summary is reported back to Gerrit (or the actual issues if there are
only a few).
See :doc:`code-formatting` for more details on code-formatting tools
and on scripts to run them.

clang-format
^^^^^^^^^^^^

This build checks and enforces code formatting, e.g.,  indentation.
Also, a second part of the build enforces the source code formatting.
As above, see :doc:`formatting` for the style guidelines.

The build runs according to :file:`admin/builds/clang-format.py`, resulting
in running ::

 admin/clang-format.sh check --rev=HEAD^

The build is marked unstable if the code formatting resulted in
any changes to the source code.

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

Updating regressiontests data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes we add new tests to the regressiontests repository. Also, as
the source code or data files change, it is sometimes necessary to
update regressiontests. This requires a particular CMake build type
and both a single and double-precision build of |Gromacs| to generate
all the data. Jenkins can automate much of the tedium here.

* Upload a regressiontests change that lacks the relevant reference
  data (either because you deleted the outdated data, or because the
  test is new). Jenkins will do the normal thing, which we ignore.
  There is now a Gerrit patch number for that change, symbolized here
  with ``MMMM``.

* Go to change ``MMMM`` on gerrit, select the patch set you want to
  update with new reference data (usually the latest one), and comment

    ``[JENKINS] Update``

  to update against the HEAD of the matching source-code branch, or

    ``[JENKINS] Cross-verify NNNN update``

  to update from builds of |Gromacs| from the latest version of
  Gerrit source-code patch ``NNNN``. You will need to do this when
  functionality changes in ``NNNN`` affect either the layout of
  the files in the reference data, or the results of the simulation,
  or the results of the subsequent analysis.

* Eventually, Jenkins will upload a new version of the regressiontests
  patch to Gerrit, which will contain the updated regressiontest data.
  That upload will again trigger Jenkins to do the normal pre-submit
  verify, which will now pass (but perhaps will only pass under
  cross-verify with patch ``NNNN``, as above).

* Later, if you later need to verify an updated version of source-code
  patch ``NNNN`` against the newly generated reference data, go to the
  source-code patch ``NNNN`` and comment

    ``[JENKINS] Cross-verify MMMM``
