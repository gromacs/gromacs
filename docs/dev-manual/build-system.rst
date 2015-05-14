.. highlight:: bash

Build system overview
=====================

The |Gromacs| build system uses CMake (version
|GMX_CMAKE_MINIMUM_REQUIRED_VERSION| or newer is required) to generate the
actual build system for the build tool choosen by the user.  See CMake
documentation for general introduction to CMake and how to use it.  This
documentation focuses on how the |Gromacs| build system is organized and
implemented, and what features it provides to developers (some of which may be
of interest to advanced users).

Most developers use ``make`` or ``ninja`` as the underlying build system, so
there can be parts of the build system that are specifically designed for
command-line consumption with these tools, and may not work ideally with other
environments, but basic building should be possible with all the environments
supported by CMake.

Also, the build system and version control is designed for out-of-source
builds.  In-source builds mostly work (there are a few custom targets that do
not), but no particular effort has been put to, e.g., having :file:`.gitignore`
files that exclude all the build outputs, or to have the ``clean`` target
remove all possible build outputs.

Build types
-----------

CMake cache variables
---------------------

This section provides an (currently incomplete) list of cache variables that
developers or advanced users can set to affect what CMake generates and/or what
will get built.

.. TODO: Figure out where to document basic variables intended for user
   consumption, and how does it relate to documentation here.

Compiler flags
^^^^^^^^^^^^^^

.. cmake:: GMX_SKIP_DEFAULT_CFLAGS

Variables affecting compilation/linking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cmake:: GMX_BUILD_FOR_COVERAGE

.. cmake:: GMX_BUILD_MDRUN_ONLY

.. cmake:: GMX_BUILD_OWN_FFTW

.. cmake:: GMX_DATA_INSTALL_DIR

.. cmake:: GMX_EXTERNAL_BLAS

.. cmake:: GMX_EXTERNAL_LAPACK

.. cmake:: GMX_EXTERNAL_BOOST

.. cmake:: GMX_EXTERNAL_TNG

.. cmake:: GMX_FFT_LIBRARY

.. cmake:: GMX_GIT_VERSION_INFO

.. cmake:: GMX_LOAD_PLUGINS

.. cmake:: GMX_MPI

.. cmake:: GMX_SIMD

.. cmake:: GMX_THREAD_MPI

.. cmake:: GMX_USE_RDTSCP

Variables affecting the ``all`` target
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cmake:: GMX_BUILD_HELP

.. cmake:: GMX_DEVELOPER_BUILD

.. cmake:: GMX_LIB_INSTALL_DIR

.. cmake:: REGRESSIONTEST_DOWNLOAD

.. cmake:: REGRESSIONTEST_PATH

Variables affecting special targets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cmake:: CPPCHECK_XML_OUTPUT

.. cmake:: GMX_BUILD_MANUAL

.. cmake:: GMX_BUILD_TARBALL

.. cmake:: GMX_BUILD_UNITTESTS

.. cmake:: GMX_COMPACT_DOXYGEN

.. cmake:: SOURCE_MD5SUM

Special targets
---------------

In addition to the default ``all`` target, the generated build system has
several custom targets that are intended to be explicitly built to perform
various tasks (some of these may also run automatically).  There are various
other targets as well used internally by these, but those are typically not
intended to be invoked directly.

check
   Builds all the binaries needed by the tests and runs the tests.  If some
   types of tests are not available, shows a note to the user.
   This is the main target intended for normal users to run the tests.
   See :doc:`testutils`.
check-source
   Runs a custom Python checker script to check for various source-level
   issues.  Uses Doxygen XML documentation as well as rudimentary parsing of
   some parts of the source files.
   This target is used as part of the Jenkins documentation job.
   All CMake code is currently in :file:`docs/doxygen/`.
   See :doc:`gmxtree`.
completion
   Runs the compiled :file:`gmx` executable to generate shell command-line
   completion definitions.  This target is only added if
   :cmake:`GMX_BUILD_HELP` is not ``OFF``, and it is run automatically as part
   of the default ``all`` target.  See :cmake:`GMX_BUILD_HELP`.
   All CMake code is in :file:`src/programs/`.
cppcheck
   Runs :command:`cppcheck` with the flags used in Jenkins for all the source
   files.  This target is directly used by the Jenkins cppcheck job.
   All CMake code is in :file:`tests/CppCheck.cmake`.
dep-graphs*
   Builds include dependency graphs for the source files using :command:`dot`
   from graphviz.
   All CMake code is in :file:`docs/doxygen/`.
   See :doc:`gmxtree`.
doxygen-*
   Targets that run Doxygen to generate the documentation.
   The ``doxygen-all`` target runs as part of the ``webpage`` target, which in
   turn runs as part of the Jenkins documentation job.
   All CMake code is in :file:`docs/doxygen/`.
   See :doc:`doxygen`.
install-guide
   Runs Sphinx to generate a plain-text INSTALL file for the source package.
   The files is generated at :file:`docs/install-guide/text/`, from where it
   gets put at the root of the source package by CPack.
   All CMake code is in :file:`docs/`.
man
   Runs Sphinx to generate man pages for the programs.
   Internally, also runs the compiled :file:`gmx` executable to generate the
   input files for Sphinx.
   All CMake code is in :file:`docs/`.
   See :cmake:`GMX_BUILD_HELP` for information on when the man pages are
   installed.
manual
   Runs LaTeX to generate the reference PDF manual.
   All CMake code is in :file:`docs/manual/`.
   See :cmake:`GMX_BUILD_MANUAL`.
package_source
   Standard target created by CPack that builds a source package.
   This target is used to generate the released source packages.
test
   Standard target created by CTest that runs all the registered tests.
   Note that this does not build the test binaries, only runs them, so you need
   to first ensure that they are up-to-date.
   See :doc:`testutils`.
tests
   Builds all the binaries needed by the tests (but not ``gmx``).
   See :doc:`testutils`.
webpage
   Collection target that runs the other documentation targets to generate the
   full set of HTML (and linked) documentaion.
   This target is used as part of the Jenkins documentation job.
   All CMake code is in :file:`docs/`.
webpage-sphinx
   Runs Sphinx to generate most content for the HTML documentation (the set of
   web pages this developer guide is also part of).
   Internally, also runs the compiled :file:`gmx` executable to generate some
   input files for Sphinx.
   All CMake code is in :file:`docs/`.

Passing information to source code
----------------------------------

The build system uses a few different mechanisms to influence the compilation:

* On the highest level, some CMake options select what files will be compiled.
* Some options are passed on the compiler command line using ``-D`` or
  equivalent, such that they are available in every compilation unit.  This
  should be used with care to keep the compiler command lines manageable.
  You can find the current set of such defines with ::

    git grep add_definitions

* A few header files are generated using CMake ``configure_file()`` and
  included in the desired source files.  These files must exist for the
  compilation to pass.  Only a few files use an ``#ifdef HAVE_CONFIG_H`` to
  protect against inclusion in case the define is not set; this is used in
  files that may get compiled outside the main build system.

  :file:`buildinfo.h`
    Contains various strings about the build environment, used mainly for
    outputting version information to log files and when requested.
  :file:`config.h`
    Contains defines for conditional compilation within source files.
  :file:`gmxpre-config.h`
    Included by :file:`gmxpre.h` as the first thing in every source file.
    Should only contain defines that are required before any other header for
    correct operation.  For example, defines that affect the behavior of system
    headers fall in this category.  See Doxygen documentation for
    :file:`gmxpre.h`.

  All the files get generated in :file:`src/`.
