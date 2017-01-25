.. highlight:: bash

Build system overview
=====================

The |Gromacs| build system uses CMake (version
|CMAKE_MINIMUM_REQUIRED_VERSION| or newer is required) to generate the
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

Build types is a CMake concept that provides overall control of how
the build tools are used on the given platform to produce executable
code. These can be set in CMake in various ways, including on a
command line such as ``cmake -DCMAKE_BUILD_TYPE=Debug``. |Gromacs|
supports the following standard CMake build types:

**Release**
  Fully optimized code intended for use in production simulation. This is the
  default.

**Debug**
  Compiled code intended for use with debugging tools, with low optimization levels
  and debug information for symbols.

**RelWithDebInfo**
  As Release, but with debug information for symbol names, which can help debugging
  issues that only emerge in optimized code.

**MinSizeRel**
  As Release, but optimized to minimize the size of the resulting executable. This
  is never a concern for |Gromacs| installations, so should not be used, but
  probably works.

Additionally, |Gromacs| provides the following build types for development and
testing. Their implementations can be found in ``cmake/gmxBuildTypeXXX.cmake``.

**Reference**
  This build type compiles a version of |Gromacs| aimed solely at correctness. All
  parallelization and optimization possibilities are disabled. This build type is
  compiled with gcc 4.7 to generate the regression test reference values, against
  which all other |Gromacs| builds are tested.

**RelWithAssert**
  As Release, but removes ``-DNDEBUG`` from compiler command lines, which makes
  all assertion statements active (and can have other safety-related side effects
  in |Gromacs| and code upon which it depends)

**Profile**
  As Release, but adds ``-pg`` for use with profiling tools. This is not
  likely to be effective for profiling the performance of :ref:`gmx mdrun`, but can
  be useful for the tools.

**TSAN**
  Builds |Gromacs| for use with ThreadSanitzer in gcc >= 4.8 and clang
  >= 3.4 (http://clang.llvm.org/docs/ThreadSanitizer.html) to detect
  data races. This disables the use of atomics in ThreadMPI,
  preferring the mutex-based implementation.

**ASAN**
  Builds |Gromacs| for use with AddressSanitzer in gcc >= 4.8 and
  clang >= 3.4 (http://clang.llvm.org/docs/AddressSanitizer.html) to
  detect many kinds of memory mis-use. By default, AddressSanitizer
  includes LeakSanitizer.

**MSAN**
  Builds |Gromacs| for use with AddressSanitzer in clang >= 3.4
  (http://clang.llvm.org/docs/MemorySanitizer.html) to detect
  reads of unitialized memory. This functionality requires that
  dependencies of the |Gromacs| build have been built in a compatible
  way (roughly, static libraries with ``-g -fsanitize=memory
  -fno-omit-frame-pointer``), which generally requires at least the C++
  standard library to have been built specially. The path where the
  includes and libraries for dependencies should be found for this
  build type is set in the CMake cache variable
  ``GMX_MSAN_PATH``. Only internal XDR and internal fftpack are
  supported at this time.

For all of the sanitizer builds, to get readable stack traces, you may
need to ensure that the ``ASAN_SYMBOLIZER_PATH`` environment variable
(or your ``PATH``) include that of the ``llvm-symbolizer`` binary.

With some generators, CMake generates the build system for more than a
single ``CMAKE_BUILD_TYPE`` from one pass over the ``CMakeLists.txt``
files, so any code that uses ``CMAKE_BUILD_TYPE`` in
``CMakeLists.txt`` directly will break. |Gromacs| does use such CMake
code, so we do not fully support all these build types in such
generators (which includes Visual Studio).

CMake cache variables
---------------------

This section provides a (currently incomplete) list of cache variables that
developers or advanced users can set to affect what CMake generates and/or what
will get built.

.. TODO: Figure out where to document basic variables intended for user
   consumption, and how does it relate to documentation here.

.. TODO: Document the remaining variables below, and identify any variables
   missing from the list.

Compiler flags
^^^^^^^^^^^^^^

Standard CMake mechanism for specifying the compiler flags is to use
``CMAKE_C_FLAGS``/``CMAKE_CXX_FLAGS`` for flags that affect all build types,
and :samp:`CMAKE_C_FLAGS_{buildtype}`/:samp:`CMAKE_CXX_FLAGS_{buildtype}` for
flags that only affect a specific build type.  CMake provides some default flags.

|Gromacs| determines its own set of default flags, grouped into two categories:

* Generic flags that are appended to the above default CMake flag variables
  (possibly for multiple build types), generally specifying optimization flags
  to use and controlling compiler warnings.
* Specific flags for certain features that the build system determines to be
  necessary for successful compilation.  One example is flags that determine
  what SIMD instruction set the compiler is allowed to use/needs to support.

All of the above flags are only added after testing that they work with the
provided compiler.

There is one cache variable to control the behavior of automatic compiler flags:

.. cmake:: GMX_SKIP_DEFAULT_CFLAGS

   If set ``ON``, the build system will not add any compiler flags
   automatically (neither generic nor specific as defined above), and will skip
   most linker flags as well.
   The default flags that would have been added are instead printed out when
   :command:`cmake` is run, and the user can set the flags themselves using the
   CMake variables.
   If ``OFF`` (the default), the flags are added as described above.

The code the determine the default generic flags is in
:file:`cmake/gmxCFlags.cmake`.
Code that sets the specific flags (e.g., SIMD flags) is in the main
:file:`CMakeLists.txt`; search for :cmake:`GMX_SKIP_DEFAULT_CFLAGS`.
The variables used there can be traced back to the locations where the actual
flags to use are determined.

Variables affecting compilation/linking
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cmake:: GMX_BROKEN_CALLOC

.. cmake:: GMX_BUILD_FOR_COVERAGE

   Special variable set ``ON`` by Jenkins when doing a build for the coverage
   job.  Allows the build system to set options to produce as useful coverage
   metrics as possible.  Currently, it disables all asserts to avoid them
   showing up as poor conditional coverage.
   Defaults to ``OFF``, and there should not be any need to change this in a
   manual build.

   .. TODO: This could likely be replaced by a (yet another) build type.

.. cmake:: GMX_BUILD_MDRUN_ONLY

   If set ``ON``, the build system is configured to only build and install a
   single :file:`mdrun` executable.  To be fully functional, the installed
   :file:`mdrun` requires a standard |Gromacs| installation (with
   ``GMX_BUILD_MDRUN_ONLY=OFF``) in the same installation prefix, as the
   mdrun-only build does not install any data files or scripts, only the
   binary.  This is intended for cases where one wants to/needs to compile one
   or more instances of :file:`mdrun` with different build options (e.g., MPI
   or SIMD) than the full installation with the other utilities.
   Defaults to ``OFF``, in which case a single :file:`gmx` executable is built
   and installed, together with all the supporting files.  :command:`mdrun` can
   be executed as :command:`gmx mdrun`.

.. cmake:: GMX_BUILD_OWN_FFTW

.. cmake:: GMX_BUILD_SHARED_EXE

.. cmake:: GMX_COMPILER_WARNINGS

   If set ``ON``, various compiler warnings are enabled for compilers that
   Jenkins uses for verification.
   Defaults to ``OFF`` when building from a source tarball so that users
   compiling with versions not tested on Jenkins are not exposed to our rather
   aggressive warning flags that can trigger a lot of warnings with, e.g., new
   versions of the compilers we use.
   When building from a git repository, defaults to ``ON``.

.. cmake:: GMX_CYCLE_SUBCOUNTERS

   If set to ``ON``, enables performance subcounters that offer more
   fine-grained mdrun performance measurement and evaluation than the default
   counters. See :doc:`/user-guide/mdrun-performance` for the description of
   subcounters which are available.
   Defaults to ``OFF``.

.. cmake:: GMX_DATA_INSTALL_DIR

   Sets the directory under :file:`share/` where data files are installed.
   The default is ``gromacs``, which puts the files under
   file:`share/gromacs/`.
   See :doc:`relocatable-binaries` for how this influences the build.

.. cmake:: GMX_DOUBLE

   Many part of |Gromacs| are implemented in terms of "real" precision,
   which is actually either a single- or double-precision type,
   according to the value of this flag. Some parts of the code
   deliberately use single- or double-precision types, and these are
   unaffected by this setting. See reference manual for further
   information.

.. cmake:: GMX_RELAXED_DOUBLE_PRECISION

   Permit a double-precision configuration to compute some quantities
   to single-precision accuracy. Particularly on architectures where
   only double-precision SIMD is available (e.g. Sparc machines such
   as the K computer), it is faster to default to ``GMX_DOUBLE=ON``
   and use SIMD than to use ``GMX_DOUBLE=OFF`` and use no
   SIMD. However, if the user does not need full double precision,
   then some optimizations can achieve the equivalent of
   single-precision results (e.g. fewer Newton-Raphson iterations for
   a reciprocal square root computation).

.. cmake:: GMX_EXTRAE

.. cmake:: GMX_EXTERNAL_BLAS

.. cmake:: GMX_EXTERNAL_LAPACK

.. cmake:: GMX_EXTERNAL_TNG

.. cmake:: GMX_FFT_LIBRARY

.. cmake:: GMX_GIT_VERSION_INFO

   Whether to generate version information dynamically from git for each build
   (e.g., HEAD commit hash).
   Defaults to ``ON`` if the build is from a git repository and :command:`git`
   is found, otherwise ``OFF``.
   If ``OFF``, static version information from
   :file:`cmake/gmxVersionInfo.cmake` is used.

.. cmake:: GMX_GPU

.. cmake:: GMX_CLANG_CUDA

   Use clang for compiling CUDA GPU code, both host and device.

.. cmake:: GMX_CUDA_CLANG_FLAGS

    Pass additional CUDA-only compiler flags to clang using this variable.

.. cmake:: GMX_LIB_INSTALL_DIR

   Sets the installation directory for libraries (default is determined by
   standard CMake package ``GNUInstallDirs``).
   See :doc:`relocatable-binaries` for how this influences the build.

.. cmake:: GMX_LOAD_PLUGINS

.. cmake:: GMX_MPI

.. cmake:: GMX_OPENMP

.. cmake:: GMX_PREFER_STATIC_LIBS

.. cmake:: GMX_SIMD

.. cmake:: GMX_SOFTWARE_INVSQRT

.. cmake:: GMX_THREAD_MPI

.. cmake:: GMX_USE_RDTSCP

.. cmake:: GMX_USE_TNG

.. cmake:: GMX_VMD_PLUGIN_PATH

.. cmake:: GMX_X11

.. cmake:: GMX_XML

   Currently, this option has no effect on the compilation or linking, since
   there is no code outside the tests that would use :file:`libxml2`.

Variables affecting the ``all`` target
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cmake:: BUILD_TESTING

   Standard variable created by CTest that enables/disables all tests.
   Defaults to ``ON``.

.. cmake:: GMX_BUILD_HELP

   Controls handling of man pages and shell completions.  Possible values:

   ``OFF`` (default for builds from release source distribution)
     Man pages and shell completions are not generated as part of the ``all``
     target, and only installed if compiling from a source package.
   ``AUTO`` (default for builds from development version)
     Shell completions are generated by executing the :file:`gmx` binary as
     part of the ``all`` target.  If it fails, a message is printed, but the
     build succeeds.
     Man pages need to be generated manually by invoking the ``man`` target.
     Man pages and shell completions are installed if they have been
     successfully generated.
   ``ON``
     Works the same as ``AUTO``, except that if invoking the :file:`gmx` binary
     fails, the build fails as well.

.. cmake:: GMX_DEVELOPER_BUILD

   If set ``ON``, the ``all`` target will include also the test binaries using
   Google Test (if :cmake:`GMX_BUILD_UNITTESTS` is ``ON``).
   Also, :cmake:`GMX_COMPILER_WARNINGS` is always enabled.
   In the future, other developer convenience features (as well as features
   inconvenient for a general user) can be added to the set controlled by this
   variable.

Variables affecting special targets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. cmake:: CPPCHECK_XML_OUTPUT

   If set ``ON``, the ``cppcheck`` target generates reports for all found
   issues in XML format.  This is used by Jenkins, which parses the XML files
   to show the issues on the web page.
   If ``OFF`` (the default), issues are reported as plain text to standard
   output and to a text file.

.. cmake:: GMX_BUILD_MANUAL

   If set ``ON``, CMake detection for LaTeX and other prerequisites for the
   reference PDF manual is done, and the ``manual`` target for building the
   manual is generated.
   If ``OFF`` (the default), all detection is skipped and the manual cannot be
   built.

   .. TODO: Consider if this is really necessary, or if we could just use
      GMX_DEVELOPER_BUILD.

.. cmake:: GMX_BUILD_TARBALL

   If set ``ON``, ``-dev`` suffix is stripped off from version strings and some
   other version info logic is adjusted such that the man pages and other
   documentation generated from this build is suitable for releasing (on the
   web page and/or in the source distribution package).
   Defaults to ``OFF``.

.. cmake:: GMX_BUILD_UNITTESTS

   If ``ON``, test binaries using Google Test are built (either as the separate
   ``tests`` targer, or also as part of the ``all`` target, depending on
   :cmake:`GMX_DEVELOPER_BUILD`).  All dependencies required for building the
   tests (Google Test and Google Mock frameworks, and tinyxml2) are
   included in :file:`src/external/`.
   Defaults to ``ON`` if :cmake:`BUILD_TESTING` is ``ON``.

.. cmake:: GMX_COMPACT_DOXYGEN

   If set ``ON``, Doxygen configuration is changed to avoid generating large
   dependency graphs, which makes it significantly faster to run Doxygen and
   reduces disk usage.  This is typically useful when developing the
   documentation to reduce the build times.
   Defaults to ``OFF``.

.. cmake:: REGRESSIONTEST_DOWNLOAD

   If set ``ON``, CMake will download the regression tests and extract them to
   a local directory.  :cmake:`REGRESSIONTEST_PATH` is set to the extracted
   tests.  Note that this happens during the configure phase, not during the
   build.
   After the download is done, the variable is automatically reset to ``OFF``
   again to avoid repeated downloads.  Can be set to ``ON`` to download again.
   Defaults to ``OFF``.

.. cmake:: REGRESSIONTEST_PATH

   Path to extracted regression test suite matching the source tree (the
   directory containing :file:`gmxtest.pl`)
   If set, CTest tests are generated to run the regression tests.
   Defaults to empty.

.. cmake:: SOURCE_MD5SUM

   Sets the MD5 sum of the release tarball when generating the HTML
   documentation.  It gets inserted into the download section of the HTML
   pages.

External libraries
------------------

.. TODO: List external libraries used (either from src/external/, or from the
   system), whether they are required or optional, what functionality they
   provide for Gromacs, and how to control their use.

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

  All the above files get generated in :file:`src/`.

  Additionally, the following file is generated by the build system:

  :file:`baseversion-gen.c`
    Provides definitions for declarations in :file:`baseversion-gen.h` for
    version info output.  The contents are generated either from Git version
    info, or from static version info if not building from a git repository.
