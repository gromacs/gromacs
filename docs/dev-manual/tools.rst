Development-time tools
======================

Several tools have their own individual pages and are listed below.

.. toctree::
   :maxdepth: 2

   doxygen
   change-management
   infrastructure
   releng/index
   gmxtree
   code-formatting
   testutils
   physical_validation

.. todo:: :issue:`3032`

   Consider what is the most reasonable structure; currently, this list
   here does not make much sense in the overall organization and creates a
   confusing TOC for the developer guide.

.. todo:: :issue:`3267`

   Add details for most of the tools, either in the form of links to wiki,
   or to a separate page that explains more details.

Change management
-----------------

|Gromacs| change management uses git and `GitLab`_ for code uploading and testing as well as issues tracking.
(For change submission guidelines, refer to :doc:`contribute`.)

git
  |Gromacs| uses `git <https://git-scm.com/>`__ as the version control system.
  Instructions for setting up git for |Gromacs|, as well as tips and tricks for
  its use, can be found in :doc:`change-management`.

  Other basic tutorial material for ``git`` can be found on the `web <https://git-scm.com/doc/ext>`__.

GitLab
  Bugs and issues, as well as some random features and discussions,
  are tracked, and all code changes go through a code review system at
  https://gitlab.com/gromacs/gromacs.

Build testing
  All changes pushed to GitLab are automatically compiled and otherwise
  checked on various platforms.
  :doc:`infrastructure` documents how builds are automated,
  providing information on how to replicate the builds (e.g., to
  diagnose issues).
  :doc:`releng/index` provides more information on the technical implementation
  of the builds.

Build system
------------

.. todo:: details, ASAN, others?

CMake
  Main tool used in the build system.

packaging for distribution (CPack)

unit testing (CTest)
  |Gromacs| uses a unit testing framework based on Google C++ Testing
  Framework (gtest) and CTest.  All unit tests are automatically run in GitLab CI
  for each commit.
  Details can be found on a separate page on :doc:`testutils`.

clang static analyzer

coverage

regression tests

floating-point exceptions
  In debug builds, floating-point exceptions (FPEs) are generated whenever one of the
  following operations is encountered: division by zero, floating-point overflow,
  invalid operation (e.g., taking sqrt of a negative number).
  Such checks are *not* performed in the following configurations:

  - release build,
  - any build by Clang with optimizations,
  - build with SYCL support.

  In these configurations, FPEs can be enabled by adding ``-fpexcept`` flag to ``gmx``
  invocation. However, FPEs are not supported on Windows and non-x86 Apple hardware.
  See ``api/legacy/include/gromacs/math/utilities.h`` for more details.

.. _dev-formatting-tools:

Code formatting and style
-------------------------

The tools and scripts listed below are used to automatically check/apply
formatting that follows |Gromacs| style guidelines described on a separate page:
:doc:`style`.

clang-format
  We use clang-format to enforce a consistent coding style, with the
  settings recorded in ``.clang-format`` in the main tree.
  See :ref:`gmx-clang-format` for details.

clang-tidy
  The source code linter clang-tidy is used to enforce common restrictions to the
  code, with the checks collected under ``.clang-tidy`` at the top of the main tree.
  See :ref:`gmx-clang-tidy` for details.

``admin/copyright.py``
  This Python script adds and formats copyright headers in source files.
  ``copyright.sh`` (see below) uses the script to check/update copyright years on
  changed files automatically.

``admin/copyright.sh``
  This ``bash`` script runs the ``copyright.py`` python script to enforce
  correct copyright information in all files that have local changes
  and checks that they conform to the prescribed
  style.  Optionally, the script can also apply changes to make the files
  conform.
  This script is automatically run by the CI to ensure that all commits adhere
  to :doc:`formatting`.  If the copyright job does not succeed, it
  means that this script has something to complain.
  See :doc:`code-formatting` for details.

``admin/clang-format.sh``
  This script enforces coding style using clang-format.
  This script is automatically run by our CI to ensure that all commits adhere
  to :doc:`formatting`.

``admin/clang-tidy.sh``
  The clang-tidy code correctness restrictions are enforced by this script.
  The script is also used by the CI to verify the code, in addition to nightly
  compilations using clang-tidy on the whole tree.

``admin/git-pre-commit``
  This sample git pre-commit hook can be used if one wants to apply
  ``clang-tidy.sh``, ``copyright.sh`` and ``clang-format.sh`` automatically
  before every commit to check for formatting
  issues.  See :doc:`code-formatting` for details.

include directive checker
  In its present form, the above include sorter script cannot be conveniently
  applied in the formatting script.  To check for issues, it is instead integrated into
  a ``check-source`` build target.  When this target is built, it also checks for
  include formatting issues.  Internally, it uses the sorter script.  This check
  is run in the CI as part of the Documentation job.
  Details for the checking mechanism are on a separate page (common for several
  checkers): :doc:`gmxtree`.

``admin/reformat_all.sh``
  This ``bash`` script runs clang-format/``copyright.py``/include sorter
  on all relevant files in the source tree (or in a particular directory).
  The script can also produce the list of files where these scripts are applied,
  for use with other scripts.  See :doc:`code-formatting` for details.

git attributes
  git attributes (specified in ``.gitattributes`` files) are used to annotate
  which files are subject to automatic formatting checks (and for automatic
  reformatting by the above scripts).  See ``man gitattributes`` for an overview of
  the mechanism.  We use the ``filter`` attribute to specify the type of automatic
  checking/formatting to apply.  Custom attributes are used for specifying some
  build system dependencies for easier processing in CMake.

