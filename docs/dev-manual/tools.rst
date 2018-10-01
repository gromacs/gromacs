Development-time tools
======================

Several tools have their own individual pages and are listed below.

.. toctree::
   :maxdepth: 2

   doxygen
   change-management
   jenkins
   releng/index
   gmxtree
   uncrustify
   testutils
   physical_validation

.. TODO: Consider what is the most reasonable structure; currently, this list
   here does not make much sense in the overall organization and creates a
   confusing TOC for the developer guide.

.. TODO: Add details for most of the tools, either in the form of links to wiki,
   or to a separate page that explains more details.

Change management
-----------------

|Gromacs| change management is supported by the following tools.
(For change submission guidelines, refer to :doc:`contribute`.)

git
  |Gromacs| uses `git <https://git-scm.com/>`__ as the version control system.
  Instructions for setting up git for |Gromacs|, as well as tips and tricks for
  its use, can be found in :doc:`change-management`.

  Other basic tutorial material for ``git`` can be found on the `web <https://git-scm.com/doc/ext>`__.

Gerrit
  All code changes go through a code review system at
  http://gerrit.gromacs.org.

Jenkins
  All changes pushed to Gerrit are automatically compiled and otherwise
  checked on various platforms using a continuous integration system at
  http://jenkins.gromacs.org.
  :doc:`jenkins` documents how Jenkins interacts with the build system,
  providing information on how to replicate the builds Jenkins does (e.g., to
  diagnose issues).
  :doc:`releng/index` provides more information on the technical implementation
  of the builds.

Redmine
  Bugs and issues, as well as some random features and discussions,
  are tracked at http://redmine.gromacs.org.

.. _Git Tips & Tricks: http://www.gromacs.org/index.php?title=Developer_Zone/Git/Git_Tips_%26_Tricks

Build system
------------

.. TODO: details, ASAN, others?

CMake
  Main tool used in the build system.

ccache
  When the `ccache <https://ccache.samba.org>`_ caching compiler wrapper is
  found on the PATH, we attempt to set up a caching compiler wrapper for CMake
  builds. Not all compilers are supported. Refer to the ``ENABLE_CCACHE``
  option in ``CMakeLists.txt`` and to ``cmake/gmxCcache.cmake``
  for details. Please submit updates if you find that the current
  configuration is too conservative.

packaging for distribution (CPack)

unit testing (CTest)
  |Gromacs| uses a unit testing framework based on Google C++ Testing
  Framework (gtest) and CTest.  All unit tests are automatically run on Jenkins
  for each commit.
  Details can be found on a separate page on :doc:`testutils`.

regression tests

clang-tidy
  `clang-tidy <http://releases.llvm.org/7.0.0/tools/clang/tools/extra/docs/clang-tidy/index.html>`_
  is used for static code analysis. clang-tidy is easy to install. It is contained in
  the llvm binary `package <http://releases.llvm.org/download.html#6.0.0>`_. Only
  version 7.0.* with libstdc++<7 or libc++ is supported. Others might miss tests or give false positives.
  It is run automatically on Jenkins for each commit. Many checks have fixes which can automatically be
  applied. To run it, the build has to be configured with
  ``cmake -DGMX_CLANG_TIDY=ON -DGMX_OPENMP=no -DCMAKE_BUILD_TYPE=Debug -DCMAKE_EXPORT_COMPILE_COMMANDS=on``.
  Any ``CMAKE_BUILD_TYPE`` which enables asserts (e.g. ASAN) works. Such a configured build will
  run both the compiler as well as clang-tidy when building. The name of the clang-tidy executable is set with
``-DCLANG_TIDY=...``, and the full path to it can be set with ``-DCLANG_TIDY_EXE=...``.
  To apply the automatic fixes to the issue identified clang-tidy should be run sepereately (running clang-tidy
  with ``-fix`` as part of the build can corrupt header files). To fix a specific file run
  ``clang-tidy -fix -header-filter '.*' {file}``, to fix all files in parallel
  ``run-clang-tidy.py -fix -header-filter '.*' '(?<!/selection/parser\.cpp|selection/scanner\.cpp)$'``,
  and to fix all modified files ``run-clang-tidy.py -fix -header-filter '.*' $(git diff HEAD --name-only)``.
  The run-clang-tidy.py script is in the
  ``share/clang/`` subfolder of the llvm distribution. ``clang-tidy`` has to be able to find the
  ``compile_commands.json`` file. Eithe run from the build folder or add a symlink to the source folder.

clang static analyzer

coverage

.. _dev-formatting-tools:

Code formatting and style
-------------------------

The tools and scripts listed below are used to automatically check/apply
formatting that follows |Gromacs| style guidelines described on a separate page:
:doc:`style`.

uncrustify
  `uncrustify <http://uncrustify.sourceforge.net>`_ is used for automatic
  indentation and other formatting of the source code to follow
  :doc:`formatting`.  All code must remain invariant under uncrustify
  with the config at ``admin/uncrustify.cfg``.  A patched version of uncrustify is
  used.  See :doc:`uncrustify` for details.

``admin/copyright.py``
  This Python script adds and formats copyright headers in source files.
  ``uncrustify.sh`` (see below) uses the script to check/update copyright years on
  changed files automatically.

``admin/uncrustify.sh``
  This ``bash`` script runs uncrustify and ``copyright.py`` for all
  files that have local changes and checks that they conform to the prescribed
  style.  Optionally, the script can also apply changes to make the files
  conform.
  This script is automatically run by Jenkins to ensure that all commits adhere
  to :doc:`formatting`.  If the uncrustify job does not succeed, it
  means that this script has something to complain.
  See :doc:`uncrustify` for details.

``admin/git-pre-commit``
  This sample git pre-commit hook can be used if one wants to apply
  ``uncrustify.sh`` automatically before every commit to check for formatting
  issues.  See :doc:`uncrustify` for details.

``docs/doxygen/includesorter.py``
  This Python script sorts and reformats #include directives according to
  the guidelines at :doc:`includestyle`.  Details are documented on a
  separate page (with the whole suite of Python scripts used for source code
  checks): :ref:`dev-include-sorter`.

include directive checker
  In its present form, the above include sorter script cannot be conveniently
  applied in ``uncrustify.sh``.  To check for issues, it is instead integrated into
  a ``check-source`` build target.  When this target is built, it also checks for
  include formatting issues.  Internally, it uses the sorter script.  This check
  is run in Jenkins as part of the Documentation job.
  Details for the checking mechanism are on a separate page (common for several
  checkers): :doc:`gmxtree`.

``admin/reformat_all.sh``
  This ``bash`` script runs uncrustify/``copyright.py``/include sorter
  on all relevant files in the source tree (or in a particular directory).
  The script can also produce the list of files where these scripts are applied,
  for use with other scripts.  See :doc:`uncrustify` for details.

git attributes
  git attributes (specified in ``.gitattributes`` files) are used to annotate
  which files are subject to automatic formatting checks (and for automatic
  reformatting by the above scripts).  See ``man gitattributes`` for an overview of
  the mechanism.  We use the ``filter`` attribute to specify the type of automatic
  checking/formatting to apply.  Custom attributes are used for specifying some
  build system dependencies for easier processing in CMake.

include-what-you-use
