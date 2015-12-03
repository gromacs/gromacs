Development-time tools
======================

.. toctree::
   :maxdepth: 2

   doxygen
   jenkins
   gmxtree
   uncrustify
   testutils

.. TODO: Consider what is the most reasonable structure; currently, this list
   here does not make much sense in the overall organization and creates a
   confusing TOC for the developer guide.

.. TODO: Add details for most of the tools, either in the form of links to wiki,
   or to a separate page that explains more details.

Change management
-----------------

git
  |Gromacs| uses git as the version control system.
  Instructions for setting up git for |Gromacs|, as well as tips and tricks for
  its use, can be found on the wiki: `Git Tips & Tricks`_

  Other basic tutorial material for ``git`` can be found on the web.

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

Redmine
  Bugs and issues, as well as some random features and discussions,
  are tracked at http://redmine.gromacs.org.

.. _Git Tips & Tricks: http://www.gromacs.org/index.php?title=Developer_Zone/Git/Git_Tips_%26_Tricks

Build system
------------

.. TODO: details, ASAN, others?

CMake
  Main tool used in the build system.

packaging for distribution (CPack)

unit testing (CTest)
  |Gromacs| uses a unit testing framework based on Google C++ Testing
  Framework (gtest) and CTest.  All unit tests are automatically run on Jenkins
  for each commit.
  Details can be found on a separate page on :doc:`testutils`.

regression tests

cppcheck
  `cppcheck <http://cppcheck.sourceforge.net>`_ is used for static code
  analysis, and is run automatically on Jenkins for each commit.  Different rules
  are used for C and C++ code (with stricter checking for C++ code, as that is
  newer).  The build system provides a ``cppcheck`` target (produced from
  ``tests/CppCheck.cmake``) to run the tool.  This target is used also by Jenkins.

clang static analyzer

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

Documentation generation
------------------------

Building the |Gromacs| documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. TODO: Move this onto a separate page

For now, there are multiple components, formats and tools for the
|Gromacs| documentation, which is aimed primarily at version-specific
deployment of the complete documentation on the website.

This is quite complex, because the dependencies for building the
documentation must not get in the way of building the code
(particularly when cross-compiling), and yet the code must build and
run in order for some documentation to be generated. Also, man page
documentation (and command-line completions) must be built from the
wrapper binary, in order to be bundled into the tarball.

The outputs of interest to most developers are generally produced in the
``docs/html/`` subdirectory of the build tree.

You need to enable at least some of the following CMake options:

``GMX_BUILD_MANUAL``
  Option needed for trying to build the PDF reference manual
  (requires LaTeX and ImageMagick)
``GMX_BUILD_HELP``
  Option that controls 1) whether shell completions are built automatically,
  and 2) whether built man pages are installed if available (the user still needs
  to build the ``man`` target manually before installing)

Some documentation cannot be built if the CMake option
``GMX_BUILD_MDRUN_ONLY`` is enabled, or when cross-compiling, as it
requires executing the ``gmx`` binary.

The following make targets are the most useful:

``manual``
  Builds the PDF reference manual
``man``
  Makes man pages from the wrapper binary with Sphinx
``doxygen-all``
  Makes the code documentation with Doxygen
``install-guide``
  Makes the INSTALL file for the tarball with Sphinx
``webpage-sphinx``
  Makes all the components of the GROMACS webpage that require Sphinx,
  including install guide and user guide.
``webpage``
  Makes the complete GROMACS webpage, requires everything. When complete,
  you can browse ``docs/html/index.html`` to find everything.

  If built from a release tarball, the ``SOURCE_MD5SUM``,
  ``SOURCE_TARBALL``, ``REGRESSIONTESTS_MD5SUM``, and
  ``REGRESSIONTESTS_TARBALL`` CMake variables can be set to pass in
  the md5sum values and names of those tarballs, for embedding into the
  final deployment to the |Gromacs| website.

The following tools are used in building parts of the documentation.

Doxygen
  `Doxygen <http://www.doxygen.org>`_ is used to extract documentation from
  source code comments.  Also some other overview
  content is laid out by Doxygen from Markdown source files.  Currently, version
  |EXPECTED_DOXYGEN_VERSION| is required for a warning-free build.  Thorough
  explanation of the Doxygen setup and instructions for documenting the source
  code can be found on a separate page: :doc:`doxygen`.

graphviz (dot)
  The Doxygen documentation uses ``dot`` from `graphviz
  <http://www.graphviz.org>`_ for building some graphs.  The tool is not
  mandatory, but the Doxygen build will produce warnings if it is not
  available, and the graphs are omitted from the documentation.

mscgen
  The Doxygen documentation uses `mscgen
  <http://www.mcternan.me.uk/mscgen/>`_ for building some graphs.  As with ``dot``,
  the tool is not mandatory, but not having it available will result in warnings
  and missing graphs.

Doxygen issue checker
  Doxygen produces warnings about some incorrect uses and wrong
  documentation, but there are many common mistakes that it does not detect.
  |Gromacs| uses an additional, custom Python script to check for such issues.
  This is most easily invoked through a ``check-source`` target in the build system.
  The script also checks that documentation for a header matches its use in the
  source code (e.g., that a header documented as internal to a module is not
  actually used from outside the module).  These checks are run in Jenkins as
  part of the Documentation job.  Details for the custom checker are on a
  separate page (common for several checkers): :doc:`gmxtree`.

module dependency graphs
  |Gromacs| uses a custom Python script to generate an annotated dependency
  graph for the code, showing #include dependencies between modules.
  The generated graph is embedded into the Doxygen documentation:
  `Module dependency graph`__
  This script shares most of its implementation with the custom checkers, and is
  documented on the same page: :doc:`gmxtree`.

__ doxygen-page-modulegraph_

Sphinx
  `Sphinx <http://sphinx-doc.org/>`_; at least version 1.2.3) is used
  for building some parts of the documentation from reStructuredText
  source files.

LaTeX
  Also requires ImageMagick for converting graphics file formats.

linkchecker

documentation exported from source files
  For man pages, HTML documentation of command-line options for executables,
  and for shell completions, the ``gmx`` binary has explicit C++ code to export
  the information required.  The build system provides targets that then invoke
  the built ``gmx`` binary to produce these documentation items.  The generated
  items are packaged into source tarballs so that this is not necessary when
  building from a source distribution (since in general, it will not work in
  cross-compilation scenarios).  To build and install these from a git
  distribution, explicit action is required.
  See `Doxygen documentation on the wrapper binary`__
  for some additional details.

__ doxygen-page-wrapperbinary_

.. include:: /fragments/doxygen-links.rst
