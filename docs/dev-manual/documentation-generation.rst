Documentation generation
========================

Building the |Gromacs| documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For now, there are multiple components, formats and tools for the
|Gromacs| documentation, which is aimed primarily at version-specific
deployment of the complete documentation on the website and in the
release tarball.

This is quite complex, because the dependencies for building the
documentation must not get in the way of building the code
(particularly when cross-compiling), and yet the code must build and
run in order for some documentation to be generated. Also, man page
documentation (and command-line completions) must be built from the
wrapper binary, in order to be bundled into the tarball. This helps
ensure that the functionality and the documentation remain in sync.

The outputs of interest to most developers are generally produced in the
``docs/html/`` subdirectory of the build tree.

You need to enable at least some of the following CMake options:

``GMX_BUILD_MANUAL``
  Option needed for trying to build the PDF reference manual
  (requires LaTeX and ImageMagick). See :cmake:`GMX_BUILD_MANUAL`.
``GMX_BUILD_HELP``
  Option that controls 1) whether shell completions are built automatically,
  and 2) whether built man pages are installed if available (the user still needs
  to build the ``man`` target manually before installing). See
  :cmake:`GMX_BUILD_HELP`.

Some documentation cannot be built if the CMake option
``GMX_BUILD_MDRUN_ONLY`` is enabled, or when cross-compiling, as it
requires executing the ``gmx`` binary.

The following make targets are the most useful:

``manual``
  Builds the PDF reference manual.
``man``
  Makes man pages from the wrapper binary with Sphinx.
``doxygen-all``
  Makes the code documentation with Doxygen.
``install-guide``
  Makes the INSTALL file for the tarball with Sphinx.
``webpage-sphinx``
  Makes all the components of the |Gromacs| webpage that require Sphinx,
  including install guide and user guide.
``webpage``
  Makes the complete |Gromacs| webpage, requires everything. When complete,
  you can browse ``docs/html/index.html`` to find everything.

  If built from a release tarball, the ``SOURCE_MD5SUM``,
  ``SOURCE_TARBALL``, ``REGRESSIONTESTS_MD5SUM``, and
  ``REGRESSIONTESTS_TARBALL`` CMake variables can be set to pass in
  the md5sum values and names of those tarballs, for embedding into the
  final deployment to the |Gromacs| website.

Needed build tools
^^^^^^^^^^^^^^^^^^

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
  `Sphinx <http://sphinx-doc.org/>`_; at least version |EXPECTED_SPHINX_VERSION| is used
  for building some parts of the documentation from reStructuredText
  source files.

LaTeX
  Also requires ImageMagick for converting graphics file formats.

linkchecker
  The linkchecker program is used together with the linkcheckerrc file to ensure
  that all the links in the documentation can be resolved correctly.

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
