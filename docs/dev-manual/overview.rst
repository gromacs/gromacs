Codebase overview
=================

The root directory of the |Gromacs| repository only contains :file:`CMakeLists.txt`
(the root file for the CMake build system), a few files supporting the build
system, and a few standard informative files (:file:`README` etc.).  The
:file:`INSTALL` is generated for source packages from
:file:`docs/install-guide/index.rst`.

All other content is in the following top-level directories:

:file:`admin/`
  Contains various scripts for developer use, as well as configuration files
  and scripts for some of the tools used.
:file:`cmake/`
  Contains code fragments and find modules for CMake.
  Some content here is copied and/or adapted from newer versions of CMake than
  the minimum currently supported.
  Default suppression file for valgrind is also included here.
  See :doc:`build-system` for details of the build system.
:file:`docs/`
  Contains the build system logic and source code for all documentation, both
  user-facing and developer-facing.  Some of the documentation is generated
  from the source code under :file:`src/`; see :ref:`dev-doc-layout`.
  This directory also contains some developer scripts that use the Doxygen
  documentation for their operation.
:file:`scripts/`
  Contains the templates for :file:`GMXRC` script, some other installed scripts,
  as well as installation rules for all these scripts.
:file:`share/`
  Contains data files that will be installed under :file:`share/`.  These
  include a template for writing C++ analysis tools, and data files used by
  |Gromacs|.
:file:`src/`
  Contains all source code.  See :ref:`dev-source-layout`.
:file:`tests/`
  Contains build system logic for some high-level tests.  Currently, only the
  regression test build system logic and cppcheck rules are in this directory,
  while other tests are under :file:`src/`.

.. _dev-source-layout:

Source code organization
------------------------

The following figure shows a high-level view of components of what gets built
from the source code under :file:`src/` and how the code is organized.
The build system is described in detail in :doc:`build-system`.
With default options, the green and white components are built as part of the
default target.  If ``GMX_BUILD_MDRUN_ONLY`` is ``ON``, then the blue and white
components are built instead; :file:`libgromacs_mdrun` is built from a subset
of the code used for :file:`libgromacs`.
The gray parts are for testing, and are by default only built as part of the
``tests`` target, but if ``GMX_DEVELOPER_BUILD`` is ``ON``, then these are
included in the default build target.
See :doc:`testutils` for details of the testing side.

.. digraph:: dev_high_level_components

   concentrate = yes
   node [ shape=box, style=filled, width=2 ]

   subgraph {
     rank = same
     externals [
       label="externals\nsrc/external/", group=common, style=rounded
     ]
     gtest [
       label="Google Test & Mock\nsrc/external/gmock-1.7.0/", group=test
       style="rounded,filled", fillcolor="0 0 0.9"
     ]
   }
   subgraph {
     rank = same
     libgromacs [
       label="libgromacs\nsrc/gromacs/", group=gmx, fillcolor="0.33 0.3 1"
     ]
     libgromacs_mdrun [
       label="libgromacs_mdrun\nsrc/gromacs/", group=mdrun, fillcolor="0.66 0.3 1"
     ]
   }
   testutils [
     label="testutils\nsrc/testutils/", group=test
     style="rounded,filled", fillcolor="0 0 0.9"
   ]
   mdrun_objlib [
     label="mdrun object lib.\nsrc/programs/mdrun/", group=common, style=rouded
   ]
   subgraph {
     rank = same
     gmx [
       label="gmx\nsrc/programs/", group=gmx, fillcolor="0.33 0.3 1"
     ]
     mdrun [
       label="mdrun\nsrc/programs/", group=mdrun, fillcolor="0.66 0.3 1"
     ]
     tests [
       label="test binaries\nsrc/.../tests/", group=test
       style="rounded,filled", fillcolor="0 0 0.9"
     ]
     template [
       label="analysis template\nshare/template/", group=common
       fillcolor="0.33 0.3 1"
     ]

     gmx -> template [ style=invis, constraint=no ]
     template -> mdrun [ style=invis, constraint=no ]
   }

   libgromacs -> externals
   libgromacs_mdrun -> externals
   mdrun_objlib -> libgromacs
   gmx -> libgromacs
   gmx -> mdrun_objlib
   mdrun -> libgromacs_mdrun
   mdrun -> mdrun_objlib
   testutils -> externals
   testutils -> gtest
   testutils -> libgromacs
   tests -> gtest
   tests -> libgromacs
   tests -> mdrun_objlib
   tests -> testutils
   template -> libgromacs

   template -> mdrun_objlib [ style=invis ]
   mdrun_objlib -> externals [ style=invis ]

All the source code (except for the analysis template) is under the
:file:`src/` directory.  Only a few files related to the build system are
included at the root level.  All actual code is in subdirectories:

:file:`src/gromacs/`
  The code under this directory is built into a single library,
  :file:`libgromacs`.  Installed headers are also located in this hierarchy.
  This is the main part of the code, and is organized into further subdirectories
  as *modules*.  See below for details.
:file:`src/programs/`
  |Gromacs| executables are built from code under this directory.
  Although some build options can change this, there is typically only a single
  binary, :file:`gmx`, built.

:file:`src/{...}/tests/`
  Various subdirectories under :file:`src/` contain a subdirectory named
  :file:`tests/`.  The code from each such directory is built into a test
  binary.  Some such directories also provide shared test code as object
  libraries that is linked into multiple test binaries from different folders.
  See :doc:`testutils` for details.
:file:`src/testutils/`
  Contains shared utility code for writing Google Test tests.
  See :doc:`testutils` for details.
:file:`src/external/`
  Contains bundled source code for various libraries and
  components that |Gromacs| uses internally.  All the code from these
  directories are built using our custom build rules into :file:`libgromacs`,
  or in some cases into the test binaries.  Some CMake options change which
  parts of this code are included in the build.
  See :doc:`build-system` for some explanation about how the code in this
  directory is used.
:file:`src/contrib/`
  Contains collection of less well maintained code that may or may
  not compile.  It is not included in the build.
:file:`src/contrib/fftw/`
  As an exception to the above, this folder contains the build system code for
  downloading and building FFTW to be included into :file:`libgromacs`.

When compiling, the include search path is set to :file:`src/`.
Some directories from under :file:`src/external/` may also be included,
depending on the compilation options.

Organization under :file:`src/gromacs/`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The :file:`libgromacs` library is built from code under :file:`src/gromacs/`.
Again, the top-level directory contains build and installation rules for the
library, and :dfn:`public API convenience headers`.  These convenience headers
provide the main installed headers that other code can use.  They do not
contain any declarations, but only include a suitable set of headers from the
subdirectories.  They typically also contain high-level Doxygen documentation
for the subdirectory with the same name: :file:`{module}.h` corresponds to
:file:`{module}/`.

The code is organized into subdirectories.  These subdirectories are denoted as
:dfn:`modules` throughout this documentation.  Each module consists of a set
of routines that do some well-defined task or a collection of tasks.

Installed headers are a subset of the headers under :file:`src/gromacs/`.
They are installed into a corresponding hierarchy under
:file:`include/gromacs/` in the installation directory.
Comments at the top of the header files contain a note about their visibility:
public (installed), intra-library (can be used from inside the library), or
intra-module/intra-file.

See :doc:`naming` for some common naming patterns for files that can help
locating declarations.

Tests, and data required for them, are in a :file:`tests/` subdirectory under
the module directory.
See :doc:`testutils` for more details.

For historical reasons, there are directories :file:`src/gromacs/gmxana/`,
:file:`src/gromacs/gmxlib/`, :file:`src/gromacs/mdlib/`, and
:file:`src/gromacs/gmxpreprocess/` that do not follow the above rules.
The installed headers for these are in :file:`src/gromacs/legacyheaders/`.
The aim is to gradually get rid of these directories and move code into proper
modules.

.. _dev-doc-layout:

Documentation organization
--------------------------

All documentation (including this developer guide) is produced from source
files under :file:`docs/`, except for some command-line help that is generated
from the source code (by executing the compiled :file:`gmx` binary).
The build system provides various custom targets that build the documentation;
see :doc:`build-system` for details.

:file:`docs/fragments/`
  Contains reStructuredText fragments used through ``.. include::`` mechanism
  from various places in the documentation.

User documentation
^^^^^^^^^^^^^^^^^^

:file:`docs/install-guide/`
  Contains reStructuredText source files for building the install guide section
  of the user documentation, as well as the :file:`INSTALL` file for the source
  package.
  The build rules are in :file:`docs/CMakeLists.txt`.
:file:`docs/manual/`
  Contains LaTeX source files and build rules for the reference (PDF) manual.
:file:`docs/user-guide/`
  Contains reStructuredText source files for building the user guide section
  of the user documentation.
  The build rules are in :file:`docs/CMakeLists.txt`.

Unix man pages
^^^^^^^^^^^^^^

Man pages for programs are generated by running the :file:`gmx` executable
after compiling it, and then using Sphinx on the reStructuredText files that
:file:`gmx` writes out.

The build rules for the man pages are in :file:`docs/CMakeLists.txt`.

Developer guide
^^^^^^^^^^^^^^^

:file:`docs/dev-manual/`
  Contains reStructuredText source files used to build the developer guide.
  The build rules are in :file:`docs/CMakeLists.txt`.

The organization of the developer guide is explained on the :doc:`front page of
the guide <index>`.

Doxygen documentation
^^^^^^^^^^^^^^^^^^^^^

:file:`docs/doxygen/`
  Contains the build rules and some overview content for the Doxygen
  documentation.
  See :doc:`doxygen` for details of how the Doxygen documentation is built and
  organized.

.. TODO: Create a separate page (at the front of the developer guide, and/or at
   the main index.rst) that describes the documentation from readers'
   perspective, and move relevant content there.  This should contain just an
   overview of how the documentation is organized in the source tree.

The Doxygen documentation is made of a few different parts.  Use the list
below as a guideline on where to look for a particular kind of content.
Since the documentation has been written over a long period of time and the
approach has evolved, not all the documentation yet follows these guidelines,
but this is where we are aiming at.

documentation pages
  These contain mainly overview content, from general-level introduction down
  into explanation of some particular areas of individual modules.
  These are generally the place to start familiarizing with the code or a new
  area of the code.
  They can be reached by links from the main page, and also through cross-links
  from places in the documentation where that information is relevant to
  understand the context.
module documentation
  These contain mainly techical content, explaining the general implementation of
  a particular module and listing the classes, functions etc. in the module.
  They complement pages that describe the concepts.
  They can be reached from the Modules tab, and also from all individual classes,
  functions etc. that make up the module.
class documentation
  These document the usage of an individual class, and in some cases that of
  closely related classes.  Where necessary (and time allowing), a broader
  overview is given on a separate page and/or in the module documentation.
method documentation
  These document the individual method.  Typically, the class documentation or
  other overview content is the place to look for how different methods interact.
file and namespace documentation
  These are generally only placeholders for links, and do not contain much else.
  The main content is the list of classes and other entities declared in that
  file.
