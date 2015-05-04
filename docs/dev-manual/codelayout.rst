General organization of the code and documentation
==================================================

Source code organization
------------------------

The source code for |Gromacs| is under the ``src/`` subdirectory
(except for an analysis tool template, which is under ``share/template/``).
The subfolders under this directory are:

``src/gromacs/``
  The code under this directory is built into a single library,
  `libgromacs`.  Installed headers are also located in this hierarchy.
  This is the main part of the code, and is organized into further subdirectories
  as *modules* (see below).
``src/programs/``
  |Gromacs| executables are built from code under this directory.
  Although some build options can change this, there is typically only a single
  binary, ``gmx``, built.

``src/testutils/``
  Shared utility code for writing unit tests is found under this directory.
``src/external/``
  This directory contains bundled source code for various libraries and
  components that |Gromacs| uses internally.
``src/contrib/``
  This directory contains collection of less well maintained code that may or may
  not compile.  It is not included in the build.

Organization under ``src/gromacs/``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is no code directly under ``src/gromacs/``, except for some public API
convenience headers.  The code is organized into subdirectories, denoted
*modules*.  Each module consists of a set of routines that do some well-defined
task or a collection of tasks.  Installed headers are a subset of the headers
in ``src/gromacs/`` and in the module subdirectories.  They are installed into a
corresponding hierarchy under ``include/gromacs/`` in the installation directory.
Comments at the top of the header files contain a note about their visibility:
public (installed), intra-library (can be used from inside the library), or
intra-module/intra-file.

Each module directory contains one or more :file:`{file}.c/.cpp` files, each of which
has a corresponding :file:`{file}.h` file that declares the external API in that file
(there are some exceptions, but this gives a general picture).
Typically, a C++ file declares a class of the same or similar name, but may
also declare some related classes.
There can also be a :file:`{file}-impl.h` file that declares classes or functions that
are not accessible outside the module.  In most cases, declarations that
are not used outside a single source file are in the source file.

Unit tests, and data required for them, are in a ``tests/`` subdirectory under
the module directory.
See :doc:`testutils` for more details.

When compiling, the include search path is set to ``src/``.  This means that
files include headers as ::

    #include "gromacs/<module>/<file>.h"

The following is also possible for intra-module headers::

    #include "<file>.h"

See :doc:`includestyle` for more details.

For historical reasons, there are directories ``src/gromacs/gmxana/``,
``src/gromacs/gmxlib/`, ``src/gromacs/mdlib/``, and ``src/gromacs/gmxpreprocess/``
that do not follow the above rules.  The installed headers for these are in
`src/gromacs/legacyheaders/``.  The aim is to gradually get rid of these
directories and move code into proper modules.

.. _dev-doc-layout:

Documentation organization
--------------------------

This Doxygen documentation is made of a few different parts.  Use the list
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
