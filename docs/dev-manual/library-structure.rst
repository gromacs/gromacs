.. _library structure:

Library structure
=================

This section discusses the source tree layout where we want to aim in
the long term. Concrete steps of how this could be reached are
discussed in the next section.

* A single library, ``libgromacs``, is built from code under ``src/gromacs/``.
* Under src/gromacs/, we have subdirectories for modules. Each module
  consists of a set of routines that do some well-defined task or a
  collection of tasks, such as basic math functions, topology
  handling, or trajectory I/O.
* Modules must have a well-defined hierarchy: cycles are not allowed
  in the (include) dependency graph for the modules, even if they
  would not create cycles when looked at the file level (separately
  forbidden below).
* Each module directory contains one or more ``file.cpp`` files, each
  of which have a corresponding ``file.h`` file that declares the external
  API in that file. There can also be a ``file-impl.h`` file that declares
  classes or functions that are not accessible outside the
  file/module. In most cases, declarations that are not used outside
  the file should go into the source file, though. If a file does not
  contain any symbols that are accessible outside the module,
  the ``-impl`` suffix can be omitted, but comments in the header file
  should indicate that it should not be used outside the module.
* Each module has a header file ``src/gromacs/module.h`` that includes
  the public API headers from the module. This header is installed
  into ``include/gromacs/``, and all headers it includes are installed
  under ``include/gromacs/module/``. Other headers are not
  installed. Comments at the top of the header files should contain a
  note what's their visibility: public (installed), intra-library (can
  be used from inside the library), or intra-module/intra-file. See
  :ref:`using Doxygen` for documenting headers.
* Unit tests, and data required for them, go into a ``tests/``
  subdirectory under the module directory.
* When compiling, the include search path is set to ``src/``. Source
  files should include headers as ``#include
  "gromacs/module/file.h"``. ``#include "file.h"`` is allowed for
  intra-module/intra-file headers. Include files should include
  headers like ``#include "../othermodule/file.h"`` (relative paths
  are mandatory for installed headers, the source-file include syntax
  is allowed for other headers). This is enforced by :doc:`gmxtree`,
  which is run by Jenkins and votes accordingly in Gerrit.
* Every installed header should compile by itself, without any
  reference to variables defined in ``config.h`` or requiring other
  headers to be included before it. Cyclic include dependencies
  prevent this, and must be avoided cause this. This is best
  guaranteed by including every header in some source file as the
  first header, even before ``config.h``. This is partly enforced by
  :doc:`gmxtree`, which is run by Jenkins and votes accordingly in Gerrit.
* Code inside the library should not unnecessarily include headers. In
  particular, headers should not include other headers if a forward
  declaration of a type is enough for the header. Within the library
  source files, include only headers from other modules that are
  necessary for that file. You can use the public API header if you
  really require everything declared in it.
* Each executable is compiled from source files under
  ``src/programs/exename/``. The main function should reside in
  ``exename.cpp``, and other code specific to that executable should
  be in files organized as in the module directories. Executables
  should only use the public API headers, and not include anything
  from the module directories directly. Also, depending on the answers
  to http://redmine.gromacs.org/issues/988, the guidelines on usage of
  public API headers only may need to be revised.
* Some larger-scale tests for whole programs or analysis modules are
  implemented in ``src/programs/mdrun/tests`` and
  ``src/gromacs/gmxana/legacytests/``. Much of the test machinery
  currently present in the ``regressiontests`` repository should
  change into similar forms, but first we need more infrastructure for
  things like checking the results of mdrun (forces and energies).


