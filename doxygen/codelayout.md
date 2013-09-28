General organization of the code and documentation {#page_codelayout}
==================================================

Source code organization
------------------------

The source code for \Gromacs is under the `src/` subdirectory
(except for an analysis tool template, which is under `share/template/`).
The subfolders under this directory are:
<dl>
<dt>`src/gromacs/`</dt>
<dd>
The code under this directory is built into a single library,
`libgromacs`.  Installed headers are also located in this hierarchy.
This is the main part of the code, and is organized into further subdirectories
as _modules_ (see below)
</dd>
<dt>`src/programs/`</dt>
<dd>
\Gromacs executables are built from code under this directory.
Although some build options can change this, there is typically only a single
binary, `gmx`, built.
</dd>
<dt>`src/testutils/`</dt>
<dd>
Shared utility code for writing unit tests is found under this directory.
</dd>
<dt>`src/contrib/`</dt>
<dd>
This directory contains collection of less well maintained code that may or may
not compile.  It is not included in the build.
</dd>
<dt>`src/external/`</dt>
<dd>
This directory contains bundled source code for various libraries and
components that \Gromacs uses internally.
</dd>
</dl>

Organization under `src/gromacs/`
---------------------------------

There is no code directly under `src/gromacs/`, except for some public API
convenience headers.  The code is organized into subdirectories, denoted
_modules_.  Each module consists of a set of routines that do some well-defined
task or a collection of tasks.  Installed headers are a subset of the headers
in `src/gromacs/` and in the module subdirectories.  They are installed into a
corresponding hierarchy under `include/gromacs/` in the installation directory.
Comments at the top of the header files contain a note about their visibility:
public (installed), intra-library (can be used from inside the library), or
intra-module/intra-file.

Each module directory contains one or more `<file>.c/.cpp` files, each of which
has a corresponding `<file>.h` file that declares the external API in that file
(there are some exceptions, but this gives a general picture).
Typically, a C++ file declares a class of the same or similar name, but may
also declare some related classes.
There can also be a `<file>-impl.h` file that declares classes or functions that
are not accessible outside the module.  In most cases, declarations that
are not used outside a single source file are in the source file.

Unit tests, and data required for them, are in a `tests/` subdirectory under
the module directory.

When compiling, the include search path is set to `src/`.  This means that
source files include headers as

    #include "gromacs/<module>/<file>.h"

The following is also possible for intra-module headers:

    #include "<file>.h"

Header files include other headers using

    #include "../<othermodule>/<file>.h"

because relative paths work best for installed headers.  For non-installed
headers, the path relative to `src/` is sometimes also used.

For historical reasons, there are directories `src/gromacs/gmxana/`,
`src/gromacs/gmxlib/`, `src/gromacs/mdlib/`, and `src/gromacs/gmxpreprocess/`
that do not follow the above rules.  The installed headers for these are in
`src/gromacs/legacyheaders/`.  The aim is to gradually get rid of these
directories and move code into proper modules.

For similar historical reasons, the include path also includes
`src/gromacs/legacyheaders/` and `src/gromacs/gmxpreprocess/` (the latter only
for part of the source).  New code should not depend on these.
