Guidelines for #include directives
==================================

The following include order is used in |Gromacs| and enforced by ``clang-format``.
An empty line should appear between each group, and headers within 
each group sorted alphabetically.

1. Each *source* file should include ``gmxpre.h`` first.
2. If a *source* file has a corresponding header, it should be included next.
   If the header is in the same directory as the source, then it is included
   without any path (i.e., relative to the source). Otherwise, the canonical
   include path of ``libraryname/modulename/header.h`` is used.
3. If the file depends on defines from ``config.h``, that comes next.
4. This is followed by standard C/C++ headers, grouped as follows:

   1. Standard C headers (e.g., ``<stdio.h>``)
   2. C++ versions of the above (e.g., ``<cstdio>``)
   3. Standard C++ headers (e.g., ``<vector>``)

   Preferably, only one of the first two groups is present, but this is not
   enforced.
5. This is followed by other system headers: platform-specific headers such as
   ``<unistd.h>``, as well as external libraries such as
   ``<gtest/gtest.h>``.
6. |Gromacs|-specific libraries from ``src/external/``, such as
   ``"thread_mpi/threads.h"``.
7. |Gromacs| headers that are not part of the including module.
8. Public |Gromacs| headers that are part of the including module.
9. Finally, |Gromacs| headers that are internal to the including module,
   executable, or test target
   (typically at the same path as the source file).

All |Gromacs| headers are included with quotes (``"gromacs/utility/path.h"``),
other headers with angle brackets (``<stdio.h>``).  Headers under ``src/external/``
are generally included with quotes (whenever the include path is relative to
``src/``, as well as for thread-MPI and TNG), but larger third-party entities are
included as if they were provided by the system.  The latter group currently
includes gtest/gmock.

In some cases, the include paths available to build targets may leak visibility
of headers inappropriately. This is usually encountered as a header that can be
used by an ``#include`` with an unusual or long path. If a header cannot be
included as described above, check that the appropriate CMake target is referenced by a
`target_link_libraries() <https://cmake.org/cmake/help/latest/command/target_link_libraries.html>`__
command. Many modules provide their own CMake target. Additionally, note

* The ``common`` CMake target provides access to :file:`gmxpre.h`, :file:`config.h`,
  :file:`gmxpre-config.h`, :file:`buildinfo.h`, and :file:`contributors.h`
* ``legacy_api`` provides access to those of the old :file:`gromacs/modulename` headers
  that are in :file:`api/legacy/include`
* ``legacy_modules`` adds :file:`src/` to the include path, exposing all headers in
  :file:`gromacs/` and :file:`gromacs/*/` for ``#include`` lines that would appear
  to comply with the guidelines above, but which may not be intended for "public" use.
  (This target was intended as a temporary measure while working towards :issue:`3288`.)

If there are any conditionally included headers (typically, only when some
#defines from ``config.h`` are set), these should be included at the end of
their respective group.  Note that the automatic checker/sorter script does not
act on such headers, nor on comments that are between #include statements; it
is up to the author of the code to put the headers in proper order in such
cases.  Trailing comments on the same line as #include statements are
preserved and do not affect the checker/sorter.

As part of the effort to build a proper API, a new scheme of separating between
public, library and module functionality in header files is planned.
See also :doc:`gmxtree` and
`API restructuring issues <https://gitlab.com/gromacs/gromacs/-/issues?label_name%5B%5D=API+restructuring>`__
for details.

Enforcing a consistent order and style has a few advantages:

* It makes it easy at a quick glance to find the dependencies of a file,
  without scanning through a long list of unorganized #includes.
* Including the header corresponding to the source file first makes most
  headers included first in some source file, revealing potential problems
  where headers would not compile unless some other header would be included
  first.  With this order, the person working on the header is most likely to
  see these problems instead of someone else seeing them later when
  refactoring unrelated code.
* Consistent usage of paths in ``#include`` directives makes it easy to use
  ``grep`` to find all uses of a header, as well as all include dependencies
  between two modules.
* An automatic script can be used to re-establish clean code after
  semi-automatic refactoring like renaming an include file with ``sed``, without
  causing other unnecessary changes.
