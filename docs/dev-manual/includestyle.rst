Guidelines for #include directives
==================================

The following include order is used in |Gromacs|. An empty line should appear
between each group, and headers within each group sorted alphabetically.

1. Each *source* file should include ``gmxpre.h`` first.
2. If a *source* file has a corresponding header, it should be included next.
   If the header is in the same directory as the source, then it is included
   without any path (i.e., relative to the source), otherwise relative to
   ``src/``.  The latter case is for headers in ``legacyheaders/`` and for tests.
3. If the file depends on defines from ``config.h``, that comes next.
4. This is followed by standard C/C++ headers, grouped as follows:

   1. Standard C headers (e.g., ``<stdio.h>``)
   2. C++ versions of the above (e.g., ``<cstdio>``)
   3. Standard C++ headers (e.g., ``<vector>``)

   Preferably, only one of the first two groups is present, but this is not
   enforced.
5. This is followed by other system headers: platform-specific headers such as
   ``<unistd.h>``, as well as external libraries such as
   ``<boost/scoped_ptr.hpp>``.
6. |Gromacs|-specific libraries from ``src/external/``, such as
   ``"thread_mpi/threads.h"``.
7. |Gromacs|-specific headers that are not internal to the including module,
   included with a path relative to ``src/``.
8. In *test* files, headers not internal to the module, but specific to
   testing code, are in a separate block at this point, paths relative to
   ``src/``.
9. Finally, |Gromacs| headers that are internal to the including module are
   included using a relative path (but never with a path starting with ``../``;
   such headers go into group 7 instead).  For test files, this group contains
   headers that are internal to tests for that module.

All |Gromacs| headers are included with quotes (``"gromacs/utility/path.h"``),
other headers with angle brackets (``<stdio.h>``).  Headers under ``src/external/``
are generally included with quotes (whenever the include path is relative to
``src/``, as well as for thread-MPI and TNG), but larger third-party entities are
included as if they were provided by the system.  The latter group currently
includes boost and gtest/gmock.

If there are any conditionally included headers (typically, only when some
#defines from ``config.h`` are set), these should be included at the end of
their respective group.  Note that the automatic checker/sorter script does not
act on such headers, nor on comments that are between #include statements; it
is up to the author of the code to put the headers in proper order in such
cases.  Trailing comments on the same line as #include statements are
preserved and do not affect the checker/sorter.

The guidelines are enforced by an automatic checker script that can also
sort/reformat include statements to follow the guidelines.
See :doc:`gmxtree` for details.

Enforcing a consistent order and style has a few advantages:

* It makes it easy at a quick glance to find the dependencies of a file,
  without scanning through a long list of unorganized #includes.
* Including the header corresponding to the source file first makes most
  headers included first in some source file, revealing potential problems
  where headers would not compile unless some other header would be included
  first.  With this order, the person working on the header is most likely to
  see these problems instead of someone else seeing them later when
  refactoring unrelated code.
* Consistent usage of paths in #include directives makes it easy to use
  ``grep`` to find all uses of a header, as well as all include dependencies
  between two modules.
* An automatic script can be used to re-establish clean code after
  semi-automatic refactoring like renaming an include file with ``sed``, without
  causing other unnecessary changes.
