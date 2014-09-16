Guidelines for \#include directives {#page_devincludes}
===================================

The following include order is used in \Gromacs. An empty line should appear
between each group, and headers within each group sorted alphabetically.

 1. Each _source_ file should include `gmxpre.h` first.
 2. If a _source_ file has a corresponding header, it should be included next.
 3. If the file depends on defines from `config.h`, that comes next.
 4. This is followed by standard C/C++ headers, grouped as follows:
     a. Standard C headers (e.g., `<stdio.h>`)
     b. C++ versions of the above (e.g., `<cstdio>`)
     c. Standard C++ headers (e.g., `<vector>`)
    Preferably, only one of the first two groups is present, but this is not
    enforced.
 5. This is followed by other system headers: platform-specific headers such as
    `<unistd.h>`, as well as external libraries such as
    `<boost/scoped_ptr.hpp>`.
 6. \Gromacs-specific libraries from `src/external/`, such as
    `"thread_mpi/threads.h"`.
 7. \Gromacs-specific headers that are not internal to the including module,
    included with a path relative to `src/`.
 8. In _test_ files, headers not internal to the module, but specific to
    testing code, are in a separate block at this point.
 9. Finally, \Gromacs headers that are internal to the including module are
    included using a relative path (but never with a path starting with `../`;
    such headers go into group 7 instead).

Enforcing a consistent order and style has a few advantages:
 * It makes it easy at a quick glance to find the dependencies of a file,
   without scanning through a long list of unorganized \#includes.
 * Including the header corresponding to the source file first makes most
   headers included first in some source file, revealing potential problems
   where headers would not compile unless some other header would be included
   first.  With this order, the person working on the header is most likely to
   see these problems instead of someone else seeing them later when
   refactoring unrelated code.
 * Consistent usage of paths in \#include directives makes it easy to use
   `grep` to find all uses of a header, as well as all include dependencies
   between two modules.

The guidelines are enforced by an automatic checker script that can also
sort/reformat include statements to follow the guidelines.
See \ref page_dev_gmxtree.
