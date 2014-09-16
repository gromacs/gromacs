Source tree checker scripts {#page_dev_gmxtree}
===========================

There is a set of Python scripts, currently under `docs/doxygen/`, that check
various aspects of the source tree for consistency.  The script is based on
producing an abstract representation of the source tree from various sources:
 * List of files in the source tree (for overall layout of the source tree)
 * List of installed headers (extracted from the generated build system)
 * git attributes (to limit the scope of some checks)
 * Doxygen XML documentation:
   * For tags about public/private nature of documented headers and other
     constructs
   * For actual documented constructs, to check them for consistency.
 * Hard-coded knowledge about the \Gromacs source tree layout.

This representation is then used for various purposes:
 * Checking Doxygen documentation elements for common mistakes: missing brief
   descriptions, mismatches in file and class visibility, etc.
 * Checking for consistent usage and documentation of headers: e.g., a header
   that is documented as internal to a module should not be used outside that
   module.
 * Checking for module-level cyclic dependencies
 * Checking for consistent style and order of \#include directives
   (see \ref page_devincludes)
 * Actually sorting and reformatting \#include directives to adhere to the
   checked style
 * Generating dependency graphs between modules and for files within modules

The checks are run as part of a single `doc-check` target, but are described
in separate sections below.  In addition to printing the issues to `stderr`,
the script also writes them into `docs/doxygen/doxygen-check.log` for later
inspection.  Jenkins runs the checks as part of the Documentation job, and the
build is marked unstable if any issues are found.

For correct functionality, the scripts depend on correct usage of Doxygen
annotations described in \ref page_doxygen, in particular the visibility and
API definitions in file-level comments.

For some false positives from the script, the suppression mechanism described
below is the easiest way to silence the script, but otherwise the goal would be
to minimize the number of suppressions.

The scripts require Python 2.7 (other versions may work, but have not been
tested).

Checker details
===============

The `doc-check` target currently checks for a few different types of issues:
* For all Doxygen documentation (currently does not apply for members that do
  not appear in the documentation):
   * If a member has documentation, it should have a brief description.
   * A note is issued for in-body documentation for functions, since this is
     ignored by our current settings.
   * If a class has documentation, it should have public documentation only if
     it appears in an installed header.
   * If a class and its containing file has documentation, the class
     documentation should not be visible if the file documentation is not.
* For all files:
   * Consistent usage of

         #include "..." // This should be used for Gromacs headers

     and

         #include <...> // This should be used for system and external headers

   * Installed headers must not include non-installed headers.
   * If the file has a git attribute to identify it as a candidate for include
     sorting, the include sorter described below should not produce any
     changes.
* For documented files:
   * Installed headers should have public documentation, and other files should
     not.
   * The API level specified for a file should not be higher than where its
     documentation is visible.  For example, only publicly documented headers
     should be specified as part of the public API.
   * If an \c \\ingroup module_foo exists, it should match the subdirectory
     that the file is actually part of in the file system.
   * If a \c \\defgroup module_foo exists for the subdirectory where the file is,
     the file should contain \c \\ingroup module_foo.
   * Files should not include other files whose documentation visibility is
     lower (if the included file is not documented, the check is skipped).
* For files that are part of documented modules
  (\c \\defgroup module_foo exists for the subdirectory):
   * Such files should not be included from outside their module if they are
     undocumented or are not specified as part of library or public API.
* For all modules:
   * There should not be cyclic include dependencies between modules.

As a side effect, the XML extraction makes Doxygen parse all comments in the
code, even if they do not appear in the documentation.  This can reveal latent
issues in the comments, like invalid Doxygen syntax.  The messages from the XML
parsing are stored in `docs/doxygen/doxygen-xml.log` in the build tree, similar to
other Doxygen runs.

Suppressing issues
------------------

The script is not currently perfect (either because of unfinished
implementation, or because Doxygen bugs or incompleteness of the Doxygen XML
output), and the current code also contains issues that the script detects, but
the authors have not fixed.  To allow the script to still be used,
`doxygen/suppressions.txt` contains a list of issues that are filtered out from
the report.  The syntax is simple:

    <file>: <text>

where `<file>` is a path to the file that reports the message, and `<text>` is
the text reported.  Both support `*` as a wildcard.  If `<file>` is empty, the
suppression matches only messages that do not have an associated file.
`<file>` is matched against the trailing portion of the file name to make it
work even though the script reports absolute paths.
Empty lines and lines starting with `#` are ignored.

To add suppression for an issue, the line that reports the issue can be copied
into `suppressions.txt`, and the line number (if any) removed.  If the
issue does not have a file name (or a pseudo-file) associated, a leading `:`
must be added.  To cover many similar issues, parts of the line can then be
replaced with wildcards.

A separate suppression mechanism is in place for cyclic dependencies: to
suppress a cycle between moduleA and moduleB, add a line with format

    moduleA -> moduleB

into `doxygen/cycle-suppressions.txt`.  This suppresses all cycles that contain
the mentioned edge.  Since a cycle contains multiple edges, the suppression
should be made for the edge that is determined to be an incorrect dependency.
This also affects the layout of the include dependency graphs (see below): the
suppressed edge is not considered when determining the dependency order, and is
shown as invalid in the graph.

Include order checking/sorting
==============================

TODO

Include dependency graphs
=========================

The build system also provides a `dep-graphs` target that generates include
dependency graphs with some additional annotations.
One graph is produced that shows all the modules under `src/gromacs/`, and
their include dependencies.  Additionally, a file-level graph is produced for
each module, showing the include dependencies within that module.  Currently,
these are mostly for eye candy, but they can also be used for analyzing
problematic dependencies to clean up the architecture.
The output is put in `docs/doxygen/depgraphs/` in the build tree.

To get `.png` versions of the graphs, `graphviz` is
additionally required.

### Module graph ###

The graph is written into `module-deps.dot`, and embedded into the Doxygen
documentation: \ref page_modulegraph.  The embedded version contains a legend
explaining the graph.

### File graph ###

The graphs are written to <em>module_name</em>`-deps.dot.png`.

Node colors:
<dl>
<dt>light blue</dt>
<dd>public API (installed) headers</dd>
<dt>dark blue</dt>
<dd>library API headers</dd>
<dt>gray</dt>
<dd>source files</dd>
<dt>light green</dt>
<dd>test files</dd>
<dt>white</dt>
<dd>other files</dd>
</dl>

Each edge signifies an include dependency; there is no additional information
currently included.
