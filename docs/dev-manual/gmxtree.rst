Source tree checker scripts
===========================

.. highlight:: bash

There is a set of Python scripts, currently under ``docs/doxygen/``, that check
various aspects of the source tree for consistency.  The script is based on
producing an abstract representation of the source tree from various sources:

* List of files in the source tree (for overall layout of the source tree)
* List of installed headers (extracted from the generated build system)
* git attributes (to limit the scope of some checks)
* Doxygen XML documentation:

  * For tags about public/private nature of documented headers and other
    constructs
  * For actual documented constructs, to check them for consistency

* Hard-coded knowledge about the |Gromacs| source tree layout

This representation is then used for various purposes:

* Checking Doxygen documentation elements for common mistakes: missing brief
  descriptions, mismatches in file and class visibility, etc.
* Checking for consistent usage and documentation of headers: e.g., a header
  that is documented as internal to a module should not be used outside that
  module.
* Checking for module-level cyclic dependencies
* Checking for consistent style and order of #include directives
  (see :doc:`includestyle`)
* Actually sorting and reformatting #include directives to adhere to the
  checked style
* Generating dependency graphs between modules and for files within modules

The checks are run as part of a single ``check-source`` target, but are described
in separate sections below.  In addition to printing the issues to ``stderr``,
the script also writes them into ``docs/doxygen/check-source.log`` for later
inspection.  Jenkins runs the checks as part of the Documentation job, and the
build is marked unstable if any issues are found.

For correct functionality, the scripts depend on correct usage of Doxygen
annotations described in :doc:`doxygen`, in particular the visibility and
API definitions in file-level comments.

For some false positives from the script, the suppression mechanism described
below is the easiest way to silence the script, but otherwise the goal would be
to minimize the number of suppressions.

The scripts require Python 2.7 (other versions may work, but have not been
tested).

To understand how the scripts work internally, see comments in the Python
source files under ``docs/doxygen/``.

Checker details
---------------

The ``check-source`` target currently checks for a few different types of issues.
These are listed in detail below, mainly related to documentation and include
dependencies.  Note in particular that the include dependency checks are much
stricter for code in modules/directories that are documented with a
``\defgroup``: all undocumented code is assumed to be internal to such modules.
The rationale is that such code has gotten some more attention, and some effort
should also have been put into defining what is the external interface of the
module and documenting it.

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

  * Consistent usage of ::

        #include "..." // This should be used for GROMACS headers

    and ::

        #include <...> // This should be used for system and external headers

  * Installed headers must not include non-installed headers.
  * All source files must include "gmxpre.h" as the first header.
  * A source/header file should include "config.h," "gromacs/simd/simd.h",
    or "gromacs/ewald/pme-simd.h" if and only if it uses a macro declared
    in such files.
  * If the file has a git attribute to identify it as a candidate for include
    sorting, the include sorter described below should not produce any
    changes (i.e., the file should follow :doc:`includestyle`).

* For documented files:

  * Installed headers should have public documentation, and other files should
    not.
  * The API level specified for a file should not be higher than where its
    documentation is visible.  For example, only publicly documented headers
    should be specified as part of the public API.
  * If an ``\ingroup module_foo`` exists, it should match the subdirectory
    that the file is actually part of in the file system.
  * If a ``\defgroup module_foo`` exists for the subdirectory where the file
    is, the file should contain ``\ingroup module_foo``.
  * Files should not include other files whose documentation visibility is
    lower (if the included file is not documented, the check is skipped).

* For files that are part of documented modules
  (``\defgroup module_foo`` exists for the subdirectory), or are explicitly
  documented to be internal or in the library API:

  * Such files should not be included from outside their module if they are
    undocumented (for documented modules) or are not specified as part of
    library or public API.

* For all modules:

  * There should not be cyclic include dependencies between modules.

As a side effect, the XML extraction makes Doxygen parse all comments in the
code, even if they do not appear in the documentation.  This can reveal latent
issues in the comments, like invalid Doxygen syntax.  The messages from the XML
parsing are stored in ``docs/doxygen/doxygen-xml.log`` in the build tree, similar to
other Doxygen runs.

Suppressing issues
^^^^^^^^^^^^^^^^^^

The script is not currently perfect (either because of unfinished
implementation, or because Doxygen bugs or incompleteness of the Doxygen XML
output), and the current code also contains issues that the script detects, but
the authors have not fixed.  To allow the script to still be used,
``doxygen/suppressions.txt`` contains a list of issues that are filtered out from
the report.  The syntax is simple::

    <file>: <text>

where ``<file>`` is a path to the file that reports the message, and ``<text>`` is
the text reported.  Both support ``*`` as a wildcard.  If ``<file>`` is empty, the
suppression matches only messages that do not have an associated file.
``<file>`` is matched against the trailing portion of the file name to make it
work even though the script reports absolute paths.
Empty lines and lines starting with ``#`` are ignored.

To add a suppression for an issue, the line that reports the issue can be
copied into ``suppressions.txt``, and the line number (if any) removed.  If the
issue does not have a file name (or a pseudo-file) associated, a leading ``:``
must be added.  To cover many similar issues, parts of the line can then be
replaced with wildcards.

A separate suppression mechanism is in place for cyclic dependencies: to
suppress a cycle between moduleA and moduleB, add a line with format ::

    moduleA -> moduleB

into ``doxygen/cycle-suppressions.txt``.  This suppresses all cycles that contain
the mentioned edge.  Since a cycle contains multiple edges, the suppression
should be made for the edge that is determined to be an incorrect dependency.
This also affects the layout of the include dependency graphs (see below): the
suppressed edge is not considered when determining the dependency order, and is
shown as invalid in the graph.

.. _dev-include-sorter:

Include order sorting
---------------------

The script checks include ordering according to :doc:`includestyle`.
If it is not obvious how the includes should be changed to make the script
happy, or bulk changes are needed in multiple files, e.g., because of a header
rename or making a previously public header private, it is possible to run a
Python script that does the sorting::

    docs/doxygen/includesorter.py -S . -B ../build <files>

The script needs to know the location of the source tree (given with ``-S``) and
the build tree (given with ``-B``), and sorts the given files.  To sort the whole
source tree, one can also use::

    admin/reformat_all.sh includesort -B=../build

For the sorter to work correctly, the build tree should contain up-to-date list
of installed files and Doxygen XML documentation.  The former is created
automatically when ``cmake`` is run, and the latter can be built using the
``doxygen-xml`` target.

Note that currently, the sorter script does not change between angle brackets
and quotes in include statements.

Include dependency graphs
-------------------------

The same set of Python scripts can also produce include dependency graphs with
some additional annotations compared to what, e.g., Doxygen produces for a
directory dependency graph.  Currently, a module-level graph is automatically
built when the Doxygen documentation is built and embedded in the documentation
(not in the public API documentation).  The graph, together with a legend, is
on a separate page: `Module dependency graph`__

__ doxygen-page-modulegraph_

The Python script produces the graphs in a format suitable for ``dot`` (from the
``graphviz`` package) to lay them out.  The build system also provides a
``dep-graphs`` target that generates PNG files from the intermediate ``dot`` files.
In addition to the module-level graph, a file-level graph is produced for
each module, showing the include dependencies within that module.
The file-level graphs can only be viewed as the PNG files, with some
explanation of the notation below.  Currently, these are mostly for eye candy,
but they can also be used for analyzing problematic dependencies to clean up
the architecture.

Both the intermediate ``.dot`` files and the final PNG files are put under
``docs/doxygen/depgraphs/`` in the build tree.

File graphs
^^^^^^^^^^^

The graphs are written to :file:`{module_name}-deps.dot.png`.

Node colors:

light blue
  public API (installed) headers
dark blue
  library API headers
gray
  source files
light green
  test files
white
  other files

Each edge signifies an include dependency; there is no additional information
currently included.

.. include:: /fragments/doxygen-links.rst
