Development-time tools {#page_devtools}
======================

TODO: Add details for most of the tools, either in the form of links to wiki,
or to a separate page that explains more details.

Change management
=================

<dl>

<dt>git</dt>
<dd>\Gromacs uses git as the version control system.
Instructions for setting up git for \Gromacs, as well as tips and tricks for
its use, can be found on the wiki:
[Git Tips & Tricks](http://www.gromacs.org/index.php?title=Developer_Zone/Git/Git_Tips_%26_Tricks) <br/>
Other basic tutorial material for `git` can be found on the web.</dd>

<dt>Gerrit</dt>
<dd>All code changes go through a code review system at
<http://gerrit.gromacs.org>.</dd>

<dt>Jenkins</dt>
<dd>All changes pushed to Gerrit are automatically compiled and otherwise
checked on various platforms using a continuous integration system at
<http://jenkins.gromacs.org>.</dd>

<dt>Redmine</dt>
<dd>Bugs and issues, as well as some random features and discussions,
are tracked at <http://redmine.gromacs.org>.</dd>

</dl>

Build system
============

TODO: details, ASAN, others?

<dl>

<dt>CMake</dt>
<dd></dd>

<dt>packaging for distribution (CPack)</dt>
<dd></dd>

<dt>unit testing (CTest)</dt>
<dd>\Gromacs uses a unit testing framework based on Google C++ Testing
Framework (gtest) and CTest.  All unit tests are automatically run on Jenkins
for each commit.
Details can be found on a separate page: \subpage page_unittesting.</dd>

<dt>regression tests</dt>
<dd></dd>

<dt>cppcheck</dt>
<dd>cppcheck (<http://cppcheck.sourceforge.net>) is used for static code
analysis, and is run automatically on Jenkins for each commit.  Different rules
are used for C and C++ code (with stricter checking for C++ code, as that is
newer).  The build system provides a `cppcheck` target (produced from
`tests/CppCheck.cmake`) to run the tool.  This target is used also by Jenkins.
</dd>

<dt>clang static analyzer</dt>
<dd></dd>

</dl>

Code formatting and style
=========================

The tools and scripts listed below are used to automatically check/apply
formatting that follows \ref page_devstyle.

<dl>

<dt>uncrustify</dt>
<dd>Uncrustify (<http://uncrustify.sourceforge.net>) is used for automatic
indentation and other formatting of the source code.  All code must remain
invariant under uncrustify with the config at `admin/uncrustify.cfg`.
A patched version of uncrustify is used.</dd>

<dt>`admin/copyright.py`</dt>
<dd>Custom Python script is used to automatically add and format copyright
headers to source files.  `uncrustify.sh` (see below) uses the script to update
copyright years on changed files automatically.</dd>

<dt>`admin/uncrustify.sh`</dt>
<dd>`bash` script is provided to run uncrustify and `copyright.py` for all
files that have local changes and check that they conform to the prescribed
style.  Optionally, the script can also apply changes to make the files
conform.
This script is automatically run by Jenkins to ensure that all commits adhere
to the formatting guidelines.  If the uncrustify job does not succeed, it means
that this script has something to complain.</dd>

<dt>git pre-commit hook</dt>
<dd>git pre-commit hook is provided at `admin/git-pre-commit` if one wants to
apply `uncrustify.sh` automatically before every commit to check for formatting
issues.</dd>

<dt>include directive checker</dt>
<dd>Custom Python script is provided to check that code follows the guidelines
from \ref page_devstyle_includes.  These checks run as part of the `doc-check`
target provided by the build system.  Details for the custom checker are
on a separate page (common for several checkers): \subpage page_dev_gmxtree.

The Python script can also be invoked to apply changes to \#include directives
to sort and reformat them according to the guidelines.  See the page above for
details.</dd>

<dt>`admin/reformat_all.sh`</dt>
<dd>`bash` script is provided to run uncrustify/`copyright.py`/include sorter
on all relevant files in the source tree (or in a particular directory).
The script can also produce the list of files where these scripts are applied,
for use with other scripts.</dd>

<dt>git attributes</dt>
<dd>git attributes (specified in `.gitattributes` files) are used to annotate
which files are subject to automatic formatting checks (and for automatic
reformatting by the above scripts).  See `man gitattributes` for an overview of
the mechanism.  We use the `filter` attribute to specify the type of automatic
checking/formatting to apply.  Custom attributes are used for specifying some
build system dependencies for easier processing in CMake.</dd>

</dl>

Documentation generation
========================

<dl>
<dt>Doxygen</dt>
<dd>Doxygen (<http://www.doxygen.org>) is used to extract documentation from
source code comments.  Also this developer manual and some other overview
content is laid out by Doxygen from Markdown source files.  Currently, version
1.8.5 is required for a warning-free build.  Thorough explanation of the
Doxygen setup and instructions for documenting the source code can be found on
a separate page: \subpage page_doxygen.</dd>

<dt>graphviz (dot)</dt>
<dd>The Doxygen documentation uses `dot` from graphviz
(<http://www.graphviz.org>) for building some graphs.  The tool is not
mandatory, but the Doxygen build will produce warnings if it is not
available, and the graphs are omitted from the documentation.</dd>

<dt>mscgen</dt>
<dd>The Doxygen documentation uses `mscgen`
(<http://www.mcternan.me.uk/mscgen/>) for building some graphs.  As with `dot`,
the tool is not mandatory, but not having it available will result in warnings
and missing graphs.</dd>

<dt>Doxygen issue checker</dt>
<dd>Doxygen produces warnings about some incorrect uses and wrong
documentation, but there are many common mistakes that it does not detect.
\Gromacs uses an additional, custom Python script to check for such issues.
This is most easily invoked through a `doc-check` target in the build system.
The script also checks that documentation for a header matches its use in the
source code (e.g., that a header documented as internal to a module is not
actually used from outside the module).  Details for the custom checker are
on a separate page (common for several checkers): \subpage page_dev_gmxtree.</dd>

<dt>module dependency graphs</dt>
<dd>\Gromacs uses a custom Python script to generate an annotated dependency
graph for the code, showing \#include dependencies between modules.
The generated graph is embedded into the Doxygen documentation:
\ref page_modulegraph.
This script shares most of its implementation with the custom checkers, and is
documented on the same page: \subpage page_dev_gmxtree.</dd>

<dt>Pandoc</dt>
<dd>Pandoc (<http://johnmacfarlane.net/pandoc/>) is used for building some
parts of the documentation from Markdown source files.</dd>

<dt>latex</dt>
<dd></dd>

<dt>linkchecker</dt>
<dd></dd>

<dt>documentation exported from source files</dt>
<dd>For man pages, HTML documentation of command-line options for executables,
and for shell completions, the `gmx` binary has explicit C++ code to export
the information required.  The build system provides targets that then invoke
the built `gmx` binary to produce these documentation items.  The generated
items are packaged into source tarballs so that this is not necessary when
building from a source distribution (since in general, it will not work in
cross-compilation scenarios).  To build and install these from a git
distribution, explicit action is required.
See \ref page_wrapperbinary for some additional details.</dd>

</dl>
