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

Code formatting and style {#section_dev_formattingtools}
=========================

The tools and scripts listed below are used to automatically check/apply
formatting that follows \Gromacs style guidelines described on a separate page:
\ref page_devstyle.

<dl>

<dt>uncrustify</dt>
<dd>Uncrustify (<http://uncrustify.sourceforge.net>) is used for automatic
indentation and other formatting of the source code to follow
\ref page_devstyle_formatting.  All code must remain invariant under uncrustify
with the config at `admin/uncrustify.cfg`.  A patched version of uncrustify is
used.  See \ref page_dev_uncrustify for details.</dd>

<dt>`admin/copyright.py`</dt>
<dd>This Python script adds and formats copyright headers in source files.
`uncrustify.sh` (see below) uses the script to check/update copyright years on
changed files automatically.</dd>

<dt>`admin/uncrustify.sh`</dt>
<dd>This `bash` script runs uncrustify and `copyright.py` for all
files that have local changes and checks that they conform to the prescribed
style.  Optionally, the script can also apply changes to make the files
conform.
This script is automatically run by Jenkins to ensure that all commits adhere
to \ref page_devstyle_formatting.  If the uncrustify job does not succeed, it
means that this script has something to complain.
See \ref page_dev_uncrustify for details.</dd>

<dt>`admin/git-pre-commit`</dt>
<dd>This sample git pre-commit hook can be used if one wants to apply
`uncrustify.sh` automatically before every commit to check for formatting
issues.  See \ref page_dev_uncrustify for details.</dd>

<dt>`docs/doxygen/includesorter.py`</dt>
<dd>This Python script sorts and reformats \#include directives according to
the guidelines at \ref page_devstyle_includes.  Details are documented on a
separate page (with the whole suite of Python scripts used for source code
checks): \ref section_dev_includesorter.</dd>

<dt>include directive checker</dt>
<dd>In its present form, the above include sorter script cannot be conveniently
applied in `uncrustify.sh`.  To check for issues, it is instead integrated into
a `check-source` build target.  When this target is built, it also checks for
include formatting issues.  Internally, it uses the sorter script.  This check
is run in Jenkins as part of the Documentation job.
Details for the checking mechanism are on a separate page (common for several
checkers): \subpage page_dev_gmxtree.</dd>

<dt>`admin/reformat_all.sh`</dt>
<dd>This `bash` script runs uncrustify/`copyright.py`/include sorter
on all relevant files in the source tree (or in a particular directory).
The script can also produce the list of files where these scripts are applied,
for use with other scripts.  See \ref page_dev_uncrustify for details.</dd>

<dt>git attributes</dt>
<dd>git attributes (specified in `.gitattributes` files) are used to annotate
which files are subject to automatic formatting checks (and for automatic
reformatting by the above scripts).  See `man gitattributes` for an overview of
the mechanism.  We use the `filter` attribute to specify the type of automatic
checking/formatting to apply.  Custom attributes are used for specifying some
build system dependencies for easier processing in CMake.</dd>

<dt>include-what-you-use</dt>
<dd></dd>

</dl>

Documentation generation
========================

Building the GROMACS documentation
----------------------------------

For now, there are multiple components, formats and tools for the
GROMACS documentation, which is aimed primarily at version-specific
deployment of the complete documentation on the website.

This is quite complex, because the dependencies for building the
documentation must not get in the way of building the code
(particularly when cross-compiling), and yet the code must build and
run in order for some documentation to be generated. Also, man page
documentation (and command-line completions) must be built from the
wrapper binary, in order to be bundled into the tarball.

The outputs of interest to most developers are generally produced in the
<tt>docs/html</tt> subdirectory of the build tree.

You need to enable at least some of the following CMake options:
<dl>
<dt><tt>GMX_BUILD_MANUAL</tt></dt>
<dd>Option needed for trying to build the PDF reference manual
(requires LaTeX and ImageMagick)</dd>
<dt><tt>GMX_BUILD_HELP</tt></dt>
<dd>Option that builds the \ref page_wrapperbinary and use it to generate
input for Sphinx</dd>
<dt><tt>GMX_BUILD_WEBPAGE</tt></dt>
<dd>Option needed for compiling all the documentation into the webpage</dd>
</dl>
Documentation cannot be built if the CMake option
<tt>GMX_BUILD_MDRUN_ONLY</tt> is enabled.

The following make targets are the most useful:
<dl>
<dt><tt>manual</tt></dt>
<dd>Builds the PDF reference manual</dd>
<dt><tt>man</tt></dt>
<dd>Makes man pages from the wrapper binary with Sphinx</dd>
<dt><tt>doxygen-all</tt></dt>
<dd>Makes the code documentation with Doxygen</dd>
<dt><tt>install-guide</tt></dt>
<dd>Makes the INSTALL file for the tarball with Sphinx</dd>
<dt><tt>webpage-sphinx</tt></dt>
<dd>Makes all the components of the GROMACS webpage that require Sphinx,
including install guide and user guide.</dd>
<dt><tt>webpage</tt></dt>
<dd>Makes the complete GROMACS webpage, requires everything. When complete,
you can browse <tt>docs/html/index.html</tt> to find everything.

If built from a release tarball, the <tt>SOURCE_MD5SUM</tt>,
<tt>SOURCE_TARBALL</tt>, <tt>REGRESSIONTESTS_MD5SUM</tt>, and
<tt>REGRESSIONTESTS_TARBALL</tt> CMake variables can be set to pass in
the md5sum values and names of those tarballs, for embedding into the
final deployment to the GROMACS website.</dd>
</dl>

The following tools are used in building parts of the documentation.

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
This is most easily invoked through a `check-source` target in the build system.
The script also checks that documentation for a header matches its use in the
source code (e.g., that a header documented as internal to a module is not
actually used from outside the module).  These checks are run in Jenkins as
part of the Documentation job.  Details for the custom checker are on a
separate page (common for several checkers): \subpage page_dev_gmxtree.</dd>

<dt>module dependency graphs</dt>
<dd>\Gromacs uses a custom Python script to generate an annotated dependency
graph for the code, showing \#include dependencies between modules.
The generated graph is embedded into the Doxygen documentation:
\ref page_modulegraph.
This script shares most of its implementation with the custom checkers, and is
documented on the same page: \subpage page_dev_gmxtree.</dd>

<dt>Sphinx</dt>
<dd>Sphinx (<http://sphinx-doc.org/>; at least version 1.23) is used
for building some parts of the documentation from reStructuredText
source files.</dd>

<dt>latex</dt>
<dd>Also requires ImageMagick for converting graphics file formats</dd>

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
