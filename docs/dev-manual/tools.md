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
<http://www.gromacs.org/index.php?title=Developer_Zone/Git/Git_Tips_%26_Tricks>
Basic tutorial material for `git` can be found on the web.</dd>

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
<dd></dd>

<dt>clang static analyzer</dt>
<dd></dd>

</dl>

Code formatting and style
=========================

The tools and scripts listed below are used to automatically check/apply
formatting that follows \ref page_devstyle.

<dl>

<dt>`uncrustify`</dt>
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

<dt>include directive checking</dt>
<dd></dd>

<dt>include sorting</dt>
<dd></dd>

<dt>reformat_all.sh</dt>
<dd></dd>

<dt>git attributes</dt>
<dd></dd>

</dl>

Documentation generation
========================

<dl>
<dt>Doxygen</dt>
<dd>Thorough explanation of the Doxygen setup and instructions for documenting
the source code can be found on a separate page: \subpage page_doxygen.</dd>

<dt>graphviz (dot)</dt>
<dd></dd>

<dt>mscgen</dt>
<dd></dd>

<dt>Pandoc</dt>
<dd></dd>

<dt>latex</dt>
<dd></dd>

<dt>linkchecker</dt>
<dd></dd>

</dl>
