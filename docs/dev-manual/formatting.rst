Guidelines for code formatting
==============================

The following list provides the general formatting/indentation rules for
|Gromacs| code (C/C++):

* Basic indentation is four spaces.
* Keep lines at a reasonable length.  There is no hard limit, but use 80
  characters as a guideline.  If you end up indenting very deeply,
  consider splitting the code into functions.
* Do not use tabs, only spaces.  Most editors can be configured to generate
  spaces even when pressing tab.  Tabs (in particular when mixed with spaces)
  easily break indentation in contexts where settings are not exactly equal
  (e.g., in ``git diff`` output).
* No trailing whitespace.
* Use braces always for delimiting blocks, even when there is only a single
  statement in an ``if`` block or similar.
* Put braces on their own lines.  The only exception is short one-line inline
  functions in C++ classes, which can be put on a single line.
* Use spaces liberally.
* ``extern "C"`` and ``namespace`` blocks are not indented, but all others
  (including ``class`` and ``switch`` bodies) are.

Additionally:

* All source files and other non-trivial scripts should contain a copyright
  header with a predetermined format and license information (check existing
  files).  Copyright holder should be "the GROMACS development team" for the
  years where the code has been in the |Gromacs| source repository, but earlier
  years can hold other copyrights.
* Whenever you update a file, you should check that the current year is listed
  as a copyright year.

Most of the above guidelines are enforced using uncrustify, an automatic source
code formatting tool.  The copyright guidelines are enforced by a separate
Python script.  See :doc:`uncrustify` for details.  Note that due to the
nature of uncrustify (it only does all-or-nothing formatting), it enforces
several additional formatting rules in addition to those above.

Enforcing a consistent formatting has a few advantages:

* No one needs to manually review code for most of these formatting issues,
  and people can focus on content.
* A separate automatic script (see below) can be applied to re-establish the
  formatting after refactoring like renaming symbols or changing some
  parameters, without needing to manually do it all.
