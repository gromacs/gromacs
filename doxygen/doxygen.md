Using Doxygen {#page_doxygen}
=============

This page documents how Doxygen is set up in the \Gromacs source tree,
as well as guidelines for adding new Doxygen comments.  Examples are included,
as well as tips and tricks for avoiding Doxygen warnings.  The guidelines focus
on C++ code and other new code that follows the new module layout.
Parts of the guidelines are still applicable to documenting older code (e.g.,
within `gmxlib/` or `mdlib/`), in particular the guidelines about formatting
the Doxygen comments and the use of \c \\internal.
\ref page_codelayout documents the overall structure of the documentation.

Documentation flavors
=====================
The \Gromacs source tree is set up to produce three different levels of Doxygen
documentation:

1. Public API documentation (suffix `-user`), which documents functions and
   classes exported from the library and intended for use outside the \Gromacs
   library.
2. Library API documentation (suffix `-lib`), which additionally includes
   functions and classes that are designed to be used from other parts of
   \Gromacs, as well as some guidelines that are mostly of interest to
   developers.
3. Full documentation (suffix `-full`), which includes (nearly) all (documented)
   functions and classes in the source tree.

Each subsequent level of documentation includes all the documentation from the
levels above it.  The suffixes above refer to the suffixes of Doxygen input and
output files, as well as the name of the output directory.  When all the
flavors have been built, the front pages of the documentation contain links to
the other flavors, and explain the differences in more detail.

As a general guideline, the public API documentation should be kept free of
anything that a user linking against an unmodified \Gromacs does not see.
In other words, the public API documentation should mainly document the
contents of installed headers, and provide the necessary overview of using
those.  Also, verbosity requirements for the public API documentation are
higher: ideally, readers of the documentation could immediately start using the
API based on the documentation, without any need to look at the implementation.

Similarly, the library API documentation should not contain things that other
modules in \Gromacs can or should never call.  In particular, anything declared
locally in source files should be only available in the full documentation.
Also, if something is documented, and is not identified to be in the library
API, then it should not be necessary to call that function from outside its
module.

Building the documentation
==========================

If you simply want to see up-to-date documentation, you can go to
http://jenkins.gromacs.org/job/Doxygen_Nightly_master/javadoc/html-lib/index.xhtml
to see the documentation for the current development version.
Jenkins also runs Doxygen for all changes pushed to Gerrit for master, and the
resulting documentation can be viewed from the link posted by Jenkins.  The
Doxygen build is marked as unstable if it introduces any Doxygen warnings.

You may need to build the documentation locally if you want to check the
results after adding/modifying a significant amount of comments.  This is
recommended in particular if you do not have much experience with Doxygen.
It is a good idea to build with all the different settings to see that the
result is what you want, and that you do not produce any warnings.
For local work, it is generally a good idea to set `GMX_COMPACT_DOXYGEN=ON`
CMake option, which removes some large generated graphs from the documentation
and speeds up the process significantly.

All files related to Doxygen reside in the `doxygen/` subdirectory in the source
and build trees.  In a freshly checked out source tree, this directory contains
various `Doxyfile-*.cmakein` files.  When you run CMake, corresponding files
`Doxyfile-user`, `Doxyfile-lib`, and `Doxyfile-full` are generated at the
corresponding location in the build tree.  There is also a
`Doxyfile-common.cmakein`, which is used to produce `Doxyfile-common`.
This file contains settings that are shared between all the input files.
`Doxyfile-compact` provides the extra settings for `GMX_COMPACT_DOXYGEN=ON`.

You can run Doxygen directly with one of the generated files (all output will
be produced under the current working directory), or built one of the
`doc-user`, `doc-lib`, and `doc-full` targets.  The targets run Doxygen in a
quieter mode and only show the warnings if there were any, and put the output
under `doxygen/` in the build tree.  The `doc-all` target builds all three
version.  The generated documentation is under `html-user/`, `html-lib/`,
and/or `html-full/`.  Open `index.xhtml` file from one of these subdirectories
to start browsing.  Log files with all Doxygen warnings are also produced as
`doygen-*.log`, so you can inspect them after the run.

You will need Doxygen 1.8.5 to build the current documentation.  Other versions
may work, but likely also produce warnings.  Additionally,
[`graphviz`][http://www.graphviz.org] and
[`mscgen`][http://www.mcternan.me.uk/mscgen/] are required for some graphs in
the documentation, and `latex` for formulas.  Working versions are likely
available through most package managers..  It is possible to build the
documentation without these tools, but you will see some errors and the related
figures will be missing from the documentation.

General guidelines for Doxygen markup
=====================================

Doxygen provides quite a few different alternative styles for documenting the
source code.  There are subtleties in how Doxygen treats the different types of
comments, and this also depends somewhat on the Doxygen configuration.  It is
possible to change the meaning of a comment by just changing the style of
comment it is enclosed in.  To avoid such issues, and to avoid needing to
manage all the alternatives, a single style throughout the source tree is
preferable.  When it comes to treatment of styles, \Gromacs uses the default
Doxygen configuration with one exception: `JAVADOC_AUTOBRIEF` is set ON to
allow more convenient one-line brief descriptions in C code.

Majority of existing comments in \Gromacs uses Qt-style comments (`/*! */` and
`//!` instead of `/** */` and `///`, \c \\brief` instead of \c \@brief etc.),
so these should be used also for new documentation.  There is a single
exception for brief comments in C code; see below.

Similarly, existing comments use `/*! */` for multiline comments in both C and
C++ code, instead of using multiple `//!` lines for C++.  The rationale is that
since the code will be a mixture of both languages for a long time, it is more
uniform to use similar style in both.  Also, since files will likely transition
from C to C++ gradually, rewriting the comments because of different style
issues should not generally be necessary.  Finally, multi-line `//!` comments
can work differently depending on Doxygen configuration, so it is better to
avoid that ambiguity.

When adding comments, ensure that a short brief description is always produced.
This is used in various listings, and should briefly explain the purpose of the
method without unnecessarily expanding those lists.
The basic guideline is to start all comment blocks with \c \\brief (possibly
after some other Doxygen commands).
If you want to avoid the \c \\brief for one-liners, you can use `//!`, but the
description must fit on a single line; otherwise, it is not interpreted as a
brief comment.  Note in particular that a simple `/*! Comment. */` does not
produce a brief description.
Also note that \c \\brief marks the whole following paragraph as a brief
description, so you should insert an empty line after the intended brief
description.

In C code, `//` comments must be avoided because some compilers do not like
them.  If you want to avoid the \c \\brief for one-liners in C code, use
`/** */` instead of `//!`.  If you do this, the brief description should not
contain unescaped periods except at the end.  Because of this, you should
prefer `//!` in C++ code.

Put the documentation comments in the header file that contains the
declaration, if such a header exists.
Implementation-specific comments that do not influence how a method
is used can go into the source file, just before the method definition, with an
\c \\internal tag in the beginning of the comment block.  Doxygen-style comments
within functions are not generally usable.

At times, you may need to exclude some part of a header or a source file such
that Doxygen does not see it at all.  In general, you should try to avoid this,
but it may be necessary to remove some functions that you do not want to appear
in the public API documentation, and which would generate warnings if left
undocumented, or to avoid Doxygen warnings from code it does not understand.
Prefer \c \\cond and \c \\endcond to do this.  If \c \\cond does not work for
you, you can also use \c \#ifndef `DOXYGEN`.  If you exclude a class method in
a header, you also need to exclude it in the source code to avoid warnings.


GROMACS specifics
=================

TODO

Documenting specific code constructs
====================================

Classes/structs
---------------

TODO

Methods/functions
-----------------

TODO

Files
-----

TODO

Modules
-------

TODO

Examples
========

Here is an example of documenting a C++ class and its containing header file.
Comments in the code and the actual documentation explain the used Doxygen
constructs.

\includelineno doxygen-example.cpp

Here is another example of documenting a C header file (so avoiding all
C++-style comments), and including free functions. It also demonstrates the use
of \c \\addtogroup to add multiple functions into a module group without repeated
\c \\ingroup tags.

\includelineno doxygen-example.c

More examples you can find by looking at existing code in the source tree.  In
particular new C++ code such as that in the `src/gromacs/analysisdata/` and
`src/gromacs/options/` subdirectories contains a large amount of code
documented along these guidelines. Some comments in `src/gromacs/selection/`
(in particular, any C-like code) predate the introduction of these guidelines,
so those are not the best examples.
