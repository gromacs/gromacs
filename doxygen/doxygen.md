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
The Gromacs source tree is set up to produce three different levels of Doxygen documentation:

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
output files, as well as the name of the output directory.  When the full
documentation is built, the front pages of the documentation contain links to
the other flavors, and explain the differences in more detail.

As a general guideline, the public API documentation should be kept free of
anything that a user linking against an unmodified \Gromacs does not see.
The public API documentation should mainly document the contents of installed
headers, and provide the necessary overview of using those.

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
`doygen-*.log`.

You will need Doxygen 1.8.5 to build the current documentation.  Other versions
may work, but likely also produce warnings.  Additionally, `graphviz` and
`mscgen` are required for some graphs in the documentation, and `latex` for
formulas.  It is possible to build the documentation without these tools, but
you will see some errors and the related figures will be missing from the
documentation.

General guidelines for Doxygen markup
=====================================

Doxygen provides quite a few different alternative styles for documenting the
source code, and the style also depends somewhat on the Doxygen configuration.
\Gromacs mainly uses Qt-style comments (`/*! */` and `//!` instead of `/** */`
and `///`, `\\brief` instead of `\@brief` etc.), so these should be preferred for
new documentation.

When adding comments, ensure that a short brief description is produced.
This is used in various listings, and should briefly explain the purpose of the
method without unnecessarily expanding those lists.
For one-liners, you can use `//!`, but the description must fit on a single
line; otherwise, it is not interpreted as a brief comment.
Note in particular that a simple `/*! Comment. */` does not produce a brief
description.
Also note that `\\brief` marks the whole following paragraph as a brief
description, so you should insert an empty line after the intended brief
description.

In C code, `//` comments must be avoided because some compilers do not like
them.  You can specify an one-line brief comment using `/** */` in C code.  The
brief description should not contain unescaped periods except at the end,
otherwise it does not work as expected.

The preferred location for the documentation comments is in header files, if
one exists.  Implementation-specific comments that do not influence how a method
is used can go into the source file, just before the method definition, with an
\c \\internal tag in the beginning of the comment block.  Doxygen-style comments
within methods are not generally usable.

If you need to exclude some part of a header or a source file from Doxygen,
prefer \c \\cond and \c \\endcond.  If you exclude a class method in a header,
you also need to exclude it in the source code to avoid warnings.
If \c \\cond does not work for you, you can also use \c \#ifndef `DOXYGEN`.

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
