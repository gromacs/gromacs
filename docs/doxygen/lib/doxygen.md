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

To get started quickly, you only need to read the first two sections to
understand the overall structure of the documentation, and take a look at the
examples at the end.  The remaining sections provide the details for
understanding why the examples are the way they are, and for more complex
situations.  They are meant more as a reference to look up solutions for
particular problems, rather than single-time reading.  To understand or find
individual Doxygen commands, you should first look at Doxygen documentation
(http://www.stack.nl/~dimitri/doxygen/manual/index.html).


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
http://jenkins.gromacs.org/job/Documentation_Nightly_master/javadoc/html-lib/index.xhtml
to see the documentation for the current development version.
Jenkins also runs Doxygen for all changes pushed to Gerrit for
release-5-0 and master branches, and the
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

All files related to Doxygen reside in the `docs/doxygen/` subdirectory in the source
and build trees.  In a freshly checked out source tree, this directory contains
various `Doxyfile-*.cmakein` files.  When you run CMake, corresponding files
`Doxyfile-user`, `Doxyfile-lib`, and `Doxyfile-full` are generated at the
corresponding location in the build tree.  There is also a
`Doxyfile-common.cmakein`, which is used to produce `Doxyfile-common`.
This file contains settings that are shared between all the input files.
`Doxyfile-compact` provides the extra settings for `GMX_COMPACT_DOXYGEN=ON`.

You can run Doxygen directly with one of the generated files (all output will
be produced under the current working directory), or build one of the
`doxygen-user`, `doxygen-lib`, and `doxygen-full` targets.  The targets run
Doxygen in a quieter mode and only show the warnings if there were any, and put
the output under `docs/html/doxygen/` in the build tree, so that the Doxygen
build cooperates with the broader `webpage` target.
The `doxygen-all` target builds all three targets with less typing.

The generated documentation is put under `html-user/`, `html-lib/`, and/or
`html-full/`.  Open `index.xhtml` file from one of
these subdirectories to start browsing (for \Gromacs developers, the
`html-lib/` is a reasonable starting point).  Log files with all Doxygen
warnings are also produced as `docs/doxygen/doxygen-*.log`, so you can inspect them after
the run.

You will need Doxygen 1.8.5 to build the current documentation.  Other versions
may work, but likely also produce warnings.  Additionally,
[graphviz](http://www.graphviz.org) and
[mscgen](http://www.mcternan.me.uk/mscgen/) are required for some graphs in
the documentation, and `latex` for formulas.  Working versions are likely
available through most package managers.  It is possible to build the
documentation without these tools, but you will see some errors and the related
figures will be missing from the documentation.


General guidelines for Doxygen markup {#section_doxygen_guidelines}
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

Majority of existing comments in \Gromacs uses Qt-style comments (`/*!` and
`//!` instead of `/*``*` and `///`, \c \\brief instead of \c \@brief etc.),
so these should be used also for new documentation.  There is a single
exception for brief comments in C code; see below.

Similarly, existing comments use `/*!` for multiline comments in both C and
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
brief comment.  Note in particular that a simple `/*!` without a \c \\brief
does not produce a brief description.
Also note that \c \\brief marks the whole following paragraph as a brief
description, so you should insert an empty line after the intended brief
description.

In C code, `//` comments must be avoided because some compilers do not like
them.  If you want to avoid the \c \\brief for one-liners in C code, use
`/*``*` instead of `//!`.  If you do this, the brief description should not
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


\Gromacs specifics
=================

The general guidelines on the style of Doxygen comments were given above.
This section introduces \Gromacs specific constructs currently used in Doxygen
documentation, as well as how \Gromacs uses Doxygen groups to organize the
documentation.

Some consistency checks are done automatically using custom scripts.
See \ref page_dev_gmxtree for details.

Controlling documentation visibility
------------------------------------

To control in which level of documentation a certain function appears, three
different mechanisms are used:
* Global Doxygen configuration.  This is mainly used to include
  declarations local to source files only in the full documentation.
  You can find the details from the `Doxyfile-*.cmakein` files, and some of
  them are also mentioned below on individual code constructs.
* The standard Doxygen command \c \\internal marks the documentation to be only
  extracted into the full documentation (`INTERNAL_DOCS` is `ON` only for the
  full documentation).  This should be used as a first command in a comment
  block to exclude all the documentation.  It is possible to use \c \\internal
  and \c \\endinternal to exclude individual paragraphs, but \c \\if `internal`
  is preferred (see below).
  In addition, \Gromacs-specific custom Doxygen command \c \\libinternal is
  provided, which should be used the same way to exclude the documentation from
  the public API documentation.  This command expands to either \c \\internal
  or to a no-op, depending on the documentation level.
* Doxygen commands \c \\if and \c \\cond can be used with section names
  `libapi` and `internal` to only include the documentation in library API and
  the full documentation, respectively.  `libapi` is also defined in the full
  documentation.  These are declared using `ENABLED_SECTIONS` in the Doxygen
  configuration files.

Examples of locations where it is necessary to use these explicit commands are
given below in the sections on individual code constructs.

Modules as Doxygen groups
-------------------------

As described in \ref page_codelayout, each subdirectory under `src/gromacs/`
represents a _module_, i.e., a somewhat coherent collection of routines.
Doxygen cannot automatically generate a list of routines in a module; it only
extracts various alphabetical indexes that contain more or less all documented
functions and classes.  To help reading the documentation, the routines for a
module should be visible in one place.

\Gromacs uses Doxygen groups to achieve this: for each documented module, there
is a \c \\defgroup definition for the module, and all the relevant classes and
functions need to be manually added to this group using \c \\ingroup and
\c \\addtogroup.
The group page also provides a natural place for overview documentation about
the module, and can be navigated to directly from the "Modules" tab in the
generated documentation.

Some notes about using \c \\addtogroup are in order:
* \c \\addtogroup only adds the elements that it directly contains into the
  group.  If it contains a namespace declaration, only the namespace is added
  to the group, but none of the namespace contents are.  For this reason,
  \c \\addtogroup should go within the innermost scope, around the members that
  should actually be added.
* If the module should not appear in the public API documentation, its
  definition (\c \\defgroup) should be prefixed with a \c \\libinternal.
  In this case, also all \c \\addtogroup commands for this module should be
  similarly prefixed.  Otherwise, they create the group in the public API
  documentation, but without any of the content from the \c \\defgroup
  definition.  This may also cause the contents of the \c \\addtogroup section
  to appear in the public API documentation, even if it otherwise would not.

Public API and library API groups
---------------------------------

In addition to the module groups, two fixed groups are provided:
`group_publicapi` and `group_libraryapi`.  Classes and files can be added to
these groups using \Gromacs specific custom \c \\inpublicapi and
\c \\inlibraryapi commands.  The generated group documentation pages are not
very useful, but annotated classes and files show the API definition under the
name, making this information more easily accessible.  These commands in
file-level comments are also used for some automatic intermodule dependency
validation (see below).

Note that functions, enumerations, and other entities that do not have a
separate page in the generated documentation can only belong to one group;
in such a case, the module group is preferred over the API group.


Documenting specific code constructs
====================================

This section describes the techical details and some tips and tricks for
documenting specific code constructs such that useful documentation is
produced.  If you are wondering where to document a certain piece of
information, see the documentation structure section on \ref page_codelayout.
The focus of the documentation should be on the overview content: Doxygen pages
and the module documentation.  An experienced developer can relatively easily
read and understand individual functions, but the documentation should help
in getting the big picture.

Doxygen pages
-------------

The pages that are accessible through navigation from the front page are
written using Markdown and are located under `docs/doxygen/`.  Each page should be
placed in the page hierarchy by making it a subpage of another page, i.e., it
should be referenced once using \c \\subpage.  `mainpage.md` is the root of the
hierarchy.

There are two subdirectories, `user/` and `lib/`, determining the highest
documentation level where the page appears.  If you add pages to `lib/`, ensure
that there are no references to the page from public API documentation.
\c \\if `libapi` can be used to add references in content that is otherwise
public.
Generally, the pages should be on a high enough level and provide overview
content that is useful enough such that it is not necessary to exclude them
from the library API documentation.

Modules
-------

For each module, decide on a header file that is the most important one for
that module (if there is no self-evident header, it may be better to designate,
e.g., module-doc.h for this purpose, but this is currently not done for any
module).  This header should contain the \c \\defgroup definition for the
module.  The name of the group should be `module_`<em>name</em>, where _name_
is the name of the subdirectory that hosts the module.

The module should be added to an appropriate group (see `docs/doxygen/misc.cpp` for
definitions) using \c \\ingroup to organize the "Modules" tab in the generated
documentation.

One or more contact persons who know about the contents of the module should be
listed using \c \\author commands.  This provides a point of contact if one
has questions.

Classes/structs
---------------

Classes and structs in header files appear always in Doxygen documentation,
even if their enclosing file is not documented.  So start the documentation
blocks of classes that are not part of the public API with \c \\internal or
\c \\libinternal.  Classes declared locally in source files or in unnamed
namespaces only appear in the full documentation.

If a whole class is not documented, this does not currently generate any
warning.  The class is simply exluded from the documentation.  But if a member
of a documented class is not documented, a warning is generated.  Guidelines for
documenting free functions apply to methods of a class as well.

For base classes, the API classification (\c \\inpublicapi or
\c \\inlibraryapi) should be based on where the class is meant to be
subclassed.  The visibility (\c \\internal or \c \\libinternal), in contrast,
should reflect the API classification of derived classes such that the base
class documentation is always generated together with the derived classes.

For classes that are meant to be subclassed and have protected members, the
protected members should only appear at the documentation level where the class
is meant to be subclassed.  For example, if a class is meant to be subclassed
only within a module, the protected members should only appear in the
full documentation.  This can be accomplished using \c \\cond (note that you
will need to add the \c \\cond command also to the source files to hide the
same methods from Doxygen, otherwise you will get confusing warnings).

Methods/functions/enums/macros
------------------------------

These items do not appear in the documentation unless their enclosing scope is
documented.  For class members, the scope is the class; otherwise, it is the
namespace if one exists, or the file.  An \c \\addtogroup can also define a
scope if the group has higher visibility than the scope outside it.
So if a function is not within a namespace (mostly applicable to C code) and
has the same visibility as its enclosing file, it is not necessary to add a
\c \\internal or \c \\libinternal.

Static functions are currently extracted for all documentation flavors to allow
headers to declare `static inline` functions (used in, for example, math code).
Functions in anonymous namespaces are only extracted into the full
documentation.  Together with the above rules, this means that you should avoid
putting a `static` function within a documented namespace, even within source
files, or it may inadvertently appear in the public API documentation.

If you want to exclude an item from the documentation, you need to put in
inside a \c \\cond block such that Doxygen does not see it.
Otherwise, a warning for an undocumented function is generated.  You need to
enclose both the declaration and the definition with \c \\cond.

Files
-----

Each documented file should start with a documentation block (right after the
copyright notice) that documents the file.  See the examples section for exact
formatting.  Things to note:
* Please do not specify the file name explicitly after \c \\file.  By default,
  a file comment applies to the file it is contained in, and an explicit file
  name only adds one more thing that can get out of date.
* \c \\brief cannot appear on the same line as the \c \\file, but should be on
  the next line.
* \c \\internal or \c \\libinternal should indicate where the header is visible.
  As a general guideline, all installed headers should appear in the public API
  documentation, i.e., not contain these commands.  If nothing else, then to
  document that it does not contain any public API functions.  Headers that
  declare anything in the library API should be marked with \c \\libinternal,
  and the rest with \c \\internal.
* All source files, as well as most test files, should be documented with
  \c \\internal, since they do not provide anything to public or library API,
  and this avoids unintentionally extracting things from the file into those
  documentations.  Shared test files used in tests from other modules should be
  marked with \c \\libinternal.
* \c \\inpublicapi or \c \\inlibraryapi should be used to indicate where the
  header is meant to be directly included.
* As with files, one or more contact persons should be listed with \c \\author.
  If you make significant modifications or additions to a file, consider adding
  an \c \\author line for yourself.

Directories
-----------

Directory documentation does not typically contain useful information beyond a
possible brief description, since they correspond very closely to modules, and
the modules themselves are documented.  A brief description is still useful to
provide a high-level overview of the source tree on the generated "Files" page.
A reference to the module is typically sufficient as a brief description for a
directory.  All directories are currently documented in
`docs/doxygen/directories.cpp`.


Examples
========

Basic C++
---------

Here is an example of documenting a C++ class and its containing header file.
Comments in the code and the actual documentation explain the used Doxygen
constructs.

\includelineno doxygen-example.cpp

Basic C
-------

Here is another example of documenting a C header file (so avoiding all
C++-style comments), and including free functions.  It also demonstrates the use
of \c \\addtogroup to add multiple functions into a module group without repeated
\c \\ingroup tags.

\includelineno doxygen-example.c

Scoping and visibility rules
----------------------------

The rules where Doxygen expects something to be documented, and when are
commands like \c \\internal needed, can be complex.  The examples below
describe some of the pitfalls.

\includelineno doxygen-example-scoping.cpp

Module documentation
--------------------

Documenting a new module should place a comment like this in a central header
for the module, such that the "Modules" tab in the generated documentation can
be used to navigate to the module.

\includelineno doxygen-example-module.cpp

Common mistakes
---------------

The most common mistake, in particular in C code, is to forget to document the
file.  This causes Doxygen to ignore most comments in the file, so it
does not validate the contents of the comments either, nor is it possible to
actually check how the generated documentation looks like.

The following example shows some other common mistakes (and some less common)
that do not produce correct documentation, as well as Doxygen "features"/bugs
that can be confusing.
The issues are explained in normal comments above each code fragment.

\includelineno doxygen-example-issues.cpp

Existing code
-------------

More examples you can find by looking at existing code in the source tree.  In
particular new C++ code such as that in the `src/gromacs/analysisdata/` and
`src/gromacs/options/` subdirectories contains a large amount of code
documented mostly along these guidelines.  Some comments in
`src/gromacs/selection/` (in particular, any C-like code) predate the
introduction of these guidelines, so those are not the best examples.
