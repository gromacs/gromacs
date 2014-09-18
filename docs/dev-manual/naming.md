Naming conventions {#page_devstyle_naming}
==================

The conventions here should be applied to all new code, and with common sense
when modifying existing code.  For example, renaming a widely used, existing
function to follow these conventions may not be justified unless the whole code
is getting a rework.

Currently, this only documents the present state of the code: no particular
attempt has been made to consolidate the naming.

Files
=====

 * Avoid having multiple files with the same name in different places within
   the same library.  In addition to making things harder to find, C++ source
   files with the same name can cause obscure problems with some compilers.
   Currently, unit tests are an exception to the rule (there is only one
   particular compiler that had problems with this, and a workaround is
   possible if/when that starts to affect more than a few of the test files).

TODO: Copy/move relevant content from \ref page_codelayout here, and add
details.

Common guidelines for C and C++ code
====================================

 * Preprocessor macros should be all upper-case.  Do not use leading
   underscores, as all such names are reserved according to the C/C++ standard.
 * Name include guards like `GMX_DIRNAME_HEADERNAME_H`.
 * Boolean variables are always named with a `b` prefix, followed by a
   CamelCase name.
 * Enum values are named with an `e` prefix.
 * Avoid abbreviations that are not obvious to a general reader.

C code
======

 * All function and variable names are lowercase, with underscores as word
   separators where needed for clarity.
 * All functions that are part of the public API should start with `gmx_`.
   Preferably, other functions should as well.
   Some parts of the code use a `_gmx_` prefix for internal functions, but
   strictly speaking, these are reserved names, so, e.g., a trailing underscore
   would be better.

C++ code
========

 * Use CamelCase for all names.  Start types (such as classes, structs, and
   typedefs) with a capital letter, other names (functions, variables) with a
   lowercase letter.
   You may use an all-lowercase name with underscores if your class closely
   resembles an external construct (e.g., a standard library construct) named
   that way.
 * C++ interfaces are named with a `Interface` suffix, and abstract base
   classes with an `Abstract` prefix.
 * Member variables are named with a trailing underscore.
 * Accessors for a variable `foo_` are named `foo()` and `setFoo()`.
 * Global variables are named with a `g_` prefix.
 * Static class variables are named with a `s_` prefix.
 * Global constants are often named with a `c_` prefix.
 * If the main responsibility of a file is to implement a particular class,
   then the name of the file should match that class, except for possible
   abbreviations to avoid repetition in file names (e.g., if all classes within
   a module start with the module name, omitting or abbreviating the module
   name is OK).  Currently, all source file names are lowercase, but this
   casing difference should be the only difference.

The rationale for the trailing underscore and the global/static prefixes is
that it is immediately clear whether a variable referenced in a method is local
to the function or has wider scope, improving the readability of the code.

Unit tests
==========

 * Test fixtures (the first parameter to `TEST`/`TEST_F`) are named with a
   `Test` suffix.
 * Classes meant as base classes for test fixtures (or as names to be typedefed
   to be fixtures) are named with a `TestBase` of `Fixture` suffix.
 * The CTest test is named with CamelCase, ending with `Tests` (e.g.,
   `OptionsUnitTests`).
 * The test binary is named with the name of the module and a `-test` suffix.
