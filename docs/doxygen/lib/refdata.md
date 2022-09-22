Using reference data in C++ tests {#page_refdata}
=================================

The \ref module_testutils module provides (among other things) utilities to
write Google Test tests that compare their results against stored reference
data.  This can either be used for

 * regression-style tests, just ensuring that the output does not change, or
 * combined with manual checking of the reference data, as a different kind of
   assertion, where the expected results would be tedious to express directly
   as C++ code (e.g., when checking complicated data structures for
   correctness).

The current reference data functionality is quite basic, but it can be extended
if/when more control over, e.g., comparison tolerances is needed.

Reference data organization
===========================

Conceptually, the reference data consists of a tree-like structure of nodes.
Each leaf node checks a single primitive value (an integer, a floating-point
value, a string etc.), and each inner node acts as a _compound_ value that
helps organizing the data.  Within each compound node (including the root of
the tree), child nodes are identified by an `id` string.  Each node within a
single compound must have a unique `id`, and it is possible to compare
multiple values produced by the test against this single node (naturally, the
test only passes if the test produces the same value in all such cases).

Each node also has a type (a string).  For leaf nodes, the type is from a
predetermined set of strings, and identifies the type of the value stored in
the node.  For compound nodes, the type is just a string provided by the test.
In all cases, the type in the reference data must match the type provided by
the test.  This provides additional safety when changing the test to detect
mismatches between the test and the reference data.  The intention is that
compound nodes whose contents have the same structure would have the same type;
this will simplify using XSLT for viewing the reference data (see below).

Some compound types are predefined, e.g., for simple sequences, but more
complicated compounds can be defined ad-hoc in tests that need them.  See below
for how to use them in the code.

As a special case, the `id` can be empty (`NULL`).  This is intended for
cases where one is checking for a sequence of items, and the only thing
distinguishing the items is their position in this sequence.  Using an empty
`id` removes the need to generate unique identifiers for the items, and makes
textual diffs of the reference data files easier to read.
Only a single sequence of nodes with an empty `id` is supported within one
parent node: if you first check some nodes with an empty `id`, followed by a
non-empty `id`, the next check for an empty `id` will again match the first
node in the sequence.
For clarity, all the nodes that have an empty `id` should be of the same
type, but this is not enforced.

Using reference data in code
============================

To use reference data in a test, the test should first create exactly one
instance of gmx::test::TestReferenceData.  It can do so as a local variable in
the test, as a member variable in its test fixture, or by subclassing a test
fixture that already contains such a variable (e.g., gmx::test::StringTestBase
or gmx::test::CommandLineTestBase).
Only use the default constructor!  The other constructor is intended for
self-testing utility code used in other tests (including self-testing the
reference data implementation itself), and behaves differently from what is
described here.

To access the root node of the data,
gmx::test::TestReferenceData::rootChecker() needs to be called.
This returns a gmx::test::TestReferenceChecker that provides various
`check*()` methods that can be used to check values against top-level nodes.
gmx::test::TestReferenceChecker::checkCompound() can be called to create custom
compound types: it returns another gmx::test::TestReferenceChecker that can be
used to check values against child nodes of the created compound.

Whenever a gmx::test::TestReferenceChecker method detects a mismatch against
reference data, it will generate a non-fatal Google Test failure in the current
test.  The test can naturally also use its own test assertions for additional
checks, but any mismatch will automatically also fail the test.

It is also possible to read values of the reference data items using
gmx::test::TestReferenceChecker, so that they can be used programmatically.
For this to work, those items should first be written in the same test.
This supports tests that want to both check data against a reference, and use
that reference as a persistence layer for storing information.  This is useful
at least for serialization tests.
This is currently not supported for all use cases, but with some caveats, it is
possible to use this for testing.

When using floating-point values in reference data, the tolerance for the
comparison can be influenced with
gmx::test::TestReferenceChecker::setDefaultTolerance().
Per-comparison tolerances would be possible to implement if necessary, but
currently you can either change the default tolerance whenever you need to, or
create copies of the gmx::test::TestReferenceChecker object and set different
tolerances in the different instances.  Note that there is an implicit
assumption that a mixed- and a double-precision build will produce the same
results (within the given tolerance).  This means that some things cannot be
tested with the reference data (e.g., multiple steps of MD integration), and
that reference data for such tests needs to be always generated in double
precision (unless the results are nice, exact binary floating-point numbers).

Just creating a gmx::test::TestReferenceData instance does not enforce using
reference data in the test; the data is loaded/used only when
gmx::test::TestReferenceData::rootChecker() is first called.  If the test never
calls this method, the gmx::test::TestReferenceData object does nothing.  This
allows using the same test fixture (e.g., CommandLineTestBase) also in tests
that do not need the reference data, but benefit from other features of the
fixture.

Running tests that use reference data
=====================================

To run a test that uses the reference data, you just execute the test binary as
you would otherwise.  However, when you first add a test, the reference data
does not exist, and the test will fail with an assertion message saying that
the reference data could not be found.  To generate the reference data, you
need to run the test binary with a `-ref-data create` command-line option
(it is also possible to use any of the `update` options below to generate the
reference data).

If you change a test (or the tested code) such that the reference data needs to
be changed, you need to run the test binary with `-ref-data update-all` or
`-ref-data update-changed`.  The first will recreate the reference data from
scratch.  The latter will retain old reference values if they are still valid.
In other words, floating-point reference values that are within the test
tolerance will be kept at their old values.  Only values that are outside the
tolerance (or otherwise do not match or do not exist) are updated.
This is useful (at least) for tests that contain floating-point data, where it
is not expected that those floating-point values would actually need to change.
This allows you to update other parts of the reference data without doing a
double-precision build, and also makes it easier to avoid spurious changes in
the last bits of other reference data values when just a single output value is
expected to change.

To create or update reference data, the test needs to pass when run with the
corresponding flag.  All comparisons against reference data will pass in these
modes, but you need to ensure that other assertions in the test also pass, and
that the test does not throw exceptions.
Note that if your test does multiple comparisons against the same `id` node,
reference data comparison can still fail during create/update if the test does
not produce the same results for each comparison.

With all the operations that create or update the reference data, you can use
the `--gtest_filter=<...>` command-line option provided by Google Test to
select the tests whose reference data you want to influence.

Persistence
===========

The reference data is stored in XML files under
`src/gromacs/`<em>module</em>`/tests/refdata/` in the source tree.
This part of the framework depends on `tinyxml2`, which is bundled in `src/external`.
One file is produced per test that uses reference data.  If you rename tests or
otherwise change the reference data, you currently need to manually manage the
files with `git`.

For inspecting the reference data in a browser, there are XSLT stylesheets that
transform the XML files into HTML.  Such custom transformations need to be
written for each type of test if the output is not easy to check otherwise.
Because of security features in browsers, the transformations may not work for
all browsers.  For the same reason, the XSLT files must be in the same folder
as the XML files.  For cases where the XSLT files are shared between multiple
modules, `src/testutils/copy_xsl.sh` takes care to synchronize the files after
a main copy is edited.
