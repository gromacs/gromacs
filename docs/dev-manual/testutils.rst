Unit testing
============

The main goal of unit tests in |Gromacs| is to help developers while developing
the code.  They focus on testing functionality of a certain module or a group
of closely related modules.  They are designed for quick execution, such that
they are easy to run after every change to check that nothing has been broken.

Finding, building and running
-----------------------------

As described in :ref:`dev-source-layout`, ``src/gromacs/`` is divided into modules,
each corresponding to a subdirectory.  If available, unit tests for that module
can be found in a ``tests/`` subdirectory under the top-level module directory.
Typically, tests for code in :file:`{file}.h` in the module is in a corresponding
:file:`tests/{file}.cpp`.  Not all files have corresponding tests, as it may not
make sense to test that individual file in isolation.  Focus of the tests is on
functionality exposed outside the module.  Some of the tests, in particular for
higher-level modules, are more like integration tests, and test the
functionality of multiple modules.
Shared code used to implement the tests is in ``src/external/gmock-1.7.0/`` and
``src/testutils/`` (see below).

The tests are built if ``BUILD_TESTING=ON`` (the default) and
``GMX_BUILD_UNITTESTS=ON`` (the default if ``libxml2`` is available) in CMake.
Each module produces a separate unit test binary (:file:`{module}-test`) under
``bin/``, which can execute all the tests for that module.

The tests can be executed in a few different ways:

- Build the ``test`` target (e.g., ``make test``):
  This runs all the tests using CTest.  This includes also the regression
  tests if CMake has been told where to find them (regression tests are not
  discussed further on this page).
  If some of the tests fail, this only prints basic summary information (only
  a pass/fail status for each test binary or regression test class).
  You can execute the failing test binaries individually to get more
  information on the failure.
  Note that ``make test`` does not rebuild the test binaries if you have changed
  the source code, so you need to separately run ``make`` or ``make tests``.
  The latter only builds the test binaries and their dependencies.
- Build the ``check`` target (e.g., ``make check``):
  This behaves the same as the ``test`` target, with a few extensions:

  1. Test binaries are rebuilt if they are outdated before the tests are run.
  2. If a test fails, the output of the test binary is shown.
  3. If unit tests and/or regression tests are not available, a message is
     printed.

- Directly executing a test binary.  This provides the most useful output for
  diagnosing failures, and allows debugging test failures.  The output
  identifies the individual test(s) that fail, and shows the results of all
  failing assertions.  Some tests also add extra information to failing
  assertions to make it easier to identify the reason.  It is possible to
  control which tests are run using command line options.  Execute the binary
  with ``-h`` to get additional information.

When executed using CTest, the tests produce XML output in
``Testing/Temporary/``, containing the result of each test as well as failure
messages.  This XML is used by Jenkins for reporting the test status for
individual tests.  Note that if a test crashes or fails because of an assert or
a gmx_fatal() call, no XML is produced for the binary, and Jenkins does not
report anything for the test binary.  The actual error is only visible in the
console output.

Unit testing framework
----------------------

The tests are written using `Google Test`_, which provides a framework for
writing unit tests and compiling them into a test binary.  Most of the command
line options provided by the test binaries are implemented by Google Test.  See
the `Google Test Primer`_ for an introduction.
Some tests also use `Google Mock`_, which provides a framework for creating
mock implementations of C++ classes.  Both components are included in the
source tree under ``src/external/gmock-1.7.0/``, and are compiled as part of the
unit test build.

``src/testutils/`` contains |Gromacs|-specific shared test code.  This includes
a few parts:

- CMake macros for declaring test binaries.  These take care of providing the
  ``main()`` method for the test executables and initializing the other parts of
  the framework, so that the test code in modules can focus on the actual
  tests.  This is the only part of the framework that you need to know to be
  able to write simple tests: you can use ``gmx_add_unit_test()`` in CMake to
  create your test binary and start writing the actual tests right away.
  See ``src/testutils/TestMacros.cmake`` and existing CMake code for examples
  how to use them.

- Generic test fixtures and helper classes.  The C++ API is documented on
  `Doxygen page for testutils`__.  Functionality here includes
  locating test input files from the source directory and constructing
  temporary files, adding custom command line
  options to the test binary, some custom test assertions
  for better exception and floating-point handling, utilities
  for constructing command line argument arrays, and
  test fixtures for tests that need to test long strings for correctness
  and for tests that execute legacy code where
  `stdin` reading etc. cannot be easily mocked.

  __ doxygen-module-testutils_

- Some classes and functions to support the above.  This code is for internal
  use of the CMake machinery to build and set up the test binaries, and to
  customize Google Test to suit our environment.

- Simple framework for building tests that check the results against reference
  data that is generated by the same test code.  This can be used if it is not
  easy to verify the results of the code with C/C++ code alone, but manual
  inspection of the results is manageable.  The C++ API for this is also
  documented on the `Doxygen page for testutils`__.

  __ doxygen-module-testutils_

  When the test using the framework is first executed, ``-ref-data create`` can
  be passed on command line to create the reference data (also options
  starting with double dashes are accepted).
  On later executions, the tests read the reference data and fail if the
  results are not the same.  It is possible to update existing reference data
  with ``-ref-data update``.

  The reference data is stored in XML files under
  :file:`src/gromacs/{module}/tests/refdata/`.  This part of the framework
  depends on ``libxml2``.  For inspecting the reference data in a browser, there
  are XSLT stylesheets that transform the XML files into HTML.  Such custom
  transformations need to be written for each type of test if the output is
  not easy to check otherwise.  Because of security features in browsers, the
  transformations may not work for all browsers.  For the same reason, the
  XSLT files must be in the same folder as the XML files.  For cases where the
  XSLT files are shared between multiple modules, ``src/testutils/copy_xsl.sh``
  takes care to synchronize the files after a master copy is edited.

  The current reference data functionality is quite basic, but it can be extended
  if/when more control over, e.g., comparison tolerances is needed.

In addition to ``src/testutils/``, some of the module test directories may
provide reusable test code that is used in higher-level tests.  For example,
the ``src/gromacs/analysisdata/tests/`` provides test fixtures, a mock
implementation for gmx::AnalysisDataModuleInterface, and some helper classes
that are also used in ``src/gromacs/trajectoryanalysis/tests/``.
These cases are handled using CMake object libraries that are linked to all the
test binaries that need them.

Getting started with new tests
------------------------------

To start working with new tests, you should first read the `Google Test`_
documentation to get a basic understanding of the testing framework, and read
the above description to understand how the tests are organized in |Gromacs|.
It is not necessary to understand all the details, but an overall understanding
helps to get started.

Writing a basic test is straightforward, and you can look at existing tests for
examples.  The existing tests have a varying level of complexity, so here are
some pointers to find tests that use certain functionality:

- ``src/gromacs/utility/tests/stringutil.cpp`` contains very simple tests for
  functions.  These do
  not use any fancy functionality, only plain Google Test assertions.
  The only thing required for these tests is the ``TEST()`` macro and the block
  following it, plus headers required to make them compile.
- The same file contains also simple tests using the reference framework to
  check line wrapping (the tests for ``gmx::TextLineWrapper``).  The test fixture
  for these tests is in ``src/testutils/stringtest.h``/``.cpp``.  The string test
  fixture also demonstrates how to add a custom command line option to the
  test binary to influence the test execution.
- ``src/gromacs/selection/tests/`` contains more complex use of the
  reference framework.  This is the code the reference framework was
  originally written for.
  ``src/gromacs/selection/tests/selectioncollection.cpp`` is the main file to
  look at.
- For more complex tests that do not use the reference framework, but instead
  do more complex verification in code, you can look at
  ``src/gromacs/selection/tests/nbsearch.cpp``.
- For complex tests with mock-up classes and the reference framework, you can
  look at ``src/gromacs/analysisdata/tests/``.

Here are some things to keep in mind when working with the unit tests:

- Try to keep the execution time for the tests as short as possible, while
  covering the most important paths in the code under test.  Generally, tests
  should take seconds instead of minutes to run, so that no one needs to
  hesitate before running the tests after they have done some changes.
  Long-running tests should go somewhere else than in the unit test set.
  Note that Jenkins runs many of the tests under Valgrind, so heavy tests are
  going to slow down also that part of the verification.
- Try to produce useful messages when a test assertion fails.  The assertion
  message should tell what went wrong, with no need to run the *test itself*
  under a debugger (e.g., if the assertion is within a loop, and the loop
  index is relevant for understanding why the assertion fails, it should be
  included in the message).  Even better if even a user can understand what
  goes wrong, but the main audience for the messages is the developer who
  caused the test to fail.

.. _Google Test: http://code.google.com/p/googletest/
.. _Google Test Primer: http://code.google.com/p/googletest/wiki/V1_7_Primer
.. _Google Mock: http://code.google.com/p/googlemock/

.. include:: /fragments/doxygen-links.rst
