.. _error handling:

Error handling
==============

To make |Gromacs| behave like a proper library, we need to handle
errors in a consistent and predictable way. In this section, "user"
refers to the end user of |Gromacs| whether via some command-line
tool, or a workflow, or a call to a public API. There are different
types of errors, and the handling reflects this. This section is a work
in progress, particularly as the broader C++ community is a long way
from consensus in these areas.

Brief summary on which method to use
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

More detailed rules and rationale are written below, but in short, when a reason
exists that code is unable to do its job:

* If the reason can be checked at compile-time, then use ``static_assert``.
* If the reason is normal in context, then express that in the types used
  (e.g. return ``std::optional``) and document that this is normal.
* If the reason is that an internal invariant or pre-condition is violated (e.g.
  unexpected null pointer passed) on a hot code path, then use ``GMX_ASSERT``.
* Otherwise, if the reason is that an internal invariant or pre-condition is violated
  then use ``GMX_RELEASE_ASSERT``.
* Otherwise, (typically an error returned from system call or GPU SDK, bad user
  input), then use ``GMX_THROW``.

Guiding principles
^^^^^^^^^^^^^^^^^^

* |Gromacs| should adopt approaches that have achieved consensus
  elsewhere, e.g. in the `C++ Core Guidelines
  <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines>`_. In
  particular, be guided by `its section on error handling
  <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#S-errors>`_
* The library should not print out anything to stdio/stderr unless it
  is part of the API specification, and even then, there should be a
  way for the user to suppress or redirect the output.
* The library should normally not terminate the program without
  the user having control over this.
* Design interfaces of functions, classes, modules, and libraries so
  that values passed at run time are valid. Pass const references or
  ``not_null`` pointers rather than raw pointers. Return objects where
  possible. Use e.g. class enums for the type of passed
  values. Consider such enums as template parameters, rather than
  passing run-time values. Refactor existing interfaces to improve
  such aspects when starting new work in an area.
* Check user input at API boundaries and establish invariants as soon
  as possible, e.g. by expressing the user's choice in the type
  system. These form the pre-conditions that error handling will rely
  on.
* Use assertions to validate invariants and pre-conditions. There is value in
  using a different technique for checking such violations in order to make 
  the reason for the check clear to the maintainer.

Specific rules
^^^^^^^^^^^^^^

* Use ``static_assert`` wherever possible to detect errors at compile
  time.
* Throw *exceptions* to indicate that a function can't do its assigned
  task, per the `C++ Core Guidelines E.2
  <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Re-throw>`_.
  In particular, constructors should throw when they cannot construct
  a valid object, per `C++ Core Guidelines C.42
  <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Re-invariant>`_.
  However, recognize that in some cases the underlying reason is that
  some other component has not set up the correct pre-condition, and
  such cases should be handled with assertions (see below).
* At API boundaries, the assigned task of some code will be to
  validate the input, and that code should express failure to validate
  by throwing.
* Many programming errors violate pre-conditions of other
  functions. Until there is language support for contracts, the best
  that can be done is to check these with *assertions*. Note that only
  one component should have the responsibility for validating any
  particular input from the user, and other components should rely
  upon that validation in their pre-conditions.
* When asserting, use ``GMX_RELEASE_ASSERT`` by default. This macro
  will run its check in all build configurations, including
  ``Release``.
* When asserting in cases where the code is called in an inner loop of
  e.g. the MD step, ``GMX_ASSERT`` can be used. This macro will run
  its check only when ``NDEBUG`` is not defined, including the
  ``RelWithAssert`` build configuration (which is the default build
  type used in CI).
* It can be appropriate to provide both checked and unchecked
  interfaces, as ``std::vector`` does with ``at()`` and
  ``operator[]``, respectively. Note that even the latter is checked
  if you build e.g. ``libstdc++`` in the right configuration!
* When calling low-level APIs (including C and C++ standard
  libraries, GPU SDKs) always check for success/failure. Generally the
  correct thing to do upon failure will be to throw, perhaps including
  a descriptive string obtained from an error code with another API
  call.
* Do catch exceptions from lower-level components
  memory or file system IO errors. As a general guideline, incorrect
  user input should not produce an untrapped exception resulting
  in execution termination telling the user an exception occured.
  Instead, you should catch exceptions in an earlier stack frame,
  make a suitable decision about diagnostic messages, and then
  decide whether execution should be terminated (if that is in the
  scope of the code making the decision) and, if so, how to terminate.
* There is a global list of possible exceptions in
  ``api/legacy/include/gromacs/utility/exceptions.h``, and the library
  should throw one of these when it fails, possibly providing a more
  detailed description of the reason for the failure. The types of
  exceptions can be extended, and currently include:

  - Out of memory (e.g. ``std::bad_alloc``)

  - File I/O error (e.g. not found)

  - Invalid user input (could not be understood)

  - Inconsistent user input (parsed correctly, but has internal conflicts)

  - Simulation instability

  - Invalid API call/value/internal error (an assertion might also be used in such cases)

  - In the internals of a module called from code that is not
    exception safe, you can use exceptions for error handling, but
    avoid propagating them to caller code.

* Avoid using exceptions to propagate errors across regions that start
  or join threads with OpenMP, since OpenMP cannot make guarantees
  about whether exceptions are caught or if the program will crash.
  Currently we catch all exceptions before we leave an OpenMP threaded
  region.  If you throw an exception, make sure that it is caught and
  handled appropriately in the same thread/OpenMP section.
* Avoid using exceptions to propagate errors within regions where
  non-blocking API calls (e.g. to MPI or GPU SDKs) have been made,
  because the possible advantage of catching at a higher level and
  continuing execution is absent when the partner in the API call
  may be left blocked.
* There are also cases where a library routine wants to report a
  warning or a non-fatal error, but is still able to continue
  processing. In this case you should try to collect all issues and
  report and report them (similar to what grompp does with notes, warnings
  and errors) instead of just returning the first error. It is irritating
  to users if they fix the reported error, but then they keep getting
  a new error message every time the rerun the program.
* A function should not fail as part of its normal operation.
  However, doing nothing can be considered normal operation. A function
  accessing data should typically also be callable when no such data is
  available, but still return through normal means. If the failure is not
  normal, it is OK to rather throw an exception.
* Error handling with ``gmx_fatal``, ``gmx_warning``, ``gmx_incons``,
  ``gmx_comm`` etc.  is deprecated and should generally be refactored
  to throw or assert according to the above guidelines.
* There is currently no attempt made to check for error states on
  other MPI ranks during the simulation and provide a coordinated
  recovery. However setup code should do such checks routinely.
* We use ``GMX_RELEASE_ASSERT`` and ``GMX_ASSERT`` rather
  than ``assert`` to ensure that non-immediate strings can be
  used to describe the problem when the error is reported.
  This is particularly useful when troubleshooting issues where
  missing test coverage leads users to uncover such errors.

  
For coding guidelines to make this all work, see :ref:`implementing exceptions`.
