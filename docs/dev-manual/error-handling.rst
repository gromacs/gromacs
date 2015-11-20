.. _error handling:

Error handling
==============

To make |Gromacs| behave like a proper library, we need to change the
way errors etc. are handled. Basically, the library should not print
out anything to stdio/stderr unless it is part of the API
specification, and even then, there should be a way for the user to
suppress the output. Also, the library should normally not terminate
the program without the user having control over this. There are
different types of errors, which also affects the handling. Different
cases are discussed separately below, split by the way they are
handled. These guidelines are starting to take their final form,
although details may still change.

* For programming errors, i.e., errors that should never occur if the
  program is correctly written, it's acceptable to assert and
  terminate the program. This applies to both errors in the library
  and errors in user code or user input that calls the library.
  Older code tends to still use ``assert()`` calls, but new
  code should prefer more expressive functionality such as
  ``GMX_RELEASE_ASSERT()``. This version of the macro will result
  in asserts that are still present when the build type is Release,
  which is what we want by default. In performance-sensitive parts
  of the code, it is acceptable to rather use ``GMX_ASSERT()`` to
  avoid the performance penalty of a branch when the code is compiled
  for production use. By default, Jenkins builds the RelWithAssert
  build type.
* For some errors it might be feasible to recover gracefully and
  continue execution. In this case, your APIs should be defined
  so that the API-user/programmer does not have to check separately
  whether the problem was due to a programming error, but it's
  better to e.g. use exceptions for recoverable errors and
  asserts for programming errors.
* Exceptions should only be used for unexpected errors, e.g., out of
  memory or file system IO errors. As a general guideline, incorrect
  user input should not produce an untrapped exception resulting
  in execution termination telling the user an exception occured.
  Instead, you should catch exceptions in an earlier stack frame,
  make a suitable decision about diagnostic messages, and then
  decide how execution should be terminated.
* There is a global list of possible exceptions in
  ``src/gromacs/utility/exceptions.h``, and the library should throw
  one of these when it fails, possibly providing a more detailed
  description of the reason for the failure. The types of exceptions
  can be extended, and currently include:

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
  or join threads with OpenMP, since OpenMP cannot make guarantees about
  whether exceptions are caught or if the program will crash.
  Currently we catch all exceptions before we leave an OpenMP threaded region.
  If you throw an exception, make sure that it is caught and handled appropriately
  in the same thread/OpenMP section.
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

For coding guidelines to make this all work, see :ref:`implementing exceptions`.
