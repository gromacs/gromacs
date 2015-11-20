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
  and errors in user code that calls the library.
* If it's feasible to recover graciously and return an error code,
  this is also possible (see the next point), but the API should in
  such cases be designed such that the users do not separately need to
  check whether they made a programming error.
* For code that is not performance-sensitive, consider using asserts
  that remain in the code even in release builds, in particular if the
  assert could fire as a result of incorrect user code.
* Performance-sensitive code should not have a lot of error-checking
  at run time in release mode, but for debugging, it is useful to be
  able to pinpoint where things start to go wrong instead of observing
  a crash that may occur much later. Using assertions is useful here.
* There is a global list of possible exceptions in
  ``src/gromacs/utility/exceptions.h``, and the library should return
  one of these when it fails, possibly providing a more detailed
  description of the reason for the failure. The exceptions include:

  - Out of memory (e.g. ``std::bad_alloc``)

  - File I/O error (e.g. not found)

  - Invalid user input (could not be understood)

  - Inconsistent user input (parsed correctly, but has internal conflicts)

  - Simulation instability

  - Invalid API call/value/internal error (an assertion might also be used in such cases)

  - In the internals of a module called from code that is not
    exception safe, you can use exceptions for error handling, but
    avoid propagating them to caller code.

  - Exceptions should only be used for unexpected errors, e.g., out of
    memory or file system IO errors. As a general guideline, incorrect
    user input should not result in the user being told that an
    exception occured, though an exception can be used to transfer
    control to an earlier stack frame that can print a better
    diagnostic message.

* Avoid exceptions in threaded code, but if you throw one, make sure
  that it will always gets caught in the same thread/OpenMP section.
* There are also cases where a library routine wants to report a
  warning or a non-fatal error, but is still able to continue
  processing. For example, what grompp does now with notes, warnings,
  and errors. If so, it should maintain an internal data structure of
  such issues and then report all of them at one time.
* A function should not "fail" as part of its normal operation,
  however doing nothing can be considered normal operation, e.g. a
  function that accesses data, that by design should also be callable
  when no data is available. If the failure is not normal, then it
  should throw an exception.

For coding guidelines for making this all work, see :ref:`implementing exceptions`.
