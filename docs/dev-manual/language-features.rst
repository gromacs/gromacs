Allowed language features
=========================

Most of these are not strict rules, but you should have a very good
reason for deviating from them.

Portability considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^

|Gromacs| uses C99 for C files and C++11 for C++ files. 
C++ has a lot of features, but to keep the source code maintainable and easy to read, 
we will avoid using some of them in |Gromacs| code. The basic principle is to keep things 
as simple as possible.
For compatiblity, certain work-arounds are required because not all compilers support 
these standards fully.

* MSVC supports only a subset of C99 and work-arounds are required in those cases.
* Before 7.0 (partial support in 6.5) CUDA didn't support C++11. Therefore any
  header file which is needed (or likly will be nedded) by CUDA should not use C++11.
* C++11 features which are not widely implemented (including in MSVC 2013 and GCC 4.6)
  should not be used.

General considerations
^^^^^^^^^^^^^^^^^^^^^^

* Use namespaces. Everything in ``libgromacs`` should be in a ``gmx``
  namespace. Don't use using in headers except possibly for aliasing
  some commonly-used names, and avoid file-level blanket ``using
  namespace gmx`` and similar. If only a small number of ``gmx``
  namespace symbols needed in a not-yet-updated file, consider
  importing just those symbols.
* Use STL, but don't use iostreams. iostreams are slow, and don't
  always play well with using C ``stdio`` routines at the same time, which
  are used extensively in the current code.
* Don't use non-const references as function parameters. They make it
  impossible to tell whether a variable passed as a parameter may
  change as a result of a function call without looking up the
  prototype.
* Don't use C-style casts; use ``const_cast``, ``static_cast`` or
  ``reinterpret_cast as appropriate``. See the point on RTTI for
  ``dynamic_cast``.
* Avoid overloading functions unless all variants really do the same
  thing, just with different types. Instead, consider making the
  function names more descriptive.
* Avoid using default function arguments.
* Don't overload operators before thorough consideration whether it
  really is the best thing to do. Never overload ``&&``, ``||``, or
  the comma operator, because it's impossible to keep their original
  behavior with respect to evaluation order.
* Don't use complex templates, complex template specialization or
  techniques like SFINAE. If nothing else, they can make the code more
  difficult to understand.
* Avoid combining complicated features in the same translation unit,
  because this multiples the risk of non-portability.
* Don't use multiple inheritance. Inheriting from multiple pure
  interfaces is OK, as long as at most one base class (which should be
  the first base class) has any code.
* Don't write excessively deep inheritance graphs. Try to not inherit
  implementation just to save a bit of coding; follow the principle
  "inherit to be reused, not to reuse."
* Don't plan to use anything from Boost without prior discussion
  (e.g. on ``gmx-developers`` mailing list). Boost is a nice library,
  but we don't want to depend on all the template magic inside, which
  may not work on more exotic compilers. Excessive template use also
  slows down compilation significantly.
* Don't include unnecessary headers (see :ref:`library structure`).
* Make liberal use of assertions to help document your intentions (but
  prefer to write the code such that no assertion is necessary).
* Prefer ``GMX_ASSERT()`` and ``GMX_RELEASE_ASSERT()`` to naked
  ``assert()`` because the former permit you to add descriptive text.
* tMPI provides basic mutexes and RAII locks for them; can be extended
  when the need arises (see
  http://redmine.gromacs.org/issues/948). TODO update code and this
  recommendation to use std::mutex.
* Use proper enums for variable whose type can only contain one of a
  limited set of values. C++ is much better than C in catching errors
  in such code.
* When writing a new class, think whether it will be necessary to make
  copies of that class. If not, declare the copy constructor and the
  assignment operator as private and don't define them, making any
  attempt to copy objects of that class fail. If you allow copies,
  either provide the copy constructor and the assignment operator, or
  write a clear comment that the compiler-generated ones will do (and
  make sure that they do what you
  want). ``src/gromacs/utility/classhelpers.h`` has some convenience
  macros for doing this well.
* Declare all constructors with one parameter as explicit unless you
  really know what you are doing. Otherwise, they can be used for
  implicit type conversions, which can make the code difficult to
  understand, or even hide bugs that would be otherwise reported by
  the compiler. For the same reason, don't declare operators for
  converting your classes to other types without thorough
  consideration.
* Write const-correct code (no ``const_cast`` unless absolutely
  necessary).
* Avoid using RTTI (run-time type information, in practice
  ``dynamic_cast`` and ``typeid``) unless you really need it. The cost
  of RTTI is not very high, neither in binary size (which you always
  pay if you compile with it) nor in execution time (which you pay
  only if you use it). If your problem seems to require RTTI, think
  about whether there would be an alternative design that
  wouldn't. Such alternative designs are often better.
* Don't depend on compiler metadata propagation (e.g. struct elements
  and captured lambda parameters tend to have ``restrict`` and
  alignment qualifiers discarded by compilers)
* Plan for code that runs in compute-sensitive kernels to have useful
  data layout for re-use, alignment for SIMD memory ops
* Recognize that some parts of the code have different requirements -
  compute kernels, mdrun setup code, high-level MD-loop code,
  simulation setup tools, and analysis tools have different needs, and
  the trade-off point between correctness vs reviewer time vs
  developer time vs compile time vs run time will differ.

.. _implementing exceptions:
   
Implementing exceptions for error handling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See :ref:`error handling` for the approach to handling run-time
errors, ie. use exceptions.

* Write exception-safe code. All new code has to offer at least the
  basic or nothrow guarantee to make this feasible.
* Use RAII for managing resources (memory, mutexes, file handles, ...).
* gmx::File provides some RAII support, but may need further work or a
  rewrite as part of http://redmine.gromacs.org/issues/950.
* Don't use ``malloc``, and try to limit the use of ``snew`` and
  ``new`` as well. Use container classes when appropriate instead of
  managing the memory everywhere manually.
* Use smart pointers for memory management. Generally use
  ``std::unique_ptr`` or ``std::shared_ptr``.
* It is not legal to call a function which might throw an exception
  from a legacy function which is not exception safe. This includes
  things like ``std::vector`` and ``std::string``, which can throw
  ``std::bad_alloc`` when out of memory!
* Functions / methods should be commented whether they are exception
  safe, whether they might throw an exception (even indirectly), and
  if so, which exception(s) they might throw.

Preprocessor considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^^
* Don't use preprocessor defines for things other than directly
  related to configuring the build. Use templates or inline functions
  to generate code, and enums or const variables for constants.
* Preprocessing variables used for configuring the build should be
  organized so that a valid value is always defined, i.e. we never
  test whether one of our preprocessor variables is defined, rather we
  test what value it has. This is much more robust under maintance,
  because a compiler can tell you that the variable is undefined.
* Avoid code with lengthy segments whose compilation depends on #if
  (or worse, #ifdef).
* Prefer to organize the definition of a const variable at the top of
  the source code file, and use that in the code.  This helps keep all
  compilation paths built in all configurations, which reduces the
  incidence of silent bugs.
* Indent nested preprocessor conditions if nesting is necessary and
  the result looks clearer than without indenting.
* Consider a comment repeating the preprocessor condition at the end
  of the region, if a lengthy region is neccessary and benefits from
  that.
