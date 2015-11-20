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
  header file which is needed (or likely will be nedded) by CUDA should not use C++11.
* We should be able to use virtually all C++ features outside of the header files
  required by CUDA code (and OpenCL kernels), since we have gradually moved to
  compilers that have full support for C++11.

C++ Standard Library
--------------------

|Gromacs| code must support the lowest common denominator of C++11 standard library
features available on supported platforms.
Some modern features are useful enough to warrant back-porting.
Consistent and forward-compatible headers are provided in ``src/gromacs/compat/``
as described in the `Library documentation <../doxygen/html-lib/group__group__compatibility.xhtml>`_

General considerations
^^^^^^^^^^^^^^^^^^^^^^
As a baseline, |Gromacs| follows the C++ Core Guidelines |linkref1|, unless
our own more specific guidelines below say otherwise. We tend to be more restrictive
in some areas, both because we depend on the code compiling with a lot of different
C++ compilers, and because we want to increase readability. However, |Gromacs| is an
advanced projects in constant development, and as our needs evolve we will both
relax and tighten many of these points. Some of these changes happen naturally as
part of agreements in code review, while major parts where we don't agree should be
pushed to a redmine thread. Large changes should be suggested early in the development
cycle for each release so we avoid being hit by last-minute compiler bugs just before
a release.

* Use namespaces. Everything in ``libgromacs`` should be in a ``gmx``
  namespace. Don't use using in headers except possibly for aliasing
  some commonly-used names, and avoid file-level blanket ``using
  namespace gmx`` and similar. If only a small number of ``gmx``
  namespace symbols needed in a not-yet-updated file, consider
  importing just those symbols. See also |linkref2|.
* Use STL, but do not use iostreams outside of the unit tests. iostreams can have
  a negative impact on performance compared to other forms 
  of string streams, depending on the use case. Also, they don't always
  play well with using C ``stdio`` routines at the same time, which
  are used extensively in the current code. However, since Google tests
  rely on iostreams, you should use it in the unit test code.
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
* Avoid using default function arguments. They can lead to the code
  being less readable than without (see |linkref3|). If you think that your specific
  case improves readability (see |linkref4|), you can justify their use.
* Don't overload operators before thorough consideration whether it
  really is the best thing to do. Never overload ``&&``, ``||``, or
  the comma operator, because it's impossible to keep their original
  behavior with respect to evaluation order.
* Try to avoid complex templates, complex template specialization or
  techniques like SFINAE as much as possible. If nothing else, they
  can make the code more difficult to understand.
* Don't use multiple inheritance. Inheriting from multiple pure
  interfaces is OK, as long as at most one base class (which should be
  the first base class) has any code. Please also refer to the
  explanation |linkref5| and |linkref6|.
* Don't write excessively deep inheritance graphs. Try to not inherit
  implementation just to save a bit of coding; follow the principle
  "inherit to be reused, not to reuse." Also, you should not
  mix implementation and interface inheritance. For explanation please
  see |linkref7|.
* Don't include unnecessary headers.
* Make liberal use of assertions to help document your intentions (but
  prefer to write the code such that no assertion is necessary).
* Prefer ``GMX_ASSERT()`` and ``GMX_RELEASE_ASSERT()`` to naked
  ``assert()`` because the former permit you to add descriptive text.
* Use gmx::Mutex rather than pthreads, std or raw thread-MPI mutexes.
* Use proper enums for variable whose type can only contain one of a
  limited set of values. C++ is much better than C in catching errors
  in such code. Ideally, all enums should be typed enums, please
  see |linkref8|. 
* When writing a new class, think whether it will be necessary to make
  copies of that class. If not, declare the copy constructor and the
  assignment operator as private and don't define them, making any
  attempt to copy objects of that class fail. If you allow copies,
  either provide the copy constructor and the assignment operator, or
  write a clear comment that the compiler-generated ones will do (and
  make sure that they do what you
  want). ``src/gromacs/utility/classhelpers.h`` has some convenience
  macros for doing this well.
  Starting from c++11, you can also use deleted functions in this case.
* Declare all constructors with one parameter as explicit unless you
  really know what you are doing. Otherwise, they can be used for
  implicit type conversions, which can make the code difficult to
  understand, or even hide bugs that would be otherwise reported by
  the compiler. For the same reason, don't declare operators for
  converting your classes to other types without thorough
  consideration. For an explanation, please see |linkref9|.
* Write const-correct code (no ``const_cast`` unless absolutely
  necessary).
* Avoid using RTTI (run-time type information, in practice
  ``dynamic_cast`` and ``typeid``) unless you really need it. The cost
  of RTTI is very high, both in binary size (which you always
  pay if you compile with it) and in execution time (which you pay
  only if you use it). If your problem seems to require RTTI, think
  about whether there would be an alternative design that
  wouldn't. Such alternative designs are often better.
* Don't depend on compiler metadata propagation. struct elements
  and captured lambda parameters tend to have ``restrict`` and
  alignment qualifiers discarded by compilers, so when you later
  define an instance of that structure or allocate memory to
  hold it, the data member might not be aligned at all.
* Plan for code that runs in compute-sensitive kernels to have useful
  data layout for re-use, alignment for SIMD memory operations
* Recognize that some parts of the code have different requirements -
  compute kernels, mdrun setup code, high-level MD-loop code,
  simulation setup tools, and analysis tools have different needs, and
  the trade-off point between correctness vs reviewer time vs
  developer time vs compile time vs run time will differ.


.. |linkref1| replace:: `c++ guidelines <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines>`__
.. |linkref2| replace:: `here <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#sf7-dont-write-using-namespace-in-a-header-file>`__
.. |linkref3| replace:: `here <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#i23-keep-the-number-of-function-arguments-low>`__
.. |linkref4| replace:: `here <https://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#f51-where-there-is-a-choice-prefer-default-arguments-over-overloading>`__
.. |linkref5| replace:: `here <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#c135-use-multiple-inheritance-to-represent-multiple-distinct-interfaces>`__
.. |linkref6| replace:: `here <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#c136-use-multiple-inheritance-to-represent-the-union-of-implementation-attributes>`__
.. |linkref7| replace:: `here <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#c129-when-designing-a-class-hierarchy-distinguish-between-implementation-inheritance-and-interface-inheritance>`__
.. |linkref8| replace:: `here <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Renum-class>`__
.. |linkref9| replace:: `here <http://isocpp.github.io/CppCoreGuidelines/CppCoreGuidelines#Rc-explicit>`__

.. _implementing exceptions:

Implementing exceptions for error handling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
See :ref:`error handling` for the approach to handling run-time
errors, ie. use exceptions.

* Write exception-safe code. All new code has to offer at least the
  basic or nothrow guarantee to make this feasible.
* Use std (or custom) containers wherever possible.
* Use smart pointers for memory management. By default, use
  ``std::unique_ptr`` and ``gmx::unique_cptr`` in assocation with any
  necessary raw ``new`` or ``snew`` calls. ``std::shared_ptr`` can be
  used wherever responsibility for lifetime must be shared.
  Never use ``malloc``.
* Use RAII for managing resources (memory, mutexes, file handles, ...).
* It is preferable to avoid calling a function which might throw an
  exception from a legacy function which is not exception safe. However,
  we make the practical exception to permit the use of features such
  as ``std::vector`` and ``std::string`` that could throw
  ``std::bad_alloc`` when out of memory. In particular, |Gromacs| has
  a lot of old C-style memory handling that checking tools continue
  to issue valid warnings about as the tools acquire more
  functionality, and fixing these with old constructs is an
  inefficient use of developer time.
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
* Please strongly consider a comment repeating the preprocessor condition at the end
  of the region, if a lengthy region is neccessary and benefits from
  that. For long regions this greatly helps in understanding 
  and debugging the code.
