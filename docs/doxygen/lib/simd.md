Single-instruction Multiple-data (SIMD) coding {#page_simd}
==============================================

Coding with SIMD instructions
=============================

One important way for \Gromacs to achieve high performance is
to use modern hardware capabilities where a single assembly
instruction operates on multiple data units, essentially short
fixed-length vectors (usually 2, 4, 8, or 16 elements). This provides
a very efficient way for the CPU to increase floating-point
performance, but it is much less versatile than general purpose
registers. For this reason it is difficult for the compiler to
generate efficient SIMD code, so the user has to organize the
data in a way where it is possible to access as vectors, and
these vectors often need to be aligned on cache boundaries.

We have supported a number of different SIMD instruction sets in
the group kernels for ages, and it is now also present in the
verlet kernels and a few other places. However, with the increased
usage and several architectures with different capabilities we now
use a vendor-agnostic \Gromacs SIMD module, as documented in
\ref module_simd.

Design of the \Gromacs SIMD module
==================================

The functions in `src/gromacs/simd` are intended to be used for writing
architecture-independent SIMD intrinsics code. Rather than making assumptions
based on architecture, we have introduced a limited number of
predefined preprocessor macros that describe the capabilities of the
current implementation - these are the ones you need to check when
writing SIMD code. As you will see, the functionality exposed by
this module as typically a small subset of general SIMD implementations,
and in particular we do not even try to expose advanced shuffling or
permute operations, simply because we haven't been able to describe those
in a generic way that can be implemented efficiently regardless of the
hardware. However, the advantage of this approach is that it is straightforward
to extend with support for new simd instruction sets in the future,
and that will instantly speed up old code too.

To support the more complex stuff in the \Gromacs nonbonded kernels and
to make it possible to use SIMD intrinsics even for some parts of the code
where the data is not in SIMD-friendly layout, we have also added about 10
higher-level utility routines. These perform gather/scatter operations
on coordinate triplets, they load table data and aligned pairs (Lennard-Jones
parameters), and sum up the forces needed in the outer loop of the nonbonded
kernels. They are very straightforward to implement, but since they are
performance-critical we want to exploit all features of each architecture,
and for this reason they are part of the SIMD implementation.

Finally, for some architectures with large or very large SIMD width (e.g. AVX
with 8 elements in single precision, or AVX-512 with 16), the nonbonded
kernels can become inefficient. Since all such architectures presently known
(AVX, AVX2, MIC, AVX512) also provide extensive support for accessing
parts of the register, we optionally define a handful of routines to
perform load, store, and reduce operations based on half-SIMD-width data,
which can improve performance. It is only useful for wide implementations,
and it can safely be ignored first when porting to new platforms - they
are only needed for the so-called 2xnn SIMD kernels.

Unfortunately there is no standard for SIMD architectures. The available
features vary a lot, but we still need to use quite a few of them to
get the best performance possible. This means some features will only
be available on certain platforms, and it is critical that we do NOT make
to many assumptions about the storage formats, their size or SIMD width.
Just to give a few examples:

- On x86, double precision (64-bit) floating-point values always convert
  to 32-bit integers, while many other platforms use 64-bit, and some cannot
  use 32-bit integers at all. This means we cannot use a mask (boolean)
  derived from integer operations to select double-precision floating-point
  values, and it could get very complex for higher-level code if all these
  decisions were exposed. Instead, we want to keep integers 32-bit since
  all algorithms anyway need to work in single precision (w. 32-bit ints).
- IBM QPX uses 4-wide SIMD both for single and double precision. Integer
  support is highly limited, and the storage format means QPX does not
  use x86-style all-ones masks (which have different widths in single/double)
  but it uses the sign bit to denote the _false_ value. In particular, this
  means we cannot use the bit contents for any fancy mask operations.
- AVX1 only supports 4-wide 128-bit integer SIMD arithmetics, but the integer
  _conversions_ can still be done 8-wide which corresponds to the single
  precision floating-point width. Similarly, with AVX1 conversions between
  double-precision and integers use the 32-bit 4-wide 128bit registers where
  we can also do integer arithmetics. AVX2 adds proper arithmetics for
  8-wide integers. We would severely limit performance if we had to say
  that integer support was not present, so instead we stick to 32-bit ints
  but limit the operations we expose (and do shuffling internally).
- For SSE2 through SSE4.1, double precision is 2-wide, but when we convert
  to integers they will be put in the first two elements of a 4-wide integer
  type. This means we cannot assume that floating-point SIMD registers and
  corresponding integer registers (after conversion) have the same width.
- The 2-wide SIMD instructions on BlueGene/L and BlueGene/P cannot do any
  floating-point logical operations (and/andnot/or/xor) whatsoever, which
  can be a pain when implementing approximations for math functions.
- Since boolean values can have different width for float/double and the
  integers corresponding to float/double, we need to use separate boolean
  types for all these values and convert between them if we e.g. want to use
  result of an integer compare to select floating-point values.

While this might sound complicated, it is actually far easier than writing
separate SIMD code for 10 architectures in both single & double. The point
is not that you need to remember the limitations above, but it is critical
that you *never assume anything about the SIMD implementation*. We
typically implement SIMD support for a new architecture in days with this
new module, and the extensions required for verlet kernels
are also very straightforward (group kernels can be more complex, but those
are gradually on their way out). For the higher-level
code, the only important thing is to never _assume_ anything about the SIMD
architecture. Our general strategy in \Gromacs is to split the SIMD coding
in three levels:

<dl>
<dt>Base level generic SIMD</dt>
<dd>
The base level SIMD module (which we get by including `gromacs/simd/simd.h`
provides the API to define and manipulate SIMD datatypes. This will be enough
for lots of cases, and it is a huge advantage that there is roughly
parity between different architectures.
</dd>
<dt>Higher-level architecture-specific SIMD utility functions</dt>
<dd>
For some parts of the code this is not enough. In particular, both the
group and Verlet kernels do insane amounts of floating-point operations,
and since we spend 85-90% of the time in these kernels it is critical that
we can optimize them as much as possible. Here, our strategy is first to
define larger high-level functions that e.g. take a number of distances
and loads the table interactions for this interaction. This way we can
move this architecture-specific implementation to the SIMD module, and
both achieve a reasonably clean kernel but still optimize a lot. This
is what we have done for the approximately 10 functions for the nonbonded
kernels, to load tables and Lennard-Jones parameters, and to sum up the
forces in the outer loop. These functions have intentionally been given
names that describe what they do with the data, rather than what their
function is in \Gromacs. By looking at the documentation for these routines,
and the reference implementation, it should be quite straightforward to
implement them for a new architecture too.
</dd>
<dt>Half-SIMD-width architecture-specific utility functions</dt>
<dd>
As described earlier, as the SIMD width increases to 8 or more elements,
the nonbonded kernels can become inefficient due to the large j-particle
cluster size. Things will still work, but if an architecture supports
efficient access to partial SIMD registers (e.g. loading half the width),
we can use this to alter the balance between memory load/store operations
and floating-point arithmetic operations by processing either e.g. 4-by-4
or 2-by-8 interactions in one iteration. When \ref
GMX_SIMD_HAVE_HSIMD_UTIL_REAL is set, a handful of routines to
use this in the nonbonded kernels is present. Avoid using these routines
outside the nonbonded kernels since they are slightly more complex, and
is is not straightforward to determine which alternative provides the best
performance.
</dd>
<dt>Architecture-specific kernels (directories/files)</dt>
<dd>
No code outside the SIMD module implementation directories should try
to execute anything hardware specific. Note that this includes even checking
for what architecture the current SIMD implementation is - you should check
for features instead, so it will work with future ports too.
</dd>
</dl>

File organization
=================

The SIMD module uses a couple of different files:

<dl>
<dt>`gromacs/simd/simd.h`</dt>
<dd>
This is the top-level wrapper that you should always include first.
It will check the settings made at configuration time and include a
suitable low-level implementation (that can be either single, double,
or both). It also contains the routines for memory alignment, and
based on the current \Gromacs precision it will set aliases to 'real'
SIMD datatypes (see further down) so the implementations do not have
to care about \Gromacs-specific details. However, note that you might
not get all SIMD support you hoped for: If you compiled \Gromacs in
double precision but the hardware only supports single-precision SIMD
there will not be any SIMD routines for default \Gromacs 'real' precision.
There are \#defines you can use to check this, as described further down.
</dd>
<dt>`gromacs/simd/impl_reference/impl_reference.h`</dt>
<dd>
This is an example of a low-level implementation. You should never, ever,
work directly with these in higher-level code. The reference implementation
contains the documentation for all SIMD wrappers, though. This file will
in turn include other separate implementation files for single, double,
simd4, etc. Since we want to be able to run the low-level SIMD implementation
in simulators for new platforms, these files are intentionally not using
the rest of the GROMACS infrastructure, e.g. for asserts().
</dd>
<dt>`gromacs/simd/simd_math.h`</dt>
<dd>
SIMD math functions. All functions in this file have to be designed
so they work no matter whether the hardware supports integer SIMD, logical
operations on integer or floating-point SIMD, or arithmetic operations
on integers. However, a few routines check for defines and use faster
algorithms if these features are present.
</dd>
<dt>`gromacs/simd/vector_operations.h`</dt>
<dd>
This file contains a few rvec-related SIMD functions, e.g. to
calculate scalar products, norms, or cross products. They obviously
cannot operate on scalar \Gromacs rvec types, but use separate SIMD
variables for X,Y, and Z vector components.
</dd>
</dl>


SIMD datatypes
==============

The SIMD module handles the challenges mentioned in the introduction
by introducing a number of datatypes;
many of these might map to the same underlying SIMD types, but we need separate
types because some architectures use different registers e.g. for boolean
types.

Floating-point data
-------------------

<dl>
<dt>`#gmx::SimdReal`</dt>
<dd>
This is the SIMD-version of \Gromacs' real type,
which is set based on the CMake configuration and internally aliased
to one of the next two types.
</dd>
<dt>`#gmx::SimdFloat`</dt>
<dd>
This is always single-precision data, but it
might not be supported on all architectures.
</dd>
<dt>`gmx::SimdDouble`</dt>
<dd>
This is always double precision when available,
and in rare cases you might want to use a specific precision.
</dd>
</dl>

Integers corresponding to floating-point values
-----------------------------------------------

For these types, 'correspond' means that it is the integer type we
get when we convert data e.g. from single (or double) precision
floating-point SIMD variables. Those need to be different, since many
common implementations only use half as many elements for double as
for single SIMD variables, and then we only get half the number of
integers too.

<dl>
<dt>`#gmx::SimdInt32`</dt>
<dd>
This is used for integers when converting to/from \Gromacs default "real" type.
</dd>
<dt>`gmx::SimdFInt32`</dt>
<dd>
Integers obtained when converting from single precision, or intended to be
converted to single precision floating-point. These are normal integers
(not a special conversion type), but since some SIMD architectures such as
SSE or AVX use different registers for integer SIMD variables having the
same width as float and double, respectively, we need to separate these
two types of integers. The actual operations you perform on the are normal
ones such as addition or multiplication.
This will also be the widest integer data type if you want to do pure
integer SIMD operations, but that will not be supported on all platforms.
If the architecture does not support any SIMD integer type at all, this
will likely be defined from the floating-point SIMD type, without support
for any integer operations apart from load/store/convert.
</dd>
<dt>`gmx::SimdDInt32`</dt>
<dd>
Integers used when converting to/from double. See the preceding item
for a detailed explanation. On many architectures,
including all x86 ones, this will be a narrower type than `gmx::SimdFInt32`.
</dd>
</dl>

Note that all integer load/stores operations defined here load/store 32-bit
integers, even when the internal register storage might be 64-bit, and we
set the "width" of the SIMD implementation based on how many float/double/
integers we load/store - even if the internal width could be larger.

Boolean values
--------------

We need a separate boolean datatype for masks and comparison results, since
we cannot assume they are identical either to integers, floats or double -
some implementations use specific predicate registers for booleans.

<dl>
<dt>`#gmx::SimdBool`</dt>
<dd>
Results from boolean operations involving reals, and the booleans we use
to select between real values. The corresponding routines have suffix `B`,
like `gmx::simdOrB()`.
</dd>
<dt>`gmx::SimdFBool`</dt>
<dd>
Booleans specifically for single precision.
</dd>
<dt>`gmx::SimdDBool`</dt>
<dd>
Operations specifically on double.
</dd>
<dt>`#gmx::SimdIBool`</dt>
<dd>
Boolean operations on integers corresponding to real (see floating-point
descriptions above).
</dd>
<dt>`gmx::SimdFIBool`</dt>
<dd>
Booleans for integers corresponding to float.
</dd>
<dt>`gmx::SimdDIBool`</dt>
<dd>
Booleans for integers corresponding to double.
</dd>
</dl>

Note: You should NOT try to store and load boolean SIMD types to memory - that
is the whole reason why there are no store or load operations provided for
them. While it will be technically possible to achieve by defining objects
inside a structure and then doing a placement new with aligned memory, this
can be a very expensive operation on platforms where special single-bit
predicate registers are used to represent booleans. You will need to find
a more portable algorithm for your code instead.

The subset you should use in practice
-------------------------------------

If this seems daunting, in practice you should only need to use these types
when you start coding:

<dl>
<dt>`#gmx::SimdReal`</dt>
<dd>
Floating-point data.
</dd>
<dt>`#gmx::SimdBool`</dt>
<dd>
Booleans.
</dd>
<dt>`#gmx::SimdInt32`</dt>
<dd>
Integer data. Might not be supported, so you must check
the preprocessor macros described below.
</dd>
</dl>

Operations on these types will be defined to either float/double (or
corresponding integers) based on the current \Gromacs precision, so the
documentation is occasionally more detailed for the lower-level actual
implementation functions.

Note that it is critical for these types to be aligned in memory. This
should always be the case when you declare variables on the stack, but
unfortunately some compilers (at least clang-3.7 on OS X) appear to be
buggy when our SIMD datatypes are placed inside a structure. Somewhere
in the processes where this structure includes our class, which in turn
includes the actual SIMD datatype, the alignment appears to be lost.
Thus, even though the compiler will not warn you, until further notice
we need to avoid putting the SIMD datatypes into other structures. This
is particular severe when allocating memory on the heap, but it occurs
for stack structures/classes too.


SIMD4 implementation
--------------------

The above should be sufficient for code that works with the full SIMD width.
Unfortunately reality is not that simple. Some algorithms like lattice
summation need quartets of elements, so even when the SIMD width is >4 we
need width-4 SIMD if it is supported. The availability of SIMD4 is indicated
by \ref GMX_SIMD4_HAVE_FLOAT and \ref GMX_SIMD4_HAVE_DOUBLE. For now we only
support a small subset of SIMD operations for SIMD4. Because SIMD4 doesn't
scale with increasingly large SIMD width it should be avoided for all new
code and SIMD4N should be used instead.

SIMD4N implementation
---------------------

Some code, like lattice summation, has inner loops which are smaller
than the full SIMD width. In GROMACS algorithms 3 and 4 iterations are common
because of PME order and three dimensions. This makes 4 an important special
case. Vectorizing such loops efficiently requires to collapse the two
most inner loops and using e.g. one 8-wide SIMD vector for 2 outer
and 4 inner iterations or one 16-wide SIMD vector for 4 outer and 4 inner
iterations. For this SIMD4N functions are
provided. The availability of these function is indicated by
\ref GMX_SIMD_HAVE_4NSIMD_UTIL_FLOAT and
\ref GMX_SIMD_HAVE_4NSIMD_UTIL_DOUBLE.
These functions return the type alias Simd4NFloat / Simd4NDouble which is
either the normal SIMD type or the SIMD4 type and thus only supports
the operations the SIMD4 type supports.

Predefined SIMD preprocessor macros
===================================

Functionality-wise, we have a small set of core set of features that we
require to be present on all platforms, while more avanced features can be
used in the code when defines like e.g. \ref GMX_SIMD_HAVE_LOADU have the
value 1.

This is a summary of the currently available preprocessor defines that
you should use to check for support when using the corresponding features.
We first list the float/double/int defines set by the _implementation_; in
most cases you do not want to check directly for float/double defines, but
you should instead use the derived "real" defines set in this file - we list
those at the end below.

Preprocessor predefined macro defines set by the low-level implementation.
These only have the value 1 if they work for all datatypes;
\ref GMX_SIMD_HAVE_LOADU thus means we can load both float, double, and
integers from unaligned memory, and that the unaligned loads are available
for SIMD4 too.

<dl>
<dt>\ref GMX_SIMD</dt>
<dd>
Some sort of SIMD architecture is enabled.
</dd>
<dt>\ref GMX_SIMD_HAVE_FLOAT</dt>
<dd>
Single-precision instructions available.
</dd>
<dt>\ref GMX_SIMD_HAVE_DOUBLE</dt>
<dd>
Double-precision instructions available.
</dd>
<dt>\ref GMX_SIMD_HAVE_LOADU</dt>
<dd>
Load from unaligned memory available.
</dd>
<dt>\ref GMX_SIMD_HAVE_STOREU</dt>
<dd>
Store to unaligned memory available.
</dd>
<dt>\ref GMX_SIMD_HAVE_LOGICAL</dt>
<dd>
Support for and/andnot/or/xor on floating-point variables.
</dd>
<dt>\ref GMX_SIMD_HAVE_FMA</dt>
<dd>
Floating-point fused multiply-add.
Note: We provide emulated FMA instructions if you do not have FMA
support, but in that case you might be able to code it more efficient w/o FMA.
</dd>
<dt>\ref GMX_SIMD_HAVE_FINT32_EXTRACT</dt>
<dd>
Support for extracting integer SIMD elements from `gmx::SimdFInt32`.
</dd>
<dt>\ref GMX_SIMD_HAVE_FINT32_LOGICAL</dt>
<dd>
Bitwise shifts on `gmx::SimdFInt32`.
</dd>
<dt>\ref GMX_SIMD_HAVE_FINT32_ARITHMETICS</dt>
<dd>
Arithmetic ops for `gmx::SimdFInt32`.
</dd>
<dt>\ref GMX_SIMD_HAVE_DINT32_EXTRACT</dt>
<dd>
Support for extracting integer SIMD elements from `gmx::SimdDInt32`.
</dd>
<dt>\ref GMX_SIMD_HAVE_DINT32_LOGICAL</dt>
<dd>
Bitwise shifts on `gmx::SimdDInt32`.
</dd>
<dt>\ref GMX_SIMD_HAVE_DINT32_ARITHMETICS</dt>
<dd>
Arithmetic ops for `gmx::SimdDInt32`.
</dd>
<dt>\ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT</dt>
<dd>
Half-SIMD-width nonbonded kernel utilities available for float SIMD.
</dd>
<dt>\ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE</dt>
<dd>
Half-SIMD-width nonbonded kernel utilities available for double SIMD.
</dd>
<dt>\ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT</dt>
<dd>
Can load pairs of unaligned floats from simd offsets (meant for linear tables).
</dd>
<dt>\ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE</dt>
<dd>
Can load pairs of unaligned doubles from simd offsets (meant for linear tables).
</dd>
</dl>

There are also two macros specific to SIMD4: \ref GMX_SIMD4_HAVE_FLOAT is set
if we can use SIMD4 in single precision, and \ref GMX_SIMD4_HAVE_DOUBLE
similarly denotes support for a double-precision SIMD4 implementation. For
generic properties (e.g. whether SIMD4 FMA is supported), you should check
the normal SIMD macros above.

Implementation properties
-------------------------

Higher-level code can use these macros to find information about the implementation,
for instance what the SIMD width is:

<dl>
<dt>\ref GMX_SIMD_FLOAT_WIDTH</dt>
<dd>
Number of elements in `gmx::SimdFloat`, and practical width of `gmx::SimdFInt32`.
</dd>
<dt>\ref GMX_SIMD_DOUBLE_WIDTH</dt>
<dd>
Number of elements in `gmx::SimdDouble`, and practical width of `gmx::SimdDInt32`</dd>
<dt>\ref GMX_SIMD_RSQRT_BITS</dt>
<dd>
Accuracy (bits) of 1/sqrt(x) lookup step.
</dd>
<dt>\ref GMX_SIMD_RCP_BITS</dt>
<dd>
Accuracy (bits) of 1/x lookup step.
</dd>
</dl>

After including the low-level architecture-specific implementation, this
header sets the following derived defines based on the current precision;
these are the ones you should check for unless you absolutely want to dig
deep into the explicit single/double precision implementations:

<dl>
<dt>\ref GMX_SIMD_HAVE_REAL</dt>
<dd>
Set to either \ref GMX_SIMD_HAVE_FLOAT or \ref GMX_SIMD_HAVE_DOUBLE
</dd>
<dt>\ref GMX_SIMD4_HAVE_REAL</dt>
<dd>
Set to either \ref GMX_SIMD4_HAVE_FLOAT or \ref GMX_SIMD4_HAVE_DOUBLE
</dd>
<dt>\ref GMX_SIMD_REAL_WIDTH</dt>
<dd>
Set to either \ref GMX_SIMD_FLOAT_WIDTH or \ref GMX_SIMD_DOUBLE_WIDTH
</dd>
<dt>\ref GMX_SIMD_HAVE_INT32_EXTRACT</dt>
<dd>
Set to either \ref GMX_SIMD_HAVE_FINT32_EXTRACT or \ref GMX_SIMD_HAVE_DINT32_EXTRACT
</dd>
<dt>\ref GMX_SIMD_HAVE_INT32_LOGICAL</dt>
<dd>
Set to either \ref GMX_SIMD_HAVE_FINT32_LOGICAL or \ref GMX_SIMD_HAVE_DINT32_LOGICAL
</dd>
<dt>\ref GMX_SIMD_HAVE_INT32_ARITHMETICS</dt>
<dd>
Set to either \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS or \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS
</dd>
<dt>\ref GMX_SIMD_HAVE_HSIMD_UTIL_REAL</dt>
<dd>
Set to either \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT or \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE
</dd>
<dt>\ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL</dt>
<dd>
Set to either \ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT or \ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE
</dd>
</dl>

For convenience we also define \ref GMX_SIMD4_WIDTH to 4. This will never vary,
but using it helps you make it clear that a loop or array refers to the
SIMD4 width rather than some other '4'.

While all these defines are available to specify the features of the
hardware, we would strongly recommend that you do NOT sprinkle your code
with defines - if nothing else it will be a debug nightmare. Instead you can
write a slower generic SIMD function that works everywhere, and then override
this with faster architecture-specific versions for some implementations. The
recommended way to do that is to add a define around the generic function
that skips it if the name is already defined. The actual implementations in
the lowest-level files are typically defined to an architecture-specific name
(such as `simdSinCosD_Sse2`) so we can override it (e.g. in SSE4) by
simply undefining and setting a new definition. Still, this is an
implementation detail you won't have to worry about until you start writing
support for a new SIMD architecture.


Function naming
---------------

We rely on C++ overloading, so the name of a function is usually identical
regardless of what datatype it operates on. There are a few exceptions to this
for functions that do not take arguments but only return a value, e.g. setZero(),
since overloading only works if the formal parameters are different. To solve this,
we use different low-level function names in these cases, but then create proxy
objects in the high-level `gromacs/simd/simd.h` so that you can still get the
functionality by simply writing setZero() in the code.

Automated checking
------------------

Having fallback implementations when SIMD is not supported can be a
performance problem if the code does not correctly include
`gromacs/simd/simd.h`, particularly after refactoring.
`make check-source` checks the whole code for the use of symbols defined
in `gromacs/simd/simd.h` and requires that files using those symbols
do the correct include. Similar checking is done for higher-level
SIMD-management headers, e.g. `gromacs/ewald/pme-simd.h`.


The SIMD math library
=====================

In addition to the low-level SIMD instructions, \Gromacs comes with a fairly
extensive SIMD math library in `gromacs/simd/simd_math.h` to support various
mathematical functions. The functions are available both in single and
double precision (overloaded on the usual math function names), and we also
provide a special version of functions that use double precision arguments,
but that only evaluate the result to single precision accuracy. This is
useful when you don’t need highly accurate results, but you want to avoid
the overhead of doing multiple single/double conversions, or if the hardware
architecture only provides a double precision SIMD implementation.

For a few functions such as the square root and exponential that are
performance-critical, we provide additional tempate parameters where the
default choice is to execute the normal function version, but it is also
possible to choose an unsafe execution path that completely bypass all
argument checking. Make absolutely sure your arguments always fulfil the
restrictions listed in the documentation of such a function before using it,
and it might even be a good idea to add a note before each call to an unsafe
function justifying why that flavor is fine to use here.
