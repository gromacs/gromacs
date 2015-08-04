Single-instruction Multiple-data (SIMD) coding {#page_simd}
==============================================

Coding with SIMD instructions
=============================

One important way for \Gromacs to achieve high performance is
to use modern hardware capabilities where a single assembly
instruction operates on multiple data units, essentially short
fixed-length vectors (usually 2,4,8, or 16 elements). This provides
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

The macros in `src/gromacs/simd` are intended to be used for writing
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
<dt>Larger architecture-specific SIMD functions</dt>
<dd>
For some parts of the code this is not enough. In particular, both the
group and Verlet kernels do insane amounts of floating-point operations,
and since we spend 85-90% of the time in these kernels it is critical that
we can optimize them as much as possible. Here, our strategy is first to
define larger high-level functions that e.g. take a number of distances
and loads the table interactions for this interaction. This way we can
move this architecture-specific implementation to the SIMD module, and
both achieve a reasonably clean kernel but still optimize a lot.
</dd>
<dt>Architecture-specific kernels (directories/files)</dt>
<dd>
When it is absolutely impossible to use a shared implementation we might
have to code SIMD (just as GPU code). When this happens, we should create
subdirectory or otherwise clearly names a file with a suffix for the
SIMD architecture, to clarify to the user that the SIMD file has a
direct non-SIMD correspondence. Since this code can be very hard to read,
it is important to be explicit and use lots of comments - this is not the
type of code where you should use smart optimization with hundreds of
preprocessor directives. Keep it simple so other developers can help you
support it. The question is not whether you can get a function 20% faster,
but whether it justifies the added complexity of the code.
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
<dt>`gromacs/simd/impl_reference.h`</dt>
<dd>
This is an example of a low-level implementation. You should never, ever,
work directly with these in higher-level code. The reference implementation
contains the documentation for all SIMD wrappers, though.
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
<dt>`#gmx_simd_real_t`</dt>
<dd>
This is the SIMD-version of \Gromacs' real type,
which is set based on the CMake configuration and internally aliased
to one of the next two types.
Operations on these variables have the suffix `_r`, e.g. `gmx_simd_add_r()`.
</dd>
<dt>`#gmx_simd_float_t`</dt>
<dd>
This is always single-precision data, but it
might not be supported on all architectures. Suffix `_f` is used for
explicit single-precision routines, e.g. `gmx_simd_mul_f()`.
</dd>
<dt>`gmx_simd_double_t`</dt>
<dd>
This is always double precision when available,
and in rare cases you might want to use a specific precision.
Suffix `_d` is used for explicit double-precision routines,
e.g. `gmx_simd_mul_d()`
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
<dt>`#gmx_simd_int32_t`</dt>
<dd>
This is used for integers when converting to/from \Gromacs default "real" type.
The corresponding routines have suffix `_i`, e.g. `gmx_simd_add_i()`.
</dd>
<dt>`gmx_simd_fint32_t`</dt>
<dd>
Integers obtained when converting from single precision, or intended to be
converted to single precision floating-point. These are normal integers
(not a special conversion type), but since some SIMD architectures such as
SSE or AVX use different registers for integer SIMD variables having the
same width as float and double, respectively, we need to separate these
two types of integers. The actual operations you perform on the are normal
ones such as addition or multiplication. The routines
operating on these variables have suffix `_fi`, like `gmx_simd_add_fi()`.
This will also be the widest integer data type if you want to do pure
integer SIMD operations, but that will not be supported on all platforms.
</dd>
<dt>`gmx_simd_dint32_t`</dt>
<dd>
Integers used when converting to/from double. See the preceding item
for a detailed explanation. On many architectures,
including all x86 ones, this will be a narrower type than `gmx_simd_fint32_t`.
The correspoding routines have suffix `_di`, like `gmx_simd_add_di()`.
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
<dt>`#gmx_simd_bool_t`</dt>
<dd>
Results from boolean operations involving reals, and the booleans we use
to select between real values. The corresponding routines have suffix `_b`,
like `gmx_simd_or_b()`.
</dd>
<dt>`gmx_simd_fbool_t`</dt>
<dd>
Booleans specifically for single precision. Corresponding function suffix
is `_fb`, like `gmx_simd_or_fb()`.
</dd>
<dt>`gmx_simd_dbool_t`</dt>
<dd>
Operations specifically on double. Operations have suffix `_db`: `gmx_simd_or_db()`
</dd>
<dt>`#gmx_simd_ibool_t`</dt>
<dd>
Boolean operations on integers corresponding to real (see floating-point
descriptions above). Operations on these booleans use suffix `_ib`,
like `gmx_simd_or_ib()`.
</dd>
<dt>`gmx_simd_fibool_t`</dt>
<dd>
Booleans for integers corresponding to float. Operation suffix is `_fib`,
like `gmx_simd_or_fib()`.
</dd>
<dt>`gmx_simd_dibool_t`</dt>
<dd>
Booleans for integers corresponding to double. Operation suffix is `_dib`,
like `gmx_simd_or_dib()`.
</dd>
</dl>

The subset you should use in practice
-------------------------------------

If this seems daunting, in practice you should only need to use these types
when you start coding:

<dl>
<dt>`#gmx_simd_real_t`</dt>
<dd>
Floating-point data.
</dd>
<dt>`#gmx_simd_bool_t`</dt>
<dd>
Booleans.
</dd>
<dt>`#gmx_simd_int32_t`</dt>
<dd>
Integer data. Might not be supported, so you must check
the preprocessor macros described below.
</dd>
</dl>

Operations on these types will be defined to either float/double (or corresponding integers) based on the current \Gromacs precision, so the documentation is occasionally more detailed for the lower-level actual implementation functions.

SIMD4 Macros
------------

The above should be sufficient for code that works with the full SIMD width.
Unfortunately reality is not that simple. Some algorithms like lattice
summation need quartets of elements, so even when the SIMD width is >4 we
need width-4 SIMD if it is supported. These datatypes and operations use the
prefix `gmx_simd4_`, and availability is indicated by `GMX_SIMD4_HAVE_FLOAT`
and `GMX_SIMD4_HAVE_DOUBLE`. For now we only support a small subset of SIMD
operations for SIMD4, but that is trivial to extend if we need to.

Predefined SIMD preprocessor macros
===================================

Functionality-wise, we have a small set of core set of features that we
require to be present on all platforms, while more avanced features can be
used in the code when defines like e.g. `GMX_SIMD_HAVE_LOADU` are set.

This is a summary of the currently available preprocessor defines that
you should use to check for support when using the corresponding features.
We first list the float/double/int defines set by the _implementation_; in
most cases you do not want to check directly for float/double defines, but
you should instead use the derived "real" defines set in this file - we list
those at the end below.

Preprocessor predefined macro defines set by the low-level implementation.
These are only set if they work for all datatypes; `GMX_SIMD_HAVE_LOADU`
thus means we can load both float, double, and integers from unaligned memory,
and that the unaligned loads are available for SIMD4 too.

<dl>
<dt>`GMX_SIMD_HAVE_FLOAT`</dt>
<dd>
Single-precision instructions available.
</dd>
<dt>`GMX_SIMD_HAVE_DOUBLE `</dt>
<dd>
Double-precision instructions available.
</dd>
<dt>`GMX_SIMD_HAVE_HARDWARE`</dt>
<dd>
Set when we are NOT emulating SIMD.
</dd>
<dt>`GMX_SIMD_HAVE_LOADU`</dt>
<dd>
Load from unaligned memory available.
</dd>
<dt>`GMX_SIMD_HAVE_STOREU`</dt>
<dd>
Store to unaligned memory available.
</dd>
<dt>`GMX_SIMD_HAVE_LOGICAL`</dt>
<dd>
Support for and/andnot/or/xor on floating-point variables.
</dd>
<dt>`GMX_SIMD_HAVE_FMA`</dt>
<dd>
Floating-point fused multiply-add.
Note: We provide emulated FMA instructions if you do not have FMA
support, but in that case you might be able to code it more efficient w/o FMA.
</dd>
<dt>`GMX_SIMD_HAVE_FRACTION`</dt>
<dd>
Instruction to get decimal fraction. Same as FMA: This denotes
hardware support, otherwise instruction will be emulated.
</dd>
<dt>`GMX_SIMD_HAVE_FINT32`</dt>
<dd>
Integer conversions to/from float available.
</dd>
<dt>`GMX_SIMD_HAVE_FINT32_EXTRACT`</dt>
<dd>
Support for extracting integer SIMD elements from `gmx_simd_fint32_t`.
</dd>
<dt>`GMX_SIMD_HAVE_FINT32_LOGICAL`</dt>
<dd>
Bitwise shifts on `gmx_simd_fint32_t`.
</dd>
<dt>`GMX_SIMD_HAVE_FINT32_ARITHMETICS`</dt>
<dd>
Arithmetic ops for `gmx_simd_fint32_t`.
</dd>
<dt>`GMX_SIMD_HAVE_DINT32`</dt>
<dd>
Integer conversions to/from double available.
</dd>
<dt>`GMX_SIMD_HAVE_DINT32_EXTRACT`</dt>
<dd>
Support for extracting integer SIMD elements from `gmx_simd_dint32_t`.
</dd>
<dt>`GMX_SIMD_HAVE_DINT32_LOGICAL`</dt>
<dd>
Bitwise shifts on `gmx_simd_dint32_t`.
</dd>
<dt>`GMX_SIMD_HAVE_DINT32_ARITHMETICS`</dt>
<dd>
Arithmetic ops for `gmx_simd_dint32_t`.
</dd>
</dl>

There are also two macros specific to SIMD4: `GMX_SIMD4_HAVE_FLOAT` is set
if we can use SIMD4 in single precision, and `GMX_SIMD4_HAVE_DOUBLE`
similarly denotes support for a double-precision SIMD4 implementation. For
generic properties (e.g. whether SIMD4 FMA is supported), you should check
the normal SIMD macros above.

Implementation properties
-------------------------

Higher-level code can use these macros to find information about the implementation,
for instance what the SIMD width is:

<dl>
<dt>`GMX_SIMD_FLOAT_WIDTH`</dt>
<dd>
Number of elements in `gmx_simd_float_t`, and practical width of `gmx_simd_fint32_t`.
</dd>
<dt>`GMX_SIMD_DOUBLE_WIDTH`</dt>
<dd>
Number of elements in `gmx_simd_double_t`, and practical width of `gmx_simd_dint32_t`</dd>
<dt>`GMX_SIMD_RSQRT_BITS`</dt>
<dd>
Accuracy (bits) of 1/sqrt(x) lookup step.
</dd>
<dt>`GMX_SIMD_RCP_BITS`</dt>
<dd>
Accuracy (bits) of 1/x lookup step.
</dd>
</dl>

After including the low-level architecture-specific implementation, this
header sets the following derived defines based on the current precision;
these are the ones you should check for unless you absolutely want to dig
deep into the explicit single/double precision implementations:

<dl>
<dt>`GMX_SIMD_HAVE_REAL`</dt>
<dd>
Set either to `GMX_SIMD_HAVE_FLOAT` or `GMX_SIMD_HAVE_DOUBLE`
</dd>
<dt>`GMX_SIMD4_HAVE_REAL`</dt>
<dd>
Set either to `GMX_SIMD4_HAVE_FLOAT` or `GMX_SIMD4_HAVE_DOUBLE`
</dd>
<dt>`GMX_SIMD_REAL_WIDTH`</dt>
<dd>
Set either to `GMX_SIMD_FLOAT_WIDTH` or `GMX_SIMD_DOUBLE_WIDTH`
</dd>
<dt>`GMX_SIMD_HAVE_INT32`</dt>
<dd>
Set either to `GMX_SIMD_HAVE_FINT32` or `GMX_SIMD_HAVE_DINT32`
</dd>
<dt>`GMX_SIMD_INT32_WIDTH`</dt>
<dd>
Set either to `GMX_SIMD_FINT32_WIDTH` or `GMX_SIMD_DINT32_WIDTH`
</dd>
<dt>`GMX_SIMD_HAVE_INT32_EXTRACT`</dt>
<dd>
Set either to `GMX_SIMD_HAVE_FINT32_EXTRACT` or `GMX_SIMD_HAVE_DINT32_EXTRACT`
</dd>
<dt>`GMX_SIMD_HAVE_INT32_LOGICAL`</dt>
<dd>
Set either to `GMX_SIMD_HAVE_FINT32_LOGICAL` or `GMX_SIMD_HAVE_DINT32_LOGICAL`
</dd>
<dt>`GMX_SIMD_HAVE_INT32_ARITHMETICS`</dt>
<dd>
Set either to `GMX_SIMD_HAVE_FINT32_ARITHMETICS` or `GMX_SIMD_HAVE_DINT32_ARITHMETICS`
</dd>
</dl>

For convenience we also define `GMX_SIMD4_WIDTH` to 4. This will never vary,
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
(such as `gmx_simd_sincos_d_sse2`) so we can override it (e.g. in SSE4) by
simply undefining and setting a new definition. Still, this is an
implementation detail you won't have to worry about until you start writing
support for a new SIMD architecture.


Automated checking
------------------

Having fallback implementations when SIMD is not supported can be a
performance problem if the code does not correctly include
`gromacs/simd/simd.h`, particularly after refactoring.
`make check-source` checks the whole code for the use of symbols defined
in `gromacs/simd/simd.h` and requires that files using those symbols
do the correct include. Similar checking is done for higher-level
SIMD-management headers, e.g. `gromacs/ewald/pme-simd.h`.
