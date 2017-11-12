/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */

/*! \libinternal
 * \defgroup module_simd SIMD intrinsics interface (simd)
 * \ingroup group_utilitymodules
 *
 * \brief Provides an architecture-independent way of doing SIMD coding.
 *
 * Overview of the SIMD implementation is provided in \ref page_simd.
 * The details are documented in gromacs/simd/simd.h and the reference
 * implementation impl_reference.h.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 */

#ifndef GMX_SIMD_SIMD_H
#define GMX_SIMD_SIMD_H

/*! \libinternal \file
 *
 * \brief Definitions, capabilities, and wrappers for SIMD module.
 *
 * The macros in this file are intended to be used for writing
 * architecture-independent SIMD intrinsics code.
 * To support a new architecture, adding a new sub-include with macros here
 * should be (nearly) all that is needed.
 *
 * The defines in this top-level file will set default Gromacs real precision
 * operations to either single or double precision based on whether
 * GMX_DOUBLE is 1. The actual implementation - including e.g.
 * conversion operations specifically between single and double - is documented
 * in impl_reference.h.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \inlibraryapi
 * \ingroup module_simd
 */

#include "config.h"

#include <cstddef>
#include <cstdint>

#include <array>

#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

//! \cond libapi


/*! \addtogroup module_simd
 * \{
 */

namespace gmx
{
/*! \libinternal \brief Tag type to select to load SimdFloat with simdLoad(U) */
struct SimdFloatTag {};
/*! \libinternal \brief Tag type to select to load SimdDouble with simdLoad(U) */
struct SimdDoubleTag {};
/*! \libinternal \brief Tag type to select to load SimdFInt32 with simdLoad(U) */
struct SimdFInt32Tag {};
/*! \libinternal \brief Tag type to select to load SimdDInt32 with simdLoad(U) */
struct SimdDInt32Tag {};
}

/*! \name SIMD predefined macros to describe high-level capabilities
 *
 *  These macros are used to describe the features available in default
 *  Gromacs real precision. They are set from the lower-level implementation
 *  files that have macros describing single and double precision individually,
 *  as well as the implementation details.
 *  \{
 */

#if GMX_SIMD_X86_SSE2
#    include "impl_x86_sse2/impl_x86_sse2.h"
#elif GMX_SIMD_X86_SSE4_1
#    include "impl_x86_sse4_1/impl_x86_sse4_1.h"
#elif GMX_SIMD_X86_AVX_128_FMA
#    include "impl_x86_avx_128_fma/impl_x86_avx_128_fma.h"
#elif GMX_SIMD_X86_AVX_256
#    include "impl_x86_avx_256/impl_x86_avx_256.h"
#elif GMX_SIMD_X86_AVX2_256
#    include "impl_x86_avx2_256/impl_x86_avx2_256.h"
#elif GMX_SIMD_X86_AVX2_128
#    include "impl_x86_avx2_128/impl_x86_avx2_128.h"
#elif GMX_SIMD_X86_MIC
#    include "impl_x86_mic/impl_x86_mic.h"
#elif GMX_SIMD_X86_AVX_512
#    include "impl_x86_avx_512/impl_x86_avx_512.h"
#elif GMX_SIMD_X86_AVX_512_KNL
#    include "impl_x86_avx_512_knl/impl_x86_avx_512_knl.h"
#elif GMX_SIMD_ARM_NEON
#    include "impl_arm_neon/impl_arm_neon.h"
#elif GMX_SIMD_ARM_NEON_ASIMD
#    include "impl_arm_neon_asimd/impl_arm_neon_asimd.h"
#elif GMX_SIMD_IBM_QPX
#    include "impl_ibm_qpx/impl_ibm_qpx.h"
#elif GMX_SIMD_IBM_VMX
#    include "impl_ibm_vmx/impl_ibm_vmx.h"
#elif GMX_SIMD_IBM_VSX
#    include "impl_ibm_vsx/impl_ibm_vsx.h"
#elif (GMX_SIMD_REFERENCE || defined DOXYGEN)
#    include "impl_reference/impl_reference.h" // Includes doxygen documentation
#else
#    include "impl_none/impl_none.h"
#endif

// The scalar SIMD-mimicking functions are always included so we can use
// templated functions even without SIMD support.
#include "gromacs/simd/scalar/scalar.h"
#include "gromacs/simd/scalar/scalar_math.h"
#include "gromacs/simd/scalar/scalar_util.h"


#if GMX_DOUBLE
#    define GMX_SIMD_HAVE_REAL                                     GMX_SIMD_HAVE_DOUBLE
#    define GMX_SIMD_REAL_WIDTH                                    GMX_SIMD_DOUBLE_WIDTH
#    define GMX_SIMD_HAVE_INT32_EXTRACT                            GMX_SIMD_HAVE_DINT32_EXTRACT
#    define GMX_SIMD_HAVE_INT32_LOGICAL                            GMX_SIMD_HAVE_DINT32_LOGICAL
#    define GMX_SIMD_HAVE_INT32_ARITHMETICS                        GMX_SIMD_HAVE_DINT32_ARITHMETICS
#    define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL    GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE
#    define GMX_SIMD_HAVE_HSIMD_UTIL_REAL                          GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE
#    define GMX_SIMD4_HAVE_REAL                                    GMX_SIMD4_HAVE_DOUBLE
#else // GMX_DOUBLE

/*! \brief 1 if SimdReal is available, otherwise 0.
 *
 *  \ref GMX_SIMD_HAVE_DOUBLE if GMX_DOUBLE is 1, otherwise \ref GMX_SIMD_HAVE_FLOAT.
 */
#    define GMX_SIMD_HAVE_REAL               GMX_SIMD_HAVE_FLOAT

/*! \brief Width of SimdReal.
 *
 *  \ref GMX_SIMD_DOUBLE_WIDTH if GMX_DOUBLE is 1, otherwise \ref GMX_SIMD_FLOAT_WIDTH.
 */
#    define GMX_SIMD_REAL_WIDTH              GMX_SIMD_FLOAT_WIDTH

/*! \brief 1 if support is available for extracting elements from SimdInt32, otherwise 0
 *
 *  \ref GMX_SIMD_HAVE_DINT32_EXTRACT if GMX_DOUBLE is 1, otherwise
 *  \ref GMX_SIMD_HAVE_FINT32_EXTRACT.
 */
#    define GMX_SIMD_HAVE_INT32_EXTRACT      GMX_SIMD_HAVE_FINT32_EXTRACT

/*! \brief 1 if logical ops are supported on SimdInt32, otherwise 0.
 *
 *  \ref GMX_SIMD_HAVE_DINT32_LOGICAL if GMX_DOUBLE is 1, otherwise
 *  \ref GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
#    define GMX_SIMD_HAVE_INT32_LOGICAL      GMX_SIMD_HAVE_FINT32_LOGICAL

/*! \brief 1 if arithmetic ops are supported on SimdInt32, otherwise 0.
 *
 *  \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS if GMX_DOUBLE is 1, otherwise
 *  \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS.
 */
#    define GMX_SIMD_HAVE_INT32_ARITHMETICS  GMX_SIMD_HAVE_FINT32_ARITHMETICS

/*! \brief 1 if gmx::simdGatherLoadUBySimdIntTranspose is present, otherwise 0
 *
 *  \ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE if GMX_DOUBLE is 1, otherwise
 *  \ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT.
 */
#    define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL    GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT

/*! \brief 1 if real half-register load/store/reduce utils present, otherwise 0
 *
 *  \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE if GMX_DOUBLE is 1, otherwise
 *  \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT.
 */
#    define GMX_SIMD_HAVE_HSIMD_UTIL_REAL    GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT

/*! \brief 1 if Simd4Real is available, otherwise 0.
 *
 *  \ref GMX_SIMD4_HAVE_DOUBLE if GMX_DOUBLE is 1, otherwise \ref GMX_SIMD4_HAVE_FLOAT.
 */
#    define GMX_SIMD4_HAVE_REAL              GMX_SIMD4_HAVE_FLOAT

#endif // GMX_DOUBLE

//! \}  end of name-group describing high-level capabilities

namespace gmx
{

template<class T, size_t N>
struct AlignedArray;

#if GMX_SIMD_HAVE_FLOAT
/*! \libinternal \brief Identical to std::array with GMX_SIMD_FLOAT_WIDTH alignment.
 *  Should not be deleted through base pointer (destructor is non-virtual).
 */
template<size_t N>
struct alignas(GMX_SIMD_FLOAT_WIDTH*sizeof(float))AlignedArray<float, N> : public std::array<float, N>
{
};
#endif

#if GMX_SIMD_HAVE_DOUBLE
/*! \libinternal \brief  Identical to std::array with GMX_SIMD_DOUBLE_WIDTH alignment.
 *  Should not be deleted through base pointer (destructor is non-virtual).
 */
template<size_t N>
struct alignas(GMX_SIMD_DOUBLE_WIDTH*sizeof(double))AlignedArray<double, N> : public std::array<double, N>
{
};
#endif

#if GMX_SIMD_HAVE_REAL

/*! \name SIMD data types
 *
 *  The actual storage of these types is implementation dependent. The
 *  documentation is generated from the reference implementation, but for
 *  normal usage this will likely not be what you are using.
 * \{
 */

/*! \brief Real precision floating-point SIMD datatype.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_REAL is 1.
 *
 * \ref SimdDouble if GMX_DOUBLE is 1, otherwise \ref SimdFloat.
 *
 * \note This variable cannot be placed inside other structures or classes, since
 *       some compilers (including at least clang-3.7) appear to lose the
 *       alignment. This is likely particularly severe when allocating such
 *       memory on the heap, but it occurs for stack structures too.
 */
#    if GMX_DOUBLE
typedef SimdDouble               SimdReal;
#    else
typedef SimdFloat                SimdReal;
#    endif


/*! \brief Boolean SIMD type for usage with \ref SimdReal.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_REAL is 1.
 *
 * If GMX_DOUBLE is 1, this will be set to \ref SimdDBool
 * internally, otherwise \ref SimdFBool. This is necessary since some
 * SIMD implementations use bitpatterns for marking truth, so single-
 * vs. double precision booleans are not necessarily exchangable.
 * As long as you just use this type you will not have to worry about precision.
 *
 * See \ref SimdIBool for an explanation of real vs. integer booleans.
 *
 * \note This variable cannot be placed inside other structures or classes, since
 *       some compilers (including at least clang-3.7) appear to lose the
 *       alignment. This is likely particularly severe when allocating such
 *       memory on the heap, but it occurs for stack structures too.
 */
#    if GMX_DOUBLE
typedef SimdDBool                SimdBool;
#    else
typedef SimdFBool                SimdBool;
#    endif


/*! \brief 32-bit integer SIMD type.
 *
 * If GMX_DOUBLE is 1, this will be set to \ref SimdDInt32
 * internally, otherwise \ref SimdFInt32. This might seem a strange
 * implementation detail, but it is because some SIMD implementations use
 * different types/widths of integers registers when converting from
 * double vs. single precision floating point. As long as you just use
 * this type you will not have to worry about precision.
 *
 * \note This variable cannot be placed inside other structures or classes, since
 *       some compilers (including at least clang-3.7) appear to lose the
 *       alignment. This is likely particularly severe when allocating such
 *       memory on the heap, but it occurs for stack structures too.
 */
#    if GMX_DOUBLE
typedef SimdDInt32               SimdInt32;
#    else
typedef SimdFInt32               SimdInt32;
#    endif

#if GMX_SIMD_HAVE_INT32_ARITHMETICS
/*! \brief Boolean SIMD type for usage with \ref SimdInt32.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_INT32_ARITHMETICS is 1.
 *
 * If GMX_DOUBLE is 1, this will be set to \ref SimdDIBool
 * internally, otherwise \ref SimdFIBool. This is necessary since some
 * SIMD implementations use bitpatterns for marking truth, so single-
 * vs. double precision booleans are not necessarily exchangable, and while
 * a double-precision boolean might be represented with a 64-bit mask, the
 * corresponding integer might only use a 32-bit mask.
 *
 * We provide conversion routines for these cases, so the only thing you need to
 * keep in mind is to use \ref SimdBool when working with
 * \ref SimdReal while you pick \ref SimdIBool when working with
 * \ref SimdInt32 .
 *
 * To convert between them, use \ref cvtB2IB and \ref cvtIB2B.
 *
 * \note This variable cannot be placed inside other structures or classes, since
 *       some compilers (including at least clang-3.7) appear to lose the
 *       alignment. This is likely particularly severe when allocating such
 *       memory on the heap, but it occurs for stack structures too.
 */
#    if GMX_DOUBLE
typedef SimdDIBool               SimdIBool;
#    else
typedef SimdFIBool               SimdIBool;
#    endif
#endif  // GMX_SIMD_HAVE_INT32_ARITHMETICS


#if GMX_DOUBLE
const int c_simdBestPairAlignment = c_simdBestPairAlignmentDouble;
#else
const int c_simdBestPairAlignment = c_simdBestPairAlignmentFloat;
#endif

#endif  // GMX_SIMD_HAVE_REAL

#if GMX_SIMD4_HAVE_REAL
/*! \brief Real precision floating-point SIMD4 datatype.
 *
 * This type is only available if \ref GMX_SIMD4_HAVE_REAL is 1.
 *
 * \ref Simd4Double if GMX_DOUBLE is 1, otherwise \ref Simd4Float.
 *
 * \note This variable cannot be placed inside other structures or classes, since
 *       some compilers (including at least clang-3.7) appear to lose the
 *       alignment. This is likely particularly severe when allocating such
 *       memory on the heap, but it occurs for stack structures too.
 */
#    if GMX_DOUBLE
typedef Simd4Double               Simd4Real;
#    else
typedef Simd4Float                Simd4Real;
#    endif


/*! \brief Boolean SIMD4 type for usage with \ref SimdReal.
 *
 * This type is only available if \ref GMX_SIMD4_HAVE_REAL is 1.
 *
 * If GMX_DOUBLE is 1, this will be set to \ref Simd4DBool
 * internally, otherwise \ref Simd4FBool. This is necessary since some
 * SIMD implementations use bitpatterns for marking truth, so single-
 * vs. double precision booleans are not necessarily exchangable.
 * As long as you just use this type you will not have to worry about precision.
 *
 * \note This variable cannot be placed inside other structures or classes, since
 *       some compilers (including at least clang-3.7) appear to lose the
 *       alignment. This is likely particularly severe when allocating such
 *       memory on the heap, but it occurs for stack structures too.
 */
#    if GMX_DOUBLE
typedef Simd4DBool                Simd4Bool;
#    else
typedef Simd4FBool                Simd4Bool;
#    endif
#endif // GMX_SIMD4_HAVE_REAL

//! \}  end of name-group describing SIMD data types

/*! \name High-level SIMD proxy objects to disambiguate load/set operations
 * \{
 */

/*! \libinternal \brief Simd traits */
template<typename T>
struct SimdTraits {};

#if GMX_SIMD_HAVE_FLOAT
template<>
struct SimdTraits<SimdFloat>
{
    using type = float;
    static constexpr int width = GMX_SIMD_FLOAT_WIDTH;
    using tag = SimdFloatTag;
};
#endif
#if GMX_SIMD_HAVE_DOUBLE
template<>
struct SimdTraits<SimdDouble>
{
    using type = double;
    static constexpr int width = GMX_SIMD_DOUBLE_WIDTH;
    using tag = SimdDoubleTag;
};
#endif
#if GMX_SIMD_HAVE_FLOAT
template<>
struct SimdTraits<SimdFInt32>
{
    using type = int;
    static constexpr int width = GMX_SIMD_FINT32_WIDTH;
    using tag = SimdFInt32Tag;
};
#endif
#if GMX_SIMD_HAVE_DOUBLE
template<>
struct SimdTraits<SimdDInt32>
{
    using type = int;
    static constexpr int width = GMX_SIMD_DINT32_WIDTH;
    using tag = SimdDInt32Tag;
};
#endif

template<typename T>
struct SimdTraits<const T>
{
    using type = const typename SimdTraits<T>::type;
    static constexpr int width = SimdTraits<T>::width;
    using tag = typename SimdTraits<T>::tag;
};

/*! \brief Load function that returns SIMD or scalar
 *
 * \tparam T Type to load (type is always mandatory)
 * \param  m Pointer to aligned memory
 * \return   Loaded value
 */
template<typename T>
static inline T
load(const typename SimdTraits<T>::type *m) //disabled by SFINAE for non-SIMD types
{
    return simdLoad(m, typename SimdTraits<T>::tag());
}

template<typename T>
static inline T
/* the enable_if serves to prevent two different type of misuse:
 * 1) load<SimdReal>(SimdReal*); should only be called on real* or int*
 * 2) load(real*); template parameter is mandatory because otherwise ambiguity is
 *    created. The dependent type disables type deduction.
 */
load(const typename std::enable_if<std::is_arithmetic<T>::value, T>::type *m)
{
    return *m;
}

template <typename T, size_t N>
static inline T gmx_simdcall
load(const AlignedArray<typename SimdTraits<T>::type, N> &m)
{
    return simdLoad(m.data(), typename SimdTraits<T>::tag());
}

/*! \brief Load function that returns SIMD or scalar based on template argument
 *
 * \tparam T Type to load (type is always mandatory)
 * \param m Pointer to unaligned memory
 * \return Loaded SimdFloat/Double/Int or basic scalar type
 */
template<typename T>
static inline T
loadU(const typename SimdTraits<T>::type *m)
{
    return simdLoadU(m, typename SimdTraits<T>::tag());
}

template<typename T>
static inline T
loadU(const typename std::enable_if<std::is_arithmetic<T>::value, T>::type *m)
{
    return *m;
}

template <typename T, size_t N>
static inline T gmx_simdcall
loadU(const AlignedArray<typename SimdTraits<T>::type, N> &m)
{
    return simdLoadU(m.data(), typename SimdTraits<T>::tag());
}

class SimdSetZeroProxyInternal;

static inline const SimdSetZeroProxyInternal gmx_simdcall
setZero();

/*! \libinternal \brief Proxy object to enable setZero() for SIMD and real types.
 *
 * This object is returned by setZero(), and depending on what type you assign
 * the result to the conversion method will call the right low-level function.
 */
class SimdSetZeroProxyInternal
{
    public:
        //!\brief Conversion method that returns 0.0 as float
        operator float() const { return 0.0f; }
        //!\brief Conversion method that returns 0.0 as double
        operator double() const { return 0.0; }
        //!\brief Conversion method that returns 0.0 as int32
        operator std::int32_t() const { return 0; }
#if GMX_SIMD_HAVE_FLOAT
        //!\brief Conversion method that will execute setZero() for SimdFloat
        operator SimdFloat() const { return setZeroF(); }
        //!\brief Conversion method that will execute setZero() for SimdFInt32
        operator SimdFInt32() const { return setZeroFI(); }
#endif
#if GMX_SIMD4_HAVE_FLOAT
        //!\brief Conversion method that will execute setZero() for Simd4Float
        operator Simd4Float() const { return simd4SetZeroF(); }
#endif
#if GMX_SIMD_HAVE_DOUBLE
        //!\brief Conversion method that will execute setZero() for SimdDouble
        operator SimdDouble() const { return setZeroD(); }
        //!\brief Conversion method that will execute setZero() for SimdDInt32
        operator SimdDInt32() const { return setZeroDI(); }
#endif
#if GMX_SIMD4_HAVE_DOUBLE
        //!\brief Conversion method that will execute setZero() for Simd4Double
        operator Simd4Double() const { return simd4SetZeroD(); }
#endif

    private:
        //! \brief Private constructor can only be called from setZero()
        SimdSetZeroProxyInternal() {}

        friend const SimdSetZeroProxyInternal gmx_simdcall
        setZero();

        GMX_DISALLOW_COPY_AND_ASSIGN(SimdSetZeroProxyInternal);
};

/*! \brief Proxy object to set any SIMD or scalar variable to zero
 *
 * \return Proxy object that will call the actual function to set a SIMD/scalar
 *         variable to zero based on the conversion function called when you
 *         assign the result.
 */
static inline const SimdSetZeroProxyInternal gmx_simdcall
setZero()
{
    return {};
}
//! \}  end of name-group proxy objects

}      // namespace gmx

// \}          end of module_simd

//! \endcond   end of condition libapi


#if GMX_SIMD_HAVE_FLOAT

/*! \brief Returns whether a pointer to float is aligned to a SIMD boundary
 *
 * \param[in] ptr  A pointer to a float
 */
static inline bool isSimdAligned(const float *ptr)
{
    return reinterpret_cast<std::size_t>(ptr) % (GMX_SIMD_FLOAT_WIDTH*sizeof(float)) == 0;
}

#endif // GMX_SIMD_HAVE_FLOAT

#if GMX_SIMD_HAVE_DOUBLE

/*! \brief Returns whether a pointer to double is aligned to a SIMD boundary
 *
 * \param[in] ptr  A pointer to a double
 */
static inline bool isSimdAligned(const double *ptr)
{
    return reinterpret_cast<std::size_t>(ptr) % (GMX_SIMD_DOUBLE_WIDTH*sizeof(double)) == 0;
}

#endif // GMX_SIMD_HAVE_DOUBLE


#if GMX_SIMD_HAVE_REAL
#if GMX_SIMD_REAL_WIDTH > GMX_REAL_MAX_SIMD_WIDTH
#error "GMX_SIMD_REAL_WIDTH > GMX_REAL_MAX_SIMD_WIDTH: increase GMX_REAL_MAX_SIMD_WIDTH in real.h"
#endif
#endif


#if 0
/* This is a hack to cover the corner case of using an
   explicit GMX_SIMD_HAVE_FLOAT or GMX_SIMD_HAVE_DOUBLE, rather than
   GMX_SIMD_HAVE_REAL.

   Such code is expected to include simd.h to get those symbols
   defined, but the actual definitions are in the implemention headers
   included by simd.h. check-source.py is not a full preprocessor, so
   it does not see the definitions in the implementation headers as
   belonging to simd.h, thus it cannot check that simd.h is being used
   correctly in the above hypothetical corner case. However, the
   checker also does not parse #if 0, so we can fool the checker into
   thinking that definition occurs here, and that will work well
   enough.

   If there's ever other kinds of SIMD code that might have the same
   problem, we might want to add other variables here.
 */
#    define GMX_SIMD_HAVE_FLOAT         1
#    define GMX_SIMD_HAVE_DOUBLE        1

#endif // end of hack

// The ArrayRef<SimdReal> specialization is always included, because compiler
// errors are confusing when template specialization aren't available.
#include "gromacs/simd/simd_memory.h"

#endif // GMX_SIMD_SIMD_H
