/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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
 * The details are documented in simd.h and the reference implementation
 * impl_reference.h.
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
 * GMX_DOUBLE is defined. The actual implementation - including e.g.
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

#include "gromacs/utility/real.h"

/*! \cond libapi */


/*! \addtogroup module_simd
 * \{
 */


/*! \name SIMD predefined macros to describe high-level capabilities
 *
 *  These macros are used to describe the features available in default
 *  Gromacs real precision. They are set from the lower-level implementation
 *  files that have macros describing single and double precision individually,
 *  as well as the implementation details.
 *  \{
 */

/* Intel MIC is a bit special since it is a co-processor. This means the rest
 * of GROMACS (which runs on the CPU) can use a default SIMD set like AVX.
 * All functions in this SIMD module are static, so it will work perfectly fine
 * to include this file with different SIMD definitions for different files.
 */
#if GMX_SIMD_X86_AVX_512ER
#    include "impl_x86_avx_512er/impl_x86_avx_512er.h"
#elif GMX_SIMD_X86_AVX_512F
#    include "impl_x86_avx_512f/impl_x86_avx_512f.h"
#elif GMX_SIMD_X86_MIC
#    include "impl_intel_mic/impl_intel_mic.h"
#elif GMX_SIMD_X86_AVX2_256
#    include "impl_x86_avx2_256/impl_x86_avx2_256.h"
#elif GMX_SIMD_X86_AVX_256
#    include "impl_x86_avx_256/impl_x86_avx_256.h"
#elif GMX_SIMD_X86_AVX_128_FMA
#    include "impl_x86_avx_128_fma/impl_x86_avx_128_fma.h"
#elif GMX_SIMD_X86_SSE4_1
#    include "impl_x86_sse4_1/impl_x86_sse4_1.h"
#elif GMX_SIMD_X86_SSE2
#    include "impl_x86_sse2/impl_x86_sse2.h"
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
#elif GMX_SIMD_SPARC64_HPC_ACE
#    include "impl_sparc64_hpc_ace/impl_sparc64_hpc_ace.h"
#elif (GMX_SIMD_REFERENCE || defined DOXYGEN)
/* Plain C SIMD reference implementation, also serves as documentation. */
#    include "impl_reference/impl_reference.h"
#else
#    include "impl_none/impl_none.h"
#endif

/* These convenience macros are ugly hacks where some source files still make
 * assumptions about the SIMD architecture. They will be removed as we implement
 * the new verlet kernels, but for now we need them, and to make sure they
 * always have values 0 or 1 we define them here rather than in the implementations.
 */
#define GMX_SIMD_X86_AVX2_256_OR_HIGHER      (GMX_SIMD_X86_AVX2_256)
#define GMX_SIMD_X86_AVX_256_OR_HIGHER       (GMX_SIMD_X86_AVX2_256_OR_HIGHER || GMX_SIMD_X86_AVX_256)
#define GMX_SIMD_X86_AVX_128_FMA_OR_HIGHER   (GMX_SIMD_X86_AVX_128_FMA)
#define GMX_SIMD_X86_SSE4_1_OR_HIGHER        (GMX_SIMD_X86_AVX_256_OR_HIGHER || GMX_SIMD_X86_AVX_128_FMA_OR_HIGHER || GMX_SIMD_X86_SSE4_1)
#define GMX_SIMD_X86_SSE2_OR_HIGHER          (GMX_SIMD_X86_SSE4_1_OR_HIGHER || GMX_SIMD_X86_SSE2)


#ifdef GMX_DOUBLE
#    define GMX_SIMD_HAVE_REAL               GMX_SIMD_HAVE_DOUBLE
#    define GMX_SIMD_REAL_WIDTH              GMX_SIMD_DOUBLE_WIDTH
#    define GMX_SIMD_HAVE_INT32              GMX_SIMD_HAVE_DINT32
#    define GMX_SIMD_INT32_WIDTH             GMX_SIMD_DINT32_WIDTH
#    define GMX_SIMD_HAVE_INT32_EXTRACT      GMX_SIMD_HAVE_DINT32_EXTRACT
#    define GMX_SIMD_HAVE_INT32_LOGICAL      GMX_SIMD_HAVE_DINT32_LOGICAL
#    define GMX_SIMD_HAVE_INT32_ARITHMETICS  GMX_SIMD_HAVE_DINT32_ARITHMETICS
#    define GMX_SIMD4_HAVE_REAL              GMX_SIMD4_HAVE_DOUBLE
#else // GMX_DOUBLE

/*! \brief 1 if SimdReal is available, otherwise 0.
 *
 *  \ref GMX_SIMD_HAVE_DOUBLE if GMX_DOUBLE is set, otherwise \ref GMX_SIMD_HAVE_FLOAT.
 */
#    define GMX_SIMD_HAVE_REAL               GMX_SIMD_HAVE_FLOAT

/*! \brief Width of SimdReal.
 *
 *  \ref GMX_SIMD_DOUBLE_WIDTH if GMX_DOUBLE is set, otherwise \ref GMX_SIMD_FLOAT_WIDTH.
 */
#    define GMX_SIMD_REAL_WIDTH              GMX_SIMD_FLOAT_WIDTH

/*! \brief 1 if SimdInt32 is available, otherwise 0.
 *
 *  \ref GMX_SIMD_HAVE_DINT32 if GMX_DOUBLE is set, otherwise \ref GMX_SIMD_HAVE_FINT32.
 */
#    define GMX_SIMD_HAVE_INT32              GMX_SIMD_HAVE_FINT32

/*! \brief Width of SimdInt32.
 *
 *  \ref GMX_SIMD_DINT32_WIDTH if GMX_DOUBLE is set, otherwise \ref GMX_SIMD_FINT32_WIDTH.
 */
#    define GMX_SIMD_INT32_WIDTH             GMX_SIMD_FINT32_WIDTH

/*! \brief 1 if simdExtractI() is available, otherwise 0.
 *
 *  \ref GMX_SIMD_HAVE_DINT32_EXTRACT if GMX_DOUBLE is set, otherwise
 *  \ref GMX_SIMD_HAVE_FINT32_EXTRACT.
 */
#    define GMX_SIMD_HAVE_INT32_EXTRACT      GMX_SIMD_HAVE_FINT32_EXTRACT

/*! \brief 1 if logical ops are supported on SimdInt32, otherwise 0.
 *
 *  \ref GMX_SIMD_HAVE_DINT32_LOGICAL if GMX_DOUBLE is set, otherwise
 *  \ref GMX_SIMD_HAVE_FINT32_LOGICAL.
 */
#    define GMX_SIMD_HAVE_INT32_LOGICAL      GMX_SIMD_HAVE_FINT32_LOGICAL

/*! \brief 1 if arithmetic ops are supported on SimdInt32, otherwise 0.
 *
 *  \ref GMX_SIMD_HAVE_DINT32_ARITHMETICS if GMX_DOUBLE is set, otherwise
 *  \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS.
 */
#    define GMX_SIMD_HAVE_INT32_ARITHMETICS  GMX_SIMD_HAVE_FINT32_ARITHMETICS

/*! \brief 1 if Simd4Real is available, otherwise 0.
 *
 *  \ref GMX_SIMD4_HAVE_DOUBLE if GMX_DOUBLE is set, otherwise \ref GMX_SIMD4_HAVE_FLOAT.
 */
#    define GMX_SIMD4_HAVE_REAL              GMX_SIMD4_HAVE_FLOAT

#endif // GMX_DOUBLE

/*! \}  end of name-group describing high-level capabilities */

namespace gmx
{

/*! \name SIMD data types
 *
 *  The actual storage of these types is implementation dependent. The
 *  documentation is generated from the reference implementation, but for
 *  normal usage this will likely not be what you are using.
 * \{
 */

#if GMX_SIMD_HAVE_REAL
/*! \brief Real precision floating-point SIMD datatype.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_REAL is 1.
 *
 * \ref SimdDouble if GMX_DOUBLE is set, otherwise \ref SimdFloat.
 */
#    ifdef GMX_DOUBLE
typedef SimdDouble               SimdReal;
#    else
typedef SimdFloat                SimdReal;
#    endif


/*! \brief Boolean SIMD type for usage with \ref SimdReal.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_REAL is 1.
 *
 * If GMX_DOUBLE is defined, this will be set to \ref SimdDBool
 * internally, otherwise \ref SimdFBool. This is necessary since some
 * SIMD implementations use bitpatterns for marking truth, so single-
 * vs. double precision booleans are not necessarily exchangable.
 * As long as you just use this type you will not have to worry about precision.
 *
 * See \ref SimdIBool for an explanation of real vs. integer booleans.
 */
#    ifdef GMX_DOUBLE
typedef SimdDBool                SimdBool;
#    else
typedef SimdFBool                SimdBool;
#    endif
#endif // GMX_SIMD_HAVE_REAL


#if GMX_SIMD_HAVE_INT32
/*! \brief 32-bit integer SIMD type.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_INT32 is 1.
 *
 * If GMX_DOUBLE is defined, this will be set to \ref SimdDInt32
 * internally, otherwise \ref SimdFInt32. This might seem a strange
 * implementation detail, but it is because some SIMD implementations use
 * different types/widths of integers registers when converting from
 * double vs. single precision floating point. As long as you just use
 * this type you will not have to worry about precision.
 */
#    ifdef GMX_DOUBLE
typedef SimdDInt32               SimdInt32;
#    else
typedef SimdFInt32               SimdInt32;
#    endif

/*! \brief Boolean SIMD type for usage with \ref gmx::SimdInt32.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_INT32 is 1.
 *
 * If GMX_DOUBLE is defined, this will be set to \ref SimdDIBool
 * internally, otherwise \ref SimdFIBool. This is necessary since some
 * SIMD implementations use bitpatterns for marking truth, so single-
 * vs. double precision booleans are not necessarily exchangable, and while
 * a double-precision boolean might be represented with a 64-bit mask, the
 * corresponding integer might only use a 32-bit mask.
 *
 * We provide conversion routines for these cases, so the only thing you need to
 * keep in mind is to use \ref SimdBool when working with
 * \ref SimdReal while you pick \ref SimdIBool when working with
 * \ref gmx::SimdInt32.
 *
 * To convert between them, use \ref simdCvtB2IB and \ref simdCvtIB2B.
 */
#    ifdef GMX_DOUBLE
typedef SimdDIBool               SimdIBool;
#    else
typedef SimdFIBool               SimdIBool;
#    endif
#endif // GMX_SIMD_HAVE_INT32

/*! \}  end of name-group describing SIMD data types */



#if GMX_SIMD_HAVE_REAL
/*  \name SIMD load/store operations on SimdReal
 *
 *  \note Unaligned load/stores are only available when
 *  \ref GMX_SIMD_HAVE_LOADU and \ref GMX_SIMD_HAVE_STOREU are set, respectively.
 *
 *  \{
 */

/*! \brief Load \ref GMX_SIMD_REAL_WIDTH values from aligned memory to \ref SimdReal
 *
 * Uses \ref simdLoadD if GMX_DOUBLE is set, otherwise \ref simdLoadF.
 *
 * \copydetails simdLoadF
 */
static inline SimdReal
simdLoad(const real *m)
{
#ifdef GMX_DOUBLE
    return simdLoadD(m);
#else
    return simdLoadF(m);
#endif
}

/*! \brief Set all elements in \ref SimdReal from single value in memory.
 *
 * Uses \ref simdLoad1D if GMX_DOUBLE is set, otherwise \ref simdLoad1F.
 *
 * \copydetails simdLoad1F
 */
static inline SimdReal
simdLoad1(const real *m)
{
#ifdef GMX_DOUBLE
    return simdLoad1D(m);
#else
    return simdLoad1F(m);
#endif
}

/*! \brief Set all elements in \ref SimdReal from a scalar.
 *
 * Uses \ref simdSet1D if GMX_DOUBLE is set, otherwise \ref simdSet1F.
 *
 * \copydetails simdSet1F
 */
static inline SimdReal
simdSet1(const real r)
{
#ifdef GMX_DOUBLE
    return simdSet1D(r);
#else
    return simdSet1F(r);
#endif
}

/*! \brief Store \ref GMX_SIMD_REAL_WIDTH values from \ref SimdReal to aligned memory.
 *
 * Uses \ref simdStoreD if GMX_DOUBLE is set, otherwise \ref simdStore F.
 *
 * \copydetails simdStoreF
 */
static inline void
simdStore(real *m, SimdReal a)
{
#ifdef GMX_DOUBLE
    simdStoreD(m, a);
#else
    simdStoreF(m, a);
#endif
}


#if GMX_SIMD_HAVE_LOADU
/*! \brief Load \ref GMX_SIMD_REAL_WIDTH values from unaligned memory to \ref SimdReal.
 *
 * Uses \ref simdLoadUD if GMX_DOUBLE is set, otherwise \ref simdLoadUF.
 *
 * \copydetails simdLoadUF
 */
static inline SimdReal
simdLoadU(const real *m)
{
#ifdef GMX_DOUBLE
    return simdLoadUD(m);
#else
    return simdLoadUF(m);
#endif
}
#endif // GMX_SIMD_HAVE_LOADU


#if GMX_SIMD_HAVE_STOREU
/*! \brief Store \ref GMX_SIMD_REAL_WIDTH values from \ref SimdReal to unaligned memory.
 *
 * Uses \ref simdStoreUD if GMX_DOUBLE is set, otherwise \ref simdStoreUF.
 *
 * \copydetails simdStoreUF
 */
static inline void
simdStoreU(real *m, SimdReal a)
{
#ifdef GMX_DOUBLE
    simdStoreUD(m, a);
#else
    simdStoreUF(m, a);
#endif
}
#endif // GMX_SIMD_HAVE_STOREU


/*! \brief Set all elements in \ref SimdReal to 0.0.
 *
 * Uses \ref simdSetZeroD if GMX_DOUBLE is set, otherwise \ref simdSetZeroF.
 *
 * \copydetails simdSetZeroF
 */
static inline SimdReal
simdSetZero()
{
#ifdef GMX_DOUBLE
    return simdSetZeroD();
#else
    return simdSetZeroF();
#endif
}

/*! \}  end of name-group describing load/store on SimdReal */

#endif // GMX_SIMD_HAVE_REAL


#if GMX_SIMD_HAVE_INT32
/*! \name SIMD load/store operations on SimdInt32
 *
 *  \note Unaligned load/stores are only available when
 *  \ref GMX_SIMD_HAVE_LOADU and \ref GMX_SIMD_HAVE_STOREU are set, respectively.
 *  \{
 */

/*! \brief Load \ref GMX_SIMD_INT32_WIDTH values from aligned memory to \ref SimdInt32 .
 *
 * Uses \ref simdLoadDI if GMX_DOUBLE is set, otherwise \ref simdLoadFI .
 *
 * \copydetails simdLoadFI
 */
static inline SimdInt32
simdLoadI(const std::int32_t *m)
{
#ifdef GMX_DOUBLE
    return simdLoadDI(m);
#else
    return simdLoadFI(m);
#endif
}

/*! \brief Set all elements in \ref SimdInt32 from a single integer.
 *
 * Uses \ref simdSet1DI if GMX_DOUBLE is set, otherwise \ref simdSet1FI .
 *
 * \copydetails simdSet1FI
 */
static inline SimdInt32
simdSet1I(std::int32_t b)
{
#ifdef GMX_DOUBLE
    return simdSet1DI(b);
#else
    return simdSet1FI(b);
#endif
}

/*! \brief Store \ref GMX_SIMD_REAL_WIDTH values from \ref SimdInt32 to aligned memory.
 *
 * Uses \ref simdStoreDI if GMX_DOUBLE is set, otherwise \ref simdStoreFI .
 *
 * \copydetails simdStoreFI
 */
static inline void
simdStoreI(std::int32_t *m, SimdInt32 a)
{
#ifdef GMX_DOUBLE
    simdStoreDI(m, a);
#else
    simdStoreFI(m, a);
#endif
}

#if GMX_SIMD_HAVE_LOADU
/*! \brief Load \ref GMX_SIMD_REAL_WIDTH values from unaligned memory to \ref SimdInt32.
 *
 * Uses \ref simdLoadUDI if GMX_DOUBLE is set, otherwise \ref simdLoadUFI .
 *
 * \copydetails simdLoadUFI
 */
static inline SimdInt32
simdLoadUI(const std::int32_t *m)
{
#ifdef GMX_DOUBLE
    return simdLoadUDI(m);
#else
    return simdLoadUFI(m);
#endif
}
#endif // GMX_SIMD_HAVE_LOADU

#if GMX_SIMD_HAVE_STOREU
/*! \brief Store \ref GMX_SIMD_REAL_WIDTH values from \ref SimdInt32 to unaligned memory.
 *
 * Uses \ref simdStoreUDI if GMX_DOUBLE is set, otherwise \ref simdStoreUFI .
 *
 * \copydetails simdStoreUFI
 */
static inline void
simdStoreUI(std::int32_t *m, SimdInt32 a)
{
#ifdef GMX_DOUBLE
    simdStoreUDI(m, a);
#else
    simdStoreUFI(m, a);
#endif
}
#endif // GMX_SIMD_HAVE_STOREU


/*! \brief Set all elements in \ref SimdInt32 to 0.
 *
 * Uses \ref simdSetZeroDI if GMX_DOUBLE is set, otherwise \ref simdSetZeroFI.
 *
 * \copydetails simdSetZeroFI
 */
static inline SimdInt32
simdSetZeroI()
{
#ifdef GMX_DOUBLE
    return simdSetZeroDI();
#else
    return simdSetZeroFI();
#endif
}

#if GMX_SIMD_HAVE_INT32_EXTRACT
/*! \brief Extract single integer from \ref gmx::SimdInt32 element.
 *
 * \tparam index Compile-time constant, position to extract (first position is 0)
 * \param  a     SIMD variable from which to extract value.
 *
 * Uses \ref simdExtractDI if GMX_DOUBLE is set, otherwise \ref simdExtractFI.
 */
template<int index>
static inline std::int32_t
simdExtractI(SimdInt32 a)
{
#ifdef GMX_DOUBLE
    return simdExtractDI(a, index);
#else
    return simdExtractFI(a, index);
#endif
}

#endif // GMX_SIMD_HAVE_INT32_EXTRACT

/*! \}  end of name-group describing load/store on SimdInt32 */

#endif // GMX_SIMD_HAVE_INT32

#if GMX_SIMD_HAVE_LOGICAL
/*! \name SIMD floating-point logical operations on SimdReal
 *
 *  These instructions are available if \ref GMX_SIMD_HAVE_LOGICAL is 1.
 *  \{
 */

/*! \brief Bitwise \a and on two \ref SimdReal.
 *
 * Uses \ref simdAndD if GMX_DOUBLE is set, otherwise \ref simdAndF.
 *
 * \copydetails simdAndF
 */
static inline SimdReal
simdAnd(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdAndD(a, b);
#else
    return simdAndF(a, b);
#endif
}

/*! \brief Bitwise \a and-not on two \ref SimdReal; 1st arg is complemented.
 *
 * Uses \ref simdAndNotD if GMX_DOUBLE is set, otherwise \ref simdAndNotF.
 *
 * \copydetails simdAndNotF
 */
static inline SimdReal
simdAndNot(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdAndNotD(a, b);
#else
    return simdAndNotF(a, b);
#endif
}

/*! \brief Bitwise \a or on two \ref SimdReal.
 *
 * Uses \ref simdOrD if GMX_DOUBLE is set, otherwise \ref simdOrF.
 *
 * \copydetails simdOrF
 */
static inline SimdReal
simdOr(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdOrD(a, b);
#else
    return simdOrF(a, b);
#endif
}

/*! \brief Bitwise \a exclusive-or on two \ref SimdReal.
 *
 * Uses \ref simdXorD if GMX_DOUBLE is set, otherwise \ref simdXorF.
 *
 * \copydetails simdXorF
 */
static inline SimdReal
simdXor(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdXorD(a, b);
#else
    return simdXorF(a, b);
#endif
}

/*! \}   end of name-group describing FP logical */
#endif // GMX_SIMD_HAVE_LOGICAL

#if GMX_SIMD_HAVE_REAL
/*! \name SIMD floating-point arithmetic operations on SimdReal
 *  \{
 */

/*! \brief SIMD a+b for two \ref SimdReal.
 *
 * Uses \ref simdAddD if GMX_DOUBLE is set, otherwise \ref simdAddF.
 *
 * \copydetails simdAddF
 */
static inline SimdReal
simdAdd(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdAddD(a, b);
#else
    return simdAddF(a, b);
#endif
}

/*! \brief SIMD a-b for two \ref SimdReal.
 *
 * Uses \ref simdSubD if GMX_DOUBLE is set, otherwise \ref simdSubF.
 *
 * \copydetails simdSubF
 */
static inline SimdReal
simdSub(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdSubD(a, b);
#else
    return simdSubF(a, b);
#endif
}

/*! \brief SIMD a*b for two \ref SimdReal.
 *
 * Uses \ref simdMulD if GMX_DOUBLE is set, otherwise \ref simdMulF.
 *
 * \copydetails simdMulF
 */
static inline SimdReal
simdMul(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdMulD(a, b);
#else
    return simdMulF(a, b);
#endif
}

/*! \brief SIMD a*b+c for three \ref SimdReal.
 *
 * Uses \ref simdFmaddD if GMX_DOUBLE is set, otherwise \ref simdFmaddF.
 *
 * \copydetails simdFmaddF
 */
static inline SimdReal
simdFmadd(SimdReal a, SimdReal b, SimdReal c)
{
#ifdef GMX_DOUBLE
    return simdFmaddD(a, b, c);
#else
    return simdFmaddF(a, b, c);
#endif
}

/*! \brief SIMD a*b-c for three \ref SimdReal.
 *
 * Uses \ref simdFmsubD if GMX_DOUBLE is set, otherwise \ref simdFmsubF.
 *
 * \copydetails simdFmsubF
 */
static inline SimdReal
simdFmsub(SimdReal a, SimdReal b, SimdReal c)
{
#ifdef GMX_DOUBLE
    return simdFmsubD(a, b, c);
#else
    return simdFmsubF(a, b, c);
#endif
}

/*! \brief SIMD -a*b+c for three \ref SimdReal.
 *
 * Uses \ref simdFnmaddD if GMX_DOUBLE is set, otherwise \ref simdFnmaddF.
 *
 * \copydetails simdFnmaddF
 */
static inline SimdReal
simdFnmadd(SimdReal a, SimdReal b, SimdReal c)
{
#ifdef GMX_DOUBLE
    return simdFnmaddD(a, b, c);
#else
    return simdFnmaddF(a, b, c);
#endif
}

/*! \brief SIMD -a*b-c for three \ref SimdReal.
 *
 * Uses \ref simdFnmsubD if GMX_DOUBLE is set, otherwise \ref simdFnmsubF.
 *
 * \copydetails simdFnmsubF
 */
static inline SimdReal
simdFnmsub(SimdReal a, SimdReal b, SimdReal c)
{
#ifdef GMX_DOUBLE
    return simdFnmsubD(a, b, c);
#else
    return simdFnmsubF(a, b, c);
#endif
}

/*! \brief SIMD table lookup for 1/sqrt(x) approximation.
 *
 * Uses \ref simdRsqrtD if GMX_DOUBLE is set, otherwise \ref simdRsqrtF.
 *
 * \copydetails simdRsqrtF
 */
static inline SimdReal
simdRsqrt(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdRsqrtD(x);
#else
    return simdRsqrtF(x);
#endif
}

/*! \brief SIMD table lookup for 1/x approximation.
 *
 * Uses \ref simdRcpD if GMX_DOUBLE is set, otherwise \ref simdRcpF.
 *
 * \copydetails simdRcpF
 */
static inline SimdReal
simdRcp(SimdReal x)
{
#ifdef GMX_DOUBLE
    return simdRcpD(x);
#else
    return simdRcpF(x);
#endif
}

/*! \brief SIMD fabs(x) for \ref SimdReal.
 *
 * Uses \ref simdAbsD if GMX_DOUBLE is set, otherwise \ref simdAbsF.
 *
 * \copydetails simdAbsF
 */
static inline SimdReal
simdAbs(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdAbsD(a);
#else
    return simdAbsF(a);
#endif
}

/*! \brief SIMD -x for \ref SimdReal.
 *
 * Uses \ref simdNegD if GMX_DOUBLE is set, otherwise \ref simdNegF.
 *
 * \copydetails simdNegF
 */
static inline SimdReal
simdNeg(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdNegD(a);
#else
    return simdNegF(a);
#endif
}

/*! \brief SIMD max(a,b) for each element in \ref SimdReal.
 *
 * Uses \ref simdMaxD if GMX_DOUBLE is set, otherwise \ref simdMaxF.
 *
 * \copydetails simdMaxF
 */
static inline SimdReal
simdMax(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdMaxD(a, b);
#else
    return simdMaxF(a, b);
#endif
}

/*! \brief SIMD min(a,b) for each element in \ref SimdReal.
 *
 * Uses \ref simdMinD if GMX_DOUBLE is set, otherwise \ref simdMinF.
 *
 * \copydetails simdMinF
 */
static inline SimdReal
simdMin(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdMinD(a, b);
#else
    return simdMinF(a, b);
#endif
}

/*! \brief Round \ref SimdReal to nearest int, return \ref SimdReal.
 *
 * Uses \ref simdRoundD if GMX_DOUBLE is set, otherwise \ref simdRoundF.
 *
 * \copydetails simdRoundF
 */
static inline SimdReal
simdRound(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdRoundD(a);
#else
    return simdRoundF(a);
#endif
}

/*! \brief Truncate \ref SimdReal towards 0, return \ref SimdReal.
 *
 * Uses \ref simdTruncD if GMX_DOUBLE is set, otherwise \ref simdTruncF.
 *
 * \copydetails simdTruncF
 */
static inline SimdReal
simdTrunc(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdTruncD(a);
#else
    return simdTruncF(a);
#endif
}

/*! \brief SIMD Fraction, i.e. x-trunc(x) for \ref SimdReal.
 *
 * Uses \ref simdFractionD if GMX_DOUBLE is set, otherwise \ref simdFractionF.
 *
 * \copydetails simdFractionF
 */
static inline SimdReal
simdFraction(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdFractionD(a);
#else
    return simdFractionF(a);
#endif
}

/*! \brief Return the FP exponent of a SIMD \ref SimdReal as a \ref SimdReal.
 *
 * Uses \ref simdGetExponentD if GMX_DOUBLE is set, otherwise \ref simdGetExponentF.
 *
 * \copydetails simdGetExponentF
 */
static inline SimdReal
simdGetExponent(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdGetExponentD(a);
#else
    return simdGetExponentF(a);
#endif
}

/*! \brief Return the FP mantissa of a SIMD \ref SimdReal as a \ref SimdReal.
 *
 * Uses \ref simdGetMantissaD if GMX_DOUBLE is set, otherwise \ref simdGetMantissaF.
 *
 * \copydetails simdGetMantissaF
 */
static inline SimdReal
simdGetMantissa(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdGetMantissaD(a);
#else
    return simdGetMantissaF(a);
#endif
}

/*! \brief Set the exponent of a SIMD \ref SimdReal from a \ref SimdReal.
 *
 * Uses \ref simdSetExponentD if GMX_DOUBLE is set, otherwise \ref simdSetExponentF.
 *
 * \copydetails simdSetExponentF
 */
static inline SimdReal
simdSetExponent(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdSetExponentD(a);
#else
    return simdSetExponentF(a);
#endif
}

/*! \} end of name-group describing FP arithmetics */


/*! \name SIMD comparison, boolean, and select operations for SimdReal
 *  \{
 */

/*! \brief SIMD a==b for \ref SimdReal. Returns a \ref SimdBool.
 *
 * Uses \ref simdCmpEqD if GMX_DOUBLE is set, otherwise \ref simdCmpEqF.
 *
 * \copydetails simdCmpEqF
 */
static inline SimdBool
simdCmpEq(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdCmpEqD(a, b);
#else
    return simdCmpEqF(a, b);
#endif
}

/*! \brief SIMD a<b for \ref SimdReal. Returns a \ref SimdBool.
 *
 * Uses \ref simdCmpLtD if GMX_DOUBLE is set, otherwise \ref simdCmpLtF.
 *
 * \copydetails simdCmpLtF
 */
static inline SimdBool
simdCmpLt(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdCmpLtD(a, b);
#else
    return simdCmpLtF(a, b);
#endif
}

/*! \brief SIMD a<=b for \ref SimdReal. Returns a \ref SimdBool.
 *
 * Uses \ref simdCmpLeD if GMX_DOUBLE is set, otherwise \ref simdCmpLeF.
 *
 * \copydetails simdCmpLeF
 */
static inline SimdBool
simdCmpLe(SimdReal a, SimdReal b)
{
#ifdef GMX_DOUBLE
    return simdCmpLeD(a, b);
#else
    return simdCmpLeF(a, b);
#endif
}

/*! \brief For each element, the result boolean is true if both arguments are true
 *
 * Uses \ref simdAndDB if GMX_DOUBLE is set, otherwise \ref simdAndFB.
 *
 * \copydetails simdAndFB
 */
static inline SimdBool
simdAndB(SimdBool a, SimdBool b)
{
#ifdef GMX_DOUBLE
    return simdAndDB(a, b);
#else
    return simdAndFB(a, b);
#endif
}

/*! \brief For each element, the result boolean is true if either argument is true
 *
 * Uses \ref simdOrDB if GMX_DOUBLE is set, otherwise \ref simdOrFB.
 *
 * \copydetails simdOrFB
 */
static inline SimdBool
simdOrB(SimdBool a, SimdBool b)
{
#ifdef GMX_DOUBLE
    return simdOrDB(a, b);
#else
    return simdOrFB(a, b);
#endif
}

/*! \brief Return nonzero if any element in SimdBool is true, otherwise 0.
 *
 * Uses \ref simdAnyTrueDB if GMX_DOUBLE is set, otherwise \ref simdAnyTrueFB.
 *
 * \copydetails simdAnyTrueFB
 */
static inline int
simdAnyTrueB(SimdBool a)
{
#ifdef GMX_DOUBLE
    return simdAnyTrueDB(a);
#else
    return simdAnyTrueFB(a);
#endif
}

/*! \brief Selects elements from \ref SimdReal where boolean is true, otherwise 0.
 *
 * Uses \ref simdMaskD if GMX_DOUBLE is set, otherwise \ref simdMaskF.
 *
 * \copydetails simdMaskF
 *
 * \sa simdMaskI
 */
static inline SimdReal
simdMask(SimdReal a, SimdBool mask)
{
#ifdef GMX_DOUBLE
    return simdMaskD(a, mask);
#else
    return simdMaskF(a, mask);
#endif
}

/*! \brief Selects elements from \ref SimdReal where boolean is false, otherwise 0.
 *
 * Uses \ref simdMaskNotD if GMX_DOUBLE is set, otherwise \ref simdMaskNotF.
 *
 * \copydetails simdMaskNotF
 */
static inline SimdReal
simdMaskNot(SimdReal a, SimdBool mask)
{
#ifdef GMX_DOUBLE
    return simdMaskNotD(a, mask);
#else
    return simdMaskNotF(a, mask);
#endif
}

/*! \brief Selects from 2nd real SIMD arg where boolean is true, otherwise 1st arg.
 *
 * Uses \ref simdBlendD if GMX_DOUBLE is set, otherwise \ref simdBlendF.
 *
 * \copydetails simdBlendF
 */
static inline SimdReal
simdBlend(SimdReal a, SimdReal b, SimdBool sel)
{
#ifdef GMX_DOUBLE
    return simdBlendD(a, b, sel);
#else
    return simdBlendF(a, b, sel);
#endif
}

/*! \brief Return sum of all elements in SIMD floating-point variable.
 *
 * Uses \ref simdReduceD if GMX_DOUBLE is set, otherwise \ref simdReduceF.
 *
 * \copydetails simdReduceF
 */
static inline real
simdReduce(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdReduceD(a);
#else
    return simdReduceF(a);
#endif
}

/*! \}  end of name-group describing FP comparison and booleans */
#endif // GMX_SIMD_HAVE_REAL

#if GMX_SIMD_HAVE_INT32_LOGICAL
/*! \name SIMD integer logical operations on SimdInt32
 *
 *  These instructions are available if \ref GMX_SIMD_HAVE_INT32_LOGICAL is 1.
 *  \{
 */

/*! \brief Shift each element in \ref SimdInt32 left by immediate
 *
 * Uses \ref simdSlliDI if GMX_DOUBLE is set, otherwise \ref simdSlliFI.
 *
 * \copydetails simdSlliFI
 */
static inline SimdInt32
simdSlliI(SimdInt32 a, int n)
{
#ifdef GMX_DOUBLE
    return simdSlliDI(a, n);
#else
    return simdSlliFI(a, n);
#endif
}

/*! \brief Shift each element in \ref SimdInt32 right by immediate
 *
 * Uses \ref simdSrliDI if GMX_DOUBLE is set, otherwise \ref simdSrliFI.
 *
 * \copydetails simdSrliFI
 */
static inline SimdInt32
simdSrliI(SimdInt32 a, int n)
{
#ifdef GMX_DOUBLE
    return simdSrliDI(a, n);
#else
    return simdSrliFI(a, n);
#endif
}

/*! \brief Bitwise \a and on two \ref SimdInt32.
 *
 * Uses \ref simdAndDI if GMX_DOUBLE is set, otherwise \ref simdAndFI.
 *
 * \copydetails simdAndFI
 */
static inline SimdInt32
simdAndI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdAndDI(a, b);
#else
    return simdAndFI(a, b);
#endif
}

/*! \brief Bitwise \a and-not on two \ref SimdInt32; 1st arg is complemented.
 *
 * Uses \ref simdAndNotDI if GMX_DOUBLE is set, otherwise \ref simdAndNotFI.
 *
 * \copydetails simdAndNotFI
 */
static inline SimdInt32
simdAndNotI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdAndNotDI(a, b);
#else
    return simdAndNotFI(a, b);
#endif
}

/*! \brief Bitwise \a or on two \ref SimdInt32.
 *
 * Uses \ref simdOrDI if GMX_DOUBLE is set, otherwise \ref simdOrFI.
 *
 * \copydetails simdOrFI
 */
static inline SimdInt32
simdOrI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdOrDI(a, b);
#else
    return simdOrFI(a, b);
#endif
}

/*! \brief Bitwise \a xor on two \ref SimdInt32.
 *
 * Uses \ref simdXorDI if GMX_DOUBLE is set, otherwise \ref simdXorFI.
 *
 * \copydetails simdXorFI
 */
static inline SimdInt32
simdXorI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdXorDI(a, b);
#else
    return simdXorFI(a, b);
#endif
}

/*! \}   end of name-group describing logical ops on SimdInt32 */
#endif // GMX_SIMD_HAVE_INT32_LOGICAL


#if GMX_SIMD_HAVE_INT32_ARITHMETICS
/*! \name SIMD integer arithmetic operations on SimdInt32
 *
 *  These instructions are available if \ref GMX_SIMD_HAVE_INT32_ARITHMETICS is 1.
 *  \{
 */

/*! \brief SIMD a+b for two \ref SimdInt32.
 *
 * Uses \ref simdAddDI if GMX_DOUBLE is set, otherwise \ref simdAddFI.
 *
 * \copydetails simdAddFI
 */
static inline SimdInt32
simdAddI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdAddDI(a, b);
#else
    return simdAddFI(a, b);
#endif
}

/*! \brief SIMD a-b for two \ref SimdInt32.
 *
 * Uses \ref simdSubDI if GMX_DOUBLE is set, otherwise \ref simdSubFI.
 *
 * \copydetails simdSubFI
 */
static inline SimdInt32
simdSubI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdSubDI(a, b);
#else
    return simdSubFI(a, b);
#endif
}

/*! \brief SIMD a*b for two \ref SimdInt32.
 *
 * Uses \ref simdMulDI if GMX_DOUBLE is set, otherwise \ref simdMulFI.
 *
 * \copydetails simdMulFI
 */
static inline SimdInt32
simdMulI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdMulDI(a, b);
#else
    return simdMulFI(a, b);
#endif
}

/*! \}  end of name-group describing SimdInt32 arithmetics */


/*! \name SIMD integer comparison, booleans, and selection on SimdInt32
 *
 *  These instructions are available if \ref GMX_SIMD_HAVE_INT32_ARITHMETICS is 1.
 *  \{
 */

/*! \brief Returns boolean describing whether a==b, for \ref SimdInt32
 *
 * Uses \ref simdCmpEqDI if GMX_DOUBLE is set, otherwise \ref simdCmpEqFI.
 *
 * \copydetails simdCmpEqFI
 */
static inline SimdIBool
simdCmpEqI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdCmpEqDI(a, b);
#else
    return simdCmpEqFI(a, b);
#endif
}

/*! \brief Returns boolean describing whether a<b, for \ref SimdInt32
 *
 * Uses \ref simdCmpLtDI if GMX_DOUBLE is set, otherwise \ref simdCmpLtFI.
 *
 * \copydetails simdCmpLtFI
 */
static inline SimdIBool
simdCmpLtI(SimdInt32 a, SimdInt32 b)
{
#ifdef GMX_DOUBLE
    return simdCmpLtDI(a, b);
#else
    return simdCmpLtFI(a, b);
#endif
}

/*! \brief For each element, the result boolean is true if both arguments are true
 *
 * Uses \ref simdAndDIB if GMX_DOUBLE is set, otherwise \ref simdAndFIB.
 *
 * \copydetails simdAndFIB
 */
static inline SimdIBool
simdAndIB(SimdIBool a, SimdIBool b)
{
#ifdef GMX_DOUBLE
    return simdAndDIB(a, b);
#else
    return simdAndFIB(a, b);
#endif
}

/*! \brief For each element, the result boolean is true if either argument is true.
 *
 * Uses \ref simdOrDIB if GMX_DOUBLE is set, otherwise \ref simdOrFIB.
 *
 * \copydetails simdOrFIB
 */
static inline SimdIBool
simdOrIB(SimdIBool a, SimdIBool b)
{
#ifdef GMX_DOUBLE
    return simdOrDIB(a, b);
#else
    return simdOrFIB(a, b);
#endif
}

/*! \brief Return nonzero if any element in SimdIBool is true, otherwise 0.
 *
 * Uses \ref simdAnyTrueDIB if GMX_DOUBLE is set, otherwise \ref simdAnyTrueFIB.
 *
 * \copydetails simdAnyTrueFIB
 */
static inline int
simdAnyTrueIB(SimdIBool a)
{
#ifdef GMX_DOUBLE
    return simdAnyTrueDIB(a);
#else
    return simdAnyTrueFIB(a);
#endif
}

/*! \brief Selects elements from \ref SimdInt32 where boolean is true, otherwise 0.
 *
 * Uses \ref simdMaskDI if GMX_DOUBLE is set, otherwise \ref simdMaskFI.
 *
 * \copydetails simdMaskFI
 */
static inline SimdInt32
simdMaskI(SimdInt32 a, SimdIBool mask)
{
#ifdef GMX_DOUBLE
    return simdMaskDI(a, mask);
#else
    return simdMaskFI(a, mask);
#endif
}

/*! \brief Selects elements from \ref SimdInt32 where boolean is false, otherwise 0.
 *
 * Uses \ref simdMaskNotDI if GMX_DOUBLE is set, otherwise \ref simdMaskNotFI.
 *
 * \copydetails simdMaskNotFI
 */
static inline SimdInt32
simdMaskNotI(SimdInt32 a, SimdIBool mask)
{
#ifdef GMX_DOUBLE
    return simdMaskNotDI(a, mask);
#else
    return simdMaskNotFI(a, mask);
#endif
}

/*! \brief Selects from 2nd int SIMD arg where boolean is true, otherwise 1st arg.
 *
 * Uses \ref simdBlendDI if GMX_DOUBLE is set, otherwise \ref simdBlendFI.
 *
 * \copydetails simdBlendFI
 */
static inline SimdInt32
simdBlendI(SimdInt32 a, SimdInt32 b, SimdIBool sel)
{
#ifdef GMX_DOUBLE
    return simdBlendDI(a, b, sel);
#else
    return simdBlendFI(a, b, sel);
#endif
}

/*! \}  end of name-group describing SimdInt32 comparisons and booleans */
#endif // GMX_SIMD_HAVE_INT32_ARITHMETICS

#if GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32
/*! \name SIMD conversion operations
 *
 *  These instructions are available when both types involved in the conversion
 *  are defined, e.g. if \ref GMX_SIMD_HAVE_REAL and \ref GMX_SIMD_HAVE_INT32
 *  are 1 for real-to-integer conversion.
 *  \{
 */

/*! \brief Convert SimdReal to SimdInt32, round to nearest integer.
 *
 * Uses \ref simdCvtD2I if GMX_DOUBLE is set, otherwise \ref simdCvtF2I.
 *
 * \copydetails simdCvtF2I
 */
static inline SimdInt32
simdCvtR2I(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdCvtD2I(a);
#else
    return simdCvtF2I(a);
#endif
}

/*! \brief Convert SimdReal to SimdInt32, truncate towards zero
 *
 * Uses \ref simdCvttD2I if GMX_DOUBLE is set, otherwise \ref simdCvttF2I.
 *
 * \copydetails simdCvttF2I
 */
static inline SimdInt32
simdCvttR2I(SimdReal a)
{
#ifdef GMX_DOUBLE
    return simdCvttD2I(a);
#else
    return simdCvttF2I(a);
#endif
}

/*! \brief Convert SimdInt32 to SimdReal
 *
 * Uses \ref simdCvtI2D if GMX_DOUBLE is set, otherwise \ref simdCvtI2F.
 *
 * \copydetails simdCvtI2F
 */
static inline SimdReal
simdCvtI2R(SimdInt32 a)
{
#ifdef GMX_DOUBLE
    return simdCvtI2D(a);
#else
    return simdCvtI2F(a);
#endif
}

/*! \brief Convert from SimdBool to SimdIBool
 *
 * Uses \ref simdCvtDB2DIB if GMX_DOUBLE is set, otherwise \ref simdCvtFB2FIB.
 *
 * \copydetails simdCvtFB2FIB
 */
static inline SimdIBool
simdCvtB2IB(SimdBool a)
{
#ifdef GMX_DOUBLE
    return simdCvtDB2DIB(a);
#else
    return simdCvtFB2FIB(a);
#endif
}

/*! \brief Convert from SimdIBool to SimdBool
 *
 * Uses \ref simdCvtDIB2DB if GMX_DOUBLE is set, otherwise \ref simdCvtFIB2FB.
 *
 * \copydetails simdCvtFIB2FB
 */
static inline SimdBool
simdCvtIB2B(SimdIBool a)
{
#ifdef GMX_DOUBLE
    return simdCvtDIB2DB(a);
#else
    return simdCvtFIB2FB(a);
#endif
}

/*! \}    end of name-group describing SIMD conversions */
#endif // GMX_SIMD_HAVE_REAL && GMX_SIMD_HAVE_INT32


#if GMX_SIMD4_HAVE_REAL
/*! \name SIMD4 - constant width-four SIMD datatypes
 *
 * These operations are only meant to be used for a few coordinate
 * manipulation and grid interpolation routines, so we only support a subset
 * of operations for SIMD4. To avoid repeating all the documentation from
 * the generic width SIMD routines, we only provide brief documentation for
 * these operations. Follow the link to the implementation documentation or the
 * reference to the corresponding generic SIMD routine. The format will be
 * exactly the same, but they have SIMD replaced with SIMD4.
 *  \{
 */

/*! \brief SIMD real datatype guaranteed to be 4 elements wide, if available.
 *
 * All the SIMD4 datatypes and operations behave like their counterparts for
 * the generic SIMD implementation, but they might be implemented with different
 * registers, or not supported at all. It is important that you check the
 * define \ref GMX_SIMD4_HAVE_REAL before using it.
 *
 * Just as the normal SIMD operations, all SIMD4 types and routines will
 * be aliased to either single or double precision ones based on whether
 * GMX_DOUBLE is defined.
 *
 * \note There is no support for integer or math operations in SIMD4.
 */
#    ifdef GMX_DOUBLE
typedef Simd4Double               Simd4Real;
#    else
typedef Simd4Float                Simd4Real;
#    endif

/*! \brief Boolean for \ref Simd4Real comparision/selection */
#    ifdef GMX_DOUBLE
typedef Simd4DBool                Simd4Bool;
#    else
typedef Simd4FBool                Simd4Bool;
#    endif


/*! \brief Load aligned data to Simd4Real.
 *
 * \copydetails simd4LoadF
 */
static inline Simd4Real
simd4Load(const real *m)
{
#ifdef GMX_DOUBLE
    return simd4LoadD(m);
#else
    return simd4LoadF(m);
#endif
}

/*! \brief Load single element to Simd4Real
 *
 * \copydetails simd4Load1F
 */
static inline Simd4Real
simd4Load1(const real *m)
{
#ifdef GMX_DOUBLE
    return simd4Load1D(m);
#else
    return simd4Load1F(m);
#endif
}

/*! \brief Set Simd4Real from scalar value
 *
 * \copydetails simd4Set1F
 */
static inline Simd4Real
simd4Set1(const real r)
{
#ifdef GMX_DOUBLE
    return simd4Set1D(r);
#else
    return simd4Set1F(r);
#endif
}

/*! \brief store aligned data from Simd4Real
 *
 * \copydetails simd4StoreF
 */
static inline void
simd4Store(real *m, Simd4Real a)
{
#ifdef GMX_DOUBLE
    simd4StoreD(m, a);
#else
    simd4StoreF(m, a);
#endif
}

#if GMX_SIMD_HAVE_LOADU
/*! \brief Load unaligned data to Simd4Real
 *
 * \copydetails simd4LoadUF
 */
static inline Simd4Real
simd4LoadU(const real *m)
{
#ifdef GMX_DOUBLE
    return simd4LoadUD(m);
#else
    return simd4LoadUF(m);
#endif
}
#endif // GMX_SIMD_HAVE_LOADU


#if GMX_SIMD_HAVE_STOREU
/*! \brief Store unaligned data from Simd4Real
 *
 * \copydetails simd4StoreUF
 */
static inline void
simd4StoreU(real *m, Simd4Real a)
{
#ifdef GMX_DOUBLE
    simd4StoreUD(m, a);
#else
    simd4StoreUF(m, a);
#endif
}
#endif // GMX_SIMD_HAVE_STOREU


/*! \brief Set all elements in Simd4Real to 0.0
 *
 * \copydetails simd4SetZeroF
 */
static inline Simd4Real
simd4SetZero()
{
#ifdef GMX_DOUBLE
    return simd4SetZeroD();
#else
    return simd4SetZeroF();
#endif
}

/*! \brief Bitwise and for two Simd4Real
 *
 * \copydetails simd4AndF
 */
static inline Simd4Real
simd4And(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4AndD(a, b);
#else
    return simd4AndF(a, b);
#endif
}

/*! \brief Bitwise and-not for two Simd4Real. 1st arg is complemented.
 *
 * \copydetails simd4AndNotF
 */
static inline Simd4Real
simd4AndNot(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4AndNotD(a, b);
#else
    return simd4AndNotF(a, b);
#endif
}

/*! \brief Bitwise or for two Simd4Real
 *
 * \copydetails simd4OrF
 */
static inline Simd4Real
simd4Or(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4OrD(a, b);
#else
    return simd4OrF(a, b);
#endif
}

/*! \brief Bitwise xor for two Simd4Real
 *
 * \copydetails simd4XorF
 */
static inline Simd4Real
simd4Xor(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4XorD(a, b);
#else
    return simd4XorF(a, b);
#endif
}

/*! \brief a+b for \ref Simd4Real
 *
 * \copydetails simd4AddF
 */
static inline Simd4Real
simd4Add(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4AddD(a, b);
#else
    return simd4AddF(a, b);
#endif
}

/*! \brief a-b for \ref Simd4Real
 *
 * \copydetails simd4SubF
 */
static inline Simd4Real
simd4Sub(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4SubD(a, b);
#else
    return simd4SubF(a, b);
#endif
}

/*! \brief a*b for \ref Simd4Real
 *
 * \copydetails simd4MulF
 */
static inline Simd4Real
simd4Mul(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4MulD(a, b);
#else
    return simd4MulF(a, b);
#endif
}

/*! \brief a*b+c for \ref Simd4Real
 *
 * \copydetails simd4FmaddF
 */
static inline Simd4Real
simd4Fmadd(Simd4Real a, Simd4Real b, Simd4Real c)
{
#ifdef GMX_DOUBLE
    return simd4FmaddD(a, b, c);
#else
    return simd4FmaddF(a, b, c);
#endif
}

/*! \brief a*b-c for \ref Simd4Real
 *
 * \copydetails simd4FmsubF
 */
static inline Simd4Real
simd4Fmsub(Simd4Real a, Simd4Real b, Simd4Real c)
{
#ifdef GMX_DOUBLE
    return simd4FmsubD(a, b, c);
#else
    return simd4FmsubF(a, b, c);
#endif
}

/*! \brief -a*b+c for \ref Simd4Real
 *
 * \copydetails simd4FnmaddF
 */
static inline Simd4Real
simd4Fnmadd(Simd4Real a, Simd4Real b, Simd4Real c)
{
#ifdef GMX_DOUBLE
    return simd4FnmaddD(a, b, c);
#else
    return simd4FnmaddF(a, b, c);
#endif
}

/*! \brief -a*b-c for \ref Simd4Real
 *
 * \copydetails simd4FnmsubF
 */
static inline Simd4Real
simd4Fnmsub(Simd4Real a, Simd4Real b, Simd4Real c)
{
#ifdef GMX_DOUBLE
    return simd4FnmsubD(a, b, c);
#else
    return simd4FnmsubF(a, b, c);
#endif
}

/*! \brief 1/sqrt(x) approximate lookup for \ref Simd4Real
 *
 * \copydetails simd4RsqrtF
 */
static inline Simd4Real
simd4Rsqrt(Simd4Real x)
{
#ifdef GMX_DOUBLE
    return simd4RsqrtD(x);
#else
    return simd4RsqrtF(x);
#endif
}

/*! \brief fabs(x) for \ref Simd4Real
 *
 * \copydetails simd4AbsF
 */
static inline Simd4Real
simd4Abs(Simd4Real a)
{
#ifdef GMX_DOUBLE
    return simd4AbsD(a);
#else
    return simd4AbsF(a);
#endif
}

/*! \brief Change sign (-x) for \ref Simd4Real
 *
 * \copydetails simd4NegF
 */
static inline Simd4Real
simd4Neg(Simd4Real a)
{
#ifdef GMX_DOUBLE
    return simd4NegD(a);
#else
    return simd4NegF(a);
#endif
}

/*! \brief Select maximum of each pair of elements from args for \ref Simd4Real
 *
 * \copydetails simd4MaxF
 */
static inline Simd4Real
simd4Max(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4MaxD(a, b);
#else
    return simd4MaxF(a, b);
#endif
}

/*! \brief Select minimum of each pair of elements from args for \ref Simd4Real
 *
 * \copydetails simd4MinF
 */
static inline Simd4Real
simd4Min(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4MinD(a, b);
#else
    return simd4MinF(a, b);
#endif
}

/*! \brief Round \ref Simd4Real to nearest integer, return \ref Simd4Real
 *
 * \copydetails simd4RoundF
 */
static inline Simd4Real
simd4Round(Simd4Real a)
{
#ifdef GMX_DOUBLE
    return simd4RoundD(a);
#else
    return simd4RoundF(a);
#endif
}

/*! \brief Truncate \ref Simd4Real towards zero, return \ref Simd4Real
 *
 * \copydetails simd4TruncF
 */
static inline Simd4Real
simd4Trunc(Simd4Real a)
{
#ifdef GMX_DOUBLE
    return simd4TruncD(a);
#else
    return simd4TruncF(a);
#endif
}

/*! \brief Scalar product of first three elements of two \ref Simd4Real
 *
 * \copydetails simd4DotProductF
 */
static inline real
simd4DotProduct(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4DotProductD(a, b);
#else
    return simd4DotProductF(a, b);
#endif
}

/*! \brief Return booleans whether a==b for each element two \ref Simd4Real
 *
 * \copydetails simd4CmpEqF
 */
static inline Simd4Bool
simd4CmpEq(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4CmpEqD(a, b);
#else
    return simd4CmpEqF(a, b);
#endif
}

/*! \brief Return booleans whether a<b for each element two \ref Simd4Real
 *
 * \copydetails simd4CmpLtF
 */
static inline Simd4Bool
simd4CmpLt(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4CmpLtD(a, b);
#else
    return simd4CmpLtF(a, b);
#endif
}

/*! \brief Return booleans whether a<=b for each element two \ref Simd4Real
 *
 * \copydetails simd4CmpLeF
 */
static inline Simd4Bool
simd4CmpLe(Simd4Real a, Simd4Real b)
{
#ifdef GMX_DOUBLE
    return simd4CmpLeD(a, b);
#else
    return simd4CmpLeF(a, b);
#endif
}

/*! \brief Logical and for two \ref gmx::Simd4Bool
 *
 * \copydetails simd4AndFB
 */
static inline Simd4Bool
simd4AndB(Simd4Bool a, Simd4Bool b)
{
#ifdef GMX_DOUBLE
    return simd4AndDB(a, b);
#else
    return simd4AndFB(a, b);
#endif
}

/*! \brief Logical or for two \ref gmx::Simd4Bool
 *
 * \copydetails simd4OrFB
 */
static inline Simd4Bool
simd4OrB(Simd4Bool a, Simd4Bool b)
{
#ifdef GMX_DOUBLE
    return simd4OrDB(a, b);
#else
    return simd4OrFB(a, b);
#endif
}

/*! \brief Return nonzero if any element in \ref gmx::Simd4Bool is true, otherwise 0
 *
 * \copydetails simd4AnyTrueFB
 */
static inline int
simd4AnyTrueB(Simd4Bool a)
{
#ifdef GMX_DOUBLE
    return simd4AnyTrueDB(a);
#else
    return simd4AnyTrueFB(a);
#endif
}

/*! \brief Selects from 2nd real SIMD4 arg where boolean is true, otherwise 1st arg
 *
 * \copydetails simd4MaskF
 */
static inline Simd4Real
simd4Mask(Simd4Real a, Simd4Bool mask)
{
#ifdef GMX_DOUBLE
    return simd4MaskD(a, mask);
#else
    return simd4MaskF(a, mask);
#endif
}

/*! \brief Selects from 2nd real SIMD4 arg where boolean is false, otherwise 1st arg
 *
 * \copydetails simd4MaskNotF
 */
static inline Simd4Real
simd4MaskNot(Simd4Real a, Simd4Bool mask)
{
#ifdef GMX_DOUBLE
    return simd4MaskNotD(a, mask);
#else
    return simd4MaskNotF(a, mask);
#endif
}

/*! \brief Selects from 2nd real SIMD4 arg where boolean is true, otherwise 1st arg
 *
 * \copydetails simd4BlendF
 */
static inline Simd4Real
simd4Blend(Simd4Real a, Simd4Real b, Simd4Bool sel)
{
#ifdef GMX_DOUBLE
    return simd4BlendD(a, b, sel);
#else
    return simd4BlendF(a, b, sel);
#endif
}

/*! \brief Return sum of all elements in SIMD4 floating-point variable.
 *
 * \copydetails simd4ReduceF
 */
static inline real
simd4Reduce(Simd4Real a)
{
#ifdef GMX_DOUBLE
    return simd4ReduceD(a);
#else
    return simd4ReduceF(a);
#endif
}

/*! \}   end of name-group describing SIMD4 */
#endif // GMX_SIMD4_HAVE_REAL

}      // namespace gmx

/*! \}  end of addtogroup module_simd */


/*! \endcond   end of condition libapi */


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

#endif /* GMX_SIMD_SIMD_H */
