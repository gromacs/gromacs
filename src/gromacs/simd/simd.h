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

// SIMD architectures other than reference are temporarily disabled in this commit
#if (GMX_SIMD_REFERENCE || defined DOXYGEN)
// Plain C SIMD reference implementation, also serves as documentation.
#    include "impl_reference/impl_reference.h"
#else
#    include "impl_none/impl_none.h"
#endif

#ifdef GMX_DOUBLE
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
 *  \ref GMX_SIMD_HAVE_DOUBLE if GMX_DOUBLE is set, otherwise \ref GMX_SIMD_HAVE_FLOAT.
 */
#    define GMX_SIMD_HAVE_REAL               GMX_SIMD_HAVE_FLOAT

/*! \brief Width of SimdReal.
 *
 *  \ref GMX_SIMD_DOUBLE_WIDTH if GMX_DOUBLE is set, otherwise \ref GMX_SIMD_FLOAT_WIDTH.
 */
#    define GMX_SIMD_REAL_WIDTH              GMX_SIMD_FLOAT_WIDTH

/*! \brief 1 if support is available for extracting elements from SimdInt32, otherwise 0
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

/*! \brief 1 if gmx::simdGatherLoadUBySimdIntTranspose is present, otherwise 0
 *
 *  \ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_DOUBLE if GMX_DOUBLE is set, otherwise
 *  \ref GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT.
 */
#    define GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_REAL    GMX_SIMD_HAVE_GATHER_LOADU_BYSIMDINT_TRANSPOSE_FLOAT

/*! \brief 1 if real half-register load/store/reduce utils present, otherwise 0
 *
 *  \ref GMX_SIMD_HAVE_HSIMD_UTIL_DOUBLE if GMX_DOUBLE is set, otherwise
 *  \ref GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT.
 */
#    define GMX_SIMD_HAVE_HSIMD_UTIL_REAL    GMX_SIMD_HAVE_HSIMD_UTIL_FLOAT

/*! \brief 1 if Simd4Real is available, otherwise 0.
 *
 *  \ref GMX_SIMD4_HAVE_DOUBLE if GMX_DOUBLE is set, otherwise \ref GMX_SIMD4_HAVE_FLOAT.
 */
#    define GMX_SIMD4_HAVE_REAL              GMX_SIMD4_HAVE_FLOAT

#endif // GMX_DOUBLE

/*! \}  end of name-group describing high-level capabilities */

namespace gmx
{

namespace
{

#if GMX_SIMD
class SimdLoadIProxy
{
    public:
        SimdLoadIProxy(const std::int32_t *m) : m_(m) {}
#if GMX_SIMD_HAVE_FLOAT
        operator SimdFInt32() const { return loadFI(m_); };
#endif
#if GMX_SIMD_HAVE_DOUBLE
        operator SimdDInt32() const { return loadDI(m_); };
#endif

    private:
        const std::int32_t * const m_;
};

const SimdLoadIProxy
load(const std::int32_t *m)
{
    return SimdLoadIProxy(m);
}

#if GMX_SIMD_HAVE_LOADU
class SimdLoadUIProxy
{
    public:
        SimdLoadUIProxy(const std::int32_t *m) : m_(m) {}
#if GMX_SIMD_HAVE_FLOAT
        operator SimdFInt32() const { return loadUFI(m_); };
#endif
#if GMX_SIMD_HAVE_DOUBLE
        operator SimdDInt32() const { return loadUDI(m_); };
#endif

    private:
        const std::int32_t * const m_;
};

const SimdLoadUIProxy
loadU(const std::int32_t *m)
{
    return SimdLoadUIProxy(m);
}
#endif // GMX_SIMD_HAVE_LOADU

#endif // GMX_SIMD


// The if-block for HAVE_REAL spans all the way until the SIMD4 definitions a long way down
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


/*! \brief 32-bit integer SIMD type.
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

#if GMX_SIMD_HAVE_INT32_ARITHMETICS
/*! \brief Boolean SIMD type for usage with \ref gmx::SimdInt32.
 *
 * This type is only available if \ref GMX_SIMD_HAVE_FINT32_ARITHMETICS is 1.
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
#endif  // GMX_SIMD_HAVE_INT32_ARITHMETICS

/*! \}  end of name-group describing SIMD data types */

#ifdef GMX_DOUBLE
const int c_simdBestPairAlignment = c_simdBestPairAlignmentD;
#else
const int c_simdBestPairAlignment = c_simdBestPairAlignmentF;
#endif

#endif // GMX_SIMD_HAVE_REAL

}      // namespace anonymous

}      // namespace gmx

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
