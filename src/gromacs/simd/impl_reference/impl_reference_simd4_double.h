/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_SIMD_IMPL_REFERENCE_SIMD4_DOUBLE_H
#define GMX_SIMD_IMPL_REFERENCE_SIMD4_DOUBLE_H

/*! \libinternal \file
 *
 * \brief Reference implementation, SIMD4 single precision.
 *
 * \author Erik Lindahl <erik.lindahl@scilifelab.se>
 *
 * \ingroup module_simd
 */

#include "config.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include <algorithm>
#include <array>

#include "impl_reference_definitions.h"

namespace gmx
{

/*! \cond libapi */
/*! \addtogroup module_simd */
/*! \{ */

/*! \name Constant width-4 double precision SIMD types and instructions
 * \{
 */

/*! \libinternal \brief SIMD4 double type.
 *
 * Available if \ref GMX_SIMD4_HAVE_DOUBLE is 1.
 *
 * \note This variable cannot be placed inside other structures or classes, since
 *       some compilers (including at least clang-3.7) appear to lose the
 *       alignment. This is likely particularly severe when allocating such
 *       memory on the heap, but it occurs for stack structures too.
 */
class Simd4Double
{
    public:
        Simd4Double() {}

        //! \brief Construct from scalar
        Simd4Double(double d) { simdInternal_.fill(d); }

        /*! \brief Internal SIMD data. Implementation dependent, don't touch.
         *
         * This has to be public to enable usage in combination with static inline
         * functions, but it should never, EVER, be accessed by any code outside
         * the corresponding implementation directory since the type will depend
         * on the architecture.
         */
        std::array<double, GMX_SIMD4_WIDTH>  simdInternal_;
};

/*! \libinternal  \brief SIMD4 variable type to use for logical comparisons on doubles.
 *
 * Available if \ref GMX_SIMD4_HAVE_DOUBLE is 1.
 *
 * \note This variable cannot be placed inside other structures or classes, since
 *       some compilers (including at least clang-3.7) appear to lose the
 *       alignment. This is likely particularly severe when allocating such
 *       memory on the heap, but it occurs for stack structures too.
 */
class Simd4DBool
{
    public:
        Simd4DBool() {}

        //! \brief Construct from scalar
        Simd4DBool(bool b) { simdInternal_.fill(b); }

        /*! \brief Internal SIMD data. Implementation dependent, don't touch.
         *
         * This has to be public to enable usage in combination with static inline
         * functions, but it should never, EVER, be accessed by any code outside
         * the corresponding implementation directory since the type will depend
         * on the architecture.
         */
        std::array<bool, GMX_SIMD4_WIDTH>  simdInternal_;
};

/*! \brief Load 4 double values from aligned memory into SIMD4 variable.
 *
 * \param m Pointer to memory aligned to 4 elements.
 * \return SIMD4 variable with data loaded.
 */
static inline Simd4Double gmx_simdcall
load4(const double *m)
{
    Simd4Double a;

    assert(std::size_t(m) % (a.simdInternal_.size()*sizeof(double)) == 0);

    std::copy(m, m+a.simdInternal_.size(), a.simdInternal_.begin());
    return a;
}

/*! \brief Store the contents of SIMD4 double to aligned memory m.
 *
 * \param[out] m Pointer to memory, aligned to 4 elements.
 * \param a SIMD4 variable to store
 */
static inline void gmx_simdcall
store4(double *m, Simd4Double a)
{
    assert(std::size_t(m) % (a.simdInternal_.size()*sizeof(double)) == 0);

    std::copy(a.simdInternal_.begin(), a.simdInternal_.end(), m);
}

/*! \brief Load SIMD4 double from unaligned memory.
 *
 * Available if \ref GMX_SIMD_HAVE_LOADU is 1.
 *
 * \param m Pointer to memory, no alignment requirement.
 * \return SIMD4 variable with data loaded.
 */
static inline Simd4Double gmx_simdcall
load4U(const double *m)
{
    Simd4Double a;
    std::copy(m, m+a.simdInternal_.size(), a.simdInternal_.begin());
    return a;
}

/*! \brief Store SIMD4 double to unaligned memory.
 *
 * Available if \ref GMX_SIMD_HAVE_STOREU is 1.
 *
 * \param[out] m Pointer to memory, no alignment requirement.
 * \param a SIMD4 variable to store.
 */
static inline void gmx_simdcall
store4U(double *m, Simd4Double a)
{
    std::copy(a.simdInternal_.begin(), a.simdInternal_.end(), m);
}

/*! \brief Set all SIMD4 double elements to 0.
 *
 * You should typically just call \ref gmx::setZero(), which uses proxy objects
 * internally to handle all types rather than adding the suffix used here.
 *
 * \return SIMD4 0.0
 */
static inline Simd4Double gmx_simdcall
simd4SetZeroD()
{
    return Simd4Double(0.0);
}


/*! \brief Bitwise and for two SIMD4 double variables.
 *
 * Supported if \ref GMX_SIMD_HAVE_LOGICAL is 1.
 *
 * \param a data1
 * \param b data2
 * \return data1 & data2
 */
static inline Simd4Double gmx_simdcall
operator&(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    union
    {
        double        r;
        std::int64_t  i;
    }
    conv1, conv2;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        conv1.r              = a.simdInternal_[i];
        conv2.r              = b.simdInternal_[i];
        conv1.i              = conv1.i & conv2.i;
        res.simdInternal_[i] = conv1.r;
    }
    return res;
}


/*! \brief Bitwise andnot for two SIMD4 double variables. c=(~a) & b.
 *
 * Available if \ref GMX_SIMD_HAVE_LOGICAL is 1.
 *
 * \param a data1
 * \param b data2
 * \return (~data1) & data2
 */
static inline Simd4Double gmx_simdcall
andNot(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    union
    {
        double        r;
        std::int64_t  i;
    }
    conv1, conv2;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        conv1.r              = a.simdInternal_[i];
        conv2.r              = b.simdInternal_[i];
        conv1.i              = ~conv1.i & conv2.i;
        res.simdInternal_[i] = conv1.r;
    }
    return res;
}


/*! \brief Bitwise or for two SIMD4 doubles.
 *
 * Available if \ref GMX_SIMD_HAVE_LOGICAL is 1.
 *
 * \param a data1
 * \param b data2
 * \return data1 | data2
 */
static inline Simd4Double gmx_simdcall
operator|(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    union
    {
        double        r;
        std::int64_t  i;
    }
    conv1, conv2;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        conv1.r              = a.simdInternal_[i];
        conv2.r              = b.simdInternal_[i];
        conv1.i              = conv1.i | conv2.i;
        res.simdInternal_[i] = conv1.r;
    }
    return res;
}

/*! \brief Bitwise xor for two SIMD4 double variables.
 *
 * Available if \ref GMX_SIMD_HAVE_LOGICAL is 1.
 *
 * \param a data1
 * \param b data2
 * \return data1 ^ data2
 */
static inline Simd4Double gmx_simdcall
operator^(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    union
    {
        double        r;
        std::int64_t  i;
    }
    conv1, conv2;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        conv1.r              = a.simdInternal_[i];
        conv2.r              = b.simdInternal_[i];
        conv1.i              = conv1.i ^ conv2.i;
        res.simdInternal_[i] = conv1.r;
    }
    return res;
}

/*! \brief Add two double SIMD4 variables.
 *
 * \param a term1
 * \param b term2
 * \return a+b
 */
static inline Simd4Double gmx_simdcall
operator+(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] + b.simdInternal_[i];
    }
    return res;
}

/*! \brief Subtract two SIMD4 variables.
 *
 * \param a term1
 * \param b term2
 * \return a-b
 */
static inline Simd4Double gmx_simdcall
operator-(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] - b.simdInternal_[i];
    }
    return res;
}

/*! \brief SIMD4 floating-point negate.
 *
 * \param a SIMD4 floating-point value
 * \return -a
 */
static inline Simd4Double gmx_simdcall
operator-(Simd4Double a)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = -a.simdInternal_[i];
    }
    return res;
}

/*! \brief Multiply two SIMD4 variables.
 *
 * \param a factor1
 * \param b factor2
 * \return a*b.
 */
static inline Simd4Double gmx_simdcall
operator*(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = a.simdInternal_[i] * b.simdInternal_[i];
    }
    return res;
}

/*! \brief SIMD4 Fused-multiply-add. Result is a*b+c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return a*b+c
 */
static inline Simd4Double gmx_simdcall
fma(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return a*b+c;
}

/*! \brief SIMD4 Fused-multiply-subtract. Result is a*b-c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return a*b-c
 */
static inline Simd4Double gmx_simdcall
fms(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return a*b-c;
}

/*! \brief SIMD4 Fused-negated-multiply-add. Result is -a*b+c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return -a*b+c
 */
static inline Simd4Double gmx_simdcall
fnma(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return c-a*b;
}

/*! \brief SIMD4 Fused-negated-multiply-subtract. Result is -a*b-c.
 *
 * \param a factor1
 * \param b factor2
 * \param c term
 * \return -a*b-c
 */
static inline Simd4Double gmx_simdcall
fnms(Simd4Double a, Simd4Double b, Simd4Double c)
{
    return -a*b-c;
}

/*! \brief SIMD4 1.0/sqrt(x) lookup.
 *
 * This is a low-level instruction that should only be called from routines
 * implementing the inverse square root in simd_math.h.
 *
 * \param x Argument, x>0
 * \return Approximation of 1/sqrt(x), accuracy is \ref GMX_SIMD_RSQRT_BITS.
 */
static inline Simd4Double gmx_simdcall
rsqrt(Simd4Double x)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        // sic - we only use single precision for the lookup
        res.simdInternal_[i] = 1.0f / std::sqrt(static_cast<float>(x.simdInternal_[i]));
    }
    return res;
};


/*! \brief SIMD4 Floating-point abs().
 *
 * \param a any floating point values
 * \return fabs(a) for each element.
 */
static inline Simd4Double gmx_simdcall
abs(Simd4Double a)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::abs(a.simdInternal_[i]);
    }
    return res;
}

/*! \brief Set each SIMD4 element to the largest from two variables.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return max(a,b) for each element.
 */
static inline Simd4Double gmx_simdcall
max(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::max(a.simdInternal_[i], b.simdInternal_[i]);
    }
    return res;
}


/*! \brief Set each SIMD4 element to the largest from two variables.
 *
 * \param a Any floating-point value
 * \param b Any floating-point value
 * \return max(a,b) for each element.
 */
static inline Simd4Double gmx_simdcall
min(Simd4Double a, Simd4Double b)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::min(a.simdInternal_[i], b.simdInternal_[i]);
    }
    return res;
}


/*! \brief SIMD4 Round to nearest integer value (in floating-point format).
 *
 * \param a Any floating-point value
 * \return The nearest integer, represented in floating-point format.
 */
static inline Simd4Double gmx_simdcall
round(Simd4Double a)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::round(a.simdInternal_[i]);
    }
    return res;
}


/*! \brief Truncate SIMD4, i.e. round towards zero - common hardware instruction.
 *
 * \param a Any floating-point value
 * \return Integer rounded towards zero, represented in floating-point format.
 *
 * \note This is truncation towards zero, not floor(). The reason for this
 * is that truncation is virtually always present as a dedicated hardware
 * instruction, but floor() frequently isn't.
 */
static inline Simd4Double gmx_simdcall
trunc(Simd4Double a)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = std::trunc(a.simdInternal_[i]);
    }
    return res;
}

/*! \brief Return dot product of two double precision SIMD4 variables.
 *
 * The dot product is calculated between the first three elements in the two
 * vectors, while the fourth is ignored. The result is returned as a scalar.
 *
 * \param a vector1
 * \param b vector2
 * \result a[0]*b[0]+a[1]*b[1]+a[2]*b[2], returned as scalar. Last element is ignored.
 */
static inline double gmx_simdcall
dotProduct(Simd4Double a, Simd4Double b)
{
    return
        (a.simdInternal_[0] * b.simdInternal_[0] +
         a.simdInternal_[1] * b.simdInternal_[1] +
         a.simdInternal_[2] * b.simdInternal_[2]);
}

/*! \brief SIMD4 double transpose
 *
 * \param[in,out] v0  Row 0 on input, column 0 on output
 * \param[in,out] v1  Row 1 on input, column 1 on output
 * \param[in,out] v2  Row 2 on input, column 2 on output
 * \param[in,out] v3  Row 3 on input, column 3 on output
 */
static inline void gmx_simdcall
transpose(Simd4Double * v0, Simd4Double * v1,
          Simd4Double * v2, Simd4Double * v3)
{
    Simd4Double t0 = *v0;
    Simd4Double t1 = *v1;
    Simd4Double t2 = *v2;
    Simd4Double t3 = *v3;
    v0->simdInternal_[0] = t0.simdInternal_[0];
    v0->simdInternal_[1] = t1.simdInternal_[0];
    v0->simdInternal_[2] = t2.simdInternal_[0];
    v0->simdInternal_[3] = t3.simdInternal_[0];
    v1->simdInternal_[0] = t0.simdInternal_[1];
    v1->simdInternal_[1] = t1.simdInternal_[1];
    v1->simdInternal_[2] = t2.simdInternal_[1];
    v1->simdInternal_[3] = t3.simdInternal_[1];
    v2->simdInternal_[0] = t0.simdInternal_[2];
    v2->simdInternal_[1] = t1.simdInternal_[2];
    v2->simdInternal_[2] = t2.simdInternal_[2];
    v2->simdInternal_[3] = t3.simdInternal_[2];
    v3->simdInternal_[0] = t0.simdInternal_[3];
    v3->simdInternal_[1] = t1.simdInternal_[3];
    v3->simdInternal_[2] = t2.simdInternal_[3];
    v3->simdInternal_[3] = t3.simdInternal_[3];
}

/*! \brief a==b for SIMD4 double
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a==b.
 */
static inline Simd4DBool gmx_simdcall
operator==(Simd4Double a, Simd4Double b)
{
    Simd4DBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] == b.simdInternal_[i]);
    }
    return res;
}

/*! \brief a!=b for SIMD4 double
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a!=b.
 */
static inline Simd4DBool gmx_simdcall
operator!=(Simd4Double a, Simd4Double b)
{
    Simd4DBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] != b.simdInternal_[i]);
    }
    return res;
}

/*! \brief a<b for SIMD4 double
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<b.
 */
static inline Simd4DBool gmx_simdcall
operator<(Simd4Double a, Simd4Double b)
{
    Simd4DBool          res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] < b.simdInternal_[i]);
    }
    return res;
}


/*! \brief a<=b for SIMD4 double.
 *
 * \param a value1
 * \param b value2
 * \return Each element of the boolean will be set to true if a<=b.
 */
static inline Simd4DBool gmx_simdcall
operator<=(Simd4Double a, Simd4Double b)
{
    Simd4DBool          res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] <= b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Logical \a and on single precision SIMD4 booleans.
 *
 * \param a logical vars 1
 * \param b logical vars 2
 * \return For each element, the result boolean is true if a \& b are true.
 *
 * \note This is not necessarily a bitwise operation - the storage format
 * of booleans is implementation-dependent.
 */
static inline Simd4DBool gmx_simdcall
operator&&(Simd4DBool a, Simd4DBool b)
{
    Simd4DBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] && b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Logical \a or on single precision SIMD4 booleans.
 *
 * \param a logical vars 1
 * \param b logical vars 2
 * \return For each element, the result boolean is true if a or b is true.
 *
 * Note that this is not necessarily a bitwise operation - the storage format
 * of booleans is implementation-dependent.
 */
static inline Simd4DBool gmx_simdcall
operator||(Simd4DBool a, Simd4DBool b)
{
    Simd4DBool         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = (a.simdInternal_[i] || b.simdInternal_[i]);
    }
    return res;
}

/*! \brief Returns non-zero if any of the boolean in SIMD4 a is True, otherwise 0.
 *
 * \param a Logical variable.
 * \return true if any element in a is true, otherwise false.
 *
 * The actual return value for truth will depend on the architecture,
 * so any non-zero value is considered truth.
 */
static inline bool gmx_simdcall
anyTrue(Simd4DBool a)
{
    bool res = false;

    for (std::size_t i = 0; i < a.simdInternal_.size(); i++)
    {
        res = res || a.simdInternal_[i];
    }
    return res;
}

/*! \brief Select from single precision SIMD4 variable where boolean is true.
 *
 * \param a Floating-point variable to select from
 * \param mask Boolean selector
 * \return  For each element, a is selected for true, 0 for false.
 */
static inline Simd4Double gmx_simdcall
selectByMask(Simd4Double a, Simd4DBool mask)
{
    Simd4Double          res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = mask.simdInternal_[i] ? a.simdInternal_[i] : 0.0;
    }
    return res;
}

/*! \brief Select from single precision SIMD4 variable where boolean is false.
 *
 * \param a Floating-point variable to select from
 * \param mask Boolean selector
 * \return  For each element, a is selected for false, 0 for true (sic).
 */
static inline Simd4Double gmx_simdcall
selectByNotMask(Simd4Double a, Simd4DBool mask)
{
    Simd4Double          res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = mask.simdInternal_[i] ? 0.0 : a.simdInternal_[i];
    }
    return res;
}


/*! \brief Vector-blend SIMD4 selection.
 *
 * \param a First source
 * \param b Second source
 * \param sel Boolean selector
 * \return For each element, select b if sel is true, a otherwise.
 */
static inline Simd4Double gmx_simdcall
blend(Simd4Double a, Simd4Double b, Simd4DBool sel)
{
    Simd4Double         res;

    for (std::size_t i = 0; i < res.simdInternal_.size(); i++)
    {
        res.simdInternal_[i] = sel.simdInternal_[i] ? b.simdInternal_[i] : a.simdInternal_[i];
    }
    return res;
}


/*! \brief Return sum of all elements in SIMD4 double variable.
 *
 * \param a SIMD4 variable to reduce/sum.
 * \return The sum of all elements in the argument variable.
 *
 */
static inline double gmx_simdcall
reduce(Simd4Double a)
{
    double sum = 0.0;

    for (std::size_t i = 0; i < a.simdInternal_.size(); i++)
    {
        sum += a.simdInternal_[i];
    }
    return sum;
}

//! \}

//! \}

//! \endcond

}      // namespace gmx

#endif // GMX_SIMD_IMPL_REFERENCE_SIMD4_DOUBLE_H
