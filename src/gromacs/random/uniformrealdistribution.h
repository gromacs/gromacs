/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2018,2019, by the GROMACS development team, led by
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

/*! \file
 * \brief The uniform real distribution
 *
 * Portable version of the uniform real that generates the same sequence
 * on all platforms. Since stdlibc++ and libc++ provide different sequences
 * we prefer this one so unit tests produce the same values on all platforms.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inpublicapi
 * \ingroup module_random
 */

#ifndef GMX_RANDOM_UNIFORMREALDISTRIBUTION_H
#define GMX_RANDOM_UNIFORMREALDISTRIBUTION_H

#include <cmath>

#include <algorithm>
#include <limits>

#include "gromacs/math/functions.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

/*
 * The portable version of the uniform real distribution (to make sure we get
 * the same values on all platforms) has been modified from the LLVM libcxx
 * headers, distributed under the MIT license:
 *
 * Copyright (c) The LLVM compiler infrastructure
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */

namespace gmx
{

/*! \brief Generate a floating-point value with specified number of random bits
 *
 * \tparam RealType  Floating-point type to generate
 * \tparam Bits      Number of random bits to generate
 * \tparam Rng       Random number generator class
 *
 * \param  g         Random number generator to use
 *
 * This implementation avoids the bug in libc++ and stdlibc++ (which is due
 * to the C++ standard being unclear) where 1.0 can be returned occasionally.
 *
 */
template<class RealType = real, unsigned int Bits, class Rng>
RealType generateCanonical(Rng& g)
{
    // No point in using more bits than fit in RealType
    const uint64_t digits   = std::numeric_limits<RealType>::digits;
    const uint64_t realBits = std::min(digits, static_cast<uint64_t>(Bits));
    const uint64_t range    = Rng::max() - Rng::min() + uint64_t(1);
    uint64_t       log2R    = (range == 0) ? std::numeric_limits<uint64_t>::digits : log2I(range);
    uint64_t       k        = realBits / log2R + (realBits % log2R != 0) + (realBits == 0);
    RealType       r        = Rng::max() - Rng::min() + RealType(1);
    RealType       s        = g() - Rng::min();
    RealType       base     = r;
    RealType       result;

    for (uint64_t i = 1; i < k; ++i)
    {
        s += RealType(g() - Rng::min()) * base;
        base *= r;
    }
    result = s / base;

    // This implementation is specified by the C++ standard, but unfortunately it
    // has a bug where 1.0 can be generated occasionally due to the limited
    // precision of floating point, while 0.0 is only generated half as often as
    // it should. We "solve" both these issues by swapping 1.0 for 0.0 when it happens.
    //
    // See:
    // https://llvm.org/bugs/show_bug.cgi?id=18767
    // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63176
    //
    // Note that we prefer not to use the gcc 'fix' of looping until the result
    // is smaller than 1.0, since that breaks the strict specification of the
    // number of times the rng will be called.
    //
    // This can only happen when we ask for the same number of bits that fit
    // in RealType, so by checking for that we avoid the extra code in all other
    // cases. If you are worried about it: Use RealType=double with 32 bits.
    //
    if (realBits == digits && result == 1.0)
    {
        result = 0.0;
    }
    return result;
}


/*! \brief Uniform real distribution
 *
 *  The C++ standard library does provide this distribution, but even
 *  though they all sample from the correct distribution different standard
 *  library implementations appear to return different sequences of numbers
 *  for the same random number generator. To make it easier to use GROMACS
 *  unit tests that depend on random numbers we have our own implementation.
 *
 * \tparam RealType Floating-point type, real by default in GROMACS.
 */
template<class RealType = real>
class UniformRealDistribution
{
public:
    /*! \brief Type of values returned */
    typedef RealType result_type;

    /*! \brief Uniform real distribution parameters */
    class param_type
    {
        /*! \brief Lower end of range (inclusive) */
        result_type a_;
        /*! \brief Upper end of range (exclusive) */
        result_type b_;

    public:
        /*! \brief Reference back to the distribution class */
        typedef UniformRealDistribution distribution_type;

        /*! \brief Construct parameter block
         *
         * \param a   Lower end of range (inclusive)
         * \param b   Upper end of range (exclusive)
         */
        explicit param_type(result_type a = 0.0, result_type b = 1.0) : a_(a), b_(b)
        {
            GMX_RELEASE_ASSERT(a < b, "The uniform real distribution requires a<b");
        }

        /*! \brief Return first parameter */
        result_type a() const { return a_; }
        /*! \brief Return second parameter */
        result_type b() const { return b_; }

        /*! \brief True if two parameter sets will return the same uniform real distribution.
         *
         * \param x  Instance to compare with.
         */
        bool operator==(const param_type& x) const { return a_ == x.a_ && b_ == x.b_; }

        /*! \brief True if two parameter sets will return different uniform real distributions
         *
         * \param x  Instance to compare with.
         */
        bool operator!=(const param_type& x) const { return !operator==(x); }
    };

public:
    /*! \brief Construct new distribution with given floating-point parameters.
     *
     * \param a   Lower end of range (inclusive)
     * \param b   Upper end of range (exclusive)
     */
    explicit UniformRealDistribution(result_type a = 0.0, result_type b = 1.0) :
        param_(param_type(a, b))
    {
    }

    /*! \brief Construct new distribution from parameter class
     *
     * \param param  Parameter class as defined inside gmx::UniformRealDistribution.
     */
    explicit UniformRealDistribution(const param_type& param) : param_(param) {}

    /*! \brief Flush all internal saved values  */
    void reset() {}

    /*! \brief Return values from uniform real distribution with internal parameters
     *
     * \tparam Rng  Random engine class
     *
     * \param  g    Random engine
     */
    template<class Rng>
    result_type operator()(Rng& g)
    {
        return (*this)(g, param_);
    }

    /*! \brief Return value from uniform real distribution with given parameters
     *
     * \tparam Rng   Random engine class
     *
     * \param  g     Random engine
     * \param  param Parameters to use
     */
    template<class Rng>
    result_type operator()(Rng& g, const param_type& param)
    {
        result_type r = generateCanonical<RealType, std::numeric_limits<RealType>::digits>(g);
        return (param.b() - param.a()) * r + param.a();
    }

    /*! \brief Return the lower range uniform real distribution */
    result_type a() const { return param_.a(); }

    /*! \brief Return the upper range of the uniform real distribution */
    result_type b() const { return param_.b(); }

    /*! \brief Return the full parameter class of the uniform real distribution */
    param_type param() const { return param_; }

    /*! \brief Smallest value that can be returned from uniform real distribution */
    result_type min() const { return a(); }

    /*! \brief Largest value that can be returned from uniform real distribution */
    result_type max() const { return b(); }

    /*! \brief True if two uniform real distributions will produce the same values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator==(const UniformRealDistribution& x) const { return param_ == x.param_; }

    /*! \brief True if two uniform real distributions will produce different values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator!=(const UniformRealDistribution& x) const { return !operator==(x); }

private:
    /*! \brief Internal value for parameters, can be overridden at generation time. */
    param_type param_;

    GMX_DISALLOW_COPY_AND_ASSIGN(UniformRealDistribution);
};

} // namespace gmx

#endif // GMX_RANDOM_UNIFORMREALDISTRIBUTION_H
