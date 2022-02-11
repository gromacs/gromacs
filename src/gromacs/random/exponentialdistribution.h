/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \file
 * \brief The exponential distribution
 *
 * Portable version of the exponential distribution that generates the same
 * sequence on all platforms. Since stdlibc++ and libc++ provide different
 * sequences we prefer this one so unit tests produce the same values on all
 * platforms.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inpublicapi
 * \ingroup module_random
 */

#ifndef GMX_RANDOM_EXPONENTIALDISTRIBUTION_H
#define GMX_RANDOM_EXPONENTIALDISTRIBUTION_H

#include <cmath>

#include <limits>
#include <memory>

#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/classhelpers.h"

/*
 * The portable version of the exponential distribution (to make sure we get the
 * same values on all platforms) has been modified from the LLVM libcxx headers,
 * distributed under the MIT license:
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

/*! \brief Exponential distribution
 *
 *  The C++ standard library does provide an exponential distribution, but even
 *  though they all sample from a correct distribution, different standard
 *  library implementations appear to return different sequences of numbers
 *  for the same random number generator. To make it easier to use GROMACS
 *  unit tests that depend on random numbers we have our own implementation.
 *
 * \tparam RealType Floating-point type, real by default in GROMACS.
 */
template<class RealType = real>
class ExponentialDistribution
{
public:
    /*! \brief Type of values returned */
    typedef RealType result_type;

    /*! \brief Exponential distribution parameters */
    class param_type
    {
        /*! \brief The lambda/decay parameter */
        result_type lambda_;

    public:
        /*! \brief Reference back to the distribution class */
        typedef ExponentialDistribution distribution_type;

        /*! \brief Construct parameter block
         *
         * \param lambda   lambda/decay parameter
         */
        explicit param_type(result_type lambda = 1.0) : lambda_(lambda) {}

        /*! \brief Return lambda parameter */
        result_type lambda() const { return lambda_; }

        /*! \brief True if two parameter sets will return the same exponential distribution.
         *
         * \param x  Instance to compare with.
         */
        bool operator==(const param_type& x) const { return lambda_ == x.lambda_; }

        /*! \brief True if two parameter sets will return different exponential distributions
         *
         * \param x  Instance to compare with.
         */
        bool operator!=(const param_type& x) const { return !operator==(x); }
    };

    /*! \brief Construct new distribution with given floating-point parameter.
     *
     * \param lambda   lambda/decay parameter

     */
    explicit ExponentialDistribution(result_type lambda = 1.0) : param_(param_type(lambda)) {}

    /*! \brief Construct new distribution from parameter class
     *
     * \param param  Parameter class as defined inside gmx::ExponentialDistribution.
     */
    explicit ExponentialDistribution(const param_type& param) : param_(param) {}

    /*! \brief Flush all internal saved values  */
    void reset() {}

    /*! \brief Return values from exponential distribution with internal parameters
     *
     *  \tparam Rng   Random engine class
     *
     *  \param  g     Random engine
     */
    template<class Rng>
    result_type operator()(Rng& g)
    {
        return (*this)(g, param_);
    }

    /*! \brief Return value from exponential distribution with given parameters
     *
     *  \tparam Rng   Random engine class
     *
     *  \param  g     Random engine
     *  \param  param Parameters to use
     */
    template<class Rng>
    result_type operator()(Rng& g, const param_type& param)
    {
        return -std::log(result_type(1)
                         - generateCanonical<result_type, std::numeric_limits<result_type>::digits>(g))
               / param.lambda();
    }

    /*! \brief Return the lambda parameter of the exponential distribution */
    result_type lambda() const { return param_.lambda(); }

    /*! \brief Return the full parameter class of exponential distribution */
    param_type param() const { return param_; }

    /*! \brief Smallest value that can be returned from exponential distribution */
    result_type min() const { return 0; }

    /*! \brief Largest value that can be returned from exponential distribution */
    result_type max() const { return std::numeric_limits<result_type>::infinity(); }

    /*! \brief True if two exponential distributions will produce the same values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator==(const ExponentialDistribution& x) const { return param_ == x.param_; }

    /*! \brief True if two exponential distributions will produce different values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator!=(const ExponentialDistribution& x) const { return !operator==(x); }

private:
    /*! \brief Internal value for parameters, can be overridden at generation time. */
    param_type param_;

    GMX_DISALLOW_COPY_AND_ASSIGN(ExponentialDistribution);
};

} // namespace gmx

#endif // GMX_MATH_RANDOM_H
