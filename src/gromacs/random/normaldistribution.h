/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2019, by the GROMACS development team, led by
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
 * \brief The normal distribution
 *
 * Portable version of the normal distribution that generates the same sequence
 * on all platforms. Since stdlibc++ and libc++ provide different sequences
 * we prefer this one so unit tests produce the same values on all platforms.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inpublicapi
 * \ingroup module_random
 */

#ifndef GMX_RANDOM_NORMALDISTRIBUTION_H
#define GMX_RANDOM_NORMALDISTRIBUTION_H

#include <cmath>

#include <limits>

#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/classhelpers.h"

/*
 * The portable version of the normal distribution (to make sure we get the same
 * values on all platforms) has been modified from the LLVM libcxx headers,
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

/*! \brief Normal distribution
 *
 *  The C++ standard library does provide a normal distribution, but even
 *  though they all sample from the normal distribution different standard
 *  library implementations appear to return different sequences of numbers
 *  for the same random number generator. To make it easier to use GROMACS
 *  unit tests that depend on random numbers we have our own implementation.
 *
 *  Be warned that the normal distribution draws values from the random engine
 *  in a loop, so you want to make sure you use a random stream with a
 *  very large margin to make sure you do not run out of random numbers
 *  in an unlucky case (which will lead to an exception with the GROMACS
 *  default random engine).
 *
 * \tparam RealType Floating-point type, real by default in GROMACS.
 */
template<class RealType = real>
class NormalDistribution
{
public:
    /*! \brief Type of values returned */
    typedef RealType result_type;

    /*! \brief Normal distribution parameters */
    class param_type
    {
        /*! \brief Mean of normal distribution */
        result_type mean_;
        /*! \brief Standard deviation of distribution */
        result_type stddev_;

    public:
        /*! \brief Reference back to the distribution class */
        typedef NormalDistribution distribution_type;

        /*! \brief Construct parameter block
         *
         * \param mean     Mean of normal distribution
         * \param stddev   Standard deviation of normal distribution
         */
        explicit param_type(result_type mean = 0.0, result_type stddev = 1.0) :
            mean_(mean),
            stddev_(stddev)
        {
        }

        /*! \brief Return first parameter */
        result_type mean() const { return mean_; }
        /*! \brief Return second parameter */
        result_type stddev() const { return stddev_; }

        /*! \brief True if two parameter sets will return the same normal distribution.
         *
         * \param x  Instance to compare with.
         */
        bool operator==(const param_type& x) const
        {
            return mean_ == x.mean_ && stddev_ == x.stddev_;
        }

        /*! \brief True if two parameter sets will return different normal distributions
         *
         * \param x  Instance to compare with.
         */
        bool operator!=(const param_type& x) const { return !operator==(x); }
    };

public:
    /*! \brief Construct new distribution with given floating-point parameters.
     *
     * \param mean     Mean of normal distribution
     * \param stddev   Standard deviation of normal distribution
     */
    explicit NormalDistribution(result_type mean = 0.0, result_type stddev = 1.0) :
        param_(param_type(mean, stddev)),
        hot_(false),
        saved_(0)
    {
    }

    /*! \brief Construct new distribution from parameter class
     *
     * \param param  Parameter class as defined inside gmx::NormalDistribution.
     */
    explicit NormalDistribution(const param_type& param) : param_(param), hot_(false), saved_(0) {}

    /*! \brief Flush all internal saved values  */
    void reset() { hot_ = false; }

    /*! \brief Return values from normal distribution with internal parameters
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

    /*! \brief Return value from normal distribution with given parameters
     *
     *  \tparam Rng   Random engine class
     *
     *  \param  g     Random engine
     *  \param  param Parameters to use
     */
    template<class Rng>
    result_type operator()(Rng& g, const param_type& param)
    {
        result_type result;

        if (hot_)
        {
            hot_   = false;
            result = saved_;
        }
        else
        {
            UniformRealDistribution<result_type> uniformDist(-1.0, 1.0);
            result_type                          u;
            result_type                          v;
            result_type                          s;

            do
            {
                u = uniformDist(g);
                v = uniformDist(g);
                s = u * u + v * v;
            } while (s > 1.0 || s == 0.0);

            s      = std::sqrt(-2.0 * std::log(s) / s);
            saved_ = v * s;
            hot_   = true;
            result = u * s;
        }
        return result * param.stddev() + param.mean();
    }

    /*! \brief Return the mean of the normal distribution */
    result_type mean() const { return param_.mean(); }

    /*! \brief Return the standard deviation of the normal distribution */
    result_type stddev() const { return param_.stddev(); }

    /*! \brief Return the full parameter class of the normal distribution */
    param_type param() const { return param_; }

    /*! \brief Smallest value that can be returned from normal distribution */
    result_type min() const { return -std::numeric_limits<result_type>::infinity(); }

    /*! \brief Largest value that can be returned from normal distribution */
    result_type max() const { return std::numeric_limits<result_type>::infinity(); }

    /*! \brief True if two normal distributions will produce the same values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator==(const NormalDistribution& x) const
    {
        /* Equal if: Params are identical, and saved-state is identical,
         * and if we have something saved, it must be identical.
         */
        return param_ == x.param_ && hot_ == x.hot_ && (!hot_ || saved_ == x.saved_);
    }

    /*! \brief True if two normal distributions will produce different values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator!=(const NormalDistribution& x) const { return !operator==(x); }

private:
    /*! \brief Internal value for parameters, can be overridden at generation time. */
    param_type param_;
    /*! \brief True if there is a saved result to return */
    bool hot_;
    /*! \brief The saved result to return - only valid if hot_ is true */
    result_type saved_;

    GMX_DISALLOW_COPY_AND_ASSIGN(NormalDistribution);
};


} // namespace gmx

#endif // GMX_RANDOM_NORMALDISTRIBUTION_H
