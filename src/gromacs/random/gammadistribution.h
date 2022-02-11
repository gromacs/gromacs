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
 * \brief The gamma distribution
 *
 * Portable version of the gamma distribution that generates the same sequence
 * on all platforms.
 *
 * \note The gamma distribution is broken in some standard library headers
 * (including those shipped with gcc-4.9), and it is not guaranteed to
 * generate the same result on stdlibc++ and libc++. Use this one instead so
 * our unit tests produce the same values on all platforms.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inpublicapi
 * \ingroup module_random
 */

#ifndef GMX_RANDOM_GAMMADISTRIBUTION_H
#define GMX_RANDOM_GAMMADISTRIBUTION_H

#include <cmath>

#include <limits>
#include <memory>

#include "gromacs/math/functions.h"
#include "gromacs/random/exponentialdistribution.h"
#include "gromacs/random/normaldistribution.h"
#include "gromacs/random/uniformrealdistribution.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/exceptions.h"

/*
 * The workaround implementation for the broken std::gamma_distribution in the
 * gcc-4.6 headers has been modified from the LLVM libcxx headers, distributed
 * under the MIT license:
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

/*! \brief Gamma distribution
 *
 *  The C++ standard library does provide a gamma distribution, but when
 *  using libstdc++-4.4.7 with at least gcc-4.6 the headers
 *  produce errors. Even for newer compilers, libstdc++ and libc++ appear to
 *  use different algorithms to generate it, which means their values differ
 *  in contrast to the uniform and normal distributions where they are
 *  identical. To avoid both compiler bugs and make it easier to use
 *  GROMACS unit tests that depend on random numbers, we have our
 *  own implementation.
 *
 *  Be warned that the gamma distribution works like the standard
 *  normal distribution and keeps drawing values from the random engine
 *  in a loop, so you want to make sure you use a random stream with a
 *  very large margin to make sure you do not run out of random numbers
 *  in an unlucky case (which will lead to an exception with the GROMACS
 *  default random engine).
 *
 *  The gamma distribution is defined as
 *
 * \f[
 *     p(x|\alpha,\theta) = \frac{1}{\Gamma(\alpha)\theta^{alpha}} x^{\alpha - 1}
 * e^{-\frac{x}{\theta}}, x\geq 0 , \alpha>0, \theta>0
 * \f]
 *
 * In this definition, the parameter &alpha; is the so-called shape, while &theta; is
 * so-called scale. This distribution will have the expectation value &alpha;&theta;, and
 *  the variance &alpha;&theta;&theta;.
 *
 * \note The gamma distribution is sometimes defined in terms of a parameter
 *  &beta; that is the inverse of &theta; (i.e., a rate rather than scale),
 *  while in other cases the parameter in the
 *  definition above is simply called &beta;. We can't do a lot about these different
 *  definitions, so make sure you look at the expectation value and variance unless
 *  you are certain you are using e.g. scale or rate for the second parameter - do not
 *  trust that it e.g. rate just because it is called &beta; in your equations.
 *
 *  \note For now, we generate the gamma distribution using the algorithm from
 *  Marsaglia G, Tsang WW (2000). ACM Trans. Math. Softw. 26(3), 363-372. DOI:10.1145/358407.358414
 *
 * \tparam RealType Floating-point type, real by default in GROMACS.
 */
template<class RealType = real>
class GammaDistribution
{
public:
    /*! \brief Type of values returned */
    typedef RealType result_type;

    /*! \brief Gamma distribution parameters */
    class param_type
    {
        /*! \brief Shape parameter of gamma distribution */
        result_type alpha_;
        /*! \brief Scale parameter of gamma distribution */
        result_type theta_;

    public:
        /*! \brief Reference back to the distribution class */
        typedef GammaDistribution distribution_type;

        /*! \brief Construct parameter block
         *
         * \param alpha  Shape parameter of gamma distribution
         * \param theta  Scale parameter of gamma distribution
         *
         *  \throws InvalidInputError if either parameter is negative or zero.
         */
        explicit param_type(result_type alpha = 1.0, result_type theta = 1.0) :
            alpha_(alpha), theta_(theta)
        {
            if (alpha <= 0 || theta <= 0)
            {
                GMX_THROW(
                        InvalidInputError("Both parameters in the gamma distribution must be >0."));
            }
        }

        /*! \brief Return shape parameter */
        result_type alpha() const { return alpha_; }
        /*! \brief Return scale parameter */
        result_type theta() const { return theta_; }

        /*! \brief True if two parameter sets will return the same distribution.
         *
         * \param x  Instance to compare with.
         */
        bool operator==(const param_type& x) const
        {
            return alpha_ == x.alpha_ && theta_ == x.theta_;
        }

        /*! \brief True if two parameter sets will return different distributions
         *
         * \param x  Instance to compare with.
         */
        bool operator!=(const param_type& x) const { return !operator==(x); }
    };

    /*! \brief Construct new distribution with given floating-point parameters.
     *
     * \param alpha  Shape parameter of gamma distribution
     * \param theta  Scale parameter of gamma distribution
     *
     *  \throws InvalidInputError if either parameter is negative or zero.
     */
    explicit GammaDistribution(result_type alpha = 1.0, result_type theta = 1.0) :
        param_(param_type(alpha, theta))
    {
    }

    /*! \brief Construct new distribution from parameter class
     *
     * \param param  Parameter class as defined inside gmx::GammaDistribution.
     */
    explicit GammaDistribution(const param_type& param) : param_(param) {}

    /*! \brief Flush all internal saved values  */
    void reset() {}

    /*! \brief Return values from gamma distribution with internal parameters
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

    /*! \brief Return value from gamma distribution with given parameters
     *
     *  \tparam Rng   Random engine class
     *
     *  \param  g     Random engine
     *  \param  param Parameters to use
     */
    template<class Rng>
    result_type operator()(Rng& g, const param_type& param)
    {
        result_type                          alpha = param.alpha();
        UniformRealDistribution<result_type> uniformDist(std::numeric_limits<result_type>::min(), 1);

        if (alpha == result_type(1))
        {
            // Special case; when alpha is unity, it is a plain exponential distribution,
            // which we can calculate faster by just using an exponential distribution.
            ExponentialDistribution<result_type> expDist;
            return expDist(g) * param.theta();
        }
        else if (alpha > result_type(1))
        {
            NormalDistribution<result_type> normDist;

            for (;;)
            {
                result_type d = alpha - result_type(1) / result_type(3);
                result_type c = result_type(1) / result_type(3) * invsqrt(d);
                result_type x, v;

                do
                {
                    x = normDist(g);
                    v = result_type(1) + c * x;
                } while (v <= result_type(0));
                v             = v * v * v;
                result_type u = uniformDist(g);
                result_type y = x * x;

                // Sieve; we first check a computationally cheaper expression that catches the majority of cases
                if (u < result_type(1) - result_type(0.0331) * y * y)
                {
                    return (d * v * param.theta());
                }
                // If we got here, we need to evaluate the two log functions to compare with the exact expression
                if (std::log(u) < result_type(0.5) * y + d * (result_type(1) - v + std::log(v)))
                {
                    return (d * v * param.theta());
                }
            }
        }
        else // alpha < 1
        {
            result_type x = this->operator()(g, param_type(alpha + result_type(1), param.theta()));
            return x * std::pow(uniformDist(g), result_type(1) / alpha);
        }
    }

    /*! \brief Return the shape parameter of gamma distribution */
    result_type alpha() const { return param_.alpha(); }

    /*! \brief Return the scale parameter of gamma distribution */
    result_type theta() const { return param_.theta(); }

    /*! \brief Return the full parameter class of gamma distribution */
    param_type param() const { return param_; }

    /*! \brief Smallest value that can be returned from gamma distribution */
    result_type min() const { return 0; }

    /*! \brief Largest value that can be returned from gamma distribution */
    result_type max() const { return std::numeric_limits<result_type>::infinity(); }

    /*! \brief True if two gamma distributions will produce the same values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator==(const GammaDistribution& x) const { return param_ == x.param_; }

    /*! \brief True if two gamma distributions will produce different values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator!=(const GammaDistribution& x) const { return !operator==(x); }

private:
    /*! \brief Internal value for parameters, can be overridden at generation time. */
    param_type param_;

    GMX_DISALLOW_COPY_AND_ASSIGN(GammaDistribution);
};

} // namespace gmx

#endif // GMX_RANDOM_GAMMADISTRIBUTION_H
