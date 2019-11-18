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
 * \brief The uniform integer distribution
 *
 * Portable version of the uniform integer that generates the same sequence
 * on all platforms. Since stdlibc++ and libc++ provide different sequences
 * we prefer this one so unit tests produce the same values on all platforms.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \inpublicapi
 * \ingroup module_random
 */

#ifndef GMX_RANDOM_UNIFORMINTDISTRIBUTION_H
#define GMX_RANDOM_UNIFORMINTDISTRIBUTION_H

#include <limits>

#include "gromacs/math/functions.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/*! \brief Uniform integer distribution
 *
 *  The C++ standard library does provide this distribution, but even
 *  though they all sample from the correct distribution different standard
 *  library implementations appear to return different sequences of numbers
 *  for the same random number generator. To make it easier to use GROMACS
 *  unit tests that depend on random numbers we have our own implementation.
 *
 * \tparam IntType Integer type, int by default.
 */
template<class IntType = int>
class UniformIntDistribution
{
public:
    /*! \brief Type of values returned */
    typedef IntType result_type;

    /*! \brief Uniform int distribution parameters */
    class param_type
    {
        /*! \brief Lower end of range (inclusive) */
        result_type a_;
        /*! \brief Upper end of range (inclusive) */
        result_type b_;

    public:
        /*! \brief Reference back to the distribution class */
        typedef UniformIntDistribution distribution_type;

        /*! \brief Construct parameter block
         *
         * \param a   Lower end of range (inclusive)
         * \param b   Upper end of range (inclusive)
         */
        explicit param_type(result_type a = 0, result_type b = std::numeric_limits<result_type>::max()) :
            a_(a),
            b_(b)
        {
            GMX_RELEASE_ASSERT(a <= b, "The uniform integer distribution requires a<=b");
        }

        /*! \brief Return lower range */
        result_type a() const { return a_; }
        /*! \brief Return upper range */
        result_type b() const { return b_; }

        /*! \brief True if two parameter sets will return the same uniform int distribution.
         *
         * \param x  Instance to compare with.
         */
        bool operator==(const param_type& x) const
        {
            // rangeBits is a function of a & b, so it does not have to be tested
            return a_ == x.a_ && b_ == x.b_;
        }

        /*! \brief True if two parameter sets will return different uniform int distributions
         *
         * \param x  Instance to compare with.
         */
        bool operator!=(const param_type& x) const { return !operator==(x); }
    };

public:
    /*! \brief Construct new distribution with given integer parameters.
     *
     * \param a   Lower end of range (inclusive)
     * \param b   Upper end of range (inclusive)
     */
    explicit UniformIntDistribution(result_type a = 0,
                                    result_type b = std::numeric_limits<result_type>::max()) :
        param_(param_type(a, b)),
        savedRandomBits_(0),
        savedRandomBitsLeft_(0)
    {
    }

    /*! \brief Construct new distribution from parameter class
     *
     * \param param  Parameter class as defined inside gmx::UniformIntDistribution.
     */
    explicit UniformIntDistribution(const param_type& param) :
        param_(param),
        savedRandomBits_(0),
        savedRandomBitsLeft_(0)
    {
    }

    /*! \brief Flush all internal saved values  */
    void reset() { savedRandomBitsLeft_ = 0; }

    /*! \brief Return values from uniform int distribution with internal parameters
     *
     * \tparam Rng  Uniform random engine class
     *
     * \param  g    Random engine
     */
    template<class Rng>
    result_type operator()(Rng& g)
    {
        return (*this)(g, param_);
    }

    /*! \brief Return value from uniform int distribution with given parameters
     *
     * \tparam Rng   Uniform random engine class
     *
     * \param  g     Random engine
     * \param  param Parameters to use
     */
    template<class Rng>
    result_type operator()(Rng& g, const param_type& param)
    {
        static_assert(sizeof(typename Rng::result_type) >= sizeof(uint32_t),
                      "The random engine result_type should be 32 or 64 bits");

        result_type  range = param.b() - param.a();
        unsigned int rangeBits;
        result_type  result;

        if (range == 0)
        {
            return param.a();
        }
        else if (range == std::numeric_limits<result_type>::max())
        {
            rangeBits = std::numeric_limits<result_type>::digits; // Use all bits in type
        }
        else
        {
            if (sizeof(result_type) == sizeof(uint32_t))
            {
                rangeBits = log2I(static_cast<uint32_t>(range));
            }
            else
            {
                rangeBits = log2I(range);
            }
            rangeBits += ((range >> rangeBits) > 0);
        }

        do
        {
            if (savedRandomBitsLeft_ < rangeBits)
            {
                savedRandomBits_     = static_cast<uint64_t>(g());
                savedRandomBitsLeft_ = std::numeric_limits<typename Rng::result_type>::digits;

                if (sizeof(typename Rng::result_type) == sizeof(uint32_t))
                {
                    savedRandomBits_ <<= std::numeric_limits<uint32_t>::digits;
                    savedRandomBits_ |= g();
                    savedRandomBitsLeft_ += std::numeric_limits<uint32_t>::digits;
                }
            }
            result = savedRandomBits_;
            savedRandomBits_ >>= rangeBits;
            result = result - (savedRandomBits_ << rangeBits);
            savedRandomBitsLeft_ -= rangeBits;
        } while (result > range);

        return result + param.a();
    }

    /*! \brief Return the lower range uniform int distribution */
    result_type a() const { return param_.a(); }

    /*! \brief Return the upper range of the uniform int distribution */
    result_type b() const { return param_.b(); }

    /*! \brief Return the full parameter class of the uniform int distribution */
    param_type param() const { return param_; }

    /*! \brief Smallest value that can be returned from uniform int distribution */
    result_type min() const { return a(); }

    /*! \brief Largest value that can be returned from uniform int distribution */
    result_type max() const { return b(); }

    /*! \brief True if two uniform int distributions will produce the same values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator==(const UniformIntDistribution& x) const { return param_ == x.param_; }

    /*! \brief True if two uniform int distributions will produce different values.
     *
     * \param  x     Instance to compare with.
     */
    bool operator!=(const UniformIntDistribution& x) const { return !operator==(x); }

private:
    /*! \brief Internal value for parameters, can be overridden at generation time. */
    param_type param_;
    /*! \brief Saved output from random engine, shifted tableBits right each time */
    uint64_t savedRandomBits_;
    /*! \brief Number of valid bits remaining i savedRandomBits_ */
    unsigned int savedRandomBitsLeft_;

    GMX_DISALLOW_COPY_AND_ASSIGN(UniformIntDistribution);
};

} // namespace gmx

#endif // GMX_RANDOM_UNIFORMINTDISTRIBUTION_H
