/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2026- The GROMACS Authors
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
 * \brief Device-side tabulated normal distribution
 *
 * Mirrors gmx::TabulatedNormalDistribution but takes a device-side pointer to
 * the distribution table instead of relying on a class-static one. The
 * annotations are written with the GMX_FUNC_ATTRIBUTE macros so the same
 * header can be used from CUDA, HIP, and SYCL device code.
 *
 * \author Erik Lindahl <erik.lindahl@gmail.com>
 * \author Magnus Lundborg <lundborg.magnus@gmail.com>
 * \inlibraryapi
 * \ingroup module_random
 */

#ifndef GMX_RANDOM_TABULATEDNORMALDISTRIBUTION_GPU_H
#define GMX_RANDOM_TABULATEDNORMALDISTRIBUTION_GPU_H

#include <cstdint>

#include <limits>

#include "gromacs/gpu_utils/gpu_device_macros.h"
#include "gromacs/utility/classhelpers.h"

#if !defined DOXYGEN || !DOXYGEN

namespace gmx
{

namespace detail
{

//! Number of bits that determines the resolution of the lookup table for the normal distribution.
constexpr int c_TabulatedNormalDistributionGpuDefaultBits = 14;

} // namespace detail

/*! \brief Tabulated normal random distribution for use on the GPU.
 *
 * Functionally equivalent to gmx::TabulatedNormalDistribution, but the
 * distribution table is provided as a device pointer rather than as a
 * class-static array. The values returned for a given sequence of random
 * bits are identical to those of the host-side class when the same table
 * data is used.
 *
 *  \tparam tableBits Size of the table, specified in bits. The storage
 *                    space required is sizeof(float)*2^tableBits. To
 *                    keep things sane this is limited to 24 bits.
 */
template<unsigned int tableBits = detail::c_TabulatedNormalDistributionGpuDefaultBits>
class TabulatedNormalDistributionGpu
{
    static_assert(tableBits <= 24,
                  "Normal distribution table is limited to 24bits (64MB in single precision)");

public:
    /*! \brief  Type of normal distribution results */
    typedef float result_type;

    /*! \brief  Normal distribution parameter class (mean and stddev) */
    class param_type
    {
    public:
        /*! \brief The type of distribution the parameters describe */
        typedef TabulatedNormalDistributionGpu distribution_type;

        /*! \brief Constructor. Default is classical distr. with mean 0, stddev 1.
         *
         * \param mean       Expectation value.
         * \param stddev     Standard deviation.
         */
        GMX_FUNC_ATTRIBUTE param_type(result_type mean = 0.0, result_type stddev = 1.0) :
            mean_(mean), stddev_(stddev)
        {
        }

        /*! \brief Return mean parameter of normal distribution */
        GMX_FUNC_ATTRIBUTE result_type mean() const { return mean_; }

        /*! \brief Return standard deviation parameter of normal distribution */
        GMX_FUNC_ATTRIBUTE result_type stddev() const { return stddev_; }

    private:
        result_type mean_;
        result_type stddev_;
    };

    /*! \brief Construct new normal distribution with specified mean & stddev.
     *
     *  \param d_tablePtr Pointer to distribution table on the device.
     *  \param mean       Mean value of tabulated normal distribution.
     *  \param stddev     Standard deviation of tabulated normal distribution.
     */
    GMX_FUNC_ATTRIBUTE explicit TabulatedNormalDistributionGpu(const float* d_tablePtr,
                                                               result_type  mean   = 0.0,
                                                               result_type  stddev = 1.0) :
        param_(param_type(mean, stddev)), table_(d_tablePtr), savedRandomBits_(0), savedRandomBitsLeft_(0)
    {
    }

    /*! \brief Clear all internal saved random bits from the random engine */
    GMX_FUNC_ATTRIBUTE void reset() { savedRandomBitsLeft_ = 0; }

    /*! \brief Return normal distribution value specified by internal parameters.
     *
     * \tparam Rng   Random engine type used to provide uniform random bits.
     * \param  g     Random engine of class Rng. For normal GROMACS usage
     *               you likely want to use ThreeFry2x64.
     */
    template<class Rng>
    GMX_FUNC_ATTRIBUTE result_type operator()(Rng& g)
    {
        return (*this)(g, param_);
    }

    /*! \brief Return normal distribution value specified by given parameters
     *
     * \tparam Rng   Random engine type used to provide uniform random bits.
     * \param  g     Random engine of class Rng.
     * \param  param Parameters used to specify normal distribution.
     */
    template<class Rng>
    GMX_FUNC_ATTRIBUTE result_type operator()(Rng& g, const param_type& param)
    {
        if (savedRandomBitsLeft_ < tableBits)
        {
            savedRandomBits_     = static_cast<uint64_t>(g());
            savedRandomBitsLeft_ = std::numeric_limits<typename Rng::result_type>::digits;
        }
        result_type value = table_[savedRandomBits_ & ((1ULL << tableBits) - 1)];
        savedRandomBits_ >>= tableBits;
        savedRandomBitsLeft_ -= tableBits;
        return param.mean() + value * param.stddev();
    }

private:
    /*! \brief Parameters of normal distribution (mean and stddev) */
    param_type param_;
    /*! \brief Pointer to tabulated values of normal distribution (on the device) */
    const float* table_;
    /*! \brief Saved output from random engine, shifted tableBits right each time */
    uint64_t savedRandomBits_;
    /*! \brief Number of valid bits remaining in savedRandomBits_ */
    unsigned int savedRandomBitsLeft_;

    GMX_DISALLOW_COPY_AND_ASSIGN(TabulatedNormalDistributionGpu);
};

} // namespace gmx

#endif

#endif // GMX_RANDOM_TABULATEDNORMALDISTRIBUTION_GPU_H
