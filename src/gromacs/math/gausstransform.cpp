/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements Gaussian function evaluations on lattices and related functionality
 *
 * \author Christian Blau <blau@kth.se>
 * \ingroup module_math
 */

#include "gmxpre.h"

#include "gausstransform.h"

#include <cmath>

#include <algorithm>

#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"

namespace gmx
{

/********************************************************************
 * GaussianOn1DLattice::Impl
 */

class GaussianOn1DLattice::Impl
{
    public:
        Impl(int numGridPointsForSpreadingHalfWidth, real sigma);
        ~Impl(){}
        Impl(const Impl &other)            = default;
        Impl &operator=(const Impl &other) = default;

        /*! \brief evaluate Gaussian function at all lattice points
         * \param[in] amplitude the amplitude of the Gaussian
         * \param[in] dx distance from the center
         */
        void spread(double amplitude, real dx);
        //! Largest distance in number of gridpoints from 0
        int numGridPointsForSpreadingHalfWidth_;
        /*! \brief Avoid overflow for E2^offset and underflow for E3(i).
         *
         * Occurs when sigma is much smaller than numGridPointsForSpreadingHalfWidth_.
         *
         * E2^offset smaller than maximum float requires
         * \f$exp(dx / (2*square(sigma))^numGridPointsForSpreadingHalfWidth_ \leq max_float \f$
         * The maximum expected distance of the Gaussian center to the next lattice point is dx = 0.5,
         * thus the maximum spread distance here is \f$4 * sigma^2 * \log(\mathrm{maxfloat})\f$ .
         *
         * E3(i) larger than minmium float requires
         * exp(i^2 / 2*(sigma)^2) > min_float
         * Thus the maximum spread distance here is \f$\sigma \sqrt(-2\log(\mathrm{minfloat}))\f$
         */
        int                    maxEvaluatedSpreadDistance_;
        //! Width of the Gaussian function
        double                 sigma_;
        //! The result of the spreading calculation
        std::vector<float>     spreadingResult_;
        //! Pre-calculated exp(-gridIndex^2/2 * (sigma^2)) named as in Greengard2004
        std::vector<float>     e3_;
        /*! \brief Equal to std::floor(std::log(std::numeric_limits<float>::max())).
         * Above expression is not constexpr and a const variable would implicitly delete default copy assignment.
         * Therefore resorting to setting number manually.
         */
        static constexpr double c_logMaxFloat =  88.72284;
        static constexpr double c_logMinFloat = -87.33654;
};

GaussianOn1DLattice::Impl::Impl(int numGridPointsForSpreadingHalfWidth, real sigma) :
    numGridPointsForSpreadingHalfWidth_(numGridPointsForSpreadingHalfWidth),
    sigma_(sigma),
    spreadingResult_(2 * numGridPointsForSpreadingHalfWidth + 1)
{
    maxEvaluatedSpreadDistance_ = std::min(numGridPointsForSpreadingHalfWidth_, static_cast<int>(std::floor(4 * square(sigma) * c_logMaxFloat)) - 1);
    maxEvaluatedSpreadDistance_ = std::min(maxEvaluatedSpreadDistance_, static_cast<int>(std::floor(sigma * sqrt(-2.0 * c_logMinFloat))) - 1);

    std::generate_n(std::back_inserter(e3_), maxEvaluatedSpreadDistance_ + 1,
                    [sigma, latticeIndex = 0]() mutable {
                        return std::exp(-0.5 * square(latticeIndex++ / sigma));
                    });

    std::fill(std::begin(spreadingResult_), std::end(spreadingResult_), 0.);
};

void GaussianOn1DLattice::Impl::spread(double amplitude, real dx)
{
    /* The spreading routine implements the fast gaussian gridding as in
     *
     * Leslie Greengard and June-Yub Lee,
     * "Accelerating the Nonuniform Fast Fourier Transform"
     * SIAM REV 2004 Vol. 46, No. 3, pp. 443-454 DOI. 10.1137/S003614450343200X
     *
     * Following the naming conventions for e1, e2 and e3, nu = 1, m = numGridPointsForSpreadingHalfWidth_.
     *
     * Speed up is achieved by factorization of the exponential that is evaluted
     * at regular lattice points i, where the distance from the
     * Gaussian center is \f$x-i\f$:
     *
     * \f[
     *      a * \exp(-(x^2-2*i*x+ i^2)/(2*\sigma^2)) =
     *      a * \exp(-x^2/2*\sigma^2) * \exp(x/\sigma^2)^i * \exp(i/2*sigma^2) =
     *      e_1(x) * e_2(x)^i * e_3(i)
     * \f]
     *
     * Requiring only two exp evaluations per spreading operation.
     *
     */
    const double e1 = amplitude * exp(-0.5 * dx * dx / square(sigma_)) / (sqrt(2 * M_PI) * sigma_);
    spreadingResult_[numGridPointsForSpreadingHalfWidth_] = e1;

    const double e2 = exp(dx / square(sigma_));

    double       e2pow = e2; //< powers of e2, e2^offset

    // Move outwards from mid-point, using e2pow value for both points simultaneously
    //      o    o    o<----O---->o    o    o
    for (int offset = 1; offset < maxEvaluatedSpreadDistance_; offset++)
    {
        const double e1_3 = e1 * e3_[offset];
        spreadingResult_[numGridPointsForSpreadingHalfWidth_ + offset] = e1_3 * e2pow;
        spreadingResult_[numGridPointsForSpreadingHalfWidth_ - offset] = e1_3 / e2pow;
        e2pow *= e2;
    }
    // seperate statement for gridpoints at the end of the range avoids
    // overflow for large sigma and saves one e2 multiplication operation
    spreadingResult_[numGridPointsForSpreadingHalfWidth_ - maxEvaluatedSpreadDistance_] = (e1 / e2pow) * e3_[maxEvaluatedSpreadDistance_];
    spreadingResult_[numGridPointsForSpreadingHalfWidth_ + maxEvaluatedSpreadDistance_] = (e1 * e2pow) * e3_[maxEvaluatedSpreadDistance_];
}

/********************************************************************
 * GaussianOn1DLattice
 */

GaussianOn1DLattice::GaussianOn1DLattice(int numGridPointsForSpreadingHalfWidth_, real sigma) : impl_(new Impl(numGridPointsForSpreadingHalfWidth_, sigma))
{
}

GaussianOn1DLattice::~GaussianOn1DLattice () {}

void GaussianOn1DLattice::spread(double amplitude, real dx)
{
    impl_->spread(amplitude, dx);
}

ArrayRef<const float> GaussianOn1DLattice::view()
{
    return impl_->spreadingResult_;
}

GaussianOn1DLattice::GaussianOn1DLattice(const GaussianOn1DLattice &other)
    : impl_(new Impl(*other.impl_))
{
}

GaussianOn1DLattice &GaussianOn1DLattice::operator=(const GaussianOn1DLattice &other)
{
    *impl_ = *other.impl_;
    return *this;
}

GaussianOn1DLattice::GaussianOn1DLattice(GaussianOn1DLattice &&) noexcept = default;

GaussianOn1DLattice &GaussianOn1DLattice::operator=(GaussianOn1DLattice &&) noexcept = default;

}    // namespace gmx
