/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2019- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements Gaussian function evaluations on lattices and related functionality
 *
 * \author Christian Blau <blau@kth.se>
 *
 * \ingroup module_math
 */
#include "gmxpre.h"

#include "gromacs/math/gausstransform.h"

#include <cmath>

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <vector>

#include "gromacs/math/functions.h"
#include "gromacs/math/multidimarray.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/extents.h"
#include "gromacs/mdspan/layouts.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

namespace gmx
{

/********************************************************************
 * GaussianOn1DLattice::Impl
 */

class GaussianOn1DLattice::Impl
{
public:
    Impl(int numGridPointsForSpreadingHalfWidth, real sigma);
    ~Impl()                 = default;
    Impl(const Impl& other) = default;
    Impl& operator=(const Impl& other) = default;

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
    int maxEvaluatedSpreadDistance_;
    //! Width of the Gaussian function
    double sigma_;
    //! The result of the spreading calculation
    std::vector<float> spreadingResult_;
    //! Pre-calculated exp(-gridIndex^2/2 * (sigma^2)) named as in Greengard2004
    std::vector<float> e3_;
    /*! \brief Equal to std::floor(std::log(std::numeric_limits<float>::max())).
     * Above expression is not constexpr and a const variable would implicitly delete default copy
     * assignment. Therefore resorting to setting number manually.
     */
    static constexpr double c_logMaxFloat = 88.72284;
    static constexpr double c_logMinFloat = -87.33654;
};

GaussianOn1DLattice::Impl::Impl(int numGridPointsForSpreadingHalfWidth, real sigma) :
    numGridPointsForSpreadingHalfWidth_(numGridPointsForSpreadingHalfWidth),
    sigma_(sigma),
    spreadingResult_(2 * numGridPointsForSpreadingHalfWidth + 1)
{
    maxEvaluatedSpreadDistance_ =
            std::min(numGridPointsForSpreadingHalfWidth_,
                     static_cast<int>(std::floor(4 * square(sigma) * c_logMaxFloat)) - 1);
    maxEvaluatedSpreadDistance_ =
            std::min(maxEvaluatedSpreadDistance_,
                     static_cast<int>(std::floor(sigma * std::sqrt(-2.0 * c_logMinFloat))) - 1);

    std::generate_n(
            std::back_inserter(e3_), maxEvaluatedSpreadDistance_ + 1, [sigma, latticeIndex = 0]() mutable {
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
    const double e1 =
            amplitude * std::exp(-0.5 * dx * dx / square(sigma_)) / (std::sqrt(2 * M_PI) * sigma_);
    spreadingResult_[numGridPointsForSpreadingHalfWidth_] = e1;

    const double e2 = std::exp(dx / square(sigma_));

    double e2pow = e2; //< powers of e2, e2^offset

    // Move outwards from mid-point, using e2pow value for both points simultaneously
    //      o    o    o<----O---->o    o    o
    for (int offset = 1; offset < maxEvaluatedSpreadDistance_; offset++)
    {
        const double e1_3                                              = e1 * e3_[offset];
        spreadingResult_[numGridPointsForSpreadingHalfWidth_ + offset] = e1_3 * e2pow;
        spreadingResult_[numGridPointsForSpreadingHalfWidth_ - offset] = e1_3 / e2pow;
        e2pow *= e2;
    }
    // separate statement for gridpoints at the end of the range avoids
    // overflow for large sigma and saves one e2 multiplication operation
    spreadingResult_[numGridPointsForSpreadingHalfWidth_ - maxEvaluatedSpreadDistance_] =
            (e1 / e2pow) * e3_[maxEvaluatedSpreadDistance_];
    spreadingResult_[numGridPointsForSpreadingHalfWidth_ + maxEvaluatedSpreadDistance_] =
            (e1 * e2pow) * e3_[maxEvaluatedSpreadDistance_];
}

/********************************************************************
 * GaussianOn1DLattice
 */

GaussianOn1DLattice::GaussianOn1DLattice(int numGridPointsForSpreadingHalfWidth_, real sigma) :
    impl_(new Impl(numGridPointsForSpreadingHalfWidth_, sigma))
{
}

GaussianOn1DLattice::~GaussianOn1DLattice() {}

void GaussianOn1DLattice::spread(double amplitude, real dx)
{
    impl_->spread(amplitude, dx);
}

ArrayRef<const float> GaussianOn1DLattice::view()
{
    return impl_->spreadingResult_;
}

GaussianOn1DLattice::GaussianOn1DLattice(const GaussianOn1DLattice& other) :
    impl_(new Impl(*other.impl_))
{
}

GaussianOn1DLattice& GaussianOn1DLattice::operator=(const GaussianOn1DLattice& other)
{
    *impl_ = *other.impl_;
    return *this;
}

GaussianOn1DLattice::GaussianOn1DLattice(GaussianOn1DLattice&&) noexcept = default;

GaussianOn1DLattice& GaussianOn1DLattice::operator=(GaussianOn1DLattice&&) noexcept = default;

namespace
{

//! rounds real-valued coordinate to the closest integer values
IVec closestIntegerPoint(const RVec& coordinate)
{
    return { roundToInt(coordinate[XX]), roundToInt(coordinate[YY]), roundToInt(coordinate[ZZ]) };
}

/*! \brief Substracts a range from a three-dimensional integer coordinate and ensures
 * the resulting coordinate is within a lattice.
 * \param[in] index point in lattice
 * \param[in] range to be shifted
 * \returns Shifted index or zero if shifted index is smaller than zero.
 */
IVec rangeBeginWithinLattice(const IVec& index, const IVec& range)
{
    return elementWiseMax({ 0, 0, 0 }, index - range);
}

/*! \brief Adds a range from a three-dimensional integer coordinate and ensures
 * the resulting coordinate is within a lattice.
 * \param[in] index point in lattice
 * \param[in] extents extent of the lattice
 * \param[in] range to be shifted
 * \returns Shifted index or the lattice extent if shifted index is larger than the extent
 */
IVec rangeEndWithinLattice(const IVec& index, const dynamicExtents3D& extents, const IVec& range)
{
    IVec extentAsIvec(static_cast<int>(extents.extent(ZZ)),
                      static_cast<int>(extents.extent(YY)),
                      static_cast<int>(extents.extent(XX)));
    return elementWiseMin(extentAsIvec, index + range);
}


} // namespace

/********************************************************************
 * OuterProductEvaluator
 */

mdspan<const float, dynamic_extent, dynamic_extent>
OuterProductEvaluator::operator()(ArrayRef<const float> x, ArrayRef<const float> y)
{
    data_.resize(ssize(x), ssize(y));
    for (gmx::Index xIndex = 0; xIndex < ssize(x); ++xIndex)
    {
        const auto xValue = x[xIndex];
        std::transform(std::begin(y), std::end(y), begin(data_.asView()[xIndex]), [xValue](float yValue) {
            return xValue * yValue;
        });
    }
    return data_.asConstView();
}

/********************************************************************
 * IntegerBox
 */

IntegerBox::IntegerBox(const IVec& begin, const IVec& end) : begin_{ begin }, end_{ end } {}

const IVec& IntegerBox::begin() const
{
    return begin_;
}
const IVec& IntegerBox::end() const
{
    return end_;
}

bool IntegerBox::empty() const
{
    return !((begin_[XX] < end_[XX]) && (begin_[YY] < end_[YY]) && (begin_[ZZ] < end_[ZZ]));
}

IntegerBox spreadRangeWithinLattice(const IVec& center, dynamicExtents3D extent, IVec range)
{
    const IVec begin = rangeBeginWithinLattice(center, range);
    const IVec end   = rangeEndWithinLattice(center, extent, range);
    return { begin, end };
}
/********************************************************************
 * GaussianSpreadKernel
 */

IVec GaussianSpreadKernelParameters::Shape::latticeSpreadRange() const
{
    DVec range(std::ceil(sigma_[XX] * spreadWidthMultiplesOfSigma_),
               std::ceil(sigma_[YY] * spreadWidthMultiplesOfSigma_),
               std::ceil(sigma_[ZZ] * spreadWidthMultiplesOfSigma_));
    return range.toIVec();
}

/********************************************************************
 * GaussTransform3D::Impl
 */

/*! \internal \brief
 * Private implementation class for GaussTransform3D.
 */
class GaussTransform3D::Impl
{
public:
    //! Construct from extent and spreading width and range
    Impl(const dynamicExtents3D& extent, const GaussianSpreadKernelParameters::Shape& kernelShapeParameters);
    ~Impl() = default;
    //! Copy constructor
    Impl(const Impl& other) = default;
    //! Copy assignment
    Impl& operator=(const Impl& other) = default;
    //! Add another gaussian
    void add(const GaussianSpreadKernelParameters::PositionAndAmplitude& localParamters);
    //! The width of the Gaussian in lattice spacing units
    BasicVector<double> sigma_;
    //! The spread range in lattice points
    IVec spreadRange_;
    //! The result of the Gauss transform
    MultiDimArray<std::vector<float>, dynamicExtents3D> data_;
    //! The outer product of a Gaussian along the z and y dimension
    OuterProductEvaluator outerProductZY_;
    //! The three one-dimensional Gaussians, whose outer product is added to the Gauss transform
    std::array<GaussianOn1DLattice, DIM> gauss1d_;
};

GaussTransform3D::Impl::Impl(const dynamicExtents3D&                      extent,
                             const GaussianSpreadKernelParameters::Shape& kernelShapeParameters) :
    sigma_{ kernelShapeParameters.sigma_ },
    spreadRange_{ kernelShapeParameters.latticeSpreadRange() },
    data_{ extent },
    gauss1d_({ GaussianOn1DLattice(spreadRange_[XX], sigma_[XX]),
               GaussianOn1DLattice(spreadRange_[YY], sigma_[YY]),
               GaussianOn1DLattice(spreadRange_[ZZ], sigma_[ZZ]) })
{
}

void GaussTransform3D::Impl::add(const GaussianSpreadKernelParameters::PositionAndAmplitude& localParameters)
{
    const IVec closestLatticePoint = closestIntegerPoint(localParameters.coordinate_);
    const auto spreadRange =
            spreadRangeWithinLattice(closestLatticePoint, data_.asView().extents(), spreadRange_);

    // do nothing if the added Gaussian will never reach the lattice
    if (spreadRange.empty())
    {
        return;
    }

    for (int dimension = XX; dimension <= ZZ; ++dimension)
    {
        // multiply with amplitude so that Gauss3D = (amplitude * Gauss_x) * Gauss_y * Gauss_z
        const float gauss1DAmplitude = dimension > XX ? 1.0 : localParameters.amplitude_;
        gauss1d_[dimension].spread(
                gauss1DAmplitude, localParameters.coordinate_[dimension] - closestLatticePoint[dimension]);
    }

    const auto spreadZY         = outerProductZY_(gauss1d_[ZZ].view(), gauss1d_[YY].view());
    const auto spreadX          = gauss1d_[XX].view();
    const IVec spreadGridOffset = spreadRange_ - closestLatticePoint;

    // \todo optimize these loops if performance critical
    // The looping strategy uses that the last, x-dimension is contiguous in the memory layout
    for (int zLatticeIndex = spreadRange.begin()[ZZ]; zLatticeIndex < spreadRange.end()[ZZ]; ++zLatticeIndex)
    {
        const auto zSlice = data_.asView()[zLatticeIndex];

        for (int yLatticeIndex = spreadRange.begin()[YY]; yLatticeIndex < spreadRange.end()[YY]; ++yLatticeIndex)
        {
            const auto  ySlice      = zSlice[yLatticeIndex];
            const float zyPrefactor = spreadZY(zLatticeIndex + spreadGridOffset[ZZ],
                                               yLatticeIndex + spreadGridOffset[YY]);

            for (int xLatticeIndex = spreadRange.begin()[XX]; xLatticeIndex < spreadRange.end()[XX];
                 ++xLatticeIndex)
            {
                const float xPrefactor = spreadX[xLatticeIndex + spreadGridOffset[XX]];
                ySlice[xLatticeIndex] += zyPrefactor * xPrefactor;
            }
        }
    }
}

/********************************************************************
 * GaussTransform3D
 */

GaussTransform3D::GaussTransform3D(const dynamicExtents3D&                      extent,
                                   const GaussianSpreadKernelParameters::Shape& kernelShapeParameters) :
    impl_(new Impl(extent, kernelShapeParameters))
{
}

void GaussTransform3D::add(const GaussianSpreadKernelParameters::PositionAndAmplitude& localParameters)
{
    impl_->add(localParameters);
}

void GaussTransform3D::setZero()
{
    std::fill(begin(impl_->data_), end(impl_->data_), 0.);
}

basic_mdspan<float, dynamicExtents3D> GaussTransform3D::view()
{
    return impl_->data_.asView();
}

basic_mdspan<const float, dynamicExtents3D> GaussTransform3D::constView() const
{
    return impl_->data_.asConstView();
}

GaussTransform3D::~GaussTransform3D() {}

GaussTransform3D::GaussTransform3D(const GaussTransform3D& other) : impl_(new Impl(*other.impl_)) {}

GaussTransform3D& GaussTransform3D::operator=(const GaussTransform3D& other)
{
    *impl_ = *other.impl_;
    return *this;
}

GaussTransform3D::GaussTransform3D(GaussTransform3D&&) noexcept = default;

GaussTransform3D& GaussTransform3D::operator=(GaussTransform3D&&) noexcept = default;

} // namespace gmx
