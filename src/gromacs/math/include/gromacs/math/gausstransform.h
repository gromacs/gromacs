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
/*! \libinternal \file
 * \brief
 * Declares Gaussian function evaluations on lattices and related functionality
 *
 * \author Christian Blau <blau@kth.se>
 * \inlibraryapi
 * \ingroup module_math
 */
#ifndef GMX_MATH_GAUSSTRANSFORM_H
#define GMX_MATH_GAUSSTRANSFORM_H

#include <memory>
#include <vector>

#include "gromacs/math/multidimarray.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdspan/extensions.h"
#include "gromacs/mdspan/mdspan.h"
#include "gromacs/utility/real.h"

namespace gmx
{
template<typename>
class ArrayRef;
/*! \internal
 * \brief Provide result of Gaussian function evaluation on a one-dimensional lattice.
 *
 * This class owns the result of the operation and provides a view on it.
 *
 * The distance between lattice points is one. Unit length is normalized by the
 * lattice spacing, thus spreading width and range are given as multiples of
 * lattice points distances.
 *
 * This works well as approximation to piece-wise integration of a Gaussian on a
 * lattice, when \f$\sigma > 1\f$. The maximal relative error to evaluating erf
 * differences for \f$\sigma = 1\f$ is 0.01602.
 *
 */
class GaussianOn1DLattice
{
public:
    /*! \brief Construct Gaussian spreader with spreading range and Gaussian width.
     *
     * Spread weights are distributed over a non-periodic lattice of length
     * 2*numGridPointsForSpreadingHalfWidth+1. The lattice comprises a center point and
     * spreadDistance points to the left and to the right.
     *
     * \note There is a maximum spreading width
     *
     * \param[in] numGridPointsForSpreadingHalfWidth maximum distance in number of gridpoints from 0
     * \param[in] sigma Gaussian width.
     */
    GaussianOn1DLattice(int numGridPointsForSpreadingHalfWidth, real sigma);
    ~GaussianOn1DLattice();
    //! Copy constructor
    GaussianOn1DLattice(const GaussianOn1DLattice& other);
    //! Copy assignment
    GaussianOn1DLattice& operator=(const GaussianOn1DLattice& other);
    //! Move constructor
    GaussianOn1DLattice(GaussianOn1DLattice&& other) noexcept;
    //! Move assignment
    GaussianOn1DLattice& operator=(GaussianOn1DLattice&& other) noexcept;
    /*! \brief Spreads weight onto grid points in one dimension.
     *
     *
     *            .            :            |            :            .
     *            o            o            o            o            o
     *                                  O---|
     *                          latticeOffset
     * O - atom position
     * o - lattice positions
     * . : | spreading value at grid points.
     *
     * \note Highest numerical accuracy is achieved when the spreading
     *       with offset to the nearest lattice coordinated < 0.5
     *
     * Spreading on lattice coordinate \f$x_i\f$
     * \f[
     *      f(x_i) = \frac{\mathrm{amplitude}}{\sigma\sqrt{2\pi}}\exp(-\frac{\mathrm{offset}-x_i^2}{2 \sigma ^2})
     * \f]
     * \param[in] amplitude of the Gaussian spread.
     * \param[in] latticeOffset The distance to the nearest grid point in lattice coordinates.
     */
    void spread(double amplitude, real latticeOffset);
    /*! \brief Returns view on spread result. */
    ArrayRef<const float> view();

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/*! \libinternal \brief Parameters for density spreading kernels.
 */
struct GaussianSpreadKernelParameters
{
    /*! \libinternal \brief Shape parameters for Gaussian spreading kernels describe
     * the kernel shape.
     */
    struct Shape
    {
        //! The width of the Gaussian function in lattice spacings
        DVec sigma_;
        //! The range of the spreading function in multiples of sigma
        double spreadWidthMultiplesOfSigma_;
        //! The spread range in lattice coordinates
        IVec latticeSpreadRange() const;
    };
    /*! \libinternal \brief Parameters that describe the kernel position and amplitude.
     */
    struct PositionAndAmplitude
    {
        //! position of the kernel to be spread onto the lattice
        const RVec& coordinate_;
        //! amplitude of the spread kernel
        real amplitude_;
    };
};

/*! \libinternal \brief Sums Gaussian values at three dimensional lattice coordinates.
 * The Gaussian is defined as \f$A \frac{1}{\sigma^3 \sqrt(2^3\pi^3)} * \exp(-\frac{(x-x0)^2}{2
 \sigma^2})\f$ \verbatim x0:              X           x
               /   \        / \
             --     --    --   --
   lattice: |    |    |    |    |    |    |
   \endverbatim
 * The lattice has spacing 1, all coordinates are given with respect to the lattice
 * coordinates.
 */
class GaussTransform3D
{
public:
    /*! \brief Construct a three-dimensional Gauss transform.
     *
     * Transform lattice values will be zero-initialized.
     *
     * \param[in] extent of the spread lattice
     * \param[in] globalParameters of the spreading kernel
     */
    GaussTransform3D(const dynamicExtents3D&                      extent,
                     const GaussianSpreadKernelParameters::Shape& globalParameters);

    ~GaussTransform3D();

    //! Copy constructor
    GaussTransform3D(const GaussTransform3D& other);

    //! Copy assignment
    GaussTransform3D& operator=(const GaussTransform3D& other);

    //! Move constructor
    GaussTransform3D(GaussTransform3D&& other) noexcept;

    //! Move assignment
    GaussTransform3D& operator=(GaussTransform3D&& other) noexcept;

    /*! \brief Add a three dimensional Gaussian with given amplitude at a coordinate.
     * \param[in] localParameters of the spreading kernel
     */
    void add(const GaussianSpreadKernelParameters::PositionAndAmplitude& localParameters);

    //! \brief Set all values on the lattice to zero.
    void setZero();

    //! Return a view on the spread lattice.
    basic_mdspan<float, dynamicExtents3D> view();

    //! Return a const view on the spread lattice.
    basic_mdspan<const float, dynamicExtents3D> constView() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

/*! \internal \brief A 3-orthotope over integer intervals.
 */
class IntegerBox
{
public:
    //! Construct from begin and end
    IntegerBox(const IVec& begin, const IVec& end);
    //! Begin indices of the box
    const IVec& begin() const;
    //! End indices of the box
    const IVec& end() const;
    //! Empty if for any dimension, end <= begin;
    bool empty() const;

private:
    const IVec begin_; //< integer indices denoting begin of box
    const IVec end_;   //< integer indices denoting one-past end of box in any dimension
};

/*! \brief Construct a box that holds all indices that are not more than a given range remote from
 * center coordinates and still within a given lattice extent.
 *
 * \param[in] center the coordinates of the center of the spread range
 * \param[in] extent the end of the lattice, number of lattice points in each dimension
 * \param[in] range the distance from the center
 * \returns box describing the range of indices
 */
IntegerBox spreadRangeWithinLattice(const IVec& center, dynamicExtents3D extent, IVec range);

/*! \internal \brief Evaluate the outer product of two number ranges.
 * Keeps the memory for the outer product allocated.
 */
class OuterProductEvaluator
{
public:
    //! Evaluate the outer product of two float number ranges.
    mdspan<const float, dynamic_extent, dynamic_extent> operator()(ArrayRef<const float> x,
                                                                   ArrayRef<const float> y);

private:
    MultiDimArray<std::vector<float>, extents<dynamic_extent, dynamic_extent>> data_;
};

} // namespace gmx

#endif /* end of include guard: GMX_MATH_GAUSSTRANSFORM_H */
