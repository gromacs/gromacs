/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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
/*!  \file
 * \brief
 * Defines volume data containers.
 *
 * \author Christian Blau <cblau@gwdg.de>
 * \inpublicapi
 */

#include "gmxpre.h"

#include <cmath>

#include "gausstransform.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/griddata/canonicalvectorbasis.h"
#include "gromacs/utility/compare.h"

namespace gmx
{

void
GaussTransform::setSigma(real sigma)
{
    sigma_ = sigma;
};

void
GaussTransform::setSpreadWidthInNSigma(real nSigma)
{
    nSigma_ = nSigma;
};

IVec
GaussTransform::minimumUsedGridIndex()
{
    return minimumUsedGridIndex_;
};

IVec
GaussTransform::maximumUsedGridIndex()
{
    return maximumUsedGridIndex_;
};

/*! \brief
 * True, if length is same in x,y,z -direction.
 */
bool GaussTransform::allVectorsSameLength(const CanonicalVectorBasis<DIM> &basis, real ftol, real abstol) const
{
    auto firstLength = basis.basisVectorLength(0);
    return equal_real(basis.basisVectorLength(YY), firstLength, ftol, abstol) && equal_real(basis.basisVectorLength(ZZ), firstLength, ftol, abstol);
}

void
GaussTransform::setGrid(GridDataReal3D * grid, bool gaussianIsIntegrated)
{
    grid_                 = grid;
    minimumUsedGridIndex_ = {GMX_INT32_MAX, GMX_INT32_MAX, GMX_INT32_MAX};
    maximumUsedGridIndex_ = {0, 0, 0};
    if (!allVectorsSameLength(grid_->getGrid().unitCell(), 1e-6, 1))
    {
        GMX_THROW(gmx::InconsistentInputError("Grid needs to be evently spaced to use the current implementation of fast gaussian gridding."));
    }
    nu_       = grid_->getGrid().unitCell().basisVectorLength(XX)/sigma_;
    m_spread_ = int(ceil(nSigma_/nu_)); // number of grid cells for spreading

    spread_2d_.resize(2*m_spread_+1);
    for (int l_z = 0; l_z < 2 * m_spread_ + 1; ++l_z)
    {
        spread_2d_[l_z].resize(2*m_spread_+1);
    }

    ceilSqrtLUT_.resize(m_spread_+1);
    for (int zDistance = 0; zDistance <= m_spread_; zDistance++)
    {
        ceilSqrtLUT_[zDistance].resize( (int)std::ceil(std::sqrt(m_spread_*m_spread_ - zDistance*zDistance))+1);
        for (int yDistance = 0; yDistance*yDistance <= m_spread_*m_spread_ - zDistance*zDistance; yDistance++)
        {
            ceilSqrtLUT_[zDistance][yDistance] = (int)std::ceil(std::sqrt(m_spread_*m_spread_ - zDistance*zDistance- yDistance*yDistance));
        }
    }
    for (auto &spreader : OneDSpreader_)
    {
        if (gaussianIsIntegrated)
        {
            spreader = std::unique_ptr<ISpreader>(new GaussIntergralOneD(m_spread_, nu_,  sigma_));
        }
        else
        {
            spreader = std::unique_ptr<ISpreader>(new GaussTransformOneD(m_spread_, nu_,  sigma_));
        }
    }
};

void
GaussTransform::tensor_product_2d_()
{
    for (int l_z = 0; l_z < 2 * m_spread_ + 1; ++l_z)
    {
        const real spread_z = spread_1d_[ZZ][l_z];
        for (int l_y = 0; l_y < 2 * m_spread_ + 1; ++l_y)
        {
            spread_2d_[l_z][l_y] = spread_z*spread_1d_[YY][l_y];
        }
    }
}


void
GaussTransform::tensor_product_()
{
    IVec minIndex;
    IVec maxIndex;
    for (size_t i = XX; i <= ZZ; ++i)
    {
        minIndex[i] = std::max(0, atomGridIndex_[i]-m_spread_);
        maxIndex[i] = std::min(atomGridIndex_[i]+m_spread_, grid_->getGrid().lattice().extend()[i]-1);

        // return if spread density never reaches grid
        if (maxIndex[i] < 0 || minIndex[i] > grid_->getGrid().lattice().extend()[i] - 1)
        {
            return;
        }

        minimumUsedGridIndex_[i] = std::min(minIndex[i], minimumUsedGridIndex_[i]);
        maximumUsedGridIndex_[i] = std::max(maxIndex[i], maximumUsedGridIndex_[i]);
    }


    for (int zzIndex = minIndex[ZZ]; zzIndex <= maxIndex[ZZ]; ++zzIndex)
    {
        int d_z              = zzIndex - atomGridIndex_[ZZ];
        int localGridIndexZZ = d_z + m_spread_;

        // spread spheres instead of cubes
        int yyIndexStart = std::max(minIndex[YY], atomGridIndex_[YY] - ceilSqrtLUT_[std::abs(d_z)][0]);
        int yyIndexEnd   = std::min(maxIndex[YY], atomGridIndex_[YY] + ceilSqrtLUT_[std::abs(d_z)][0]);

        for (int yyIndex = yyIndexStart; yyIndex <= yyIndexEnd; ++yyIndex)
        {
            int        d_y              = yyIndex - atomGridIndex_[YY];
            int        localGridIndexYY = d_y + m_spread_;
            const real spread_zy        = spread_2d_[localGridIndexZZ][localGridIndexYY];

            int        xxIndexStart = std::max(minIndex[XX], atomGridIndex_[XX] - ceilSqrtLUT_[std::abs(d_z)][std::abs(d_y)]);
            xxIndexStart = std::min(maxIndex[XX], xxIndexStart);

            int  xxIndexEnd   = std::min(maxIndex[XX], atomGridIndex_[XX] + ceilSqrtLUT_[std::abs(d_z)][std::abs(d_y)]);
            xxIndexEnd = std::max(minIndex[XX], xxIndexEnd);

            const int  xxIndexLocalStart      = xxIndexStart - atomGridIndex_[XX]+m_spread_;
            const int  numberSpreadVoxelsXX   = xxIndexEnd-xxIndexStart;
            auto       voxel                  = grid_->iteratorAtMultiIndex({{xxIndexStart, yyIndex, zzIndex}});
            auto       spread_1d_XX           = &(spread_1d_[XX][xxIndexLocalStart]);

            for (int l_x = 0; l_x <= numberSpreadVoxelsXX; ++l_x)
            {
                *voxel += spread_zy * (*spread_1d_XX);
                ++spread_1d_XX;
                ++voxel;
            }
        }
    }
};

void GaussTransform::prepare_2d_grid(const rvec x, const real weight)
{
    atomGridIndex_ = grid_->getGrid().coordinateToFloorMultiIndex({{ x[XX], x[YY], x[ZZ] }});
    RVec dx; // (x-nearest voxel)/sigma
    for (size_t i = XX; i <= ZZ; ++i)
    {
        dx[i]  = (x[i]-grid_->getGrid().multiIndexToCoordinate(atomGridIndex_)[i])/sigma_;
    }
    spread_1d_[XX] = OneDSpreader_[XX]->spread(weight, dx[XX]);
    spread_1d_[YY] = OneDSpreader_[YY]->spread(1, dx[YY]);
    spread_1d_[ZZ] = OneDSpreader_[ZZ]->spread(1, dx[ZZ]);
    tensor_product_2d_();
}

void
GaussTransform::addTransform(const rvec x, real weight)
{
    prepare_2d_grid(x, weight);
    tensor_product_();
}


GaussTransformOneD::GaussTransformOneD(int m_spread, real nu, real sigma) : m_spread_(m_spread), nu_(nu), sigma_(sigma), spread_(2*m_spread+1), E3_(2*m_spread+1)
{
    for (int i = -m_spread_; i <= m_spread_; ++i)
    {
        E3_[i+m_spread_] = exp( -i * i * nu_ * nu_ / 2. );
    }
}

ArrayRef<real>  GaussTransformOneD::spread(real weight, real dx)
{
    real E1 = weight* 1./(sqrt(2*M_PI*sigma_))*exp(-dx*dx/2.0); //< weight * nu / sqrt(2*pi) * exp(-dx_*dx_/2) , following the naming convention of Greengard et al., ;
    // avoid over / underflow for E2_power operation
    // this only occurs when the Gaussian width is much smaller than the grid
    // maximum expected E2 ^ m_spread may never be larger than max float
    // which in turn requires exp(dx * nu)^m_spread <= max_float, or
    // m_spread * dx * nu < log(max_float)
    if (m_spread_ * dx * nu_ > 86)
    {
        spread_[m_spread_] = E1;
        for (int l = 1; l <= m_spread_; l++)
        {
            spread_[m_spread_-l] = 0;
            spread_[m_spread_+l] = 0;
        }
        return spread_;
    }

    real E2         = exp(dx*nu_); //< following the naming convention of Greengard et al., ;
    real E2_power_l = E2;

    spread_[m_spread_] = E1;
    for (int l = 1; l < m_spread_; l++)
    {
        spread_[m_spread_-l] = (E1 / E2_power_l) * E3_[m_spread_-l];
        spread_[m_spread_+l] = (E1 * E2_power_l) * E3_[m_spread_+l];

        E2_power_l *= E2;
    }
    // seperate statement for l = 0 or l= 2*m_spread_ avoids
    // E2 multiplication operation and overflow for very small nu
    spread_[0]           = (E1 / E2_power_l) * E3_[0];
    spread_[2*m_spread_] = (E1 * E2_power_l) * E3_[2*m_spread_];

    return spread_;
}

GaussIntergralOneD::GaussIntergralOneD(int m_spread, real nu, real sigma) :
    m_spread_(m_spread), nu_(nu), sigma_(sigma)
{
    int bufferSize = 2*m_spread_+2;
  #if GMX_SIMD_HAVE_REAL
    bufferSize = GMX_SIMD_REAL_WIDTH * (1+(bufferSize-1) / GMX_SIMD_REAL_WIDTH);
  #endif

    voxelBoundaryBuffer_.resize(bufferSize);
    erfBuffer_.resize(bufferSize);

    // coordinates of the left and right boundaries of the spread intervall,
    // note that the right boundary of a spread intervall is the left boundary
    // of the following interval
    real voxelCoordinate = nu * (-m_spread_- 0.5 - 1);
    std::generate(std::begin(voxelBoundaryBuffer_), std::end(voxelBoundaryBuffer_), [&voxelCoordinate, nu](){voxelCoordinate += nu; return voxelCoordinate; });
}

ArrayRef<real>  GaussIntergralOneD::spread( real weight, real x)
{
    // shift the voxel boundaries by the normalised distance of spread atom to neares voxel
    std::transform(voxelBoundaryBuffer_.begin(), voxelBoundaryBuffer_.end(), erfBuffer_.begin(), [x](real voxelBoundary) {return (voxelBoundary-x)/std::sqrt(2); });
#if GMX_SIMD_HAVE_REAL
    for (size_t i = 0; i < erfBuffer_.size(); i += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal sx = load<SimdReal>(erfBuffer_.data() + i);
        sx = weight * 0.5 * gmx::erf(sx);
        store(erfBuffer_.data() + i, sx);
    }
#else       // GMX_SIMD_HAVE_REAL
    std::transform(std::begin(erfBuffer_), std::end(erfBuffer_), std::begin(erfBuffer_), [weight](real x){return weight*0.5*std::erf(x); });
#endif      // GMX_SIMD_HAVE_REAL
    // the differences of the error function values at the voxel boundaries
    // are the integrals of the gaussian over the respective integral
    std::adjacent_difference(std::begin(erfBuffer_), std::end(erfBuffer_), std::begin(erfBuffer_));
    return ArrayRef<real>(erfBuffer_.data()+1, &(*erfBuffer_.end()));
}

}    // namespace gmx
