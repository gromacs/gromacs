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

#include <cmath>

#include "gausstransform.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/griddata/canonicalvectorbasis.h"
#include "gromacs/utility/compare.h"

namespace gmx
{
/****************************GaussTransform**********************************/
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

/****************************FastGaussianGridding**********************************/


/*! \brief
 * True, if length is same in x,y,z -direction.
 */
bool FastGaussianGridding::allVectorsSameLength(const CanonicalVectorBasis<DIM> &basis, real ftol, real abstol) const
{
    auto firstLength = basis.basisVectorLength(0);
    return equal_real(basis.basisVectorLength(YY), firstLength, ftol, abstol) && equal_real(basis.basisVectorLength(ZZ), firstLength, ftol, abstol);
}

void
FastGaussianGridding::setGrid(GridDataReal3D * grid)
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
    E3_.resize(2*m_spread_+1);

    spread_2d_.resize(2*m_spread_+1);
    for (int l_z = 0; l_z < 2 * m_spread_ + 1; ++l_z)
    {
        spread_2d_[l_z].resize(2*m_spread_+1);
    }

    for (int i = -m_spread_; i <= m_spread_; ++i)
    {
        E3_[i+m_spread_] = exp( -i * i * nu_ * nu_ / 2. );
    }

    ceilSqrtLUT_.resize(m_spread_+1);
    for (int zDistance = 0; zDistance <= m_spread_; zDistance++)
    {
        ceilSqrtLUT_[zDistance].resize( (int)std::ceil(sqrt(m_spread_*m_spread_ - zDistance*zDistance))+1);
        for (int yDistance = 0; yDistance*yDistance <= m_spread_*m_spread_ - zDistance*zDistance; yDistance++)
        {
            ceilSqrtLUT_[zDistance][yDistance] = (int)std::ceil(sqrt(m_spread_*m_spread_ - zDistance*zDistance- yDistance*yDistance));
        }
    }


};

void
FastGaussianGridding::spread_1d(real weight, int m_spread, rvec dx, real nu, const std::vector<real> &E3, int dimension)
{
    spread_1d_[dimension].resize(2*m_spread+1);

    real E1         = weight* 1./(sqrt(2*M_PI*sigma_))*exp(-dx[dimension]*dx[dimension]/2.0); //< weight * nu / sqrt(2*pi) * exp(-dx_*dx_/2) , following the naming convention of Greengard et al., ;
    // avoid over / underflow for E2_power operation
    // this only occurs when the Gaussian width is much smaller than the grid
    // maximum expected E2 ^ m_spread may never be larger than max float
    // which in turn requires exp(dx * nu)^m_spread <= max_float, or
    // m_spread * dx * nu < log(max_float)
    if (m_spread * dx[dimension] * nu > 86)
    {
        spread_1d_[dimension][m_spread] = E1;
        for (int l = 1; l <= m_spread; l++)
        {
            spread_1d_[dimension][m_spread-l] = 0;
            spread_1d_[dimension][m_spread+l] = 0;
        }
        return;
    }

    real E2         = exp(dx[dimension]*nu);                                                 //< following the naming convention of Greengard et al., ;
    real E2_power_l = E2;

    spread_1d_[dimension][m_spread] = E1;
    for (int l = 1; l < m_spread; l++)
    {
        spread_1d_[dimension][m_spread-l] = (E1 / E2_power_l) * E3[m_spread-l];
        spread_1d_[dimension][m_spread+l] = (E1 * E2_power_l) * E3[m_spread+l];

        E2_power_l *= E2;
    }
    // repeat statement to avoid unnecessary multiplication operation above
    // avoids overflow issues for very small Gaussian widths
    spread_1d_[dimension][0]          = (E1 / E2_power_l) * E3[0];
    spread_1d_[dimension][2*m_spread] = (E1 * E2_power_l) * E3[2*m_spread];

};

void
FastGaussianGridding::tensor_product_2d_()
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
FastGaussianGridding::tensor_product_()
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

void FastGaussianGridding::prepare_2d_grid(const rvec x, const real weight)
{
    atomGridIndex_ = grid_->getGrid().coordinateToFloorMultiIndex({{ x[XX], x[YY], x[ZZ] }});
    RVec dx; // (x-nearest voxel)/sigma
    for (size_t i = XX; i <= ZZ; ++i)
    {
        dx[i]  = (x[i]-grid_->getGrid().multiIndexToCoordinate(atomGridIndex_)[i])/sigma_;
    }
    spread_1d(weight, m_spread_, dx, nu_, E3_, XX);
    spread_1d(1, m_spread_, dx, nu_, E3_, YY);
    spread_1d(1, m_spread_, dx, nu_, E3_, ZZ);
    tensor_product_2d_();
}

void
FastGaussianGridding::addTransform(const rvec x, real weight)
{
    prepare_2d_grid(x, weight);
    tensor_product_();
}


/****************************IntegratedGaussian**********************************/

void IntegratedGaussian::setGrid(GridDataReal3D *grid)
{

}

void IntegratedGaussian::addTransform(const real *x, real weight)
{

}



}    // namespace gmx
