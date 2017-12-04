/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016, by the GROMACS development team, led by
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
#include "gausstransform.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/math/griddata/canonicalvectorbasis.h"
#include "gromacs/utility/compare.h"

namespace gmx
{

void
GaussTransform::set_sigma(real sigma)
{
    sigma_ = sigma;
};

void
GaussTransform::set_n_sigma(real n_sigma)
{
    n_sigma_ = n_sigma;
};

IVec
GaussTransform::getMinimumUsedGridIndex()
{
    return minimumUsedGridIndex_;
};

IVec
GaussTransform::getMaximumUsedGridIndex()
{
    return maximumUsedGridIndex_;
};


/*! \brief
 * True, if length is same in x,y,z -direction.
 */
bool FastGaussianGridding::allVectorsSameLength(const CanonicalVectorBasis<DIM> &basis, real ftol, real abstol) const
{
    auto firstLength           = basis.basisVectorLength(0);
    for (int dimension = 1; dimension < DIM; ++dimension)
    {
        if (!equal_real(basis.basisVectorLength(dimension), firstLength, ftol, abstol))
        {
            return false;
        }
    }
    ;
    return true;
}

void
FastGaussianGridding::set_grid(std::unique_ptr < GridDataReal3D> grid)
{
    grid_ = std::move((grid));
    if (!allVectorsSameLength(grid_->getGrid().unitCell(), 1e-6, 1))
    {
        GMX_THROW(gmx::InconsistentInputError("Grid needs to be evently spaced to use the current implementation of fast gaussian gridding."));
    }
    nu_       = grid_->getGrid().unitCell().basisVectorLength(XX)/sigma_;
    m_spread_ = int(ceil(n_sigma_/nu_)); // number of grid cells for spreading
    E3_.resize(2*m_spread_+1);

    for (int i = -m_spread_; i <= m_spread_; ++i)
    {
        E3_[i+m_spread_] = exp( -i * i * nu_ * nu_ / 2. );
    }
    minimumUsedGridIndex_ = {GMX_INT32_MAX, GMX_INT32_MAX, GMX_INT32_MAX};
    maximumUsedGridIndex_ = {0, 0, 0};
};

void
FastGaussianGridding::spread_1d(real weight, int m_spread, rvec dx, real nu, const std::vector<real> &E3, int dimension)
{
    spread_1d_[dimension].resize(2*m_spread+1);

    real E1         = weight*1./(sqrt(2*M_PI)*sigma_)*exp(-dx[dimension]*dx[dimension]/2.0); //< weight * nu / sqrt(2*pi) * exp(-dx_*dx_/2) , following the naming convention of Greengard et al., ;
    real E2         = exp(dx[dimension]*nu);                                                 //< exp(dx_*nu_) , following the naming convention of Greengard et al., ;
    real E2_power_l = E2;

    spread_1d_[dimension][m_spread] = E1;
    for (int l = 1; l <= m_spread; l++)
    {
        spread_1d_[dimension][m_spread-l] = (E1 / E2_power_l) * E3[m_spread-l];
        spread_1d_[dimension][m_spread+l] = (E1 * E2_power_l) * E3[m_spread+l];

        E2_power_l *= E2;
    }
};

void
FastGaussianGridding::tensor_product_2d_()
{
    spread_2d_.resize(2*m_spread_+1);
    real spread_z;
    for (int l_z = 0; l_z < 2 * m_spread_ + 1; ++l_z)
    {
        spread_2d_[l_z].resize(2*m_spread_+1);
        spread_z = spread_1d_[ZZ][l_z];
        for (int l_y = 0; l_y < 2 * m_spread_ + 1; ++l_y)
        {
            spread_2d_[l_z][l_y] = spread_z*spread_1d_[YY][l_y];
        }
    }
}


void
FastGaussianGridding::tensor_product_()
{

    IVec minimumGlobalGridIndex;
    IVec maximumGlobalGridIndex;
    for (size_t i = XX; i <= ZZ; ++i)
    {
        minimumGlobalGridIndex[i] = std::max(0, grid_index_of_spread_atom_[i]-m_spread_);
        minimumUsedGridIndex_[i]  = std::min(minimumGlobalGridIndex[i], minimumUsedGridIndex_[i]);
        maximumGlobalGridIndex[i] = std::min(grid_index_of_spread_atom_[i]+m_spread_, grid_->getGrid().lattice().extend()[i]-1);
        maximumUsedGridIndex_[i]  = std::max(maximumGlobalGridIndex[i], maximumUsedGridIndex_[i]);
    }

    real   spread_zy;
    std::vector<real>::iterator voxel;
    real * spread_1d_XX;
    std::vector < std::vector < int>> ceilSqrtLUT(m_spread_+1);
    for (int d_z = 0; d_z <= m_spread_; d_z++)
    {
        ceilSqrtLUT[d_z].resize( (int)std::ceil(sqrt(m_spread_*m_spread_ - d_z*d_z))+1);
        for (int d_y = 0; d_y*d_y <= m_spread_*m_spread_ - d_z*d_z; d_y++)
        {
            ceilSqrtLUT[d_z][d_y] = (int)std::ceil(sqrt(m_spread_*m_spread_ - d_z*d_z- d_y*d_y));
        }
    }
    #pragma omp simd
    for (int globalGridIndexZZ = minimumGlobalGridIndex[ZZ]; globalGridIndexZZ <= maximumGlobalGridIndex[ZZ]; ++globalGridIndexZZ)
    {
        int d_z              = globalGridIndexZZ - grid_index_of_spread_atom_[ZZ];
        int localGridIndexZZ = d_z + m_spread_;
        // spread spheres instead of cubes
        int globalGridIndexYYStart = std::max(minimumGlobalGridIndex[YY], grid_index_of_spread_atom_[YY] - ceilSqrtLUT[std::abs(d_z)][0]);
        int globalGridIndexYYEnd   = std::min(maximumGlobalGridIndex[YY], grid_index_of_spread_atom_[YY] + ceilSqrtLUT[std::abs(d_z)][0]);

        for (int globalGridIndexYY = globalGridIndexYYStart; globalGridIndexYY <= globalGridIndexYYEnd; ++globalGridIndexYY)
        {
            int d_y              = globalGridIndexYY - grid_index_of_spread_atom_[YY];
            int localGridIndexYY = d_y + m_spread_;
            spread_zy    = spread_2d_[localGridIndexZZ][localGridIndexYY];
            int globalGridIndexXXStart = std::max(minimumGlobalGridIndex[XX], grid_index_of_spread_atom_[XX] - ceilSqrtLUT[std::abs(d_z)][std::abs(d_y)]);
            int globalGridIndexXXEnd   = std::min(maximumGlobalGridIndex[XX], grid_index_of_spread_atom_[XX] + ceilSqrtLUT[std::abs(d_z)][std::abs(d_y)]);
            int localGridIndexXXStart  = globalGridIndexXXStart - grid_index_of_spread_atom_[XX]+m_spread_;
            int numberSpreadVoxelsXX   = globalGridIndexXXEnd-globalGridIndexXXStart;
            voxel        = grid_->iteratorAtMultiIndex({{globalGridIndexXXStart, globalGridIndexYY, globalGridIndexZZ}});
            spread_1d_XX = &(spread_1d_[XX][localGridIndexXXStart]);

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
    grid_index_of_spread_atom_ = grid_->getGrid().coordinateToFloorMultiIndex({{ x[XX], x[YY], x[ZZ] }});
    RVec dx; // (x-nearest voxel)/sigma
    for (size_t i = XX; i <= ZZ; ++i)
    {
        dx[i]  = (x[i]-grid_->getGrid().multiIndexToCoordinate(grid_index_of_spread_atom_)[i])/sigma_;
    }
    spread_1d(weight, m_spread_, dx, nu_, E3_, XX);
    spread_1d(1, m_spread_, dx, nu_, E3_, YY);
    spread_1d(1, m_spread_, dx, nu_, E3_, ZZ);
    tensor_product_2d_();
}

void
FastGaussianGridding::transform(const rvec x, real weight)
{
    prepare_2d_grid(x, weight);
    tensor_product_();
}

std::unique_ptr < GridDataReal3D>
FastGaussianGridding::finish_and_return_grid()
{
    return std::move(grid_);
}

}    // namespace gmx
