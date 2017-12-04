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

#include "gmxpre.h"

#include "densityspreader.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/math/griddata/operations/gausstransform.h"
#include "gromacs/math/vectypes.h"
#include <algorithm>

namespace gmx
{

DensitySpreader::DensitySpreader(const IGrid<DIM> &grid, int numberOfThreads, int n_sigma, real sigma, bool integrate) :
    simulated_density_(grid.duplicate()), number_of_threads_(numberOfThreads), integrate_(integrate)
{

#pragma omp parallel for ordered num_threads(number_of_threads_)
    for (int thread = 0; thread < number_of_threads_; ++thread)
    {
    #pragma omp critical
        {
            simulated_density_buffer_.emplace_back(GridDataReal3D(grid.duplicate()));
            gauss_transform_.emplace_back(GaussTransform());
        }
    }
    for (auto &transform : gauss_transform_)
    {
        transform.setSigma(sigma);
        transform.setSpreadWidthInNSigma(n_sigma);
    }
}

DensitySpreader::~DensitySpreader()
{
};


void DensitySpreader::setSigma(real sigma)
{
    for (auto &transform : gauss_transform_)
    {
        transform.setSigma(sigma);
    }
}


const GridDataReal3D &
DensitySpreader::getSpreadGrid() const
{
    return simulated_density_;
}

void DensitySpreader::zero()
{
    containeroperation::setZero(&simulated_density_);
}

void DensitySpreader::zeroBuffers(int thread)
{
    auto     &densitybuffer = simulated_density_buffer_[thread];
    auto      minIndex      = gauss_transform_[thread].minimumUsedGridIndex();
    auto      maxIndex      = gauss_transform_[thread].maximumUsedGridIndex();
    for (int iz = minIndex[ZZ]; iz < maxIndex[ZZ]; ++iz)
    {
        for (int iy = minIndex[YY]; iy < maxIndex[YY]; iy++)
        {
            std::fill(densitybuffer.iteratorAtMultiIndex({{minIndex[XX], iy, iz}}), densitybuffer.iteratorAtMultiIndex({{maxIndex[XX], iy, iz}}), 0.);
        }
    }
}


const GridDataReal3D &
DensitySpreader::spreadLocalAtoms(PaddedArrayRef<const RVec> x, ArrayRef<const real> weights)
{
    std::vector<IVec> minimumUsedGridIndex(number_of_threads_);
    std::vector<IVec> maximumUsedGridIndex(number_of_threads_);

    const size_t      nAtoms = x.size();

#pragma omp parallel num_threads(number_of_threads_)
    {
        const int thread        = gmx_omp_get_thread_num();

        zeroBuffers(thread);
        gauss_transform_[thread].setGrid(&simulated_density_buffer_[thread], integrate_);

        const size_t beginThreadAtoms = thread * (nAtoms / number_of_threads_ );
        const size_t endThreadAtoms   = (thread == number_of_threads_-1) ?  nAtoms : (thread+1) * ( nAtoms / number_of_threads_);

        for (size_t atomIndex = beginThreadAtoms; atomIndex != endThreadAtoms; ++atomIndex)
        {
            gauss_transform_[thread].addTransform(x[atomIndex], weights[atomIndex]);
        }

        minimumUsedGridIndex[thread] = gauss_transform_[thread].minimumUsedGridIndex();
        maximumUsedGridIndex[thread] = gauss_transform_[thread].maximumUsedGridIndex();
    }
    return sumThreadLocalGrids_(minimumUsedGridIndex, maximumUsedGridIndex);
}

const GridDataReal3D &
DensitySpreader::sumThreadLocalGrids_(const std::vector<IVec> &minimumUsedGridIndex, const std::vector<IVec> &maximumUsedGridIndex)
{

    std::array<int, 3> gridStart;
    std::array<int, 3> gridEnd;
    for (int i = XX; i <= ZZ; ++i)
    {
        auto ithElementLarger = [i](IVec a, IVec b){
                return a[i] < b[i];
            };
        gridStart[i] = (*std::min_element(std::begin(minimumUsedGridIndex), std::end(minimumUsedGridIndex), ithElementLarger))[i];
        gridEnd[i]   = (*std::max_element(std::begin(maximumUsedGridIndex), std::end(maximumUsedGridIndex), ithElementLarger))[i];
    }
    // add together all thread-local grids, using only density values, where there is
    // actual spread density
    int  nGridPointsXX        = gridEnd[XX]-gridStart[XX];
    for (int gridIndexZZ = gridStart[ZZ]; gridIndexZZ <= gridEnd[ZZ]; ++gridIndexZZ)
    {
        for (int gridIndexYY = gridStart[YY]; gridIndexYY <= gridEnd[YY]; ++gridIndexYY)
        {
            const std::array<int, DIM>            xRowStartIndex({{gridStart[XX], gridIndexYY, gridIndexZZ}});
            std::vector<GridDataReal3D::iterator> contributingThreadVoxels;
            // access rows in thread local grids if they contribute to the density
            for (int thread = 0; thread < number_of_threads_; ++thread)
            {
                if ((minimumUsedGridIndex[thread][ZZ] <= gridIndexZZ) &&
                    ( gridIndexZZ <= maximumUsedGridIndex[thread][ZZ] ) &&
                    (minimumUsedGridIndex[thread][YY] <= gridIndexYY) &&
                    (gridIndexYY <= maximumUsedGridIndex[thread][YY]))
                {
                    contributingThreadVoxels.push_back(simulated_density_buffer_[thread].iteratorAtMultiIndex(xRowStartIndex));
                }
            }

            auto simulatedDensityVoxelIterator = simulated_density_.iteratorAtMultiIndex(xRowStartIndex);
            // step though grid row by row
            for (int gridIndexXX = 0; gridIndexXX <= nGridPointsXX; ++gridIndexXX)
            {
                // loop over the local threads, collect voxel contributions and advance all threadlocal iterators to next voxel in row
                for (auto &threadLocalVoxel : contributingThreadVoxels)
                {
                    *simulatedDensityVoxelIterator += *threadLocalVoxel;
                    ++threadLocalVoxel;
                }
                ++simulatedDensityVoxelIterator;
            }
        }
    }
    return simulated_density_;
};



} /* gmx */
