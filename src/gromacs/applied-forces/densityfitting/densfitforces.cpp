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
#include "gmxpre.h"

#include "densfitforces.h"

#include <functional>

#include "gromacs/math/griddata/griddata.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"

namespace gmx
{

namespace
{

struct erfFunctor {
#if GMX_SIMD_HAVE_REAL
    gmx::SimdReal operator()(gmx::SimdReal x) { return gmx::erf(x); };
#else  // GMX_SIMD_HAVE_REAL
    real operator()(real x) { return gmx::erf(x); };
#endif // GMX_SIMD_HAVE_REAL
};

struct expFunctor {
#if GMX_SIMD_HAVE_REAL
    gmx::SimdReal operator()(gmx::SimdReal x) { return gmx::exp(x); };
#else  // GMX_SIMD_HAVE_REAL
    real operator()(real x) { return gmx::exp(x); };
#endif // GMX_SIMD_HAVE_REAL
};

struct multipliesFunctor {
#if GMX_SIMD_HAVE_REAL
    SimdReal operator()(SimdReal x, SimdReal y) { return x * y; };
#else  // GMX_SIMD_HAVE_REAL
    real operator()(real x, real y) { return x * y; };
#endif // GMX_SIMD_HAVE_REAL
};

struct dividesFunctor {
#if GMX_SIMD_HAVE_REAL
    SimdReal operator()(SimdReal x, SimdReal y) { return x / y; };
#else  // GMX_SIMD_HAVE_REAL
    real operator()(real x, real y) { return x / y; };
#endif // GMX_SIMD_HAVE_REAL
};

template <class F>
void inline inPlaceTransform(PaddedArrayRef<real> valueArray)
{
#if GMX_SIMD_HAVE_REAL
    const auto simdAlignedSize =
        (valueArray.size() / GMX_SIMD_REAL_WIDTH) * GMX_SIMD_REAL_WIDTH;
    for (size_t i = 0; i < simdAlignedSize; i += GMX_SIMD_REAL_WIDTH)
    {
        auto sx = F() (load<SimdReal>(valueArray.data() + i));
        store(valueArray.data() + i, sx);
    }
    std::vector < real, AlignedAllocator < real>> tmp(GMX_SIMD_REAL_WIDTH, 0.);
    for (size_t i = simdAlignedSize; i < valueArray.size(); ++i)
    {
        tmp[i - simdAlignedSize] = valueArray[i];
    }

    std::vector < real, AlignedAllocator < real>> lastvalues(GMX_SIMD_REAL_WIDTH, 0.);
    store(lastvalues.data(), F() (load<SimdReal>(tmp.data())));
    for (size_t i = simdAlignedSize; i < valueArray.size(); ++i)
    {
        valueArray[i] = lastvalues[i - simdAlignedSize];
    }
#else       // GMX_SIMD_HAVE_REAL
    std::transform(std::begin(valueArray), std::end(valueArray),
                   std::begin(valueArray), F());
#endif      // GMX_SIMD_HAVE_REAL
}
}

DensfitForces::DensfitForces(const IGrid<DIM> &grid, real sigma, real nSigma) :
    grid_(grid.duplicate())
{
    setSigma(sigma, nSigma);
}

void DensfitForces::setSigma(real sigma, real nSigma)
{
    sigma_    = sigma;
    voxrange_ = ceil(0.5 * sigma * nSigma / grid_->unitCell().basisVectorLength(0));
    nu_       = (grid_->unitCell().basisVectorLength(0) + grid_->unitCell().basisVectorLength(1) + grid_->unitCell().basisVectorLength(2)) / (3.0 * sigma);

    const int voxPerAtom = 3 * (2 * voxrange_ + 2); /* Total number of voxels for one atom (x,y,z) */
    const int ompthreads = std::max(1, gmx_omp_nthreads_get(emntDefault));

    #pragma omp parallel for ordered num_threads(ompthreads) \
    schedule(static) default(none)
    for (int i = 0; i < ompthreads; i++)
    {
    #pragma omp ordered
        {
            erfVector.emplace_back(AlignedRealVector(voxPerAtom));
            expVector.emplace_back(AlignedRealVector(voxPerAtom));
        }
    }
}

RVec DensfitForces::force(const RVec           &x,
                          const GridDataReal3D &densityDensityDerivative)
{
    const IGrid<DIM> &grid           = *grid_;
    auto              atomGridIndex_ =
        grid.coordinateToFloorMultiIndex({{x[XX], x[YY], x[ZZ]}});

    // continue to next atom if far outside the grid
    if (atomGridIndex_[XX] + voxrange_ < 0 ||
        atomGridIndex_[YY] + voxrange_ < 0 ||
        atomGridIndex_[ZZ] + voxrange_ < 0)
    {
        return {};
    }
    if (atomGridIndex_[XX] - voxrange_ >= grid.lattice().extend()[XX] ||
        atomGridIndex_[YY] - voxrange_ >= grid.lattice().extend()[YY] ||
        atomGridIndex_[ZZ] - voxrange_ >= grid.lattice().extend()[ZZ])
    {
        return {};
    }

    // (x-nearest voxel)/sigma
    RVec dx;
    for (size_t i = XX; i <= ZZ; ++i)
    {
        dx[i] = (x[i] - grid.multiIndexToCoordinate(atomGridIndex_)[i]) / sigma_;
    }

    /* For erf() performance, we want all values for one atom in one linear array
     * (not in three).
     * Each OpenMP thread th accesses its own temporary storage vector
     */
    int                  th = gmx_omp_get_thread_num();
    PaddedArrayRef<real> erfThreadLocal(erfVector[th]);
    PaddedArrayRef<real> expThreadLocal(expVector[th]);

    /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
     * boundaries.*/
    int i = 0;
    for (int dim = XX; dim <= ZZ; ++dim)
    {
        for (int ix = -voxrange_; ix <= voxrange_ + 1; ix++)
        {
            const double left = (ix - 0.5) * nu_ - dx[dim];
            expThreadLocal[i]   = -0.5 * left * left;
            erfThreadLocal[i++] = left / sqrtf(2.);
        }
    }

    inPlaceTransform<erfFunctor>(erfThreadLocal);
    inPlaceTransform<expFunctor>(expThreadLocal);

    /* Transform in place to differences of adjacent valules. */
    std::adjacent_difference(std::begin(erfThreadLocal), std::end(erfThreadLocal),
                             std::begin(erfThreadLocal));
    std::adjacent_difference(std::begin(expThreadLocal), std::end(expThreadLocal),
                             std::begin(expThreadLocal));

    const auto iStart = std::max(0, atomGridIndex_[XX] - voxrange_);
    const auto iEnd   =
        std::min(grid.lattice().extend()[XX], atomGridIndex_[XX] + voxrange_);
    const int  iyOffset = 2 * voxrange_ + 2;
    const int  izOffset = 2 * (2 * voxrange_ + 2);
    real       forceX   = 0.;
    real       forceY   = 0.;
    real       forceZ   = 0.;

    for (int iz = -voxrange_; iz <= voxrange_; iz++)
    {
        Grid<DIM>::MultiIndex currentIndex {
            {
                atomGridIndex_[XX], atomGridIndex_[YY], atomGridIndex_[ZZ] + iz
            }
        };
        if (!grid.lattice().inLattice(currentIndex))
        {
            continue;
        }
        real erfz = erfThreadLocal[iz + voxrange_ + 1 + izOffset];
        real expz = expThreadLocal[iz + voxrange_ + 1 + izOffset];

        for (int iy = -voxrange_; iy <= voxrange_; iy++)
        {
            currentIndex[YY] = atomGridIndex_[YY] + iy;
            if (!grid.lattice().inLattice(currentIndex))
            {
                continue;
            }

            const real erfy = erfThreadLocal[iy + voxrange_ + 1 + iyOffset];
            const real expy = expThreadLocal[iy + voxrange_ + 1 + iyOffset];

            const real erfy_erfz             = erfy * erfz;
            const real expy_erfz             = expy * erfz;
            const real erfy_expz             = erfy * expz;
            auto       derivativeMapIterator =
                densityDensityDerivative.iteratorAtMultiIndex(
                        {{iStart, atomGridIndex_[YY] + iy, atomGridIndex_[ZZ] + iz}});

/* Calculate d/dx_l rho_sim_ijk */
#pragma omp simd
            for (int ix = -voxrange_; ix <= voxrange_; ix++)
            {
                if (ix + atomGridIndex_[XX] < iStart ||
                    ix + atomGridIndex_[XX] >= iEnd)
                {
                    continue;
                }

                const real erfx = erfThreadLocal[ix + voxrange_ + 1];
                const real expx = expThreadLocal[ix + voxrange_ + 1];

                const auto factor = *derivativeMapIterator * dRhoDxPrefactor_;

                forceX += factor * expx * erfy_erfz;
                forceY += factor * erfx * expy_erfz;
                forceZ += factor * erfx * erfy_expz;

                derivativeMapIterator++;
            }
        }
    }   /* end of loop over voxels of this atom */
    return {forceX, forceY, forceZ};
}

} // namespace gmx
