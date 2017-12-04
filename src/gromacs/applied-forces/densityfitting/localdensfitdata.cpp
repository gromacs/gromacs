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

#include "localdensfitdata.h"
#include "densityspreader.h"

#include <cassert>
#include <numeric>
#include <functional>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/mdlib/sim_util.h"


namespace gmx
{

namespace
{

struct erfFunctor{
    #if GMX_SIMD_HAVE_REAL
    gmx::SimdReal operator()(gmx::SimdReal x) {return gmx::erf(x); };
    #else   // GMX_SIMD_HAVE_REAL
    real operator()(real x) {return gmx::erf(x); };
    #endif  // GMX_SIMD_HAVE_REAL
};

struct expFunctor{
    #if GMX_SIMD_HAVE_REAL
    gmx::SimdReal operator()(gmx::SimdReal x) {return gmx::exp(x); };
    #else   // GMX_SIMD_HAVE_REAL
    real operator()(real x) {return gmx::exp(x); };
    #endif  // GMX_SIMD_HAVE_REAL
};

struct multipliesFunctor{
    #if GMX_SIMD_HAVE_REAL
    gmx::SimdReal operator()(gmx::SimdReal x, gmx::SimdReal y) {return x*y; };
    #else   // GMX_SIMD_HAVE_REAL
    real operator()(real x, real y) {return x*y; };
    #endif  // GMX_SIMD_HAVE_REAL
};

struct crossCorrelationFunctor{
    #if GMX_SIMD_HAVE_REAL
    gmx::SimdReal operator()(gmx::SimdReal x, gmx::SimdReal y) {return 0.5*(x*x - 2* x*y + y*y); };
    #else   // GMX_SIMD_HAVE_REAL
    real operator()(real x, real y) {return 0.5*(x*x - 2* x*y + y*y); };
    #endif  // GMX_SIMD_HAVE_REAL
};

struct crossCorrelationDerivative{
    // #if GMX_SIMD_HAVE_REAL
    // crossCorrelationDerivative(const real forceConstant)
    // {
    //     alignas(GMX_SIMD_ALIGNMENT) real        mem[GMX_SIMD_REAL_WIDTH];
    //
    //     for (int i = 0; i < GMX_SIMD_REAL_WIDTH; i++)
    //     {
    //         mem[i] = forceConstant;  // repeat vector contents to fill simd width
    //     }
    //     forceConstant_= load<SimdReal>(mem);
    // };
    // gmx::SimdReal operator()(gmx::SimdReal x, gmx::SimdReal y) {return forceConstant_ *(y-x); };
    // gmx::SimdReal forceConstant_;
    // #else   // GMX_SIMD_HAVE_REAL
    crossCorrelationDerivative(const real forceConstant) : forceConstant_ {forceConstant}
    {};
    real operator()(real x, real y) {return forceConstant_ *(y-x); };
    real forceConstant_;
    // #endif  // GMX_SIMD_HAVE_REAL
};

struct dividesFunctor{
    #if GMX_SIMD_HAVE_REAL
    gmx::SimdReal operator()(gmx::SimdReal x, gmx::SimdReal y) {return x/y; };
    #else   // GMX_SIMD_HAVE_REAL
    real operator()(real x, real y) {return x/y; };
    #endif  // GMX_SIMD_HAVE_REAL
};

template<class F>
void inline inPlaceTransform(PaddedArrayRef<real> valueArray)
{
#if GMX_SIMD_HAVE_REAL
    const auto simdAlignedSize = (valueArray.size()/GMX_SIMD_REAL_WIDTH) * GMX_SIMD_REAL_WIDTH;
    for (size_t i = 0; i < simdAlignedSize; i += GMX_SIMD_REAL_WIDTH)
    {
        auto sx = F() (load<SimdReal>(valueArray.data() + i));
        store(valueArray.data() + i, sx);
    }
    std::vector < real, AlignedAllocator < real>> tmp(GMX_SIMD_REAL_WIDTH, 0.);
    for (size_t i = simdAlignedSize; i < valueArray.size(); ++i)
    {
        tmp[i-simdAlignedSize] = valueArray[i];
    }

    std::vector < real, AlignedAllocator < real>> lastvalues(GMX_SIMD_REAL_WIDTH, 0.);
    store(lastvalues.data(), F() (load<SimdReal>(tmp.data())));
    for (size_t i = simdAlignedSize; i < valueArray.size(); ++i)
    {
        valueArray[i] = lastvalues[i-simdAlignedSize];
    }
#else       // GMX_SIMD_HAVE_REAL
    std::transform(std::begin(valueArray), std::end(valueArray), std::begin(valueArray), F());
#endif      // GMX_SIMD_HAVE_REAL
}

template<class F>
real inline innerProduct(PaddedArrayRef<const real> reference, PaddedArrayRef<const real> comparand)
{
#if GMX_SIMD_HAVE_REAL
    SimdReal   simdSum = 0.;

    const auto simdAlignedSize = (reference.size()/GMX_SIMD_REAL_WIDTH) * GMX_SIMD_REAL_WIDTH;
    for (size_t i = 0; i < simdAlignedSize; i += GMX_SIMD_REAL_WIDTH)
    {
        const SimdReal sx = load<SimdReal>(reference.data() + i);
        const SimdReal sy = load<SimdReal>(comparand.data() + i);
        simdSum = simdSum + F() (sx, sy);
    }
    std::vector < real, AlignedAllocator < real>> summationArray(GMX_SIMD_REAL_WIDTH);
    store(summationArray.data(), simdSum);
    real result = std::accumulate(std::begin(summationArray), std::end(summationArray), 0.);

    std::vector < real, AlignedAllocator < real>> referenceExtra(GMX_SIMD_REAL_WIDTH, 0.);
    std::vector < real, AlignedAllocator < real>> comparandExtra(GMX_SIMD_REAL_WIDTH, 0.);
    for (size_t i = simdAlignedSize; i < reference.size(); ++i)
    {
        referenceExtra[i-simdAlignedSize] = reference[i];
        comparandExtra[i-simdAlignedSize] = comparand[i];
    }
    std::vector < real, AlignedAllocator < real>> lastvalues(GMX_SIMD_REAL_WIDTH, 0.);
    store(lastvalues.data(), F() (load<SimdReal>(referenceExtra.data()), load<SimdReal>(comparandExtra.data())));
    for (size_t i = simdAlignedSize; i < reference.size(); ++i)
    {
        result += lastvalues[i-simdAlignedSize];
    }
    return result;
#else       // GMX_SIMD_HAVE_REAL
    return std::inner_product(std::begin(reference), std::end(reference), std::begin(comparand), 0., std::plus<real>(), F());
#endif      // GMX_SIMD_HAVE_REAL
}

}

//! \brief Allocate and initialize the arrays for assembled positions, forces
//! and atomic weights
LocalDensfitData::LocalDensfitData(LocalAtomSet localAtomSet, const DensfitData &parameters,
                                   const std::vector<RVec> &x_densfit_whole) :
    localAtoms(localAtomSet), densitySpreader_(
            parameters.referenceMap().getGrid(),
            std::max(1, gmx_omp_nthreads_get(emntDefault)),
            parameters.dist(),
            parameters.timeDependent().currentSigma(0),
            false)
{
    const auto &grid = parameters.referenceMap().getGrid();
    map_sim_.setGrid(grid.duplicate());
    derivativeMap_.setGrid(grid.duplicate());

    /* Initialize the possibly time-dependent parameters */
    updateTimeDependentData(0., parameters);

    const auto nAtoms = parameters.fittingGroup().ind_.size();

    const RVec zeroRVecVector({0, 0, 0});
    x_assembled.resize(nAtoms, zeroRVecVector);
    f_loc.resize(nAtoms, zeroRVecVector);
    xWholeMoleculeReference_.resize(nAtoms, zeroRVecVector);

    const IVec zeroIVecVector({0, 0, 0});
    x_shifts.resize(nAtoms, zeroIVecVector);
    extra_shifts.resize(nAtoms, zeroIVecVector);

    bUpdateShifts  = true;

    xWholeMoleculeReference_ = x_densfit_whole;

    weights.resize(nAtoms, 1.);

    voxrange =
        ceil(0.5 * sigma * parameters.dist() / grid.unitCell().basisVectorLength(0));
    vox_per_atom = 3 * (2 * voxrange + 2);
    const int ompthreads      = std::max(1, gmx_omp_nthreads_get(emntDefault));
#pragma omp parallel for ordered num_threads(ompthreads) schedule(static) \
    default(none)
    for (int i = 0; i < ompthreads; i++)
    {
#pragma omp ordered
        {
            erfVector.emplace_back(AlignedRealVector(vox_per_atom));
            expVector.emplace_back(AlignedRealVector(vox_per_atom));
        }
    }
}

void LocalDensfitData::updateTimeDependentData(real time, const DensfitData &parameters)
{
    sigma          = parameters.timeDependent().currentSigma(time);
    forceConstant_ = parameters.nStFit() * parameters.timeDependent().currentForceConstant(time);
}

void LocalDensfitData::do_forces()
{

    if (localAtoms.numAtomsLocal() <= 0)
    {
        return;
    }

    // cppcheck-suppress unreadVariable
    int gmx_unused nth = std::max(1, gmx_omp_nthreads_get(emntDefault)); // number of threads

/* Zero out the forces array TODO */
#pragma omp parallel for num_threads(nth) schedule(static)
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        clear_rvec(f_loc[l]);
    }
    const auto  &grid            = map_sim_.getGrid();
    const auto   nu              = (grid.unitCell().basisVectorLength(0)+grid.unitCell().basisVectorLength(1)+grid.unitCell().basisVectorLength(2))/(3.0 * sigma);
    const double dRhoDxPrefactor = -sqrt(2./M_PI);

#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(grid) default(none)
    /* Loop over all atoms of the density fitting group */
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        const auto &x              = x_assembled[localAtoms.collectiveIndex()[l]];
        auto        atomGridIndex_ = grid.coordinateToFloorMultiIndex({{ x[XX], x[YY], x[ZZ] }});

        // continue to next atom if far outside the grid
        if (atomGridIndex_[XX]+voxrange < 0 || atomGridIndex_[YY]+voxrange < 0 || atomGridIndex_[ZZ]+voxrange < 0)
        {
            continue;
        }
        if (atomGridIndex_[XX]-voxrange >= grid.lattice().extend()[XX] || atomGridIndex_[YY]-voxrange >= grid.lattice().extend()[YY] || atomGridIndex_[ZZ]-voxrange >= grid.lattice().extend()[ZZ])
        {
            continue;
        }

        // (x-nearest voxel)/sigma
        RVec        dx;
        for (size_t i = XX; i <= ZZ; ++i)
        {
            dx[i]  = (x[i]-grid.multiIndexToCoordinate(atomGridIndex_)[i])/sigma;
        }

        /* For erf() performance, we want all values for one atom in one linear array (not in three).
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
            for (int ix = -voxrange; ix <= voxrange + 1; ix++)
            {
                const double left = (ix - 0.5) * nu - dx[dim];
                expThreadLocal[i]   = -0.5 * left * left;
                erfThreadLocal[i++] = left / sqrtf(2);
            }
        }

        inPlaceTransform<erfFunctor>(erfThreadLocal);
        inPlaceTransform<expFunctor>(expThreadLocal);

        /* Transform in place to differences of adjacent valules. */
        std::adjacent_difference(std::begin(erfThreadLocal), std::end(erfThreadLocal), std::begin(erfThreadLocal));
        std::adjacent_difference(std::begin(expThreadLocal), std::end(expThreadLocal), std::begin(expThreadLocal));

        const auto iStart   = std::max(0, atomGridIndex_[XX]-voxrange);
        const auto iEnd     = std::min(grid.lattice().extend()[XX], atomGridIndex_[XX]+voxrange);
        const int  iyOffset = 2*(2*voxrange+2);
        const int  izOffset = 2*voxrange+2;
        RVec       force    = {0., 0., 0.};

        for (int iz = -voxrange; iz <= voxrange; iz++)
        {
            Grid<DIM>::MultiIndex currentIndex {{
                                                    atomGridIndex_[XX], atomGridIndex_[YY], atomGridIndex_[ZZ]+iz
                                                }};
            if (!grid.lattice().inLattice(currentIndex))
            {
                continue;
            }
            real erfz = erfThreadLocal[iz+voxrange+1+ izOffset];
            real expz = expThreadLocal[iz+voxrange+1+ izOffset];

            for (int iy = -voxrange; iy <= voxrange; iy++)
            {
                currentIndex[YY] = atomGridIndex_[YY]+iy;
                if (!grid.lattice().inLattice(currentIndex))
                {
                    continue;
                }

                const real erfy = erfThreadLocal[iy+voxrange+1+ iyOffset];
                const real expy = expThreadLocal[iy+voxrange+1+ iyOffset];

                const real erfy_erfz             = erfy * erfz;
                const real expy_erfz             = expy * erfz;
                const real erfy_expz             = erfy * expz;
                auto       derivativeMapIterator = derivativeMap_.iteratorAtMultiIndex({{iStart, atomGridIndex_[YY]+iy, atomGridIndex_[ZZ]+iz}});

                /* Calculate d/dx_l rho_sim_ijk */
                for (int ix = -voxrange; ix <= voxrange; ix++)
                {
                    if (ix+atomGridIndex_[XX] < iStart || ix+atomGridIndex_[XX] >= iEnd)
                    {
                        continue;
                    }

                    const double densityDerivative = (double)*derivativeMapIterator;

                    const real   erfx = erfThreadLocal[ix+voxrange+1];
                    const real   expx = expThreadLocal[ix+voxrange+1];

                    const auto   factor =  densityDerivative * dRhoDxPrefactor;

                    force[XX] += factor * expx * erfy_erfz;
                    force[YY] += factor * erfx * expy_erfz;
                    force[ZZ] += factor * erfx * erfy_expz;

                    derivativeMapIterator++;
                }
            }
        } /* end of loop over voxels of this atom */
        rvec_inc(f_loc[l], force);
    }     /* end of loop over atoms */
}

const GridDataReal3D &LocalDensfitData::simulatedMap() const
{
    return map_sim_;
}

std::string LocalDensfitData::infoString() const
{
    char buffer[STRLEN];
    sprintf(buffer, "%12g %12.5e %12.5e %12.7f %12.5e", 0., forceConstant_, sigma, goodnessOfFit_, -forceConstant_ * goodnessOfFit_);
    return buffer;
}

void LocalDensfitData::add_forces(rvec * f) const
{
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        /* Get the right index of the local force, since typically not all local
         * atoms are subject to density fitting forces */
        rvec_inc(f[localAtoms.localIndex()[l]], f_loc[l]);
    }
}

bool LocalDensfitData::afterDomainDecomposition() const
{
    return bUpdateShifts;
}

void LocalDensfitData::communicate( PaddedArrayRef<RVec> x, const t_commrec * cr, matrix box)
{
    communicate_group_positions(
            cr, as_vec_array(x_assembled.data()), as_vec_array(x_shifts.data()), as_vec_array(extra_shifts.data()),
            bUpdateShifts, as_rvec_array(x.data()), localAtoms.numAtomsGlobal(),
            localAtoms.numAtomsLocal(), localAtoms.localIndex().data(), localAtoms.collectiveIndex().data(), as_vec_array(xWholeMoleculeReference_.data()), box);

/* If bUpdateShifts was TRUE then the shifts have just been updated in
 * communicate_group_positions. We do not need to update the shifts until
 * the next NS step */
    bUpdateShifts = false;

}

void LocalDensfitData::triggerShiftUpdate()
{
    bUpdateShifts = true;
}


void LocalDensfitData::spreadAtoms(const t_commrec * cr, bool normalize)
{

    densitySpreader_.zero();
    densitySpreader_.setSigma(sigma);

    if (cr == nullptr)
    {
        ArrayRef<real> weightsRef {
            weights.data(), weights.data()+weights.size()
        };
        PaddedArrayRef<RVec> xRef {
            x_assembled.data(), x_assembled.data()+x_assembled.size()
        };
        densitySpreader_.spreadLocalAtoms(x_assembled, weightsRef);
        map_sim_ = densitySpreader_.getSpreadGrid();
        return;
    }

    int  istart;
    int  nspread;
    int  nnodes;

    if (PAR(cr))
    {
        nnodes = cr->nnodes - cr->npmenodes;

        nspread = ceil((real) localAtoms.numAtomsGlobal() / nnodes);
        istart  = cr->nodeid * nspread;

        if ((nnodes - 1) == cr->nodeid)
        {
            nspread = localAtoms.numAtomsGlobal()  - nspread * (nnodes - 1);
        }
    }
    else
    {
        istart  = 0;
        nspread = localAtoms.numAtomsGlobal();
    }

    assert(istart >= 0);
    assert(nspread >= 0);

    ArrayRef<real> weightsRef {
        weights.data()+istart, weights.data()+istart+nspread
    };
    PaddedArrayRef<RVec> xRef {
        x_assembled.data()+istart, x_assembled.data()+istart+nspread
    };
    densitySpreader_.spreadLocalAtoms(xRef, weightsRef);
    map_sim_ = densitySpreader_.getSpreadGrid();
    if (normalize)
    {
        containeroperation::divide(&map_sim_, containermeasure::norm(weights));
    }
}

void LocalDensfitData::sumSimulatedMap(const t_commrec * cr)
{
    gmx_sumf(map_sim_.size(), map_sim_.data(), cr);
}

void LocalDensfitData::densityDensityDerivative(const GridDataReal3D &referenceMap)
{
    std::transform(std::begin(referenceMap), std::end(referenceMap), std::begin(map_sim_), std::begin(derivativeMap_),
                   crossCorrelationDerivative(forceConstant_));
}

void LocalDensfitData::goodnessOfFit(const GridDataReal3D &referenceMap)
{
    goodnessOfFit_ = innerProduct<crossCorrelationFunctor>(referenceMap, map_sim_);
}

}
