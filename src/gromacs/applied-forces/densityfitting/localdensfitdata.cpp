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

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/simd/simd.h"
#include "gromacs/simd/simd_math.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/math/griddata/griddata.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/mdlib/sim_util.h"


namespace gmx
{

namespace
{
//! \brief Exponential function exp for a vector of n elements
void gmx_exp_vector(int n, real x[])
{
    //    gmx_cycles_t c1 = gmx_cycles_read();

    #if GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal sx = load<SimdReal>(x + i);
        sx = gmx::exp(sx);
        store(x + i, sx);
    }
    #else       // GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i++)
    {
        x[i] = std::exp(x[i]);
    }
    #endif      // GMX_SIMD_HAVE_REAL
    //    gmx_cycles_t c2 = gmx_cycles_read();
    //    fprintf(stderr, "exp vector %15d elements %15lld cycles  %15.3f per
    //    exp\n", n, c2-c1, (1.0*(c2-c1))/n);
}

//! \brief The error function erf for a vector of n elements
void gmx_erf_vector(int n, real x[])
{
    //    gmx_cycles_t c1 = gmx_cycles_read();

 #if GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i += GMX_SIMD_REAL_WIDTH)
    {
        SimdReal sx = load<SimdReal>(x + i);
        sx = gmx::erf(sx);
        store(x + i, sx);
    }
 #else       // GMX_SIMD_HAVE_REAL
    for (int i = 0; i < n; i++)
    {
        x[i] = std::erf(x[i]);
    }
 #endif      // GMX_SIMD_HAVE_REAL
    //    gmx_cycles_t c2 = gmx_cycles_read();
    //    fprintf(stderr, "erf vector %15d elements %15lld cycles  %15.3f per erf
    //    (SIMD width is %d)\n", n, c2-c1, (1.0*(c2-c1))/n, GMX_SIMD_REAL_WIDTH);
}

}
/*! \brief Calculate one of the terms of the correlation coefficient
 *
 * Multiply each voxel of the map with the corresponding voxel of the other
 * map and sum up everything, i.e. compute
 *
 *  sum_ijk [ rhoA_ijk * rhoB_ijk ]
 *
 */
double calc_sum_rhoA_rhoB(const std::vector<float> &vecA,
                          const std::vector<float> &vecB)
{
    double sum_ijk = 0.;
    for (size_t i = 0; i < vecA.size(); i++)
    {
        sum_ijk += vecA[i] * vecB[i];
    }
    return sum_ijk;
}


//! \brief Allocate and initialize the arrays for assembled positions, forces
//! and atomic weights
LocalDensfitData::LocalDensfitData(LocalAtomSet localAtomSet, const DensfitData &parameters,
                                   rvec *x_densfit_whole) :
    localAtoms(localAtomSet), densitySpreader_(
            parameters.referenceMap().getGrid(),
            std::max(1, gmx_omp_nthreads_get(emntDefault)),
            parameters.dist(),
            parameters.timeDependent().currentSigma(0),
            false)
{
    const auto &grid = parameters.referenceMap().getGrid();
    map_sim_.setGrid(grid.duplicate());

    /* Initialize the possibly time-dependent parameters */
    updateTimeDependentData(0., parameters);

    const auto nAtoms = parameters.fittingGroup().ind_.size();

    const RVec zeroRVecVector({0, 0, 0});
    x_assembled.resize(nAtoms, zeroRVecVector);
    f_loc.resize(nAtoms, zeroRVecVector);
    x_old.resize(nAtoms, zeroRVecVector);

    const IVec zeroIVecVector({0, 0, 0});
    x_shifts.resize(nAtoms, zeroIVecVector);
    extra_shifts.resize(nAtoms, zeroIVecVector);

    bUpdateShifts  = true;

    /* In mdruns we have to make the density fitting structure whole */
    if (x_densfit_whole != nullptr)
    {
        for (size_t i = 0; i < nAtoms; i++)
        {
            copy_rvec(x_densfit_whole[i], x_old[i]);
        }
    }

    sum_rho_exp2 = calc_sum_rhoA_rhoB(parameters.referenceMap(), parameters.referenceMap());
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
            int extra = 0;
#if GMX_SIMD_HAVE_REAL
            extra = GMX_SIMD_REAL_WIDTH;
#endif
            erfVector.emplace_back(AlignedRealVector(vox_per_atom+extra));
            expVector.emplace_back(AlignedRealVector(vox_per_atom+extra));
        }
    }

}

void LocalDensfitData::updateTimeDependentData(real time, const DensfitData &parameters)
{
    sigma = parameters.timeDependent().currentSigma(time);
    k     = parameters.timeDependent().currentForceConstant(time);
}

void LocalDensfitData::do_forces(const GridDataReal3D &referenceMap)
{

    if (localAtoms.numAtomsLocal() <= 0)
    {
        return;
    }

    // cppcheck-suppress unreadVariable
    int gmx_unused nth =
        std::max(1, gmx_omp_nthreads_get(emntDefault)); // number of threads

/* Zero out the forces array TODO */
#pragma omp parallel for num_threads(nth) schedule(static)
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        clear_rvec(f_loc[l]);
    }
    const auto  &grid = map_sim_.getGrid();
    /* Calculate various prefactors */
    const double sum_rho_sim2 =
        calc_sum_rhoA_rhoB(map_sim_, map_sim_);
    const double sum_rho_exp_rho_sim =
        calc_sum_rhoA_rhoB(referenceMap, map_sim_);

    const double term1_prefactor = k / sqrt(sum_rho_exp2 * sum_rho_sim2);
    const double term2_prefactor = -k * sum_rho_exp_rho_sim /
        (sqrt(sum_rho_exp2) * pow(sum_rho_sim2, 1.5));
    const auto   spacing          = (grid.unitCell().basisVectorLength(0)+grid.unitCell().basisVectorLength(1)+grid.unitCell().basisVectorLength(2))/3.0;
    const double V_vox            = spacing * spacing * spacing;
    const double prefactor        = -M_PI * sigma * sigma / (6.0 * V_vox);

#pragma omp parallel for num_threads(nth) schedule(static) \
    shared(referenceMap,grid) default(none)
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
        RVec        dx; // (x-nearest voxel)/sigma
        for (size_t i = XX; i <= ZZ; ++i)
        {
            dx[i]  = (x[i]-grid.multiIndexToCoordinate(atomGridIndex_)[i])/sigma;
        }

        real *aarr[3];
        real *garr[3];

        int   th = gmx_omp_get_thread_num();

        dvec  term1_sum;
        dvec  term2_sum;

        /* Loop over all voxels and calculate term1 and term2 for the force
         * simultaneously */
        clear_dvec(term1_sum);
        clear_dvec(term2_sum);

        /* For erf() performance, we want all values for one atom in one linear
         * array (not
         * int three). However, we build pointers to access it as x, y, and z
         * arrays.
         *
         * Each OpenMP thread th accesses its own temporary storage vector
         */
        aarr[XX] = &erfVector[th][0];
        aarr[YY] = &erfVector[th][0 + 2 * voxrange + 2];
        aarr[ZZ] = &erfVector[th][0 + 2 * (2* voxrange + 2)];

        /* Same for the temporary array that keeps the exp values */
        garr[XX] = &expVector[th][0];
        garr[YY] = &expVector[th][0 + 2 * voxrange + 2];
        garr[ZZ] = &expVector[th][0 + 2 * (2* voxrange + 2)];

        /* Store the x,y,z-arguments for the erf(x,y,z) evaluations at the voxel
         * boundaries in the tmpgrid[] array. This way we can compute the erf's
         * for this atom all in one (vector) call. */
        for (int dim = XX; dim <= ZZ; ++dim)
        {
            int i = 0;
            for (int ix = -voxrange; ix <= voxrange + 1; ix++)
            {
                double left = (ix - 0.5) * spacing/sigma - dx[dim];
                garr[dim][i]   = -1.5 * left * left;
                aarr[dim][i++] = sqrtf(1.5) * left;
            }
        }

        /* Call erf() and exp() for all input values for this atom in one go */
        gmx_erf_vector(3*(2*voxrange +2 ), aarr[XX]);
        gmx_exp_vector(3*(2*voxrange +2 ), garr[XX]);

        /* Transform (in place) into a new array that directly stores the
         * differences of
         * the error functions, i.e. the area under the Gaussian curve */
        std::adjacent_difference(aarr[0], aarr[0] + 3*(2*voxrange +2), aarr[0]);
        std::adjacent_difference(garr[0], garr[0] + 3*(2*voxrange +2), garr[0]);

        const auto iStart = std::max(0, atomGridIndex_[XX]-voxrange);
        const auto iEnd   = std::min(grid.lattice().extend()[XX], atomGridIndex_[XX]+voxrange);

        for (int iz = -voxrange; iz <= voxrange; iz++)
        {
            if (!grid.lattice().inLattice({atomGridIndex_[XX], atomGridIndex_[YY], atomGridIndex_[ZZ]+iz}))
            {
                continue;
            }
            real erfz = aarr[ZZ][iz+voxrange+1];
            real expz = garr[ZZ][iz+voxrange+1];

            for (int iy = -voxrange; iy <= voxrange; iy++)
            {
                if (!grid.lattice().inLattice({atomGridIndex_[XX], atomGridIndex_[YY]+iy, atomGridIndex_[ZZ]+iz}))
                {
                    continue;
                }

                real erfy = aarr[YY][iy+voxrange+1];
                real expy = garr[YY][iy+voxrange+1];

                real erfy_erfz = erfy * erfz;
                real expy_erfz = expy * erfz;
                real erfy_expz = erfy * expz;
                auto ptr_ref   = referenceMap.iteratorAtMultiIndex({iStart, atomGridIndex_[YY]+iy, atomGridIndex_[ZZ]+iz});
                auto ptr_sim   = map_sim_.iteratorAtMultiIndex({iStart, atomGridIndex_[YY]+iy, atomGridIndex_[ZZ]+iz});

                /* Calculate d/dx_l rho_sim_ijk */
                for (int ix = -voxrange; ix <= voxrange; ix++)
                {
                    if (ix+atomGridIndex_[XX] < iStart || ix+atomGridIndex_[XX] >= iEnd)
                    {
                        continue;
                    }

                    const double ref_dens = (double)*ptr_ref;
                    const double sim_dens = (double)*ptr_sim;

                    const real   erfx = aarr[XX][ix+voxrange+1];
                    const real   expx = garr[XX][ix+voxrange+1];

                    rvec         drho_dxl;

                    drho_dxl[XX] = expx * erfy_erfz;
                    drho_dxl[YY] = erfx * expy_erfz;
                    drho_dxl[ZZ] = erfx * erfy_expz;

                    term1_sum[XX] += ref_dens * drho_dxl[XX];
                    term1_sum[YY] += ref_dens * drho_dxl[YY];
                    term1_sum[ZZ] += ref_dens * drho_dxl[ZZ];

                    term2_sum[XX] += sim_dens * drho_dxl[XX];
                    term2_sum[YY] += sim_dens * drho_dxl[YY];
                    term2_sum[ZZ] += sim_dens * drho_dxl[ZZ];

                    ptr_ref++;
                    ptr_sim++;
                }
            }
        }     /* end of loop over voxels of this atom */

        dvec tmp1;
        dvec tmp2;
        dsvmul(prefactor * term1_prefactor, term1_sum, tmp1);
        dsvmul(prefactor * term2_prefactor, term2_sum, tmp2);

        dvec force;
        dvec_add(tmp1, tmp2, force);

        f_loc[l][XX] = force[XX];
        f_loc[l][YY] = force[YY];
        f_loc[l][ZZ] = force[ZZ];

    }   /* end of loop over atoms */
}

const GridDataReal3D &LocalDensfitData::simulatedMap() const
{
    return map_sim_;
}


std::string LocalDensfitData::infoString() const
{
    char buffer[STRLEN];
    sprintf(buffer, "%12g%12.5e%12.5e%12.7f%12.5e", 0., k, sigma, cc, k * (1.0 - cc));
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

void LocalDensfitData::communicate( PaddedArrayRef<RVec> x, const t_commrec * cr, matrix box, int ePBC)
{
    communicate_group_positions(
            cr, as_vec_array(x_assembled.data()), as_vec_array(x_shifts.data()), as_vec_array(extra_shifts.data()),
            bUpdateShifts, as_rvec_array(x.data()), localAtoms.numAtomsGlobal(),
            localAtoms.numAtomsLocal(), localAtoms.localIndex().data(), localAtoms.collectiveIndex().data(), as_vec_array(x_old.data()), box);

/* If bUpdateShifts was TRUE then the shifts have just been updated in
 * communicate_group_positions. We do not need to update the shifts until
 * the next NS step */
    bUpdateShifts = false;

/* Put all atoms in the box (TODO: if we do that, we do not need to construct
 * a whole DF group before with communicate_group_positions!)
 */
    if (ePBC != epbcNONE)
    {
        auto xAssembledArrayRef = arrayRefFromArray(
                    x_assembled.data(), localAtoms.numAtomsGlobal());
        put_atoms_in_box_omp(ePBC, box, xAssembledArrayRef);
    }
}

void LocalDensfitData::triggerShiftUpdate()
{
    bUpdateShifts = true;
}


void LocalDensfitData::spreadAtoms(const t_commrec * cr)
{

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

}

void LocalDensfitData::sumSimulatedMap(const t_commrec * cr)
{
    gmx_sumf(map_sim_.size(), map_sim_.data(), cr);
}


void LocalDensfitData::calculateCorrelationCoefficient(const GridDataReal3D &referenceMap)
{
    cc = calc_sum_rhoA_rhoB(referenceMap, map_sim_) / sqrt(sum_rho_exp2 *
                                                           calc_sum_rhoA_rhoB(map_sim_, map_sim_));
}

}
