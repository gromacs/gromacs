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

#include "gromacs/math/griddata/operations/densityspreader.h"
#include "localdensfitdata.h"
#include "measures.h"

#include <cassert>
#include <numeric>
#include <functional>

#include "densfitforces.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/math/container/containeroperation.h"
#include "gromacs/mdlib/gmx_omp_nthreads.h"
#include "gromacs/mdlib/groupcoord.h"

namespace gmx
{

//! \brief Allocate and initialize the arrays for assembled positions, forces
//! and atomic weights
LocalDensfitData::LocalDensfitData(LocalAtomSet localAtomSet, const DensfitData &parameters,
                                   const std::vector<RVec> &x_densfit_whole) :
    localAtoms(localAtomSet), densitySpreader_(
            parameters.referenceMap().getGrid(),
            std::max(1, gmx_omp_nthreads_get(emntDefault)),
            parameters.dist(),
            parameters.timeDependent().currentSigma(0),
            true)
{
    /* Initialize the possibly time-dependent parameters */
    updateTimeDependentData(0., parameters);

    measure_ = createDensityDensityMeasure(parameters.densityPotential(), parameters.referenceMap());

    const auto &grid = parameters.referenceMap().getGrid();
    densfitForces_ = std::unique_ptr<DensfitForces>(new DensfitForces(grid, sigma, parameters.dist()));
    map_sim_.setGrid(grid.duplicate());

    const auto nAtoms = parameters.fittingGroup().ind_.size();

    const RVec zeroRVecVector({0, 0, 0});
    x_assembled.resize(nAtoms, zeroRVecVector);
    f_loc.resize(nAtoms, zeroRVecVector);
    xWholeMoleculeReference_.resize(nAtoms, zeroRVecVector);

    const IVec zeroIVecVector({0, 0, 0});
    x_shifts.resize(nAtoms, zeroIVecVector);
    extra_shifts.resize(nAtoms, zeroIVecVector);

    bUpdateShifts            = true;
    xWholeMoleculeReference_ = x_densfit_whole;
    weights.resize(nAtoms, 1.);
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

    int gmx_unused nth = std::max(1, gmx_omp_nthreads_get(emntDefault)); // number of threads

/* Zero out the forces array TODO */
#pragma omp parallel for num_threads(nth) schedule(static)
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        clear_rvec(f_loc[l]);
    }
#pragma omp parallel for num_threads(nth) schedule(static) default(none)
    for (size_t l = 0; l < localAtoms.numAtomsLocal(); l++)
    {
        const RVec force = densfitForces_->force(x_assembled[localAtoms.collectiveIndex()[l]], derivativeMap_);
        rvec_inc(f_loc[l], force);
    }
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

    int  istart  = 0;
    int  nspread =  localAtoms.numAtomsGlobal();

    if (PAR(cr))
    {
        int nnodes = cr->nnodes - cr->npmenodes;

        nspread = ceil((real) localAtoms.numAtomsGlobal() / nnodes);
        istart  = cr->nodeid * nspread;

        if ((nnodes - 1) == cr->nodeid)
        {
            nspread = localAtoms.numAtomsGlobal()  - nspread * (nnodes - 1);
        }
    }

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

void LocalDensfitData::densityDensityDerivative()
{
    measure_->evaluateDensityDensityDerivative(map_sim_, forceConstant_);
}

void LocalDensfitData::goodnessOfFit()
{
    goodnessOfFit_ = measure_->goodnessOfFit(map_sim_);
}

}
