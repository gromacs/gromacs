/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016,2017, by the GROMACS development team, led by
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

/*! \internal \file
 * \brief
 * Implements the Awh class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "awh.h"

#include <assert.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"

#include "bias.h"
#include "grid.h"
#include "internal.h"
#include "pointstate.h"

/*! \internal
 * \brief A bias and its coupling to the system.
 *
 * Currently the AWH couples to the system by mapping each
 * AWH bias to a pull coordinate. This could be more general.
 */
class BiasCoupledToSystem
{
    public:
        /*! \brief Couple a bias to a set of pull coordinates.
         *
         * \param[in] bias            An initialized bias.
         * \param[in] pullCoordIndex  The pull coordinate indices.
         */
        BiasCoupledToSystem(Bias *bias, const std::vector<int> &pullCoordIndex);

        /*! \brief Returns the bias.
         */
        const Bias *bias() const
        {
            return bias_.get();
        }

        /*! \brief Returns the bias.
         */
        Bias *bias()
        {
            return bias_.get();
        }

        /*! \brief Returns the pull coordinate index of the bias dimension.
         *
         * \param[in] dim  The dimension.
         */
        int pullCoordIndex(int dim)
        {
            return pullCoordIndex_[dim];
        }

    private:
        std::unique_ptr<Bias> bias_;           /**< The bias. */
        std::vector<int>      pullCoordIndex_; /**< The pull coordinates this bias acts on. */

        /* Here AWH can be extended to work on other coordinates than pull. */
};

/* Construct a bias with coupling to a pull coordinate. */
BiasCoupledToSystem::BiasCoupledToSystem(Bias                   *bias,
                                         const std::vector<int> &pullCoordIndex) :
    bias_(bias),
    pullCoordIndex_(pullCoordIndex)
{
    GMX_RELEASE_ASSERT(static_cast<size_t>(bias->ndim()) == pullCoordIndex.size(), "The bias dimensionality should match the number of pull coordinates.");
}

/* Construct an AWH at the start of a simulation. */
Awh::Awh(FILE                *fplog,
         const t_inputrec    *ir,
         const t_commrec     *cr,
         const awh_params_t  *awhParams,
         t_state             *state_global,
         struct pull_t       *pull_work,
         bool                 startingFromCheckpoint) :
    seed_(awhParams->seed),
    potentialOffset_(0)
{
    GMX_RELEASE_ASSERT(ir->pull != nullptr, "With AWH we should have pull parameters");

    if (fplog != nullptr)
    {
        please_cite(fplog, "Lindahl2014");
    }

    /* Initialize all the biases */
    double beta = 1/(BOLTZ*ir->opts.ref_t[0]);
    for (int k = 0; k < awhParams->nbias; k++)
    {
        const awh_bias_params_t &awhBiasParams   = awhParams->awh_bias_params[k];

        std::vector<int>         pullCoordIndex;
        std::vector<DimParams>   dimParams;
        for (int d = 0; d < awhBiasParams.ndim; d++)
        {
            int                        pullCoordIndexDim = awhBiasParams.dim_params[d].pull_coord_index;
            const t_pull_coord        &pullCoord         = ir->pull->coord[pullCoordIndexDim];
            double                     conversionFactor  = pull_coordinate_is_angletype(&pullCoord) ? DEG2RAD : 1;
            dimParams.emplace_back(DimParams(conversionFactor, pullCoord.k, beta));

            pullCoordIndex.push_back(pullCoordIndexDim);
        }

        /* Initialize the bias */
        Bias *bias = new Bias(fplog, cr, k, *awhParams, awhParams->awh_bias_params[k], dimParams, beta, ir->delta_t);

        biasCoupledToSystem_.emplace_back(BiasCoupledToSystem(bias, pullCoordIndex));
    }

    /* Need to register the AWH coordinates to be allowed to apply forces to the pull coordinates. */
    registerAwhWithPull(awhParams, pull_work);

    if (startingFromCheckpoint)
    {
        GMX_RELEASE_ASSERT(!MASTER(cr) || state_global->awhHistory != nullptr, "The master rank should have the history when starting from checkpoint.");

        restoreStateFromHistory(*state_global->awhHistory.get(), cr);
    }
    else if (MASTER(cr))
    {
        state_global->awhHistory = std::shared_ptr<AwhHistory>(new AwhHistory());
        initHistoryFromState(state_global->awhHistory.get());
    }
}

/* Peform an AWH update every MD step. */
real Awh::applyBiasForcesAndUpdateBias(struct pull_t          *pull_work,
                                       int                     ePBC,
                                       const t_mdatoms        *mdatoms,
                                       const matrix            box,
                                       rvec                   *force,
                                       tensor                  virial,
                                       const gmx_multisim_t   *ms,
                                       double                  t,
                                       gmx_int64_t             step,
                                       gmx_wallcycle          *wallcycle,
                                       FILE                   *fplog)
{
    wallcycle_start(wallcycle, ewcAWH);

    t_pbc  pbc;
    set_pbc(&pbc, ePBC, box);

    /* During the AWH update the potential can instantaneously jump due to either
       an bias update or moving the umbrella. The jumps are kept track of and
       subtracted from the potential in order to get a useful conserved energy quantity. */
    double potentialJump = 0;
    double biasPotential = potentialOffset_;

    for (auto &biasCTS : biasCoupledToSystem_)
    {
        Bias *bias = biasCTS.bias();

        /* Update the AWH coordinate values with those of the corresponding
         * pull coordinates.
         */
        for (int d = 0; d < bias->ndim(); d++)
        {
            double coordValue;
            get_pull_coord_value(pull_work, biasCTS.pullCoordIndex(d), &pbc,
                                 &coordValue);
            bias->setCoordValue(d, coordValue);
        }

        /* Perform an AWH biasing step: this means, at regular intervals,
         * sampling observables based on the input pull coordinate value,
         * setting the bias force and/or updating the AWH bias state.
         */
        awh_dvec biasForce;
        bias->calcForceAndUpdateBias(biasForce, &biasPotential, &potentialJump,
                                     ms, t, step, seed_, fplog);

        /* Keep track of the total potential shift needed to remove the potential jumps. */
        potentialOffset_ -= potentialJump;

        /* Communicate the bias force to the pull struct.
         * The bias potential is returned at the end of this function,
         * so that it can be added externally to the correct energy data block.
         */
        for (int d = 0; d < bias->ndim(); d++)
        {
            apply_external_pull_coord_force(pull_work, biasCTS.pullCoordIndex(d),
                                            biasForce[d], mdatoms, force, virial);
        }
    }

    wallcycle_stop(wallcycle, ewcAWH);

    return static_cast<real>(biasPotential);
}

/* Allocate and initialize an AWH history with the given AWH state. */
void Awh::initHistoryFromState(AwhHistory *awhHistory) const
{
    GMX_RELEASE_ASSERT(awhHistory != nullptr, "Should be called with a valid history struct (and only on the master rank)");

    awhHistory->bias.clear();
    awhHistory->bias.resize(biasCoupledToSystem_.size());

    for (size_t k = 0; k < awhHistory->bias.size(); k++)
    {
        AwhBiasHistory *biasHistory = &awhHistory->bias[k];
        biasHistory->pointState.resize(biasCoupledToSystem_[k].bias()->pointState().size());
    }
}

/* Restore the AWH state from the given history. */
void Awh::restoreStateFromHistory(const AwhHistory &awhHistory,
                                  const t_commrec  *cr)
{
    /* Restore the history to the current state */
    if (MASTER(cr))
    {
        GMX_RELEASE_ASSERT(awhHistory.bias.size() == biasCoupledToSystem_.size(), "AWH state and history bias count should match");

        potentialOffset_ = awhHistory.potentialOffset;
    }
    if (PAR(cr))
    {
        gmx_bcast(sizeof(potentialOffset_), &potentialOffset_, cr);
    }

    for (size_t k = 0; k < biasCoupledToSystem_.size(); k++)
    {
        Bias *bias = biasCoupledToSystem_[k].bias();

        if (MASTER(cr))
        {
            bias->restoreStateFromHistory(&awhHistory.bias[k]);
        }

        if (PAR(cr))
        {
            bias->broadcast(cr);
        }
    }
}

/* Update the AWH history for checkpointing. */
void Awh::updateHistory(AwhHistory *awhHistory) const
{
    /* This assert will also catch a non-master rank calling this function. */
    GMX_RELEASE_ASSERT(awhHistory->bias.size() == biasCoupledToSystem_.size(), "AWH state and history bias count should match");

    awhHistory->potentialOffset = potentialOffset_;

    for (size_t k = 0; k < awhHistory->bias.size(); k++)
    {
        biasCoupledToSystem_[k].bias()->updateHistory(&awhHistory->bias[k]);
    }
}

/* Register the AWH biased coordinates with pull. */
void Awh::registerAwhWithPull(const awh_params_t *awhParams,
                              struct pull_t      *pull_work)
{
    for (int k = 0; k < awhParams->nbias; k++)
    {
        const awh_bias_params_t &biasParams = awhParams->awh_bias_params[k];

        for (int d = 0; d < biasParams.ndim; d++)
        {
            register_external_pull_potential(pull_work, biasParams.dim_params[d].pull_coord_index, Awh::externalPotentialString());
        }
    }
}
