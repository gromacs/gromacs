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
 * Implements the Bias class.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include "bias.h"

#include <assert.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"

#include "grid.h"
#include "internal.h"
#include "math.h"
#include "pointstate.h"


/* Makes checks for the collected histograms and warns if issues are detected. */
void Bias::checkHistograms(double t, gmx_int64_t step, FILE *fplog)
{
    const int    maxNumWarningsInCheck = 1;   /* The maximum number of warnings to print per check */
    const int    maxNumWarningsInRun   = 10;  /* The maximum number of warnings to print in a run */

    if (fplog == nullptr || numWarningsIssued_ >= maxNumWarningsInRun || state_.inInitialStage() ||
        (step == 0) || isCheckStep(params_, state_.points(), step))
    {
        return;
    }

    numWarningsIssued_ +=
        state_.checkHistograms(grid(), biasIndex(), t, fplog,
                               maxNumWarningsInCheck);

    if (numWarningsIssued_ >= maxNumWarningsInRun)
    {
        fprintf(fplog, "\nawh%d: suppressing future AWH warnings.\n", biasIndex() + 1);
    }
}

/* Do all previously skipped updates. */
void Bias::doSkippedUpdatesForAllPoints()
{
    state_.doSkippedUpdatesForAllPoints(params_);
}

/* Do what the AWH bias needs to do every step. */
void Bias::calcForceAndUpdateBias(awh_dvec biasForce,
                                  double *awhPotential, double *potential_jump,
                                  const gmx_multisim_t *ms,
                                  double t, gmx_int64_t step, int seed,
                                  FILE *fplog)
{
    if (step < 0)
    {
        gmx_fatal(FARGS, "The step number is negative which is not supported by the AWH code.");
    }

    const BiasParams    &params = params_;

    std::vector<double> *probWeightNeighbor = &tempWorkSpace_;

    /* If the convolved force is needed or this is a sampling step,
     * the bias in the current neighborhood needs to be up-to-date
     * and the probablity weights need to be calculated.
     */
    const bool sampleCoord   = doAtStep(params.numStepsSampleCoord, step);
    double     convolvedBias = 0;
    if (params.convolveForce || sampleCoord)
    {
        if (params.skipUpdates())
        {
            state_.doSkippedUpdatesInNeighborhood(params, grid());
        }

        convolvedBias = state_.updateProbabilityWeightsAndConvolvedBias(dimParams_, grid(), probWeightNeighbor);

        if (step > 0 && sampleCoord)
        {
            state_.sampleCoordAndPmf(grid(), *probWeightNeighbor, convolvedBias);
        }
    }

    /* Set the bias force and get the potential contribution from this bias.
     * The potential jump occurs at different times depending on how
     * the force is applied (and how the potential is normalized).
     * For the convolved force it happens when the bias is updated,
     * for the umbrella when the umbrella is moved.
     */
    double potential, newPotential;
    if (params.convolveForce)
    {
        state_.calcConvolvedForce(dimParams_, grid(), *probWeightNeighbor,
                                  biasForce);

        potential    = -convolvedBias*params.invBeta;
        newPotential = potential; /* Assume no jump */
    }
    else
    {
        /* Umbrella force */
        GMX_RELEASE_ASSERT(state_.points()[state_.refGridpoint()].inTargetRegion(),
                           "AWH bias grid point reference value is outside of the target region.");
        potential =
            state_.calcUmbrellaForceAndPotential(dimParams_, grid(), state_.refGridpoint(), biasForce);

        /* Moving the umbrella results in a force correction and
         * a new potential. The umbrella center is sampled as often as
         * the coordinate so we know the probability weights needed
         * for moving the umbrella are up-to-date.
         */
        if (sampleCoord)
        {
            newPotential = state_.moveUmbrella(dimParams_, grid(), *probWeightNeighbor, biasForce, step, seed, params.biasIndex);
        }
        else
        {
            newPotential = potential;
        }
    }

    /* Add this bias potential to the total sum. */
    *awhPotential += potential;

    /* Update the free energy estimates and bias and other history dependent method parameters */
    const int stepIntervalUpdateFreeEnergy = params.numSamplesUpdateFreeEnergy*params.numStepsSampleCoord;
    if (step > 0 && doAtStep(stepIntervalUpdateFreeEnergy, step))
    {
        state_.updateFreeEnergyAndAddSamplesToHistogram(dimParams_, grid(),
                                                        params_,
                                                        ms, t, step, fplog,
                                                        &updateList_);

        if (params.convolveForce)
        {
            /* The update results in a potential jump, so we need the new convolved potential. */
            newPotential = -calcConvolvedBias(dimParams_, grid(), state_.points(), state_.coordValue())*params.invBeta;
        }
    }

    /* Add the potential jump of this bias to the total sum. */
    *potential_jump += (newPotential - potential);

    /* Check the sampled histograms and potentially warn user if something is suspicious */
    checkHistograms(t, step, fplog);
}

/* Update the coordinate value with coordValue. */
void Bias::setCoordValue(int dim, double coordValue)
{
    GMX_RELEASE_ASSERT(dim >= 0 && dim < ndim(), "The dimension should be in range");

    state_.setCoordValue(grid(), dim, coordValue);
}

/* Restore the bias state from history on the master rank and broadcast it. */
void Bias::restoreStateFromHistory(const AwhBiasHistory &biasHistory,
                                   const t_commrec      *cr)
{
    if (MASTER(cr))
    {
        state_.restoreFromHistory(biasHistory, grid());
    }

    if (PAR(cr))
    {
        state_.broadcast(cr);
    }
}

/*! \brief
 * Get the number of domains.
 *
 * \param[in] awhBiasParams  Bias parameters.
 * \returns the number of domains.
 */
static int getNumDomains(const awh_bias_params_t &awhBiasParams)
{
    int numDomains = 1;
    for (int d = 0; d < awhBiasParams.ndim; d++)
    {
        numDomains *= awhBiasParams.dim_params[d].ninterval;
    }
    return numDomains;
}

/*! \brief
 * Get the domain id that the given point maps to.
 *
 * Each point maps to the domain that it is closest to the center of (for non-overlapping domains).
 *
 * \param[in] awhBiasParams  Bias parameters.
 * \param[in] grid           Grid.
 * \param[in] point          Point index.
 * \returns domain id.
 */
static int getDomainId(const awh_bias_params_t &awhBiasParams, const Grid &grid, int point)
{
    const awh_dim_params_t *awh_dim_params = awhBiasParams.dim_params;

    /* First get the multidimensional id. */
    awh_ivec domainIdDim;

    for (int d = 0; d < grid.ndim(); d++)
    {
        domainIdDim[d] = static_cast<int>(std::floor(1.*grid.point(point).index[d]/grid.axis(d).numPoints()*awh_dim_params[d].ninterval));
    }

    /* Convert to a linear domain id for returning. */
    awh_ivec numIntervalDim;

    for (int d = 0; d < awhBiasParams.ndim; d++)
    {
        numIntervalDim[d] = awh_dim_params[d].ninterval;
    }

    return multidim_array_index_to_linear(domainIdDim, grid.ndim(), numIntervalDim);
}


/*! \brief
 * Get the multidim boundary corner points for the given domain id.
 *
 * \param[in] awhBiasParams  Bias parameters.
 * \param[in] grid           Grid.
 * \param[in] domainId       Domain id.
 * \param[in,out] pointMin   Mininmum boundary point index.
 * \param[in,out] pointMax   Maximum boundary point index.
 */
static void getDomainBoundary(const awh_bias_params_t &awhBiasParams, const Grid &grid,
                              int domainId,
                              int *pointMin, int *pointMax)
{
    const awh_dim_params_t *awh_dim_params = awhBiasParams.dim_params;

    /* Convert the given linear domain id to a multidimensional one. */
    awh_ivec numInterval, domainIdDim = {0};

    for (int d = 0; d < awhBiasParams.ndim; d++)
    {
        numInterval[d] = awh_dim_params[d].ninterval;
    }
    linear_array_index_to_multidim(domainId, grid.ndim(), numInterval, domainIdDim);

    /* Get the multidimensional min and max coordinates that defines the domain */
    awh_ivec         numPoints, domainIMin, domainIMax;

    for (int d = 0; d < grid.ndim(); d++)
    {
        const double nonOverlap       = 1 - awh_dim_params[d].interval_overlap;
        const double intervalFraction = 1./(nonOverlap*(awh_dim_params[d].ninterval - 1) + 1);
        const double intervalSize     = grid.axis(d).numPoints()*intervalFraction;

        numPoints[d] = grid.axis(d).numPoints();

        /* The domain end points are given by the domain id scaled by the interval size and a non-overlap factor.
           For overlap > 0, i.e. non-overlap < 1, the interval is extended in both directions. */
        domainIMin[d] = static_cast<int>(std::floor(intervalSize*nonOverlap*domainIdDim[d]));

        /* Note that the max is just the min + the interval size. */
        domainIMax[d] = static_cast<int>(std::floor(intervalSize*(nonOverlap*domainIdDim[d] + 1)));
        domainIMax[d] = std::min(domainIMax[d], numPoints[d] - 1);
    }

    /* Finally, convert the multidimensional point indices to linear ones */
    *pointMin = multidim_array_index_to_linear(domainIMin, grid.ndim(), numPoints);
    *pointMax = multidim_array_index_to_linear(domainIMax, grid.ndim(), numPoints);
}

/*! \brief
 * Get the min and max boundary points of the rectangular domain that the given point maps to.
 *
 * \param[in] awhBiasParams  Bias parameters.
 * \param[in] grid           Grid.
 * \param[in] point          Point index.
 * \param[in,out] point_min  Mininmum boundary point index.
 * \param[in,out] point_max  Maximum boundary point index.
 */
static void getDomainBoundaryForPoint(const awh_bias_params_t &awhBiasParams, const Grid &grid, int point,
                                      int *point_min, int *point_max)
{
    getDomainBoundary(awhBiasParams, grid, getDomainId(awhBiasParams, grid, point),
                      point_min, point_max);

    /* Make sure that point is inside of the domain */
    for (int d = 0; d < grid.ndim(); d++)
    {
        int index_d = grid.point(point).index[d];
        GMX_RELEASE_ASSERT(grid.point(*point_min).index[d] <= index_d &&
                           grid.point(*point_max).index[d] >= index_d,
                           "AWH coord point outside of domain while partitioning.");
    }
}

/*! \brief
 * Prin information about domain partitioning.
 *
 * \param[in] awhPrefix      Prefix string.
 * \param[in] bias           The AWH bias.
 * \param[in] awhBiasParams  AWH bias parameters.
 * \param[in,out] fplog      Log file.
 */
static void printPartitioningDomainInit(const char              *awhPrefix,
                                        const Bias              &bias,
                                        const awh_bias_params_t &awhBiasParams,
                                        FILE                    *fplog)
{
    int numDomains = getNumDomains(awhBiasParams);
    if (numDomains == 1)
    {
        return;
    }

    const Grid &grid  = bias.grid();

    int         my_id = getDomainId(awhBiasParams, grid, bias.state().gridpointIndex());

    if (fplog != nullptr)
    {
        /* Print partitioning info for this simulation. */
        int domain_imin, domain_imax;
        getDomainBoundary(awhBiasParams, grid, my_id, &domain_imin, &domain_imax);

        fprintf(fplog, "%s partitioned the AWH domain into %d subdomains. This sim has target domain: ",
                awhPrefix, numDomains);
        for (int d = 0; d < bias.ndim(); d++)
        {
            fprintf(fplog, "[%g, %g]",
                    bias.dimParams()[d].scaleInternalToUserInput(grid.point(domain_imin).coordValue[d]),
                    bias.dimParams()[d].scaleInternalToUserInput(grid.point(domain_imax).coordValue[d]));
            if (d < bias.ndim() - 1)
            {
                fprintf(fplog, ", ");
            }
        }
        fprintf(fplog, "\n");
    }
    else
    {
        /* Print partitioning info about partitioning stderr if there is no log file (at preprocessing). */
        char output[STRLEN];

        sprintf(output, "%s will partition the AWH domain into %d subdomains. "
                "The coordinate value to subdomain mapping ('-->') for each dimensions is "
                "approximately as follows. This simulation will be assigned target domain id %d.\n",
                awhPrefix, numDomains, my_id);
        fprintf(stderr, "%s", wrap_lines(output, linewidth, indent, FALSE));

        for (int id = 0; id < numDomains; id++)
        {
            /* Get the multidim id of this domain */
            awh_ivec ninterval_dim, id_dim;
            for (int d = 0; d < awhBiasParams.ndim; d++)
            {
                ninterval_dim[d] = awhBiasParams.dim_params[d].ninterval;
            }
            linear_array_index_to_multidim(id, grid.ndim(), ninterval_dim, id_dim);

            /* Get the boundary points of each domain and the coordinate values that map to it. */
            int domain_imin, domain_imax;
            getDomainBoundary(awhBiasParams, grid, id, &domain_imin, &domain_imax);

            fprintf(stderr, "\n");
            for (int d = 0; d < bias.ndim(); d++)
            {
                int coord_imin = static_cast<int>(std::floor(1.0*id_dim[d]/ninterval_dim[d]*(bias.state().points().size() - 1)));
                int coord_imax = static_cast<int>(std::floor(1.0*(id_dim[d] + 1)/ninterval_dim[d]*(bias.state().points().size() - 1)));

                fprintf(stderr, "[%.4g, %.4g]",
                        bias.dimParams()[d].scaleInternalToUserInput(grid.point(coord_imin).coordValue[d]),
                        bias.dimParams()[d].scaleInternalToUserInput(grid.point(coord_imax).coordValue[d]));

                if (d < bias.ndim() - 1)
                {
                    fprintf(stderr, " x ");
                }
            }
            fprintf(stderr, " --> domain id %d: ", id);
            for (int d = 0; d < bias.ndim(); d++)
            {
                fprintf(stderr, "[%.4g, %.4g]",
                        bias.dimParams()[d].scaleInternalToUserInput(grid.point(domain_imin).coordValue[d]),
                        bias.dimParams()[d].scaleInternalToUserInput(grid.point(domain_imax).coordValue[d]));
                if (d < bias.ndim() - 1)
                {
                    fprintf(stderr, " x ");
                }
            }
        }
        fprintf(stderr, "\n\n");
    }
}

/*! \brief
 * Print information about initialization to log file.
 *
 * \param[in] bias           The AWH bias.
 * \param[in] awhBiasParams  AWH bias parameters.
 * \param[in,out] fplog      Log file.
 */
static void printInitializationToLog(const Bias              &bias,
                                     const awh_bias_params_t &awhBiasParams,
                                     FILE                    *fplog)
{
    char           awhstr[STRLEN];
    sprintf(awhstr, "\nawh%d:", bias.params().biasIndex + 1);

    printPartitioningDomainInit(awhstr, bias, awhBiasParams, fplog);
}

/* Constructor. */
Bias::Bias(FILE                          *fplog,
           const t_commrec               *cr,
           int                            biasIndexInCollection,
           const awh_params_t            &awhParams,
           const awh_bias_params_t       &awhBiasParams,
           const std::vector<DimParams>  &dimParamsInit,
           double                         beta,
           double                         mdTimeStep,
           BiasParams::DisableUpdateSkips disableUpdateSkips) :
    dimParams_(dimParamsInit),
    grid_(new Grid(dimParamsInit, awhBiasParams.dim_params)),
    params_(awhParams, awhBiasParams, dimParams_, beta, mdTimeStep, disableUpdateSkips, cr, grid_->axis(), biasIndexInCollection),
    state_(awhBiasParams, params_.histSizeInitial, dimParams_, grid()),
    tempWorkSpace_(),
    numWarningsIssued_(0)
{
    /* For a global update updateList covers all points, so reserve that */
    updateList_.reserve(grid_->numPoints());

    state_.initGridPointState(awhBiasParams, dimParams_, grid(), params_, awhParams.nbias, (cr != nullptr ? cr->ms : nullptr));

    /* Partition the AWH domain if any of the dimensions is set to be divided into more than 1 interval. */
    if (getNumDomains(awhBiasParams) > 1)
    {
        /* Partitioning  modifies the target distribution and the weighthistogram outside of the target domain.
           The target domain is determined based on the current coordinate reference value. */
        int pointMin, pointMax;
        getDomainBoundaryForPoint(awhBiasParams, grid(), state_.gridpointIndex(),
                                  &pointMin, &pointMax);

        params_.histSizeInitial = state_.partitionDomain(grid(), pointMin, pointMax);
    }

    /* Translate the point cover diameter fraction to a radius in units of number of points. */
    for (int d = 0; d < ndim(); d++)
    {
        double coverRadius     = 0.5*awhBiasParams.dim_params[d].coverDiameter;
        double spacing         = grid_->axis(d).spacing();
        params_.coverRadius[d] = spacing > 0 ?  static_cast<int>(std::round(coverRadius/spacing)) : 0;
    }

    /* Print information about AWH variables that are set internally but might be of interest to the user. */
    if ((cr == nullptr) || (MASTER(cr)))
    {
        printInitializationToLog(*this, awhBiasParams, fplog);
    }
}
