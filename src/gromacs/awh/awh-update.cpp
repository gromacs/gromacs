/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015,2016, by the GROMACS development team, led by
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

#include <assert.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/awh/awh.h"
#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/utilities.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "grid.h"
#include "history.h"
#include "internal.h"
#include "math.h"
#include "types.h"


//! A value that can be passed to exp() with result 0, also with SIMD
static const double c_largeNegativeExponent = -10000.0;

//! The largest acceptable positive exponent for variables that are passed to exp().
static const double c_largePositiveExponent =  700.0;

/* Convert internal coordinate units to external, user coordinate units. */
double scaleInternalToUserInput(const awh_bias_t *awh_bias, int dimIndex, double value)
{
    return value/awh_bias->userCoordUnitsToInternal[dimIndex];
}

/* Convert external, user coordinate units to internal coordinate units. */
double scaleUserInputToInternal(const awh_bias_t *awh_bias, int dimIndex, double value)
{
    return value*awh_bias->userCoordUnitsToInternal[dimIndex];
}

/*! \brief
 * Query if the coordinate point is in the target region.
 *
 * \param[in] coordpoint   Coordinate point to query about.
 * \returns true if the point is in the target region.
 */
static bool in_target_region(const coordpoint_t *coordpoint)
{
    return coordpoint->target > 0;
}

/* Sets the given array with PMF values. */
void getPmf(const awh_bias_t *awh_bias, const gmx_multisim_t *ms, double *pmf)
{
    /* The PMF is just the negative of the log of the sampled PMF histogram.
       Points with zero target weight are ignored, they will mostly contain noise. */

    if (awh_bias->numSharedUpdate > 1)
    {
        /* Need to temporarily exponentiate the log weights to sum over simulations */
        for (int i = 0; i < awh_bias->npoints; i++)
        {
            pmf[i] = in_target_region(&awh_bias->coordpoint[i]) ? std::exp(awh_bias->coordpoint[i].log_pmfsum) : 0;
        }

        gmx_sumd_sim(awh_bias->npoints, pmf, ms);

        /* Take log again to get (non-normalized) PMF */
        for (int i = 0; i < awh_bias->npoints; i++)
        {
            pmf[i] = in_target_region(&awh_bias->coordpoint[i]) ? -std::log(pmf[i]) : GMX_DOUBLE_MAX;
        }
    }
    else
    {
        for (int i = 0; i < awh_bias->npoints; i++)
        {
            pmf[i] = in_target_region(&awh_bias->coordpoint[i]) ? -awh_bias->coordpoint[i].log_pmfsum : GMX_DOUBLE_MAX;
        }
    }
}

/*! \brief
 * Find the minimum free energy value.
 *
 * \param[in] awh_bias  AWH bias.
 * \returns the minimum free energy value.
 */
static double freeEnergyMinimumValue(const awh_bias_t *awh_bias)
{
    double        f_min = GMX_DOUBLE_MAX;

    for (int m = 0; m < awh_bias->npoints; m++)
    {
        if (in_target_region(&awh_bias->coordpoint[m]) && awh_bias->coordpoint[m].free_energy < f_min)
        {
            f_min = awh_bias->coordpoint[m].free_energy;
        }
    }

    return f_min;
}

/*! \brief
 * Find the probability weight of a point given a coordinate value.
 *
 * The unnormalized weight is given by
 * w(point|value) = exp(bias(point) - U(value,point)),
 * where U is an harmonic umbrella potential.
 *
 * \param[in] awh_bias      AWH bias.
 * \param[in] point_index   Point to evaluate probability weight for.
 * \param[in] pointBias     Bias for the point (as a log weight).
 * \param[in] value         Coordinate value.
 * \returns the biased probability weight.
 */
static double biasedWeightFromPoint(const awh_bias_t *awh_bias, int point_index, double pointBias, const awh_dvec value)
{
    double weight = 0;

    /* Only points in the target reigon have non-zero weight */
    if (in_target_region(&awh_bias->coordpoint[point_index]))
    {
        double log_weight = pointBias;

        /* Add potential for all parameter dimensions */
        for (int d = 0; d < awh_bias->ndim; d++)
        {
            double dev = get_deviation_from_point_along_gridaxis(awh_bias->grid, d, point_index, value[d]);
            log_weight -= 0.5*awh_bias->betak[d]*dev*dev;
        }

        weight = std::exp(log_weight);
    }

    return weight;
}

/* Convolves the given PMF using the given AWH bias. */
void getConvolvedPmf(const awh_bias_t *awh_bias, const gmx_multisim_t *ms, double *convolvedPmf)
{
    double *pmf;

    snew(pmf, awh_bias->npoints);

    /* Get the PMF to convolve. */
    getPmf(awh_bias, ms, pmf);

    grid_point_t      *gridpoints  = awh_bias->grid->point;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        double freeEnergyWeights = 0;
        for (int n = 0; n < gridpoints[m].nneighbors; n++)
        {
            int neighbor = gridpoints[m].neighbor[n];

            /* The negative PMF is a positive bias. */
            double biasNeighbor = -pmf[neighbor];

            /* Add the convolved PMF weights for the neighbors of this point.
               Note that this function only adds point within the target > 0 region.
               Sum weights, take the logarithm last to get the free energy. */
            freeEnergyWeights += biasedWeightFromPoint(awh_bias, neighbor, biasNeighbor,
                                                       gridpoints[m].value);
        }

        GMX_RELEASE_ASSERT(freeEnergyWeights > 0, "Attempting to do log(<= 0) in AWH convolved PMF calculation.");
        convolvedPmf[m] =  -std::log(freeEnergyWeights);
    }

    sfree(pmf);
}

/*! \brief
 * Updates the target distribution of the AWH bias.
 *
 * The target distribution is always updated for all points
 * at the same time.
 *
 * \param[in,out] awh_bias      AWH bias.
 */
static void update_target(awh_bias_t *awh_bias)
{
    coordpoint_t       *coordpoint  = awh_bias->coordpoint;
    double              f_max       = 0, sum_target = 0, inv_sum;

    if (awh_bias->eTarget == eawhtargetCUTOFF)
    {
        f_max = freeEnergyMinimumValue(awh_bias) + awh_bias->target_param;
    }

    for (int m = 0; m < awh_bias->npoints; m++)
    {
        if (awh_bias->eTarget == eawhtargetCONSTANT)
        {
            coordpoint[m].target = 1;
        }
        else if (awh_bias->eTarget == eawhtargetCUTOFF)
        {
            double df =  coordpoint[m].free_energy - f_max;
            coordpoint[m].target = df > 0 ? std::exp(-df) : 1;
        }
        else if (awh_bias->eTarget == eawhtargetBOLTZMANN)
        {
            coordpoint[m].target = std::exp(-awh_bias->target_param*coordpoint[m].free_energy);
        }
        else if (awh_bias->eTarget == eawhtargetLOCALBOLTZMANN)
        {
            coordpoint[m].target = coordpoint[m].weightsum_ref;
        }

        /* All target types can be modulated by a constant factor. */
        coordpoint[m].target *= coordpoint[m].target_constant_weight;

        sum_target += coordpoint[m].target;
    }

    /* Normalize to 1 */
    inv_sum = 1./sum_target;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        coordpoint[m].target *= inv_sum;
    }
}

/*! \brief
 * Updates the bias of a coordinate point.
 *
 * \param[in,out] coordpoint  A coordinate point.
 */
static void update_bias_point(coordpoint_t *coordpoint)
{
    GMX_RELEASE_ASSERT(coordpoint->target > 0, "AWH target distribution must be > 0 to calculate the point bias.");
    coordpoint->bias = coordpoint->free_energy + std::log(coordpoint->target);
}

/* Initialize the target and bias values for the given AWH bias. */
void initTargetAndBias(awh_bias_t *awh_bias)
{

    /* The local Boltzmann distribution is special because the target distribution is updated as a function of the reference weighthistogram. */
    GMX_ASSERT((awh_bias->eTarget != eawhtargetLOCALBOLTZMANN) ||
               (awh_bias->eTarget == eawhtargetLOCALBOLTZMANN && awh_bias->coordpoint[0].weightsum_ref != 0),
               "AWH reference weight histogram not initialized properly with local Boltzmann target distribution." );

    update_target(awh_bias);

    for (int m = 0; m < awh_bias->npoints; m++)
    {
        coordpoint_t *coordpoint = &awh_bias->coordpoint[m];
        if (in_target_region(coordpoint))
        {
            update_bias_point(coordpoint);
        }
        else
        {
            /* Note that for zero target this is a value that represents -infinity but should not be used for biasing. */
            coordpoint->bias = c_largeNegativeExponent;
        }
    }
}

/*! \brief
 * Query if something should be done at this step.
 *
 * \param[in] stepInterval   Step interval for doing something.
 * \param[in] step           Current step.
 * \returns true if something should be done at this step.
 */
static bool do_at_step(int stepInterval, gmx_int64_t step)
{
    GMX_ASSERT(stepInterval > 0, "All step intervals in AWH should be > 0");

    return (step % stepInterval == 0);
}

/*! \brief
 * Set pull forces for the pull coordinates AWH is responsible for.
 *
 * \param[in] awh          AWH struct.
 * \param[in,out] pull     Pull struct.
 * \param[in] mdatoms      Atom properties.
 * \param[in,out] force    Forces.
 * \param[in,out] force    Virial tensor.
 */
static void apply_bias_force_to_pull_coords(const awh_t        *awh,
                                            struct pull_t      *pull_work,
                                            const t_mdatoms    *mdatoms,
                                            rvec               *force,
                                            tensor              virial)
{
    for (int k = 0; k < awh->nbias; k++)
    {
        awh_bias_t *awh_bias = &awh->awh_bias[k];

        for (int d = 0; d < awh_bias->ndim; d++)
        {
            apply_external_pull_coord_force(pull_work, awh_bias->pull_coord_index[d],
                                            awh_bias->bias_force[d], mdatoms, force, virial);
        }
    }
}

/*! \brief
 * Puts together a string describing a grid point.
 *
 * \param[in] grid         The grid.
 * \param[in] point        Grid point index.
 * \param[in,out] pointstr String for the point.
 */
static void set_grid_point_value_string(const grid_t *grid, int point, char *pointstr)
{
    char buf[STRLEN];

    strcpy(pointstr, "(");

    for (int d = 0; d < grid->ndim; d++)
    {
        sprintf(buf, "%g", grid->point[point].value[d]);
        strcat(pointstr, buf);
        if (d < grid->ndim - 1)
        {
            sprintf(buf, ",");
            strcat(pointstr, buf);
        }
        else
        {
            sprintf(buf, ")");
            strcat(pointstr, buf);
        }
    }
}

/*! \brief
 * Makes checks for the collected data (histograms) and warns if issues are detected.
 *
 * \param[in] awh_bias     AWH bias.
 * \param[in] t            Time.
 * \param[in] step         Time step.
 * \param[in,out] fplog    Output file for warnings.
 */
static void checkOnData(const awh_bias_t *awh_bias, double t, gmx_int64_t step, FILE *fplog)
{
    const int          max_nwarnings          = 1;   /* The maximum number of warnings to print per check */
    const int          max_alltime_nwarnings  = 5;   /* The maximum number of warnings to print ever */
    static int         alltime_nwarnings      = 0;   /* The number of warnings printed ever */
    const double       max_histogram_ratio    = 0.5; /* Tolerance for printing a warning about the histogram ratios */
    const int          nstcheck_on_data       = std::max(awh_bias->nstupdate_free_energy,
                                                         (50000*awh_bias->nstsample_coord/awh_bias->nstupdate_free_energy)*awh_bias->nstupdate_free_energy);
    int                nwarnings;
    double             sum_W = 0, sum_visits = 0;
    double             inv_norm_visits, inv_norm_W;
    coordpoint_t      *coordpoint  = awh_bias->coordpoint;
    grid_t            *grid        = awh_bias->grid;

    if (fplog == NULL || alltime_nwarnings >= max_alltime_nwarnings || awh_bias->in_initial ||
        (step == 0) || (step % nstcheck_on_data != 0))
    {
        return;
    }

    /* Sum up the histograms and get their normalization */
    for (int m = 0; m < grid->npoints; m++)
    {
        if (in_target_region(&coordpoint[m]))
        {
            sum_visits += coordpoint[m].visits_tot;
            sum_W      += coordpoint[m].weightsum_tot;
        }
    }
    inv_norm_visits = 1./sum_visits;
    inv_norm_W      = 1./sum_W;

    /* Check all points for warnings */
    nwarnings = 0;
    for (int m = 0; m < grid->npoints; m++)
    {
        bool bSkip_point = FALSE;

        /* Skip points close to boundary or non-target region */
        for (int n = 0; (n < grid->point[m].nneighbors) && !bSkip_point; n++)
        {
            int ineighbor = grid->point[m].neighbor[n];
            bSkip_point = !in_target_region(&coordpoint[ineighbor]);
            for (int d = 0; (d < grid->ndim) && !bSkip_point; d++)
            {
                bSkip_point = (grid->point[ineighbor].index[d] == 0) ||  (grid->point[ineighbor].index[d] ==  grid->axis[d].npoints - 1);
            }
        }

        /* Warn if the coordinate distribution less than the target distribution with a certain fraction somewhere */
        if (!bSkip_point &&
            (coordpoint[m].weightsum_tot*inv_norm_W*max_histogram_ratio > coordpoint[m].visits_tot*inv_norm_visits ))
        {
            char pointvaluestr[STRLEN], warningmsg[STRLEN];

            set_grid_point_value_string(grid, m, pointvaluestr);
            sprintf(warningmsg, "\nawh%d warning: "
                    "at t = %g ps the obtained coordinate distribution at coordinate value %s "
                    "is less than a fraction %g of the reference distribution at that point. "
                    "If you are not certain about your settings you might want to increase your pull force constant or "
                    "modify your sampling region.\n",
                    awh_bias->biasIndex + 1, t, pointvaluestr, max_histogram_ratio);
            fprintf(fplog, "%s", wrap_lines(warningmsg, linewidth, indent, FALSE));

            nwarnings++;
            alltime_nwarnings++;
        }
        if (nwarnings >= max_nwarnings)
        {
            break;
        }
    }
    if (alltime_nwarnings >= max_alltime_nwarnings)
    {
        fprintf(fplog, "\nawh%d: suppressing future AWH warnings since the maximum number has .\n", awh_bias->biasIndex + 1);
    }
}

/*! \brief
 * Calculates and sets the force the coordinate experiences from an umbrella centered at the given point.
 *
 * The umbrella potential is an harmonic potential given by 0.5k(coord value - point value)^2. This
 * value is also returned.
 *
 * \param[in] awh_bias     AWH bias.
 * \param[in] k            Force constant for each dimension.
 * \param[in] point        Point for umbrella center.
 * \param[in,out] force    Force vector to set.
 * Returns the umbrella potential.
 */
static double calcUmbrellaForceAndPotential(const awh_bias_t *awh_bias, const awh_dvec k, int point, awh_dvec force)
{
    double potential = 0;
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        double dev = get_deviation_from_point_along_gridaxis(awh_bias->grid, d, point, awh_bias->coord_value[d]);

        /* Force from harmonic potential 0.5*k*dev^2 */
        force[d]   = -k[d]*dev;
        potential += 0.5*k[d]*dev*dev;
    }

    return potential;
}

/*! \brief
 * Calculates and sets the convolved force acting on the coordinate.
 *
 * The convolved force is the weighted sum of forces from umbrellas
 * located at each point in the grid.
 *
 * \param[in] awh_bias     AWH bias.
 * \param[in,out] force    Bias force vector to set.
 */
static void calc_convolved_force(const awh_bias_t *awh_bias, awh_dvec force)
{
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        force[d] = 0;
    }

    /* Only neighboring points have non-negligible contribution. */
    for (int n = 0; n < awh_bias->grid->point[awh_bias->coord_value_index].nneighbors; n++)
    {
        int      m_neighbor;
        double   weight_neighbor;
        awh_dvec force_from_neighbor;

        weight_neighbor = awh_bias->prob_weight_neighbor[n];
        m_neighbor      = awh_bias->grid->point[awh_bias->coord_value_index].neighbor[n];

        /* Get the umbrella force from this point. The returned potential is ignored here. */
        calcUmbrellaForceAndPotential(awh_bias, awh_bias->k, m_neighbor, force_from_neighbor);

        /* Add the weighted umbrella force to the convolved force. */
        for (int d = 0; d < awh_bias->ndim; d++)
        {
            force[d] += force_from_neighbor[d]*weight_neighbor;
        }
    }
}

/*! \brief
 * Sample a new reference point given the current coordinate value.
 *
 * It is assumed that the probability distribution has been updated.
 *
 * \param[in] awh_bias     AWH bias.
 * \param[in] step         Time step needed for the random number generator.
 * \param[in] seed         Random seed.
 * \returns the index of the sampled point.
 */
static int sampleReferenceCoordpoint(const awh_bias_t *awh_bias, gmx_int64_t step, int seed)
{
    /* Sample new reference value from the probability distribution which is defined for the neighboring
       points of the current coordinate value.*/

    /* In order to use the same seed for all AWH biases and get independent
       samples we use the index of the bias. */
    int n_sampled  = get_sample_from_distribution(awh_bias->prob_weight_neighbor, awh_bias->grid->point[awh_bias->coord_value_index].nneighbors,
                                                  step, seed, awh_bias->biasIndex);

    return awh_bias->grid->point[awh_bias->coord_value_index].neighbor[n_sampled];
}

/*! \brief
 * Move the center point of the umbrella potential.
 *
 * A new umbrella center is sampled from the biased distibution. Also, the bias
 * force is updated and the new potential is return.
 *
 * This function should only be called when the bias force is not being convolved.
 * It is assumed that the probability distribution has been updated.
 *
 * \param[in,out] awh_bias     AWH bias.
 * \param[in] step             Time step needed for the random number generator.
 * \param[in] seed             Random seed.
 * \returns the new potential value.
 */
double moveUmbrella(awh_bias_t *awh_bias, gmx_int64_t step, int seed)
{
    /* Generate and set a new coordinate reference value */
    awh_bias->refCoordpoint = sampleReferenceCoordpoint(awh_bias, step, seed);

    awh_dvec newForce;
    double   newPotential =
        calcUmbrellaForceAndPotential(awh_bias, awh_bias->k, awh_bias->refCoordpoint, newForce);

    /*  A modification of the reference value at time t will lead to a different
        force over t-dt/2 to t and over t to t+dt/2. For high switching rates
        this means the force and velocity will change signs roughly as often.
        To avoid any issues we take the average of the previous and new force
        at steps when the reference value has been moved. E.g. if the ref. value
        is set every step to (coord dvalue +/- delta) would give zero force.
     */
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        /* Average of the current and new force */
        awh_bias->bias_force[d] = 0.5*(awh_bias->bias_force[d] + newForce[d]);
    }

    return newPotential;
}

/*! \brief
 * Sets the histogram rescaling factors needed to control the histogram size.
 *
 * For sake of robustness, the reference weight histogram can grow at a rate
 * different from the actual sampling rate. Typically this happens for a limited
 * initial time, alternatively growth is scaled down by a constant factor for all
 * times. Since the size of the reference histogram sets the size of the free
 * energy update this should be reflected also in the PMF. Thus the PMF histogram
 * needs to be rescaled too.
 *
 * This function should only be called by update() or wrapped by a function that
 * knows what scale factors should be applied when, e.g,
 * setSkippedUpdateHistogramScaleFactors().
 *
 * \param[in,out] awh_bias             AWH bias.
 * \param[in] histsizeNew              New reference weight histogram size.
 * \param[in] histsizeOld              Previous reference weight histogram size (before adding new samples).
 * \param[in,out] weighthistScaling    Scaling factor for the reference weight histogram.
 * \param[in,out] logPmfsumScaling     Log of the scaling factor for the PMF histogram.
 */
static void setHistogramUpdateScaleFactors(const awh_bias_t *awh_bias, double histsizeNew, double histsizeOld,
                                           double *weighthistScaling, double *logPmfsumScaling)
{

    /* The two scaling factors below are slightly different (ignoring the log factor) because the
       reference and the PMF histogram apply weight scaling differently. The weight histogram
       applies is  locally, i.e. each sample is scaled down meaning all samples get equal weight.
       It is done this way because that is what target type local Boltzmann (for which
       target = weight histogram) needs. In contrast, the PMF histogram is rescaled globally
       by repeatedly scaling down the whole histogram. The reasons for doing it this way are:
       1) empirically this is necessary for converging the PMF; 2) since the extraction of
       the PMF is theoretically only valid for a constant bias, new samples should get more
       weight than old ones for which the bias is fluctuating more. */
    *weighthistScaling = histsizeNew/(histsizeOld + awh_bias->update_weight*awh_bias->localWeightScaling);
    *logPmfsumScaling  = std::log(histsizeNew/(histsizeOld + awh_bias->update_weight));
}

/*! \brief
 * Check if the AWH constant parameters are such that we attempt to skip updates.
 *
 * Generally, we can skip (postpone) updates of points that are non-local at the
 * time of the update if we for later times, when the points with skipped updates
 * have become local, know exactly how to apply the previous updates. The free
 * energy updates only depend on local sampling, but the histogram rescaling factors
 * generally depend on the histogram size (all samples). If the histogram size is
 * kept constant or the scaling factors are trivial, this is not a problem. However,
 * if the histogram growth is scaled down by some factor the size at the time of the
 * update needs to be known. It would be fairly simple to, for a deterministically
 * growing histogram, backtrack and calculate this value, but currently we just disallow
 * this case. This is not a restriction because it only affects the local Boltzmann
 * target type for which every update is currently anyway global because the target
 * is always updated globally.
 *
 * \param[in] awh_bias             AWH bias.
 */
bool canSkipUpdates(const awh_bias_t *awh_bias)
{
    return (awh_bias->localWeightScaling == 1);
}

/*! \brief
 * Sets the histogram rescaling factors needed for skipped updates.
 *
 * \param[in,out] awh_bias             AWH bias.
 * \param[in,out] weighthistScaling    Scaling factor for the reference weight histogram.
 * \param[in,out] logPmfsumScaling     Log of the scaling factor for the PMF histogram.
 */
static void setSkippedUpdateHistogramScaleFactors(const awh_bias_t *awh_bias, double *weighthistScaling, double *logPmfsumScaling)
{
    GMX_ASSERT(canSkipUpdates(awh_bias), "Calling function for skipped updates when skipping updates is not allowed");

    if (awh_bias->in_initial)
    {
        /* In between global updates the reference histogram size is kept constant so we trivially know what the
            histogram size was at the time of the skipped update. */
        setHistogramUpdateScaleFactors(awh_bias, awh_bias->histsize, awh_bias->histsize,
                                       weighthistScaling, logPmfsumScaling);
    }
    else
    {
        /* In the final stage, the reference histogram grows at the sampling rate which gives trivial scale factors. */
        *weighthistScaling   = 1;
        *logPmfsumScaling    = 0;
    }
}

/*! \brief
 * Update the free energy estimate of a point.
 *
 * The free energy update here is inherently local, i.e. it just depends on local sampling and on constant
 * AWH parameters. This assumes that the target distribution is a constant, or at least constant in between
 * global updates.
 *
 * \param[in,out] coordpoint    Coordinate point to update free energy for.
 * \param[in] awh_bias          AWH bias.
 * \param[in] weight_point      Sampled probability weight at this point.
 */
static void update_free_energy_point(coordpoint_t *coordpoint,  const awh_bias_t *awh_bias, double weight_point)
{
    double        weighthist_sampled, weighthist_target, df;

    weighthist_sampled              = coordpoint->weightsum_ref + weight_point;
    weighthist_target               = coordpoint->weightsum_ref + awh_bias->update_weight*coordpoint->target;
    df                              = -std::log(weighthist_sampled/weighthist_target);
    coordpoint->free_energy        += df;

    GMX_RELEASE_ASSERT(coordpoint->free_energy < c_largePositiveExponent,
                       "Very large free energy differences or badly normalized free energy in AWH update.");
    /* TODO: potentially we need to add a normalization of the free energy here. Non-local points
       (with sampling weight = 0) keep accumulating a positive update ~ log(1 + 1/histsize). The free energy
       also sets the normalization of the bias. We can run into problems with very large numbers since we
       need to be able to expoentiate the bias. Usually, if histsize increases fast enough this does not
       become a problem. But in the extreme case, histsize = constant and then failure will surely happen at
       some point. In the initial stage this can happen if for one reason or another the sampling interval is
       not covered (fast enough). Not that adding a normalization here, e.g. subtracting log(1 - 1/histsize),
       makes the free energy update explicitly dependent on histsize. When skipping updates the correct value
       would need to be provided. */
}

/*! \brief
 * Update the reference weight histogram of a point.
 *
 * \param[in,out] coordpoint    Coordinate point to update.
 * \param[in] awh_bias          AWH bias.
 * \param[in] weight_point      Sampled probability weight at this point.
 * \param[in] scaleFactor       Factor to rescale the histogram with.
 */
static void update_weighthistogram_point(coordpoint_t *coordpoint, const awh_bias_t *awh_bias, double weight_point, double scaleFactor)
{
    if (awh_bias->idealWeighthistUpdate)
    {
        /* Grow histogram using the target distribution. */
        coordpoint->weightsum_ref += coordpoint->target*awh_bias->update_weight*awh_bias->localWeightScaling;
    }
    else
    {
        /* Grow using the actual samples (which are distributed ~ as target). */
        coordpoint->weightsum_ref += weight_point*awh_bias->localWeightScaling;
    }

    coordpoint->weightsum_ref *= scaleFactor;
}

/*! \brief
 * Apply a point update.
 *
 * This updates local properties that can be updated without accessing or affecting
 * all points. This excludes updating the size of reference weight histogram and
 * the target distribution. The bias update is excluded only because if updates
 * have been skipped this function will be called multiple times, while the bias
 * only needs to be updated once (last).
 *
 * Since this function only performs the update with the given arguments and does
 * not know anything about the time of the update, the last update index is not
 * updated here. The caller should take care of updating the update index.
 *
 * \param[in,out] coordpoint       Coordinate point to update.
 * \param[in] awh_bias             AWH bias.
 * \param[in] weight_point         Sampled probability weight at this point.
 * \param[in] weighthistScaling    Scaling factor for the reference weight histogram.
 * \param[in] logPmfsumScaling     Log of the scaling factor for the PMF histogram.
 */
static void update_point(coordpoint_t *coordpoint, const awh_bias_t *awh_bias,
                         double weight_point, double weighthistScaling, double logPmfsumScaling)
{
    update_free_energy_point(coordpoint, awh_bias, weight_point);
    update_weighthistogram_point(coordpoint, awh_bias, weight_point, weighthistScaling);
    coordpoint->log_pmfsum += logPmfsumScaling;
}

/*! \brief
 * Apply previous updates that were skipped.
 *
 * An update can only be skipped if the parameters needed for the update are constant or
 * deterministic so that the same update can be performed at a later time.
 * Here, the necessary parameters are the sampled weight and scaling factors for the
 * histograms. The scaling factors are provided as arguments only to avoid recalculating
 * them for each point
 *
 * The last update index is also updated here.
 *
 * \param[in,out] coordpoint    Coordinate point to update.
 * \param[in] awh_bias          AWH bias.
 * \param[in] step              Current time step.
 * \param[in] weighthistScaling  Scale factor for the reference weight histogram.
 * \param[in] logPmfsumScaling      Scale factor for the reference PMF histogram.
 * \return true if at least one update was applied.
 */
static bool update_point_skipped(coordpoint_t *coordpoint, const awh_bias_t *awh_bias,
                                 gmx_int64_t step, double weighthistScaling, double logPmfsumScaling)
{
    int last_update_index, nupdates_skipped;

    GMX_ASSERT(canSkipUpdates(awh_bias), "Calling function for skipped updates when skipping updates is not allowed");

    if (!in_target_region(coordpoint))
    {
        return false;
    }

    /* The most current past update */
    last_update_index = (step - 1)/awh_bias->nstupdate_free_energy;
    nupdates_skipped  = last_update_index - coordpoint->last_update_index;

    if (nupdates_skipped == 0)
    {
        /* Was not updated */
        return false;
    }

    for (int i = 0; i < nupdates_skipped; i++)
    {
        /* This point was non-local at the time of the update meaning no weight */
        update_point(coordpoint, awh_bias, 0, weighthistScaling, logPmfsumScaling);
    }

    /* Only past updates are applied here. */
    coordpoint->last_update_index = last_update_index;

    /* Was updated */
    return true;
}

/* TODO: should probably merge do_skipped_updates_for_all_points and
   do_skipped_updates_in_neighborhood into one function that takes as
   input a vector with points to be updated. */

/*! \brief
 * Do all previously skipped updates.
 *
 * \param[in,out] awh_bias          AWH bias.
 * \param[in] step                  Current time step.
 */
static void do_skipped_updates_for_all_points(awh_bias_t *awh_bias, gmx_int64_t step)
{
    double weighthistScaling, logPmfsumScaling;

    setSkippedUpdateHistogramScaleFactors(awh_bias, &weighthistScaling, &logPmfsumScaling);

    for (int m = 0; m < awh_bias->npoints; m++)
    {
        bool bUpdated =  update_point_skipped(&awh_bias->coordpoint[m], awh_bias, step, weighthistScaling, logPmfsumScaling);

        /* Update the bias for this point only if there were skipped updates in the past to avoid calculating the log unneccessarily */
        if (bUpdated)
        {
            update_bias_point(&awh_bias->coordpoint[m]);
        }
    }
}

/*! \brief
 * Do previously skipped updates in this neighborhood.
 *
 * \param[in,out] awh_bias          AWH bias.
 * \param[in] step                  Current time step.
 */
static void do_skipped_updates_in_neighborhood(awh_bias_t *awh_bias, gmx_int64_t step)
{
    double weighthistScaling, logPmfsumScaling;
    int    nneighbors;
    int   *neighbor;

    setSkippedUpdateHistogramScaleFactors(awh_bias, &weighthistScaling, &logPmfsumScaling);

    /* For each neighbor point of the center point, refresh its state by adding the results of all past, skipped updates. */
    nneighbors = awh_bias->grid->point[awh_bias->coord_value_index].nneighbors;
    neighbor   = awh_bias->grid->point[awh_bias->coord_value_index].neighbor;
    for (int n = 0; n < nneighbors; n++)
    {
        bool bUpdated = update_point_skipped(&awh_bias->coordpoint[neighbor[n]], awh_bias, step, weighthistScaling, logPmfsumScaling);

        if (bUpdated)
        {
            update_bias_point(&awh_bias->coordpoint[neighbor[n]]);
        }
    }
}

/*! \brief
 * Merge update lists from multiple sharing simulations.
 *
 * \param[in,out] updateList        Update list for this simulation (assumed >= npoints long).
 * \param[in,out] nupdate_ptr       Number of points in the update list.
 * \param[in] npoints               Total number of points.
 * \param[in] ms                    Struct for multi-simulation communication.
 */
static void merge_shared_update_lists(int *updateList, int *nupdate_ptr, int npoints, const gmx_multisim_t *ms)
{
    int *nupdate_point;
    int  nupdate_merged;

    /* Flag the update points of this sim.
       TODO: we can probably avoid allocating this array and just use the input array. */
    snew(nupdate_point, npoints);
    for (int i = 0; i < *nupdate_ptr; i++)
    {
        nupdate_point[updateList[i]] = 1;
    }

    /* Sum over the sims to get all the flagged points (snew initializes array elements to 0) */
    gmx_sumi_sim(npoints, nupdate_point, ms);

    /* Collect the indices of the flagged points in place. The resulting array will be the merged update list.*/
    nupdate_merged = 0;
    for (int m = 0; m < npoints; m++)
    {
        if (nupdate_point[m] > 0)
        {
            updateList[nupdate_merged++] = m;
        }
    }

    sfree(nupdate_point);

    *nupdate_ptr = nupdate_merged;
}

/*! \brief
 * Generate an update list of points sampled since the last update.
 *
 * \param[in] awh_bias              AWH bias.
 * \param[in] ms                    Struct for multi-simulation communication.
 * \param[in,out] updateList        Local update list to set (assumed >= npoints long).
 * \param[in,out] nupdate_ptr       Number of points in the update list.
 */
static void make_local_update_list(const awh_bias_t *awh_bias, const gmx_multisim_t *ms, int *updateList, int *nupdate_ptr)
{
    int      point_index, nupdate;
    awh_ivec origin_dim, npoints_dim;
    bool     pointExists;

    /* Define the update search grid */
    for (int d = 0; d < awh_bias->grid->ndim; d++)
    {
        origin_dim[d]  = awh_bias->origin_updatelist[d];
        npoints_dim[d] = awh_bias->end_updatelist[d] - awh_bias->origin_updatelist[d] + 1;

        /* Because the end_updatelist is unwrapped it can be > (npoints - 1) so that npoints_dim can be > npoints in grid.
           This helps for calculating the distance/number of points but should be removed and fixed when the way of
           updating origin/end updatelist is changed (see sampleProbabilityWeights). */
        npoints_dim[d] = std::min(awh_bias->grid->axis[d].npoints, npoints_dim[d]);
    }

    /* Make the update list */
    point_index   = -1;
    nupdate       = 0;
    pointExists   = true;
    while (pointExists)
    {
        pointExists = get_next_point_in_local_grid(awh_bias->grid, origin_dim, npoints_dim, &point_index);

        if (pointExists && in_target_region(&awh_bias->coordpoint[point_index]))
        {
            updateList[nupdate] = point_index;
            nupdate++;
        }
    }

    if (awh_bias->numSharedUpdate > 1)
    {
        merge_shared_update_lists(updateList, &nupdate, awh_bias->npoints, ms);
    }

    /* Return the number of points to update */
    *nupdate_ptr = nupdate;
}

/*! \brief
 * Reset the range used to make the local update list.
 *
 * \param[in,out] awh_bias              AWH bias.
 */
static void reset_local_update_range(awh_bias_t *awh_bias)
{
    for (int d = 0; d < awh_bias->grid->ndim; d++)
    {
        /* This gives the  minimum range consisting only of the current closest point. */
        awh_bias->origin_updatelist[d] = awh_bias->grid->point[awh_bias->coord_value_index].index[d];
        awh_bias->end_updatelist[d]    = awh_bias->grid->point[awh_bias->coord_value_index].index[d];
    }
}

/*! \brief
 * Apply a point update with new sampling data.
 *
 * The last update index is also updated here.
 *
 * \param[in,out] coordpoint       Coordinate point to update.
 * \param[in] awh_bias             AWH bias.
 * \param[in] step                 Time step for updating the last update index.
 * \param[in] weighthistScaling   Scaling factor for the reference weight histogram.
 * \param[in] logPmfsumScaling     Log of the scaling factor for the PMF histogram.
 */
static void update_point_new(coordpoint_t *coordpoint, const awh_bias_t *awh_bias,
                             gmx_int64_t step, double weighthistScaling, double logPmfsumScaling)
{
    update_point(coordpoint, awh_bias, coordpoint->weightsum_iteration, weighthistScaling, logPmfsumScaling);

    /* Remember that this was most recent update performed on this point. */
    coordpoint->last_update_index   = step/awh_bias->nstupdate_free_energy;
}


/*! \brief
 * Add partial histograms (accumulating between updates) to accumulating histograms.
 *
 * \param[in] awh_bias                 AWH bias.
 * \param[in] ms                       Struct for multi-simulation communication.
 * \param[in] local_update_list        List of points with data.
 * \param[in] nlocal_update            Number of points in the list.
 */
static void sum_histograms(const awh_bias_t *awh_bias, const gmx_multisim_t *ms, const int *local_update_list, int nlocal_update)
{
    /* Sum histograms over multiple simulations if needed. */
    if (awh_bias->numSharedUpdate > 1)
    {
        /* Collect the weights and counts in linear arrays to be able to use gmx_sumd_sim. */
        double               *weight_distr, *coord_visits;

        snew(weight_distr, nlocal_update);
        snew(coord_visits, nlocal_update);

        for (int ilocal = 0; ilocal < nlocal_update; ilocal++)
        {
            int iglobal = local_update_list[ilocal];

            weight_distr[ilocal] = awh_bias->coordpoint[iglobal].weightsum_iteration;
            coord_visits[ilocal] = awh_bias->coordpoint[iglobal].visits_iteration;
        }

        gmx_sumd_sim(nlocal_update, weight_distr, ms);
        gmx_sumd_sim(nlocal_update, coord_visits, ms);

        /* Transfer back the result */
        for (int ilocal = 0; ilocal < nlocal_update; ilocal++)
        {
            int iglobal = local_update_list[ilocal];

            awh_bias->coordpoint[iglobal].weightsum_iteration = weight_distr[ilocal];
            awh_bias->coordpoint[iglobal].visits_iteration    = coord_visits[ilocal];
        }

        sfree(coord_visits);
        sfree(weight_distr);
    }

    /* Now add the partial counts and weights to the accumulating histograms.
       Note: we still need to use the weights for the update so we wait
       with resetting them until the end of the update. */
    for (int ilocal = 0; ilocal < nlocal_update; ilocal++)
    {
        int    iglobal = local_update_list[ilocal];

        double weights = awh_bias->coordpoint[iglobal].weightsum_iteration;
        double visits  = awh_bias->coordpoint[iglobal].visits_iteration;
        awh_bias->coordpoint[iglobal].weightsum_covering += weights;
        awh_bias->coordpoint[iglobal].weightsum_tot      += weights;
        awh_bias->coordpoint[iglobal].visits_tot         += visits;
    }
}

/*! \brief
 * Returns the new size of the reference weight histogram in the initial stage.
 *
 * This function also takes care resetting the histogram used for covering checks
 * and for exiting the initial stage.
 *
 * \param[in,out] awh_bias    AWH bias.
 * \param[in] t               Time.
 * \param[in] bCovered        True if the sampling interval has been covered enough.
 * \param[in,out] fplog       Log file.
 * \returns the new histogram size.
 */
static double newHistsizeInitialStage(awh_bias_t *awh_bias, double t, bool bCovered, FILE *fplog)
{
    /* The histogram size is kept constant until the sampling region has been covered
       and the the current sample weight is large enough. */
    if (!bCovered ||
        (awh_bias->scaledSampleWeight < awh_bias->maxScaledSampleWeight))
    {
        return awh_bias->histsize;
    }

    /* Reset the covering weight histogram. If we got this far we are either entering a
       new covering stage with a new covering histogram or exiting the initial stage
       altogether. */
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        awh_bias->coordpoint[m].weightsum_covering = 0.;
    }

    /*  The current sample weigth is now the maximum. */
    double prevMaxScaledSampleWeight = awh_bias->maxScaledSampleWeight;
    awh_bias->maxScaledSampleWeight = awh_bias->scaledSampleWeight;

    /* Increase the histogram size by a constant scale factor if we can, i.e. if the sample weight
       resulting from such a scaling is still larger than the previous maximum sample weight
       (ensuring that the sample weights at the end of each covering stage are monotonically
       increasing). If we cannot, exit the initial stage without changing the histogram size. */

    /* The scale factor. The value is not very critical but should obviously be > 1 (or the exit
       will happen very late) and probably < 5 or so (or there will be no initial stage). */
    static const double scaleFactor = 3;

    bool                bExit        = (awh_bias->scaledSampleWeight - std::log(scaleFactor) <= prevMaxScaledSampleWeight);
    double              histsize_new = bExit ? awh_bias->histsize : awh_bias->histsize*scaleFactor;

    /* Update the AWH bias about the exit. */
    awh_bias->in_initial = !bExit;

    /* Print information about coverings and if there was an exit. */
    if (fplog != NULL)
    {
        char                buf[STRLEN];
        sprintf(buf, "\nawh%d:", awh_bias->biasIndex + 1);
        fprintf(fplog, "%s covering at t = %g ps. Decreased the update size.\n", buf, t);

        if (bExit)
        {
            fprintf(fplog, "%s out of the initial stage at t = %g.\n", buf, t);
            /* It would be nice to have a way of estimating a minimum time until exit but it
               is difficult because the exit time is determined by how long it takes to cover
               relative to the time it takes to "regaining" enough sample weight. The latter
               is easy to calculate, but how the former depends on the histogram size
               is not known. */
        }
        fflush(fplog);
    }
    return histsize_new;
}

/*! \brief
 * Check if the sampling region has been covered "enough" or not.
 *
 * A one-dimensional interval is defined as covered if each point has
 * accumulated the same weight as is in the peak of a discretized normal
 * distribution. For multiple dimensions, the weights are simply projected
 * onto each dimension and the multidimensional space is covered if each
 * dimension is.
 *
 * \note The covering criterion for multiple dimensions could improved, e.g.
 * by using a path finding algorithm.
 *
 * \param[in] awh_bias        AWH bias.
 * \returns true if covered.
 */
static bool isCovered(const awh_bias_t *awh_bias)
{
    int           dimindex;
    double        weight_peak, f_max;
    bool        **bCovered, **bCheck;
    bool          bAll_covered;
    grid_t       *grid = awh_bias->grid;

    /* Allocate arrays: one for checking that each coordinate value along each dimension as been covered,
       and one for keeping track of which points to check. */
    snew(bCovered, grid->ndim);
    snew(bCheck, grid->ndim);
    for (int d = 0; d < grid->ndim; d++)
    {
        snew(bCovered[d], grid->axis[d].npoints);
        snew(bCheck[d], grid->axis[d].npoints);
    }

    for (int d = 0; d < grid->ndim; d++)
    {
        for (int n = 0; n < grid->axis[d].npoints; n++)
        {
            bCovered[d][n] = FALSE;
            bCheck[d][n]   = FALSE;
        }
    }

    /* Get free energy cutoff if there is one. Points above the cutoff should be ignored. */
    f_max = GMX_DOUBLE_MAX;

    if (awh_bias->eTarget == eawhtargetCUTOFF)
    {
        f_max = freeEnergyMinimumValue(awh_bias) + awh_bias->target_param;
    }

    /* Covering weight cutoff: weight_peak = weight at peak from mean assuming Gaussian transition prob. distr. */
    weight_peak = 1;
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        weight_peak *= grid->axis[d].spacing*std::sqrt(awh_bias->betak[d]*0.5*M_1_PI);
    }

    /* Project the covering data of all points onto each dimension */
    for (int m = 0; m < grid->npoints; m++)
    {
        coordpoint_t *coordpoint = &awh_bias->coordpoint[m];

        for (int d = 0; d < grid->ndim; d++)
        {
            int n = grid->point[m].index[d];

            /* Is covered if it was already covered or if there is enough weight at the current point */
            bCovered[d][n] = bCovered[d][n] || coordpoint->weightsum_covering > weight_peak;

            /* Check for covering if there is at least point in this slice that is in the target region and within the cutoff */
            bCheck[d][n] = bCheck[d][n] || (in_target_region(coordpoint) && coordpoint->free_energy < f_max);
        }
    }

    /* Check for global covering, i.e. all points covered. Break if not covered somewhere. */
    bAll_covered        = TRUE;
    dimindex            = 0;
    while (bAll_covered && dimindex < grid->ndim)
    {
        int n = 0;
        while (bAll_covered && n < grid->axis[dimindex].npoints)
        {
            bAll_covered = !bCheck[dimindex][n] || bCovered[dimindex][n];
            n++;
        }
        dimindex++;
    }

    for (int d = 0; d < grid->ndim; d++)
    {
        sfree(bCovered[d]);
    }
    sfree(bCovered);

    return bAll_covered;
}

/*! \brief
 * Return the new reference weight histogram size for the current update.
 *
 * This function also takes care of checking for covering in the initial stage.
 *
 * \param[in,out] awh_bias       AWH bias.
 * \param[in] t                  Time.
 * \param[in] step               Time step.
 * \param[in,out] bCovered_ptr   Set to true if the sampling interval has been covered enough.
 * \param[in,out] fplog          Log file.
 * \returns the new histogram size.
 */
static double newHistsize(awh_bias_t *awh_bias, double t, gmx_int64_t step, bool *bCovered_ptr, FILE *fplog)
{
    double    histsize_new;
    bool      bCovered = FALSE;

    if (awh_bias->in_initial)
    {
        int nstcheck_covered;

        /* In the initial stage, the histogram grows dynamically as a function of the number of coverings. */
        nstcheck_covered = std::max((500*awh_bias->nstsample_coord/awh_bias->nstupdate_free_energy)*awh_bias->nstupdate_free_energy,
                                    awh_bias->nstupdate_free_energy);
        bCovered         = (step % nstcheck_covered == 0) && isCovered(awh_bias);
        histsize_new     = newHistsizeInitialStage(awh_bias, t, bCovered, fplog);
    }
    else
    {
        /* If not in the initial stage, the histogram grows at a linear, possibly scaled down, rate. */
        histsize_new = awh_bias->histsize + awh_bias->update_weight*awh_bias->localWeightScaling;
    }

    *bCovered_ptr = bCovered;
    return histsize_new;
}

/*! \brief
 * Performs an AWH update.
 *
 * The objective of the update is to use collected samples (probability weights)
 * to improve the free energy estimate. For sake of efficiency, the update is
 * local whenever possible meaning that only points that have actually been sampled
 * are accessed and updated here. For certain AWH settings or at certain steps
 * however, global need to be performed. Besides the actual free energy update, this
 * function takes care of ensuring future convergence of the free energy. Convergence
 * is obtained by increasing the size of the reference weight histogram in a controlled
 * (sometimes dynamic) manner. Also, there are AWH variables that are direct functions
 * of the free energy or sampling history that need to be updated here, namely the target
 * distribution and the bias function.
 *
 * \param[in,out] awh_bias       AWH bias.
 * \param[in] ms                 Struct for multi-simulation communication.
 * \param[in] t                  Time.
 * \param[in step                Time step.
 * \param[in,out] fplog           Log file.
 */
static void update(awh_bias_t *awh_bias, const gmx_multisim_t *ms, double t, gmx_int64_t step, FILE *fplog)
{
    /* This array is only used in this scope and is always re-initialized */
    int     *updateList = awh_bias->updateList;

    /* Make a list of all local points, i.e. those that could have been touched since
       the last update. These are the points needed for summing histograms below
       (non-local points only add zeros). For local updates, this will also be the
       final update list. */
    int      nlocal_update;
    make_local_update_list(awh_bias, ms, updateList, &nlocal_update);

    /* Reset the range for the next update */
    reset_local_update_range(awh_bias);

    /* Add samples to histograms for all local points and sync simulations if needed */
    sum_histograms(awh_bias, ms, updateList, nlocal_update);

    /* Update target distribution? */
    bool updateTarget = (awh_bias->eTarget != eawhtargetCONSTANT &&
                         do_at_step(awh_bias->nstupdate_target, step));

    /* The weighthistogram size after this update. */
    bool   haveCovered;
    double histsize_new = newHistsize(awh_bias, t, step, &haveCovered, fplog);

    /* Make the update list. Usually we try to only update local points but if the update has
       non-trivial or non-deterministic affect on non-local points a global update is needed.
       This is the case when: 1) a covering occurred in the initial stage, leading to non-trivial
       histogram rescaling factors; or 2) the target distribution will be updated since we don't
       make any assumption on its form; or 3) the AWH parameters are such that we never attempt
       to skip non-local updates. */
    int nupdate;
    if (updateTarget || haveCovered || !canSkipUpdates(awh_bias))
    {
        /* Global update, just add all points. */
        nupdate = 0;
        for (int m = 0; m < awh_bias->npoints; m++)
        {
            if (in_target_region(&awh_bias->coordpoint[m]))
            {
                updateList[nupdate] = m;
                nupdate++;
            }
        }
    }
    else
    {
        /* Local update. The update list equals the local update list of touched points (which
           has already been set). */
        nupdate = nlocal_update;
    }

    /* Set histogram scale factors. */
    double   weighthistScalingSkipped, weighthistScalingNew, logPmfsumScalingSkipped, logPmfsumScalingNew;
    if (canSkipUpdates(awh_bias))
    {
        setSkippedUpdateHistogramScaleFactors(awh_bias, &weighthistScalingSkipped, &logPmfsumScalingSkipped);
    }
    setHistogramUpdateScaleFactors(awh_bias, histsize_new, awh_bias->histsize,
                                   &weighthistScalingNew, &logPmfsumScalingNew);

    /* Update free energy and reference weight histogram for points in the update list. */
    for (int iupdate = 0; iupdate < nupdate; iupdate++)
    {
        coordpoint_t *coordpoint_to_update = &awh_bias->coordpoint[updateList[iupdate]];

        /* Do updates from previous update steps that were skipped because this point was at that time non-local. */
        if (canSkipUpdates(awh_bias))
        {
            update_point_skipped(coordpoint_to_update, awh_bias, step, weighthistScalingSkipped, logPmfsumScalingSkipped);
        }

        /* Now do an update with new sampling data. */
        update_point_new(coordpoint_to_update, awh_bias, step, weighthistScalingNew, logPmfsumScalingNew);

        /* Reset histograms collected for this update. */
        coordpoint_to_update->weightsum_iteration = 0;
        coordpoint_to_update->visits_iteration    = 0;
    }

    /* Only update the histogram size after we are done with the local point updates */
    awh_bias->histsize = histsize_new;

    /* The weight of new samples relative to previous ones change when the histogram is
       rescaled. We keep the log since this number can become very large. */
    awh_bias->scaledSampleWeight -= std::log(weighthistScalingNew);

    if (updateTarget)
    {
        /* The target distribution is always updated for all points at once. */
        update_target(awh_bias);
    }

    /* Update the bias. The bias is updated separately and last since it simply a function of
       the free energy and the target distribution and we want to avoid doing extra work. */
    for (int iupdate = 0; iupdate < nupdate; iupdate++)
    {
        update_bias_point(&awh_bias->coordpoint[updateList[iupdate]]);
    }
}

/*! \brief
 * Update the probability weights and the convolved bias.
 *
 * Given a coordinate value, each coordinate point is assigned a probability
 * weight, w(point|value), that depends on the current bias function. The sum
 * of these weights is needed for normalizing the probability sum to 1 but
 * also equals the effective, or convolved, biasing weight for this coordinate
 * value. The convolved bias is needed e.g. for extracting the PMF, so we save
 * it here since this saves us from doing extra exponential function evaluations
 * later on.
 *
 * \param[in,out] awh_bias       AWH bias.
 */
static void updateProbabilityWeightsAndConvolvedBias(awh_bias_t *awh_bias)
{
    double      weight_sum, inv_weight_sum;
    double     *weight = awh_bias->prob_weight_neighbor;
    grid_t     *grid   = awh_bias->grid;

    /* Sum of probability weights */
    weight_sum = 0;

    /* Only neighbors of the current coordinate value will have a non-negligible chance of getting sampled */
    for (int n = 0; n < grid->point[awh_bias->coord_value_index].nneighbors; n++)
    {
        int neighbor =  grid->point[awh_bias->coord_value_index].neighbor[n];
        weight[n]   = biasedWeightFromPoint(awh_bias, neighbor, awh_bias->coordpoint[neighbor].bias,
                                            awh_bias->coord_value);
        weight_sum += weight[n];
    }
    GMX_RELEASE_ASSERT(weight_sum > 0, "zero probability weight when updating AWH probability weights.");
    inv_weight_sum = 1./weight_sum;

    /* Normalize probabilities to 1 */
    for (int n = 0; n < grid->point[awh_bias->coord_value_index].nneighbors; n++)
    {
        weight[n] *= inv_weight_sum;
    }

    awh_bias->convolved_bias = std::log(weight_sum);
}

/* Calculates the convolved bias for the given coordinate value. */
double calc_convolved_bias(const awh_bias_t *awh_bias, const awh_dvec coord_value)
{
    grid_point_t      *gridpoints  = awh_bias->grid->point;
    coordpoint_t      *coordpoints = awh_bias->coordpoint;
    int                point       = get_closest_index_in_grid(coord_value, awh_bias->grid);
    double             weight_sum  = 0;

    /* Sum the probability weights from the neighborhood of the given point */
    for (int n = 0; n < gridpoints[point].nneighbors; n++)
    {
        int neighbor =  gridpoints[point].neighbor[n];
        weight_sum += biasedWeightFromPoint(awh_bias, neighbor, coordpoints[neighbor].bias,
                                            coord_value);
    }

    /* Returns -GMX_DOUBLE_MAX if no neighboring points where in the target region. */
    return (weight_sum > 0) ? std::log(weight_sum) : -GMX_DOUBLE_MAX;
}

/*! \brief
 * Save the current probability weights for future updates and analysis.
 *
 * Points in the current neighborhood will now have data meaning they
 * need to be included in the local update list of the next update.
 * Therefore, the local update range is also update here.
 *
 * \param[in,out] awh_bias       AWH bias.
 */
static void sampleProbabilityWeights(awh_bias_t *awh_bias)
{
    int m_neighbor_origin, m_neighbor_end;
    int nneighbors = awh_bias->grid->point[awh_bias->coord_value_index].nneighbors;

    /* Save weights for next update */
    for (int n = 0; n < nneighbors; n++)
    {
        int m_neighbor = awh_bias->grid->point[awh_bias->coord_value_index].neighbor[n];
        awh_bias->coordpoint[m_neighbor].weightsum_iteration += awh_bias->prob_weight_neighbor[n];
    }

    /* Update the local update range. Two corner points define this rectangular domain. We need to
       choose two new corner points such that the new domain contains both the old update range and
       the current neighborhood. In the simplest case when an update is performed every sample, the
        update range would simply equal the current neighborhood. */
    m_neighbor_origin = awh_bias->grid->point[awh_bias->coord_value_index].neighbor[0];
    m_neighbor_end    = awh_bias->grid->point[awh_bias->coord_value_index].neighbor[nneighbors - 1];
    for (int d = 0; d < awh_bias->grid->ndim; d++)
    {
        int origin_d, end_d;

        origin_d = awh_bias->grid->point[m_neighbor_origin].index[d];
        end_d    = awh_bias->grid->point[m_neighbor_end].index[d];

        if (origin_d > end_d)
        {
            /* Unwrap if wrapped around the boundary (only happens for periodic boundaries).
               This has been already for the stored index interval. */

            /* TODO: what we want to do is to find the smallest the update interval that contains all points that need to be updated.
               This amounts to combining two intervals, the current [origin, end] update interval and the new touched neighborhood
               into a new interval that contains all points from both the old intervals.

               For periodic boundaries it becomes slightly more complicated than for closed boundaries because the it needs not be
               true that origin < end (so one can't simply relate the origin/end in the min()/max() below). The strategy here is to choose the
               origin closest to a reference point (index 0) and then unwrap the end index if needed and choose the largest end index.
               This ensures that both intervals are in the new interval but it's not necessarily the smallest.
               I can't think of a better way of solving this than going through each possibility and checking them.
             */
            end_d += awh_bias->grid->axis[d].npoints_period;
        }

        awh_bias->origin_updatelist[d] = std::min(awh_bias->origin_updatelist[d], origin_d);
        awh_bias->end_updatelist[d]    = std::max(awh_bias->end_updatelist[d], end_d);
    }
}

/*! \brief
 * Update the PMF histogram with the current coordinate value.
 *
 * Because of the finite width of the harmonic potential, the free energy
 * defined for each coordinate point does not exactly equal that of the
 * actual coordinate, the PMF. However, the PMF can be estimated by applying
 * the relation exp(-PMF) = exp(-bias_convolved)*P_biased/Z, i.e. by keeping a
 * reweighted histogram of the coordinate value. Strictly, this relies on
 * the unknown normalization constant Z being either constant or known. Here,
 * neither is true except in the long simulation time limit. Empirically however,
 * it works (mainly because how the PMF histogram is rescaled).
 *
 * \param[in,out] awh_bias       AWH bias.
 */
static void sample_pmf(awh_bias_t *awh_bias)
{
    bool               bIn_array;
    coordpoint_t      *coordpoint;

    /* Only save coordinate data that is in range (the given index is always
       in range even if the coordinate value is not). */
    bIn_array   = value_is_in_grid(awh_bias->coord_value, awh_bias->grid);
    coordpoint  = &awh_bias->coordpoint[awh_bias->coord_value_index];

    /* Save PMF sum and keep a histogram of the sampled coordinate values */
    if (bIn_array  && in_target_region(coordpoint))
    {
        coordpoint->log_pmfsum        = expsum(coordpoint->log_pmfsum, -awh_bias->convolved_bias);
        coordpoint->visits_iteration += 1;
    }
}

/*! \brief
 * Sample observables for future updates or analysis.
 *
 * \param[in,out] awh_bias       AWH bias.
 */
static void do_sampling(awh_bias_t *awh_bias)
{
    /* Sampling-based deconvolution extracting the PMF */
    sample_pmf(awh_bias);

    /* Save probability weights for the update */
    sampleProbabilityWeights(awh_bias);
}

/*! \brief
 * Do what the AWH bias needs to do every step.
 *
 * The main objectives of the AWH bias step is to 1) set the bias force
 * and potential; 2) update the free energy and bias; and 3) sample
 * observables (e.g. the PMF).
 *
 * \param[in,out] awh_bias       AWH bias.
 * \param[in] convolveForce      True if the bias force is convolved, otherwise apply single umbrella force.
 * \param[in,out] awhPotential   Bias potential from all AWH biases.
 * \param[in,out] potentialJump  Change in bias potential for all AWH biases.
 * \param[in] ms                 Struct for multi-simulation communication.
 * \param[in] t                  Time.
 * \param[in] step               Time step.
 * \param[in] seed               Random seed.
 * \param[in,out] fplog          Log file.
 */
static void do_awh_bias_step(awh_bias_t *awh_bias,
                             bool convolveForce,
                             double *awhPotential, double *potential_jump,
                             const gmx_multisim_t *ms, double t, gmx_int64_t step, int seed, FILE *fplog)
{
    /* The grid point closest to the coordinate value defines the current neighborhood of points. Besides
       at steps when global updates and/or checks are performed, only the neighborhood will be touched. */
    awh_bias->coord_value_index = get_closest_index_in_grid(awh_bias->coord_value, awh_bias->grid);

    /* If the convolved force is needed or this is a sampling step, the bias in the current neighborhood
       needs to be up-to-date and the probablity weights need to be calculated. */
    if (convolveForce || do_at_step(awh_bias->nstsample_coord, step))
    {
        if (canSkipUpdates(awh_bias))
        {
            do_skipped_updates_in_neighborhood(awh_bias, step);
        }

        updateProbabilityWeightsAndConvolvedBias(awh_bias);

        if (do_at_step(awh_bias->nstsample_coord, step))
        {
            do_sampling(awh_bias);
        }
    }

    /* Set the bias force and get the potential contribution from this bias.
       The potential jump occurs at different times depending on how the force is applied
       (and how the potential is normalized). For the convolved force it happens when
       the bias is updated, for the umbrella when the umbrella is moved. */
    double potential, newPotential;
    if (convolveForce)
    {
        calc_convolved_force(awh_bias, awh_bias->bias_force);
        potential    =  -awh_bias->convolved_bias*awh_bias->invBeta;
        newPotential = potential; /* Assume no jump */
    }
    else
    {
        /* Umbrella force */
        GMX_RELEASE_ASSERT(in_target_region(&awh_bias->coordpoint[awh_bias->refCoordpoint]),
                           "AWH bias coordinate reference value is outside of the target region.");
        potential =
            calcUmbrellaForceAndPotential(awh_bias, awh_bias->k, awh_bias->refCoordpoint, awh_bias->bias_force);

        /* Moving the umbrella results in a force correction and a new potential. The umbrella center is sampled as
           often as the coordinate so we know the probability weights needed for moving the umbrella are up-to-date. */
        if (do_at_step(awh_bias->nstsample_coord, step))
        {
            newPotential = moveUmbrella(awh_bias, step, seed);
        }
        else
        {
            newPotential = potential;
        }
    }

    /* Add this bias potential to the total sum. */
    *awhPotential += potential;

    /* Update the free energy estimates and bias and other history dependent method parameters */
    if (do_at_step(awh_bias->nstupdate_free_energy, step))
    {
        update(awh_bias, ms, t, step, fplog);

        if (convolveForce)
        {
            /* The update results in a potential jump, so we need the new convolved potential. */
            newPotential = -calc_convolved_bias(awh_bias, awh_bias->coord_value)*awh_bias->invBeta;
        }
    }

    /* Add the the potential jump of this bias to the total sum. */
    *potential_jump += (newPotential - potential);

    /* Check the sampled data (histograms) and potentially warn user if something is suspicious */
    checkOnData(awh_bias, t, step, fplog);
}

/*! \brief
 * Do what AWH needs to do every step.
 *
 * The main objectives of the AWH step is to sum up the total bias potential and
 * to make sure each bias does what it needs to do every step (set forces, sample,
 * do updates; see do_awh_bias_step). The bias potential offset, resulting from
 * the accumulated change of all biases, is also updated here.
 *
 * \param[in,out] awh            AWH struct.
 * \param[in,out] potential      Potential from all AWH biases.
 * \param[in] ms                 Struct for multi-simulation communication.
 * \param[in] t                  Time.
 * \param[in] step               Time step.
 * \param[in,out] fplog          Log file.
 */
static void do_awh_step(awh_t                *awh,
                        double               *potential,
                        const gmx_multisim_t *ms,
                        double                t,
                        gmx_int64_t           step,
                        FILE                 *fplog)
{
    /* During the AWH step the potential can instantaneously jump due to either
       an bias update or moving the umbrella. The jumps are kept track of and
       subtracted from the potential in order to get a useful conserved energy quantity. */
    double potentialJump = 0;
    *potential      = awh->potential_offset;

    for (int k = 0; k < awh->nbias; k++)
    {
        do_awh_bias_step(&awh->awh_bias[k],
                         awh->convolveForce,
                         potential, &potentialJump,
                         ms, t, step, awh->seed, fplog);
    }

    /* Keep track of the total potential shift needed to remove the potential jumps. */
    awh->potential_offset -= potentialJump;
}

/*! \brief
 * Update the coordinate value with the corresponding pull coordinate value.
 *
 * \param[in,out] awh_bias       AWH bias.
 * \param[in, out] pull_work     Pull struct.
 * \param[in] pbc                Struct with information about periodic boundary conditions.
 */
static void setAwhBiasCoordValue(awh_bias_t *awh_bias, struct pull_t *pull_work, const t_pbc *pbc)
{
    /* Keep own copy of current coordinate value. */
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        get_pull_coord_value(pull_work, awh_bias->pull_coord_index[d], pbc, &awh_bias->coord_value[d]);
    }
}

/*! \brief
 * Update the coordinate values for all AWH biases.
 *
 * \param[in,out] awh            AWH struct.
 * \param[in, out] pull_work     Pull struct
 * \param[in] ePBC               Type of periodic boundary conditions
 * \param[in] box                Box vectors.
 */
static void setAwhCoordValues(awh_t *awh, struct pull_t *pull_work,
                              int ePBC, const matrix box)
{
    t_pbc    pbc;

    set_pbc(&pbc, ePBC, box);

    for (int k = 0; k < awh->nbias; k++)
    {
        setAwhBiasCoordValue(&awh->awh_bias[k], pull_work, &pbc);
    }
}

/* Do an AWH biasing update. */
real update_awh(awh_t                  *awh,
                struct pull_t          *pull_work,
                int                     ePBC,
                const t_mdatoms        *mdatoms,
                const matrix            box,
                rvec                   *force,
                tensor                  virial,
                const gmx_multisim_t   *ms,
                double                  t,
                gmx_int64_t             step,
                struct gmx_wallcycle   *wallcycle,
                FILE                   *fplog)
{
    double   bias_potential = 0;

    wallcycle_start(wallcycle, ewcAWH);

    /* Update the AWH coordinate values with those of the corresponding pull coordinates. */
    setAwhCoordValues(awh, pull_work, ePBC, box);

    /* Perform an AWH biasing step: this means, at regular intervals, sampling observables
       based on the input pull coordinate value, setting the bias force and/or updating the AWH bias state. */
    do_awh_step(awh, &bias_potential, ms, t, step, fplog);

    /* Communicate the bias force to the pull struct.
     * The bias potential is returned at the end of this function,
     * so that it can be added externally to the correct energy data block.
     */
    apply_bias_force_to_pull_coords(awh, pull_work, mdatoms, force, virial);

    wallcycle_stop(wallcycle, ewcAWH);

    return static_cast<real>(bias_potential);
}
