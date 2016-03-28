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

#include "data-writer.h"
#include "grid.h"
#include "history.h"
#include "internal.h"
#include "math.h"
#include "types.h"


//! A value that can be passed to exp() with result 0, also with SIMD
static const double c_largeNegativeExponent = -10000.0;

static bool in_target_region(const coordpoint_t *coordpoint)
{
    return coordpoint->target > 0;
}

bool write_point_to_output(const awh_bias_t *awh_bias, int point_index)
{
    return in_target_region(&awh_bias->coordpoint[point_index]);
}

void getPmf(const awh_bias_t *awh_bias, const gmx_multisim_t *ms, float *pmf)
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

        gmx_sumf_sim(awh_bias->npoints, pmf, ms);

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

/* Find max frelative to global min */
static double get_f_min(const awh_bias_t *awh_bias)
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
static double get_biased_weight_from_point(const awh_bias_t *awh_bias, int point_index, double bias, const awh_dvec value)
{
    double weight = 0;

    /* Only points in the target reigon have non-zero weight */
    if (in_target_region(&awh_bias->coordpoint[point_index]))
    {
        double log_weight = bias;

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

void getConvolvedPmf(const awh_bias_t *awh_bias, const gmx_multisim_t *ms, float *convolvedPmf)
{
    float *pmf;

    snew(pmf, awh_bias->npoints);

    /* Get the PMF to convolve. */
    getPmf(awh_bias, ms, pmf);

    grid_point_t      *gridpoints  = awh_bias->grid->point;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        double freeEnergyWeights = 0;
        for (int n = 0; n < gridpoints->nneighbors; n++)
        {
            int neighbor = gridpoints[m].neighbor[n];

            /* The negative PMF is a positive bias. */
            double biasNeighbor = -pmf[neighbor];

            /* Add the convolved PMF weights for the neighbors of this point.
               Note that this function only adds point within the target > 0 region.
               Sum weights, take the logarithm last to get the free energy. */
            freeEnergyWeights += get_biased_weight_from_point(awh_bias, neighbor, biasNeighbor,
                                                              gridpoints[m].value);
        }

        GMX_RELEASE_ASSERT(freeEnergyWeights > 0, "Attempting to do log(< 0) in AWH convolved PMF calculation.");
        convolvedPmf[m] =  -std::log(freeEnergyWeights);
    }

    sfree(pmf);
}

void update_target(awh_bias_t *awh_bias)
{
    coordpoint_t       *coordpoint  = awh_bias->coordpoint;
    double              f_max       = 0, sum_target = 0, inv_sum;

    if (awh_bias->eTarget == eawhtargetCUTOFF)
    {
        f_max = get_f_min(awh_bias) + awh_bias->target_param;
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
        else if (awh_bias->eTarget == eawhtargetWEIGHTHIST)
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

static void update_bias_point(coordpoint_t *coordpoint)
{
    GMX_RELEASE_ASSERT(coordpoint->target > 0, "AWH target distribution must be > 0 to calculate the point bias.");
    coordpoint->bias = coordpoint->free_energy + std::log(coordpoint->target);
}

void initTargetAndBias(awh_bias_t *awh_bias)
{

    /* The local Boltzmann distribution is special because the target distribution is updated as a function of the reference weighthistogram. */
    GMX_ASSERT((awh_bias->eTarget != eawhtargetWEIGHTHIST) ||
               (awh_bias->eTarget == eawhtargetWEIGHTHIST && awh_bias->coordpoint[0].weightsum_ref != 0),
               "AWH reference weight histogram not initialized properly with local Boltzmann target distribution." );

    update_target(awh_bias);

    coordpoint_t *coordpoint = awh_bias->coordpoint;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        /* Note that for zero target this is a value that represents -infinity but should not be used for biasing. */
        coordpoint->bias = coordpoint->target > 0 ? coordpoint->free_energy + std::log(coordpoint->target) : c_largeNegativeExponent;
    }
}

static bool do_at_step(int stepInterval, gmx_int64_t step)
{
    GMX_ASSERT(stepInterval > 0, "All step intervals in AWH should be > 0");

    return (step % stepInterval == 0);
}

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

static void set_grid_point_value_string(const grid_t *grid, int ipoint, char *coordstr)
{
    char buf[STRLEN];

    strcpy(coordstr, "(");

    for (int d = 0; d < grid->ndim; d++)
    {
        sprintf(buf, "%g", grid->point[ipoint].value[d]);
        strcat(coordstr, buf);
        if (d < grid->ndim - 1)
        {
            sprintf(buf, ",");
            strcat(coordstr, buf);
        }
        else
        {
            sprintf(buf, ")");
            strcat(coordstr, buf);
        }
    }
}

static void check_on_data(const awh_bias_t *awh_bias, int awh_id, double t, FILE *fplog, gmx_int64_t step)
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
                    awh_id + 1, t, pointvaluestr, max_histogram_ratio);
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
        fprintf(fplog, "\nawh%d: suppressing future AWH warnings since the maximum number has .\n", awh_id + 1);
    }
}

/* Calculates and sets the force on the coordinate from a single grid point.
 * Returns the potential.
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

/* Calculates and sets the convolved force on the coordinate,
   i.e. the sum of force contributions from all grid points. */
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

        /* Add the force data of this point */
        calcUmbrellaForceAndPotential(awh_bias, awh_bias->k, m_neighbor, force_from_neighbor);

        for (int d = 0; d < awh_bias->ndim; d++)
        {
            force[d] += force_from_neighbor[d]*weight_neighbor;
        }
    }
}

static int sample_coord_refvalue_index(const awh_bias_t *awh_bias, gmx_int64_t step, int seed)
{
    /* Sample new reference value from the probability distribution which is defined for the neighboring
       points of the current coordinate value. */

    /* In order to use the same seed for all AWH biases and get independent
     * samples we use the index of the bias.
     */
    int n_sampled  = get_sample_from_distribution(awh_bias->prob_weight_neighbor, awh_bias->grid->point[awh_bias->coord_value_index].nneighbors,
                                                  step, seed, awh_bias->biasIndex);

    return awh_bias->grid->point[awh_bias->coord_value_index].neighbor[n_sampled];
}

double moveUmbrella(awh_bias_t *awh_bias, gmx_int64_t step, int seed)
{
    /* Generate and set a new coordinate reference value */
    awh_bias->coord_refvalue_index =
        sample_coord_refvalue_index(awh_bias, step, seed);

    awh_dvec newForce;
    double   newPotential =
        calcUmbrellaForceAndPotential(awh_bias, awh_bias->k, awh_bias->coord_refvalue_index, newForce);

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

static void get_histogram_scaling_skipped_update(const awh_bias_t *awh_bias, double *scale_weighthist_ptr, double *scale_pmfsum_ptr)
{
    if (awh_bias->in_initial || awh_bias->eGrowth == eawhgrowthNONE)
    {
        /* These scaling factors account for keeping the histogram size constant. */
        *scale_weighthist_ptr = awh_bias->histsize/(awh_bias->histsize + awh_bias->update_weight*awh_bias->weight_scaling);
        *scale_pmfsum_ptr     = std::log(awh_bias->histsize/(awh_bias->histsize + awh_bias->update_weight));
    }
    else
    {
        /* Linear growth */
        *scale_weighthist_ptr     = 1;
        *scale_pmfsum_ptr         = 0;
    }
}

static void get_histogram_scaling_new_update(const awh_bias_t *awh_bias, double histsize_new, double histsize_old,
                                             bool bCovered, double *scale_weighthist_ptr, double *scale_pmfsum_ptr)
{
    if (bCovered)
    {
        /* The histsize can change non-deterministically if we just covered in the intial stage */
        *scale_weighthist_ptr = histsize_new/(histsize_old + awh_bias->update_weight*awh_bias->weight_scaling);
        *scale_pmfsum_ptr     = std::log(histsize_new/(histsize_old + awh_bias->update_weight));
    }
    else
    {
        /* For all other cases, the scaling factor is constant. */
        get_histogram_scaling_skipped_update(awh_bias, scale_weighthist_ptr, scale_pmfsum_ptr);
    }
}

static void scale_pmfsum_point(coordpoint_t *coordpoint, double log_scale)
{
    coordpoint->log_pmfsum += log_scale;
}

static void update_free_energy_point(coordpoint_t *coordpoint,  const awh_bias_t *awh_bias, double weight_point)
{
    double        weighthist_sampled, weighthist_target, df;

    weighthist_sampled              = coordpoint->weightsum_ref + weight_point;
    weighthist_target               = coordpoint->weightsum_ref + awh_bias->update_weight*coordpoint->target;
    df                              = -std::log(weighthist_sampled/weighthist_target);
    coordpoint->free_energy        += df;

    /* Note: potentially it could be useful to renormalize f since this determines the bias normalization
       which in turn is used for the extracting the PMF. */
}

static void update_weighthistogram_point(coordpoint_t *coordpoint, const awh_bias_t *awh_bias, double weight_point, double rescale)
{
    if (awh_bias->idealWeighthistUpdate)
    {
        /* Grow histogram using the target distribution. */
        coordpoint->weightsum_ref += coordpoint->target*awh_bias->update_weight*awh_bias->weight_scaling;
    }
    else
    {
        /* Grow using the actual samples (which are distributed ~ as target). */
        coordpoint->weightsum_ref += weight_point*awh_bias->weight_scaling;
    }

    coordpoint->weightsum_ref *= rescale;
}

static void update_point(coordpoint_t *coordpoint, const awh_bias_t *awh_bias,
                         double weight_point, double scale_weighthist, double scale_pmfsum)
{
    update_free_energy_point(coordpoint, awh_bias, weight_point);
    update_weighthistogram_point(coordpoint, awh_bias, weight_point, scale_weighthist);
    scale_pmfsum_point(coordpoint, scale_pmfsum);
}

/* The previously post-poned non-local updates */
static bool update_point_skipped(coordpoint_t *coordpoint, const awh_bias_t *awh_bias,
                                 gmx_int64_t step, double scale_weighthist, double scale_pmfsum)
{
    int last_update_index, nupdates_skipped;

    if (!in_target_region(coordpoint))
    {
        return FALSE;
    }

    /* The most current past update */
    last_update_index = (step - 1)/awh_bias->nstupdate_free_energy;
    nupdates_skipped  = last_update_index - coordpoint->last_update_index;

    if (nupdates_skipped == 0)
    {
        /* Was not updated */
        return FALSE;
    }

    for (int i = 0; i < nupdates_skipped; i++)
    {
        /* This point was non-local at the time of the update meaning no weight */
        update_point(coordpoint, awh_bias, 0, scale_weighthist, scale_pmfsum);
    }

    coordpoint->last_update_index = last_update_index;

    /* Was updated */
    return TRUE;
}

/* Make sure the bias is up-to-date.
 * Note: this only deals with past updates, i.e. not including the new update at an step.
 * This means that if we want the point to be up-to-date we need to call this function _before_ the update.
 * After the update, the point may again be outdated if the update was non-global. The only way to get around
 * from this limitation is to either always do global updates (no skipped updates) or for each step flag when
 *  an AWH update has been performed, but currently there's not much point to adding more complexity.
 */
static void do_skipped_updates_for_all_points(awh_bias_t *awh_bias, gmx_int64_t step)
{
    double scale_weighthist, scale_pmfsum;

    get_histogram_scaling_skipped_update(awh_bias, &scale_weighthist, &scale_pmfsum);

    for (int m = 0; m < awh_bias->npoints; m++)
    {
        bool bUpdated =  update_point_skipped(&awh_bias->coordpoint[m], awh_bias, step, scale_weighthist, scale_pmfsum);

        /* Update the bias for this point only if there were skipped updates in the past to avoid calculating the log unneccessarily */
        if (bUpdated)
        {
            update_bias_point(&awh_bias->coordpoint[m]);
        }
    }
}
/* The global and neighborhood updates are the same except for the points touched. Could have one function
 * and a list of points to update instead. */
static void do_skipped_updates_in_neighborhood(awh_bias_t *awh_bias, gmx_int64_t step)
{
    double scale_weighthist, scale_pmfsum;
    int    nneighbors;
    int   *neighbor;

    get_histogram_scaling_skipped_update(awh_bias, &scale_weighthist, &scale_pmfsum);

    /* For each neighbor point of the center point, refresh its state by adding the results of all past, skipped updates. */
    nneighbors = awh_bias->grid->point[awh_bias->coord_value_index].nneighbors;
    neighbor   = awh_bias->grid->point[awh_bias->coord_value_index].neighbor;
    for (int n = 0; n < nneighbors; n++)
    {
        bool bUpdated = update_point_skipped(&awh_bias->coordpoint[neighbor[n]], awh_bias, step, scale_weighthist, scale_pmfsum);

        if (bUpdated)
        {
            update_bias_point(&awh_bias->coordpoint[neighbor[n]]);
        }
    }
}

static void merge_shared_update_lists(int *updateList, int *nupdate_ptr, int npoints, const gmx_multisim_t *ms)
{
    int *nupdate_point;
    int  nupdate_merged;

    /* Flag the update points of this sim */
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

/* Make an update list of all points that where touched since the last update */
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
           updating origin/end updatelist is changed (see sample_transition_probability_weights). */
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

static void reset_local_update_range(awh_bias_t *awh_bias)
{
    for (int d = 0; d < awh_bias->grid->ndim; d++)
    {
        awh_bias->origin_updatelist[d] = awh_bias->grid->point[awh_bias->coord_value_index].index[d];
        awh_bias->end_updatelist[d]    = awh_bias->grid->point[awh_bias->coord_value_index].index[d];
    }
}

/* Make a new point update */
static void update_point_new(coordpoint_t *coordpoint, const awh_bias_t *awh_bias,
                             gmx_int64_t step, double scale_weighthist, double scale_pmfsum)
{
    update_point(coordpoint, awh_bias, coordpoint->weightsum_iteration, scale_weighthist, scale_pmfsum);

    /* Remember that this was most recent update performed on this point. */
    coordpoint->last_update_index   = step/awh_bias->nstupdate_free_energy;
}


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

    /* Now add the counts and weights to the accumulating histograms.
       Note: we still need to use the weights for the update so we wait
       with resetting the histograms until the end of the update. */
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

static double get_new_histsize_initial(awh_bias_t *awh_bias, int awh_id, double t, bool bCovered, FILE *fplog)
{
    static const double exp_scale = 2;
    double              histsize_new, histsize_ref, nsamples_tot;
    bool                bExit;
    char                buf[STRLEN];

    if (!bCovered)
    {
        return awh_bias->histsize;
    }

    nsamples_tot = 0.;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        awh_bias->coordpoint[m].weightsum_covering = 0.;
        nsamples_tot                              += awh_bias->coordpoint[m].weightsum_tot;
    }

    histsize_new = awh_bias->histsize*exp_scale;

    /* Check if out of initial stage and reset N if so */
    histsize_ref  =  awh_bias->histsize_initial + nsamples_tot;
    bExit         = histsize_new >= histsize_ref;

    if (bExit)
    {
        histsize_new               = histsize_ref;
        awh_bias->in_initial       = FALSE;
    }

    if (fplog != NULL)
    {
        sprintf(buf, "\nawh%d:", awh_id + 1);
        fprintf(fplog, "%s covering at t = %g ps. Decreased the update size.\n",
                buf, t);

        if (bExit)
        {
            fprintf(fplog, "%s out of the initial stage. Update size will now grow continously.\n",
                    buf);
        }
        else
        {
            int ndoublings = ceil_log2(histsize_ref/histsize_new);
            fprintf(fplog, "%s at least %d more coverings needed to exit the initial stage\n",
                    buf, ndoublings);
        }
    }
    return histsize_new;
}

/* Return if the coordinate grid has been covered "enough" or not.
   Note: this could be improved, e.g. by using  path finding algorithm, e.g. Djikstra's algorithm. */
static bool is_covered(awh_bias_t *awh_bias)
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
        f_max = get_f_min(awh_bias) + awh_bias->target_param;
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

static double get_new_histsize(awh_bias_t *awh_bias, int awh_id, double t, gmx_int64_t step, bool *bCovered_ptr, FILE *fplog)
{
    double    histsize_new;
    bool      bCovered = FALSE;

    switch (awh_bias->eGrowth)
    {
        case eawhgrowthNONE:
            histsize_new = awh_bias->histsize;
            break;
        case eawhgrowthLINEAR:
            histsize_new = awh_bias->histsize + awh_bias->update_weight*awh_bias->weight_scaling;
            break;
        default:
            if (awh_bias->in_initial)
            {
                int nstcheck_covered;

                /* If covered enough, increase the histsize */
                nstcheck_covered = std::max((500*awh_bias->nstsample_coord/awh_bias->nstupdate_free_energy)*awh_bias->nstupdate_free_energy,
                                            awh_bias->nstupdate_free_energy);
                bCovered         = step % nstcheck_covered == 0 && is_covered(awh_bias);
                histsize_new     = get_new_histsize_initial(awh_bias, awh_id, t, bCovered, fplog);
            }
            else
            {
                histsize_new = awh_bias->histsize + awh_bias->update_weight*awh_bias->weight_scaling;
            }
    }

    *bCovered_ptr = bCovered;
    return histsize_new;
}

static void update(awh_bias_t *awh_bias, int awh_id, const gmx_multisim_t *ms, double t, gmx_int64_t step, FILE *fplog)
{
    int      nlocal_update, nupdate;
    bool     updateTarget, haveCovered;
    double   histsize_new;
    double   scale_weighthist_skipped, scale_weighthist_new, scale_pmfsum_skipped, scale_pmfsum_new;

    /* This array is only used in this scope and is always re-initialized */
    int     *updateList = awh_bias->updateList;

    /* Make a list of all points that could have been touched since the last update.
       For non-global updates, this is the final update list. For global updates this
       list is used for summing histograms below. */
    make_local_update_list(awh_bias, ms, updateList, &nlocal_update);

    /* Reset the range for the next update */
    reset_local_update_range(awh_bias);

    /* Add samples to histograms and sync sims if needed */
    sum_histograms(awh_bias, ms, updateList, nlocal_update);

    /* Update target distribution? */
    updateTarget = (awh_bias->eTarget != eawhtargetCONSTANT &&
                    do_at_step(awh_bias->nstupdate_target, step));

    /* The weighthistogram size after this update. */
    histsize_new = get_new_histsize(awh_bias, awh_id, t, step,
                                    &haveCovered, fplog);

    /* We update all points simultaneously either if 1) we covered (in the initial stage) since
       we don't want to keep track of when this happened and how it affects non-local points in future updates,
       or 2) we are updating the target distribution since its normalization has a global affect,
     */

    /*  Make the update list. For a global update, just add all points. For a local update,
        the update list equals the local update list of touched points. */
    if (updateTarget || haveCovered)
    {
        /* Global update. */
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
        /* Local update. The update list was set before summing the histograms. */
        nupdate = nlocal_update;
    }

    get_histogram_scaling_skipped_update(awh_bias, &scale_weighthist_skipped, &scale_pmfsum_skipped);
    get_histogram_scaling_new_update(awh_bias, histsize_new, awh_bias->histsize, haveCovered,
                                     &scale_weighthist_new, &scale_pmfsum_new);

    /* Update free energy and weight histogram for points in the list.
       Note: the bias is updated separately since it simply a function of the free energy and the target distribution. */
    for (int iupdate = 0; iupdate < nupdate; iupdate++)
    {
        coordpoint_t *coordpoint_to_update = &awh_bias->coordpoint[updateList[iupdate]];

        /* The previously skipped non-local updates */
        update_point_skipped(coordpoint_to_update, awh_bias, step, scale_weighthist_skipped, scale_pmfsum_skipped);

        /* The current local update. If we covered the update is different. */
        update_point_new(coordpoint_to_update, awh_bias, step, scale_weighthist_new, scale_pmfsum_new);

        /* Reset histograms collected for this update. */
        coordpoint_to_update->weightsum_iteration = 0;
        coordpoint_to_update->visits_iteration    = 0;
    }

    /* Only update the histogram size after we are done with the local point updates */
    awh_bias->histsize = histsize_new;

    /* The weight of new samples change if we rescale the weight of previous samples */
    awh_bias->log_relative_sampleweight -= std::log(scale_weighthist_new);

    if (updateTarget)
    {
        update_target(awh_bias);
    }

    /* Update the bias */
    for (int iupdate = 0; iupdate < nupdate; iupdate++)
    {
        update_bias_point(&awh_bias->coordpoint[updateList[iupdate]]);
    }
}

static void update_transition_probability_and_convolved_bias(awh_bias_t *awh_bias)
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
        weight[n]   = get_biased_weight_from_point(awh_bias, neighbor, awh_bias->coordpoint[neighbor].bias,
                                                   awh_bias->coord_value);
        weight_sum += weight[n];
    }
    inv_weight_sum = 1./weight_sum;

    /* Normalize probabilities to 1 */
    for (int n = 0; n < grid->point[awh_bias->coord_value_index].nneighbors; n++)
    {
        weight[n] *= inv_weight_sum;
    }

    /* The integral of the transition probabilities normalizes the distribution
       but is also the integrated bias at the current coordinate value which is needed
       later on for extracting the PMF. We avoid recalculating the exponentials below by
       returning this integral so that we can re-use it. */
    awh_bias->convolved_bias = std::log(weight_sum);
}

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
        weight_sum += get_biased_weight_from_point(awh_bias, neighbor, coordpoints[neighbor].bias,
                                                   coord_value);
    }

    /* Returns -GMX_DOUBLE_MAX if no neighboring points where in the target region. */
    return (weight_sum > 0) ? std::log(weight_sum) : -GMX_DOUBLE_MAX;
}

static void sample_transition_probability_weights(awh_bias_t *awh_bias)
{
    int m_neighbor_origin, m_neighbor_end;
    int nneighbors = awh_bias->grid->point[awh_bias->coord_value_index].nneighbors;

    /* Save weights for next update */
    for (int n = 0; n < nneighbors; n++)
    {
        int m_neighbor = awh_bias->grid->point[awh_bias->coord_value_index].neighbor[n];
        awh_bias->coordpoint[m_neighbor].weightsum_iteration += awh_bias->prob_weight_neighbor[n];
    }

    /* Keep track of which points will be affected by the next update. Only need to check the neighbors at the
       corners of the neighborhood grid.
       Note: if we always use one sample per update the update grid would be the same as the neighborhood grid.
       If the local updates are made efficient enpough there might not be any point to updating less often. */
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

static void do_sampling(awh_bias_t *awh_bias)
{
    /* Sampling-based deconvolution extracting the PMF */
    sample_pmf(awh_bias);

    /* Save probability weights for the update */
    sample_transition_probability_weights(awh_bias);
}


static void do_awh_bias_step(awh_bias_t *awh_bias, int awh_id,
                             bool convolveForce,
                             double *awhPotential, double *potential_jump,
                             const gmx_multisim_t *ms, double t, gmx_int64_t step, int seed, FILE *fplog)
{
    /* The grid point closest to the coordinate value defines the current neighborhood of points. Besides
       at steps when global updates and/or checks are performed, only the neighborhood will be touched. */
    awh_bias->coord_value_index = get_closest_index_in_grid(awh_bias->coord_value, awh_bias->grid);

    if (convolveForce || do_at_step(awh_bias->nstsample_coord, step))
    {
        /* In between updates we do things that require the current neighborhood of points to be up-to-date */
        do_skipped_updates_in_neighborhood(awh_bias, step);

        /* Update the transition probabilities for the neighborhood and save the convoluted coordinate bias which is
           used for extracting the PMF and equals convolved total potential energy (for all dimensions of the AWH coordinate). */
        update_transition_probability_and_convolved_bias(awh_bias);

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
        GMX_RELEASE_ASSERT(in_target_region(&awh_bias->coordpoint[awh_bias->coord_refvalue_index]),
                           "AWH bias coordinate reference value is outside of the target region.");
        potential =
            calcUmbrellaForceAndPotential(awh_bias, awh_bias->k, awh_bias->coord_refvalue_index, awh_bias->bias_force);

        /* Moving the umbrella results in a force correction and a new potential. */
        if (do_at_step(awh_bias->nstmove_refvalue, step))
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
        update(awh_bias, awh_id, ms, t, step, fplog);

        if (convolveForce)
        {
            /* The update results in a potential jump, so we need the new convolved potential. */
            newPotential = -calc_convolved_bias(awh_bias, awh_bias->coord_value)*awh_bias->invBeta;
        }
    }

    /* Add the the potential jump of this bias to the total sum. */
    *potential_jump += (newPotential - potential);

    /* Check the sampled data (histograms) and potentially warn user if something is suspicious */
    check_on_data(awh_bias, awh_id, t, fplog, step);
}

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
        do_awh_bias_step(&awh->awh_bias[k], k,
                         awh->convolveForce,
                         potential, &potentialJump,
                         ms, t, step, awh->seed, fplog);
    }

    /* Keep track of the total potential shift needed to remove the potential jumps. */
    awh->potential_offset -= potentialJump;
}

static void set_awh_coord_value(awh_bias_t *awh_bias, struct pull_t *pull_work, const t_pbc *pbc)
{
    /* Keep own copy of current coordinate value. */
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        get_pull_coord_value(pull_work, awh_bias->pull_coord_index[d], pbc, &awh_bias->coord_value[d]);
    }
}

static void set_awh_coord_value(awh_t *awh, struct pull_t *pull_work,
                                int ePBC, const matrix box)
{
    t_pbc    pbc;

    set_pbc(&pbc, ePBC, box);

    for (int k = 0; k < awh->nbias; k++)
    {
        set_awh_coord_value(&awh->awh_bias[k], pull_work, &pbc);
    }
}

real update_awh(awh_t                  *awh,
                const awh_params_t     *awh_params,
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

    /* Prepare AWH output data to later print to the energy file */
    if (time_to_write(step, awh->writer))
    {
        /* Make sure bias is up to date globally. This will also update the free energy and weight histogram. */
        for (int k = 0; k < awh->nbias; k++)
        {
            do_skipped_updates_for_all_points(awh->awh_bias, step);
        }

        prep_awh_output(awh->writer, awh_params, awh, ms);
    }

    /* Update the AWH coordinate values with those of the corresponding pull coordinates. */
    set_awh_coord_value(awh, pull_work, ePBC, box);

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

void write_awh_to_energyframe(t_enxframe *fr, const awh_t *awh)
{
    write_awh_to_frame(fr, awh->writer);
}
