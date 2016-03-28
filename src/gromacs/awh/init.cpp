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

#include "awh.h"

#include <assert.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/enxio.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "grid.h"
#include "history.h"
#include "internal.h"
#include "math.h"
#include "types.h"

//! String for registering AWH as an external potential with pull.
static const char *external_potential_string = "AWH";

/* Get the multidimensional id of the domain that contains the current point */
static void get_domain_multidim_id_for_point(const awh_dim_params_t *awh_dim_params, const grid_t *grid, int ipoint, awh_ivec domain_id_dim)
{
    for (int d = 0; d < grid->ndim; d++)
    {
        domain_id_dim[d] = static_cast<int>(std::floor(1.*grid->point[ipoint].index[d]/grid->axis[d].npoints*awh_dim_params[d].ninterval));
    }
}

/* Get the min and max boundary points of the rectangular domain that the given point maps to.
   Each point maps to the domain that it is closest to the center of. */
static void get_domain_boundary_points(const awh_dim_params_t *awh_dim_params, const grid_t *grid, int point,
                                       int *point_min_ptr, int *point_max_ptr)
{
    awh_ivec          npoints_dim, domain_imin_dim, domain_imax_dim, domain_id_dim;

    get_domain_multidim_id_for_point(awh_dim_params, grid, point, domain_id_dim);

    /* Get the multidimensional min and max coordinates that defines the domain */
    for (int d = 0; d < grid->ndim; d++)
    {
        const double  non_overlap       = 1 - awh_dim_params[d].interval_overlap;
        const double  interval_fraction = 1./(non_overlap*(awh_dim_params[d].ninterval - 1) + 1);
        const double  interval_size     = grid->axis[d].npoints*interval_fraction;

        npoints_dim[d] = grid->axis[d].npoints;

        domain_imin_dim[d] = static_cast<int>(std::floor(interval_size*non_overlap*domain_id_dim[d]));
        domain_imax_dim[d] = static_cast<int>(std::floor(interval_size*(non_overlap*domain_id_dim[d] + 1)));
        domain_imax_dim[d] = std::min(domain_imax_dim[d], npoints_dim[d] - 1);

        /* Make sure that current point is inside of the domain */
        domain_imin_dim[d] = std::min(domain_imin_dim[d], grid->point[point].index[d]);
        domain_imax_dim[d] = std::max(domain_imax_dim[d], grid->point[point].index[d]);
    }

    /* Finally, convert the multidimensional indices to linear ones */
    (*point_min_ptr) = multidim_array_index_to_linear(domain_imin_dim, grid->ndim, npoints_dim);
    (*point_max_ptr) = multidim_array_index_to_linear(domain_imax_dim, grid->ndim, npoints_dim);
}

double coord_value_conversion_factor_internal2userinput(const awh_bias_t *awh_bias, const pull_params_t *pull_params, int awh_dimindex)
{
    return pull_coordinate_is_angletype(&pull_params->coord[awh_bias->pull_coord_index[awh_dimindex]]) ? RAD2DEG : 1;
}

double coord_value_conversion_factor_userinput2internal(const awh_bias_t *awh_bias, const pull_params_t *pull_params, int awh_dimindex)
{
    return pull_coordinate_is_angletype(&pull_params->coord[awh_bias->pull_coord_index[awh_dimindex]]) ? DEG2RAD : 1;
}

static void print_log_init(const awh_bias_t       *awh_bias,
                           const awh_dim_params_t *awh_dim_params,
                           FILE                   *fplog,
                           int                     ndomains)
{
    if (fplog != NULL)
    {
        char           awhstr[STRLEN];

        sprintf(awhstr, "\nawh%d:", awh_bias->biasIndex + 1);

        if (ndomains > 1)
        {
            int domain_imin, domain_imax;

            get_domain_boundary_points(awh_dim_params, awh_bias->grid, awh_bias->coord_value_index, &domain_imin, &domain_imax);

            fprintf(fplog, "%s partitioned the AWH domain into %d subdomains. This sim has target domain: ",
                    awhstr, ndomains);
            for (int d = 0; d < awh_bias->ndim; d++)
            {
                fprintf(fplog, "[%g, %g]",
                        awh_bias->grid->point[domain_imin].value[d], awh_bias->grid->point[domain_imax].value[d]);
                if (d < awh_bias->ndim - 1)
                {
                    fprintf(fplog, ", ");
                }
            }
        }
        fprintf(fplog, "\n");
    }
}

static void partition_domain(awh_bias_t *awh_bias, const awh_dim_params_t *awh_dim_params, int point_index)
{
    int                domain_imin, domain_imax;
    awh_ivec           domain_imin_dim, domain_imax_dim, npoints_dim;
    double             target_sum, inv_target_sum;
    coordpoint_t      *coordpoint         = awh_bias->coordpoint;

    /* Partition by zeroing the target distribution outside of the subinterval */
    get_domain_boundary_points(awh_dim_params, awh_bias->grid, point_index, &domain_imin, &domain_imax);

    /* Convert linear index to multidimensional index */
    for (int d = 0; d < awh_bias->grid->ndim; d++)
    {
        npoints_dim[d] = awh_bias->grid->axis[d].npoints;
    }
    linear_array_index_to_multidim(domain_imin, awh_bias->grid->ndim, npoints_dim, domain_imin_dim);
    linear_array_index_to_multidim(domain_imax, awh_bias->grid->ndim, npoints_dim, domain_imax_dim);

    target_sum = 0.;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        for (int d = 0; d < awh_bias->ndim; d++)
        {
            if (awh_bias->grid->point[m].index[d] < domain_imin_dim[d] || awh_bias->grid->point[m].index[d] > domain_imax_dim[d])
            {
                coordpoint[m].target = 0;
                coordpoint[m].bias   = -GMX_DOUBLE_MAX; /* the bias = log(target) + const = -infty */
            }
        }
        target_sum += coordpoint[m].target;
    }

    /* Renormalize target distribution to 1 and scale down the histogram size */
    awh_bias->histsize        *= target_sum;
    awh_bias->histsize_initial = awh_bias->histsize;
    inv_target_sum             = 1./target_sum;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        coordpoint[m].target *= inv_target_sum;
    }
}

static double get_initial_histsize_estimate(const awh_dvec betak, const awh_bias_params_t *awh_bias_params,
                                            const grid_t *grid, double dt_sample)
{
    int      ndim_skip;
    double   s, error_initial2, L2invD, histsize;
    awh_dvec x;

    /* Get diffusion factor */
    ndim_skip = 0;
    L2invD    = 0.;
    for (int d = 0; d < grid->ndim; d++)
    {
        double L = get_gridaxis_length(grid, d);
        if (L > 0)
        {
            L2invD += awh_bias_params->dim_params[d].diffusion/(L*L);
            s       = 1./std::sqrt(betak[d]);
            x[d]    = s/L;
        }
        else
        {
            /* Avoid division by zero for the rare case when there is only one point along one dimension */
            ndim_skip += 1;
            x[d]       = 1;
        }
    }
    L2invD         = 1./L2invD;
    error_initial2 = awh_bias_params->error_initial*awh_bias_params->error_initial;
    histsize       = L2invD*gaussian_geometry_factor(x, grid->ndim - ndim_skip)/(error_initial2*dt_sample);

    return histsize;
}

static void set_free_energy_to_convolved_pmf(awh_bias_t *awh_bias, const gmx_multisim_t *ms)
{
    double *convolvedPmf;

    snew(convolvedPmf, awh_bias->npoints);

    getConvolvedPmf(awh_bias, ms, convolvedPmf);

    for (int m = 0; m < awh_bias->npoints; m++)
    {
        awh_bias->coordpoint[m].free_energy = convolvedPmf[m];
    }

    sfree(convolvedPmf);
}

static void read_user_pmf_and_target(awh_bias_t *awh_bias, int nbias)
{
    char          filename[STRLEN];
    int           colindex_pmf, colindex_target;
    int           nrows, nrows_valid, ncols, ncols_min;
    int          *grid2data_index;
    double        val_prev;
    double      **data;

    /* Read the PMF and target distribution.
       From the PMF, the convolved PMF, or the reference value free energy, can be calculated
       base on the force constant. The free energy and target together determine the bias.
     */
    if (nbias == 1)
    {
        sprintf(filename, "awh-init.xvg");
    }
    else
    {
        sprintf(filename, "awh%d-init.xvg", awh_bias->biasIndex + 1);
    }

    nrows = read_xvg(filename, &data, &ncols);

    if (nrows <= 0)
    {
        gmx_fatal(FARGS, "%s is empty!", filename);
    }

    /* Make sure there are enough columns of data.

       Two formats are allowed. Either with columns  {coords, PMF, target} or
       {coords, PMF, x, y, z, target, ...}. The latter format is allowed since that
       is how AWH output is written (x, y, z being other AWH variables). For this format,
       trailing columns are ignored.
     */
    ncols_min    = awh_bias->ndim + 2;
    colindex_pmf = awh_bias->ndim;
    if (ncols == ncols_min)
    {
        colindex_target = colindex_pmf + 1;
    }
    else
    {
        colindex_target = colindex_pmf + 4;
    }

    /* Check the data */
    if (ncols < ncols_min)
    {
        gmx_fatal(FARGS, "The number of columns in %s (%d) should be at least %d. "
                  "The allowed formats are: [coords x ndim | PMF | target ] or "
                  "the format of an output AWH xvg-file.",
                  filename, ncols, ncols_min);
    }

    /* Make sure the coordinate values in the first dimension is monotonically increasing */
    val_prev    = data[0][0];
    nrows_valid = nrows;
    for (int m = 1; m < nrows; m++)
    {
        double val = data[0][m];
        if (val < val_prev)
        {
            /* read_xvg might return columns with zeros at the end of file if there were empty lines.
             * Check the remaining lines and ignore them if this is the case.
             */
            int      nskip       = 0;
            bool     bIgnore_row = TRUE;

            for (int n = m; n < nrows; n++)
            {
                for (int c = 0; c < ncols; c++)
                {
                    if (data[c][n] != 0)
                    {
                        bIgnore_row = FALSE;
                        break;
                    }
                }
                nskip++;
            }
            if (bIgnore_row)
            {
                char        warn_msg[STRLEN];
                nrows_valid -= nskip;
                sprintf(warn_msg, "Found %d empty lines in %s.\n", nskip, filename);
                gmx_warning(wrap_lines(warn_msg, linewidth, indent, FALSE));
                break;
            }
            else
            {
                gmx_fatal(FARGS, "The coordinate value of point %d, dimension 0 (%g) in %s smaller than"
                          " that of the previous point %d (%g)."
                          " The coordinate values should increase monotonically.",
                          m, val, filename, m-1, val_prev);
            }
        }
        else
        {
            val_prev = val;
        }
    }

    /* Less than 2 points is not useful, */
    if (nrows_valid < 2)
    {
        gmx_fatal(FARGS, "%s contains too few data points (%d)."
                  "The minimum number of points is 2.",
                  filename, nrows_valid);
    }

    /* Get a data point for each AWH grid point */
    snew(grid2data_index, awh_bias->grid->npoints);
    map_grid_to_datagrid(grid2data_index, data, nrows_valid, filename, awh_bias->grid);

    /* Extract the data for each grid point */
    bool target_is_zero = true;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        double target;

        awh_bias->coordpoint[m].log_pmfsum = -data[colindex_pmf][grid2data_index[m]];
        target = data[colindex_target][grid2data_index[m]];

        /* Check if the values are allowed. */
        if (target < 0)
        {
            gmx_fatal(FARGS, "Target distribution weight at point %d (%g) in %s is negative.",
                      m, target, filename);
        }
        if (target > 0)
        {
            target_is_zero = false;
        }
        awh_bias->coordpoint[m].target_constant_weight = target;
    }

    if (target_is_zero)
    {
        gmx_fatal(FARGS, "The target weights given in column %d in %s are all 0",
                  filename, colindex_target);
    }

    /* Free the arrays. */
    for (int m = 0; m < ncols; m++)
    {
        sfree(data[m]);
    }
    sfree(data);
    sfree(grid2data_index);
}

static void normalize_pmf(awh_bias_t *awh_bias,  int nsharing_sims)
{
    double expsum_pmf, expsum_f, nsamples, log_renorm;

    /* The normalization of the PMF estimate matters because it determines how big effect the next sample has.
       Approximately (for large enough force constant) we should have:
       sum_x(exp(-pmf(x)) = nsamples*sum_xref(exp(-f(xref)).
     */

    /* Calculate the normalization factor, i.e. divide by the pmf sum, multiply by the number of samples and the f sum */
    expsum_pmf = 0;
    expsum_f   = 0;
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        expsum_pmf = std::exp(awh_bias->coordpoint[m].log_pmfsum);
        expsum_f   = std::exp(-awh_bias->coordpoint[m].free_energy);
    }
    nsamples = awh_bias->histsize/nsharing_sims;

    /* Renormalize */
    log_renorm = std::log(nsamples*expsum_f/expsum_pmf);
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        awh_bias->coordpoint[m].log_pmfsum += log_renorm;
    }
}

static void init_coordpoints(awh_bias_t *awh_bias, const awh_bias_params_t *awh_bias_params,
                             int nbias, int nsharing_sims, const gmx_multisim_t *ms)
{
    coordpoint_t *coordpoint = awh_bias->coordpoint;

    /* Apply default settings. These may later be modified because of user input data,
       normalization, target type, etc. Note that later initializations may depend on
       these values being set here. */
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        coordpoint[m].bias                                    = 0;
        coordpoint[m].free_energy                             = 0;

        coordpoint[m].target                                  = 1;
        coordpoint[m].target_constant_weight                  = 1;

        coordpoint[m].weightsum_iteration                     = 0;
        coordpoint[m].weightsum_covering                      = 0;
        coordpoint[m].weightsum_tot                           = 0;
        coordpoint[m].weightsum_ref                           = 1;

        coordpoint[m].last_update_index                       = 0;

        coordpoint[m].log_pmfsum                              = 0;

        coordpoint[m].visits_iteration                        = 0;
        coordpoint[m].visits_tot                              = 0;
    }

    /* Modify PMF, free energy and the constant target distribution factor to user input
       values if there is data given. */
    if (awh_bias_params->bUser_data)
    {
        read_user_pmf_and_target(awh_bias, nbias);
        set_free_energy_to_convolved_pmf(awh_bias, ms);
    }

    /* We need to borrow update functions to initialize the target and bias.
       Note that free_energy needs to be set here. */
    initTargetAndBias(awh_bias);

    /* Set the initial reference weighthistogram. */
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        coordpoint[m].weightsum_ref = awh_bias->histsize*coordpoint[m].target;
    }

    /* Make sure the pmf is normalized consistently with the histogram size.
       Note: the target distribution and free energy need to be set here. */
    normalize_pmf(awh_bias, nsharing_sims);
}

static void init_awh_bias(FILE                       *fplog,
                          const t_inputrec           *ir,
                          const t_commrec            *cr,
                          int                         biasIndex,
                          int                         nbias,
                          awh_bias_t                 *awh_bias,
                          awh_bias_params_t          *awh_bias_params)
{
    int                  nneighbors_alloc, ndomains;
    awh_dvec             grid_origin, grid_end, grid_period;
    awh_dvec             coord_value_init;
    double               beta;
    const pull_params_t *pull_params = ir->pull;

    awh_bias->biasIndex       = biasIndex;

    /* Directly transfer some of the input parameters. Currently, we don't keep an a direct copy
       of the input data structure so that we can be flexible with the layout of the internal data structure. */
    awh_bias->eTarget         = awh_bias_params->eTarget;
    awh_bias->target_param    = awh_bias_params->target_param;
    awh_bias->eGrowth         = awh_bias_params->eGrowth;
    awh_bias->ndim            = awh_bias_params->ndim;
    awh_bias->numSharedUpdate = (awh_bias_params->bShare && (cr != NULL) && MULTISIM(cr) ? cr->ms->nsim : 1);
    awh_bias->update_weight   = awh_bias->nstupdate_free_energy/awh_bias->nstsample_coord*awh_bias->numSharedUpdate;
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        double user2internal_units = coord_value_conversion_factor_userinput2internal(awh_bias, pull_params, d);
        awh_bias->pull_coord_index[d]    = awh_bias_params->dim_params[d].pull_coord_index;
        grid_origin[d]                   = user2internal_units*awh_bias_params->dim_params[d].origin;
        grid_end[d]                      = user2internal_units*awh_bias_params->dim_params[d].end;
        grid_period[d]                   = user2internal_units*awh_bias_params->dim_params[d].period;
        coord_value_init[d]              = user2internal_units*awh_bias_params->dim_params[d].coord_value_init;
    }
    awh_bias->log_relative_sampleweight = 0;
    awh_bias->in_initial                = awh_bias->eGrowth == eawhgrowthEXP_LINEAR;

    /* When the target distribtution equals the reference histogram there are special settings needed for it to be
       useful:

       1) The target parameter equals the effective temperature factor T/T_eff which is fed into the biasing by scaling
       down the sampled weight.

       2) The reference weight histogram is updated by the adding the actually sampled weights (scaled down by weight_scaling).
     */
    awh_bias->weight_scaling            = awh_bias->eTarget == eawhtargetWEIGHTHIST ? awh_bias->target_param : 1;
    awh_bias->idealWeighthistUpdate     = awh_bias->eTarget != eawhtargetWEIGHTHIST;

    /* Set the target update frequency based on the target distrbution type (this could be made a user-option but
       there is most likely no big need for tweaking this for most users. */
    switch (awh_bias->eTarget)
    {
        case eawhtargetCONSTANT:
            awh_bias->nstupdate_target = 0;
            break;
        case eawhtargetCUTOFF:
        case eawhtargetBOLTZMANN:
            /* Updating the target generally requires updating the whole grid so to keep the cost down
               we generall update the target less often than the free energy (unless the free energy
               update step is set to > 100 samples). */
            awh_bias->nstupdate_target = std::max((100*awh_bias->nstsample_coord/awh_bias->nstupdate_free_energy)*awh_bias->nstupdate_free_energy,
                                                  awh_bias->nstupdate_free_energy);
            break;
        case eawhtargetWEIGHTHIST:
            /* If the target distribution is set equal to the reference histogram which is updated every free energy update,
               then the target has to be updated at the same times. */
            awh_bias->nstupdate_target = awh_bias->nstupdate_free_energy;
            break;
        default:
            gmx_incons("Unknown AWH target type");
    }

    /* The force constant k is needed for calculating the force. beta*k is
     * basically the inverse variance of the coordinate and is e.g. used for defining the grid spacing */
    beta = 1./(BOLTZ*ir->opts.ref_t[0]);
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        awh_bias->k[d]     = pull_params->coord[awh_bias->pull_coord_index[d]].k;
        awh_bias->betak[d] = beta*awh_bias->k[d];
    }

    /* Initialize grid and the functions of the grid */
    awh_bias->grid        = init_grid(awh_bias->ndim, awh_bias->betak, grid_origin, grid_end, grid_period, &nneighbors_alloc);
    awh_bias->npoints     = awh_bias->grid->npoints;

    /* The transition probability weights in a neighborhood of the current coordinate */
    snew(awh_bias->prob_weight_neighbor, nneighbors_alloc);

    /* Estimate and initialize histsize_initial. The estimation depends on the grid. */
    awh_bias->histsize_initial = get_initial_histsize_estimate(awh_bias->betak, awh_bias_params, awh_bias->grid, awh_bias->nstsample_coord*ir->delta_t);
    awh_bias->histsize         = awh_bias->histsize_initial;

    /* Set initial coordinate reference value to the one closest to the initial reference value given in pull.
       More correctly one would sample from the biased distribution, but it doesn't really matter. */
    awh_bias->coord_value_index    = get_closest_index_in_grid(coord_value_init, awh_bias->grid);
    awh_bias->coord_refvalue_index = awh_bias->coord_value_index;

    /* The minimum and maximum multidimensional point indices that are affected by the next update */
    for (int d = 0; d < awh_bias->ndim; d++)
    {
        awh_bias->origin_updatelist[d] = awh_bias->grid->point[awh_bias->coord_value_index].index[d];
        awh_bias->end_updatelist[d]    = awh_bias->grid->point[awh_bias->coord_value_index].index[d];
    }

    /* Initialize free energy functions and biases */
    snew(awh_bias->coordpoint, awh_bias->npoints);
    snew(awh_bias->updateList, awh_bias->npoints);
    init_coordpoints(awh_bias, awh_bias_params, nbias, awh_bias->numSharedUpdate, (cr != NULL ? cr->ms : NULL));

    /* Partition the AWH domain if any of the dimensions is set to be divided into more than 1 interval. */
    ndomains = 1;
    for (int d = 0; d < awh_bias_params->ndim; d++)
    {
        ndomains *= awh_bias_params->dim_params[d].ninterval;
    }
    if (ndomains > 1)
    {
        /* Partitioning  modifies the target distribution and the weighthistogram outside of the target domain.
           The target domain is determined based on the current coordinate reference value. */
        partition_domain(awh_bias, awh_bias_params->dim_params, awh_bias->coord_value_index);
    }

    /* Print information about AWH variables that are set internally but might be of interest to the user. */
    print_log_init(awh_bias, awh_bias_params->dim_params, fplog, ndomains);
}



void register_bias_with_pull(const awh_t *awh, struct pull_t *pull_work)
{
    for (int k = 0; k < awh->nbias; k++)
    {
        awh_bias_t *awh_bias = &awh->awh_bias[k];

        for (int d = 0; d < awh_bias->ndim; d++)
        {
            register_external_pull_potential(pull_work, awh_bias->pull_coord_index[d], external_potential_string);
        }
    }
}

awh_t *init_awh(FILE                    *fplog,
                const t_inputrec        *ir,
                const t_commrec         *cr,
                const awh_params_t      *awh_params)
{
    awh_t         *awh;
    const int      nstsample_coord        = awh_params->nstsample_coord;
    const int      nstmove_refvalue       = awh_params->nsamples_move_refvalue*nstsample_coord;
    const int      nstupdate_free_energy  = awh_params->nsamples_update_free_energy*nstsample_coord;

    snew(awh, 1);

    /* Copy over the parameters */
    awh->nbias            = awh_params->nbias;
    awh->convolveForce    = (awh_params->ePotential == eawhpotentialCONVOLVED);
    awh->invBeta          = BOLTZ*ir->opts.ref_t[0];
    awh->seed             = awh_params->seed;
    awh->potential_offset = 0;

    snew(awh->awh_bias, awh->nbias);
    for (int k = 0; k < awh->nbias; k++)
    {
        /* Some things are common for all biases */
        awh->awh_bias[k].nstsample_coord          = nstsample_coord;
        awh->awh_bias[k].nstmove_refvalue         = nstmove_refvalue;
        awh->awh_bias[k].nstupdate_free_energy    = nstupdate_free_energy;

        init_awh_bias(fplog, ir, cr, k, awh->nbias, &awh->awh_bias[k], &awh_params->awh_bias_params[k]);
    }

    return awh;
}

awh_t *init_awh_md(FILE                    *fplog,
                   const t_inputrec        *ir,
                   const t_commrec         *cr,
                   const awh_params_t      *awh_params,
                   awh_history_t           *awh_history,
                   struct pull_t           *pull_work,
                   bool                     startingFromCheckpoint)
{
    awh_t *awh = init_awh(fplog, ir, cr, awh_params);

    /* Need to register the AWH coordinates to be allowed to apply forces to the pull coordinates. */
    register_bias_with_pull(awh, pull_work);

    if (startingFromCheckpoint)
    {
        init_awh_history_from_checkpoint(awh_history, cr);
        restore_awh_state_from_history(awh_history, awh);
    }
    else
    {
        init_awh_history_from_state(awh_history, awh);
    }

    return awh;
}
