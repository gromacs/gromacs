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

#include "correlation.h"
#include "data-writer.h"
#include "grid.h"
#include "history.h"
#include "internal.h"
#include "math.h"
#include "types.h"

//! String for registering AWH as an external potential with pull.
static const char *external_potential_string = "AWH";

//! AWH weight histogram type enum
enum {
    eawhweighthistREAL, eawhweighthistIDEAL, eawhweighthistNR
};

//! Linewidth used for warning output
static const int linewidth = 78;

//! Indent used for warning output
static const int indent    = 0;

static bool in_target_region(const coordpoint_t *coordpoint)
{
    return coordpoint->target > 0;
}

bool write_point_to_output(const awh_bias_t *awh_bias, int point_index)
{
    return in_target_region(&awh_bias->coordpoint[point_index]);
}

void get_pmf(const awh_bias_t *awh_bias, const gmx_multisim_t *ms, double *pmf)
{
    /* Need to temporarily exponentiate (~probability) to sum over simulations */
    for (int i = 0; i < awh_bias->npoints; i++)
    {
        if (in_target_region(&awh_bias->coordpoint[i]))
        {
            pmf[i] = std::exp(awh_bias->coordpoint[i].log_pmfsum);
        }
    }

    if (awh_bias->numSharedUpdate > 1)
    {
        gmx_sumd_sim(awh_bias->npoints, pmf, ms);
    }

    /* Take log again to get (non-normalized) PMF */
    for (int i = 0; i < awh_bias->npoints; i++)
    {
        if (in_target_region(&awh_bias->coordpoint[i]))
        {
            pmf[i] = -std::log(pmf[i]);
        }
        else
        {
            pmf[i] = 0;
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
                           bool                    bBlocklength_in_weight,
                           int                     ndomains)
{
    if (fplog != NULL)
    {
        char           awhstr[STRLEN];

        sprintf(awhstr, "\nawh%d:", awh_bias->biasIndex + 1);

        fprintf(fplog,
                "%s initial force correlation block length = %g %s"
                "%s force correlation number of blocks = %d",
                awhstr, get_blocklength(awh_bias->forcecorr),
                bBlocklength_in_weight ? "" : "ps",
                awhstr, get_nblocks(awh_bias->forcecorr));

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

static void set_free_energy_to_convolved_pmf(awh_bias_t *awh_bias)
{
    grid_point_t      *gridpoints  = awh_bias->grid->point;
    coordpoint_t      *coordpoints = awh_bias->coordpoint;

    for (int m = 0; m < awh_bias->npoints; m++)
    {
        double free_energy_weights = 0;
        for (int n = 0; n < gridpoints->nneighbors; n++)
        {
            int neighbor = gridpoints[m].neighbor[n];

            /* The log of the PMF histogram is the negative of the PMF,
               which is exactly the (positive) bias needed. */
            double neg_pmf_neighbor = coordpoints[neighbor].log_pmfsum;

            /* Add the convolved PMF weights for the neighbors of this point.
               Note that this function only adds point within the target > 0 reggon.
               Sum weights, take the logarithm last to get the free energy. */
            free_energy_weights += get_biased_weight_from_point(awh_bias, neighbor, neg_pmf_neighbor,
                                                                gridpoints[m].value);
        }

        GMX_RELEASE_ASSERT(free_energy_weights > 0, "Attempting to do log(< 0) in AWH convolved PMF calculation.");
        coordpoints[m].free_energy = -std::log(free_energy_weights);
    }
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

static void update_target(awh_bias_t *awh_bias, bool bInit)
{
    coordpoint_t       *coordpoint  = awh_bias->coordpoint;
    double              f_max       = 0, sum_target = 0, inv_sum;

    if (awh_bias->eTarget == eawhtargetCONSTANT && !bInit)
    {
        return;
    }

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
        else if (awh_bias->eTarget == eawhtargetWEIGHTHIST && !bInit)
        {
            coordpoint[m].target = bInit ? 1 : coordpoint[m].weightsum_ref;
        }

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
                             int nbias, int nsharing_sims)
{
    coordpoint_t *coordpoint = awh_bias->coordpoint;

    /* Apply default settings */
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        coordpoint[m].bias                                    = 0;
        coordpoint[m].free_energy                             = 0;
        coordpoint[m].target                                  = 1;
        coordpoint[m].target_constant_weight                  = 1;
        coordpoint[m].weightsum_iteration                     = 0;
        coordpoint[m].weightsum_covering                      = 0;
        coordpoint[m].weightsum_tot                           = 0;
        coordpoint[m].weightsum_ref                           = 0;
        coordpoint[m].last_update_index                       = 0;
        coordpoint[m].log_pmfsum                              = 0;
        coordpoint[m].visits_iteration                        = 0;
        coordpoint[m].visits_tot                              = 0;
    }

    /* Modify PMF, free energy and target distribution to user inputs values if there are any */
    if (awh_bias_params->bUser_data)
    {
        read_user_pmf_and_target(awh_bias, nbias);
        set_free_energy_to_convolved_pmf(awh_bias);
    }

    /* Initialize and normalize target. Note: depends on the free energy! */
    update_target(awh_bias, TRUE);

    /* Set the initial reference weighthistogram. Note: depends on the target distribution! */
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        coordpoint[m].weightsum_ref = awh_bias->histsize*coordpoint[m].target;
    }

    /* Note: the two initializations below depend on the target distribution and the free energy! */

    /* Set the bias */
    for (int m = 0; m < awh_bias->npoints; m++)
    {
        if (in_target_region(&coordpoint[m]))
        {
            update_bias_point(&coordpoint[m]);
        }
    }

    /* Make sure the pmf is normalized consistently with the histogram size. */
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
    awh_bias->eWeighthist               = awh_bias->eTarget == eawhtargetWEIGHTHIST ? eawhweighthistREAL : eawhweighthistIDEAL;

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

    /* The force constant k is needed for calculating the force correlation and the convoluted umbrella force. beta*k is
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
    init_coordpoints(awh_bias, awh_bias_params, nbias, awh_bias->numSharedUpdate);

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

    bool blocklength_in_weight   = false;

    /* We let the correlation init function set its parameters to something useful for now. */
    int    nblocks               = 0;
    double blocklength           = 0;

    awh_bias->forcecorr = init_correlation_grid(awh_bias->npoints, awh_bias->ndim,
                                                nblocks, blocklength,
                                                blocklength_in_weight,
                                                ir->delta_t,
                                                awh_bias->nstsample_coord);

    /* Print information about AWH variables that are set internally but might be of interest to the user. */
    print_log_init(awh_bias, awh_bias_params->dim_params, fplog,
                   blocklength_in_weight, ndomains);
}

static bool do_at_step(int nst, gmx_int64_t step)
{
    GMX_ASSERT(nst > 0, "All nst intervals in AWH should be > 0");

    return (step % nst == 0);
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


    /* Keep an array with the data to print to the energy file */
    awh->writer = init_awh_energywriter(awh_params->nstout, awh, ir->pull);

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
static double calc_umbrella_force(const awh_bias_t *awh_bias, const awh_dvec k, int point, awh_dvec force)
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
        calc_umbrella_force(awh_bias, awh_bias->k, m_neighbor, force_from_neighbor);

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

/* Calculates the umbrella force and increases the potential.
 * If this is a move step, also (attempt to) move the umbrella
 * and add the jump in the potential to potential_jump.
 */
static void set_umbrella_force(awh_bias_t *awh_bias, awh_dvec force,
                               double *potential, double *potential_jump,
                               gmx_int64_t step, int seed)
{
    double pot_cur =
        calc_umbrella_force(awh_bias, awh_bias->k, awh_bias->coord_refvalue_index, force);

    *potential += pot_cur;

    if (do_at_step(awh_bias->nstmove_refvalue, step))
    {
        awh_dvec new_force;

        /* Generate and set a new coordinate reference value */
        awh_bias->coord_refvalue_index =
            sample_coord_refvalue_index(awh_bias, step, seed);

        double pot_new =
            calc_umbrella_force(awh_bias, awh_bias->k, awh_bias->coord_refvalue_index, new_force);

        /*  A modification of the reference value at time t will lead to a different
            force over t-dt/2 to t and over t to t+dt/2. For high switching rates
            this means the force and velocity will change signs roughly as often.
            To avoid any issues we take the average of the previous and new force
            at steps when the reference value has been moved. E.g. if the ref. value
            is set every step to (coord dvalue +/- delta) would give zero force.
         */
        for (int d = 0; d < awh_bias->ndim; d++)
        {
            /* Average of the old and new force */
            force[d] = 0.5*(force[d] + new_force[d]);
        }

        *potential_jump += pot_new - pot_cur;
    }
}

static void update_force_correlation(correlation_grid_t *forcecorr, const awh_bias_t *awh_bias, double t)
{
    for (int n = 0; n < awh_bias->grid->point[awh_bias->coord_value_index].nneighbors; n++)
    {
        int      m_neighbor;
        double   weight_neighbor;
        awh_dvec force_from_neighbor;

        weight_neighbor = awh_bias->prob_weight_neighbor[n];
        m_neighbor      = awh_bias->grid->point[awh_bias->coord_value_index].neighbor[n];

        /* Add the force data of this neighbor point. Note: the sum of these forces is the convolved force.

           We actually add the force normalized by beta which has the units of 1/length. This means that the
           resulting correlation time integral is directly in units of friction time/length^2 which is really what
           we're interested in. */
        calc_umbrella_force(awh_bias, awh_bias->betak, m_neighbor, force_from_neighbor);

        /* Note: we might want to give a whole list of data to add instead and have this loop in the data adding function */
        add_data_to_correlation_matrix(forcecorr, m_neighbor, weight_neighbor, force_from_neighbor, t);
    }
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
    if (awh_bias->eWeighthist == eawhweighthistREAL)
    {
        coordpoint->weightsum_ref += weight_point*awh_bias->weight_scaling;
    }
    else
    {
        coordpoint->weightsum_ref += coordpoint->target*awh_bias->update_weight*awh_bias->weight_scaling;
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
    pointExists   = TRUE;
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
        awh_bias->origin_updatelist[d] = awh_bias->grid->point[awh_bias->coord_refvalue_index].index[d];
        awh_bias->end_updatelist[d]    = awh_bias->grid->point[awh_bias->coord_refvalue_index].index[d];
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

    /* Make a list of all points that could have been touched since the last update */
    make_local_update_list(awh_bias, ms, updateList, &nlocal_update);
    /* Reset the range for the next update */
    reset_local_update_range(awh_bias);

    /* Add samples to histograms and sync sims if needed */
    sum_histograms(awh_bias, ms, updateList, nlocal_update);

    /* Update target? */
    updateTarget = (awh_bias->eTarget != eawhtargetCONSTANT &&
                    do_at_step(awh_bias->nstupdate_target, step));

    /* The weighthistogram size after this update. */
    histsize_new = get_new_histsize(awh_bias, awh_id, t, step,
                                    &haveCovered, fplog);

    /* We update all points simultaneously either if 1) we covered (in the initial stage) since
       we don't want to keep track of when this happened and how it affects non-local points in future updates,
       or 2) we are updating the target distribution since its normalization has a global affect,
     */

    /*  Make the update list */
    if (updateTarget || haveCovered)
    {
        /* Perform a global update */
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
        update_target(awh_bias, FALSE);
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

static void do_sampling(awh_bias_t *awh_bias, double t)
{
    /* Force correlation */
    update_force_correlation(awh_bias->forcecorr, awh_bias, t);

    /* Sampling-based deconvolution extracting the PMF */
    sample_pmf(awh_bias);

    /* Save probability weights for the update */
    sample_transition_probability_weights(awh_bias);
}


static void do_awh_step(awh_bias_t *awh_bias, int awh_id,
                        bool convolveForce,
                        double *potential, double *potential_jump,
                        const gmx_multisim_t *ms, double t, gmx_int64_t step, int seed, FILE *fplog)
{
    GMX_RELEASE_ASSERT(in_target_region(&awh_bias->coordpoint[awh_bias->coord_refvalue_index]),
                       "The AWH bias coordinate reference value is outside of the target region.");

    if (convolveForce || do_at_step(awh_bias->nstsample_coord, step))
    {
        awh_bias->coord_value_index = get_closest_index_in_grid(awh_bias->coord_value, awh_bias->grid);

        /* In between updates we do things that require the current neighborhood of points to be up-to-date */
        do_skipped_updates_in_neighborhood(awh_bias, step);

        /* Update the transition probabilities for the neighborhood and save the convoluted coordinate bias which is
           used for extracting the PMF and equals convolved total potential energy (for all dimensions of the AWH coordinate). */
        update_transition_probability_and_convolved_bias(awh_bias);

        if (do_at_step(awh_bias->nstsample_coord, step))
        {
            do_sampling(awh_bias, t);
        }
    }
    /* The force on the coordinate resulting from the bias. */
    if (convolveForce)
    {
        calc_convolved_force(awh_bias, awh_bias->bias_force);
    }
    else
    {
        set_umbrella_force(awh_bias, awh_bias->bias_force,
                           potential, potential_jump,
                           step, seed);
    }

    /* Update the free energy estimates and bias and other history dependent method parameters */
    if (do_at_step(awh_bias->nstupdate_free_energy, step))
    {
        update(awh_bias, awh_id, ms, t, step, fplog);
    }

    /* Check the sampled data (histograms) and potentially warn user if something is suspicious */
    check_on_data(awh_bias, awh_id, t, fplog, step);
}

static void do_awh_step(awh_t                *awh,
                        double               *bias_potential,
                        const gmx_multisim_t *ms,
                        double                t,
                        gmx_int64_t           step,
                        FILE                 *fplog)
{
    double potential      = awh->potential_offset;
    double potential_jump = 0;

    for (int k = 0; k < awh->nbias; k++)
    {
        do_awh_step(&awh->awh_bias[k], k,
                    awh->convolveForce,
                    &potential, &potential_jump,
                    ms, t, step, awh->seed, fplog);
    }

    if (!awh->convolveForce)
    {
        /* Add the current potential jump to the shift to remove the jump */
        awh->potential_offset -= potential_jump;

        *bias_potential        = potential;
    }
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

static double get_bias_potential(const awh_t *awh)
{
    /* If there was an update the bias is modified which instantaneously changes the
       convolved bias potential at the current coordinate value. This shift is subtracted before
       adding the bias potential to the energy data. Note that a positive bias change corresponds
       to a negative potential change. */
    double bias_potential = awh->potential_offset;

    for (int k = 0; k < awh->nbias; k++)
    {
        awh_bias_t *awh_bias = &awh->awh_bias[k];

        /* Convert the bias potential from units of kT to energy units and add to
           the potential sum. */
        bias_potential += -awh_bias->convolved_bias*awh->invBeta;
    }

    return bias_potential;
}

static double convolved_bias_shift(const awh_t *awh, gmx_int64_t step)
{
    double bias_shift = 0;

    for (int k = 0; k < awh->nbias; k++)
    {
        awh_bias_t *awh_bias = &awh->awh_bias[k];

        /* If there was a free energy/bias update, update the bias potential shift. */
        if (do_at_step(awh_bias->nstupdate_free_energy, step))
        {
            /* The shift in the convolved bias for the same coordinate value. */
            bias_shift += (calc_convolved_bias(awh_bias, awh_bias->coord_value) - awh_bias->convolved_bias);
        }
    }

    return bias_shift;
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

    if (awh->convolveForce)
    {
        /* Get the bias potential, shifted to account for bias updates. */
        bias_potential = get_bias_potential(awh);

        /* Add the bias shift of this update if there was one. */
        awh->potential_offset += convolved_bias_shift(awh, step)*awh->invBeta;
    }

    wallcycle_stop(wallcycle, ewcAWH);

    return static_cast<real>(bias_potential);
}

void write_awh_to_energyframe(t_enxframe *fr, const awh_t *awh)
{
    write_awh_to_frame(fr, awh->writer);
}
