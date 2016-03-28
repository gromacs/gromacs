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
 * Implements the initialization of the Bias class and its helpers.
 *
 * \author Viveca Lindahl
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_awh
 */

#include "gmxpre.h"

#include <assert.h>

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <algorithm>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/math/units.h"
#include "gromacs/mdtypes/awh-history.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/smalloc.h"

#include "bias.h"
#include "biaswriter.h"
#include "correlation.h"
#include "grid.h"
#include "internal.h"
#include "math.h"
#include "pointstate.h"

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
        domainIdDim[d] = static_cast<int>(std::floor(1.*grid.point(point).index[d]/grid.axis(d).npoints*awh_dim_params[d].ninterval));
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
 * \param[in] domain_id      Domain id.
 * \param[in,out] point_min  Mininmum boundary point index.
 * \param[in,out] point_max  Maximum boundary point index.
 */
static void getDomainBoundary(const awh_bias_params_t &awhBiasParams, const Grid &grid,
                              int domain_id,
                              int *point_min, int *point_max)
{
    const awh_dim_params_t *awh_dim_params = awhBiasParams.dim_params;

    /* Convert the given linear domain id to a multidimensional one. */
    awh_ivec ninterval_dim, domain_id_dim = {0};

    for (int d = 0; d < awhBiasParams.ndim; d++)
    {
        ninterval_dim[d] = awh_dim_params[d].ninterval;
    }
    linear_array_index_to_multidim(domain_id, grid.ndim(), ninterval_dim, domain_id_dim);

    /* Get the multidimensional min and max coordinates that defines the domain */
    awh_ivec          npoints_dim, domain_imin_dim, domain_imax_dim;

    for (int d = 0; d < grid.ndim(); d++)
    {
        const double  non_overlap       = 1 - awh_dim_params[d].interval_overlap;
        const double  interval_fraction = 1./(non_overlap*(awh_dim_params[d].ninterval - 1) + 1);
        const double  interval_size     = grid.axis(d).npoints*interval_fraction;

        npoints_dim[d] = grid.axis(d).npoints;

        /* The domain end points are given by the domain id scaled by the interval size and a non-overlap factor.
           For overlap > 0, i.e. non-overlap < 1, the interval is extended in both directions. */
        domain_imin_dim[d] = static_cast<int>(std::floor(interval_size*non_overlap*domain_id_dim[d]));

        /* Note that the max is just the min + the interval size. */
        domain_imax_dim[d] = static_cast<int>(std::floor(interval_size*(non_overlap*domain_id_dim[d] + 1)));
        domain_imax_dim[d] = std::min(domain_imax_dim[d], npoints_dim[d] - 1);
    }

    /* Finally, convert the multidimensional point indices to linear ones */
    (*point_min) = multidim_array_index_to_linear(domain_imin_dim, grid.ndim(), npoints_dim);
    (*point_max) = multidim_array_index_to_linear(domain_imax_dim, grid.ndim(), npoints_dim);
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

    const Grid &grid  = *bias.grid();

    int         my_id = getDomainId(awhBiasParams, grid, bias.state().gridpointIndex);

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
                    bias.dimParams(d).scaleInternalToUserInput(grid.point(domain_imin).coordValue[d]),
                    bias.dimParams(d).scaleInternalToUserInput(grid.point(domain_imax).coordValue[d]));
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
                int coord_imin = static_cast<int>(std::floor(1.0*id_dim[d]/ninterval_dim[d]*(bias.pointState().size() - 1)));
                int coord_imax = static_cast<int>(std::floor(1.0*(id_dim[d] + 1)/ninterval_dim[d]*(bias.pointState().size() - 1)));

                fprintf(stderr, "[%.4g, %.4g]",
                        bias.dimParams(d).scaleInternalToUserInput(grid.point(coord_imin).coordValue[d]),
                        bias.dimParams(d).scaleInternalToUserInput(grid.point(coord_imax).coordValue[d]));

                if (d < bias.ndim() - 1)
                {
                    fprintf(stderr, " x ");
                }
            }
            fprintf(stderr, " --> domain id %d: ", id);
            for (int d = 0; d < bias.ndim(); d++)
            {
                fprintf(stderr, "[%.4g, %.4g]",
                        bias.dimParams(d).scaleInternalToUserInput(grid.point(domain_imin).coordValue[d]),
                        bias.dimParams(d).scaleInternalToUserInput(grid.point(domain_imax).coordValue[d]));
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
 * \param[in] bias                        The AWH bias.
 * \param[in] awhBiasParams               AWH bias parameters.
 * \param[in] measureBlocklengthInWeight  True if correlation block length should be measured in weight, otherwise time.
 * \param[in,out] fplog                   Log file.
 */
static void printInitializationToLog(const Bias              &bias,
                                     const awh_bias_params_t &awhBiasParams,
                                     bool                     measureBlocklengthInWeight,
                                     FILE                    *fplog)
{
    char           awhstr[STRLEN];
    sprintf(awhstr, "\nawh%d:", bias.biasIndex + 1);

    if (fplog != nullptr)
    {
        fprintf(fplog,
                "%s initial force correlation block length = %g %s"
                "%s force correlation number of blocks = %d",
                awhstr, getBlockLength(bias.forceCorr()),
                measureBlocklengthInWeight ? "" : "ps",
                awhstr, getNumBlocks(bias.forceCorr()));
    }


    printPartitioningDomainInit(awhstr, bias, awhBiasParams, fplog);
}

/*! \brief
 * Partition sampling domain.
 *
 * \param[in] awhBiasParams  AWH bias parameters.
 * \param[in] grid           The grid setup.
 * \param[in] point_index    Point determining the domain id after partitioning.
 * \param[in,out] pointState The state of the points in the bias.
 */
static double partitionDomain(const awh_bias_params_t &awhBiasParams, const Grid &grid, int point_index,
                              std::vector<PointState> *pointState)
{
    int         domain_imin, domain_imax;
    awh_ivec    domain_imin_dim, domain_imax_dim, npoints_dim;

    /* Partition by zeroing the target distribution outside of the subinterval */
    getDomainBoundaryForPoint(awhBiasParams, grid, point_index, &domain_imin, &domain_imax);

    /* Convert linear index to multidimensional index */
    for (int d = 0; d < grid.ndim(); d++)
    {
        npoints_dim[d] = grid.axis(d).npoints;
    }
    linear_array_index_to_multidim(domain_imin, grid.ndim(), npoints_dim, domain_imin_dim);
    linear_array_index_to_multidim(domain_imax, grid.ndim(), npoints_dim, domain_imax_dim);

    double targetSum = 0.;
    for (size_t m = 0; m < pointState->size(); m++)
    {
        PointState &ps = (*pointState)[m];
        for (int d = 0; d < grid.ndim(); d++)
        {
            int index_d = grid.point(m).index[d];
            if (index_d < domain_imin_dim[d] ||
                index_d > domain_imax_dim[d])
            {
                ps.setTargetToZero();
            }
        }
        targetSum += ps.target();
    }

    /* Renormalize target distribution to 1.
     * NOTE: The total histogram size in state still needs to be scaled down.
     */
    double invTargetSum = 1/targetSum;
    for (auto &ps : *pointState)
    {
        ps.scaleTarget(invTargetSum);
    }

    return targetSum;
}

/*! \brief
 * Estimate a reasonable initial reference weight histogram size.
 *
 * \param[in] dimParams      Parameters for the dimensions of the coordinate.
 * \param[in] awhBiasParams  Bias parameters.
 * \param[in] gridAxis       The Grid axes.
 * \param[in] dt_sample      Sampling frequency of probability weights.
 * \returns estimate of initial histogram size.
 */
static double getInitialHistSizeEstimate(const std::vector<DimParams> &dimParams, const awh_bias_params_t &awhBiasParams,
                                         const std::vector<GridAxis> &gridAxis, double dt_sample)
{
    int      ndim_skip;
    double   s, error_initial2, L2invD, histSize;
    awh_dvec x;

    /* Get diffusion factor */
    ndim_skip = 0;
    L2invD    = 0.;
    for (size_t d = 0; d < gridAxis.size(); d++)
    {
        double L = getGridaxisLength(gridAxis[d]);
        if (L > 0)
        {
            L2invD += awhBiasParams.dim_params[d].diffusion/(L*L);
            s       = 1./std::sqrt(dimParams[d].betak);
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
    error_initial2 = awhBiasParams.error_initial*awhBiasParams.error_initial;
    histSize       = L2invD*gaussian_geometry_factor(x, gridAxis.size() - ndim_skip)/(error_initial2*dt_sample);

    return histSize;
}

/*! \brief
 * Convolves the PMF and sets the initial free energy to its convolution.
 *
 * \param[in] bias            The AWH bias.
 * \param[in] ms              Struct for multi-simulation communication.
 * \param[in,out] pointState  The state of all points.
 */
static void setFreeEnergyToConvolvedPmf(const Bias &bias, const gmx_multisim_t *ms,
                                        std::vector<PointState> *pointState)
{
    std::vector<float> convolvedPmf;

    getConvolvedPmf(bias, ms, &convolvedPmf);

    for (size_t m = 0; m < pointState->size(); m++)
    {
        (*pointState)[m].setFreeEnergy(convolvedPmf[m]);
    }
}

/*! \brief
 * Find trailing data rows containing only zeros.
 *
 * \param[in] data    2D data array.
 * \param[in] nrows   Number of rows in array.
 * \param[in] ncols   Number of cols in array.
 * \returns the number of trailing zero rows.
 */
static int findTrailingZeroRows(const double* const *data, int nrows, int ncols)
{
    int nZeroRows = 0;
    for (int m = nrows - 1; m >= 0; m--)
    {
        bool rowIsZero = true;
        for (int d = 0; d < ncols; d++)
        {
            if (data[d][m] != 0)
            {
                rowIsZero = false;
                break;
            }
        }

        if (!rowIsZero)
        {
            /* At a row with non-zero data */
            break;
        }
        else
        {
            /* Still at a zero data row, keep checking rows higher up. */
            nZeroRows++;
        }
    }

    return nZeroRows;
}

/*! \brief
 * Initializes the PMF and target with data read from an input table.
 *
 * \param[in] bias            The AWH Bias.
 * \param[in] nbias           Number of biases.
 * \param[in,out] pointState  The state of the points in this bias.
 */
static void readUserPmfAndTargetDistribution(const Bias &bias, int nbias,
                                             std::vector<PointState> *pointState)
{
    char          filename[STRLEN];

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
        sprintf(filename, "awh%d-init.xvg", bias.biasIndex + 1);
    }

    char buf[STRLEN];
    sprintf(buf,
            "%s is expected in the following format. "
            "The first ndim column(s) should contain the coordinate values for each point, "
            " each column containing values of one dimension (in ascending order). "
            "For a multidimensional coordinate, points should be listed "
            "in the order obtained by traversing lower dimensions first. "
            "E.g. for two-dimensional grid of size nxn: "
            "(1, 1), (1, 2),..., (1, n), (2, 1), (2, 2), ..., , (n, n - 1), (n, n). "
            "Column ndim +  1 should contain the PMF value for each coordinate value. "
            "The target distribution values should be in column ndim + 2  or column ndim + 5. "
            "Make sure there input file ends with a new line but has no trailing new lines.",
            filename);
    char correctFormatMessage[STRLEN];
    sprintf(correctFormatMessage, "%s", wrap_lines(buf, linewidth, indent, FALSE));

    double  **data;
    int       nrows, ncols;
    nrows = read_xvg(filename, &data, &ncols);

    /* Check basic data properties here. Grid takes care of more complicated things. */

    if (nrows <= 0)
    {
        gmx_fatal(FARGS, "%s is empty!.\n\n%s", filename, correctFormatMessage);
    }

    /* Less than 2 points is not useful for PMF or target. */
    if (nrows <  2)
    {
        gmx_fatal(FARGS, "%s contains too few data points (%d)."
                  "The minimum number of points is 2.",
                  filename, nrows);
    }

    /* Make sure there are enough columns of data.

       Two formats are allowed. Either with columns  {coords, PMF, target} or
       {coords, PMF, x, y, z, target, ...}. The latter format is allowed since that
       is how AWH output is written (x, y, z being other AWH variables). For this format,
       trailing columns are ignored.
     */
    int colindex_target;
    int ncols_min    = bias.ndim() + 2;
    int colindex_pmf = bias.ndim();
    if (ncols == ncols_min)
    {
        colindex_target = colindex_pmf + 1;
    }
    else
    {
        colindex_target = colindex_pmf + 4;
    }

    if (ncols < ncols_min)
    {
        gmx_fatal(FARGS, "The number of columns in %s (%d) should be at least %d."
                  "\n\n%s",
                  filename, correctFormatMessage);
    }

    /* read_xvg can give trailing zero data rows for trailing new lines in the input. We allow 1 zero row,
       since this could be real data. But multiple trailing zero rows cannot correspond to valid data. */
    int nZeroRows = findTrailingZeroRows(data, nrows, ncols);
    if (nZeroRows > 1)
    {
        gmx_fatal(FARGS, "Found %d trailing zero data rows in %s. Please remove trailing empty lines and try again.",
                  nZeroRows, filename);
    }

    /* Convert from user units to internal units before sending the data of to grid. */
    for (int d = 0; d < bias.ndim(); d++)
    {
        double scalingFactor = bias.dimParams(d).scaleUserInputToInternal(1);
        if (scalingFactor == 1)
        {
            continue;
        }
        for (size_t m = 0; m < pointState->size(); m++)
        {
            data[d][m] *= scalingFactor;
        }
    }

    /* Get a data point for each AWH grid point so that they all get data. */
    std::vector<int> grid2data_index;
    grid2data_index.resize(bias.grid()->numPoints());
    mapGridToDatagrid(&grid2data_index, data, nrows, filename, bias.grid(), correctFormatMessage);

    /* Extract the data for each grid point */
    bool target_is_zero = true;
    for (size_t m = 0; m < pointState->size(); m++)
    {
        double target;

        (*pointState)[m].setLogPmfsum(-data[colindex_pmf][grid2data_index[m]]);
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
        (*pointState)[m].setTargetConstantWeight(target);
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
}

/*! \brief
 * Normalize the PMF histogram.
 *
 * \param[in,out] pointState  The bias points to normalize.
 * \param[in] state           The global state of the bias.
 * \param[in] numSharingSims  The number of simulations sharing the bias.
 */
static void normalizePmf(std::vector<PointState> *pointState,
                         const BiasState &state, int numSharingSims)
{
    /* The normalization of the PMF estimate matters because it determines how big effect the next sample has.
       Approximately (for large enough force constant) we should have:
       sum_x(exp(-pmf(x)) = nsamples*sum_xref(exp(-f(xref)).
     */

    /* Calculate the normalization factor, i.e. divide by the pmf sum, multiply by the number of samples and the f sum */
    double expsum_pmf = 0;
    double expsum_f   = 0;
    for (auto &ps : *pointState)
    {
        if (ps.inTargetRegion())
        {
            expsum_pmf += std::exp( ps.logPmfsum());
            expsum_f   += std::exp(-ps.freeEnergy());
        }
    }
    double numSamples = state.histSize/numSharingSims;

    /* Renormalize */
    double logRenorm = std::log(numSamples*expsum_f/expsum_pmf);
    for (auto &ps : *pointState)
    {
        if (ps.inTargetRegion())
        {
            ps.setLogPmfsum(ps.logPmfsum() + logRenorm);
        }
    }
}

/*! \brief
 * Initialize the state of grid coordinate points.
 *
 * \param[in] bias            The AWH bias.
 * \param[in] awhBiasParams   Bias parameters.
 * \param[in] nbias           Number of biases.
 * \param[in] numSharingSims  The number of simulations sharing the bias.
 * \param[in] ms              Struct for multi-simulation communication.
 * \param[in,out] pointState  The state of the points in the bias.
 */
static void initGridPointState(const Bias &bias, const awh_bias_params_t &awhBiasParams,
                               int nbias, int numSharingSims, const gmx_multisim_t *ms,
                               std::vector<PointState> *pointState)
{
    /* Modify PMF, free energy and the constant target distribution factor
     * to user input values if there is data given.
     */
    if (awhBiasParams.bUser_data)
    {
        readUserPmfAndTargetDistribution(bias, nbias, pointState);
        setFreeEnergyToConvolvedPmf(bias, ms, pointState);
    }

/* The local Boltzmann distribution is special because the target distribution is updated as a function of the reference weighthistogram. */
    GMX_ASSERT((bias.params().eTarget != eawhtargetLOCALBOLTZMANN) ||
               (bias.params().eTarget == eawhtargetLOCALBOLTZMANN && (*pointState)[0].weightsumRef() != 0),
               "AWH reference weight histogram not initialized properly with local Boltzmann target distribution." );

    updateTarget(pointState, bias.params());

    for (auto &ps : *pointState)
    {
        if (ps.inTargetRegion())
        {
            ps.updateBias();
        }
        else
        {
            /* Note that for zero target this is a value that represents -infinity but should not be used for biasing. */
            ps.setTargetToZero();
        }
    }

    /* Set the initial reference weighthistogram. */
    for (auto &ps : *pointState)
    {
        ps.setInitialReferenceWeightHistogram(bias.state());
    }

    /* Make sure the pmf is normalized consistently with the histogram size.
       Note: the target distribution and free energy need to be set here. */
    normalizePmf(pointState, bias.state(), numSharingSims);
}

BiasParams::BiasParams(const awh_params_t           &awhParams,
                       const awh_bias_params_t      &awhBiasParams,
                       const std::vector<DimParams> &dimParams,
                       double                        beta,
                       double                        mdTimeStep,
                       DisableUpdateSkips            disableUpdateSkips,
                       const t_commrec              *cr,
                       const std::vector<GridAxis>  &gridAxis) :
    numStepsSampleCoord(awhParams.nstsample_coord),
    numSamplesUpdateFreeEnergy(awhParams.nsamples_update_free_energy),
    eTarget(awhBiasParams.eTarget),
    convolveForce(awhParams.ePotential == eawhpotentialCONVOLVED),
    disableUpdateSkips_(disableUpdateSkips)
{
    switch (eTarget)
    {
        case eawhtargetCUTOFF:
            targetParam = awhBiasParams.targetCutoff;
            break;
        case eawhtargetBOLTZMANN:
        case eawhtargetLOCALBOLTZMANN:
            targetParam = awhBiasParams.targetBetaScaling;
            break;
        default:
            targetParam = 0;
            break;
    }

    int numMultiSims = ((cr != NULL) && MULTISIM(cr)) ? cr->ms->nsim : 1;
    numSharedUpdate  = (awhBiasParams.bShare ? numMultiSims : 1);
    update_weight    = numSamplesUpdateFreeEnergy*numSharedUpdate;

    invBeta          = 1./beta;

    /* The local Boltzmann target distibution is defined by
       1) Adding the sampled weights instead of the target weights to the reference weight histogram.
       2) Scaling the weights of these samples by the beta scaling factor.
       3) Setting the target distribution equal the reference weight histogram.
       This requires the following special update settings.
     */
    localWeightScaling      = (eTarget == eawhtargetLOCALBOLTZMANN ? targetParam : 1);
    idealWeighthistUpdate   = (eTarget != eawhtargetLOCALBOLTZMANN);
    /* Note: these variables could in principle be set to something else also for other target distribution types.
       However, localWeightScaling < 1  is in general expected to give lower efficiency and, except for local Boltzmann,
       idealWeightHistUpdate = false gives (in my experience) unstable, non-converging results. */

    /* Set the target update frequency based on the target distrbution type (this could be made a user-option but
       there is most likely no big need for tweaking this for most users. */
    switch (eTarget)
    {
        case eawhtargetCONSTANT:
            nstupdate_target = 0;
            break;
        case eawhtargetCUTOFF:
        case eawhtargetBOLTZMANN:
            /* Updating the target generally requires updating the whole grid so to keep the cost down
               we generally update the target less often than the free energy (unless the free energy
               update step is set to > 100 samples). */
            nstupdate_target = std::max(100 % numSamplesUpdateFreeEnergy,
                                        numSamplesUpdateFreeEnergy)*numStepsSampleCoord;
            break;
        case eawhtargetLOCALBOLTZMANN:
            /* The target distribution is set equal to the reference histogram which is updated every free energy update.
               So the target has to be updated at the same time. This leads to a global update each time because it is
               assumed that a target distribution update can take any form. This is a bit unfortunate for a "local"
               target distribution. One could avoid the global update by making a local target update function (and
               postponing target updates for non-local points as for the free energy update). We avoid such additions
               for now and accept that this target type always does global updates. */
            nstupdate_target = numSamplesUpdateFreeEnergy*numStepsSampleCoord;
            break;
        default:
            gmx_incons("Unknown AWH target type");
    }

    for (int d = 0; d < awhBiasParams.ndim; d++)
    {
        double coverRadiusInNm = 0.5*awhBiasParams.dim_params[d].coverDiameter;
        double spacing         = gridAxis[d].spacing;
        coverRadius[d]         = spacing > 0 ?  static_cast<int>(std::round(coverRadiusInNm/spacing)) : 0;
    }

    errorInitial    = awhBiasParams.error_initial;

    /* Estimate and initialize histSizeInitial. The estimation depends on the grid. */
    histSizeInitial = getInitialHistSizeEstimate(dimParams, awhBiasParams, gridAxis, numStepsSampleCoord*mdTimeStep);
}

DimParams::DimParams(double userCoordUnitsToInternalFactor,
                     double forceConstant,
                     double beta) :
    k(forceConstant),
    betak(beta*forceConstant),
    userCoordUnitsToInternal(userCoordUnitsToInternalFactor)
{
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
    biasIndex(biasIndexInCollection),
    dimParams_(dimParamsInit),
    grid_(new Grid(dimParamsInit, awhBiasParams.dim_params)),
    params_(awhParams, awhBiasParams, dimParams_, beta, mdTimeStep, disableUpdateSkips, cr, grid_->axis_),
    tempWorkSpace_()
{
    awh_dvec             coordValueInit;

    for (int d = 0; d < ndim(); d++)
    {
        coordValueInit[d] = dimParams_[d].scaleUserInputToInternal(awhBiasParams.dim_params[d].coord_value_init);
    }

    /* The initial sample weight is set to 1 and we keep the logarithm. */
    state_.scaledSampleWeight    = 0;
    state_.maxScaledSampleWeight = 0;

    state_.inInitialStage        = awhBiasParams.eGrowth == eawhgrowthEXP_LINEAR;
    state_.equilibrateHistogram  = awhBiasParams.equilibrateHistogram;
    state_.histSize              = params_.histSizeInitial;

    state_.numUpdates            = 0;

    /* Set initial coordinate reference value to the one closest to the initial reference value given in pull.
       More correctly one would sample from the biased distribution, but it doesn't really matter. */
    state_.gridpointIndex = getClosestIndexInGrid(coordValueInit, grid_.get());
    state_.refGridpoint   = state_.gridpointIndex;

    /* The minimum and maximum multidimensional point indices that are affected by the next update */
    for (int d = 0; d < ndim(); d++)
    {
        int index_d                = grid_->point(state_.gridpointIndex).index[d];
        state_.originUpdatelist[d] = index_d;
        state_.endUpdatelist[d]    = index_d;
    }

    /* Initialize free energy functions and biases */
    pointState_.resize(grid_->numPoints());
    updateList_.resize(grid_->numPoints());

    initGridPointState(*this, awhBiasParams, awhParams.nbias, params_.numSharedUpdate, (cr != NULL ? cr->ms : NULL),
                       &pointState_);

    /* Partition the AWH domain if any of the dimensions is set to be divided into more than 1 interval. */
    if (getNumDomains(awhBiasParams) > 1)
    {
        /* Partitioning  modifies the target distribution and the weighthistogram outside of the target domain.
           The target domain is determined based on the current coordinate reference value. */
        double targetSum = partitionDomain(awhBiasParams, *grid(), state_.gridpointIndex,
                                           &pointState_);
        /* Scale down the histogram size */
        state_.histSize         *= targetSum;
        params_.histSizeInitial  = state_.histSize;
    }

    /* Translate the point cover diameter fraction to a radius in units of number of points. */
    for (int d = 0; d < ndim(); d++)
    {
        double coverRadius     = 0.5*awhBiasParams.dim_params[d].coverDiameter;
        double spacing         = grid_->axis(d).spacing;
        params_.coverRadius[d] = spacing > 0 ?  static_cast<int>(std::round(coverRadius/spacing)) : 0;
    }

    bool blocklengthInWeight   = false;

    /* We let the correlation init function set its parameters to something useful for now. */
    double blockLength         = 0;

    if (MASTER(cr))
    {
        forceCorr_ = std::unique_ptr<CorrelationGrid>(new CorrelationGrid(pointState_.size(), ndim(),
                                                                         blockLength, blocklengthInWeight,
                                                                         params_.numStepsSampleCoord*mdTimeStep));

        writer_ = std::unique_ptr<BiasWriter>(new BiasWriter(*this));
    }

    /* Print information about AWH variables that are set internally but might be of interest to the user. */
    if ((cr == NULL) || (MASTER(cr)))
    {
        printInitializationToLog(*this, awhBiasParams, blocklengthInWeight, fplog);
    }
}
