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

/* Sets the given array with PMF values. */
void getPmf(const Bias &bias, const gmx_multisim_t *ms,
            std::vector<float> *pmf)
{
    /* The PMF is just the negative of the log of the sampled PMF histogram.
       Points with zero target weight are ignored, they will mostly contain noise. */

    const std::vector<CoordPoint> &coordPoint = bias.coordPoint;

    pmf->resize(coordPoint.size());

    if (bias.params.numSharedUpdate > 1)
    {
        /* Need to temporarily exponentiate the log weights to sum over simulations */
        for (size_t i = 0; i < coordPoint.size(); i++)
        {
            (*pmf)[i] = coordPoint[i].inTargetRegion() ? std::exp(static_cast<float>(coordPoint[i].log_pmfsum)) : 0;
        }

        gmx_sumf_sim(bias.coordPoint.size(), pmf->data(), ms);

        /* Take log again to get (non-normalized) PMF */
        for (size_t i = 0; i < coordPoint.size(); i++)
        {
            (*pmf)[i] = coordPoint[i].inTargetRegion() ? -std::log((*pmf)[i]) : GMX_FLOAT_MAX;
        }
    }
    else
    {
        for (size_t i = 0; i < coordPoint.size(); i++)
        {
            (*pmf)[i] = coordPoint[i].inTargetRegion() ? -coordPoint[i].log_pmfsum : GMX_FLOAT_MAX;
        }
    }
}

/*! \brief
 * Find the minimum free energy value.
 *
 * \param[in] bias      The AWH bias.
 * \returns the minimum free energy value.
 */
static double freeEnergyMinimumValue(const Bias &bias)
{
    double        f_min = GMX_DOUBLE_MAX;

    for (auto &coordPoint : bias.coordPoint)
    {
        if (coordPoint.inTargetRegion() && coordPoint.free_energy < f_min)
        {
            f_min = coordPoint.free_energy;
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
 * \param[in] bias          The AWH bias.
 * \param[in] point_index   Point to evaluate probability weight for.
 * \param[in] pointBias     Bias for the point (as a log weight).
 * \param[in] value         Coordinate value.
 * \returns the biased probability weight.
 */
static double biasedWeightFromPoint(const Bias &bias, int point_index, double pointBias, const awh_dvec value)
{
    double weight = 0;

    /* Only points in the target reigon have non-zero weight */
    if (bias.coordPoint[point_index].inTargetRegion())
    {
        double log_weight = pointBias;

        /* Add potential for all parameter dimensions */
        for (int d = 0; d < bias.ndim; d++)
        {
            double dev = get_deviation_from_point_along_gridaxis(bias.grid.get(), d, point_index, value[d]);
            log_weight -= 0.5*bias.dimParams[d].betak*dev*dev;
        }

        weight = std::exp(log_weight);
    }

    return weight;
}

/* Convolves the given PMF using the given AWH bias. */
void getConvolvedPmf(const Bias &bias, const gmx_multisim_t *ms,
                     std::vector<float> *convolvedPmf)
{
    size_t             numPoints = bias.grid.get()->numPoints();

    std::vector<float> pmf;

    convolvedPmf->resize(numPoints);

    /* Get the PMF to convolve. */
    getPmf(bias, ms, &pmf);

    for (size_t m = 0; m < numPoints; m++)
    {
        double           freeEnergyWeights = 0;
        const GridPoint &point             = bias.grid.get()->point(m);
        for (auto &neighbor : point.neighbor)
        {
            /* The negative PMF is a positive bias. */
            double biasNeighbor = -pmf[neighbor];

            /* Add the convolved PMF weights for the neighbors of this point.
               Note that this function only adds point within the target > 0 region.
               Sum weights, take the logarithm last to get the free energy. */
            freeEnergyWeights += biasedWeightFromPoint(bias, neighbor, biasNeighbor,
                                                       point.value);
        }

        GMX_RELEASE_ASSERT(freeEnergyWeights > 0, "Attempting to do log(<= 0) in AWH convolved PMF calculation.");
        (*convolvedPmf)[m] = -std::log(static_cast<float>(freeEnergyWeights));
    }
}

/*! \brief
 * Updates the target distribution of the AWH bias.
 *
 * The target distribution is always updated for all points
 * at the same time.
 *
 * \param[in,out] bias          The AWH bias.
 */
static void updateTarget(Bias *bias)
{
    const BiasParams &params = bias->params;
    double            f_max  = 0;

    if (bias->params.eTarget == eawhtargetCUTOFF)
    {
        f_max = freeEnergyMinimumValue(*bias) + params.target_param;
    }

    double sum_target = 0;
    for (auto &coordPoint : bias->coordPoint)
    {
        if (params.eTarget == eawhtargetCONSTANT)
        {
            coordPoint.target = 1;
        }
        else if (params.eTarget == eawhtargetCUTOFF)
        {
            double df =  coordPoint.free_energy - f_max;
            coordPoint.target = df > 0 ? std::exp(-df) : 1;
        }
        else if (params.eTarget == eawhtargetBOLTZMANN)
        {
            coordPoint.target = std::exp(-params.target_param*coordPoint.free_energy);
        }
        else if (params.eTarget == eawhtargetLOCALBOLTZMANN)
        {
            coordPoint.target = coordPoint.weightsum_ref;
        }

        /* All target types can be modulated by a constant factor. */
        coordPoint.target *= coordPoint.target_constant_weight;

        sum_target += coordPoint.target;
    }

    /* Normalize to 1 */
    double inv_sum = 1./sum_target;
    for (auto &coordPoint : bias->coordPoint)
    {
        coordPoint.target *= inv_sum;
    }
}

/*! \brief
 * Updates the bias of a coordinate point.
 *
 * \param[in,out] coordPoint  A coordinate point.
 */
static void updateBiasPoint(CoordPoint *coordPoint)
{
    GMX_RELEASE_ASSERT(coordPoint->target > 0, "AWH target distribution must be > 0 to calculate the point bias.");
    coordPoint->bias = coordPoint->free_energy + std::log(coordPoint->target);
}

/* Initialize the target and bias values for the given AWH bias. */
void initTargetAndBias(Bias *bias)
{

    /* The local Boltzmann distribution is special because the target distribution is updated as a function of the reference weighthistogram. */
    GMX_ASSERT((bias->params.eTarget != eawhtargetLOCALBOLTZMANN) ||
               (bias->params.eTarget == eawhtargetLOCALBOLTZMANN && bias->coordPoint[0].weightsum_ref != 0),
               "AWH reference weight histogram not initialized properly with local Boltzmann target distribution." );

    updateTarget(bias);

    for (auto &coordPoint : bias->coordPoint)
    {
        if (coordPoint.inTargetRegion())
        {
            updateBiasPoint(&coordPoint);
        }
        else
        {
            /* Note that for zero target this is a value that represents -infinity but should not be used for biasing. */
            coordPoint.bias = c_largeNegativeExponent;
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
 * \param[in] awh               AWH struct.
 * \param[in] biasForce         The AWH bias force.
 * \param[in,out] pull_work     Pull struct.
 * \param[in] mdatoms           Atom properties.
 * \param[in,out] force         Forces.
 * \param[in,out] virial        Virial tensor.
 */
static void applyBiasForceToPullCoords(const Bias       &Bias,
                                       const awh_dvec    biasForce,
                                       struct pull_t    *pull_work,
                                       const t_mdatoms  *mdatoms,
                                       rvec             *force,
                                       tensor            virial)
{
    for (int d = 0; d < Bias.ndim; d++)
    {
        apply_external_pull_coord_force(pull_work, Bias.dimParams[d].pullCoordIndex,
                                        biasForce[d], mdatoms, force, virial);
    }
}

/*! \brief
 * Puts together a string describing a grid point.
 *
 * \param[in] grid         The grid.
 * \param[in] point        Grid point index.
 * \param[in,out] pointstr String for the point.
 */
static void set_grid_point_value_string(const Grid *grid, int point, char *pointstr)
{
    char buf[STRLEN];

    strcpy(pointstr, "(");

    for (int d = 0; d < grid->ndim(); d++)
    {
        sprintf(buf, "%g", grid->point(point).value[d]);
        strcat(pointstr, buf);
        if (d < grid->ndim() - 1)
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
 * \param[in] bias         The AWH bias.
 * \param[in] t            Time.
 * \param[in] step         Time step.
 * \param[in,out] fplog    Output file for warnings. */
static void checkOnData(const Bias &bias, double t, gmx_int64_t step, FILE *fplog)
{
    const BiasParams  &params                 = bias.params;

    const int          max_nwarnings          = 1;   /* The maximum number of warnings to print per check */
    const int          max_alltime_nwarnings  = 5;   /* The maximum number of warnings to print ever */
    static int         alltime_nwarnings      = 0;   /* The number of warnings printed ever */
    const double       max_histogram_ratio    = 0.5; /* Tolerance for printing a warning about the histogram ratios */
    const int          nstcheck_on_data       = std::max(params.nstupdate_free_energy,
                                                         (static_cast<int>(bias.coordPoint.size())*100*params.nstsample_coord/params.nstupdate_free_energy)*params.nstupdate_free_energy);
    double             sum_W = 0, sum_visits = 0;
    double             inv_norm_visits, inv_norm_W;

    if (fplog == NULL || alltime_nwarnings >= max_alltime_nwarnings || bias.state.in_initial ||
        (step == 0) || (step % nstcheck_on_data != 0))
    {
        return;
    }

    /* Sum up the histograms and get their normalization */
    for (auto &coordPoint : bias.coordPoint)
    {
        if (coordPoint.inTargetRegion())
        {
            sum_visits += coordPoint.visits_tot;
            sum_W      += coordPoint.weightsum_tot;
        }
    }
    inv_norm_visits = 1./sum_visits;
    inv_norm_W      = 1./sum_W;

    /* Check all points for warnings */
    int                            nwarnings  = 0;
    const Grid                    *grid       = bias.grid.get();
    size_t                         numPoints  = grid->numPoints();
    const std::vector<CoordPoint> &coordPoint = bias.coordPoint;
    for (size_t m = 0; m < numPoints; m++)
    {
        bool bSkip_point = FALSE;

        /* Skip points close to boundary or non-target region */
        const GridPoint &gridPoint = grid->point(m);
        for (size_t n = 0; (n < gridPoint.neighbor.size()) && !bSkip_point; n++)
        {
            int ineighbor = gridPoint.neighbor[n];
            bSkip_point   = !coordPoint[ineighbor].inTargetRegion();
            for (int d = 0; (d < grid->ndim()) && !bSkip_point; d++)
            {
                const GridPoint &neighborPoint = grid->point(ineighbor);
                bSkip_point =
                    neighborPoint.index[d] == 0 ||
                    neighborPoint.index[d] == grid->axis(d).npoints - 1;
            }
        }

        /* Warn if the coordinate distribution less than the target distribution with a certain fraction somewhere */
        if (!bSkip_point &&
            (coordPoint[m].weightsum_tot*inv_norm_W*max_histogram_ratio > coordPoint[m].visits_tot*inv_norm_visits ))
        {
            char pointvaluestr[STRLEN], warningmsg[STRLEN];

            set_grid_point_value_string(grid, m, pointvaluestr);
            sprintf(warningmsg, "\nawh%d warning: "
                    "at t = %g ps the obtained coordinate distribution at coordinate value %s "
                    "is less than a fraction %g of the reference distribution at that point. "
                    "If you are not certain about your settings you might want to increase your pull force constant or "
                    "modify your sampling region.\n",
                    bias.biasIndex + 1, t, pointvaluestr, max_histogram_ratio);
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
        fprintf(fplog, "\nawh%d: suppressing future AWH warnings.\n", bias.biasIndex + 1);
    }
}

/*! \brief
 * Calculates and sets the force the coordinate experiences from an umbrella centered at the given point.
 *
 * The umbrella potential is an harmonic potential given by 0.5k(coord value - point value)^2. This
 * value is also returned.
 *
 * \param[in] bias         The AWH bias.
 * \param[in] point        Point for umbrella center.
 * \param[in,out] force    Force vector to set.
 * Returns the umbrella potential.
 */
static double calcUmbrellaForceAndPotential(const Bias &bias, int point, awh_dvec force)
{
    double potential = 0;
    for (int d = 0; d < bias.ndim; d++)
    {
        double dev = get_deviation_from_point_along_gridaxis(bias.grid.get(), d, point, bias.state.coordValue[d]);

        double k   = bias.dimParams[d].k;

        /* Force from harmonic potential 0.5*k*dev^2 */
        force[d]   = -k*dev;
        potential += 0.5*k*dev*dev;
    }

    return potential;
}

/*! \brief
 * Calculates and sets the convolved force acting on the coordinate.
 *
 * The convolved force is the weighted sum of forces from umbrellas
 * located at each point in the grid.
 *
 * \param[in] bias                The AWH bias.
 * \param[in] probWeightNeighbor  Probability weights of the neighbors.
 * \param[in,out] force           Bias force vector to set.
 */
static void calcConvolvedForce(const Bias                &bias,
                               const std::vector<double> &probWeightNeighbor,
                               awh_dvec                   force)
{
    for (int d = 0; d < bias.ndim; d++)
    {
        force[d] = 0;
    }

    /* Only neighboring points have non-negligible contribution. */
    const std::vector<int> &neighbor = bias.grid->point(bias.state.coordValueIndex).neighbor;
    for (size_t n = 0; n < neighbor.size(); n++)
    {
        double weightNeighbor = probWeightNeighbor[n];
        int    indexNeighbor  = neighbor[n];

        /* Get the umbrella force from this point. The returned potential is ignored here. */
        awh_dvec forceFromNeighbor;
        calcUmbrellaForceAndPotential(bias, indexNeighbor,
                                      forceFromNeighbor);

        /* Add the weighted umbrella force to the convolved force. */
        for (int d = 0; d < bias.ndim; d++)
        {
            force[d] += forceFromNeighbor[d]*weightNeighbor;
        }
    }
}

/*! \brief
 * Sample a new reference point given the current coordinate value.
 *
 * It is assumed that the probability distribution has been updated.
 *
 * \param[in] bias                The AWH bias.
 * \param[in] probWeightNeighbor  Probability weights of the neighbors.
 * \param[in] step                Step number, needed for the random number generator.
 * \param[in] seed                Random seed.
 * \returns the index of the sampled point.
 */
static int sampleReferenceCoordpoint(const Bias *bias,
                                     const std::vector<double> &probWeightNeighbor,
                                     gmx_int64_t step, int seed)
{
    /* Sample new reference value from the probability distribution which is defined for the neighboring
       points of the current coordinate value.*/

    /* In order to use the same seed for all AWH biases and get independent
       samples we use the index of the bias. */
    int n_sampled  = get_sample_from_distribution(probWeightNeighbor, bias->grid->point(bias->state.coordValueIndex).neighbor.size(),
                                                  step, seed, bias->biasIndex);

    return bias->grid->point(bias->state.coordValueIndex).neighbor[n_sampled];
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
 * \param[in,out] bias                The AWH bias.
 * \param[in]     probWeightNeighbor  Probability weights of the neighbors.

 * \param[in,out] biasForce           The AWH bias force.
 * \param[in] step                    Step number, needed for the random number generator.
 * \param[in] seed                    Random seed.
 * \returns the new potential value.
 */
double moveUmbrella(Bias *bias, const std::vector<double> &probWeightNeighbor,
                    awh_dvec biasForce,
                    gmx_int64_t step, int seed)
{
    /* Generate and set a new coordinate reference value */
    bias->state.refCoordpoint = sampleReferenceCoordpoint(bias, probWeightNeighbor, step, seed);

    awh_dvec newForce;
    double   newPotential =
        calcUmbrellaForceAndPotential(*bias, bias->state.refCoordpoint, newForce);

    /*  A modification of the reference value at time t will lead to a different
        force over t-dt/2 to t and over t to t+dt/2. For high switching rates
        this means the force and velocity will change signs roughly as often.
        To avoid any issues we take the average of the previous and new force
        at steps when the reference value has been moved. E.g. if the ref. value
        is set every step to (coord dvalue +/- delta) would give zero force.
     */
    for (int d = 0; d < bias->ndim; d++)
    {
        /* Average of the current and new force */
        biasForce[d] = 0.5*(biasForce[d] + newForce[d]);
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
 * \param[in,out] bias                 The AWH bias.
 * \param[in] histSizeNew              New reference weight histogram size.
 * \param[in] histSizeOld              Previous reference weight histogram size (before adding new samples).
 * \param[in,out] weighthistScaling    Scaling factor for the reference weight histogram.
 * \param[in,out] logPmfsumScaling     Log of the scaling factor for the PMF histogram.
 */
static void setHistogramUpdateScaleFactors(const Bias *bias, double histSizeNew, double histSizeOld,
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
    *weighthistScaling = histSizeNew/(histSizeOld + bias->params.update_weight*bias->params.localWeightScaling);
    *logPmfsumScaling  = std::log(histSizeNew/(histSizeOld + bias->params.update_weight));
}

/*! \brief
 * Sets the histogram rescaling factors needed for skipped updates.
 *
 * \param[in,out] bias                 The AWH bias.
 * \param[in,out] weighthistScaling    Scaling factor for the reference weight histogram.
 * \param[in,out] logPmfsumScaling     Log of the scaling factor for the PMF histogram.
 */
static void setSkippedUpdateHistogramScaleFactors(const Bias *bias, double *weighthistScaling, double *logPmfsumScaling)
{
    GMX_ASSERT(bias->params.canSkipUpdates(), "Calling function for skipped updates when skipping updates is not allowed");

    if (bias->state.in_initial)
    {
        /* In between global updates the reference histogram size is kept constant so we trivially know what the
            histogram size was at the time of the skipped update. */
        setHistogramUpdateScaleFactors(bias, bias->state.histSize, bias->state.histSize,
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
 * AWH parameters. This assumes that the variables used here are kept constant, at least in between
 * global updates.
 *
 * \param[in,out] coordPoint    Coordinate point to update free energy for.
 * \param[in] params            The AWH bias parameters.
 * \param[in] weight_point      Sampled probability weight at this point.
 */
static void update_free_energy_point(CoordPoint *coordPoint, const BiasParams &params, double weight_point)
{
    double weighthist_sampled       = coordPoint->weightsum_ref + weight_point;
    double weighthist_target        = coordPoint->weightsum_ref + params.update_weight*coordPoint->target;

    double df                       = -std::log(weighthist_sampled/weighthist_target);
    coordPoint->free_energy        += df;

    GMX_RELEASE_ASSERT(std::fabs(coordPoint->free_energy) < c_largePositiveExponent,
                       "Very large free energy differences or badly normalized free energy in AWH update.");
}

/*! \brief
 * Update the reference weight histogram of a point.
 *
 * \param[in,out] coordPoint    Coordinate point to update.
 * \param[in] params            The AWH bias parameters.
 * \param[in] weight_point      Sampled probability weight at this point.
 * \param[in] scaleFactor       Factor to rescale the histogram with.
 */
static void update_weighthistogram_point(CoordPoint *coordPoint, const BiasParams &params, double weight_point, double scaleFactor)
{
    if (params.idealWeighthistUpdate)
    {
        /* Grow histogram using the target distribution. */
        coordPoint->weightsum_ref += coordPoint->target*params.update_weight*params.localWeightScaling;
    }
    else
    {
        /* Grow using the actual samples (which are distributed ~ as target). */
        coordPoint->weightsum_ref += weight_point*params.localWeightScaling;
    }

    coordPoint->weightsum_ref *= scaleFactor;
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
 * \param[in,out] coordPoint       Coordinate point to update.
 * \param[in] params               The AWH bias parameters.
 * \param[in] weight_point         Sampled probability weight at this point.
 * \param[in] weighthistScaling    Scaling factor for the reference weight histogram.
 * \param[in] logPmfsumScaling     Log of the scaling factor for the PMF histogram.
 */
static void update_point(CoordPoint *coordPoint, const BiasParams &params,
                         double weight_point, double weighthistScaling, double logPmfsumScaling)
{
    update_free_energy_point(coordPoint, params, weight_point);
    update_weighthistogram_point(coordPoint, params, weight_point, weighthistScaling);
    coordPoint->log_pmfsum += logPmfsumScaling;
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
 * \param[in,out] coordPoint     Coordinate point to update.
 * \param[in] params             The AWH bias parameters.
 * \param[in] step               Current time step.
 * \param[in] weighthistScaling  Scale factor for the reference weight histogram.
 * \param[in] logPmfsumScaling      Scale factor for the reference PMF histogram.
 * \return true if at least one update was applied.
 */
static bool updateSkippedPoint(CoordPoint *coordPoint, const BiasParams &params,
                               gmx_int64_t step, double weighthistScaling, double logPmfsumScaling)
{
    int last_update_index, nupdates_skipped;

    GMX_ASSERT(params.canSkipUpdates(), "Calling function for skipped updates when skipping updates is not allowed");

    if (coordPoint->inTargetRegion())
    {
        return false;
    }

    /* The most current past update */
    last_update_index = step > 0 ? (step - 1)/params.nstupdate_free_energy : 0;
    nupdates_skipped  = last_update_index - coordPoint->last_update_index;

    if (nupdates_skipped == 0)
    {
        /* Was not updated */
        return false;
    }

    for (int i = 0; i < nupdates_skipped; i++)
    {
        /* This point was non-local at the time of the update meaning no weight */
        update_point(coordPoint, params, 0, weighthistScaling, logPmfsumScaling);
    }

    /* Only past updates are applied here. */
    coordPoint->last_update_index = last_update_index;

    /* Was updated */
    return true;
}

/* TODO: should probably merge do_skipped_updates_for_all_points and
   do_skipped_updates_in_neighborhood into one function that takes as
   input a vector with points to be updated. */

/*! \brief
 * Do all previously skipped updates.
 *
 * \param[in,out] bias              The AWH bias.
 * \param[in] step                  Current time step.
 */
static void doSkippedUpdatesForAllPoints(Bias *bias, gmx_int64_t step)
{
    double weighthistScaling, logPmfsumScaling;

    setSkippedUpdateHistogramScaleFactors(bias, &weighthistScaling, &logPmfsumScaling);

    for (auto &coordPoint : bias->coordPoint)
    {
        bool bUpdated = updateSkippedPoint(&coordPoint, bias->params, step, weighthistScaling, logPmfsumScaling);

        /* Update the bias for this point only if there were skipped updates in the past to avoid calculating the log unneccessarily */
        if (bUpdated)
        {
            updateBiasPoint(&coordPoint);
        }
    }
}

/*! \brief
 * Do previously skipped updates in this neighborhood.
 *
 * \param[in,out] bias              The AWH bias.
 * \param[in] step                  Current time step.
 */
static void doSkippedUpdatesInNeighborhood(Bias *bias, gmx_int64_t step)
{
    double weighthistScaling, logPmfsumScaling;

    setSkippedUpdateHistogramScaleFactors(bias, &weighthistScaling, &logPmfsumScaling);

    /* For each neighbor point of the center point, refresh its state by adding the results of all past, skipped updates. */
    const std::vector<int> &neighbors = bias->grid->point(bias->state.coordValueIndex).neighbor;
    for (auto &neighbor : neighbors)
    {
        bool bUpdated = updateSkippedPoint(&bias->coordPoint[neighbor], bias->params, step, weighthistScaling, logPmfsumScaling);

        if (bUpdated)
        {
            updateBiasPoint(&bias->coordPoint[neighbor]);
        }
    }
}

/*! \brief
 * Merge update lists from multiple sharing simulations.
 *
 * \param[in,out] updateList        Update list for this simulation (assumed >= npoints long).
 * \param[in,out] numUpdatePtr      Number of points in the update list.
 * \param[in] npoints               Total number of points.
 * \param[in] ms                    Struct for multi-simulation communication.
 */
static void merge_shared_update_lists(std::vector<int> *updateList, int *numUpdatePtr, int npoints, const gmx_multisim_t *ms)
{
    std::vector<int> nupdate_point;

    /* Flag the update points of this sim.
       TODO: we can probably avoid allocating this array and just use the input array. */
    nupdate_point.resize(npoints, 0);
    for (int i = 0; i < *numUpdatePtr; i++)
    {
        nupdate_point[(*updateList)[i]] = 1;
    }

    /* Sum over the sims to get all the flagged points */
    gmx_sumi_sim(npoints, nupdate_point.data(), ms);

    /* Collect the indices of the flagged points in place. The resulting array will be the merged update list.*/
    int numUpdateMerged = 0;
    for (int m = 0; m < npoints; m++)
    {
        if (nupdate_point[m] > 0)
        {
            (*updateList)[numUpdateMerged++] = m;
        }
    }

    *numUpdatePtr = numUpdateMerged;
}

/*! \brief
 * Generate an update list of points sampled since the last update.
 *
 * \param[in] bias                  The AWH bias.
 * \param[in] ms                    Struct for multi-simulation communication.
 * \param[in,out] updateList        Local update list to set (assumed >= npoints long).
 * \param[in,out] nupdate_ptr       Number of points in the update list.
 */
static void make_local_update_list(const Bias *bias, const gmx_multisim_t *ms, std::vector<int> *updateList, int *nupdate_ptr)
{
    awh_ivec origin_dim, npoints_dim;

    /* Define the update search grid */
    for (int d = 0; d < bias->grid->ndim(); d++)
    {
        origin_dim[d]  = bias->state.origin_updatelist[d];
        npoints_dim[d] = bias->state.end_updatelist[d] - bias->state.origin_updatelist[d] + 1;

        /* Because the end_updatelist is unwrapped it can be > (npoints - 1) so that npoints_dim can be > npoints in grid.
           This helps for calculating the distance/number of points but should be removed and fixed when the way of
           updating origin/end updatelist is changed (see sampleProbabilityWeights). */
        npoints_dim[d] = std::min(bias->grid->axis(d).npoints, npoints_dim[d]);
    }

    /* Make the update list */
    int  point_index   = -1;
    int  nupdate       = 0;
    bool pointExists   = true;
    while (pointExists)
    {
        pointExists = get_next_point_in_local_grid(bias->grid.get(), origin_dim, npoints_dim, &point_index);

        if (pointExists && bias->coordPoint[point_index].inTargetRegion())
        {
            (*updateList)[nupdate] = point_index;
            nupdate++;
        }
    }

    if (bias->params.numSharedUpdate > 1)
    {
        merge_shared_update_lists(updateList, &nupdate, bias->coordPoint.size(), ms);
    }

    /* Return the number of points to update */
    *nupdate_ptr = nupdate;
}

/*! \brief
 * Reset the range used to make the local update list.
 *
 * \param[in,out] bias  The AWH bias.
 */
static void reset_local_update_range(Bias *bias)
{
    BiasState *state = &bias->state;

    for (int d = 0; d < bias->grid->ndim(); d++)
    {
        /* This gives the  minimum range consisting only of the current closest point. */
        state->origin_updatelist[d] = bias->grid->point(state->coordValueIndex).index[d];
        state->end_updatelist[d]    = bias->grid->point(state->coordValueIndex).index[d];
    }
}

/*! \brief
 * Apply a point update with new sampling data.
 *
 * The last update index is also updated here.
 *
 * \param[in,out] coordPoint      Coordinate point to update.
 * \param[in] params              The AWH bias parameters.
 * \param[in] step                Time step for updating the last update index.
 * \param[in] weighthistScaling   Scaling factor for the reference weight histogram.
 * \param[in] logPmfsumScaling    Log of the scaling factor for the PMF histogram.
 */
static void updatePointNew(CoordPoint *coordPoint, const BiasParams &params,
                           gmx_int64_t step, double weighthistScaling, double logPmfsumScaling)
{
    update_point(coordPoint, params, coordPoint->weightsum_iteration, weighthistScaling, logPmfsumScaling);

    /* Remember that this was most recent update performed on this point. */
    coordPoint->last_update_index   = step/params.nstupdate_free_energy;
}

/*! \brief
 * Add partial histograms (accumulating between updates) to accumulating histograms.
 *
 * \param[in] bias               The AWH bias.
 * \param[in] ms                 Struct for multi-simulation communication.
 * \param[in] local_update_list  List of points with data.
 * \param[in] nlocal_update      Number of points in the list.
 */
static void sum_histograms(Bias *bias, const gmx_multisim_t *ms, const std::vector<int> &local_update_list, int nlocal_update)
{
    /* The covering checking histograms are added before summing over simulations, so that the weights from different
       simulations are kept distinguishable. */
    for (int ilocal = 0; ilocal < nlocal_update; ilocal++)
    {
        int    iglobal = local_update_list[ilocal];
        bias->coordPoint[iglobal].weightsum_covering +=  bias->coordPoint[iglobal].weightsum_iteration;
    }

    /* Sum histograms over multiple simulations if needed. */
    if (bias->params.numSharedUpdate > 1)
    {
        /* Collect the weights and counts in linear arrays to be able to use gmx_sumd_sim. */
        std::vector<double> weight_distr, coord_visits;

        weight_distr.resize(nlocal_update);
        coord_visits.resize(nlocal_update);

        for (int ilocal = 0; ilocal < nlocal_update; ilocal++)
        {
            int iglobal = local_update_list[ilocal];

            weight_distr[ilocal] = bias->coordPoint[iglobal].weightsum_iteration;
            coord_visits[ilocal] = bias->coordPoint[iglobal].visits_iteration;
        }

        gmx_sumd_sim(nlocal_update, weight_distr.data(), ms);
        gmx_sumd_sim(nlocal_update, coord_visits.data(), ms);

        /* Transfer back the result */
        for (int ilocal = 0; ilocal < nlocal_update; ilocal++)
        {
            int iglobal = local_update_list[ilocal];

            bias->coordPoint[iglobal].weightsum_iteration = weight_distr[ilocal];
            bias->coordPoint[iglobal].visits_iteration    = coord_visits[ilocal];
        }
    }

    /* Now add the partial counts and weights to the accumulating histograms.
       Note: we still need to use the weights for the update so we wait
       with resetting them until the end of the update. */
    for (int ilocal = 0; ilocal < nlocal_update; ilocal++)
    {
        int    iglobal = local_update_list[ilocal];

        double weights = bias->coordPoint[iglobal].weightsum_iteration;
        double visits  = bias->coordPoint[iglobal].visits_iteration;
        bias->coordPoint[iglobal].weightsum_tot      += weights;
        bias->coordPoint[iglobal].visits_tot         += visits;
    }
}

/*! \brief
 * Returns the new size of the reference weight histogram in the initial stage.
 *
 * This function also takes care resetting the histogram used for covering checks
 * and for exiting the initial stage.
 *
 * \param[in,out] bias   The AWH bias.
 * \param[in] t          Time.
 * \param[in] bCovered   True if the sampling interval has been covered enough.
 * \param[in,out] fplog  Log file.
 * \returns the new histogram size.
 */
static double newHistSizeInitialStage(Bias *bias, double t, bool bCovered, FILE *fplog)
{
    BiasState *state = &bias->state;

    /* The histogram size is kept constant until the sampling region has been covered
       and the the current sample weight is large enough and the histogram is ready. */
    if (!bCovered ||
        (state->scaledSampleWeight < state->maxScaledSampleWeight) ||
        state->equilibrateHistogram)
    {
        return state->histSize;
    }

    /* Reset the covering weight histogram. If we got this far we are either entering a
       new covering stage with a new covering histogram or exiting the initial stage
       altogether. */
    for (auto &coordPoint : bias->coordPoint)
    {
        coordPoint.weightsum_covering = 0.;
    }

    /*  The current sample weigth is now the maximum. */
    double prevMaxScaledSampleWeight = state->maxScaledSampleWeight;
    state->maxScaledSampleWeight = state->scaledSampleWeight;

    /* Increase the histogram size by a constant scale factor if we can, i.e. if the sample weight
       resulting from such a scaling is still larger than the previous maximum sample weight
       (ensuring that the sample weights at the end of each covering stage are monotonically
       increasing). If we cannot, exit the initial stage without changing the histogram size. */

    /* The scale factor. The value is not very critical but should obviously be > 1 (or the exit
       will happen very late) and probably < 5 or so (or there will be no initial stage). */
    static const double growthFactor = 3;

    /* The scale factor is in most cases very close to the histogram growth factor. */
    double              scaleFactor = growthFactor/(1. + bias->params.update_weight*bias->params.localWeightScaling/state->histSize);

    bool                bExit       = (state->scaledSampleWeight - std::log(scaleFactor) <= prevMaxScaledSampleWeight);
    double              histSizeNew = bExit ? state->histSize : state->histSize*growthFactor;

    /* Update the AWH bias about the exit. */
    state->in_initial = !bExit;

    /* Print information about coverings and if there was an exit. */
    if (fplog != NULL)
    {
        char                buf[STRLEN];
        sprintf(buf, "\nawh%d:", bias->biasIndex + 1);
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
    return histSizeNew;
}

/*! \brief
 * Label points along an axis as covered or not.
 *
 * A point is covered if it is surrounded by visited points up to a radius = coverRadius.
 *
 * \param[in] visited           Visited? For each point.
 * \param[in] checkCovering     Check for covering? For each point.
 * \param[in] npoints           Number of points.
 * \param[in] period            Period in number of points.
 * \param[in] coverRadius       Cover radius, in points, needed for defining a point as covered.
 * \param[in,out] covered       Covered? For each point.
 */
static void labelCoveredPoints(const std::vector<bool> &visited, const std::vector<bool> &checkCovering, int npoints, int period,
                               int coverRadius, std::vector<int> *covered)
{

    bool haveFirstNotVisited = false;
    int  firstNotVisited     = -1, notVisitedLow = -1, notVisitedHigh = -1;

    for (int n = 0; n < npoints; n++)
    {
        if (checkCovering[n] && !visited[n])
        {
            if (!haveFirstNotVisited)
            {
                notVisitedLow       = n;
                firstNotVisited     = n;
                haveFirstNotVisited = true;
            }
            else
            {
                notVisitedHigh = n;

                /* Have now an interval I = [notVisitedLow,notVisitedHigh] of visited points bounded
                   by unvisited points. The unvisted end points affect the coveredness of the visited
                   with a reach equal to the cover radius. */
                int notCoveredLow  = notVisitedLow + coverRadius;
                int notCoveredHigh = notVisitedHigh - coverRadius;
                for (int i = notVisitedLow; i <= notVisitedHigh; i++)
                {
                    (*covered)[i] = (i > notCoveredLow) && (i < notCoveredHigh);
                }

                /* Find a new interval to set covering for. Make the notVisitedHigh of this interval the
                   notVisitedLow of the next. */
                notVisitedLow = notVisitedHigh;
            }
        }
    }

    /* Have labelled all the internal points. Now take care of the boundary regions. */
    if (!haveFirstNotVisited)
    {
        /* No non-visited points <=> all points visited => all points covered. */

        for (int n = 0; n < npoints; n++)
        {
            (*covered)[n] = TRUE;
        }
    }
    else
    {
        int lastNotVisited = notVisitedLow;

        /* For periodic boundaries, non-visited points can influence points
           on the other side of the boundary so we need to wrap around. */

        /* Lower end. For periodic boundaries the last upper end not visited point becomes the low-end not visited point.
           For non-periodic boundaries there is no lower end point so a dummy value is used. */
        int notVisitedHigh = firstNotVisited;
        int notVisitedLow  = period > 0 ? (lastNotVisited - period) : -(coverRadius + 1);

        int notCoveredLow  =  notVisitedLow + coverRadius;
        int notCoveredHigh = notVisitedHigh - coverRadius;

        for (int i = 0; i <= notVisitedHigh; i++)
        {
            /* For non-periodic boundaries notCoveredLow = -1 will impose no restriction. */
            (*covered)[i] = (i > notCoveredLow) && (i < notCoveredHigh);
        }

        /* Upper end. Same as for lower end but in the other direction. */
        notVisitedHigh = period > 0 ? (firstNotVisited + period) : (npoints + coverRadius);
        notVisitedLow  = lastNotVisited;

        notCoveredLow  = notVisitedLow + coverRadius;
        notCoveredHigh = notVisitedHigh - coverRadius;

        for (int i = notVisitedLow; i <= npoints - 1; i++)
        {
            /* For non-periodic boundaries notCoveredHigh = npoints will impose no restriction. */
            (*covered)[i] = (i > notCoveredLow) && (i < notCoveredHigh);
        }
    }
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
 * \param[in] bias  The AWH bias.
 * \param[in] ms    Struct for multi-simulation communication.
 * \returns true if covered.
 */
static bool isCovered(const Bias &bias, const gmx_multisim_t *ms)
{
    const Grid   *grid = bias.grid.get();

    /* Allocate and initialize arrays: one for checking visits along each dimension,
       one for keeping track of which points to check and one for the covered points.
       Possibly these could be kept as AWH variables to avoid these allocations. */
    struct CheckDim
    {
        std::vector<bool> visited;
        std::vector<bool> checkCovering;
        // We use int for the covering array since we might use gmx_sumi_sim.
        std::vector<int>  covered;
    };

    std::vector<CheckDim> checkDim;
    checkDim.resize(grid->ndim());

    for (int d = 0; d < grid->ndim(); d++)
    {
        size_t npoints = grid->axis(d).npoints;
        checkDim[d].visited.resize(npoints, false);
        checkDim[d].checkCovering.resize(npoints, false);
        checkDim[d].covered.resize(npoints, 0);
    }

    /* Set visited points along each dimension and which points should be checked for covering.
       Specifically, points above the free energy cutoff (if there is one) or points outside
       of the target region are ignored. */

    /* Set the free energy cutoff */
    double maxFreeEnergy = GMX_DOUBLE_MAX;

    if (bias.params.eTarget == eawhtargetCUTOFF)
    {
        maxFreeEnergy = freeEnergyMinimumValue(bias) + bias.params.target_param;
    }

    /* Set the cutoff weight for a point to be considered visited. */
    double weight_peak = 1;
    for (int d = 0; d < bias.ndim; d++)
    {
        weight_peak *= grid->axis(d).spacing*std::sqrt(bias.dimParams[d].betak*0.5*M_1_PI);
    }

    /* Project the sampling weights onto each dimension */
    for (size_t m = 0; m < grid->numPoints(); m++)
    {
        const CoordPoint &coordPoint = bias.coordPoint[m];

        for (int d = 0; d < grid->ndim(); d++)
        {
            int n = grid->point(m).index[d];

            /* Is visited if it was already visited or if there is enough weight at the current point */
            checkDim[d].visited[n] = checkDim[d].visited[n] || (coordPoint.weightsum_covering > weight_peak);

            /* Check for covering if there is at least point in this slice that is in the target region and within the cutoff */
            checkDim[d].checkCovering[n] = checkDim[d].checkCovering[n] || (coordPoint.inTargetRegion() && coordPoint.free_energy < maxFreeEnergy);
        }
    }

    /* Label each point along each dimension as covered or not. */
    for (int d = 0; d < grid->ndim(); d++)
    {
        labelCoveredPoints(checkDim[d].visited, checkDim[d].checkCovering, grid->axis(d).npoints, grid->axis(d).npoints_period,
                           bias.params.coverRadius[d], &checkDim[d].covered);
    }

    /* Now check for global covering. Each dimension needs to be covered separately.
       A dimension is covered if each point is covered.  Multiple simulations collectively
       cover the points, i.e. a point is covered if any of the simulations covered it.
       However, visited points are not shared, i.e. if a point is covered or not is
       determined by the visits of a single simulation. In general the covering criterion is
       all points covered => all points are surrounded by visited points up to a radius = coverRadius.
       For 1 simulation, all points covered <=> all points visited. For multiple simulations
       however, all points visited collectively !=> all points covered, except for coverRadius = 0.
       In the limit of large coverRadius, all points covered => all points visited by at least one
       simulation (since no point will be covered until all points have been visited by a
       single simulation). Basically coverRadius sets how much "connectedness" (or mixing) a point
       needs with surrounding points before sharing covering information with other simulations. */

    /* Communicate the covered points between sharing simulations if needed. */
    if (bias.params.numSharedUpdate > 1)
    {
        /* For multiple dimensions this may not be the best way to do it. */
        for (int d = 0; d < grid->ndim(); d++)
        {
            gmx_sumi_sim(grid->axis(d).npoints, checkDim[d].covered.data(), ms);
        }
    }

    /* Now check if for each dimension all points are covered. Break if not true. */
    bool allPointsCovered = true;
    for (int d = 0; d < grid->ndim() && allPointsCovered; d++)
    {
        for (int n = 0; n < grid->axis(d).npoints && allPointsCovered; n++)
        {
            allPointsCovered = checkDim[d].covered[n];
        }
    }

    return allPointsCovered;
}

/*! \brief
 * Checks if the histogram has equilibrated to the target distribution.
 *
 * The histogram is considered equilibrated if, for a minimum fraction of
 * the target region, the relative error of the sampled weight relative
 * to the target is less than a tolerance value.
 *
 * \param[in] bias  The AWH bias.
 * \returns true if the histogram is equilibrated.
 */
static bool histogramIsEquilibrated(const Bias &bias)
{
    if (!bias.state.equilibrateHistogram)
    {
        return true;
    }

    /* Get the total weight of the total weight histogram; needed for normalization. */
    double totalWeight     = 0;
    int    numTargetPoints = 0;
    for (auto &coordPoint : bias.coordPoint)
    {
        if (!coordPoint.inTargetRegion())
        {
            continue;
        }
        totalWeight += coordPoint.weightsum_tot;
        numTargetPoints++;
    }
    GMX_ASSERT(totalWeight > 0, "No samples when normalizing AWH histogram.");
    double              inverseTotalWeight   = 1./totalWeight;

    /* Points with target weight below a certain cutoff are ignored. */
    static const double minTargetCutoff  = 0.05;
    double              minTargetWeight  = 1./numTargetPoints*minTargetCutoff;

    /* Points with error less than this tolerance pass the check.*/
    static const double errorTolerance   = 0.2;

    /* Sum up weight of points that do or don't pass the check. */
    double equilibratedWeight    = 0;
    double notEquilibratedWeight = 0;
    for (auto &coordPoint : bias.coordPoint)
    {
        double targetWeight  = coordPoint.target;
        double sampledWeight = coordPoint.weightsum_tot*inverseTotalWeight;

        /* Ignore these points. */
        if (!coordPoint.inTargetRegion() || targetWeight < minTargetWeight)
        {
            continue;
        }

        if (std::fabs(sampledWeight/targetWeight - 1) > errorTolerance)
        {
            notEquilibratedWeight += targetWeight;
        }
        else
        {
            equilibratedWeight += targetWeight;
        }
    }

    /* It is enough if sampling in at least a fraction of the target region follows the target
       distribution. Boundaries will in general fail and this should be ignored (to some extent). */
    static const double minFraction = 0.8;

    return equilibratedWeight/(equilibratedWeight + notEquilibratedWeight) > minFraction;;
}

/*! \brief
 * Return the new reference weight histogram size for the current update.
 *
 * This function also takes care of checking for covering in the initial stage.
 *
 * \param[in,out] bias   The AWH bias.
 * \param[in] t          Time.
 * \param[in] covered    True if the sampling interval has been covered enough.
 * \param[in,out] fplog  Log file.
 * \returns the new histogram size.
 */
static double newHistSize(Bias *bias, double t, bool covered, FILE *fplog)
{
    BiasState *state = &bias->state;

    double     histSizeNew;
    if (state->in_initial)
    {
        /* Only bother with checking equilibration if we have covered already. */
        if (state->equilibrateHistogram && covered)
        {
            /* The histogram is equilibrated at most once. */
            state->equilibrateHistogram = !histogramIsEquilibrated(*bias);

            static bool printAboutCovering = true;
            if (!state->equilibrateHistogram)
            {
                char buf[STRLEN];
                sprintf(buf, "\nawh%d:", bias->biasIndex + 1);
                fprintf(fplog, "%s equilibrated histogram at t = %g ps.\n", buf, t);
            }
            else if (printAboutCovering)
            {
                char buf[STRLEN];
                sprintf(buf, "\nawh%d:", bias->biasIndex + 1);
                fprintf(fplog, "%s covered but histogram not equilibrated at t = %g ps.\n", buf, t);
                printAboutCovering = false; /* Just print once. */
            }
        }

        /* In the initial stage, the histogram grows dynamically as a function of the number of coverings. */
        histSizeNew = newHistSizeInitialStage(bias, t, covered, fplog);
    }
    else
    {
        /* If not in the initial stage, the histogram grows at a linear, possibly scaled down, rate. */
        histSizeNew = state->histSize + bias->params.update_weight*bias->params.localWeightScaling;
    }

    return histSizeNew;
}

/*! \brief
 * Normalizes the free energy and PMF sum.
 *
 * The choice of free energy normalization is quite arbitrary. Here, the
 * minimum free energy is set to 0. The normalization of the free energy
 * and PMF sum are not independent because the PMF is sampled by reweighting
 * with the bias = free energy + log(target).
 *
 * \param[in,out] bias  The AWH bias.
 */
static void normalizeFreeEnergyAndPmfSum(Bias *bias)
{
    double minF = freeEnergyMinimumValue(*bias);

    for (auto &coordPoint : bias->coordPoint)
    {
        if (coordPoint.inTargetRegion())
        {
            /* The sign of the free energy and PMF constants are opposite
               because the PMF samples are reweighted with the negative
               bias e^(-bias) ~ e^(-free energy). */
            coordPoint.free_energy -= minF;
            coordPoint.log_pmfsum  += minF;
        }
    }
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
 * \param[in,out] bias   The AWH bias.
 * \param[in] ms         Struct for multi-simulation communication.
 * \param[in] t          Time.
 * \param[in] step       Time step.
 * \param[in,out] fplog  Log file.
 */
static void update(Bias *bias, const gmx_multisim_t *ms, double t, gmx_int64_t step, FILE *fplog)
{
    /* This array is only used in this scope and is always re-initialized */
    std::vector<int> &updateList = bias->updateList;

    /* Make a list of all local points, i.e. those that could have been touched since
       the last update. These are the points needed for summing histograms below
       (non-local points only add zeros). For local updates, this will also be the
       final update list. */
    int      nlocal_update;
    make_local_update_list(bias, ms, &updateList, &nlocal_update);

    /* Reset the range for the next update */
    reset_local_update_range(bias);

    /* Add samples to histograms for all local points and sync simulations if needed */
    sum_histograms(bias, ms, updateList, nlocal_update);

    /* Renormalize the free energy if values are too large. */
    bool needToNormalizeFreeEnergy = false;
    for (int ilocal = 0; ilocal < nlocal_update; ilocal++)
    {
        int    iglobal = updateList[ilocal];
        /* We want to keep the absolute value of the free energies to be less c_largePositiveExponent
           to be able to safely pass these values to exp(). The check below ensures this as long as
           the free energy values grow less than 0.5*c_largePositiveExponent in a return time to this
           neighborhood. For reasonable update sizes it's unlikely that this requirement would be
           broken. */
        if (std::fabs(bias->coordPoint[iglobal].free_energy) > 0.5*c_largePositiveExponent)
        {
            needToNormalizeFreeEnergy = true;
            break;
        }
    }

    const BiasParams &params = bias->params;

    /* Update target distribution? */
    bool doUpdateTarget = (params.eTarget != eawhtargetCONSTANT &&
                           do_at_step(params.nstupdate_target, step));


    bool haveCovered  = false;
    if (bias->state.in_initial)
    {
        /* In the initial stage, the histogram grows dynamically as a function of the number of coverings. */
        int nstcheck_covered = std::max((static_cast<int>(bias->coordPoint.size())*params.nstsample_coord/params.nstupdate_free_energy)*params.nstupdate_free_energy,
                                        params.nstupdate_free_energy);
        haveCovered          = (step % nstcheck_covered == 0) && isCovered(*bias, ms);
    }

    /* The weighthistogram size after this update. */
    double histSizeNew = newHistSize(bias, t, haveCovered, fplog);

    /* Make the update list. Usually we try to only update local points but if the update has
       non-trivial or non-deterministic affect on non-local points a global update is needed.
       This is the case when: 1) a covering occurred in the initial stage, leading to non-trivial
       histogram rescaling factors; or 2) the target distribution will be updated since we don't
       make any assumption on its form; or 3) the AWH parameters are such that we never attempt
       to skip non-local updates; or 4) the free energy values have grown so large that the
       a renormalization is needed. */
    int nupdate;
    if (doUpdateTarget || haveCovered || !params.canSkipUpdates() || needToNormalizeFreeEnergy)
    {
        /* Global update, just add all points. */
        nupdate = 0;
        for (size_t m = 0; m < bias->coordPoint.size(); m++)
        {
            if (bias->coordPoint[m].inTargetRegion())
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
    double   weighthistScalingSkipped = 0, logPmfsumScalingSkipped = 0;
    double   weighthistScalingNew, logPmfsumScalingNew;
    if (params.canSkipUpdates())
    {
        setSkippedUpdateHistogramScaleFactors(bias, &weighthistScalingSkipped, &logPmfsumScalingSkipped);
    }
    setHistogramUpdateScaleFactors(bias, histSizeNew, bias->state.histSize,
                                   &weighthistScalingNew, &logPmfsumScalingNew);

    /* Update free energy and reference weight histogram for points in the update list. */
    for (int iupdate = 0; iupdate < nupdate; iupdate++)
    {
        CoordPoint *coordPointToUpdate = &bias->coordPoint[updateList[iupdate]];

        /* Do updates from previous update steps that were skipped because this point was at that time non-local. */
        if (params.canSkipUpdates())
        {
            updateSkippedPoint(coordPointToUpdate, params, step, weighthistScalingSkipped, logPmfsumScalingSkipped);
        }

        /* Now do an update with new sampling data. */
        updatePointNew(coordPointToUpdate, params, step, weighthistScalingNew, logPmfsumScalingNew);

        /* Reset histograms collected for this update. */
        coordPointToUpdate->weightsum_iteration = 0;
        coordPointToUpdate->visits_iteration    = 0;
    }

    /* Only update the histogram size after we are done with the local point updates */
    bias->state.histSize = histSizeNew;

    /* The weight of new samples relative to previous ones change when the histogram is
       rescaled. We keep the log since this number can become very large. */
    bias->state.scaledSampleWeight -= std::log(weighthistScalingNew);

    if (needToNormalizeFreeEnergy)
    {
        normalizeFreeEnergyAndPmfSum(bias);
    }

    if (doUpdateTarget)
    {
        /* The target distribution is always updated for all points at once. */
        updateTarget(bias);
    }

    /* Update the bias. The bias is updated separately and last since it simply a function of
       the free energy and the target distribution and we want to avoid doing extra work. */
    for (int iupdate = 0; iupdate < nupdate; iupdate++)
    {
        updateBiasPoint(&bias->coordPoint[updateList[iupdate]]);
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
 * \param[in,out] bias           The AWH bias.
 * \param[out]    weight         Probability weights of the neighbors.
 * \param[out]    convolvedBias  The convolved bias.
 */
static void updateProbabilityWeightsAndConvolvedBias(Bias                *bias,
                                                     std::vector<double> *weight,
                                                     double              *convolvedBias)
{
    const Grid *grid = bias->grid.get();

    /* Only neighbors of the current coordinate value will have a non-negligible chance of getting sampled */
    const std::vector<int> &neighbors = grid->point(bias->state.coordValueIndex).neighbor;

    weight->resize(neighbors.size());

    double weightSum = 0;
    for (size_t n = 0; n < neighbors.size(); n++)
    {
        int neighbor = neighbors[n];
        (*weight)[n] = biasedWeightFromPoint(*bias, neighbor, bias->coordPoint[neighbor].bias,
                                             bias->state.coordValue);
        weightSum   += (*weight)[n];
    }
    GMX_RELEASE_ASSERT(weightSum > 0, "zero probability weight when updating AWH probability weights.");

    /* Normalize probabilities to sum to 1 */
    double invWeightSum = 1/weightSum;
    for (auto &w : *weight)
    {
        w *= invWeightSum;
    }

    *convolvedBias = std::log(weightSum);
}

/* Calculates the convolved bias for the given coordinate value. */
double calcConvolvedBias(const Bias &bias, const awh_dvec coordValue)
{
    int                            point      = get_closest_index_in_grid(coordValue, bias.grid.get());
    const GridPoint               &gridPoint  = bias.grid->point(point);
    const std::vector<CoordPoint> &coordPoint = bias.coordPoint;

    /* Sum the probability weights from the neighborhood of the given point */
    double weight_sum = 0;
    for (auto &neighbor : gridPoint.neighbor)
    {
        weight_sum += biasedWeightFromPoint(bias, neighbor, coordPoint[neighbor].bias,
                                            coordValue);
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
 * \param[in,out] bias            The AWH bias.
 * \param[in] probWeightNeighbor  Probability weights of the neighbors.
 */
static void sampleProbabilityWeights(Bias                      *bias,
                                     const std::vector<double> &probWeightNeighbor)
{
    const std::vector<int> &neighbor = bias->grid->point(bias->state.coordValueIndex).neighbor;

    /* Save weights for next update */
    for (size_t n = 0; n < neighbor.size(); n++)
    {
        bias->coordPoint[neighbor[n]].weightsum_iteration += probWeightNeighbor[n];
    }

    /* Update the local update range. Two corner points define this rectangular domain. We need to
       choose two new corner points such that the new domain contains both the old update range and
       the current neighborhood. In the simplest case when an update is performed every sample, the
        update range would simply equal the current neighborhood. */
    int m_neighbor_origin = neighbor[0];
    int m_neighbor_end    = neighbor[neighbor.size() - 1];
    for (int d = 0; d < bias->grid->ndim(); d++)
    {
        int origin_d, end_d;

        origin_d = bias->grid->point(m_neighbor_origin).index[d];
        end_d    = bias->grid->point(m_neighbor_end).index[d];

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
            end_d += bias->grid->axis(d).npoints_period;
        }

        bias->state.origin_updatelist[d] = std::min(bias->state.origin_updatelist[d], origin_d);
        bias->state.end_updatelist[d]    = std::max(bias->state.end_updatelist[d], end_d);
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
 * \param[in,out] bias           The AWH bias.
 * \param[in]     convolvedBias  The convolved bias.
 */
static void sample_pmf(Bias *bias, double convolvedBias)
{
    /* Only save coordinate data that is in range (the given index is always
       in range even if the coordinate value is not). */
    bool        bIn_array  = value_is_in_grid(bias->state.coordValue, bias->grid.get());
    CoordPoint *coordPoint = &bias->coordPoint[bias->state.coordValueIndex];

    /* Save PMF sum and keep a histogram of the sampled coordinate values */
    if (bIn_array && coordPoint->inTargetRegion())
    {
        coordPoint->log_pmfsum        = expsum(coordPoint->log_pmfsum, -convolvedBias);
        coordPoint->visits_iteration += 1;
    }
}

/*! \brief
 * Sample observables for future updates or analysis.
 *
 * \param[in,out] bias            The AWH bias.
 * \param[in] probWeightNeighbor  Probability weights of the neighbors.
 * \param[in]     convolvedBias   The convolved bias.
 */
static void doObservableSampling(Bias                      *bias,
                                 const std::vector<double> &probWeightNeighbor,
                                 double                     convolvedBias)
{
    /* Sampling-based deconvolution extracting the PMF */
    sample_pmf(bias, convolvedBias);

    /* Save probability weights for the update */
    sampleProbabilityWeights(bias, probWeightNeighbor);
}

/*! \brief
 * Do what the AWH bias needs to do every step.
 *
 * The main objectives of the AWH bias step is to 1) set the bias force
 * and potential; 2) update the free energy and bias; and 3) sample
 * observables (e.g. the PMF).
 *
 * \param[in,out] bias           The AWH bias.
 * \param[in] convolveForce      True if the bias force is convolved, otherwise apply single umbrella force.
 * \param[in] probWeightNeighbor Working arrary to store weights of neighbors.
 * \param[out]    awhForce,      The AWH bias force.
 * \param[in,out] awhPotential   Bias potential from all AWH biases.
 * \param[in,out] potential_jump Change in bias potential for all AWH biases.
 * \param[in] ms                 Struct for multi-simulation communication.
 * \param[in] t                  Time.
 * \param[in] step               Time step.
 * \param[in] seed               Random seed.
 * \param[in,out] fplog          Log file.
 */
static void doBiasStep(Bias *bias,
                       bool convolveForce,
                       std::vector<double> *probWeightNeighbor,
                       awh_dvec biasForce,
                       double *awhPotential, double *potential_jump,
                       const gmx_multisim_t *ms,
                       double t, gmx_int64_t step, int seed,
                       FILE *fplog)
{
    BiasState *state = &bias->state;

    /* The grid point closest to the coordinate value defines the current neighborhood of points. Besides
       at steps when global updates and/or checks are performed, only the neighborhood will be touched. */
    state->coordValueIndex = get_closest_index_in_grid(state->coordValue, bias->grid.get());

    const BiasParams &params = bias->params;

    /* If the convolved force is needed or this is a sampling step, the bias in the current neighborhood
       needs to be up-to-date and the probablity weights need to be calculated. */
    double convolvedBias;
    if (convolveForce || do_at_step(bias->params.nstsample_coord, step))
    {
        if (params.canSkipUpdates())
        {
            doSkippedUpdatesInNeighborhood(bias, step);
        }

        updateProbabilityWeightsAndConvolvedBias(bias, probWeightNeighbor, &convolvedBias);

        if (step > 0 && do_at_step(params.nstsample_coord, step))
        {
            doObservableSampling(bias, *probWeightNeighbor, convolvedBias);
        }
    }

    /* Set the bias force and get the potential contribution from this bias.
       The potential jump occurs at different times depending on how the force is applied
       (and how the potential is normalized). For the convolved force it happens when
       the bias is updated, for the umbrella when the umbrella is moved. */
    double potential, newPotential;
    if (convolveForce)
    {
        calcConvolvedForce(*bias, *probWeightNeighbor, biasForce);
        potential    =  -convolvedBias*params.invBeta;
        newPotential = potential; /* Assume no jump */
    }
    else
    {
        /* Umbrella force */
        GMX_RELEASE_ASSERT(bias->coordPoint[state->refCoordpoint].inTargetRegion(),
                           "AWH bias coordinate reference value is outside of the target region.");
        potential =
            calcUmbrellaForceAndPotential(*bias, state->refCoordpoint, biasForce);

        /* Moving the umbrella results in a force correction and a new potential. The umbrella center is sampled as
           often as the coordinate so we know the probability weights needed for moving the umbrella are up-to-date. */
        if (do_at_step(params.nstsample_coord, step))
        {
            newPotential = moveUmbrella(bias, *probWeightNeighbor, biasForce, step, seed);
        }
        else
        {
            newPotential = potential;
        }
    }

    /* Add this bias potential to the total sum. */
    *awhPotential += potential;

    /* Update the free energy estimates and bias and other history dependent method parameters */
    if (step > 0 && do_at_step(params.nstupdate_free_energy, step))
    {
        update(bias, ms, t, step, fplog);

        if (convolveForce)
        {
            /* The update results in a potential jump, so we need the new convolved potential. */
            newPotential = -calcConvolvedBias(*bias, state->coordValue)*params.invBeta;
        }
    }

    /* Add the potential jump of this bias to the total sum. */
    *potential_jump += (newPotential - potential);

    /* Check the sampled data (histograms) and potentially warn user if something is suspicious */
    checkOnData(*bias, t, step, fplog);
}

/*! \brief
 * Update the coordinate value with the corresponding pull coordinate value.
 *
 * \param[in,out] bias        The AWH bias
 * \param[in, out] pull_work  Pull struct.
 * \param[in] pbc             Struct with information about periodic boundary conditions.
 */
static void setBiasCoordValue(Bias *bias, struct pull_t *pull_work, const t_pbc *pbc)
{
    /* Keep own copy of current coordinate value. */
    for (int d = 0; d < bias->ndim; d++)
    {
        get_pull_coord_value(pull_work, bias->dimParams[d].pullCoordIndex, pbc, &bias->state.coordValue[d]);
    }
}

/* Do an AWH biasing update. */
real update_awh(AwhBiasCollection      *awh,
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
    wallcycle_start(wallcycle, ewcAWH);

    t_pbc  pbc;
    set_pbc(&pbc, ePBC, box);

    /* During the AWH step the potential can instantaneously jump due to either
       an bias update or moving the umbrella. The jumps are kept track of and
       subtracted from the potential in order to get a useful conserved energy quantity. */
    double potentialJump = 0;
    double biasPotential = awh->potential_offset;

    for (auto &bias : awh->bias)
    {
        /* Update the AWH coordinate values with those of the corresponding
         * pull coordinates.
         */
        setBiasCoordValue(&bias, pull_work, &pbc);

        /* Perform an AWH biasing step: this means, at regular intervals,
         * sampling observables based on the input pull coordinate value,
         * setting the bias force and/or updating the AWH bias state.
         */
        awh_dvec biasForce;
        doBiasStep(&bias,
                   awh->convolveForce,
                   &awh->work,
                   biasForce, &biasPotential, &potentialJump,
                   ms, t, step, awh->seed, fplog);

        /* Keep track of the total potential shift needed to remove the potential jumps. */
        awh->potential_offset -= potentialJump;

        /* Communicate the bias force to the pull struct.
         * The bias potential is returned at the end of this function,
         * so that it can be added externally to the correct energy data block.
         */
        applyBiasForceToPullCoords(bias, biasForce,
                                   pull_work, mdatoms, force, virial);
    }

    wallcycle_stop(wallcycle, ewcAWH);

    return static_cast<real>(biasPotential);
}
