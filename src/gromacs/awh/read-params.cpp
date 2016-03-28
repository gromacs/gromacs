/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
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

#include "read-params.h"

#include "gromacs/awh/awh.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/random/seed.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "internal.h"
#include "math.h"

const char *eawhtarget_names[eawhtargetNR+1] = {
    "constant", "cutoff", "boltzmann", "weighthistogram", nullptr
};

const char *eawhgrowth_names[eawhgrowthNR+1] = {
    "exp-linear", "none", "linear", nullptr
};

const char *eawhpotential_names[eawhpotentialNR+1] = {
    "convolved", "umbrella", nullptr
};

static void read_dim_params(int *ninp_p, t_inpfile **inp_p, const char *prefix,
                            awh_dim_params_t *dim_params, const pull_params_t *pull_params,
                            warninp_t wi, bool bComment)
{
    int        icoord_ext, icoord, eGeom;
    char       warningmsg[STRLEN], opt[STRLEN];
    int        ninp;
    t_inpfile *inp;

    ninp   = *ninp_p;
    inp    = *inp_p;

    if (bComment)
    {
        CTYPE("The pull coordinate for each AWH dimension");
    }
    sprintf(opt, "%s-pull-coord", prefix);
    ITYPE(opt, icoord_ext, 0);
    if (icoord_ext <  1)
    {
        gmx_fatal(FARGS, "Failed to read a valid pull coordinate index for %s. "
                  "Note that the pull coordinate indexing starts at 1.", opt);
    }

    /* The pull coordinate indices start at 1 in the input file, at 0 internally */
    icoord                       = icoord_ext - 1;
    dim_params->pull_coord_index = icoord;

    /* The pull settings need to be consistent with the AWH settings */
    if (!(pull_params->coord[icoord].eType == epullEXTERNAL) )
    {
        gmx_fatal(FARGS, "AWH biasing can only be  applied to pull type %s",
                  EPULLTYPE(epullEXTERNAL));
    }

    if (icoord >= pull_params->ncoord)
    {
        gmx_fatal(FARGS, "The given AWH coordinate index (%d) is larger than the number of pull coordinates (%d)",
                  icoord_ext, pull_params->ncoord);
    }
    if (pull_params->coord[icoord].rate != 0)
    {
        gmx_fatal(FARGS, "Setting pull-coord%d-rate (%g) is incompatible with AWH biasing this coordinate",
                  icoord_ext, pull_params->coord[icoord].rate);
    }

    if (bComment)
    {
        CTYPE("Estimated diffusion constant (nm^2/ps or rad^2/ps). Hint: a larger value gives a faster initial biasing rate.");
    }
    sprintf(opt, "%s-diffusion", prefix);
    RTYPE(opt, dim_params->diffusion, 0);

    if (dim_params->diffusion <= 0)
    {
        const double diffusion_default = 1e-5;
        sprintf(warningmsg, "%s not explicitly set by user."
                " You can choose to use a default value (%g nm^2/ps or rad^2/ps) but this may very well be non-optimal for your system!",
                opt, diffusion_default);
        warning(wi, warningmsg);
        dim_params->diffusion = diffusion_default;
    }

    /* Grid params for each axis */
    eGeom = pull_params->coord[dim_params->pull_coord_index].eGeom;

    if (bComment)
    {
        CTYPE("Minimum and maximum reference points for each dimension");
    }

    sprintf(opt, "%s-min", prefix);
    RTYPE(opt, dim_params->origin, 0.);

    sprintf(opt, "%s-max", prefix);
    RTYPE(opt, dim_params->end, 0.);

    if (gmx_within_tol(dim_params->end - dim_params->origin, 0, GMX_REAL_EPS))
    {
        sprintf(warningmsg, "The given interval length given by %s-min (%g) and %s-max (%g) is zero. "
                "This will result in only one point along this axis in the AWH grid.",
                prefix, dim_params->origin, prefix, dim_params->end);
        warning(wi, warningmsg);
    }
    /* Check that the requested interval is in allowed range */
    if (eGeom == epullgDIST)
    {
        if (dim_params->origin < 0 || dim_params->end < 0)
        {
            gmx_fatal(FARGS, "%s-min (%g) or %s-max (%g) set to a negative value. With pull geometry distance coordinate values are non-negative. "
                      "Perhaps you want to use geometry %s instead?",
                      prefix, dim_params->origin, prefix, dim_params->end, EPULLGEOM(epullgDIR));
        }
    }
    else if (eGeom == epullgANGLE || eGeom == epullgANGLEAXIS)
    {
        if (dim_params->origin < 0 || dim_params->end > 180)
        {
            gmx_fatal(FARGS, "%s-min (%g) and %s-max (%g) are outside of the allowed range 0 to 180 deg for pull geometries %s and %s ",
                      prefix, dim_params->origin, prefix, dim_params->end, EPULLGEOM(epullgANGLE), EPULLGEOM(epullgANGLEAXIS));
        }
    }
    else if (eGeom == epullgDIHEDRAL)
    {
        if (dim_params->origin < -180 || dim_params->end > 180)
        {
            gmx_fatal(FARGS, "%s-min (%g) and %s-max (%g) are outside of the allowed range -180 to 180 deg for pull geometry %s. ",
                      prefix, dim_params->origin, prefix, dim_params->end, EPULLGEOM(epullgDIHEDRAL));
        }
    }
    if (bComment)
    {
        CTYPE("The number of subintervals to partition the AWH coordinate interval into");
    }
    sprintf(opt, "%s-ninterval", prefix);
    ITYPE(opt, dim_params->ninterval, 1);
    if (dim_params->ninterval > 1)
    {
    }
    if (dim_params->ninterval <= 0)
    {
        gmx_fatal(FARGS, "%s (%d) needs to a positive number", opt, dim_params->ninterval);
    }

    CTYPE("Overlap of the AWH coordinate subintervals when partitioning the interval");
    sprintf(opt, "%s-interval-overlap", prefix);
    RTYPE(opt, dim_params->interval_overlap, 1.);
    if (dim_params->interval_overlap < 0 || dim_params->interval_overlap > 1)
    {
        gmx_fatal(FARGS, "%s (%g) only takes values in the range of 0 to 1",
                  opt, dim_params->interval_overlap);
    }

    *ninp_p   = ninp;
    *inp_p    = inp;
}

static void check_input_consistency_awh_bias(const awh_bias_params_t *awh_bias_params, warninp_t wi)
{
    char warningmsg[STRLEN];

    /* Settings for partitioning the AWH domain */
    for (int d = 0; d < awh_bias_params->ndim; d++)
    {
        if (awh_bias_params->dim_params[d].ninterval > 1)
        {
            if (awh_bias_params->bShare)
            {
                gmx_fatal(FARGS, "Partitioning the AWH dimension %d into %d intervals is not compatible "
                          "with having sharing simulations since all intervals need to be identical.",
                          d + 1, awh_bias_params->dim_params[d].ninterval);
            }
            else if (awh_bias_params->dim_params[d].interval_overlap == 1)
            {
                sprintf(warningmsg, "Splitting the interval of independent simulations "
                        "with an overlap of 1 has no effect. "
                        "Perhaps you forgot to set the interval overlap in the mdp file?" );
                warning(wi, warningmsg);
            }
        }
        else if ((awh_bias_params->dim_params[d].ninterval == 1) && (awh_bias_params->dim_params[d].interval_overlap < 1))
        {
            gmx_fatal(FARGS, "With a single AWH subinterval the AWH overlap can only equal 1. "
                      "Perhaps you forgot to set the number of AWH subintervals in the mdp file?");
        }
    }
}

static void read_bias_params(int *ninp_p, t_inpfile **inp_p, awh_bias_params_t *awh_bias_params, const char *prefix,
                             const t_inputrec *ir, warninp_t wi, bool bComment)
{
    int        ninp;
    t_inpfile *inp;
    char       opt[STRLEN], prefixdim[STRLEN];
    char       warningmsg[STRLEN];

    /* These are assumed to be declared by the gromacs reading functions */
    ninp   = *ninp_p;
    inp    = *inp_p;

    if (bComment)
    {
        CTYPE("Target distribution type");
    }
    sprintf(opt, "%s-target", prefix);
    EETYPE(opt, awh_bias_params->eTarget, eawhtarget_names);

    if (bComment)
    {
        CTYPE("Parameter value used for certain target distributions");
    }
    sprintf(opt, "%s-target-paramvalue", prefix);
    RTYPE(opt, awh_bias_params->target_param, 0);

    switch (awh_bias_params->eTarget)
    {
        case eawhtargetCONSTANT:
            if (awh_bias_params->target_param != 0)
            {
                sprintf(warningmsg, "%s was set explicitly but will not be used for target type %s.",
                        opt, EAWHTARGET(awh_bias_params->eTarget));
                warning(wi, warningmsg);
            }
            break;
        case eawhtargetCUTOFF:
            if (awh_bias_params->target_param <= 0)
            {
                gmx_fatal(FARGS, "%s = %g is not useful for target type %s.",
                          opt, awh_bias_params->target_param, EAWHTARGET(awh_bias_params->eTarget));
            }
            break;
        case eawhtargetBOLTZMANN:
        case eawhtargetWEIGHTHIST:
            if (awh_bias_params->target_param < 0 || awh_bias_params->target_param > 1)
            {
                gmx_fatal(FARGS, "%s = %g is not useful for target type %s.",
                          opt, awh_bias_params->target_param, EAWHTARGET(awh_bias_params->eTarget));
            }
            break;
    }

    if (bComment)
    {
        CTYPE("Growth rate of the AWH biasing histogram");
    }
    sprintf(opt, "%s-growth", prefix);
    EETYPE(opt, awh_bias_params->eGrowth, eawhgrowth_names);

    if (bComment)
    {
        CTYPE("Initialize PMF and target with user data: yes or no");
    }
    sprintf(opt, "%s-user-data", prefix);
    EETYPE(opt, awh_bias_params->bUser_data, yesno_names);

    if (bComment)
    {
        CTYPE("Estimated initial free energy error (kT). Hint: a larger value gives a faster initial biasing rate.");
    }
    sprintf(opt, "%s-error-init", prefix);

    /* We allow using a default value here (and just warn the user if the diffusion constant is not set) */
    RTYPE(opt, awh_bias_params->error_initial, 4);

    if (bComment)
    {
        CTYPE("Dimensionality of the coordinate");
    }
    sprintf(opt, "%s-ndim", prefix);
    ITYPE(opt, awh_bias_params->ndim, 0);

    if (awh_bias_params->ndim <= 0)
    {
        gmx_fatal(FARGS, "%s-ndim (%d) should be > 0\n", prefix,  awh_bias_params->ndim);
    }
    if (awh_bias_params->ndim > 2)
    {
        warning_note(wi, "For awh-dim > 2 the estimate based on the diffusion and the initial error is currently only a rough guideline."
                     " You should verify its usefulness for your system before production runs!");
    }

    CTYPE("Share the AWH bias across multiple simulations: yes or no");
    sprintf(opt, "%s-share", prefix);
    EETYPE(opt, awh_bias_params->bShare, yesno_names);

    snew(awh_bias_params->dim_params, awh_bias_params->ndim);
    for (int d = 0; d < awh_bias_params->ndim; d++)
    {
        bComment = bComment && d == 0;
        sprintf(prefixdim, "%s-dim%d", prefix, d + 1);
        read_dim_params(&ninp, &inp, prefixdim, &awh_bias_params->dim_params[d], ir->pull, wi, bComment);
    }

    /* Check consistencies here that cannot be checked at read time at a lower level. */
    check_input_consistency_awh_bias(awh_bias_params, wi);

    *ninp_p   = ninp;
    *inp_p    = inp;
}

static void check_input_consistency_awh(const awh_params_t *awh_params)
{
    /* Each pull coord can map to at most 1 AWH coord */
    for (int k1 = 0; k1 < awh_params->nbias; k1++)
    {
        awh_bias_params_t *awh_bias_params1 = &awh_params->awh_bias_params[k1];

        /* k1 is the reference AWH, k2 is the AWH we compare with (can be equal to k1) */
        for (int k2 = k1; k2 < awh_params->nbias; k2++)
        {
            for (int d1 = 0; d1 < awh_bias_params1->ndim; d1++)
            {
                awh_bias_params_t *awh_bias_params2 = &awh_params->awh_bias_params[k2];

                /* d1 is the reference dimension of the reference AWH. d2 is the dim index of the AWH to compare with. */
                for (int d2 = 0; d2 < awh_bias_params2->ndim; d2++)
                {
                    /* Give an error if (d1, k1) is different from (d2, k2) but the pull coordinate is the same */
                    if ( (d1 != d2 || k1 != k2) && (awh_bias_params1->dim_params[d1].pull_coord_index == awh_bias_params2->dim_params[d2].pull_coord_index) )
                    {
                        char errormsg[STRLEN];
                        sprintf(errormsg, "One pull coordinate (%d) cannot be mapped to two separate AWH dimensions (awh%d-dim%d and awh%d-dim%d). "
                                "If this is really what you want to do you will have to duplicate this pull coordinate.",
                                awh_bias_params1->dim_params[d1].pull_coord_index + 1, k1 + 1, d1 + 1, k2 + 1, d2 + 1);
                        gmx_fatal(FARGS, errormsg);
                    }
                }
            }
        }
    }
}

awh_params_t *read_awh_params(int *ninp_p, t_inpfile **inp_p, const t_inputrec *ir, warninp_t wi)
{
    int               ninp;
    t_inpfile        *inp;
    char              opt[STRLEN], prefix[STRLEN], prefixawh[STRLEN];

    awh_params_t     *awh_params;
    snew(awh_params, 1);

    ninp   = *ninp_p;
    inp    = *inp_p;

    sprintf(prefix, "%s", "awh");

    CTYPE("The number of independent AWH biases");
    sprintf(opt, "%s-nbias", prefix);
    ITYPE(opt, awh_params->nbias, 1);
    if (awh_params->nbias <= 0)
    {
        gmx_fatal(FARGS, "%s needs to be an integer > 0", opt);
    }

    CTYPE("The AWH random seed");
    sprintf(opt, "%s-seed", prefix);
    ITYPE(opt, awh_params->seed, -1);
    if (awh_params->seed == -1)
    {
        /* TODO: I'm not sure why for the ld and expanded (in grompp.cpp) the check was changed from -1 to 0
           although those default values are still -1 ? */
        awh_params->seed = static_cast<int>(gmx::makeRandomSeed());
        fprintf(stderr, "Setting the AWH bias MC random seed to %" GMX_PRId64 "\n", awh_params->seed);
    }

    CTYPE("Number of steps per AWH printing");
    sprintf(opt, "%s-nstout", prefix);
    ITYPE(opt, awh_params->nstout, ir->nstenergy);
    if (awh_params->nstout > 0 && (ir->nstenergy <= 0 || (awh_params->nstout % ir->nstenergy != 0)))
    {
        gmx_fatal(FARGS, "%s (%d) needs to be a multiple of nstenergy (%d)",
                  opt, awh_params->nstout, ir->nstenergy);
    }

    snew(awh_params->awh_bias_params, awh_params->nbias);

    CTYPE("Number of steps per coordinate sample");
    sprintf(opt, "%s-nstsample", prefix);
    ITYPE(opt, awh_params->nstsample_coord, 10);

    CTYPE("Number of samples per moving the coordinate reference value");
    sprintf(opt, "%s-nsamples-refvalue", prefix);
    ITYPE(opt, awh_params->nsamples_move_refvalue, 5);

    CTYPE("Number of coordinate samples per free energy update and bias update");
    sprintf(opt, "%s-nsamples-update", prefix);
    ITYPE(opt, awh_params->nsamples_update_free_energy, 5);

    CTYPE("The way to apply the biasing potential: convolved or umbrella");
    sprintf(opt, "%s-potential", prefix);
    EETYPE(opt, awh_params->ePotential, eawhpotential_names);

    /* Read the parameters specific to each AWH bias */
    for (int k = 0; k < awh_params->nbias; k++)
    {
        bool bComment = (k == 0);
        sprintf(prefixawh, "%s%d", prefix, k + 1);
        read_bias_params(&ninp, &inp, &awh_params->awh_bias_params[k], prefixawh, ir, wi, bComment);
    }

    /* Do a final consistency check before returning */
    check_input_consistency_awh(awh_params);

    *ninp_p   = ninp;
    *inp_p    = inp;

    return awh_params;
}

static double get_pull_coord_period(const pull_params_t *pull_params,
                                    int                  coord_ind,
                                    const matrix         box)
{
    double        period;
    t_pull_coord *pcrd_params = &pull_params->coord[coord_ind];

    if (pcrd_params->eGeom == epullgDIRPBC)
    {
        /* For direction periodic, we need the pull vector to be one of the box vectors
           (or more generally I guess it could be a integer combination of boxvectors).
           This boxvector should to be orthogonal to the (periodic) plane spanned by the other two box vectors.
           Here we assume that the pull vector is either x, y or z.
         * E.g. for pull vec = (1, 0, 0) the box vector tensor should look like:
         * | x 0 0 |
         * | 0 a c |
         * | 0 b d |
         *
           The period is then given by the box length x.

           Note: we make these checks here for AWH and not in pull because we allow pull to be more general.
         */
        int m_pullvec, count_nonzeros = 0;

        /* Check that pull vec has only one component and which component it is. This component gives the relevant box vector */
        for (int m = 0; m < DIM; m++)
        {
            if (pcrd_params->vec[m] != 0)
            {
                m_pullvec = m;
                count_nonzeros++;
            }
        }
        if (count_nonzeros != 1)
        {
            gmx_fatal(FARGS, "For AWH biasing pull coordinate %d with pull geometry %s, the pull vector needs to be parallel to "
                      "a box vector that is parallel to either the x, y or z axis and is orthogonal to the other box vectors.",
                      coord_ind + 1, EPULLGEOM(epullgDIRPBC));
        }

        /* Check that there is a box vec parallel to pull vec and that this boxvec is orthogonal to the other box vectors */
        for (int m = 0; m < DIM; m++)
        {
            for (int n = 0; n < DIM; n++)
            {
                if ((n != m) && (n == m_pullvec || m == m_pullvec) && box[m][n] > 0)
                {
                    gmx_fatal(FARGS, "For AWH biasing pull coordinate %d with pull geometry %s, there needs to be a box vector parallel to the pull vector that is "
                              "orthogonal to the other box vectors.",
                              coord_ind + 1, EPULLGEOM(epullgDIRPBC));
                }
            }
        }

        /* If this box vector only has one component as we assumed the norm should be equal to the absolute value of that component */
        period = static_cast<double>(norm(box[m_pullvec]));
    }
    else if (pcrd_params->eGeom == epullgDIHEDRAL)
    {
        /* The dihedral angle is periodic in -180 to 180 deg */
        period = 360;
    }
    else
    {
        period = 0;
    }

    return period;
}

static bool interval_is_in_periodic_interval(double origin, double end, double period)
{
    return (period == 0) || (std::fabs(origin) <= 0.5*period && std::fabs(end) <= 0.5*period);
}

static bool value_is_in_interval(double origin, double end, double period, double value)
{
    bool bIn_interval;

    if (period > 0)
    {
        if (origin < end)
        {
            /* The interval closes within the periodic interval */
            bIn_interval = (value >= origin) && (value <= end);
        }
        else
        {
            /* The interval wraps around the periodic boundary */
            bIn_interval = ((value >= origin) && (value <= 0.5*period)) || ((value >= -0.5*period) && (value <= end));
        }
    }
    else
    {
        bIn_interval = (value >= origin) && (value <= end);
    }

    return bIn_interval;
}

/* Checks that the pull coordinate input is consistent with the requested AWH interval */
static void check_input_consistency_interval(const awh_params_t *awh_params, warninp_t wi)
{
    for (int k = 0; k < awh_params->nbias; k++)
    {
        awh_bias_params_t *awh_bias_params = &awh_params->awh_bias_params[k];
        for (int d = 0; d < awh_bias_params->ndim; d++)
        {
            awh_dim_params_t *dim_params            = &awh_bias_params->dim_params[d];
            int               pull_coord_index      = dim_params->pull_coord_index;
            double            origin                = dim_params->origin, end = dim_params->end, period = dim_params->period;
            double            coord_value_init      = dim_params->coord_value_init;

            if ((period == 0) && (origin > end))
            {
                gmx_fatal(FARGS, "For the non-periodic pull coordinates awh%d-dim%d-min cannot be larger than awh%d-dim%d-max",
                          k + 1, d + 1, origin, k + 1, d + 1, end);
            }

            /* Currently we assume symmetric periodic intervals, meaning we use [-period/2, period/2] as the reference interval.
               Make sure the AWH interval is within this reference interval.

               Note: we could fairly simply allow using a  more general interval (e.g. [x, x + period]) but it complicates
               things slightly and I don't see that there is a great need for it. It would also mean that the interval would
               depend on AWH input. Also, for dihedral angles you would always want the reference interval to be -180, +180,
               independent of AWH parameters.
             */
            if (!interval_is_in_periodic_interval(origin, end, period))
            {
                gmx_fatal(FARGS, "When using AWH with periodic pull coordinate geometries awh%d-dim%d-min (%.8g) and "
                          "awh%d-dim%d-max (%.8g) should cover at most one period (%.8g) and take values in between "
                          "minus half a period and plus half a period, i.e. in the interval [%.8g, %.8g].",
                          k + 1, d + 1, origin, k + 1, d + 1, end,
                          period, -0.5*period, 0.5*period);

            }

            /* Warn if the pull initial coordinate value is not in the grid */
            if (!value_is_in_interval(origin, end, period, coord_value_init))
            {
                char       warningmsg[STRLEN];
                sprintf(warningmsg, "The initial pull coordinate value (%.8g) for pull coordinate index (%d) falls outside "
                        "of the AWH interval given by awh%d-dim%d-min (%.8g) to awh%d-dim%d-max (%.8g). "
                        "The coordinate reference value will be set the the closest grid point inside of the interval.",
                        coord_value_init, pull_coord_index + 1,
                        k + 1, d + 1, origin, k + 1, d + 1, end);
                warning(wi, warningmsg);
            }
        }
    }
}

/* Some AWH parameters are dependent on state parameters, e.g. the initial configuration which we don't have at the time we read the
   mdp file. We could read this at mdrun time, but there might be inconsistencies giving errors and we would rather throw those errors
   before earlier. */
void set_state_dependent_awh_params(awh_params_t *awh_params,
                                    const pull_params_t *pull_params, pull_t *pull_work,
                                    const t_inputrec *inputrec, const matrix box,  int ePBC,
                                    const t_grpopts *inputrec_group_options, warninp_t wi)
{
    /* The temperature is not really state depenendent but is not know when read_awh_params is called (in get ir).
       It is known first after do_index has been called in grompp.cpp. */
    if (inputrec_group_options->ref_t == NULL || inputrec_group_options->ref_t[0] <= 0)
    {
        gmx_fatal(FARGS, "AWH biasing is only supported for temperatures > 0");
    }
    for (int i = 1; i < inputrec_group_options->ngtc; i++)
    {
        if (inputrec_group_options->ref_t[i] != inputrec_group_options->ref_t[0])
        {
            gmx_fatal(FARGS, "AWH biasing is currently only supported for identical temperatures for all temperature coupling groups");
        }
    }

    t_pbc          pbc;
    set_pbc(&pbc, ePBC, box);

    for (int k = 0; k < awh_params->nbias; k++)
    {
        awh_bias_params_t *awh_bias_params = &awh_params->awh_bias_params[k];
        for (int d = 0; d < awh_bias_params->ndim; d++)
        {
            awh_dim_params_t *dim_params = &awh_bias_params->dim_params[d];

            /* The periodiciy of the AWH grid in certain cases depends on the simulation box */
            dim_params->period = get_pull_coord_period(pull_params, dim_params->pull_coord_index, box);

            /* The initial coordinate value */
            get_pull_coord_value(pull_work, dim_params->pull_coord_index, &pbc, &dim_params->coord_value_init);
        }
    }
    check_input_consistency_interval(awh_params, wi);

    /* Make a temporary working awh (with NULL log file and commrec) and register
       AWH as external potential with pull. */
    awh_t     *awh = init_awh(NULL, inputrec, NULL, awh_params);
    register_bias_with_pull(awh, pull_work);
}
