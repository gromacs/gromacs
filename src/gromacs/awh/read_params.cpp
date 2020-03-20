/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017 by the GROMACS development team.
 * Copyright (c) 2018,2019,2020, by the GROMACS development team, led by
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

#include "read_params.h"

#include "gromacs/awh/awh.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/fileio/warninp.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/awh_params.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull_params.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pulling/pull.h"
#include "gromacs/random/seed.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "biasparams.h"
#include "biassharing.h"

namespace gmx
{

const char* eawhtarget_names[eawhtargetNR + 1] = { "constant", "cutoff", "boltzmann",
                                                   "local-boltzmann", nullptr };

const char* eawhgrowth_names[eawhgrowthNR + 1] = { "exp-linear", "linear", nullptr };

const char* eawhpotential_names[eawhpotentialNR + 1] = { "convolved", "umbrella", nullptr };

const char* eawhcoordprovider_names[eawhcoordproviderNR + 1] = { "pull", nullptr };

/*! \brief
 * Read parameters of an AWH bias dimension.
 *
 * \param[in,out] inp        Input file entries.
 * \param[in] prefix         Prefix for dimension parameters.
 * \param[in,out] dimParams  AWH dimensional parameters.
 * \param[in,out] wi         Struct for bookeeping warnings.
 * \param[in] bComment       True if comments should be printed.
 */
static void readDimParams(std::vector<t_inpfile>* inp,
                          const std::string&      prefix,
                          AwhDimParams*           dimParams,
                          warninp_t               wi,
                          bool                    bComment)
{
    std::string opt;
    if (bComment)
    {
        printStringNoNewline(
                inp, "The provider of the reaction coordinate, currently only pull is supported");
    }

    opt                       = prefix + "-coord-provider";
    dimParams->eCoordProvider = get_eeenum(inp, opt, eawhcoordprovider_names, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "The coordinate index for this dimension");
    }
    opt = prefix + "-coord-index";
    int coordIndexInput;
    coordIndexInput = get_eint(inp, opt, 1, wi);
    if (coordIndexInput < 1)
    {
        gmx_fatal(FARGS,
                  "Failed to read a valid coordinate index for %s. "
                  "Note that the pull coordinate indexing starts at 1.",
                  opt.c_str());
    }

    /* The pull coordinate indices start at 1 in the input file, at 0 internally */
    dimParams->coordIndex = coordIndexInput - 1;

    if (bComment)
    {
        printStringNoNewline(inp, "Start and end values for each coordinate dimension");
    }

    opt               = prefix + "-start";
    dimParams->origin = get_ereal(inp, opt, 0., wi);

    opt            = prefix + "-end";
    dimParams->end = get_ereal(inp, opt, 0., wi);

    if (bComment)
    {
        printStringNoNewline(
                inp, "The force constant for this coordinate (kJ/mol/nm^2 or kJ/mol/rad^2)");
    }
    opt                      = prefix + "-force-constant";
    dimParams->forceConstant = get_ereal(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Estimated diffusion constant (nm^2/ps or rad^2/ps)");
    }
    opt                  = prefix + "-diffusion";
    dimParams->diffusion = get_ereal(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Diameter that needs to be sampled around a point before it is "
                             "considered covered.");
    }
    opt                      = prefix + "-cover-diameter";
    dimParams->coverDiameter = get_ereal(inp, opt, 0, wi);
}

/*! \brief
 * Check the parameters of an AWH bias dimension.
 *
 * \param[in] prefix         Prefix for dimension parameters.
 * \param[in,out] dimParams  AWH dimensional parameters.
 * \param[in] pull_params    Pull parameters.
 * \param[in,out] wi         Struct for bookeeping warnings.
 */
static void checkDimParams(const std::string&   prefix,
                           AwhDimParams*        dimParams,
                           const pull_params_t* pull_params,
                           warninp_t            wi)
{
    std::string opt;

    /* The pull settings need to be consistent with the AWH settings */
    if (!(pull_params->coord[dimParams->coordIndex].eType == epullEXTERNAL))
    {
        gmx_fatal(FARGS, "AWH biasing can only be  applied to pull type %s", EPULLTYPE(epullEXTERNAL));
    }

    if (dimParams->coordIndex >= pull_params->ncoord)
    {
        gmx_fatal(FARGS,
                  "The given AWH coordinate index (%d) is larger than the number of pull "
                  "coordinates (%d)",
                  dimParams->coordIndex + 1, pull_params->ncoord);
    }
    if (pull_params->coord[dimParams->coordIndex].rate != 0)
    {
        auto message = formatString(
                "Setting pull-coord%d-rate (%g) is incompatible with AWH biasing this coordinate",
                dimParams->coordIndex + 1, pull_params->coord[dimParams->coordIndex].rate);
        warning_error(wi, message);
    }

    /* BiasGrid params for each axis */
    int eGeom = pull_params->coord[dimParams->coordIndex].eGeom;

    if (gmx_within_tol(dimParams->end - dimParams->origin, 0, GMX_REAL_EPS))
    {
        auto message = formatString(
                "The given interval length given by %s-start (%g) and %s-end (%g) is zero. "
                "This will result in only one point along this axis in the coordinate value grid.",
                prefix.c_str(), dimParams->origin, prefix.c_str(), dimParams->end);
        warning(wi, message);
    }
    /* Check that the requested interval is in allowed range */
    if (eGeom == epullgDIST)
    {
        if (dimParams->origin < 0 || dimParams->end < 0)
        {
            gmx_fatal(FARGS,
                      "%s-start (%g) or %s-end (%g) set to a negative value. With pull geometry "
                      "distance coordinate values are non-negative. "
                      "Perhaps you want to use geometry %s instead?",
                      prefix.c_str(), dimParams->origin, prefix.c_str(), dimParams->end,
                      EPULLGEOM(epullgDIR));
        }
    }
    else if (eGeom == epullgANGLE || eGeom == epullgANGLEAXIS)
    {
        if (dimParams->origin < 0 || dimParams->end > 180)
        {
            gmx_fatal(FARGS,
                      "%s-start (%g) and %s-end (%g) are outside of the allowed range 0 to 180 deg "
                      "for pull geometries %s and %s ",
                      prefix.c_str(), dimParams->origin, prefix.c_str(), dimParams->end,
                      EPULLGEOM(epullgANGLE), EPULLGEOM(epullgANGLEAXIS));
        }
    }
    else if (eGeom == epullgDIHEDRAL)
    {
        if (dimParams->origin < -180 || dimParams->end > 180)
        {
            gmx_fatal(FARGS,
                      "%s-start (%g) and %s-end (%g) are outside of the allowed range -180 to 180 "
                      "deg for pull geometry %s. ",
                      prefix.c_str(), dimParams->origin, prefix.c_str(), dimParams->end,
                      EPULLGEOM(epullgDIHEDRAL));
        }
    }

    opt = prefix + "-force-constant";
    if (dimParams->forceConstant <= 0)
    {
        warning_error(wi, "The force AWH bias force constant should be > 0");
    }

    opt = prefix + "-diffusion";
    if (dimParams->diffusion <= 0)
    {
        const double diffusion_default = 1e-5;
        auto         message           = formatString(
                "%s not explicitly set by user. You can choose to use a default "
                "value (%g nm^2/ps or rad^2/ps) but this may very well be "
                "non-optimal for your system!",
                opt.c_str(), diffusion_default);
        warning(wi, message);
        dimParams->diffusion = diffusion_default;
    }

    opt = prefix + "-cover-diameter";
    if (dimParams->coverDiameter < 0)
    {
        gmx_fatal(FARGS, "%s (%g) cannot be negative.", opt.c_str(), dimParams->coverDiameter);
    }
}

/*! \brief
 * Check consistency of input at the AWH bias level.
 *
 * \param[in]     awhBiasParams  AWH bias parameters.
 * \param[in,out] wi             Struct for bookkeeping warnings.
 */
static void checkInputConsistencyAwhBias(const AwhBiasParams& awhBiasParams, warninp_t wi)
{
    /* Covering diameter and sharing warning. */
    for (int d = 0; d < awhBiasParams.ndim; d++)
    {
        double coverDiameter = awhBiasParams.dimParams[d].coverDiameter;
        if (awhBiasParams.shareGroup <= 0 && coverDiameter > 0)
        {
            warning(wi,
                    "The covering diameter is only relevant to set for bias sharing simulations.");
        }
    }
}

/*! \brief
 * Read parameters of an AWH bias.
 *
 * \param[in,out] inp            Input file entries.
 * \param[in,out] awhBiasParams  AWH dimensional parameters.
 * \param[in]     prefix         Prefix for bias parameters.
 * \param[in,out] wi             Struct for bookeeping warnings.
 * \param[in]     bComment       True if comments should be printed.
 */
static void readBiasParams(std::vector<t_inpfile>* inp,
                           AwhBiasParams*          awhBiasParams,
                           const std::string&      prefix,
                           warninp_t               wi,
                           bool                    bComment)
{
    if (bComment)
    {
        printStringNoNewline(inp, "Estimated initial PMF error (kJ/mol)");
    }

    std::string opt             = prefix + "-error-init";
    awhBiasParams->errorInitial = get_ereal(inp, opt, 10, wi);

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Growth rate of the reference histogram determining the bias update "
                             "size: exp-linear or linear");
    }
    opt                    = prefix + "-growth";
    awhBiasParams->eGrowth = get_eeenum(inp, opt, eawhgrowth_names, wi);

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Start the simulation by equilibrating histogram towards the target "
                             "distribution: no or yes");
    }
    opt                                 = prefix + "-equilibrate-histogram";
    awhBiasParams->equilibrateHistogram = (get_eeenum(inp, opt, yesno_names, wi) != 0);

    if (bComment)
    {
        printStringNoNewline(
                inp, "Target distribution type: constant, cutoff, boltzmann or local-boltzmann");
    }
    opt                    = prefix + "-target";
    awhBiasParams->eTarget = get_eeenum(inp, opt, eawhtarget_names, wi);

    if (bComment)
    {
        printStringNoNewline(inp,
                             "Boltzmann beta scaling factor for target distribution types "
                             "'boltzmann' and 'boltzmann-local'");
    }
    opt                              = prefix + "-target-beta-scaling";
    awhBiasParams->targetBetaScaling = get_ereal(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Free energy cutoff value for target distribution type 'cutoff'");
    }
    opt                         = prefix + "-target-cutoff";
    awhBiasParams->targetCutoff = get_ereal(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Initialize PMF and target with user data: no or yes");
    }
    opt                      = prefix + "-user-data";
    awhBiasParams->bUserData = get_eeenum(inp, opt, yesno_names, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Group index to share the bias with, 0 means not shared");
    }
    opt                       = prefix + "-share-group";
    awhBiasParams->shareGroup = get_eint(inp, opt, 0, wi);

    if (bComment)
    {
        printStringNoNewline(inp, "Dimensionality of the coordinate");
    }
    opt                 = prefix + "-ndim";
    awhBiasParams->ndim = get_eint(inp, opt, 0, wi);

    /* Check this before starting to read the AWH dimension parameters. */
    if (awhBiasParams->ndim <= 0 || awhBiasParams->ndim > c_biasMaxNumDim)
    {
        gmx_fatal(FARGS, "%s (%d) needs to be > 0 and at most %d\n", opt.c_str(),
                  awhBiasParams->ndim, c_biasMaxNumDim);
    }
    snew(awhBiasParams->dimParams, awhBiasParams->ndim);
    for (int d = 0; d < awhBiasParams->ndim; d++)
    {
        bComment              = bComment && d == 0;
        std::string prefixdim = prefix + formatString("-dim%d", d + 1);
        readDimParams(inp, prefixdim, &awhBiasParams->dimParams[d], wi, bComment);
    }
}

/*! \brief
 * Check the parameters of an AWH bias.
 *
 * \param[in]     awhBiasParams  AWH dimensional parameters.
 * \param[in]     prefix         Prefix for bias parameters.
 * \param[in]     ir             Input parameter struct.
 * \param[in,out] wi             Struct for bookeeping warnings.
 */
static void checkBiasParams(const AwhBiasParams* awhBiasParams,
                            const std::string&   prefix,
                            const t_inputrec*    ir,
                            warninp_t            wi)
{
    std::string opt = prefix + "-error-init";
    if (awhBiasParams->errorInitial <= 0)
    {
        gmx_fatal(FARGS, "%s needs to be > 0.", opt.c_str());
    }

    opt = prefix + "-equilibrate-histogram";
    if (awhBiasParams->equilibrateHistogram && awhBiasParams->eGrowth != eawhgrowthEXP_LINEAR)
    {
        auto message =
                formatString("Option %s will only have an effect for histogram growth type '%s'.",
                             opt.c_str(), EAWHGROWTH(eawhgrowthEXP_LINEAR));
        warning(wi, message);
    }

    if ((awhBiasParams->eTarget == eawhtargetLOCALBOLTZMANN)
        && (awhBiasParams->eGrowth == eawhgrowthEXP_LINEAR))
    {
        auto message = formatString(
                "Target type '%s' combined with histogram growth type '%s' is not "
                "expected to give stable bias updates. You probably want to use growth type "
                "'%s' instead.",
                EAWHTARGET(eawhtargetLOCALBOLTZMANN), EAWHGROWTH(eawhgrowthEXP_LINEAR),
                EAWHGROWTH(eawhgrowthLINEAR));
        warning(wi, message);
    }

    opt = prefix + "-target-beta-scaling";
    switch (awhBiasParams->eTarget)
    {
        case eawhtargetBOLTZMANN:
        case eawhtargetLOCALBOLTZMANN:
            if (awhBiasParams->targetBetaScaling < 0 || awhBiasParams->targetBetaScaling > 1)
            {
                gmx_fatal(FARGS, "%s = %g is not useful for target type %s.", opt.c_str(),
                          awhBiasParams->targetBetaScaling, EAWHTARGET(awhBiasParams->eTarget));
            }
            break;
        default:
            if (awhBiasParams->targetBetaScaling != 0)
            {
                gmx_fatal(
                        FARGS,
                        "Value for %s (%g) set explicitly but will not be used for target type %s.",
                        opt.c_str(), awhBiasParams->targetBetaScaling,
                        EAWHTARGET(awhBiasParams->eTarget));
            }
            break;
    }

    opt = prefix + "-target-cutoff";
    switch (awhBiasParams->eTarget)
    {
        case eawhtargetCUTOFF:
            if (awhBiasParams->targetCutoff <= 0)
            {
                gmx_fatal(FARGS, "%s = %g is not useful for target type %s.", opt.c_str(),
                          awhBiasParams->targetCutoff, EAWHTARGET(awhBiasParams->eTarget));
            }
            break;
        default:
            if (awhBiasParams->targetCutoff != 0)
            {
                gmx_fatal(
                        FARGS,
                        "Value for %s (%g) set explicitly but will not be used for target type %s.",
                        opt.c_str(), awhBiasParams->targetCutoff, EAWHTARGET(awhBiasParams->eTarget));
            }
            break;
    }

    opt = prefix + "-share-group";
    if (awhBiasParams->shareGroup < 0)
    {
        warning_error(wi, "AWH bias share-group should be >= 0");
    }

    opt = prefix + "-ndim";
    if (awhBiasParams->ndim <= 0 || awhBiasParams->ndim > c_biasMaxNumDim)
    {
        gmx_fatal(FARGS, "%s (%d) needs to be > 0 and at most %d\n", opt.c_str(),
                  awhBiasParams->ndim, c_biasMaxNumDim);
    }
    if (awhBiasParams->ndim > 2)
    {
        warning_note(wi,
                     "For awh-dim > 2 the estimate based on the diffusion and the initial error is "
                     "currently only a rough guideline."
                     " You should verify its usefulness for your system before production runs!");
    }
    for (int d = 0; d < awhBiasParams->ndim; d++)
    {
        std::string prefixdim = prefix + formatString("-dim%d", d + 1);
        checkDimParams(prefixdim, &awhBiasParams->dimParams[d], ir->pull, wi);
    }

    /* Check consistencies here that cannot be checked at read time at a lower level. */
    checkInputConsistencyAwhBias(*awhBiasParams, wi);
}

/*! \brief
 * Check consistency of input at the AWH level.
 *
 * \param[in]     awhParams  AWH parameters.
 * \param[in,out] wi         Struct for bookkeeping warnings.
 */
static void checkInputConsistencyAwh(const AwhParams& awhParams, warninp_t wi)
{
    /* Each pull coord can map to at most 1 AWH coord.
     * Check that we have a shared bias when requesting multisim sharing.
     */
    bool haveSharedBias = false;
    for (int k1 = 0; k1 < awhParams.numBias; k1++)
    {
        const AwhBiasParams& awhBiasParams1 = awhParams.awhBiasParams[k1];

        if (awhBiasParams1.shareGroup > 0)
        {
            haveSharedBias = true;
        }

        /* k1 is the reference AWH, k2 is the AWH we compare with (can be equal to k1) */
        for (int k2 = k1; k2 < awhParams.numBias; k2++)
        {
            for (int d1 = 0; d1 < awhBiasParams1.ndim; d1++)
            {
                const AwhBiasParams& awhBiasParams2 = awhParams.awhBiasParams[k2];

                /* d1 is the reference dimension of the reference AWH. d2 is the dim index of the AWH to compare with. */
                for (int d2 = 0; d2 < awhBiasParams2.ndim; d2++)
                {
                    /* Give an error if (d1, k1) is different from (d2, k2) but the pull coordinate is the same */
                    if ((d1 != d2 || k1 != k2)
                        && (awhBiasParams1.dimParams[d1].coordIndex == awhBiasParams2.dimParams[d2].coordIndex))
                    {
                        char errormsg[STRLEN];
                        sprintf(errormsg,
                                "One pull coordinate (%d) cannot be mapped to two separate AWH "
                                "dimensions (awh%d-dim%d and awh%d-dim%d). "
                                "If this is really what you want to do you will have to duplicate "
                                "this pull coordinate.",
                                awhBiasParams1.dimParams[d1].coordIndex + 1, k1 + 1, d1 + 1, k2 + 1,
                                d2 + 1);
                        gmx_fatal(FARGS, "%s", errormsg);
                    }
                }
            }
        }
    }

    if (awhParams.shareBiasMultisim && !haveSharedBias)
    {
        warning(wi,
                "Sharing of biases over multiple simulations is requested, but no bias is marked "
                "as shared (share-group > 0)");
    }

    /* mdrun does not support this (yet), but will check again */
    if (haveBiasSharingWithinSimulation(awhParams))
    {
        warning(wi,
                "You have shared biases within a single simulation, but mdrun does not support "
                "this (yet)");
    }
}

AwhParams* readAwhParams(std::vector<t_inpfile>* inp, warninp_t wi)
{
    AwhParams* awhParams;
    snew(awhParams, 1);
    std::string opt;

    /* Parameters common for all biases */

    printStringNoNewline(inp, "The way to apply the biasing potential: convolved or umbrella");
    opt                   = "awh-potential";
    awhParams->ePotential = get_eeenum(inp, opt, eawhpotential_names, wi);

    printStringNoNewline(inp,
                         "The random seed used for sampling the umbrella center in the case of "
                         "umbrella type potential");
    opt             = "awh-seed";
    awhParams->seed = get_eint(inp, opt, -1, wi);
    if (awhParams->seed == -1)
    {
        awhParams->seed = static_cast<int>(gmx::makeRandomSeed());
        fprintf(stderr, "Setting the AWH bias MC random seed to %" PRId64 "\n", awhParams->seed);
    }

    printStringNoNewline(inp, "Data output interval in number of steps");
    opt               = "awh-nstout";
    awhParams->nstOut = get_eint(inp, opt, 100000, wi);

    printStringNoNewline(inp, "Coordinate sampling interval in number of steps");
    opt                       = "awh-nstsample";
    awhParams->nstSampleCoord = get_eint(inp, opt, 10, wi);

    printStringNoNewline(inp, "Free energy and bias update interval in number of samples");
    opt                                   = "awh-nsamples-update";
    awhParams->numSamplesUpdateFreeEnergy = get_eint(inp, opt, 10, wi);

    printStringNoNewline(
            inp, "When true, biases with share-group>0 are shared between multiple simulations");
    opt                          = "awh-share-multisim";
    awhParams->shareBiasMultisim = (get_eeenum(inp, opt, yesno_names, wi) != 0);

    printStringNoNewline(inp, "The number of independent AWH biases");
    opt                = "awh-nbias";
    awhParams->numBias = get_eint(inp, opt, 1, wi);
    /* Check this before starting to read the AWH biases. */
    if (awhParams->numBias <= 0)
    {
        gmx_fatal(FARGS, "%s needs to be an integer > 0", opt.c_str());
    }

    /* Read the parameters specific to each AWH bias */
    snew(awhParams->awhBiasParams, awhParams->numBias);

    for (int k = 0; k < awhParams->numBias; k++)
    {
        bool        bComment  = (k == 0);
        std::string prefixawh = formatString("awh%d", k + 1);
        readBiasParams(inp, &awhParams->awhBiasParams[k], prefixawh, wi, bComment);
    }

    return awhParams;
}

void checkAwhParams(const AwhParams* awhParams, const t_inputrec* ir, warninp_t wi)
{
    std::string opt;

    if (!ir->bPull)
    {
        gmx_fatal(FARGS, "AWH biasing is only compatible with COM pulling turned on");
    }

    opt = "awh-nstout";
    if (awhParams->nstOut <= 0)
    {
        auto message = formatString("Not writing AWH output with AWH (%s = %d) does not make sense",
                                    opt.c_str(), awhParams->nstOut);
        warning_error(wi, message);
    }
    /* This restriction can be removed by changing a flag of print_ebin() */
    if (ir->nstenergy == 0 || awhParams->nstOut % ir->nstenergy != 0)
    {
        auto message = formatString("%s (%d) should be a multiple of nstenergy (%d)", opt.c_str(),
                                    awhParams->nstOut, ir->nstenergy);
        warning_error(wi, message);
    }

    opt = "awh-nsamples-update";
    if (awhParams->numSamplesUpdateFreeEnergy <= 0)
    {
        warning_error(wi, opt + " needs to be an integer > 0");
    }

    for (int k = 0; k < awhParams->numBias; k++)
    {
        std::string prefixawh = formatString("awh%d", k + 1);
        checkBiasParams(&awhParams->awhBiasParams[k], prefixawh, ir, wi);
    }

    /* Do a final consistency check before returning */
    checkInputConsistencyAwh(*awhParams, wi);

    if (ir->init_step != 0)
    {
        warning_error(wi, "With AWH init-step should be 0");
    }
}

/*! \brief
 * Gets the period of a pull coordinate.
 *
 * \param[in] pullCoordParams   The parameters for the pull coordinate.
 * \param[in] pbc               The PBC setup
 * \param[in] intervalLength    The length of the AWH interval for this pull coordinate
 * \returns the period (or 0 if not periodic).
 */
static double get_pull_coord_period(const t_pull_coord& pullCoordParams, const t_pbc& pbc, const real intervalLength)
{
    double period = 0;

    if (pullCoordParams.eGeom == epullgDIR)
    {
        const real margin = 0.001;
        // Make dims periodic when the interval covers > 95%
        const real periodicFraction = 0.95;

        // Check if the pull direction is along a box vector
        for (int dim = 0; dim < pbc.ndim_ePBC; dim++)
        {
            const real boxLength    = norm(pbc.box[dim]);
            const real innerProduct = iprod(pullCoordParams.vec, pbc.box[dim]);
            if (innerProduct >= (1 - margin) * boxLength && innerProduct <= (1 + margin) * boxLength)
            {
                if (intervalLength > (1 + margin) * boxLength)
                {
                    gmx_fatal(FARGS,
                              "The AWH interval (%f nm) for a pull coordinate is larger than the "
                              "box size (%f nm)",
                              intervalLength, boxLength);
                }

                if (intervalLength > periodicFraction * boxLength)
                {
                    period = boxLength;
                }
            }
        }
    }
    else if (pullCoordParams.eGeom == epullgDIHEDRAL)
    {
        /* The dihedral angle is periodic in -180 to 180 deg */
        period = 360;
    }

    return period;
}

/*! \brief
 * Checks if the given interval is defined in the correct periodic interval.
 *
 * \param[in] origin      Start value of interval.
 * \param[in] end         End value of interval.
 * \param[in] period      Period (or 0 if not periodic).
 * \returns true if the end point values are in the correct periodic interval.
 */
static bool intervalIsInPeriodicInterval(double origin, double end, double period)
{
    return (period == 0) || (std::fabs(origin) <= 0.5 * period && std::fabs(end) <= 0.5 * period);
}

/*! \brief
 * Checks if a value is within an interval.
 *
 * \param[in] origin      Start value of interval.
 * \param[in] end         End value of interval.
 * \param[in] period      Period (or 0 if not periodic).
 * \param[in] value       Value to check.
 * \returns true if the value is within the interval.
 */
static bool valueIsInInterval(double origin, double end, double period, double value)
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
            bIn_interval = ((value >= origin) && (value <= 0.5 * period))
                           || ((value >= -0.5 * period) && (value <= end));
        }
    }
    else
    {
        bIn_interval = (value >= origin) && (value <= end);
    }

    return bIn_interval;
}

/*! \brief
 * Check if the starting configuration is consistent with the given interval.
 *
 * \param[in]     awhParams  AWH parameters.
 * \param[in,out] wi         Struct for bookeeping warnings.
 */
static void checkInputConsistencyInterval(const AwhParams* awhParams, warninp_t wi)
{
    for (int k = 0; k < awhParams->numBias; k++)
    {
        AwhBiasParams* awhBiasParams = &awhParams->awhBiasParams[k];
        for (int d = 0; d < awhBiasParams->ndim; d++)
        {
            AwhDimParams* dimParams  = &awhBiasParams->dimParams[d];
            int           coordIndex = dimParams->coordIndex;
            double origin = dimParams->origin, end = dimParams->end, period = dimParams->period;
            double coordValueInit = dimParams->coordValueInit;

            if ((period == 0) && (origin > end))
            {
                gmx_fatal(FARGS,
                          "For the non-periodic pull coordinates awh%d-dim%d-start (%f) cannot be "
                          "larger than awh%d-dim%d-end (%f)",
                          k + 1, d + 1, origin, k + 1, d + 1, end);
            }

            /* Currently we assume symmetric periodic intervals, meaning we use [-period/2, period/2] as the reference interval.
               Make sure the AWH interval is within this reference interval.

               Note: we could fairly simply allow using a  more general interval (e.g. [x, x + period]) but it complicates
               things slightly and I don't see that there is a great need for it. It would also mean that the interval would
               depend on AWH input. Also, for dihedral angles you would always want the reference interval to be -180, +180,
               independent of AWH parameters.
             */
            if (!intervalIsInPeriodicInterval(origin, end, period))
            {
                gmx_fatal(FARGS,
                          "When using AWH with periodic pull coordinate geometries "
                          "awh%d-dim%d-start (%.8g) and "
                          "awh%d-dim%d-end (%.8g) should cover at most one period (%.8g) and take "
                          "values in between "
                          "minus half a period and plus half a period, i.e. in the interval [%.8g, "
                          "%.8g].",
                          k + 1, d + 1, origin, k + 1, d + 1, end, period, -0.5 * period, 0.5 * period);
            }

            /* Warn if the pull initial coordinate value is not in the grid */
            if (!valueIsInInterval(origin, end, period, coordValueInit))
            {
                auto message = formatString(
                        "The initial coordinate value (%.8g) for pull coordinate index %d falls "
                        "outside "
                        "of the sampling nterval awh%d-dim%d-start (%.8g) to awh%d-dim%d-end "
                        "(%.8g). "
                        "This can lead to large initial forces pulling the coordinate towards the "
                        "sampling interval.",
                        coordValueInit, coordIndex + 1, k + 1, d + 1, origin, k + 1, d + 1, end);
                warning(wi, message);
            }
        }
    }
}

void setStateDependentAwhParams(AwhParams*           awhParams,
                                const pull_params_t* pull_params,
                                pull_t*              pull_work,
                                const matrix         box,
                                PbcType              pbcType,
                                const tensor&        compressibility,
                                const t_grpopts*     inputrecGroupOptions,
                                warninp_t            wi)
{
    /* The temperature is not really state depenendent but is not known
     * when read_awhParams is called (in get ir).
     * It is known first after do_index has been called in grompp.cpp.
     */
    if (inputrecGroupOptions->ref_t == nullptr || inputrecGroupOptions->ref_t[0] <= 0)
    {
        gmx_fatal(FARGS, "AWH biasing is only supported for temperatures > 0");
    }
    for (int i = 1; i < inputrecGroupOptions->ngtc; i++)
    {
        if (inputrecGroupOptions->ref_t[i] != inputrecGroupOptions->ref_t[0])
        {
            gmx_fatal(FARGS,
                      "AWH biasing is currently only supported for identical temperatures for all "
                      "temperature coupling groups");
        }
    }

    t_pbc pbc;
    set_pbc(&pbc, pbcType, box);

    for (int k = 0; k < awhParams->numBias; k++)
    {
        AwhBiasParams* awhBiasParams = &awhParams->awhBiasParams[k];
        for (int d = 0; d < awhBiasParams->ndim; d++)
        {
            AwhDimParams*       dimParams       = &awhBiasParams->dimParams[d];
            const t_pull_coord& pullCoordParams = pull_params->coord[dimParams->coordIndex];

            if (pullCoordParams.eGeom == epullgDIRPBC)
            {
                gmx_fatal(FARGS,
                          "AWH does not support pull geometry '%s'. "
                          "If the maximum distance between the groups is always less than half the "
                          "box size, "
                          "you can use geometry '%s' instead.",
                          EPULLGEOM(epullgDIRPBC), EPULLGEOM(epullgDIR));
            }

            dimParams->period =
                    get_pull_coord_period(pullCoordParams, pbc, dimParams->end - dimParams->origin);
            // We would like to check for scaling, but we don't have the full inputrec available here
            if (dimParams->period > 0
                && !(pullCoordParams.eGeom == epullgANGLE || pullCoordParams.eGeom == epullgDIHEDRAL))
            {
                bool coordIsScaled = false;
                for (int d2 = 0; d2 < DIM; d2++)
                {
                    if (pullCoordParams.vec[d2] != 0 && norm2(compressibility[d2]) != 0)
                    {
                        coordIsScaled = true;
                    }
                }
                if (coordIsScaled)
                {
                    std::string mesg = gmx::formatString(
                            "AWH dimension %d of bias %d is periodic with pull geometry '%s', "
                            "while you should are applying pressure scaling to the corresponding "
                            "box vector, this is not supported.",
                            d + 1, k + 1, EPULLGEOM(pullCoordParams.eGeom));
                    warning(wi, mesg.c_str());
                }
            }

            /* The initial coordinate value, converted to external user units. */
            dimParams->coordValueInit = get_pull_coord_value(pull_work, dimParams->coordIndex, &pbc);

            dimParams->coordValueInit *= pull_conversion_factor_internal2userinput(&pullCoordParams);
        }
    }
    checkInputConsistencyInterval(awhParams, wi);

    /* Register AWH as external potential with pull to check consistency. */
    Awh::registerAwhWithPull(*awhParams, pull_work);
}

} // namespace gmx
