/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/*! \internal \file
 *
 * \brief Implements functions for tuning adjustable parameters for the nbnxn non-bonded search and interaction kernels
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup module_nbnxm
 */

#include "gmxpre.h"

#include "pairlist_tuning.h"

#include <cassert>
#include <cmath>
#include <cstdlib>

#include <algorithm>
#include <filesystem>
#include <string>

#include "gromacs/domdec/domdec.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/multipletimestepping.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/nbnxm/pairlistparams.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"

#include "nbnxm_geometry.h"
#include "pairlistsets.h"

namespace gmx
{

/*! \brief Returns if we can (heuristically) change nstlist and rlist
 *
 * \param [in] ir  The input parameter record
 */
static bool supportsDynamicPairlistGenerationInterval(const t_inputrec& ir)
{
    return ir.cutoff_scheme == CutoffScheme::Verlet && EI_DYNAMICS(ir.eI)
           && !(EI_MD(ir.eI) && ir.etc == TemperatureCoupling::No) && ir.verletbuf_tol > 0;
}

/*! \brief Cost of non-bonded kernels
 *
 * We determine the extra cost of the non-bonded kernels compared to
 * a reference nstlist value of 10 (which is the default in grompp).
 */
static const int nbnxnReferenceNstlist = 10;
//! The values to try when switching
const int nstlist_try[] = { 20, 25, 40, 50, 80, 100 };
//! Number of elements in the neighborsearch list trials.
#define NNSTL (sizeof(nstlist_try) / sizeof(nstlist_try[0]))
/* Increase nstlist until the size of the pair-list increased by
 * \p c_nbnxnListSizeFactor??? or more, but never more than
 * \p c_nbnxnListSizeFactor??? + \p c_nbnxnListSizeFactorMargin.
 * Since we have dynamic pair list pruning, the force kernel cost depends
 * only very weakly on nstlist. It depends strongly on nstlistPrune.
 * Increasing nstlist mainly affects the cost of the pair search (down due
 * to lower frequency, up due to larger list) and the list pruning kernel.
 * We increase nstlist conservatively with regard to kernel performance.
 * In serial the search cost is not high and thus we don't gain much by
 * increasing nstlist a lot. In parallel the MPI and CPU-GPU communication
 * volume as well as the communication buffer preparation and reduction time
 * increase quickly with rlist and thus nslist. Therefore we should avoid
 * large nstlist, even if that also reduces the domain decomposition cost.
 * With GPUs we perform the dynamic pruning in a rolling fashion and this
 * overlaps with the update on the CPU, which allows even larger nstlist.
 */
// CPU: pair-search is a factor ~1.5 slower than the non-bonded kernel.
//! Target pair-list size increase ratio for CPU
static const float c_nbnxnListSizeFactorCpu = 1.25;
// Intel KNL: pair-search is a factor ~2-3 slower than the non-bonded kernel.
//! Target pair-list size increase ratio for Intel KNL
static const float c_nbnxnListSizeFactorIntelXeonPhi = 1.4;
// GPU: pair-search is a factor 1.5-3 slower than the non-bonded kernel.
//! Target pair-list size increase ratio for GPU
static const float c_nbnxnListSizeFactorGPU = 1.4;
//! Never increase the size of the pair-list more than the factor above plus this margin
static const float c_nbnxnListSizeFactorMargin = 0.1;

//! Returns the Verlet buffer pressure tolerance set by an env.var. or from input
static real getPressureTolerance(const real inputrecVerletBufferPressureTolerance)
{
    const char* pressureToleranceString = getenv("GMX_VERLET_BUFFER_PRESSURE_TOLERANCE");
    real        pressureTolerance       = -1;
    if (pressureToleranceString != nullptr)
    {
        if (inputrecVerletBufferPressureTolerance > 0)
        {
            GMX_THROW(
                    InvalidInputError("GMX_VERLET_BUFFER_PRESSURE_TOLERANCE cannot be used when "
                                      "verlet-buffer-pressure-tolerance is set in the tpr file"));
        }

        pressureTolerance = std::stod(pressureToleranceString);
        if (pressureTolerance <= 0)
        {
            GMX_THROW(InvalidInputError("Max pressure error should be positive"));
        }
    }
    else
    {
        pressureTolerance = inputrecVerletBufferPressureTolerance;
    }

    return pressureTolerance;
}

void increaseNstlist(FILE*             fp,
                     t_commrec*        cr,
                     t_inputrec*       ir,
                     int               nstlist_cmdline,
                     const gmx_mtop_t* mtop,
                     const matrix      box,
                     const real        effectiveAtomDensity,
                     bool              useOrEmulateGpuForNonbondeds,
                     const CpuInfo&    cpuinfo)
{
    if (!EI_DYNAMICS(ir->eI))
    {
        /* Can only increase nstlist with dynamics */
        return;
    }

    size_t      nstlist_ind = 0;
    const char* nstl_gpu =
            "\nFor optimal performance with a GPU nstlist (now %d) should be larger.\nThe "
            "optimum depends on your CPU and GPU resources.\nYou might want to try several "
            "nstlist values.\n";
    const char* nve_err = "Can not increase nstlist because an NVE ensemble is used";
    const char* vbd_err =
            "Can not increase nstlist because verlet-buffer-tolerance is not set or used";
    const char* box_err = "Can not increase nstlist because the box is too small";
    const char* dd_err  = "Can not increase nstlist because of domain decomposition limitations";
    char        buf[STRLEN];

    /* When most of the computation, and in particular the non-bondeds is only
     * performed every ir->mtsFactor steps due to multiple time stepping,
     * we scale all nstlist values by this factor.
     */
    const int mtsFactor = nonbondedMtsFactor(*ir);

    if (nstlist_cmdline <= 0)
    {
        if (ir->nstlist <= mtsFactor)
        {
            /* The user probably set nstlist<=mtsFactor for a reason,
             * don't mess with the settings, except when < mtsFactor.
             */
            ir->nstlist = mtsFactor;

            return;
        }

        /* With a GPU and fixed nstlist suggest tuning nstlist */
        if (fp != nullptr && useOrEmulateGpuForNonbondeds && ir->nstlist < nstlist_try[0] * mtsFactor
            && !supportsDynamicPairlistGenerationInterval(*ir))
        {
            fprintf(fp, nstl_gpu, ir->nstlist);
        }

        nstlist_ind = 0;
        while (nstlist_ind < NNSTL && ir->nstlist >= nstlist_try[nstlist_ind] * mtsFactor)
        {
            nstlist_ind++;
        }
        if (nstlist_ind == NNSTL)
        {
            /* There are no larger nstlist value to try */
            return;
        }
    }

    if (EI_MD(ir->eI) && ir->etc == TemperatureCoupling::No)
    {
        if (MAIN(cr))
        {
            fprintf(stderr, "%s\n", nve_err);
        }
        if (fp != nullptr)
        {
            fprintf(fp, "%s\n", nve_err);
        }

        return;
    }

    if (ir->verletbuf_tol == 0 && useOrEmulateGpuForNonbondeds)
    {
        gmx_fatal(FARGS,
                  "You are using an old tpr file with a GPU, please generate a new tpr file with "
                  "an up to date version of grompp");
    }

    if (ir->verletbuf_tol < 0)
    {
        if (MAIN(cr))
        {
            fprintf(stderr, "%s\n", vbd_err);
        }
        if (fp != nullptr)
        {
            fprintf(fp, "%s\n", vbd_err);
        }

        return;
    }

    GMX_RELEASE_ASSERT(supportsDynamicPairlistGenerationInterval(*ir),
                       "In all cases that do not support dynamic nstlist, we should have returned "
                       "with an appropriate message above");

    const bool  runningOnXeonPhi = (cpuinfo.brandString().find("Xeon Phi") != std::string::npos);
    const float listfac_ok       = useOrEmulateGpuForNonbondeds ? c_nbnxnListSizeFactorGPU
                                   : runningOnXeonPhi           ? c_nbnxnListSizeFactorIntelXeonPhi
                                                                : c_nbnxnListSizeFactorCpu;
    float       listfac_max      = listfac_ok + c_nbnxnListSizeFactorMargin;

    const int nstlist_orig = ir->nstlist;
    if (nstlist_cmdline > 0)
    {
        if (fp)
        {
            sprintf(buf, "Getting nstlist=%d from command line option", nstlist_cmdline);
        }
        ir->nstlist = nstlist_cmdline;
    }

    ListSetupType listType =
            (useOrEmulateGpuForNonbondeds ? ListSetupType::Gpu : ListSetupType::CpuSimdWhenSupported);
    VerletbufListSetup listSetup = verletbufGetSafeListSetup(listType);

    const real pressureTolerance = getPressureTolerance(ir->verletBufferPressureTolerance);

    /* Allow rlist to make the list a given factor larger than the list
     * would be with the reference value for nstlist (10*mtsFactor).
     *
     * We use half the pressure tolerance to have margin for a, potential, inner pair list.
     */
    int nstlist_prev                     = ir->nstlist;
    ir->nstlist                          = nbnxnReferenceNstlist * mtsFactor;
    const real rlistWithReferenceNstlist = calcVerletBufferSize(
            *mtop, effectiveAtomDensity, *ir, 0.5 * pressureTolerance, ir->nstlist, ir->nstlist - 1, -1, listSetup);
    ir->nstlist = nstlist_prev;

    /* Determine the pair list size increase due to zero interactions */
    real rlist_inc = nbnxmPairlistVolumeRadiusIncrease(useOrEmulateGpuForNonbondeds, effectiveAtomDensity);
    real rlist_ok  = (rlistWithReferenceNstlist + rlist_inc) * std::cbrt(listfac_ok) - rlist_inc;
    real rlist_max = (rlistWithReferenceNstlist + rlist_inc) * std::cbrt(listfac_max) - rlist_inc;
    if (debug)
    {
        fprintf(debug, "nstlist tuning: rlist_inc %.3f rlist_ok %.3f rlist_max %.3f\n", rlist_inc, rlist_ok, rlist_max);
    }

    nstlist_prev    = nstlist_orig;
    real rlist_prev = ir->rlist;
    real rlist_new  = 0;
    bool bBox = false, bDD = false, bCont = false;
    do
    {
        if (nstlist_cmdline <= 0)
        {
            ir->nstlist = nstlist_try[nstlist_ind] * mtsFactor;
        }

        /* Set the pair-list buffer size in ir.
         * We use half the pressure tolerance to have margin for a, potential, inner pair list.
         */
        rlist_new = calcVerletBufferSize(*mtop,
                                         effectiveAtomDensity,
                                         *ir,
                                         0.5 * pressureTolerance,
                                         ir->nstlist,
                                         ir->nstlist - mtsFactor,
                                         -1,
                                         listSetup);

        /* Does rlist fit in the box? */
        bBox = (square(rlist_new) < max_cutoff2(ir->pbcType, box));
        bDD  = true;
        if (bBox && haveDDAtomOrdering(*cr))
        {
            /* Currently (as of July 2020), the code in this if clause is never executed.
             * increaseNstlist(...) is only called from prepare_verlet_scheme, which in turns
             * gets called by the runner _before_ setting up DD. haveDDAtomOrdering(*cr) will
             * therefore always be false here. See #3334.
             */
            /* Check if rlist fits in the domain decomposition */
            if (inputrec2nboundeddim(ir) < DIM)
            {
                gmx_incons(
                        "Changing nstlist with domain decomposition and unbounded dimensions is "
                        "not implemented yet");
            }
            // nstlist tuning happens before GPU DD is initialized so we can't check
            // whether the new cutoff would conflict with direct GPU communication.
            const bool checkGpuDdLimitation = false;
            bDD = change_dd_cutoff(cr, box, ArrayRef<const RVec>(), rlist_new, checkGpuDdLimitation);
        }

        if (debug)
        {
            fprintf(debug,
                    "nstlist %d rlist %.3f bBox %s bDD %s\n",
                    ir->nstlist,
                    rlist_new,
                    boolToString(bBox),
                    boolToString(bDD));
        }

        bCont = false;

        if (nstlist_cmdline <= 0)
        {
            if (bBox && bDD && rlist_new <= rlist_max)
            {
                /* Increase nstlist */
                nstlist_prev = ir->nstlist;
                rlist_prev   = rlist_new;
                bCont        = (nstlist_ind + 1 < NNSTL && rlist_new < rlist_ok);
            }
            else
            {
                /* Stick with the previous nstlist */
                ir->nstlist = nstlist_prev;
                rlist_new   = rlist_prev;
                bBox        = true;
                bDD         = true;
            }
        }

        nstlist_ind++;
    } while (bCont);

    if (!bBox || !bDD)
    {
        gmx_warning("%s", !bBox ? box_err : dd_err);
        if (fp != nullptr)
        {
            fprintf(fp, "\n%s\n", !bBox ? box_err : dd_err);
        }
        ir->nstlist = nstlist_orig;
    }
    else if (ir->nstlist != nstlist_orig || rlist_new != ir->rlist)
    {
        sprintf(buf,
                "Changing nstlist from %d to %d, rlist from %g to %g",
                nstlist_orig,
                ir->nstlist,
                ir->rlist,
                rlist_new);
        if (MAIN(cr))
        {
            fprintf(stderr, "%s\n\n", buf);
        }
        if (fp != nullptr)
        {
            fprintf(fp, "%s\n\n", buf);
        }
        ir->rlist = rlist_new;
    }
}

/*! \brief The interval in steps at which we perform dynamic, rolling pruning on a GPU.
 *
 * Ideally we should auto-tune this value.
 * Not considering overheads, 1 would be the ideal value. But 2 seems
 * a reasonable compromise that reduces GPU kernel launch overheads and
 * also avoids inefficiency on large GPUs when pruning small lists.
 * Because with domain decomposition we alternate local/non-local pruning
 * at even/odd steps, which gives a period of 2, this value currenly needs
 * to be 2, which is indirectly asserted when the GPU pruning is dispatched
 * during the force evaluation.
 */
static const int c_nbnxnGpuRollingListPruningInterval = 2;

/*! \brief The minimum nstlist for dynamic pair list pruning on CPUs.
 *
 * In most cases going lower than 5 will lead to a too high pruning cost.
 */
static const int c_nbnxnCpuDynamicListPruningMinLifetime = 5;

/*! \brief The minimum nstlist for dynamic pair list pruning om GPUs.
 *
 * In most cases going lower than 4 will lead to a too high pruning cost.
 * This value should be a multiple of \p c_nbnxnGpuRollingListPruningInterval
 */
static const int c_nbnxnGpuDynamicListPruningMinLifetime = 4;

//! Struct with references for most parameters for calling calcVerletBufferSize()
struct CalcVerletBufferParameters
{
    const gmx_mtop_t&         mtop;
    const real                effectiveAtomDensity;
    const t_inputrec&         inputrec;
    const real                pressureTolerance;
    const VerletbufListSetup& listSetup;
    const bool                useGpuList;
    const int                 mtsFactor;
};

/*! \brief Wrapper for calcVerletBufferSize() for determining the pruning cut-off
 *
 * \param[in] params   References to most parameters for calcVerletBufferSize()
 * \param[in] nstlist  The pruning interval, also used for setting the list lifetime
 * \return The cut-off for pruning the pairlist
 */
static real calcPruneVerletBufferSize(const CalcVerletBufferParameters& params, const int nstlist)
{
    const int listLifetime = nstlist - (params.useGpuList ? 0 : params.mtsFactor);

    return calcVerletBufferSize(params.mtop,
                                params.effectiveAtomDensity,
                                params.inputrec,
                                params.pressureTolerance,
                                nstlist,
                                listLifetime,
                                -1,
                                params.listSetup);
}

/*! \brief Set the dynamic pairlist pruning parameters in \p ic
 *
 * \param[in]     inputrec          The input parameter record
 * \param[in]     mtop        The global topology
 * \param[in]     effectiveAtomDensity  The effective atom density of the system
 * \param[in]     useGpuList  Tells if we are using a GPU type pairlist
 * \param[in]     listSetup   The nbnxn pair list setup
 * \param[in]     userSetNstlistPrune  The user set ic->nstlistPrune (using an env.var.)
 * \param[in]     interactionConst              The nonbonded interactions constants
 * \param[in,out] listParams  The list setup parameters
 */
static void setDynamicPairlistPruningParameters(const t_inputrec&          inputrec,
                                                const gmx_mtop_t&          mtop,
                                                const real                 effectiveAtomDensity,
                                                const bool                 useGpuList,
                                                const VerletbufListSetup&  listSetup,
                                                const bool                 userSetNstlistPrune,
                                                const interaction_const_t& interactionConst,
                                                PairlistParams*            listParams)
{
    /* Note that we do not treat the energy drift and pressure error consistently here.
     * The contributions to the pressure errors are added up for the outer and inner lists.
     * The energy drift of the outer and inner list are required to independently obey
     * the tolerance. This can lead to a slight underestimate of the drift, but the effect
     * is very small as the energy increases linearly with the distance from the cut-off.
     * Summing the drift estimates from the outer and inner list would lead to significant
     * double counting. This is different for the LJ pressure error, as the LJ force is
     * usually a delta function at the cut-off and thus outer and inner list contributions
     * do add up in practice, although not completely.
     */

    real pressureTolerance = getPressureTolerance(inputrec.verletBufferPressureTolerance);
    if (pressureTolerance > 0)
    {
        // The tolerance for the inner list is the total minus the contribution from the outer list
        pressureTolerance -= verletBufferPressureError(
                mtop, effectiveAtomDensity, inputrec, inputrec.nstlist, false, listParams->rlistOuter, listSetup);
    }

    /* When applying multiple time stepping to the non-bonded forces,
     * we only compute them every mtsFactor steps, so all parameters here
     * should be a multiple of mtsFactor.
     */
    listParams->mtsFactor = nonbondedMtsFactor(inputrec);

    const int mtsFactor = listParams->mtsFactor;

    GMX_RELEASE_ASSERT(inputrec.nstlist % mtsFactor == 0,
                       "nstlist should be a multiple of mtsFactor");

    listParams->lifetime = inputrec.nstlist - mtsFactor;

    /* We are now left with determining the following parameters of listParams:
     * - useDynamicPruning
     * - nstlistPrune (left untouched when set by user)
     * - rlistInner
     */

    // Gather (references to) all constant parameters to calcVerletBufferSize() in a struct
    CalcVerletBufferParameters calcBufferParams(
            { mtop, effectiveAtomDensity, inputrec, pressureTolerance, listSetup, useGpuList, mtsFactor });

    if (userSetNstlistPrune)
    {
        listParams->useDynamicPruning = true;
        listParams->rlistInner = calcPruneVerletBufferSize(calcBufferParams, listParams->nstlistPrune);

        return;
    }

    /* We compute rlistInner and increase nstlist as long as we have
     * a pairlist buffer of length 0 (i.e. rlistInner == cutoff).
     */
    const real interactionCutoff = std::max(interactionConst.rcoulomb, interactionConst.rvdw);
    int        nstlistPrune;
    real       rlistInner;
    int        tunedNstlistPrune = listParams->nstlistPrune;
    do
    {
        /* Dynamic pruning on the GPU is performed on the list for
         * the next step on the coordinates of the current step,
         * so the list lifetime is nstlistPrune (not the usual nstlist-mtsFactor).
         */
        nstlistPrune = tunedNstlistPrune;
        rlistInner   = calcPruneVerletBufferSize(calcBufferParams, nstlistPrune);

        /* On the GPU we apply the dynamic pruning in a rolling fashion
         * every c_nbnxnGpuRollingListPruningInterval steps,
         * so keep nstlistPrune a multiple of the interval.
         */
        tunedNstlistPrune += (useGpuList ? c_nbnxnGpuRollingListPruningInterval : 1) * mtsFactor;
    } while (tunedNstlistPrune < inputrec.nstlist && rlistInner == interactionCutoff);

    /* The current nstlistPrune in listParams is in most cases sub-optimal,
     * as it often just increases the buffer from zero to a (small) non-zero value.
     * When pruning on GPU this doesn't matter much as we prune parts every (second) step.
     */
    if (!useGpuList)
    {
        /* Generally nstlistPrune can be lowered without generating more pruning events.
         * Thus here we decrease nstlistPrune to the lowest value that has the same number
         * of pruning events and therefore the same pruning cost.
         */
        const int numPrunings = (inputrec.nstlist + nstlistPrune - 1) / nstlistPrune;
        // Compute the lowest nstlistPrune that has numPrunings pruning steps
        const int lowerNstlistPrune = (inputrec.nstlist + numPrunings - 1) / numPrunings;
        if (lowerNstlistPrune < nstlistPrune)
        {
            nstlistPrune = lowerNstlistPrune;
            rlistInner   = calcPruneVerletBufferSize(calcBufferParams, nstlistPrune);
        }
    }

    /* Determine the pair list size increase due to zero interactions */
    real rlistInc = nbnxmPairlistVolumeRadiusIncrease(useGpuList, effectiveAtomDensity);

    /* Dynamic pruning is only useful when the inner list is smaller than
     * the outer. The factor 0.99 ensures at least 3% list size reduction.
     *
     * With dynamic pruning on the CPU we prune after updating,
     * so nstlistPrune=nstlist-1 would add useless extra work.
     * With the GPU there will probably be more overhead than gain
     * with nstlistPrune=nstlist-1, so we disable dynamic pruning.
     * Note that in such cases the first sub-condition is likely also false.
     */
    listParams->useDynamicPruning = (rlistInner + rlistInc < 0.99 * (listParams->rlistOuter + rlistInc)
                                     && nstlistPrune < listParams->lifetime);

    if (listParams->useDynamicPruning)
    {
        listParams->nstlistPrune = nstlistPrune;
        listParams->rlistInner   = rlistInner;
    }
    else
    {
        /* These parameters should not be used, but set them to useful values */
        listParams->nstlistPrune = -1;
        listParams->rlistInner   = listParams->rlistOuter;
    }
}

/*! \brief Returns a string describing the setup of a single pair-list
 *
 * \param[in] listName           Short name of the list, can be ""
 * \param[in] nstList            The list update interval in steps
 * \param[in] nstListForSpacing  Update interval for setting the number characters for printing \p nstList
 * \param[in] rList              List cut-off radius
 * \param[in] interactionCutoff  The interaction cut-off, use for printing the list buffer size
 */
static std::string formatListSetup(const std::string& listName,
                                   int                nstList,
                                   int                nstListForSpacing,
                                   real               rList,
                                   real               interactionCutoff)
{
    std::string listSetup = "  ";
    if (!listName.empty())
    {
        listSetup += listName + " list: ";
    }
    listSetup += "updated every ";
    // Make the shortest int format string that fits nstListForSpacing
    std::string nstListFormat =
            "%" + formatString("%zu", formatString("%d", nstListForSpacing).size()) + "d";
    listSetup += formatString(nstListFormat.c_str(), nstList);
    listSetup += formatString(" steps, buffer %.3f nm, rlist %.3f nm\n", rList - interactionCutoff, rList);

    return listSetup;
}

void setupDynamicPairlistPruning(const MDLogger&            mdlog,
                                 const t_inputrec&          inputrec,
                                 const gmx_mtop_t&          mtop,
                                 const real                 effectiveAtomDensity,
                                 const interaction_const_t& interactionConst,
                                 PairlistParams*            listParams)
{
    GMX_RELEASE_ASSERT(listParams->rlistOuter > 0, "With the nbnxn setup rlist should be > 0");

    /* Initialize the parameters to no dynamic list pruning */
    listParams->useDynamicPruning = false;

    const VerletbufListSetup ls = { IClusterSizePerListType[listParams->pairlistType],
                                    JClusterSizePerListType[listParams->pairlistType] };

    /* Currently emulation mode does not support dual pair-lists */
    const bool useGpuList = sc_isGpuPairListType[listParams->pairlistType];

    if (supportsDynamicPairlistGenerationInterval(inputrec) && getenv("GMX_DISABLE_DYNAMICPRUNING") == nullptr)
    {
        /* Note that nstlistPrune can have any value independently of nstlist.
         * Actually applying rolling pruning is only useful when
         * nstlistPrune < nstlist -1
         */
        char* env                 = getenv("GMX_NSTLIST_DYNAMICPRUNING");
        bool  userSetNstlistPrune = (env != nullptr);

        if (userSetNstlistPrune)
        {
            char* end                = nullptr;
            listParams->nstlistPrune = strtol(env, &end, 10);
            if (!end || (*end != 0)
                || !(listParams->nstlistPrune > 0 && listParams->nstlistPrune < inputrec.nstlist))
            {
                gmx_fatal(FARGS,
                          "Invalid value passed in GMX_NSTLIST_DYNAMICPRUNING=%s, should be > 0 "
                          "and < nstlist",
                          env);
            }
        }
        else
        {
            static_assert(c_nbnxnGpuDynamicListPruningMinLifetime % c_nbnxnGpuRollingListPruningInterval == 0,
                          "c_nbnxnGpuDynamicListPruningMinLifetime sets the starting value for "
                          "nstlistPrune, which should be divisible by the rolling pruning interval "
                          "for efficiency reasons.");

            // TODO: Use auto-tuning to determine nstlistPrune
            listParams->nstlistPrune = (useGpuList ? c_nbnxnGpuDynamicListPruningMinLifetime
                                                   : c_nbnxnCpuDynamicListPruningMinLifetime);
        }

        setDynamicPairlistPruningParameters(
                inputrec, mtop, effectiveAtomDensity, useGpuList, ls, userSetNstlistPrune, interactionConst, listParams);

        if (listParams->useDynamicPruning && useGpuList)
        {
            /* Note that we can round down here. This makes the effective
             * rolling pruning interval slightly shorter than nstlistTune,
             * thus giving correct results, but a slightly lower efficiency.
             */
            GMX_RELEASE_ASSERT(listParams->nstlistPrune >= c_nbnxnGpuRollingListPruningInterval,
                               ("With dynamic list pruning on GPUs pruning frequency must be at "
                                "least as large as the rolling pruning interval ("
                                + std::to_string(c_nbnxnGpuRollingListPruningInterval) + ").")
                                       .c_str());
            listParams->numRollingPruningParts =
                    listParams->nstlistPrune / c_nbnxnGpuRollingListPruningInterval;
        }
        else
        {
            listParams->numRollingPruningParts = 1;
        }
    }

    std::string mesg;

    const real interactionCutoff = std::max(interactionConst.rcoulomb, interactionConst.rvdw);
    if (listParams->useDynamicPruning)
    {
        mesg += formatString("Using a dual %dx%d pair-list setup updated with dynamic%s pruning:\n",
                             ls.cluster_size_i,
                             ls.cluster_size_j,
                             listParams->numRollingPruningParts > 1 ? ", rolling" : "");
        mesg += formatListSetup(
                "outer", inputrec.nstlist, inputrec.nstlist, listParams->rlistOuter, interactionCutoff);
        mesg += formatListSetup(
                "inner", listParams->nstlistPrune, inputrec.nstlist, listParams->rlistInner, interactionCutoff);
    }
    else
    {
        mesg += formatString("Using a %dx%d pair-list setup:\n", ls.cluster_size_i, ls.cluster_size_j);
        mesg += formatListSetup(
                "", inputrec.nstlist, inputrec.nstlist, listParams->rlistOuter, interactionCutoff);
    }
    if (supportsDynamicPairlistGenerationInterval(inputrec))
    {
        const real pressureTolerance = getPressureTolerance(inputrec.verletBufferPressureTolerance);

        const VerletbufListSetup listSetup1x1 = { 1, 1 };
        const real               rlistOuter   = calcVerletBufferSize(mtop,
                                                     effectiveAtomDensity,
                                                     inputrec,
                                                     pressureTolerance,
                                                     inputrec.nstlist,
                                                     inputrec.nstlist - 1,
                                                     -1,
                                                     listSetup1x1);
        real                     rlistInner   = rlistOuter;
        if (listParams->useDynamicPruning)
        {
            real pressureToleranceInner = -1;
            if (pressureTolerance > 0)
            {
                pressureToleranceInner =
                        pressureTolerance
                        - verletBufferPressureError(
                                mtop, effectiveAtomDensity, inputrec, inputrec.nstlist, false, rlistOuter, listSetup1x1);
            }
            int listLifeTime = listParams->nstlistPrune - (useGpuList ? 0 : 1);
            rlistInner       = calcVerletBufferSize(mtop,
                                              effectiveAtomDensity,
                                              inputrec,
                                              pressureToleranceInner,
                                              listParams->nstlistPrune,
                                              listLifeTime,
                                              -1,
                                              listSetup1x1);
        }

        mesg += formatString(
                "At tolerance %g kJ/mol/ps per atom, equivalent classical 1x1 list would be:\n",
                inputrec.verletbuf_tol);
        if (listParams->useDynamicPruning)
        {
            mesg += formatListSetup(
                    "outer", inputrec.nstlist, inputrec.nstlist, rlistOuter, interactionCutoff);
            mesg += formatListSetup(
                    "inner", listParams->nstlistPrune, inputrec.nstlist, rlistInner, interactionCutoff);
        }
        else
        {
            mesg += formatListSetup("", inputrec.nstlist, inputrec.nstlist, rlistOuter, interactionCutoff);
        }
    }

    GMX_LOG(mdlog.info).asParagraph().appendText(mesg);
}

void printNbnxmPressureError(const MDLogger&       mdlog,
                             const t_inputrec&     inputrec,
                             const gmx_mtop_t&     mtop,
                             const real            effectiveAtomDensity,
                             const PairlistParams& listParams)
{
    const VerletbufListSetup ls = { IClusterSizePerListType[listParams.pairlistType],
                                    JClusterSizePerListType[listParams.pairlistType] };

    real pressureError = verletBufferPressureError(
            mtop, effectiveAtomDensity, inputrec, inputrec.nstlist, false, listParams.rlistOuter, ls);
    if (listParams.useDynamicPruning)
    {
        /* Currently emulation mode does not support dual pair-lists */
        const bool useGpuList = sc_isGpuPairListType[listParams.pairlistType];

        // Add the error due to the pruning of the inner list.
        // The errors are not completely independent, so this results
        // in a (slight) overestimate.
        pressureError += verletBufferPressureError(mtop,
                                                   effectiveAtomDensity,
                                                   inputrec,
                                                   listParams.nstlistPrune,
                                                   useGpuList,
                                                   listParams.rlistInner,
                                                   ls);
    }

    GMX_LOG(mdlog.info)
            .asParagraph()
            .appendText(formatString("The average pressure is off by at most %.2f bar due to "
                                     "missing LJ interactions",
                                     pressureError));
}

} // namespace gmx
