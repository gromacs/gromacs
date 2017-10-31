/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
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
 *
 * \brief Implements functions for tuning adjustable parameters for the nbnxn non-bonded search and interaction kernels
 *
 * \author Berk Hess <hess@kth.se>
 * \ingroup __module_nb_verlet
 */

#include "gmxpre.h"

#include "nbnxn_tuning.h"

#include <assert.h>
#include <stdlib.h>

#include <cmath>

#include <algorithm>
#include <string>

#include "gromacs/domdec/domdec.h"
#include "gromacs/hardware/cpuinfo.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/calc_verletbuf.h"
#include "gromacs/mdlib/nb_verlet.h"
#include "gromacs/mdlib/nbnxn_search.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/stringutil.h"

/*! \brief Returns if we can (heuristically) change nstlist and rlist
 *
 * \param [in] ir  The input parameter record
 */
static bool supportsDynamicPairlistGenerationInterval(const t_inputrec &ir)
{
    return
        ir.cutoff_scheme == ecutsVERLET &&
        EI_DYNAMICS(ir.eI) &&
        !(EI_MD(ir.eI) && ir.etc == etcNO) &&
        ir.verletbuf_tol > 0;
}

/*! \brief Cost of non-bonded kernels
 *
 * We determine the extra cost of the non-bonded kernels compared to
 * a reference nstlist value of 10 (which is the default in grompp).
 */
static const int    nbnxnReferenceNstlist = 10;
//! The values to try when switching
const int           nstlist_try[] = { 20, 25, 40, 50, 80, 100 };
//! Number of elements in the neighborsearch list trials.
#define NNSTL  sizeof(nstlist_try)/sizeof(nstlist_try[0])
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
static const float c_nbnxnListSizeFactorCpu           = 1.25;
// Intel KNL: pair-search is a factor ~2-3 slower than the non-bonded kernel.
//! Target pair-list size increase ratio for Intel KNL
static const float c_nbnxnListSizeFactorIntelXeonPhi  = 1.4;
// GPU: pair-search is a factor 1.5-3 slower than the non-bonded kernel.
//! Target pair-list size increase ratio for GPU
static const float c_nbnxnListSizeFactorGPU           = 1.4;
//! Never increase the size of the pair-list more than the factor above plus this margin
static const float c_nbnxnListSizeFactorMargin        = 0.1;

void increaseNstlist(FILE *fp, t_commrec *cr,
                     t_inputrec *ir, int nstlist_cmdline,
                     const gmx_mtop_t *mtop,
                     const matrix box,
                     bool useOrEmulateGpuForNonbondeds,
                     const gmx::CpuInfo &cpuinfo)
{
    if (!EI_DYNAMICS(ir->eI))
    {
        /* Can only increase nstlist with dynamics */
        return;
    }

    float                  listfac_ok, listfac_max;
    int                    nstlist_orig, nstlist_prev;
    real                   rlistWithReferenceNstlist, rlist_inc, rlist_ok, rlist_max;
    real                   rlist_new, rlist_prev;
    size_t                 nstlist_ind = 0;
    gmx_bool               bBox, bDD, bCont;
    const char            *nstl_gpu = "\nFor optimal performance with a GPU nstlist (now %d) should be larger.\nThe optimum depends on your CPU and GPU resources.\nYou might want to try several nstlist values.\n";
    const char            *nve_err  = "Can not increase nstlist because an NVE ensemble is used";
    const char            *vbd_err  = "Can not increase nstlist because verlet-buffer-tolerance is not set or used";
    const char            *box_err  = "Can not increase nstlist because the box is too small";
    const char            *dd_err   = "Can not increase nstlist because of domain decomposition limitations";
    char                   buf[STRLEN];

    if (nstlist_cmdline <= 0)
    {
        if (ir->nstlist == 1)
        {
            /* The user probably set nstlist=1 for a reason,
             * don't mess with the settings.
             */
            return;
        }

        /* With a GPU and fixed nstlist suggest tuning nstlist */
        if (fp != nullptr &&
            useOrEmulateGpuForNonbondeds &&
            ir->nstlist < nstlist_try[0] &&
            !supportsDynamicPairlistGenerationInterval(*ir))
        {
            fprintf(fp, nstl_gpu, ir->nstlist);
        }

        nstlist_ind = 0;
        while (nstlist_ind < NNSTL && ir->nstlist >= nstlist_try[nstlist_ind])
        {
            nstlist_ind++;
        }
        if (nstlist_ind == NNSTL)
        {
            /* There are no larger nstlist value to try */
            return;
        }
    }

    if (EI_MD(ir->eI) && ir->etc == etcNO)
    {
        if (MASTER(cr))
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
        gmx_fatal(FARGS, "You are using an old tpr file with a GPU, please generate a new tpr file with an up to date version of grompp");
    }

    if (ir->verletbuf_tol < 0)
    {
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n", vbd_err);
        }
        if (fp != nullptr)
        {
            fprintf(fp, "%s\n", vbd_err);
        }

        return;
    }

    GMX_RELEASE_ASSERT(supportsDynamicPairlistGenerationInterval(*ir), "In all cases that do not support dynamic nstlist, we should have returned with an appropriate message above");

    if (useOrEmulateGpuForNonbondeds)
    {
        listfac_ok  = c_nbnxnListSizeFactorGPU;
    }
    else if (cpuinfo.brandString().find("Xeon Phi") != std::string::npos)
    {
        listfac_ok  = c_nbnxnListSizeFactorIntelXeonPhi;
    }
    else
    {
        listfac_ok  = c_nbnxnListSizeFactorCpu;
    }
    listfac_max     = listfac_ok + c_nbnxnListSizeFactorMargin;

    nstlist_orig    = ir->nstlist;
    if (nstlist_cmdline > 0)
    {
        if (fp)
        {
            sprintf(buf, "Getting nstlist=%d from command line option",
                    nstlist_cmdline);
        }
        ir->nstlist = nstlist_cmdline;
    }

    ListSetupType      listType  = (useOrEmulateGpuForNonbondeds ? ListSetupType::Gpu : ListSetupType::CpuSimdWhenSupported);
    VerletbufListSetup listSetup = verletbufGetSafeListSetup(listType);

    /* Allow rlist to make the list a given factor larger than the list
     * would be with the reference value for nstlist (10).
     */
    nstlist_prev = ir->nstlist;
    ir->nstlist  = nbnxnReferenceNstlist;
    calc_verlet_buffer_size(mtop, det(box), ir, ir->nstlist, ir->nstlist - 1,
                            -1, &listSetup,
                            nullptr, &rlistWithReferenceNstlist);
    ir->nstlist  = nstlist_prev;

    /* Determine the pair list size increase due to zero interactions */
    rlist_inc = nbnxn_get_rlist_effective_inc(listSetup.cluster_size_j,
                                              mtop->natoms/det(box));
    rlist_ok  = (rlistWithReferenceNstlist + rlist_inc)*std::cbrt(listfac_ok) - rlist_inc;
    rlist_max = (rlistWithReferenceNstlist + rlist_inc)*std::cbrt(listfac_max) - rlist_inc;
    if (debug)
    {
        fprintf(debug, "nstlist tuning: rlist_inc %.3f rlist_ok %.3f rlist_max %.3f\n",
                rlist_inc, rlist_ok, rlist_max);
    }

    nstlist_prev = nstlist_orig;
    rlist_prev   = ir->rlist;
    do
    {
        if (nstlist_cmdline <= 0)
        {
            ir->nstlist = nstlist_try[nstlist_ind];
        }

        /* Set the pair-list buffer size in ir */
        calc_verlet_buffer_size(mtop, det(box), ir, ir->nstlist, ir->nstlist - 1, -1, &listSetup, nullptr, &rlist_new);

        /* Does rlist fit in the box? */
        bBox = (gmx::square(rlist_new) < max_cutoff2(ir->ePBC, box));
        bDD  = TRUE;
        if (bBox && DOMAINDECOMP(cr))
        {
            /* Check if rlist fits in the domain decomposition */
            if (inputrec2nboundeddim(ir) < DIM)
            {
                gmx_incons("Changing nstlist with domain decomposition and unbounded dimensions is not implemented yet");
            }
            t_state state_tmp;
            copy_mat(box, state_tmp.box);
            bDD = change_dd_cutoff(cr, &state_tmp, ir, rlist_new);
        }

        if (debug)
        {
            fprintf(debug, "nstlist %d rlist %.3f bBox %d bDD %d\n",
                    ir->nstlist, rlist_new, bBox, bDD);
        }

        bCont = FALSE;

        if (nstlist_cmdline <= 0)
        {
            if (bBox && bDD && rlist_new <= rlist_max)
            {
                /* Increase nstlist */
                nstlist_prev = ir->nstlist;
                rlist_prev   = rlist_new;
                bCont        = (nstlist_ind+1 < NNSTL && rlist_new < rlist_ok);
            }
            else
            {
                /* Stick with the previous nstlist */
                ir->nstlist = nstlist_prev;
                rlist_new   = rlist_prev;
                bBox        = TRUE;
                bDD         = TRUE;
            }
        }

        nstlist_ind++;
    }
    while (bCont);

    if (!bBox || !bDD)
    {
        gmx_warning(!bBox ? box_err : dd_err);
        if (fp != nullptr)
        {
            fprintf(fp, "\n%s\n", !bBox ? box_err : dd_err);
        }
        ir->nstlist = nstlist_orig;
    }
    else if (ir->nstlist != nstlist_orig || rlist_new != ir->rlist)
    {
        sprintf(buf, "Changing nstlist from %d to %d, rlist from %g to %g",
                nstlist_orig, ir->nstlist,
                ir->rlist, rlist_new);
        if (MASTER(cr))
        {
            fprintf(stderr, "%s\n\n", buf);
        }
        if (fp != nullptr)
        {
            fprintf(fp, "%s\n\n", buf);
        }
        ir->rlist     = rlist_new;
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

/*! \brief The minimum nstlist for dynamic pair list pruning.
 *
 * In most cases going lower than 4 will lead to a too high pruning cost.
 * This value should be a multiple of \p c_nbnxnGpuRollingListPruningInterval
 */
static const int c_nbnxnDynamicListPruningMinLifetime = 4;

/*! \brief Set the dynamic pairlist pruning parameters in \p ic
 *
 * \param[in]     ir          The input parameter record
 * \param[in]     mtop        The global topology
 * \param[in]     box         The unit cell
 * \param[in]     useGpu      Tells if we are using a GPU for non-bondeds
 * \param[in]     listSetup   The nbnxn pair list setup
 * \param[in]     userSetNstlistPrune  The user set ic->nstlistPrune (using an env.var.)
 * \param[in] ic              The nonbonded interactions constants
 * \param[in,out] listParams  The list setup parameters
 */
static void
setDynamicPairlistPruningParameters(const t_inputrec             *ir,
                                    const gmx_mtop_t             *mtop,
                                    matrix                        box,
                                    gmx_bool                      useGpu,
                                    const VerletbufListSetup     &listSetup,
                                    bool                          userSetNstlistPrune,
                                    const interaction_const_t    *ic,
                                    NbnxnListParameters          *listParams)
{
    /* When nstlistPrune was set by the user, we need to execute one loop
     * iteration to determine rlistInner.
     * Otherwise we compute rlistInner and increase nstlist as long as
     * we have a pairlist buffer of length 0 (i.e. rlistInner == cutoff).
     */
    const real interactionCutoff = std::max(ic->rcoulomb, ic->rvdw);
    int        tunedNstlistPrune = listParams->nstlistPrune;
    do
    {
        /* Dynamic pruning on the GPU is performed on the list for
         * the next step on the coordinates of the current step,
         * so the list lifetime is nstlistPrune (not the usual nstlist-1).
         */
        int listLifetime         = tunedNstlistPrune - (useGpu ? 0 : 1);
        listParams->nstlistPrune = tunedNstlistPrune;
        calc_verlet_buffer_size(mtop, det(box), ir,
                                tunedNstlistPrune, listLifetime,
                                -1, &listSetup, NULL,
                                &listParams->rlistInner);

        /* On the GPU we apply the dynamic pruning in a rolling fashion
         * every c_nbnxnGpuRollingListPruningInterval steps,
         * so keep nstlistPrune a multiple of the interval.
         */
        tunedNstlistPrune += useGpu ? c_nbnxnGpuRollingListPruningInterval : 1;
    }
    while (!userSetNstlistPrune &&
           tunedNstlistPrune < ir->nstlist &&
           listParams->rlistInner == interactionCutoff);

    if (userSetNstlistPrune)
    {
        listParams->useDynamicPruning = true;
    }
    else
    {
        /* Determine the pair list size increase due to zero interactions */
        real rlistInc = nbnxn_get_rlist_effective_inc(listSetup.cluster_size_j,
                                                      mtop->natoms/det(box));

        /* Dynamic pruning is only useful when the inner list is smaller than
         * the outer. The factor 0.99 ensures at least 3% list size reduction.
         *
         * With dynamic pruning on the CPU we prune after updating,
         * so nstlistPrune=nstlist-1 would add useless extra work.
         * With the GPU there will probably be more overhead than gain
         * with nstlistPrune=nstlist-1, so we disable dynamic pruning.
         * Note that in such cases the first sub-condition is likely also false.
         */
        listParams->useDynamicPruning =
            (listParams->rlistInner + rlistInc < 0.99*(listParams->rlistOuter + rlistInc) &&
             listParams->nstlistPrune < ir->nstlist - 1);
    }

    if (!listParams->useDynamicPruning)
    {
        /* These parameters should not be used, but set them to useful values */
        listParams->nstlistPrune  = -1;
        listParams->rlistInner    = listParams->rlistOuter;
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
static std::string formatListSetup(const std::string &listName,
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
    std::string nstListFormat = "%" + gmx::formatString("%d", gmx::formatString("%zu", nstListForSpacing).size()) + "d";
    listSetup += gmx::formatString(nstListFormat.c_str(), nstList);
    listSetup += gmx::formatString(" steps, buffer %.3f nm, rlist %.3f nm\n",
                                   rList - interactionCutoff, rList);

    return listSetup;
}

void setupDynamicPairlistPruning(const gmx::MDLogger       &mdlog,
                                 const t_inputrec          *ir,
                                 const gmx_mtop_t          *mtop,
                                 matrix                     box,
                                 int                        nbnxnKernelType,
                                 const interaction_const_t *ic,
                                 NbnxnListParameters       *listParams)
{
    GMX_RELEASE_ASSERT(listParams->rlistOuter > 0, "With the nbnxn setup rlist should be > 0");

    /* Initialize the parameters to no dynamic list pruning */
    listParams->useDynamicPruning = false;

    const VerletbufListSetup ls   = verletbufGetListSetup(nbnxnKernelType);

    /* Currently emulation mode does not support dual pair-lists */
    const bool useGpu             = (nbnxnKernelType == nbnxnk8x8x8_GPU);

    if (supportsDynamicPairlistGenerationInterval(*ir) &&
        getenv("GMX_DISABLE_DYNAMICPRUNING") == NULL)
    {
        /* Note that nstlistPrune can have any value independently of nstlist.
         * Actually applying rolling pruning is only useful when
         * nstlistPrune < nstlist -1
         */
        char *env                 = getenv("GMX_NSTLIST_DYNAMICPRUNING");
        bool  userSetNstlistPrune = (env != NULL);

        if (userSetNstlistPrune)
        {
            char *end;
            listParams->nstlistPrune = strtol(env, &end, 10);
            if (!end || (*end != 0) ||
                !(listParams->nstlistPrune > 0 && listParams->nstlistPrune < ir->nstlist))
            {
                gmx_fatal(FARGS, "Invalid value passed in GMX_NSTLIST_DYNAMICPRUNING=%s, should be > 0 and < nstlist", env);
            }
        }
        else
        {
            static_assert(c_nbnxnDynamicListPruningMinLifetime % c_nbnxnGpuRollingListPruningInterval == 0,
                          "c_nbnxnDynamicListPruningMinLifetime sets the starting value for nstlistPrune, which should be divisible by the rolling pruning interval for efficiency reasons.");

            // TODO: Use auto-tuning to determine nstlistPrune
            listParams->nstlistPrune = c_nbnxnDynamicListPruningMinLifetime;
        }

        setDynamicPairlistPruningParameters(ir, mtop, box, useGpu, ls,
                                            userSetNstlistPrune, ic,
                                            listParams);

        if (listParams->useDynamicPruning && useGpu)
        {
            /* Note that we can round down here. This makes the effective
             * rolling pruning interval slightly shorter than nstlistTune,
             * thus giving correct results, but a slightly lower efficiency.
             */
            GMX_RELEASE_ASSERT(listParams->nstlistPrune >= c_nbnxnGpuRollingListPruningInterval,
                               ( "With dynamic list pruning on GPUs pruning frequency must be at least as large as the rolling pruning interval (" +
                                 std::to_string(c_nbnxnGpuRollingListPruningInterval) +
                                 ").").c_str() );
            listParams->numRollingParts = listParams->nstlistPrune/c_nbnxnGpuRollingListPruningInterval;
        }
        else
        {
            listParams->numRollingParts = 1;
        }
    }

    std::string mesg;

    const real  interactionCutoff = std::max(ic->rcoulomb, ic->rvdw);
    if (listParams->useDynamicPruning)
    {
        mesg += gmx::formatString("Using a dual %dx%d pair-list setup updated with dynamic%s pruning:\n",
                                  ls.cluster_size_i, ls.cluster_size_j,
                                  listParams->numRollingParts > 1 ? ", rolling" : "");
        mesg += formatListSetup("outer", ir->nstlist, ir->nstlist, listParams->rlistOuter, interactionCutoff);
        mesg += formatListSetup("inner", listParams->nstlistPrune, ir->nstlist, listParams->rlistInner, interactionCutoff);
    }
    else
    {
        mesg += gmx::formatString("Using a %dx%d pair-list setup:\n",
                                  ls.cluster_size_i, ls.cluster_size_j);
        mesg += formatListSetup("", ir->nstlist, ir->nstlist, listParams->rlistOuter, interactionCutoff);
    }
    if (supportsDynamicPairlistGenerationInterval(*ir))
    {
        VerletbufListSetup listSetup1x1 = { 1, 1 };
        real               rlistOuter;
        real               rlistInner;
        calc_verlet_buffer_size(mtop, det(box), ir, ir->nstlist, ir->nstlist - 1,
                                -1, &listSetup1x1, nullptr, &rlistOuter);
        if (listParams->useDynamicPruning)
        {
            int listLifeTime = listParams->nstlistPrune - (useGpu ? 0 : 1);
            calc_verlet_buffer_size(mtop, det(box), ir, listParams->nstlistPrune, listLifeTime,
                                    -1, &listSetup1x1, nullptr, &rlistInner);
        }

        mesg += gmx::formatString("At tolerance %g kJ/mol/ps per atom, equivalent classical 1x1 list would be:\n",
                                  ir->verletbuf_tol);
        if (listParams->useDynamicPruning)
        {
            mesg += formatListSetup("outer", ir->nstlist, ir->nstlist, rlistOuter, interactionCutoff);
            mesg += formatListSetup("inner", listParams->nstlistPrune, ir->nstlist, rlistInner, interactionCutoff);
        }
        else
        {
            mesg += formatListSetup("", ir->nstlist, ir->nstlist, rlistOuter, interactionCutoff);
        }
    }

    GMX_LOG(mdlog.info).asParagraph().appendText(mesg);
}
