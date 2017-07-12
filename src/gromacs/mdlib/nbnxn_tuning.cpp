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
const int           nstlist_try[] = { 20, 25, 40 };
//! Number of elements in the neighborsearch list trials.
#define NNSTL  sizeof(nstlist_try)/sizeof(nstlist_try[0])
/* Increase nstlist until the non-bonded cost increases more than listfac_ok,
 * but never more than listfac_max.
 * A standard (protein+)water system at 300K with PME ewald_rtol=1e-5
 * needs 1.28 at rcoulomb=0.9 and 1.24 at rcoulomb=1.0 to get to nstlist=40.
 * Note that both CPU and GPU factors are conservative. Performance should
 * not go down due to this tuning, except with a relatively slow GPU.
 * On the other hand, at medium/high parallelization or with fast GPUs
 * nstlist will not be increased enough to reach optimal performance.
 */
/* CPU: pair-search is about a factor 1.5 slower than the non-bonded kernel */
//! Max OK performance ratio beween force calc and neighbor searching
static const float  nbnxn_cpu_listfac_ok    = 1.05;
//! Too high performance ratio beween force calc and neighbor searching
static const float  nbnxn_cpu_listfac_max   = 1.09;
/* CPU: pair-search is about a factor 2-3 slower than the non-bonded kernel */
//! Max OK performance ratio beween force calc and neighbor searching
static const float  nbnxn_knl_listfac_ok    = 1.22;
//! Too high performance ratio beween force calc and neighbor searching
static const float  nbnxn_knl_listfac_max   = 1.3;
/* GPU: pair-search is a factor 1.5-3 slower than the non-bonded kernel */
//! Max OK performance ratio beween force calc and neighbor searching
static const float  nbnxn_gpu_listfac_ok    = 1.20;
//! Too high performance ratio beween force calc and neighbor searching
static const float  nbnxn_gpu_listfac_max   = 1.30;

void increase_nstlist(FILE *fp, t_commrec *cr,
                      t_inputrec *ir, int nstlist_cmdline,
                      const gmx_mtop_t *mtop, matrix box,
                      bool makeGpuPairList, const gmx::CpuInfo &cpuinfo)
{
    float                  listfac_ok, listfac_max;
    int                    nstlist_orig, nstlist_prev;
    verletbuf_list_setup_t ls;
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

        if (fp != nullptr && makeGpuPairList && ir->nstlist < nstlist_try[0])
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

    if (ir->verletbuf_tol == 0 && makeGpuPairList)
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

    if (makeGpuPairList)
    {
        listfac_ok  = nbnxn_gpu_listfac_ok;
        listfac_max = nbnxn_gpu_listfac_max;
    }
    else if (cpuinfo.feature(gmx::CpuInfo::Feature::X86_Avx512ER))
    {
        listfac_ok  = nbnxn_knl_listfac_ok;
        listfac_max = nbnxn_knl_listfac_max;
    }
    else
    {
        listfac_ok  = nbnxn_cpu_listfac_ok;
        listfac_max = nbnxn_cpu_listfac_max;
    }

    nstlist_orig = ir->nstlist;
    if (nstlist_cmdline > 0)
    {
        if (fp)
        {
            sprintf(buf, "Getting nstlist=%d from command line option",
                    nstlist_cmdline);
        }
        ir->nstlist = nstlist_cmdline;
    }

    verletbuf_get_list_setup(true, makeGpuPairList, &ls);

    /* Allow rlist to make the list a given factor larger than the list
     * would be with the reference value for nstlist (10).
     */
    nstlist_prev = ir->nstlist;
    ir->nstlist  = nbnxnReferenceNstlist;
    calc_verlet_buffer_size(mtop, det(box), ir, ir->nstlist, ir->nstlist - 1,
                            -1, &ls, nullptr, &rlistWithReferenceNstlist);
    ir->nstlist  = nstlist_prev;

    /* Determine the pair list size increase due to zero interactions */
    rlist_inc = nbnxn_get_rlist_effective_inc(ls.cluster_size_j,
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
        calc_verlet_buffer_size(mtop, det(box), ir, ir->nstlist, ir->nstlist - 1, -1, &ls, nullptr, &rlist_new);

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
            fprintf(fp, "\n%s\n", bBox ? box_err : dd_err);
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
                                    const verletbuf_list_setup_t &listSetup,
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

void setupDynamicPairlistPruning(FILE                      *fplog,
                                 const t_inputrec          *ir,
                                 const gmx_mtop_t          *mtop,
                                 matrix                     box,
                                 bool                       useGpu,
                                 const interaction_const_t *ic,
                                 NbnxnListParameters       *listParams)
{
    GMX_RELEASE_ASSERT(listParams->rlistOuter > 0, "With the nbnxn setup rlist should be > 0");

    /* Initialize the parameters to no dynamic list pruning */
    listParams->useDynamicPruning = false;

    if (supportsDynamicPairlistGenerationInterval(*ir) &&
        getenv("GMX_DISABLE_DYNAMICPRUNING") == NULL)
    {
        verletbuf_list_setup_t ls;
        verletbuf_get_list_setup(TRUE, TRUE, &ls);

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

        // Dynamic pruning disabled until the kernels are present
        listParams->useDynamicPruning = false;

        if (listParams->useDynamicPruning && useGpu)
        {
            /* Note that we can round down here. This makes the effective
             * rolling pruning interval slightly shorter than nstlistTune,
             * thus giving correct results, but a slightly lower efficiency.
             */
            listParams->numRollingParts = listParams->nstlistPrune/c_nbnxnGpuRollingListPruningInterval;
        }
        else
        {
            listParams->numRollingParts = 1;
        }

        if (fplog && listParams->useDynamicPruning)
        {
            const real interactionCutoff = std::max(ic->rcoulomb, ic->rvdw);
            fprintf(fplog,
                    "Using a dual pair-list setup updated with dynamic%s pruning:\n"
                    "  outer list: updated every %3d steps, buffer %.3f nm, rlist %.3f nm\n"
                    "  inner list: updated every %3d steps, buffer %.3f nm, rlist %.3f nm\n",
                    listParams->numRollingParts > 1 ? ", rolling" : "",
                    ir->nstlist, listParams->rlistOuter - interactionCutoff, listParams->rlistOuter,
                    listParams->nstlistPrune, listParams->rlistInner - interactionCutoff, listParams->rlistInner);
        }
    }
}
