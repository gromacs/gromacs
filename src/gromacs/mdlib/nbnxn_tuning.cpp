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
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
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

        /* Dynamic pruning is only useful when the inner list is smaller
         * than the outer. The factor 0.99 ensures at least 3% list reduction.
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

/* Set up the dynamic pairlist pruning */
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
            if (useGpu && listParams->nstlistPrune % c_nbnxnGpuRollingListPruningInterval != 0)
            {
                gmx_fatal(FARGS, "With a GPU, GMX_NSTLIST_DYNAMICPRUNING should be divisble by %d", c_nbnxnGpuRollingListPruningInterval);
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
            GMX_RELEASE_ASSERT(listParams->nstlistPrune % c_nbnxnGpuRollingListPruningInterval == 0, "For efficiency reasons, and correctness of numRollingParts below, nstlistPrune should be divisible by c_nbnxnGpuRollingListPruningInterval");
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
