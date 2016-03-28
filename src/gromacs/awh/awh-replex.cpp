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

#include "awh-replex.h"

#include <assert.h>
#include <cmath>
#include <algorithm>

#include "gromacs/awh/awh-grid.h"
#include "gromacs/awh/awh-internal.h"
#include "gromacs/awh/awh-params.h"
#include "gromacs/awh/awh-types.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/qsort_threadsafe.h"
#include "gromacs/utility/smalloc.h"

typedef struct awh_replex_t
{
    int      replica_id;            /* The id of this replica */
    int      nreplica;              /* The number of replicas */
    int      npoints;               /* Number of points in bias array. Currently all replicas need to have the same number of points. */
    bool     all_replicas_exchange; /* True if all replicas are allowed to exchange */
    int     *domain_id;             /* AWH domain id's. Replicas in the same interval do not exchange */
    int     *point_index;           /* Current point index for each replica */
    double **point_bias;            /* Current bias for each point and replica */
} awh_replex_t;

static bool replicas_are_equivalent(const awh_replex_t *replex, int replica_id1, int replica_id2)
{
    return (replex->domain_id[replica_id1] == replex->domain_id[replica_id2]);
}

gmx_bool awhbias_replicas_can_exchange(const awh_replex_t *replex, int replica_id1, int replica_id2)
{
    /* Unless exchange between all replicas has been set, there are no exchanges between replicas with the same domain id.
       This way, the non-overlap region of each domain is ensured to have continuous trajectories. */
    return replex->all_replicas_exchange || !replicas_are_equivalent(replex, replica_id1, replica_id2);
}

/* Comparison function for integers needed by qsort. */
int int_compare(const void *a, const void *b)
{
    return (*(int *)a) - (*(int *)b);
}

/* Find the the first point with target > 0 of this replica. Assuming this uniquely defines the target range, this
   defines the domain id of the replica. The replica ids are sorted by the starting index. The domain id is then given
   by the order of the replica ids. */
static int *init_domain_id(const awh_replex_t *replex, const t_awh *awh, const gmx_multisim_t *ms)
{
    int      start              = -1;
    gmx_bool bFound_first_point = FALSE;

    for (int m = 0; m < awh->npoints && !bFound_first_point; m++)
    {
        bFound_first_point = (awh->coord_point[m].target > 0);
        start              = m;
    }

    /* Add the starting index, replica id pair to an helper array for sorting */
    int *start_replicaid;
    snew(start_replicaid, replex->nreplica*2);
    start_replicaid[replex->replica_id*2]     = start;
    start_replicaid[replex->replica_id*2 + 1] = replex->replica_id;

    /* Sort the start index, replica id pair */
    gmx_sumi_sim(replex->nreplica*2, start_replicaid, ms);
    qsort(start_replicaid, replex->nreplica, 2*sizeof(int), &int_compare);

    /* Allocate memory and set the domain id for each replica */
    int     *domain_id;
    int      domain_id_count = -1;
    int      start_prev      = -1;

    snew(domain_id, replex->nreplica);

    for (int i = 0; i < replex->nreplica*2; i += 2)
    {
        int start      = start_replicaid[i];
        int replica_id = start_replicaid[i + 1];

        if (start > start_prev)
        {
            domain_id_count++;
        }
        domain_id[replica_id] = domain_id_count;
        start_prev            = start;
    }

    sfree(start_replicaid);
    return domain_id;
}

static awh_replex_t *init_awh_replica_exchange(const t_awh *awh, const gmx_multisim_t *ms, bool all_replicas_exchange)
{
    awh_replex_t *replex;

    snew(replex, 1);

    replex->replica_id  = ms->sim;
    replex->nreplica    = ms->nsim;
    replex->npoints     = awh->npoints;

    /* Allocate memory */
    snew(replex->point_index, replex->nreplica);
    snew(replex->point_bias, replex->nreplica);
    for (int i = 0; i < replex->nreplica; i++)
    {
        snew(replex->point_bias[i], replex->npoints);
    }


    replex->all_replicas_exchange = all_replicas_exchange;
    if (all_replicas_exchange)
    {
        replex->domain_id = NULL;
    }
    else
    {
        replex->domain_id = init_domain_id(replex, awh, ms);
    }

    return replex;
}

awh_replex_t *init_awhbias_replica_exchange(const awhbias_params_t *awhbias_params, const t_awhbias *awhbias, const gmx_multisim_t *ms,
                                            const int *exchange_indices, int bMulti_exchange)
{
    awh_replex_t *replex;

    if (awhbias->nawh > 1)
    {
        gmx_fatal(FARGS, "AWH bias replica exchange only supports one AWH bias per replica, but the number of AWH biases is %d.", awhbias->nawh);
    }
    /* The target distribution can generally be of any form and be history dependent. This makes it less straightforward to define
       equivalence between AWH domains (target distributions exactly equal or > 0 in the same region?). So currently we just
       support this equivalence check for constant, non-user defined target distributions. */
    bool all_replicas_exchange = !((awhbias_params->awh_params[0].eTarget == eawhtargetCONSTANT) && (!awhbias_params->awh_params[0].bUser_data));

    replex = init_awh_replica_exchange(&awhbias->awh[0], ms, all_replicas_exchange);

    if (bMulti_exchange && !replex->all_replicas_exchange)
    {
        /* When we have multiple exchanges per exchange step (instead of nearest neighbors) we don't allow all replicas to be equivalent since
           this leads to an infinite loop when looking for non-equivalent exchange partners. */
        int nnot_equivalent = 0;

        for (int i = 1; i < replex->nreplica; i++)
        {
            if (!replicas_are_equivalent(replex, i - 1, i))
            {
                nnot_equivalent++;
            }
        }
        if (nnot_equivalent == 0)
        {
            gmx_fatal(FARGS, "Replica AWH bias exchange with multi-exchange and all AWH domains identical is not allowed. "
                      "You can try modifying your initial configurations or modifying your AWH interval mdp settings.");
        }
    }
    else
    {
        /* Nearest neighbor exchange only makes much sense in the one-dimensional coordinate case with one replica per domain. If there
         * are several replicas per domain these will not exchange, preventing mixing. The exchange indices should order the replicas
         * by increasing domain id so that there is overlap between neighboring replicas.
         */
        if (awhbias_params->awh_params[0].ndim > 1)
        {
            gmx_fatal(FARGS, "Replica AWH bias exchange with nearest neighbor exhanges is not useful for an AWH coordinate of dimension > 1 since "
                      "multidimensional indexing is currently not handled. Try using multi-exchange with the -nex option.");
        }
        for (int i = 1; i < replex->nreplica; i++)
        {
            if (replex->domain_id[exchange_indices[i]] - replex->domain_id[exchange_indices[i - 1]] != 1)
            {
                gmx_fatal(FARGS, "Replica %d has domain id %d and replica %d has domain id %d. "
                          "With nearest neighbor AWH bias replica exchange the replicas should be given in order of increasing "
                          "domain id with a domain id step of 1. You can try modifying your initial configurations or modifying your "
                          "AWH interval mdp settings. You can also use multi-exchange instead with the -nex option.",
                          exchange_indices[i - 1], replex->domain_id[exchange_indices[i - 1]],
                          exchange_indices[i], replex->domain_id[exchange_indices[i]]);
            }
        }
    }

    /* Note/TODO: obviously there is a lot of trust in the user here that the multiple simulations are compatible. Only things
       directly related to replica exchange are checked. */
    return replex;
}


static gmx_bool update_awh_replica_exchange_state(const t_awh *awh, awh_replex_t *replex, const gmx_multisim_t *ms)
{
    int      nreplicas_ok;
    gmx_bool bCoord_value_in_range, bUpdated;

    /* Only fill in the data of this replica. Zero the rest. These arrays will later be summed over all replicas. */
    for (int i = 0; i < replex->nreplica; i++)
    {
        replex->point_index[i] = 0;
        for (int j = 0; j < replex->npoints; j++)
        {
            replex->point_bias[i][j] = 0.;
        }
    }

    /* Only do exchange if the coordinate of this sim is in range of the AWH grid (it should be quite rare that it is not) */
    bCoord_value_in_range = coord_value_is_in_grid(awh);

    if (bCoord_value_in_range)
    {
        replex->point_index[replex->replica_id] = get_coord_value_index(awh);

        /* Get this replica's coordinate bias evaluated for all grid points.
         * Note: points outside of the target domain will have maximally low bias = -GMX_DOUBLE_MAX.
         * We really only need to update points close to the target domain since the others should not change
         * and for many points, looping over the whole grid is inefficient. */
        for (int i = 0; i < replex->npoints; i++)
        {
            replex->point_bias[replex->replica_id][i] = calc_convolved_bias(awh, awh->grid->point[i].value);
        }
    }

    /* Now do communication */

    /* First make sure all replicas are up for exchanging. */
    nreplicas_ok = bCoord_value_in_range ? 1 : 0;
    gmx_sumi_sim(1, &nreplicas_ok, ms);

    if (nreplicas_ok == replex->nreplica)
    {
        gmx_sumi_sim(replex->nreplica, replex->point_index, ms);
        for (int i = 0; i < replex->nreplica; i++)
        {
            gmx_sumd_sim(replex->npoints, replex->point_bias[i], ms);
        }
        bUpdated = TRUE;
    }
    else
    {
        bUpdated = FALSE;

        /* If there has been no update, replica exchange should not be attempted. */
        for (int i = 0; i < replex->nreplica; i++)
        {
            replex->point_index[i] = -1;
        }
    }

    return bUpdated;
}

gmx_bool update_awhbias_replica_exchange_state(const t_awhbias *awhbias, awh_replex_t *replex, const gmx_multisim_t *ms)
{
    /* Only a single bias can currently do replica exchange */
    return update_awh_replica_exchange_state(&awhbias->awh[0], replex, ms);
}

real calc_awhbias_replica_exchange_delta(const awh_replex_t *replex, int config_A, int config_B, int state_A, int state_B)
{
    /* Note about notation: Note that in the call to this function
       from calc_delta in repl_ex.cpp the state indices state_A, state_B are called ap, bp (p labelling permuted indices).
       This is confusing but the notation of that function assumes states are permuted (currently in mdrun it is however
       the configurations that are exchanged). The input arguments to calc_delta
       are however themselves swapped so that A, B and not ap, bp corresponds to what is actually being swapped. */
    double delta;
    int    point_A = replex->point_index[config_A], point_B = replex->point_index[config_B];

    GMX_RELEASE_ASSERT(point_A >= 0 && point_B >= 0, "Attempted to do AWH replica exchange with unset AWH coordinte points");

    /* If the bias is extremely low for one of the considered point indices we directly set delta to
     *  something large to avoid over-/underflow in the delta calculation. The bias can be this low if we are
     *  trying to evaluate the bias outside of it's domain (low probability). Most often this would happen when evaluating
     *  point B for bias A or vice versa. The bias  would not become extremely high however unless something has gone wrong
     *  in AWH making the bias values very large.
     */
    if ((replex->point_bias[state_A][point_B] < -GMX_FLOAT_MAX) || (replex->point_bias[state_A][point_A] < -GMX_FLOAT_MAX) ||
        (replex->point_bias[state_B][point_A] < -GMX_FLOAT_MAX) || (replex->point_bias[state_B][point_B] < -GMX_FLOAT_MAX))
    {
        delta = GMX_FLOAT_MAX;

        if ((replex->point_bias[state_A][point_A] < -GMX_FLOAT_MAX) || (replex->point_bias[state_B][point_B] < -GMX_FLOAT_MAX))
        {
            /* This should rarely (not ever?) happen. If it did it would mean that a replica coordinate is not following the current bias but is
               for some reason wandering out in low target probability regions. */
            char warningmsg[STRLEN];
            sprintf(warningmsg, "Between replica %d at point index %d with bias %g and replica %d at point index %d with bias %g. "
                    " A very low bias at the current point likely means that the coordinate is not following the AWH bias.",
                    state_A, point_A, replex->point_bias[state_A][point_A], state_B, point_B, replex->point_bias[state_B][point_B]);
            gmx_warning(warningmsg);
        }
    }
    else
    {
        delta = -(replex->point_bias[state_A][point_B] - replex->point_bias[state_A][point_A]) - (replex->point_bias[state_B][point_A] - replex->point_bias[state_B][point_B]);
        /* delta = -(change in bias state A  due to coordinate swap) -(change in bias state B due to coordinate swap)
         * I.e., the more positive the overall bias change, the more negative delta and the higher exchange probability.
         */
    }

    /*  The replica exchange code uses reals */
    return static_cast<real>(delta);
}
