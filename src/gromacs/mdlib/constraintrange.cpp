/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018,2019, by the GROMACS development team, led by
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
/* This file is completely threadsafe - keep it that way! */
#include "gmxpre.h"

#include "constraintrange.h"

#include <cmath>

#include <algorithm>

#include "gromacs/mdlib/constr.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

//! Recursing function to help find all adjacent constraints.
static void constr_recur(const t_blocka*                at2con,
                         const InteractionLists&        ilist,
                         gmx::ArrayRef<const t_iparams> iparams,
                         gmx_bool                       bTopB,
                         int                            at,
                         int                            depth,
                         int                            nc,
                         int*                           path,
                         real                           r0,
                         real                           r1,
                         real*                          r2max,
                         int*                           count)
{
    int      c, con, a1;
    gmx_bool bUse;
    real     len, rn0, rn1;

    (*count)++;

    gmx::ArrayRef<const int> ia1 = ilist[F_CONSTR].iatoms;
    gmx::ArrayRef<const int> ia2 = ilist[F_CONSTRNC].iatoms;

    /* Loop over all constraints connected to this atom */
    for (c = at2con->index[at]; c < at2con->index[at + 1]; c++)
    {
        con = at2con->a[c];
        /* Do not walk over already used constraints */
        bUse = TRUE;
        for (a1 = 0; a1 < depth; a1++)
        {
            if (con == path[a1])
            {
                bUse = FALSE;
            }
        }
        if (bUse)
        {
            const int* ia = constr_iatomptr(ia1, ia2, con);
            /* Flexible constraints currently have length 0, which is incorrect */
            if (!bTopB)
            {
                len = iparams[ia[0]].constr.dA;
            }
            else
            {
                len = iparams[ia[0]].constr.dB;
            }
            /* In the worst case the bond directions alternate */
            if (nc % 2 == 0)
            {
                rn0 = r0 + len;
                rn1 = r1;
            }
            else
            {
                rn0 = r0;
                rn1 = r1 + len;
            }
            /* Assume angles of 120 degrees between all bonds */
            if (rn0 * rn0 + rn1 * rn1 + rn0 * rn1 > *r2max)
            {
                *r2max = rn0 * rn0 + rn1 * rn1 + r0 * rn1;
                if (debug)
                {
                    fprintf(debug,
                            "Found longer constraint distance: r0 %5.3f r1 %5.3f rmax %5.3f\n", rn0,
                            rn1, sqrt(*r2max));
                    for (a1 = 0; a1 < depth; a1++)
                    {
                        fprintf(debug, " %d %5.3f", path[a1],
                                iparams[constr_iatomptr(ia1, ia2, con)[0]].constr.dA);
                    }
                    fprintf(debug, " %d %5.3f\n", con, len);
                }
            }
            /* Limit the number of recursions to 1000*nc,
             * so a call does not take more than a second,
             * even for highly connected systems.
             */
            if (depth + 1 < nc && *count < 1000 * nc)
            {
                if (ia[1] == at)
                {
                    a1 = ia[2];
                }
                else
                {
                    a1 = ia[1];
                }
                /* Recursion */
                path[depth] = con;
                constr_recur(at2con, ilist, iparams, bTopB, a1, depth + 1, nc, path, rn0, rn1, r2max, count);
                path[depth] = -1;
            }
        }
    }
}

//! Find the interaction radius needed for constraints for this molecule type.
static real constr_r_max_moltype(const gmx_moltype_t*           molt,
                                 gmx::ArrayRef<const t_iparams> iparams,
                                 const t_inputrec*              ir)
{
    int natoms, *path, at, count;

    t_blocka at2con;
    real     r0, r1, r2maxA, r2maxB, rmax, lam0, lam1;

    if (molt->ilist[F_CONSTR].size() == 0 && molt->ilist[F_CONSTRNC].size() == 0)
    {
        return 0;
    }

    natoms = molt->atoms.nr;

    at2con = make_at2con(*molt, iparams, flexibleConstraintTreatment(EI_DYNAMICS(ir->eI)));
    snew(path, 1 + ir->nProjOrder);
    for (at = 0; at < 1 + ir->nProjOrder; at++)
    {
        path[at] = -1;
    }

    r2maxA = 0;
    for (at = 0; at < natoms; at++)
    {
        r0 = 0;
        r1 = 0;

        count = 0;
        constr_recur(&at2con, molt->ilist, iparams, FALSE, at, 0, 1 + ir->nProjOrder, path, r0, r1,
                     &r2maxA, &count);
    }
    if (ir->efep == efepNO)
    {
        rmax = sqrt(r2maxA);
    }
    else
    {
        r2maxB = 0;
        for (at = 0; at < natoms; at++)
        {
            r0    = 0;
            r1    = 0;
            count = 0;
            constr_recur(&at2con, molt->ilist, iparams, TRUE, at, 0, 1 + ir->nProjOrder, path, r0,
                         r1, &r2maxB, &count);
        }
        lam0 = ir->fepvals->init_lambda;
        if (EI_DYNAMICS(ir->eI))
        {
            lam0 += ir->init_step * ir->fepvals->delta_lambda;
        }
        rmax = (1 - lam0) * sqrt(r2maxA) + lam0 * sqrt(r2maxB);
        if (EI_DYNAMICS(ir->eI))
        {
            lam1 = ir->fepvals->init_lambda + (ir->init_step + ir->nsteps) * ir->fepvals->delta_lambda;
            rmax = std::max(rmax, (1 - lam1) * std::sqrt(r2maxA) + lam1 * std::sqrt(r2maxB));
        }
    }

    done_blocka(&at2con);
    sfree(path);

    return rmax;
}

real constr_r_max(const MDLogger& mdlog, const gmx_mtop_t* mtop, const t_inputrec* ir)
{
    real rmax = 0;
    for (const gmx_moltype_t& molt : mtop->moltype)
    {
        rmax = std::max(rmax, constr_r_max_moltype(&molt, mtop->ffparams.iparams, ir));
    }

    GMX_LOG(mdlog.info)
            .appendTextFormatted(
                    "Maximum distance for %d constraints, at 120 deg. angles, all-trans: %.3f nm",
                    1 + ir->nProjOrder, rmax);

    return rmax;
}

} // namespace gmx
