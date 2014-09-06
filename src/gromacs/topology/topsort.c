/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#include "topsort.h"

#include <stdio.h>

#include "gromacs/legacyheaders/types/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static gmx_bool ip_pert(int ftype, const t_iparams *ip)
{
    gmx_bool bPert;
    int      i;

    if (NRFPB(ftype) == 0)
    {
        return FALSE;
    }

    switch (ftype)
    {
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
        case F_ANGLES:
        case F_G96ANGLES:
        case F_IDIHS:
            bPert = (ip->harmonic.rA  != ip->harmonic.rB ||
                     ip->harmonic.krA != ip->harmonic.krB);
            break;
        case F_MORSE:
            bPert = (ip->morse.b0A  != ip->morse.b0B ||
                     ip->morse.cbA  != ip->morse.cbB ||
                     ip->morse.betaA  != ip->morse.betaB);
            break;
        case F_RESTRBONDS:
            bPert = (ip->restraint.lowA  != ip->restraint.lowB ||
                     ip->restraint.up1A  != ip->restraint.up1B ||
                     ip->restraint.up2A  != ip->restraint.up2B ||
                     ip->restraint.kA    != ip->restraint.kB);
            break;
        case F_UREY_BRADLEY:
            bPert = (ip->u_b.thetaA  != ip->u_b.thetaB  ||
                     ip->u_b.kthetaA != ip->u_b.kthetaB ||
                     ip->u_b.r13A    != ip->u_b.r13B    ||
                     ip->u_b.kUBA    != ip->u_b.kUBB);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            bPert = (ip->pdihs.phiA != ip->pdihs.phiB ||
                     ip->pdihs.cpA  != ip->pdihs.cpB);
            break;
        case F_RBDIHS:
            bPert = FALSE;
            for (i = 0; i < NR_RBDIHS; i++)
            {
                if (ip->rbdihs.rbcA[i] != ip->rbdihs.rbcB[i])
                {
                    bPert = TRUE;
                }
            }
            break;
        case F_TABBONDS:
        case F_TABBONDSNC:
        case F_TABANGLES:
        case F_TABDIHS:
            bPert = (ip->tab.kA != ip->tab.kB);
            break;
        case F_POSRES:
            bPert = FALSE;
            for (i = 0; i < DIM; i++)
            {
                if (ip->posres.pos0A[i] != ip->posres.pos0B[i] ||
                    ip->posres.fcA[i]   != ip->posres.fcB[i])
                {
                    bPert = TRUE;
                }
            }
            break;
        case F_DIHRES:
            bPert = ((ip->dihres.phiA != ip->dihres.phiB) ||
                     (ip->dihres.dphiA != ip->dihres.dphiB) ||
                     (ip->dihres.kfacA != ip->dihres.kfacB));
            break;
        case F_LJ14:
            bPert = (ip->lj14.c6A  != ip->lj14.c6B ||
                     ip->lj14.c12A != ip->lj14.c12B);
            break;
        case F_CMAP:
            bPert = FALSE;
            break;
        case F_RESTRANGLES:
        case F_RESTRDIHS:
        case F_CBTDIHS:
            bPert = FALSE;
            gmx_fatal(FARGS, "Function type %s does not support currentely free energy calculations",
                      interaction_function[ftype].longname);
        default:
            bPert = FALSE;
            gmx_fatal(FARGS, "Function type %s not implemented in ip_pert",
                      interaction_function[ftype].longname);
    }

    return bPert;
}

static gmx_bool ip_q_pert(int ftype, const t_iatom *ia,
                          const t_iparams *ip, const real *qA, const real *qB)
{
    /* 1-4 interactions do not have the charges stored in the iparams list,
     * so we need a separate check for those.
     */
    return (ip_pert(ftype, ip+ia[0]) ||
            (ftype == F_LJ14 && (qA[ia[1]] != qB[ia[1]] ||
                                 qA[ia[2]] != qB[ia[2]])));
}

gmx_bool gmx_mtop_bondeds_free_energy(const gmx_mtop_t *mtop)
{
    const gmx_ffparams_t *ffparams;
    int                   i, ftype;
    int                   mb;
    t_atom               *atom;
    t_ilist              *il;
    t_iatom              *ia;
    gmx_bool              bPert;

    ffparams = &mtop->ffparams;

    /* Loop over all the function types and compare the A/B parameters */
    bPert = FALSE;
    for (i = 0; i < ffparams->ntypes; i++)
    {
        ftype = ffparams->functype[i];
        if (interaction_function[ftype].flags & IF_BOND)
        {
            if (ip_pert(ftype, &ffparams->iparams[i]))
            {
                bPert = TRUE;
            }
        }
    }

    /* Check perturbed charges for 1-4 interactions */
    for (mb = 0; mb < mtop->nmolblock; mb++)
    {
        atom = mtop->moltype[mtop->molblock[mb].type].atoms.atom;
        il   = &mtop->moltype[mtop->molblock[mb].type].ilist[F_LJ14];
        ia   = il->iatoms;
        for (i = 0; i < il->nr; i += 3)
        {
            if (atom[ia[i+1]].q != atom[ia[i+1]].qB ||
                atom[ia[i+2]].q != atom[ia[i+2]].qB)
            {
                bPert = TRUE;
            }
        }
    }

    return bPert;
}

void gmx_sort_ilist_fe(t_idef *idef, const real *qA, const real *qB)
{
    int        ftype, nral, i, ic, ib, a;
    t_iparams *iparams;
    t_ilist   *ilist;
    t_iatom   *iatoms;
    gmx_bool   bPert;
    t_iatom   *iabuf;
    int        iabuf_nalloc;

    if (qB == NULL)
    {
        qB = qA;
    }

    iabuf_nalloc = 0;
    iabuf        = NULL;

    iparams = idef->iparams;

    for (ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_BOND)
        {
            ilist  = &idef->il[ftype];
            iatoms = ilist->iatoms;
            nral   = NRAL(ftype);
            ic     = 0;
            ib     = 0;
            i      = 0;
            while (i < ilist->nr)
            {
                /* Check if this interaction is perturbed */
                if (ip_q_pert(ftype, iatoms+i, iparams, qA, qB))
                {
                    /* Copy to the perturbed buffer */
                    if (ib + 1 + nral > iabuf_nalloc)
                    {
                        iabuf_nalloc = over_alloc_large(ib+1+nral);
                        srenew(iabuf, iabuf_nalloc);
                    }
                    for (a = 0; a < 1+nral; a++)
                    {
                        iabuf[ib++] = iatoms[i++];
                    }
                }
                else
                {
                    /* Copy in place */
                    for (a = 0; a < 1+nral; a++)
                    {
                        iatoms[ic++] = iatoms[i++];
                    }
                }
            }
            /* Now we now the number of non-perturbed interactions */
            ilist->nr_nonperturbed = ic;

            /* Copy the buffer with perturbed interactions to the ilist */
            for (a = 0; a < ib; a++)
            {
                iatoms[ic++] = iabuf[a];
            }

            if (debug)
            {
                fprintf(debug, "%s non-pert %d pert %d\n",
                        interaction_function[ftype].longname,
                        ilist->nr_nonperturbed,
                        ilist->nr-ilist->nr_nonperturbed);
            }
        }
    }

    sfree(iabuf);

    idef->ilsort = ilsortFE_SORTED;
}
