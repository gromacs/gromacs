/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/topology/topsort.h"

#include <cstdint>
#include <cstdio>

#include <array>
#include <filesystem>
#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/atominfo.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

static gmx_bool ip_pert(int ftype, const t_iparams* ip)
{
    if (NRFPB(ftype) == 0)
    {
        return FALSE;
    }

    bool bPert = false;
    switch (ftype)
    {
        case F_BONDS:
        case F_G96BONDS:
        case F_HARMONIC:
        case F_ANGLES:
        case F_G96ANGLES:
        case F_IDIHS:
            bPert = (ip->harmonic.rA != ip->harmonic.rB || ip->harmonic.krA != ip->harmonic.krB);
            break;
        case F_MORSE:
            bPert = (ip->morse.b0A != ip->morse.b0B || ip->morse.cbA != ip->morse.cbB
                     || ip->morse.betaA != ip->morse.betaB);
            break;
        case F_RESTRBONDS:
            bPert = (ip->restraint.lowA != ip->restraint.lowB || ip->restraint.up1A != ip->restraint.up1B
                     || ip->restraint.up2A != ip->restraint.up2B
                     || ip->restraint.kA != ip->restraint.kB);
            break;
        case F_UREY_BRADLEY:
            bPert = (ip->u_b.thetaA != ip->u_b.thetaB || ip->u_b.kthetaA != ip->u_b.kthetaB
                     || ip->u_b.r13A != ip->u_b.r13B || ip->u_b.kUBA != ip->u_b.kUBB);
            break;
        case F_LINEAR_ANGLES:
            bPert = (ip->linangle.klinA != ip->linangle.klinB || ip->linangle.aA != ip->linangle.aB);
            break;
        case F_PDIHS:
        case F_PIDIHS:
        case F_ANGRES:
        case F_ANGRESZ:
            bPert = (ip->pdihs.phiA != ip->pdihs.phiB || ip->pdihs.cpA != ip->pdihs.cpB);
            break;
        /* Fourier dihedrals have been converted to Ryckaert-Bellemans by now. Treat them the same way. */
        case F_RBDIHS:
        case F_FOURDIHS:
            bPert = FALSE;
            for (int i = 0; i < NR_RBDIHS; i++)
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
        case F_TABDIHS: bPert = (ip->tab.kA != ip->tab.kB); break;
        case F_POSRES:
            bPert = FALSE;
            for (int i = 0; i < DIM; i++)
            {
                if (ip->posres.pos0A[i] != ip->posres.pos0B[i] || ip->posres.fcA[i] != ip->posres.fcB[i])
                {
                    bPert = TRUE;
                }
            }
            break;
        case F_DIHRES:
            bPert = ((ip->dihres.phiA != ip->dihres.phiB) || (ip->dihres.dphiA != ip->dihres.dphiB)
                     || (ip->dihres.kfacA != ip->dihres.kfacB));
            break;
        case F_LJ14:
            bPert = (ip->lj14.c6A != ip->lj14.c6B || ip->lj14.c12A != ip->lj14.c12B);
            break;
        case F_CMAP: bPert = FALSE; break;
        case F_RESTRANGLES:
            bPert = (ip->harmonic.rA != ip->harmonic.rB) || (ip->harmonic.krA != ip->harmonic.krB);
            break;
        case F_RESTRDIHS:
            bPert = (ip->pdihs.phiA != ip->pdihs.phiB) || (ip->pdihs.cpA != ip->pdihs.cpB);
            break;
        case F_CBTDIHS:
            bPert = false;
            for (int i = 0; i < NR_CBTDIHS && !bPert; i++)
            {
                bPert = ip->cbtdihs.cbtcA[i] != ip->cbtdihs.cbtcB[i];
            }
            break;
        default:
            gmx_fatal(FARGS,
                      "Function type %s not implemented in ip_pert",
                      interaction_function[ftype].longname);
    }

    if (bPert && (ftype == F_RESTRANGLES || ftype == F_RESTRDIHS || ftype == F_CBTDIHS))
    {
        gmx_fatal(FARGS,
                  "Function type %s does not currently support being perturbed in free energy "
                  "calculations",
                  interaction_function[ftype].longname);
    }
    return bPert;
}


//! Return whether the atom has its charge perturbed per \c atomInfo
static bool hasPerturbedCharge(int atom, gmx::ArrayRef<const int32_t> atomInfo)
{
    return bool(atomInfo[atom] & gmx::sc_atomInfo_HasPerturbedCharge);
}

gmx_bool gmx_mtop_bondeds_free_energy(const gmx_mtop_t* mtop)
{
    const gmx_ffparams_t* ffparams = &mtop->ffparams;

    /* Loop over all the function types and compare the A/B parameters */
    gmx_bool bPert = FALSE;
    for (int i = 0; i < ffparams->numTypes(); i++)
    {
        int ftype = ffparams->functype[i];
        if (interaction_function[ftype].flags & IF_BOND)
        {
            if (ip_pert(ftype, &ffparams->iparams[i]))
            {
                bPert = TRUE;
            }
        }
    }

    /* Check perturbed charges for 1-4 interactions */
    for (const gmx_molblock_t& molb : mtop->molblock)
    {
        const t_atom*            atom = mtop->moltype[molb.type].atoms.atom;
        const InteractionList&   il   = mtop->moltype[molb.type].ilist[F_LJ14];
        gmx::ArrayRef<const int> ia   = il.iatoms;
        for (int i = 0; i < il.size(); i += 3)
        {
            if (atomHasPerturbedCharge(atom[ia[i + 1]]) || atomHasPerturbedCharge(atom[ia[i + 2]]))
            {
                bPert = TRUE;
            }
        }
    }

    return bPert;
}

void gmx_sort_ilist_fe(InteractionDefinitions* idef, gmx::ArrayRef<const int32_t> atomInfo)
{
    bool havePerturbedInteractions = false;

    int      iabuf_nalloc = 0;
    t_iatom* iabuf        = nullptr;

    for (int ftype = 0; ftype < F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_BOND)
        {
            InteractionList* ilist  = &idef->il[ftype];
            int*             iatoms = ilist->iatoms.data();
            const int        nral   = NRAL(ftype);
            int              ic     = 0;
            int              ib     = 0;
            int              i      = 0;
            while (i < ilist->size())
            {
                /* Check if this interaction is perturbed */
                if (ip_pert(ftype, idef->iparams.data() + iatoms[i])
                    || (ftype == F_LJ14
                        && (hasPerturbedCharge(iatoms[i + 1], atomInfo)
                            || hasPerturbedCharge(iatoms[i + 2], atomInfo))))
                {
                    /* Copy to the perturbed buffer */
                    if (ib + 1 + nral > iabuf_nalloc)
                    {
                        iabuf_nalloc = over_alloc_large(ib + 1 + nral);
                        srenew(iabuf, iabuf_nalloc);
                    }
                    for (int a = 0; a < 1 + nral; a++)
                    {
                        iabuf[ib++] = iatoms[i++];
                    }

                    havePerturbedInteractions = true;
                }
                else
                {
                    /* Copy in place */
                    for (int a = 0; a < 1 + nral; a++)
                    {
                        iatoms[ic++] = iatoms[i++];
                    }
                }
            }
            /* Now we know the number of non-perturbed interactions */
            idef->numNonperturbedInteractions[ftype] = ic;

            /* Copy the buffer with perturbed interactions to the ilist */
            for (int a = 0; a < ib; a++)
            {
                iatoms[ic++] = iabuf[a];
            }

            if (debug)
            {
                const int numNonperturbed = idef->numNonperturbedInteractions[ftype];
                fprintf(debug,
                        "%s non-pert %d pert %d\n",
                        interaction_function[ftype].longname,
                        numNonperturbed,
                        ilist->size() - numNonperturbed);
            }
        }
    }

    sfree(iabuf);

    idef->ilsort = (havePerturbedInteractions ? ilsortFE_SORTED : ilsortNO_FE);
}
