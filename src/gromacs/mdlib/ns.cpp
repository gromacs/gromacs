/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2017, by the GROMACS development team, led by
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

#include "ns.h"

#include <stdlib.h>
#include <string.h>

#include <cmath>

#include <algorithm>

#include "gromacs/domdec/domdec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxlib/nonbonded/nonbonded.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/nsgrid.h"
#include "gromacs/mdlib/qmmm.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/group.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/ishift.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

/*
 *    E X C L U S I O N   H A N D L I N G
 */

#ifdef DEBUG
static void SETEXCL_(t_excl e[], int i, int j)
{
    e[j] = e[j] | (1<<i);
}
static void RMEXCL_(t_excl e[], int i, int j)
{
    e[j] = e[j] & ~(1<<i);
}
static gmx_bool ISEXCL_(t_excl e[], int i, int j)
{
    return (gmx_bool)(e[j] & (1<<i));
}
static gmx_bool NOTEXCL_(t_excl e[], int i, int j)
{
    return !(ISEXCL(e, i, j));
}
#else
#define SETEXCL(e, i, j) (e)[((int) (j))] |= (1<<((int) (i)))
#define RMEXCL(e, i, j)  (e)[((int) (j))] &= (~(1<<((int) (i))))
#define ISEXCL(e, i, j)  (gmx_bool) ((e)[((int) (j))] & (1<<((int) (i))))
#define NOTEXCL(e, i, j) !(ISEXCL(e, i, j))
#endif

static int
round_up_to_simd_width(int length, int simd_width)
{
    int offset;

    offset = (simd_width > 0) ? length % simd_width : 0;

    return (offset == 0) ? length : length-offset+simd_width;
}
/************************************************
 *
 *  U T I L I T I E S    F O R    N S
 *
 ************************************************/

void reallocate_nblist(t_nblist *nl)
{
    if (gmx_debug_at)
    {
        fprintf(debug, "reallocating neigborlist (ielec=%d, ivdw=%d, igeometry=%d, type=%d), maxnri=%d\n",
                nl->ielec, nl->ivdw, nl->igeometry, nl->type, nl->maxnri);
    }
    srenew(nl->iinr,   nl->maxnri);
    if (nl->igeometry == GMX_NBLIST_GEOMETRY_CG_CG)
    {
        srenew(nl->iinr_end, nl->maxnri);
    }
    srenew(nl->gid,    nl->maxnri);
    srenew(nl->shift,  nl->maxnri);
    srenew(nl->jindex, nl->maxnri+1);
}


static void init_nblist(FILE *log, t_nblist *nl_sr,
                        int maxsr,
                        int ivdw, int ivdwmod,
                        int ielec, int ielecmod,
                        int igeometry, int type,
                        gmx_bool bElecAndVdwSwitchDiffers)
{
    t_nblist *nl;
    int       homenr;

    {
        nl     = nl_sr;
        homenr = maxsr;

        if (nl == nullptr)
        {
            return;
        }


        /* Set coul/vdw in neighborlist, and for the normal loops we determine
         * an index of which one to call.
         */
        nl->ivdw        = ivdw;
        nl->ivdwmod     = ivdwmod;
        nl->ielec       = ielec;
        nl->ielecmod    = ielecmod;
        nl->type        = type;
        nl->igeometry   = igeometry;

        if (nl->type == GMX_NBLIST_INTERACTION_FREE_ENERGY)
        {
            nl->igeometry  = GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE;
        }

        /* This will also set the simd_padding_width field */
        gmx_nonbonded_set_kernel_pointers(log, nl, bElecAndVdwSwitchDiffers);

        /* maxnri is influenced by the number of shifts (maximum is 8)
         * and the number of energy groups.
         * If it is not enough, nl memory will be reallocated during the run.
         * 4 seems to be a reasonable factor, which only causes reallocation
         * during runs with tiny and many energygroups.
         */
        nl->maxnri      = homenr*4;
        nl->maxnrj      = 0;
        nl->nri         = -1;
        nl->nrj         = 0;
        nl->iinr        = nullptr;
        nl->gid         = nullptr;
        nl->shift       = nullptr;
        nl->jindex      = nullptr;
        nl->jjnr        = nullptr;
        nl->excl_fep    = nullptr;
        reallocate_nblist(nl);
        nl->jindex[0] = 0;

        if (debug)
        {
            fprintf(debug, "Initiating neighbourlist (ielec=%d, ivdw=%d, type=%d) for %s interactions,\nwith %d SR atoms.\n",
                    nl->ielec, nl->ivdw, nl->type, gmx_nblist_geometry_names[nl->igeometry], maxsr);
        }
    }
}

void init_neighbor_list(FILE *log, t_forcerec *fr, int homenr)
{
    int        maxsr, maxsr_wat;
    int        ielec, ivdw, ielecmod, ivdwmod, type;
    int        igeometry_def, igeometry_w, igeometry_ww;
    int        i;
    gmx_bool   bElecAndVdwSwitchDiffers;
    t_nblists *nbl;

    /* maxsr     = homenr-fr->nWatMol*3; */
    maxsr     = homenr;

    if (maxsr < 0)
    {
        gmx_fatal(FARGS, "%s, %d: Negative number of short range atoms.\n"
                  "Call your GROMACS dealer for assistance.", __FILE__, __LINE__);
    }
    /* This is just for initial allocation, so we do not reallocate
     * all the nlist arrays many times in a row.
     * The numbers seem very accurate, but they are uncritical.
     */
    maxsr_wat = std::min(fr->nWatMol, (homenr+2)/3);

    /* Determine the values for ielec/ivdw. */
    ielec                    = fr->nbkernel_elec_interaction;
    ivdw                     = fr->nbkernel_vdw_interaction;
    ielecmod                 = fr->nbkernel_elec_modifier;
    ivdwmod                  = fr->nbkernel_vdw_modifier;
    type                     = GMX_NBLIST_INTERACTION_STANDARD;
    bElecAndVdwSwitchDiffers = ( (fr->ic->rcoulomb_switch != fr->ic->rvdw_switch) || (fr->ic->rcoulomb != fr->ic->rvdw));

    fr->ns->bCGlist = (getenv("GMX_NBLISTCG") != nullptr);
    if (!fr->ns->bCGlist)
    {
        igeometry_def = GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE;
    }
    else
    {
        igeometry_def = GMX_NBLIST_GEOMETRY_CG_CG;
        if (log != nullptr)
        {
            fprintf(log, "\nUsing charge-group - charge-group neighbor lists and kernels\n\n");
        }
    }

    if (fr->solvent_opt == esolTIP4P)
    {
        igeometry_w  = GMX_NBLIST_GEOMETRY_WATER4_PARTICLE;
        igeometry_ww = GMX_NBLIST_GEOMETRY_WATER4_WATER4;
    }
    else
    {
        igeometry_w  = GMX_NBLIST_GEOMETRY_WATER3_PARTICLE;
        igeometry_ww = GMX_NBLIST_GEOMETRY_WATER3_WATER3;
    }

    for (i = 0; i < fr->nnblists; i++)
    {
        nbl = &(fr->nblists[i]);

        init_nblist(log, &nbl->nlist_sr[eNL_VDWQQ],
                    maxsr, ivdw, ivdwmod, ielec, ielecmod, igeometry_def, type, bElecAndVdwSwitchDiffers);
        init_nblist(log, &nbl->nlist_sr[eNL_VDW],
                    maxsr, ivdw, ivdwmod, GMX_NBKERNEL_ELEC_NONE, eintmodNONE, igeometry_def, type, bElecAndVdwSwitchDiffers);
        init_nblist(log, &nbl->nlist_sr[eNL_QQ],
                    maxsr, GMX_NBKERNEL_VDW_NONE, eintmodNONE, ielec, ielecmod, igeometry_def, type, bElecAndVdwSwitchDiffers);
        init_nblist(log, &nbl->nlist_sr[eNL_VDWQQ_WATER],
                    maxsr_wat, ivdw, ivdwmod, ielec, ielecmod, igeometry_w, type, bElecAndVdwSwitchDiffers);
        init_nblist(log, &nbl->nlist_sr[eNL_QQ_WATER],
                    maxsr_wat, GMX_NBKERNEL_VDW_NONE, eintmodNONE, ielec, ielecmod, igeometry_w, type, bElecAndVdwSwitchDiffers);
        init_nblist(log, &nbl->nlist_sr[eNL_VDWQQ_WATERWATER],
                    maxsr_wat, ivdw, ivdwmod, ielec, ielecmod, igeometry_ww, type, bElecAndVdwSwitchDiffers);
        init_nblist(log, &nbl->nlist_sr[eNL_QQ_WATERWATER],
                    maxsr_wat, GMX_NBKERNEL_VDW_NONE, eintmodNONE, ielec, ielecmod, igeometry_ww, type, bElecAndVdwSwitchDiffers);

        /* Did we get the solvent loops so we can use optimized water kernels? */
        if (nbl->nlist_sr[eNL_VDWQQ_WATER].kernelptr_vf == nullptr
            || nbl->nlist_sr[eNL_QQ_WATER].kernelptr_vf == nullptr
            || nbl->nlist_sr[eNL_VDWQQ_WATERWATER].kernelptr_vf == nullptr
            || nbl->nlist_sr[eNL_QQ_WATERWATER].kernelptr_vf == nullptr)
        {
            fr->solvent_opt = esolNO;
            if (log != nullptr)
            {
                fprintf(log, "Note: The available nonbonded kernels do not support water optimization - disabling.\n");
            }
        }

        if (fr->efep != efepNO)
        {
            init_nblist(log, &nbl->nlist_sr[eNL_VDWQQ_FREE],
                        maxsr, ivdw, ivdwmod, ielec, ielecmod, GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE, GMX_NBLIST_INTERACTION_FREE_ENERGY, bElecAndVdwSwitchDiffers);
            init_nblist(log, &nbl->nlist_sr[eNL_VDW_FREE],
                        maxsr, ivdw, ivdwmod, GMX_NBKERNEL_ELEC_NONE, eintmodNONE, GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE, GMX_NBLIST_INTERACTION_FREE_ENERGY, bElecAndVdwSwitchDiffers);
            init_nblist(log, &nbl->nlist_sr[eNL_QQ_FREE],
                        maxsr, GMX_NBKERNEL_VDW_NONE, eintmodNONE, ielec, ielecmod, GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE, GMX_NBLIST_INTERACTION_FREE_ENERGY, bElecAndVdwSwitchDiffers);
        }
    }
    /* QMMM MM list */
    if (fr->bQMMM && fr->qr->QMMMscheme != eQMMMschemeoniom)
    {
        if (nullptr == fr->QMMMlist)
        {
            snew(fr->QMMMlist, 1);
        }
        init_nblist(log, fr->QMMMlist,
                    maxsr, 0, 0, ielec, ielecmod, GMX_NBLIST_GEOMETRY_PARTICLE_PARTICLE, GMX_NBLIST_INTERACTION_STANDARD, bElecAndVdwSwitchDiffers);
    }

    if (log != nullptr)
    {
        fprintf(log, "\n");
    }

    fr->ns->nblist_initialized = TRUE;
}

static void reset_nblist(t_nblist *nl)
{
    nl->nri       = -1;
    nl->nrj       = 0;
    if (nl->jindex)
    {
        nl->jindex[0] = 0;
    }
}

static void reset_neighbor_lists(t_forcerec *fr)
{
    int n, i;

    if (fr->bQMMM)
    {
        /* only reset the short-range nblist */
        reset_nblist(fr->QMMMlist);
    }

    for (n = 0; n < fr->nnblists; n++)
    {
        for (i = 0; i < eNL_NR; i++)
        {
            reset_nblist( &(fr->nblists[n].nlist_sr[i]) );
        }
    }
}




static gmx_inline void new_i_nblist(t_nblist *nlist, int i_atom, int shift, int gid)
{
    int    nri = nlist->nri;

    /* Check whether we have to increase the i counter */
    if ((nri == -1) ||
        (nlist->iinr[nri]  != i_atom) ||
        (nlist->shift[nri] != shift) ||
        (nlist->gid[nri]   != gid))
    {
        /* This is something else. Now see if any entries have
         * been added in the list of the previous atom.
         */
        if ((nri == -1) ||
            ((nlist->jindex[nri+1] > nlist->jindex[nri]) &&
             (nlist->gid[nri] != -1)))
        {
            /* If so increase the counter */
            nlist->nri++;
            nri++;
            if (nlist->nri >= nlist->maxnri)
            {
                nlist->maxnri += over_alloc_large(nlist->nri);
                reallocate_nblist(nlist);
            }
        }
        /* Set the number of neighbours and the atom number */
        nlist->jindex[nri+1] = nlist->jindex[nri];
        nlist->iinr[nri]     = i_atom;
        nlist->gid[nri]      = gid;
        nlist->shift[nri]    = shift;
    }
    else
    {
        /* Adding to previous list. First remove possible previous padding */
        if (nlist->simd_padding_width > 1)
        {
            while (nlist->nrj > 0 && nlist->jjnr[nlist->nrj-1] < 0)
            {
                nlist->nrj--;
            }
        }
    }
}

static gmx_inline void close_i_nblist(t_nblist *nlist)
{
    int nri = nlist->nri;
    int len;

    if (nri >= 0)
    {
        /* Add elements up to padding. Since we allocate memory in units
         * of the simd_padding width, we do not have to check for possible
         * list reallocation here.
         */
        while ((nlist->nrj % nlist->simd_padding_width) != 0)
        {
            /* Use -4 here, so we can write forces for 4 atoms before real data */
            nlist->jjnr[nlist->nrj++] = -4;
        }
        nlist->jindex[nri+1] = nlist->nrj;

        len = nlist->nrj -  nlist->jindex[nri];
        /* If there are no j-particles we have to reduce the
         * number of i-particles again, to prevent errors in the
         * kernel functions.
         */
        if ((len == 0) && (nlist->nri > 0))
        {
            nlist->nri--;
        }
    }
}

static gmx_inline void close_nblist(t_nblist *nlist)
{
    /* Only close this nblist when it has been initialized.
     * Avoid the creation of i-lists with no j-particles.
     */
    if (nlist->nrj == 0)
    {
        /* Some assembly kernels do not support empty lists,
         * make sure here that we don't generate any empty lists.
         * With the current ns code this branch is taken in two cases:
         * No i-particles at all: nri=-1 here
         * There are i-particles, but no j-particles; nri=0 here
         */
        nlist->nri = 0;
    }
    else
    {
        /* Close list number nri by incrementing the count */
        nlist->nri++;
    }
}

static gmx_inline void close_neighbor_lists(t_forcerec *fr, gmx_bool bMakeQMMMnblist)
{
    int n, i;

    if (bMakeQMMMnblist)
    {
        close_nblist(fr->QMMMlist);
    }

    for (n = 0; n < fr->nnblists; n++)
    {
        for (i = 0; (i < eNL_NR); i++)
        {
            close_nblist(&(fr->nblists[n].nlist_sr[i]));
        }
    }
}


static gmx_inline void add_j_to_nblist(t_nblist *nlist, int j_atom)
{
    int nrj = nlist->nrj;

    if (nlist->nrj >= nlist->maxnrj)
    {
        nlist->maxnrj = round_up_to_simd_width(over_alloc_small(nlist->nrj + 1), nlist->simd_padding_width);

        if (gmx_debug_at)
        {
            fprintf(debug, "Increasing SR nblist (ielec=%d,ivdw=%d,type=%d,igeometry=%d) j size to %d\n",
                    nlist->ielec, nlist->ivdw, nlist->type, nlist->igeometry, nlist->maxnrj);
        }

        srenew(nlist->jjnr, nlist->maxnrj);
    }

    nlist->jjnr[nrj] = j_atom;
    nlist->nrj++;
}

static gmx_inline void add_j_to_nblist_cg(t_nblist *nlist,
                                          int j_start, int j_end,
                                          t_excl *bexcl, gmx_bool i_is_j)
{
    int nrj = nlist->nrj;
    int j;

    if (nlist->nrj >= nlist->maxnrj)
    {
        nlist->maxnrj = over_alloc_small(nlist->nrj + 1);
        if (gmx_debug_at)
        {
            fprintf(debug, "Increasing SR nblist (ielec=%d,ivdw=%d,type=%d,igeometry=%d) j size to %d\n",
                    nlist->ielec, nlist->ivdw, nlist->type, nlist->igeometry, nlist->maxnrj);
        }

        srenew(nlist->jjnr, nlist->maxnrj);
        srenew(nlist->jjnr_end, nlist->maxnrj);
        srenew(nlist->excl, nlist->maxnrj*MAX_CGCGSIZE);
    }

    nlist->jjnr[nrj]     = j_start;
    nlist->jjnr_end[nrj] = j_end;

    if (j_end - j_start > MAX_CGCGSIZE)
    {
        gmx_fatal(FARGS, "The charge-group - charge-group neighborlist do not support charge groups larger than %d, found a charge group of size %d", MAX_CGCGSIZE, j_end-j_start);
    }

    /* Set the exclusions */
    for (j = j_start; j < j_end; j++)
    {
        nlist->excl[nrj*MAX_CGCGSIZE + j - j_start] = bexcl[j];
    }
    if (i_is_j)
    {
        /* Avoid double counting of intra-cg interactions */
        for (j = 1; j < j_end-j_start; j++)
        {
            nlist->excl[nrj*MAX_CGCGSIZE + j] |= (1<<j) - 1;
        }
    }

    nlist->nrj++;
}

typedef void
    put_in_list_t (gmx_bool              bHaveVdW[],
                   int                   ngid,
                   t_mdatoms     *       md,
                   int                   icg,
                   int                   jgid,
                   int                   nj,
                   int                   jjcg[],
                   int                   index[],
                   t_excl                bExcl[],
                   int                   shift,
                   t_forcerec     *      fr,
                   gmx_bool              bDoVdW,
                   gmx_bool              bDoCoul,
                   int                   solvent_opt);

static void
put_in_list_at(gmx_bool              bHaveVdW[],
               int                   ngid,
               t_mdatoms     *       md,
               int                   icg,
               int                   jgid,
               int                   nj,
               int                   jjcg[],
               int                   index[],
               t_excl                bExcl[],
               int                   shift,
               t_forcerec     *      fr,
               gmx_bool              bDoVdW,
               gmx_bool              bDoCoul,
               int                   solvent_opt)
{
    /* The a[] index has been removed,
     * to put it back in i_atom should be a[i0] and jj should be a[jj].
     */
    t_nblist  *   vdwc;
    t_nblist  *   vdw;
    t_nblist  *   coul;
    t_nblist  *   vdwc_free  = nullptr;
    t_nblist  *   vdw_free   = nullptr;
    t_nblist  *   coul_free  = nullptr;
    t_nblist  *   vdwc_ww    = nullptr;
    t_nblist  *   coul_ww    = nullptr;

    int           i, j, jcg, igid, gid, nbl_ind;
    int           jj, jj0, jj1, i_atom;
    int           i0, nicg;

    int          *cginfo;
    int          *type, *typeB;
    real         *charge, *chargeB;
    real          qi, qiB;
    gmx_bool      bFreeEnergy, bFree, bFreeJ, bNotEx, *bPert;
    gmx_bool      bDoVdW_i, bDoCoul_i, bDoCoul_i_sol;
    int           iwater, jwater;
    t_nblist     *nlist;

    /* Copy some pointers */
    cginfo  = fr->cginfo;
    charge  = md->chargeA;
    chargeB = md->chargeB;
    type    = md->typeA;
    typeB   = md->typeB;
    bPert   = md->bPerturbed;

    /* Get atom range */
    i0     = index[icg];
    nicg   = index[icg+1]-i0;

    /* Get the i charge group info */
    igid   = GET_CGINFO_GID(cginfo[icg]);

    iwater = (solvent_opt != esolNO) ? GET_CGINFO_SOLOPT(cginfo[icg]) : esolNO;

    bFreeEnergy = FALSE;
    if (md->nPerturbed)
    {
        /* Check if any of the particles involved are perturbed.
         * If not we can do the cheaper normal put_in_list
         * and use more solvent optimization.
         */
        for (i = 0; i < nicg; i++)
        {
            bFreeEnergy |= bPert[i0+i];
        }
        /* Loop over the j charge groups */
        for (j = 0; (j < nj && !bFreeEnergy); j++)
        {
            jcg = jjcg[j];
            jj0 = index[jcg];
            jj1 = index[jcg+1];
            /* Finally loop over the atoms in the j-charge group */
            for (jj = jj0; jj < jj1; jj++)
            {
                bFreeEnergy |= bPert[jj];
            }
        }
    }

    /* Unpack pointers to neighbourlist structs */
    if (fr->nnblists == 1)
    {
        nbl_ind = 0;
    }
    else
    {
        nbl_ind = fr->gid2nblists[GID(igid, jgid, ngid)];
    }
    nlist = fr->nblists[nbl_ind].nlist_sr;

    if (iwater != esolNO)
    {
        vdwc    = &nlist[eNL_VDWQQ_WATER];
        vdw     = &nlist[eNL_VDW];
        coul    = &nlist[eNL_QQ_WATER];
        vdwc_ww = &nlist[eNL_VDWQQ_WATERWATER];
        coul_ww = &nlist[eNL_QQ_WATERWATER];
    }
    else
    {
        vdwc = &nlist[eNL_VDWQQ];
        vdw  = &nlist[eNL_VDW];
        coul = &nlist[eNL_QQ];
    }

    if (!bFreeEnergy)
    {
        if (iwater != esolNO)
        {
            /* Loop over the atoms in the i charge group */
            i_atom  = i0;
            gid     = GID(igid, jgid, ngid);
            /* Create new i_atom for each energy group */
            if (bDoCoul && bDoVdW)
            {
                new_i_nblist(vdwc, i_atom, shift, gid);
                new_i_nblist(vdwc_ww, i_atom, shift, gid);
            }
            if (bDoVdW)
            {
                new_i_nblist(vdw, i_atom, shift, gid);
            }
            if (bDoCoul)
            {
                new_i_nblist(coul, i_atom, shift, gid);
                new_i_nblist(coul_ww, i_atom, shift, gid);
            }
            /* Loop over the j charge groups */
            for (j = 0; (j < nj); j++)
            {
                jcg = jjcg[j];

                if (jcg == icg)
                {
                    continue;
                }

                jj0    = index[jcg];
                jwater = GET_CGINFO_SOLOPT(cginfo[jcg]);

                if (iwater == esolSPC && jwater == esolSPC)
                {
                    /* Interaction between two SPC molecules */
                    if (!bDoCoul)
                    {
                        /* VdW only - only first atoms in each water interact */
                        add_j_to_nblist(vdw, jj0);
                    }
                    else
                    {
                        /* One entry for the entire water-water interaction */
                        if (!bDoVdW)
                        {
                            add_j_to_nblist(coul_ww, jj0);
                        }
                        else
                        {
                            add_j_to_nblist(vdwc_ww, jj0);
                        }
                    }
                }
                else if (iwater == esolTIP4P && jwater == esolTIP4P)
                {
                    /* Interaction between two TIP4p molecules */
                    if (!bDoCoul)
                    {
                        /* VdW only - only first atoms in each water interact */
                        add_j_to_nblist(vdw, jj0);
                    }
                    else
                    {
                        /* One entry for the entire water-water interaction */
                        if (!bDoVdW)
                        {
                            add_j_to_nblist(coul_ww, jj0);
                        }
                        else
                        {
                            add_j_to_nblist(vdwc_ww, jj0);
                        }
                    }
                }
                else
                {
                    /* j charge group is not water, but i is.
                     * Add entries to the water-other_atom lists; the geometry of the water
                     * molecule doesn't matter - that is taken care of in the nonbonded kernel,
                     * so we don't care if it is SPC or TIP4P...
                     */

                    jj1 = index[jcg+1];

                    if (!bDoVdW)
                    {
                        for (jj = jj0; (jj < jj1); jj++)
                        {
                            if (charge[jj] != 0)
                            {
                                add_j_to_nblist(coul, jj);
                            }
                        }
                    }
                    else if (!bDoCoul)
                    {
                        for (jj = jj0; (jj < jj1); jj++)
                        {
                            if (bHaveVdW[type[jj]])
                            {
                                add_j_to_nblist(vdw, jj);
                            }
                        }
                    }
                    else
                    {
                        /* _charge_ _groups_ interact with both coulomb and LJ */
                        /* Check which atoms we should add to the lists!       */
                        for (jj = jj0; (jj < jj1); jj++)
                        {
                            if (bHaveVdW[type[jj]])
                            {
                                if (charge[jj] != 0)
                                {
                                    add_j_to_nblist(vdwc, jj);
                                }
                                else
                                {
                                    add_j_to_nblist(vdw, jj);
                                }
                            }
                            else if (charge[jj] != 0)
                            {
                                add_j_to_nblist(coul, jj);
                            }
                        }
                    }
                }
            }
            close_i_nblist(vdw);
            close_i_nblist(coul);
            close_i_nblist(vdwc);
            close_i_nblist(coul_ww);
            close_i_nblist(vdwc_ww);
        }
        else
        {
            /* no solvent as i charge group */
            /* Loop over the atoms in the i charge group */
            for (i = 0; i < nicg; i++)
            {
                i_atom  = i0+i;
                gid     = GID(igid, jgid, ngid);
                qi      = charge[i_atom];

                /* Create new i_atom for each energy group */
                if (bDoVdW && bDoCoul)
                {
                    new_i_nblist(vdwc, i_atom, shift, gid);
                }
                if (bDoVdW)
                {
                    new_i_nblist(vdw, i_atom, shift, gid);
                }
                if (bDoCoul)
                {
                    new_i_nblist(coul, i_atom, shift, gid);
                }
                bDoVdW_i  = (bDoVdW  && bHaveVdW[type[i_atom]]);
                bDoCoul_i = (bDoCoul && qi != 0);

                if (bDoVdW_i || bDoCoul_i)
                {
                    /* Loop over the j charge groups */
                    for (j = 0; (j < nj); j++)
                    {
                        jcg = jjcg[j];

                        /* Check for large charge groups */
                        if (jcg == icg)
                        {
                            jj0 = i0 + i + 1;
                        }
                        else
                        {
                            jj0 = index[jcg];
                        }

                        jj1 = index[jcg+1];
                        /* Finally loop over the atoms in the j-charge group */
                        for (jj = jj0; jj < jj1; jj++)
                        {
                            bNotEx = NOTEXCL(bExcl, i, jj);

                            if (bNotEx)
                            {
                                if (!bDoVdW_i)
                                {
                                    if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul, jj);
                                    }
                                }
                                else if (!bDoCoul_i)
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        add_j_to_nblist(vdw, jj);
                                    }
                                }
                                else
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        if (charge[jj] != 0)
                                        {
                                            add_j_to_nblist(vdwc, jj);
                                        }
                                        else
                                        {
                                            add_j_to_nblist(vdw, jj);
                                        }
                                    }
                                    else if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul, jj);
                                    }
                                }
                            }
                        }
                    }
                }
                close_i_nblist(vdw);
                close_i_nblist(coul);
                close_i_nblist(vdwc);
            }
        }
    }
    else
    {
        /* we are doing free energy */
        vdwc_free = &nlist[eNL_VDWQQ_FREE];
        vdw_free  = &nlist[eNL_VDW_FREE];
        coul_free = &nlist[eNL_QQ_FREE];
        /* Loop over the atoms in the i charge group */
        for (i = 0; i < nicg; i++)
        {
            i_atom  = i0+i;
            gid     = GID(igid, jgid, ngid);
            qi      = charge[i_atom];
            qiB     = chargeB[i_atom];

            /* Create new i_atom for each energy group */
            if (bDoVdW && bDoCoul)
            {
                new_i_nblist(vdwc, i_atom, shift, gid);
            }
            if (bDoVdW)
            {
                new_i_nblist(vdw, i_atom, shift, gid);
            }
            if (bDoCoul)
            {
                new_i_nblist(coul, i_atom, shift, gid);
            }

            new_i_nblist(vdw_free, i_atom, shift, gid);
            new_i_nblist(coul_free, i_atom, shift, gid);
            new_i_nblist(vdwc_free, i_atom, shift, gid);

            bDoVdW_i  = (bDoVdW  &&
                         (bHaveVdW[type[i_atom]] || bHaveVdW[typeB[i_atom]]));
            bDoCoul_i = (bDoCoul && (qi != 0 || qiB != 0));
            /* For TIP4P the first atom does not have a charge,
             * but the last three do. So we should still put an atom
             * without LJ but with charge in the water-atom neighborlist
             * for a TIP4p i charge group.
             * For SPC type water the first atom has LJ and charge,
             * so there is no such problem.
             */
            if (iwater == esolNO)
            {
                bDoCoul_i_sol = bDoCoul_i;
            }
            else
            {
                bDoCoul_i_sol = bDoCoul;
            }

            if (bDoVdW_i || bDoCoul_i_sol)
            {
                /* Loop over the j charge groups */
                for (j = 0; (j < nj); j++)
                {
                    jcg = jjcg[j];

                    /* Check for large charge groups */
                    if (jcg == icg)
                    {
                        jj0 = i0 + i + 1;
                    }
                    else
                    {
                        jj0 = index[jcg];
                    }

                    jj1 = index[jcg+1];
                    /* Finally loop over the atoms in the j-charge group */
                    bFree = bPert[i_atom];
                    for (jj = jj0; (jj < jj1); jj++)
                    {
                        bFreeJ = bFree || bPert[jj];
                        /* Complicated if, because the water H's should also
                         * see perturbed j-particles
                         */
                        if (iwater == esolNO || i == 0 || bFreeJ)
                        {
                            bNotEx = NOTEXCL(bExcl, i, jj);

                            if (bNotEx)
                            {
                                if (bFreeJ)
                                {
                                    if (!bDoVdW_i)
                                    {
                                        if (charge[jj] != 0 || chargeB[jj] != 0)
                                        {
                                            add_j_to_nblist(coul_free, jj);
                                        }
                                    }
                                    else if (!bDoCoul_i)
                                    {
                                        if (bHaveVdW[type[jj]] || bHaveVdW[typeB[jj]])
                                        {
                                            add_j_to_nblist(vdw_free, jj);
                                        }
                                    }
                                    else
                                    {
                                        if (bHaveVdW[type[jj]] || bHaveVdW[typeB[jj]])
                                        {
                                            if (charge[jj] != 0 || chargeB[jj] != 0)
                                            {
                                                add_j_to_nblist(vdwc_free, jj);
                                            }
                                            else
                                            {
                                                add_j_to_nblist(vdw_free, jj);
                                            }
                                        }
                                        else if (charge[jj] != 0 || chargeB[jj] != 0)
                                        {
                                            add_j_to_nblist(coul_free, jj);
                                        }
                                    }
                                }
                                else if (!bDoVdW_i)
                                {
                                    /* This is done whether or not bWater is set */
                                    if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul, jj);
                                    }
                                }
                                else if (!bDoCoul_i_sol)
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        add_j_to_nblist(vdw, jj);
                                    }
                                }
                                else
                                {
                                    if (bHaveVdW[type[jj]])
                                    {
                                        if (charge[jj] != 0)
                                        {
                                            add_j_to_nblist(vdwc, jj);
                                        }
                                        else
                                        {
                                            add_j_to_nblist(vdw, jj);
                                        }
                                    }
                                    else if (charge[jj] != 0)
                                    {
                                        add_j_to_nblist(coul, jj);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            close_i_nblist(vdw);
            close_i_nblist(coul);
            close_i_nblist(vdwc);
            close_i_nblist(vdw_free);
            close_i_nblist(coul_free);
            close_i_nblist(vdwc_free);
        }
    }
}

static void
put_in_list_qmmm(gmx_bool gmx_unused              bHaveVdW[],
                 int                              ngid,
                 t_mdatoms gmx_unused     *       md,
                 int                              icg,
                 int                              jgid,
                 int                              nj,
                 int                              jjcg[],
                 int                              index[],
                 t_excl                           bExcl[],
                 int                              shift,
                 t_forcerec                *      fr,
                 gmx_bool  gmx_unused             bDoVdW,
                 gmx_bool  gmx_unused             bDoCoul,
                 int       gmx_unused             solvent_opt)
{
    t_nblist  *   coul;
    int           i, j, jcg, igid, gid;
    int           jj, jj0, jj1, i_atom;
    int           i0, nicg;
    gmx_bool      bNotEx;

    /* Get atom range */
    i0     = index[icg];
    nicg   = index[icg+1]-i0;

    /* Get the i charge group info */
    igid   = GET_CGINFO_GID(fr->cginfo[icg]);

    coul = fr->QMMMlist;

    /* Loop over atoms in the ith charge group */
    for (i = 0; i < nicg; i++)
    {
        i_atom = i0+i;
        gid    = GID(igid, jgid, ngid);
        /* Create new i_atom for each energy group */
        new_i_nblist(coul, i_atom, shift, gid);

        /* Loop over the j charge groups */
        for (j = 0; j < nj; j++)
        {
            jcg = jjcg[j];

            /* Charge groups cannot have QM and MM atoms simultaneously */
            if (jcg != icg)
            {
                jj0 = index[jcg];
                jj1 = index[jcg+1];
                /* Finally loop over the atoms in the j-charge group */
                for (jj = jj0; jj < jj1; jj++)
                {
                    bNotEx = NOTEXCL(bExcl, i, jj);
                    if (bNotEx)
                    {
                        add_j_to_nblist(coul, jj);
                    }
                }
            }
        }
        close_i_nblist(coul);
    }
}

static void
put_in_list_cg(gmx_bool  gmx_unused             bHaveVdW[],
               int                              ngid,
               t_mdatoms  gmx_unused    *       md,
               int                              icg,
               int                              jgid,
               int                              nj,
               int                              jjcg[],
               int                              index[],
               t_excl                           bExcl[],
               int                              shift,
               t_forcerec                *      fr,
               gmx_bool   gmx_unused            bDoVdW,
               gmx_bool   gmx_unused            bDoCoul,
               int        gmx_unused            solvent_opt)
{
    int          cginfo;
    int          igid, gid, nbl_ind;
    t_nblist *   vdwc;
    int          j, jcg;

    cginfo = fr->cginfo[icg];

    igid = GET_CGINFO_GID(cginfo);
    gid  = GID(igid, jgid, ngid);

    /* Unpack pointers to neighbourlist structs */
    if (fr->nnblists == 1)
    {
        nbl_ind = 0;
    }
    else
    {
        nbl_ind = fr->gid2nblists[gid];
    }
    vdwc = &fr->nblists[nbl_ind].nlist_sr[eNL_VDWQQ];

    /* Make a new neighbor list for charge group icg.
     * Currently simply one neighbor list is made with LJ and Coulomb.
     * If required, zero interactions could be removed here
     * or in the force loop.
     */
    new_i_nblist(vdwc, index[icg], shift, gid);
    vdwc->iinr_end[vdwc->nri] = index[icg+1];

    for (j = 0; (j < nj); j++)
    {
        jcg = jjcg[j];
        /* Skip the icg-icg pairs if all self interactions are excluded */
        if (!(jcg == icg && GET_CGINFO_EXCL_INTRA(cginfo)))
        {
            /* Here we add the j charge group jcg to the list,
             * exclusions are also added to the list.
             */
            add_j_to_nblist_cg(vdwc, index[jcg], index[jcg+1], bExcl, icg == jcg);
        }
    }

    close_i_nblist(vdwc);
}

static void setexcl(int start, int end, t_blocka *excl, gmx_bool b,
                    t_excl bexcl[])
{
    int i, k;

    if (b)
    {
        for (i = start; i < end; i++)
        {
            for (k = excl->index[i]; k < excl->index[i+1]; k++)
            {
                SETEXCL(bexcl, i-start, excl->a[k]);
            }
        }
    }
    else
    {
        for (i = start; i < end; i++)
        {
            for (k = excl->index[i]; k < excl->index[i+1]; k++)
            {
                RMEXCL(bexcl, i-start, excl->a[k]);
            }
        }
    }
}

int calc_naaj(int icg, int cgtot)
{
    int naaj;

    if ((cgtot % 2) == 1)
    {
        /* Odd number of charge groups, easy */
        naaj = 1 + (cgtot/2);
    }
    else if ((cgtot % 4) == 0)
    {
        /* Multiple of four is hard */
        if (icg < cgtot/2)
        {
            if ((icg % 2) == 0)
            {
                naaj = 1+(cgtot/2);
            }
            else
            {
                naaj = cgtot/2;
            }
        }
        else
        {
            if ((icg % 2) == 1)
            {
                naaj = 1+(cgtot/2);
            }
            else
            {
                naaj = cgtot/2;
            }
        }
    }
    else
    {
        /* cgtot/2 = odd */
        if ((icg % 2) == 0)
        {
            naaj = 1+(cgtot/2);
        }
        else
        {
            naaj = cgtot/2;
        }
    }
#ifdef DEBUG
    fprintf(log, "naaj=%d\n", naaj);
#endif

    return naaj;
}

/************************************************
 *
 *  S I M P L E      C O R E     S T U F F
 *
 ************************************************/

static real calc_image_tric(rvec xi, rvec xj, matrix box,
                            rvec b_inv, int *shift)
{
    /* This code assumes that the cut-off is smaller than
     * a half times the smallest diagonal element of the box.
     */
    const real h25 = 2.5;
    real       dx, dy, dz;
    real       r2;
    int        tx, ty, tz;

    /* Compute diff vector */
    dz = xj[ZZ] - xi[ZZ];
    dy = xj[YY] - xi[YY];
    dx = xj[XX] - xi[XX];

    /* Perform NINT operation, using trunc operation, therefore
     * we first add 2.5 then subtract 2 again
     */
    tz  = static_cast<int>(dz*b_inv[ZZ] + h25);
    tz -= 2;
    dz -= tz*box[ZZ][ZZ];
    dy -= tz*box[ZZ][YY];
    dx -= tz*box[ZZ][XX];

    ty  = static_cast<int>(dy*b_inv[YY] + h25);
    ty -= 2;
    dy -= ty*box[YY][YY];
    dx -= ty*box[YY][XX];

    tx  = static_cast<int>(dx*b_inv[XX]+h25);
    tx -= 2;
    dx -= tx*box[XX][XX];

    /* Distance squared */
    r2 = (dx*dx) + (dy*dy) + (dz*dz);

    *shift = XYZ2IS(tx, ty, tz);

    return r2;
}

static real calc_image_rect(rvec xi, rvec xj, rvec box_size,
                            rvec b_inv, int *shift)
{
    const real h15 = 1.5;
    real       ddx, ddy, ddz;
    real       dx, dy, dz;
    real       r2;
    int        tx, ty, tz;

    /* Compute diff vector */
    dx = xj[XX] - xi[XX];
    dy = xj[YY] - xi[YY];
    dz = xj[ZZ] - xi[ZZ];

    /* Perform NINT operation, using trunc operation, therefore
     * we first add 1.5 then subtract 1 again
     */
    tx = static_cast<int>(dx*b_inv[XX] + h15);
    ty = static_cast<int>(dy*b_inv[YY] + h15);
    tz = static_cast<int>(dz*b_inv[ZZ] + h15);
    tx--;
    ty--;
    tz--;

    /* Correct diff vector for translation */
    ddx = tx*box_size[XX] - dx;
    ddy = ty*box_size[YY] - dy;
    ddz = tz*box_size[ZZ] - dz;

    /* Distance squared */
    r2 = (ddx*ddx) + (ddy*ddy) + (ddz*ddz);

    *shift = XYZ2IS(tx, ty, tz);

    return r2;
}

static void add_simple(t_ns_buf * nsbuf, int nrj, int cg_j,
                       gmx_bool bHaveVdW[], int ngid, t_mdatoms *md,
                       int icg, int jgid, t_block *cgs, t_excl bexcl[],
                       int shift, t_forcerec *fr, put_in_list_t *put_in_list)
{
    if (nsbuf->nj + nrj > MAX_CG)
    {
        put_in_list(bHaveVdW, ngid, md, icg, jgid, nsbuf->ncg, nsbuf->jcg,
                    cgs->index, bexcl, shift, fr, TRUE, TRUE, fr->solvent_opt);
        /* Reset buffer contents */
        nsbuf->ncg = nsbuf->nj = 0;
    }
    nsbuf->jcg[nsbuf->ncg++] = cg_j;
    nsbuf->nj               += nrj;
}

static void ns_inner_tric(rvec x[], int icg, int *i_egp_flags,
                          int njcg, int jcg[],
                          matrix box, rvec b_inv, real rcut2,
                          t_block *cgs, t_ns_buf **ns_buf,
                          gmx_bool bHaveVdW[], int ngid, t_mdatoms *md,
                          t_excl bexcl[], t_forcerec *fr,
                          put_in_list_t *put_in_list)
{
    int       shift;
    int       j, nrj, jgid;
    int      *cginfo = fr->cginfo;
    int       cg_j, *cgindex;

    cgindex = cgs->index;
    shift   = CENTRAL;
    for (j = 0; (j < njcg); j++)
    {
        cg_j   = jcg[j];
        nrj    = cgindex[cg_j+1]-cgindex[cg_j];
        if (calc_image_tric(x[icg], x[cg_j], box, b_inv, &shift) < rcut2)
        {
            jgid  = GET_CGINFO_GID(cginfo[cg_j]);
            if (!(i_egp_flags[jgid] & EGP_EXCL))
            {
                add_simple(&ns_buf[jgid][shift], nrj, cg_j,
                           bHaveVdW, ngid, md, icg, jgid, cgs, bexcl, shift, fr,
                           put_in_list);
            }
        }
    }
}

static void ns_inner_rect(rvec x[], int icg, int *i_egp_flags,
                          int njcg, int jcg[],
                          gmx_bool bBox, rvec box_size, rvec b_inv, real rcut2,
                          t_block *cgs, t_ns_buf **ns_buf,
                          gmx_bool bHaveVdW[], int ngid, t_mdatoms *md,
                          t_excl bexcl[], t_forcerec *fr,
                          put_in_list_t *put_in_list)
{
    int       shift;
    int       j, nrj, jgid;
    int      *cginfo = fr->cginfo;
    int       cg_j, *cgindex;

    cgindex = cgs->index;
    if (bBox)
    {
        shift = CENTRAL;
        for (j = 0; (j < njcg); j++)
        {
            cg_j   = jcg[j];
            nrj    = cgindex[cg_j+1]-cgindex[cg_j];
            if (calc_image_rect(x[icg], x[cg_j], box_size, b_inv, &shift) < rcut2)
            {
                jgid  = GET_CGINFO_GID(cginfo[cg_j]);
                if (!(i_egp_flags[jgid] & EGP_EXCL))
                {
                    add_simple(&ns_buf[jgid][shift], nrj, cg_j,
                               bHaveVdW, ngid, md, icg, jgid, cgs, bexcl, shift, fr,
                               put_in_list);
                }
            }
        }
    }
    else
    {
        for (j = 0; (j < njcg); j++)
        {
            cg_j   = jcg[j];
            nrj    = cgindex[cg_j+1]-cgindex[cg_j];
            if ((rcut2 == 0) || (distance2(x[icg], x[cg_j]) < rcut2))
            {
                jgid  = GET_CGINFO_GID(cginfo[cg_j]);
                if (!(i_egp_flags[jgid] & EGP_EXCL))
                {
                    add_simple(&ns_buf[jgid][CENTRAL], nrj, cg_j,
                               bHaveVdW, ngid, md, icg, jgid, cgs, bexcl, CENTRAL, fr,
                               put_in_list);
                }
            }
        }
    }
}

/* ns_simple_core needs to be adapted for QMMM still 2005 */

static int ns_simple_core(t_forcerec *fr,
                          gmx_localtop_t *top,
                          t_mdatoms *md,
                          matrix box, rvec box_size,
                          t_excl bexcl[], int *aaj,
                          int ngid, t_ns_buf **ns_buf,
                          put_in_list_t *put_in_list, gmx_bool bHaveVdW[])
{
    int          naaj, k;
    real         rlist2;
    int          nsearch, icg, igid, nn;
    int         *cginfo;
    t_ns_buf    *nsbuf;
    /* int  *i_atoms; */
    t_block     *cgs  = &(top->cgs);
    t_blocka    *excl = &(top->excls);
    rvec         b_inv;
    int          m;
    gmx_bool     bBox, bTriclinic;
    int         *i_egp_flags;

    rlist2 = gmx::square(fr->rlist);

    bBox = (fr->ePBC != epbcNONE);
    if (bBox)
    {
        for (m = 0; (m < DIM); m++)
        {
            if (gmx_numzero(box_size[m]))
            {
                gmx_fatal(FARGS, "Dividing by zero box size!");
            }
            b_inv[m] = 1.0/box_size[m];
        }
        bTriclinic = TRICLINIC(box);
    }
    else
    {
        bTriclinic = FALSE;
    }

    cginfo = fr->cginfo;

    nsearch = 0;
    for (icg = fr->cg0; (icg < fr->hcg); icg++)
    {
        /*
           i0        = cgs->index[icg];
           nri       = cgs->index[icg+1]-i0;
           i_atoms   = &(cgs->a[i0]);
           i_eg_excl = fr->eg_excl + ngid*md->cENER[*i_atoms];
           setexcl(nri,i_atoms,excl,TRUE,bexcl);
         */
        igid        = GET_CGINFO_GID(cginfo[icg]);
        i_egp_flags = fr->egp_flags + ngid*igid;
        setexcl(cgs->index[icg], cgs->index[icg+1], excl, TRUE, bexcl);

        naaj = calc_naaj(icg, cgs->nr);
        if (bTriclinic)
        {
            ns_inner_tric(fr->cg_cm, icg, i_egp_flags, naaj, &(aaj[icg]),
                          box, b_inv, rlist2, cgs, ns_buf,
                          bHaveVdW, ngid, md, bexcl, fr, put_in_list);
        }
        else
        {
            ns_inner_rect(fr->cg_cm, icg, i_egp_flags, naaj, &(aaj[icg]),
                          bBox, box_size, b_inv, rlist2, cgs, ns_buf,
                          bHaveVdW, ngid, md, bexcl, fr, put_in_list);
        }
        nsearch += naaj;

        for (nn = 0; (nn < ngid); nn++)
        {
            for (k = 0; (k < SHIFTS); k++)
            {
                nsbuf = &(ns_buf[nn][k]);
                if (nsbuf->ncg > 0)
                {
                    put_in_list(bHaveVdW, ngid, md, icg, nn, nsbuf->ncg, nsbuf->jcg,
                                cgs->index, bexcl, k, fr, TRUE, TRUE, fr->solvent_opt);
                    nsbuf->ncg = nsbuf->nj = 0;
                }
            }
        }
        /* setexcl(nri,i_atoms,excl,FALSE,bexcl); */
        setexcl(cgs->index[icg], cgs->index[icg+1], excl, FALSE, bexcl);
    }
    close_neighbor_lists(fr, FALSE);

    return nsearch;
}

/************************************************
 *
 *    N S 5     G R I D     S T U F F
 *
 ************************************************/

static gmx_inline void get_dx_dd(int Nx, real gridx, real rc2, int xgi, real x,
                                 int ncpddc, int shift_min, int shift_max,
                                 int *g0, int *g1, real *dcx2)
{
    real dcx, tmp;
    int  g_min, g_max, shift_home;

    if (xgi < 0)
    {
        g_min = 0;
        g_max = Nx - 1;
        *g0   = 0;
        *g1   = -1;
    }
    else if (xgi >= Nx)
    {
        g_min = 0;
        g_max = Nx - 1;
        *g0   = Nx;
        *g1   = Nx - 1;
    }
    else
    {
        if (ncpddc == 0)
        {
            g_min = 0;
            g_max = Nx - 1;
        }
        else
        {
            if (xgi < ncpddc)
            {
                shift_home = 0;
            }
            else
            {
                shift_home = -1;
            }
            g_min = (shift_min == shift_home ? 0          : ncpddc);
            g_max = (shift_max == shift_home ? ncpddc - 1 : Nx - 1);
        }
        if (shift_min > 0)
        {
            *g0 = g_min;
            *g1 = g_min - 1;
        }
        else if (shift_max < 0)
        {
            *g0 = g_max + 1;
            *g1 = g_max;
        }
        else
        {
            *g0       = xgi;
            *g1       = xgi;
            dcx2[xgi] = 0;
        }
    }

    while (*g0 > g_min)
    {
        /* Check one grid cell down */
        dcx = ((*g0 - 1) + 1)*gridx - x;
        tmp = dcx*dcx;
        if (tmp >= rc2)
        {
            break;
        }
        (*g0)--;
        dcx2[*g0] = tmp;
    }

    while (*g1 < g_max)
    {
        /* Check one grid cell up */
        dcx = (*g1 + 1)*gridx - x;
        tmp = dcx*dcx;
        if (tmp >= rc2)
        {
            break;
        }
        (*g1)++;
        dcx2[*g1] = tmp;
    }
}


#define calc_dx2(XI, YI, ZI, y) (gmx::square(XI-y[XX]) + gmx::square(YI-y[YY]) + gmx::square(ZI-y[ZZ]))
#define calc_cyl_dx2(XI, YI, y) (gmx::square(XI-y[XX]) + gmx::square(YI-y[YY]))
/****************************************************
 *
 *    F A S T   N E I G H B O R  S E A R C H I N G
 *
 *    Optimized neighboursearching routine using grid
 *    at least 1x1x1, see GROMACS manual
 *
 ****************************************************/


static void get_cutoff2(t_forcerec *fr, real *rs2)
{
    *rs2 = gmx::square(fr->rlist);
}

static void init_nsgrid_lists(t_forcerec *fr, int ngid, gmx_ns_t *ns)
{
    real rs2;
    int  j;

    get_cutoff2(fr, &rs2);

    /* Short range buffers */
    snew(ns->nl_sr, ngid);
    /* Counters */
    snew(ns->nsr, ngid);

    for (j = 0; (j < ngid); j++)
    {
        snew(ns->nl_sr[j], MAX_CG);
    }
    if (debug)
    {
        fprintf(debug,
                "ns5_core: rs2 = %g (nm^2)\n",
                rs2);
    }
}

static int nsgrid_core(t_commrec *cr, t_forcerec *fr,
                       matrix box, int ngid,
                       gmx_localtop_t *top,
                       t_grid *grid,
                       t_excl bexcl[], gmx_bool *bExcludeAlleg,
                       t_mdatoms *md,
                       put_in_list_t *put_in_list,
                       gmx_bool bHaveVdW[],
                       gmx_bool bMakeQMMMnblist)
{
    gmx_ns_t     *ns;
    int         **nl_sr;
    int          *nsr;
    gmx_domdec_t *dd;
    t_block      *cgs    = &(top->cgs);
    int          *cginfo = fr->cginfo;
    /* int *i_atoms,*cgsindex=cgs->index; */
    ivec          sh0, sh1, shp;
    int           cell_x, cell_y, cell_z;
    int           d, tx, ty, tz, dx, dy, dz, cj;
#ifdef ALLOW_OFFDIAG_LT_HALFDIAG
    int           zsh_ty, zsh_tx, ysh_tx;
#endif
    int           dx0, dx1, dy0, dy1, dz0, dz1;
    int           Nx, Ny, Nz, shift = -1, j, nrj, nns, nn = -1;
    real          gridx, gridy, gridz, grid_x, grid_y;
    real         *dcx2, *dcy2, *dcz2;
    int           zgi, ygi, xgi;
    int           cg0, cg1, icg = -1, cgsnr, i0, igid, naaj, max_jcg;
    int           jcg0, jcg1, jjcg, cgj0, jgid;
    int          *grida, *gridnra, *gridind;
    rvec         *cgcm, grid_offset;
    real          r2, rs2, XI, YI, ZI, tmp1, tmp2;
    int          *i_egp_flags;
    gmx_bool      bDomDec, bTriclinicX, bTriclinicY;
    ivec          ncpddc;

    ns = fr->ns;

    bDomDec = DOMAINDECOMP(cr);
    dd      = cr->dd;

    bTriclinicX = ((YY < grid->npbcdim &&
                    (!bDomDec || dd->nc[YY] == 1) && box[YY][XX] != 0) ||
                   (ZZ < grid->npbcdim &&
                    (!bDomDec || dd->nc[ZZ] == 1) && box[ZZ][XX] != 0));
    bTriclinicY =  (ZZ < grid->npbcdim &&
                    (!bDomDec || dd->nc[ZZ] == 1) && box[ZZ][YY] != 0);

    cgsnr    = cgs->nr;

    get_cutoff2(fr, &rs2);

    nl_sr     = ns->nl_sr;
    nsr       = ns->nsr;

    /* Unpack arrays */
    cgcm    = fr->cg_cm;
    Nx      = grid->n[XX];
    Ny      = grid->n[YY];
    Nz      = grid->n[ZZ];
    grida   = grid->a;
    gridind = grid->index;
    gridnra = grid->nra;
    nns     = 0;

    gridx      = grid->cell_size[XX];
    gridy      = grid->cell_size[YY];
    gridz      = grid->cell_size[ZZ];
    grid_x     = 1/gridx;
    grid_y     = 1/gridy;
    copy_rvec(grid->cell_offset, grid_offset);
    copy_ivec(grid->ncpddc, ncpddc);
    dcx2       = grid->dcx2;
    dcy2       = grid->dcy2;
    dcz2       = grid->dcz2;

#ifdef ALLOW_OFFDIAG_LT_HALFDIAG
    zsh_ty = floor(-box[ZZ][YY]/box[YY][YY]+0.5);
    zsh_tx = floor(-box[ZZ][XX]/box[XX][XX]+0.5);
    ysh_tx = floor(-box[YY][XX]/box[XX][XX]+0.5);
    if (zsh_tx != 0 && ysh_tx != 0)
    {
        /* This could happen due to rounding, when both ratios are 0.5 */
        ysh_tx = 0;
    }
#endif

    if (fr->n_tpi)
    {
        /* We only want a list for the test particle */
        cg0 = cgsnr - 1;
    }
    else
    {
        cg0 = grid->icg0;
    }
    cg1 = grid->icg1;

    /* Set the shift range */
    for (d = 0; d < DIM; d++)
    {
        sh0[d] = -1;
        sh1[d] = 1;
        /* Check if we need periodicity shifts.
         * Without PBC or with domain decomposition we don't need them.
         */
        if (d >= ePBC2npbcdim(fr->ePBC) || (bDomDec && dd->nc[d] > 1))
        {
            shp[d] = 0;
        }
        else
        {
            if (d == XX &&
                box[XX][XX] - std::abs(box[YY][XX]) - std::abs(box[ZZ][XX]) < std::sqrt(rs2))
            {
                shp[d] = 2;
            }
            else
            {
                shp[d] = 1;
            }
        }
    }

    /* Loop over charge groups */
    for (icg = cg0; (icg < cg1); icg++)
    {
        igid = GET_CGINFO_GID(cginfo[icg]);
        /* Skip this charge group if all energy groups are excluded! */
        if (bExcludeAlleg[igid])
        {
            continue;
        }

        i0   = cgs->index[icg];

        if (bMakeQMMMnblist)
        {
            /* Skip this charge group if it is not a QM atom while making a
             * QM/MM neighbourlist
             */
            if (md->bQM[i0] == FALSE)
            {
                continue; /* MM particle, go to next particle */
            }

            /* Compute the number of charge groups that fall within the control
             * of this one (icg)
             */
            naaj    = calc_naaj(icg, cgsnr);
            jcg0    = icg;
            jcg1    = icg + naaj;
            max_jcg = cgsnr;
        }
        else
        {
            /* make a normal neighbourlist */

            if (bDomDec)
            {
                /* Get the j charge-group and dd cell shift ranges */
                dd_get_ns_ranges(cr->dd, icg, &jcg0, &jcg1, sh0, sh1);
                max_jcg = 0;
            }
            else
            {
                /* Compute the number of charge groups that fall within the control
                 * of this one (icg)
                 */
                naaj = calc_naaj(icg, cgsnr);
                jcg0 = icg;
                jcg1 = icg + naaj;

                if (fr->n_tpi)
                {
                    /* The i-particle is awlways the test particle,
                     * so we want all j-particles
                     */
                    max_jcg = cgsnr - 1;
                }
                else
                {
                    max_jcg  = jcg1 - cgsnr;
                }
            }
        }

        i_egp_flags = fr->egp_flags + igid*ngid;

        /* Set the exclusions for the atoms in charge group icg using a bitmask */
        setexcl(i0, cgs->index[icg+1], &top->excls, TRUE, bexcl);

        ci2xyz(grid, icg, &cell_x, &cell_y, &cell_z);

        /* Changed iicg to icg, DvdS 990115
         * (but see consistency check above, DvdS 990330)
         */
#ifdef NS5DB
        fprintf(log, "icg=%5d, naaj=%5d, cell %d %d %d\n",
                icg, naaj, cell_x, cell_y, cell_z);
#endif
        /* Loop over shift vectors in three dimensions */
        for (tz = -shp[ZZ]; tz <= shp[ZZ]; tz++)
        {
            ZI = cgcm[icg][ZZ]+tz*box[ZZ][ZZ];
            /* Calculate range of cells in Z direction that have the shift tz */
            zgi = cell_z + tz*Nz;
            get_dx_dd(Nz, gridz, rs2, zgi, ZI-grid_offset[ZZ],
                      ncpddc[ZZ], sh0[ZZ], sh1[ZZ], &dz0, &dz1, dcz2);
            if (dz0 > dz1)
            {
                continue;
            }
            for (ty = -shp[YY]; ty <= shp[YY]; ty++)
            {
                YI = cgcm[icg][YY]+ty*box[YY][YY]+tz*box[ZZ][YY];
                /* Calculate range of cells in Y direction that have the shift ty */
                if (bTriclinicY)
                {
                    ygi = (int)(Ny + (YI - grid_offset[YY])*grid_y) - Ny;
                }
                else
                {
                    ygi = cell_y + ty*Ny;
                }
                get_dx_dd(Ny, gridy, rs2, ygi, YI-grid_offset[YY],
                          ncpddc[YY], sh0[YY], sh1[YY], &dy0, &dy1, dcy2);
                if (dy0 > dy1)
                {
                    continue;
                }
                for (tx = -shp[XX]; tx <= shp[XX]; tx++)
                {
                    XI = cgcm[icg][XX]+tx*box[XX][XX]+ty*box[YY][XX]+tz*box[ZZ][XX];
                    /* Calculate range of cells in X direction that have the shift tx */
                    if (bTriclinicX)
                    {
                        xgi = (int)(Nx + (XI - grid_offset[XX])*grid_x) - Nx;
                    }
                    else
                    {
                        xgi = cell_x + tx*Nx;
                    }
                    get_dx_dd(Nx, gridx, rs2, xgi, XI-grid_offset[XX],
                              ncpddc[XX], sh0[XX], sh1[XX], &dx0, &dx1, dcx2);
                    if (dx0 > dx1)
                    {
                        continue;
                    }
                    /* Get shift vector */
                    shift = XYZ2IS(tx, ty, tz);
#ifdef NS5DB
                    range_check(shift, 0, SHIFTS);
#endif
                    for (nn = 0; (nn < ngid); nn++)
                    {
                        nsr[nn]      = 0;
                    }
#ifdef NS5DB
                    fprintf(log, "shift: %2d, dx0,1: %2d,%2d, dy0,1: %2d,%2d, dz0,1: %2d,%2d\n",
                            shift, dx0, dx1, dy0, dy1, dz0, dz1);
                    fprintf(log, "cgcm: %8.3f  %8.3f  %8.3f\n", cgcm[icg][XX],
                            cgcm[icg][YY], cgcm[icg][ZZ]);
                    fprintf(log, "xi:   %8.3f  %8.3f  %8.3f\n", XI, YI, ZI);
#endif
                    for (dx = dx0; (dx <= dx1); dx++)
                    {
                        tmp1 = rs2 - dcx2[dx];
                        for (dy = dy0; (dy <= dy1); dy++)
                        {
                            tmp2 = tmp1 - dcy2[dy];
                            if (tmp2 > 0)
                            {
                                for (dz = dz0; (dz <= dz1); dz++)
                                {
                                    if (tmp2 > dcz2[dz])
                                    {
                                        /* Find grid-cell cj in which possible neighbours are */
                                        cj   = xyz2ci(Ny, Nz, dx, dy, dz);

                                        /* Check out how many cgs (nrj) there in this cell */
                                        nrj  = gridnra[cj];

                                        /* Find the offset in the cg list */
                                        cgj0 = gridind[cj];

                                        /* Check if all j's are out of range so we
                                         * can skip the whole cell.
                                         * Should save some time, especially with DD.
                                         */
                                        if (nrj == 0 ||
                                            (grida[cgj0] >= max_jcg &&
                                             (grida[cgj0] >= jcg1 || grida[cgj0+nrj-1] < jcg0)))
                                        {
                                            continue;
                                        }

                                        /* Loop over cgs */
                                        for (j = 0; (j < nrj); j++)
                                        {
                                            jjcg = grida[cgj0+j];

                                            /* check whether this guy is in range! */
                                            if ((jjcg >= jcg0 && jjcg < jcg1) ||
                                                (jjcg < max_jcg))
                                            {
                                                r2 = calc_dx2(XI, YI, ZI, cgcm[jjcg]);
                                                if (r2 < rs2)
                                                {
                                                    /* jgid = gid[cgsatoms[cgsindex[jjcg]]]; */
                                                    jgid = GET_CGINFO_GID(cginfo[jjcg]);
                                                    /* check energy group exclusions */
                                                    if (!(i_egp_flags[jgid] & EGP_EXCL))
                                                    {
                                                        if (nsr[jgid] >= MAX_CG)
                                                        {
                                                            /* Add to short-range list */
                                                            put_in_list(bHaveVdW, ngid, md, icg, jgid,
                                                                        nsr[jgid], nl_sr[jgid],
                                                                        cgs->index, /* cgsatoms, */ bexcl,
                                                                        shift, fr, TRUE, TRUE, fr->solvent_opt);
                                                            nsr[jgid] = 0;
                                                        }
                                                        nl_sr[jgid][nsr[jgid]++] = jjcg;
                                                    }
                                                }
                                                nns++;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    /* CHECK whether there is anything left in the buffers */
                    for (nn = 0; (nn < ngid); nn++)
                    {
                        if (nsr[nn] > 0)
                        {
                            put_in_list(bHaveVdW, ngid, md, icg, nn, nsr[nn], nl_sr[nn],
                                        cgs->index, /* cgsatoms, */ bexcl,
                                        shift, fr, TRUE, TRUE, fr->solvent_opt);
                        }
                    }
                }
            }
        }
        setexcl(cgs->index[icg], cgs->index[icg+1], &top->excls, FALSE, bexcl);
    }
    /* No need to perform any left-over force calculations anymore (as we used to do here)
     * since we now save the proper long-range lists for later evaluation.
     */

    /* Close neighbourlists */
    close_neighbor_lists(fr, bMakeQMMMnblist);

    return nns;
}

static void ns_realloc_natoms(gmx_ns_t *ns, int natoms)
{
    int i;

    if (natoms > ns->nra_alloc)
    {
        ns->nra_alloc = over_alloc_dd(natoms);
        srenew(ns->bexcl, ns->nra_alloc);
        for (i = 0; i < ns->nra_alloc; i++)
        {
            ns->bexcl[i] = 0;
        }
    }
}

void init_ns(FILE *fplog, const t_commrec *cr,
             gmx_ns_t *ns, t_forcerec *fr,
             const gmx_mtop_t *mtop)
{
    int      mt, icg, nr_in_cg, maxcg, i, j, jcg, ngid, ncg;
    t_block *cgs;

    /* Compute largest charge groups size (# atoms) */
    nr_in_cg = 1;
    for (mt = 0; mt < mtop->nmoltype; mt++)
    {
        cgs = &mtop->moltype[mt].cgs;
        for (icg = 0; (icg < cgs->nr); icg++)
        {
            nr_in_cg = std::max(nr_in_cg, (int)(cgs->index[icg+1]-cgs->index[icg]));
        }
    }

    /* Verify whether largest charge group is <= max cg.
     * This is determined by the type of the local exclusion type
     * Exclusions are stored in bits. (If the type is not large
     * enough, enlarge it, unsigned char -> unsigned short -> unsigned long)
     */
    maxcg = sizeof(t_excl)*8;
    if (nr_in_cg > maxcg)
    {
        gmx_fatal(FARGS, "Max #atoms in a charge group: %d > %d\n",
                  nr_in_cg, maxcg);
    }

    ngid = mtop->groups.grps[egcENER].nr;
    snew(ns->bExcludeAlleg, ngid);
    for (i = 0; i < ngid; i++)
    {
        ns->bExcludeAlleg[i] = TRUE;
        for (j = 0; j < ngid; j++)
        {
            if (!(fr->egp_flags[i*ngid+j] & EGP_EXCL))
            {
                ns->bExcludeAlleg[i] = FALSE;
            }
        }
    }

    if (fr->bGrid)
    {
        /* Grid search */
        ns->grid = init_grid(fplog, fr);
        init_nsgrid_lists(fr, ngid, ns);
    }
    else
    {
        /* Simple search */
        snew(ns->ns_buf, ngid);
        for (i = 0; (i < ngid); i++)
        {
            snew(ns->ns_buf[i], SHIFTS);
        }
        ncg = ncg_mtop(mtop);
        snew(ns->simple_aaj, 2*ncg);
        for (jcg = 0; (jcg < ncg); jcg++)
        {
            ns->simple_aaj[jcg]     = jcg;
            ns->simple_aaj[jcg+ncg] = jcg;
        }
    }

    /* Create array that determines whether or not atoms have VdW */
    snew(ns->bHaveVdW, fr->ntype);
    for (i = 0; (i < fr->ntype); i++)
    {
        for (j = 0; (j < fr->ntype); j++)
        {
            ns->bHaveVdW[i] = (ns->bHaveVdW[i] ||
                               (fr->bBHAM ?
                                ((BHAMA(fr->nbfp, fr->ntype, i, j) != 0) ||
                                 (BHAMB(fr->nbfp, fr->ntype, i, j) != 0) ||
                                 (BHAMC(fr->nbfp, fr->ntype, i, j) != 0)) :
                                ((C6(fr->nbfp, fr->ntype, i, j) != 0) ||
                                 (C12(fr->nbfp, fr->ntype, i, j) != 0))));
        }
    }
    if (debug)
    {
        pr_bvec(debug, 0, "bHaveVdW", ns->bHaveVdW, fr->ntype, TRUE);
    }

    ns->nra_alloc = 0;
    ns->bexcl     = nullptr;
    if (!DOMAINDECOMP(cr))
    {
        ns_realloc_natoms(ns, mtop->natoms);
    }

    ns->nblist_initialized = FALSE;

    /* nbr list debug dump */
    {
        char *ptr = getenv("GMX_DUMP_NL");
        if (ptr)
        {
            ns->dump_nl = strtol(ptr, nullptr, 10);
            if (fplog)
            {
                fprintf(fplog, "GMX_DUMP_NL = %d", ns->dump_nl);
            }
        }
        else
        {
            ns->dump_nl = 0;
        }
    }
}


int search_neighbours(FILE *log, t_forcerec *fr,
                      matrix box,
                      gmx_localtop_t *top,
                      gmx_groups_t *groups,
                      t_commrec *cr,
                      t_nrnb *nrnb, t_mdatoms *md,
                      gmx_bool bFillGrid)
{
    t_block            *cgs = &(top->cgs);
    rvec                box_size, grid_x0, grid_x1;
    int                 m, ngid;
    real                min_size, grid_dens;
    int                 nsearch;
    gmx_bool            bGrid;
    int                 start, end;
    gmx_ns_t           *ns;
    t_grid             *grid;
    gmx_domdec_zones_t *dd_zones;
    put_in_list_t      *put_in_list;

    ns = fr->ns;

    /* Set some local variables */
    bGrid = fr->bGrid;
    ngid  = groups->grps[egcENER].nr;

    for (m = 0; (m < DIM); m++)
    {
        box_size[m] = box[m][m];
    }

    if (fr->ePBC != epbcNONE)
    {
        if (gmx::square(fr->rlist) >= max_cutoff2(fr->ePBC, box))
        {
            gmx_fatal(FARGS, "One of the box vectors has become shorter than twice the cut-off length or box_yy-|box_zy| or box_zz has become smaller than the cut-off.");
        }
        if (!bGrid)
        {
            min_size = std::min(box_size[XX], std::min(box_size[YY], box_size[ZZ]));
            if (2*fr->rlist >= min_size)
            {
                gmx_fatal(FARGS, "One of the box diagonal elements has become smaller than twice the cut-off length.");
            }
        }
    }

    if (DOMAINDECOMP(cr))
    {
        ns_realloc_natoms(ns, cgs->index[cgs->nr]);
    }

    /* Reset the neighbourlists */
    reset_neighbor_lists(fr);

    if (bGrid && bFillGrid)
    {

        grid = ns->grid;
        if (DOMAINDECOMP(cr))
        {
            dd_zones = domdec_zones(cr->dd);
        }
        else
        {
            dd_zones = nullptr;

            get_nsgrid_boundaries(grid->nboundeddim, box, nullptr, nullptr, nullptr, nullptr,
                                  cgs->nr, fr->cg_cm, grid_x0, grid_x1, &grid_dens);

            grid_first(log, grid, nullptr, nullptr, box, grid_x0, grid_x1,
                       fr->rlist, grid_dens);
        }

        start = 0;
        end   = cgs->nr;

        if (DOMAINDECOMP(cr))
        {
            end = cgs->nr;
            fill_grid(dd_zones, grid, end, -1, end, fr->cg_cm);
            grid->icg0 = 0;
            grid->icg1 = dd_zones->izone[dd_zones->nizone-1].cg1;
        }
        else
        {
            fill_grid(nullptr, grid, cgs->nr, fr->cg0, fr->hcg, fr->cg_cm);
            grid->icg0 = fr->cg0;
            grid->icg1 = fr->hcg;
        }

        calc_elemnr(grid, start, end, cgs->nr);
        calc_ptrs(grid);
        grid_last(grid, start, end, cgs->nr);

        if (gmx_debug_at)
        {
            check_grid(grid);
            print_grid(debug, grid);
        }
    }
    else if (fr->n_tpi)
    {
        /* Set the grid cell index for the test particle only.
         * The cell to cg index is not corrected, but that does not matter.
         */
        fill_grid(nullptr, ns->grid, fr->hcg, fr->hcg-1, fr->hcg, fr->cg_cm);
    }

    if (!fr->ns->bCGlist)
    {
        put_in_list = put_in_list_at;
    }
    else
    {
        put_in_list = put_in_list_cg;
    }

    /* Do the core! */
    if (bGrid)
    {
        grid    = ns->grid;
        nsearch = nsgrid_core(cr, fr, box, ngid, top,
                              grid, ns->bexcl, ns->bExcludeAlleg,
                              md, put_in_list, ns->bHaveVdW,
                              FALSE);

        /* neighbour searching withouth QMMM! QM atoms have zero charge in
         * the classical calculation. The charge-charge interaction
         * between QM and MM atoms is handled in the QMMM core calculation
         * (see QMMM.c). The VDW however, we'd like to compute classically
         * and the QM MM atom pairs have just been put in the
         * corresponding neighbourlists. in case of QMMM we still need to
         * fill a special QMMM neighbourlist that contains all neighbours
         * of the QM atoms. If bQMMM is true, this list will now be made:
         */
        if (fr->bQMMM && fr->qr->QMMMscheme != eQMMMschemeoniom)
        {
            nsearch += nsgrid_core(cr, fr, box, ngid, top,
                                   grid, ns->bexcl, ns->bExcludeAlleg,
                                   md, put_in_list_qmmm, ns->bHaveVdW,
                                   TRUE);
        }
    }
    else
    {
        nsearch = ns_simple_core(fr, top, md, box, box_size,
                                 ns->bexcl, ns->simple_aaj,
                                 ngid, ns->ns_buf, put_in_list, ns->bHaveVdW);
    }

    inc_nrnb(nrnb, eNR_NS, nsearch);

    return nsearch;
}
