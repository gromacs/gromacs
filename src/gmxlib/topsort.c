/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-  */

/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "smalloc.h"

static bool ip_pert(int ftype,t_iparams *ip)
{
    bool bPert;
    int  i;

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
    case F_PDIHS:
    case F_ANGRES:
    case F_ANGRESZ:
        bPert = (ip->pdihs.phiA != ip->pdihs.phiB ||
                 ip->pdihs.cpA  != ip->pdihs.cpB);
        break;
    case F_RBDIHS:
        bPert = FALSE;
        for(i=0; i<NR_RBDIHS; i++)
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
        for(i=0; i<DIM; i++)
        {
            if (ip->posres.pos0A[i] != ip->posres.pos0B[i] ||
                ip->posres.fcA[i]   != ip->posres.fcB[i])
            {
                bPert = TRUE;
            }
        }
        break;
    case F_LJ14:
        bPert = (ip->lj14.c6A  != ip->lj14.c6B ||
                 ip->lj14.c12A != ip->lj14.c12B);
        break;
    default:
        bPert = FALSE;
        gmx_fatal(FARGS,"Function type %s not implemented in ip_pert",
                  interaction_function[ftype].longname);
    }

    return bPert;
}

bool gmx_mtop_bondeds_free_energy(gmx_mtop_t *mtop)
{
    gmx_ffparams_t *ffparams;
    int  i,ftype;
    bool bPert;

    ffparams = &mtop->ffparams;
    
    /* Loop over all the function types and compare the A/B parameters */
    bPert = FALSE;
    for(i=0; i<ffparams->ntypes; i++)
    {
        ftype = ffparams->functype[i];
        if (interaction_function[ftype].flags & IF_BOND)
        {
            if (ip_pert(ftype,&ffparams->iparams[i]))
            {
                bPert = TRUE;
            }
        }
    }

    return (bPert ? ilsortFE_UNSORTED : ilsortNO_FE);
}

void gmx_sort_ilist_fe(t_idef *idef)
{
    int  ftype,nral,i,ic,ib,a;
    t_iparams *iparams;
    t_ilist *ilist;
    t_iatom *iatoms;
    bool bPert;
    t_iatom *iabuf;
    int  iabuf_nalloc;

    iabuf_nalloc = 0;
    iabuf        = NULL;
    
    iparams = idef->iparams;

    for(ftype=0; ftype<F_NRE; ftype++)
    {
        if (interaction_function[ftype].flags & IF_BOND)
        {
            ilist = &idef->il[ftype];
            iatoms = ilist->iatoms;
            nral  = NRAL(ftype);
            ic = 0;
            ib = 0;
            i  = 0;
            while (i < ilist->nr)
            {
                /* The first element of ia gives the type */
                if (ip_pert(ftype,&iparams[iatoms[i]]))
                {
                    /* Copy to the perturbed buffer */
                    if (ib + 1 + nral > iabuf_nalloc)
                    {
                        iabuf_nalloc = over_alloc_large(ib+1+nral);
                        srenew(iabuf,iabuf_nalloc);
                    }
                    for(a=0; a<1+nral; a++)
                    {
                        iabuf[ib++] = iatoms[i++];
                    }
                }
                else
                {
                    /* Copy in place */
                    for(a=0; a<1+nral; a++)
                    {
                        iatoms[ic++] = iatoms[i++];
                    }
                }
            }
            /* Now we now the number of non-perturbed interactions */
            ilist->nr_nonperturbed = ic;
            
            /* Copy the buffer with perturbed interactions to the ilist */
            for(a=0; a<ib; a++)
            {
                iatoms[ic++] = iabuf[a];
            }

            if (debug)
            {
                fprintf(debug,"%s non-pert %d pert %d\n",
                        interaction_function[ftype].longname,
                        ilist->nr_nonperturbed,
                        ilist->nr-ilist->nr_nonperturbed);
            }
        }
    }

    sfree(iabuf);

    idef->ilsort = ilsortFE_SORTED;
}
