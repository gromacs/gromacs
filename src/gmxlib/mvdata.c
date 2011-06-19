/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sysstuff.h>
#include <string.h>
#include "typedefs.h"
#include "main.h"
#include "mvdata.h"
#include "network.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "symtab.h"
#include "vec.h"
#include "tgroup.h"

#ifndef GMX_MPI

#else

typedef enum {eptPackSize, eptPack, eptUnpack} ePackType;

/* Helper function */
void *snew_pack(ePackType packtype, void *d, int nr)
{
    if ((packtype) == bcUnpack)
    {
        snew((d),(nr));
    }
}

typedef struct
{
    const t_commrec *cr;
    char *buffer;
    int size_of_buffer;
    int position;
} t_gmx_packed;
typedef t_gmx_packed *gmx_packed_t;

gmx_packed_t gmx_packed_init(const t_commrec *cr)
{
    gmx_packed_t packed;
    snew(packed, 1);
    packed->cr = cr;
    return packed;
}

int gmx_packed_allocate(gmx_packed_t packed)
{
    snew(packed->buffer, packed->size_of_buffer);
    packed->position = 0;
    return GMX_SUCCESS;
}

int gmx_packed_refresh(gmx_packed_t packed)
{
    sfree(packed->buffer);
    packed->size_of_buffer = 0;
    packed->position = 0;
    return GMX_SUCCESS;
}

int gmx_packed_destroy(gmx_packed_t packed)
{
    sfree(packed->buffer);
    sfree(packed);
    return GMX_SUCCESS;
}

typedef int t_pack_function(gmx_packed_t, void *, int, MPI_Datatype);

int gmx_pack_size(gmx_packed_t packed, void *data, int nr, MPI_Datatype mpitype)
{
    int packed_size = 0;
    MPI_Pack_size(nr, mpitype, packed->cr->mpi_comm_group, &packed_size);
    packed->size_of_buffer += packed_size;
    return GMX_SUCCESS;
}

int gmx_pack(gmx_packed_t packed, void *data, int nr, MPI_Datatype mpitype)
{
    if (0 < nr)
    {
        MPI_Pack(data, nr, mpitype, packed->buffer, packed->size_of_buffer, &packed->position, (cr)->mpi_comm_group);
    }
    return GMX_SUCCESS;
}

int gmx_unpack(gmx_packed_t packed, void *data, int nr, MPI_Datatype mpitype)
{
    if (0 < nr)
    {
        MPI_Unpack(packed->buffer, packed->size_of_buffer, &packed->position, data, nr, mpitype, packed->cr->mpi_comm_group);
    }
    return GMX_SUCCESS;
}

int do_packed_bcast(gmx_packed_t packed)
{
    /* Probably the test for (nr) > 0 is only needed on BlueGene(/L),
     * where IBM's MPI_Bcast will segfault after dereferencing a null
     * pointer, even when no data is to be transferred. */
    if (0 < packed->position)
    {
        MPI_Bcast(packed->buffer, packed->position, MPI_PACKED, MASTERRANK(packed->cr), packed->cr->mpi_comm_group);
    }
    return GMX_SUCCESS;
}

/* Broadcast a t_state struct to the other processors in the
 * simulation. */
void bcast_state(const t_commrec *cr, t_state *state, gmx_bool bAlloc)
{
    int i, nnht, nnhtp;
    gmx_packed_t packed = gmx_packed_init(cr);
    int num_to_pack = 7;
    t_pack_function generic_pack_function = MASTER(cr) ? gmx_pack : gmx_unpack;
    ePackType packtype = eptPackSize;
    /* First, send around various integers that determine sizes of
     * arrays, etc. */

    /* How much will we pack? */
    if (state->flags & estDISRE_RM3TAV)
    {
        num_to_pack++;
    }
    if (state->flags & estORIRE_DTAV)
    {
        num_to_pack++;
    }
    /* The next line has to be kept in sync with number of function
     * MPI_INT lines below */
    gmx_pack_size(packed, num_to_pack, MPI_INT);
    
    /* Allocate */
    gmx_packed_allocate(packed);

    if (MASTER(cr))
    {
        packtype = eptPack;
    }
    else
    {
        packtype = eptUnpack;
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Do packing or unpacking, as appropriate */
    generic_pack_function(packed, state->natoms, 1, MPI_INT);
    generic_pack_function(packed, state->ngtc, 1, MPI_INT);
    generic_pack_function(packed, state->nnhpres, 1, MPI_INT);
    generic_pack_function(packed, state->nhchainlength, 1, MPI_INT);
    generic_pack_function(packed, state->nrng, 1, MPI_INT);
    generic_pack_function(packed, state->nrngi, 1, MPI_INT);
    generic_pack_function(packed, state->flags, 1, MPI_INT);
    if (state->flags & estDISRE_RM3TAV)
    {
        generic_pack_function(packed, state->hist.ndisrepairs, MPI_INT);
    }
    if (state->flags & estORIRE_DTAV)
    {
        generic_pack_function(packed, state->hist.norire_Dtav, MPI_INT);
    }

    if (MASTER(cr))
    {
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Clean up, prepare for next round of packing */
    gmx_packed_refresh(packed);
    packtype = eptPackSize;

    /* Now do some allocation and send the remainder of the state
     * data. */

    nnht = (state->ngtc)*(state->nhchainlength); 
    nnhtp = (state->nnhpres)*(state->nhchainlength); 

    if (MASTER(cr))
    {
        bAlloc = FALSE;
    }
    if (bAlloc)
    {
        state->nalloc = state->natoms;
    }

    /* How much will we pack? */
    for(i = 0; i < estNR; i++)
    {
        if (state->flags & (1<<i))
        {
            switch (i)
            {
            case estLAMBDA:      gmx_pack_size(packed, 1, GMX_MPI_REAL); break;
            case estBOX:         gmx_pack_size(packed, DIM*DIM, GMX_MPI_REAL); break;
            case estBOX_REL:     gmx_pack_size(packed, DIM*DIM, GMX_MPI_REAL); break;
            case estBOXV:        gmx_pack_size(packed, DIM*DIM, GMX_MPI_REAL); break;
            case estPRES_PREV:   gmx_pack_size(packed, DIM*DIM, GMX_MPI_REAL); break;
            case estSVIR_PREV:   gmx_pack_size(packed, DIM*DIM, GMX_MPI_REAL); break;
            case estFVIR_PREV:   gmx_pack_size(packed, DIM*DIM, GMX_MPI_REAL); break;
            case estNH_XI:       gmx_pack_size(packed, nnht, MPI_DOUBLE); break;
            case estNH_VXI:      gmx_pack_size(packed, nnht, MPI_DOUBLE); break;
            case estNHPRES_XI:   gmx_pack_size(packed, nnhtp, MPI_DOUBLE); break;
            case estNHPRES_VXI:  gmx_pack_size(packed, nnhtp, MPI_DOUBLE); break;
            case estTC_INT:      gmx_pack_size(packed, state->ngtc, MPI_DOUBLE); break;
            case estVETA:        gmx_pack_size(packed, 1, GMX_MPI_REAL); break;
            case estVOL0:        gmx_pack_size(packed, 1, GMX_MPI_REAL); break;
            case estX:           break;
            case estV:           break;
            case estSDX:         break;
            case estCGP:         break;
            case estLD_RNG:      if(1 == state->nrngi) gmx_pack_size(packed, state->nrng, MPI_UNSIGNED); break;
            case estLD_RNGI:     if(1 == state->nrngi) gmx_pack_size(packed, state->nrngi, MPI_INT); break;
            case estDISRE_INITF: gmx_pack_size(packed, 1, GMX_MPI_REAL); break;
            case estDISRE_RM3TAV:gmx_pack_size(packed, state->hist.ndisrepairs, GMX_MPI_REAL); break;
            case estORIRE_INITF: gmx_pack_size(packed, 1, GMX_MPI_REAL); break;
            case estORIRE_DTAV:  gmx_pack_size(packed, state->hist.norire_Dtav, GMX_MPI_REAL);
                break;
            default:
                gmx_fatal(FARGS,
                          "Communication is not implemented for %s in bcast_state",
                          est_names[i]);
            }
        }
    }

    /* Allocate */
    gmx_packed_allocate(packed);

    if (MASTER(cr))
    {
        packtype = eptPack;
    }
    else
    {
        packtype = eptUnpack;
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Do packing or unpacking, as appropriate */
    generic_pack_function(packed, state->natoms, 1, MPI_INT);
    for(i = 0; i < estNR; i++)
    {
        if (state->flags & (1<<i))
        {
            switch (i)
            {
            case estLAMBDA:
                generic_pack_function(packed, state->lambda, 1, GMX_MPI_REAL);
                break;
            case estBOX:
                generic_pack_function(packed, state->box, DIM*DIM, GMX_MPI_REAL);
                break;
            case estBOX_REL:
                generic_pack_function(packed, state->box_rel, DIM*DIM, GMX_MPI_REAL);
                break;
            case estBOXV:
                generic_pack_function(packed, state->boxv, DIM*DIM, GMX_MPI_REAL);
                break;
            case estPRES_PREV:
                generic_pack_function(packed, state->pres_prev, DIM*DIM, GMX_MPI_REAL);
                break;
            case estSVIR_PREV:
                generic_pack_function(packed, state->svir_prev, DIM*DIM, GMX_MPI_REAL);
                break;
            case estFVIR_PREV:
                generic_pack_function(packed, state->fvir_prev, DIM*DIM, GMX_MPI_REAL);
                break;
            case estNH_XI:
                if (bAlloc)
                {
                    snew(state->nosehoover_xi, nnht);
                }
                generic_pack_function(packed, state->nosehoover_xi, nnht, MPI_DOUBLE);
                break;
            case estNH_VXI:
                if (bAlloc)
                {
                    snew(state->nosehoover_vxi, nnht);
                }
                generic_pack_function(packed, state->nosehoover_vxi, nnht, MPI_DOUBLE);
                break;
            case estNHPRES_XI:
                if (bAlloc)
                {
                    snew(state->nhpres_xi, nnhtp);
                }
                generic_pack_function(packed, state->nhpres_xi, nnhtp, MPI_DOUBLE);
                break;
            case estNHPRES_VXI:
                if (bAlloc)
                {
                    snew(state->nhpres_vxi, nnhtp);
                }
                generic_pack_function(packed, state->nhpres_vxi, nnhtp, MPI_DOUBLE);
                break;
            case estTC_INT:
                if (bAlloc)
                {
                    snew(state->them_integral, state->ngtc);
                }
                generic_pack_function(packed, state->ngtc, state->therm_integral, 1, MPI_DOUBLE);
                break;
            case estVETA:
                generic_pack_function(packed, state->veta, 1, GMX_MPI_REAL);
                break;
            case estVOL0:
                generic_pack_function(packed, state->vol0, 1, GMX_MPI_REAL);
                break;
                /* The next four items are generally large enough to be communicated by themselves */
            case estX:
                if (bAlloc)
                {
                    snew(state->x, state->natoms);
                }
                generic_pack_function(packed, state->x, state->natoms * DIM, GMX_MPI_REAL);
                break;
            case estV:
                if (bAlloc)
                {
                    snew(state->v, state->natoms);
                }
                generic_pack_function(packed, state->v, state->natoms * DIM, GMX_MPI_REAL);
                break;
            case estSDX:
                if (bAlloc)
                {
                    snew(state->sd_X, state->natoms);
                }
                generic_pack_function(packed, state->sd_X, state->natoms * DIM, GMX_MPI_REAL);
                break;
            case estCGP:
                if (bAlloc)
                {
                    snew(state->cg_p, state->natoms);
                }
                generic_pack_function(packed, state->cg_p, state->natoms * DIM, GMX_MPI_REAL);
                break;
            case estLD_RNG:
                if(state->nrngi == 1)
                {
                    if (bAlloc)
                    {
                        snew(state->ld_rng, state->nrng);
                    }
                    generic_pack_function(packed, state->ld_rng, state->nrng, MPI_UNSIGNED);
                    break;
                }
            case estLD_RNGI:
                if(state->nrngi == 1)
                {
                    if (bAlloc)
                    {
                        snew(state->ld_rngi, state->nrngi);
                    }
                    generic_pack_function(packed, state->ld_rngi, state->nrngi, MPI_INT); 
                    break;
                }
            case estDISRE_INITF:
                generic_pack_function(packed, state->hist.disre_initf, 1, GMX_MPI_REAL);
                break;
            case estDISRE_RM3TAV:
                if (bAlloc)
                {
                    snew(state->hist.disre_rm3tav, state->hist.ndisrepairs);
                }
                generic_pack_function(packed, state->hist.disre_rm3tav, state->hist.ndisrepairs, GMX_MPI_REAL);
                break;
            case estORIRE_INITF:
                generic_pack_function(packed, state->hist.orire_initf, 1, GMX_MPI_REAL);
                break;
            case estORIRE_DTAV:
                if (bAlloc)
                {
                    snew(state->hist.orire_Dtav, state->hist.norire_Dtav);
                }
                generic_pack_function(packed, state->hist.orire_Dtav, state->hist.norire_Dtav, GMX_MPI_REAL);
                break;
            default:
                gmx_fatal(FARGS,
                          "Communication is not implemented for %s in bcast_state",
                          est_names[i]);
            }
        }
    }

    if (MASTER(cr))
    {
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Clean up */
    gmx_packed_destroy(packed);

/*
  for(i=0; i<estNR; i++)
  {
  if (state->flags & (1<<i))
  {
  switch (i)
  {
  case estLAMBDA:  block_bc(cr,state->lambda); break;
  case estBOX:     block_bc(cr,state->box); break;
  case estBOX_REL: block_bc(cr,state->box_rel); break;
  case estBOXV:    block_bc(cr,state->boxv); break;
  case estPRES_PREV: block_bc(cr,state->pres_prev); break;
  case estSVIR_PREV: block_bc(cr,state->svir_prev); break;
  case estFVIR_PREV: block_bc(cr,state->fvir_prev); break;
  case estNH_XI:   nblock_abc(cr,nnht,state->nosehoover_xi); break;
  case estNH_VXI:  nblock_abc(cr,nnht,state->nosehoover_vxi); break;
  case estNHPRES_XI:   nblock_abc(cr,nnhtp,state->nhpres_xi); break;
  case estNHPRES_VXI:  nblock_abc(cr,nnhtp,state->nhpres_vxi); break;
  case estTC_INT:  nblock_abc(cr,state->ngtc,state->therm_integral); break;
  case estVETA:    block_bc(cr,state->veta); break;
  case estVOL0:    block_bc(cr,state->vol0); break;
  case estX:       nblock_abc(cr,state->natoms,state->x); break;
  case estV:       nblock_abc(cr,state->natoms,state->v); break;
  case estSDX:     nblock_abc(cr,state->natoms,state->sd_X); break;
  case estCGP:     nblock_abc(cr,state->natoms,state->cg_p); break;
  case estLD_RNG:  if(state->nrngi == 1) nblock_abc(cr,state->nrng,state->ld_rng); break;
  case estLD_RNGI: if(state->nrngi == 1) nblock_abc(cr,state->nrngi,state->ld_rngi); break;
  case estDISRE_INITF: block_bc(cr,state->hist.disre_initf); break;
  case estDISRE_RM3TAV:
  block_bc(cr,state->hist.ndisrepairs);
  nblock_abc(cr,state->hist.ndisrepairs,state->hist.disre_rm3tav);
  break;
  case estORIRE_INITF: block_bc(cr,state->hist.orire_initf); break;
  case estORIRE_DTAV:
  block_bc(cr,state->hist.norire_Dtav);
  nblock_abc(cr,state->hist.norire_Dtav,state->hist.orire_Dtav);
  break;
  default:
  gmx_fatal(FARGS,
  "Communication is not implemented for %s in bcast_state",
  est_names[i]);
  }
  }
  }
*/
}

static void generic_pack_string(gmx_packed_t packed, ePackType packtype, t_symtab *symtab, char ***s)
{
    int handle;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    if (eptPack == packtype)
    {
        handle = lookup_symtab(symtab, *s);
    }
    generic_pack_function(packed, handle, 1, MPI_INT);
    if (eptUnpack == packtype)
    {
        *s = get_symtab_handle(symtab, handle);
    }
/*
  if (MASTER(cr)) {
    handle = lookup_symtab(symtab,*s);
  }
  block_bc(cr,handle);
  if (!MASTER(cr)) {
    *s = get_symtab_handle(symtab,handle);
  }
*/
}

static void generic_pack_strings(gmx_packed_t packed, ePackType packtype, t_symtab *symtab, int nr, char ****nm)
{
    int  i;
    int  *handle;
    char ***NM;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    if (eptPack == packtype)
    {
        snew(handle, nr);
        for(i = 0; i < nr; i++)
        {
            handle[i] = lookup_symtab(symtab, NM[i]);
        }
    }
    generic_pack_function(packed, handle, nr, MPI_INT);
    if (eptUnpack == packtype)
    {
        snew_pack(cr, *nm, nr);
        for (i = 0; i < nr; i++)
        {
            (*nm)[i] = get_symtab_handle(symtab, handle[i]);
        }
    }
    if (eptPackSize != packtype)
    {
        sfree(handle);
    }
/*
  snew(handle,nr);
  if (MASTER(cr)) {
    NM = *nm;
    for(i=0; (i<nr); i++)
      handle[i] = lookup_symtab(symtab,NM[i]);
  }
  nblock_bc(cr,nr,handle);

  if (!MASTER(cr)) {
    snew_pack(cr,*nm,nr);
    NM = *nm;
    for (i=0; (i<nr); i++) 
      (*nm)[i] = get_symtab_handle(symtab,handle[i]);
  }
  sfree(handle);
*/
}

static void generic_pack_strings_resinfo(gmx_packed_t packed, ePackType packtype, t_symtab *symtab,
                               int nr, t_resinfo *resinfo)
{
    int  i;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    for(i = 0; i < nr; i++)
    {
        int  handle = -1;
        if (eptPack == packtype)
        {
            handle = lookup_symtab(symtab, resinfo[i].name);
        }
        generic_pack_function(packed, &handle, 1, MPI_INT);
        if (eptUnpack == packtype)
        {
            resinfo[i].name = get_symtab_handle(symtab, handle);
        }
    }
/*
  snew(handle,nr);
  if (MASTER(cr)) {
    for(i=0; (i<nr); i++)
      handle[i] = lookup_symtab(symtab,resinfo[i].name);
  }
  nblock_bc(cr,nr,handle);

  if (!MASTER(cr)) {
    for (i=0; (i<nr); i++) 
      resinfo[i].name = get_symtab_handle(symtab,handle[i]);
  }
  sfree(handle);
*/
}

static void generic_pack_symtab(gmx_packed_t packed, ePackType packtype, t_symtab *symtab, int *buf_lengths)
{
    int i, len;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    for (i = 0; i < symtab->nr; i++)
    {
        snew_pack(packtype, mtop->symtab->symbuf->buf[i], buf_lengths[i]);
        generic_pack_function(packed, symtab->symbuf->buf[i], buf_lengths[i], MPI_CHAR);
    }
/*
    block_bc(cr,symtab->nr);
    snew_pack(cr,symtab->symbuf,1);
    symtab->symbuf->bufsize = symtab->nr;
    snew_pack(cr,symtab->symbuf->buf,symtab->nr);
    for (i=0; i<nr; i++)
    {
        if (MASTER(cr))
        {
            len = strlen(symtab->symbuf->buf[i]) + 1;
        }
        block_bc(cr,len);
        snew_pack(cr,symtab->symbuf->buf[i],len);
        nblock_bc(cr,len,symtab->symbuf->buf[i]);
    }
*/
}

static void generic_pack_block(gmx_packed_t packed, ePackType packtype, t_block *block)
{
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    snew_pack(packtype, block->index, block->nr+1);
    generic_pack_function(packed, block->index, block->nr+1, MPI_INT);
/*
    block_bc(cr,block->nr);
    snew_pack(cr,block->index,block->nr+1);
    nblock_bc(cr,block->nr+1,block->index);
*/
}

static void generic_pack_blocka(gmx_packed_t packed, ePackType packtype, t_blocka *block)
{
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    snew_pack(packtype, block->index, block->nr+1);
    generic_pack_function(packed, block->index, block->nr+1, MPI_INT);
    if (block->nra)
    {
        snew_pack(packtype, block->a, block->nra);
        generic_pack_function(packed, block->a, block->nra, MPI_INT);
    }
/*
    block_bc(cr,block->nr);
    snew_pack(cr,block->index,block->nr+1);
    nblock_bc(cr,block->nr+1,block->index);
    block_bc(cr,block->nra);
    if (block->nra)
    {
        snew_pack(cr,block->a,block->nra);
        nblock_bc(cr,block->nra,block->a);
    }
*/
}

static void generic_pack_grps(gmx_packed_t packed, ePackType packtype, t_grps grps[])
{
    int i;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);
  
    for(i = 0; i < egcNR; i++)
    {
        snew_pack(packtype, grps[i].nm_ind, grps[i].nr);
        generic_pack_function(packed, grps[i].nm_ind, grps[i].nr, MPI_INT);
    }
/*
    for(i=0; (i<egcNR); i++)
    {
        block_bc(cr,grps[i].nr);
        snew_pack(cr,grps[i].nm_ind,grps[i].nr);
        nblock_bc(cr,grps[i].nr,grps[i].nm_ind);
    }
*/
}

static void generic_pack_atoms(gmx_packed_t packed, ePackType packtype, t_symtab *symtab, t_atoms *atoms)
{
    int dummy;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    snew_pack(packtype, atoms->atom, atoms->nr);
    generic_pack_function(packed, atoms->atom, sizeof(atoms->atom[0]) * atoms->nr, MPI_BYTE);
    generic_pack_strings(cr, symtab, atoms->nr, &atoms->atomname);
    snew_pack(packtype, atoms->resinfo, atoms->nres);
    generic_pack_function(packed, atoms->resinfo, sizeof(atoms->resinfo[0]) * atoms->nres, MPI_BYTE);
    generic_pack_strings_resinfo(cr, symtab, atoms->nres, atoms->resinfo);
    /* QMMM requires atomtypes to be known on all nodes as well */
    generic_pack_strings(cr, symtab, atoms->nr, &atoms->atomtype);
    generic_pack_strings(cr, symtab, atoms->nr, &atoms->atomtypeB);
/*
    block_bc(cr,atoms->nr);
    snew_pack(cr,atoms->atom,atoms->nr);
    nblock_bc(cr,atoms->nr,atoms->atom);
    bc_strings(cr,symtab,atoms->nr,&atoms->atomname);
    block_bc(cr,atoms->nres);
    snew_pack(cr,atoms->resinfo,atoms->nres);
    nblock_bc(cr,atoms->nres,atoms->resinfo);
    bc_strings_resinfo(cr,symtab,atoms->nres,atoms->resinfo);
    /* QMMM requires atomtypes to be known on all nodes as well * /
    bc_strings(cr,symtab,atoms->nr,&atoms->atomtype);
    bc_strings(cr,symtab,atoms->nr,&atoms->atomtypeB);
*/
}

static void generic_pack_groups(gmx_packed_t packed, ePackType packtype, t_symtab *symtab,
                      int natoms, gmx_groups_t *groups, int *grpnr_sizes)
{
    int dummy;
    int g, n;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    generic_pack_grps(packed, packtype, groups->grps);
    generic_pack_strings(packed, packtype, symtab, groups->ngrpname, &groups->grpname);
    for(g = 0; g < egcNR; g++)
    {
        if (0 == grpnr_sizes[g])
        {
            mtop->groups->grpnr[g] = NULL;
        }
        else
        {
            snew_pack(packtype, groups->grpnr[g], grpnr_sizes[g]);
            generic_pack_function(packed, groups->grpnr[g], grpnr_sizes[g], MPI_CHAR);
        }
    }
    if (debug) fprintf(debug,"after generic_pack_groups\n");
/*
    bc_grps(cr,groups->grps);
    block_bc(cr,groups->ngrpname);
    bc_strings(cr,symtab,groups->ngrpname,&groups->grpname);
    for(g=0; g<egcNR; g++)
    {
        if (MASTER(cr))
        {
            if (groups->grpnr[g])
            {
                n = natoms;
            }
            else
            {
                n = 0;
            }
        }
        block_bc(cr,n);
        if (n == 0)
        {
            groups->grpnr[g] = NULL;
        }
        else
        {
            snew_pack(cr,groups->grpnr[g],n);
            nblock_bc(cr,n,groups->grpnr[g]);
        }
    }
    if (debug) fprintf(debug,"after bc_groups\n");
*/
}

static void generic_pack_ilists(gmx_packed_t packed, ePackType packtype, t_ilist *ilist)
{
    int ftype;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    for(ftype = 0; ftype < F_NRE; ftype++)
    {
        snew_pack(packtype, ilist[ftype].iatoms, ilist[ftype].nr);
        generic_pack_function(packed, ilist[ftype].iatoms, ilist[ftype].nr, MPI_INT);
    }

    if (debug) fprintf(debug,"after generic_pack_ilists\n");
/*
    /* Here we only communicate the non-zero length ilists * /
    if (MASTER(cr))
    {
        for(ftype=0; ftype<F_NRE; ftype++)
        {
            if (ilist[ftype].nr > 0)
            {
                block_bc(cr,ftype); //omit
                block_bc(cr,ilist[ftype].nr);
                nblock_bc(cr,ilist[ftype].nr,ilist[ftype].iatoms);
            }
        }
        ftype = -1;
        block_bc(cr,ftype);
    } else {
        for(ftype=0; ftype<F_NRE; ftype++)
        {
            ilist[ftype].nr = 0;
        }
        do {
            block_bc(cr,ftype);//omit
            if (ftype >= 0)
            {
                block_bc(cr,ilist[ftype].nr);
                snew_pack(cr,ilist[ftype].iatoms,ilist[ftype].nr);
                nblock_bc(cr,ilist[ftype].nr,ilist[ftype].iatoms);
            }
        } while (ftype >= 0);
    }

    if (debug) fprintf(debug,"after bc_ilists\n");
*/
}

/*
static void bc_idef(gmx_packed_t packed, ept packtype, t_idef *idef)
{
    block_bc(cr,idef->ntypes); //send early
    block_bc(cr,idef->atnr);
    snew_pack(cr,idef->functype,idef->ntypes);
    snew_pack(cr,idef->iparams,idef->ntypes);
    nblock_bc(cr,idef->ntypes,idef->functype);
    nblock_bc(cr,idef->ntypes,idef->iparams);
    block_bc(cr,idef->fudgeQQ);
    bc_ilists(cr,idef->il);
    block_bc(cr,idef->ilsort);
}
*/

static void generic_pack_cmap(gmx_packed_t packed, ePackType packtype, gmx_cmap_t *cmap_grid)
{
    int i, j,nelem, ngrid;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);
	
    ngrid = cmap_grid->ngrid;
    nelem = cmap_grid->grid_spacing * cmap_grid->grid_spacing;
	
    if(ngrid > 0)
    {
        snew_pack(packtype, cmap_grid->cmapdata, ngrid);
		
        for(i = 0; i < ngrid; i++)
        {
            snew_pack(packtype, cmap_grid->cmapdata[i].cmap, 4*nelem);
            generic_pack_function(packed, cmap_grid->cmapdata[i].cmap, 4*nelem, GMX_MPI_REAL);
        }
    }
/*
    block_bc(cr,cmap_grid->ngrid);//send early DONE
    block_bc(cr,cmap_grid->grid_spacing);//send early DONE
	
    ngrid = cmap_grid->ngrid;
    nelem = cmap_grid->grid_spacing * cmap_grid->grid_spacing;
	
    if(ngrid>0)
    {
        snew_pack(cr,cmap_grid->cmapdata,ngrid);
		
        for(i=0;i<ngrid;i++)
        {
            snew_pack(cr,cmap_grid->cmapdata[i].cmap,4*nelem);
            nblock_bc(cr,4*nelem,cmap_grid->cmapdata[i].cmap);
        }
    }
*/
}

static void generic_pack_ffparams(gmx_packed_t packed, ePackType packtype, gmx_ffparams_t *ffp)
{
    int i;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);
  
    generic_pack_function(packed, ffp->atnr, 1, MPI_INT);
    snew_pack(packtype, ffp->functype, ffp->ntypes);
    snew_pack(packtype, ffp->iparams, ffp->ntypes);
    generic_pack_function(packed, ffp->functype, ffp->ntypes, MPI_INT);
    generic_pack_function(packed, ffp->iparams, sizeof(ffp->iparams[0]) * ffp->ntypes, MPI_BYTE);
    generic_pack_function(packed, ffp->reppow, 1, MPI_DOUBLE);
    generic_pack_function(packed, ffp->fudgeQQ, 1, GMX_MPI_REAL);
    generic_pack_cmap(packed, packtype, &ffp->cmap_grid);
/*
    block_bc(cr,ffp->ntypes);
    block_bc(cr,ffp->atnr);
    snew_pack(cr,ffp->functype,ffp->ntypes);
    snew_pack(cr,ffp->iparams,ffp->ntypes);
    nblock_bc(cr,ffp->ntypes,ffp->functype);
    nblock_bc(cr,ffp->ntypes,ffp->iparams);
    block_bc(cr,ffp->reppow);
    block_bc(cr,ffp->fudgeQQ);
    bc_cmap(cr,&ffp->cmap_grid);
*/
}

static void generic_pack_grpopts(gmx_packed_t packed, ePackType packtype, t_grpopts *g)
{
    int i, n;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    snew_pack(packtype, g->nrdf, g->ngtc);
    snew_pack(packtype, g->tau_t, g->ngtc);
    snew_pack(packtype, g->ref_t, g->ngtc);
    snew_pack(packtype, g->acc, g->ngacc);
    snew_pack(packtype, g->nFreeze, g->ngfrz);
    snew_pack(packtype, g->egp_flags, g->ngener*g->ngener);
        
    generic_pack_function(packed, g->nrdf, g->ngtc, GMX_MPI_REAL);
    generic_pack_function(packed, g->tau_t, g->ngtc, GMX_MPI_REAL);
    generic_pack_function(packed, g->ref_t, g->ngtc, GMX_MPI_REAL);
    generic_pack_function(packed, g->acc, g->ngacc * DIM, GMX_MPI_REAL);
    generic_pack_function(packed, g->nFreeze, g->ngfrz * DIM, MPI_INT);
    generic_pack_function(packed, g->egp_flags, g->ngener*g->ngener, MPI_INT);
    snew_pack(packtype, g->annealing, g->ngtc);
    snew_pack(packtype, g->anneal_npoints, g->ngtc);
    snew_pack(packtype, g->anneal_time, g->ngtc);
    snew_pack(packtype, g->anneal_temp, g->ngtc);
    generic_pack_function(packed, g->annealing, g->ngtc, MPI_INT);
    for(i = 0;i < g->ngtc; i++)
    {
        n = g->anneal_npoints[i];
        if (n > 0)
        {
            snew_pack(packtype, g->anneal_time[i], n);
            snew_pack(packtype, g->anneal_temp[i], n);
            generic_pack_function(packed, g->anneal_time[i], n, GMX_MPI_REAL);
            generic_pack_function(packed, g->anneal_temp[i], n, GMX_MPI_REAL);
        }
    }
        
    /* QMMM stuff, see inputrec */
    snew_pack(packtype, g->QMmethod, g->ngQM);
    snew_pack(packtype, g->QMbasis, g->ngQM);
    snew_pack(packtype, g->QMcharge, g->ngQM);
    snew_pack(packtype, g->QMmult, g->ngQM);
    snew_pack(packtype, g->bSH, g->ngQM);
    snew_pack(packtype, g->CASorbitals, g->ngQM);
    snew_pack(packtype, g->CASelectrons, g->ngQM);
    snew_pack(packtype, g->SAon, g->ngQM);
    snew_pack(packtype, g->SAoff, g->ngQM);
    snew_pack(packtype, g->SAsteps, g->ngQM);
    
    if (g->ngQM)
    {
        generic_pack_function(packed, g->QMmethod, g->ngQM, MPI_INT);
        generic_pack_function(packed, g->QMbasis, g->ngQM, MPI_INT);
        generic_pack_function(packed, g->QMcharge, g->ngQM, MPI_INT);
        generic_pack_function(packed, g->QMmult, g->ngQM, MPI_INT);
        generic_pack_function(packed, g->bSH, g->ngQM, MPI_INT);
        generic_pack_function(packed, g->CASorbitals, g->ngQM, MPI_INT);
        generic_pack_function(packed, g->CASelectrons, g->ngQM, MPI_INT);
        generic_pack_function(packed, g->SAon, g->ngQM, GMX_MPI_REAL);
        generic_pack_function(packed, g->SAoff, g->ngQM, GMX_MPI_REAL);
        generic_pack_function(packed, g->SAsteps, g->ngQM, MPI_INT);
        /* end of QMMM stuff */
    }
/*
    block_bc(cr,g->ngtc);//send early DONE
    block_bc(cr,g->ngacc);//send early DONE
    block_bc(cr,g->ngfrz);//send early DONE
    block_bc(cr,g->ngener);//send early DONE
    snew_pack(cr,g->nrdf,g->ngtc);
    snew_pack(cr,g->tau_t,g->ngtc);
    snew_pack(cr,g->ref_t,g->ngtc);
    snew_pack(cr,g->acc,g->ngacc);
    snew_pack(cr,g->nFreeze,g->ngfrz);
    snew_pack(cr,g->egp_flags,g->ngener*g->ngener);
    
    nblock_bc(cr,g->ngtc,g->nrdf);
    nblock_bc(cr,g->ngtc,g->tau_t);
    nblock_bc(cr,g->ngtc,g->ref_t);
    nblock_bc(cr,g->ngacc,g->acc);
    nblock_bc(cr,g->ngfrz,g->nFreeze);
    nblock_bc(cr,g->ngener*g->ngener,g->egp_flags);
    snew_pack(cr,g->annealing,g->ngtc);
    snew_pack(cr,g->anneal_npoints,g->ngtc);
    snew_pack(cr,g->anneal_time,g->ngtc);
    snew_pack(cr,g->anneal_temp,g->ngtc);
    nblock_bc(cr,g->ngtc,g->annealing);
    nblock_bc(cr,g->ngtc,g->anneal_npoints); //send not quite so early DONE
    for(i=0;(i<g->ngtc); i++)
    {
        n = g->anneal_npoints[i];
        if (n > 0)
        {
            snew_pack(cr,g->anneal_time[i],n);
            snew_pack(cr,g->anneal_temp[i],n);
            nblock_bc(cr,n,g->anneal_time[i]);
            nblock_bc(cr,n,g->anneal_temp[i]);
        }
    }
    
    /* QMMM stuff, see inputrec * /
    block_bc(cr,g->ngQM);//send early DONE
    snew_pack(cr,g->QMmethod,g->ngQM);
    snew_pack(cr,g->QMbasis,g->ngQM);
    snew_pack(cr,g->QMcharge,g->ngQM);
    snew_pack(cr,g->QMmult,g->ngQM);
    snew_pack(cr,g->bSH,g->ngQM);
    snew_pack(cr,g->CASorbitals,g->ngQM);
    snew_pack(cr,g->CASelectrons,g->ngQM);
    snew_pack(cr,g->SAon,g->ngQM);
    snew_pack(cr,g->SAoff,g->ngQM);
    snew_pack(cr,g->SAsteps,g->ngQM);
    
    if (g->ngQM)
    {
        nblock_bc(cr,g->ngQM,g->QMmethod);
        nblock_bc(cr,g->ngQM,g->QMbasis);
        nblock_bc(cr,g->ngQM,g->QMcharge);
        nblock_bc(cr,g->ngQM,g->QMmult);
        nblock_bc(cr,g->ngQM,g->bSH);
        nblock_bc(cr,g->ngQM,g->CASorbitals);
        nblock_bc(cr,g->ngQM,g->CASelectrons);
        nblock_bc(cr,g->ngQM,g->SAon);
        nblock_bc(cr,g->ngQM,g->SAoff);
        nblock_bc(cr,g->ngQM,g->SAsteps);
        /* end of QMMM stuff * /
    }
*/
}

static void generic_pack_cosines(gmx_packed_t packed, ePackType packtype, t_cosines *cs)
{
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    snew_pack(packtype, cs->a, cs->n);
    snew_pack(packtype, cs->phi, cs->n);
    if (cs->n > 0)
    {
        generic_pack_function(packed, cs->a, cs->n);
        generic_pack_function(packed, cs->phi, cs->n);
    }
/*
    block_bc(cr,cs->n);
    snew_pack(cr,cs->a,cs->n);
    snew_pack(cr,cs->phi,cs->n);
    if (cs->n > 0)
    {
        nblock_bc(cr,cs->n,cs->a);
        nblock_bc(cr,cs->n,cs->phi);
    }
*/
}

static void generic_pack_pullgrp(gmx_packed_t packed, ePackType packtype, t_pullgrp *pgrp)
{
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    if (pgrp->nat > 0)
    {
        snew_pack(packtype, pgrp->ind, pgrp->nat);
        generic_pack_function(packed, pgrp->ind, pgrp->nat, MPI_INT);
    }
    if (pgrp->nweight > 0)
    {
        snew_pack(packtype, pgrp->weight, pgrp->nweight);
        generic_pack_function(packed, pgrp->weight, pgrp->nweight, GMX_MPI_REAL);
    }
/*
    block_bc(cr,*pgrp);
    if (pgrp->nat > 0)
    {
        snew_pack(cr,pgrp->ind,pgrp->nat);
        nblock_bc(cr,pgrp->nat,pgrp->ind);
    }
    if (pgrp->nweight > 0)
    {
        snew_pack(cr,pgrp->weight,pgrp->nweight);
        nblock_bc(cr,pgrp->nweight,pgrp->weight);
    }
*/
}

static void generic_pack_pull(gmx_packed_t packed, ePackType packtype, t_pull *pull)
{
    int g;

    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    snew_pack(packtype, pull->grp, pull->ngrp+1);
    for(g = 0; g < pull->ngrp+1; g++)
    {
        generic_pack_pullgrp(cr, &pull->grp[g]);
    }
/*
    block_bc(cr,*pull);
    snew_pack(cr,pull->grp,pull->ngrp+1);
    for(g=0; g<pull->ngrp+1; g++)
    {
        bc_pullgrp(cr,&pull->grp[g]);
    }
*/
}

static void generic_pack_inputrec(gmx_packed_t packed, ePackType packtype, t_inputrec *inputrec)
{
    gmx_bool bAlloc = TRUE;
    int i;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    snew_pack(packtype, inputrec->flambda, inputrec->n_flambda);
    generic_pack_function(packed, inputrec->flambda, inputrec->n_flambda, MPI_DOUBLE);

    generic_pack_grpopts(packed, packtype, &(inputrec->opts));
    if (inputrec->ePull != epullNO)
    {
        snew_pack(packtype, inputrec->pull, 1);
        generic_pack_pull(packed, packtype, inputrec->pull);
    }
    for(i = 0; i < DIM; i++)
    {
        generic_pack_cosines(packed, packtype, &(inputrec->ex[i]));
        generic_pack_cosines(packed, packtype, &(inputrec->et[i]));
    }
/*
    block_bc(cr,*inputrec);
    snew_pack(cr,inputrec->flambda,inputrec->n_flambda);
    nblock_bc(cr,inputrec->n_flambda,inputrec->flambda);
    bc_grpopts(cr,&(inputrec->opts));
    if (inputrec->ePull != epullNO)
    {
        snew_pack(cr,inputrec->pull,1);
        bc_pull(cr,inputrec->pull);
    }
    for(i=0; (i<DIM); i++)
    {
        bc_cosines(cr,&(inputrec->ex[i]));
        bc_cosines(cr,&(inputrec->et[i]));
    }
*/
}

static void generic_pack_moltype(gmx_packed_t packed, ePackType packtype, t_symtab *symtab,
                       gmx_moltype_t *moltype)
{
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    generic_pack_string(cr, symtab, &moltype->name);
    generic_pack_atoms(cr, symtab, &moltype->atoms); // stuff inside not quite so early
    if (debug) fprintf(debug,"after generic_pack_atoms\n");

    generic_pack_ilists(cr, moltype->ilist); //stuff inside not quite so early
    generic_pack_block(cr, &moltype->cgs); //stuff inside not quite so early DONE
    generic_pack_blocka(cr, &moltype->excls);
}

static void generic_pack_molblock(gmx_packed_t packed, ePackType packtype, gmx_molblock_t *molb)
{
    gmx_bool bAlloc = TRUE;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);
  
    generic_pack_function(packed, molb->type, 1, MPI_INT);
    generic_pack_function(packed, molb->nmol, 1, MPI_INT);
    generic_pack_function(packed, molb->natoms_mol, 1, MPI_INT);
    if (molb->nposres_xA > 0)
    {
        snew_pack(packtype, molb->posres_xA, molb->nposres_xA);
        generic_pack_function(packed, molb->posres_xA[0], molb->nposres_xA * DIM, GMX_MPI_REAL);
    }
    if (molb->nposres_xB > 0)
    {
        snew_pack(packtype, molb->posres_xB, molb->nposres_xB);
        generic_pack_function(packed, molb->posres_xB[0], molb->nposres_xB * DIM, GMX_MPI_REAL);
    }
    if (debug) fprintf(debug,"after generic_pack_molblock\n");
/*
    block_bc(cr,molb->type);
    block_bc(cr,molb->nmol);
    block_bc(cr,molb->natoms_mol);
    block_bc(cr,molb->nposres_xA);
    if (molb->nposres_xA > 0)
    {
        snew_pack(cr,molb->posres_xA,molb->nposres_xA);
        nblock_bc(cr,molb->nposres_xA*DIM,molb->posres_xA[0]);
    }
    block_bc(cr,molb->nposres_xB);
    if (molb->nposres_xB > 0)
    {
        snew_pack(cr,molb->posres_xB,molb->nposres_xB);
        nblock_bc(cr,molb->nposres_xB*DIM,molb->posres_xB[0]);
    }
    if (debug) fprintf(debug,"after bc_molblock\n");
*/
}

static void generic_pack_atomtypes(gmx_packed_t packed, ePackType packtype, t_atomtypes *atomtypes)
{
    int nr;
    t_pack_function generic_pack_function = (eptPackSize == packtype ) ? gmx_pack_size : (eptPack == packtype ? gmx_pack : gmx_unpack);

    nr = atomtypes->nr;

    snew_pack(packtype, atomtypes->radius, nr);
    snew_pack(packtype, atomtypes->vol, nr);
    snew_pack(packtype, atomtypes->surftens, nr);
    snew_pack(packtype, atomtypes->gb_radius, nr);
    snew_pack(packtype, atomtypes->S_hct, nr);

    generic_pack_function(packed, atomtypes->radius, nr, GMX_MPI_REAL);
    generic_pack_function(packed, atomtypes->vol, nr, GMX_MPI_REAL);
    generic_pack_function(packed, atomtypes->surftens, nr, GMX_MPI_REAL);
    generic_pack_function(packed, atomtypes->gb_radius, nr, GMX_MPI_REAL);
    generic_pack_function(packed, atomtypes->S_hct, nr, GMX_MPI_REAL);
/*
    block_bc(cr,atomtypes->nr);

    nr = atomtypes->nr;

    snew_pack(cr,atomtypes->radius,nr);
    snew_pack(cr,atomtypes->vol,nr);
    snew_pack(cr,atomtypes->surftens,nr);
    snew_pack(cr,atomtypes->gb_radius,nr);
    snew_pack(cr,atomtypes->S_hct,nr);

    nblock_bc(cr,nr,atomtypes->radius);
    nblock_bc(cr,nr,atomtypes->vol);
    nblock_bc(cr,nr,atomtypes->surftens);
    nblock_bc(cr,nr,atomtypes->gb_radius);
    nblock_bc(cr,nr,atomtypes->S_hct);
*/
}

/* Broadcast a t_inputrec and gmx_mtop_t to the other processors in
 * the simulation. */
void bcast_ir_mtop(const t_commrec *cr, t_inputrec *inputrec, gmx_mtop_t *mtop)
{
    int i, g;
    int *buf_lengths, *grpnr_sizes;
    gmx_packed_t packed = gmx_packed_init(cr);
    t_pack_function generic_pack_function = MASTER(cr) ? gmx_pack : gmx_unpack;
    ePackType packtype;

    /* First, send around various integers that determine sizes of
     * arrays, etc. */

    /* How much will we pack? This has to be kept in sync with the
     * function lines below. */
    gmx_pack_size(packed, sizeof(*inputrec), MPI_BYTE);
    gmx_pack_size(packed, 5, MPI_INT);
    gmx_pack_size(packed, sizeof(*pull), MPI_BYTE);
    gmx_pack_size(packed, 2 * DIM + 10 + 2 * egcNR, MPI_INT);

    /* Allocate */
    gmx_packed_allocate(packed);

    if (MASTER(cr))
    {
        packtype = bcPack;
    }
    else
    {
        packtype = bcUnpack;
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Do packing or unpacking, as appropriate */
    generic_pack_function(packed, *inputrec, sizeof(*inputrec), MPI_BYTE);
    generic_pack_function(packed, inputrec->opts->ngtc, 1, MPI_INT);
    generic_pack_function(packed, inputrec->opts->ngacc, 1, MPI_INT);
    generic_pack_function(packed, inputrec->opts->ngfrz, 1, MPI_INT);
    generic_pack_function(packed, inputrec->opts->ngener, 1, MPI_INT);
    generic_pack_function(packed, inputrec->opts->ngQM, 1, MPI_INT);
    generic_pack_function(packed, *inputrec->pull, sizeof(*inputrec->pull), MPI_BYTE);
    for (i = 0; i < DIM; i++)
    {
        generic_pack_function(packed, inputrec->ex[i]->n, 1, MPI_INT);
        generic_pack_function(packed, inputrec->et[i]->n, 1, MPI_INT);
    }
    generic_pack_function(packed, mtop->symtab->nr, 1, MPI_INT);
    generic_pack_function(packed, mtop->ffparams->ntypes, 1, MPI_INT);
    generic_pack_function(packed, mtop->ffparams->cmap->ngrid, 1, MPI_INT);
    generic_pack_function(packed, mtop->ffparams->cmap->ngrid_spacing, 1, MPI_INT);
    generic_pack_function(packed, mtop->nmoltype, 1, MPI_INT);
    generic_pack_function(packed, mtop->nmolblock, 1, MPI_INT);
    generic_pack_function(packed, mtop->natoms, 1, MPI_INT);
    generic_pack_function(packed, mtop->atomtypes->nr, 1, MPI_INT);
    generic_pack_function(packed, mtop->mols->nr, 1, MPI_INT);
    for(g = 0; g < egcNR; g++)
    {
        generic_pack_function(packed, mtop->groups->grps[g].nr, 1, MPI_INT);
    }
    generic_pack_function(packed, mtop->groups->ngrpname, 1, MPI_INT);
    snew(grpnr_sizes, egcNR);
    for(g = 0; g < egcNR; g++)
    {
        grpnr_sizes[j] = mtop->groups->grpnr[g] ? mtop->natoms : 0;
    }
    generic_pack_function(packed, grpnr_sizes, egcNR, MPI_INT);

    if (MASTER(cr))
    {
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Clean up, prepare for next round of packing */
    gmx_packed_refresh(packed);
    packtype = bcPackSize;

    /* Now do some allocation for the next stage of data*/
    snew_pack(packtype, inputrec->opts->anneal_npoints, g->ngtc);
    snew_pack(packtype, inputrec->pull->grp, inputrec->pull->ngrp+1);
    snew_pack(packtype, mtop->symtab->symbuf, 1);
    mtop->symtab->symbuf->bufsize = mtop->symtab->nr;
    snew_pack(packtype, mtop->symtab->symbuf->buf, symtab->nr);
    snew_pack(packtype, mtop->moltype, mtop->nmoltype);
    snew_pack(packtype, mtop->molblock, mtop->nmolblock);

    /* Now send the next stage of the mtop data. How much will we
     * pack? This has to be kept in sync with the function lines
     * below. */
    gmx_pack_size(packed, inputrec->opts->ngtc, MPI_INT);
    gmx_pack_size(packed, sizeof(inputrec->pull->grp[0]) * inputrec->pull->ngrp+1, MPI_BYTE);
    gmx_pack_size(packed, mtop->symtab->nr, MPI_INT);
    gmx_pack_size(packed, 4 * mtop->nmoltype, MPI_INT);
    gmx_pack_size(packed, 2 * mtop->nmolblock, MPI_INT);

    /* Allocate */
    gmx_packed_allocate(packed);

    if (!MASTER(cr))
    {
        packtype = bcPack;
    }
    else
    {
        packtype = bcUnpack;
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Do packing or unpacking, as appropriate */
    generic_pack_function(packed, inputrec->opts->anneal_npoints, inputrec->opts->ngtc, MPI_INT);
    for (g = 0; g < inputrec->pull->ngrp+1; g++)
    {
        generic_pack_function(packed, &inputrec->pull->grp[g], sizeof(inputrec->pull->grp[g]), MPI_BYTE);
    }
    for (i = 0; i < mtop->symtab->nr; i++)
    {
        if (MASTER(cr))
        {
            buf_lengths[i] = strlen(mtop->symtab->symbuf->buf[i]) + 1;
        }
    }
    generic_pack_function(packed, buf_lengths, mtop->symtab->nr, MPI_INT);
    for (i = 0; i < mtop->nmoltype; i++)
    {
        generic_pack_function(packed, mtop->moltype[i]->atoms.nr, 1, MPI_INT);
        generic_pack_function(packed, mtop->moltype[i]->atoms.nres, 1, MPI_INT);
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            generic_pack_function(packed, mtop->moltype[i]->ilist[ftype].nr, 1, MPI_INT);
        }
        generic_pack_function(packed, mtop->moltype[i]->cgs.nr, 1, MPI_INT);
        generic_pack_function(packed, mtop->moltype[i]->cgs.nra, 1, MPI_INT);
    }
    for(i = 0; i < mtop->nmolblock; i++)
    {
        generic_pack_function(packed, mtop->molblock[i]->nposres_xA, 1, MPI_INT);
        generic_pack_function(packed, mtop->molblock[i]->nposres_xB, 1, MPI_INT);
    }

    if (MASTER(cr))
    {
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Clean up, prepare for final round of packing */
    gmx_packed_refresh(packed);
    packtype = bcPackSize;

    generic_pack_inputrec(packed, packtype, eptPackSize, inputrec);
    generic_pack_symtab(packed, packtype, eptPackSize, &mtop->symtab, buf_lengths);
    generic_pack_string(packed, packtype, eptPackSize, &mtop->symtab, &mtop->name);
    generic_pack_ffparams(packed, packtype, eptPackSize, &mtop->ffparams);
    for(i = 0; i < mtop->nmoltype; i++)
    {
        generic_pack_moltype(packed, packtype, eptPackSize, &mtop->symtab, &mtop->moltype[i]);
    }
    for(i = 0; i < mtop->nmolblock; i++)
    {
        generic_pack_molblock(packed, packtype, eptPackSize, &mtop->molblock[i]);
    }
    block_bc(packed, packtype, eptPackSize, mtop->natoms);
    generic_pack_atomtypes(packed, packtype, eptPackSize, &mtop->atomtypes);
    generic_pack_block(packed, packtype, eptPackSize, &mtop->mols);
    generic_pack_groups(packed, packtype, eptPackSize, &mtop->symtab, mtop->natoms, &mtop->groups);

    /* Allocate */
    gmx_packed_allocate(packed);

    if (!MASTER(cr))
    {
        packtype = bcPack;
    }
    else
    {
        packtype = bcUnpack;
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Do packing or unpacking, as appropriate */
    generic_pack_inputrec(packed, packtype, inputrec);
    generic_pack_symtab(packed, packtype, &mtop->symtab, buf_lengths);
    generic_pack_string(packed, packtype, &mtop->symtab, &mtop->name);

    generic_pack_ffparams(packed, packtype, &mtop->ffparams);

    for(i = 0; i < mtop->nmoltype; i++)
    {
        generic_pack_moltype(packed, packtype, &mtop->symtab, &mtop->moltype[i]);
    }

    for(i = 0; i < mtop->nmolblock; i++)
    {
        generic_pack_molblock(packed, packtype, &mtop->molblock[i]);
    }

    generic_pack_atomtypes(packed, packtype, &mtop->atomtypes);

    generic_pack_block(packed, packtype, &mtop->mols);
    generic_pack_groups(packed, packtype, &mtop->symtab, mtop->natoms, &mtop->groups);

    if (MASTER(cr))
    {
        /* communicate */
        do_packed_bcast(packed);
    }

    /* Clean up */
    gmx_packed_destroy(packed);
    sfree(buf_lengths);
    sfree(grpnr_sizes);
}
#endif
