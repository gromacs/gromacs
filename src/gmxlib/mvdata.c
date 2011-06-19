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

/*! \brief Code for packed MPI communication of t_state, t_inputrec
 *         and gmx_mtop_t variables.
 *
 * This code uses proper MPI packing and unpacking to communicate
 * dynamically-allocated data structures with minimal calls to
 * functions that actually communicate between processors.
 *
 * It does so by abstracting the necessary calls to MPI_Pack_size,
 * MPI_Pack and MPI_Unpack through a common umbrella function. This
 * umbrella function is called from a callback function that is unique
 * to each kind of data. The necessary house-keeping is performed
 * behind the scenes. The callback is first invoked in a mode that
 * determines the (upper bound of the) size of the packed data
 * structure, and memory is allocated accordingly on all nodes. The
 * sending (master) node then invokes the callback function again in
 * "pack" mode. All nodes execute the MPI communication call. Then the
 * receiving (non-master) nodes then invoke the callback function
 * again in "unpack" mode. This way, only the callback function needs
 * to contain code that maps normal data types to MPI_Datatypes,
 * rather than in three places in a simple MPI_Pack implementation.
 *
 * Clearly, dynamically allocated arrays of dynamically allocated
 * arrays (etc.) present a difficulty, because all nodes need to know
 * in advance how much data they will communicate. So multiple phases
 * of packing and unpacking can be required. The packaging of the
 * house-keeping makes this very simple. All that is required is
 * multiple callback functions. An example exists in bcast_ir_mtop().
 *
 * Large arrays (e.g. of size natoms) don't benefit from reducing
 * calls to the communication library by doing such packing, anda
 * developer can feel free to include a normal MPI communcation call
 * in a callback function, so long as they take care to see that it is
 * called suitably from the different kinds of nodes. An example
 * exists in bcast_state().
 *
 * The usual trick of defining a struct to contain the necessary data
 * that must be preserved between calls to the callback function can
 * be used, via the usual void pointer argument to the callback and
 * the house-keeping function.
 *
 * The 4.5.x (and earlier) code used dozens of calls to MPI_Bcast to
 * deal with the fact that propagating dynamically allocated structs
 * is awkward in MPI. Since this process only really occured during
 * system setup and replica exchange, it was not a big deal. However,
 * the nature of the solution meant that
 *
 */

#ifndef GMX_MPI
/* Calling any of these entry points is illegal without MPI, but the
 * non-MPI code contains calls to them them, so they have to be
 * defined and able to be compiled. */
void bcast_state_setup(const t_commrec *cr, t_state *state)
{
    gmx_call("bcast_state_setup");
}

void bcast_state(const t_commrec *cr, t_state *state, gmx_bool bAlloc)
{
    gmx_call("bcast_state");
}

void bcast_ir_mtop(const t_commrec *cr, t_inputrec *inputrec, gmx_mtop_t *mtop)
{
    gmx_call("bcast_ir_mtop");
}

#else

/*! Determines the mode in which a function is invoked during MPI
 *  packing. */
typedef enum {eptPackSize, eptPack, eptUnpack, eptPackTypeNR} ePackType;

/*! /brief Manage packing and unpacking transparently.
 */
typedef struct
{
    const t_commrec *cr;
    ePackType packtype;
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
    packed->packtype = eptPackSize;
    return packed;
}

void gmx_packed_allocate(gmx_packed_t packed)
{
    snew(packed->buffer, packed->size_of_buffer);
    packed->position = 0;
}

void gmx_packed_destroy(gmx_packed_t packed)
{
    sfree(packed->buffer);
    sfree(packed);
}

/*! /brief Do the house-keeping around the calls to the callback
 *         function.
 *
 * It would be nice if this function could allocate the data pointer
 * upon request (which is wanted by bcast_state()), but without proper
 * function overloading, that's possible only by writing ugly code.
 */
void gmx_pack_function(gmx_packed_t packed, void *data, int nr, MPI_Datatype mpitype)
{
    int packed_size = 0;
    switch(packed->packtype)
    {
    case eptPackSize:
        MPI_Pack_size(nr, mpitype, packed->cr->mpi_comm_mygroup, &packed_size);
        packed->size_of_buffer += packed_size;
        break;
    case eptPack:
        if (0 < nr)
        {
            MPI_Pack(data, nr, mpitype, packed->buffer, packed->size_of_buffer, &packed->position, packed->cr->mpi_comm_mygroup);
        }
        break;
    case eptUnpack:
        if (0 < nr)
        {
            MPI_Unpack(packed->buffer, packed->size_of_buffer, &packed->position, data, nr, mpitype, packed->cr->mpi_comm_mygroup);
        }
        break;
    }
}

/*! /brief Actually do the broadcast of the packed data.
 */
void do_packed_bcast(gmx_packed_t packed)
{
    /* Probably the test for 0 < packed->position is only needed on
     * BlueGene(/L), where IBM's MPI_Bcast will segfault after
     * dereferencing a null pointer, even when no data is to be
     * transferred. */
    if (0 < packed->size_of_buffer)
    {
        MPI_Bcast(packed->buffer, packed->size_of_buffer, MPI_PACKED, MASTERRANK(packed->cr), packed->cr->mpi_comm_mygroup);
    }
}

/*! /brief Helper allocation macro
 *
 * The idea is that we want to be able allocate sub-structures before
 * any dereferences in the callbacks, and never repeat the same
 * allocation.
 */
#define snew_pack(packed,d,nr) if (!MASTER((packed)->cr) && ((packed)->packtype == eptPackSize)) snew((d),(nr))

/*! /brief Helper allocation macro for bcast_state
 */
#define snew_bAlloc(packed,bAlloc,d,nr) if ((bAlloc) && ((packed)->packtype == eptUnpack)) snew((d),(nr))

/*! /brief Data type for callback functions.
 */
typedef void (*t_packed_callback)(gmx_packed_t packed, void *data);

/*!/brief Do the house-keeping for broadcasting structures from master
 *        to other nodes.
 */
static
void pack_and_bcast(const t_commrec *cr, void *data,
                    t_packed_callback callback) {
    gmx_packed_t packed = gmx_packed_init(cr);
    /* How much data will be packed and sent later? */
    callback(packed, data);
    
    /* Allocate */
    gmx_packed_allocate(packed);
    if (MASTER(cr))
    {
        packed->packtype = eptPack;
    }
    else
    {
        packed->packtype = eptUnpack;
        /* receive broadcast before unpacking */
        do_packed_bcast(packed);
    }

    /* Do pack or unpack, as applicable */
    callback(packed, data);

    if (MASTER(cr))
    {
        /* send broadcast after packing */
        do_packed_bcast(packed);
    }
    /* Clean up */
    gmx_packed_destroy(packed);
}

/* The functions for the two stages of broadcasting a t_state follow.
 */

typedef struct state_helper
{
    t_state *state;
    gmx_bool bAlloc;
} t_state_helper;

void state_setup_callback(gmx_packed_t packed, void *data)
{
    t_state_helper *helper = (t_state_helper *) data;
    t_state *state = helper->state;

    gmx_pack_function(packed, &state->natoms, 1, MPI_INT);
    gmx_pack_function(packed, &state->ngtc, 1, MPI_INT);
    gmx_pack_function(packed, &state->nnhpres, 1, MPI_INT);
    gmx_pack_function(packed, &state->nhchainlength, 1, MPI_INT);
    gmx_pack_function(packed, &state->nrng, 1, MPI_INT);
    gmx_pack_function(packed, &state->nrngi, 1, MPI_INT);
    gmx_pack_function(packed, &state->flags, 1, MPI_INT);
    if (state->flags & estDISRE_RM3TAV)
    {
        gmx_pack_function(packed, &state->hist.ndisrepairs, 1, MPI_INT);
    }
    if (state->flags & estORIRE_DTAV)
    {
        gmx_pack_function(packed, &state->hist.norire_Dtav, 1, MPI_INT);
    }
}

/*! /brief Broadcast various integers that determine sizes of arrays,
 *         etc. in a t_state.
 */
void bcast_state_setup(const t_commrec *cr, t_state *state)
{
    t_state_helper helper;

    helper.state = state;
    helper.bAlloc = FALSE;
    pack_and_bcast(cr, &helper, state_setup_callback);
}

void state_callback(gmx_packed_t packed, void *data)
{
    int i, nnht, nnhtp;
    t_state_helper *helper = (t_state_helper *) data;
    t_state *state = helper->state;
    gmx_bool bAlloc = helper->bAlloc;

    nnht = (state->ngtc)*(state->nhchainlength);
    nnhtp = (state->nnhpres)*(state->nhchainlength);

    for(i = 0; i < estNR; i++)
    {
        if (state->flags & (1<<i))
        {
            switch (i)
            {
            case estLAMBDA:      gmx_pack_function(packed, &state->lambda, 1, GMX_MPI_REAL); break;
            case estBOX:         gmx_pack_function(packed, state->box, DIM*DIM, GMX_MPI_REAL); break;
            case estBOX_REL:     gmx_pack_function(packed, state->box_rel, DIM*DIM, GMX_MPI_REAL); break;
            case estBOXV:        gmx_pack_function(packed, state->boxv, DIM*DIM, GMX_MPI_REAL); break;
            case estPRES_PREV:   gmx_pack_function(packed, state->pres_prev, DIM*DIM, GMX_MPI_REAL); break;
            case estSVIR_PREV:   gmx_pack_function(packed, state->svir_prev, DIM*DIM, GMX_MPI_REAL); break;
            case estFVIR_PREV:   gmx_pack_function(packed, state->fvir_prev, DIM*DIM, GMX_MPI_REAL); break;
            case estNH_XI:
                snew_bAlloc(packed, bAlloc, state->nosehoover_xi, nnht);
                gmx_pack_function(packed, state->nosehoover_xi, nnht, MPI_DOUBLE);
                break;
            case estNH_VXI:
                snew_bAlloc(packed, bAlloc, state->nosehoover_vxi, nnht);
                gmx_pack_function(packed, state->nosehoover_vxi, nnht, MPI_DOUBLE);
                break;
            case estNHPRES_XI:
                snew_bAlloc(packed, bAlloc, state->nhpres_xi, nnhtp);
                gmx_pack_function(packed, state->nhpres_xi, nnhtp, MPI_DOUBLE);
                break;
            case estNHPRES_VXI:
                snew_bAlloc(packed, bAlloc, state->nhpres_vxi, nnhtp);
                gmx_pack_function(packed, state->nhpres_vxi, nnhtp, MPI_DOUBLE);
                break;
            case estTC_INT:
                snew_bAlloc(packed, bAlloc, state->therm_integral, state->ngtc);
                gmx_pack_function(packed, state->therm_integral, state->ngtc, MPI_DOUBLE);
                break;
            case estVETA:        gmx_pack_function(packed, &state->veta, 1, GMX_MPI_REAL); break;
            case estVOL0:        gmx_pack_function(packed, &state->vol0, 1, GMX_MPI_REAL); break;
                /* The next four items are generally large enough
                 * to broadcast on their own, rather than take time
                 * packing. */
            case estX:
                snew_bAlloc(packed, bAlloc, state->x, state->natoms);
                if (eptPackSize != packed->packtype)
                {
                    gmx_bcast(state->natoms * sizeof(state->x[0]), state->x, packed->cr);
                }
                break;
            case estV:
                snew_bAlloc(packed, bAlloc, state->v, state->natoms);
                if (eptPackSize != packed->packtype)
                {
                    gmx_bcast(state->natoms * sizeof(state->v[0]), state->v, packed->cr);
                }
                break;
            case estSDX:
                snew_bAlloc(packed, bAlloc, state->sd_X, state->natoms);
                if (eptPackSize != packed->packtype)
                {
                    gmx_bcast(state->natoms * sizeof(state->sd_X[0]), state->sd_X, packed->cr);
                }
                break;
            case estCGP:
                snew_bAlloc(packed, bAlloc, state->cg_p, state->natoms);
                if (eptPackSize != packed->packtype)
                {
                    gmx_bcast(state->natoms * sizeof(state->cg_p[0]), state->cg_p, packed->cr);
                }
                break;
            case estLD_RNG:
                if(1 == state->nrngi)
                {
                    snew_bAlloc(packed, bAlloc, state->ld_rng, state->nrng);
                    gmx_pack_function(packed, state->ld_rng, state->nrng, MPI_UNSIGNED);
                }
                break;
            case estLD_RNGI:
                if(1 == state->nrngi)
                {
                    snew_bAlloc(packed, bAlloc, state->ld_rngi, state->nrngi);
                    gmx_pack_function(packed, state->ld_rngi, state->nrngi, MPI_INT);
                }
                break;
            case estDISRE_INITF:
                gmx_pack_function(packed, &state->hist.disre_initf, 1, GMX_MPI_REAL);
                break;
            case estDISRE_RM3TAV:
                snew_bAlloc(packed, bAlloc, state->hist.disre_rm3tav, state->hist.ndisrepairs);
                gmx_pack_function(packed, state->hist.disre_rm3tav, state->hist.ndisrepairs, GMX_MPI_REAL);
                break;
            case estORIRE_INITF:
                gmx_pack_function(packed, &state->hist.orire_initf, 1, GMX_MPI_REAL);
                break;
            case estORIRE_DTAV:
                snew_bAlloc(packed, bAlloc, state->hist.orire_Dtav, state->hist.norire_Dtav);
                gmx_pack_function(packed, state->hist.orire_Dtav, state->hist.norire_Dtav, GMX_MPI_REAL);
                break;
            default:
                gmx_fatal(FARGS,
                          "Communication is not implemented for %s in bcast_state",
                          est_names[i]);
            }
        }
    }
}

/*! /brief Broadcast a t_state struct to the other processors in the
 *         simulation.
 */
void bcast_state(const t_commrec *cr, t_state *state, gmx_bool bAlloc)
{
    t_state_helper helper;

    if (MASTER(cr))
    {
        bAlloc = FALSE;
    }
    if (bAlloc)
    {
        state->nalloc = state->natoms;
    }

    helper.state = state;
    helper.bAlloc = bAlloc;
    pack_and_bcast(cr, &helper, state_callback);
}

/* Helper functions for doing the bcast of inputrec and mtop follow.
 */

static
void gmx_pack_string(gmx_packed_t packed, t_symtab *symtab, char ***s)
{
    int handle;

    if (eptPack == packed->packtype)
    {
        handle = lookup_symtab(symtab, *s);
    }
    gmx_pack_function(packed, &handle, 1, MPI_INT);
    if (eptUnpack == packed->packtype)
    {
        *s = get_symtab_handle(symtab, handle);
    }
}

static
void gmx_pack_strings(gmx_packed_t packed, t_symtab *symtab, int nr, char ****nm)
{
    int  i;
    int  *handle;

    /* Allocate and pre-process */
    switch(packed->packtype)
    {
    case eptPackSize:
        snew_pack(packed, *nm, nr);
        break;
    case eptPack:
        snew(handle, nr);
        for(i = 0; i < nr; i++)
        {
            handle[i] = lookup_symtab(symtab,(*nm)[i]);
        }
        break;
    case eptUnpack:
        snew(handle, nr);
        break;
    }

    /* Do packing */
    gmx_pack_function(packed, handle, nr, MPI_INT);

    /* Deallocate and post-process */
    switch(packed->packtype)
    {
    case eptPackSize:
        break;
    case eptPack:
        sfree(handle);
        break;
    case eptUnpack:
        for (i = 0; i < nr; i++)
        {
            (*nm)[i] = get_symtab_handle(symtab, handle[i]);
        }
        sfree(handle);
        break;
    }
}

static
void gmx_pack_strings_resinfo(gmx_packed_t packed, t_symtab *symtab,
                                         int nr, t_resinfo *resinfo)
{
    int  i;

    for(i = 0; i < nr; i++)
    {
        int  handle = -1;
        if (eptPack == packed->packtype)
        {
            handle = lookup_symtab(symtab, resinfo[i].name);
        }
        gmx_pack_function(packed, &handle, 1, MPI_INT);
        if (eptUnpack == packed->packtype)
        {
            resinfo[i].name = get_symtab_handle(symtab, handle);
        }
    }
}

static
void gmx_pack_symtab(gmx_packed_t packed, t_symtab *symtab, int *buf_lengths)
{
    int i, len;

    for (i = 0; i < symtab->nr; i++)
    {
        snew_pack(packed, symtab->symbuf->buf[i], buf_lengths[i]);
        gmx_pack_function(packed, symtab->symbuf->buf[i], buf_lengths[i], MPI_CHAR);
    }
}

static
void gmx_pack_block(gmx_packed_t packed, t_block *block)
{

    snew_pack(packed, block->index, block->nr+1);
    gmx_pack_function(packed, block->index, block->nr+1, MPI_INT);
}

static
void gmx_pack_blocka(gmx_packed_t packed, t_blocka *block)
{

    snew_pack(packed, block->index, block->nr+1);
    gmx_pack_function(packed, block->index, block->nr+1, MPI_INT);
    if (block->nra)
    {
        snew_pack(packed, block->a, block->nra);
        gmx_pack_function(packed, block->a, block->nra, MPI_INT);
    }
}

static
void gmx_pack_grps(gmx_packed_t packed, t_grps grps[])
{
    int i;
  
    for(i = 0; i < egcNR; i++)
    {
        snew_pack(packed, grps[i].nm_ind, grps[i].nr);
        gmx_pack_function(packed, grps[i].nm_ind, grps[i].nr, MPI_INT);
    }
}

static
void gmx_pack_atoms(gmx_packed_t packed, t_symtab *symtab, t_atoms *atoms)
{
    int dummy;

    snew_pack(packed, atoms->atom, atoms->nr);
    gmx_pack_function(packed, atoms->atom, sizeof(atoms->atom[0]) * atoms->nr, MPI_BYTE);
    gmx_pack_strings(packed, symtab, atoms->nr, &atoms->atomname);
    snew_pack(packed, atoms->resinfo, atoms->nres);
    gmx_pack_function(packed, atoms->resinfo, sizeof(atoms->resinfo[0]) * atoms->nres, MPI_BYTE);
    gmx_pack_strings_resinfo(packed, symtab, atoms->nres, atoms->resinfo);
    /* QMMM requires atomtypes to be known on all nodes as well */
    gmx_pack_strings(packed, symtab, atoms->nr, &atoms->atomtype);
    gmx_pack_strings(packed, symtab, atoms->nr, &atoms->atomtypeB);
}

static
void gmx_pack_groups(gmx_packed_t packed, t_symtab *symtab,
                                gmx_groups_t *groups, int *grpnr_sizes)
{
    int dummy;
    int g, n;

    gmx_pack_grps(packed, groups->grps);
    gmx_pack_strings(packed, symtab, groups->ngrpname, &groups->grpname);
    for(g = 0; g < egcNR; g++)
    {
        if (0 == grpnr_sizes[g])
        {
            groups->grpnr[g] = NULL;
        }
        else
        {
            snew_pack(packed, groups->grpnr[g], grpnr_sizes[g]);
            gmx_pack_function(packed, groups->grpnr[g], grpnr_sizes[g], MPI_CHAR);
        }
    }
    if (debug) fprintf(debug,"after gmx_pack_groups\n");
}

static
void gmx_pack_ilists(gmx_packed_t packed, t_ilist *ilist)
{
    int ftype;

    for(ftype = 0; ftype < F_NRE; ftype++)
    {
        snew_pack(packed, ilist[ftype].iatoms, ilist[ftype].nr);
        gmx_pack_function(packed, ilist[ftype].iatoms, ilist[ftype].nr, MPI_INT);
    }

    if (debug) fprintf(debug,"after gmx_pack_ilists\n");
}

/* This was in the old code, but not used.
static
void bc_idef(gmx_packed_t packed, ept packtype, t_idef *idef)
{
    block_bc(cr,idef->ntypes);
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

static
void gmx_pack_cmap(gmx_packed_t packed, gmx_cmap_t *cmap_grid)
{
    int i, j,nelem, ngrid;
	
    ngrid = cmap_grid->ngrid;
    nelem = cmap_grid->grid_spacing * cmap_grid->grid_spacing;
	
    if(ngrid > 0)
    {
        snew_pack(packed, cmap_grid->cmapdata, ngrid);
		
        for(i = 0; i < ngrid; i++)
        {
            snew_pack(packed, cmap_grid->cmapdata[i].cmap, 4*nelem);
            gmx_pack_function(packed, cmap_grid->cmapdata[i].cmap, 4*nelem, GMX_MPI_REAL);
        }
    }
}

static
void gmx_pack_ffparams(gmx_packed_t packed, gmx_ffparams_t *ffp)
{
    int i;
  
    gmx_pack_function(packed, &ffp->atnr, 1, MPI_INT);
    snew_pack(packed, ffp->functype, ffp->ntypes);
    snew_pack(packed, ffp->iparams, ffp->ntypes);
    gmx_pack_function(packed, ffp->functype, ffp->ntypes, MPI_INT);
    gmx_pack_function(packed, ffp->iparams, sizeof(ffp->iparams[0]) * ffp->ntypes, MPI_BYTE);
    gmx_pack_function(packed, &ffp->reppow, 1, MPI_DOUBLE);
    gmx_pack_function(packed, &ffp->fudgeQQ, 1, GMX_MPI_REAL);
    gmx_pack_cmap(packed, &ffp->cmap_grid);
}

static
void gmx_pack_grpopts(gmx_packed_t packed, t_grpopts *g)
{
    int i, n;

    snew_pack(packed, g->nrdf, g->ngtc);
    snew_pack(packed, g->tau_t, g->ngtc);
    snew_pack(packed, g->ref_t, g->ngtc);
    snew_pack(packed, g->acc, g->ngacc);
    snew_pack(packed, g->nFreeze, g->ngfrz);
    snew_pack(packed, g->egp_flags, g->ngener*g->ngener);
    gmx_pack_function(packed, g->nrdf, g->ngtc, GMX_MPI_REAL);
    gmx_pack_function(packed, g->tau_t, g->ngtc, GMX_MPI_REAL);
    gmx_pack_function(packed, g->ref_t, g->ngtc, GMX_MPI_REAL);
    gmx_pack_function(packed, g->acc, g->ngacc * DIM, GMX_MPI_REAL);
    gmx_pack_function(packed, g->nFreeze, g->ngfrz * DIM, MPI_INT);
    gmx_pack_function(packed, g->egp_flags, g->ngener*g->ngener, MPI_INT);
    snew_pack(packed, g->annealing, g->ngtc);
    snew_pack(packed, g->anneal_npoints, g->ngtc);
    snew_pack(packed, g->anneal_time, g->ngtc);
    snew_pack(packed, g->anneal_temp, g->ngtc);
    gmx_pack_function(packed, g->annealing, g->ngtc, MPI_INT);
    for(i = 0;i < g->ngtc; i++)
    {
        n = g->anneal_npoints[i];
        if (n > 0)
        {
            snew_pack(packed, g->anneal_time[i], n);
            snew_pack(packed, g->anneal_temp[i], n);
            gmx_pack_function(packed, g->anneal_time[i], n, GMX_MPI_REAL);
            gmx_pack_function(packed, g->anneal_temp[i], n, GMX_MPI_REAL);
        }
    }
        
    /* QMMM stuff, see inputrec */
    snew_pack(packed, g->QMmethod, g->ngQM);
    snew_pack(packed, g->QMbasis, g->ngQM);
    snew_pack(packed, g->QMcharge, g->ngQM);
    snew_pack(packed, g->QMmult, g->ngQM);
    snew_pack(packed, g->bSH, g->ngQM);
    snew_pack(packed, g->CASorbitals, g->ngQM);
    snew_pack(packed, g->CASelectrons, g->ngQM);
    snew_pack(packed, g->SAon, g->ngQM);
    snew_pack(packed, g->SAoff, g->ngQM);
    snew_pack(packed, g->SAsteps, g->ngQM);
    
    if (g->ngQM)
    {
        gmx_pack_function(packed, g->QMmethod, g->ngQM, MPI_INT);
        gmx_pack_function(packed, g->QMbasis, g->ngQM, MPI_INT);
        gmx_pack_function(packed, g->QMcharge, g->ngQM, MPI_INT);
        gmx_pack_function(packed, g->QMmult, g->ngQM, MPI_INT);
        gmx_pack_function(packed, g->bSH, g->ngQM, MPI_INT);
        gmx_pack_function(packed, g->CASorbitals, g->ngQM, MPI_INT);
        gmx_pack_function(packed, g->CASelectrons, g->ngQM, MPI_INT);
        gmx_pack_function(packed, g->SAon, g->ngQM, GMX_MPI_REAL);
        gmx_pack_function(packed, g->SAoff, g->ngQM, GMX_MPI_REAL);
        gmx_pack_function(packed, g->SAsteps, g->ngQM, MPI_INT);
        /* end of QMMM stuff */
    }
}

static
void gmx_pack_cosines(gmx_packed_t packed, t_cosines *cs)
{

    snew_pack(packed, cs->a, cs->n);
    snew_pack(packed, cs->phi, cs->n);
    if (cs->n > 0)
    {
        gmx_pack_function(packed, cs->a, cs->n, GMX_MPI_REAL);
        gmx_pack_function(packed, cs->phi, cs->n, GMX_MPI_REAL);
    }
}

static
void gmx_pack_pullgrp(gmx_packed_t packed, t_pullgrp *pgrp)
{

    if (pgrp->nat > 0)
    {
        snew_pack(packed, pgrp->ind, pgrp->nat);
        gmx_pack_function(packed, pgrp->ind, pgrp->nat, MPI_INT);
    }
    if (pgrp->nweight > 0)
    {
        snew_pack(packed, pgrp->weight, pgrp->nweight);
        gmx_pack_function(packed, pgrp->weight, pgrp->nweight, GMX_MPI_REAL);
    }
}

static
void gmx_pack_pull(gmx_packed_t packed, t_pull *pull)
{
    int g;

    for(g = 0; g < pull->ngrp+1; g++)
    {
        gmx_pack_pullgrp(packed, &pull->grp[g]);
    }
}

static
void gmx_pack_inputrec(gmx_packed_t packed, t_inputrec *inputrec)
{
    int i;

    snew_pack(packed, inputrec->flambda, inputrec->n_flambda);
    gmx_pack_function(packed, inputrec->flambda, inputrec->n_flambda, MPI_DOUBLE);

    gmx_pack_grpopts(packed, &(inputrec->opts));
    if (inputrec->ePull != epullNO)
    {
        gmx_pack_pull(packed, inputrec->pull);
    }
    for(i = 0; i < DIM; i++)
    {
        gmx_pack_cosines(packed, &(inputrec->ex[i]));
        gmx_pack_cosines(packed, &(inputrec->et[i]));
    }
}

static
void gmx_pack_moltype(gmx_packed_t packed, t_symtab *symtab,
                                 gmx_moltype_t *moltype)
{
    gmx_pack_string(packed, symtab, &moltype->name);
    gmx_pack_atoms(packed, symtab, &moltype->atoms);
    if (debug) fprintf(debug,"after gmx_pack_atoms\n");

    gmx_pack_ilists(packed, moltype->ilist);
    gmx_pack_block(packed, &moltype->cgs);
    gmx_pack_blocka(packed, &moltype->excls);
}

static
void gmx_pack_molblock(gmx_packed_t packed, gmx_molblock_t *molb)
{
    gmx_pack_function(packed, &molb->type, 1, MPI_INT);
    gmx_pack_function(packed, &molb->nmol, 1, MPI_INT);
    gmx_pack_function(packed, &molb->natoms_mol, 1, MPI_INT);
    if (molb->nposres_xA > 0)
    {
        snew_pack(packed, molb->posres_xA, molb->nposres_xA);
        gmx_pack_function(packed, molb->posres_xA[0], molb->nposres_xA * DIM, GMX_MPI_REAL);
    }
    if (molb->nposres_xB > 0)
    {
        snew_pack(packed, molb->posres_xB, molb->nposres_xB);
        gmx_pack_function(packed, molb->posres_xB[0], molb->nposres_xB * DIM, GMX_MPI_REAL);
    }
    if (debug) fprintf(debug,"after gmx_pack_molblock\n");
}

static
void gmx_pack_atomtypes(gmx_packed_t packed, t_atomtypes *atomtypes)
{
    int nr;

    nr = atomtypes->nr;

    snew_pack(packed, atomtypes->radius, nr);
    snew_pack(packed, atomtypes->vol, nr);
    snew_pack(packed, atomtypes->surftens, nr);
    snew_pack(packed, atomtypes->gb_radius, nr);
    snew_pack(packed, atomtypes->S_hct, nr);

    gmx_pack_function(packed, atomtypes->radius, nr, GMX_MPI_REAL);
    gmx_pack_function(packed, atomtypes->vol, nr, GMX_MPI_REAL);
    gmx_pack_function(packed, atomtypes->surftens, nr, GMX_MPI_REAL);
    gmx_pack_function(packed, atomtypes->gb_radius, nr, GMX_MPI_REAL);
    gmx_pack_function(packed, atomtypes->S_hct, nr, GMX_MPI_REAL);
}

typedef struct ir_mtop_helper
{
    t_inputrec *inputrec;
    gmx_mtop_t *mtop;
    int *buf_lengths, *grpnr_sizes;
} t_ir_mtop_helper;

void ir_mtop_first_callback(gmx_packed_t packed, void *data)
{
    t_ir_mtop_helper *helper = (t_ir_mtop_helper *) data;
    t_inputrec *inputrec = helper->inputrec;
    gmx_mtop_t *mtop = helper->mtop;
    int *buf_lengths = helper->buf_lengths;
    int *grpnr_sizes = helper->grpnr_sizes;
    int i, g;

    gmx_pack_function(packed, inputrec, sizeof(*inputrec), MPI_BYTE);
    for (i = 0; i < DIM; i++)
    {
        gmx_pack_function(packed, &inputrec->ex[i].n, 1, MPI_INT);
        gmx_pack_function(packed, &inputrec->et[i].n, 1, MPI_INT);
    }
    gmx_pack_function(packed, mtop, sizeof(*mtop), MPI_BYTE);
}

void ir_mtop_second_callback(gmx_packed_t packed, void *data)
{
    t_ir_mtop_helper *helper = (t_ir_mtop_helper *) data;
    t_inputrec *inputrec = helper->inputrec;
    gmx_mtop_t *mtop = helper->mtop;
    int *buf_lengths = helper->buf_lengths;
    int *grpnr_sizes = helper->grpnr_sizes;
    int i, g, ftype;

    /* Now do some allocation for the next stage of data */
    if (packed->packtype == eptUnpack)
    {
        snew(inputrec->opts.anneal_npoints, inputrec->opts.ngtc);
        snew(mtop->symtab.symbuf, 1);
        mtop->symtab.symbuf->bufsize = mtop->symtab.nr;
        snew(mtop->symtab.symbuf->buf, mtop->symtab.nr);
        snew(mtop->moltype, mtop->nmoltype);
        snew(mtop->molblock, mtop->nmolblock);
        if (inputrec->ePull != epullNO)
        {
            snew(inputrec->pull, 1);
        }
    }

    gmx_pack_function(packed, inputrec->opts.anneal_npoints, inputrec->opts.ngtc, MPI_INT);

    if (inputrec->ePull != epullNO)
    {
        gmx_pack_function(packed, inputrec->pull, sizeof(*inputrec->pull), MPI_BYTE);
    }

    for (i = 0; i < mtop->symtab.nr; i++)
    {
        if (eptPack == packed->packtype)
        {
            buf_lengths[i] = strlen(mtop->symtab.symbuf->buf[i]) + 1;
        }
    }
    gmx_pack_function(packed, buf_lengths, mtop->symtab.nr, MPI_INT);

    gmx_pack_function(packed, &mtop->groups.ngrpname, 1, MPI_INT);
    for(g = 0; g < egcNR; g++)
    {
        gmx_pack_function(packed, &mtop->groups.grps[g].nr, 1, MPI_INT);
    }
    for(g = 0; g < egcNR; g++)
    {
        if (eptPack == packed->packtype)
        {
            grpnr_sizes[g] = mtop->groups.grpnr[g] ? mtop->natoms : 0;
        }
    }
    gmx_pack_function(packed, grpnr_sizes, egcNR, MPI_INT);

    for (i = 0; i < mtop->nmoltype; i++)
    {
        gmx_pack_function(packed, &mtop->moltype[i].atoms.nr, 1, MPI_INT);
        gmx_pack_function(packed, &mtop->moltype[i].atoms.nres, 1, MPI_INT);
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            gmx_pack_function(packed, &mtop->moltype[i].ilist[ftype].nr, 1, MPI_INT);
        }
        gmx_pack_function(packed, &mtop->moltype[i].cgs.nr, 1, MPI_INT);
        gmx_pack_function(packed, &mtop->moltype[i].excls.nr, 1, MPI_INT);
        gmx_pack_function(packed, &mtop->moltype[i].excls.nra, 1, MPI_INT);
    }
    for(i = 0; i < mtop->nmolblock; i++)
    {
        gmx_pack_function(packed, &mtop->molblock[i].nposres_xA, 1, MPI_INT);
        gmx_pack_function(packed, &mtop->molblock[i].nposres_xB, 1, MPI_INT);
    }
}

void ir_mtop_third_callback(gmx_packed_t packed, void *data)
{
    t_ir_mtop_helper *helper = (t_ir_mtop_helper *) data;
    t_inputrec *inputrec = helper->inputrec;
    gmx_mtop_t *mtop = helper->mtop;
    int *buf_lengths = helper->buf_lengths;
    int *grpnr_sizes = helper->grpnr_sizes;
    int g;

    snew_pack(packed, inputrec->pull->grp, inputrec->pull->ngrp+1);
    for (g = 0; g < inputrec->pull->ngrp+1; g++)
    {
        gmx_pack_function(packed, &inputrec->pull->grp[g], sizeof(inputrec->pull->grp[g]), MPI_BYTE);
    }
}

void ir_mtop_fourth_callback(gmx_packed_t packed, void *data)
{
    t_ir_mtop_helper *helper = (t_ir_mtop_helper *) data;
    t_inputrec *inputrec = helper->inputrec;
    gmx_mtop_t *mtop = helper->mtop;
    int *buf_lengths = helper->buf_lengths;
    int *grpnr_sizes = helper->grpnr_sizes;
    int i;

    gmx_pack_inputrec(packed, inputrec);
    gmx_pack_symtab(packed, &mtop->symtab, buf_lengths);
    gmx_pack_string(packed, &mtop->symtab, &mtop->name);
    gmx_pack_ffparams(packed, &mtop->ffparams);
    for(i = 0; i < mtop->nmoltype; i++)
    {
        gmx_pack_moltype(packed, &mtop->symtab, &mtop->moltype[i]);
    }
    for(i = 0; i < mtop->nmolblock; i++)
    {
        gmx_pack_molblock(packed, &mtop->molblock[i]);
    }
    gmx_pack_atomtypes(packed, &mtop->atomtypes);
    gmx_pack_block(packed, &mtop->mols);
    gmx_pack_groups(packed, &mtop->symtab, &mtop->groups, grpnr_sizes);
}

/*! /brief Broadcast t_inputrec and gmx_mtop_t structs to the other
 *         processors in the simulation.
 */
void bcast_ir_mtop(const t_commrec *cr, t_inputrec *inputrec, gmx_mtop_t *mtop)
{
    int i, g, ftype;
    t_ir_mtop_helper helper;

    helper.inputrec = inputrec;
    helper.mtop = mtop;
    snew(helper.grpnr_sizes, egcNR);
    helper.buf_lengths = NULL;

    /* First, send around various integers that determine sizes of
     * arrays. */
    pack_and_bcast(cr, &helper, ir_mtop_first_callback);

    snew(helper.buf_lengths, mtop->symtab.nr);

    /* Then, send around some arrays of more sizes, whose lengths were
     * found in the first stage. */
    pack_and_bcast(cr, &helper, ir_mtop_second_callback);

    /* Pull groups require another round of communication */
    if (inputrec->ePull != epullNO)
    {
        pack_and_bcast(cr, &helper, ir_mtop_third_callback);
    }

    /* Finally, flesh it all out */
    pack_and_bcast(cr, &helper, ir_mtop_fourth_callback);

    /* Clean up */
    sfree(helper.buf_lengths);
    sfree(helper.grpnr_sizes);
}
#endif
