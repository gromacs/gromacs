/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "smalloc.h"
#include "symtab.h"
#include "vec.h"
#include "pbc.h"
#include <string.h>

#ifdef GMX_THREAD_MPI
#include "thread_mpi.h"
#endif

/* The source code in this file should be thread-safe.
      Please keep it that way. */



static gmx_bool            bOverAllocDD = FALSE;
#ifdef GMX_THREAD_MPI
static tMPI_Thread_mutex_t over_alloc_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif


void set_over_alloc_dd(gmx_bool set)
{
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&over_alloc_mutex);
    /* we just make sure that we don't set this at the same time.
       We don't worry too much about reading this rarely-set variable */
#endif
    bOverAllocDD = set;
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_unlock(&over_alloc_mutex);
#endif
}

int over_alloc_dd(int n)
{
    if (bOverAllocDD)
    {
        return OVER_ALLOC_FAC*n + 100;
    }
    else
    {
        return n;
    }
}

int gmx_large_int_to_int(gmx_large_int_t step, const char *warn)
{
    int i;

    i = (int)step;

    if (warn != NULL && (step < INT_MIN || step > INT_MAX))
    {
        fprintf(stderr, "\nWARNING during %s:\n", warn);
        fprintf(stderr, "step value ");
        fprintf(stderr, gmx_large_int_pfmt, step);
        fprintf(stderr, " does not fit in int, converted to %d\n\n", i);
    }

    return i;
}

char *gmx_step_str(gmx_large_int_t i, char *buf)
{
    sprintf(buf, gmx_large_int_pfmt, i);

    return buf;
}

void init_block(t_block *block)
{
    int i;

    block->nr           = 0;
    block->nalloc_index = 1;
    snew(block->index, block->nalloc_index);
    block->index[0]     = 0;
}

void init_blocka(t_blocka *block)
{
    int i;

    block->nr           = 0;
    block->nra          = 0;
    block->nalloc_index = 1;
    snew(block->index, block->nalloc_index);
    block->index[0]     = 0;
    block->nalloc_a     = 0;
    block->a            = NULL;
}

void init_atom(t_atoms *at)
{
    int i;

    at->nr        = 0;
    at->nres      = 0;
    at->atom      = NULL;
    at->resinfo   = NULL;
    at->atomname  = NULL;
    at->atomtype  = NULL;
    at->atomtypeB = NULL;
    at->pdbinfo   = NULL;
}

void init_atomtypes(t_atomtypes *at)
{
    at->nr         = 0;
    at->radius     = NULL;
    at->vol        = NULL;
    at->atomnumber = NULL;
    at->gb_radius  = NULL;
    at->S_hct      = NULL;
}

void init_groups(gmx_groups_t *groups)
{
    int g;

    groups->ngrpname = 0;
    groups->grpname  = NULL;
    for (g = 0; (g < egcNR); g++)
    {
        groups->grps[g].nm_ind = NULL;
        groups->ngrpnr[g]      = 0;
        groups->grpnr[g]       = NULL;
    }

}

void init_mtop(gmx_mtop_t *mtop)
{
    mtop->name         = NULL;
    mtop->nmoltype     = 0;
    mtop->moltype      = NULL;
    mtop->nmolblock    = 0;
    mtop->molblock     = NULL;
    mtop->maxres_renum = 0;
    mtop->maxresnr     = -1;
    init_groups(&mtop->groups);
    init_block(&mtop->mols);
    open_symtab(&mtop->symtab);
}

void init_top (t_topology *top)
{
    int i;

    top->name = NULL;
    init_atom (&(top->atoms));
    init_atomtypes(&(top->atomtypes));
    init_block(&top->cgs);
    init_block(&top->mols);
    init_blocka(&top->excls);
    open_symtab(&top->symtab);
}

void init_inputrec(t_inputrec *ir)
{
    memset(ir, 0, (size_t)sizeof(*ir));
    snew(ir->fepvals, 1);
    snew(ir->expandedvals, 1);
    snew(ir->simtempvals, 1);
}

void stupid_fill_block(t_block *grp, int natom, gmx_bool bOneIndexGroup)
{
    int i;

    if (bOneIndexGroup)
    {
        grp->nalloc_index = 2;
        snew(grp->index, grp->nalloc_index);
        grp->index[0] = 0;
        grp->index[1] = natom;
        grp->nr       = 1;
    }
    else
    {
        grp->nalloc_index = natom+1;
        snew(grp->index, grp->nalloc_index);
        snew(grp->index, natom+1);
        for (i = 0; (i <= natom); i++)
        {
            grp->index[i] = i;
        }
        grp->nr = natom;
    }
}

void stupid_fill_blocka(t_blocka *grp, int natom)
{
    int i;

    grp->nalloc_a = natom;
    snew(grp->a, grp->nalloc_a);
    for (i = 0; (i < natom); i++)
    {
        grp->a[i] = i;
    }
    grp->nra = natom;

    grp->nalloc_index = natom + 1;
    snew(grp->index, grp->nalloc_index);
    for (i = 0; (i <= natom); i++)
    {
        grp->index[i] = i;
    }
    grp->nr = natom;
}

void copy_blocka(const t_blocka *src, t_blocka *dest)
{
    int i;

    dest->nr           = src->nr;
    dest->nalloc_index = dest->nr + 1;
    snew(dest->index, dest->nalloc_index);
    for (i = 0; i < dest->nr+1; i++)
    {
        dest->index[i] = src->index[i];
    }
    dest->nra      = src->nra;
    dest->nalloc_a = dest->nra + 1;
    snew(dest->a, dest->nalloc_a);
    for (i = 0; i < dest->nra+1; i++)
    {
        dest->a[i] = src->a[i];
    }
}

void done_block(t_block *block)
{
    block->nr    = 0;
    sfree(block->index);
    block->nalloc_index = 0;
}

void done_blocka(t_blocka *block)
{
    block->nr    = 0;
    block->nra   = 0;
    sfree(block->index);
    if (block->a)
    {
        sfree(block->a);
    }
    block->nalloc_index = 0;
    block->nalloc_a     = 0;
}

void done_atom (t_atoms *at)
{
    at->nr       = 0;
    at->nres     = 0;
    sfree(at->atom);
    sfree(at->resinfo);
    sfree(at->atomname);
    sfree(at->atomtype);
    sfree(at->atomtypeB);
    if (at->pdbinfo)
    {
        sfree(at->pdbinfo);
    }
}

void done_atomtypes(t_atomtypes *atype)
{
    atype->nr = 0;
    sfree(atype->radius);
    sfree(atype->vol);
    sfree(atype->surftens);
    sfree(atype->atomnumber);
    sfree(atype->gb_radius);
    sfree(atype->S_hct);
}

void done_moltype(gmx_moltype_t *molt)
{
    int f;

    done_atom(&molt->atoms);
    done_block(&molt->cgs);
    done_blocka(&molt->excls);

    for (f = 0; f < F_NRE; f++)
    {
        sfree(molt->ilist[f].iatoms);
        molt->ilist[f].nalloc = 0;
    }
}

void done_molblock(gmx_molblock_t *molb)
{
    if (molb->nposres_xA > 0)
    {
        molb->nposres_xA = 0;
        free(molb->posres_xA);
    }
    if (molb->nposres_xB > 0)
    {
        molb->nposres_xB = 0;
        free(molb->posres_xB);
    }
}

void done_mtop(gmx_mtop_t *mtop, gmx_bool bDoneSymtab)
{
    int i;

    if (bDoneSymtab)
    {
        done_symtab(&mtop->symtab);
    }

    sfree(mtop->ffparams.functype);
    sfree(mtop->ffparams.iparams);

    for (i = 0; i < mtop->nmoltype; i++)
    {
        done_moltype(&mtop->moltype[i]);
    }
    sfree(mtop->moltype);
    for (i = 0; i < mtop->nmolblock; i++)
    {
        done_molblock(&mtop->molblock[i]);
    }
    sfree(mtop->molblock);
    done_block(&mtop->mols);
}

void done_top(t_topology *top)
{
    int f;

    sfree(top->idef.functype);
    sfree(top->idef.iparams);
    for (f = 0; f < F_NRE; ++f)
    {
        sfree(top->idef.il[f].iatoms);
        top->idef.il[f].iatoms = NULL;
        top->idef.il[f].nalloc = 0;
    }

    done_atom (&(top->atoms));

    /* For GB */
    done_atomtypes(&(top->atomtypes));

    done_symtab(&(top->symtab));
    done_block(&(top->cgs));
    done_block(&(top->mols));
    done_blocka(&(top->excls));
}

static void done_pullgrp(t_pullgrp *pgrp)
{
    sfree(pgrp->ind);
    sfree(pgrp->ind_loc);
    sfree(pgrp->weight);
    sfree(pgrp->weight_loc);
}

static void done_pull(t_pull *pull)
{
    int i;

    for (i = 0; i < pull->ngrp+1; i++)
    {
        done_pullgrp(pull->grp);
        done_pullgrp(pull->dyna);
    }
}

void done_inputrec(t_inputrec *ir)
{
    int m;

    for (m = 0; (m < DIM); m++)
    {
        if (ir->ex[m].a)
        {
            sfree(ir->ex[m].a);
        }
        if (ir->ex[m].phi)
        {
            sfree(ir->ex[m].phi);
        }
        if (ir->et[m].a)
        {
            sfree(ir->et[m].a);
        }
        if (ir->et[m].phi)
        {
            sfree(ir->et[m].phi);
        }
    }

    sfree(ir->opts.nrdf);
    sfree(ir->opts.ref_t);
    sfree(ir->opts.annealing);
    sfree(ir->opts.anneal_npoints);
    sfree(ir->opts.anneal_time);
    sfree(ir->opts.anneal_temp);
    sfree(ir->opts.tau_t);
    sfree(ir->opts.acc);
    sfree(ir->opts.nFreeze);
    sfree(ir->opts.QMmethod);
    sfree(ir->opts.QMbasis);
    sfree(ir->opts.QMcharge);
    sfree(ir->opts.QMmult);
    sfree(ir->opts.bSH);
    sfree(ir->opts.CASorbitals);
    sfree(ir->opts.CASelectrons);
    sfree(ir->opts.SAon);
    sfree(ir->opts.SAoff);
    sfree(ir->opts.SAsteps);
    sfree(ir->opts.bOPT);
    sfree(ir->opts.bTS);

    if (ir->pull)
    {
        done_pull(ir->pull);
        sfree(ir->pull);
    }
}

static void zero_history(history_t *hist)
{
    hist->disre_initf = 0;
    hist->ndisrepairs = 0;
    hist->disre_rm3tav = NULL;
    hist->orire_initf = 0;
    hist->norire_Dtav = 0;
    hist->orire_Dtav = NULL;
}

static void zero_ekinstate(ekinstate_t *eks)
{
    eks->ekin_n         = 0;
    eks->ekinh          = NULL;
    eks->ekinf          = NULL;
    eks->ekinh_old      = NULL;
    eks->ekinscalef_nhc = NULL;
    eks->ekinscaleh_nhc = NULL;
    eks->vscale_nhc     = NULL;
    eks->dekindl        = 0;
    eks->mvcos          = 0;
}

void init_energyhistory(energyhistory_t * enerhist)
{
    enerhist->nener = 0;

    enerhist->ener_ave     = NULL;
    enerhist->ener_sum     = NULL;
    enerhist->ener_sum_sim = NULL;
    enerhist->dht          = NULL;

    enerhist->nsteps     = 0;
    enerhist->nsum       = 0;
    enerhist->nsteps_sim = 0;
    enerhist->nsum_sim   = 0;

    enerhist->dht = NULL;
}

static void done_delta_h_history(delta_h_history_t *dht)
{
    int i;

    for (i = 0; i < dht->nndh; i++)
    {
        sfree(dht->dh[i]);
    }
    sfree(dht->dh);
    sfree(dht->ndh);
}

void done_energyhistory(energyhistory_t * enerhist)
{
    sfree(enerhist->ener_ave);
    sfree(enerhist->ener_sum);
    sfree(enerhist->ener_sum_sim);

    if (enerhist->dht != NULL)
    {
        done_delta_h_history(enerhist->dht);
        sfree(enerhist->dht);
    }
}

void init_gtc_state(t_state *state, int ngtc, int nnhpres, int nhchainlength)
{
    int i, j;

    state->ngtc          = ngtc;
    state->nnhpres       = nnhpres;
    state->nhchainlength = nhchainlength;
    if (state->ngtc > 0)
    {
        snew(state->nosehoover_xi, state->nhchainlength*state->ngtc);
        snew(state->nosehoover_vxi, state->nhchainlength*state->ngtc);
        snew(state->therm_integral, state->ngtc);
        for (i = 0; i < state->ngtc; i++)
        {
            for (j = 0; j < state->nhchainlength; j++)
            {
                state->nosehoover_xi[i*state->nhchainlength + j]   = 0.0;
                state->nosehoover_vxi[i*state->nhchainlength + j]  = 0.0;
            }
        }
        for (i = 0; i < state->ngtc; i++)
        {
            state->therm_integral[i]  = 0.0;
        }
    }
    else
    {
        state->nosehoover_xi  = NULL;
        state->nosehoover_vxi = NULL;
        state->therm_integral = NULL;
    }

    if (state->nnhpres > 0)
    {
        snew(state->nhpres_xi, state->nhchainlength*nnhpres);
        snew(state->nhpres_vxi, state->nhchainlength*nnhpres);
        for (i = 0; i < nnhpres; i++)
        {
            for (j = 0; j < state->nhchainlength; j++)
            {
                state->nhpres_xi[i*nhchainlength + j]   = 0.0;
                state->nhpres_vxi[i*nhchainlength + j]  = 0.0;
            }
        }
    }
    else
    {
        state->nhpres_xi  = NULL;
        state->nhpres_vxi = NULL;
    }
}


void init_state(t_state *state, int natoms, int ngtc, int nnhpres, int nhchainlength, int nlambda)
{
    int i;

    state->natoms = natoms;
    state->nrng   = 0;
    state->flags  = 0;
    state->lambda = 0;
    snew(state->lambda, efptNR);
    for (i = 0; i < efptNR; i++)
    {
        state->lambda[i] = 0;
    }
    state->veta   = 0;
    clear_mat(state->box);
    clear_mat(state->box_rel);
    clear_mat(state->boxv);
    clear_mat(state->pres_prev);
    clear_mat(state->svir_prev);
    clear_mat(state->fvir_prev);
    init_gtc_state(state, ngtc, nnhpres, nhchainlength);
    state->nalloc = state->natoms;
    if (state->nalloc > 0)
    {
        snew(state->x, state->nalloc);
        snew(state->v, state->nalloc);
    }
    else
    {
        state->x = NULL;
        state->v = NULL;
    }
    state->sd_X = NULL;
    state->cg_p = NULL;
    zero_history(&state->hist);
    zero_ekinstate(&state->ekinstate);
    init_energyhistory(&state->enerhist);
    init_df_history(&state->dfhist, nlambda, 0);
    state->ddp_count       = 0;
    state->ddp_count_cg_gl = 0;
    state->cg_gl           = NULL;
    state->cg_gl_nalloc    = 0;
}

void done_state(t_state *state)
{
    if (state->nosehoover_xi)
    {
        sfree(state->nosehoover_xi);
    }
    if (state->x)
    {
        sfree(state->x);
    }
    if (state->v)
    {
        sfree(state->v);
    }
    if (state->sd_X)
    {
        sfree(state->sd_X);
    }
    if (state->cg_p)
    {
        sfree(state->cg_p);
    }
    state->nalloc = 0;
    if (state->cg_gl)
    {
        sfree(state->cg_gl);
    }
    state->cg_gl_nalloc = 0;
}

static void do_box_rel(t_inputrec *ir, matrix box_rel, matrix b, gmx_bool bInit)
{
    int d, d2;

    for (d = YY; d <= ZZ; d++)
    {
        for (d2 = XX; d2 <= (ir->epct == epctSEMIISOTROPIC ? YY : ZZ); d2++)
        {
            /* We need to check if this box component is deformed
             * or if deformation of another component might cause
             * changes in this component due to box corrections.
             */
            if (ir->deform[d][d2] == 0 &&
                !(d == ZZ && d2 == XX && ir->deform[d][YY] != 0 &&
                  (b[YY][d2] != 0 || ir->deform[YY][d2] != 0)))
            {
                if (bInit)
                {
                    box_rel[d][d2] = b[d][d2]/b[XX][XX];
                }
                else
                {
                    b[d][d2] = b[XX][XX]*box_rel[d][d2];
                }
            }
        }
    }
}

void set_box_rel(t_inputrec *ir, t_state *state)
{
    /* Make sure the box obeys the restrictions before we fix the ratios */
    correct_box(NULL, 0, state->box, NULL);

    clear_mat(state->box_rel);

    if (PRESERVE_SHAPE(*ir))
    {
        do_box_rel(ir, state->box_rel, state->box, TRUE);
    }
}

void preserve_box_shape(t_inputrec *ir, matrix box_rel, matrix b)
{
    if (PRESERVE_SHAPE(*ir))
    {
        do_box_rel(ir, box_rel, b, FALSE);
    }
}

void add_t_atoms(t_atoms *atoms, int natom_extra, int nres_extra)
{
    int i;

    if (natom_extra > 0)
    {
        srenew(atoms->atomname, atoms->nr+natom_extra);
        srenew(atoms->atom, atoms->nr+natom_extra);
        if (NULL != atoms->pdbinfo)
        {
            srenew(atoms->pdbinfo, atoms->nr+natom_extra);
        }
        if (NULL != atoms->atomtype)
        {
            srenew(atoms->atomtype, atoms->nr+natom_extra);
        }
        if (NULL != atoms->atomtypeB)
        {
            srenew(atoms->atomtypeB, atoms->nr+natom_extra);
        }
        for (i = atoms->nr; (i < atoms->nr+natom_extra); i++)
        {
            atoms->atomname[i] = NULL;
            memset(&atoms->atom[i], 0, sizeof(atoms->atom[i]));
            if (NULL != atoms->pdbinfo)
            {
                memset(&atoms->pdbinfo[i], 0, sizeof(atoms->pdbinfo[i]));
            }
            if (NULL != atoms->atomtype)
            {
                atoms->atomtype[i] = NULL;
            }
            if (NULL != atoms->atomtypeB)
            {
                atoms->atomtypeB[i] = NULL;
            }
        }
        atoms->nr += natom_extra;
    }
    if (nres_extra > 0)
    {
        srenew(atoms->resinfo, atoms->nres+nres_extra);
        for (i = atoms->nres; (i < atoms->nres+nres_extra); i++)
        {
            memset(&atoms->resinfo[i], 0, sizeof(atoms->resinfo[i]));
        }
        atoms->nres += nres_extra;
    }
}

void init_t_atoms(t_atoms *atoms, int natoms, gmx_bool bPdbinfo)
{
    atoms->nr   = natoms;
    atoms->nres = 0;
    snew(atoms->atomname, natoms);
    atoms->atomtype  = NULL;
    atoms->atomtypeB = NULL;
    snew(atoms->resinfo, natoms);
    snew(atoms->atom, natoms);
    if (bPdbinfo)
    {
        snew(atoms->pdbinfo, natoms);
    }
    else
    {
        atoms->pdbinfo = NULL;
    }
}

t_atoms *copy_t_atoms(t_atoms *src)
{
    t_atoms *dst;
    int      i;

    snew(dst, 1);
    init_t_atoms(dst, src->nr, (NULL != src->pdbinfo));
    dst->nr = src->nr;
    if (NULL != src->atomname)
    {
        snew(dst->atomname, src->nr);
    }
    if (NULL != src->atomtype)
    {
        snew(dst->atomtype, src->nr);
    }
    if (NULL != src->atomtypeB)
    {
        snew(dst->atomtypeB, src->nr);
    }
    for (i = 0; (i < src->nr); i++)
    {
        dst->atom[i] = src->atom[i];
        if (NULL != src->pdbinfo)
        {
            dst->pdbinfo[i] = src->pdbinfo[i];
        }
        if (NULL != src->atomname)
        {
            dst->atomname[i]  = src->atomname[i];
        }
        if (NULL != src->atomtype)
        {
            dst->atomtype[i] = src->atomtype[i];
        }
        if (NULL != src->atomtypeB)
        {
            dst->atomtypeB[i] = src->atomtypeB[i];
        }
    }
    dst->nres = src->nres;
    for (i = 0; (i < src->nres); i++)
    {
        dst->resinfo[i] = src->resinfo[i];
    }
    return dst;
}

void t_atoms_set_resinfo(t_atoms *atoms, int atom_ind, t_symtab *symtab,
                         const char *resname, int resnr, unsigned char ic,
                         int chainnum, char chainid)
{
    t_resinfo *ri;

    ri           = &atoms->resinfo[atoms->atom[atom_ind].resind];
    ri->name     = put_symtab(symtab, resname);
    ri->rtp      = NULL;
    ri->nr       = resnr;
    ri->ic       = ic;
    ri->chainnum = chainnum;
    ri->chainid  = chainid;
}

void free_t_atoms(t_atoms *atoms, gmx_bool bFreeNames)
{
    int i;

    if (bFreeNames)
    {
        for (i = 0; i < atoms->nr; i++)
        {
            sfree(*atoms->atomname[i]);
            *atoms->atomname[i] = NULL;
        }
        for (i = 0; i < atoms->nres; i++)
        {
            sfree(*atoms->resinfo[i].name);
            *atoms->resinfo[i].name = NULL;
        }
    }
    sfree(atoms->atomname);
    /* Do we need to free atomtype and atomtypeB as well ? */
    sfree(atoms->resinfo);
    sfree(atoms->atom);
    if (atoms->pdbinfo)
    {
        sfree(atoms->pdbinfo);
    }
    atoms->nr   = 0;
    atoms->nres = 0;
}

real max_cutoff(real cutoff1, real cutoff2)
{
    if (cutoff1 == 0 || cutoff2 == 0)
    {
        return 0;
    }
    else
    {
        return max(cutoff1, cutoff2);
    }
}

extern void init_df_history(df_history_t *dfhist, int nlambda, real wl_delta)
{
    int i;

    dfhist->bEquil   = 0;
    dfhist->nlambda  = nlambda;
    dfhist->wl_delta = wl_delta;
    snew(dfhist->sum_weights, dfhist->nlambda);
    snew(dfhist->sum_dg, dfhist->nlambda);
    snew(dfhist->sum_minvar, dfhist->nlambda);
    snew(dfhist->sum_variance, dfhist->nlambda);
    snew(dfhist->n_at_lam, dfhist->nlambda);
    snew(dfhist->wl_histo, dfhist->nlambda);

    /* allocate transition matrices here */
    snew(dfhist->Tij, dfhist->nlambda);
    snew(dfhist->Tij_empirical, dfhist->nlambda);

    for (i = 0; i < dfhist->nlambda; i++)
    {
        snew(dfhist->Tij[i], dfhist->nlambda);
        snew(dfhist->Tij_empirical[i], dfhist->nlambda);
    }

    snew(dfhist->accum_p, dfhist->nlambda);
    snew(dfhist->accum_m, dfhist->nlambda);
    snew(dfhist->accum_p2, dfhist->nlambda);
    snew(dfhist->accum_m2, dfhist->nlambda);

    for (i = 0; i < dfhist->nlambda; i++)
    {
        snew((dfhist->accum_p)[i], dfhist->nlambda);
        snew((dfhist->accum_m)[i], dfhist->nlambda);
        snew((dfhist->accum_p2)[i], dfhist->nlambda);
        snew((dfhist->accum_m2)[i], dfhist->nlambda);
    }
}

extern void copy_df_history(df_history_t *df_dest, df_history_t *df_source)
{
    int i, j;

    init_df_history(df_dest, df_source->nlambda, df_source->wl_delta);
    df_dest->nlambda = df_source->nlambda;
    df_dest->bEquil  = df_source->bEquil;
    for (i = 0; i < df_dest->nlambda; i++)
    {
        df_dest->sum_weights[i]  = df_source->sum_weights[i];
        df_dest->sum_dg[i]       = df_source->sum_dg[i];
        df_dest->sum_minvar[i]   = df_source->sum_minvar[i];
        df_dest->sum_variance[i] = df_source->sum_variance[i];
        df_dest->n_at_lam[i]     = df_source->n_at_lam[i];
        df_dest->wl_histo[i]     = df_source->wl_histo[i];
        df_dest->accum_p[i]      = df_source->accum_p[i];
        df_dest->accum_m[i]      = df_source->accum_m[i];
        df_dest->accum_p2[i]     = df_source->accum_p2[i];
        df_dest->accum_m2[i]     = df_source->accum_m2[i];
    }

    for (i = 0; i < df_dest->nlambda; i++)
    {
        for (j = 0; j < df_dest->nlambda; j++)
        {
            df_dest->Tij[i][j]            = df_source->Tij[i][j];
            df_dest->Tij_empirical[i][j]  = df_source->Tij_empirical[i][j];
        }
    }
}
