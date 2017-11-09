/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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

#include "broadcaststructs.h"

#include <string.h>

#include "gromacs/gmxlib/network.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/mdrun.h"
#include "gromacs/mdlib/tgroup.h"
#include "gromacs/mdtypes/awh-params.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/pull-params.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/inmemoryserializer.h"
#include "gromacs/utility/keyvaluetree.h"
#include "gromacs/utility/keyvaluetreeserializer.h"
#include "gromacs/utility/smalloc.h"

static void bc_cstring(const t_commrec *cr, char **s)
{
    int size = 0;

    if (MASTER(cr) && *s != nullptr)
    {
        /* Size of the char buffer is string length + 1 for '\0' */
        size = strlen(*s) + 1;
    }
    block_bc(cr, size);
    if (size > 0)
    {
        if (!MASTER(cr))
        {
            srenew(*s, size);
        }
        nblock_bc(cr, size, *s);
    }
    else if (!MASTER(cr) && *s != nullptr)
    {
        sfree(*s);
        *s = nullptr;
    }
}

static void bc_string(const t_commrec *cr, t_symtab *symtab, char ***s)
{
    int handle;

    if (MASTER(cr))
    {
        handle = lookup_symtab(symtab, *s);
    }
    block_bc(cr, handle);
    if (!MASTER(cr))
    {
        *s = get_symtab_handle(symtab, handle);
    }
}

static void bc_strings(const t_commrec *cr, t_symtab *symtab, int nr, char ****nm)
{
    int     i;
    int    *handle;

    snew(handle, nr);
    if (MASTER(cr))
    {
        for (i = 0; (i < nr); i++)
        {
            handle[i] = lookup_symtab(symtab, (*nm)[i]);
        }
    }
    nblock_bc(cr, nr, handle);

    if (!MASTER(cr))
    {
        snew_bc(cr, *nm, nr);
        for (i = 0; (i < nr); i++)
        {
            (*nm)[i] = get_symtab_handle(symtab, handle[i]);
        }
    }
    sfree(handle);
}

static void bc_strings_resinfo(const t_commrec *cr, t_symtab *symtab,
                               int nr, t_resinfo *resinfo)
{
    int   i;
    int  *handle;

    snew(handle, nr);
    if (MASTER(cr))
    {
        for (i = 0; (i < nr); i++)
        {
            handle[i] = lookup_symtab(symtab, resinfo[i].name);
        }
    }
    nblock_bc(cr, nr, handle);

    if (!MASTER(cr))
    {
        for (i = 0; (i < nr); i++)
        {
            resinfo[i].name = get_symtab_handle(symtab, handle[i]);
        }
    }
    sfree(handle);
}

static void bc_symtab(const t_commrec *cr, t_symtab *symtab)
{
    int       i, nr, len;
    t_symbuf *symbuf;

    block_bc(cr, symtab->nr);
    nr = symtab->nr;
    snew_bc(cr, symtab->symbuf, 1);
    symbuf          = symtab->symbuf;
    symbuf->bufsize = nr;
    snew_bc(cr, symbuf->buf, nr);
    for (i = 0; i < nr; i++)
    {
        if (MASTER(cr))
        {
            len = strlen(symbuf->buf[i]) + 1;
        }
        block_bc(cr, len);
        snew_bc(cr, symbuf->buf[i], len);
        nblock_bc(cr, len, symbuf->buf[i]);
    }
}

static void bc_block(const t_commrec *cr, t_block *block)
{
    block_bc(cr, block->nr);
    snew_bc(cr, block->index, block->nr+1);
    nblock_bc(cr, block->nr+1, block->index);
}

static void bc_blocka(const t_commrec *cr, t_blocka *block)
{
    block_bc(cr, block->nr);
    snew_bc(cr, block->index, block->nr+1);
    nblock_bc(cr, block->nr+1, block->index);
    block_bc(cr, block->nra);
    if (block->nra)
    {
        snew_bc(cr, block->a, block->nra);
        nblock_bc(cr, block->nra, block->a);
    }
}

static void bc_grps(const t_commrec *cr, t_grps grps[])
{
    int i;

    for (i = 0; (i < egcNR); i++)
    {
        block_bc(cr, grps[i].nr);
        snew_bc(cr, grps[i].nm_ind, grps[i].nr);
        nblock_bc(cr, grps[i].nr, grps[i].nm_ind);
    }
}

static void bc_atoms(const t_commrec *cr, t_symtab *symtab, t_atoms *atoms)
{
    block_bc(cr, atoms->nr);
    snew_bc(cr, atoms->atom, atoms->nr);
    nblock_bc(cr, atoms->nr, atoms->atom);
    bc_strings(cr, symtab, atoms->nr, &atoms->atomname);
    block_bc(cr, atoms->nres);
    snew_bc(cr, atoms->resinfo, atoms->nres);
    nblock_bc(cr, atoms->nres, atoms->resinfo);
    bc_strings_resinfo(cr, symtab, atoms->nres, atoms->resinfo);
    /* QMMM requires atomtypes to be known on all nodes as well */
    bc_strings(cr, symtab, atoms->nr, &atoms->atomtype);
    bc_strings(cr, symtab, atoms->nr, &atoms->atomtypeB);
}

static void bc_groups(const t_commrec *cr, t_symtab *symtab,
                      int natoms, gmx_groups_t *groups)
{
    int g, n;

    bc_grps(cr, groups->grps);
    block_bc(cr, groups->ngrpname);
    bc_strings(cr, symtab, groups->ngrpname, &groups->grpname);
    for (g = 0; g < egcNR; g++)
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
        block_bc(cr, n);
        if (n == 0)
        {
            groups->grpnr[g] = nullptr;
        }
        else
        {
            snew_bc(cr, groups->grpnr[g], n);
            nblock_bc(cr, n, groups->grpnr[g]);
        }
    }
    if (debug)
    {
        fprintf(debug, "after bc_groups\n");
    }
}

template <typename AllocatorType>
static void bcastPaddedRVecVector(const t_commrec *cr, std::vector<gmx::RVec, AllocatorType> *v, int numAtoms)
{
    v->resize(gmx::paddedRVecVectorSize(numAtoms));
    nblock_bc(cr, numAtoms, as_rvec_array(v->data()));
}

void broadcastStateWithoutDynamics(const t_commrec *cr, t_state *state)
{
    GMX_RELEASE_ASSERT(!DOMAINDECOMP(cr), "broadcastStateWithoutDynamics should only be used for special cases without domain decomposition");

    if (!PAR(cr))
    {
        return;
    }

    /* Broadcasts the state sizes and flags from the master to all ranks
     * in cr->mpi_comm_mygroup.
     */
    block_bc(cr, state->natoms);
    block_bc(cr, state->flags);

    for (int i = 0; i < estNR; i++)
    {
        if (state->flags & (1 << i))
        {
            switch (i)
            {
                case estLAMBDA:
                    nblock_bc(cr, efptNR, state->lambda.data());
                    break;
                case estFEPSTATE:
                    block_bc(cr, state->fep_state);
                    break;
                case estBOX:
                    block_bc(cr, state->box);
                    break;
                case estX:
                    bcastPaddedRVecVector(cr, &state->x, state->natoms);
                    break;
                default:
                    GMX_RELEASE_ASSERT(false, "The state has a dynamic entry, while no dynamic entries should be present");
                    break;
            }
        }
    }
}

static void bc_ilists(const t_commrec *cr, t_ilist *ilist)
{
    int ftype;

    /* Here we only communicate the non-zero length ilists */
    if (MASTER(cr))
    {
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            if (ilist[ftype].nr > 0)
            {
                block_bc(cr, ftype);
                block_bc(cr, ilist[ftype].nr);
                nblock_bc(cr, ilist[ftype].nr, ilist[ftype].iatoms);
            }
        }
        ftype = -1;
        block_bc(cr, ftype);
    }
    else
    {
        for (ftype = 0; ftype < F_NRE; ftype++)
        {
            ilist[ftype].nr = 0;
        }
        do
        {
            block_bc(cr, ftype);
            if (ftype >= 0)
            {
                block_bc(cr, ilist[ftype].nr);
                snew_bc(cr, ilist[ftype].iatoms, ilist[ftype].nr);
                nblock_bc(cr, ilist[ftype].nr, ilist[ftype].iatoms);
            }
        }
        while (ftype >= 0);
    }

    if (debug)
    {
        fprintf(debug, "after bc_ilists\n");
    }
}

static void bc_cmap(const t_commrec *cr, gmx_cmap_t *cmap_grid)
{
    int i, nelem, ngrid;

    block_bc(cr, cmap_grid->ngrid);
    block_bc(cr, cmap_grid->grid_spacing);

    ngrid = cmap_grid->ngrid;
    nelem = cmap_grid->grid_spacing * cmap_grid->grid_spacing;

    if (ngrid > 0)
    {
        snew_bc(cr, cmap_grid->cmapdata, ngrid);

        for (i = 0; i < ngrid; i++)
        {
            snew_bc(cr, cmap_grid->cmapdata[i].cmap, 4*nelem);
            nblock_bc(cr, 4*nelem, cmap_grid->cmapdata[i].cmap);
        }
    }
}

static void bc_ffparams(const t_commrec *cr, gmx_ffparams_t *ffp)
{
    block_bc(cr, ffp->ntypes);
    block_bc(cr, ffp->atnr);
    snew_bc(cr, ffp->functype, ffp->ntypes);
    snew_bc(cr, ffp->iparams, ffp->ntypes);
    nblock_bc(cr, ffp->ntypes, ffp->functype);
    nblock_bc(cr, ffp->ntypes, ffp->iparams);
    block_bc(cr, ffp->reppow);
    block_bc(cr, ffp->fudgeQQ);
    bc_cmap(cr, &ffp->cmap_grid);
}

static void bc_grpopts(const t_commrec *cr, t_grpopts *g)
{
    int i, n;

    block_bc(cr, g->ngtc);
    block_bc(cr, g->ngacc);
    block_bc(cr, g->ngfrz);
    block_bc(cr, g->ngener);
    snew_bc(cr, g->nrdf, g->ngtc);
    snew_bc(cr, g->tau_t, g->ngtc);
    snew_bc(cr, g->ref_t, g->ngtc);
    snew_bc(cr, g->acc, g->ngacc);
    snew_bc(cr, g->nFreeze, g->ngfrz);
    snew_bc(cr, g->egp_flags, g->ngener*g->ngener);

    nblock_bc(cr, g->ngtc, g->nrdf);
    nblock_bc(cr, g->ngtc, g->tau_t);
    nblock_bc(cr, g->ngtc, g->ref_t);
    nblock_bc(cr, g->ngacc, g->acc);
    nblock_bc(cr, g->ngfrz, g->nFreeze);
    nblock_bc(cr, g->ngener*g->ngener, g->egp_flags);
    snew_bc(cr, g->annealing, g->ngtc);
    snew_bc(cr, g->anneal_npoints, g->ngtc);
    snew_bc(cr, g->anneal_time, g->ngtc);
    snew_bc(cr, g->anneal_temp, g->ngtc);
    nblock_bc(cr, g->ngtc, g->annealing);
    nblock_bc(cr, g->ngtc, g->anneal_npoints);
    for (i = 0; (i < g->ngtc); i++)
    {
        n = g->anneal_npoints[i];
        if (n > 0)
        {
            snew_bc(cr, g->anneal_time[i], n);
            snew_bc(cr, g->anneal_temp[i], n);
            nblock_bc(cr, n, g->anneal_time[i]);
            nblock_bc(cr, n, g->anneal_temp[i]);
        }
    }

    /* QMMM stuff, see inputrec */
    block_bc(cr, g->ngQM);
    snew_bc(cr, g->QMmethod, g->ngQM);
    snew_bc(cr, g->QMbasis, g->ngQM);
    snew_bc(cr, g->QMcharge, g->ngQM);
    snew_bc(cr, g->QMmult, g->ngQM);
    snew_bc(cr, g->bSH, g->ngQM);
    snew_bc(cr, g->CASorbitals, g->ngQM);
    snew_bc(cr, g->CASelectrons, g->ngQM);
    snew_bc(cr, g->SAon, g->ngQM);
    snew_bc(cr, g->SAoff, g->ngQM);
    snew_bc(cr, g->SAsteps, g->ngQM);

    if (g->ngQM)
    {
        nblock_bc(cr, g->ngQM, g->QMmethod);
        nblock_bc(cr, g->ngQM, g->QMbasis);
        nblock_bc(cr, g->ngQM, g->QMcharge);
        nblock_bc(cr, g->ngQM, g->QMmult);
        nblock_bc(cr, g->ngQM, g->bSH);
        nblock_bc(cr, g->ngQM, g->CASorbitals);
        nblock_bc(cr, g->ngQM, g->CASelectrons);
        nblock_bc(cr, g->ngQM, g->SAon);
        nblock_bc(cr, g->ngQM, g->SAoff);
        nblock_bc(cr, g->ngQM, g->SAsteps);
        /* end of QMMM stuff */
    }
}

static void bc_awhBias(const t_commrec *cr, gmx::AwhBiasParams *awhBiasParams)
{
    block_bc(cr, *awhBiasParams);

    snew_bc(cr, awhBiasParams->dimParams, awhBiasParams->ndim);
    nblock_bc(cr, awhBiasParams->ndim, awhBiasParams->dimParams);
}

static void bc_awh(const t_commrec *cr, gmx::AwhParams *awhParams)
{
    int k;

    block_bc(cr, *awhParams);
    snew_bc(cr, awhParams->awhBiasParams, awhParams->numBias);
    for (k = 0; k < awhParams->numBias; k++)
    {
        bc_awhBias(cr, &awhParams->awhBiasParams[k]);
    }
}

static void bc_pull_group(const t_commrec *cr, t_pull_group *pgrp)
{
    block_bc(cr, *pgrp);
    if (pgrp->nat > 0)
    {
        snew_bc(cr, pgrp->ind, pgrp->nat);
        nblock_bc(cr, pgrp->nat, pgrp->ind);
    }
    if (pgrp->nweight > 0)
    {
        snew_bc(cr, pgrp->weight, pgrp->nweight);
        nblock_bc(cr, pgrp->nweight, pgrp->weight);
    }
}

static void bc_pull(const t_commrec *cr, pull_params_t *pull)
{
    int g;

    block_bc(cr, *pull);
    snew_bc(cr, pull->group, pull->ngroup);
    for (g = 0; g < pull->ngroup; g++)
    {
        bc_pull_group(cr, &pull->group[g]);
    }
    snew_bc(cr, pull->coord, pull->ncoord);
    nblock_bc(cr, pull->ncoord, pull->coord);
    for (int c = 0; c < pull->ncoord; c++)
    {
        if (!MASTER(cr))
        {
            pull->coord[c].externalPotentialProvider = nullptr;
        }
        if (pull->coord[c].eType == epullEXTERNAL)
        {
            bc_cstring(cr, &pull->coord[c].externalPotentialProvider);
        }
    }
}

static void bc_rotgrp(const t_commrec *cr, t_rotgrp *rotg)
{
    block_bc(cr, *rotg);
    if (rotg->nat > 0)
    {
        snew_bc(cr, rotg->ind, rotg->nat);
        nblock_bc(cr, rotg->nat, rotg->ind);
        snew_bc(cr, rotg->x_ref, rotg->nat);
        nblock_bc(cr, rotg->nat, rotg->x_ref);
    }
}

static void bc_rot(const t_commrec *cr, t_rot *rot)
{
    int g;

    block_bc(cr, *rot);
    snew_bc(cr, rot->grp, rot->ngrp);
    for (g = 0; g < rot->ngrp; g++)
    {
        bc_rotgrp(cr, &rot->grp[g]);
    }
}

static void bc_imd(const t_commrec *cr, t_IMD *imd)
{
    block_bc(cr, *imd);
    snew_bc(cr, imd->ind, imd->nat);
    nblock_bc(cr, imd->nat, imd->ind);
}

static void bc_fepvals(const t_commrec *cr, t_lambda *fep)
{
    int      i;

    block_bc(cr, fep->nstdhdl);
    block_bc(cr, fep->init_lambda);
    block_bc(cr, fep->init_fep_state);
    block_bc(cr, fep->delta_lambda);
    block_bc(cr, fep->edHdLPrintEnergy);
    block_bc(cr, fep->n_lambda);
    if (fep->n_lambda > 0)
    {
        snew_bc(cr, fep->all_lambda, efptNR);
        nblock_bc(cr, efptNR, fep->all_lambda);
        for (i = 0; i < efptNR; i++)
        {
            snew_bc(cr, fep->all_lambda[i], fep->n_lambda);
            nblock_bc(cr, fep->n_lambda, fep->all_lambda[i]);
        }
    }
    block_bc(cr, fep->sc_alpha);
    block_bc(cr, fep->sc_power);
    block_bc(cr, fep->sc_r_power);
    block_bc(cr, fep->sc_sigma);
    block_bc(cr, fep->sc_sigma_min);
    block_bc(cr, fep->bScCoul);
    nblock_bc(cr, efptNR, &(fep->separate_dvdl[0]));
    block_bc(cr, fep->dhdl_derivatives);
    block_bc(cr, fep->dh_hist_size);
    block_bc(cr, fep->dh_hist_spacing);
    if (debug)
    {
        fprintf(debug, "after bc_fepvals\n");
    }
}

static void bc_expandedvals(const t_commrec *cr, t_expanded *expand, int n_lambda)
{
    block_bc(cr, expand->nstexpanded);
    block_bc(cr, expand->elamstats);
    block_bc(cr, expand->elmcmove);
    block_bc(cr, expand->elmceq);
    block_bc(cr, expand->equil_n_at_lam);
    block_bc(cr, expand->equil_wl_delta);
    block_bc(cr, expand->equil_ratio);
    block_bc(cr, expand->equil_steps);
    block_bc(cr, expand->equil_samples);
    block_bc(cr, expand->lmc_seed);
    block_bc(cr, expand->minvar);
    block_bc(cr, expand->minvar_const);
    block_bc(cr, expand->c_range);
    block_bc(cr, expand->bSymmetrizedTMatrix);
    block_bc(cr, expand->nstTij);
    block_bc(cr, expand->lmc_repeats);
    block_bc(cr, expand->lmc_forced_nstart);
    block_bc(cr, expand->gibbsdeltalam);
    block_bc(cr, expand->wl_scale);
    block_bc(cr, expand->wl_ratio);
    block_bc(cr, expand->init_wl_delta);
    block_bc(cr, expand->bInit_weights);
    snew_bc(cr, expand->init_lambda_weights, n_lambda);
    nblock_bc(cr, n_lambda, expand->init_lambda_weights);
    block_bc(cr, expand->mc_temp);
    if (debug)
    {
        fprintf(debug, "after bc_expandedvals\n");
    }
}

static void bc_simtempvals(const t_commrec *cr, t_simtemp *simtemp, int n_lambda)
{
    block_bc(cr, simtemp->simtemp_low);
    block_bc(cr, simtemp->simtemp_high);
    block_bc(cr, simtemp->eSimTempScale);
    snew_bc(cr, simtemp->temperatures, n_lambda);
    nblock_bc(cr, n_lambda, simtemp->temperatures);
    if (debug)
    {
        fprintf(debug, "after bc_simtempvals\n");
    }
}


static void bc_swapions(const t_commrec *cr, t_swapcoords *swap)
{
    block_bc(cr, *swap);

    /* Broadcast atom indices for split groups, solvent group, and for all user-defined swap groups */
    snew_bc(cr, swap->grp, swap->ngrp);
    for (int i = 0; i < swap->ngrp; i++)
    {
        t_swapGroup *g = &swap->grp[i];

        block_bc(cr, *g);
        snew_bc(cr, g->ind, g->nat);
        nblock_bc(cr, g->nat, g->ind);

        int len = 0;
        if (MASTER(cr))
        {
            len = strlen(g->molname);
        }
        block_bc(cr, len);
        snew_bc(cr, g->molname, len);
        nblock_bc(cr, len, g->molname);
    }
}


static void bc_inputrec(const t_commrec *cr, t_inputrec *inputrec)
{
    // Note that this overwrites pointers in inputrec, so all pointer fields
    // Must be initialized separately below.
    block_bc(cr, *inputrec);
    if (SIMMASTER(cr))
    {
        gmx::InMemorySerializer serializer;
        gmx::serializeKeyValueTree(*inputrec->params, &serializer);
        std::vector<char>       buffer = serializer.finishAndGetBuffer();
        size_t                  size   = buffer.size();
        block_bc(cr, size);
        nblock_bc(cr, size, buffer.data());
    }
    else
    {
        // block_bc() above overwrites the old pointer, so set it to a
        // reasonable value in case code below throws.
        // cppcheck-suppress redundantAssignment
        inputrec->params = nullptr;
        std::vector<char> buffer;
        size_t            size;
        block_bc(cr, size);
        nblock_abc(cr, size, &buffer);
        gmx::InMemoryDeserializer serializer(buffer);
        // cppcheck-suppress redundantAssignment
        inputrec->params = new gmx::KeyValueTreeObject(
                    gmx::deserializeKeyValueTree(&serializer));
    }

    bc_grpopts(cr, &(inputrec->opts));

    /* even if efep is efepNO, we need to initialize to make sure that
     * n_lambda is set to zero */

    snew_bc(cr, inputrec->fepvals, 1);
    if (inputrec->efep != efepNO || inputrec->bSimTemp)
    {
        bc_fepvals(cr, inputrec->fepvals);
    }
    /* need to initialize this as well because of data checked for in the logic */
    snew_bc(cr, inputrec->expandedvals, 1);
    if (inputrec->bExpanded)
    {
        bc_expandedvals(cr, inputrec->expandedvals, inputrec->fepvals->n_lambda);
    }
    snew_bc(cr, inputrec->simtempvals, 1);
    if (inputrec->bSimTemp)
    {
        bc_simtempvals(cr, inputrec->simtempvals, inputrec->fepvals->n_lambda);
    }
    if (inputrec->bPull)
    {
        snew_bc(cr, inputrec->pull, 1);
        bc_pull(cr, inputrec->pull);
    }
    if (inputrec->bDoAwh)
    {
        snew_bc(cr, inputrec->awhParams, 1);
        bc_awh(cr, inputrec->awhParams);
    }

    if (inputrec->bRot)
    {
        snew_bc(cr, inputrec->rot, 1);
        bc_rot(cr, inputrec->rot);
    }
    if (inputrec->bIMD)
    {
        snew_bc(cr, inputrec->imd, 1);
        bc_imd(cr, inputrec->imd);
    }
    if (inputrec->eSwapCoords != eswapNO)
    {
        snew_bc(cr, inputrec->swap, 1);
        bc_swapions(cr, inputrec->swap);
    }
}

static void bc_moltype(const t_commrec *cr, t_symtab *symtab,
                       gmx_moltype_t *moltype)
{
    bc_string(cr, symtab, &moltype->name);
    bc_atoms(cr, symtab, &moltype->atoms);
    if (debug)
    {
        fprintf(debug, "after bc_atoms\n");
    }

    bc_ilists(cr, moltype->ilist);
    bc_block(cr, &moltype->cgs);
    bc_blocka(cr, &moltype->excls);
}

static void bc_molblock(const t_commrec *cr, gmx_molblock_t *molb)
{
    block_bc(cr, *molb);
    if (molb->nposres_xA > 0)
    {
        snew_bc(cr, molb->posres_xA, molb->nposres_xA);
        nblock_bc(cr, molb->nposres_xA*DIM, molb->posres_xA[0]);
    }
    if (molb->nposres_xB > 0)
    {
        snew_bc(cr, molb->posres_xB, molb->nposres_xB);
        nblock_bc(cr, molb->nposres_xB*DIM, molb->posres_xB[0]);
    }
    if (debug)
    {
        fprintf(debug, "after bc_molblock\n");
    }
}

static void bc_atomtypes(const t_commrec *cr, t_atomtypes *atomtypes)
{
    int nr;

    block_bc(cr, atomtypes->nr);

    nr = atomtypes->nr;

    snew_bc(cr, atomtypes->radius, nr);
    snew_bc(cr, atomtypes->vol, nr);
    snew_bc(cr, atomtypes->surftens, nr);
    snew_bc(cr, atomtypes->gb_radius, nr);
    snew_bc(cr, atomtypes->S_hct, nr);

    nblock_bc(cr, nr, atomtypes->radius);
    nblock_bc(cr, nr, atomtypes->vol);
    nblock_bc(cr, nr, atomtypes->surftens);
    nblock_bc(cr, nr, atomtypes->gb_radius);
    nblock_bc(cr, nr, atomtypes->S_hct);
}

/*! \brief Broadcasts ir and mtop from the master to all nodes in
 * cr->mpi_comm_mygroup. */
static
void bcast_ir_mtop(const t_commrec *cr, t_inputrec *inputrec, gmx_mtop_t *mtop)
{
    int i;
    if (debug)
    {
        fprintf(debug, "in bc_data\n");
    }
    bc_inputrec(cr, inputrec);
    if (debug)
    {
        fprintf(debug, "after bc_inputrec\n");
    }
    bc_symtab(cr, &mtop->symtab);
    if (debug)
    {
        fprintf(debug, "after bc_symtab\n");
    }
    bc_string(cr, &mtop->symtab, &mtop->name);
    if (debug)
    {
        fprintf(debug, "after bc_name\n");
    }

    bc_ffparams(cr, &mtop->ffparams);

    block_bc(cr, mtop->nmoltype);
    snew_bc(cr, mtop->moltype, mtop->nmoltype);
    for (i = 0; i < mtop->nmoltype; i++)
    {
        bc_moltype(cr, &mtop->symtab, &mtop->moltype[i]);
    }

    block_bc(cr, mtop->bIntermolecularInteractions);
    if (mtop->bIntermolecularInteractions)
    {
        snew_bc(cr, mtop->intermolecular_ilist, F_NRE);
        bc_ilists(cr, mtop->intermolecular_ilist);
    }

    block_bc(cr, mtop->nmolblock);
    snew_bc(cr, mtop->molblock, mtop->nmolblock);
    for (i = 0; i < mtop->nmolblock; i++)
    {
        bc_molblock(cr, &mtop->molblock[i]);
    }

    block_bc(cr, mtop->natoms);

    bc_atomtypes(cr, &mtop->atomtypes);

    bc_block(cr, &mtop->mols);
    bc_groups(cr, &mtop->symtab, mtop->natoms, &mtop->groups);
}

void init_parallel(t_commrec *cr, t_inputrec *inputrec,
                   gmx_mtop_t *mtop)
{
    bcast_ir_mtop(cr, inputrec, mtop);
}
