/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
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

#include "topology.h"

#include "gromacs/topology/symtab.h"
#include "gromacs/utility/smalloc.h"

static void init_groups(gmx_groups_t *groups)
{
    groups->ngrpname = 0;
    groups->grpname  = NULL;
    for (int g = 0; g < egcNR; g++)
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

void init_top(t_topology *top)
{
    top->name = NULL;
    init_atom(&(top->atoms));
    init_atomtypes(&(top->atomtypes));
    init_block(&top->cgs);
    init_block(&top->mols);
    init_blocka(&top->excls);
    open_symtab(&top->symtab);
}


void done_moltype(gmx_moltype_t *molt)
{
    done_atom(&molt->atoms);
    done_block(&molt->cgs);
    done_blocka(&molt->excls);

    for (int f = 0; f < F_NRE; f++)
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
        sfree(molb->posres_xA);
    }
    if (molb->nposres_xB > 0)
    {
        molb->nposres_xB = 0;
        sfree(molb->posres_xB);
    }
}

void done_mtop(gmx_mtop_t *mtop, gmx_bool bDoneSymtab)
{
    if (bDoneSymtab)
    {
        done_symtab(&mtop->symtab);
    }

    sfree(mtop->ffparams.functype);
    sfree(mtop->ffparams.iparams);

    for (int i = 0; i < mtop->nmoltype; i++)
    {
        done_moltype(&mtop->moltype[i]);
    }
    sfree(mtop->moltype);
    for (int i = 0; i < mtop->nmolblock; i++)
    {
        done_molblock(&mtop->molblock[i]);
    }
    sfree(mtop->molblock);
    done_block(&mtop->mols);
}

void done_top(t_topology *top)
{
    sfree(top->idef.functype);
    sfree(top->idef.iparams);
    for (int f = 0; f < F_NRE; ++f)
    {
        sfree(top->idef.il[f].iatoms);
        top->idef.il[f].iatoms = NULL;
        top->idef.il[f].nalloc = 0;
    }

    done_atom(&(top->atoms));

    /* For GB */
    done_atomtypes(&(top->atomtypes));

    done_symtab(&(top->symtab));
    done_block(&(top->cgs));
    done_block(&(top->mols));
    done_blocka(&(top->excls));
}
