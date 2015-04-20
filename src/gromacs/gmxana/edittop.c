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

#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

void replace_atom(t_topology *top, int inr, char *anm, char *resnm,
                  real q, real m, int type)
{
    t_atoms *atoms;

    atoms = &(top->atoms);

    /* Replace important properties of an atom by other properties */
    if ((inr < 0) || (inr > atoms->nr))
    {
        gmx_fatal(FARGS, "Replace_atom: inr (%d) not in %d .. %d", inr, 0, atoms->nr);
    }
    if (debug)
    {
        fprintf(debug, "Replacing atom %d ... ", inr);
    }
    /* Charge, mass and type */
    atoms->atom[inr].q    = atoms->atom[inr].qB    = q;
    atoms->atom[inr].m    = atoms->atom[inr].mB    = m;
    atoms->atom[inr].type = atoms->atom[inr].typeB = type;

    /* Residue name */
    atoms->resinfo[atoms->atom[inr].resind].name = put_symtab(&top->symtab, resnm);
    /* Atom name */
    atoms->atomname[inr] = put_symtab(&top->symtab, anm);
    if (debug)
    {
        fprintf(debug, " done\n");
    }
}

static void delete_from_interactions(t_idef *idef, int inr)
{
    int      i, j, k, nra, nnr;
    t_iatom *niatoms;
    gmx_bool bDel;

    /* Delete interactions including atom inr from lists */
    for (i = 0; (i < F_NRE); i++)
    {
        nra = interaction_function[i].nratoms;
        nnr = 0;
        snew(niatoms, idef->il[i].nr);
        for (j = 0; (j < idef->il[i].nr); j += nra+1)
        {
            bDel = FALSE;
            for (k = 0; (k < nra); k++)
            {
                if (idef->il[i].iatoms[j+k+1] == inr)
                {
                    bDel = TRUE;
                }
            }
            if (!bDel)
            {
                /* If this does not need to be deleted, then copy it to temp array */
                for (k = 0; (k < nra+1); k++)
                {
                    niatoms[nnr+k] = idef->il[i].iatoms[j+k];
                }
                nnr += nra+1;
            }
        }
        /* Copy temp array back */
        for (j = 0; (j < nnr); j++)
        {
            idef->il[i].iatoms[j] = niatoms[j];
        }
        idef->il[i].nr = nnr;
        sfree(niatoms);
    }
}

static void delete_from_block(t_block *block, int inr)
{
    /* Update block data structure */
    int i, i1, j1, j, k;

    for (i = 0; (i < block->nr); i++)
    {
        for (j = block->index[i]; (j < block->index[i+1]); j++)
        {
            if (j == inr)
            {
                /* This atom has to go */
                /* Change indices too */
                for (i1 = i+1; (i1 <= block->nr); i1++)
                {
                    block->index[i1]--;
                }
            }
        }
    }
}

static void delete_from_blocka(t_blocka *block, int inr)
{
    /* Update block data structure */
    int i, i1, j1, j, k;

    for (i = 0; (i < block->nr); i++)
    {
        for (j = block->index[i]; (j < block->index[i+1]); j++)
        {
            k = block->a[j];
            if (k == inr)
            {
                /* This atom has to go */
                for (j1 = j; (j1 < block->nra-1); j1++)
                {
                    block->a[j1] = block->a[j1+1];
                }
                block->nra--;
                /* Change indices too */
                for (i1 = i+1; (i1 <= block->nr); i1++)
                {
                    block->index[i1]--;
                }
            }
        }
    }
}

static void delete_from_atoms(t_atoms *atoms, int inr)
{
    int i;

    /* Shift the atomnames down */
    for (i = inr; (i < atoms->nr-1); i++)
    {
        atoms->atomname[i] = atoms->atomname[i+1];
    }

    /* Shift the atom struct down */
    for (i = inr; (i < atoms->nr-1); i++)
    {
        atoms->atom[i] = atoms->atom[i+1];
    }

    if (atoms->pdbinfo)
    {
        /* Shift the pdbatom struct down */
        for (i = inr; (i < atoms->nr-1); i++)
        {
            atoms->pdbinfo[i] = atoms->pdbinfo[i+1];
        }
    }
    atoms->nr--;
}

void delete_atom(t_topology *top, int inr)
{
    int k;

    if ((inr < 0) || (inr >= top->atoms.nr))
    {
        gmx_fatal(FARGS, "Delete_atom: inr (%d) not in %d .. %d", inr, 0,
                  top->atoms.nr);
    }
    if (debug)
    {
        fprintf(debug, "Deleting atom %d ...", inr);
    }

    /* First remove bonds etc. */
    delete_from_interactions(&top->idef, inr);
    /* Now charge groups etc. */
    delete_from_block(&(top->cgs), inr);
    delete_from_block(&(top->mols), inr);
    delete_from_blocka(&(top->excls), inr);
    /* Now from the atoms struct */
    delete_from_atoms(&top->atoms, inr);
    if (debug)
    {
        fprintf(debug, " done\n");
    }
}
