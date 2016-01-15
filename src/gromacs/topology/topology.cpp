/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
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

#include <cstdio>

#include <algorithm>

#include "gromacs/math/vecdump.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/txtdump.h"

const char *gtypes[egcNR+1] = {
    "T-Coupling", "Energy Mon.", "Acceleration", "Freeze",
    "User1", "User2", "VCM", "Compressed X", "Or. Res. Fit", "QMMM", NULL
};

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

static void pr_grps(FILE *fp, const char *title, const t_grps grps[], char **grpname[])
{
    int i, j;

    for (i = 0; (i < egcNR); i++)
    {
        fprintf(fp, "%s[%-12s] nr=%d, name=[", title, gtypes[i], grps[i].nr);
        for (j = 0; (j < grps[i].nr); j++)
        {
            fprintf(fp, " %s", *(grpname[grps[i].nm_ind[j]]));
        }
        fprintf(fp, "]\n");
    }
}

static void pr_groups(FILE *fp, int indent,
                      const gmx_groups_t *groups,
                      gmx_bool bShowNumbers)
{
    int nat_max, i, g;

    pr_grps(fp, "grp", groups->grps, groups->grpname);
    pr_strings(fp, indent, "grpname", groups->grpname, groups->ngrpname, bShowNumbers);

    pr_indent(fp, indent);
    fprintf(fp, "groups          ");
    for (g = 0; g < egcNR; g++)
    {
        printf(" %5.5s", gtypes[g]);
    }
    printf("\n");

    pr_indent(fp, indent);
    fprintf(fp, "allocated       ");
    nat_max = 0;
    for (g = 0; g < egcNR; g++)
    {
        printf(" %5d", groups->ngrpnr[g]);
        nat_max = std::max(nat_max, groups->ngrpnr[g]);
    }
    printf("\n");

    if (nat_max == 0)
    {
        pr_indent(fp, indent);
        fprintf(fp, "groupnr[%5s] =", "*");
        for (g = 0; g < egcNR; g++)
        {
            fprintf(fp, "  %3d ", 0);
        }
        fprintf(fp, "\n");
    }
    else
    {
        for (i = 0; i < nat_max; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "groupnr[%5d] =", i);
            for (g = 0; g < egcNR; g++)
            {
                fprintf(fp, "  %3d ",
                        groups->grpnr[g] ? groups->grpnr[g][i] : 0);
            }
            fprintf(fp, "\n");
        }
    }
}

static void pr_moltype(FILE *fp, int indent, const char *title,
                       const gmx_moltype_t *molt, int n,
                       const gmx_ffparams_t *ffparams,
                       gmx_bool bShowNumbers)
{
    int j;

    indent = pr_title_n(fp, indent, title, n);
    pr_indent(fp, indent);
    fprintf(fp, "name=\"%s\"\n", *(molt->name));
    pr_atoms(fp, indent, "atoms", &(molt->atoms), bShowNumbers);
    pr_block(fp, indent, "cgs", &molt->cgs, bShowNumbers);
    pr_blocka(fp, indent, "excls", &molt->excls, bShowNumbers);
    for (j = 0; (j < F_NRE); j++)
    {
        pr_ilist(fp, indent, interaction_function[j].longname,
                 ffparams->functype, &molt->ilist[j], bShowNumbers);
    }
}

static void pr_molblock(FILE *fp, int indent, const char *title,
                        const gmx_molblock_t *molb, int n,
                        const gmx_moltype_t *molt)
{
    indent = pr_title_n(fp, indent, title, n);
    pr_indent(fp, indent);
    fprintf(fp, "%-20s = %d \"%s\"\n",
            "moltype", molb->type, *(molt[molb->type].name));
    pr_int(fp, indent, "#molecules", molb->nmol);
    pr_int(fp, indent, "#atoms_mol", molb->natoms_mol);
    pr_int(fp, indent, "#posres_xA", molb->nposres_xA);
    if (molb->nposres_xA > 0)
    {
        pr_rvecs(fp, indent, "posres_xA", molb->posres_xA, molb->nposres_xA);
    }
    pr_int(fp, indent, "#posres_xB", molb->nposres_xB);
    if (molb->nposres_xB > 0)
    {
        pr_rvecs(fp, indent, "posres_xB", molb->posres_xB, molb->nposres_xB);
    }
}

void pr_mtop(FILE *fp, int indent, const char *title, const gmx_mtop_t *mtop,
             gmx_bool bShowNumbers)
{
    int mt, mb, j;

    if (available(fp, mtop, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "name=\"%s\"\n", *(mtop->name));
        pr_int(fp, indent, "#atoms", mtop->natoms);
        pr_int(fp, indent, "#molblock", mtop->nmolblock);
        for (mb = 0; mb < mtop->nmolblock; mb++)
        {
            pr_molblock(fp, indent, "molblock", &mtop->molblock[mb], mb, mtop->moltype);
        }
        pr_str(fp, indent, "bIntermolecularInteractions",
               gmx::boolToString(mtop->bIntermolecularInteractions));
        if (mtop->bIntermolecularInteractions)
        {
            for (j = 0; (j < F_NRE); j++)
            {
                pr_ilist(fp, indent, interaction_function[j].longname,
                         mtop->ffparams.functype,
                         &mtop->intermolecular_ilist[j], bShowNumbers);
            }
        }
        pr_ffparams(fp, indent, "ffparams", &(mtop->ffparams), bShowNumbers);
        pr_atomtypes(fp, indent, "atomtypes", &(mtop->atomtypes), bShowNumbers);
        for (mt = 0; mt < mtop->nmoltype; mt++)
        {
            pr_moltype(fp, indent, "moltype", &mtop->moltype[mt], mt,
                       &mtop->ffparams, bShowNumbers);
        }
        pr_groups(fp, indent, &mtop->groups, bShowNumbers);
    }
}

void pr_top(FILE *fp, int indent, const char *title, const t_topology *top, gmx_bool bShowNumbers)
{
    if (available(fp, top, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "name=\"%s\"\n", *(top->name));
        pr_atoms(fp, indent, "atoms", &(top->atoms), bShowNumbers);
        pr_atomtypes(fp, indent, "atomtypes", &(top->atomtypes), bShowNumbers);
        pr_block(fp, indent, "cgs", &top->cgs, bShowNumbers);
        pr_block(fp, indent, "mols", &top->mols, bShowNumbers);
        pr_str(fp, indent, "bIntermolecularInteractions",
               gmx::boolToString(top->bIntermolecularInteractions));
        pr_blocka(fp, indent, "excls", &top->excls, bShowNumbers);
        pr_idef(fp, indent, "idef", &top->idef, bShowNumbers);
    }
}
