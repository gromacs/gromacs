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
#include "gmxpre.h"

#include "topology.h"

#include <cstdio>

#include <algorithm>

#include "gromacs/math/vecdump.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/txtdump.h"

const char *gtypes[egcNR+1] = {
    "T-Coupling", "Energy Mon.", "Acceleration", "Freeze",
    "User1", "User2", "VCM", "Compressed X", "Or. Res. Fit", "QMMM", nullptr
};

static void init_groups(gmx_groups_t *groups)
{
    groups->ngrpname = 0;
    groups->grpname  = nullptr;
    for (int g = 0; g < egcNR; g++)
    {
        groups->grps[g].nm_ind = nullptr;
        groups->ngrpnr[g]      = 0;
        groups->grpnr[g]       = nullptr;
    }

}

void init_mtop(gmx_mtop_t *mtop)
{
    mtop->name         = nullptr;
    mtop->nmoltype     = 0;
    mtop->moltype      = nullptr;
    mtop->nmolblock    = 0;
    mtop->molblock     = nullptr;
    mtop->maxres_renum = 0;
    mtop->maxresnr     = -1;
    init_groups(&mtop->groups);
    init_block(&mtop->mols);
    open_symtab(&mtop->symtab);
}

void init_top(t_topology *top)
{
    top->name = nullptr;
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

void done_gmx_groups_t(gmx_groups_t *g)
{
    int i;

    for (i = 0; (i < egcNR); i++)
    {
        if (nullptr != g->grps[i].nm_ind)
        {
            sfree(g->grps[i].nm_ind);
            g->grps[i].nm_ind = nullptr;
        }
        if (nullptr != g->grpnr[i])
        {
            sfree(g->grpnr[i]);
            g->grpnr[i] = nullptr;
        }
    }
    /* The contents of this array is in symtab, don't free it here */
    sfree(g->grpname);
}

void done_mtop(gmx_mtop_t *mtop)
{
    done_symtab(&mtop->symtab);

    sfree(mtop->ffparams.functype);
    sfree(mtop->ffparams.iparams);
    for (int i = 0; i < mtop->ffparams.cmap_grid.ngrid; i++)
    {
        sfree(mtop->ffparams.cmap_grid.cmapdata[i].cmap);
    }
    sfree(mtop->ffparams.cmap_grid.cmapdata);

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
    done_atomtypes(&mtop->atomtypes);
    done_gmx_groups_t(&mtop->groups);
    done_block(&mtop->mols);
}

void done_top(t_topology *top)
{
    sfree(top->idef.functype);
    sfree(top->idef.iparams);
    for (int f = 0; f < F_NRE; ++f)
    {
        sfree(top->idef.il[f].iatoms);
        top->idef.il[f].iatoms = nullptr;
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

void done_top_mtop(t_topology *top, gmx_mtop_t *mtop)
{
    if (mtop != nullptr)
    {
        if (top != nullptr)
        {
            for (int f = 0; f < F_NRE; ++f)
            {
                sfree(top->idef.il[f].iatoms);
                top->idef.il[f].iatoms = nullptr;
                top->idef.il[f].nalloc = 0;
            }
            done_atom(&top->atoms);
            done_block(&top->cgs);
            done_blocka(&top->excls);
            done_symtab(&top->symtab);
            open_symtab(&mtop->symtab);
        }
        done_mtop(mtop);
    }
}

bool gmx_mtop_has_masses(const gmx_mtop_t *mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->nmoltype == 0 || mtop->moltype[0].atoms.haveMass;
}

bool gmx_mtop_has_charges(const gmx_mtop_t *mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->nmoltype == 0 || mtop->moltype[0].atoms.haveCharge;
}

bool gmx_mtop_has_atomtypes(const gmx_mtop_t *mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->nmoltype == 0 || mtop->moltype[0].atoms.haveType;
}

bool gmx_mtop_has_pdbinfo(const gmx_mtop_t *mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->nmoltype == 0 || mtop->moltype[0].atoms.havePdbInfo;
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
                       gmx_bool bShowNumbers, gmx_bool bShowParameters)
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
                 ffparams->functype, &molt->ilist[j],
                 bShowNumbers, bShowParameters, ffparams->iparams);
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
             gmx_bool bShowNumbers, gmx_bool bShowParameters)
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
                         &mtop->intermolecular_ilist[j],
                         bShowNumbers, bShowParameters, mtop->ffparams.iparams);
            }
        }
        pr_ffparams(fp, indent, "ffparams", &(mtop->ffparams), bShowNumbers);
        pr_atomtypes(fp, indent, "atomtypes", &(mtop->atomtypes), bShowNumbers);
        for (mt = 0; mt < mtop->nmoltype; mt++)
        {
            pr_moltype(fp, indent, "moltype", &mtop->moltype[mt], mt,
                       &mtop->ffparams, bShowNumbers, bShowParameters);
        }
        pr_groups(fp, indent, &mtop->groups, bShowNumbers);
    }
}

void pr_top(FILE *fp, int indent, const char *title, const t_topology *top,
            gmx_bool bShowNumbers, gmx_bool bShowParameters)
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
        pr_idef(fp, indent, "idef", &top->idef, bShowNumbers, bShowParameters);
    }
}

static void cmp_ilist(FILE *fp, int ftype, const t_ilist *il1, const t_ilist *il2)
{
    int  i;
    char buf[256];

    fprintf(fp, "comparing ilist %s\n", interaction_function[ftype].name);
    sprintf(buf, "%s->nr", interaction_function[ftype].name);
    cmp_int(fp, buf, -1, il1->nr, il2->nr);
    sprintf(buf, "%s->iatoms", interaction_function[ftype].name);
    if (((il1->nr > 0) && (!il1->iatoms)) ||
        ((il2->nr > 0) && (!il2->iatoms)) ||
        ((il1->nr != il2->nr)))
    {
        fprintf(fp, "Comparing radically different topologies - %s is different\n",
                buf);
    }
    else
    {
        for (i = 0; (i < il1->nr); i++)
        {
            cmp_int(fp, buf, i, il1->iatoms[i], il2->iatoms[i]);
        }
    }
}

static void cmp_iparm(FILE *fp, const char *s, t_functype ft,
                      const t_iparams &ip1, const t_iparams &ip2, real ftol, real abstol)
{
    int      i;
    gmx_bool bDiff;

    bDiff = FALSE;
    for (i = 0; i < MAXFORCEPARAM && !bDiff; i++)
    {
        bDiff = !equal_real(ip1.generic.buf[i], ip2.generic.buf[i], ftol, abstol);
    }
    if (bDiff)
    {
        fprintf(fp, "%s1: ", s);
        pr_iparams(fp, ft, &ip1);
        fprintf(fp, "%s2: ", s);
        pr_iparams(fp, ft, &ip2);
    }
}

static void cmp_iparm_AB(FILE *fp, const char *s, t_functype ft,
                         const t_iparams &ip1, real ftol, real abstol)
{
    int      nrfpA, nrfpB, p0, i;
    gmx_bool bDiff;

    /* Normally the first parameter is perturbable */
    p0    = 0;
    nrfpA = interaction_function[ft].nrfpA;
    nrfpB = interaction_function[ft].nrfpB;
    if (ft == F_PDIHS)
    {
        nrfpB = 2;
    }
    else if (interaction_function[ft].flags & IF_TABULATED)
    {
        /* For tabulated interactions only the second parameter is perturbable */
        p0    = 1;
        nrfpB = 1;
    }
    bDiff = FALSE;
    for (i = 0; i < nrfpB && !bDiff; i++)
    {
        bDiff = !equal_real(ip1.generic.buf[p0+i], ip1.generic.buf[nrfpA+i], ftol, abstol);
    }
    if (bDiff)
    {
        fprintf(fp, "%s: ", s);
        pr_iparams(fp, ft, &ip1);
    }
}

static void cmp_cmap(FILE *fp, const gmx_cmap_t *cmap1, const gmx_cmap_t *cmap2, real ftol, real abstol)
{
    cmp_int(fp, "cmap ngrid", -1, cmap1->ngrid, cmap2->ngrid);
    cmp_int(fp, "cmap grid_spacing", -1, cmap1->grid_spacing, cmap2->grid_spacing);
    if (cmap1->ngrid == cmap2->ngrid &&
        cmap1->grid_spacing == cmap2->grid_spacing)
    {
        int g;

        for (g = 0; g < cmap1->ngrid; g++)
        {
            int i;

            fprintf(fp, "comparing cmap %d\n", g);

            for (i = 0; i < 4*cmap1->grid_spacing*cmap1->grid_spacing; i++)
            {
                cmp_real(fp, "", i, cmap1->cmapdata[g].cmap[i], cmap2->cmapdata[g].cmap[i], ftol, abstol);
            }
        }
    }
}

static void cmp_idef(FILE *fp, const t_idef *id1, const t_idef *id2, real ftol, real abstol)
{
    int  i;
    char buf1[64], buf2[64];

    fprintf(fp, "comparing idef\n");
    if (id2)
    {
        cmp_int(fp, "idef->ntypes", -1, id1->ntypes, id2->ntypes);
        cmp_int(fp, "idef->atnr",  -1, id1->atnr, id2->atnr);
        for (i = 0; (i < std::min(id1->ntypes, id2->ntypes)); i++)
        {
            sprintf(buf1, "idef->functype[%d]", i);
            sprintf(buf2, "idef->iparam[%d]", i);
            cmp_int(fp, buf1, i, (int)id1->functype[i], (int)id2->functype[i]);
            cmp_iparm(fp, buf2, id1->functype[i],
                      id1->iparams[i], id2->iparams[i], ftol, abstol);
        }
        cmp_real(fp, "fudgeQQ", -1, id1->fudgeQQ, id2->fudgeQQ, ftol, abstol);
        cmp_cmap(fp, &id1->cmap_grid, &id2->cmap_grid, ftol, abstol);
        for (i = 0; (i < F_NRE); i++)
        {
            cmp_ilist(fp, i, &(id1->il[i]), &(id2->il[i]));
        }
    }
    else
    {
        for (i = 0; (i < id1->ntypes); i++)
        {
            cmp_iparm_AB(fp, "idef->iparam", id1->functype[i], id1->iparams[i], ftol, abstol);
        }
    }
}

static void cmp_block(FILE *fp, const t_block *b1, const t_block *b2, const char *s)
{
    char buf[32];

    fprintf(fp, "comparing block %s\n", s);
    sprintf(buf, "%s.nr", s);
    cmp_int(fp, buf, -1, b1->nr, b2->nr);
}

static void cmp_blocka(FILE *fp, const t_blocka *b1, const t_blocka *b2, const char *s)
{
    char buf[32];

    fprintf(fp, "comparing blocka %s\n", s);
    sprintf(buf, "%s.nr", s);
    cmp_int(fp, buf, -1, b1->nr, b2->nr);
    sprintf(buf, "%s.nra", s);
    cmp_int(fp, buf, -1, b1->nra, b2->nra);
}

void cmp_top(FILE *fp, const t_topology *t1, const t_topology *t2, real ftol, real abstol)
{
    fprintf(fp, "comparing top\n");
    if (t2)
    {
        cmp_idef(fp, &(t1->idef), &(t2->idef), ftol, abstol);
        cmp_atoms(fp, &(t1->atoms), &(t2->atoms), ftol, abstol);
        cmp_block(fp, &t1->cgs, &t2->cgs, "cgs");
        cmp_block(fp, &t1->mols, &t2->mols, "mols");
        cmp_bool(fp, "bIntermolecularInteractions", -1, t1->bIntermolecularInteractions, t2->bIntermolecularInteractions);
        cmp_blocka(fp, &t1->excls, &t2->excls, "excls");
    }
    else
    {
        cmp_idef(fp, &(t1->idef), nullptr, ftol, abstol);
        cmp_atoms(fp, &(t1->atoms), nullptr, ftol, abstol);
    }
}

void cmp_groups(FILE *fp, const gmx_groups_t *g0, const gmx_groups_t *g1,
                int natoms0, int natoms1)
{
    int  i, j;
    char buf[32];

    fprintf(fp, "comparing groups\n");

    for (i = 0; i < egcNR; i++)
    {
        sprintf(buf, "grps[%d].nr", i);
        cmp_int(fp, buf, -1, g0->grps[i].nr, g1->grps[i].nr);
        if (g0->grps[i].nr == g1->grps[i].nr)
        {
            for (j = 0; j < g0->grps[i].nr; j++)
            {
                sprintf(buf, "grps[%d].name[%d]", i, j);
                cmp_str(fp, buf, -1,
                        *g0->grpname[g0->grps[i].nm_ind[j]],
                        *g1->grpname[g1->grps[i].nm_ind[j]]);
            }
        }
        cmp_int(fp, "ngrpnr", i, g0->ngrpnr[i], g1->ngrpnr[i]);
        if (g0->ngrpnr[i] == g1->ngrpnr[i] && natoms0 == natoms1 &&
            (g0->grpnr[i] != nullptr || g1->grpnr[i] != nullptr))
        {
            for (j = 0; j < natoms0; j++)
            {
                cmp_int(fp, gtypes[i], j, ggrpnr(g0, i, j), ggrpnr(g1, i, j));
            }
        }
    }
    /* We have compared the names in the groups lists,
     * so we can skip the grpname list comparison.
     */
}
