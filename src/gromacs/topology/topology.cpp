/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/gmxassert.h"
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
        groups->grps[g].nr     = 0;
        groups->grps[g].nm_ind = nullptr;
        groups->ngrpnr[g]      = 0;
        groups->grpnr[g]       = nullptr;
    }

}

void init_top(t_topology *top)
{
    top->name = nullptr;
    init_idef(&top->idef);
    init_atom(&(top->atoms));
    init_atomtypes(&(top->atomtypes));
    init_block(&top->cgs);
    init_block(&top->mols);
    init_blocka(&top->excls);
    open_symtab(&top->symtab);
}


gmx_moltype_t::gmx_moltype_t() :
    name(nullptr),
    cgs(),
    excls()
{
    init_t_atoms(&atoms, 0, FALSE);
}

gmx_moltype_t::~gmx_moltype_t()
{
    done_atom(&atoms);
    done_block(&cgs);
    done_blocka(&excls);
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

gmx_mtop_t::gmx_mtop_t()
{
    init_atomtypes(&atomtypes);
    init_groups(&groups);
    open_symtab(&symtab);
}

gmx_mtop_t::~gmx_mtop_t()
{
    done_symtab(&symtab);

    moltype.clear();
    molblock.clear();
    done_atomtypes(&atomtypes);
    done_gmx_groups_t(&groups);
}

void done_top(t_topology *top)
{
    done_idef(&top->idef);
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
            done_idef(&top->idef);
            done_atom(&top->atoms);
            done_block(&top->cgs);
            done_blocka(&top->excls);
            done_block(&top->mols);
            done_symtab(&top->symtab);
            open_symtab(&mtop->symtab);
            done_atomtypes(&(top->atomtypes));
        }

        // Note that the rest of mtop will be freed by the destructor
    }
}

void init_localtop(gmx_localtop_t *top)
{
    init_block(&top->cgs);
    init_blocka(&top->excls);
    init_idef(&top->idef);
    init_atomtypes(&top->atomtypes);
}

void done_localtop(gmx_localtop_t *top)
{
    if (top == nullptr)
    {
        return;
    }
    done_idef(&top->idef);
    done_block(&top->cgs);
    done_blocka(&top->excls);
    done_atomtypes(&top->atomtypes);
}

void done_and_sfree_localtop(gmx_localtop_t *top)
{
    done_localtop(top);
    sfree(top);
}

bool gmx_mtop_has_masses(const gmx_mtop_t *mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->moltype.empty() || mtop->moltype[0].atoms.haveMass;
}

bool gmx_mtop_has_charges(const gmx_mtop_t *mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->moltype.empty() || mtop->moltype[0].atoms.haveCharge;
}

bool gmx_mtop_has_perturbed_charges(const gmx_mtop_t &mtop)
{
    for (const gmx_moltype_t &moltype : mtop.moltype)
    {
        const t_atoms &atoms = moltype.atoms;
        if (atoms.haveBState)
        {
            for (int a = 0; a < atoms.nr; a++)
            {
                if (atoms.atom[a].q != atoms.atom[a].qB)
                {
                    return true;
                }
            }
        }
    }
    return false;
}

bool gmx_mtop_has_atomtypes(const gmx_mtop_t *mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->moltype.empty() || mtop->moltype[0].atoms.haveType;
}

bool gmx_mtop_has_pdbinfo(const gmx_mtop_t *mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->moltype.empty() || mtop->moltype[0].atoms.havePdbInfo;
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
                 ffparams->functype.data(), molt->ilist[j],
                 bShowNumbers, bShowParameters, ffparams->iparams.data());
    }
}

static void pr_molblock(FILE *fp, int indent, const char *title,
                        const gmx_molblock_t *molb, int n,
                        const std::vector<gmx_moltype_t> &molt)
{
    indent = pr_title_n(fp, indent, title, n);
    pr_indent(fp, indent);
    fprintf(fp, "%-20s = %d \"%s\"\n",
            "moltype", molb->type, *(molt[molb->type].name));
    pr_int(fp, indent, "#molecules", molb->nmol);
    pr_int(fp, indent, "#posres_xA", molb->posres_xA.size());
    if (!molb->posres_xA.empty())
    {
        pr_rvecs(fp, indent, "posres_xA", as_rvec_array(molb->posres_xA.data()), molb->posres_xA.size());
    }
    pr_int(fp, indent, "#posres_xB", molb->posres_xB.size());
    if (!molb->posres_xB.empty())
    {
        pr_rvecs(fp, indent, "posres_xB", as_rvec_array(molb->posres_xB.data()), molb->posres_xB.size());
    }
}

void pr_mtop(FILE *fp, int indent, const char *title, const gmx_mtop_t *mtop,
             gmx_bool bShowNumbers, gmx_bool bShowParameters)
{
    if (available(fp, mtop, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "name=\"%s\"\n", *(mtop->name));
        pr_int(fp, indent, "#atoms", mtop->natoms);
        pr_int(fp, indent, "#molblock", mtop->molblock.size());
        for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
        {
            pr_molblock(fp, indent, "molblock", &mtop->molblock[mb], mb, mtop->moltype);
        }
        pr_str(fp, indent, "bIntermolecularInteractions",
               gmx::boolToString(mtop->bIntermolecularInteractions));
        if (mtop->bIntermolecularInteractions)
        {
            for (int j = 0; j < F_NRE; j++)
            {
                pr_ilist(fp, indent, interaction_function[j].longname,
                         mtop->ffparams.functype.data(),
                         (*mtop->intermolecular_ilist)[j],
                         bShowNumbers, bShowParameters, mtop->ffparams.iparams.data());
            }
        }
        pr_ffparams(fp, indent, "ffparams", &(mtop->ffparams), bShowNumbers);
        pr_atomtypes(fp, indent, "atomtypes", &(mtop->atomtypes), bShowNumbers);
        for (size_t mt = 0; mt < mtop->moltype.size(); mt++)
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
    int cmap1_ngrid = (cmap1 ? cmap1->cmapdata.size() : 0);
    int cmap2_ngrid = (cmap2 ? cmap2->cmapdata.size() : 0);

    cmp_int(fp, "cmap ngrid", -1, cmap1_ngrid, cmap2_ngrid);

    if (cmap1 == nullptr || cmap2 == nullptr)
    {
        return;
    }

    cmp_int(fp, "cmap grid_spacing", -1, cmap1->grid_spacing, cmap2->grid_spacing);
    if (cmap1->cmapdata.size() == cmap2->cmapdata.size() &&
        cmap1->grid_spacing == cmap2->grid_spacing)
    {
        for (size_t g = 0; g < cmap1->cmapdata.size(); g++)
        {
            int i;

            fprintf(fp, "comparing cmap %zu\n", g);

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
            cmp_int(fp, buf1, i, static_cast<int>(id1->functype[i]), static_cast<int>(id2->functype[i]));
            cmp_iparm(fp, buf2, id1->functype[i],
                      id1->iparams[i], id2->iparams[i], ftol, abstol);
        }
        cmp_real(fp, "fudgeQQ", -1, id1->fudgeQQ, id2->fudgeQQ, ftol, abstol);
        cmp_cmap(fp, id1->cmap_grid, id2->cmap_grid, ftol, abstol);
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

static void compareFfparams(FILE *fp, const gmx_ffparams_t &ff1, const gmx_ffparams_t &ff2, real ftol, real abstol)
{
    fprintf(fp, "comparing ffparams\n");
    cmp_int(fp, "numTypes", -1, ff1.numTypes(), ff2.numTypes());
    cmp_int(fp, "atnr", -1, ff1.atnr, ff1.atnr);
    cmp_double(fp, "reppow", -1, ff1.reppow, ff2.reppow, ftol, abstol);
    cmp_real(fp, "fudgeQQ", -1, ff1.fudgeQQ, ff2.fudgeQQ, ftol, abstol);
    cmp_cmap(fp, &ff1.cmap_grid, &ff2.cmap_grid, ftol, abstol);
    int nr = std::min(ff1.numTypes(), ff2.numTypes());
    for (int i = 0; i < nr; i++)
    {
        std::string buf = gmx::formatString("ffparams->functype[%d]", i);
        cmp_int(fp, buf.c_str(), i, ff1.functype[i], ff2.functype[i]);
        buf = gmx::formatString("ffparams->iparams[%d]", i);
        cmp_iparm(fp, buf.c_str(), ff1.functype[i], ff1.iparams[i], ff2.iparams[i], ftol, abstol);
    }

}

static void compareAtom(FILE *fp, int atom, const t_atom &a1, const t_atom &a2, real ftol, real abstol)
{
    fprintf(fp, "comparing t_atom\n");
    cmp_real(fp, "m", atom, a1.m, a2.m, ftol, abstol);
    cmp_real(fp, "q", atom, a1.q, a2.q, ftol, abstol);
    cmp_us(fp, "type", atom, a1.type, a2.type);
    cmp_us(fp, "typeB", atom, a1.typeB, a2.typeB);
    cmp_int(fp, "ptype", atom, a1.ptype, a2.ptype);
    cmp_int(fp, "resind", atom, a1.resind, a2.resind);
    cmp_int(fp, "atomnumber", atom, a1.atomnumber, a2.atomnumber);
    cmp_str(fp, "elem", atom, a1.elem, a2.elem);
}

static void compareResinfo(FILE *fp, int residue, const t_resinfo &r1, const t_resinfo &r2)
{
    fprintf(fp, "comparing t_resinfo\n");
    cmp_str(fp, "name", residue, *r1.name, *r2.name);
    cmp_int(fp, "nr", residue, r1.nr, r2.nr);
    cmp_uc(fp, "ic", residue, r1.ic, r2.ic);
    cmp_int(fp, "chainnum", residue, r1.chainnum, r2.chainnum);
    cmp_uc(fp, "chainid", residue, r1.chainid, r2.chainid);
    if ((r1.rtp || r2.rtp ) && (!r1.rtp || !r2.rtp))
    {
        fprintf(fp, "rtp info is present in topology %d but not in the other\n", r1.rtp ? 1 : 2);
    }
    if (r1.rtp && r2.rtp)
    {
        cmp_str(fp, "rtp", residue, *r1.rtp, *r2.rtp);
    }
}

static void comparePdbinfo(FILE *fp, int pdb, const t_pdbinfo &pdb1, const t_pdbinfo &pdb2, real ftol, real abstol)
{
    fprintf(fp, "comparing t_pdbinfo\n");
    cmp_int(fp, "type", pdb, pdb1.type, pdb2.type);
    cmp_int(fp, "atomnr", pdb, pdb1.atomnr, pdb2.atomnr);
    cmp_uc(fp, "altloc", pdb, pdb1.altloc, pdb2.altloc);
    cmp_str(fp, "atomnm", pdb, pdb1.atomnm, pdb2.atomnm);
    cmp_real(fp, "occup", pdb, pdb1.occup, pdb2.occup, ftol, abstol);
    cmp_real(fp, "bfac", pdb, pdb1.bfac, pdb2.bfac, ftol, abstol);
    cmp_bool(fp, "bAnistropic", pdb, pdb1.bAnisotropic, pdb2.bAnisotropic);
    for (int i = 0; i < 6; i++)
    {
        std::string buf = gmx::formatString("uij[%d]", i);
        cmp_int(fp, buf.c_str(), pdb, pdb1.uij[i], pdb2.uij[i]);
    }
}

static void compareAtoms(FILE *fp, const t_atoms &a1, const t_atoms &a2, real ftol, real abstol)
{
    fprintf(fp, "comparing t_atoms\n");
    cmp_int(fp, "nr", -1, a1.nr, a2.nr);
    cmp_int(fp, "nres", -1, a1.nres, a2.nres);
    cmp_bool(fp, "haveMass", -1, a1.haveMass, a2.haveMass);
    cmp_bool(fp, "haveCharge", -1, a1.haveCharge, a2.haveCharge);
    cmp_bool(fp, "haveType", -1, a1.haveType, a2.haveType);
    cmp_bool(fp, "haveBState", -1, a1.haveBState, a2.haveBState);
    cmp_bool(fp, "havePdbInfo", -1, a1.havePdbInfo, a2.havePdbInfo);

    if ((a1.atom || a2.atom) && (!a1.atom || !a2.atom))
    {
        fprintf(fp, "t_atoms are present in topology %d but not in the other\n", a1.atom ? 1 : 2);
    }
    int nr = std::min(a1.nr, a2.nr);
    for (int i = 0; i < nr; i++)
    {
        if (a1.atom && a2.atom)
        {
            compareAtom(fp, i, a1.atom[i], a2.atom[i], ftol, abstol);
        }
        if (a1.atomname && a2.atomname)
        {
            cmp_str(fp, "atomname", i, *a1.atomname[i], *a2.atomname[i]);
        }
        if (a1.havePdbInfo && a2.havePdbInfo)
        {
            comparePdbinfo(fp, i, a1.pdbinfo[i], a2.pdbinfo[i], ftol, abstol);
        }
        if (a1.haveType && a2.haveType)
        {
            cmp_str(fp, "atomtype", i, *a1.atomtype[i], *a2.atomtype[i]);
        }
        if (a1.haveBState && a2.haveBState)
        {
            cmp_str(fp, "atomtypeB", i, *a1.atomtypeB[i], *a2.atomtypeB[i]);
        }
    }
    nr = std::min(a1.nres, a2.nres);
    for (int i = 0; i < nr; i++)
    {
        compareResinfo(fp, i, a1.resinfo[i], a2.resinfo[i]);
    }
}

static void compareIntermolecular(FILE *fp, const InteractionLists *il1, const InteractionLists *il2)
{
    fprintf(fp, "comparing InteractionLists\n");
    if ((il1 || il2) && (!il1 || !il2))
    {
        fprintf(fp, "InteractionLists is present in topology %d but not in the other\n", il1 ? 1 : 2);
    }
    if (il1 && il2)
    {
        for (int i = 0; i < F_NRE; i++)
        {
            cmp_int(fp, "InteractionList size", i, il1->at(i).size(), il2->at(i).size());
            int nr = std::min(il1->at(i).size(), il2->at(i).size());
            for (int j = 0; j < nr; j++)
            {
                cmp_int(fp, "InteractionList entry", j, il1->at(i).iatoms.at(j), il2->at(i).iatoms.at(j));
            }
        }
    }
}

static void compareMoltypes(FILE *fp, gmx::ArrayRef<const gmx_moltype_t> mt1, gmx::ArrayRef<const gmx_moltype_t> mt2, real ftol, real abstol)
{
    fprintf(fp, "comparing moltypes\n");
    cmp_int(fp, "moltype size", -1, mt1.size(), mt2.size());
    int nr = std::min(mt1.size(), mt2.size());
    for (int i = 0; i < nr; i++)
    {
        cmp_str(fp, "Name", i, *mt1[i].name, *mt2[i].name);
        compareAtoms(fp, mt1[i].atoms, mt2[i].atoms, ftol, abstol);
        compareIntermolecular(fp, &mt1[i].ilist, &mt2[i].ilist);
        std::string buf = gmx::formatString("cgs[%d]", i);
        cmp_block(fp, &mt1[i].cgs, &mt2[i].cgs, buf.c_str());
        buf = gmx::formatString("excls[%d]", i);
        cmp_blocka(fp, &mt1[i].excls, &mt2[i].excls, buf.c_str());
    }
}

static void compareMolblocks(FILE *fp, gmx::ArrayRef<const gmx_molblock_t> mb1, gmx::ArrayRef<const gmx_molblock_t> mb2)
{
    fprintf(fp, "comparing molblocks\n");
    cmp_int(fp, "molblock size", -1, mb1.size(), mb2.size());
    int nr = std::min(mb1.size(), mb2.size());
    for (int i = 0; i < nr; i++)
    {
        cmp_int(fp, "type", i, mb1[i].type, mb2[i].type);
        cmp_int(fp, "nmol", i, mb1[i].nmol, mb2[i].nmol);
        // Only checking size of restraint vectors for now
        cmp_int(fp, "posres_xA size", i, mb1[i].posres_xA.size(), mb2[i].posres_xA.size());
        cmp_int(fp, "posres_xB size", i, mb1[i].posres_xB.size(), mb2[i].posres_xB.size());
    }

}

static void compareAtomtypes(FILE *fp, const t_atomtypes &at1, const t_atomtypes &at2)
{
    fprintf(fp, "comparing atomtypes\n");
    cmp_int(fp, "nr", -1, at1.nr, at2.nr);
    int nr = std::min(at1.nr, at2.nr);
    if (nr > 0)
    {
        for (int i = 0; i < nr; i++)
        {
            cmp_int(fp, "atomtype", i, at1.atomnumber[i], at2.atomnumber[i]);
        }
    }

}

static void compareIntermolecularExclusions(FILE *fp, gmx::ArrayRef<const int> ime1, gmx::ArrayRef<const int> ime2)
{
    fprintf(fp, "comparing intermolecular exclusions\n");
    cmp_int(fp, "exclusion number", -1, ime1.size(), ime2.size());
    int nr = std::min(ime1.size(), ime2.size());
    if (nr > 0)
    {
        for (int i = 0; i < nr; i++)
        {
            cmp_int(fp, "exclusion", i, ime1[i], ime2[i]);
        }
    }

}

static void compareBlockIndices(FILE *fp, gmx::ArrayRef<const MoleculeBlockIndices> mbi1, gmx::ArrayRef<const MoleculeBlockIndices> mbi2)
{
    fprintf(fp, "comparing moleculeBlockIndices\n");
    cmp_int(fp, "size", -1, mbi1.size(), mbi2.size());
    int nr = std::min(mbi1.size(), mbi2.size());
    for (int i = 0; i < nr; i++)
    {
        cmp_int(fp, "numAtomsPerMolecule", i, mbi1[i].numAtomsPerMolecule, mbi2[i].numAtomsPerMolecule);
        cmp_int(fp, "globalAtomStart", i, mbi1[i].globalAtomStart, mbi2[i].globalAtomStart);
        cmp_int(fp, "globalAtomEnd", i, mbi1[i].globalAtomEnd, mbi2[i].globalAtomEnd);
        cmp_int(fp, "globalResidueStart", i, mbi1[i].globalResidueStart, mbi2[i].globalResidueStart);
        cmp_int(fp, "moleculeIndexStart", i, mbi1[i].moleculeIndexStart, mbi2[i].moleculeIndexStart);
    }
}

void compareMtop(FILE *fp, const gmx_mtop_t *mtop1, const gmx_mtop_t *mtop2, real ftol, real abstol)
{
    fprintf(fp, "comparing top\n");
    if (mtop2)
    {
        cmp_str(fp, "Name", -1, *mtop1->name, *mtop2->name);
        cmp_int(fp, "natoms", -1, mtop1->natoms, mtop2->natoms);
        cmp_int(fp, "maxres_renum", -1, mtop1->maxres_renum, mtop2->maxres_renum);
        cmp_int(fp, "maxresnr", -1, mtop1->maxresnr, mtop2->maxresnr);
        cmp_bool(fp, "bIntermolecularInteractions", -1, mtop1->bIntermolecularInteractions, mtop2->bIntermolecularInteractions);
        cmp_bool(fp, "haveMoleculeIndices", -1, mtop1->haveMoleculeIndices, mtop2->haveMoleculeIndices);

        compareFfparams(fp, mtop1->ffparams, mtop2->ffparams, ftol, abstol);
        compareMoltypes(fp, mtop1->moltype, mtop2->moltype, ftol, abstol);
        compareMolblocks(fp, mtop1->molblock, mtop2->molblock);
        compareIntermolecular(fp, mtop1->intermolecular_ilist.get(), mtop2->intermolecular_ilist.get());
        compareAtomtypes(fp, mtop1->atomtypes, mtop2->atomtypes);
        compareGmxGroups(fp, mtop1->groups, mtop2->groups, mtop1->natoms, mtop2->natoms);
        compareIntermolecularExclusions(fp, mtop1->intermolecularExclusionGroup, mtop2->intermolecularExclusionGroup);
        compareBlockIndices(fp, mtop1->moleculeBlockIndices, mtop2->moleculeBlockIndices);

    }

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

void compareGmxGroups(FILE *fp, const gmx_groups_t &g0, const gmx_groups_t &g1,
                      int natoms0, int natoms1)
{
    fprintf(fp, "comparing groups\n");

    for (int i = 0; i < egcNR; i++)
    {
        std::string buf = gmx::formatString("grps[%d].nr", i);
        cmp_int(fp, buf.c_str(), -1, g0.grps[i].nr, g1.grps[i].nr);
        if (g0.grps[i].nr == g1.grps[i].nr)
        {
            for (int j = 0; j < g0.grps[i].nr; j++)
            {
                buf = gmx::formatString("grps[%d].name[%d]", i, j);
                cmp_str(fp, buf.c_str(), -1,
                        *g0.grpname[g0.grps[i].nm_ind[j]],
                        *g1.grpname[g1.grps[i].nm_ind[j]]);
            }
        }
        cmp_int(fp, "ngrpnr", i, g0.ngrpnr[i], g1.ngrpnr[i]);
        if (g0.ngrpnr[i] == g1.ngrpnr[i] && natoms0 == natoms1 &&
            (g0.grpnr[i] != nullptr || g1.grpnr[i] != nullptr))
        {
            for (int j = 0; j < natoms0; j++)
            {
                cmp_int(fp, gtypes[i], j, getGroupType(g0, i, j), getGroupType(g1, i, j));
            }
        }
    }
    /* We have compared the names in the groups lists,
     * so we can skip the grpname list comparison.
     */
}

int getGroupType(const gmx_groups_t &group, int type, int atom)
{
    return (group.grpnr[type] ? group.grpnr[type][atom] : 0);
}

void copy_moltype(const gmx_moltype_t *src, gmx_moltype_t *dst)
{
    dst->name = src->name;
    copy_blocka(&src->excls, &dst->excls);
    copy_block(&src->cgs, &dst->cgs);
    t_atoms *atomsCopy = copy_t_atoms(&src->atoms);
    dst->atoms = *atomsCopy;
    sfree(atomsCopy);

    for (int i = 0; i < F_NRE; ++i)
    {
        dst->ilist[i] = src->ilist[i];
    }
}
