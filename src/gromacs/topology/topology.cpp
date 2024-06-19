/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
#include "gmxpre.h"

#include "gromacs/topology/topology.h"

#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/forcefieldparameters.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/compare.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/txtdump.h"

const char* shortName(SimulationAtomGroupType type)
{
    constexpr gmx::EnumerationArray<SimulationAtomGroupType, const char*> sc_simulationAtomGroupTypeShortNames = {
        { "T-Coupling",
          "Energy Mon.",
          "Acc. not used",
          "Freeze",
          "User1",
          "User2",
          "VCM",
          "Compressed X",
          "Or. Res. Fit",
          "QMMM" }
    };

    return sc_simulationAtomGroupTypeShortNames[type];
}

gmx_moltype_t::gmx_moltype_t() : name(nullptr)
{
    init_t_atoms(&atoms, 0, FALSE);
}

gmx_moltype_t::~gmx_moltype_t()
{
    done_atom(&atoms);
}

static int gmx_mtop_maxresnr(const gmx::ArrayRef<const gmx_moltype_t> moltypes, int maxres_renum)
{
    int maxresnr = 0;

    for (const gmx_moltype_t& moltype : moltypes)
    {
        const t_atoms& atoms = moltype.atoms;
        if (atoms.nres > maxres_renum)
        {
            for (int r = 0; r < atoms.nres; r++)
            {
                if (atoms.resinfo[r].nr > maxresnr)
                {
                    maxresnr = atoms.resinfo[r].nr;
                }
            }
        }
    }

    return maxresnr;
}

gmx_mtop_t::gmx_mtop_t()
{
    open_symtab(&symtab);
}

gmx_mtop_t::~gmx_mtop_t()
{
    done_symtab(&symtab);

    moltype.clear();
    molblock.clear();
}

void gmx_mtop_t::finalize()
{
    if (molblock.size() == 1 && molblock[0].nmol == 1)
    {
        /* We have a single molecule only, no renumbering needed.
         * This case also covers an mtop converted from pdb/gro/... input,
         * so we retain the original residue numbering.
         */
        maxResiduesPerMoleculeToTriggerRenumber_ = 0;
    }
    else
    {
        /* We only renumber single residue molecules. Their intra-molecular
         * residue numbering is anyhow irrelevant.
         */
        maxResiduesPerMoleculeToTriggerRenumber_ = 1;
    }

    const char* env = getenv("GMX_MAXRESRENUM");
    if (env != nullptr)
    {
        sscanf(env, "%d", &maxResiduesPerMoleculeToTriggerRenumber_);
    }
    if (maxResiduesPerMoleculeToTriggerRenumber_ == -1)
    {
        /* -1 signals renumber residues in all molecules */
        maxResiduesPerMoleculeToTriggerRenumber_ = std::numeric_limits<int>::max();
    }

    maxResNumberNotRenumbered_ = gmx_mtop_maxresnr(moltype, maxResiduesPerMoleculeToTriggerRenumber_);

    buildMolblockIndices();
}

void gmx_mtop_t::buildMolblockIndices()
{
    moleculeBlockIndices.resize(molblock.size());

    int atomIndex          = 0;
    int residueIndex       = 0;
    int residueNumberStart = maxResNumberNotRenumbered_ + 1;
    int moleculeIndexStart = 0;
    for (size_t mb = 0; mb < molblock.size(); mb++)
    {
        const gmx_molblock_t& molb         = molblock[mb];
        MoleculeBlockIndices& indices      = moleculeBlockIndices[mb];
        const int             numResPerMol = moltype[molb.type].atoms.nres;

        indices.numAtomsPerMolecule = moltype[molb.type].atoms.nr;
        indices.globalAtomStart     = atomIndex;
        indices.globalResidueStart  = residueIndex;
        atomIndex += molb.nmol * indices.numAtomsPerMolecule;
        residueIndex += molb.nmol * numResPerMol;
        indices.globalAtomEnd      = atomIndex;
        indices.residueNumberStart = residueNumberStart;
        if (numResPerMol <= maxResiduesPerMoleculeToTriggerRenumber_)
        {
            residueNumberStart += molb.nmol * numResPerMol;
        }
        indices.moleculeIndexStart = moleculeIndexStart;
        moleculeIndexStart += molb.nmol;
    }
}

void done_top(t_topology* top)
{
    done_idef(&top->idef);
    done_atom(&(top->atoms));

    done_symtab(&(top->symtab));
    done_block(&(top->mols));
}

void done_top_mtop(t_topology* top, gmx_mtop_t* mtop)
{
    if (mtop != nullptr)
    {
        if (top != nullptr)
        {
            done_idef(&top->idef);
            done_atom(&top->atoms);
            done_block(&top->mols);
            done_symtab(&top->symtab);
            open_symtab(&mtop->symtab);
        }

        // Note that the rest of mtop will be freed by the destructor
    }
}

gmx_localtop_t::gmx_localtop_t(const gmx_ffparams_t& ffparams) : idef(ffparams) {}

bool gmx_mtop_has_masses(const gmx_mtop_t* mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->moltype.empty() || mtop->moltype[0].atoms.haveMass;
}

bool gmx_mtop_has_charges(const gmx_mtop_t* mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->moltype.empty() || mtop->moltype[0].atoms.haveCharge;
}

bool gmx_mtop_has_atomtypes(const gmx_mtop_t* mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->moltype.empty() || mtop->moltype[0].atoms.haveType;
}

bool gmx_mtop_has_pdbinfo(const gmx_mtop_t* mtop)
{
    if (mtop == nullptr)
    {
        return false;
    }
    return mtop->moltype.empty() || mtop->moltype[0].atoms.havePdbInfo;
}

static void pr_grps(FILE*                                 fp,
                    const char*                           title,
                    gmx::ArrayRef<const AtomGroupIndices> grps,
                    const char* const* const*             grpname)
{
    int index = 0;
    for (const auto& group : grps)
    {
        fprintf(fp,
                "%s[%-12s] nr=%zu, name=[",
                title,
                shortName(static_cast<SimulationAtomGroupType>(index)),
                group.size());
        for (const auto& entry : group)
        {
            fprintf(fp, " %s", *(grpname[entry]));
        }
        fprintf(fp, "]\n");
        index++;
    }
}

static void pr_groups(FILE* fp, int indent, const SimulationGroups& groups, gmx_bool bShowNumbers)
{
    pr_grps(fp, "grp", groups.groups, groups.groupNames.data());
    pr_strings(fp, indent, "grpname", groups.groupNames.data(), groups.groupNames.size(), bShowNumbers);

    pr_indent(fp, indent);
    fprintf(fp, "groups          ");
    for (const auto group : gmx::EnumerationWrapper<SimulationAtomGroupType>{})
    {
        printf(" %5.5s", shortName(group));
    }
    printf("\n");

    pr_indent(fp, indent);
    fprintf(fp, "allocated       ");
    int nat_max = 0;
    for (auto group : keysOf(groups.groups))
    {
        printf(" %5d", groups.numberOfGroupNumbers(group));
        nat_max = std::max(nat_max, groups.numberOfGroupNumbers(group));
    }
    printf("\n");

    if (nat_max == 0)
    {
        pr_indent(fp, indent);
        fprintf(fp, "groupnr[%5s] =", "*");
        for (auto gmx_unused group : keysOf(groups.groups))
        {
            fprintf(fp, "  %3d ", 0);
        }
        fprintf(fp, "\n");
    }
    else
    {
        for (int i = 0; i < nat_max; i++)
        {
            pr_indent(fp, indent);
            fprintf(fp, "groupnr[%5d] =", i);
            for (auto group : keysOf(groups.groups))
            {
                fprintf(fp,
                        "  %3d ",
                        !groups.groupNumbers[group].empty() ? groups.groupNumbers[group][i] : 0);
            }
            fprintf(fp, "\n");
        }
    }
}

static void pr_moltype(FILE*                 fp,
                       int                   indent,
                       const char*           title,
                       const gmx_moltype_t*  molt,
                       int                   n,
                       const gmx_ffparams_t* ffparams,
                       gmx_bool              bShowNumbers,
                       gmx_bool              bShowParameters)
{
    indent = pr_title_n(fp, indent, title, n);
    pr_indent(fp, indent);
    fprintf(fp, "name=\"%s\"\n", *(molt->name));
    pr_atoms(fp, indent, "atoms", &(molt->atoms), bShowNumbers);
    pr_listoflists(fp, indent, "excls", &molt->excls, bShowNumbers);
    for (int j = 0; (j < F_NRE); j++)
    {
        pr_ilist(fp,
                 indent,
                 interaction_function[j].longname,
                 ffparams->functype.data(),
                 molt->ilist[j],
                 bShowNumbers,
                 bShowParameters,
                 ffparams->iparams.data());
    }
}

static void pr_molblock(FILE*                             fp,
                        int                               indent,
                        const char*                       title,
                        const gmx_molblock_t*             molb,
                        int                               n,
                        const std::vector<gmx_moltype_t>& molt)
{
    indent = pr_title_n(fp, indent, title, n);
    pr_indent(fp, indent);
    fprintf(fp, "%-20s = %d \"%s\"\n", "moltype", molb->type, *(molt[molb->type].name));
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

void pr_mtop(FILE* fp, int indent, const char* title, const gmx_mtop_t* mtop, gmx_bool bShowNumbers, gmx_bool bShowParameters)
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
        pr_str(fp, indent, "bIntermolecularInteractions", gmx::boolToString(mtop->bIntermolecularInteractions));
        if (mtop->bIntermolecularInteractions)
        {
            for (int j = 0; j < F_NRE; j++)
            {
                pr_ilist(fp,
                         indent,
                         interaction_function[j].longname,
                         mtop->ffparams.functype.data(),
                         (*mtop->intermolecular_ilist)[j],
                         bShowNumbers,
                         bShowParameters,
                         mtop->ffparams.iparams.data());
            }
        }
        pr_ffparams(fp, indent, "ffparams", &(mtop->ffparams), bShowNumbers);
        for (size_t mt = 0; mt < mtop->moltype.size(); mt++)
        {
            pr_moltype(fp, indent, "moltype", &mtop->moltype[mt], mt, &mtop->ffparams, bShowNumbers, bShowParameters);
        }
        pr_groups(fp, indent, mtop->groups, bShowNumbers);
    }
}

void pr_top(FILE* fp, int indent, const char* title, const t_topology* top, gmx_bool bShowNumbers, gmx_bool bShowParameters)
{
    if (available(fp, top, indent, title))
    {
        indent = pr_title(fp, indent, title);
        pr_indent(fp, indent);
        fprintf(fp, "name=\"%s\"\n", *(top->name));
        pr_atoms(fp, indent, "atoms", &(top->atoms), bShowNumbers);
        pr_block(fp, indent, "mols", &top->mols, bShowNumbers);
        pr_str(fp, indent, "bIntermolecularInteractions", gmx::boolToString(top->bIntermolecularInteractions));
        pr_idef(fp, indent, "idef", &top->idef, bShowNumbers, bShowParameters);
    }
}

static void cmp_iparm(FILE*            fp,
                      const char*      s,
                      t_functype       ft,
                      const t_iparams& ip1,
                      const t_iparams& ip2,
                      real             relativeTolerance,
                      real             absoluteTolerance)
{
    bool bDiff = false;
    for (int i = 0; i < MAXFORCEPARAM && !bDiff; i++)
    {
        bDiff = !equal_real(ip1.generic.buf[i], ip2.generic.buf[i], relativeTolerance, absoluteTolerance);
    }
    if (bDiff)
    {
        fprintf(fp, "%s1: ", s);
        pr_iparams(fp, ft, ip1);
        fprintf(fp, "%s2: ", s);
        pr_iparams(fp, ft, ip2);
    }
}

static void cmp_iparm_AB(FILE* fp, const char* s, t_functype ft, const t_iparams& ip1, real relativeTolerance, real absoluteTolerance)
{
    /* Normally the first parameter is perturbable */
    int p0    = 0;
    int nrfpA = interaction_function[ft].nrfpA;
    int nrfpB = interaction_function[ft].nrfpB;
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
    bool bDiff = false;
    for (int i = 0; i < nrfpB && !bDiff; i++)
    {
        bDiff = !equal_real(
                ip1.generic.buf[p0 + i], ip1.generic.buf[nrfpA + i], relativeTolerance, absoluteTolerance);
    }
    if (bDiff)
    {
        fprintf(fp, "%s: ", s);
        pr_iparams(fp, ft, ip1);
    }
}

static void cmp_cmap(FILE* fp, const gmx_cmap_t* cmap1, const gmx_cmap_t* cmap2, real relativeTolerance, real absoluteTolerance)
{
    int cmap1_ngrid = (cmap1 ? cmap1->cmapdata.size() : 0);
    int cmap2_ngrid = (cmap2 ? cmap2->cmapdata.size() : 0);

    cmp_int(fp, "cmap ngrid", -1, cmap1_ngrid, cmap2_ngrid);

    if (cmap1 == nullptr || cmap2 == nullptr)
    {
        return;
    }

    cmp_int(fp, "cmap grid_spacing", -1, cmap1->grid_spacing, cmap2->grid_spacing);
    if (cmap1->cmapdata.size() == cmap2->cmapdata.size() && cmap1->grid_spacing == cmap2->grid_spacing)
    {
        for (size_t g = 0; g < cmap1->cmapdata.size(); g++)
        {
            fprintf(fp, "comparing cmap %zu\n", g);

            for (int i = 0; i < 4 * cmap1->grid_spacing * cmap1->grid_spacing; i++)
            {
                cmp_real(fp, "", i, cmap1->cmapdata[g].cmap[i], cmap2->cmapdata[g].cmap[i], relativeTolerance, absoluteTolerance);
            }
        }
    }
}

static void cmp_listoflists(FILE*                        fp,
                            const gmx::ListOfLists<int>& list1,
                            const gmx::ListOfLists<int>& list2,
                            const char*                  s)
{
    char buf[32];

    fprintf(fp, "comparing blocka %s\n", s);
    sprintf(buf, "%s.numLists", s);
    cmp_int(fp, buf, -1, list1.ssize(), list2.ssize());
    sprintf(buf, "%s.numElements", s);
    cmp_int(fp, buf, -1, list1.numElements(), list2.numElements());
}

static void compareFfparams(FILE*                 fp,
                            const gmx_ffparams_t& ff1,
                            const gmx_ffparams_t& ff2,
                            real                  relativeTolerance,
                            real                  absoluteTolerance)
{
    fprintf(fp, "comparing force field parameters\n");
    cmp_int(fp, "numTypes", -1, ff1.numTypes(), ff2.numTypes());
    cmp_int(fp, "atnr", -1, ff1.atnr, ff1.atnr);
    cmp_double(fp, "reppow", -1, ff1.reppow, ff2.reppow, relativeTolerance, absoluteTolerance);
    cmp_real(fp, "fudgeQQ", -1, ff1.fudgeQQ, ff2.fudgeQQ, relativeTolerance, absoluteTolerance);
    cmp_cmap(fp, &ff1.cmap_grid, &ff2.cmap_grid, relativeTolerance, absoluteTolerance);
    for (int i = 0; i < std::min(ff1.numTypes(), ff2.numTypes()); i++)
    {
        std::string buf = gmx::formatString("ffparams->functype[%d]", i);
        cmp_int(fp, buf.c_str(), i, ff1.functype[i], ff2.functype[i]);
        buf = gmx::formatString("ffparams->iparams[%d]", i);
        cmp_iparm(fp, buf.c_str(), ff1.functype[i], ff1.iparams[i], ff2.iparams[i], relativeTolerance, absoluteTolerance);
    }
}

static void compareFfparamAB(FILE* fp, const gmx_ffparams_t& ff1, real relativeTolerance, real absoluteTolerance)
{
    fprintf(fp, "comparing free energy parameters\n");
    for (int i = 0; i < ff1.numTypes(); i++)
    {
        std::string buf = gmx::formatString("ffparams->iparams[%d]", i);
        cmp_iparm_AB(fp, buf.c_str(), ff1.functype[i], ff1.iparams[i], relativeTolerance, absoluteTolerance);
    }
}
static void compareInteractionLists(FILE* fp, const InteractionLists* il1, const InteractionLists* il2)
{
    fprintf(fp, "comparing InteractionLists\n");
    if ((il1 || il2) && (!il1 || !il2))
    {
        fprintf(fp, "InteractionLists are present in topology %d but not in the other\n", il1 ? 1 : 2);
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

static void compareMoltypes(FILE*                              fp,
                            gmx::ArrayRef<const gmx_moltype_t> mt1,
                            gmx::ArrayRef<const gmx_moltype_t> mt2,
                            real                               relativeTolerance,
                            real                               absoluteTolerance)
{
    fprintf(fp, "comparing molecule types\n");
    cmp_int(fp, "moltype size", -1, mt1.size(), mt2.size());
    for (int i = 0; i < std::min(mt1.ssize(), mt2.ssize()); i++)
    {
        cmp_str(fp, "Name", i, *mt1[i].name, *mt2[i].name);
        compareAtoms(fp, &mt1[i].atoms, &mt2[i].atoms, relativeTolerance, absoluteTolerance);
        compareInteractionLists(fp, &mt1[i].ilist, &mt2[i].ilist);
        std::string buf = gmx::formatString("excls[%d]", i);
        cmp_listoflists(fp, mt1[i].excls, mt2[i].excls, buf.c_str());
    }
}

static void compareMoletypeAB(FILE* fp, gmx::ArrayRef<const gmx_moltype_t> mt1, real relativeTolerance, real absoluteTolerance)
{
    fprintf(fp, "comparing free energy molecule types\n");
    for (gmx::Index i = 0; i < mt1.ssize(); i++)
    {
        compareAtoms(fp, &mt1[i].atoms, nullptr, relativeTolerance, absoluteTolerance);
    }
}
static void compareMolblocks(FILE*                               fp,
                             gmx::ArrayRef<const gmx_molblock_t> mb1,
                             gmx::ArrayRef<const gmx_molblock_t> mb2)
{
    fprintf(fp, "comparing molecule blocks\n");
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

static void compareIntermolecularExclusions(FILE*                    fp,
                                            gmx::ArrayRef<const int> ime1,
                                            gmx::ArrayRef<const int> ime2)
{
    fprintf(fp, "comparing intermolecular exclusions\n");
    cmp_int(fp, "exclusion number", -1, ime1.size(), ime2.size());
    int nr = std::min(ime1.size(), ime2.size());
    for (int i = 0; i < nr; i++)
    {
        cmp_int(fp, "exclusion", i, ime1[i], ime2[i]);
    }
}

static void compareBlockIndices(FILE*                                     fp,
                                gmx::ArrayRef<const MoleculeBlockIndices> mbi1,
                                gmx::ArrayRef<const MoleculeBlockIndices> mbi2)
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

void compareMtop(FILE* fp, const gmx_mtop_t& mtop1, const gmx_mtop_t& mtop2, real relativeTolerance, real absoluteTolerance)
{
    fprintf(fp, "comparing mtop topology\n");
    cmp_str(fp, "Name", -1, *mtop1.name, *mtop2.name);
    cmp_int(fp, "natoms", -1, mtop1.natoms, mtop2.natoms);
    cmp_int(fp,
            "maxres_renum",
            -1,
            mtop1.maxResiduesPerMoleculeToTriggerRenumber(),
            mtop2.maxResiduesPerMoleculeToTriggerRenumber());
    cmp_int(fp, "maxresnr", -1, mtop1.maxResNumberNotRenumbered(), mtop2.maxResNumberNotRenumbered());
    cmp_bool(fp, "bIntermolecularInteractions", -1, mtop1.bIntermolecularInteractions, mtop2.bIntermolecularInteractions);
    cmp_bool(fp, "haveMoleculeIndices", -1, mtop1.haveMoleculeIndices, mtop2.haveMoleculeIndices);

    compareFfparams(fp, mtop1.ffparams, mtop2.ffparams, relativeTolerance, absoluteTolerance);
    compareMoltypes(fp, mtop1.moltype, mtop2.moltype, relativeTolerance, absoluteTolerance);
    compareMolblocks(fp, mtop1.molblock, mtop2.molblock);
    compareInteractionLists(fp, mtop1.intermolecular_ilist.get(), mtop2.intermolecular_ilist.get());
    compareAtomGroups(fp, mtop1.groups, mtop2.groups, mtop1.natoms, mtop2.natoms);
    compareIntermolecularExclusions(
            fp, mtop1.intermolecularExclusionGroup, mtop2.intermolecularExclusionGroup);
    compareBlockIndices(fp, mtop1.moleculeBlockIndices, mtop2.moleculeBlockIndices);
}

void compareMtopAB(FILE* fp, const gmx_mtop_t& mtop1, real relativeTolerance, real absoluteTolerance)
{
    fprintf(fp, "comparing topAB\n");
    compareFfparamAB(fp, mtop1.ffparams, relativeTolerance, absoluteTolerance);
    compareMoletypeAB(fp, mtop1.moltype, relativeTolerance, absoluteTolerance);
}

void compareAtomGroups(FILE* fp, const SimulationGroups& g0, const SimulationGroups& g1, int natoms0, int natoms1)
{
    fprintf(fp, "comparing groups\n");

    for (auto group : keysOf(g0.groups))
    {
        std::string buf = gmx::formatString("grps[%d].nr", static_cast<int>(group));
        cmp_int(fp, buf.c_str(), -1, g0.groups[group].size(), g1.groups[group].size());
        if (g0.groups[group].size() == g1.groups[group].size())
        {
            for (gmx::Index j = 0; j < gmx::ssize(g0.groups[group]); j++)
            {
                buf = gmx::formatString("grps[%d].name[%zd]", static_cast<int>(group), j);
                cmp_str(fp,
                        buf.c_str(),
                        -1,
                        *g0.groupNames[g0.groups[group][j]],
                        *g1.groupNames[g1.groups[group][j]]);
            }
        }
        cmp_int(fp,
                "ngrpnr",
                static_cast<int>(group),
                g0.numberOfGroupNumbers(group),
                g1.numberOfGroupNumbers(group));
        if (g0.numberOfGroupNumbers(group) == g1.numberOfGroupNumbers(group) && natoms0 == natoms1
            && (!g0.groupNumbers[group].empty() || !g1.groupNumbers[group].empty()))
        {
            for (int j = 0; j < natoms0; j++)
            {
                cmp_int(fp, shortName(group), j, getGroupType(g0, group, j), getGroupType(g1, group, j));
            }
        }
    }
    /* We have compared the names in the groups lists,
     * so we can skip the grpname list comparison.
     */
}

int getGroupType(const SimulationGroups& group, SimulationAtomGroupType type, int atom)
{
    return (group.groupNumbers[type].empty() ? 0 : group.groupNumbers[type][atom]);
}

void copy_moltype(const gmx_moltype_t* src, gmx_moltype_t* dst)
{
    dst->name          = src->name;
    dst->excls         = src->excls;
    t_atoms* atomsCopy = copy_t_atoms(&src->atoms);
    dst->atoms         = *atomsCopy;
    sfree(atomsCopy);

    for (int i = 0; i < F_NRE; ++i)
    {
        dst->ilist[i] = src->ilist[i];
    }
}
