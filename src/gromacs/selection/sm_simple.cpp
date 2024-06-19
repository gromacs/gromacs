/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2009- The GROMACS Authors
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
/*! \internal \file
 * \brief
 * Implements simple keyword selection methods.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include <cctype>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/selection/indexutil.h"
#include "gromacs/selection/position.h"
#include "gromacs/selection/selparam.h"
#include "gromacs/selection/selvalue.h"
#include "gromacs/topology/mtop_lookup.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/real.h"

#include "selmethod.h"
#include "selmethod_impl.h"

/** Evaluates the \p all selection keyword. */
static void evaluate_all(const gmx::SelMethodEvalContext& context,
                         gmx_ana_index_t*                 g,
                         gmx_ana_selvalue_t*              out,
                         void*                            data);
/** Evaluates the \p none selection keyword. */
static void evaluate_none(const gmx::SelMethodEvalContext& context,
                          gmx_ana_index_t*                 g,
                          gmx_ana_selvalue_t*              out,
                          void*                            data);
/** Evaluates the \p atomnr selection keyword. */
static void evaluate_atomnr(const gmx::SelMethodEvalContext& context,
                            gmx_ana_index_t*                 g,
                            gmx_ana_selvalue_t*              out,
                            void*                            data);
/** Evaluates the \p resnr selection keyword. */
static void evaluate_resnr(const gmx::SelMethodEvalContext& context,
                           gmx_ana_index_t*                 g,
                           gmx_ana_selvalue_t*              out,
                           void*                            data);
/** Evaluates the \p resindex selection keyword. */
static void evaluate_resindex(const gmx::SelMethodEvalContext& context,
                              gmx_ana_index_t*                 g,
                              gmx_ana_selvalue_t*              out,
                              void*                            data);
/*! \brief
 * Checks whether molecule information is present in the topology.
 *
 * \param[in] top  Topology structure.
 * \param     npar Not used.
 * \param     param Not used.
 * \param     data Not used.
 * \returns   0 if molecule info is present in the topology, -1 otherwise.
 *
 * If molecule information is not found, also prints an error message.
 */
static void check_molecules(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/** Evaluates the \p molindex selection keyword. */
static void evaluate_molindex(const gmx::SelMethodEvalContext& context,
                              gmx_ana_index_t*                 g,
                              gmx_ana_selvalue_t*              out,
                              void*                            data);
/** Evaluates the \p atomname selection keyword. */
static void evaluate_atomname(const gmx::SelMethodEvalContext& context,
                              gmx_ana_index_t*                 g,
                              gmx_ana_selvalue_t*              out,
                              void*                            data);
/** Evaluates the \p pdbatomname selection keyword. */
static void evaluate_pdbatomname(const gmx::SelMethodEvalContext& context,
                                 gmx_ana_index_t*                 g,
                                 gmx_ana_selvalue_t*              out,
                                 void*                            data);
/*! \brief
 * Checks whether atom types are present in the topology.
 *
 * \param[in] top  Topology structure.
 * \param     npar Not used.
 * \param     param Not used.
 * \param     data Not used.
 */
static void check_atomtype(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/** Evaluates the \p atomtype selection keyword. */
static void evaluate_atomtype(const gmx::SelMethodEvalContext& context,
                              gmx_ana_index_t*                 g,
                              gmx_ana_selvalue_t*              out,
                              void*                            data);
/** Evaluates the \p insertcode selection keyword. */
static void evaluate_insertcode(const gmx::SelMethodEvalContext& context,
                                gmx_ana_index_t*                 g,
                                gmx_ana_selvalue_t*              out,
                                void*                            data);
/** Evaluates the \p chain selection keyword. */
static void evaluate_chain(const gmx::SelMethodEvalContext& context,
                           gmx_ana_index_t*                 g,
                           gmx_ana_selvalue_t*              out,
                           void*                            data);
/** Evaluates the \p mass selection keyword. */
static void evaluate_mass(const gmx::SelMethodEvalContext& context,
                          gmx_ana_index_t*                 g,
                          gmx_ana_selvalue_t*              out,
                          void*                            data);
/*! \brief
 * Checks whether charges are present in the topology.
 *
 * \param[in] top  Topology structure.
 * \param     npar Not used.
 * \param     param Not used.
 * \param     data Not used.
 */
static void check_charge(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/** Evaluates the \p charge selection keyword. */
static void evaluate_charge(const gmx::SelMethodEvalContext& context,
                            gmx_ana_index_t*                 g,
                            gmx_ana_selvalue_t*              out,
                            void*                            data);
/*! \brief
 * Checks whether PDB info is present in the topology.
 *
 * \param[in] top  Topology structure.
 * \param     npar Not used.
 * \param     param Not used.
 * \param     data Not used.
 * \returns   0 if PDB info is present in the topology, -1 otherwise.
 *
 * If PDB info is not found, also prints an error message.
 */
static void check_pdbinfo(const gmx_mtop_t* top, int npar, gmx_ana_selparam_t* param, void* data);
/** Evaluates the \p altloc selection keyword. */
static void evaluate_altloc(const gmx::SelMethodEvalContext& context,
                            gmx_ana_index_t*                 g,
                            gmx_ana_selvalue_t*              out,
                            void*                            data);
/** Evaluates the \p occupancy selection keyword. */
static void evaluate_occupancy(const gmx::SelMethodEvalContext& context,
                               gmx_ana_index_t*                 g,
                               gmx_ana_selvalue_t*              out,
                               void*                            data);
/** Evaluates the \p betafactor selection keyword. */
static void evaluate_betafactor(const gmx::SelMethodEvalContext& context,
                                gmx_ana_index_t*                 g,
                                gmx_ana_selvalue_t*              out,
                                void*                            data);
/** Evaluates the \p resname selection keyword. */
static void evaluate_resname(const gmx::SelMethodEvalContext& context,
                             gmx_ana_index_t*                 g,
                             gmx_ana_selvalue_t*              out,
                             void*                            data);

/** Evaluates the \p x selection keyword. */
static void evaluate_x(const gmx::SelMethodEvalContext& context,
                       gmx_ana_pos_t*                   pos,
                       gmx_ana_selvalue_t*              out,
                       void*                            data);
/** Evaluates the \p y selection keyword. */
static void evaluate_y(const gmx::SelMethodEvalContext& context,
                       gmx_ana_pos_t*                   pos,
                       gmx_ana_selvalue_t*              out,
                       void*                            data);
/** Evaluates the \p z selection keyword. */
static void evaluate_z(const gmx::SelMethodEvalContext& context,
                       gmx_ana_pos_t*                   pos,
                       gmx_ana_selvalue_t*              out,
                       void*                            data);

//! Help title for atom name selection keywords.
static const char helptitle_atomname[] = "Selecting atoms by name";
//! Help text for atom name selection keywords.
static const char* const help_atomname[] = {
    "::",
    "",
    "  name",
    "  pdbname",
    "  atomname",
    "  pdbatomname",
    "",
    "These keywords select atoms by name. [TT]name[tt] selects atoms using",
    "the GROMACS atom naming convention.",
    "For input formats other than PDB, the atom names are matched exactly",
    "as they appear in the input file. For PDB files, 4 character atom names",
    "that start with a digit are matched after moving the digit to the end",
    "(e.g., to match 3HG2 from a PDB file, use [TT]name HG23[tt]).",
    "[TT]pdbname[tt] can only be used with a PDB input file, and selects",
    "atoms based on the exact name given in the input file, without the",
    "transformation described above.[PAR]",

    "[TT]atomname[tt] and [TT]pdbatomname[tt] are synonyms for the above two",
    "keywords."
};

//! Help title for residue index selection keywords.
static const char helptitle_resindex[] = "Selecting atoms by residue number";
//! Help text for residue index selection keywords.
static const char* const help_resindex[] = {
    "::",
    "",
    "  resnr",
    "  resid",
    "  resindex",
    "  residue",
    "",
    "[TT]resnr[tt] selects atoms using the residue numbering in the input",
    "file. [TT]resid[tt] is synonym for this keyword for VMD compatibility.",
    "",
    "[TT]resindex N[tt] selects the [TT]N[tt] th residue starting from the",
    "beginning of the input file. This is useful for uniquely identifying",
    "residues if there are duplicate numbers in the input file (e.g., in",
    "multiple chains).",
    "[TT]residue[tt] is a synonym for [TT]resindex[tt]. This allows",
    "[TT]same residue as[tt] to work as expected."
};

/** Selection method data for \p all selection keyword. */
gmx_ana_selmethod_t sm_all = {
    "all",   GROUP_VALUE, 0,       0,       nullptr,       nullptr, nullptr,
    nullptr, nullptr,     nullptr, nullptr, &evaluate_all, nullptr,
};

/** Selection method data for \p none selection keyword. */
gmx_ana_selmethod_t sm_none = {
    "none",  GROUP_VALUE, 0,       0,       nullptr,        nullptr, nullptr,
    nullptr, nullptr,     nullptr, nullptr, &evaluate_none, nullptr,
};

/** Selection method data for \p atomnr selection keyword. */
gmx_ana_selmethod_t sm_atomnr = {
    "atomnr", INT_VALUE, 0,       0,       nullptr,          nullptr, nullptr,
    nullptr,  nullptr,   nullptr, nullptr, &evaluate_atomnr, nullptr,
};

/** Selection method data for \p resnr selection keyword. */
gmx_ana_selmethod_t sm_resnr = {
    "resnr",      INT_VALUE,
    SMETH_REQTOP, 0,
    nullptr,      nullptr,
    nullptr,      nullptr,
    nullptr,      nullptr,
    nullptr,      &evaluate_resnr,
    nullptr,      { nullptr, helptitle_resindex, asize(help_resindex), help_resindex }
};

/** Selection method data for \p resindex selection keyword. */
gmx_ana_selmethod_t sm_resindex = {
    "resindex",   INT_VALUE,
    SMETH_REQTOP, 0,
    nullptr,      nullptr,
    nullptr,      nullptr,
    nullptr,      nullptr,
    nullptr,      &evaluate_resindex,
    nullptr,      { nullptr, helptitle_resindex, asize(help_resindex), help_resindex }
};

/** Selection method data for \p molindex selection keyword. */
gmx_ana_selmethod_t sm_molindex = {
    "molindex", INT_VALUE,        SMETH_REQTOP, 0,       nullptr, nullptr,
    nullptr,    &check_molecules, nullptr,      nullptr, nullptr, &evaluate_molindex,
    nullptr,
};

/** Selection method data for \p atomname selection keyword. */
gmx_ana_selmethod_t sm_atomname = {
    "atomname",   STR_VALUE,
    SMETH_REQTOP, 0,
    nullptr,      nullptr,
    nullptr,      nullptr,
    nullptr,      nullptr,
    nullptr,      &evaluate_atomname,
    nullptr,      { nullptr, helptitle_atomname, asize(help_atomname), help_atomname }
};

/** Selection method data for \p pdbatomname selection keyword. */
gmx_ana_selmethod_t sm_pdbatomname = {
    "pdbatomname", STR_VALUE,
    SMETH_REQTOP,  0,
    nullptr,       nullptr,
    nullptr,       &check_pdbinfo,
    nullptr,       nullptr,
    nullptr,       &evaluate_pdbatomname,
    nullptr,       { nullptr, helptitle_atomname, asize(help_atomname), help_atomname }
};

/** Selection method data for \p atomtype selection keyword. */
gmx_ana_selmethod_t sm_atomtype = {
    "atomtype", STR_VALUE,       SMETH_REQTOP, 0,       nullptr, nullptr,
    nullptr,    &check_atomtype, nullptr,      nullptr, nullptr, &evaluate_atomtype,
    nullptr,
};

/** Selection method data for \p resname selection keyword. */
gmx_ana_selmethod_t sm_resname = {
    "resname", STR_VALUE, SMETH_REQTOP, 0,       nullptr,           nullptr, nullptr,
    nullptr,   nullptr,   nullptr,      nullptr, &evaluate_resname, nullptr,
};

/** Selection method data for \p chain selection keyword. */
gmx_ana_selmethod_t sm_insertcode = {
    "insertcode",
    STR_VALUE,
    SMETH_REQTOP | SMETH_CHARVAL,
    0,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    nullptr,
    &evaluate_insertcode,
    nullptr,
};

/** Selection method data for \p chain selection keyword. */
gmx_ana_selmethod_t sm_chain = {
    "chain", STR_VALUE, SMETH_REQTOP | SMETH_CHARVAL,
    0,       nullptr,   nullptr,
    nullptr, nullptr,   nullptr,
    nullptr, nullptr,   &evaluate_chain,
    nullptr,
};

/** Selection method data for \p mass selection keyword. */
gmx_ana_selmethod_t sm_mass = {
    "mass",  REAL_VALUE, SMETH_REQMASS, 0,       nullptr,        nullptr, nullptr,
    nullptr, nullptr,    nullptr,       nullptr, &evaluate_mass, nullptr,
};

/** Selection method data for \p charge selection keyword. */
gmx_ana_selmethod_t sm_charge = {
    "charge",      REAL_VALUE, SMETH_REQTOP, 0,       nullptr,          nullptr, nullptr,
    &check_charge, nullptr,    nullptr,      nullptr, &evaluate_charge, nullptr,
};

/** Selection method data for \p chain selection keyword. */
gmx_ana_selmethod_t sm_altloc = {
    "altloc", STR_VALUE,      SMETH_REQTOP | SMETH_CHARVAL,
    0,        nullptr,        nullptr,
    nullptr,  &check_pdbinfo, nullptr,
    nullptr,  nullptr,        &evaluate_altloc,
    nullptr,
};

/** Selection method data for \p occupancy selection keyword. */
gmx_ana_selmethod_t sm_occupancy = {
    "occupancy", REAL_VALUE,     SMETH_REQTOP, 0,       nullptr, nullptr,
    nullptr,     &check_pdbinfo, nullptr,      nullptr, nullptr, &evaluate_occupancy,
    nullptr,
};

/** Selection method data for \p betafactor selection keyword. */
gmx_ana_selmethod_t sm_betafactor = {
    "betafactor", REAL_VALUE,     SMETH_REQTOP, 0,       nullptr, nullptr,
    nullptr,      &check_pdbinfo, nullptr,      nullptr, nullptr, &evaluate_betafactor,
    nullptr,
};

/** Selection method data for \p x selection keyword. */
gmx_ana_selmethod_t sm_x = {
    "x",     REAL_VALUE, SMETH_DYNAMIC, 0,       nullptr, nullptr,     nullptr,
    nullptr, nullptr,    nullptr,       nullptr, nullptr, &evaluate_x,
};

/** Selection method data for \p y selection keyword. */
gmx_ana_selmethod_t sm_y = {
    "y",     REAL_VALUE, SMETH_DYNAMIC, 0,       nullptr, nullptr,     nullptr,
    nullptr, nullptr,    nullptr,       nullptr, nullptr, &evaluate_y,
};

/** Selection method data for \p z selection keyword. */
gmx_ana_selmethod_t sm_z = {
    "z",     REAL_VALUE, SMETH_DYNAMIC, 0,       nullptr, nullptr,     nullptr,
    nullptr, nullptr,    nullptr,       nullptr, nullptr, &evaluate_z,
};

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Copies \p g to \p out->u.g.
 */
static void evaluate_all(const gmx::SelMethodEvalContext& /*context*/,
                         gmx_ana_index_t*    g,
                         gmx_ana_selvalue_t* out,
                         void* /* data */)
{
    gmx_ana_index_copy(out->u.g, g, false);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns an empty \p out->u.g.
 */
static void evaluate_none(const gmx::SelMethodEvalContext& /*context*/,
                          gmx_ana_index_t* /* g */,
                          gmx_ana_selvalue_t* out,
                          void* /* data */)
{
    out->u.g->isize = 0;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the indices for each atom in \p out->u.i.
 */
static void evaluate_atomnr(const gmx::SelMethodEvalContext& /*context*/,
                            gmx_ana_index_t*    g,
                            gmx_ana_selvalue_t* out,
                            void* /* data */)
{
    int i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.i[i] = g->index[i] + 1;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the residue numbers for each atom in \p out->u.i.
 */
static void evaluate_resnr(const gmx::SelMethodEvalContext& context,
                           gmx_ana_index_t*                 g,
                           gmx_ana_selvalue_t*              out,
                           void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        mtopGetAtomAndResidueName(*context.top_, g->index[i], &molb, nullptr, &out->u.i[i], nullptr, nullptr);
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the residue indices for each atom in \p out->u.i.
 */
static void evaluate_resindex(const gmx::SelMethodEvalContext& context,
                              gmx_ana_index_t*                 g,
                              gmx_ana_selvalue_t*              out,
                              void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        int resind;
        mtopGetAtomAndResidueName(*context.top_, g->index[i], &molb, nullptr, nullptr, nullptr, &resind);
        out->u.i[i] = resind + 1;
    }
}

static void check_molecules(const gmx_mtop_t* top, int /* npar */, gmx_ana_selparam_t* /* param */, void* /* data */)
{
    bool bOk;

    bOk = (top != nullptr && top->haveMoleculeIndices);
    if (!bOk)
    {
        GMX_THROW(gmx::InconsistentInputError("Molecule information not available in topology"));
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the molecule indices for each atom in \p out->u.i.
 */
static void evaluate_molindex(const gmx::SelMethodEvalContext& context,
                              gmx_ana_index_t*                 g,
                              gmx_ana_selvalue_t*              out,
                              void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.i[i] = mtopGetMoleculeIndex(*context.top_, g->index[i], &molb) + 1;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the atom name for each atom in \p out->u.s.
 */
static void evaluate_atomname(const gmx::SelMethodEvalContext& context,
                              gmx_ana_index_t*                 g,
                              gmx_ana_selvalue_t*              out,
                              void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        const char* atom_name;
        mtopGetAtomAndResidueName(*context.top_, g->index[i], &molb, &atom_name, nullptr, nullptr, nullptr);
        out->u.s[i] = const_cast<char*>(atom_name);
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the PDB atom name for each atom in \p out->u.s.
 */
static void evaluate_pdbatomname(const gmx::SelMethodEvalContext& context,
                                 gmx_ana_index_t*                 g,
                                 gmx_ana_selvalue_t*              out,
                                 void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        const char* s = mtopGetAtomPdbInfo(*context.top_, g->index[i], &molb).atomnm;
        while (std::isspace(*s))
        {
            ++s;
        }
        out->u.s[i] = const_cast<char*>(s);
    }
}

static void check_atomtype(const gmx_mtop_t* top, int /* npar */, gmx_ana_selparam_t* /* param */, void* /* data */)
{
    if (!gmx_mtop_has_atomtypes(top))
    {
        GMX_THROW(gmx::InconsistentInputError("Atom types not available in topology"));
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the atom type for each atom in \p out->u.s.
 * Segfaults if atom types are not found in the topology.
 */
static void evaluate_atomtype(const gmx::SelMethodEvalContext& context,
                              gmx_ana_index_t*                 g,
                              gmx_ana_selvalue_t*              out,
                              void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        int atomIndexInMolecule;
        mtopGetMolblockIndex(*context.top_, g->index[i], &molb, nullptr, &atomIndexInMolecule);
        const gmx_moltype_t& moltype = context.top_->moltype[context.top_->molblock[molb].type];
        out->u.s[i]                  = *moltype.atoms.atomtype[atomIndexInMolecule];
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the residue name for each atom in \p out->u.s.
 */
static void evaluate_resname(const gmx::SelMethodEvalContext& context,
                             gmx_ana_index_t*                 g,
                             gmx_ana_selvalue_t*              out,
                             void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.s[i] = *mtopGetResidueInfo(*context.top_, g->index[i], &molb).name;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the insertion code for each atom in \p out->u.s.
 */
static void evaluate_insertcode(const gmx::SelMethodEvalContext& context,
                                gmx_ana_index_t*                 g,
                                gmx_ana_selvalue_t*              out,
                                void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.s[i][0] = mtopGetResidueInfo(*context.top_, g->index[i], &molb).ic;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the chain for each atom in \p out->u.s.
 */
static void evaluate_chain(const gmx::SelMethodEvalContext& context,
                           gmx_ana_index_t*                 g,
                           gmx_ana_selvalue_t*              out,
                           void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.s[i][0] = mtopGetResidueInfo(*context.top_, g->index[i], &molb).chainid;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the mass for each atom in \p out->u.r.
 */
static void evaluate_mass(const gmx::SelMethodEvalContext& context,
                          gmx_ana_index_t*                 g,
                          gmx_ana_selvalue_t*              out,
                          void* /* data */)
{
    GMX_RELEASE_ASSERT(gmx_mtop_has_masses(context.top_), "Masses not available for evaluation");
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.r[i] = mtopGetAtomMass(*context.top_, g->index[i], &molb);
    }
}


static void check_charge(const gmx_mtop_t* top, int /* npar */, gmx_ana_selparam_t* /* param */, void* /* data */)
{
    if (!gmx_mtop_has_charges(top))
    {
        GMX_THROW(gmx::InconsistentInputError("Charges not available in topology"));
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the charge for each atom in \p out->u.r.
 */
static void evaluate_charge(const gmx::SelMethodEvalContext& context,
                            gmx_ana_index_t*                 g,
                            gmx_ana_selvalue_t*              out,
                            void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.r[i] = mtopGetAtomParameters(*context.top_, g->index[i], &molb).q;
    }
}

static void check_pdbinfo(const gmx_mtop_t* top, int /* npar */, gmx_ana_selparam_t* /* param */, void* /* data */)
{
    if (!gmx_mtop_has_pdbinfo(top))
    {
        GMX_THROW(gmx::InconsistentInputError("PDB info not available in topology"));
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the alternate location identifier for each atom in \p out->u.s.
 */
static void evaluate_altloc(const gmx::SelMethodEvalContext& context,
                            gmx_ana_index_t*                 g,
                            gmx_ana_selvalue_t*              out,
                            void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.s[i][0] = mtopGetAtomPdbInfo(*context.top_, g->index[i], &molb).altloc;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the occupancy numbers for each atom in \p out->u.r.
 * Segfaults if PDB info is not found in the topology.
 */
static void evaluate_occupancy(const gmx::SelMethodEvalContext& context,
                               gmx_ana_index_t*                 g,
                               gmx_ana_selvalue_t*              out,
                               void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.r[i] = mtopGetAtomPdbInfo(*context.top_, g->index[i], &molb).occup;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the B-factors for each atom in \p out->u.r.
 * Segfaults if PDB info is not found in the topology.
 */
static void evaluate_betafactor(const gmx::SelMethodEvalContext& context,
                                gmx_ana_index_t*                 g,
                                gmx_ana_selvalue_t*              out,
                                void* /* data */)
{
    out->nr  = g->isize;
    int molb = 0;
    for (int i = 0; i < g->isize; ++i)
    {
        out->u.r[i] = mtopGetAtomPdbInfo(*context.top_, g->index[i], &molb).bfac;
    }
}

/*! \brief
 * Internal utility function for position keyword evaluation.
 *
 * \param[out] out  Output array.
 * \param[in]  pos  Position data to use instead of atomic coordinates.
 * \param[in]  d    Coordinate index to evaluate (\p XX, \p YY or \p ZZ).
 *
 * This function is used internally by evaluate_x(), evaluate_y() and
 * evaluate_z() to do the actual evaluation.
 */
static void evaluate_coord(real out[], gmx_ana_pos_t* pos, int d)
{
    for (int i = 0; i < pos->count(); ++i)
    {
        out[i] = pos->x[i][d];
    }
    // TODO: Make this more efficient by directly extracting the coordinates
    // from the frame coordinates for atomic positions instead of going through
    // a position calculation.
}

/*!
 * See sel_updatefunc_pos() for description of the parameters.
 * \p data is not used.
 *
 * Returns the \p x coordinate for each position in \p out->u.r.
 */
static void evaluate_x(const gmx::SelMethodEvalContext& /*context*/,
                       gmx_ana_pos_t*      pos,
                       gmx_ana_selvalue_t* out,
                       void* /*data*/)
{
    out->nr = pos->count();
    evaluate_coord(out->u.r, pos, XX);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the \p y coordinate for each position in \p out->u.r.
 */
static void evaluate_y(const gmx::SelMethodEvalContext& /*context*/,
                       gmx_ana_pos_t*      pos,
                       gmx_ana_selvalue_t* out,
                       void* /*data*/)
{
    out->nr = pos->count();
    evaluate_coord(out->u.r, pos, YY);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the \p z coordinate for each position in \p out->u.r.
 */
static void evaluate_z(const gmx::SelMethodEvalContext& /*context*/,
                       gmx_ana_pos_t*      pos,
                       gmx_ana_selvalue_t* out,
                       void* /*data*/)
{
    out->nr = pos->count();
    evaluate_coord(out->u.r, pos, ZZ);
}
