/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2009,2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements simple keyword selection methods.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_selection
 */
#include "gmxpre.h"

#include <cctype>

#include "gromacs/legacyheaders/macros.h"
#include "gromacs/selection/position.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/exceptions.h"

#include "selmethod.h"

/** Evaluates the \p all selection keyword. */
static void
evaluate_all(t_topology *top, t_trxframe *fr, t_pbc *pbc,
             gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p none selection keyword. */
static void
evaluate_none(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p atomnr selection keyword. */
static void
evaluate_atomnr(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p resnr selection keyword. */
static void
evaluate_resnr(t_topology *top, t_trxframe *fr, t_pbc *pbc,
               gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p resindex selection keyword. */
static void
evaluate_resindex(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
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
static void
check_molecules(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Evaluates the \p molindex selection keyword. */
static void
evaluate_molindex(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p atomname selection keyword. */
static void
evaluate_atomname(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p pdbatomname selection keyword. */
static void
evaluate_pdbatomname(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/*! \brief
 * Checks whether atom types are present in the topology.
 *
 * \param[in] top  Topology structure.
 * \param     npar Not used.
 * \param     param Not used.
 * \param     data Not used.
 * \returns   0 if atom types are present in the topology, -1 otherwise.
 *
 * If the atom types are not found, also prints an error message.
 */
static void
check_atomtype(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Evaluates the \p atomtype selection keyword. */
static void
evaluate_atomtype(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p insertcode selection keyword. */
static void
evaluate_insertcode(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                    gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p chain selection keyword. */
static void
evaluate_chain(t_topology *top, t_trxframe *fr, t_pbc *pbc,
               gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p mass selection keyword. */
static void
evaluate_mass(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p charge selection keyword. */
static void
evaluate_charge(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
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
static void
check_pdbinfo(t_topology *top, int npar, gmx_ana_selparam_t *param, void *data);
/** Evaluates the \p altloc selection keyword. */
static void
evaluate_altloc(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p occupancy selection keyword. */
static void
evaluate_occupancy(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                   gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p betafactor selection keyword. */
static void
evaluate_betafactor(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                    gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p resname selection keyword. */
static void
evaluate_resname(t_topology *top, t_trxframe *fr, t_pbc *pbc,
                 gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void *data);

/** Evaluates the \p x selection keyword. */
static void
evaluate_x(t_topology *top, t_trxframe *fr, t_pbc *pbc,
           gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p y selection keyword. */
static void
evaluate_y(t_topology *top, t_trxframe *fr, t_pbc *pbc,
           gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data);
/** Evaluates the \p z selection keyword. */
static void
evaluate_z(t_topology *top, t_trxframe *fr, t_pbc *pbc,
           gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void *data);

//! Help title for atom name selection keywords.
static const char        helptitle_atomname[] = "Selecting atoms by name";
//! Help text for atom name selection keywords.
static const char *const help_atomname[] = {
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

/** Selection method data for \p all selection keyword. */
gmx_ana_selmethod_t sm_all = {
    "all", GROUP_VALUE, 0,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_all,
    NULL,
};

/** Selection method data for \p none selection keyword. */
gmx_ana_selmethod_t sm_none = {
    "none", GROUP_VALUE, 0,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_none,
    NULL,
};

/** Selection method data for \p atomnr selection keyword. */
gmx_ana_selmethod_t sm_atomnr = {
    "atomnr", INT_VALUE, 0,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_atomnr,
    NULL,
};

/** Selection method data for \p resnr selection keyword. */
gmx_ana_selmethod_t sm_resnr = {
    "resnr", INT_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_resnr,
    NULL,
};

/** Selection method data for \p resindex selection keyword. */
gmx_ana_selmethod_t sm_resindex = {
    "resindex", INT_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_resindex,
    NULL,
};

/** Selection method data for \p molindex selection keyword. */
gmx_ana_selmethod_t sm_molindex = {
    "molindex", INT_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    &check_molecules,
    NULL,
    NULL,
    NULL,
    &evaluate_molindex,
    NULL,
};

/** Selection method data for \p atomname selection keyword. */
gmx_ana_selmethod_t sm_atomname = {
    "atomname", STR_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_atomname,
    NULL,
    {NULL, helptitle_atomname, asize(help_atomname), help_atomname}
};

/** Selection method data for \p pdbatomname selection keyword. */
gmx_ana_selmethod_t sm_pdbatomname = {
    "pdbatomname", STR_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    &check_pdbinfo,
    NULL,
    NULL,
    NULL,
    &evaluate_pdbatomname,
    NULL,
    {NULL, helptitle_atomname, asize(help_atomname), help_atomname}
};

/** Selection method data for \p atomtype selection keyword. */
gmx_ana_selmethod_t sm_atomtype = {
    "atomtype", STR_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    &check_atomtype,
    NULL,
    NULL,
    NULL,
    &evaluate_atomtype,
    NULL,
};

/** Selection method data for \p resname selection keyword. */
gmx_ana_selmethod_t sm_resname = {
    "resname", STR_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_resname,
    NULL,
};

/** Selection method data for \p chain selection keyword. */
gmx_ana_selmethod_t sm_insertcode = {
    "insertcode", STR_VALUE, SMETH_REQTOP | SMETH_CHARVAL,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_insertcode,
    NULL,
};

/** Selection method data for \p chain selection keyword. */
gmx_ana_selmethod_t sm_chain = {
    "chain", STR_VALUE, SMETH_REQTOP | SMETH_CHARVAL,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_chain,
    NULL,
};

/** Selection method data for \p mass selection keyword. */
gmx_ana_selmethod_t sm_mass = {
    "mass", REAL_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_mass,
    NULL,
};

/** Selection method data for \p charge selection keyword. */
gmx_ana_selmethod_t sm_charge = {
    "charge", REAL_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_charge,
    NULL,
};

/** Selection method data for \p chain selection keyword. */
gmx_ana_selmethod_t sm_altloc = {
    "altloc", STR_VALUE, SMETH_REQTOP | SMETH_CHARVAL,
    0, NULL,
    NULL,
    NULL,
    &check_pdbinfo,
    NULL,
    NULL,
    NULL,
    &evaluate_altloc,
    NULL,
};

/** Selection method data for \p occupancy selection keyword. */
gmx_ana_selmethod_t sm_occupancy = {
    "occupancy", REAL_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    &check_pdbinfo,
    NULL,
    NULL,
    NULL,
    &evaluate_occupancy,
    NULL,
};

/** Selection method data for \p betafactor selection keyword. */
gmx_ana_selmethod_t sm_betafactor = {
    "betafactor", REAL_VALUE, SMETH_REQTOP,
    0, NULL,
    NULL,
    NULL,
    &check_pdbinfo,
    NULL,
    NULL,
    NULL,
    &evaluate_betafactor,
    NULL,
};

/** Selection method data for \p x selection keyword. */
gmx_ana_selmethod_t sm_x = {
    "x", REAL_VALUE, SMETH_DYNAMIC,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_x,
};

/** Selection method data for \p y selection keyword. */
gmx_ana_selmethod_t sm_y = {
    "y", REAL_VALUE, SMETH_DYNAMIC,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_y,
};

/** Selection method data for \p z selection keyword. */
gmx_ana_selmethod_t sm_z = {
    "z", REAL_VALUE, SMETH_DYNAMIC,
    0, NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    NULL,
    &evaluate_z,
};

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Copies \p g to \p out->u.g.
 */
static void
evaluate_all(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
             gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    gmx_ana_index_copy(out->u.g, g, false);
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns an empty \p out->u.g.
 */
static void
evaluate_none(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
              gmx_ana_index_t * /* g */, gmx_ana_selvalue_t *out, void * /* data */)
{
    out->u.g->isize = 0;
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the indices for each atom in \p out->u.i.
 */
static void
evaluate_atomnr(t_topology * /* top */, t_trxframe * /* fr */, t_pbc * /* pbc */,
                gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

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
static void
evaluate_resnr(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
               gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;
    int  resind;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        resind      = top->atoms.atom[g->index[i]].resind;
        out->u.i[i] = top->atoms.resinfo[resind].nr;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the residue indices for each atom in \p out->u.i.
 */
static void
evaluate_resindex(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.i[i] = top->atoms.atom[g->index[i]].resind + 1;
    }
}

static void
check_molecules(t_topology *top, int /* npar */, gmx_ana_selparam_t * /* param */, void * /* data */)
{
    bool bOk;

    bOk = (top != NULL && top->mols.nr > 0);
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
static void
evaluate_molindex(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i, j;

    out->nr = g->isize;
    for (i = j = 0; i < g->isize; ++i)
    {
        while (top->mols.index[j + 1] <= g->index[i])
        {
            ++j;
        }
        out->u.i[i] = j + 1;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the atom name for each atom in \p out->u.s.
 */
static void
evaluate_atomname(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.s[i] = *top->atoms.atomname[g->index[i]];
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the PDB atom name for each atom in \p out->u.s.
 */
static void
evaluate_pdbatomname(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                     gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        char *s = top->atoms.pdbinfo[g->index[i]].atomnm;
        while (std::isspace(*s))
        {
            ++s;
        }
        out->u.s[i] = s;
    }
}

static void
check_atomtype(t_topology *top, int /* npar */, gmx_ana_selparam_t * /* param */, void * /* data */)
{
    bool bOk;

    bOk = (top != NULL && top->atoms.atomtype != NULL);
    if (!bOk)
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
static void
evaluate_atomtype(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                  gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.s[i] = *top->atoms.atomtype[g->index[i]];
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the residue name for each atom in \p out->u.s.
 */
static void
evaluate_resname(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                 gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;
    int  resind;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        resind      = top->atoms.atom[g->index[i]].resind;
        out->u.s[i] = *top->atoms.resinfo[resind].name;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the insertion code for each atom in \p out->u.s.
 */
static void
evaluate_insertcode(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                    gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;
    int  resind;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        resind         = top->atoms.atom[g->index[i]].resind;
        out->u.s[i][0] = top->atoms.resinfo[resind].ic;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the chain for each atom in \p out->u.s.
 */
static void
evaluate_chain(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
               gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;
    int  resind;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        resind         = top->atoms.atom[g->index[i]].resind;
        out->u.s[i][0] = top->atoms.resinfo[resind].chainid;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the mass for each atom in \p out->u.r.
 */
static void
evaluate_mass(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
              gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.r[i] = top->atoms.atom[g->index[i]].m;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the charge for each atom in \p out->u.r.
 */
static void
evaluate_charge(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.r[i] = top->atoms.atom[g->index[i]].q;
    }
}

static void
check_pdbinfo(t_topology *top, int /* npar */, gmx_ana_selparam_t * /* param */, void * /* data */)
{
    bool bOk;

    bOk = (top != NULL && top->atoms.pdbinfo != NULL);
    if (!bOk)
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
static void
evaluate_altloc(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.s[i][0] = top->atoms.pdbinfo[g->index[i]].altloc;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the occupancy numbers for each atom in \p out->u.r.
 * Segfaults if PDB info is not found in the topology.
 */
static void
evaluate_occupancy(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                   gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.r[i] = top->atoms.pdbinfo[g->index[i]].occup;
    }
}

/*!
 * See sel_updatefunc() for description of the parameters.
 * \p data is not used.
 *
 * Returns the B-factors for each atom in \p out->u.r.
 * Segfaults if PDB info is not found in the topology.
 */
static void
evaluate_betafactor(t_topology *top, t_trxframe * /* fr */, t_pbc * /* pbc */,
                    gmx_ana_index_t *g, gmx_ana_selvalue_t *out, void * /* data */)
{
    int  i;

    out->nr = g->isize;
    for (i = 0; i < g->isize; ++i)
    {
        out->u.r[i] = top->atoms.pdbinfo[g->index[i]].bfac;
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
static void
evaluate_coord(real out[], gmx_ana_pos_t *pos, int d)
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
static void
evaluate_x(t_topology * /*top*/, t_trxframe * /*fr*/, t_pbc * /*pbc*/,
           gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void * /*data*/)
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
static void
evaluate_y(t_topology * /*top*/, t_trxframe * /*fr*/, t_pbc * /*pbc*/,
           gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void * /*data*/)
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
static void
evaluate_z(t_topology * /*top*/, t_trxframe * /*fr*/, t_pbc * /*pbc*/,
           gmx_ana_pos_t *pos, gmx_ana_selvalue_t *out, void * /*data*/)
{
    out->nr = pos->count();
    evaluate_coord(out->u.r, pos, ZZ);
}
