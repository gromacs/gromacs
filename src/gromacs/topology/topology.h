/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2011,2014,2015,2016,2018, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_TOPOLOGY_H
#define GMX_TOPOLOGY_TOPOLOGY_H

#include <cstdio>

#include <vector>

#include "gromacs/math/vectypes.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/symtab.h"

enum {
    egcTC,    egcENER,   egcACC, egcFREEZE,
    egcUser1, egcUser2,  egcVCM, egcCompressedX,
    egcORFIT, egcQMMM,
    egcNR
};
/* Names corresponding to groups */
extern const char *gtypes[egcNR+1];

/*! \brief Molecules type data: atoms, interactions and exclusions */
struct gmx_moltype_t
{
    /*! \brief Constructor */
    gmx_moltype_t();

    /*! \brief Destructor */
    ~gmx_moltype_t();

    /*! \brief Deleted copy assignment operator to avoid (not) freeing pointers */
    gmx_moltype_t &operator=(const gmx_moltype_t &) = delete;

    /*! \brief Default copy constructor */
    gmx_moltype_t(const gmx_moltype_t &) = default;

    char          **name;         /**< Name of the molecule type            */
    t_atoms         atoms;        /**< The atoms in this molecule           */
    t_ilist         ilist[F_NRE]; /**< Interaction list with local indices  */
    t_block         cgs;          /**< The charge groups                    */
    t_blocka        excls;        /**< The exclusions                       */
};

/*! \brief Block of molecules of the same type, used in gmx_mtop_t */
struct gmx_molblock_t
{
    int                    type = -1; /**< The molecule type index in mtop.moltype  */
    int                    nmol = 0;  /**< The number of molecules in this block    */
    std::vector<gmx::RVec> posres_xA; /**< Position restraint coordinates for top A */
    std::vector<gmx::RVec> posres_xB; /**< Position restraint coordinates for top B */
};

/*! \brief Indices for a gmx_molblock_t, derived from other gmx_mtop_t contents */
struct MoleculeBlockIndices
{
    int     numAtomsPerMolecule; /**< Number of atoms in a molecule in the block */
    int     globalAtomStart;     /**< Global atom index of the first atom in the block */
    int     globalAtomEnd;       /**< Global atom index + 1 of the last atom in the block */
    int     globalResidueStart;  /**< Global residue index of the first residue in the block */
    int     residueNumberStart;  /**< Residue numbers start from this value if the number of residues per molecule is <= maxres_renum */
    int     moleculeIndexStart;  /**< Global molecule indexing starts from this value */
};

typedef struct gmx_groups_t
{
    t_grps            grps[egcNR];  /* Groups of things                     */
    int               ngrpname;     /* Number of groupnames                 */
    char           ***grpname;      /* Names of the groups                  */
    int               ngrpnr[egcNR];
    unsigned char    *grpnr[egcNR]; /* Group numbers or NULL                */
} gmx_groups_t;

/* This macro gives the group number of group type egc for atom i.
 * This macro is useful, since the grpnr pointers are NULL
 * for group types that have all entries 0.
 */
#define ggrpnr(groups, egc, i) ((groups)->grpnr[egc] ? (groups)->grpnr[egc][i] : 0)

/* The global, complete system topology struct, based on molecule types.
 * This structure should contain no data that is O(natoms) in memory.
 *
 * TODO: Find a solution for ensuring that the derived data is in sync
 *       with the primary data, possibly by converting to a class.
 */
struct gmx_mtop_t
{
    /* Constructor */
    gmx_mtop_t();

    /* Destructor */
    ~gmx_mtop_t();

    char                      **name; /* Name of the topology                 */
    gmx_ffparams_t              ffparams;
    std::vector<gmx_moltype_t>  moltype;
    std::vector<gmx_molblock_t> molblock;
    gmx_bool                    bIntermolecularInteractions; /* Are there intermolecular
                                                              * interactions?            */
    t_ilist                    *intermolecular_ilist;        /* List of intermolecular interactions
                                                              * using system wide atom indices,
                                                              * either NULL or size F_NRE           */
    int              natoms;
    int              maxres_renum;                           /* Parameter for residue numbering      */
    int              maxresnr;                               /* The maximum residue number in moltype */
    t_atomtypes      atomtypes;                              /* Atomtype properties                  */
    gmx_groups_t     groups;                                 /* Groups of atoms for different purposes */
    t_symtab         symtab;                                 /* The symbol table                     */
    bool             haveMoleculeIndices;                    /* Tells whether we have valid molecule indices */

    /* Derived data */
    std::vector<MoleculeBlockIndices> moleculeBlockIndices;  /* Indices for each molblock entry for fast lookup of atom properties */
};

/* The mdrun node-local topology struct, completely written out */
typedef struct gmx_localtop_t
{
    t_idef        idef;         /* The interaction function definition  */
    t_atomtypes   atomtypes;    /* Atomtype properties                  */
    t_block       cgs;          /* The charge groups                    */
    t_blocka      excls;        /* The exclusions                       */
} gmx_localtop_t;

/* The old topology struct, completely written out, used in analysis tools */
typedef struct t_topology
{
    char          **name;                        /* Name of the topology                 */
    t_idef          idef;                        /* The interaction function definition  */
    t_atoms         atoms;                       /* The atoms                            */
    t_atomtypes     atomtypes;                   /* Atomtype properties                  */
    t_block         cgs;                         /* The charge groups                    */
    t_block         mols;                        /* The molecules                        */
    gmx_bool        bIntermolecularInteractions; /* Inter.mol. int. ?   */
    t_blocka        excls;                       /* The exclusions                       */
    t_symtab        symtab;                      /* The symbol table                     */
} t_topology;

void init_mtop(gmx_mtop_t *mtop);
void init_top(t_topology *top);
void done_gmx_groups_t(gmx_groups_t *g);
void done_top(t_topology *top);
// Frees both t_topology and gmx_mtop_t when the former has been created from
// the latter.
void done_top_mtop(t_topology *top, gmx_mtop_t *mtop);

bool gmx_mtop_has_masses(const gmx_mtop_t *mtop);
bool gmx_mtop_has_charges(const gmx_mtop_t *mtop);
bool gmx_mtop_has_atomtypes(const gmx_mtop_t *mtop);
bool gmx_mtop_has_pdbinfo(const gmx_mtop_t *mtop);

void pr_mtop(FILE *fp, int indent, const char *title, const gmx_mtop_t *mtop,
             gmx_bool bShowNumbers, gmx_bool bShowParameters);
void pr_top(FILE *fp, int indent, const char *title, const t_topology *top,
            gmx_bool bShowNumbers, gmx_bool bShowParameters);

void cmp_top(FILE *fp, const t_topology *t1, const t_topology *t2, real ftol, real abstol);
void cmp_groups(FILE *fp, const gmx_groups_t *g0, const gmx_groups_t *g1,
                int natoms0, int natoms1);

#endif
