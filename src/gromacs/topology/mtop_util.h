/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015,2016, by the GROMACS development team, led by
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
#ifndef GMX_TOPOLOGY_MTOP_UTIL_H
#define GMX_TOPOLOGY_MTOP_UTIL_H

#include <cstddef>

#include <vector>

#include "gromacs/utility/basedefinitions.h"

struct gmx_localtop_t;
struct gmx_moltype_t;
struct gmx_mtop_t;
struct t_atom;
struct t_atoms;
struct t_block;
struct t_ilist;
struct t_topology;

/* Should be called after generating or reading mtop,
 * to set some compute intesive variables to avoid
 * N^2 operations later on.
 */
void
gmx_mtop_finalize(gmx_mtop_t *mtop);

/* Counts the number of atoms of each type. State should be 0 for
 * state A and 1 for state B types.  typecount should have at
 * least mtop->ffparams.atnr elements.
 */
void
gmx_mtop_count_atomtypes(const gmx_mtop_t *mtop, int state, int typecount[]);

/* Returns the total number of charge groups in mtop */
int
ncg_mtop(const gmx_mtop_t *mtop);

/* Removes the charge groups, i.e. makes single atom charge groups, in mtop */
void gmx_mtop_remove_chargegroups(gmx_mtop_t *mtop);


/* Abstract data type for looking up atoms by global atom number */
typedef struct gmx_mtop_atomlookup *gmx_mtop_atomlookup_t;

/* Initialize atom lookup by global atom number */
gmx_mtop_atomlookup_t
gmx_mtop_atomlookup_init(const gmx_mtop_t *mtop);

/* As gmx_mtop_atomlookup_init, but optimized for atoms involved in settle */
gmx_mtop_atomlookup_t
gmx_mtop_atomlookup_settle_init(const gmx_mtop_t *mtop);

/* Destroy a gmx_mtop_atomlookup_t data structure */
void
gmx_mtop_atomlookup_destroy(gmx_mtop_atomlookup_t alook);


/* Returns a pointer to the t_atom struct belonging to atnr_global.
 * This can be an expensive operation, so if possible use
 * one of the atom loop constructs below.
 */
void
gmx_mtop_atomnr_to_atom(const gmx_mtop_atomlookup_t alook,
                        int                         atnr_global,
                        t_atom                    **atom);


/* Returns a pointer to the molecule interaction array ilist_mol[F_NRE]
 * and the local atom number in the molecule belonging to atnr_global.
 */
void
gmx_mtop_atomnr_to_ilist(const gmx_mtop_atomlookup_t alook,
                         int atnr_global,
                         t_ilist **ilist_mol, int *atnr_offset);


/* Returns the molecule block index
 * and the molecule number in the block
 * and the atom number offset for the atom indices in moltype
 * belonging to atnr_global.
 */
void
gmx_mtop_atomnr_to_molblock_ind(const gmx_mtop_atomlookup_t alook,
                                int atnr_global,
                                int *molb, int *molnr, int *atnr_mol);


/* Returns atom name, global resnr and residue name  of atom atnr_global */
void
gmx_mtop_atominfo_global(const gmx_mtop_t *mtop, int atnr_global,
                         char **atomname, int *resnr, char **resname);


/* Abstract type for atom loop over all atoms */
typedef struct gmx_mtop_atomloop_all *gmx_mtop_atomloop_all_t;

/* Initialize an atom loop over all atoms in the system.
 * The order of the atoms will be as in the state struct.
 * Only use this when you really need to loop over all atoms,
 * i.e. when you use groups which might differ per molecule,
 * otherwise use gmx_mtop_atomloop_block.
 */
gmx_mtop_atomloop_all_t
gmx_mtop_atomloop_all_init(const gmx_mtop_t *mtop);

/* Loop to the next atom.
 * When not at the end:
 *   returns TRUE and at_global,
 *   writes the global atom number in *at_global
 *   and sets the pointer atom to the t_atom struct of that atom.
 * When at the end, destroys aloop and returns FALSE.
 * Use as:
 * gmx_mtop_atomloop_all_t aloop;
 * aloop = gmx_mtop_atomloop_all_init(mtop)
 * while (gmx_mtop_atomloop_all_next(aloop,&at_global,&atom)) {
 *     ...
 * }
 */
gmx_bool
gmx_mtop_atomloop_all_next(gmx_mtop_atomloop_all_t aloop,
                           int *at_global, t_atom **atom);

/* Return the atomname, the residue number and residue name
 * of the current atom in the loop.
 */
void
gmx_mtop_atomloop_all_names(gmx_mtop_atomloop_all_t aloop,
                            char **atomname, int *resnr, char **resname);

/* Return the a pointer to the moltype struct of the current atom
 * in the loop and the atom number in the molecule.
 */
void
gmx_mtop_atomloop_all_moltype(gmx_mtop_atomloop_all_t aloop,
                              gmx_moltype_t **moltype, int *at_mol);


/* Abstract type for atom loop over atoms in all molecule blocks */
typedef struct gmx_mtop_atomloop_block *gmx_mtop_atomloop_block_t;

/* Initialize an atom loop over atoms in all molecule blocks the system.
 */
gmx_mtop_atomloop_block_t
gmx_mtop_atomloop_block_init(const gmx_mtop_t *mtop);

/* Loop to the next atom.
 * When not at the end:
 *   returns TRUE
 *   sets the pointer atom to the t_atom struct of that atom
 *   and return the number of molecules corresponding to this atom.
 * When at the end, destroys aloop and returns FALSE.
 * Use as:
 * gmx_mtop_atomloop_block_t aloop;
 * aloop = gmx_mtop_atomloop_block_init(mtop)
 * while (gmx_mtop_atomloop_block_next(aloop,&atom,&nmol)) {
 *     ...
 * }
 */
gmx_bool
gmx_mtop_atomloop_block_next(gmx_mtop_atomloop_block_t aloop,
                             t_atom **atom, int *nmol);


/* Abstract type for ilist loop over all ilists */
typedef struct gmx_mtop_ilistloop *gmx_mtop_ilistloop_t;

/* Initialize an ilist loop over all molecule types in the system. */
gmx_mtop_ilistloop_t
gmx_mtop_ilistloop_init(const gmx_mtop_t *mtop);


/* Loop to the next molecule,
 * When not at the end:
 *   returns TRUE and a pointer to the next array ilist_mol[F_NRE],
 *   writes the number of molecules for this ilist in *nmol.
 * When at the end, destroys iloop and returns FALSE.
 */
gmx_bool
gmx_mtop_ilistloop_next(gmx_mtop_ilistloop_t iloop,
                        t_ilist **ilist_mol, int *nmol);


/* Abstract type for ilist loop over all ilists of all molecules */
typedef struct gmx_mtop_ilistloop_all *gmx_mtop_ilistloop_all_t;

/* Initialize an ilist loop over all molecule types in the system.
 * Only use this when you really need to loop over all molecules,
 * i.e. when you use groups which might differ per molecule,
 * otherwise use gmx_mtop_ilistloop.
 */
gmx_mtop_ilistloop_all_t
gmx_mtop_ilistloop_all_init(const gmx_mtop_t *mtop);

/* Loop to the next molecule,
 * When not at the end:
 *   returns TRUE and a pointer to the next array ilist_mol[F_NRE],
 *   writes the atom offset which should be added to iatoms in atnr_offset.
 * When at the end, destroys iloop and returns FALSE.
 */
gmx_bool
gmx_mtop_ilistloop_all_next(gmx_mtop_ilistloop_all_t iloop,
                            t_ilist **ilist_mol, int *atnr_offset);


/* Returns the total number of interactions in the system of type ftype */
int
gmx_mtop_ftype_count(const gmx_mtop_t *mtop, int ftype);


/* Returns a charge group index for the whole system */
t_block
gmx_mtop_global_cgs(const gmx_mtop_t *mtop);


/* Returns a single t_atoms struct for the whole system */
t_atoms
gmx_mtop_global_atoms(const gmx_mtop_t *mtop);


/* Generate a 'local' topology for the whole system.
 * When feeEnergyInteractionsAtEnd == true, the free energy interactions will
 * be sorted to the end.
 */
gmx_localtop_t *
gmx_mtop_generate_local_top(const gmx_mtop_t *mtop, bool freeEnergyInteractionsAtEnd);


/* Converts a gmx_mtop_t struct to t_topology.
 * All memory relating only to mtop will be freed.
 */
t_topology
gmx_mtop_t_to_t_topology(gmx_mtop_t *mtop);

/*! \brief Get vector of atoms indices from topology
 *
 * This function returns the indices of all particles with type
 * eptAtom, that is shells, vsites etc. are left out.
 * \param[in]  mtop Molecular topology
 * \returns Vector that will be filled with the atom indices
 */
std::vector<size_t> get_atom_index(const gmx_mtop_t *mtop);

#endif
