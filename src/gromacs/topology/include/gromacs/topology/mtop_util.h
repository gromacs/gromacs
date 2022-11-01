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
#ifndef GMX_TOPOLOGY_MTOP_UTIL_H
#define GMX_TOPOLOGY_MTOP_UTIL_H

#include <cstddef>

#include <array>
#include <vector>

#include "gromacs/topology/block.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/topology_enums.h"
#include "gromacs/utility/enumerationhelpers.h"

struct gmx_mtop_t;
struct gmx_moltype_t;
enum class ParticleType : int;
struct gmx_localtop_t;
struct t_atoms;
struct t_symtab;
struct t_topology;

namespace gmx
{
template<typename>
class ArrayRef;
} // namespace gmx

// TODO All of the functions taking a const gmx_mtop * are deprecated
// and should be replaced by versions taking const gmx_mtop & when
// their callers are refactored similarly.

/* Counts the number of atoms of each type. State should be 0 for
 * state A and 1 for state B types.  typecount should have at
 * least mtop->ffparams.atnr elements.
 */
void gmx_mtop_count_atomtypes(const gmx_mtop_t& mtop, int state, int typecount[]);

/*!\brief Returns the total number of molecules in mtop
 *
 * \param[in] mtop  The global topology
 */
int gmx_mtop_num_molecules(const gmx_mtop_t& mtop);

/* Returns the total number of residues in mtop. */
int gmx_mtop_nres(const gmx_mtop_t& mtop);

/* Returns the total number of interactions in the system of type ftype */
int gmx_mtop_ftype_count(const gmx_mtop_t& mtop, int ftype);

/* Returns the total number of interactions in the system with all interaction flags that are set in \p if_flags set */
int gmx_mtop_interaction_count(const gmx_mtop_t& mtop, int unsigned if_flags);

/* Returns the count of atoms for each particle type */
gmx::EnumerationArray<ParticleType, int> gmx_mtop_particletype_count(const gmx_mtop_t& mtop);

/* Returns a single t_atoms struct for the whole system */
t_atoms gmx_mtop_global_atoms(const gmx_mtop_t& mtop);

/*! \brief
 * Populate a 'local' topology for the whole system.
 *
 * When freeEnergyInteractionsAtEnd == true, the free energy interactions will
 * be sorted to the end.
 *
 * \param[in]     mtop                        The global topology used to populate the local one.
 * \param[in,out] top                         New local topology populated from global \p mtop.
 * \param[in]     freeEnergyInteractionsAtEnd If free energy interactions will be sorted.
 */
void gmx_mtop_generate_local_top(const gmx_mtop_t& mtop, gmx_localtop_t* top, bool freeEnergyInteractionsAtEnd);


/*!\brief Creates and returns a struct with begin/end atom indices of all molecules
 *
 * \param[in] mtop  The global topology
 * \returns A RangePartitioning object with numBlocks() equal to the number
 * of molecules and atom indices such that molecule m contains atoms a with:
 * index[m] <= a < index[m+1].
 */
gmx::RangePartitioning gmx_mtop_molecules(const gmx_mtop_t& mtop);

/*! \brief
 * Returns the index range from residue begin to end for each residue in a molecule block.
 *
 * Note that residues will always have consecutive atoms numbers internally.
 *
 * \param[in] moltype  Molecule Type to parse for start and end.
 * \returns Vector of ranges for all residues.
 */
std::vector<gmx::Range<int>> atomRangeOfEachResidue(const gmx_moltype_t& moltype);

/* Converts a gmx_mtop_t struct to t_topology.
 *
 * If the lifetime of the returned topology should be longer than that
 * of mtop, your need to pass freeMtop==true.
 * If freeMTop == true, memory related to mtop will be freed so that done_top()
 * on the result value will free all memory.
 * If freeMTop == false, mtop and the return value will share some of their
 * memory, and there is currently no way to consistently free all the memory.
 */
t_topology gmx_mtop_t_to_t_topology(gmx_mtop_t* mtop, bool freeMTop);

/*! \brief Get vector of atoms indices from topology
 *
 * This function returns the indices of all particles with type
 * eptAtom, that is shells, vsites etc. are left out.
 * \param[in]  mtop Molecular topology
 * \returns Vector that will be filled with the atom indices
 */
std::vector<int> get_atom_index(const gmx_mtop_t& mtop);

/*! \brief Converts a t_atoms struct to an mtop struct
 *
 * All pointers contained in \p atoms will be copied into \p mtop.
 * Note that this will produce one moleculetype encompassing the whole system.
 *
 * \param[in]  symtab  The symbol table
 * \param[in]  name    Pointer to the name for the topology
 * \param[in]  atoms   The atoms to convert
 * \param[out] mtop    The molecular topology output containing atoms.
 */
void convertAtomsToMtop(t_symtab* symtab, char** name, t_atoms* atoms, gmx_mtop_t* mtop);

//! Checks and returns whether non-bonded interactions are perturbed for free-energy calculations
bool haveFepPerturbedNBInteractions(const gmx_mtop_t& mtop);

//! Checks whether masses are perturbed for free-energy calculations
bool haveFepPerturbedMasses(const gmx_mtop_t& mtop);

//! Checks whether masses are perturbed for free-energy calculations in SETTLE interactions
bool haveFepPerturbedMassesInSettles(const gmx_mtop_t& mtop);

//! Checks whether constraints are perturbed for free-energy calculations
bool havePerturbedConstraints(const gmx_mtop_t& mtop);

#endif
