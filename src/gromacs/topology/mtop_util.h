/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014,2015,2016,2018,2019,2020, by the GROMACS development team, led by
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

#include <array>
#include <vector>

#include "gromacs/topology/topology.h"
#include "gromacs/utility/basedefinitions.h"

struct gmx_localtop_t;
struct t_atom;
struct t_atoms;
struct t_block;
struct t_symtab;
enum struct GmxQmmmMode;

// TODO All of the functions taking a const gmx_mtop * are deprecated
// and should be replaced by versions taking const gmx_mtop & when
// their callers are refactored similarly.

/* Should be called after generating or reading mtop,
 * to set some compute intesive variables to avoid
 * N^2 operations later on.
 */
void gmx_mtop_finalize(gmx_mtop_t* mtop);

/* Counts the number of atoms of each type. State should be 0 for
 * state A and 1 for state B types.  typecount should have at
 * least mtop->ffparams.atnr elements.
 */
void gmx_mtop_count_atomtypes(const gmx_mtop_t* mtop, int state, int typecount[]);

/*!\brief Returns the total number of molecules in mtop
 *
 * \param[in] mtop  The global topology
 */
int gmx_mtop_num_molecules(const gmx_mtop_t& mtop);

/* Returns the total number of residues in mtop. */
int gmx_mtop_nres(const gmx_mtop_t* mtop);

class AtomIterator;

//! Proxy object returned from AtomIterator
class AtomProxy
{
public:
    //! Default constructor.
    AtomProxy(const AtomIterator* it) : it_(it) {}
    //! Access current global atom number.
    int globalAtomNumber() const;
    //! Access current t_atom struct.
    const t_atom& atom() const;
    //! Access current name of the atom.
    const char* atomName() const;
    //! Access current name of the residue the atom is in.
    const char* residueName() const;
    //! Access current residue number.
    int residueNumber() const;
    //! Access current molecule type.
    const gmx_moltype_t& moleculeType() const;
    //! Access the position of the current atom in the molecule.
    int atomNumberInMol() const;

private:
    const AtomIterator* it_;
};

//! Wrapper around proxy object to implement operator->
template<typename T>
class ProxyPtr
{
public:
    //! Construct with proxy object.
    ProxyPtr(T t) : t_(t) {}
    //! Member of pointer operator.
    T* operator->() { return &t_; }

private:
    T t_;
};

/*! \brief
 * Object that allows looping over all atoms in an mtop.
 */
class AtomIterator
{
public:
    //! Construct from topology and optionalally a global atom number.
    explicit AtomIterator(const gmx_mtop_t& mtop, int globalAtomNumber = 0);

    //! Prefix increment.
    AtomIterator& operator++();
    //! Postfix increment.
    AtomIterator operator++(int);

    //! Equality comparison.
    bool operator==(const AtomIterator& o) const;
    //! Non-equal comparison.
    bool operator!=(const AtomIterator& o) const;

    //! Dereference operator. Returns proxy.
    AtomProxy operator*() const { return { this }; }
    //! Member of pointer operator.
    ProxyPtr<AtomProxy> operator->() const { return { this }; }

private:
    //! Global topology.
    const gmx_mtop_t* mtop_;
    //! Current molecule block.
    size_t mblock_;
    //! The atoms of the current molecule.
    const t_atoms* atoms_;
    //! The current molecule.
    int currentMolecule_;
    //! Current highest number for residues.
    int highestResidueNumber_;
    //! Current local atom number.
    int localAtomNumber_;
    //! Global current atom number.
    int globalAtomNumber_;

    friend class AtomProxy;
};

//! Range over all atoms of topology.
class AtomRange
{
public:
    //! Default constructor.
    explicit AtomRange(const gmx_mtop_t& mtop) : begin_(mtop), end_(mtop, mtop.natoms) {}
    //! Iterator to begin of range.
    AtomIterator& begin() { return begin_; }
    //! Iterator to end of range.
    AtomIterator& end() { return end_; }

private:
    AtomIterator begin_, end_;
};

/* Abstract type for atom loop over atoms in all molecule blocks */
typedef struct gmx_mtop_atomloop_block* gmx_mtop_atomloop_block_t;

/* Initialize an atom loop over atoms in all molecule blocks the system.
 */
gmx_mtop_atomloop_block_t gmx_mtop_atomloop_block_init(const gmx_mtop_t* mtop);

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
gmx_bool gmx_mtop_atomloop_block_next(gmx_mtop_atomloop_block_t aloop, const t_atom** atom, int* nmol);


/* Abstract type for ilist loop over all ilists */
typedef struct gmx_mtop_ilistloop* gmx_mtop_ilistloop_t;

/* Initialize an ilist loop over all molecule types in the system. */
gmx_mtop_ilistloop_t gmx_mtop_ilistloop_init(const gmx_mtop_t* mtop);

/* Initialize an ilist loop over all molecule types in the system. */
gmx_mtop_ilistloop_t gmx_mtop_ilistloop_init(const gmx_mtop_t& mtop);

/* Loop to the next molecule,
 * When not at the end:
 *   returns a valid pointer to the next array ilist_mol[F_NRE],
 *   writes the number of molecules for this ilist in *nmol.
 * When at the end, destroys iloop and returns nullptr.
 */
const InteractionLists* gmx_mtop_ilistloop_next(gmx_mtop_ilistloop_t iloop, int* nmol);

/* Abstract type for ilist loop over all ilists of all molecules */
typedef struct gmx_mtop_ilistloop_all* gmx_mtop_ilistloop_all_t;

/* Initialize an ilist loop over all molecule types in the system.
 * Only use this when you really need to loop over all molecules,
 * i.e. when you use groups which might differ per molecule,
 * otherwise use gmx_mtop_ilistloop.
 */
gmx_mtop_ilistloop_all_t gmx_mtop_ilistloop_all_init(const gmx_mtop_t* mtop);

/* Loop to the next molecule,
 * When not at the end:
 *   returns a valid pointer to the next array ilist_mol[F_NRE],
 *   writes the atom offset which should be added to iatoms in atnr_offset.
 * When at the end, destroys iloop and returns nullptr.
 */
const InteractionLists* gmx_mtop_ilistloop_all_next(gmx_mtop_ilistloop_all_t iloop, int* atnr_offset);


/* Returns the total number of interactions in the system of type ftype */
int gmx_mtop_ftype_count(const gmx_mtop_t* mtop, int ftype);

/* Returns the total number of interactions in the system of type ftype */
int gmx_mtop_ftype_count(const gmx_mtop_t& mtop, int ftype);

/* Returns the total number of interactions in the system with all interaction flags that are set in \p if_flags set */
int gmx_mtop_interaction_count(const gmx_mtop_t& mtop, int unsigned if_flags);

/* Returns the count of atoms for each particle type */
std::array<int, eptNR> gmx_mtop_particletype_count(const gmx_mtop_t& mtop);

/* Returns a single t_atoms struct for the whole system */
t_atoms gmx_mtop_global_atoms(const gmx_mtop_t* mtop);


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
std::vector<int> get_atom_index(const gmx_mtop_t* mtop);

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

#endif
