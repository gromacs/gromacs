/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2022- The GROMACS Authors
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
/*! \libinternal \file
 *
 * \brief This file contains functions to loop over topology contents.
 *
 * \author Berk Hess <hess@kth.se>
 * \author Paul Bauer <paul.bauer.q@gmail.com>
 * \inlibraryapi
 * \ingroup module_mtop
 */
#ifndef GMX_TOPOLOGY_MTOP_ATOMLOOPS_H
#define GMX_TOPOLOGY_MTOP_ATOMLOOPS_H

#include "external/boost/stl_interfaces/iterator_interface.hpp"

#include "gromacs/topology/topology.h"
#include "gromacs/topology/topology_enums.h"

namespace gmx
{
template<typename>
class ArrayRef;
} // namespace gmx

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

/*! \brief
 * Object that allows looping over all atoms in an mtop.
 */
class AtomIterator :
    public gmx::boost::stl_interfaces::proxy_iterator_interface<AtomIterator, std::forward_iterator_tag, t_atom, AtomProxy>
{
    using Base =
            gmx::boost::stl_interfaces::proxy_iterator_interface<AtomIterator, std::forward_iterator_tag, t_atom, AtomProxy>;

public:
    //! Construct from topology and optionalally a global atom number.
    explicit AtomIterator(const gmx_mtop_t& mtop, int globalAtomNumber = 0);

    //! Prefix increment.
    AtomIterator& operator++();
    using Base::  operator++;

    //! Equality comparison.
    bool operator==(const AtomIterator& o) const;

    //! Dereference operator. Returns proxy.
    AtomProxy operator*() const { return { this }; }

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

class IListIterator;

//! Proxy object returned from IListIterator
class IListProxy
{
public:
    //! Default constructor.
    IListProxy(const IListIterator* it) : it_(it) {}
    //! Access current global atom number.
    const InteractionLists& list() const;
    //! Access current molecule.
    int nmol() const;

private:
    const IListIterator* it_;
};

/*! \brief
 * Object that allows looping over all atoms in an mtop.
 */
class IListIterator :
    public gmx::boost::stl_interfaces::proxy_iterator_interface<IListIterator, std::forward_iterator_tag, InteractionLists, IListProxy>
{
    using Base =
            gmx::boost::stl_interfaces::proxy_iterator_interface<IListIterator, std::forward_iterator_tag, InteractionLists, IListProxy>;

public:
    //! Construct from topology.
    explicit IListIterator(const gmx_mtop_t& mtop, size_t mblock = 0);

    //! Prefix increment.
    IListIterator& operator++();
    using Base::   operator++;

    //! Equality comparison.
    bool operator==(const IListIterator& o) const;

    //! Dereference operator. Returns proxy.
    IListProxy operator*() const { return { this }; }

private:
    //! Global topology.
    const gmx_mtop_t* mtop_;
    //! Index of molecule block corresponding to the current location.
    size_t mblock_;

    friend class IListProxy;
};


/*! \brief
 * Range over all interaction lists of topology.
 *
 * Includes the intermolecular interactions as the final element in the
 * range if present.
 */
class IListRange
{
public:
    //! Default constructor.
    explicit IListRange(const gmx_mtop_t& mtop);
    //! Iterator to begin of range.
    IListIterator& begin() { return begin_; }
    //! Iterator to end of range.
    IListIterator& end() { return end_; }

private:
    IListIterator begin_, end_;
};

/* Abstract type for atom loop over atoms in all molecule blocks */
typedef struct gmx_mtop_atomloop_block* gmx_mtop_atomloop_block_t;

/* Initialize an atom loop over atoms in all molecule blocks the system.
 */
gmx_mtop_atomloop_block_t gmx_mtop_atomloop_block_init(const gmx_mtop_t& mtop);

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
bool gmx_mtop_atomloop_block_next(gmx_mtop_atomloop_block_t aloop, const t_atom** atom, int* nmol);

#endif
