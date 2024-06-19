/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2021- The GROMACS Authors
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
 * \brief This file makes declarations used for building
 * the reverse topology
 *
 * \inlibraryapi
 * \ingroup module_domdec
 */

#ifndef GMX_DOMDEC_REVERSETOPOLOGY_H
#define GMX_DOMDEC_REVERSETOPOLOGY_H

#include <cstdio>

#include <memory>
#include <vector>

#include "gromacs/mdlib/vsite.h"
#include "gromacs/topology/idef.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/listoflists.h"

struct gmx_domdec_t;
struct gmx_ffparams_t;
struct gmx_mtop_t;
struct t_atoms;
struct t_inputrec;
struct ReverseTopOptions;

namespace gmx
{
class VirtualSitesHandler;
enum class DDBondedChecking : bool;
} // namespace gmx

//! Options for linking atoms in make_reverse_ilist
enum class AtomLinkRule
{
    FirstAtom,        //!< Link all interactions to the first atom in the atom list
    AllAtomsInBondeds //!< Link bonded interactions to all atoms involved, don't link vsites
};

struct MolblockIndices
{
    int a_start;
    int a_end;
    int natoms_mol;
    int type;
};

struct reverse_ilist_t
{
    std::vector<int> index;              /* Index for each atom into il          */
    std::vector<int> il;                 /* ftype|type|a0|...|an|ftype|...       */
    int              numAtomsInMolecule; /* The number of atoms in this molecule */
};

/*! \internal \brief Struct for thread local work data for local topology generation */
struct thread_work_t
{
    /*! \brief Constructor
     *
     * \param[in] ffparams  The interaction parameters, the lifetime of the created object should not exceed the lifetime of the passed parameters
     */
    thread_work_t(const gmx_ffparams_t& ffparams) : idef(ffparams) {}

    InteractionDefinitions         idef;               /**< Partial local topology */
    std::unique_ptr<gmx::VsitePbc> vsitePbc = nullptr; /**< vsite PBC structure */
    int numBondedInteractions               = 0; /**< The number of bonded interactions observed */
    gmx::ListOfLists<int> excl;                  /**< List of exclusions */
};

/*! \internal \brief Options for setting up gmx_reverse_top_t */
struct ReverseTopOptions
{
    //! Constructor, constraints and settles are not including with a single argument
    ReverseTopOptions(gmx::DDBondedChecking ddBondedChecking,
                      bool                  includeConstraints = false,
                      bool                  includeSettles     = false) :
        ddBondedChecking_(ddBondedChecking),
        includeConstraints_(includeConstraints),
        includeSettles_(includeSettles)
    {
    }

    //! \brief For which bonded interactions to check assignments
    const gmx::DDBondedChecking ddBondedChecking_;
    //! \brief Whether constraints are stored in this reverse top
    const bool includeConstraints_;
    //! \brief Whether settles are stored in this reverse top
    const bool includeSettles_;
};

//! \internal \brief Reverse topology class
class gmx_reverse_top_t
{
public:
    //! Constructor
    gmx_reverse_top_t(const gmx_mtop_t& mtop, bool useFreeEnergy, const ReverseTopOptions& reverseTopOptions);
    //! Destructor
    ~gmx_reverse_top_t();

    //! Gets the options that configured the construction
    const ReverseTopOptions& options() const;

    //! Gets the interaction list for the given molecule type
    const reverse_ilist_t& interactionListForMoleculeType(int moleculeType) const;

    //! Returns the molecule block indices
    gmx::ArrayRef<const MolblockIndices> molblockIndices() const;
    //! Returns whether the reverse topology describes intermolecular interactions
    bool hasIntermolecularInteractions() const;
    //! Gets the interaction list for any intermolecular interactions
    const reverse_ilist_t& interactionListForIntermolecularInteractions() const;
    //! Returns whether the reverse topology describes interatomic interactions
    bool hasInterAtomicInteractions() const;
    //! Returns whether there are interactions of type F_POSRES and/or F_FBPOSRES
    bool hasPositionRestraints() const;
    //! Returns the per-thread working structures for making the local topology
    gmx::ArrayRef<thread_work_t> threadWorkObjects() const;
    //! Returns whether the local topology listed-forces interactions should be sorted
    bool doListedForcesSorting() const;

    //! Private implementation definition
    struct Impl;
    //! Private implementation declaration
    std::unique_ptr<Impl> impl_;
};

/*! \brief Returns the number of atom entries for il in gmx_reverse_top_t */
int nral_rt(int ftype);

/*! \brief Return whether interactions of type \p ftype need to be assigned exactly once */
bool dd_check_ftype(int ftype, const ReverseTopOptions& rtOptions);

/*! \internal
 *  \brief Molecular topology indices of a global molecule a global atom belongs to
 */
struct MolecularTopologyAtomIndices
{
    //! The index of the molecule block
    int blockIndex;
    //! The molecule type
    int moleculeType;
    //! The index of the molecule in the block
    int moleculeIndex;
    //! The index of the atom in the molecule
    int atomIndex;
};

//! Return global topology molecule information for global atom index \p globalAtomIndex
MolecularTopologyAtomIndices globalAtomIndexToMoltypeIndices(gmx::ArrayRef<const MolblockIndices> molblockIndices,
                                                             int globalAtomIndex);

/*! \brief Make the reverse ilist: a list of bonded interactions linked to atoms */
void make_reverse_ilist(const InteractionLists&  ilist,
                        const t_atoms*           atoms,
                        const ReverseTopOptions& rtOptions,
                        AtomLinkRule             atomLinkRule,
                        reverse_ilist_t*         ril_mt);

/*! \brief Generate and store the reverse topology */
void dd_make_reverse_top(FILE*                           fplog,
                         gmx_domdec_t*                   dd,
                         const gmx_mtop_t&               mtop,
                         const gmx::VirtualSitesHandler* vsite,
                         const t_inputrec&               inputrec,
                         gmx::DDBondedChecking           ddBondedChecking);

#endif
