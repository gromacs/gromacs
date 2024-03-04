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
#ifndef GMX_MDLIB_CONSTRAINT_GPU_HELPERS_H
#define GMX_MDLIB_CONSTRAINT_GPU_HELPERS_H

#include <numeric>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/listoflists.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
class InteractionDefinitions;
template<typename>
class ListOfLists;

//! Check if the coupled contraints work for LINCS for specific number of threads.
bool isNumCoupledConstraintsSupported(const gmx_mtop_t& mtop, int threadsPerBlock);

//! Helper type for discovering coupled constraints
struct AtomsAdjacencyListElement
{
    AtomsAdjacencyListElement(int indexOfSecondConstrainedAtom, int indexOfConstraint, int signFactor) :
        indexOfSecondConstrainedAtom_(indexOfSecondConstrainedAtom),
        indexOfConstraint_(indexOfConstraint),
        signFactor_(signFactor)
    {
    }
    //! The index of the other atom constrained to this atom.
    int indexOfSecondConstrainedAtom_;
    //! The index of this constraint in the container of constraints.
    int indexOfConstraint_;
    /*! \brief A multiplicative factor that indicates the relative
     * order of the atoms in the atom list.
     *
     * Used for computing the mass factor of this constraint
     * relative to any coupled constraints. */
    int signFactor_;
};


//! Constructs and returns an atom constraint adjacency list
gmx::ListOfLists<AtomsAdjacencyListElement> constructAtomsAdjacencyList(int numAtoms,
                                                                        gmx::ArrayRef<const int> iatoms);

/*! \brief Computes and returns how many constraints are coupled to each constraint
 *
 * Needed to introduce splits in data so that all coupled constraints will be computed in a
 * single GPU block. The position \p c of the vector \p numCoupledConstraints should have the number
 * of constraints that are coupled to a constraint \p c and are after \p c in the vector. Only
 * first index of the connected group of the constraints is needed later in the code, hence the
 * numCoupledConstraints vector is also used to keep track if the constrain was already counted.
 */
std::vector<int> countNumCoupledConstraints(gmx::ArrayRef<const int> iatoms,
                                            const gmx::ListOfLists<AtomsAdjacencyListElement>& atomAdjacencyList);

/*! \brief Helper function to go through constraints recursively.
 *
 *  For each constraint, counts the number of coupled constraints and stores the value in \p numCoupledConstraints array.
 *  This information is used to split the array of constraints between thread blocks on a GPU so there is no
 *  coupling between constraints from different thread blocks. After the \p numCoupledConstraints array is filled, the
 *  value \p numCoupledConstraints[c] should be equal to the number of constraints that are coupled to \p c and located
 *  after it in the constraints array.
 *
 * \param[in]     a                   Atom index.
 * \param[in,out] numCoupledConstraints  Indicates if the constraint was already counted and stores
 *                                    the number of constraints (i) connected to it and (ii) located
 *                                    after it in memory. This array is filled by this recursive function.
 *                                    For a set of coupled constraints, only for the first one in this list
 *                                    the number of consecutive coupled constraints is needed: if there is
 *                                    not enough space for this set of constraints in the thread block,
 *                                    the group has to be moved to the next one.
 * \param[in]     atomsAdjacencyList  Stores information about connections between atoms.
 */
int countCoupled(int                                                a,
                 gmx::ArrayRef<int>                                 numCoupledConstraints,
                 const gmx::ListOfLists<AtomsAdjacencyListElement>& atomsAdjacencyList);

//! \brief Indices of atoms in a water molecule
struct WaterMolecule
{
    //! Oxygen atom
    int ow1;
    //! First hydrogen atom
    int hw2;
    //! Second hydrogen atom
    int hw3;
};

//! Calculate total number of molecules to settle.
int computeTotalNumSettles(const gmx_mtop_t& mtop);

//! Composite data for settle initialization.
struct SettleWaterTopology
{
    //! Oxygen mass.
    real mO;
    //! Hydrogen mass.
    real mH;
    //! Distance H-O.
    real dOH;
    //! Distance H-H.
    real dHH;
};

//! Validates that settle can run and returns setup data.
SettleWaterTopology getSettleTopologyData(const gmx_mtop_t& mtop);

//! Local settle data.
struct LocalSettleData
{
    //! Number of settle interactions.
    int numSettle;
    //! NRAL
    int nral;
};

//! Compute numsettles from interaction def
LocalSettleData computeNumSettles(const InteractionDefinitions& idef);

//! Access to atom array from local settles.
gmx::ArrayRef<const int> localSettleAtoms(const InteractionDefinitions& idef);

#endif
