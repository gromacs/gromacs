/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2025- The GROMACS Authors
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
 * Collection of topology preprocessing functions responsible for
 * all modifications of the topology during input pre-processing
 * when a part of the system is handled outside the classical
 * MM approximation (e.g. QM/MM or neural-network/MM).
 *
 * \author Lukas Müllender <lukas.muellender@gmail.com>
 */
#ifndef GMX_TOPOLOGY_EMBEDDED_SYSTEM_PREPROCESSING_H
#define GMX_TOPOLOGY_EMBEDDED_SYSTEM_PREPROCESSING_H

#include <set>
#include <vector>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_mtop_t;
class WarningHandler;

namespace gmx
{

class MDLogger;

/*! \internal
 * \brief Helper structure with indexes of broken bonds between embedded and MM
 * Used to determine and store pair of embedded and MM atoms between which chemical bond is broken
 */
struct LinkFrontier
{
    //! Global index of embedded atom at Frontier
    Index embedded;
    //! Global index of MM atom at Frontier
    Index mm;
};

/*! \brief Splits embedded atom containing molecules out of MM blocks in topology
 *
 * Modifies molblocks in topology \p mtop
 * \param[in,out] mtop topology to be modified
 * \param[in] embeddedIndices set with global indices of embedded atoms
 * \returns vector of flags for embedded atom-containing blocks in modified mtop
 */
std::vector<bool> splitEmbeddedBlocks(gmx_mtop_t* mtop, const std::set<int>& embeddedIndices);

/*! \brief Removes classical charges from embedded atoms and virtual sites
 *
 * Also removes charges from virtual sites built from embedded atoms only.
 * \param[in,out] mtop topology to be modified
 * \param[in] embeddedIndices set with global indices of embedded atoms
 * \param[in] isEmbeddedBlock vector with flags for embedded atom-containing blocks
 * \param[in] refQ reference total charge of the system, used for warning messages
 * \param[in] logger MDLogger for logging info about modifications
 * \param[in] wi WarningHandler for handling warnings
 * \returns vector of point charges for all atoms in the modified topology
 */
std::vector<real> removeEmbeddedClassicalCharges(gmx_mtop_t*              mtop,
                                                 const std::set<int>&     embeddedIndices,
                                                 const std::vector<bool>& isEmbeddedBlock,
                                                 real                     refQ,
                                                 const MDLogger&          logger,
                                                 WarningHandler*          wi);

/*! \brief Build exclusion list for non-bonded interactions between embedded atoms
 *
 * Adds embedded atoms to \c mtop->intermolecularExclusionGroup
 * \param[in,out] mtop topology to be modified
 * \param[in] embeddedIndices set with global indices of embedded atoms
 * \param[in] logger MDLogger for logging info about modifications
 */
void addEmbeddedNBExclusions(gmx_mtop_t* mtop, const std::set<int>& embeddedIndices, const MDLogger& logger);

/*! \brief Builds and returns a vector of atom numbers for all atoms in \p mtop.
 *
 * \param[in] mtop topology to be processed
 * \returns vector of atom numbers for all atoms
 */
std::vector<int> buildEmbeddedAtomNumbers(const gmx_mtop_t& mtop);

/*! \brief Modifies pairwise bonded interactions
 *
 * Removes any other pairwise bonded interactions between embedded atoms
 * Creates F_CONNBOND between embedded atoms
 * Any restraints and constraints will be kept
 * \param[in,out] mtop topology to be modified
 * \param[in] embeddedIndices set with global indices of embedded atoms
 * \param[in] isEmbeddedBlock vector with flags for embedded atom-containing blocks
 * \param[in] logger MDLogger for logging info about modifications
 */
void modifyEmbeddedTwoCenterInteractions(gmx_mtop_t*              mtop,
                                         const std::set<int>&     embeddedIndices,
                                         const std::vector<bool>& isEmbeddedBlock,
                                         const MDLogger&          logger);

/*! \brief Modifies three-centers interactions (i.e. Angles, Settles)
 *
 * Removes any other three-centers bonded interactions including 2 or more embedded atoms
 * Any restraints and constraints will be kept
 * Any F_SETTLE containing embedded atoms will be converted to the pair of F_CONNBONDS
 * \param[in,out] mtop topology to be modified
 * \param[in] embeddedIndices set with global indices of embedded atoms
 * \param[in] isEmbeddedBlock vector with flags for embedded atom-containing blocks
 * \param[in] logger MDLogger for logging info about modifications
 */
void modifyEmbeddedThreeCenterInteractions(gmx_mtop_t*              mtop,
                                           const std::set<int>&     embeddedIndices,
                                           const std::vector<bool>& isEmbeddedBlock,
                                           const MDLogger&          logger);

/*! \brief Modifies four-centers interactions
 *
 * Removes any other four-centers bonded interactions including 3 or more embedded atoms
 * Any restraints and constraints will be kept
 * \param[in,out] mtop topology to be modified
 * \param[in] embeddedIndices set with global indices of embedded atoms
 * \param[in] isEmbeddedBlock vector with flags for embedded atom-containing blocks
 * \param[in] logger MDLogger for logging info about modifications
 */
void modifyEmbeddedFourCenterInteractions(gmx_mtop_t*              mtop,
                                          const std::set<int>&     embeddedIndices,
                                          const std::vector<bool>& isEmbeddedBlock,
                                          const MDLogger&          logger);

/*! \brief Checks for constrained bonds within embedded subsystem
 *
 * Provides warnings via WarningHandler if any are found.
 * Separated from buildLinkFrontier for better modularity.
 * \param[in,out] mtop topology to be modified
 * \param[in] embeddedIndices set with global indices of embedded atoms
 * \param[in] isEmbeddedBlock vector with flags for embedded atom-containing blocks
 * \param[in] wi WarningHandler for handling warnings
 */
void checkConstrainedBonds(gmx_mtop_t*              mtop,
                           const std::set<int>&     embeddedIndices,
                           const std::vector<bool>& isEmbeddedBlock,
                           WarningHandler*          wi);

/*! \brief Builds link frontier vector with pairs of atoms indicting broken embedded - MM chemical bonds.
 *
 * Also performs search of constrained bonds within embedded subsystem.
 * \param[in,out] mtop topology to be modified
 * \param[in] embeddedIndices set with global indices of embedded atoms
 * \param[in] isEmbeddedBlock vector with flags for embedded atom-containing blocks
 * \param[in] logger MDLogger for logging info about modifications
 * \returns vector of link atom pairs
 */
std::vector<LinkFrontier> buildLinkFrontier(gmx_mtop_t*              mtop,
                                            const std::set<int>&     embeddedIndices,
                                            const std::vector<bool>& isEmbeddedBlock,
                                            const MDLogger&          logger);


} // namespace gmx

#endif // GMX_TOPOLOGY_EMBEDDED_SYSTEM_PREPROCESSING_H
