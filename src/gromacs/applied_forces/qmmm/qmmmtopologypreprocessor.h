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
/*! \internal \file
 * \brief
 * QMMMTopologyPrepocessor class responsible for
 * all modificatios of the topology during input pre-processing
 *
 * \author Dmitry Morozov <dmitry.morozov@jyu.fi>
 * \ingroup module_applied_forces
 */
#ifndef GMX_APPLIED_FORCES_QMMMTOPOLOGYPREPROCESSOR_H
#define GMX_APPLIED_FORCES_QMMMTOPOLOGYPREPROCESSOR_H

#include <set>
#include <string>
#include <vector>

#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

#include "qmmmtypes.h"

struct gmx_mtop_t;
class WarningHandler;

namespace gmx
{

class MDLogger;

/*! \internal
 * \brief Contains various information about topology modifications
 * Used for statistics during topology pre-processing within QMMMTopologyPreprocessor class
 */
struct QMMMTopologyInfo
{
    //! Total number of MM atoms
    int numMMAtoms = 0;
    //! Total number of QM atoms
    int numQMAtoms = 0;
    //! Total remaining charge of MM part
    real remainingMMCharge = 0.0;
    //! Total classical charge removed from QM atoms
    real totalClassicalChargeOfQMAtoms = 0.0;
    //! Total number of Non-bonded (LJ) exclusions made for QM-QM interactions
    int numExclusionsMade = 0;
    //! Total number of removed classical Bonds between QM-QM atoms
    int numBondsRemoved = 0;
    //! Total number of removed classical Angles between QM-QM atoms
    int numAnglesRemoved = 0;
    //! Total number of removed classical Dihedrals between QM-QM atoms
    int numDihedralsRemoved = 0;
    //! Total number of removed F_SETTLE between QM-QM atoms
    int numSettleRemoved = 0;
    //! Total number of empty chemical bonds (F_CONNBONDS) added between QM-QM atoms
    int numConnBondsAdded = 0;
    //! Total number of virtual sites, that consisting of QM atoms only, which charge has been removed
    int numVirtualSitesModified = 0;
    //! Total number of constrained bonds within QM subsystem
    int numConstrainedBondsInQMSubsystem = 0;
    //! Total number of broken bonds between QM and MM atoms (Link Frontier)
    int numLinkBonds = 0;
};

/*! \internal
 * \brief Class implementing gmx_mtop_t QMMM modifications during preprocessing
 */
class QMMMTopologyPreprocessor
{
public:
    /*! \brief Constructor for QMMMTopologyPreprocessor from its parameters
     *
     * \param[in] qmIndices Array with global indicies of QM atoms
     */
    QMMMTopologyPreprocessor(ArrayRef<const Index> qmIndices);

    /*! \brief Pocesses mtop topology and prepares atomNumbers_ and linkFrontier_ vectors
     * Builds topInfo_ containing information about topology modifications
     *
     * \param[in,out] mtop Topology that needs to be modified
     * \param[in] refQ Reference total charge of the system, used for warnings
     * \param[in] logger MDLogger for logging info about modifications
     * \param[in] wi WarningHandler for handling warnings
     */
    void preprocess(gmx_mtop_t* mtop, real refQ, const MDLogger& logger, WarningHandler* wi);

    //! \brief Returns data about modifications made via QMMMTopologyInfo
    const QMMMTopologyInfo& topInfo() const;

    //! \brief Returns view of atomic numbers for all atoms in the processed topology
    ArrayRef<const int> atomNumbers() const;

    //! \brief Returns view of point charges for all atoms in the processed topology
    ArrayRef<const real> atomCharges() const;

    //! \brief Returns view of the whole Link Frontier for the processed topology
    ArrayRef<const LinkFrontier> linkFrontier() const;

private:
    /*! \brief Global indices of QM atoms;
     * The dominant operation is search and we also expect the set of qm atoms to be very small
     * relative to the rest, so set should outperform unordered set, i.e. unsorted std::vector.
     */
    std::set<int> qmIndices_;
    //! Vector with atom numbers for the whole system
    std::vector<int> atomNumbers_;
    //! Vector with atom point charges for the whole system
    std::vector<real> atomCharges_;
    //! Vector with pairs of indices defining broken bonds in QMMM
    std::vector<LinkFrontier> linkFrontier_;
    //! Structure with information about modifications made
    QMMMTopologyInfo topInfo_;
};

/*! \brief Splits QM containing molecules out of MM blocks in topology
 *
 * Modifies molblocks in topology \p mtop
 * \param[in,out] mtop topology to be modified
 * \param[in] qmIndices set with global indices of QM atoms
 * \returns vector of flags for QM-containing blocks in modified mtop
 */
std::vector<bool> splitQMBlocks(gmx_mtop_t* mtop, const std::set<int>& qmIndices);

/*! \brief Removes classical charges from QM atoms and virtual sites
 *
 * Provides data about removed charge via \p topInfo.
 * Also removes charges from virtual sites built from QM atoms only.
 * \param[in,out] mtop topology to be modified
 * \param[in] qmIndices set with global indices of QM atoms
 * \param[in] bQMBlock vector with flags for QM-containing blocks
 * \param[in,out] topInfo structure with information about topology modifications
 * \param[in] refQ reference total charge of the system, used for warning messages
 * \param[in] logger MDLogger for logging info about modifications
 * \param[in] wi WarningHandler for handling warnings
 * \returns vector of point charges for all atoms in the modified topology
 */
std::vector<real> removeQMClassicalCharges(gmx_mtop_t*              mtop,
                                           const std::set<int>&     qmIndices,
                                           const std::vector<bool>& bQMBlock,
                                           QMMMTopologyInfo&        topInfo,
                                           real                     refQ,
                                           const MDLogger&          logger,
                                           WarningHandler*          wi);

/*! \brief Build exclusion list for non-bonded interactions between QM atoms
 *
 * Adds QM atoms to \c mtop->intermolecularExclusionGroup
 * \param[in,out] mtop topology to be modified
 * \param[in] qmIndices set with global indices of QM atoms
 * \param[in,out] topInfo structure with information about topology modifications
 * \param[in] logger MDLogger for logging info about modifications
 */
void addQMLJExclusions(gmx_mtop_t*          mtop,
                       const std::set<int>& qmIndices,
                       QMMMTopologyInfo&    topInfo,
                       const MDLogger&      logger);

/*! \brief Builds and returns a vector of atom numbers for all atoms in \p mtop.
 *
 * Provides data about total number of QM and MM atoms via topInfo_
 * \param[in] mtop topology to be processed
 * \returns vector of atom numbers for all atoms
 */
std::vector<int> buildQMMMAtomNumbers(const gmx_mtop_t& mtop);

/*! \brief Modifies pairwise bonded interactions
 *
 * Removes any other pairwise bonded interactions between QM-QM atoms
 * Creates F_CONNBOND between QM atoms
 * Any restraints and constraints will be kept
 * Provides data about modifications via \p topInfo
 * \param[in,out] mtop topology to be modified
 * \param[in] qmIndices set with global indices of QM atoms
 * \param[in] bQMBlock vector with flags for QM-containing blocks
 * \param[in,out] topInfo structure with information about topology modifications
 * \param[in] logger MDLogger for logging info about modifications
 */
void modifyQMMMTwoCenterInteractions(gmx_mtop_t*              mtop,
                                     const std::set<int>&     qmIndices,
                                     const std::vector<bool>& bQMBlock,
                                     QMMMTopologyInfo&        topInfo,
                                     const MDLogger&          logger);

/*! \brief Modifies three-centers interactions (i.e. Angles, Settles)
 *
 * Removes any other three-centers bonded interactions including 2 or more QM atoms
 * Any restraints and constraints will be kept
 * Any F_SETTLE containing QM atoms will be converted to the pair of F_CONNBONDS
 * Provides data about modifications via \p topInfo
 * \param[in,out] mtop topology to be modified
 * \param[in] qmIndices set with global indices of QM atoms
 * \param[in] bQMBlock vector with flags for QM-containing blocks
 * \param[in,out] topInfo structure with information about topology modifications
 * \param[in] logger MDLogger for logging info about modifications
 */
void modifyQMMMThreeCenterInteractions(gmx_mtop_t*              mtop,
                                       const std::set<int>&     qmIndices,
                                       const std::vector<bool>& bQMBlock,
                                       QMMMTopologyInfo&        topInfo,
                                       const MDLogger&          logger);

/*! \brief Modifies four-centers interactions
 *
 * Removes any other four-centers bonded interactions including 3 or more QM atoms
 * Any restraints and constraints will be kept
 * Provides data about modifications via \p topInfo
 * \param[in,out] mtop topology to be modified
 * \param[in] qmIndices set with global indices of QM atoms
 * \param[in] bQMBlock vector with flags for QM-containing blocks
 * \param[in,out] topInfo structure with information about topology modifications
 * \param[in] logger MDLogger for logging info about modifications
 */
void modifyQMMMFourCenterInteractions(gmx_mtop_t*              mtop,
                                      const std::set<int>&     qmIndices,
                                      const std::vector<bool>& bQMBlock,
                                      QMMMTopologyInfo&        topInfo,
                                      const MDLogger&          logger);

/*! \brief Checks for constrained bonds within QM subsystem
 *
 * Provides warnings via WarningHandler if any are found.
 * Separated from buildQMMMLink for better modularity.
 * \param[in,out] mtop topology to be modified
 * \param[in] qmIndices set with global indices of QM atoms
 * \param[in] bQMBlock vector with flags for QM-containing blocks
 * \param[in,out] topInfo structure with information about topology modifications
 * \param[in] wi WarningHandler for handling warnings
 */
void checkConstrainedBonds(gmx_mtop_t*              mtop,
                           const std::set<int>&     qmIndices,
                           const std::vector<bool>& bQMBlock,
                           QMMMTopologyInfo&        topInfo,
                           WarningHandler*          wi);

/*! \brief Builds link frontier vector with pairs of atoms indicting broken QM - MM chemical bonds.
 *
 * Also performs search of constrained bonds within QM subsystem.
 * Provides data about links via \p topInfo
 * \param[in,out] mtop topology to be modified
 * \param[in] qmIndices set with global indices of QM atoms
 * \param[in] bQMBlock vector with flags for QM-containing blocks
 * \param[in,out] topInfo structure with information about topology modifications
 * \param[in] logger MDLogger for logging info about modifications
 * \returns vector of link atom pairs
 */
std::vector<LinkFrontier> buildQMMMLink(gmx_mtop_t*              mtop,
                                        const std::set<int>&     qmIndices,
                                        const std::vector<bool>& bQMBlock,
                                        QMMMTopologyInfo&        topInfo,
                                        const MDLogger&          logger);

} // namespace gmx

#endif
