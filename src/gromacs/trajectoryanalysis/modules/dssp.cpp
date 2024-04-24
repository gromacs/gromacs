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
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::Dssp.
 *
 * \author Sergey Gorelov <gorelov_sv@pnpi.nrcki.ru>
 * \author Anatoly Titov <titov_ai@pnpi.nrcki.ru>
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "dssp.h"

#include <algorithm>
#include <bitset>
#include <fstream>
#include <iostream>
#include <set>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/units.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

//! Structure that contains storage information from different frames.
struct DsspStorageFrame
{
    //! Frame number.
    int frameNumber_ = 0;
    //! Frame dssp data.
    std::string dsspData_;
};

/*! \brief
 * Class that stores frame information in storage and, upon request, can return it.
 */
class DsspStorage
{
public:
    /*! \brief
     * Function that stores frame information in storage.
     */
    void addData(int frnr, const std::string& data);
    /*! \brief
     * Function that returns frame information from storage.
     */
    const std::vector<DsspStorageFrame>& getData();

private:
    /*! \brief
     * Vector that contains information from different frames.
     */
    std::vector<DsspStorageFrame> data_;
};

void DsspStorage::addData(int frnr, const std::string& data)
{
    DsspStorageFrame dsspData;
    dsspData.frameNumber_ = frnr;
    dsspData.dsspData_    = data;
    data_.push_back(dsspData);
}

const std::vector<DsspStorageFrame>& DsspStorage::getData()
{
    return data_;
}

//! Enum of backbone atoms' types.
enum class BackboneAtomTypes : std::size_t
{
    AtomCA,
    AtomC,
    AtomO,
    AtomN,
    AtomH,
    Count
};

//! String values corresponding to backbone atom types.
const gmx::EnumerationArray<BackboneAtomTypes, const char*> c_backboneAtomTypeNames = {
    { "CA", "C", "O", "N", "H" }
};

/*! \brief
 * Structure of residues' information that can operate with atoms' indices.
 */
struct ResInfo
{
    /*! \brief
     * Size_t array of atoms' indices corresponding to backbone atom types.
     */
    gmx::EnumerationArray<BackboneAtomTypes, size_t> backboneIndices_ = { 0, 0, 0, 0, 0 };
    /*! \brief
     * Bitset of atoms' status. Used to minimize misinterpretation of data.
     */
    std::bitset<static_cast<size_t>(BackboneAtomTypes::Count)> backboneIndicesStatus_{ 0x0 };
    /*! \brief
     * Function that returns atom's index based on specific atom type.
     */
    std::size_t getIndex(BackboneAtomTypes atomTypeName) const;
    /*! \brief
     * Function that returns atom's status based on specific atom type. Used to minimize misinterpretation of data.
     */
    bool hasIndex(BackboneAtomTypes atomTypeName) const;
    /*! \brief
     * Function that sets atom's index and status based on specific atom type.
     */
    void setIndex(BackboneAtomTypes atomTypeName, std::size_t index);
    //! Constant value that defines maximum amount of donors and acceptors per residue.
    static constexpr std::size_t sc_maxDonorsPerResidue = 2;
    //! Pointer to t_resinfo which contains full information about this specific residue.
    t_resinfo* info_ = nullptr;
    //! Pointer to t_resinfo which contains full information about this residue's h-bond donors.
    t_resinfo* donor_[sc_maxDonorsPerResidue] = { nullptr, nullptr };
    //! Pointer to t_resinfo which contains full information about this residue's h-bond acceptors.
    t_resinfo* acceptor_[sc_maxDonorsPerResidue] = { nullptr, nullptr };
    //! Pointer to previous residue in list.
    ResInfo* prevResi_ = nullptr;
    //! Pointer to next residue in list.
    ResInfo* nextResi_ = nullptr;
    //! Float value of h-bond energy with this residue's donors.
    float donorEnergy_[sc_maxDonorsPerResidue] = { 0, 0 };
    //! Float value of h-bond energy with this residue's accpetors.
    float acceptorEnergy_[sc_maxDonorsPerResidue] = { 0, 0 };
    //! Bool value that defines either this residue is proline (PRO) or not.
    bool isProline_ = false;
};

std::size_t ResInfo::getIndex(BackboneAtomTypes atomTypeName) const
{
    return backboneIndices_[atomTypeName];
}

bool ResInfo::hasIndex(BackboneAtomTypes atomTypeName) const
{
    return backboneIndicesStatus_.test(static_cast<std::size_t>(atomTypeName));
}

void ResInfo::setIndex(BackboneAtomTypes atomTypeName, std::size_t index)
{
    backboneIndices_[atomTypeName] = index;
    backboneIndicesStatus_.set(static_cast<std::size_t>(atomTypeName));
}

//! Enum of secondary structures' types.
enum class SecondaryStructureTypes : std::size_t
{
    Loop = 0, //! ~
    Break,    //! =
    Bend,     //! S
    Turn,     //! T
    Helix_PP, //! P
    Helix_5,  //! I
    Helix_3,  //! G
    Strand,   //! E
    Bridge,   //! B
    Helix_4,  //! H
    Count

};

//! String values (symbols) corresponding to secondary structures' types.
const gmx::EnumerationArray<SecondaryStructureTypes, const char> c_secondaryStructureTypeNames = {
    { '~', '=', 'S', 'T', 'P', 'I', 'G', 'E', 'B', 'H' }
};

//! String values (full) corresponding to secondary structures' types.
const gmx::EnumerationArray<SecondaryStructureTypes, const char*> c_secondaryStructureTypeNamesFull = {
    { "Loops",
      "Breaks",
      "Bends",
      "Turns",
      "PP_Helices",
      "π-Helices",
      "3⏨-Helices",
      "β-Strands",
      "β-Bridges",
      "α-Helices" }
};

//! Enum of turns' types.
enum class TurnsTypes : std::size_t
{
    Turn_3 = 0,
    Turn_4,
    Turn_5,
    Turn_PP,
    Count
};

//! Enum of different possible helix positions' types.
enum class HelixPositions : std::size_t
{
    None = 0,
    Start,
    Middle,
    End,
    StartAndEnd,
    Count
};

//! Enum of bridges' types.
enum class BridgeTypes : std::size_t
{
    None = 0,
    AntiParallelBridge,
    ParallelBridge,
    Count
};

/*! \brief
 * Enum of various modes of use of hydrogen atoms. Gromacs mode strictly uses hydrogen atoms from the protein structure,
 * while Dssp mode exclusively uses non-existent hydrogen atoms with coordinates calculated from the positions of carbon and oxygen atoms in residues.
 */
enum class HydrogenMode : std::size_t
{
    Gromacs = 0,
    Dssp,
    Count
};

//! String values corresponding to hydrogen-assignment modes.
const gmx::EnumerationArray<HydrogenMode, const char*> c_HydrogenModeNames = { { "gromacs",
                                                                                 "dssp" } };


/*! \brief
 * Enum of various size of stretch of polyproline helices.
 */
enum class PPStretches : std::size_t
{
    Shortened = 0,
    Default,
    Count
};

//! String values corresponding to neighbor-search modes.
const gmx::EnumerationArray<PPStretches, const char*> c_PPStretchesNames = { { "shortened",
                                                                               "default" } };

/*! \brief
 * Enum of various modes of hydrogen bond definition.
 */
enum class HBondDefinition : std::size_t
{
    Energy = 0,
    Geometry,
    Count
};

//! String values corresponding to hydrogen bond definition modes.
const gmx::EnumerationArray<HBondDefinition, const char*> c_HBondDefinition = { { "energy",
                                                                                  "geometry" } };

/*! \brief
 * Describes and manipulates secondary structure attributes of a residue.
 */
class SecondaryStructuresData
{
public:
    /*! \brief
     * Function that sets status of specific secondary structure to a residue.
     */
    void setSecondaryStructureType(SecondaryStructureTypes secondaryStructureTypeName);
    /*! \brief
     * Function that sets status of specific helix position of specific turns' type to a residue.
     */
    void setHelixPosition(HelixPositions helixPosition, TurnsTypes turn);
    /*! \brief
     * Function that sets status "Break" to a residue and its break partner.
     */
    void setBreak(SecondaryStructuresData* breakPartner);
    /*! \brief
     * Function that sets status "Bridge" or "Anti-Bridge" to a residue and its bridge partner.
     */
    void setBridge(std::size_t bridgePartnerIndex, BridgeTypes bridgeType);
    /*! \brief
     * Function that returns array of residue's bridges indexes.
     */
    const std::vector<std::size_t>& getBridges(BridgeTypes bridgeType);
    /*! \brief
     * Function that returns boolean status of break existence with another specific residue.
     */
    bool isBreakPartnerWith(const SecondaryStructuresData* partner) const;
    /*! \brief
     * Function that returns boolean status of break existence within residue.
     */
    bool hasBreaks() const;
    /*! \brief
     * Function that returns boolean status of bridge existence with another specific residue.
     */
    bool hasBridges(BridgeTypes bridgeType) const;
    /*! \brief
     * Returns which part (None/Start/Middle/End/StartAndEnd) of a helix is represented by \c turn.
     */
    HelixPositions getHelixPosition(TurnsTypes turn) const;
    /*! \brief
     * Function that returns status of specific secondary structure in a residue.
     */
    SecondaryStructureTypes getSecondaryStructure() const;

private:
    //! Constant value that defines maximum amount of breaks (between residues) per residue.
    static const std::size_t sc_maxBreaksPerResidue = 2;
    //! Array of pointers to other residues that forms breaks with this residue.
    SecondaryStructuresData* breakPartners_[sc_maxBreaksPerResidue] = { nullptr, nullptr };
    //! Array of other residues indexes that forms parralel bridges with this residue.
    std::vector<std::size_t> parallelBridgePartners_;
    //! Array of other residues indexes that forms antiparallel bridges with this residue.
    std::vector<std::size_t> antiBridgePartners_;
    //! Secondary structure's status of this residue.
    SecondaryStructureTypes secondaryStructure_ = SecondaryStructureTypes::Loop;
    //! Break status of this residue.
    bool isBreak_ = false;
    //! Helix positions (None/Start/Middle/End/StartAndEnd) corresponding to turns types.
    gmx::EnumerationArray<TurnsTypes, HelixPositions> turnsStatusArray_{ HelixPositions::None,
                                                                         HelixPositions::None,
                                                                         HelixPositions::None,
                                                                         HelixPositions::None };
};

void SecondaryStructuresData::setSecondaryStructureType(const SecondaryStructureTypes secondaryStructureTypeName)
{
    secondaryStructure_ = secondaryStructureTypeName;
}

void SecondaryStructuresData::setHelixPosition(const HelixPositions helixPosition, const TurnsTypes turn)
{
    turnsStatusArray_[turn] = helixPosition;
}

bool SecondaryStructuresData::isBreakPartnerWith(const SecondaryStructuresData* partner) const
{
    return breakPartners_[0] == partner || breakPartners_[1] == partner;
}

HelixPositions SecondaryStructuresData::getHelixPosition(const TurnsTypes turn) const
{
    return turnsStatusArray_[turn];
}

SecondaryStructureTypes SecondaryStructuresData::getSecondaryStructure() const
{
    return secondaryStructure_;
}

void SecondaryStructuresData::setBreak(SecondaryStructuresData* breakPartner)
{
    if (breakPartners_[0] != nullptr)
    {
        breakPartners_[1] = breakPartner;
    }
    else
    {
        breakPartners_[0] = breakPartner;
    }
    isBreak_ = true;
}

void SecondaryStructuresData::setBridge(std::size_t bridgePartnerIndex, BridgeTypes bridgeType)
{
    GMX_RELEASE_ASSERT((bridgeType == BridgeTypes::ParallelBridge
                        || bridgeType == BridgeTypes::AntiParallelBridge),
                       "Unsupported bridge type.");

    if (bridgeType == BridgeTypes::ParallelBridge)
    {
        parallelBridgePartners_.emplace_back(bridgePartnerIndex);
    }
    else
    {
        antiBridgePartners_.emplace_back(bridgePartnerIndex);
    }
}

const std::vector<std::size_t>& SecondaryStructuresData::getBridges(BridgeTypes bridgeType)
{
    GMX_RELEASE_ASSERT((bridgeType == BridgeTypes::ParallelBridge
                        || bridgeType == BridgeTypes::AntiParallelBridge),
                       "Unsupported bridge type.");

    if (bridgeType == BridgeTypes::ParallelBridge)
    {
        return parallelBridgePartners_;
    }
    else
    {
        return antiBridgePartners_;
    }
}

bool SecondaryStructuresData::hasBreaks() const
{
    return isBreak_;
}

bool SecondaryStructuresData::hasBridges(BridgeTypes bridgeType) const
{

    GMX_RELEASE_ASSERT((bridgeType == BridgeTypes::ParallelBridge
                        || bridgeType == BridgeTypes::AntiParallelBridge),
                       "Unsupported bridge type.");

    if (bridgeType == BridgeTypes::ParallelBridge)
    {
        return !parallelBridgePartners_.empty();
    }
    else
    {
        return !antiBridgePartners_.empty();
    }
}

/*! \brief
 * Class that provides search of specific h-bond patterns within residues.
 */
class SecondaryStructures
{
public:
    /*! \brief
     * Function that parses topology to construct vector containing information about the residues.
     */
    void analyseTopology(const TopologyInformation& top,
                         const Selection&           sel,
                         const HydrogenMode&        transferredHMode,
                         bool                       clearStructure);
    /*! \brief
     * Function that checks if topologyVector_ is empty. Used after parsing topology data. If it is
     * empty after running analyseTopology(), then some error has occurred.
     */
    bool topologyIsIncorrect() const;
    /*! \brief
     * Complex function that provides h-bond patterns search and returns string of one-letter secondary structure definitions.
     */
    std::string performPatternSearch(const t_trxframe& fr,
                                     const t_pbc*      pbc,
                                     bool              transferredNbsMode,
                                     real              transferredCutoff,
                                     bool              transferredPiHelicesPreference,
                                     PPStretches       transferredPolyProStretch,
                                     HBondDefinition   transferredHbDef);

private:
    //! Function that parses information from a frame to determine hydrogen bonds (via energy or geometry calculation) patterns.
    void analyzeHydrogenBondsInFrame(const t_trxframe& fr, const t_pbc* pbc, bool nBSmode, real cutoff);
    /*! \brief
     * Function that provides a simple test if a h-bond exists within two residues of specific indices.
     */
    bool hasHBondBetween(std::size_t Donor, std::size_t Acceptor) const;
    /*! \brief
     * Function that provides a simple test if a chain break exists within two residues of specific indices.
     */
    bool noChainBreaksBetween(std::size_t residueA, std::size_t residueB) const;
    /*! \brief
     * Function that calculates if bridge or anti-bridge exists within two residues of specific indices.
     */
    BridgeTypes calculateBridge(std::size_t residueA, std::size_t residueB) const;
    /*! \brief
     * Complex function that provides h-bond patterns search of bridges and strands. Part of patternSearch() complex function.
     */
    void analyzeBridgesAndStrandsPatterns();
    /*! \brief
     * Complex function that provides h-bond patterns search of turns and helices. Part of patternSearch() complex function.
     */
    void analyzeTurnsAndHelicesPatterns();
    /*! \brief
     * Function that calculates atomic distances between atoms A and B based on atom indices.
     */
    static float calculateAtomicDistances(std::size_t       atomA,
                                          std::size_t       atomB,
                                          const t_trxframe& fr,
                                          const t_pbc*      pbc);
    /*! \brief
     * Function that calculates atomic distances between atoms A and B based on atom indices (for atom B) and atom coordinates (for atom A).
     */
    static float calculateAtomicDistances(rvec atomA, std::size_t atomB, const t_trxframe& fr, const t_pbc* pbc);
    /*! \brief
     * Function that calculates Dihedral Angles based on atom indices.
     */
    static float
    calculateDihedralAngle(int atomA, int atomB, int atomC, int atomD, const t_trxframe& fr, const t_pbc* pbc);
    /*! \brief
     * Function that calculates dihedral angles in secondary structure map.
     */
    void calculateDihedrals(const t_trxframe& fr, const t_pbc* pbc);
    /*! \brief
     * Function that calculates bends and breaks in secondary structure map.
     */
    void calculateBends(const t_trxframe& fr, const t_pbc* pbc);

    /*! \brief
     * Function that checks if H-Bond exist according to DSSP algorithm
     * kCouplingConstant = 27.888,  //  = 332 * 0.42 * 0.2
     * E = k * (1/rON + 1/rCH - 1/rOH - 1/rCN) where CO comes from one AA and NH from another
     * if R is in A
     * Hbond exists if E < -0.5
     */
    void calculateHBondEnergy(ResInfo* Donor, ResInfo* Acceptor, const t_trxframe& fr, const t_pbc* pbc);
    /*! \brief
     * Function that checks if H-Bond exist according to HBOND algorithm
     * H-Bond exists if distance between Donor and Acceptor
     * d <= 0.35 nm
     * and Hydrogen-Donor-Acceptor angle
     * α is < 30°.
     */
    void calculateHBondGeometry(ResInfo* Donor, ResInfo* Acceptor, const t_trxframe& fr, const t_pbc* pbc);
    //! Vector that contains h-bond pattern information-manipulating class for each residue in selection.
    std::vector<SecondaryStructuresData> secondaryStructuresStatusVector_;
    //! Vector of ResInfo struct that contains all important information from topology about residues in the protein structure.
    std::vector<ResInfo> topologyVector_;
    //! Vector of ResInfo. Each new frame information from topologyVector_ is copied to frameVector_ to provide frame information independency.
    std::vector<ResInfo> frameVector_;
    //! String that contains result of dssp calculations for output.
    std::string secondaryStructuresStringLine_;
    //! Constant float value of h-bond energy. If h-bond energy within residues is smaller than that value, then h-bond exists.
    const float hBondEnergyCutOff_ = -0.5F;
    //! Constant float value that determines the minimum possible distance (in Å) between two Ca atoms of amino acids of the protein,
    //! exceeding which a hydrogen bond between these two residues will be impossible.
    const float minimalCAdistance_ = 9.0F;
    //! Boolean value that indicates the priority of calculating pi-helices.
    bool piHelicesPreference_ = false;
    //! Enum value for creating hydrogen atoms mode. Very useful for structures without hydrogen atoms. Set in initial options.
    HydrogenMode hMode_ = HydrogenMode::Gromacs;
    //! Enum value that defines polyproline helix stretch. Can be only equal to 2 or 3. Set in initial options.
    PPStretches polyProStretch_ = PPStretches::Default;
    //! Enum value that defines hydrogen bond definition. Set in initial options.
    HBondDefinition hbDef_ = HBondDefinition::Energy;
};

void SecondaryStructures::analyseTopology(const TopologyInformation& top,
                                          const Selection&           sel,
                                          const HydrogenMode&        transferredHMode,
                                          bool                       clearStructure)
{
    hMode_ = transferredHMode;
    int resicompare =
            top.atoms()->atom[static_cast<std::size_t>(*(sel.atomIndices().begin()))].resind - 1;
    for (const auto& ai : sel.atomIndices())
    {
        if (resicompare != top.atoms()->atom[static_cast<std::size_t>(ai)].resind)
        {
            resicompare = top.atoms()->atom[static_cast<std::size_t>(ai)].resind;
            topologyVector_.emplace_back();
            topologyVector_.back().info_ = &(top.atoms()->resinfo[resicompare]);
            std::string residueName      = *(topologyVector_.back().info_->name);
            if (residueName == "PRO")
            {
                topologyVector_.back().isProline_ = true;
            }
        }
        std::string atomName(*(top.atoms()->atomname[static_cast<std::size_t>(ai)]));
        if (atomName == c_backboneAtomTypeNames[BackboneAtomTypes::AtomCA])
        {
            topologyVector_.back().setIndex(BackboneAtomTypes::AtomCA, ai);
        }
        else if (atomName == c_backboneAtomTypeNames[BackboneAtomTypes::AtomC])
        {
            topologyVector_.back().setIndex(BackboneAtomTypes::AtomC, ai);
        }
        else if (atomName == c_backboneAtomTypeNames[BackboneAtomTypes::AtomO])
        {
            topologyVector_.back().setIndex(BackboneAtomTypes::AtomO, ai);
        }
        else if (atomName == c_backboneAtomTypeNames[BackboneAtomTypes::AtomN])
        {
            topologyVector_.back().setIndex(BackboneAtomTypes::AtomN, ai);
            if (hMode_ == HydrogenMode::Dssp)
            {
                topologyVector_.back().setIndex(BackboneAtomTypes::AtomH, ai);
            }
        }
        else if (hMode_ == HydrogenMode::Gromacs
                 && atomName == c_backboneAtomTypeNames[BackboneAtomTypes::AtomH])
        {
            topologyVector_.back().setIndex(BackboneAtomTypes::AtomH, ai);
        }
    }
    if (clearStructure)
    {
        auto isCorrupted = [](const ResInfo& Res) -> bool {
            return !Res.hasIndex(BackboneAtomTypes::AtomCA) || !Res.hasIndex(BackboneAtomTypes::AtomC)
                   || !Res.hasIndex(BackboneAtomTypes::AtomO) || !Res.hasIndex(BackboneAtomTypes::AtomN)
                   || !Res.hasIndex(BackboneAtomTypes::AtomH);
        };
        auto corruptedResidues =
                std::remove_if(topologyVector_.begin(), topologyVector_.end(), isCorrupted);
        topologyVector_.erase(corruptedResidues, topologyVector_.end());
    }
    for (std::size_t i = 1; i < topologyVector_.size(); ++i)
    {
        topologyVector_[i].prevResi_     = &(topologyVector_[i - 1]);
        topologyVector_[i - 1].nextResi_ = &(topologyVector_[i]);
    }
}

bool SecondaryStructures::topologyIsIncorrect() const
{
    return (topologyVector_.empty());
}

void SecondaryStructures::analyzeHydrogenBondsInFrame(const t_trxframe& fr, const t_pbc* pbc, bool nBSmode, real cutoff)
{
    if (nBSmode)
    {
        std::vector<gmx::RVec> positionsCA;
        for (std::size_t i = 0; i < frameVector_.size(); ++i)
        {
            positionsCA.emplace_back(fr.x[frameVector_[i].getIndex(BackboneAtomTypes::AtomCA)]);
        }
        AnalysisNeighborhood nb;
        nb.setCutoff(cutoff);
        AnalysisNeighborhoodPositions       nbPos(positionsCA);
        gmx::AnalysisNeighborhoodSearch     start      = nb.initSearch(pbc, nbPos);
        gmx::AnalysisNeighborhoodPairSearch pairSearch = start.startPairSearch(nbPos);
        gmx::AnalysisNeighborhoodPair       pair;
        ResInfo*                            donor;
        ResInfo*                            acceptor;
        while (pairSearch.findNextPair(&pair))
        {
            if (pair.refIndex() < pair.testIndex())
            {
                donor    = &frameVector_[pair.refIndex()];
                acceptor = &frameVector_[pair.testIndex()];
            }
            else
            {
                continue;
            }
            switch (hbDef_)
            {
                case HBondDefinition::Energy:
                    calculateHBondEnergy(donor, acceptor, fr, pbc);
                    if (acceptor->info_ != donor->nextResi_->info_)
                    {
                        calculateHBondEnergy(acceptor, donor, fr, pbc);
                    }
                    break;
                case HBondDefinition::Geometry:
                    calculateHBondGeometry(donor, acceptor, fr, pbc);
                    if (acceptor->info_ != donor->nextResi_->info_)
                    {
                        calculateHBondGeometry(acceptor, donor, fr, pbc);
                    }
                    break;
                default: continue;
            }
        }
    }
    else
    {
        for (std::size_t donor = 0; donor + 1 < frameVector_.size(); ++donor)
        {
            for (std::size_t acceptor = donor + 1; acceptor < frameVector_.size(); ++acceptor)
            {
                switch (hbDef_)
                {
                    case HBondDefinition::Energy:
                        calculateHBondEnergy(&frameVector_[donor], &frameVector_[acceptor], fr, pbc);
                        if (acceptor != donor + 1)
                        {
                            calculateHBondEnergy(&frameVector_[acceptor], &frameVector_[donor], fr, pbc);
                        }
                        break;
                    case HBondDefinition::Geometry:
                        calculateHBondGeometry(&frameVector_[donor], &frameVector_[acceptor], fr, pbc);
                        if (acceptor != donor + 1)
                        {
                            calculateHBondGeometry(&frameVector_[acceptor], &frameVector_[donor], fr, pbc);
                        }
                        break;
                    default: continue;
                }
            }
        }
    }
}


bool SecondaryStructures::hasHBondBetween(std::size_t donor, std::size_t acceptor) const
{
    for (std::size_t i = 0; i < ResInfo::sc_maxDonorsPerResidue; ++i)
    {
        if (frameVector_[donor].acceptor_[i] == frameVector_[acceptor].info_
            && (frameVector_[donor].acceptorEnergy_[i] < hBondEnergyCutOff_
                || hbDef_ == HBondDefinition::Geometry))
        {
            return true;
        }
    }
    return false;
}

bool SecondaryStructures::noChainBreaksBetween(std::size_t residueA, std::size_t residueB) const
{
    if (residueA > residueB)
    {
        std::swap(residueA, residueB);
    }
    for (; residueA != residueB; ++residueA)
    {
        if (secondaryStructuresStatusVector_[residueA].isBreakPartnerWith(
                    &secondaryStructuresStatusVector_[residueA + 1])
            && secondaryStructuresStatusVector_[residueA + 1].isBreakPartnerWith(
                    &secondaryStructuresStatusVector_[residueA]))
        {
            return false;
        }
    }
    return true;
}

BridgeTypes SecondaryStructures::calculateBridge(std::size_t residueA, std::size_t residueB) const
{
    if (residueA < 1 || residueB < 1 || residueA + 1 >= frameVector_.size()
        || residueB + 1 >= frameVector_.size())
    {
        return BridgeTypes::None;
    }
    if (noChainBreaksBetween(residueA - 1, residueA + 1)
        && noChainBreaksBetween(residueB - 1, residueB + 1) && frameVector_[residueA].prevResi_
        && frameVector_[residueA].nextResi_ && frameVector_[residueB].prevResi_
        && frameVector_[residueB].nextResi_)
    {
        if ((hasHBondBetween(residueA + 1, residueB) && hasHBondBetween(residueB, residueA - 1))
            || (hasHBondBetween(residueB + 1, residueA) && hasHBondBetween(residueA, residueB - 1)))
        {
            return BridgeTypes::ParallelBridge;
        }
        else if ((hasHBondBetween(residueA + 1, residueB - 1) && hasHBondBetween(residueB + 1, residueA - 1))
                 || (hasHBondBetween(residueB, residueA) && hasHBondBetween(residueA, residueB)))
        {
            return BridgeTypes::AntiParallelBridge;
        }
        else
        {
            return BridgeTypes::None;
        }
    }
    return BridgeTypes::None;
}

void SecondaryStructures::analyzeBridgesAndStrandsPatterns()
{
    for (std::size_t i = 1; i + 4 < secondaryStructuresStatusVector_.size(); ++i)
    {
        for (std::size_t j = i + 3; j + 1 < secondaryStructuresStatusVector_.size(); ++j)
        {
            switch (calculateBridge(i, j))
            {
                case BridgeTypes::ParallelBridge:
                {
                    secondaryStructuresStatusVector_[i].setBridge(j, BridgeTypes::ParallelBridge);
                    secondaryStructuresStatusVector_[j].setBridge(i, BridgeTypes::ParallelBridge);
                    break;
                }
                case BridgeTypes::AntiParallelBridge:
                {
                    secondaryStructuresStatusVector_[i].setBridge(j, BridgeTypes::AntiParallelBridge);
                    secondaryStructuresStatusVector_[j].setBridge(i, BridgeTypes::AntiParallelBridge);
                    break;
                }
                default: continue;
            }
        }
    }
    for (std::size_t i = 1; i + 1 < secondaryStructuresStatusVector_.size(); ++i)
    {
        for (std::size_t j = 1; j < 3 and i + j < secondaryStructuresStatusVector_.size(); ++j)
        {
            for (const BridgeTypes& bridgeType :
                 { BridgeTypes::ParallelBridge, BridgeTypes::AntiParallelBridge })
            {
                if (secondaryStructuresStatusVector_[i].hasBridges(bridgeType)
                    && secondaryStructuresStatusVector_[i + j].hasBridges(bridgeType)
                    && (noChainBreaksBetween(i - 1, i + 1) && noChainBreaksBetween(i + j - 1, i + j + 1)))
                {
                    std::vector<std::size_t> iPartners =
                            secondaryStructuresStatusVector_[i].getBridges(bridgeType);
                    std::vector<std::size_t> jPartners =
                            secondaryStructuresStatusVector_[i + j].getBridges(bridgeType);
                    for (const size_t iPartner : iPartners)
                    {
                        for (const size_t jPartner : jPartners)
                        {

                            int delta = abs(static_cast<int>(iPartner) - static_cast<int>(jPartner));
                            if (delta < 6)
                            {
                                int secondStrandStart = iPartner;
                                int secondStrandEnd   = jPartner;
                                if (secondStrandStart > secondStrandEnd)
                                {
                                    std::swap(secondStrandStart, secondStrandEnd);
                                }
                                for (std::size_t k = secondStrandStart;
                                     k <= static_cast<std::size_t>(secondStrandEnd);
                                     ++k)
                                {
                                    secondaryStructuresStatusVector_[k].setSecondaryStructureType(
                                            SecondaryStructureTypes::Strand);
                                }
                                for (std::size_t k = 0; k <= j; ++k)
                                {
                                    secondaryStructuresStatusVector_[i + k].setSecondaryStructureType(
                                            SecondaryStructureTypes::Strand);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    for (std::size_t i = 1; i + 1 < secondaryStructuresStatusVector_.size(); ++i)
    {
        if (!(secondaryStructuresStatusVector_[i].getSecondaryStructure() == SecondaryStructureTypes::Strand)
            && (secondaryStructuresStatusVector_[i].hasBridges(BridgeTypes::ParallelBridge)
                || secondaryStructuresStatusVector_[i].hasBridges(BridgeTypes::AntiParallelBridge)))
        {
            secondaryStructuresStatusVector_[i].setSecondaryStructureType(SecondaryStructureTypes::Bridge);
        }
    }
}

void SecondaryStructures::analyzeTurnsAndHelicesPatterns()
{
    for (const TurnsTypes& i : { TurnsTypes::Turn_3, TurnsTypes::Turn_4, TurnsTypes::Turn_5 })
    {
        std::size_t stride = static_cast<std::size_t>(i) + 3;
        for (std::size_t j = 0; j + stride < secondaryStructuresStatusVector_.size(); ++j)
        {
            if (hasHBondBetween(j + stride, j) && noChainBreaksBetween(j, j + stride))
            {
                secondaryStructuresStatusVector_[j + stride].setHelixPosition(HelixPositions::End, i);

                for (std::size_t k = 1; k < stride; ++k)
                {
                    if (secondaryStructuresStatusVector_[j + k].getHelixPosition(i) == HelixPositions::None)
                    {
                        secondaryStructuresStatusVector_[j + k].setHelixPosition(HelixPositions::Middle, i);
                    }
                }
                if (secondaryStructuresStatusVector_[j].getHelixPosition(i) == HelixPositions::End)
                {
                    secondaryStructuresStatusVector_[j].setHelixPosition(HelixPositions::StartAndEnd, i);
                }
                else
                {
                    secondaryStructuresStatusVector_[j].setHelixPosition(HelixPositions::Start, i);
                }
            }
        }
    }

    for (const TurnsTypes& i : { TurnsTypes::Turn_4, TurnsTypes::Turn_3, TurnsTypes::Turn_5 })
    {
        std::size_t stride = static_cast<std::size_t>(i) + 3;
        for (std::size_t j = 1; j + stride < secondaryStructuresStatusVector_.size(); ++j)
        {
            if ((secondaryStructuresStatusVector_[j - 1].getHelixPosition(i) == HelixPositions::Start
                 || secondaryStructuresStatusVector_[j - 1].getHelixPosition(i) == HelixPositions::StartAndEnd)
                && (secondaryStructuresStatusVector_[j].getHelixPosition(i) == HelixPositions::Start
                    || secondaryStructuresStatusVector_[j].getHelixPosition(i) == HelixPositions::StartAndEnd))
            {
                bool                    empty = true;
                SecondaryStructureTypes helix;
                switch (i)
                {
                    case TurnsTypes::Turn_3:
                        for (std::size_t k = 0; empty && k < stride; ++k)
                        {
                            empty = secondaryStructuresStatusVector_[j + k].getSecondaryStructure()
                                    <= SecondaryStructureTypes::Helix_3;
                        }
                        helix = SecondaryStructureTypes::Helix_3;
                        break;
                    case TurnsTypes::Turn_5:
                        for (std::size_t k = 0; empty && k < stride; ++k)
                        {
                            empty = secondaryStructuresStatusVector_[j + k].getSecondaryStructure()
                                            <= SecondaryStructureTypes::Helix_5
                                    || (piHelicesPreference_
                                        && secondaryStructuresStatusVector_[j + k].getSecondaryStructure()
                                                   == SecondaryStructureTypes::Helix_4);
                        }
                        helix = SecondaryStructureTypes::Helix_5;
                        break;
                    default: helix = SecondaryStructureTypes::Helix_4; break;
                }
                if (empty || helix == SecondaryStructureTypes::Helix_4)
                {
                    for (std::size_t k = 0; k < stride; ++k)
                    {
                        secondaryStructuresStatusVector_[j + k].setSecondaryStructureType(helix);
                    }
                }
            }
        }
    }
    for (std::size_t i = 1; i + 1 < secondaryStructuresStatusVector_.size(); ++i)
    {
        if (secondaryStructuresStatusVector_[i].getSecondaryStructure() <= SecondaryStructureTypes::Turn)
        {
            bool isTurn = false;
            for (const TurnsTypes& j : { TurnsTypes::Turn_3, TurnsTypes::Turn_4, TurnsTypes::Turn_5 })
            {
                std::size_t stride = static_cast<std::size_t>(j) + 3;
                for (std::size_t k = 1; k < stride and !isTurn; ++k)
                {
                    isTurn = (i >= k)
                             && (secondaryStructuresStatusVector_[i - k].getHelixPosition(j)
                                         == HelixPositions::Start
                                 || secondaryStructuresStatusVector_[i - k].getHelixPosition(j)
                                            == HelixPositions::StartAndEnd);
                }
            }
            if (isTurn)
            {
                secondaryStructuresStatusVector_[i].setSecondaryStructureType(SecondaryStructureTypes::Turn);
            }
        }
    }
}

std::string SecondaryStructures::performPatternSearch(const t_trxframe& fr,
                                                      const t_pbc*      pbc,
                                                      bool              transferredNbsMode,
                                                      real              transferredCutoff,
                                                      bool        transferredPiHelicesPreference,
                                                      PPStretches transferredPolyProStretch,
                                                      HBondDefinition transferredHbDef)
{
    GMX_RELEASE_ASSERT(
            !topologyVector_.empty(),
            "Invalid usage of this function. You have to load topology information before. Run "
            "analyseTopology(...) first.");
    frameVector_         = topologyVector_;
    piHelicesPreference_ = transferredPiHelicesPreference;
    polyProStretch_      = transferredPolyProStretch;
    hbDef_               = transferredHbDef;
    secondaryStructuresStatusVector_.resize(0);
    secondaryStructuresStatusVector_.resize(frameVector_.size());
    secondaryStructuresStringLine_.resize(0);
    analyzeHydrogenBondsInFrame(fr, pbc, transferredNbsMode, transferredCutoff);
    secondaryStructuresStringLine_.resize(frameVector_.size(), '~');
    calculateBends(fr, pbc);
    analyzeBridgesAndStrandsPatterns();
    analyzeTurnsAndHelicesPatterns();
    calculateDihedrals(fr, pbc);
    for (auto i = static_cast<std::size_t>(SecondaryStructureTypes::Bend);
         i != static_cast<std::size_t>(SecondaryStructureTypes::Count);
         ++i)
    {
        for (std::size_t j = 0; j < secondaryStructuresStatusVector_.size(); ++j)
        {
            if (secondaryStructuresStatusVector_[j].getSecondaryStructure()
                == static_cast<SecondaryStructureTypes>(i))
            {
                secondaryStructuresStringLine_[j] = c_secondaryStructureTypeNames[i];
            }
        }
    }
    if (secondaryStructuresStatusVector_.size() > 1)
    {
        for (std::size_t i = 0, lineFactor = 1; i + 1 < secondaryStructuresStatusVector_.size(); ++i)
        {
            if (secondaryStructuresStatusVector_[i].hasBreaks()
                && secondaryStructuresStatusVector_[i + 1].hasBreaks())
            {
                if (secondaryStructuresStatusVector_[i].isBreakPartnerWith(
                            &secondaryStructuresStatusVector_[i + 1])
                    && secondaryStructuresStatusVector_[i + 1].isBreakPartnerWith(
                            &secondaryStructuresStatusVector_[i]))
                {
                    secondaryStructuresStringLine_.insert(
                            secondaryStructuresStringLine_.begin() + i + lineFactor,
                            c_secondaryStructureTypeNames[SecondaryStructureTypes::Break]);
                    ++lineFactor;
                }
            }
        }
    }
    return secondaryStructuresStringLine_;
}

float SecondaryStructures::calculateAtomicDistances(std::size_t       atomA,
                                                    std::size_t       atomB,
                                                    const t_trxframe& fr,
                                                    const t_pbc*      pbc)
{
    gmx::RVec vectorBA = { 0, 0, 0 };
    pbc_dx(pbc, fr.x[atomA], fr.x[atomB], vectorBA.as_vec());
    return vectorBA.norm() * gmx::c_nm2A;
}

float SecondaryStructures::calculateAtomicDistances(rvec              atomA,
                                                    std::size_t       atomB,
                                                    const t_trxframe& fr,
                                                    const t_pbc*      pbc)
{
    gmx::RVec vectorBA = { 0, 0, 0 };
    pbc_dx(pbc, atomA, fr.x[atomB], vectorBA.as_vec());
    return vectorBA.norm() * gmx::c_nm2A;
}

float SecondaryStructures::calculateDihedralAngle(int               atomA,
                                                  int               atomB,
                                                  int               atomC,
                                                  int               atomD,
                                                  const t_trxframe& fr,
                                                  const t_pbc*      pbc)
{
    float     result               = 360;
    float     vdot1                = 0;
    float     vdot2                = 0;
    gmx::RVec vectorBA             = { 0, 0, 0 };
    gmx::RVec vectorCD             = { 0, 0, 0 };
    gmx::RVec vectorCB             = { 0, 0, 0 };
    gmx::RVec vectorCBxBA          = { 0, 0, 0 };
    gmx::RVec vectorCBxCD          = { 0, 0, 0 };
    gmx::RVec vectorCBxvectorCBxCD = { 0, 0, 0 };
    pbc_dx(pbc, fr.x[atomA], fr.x[atomB], vectorBA.as_vec());
    pbc_dx(pbc, fr.x[atomD], fr.x[atomC], vectorCD.as_vec());
    pbc_dx(pbc, fr.x[atomB], fr.x[atomC], vectorCB.as_vec());
    vectorBA *= gmx::c_nm2A;
    vectorCD *= gmx::c_nm2A;
    vectorCB *= gmx::c_nm2A;
    vectorCBxBA          = vectorCB.cross(vectorBA);
    vectorCBxCD          = vectorCB.cross(vectorCD);
    vectorCBxvectorCBxCD = vectorCB.cross(vectorCBxCD);
    vdot1                = vectorCBxCD.dot(vectorCBxCD);
    vdot2                = vectorCBxvectorCBxCD.dot(vectorCBxvectorCBxCD);
    if (vdot1 > 0 and vdot2 > 0)
    {
        vdot1 = vectorCBxBA.dot(vectorCBxCD) / std::sqrt(vdot1);
        vdot2 = vectorCBxBA.dot(vectorCBxvectorCBxCD) / std::sqrt(vdot2);
        if (vdot1 != 0 or vdot2 != 0)
        {
            result = std::atan2(vdot2, vdot1) * gmx::c_rad2Deg;
        }
    }
    return result;
}

void SecondaryStructures::calculateDihedrals(const t_trxframe& fr, const t_pbc* pbc)
{
    // Values are taken from original DSSP algorithm, file Secondary.cpp from https://github.com/PDB-REDO/libcifpp/releases/tag/v3.0.0
    const float        epsilon = 29;
    const float        phiMin  = -75 - epsilon;
    const float        phiMax  = -75 + epsilon;
    const float        psiMin  = 145 - epsilon;
    const float        psiMax  = 145 + epsilon;
    std::vector<float> phi(frameVector_.size(), 360);
    std::vector<float> psi(frameVector_.size(), 360);
    for (std::size_t i = 1; i + 1 < frameVector_.size(); ++i)
    {
        if (frameVector_[i - 1].hasIndex(BackboneAtomTypes::AtomC)
            && frameVector_[i].hasIndex(BackboneAtomTypes::AtomN)
            && frameVector_[i].hasIndex(BackboneAtomTypes::AtomCA)
            && frameVector_[i].hasIndex(BackboneAtomTypes::AtomC))
        {
            phi[i] = calculateDihedralAngle(frameVector_[i - 1].getIndex(BackboneAtomTypes::AtomC),
                                            frameVector_[i].getIndex(BackboneAtomTypes::AtomN),
                                            frameVector_[i].getIndex(BackboneAtomTypes::AtomCA),
                                            frameVector_[i].getIndex(BackboneAtomTypes::AtomC),
                                            fr,
                                            pbc);
        }
        if (frameVector_[i].hasIndex(BackboneAtomTypes::AtomN)
            && frameVector_[i].hasIndex(BackboneAtomTypes::AtomCA)
            && frameVector_[i].hasIndex(BackboneAtomTypes::AtomC)
            && frameVector_[i + 1].hasIndex(BackboneAtomTypes::AtomN))
        {
            psi[i] = calculateDihedralAngle(frameVector_[i].getIndex(BackboneAtomTypes::AtomN),
                                            frameVector_[i].getIndex(BackboneAtomTypes::AtomCA),
                                            frameVector_[i].getIndex(BackboneAtomTypes::AtomC),
                                            frameVector_[i + 1].getIndex(BackboneAtomTypes::AtomN),
                                            fr,
                                            pbc);
        }
    }
    for (std::size_t i = 1; i + 3 < frameVector_.size(); ++i)
    {
        switch (polyProStretch_)
        {
            case PPStretches::Shortened:
            {
                if (phiMin > phi[i] or phi[i] > phiMax or phiMin > phi[i + 1] or phi[i + 1] > phiMax)
                {
                    continue;
                }

                if (psiMin > psi[i] or psi[i] > psiMax or psiMin > psi[i + 1] or psi[i + 1] > psiMax)
                {
                    continue;
                }

                switch (secondaryStructuresStatusVector_[i].getHelixPosition(TurnsTypes::Turn_PP))
                {
                    case HelixPositions::None:
                        secondaryStructuresStatusVector_[i].setHelixPosition(HelixPositions::Start,
                                                                             TurnsTypes::Turn_PP);
                        break;

                    case HelixPositions::End:
                        secondaryStructuresStatusVector_[i].setHelixPosition(
                                HelixPositions::StartAndEnd, TurnsTypes::Turn_PP);
                        break;

                    default: break;
                }
                secondaryStructuresStatusVector_[i + 1].setHelixPosition(HelixPositions::End,
                                                                         TurnsTypes::Turn_PP);
                if (secondaryStructuresStatusVector_[i].getSecondaryStructure()
                    == SecondaryStructureTypes::Loop)
                {
                    secondaryStructuresStatusVector_[i].setSecondaryStructureType(
                            SecondaryStructureTypes::Helix_PP);
                }
                if (secondaryStructuresStatusVector_[i + 1].getSecondaryStructure()
                    == SecondaryStructureTypes::Loop)
                {
                    secondaryStructuresStatusVector_[i + 1].setSecondaryStructureType(
                            SecondaryStructureTypes::Helix_PP);
                }
                break;
            }
            case PPStretches::Default:
            {
                if (phiMin > phi[i] or phi[i] > phiMax or phiMin > phi[i + 1] or phi[i + 1] > phiMax
                    or phiMin > phi[i + 2] or phi[i + 2] > phiMax)
                {
                    continue;
                }

                if (psiMin > psi[i] or psi[i] > psiMax or psiMin > psi[i + 1] or psi[i + 1] > psiMax
                    or psiMin > psi[i + 2] or psi[i + 2] > psiMax)
                {
                    continue;
                }
                switch (secondaryStructuresStatusVector_[i].getHelixPosition(TurnsTypes::Turn_PP))
                {
                    case HelixPositions::None:
                        secondaryStructuresStatusVector_[i].setHelixPosition(HelixPositions::Start,
                                                                             TurnsTypes::Turn_PP);
                        break;

                    case HelixPositions::End:
                        secondaryStructuresStatusVector_[i].setHelixPosition(
                                HelixPositions::StartAndEnd, TurnsTypes::Turn_PP);
                        break;

                    default: break;
                }
                secondaryStructuresStatusVector_[i + 1].setHelixPosition(HelixPositions::Middle,
                                                                         TurnsTypes::Turn_PP);
                secondaryStructuresStatusVector_[i + 2].setHelixPosition(HelixPositions::End,
                                                                         TurnsTypes::Turn_PP);
                if (secondaryStructuresStatusVector_[i].getSecondaryStructure()
                    == SecondaryStructureTypes::Loop)
                {
                    secondaryStructuresStatusVector_[i].setSecondaryStructureType(
                            SecondaryStructureTypes::Helix_PP);
                }
                if (secondaryStructuresStatusVector_[i + 1].getSecondaryStructure()
                    == SecondaryStructureTypes::Loop)
                {
                    secondaryStructuresStatusVector_[i + 1].setSecondaryStructureType(
                            SecondaryStructureTypes::Helix_PP);
                }
                if (secondaryStructuresStatusVector_[i + 2].getSecondaryStructure()
                    == SecondaryStructureTypes::Loop)
                {
                    secondaryStructuresStatusVector_[i + 2].setSecondaryStructureType(
                            SecondaryStructureTypes::Helix_PP);
                }
                break;
            }
            default: GMX_RELEASE_ASSERT(false, "Unsupported stretch length.");
        }
    }
}

void SecondaryStructures::calculateBends(const t_trxframe& fr, const t_pbc* pbc)
{
    // Values are taken from original DSSP algorithm, file Secondary.cpp from https://github.com/PDB-REDO/libcifpp/releases/tag/v3.0.0
    const float bendDegreeMin = 70.0;
    const float bendDegreeMax = 360.0;
    const float maxDist       = 2.5; // note, in Angstrom
    float       degree        = 0;
    gmx::RVec   vecAB{ 0, 0, 0 };
    gmx::RVec   vecAC{ 0, 0, 0 };
    for (std::size_t i = 0; i + 1 < frameVector_.size(); ++i)
    {
        if (frameVector_[i].hasIndex(BackboneAtomTypes::AtomC)
            && frameVector_[i + 1].hasIndex(BackboneAtomTypes::AtomN))
        {
            if (calculateAtomicDistances(frameVector_[i].getIndex(BackboneAtomTypes::AtomC),
                                         frameVector_[i + 1].getIndex(BackboneAtomTypes::AtomN),
                                         fr,
                                         pbc)
                > maxDist)
            {
                secondaryStructuresStatusVector_[i].setBreak(&secondaryStructuresStatusVector_[i + 1]);
                secondaryStructuresStatusVector_[i + 1].setBreak(&secondaryStructuresStatusVector_[i]);
            }
        }
        else
        {
            secondaryStructuresStatusVector_[i].setBreak(&secondaryStructuresStatusVector_[i + 1]);
            secondaryStructuresStatusVector_[i + 1].setBreak(&secondaryStructuresStatusVector_[i]);
        }
    }
    for (std::size_t i = 2; i + 2 < frameVector_.size(); ++i)
    {
        if (secondaryStructuresStatusVector_[i - 2].isBreakPartnerWith(
                    &(secondaryStructuresStatusVector_[i - 1]))
            || secondaryStructuresStatusVector_[i - 1].isBreakPartnerWith(
                    &(secondaryStructuresStatusVector_[i]))
            || secondaryStructuresStatusVector_[i].isBreakPartnerWith(
                    &(secondaryStructuresStatusVector_[i + 1]))
            || secondaryStructuresStatusVector_[i + 1].isBreakPartnerWith(
                    &(secondaryStructuresStatusVector_[i + 2])))
        {
            continue;
        }
        pbc_dx(pbc,
               fr.x[frameVector_[i].getIndex(BackboneAtomTypes::AtomCA)],
               fr.x[frameVector_[i - 2].getIndex(BackboneAtomTypes::AtomCA)],
               vecAB.as_vec());
        pbc_dx(pbc,
               fr.x[frameVector_[i + 2].getIndex(BackboneAtomTypes::AtomCA)],
               fr.x[frameVector_[i].getIndex(BackboneAtomTypes::AtomCA)],
               vecAC.as_vec());
        degree = gmx_angle(vecAB, vecAC) * gmx::c_rad2Deg;
        if (degree != bendDegreeMax and degree > bendDegreeMin)
        {
            secondaryStructuresStatusVector_[i].setSecondaryStructureType(SecondaryStructureTypes::Bend);
        }
    }
}

void SecondaryStructures::calculateHBondEnergy(ResInfo*          donor,
                                               ResInfo*          acceptor,
                                               const t_trxframe& fr,
                                               const t_pbc*      pbc)
{
    if (!(donor->isProline_)
        && (acceptor->hasIndex(BackboneAtomTypes::AtomC) && acceptor->hasIndex(BackboneAtomTypes::AtomO)
            && donor->hasIndex(BackboneAtomTypes::AtomN) && donor->hasIndex(BackboneAtomTypes::AtomH)))
    {
        if (calculateAtomicDistances(donor->getIndex(BackboneAtomTypes::AtomCA),
                                     acceptor->getIndex(BackboneAtomTypes::AtomCA),
                                     fr,
                                     pbc)
            < minimalCAdistance_)
        {
            float distanceHO = 0;
            float distanceHC = 0;
            float distanceNO = calculateAtomicDistances(donor->getIndex(BackboneAtomTypes::AtomN),
                                                        acceptor->getIndex(BackboneAtomTypes::AtomO),
                                                        fr,
                                                        pbc);
            float distanceNC = calculateAtomicDistances(donor->getIndex(BackboneAtomTypes::AtomN),
                                                        acceptor->getIndex(BackboneAtomTypes::AtomC),
                                                        fr,
                                                        pbc);
            if (hMode_ == HydrogenMode::Dssp)
            {
                if (donor->prevResi_ != nullptr && donor->prevResi_->getIndex(BackboneAtomTypes::AtomC)
                    && donor->prevResi_->getIndex(BackboneAtomTypes::AtomO))
                {
                    gmx::RVec atomH  = fr.x[donor->getIndex(BackboneAtomTypes::AtomH)];
                    gmx::RVec prevCO = fr.x[donor->prevResi_->getIndex(BackboneAtomTypes::AtomC)];
                    prevCO -= fr.x[donor->prevResi_->getIndex(BackboneAtomTypes::AtomO)];
                    float prevCODist = calculateAtomicDistances(
                            donor->prevResi_->getIndex(BackboneAtomTypes::AtomC),
                            donor->prevResi_->getIndex(BackboneAtomTypes::AtomO),
                            fr,
                            pbc);
                    atomH += prevCO / prevCODist;
                    distanceHO = calculateAtomicDistances(
                            atomH, acceptor->getIndex(BackboneAtomTypes::AtomO), fr, pbc);
                    distanceHC = calculateAtomicDistances(
                            atomH, acceptor->getIndex(BackboneAtomTypes::AtomC), fr, pbc);
                }
                else
                {
                    distanceHO = distanceNO;
                    distanceHC = distanceNC;
                }
            }
            else
            {
                distanceHO = calculateAtomicDistances(donor->getIndex(BackboneAtomTypes::AtomH),
                                                      acceptor->getIndex(BackboneAtomTypes::AtomO),
                                                      fr,
                                                      pbc);
                distanceHC = calculateAtomicDistances(donor->getIndex(BackboneAtomTypes::AtomH),
                                                      acceptor->getIndex(BackboneAtomTypes::AtomC),
                                                      fr,
                                                      pbc);
            }
            // Values are taken from original DSSP algorithm, file Secondary.cpp from https://github.com/PDB-REDO/libcifpp/releases/tag/v3.0.0
            float       HbondEnergy         = 0;
            const float minEnergy           = -9.9;
            const float minimalAtomDistance = 0.5;
            const float kCouplingConstant   = 27.888;
            if ((distanceNO < minimalAtomDistance) || (distanceHC < minimalAtomDistance)
                || (distanceHO < minimalAtomDistance) || (distanceNC < minimalAtomDistance))
            {
                HbondEnergy = minEnergy;
            }
            else
            {
                HbondEnergy =
                        kCouplingConstant
                        * ((1 / distanceNO) + (1 / distanceHC) - (1 / distanceHO) - (1 / distanceNC));
            }
            if (HbondEnergy < donor->acceptorEnergy_[0])
            {
                donor->acceptor_[1]       = donor->acceptor_[0];
                donor->acceptorEnergy_[1] = donor->acceptorEnergy_[0];
                donor->acceptor_[0]       = acceptor->info_;
                donor->acceptorEnergy_[0] = HbondEnergy;
            }
            else if (HbondEnergy < donor->acceptorEnergy_[1])
            {
                donor->acceptor_[1]       = acceptor->info_;
                donor->acceptorEnergy_[1] = HbondEnergy;
            }

            if (HbondEnergy < acceptor->donorEnergy_[0])
            {
                acceptor->donor_[1]       = acceptor->donor_[0];
                acceptor->donorEnergy_[1] = acceptor->donorEnergy_[0];
                acceptor->donor_[0]       = donor->info_;
                acceptor->donorEnergy_[0] = HbondEnergy;
            }
            else if (HbondEnergy < acceptor->donorEnergy_[1])
            {
                acceptor->donor_[1]       = donor->info_;
                acceptor->donorEnergy_[1] = HbondEnergy;
            }
        }
    }
}

void SecondaryStructures::calculateHBondGeometry(ResInfo*          donor,
                                                 ResInfo*          acceptor,
                                                 const t_trxframe& fr,
                                                 const t_pbc*      pbc)
{
    if (!(donor->isProline_)
        && (acceptor->hasIndex(BackboneAtomTypes::AtomC) && acceptor->hasIndex(BackboneAtomTypes::AtomO)
            && donor->hasIndex(BackboneAtomTypes::AtomN) && donor->hasIndex(BackboneAtomTypes::AtomH)))
    {
        gmx::RVec vectorNO = { 0, 0, 0 };
        pbc_dx(pbc,
               fr.x[acceptor->getIndex(BackboneAtomTypes::AtomO)],
               fr.x[donor->getIndex(BackboneAtomTypes::AtomN)],
               vectorNO.as_vec());
        // Value is taken from the HBOND algorithm.
        const float c_rMaxDistanceNM_ = 0.35;
        if (vectorNO.norm() <= c_rMaxDistanceNM_)
        {
            gmx::RVec vectorH = fr.x[donor->getIndex(BackboneAtomTypes::AtomH)];
            if (hMode_ == HydrogenMode::Dssp)
            {
                if (donor->prevResi_ != nullptr && donor->prevResi_->getIndex(BackboneAtomTypes::AtomC)
                    && donor->prevResi_->getIndex(BackboneAtomTypes::AtomO))
                {

                    gmx::RVec prevCO = fr.x[donor->prevResi_->getIndex(BackboneAtomTypes::AtomC)];
                    prevCO -= fr.x[donor->prevResi_->getIndex(BackboneAtomTypes::AtomO)];
                    float prevCODist = calculateAtomicDistances(
                            donor->prevResi_->getIndex(BackboneAtomTypes::AtomC),
                            donor->prevResi_->getIndex(BackboneAtomTypes::AtomO),
                            fr,
                            pbc);
                    vectorH += prevCO / prevCODist;
                }
            }
            gmx::RVec vectorNH = { 0, 0, 0 };
            pbc_dx(pbc, vectorH, fr.x[donor->getIndex(BackboneAtomTypes::AtomN)], vectorNH.as_vec());
            // Values are taken from the HBOND algorithm.
            float       degree            = 0;
            const float c_angleMaxDegree_ = 30;
            degree                        = gmx_angle(vectorNO, vectorNH) * gmx::c_rad2Deg;
            if (degree <= c_angleMaxDegree_)
            {
                if (donor->acceptor_[0] == nullptr)
                {
                    donor->acceptor_[0] = acceptor->info_;
                }
                else if (donor->acceptor_[1] == nullptr)
                {
                    donor->acceptor_[1] = donor->acceptor_[0];
                    donor->acceptor_[0] = acceptor->info_;
                }

                if (acceptor->donor_[0] == nullptr)
                {
                    acceptor->donor_[0] = donor->info_;
                }
                else if (acceptor->donor_[1] == nullptr)
                {
                    acceptor->donor_[1] = acceptor->donor_[0];
                    acceptor->donor_[0] = donor->info_;
                }
            }
        }
    }
}

class Dssp : public TrajectoryAnalysisModule
{
public:
    Dssp();
    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void initAfterFirstFrame(const TrajectoryAnalysisSettings& settings, const t_trxframe& fr) override;
    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;
    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    //! Selections for DSSP output. Set in initial options.
    Selection sel_;
    //! Boolean value for Preferring P-Helices mode. Set in initial options.
    bool piHelicesPreference_ = false;
    //! Enum value for creating hydrogen atoms mode. Very useful for structures without hydrogen atoms. Set in initial options.
    HydrogenMode hMode_ = HydrogenMode::Gromacs;
    //! Boolean value determines different calculation methods for searching neighbor residues. Set in initial options.
    bool nBSmode_ = true;
    //! Boolean value determines the removal of defective residues from the structure. Set in initial options.
    bool clearStructure_ = false;
    //! Real value that defines maximum distance from residue to its neighbor residue.
    real cutoff_ = 0.9;
    //! Enum value that defines polyproline helix stretch. Set in initial options.
    PPStretches polyProStretch_ = PPStretches::Default;
    //! Enum value that defines hydrogen bond definition. Set in initial options.
    HBondDefinition hbDef_ = HBondDefinition::Energy;
    //! String value that defines output filename. Set in initial options.
    std::string fnmDSSPOut_ = "dssp.dat";
    //! String value that defines plot output filename. Set in initial options.
    std::string fnmPlotOut_;
    //! Class that calculates h-bond patterns in secondary structure map based on original DSSP algorithm.
    SecondaryStructures patternSearch_;
    //! A storage that contains DSSP info_ from different frames.
    DsspStorage storage_;
    //! Data container for raw data of number of secondary structures per frame.
    AnalysisData ssNumPerFrame_;
};

Dssp::Dssp()
{
    registerAnalysisDataset(&ssNumPerFrame_, "secondaryStructuresNum");
}

void Dssp::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] allows using the DSSP algorithm (namely, by detecting specific patterns of "
        "hydrogen bonds between amino acid residues) "
        "to determine the secondary structure of a protein.[PAR]"
        "One-symbol secondary structure designations that are used in the output file:[PAR]"
        "[TT]H[tt] — [GRK]alpha[grk]-helix;[PAR]"
        "[TT]B[tt] — residue in isolated [GRK]beta[grk]-bridge;[PAR]"
        "[TT]E[tt] — extended strand that participates in [GRK]beta[grk]-ladder;[PAR]"
        "[TT]G[tt] — 3[SUB]10[sub]-helix;[PAR]"
        "[TT]I[tt] — [GRK]pi[grk]-helix;[PAR]"
        "[TT]P[tt] — [GRK]kappa[grk]-helix (poly-proline II helix);[PAR]"
        "[TT]S[tt] — bend;[PAR]"
        "[TT]T[tt] — hydrogen-bonded turn;[PAR]"
        "[TT]=[tt] — break;[PAR]"
        "[TT]~[tt] — loop (no special secondary structure designation).[PAR]"
        "[TT]-num[tt] allows you to get a plot of the number of secondary structures of each type "
        "as a function of time at the output.[PAR]"
        "[TT]-hmode[tt] selects between using hydrogen atoms directly from the structure "
        "(\"gromacs\" option) and using hydrogen pseudo-atoms based on C and O atom coordinates of "
        "previous residue (\"dssp\" option). You should always use the \"dssp\" option for "
        "structures "
        "with absent hydrogen atoms![PAR]"
        "[TT]-hbond[tt] selects between different definitions of hydrogen bond. \"energy\" means "
        "the calculation of a hydrogen bond using the electrostatic interaction energy and "
        "\"geometry\" means the calculation of the hydrogen bond using geometric criterion "
        "for the existence of a hydrogen bond.[PAR]"
        "[TT]-nb[tt] allows using GROMACS neighbor-search method to find residue pairs that may "
        "have a hydrogen bond instead of simply iterating over the residues among themselves.[PAR]"
        "[TT]-cutoff[tt] is a real value that defines maximum distance from residue to its "
        "neighbor residue used in [TT]-nb[tt]. Minimum (and also recommended) value is 0.9.[PAR]"
        "[TT]-clear[tt] allows you to ignore the analysis of the secondary structure residues "
        "that are missing one or more critical atoms (CA, C, N, O or H). Always use this option "
        "together with [TT]-hmode dssp[tt] for structures that lack hydrogen atoms![PAR]"
        "[TT]-pihelix[tt] changes pattern-search algorithm towards preference of pi-helices.[PAR]"
        "[TT]-ppstretch[tt] defines stretch value of polyproline-helices. \"shortened\" means "
        "stretch with size 2 and \"default\" means stretch with size 3.[PAR]"
        "Note that [THISMODULE] currently is not capable of reproducing "
        "the secondary structure of proteins whose structure is determined by methods other than "
        "X-ray crystallography (structures in .pdb format with "
        "incorrect values in the CRYST1 line) due to the incorrect cell size in such "
        "structures.[PAR]"
        "Please note that the computation is always done in single precision, regardless of the "
        "precision for which GROMACS was configured."
    };
    options->addOption(FileNameOption("o")
                               .outputFile()
                               .store(&fnmDSSPOut_)
                               .required()
                               .defaultBasename("dssp")
                               .filetype(OptionFileType::GenericData)
                               .description("Filename for DSSP output"));
    options->addOption(FileNameOption("num")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnmPlotOut_)
                               .defaultBasename("num")
                               .description("Output file name for secondary structures statistics "
                                            "for the trajectory"));
    options->addOption(SelectionOption("sel").store(&sel_).defaultSelectionText("Protein").description(
            "Group for DSSP"));
    options->addOption(EnumOption<HydrogenMode>("hmode")
                               .store(&hMode_)
                               .defaultValue(HydrogenMode::Gromacs)
                               .enumValue(c_HydrogenModeNames)
                               .description("Hydrogens pseudoatoms creating mode"));
    options->addOption(
            EnumOption<HBondDefinition>("hbond")
                    .store(&hbDef_)
                    .defaultValue(HBondDefinition::Energy)
                    .enumValue(c_HBondDefinition)
                    .description("Selects between different definitions of hydrogen bond"));
    options->addOption(BooleanOption("nb").store(&nBSmode_).defaultValue(true).description(
            "Use GROMACS neighbor-search method"));
    options->addOption(RealOption("cutoff").store(&cutoff_).required().defaultValue(0.9).description(
            "Distance from residue to its neighbor residue in neighbor search. Must be >= 0.9"));
    options->addOption(BooleanOption("clear")
                               .store(&clearStructure_)
                               .defaultValue(false)
                               .description("Clear defective residues from the structure"));
    options->addOption(
            BooleanOption("pihelix").store(&piHelicesPreference_).defaultValue(false).description("Prefer Pi Helices"));
    options->addOption(EnumOption<PPStretches>("ppstretch")
                               .store(&polyProStretch_)
                               .defaultValue(PPStretches::Default)
                               .enumValue(c_PPStretchesNames)
                               .description("Stretch value for PP-helices"));
    settings->setHelpText(desc);
}

void Dssp::optionsFinished(TrajectoryAnalysisSettings* /* settings */)
{
    if (cutoff_ < real(0.9))
    {
        GMX_THROW(InconsistentInputError("Invalid cutoff value. It must be >= 0.9."));
    }
}

void Dssp::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top)
{
    patternSearch_.analyseTopology(top, sel_, hMode_, clearStructure_);

    if (patternSearch_.topologyIsIncorrect())
    {
        std::string errorDesc =
                "From these inputs, it is not possible to obtain proper information about the "
                "patterns of hydrogen bonds.";
        if (hMode_ != HydrogenMode::Dssp)
        {
            errorDesc += " Maybe you should add the \"-hmode dssp\" option?";
        }
        GMX_THROW(InconsistentInputError(errorDesc));
    }

    if (!fnmPlotOut_.empty())
    {
        ssNumPerFrame_.setColumnCount(0, static_cast<std::size_t>(SecondaryStructureTypes::Count));
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnmPlotOut_);
        plotm->setTitle("Number of Secondary Structures");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Secondary Structures");
        for (const auto& i : c_secondaryStructureTypeNamesFull)
        {
            plotm->appendLegend(i);
        }
        plotm->setYFormat(10, 0);
        ssNumPerFrame_.addModule(plotm);
    }
}

void Dssp::initAfterFirstFrame(const TrajectoryAnalysisSettings& /* settings */, const t_trxframe& /* fr */)
{
}

void Dssp::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata)
{
    AnalysisDataHandle dhNum_ = pdata->dataHandle(ssNumPerFrame_);
    storage_.addData(frnr,
                     patternSearch_.performPatternSearch(
                             fr, pbc, nBSmode_, cutoff_, piHelicesPreference_, polyProStretch_, hbDef_));
    if (!fnmPlotOut_.empty())
    {
        dhNum_.startFrame(frnr, fr.time);
        const std::string& temp = storage_.getData().back().dsspData_;
        for (std::size_t i = 0; i < static_cast<std::size_t>(SecondaryStructureTypes::Count); ++i)
        {
            dhNum_.setPoint(i, std::count(temp.begin(), temp.end(), c_secondaryStructureTypeNames[i]));
        }
        dhNum_.finishFrame();
    }
}

void Dssp::finishAnalysis(int /*nframes*/)
{
    please_cite(stdout, "Kabsch1983");
    please_cite(stdout, "Gorelov2024");
}

void Dssp::writeOutput()
{
    std::vector<DsspStorageFrame> dataOut;
    FILE*                         fp;
    fp      = gmx_ffopen(fnmDSSPOut_, "w");
    dataOut = storage_.getData();
    for (auto& i : dataOut)
    {
        std::fprintf(fp, "%s\n", i.dsspData_.c_str());
    }
    gmx_ffclose(fp);
}

} // namespace

const char DsspInfo::name[] = "dssp";
const char DsspInfo::shortDescription[] =
        "Calculate protein secondary structure via DSSP algorithm";

TrajectoryAnalysisModulePointer DsspInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Dssp);
}

} // namespace analysismodules

} // namespace gmx
