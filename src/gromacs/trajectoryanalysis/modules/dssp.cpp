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
#include <set>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

/********************************************************************
 * Dssp
 */

//! BackBone atom types
enum class BackboneAtomTypes : size_t
{
    AtomCA = 0,
    AtomC  = 1,
    AtomO  = 2,
    AtomN  = 3,
    AtomH  = 4,
    Count
};
//! String values corresponding to backbone atom types
const gmx::EnumerationArray<BackboneAtomTypes, const char*> BackboneAtomTypeNames = {
    { "CA", "C", "O", "N", "H" }
};

//! Class for backbone Atoms
class BackboneAtomIndexes
{
public:
    void   addAtomIndex(BackboneAtomTypes atomTypeName, size_t atomIndex);
    size_t getIndex(BackboneAtomTypes atomTypeName) const;

private:
    std::array<size_t, 5> backBoneAtomIndexes_{ 0, 0, 0, 0, 0 };
};


void BackboneAtomIndexes::addAtomIndex(BackboneAtomTypes atomTypeName, size_t atomIndex)
{
    backBoneAtomIndexes_.at(static_cast<size_t>(atomTypeName)) = atomIndex;
}
size_t BackboneAtomIndexes::getIndex(BackboneAtomTypes atomTypeName) const
{
    return backBoneAtomIndexes_[static_cast<size_t>(atomTypeName)];
}

struct DsspStorageFrame
{
    int         frnr = 0;
    std::string dsspData;
};

class DsspStorage
{
public:
    void                          addData(int frnr, const std::string& data);
    std::vector<DsspStorageFrame> getData();

private:
    std::vector<DsspStorageFrame> data_;
};

void DsspStorage::addData(int frnr, const std::string& data)
{
    DsspStorageFrame dsspData;
    dsspData.frnr     = frnr;
    dsspData.dsspData = data;
    data_.push_back(dsspData);
}

std::vector<DsspStorageFrame> DsspStorage::getData()
{
    return data_;
}


class Dssp : public TrajectoryAnalysisModule
{
public:
    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;
    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    real                             cutoff_ = 1.0;
    Selection                        sel_;
    AtomsDataPtr                     atoms_;
    std::string                      fnmDSSPOut_;
    DsspStorage                      DsspStorage_;
    std::vector<std::size_t>         AtomResi_, SecondaryStructuresMap_;
    std::vector<BackboneAtomIndexes> IndexMap_;
    std::vector<std::vector<bool>>   HBondsMap_;
    std::vector<std::size_t>         bendmap_, breakmap_;
    /*! \brief
     * Function that calculate Bends and Breakes in Secondary Structure Map
     */
    void calculateBends(const t_trxframe& fr, const t_pbc* pbc);
    /*! \brief
     * Function that calculate search patterns in Secondary Structure Map
     *
     * Assignment done based on original DSSP algo
     */
    void PatternSearch();
    /*! \brief
     * Function that Checks if H-Bond exist
     *
     * Check done according to DSSP algo
     */
    static bool  isHbondExist(const BackboneAtomIndexes& resA,
                              const BackboneAtomIndexes& resB,
                              const t_trxframe&          fr,
                              const t_pbc*               pbc);
    static float CalculateAtomicDistances(int A, int B, const t_trxframe& fr, const t_pbc* pbc);
};

void Dssp::calculateBends(const t_trxframe& fr, const t_pbc* pbc)
{
    const float benddegree{ 70.0 }, maxdist{ 2.5 };
    float       degree{ 0 }, vdist{ 0 }, vprod{ 0 };
    gmx::RVec   a{ 0, 0, 0 }, b{ 0, 0, 0 };
    bendmap_.resize(0);
    breakmap_.resize(0);
    bendmap_.resize(IndexMap_.size(), 0);
    breakmap_.resize(IndexMap_.size(), 0);
    for (size_t i = 0; i < IndexMap_.size() - 1; ++i)
    {
        if (CalculateAtomicDistances(IndexMap_[i].getIndex(BackboneAtomTypes::AtomC),
                                     IndexMap_[i + 1].getIndex(BackboneAtomTypes::AtomN),
                                     fr,
                                     pbc)
            > maxdist)
        {
            breakmap_[i]     = 1;
            breakmap_[i + 1] = 1;
        }
    }
    for (size_t i = 2; i < IndexMap_.size() - 2; ++i)
    {
        if (breakmap_[i - 1] || breakmap_[i] || breakmap_[i + 1])
        {
            continue;
        }
        for (int j = 0; j < 3; ++j)
        {
            a[j] = fr.x[IndexMap_[i].getIndex(BackboneAtomTypes::AtomCA)][j]
                   - fr.x[IndexMap_[i - 2].getIndex(BackboneAtomTypes::AtomCA)][j];
            b[j] = fr.x[IndexMap_[i + 2].getIndex(BackboneAtomTypes::AtomCA)][j]
                   - fr.x[IndexMap_[i].getIndex(BackboneAtomTypes::AtomCA)][j];
        }
        vdist = (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2]);
        vprod = CalculateAtomicDistances(IndexMap_[i - 2].getIndex(BackboneAtomTypes::AtomCA),
                                         IndexMap_[i].getIndex(BackboneAtomTypes::AtomCA),
                                         fr,
                                         pbc)
                * gmx::c_angstrom / gmx::c_nano
                * CalculateAtomicDistances(IndexMap_[i].getIndex(BackboneAtomTypes::AtomCA),
                                           IndexMap_[i + 2].getIndex(BackboneAtomTypes::AtomCA),
                                           fr,
                                           pbc)
                * gmx::c_angstrom / gmx::c_nano;
        degree = std::acos(vdist / vprod) * gmx::c_rad2Deg;
        if (degree > benddegree)
        {
            bendmap_[i] = 1;
        }
    }
}

void Dssp::PatternSearch()
{
    std::vector<std::size_t>              v1_;
    std::vector<std::vector<std::size_t>> nturnsmap_;
    std::vector<std::vector<std::size_t>> Bridges_, AntiBridges_;
    v1_.resize(0);
    v1_.resize(HBondsMap_.size(), 0);
    SecondaryStructuresMap_.resize(0);
    SecondaryStructuresMap_.resize(HBondsMap_.size(), 0);
    nturnsmap_.resize(0);
    nturnsmap_.resize(6, v1_);
    v1_.resize(0);
    Bridges_.resize(0);
    Bridges_.resize(HBondsMap_.size(), v1_);
    AntiBridges_.resize(0);
    AntiBridges_.resize(HBondsMap_.size(), v1_);
    for (size_t i = 0; i < HBondsMap_.front().size(); ++i)
    {
        if (bendmap_[i])
        {
            SecondaryStructuresMap_[i] = 7;
        }
    }
    for (size_t n = 3; n <= 5; ++n)
    {
        for (size_t i = 0; i + n < HBondsMap_.front().size(); ++i)
        {
            if (HBondsMap_[i][i + n])
            {
                nturnsmap_[n - 3][i] = n;
                for (size_t j{ 1 }; j < n; ++j)
                {
                    if ((i + j) < nturnsmap_.front().size())
                    {
                        SecondaryStructuresMap_[i + j] = 6;
                    }
                }
            }
            // Go back
            if (HBondsMap_[i + n][i])
            {
                nturnsmap_[n][i + n] = n;
                for (size_t j{ 1 }; j < n; ++j)
                {
                    if ((i + n) >= j)
                    {
                        SecondaryStructuresMap_[i + n - j] = 6;
                    }
                }
            }
        }
    }
    for (size_t i = 1; i < HBondsMap_.front().size() - 1; ++i)
    {
        for (size_t j = 1; j < HBondsMap_.front().size() - 1; ++j)
        {
            if (std::abs(static_cast<int>(i) - static_cast<int>(j)) > 2)
            {
                if ((HBondsMap_[i - 1][j] && HBondsMap_[j][i + 1])
                    || (HBondsMap_[j - 1][i] && HBondsMap_[i][j + 1]))
                {
                    Bridges_[i].push_back(j);
                }
                if ((HBondsMap_[i][j] && HBondsMap_[j][i])
                    || (HBondsMap_[i - 1][j + 1] && HBondsMap_[j - 1][i + 1]))
                {
                    AntiBridges_[i].push_back(j);
                }
            }
        }
    }
    for (int n = 3; n <= 5; ++n)
    {
        for (size_t i = 1; i < HBondsMap_.front().size() - 1; ++i)
        {
            if (nturnsmap_[n - 3][i - 1] && nturnsmap_[n - 3][i])
            {
                for (int j = 0; j < n; ++j)
                {
                    if ((j + i) < SecondaryStructuresMap_.size())
                    {
                        SecondaryStructuresMap_[j + i] = n;
                    }
                }
            }
        }
        for (int i = (HBondsMap_.front().size() - 2); i >= 0; --i)
        {
            if (nturnsmap_[n][i] && nturnsmap_[n][i + 1])
            {
                for (int j = 0; j < n; ++j)
                {
                    if ((i - j) >= 0)
                    {
                        SecondaryStructuresMap_[i - j] = n;
                    }
                }
            }
        }

        if (n == 3)
        {
            for (size_t i = 0; i < HBondsMap_.front().size(); ++i)
            {
                if ((!Bridges_[i].empty() || !AntiBridges_[i].empty()))
                {
                    SecondaryStructuresMap_[i] = 1;
                }
            }
            for (size_t i = 2; i < HBondsMap_.front().size() - 3; ++i)
            {
                for (size_t j = i - 2; j <= (i + 2); ++j)
                {
                    if (j == i)
                    {
                        continue;
                    }
                    else
                    {
                        if (!Bridges_[i].empty() || !Bridges_[j].empty())
                        {
                            for (size_t i_resi = 0; i_resi < Bridges_[i].size(); ++i_resi)
                            {
                                for (size_t j_resi = 0; j_resi < Bridges_[j].size(); ++j_resi)
                                {
                                    if (abs(static_cast<int>(Bridges_[i][i_resi])
                                            - static_cast<int>(Bridges_[j][j_resi]))
                                        && (abs(static_cast<int>(Bridges_[i][i_resi])
                                                - static_cast<int>(Bridges_[j][j_resi]))
                                            < 5))
                                    {
                                        if (j < i)
                                        {
                                            for (size_t k = 0; k <= i - j; ++k)
                                            {
                                                SecondaryStructuresMap_[j + k] = 2;
                                            }
                                        }
                                        else
                                        {
                                            for (size_t k = 0; k <= j - i; ++k)
                                            {
                                                SecondaryStructuresMap_[i + k] = 2;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        if (!AntiBridges_[i].empty() || !AntiBridges_[j].empty())
                        {
                            for (size_t i_resi = 0; i_resi < AntiBridges_[i].size(); ++i_resi)
                            {
                                for (size_t j_resi = 0; j_resi < AntiBridges_[j].size(); ++j_resi)
                                {
                                    if (abs(static_cast<int>(AntiBridges_[i][i_resi])
                                            - static_cast<int>(AntiBridges_[j][j_resi]))
                                        && (abs(static_cast<int>(AntiBridges_[i][i_resi])
                                                - static_cast<int>(AntiBridges_[j][j_resi]))
                                            < 5))
                                    {
                                        if (j < i)
                                        {
                                            for (size_t k = 0; k <= i - j; ++k)
                                            {
                                                SecondaryStructuresMap_[j + k] = 2;
                                            }
                                        }
                                        else
                                        {
                                            for (size_t k = 0; k <= j - i; ++k)
                                            {
                                                SecondaryStructuresMap_[i + k] = 2;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

bool Dssp::isHbondExist(const BackboneAtomIndexes& resA,
                        const BackboneAtomIndexes& resB,
                        const t_trxframe&          fr,
                        const t_pbc*               pbc)
{
    /*
     * DSSP uses eq from dssp 2.x
     * kCouplingConstant = 27.888,  //  = 332 * 0.42 * 0.2
     * E = k * (1/rON + 1/rCH - 1/rOH - 1/rCN) where CO comes from one AA and NH from another
     * if R is in A
     * Hbond if E < -0.5
     */
    const float kCouplingConstant = 27.888,
                HBondEnergyCutOff = -0.5, // from dssp
            minimalAtomDistance   = 0.5,  // from original dssp
            minimalCAdistance     = 9.0,  // from original dssp
            minEnergy             = -9.9; // from original dssp algo in A
    float HbondEnergy             = 0.;
    float distanceON = 0., distanceCH = 0., distanceOH = 0., distanceCN = 0.;
    distanceON = CalculateAtomicDistances(
            resA.getIndex(BackboneAtomTypes::AtomO), resB.getIndex(BackboneAtomTypes::AtomN), fr, pbc);
    distanceCH = CalculateAtomicDistances(
            resA.getIndex(BackboneAtomTypes::AtomC), resB.getIndex(BackboneAtomTypes::AtomH), fr, pbc);
    distanceOH = CalculateAtomicDistances(
            resA.getIndex(BackboneAtomTypes::AtomO), resB.getIndex(BackboneAtomTypes::AtomH), fr, pbc);
    distanceCN = CalculateAtomicDistances(
            resA.getIndex(BackboneAtomTypes::AtomC), resB.getIndex(BackboneAtomTypes::AtomN), fr, pbc);

    if (resA.getIndex(BackboneAtomTypes::AtomC) && resA.getIndex(BackboneAtomTypes::AtomO)
        && resB.getIndex(BackboneAtomTypes::AtomN) && resB.getIndex(BackboneAtomTypes::AtomH))
    {
        if (CalculateAtomicDistances(
                    resA.getIndex(BackboneAtomTypes::AtomCA), resB.getIndex(BackboneAtomTypes::AtomCA), fr, pbc)
            < minimalCAdistance)
        {
            if ((distanceON < minimalAtomDistance) || (distanceCH < minimalAtomDistance)
                || (distanceOH < minimalAtomDistance) || (distanceCN < minimalAtomDistance))
            {
                HbondEnergy = minEnergy;
            }
            else
            {
                HbondEnergy =
                        kCouplingConstant
                        * ((1 / distanceON) + (1 / distanceCH) - (1 / distanceOH) - (1 / distanceCN));
            }

            return HbondEnergy < HBondEnergyCutOff;
        }
        else
        {
            return false;
        }
    }
    else
    {
        return false;
    }
}

float Dssp::CalculateAtomicDistances(int A, int B, const t_trxframe& fr, const t_pbc* pbc)
{
    gmx::RVec r{ 0, 0, 0 };
    pbc_dx(pbc, fr.x[B], fr.x[A], r.as_vec());
    return r.norm() * gmx::c_nm2A;
}

void Dssp::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] calculate protein secondary structure via DSSP algo",
    };
    options->addOption(FileNameOption("o")
                               .outputFile()
                               .store(&fnmDSSPOut_)
                               .required()
                               .defaultBasename("dssp")
                               .filetype(OptionFileType::GenericData)
                               .description("Filename for DSSP output"));
    options->addOption(RealOption("cutoff").store(&cutoff_).required().defaultValue(1.0).description(
            "cutoff for neighbour search"));
    options->addOption(SelectionOption("sel").store(&sel_).defaultSelectionText("Protein").description(
            "Group for DSSP"));
    settings->setHelpText(desc);
}


void Dssp::optionsFinished(TrajectoryAnalysisSettings* /* settings */) {}


void Dssp::initAnalysis(const TrajectoryAnalysisSettings& /* settings */, const TopologyInformation& top)
{
    BackboneAtomIndexes backboneAtoms;
    int                 i = 0;
    int resicompare{ top.atoms()->atom[static_cast<std::size_t>(*(sel_.atomIndices().begin()))].resind };
    AtomResi_.resize(0);
    IndexMap_.resize(0);
    IndexMap_.push_back(backboneAtoms);
    for (gmx::ArrayRef<const int>::iterator ai{ sel_.atomIndices().begin() };
         (ai < sel_.atomIndices().end());
         ++ai)
    {

        if (resicompare != top.atoms()->atom[static_cast<std::size_t>(*ai)].resind)
        {
            ++i;
            resicompare = top.atoms()->atom[static_cast<std::size_t>(*ai)].resind;
            IndexMap_.push_back(backboneAtoms);
        }
        AtomResi_.push_back(i);
        std::string atomname(*(top.atoms()->atomname[static_cast<std::size_t>(*ai)]));
        if (atomname == BackboneAtomTypeNames[BackboneAtomTypes::AtomCA])
        {
            IndexMap_[i].addAtomIndex(BackboneAtomTypes::AtomCA, static_cast<std::size_t>(*ai));
        }
        else if (atomname == BackboneAtomTypeNames[BackboneAtomTypes::AtomC])
        {
            IndexMap_[i].addAtomIndex(BackboneAtomTypes::AtomC, static_cast<std::size_t>(*ai));
        }
        else if (atomname == BackboneAtomTypeNames[BackboneAtomTypes::AtomO])
        {
            IndexMap_[i].addAtomIndex(BackboneAtomTypes::AtomO, static_cast<std::size_t>(*ai));
        }
        else if (atomname == BackboneAtomTypeNames[BackboneAtomTypes::AtomN])
        {
            IndexMap_[i].addAtomIndex(BackboneAtomTypes::AtomN, static_cast<std::size_t>(*ai));
        }
        else if (atomname == BackboneAtomTypeNames[BackboneAtomTypes::AtomH])
        {
            IndexMap_[i].addAtomIndex(BackboneAtomTypes::AtomH, static_cast<std::size_t>(*ai));
        }
    }
}

void Dssp::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* /* pdata */)
{
    // store positions of CA atoms to use them for nbSearch
    std::vector<gmx::RVec> positionsCA_;
    std::string            sspPattern_;
    for (auto& i : IndexMap_)
    {
        positionsCA_.emplace_back(fr.x[i.getIndex(BackboneAtomTypes::AtomCA)]);
    }
    // resize HBondsMap_
    HBondsMap_.resize(0);
    HBondsMap_.resize(IndexMap_.size(), std::vector<bool>(IndexMap_.size(), false));

    // Init nbSearch
    AnalysisNeighborhood nb_;
    nb_.setCutoff(cutoff_);
    AnalysisNeighborhoodPositions       nbPos_(positionsCA_);
    gmx::AnalysisNeighborhoodSearch     start      = nb_.initSearch(pbc, nbPos_);
    gmx::AnalysisNeighborhoodPairSearch pairSearch = start.startPairSearch(nbPos_);
    gmx::AnalysisNeighborhoodPair       pair;
    while (pairSearch.findNextPair(&pair))
    {
        HBondsMap_[pair.refIndex()][pair.testIndex()] =
                isHbondExist(IndexMap_[pair.refIndex()], IndexMap_[pair.testIndex()], fr, pbc);
    }

    calculateBends(fr, pbc);
    PatternSearch();
    sspPattern_.resize(0);
    for (size_t i = 0; i < SecondaryStructuresMap_.size(); ++i)
    {
        switch (SecondaryStructuresMap_[i])
        {
            case 0: sspPattern_.push_back('~'); break;
            case 1: sspPattern_.push_back('B'); break;
            case 2: sspPattern_.push_back('E'); break;
            case 3: sspPattern_.push_back('G'); break;
            case 4: sspPattern_.push_back('H'); break;
            case 5: sspPattern_.push_back('I'); break;
            case 6: sspPattern_.push_back('T'); break;
            case 7: sspPattern_.push_back('S'); break;
            default: sspPattern_.push_back('?'); break;
        }

        if (i != SecondaryStructuresMap_.size() - 1)
        {
            if (breakmap_[i] && breakmap_[i + 1])
            {
                sspPattern_.push_back('=');
            }
        }
    }
    DsspStorage_.addData(frnr, sspPattern_);
}

void Dssp::finishAnalysis(int /*nframes*/) {}


void Dssp::writeOutput()
{
    std::vector<DsspStorageFrame> dataOut_;
    FILE*                         fp_;
    fp_      = gmx_ffopen(fnmDSSPOut_, "w");
    dataOut_ = DsspStorage_.getData();
    for (auto& i : dataOut_)
    {
        std::fprintf(fp_, "%s\n", i.dsspData.c_str());
    }
    gmx_ffclose(fp_);
}

} // namespace

const char DsspInfo::name[]             = "dssp";
const char DsspInfo::shortDescription[] = "Calculate protein secondary structure via DSSP algo";

TrajectoryAnalysisModulePointer DsspInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Dssp);
}

} // namespace analysismodules

} // namespace gmx
