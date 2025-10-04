/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Implements gmx::analysismodules::Hbond2.
 *
 * \author Sergey Gorelov <gorelov_sv@pnpi.nrcki.ru>
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "hbond.h"

#include <cstdio>

#include <algorithm>
#include <array>
#include <filesystem>
#include <initializer_list>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/math/units.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/options/options.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/nbsearch.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/idef.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/vec.h"
#include "gromacs/utility/vectypes.h"

#include "isotope.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

struct HBond
{
    HBond(int d, int a, bool aiad, int h = -1) :
        donor(d), acceptor(a), acceptorIsAlsoDonor(aiad), hydrogen(h)
    {
    }
    int   donor;
    int   acceptor;
    bool  acceptorIsAlsoDonor;
    int   hydrogen;
    bool  operator<(const HBond& other) const;
    bool  operator==(const HBond& other) const;
    HBond operator!() const;
};

bool HBond::operator<(const HBond& other) const
{
    return (this->donor < other.donor) || (this->donor == other.donor && this->acceptor < other.acceptor)
           || (this->donor == other.donor && this->acceptor == other.acceptor
               && this->hydrogen < other.hydrogen);
}

bool HBond::operator==(const HBond& other) const
{
    return this->donor == other.donor && this->acceptor == other.acceptor && this->hydrogen == other.hydrogen;
}

HBond HBond::operator!() const
{
    return HBond(this->acceptor, this->donor, true);
}

//! Structure that contains storage information from different frames.
struct HbondStorageFrame
{
    HbondStorageFrame(int frnr, std::vector<HBond> data) :
        frameNumber_(frnr), hbondData_(std::move(data))
    {
    }
    //! Frame number.
    int frameNumber_ = 0;
    //! Frame hbond data.
    std::vector<HBond> hbondData_;
};

/*! \brief
 * Class that stores frame information in storage and, upon request, can return it.
 */
class HbondStorage
{
public:
    /*! \brief
     * Function that stores frame information in storage.
     */
    void addData(int frnr, const std::vector<HBond>& data);
    /*! \brief
     * Function that returns frame information from storage.
     */
    const std::vector<HbondStorageFrame>& getData() const;

private:
    /*! \brief
     * Vector that contains information from different frames.
     */
    std::vector<HbondStorageFrame> data_;
};

void HbondStorage::addData(int frnr, const std::vector<HBond>& data)
{
    data_.emplace_back(frnr, data);
}

const std::vector<HbondStorageFrame>& HbondStorage::getData() const
{
    return data_;
}

struct t_acceptor
{
    t_acceptor(int atomIndex) : ai(atomIndex) {}
    // Acceptor Atom Qndex
    int  ai;
    bool isAlsoDonor = false;
};

const std::size_t c_maxHydrogensWithOxygen = 3;

const std::size_t c_maxHydrogensWithNitrogen = 4;

const std::size_t c_maxHydrogenBonds = 2;

struct t_donor
{
    t_donor(int atomIndex, std::string donorElement) : ai(atomIndex), de(std::move(donorElement)) {}
    // Donor Atom Index
    int ai;
    // Associated Hydrogen Atoms Indices
    std::vector<int> h_atoms;
    // Donor Element
    std::string de;
    void        addHydrogen(int hi);
};

void t_donor::addHydrogen(const int hi)
{
    if (this->de == "O" && h_atoms.size() > c_maxHydrogensWithOxygen)
    {
        GMX_THROW(InconsistentInputError("Donor " + std::to_string(this->ai) + " has more than "
                                         + std::to_string(c_maxHydrogensWithOxygen)
                                         + " covalent bonds with hydrogens!"));
    }
    else if (this->de == "N" && h_atoms.size() > c_maxHydrogensWithNitrogen)
    {
        GMX_THROW(InconsistentInputError("Donor " + std::to_string(this->ai) + " has more than "
                                         + std::to_string(c_maxHydrogensWithNitrogen)
                                         + " covalent bonds with hydrogens!"));
    }
    else
    {
        h_atoms.push_back(hi);
    }
}

struct t_info
{
    std::vector<t_acceptor> acceptors;
    std::vector<t_donor>    donors;
    const Selection*        selectionPtr  = nullptr;
    std::string             selectionType = "Selection";
    std::vector<int>        atomIndices;
};

class Hbond : public TrajectoryAnalysisModule
{
public:
    Hbond();
    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;
    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    static bool isInSelection(int ai, const std::vector<int>& selection);
    bool        isValidDonorElement(const std::string& element) const;
    bool        isValidAcceptorElement(const std::string& element) const;
    /*! \brief Check partial overlapping between two parsed selections
     * and throw an error when groups are overlapping but are not equal. */
    void checkOverlap() const;
    void searchAcceptors(const TopologyInformation& top,
                         t_info*                    selectionTool,
                         const std::vector<int>*    selection) const;
    void searchDonors(const TopologyInformation& top, t_info* selectionTool, const std::vector<int>* selection);
    void                       searchDonorsAndAcceptors();
    static void                linkDA(t_info* selectionTool);
    void                       prepareForAnalysis(const TrajectoryAnalysisSettings& settings);
    std::vector<HBond>         prepareFrameData(const std::vector<HbondStorageFrame>& data) const;
    std::vector<HBond>         prepareFrameData(const std::vector<HBond>& data) const;
    Selection                  refSelection_, targetSelection_;
    std::string                fnmHbondOut_ = "hbond";
    std::string                fnmHbnumOut_;
    std::string                fnmHbdistOut_;
    std::string                fnmHbangOut_;
    std::string                fnmHbdanOut_;
    std::string                refSelectionType    = "Selection";
    std::string                targetSelectionType = "Selection";
    t_info                     refInfo_, targetInfo_;
    const TopologyInformation* localTopology                    = nullptr;
    bool                       isTwoDiffGroups_                 = false;
    bool                       perFrame_                        = false;
    bool                       mergeHydrogens_                  = false;
    bool                       staticSelectionRefInitialized    = false;
    bool                       staticSelectionTargetInitialized = false;
    real                       nbsearchCutoff_                  = 0.35;
    real                       hbDistCutoff_                    = 0.35;
    real                       hbAngleCutoff_                   = 30;
    std::vector<std::string>   donorElements                    = { "N", "O" };
    std::vector<std::string>   acceptorElements                 = { "N", "O" };
    HbondStorage               storage_;

    AnalysisData                             hbnum_;
    AnalysisData                             distances_;
    AnalysisData                             angles_;
    AnalysisData                             dan_;
    AnalysisDataSimpleHistogramModulePointer histogramModuleDist_;
    AnalysisDataSimpleHistogramModulePointer histogramModuleAng_;
};

Hbond::Hbond() :
    histogramModuleDist_(new AnalysisDataSimpleHistogramModule),
    histogramModuleAng_(new AnalysisDataSimpleHistogramModule)
{
    registerAnalysisDataset(&hbnum_, "hbnum");

    registerAnalysisDataset(&distances_, "hbdist");
    distances_.addModule(histogramModuleDist_);
    registerBasicDataset(&histogramModuleDist_->averager(), "histogram_dist");

    registerAnalysisDataset(&angles_, "hbang");
    angles_.addModule(histogramModuleAng_);
    registerBasicDataset(&histogramModuleAng_->averager(), "histogram_ang");

    registerAnalysisDataset(&dan_, "hbdan");
}

void Hbond::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] allows using geometric definition of hydrogen bonds to define them "
        "throughout the structure.[PAR]"
        "[TT]-r[tt] specifies reference selection, relative to which the search for hydrogen bonds "
        "in target selection will develop. Note that all atoms in reference and target selections "
        "should be either absolutely identical or non-overlapping at all. Accepts dynamic "
        "selection.[PAR]"
        "[TT]-t[tt] specifies target selection, relative to which the search for hydrogen bonds in "
        "reference selection will develop. Note that all atoms in reference and target selections "
        "should be either absolutely identical or non-overlapping at all. Accepts dynamic "
        "selection.[PAR]"
        "[TT]-m[tt] forces to merge together information in output index file about hydrogen bonds "
        "if they differ only in hydrogen indices. This also means that information about hydrogen "
        "atoms in the hydrogen bonds would not be written in output index file at all.[PAR]"
        "[TT]-pf[tt] forces to write hydrogen bonds for each frame separately instead of writing "
        "hydrogen bonds for the whole system. Each information about hydrogen bonds in new frame "
        "will be stored in its own section of the output index file.[PAR]"
        "[TT]-cutoff[tt] is a real value that defines distance from donor to acceptor "
        "(and vice versa) that used in neighbor search. Minimum (and also recommended) value is "
        "0.35.[PAR]"
        "[TT]-hbr[tt] Sets the cutoff that is used when calculating hydrogen bond distances. "
        "Recommended value: 0.35. [PAR]"
        "[TT]-hba[tt] Sets the cutoff that is used when calculating hydrogen bond angles. "
        "Recommended value: 30. [PAR]"
        "[TT]-de[tt] Specifies the atomic elements that will be selected from the topology to "
        "check if a given element is a potential hydrogen bond donor. [PAR]"
        "[TT]-ae[tt] Specifies the atomic elements that will be selected from the topology to "
        "check if a given element is a potential hydrogen bond acceptor. [PAR]"
        "[TT]-num[tt] allows you to get a plot of the number of hydrogen bonds as a function of "
        "time at the output.[PAR]"
        "[TT]-dist[tt] allows you to get a plot of the distance distribution of all hydrogen bonds "
        "at the output.[PAR]"
        "[TT]-ang[tt] allows you to get a plot of the angular distribution of all hydrogen bonds "
        "at the output.[PAR]"
        "[TT]-dan[tt] allows you to get a plot of the number of analyzed donors and acceptors for "
        "each frame at the output.[PAR]",
        "[PAR]",
        "Note that this is a new implementation of the hbond utility added in",
        "GROMACS 2024. If you need the old one, use [TT]gmx hbond-legacy[tt]."
    };


    options->addOption(
            FileNameOption("o")
                    .filetype(OptionFileType::AtomIndex)
                    .outputFile()
                    .required()
                    .store(&fnmHbondOut_)
                    .defaultBasename("hbond")
                    .description("Index file that contains selected groups', acceptors', "
                                 "donors' and hydrogens' indices and hydrogen bond "
                                 "pairs between or within selected groups."));
    options->addOption(FileNameOption("num")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnmHbnumOut_)
                               .defaultBasename("hbnum")
                               .description("Number of hydrogen bonds as a function of time."));
    options->addOption(FileNameOption("dist")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnmHbdistOut_)
                               .defaultBasename("hbdist")
                               .description("Distance distribution of all hydrogen bonds."));
    options->addOption(FileNameOption("ang")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnmHbangOut_)
                               .defaultBasename("hbang")
                               .description("Angle distribution of all hydrogen bonds."));
    options->addOption(
            FileNameOption("dan")
                    .filetype(OptionFileType::Plot)
                    .outputFile()
                    .store(&fnmHbdanOut_)
                    .defaultBasename("hbdan")
                    .description("Number of donors and acceptors analyzed for each frame."));
    options->addOption(SelectionOption("r")
                               .store(&refSelection_)
                               .required()
                               .description("Reference selection, relative to which the search for "
                                            "hydrogen bonds in target "
                                            "selection will develop."));
    options->addOption(
            SelectionOption("t").store(&targetSelection_).required().description("Target selection, relative to which the search for hydrogen bonds in reference selection will develop."));
    options->addOption(
            BooleanOption("m").store(&mergeHydrogens_).description("Merge together information about hydrogen bonds if they differ only in hydrogen indices."));
    options->addOption(
            BooleanOption("pf").store(&perFrame_).description("Write hydrogen bonds for each frame separately instead of writing hydrogen bonds for the whole system."));
    options->addOption(
            RealOption("cutoff").store(&nbsearchCutoff_).description("Distance from donor to acceptor (and vice versa) that used in neighbor search (nm). Must be > 0."));
    options->addOption(
            RealOption("hbr").store(&hbDistCutoff_).description("Hydrogen bond cutoff distance, between donor and acceptor (nm). The value must not exceed the neighbor search cutoff and must be > 0."));
    options->addOption(
            RealOption("hba").store(&hbAngleCutoff_).description("A-D-H hydrogen bond cutoff angle (degrees). Must be > 0."));
    options->addOption(
            StringOption("de").storeVector(&donorElements).multiValue().description("Donor elements. Default elements: N, O."));
    options->addOption(
            StringOption("ae").storeVector(&acceptorElements).multiValue().description("Acceptor elements. Default elements: N, O."));
    settings->setHelpText(desc);
    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}

void Hbond::optionsFinished(TrajectoryAnalysisSettings* /* settings */)
{
    if (nbsearchCutoff_ <= 0)
    {
        GMX_THROW(InconsistentInputError("Invalid cutoff value. It must be > 0."));
    }
    if (hbDistCutoff_ <= 0)
    {
        GMX_THROW(InconsistentInputError("Invalid hbond distance cutoff value. It must be > 0."));
    }
    if (hbAngleCutoff_ <= 0)
    {
        GMX_THROW(InconsistentInputError("Invalid hbond angle cutoff value. It must be > 0."));
    }
    if (nbsearchCutoff_ < hbDistCutoff_)
    {
        GMX_THROW(
                InconsistentInputError("Hydrogen bond cutoff distance should not exceed neighbor "
                                       "search cutoff distance."));
    }
}

bool Hbond::isInSelection(int ai, const std::vector<int>& selection)
{
    return std::find(selection.begin(), selection.end(), ai) != selection.end();
}

bool Hbond::isValidDonorElement(const std::string& element) const
{
    return std::find(donorElements.begin(), donorElements.end(), element) != donorElements.end();
}


bool Hbond::isValidAcceptorElement(const std::string& element) const
{
    return std::find(acceptorElements.begin(), acceptorElements.end(), element) != acceptorElements.end();
}


void Hbond::checkOverlap() const
{
    printf("\nChecking for overlap in atoms between %s (%d atoms) and %s (%d atoms)\n",
           refSelection_.name(),
           refSelection_.atomCount(),
           targetSelection_.name(),
           targetSelection_.atomCount());
    bool isIdentical = false;
    bool isDifferent = false;

    auto i = refInfo_.atomIndices.begin(), j = targetInfo_.atomIndices.begin();
    while ((i != refInfo_.atomIndices.end() && j != targetInfo_.atomIndices.end())
           && !(isIdentical && isDifferent))
    {
        if (*i == *j)
        {
            isIdentical = true;
            ++i;
            ++j;
        }
        else if (*i < *j)
        {
            isDifferent = true;
            ++i;
        }
        else
        {
            isDifferent = true;
            ++j;
        }
    }

    GMX_RELEASE_ASSERT(not(isIdentical && isDifferent),
                       (std::string("Partial overlap between groups '") + refSelection_.name()
                        + "' and '" + targetSelection_.name() + "'")
                               .c_str());
}

void Hbond::searchAcceptors(const TopologyInformation& top,
                            t_info*                    selectionTool,
                            const std::vector<int>*    selection) const
{
    for (auto ai = selection->begin(); ai != selection->end(); ++ai)
    {
        std::string element                = top.atoms()->atom[static_cast<std::size_t>(*ai)].elem;
        bool        acceptorElementIsValid = isValidAcceptorElement(element);
        if (acceptorElementIsValid && isInSelection(*ai, *selection))
        {
            selectionTool->acceptors.emplace_back(*ai);
        }
    }

    if (!selectionTool->acceptors.empty())
    {
        std::sort(selectionTool->acceptors.begin(),
                  selectionTool->acceptors.end(),
                  [](t_acceptor a1, t_acceptor a2) { return a1.ai < a2.ai; });
    }
    printf("%s '%s' has %zu acceptors",
           selectionTool->selectionType.c_str(),
           selectionTool->selectionPtr->name(),
           selectionTool->acceptors.size());
}

void Hbond::searchDonors(const TopologyInformation& top, t_info* selectionTool, const std::vector<int>* selection)
{
    std::vector<std::pair<int, int>> dhUnsorted;
    for (const auto func_type : gmx::EnumerationWrapper<InteractionFunction>{})
    {
        const InteractionList* interaction = &(top.expandedTopology()->idef.il[func_type]);
        if (func_type == InteractionFunction::PositionRestraints
            || func_type == InteractionFunction::FlatBottomedPositionRestraints)
        {
            /* We don't need posre interactions for hbonds anyway.*/
            continue;
        }
        for (std::size_t i = 0; i < interaction->iatoms.size();
             i +=
             interaction_function[top.expandedTopology()->idef.functype[interaction->iatoms[i]]].nratoms + 1)
        {
            /* next function */
            if (func_type != top.expandedTopology()->idef.functype[interaction->iatoms[i]])
            {
                fprintf(stderr, "Error in func_type '%s'", interaction_function[func_type].longname);
                continue;
            }
            if (func_type == InteractionFunction::SETTLE)
            {
                int nr1 = interaction->iatoms[i + 1]; // DONOR
                int nr2 = interaction->iatoms[i + 2]; // HYDRO 1
                int nr3 = interaction->iatoms[i + 3]; // HYDRO 2
                if (isInSelection(nr1, *selection))
                {
                    if (isInSelection(nr2, *selection))
                    {
                        dhUnsorted.emplace_back(nr1, nr2);
                    }
                    if (isInSelection(nr3, *selection))
                    {
                        dhUnsorted.emplace_back(nr1, nr3);
                    }
                }
            }
            else if (interaction_function[(func_type)].nratoms == 2
                     && (interaction_function[(func_type)].flags & IF_CHEMBOND))
            {
                for (std::size_t j = 0; j < 2; ++j)
                {
                    int         nr1      = interaction->iatoms[i + 1 + j]; // POTENTIAL HYDRO
                    int         nr2      = interaction->iatoms[i + 2 - j]; // POTENTIAL DONOR
                    std::string element1 = top.atoms()->atom[static_cast<std::size_t>(nr1)].elem;
                    std::string element2 = top.atoms()->atom[static_cast<std::size_t>(nr2)].elem;
                    if (element1 == "H" && (isValidDonorElement(element2)))
                    {
                        if (isInSelection(nr1, *selection) && isInSelection(nr2, *selection))
                        {
                            dhUnsorted.emplace_back(nr2, nr1);
                        }
                    }
                }
            }
        }
    }
    if (!dhUnsorted.empty())
    {
        std::sort(dhUnsorted.begin(),
                  dhUnsorted.end(),
                  [](std::pair<int, int> d1, std::pair<int, int> d2)
                  { return d1.first < d2.first || (d1.first == d2.first && d1.second < d2.second); });

        int memoryIndex = dhUnsorted.begin()->first - 1;

        for (const auto& pi : dhUnsorted)
        {
            if (memoryIndex != pi.first)
            {
                std::string donorElement = top.atoms()->atom[static_cast<std::size_t>(pi.first)].elem;
                selectionTool->donors.emplace_back(pi.first, donorElement);
                memoryIndex = pi.first;
            }
            selectionTool->donors.back().addHydrogen(pi.second);
        }
    }
    printf(" and %zu donors.\n", selectionTool->donors.size());
}

void Hbond::linkDA(t_info* selectionTool)
{
    if (!selectionTool->acceptors.empty() || !selectionTool->donors.empty())
    {

        auto i = selectionTool->acceptors.begin();
        auto j = selectionTool->donors.begin();

        while (i != selectionTool->acceptors.end() && j != selectionTool->donors.end())
        {
            if (i->ai == j->ai)
            {
                i->isAlsoDonor = true;
            }

            if (i->ai <= j->ai)
            {
                ++i;
            }
            else
            {
                ++j;
            }
        }
    }
    else
    {
        GMX_THROW(InconsistentInputError(
                selectionTool->selectionType + " " + std::string(selectionTool->selectionPtr->name())
                + "' has no donors AND has no acceptors! Nothing to be done."));
    }
}


void Hbond::prepareForAnalysis(const TrajectoryAnalysisSettings& settings)
{
    if (!fnmHbnumOut_.empty())
    {
        hbnum_.setColumnCount(0, 1);
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnmHbnumOut_);
        plotm->setTitle("Number of hydrogen bonds");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Hbonds");
        plotm->appendLegend("Hydrogen bonds");
        plotm->setYFormat(10, 0);
        hbnum_.addModule(plotm);
    }
    if (!fnmHbdistOut_.empty())
    {
        const double histogramMin = 0;
        const double histogramMax = 0.35;
        histogramModuleDist_->init(
                histogramFromRange(histogramMin, histogramMax).binWidth(0.005).includeAll());

        distances_.setColumnCount(
                0, ((refInfo_.donors.size() + targetInfo_.donors.size()) * c_maxHydrogenBonds));

        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnmHbdistOut_);
        plotm->setTitle("Hydrogen Bond Distance Distribution Histogram");
        plotm->setXLabel("Distance (nm)");
        plotm->setYLabel("Probability");
        plotm->setXFormat(10, 4);
        plotm->setYFormat(10, 4);
        histogramModuleDist_->averager().addModule(plotm);
    }
    if (!fnmHbangOut_.empty())
    {
        const double histogramMin = 0;
        const double histogramMax = 30;
        histogramModuleAng_->init(histogramFromRange(histogramMin, histogramMax).binWidth(1.0).includeAll());

        angles_.setColumnCount(
                0, ((refInfo_.donors.size() + targetInfo_.donors.size()) * c_maxHydrogenBonds));

        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnmHbangOut_);
        plotm->setTitle("Hydrogen Bond Angle Distribution Histogram");
        plotm->setXLabel("Hydrogen - Donor - Acceptor Angle (degrees)");
        plotm->setYLabel("Probability");
        plotm->setXFormat(10, 1);
        plotm->setYFormat(10, 7);
        histogramModuleAng_->averager().addModule(plotm);
    }
    if (!fnmHbdanOut_.empty())
    {
        dan_.setColumnCount(0, 2);
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnmHbdanOut_);
        plotm->setTitle("Number of Donors and Acceptors");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Donors/Acceptors");
        plotm->setXFormat(10, 0);
        plotm->setYFormat(10, 0);
        plotm->appendLegend("Number of Donors");
        plotm->appendLegend("Number of Acceptors");
        dan_.addModule(plotm);
    }
}


void Hbond::searchDonorsAndAcceptors()
{
    FILE* fp;
    if (refSelection_.isDynamic())
    {
        refInfo_.atomIndices.clear();
        refInfo_.acceptors.clear();
        refInfo_.donors.clear();
        refInfo_.atomIndices.resize(refSelection_.atomIndices().size());
        std::copy(refSelection_.atomIndices().begin(),
                  refSelection_.atomIndices().end(),
                  refInfo_.atomIndices.begin());
        if (!refSelection_.hasSortedAtomIndices())
        {
            std::sort(refInfo_.atomIndices.begin(), refInfo_.atomIndices.end());
        }
    }

    if (targetSelection_.isDynamic())
    {
        targetInfo_.atomIndices.clear();
        targetInfo_.acceptors.clear();
        targetInfo_.donors.clear();
        targetInfo_.atomIndices.resize(targetSelection_.atomIndices().size());

        std::copy(targetSelection_.atomIndices().begin(),
                  targetSelection_.atomIndices().end(),
                  targetInfo_.atomIndices.begin());

        if (!targetSelection_.hasSortedAtomIndices())
        {
            std::sort(targetInfo_.atomIndices.begin(), targetInfo_.atomIndices.end());
        }
    }

    checkOverlap();

    if (refInfo_.atomIndices[0] != targetInfo_.atomIndices[0])
    {
        isTwoDiffGroups_ = true;
    }

    if (refSelection_.isDynamic() || !staticSelectionRefInitialized)
    {
        searchAcceptors(*localTopology, &refInfo_, &(refInfo_.atomIndices));
        searchDonors(*localTopology, &refInfo_, &(refInfo_.atomIndices));
        linkDA(&refInfo_);
    }
    if (isTwoDiffGroups_ && (targetSelection_.isDynamic() || !staticSelectionTargetInitialized))
    {
        searchAcceptors(*localTopology, &targetInfo_, &(targetInfo_.atomIndices));
        searchDonors(*localTopology, &targetInfo_, &(targetInfo_.atomIndices));
        linkDA(&targetInfo_);
    }

    if ((!refSelection_.isDynamic() && !staticSelectionRefInitialized)
        || (!targetSelection_.isDynamic() && !staticSelectionTargetInitialized && isTwoDiffGroups_))
    {
        fp = gmx_ffopen(fnmHbondOut_, "w");
        gmx_ffclose(fp);
    }

    if (!refSelection_.isDynamic() && !staticSelectionRefInitialized)
    {
        fp = gmx_ffopen(fnmHbondOut_, "a");
        fprintf(fp, "[ %s ]", refSelection_.name());
        for (std::size_t i = 0; i < refInfo_.selectionPtr->atomIndices().size(); ++i)
        {
            fprintf(fp, "%c %4i", (i % 15 == 0 ? '\n' : ' '), refSelection_.atomIndices()[i] + 1);
        }
        fprintf(fp, "\n[ donors_hydrogens_%s ]", refSelection_.name());
        for (const auto& donor : refInfo_.donors)
        {
            fprintf(fp, "\n %4i", donor.ai + 1);
            for (std::size_t j = 0; j < donor.h_atoms.size(); ++j)
            {
                fprintf(fp, " %4i", donor.h_atoms[j] + 1);
            }
        }
        fprintf(fp, "\n[ acceptors_%s ]", refSelection_.name());
        for (std::size_t i = 0; i < refInfo_.acceptors.size(); ++i)
        {
            fprintf(fp, "%c %4i", (i % 15 == 0 ? '\n' : ' '), refInfo_.acceptors[i].ai + 1);
        }
        fprintf(fp, "\n");
        gmx_ffclose(fp);
        staticSelectionRefInitialized = true;
    }
    if (!targetSelection_.isDynamic() && !staticSelectionTargetInitialized && isTwoDiffGroups_)
    {
        fp = gmx_ffopen(fnmHbondOut_, "a");
        fprintf(fp, "[ %s ]", targetSelection_.name());
        for (std::size_t i = 0; i < targetInfo_.selectionPtr->atomIndices().size(); ++i)
        {
            fprintf(fp, "%c %4i", (i % 15 == 0 ? '\n' : ' '), targetSelection_.atomIndices()[i] + 1);
        }
        fprintf(fp, "\n[ donors_hydrogens_%s ]", targetSelection_.name());
        for (const auto& donor : targetInfo_.donors)
        {
            fprintf(fp, "\n %4i", donor.ai + 1);
            for (std::size_t j = 0; j < donor.h_atoms.size(); ++j)
            {
                fprintf(fp, " %4i", donor.h_atoms[j] + 1);
            }
        }
        fprintf(fp, "\n[ acceptors_%s ]", targetSelection_.name());
        for (std::size_t i = 0; i < targetInfo_.acceptors.size(); ++i)
        {
            fprintf(fp, "%c %4i", (i % 15 == 0 ? '\n' : ' '), targetInfo_.acceptors[i].ai + 1);
        }
        fprintf(fp, "\n");
        gmx_ffclose(fp);
        staticSelectionTargetInitialized = true;
    }
}

void Hbond::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top)
{
    if (!refSelection_.isValid() || !targetSelection_.isValid())
    {
        GMX_THROW(InvalidInputError("Invalid selection(s). Nothing to be done."));
    }

    localTopology = &top;

    refInfo_.selectionPtr    = &refSelection_;
    targetInfo_.selectionPtr = &targetSelection_;

    if (!refSelection_.isDynamic())
    {
        refInfo_.atomIndices.resize(refSelection_.atomIndices().size());
        std::copy(refSelection_.atomIndices().begin(),
                  refSelection_.atomIndices().end(),
                  refInfo_.atomIndices.begin());
        if (!refSelection_.hasSortedAtomIndices())
        {
            std::sort(refInfo_.atomIndices.begin(), refInfo_.atomIndices.end());
        }
    }
    else
    {
        refInfo_.selectionType = "Dynamic selection";
    }

    if (!targetSelection_.isDynamic())
    {
        targetInfo_.atomIndices.resize(targetSelection_.atomIndices().size());
        std::copy(targetSelection_.atomIndices().begin(),
                  targetSelection_.atomIndices().end(),
                  targetInfo_.atomIndices.begin());

        if (!targetSelection_.hasSortedAtomIndices())
        {
            std::sort(targetInfo_.atomIndices.begin(), targetInfo_.atomIndices.end());
        }
    }
    else
    {
        targetInfo_.selectionType = "Dynamic selection";
    }

    if (!refSelection_.isDynamic() && !targetSelection_.isDynamic())
    {
        searchDonorsAndAcceptors();
    }

    prepareForAnalysis(settings);
}


void Hbond::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata)
{
    if (refSelection_.isDynamic() || targetSelection_.isDynamic())
    {
        searchDonorsAndAcceptors();
        FILE* fp;
        fp = gmx_ffopen(fnmHbondOut_, "a");
        if (refSelection_.isDynamic())
        {
            // Print Selection Indices
            fprintf(fp, "[ %s_frame_%i ]", refSelection_.name(), frnr);
            for (std::size_t i = 0; i < refSelection_.atomIndices().size(); ++i)
            {
                fprintf(fp, "%c %4i", (i % 15 == 0 ? '\n' : ' '), refSelection_.atomIndices()[i] + 1);
            }
            // Print Donors Indices
            fprintf(fp, "\n[ donors_hydrogens_%s_frame_%i ]", refSelection_.name(), frnr);
            for (const auto& donor : refInfo_.donors)
            {
                fprintf(fp, "\n %4i", donor.ai + 1);
                for (std::size_t j = 0; j < donor.h_atoms.size(); ++j)
                {
                    fprintf(fp, " %4i", donor.h_atoms[j] + 1);
                }
            }
            // Print Acceptors Indices
            fprintf(fp, "\n[ acceptors_%s_frame_%i ]", refSelection_.name(), frnr);
            for (std::size_t i = 0; i < refInfo_.acceptors.size(); ++i)
            {
                fprintf(fp, "%c %4i", (i % 15 == 0 ? '\n' : ' '), refInfo_.acceptors[i].ai + 1);
            }
            fprintf(fp, "\n");
        }
        if (targetSelection_.isDynamic() && isTwoDiffGroups_)
        {
            // Print Selection Indices
            fprintf(fp, "[ %s_frame_%i ]", targetSelection_.name(), frnr);
            for (std::size_t i = 0; i < targetSelection_.atomIndices().size(); ++i)
            {
                fprintf(fp, "%c %4i", (i % 15 == 0 ? '\n' : ' '), targetSelection_.atomIndices()[i] + 1);
            }
            // Print Donors Indices
            fprintf(fp, "\n[ donors_hydrogens_%s_frame_%i ]", targetSelection_.name(), frnr);
            for (const auto& donor : targetInfo_.donors)
            {
                fprintf(fp, "\n %4i", donor.ai + 1);
                for (std::size_t j = 0; j < donor.h_atoms.size(); ++j)
                {
                    fprintf(fp, " %4i", donor.h_atoms[j] + 1);
                }
            }
            // Print Acceptors Indices
            fprintf(fp, "\n[ acceptors_%s_frame_%i ]", targetSelection_.name(), frnr);
            for (std::size_t i = 0; i < targetInfo_.acceptors.size(); ++i)
            {
                fprintf(fp, "%c %4i", (i % 15 == 0 ? '\n' : ' '), targetInfo_.acceptors[i].ai + 1);
            }
            fprintf(fp, "\n");
        }
        gmx_ffclose(fp);
    }

    AnalysisDataHandle              dhHbnum       = pdata->dataHandle(hbnum_);
    AnalysisDataHandle              dhDist        = pdata->dataHandle(distances_);
    AnalysisDataHandle              dhAng         = pdata->dataHandle(angles_);
    AnalysisDataHandle              dhDan         = pdata->dataHandle(dan_);
    std::size_t                     distIterator  = 0;
    std::size_t                     angleIterator = 0;
    std::unordered_set<std::size_t> donors;
    std::unordered_set<std::size_t> acceptors;


    if (!fnmHbnumOut_.empty())
    {
        dhHbnum.startFrame(frnr, fr.time);
    }
    if (!fnmHbdistOut_.empty())
    {
        dhDist.startFrame(frnr, fr.time);
    }
    if (!fnmHbangOut_.empty())
    {
        dhAng.startFrame(frnr, fr.time);
    }
    if (!fnmHbdanOut_.empty())
    {
        dhDan.startFrame(frnr, fr.time);
    }


    std::vector<HBond>   daMap;
    AnalysisNeighborhood nb;
    nb.setCutoff(nbsearchCutoff_);
    const t_info* infoTool1 = &refInfo_;
    const t_info* infoTool2;
    if (isTwoDiffGroups_)
    {
        infoTool2 = &targetInfo_;
    }
    else
    {
        infoTool2 = &refInfo_;
    }
    for (std::size_t cycleIterator = 0; cycleIterator < 2; ++cycleIterator)
    {
        std::vector<gmx::RVec> positionsDonor, positionsAcceptor;
        for (const auto& acceptor : infoTool1->acceptors)
        {
            positionsAcceptor.emplace_back(fr.x[acceptor.ai]);
        }
        for (const auto& donor : infoTool2->donors)
        {
            positionsDonor.emplace_back(fr.x[donor.ai]);
        }
        AnalysisNeighborhoodPositions       nbPos(positionsAcceptor);
        gmx::AnalysisNeighborhoodSearch     start      = nb.initSearch(pbc, nbPos);
        gmx::AnalysisNeighborhoodPairSearch pairSearch = start.startPairSearch(positionsDonor);
        gmx::AnalysisNeighborhoodPair       pair;
        while (pairSearch.findNextPair(&pair))
        {
            if (infoTool1->acceptors[pair.refIndex()].ai == infoTool2->donors[pair.testIndex()].ai)
            {
                continue;
            }
            t_acceptor acceptor = infoTool1->acceptors[pair.refIndex()];
            t_donor    donor    = infoTool2->donors[pair.testIndex()];
            gmx::RVec  vectorDA = { 0, 0, 0 };
            pbc_dx(pbc, fr.x[acceptor.ai], fr.x[donor.ai], vectorDA.as_vec());
            if (vectorDA.norm() > hbDistCutoff_)
            {
                continue;
            }
            else
            {
                for (const auto hIndex : donor.h_atoms)
                {
                    gmx::RVec vectorDH = { 0, 0, 0 };
                    float     degree   = 0;
                    pbc_dx(pbc, fr.x[hIndex], fr.x[donor.ai], vectorDH.as_vec());
                    degree = gmx_angle(vectorDA, vectorDH) * gmx::c_rad2Deg;
                    if (degree <= hbAngleCutoff_)
                    {
                        if (!isTwoDiffGroups_ && (acceptor.isAlsoDonor && acceptor.ai < donor.ai)
                            && mergeHydrogens_)
                        {
                            daMap.emplace_back(acceptor.ai, donor.ai, true);
                        }
                        else if (mergeHydrogens_)
                        {
                            daMap.emplace_back(donor.ai, acceptor.ai, acceptor.isAlsoDonor);
                        }
                        else
                        {
                            daMap.emplace_back(donor.ai, acceptor.ai, acceptor.isAlsoDonor, hIndex);
                        }

                        if (!fnmHbdistOut_.empty())
                        {
                            dhDist.setPoint(distIterator++, vectorDA.norm());
                        }

                        if (!fnmHbangOut_.empty())
                        {
                            dhAng.setPoint(angleIterator++, degree);
                        }
                        if (!fnmHbdanOut_.empty())
                        {
                            if (!isTwoDiffGroups_ && (acceptor.isAlsoDonor && acceptor.ai < donor.ai)
                                && mergeHydrogens_)
                            {
                                donors.insert(acceptor.ai);
                                acceptors.insert(donor.ai);
                            }
                            else
                            {
                                donors.insert(donor.ai);
                                acceptors.insert(acceptor.ai);
                            }
                        }
                    }
                }
            }
        }
        if (isTwoDiffGroups_)
        {
            std::swap(infoTool1, infoTool2);
        }
        else
        {
            break;
        }
    }
    storage_.addData(frnr, daMap);

    if (!fnmHbnumOut_.empty())
    {

        dhHbnum.setPoint(0, prepareFrameData(daMap).size());
        dhHbnum.finishFrame();
    }
    if (!fnmHbdistOut_.empty())
    {
        dhDist.finishFrame();
    }
    if (!fnmHbangOut_.empty())
    {
        dhAng.finishFrame();
    }
    if (!fnmHbdanOut_.empty())
    {
        dhDan.setPoint(0, donors.size());
        dhDan.setPoint(1, acceptors.size());
        dhDan.finishFrame();
    }
}


void Hbond::finishAnalysis(int /*nframes*/)
{
    if (!fnmHbdistOut_.empty())
    {
        AbstractAverageHistogram& averageHistogramDist = histogramModuleDist_->averager();
        averageHistogramDist.normalizeProbability();
        averageHistogramDist.done();
    }
    if (!fnmHbangOut_.empty())
    {
        AbstractAverageHistogram& averageHistogramAng = histogramModuleAng_->averager();
        averageHistogramAng.normalizeProbability();
        averageHistogramAng.done();
    }
    please_cite(stdout, "Gorelov2024b");
}

std::vector<HBond> Hbond::prepareFrameData(const std::vector<HbondStorageFrame>& data) const
{
    std::vector<HBond> temp;
    for (const auto& i : data)
    {
        for (const auto& j : i.hbondData_)
        {
            if (!isTwoDiffGroups_ && mergeHydrogens_ && j.acceptorIsAlsoDonor && j.acceptor < j.donor)
            {
                temp.emplace_back(!j);
            }
            else
            {
                temp.emplace_back(j);
            }
        }
    }
    std::sort(temp.begin(), temp.end());
    temp.erase(std::unique(temp.begin(), temp.end()), temp.end());
    if (isTwoDiffGroups_ && mergeHydrogens_)
    {
        printf("Merging hbonds with Acceptor and Donor swapped\n");
        for (std::vector<HBond>::iterator i = temp.begin(); i != temp.end();)
        {
            auto find = std::equal_range(i, temp.end(), !(*i));
            if (find.first != find.second)
            {
                i = temp.erase(find.first, find.second);
            }
            else
            {
                ++i;
            }
        }
    }
    return temp;
}

std::vector<HBond> Hbond::prepareFrameData(const std::vector<HBond>& data) const
{
    std::vector<HBond> temp;
    for (const auto& i : data)
    {
        if (!isTwoDiffGroups_ && mergeHydrogens_ && i.acceptorIsAlsoDonor && i.acceptor < i.donor)
        {
            temp.emplace_back(!i);
        }
        else
        {
            temp.emplace_back(i);
        }
    }
    std::sort(temp.begin(), temp.end());
    temp.erase(std::unique(temp.begin(), temp.end()), temp.end());
    if (isTwoDiffGroups_ && mergeHydrogens_)
    {
        for (std::vector<HBond>::iterator i = temp.begin(); i != temp.end();)
        {
            auto find = std::equal_range(i, temp.end(), !(*i));
            if (find.first != find.second)
            {
                i = temp.erase(find.first, find.second);
            }
            else
            {
                ++i;
            }
        }
    }
    return temp;
}

void Hbond::writeOutput()
{
    std::vector<HbondStorageFrame> dataOut;
    FILE*                          fp;
    fp      = gmx_ffopen(fnmHbondOut_, "a");
    dataOut = storage_.getData();
    // Print Hbond Indices
    std::string selects;
    if (isTwoDiffGroups_)
    {
        selects = std::string(refInfo_.selectionPtr->name()) + "-"
                  + std::string(targetInfo_.selectionPtr->name());
    }
    else
    {
        selects = refInfo_.selectionPtr->name();
    }
    if (perFrame_)
    {
        if (isTwoDiffGroups_ && mergeHydrogens_)
        {
            printf("Merging hbonds with Acceptor and Donor swapped\n");
        }
        for (const auto& i : dataOut)
        {
            std::vector<HBond> tempFrameData = prepareFrameData(i.hbondData_);
            fprintf(fp, "[ hbonds_%s_frame_%i ]", selects.c_str(), i.frameNumber_);
            for (const auto& j : tempFrameData)
            {
                if (mergeHydrogens_)
                {
                    fprintf(fp, "\n %6i %6i", j.donor + 1, j.acceptor + 1);
                }
                else
                {
                    fprintf(fp, "\n %6i %6i %6i", j.donor + 1, j.hydrogen + 1, j.acceptor + 1);
                }
            }
            fprintf(fp, "\n");
        }
    }
    else
    {
        std::vector<HBond> mergedFrameData = prepareFrameData(dataOut);
        fprintf(fp, "[ hbonds_%s ]", selects.c_str());
        if (mergeHydrogens_)
        {
            for (const auto& i : mergedFrameData)
            {
                fprintf(fp, "\n %6i %6i", i.donor + 1, i.acceptor + 1);
            }
        }
        else
        {
            for (const auto& i : mergedFrameData)
            {
                fprintf(fp, "\n %6i %6i %6i", i.donor + 1, i.hydrogen + 1, i.acceptor + 1);
            }
        }
    }
    fprintf(fp, "\n");

    gmx_ffclose(fp);
}

} // namespace

const char HbondInfo::name[]             = "hbond";
const char HbondInfo::shortDescription[] = "Compute and analyze hydrogen bonds.";

TrajectoryAnalysisModulePointer HbondInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Hbond);
}

} // namespace analysismodules

} // namespace gmx
