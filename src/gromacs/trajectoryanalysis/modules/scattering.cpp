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
 * Implements gmx::analysismodules::scattering.
 *
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */

#include "gmxpre.h"

#include "scattering.h"

#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/pleasecite.h"

#include "scattering-debye-sans.h"
#include "scattering-debye-saxs.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

//! Type of scattering calculation to perform
enum class ScatterType
{
    SAXS,
    SANS,
    Count
};

//! String values corresponding to ScatterType
const gmx::EnumerationArray<ScatterType, const char*> ScatterTypeNames = { { "saxs", "sans" } };

class Scattering : public TrajectoryAnalysisModule
{
public:
    Scattering();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;
    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    SelectionList                                        sel_;
    std::string                                          fnScattering_;
    double                                               startQ_;
    double                                               endQ_;
    double                                               qSpacing_;
    double                                               monteCarloCoverage_;
    double                                               binwidth_;
    int                                                  seed_;
    bool                                                 bNormalize_;
    bool                                                 bMonteCarloScattering_;
    bool                                                 bSeedSet_;
    bool                                                 bMonteCarloCoverageSet_;
    ScatterType                                          scatterEnum_;
    AnalysisData                                         intensity_;
    AnalysisDataAverageModulePointer                     averageModule_;
    std::vector<double>                                  qList_;
    std::vector<double>                                  scatteringAtQ0_;
    std::vector<std::unique_ptr<ComputeDebyeScattering>> computeScattering_;
};

Scattering::Scattering() :
    startQ_(0.0),
    endQ_(2.0),
    qSpacing_(0.01),
    monteCarloCoverage_(0.2),
    binwidth_(0.1),
    seed_(2023),
    bNormalize_(true),
    bMonteCarloScattering_(false),
    bSeedSet_(false),
    bMonteCarloCoverageSet_(false),
    scatterEnum_(ScatterType::SANS),
    averageModule_(std::make_unique<AnalysisDataAverageModule>())
{
    intensity_.addModule(averageModule_);
    intensity_.setMultipoint(true);
    registerAnalysisDataset(&intensity_, "scattering");
}

void Scattering::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] calculates SANS and SAXS scattering curves using Debye method.[PAR]",
        "The scattering intensity, I(q), as a function of scattering angle q",
        "with averaging over frames. [PAR]",
        "[PAR]",
        "Note that this is a new implementation of the SANS/SAXS utilities added in",
        "GROMACS 2024. If you need the old ones,",
        "use [TT]gmx sans-legacy[tt] or [TT]gmx saxs-legacy[tt]."
    };

    settings->setHelpText(desc);
    options->addOption(FileNameOption("o")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnScattering_)
                               .defaultBasename("scattering")
                               .description("scattering intensity as a function of q"));
    options->addOption(SelectionOption("sel").storeVector(&sel_).required().multiValue().description(
            "Selection for Scattering calculation"));
    options->addOption(DoubleOption("startq").store(&startQ_).defaultValue(0.0).description(
            "smallest q value (1/nm)"));
    options->addOption(DoubleOption("endq").store(&endQ_).defaultValue(2).description(
            "largest q value (1/nm)"));
    options->addOption(
            DoubleOption("qspacing").store(&qSpacing_).defaultValue(0.01).description("spacing of q values (1/nm)"));
    options->addOption(
            DoubleOption("binwidth").store(&binwidth_).defaultValue(0.1).description("Bin width (nm) for P(r)"));
    options->addOption(DoubleOption("mc-coverage")
                               .store(&monteCarloCoverage_)
                               .storeIsSet(&bMonteCarloCoverageSet_)
                               .defaultValue(0.2)
                               .description("coverage of Monte Carlo (%)"));
    options->addOption(
            IntegerOption("seed").store(&seed_).storeIsSet(&bSeedSet_).defaultValue(2023).description("random seed for Monte Carlo"));
    options->addOption(
            BooleanOption("norm").store(&bNormalize_).defaultValue(false).description("normalize scattering intensities"));
    options->addOption(BooleanOption("mc")
                               .store(&bMonteCarloScattering_)
                               .defaultValue(true)
                               .description("use Monte Carlo to scattering intensities"));
    options->addOption(EnumOption<ScatterType>("scattering-type")
                               .enumValue(ScatterTypeNames)
                               .store(&scatterEnum_)
                               .defaultValue(ScatterType::SANS)
                               .description("Scattering type"));
}

void Scattering::optionsFinished(TrajectoryAnalysisSettings* /*settings*/)
{
    if (!bMonteCarloScattering_ && (bMonteCarloCoverageSet_ || bSeedSet_))
    {
        GMX_THROW(InvalidInputError("You cannot set seed or coverage unless you specify -mc"));
    }
    if (monteCarloCoverage_ <= 0 || monteCarloCoverage_ > 1)
    {
        GMX_THROW(InvalidInputError("You must specify coverage in (0, 1]"));
    }

    if (startQ_ < 0.0)
    {
        GMX_THROW(InconsistentInputError("startq cannot be < 0"));
    }

    if (endQ_ < startQ_)
    {
        GMX_THROW(InconsistentInputError("endq cannot be < startq"));
    }

    if (qSpacing_ <= 0)
    {
        GMX_THROW(InconsistentInputError("qspacing cannot be <= 0"));
    }

    // actually set the q values to be used
    double qVal = startQ_;
    while (qVal < endQ_)
    {
        qList_.push_back(qVal);
        qVal += qSpacing_;
    }
}

void Scattering::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top)
{
    const auto*          atoms    = top.atoms();
    std::vector<Isotope> isotopes = getIsotopes(atoms);
    intensity_.setDataSetCount(sel_.size());
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        intensity_.setColumnCount(g, qList_.size());
        switch (scatterEnum_)
        {
            case ScatterType::SANS:
                computeScattering_.push_back(std::make_unique<SansDebye>(isotopes));
                break;
            case ScatterType::SAXS:
                computeScattering_.push_back(std::make_unique<SaxsDebye>(isotopes, qList_));
                break;
            default: GMX_THROW(InconsistentInputError("Unknown scattering type"));
        }
        computeScattering_[g]->setBinWidth(binwidth_);
        computeScattering_[g]->addQList(qList_);
    }

    if (!fnScattering_.empty())
    {
        averageModule_->setXAxis(startQ_, qSpacing_);
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnScattering_);
        plotm->setXFormat(2, 8);
        plotm->setYFormat(1, 8);
        plotm->setTitle("Scattering intensity: I(q)");
        plotm->setXLabel("q (1/nm)");
        plotm->setYLabel("Intensity");
        averageModule_->addModule(plotm);
    }
}

void Scattering::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata)
{
    AnalysisDataHandle   scatterHandle = pdata->dataHandle(intensity_);
    const SelectionList& sel = gmx::TrajectoryAnalysisModuleData::parallelSelections(sel_);
    scatterHandle.startFrame(frnr, fr.time);
    matrix fBox;
    copy_mat(fr.box, fBox);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        scatterHandle.selectDataSet(g);
        computeScattering_[g]->getMaxDist(fBox);
        computeScattering_[g]->initPairDistHist();
        if (bMonteCarloScattering_)
        {
            computeScattering_[g]->computeMonteCarloPairDistancesHistogram(
                    pbc, sel[g], monteCarloCoverage_, seed_);
        }
        else
        {
            computeScattering_[g]->computeDirectPairDistancesHistogram(pbc, sel[g]);
        }
        if (scatteringAtQ0_.empty())
        {
            scatteringAtQ0_.push_back(computeScattering_[g]->computeIntensityZeroQ());
        }
        computeScattering_[g]->computeIntensity();
        for (size_t qi = 0; qi != qList_.size(); ++qi)
        {
            double DebyeScattering = computeScattering_[g]->getIntensity(qi);

            if (bNormalize_)
            {
                DebyeScattering /= scatteringAtQ0_[g];
            }
            scatterHandle.setPoint(qi, DebyeScattering, true);
            scatterHandle.finishPointSet();
        }
        computeScattering_[g]->clearHist();
    }
    scatterHandle.finishFrame();
}


void Scattering::finishAnalysis(int /* nframes */)
{
    please_cite(stdout, "Shvetsov2013");
    if (scatterEnum_ == ScatterType::SAXS)
    {
        please_cite(stdout, "Cromer1968a");
    }
}


void Scattering::writeOutput() {}

} // namespace

const char ScatteringInfo::name[] = "scattering";
const char ScatteringInfo::shortDescription[] =
        "Calculate small angle scattering profiles for SANS or SAXS";

TrajectoryAnalysisModulePointer ScatteringInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Scattering);
}

} // namespace analysismodules

} // namespace gmx
