/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
 * Implements gmx::analysismodules::Distance.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "distance.h"

#include <cstdio>

#include <memory>
#include <string>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class TopologyInformation;

namespace analysismodules
{

namespace
{

class Distance : public TrajectoryAnalysisModule
{
public:
    Distance();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;

    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;

    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    SelectionList sel_;
    std::string   fnAverage_;
    std::string   fnAll_;
    std::string   fnXYZ_;
    std::string   fnHistogram_;
    std::string   fnAllStats_;
    double        meanLength_;
    double        lengthDev_;
    double        binWidth_;

    AnalysisData                             distances_;
    AnalysisData                             xyz_;
    AnalysisDataAverageModulePointer         summaryStatsModule_;
    AnalysisDataAverageModulePointer         allStatsModule_;
    AnalysisDataFrameAverageModulePointer    averageModule_;
    AnalysisDataSimpleHistogramModulePointer histogramModule_;

    // Copy and assign disallowed by base.
};

Distance::Distance() :
    meanLength_(0.1),
    lengthDev_(1.0),
    binWidth_(0.001),
    summaryStatsModule_(std::make_unique<AnalysisDataAverageModule>()),
    allStatsModule_(std::make_unique<AnalysisDataAverageModule>()),
    averageModule_(std::make_unique<AnalysisDataFrameAverageModule>()),
    histogramModule_(std::make_unique<AnalysisDataSimpleHistogramModule>())
{
    summaryStatsModule_->setAverageDataSets(true);
    distances_.addModule(summaryStatsModule_);
    distances_.addModule(allStatsModule_);
    distances_.addModule(averageModule_);
    distances_.addModule(histogramModule_);

    registerAnalysisDataset(&distances_, "dist");
    registerAnalysisDataset(&xyz_, "xyz");
    registerBasicDataset(summaryStatsModule_.get(), "stats");
    registerBasicDataset(allStatsModule_.get(), "allstats");
    registerBasicDataset(averageModule_.get(), "average");
    registerBasicDataset(&histogramModule_->averager(), "histogram");
}


void Distance::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] calculates distances between pairs of positions",
        "as a function of time. Each selection specifies an independent set",
        "of distances to calculate. Each selection should consist of pairs",
        "of positions, and the distances are computed between positions 1-2,",
        "3-4, etc.[PAR]",
        "[TT]-oav[tt] writes the average distance as a function of time for",
        "each selection.",
        "[TT]-oall[tt] writes all the individual distances.",
        "[TT]-oxyz[tt] does the same, but the x, y, and z components of the",
        "distance are written instead of the norm.",
        "[TT]-oh[tt] writes a histogram of the distances for each selection.",
        "The location of the histogram is set with [TT]-len[tt] and",
        "[TT]-tol[tt]. Bin width is set with [TT]-binw[tt].",
        "[TT]-oallstat[tt] writes out the average and standard deviation for",
        "each individual distance, calculated over the frames.[PAR]",
        "Note that [THISMODULE] calculates distances between fixed pairs",
        "(1-2, 3-4, etc.) within a single selection.  To calculate distances",
        "between two selections, including minimum, maximum, and pairwise",
        "distances, use [gmx-pairdist]."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("oav")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnAverage_)
                               .defaultBasename("distave")
                               .description("Average distances as function of time"));
    options->addOption(FileNameOption("oall")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnAll_)
                               .defaultBasename("dist")
                               .description("All distances as function of time"));
    options->addOption(FileNameOption("oxyz")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnXYZ_)
                               .defaultBasename("distxyz")
                               .description("Distance components as function of time"));
    options->addOption(FileNameOption("oh")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnHistogram_)
                               .defaultBasename("disthist")
                               .description("Histogram of the distances"));
    options->addOption(FileNameOption("oallstat")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnAllStats_)
                               .defaultBasename("diststat")
                               .description("Statistics for individual distances"));
    options->addOption(
            SelectionOption("select").storeVector(&sel_).required().dynamicMask().multiValue().description(
                    "Position pairs to calculate distances for"));
    // TODO: Extend the histogramming implementation to allow automatic
    // extension of the histograms to cover the data, removing the need for
    // the first two options.
    options->addOption(
            DoubleOption("len").store(&meanLength_).description("Mean distance for histogramming"));
    options->addOption(
            DoubleOption("tol").store(&lengthDev_).description("Width of full distribution as fraction of [TT]-len[tt]"));
    options->addOption(
            DoubleOption("binw").store(&binWidth_).description("Bin width for histogramming"));
}


/*! \brief
 * Checks that selections conform to the expectations of the tool.
 */
void checkSelections(const SelectionList& sel)
{
    for (size_t g = 0; g < sel.size(); ++g)
    {
        if (sel[g].posCount() % 2 != 0)
        {
            std::string message = formatString(
                    "Selection '%s' does not evaluate into an even number of positions "
                    "(there are %d positions)",
                    sel[g].name(),
                    sel[g].posCount());
            GMX_THROW(InconsistentInputError(message));
        }
        if (sel[g].isDynamic())
        {
            for (int i = 0; i < sel[g].posCount(); i += 2)
            {
                if (sel[g].position(i).selected() != sel[g].position(i + 1).selected())
                {
                    std::string message = formatString(
                            "Dynamic selection %d does not select "
                            "a consistent set of pairs over the frames",
                            static_cast<int>(g + 1));
                    GMX_THROW(InconsistentInputError(message));
                }
            }
        }
    }
}


void Distance::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& /*top*/)
{
    checkSelections(sel_);

    distances_.setDataSetCount(sel_.size());
    xyz_.setDataSetCount(sel_.size());
    for (size_t i = 0; i < sel_.size(); ++i)
    {
        const int distCount = sel_[i].posCount() / 2;
        distances_.setColumnCount(i, distCount);
        xyz_.setColumnCount(i, distCount * 3);
    }
    const double histogramMin = (1.0 - lengthDev_) * meanLength_;
    const double histogramMax = (1.0 + lengthDev_) * meanLength_;
    histogramModule_->init(histogramFromRange(histogramMin, histogramMax).binWidth(binWidth_).includeAll());

    if (!fnAverage_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnAverage_);
        plotm->setTitle("Average distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plotm->appendLegend(sel_[g].name());
        }
        averageModule_->addModule(plotm);
    }

    if (!fnAll_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnAll_);
        plotm->setTitle("Distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        // TODO: Add legends? (there can be a massive amount of columns)
        distances_.addModule(plotm);
    }

    if (!fnXYZ_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnXYZ_);
        plotm->setTitle("Distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        // TODO: Add legends? (there can be a massive amount of columns)
        xyz_.addModule(plotm);
    }

    if (!fnHistogram_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnHistogram_);
        plotm->setTitle("Distance histogram");
        plotm->setXLabel("Distance (nm)");
        plotm->setYLabel("Probability");
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plotm->appendLegend(sel_[g].name());
        }
        histogramModule_->averager().addModule(plotm);
    }

    if (!fnAllStats_.empty())
    {
        AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnAllStats_);
        plotm->setErrorsAsSeparateColumn(true);
        plotm->setTitle("Statistics for individual distances");
        plotm->setXLabel("Distance index");
        plotm->setYLabel("Average/standard deviation (nm)");
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plotm->appendLegend(std::string(sel_[g].name()) + " avg");
            plotm->appendLegend(std::string(sel_[g].name()) + " std.dev.");
        }
        // TODO: Consider whether this output format is the best possible.
        allStatsModule_->addModule(plotm);
    }
}


void Distance::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata)
{
    AnalysisDataHandle   distHandle = pdata->dataHandle(distances_);
    AnalysisDataHandle   xyzHandle  = pdata->dataHandle(xyz_);
    const SelectionList& sel        = TrajectoryAnalysisModuleData::parallelSelections(sel_);

    checkSelections(sel);

    distHandle.startFrame(frnr, fr.time);
    xyzHandle.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        distHandle.selectDataSet(g);
        xyzHandle.selectDataSet(g);
        for (int i = 0, n = 0; i < sel[g].posCount(); i += 2, ++n)
        {
            const SelectionPosition& p1 = sel[g].position(i);
            const SelectionPosition& p2 = sel[g].position(i + 1);
            rvec                     dx;
            if (pbc != nullptr)
            {
                pbc_dx(pbc, p2.x(), p1.x(), dx);
            }
            else
            {
                rvec_sub(p2.x(), p1.x(), dx);
            }
            real dist     = norm(dx);
            bool bPresent = p1.selected() && p2.selected();
            distHandle.setPoint(n, dist, bPresent);
            xyzHandle.setPoints(n * 3, 3, dx, bPresent);
        }
    }
    distHandle.finishFrame();
    xyzHandle.finishFrame();
}


void Distance::finishAnalysis(int /*nframes*/)
{
    AbstractAverageHistogram& averageHistogram = histogramModule_->averager();
    averageHistogram.normalizeProbability();
    averageHistogram.done();
}


void Distance::writeOutput()
{
    SelectionList::const_iterator sel;
    int                           index = 0;
    for (sel = sel_.begin(); sel != sel_.end(); ++sel, ++index)
    {
        printf("%s:\n", sel->name());
        printf("  Number of samples:  %d\n", summaryStatsModule_->sampleCount(index, 0));
        printf("  Average distance:   %-8.5f nm\n", summaryStatsModule_->average(index, 0));
        printf("  Standard deviation: %-8.5f nm\n", summaryStatsModule_->standardDeviation(index, 0));
    }
}

} // namespace

const char DistanceInfo::name[]             = "distance";
const char DistanceInfo::shortDescription[] = "Calculate distances between pairs of positions";

TrajectoryAnalysisModulePointer DistanceInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Distance);
}

} // namespace analysismodules

} // namespace gmx
