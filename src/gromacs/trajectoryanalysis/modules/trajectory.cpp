/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
 * Implements gmx::analysismodules::Trajectory.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "trajectory.h"

#include <cstddef>

#include <algorithm>
#include <array>
#include <iterator>
#include <memory>
#include <string>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/optionfiletype.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"

struct t_pbc;

namespace gmx
{
class TopologyInformation;

namespace analysismodules
{

namespace
{

/********************************************************************
 * Trajectory
 */

class Trajectory : public TrajectoryAnalysisModule
{
public:
    Trajectory();

    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void optionsFinished(TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;

    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;

    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    SelectionList sel_;

    std::string         fnX_;
    std::string         fnV_;
    std::string         fnF_;
    std::array<bool, 4> dimMask_;
    std::array<bool, 4> maskSet_;

    AnalysisData xdata_;
    AnalysisData vdata_;
    AnalysisData fdata_;
};

Trajectory::Trajectory() : dimMask_{ true, true, true, false }, maskSet_{}
{
    registerAnalysisDataset(&xdata_, "x");
    registerAnalysisDataset(&vdata_, "v");
    registerAnalysisDataset(&fdata_, "f");
}


void Trajectory::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] plots coordinates, velocities, and/or forces for",
        "provided selections. By default, the X, Y, and Z components for",
        "the requested vectors are plotted, but specifying one or more of",
        "[TT]-len[tt], [TT]-x[tt], [TT]-y[tt], and [TT]-z[tt] overrides this.",
        "",
        "For dynamic selections, currently the values are written out for",
        "all positions that the selection could select."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("ox")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnX_)
                               .defaultBasename("coord")
                               .description("Coordinates for each position as a function of time"));
    options->addOption(FileNameOption("ov")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnV_)
                               .defaultBasename("veloc")
                               .description("Velocities for each position as a function of time"));
    options->addOption(FileNameOption("of")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnF_)
                               .defaultBasename("force")
                               .description("Forces for each position as a function of time"));

    options->addOption(
            SelectionOption("select").storeVector(&sel_).required().dynamicMask().multiValue().description(
                    "Selections to analyze"));

    options->addOption(
            BooleanOption("x").store(&dimMask_[XX]).storeIsSet(&maskSet_[XX]).description("Plot X component"));
    options->addOption(
            BooleanOption("y").store(&dimMask_[YY]).storeIsSet(&maskSet_[YY]).description("Plot Y component"));
    options->addOption(
            BooleanOption("z").store(&dimMask_[ZZ]).storeIsSet(&maskSet_[ZZ]).description("Plot Z component"));
    options->addOption(
            BooleanOption("len").store(&dimMask_[DIM]).storeIsSet(&maskSet_[DIM]).description("Plot vector length"));
}


void Trajectory::optionsFinished(TrajectoryAnalysisSettings* settings)
{
    int frameFlags = TRX_NEED_X;
    if (!fnV_.empty())
    {
        frameFlags |= TRX_READ_V;
    }
    if (!fnF_.empty())
    {
        frameFlags |= TRX_READ_F;
    }
    settings->setFrameFlags(frameFlags);
    if (std::count(std::begin(maskSet_), std::end(maskSet_), true) > 0)
    {
        for (int i = 0; i <= DIM; ++i)
        {
            if (!maskSet_[i])
            {
                dimMask_[i] = false;
            }
        }
    }
}


void Trajectory::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& /*top*/)
{
    if (!fnX_.empty())
    {
        xdata_.setDataSetCount(sel_.size());
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            xdata_.setColumnCount(g, 3 * sel_[g].posCount());
        }
        AnalysisDataVectorPlotModulePointer plot(new AnalysisDataVectorPlotModule(settings.plotSettings()));
        plot->setWriteMask(dimMask_.data());
        plot->setFileName(fnX_);
        plot->setTitle("Coordinates");
        plot->setXAxisIsTime();
        plot->setYLabel("Value [nm]");
        xdata_.addModule(plot);
    }
    if (!fnV_.empty())
    {
        vdata_.setDataSetCount(sel_.size());
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            sel_[g].setEvaluateVelocities(true);
            vdata_.setColumnCount(g, 3 * sel_[g].posCount());
        }
        AnalysisDataVectorPlotModulePointer plot(new AnalysisDataVectorPlotModule(settings.plotSettings()));
        plot->setWriteMask(dimMask_.data());
        plot->setFileName(fnV_);
        plot->setTitle("Velocities");
        plot->setXAxisIsTime();
        plot->setYLabel("Value [nm/ps]");
        vdata_.addModule(plot);
    }
    if (!fnF_.empty())
    {
        fdata_.setDataSetCount(sel_.size());
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            sel_[g].setEvaluateForces(true);
            fdata_.setColumnCount(g, 3 * sel_[g].posCount());
        }
        AnalysisDataVectorPlotModulePointer plot(new AnalysisDataVectorPlotModule(settings.plotSettings()));
        plot->setWriteMask(dimMask_.data());
        plot->setFileName(fnF_);
        plot->setTitle("Forces");
        plot->setXAxisIsTime();
        plot->setYLabel("Value [kJ mol\\S-1\\N nm\\S-1\\N]");
        fdata_.addModule(plot);
    }
}


//! Helper function for Trajectory::analyzeFrame
template<typename T>
void analyzeFrameImpl(int frnr, const t_trxframe& fr, AnalysisDataHandle* dh, const SelectionList& sel, T getField)
{
    if (dh->isValid())
    {
        dh->startFrame(frnr, fr.time);
        for (size_t g = 0; g < sel.size(); ++g)
        {
            dh->selectDataSet(g);
            for (int i = 0; i < sel[g].posCount(); ++i)
            {
                const SelectionPosition& pos = sel[g].position(i);
                dh->setPoints(i * 3, 3, getField(pos), pos.selected());
            }
        }
        dh->finishFrame();
    }
}

void Trajectory::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* /* pbc */, TrajectoryAnalysisModuleData* pdata)
{
    AnalysisDataHandle   dh  = pdata->dataHandle(xdata_);
    const SelectionList& sel = TrajectoryAnalysisModuleData::parallelSelections(sel_);
    analyzeFrameImpl(frnr, fr, &dh, sel, [](const SelectionPosition& pos) { return pos.x(); });
    if (fr.bV)
    {
        analyzeFrameImpl(frnr, fr, &dh, sel, [](const SelectionPosition& pos) { return pos.v(); });
    }
    if (fr.bF)
    {
        analyzeFrameImpl(frnr, fr, &dh, sel, [](const SelectionPosition& pos) { return pos.f(); });
    }
}


void Trajectory::finishAnalysis(int /*nframes*/) {}


void Trajectory::writeOutput() {}

} // namespace

const char TrajectoryInfo::name[] = "trajectory";
const char TrajectoryInfo::shortDescription[] =
        "Print coordinates, velocities, and/or forces for selections";

TrajectoryAnalysisModulePointer TrajectoryInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Trajectory);
}

} // namespace analysismodules

} // namespace gmx
