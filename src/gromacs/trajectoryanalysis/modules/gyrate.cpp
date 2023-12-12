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
 * Implements gmx::analysismodules::Gyrate.
 *
 * \author Vladimir Basov <vovabasov830@gmail.com>
 * \author Alexey Shvetsov <alexxyum@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gyrate.h"

#include <algorithm>
#include <memory>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

//! Enum of different ways to use the weight of atoms.
enum class GyrateMode : std::size_t
{
    Mass = 0,
    Charge,
    Geometry,
    Count
};

//! String values corresponding to weight-assignment modes.
const gmx::EnumerationArray<GyrateMode, const char*> c_GyrateModeNames = {
    { "mass", "charge", "geometry" }
};

class Gyrate : public TrajectoryAnalysisModule
{
public:
    Gyrate();
    void initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings) override;
    void initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& top) override;
    void analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata) override;
    void finishAnalysis(int nframes) override;
    void writeOutput() override;

private:
    //! String value that defines plot output filename. Set in initial options.
    std::string fnGyrate_;
    //! Selections for Gyrate output. Set in initial options
    Selection sel_;
    //! Container in which the gyrate data is written
    AnalysisData gyrate_;
    //! Enum value for creating weight mode. Set in initial options.
    GyrateMode gMode_;
    //! Function that receives the value of the weights of atoms(mass, charge).
    static real getWeighFactor(SelectionPosition position, GyrateMode mode);
};

Gyrate::Gyrate() : gMode_(GyrateMode::Mass)
{
    registerAnalysisDataset(&gyrate_, "gyrate");
}

real Gyrate::getWeighFactor(const SelectionPosition position, const GyrateMode mode)
{
    switch (mode)
    {
        case GyrateMode::Geometry: return 1.;
        case GyrateMode::Mass: return position.mass();
        case GyrateMode::Charge: return position.charge();
        default: GMX_RELEASE_ASSERT(false, "Invalid value of GyrateMode"); return 0.;
    }
}

void Gyrate::initOptions(IOptionsContainer* options, TrajectoryAnalysisSettings* settings)
{
    static const char* const desc[] = {
        "[THISMODULE] computes the radius of gyration of a molecule",
        "and the radii of gyration about the [IT]x[it]-, [IT]y[it]- and [IT]z[it]-axes,",
        "as a function of time. The atoms are explicitly mass weighted.[PAR]",
        "The axis components corresponds to the mass-weighted root-mean-square",
        "of the radii components orthogonal to each axis, for example:[PAR]",
        "Rg(x) = sqrt((sum_i w_i (R_i(y)^2 + R_i(z)^2))/(sum_i w_i)).[PAR]",
        "where w_i is the weight value in the given situation (mass, charge, unit)",
        "[PAR]",
        "Note that this is a new implementation of the gyrate utility added in",
        "GROMACS 2024. If you need the old one, use [TT]gmx gyrate-legacy[tt]."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("o")
                               .filetype(OptionFileType::Plot)
                               .outputFile()
                               .store(&fnGyrate_)
                               .defaultBasename("gyrate-taf")
                               .required(true)
                               .description("Filename for gyrate plot output"));
    options->addOption(SelectionOption("sel").store(&sel_).required().dynamicMask().description(
            "Select group to compute gyrate radius"));
    options->addOption(EnumOption<GyrateMode>("mode")
                               .store(&gMode_)
                               .defaultValue(GyrateMode::Mass)
                               .enumValue(c_GyrateModeNames)
                               .description("Atom weighting mode"));
    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efNoUserPBC);
    settings->setPBC(true);
}

void Gyrate::initAnalysis(const TrajectoryAnalysisSettings& settings, const TopologyInformation& /*top*/)
{
    gyrate_.setColumnCount(0, 4);
    AnalysisDataPlotModulePointer plotm(new AnalysisDataPlotModule(settings.plotSettings()));
    plotm->setFileName(fnGyrate_);
    plotm->setTitle("Radius of gyration (total and around axes)");
    plotm->setXAxisIsTime();
    plotm->setYFormat(1, 6);
    plotm->setYLabel("Radius (nm)");
    plotm->appendLegend("Rg");
    plotm->appendLegend("Rg/sX/N");
    plotm->appendLegend("Rg/sY/N");
    plotm->appendLegend("Rg/sZ/N");
    gyrate_.addModule(plotm);
}

void Gyrate::analyzeFrame(int frnr, const t_trxframe& fr, t_pbc* pbc, TrajectoryAnalysisModuleData* pdata)
{
    const Selection&   sel       = TrajectoryAnalysisModuleData::parallelSelection(sel_);
    AnalysisDataHandle gyrHandle = pdata->dataHandle(gyrate_);

    real weighTotal           = 0.;
    real totalMomentOfInertia = 0.;
    real gyrationRadius       = 0.;
    RVec centerOfMass         = { 0., 0., 0. };
    RVec momentOfInertia      = { 0., 0., 0. };
    RVec gyrationRadiusOfAxes = { 0., 0., 0. };

    for (int i = 0; i < sel.posCount(); ++i)
    {
        const SelectionPosition& p  = sel.position(i);
        const real               wf = getWeighFactor(p, gMode_);
        RVec                     x(p.x());
        centerOfMass += x * wf;
        weighTotal += wf;
    }
    centerOfMass /= weighTotal;

    for (int i = 0; i < sel.posCount(); ++i)
    {
        RVec                     dx;
        const SelectionPosition& p = sel.position(i);

        if (pbc != nullptr)
        {
            pbc_dx(pbc, centerOfMass.as_vec(), p.x(), dx.as_vec());
        }
        else
        {
            rvec_sub(centerOfMass.as_vec(), p.x(), dx.as_vec());
        }


        totalMomentOfInertia += dx.norm2() * getWeighFactor(p, gMode_);

        momentOfInertia[XX] += (dx[YY] * dx[YY] + dx[ZZ] * dx[ZZ]) * getWeighFactor(p, gMode_);
        momentOfInertia[YY] += (dx[XX] * dx[XX] + dx[ZZ] * dx[ZZ]) * getWeighFactor(p, gMode_);
        momentOfInertia[ZZ] += (dx[XX] * dx[XX] + dx[YY] * dx[YY]) * getWeighFactor(p, gMode_);
    }
    gyrationRadius           = std::sqrt(totalMomentOfInertia / weighTotal);
    gyrationRadiusOfAxes[XX] = std::sqrt(momentOfInertia[XX] / weighTotal);
    gyrationRadiusOfAxes[YY] = std::sqrt(momentOfInertia[YY] / weighTotal);
    gyrationRadiusOfAxes[ZZ] = std::sqrt(momentOfInertia[ZZ] / weighTotal);

    gyrHandle.startFrame(frnr, fr.time);

    gyrHandle.setPoint(0, gyrationRadius);
    gyrHandle.setPoint(1, gyrationRadiusOfAxes[XX]);
    gyrHandle.setPoint(2, gyrationRadiusOfAxes[YY]);
    gyrHandle.setPoint(3, gyrationRadiusOfAxes[ZZ]);

    gyrHandle.finishFrame();
}


void Gyrate::finishAnalysis(int /*nframes*/) {}


void Gyrate::writeOutput() {}

} // namespace

const char GyrateInfo::name[]             = "gyrate";
const char GyrateInfo::shortDescription[] = "Calculate radius of gyration of a molecule";

TrajectoryAnalysisModulePointer GyrateInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Gyrate);
}

} // namespace analysismodules

} // namespace gmx
