/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018 by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::scattering.
 *
 * \author Joe Jordan <e.jjordan12@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "scattering.h"

#include <map>
#include <string>
#include <vector>
#include <iostream>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/applied-forces/scattering/scattering.h"

namespace gmx
{

namespace analysismodules
{

namespace
{


class Scattering : public TrajectoryAnalysisModule
{
    public:
        Scattering();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void optionsFinished(TrajectoryAnalysisSettings *settings);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);
        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();
        void initScattering(const SelectionList       &sel,
                            const TopologyInformation &top);
    private:
        SelectionList                            sel_;
        std::string                              fnScattering_;
        std::string                              fnAll_;
        std::string                              fnStructureFactor_;
        double                                   startQ_;    
        double                                   endQ_;
        double                                   qSpacing_;
        int                                      numQ_;    
        AnalysisData                             distances_;
        AnalysisDataAverageModulePointer         averageModule_;
        std::vector <double>                     qList_;
        std::vector<double>                      scatteringAtQ0_;

        friend class ComputeScattering;
        std::vector<ComputeScattering>           computeScattering_;
        // Copy and assign disallowed by base.
};

Scattering::Scattering()
{
    averageModule_.reset(new AnalysisDataAverageModule());
    distances_.addModule(averageModule_);

    distances_.setMultipoint(true);
    registerAnalysisDataset(&distances_, "dist");
}

void
Scattering::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] calculates neutron scatting curves.[PAR]",
        "[TT]-o[tt] writes scattering intensity, I(q), as a function of ",
        "scattering angle, q.",
    };

    settings->setHelpText(desc);
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnScattering_).defaultBasename("scattering")
                           .description("Scattering intensity as a function of q"));
    options->addOption(SelectionOption("sel").storeVector(&sel_).required().multiValue()
                       .description("Positions to calculate distances for"));
    options->addOption(DoubleOption("startq").store(&startQ_)
                       .defaultValue(0.0)
                       .description("smallest q value"));
    options->addOption(DoubleOption("endq").store(&endQ_)
                       .defaultValue(2.0)
                       .description("largest q value"));
    options->addOption(DoubleOption("qspacing").store(&qSpacing_)
                       .defaultValue(0.0)
                       .description("spacing of q values"));
    options->addOption(IntegerOption("numq").store(&numQ_)
                       .defaultValue(200)
                       .description("number of q values to evaluate"));
}

void
Scattering::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    if (startQ_ < 0.0)
    {
        GMX_THROW(InconsistentInputError("-startq cannot be less than 0"));
    }

    if (endQ_ < startQ_)
    {
        GMX_THROW(InconsistentInputError("-endq cannot be less than sQ"));
    }

    if (numQ_>0 && qSpacing_>0)
    {
        GMX_THROW(InconsistentInputError("-numq and qspacing cannot both be set"));
    }

    //actually set the q values to be used
    if ((numQ_>0) && ! (qSpacing_>0))
    {
        qSpacing_ = (endQ_ - startQ_) / numQ_;
        for (int i = 0; i < numQ_; ++i)
        {
            qList_.push_back(startQ_ + qSpacing_ * i);
        }
    }
    // need:  if (qSpacing_ and not numQ_)
}

/*! \brief
 * sets of the list of scattering factors, indexed with selections
 */
void Scattering::initScattering(const SelectionList       &sel,
                                const TopologyInformation &top)
{
    for (size_t g = 0; g < sel.size(); ++g)
    {
        auto ComputeScattering = makeScattering(sel_[g], top);
        computeScattering_.push_back(ComputeScattering);
    }
}

void
Scattering::initAnalysis(const TrajectoryAnalysisSettings &settings,
                         const TopologyInformation        &top)
{
    initScattering(sel_, top);

    distances_.setDataSetCount(sel_.size());
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        distances_.setColumnCount(g, qList_.size());
        scatteringAtQ0_.push_back(computeScattering_[g].computeS0());
    }

    if (!fnScattering_.empty())
    {
        averageModule_->setXAxis(startQ_, qSpacing_);
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnScattering_);
        plotm->setXFormat(2, 4); //is this precision reasonable?
        plotm->setYFormat(1, 6); //is this precision reasonable?
        plotm->setTitle("Scattering intensity: I(q)");
        plotm->setXLabel("q (1/nm)");
        plotm->setYLabel("Intensity");
        averageModule_->addModule(plotm);
    }
}

void
Scattering::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                       TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle    scatterHandle   = pdata->dataHandle(distances_);
    const SelectionList   &sel         = pdata->parallelSelections(sel_);

    scatterHandle.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        scatterHandle.selectDataSet(g);
        //loop over q's and just store the debeye evaluations in the histogram
        for (size_t qi = 0; qi != qList_.size(); ++qi)
        {
            if (qList_[qi]==0.0)
            {
                scatterHandle.setPoint(qi, 1.0, true);
            }
            else
            {
                double Debeye = computeScattering_[g].compute_scattering(pbc, qList_[qi], sel[g]);
                Debeye /= scatteringAtQ0_[g];
                scatterHandle.setPoint(qi, Debeye, true);
            }
            scatterHandle.finishPointSet();
        }
    }    
    scatterHandle.finishFrame();
}


void
Scattering::finishAnalysis(int /*nframes*/)
{
    //nothing for now
}


void
Scattering::writeOutput()
{
    //nothing for now
}

}       // namespace

const char ScatteringInfo::name[]             = "scattering";
const char ScatteringInfo::shortDescription[] =
    "Calculate small angle scattering profiles";

TrajectoryAnalysisModulePointer ScatteringInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Scattering);
}

} // namespace analysismodules

} // namespace gmx
