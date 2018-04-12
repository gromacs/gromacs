/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include <algorithm>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/histogram.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/applied-forces/scattering/scattering-forces.h"
#include "gromacs/compat/make_unique.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/math/vec.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/topologyinformation.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

//! Type of scattering calculation to perform
enum ScatterType {
    enXray, enNeutron
};

//! String values corresponding to ScatterType
const char               *ScatterEnum[] = {
    "xray", "neutron"
};

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

    private:
        SelectionList                    sel_;
        std::string                      fnScattering_;
        std::string                      fnAll_;
        std::string                      fnReferenceScattering_;
        std::string                      fnStructureFactor_;
        double                           startQ_;
        double                           endQ_;
        double                           qSpacing_;
        int                              numQ_;
        bool                             bNormalize_;
        bool                             bLogarithm_;
        ScatterType                      scatterEnum_;
        AnalysisData                     distances_;
        AnalysisDataAverageModulePointer averageModule_;
        std::vector<double>              qList_;
        std::vector<double>              refIntensities_;
        std::vector<double>              refErrors_;
        std::vector<double>              scatteringDifferences_;
        std::vector<double>              scatteringAtQ0_;

        friend class ComputeScattering;
        std::vector < std::unique_ptr < ComputeScattering>>  computeScattering_;
        void readReferenceScattering();
        // Copy and assign disallowed by base.
};

Scattering::Scattering()
    : startQ_(0.0), endQ_(0.0),
      qSpacing_(0.0), numQ_(0),
      bNormalize_(true), bLogarithm_(false),
      scatterEnum_(enXray),
      averageModule_(compat::make_unique<AnalysisDataAverageModule>())
{
    distances_.addModule(averageModule_);
    distances_.setMultipoint(true);
    registerAnalysisDataset(&distances_, "dist");
}

void
Scattering::initOptions(IOptionsContainer *options, TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] calculates neutron and x-ray scattering curves.[PAR]",

        "The scattering intensity, I(q), as a function of scattering angle, ",
        "q can be written out with [TT]-o[TT]. It is possible to set the ",
        "start of the q-value range with [TT]-startq[TT] and the end of the ",
        "q-value range with [TT]-endq[TT]. The number of scattering points ",
        "can be set with [TT]-numq[TT] or the spacing between q-values can ",
        "be set with [TT]-qspacing[TT]. If one of these values is set the ",
        "other must explicitly be set to zero. Comparison to a reference ",
        "scattering curve can also be performed, but the q-values must be ",
        "equally spaced in the reference file provided by option ",
        "[TT]-reference[TT].[PAR]",

    };

    settings->setHelpText(desc);
    options->addOption(FileNameOption("o").filetype(eftPlot).outputFile()
                           .store(&fnScattering_).defaultBasename("scattering")
                           .description("Scattering intensity as a function of q"));
    options->addOption(SelectionOption("sel").storeVector(&sel_).required().multiValue()
                           .description("Positions to calculate distances for"));
    options->addOption(DoubleOption("startq").store(&startQ_)
                           .defaultValue(0.0)
                           .description("smallest q value (nm)"));
    options->addOption(DoubleOption("endq").store(&endQ_)
                           .defaultValue(0.5)
                           .description("largest q value (nm)"));
    options->addOption(DoubleOption("qspacing").store(&qSpacing_)
                           .defaultValue(0.0)
                           .description("spacing of q values (nm)"));
    options->addOption(IntegerOption("numq").store(&numQ_)
                           .defaultValue(256)
                           .description("number of q values to evaluate"));
    options->addOption(BooleanOption("norm").store(&bNormalize_)
                           .defaultValue(true)
                           .description("normalize scattering intensities"));
    options->addOption(BooleanOption("log").store(&bLogarithm_)
                           .defaultValue(false)
                           .description("take the log of scattering intensities"));
    options->addOption(EnumOption<ScatterType>("scatter-type").enumValue(ScatterEnum)
                           .store(&scatterEnum_)
                           .defaultValue(enXray)
                           .description("perform x-ray or neutron scattering"));
    options->addOption(FileNameOption("ref").filetype(eftGenericData)
                           .inputFile().store(&fnReferenceScattering_)
                           .defaultBasename("reference")
                           .description("reference scattering values"));
}

void
Scattering::optionsFinished(TrajectoryAnalysisSettings * /*settings*/)
{
    // only use q's from reference file if provided
    if (!fnReferenceScattering_.empty())
    {
        readReferenceScattering();
        startQ_   = qList_[0];
        qSpacing_ = qList_[1] - startQ_;
        auto spacingLambda = [this](double a, double b){
                return a-b != qSpacing_;
            };
        if (adjacent_find(qList_.begin(), qList_.end(), spacingLambda) != qList_.end())
        {
            GMX_THROW(InconsistentInputError("q spacing in ref must be constant"));
        }
    }

    if (startQ_ < 0.0)
    {
        GMX_THROW(InconsistentInputError("-startq cannot be less than 0"));
    }

    if (endQ_ < startQ_)
    {
        GMX_THROW(InconsistentInputError("-endq cannot be less than sQ"));
    }

    if (numQ_ > 0 && qSpacing_ > 0)
    {
        GMX_THROW(InconsistentInputError("-numq and -qspacing cannot both be set"));
    }
    if (qSpacing_ <= 0 && numQ_ <= 0)
    {
        GMX_THROW(InconsistentInputError("-numq and -qspacing cannot both be set to 0"));
    }

    //actually set the q values to be used
    if ((numQ_ > 0) && !(qSpacing_ > 0))
    {
        qSpacing_ = (endQ_ - startQ_) / (numQ_ - 1);
        for (int i = 0; i < numQ_; ++i)
        {
            qList_.push_back(startQ_ + qSpacing_ * i);
        }
    }

    if (qSpacing_  > 0 && !(numQ_ > 0))
    {
        double qVal = startQ_;
        while (qVal < endQ_)
        {
            qList_.push_back(qVal);
            qVal += qSpacing_;
        }
        if (qVal != endQ_)
        {
            GMX_THROW(gmx::InvalidInputError("The requested -qspacing did not result in an integer number of q values"));
        }
        numQ_ = qList_.size();
    }
}

void Scattering::readReferenceScattering()
{
    // read scattering parameters from file if provided, else use defaults

    if (fnReferenceScattering_.empty())
    {
        GMX_THROW(gmx::InvalidInputError("The reference scattering file was not properly read\n"));
    }

    FILE *fp = gmx_fio_fopen(fnReferenceScattering_.c_str(), "r");
    char  line[1000];
    // all non-header lines
    while (get_a_line(fp, line, 1000))
    {
        float scatteringAngle;
        float intensity;
        float experimentalError;
        if (sscanf(line, "%f %f %f", &scatteringAngle, &intensity, &experimentalError) == 3)
        {
            qList_.push_back(scatteringAngle);
            refIntensities_.push_back(intensity);
            refErrors_.push_back(experimentalError);
        }
        else if (sscanf(line, "%f %f", &scatteringAngle, &intensity) == 2)
        {
            qList_.push_back(scatteringAngle);
            refIntensities_.push_back(intensity);
        }
    }
    gmx_fio_fclose(fp);
}

void
Scattering::initAnalysis(const TrajectoryAnalysisSettings &settings,
                         const TopologyInformation        &top)
{
    const gmx_mtop_t *mtop  = top.mtop();
    const t_atoms     atoms = gmx_mtop_global_atoms(mtop);
    distances_.setDataSetCount(sel_.size());
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        distances_.setColumnCount(g, qList_.size());
        if (scatterEnum_ == enXray)
        {
            computeScattering_.push_back(compat::make_unique<XrayInfo>(atoms));
        }
        if (scatterEnum_ == enNeutron)
        {
            computeScattering_.push_back(compat::make_unique<NeutronInfo>(atoms));
        }
        scatteringAtQ0_.push_back(computeScattering_[g]->computeS0(sel_[g]));
    }

    if (!fnScattering_.empty())
    {
        averageModule_->setXAxis(startQ_, qSpacing_);
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnScattering_);
        plotm->setXFormat(2, 10); //is this precision reasonable?
        plotm->setYFormat(1, 10); //is this precision reasonable?
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
    AnalysisDataHandle     scatterHandle   = pdata->dataHandle(distances_);
    const SelectionList   &sel             = pdata->parallelSelections(sel_);
    scatterHandle.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        scatterHandle.selectDataSet(g);
        //loop over q's and just store the debeye evaluations in the histogram
        for (size_t qi = 0; qi != qList_.size(); ++qi)
        {
            double Debeye;
            if (qList_[qi] == 0.0)
            {
                Debeye = scatteringAtQ0_[g];
            }
            else
            {
                Debeye = computeScattering_[g]->computeIntensity(pbc, qList_[qi], sel[g]);
            }
            if (!fnReferenceScattering_.empty())
            {
                Debeye -= refIntensities_[qi];
            }
            if (bNormalize_)
            {
                Debeye /= scatteringAtQ0_[g];
            }
            if (bLogarithm_)
            {
                Debeye = log(Debeye);
            }
            scatterHandle.setPoint(qi, Debeye, true);
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
