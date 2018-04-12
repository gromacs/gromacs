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

namespace gmx
{

namespace analysismodules
{

namespace
{

struct xRayStructureFactors {
    int   p;      /* proton number */
    int   n;      /* neutron number */
    /* Parameters for the Cromer Mann fit */
    real  a;      /* parameter a */
    real  b;      /* parameter b */
    real  c;      /* parameter c */
    char  atomnm; /* atomname */

};

struct neutronStructureFactor {
    int    p;              /* proton number */
    int    n;              /* neutron number */
    float  scatterLength;  /* neutron scattering length*/
};

//! \brief Set the scattering lengths for the six atom types H, C, N, O, P, and S
std::map <std::string, neutronStructureFactor> set_scattering_coeff( std::string parameterFileName )
{
    // read gaussian parameters from file if provided, else use defaults
    FILE *parameterFile;
    if (!parameterFileName.empty())
    {
        fprintf(stdout, "Reading file with neutron scattering factors\n");
        parameterFile = fopen(parameterFileName.c_str(), "r");
    }
    else
    {
        GMX_THROW(FileIOError("Could not open scattering file " + parameterFileName + " for reading"));
    }
    std::map < std::string, neutronStructureFactor> neutronStructureFactorLookup;
    char line[1000];
    // all non-header lines
    while (get_a_line(parameterFile, line, 1000))
    {

        char  currentAtomType[8];
        int   p;
        int   n;
        float  scatterLength;

        if (sscanf(line, "%s %d %d %f", currentAtomType, &p, &n, &scatterLength) == 4)
        {
            neutronStructureFactor nsf = {p, n, scatterLength};
            neutronStructureFactorLookup[currentAtomType] = nsf;
        }
    }
    fclose(parameterFile);
    return neutronStructureFactorLookup;
}

class Scattering : public TrajectoryAnalysisModule
{
    public:
        Scattering();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        SelectionList                            sel_;
        std::string                              fnScattering_;
        std::string                              fnAll_;
        std::string                              fnStructureFactor_;
        double                                   startQ_;    
        double                                   endQ_;
        double                                   qSpacing_;
        int                                      numQ_;    

        AnalysisData                             scattering_;
        AnalysisDataAverageModulePointer         averageModule_;

        std::map <float, float> mass2scatterLength_;
        std::vector <real>      qList_;
        real                    scatteringAtQ0_;
        // Copy and assign disallowed by base.
};

Scattering::Scattering()
    : numQ_(100), startQ_(0.01), endQ_(2.1), qSpacing_(0.0)
{
    averageModule_.reset(new AnalysisDataAverageModule());
    scattering_.addModule(averageModule_);

    scattering_.setMultipoint(true);
    registerAnalysisDataset(&scattering_, "dist");
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
    options->addOption(DoubleOption("sQ").store(&startQ_)
                           .description("smallest q value"));
    options->addOption(DoubleOption("eQ").store(&endQ_)
                           .description("largest q value"));
    options->addOption(DoubleOption("qs").store(&qSpacing_)
                           .description("spacing of q values"));
    options->addOption(IntegerOption("nQ").store(&numQ_)
                           .description("number of q values to evaluate"));

    // optional input file
    options->addOption(FileNameOption("sfactor").filetype(eftGenericData)
                       .inputFile().required()
                       .store(&fnStructureFactor_).defaultBasename("nsfactor")
                       .description("Structure factors")); 
}

/*! \brief
 * Checks that selections conform to the expectations of the tool.
 */
void checkSelections(const SelectionList &sel)
{
    for (size_t g = 0; g < sel.size(); ++g)
    {
        // do nothing for now
        // at some point will need checks on option consistency here
    }
}

void checkQs(int numQ_, double startQ_, double endQ_,
             double &qSpacing_,  std::vector <float> &qList_)
{
    GMX_ASSERT(numQ_>0 and qSpacing_>0,
               "You cannot set both the number of q points and the q spacing");
    if (numQ_>0 and not qSpacing_>0)
    {
        qSpacing_ = (endQ_ - startQ_) / numQ_;
        for (int i = 0; i < numQ_; ++i)
        {
            qList_.push_back(startQ_ + qSpacing_ * i);
        }
    }
    // need:  if (qSpacing_ and not numQ_)
}

void
Scattering::initAnalysis(const TrajectoryAnalysisSettings &settings,
                         const TopologyInformation        &top)
{
    checkSelections(sel_);

    const t_atoms &atoms = top.topology()->atoms;
    std::map <float, std::string> mass2element;
    std::map <std::string, neutronStructureFactor> scatteringPerAtomType;

    checkQs(numQ_, startQ_, endQ_, qSpacing_, qList_);

    for (int i = 0; i < atoms.nr; ++i)
    {
        const t_atom atom = atoms.atom[i];
        if (atom.atomnumber > 0 && atom.m > 0) //atomnumber==0 and mass<-m==0 are virtual sites
        {
            mass2element.insert( std::pair<float, std::string> (atom.m, atom.elem));
        }
    }
    scatteringPerAtomType = set_scattering_coeff(fnStructureFactor_);
    for (auto j = mass2element.begin(); j != mass2element.end(); ++j)
    {
        auto mass = j->first;
        auto element = j->second;
        auto scatterLength = scatteringPerAtomType.find(element)->second.scatterLength;
        mass2scatterLength_.insert( std::pair<float, float> (mass, scatterLength));
    }
    scattering_.setDataSetCount(sel_.size());
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        scattering_.setColumnCount(g, qList_.size());
        const int posCount = sel_[g].posCount();
        for (int i = 0; i < posCount; ++i)
        {
            const SelectionPosition &pi = sel_[g].position(i);
            auto mass_i = pi.mass();
            auto scatterLength_i = mass2scatterLength_.find(mass_i)->second;
            for (int j = 0; j < i; ++j)
            {
                const SelectionPosition &pj = sel_[g].position(j);
                auto  mass_j = pj.mass();
                auto  scatterLength_j = mass2scatterLength_.find(mass_j)->second;
                scatteringAtQ0_ += scatterLength_i * scatterLength_j;
            }
        }
    }

    if (!fnScattering_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnScattering_);
        plotm->setYFormat(6, 6); //is this precision reasonable?
        plotm->setTitle("Scattering intensity: I(q)");
        plotm->setXLabel("q (1/nm)");
        plotm->setYLabel("Intensity");
        for (size_t g = 0; g < sel_.size(); ++g)
        {
            plotm->appendLegend(sel_[g].name());
        }
        averageModule_->addModule(plotm);
    }
}

void
Scattering::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                       TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle    scatterHandle   = pdata->dataHandle(scattering_);
    //may want to put in a normalization factor like in rdf to reduce overhead
    const SelectionList   &sel         = pdata->parallelSelections(sel_);

    checkSelections(sel);
    scatterHandle.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel.size(); ++g)
    {
        scatterHandle.selectDataSet(g);
        const int posCount = sel[g].posCount();
        //loop over q's and just store the debeye evaluations in the histogram
        for (int qi = 0; qi != qList_.size(); ++qi)
        {
            real Debeye = 0;
            for (int i = 0; i < posCount; ++i)
            {
                const SelectionPosition &pi = sel[g].position(i);
                auto mass_i = pi.mass();
                auto scatterLength_i = mass2scatterLength_.find(mass_i)->second;
                for (int j = 0; j < i; ++j)
                {
                    const SelectionPosition &pj = sel[g].position(j);
                    auto  mass_j = pj.mass();
                    auto  scatterLength_j = mass2scatterLength_.find(mass_j)->second;
                    real  formFactor = scatterLength_i * scatterLength_j;
                    rvec  dx;
                    if (pbc != nullptr)
                    {
                        pbc_dx(pbc, pi.x(), pj.x(), dx);
                    }
                    else
                    {
                        rvec_sub(pi.x(), pj.x(), dx);
                    }
                    real qDist = norm(dx) * qList_[qi];
                    real sinQDist = sin(qDist);
                    Debeye += formFactor * (sinQDist / qDist);
                }
            }
            Debeye /= scatteringAtQ0_;
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
