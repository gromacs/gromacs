/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

using namespace gmx;

/*! \brief
 * Template class to serve as a basis for user analysis tools.
 */
class AnalysisTemplate : public TrajectoryAnalysisModule
{
    public:
        AnalysisTemplate();

        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual TrajectoryAnalysisModuleDataPointer startFrames(
            const AnalysisDataParallelOptions &opt,
            const SelectionCollection         &selections);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        class ModuleData;

        std::string                      fnDist_;
        double                           cutoff_;
        Selection                        refsel_;
        SelectionList                    sel_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;
};

/*! \brief
 * Frame-local data needed in analysis.
 */
class AnalysisTemplate::ModuleData : public TrajectoryAnalysisModuleData
{
    public:
        /*! \brief
         * Initializes frame-local data.
         *
         * \param[in] module     Analysis module to use for data objects.
         * \param[in] opt        Data parallelization options.
         * \param[in] selections Thread-local selection collection.
         * \param[in] cutoff     Cutoff distance for the search
         *   (<=0 stands for no cutoff).
         * \param[in] posCount   Maximum number of reference particles.
         */
        ModuleData(TrajectoryAnalysisModule *module,
                   const AnalysisDataParallelOptions &opt,
                   const SelectionCollection &selections,
                   double cutoff, int posCount)
            : TrajectoryAnalysisModuleData(module, opt, selections),
              nb_(cutoff, posCount)
        {
        }

        virtual void finish()
        {
            finishDataHandles();
        }

        //! Neighborhood search data for distance calculation.
        NeighborhoodSearch      nb_;
};


AnalysisTemplate::AnalysisTemplate()
    : TrajectoryAnalysisModule("template", "Template analysis tool"),
      cutoff_(0.0)
{
    registerAnalysisDataset(&data_, "avedist");
}


void
AnalysisTemplate::initOptions(Options                    *options,
                              TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "This is a template for writing your own analysis tools for",
        "Gromacs. The advantage of using Gromacs for this is that you",
        "have access to all information in the topology, and your",
        "program will be able to handle all types of coordinates and",
        "trajectory files supported by Gromacs. In addition,",
        "you get a lot of functionality for free from the trajectory",
        "analysis library, including support for flexible dynamic",
        "selections. Go ahead an try it![PAR]",
        "To get started with implementing your own analysis program,",
        "follow the instructions in the README file provided.",
        "This template implements a simple analysis programs that calculates",
        "average distances from a reference group to one or more",
        "analysis groups."
    };

    options->setDescription(concatenateStrings(desc));

    options->addOption(FileNameOption("o")
                           .filetype(eftPlot).outputFile()
                           .store(&fnDist_).defaultBasename("avedist")
                           .description("Average distances from reference group"));

    options->addOption(SelectionOption("reference")
                           .store(&refsel_).required()
                           .description("Reference group to calculate distances from"));
    options->addOption(SelectionOption("select")
                           .storeVector(&sel_).required().multiValue()
                           .description("Groups to calculate distances to"));

    options->addOption(DoubleOption("cutoff").store(&cutoff_)
                           .description("Cutoff for distance calculation (0 = no cutoff)"));

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


void
AnalysisTemplate::initAnalysis(const TrajectoryAnalysisSettings &settings,
                               const TopologyInformation         & /*top*/)
{
    data_.setColumnCount(sel_.size());

    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    if (!fnDist_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnDist_);
        plotm->setTitle("Average distance");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Distance (nm)");
        data_.addModule(plotm);
    }
}


TrajectoryAnalysisModuleDataPointer
AnalysisTemplate::startFrames(const AnalysisDataParallelOptions &opt,
                              const SelectionCollection         &selections)
{
    return TrajectoryAnalysisModuleDataPointer(
            new ModuleData(this, opt, selections, cutoff_, refsel_.posCount()));
}


void
AnalysisTemplate::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle  dh     = pdata->dataHandle(data_);
    NeighborhoodSearch &nb     = static_cast<ModuleData *>(pdata)->nb_;
    const Selection    &refsel = pdata->parallelSelection(refsel_);

    nb.init(pbc, refsel.positions());
    dh.startFrame(frnr, fr.time);
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        const Selection &sel   = pdata->parallelSelection(sel_[g]);
        int              nr    = sel.posCount();
        real             frave = 0.0;
        for (int i = 0; i < nr; ++i)
        {
            SelectionPosition p = sel.position(i);
            frave += nb.minimumDistance(p.x());
        }
        frave /= nr;
        dh.setPoint(g, frave);
    }
    dh.finishFrame();
}


void
AnalysisTemplate::finishAnalysis(int /*nframes*/)
{
}


void
AnalysisTemplate::writeOutput()
{
    // We print out the average of the mean distances for each group.
    for (size_t g = 0; g < sel_.size(); ++g)
    {
        fprintf(stderr, "Average mean distance for '%s': %.3f nm\n",
                sel_[g].name(), avem_->average(g));
    }
}

/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[])
{
    ProgramInfo::init(argc, argv);
    try
    {
        AnalysisTemplate                    module;
        TrajectoryAnalysisCommandLineRunner runner(&module);
        return runner.run(argc, argv);
    }
    catch (const std::exception &ex)
    {
        gmx::printFatalErrorMessage(stderr, ex);
        return 1;
    }
}
