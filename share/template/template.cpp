/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
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

        virtual Options *initOptions(TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TopologyInformation &top);

        virtual TrajectoryAnalysisModuleData *startFrames(
                    AnalysisDataParallelOptions opt,
                    const SelectionCollection &selections);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        class ModuleData;

        Options                      _options;
        std::string                  _fnDist;
        double                       _cutoff;
        Selection                   *_refsel;
        std::vector<Selection *>     _sel;
        AnalysisData                 _data;
        AnalysisDataAverageModule   *_avem;
};

/*! \brief
 * Frame-local data needed in analysis.
 */
class AnalysisTemplate::ModuleData : public TrajectoryAnalysisModuleData
{
    public:
        ModuleData(TrajectoryAnalysisModule *module,
                   AnalysisDataParallelOptions opt,
                   const SelectionCollection &selections)
            : TrajectoryAnalysisModuleData(module, opt, selections),
              _nb(NULL)
        {
        }

        virtual ~ModuleData()
        {
            delete _nb;
        }

        virtual void finish()
        {
            finishDataHandles();
        }

        NeighborhoodSearch          *_nb;
};


AnalysisTemplate::AnalysisTemplate()
    : _options("template", "Template options"), _cutoff(0.0),
      _refsel(NULL), _avem(NULL)
{
}


Options *
AnalysisTemplate::initOptions(TrajectoryAnalysisSettings *settings)
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
        "average distances from the a reference group to one or more",
        "analysis groups.",
        NULL
    };

    _options.setDescription(desc);

    _options.addOption(FileNameOption("o")
        .filetype(eftPlot).writeOnly()
        .store(&_fnDist).defaultValueIfSet("avedist")
        .description("Average distances from reference group"));

    _options.addOption(SelectionOption("reference")
        .store(&_refsel).required()
        .description("Reference group to calculate distances from"));
    _options.addOption(SelectionOption("select")
        .storeVector(&_sel).required().multiValue()
        .description("Groups to calculate distances to"));

    _options.addOption(DoubleOption("cutoff").store(&_cutoff)
        .description("Cutoff for distance calculation (0 = no cutoff)"));

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);

    return &_options;
}


void
AnalysisTemplate::initAnalysis(const TopologyInformation & /*top*/)
{
    _data.setColumns(_sel.size());
    registerAnalysisDataset(&_data, "avedist");

    _avem = new AnalysisDataAverageModule();
    _data.addModule(_avem);

    if (!_fnDist.empty())
    {
        AnalysisDataPlotModule *plotm = new AnalysisDataPlotModule(_options);
        plotm->setFileName(_fnDist);
        plotm->setTitle("Average distance");
        plotm->setXLabel("Time (ps)");
        plotm->setYLabel("Distance (nm)");
        _data.addModule(plotm);
    }
}


TrajectoryAnalysisModuleData *
AnalysisTemplate::startFrames(AnalysisDataParallelOptions opt,
                              const SelectionCollection &selections)
{
    ModuleData *pdata = new ModuleData(this, opt, selections);
    int rc = NeighborhoodSearch::create(&pdata->_nb, _cutoff, _refsel->posCount());
    if (rc != 0)
    {
        delete pdata;
        // FIXME: Use exceptions in the neighborhood search API
        GMX_THROW(InternalError("Neighborhood search initialization failed"));
    }
    return pdata;
}


void
AnalysisTemplate::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                               TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle *dh = pdata->dataHandle("avedist");
    NeighborhoodSearch *nb = static_cast<ModuleData *>(pdata)->_nb;

    int rc = nb->init(pbc, _refsel->positions());
    if (rc != 0)
    {
        // FIXME: Use exceptions in the neighborhood search API
        GMX_THROW(InternalError("Neighborhood search frame initialization failed"));
    }
    dh->startFrame(frnr, fr.time);
    for (size_t g = 0; g < _sel.size(); ++g)
    {
        Selection *sel = pdata->parallelSelection(_sel[g]);
        int   nr = sel->posCount();
        real  frave = 0.0;
        for (int i = 0; i < nr; ++i)
        {
            frave += nb->minimumDistance(sel->x(i));
        }
        frave /= nr;
        dh->addPoint(g, frave);
    }
    dh->finishFrame();
}


void
AnalysisTemplate::finishAnalysis(int /*nframes*/)
{
}


void
AnalysisTemplate::writeOutput()
{
    // We print out the average of the mean distances for each group.
    for (size_t g = 0; g < _sel.size(); ++g)
    {
        fprintf(stderr, "Average mean distance for '%s': %.3f nm\n",
                _sel[g]->name(), _avem->average(g));
    }
}

/*! \brief
 * The main function for the analysis template.
 */
int
main(int argc, char *argv[])
{
    try
    {
        AnalysisTemplate module;
        TrajectoryAnalysisCommandLineRunner runner(&module);
        return runner.run(argc, argv);
    }
    catch (std::exception &ex)
    {
        fprintf(stderr, "%s", gmx::formatErrorMessage(ex).c_str());
        return 1;
    }
}
