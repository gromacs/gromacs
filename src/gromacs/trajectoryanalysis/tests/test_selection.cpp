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
/*! \internal \file
 * \brief Testing/debugging tool for the selection engine.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/format.h"
#include "gromacs/utility/programinfo.h"

namespace gmx
{

class SelectionTester : public TrajectoryAnalysisModule
{
    public:
        SelectionTester();
        virtual ~SelectionTester();

        virtual Options &initOptions(TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        void printSelections();

        Options                  _options;
        SelectionList            _selections;
        int                      _nmaxind;
};

SelectionTester::SelectionTester()
    : _options("testing", "Selection testing and debugging"),
      _nmaxind(20)
{
}

SelectionTester::~SelectionTester()
{
}

void
SelectionTester::printSelections()
{
    fprintf(stderr, "\nSelections:\n");
    for (size_t g = 0; g < _selections.size(); ++g)
    {
        _selections[g].printDebugInfo(stderr, _nmaxind);
    }
    fprintf(stderr, "\n");
}

Options &
SelectionTester::initOptions(TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "This is a test program for selections."
    };

    _options.setDescription(concatenateStrings(desc));

    _options.addOption(SelectionOption("select").storeVector(&_selections)
                           .required().multiValue()
                           .description("Selections to test"));
    _options.addOption(IntegerOption("pmax").store(&_nmaxind)
                           .description("Maximum number of indices to print in lists (-1 = print all)"));

    return _options;
}

void
SelectionTester::initAnalysis(const TrajectoryAnalysisSettings &/*settings*/,
                              const TopologyInformation &/*top*/)
{
    printSelections();
}

void
SelectionTester::analyzeFrame(int /*frnr*/, const t_trxframe &/*fr*/, t_pbc * /*pbc*/,
                              TrajectoryAnalysisModuleData * /*pdata*/)
{
    fprintf(stderr, "\n");
    for (size_t g = 0; g < _selections.size(); ++g)
    {
        const Selection &sel = _selections[g];
        int n;

        fprintf(stderr, "  Atoms (%d pcs):", sel.atomCount());
        n = sel.atomCount();
        if (_nmaxind >= 0 && n > _nmaxind)
        {
            n = _nmaxind;
        }
        ConstArrayRef<int> atoms = sel.atomIndices();
        for (int i = 0; i < n; ++i)
        {
            fprintf(stderr, " %d", atoms[i]+1);
        }
        if (n < sel.atomCount())
        {
            fprintf(stderr, " ...");
        }
        fprintf(stderr, "\n");

        fprintf(stderr, "  Positions (%d pcs):\n", sel.posCount());
        n = sel.posCount();
        if (_nmaxind >= 0 && n > _nmaxind)
        {
            n = _nmaxind;
        }
        for (int i = 0; i < n; ++i)
        {
            const SelectionPosition &p = sel.position(i);
            fprintf(stderr, "    (%.2f,%.2f,%.2f) r=%d, m=%d, n=%d\n",
                    p.x()[XX], p.x()[YY], p.x()[ZZ],
                    p.refId(), p.mappedId(), p.atomCount());
        }
        if (n < sel.posCount())
        {
            fprintf(stderr, "    ...\n");
        }
    }
    fprintf(stderr, "\n");
}

void
SelectionTester::finishAnalysis(int /*nframes*/)
{
    printSelections();
}

void
SelectionTester::writeOutput()
{
}

}

/*! \internal \brief
 * The main function for the selection testing tool.
 */
int
main(int argc, char *argv[])
{
    gmx::ProgramInfo::init(argc, argv);
    try
    {
        gmx::SelectionTester module;
        gmx::TrajectoryAnalysisCommandLineRunner runner(&module);
        runner.setSelectionDebugLevel(1);
        return runner.run(argc, argv);
    }
    catch (const std::exception &ex)
    {
        gmx::printFatalErrorMessage(stderr, ex);
        return 1;
    }
}
