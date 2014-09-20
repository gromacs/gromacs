/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
 * \brief Testing/debugging tool for the selection engine.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/options/basicoptions.h"
#include "gromacs/options/options.h"
#include "gromacs/selection/selection.h"
#include "gromacs/selection/selectionoption.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/exceptions.h"

namespace gmx
{

class SelectionTester : public TrajectoryAnalysisModule
{
    public:
        SelectionTester();
        virtual ~SelectionTester();

        virtual void initOptions(Options                    *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        void printSelections();

        SelectionList            selections_;
        int                      nmaxind_;
};

SelectionTester::SelectionTester()
    : TrajectoryAnalysisModule("testing", "Selection testing and debugging"),
      nmaxind_(20)
{
}

SelectionTester::~SelectionTester()
{
}

void
SelectionTester::printSelections()
{
    fprintf(stderr, "\nSelections:\n");
    for (size_t g = 0; g < selections_.size(); ++g)
    {
        selections_[g].printDebugInfo(stderr, nmaxind_);
    }
    fprintf(stderr, "\n");
}

void
SelectionTester::initOptions(Options                   *options,
                             TrajectoryAnalysisSettings * /*settings*/)
{
    static const char *const desc[] = {
        "This is a test program for selections."
    };

    options->setDescription(desc);

    options->addOption(SelectionOption("select").storeVector(&selections_)
                           .required().multiValue()
                           .description("Selections to test"));
    options->addOption(IntegerOption("pmax").store(&nmaxind_)
                           .description("Maximum number of indices to print in lists (-1 = print all)"));
}

void
SelectionTester::initAnalysis(const TrajectoryAnalysisSettings & /*settings*/,
                              const TopologyInformation        & /*top*/)
{
    printSelections();
}

void
SelectionTester::analyzeFrame(int /*frnr*/, const t_trxframe & /*fr*/, t_pbc * /*pbc*/,
                              TrajectoryAnalysisModuleData * /*pdata*/)
{
    fprintf(stderr, "\n");
    for (size_t g = 0; g < selections_.size(); ++g)
    {
        const Selection &sel = selections_[g];
        int              n;

        fprintf(stderr, "  Atoms (%d pcs):", sel.atomCount());
        n = sel.atomCount();
        if (nmaxind_ >= 0 && n > nmaxind_)
        {
            n = nmaxind_;
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
        if (nmaxind_ >= 0 && n > nmaxind_)
        {
            n = nmaxind_;
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
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<gmx::SelectionTester>(argc, argv);
}
