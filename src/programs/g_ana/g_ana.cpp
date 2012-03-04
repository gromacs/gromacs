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
/*! \internal \brief
 * Implements the g_ana tool.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
#include <copyrite.h>

#include "gromacs/fatalerror.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/trajectoryanalysis/modules.h"

int
main(int argc, char *argv[])
{
    bool bPrintCopyrightOnError = true;

    try
    {
        if (argc < 2)
        {
            GMX_THROW(gmx::InvalidInputError("Not enough command-line arguments"));
        }

        gmx::TrajectoryAnalysisModulePointer
            mod(gmx::createTrajectoryAnalysisModule(argv[1]));
        --argc;
        ++argv;

        gmx::TrajectoryAnalysisCommandLineRunner runner(mod.get());
        bPrintCopyrightOnError = false;
        return runner.run(argc, argv);
    }
    catch (const std::exception &ex)
    {
        if (bPrintCopyrightOnError)
        {
            CopyRight(stderr, argv[0]);
        }
        fprintf(stderr, "%s", gmx::formatErrorMessage(ex).c_str());
        return 1;
    }
}
