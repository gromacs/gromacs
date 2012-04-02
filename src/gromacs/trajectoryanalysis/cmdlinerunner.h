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
/*! \file
 * \brief
 * Declares gmx::TrajectoryAnalysisCommandLineRunner.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_CMDLINERUNNER_H
#define GMX_TRAJECTORYANALYSIS_CMDLINERUNNER_H

#include "../utility/common.h"

namespace gmx
{

class TrajectoryAnalysisModule;

/*! \brief
 * Runner class for command-line analysis tools.
 *
 * This class implements a command-line analysis program, given a
 * TrajectoryAnalysisModule object.  It takes care of converting command-line
 * parameters to a form understood by the module, as well as parsing common
 * options, initializing and evaluating selections, and looping over trajectory
 * frames.
 *
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
class TrajectoryAnalysisCommandLineRunner
{
    public:
        /*! \brief
         * Create a new runner with the provided module.
         *
         * \param  module  Analysis module to run using the runner.
         * \throws std::bad_alloc if out of memory.
         *
         * The caller should ensure that the provided module is not destroyed
         * while the runner exists.
         */
        TrajectoryAnalysisCommandLineRunner(TrajectoryAnalysisModule *module);
        ~TrajectoryAnalysisCommandLineRunner();

        /*! \brief
         * Sets the default debugging level for selections.
         *
         * This is intended only for use by internal debugging tools.
         *
         * Does not throw.
         *
         * \see SelectionCollection::setDebugLevel()
         */
        void setSelectionDebugLevel(int debuglevel);
        /*! \brief
         * Parses options from the given command line and runs the analysis.
         *
         * \throws  multiple  Exceptions are used to indicate errors.
         * \returns Zero on success.
         */
        int run(int argc, char *argv[]);

    private:
        class Impl;

        PrivateImplPointer<Impl> _impl;
};

} // namespace gmx

#endif
