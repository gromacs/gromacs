/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares gmx::TrajectoryAnalysisCommandLineRunner.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_trajectoryanalysis
 */
#ifndef GMX_TRAJECTORYANALYSIS_CMDLINERUNNER_H
#define GMX_TRAJECTORYANALYSIS_CMDLINERUNNER_H

#include "../utility/common.h"

namespace gmx
{

class HelpWriterContext;
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
         * Sets whether the runner will print the copyright header.
         *
         * \param[in] bPrint  Whether to print the copyright header.
         *
         * By default, the copyright header is printed.
         * This is used internally when executing the runner in a context where
         * the copyright has already been printed at a higher level.
         *
         * Does not throw.
         */
        void setPrintCopyright(bool bPrint);
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
        /*! \brief
         * Prints help for the module, including common options from the runner.
         *
         * \param[in] context  Context object for writing the help.
         * \throws    std::bad_alloc if out of memory.
         * \throws    FileIOError on any I/O error.
         */
        void writeHelp(const HelpWriterContext &context);

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
