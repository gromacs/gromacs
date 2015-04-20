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

#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/utility/classhelpers.h"

namespace gmx
{

class CommandLineModuleManager;
class CommandLineHelpContext;

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
         * Factory method type for creating a trajectory analysis module.
         *
         * This method allows the module creation to be postponed to be inside
         * the try/catch block in runAsMain()/registerModule() implementation
         * methods and still keep the implementation out of the header, making
         * the ABI more stable.
         */
        typedef TrajectoryAnalysisModulePointer (*ModuleFactoryMethod)();

        /*! \brief
         * Implements a main() method that runs a given module.
         *
         * \tparam ModuleType  Trajectory analysis module.
         * \param  argc        \c argc passed to main().
         * \param  argv        \c argv passed to main().
         *
         * This method abstracts away all the logic required to implement a
         * main() method in user tools, allowing that to be changed without
         * requiring changes to the tools themselves.
         *
         * \p ModuleType should be default-constructible and derive from
         * TrajectoryAnalysisModule.
         *
         * Does not throw.  All exceptions are caught and handled internally.
         */
        template <class ModuleType>
        static int runAsMain(int argc, char *argv[])
        {
            return runAsMain(argc, argv, &createModule<ModuleType>);
        }
        /*! \brief
         * Registers a command-line module that runs a given module.
         *
         * \tparam ModuleType  Trajectory analysis module.
         * \param  manager     Manager to register the module to.
         * \param  name        Name of the module to register.
         * \param  description One-line description for the module to register.
         *
         * \p ModuleType should be default-constructible and derive from
         * TrajectoryAnalysisModule.
         *
         * \p name and \p descriptions must be string constants or otherwise
         * stay valid for the duration of the program execution.
         */
        template <class ModuleType>
        static void registerModule(CommandLineModuleManager *manager,
                                   const char *name, const char *description)
        {
            registerModule(manager, name, description, &createModule<ModuleType>);
        }
        /*! \brief
         * Registers a command-line module that runs a given module.
         *
         * \tparam ModuleType  Trajectory analysis module.
         * \param  manager     Manager to register the module to.
         * \param  name        Name of the module to register.
         * \param  description One-line description for the module to register.
         * \param  factory     Function that creates the module on demand.
         *
         * \p name and \p descriptions must be string constants or otherwise
         * stay valid for the duration of the program execution.
         *
         * Implements the template registerModule() method, but can also be
         * used independently.
         */
        static void registerModule(CommandLineModuleManager *manager,
                                   const char *name, const char *description,
                                   ModuleFactoryMethod factory);

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
         * Sets whether default index groups are initialized.
         *
         * This is intended only for internal unit testing purposes to avoid
         * repeated, unnecessary initialization of the default groups, which
         * can be expensive under, e.g., valgrind.
         *
         * Does not throw.
         */
        void setUseDefaultGroups(bool bUseDefaults);
        /*! \brief
         * Sets the default debugging level for selections.
         *
         * \param[in] debuglevel  Level of debugging verbosity.
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
        void writeHelp(const CommandLineHelpContext &context);

    private:
        /*! \brief
         * Creates a trajectory analysis module of a given type.
         *
         * \tparam ModuleType  Module to create.
         */
        template <class ModuleType>
        static TrajectoryAnalysisModulePointer createModule()
        {
            return TrajectoryAnalysisModulePointer(new ModuleType());
        }

        //! Implements the template runAsMain() method.
        static int runAsMain(int argc, char *argv[],
                             ModuleFactoryMethod factory);

        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
