/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares gmx::energyanalysis::EnergyHandler
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#ifndef GMX_ENERGYANALYSIS_ENERGYANALYSISRUNNER_H
#define GMX_ENERGYANALYSIS_ENERGYANALYSISRUNNER_H

#include <functional>
#include <memory>

#include "analysismodule.h"

namespace gmx
{

class CommandLineModuleManager;
class ICommandLineOptionsModule;

namespace energyanalysis
{

/*! \libinternal
 * \brief
 * Runner for command-line energy analysis tools.
 *
 * This class provides static methods to implement a command-line analysis
 * program, given an EnergyAnalysisModule object (or a factory of such).
 * It takes care of common command-line parameters and reading the energy
 * files.
 */
class EnergyAnalysisRunner
{
    public:
        /*! \brief
         * Factory method type for creating an energy analysis module.
         *
         * This method allows the module creation to be postponed to the point
         * where the module is needed, reducing initialization costs in, e.g.,
         * the `gmx` binary, and simplifying exception handling.
         */
        typedef std::function<EnergyAnalysisModulePointer()>
            ModuleFactoryMethod;

        /*! \brief
         * Registers a command-line module that runs a given module.
         *
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
         * Create a command-line module that runs the provided analysis module.
         *
         * \param[in]  module     Module to run.
         * \returns    Command-line module that runs the provided analysis
         *      module.
         * \throws std::bad_alloc if out of memory.
         *
         * This is mainly provided for testing purposes that want to bypass
         * CommandLineModuleManager.
         */
        static std::unique_ptr<ICommandLineOptionsModule>
        createModule(EnergyAnalysisModulePointer module);

    private:
        // Prevent instantiation.
        EnergyAnalysisRunner() {}
};

} // namespace energyanalysis

} // namespace gmx

#endif
