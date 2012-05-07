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
 * \brief
 * Implements gmx::CommandLineModuleManager.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_commandline
 */
#include "cmdlinemodulemanager.h"

#include <cstdio>

#include <map>
#include <string>
#include <utility>

#include "gromacs/legacyheaders/statutil.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/utility/gmxassert.h"

namespace gmx
{

/********************************************************************
 * CommandLineModuleManager::Impl
 */

/*! \internal \brief
 * Private implementation class for CommandLineModuleManager.
 *
 * \ingroup module_commandline
 */
class CommandLineModuleManager::Impl
{
    public:
        //! Container for mapping module names to module objects.
        typedef std::map<std::string, CommandLineModulePointer> ModuleMap;

        //! Prints usage message to stderr.
        void printUsage() const;

        /*! \brief
         * Maps module names to module objects.
         *
         * Owns the contained modules.
         */
        ModuleMap               modules_;
};

void CommandLineModuleManager::Impl::printUsage() const
{
    fprintf(stderr, "Usage: %s <command> [<args>]\n",
            ShortProgram());
}

/********************************************************************
 * CommandLineModuleManager
 */

CommandLineModuleManager::CommandLineModuleManager()
    : impl_(new Impl)
{
}

CommandLineModuleManager::~CommandLineModuleManager()
{
}

void CommandLineModuleManager::addModule(CommandLineModulePointer module)
{
    GMX_ASSERT(impl_->modules_.find(module->name()) == impl_->modules_.end(),
               "Attempted to register a duplicate module name");
    impl_->modules_.insert(std::make_pair(std::string(module->name()), module));
}

int CommandLineModuleManager::run(int argc, char *argv[])
{
    if (argc < 2)
    {
        fprintf(stderr, "\n");
        impl_->printUsage();
        return 2;
    }
    // TODO: Accept unambiguous prefixes?
    Impl::ModuleMap::const_iterator module = impl_->modules_.find(argv[1]);
    if (module == impl_->modules_.end())
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Unknown command: %s\n", argv[1]);
        impl_->printUsage();
        return 2;
    }
    return module->second->run(argc - 1, argv + 1);
}

} // namespace gmx
