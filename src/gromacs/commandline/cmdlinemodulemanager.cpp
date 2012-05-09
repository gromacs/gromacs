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

// For GMX_BINARY_SUFFIX
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdio>
#include <cstring>

#include <map>
#include <string>
#include <utility>

#include "gromacs/legacyheaders/statutil.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"

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

        /*! \brief
         * Initializes the implementation class.
         *
         * \param[in] realBinaryName  Name of the binary that this manager runs.
         */
        explicit Impl(const char *realBinaryName);

        /*! \brief
         * Finds a module that matches a name.
         *
         * \param[in] name  Module name to find.
         * \returns   Iterator to the found module, or
         *      \c modules_.end() if not found.
         *
         * Does not throw.
         */
        ModuleMap::const_iterator findModuleByName(const std::string &name) const;
        /*! \brief
         * Finds a module that the name of the binary.
         *
         * \param[in] argv0  argv[0] passed to the program.
         * \throws    std::bad_alloc if out of memory.
         * \returns   Iterator to the found module, or
         *      \c modules_.end() if not found.
         *
         * Checks whether the program is invoked through a symlink whose name
         * is different from \a realBinaryName_, and if so, checks if a module
         * name matches the name of the symlink.
         */
        ModuleMap::const_iterator findModuleFromBinaryName(const std::string &argv0) const;

        //! Prints usage message to stderr.
        void printUsage(bool bModuleList) const;
        //! Prints the list of modules to stderr.
        void printModuleList() const;

        /*! \brief
         * Maps module names to module objects.
         *
         * Owns the contained modules.
         */
        ModuleMap               modules_;
        //! Real name of the binary that is running (without suffixes).
        std::string             realBinaryName_;
};

CommandLineModuleManager::Impl::Impl(const char *realBinaryName)
    : realBinaryName_(realBinaryName)
{
}

CommandLineModuleManager::Impl::ModuleMap::const_iterator
CommandLineModuleManager::Impl::findModuleByName(const std::string &name) const
{
    // TODO: Accept unambiguous prefixes?
    return modules_.find(name);
}

CommandLineModuleManager::Impl::ModuleMap::const_iterator
CommandLineModuleManager::Impl::findModuleFromBinaryName(const std::string &argv0) const
{
    // TODO: Move this logic into a common place in utility/ and remove
    // dependency on config.h from this file.
    // (most natural place would be in a location that wraps Program() etc.)
    std::string binaryName = Path::splitToPathAndFilename(argv0).second;
    if (binaryName.length() >= 4
        && binaryName.compare(binaryName.length() - 4, 4, ".exe") == 0)
    {
        binaryName.erase(binaryName.length() - 4);
    }
#ifdef GMX_BINARY_SUFFIX
    size_t suffixLength = std::strlen(GMX_BINARY_SUFFIX);
    if (suffixLength > 0 && binaryName.length() >= suffixLength
        && binaryName.compare(binaryName.length() - suffixLength, suffixLength,
                              GMX_BINARY_SUFFIX) == 0)
    {
        binaryName.erase(binaryName.length() - suffixLength);
    }
#endif
    if (binaryName == realBinaryName_)
    {
        return modules_.end();
    }
    if (binaryName.compare(0, 2, "g_") == 0)
    {
        binaryName.erase(0, 2);
    }
    return findModuleByName(binaryName);
}

void CommandLineModuleManager::Impl::printUsage(bool bModuleList) const
{
    const char *program = ShortProgram();
    fprintf(stderr, "Usage: %s <command> [<args>]\n\n", program);
    if (bModuleList)
    {
        printModuleList();
    }
    else
    {
        fprintf(stderr, "See '%s help' for list of commands.\n", program);
    }
}

void CommandLineModuleManager::Impl::printModuleList() const
{
    int maxNameLength = 0;
    ModuleMap::const_iterator module;
    for (module = modules_.begin(); module != modules_.end(); ++module)
    {
        int nameLength = static_cast<int>(module->first.length());
        if (nameLength > maxNameLength)
        {
            maxNameLength = nameLength;
        }
    }
    fprintf(stderr, "The following commands are available:\n");
    for (module = modules_.begin(); module != modules_.end(); ++module)
    {
        const char *name = module->first.c_str();
        const char *description = module->second->shortDescription();
        fprintf(stderr, "    %*s  %s\n", -maxNameLength, name, description);
    }
}


/********************************************************************
 * CommandLineHelpModule
 */

namespace internal
{

/*! \internal \brief
 * Command-line module for producing help.
 *
 * This module implements the 'help' subcommand that is automatically added by
 * CommandLineModuleManager.
 *
 * \ingroup module_commandline
 */
class CommandLineHelpModule : public CommandLineModuleInterface
{
    public:
        /*! \brief
         * Creates a help module for the given module manager.
         *
         * \param[in] manager  Manager for which this module provides help.
         *
         * Does not throw.
         */
        explicit CommandLineHelpModule(const CommandLineModuleManager &manager);

        virtual const char *name() const { return "help"; }
        virtual const char *shortDescription() const
        {
            return "Print help information";
        }

        virtual int run(int argc, char *argv[]);

    private:
        const CommandLineModuleManager &manager_;

        GMX_DISALLOW_COPY_AND_ASSIGN(CommandLineHelpModule);
};

CommandLineHelpModule::CommandLineHelpModule(const CommandLineModuleManager &manager)
    : manager_(manager)
{
}

int CommandLineHelpModule::run(int argc, char *argv[])
{
    manager_.impl_->printUsage(true);
    return 0;
}

} // namespace internal


/********************************************************************
 * CommandLineModuleManager
 */

CommandLineModuleManager::CommandLineModuleManager(const char *realBinaryName)
    : impl_(new Impl(realBinaryName))
{
    addModule(CommandLineModulePointer(new internal::CommandLineHelpModule(*this)));
}

CommandLineModuleManager::~CommandLineModuleManager()
{
}

void CommandLineModuleManager::addModule(CommandLineModulePointer module)
{
    GMX_ASSERT(impl_->modules_.find(module->name()) == impl_->modules_.end(),
               "Attempted to register a duplicate module name");
    impl_->modules_.insert(std::make_pair(std::string(module->name()),
                                          move(module)));
}

int CommandLineModuleManager::run(int argc, char *argv[])
{
    int argOffset = 0;
    Impl::ModuleMap::const_iterator module
        = impl_->findModuleFromBinaryName(argv[0]);
    if (module == impl_->modules_.end())
    {
        if (argc < 2)
        {
            fprintf(stderr, "\n");
            impl_->printUsage(false);
            return 2;
        }
        module = impl_->findModuleByName(argv[1]);
        argOffset = 1;
    }
    if (module == impl_->modules_.end())
    {
        fprintf(stderr, "\n");
        fprintf(stderr, "Unknown command: '%s'\n", argv[1]);
        impl_->printUsage(true);
        return 2;
    }
    return module->second->run(argc - argOffset, argv + argOffset);
}

} // namespace gmx
