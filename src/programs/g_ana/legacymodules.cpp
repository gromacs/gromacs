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
 * Implements command-line modules for pre-5.0 binaries.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
#include "legacymodules.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/onlinehelp/helpwritercontext.h"
#include "gromacs/utility/exceptions.h"

#include "gromacs/legacyheaders/gmx_ana.h"

namespace
{

/*! \internal \brief
 * Command-line module for wrapping pre-5.0 binaries.
 *
 * Implements a gmx::CommandLineModuleInterface, given a function with
 * C/C++ main signature.
 */
class LegacyCmdLineWrapper : public gmx::CommandLineModuleInterface
{
    public:
        //! Function pointer type for the main function of the module.
        typedef int (*MainFunction)(int argc, char *argv[]);

        /*! \brief
         * Convenience function for creating and registering a module.
         *
         * \param[in] manager  Module manager to which to register the module.
         * \param[in] main     Main function to wrap.
         * \param[in] name     Name for the new module.
         * \param[in] shortDescription One-line description for the new module.
         */
        static void registerModule(gmx::CommandLineModuleManager *manager,
                                   MainFunction main, const char *name,
                                   const char *shortDescription)
        {
            gmx::CommandLineModulePointer module(
                    new LegacyCmdLineWrapper(main, name, shortDescription));
            manager->addModule(gmx::move(module));
        }

        /*! \brief
         * Creates a wrapper module for the given main function.
         *
         * \see registerModule()
         */
        LegacyCmdLineWrapper(MainFunction main, const char *name,
                             const char *shortDescription)
            : main_(main), name_(name), shortDescription_(shortDescription)
        {
        }

        virtual const char *name() const
        {
            return name_;
        }
        virtual const char *shortDescription() const
        {
            return shortDescription_;
        }

        virtual int run(int argc, char *argv[]);
        virtual void writeHelp(const gmx::HelpWriterContext &context) const;

    private:
        MainFunction            main_;
        const char             *name_;
        const char             *shortDescription_;

};

int LegacyCmdLineWrapper::run(int argc, char *argv[])
{
    return main_(argc, argv);
}

void LegacyCmdLineWrapper::writeHelp(const gmx::HelpWriterContext &context) const
{
    if (context.outputFormat() != gmx::eHelpOutputFormat_Console)
    {
        GMX_THROW(gmx::NotImplementedError(
                    "Command-line help is not implemented for this output format"));
    }
    char *argv[2];
    // TODO: The constness should not be cast away.
    argv[0] = const_cast<char *>(name_);
    argv[1] = const_cast<char *>("-h");
    main_(2, argv);
}

} // namespace

void registerLegacyModules(gmx::CommandLineModuleManager *manager)
{
    LegacyCmdLineWrapper::registerModule(manager, &gmx_dist, "dist",
            "Calculates distances between centers of mass of two groups");
}
