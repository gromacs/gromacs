/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/energyanalysis/energyanalysis.h"
#include "gromacs/energyanalysis/modules.h"
#include "gromacs/energyanalysis/viscosity.h"
#include "gromacs/energyanalysis/fluctprops.h"
#include "gromacs/energyanalysis/dhdl.h"
#include "gromacs/energyanalysis/freeenergydifference.h"
#include "gromacs/energyanalysis/simple.h"
#include "gromacs/energyanalysis/handler.h"

namespace gmx
{

/*! \brief
 * Factory method type for creating an energy analysis module.
 *
 * This method allows the module creation to be postponed to be inside
 * the try/catch block in runAsMain()/registerModule() implementation
 * methods and still keep the implementation out of the header, making
 * the ABI more stable.
 */
typedef EnergyAnalysisPointer (*ModuleFactoryMethod)();

/*! \brief
 * Class to interface between the energy handler and the command line tools
 */
class EnergyHandlerInterface : public CommandLineModuleInterface
{
    public:
        /*! \brief
         * Constructs a module.
         *
         * \param[in] name         Name for the module.
         * \param[in] description  One-line description for the module.
         * \param[in] factory      Factory method to create the analysis module.
         *
         * Does not throw.  This is important for correct implementation of
         * runAsMain().
         * TODO: check that the above holds!
         */
        EnergyHandlerInterface(const char *name, const char *description,
                               ModuleFactoryMethod factory)
            : name_(name), description_(description), factory_(factory)
        {
        }

        //! Return the name of this module
        virtual const char *name() const { return name_; }

        //! Return a short description as seen with gmx help
        virtual const char *shortDescription() const { return description_; };

        /*! \brief
         * Run the modules
         *
         * \param[in] argc Number of command line arguments
         * \param[in] argv The command line arguments
         * \return 0 if all OK , 1 otherwise
         */
        virtual int run(int argc, char *argv[]);

        /*! \brief
         * Print a help text on the terminal
         * \param[in] context Something I do not quite understand
         */
        virtual void writeHelp(const CommandLineHelpContext &context) const;

    private:
        //! Name corresponding to the module
        const char             *name_;
        //! Short description (see above)
        const char             *description_;
        //! Function pointer that creates an instance of this module
        ModuleFactoryMethod     factory_;
        //! Hocus pocus
        GMX_DISALLOW_COPY_AND_ASSIGN(EnergyHandlerInterface);
};

class EnergyAnalysisCommandLineRunner
{
    public:
        /*! \brief
         * Implements a main() method that runs a given module.
         *
         * \tparam ModuleType  Energy analysis module.
         * \param  argc        \c argc passed to main().
         * \param  argv        \c argv passed to main().
         *
         * This method abstracts away all the logic required to implement a
         * main() method in user tools, allowing that to be changed without
         * requiring changes to the tools themselves.
         *
         * \p ModuleType should be default-constructible and derive from
         * Energyanalysis.
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
         * \tparam ModuleType  Energy analysis module.
         * \param  manager     Manager to register the module to.
         * \param  name        Name of the module to register.
         * \param  description One-line description for the module to register.
         *
         * \p ModuleType should be default-constructible and derive from
         * Energyanalysis.
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
         * \tparam ModuleType  Energy analysis module.
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
        EnergyAnalysisCommandLineRunner(EnergyAnalysis *module);
        ~EnergyAnalysisCommandLineRunner();

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
        static EnergyAnalysisPointer createModule()
        {
            return EnergyAnalysisPointer(new ModuleType());
        }

        //! Implements the template runAsMain() method.
        static int runAsMain(int argc, char *argv[],
                             ModuleFactoryMethod factory);
};

void EnergyAnalysisCommandLineRunner::registerModule(CommandLineModuleManager *manager,
                                                     const char               *name,
                                                     const char               *description,
                                                     ModuleFactoryMethod       factory)
{
    //    CommandLineModulePointer module(
    //      new Impl::RunnerCommandLineModule(name, description, factory));
    printf("Hallo!\n");
}


/*! \brief
 * Convenience method for registering a command-line module for energy
 * analysis.
 *
 * \tparam ModuleInfo  Info about energy analysis module to wrap.
 *
 * \p ModuleInfo should have static public members
 * `const char name[]`, `const char shortDescription[]`, and
 * `gmx::TrajectoryAnalysisModulePointer create()`.
 *
 * \ingroup module_energyyanalysis
 */
template <class ModuleInfo>
void registerModule(CommandLineModuleManager *manager,
                    CommandLineModuleGroup    group)
{
    EnergyAnalysisCommandLineRunner::registerModule(manager,
                                                    ModuleInfo::name,
                                                    ModuleInfo::shortDescription,
                                                    &ModuleInfo::create);
    // Need to add module before adding it to group
    group.addModule(ModuleInfo::name);
}

//! \cond libapi
void registerEnergyAnalysisModules(CommandLineModuleManager *manager)
{
    CommandLineModuleGroup group = manager->addModuleGroup("Energy analysis");
    registerModule<SimpleInfo>(manager, group);
    registerModule<DhdlInfo>(manager, group);
    registerModule<FreeEnergyDifferenceInfo>(manager, group);
    registerModule<FluctPropsInfo>(manager, group);
    registerModule<ViscosityInfo>(manager, group);
}
//! \endcond

int EnergyHandlerInterface::run(int   argc,
                                char *argv[])
{
    // Energy handler utility
    EnergyHandler         eh;
    EnergyAnalysisPointer ptr = factory_();
    eh.addAnalysisTool(ptr);

    eh.prepare(&argc, argv);

    return eh.readFiles();
}

void EnergyHandlerInterface::writeHelp(const CommandLineHelpContext gmx_unused &context) const
{
    printf("Please implement proper help writing\n");
}

}
