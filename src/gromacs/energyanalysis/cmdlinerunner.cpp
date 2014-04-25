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
#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/selection/selectioncollection.h"
#include "gromacs/selection/selectionoptionmanager.h"
#include "gromacs/utility/file.h"
#include "analysismodule.h"
#include "cmdlinerunner.h"
#include "handler.h"

namespace gmx
{

/********************************************************************
 * EnergyAnalysisCommandLineRunner::Impl
 */

class EnergyAnalysisCommandLineRunner::Impl
{
    public:
        class RunnerCommandLineModule;

        Impl(EnergyAnalysisModule *module);
        ~Impl();

        void parseOptions(SelectionCollection *selections,
                          int *argc, char *argv[]);

        EnergyAnalysisModule *module_;
};


EnergyAnalysisCommandLineRunner::Impl::Impl(
        EnergyAnalysisModule *module)
    : module_(module)
{
}


EnergyAnalysisCommandLineRunner::Impl::~Impl()
{
}


void
EnergyAnalysisCommandLineRunner::Impl::parseOptions(
        //EnergyAnalysisSettings *settings,
        //EnergyAnalysisRunnerCommon *common,
        SelectionCollection *selections,
        int *argc, char *argv[])
{
    FileNameOptionManager  fileoptManager;
    SelectionOptionManager seloptManager(selections);
    Options                options(NULL, NULL);
    Options                moduleOptions(module_->name(), module_->description());
    //Options                commonOptions("common", "Common analysis control");
    Options                selectionOptions("selection", "Common selection control");

    options.addManager(&fileoptManager);
    options.addManager(&seloptManager);
    //options.addSubSection(&commonOptions);
    //options.addSubSection(&selectionOptions);
    options.addSubSection(&moduleOptions);

    module_->initOptions(&moduleOptions);
    //common->initOptions(&commonOptions);
    selections->initOptions(&selectionOptions);

    {
        CommandLineParser  parser(&options);
        // TODO: Print the help if user provides an invalid option?
        // Or just add a message advising the user to invoke the help?
        parser.parse(argc, argv);
        //common->scaleTimeOptions(&options);
        options.finish();
    }

    //common->optionsFinished(&commonOptions);
    //module_->optionsFinished(&moduleOptions, settings);

    //common->initIndexGroups(selections, bUseDefaultGroups_);

    const bool bInteractive = File::standardInput().isInteractive();
    seloptManager.parseRequestedFromStdin(bInteractive);
    //common->doneIndexGroups(selections);

    //common->initTopology(selections);
    selections->compile();
}


/*! \brief
 * Class to interface between the energy handler and the command line tools
 */
class EnergyAnalysisCommandLineRunner::Impl::RunnerCommandLineModule
    : public CommandLineModuleInterface
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
        RunnerCommandLineModule(const char *name, const char *description,
                                ModuleFactoryMethod factory)
            : name_(name), description_(description), factory_(factory)
        {
        }
        //! Destructor
        virtual ~RunnerCommandLineModule();

        //! Return the name of this module
        const char *name() const { return name_; }

        //! Return a short description as seen with gmx help
        const char *shortDescription() const { return description_; };

        /*! \brief
         * Initializes the module and provides settings for the runner.
         *
         * This will be called before run(), and can be used to adjust
         * initialization that the runner does.
         */
        virtual void init(CommandLineModuleSettings *settings);
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
        GMX_DISALLOW_COPY_AND_ASSIGN(RunnerCommandLineModule);
};

void EnergyAnalysisCommandLineRunner::Impl::RunnerCommandLineModule::init(CommandLineModuleSettings * /*settings*/)
{

}

void EnergyAnalysisCommandLineRunner::registerModule(CommandLineModuleManager *manager,
                                                     const char               *name,
                                                     const char               *description,
                                                     ModuleFactoryMethod       factory)
{
    CommandLineModulePointer module
        (new Impl::RunnerCommandLineModule(name, description, factory));
    manager->addModule(move(module));
}

EnergyAnalysisCommandLineRunner::EnergyAnalysisCommandLineRunner(EnergyAnalysisModule *module)
    : impl_(new Impl(module))
{

}

EnergyAnalysisCommandLineRunner::~EnergyAnalysisCommandLineRunner()
{
}

EnergyAnalysisCommandLineRunner::Impl::RunnerCommandLineModule::~RunnerCommandLineModule()
{
}

void EnergyAnalysisCommandLineRunner::Impl::RunnerCommandLineModule::writeHelp(const CommandLineHelpContext &context) const
{
    // TODO: This method duplicates some code from run().
    // See how to best refactor it to share the common code.
    SelectionCollection             selections;
    //EnergyAnalysisSettings      settings;
    //EnergyAnalysisRunnerCommon  common(&settings);

    SelectionOptionManager          seloptManager(&selections);
    Options                         options(NULL, NULL);
    Options                         moduleOptions(name(), shortDescription());
    //Options                         commonOptions("common", "Common analysis control");
    Options                         selectionOptions("selection", "Common selection control");

    options.addManager(&seloptManager);
    //options.addSubSection(&commonOptions);
    options.addSubSection(&selectionOptions);
    options.addSubSection(&moduleOptions);

    factory_()->initOptions(&moduleOptions);
    //common.initOptions(&commonOptions);
    selections.initOptions(&selectionOptions);

    CommandLineHelpWriter(options)
        .setShowDescriptions(true)
        .setTimeUnitString("ps") //settings.timeUnitManager().timeUnitAsString())
        .writeHelp(context);
}

int EnergyAnalysisCommandLineRunner::Impl::RunnerCommandLineModule::run(int   argc,
                                                                        char *argv[])
{
    // Energy handler utility
    EnergyHandler         eh;
    //EnergyAnalysisModulePointer ptr = (EnergyAnalysisModulePointer)impl_->module_;
    eh.addAnalysisTool(factory_());

    eh.prepare(&argc, argv);

    return eh.readFiles();
}

void EnergyAnalysisCommandLineRunner::writeHelp(const CommandLineHelpContext gmx_unused &context)
{
    printf("Fix me %s %d\n", __FILE__, __LINE__);
    //impl_->module_->writeHelp(context);
}

}
