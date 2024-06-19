/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2014- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements supporting routines for gmx::ICommandLineOptionsModule.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlineoptionsmodule.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/options/behaviorcollection.h"
#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/options/ioptionsbehavior.h"
#include "gromacs/options/options.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{
class CommandLineHelpContext;

namespace
{

/********************************************************************
 * CommandLineOptionsModuleSettings
 */

class CommandLineOptionsModuleSettings : public ICommandLineOptionsModuleSettings
{
public:
    explicit CommandLineOptionsModuleSettings(OptionsBehaviorCollection* behaviors) :
        behaviors_(*behaviors)
    {
    }

    const std::string& helpText() const { return helpText_; }

    ArrayRef<const std::string> bugText() const { return bugText_; }

    void setHelpText(const ArrayRef<const char* const>& help) override
    {
        helpText_ = joinStrings(help, "\n");
    }

    void setBugText(const ArrayRef<const char* const>& bug) override
    {
        bugText_ = std::vector<std::string>(bug.begin(), bug.end());
    }
    void addOptionsBehavior(const OptionsBehaviorPointer& behavior) override
    {
        behaviors_.addBehavior(behavior);
    }

private:
    std::string                helpText_;
    std::vector<std::string>   bugText_;
    OptionsBehaviorCollection& behaviors_;
};

/********************************************************************
 * CommandLineOptionsModule
 */

class CommandLineOptionsModule : public ICommandLineModule
{
public:
    //! Shorthand for the factory function pointer type.
    typedef ICommandLineOptionsModule::FactoryMethod FactoryMethod;

    CommandLineOptionsModule(const char* name, const char* description, FactoryMethod factory) :
        name_(name), description_(description), factory_(std::move(factory))
    {
    }
    CommandLineOptionsModule(const char* name, const char* description, ICommandLineOptionsModulePointer module) :
        name_(name), description_(description), module_(std::move(module))
    {
    }
    const char* name() const override { return name_; }
    const char* shortDescription() const override { return description_; }

    void init(CommandLineModuleSettings* settings) override;
    int  run(int argc, char* argv[]) override;
    void writeHelp(const CommandLineHelpContext& context) const override;

private:
    void parseOptions(int argc, char* argv[]);

    const char*                      name_;
    const char*                      description_;
    FactoryMethod                    factory_;
    ICommandLineOptionsModulePointer module_;
};

void CommandLineOptionsModule::init(CommandLineModuleSettings* settings)
{
    if (!module_)
    {
        GMX_RELEASE_ASSERT(factory_ != nullptr, "Neither factory nor module provided");
        module_ = factory_();
    }
    module_->init(settings);
}

int CommandLineOptionsModule::run(int argc, char* argv[])
{
    GMX_RELEASE_ASSERT(module_, "init() has not been called");
    parseOptions(argc, argv);
    return module_->run();
}

void CommandLineOptionsModule::writeHelp(const CommandLineHelpContext& context) const
{
    ICommandLineOptionsModulePointer moduleGuard;
    ICommandLineOptionsModule*       module = module_.get();
    if (!module)
    {
        GMX_RELEASE_ASSERT(factory_ != nullptr, "Neither factory nor module provided");
        moduleGuard = factory_();
        module      = moduleGuard.get();
    }
    Options                          options;
    OptionsBehaviorCollection        behaviors(&options);
    CommandLineOptionsModuleSettings settings(&behaviors);
    module->initOptions(&options, &settings);
    CommandLineHelpWriter(options)
            .setHelpText(settings.helpText())
            .setKnownIssues(settings.bugText())
            .writeHelp(context);
}

void CommandLineOptionsModule::parseOptions(int argc, char* argv[])
{
    FileNameOptionManager fileoptManager;
    Options               options;

    options.addManager(&fileoptManager);

    OptionsBehaviorCollection        behaviors(&options);
    CommandLineOptionsModuleSettings settings(&behaviors);
    module_->initOptions(&options, &settings);
    {
        CommandLineParser parser(&options);
        parser.parse(&argc, argv);
        behaviors.optionsFinishing();
        options.finish();
    }
    module_->optionsFinished();
    behaviors.optionsFinished();
}

} // namespace

/********************************************************************
 * ICommandLineOptionsModuleSettings
 */

ICommandLineOptionsModuleSettings::~ICommandLineOptionsModuleSettings() {}

/********************************************************************
 * ICommandLineOptionsModule
 */

ICommandLineOptionsModule::~ICommandLineOptionsModule() {}

// static
std::unique_ptr<ICommandLineModule> ICommandLineOptionsModule::createModule(const char* name,
                                                                            const char* description,
                                                                            ICommandLineOptionsModulePointer module)
{
    return std::unique_ptr<ICommandLineModule>(
            new CommandLineOptionsModule(name, description, std::move(module)));
}

// static
int ICommandLineOptionsModule::runAsMain(int           argc,
                                         char*         argv[],
                                         const char*   name,
                                         const char*   description,
                                         FactoryMethod factory)
{
    CommandLineOptionsModule module(name, description, std::move(factory));
    return CommandLineModuleManager::runAsMainSingleModule(argc, argv, &module);
}

// static
void ICommandLineOptionsModule::registerModuleFactory(CommandLineModuleManager* manager,
                                                      const char*               name,
                                                      const char*               description,
                                                      FactoryMethod             factory)
{
    CommandLineModulePointer module(new CommandLineOptionsModule(name, description, std::move(factory)));
    manager->addModule(std::move(module));
}

// static
void ICommandLineOptionsModule::registerModuleDirect(CommandLineModuleManager*        manager,
                                                     const char*                      name,
                                                     const char*                      description,
                                                     ICommandLineOptionsModulePointer module)
{
    CommandLineModulePointer wrapperModule(createModule(name, description, std::move(module)));
    manager->addModule(std::move(wrapperModule));
}

} // namespace gmx
