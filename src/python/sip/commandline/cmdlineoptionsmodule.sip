/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016 by the GROMACS development team, led by
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

class ICommandLineOptionsModuleSettings /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/utility/arrayref.h>
#include <gromacs/commandline/cmdlineoptionsmodule.h>
using gmx::ICommandLineOptionsModuleSettings;
%End
public:
    virtual void setHelpText(const char *help) = 0;
    %MethodCode
        const char *const helpText[] = {a0};
        sipCpp->setHelpText(helpText);
    %End
    virtual void addOptionsBehavior(
        IOptionsBehavior *behavior /Transfer/) = 0;
    %MethodCode
        sipCpp->addOptionsBehavior(std::shared_ptr<IOptionsBehavior>(a0));
    %End
private:
    virtual ~ICommandLineOptionsModuleSettings();
};

// This class is needed by Pipeline python implementation.
// ICommandLineOptionsModuleSettings can't be subclassed in Python
// due to private destructor, so subclass it here.

%ModuleHeaderCode
#include <gromacs/commandline/cmdlineoptionsmodule.h>
#include <memory>

namespace pygmx {
    class DummyOptionsModuleSettings : public gmx::ICommandLineOptionsModuleSettings {
    public:
        virtual void setHelpText(const gmx::ConstArrayRef<const char*>&) {};
        virtual void addOptionsBehavior(const std::shared_ptr<gmx::IOptionsBehavior>&) {};
    };
}
%End

class DummyOptionsModuleSettings : public ICommandLineOptionsModuleSettings {
%TypeHeaderCode
#include <gromacs/commandline/cmdlineoptionsmodule.h>
using gmx::ICommandLineOptionsModuleSettings;
using pygmx::DummyOptionsModuleSettings;
%End
private:
    virtual void setHelpText(const char *);
    virtual void addOptionsBehavior(IOptionsBehavior* /Transfer/);
};

class CommandLineModuleSettings /NoDefaultCtors/ {
%TypeHeaderCode
#include <gromacs/commandline/cmdlinemodule.h>
using gmx::CommandLineModuleSettings;
%End
};

class ICommandLineOptionsModule /VirtualErrorHandler=vehandler/ {
%TypeHeaderCode
#include <gromacs/commandline/cmdlineoptionsmodule.h>
using gmx::ICommandLineOptionsModule;
%End
public:
    virtual void init(CommandLineModuleSettings *settings) = 0;
    virtual void initOptions(IOptionsContainer *options, ICommandLineOptionsModuleSettings *settings) = 0;
    virtual void optionsFinished() = 0;
    virtual int run() = 0;
    static int runAsMain(SIP_PYOBJECT mod, SIP_PYLIST argv, const char *name, const char *description);
    %MethodCode
        auto runner = [a0, a2, a3] (int argc, char **argv) -> int {
            return ICommandLineOptionsModule::runAsMain(argc, argv, a2, a3, PyFactory<gmx::ICommandLineOptionsModule>(sipType_ICommandLineOptionsModule, a0));
        };
        sipRes = argv_wrapper(a1, runner);
    %End
};
