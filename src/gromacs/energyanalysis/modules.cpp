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

#include "gromacs/energyanalysis/modules.h"
#include "gromacs/energyanalysis/viscosity.h"
#include "gromacs/energyanalysis/fluctprops.h"
#include "gromacs/energyanalysis/dhdl.h"
#include "gromacs/energyanalysis/freeenergydifference.h"
#include "gromacs/energyanalysis/simple.h"
#include "gromacs/energyanalysis/handler.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"

namespace gmx
{

typedef EnergyAnalysisPointer (*ModuleFactoryMethod)();

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
         */
        EnergyHandlerInterface(const char *name, const char *description,
                               ModuleFactoryMethod factory)
            : name_(name), description_(description), factory_(factory)
        {
        }

        virtual const char *name() const { return name_; }
        virtual const char *shortDescription() const { return description_; };

        virtual int run(int argc, char *argv[]);
        virtual void writeHelp(const CommandLineHelpContext &context) const;

    private:
        const char             *name_;
        const char             *description_;
        ModuleFactoryMethod     factory_;

        GMX_DISALLOW_COPY_AND_ASSIGN(EnergyHandlerInterface);
};

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
    CommandLineModulePointer module(
            new EnergyHandlerInterface(ModuleInfo::name,
                                       ModuleInfo::shortDescription,
                                       &ModuleInfo::create));
    manager->addModule(move(module));
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

void EnergyHandlerInterface::writeHelp(const CommandLineHelpContext &context) const
{
}

}
