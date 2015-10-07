/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements classes in modules.h.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_energyanalysis
 */
#include "gmxpre.h"

#include "modules.h"

#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"

#include "dhdl.h"
#include "fluctprops.h"
#include "simple.h"
#include "viscosity.h"

namespace gmx
{

namespace energyanalysis
{

namespace
{
/*! \brief
 * Convenience method for registering a command-line module for energy
 * analysis.
 *
 * \tparam ModuleInfo  Info about energy analysis module to wrap.
 *
 * \p ModuleInfo should have static public members
 * `const char name[]`, `const char shortDescription[]`, and
 * `gmx::EnergyAnalysisModulePointer create()`.
 *
 * \ingroup module_energyanalysis
 */
template <class ModuleInfo>
void registerModule(CommandLineModuleManager *manager,
                    CommandLineModuleGroup    group)
{
    ICommandLineOptionsModule::registerModuleFactory(manager,
                                                     ModuleInfo::name,
                                                     ModuleInfo::shortDescription,
                                                     &ModuleInfo::create);
    group.addModule(ModuleInfo::name);
}

}       // namespace

//! \cond libapi
void registerEnergyAnalysisModules(CommandLineModuleManager *manager)
{
    CommandLineModuleGroup group = manager->addModuleGroup("Energy analysis");
    registerModule<SimpleInfo>(manager, group);
    registerModule<DhdlInfo>(manager, group);
    registerModule<FluctPropsInfo>(manager, group);
    registerModule<ViscosityInfo>(manager, group);
}
//! \endcond

} // namespace energyanalysis

} // namespace gmx
