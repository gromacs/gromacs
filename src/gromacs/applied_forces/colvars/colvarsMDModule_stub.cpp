/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2023- The GROMACS Authors
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
 * Stub implementation of ColvarsMDModule
 * Compiled in case if Colvars is not activated with CMake.
 *
 * \author Hubert Santuz <hubert.santuz@gmail.com>
 * \ingroup module_applied_forces
 */
#include "gmxpre.h"

#include <memory>
#include <string>

#include "gromacs/domdec/localatomsetmanager.h"
#include "gromacs/fileio/checkpoint.h"
#include "gromacs/mdtypes/imdmodule.h"
#include "gromacs/mdtypes/imdpoptionprovider.h"
#include "gromacs/utility/keyvaluetreebuilder.h"

#include "colvarsMDModule.h"
#include "colvarsoptions.h"


namespace gmx
{

namespace
{

/*! \internal
 * \brief Colvars module
 *
 * Stub Implementation in case Colvars library is not compiled
 */
class ColvarsMDModule final : public IMDModule
{
public:
    //! \brief Construct the colvars module.
    explicit ColvarsMDModule() = default;


    void subscribeToPreProcessingNotifications(MDModulesNotifiers* /*notifier*/) override
    {
        // Proper exit when Colvars is activated from the mdp but has not been compiled along with GROMACS.
        // In this case, Colvars config file cannot be parsed & verified and the tpr created could be flawed.
        if (colvarsOptionsStub_.isActive())
        {
            GMX_THROW(InternalError(
                    "Colvars module is activated but GROMACS has not been compiled with Colvars, "
                    "Colvars simulation is not possible.\n"
                    "Please, reconfigure GROMACS with -DGMX_USE_COLVARS=internal.\n"));
        }
    }

    void subscribeToSimulationSetupNotifications(MDModulesNotifiers* /*notifier*/) override {}

    void subscribeToSimulationRunNotifications(MDModulesNotifiers* /*notifier*/) override {}

    IMdpOptionProvider* mdpOptionProvider() override { return &colvarsOptionsStub_; }

    IMDOutputProvider* outputProvider() override { return nullptr; }

    void initForceProviders(ForceProviders* /*forceProviders*/) override
    {
        // Proper exit when Colvars is activated from the tpr but has not been compiled along with GROMACS.
        if (colvarsOptionsStub_.isActive())
        {
            GMX_THROW(InternalError(
                    "Colvars module is activated but GROMACS has not been compiled with Colvars, "
                    "Colvars simulation is not possible.\n"
                    "Please, reconfigure GROMACS with -DGMX_USE_COLVARS=internal.\n"));
        }
    }

private:
    ColvarsOptions colvarsOptionsStub_;
};

} // namespace

std::unique_ptr<IMDModule> ColvarsModuleInfo::create()
{
    return std::make_unique<ColvarsMDModule>();
}

std::string colvarsDescription()
{
    return "disabled";
}

} // namespace gmx
