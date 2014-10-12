/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2016, by the GROMACS development team, led by
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
 * Implements registerTrajectoryAnalysisModules().
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "modules.h"

#include "gromacs/commandline/cmdlinemodulemanager.h"
#include "gromacs/trajectoryanalysis/analysismodule.h"
#include "gromacs/trajectoryanalysis/cmdlinerunner.h"
#include "gromacs/utility/exceptions.h"

#include "modules/angle.h"
#include "modules/distance.h"
#include "modules/freevolume.h"
#include "modules/pairdist.h"
#include "modules/rdf.h"
#include "modules/sasa.h"
#include "modules/select.h"

#include <string>
#include <map>

namespace gmx
{

namespace
{

using namespace gmx::analysismodules;
//! \brief Storage for module factories
const std::pair < std::string, std::pair < const char *, TrajectoryAnalysisCommandLineRunner::ModuleFactoryMethod>>
registeredModules[] {
    {
        AngleInfo::name, {
            AngleInfo::shortDescription, AngleInfo::create
        }
    },
    {
        DistanceInfo::name, {
            DistanceInfo::shortDescription, DistanceInfo::create
        }
    },
    {
        FreeVolumeInfo::name, {
            FreeVolumeInfo::shortDescription, FreeVolumeInfo::create
        }
    },
    {
        PairDistanceInfo::name, {
            PairDistanceInfo::shortDescription, PairDistanceInfo::create
        }
    },
    {
        RdfInfo::name, {
            RdfInfo::shortDescription, RdfInfo::create
        }
    },
    {
        SasaInfo::name, {
            SasaInfo::shortDescription, SasaInfo::create
        }
    },
    {
        SelectInfo::name, {
            SelectInfo::shortDescription, SelectInfo::create
        }
    },
};
}   // namespace

//! \cond libapi
void registerTrajectoryAnalysisModules(CommandLineModuleManager *manager)
{
    CommandLineModuleGroup group = manager->addModuleGroup("Trajectory analysis");

    for (const auto &item : registeredModules)
    {
        TrajectoryAnalysisCommandLineRunner::registerModule(
                manager, item.first.c_str(), item.second.first, item.second.second);
        group.addModule(item.first.c_str());
    }
}

std::unique_ptr<TrajectoryAnalysisModule>
createTrajectoryAnalysisModuleByName(const char *name)
{
    for (const auto &item : registeredModules)
    {
        if (item.first == name)
        {
            return item.second.second();
        }
    }

    GMX_THROW(APIError("Module has not been registered"));
}
//! \endcond

} // namespace gmx
