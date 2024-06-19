/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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
 * Implements classes in analysissettings.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "gromacs/trajectoryanalysis/analysissettings.h"

#include <memory>

#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/gmxassert.h"

#include "analysissettings_impl.h"

namespace gmx
{
class AnalysisDataPlotSettings;
enum class TimeUnit : int;


/********************************************************************
 * TrajectoryAnalysisSettings
 */

TrajectoryAnalysisSettings::TrajectoryAnalysisSettings() : impl_(new Impl)
{
    impl_->frflags |= TRX_NEED_X;
}


TrajectoryAnalysisSettings::~TrajectoryAnalysisSettings() {}


void TrajectoryAnalysisSettings::setOptionsModuleSettings(ICommandLineOptionsModuleSettings* settings)
{
    impl_->optionsModuleSettings_ = settings;
}


TimeUnit TrajectoryAnalysisSettings::timeUnit() const
{
    return impl_->timeUnit;
}


const AnalysisDataPlotSettings& TrajectoryAnalysisSettings::plotSettings() const
{
    return impl_->plotSettings;
}


unsigned long TrajectoryAnalysisSettings::flags() const
{
    return impl_->flags;
}


bool TrajectoryAnalysisSettings::hasFlag(unsigned long flag) const
{
    return (impl_->flags & flag) != 0U;
}


bool TrajectoryAnalysisSettings::hasPBC() const
{
    return impl_->bPBC;
}


bool TrajectoryAnalysisSettings::hasRmPBC() const
{
    return impl_->bRmPBC;
}


int TrajectoryAnalysisSettings::frflags() const
{
    return impl_->frflags;
}


void TrajectoryAnalysisSettings::setFlags(unsigned long flags)
{
    impl_->flags = flags;
}


void TrajectoryAnalysisSettings::setFlag(unsigned long flag, bool bSet)
{
    if (bSet)
    {
        impl_->flags |= flag;
    }
    else
    {
        impl_->flags &= ~flag;
    }
}


void TrajectoryAnalysisSettings::setPBC(bool bPBC)
{
    impl_->bPBC = bPBC;
}


void TrajectoryAnalysisSettings::setRmPBC(bool bRmPBC)
{
    impl_->bRmPBC = bRmPBC;
}


void TrajectoryAnalysisSettings::setFrameFlags(int frflags)
{
    impl_->frflags = frflags;
}

void TrajectoryAnalysisSettings::setHelpText(const ArrayRef<const char* const>& help)
{
    GMX_RELEASE_ASSERT(impl_->optionsModuleSettings_ != nullptr,
                       "setHelpText() called in invalid context");
    impl_->optionsModuleSettings_->setHelpText(help);
}

} // namespace gmx
