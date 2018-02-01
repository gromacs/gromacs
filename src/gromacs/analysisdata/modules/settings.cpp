/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014,2015, by the GROMACS development team, led by
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
 * Implements classes in settings.h to configure trajectory data writing into a file.
 *
 * \ingroup module_analysisdata
 * \author 
 */
#include "gmxpre.h"

#include "settings.h"

#include <memory>
#include <string>

#include "gromacs/analysisdata/datamodule.h"
#include "gromacs/fileio/filetypes.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/timeunitmanager.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/classhelpers.h"
#include "gromacs/utility/real.h"

namespace gmx
{

class SelectionCollection;
class Selection;

/*! \brief
 * Common settings for data plots.
 *
 */
TrajectoryDataWriteSettings::TrajectoryDataWriteSettings() : selections_(nullptr),
    sel_(nullptr), filetype_(efNR), mtop_(nullptr), top_(nullptr), bgenCon_(false), prec_(1000),
    bPrec_(false), time_(0), bTime_(false), bVel_(false), bForce_(false), bAtoms_(false)
{
    init_atom(&atoms_);
}

TrajectoryDataWriteSettings::~TrajectoryDataWriteSettings()
{
    if (bgenCon_)
    {
        gmx_conect_done(connections_);
    }
    if (bAtoms_)
    {
        done_atom(&atoms_);
    }
}


void
TrajectoryDataWriteSettings::initOptions(IOptionsContainer *options)
{
    options->addOption(BooleanOption("vel").store(&bVel_)
            .description("Write velocities to file if possible"));
    options->addOption(BooleanOption("force").store(&bForce_)
            .description("Write forces to file if possible"));
    options->addOption(BooleanOption("conect").store(&bgenCon_)
            .description("Write connection information to file"));
    options->addOption(DoubleOption("prec").store(&prec_)
            .description("Output precision for compressed files")
            .storeIsSet(&bPrec_));
    options->addOption(DoubleOption("t0").store(&time_)
            .description("Set custom start time")
            .storeIsSet(&bTime_));
}

} // namespace gmx

