/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2015, by the GROMACS development team, led by
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
 *
 * \brief This file defines gmx::MpiHandler that provides behaviour
 * that depends on MPI functionality
 *
 * \todo Gradually move data and functionality from t_commrec here.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 *
 * \ingroup module_utility
 */

#include "gmxpre.h"

#include "mpihandler.h"

#include <stdio.h>

#include "gromacs/legacyheaders/main.h"
#include "gromacs/legacyheaders/network.h"
#include "gromacs/legacyheaders/types/commrec.h"

namespace gmx
{

MpiHandler::MpiHandler(const t_commrec *cr)
    : cr_(cr)
{
}

MpiHandler::~MpiHandler()
{
}

bool
MpiHandler::hasMultipleRanks() const
{
    return PAR(cr_);
}

bool
MpiHandler::isMaster() const
{
    return MASTER(cr_);
}

bool
MpiHandler::isSimMaster() const
{
    return SIMMASTER(cr_);
}

bool
MpiHandler::isMultiMaster() const
{
    return MULTIMASTER(cr_);
}

bool
MpiHandler::isMultiSim() const
{
    return MULTISIM(cr_);
}

void
MpiHandler::broadcast(int byteSize, void *data) const
{
    gmx_bcast(byteSize, data, cr_);
}

void
MpiHandler::checkAcrossMultiSim(FILE *fp, int theInteger, const char *description, bool bQuiet) const
{
    if (isMultiSim() && isMaster())
    {
        check_multi_int(fp, cr_->ms, theInteger, description, bQuiet);
    }
}

} // namespace
