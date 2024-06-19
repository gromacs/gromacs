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
 * Implements helper class for autocorrelation tests
 *
 * \author Anders G&auml;rden&auml;s <anders.gardenas@gmail.com>
 * \ingroup module_correlationfunctions
 */
#include "gmxpre.h"

#include "correlationdataset.h"

#include <cmath>

#include <filesystem>
#include <sstream>

#include "gromacs/fileio/xvgr.h"
#include "gromacs/utility/smalloc.h"

#include "testutils/testfilemanager.h"

CorrelationDataSet::CorrelationDataSet(const std::string& fileName)
{
    std::filesystem::path fileNm = gmx::test::TestFileManager::getInputFilePath(fileName);
    nrLines_                     = read_xvg(fileNm, &tempValues_, &nrColumns_);

    dt_        = tempValues_[0][1] - tempValues_[0][0];
    startTime_ = tempValues_[0][0];
    endTime_   = tempValues_[0][nrLines_ - 1];
}

CorrelationDataSet::~CorrelationDataSet()
{
    // Allocated in read_xvg, destroyed here.
    for (int i = 0; i < nrColumns_; i++)
    {
        sfree(tempValues_[i]);
        tempValues_[i] = nullptr;
    }
    sfree(tempValues_);
    tempValues_ = nullptr;
}

real CorrelationDataSet::getValue(int set, int time) const
{
    if (set + 1 < nrColumns_)
    {
        return tempValues_[set + 1][time];
    }
    else
    {
        return 0;
    }
}
