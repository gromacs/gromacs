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

#include <cmath>
#include <sstream>
#include "gromacs/fileio/xvgr.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/utility/smalloc.h"
#include "correlationDataSet.h"
#include "testutils/testfilemanager.h"





CorrelationDataSet::CorrelationDataSet(std::string fileName)
{
    double ** tempValues;
    fileName    = gmx::test::TestFileManager::getInputFilePath(fileName.c_str());
    nrLines     = read_xvg(fileName.c_str(), &tempValues, &nrColums );

    values.reserve(nrLines);

    //Comverting double to real
    for (int i = 0; i  < nrLines; i++)
    {
        values[i]   = (real)tempValues[1][i];

    }
    dt        =  tempValues[0][1] - tempValues[0][0];
    startTime = tempValues[0][0];
    endTime   = tempValues[0][nrLines-1];



    //alocated in read_xvg
    for (int i = 0; i < nrColums; i++)
    {
        sfree(tempValues[i]);
        tempValues[i] = NULL;
    }
    sfree(tempValues);
    tempValues = NULL;
}


real CorrelationDataSet::getValue(int x)
{
    return values[x];
}

int CorrelationDataSet::getNrColums()
{
    return nrColums;
}

int CorrelationDataSet::getNrLines()
{
    return nrLines;
}

real CorrelationDataSet::getStartTime()
{
    return startTime;
}

real CorrelationDataSet::getEndTime()
{
    return endTime;
}

real CorrelationDataSet::getDt()
{
    return dt;
}
