/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2018, by the GROMACS development team, led by
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

#include "gmxpre.h"

#include "densfitoutput.h"

#include "gromacs/fileio/xvgr.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/pleasecite.h"

namespace gmx
{

DensfitOutput::DensfitOutput(FILE * fplog, const char *fnLogFile, const char * fnDensity, const gmx_output_env_t *oenv, bool bAppend)
{
    if (bAppend)
    {
        logFile_ = gmx_fio_fopen(fnLogFile, "a");
    }
    else
    {
        please_cite(fplog, "Tama2008");
        logFile_ = xvgropen(fnLogFile, "Density fitting correlation coefficient", "Time (ps)",
                            "correlation coefficient", oenv);
        fflush(logFile_);
    }
    fnDensity_ = fnDensity;
}

FILE* DensfitOutput::logFile()
{
    return logFile_;
}

DensfitOutput::~DensfitOutput()
{
    if (logFile_) // Close the density fitting output file
    {
        if (bAppend_)
        {
            gmx_fio_fclose(logFile_);
        }
        else
        {
            xvgrclose(logFile_);
        }
    }
}

std::string DensfitOutput::outputMapFileNameWithStep(gmx_int64_t step) const
{
    return fnDensity_.substr(0, fnDensity_.find_last_of("."))+std::to_string(step)+".ccp4";
}
}
