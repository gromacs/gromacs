/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2015- The GROMACS Authors
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
 * Declares data structure for RAMD output
 *
 * \author Bernd Doser <bernd.doser@h-its.org>
 * \ingroup module_applied_forces
 */

#include "gmxpre.h"

#include "ramdoutputprovider.h"

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/xvgr.h"

struct gmx_output_env_t;

namespace gmx
{

void RAMDOutputProvider::initOutput(FILE* /*fplog*/,
                                    int                     nfile,
                                    const t_filenm          fnm[],
                                    bool                    bAppendFiles,
                                    const gmx_output_env_t* oenv)
{
    if (opt2bSet("-ramd", nfile, fnm))
    {
        std::string filename = opt2fn("-ramd", nfile, fnm);
        if (bAppendFiles)
        {
            fpRAMD_ = gmx_fio_fopen(filename.c_str(), "a+");
        }
        else
        {
            fpRAMD_ = xvgropen(filename.c_str(),
                               "RAMD Ligand-receptor COM distance",
                               "Time (ps)",
                               "Distance (nm)",
                               oenv);
            // std::vector<std::string> setnames;
            // for (int g = 1; g <= params.ngroup; ++g)
            // {
            //     setnames.push_back(std::to_string(g));
            // }
            // xvgrLegend(fpRAMD_, setnames, oenv);
        }
    }
}

void RAMDOutputProvider::addTime(real time)
{
    if (fpRAMD_)
    {
        fprintf(fpRAMD_, "%.4f", time);
    }
}

void RAMDOutputProvider::addDistance(real distance)
{
    if (fpRAMD_)
    {
        fprintf(fpRAMD_, "\t%10.6f", distance);
    }
}

void RAMDOutputProvider::newLine()
{
    if (fpRAMD_)
    {
        fprintf(fpRAMD_, "\n");
    }
}

void RAMDOutputProvider::flush()
{
    if (fpRAMD_)
    {
        fflush(fpRAMD_);
    }
}

void RAMDOutputProvider::finishOutput() {}

} // namespace gmx
