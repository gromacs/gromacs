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
/*! \internal \file
 *
 * \brief Implements the MD log file handling routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrun
 */
#include "gmxpre.h"

#include "logging.h"

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"

//! Implements aspects of logfile handling common to opening either for writing or appending.
static void gmx_log_setup(gmx::BinaryInformationSettings settings,
                          const int                      rankIndex,
                          const int                      numRanks,
                          FILE                         * fplog)
{
    int    pid;
    char   host[256];

    GMX_RELEASE_ASSERT(fplog != nullptr, "Log file must be already open");
    gmx_fatal_set_log_file(fplog);

    /* Get some machine parameters */
    gmx_gethostname(host, 256);
    pid = gmx_getpid();

    fprintf(fplog,
            "Log file opened on %s"
            "Host: %s  pid: %d  rank ID: %d  number of ranks:  %d\n",
            gmx_format_current_time().c_str(), host, pid, rankIndex, numRanks);
    try
    {
        settings.extendedInfo(true);
        gmx::printBinaryInformation(fplog, gmx::getProgramContext(), settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    fprintf(fplog, "\n");

    fflush(fplog);
}

void gmx_log_open(const char *lognm,
                  const int   rankIndex,
                  const int   numRanks,
                  FILE     ** fplog)
{
    *fplog = gmx_fio_fopen(lognm, "w+");
    if (*fplog == nullptr)
    {
        GMX_THROW(gmx::FileIOError("Could not open log file" + std::string(lognm)));
    }
    gmx::BinaryInformationSettings settings;
    settings.copyright(true);
    gmx_log_setup(settings, rankIndex, numRanks, *fplog);
}

void gmx_log_append(const int rankIndex,
                    const int numRanks,
                    FILE     *fplog)
{
    GMX_RELEASE_ASSERT(fplog != nullptr, "Log file must be already open");
    fprintf(fplog,
            "\n"
            "\n"
            "-----------------------------------------------------------\n"
            "Restarting from checkpoint, appending to previous log file.\n"
            "\n"
            );
    gmx::BinaryInformationSettings settings;
    settings.copyright(false);
    gmx_log_setup(settings, rankIndex, numRanks, fplog);
}

void gmx_log_close(FILE *fp)
{
    if (fp)
    {
        gmx_fatal_set_log_file(nullptr);
        gmx_fio_fclose(fp);
    }
}
