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
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/sysinfo.h"

void gmx_log_open(const char *lognm, const t_commrec *cr,
                  gmx_bool bAppendFiles, FILE** fplog)
{
    int    pid;
    char   host[256];
    char   timebuf[STRLEN];
    FILE  *fp = *fplog;

    if (!bAppendFiles)
    {
        fp = gmx_fio_fopen(lognm, bAppendFiles ? "a+" : "w+" );
    }

    gmx_fatal_set_log_file(fp);

    /* Get some machine parameters */
    gmx_gethostname(host, 256);
    pid = gmx_getpid();
    gmx_format_current_time(timebuf, STRLEN);

    if (bAppendFiles)
    {
        fprintf(fp,
                "\n"
                "\n"
                "-----------------------------------------------------------\n"
                "Restarting from checkpoint, appending to previous log file.\n"
                "\n"
                );
    }

    fprintf(fp,
            "Log file opened on %s"
            "Host: %s  pid: %d  rank ID: %d  number of ranks:  %d\n",
            timebuf, host, pid, cr->nodeid, cr->nnodes);
    try
    {
        gmx::BinaryInformationSettings settings;
        settings.extendedInfo(true);
        settings.copyright(!bAppendFiles);
        gmx::printBinaryInformation(fp, gmx::getProgramContext(), settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
    fprintf(fp, "\n");

    fflush(fp);

    *fplog = fp;
}

void gmx_log_close(FILE *fp)
{
    if (fp)
    {
        gmx_fatal_set_log_file(nullptr);
        gmx_fio_fclose(fp);
    }
}
