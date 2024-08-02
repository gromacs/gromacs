/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2018- The GROMACS Authors
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
 *
 * \brief Implements the MD log file handling routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_mdrunutility
 */
#include "gmxpre.h"

#include "logging.h"

#include <filesystem>
#include <memory>
#include <string>

#include "gromacs/fileio/gmxfio.h"
#include "gromacs/utility/binaryinformation.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"

namespace gmx
{

//! Implements aspects of logfile handling common to opening either for writing or appending.
static void prepareLogFile(BinaryInformationSettings settings, FILE* fplog)
{
    GMX_RELEASE_ASSERT(fplog != nullptr, "Log file must be already open");
    // TODO This function is writing initial content to the log
    // file. Preparing the error output handling should happen at some
    // later point, using this log file, but should not be done at the
    // same time as writing content. Move this call there.
    gmx_fatal_set_log_file(fplog);

    try
    {
        settings.extendedInfo(true).processId(true);
        printBinaryInformation(fplog, getProgramContext(), settings);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
    fprintf(fplog, "\n");

    fflush(fplog);
}

LogFilePtr openLogFile(const char* lognm, bool appendFiles)
{
    const char* fileOpeningMode = "w+";
    if (appendFiles)
    {
        fileOpeningMode = GMX_FAHCORE ? "a" : "r+";
    }

    LogFilePtr logfio(gmx_fio_open(lognm, fileOpeningMode));
    if (!logfio)
    {
        GMX_THROW(FileIOError("Could not open log file" + std::string(lognm)));
    }
    // If appending, then there is no need to write this header
    // information, and we don't want to change the file until the
    // checksum has been computed.
    if (!appendFiles)
    {
        FILE*                          fplog = gmx_fio_getfp(logfio.get());
        gmx::BinaryInformationSettings settings;
        settings.copyright(true);
        prepareLogFile(settings, fplog);
    }
    return logfio;
}

void prepareLogAppending(FILE* fplog)
{
    GMX_RELEASE_ASSERT(fplog != nullptr, "Log file must be already open");
    fprintf(fplog,
            "\n"
            "\n"
            "-----------------------------------------------------------\n"
            "Restarting from checkpoint, appending to previous log file.\n"
            "\n");
    gmx::BinaryInformationSettings settings;
    settings.copyright(false);
    prepareLogFile(settings, fplog);
}

void closeLogFile(t_fileio* logfio)
{
    if (logfio)
    {
        gmx_fatal_set_log_file(nullptr);
        gmx_fio_close(logfio);
    }
}

} // namespace gmx
