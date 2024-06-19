/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2016- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/utility/logger.h"

#include <cstdarg>

#include <string>

#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

//! Helper method for reading logging targets from an array.
ILogTarget* getTarget(ILogTarget* targets[MDLogger::LogLevelCount], MDLogger::LogLevel level)
{
    return targets[static_cast<int>(level)];
}

} // namespace

ILogTarget::~ILogTarget() {}


LogEntryWriter& LogEntryWriter::appendTextFormatted(gmx_fmtstr const char* fmt, ...)
{
    std::va_list ap;

    va_start(ap, fmt);
    entry_.text.append(formatStringV(fmt, ap));
    va_end(ap);
    return *this;
}

MDLogger::MDLogger() :
    warning(nullptr), error(nullptr), debug(nullptr), verboseDebug(nullptr), info(nullptr)
{
}

MDLogger::MDLogger(ILogTarget* targets[LogLevelCount]) :
    warning(getTarget(targets, LogLevel::Warning)),
    error(getTarget(targets, LogLevel::Error)),
    debug(getTarget(targets, LogLevel::Debug)),
    verboseDebug(getTarget(targets, LogLevel::VerboseDebug)),
    info(getTarget(targets, LogLevel::Info))
{
}

} // namespace gmx
