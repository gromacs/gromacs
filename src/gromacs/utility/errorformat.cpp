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
 * Implements functions declared in errorformat.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "errorformat.h"

#include <cctype>
#include <cstdio>
#include <cstring>

#include "buildinfo.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! \cond internal */
namespace internal
{

void printFatalErrorHeader(FILE *fp, const char *title,
                           const char *func, const char *file, int line)
{
    // In case ProgramInfo is not initialized and there is an issue with the
    // initialization, fall back to "GROMACS".
    const char *programName = "GROMACS";
    try
    {
        programName = getProgramContext().displayName();
    }
    catch (const std::exception &)
    {
    }

    std::fprintf(fp, "\n-------------------------------------------------------\n");
    std::fprintf(fp, "Program:     %s, %s\n", programName, gmx_version());
    if (file != NULL)
    {
        // TODO: Check whether this works on Windows. If it doesn't, perhaps
        // add Path::startsWith().
        if (startsWith(file, CMAKE_SOURCE_DIR))
        {
            file += std::strlen(CMAKE_SOURCE_DIR);
            if (file[0] == '/' || file[0] == '\\')
            {
                ++file;
            }
        }
        std::fprintf(fp, "Source file: %s (line %d)\n", file, line);
    }
    if (func != NULL)
    {
        std::fprintf(fp, "Function:    %s\n", func);
    }
    std::fprintf(fp, "\n");
    std::fprintf(fp, "%s:\n", title);
}

void printFatalErrorMessageLine(FILE *fp, const char *text, int indent)
{
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setLineLength(78 - indent);
    size_t               lineStart = 0;
    size_t               length    = std::strlen(text);
    while (lineStart < length)
    {
        size_t nextLineStart = wrapper.findNextLine(text, lineStart);
        int    lineLength    = static_cast<int>(nextLineStart - lineStart);
        while (lineLength > 0 && std::isspace(text[lineStart + lineLength - 1]))
        {
            --lineLength;
        }
        std::fprintf(fp, "%*s%.*s\n", indent, "", lineLength, text + lineStart);
        lineStart = nextLineStart;
    }
}

void printFatalErrorFooter(FILE *fp)
{
    std::fprintf(fp, "\n");
    std::fprintf(fp, "For more information and tips for troubleshooting, please check the GROMACS\n"
                 "website at http://www.gromacs.org/Documentation/Errors");
    std::fprintf(fp, "\n-------------------------------------------------------\n");
}

}   // namespace internal
//! \endcond

} // namespace gmx
