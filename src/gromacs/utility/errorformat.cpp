/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2010- The GROMACS Authors
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

#include <string>

#include "gromacs/utility/basenetwork.h"
#include "gromacs/utility/baseversion.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

/*! \cond internal */
namespace internal
{

void printFatalErrorHeader(FILE* fp, const char* title, const char* func, const char* file, int line)
{
    // In case ProgramInfo is not initialized and there is an issue with the
    // initialization, fall back to "GROMACS".
    const char* programName = "GROMACS";
    try
    {
        programName = getProgramContext().displayName();
    }
    catch (const std::exception&)
    {
    }

    std::fprintf(fp, "\n-------------------------------------------------------\n");
    std::fprintf(fp, "Program:     %s, version %s\n", programName, gmx_version());
    if (file)
    {
        std::fprintf(fp, "Source file: %s (line %d)\n", stripSourcePrefix(file).c_str(), line);
    }
    if (func != nullptr)
    {
        std::fprintf(fp, "Function:    %s\n", func);
    }
    if (gmx_node_num() > 1)
    {
        std::fprintf(fp, "MPI rank:    %d (out of %d)\n", gmx_node_rank(), gmx_node_num());
    }
    std::fprintf(fp, "\n");
    std::fprintf(fp, "%s:\n", title);
}

void printFatalErrorMessageLine(FILE* fp, const char* text, int indent)
{
    gmx::TextLineWrapper wrapper;
    wrapper.settings().setLineLength(78 - indent);
    size_t lineStart = 0;
    size_t length    = std::strlen(text);
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

void printFatalErrorFooter(FILE* fp)
{
    std::fprintf(fp, "\n");
    std::fprintf(fp,
                 "For more information and tips for troubleshooting, please check the GROMACS\n"
                 "website at https://manual.gromacs.org/current/user-guide/run-time-errors.html");
    std::fprintf(fp, "\n-------------------------------------------------------\n");
}

} // namespace internal
//! \endcond

} // namespace gmx
