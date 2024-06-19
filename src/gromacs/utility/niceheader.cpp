/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2017- The GROMACS Authors
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
 * Implements functions from niceheader.h.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "gromacs/utility/niceheader.h"

#include <string>

#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/sysinfo.h"
#include "gromacs/utility/textwriter.h"

namespace gmx
{

void niceHeader(TextWriter* writer, const char* fn, char commentChar)
{
    char userbuf[256];
    char hostbuf[256];

    /* Write a nice header for an output file */
    writer->writeLine(formatString("%c", commentChar));
    writer->writeLine(formatString("%c\tFile '%s' was generated", commentChar, fn ? fn : "unknown"));

    int uid = gmx_getuid();
    gmx_getusername(userbuf, 256);
    gmx_gethostname(hostbuf, 256);

    writer->writeLine(formatString("%c\tBy user: %s (%d)", commentChar, userbuf, uid));
    writer->writeLine(formatString("%c\tOn host: %s", commentChar, hostbuf));
    writer->writeLine(formatString("%c\tAt date: %s", commentChar, gmx_format_current_time().c_str()));
    writer->writeLine(formatString("%c", commentChar));
}

} // namespace gmx
