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
/*! \libinternal \file
 *
 * \brief Declares the MD log file handling routines.
 *
 * \author Mark Abraham <mark.j.abraham@gmail.com>
 * \inlibraryapi
 * \ingroup module_mdrunutility
 */
#ifndef GMX_MDRUNUTILITY_LOGGING_H
#define GMX_MDRUNUTILITY_LOGGING_H

#include <cstdio>

#include <memory>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/unique_cptr.h"

struct t_fileio;

namespace gmx
{

/*! \brief Close the log file */
void closeLogFile(t_fileio* logfio);

//! Simple guard pointer See unique_cptr for details.
using LogFilePtr = std::unique_ptr<t_fileio, functor_wrapper<t_fileio, closeLogFile>>;

/*! \brief Open the log file for writing/appending.
 *
 * \throws FileIOError when the log file cannot be opened. */
LogFilePtr openLogFile(const char* lognm, bool appendFiles);

/*! \brief Prepare to use the open log file when appending.
 *
 * Does not throw.
 */
void prepareLogAppending(FILE* fplog);

} // namespace gmx

#endif
