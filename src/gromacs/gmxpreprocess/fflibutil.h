/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2014,2015, by the GROMACS development team, led by
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

#ifndef GMX_GMXPREPROCESS_FFLIBUTIL_H
#define GMX_GMXPREPROCESS_FFLIBUTIL_H

#include <stdio.h>

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
#include <vector>

#include "gromacs/utility/datafilefinder.h"

/*! \brief
 * Enumerates forcefields in the data directories.
 */
std::vector<gmx::DataFileInfo> fflib_enumerate_forcefields();

extern "C" {
#endif

const char *fflib_forcefield_dir_ext();
/* Returns the name of the force field directory extension */

const char *fflib_forcefield_itp();
/* Returns the name of the main forcefield itp file */

const char *fflib_forcefield_doc();
/* Returns the name of the forcefield documentation file */

void fflib_filename_base(const char *filename, char *filebase, int maxlen);
/* Return the base file name of filename in base,
 * i.e. remove path and extension, if present.
 * base should be at least of size maxlen.
 */

int fflib_search_file_end(const char *ffdir,
                          const char *file_end,
                          gmx_bool    bFatalError,
                          char     ***filenames);
/* Search for files ending on file_end in the force field directory fflib.
 * fflib should be in the GROMACS lib.path.
 * Return the number of files and the file names in filenames.
 */

gmx_bool fflib_fexist(const char *file);
/* Check if a file exists in the force field library */

FILE *fflib_open(const char *file);
/* Open force field library file "file" for reading.
 * "file" should contain the whole path to the force field library,
 * either absolute or relative to the current dir.
 */

#ifdef __cplusplus
}
#endif

#endif
