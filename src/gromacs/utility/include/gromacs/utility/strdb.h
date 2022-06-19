/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
 * \brief
 * Declares C functions for reading files with a list of strings.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_STRDB_H
#define GMX_UTILITY_STRDB_H

#include <cstdio>

#include "gromacs/utility/basedefinitions.h"

/*! \brief
 * Reads a line of at most n characters from *fp to line.
 *
 * Comment ';...' and leading spaces are removed, empty lines are skipped.
 * Return FALSE when eof.
 */
gmx_bool get_a_line(FILE* fp, char line[], int n);

/*! \brief
 * Read a header between '[' and ']' from line to header.
 *
 * Returns FALSE if no header is found.
 */
gmx_bool get_header(char line[], char header[]);

/*! \brief
 * Opens file db, or if non-existant file $GMXLIB/db and read strings.
 *
 * First line in the file needs to specify the number of strings following.
 * Returns the number of strings.
 */
int get_lines(const char* db, char*** strings);

/*! \brief
 * Searches an array of strings for key, return the index if found.
 *
 * Returns -1 if not found.
 */
int search_str(int nstr, char** str, char* key);

#endif
