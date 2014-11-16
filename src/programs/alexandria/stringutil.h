/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014, by the GROMACS development team, led by
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
#ifndef STRINGUTIL_H
#define STRINGUTIL_H

#include <string>
#include <vector>

/*! \brief
 * Return the index of the selected term
 * \param[in] opts List of strings to select from
 * \returns the index or -1 if not found
 */
extern int get_option(const char **opts);

/*! \brief
 * Split a string into substrings separated by a delimiter character
 * \param[in] s The string to be split
 * \parin[in] delim The delimiting character
 * \param[in] elems A vector of substring elements
 * \returns The vector of substring elements (same as input elems)
 */
extern std::vector<std::string> &split(const std::string        &s,
                                       char                      delim,
                                       std::vector<std::string> &elems);

/*! \brief
 * Split a string into substrings separated by a delimiter character
 * \param[in] s The string to be split
 * \parin[in] delim The delimiting character
 * \returns A vector of substring elements
 */
extern std::vector<std::string> split(const std::string &s,
                                      char               delim);

/** return new string with f printed. */
extern std::string gmx_ftoa(double f);

/** return new string with f printed. */
extern std::string gmx_itoa(int f);

#endif
