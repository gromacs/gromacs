/*
 * This source file is part of the Alexandria Chemistry Toolkit.
 *
 * Copyright (C) 2014-2020 
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour, 
 *             Paul J. van Maaren, 
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, 
 * Boston, MA  02110-1301, USA.
 */
 
/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
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
int get_option(const char **opts);

/*! \brief
 * Split a string into substrings separated by a delimiter character
 * \param[in] s The string to be split
 * \parin[in] delim The delimiting character
 * \param[in] elems A vector of substring elements
 * \returns The vector of substring elements (same as input elems)
 */
std::vector<std::string> split(const std::string        &s,
                               char                      delim,
                               std::vector<std::string> &elems);

/*! \brief
 * Split a string into substrings separated by a delimiter character
 * \param[in] s The string to be split
 * \parin[in] delim The delimiting character
 * \returns A vector of substring elements
 */
std::vector<std::string> split(const std::string &s,
                               char               delim);

/** return new string with f printed. */
std::string gmx_ftoa(double f);

/** return new string with f printed. */
std::string gmx_dtoa(double f);

/** return new string with f printed. */
std::string gmx_itoa(int f);

/*! \brief
 *
 * Read a double from a string and apply error checking.
 * \param[in] str         The string to read a double from
 * \param[in] description If there is any issue this will be used
 *                        to describe the variable being read.
 *                        May be nullptr.
 * \return The value read or -1 in case of issues.
 */
double my_atof(const char *str, const char *description);

#endif
