/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2016,2017, by the GROMACS development team, led by
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
/*! \libinternal \file
 * \brief
 * Declares functions for using keyvaluetree.h with Options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_TREESUPPORT_H
#define GMX_OPTIONS_TREESUPPORT_H

namespace gmx
{

class IKeyValueTreeErrorHandler;
class KeyValueTreeObject;
class Options;

//! \cond libapi

/*! \libinternal \brief
 * Assigns option values from a given KeyValueTreeObject.
 *
 * Each property with a simple value (or an array of simple values) is assigned
 * to an option with the same name.  Objects and arrays of objects are assigned
 * to sections with the same name.
 *
 * \ingroup module_options
 */
void assignOptionsFromKeyValueTree(Options                   *options,
                                   const KeyValueTreeObject  &tree,
                                   IKeyValueTreeErrorHandler *errorHandler);
/*! \libinternal \brief
 * Checks that a given KeyValueTreeObject can be assigned to given Options.
 *
 * Throws an exception if `tree` contains any values that are not recognized by
 * `options`.  Does not verify the type of the values, only that an option with
 * the correct names exists.
 *
 * \ingroup module_options
 */
void checkForUnknownOptionsInKeyValueTree(const KeyValueTreeObject &tree,
                                          const Options            &options);
/*! \libinternal \brief
 * Adjusts a KeyValueTreeObject to the structure of given Options.
 *
 * Assumes that all values in the input KeyValueTreeObject are valid values for
 * the options.  The output has all the values in the input, but in the order
 * they are in the options.  Values are also converted to the native type for
 * the underlying option (e.g., strings are parsed to integers if the option
 * accepts those).  For any option that does not have a corresponding value in
 * the input, the output has it with a default value (if one exists for the
 * option).
 *
 * Any values in `tree` that do not have matching options are not present in
 * the output.  If this is not desirable, call
 * checkForUnknownOptionsInKeyValueTree() before calling this function to
 * ensure that no such values are present.
 *
 * Does not currently work for option sections in an array.
 *
 * \ingroup module_options
 */
KeyValueTreeObject
adjustKeyValueTreeFromOptions(const KeyValueTreeObject &tree,
                              const Options            &options);

//! \endcond

} // namespace gmx

#endif
