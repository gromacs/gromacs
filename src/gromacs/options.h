/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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
/*! \defgroup module_options Extensible Handling of Options
 * \ingroup group_utilitymodules
 * \brief
 * Provides functionality for handling options.
 *
 * <H3>Basic Use</H3>
 *
 * Basic interface for providing options is implemented by the Options class
 * and classes defined in basicoptions.h for specifying individual options.
 * Only these are needed if a class wants to provide a set of standard options.
 * When creating an Options object and adding options, it is possible to add
 * descriptions for individual options as well as for the whole set of options.
 * These can then be used to write out help text.
 *
 * Module \ref module_commandline implements classes that assign option values
 * from command line and produce help for programs that use the command line
 * parser.
 *
 * \if libapi
 * <H3>Advanced Use (in library API)</H3>
 *
 * It is possible to extend the module with new option types and/or parsers for
 * option values.
 *
 * To implement new option types, it is necessary to subclass the templates
 * OptionTemplate and OptionStorageTemplate with the type of the values that
 * the option should provide as the template argument.  After this is done, it
 * is possible to add options of this new type using Options::addOption().
 *
 * To implement new parsers, one can use OptionsAssigner, which provides an
 * interface to set values in an Options object.
 *
 * There is also an interface to iterate over all options in an Options object.
 * One should implement the OptionsVisitor interface, and then use
 * OptionsIterator to apply this visitor to the Options object.
 * \endif
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for handling of options.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_H
#define GMX_OPTIONS_H

#include "options/basicoptions.h"
#include "options/filenameoption.h"
#include "options/options.h"

#endif
