/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2010,2011,2012,2013,2014, by the GROMACS development team, led by
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
/*! \defgroup module_options Extensible Handling of Options (options)
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
 * The sequence charts below provides an overview of how the options work from
 * usage perspective.  They include two fictional modules, A and B, that provide
 * options, and a main routine that manages these.  The first chart shows a
 * typical initialization sequence, where the main routine creates an options
 * object, and calls an initOptions() method in each module that can provide
 * options (the modules may also request their submodules to add their own
 * options).  Each module uses gmx::Options::addOption() to add the options
 * they require, and specify output variables into which the options values are
 * stored.
 * \msc
 *     main,
 *     options [ label="Options", URL="\ref gmx::Options" ],
 *     A [ label="module A" ],
 *     B [ label="module B" ];
 *
 *     main box B [ label="main owns all objects" ];
 *     main => options [ label="create", URL="\ref gmx::Options::Options()" ];
 *     main => A [ label="initOptions()" ];
 *     A => options [ label="addOption()", URL="\ref gmx::Options::addOption()" ];
 *     ...;
 *     main << A;
 *     main => B [ label="initOptions()" ];
 *     B => options [ label="addOption()", URL="\ref gmx::Options::addOption()" ];
 *     ...;
 *     main << B;
 * \endmsc
 *
 * After all options have been specified, they can be parsed.  A parser splits
 * the input into option-value pairs (one option may have multiple values), and
 * passes these into the gmx::Options object, which is responsible for
 * converting them into the appropriate types and storing the values into the
 * variables provided in the calls to gmx::Options::addOption().
 * \msc
 *     main,
 *     parser [ label="parser" ],
 *     options [ label="Options", URL="\ref gmx::Options" ],
 *     A [ label="module A" ],
 *     B [ label="module B" ];
 *
 *     main => parser [ label="parse()" ];
 *     parser => options [ label="assign(string)" ];
 *     options -> A [ label="set variable" ];
 *     parser => options [ label="assign(string)" ];
 *     options -> B [ label="set variable" ];
 *     ...;
 * \endmsc
 *
 * After all options have been parsed (possibly using multiple different
 * parsers), gmx::Options::finish() is called.  This performs final
 * validation of the options and may further adjust the values stored in the
 * output variables (see documentation on individual option types on when this
 * may happen).
 * \msc
 *     main,
 *     options [ label="Options", URL="\ref gmx::Options" ],
 *     A [ label="module A" ],
 *     B [ label="module B" ];
 *
 *     main => options [ label="finish()", URL="\ref gmx::Options::finish()" ];
 *     options -> A [ label="set variable" ];
 *     options -> B [ label="set variable" ];
 *     ...;
 * \endmsc
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

#include "gromacs/fileio/filenm.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/options/options.h"

#endif
