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
/*! \defgroup module_options Extensible Handling of Options (options)
 * \ingroup group_utilitymodules
 * \brief
 * Provides functionality for handling options.
 *
 * <H3>Basic Use</H3>
 *
 * Code that provides options does so using methods in gmx::IOptionsContainer
 * and classes defined in basicoptions.h.
 * Only these are needed if a class wants to provide a set of standard options
 * (other modules can provide additional option types, such as
 * gmx::SelectionOption).
 * For each option, the caller provides an output variable that will receive
 * the final value of the option once user input has been parsed.
 * When adding options, it is possible to also provide descriptions for the
 * options for use in generated help text.
 *
 * Generic code that handles the user input does so by creating a gmx::Options
 * instance and passing it (as gmx::IOptionsContainer) to the classes that add
 * the actual options.  It can then use a parser to set values to the options.
 * Final values for the options can be inspected in the code that added the
 * individual options, from the provided output variables.
 *
 * The sequence charts below provides an overview of how the options work from
 * usage perspective.  They include two fictional modules, A and B, that provide
 * options, and a main routine that manages these.  The first chart shows a
 * typical initialization sequence, where the main routine creates an options
 * object, and calls an initOptions() method in each module that can provide
 * options (the modules may also request their submodules to add their own
 * options).  Each module uses gmx::IOptionsContainer::addOption() to add the
 * options they require, and specify output variables into which the options
 * values are stored.
 * \msc
 *     main,
 *     options [ label="Options", URL="\ref gmx::Options" ],
 *     A [ label="module A" ],
 *     B [ label="module B" ];
 *
 *     main box B [ label="main owns all objects" ];
 *     main => options [ label="create", URL="\ref gmx::Options::Options()" ];
 *     main => A [ label="initOptions()" ];
 *     A => options [ label="addOption()", URL="\ref gmx::IOptionsContainer::addOption()" ];
 *     ...;
 *     main << A;
 *     main => B [ label="initOptions()" ];
 *     B => options [ label="addOption()", URL="\ref gmx::IOptionsContainer::addOption()" ];
 *     ...;
 *     main << B;
 * \endmsc
 *
 * After all options have been specified, they can be parsed.  A parser splits
 * the input into option-value pairs (one option may have multiple values), and
 * passes these into the gmx::Options object, which is responsible for
 * converting them into the appropriate types and storing the values into the
 * variables provided in the calls to gmx::IOptionsContainer::addOption().
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
 * is possible to add options of this new type using IOptionsContainer::addOption().
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

#include "gromacs/fileio/filetypes.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/filenameoptionmanager.h"
#include "gromacs/options/ioptionsbehavior.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/options/options.h"

#endif
