/*
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 *
 * For more info, check our website at http://www.gromacs.org
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
 *
 * Values for options can be set using option parsers.
 * Currently, only one is defined: CommandLineParser.
 * As the name suggests, this parser gets values from command-line arguments.
 * Multiple parsers can be used in sequence to provide option values from
 * multiple sources; in such cases, if an option is provided in multiple
 * sources, new values overwrite any values from previous sources.
 *
 * In most cases, it is desirable to be able to provide a help describing the
 * available options.  When creating an Options object and adding options, it
 * is possible to add descriptions for individual options as well as for the
 * whole set of options.  These can then be used to write out a help using one
 * of the provided help writers.
 * Currently, only one is defined: CommandLineHelpWriter
 * (implementation is not yet complete).
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
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 */
/*! \file
 * \brief
 * Public API convenience header for handling of options.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \inpublicapi
 * \ingroup module_options
 */
#ifndef GMX_OPTIONS_H
#define GMX_OPTIONS_H

#include "options/basicoptioninfo.h"
#include "options/basicoptions.h"
#include "options/filenameoptioninfo.h"
#include "options/filenameoption.h"
#include "options/options.h"

#endif
