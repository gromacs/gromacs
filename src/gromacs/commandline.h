/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
/*! \defgroup module_commandline Command Line Program Management (commandline)
 * \ingroup group_utilitymodules
 * \brief
 * Provides functionality for managing command line programs.
 *
 * This module provides utility classes and functions for implementing command
 * line programs.  They are mainly used within \Gromacs, but can also be used
 * from external programs if they want to get a similar user experience to
 * \Gromacs tools.
 *
 * The classes exposed from this module can be roughly divided into two groups:
 *
 *  - Helper classes/functions for implementing the %main() function.
 *    See \ref page_usinglibrary for an overview of those available for user
 *    programs.  These are declared in cmdlineinit.h
 *    (gmx::CommandLineModuleInterface is declared in cmdlinemodule.h and
 *    gmx::CommandLineOptionsInterface in cmdlineoptionsmodule.h).
 *    \if libapi
 *
 *    Additionally, for internal \Gromacs use, gmx::CommandLineModuleManager
 *    provides the functionality to implement the `gmx` wrapper binary, as well
 *    as command line options common to all \Gromacs programs (such as
 *    `-version`).
 *    This is described in more detail at \ref page_wrapperbinary.
 *    \endif
 *
 *  - Helper classes for particular command line tasks:
 *     - gmx::CommandLineParser implements command line parsing to assign
 *       values to gmx::Options (see \ref module_options).
 *     - gmx::CommandLineHelpWriter writes help text for a program that uses
 *       the parser.
 *     - parse_common_args() is an old interface to \Gromacs command line
 *       parsing.  This is still used by many parts of \Gromacs.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 */
/*! \file
 * \brief
 * Public API convenience header for managing command line programs.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_H
#define GMX_COMMANDLINE_H

#include "gromacs/commandline/cmdlinehelpwriter.h"
#include "gromacs/commandline/cmdlineinit.h"
#include "gromacs/commandline/cmdlinemodule.h"
#include "gromacs/commandline/cmdlineoptionsmodule.h"
#include "gromacs/commandline/cmdlineparser.h"
#include "gromacs/commandline/pargs.h"

#endif
