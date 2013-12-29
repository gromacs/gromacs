/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares functions for initializing the \Gromacs library for command-line use.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEINIT_H
#define GMX_COMMANDLINE_CMDLINEINIT_H

// Forward declaration of class CommandLineProgramContext is not sufficient for
// MSVC if the return value of initForCommandLine() is ignored(!)
#include "cmdlineprogramcontext.h"

namespace gmx
{

/*! \brief
 * Initializes the \Gromacs library for command-line use.
 *
 * \param[in] argc  argc value passed to main().
 * \param[in] argv  argv array passed to main().
 * \returns   Reference to initialized program context object.
 *
 * This function is tailored for use in command line applications.
 * For other usage, combination of gmx::init() and gmx::setProgramContext()
 * provides more flexible initialization alternatives.
 * Unlike gmx::init(), calls to this method cannot be nested.
 *
 * The command line arguments are communicated so that they can be
 * parsed on each processor.
 * \p argc and \p argv are passed to gmx::init(); see there for additional
 * discussion.  This method does not place any additional limitations, but
 * generally there should be no need to pass NULL values.
 *
 * Does not throw. Terminates the program on out-of-memory error.
 *
 * This method is not thread-safe, since it is intended to be the first method
 * called.  See setProgramContext() for additional discussion.
 *
 * \see gmx::init()
 * \see setProgramContext()
 * \ingroup module_commandline
 */
CommandLineProgramContext &initForCommandLine(int *argc, char ***argv);
/*! \brief
 * Deinitializes the \Gromacs library after initForCommandLine().
 *
 * Calls gmx::finalize() and additionally undoes the work done by
 * initForCommandLine().
 *
 * \see gmx::finalize()
 * \ingroup module_commandline
 */
void finalizeForCommandLine();

} // namespace gmx

#endif
