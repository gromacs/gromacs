/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
/*! \file
 * \brief
 * Declares functions for initializing the \Gromacs library.
 *
 * Currently, only MPI initialization/finalization management is
 * required, and only if external MPI support is enabled.
 *
 * If MPI is already initialized, we should not call MPI_Init() or
 * MPI_Finalize(). This management object permits \Gromacs test code to
 * nest calls to functions that might normally implement a stand-alone
 * MPI-using tool. It also permits \Gromacs code to be called from code
 * that has already initialized MPI and needs that environment to work
 * and persist after \Gromacs code returns (e.g. \Gromacs tests,
 * external libraries that call \Gromacs code).
 *
 * It does so by maintaining a counter of the number of MPI
 * initializations, and only calling MPI_Init() or MPI_Finalize when
 * it is safe (ie. when the counter is at zero).
 *
 * Thread-MPI initialization and finalization for mdrun is all managed
 * in runner.c.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_INIT_H
#define GMX_UTILITY_INIT_H

// Forward declaration of class ProgramInfo is not sufficient for MSVC if
// the return value of init() is ignored(!)
#include "programinfo.h"

namespace gmx
{

/*! \brief
 * Initializes the \Gromacs library with explicit binary name.
 *
 * \param[in] realBinaryName  Name of the binary
 *     (without Gromacs binary suffix or .exe on Windows).
 * \param[in] argc  argc value passed to main().
 * \param[in] argv  argv array passed to main().
 * \returns   Reference to initialized program information object.
 *
 * This overload is provided for cases where the program may be invoked
 * through a symlink, and it is necessary to know the real name of the
 * binary.
 *
 * Currently, this is tailored for use in command-line/standalone applications.
 * Some additional thought may be required to make it generally usable.
 *
 * The command line arguments are communicated so that they can be
 * parsed on each processor.
 * Arguments are the number of command line arguments, and a pointer to the
 * array of argument strings. Both are allowed to be NULL.
 *
 * Does not throw. Terminates the program on out-of-memory error.
 *
 * \ingroup module_utility
 */
ProgramInfo &init(const char *realBinaryName, int *argc, char ***argv);
/*! \brief
 * Initializes the \Gromacs library.
 *
 * \param[in] argc  argc value passed to main().
 * \param[in] argv  argv array passed to main().
 * \returns   Reference to initialized program information object.
 *
 * Does not throw. Terminates the program on out-of-memory error.
 *
 * \ingroup module_utility
 */
ProgramInfo &init(int *argc, char ***argv);
/*! \brief
 * Deinitializes the \Gromacs library.
 *
 * Decrements the initialization counter, and calls MPI_Finalize()
 * if \Gromacs is compiled with MPI support and the counter has
 * reached zero.  In that case, it is not possible to reinitialize
 * \Gromacs after calling this function.  Instead, call gmx::init() at
 * a higher level, and note that calls to init can be nested safely.
 *
 * \ingroup module_utility
 */
void finalize();

} // namespace gmx

#endif
