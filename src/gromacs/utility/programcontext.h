/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013, by the GROMACS development team, led by
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
 * Declares gmx::ProgramContextInterface and related methods.
 *
 * This header is installed to support init.h because some compilers don't
 * allow returning a reference to an incomplete type from a function.
 * It should not be necessary to use gmx::ProgramInfo outside the Gromacs
 * library.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_PROGRAMCONTEXT_H
#define GMX_UTILITY_PROGRAMCONTEXT_H

namespace gmx
{

class ProgramContextInterface
{
    public:
        virtual const char *programName() const = 0;
        virtual const char *displayName() const = 0;
        virtual const char *fullBinaryPath() const = 0;
        virtual const char *commandLine() const = 0;

    protected:
        virtual ~ProgramContextInterface() {}
};

/*! \brief
 * Returns the singleton ProgramInfo object.
 *
 * \returns The same object as initialized with the last call to init().
 * \throws  std::bad_alloc if out of memory (only if this is the first
 *      call and init() has not been called either).
 * \throws  tMPI::system_error on thread synchronization errors.
 */
const ProgramContextInterface &getProgramContext();
void setProgramContext(const ProgramContextInterface *context);

} // namespace gmx

#endif
