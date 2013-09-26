/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
 * Declares gmx::ProgramInfo.
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
#ifndef GMX_UTILITY_PROGRAMINFO_H
#define GMX_UTILITY_PROGRAMINFO_H

#include <string>

#include "common.h"

namespace gmx
{

/*! \libinternal \brief
 * Helper class for managing information about the running binary.
 *
 * This class provides access to the name of the binary currently running, as
 * well as other information derived from it.
 *
 * ProgramInfo::init() should be called before any other (C++) Gromacs calls in
 * a command-line program, as the information is used for printing error
 * messages.
 *
 * Constructors are provided mostly for unit testing purposes; in normal usage,
 * a single ProgramInfo object is constructed with init() in the beginning of
 * the program.  The returned object can be explicitly passed to other methods,
 * or accessed through getInstance().
 *
 * Unless explicitly noted otherwise, methods in this class may throw
 * std::bad_alloc on out-of-memory conditions, but do not throw other
 * exceptions.
 *
 * \inlibraryapi
 * \ingroup module_utility
 */
class ProgramInfo
{
    public:
        /*! \brief
         * Returns the singleton ProgramInfo object.
         *
         * \returns The same object as initialized with the last call to init().
         * \throws  std::bad_alloc if out of memory (only if this is the first
         *      call and init() has not been called either).
         * \throws  tMPI::system_error on thread synchronization errors.
         */
        static const ProgramInfo &getInstance();
        /*! \brief
         * Initializes global program information.
         *
         * \param[in] argc  argc value passed to main().
         * \param[in] argv  argv array passed to main().
         * \returns   Reference to initialized program information object.
         *
         * The return value of realBinaryName() is the same as
         * invariantProgramName().
         *
         * Does not throw. Terminates the program on out-of-memory error.
         */
        static ProgramInfo &init(int argc, const char *const argv[]);
        /*! \brief
         * Initializes global program information with explicit binary name.
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
         * Does not throw. Terminates the program on out-of-memory error.
         */
        static ProgramInfo &init(const char *realBinaryName,
                                 int argc, const char *const argv[]);

        /*! \brief
         * Constructs an empty program info objects.
         *
         * All methods in the constructed object return dummy values.
         */
        ProgramInfo();
        /*! \brief
         * Initializes a program information object with binary name only.
         *
         * \param[in] realBinaryName  Name of the binary
         *     (without Gromacs binary suffix or .exe on Windows).
         *
         * This is needed for unit testing purposes.
         * The constructed object works as if the command line consisted of
         * only of the binary name.
         */
        explicit ProgramInfo(const char *realBinaryName);
        /*! \brief
         * Initializes a program information object based on command line.
         *
         * \param[in] argc  argc value passed to main().
         * \param[in] argv  argv array passed to main().
         */
        ProgramInfo(int argc, const char *const argv[]);
        /*! \brief
         * Initializes a program information object based on binary name and
         * command line.
         *
         * \param[in] realBinaryName  Name of the binary
         *     (without Gromacs binary suffix or .exe on Windows).
         * \param[in] argc  argc value passed to main().
         * \param[in] argv  argv array passed to main().
         */
        ProgramInfo(const char *realBinaryName,
                    int argc, const char *const argv[]);
        ~ProgramInfo();

        /*! \brief
         * Sets a display name for the binary.
         *
         * \throws std::bad_alloc if out of memory.
         * \throws tMPI::system_error on thread synchronization errors.
         *
         * This is used with the wrapper binary to add the name of the invoked
         * module to the name of the binary shown.
         */
        void setDisplayName(const std::string &name);

        /*! \brief
         * Returns the real name of the binary.
         *
         * The returned value is comparable to invariantProgramName(), i.e., it
         * has suffixes and OS-specific extensions removed.
         *
         * Does not throw.
         */
        const std::string &realBinaryName() const;
        /*! \brief
         * Returns the path and name of the binary as it was invoked.
         *
         * Does not throw.
         */
        const std::string &programNameWithPath() const;
        /*! \brief
         * Returns the name of the binary as it was invoked without any path.
         *
         * Does not throw.
         */
        const std::string &programName() const;
        /*! \brief
         * Returns an invariant name of the binary.
         *
         * The returned value has OS-specific extensions (.exe on Windows)
         * removed, as well as any binary suffix that was configured.
         *
         * Does not throw.
         */
        const std::string &invariantProgramName() const;
        /*! \brief
         * Returns a display name of the current module.
         *
         * \throws  tMPI::system_error on thread synchronization errors.
         *
         * The returned value equals programName(), unless a separate display
         * name has been set with setDisplayName().
         */
        const std::string &displayName() const;
        /*! \brief
         * Returns the full command line used to invoke the binary.
         *
         * Does not throw.
         */
        const std::string &commandLine() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
