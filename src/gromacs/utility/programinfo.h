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
 *
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
/*! \libinternal \file
 * \brief
 * Declares gmx::ProgramInfo.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
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
        static const ProgramInfo &init(int argc, const char *const argv[]);
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
        static const ProgramInfo &init(const char *realBinaryName,
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
         * Returns the real name of the binary.
         *
         * The returned value is comparable to invariantProgramName(), i.e., it
         * has suffixes and OS-specific extensions removed.
         */
        std::string realBinaryName() const;
        /*! \brief
         * Returns the path and name of the binary as it was invoked.
         */
        std::string programNameWithPath() const;
        /*! \brief
         * Returns the name of the binary as it was invoked without any path.
         */
        std::string programName() const;
        /*! \brief
         * Returns an invariant name of the binary.
         *
         * The returned value has OS-specific extensions (.exe on Windows)
         * removed, as well as any binary suffix that was configured.
         */
        std::string invariantProgramName() const;

    private:
        class Impl;

        PrivateImplPointer<Impl> impl_;
};

} // namespace gmx

#endif
