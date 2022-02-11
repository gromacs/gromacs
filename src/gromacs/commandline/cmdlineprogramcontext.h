/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
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
/*! \file
 * \brief
 * Declares gmx::CommandLineProgramContext.
 *
 * This header is installed to support cmdlineinit.h because some compilers
 * don't allow returning a reference to an incomplete type from a function.
 * It should not be necessary to use gmx::CommandLineProgramContext outside the
 * \Gromacs library.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \inlibraryapi
 * \ingroup module_commandline
 */
#ifndef GMX_COMMANDLINE_CMDLINEPROGRAMCONTEXT_H
#define GMX_COMMANDLINE_CMDLINEPROGRAMCONTEXT_H

#include <memory>
#include <string>
#include <vector>

#include "gromacs/utility/programcontext.h"

namespace gmx
{

//! \addtogroup module_commandline
//! \{

/*! \libinternal \brief
 * Allows customization of the way various directories are found by
 * CommandLineProgramContext.
 *
 * For the CommandLineProgramContext constructors that do not take this
 * interface as a parameter, a default implementation is used that forwards
 * the calls to the corresponding methods in gmx::Path.
 *
 * \inlibraryapi
 */
class IExecutableEnvironment
{
public:
    virtual ~IExecutableEnvironment() {}

    /*! \brief
     * Returns the working directory when the program was launched.
     */
    virtual std::string getWorkingDirectory() const = 0;
    /*! \brief
     * Returns list of paths where executables are searched for.
     *
     * The returned list should be in priority order.  An empty string in
     * the returned list corresponds to getWorkindDirectory().
     */
    virtual std::vector<std::string> getExecutablePaths() const = 0;
};

//! Shorthand for a smart pointer to IExecutableEnvironment.
typedef std::unique_ptr<IExecutableEnvironment> ExecutableEnvironmentPointer;

/*! \libinternal \brief
 * Program context implementation for command line programs.
 *
 * Constructors are provided mostly for unit testing purposes; in normal usage,
 * a single CommandLineProgramContext object is constructed with
 * initForCommandLine() in the beginning of the program.
 * The returned object can be explicitly passed to other methods, or accessed
 * through getProgramContext().
 *
 * Unless explicitly noted otherwise, methods in this class may throw
 * std::bad_alloc on out-of-memory conditions, but do not throw other
 * exceptions.
 *
 * \inlibraryapi
 */
class CommandLineProgramContext : public IProgramContext
{
public:
    /*! \brief
     * Constructs an empty context object.
     *
     * All methods in the constructed object return dummy values.
     */
    CommandLineProgramContext();
    /*! \brief
     * Initializes a program context object with binary name only.
     *
     * \param[in] binaryName  Name of the binary.
     *
     * This is needed for unit testing purposes.
     * The constructed object works as if the command line consisted of
     * only of the binary name.
     */
    explicit CommandLineProgramContext(const char* binaryName);
    /*! \brief
     * Initializes a program context object based on command line.
     *
     * \param[in] argc  argc value passed to main().
     * \param[in] argv  argv array passed to main().
     */
    CommandLineProgramContext(int argc, const char* const argv[]);
    /*! \brief
     * Initializes a program context object based on command line.
     *
     * \param[in] argc  argc value passed to main().
     * \param[in] argv  argv array passed to main().
     * \param[in] env   Customizes the way the binary name is handled.
     *
     * This overload allows one to customize the way the binary is located
     * by providing a custom IExecutableEnvironment implementation.
     * This is mainly useful for testing purposes to make it possible to
     * test different paths without setting environment variables, changing
     * the working directory or doing other process-wide operations.
     * It may also be useful for making Gromacs behave better when linked
     * into a non-Gromacs executable (with possible extensions in
     * IExecutableEnvironment).
     */
    CommandLineProgramContext(int argc, const char* const argv[], ExecutableEnvironmentPointer env);
    ~CommandLineProgramContext() override;

    /*! \brief
     * Sets a display name for the binary.
     *
     * \throws std::bad_alloc if out of memory.
     *
     * This is used with the wrapper binary to add the name of the invoked
     * module to the name of the binary shown.
     *
     * It is not threadsafe if there are concurrent calls to displayName()
     * before this method has returned.  Thread safety is not required for
     * the normal initialization sequence of command line programs; it is
     * called in the wrapper binary before the control passes to the actual
     * module which may create threads.
     */
    void setDisplayName(const std::string& name);

    /*! \brief
     * Returns the name of the binary as it was invoked without any path.
     *
     * Does not throw.
     */
    const char* programName() const override;
    /*! \brief
     * Returns a display name of the current module.
     *
     * The returned value equals programName(), unless a separate display
     * name has been set with setDisplayName().
     *
     * Does not throw.
     */
    const char* displayName() const override;
    /*! \brief
     * Returns the full path of the running binary.
     *
     * \throws std::bad_alloc if out of memory.
     *
     * Returns argv[0] if there was an error in finding the absolute path.
     */
    const char* fullBinaryPath() const override;
    /*! \brief
     * Returns the installation prefix (for finding \Gromacs data files).
     *
     * \throws std::bad_alloc if out of memory.
     *
     * Returns a hardcoded path set during configuration time if there is
     * an error in finding the library data files.
     */
    InstallationPrefixInfo installationPrefix() const override;
    /*! \brief
     * Returns the full command line used to invoke the binary.
     *
     * Does not throw.
     */
    const char* commandLine() const override;

private:
    class Impl;

    std::unique_ptr<Impl> impl_;
};

} // namespace gmx

#endif
