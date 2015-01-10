/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013,2014,2015, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements gmx::CommandLineProgramContext.
 *
 * See \linktodevmanual{relocatable-binaries,developer guide section on
 * relocatable binaries} for explanation of the searching logic.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "gmxpre.h"

#include "cmdlineprogramcontext.h"

#include "config.h"

#include <cstdlib>
#include <cstring>

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "thread_mpi/mutex.h"

#include "buildinfo.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/file.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{

//! \addtogroup module_commandline
//! \{

/*! \brief
 * Quotes a string if it contains spaces.
 */
std::string quoteIfNecessary(const char *str)
{
    const bool bSpaces = (std::strchr(str, ' ') != NULL);
    if (bSpaces)
    {
        return formatString("'%s'", str);
    }
    return str;
}

/*! \brief
 * Default implementation for ExecutableEnvironmentInterface.
 *
 * Used if ExecutableEnvironmentInterface is not explicitly provided when
 * constructing CommandLineProgramContext.
 */
class DefaultExecutableEnvironment : public ExecutableEnvironmentInterface
{
    public:
        //! Allocates a default environment.
        static ExecutableEnvironmentPointer create()
        {
            return ExecutableEnvironmentPointer(new DefaultExecutableEnvironment());
        }

        DefaultExecutableEnvironment()
            : initialWorkingDirectory_(Path::getWorkingDirectory())
        {
        }

        virtual std::string getWorkingDirectory() const
        {
            return initialWorkingDirectory_;
        }
        virtual std::vector<std::string> getExecutablePaths() const
        {
            return Path::getExecutablePaths();
        }

    private:
        std::string   initialWorkingDirectory_;
};

/*! \brief
 * Finds the absolute path of the binary from \c argv[0].
 *
 * \param[in] invokedName \c argv[0] the binary was invoked with.
 * \param[in] env         Executable environment.
 * \returns   The full path of the binary.
 *
 * If a binary with the given name cannot be located, \p invokedName is
 * returned.
 */
std::string findFullBinaryPath(const std::string                    &invokedName,
                               const ExecutableEnvironmentInterface &env)
{
    std::string searchName = invokedName;
    // On Windows & Cygwin we need to add the .exe extension,
    // or we wont be able to detect that the file exists.
#if (defined GMX_NATIVE_WINDOWS || defined GMX_CYGWIN)
    if (!endsWith(searchName, ".exe"))
    {
        searchName.append(".exe");
    }
#endif
    if (!Path::containsDirectory(searchName))
    {
        // No directory in name means it must be in the path - search it!
        std::vector<std::string>                 pathEntries = env.getExecutablePaths();
        std::vector<std::string>::const_iterator i;
        for (i = pathEntries.begin(); i != pathEntries.end(); ++i)
        {
            const std::string &dir      = i->empty() ? env.getWorkingDirectory() : *i;
            std::string        testPath = Path::join(dir, searchName);
            if (File::exists(testPath))
            {
                return testPath;
            }
        }
    }
    else if (!Path::isAbsolute(searchName))
    {
        // Name contains directories, but is not absolute, i.e.,
        // it is relative to the current directory.
        std::string cwd      = env.getWorkingDirectory();
        std::string testPath = Path::join(cwd, searchName);
        return testPath;
    }
    return searchName;
}

/*! \brief
 * Returns whether given path contains files from `share/top/`.
 *
 * Only checks for a single file that has an uncommon enough name.
 */
bool isAcceptableLibraryPath(const std::string &path)
{
    return Path::exists(Path::join(path, "gurgle.dat"));
}

/*! \brief
 * Returns whether given path prefix contains files from `share/top/`.
 *
 * \param[in]  path   Path prefix to check.
 * \returns  `true` if \p path contains the data files.
 *
 * Checks whether \p path could be the installation prefix where `share/top/`
 * files have been installed:  appends the relative installation path of the
 * data files and calls isAcceptableLibraryPath().
 */
bool isAcceptableLibraryPathPrefix(const std::string &path)
{
    std::string testPath = Path::join(path, DATA_INSTALL_DIR, "top");
    if (isAcceptableLibraryPath(testPath))
    {
        return true;
    }
    return false;
}

/*! \brief
 * Returns a fallback installation prefix path.
 *
 * Checks a few standard locations for the data files before returning a
 * configure-time hard-coded path.  The hard-coded path is preferred if it
 * actually contains the data files, though.
 */
std::string findFallbackInstallationPrefixPath()
{
#ifndef GMX_NATIVE_WINDOWS
    if (!isAcceptableLibraryPathPrefix(CMAKE_INSTALL_PREFIX))
    {
        if (isAcceptableLibraryPathPrefix("/usr/local"))
        {
            return "/usr/local";
        }
        if (isAcceptableLibraryPathPrefix("/usr"))
        {
            return "/usr";
        }
        if (isAcceptableLibraryPathPrefix("/opt"))
        {
            return "/opt";
        }
    }
#endif
    return CMAKE_INSTALL_PREFIX;
}

/*! \brief
 * Generic function to find data files based on path of the binary.
 *
 * \param[in]  binaryPath     Absolute path to the binary.
 * \param[out] bSourceLayout  Set to `true` if the binary is run from
 *     the build tree and the original source directory can be found.
 * \returns  Path to the `share/top/` data files.
 *
 * The search based on the path only works if the binary is in the same
 * relative path as the installed \Gromacs binaries.  If the binary is
 * somewhere else, a hard-coded fallback is used.  This doesn't work if the
 * binaries are somewhere else than the path given during configure time...
 *
 * Extra logic is present to allow running binaries from the build tree such
 * that they use up-to-date data files from the source tree.
 */
std::string findInstallationPrefixPath(const std::string &binaryPath,
                                       bool              *bSourceLayout)
{
    *bSourceLayout = false;
    // Don't search anything if binary cannot be found.
    if (Path::exists(binaryPath))
    {
        // Remove the executable name.
        std::string searchPath = Path::getParentPath(binaryPath);
        // If running directly from the build tree, try to use the source
        // directory.
#if (defined CMAKE_SOURCE_DIR && defined CMAKE_BINARY_DIR)
        std::string buildBinPath;
#ifdef CMAKE_INTDIR /*In multi-configuration build systems the output subdirectory*/
        buildBinPath = Path::join(CMAKE_BINARY_DIR, "bin", CMAKE_INTDIR);
#else
        buildBinPath = Path::join(CMAKE_BINARY_DIR, "bin");
#endif
        if (Path::isEquivalent(searchPath, buildBinPath))
        {
            std::string testPath = Path::join(CMAKE_SOURCE_DIR, "share/top");
            if (isAcceptableLibraryPath(testPath))
            {
                *bSourceLayout = true;
                return CMAKE_SOURCE_DIR;
            }
        }
#endif

        // Use the executable path to (try to) find the library dir.
        // TODO: Consider only going up exactly the required number of levels.
        while (!searchPath.empty())
        {
            if (isAcceptableLibraryPathPrefix(searchPath))
            {
                return searchPath;
            }
            searchPath = Path::getParentPath(searchPath);
        }
    }

    // End of smart searching. If we didn't find it in our parent tree,
    // or if the program name wasn't set, return a fallback.
    return findFallbackInstallationPrefixPath();
}

//! \}

}   // namespace

/********************************************************************
 * CommandLineProgramContext::Impl
 */

class CommandLineProgramContext::Impl
{
    public:
        Impl();
        Impl(int argc, const char *const argv[],
             ExecutableEnvironmentPointer env);

        /*! \brief
         * Finds the full binary path if it isn't searched yet.
         *
         * Sets \a fullBinaryPath_ if it isn't set yet.
         *
         * The \a binaryPathMutex_ should be locked by the caller before
         * calling this function.
         */
        void findBinaryPath() const;

        ExecutableEnvironmentPointer  executableEnv_;
        std::string                   invokedName_;
        std::string                   programName_;
        std::string                   displayName_;
        std::string                   commandLine_;
        mutable std::string           fullBinaryPath_;
        mutable std::string           installationPrefix_;
        mutable bool                  bSourceLayout_;
        mutable tMPI::mutex           binaryPathMutex_;
};

CommandLineProgramContext::Impl::Impl()
    : programName_("GROMACS"), bSourceLayout_(false)
{
}

CommandLineProgramContext::Impl::Impl(int argc, const char *const argv[],
                                      ExecutableEnvironmentPointer env)
    : executableEnv_(env), bSourceLayout_(false)
{
    invokedName_ = (argc != 0 ? argv[0] : "");
    programName_ = Path::getFilename(invokedName_);
    programName_ = stripSuffixIfPresent(programName_, ".exe");

    commandLine_ = quoteIfNecessary(programName_.c_str());
    for (int i = 1; i < argc; ++i)
    {
        commandLine_.append(" ");
        commandLine_.append(quoteIfNecessary(argv[i]));
    }
}

void CommandLineProgramContext::Impl::findBinaryPath() const
{
    if (fullBinaryPath_.empty())
    {
        fullBinaryPath_ = findFullBinaryPath(invokedName_, *executableEnv_);
        fullBinaryPath_ = Path::normalize(Path::resolveSymlinks(fullBinaryPath_));
        // TODO: Investigate/Consider using a dladdr()-based solution.
        // Potentially less portable, but significantly simpler, and also works
        // with user binaries even if they are located in some arbitrary location,
        // as long as shared libraries are used.
    }
}

/********************************************************************
 * CommandLineProgramContext
 */

CommandLineProgramContext::CommandLineProgramContext()
    : impl_(new Impl)
{
}

CommandLineProgramContext::CommandLineProgramContext(const char *binaryName)
    : impl_(new Impl(1, &binaryName, DefaultExecutableEnvironment::create()))
{
}

CommandLineProgramContext::CommandLineProgramContext(
        int argc, const char *const argv[])
    : impl_(new Impl(argc, argv, DefaultExecutableEnvironment::create()))
{
}

CommandLineProgramContext::CommandLineProgramContext(
        int argc, const char *const argv[], ExecutableEnvironmentPointer env)
    : impl_(new Impl(argc, argv, env))
{
}

CommandLineProgramContext::~CommandLineProgramContext()
{
}

void CommandLineProgramContext::setDisplayName(const std::string &name)
{
    GMX_RELEASE_ASSERT(impl_->displayName_.empty(),
                       "Can only set display name once");
    impl_->displayName_ = name;
}

const char *CommandLineProgramContext::programName() const
{
    return impl_->programName_.c_str();
}

const char *CommandLineProgramContext::displayName() const
{
    return impl_->displayName_.empty()
           ? impl_->programName_.c_str()
           : impl_->displayName_.c_str();
}

const char *CommandLineProgramContext::commandLine() const
{
    return impl_->commandLine_.c_str();
}

const char *CommandLineProgramContext::fullBinaryPath() const
{
    tMPI::lock_guard<tMPI::mutex> lock(impl_->binaryPathMutex_);
    impl_->findBinaryPath();
    return impl_->fullBinaryPath_.c_str();
}

InstallationPrefixInfo CommandLineProgramContext::installationPrefix() const
{
    tMPI::lock_guard<tMPI::mutex> lock(impl_->binaryPathMutex_);
    if (impl_->installationPrefix_.empty())
    {
        impl_->findBinaryPath();
        impl_->installationPrefix_ =
            Path::normalize(findInstallationPrefixPath(impl_->fullBinaryPath_,
                                                       &impl_->bSourceLayout_));
    }
    return InstallationPrefixInfo(
            impl_->installationPrefix_.c_str(),
            impl_->bSourceLayout_);
}

} // namespace gmx
