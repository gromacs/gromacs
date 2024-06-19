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

#include <filesystem>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/gmxassert.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

#include "buildinfo.h"

namespace gmx
{

namespace
{

//! \addtogroup module_commandline
//! \{

/*! \brief
 * Quotes a string if it contains spaces.
 */
std::string quoteIfNecessary(const char* str)
{
    const bool bSpaces = (std::strchr(str, ' ') != nullptr);
    if (bSpaces)
    {
        return formatString("'%s'", str);
    }
    return str;
}

/*! \brief
 * Default implementation for IExecutableEnvironment.
 *
 * Used if IExecutableEnvironment is not explicitly provided when
 * constructing CommandLineProgramContext.
 */
class DefaultExecutableEnvironment : public IExecutableEnvironment
{
public:
    //! Allocates a default environment.
    static ExecutableEnvironmentPointer create()
    {
        return ExecutableEnvironmentPointer(new DefaultExecutableEnvironment());
    }

    DefaultExecutableEnvironment() : initialWorkingDirectory_(std::filesystem::current_path()) {}

    std::filesystem::path getWorkingDirectory() const override { return initialWorkingDirectory_; }
    std::vector<std::filesystem::path> getExecutablePaths() const override
    {
        return getSystemExecutablePaths();
    }

private:
    std::filesystem::path initialWorkingDirectory_;
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
std::filesystem::path findFullBinaryPath(const std::filesystem::path&  invokedName,
                                         const IExecutableEnvironment& env)
{
    std::filesystem::path searchName{ invokedName };
    // On Windows & Cygwin we need to add the .exe extension,
    // or we wont be able to detect that the file exists.
#if GMX_NATIVE_WINDOWS || GMX_CYGWIN
    if (!searchName.has_extension() || searchName.extension() != ".exe")
    {
        searchName.replace_extension(".exe");
    }
#endif
    if (!searchName.has_parent_path())
    {
        // No directory in name means it must be in the path - search it!
        std::vector<std::filesystem::path> pathEntries = env.getExecutablePaths();
        for (const auto& path : pathEntries)
        {
            auto dir = path.empty() ? env.getWorkingDirectory() : path;
            dir.append(searchName.string());
            if (File::exists(dir, File::returnFalseOnError))
            {
                return dir;
            }
        }
    }
    else if (!searchName.is_absolute())
    {
        // Name contains directories, but is not absolute, i.e.,
        // it is relative to the current directory.
        return env.getWorkingDirectory().append(searchName.string());
    }
    return searchName;
}

/*! \brief
 * Returns whether given path contains files from `share/top/`.
 *
 * Only checks for a single file that has an uncommon enough name.
 */
bool isAcceptableLibraryPath(const std::filesystem::path& path)
{
    return std::filesystem::exists(std::filesystem::path(path).append("residuetypes.dat"));
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
bool isAcceptableLibraryPathPrefix(const std::filesystem::path& path)
{
    return isAcceptableLibraryPath(
            std::filesystem::path(path).append(GMX_INSTALL_GMXDATADIR).append("top"));
}

/*! \brief
 * Returns a fallback installation prefix path.
 *
 * Checks a few standard locations for the data files before returning a
 * configure-time hard-coded path.  The hard-coded path is preferred if it
 * actually contains the data files, though.
 */
std::filesystem::path findFallbackInstallationPrefixPath()
{
#if !GMX_NATIVE_WINDOWS
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
std::filesystem::path findInstallationPrefixPath(const std::filesystem::path& binaryPath, bool* bSourceLayout)
{
    *bSourceLayout = false;
    // Don't search anything if binary cannot be found.
    if (std::filesystem::exists(binaryPath))
    {
        // Remove the executable name.
        auto searchPath = binaryPath.parent_path();
        // If running directly from the build tree, try to use the source
        // directory.
#if (defined CMAKE_SOURCE_DIR && defined CMAKE_BINARY_DIR)
        std::filesystem::path buildBinPath;
#    ifdef CMAKE_INTDIR /*In multi-configuration build systems the output subdirectory*/
        buildBinPath = std::filesystem::path(CMAKE_BINARY_DIR).append("bin").append(CMAKE_INTDIR);
#    else
        buildBinPath = std::filesystem::path(CMAKE_BINARY_DIR).append("bin");
#    endif
        // If either path does not exist than an error is produced by
        // equivalent(). When an error is produced, the non-throwing
        // call to equivalent() returns false, which is what we want.
        // We don't care what the error when we are merely finding a
        // valid GROMACS installation prefix.
        if (std::error_code c; std::filesystem::equivalent(searchPath, buildBinPath, c))
        {
            if (isAcceptableLibraryPath(std::filesystem::path(CMAKE_SOURCE_DIR).append("share/top")))
            {
                *bSourceLayout = true;
                return CMAKE_SOURCE_DIR;
            }
        }
#endif

        // Use the executable path to (try to) find the library dir.
        // TODO: Consider only going up exactly the required number of levels.
        while (searchPath != searchPath.root_path()) // While we can go one level up
        {
            if (isAcceptableLibraryPathPrefix(searchPath))
            {
                return searchPath;
            }
            searchPath = searchPath.parent_path();
        }
    }

    // End of smart searching. If we didn't find it in our parent tree,
    // or if the program name wasn't set, return a fallback.
    return findFallbackInstallationPrefixPath();
}

//! \}

} // namespace

/********************************************************************
 * CommandLineProgramContext::Impl
 */

class CommandLineProgramContext::Impl
{
public:
    Impl();
    Impl(int argc, const char* const argv[], ExecutableEnvironmentPointer env);

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
    std::filesystem::path         invokedName_;
    std::string                   programName_;
    std::string                   displayName_;
    std::string                   commandLine_;
    mutable std::filesystem::path fullBinaryPath_;
    mutable std::filesystem::path installationPrefix_;
    mutable bool                  bSourceLayout_;
    mutable std::mutex            binaryPathMutex_;
};

CommandLineProgramContext::Impl::Impl() : programName_("GROMACS"), bSourceLayout_(false) {}

CommandLineProgramContext::Impl::Impl(int argc, const char* const argv[], ExecutableEnvironmentPointer env) :
    executableEnv_(std::move(env)), invokedName_(argc != 0 ? argv[0] : ""), bSourceLayout_(false)
{
    programName_ = invokedName_.filename().string();
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
        if (std::filesystem::is_symlink(fullBinaryPath_))
        {
            auto tempPath = std::filesystem::read_symlink(fullBinaryPath_);
            // if the path from the symlink is not an absolute path, we need
            // to later combine the path information again from the binaryPath
            // and the relative position of the resolved symlink
            if (!tempPath.is_absolute())
            {
                fullBinaryPath_ = fullBinaryPath_.parent_path().append(tempPath.string());
            }
            else
            {
                fullBinaryPath_ = tempPath;
            }
        }
        fullBinaryPath_ = fullBinaryPath_.make_preferred();
        // TODO: Investigate/Consider using a dladdr()-based solution.
        // Potentially less portable, but significantly simpler, and also works
        // with user binaries even if they are located in some arbitrary location,
        // as long as shared libraries are used.
    }
}

/********************************************************************
 * CommandLineProgramContext
 */

CommandLineProgramContext::CommandLineProgramContext() : impl_(new Impl) {}

CommandLineProgramContext::CommandLineProgramContext(const char* binaryName) :
    impl_(new Impl(1, &binaryName, DefaultExecutableEnvironment::create()))
{
}

CommandLineProgramContext::CommandLineProgramContext(int argc, const char* const argv[]) :
    impl_(new Impl(argc, argv, DefaultExecutableEnvironment::create()))
{
}

CommandLineProgramContext::CommandLineProgramContext(int                          argc,
                                                     const char* const            argv[],
                                                     ExecutableEnvironmentPointer env) :
    impl_(new Impl(argc, argv, std::move(env)))
{
}

CommandLineProgramContext::~CommandLineProgramContext() {}

void CommandLineProgramContext::setDisplayName(const std::string& name)
{
    GMX_RELEASE_ASSERT(impl_->displayName_.empty(), "Can only set display name once");
    impl_->displayName_ = name;
}

const char* CommandLineProgramContext::programName() const
{
    return impl_->programName_.c_str();
}

const char* CommandLineProgramContext::displayName() const
{
    return impl_->displayName_.empty() ? impl_->programName_.c_str() : impl_->displayName_.c_str();
}

const char* CommandLineProgramContext::commandLine() const
{
    return impl_->commandLine_.c_str();
}

std::filesystem::path CommandLineProgramContext::fullBinaryPath() const
{
    std::lock_guard<std::mutex> lock(impl_->binaryPathMutex_);
    impl_->findBinaryPath();
    return impl_->fullBinaryPath_.c_str();
}

InstallationPrefixInfo CommandLineProgramContext::installationPrefix() const
{
    std::lock_guard<std::mutex> lock(impl_->binaryPathMutex_);
    if (impl_->installationPrefix_.empty())
    {
        impl_->findBinaryPath();
        impl_->installationPrefix_ =
                findInstallationPrefixPath(impl_->fullBinaryPath_, &impl_->bSourceLayout_).make_preferred();
    }
    return InstallationPrefixInfo(impl_->installationPrefix_.c_str(), impl_->bSourceLayout_);
}

} // namespace gmx
