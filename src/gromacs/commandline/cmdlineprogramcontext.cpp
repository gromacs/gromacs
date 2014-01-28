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
/*! \internal \file
 * \brief
 * Implements gmx::ProgramInfo.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_commandline
 */
#include "cmdlineprogramcontext.h"

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdlib>
#include <cstring>

#include <string>
#include <vector>

#include <boost/scoped_ptr.hpp>

#include "gromacs/legacyheaders/thread_mpi/mutex.h"

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
 * constructing ProgramInfo.
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

//! \}

}   // namespace

/********************************************************************
 * ProgramInfo::Impl
 */

class ProgramInfo::Impl
{
    public:
        Impl();
        Impl(int argc, const char *const argv[],
             ExecutableEnvironmentPointer env);

        ExecutableEnvironmentPointer  executableEnv_;
        std::string                   invokedName_;
        std::string                   programName_;
        std::string                   displayName_;
        std::string                   commandLine_;
        mutable std::string           fullBinaryPath_;
        mutable tMPI::mutex           binaryPathMutex_;
};

ProgramInfo::Impl::Impl()
    : programName_("GROMACS")
{
}

ProgramInfo::Impl::Impl(int argc, const char *const argv[],
                        ExecutableEnvironmentPointer env)
    : executableEnv_(move(env))
{
    invokedName_          = (argc != 0 ? argv[0] : "");
    programName_          = Path::splitToPathAndFilename(invokedName_).second;
    programName_          = stripSuffixIfPresent(programName_, ".exe");

    commandLine_ = quoteIfNecessary(programName_.c_str());
    for (int i = 1; i < argc; ++i)
    {
        commandLine_.append(" ");
        commandLine_.append(quoteIfNecessary(argv[i]));
    }
}

/********************************************************************
 * ProgramInfo
 */

ProgramInfo::ProgramInfo()
    : impl_(new Impl)
{
}

ProgramInfo::ProgramInfo(const char *binaryName)
    : impl_(new Impl(1, &binaryName,
                     DefaultExecutableEnvironment::create()))
{
}

ProgramInfo::ProgramInfo(int argc, const char *const argv[])
    : impl_(new Impl(argc, argv,
                     DefaultExecutableEnvironment::create()))
{
}

ProgramInfo::ProgramInfo(int argc, const char *const argv[],
                         ExecutableEnvironmentPointer env)
    : impl_(new Impl(argc, argv, move(env)))
{
}

ProgramInfo::~ProgramInfo()
{
}

void ProgramInfo::setDisplayName(const std::string &name)
{
    GMX_RELEASE_ASSERT(impl_->displayName_.empty(),
                       "Can only set display name once");
    impl_->displayName_ = name;
}

const char *ProgramInfo::programName() const
{
    return impl_->programName_.c_str();
}

const char *ProgramInfo::displayName() const
{
    return impl_->displayName_.empty()
           ? impl_->programName_.c_str()
           : impl_->displayName_.c_str();
}

const char *ProgramInfo::commandLine() const
{
    return impl_->commandLine_.c_str();
}

const char *ProgramInfo::fullBinaryPath() const
{
    tMPI::lock_guard<tMPI::mutex> lock(impl_->binaryPathMutex_);
    if (impl_->fullBinaryPath_.empty())
    {
        impl_->fullBinaryPath_ =
            Path::normalize(
                    Path::resolveSymlinks(
                            findFullBinaryPath(impl_->invokedName_,
                                               *impl_->executableEnv_)));
        // TODO: Investigate/Consider using a dladdr()-based solution.
        // Potentially less portable, but significantly simpler, and also works
        // with user binaries even if they are located in some arbitrary location,
        // as long as shared libraries are used.
    }
    return impl_->fullBinaryPath_.c_str();
}

} // namespace gmx
