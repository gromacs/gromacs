/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012, by the GROMACS development team, led by
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
/*! \internal \file
 * \brief
 * Implements gmx::ProgramInfo.
 *
 * \author Teemu Murtola <teemu.murtola@cbr.su.se>
 * \ingroup module_utility
 */
#include "programinfo.h"

// For GMX_BINARY_SUFFIX
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <cstdlib>
#include <cstring>

#include <algorithm>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "gromacs/legacyheaders/futil.h"
#include "gromacs/legacyheaders/thread_mpi/mutex.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{
//! Mutex for updates to the global program info objects.
tMPI::mutex                    g_programInfoMutex;
//! Partially filled program info, needed to support set_program_name().
boost::scoped_ptr<ProgramInfo> g_partialProgramInfo;
//! Global program info; stores the object initialized with ProgramInfo::init().
boost::scoped_ptr<ProgramInfo> g_programInfo;
}   // namespace

/********************************************************************
 * ProgramInfo::Impl
 */

class ProgramInfo::Impl
{
    public:
        Impl();
        Impl(const char *realBinaryName, int argc, const char *const argv[]);

        std::string realBinaryName_;
        std::string fullInvokedProgram_;
        std::string programName_;
        std::string invariantProgramName_;
        std::string commandLine_;
};

ProgramInfo::Impl::Impl()
    : realBinaryName_("GROMACS"), fullInvokedProgram_("GROMACS"),
      programName_("GROMACS"), invariantProgramName_("GROMACS")
{
}

ProgramInfo::Impl::Impl(const char *realBinaryName,
                        int argc, const char *const argv[])
    : realBinaryName_(realBinaryName != NULL ? realBinaryName : ""),
      fullInvokedProgram_(argc != 0 ? argv[0] : ""),
      programName_(Path::splitToPathAndFilename(fullInvokedProgram_).second)
{
    // Temporary hack to make things work on Windows while waiting for #950.
    // Some places in the existing code expect to have DIR_SEPARATOR in all
    // input paths, but Windows may also give '/' (and does that, e.g., for
    // tests invoked through CTest).
    // When removing this, remove also the #include "futil.h".
    if (DIR_SEPARATOR == '\\')
    {
        std::replace(fullInvokedProgram_.begin(), fullInvokedProgram_.end(),
                     '/', '\\');
    }
    programName_          = stripSuffixIfPresent(programName_, ".exe");
    invariantProgramName_ = programName_;
#ifdef GMX_BINARY_SUFFIX
    invariantProgramName_ =
        stripSuffixIfPresent(invariantProgramName_, GMX_BINARY_SUFFIX);
#endif
    if (realBinaryName == NULL)
    {
        realBinaryName_ = invariantProgramName_;
    }

    for (int i = 0; i < argc; ++i)
    {
        if (i > 0)
        {
            commandLine_.append(" ");
        }
        const char *arg     = argv[i];
        bool        bSpaces = (std::strchr(arg, ' ') != NULL);
        if (bSpaces)
        {
            commandLine_.append("'");
        }
        commandLine_.append(arg);
        if (bSpaces)
        {
            commandLine_.append("'");
        }
    }
}

/********************************************************************
 * ProgramInfo
 */

// static
const ProgramInfo &ProgramInfo::getInstance()
{
    tMPI::lock_guard<tMPI::mutex> lock(g_programInfoMutex);
    if (g_programInfo.get() == NULL)
    {
        if (g_partialProgramInfo.get() != NULL)
        {
            return *g_partialProgramInfo;
        }
        static ProgramInfo fallbackInfo;
        return fallbackInfo;
    }
    return *g_programInfo;
}

// static
const ProgramInfo &ProgramInfo::init(int argc, const char *const argv[])
{
    return init(NULL, argc, argv);
}

// static
const ProgramInfo &ProgramInfo::init(const char *realBinaryName,
                                     int argc, const char *const argv[])
{
    try
    {
        tMPI::lock_guard<tMPI::mutex> lock(g_programInfoMutex);
        if (g_programInfo.get() == NULL)
        {
            // TODO: Remove this hack with negative argc once there is no need for
            // set_program_name().
            if (argc < 0)
            {
                if (g_partialProgramInfo.get() == NULL)
                {
                    g_partialProgramInfo.reset(
                            new ProgramInfo(realBinaryName, -argc, argv));
                }
                return *g_partialProgramInfo;
            }
            g_programInfo.reset(new ProgramInfo(realBinaryName, argc, argv));
        }
        return *g_programInfo;
    }
    catch (const std::exception &ex)
    {
        printFatalErrorMessage(stderr, ex);
        std::exit(1);
    }
}

ProgramInfo::ProgramInfo()
    : impl_(new Impl)
{
}

ProgramInfo::ProgramInfo(const char *realBinaryName)
    : impl_(new Impl(realBinaryName, 1, &realBinaryName))
{
}

ProgramInfo::ProgramInfo(int argc, const char *const argv[])
    : impl_(new Impl(NULL, argc, argv))
{
}

ProgramInfo::ProgramInfo(const char *realBinaryName,
                         int argc, const char *const argv[])
    : impl_(new Impl(realBinaryName, argc, argv))
{
}

ProgramInfo::~ProgramInfo()
{
}

const std::string &ProgramInfo::realBinaryName() const
{
    return impl_->realBinaryName_;
}

const std::string &ProgramInfo::programNameWithPath() const
{
    return impl_->fullInvokedProgram_;
}

const std::string &ProgramInfo::programName() const
{
    return impl_->programName_;
}

const std::string &ProgramInfo::invariantProgramName() const
{
    return impl_->invariantProgramName_;
}

const std::string &ProgramInfo::commandLine() const
{
    return impl_->commandLine_;
}

} // namespace gmx
