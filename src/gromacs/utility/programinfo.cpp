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
tMPI::mutex g_programInfoMutex;
boost::scoped_ptr<ProgramInfo> g_programInfo;
} // namespace

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
        bool        bIncomplete_;
};

ProgramInfo::Impl::Impl()
    : realBinaryName_("GROMACS"), fullInvokedProgram_("GROMACS"),
      bIncomplete_(true)
{
}

ProgramInfo::Impl::Impl(const char *realBinaryName,
                        int argc, const char *const argv[])
    : realBinaryName_(realBinaryName != NULL ? realBinaryName : ""),
      fullInvokedProgram_(argc != 0 ? argv[0] : ""),
      programName_(Path::splitToPathAndFilename(fullInvokedProgram_).second),
      invariantProgramName_(programName_),
      bIncomplete_(argc <= 0)
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
    invariantProgramName_ = stripSuffixIfPresent(invariantProgramName_, ".exe");
#ifdef GMX_BINARY_SUFFIX
    invariantProgramName_ =
        stripSuffixIfPresent(invariantProgramName_, GMX_BINARY_SUFFIX);
#endif
    if (realBinaryName == NULL)
    {
        realBinaryName_ = invariantProgramName_;
    }

    // TODO: Remove this hack with negative argc once there is no need for
    // set_program_name().
    if (argc < 0)
    {
        argc = -argc;
    }
    for (int i = 0; i < argc; ++i)
    {
        if (i > 0)
        {
            commandLine_.append(" ");
        }
        const char *arg = argv[i];
        bool bSpaces = (std::strchr(arg, ' ') != NULL);
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
        g_programInfo.reset(new ProgramInfo());
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
        if (g_programInfo.get() == NULL || g_programInfo->impl_->bIncomplete_)
        {
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

std::string ProgramInfo::realBinaryName() const
{
    return impl_->realBinaryName_;
}

std::string ProgramInfo::programNameWithPath() const
{
    return impl_->fullInvokedProgram_;
}

std::string ProgramInfo::programName() const
{
    return impl_->programName_;
}

std::string ProgramInfo::invariantProgramName() const
{
    return impl_->invariantProgramName_;
}

std::string ProgramInfo::commandLine() const
{
    return impl_->commandLine_;
}

} // namespace gmx
