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

#include <new>
#include <string>

#include <boost/scoped_ptr.hpp>

#include "gromacs/legacyheaders/statutil.h"

#include "gromacs/utility/path.h"
#include "gromacs/utility/stringutil.h"

namespace gmx
{

namespace
{
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
};

ProgramInfo::Impl::Impl()
    : realBinaryName_("GROMACS"), fullInvokedProgram_("GROMACS")
{
}

ProgramInfo::Impl::Impl(const char *realBinaryName,
                        int argc, const char *const argv[])
    : realBinaryName_(realBinaryName != NULL ? realBinaryName : ""),
      fullInvokedProgram_(argv[0]),
      programName_(Path::splitToPathAndFilename(fullInvokedProgram_).second),
      invariantProgramName_(programName_)
{
    invariantProgramName_ = stripSuffixIfPresent(invariantProgramName_, ".exe");
#ifdef GMX_BINARY_SUFFIX
    invariantProgramName_ =
        stripSuffixIfPresent(invariantProgramName_, GMX_BINARY_SUFFIX);
#endif
    if (realBinaryName == NULL)
    {
        realBinaryName_ = invariantProgramName_;
    }
}

/********************************************************************
 * ProgramInfo
 */

// static
const ProgramInfo &ProgramInfo::getInstance()
{
    // TODO: For completeness, this should be thread-safe.
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
        set_program_name(argv[0]);
        // TODO: Constness should not be cast away.
        set_command_line(argc, const_cast<char **>(argv));
        // TODO: For completeness, this should be thread-safe.
        g_programInfo.reset(new ProgramInfo(realBinaryName, argc, argv));
        return *g_programInfo;
    }
    catch (const std::bad_alloc &)
    {
        // TODO: Report error.
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

} // namespace gmx
