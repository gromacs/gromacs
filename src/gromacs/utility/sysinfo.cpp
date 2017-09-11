/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2014,2015,2016,2017, by the GROMACS development team, led by
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
 * Implements functions from sysinfo.h.
 *
 * \author Teemu Murtola <teemu.murtola@gmail.com>
 * \ingroup module_utility
 */
#include "gmxpre.h"

#include "sysinfo.h"

#include "config.h"

#include <string.h>
#include <time.h>

#include <sys/types.h>
#ifdef HAVE_SYS_TIME_H
#include <sys/time.h>
#endif
#if GMX_NATIVE_WINDOWS
#include <Windows.h>
#include <process.h>
#endif
#if HAVE_PWD_H
#include <pwd.h>
#endif
#ifdef HAVE_UNISTD_H
#include <unistd.h>
#endif

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/gmxassert.h"

namespace
{
//! Static return value for cases when a string value is not available.
const char c_unknown[] = "unknown";
} // namespace

int gmx_gethostname(char *buf, size_t len)
{
    GMX_RELEASE_ASSERT(len >= 8, "Input buffer is too short");
#if GMX_NATIVE_WINDOWS
    DWORD  dlen = len;
    if (GetComputerName(buf, &dlen))
    {
        return 0;
    }
#elif defined(HAVE_UNISTD_H) && !defined(__native_client__)
    if (gethostname(buf, len-1) == 0)
    {
        buf[len-1] = '\0';
        return 0;
    }
#endif
    strcpy(buf, c_unknown);
    return -1;
}

int gmx_getpid()
{
#if GMX_NATIVE_WINDOWS
    return _getpid();
#else
    return getpid();
#endif
}

int gmx_getuid()
{
#if defined(HAVE_UNISTD_H) && !defined(__MINGW32__)
    return getuid();
#else
    return -1;
#endif
}

int gmx_getusername(char *buf, size_t len)
{
    GMX_RELEASE_ASSERT(len >= 8, "Input buffer is too short");
    // TODO: nice_header() used getpwuid() instead; consider using getpwuid_r()
    // here.  If not, get rid of HAVE_PWD_H completely.
#if GMX_NATIVE_WINDOWS
    DWORD  dlen = len;
    if (GetUserName(buf, &dlen))
    {
        return 0;
    }
#elif defined(HAVE_UNISTD_H) && !__has_feature(memory_sanitizer) //MSAN Issue 83
    if (!getlogin_r(buf, len))
    {
        buf[len-1] = '\0';
        return 0;
    }
#endif
    strcpy(buf, c_unknown);
    return -1;
}

char *
gmx_ctime_r(const time_t *clock, char *buf, size_t len)
{
#ifdef _MSC_VER
    /* Windows */
    ctime_s(buf, len, clock);
#elif GMX_NATIVE_WINDOWS
    char *tmpbuf = ctime(clock);
    strncpy(buf, tmpbuf, len-1);
    buf[len-1] = '\0';
#elif (defined(__sun))
    /*Solaris*/
    ctime_r(clock, buf);
#else
    char tmpbuf[30];
    ctime_r(clock, tmpbuf);
    strncpy(buf, tmpbuf, len-1);
    buf[len-1] = '\0';
#endif
    return buf;
}

void gmx_format_current_time(char *buf, size_t len)
{
    time_t clock = time(nullptr);
    gmx_ctime_r(&clock, buf, len);
}

int gmx_set_nice(int level)
{
#if GMX_USE_NICE
    // TODO: This may not be reliable, but currently the return value is not
    // used.
    if (nice(level) != -1)
    {
        return 0;
    }
#else
    GMX_UNUSED_VALUE(level);
#endif
    return -1;
}
