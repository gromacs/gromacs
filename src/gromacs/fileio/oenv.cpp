/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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
#include "gmxpre.h"

#include "gromacs/fileio/oenv.h"

#include <cstdlib>

#include <string>

#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/enumerationhelpers.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t
{
    explicit gmx_output_env_t(const gmx::IProgramContext& context) :
        programContext(context),
        timeUnit(gmx::TimeUnit::Picoseconds),
        view(FALSE),
        xvgFormat(XvgFormat::None),
        verbosity(0),
        trajectory_io_verbosity(0)
    {
    }


    const gmx::IProgramContext& programContext;

    /* the time unit, enum defined in timeunitmanager.h */
    gmx::TimeUnit timeUnit;
    /* view of file requested */
    gmx_bool view;
    /* xvg output format, enum defined in oenv.h */
    XvgFormat xvgFormat;
    /* The level of verbosity for this program */
    int verbosity;
    /* The level of verbosity during trajectory I/O. Default=1, quiet=0. */
    int trajectory_io_verbosity;
};

/* The source code in this file should be thread-safe.
      Please keep it that way. */

/******************************************************************
 *
 *             T R A J E C T O R Y   S T U F F
 *
 ******************************************************************/

/* read only time names */
static const gmx::EnumerationArray<gmx::TimeUnit, real> c_picosecondsInTimeUnit = {
    { real(1e3), real(1), real(1e-3), real(1e-6), real(1e-9), real(1e-12) }
};
static const gmx::EnumerationArray<gmx::TimeUnit, real> c_timeUnitsInPicoseconds = {
    { real(1e-3), real(1), real(1e3), real(1e6), real(1e9), real(1e12) }
};
static const gmx::EnumerationArray<gmx::TimeUnit, const char*> c_timeUnitNames = {
    { "fs", "ps", "ns", "us", "ms", "s" }
};
static const gmx::EnumerationArray<gmx::TimeUnit, const char*> c_timeUnitNamesForXvgr = {
    { "fs", "ps", "ns", "\\mus", "ms", "s" }
};


/***** OUTPUT_ENV MEMBER FUNCTIONS ******/

void output_env_init(gmx_output_env_t**          oenvp,
                     const gmx::IProgramContext& context,
                     gmx::TimeUnit               tmu,
                     gmx_bool                    view,
                     XvgFormat                   xvgFormat,
                     int                         verbosity)
{
    try
    {
        gmx_output_env_t* oenv        = new gmx_output_env_t(context);
        *oenvp                        = oenv;
        oenv->timeUnit                = tmu;
        oenv->view                    = view;
        oenv->xvgFormat               = xvgFormat;
        oenv->verbosity               = verbosity;
        const char* env               = getenv("GMX_TRAJECTORY_IO_VERBOSITY");
        oenv->trajectory_io_verbosity = (env != nullptr ? strtol(env, nullptr, 10) : 1);
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
}

void output_env_init_default(gmx_output_env_t** oenvp)
{
    try
    {
        gmx_output_env_t* oenv = new gmx_output_env_t(gmx::getProgramContext());
        *oenvp                 = oenv;
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR
}

void output_env_done(gmx_output_env_t* oenv)
{
    delete oenv;
}


int output_env_get_verbosity(const gmx_output_env_t* oenv)
{
    return oenv->verbosity;
}

int output_env_get_trajectory_io_verbosity(const gmx_output_env_t* oenv)
{
    return oenv->trajectory_io_verbosity;
}

std::string output_env_get_time_unit(const gmx_output_env_t* oenv)
{
    return c_timeUnitNames[oenv->timeUnit];
}

std::string output_env_get_time_label(const gmx_output_env_t* oenv)
{
    return gmx::formatString("Time (%s)", c_timeUnitNames[oenv->timeUnit]);
}

std::string output_env_get_xvgr_tlabel(const gmx_output_env_t* oenv)
{
    return gmx::formatString("Time (%s)", c_timeUnitNamesForXvgr[oenv->timeUnit]);
}

real output_env_get_time_factor(const gmx_output_env_t* oenv)
{
    return c_picosecondsInTimeUnit[oenv->timeUnit];
}

real output_env_get_time_invfactor(const gmx_output_env_t* oenv)
{
    return c_timeUnitsInPicoseconds[oenv->timeUnit];
}

real output_env_conv_time(const gmx_output_env_t* oenv, real time)
{
    return time * c_picosecondsInTimeUnit[oenv->timeUnit];
}

void output_env_conv_times(const gmx_output_env_t* oenv, int n, real* time)
{
    int    i;
    double fact = c_picosecondsInTimeUnit[oenv->timeUnit];

    if (fact != 1.)
    {
        for (i = 0; i < n; i++)
        {
            time[i] *= fact;
        }
    }
}

gmx_bool output_env_get_view(const gmx_output_env_t* oenv)
{
    return oenv->view;
}

XvgFormat output_env_get_xvg_format(const gmx_output_env_t* oenv)
{
    return oenv->xvgFormat;
}

const char* output_env_get_program_display_name(const gmx_output_env_t* oenv)
{
    const char* displayName = nullptr;

    try
    {
        displayName = oenv->programContext.displayName();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR

    return displayName;
}

const gmx::IProgramContext& output_env_get_program_context(const gmx_output_env_t* oenv)
{
    return oenv->programContext;
}
