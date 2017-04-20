/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "oenv.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

struct gmx_output_env_t
{
    explicit gmx_output_env_t(const gmx::IProgramContext &context)
        : programContext(context),
          time_unit(time_ps),
          view(FALSE),
          xvg_format(exvgNONE),
          verbosity(0) {}

    const gmx::IProgramContext  &programContext;

    /* the time unit, enum defined in oenv.h */
    time_unit_t                          time_unit;
    /* view of file requested */
    gmx_bool                             view;
    /* xvg output format, enum defined in oenv.h */
    xvg_format_t                         xvg_format;
    /* The level of verbosity for this program */
    int                                  verbosity;
};

/* The source code in this file should be thread-safe.
      Please keep it that way. */

/******************************************************************
 *
 *             T R A J E C T O R Y   S T U F F
 *
 ******************************************************************/

/* read only time names */
/* These must correspond to the time units type time_unit_t in oenv.h */
static const real  timefactors[] =   { real(0),  real(1e3),  real(1), real(1e-3), real(1e-6), real(1e-9), real(1e-12), real(0) };
static const real  timeinvfactors[] = { real(0), real(1e-3),  real(1),  real(1e3),  real(1e6),  real(1e9),  real(1e12), real(0) };
static const char *time_units_str[] = {
    nullptr, "fs", "ps", "ns", "us",
    "\\mus", "ms", "s"
};
static const char *time_units_xvgr[] = {
    nullptr, "fs", "ps", "ns",
    "ms", "s", nullptr
};


/***** OUTPUT_ENV MEMBER FUNCTIONS ******/

void output_env_init(gmx_output_env_t **oenvp,
                     const gmx::IProgramContext &context,
                     time_unit_t tmu, gmx_bool view, xvg_format_t xvg_format,
                     int verbosity)
{
    try
    {
        gmx_output_env_t *oenv = new gmx_output_env_t(context);
        *oenvp            = oenv;
        oenv->time_unit   = tmu;
        oenv->view        = view;
        oenv->xvg_format  = xvg_format;
        oenv->verbosity   = verbosity;
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void output_env_init_default(gmx_output_env_t **oenvp)
{
    try
    {
        gmx_output_env_t *oenv = new gmx_output_env_t(gmx::getProgramContext());
        *oenvp = oenv;
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void output_env_done(gmx_output_env_t *oenv)
{
    delete oenv;
}


int output_env_get_verbosity(const gmx_output_env_t *oenv)
{
    return oenv->verbosity;
}

std::string output_env_get_time_unit(const gmx_output_env_t *oenv)
{
    return time_units_str[oenv->time_unit];
}

std::string output_env_get_time_label(const gmx_output_env_t *oenv)
{
    return gmx::formatString("Time (%s)", time_units_str[oenv->time_unit] ?
                             time_units_str[oenv->time_unit] : "ps");
}

std::string output_env_get_xvgr_tlabel(const gmx_output_env_t *oenv)
{
    return gmx::formatString("Time (%s)", time_units_xvgr[oenv->time_unit] ?
                             time_units_xvgr[oenv->time_unit] : "ps");
}

real output_env_get_time_factor(const gmx_output_env_t *oenv)
{
    return timefactors[oenv->time_unit];
}

real output_env_get_time_invfactor(const gmx_output_env_t *oenv)
{
    return timeinvfactors[oenv->time_unit];
}

real output_env_conv_time(const gmx_output_env_t *oenv, real time)
{
    return time*timefactors[oenv->time_unit];
}

void output_env_conv_times(const gmx_output_env_t *oenv, int n, real *time)
{
    int    i;
    double fact = timefactors[oenv->time_unit];

    if (fact != 1.)
    {
        for (i = 0; i < n; i++)
        {
            time[i] *= fact;
        }
    }
}

gmx_bool output_env_get_view(const gmx_output_env_t *oenv)
{
    return oenv->view;
}

xvg_format_t output_env_get_xvg_format(const gmx_output_env_t *oenv)
{
    return oenv->xvg_format;
}

const char *output_env_get_program_display_name(const gmx_output_env_t *oenv)
{
    const char *displayName = nullptr;

    try
    {
        displayName = oenv->programContext.displayName();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;

    return displayName;
}

const gmx::IProgramContext &
output_env_get_program_context(const gmx_output_env_t *oenv)
{
    return oenv->programContext;
}
