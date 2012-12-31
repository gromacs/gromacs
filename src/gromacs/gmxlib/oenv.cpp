/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 *
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 *          GROningen MAchine for Chemical Simulations
 *
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.

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
 *
 * And Hey:
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#include "oenv.h"

#include "smalloc.h"

#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/programinfo.h"

struct output_env
{
    output_env()
    {
        setDefaults();
    }
    output_env(int argc, const char *const argv[])
        : programInfo(argc, argv)
    {
        setDefaults();
    }

    void setDefaults()
    {
        time_unit   = time_ps;
        view        = FALSE;
        xvg_format  = exvgNONE;
        verbosity   = 0;
        debug_level = 0;
    }

    gmx::ProgramInfo programInfo;

    time_unit_t time_unit; /* the time unit, enum defined in oenv.h */
    gmx_bool view;  /* view of file requested */
    xvg_format_t xvg_format; /* xvg output format, enum defined in oenv.h */
    int  verbosity; /* The level of verbosity for this program */
    int debug_level; /* the debug level */
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
static const real timefactors[] =   { 0,  1e3,  1, 1e-3, 1e-6, 1e-9, 1e-12, 0 };
static const real timeinvfactors[] = { 0, 1e-3,  1,  1e3,  1e6,  1e9,  1e12, 0 };
static const char *time_units_str[] = { NULL, "fs", "ps", "ns", "us",
                                        "\\mus", "ms", "s"
                                      };
static const char *time_units_xvgr[] = { NULL, "fs", "ps", "ns",
                                         "ms", "s", NULL
                                       };


/***** OUTPUT_ENV MEMBER FUNCTIONS ******/

void output_env_init(output_env_t *oenvp, int argc, char *argv[],
                     time_unit_t tmu, gmx_bool view, xvg_format_t xvg_format,
                     int verbosity, int debug_level)
{
    try
    {
        output_env_t oenv = new output_env(argc, argv);
        *oenvp = oenv;
        oenv->time_unit   = tmu;
        oenv->view        = view;
        oenv->xvg_format  = xvg_format;
        oenv->verbosity   = verbosity;
        oenv->debug_level = debug_level;
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void output_env_init_default(output_env_t *oenvp)
{
    try
    {
        output_env_t oenv = new output_env();
        *oenvp = oenv;
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

void output_env_done(output_env_t oenv)
{
    delete oenv;
}


int output_env_get_verbosity(const output_env_t oenv)
{
    return oenv->verbosity;
}

int output_env_get_debug_level(const output_env_t oenv)
{
    return oenv->debug_level;
}

const char *output_env_get_time_unit(const output_env_t oenv)
{
    return time_units_str[oenv->time_unit];
}

const char *output_env_get_time_label(const output_env_t oenv)
{
    char *label;
    snew(label, 20);

    sprintf(label,"Time (%s)",time_units_str[oenv->time_unit] ?
            time_units_str[oenv->time_unit]: "ps");

    return label;
}

const char *output_env_get_xvgr_tlabel(const output_env_t oenv)
{
    char *label;
    snew(label, 20);

    sprintf(label,"Time (%s)", time_units_xvgr[oenv->time_unit] ?
            time_units_xvgr[oenv->time_unit] : "ps");

    return label;
}

real output_env_get_time_factor(const output_env_t oenv)
{
    return timefactors[oenv->time_unit];
}

real output_env_get_time_invfactor(const output_env_t oenv)
{
    return timeinvfactors[oenv->time_unit];
}

real output_env_conv_time(const output_env_t oenv, real time)
{
    return time*timefactors[oenv->time_unit];
}

void output_env_conv_times(const output_env_t oenv, int n, real *time)
{
    int i;
    double fact=timefactors[oenv->time_unit];

    if (fact!=1.)
        for (i=0; i<n; i++)
        {
            time[i] *= fact;
        }
}

gmx_bool output_env_get_view(const output_env_t oenv)
{
    return oenv->view;
}

xvg_format_t output_env_get_xvg_format(const output_env_t oenv)
{
    return oenv->xvg_format;
}

const char *output_env_get_program_name(const output_env_t oenv)
{
    try
    {
        return oenv->programInfo.programNameWithPath().c_str();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

const char *output_env_get_short_program_name(const output_env_t oenv)
{
    try
    {
        return oenv->programInfo.programName().c_str();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}

const char *output_env_get_cmd_line(const output_env_t oenv)
{
    try
    {
        return oenv->programInfo.commandLine().c_str();
    }
    GMX_CATCH_ALL_AND_EXIT_WITH_FATAL_ERROR;
}
