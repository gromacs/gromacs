/*
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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _oenv_h
#define _oenv_h

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0 /* avoid screwing up indentation */
}
#endif


/* output_env member functions */

/* The output_env structure holds information about program name, cmd line, 
   default times, etc. 

   There are still legacy functions for the program name, and the command 
   line, but the output_env versions are now preferred.*/

typedef enum 
{ 
    timeNULL, time_fs, time_ps, time_ns, time_us, time_ms, time_s 
} time_unit_t;
/* the time units. For the time being, ps means no conversion. */

typedef enum { exvgNULL, exvgXMGRACE, exvgXMGR, exvgNONE } xvg_format_t;
/* the xvg output formattings */


struct output_env
{
    time_unit_t time_unit; /* the time unit, enum defined in statuti.h */
    gmx_bool view;  /* view of file requested */
    xvg_format_t xvg_format; /* xvg output format, enum defined in statutil.h */
    int  verbosity; /* The level of verbosity for this program */
    int debug_level; /* the debug level */

    char *program_name; /* the program name */
    char *cmd_line; /* the re-assembled command line */
};


void output_env_init(output_env_t oenv,  int argc, char *argv[],
                     time_unit_t tmu, gmx_bool view, xvg_format_t xvg_format,
                     int verbosity, int debug_level);
/* initialize an output_env structure, setting the command line, 
   the default time value a gmx_boolean view that is set to TRUE when the 
   user requests direct viewing of graphs, 
   the graph formatting type, the verbosity, and debug level */

void output_env_init_default(output_env_t oenv);
/* initialize an output_env structure, with reasonable default settings.
    (the time unit is set to time_ps, which means no conversion).  */

void output_env_done(output_env_t oenv);
/* free memory allocated for an output_env structure. */


int output_env_get_verbosity(const output_env_t oenv);
/* return the verbosity */

int output_env_get_debug_level(const output_env_t oenv);
/* return the debug level */

const char *output_env_get_time_unit(const output_env_t oenv);
/* return time unit (e.g. ps or ns) */

const char *output_env_get_time_label(const output_env_t oenv);
/* return time unit label (e.g. "Time (ps)") */

const char *output_env_get_xvgr_tlabel(const output_env_t oenv);
/* retrun x-axis time label for xmgr */

real output_env_get_time_factor(const output_env_t oenv);
/* return time conversion factor from ps (i.e. 1e-3 for ps->ns) */

real output_env_get_time_invfactor(const output_env_t oenv);
/* return inverse time conversion factor from ps (i.e. 1e3 for ps->ns) */

real output_env_conv_time(const output_env_t oenv, real time);
/* return converted time */

void output_env_conv_times(const output_env_t oenv, int n, real *time);
/* convert array of times */

gmx_bool output_env_get_view(const output_env_t oenv);
/* Return TRUE when user requested viewing of the file */


xvg_format_t output_env_get_xvg_format(const output_env_t oenv);
/* Returns enum (see above) for xvg output formatting */

const char *output_env_get_program_name(const output_env_t oenv);
/* return the program name */

const char *output_env_get_cmd_line(const output_env_t oenv);
/* return the command line */

const char *output_env_get_short_program_name(const output_env_t oenv);
/* get the short version (without path component) of the program name */



#ifdef __cplusplus
}
#endif

#endif
