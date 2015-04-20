/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2014, by the GROMACS development team, led by
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

#ifndef _oenv_h
#define _oenv_h

#include "gromacs/legacyheaders/types/oenv.h"
#include "gromacs/legacyheaders/types/simple.h"

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

typedef enum {
    exvgNULL, exvgXMGRACE, exvgXMGR, exvgNONE
} xvg_format_t;
/* the xvg output formattings */


void output_env_init_default(output_env_t *oenvp);
/* initialize an output_env structure, with reasonable default settings.
    (the time unit is set to time_ps, which means no conversion).  */

void output_env_done(output_env_t oenv);
/* free memory allocated for an output_env structure. */


int output_env_get_verbosity(const output_env_t oenv);
/* return the verbosity */

const char *output_env_get_time_unit(const output_env_t oenv);
/* return time unit (e.g. ps or ns) */

const char *output_env_get_time_label(const output_env_t oenv);
/* return time unit label (e.g. "Time (ps)") */

const char *output_env_get_xvgr_tlabel(const output_env_t oenv);
/* retrun x-axis time label for xmgr */

real output_env_get_time_factor(const output_env_t oenv);
/* return time conversion factor from ps (i.e. 1e-3 for ps->ns) */

real output_env_get_time_invfactor(const output_env_t oenv);
/* return inverse time conversion factor to ps (i.e. 1e3 for ns->ps) */

real output_env_conv_time(const output_env_t oenv, real time);
/* return converted time */

void output_env_conv_times(const output_env_t oenv, int n, real *time);
/* convert array of times */

gmx_bool output_env_get_view(const output_env_t oenv);
/* Return TRUE when user requested viewing of the file */

xvg_format_t output_env_get_xvg_format(const output_env_t oenv);
/* Returns enum (see above) for xvg output formatting */

/*! \brief
 * Returns display name for the currently running program.
 */
const char *output_env_get_program_display_name(const output_env_t oenv);

#ifdef __cplusplus
}

namespace gmx
{
class ProgramContextInterface;
} // namespace gmx

void output_env_init(output_env_t *oenvp,
                     const gmx::ProgramContextInterface &context,
                     time_unit_t tmu, gmx_bool view, xvg_format_t xvg_format,
                     int verbosity);
/* initialize an output_env structure, setting the command line,
   the default time value a gmx_boolean view that is set to TRUE when the
   user requests direct viewing of graphs,
   the graph formatting type, the verbosity, and debug level */

/*! \brief
 * Returns gmx::ProgramContextInterface from an output_env structure.
 */
const gmx::ProgramContextInterface &
output_env_get_program_context(const output_env_t oenv);

#endif

#endif
