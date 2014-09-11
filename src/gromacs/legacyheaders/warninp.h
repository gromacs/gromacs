/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2010,2014, by the GROMACS development team, led by
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

#ifndef _warninp_h
#define _warninp_h

#include "gromacs/utility/basedefinitions.h"

#ifdef __cplusplus
extern "C" {
#endif


/* Abstract type for warning bookkeeping */
typedef struct warninp *warninp_t;


warninp_t
init_warning(gmx_bool bAllowWarnings, int maxwarning);
/* Initialize the warning data structure.
 * If bAllowWarnings=FALSE, all warnings (calls to warning()) will be
 * transformed into errors, calls to warning_note still produce notes.
 * maxwarning determines the maximum number of warnings that are allowed
 * for proceeding. When this number is exceeded check_warning_error
 * and done_warning will generate a fatal error.
 * bAllowWarnings=TRUE should only be used by programs that have
 * a -maxwarn command line option.
 */

void
set_warning_line(warninp_t wi, const char *fn, int line);
/* Set filename and linenumber for the warning */

int
get_warning_line(warninp_t wi);
/* Get linenumber for the warning */


const char *
get_warning_file(warninp_t wi);
/* Get filename for the warning */

void
warning(warninp_t wi, const char *s);
/* Issue a warning, with the string s. If s == NULL, then warn_buf
 * will be printed instead. The file and line set by set_warning_line
 * are printed, nwarn_warn (local) is incremented.
 * A fatal error will be generated after processing the input
 * when nwarn_warn is larger than maxwarning passed to init_warning.
 * So warning should only be called for issues that should be resolved,
 * otherwise warning_note should be called.
 */

void
warning_note(warninp_t wi, const char *s);
/* Issue a note, with the string s. If s == NULL, then warn_buf
 * will be printed instead. The file and line set by set_warning_line
 * are printed, nwarn_note (local) is incremented.
 * This is for issues which could be a problem for some systems,
 * but 100% ok for other systems.
 */

void
warning_error(warninp_t wi, const char *s);
/* Issue an error, with the string s. If s == NULL, then warn_buf
 * will be printed instead. The file and line set by set_warning_line
 * are printed, nwarn_error (local) is incremented.
 */

void
check_warning_error(warninp_t wi, int f_errno, const char *file, int line);
/* When warning_error has been called at least once gmx_fatal is called,
 * otherwise does nothing.
 */

void
done_warning(warninp_t wi, int f_errno, const char *file, int line);
/* Should be called when finished processing the input file.
 * Prints the number of notes and warnings
 * and generates a fatal error when errors were found or too many
 * warnings were generatesd.
 * Frees the data structure pointed to by wi.
 */

void
_too_few(warninp_t wi, const char *fn, int line);
#define too_few(wi) _too_few(wi, __FILE__, __LINE__)
/* Issue a warning stating 'Too few parameters' */

void
_incorrect_n_param(warninp_t wi, const char *fn, int line);
#define incorrect_n_param(wi) _incorrect_n_param(wi, __FILE__, __LINE__)
/* Issue a warning stating 'Incorrect number of parameters' */

#ifdef __cplusplus
}
#endif

#endif  /* _warninp_h */
