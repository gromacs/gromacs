/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
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
/*! \file
 * \brief Generic string handling functions.
 */
#ifndef _string2_h
#define _string2_h

/*
 *
 * string2.h
 * David van der Spoel
 *
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>
#include "visibility.h"
#include "gmx_header_config.h"

/*#include "typedefs.h"*/
#include "types/simple.h"

/* Suppress Cygwin compiler warnings from using newlib version of
 * ctype.h */
#ifdef GMX_CYGWIN
#undef isdigit
#undef isstring
#undef isspace
#undef isalnum
#undef isalpha
#undef ispunct
#undef isxdigit
#undef isupper
#undef islower
#undef toupper
#undef tolower
#endif

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

#define CONTINUE    '\\'
#define COMMENTSIGN ';'

GMX_LIBGMX_EXPORT
int continuing(char *s);

GMX_LIBGMX_EXPORT
char *fgets2(char *s, int n, FILE *stream);

GMX_LIBGMX_EXPORT
void strip_comment (char *line);

int break_line (char *line,
                char *variable,
                char *value);

GMX_LIBGMX_EXPORT
void upstring (char *str);

GMX_LIBGMX_EXPORT
void ltrim (char *str);

GMX_LIBGMX_EXPORT
void rtrim (char *str);

GMX_LIBGMX_EXPORT
void trim (char *str);

GMX_LIBGMX_EXPORT
void nice_header (FILE *out, const char *fn);

GMX_LIBGMX_EXPORT
int gmx_strcasecmp_min(const char *str1, const char *str2);
GMX_LIBGMX_EXPORT
int gmx_strncasecmp_min(const char *str1, const char *str2, int n);
/* This funny version of strcasecmp, is not only case-insensitive,
 * but also ignores '-' and '_'.
 */

GMX_LIBGMX_EXPORT
int gmx_strcasecmp(const char *str1, const char *str2);
GMX_LIBGMX_EXPORT
int gmx_strncasecmp(const char *str1, const char *str2, int n);

GMX_LIBGMX_EXPORT
char *gmx_strdup(const char *src);
char *gmx_strndup(const char *src, int n);

/* Magic hash initialization number from Dan J. Bernstein. */
extern const unsigned int
    gmx_string_hash_init;

/* Return a hash of the string according to Dan J. Bernsteins algorithm.
 * On the first invocation for a new string, use the constant
 * gmx_string_hash_init for the second argument. If you want to create a hash
 * corresponding to several concatenated strings, provide the returned hash
 * value as hash_init for the second string, etc.
 */
unsigned int
gmx_string_fullhash_func(const char *s, unsigned int hash_init);

/* Identical to gmx_string_fullhash_func, except that
 * this routine only uses characters for which isalnum(c) is true,
 * and all characters are converted to upper case.
 */
unsigned int
gmx_string_hash_func(const char *s, unsigned int hash_init);

/** Pattern matcing with wildcards. */
GMX_LIBGMX_EXPORT
int gmx_wcmatch(const char *pattern, const char *src);

/** Return value for gmx_wcmatch() when there is no match. */
#define GMX_NO_WCMATCH 1


/* this is our implementation of strsep, the thread-safe replacement for
   strtok */
GMX_LIBGMX_EXPORT
char *gmx_strsep(char **stringp, const char *delim);


GMX_LIBGMX_EXPORT
char *wrap_lines(const char *buf, int line_width, int indent,
                 gmx_bool bIndentFirst);
/* wraps lines at 'linewidth', indenting all following
 * lines by 'indent' spaces. A temp buffer is allocated and returned,
 * which can be disposed of if no longer needed.
 * If !bIndentFirst, then the first line will not be indented, only
 * the lines that are created due to wapping.
 */


char **split(char sep, char *str);
/* Implementation of the well-known Perl function split */

gmx_large_int_t str_to_large_int_t(const char *str, char **endptr);

#ifdef GMX_NATIVE_WINDOWS
#define snprintf _snprintf
#endif

#ifdef __cplusplus
}
#endif

#endif  /* _string2_h */
