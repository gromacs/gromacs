/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2012,2013,2014, by the GROMACS development team, led by
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
/*! \file
 * \brief Generic C string handling functions.
 *
 * \inpublicapi
 * \ingroup module_utility
 */
#ifndef GMX_UTILITY_CSTRINGUTIL_H
#define GMX_UTILITY_CSTRINGUTIL_H

#include <stdio.h>
#include <time.h>

#include "../legacyheaders/types/simple.h"

#include "gmx_header_config.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/** Continuation character. */
#define CONTINUE    '\\'
/** Comment sign to use. */
#define COMMENTSIGN ';'

/*! \brief
 * Strip trailing spaces and if s ends with a ::CONTINUE remove that too.
 *
 * \returns TRUE if s ends with a CONTINUE, FALSE otherwise.
 */
int continuing(char *s);

/*! \brief
 * Reads a line from a stream.
 *
 * This routine reads a string from stream of max length n
 * and zero terminated, without newlines.
 * \p s should be long enough (>= \p n)
 */
char *fgets2(char *s, int n, FILE *stream);

/** Remove portion of a line after a ::COMMENTSIGN.  */
void strip_comment(char *line);

/** Make a string uppercase. */
void upstring(char *str);

/** Remove leading whitespace from a string. */
void ltrim(char *str);

/** Remove trailing whitespace from a string. */
void rtrim(char *str);

/** Remove leading and trailing whitespace from a string. */
void trim(char *str);

/** Portable version of ctime_r. */
char *gmx_ctime_r(const time_t *clock, char *buf, int n);

/** Prints creation time stamp and user information into a file as comments. */
void nice_header(FILE *out, const char *fn);

/** Version of gmx_strcasecmp() that also ignores '-' and '_'. */
int gmx_strcasecmp_min(const char *str1, const char *str2);
/** Version of gmx_strncasecmp() that also ignores '-' and '_'. */
int gmx_strncasecmp_min(const char *str1, const char *str2, int n);

/** Case-insensitive strcmp(). */
int gmx_strcasecmp(const char *str1, const char *str2);
/** Case-insensitive strncmp(). */
int gmx_strncasecmp(const char *str1, const char *str2, int n);

/** Creates a duplicate of \p src. */
char *gmx_strdup(const char *src);
/** Duplicates first \p n characters of \p src. */
char *gmx_strndup(const char *src, int n);

/*! \brief
 * Pattern matching with wildcards.
 *
 * \param[in] pattern  Pattern to match against.
 * \param[in] str      String to match.
 * \returns   0 on match, GMX_NO_WCMATCH if there is no match.
 *
 * Matches \p str against \p pattern, which may contain * and ? wildcards.
 * All other characters are matched literally.
 * Currently, it is not possible to match literal * or ?.
 */
int gmx_wcmatch(const char *pattern, const char *str);

/** Return value for gmx_wcmatch() when there is no match. */
#define GMX_NO_WCMATCH 1

/** Magic hash initialization number from Dan J. Bernstein. */
extern const unsigned int
    gmx_string_hash_init;

/*! \brief
 * Return a hash of the string according to Dan J. Bernsteins algorithm.
 *
 * \param[in] s          String to calculate hash for.
 * \param[in] hash_init  Initial (or previous) hash value.
 * \returns   Updated hash value (hash_init combined with string hash).
 *
 * On the first invocation for a new string, use the constant
 * gmx_string_hash_init for the second argument. If you want to create a hash
 * corresponding to several concatenated strings, provide the returned hash
 * value as hash_init for the second string, etc.
 */
unsigned int
gmx_string_fullhash_func(const char *s, unsigned int hash_init);

/*! \brief
 * Return a hash of the string according to Dan J. Bernsteins algorithm.
 *
 * \param[in] s          String to calculate hash for.
 * \param[in] hash_init  Initial (or previous) hash value.
 * \returns   Updated hash value (hash_init combined with string hash).
 *
 * Identical to gmx_string_fullhash_func, except that
 * this routine only uses characters for which isalnum(c) is true,
 * and all characters are converted to upper case.
 */
unsigned int
gmx_string_hash_func(const char *s, unsigned int hash_init);

/*! \brief
 * Wraps lines, optionally indenting lines.
 *
 * Wraps lines at \p linewidth, indenting all following lines by \p indent
 * spaces.  A temp buffer is allocated and returned, which can be disposed of
 * if no longer needed.
 * If \p bIndentFirst is FALSE, then the first line will not be indented, only
 * the lines that are created due to wapping.
 */
char *wrap_lines(const char *buf, int line_width, int indent,
                 gmx_bool bIndentFirst);

/*! \brief
 * Convert a string to gmx_int64_t.
 *
 * This method works as the standard library function strtol(), except that it
 * does not support different bases.
 */
gmx_int64_t str_to_int64_t(const char *str, char **endptr);

#ifdef GMX_NATIVE_WINDOWS
#define snprintf _snprintf
#endif

/*! \brief Construct an array of digits found in the input string
 *
 * \param[in]  digitstring  String that must contain only digits
 * \param[out] ndigits      Size of return array with the values of the digits
 * \param[out] digitlist    Array of digits found in
 *                          digitstring. Allocated by this function
 *                          with size *ndigits. Calling code is
 *                          responsible for deallocation.
 *
 * If digitstring is NULL, then ndigits is set to zero and digitlist
 * to NULL. If digitstring contains a non-digit character, a fatal
 * error results.
 */
void parse_digits_from_plain_string(const char *digitstring, int *ndigits, int **digitlist);

#ifdef __cplusplus
}
#endif

#endif
