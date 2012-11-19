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
#include "../utility/gmx_header_config.h"

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
 * This routine only uses characters for which isalnum(c) is true,
 * and all characters are converted to upper case.
 * On the first invocation for a new string, use the constant
 * gmx_string_hash_init for the second argument. If you want to create a hash
 * corresponding to several concatenated strings, provide the returned hash
 * value as hash_init for the second string, etc.
 */
unsigned int
gmx_string_hash_func(const char *s, unsigned int hash_init);
    
/** Return value for gmx_wcmatch() when there is no match. */
#define GMX_NO_WCMATCH 1

/** Our implementation of strsep, the thread-safe replacement for strtok. */
char *gmx_strsep(char **stringp, const char *delim);

/*! \brief
 * Wraps lines, optionally indenting lines.
 *
 * Wraps lines at \p linewidth, indenting all following lines by \p indent
 * spaces.  A temp buffer is allocated and returned, which can be disposed of
 * if no longer needed.
 * If \p bIndentFirst is FALSE, then the first line will not be indented, only
 * the lines that are created due to wapping.
 */
char *wrap_lines(const char *buf,int line_width, int indent,
                 gmx_bool bIndentFirst);

/** Implementation of the well-known Perl function split. */
char **split(char sep,const char *str);

/** return new pointer to a string with f printed. Pointer can be discarded after use. */ 
extern char *gmx_ftoa(double f);
 
/** return new pointer to a string with f printed. Pointer can be discarded after use. */  
extern char *gmx_itoa(int f);

/*! \brief
 * Convert a string to gmx_large_int_t.
 *
 * This method works as the standard library function strtol(), except that it
 * does not support different bases.
 *
 * \attention
 * The following differences are present from the standard behavior:
 *  - \p endptr cannot be NULL.
 *  - If an overflow occurs, returns zero and \p *endptr will equal \p str.
 *    errno is still set to ERANGE.
 */
gmx_large_int_t str_to_large_int_t(const char *str, char **endptr);
    
#ifdef GMX_NATIVE_WINDOWS
#define snprintf _snprintf
#endif

#ifdef __cplusplus
}
#endif

#endif	/* _string2_h */
