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

#define CONTINUE    '\\'
#define COMMENTSIGN ';'

int continuing(char *s);

char *fgets2(char *s, int n, FILE *stream);

void strip_comment (char *line);

int break_line (char *line,
		       char *variable,
		       char *value);

void upstring (char *str);

void ltrim (char *str);

void rtrim (char *str);

void trim (char *str);

void nice_header (FILE *out,const char *fn);

int gmx_strcasecmp_min(const char *str1, const char *str2);
int gmx_strncasecmp_min(const char *str1, const char *str2, int n);
/* This funny version of strcasecmp, is not only case-insensitive,
 * but also ignores '-' and '_'.
 */

int gmx_strcasecmp(const char *str1, const char *str2);
int gmx_strncasecmp(const char *str1, const char *str2, int n);

char *gmx_strdup(const char *src);
char *gmx_strndup(const char *src, int n);
    
/** Pattern matcing with wildcards. */
int gmx_wcmatch(const char *pattern, const char *src);

/** Return value for gmx_wcmatch() when there is no match. */
#define GMX_NO_WCMATCH 1


/* this is our implementation of strsep, the thread-safe replacement for
   strtok */
char *gmx_strsep(char **stringp, const char *delim);


char *wrap_lines(const char *buf,int line_width, int indent,
			gmx_bool bIndentFirst);
/* wraps lines at 'linewidth', indenting all following
 * lines by 'indent' spaces. A temp buffer is allocated and returned,
 * which can be disposed of if no longer needed.
 * If !bIndentFirst, then the first line will not be indented, only 
 * the lines that are created due to wapping.
 */


char **split(char sep,char *str);
/* Implementation of the well-known Perl function split */

gmx_large_int_t str_to_large_int_t(const char *str, char **endptr);

#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#define snprintf _snprintf
#endif

#ifdef __cplusplus
}
#endif

#endif	/* _string2_h */
