/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef _string2_h
#define _string2_h

static char *SRCID_string2_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) string2.h 1.13 11/23/92"
#endif /* HAVE_IDENT */

/*
 *
 * string2.h
 * David van der Spoel
 *
 */

#ifdef CPLUSPLUS
extern "C" {
#endif

#include <string.h>
#include <stdio.h>

#define CONTINUE    '\\'
#define COMMENTSIGN ';'
#define STRLEN      4096

int continuing(char *s);

char *fgets2(char *s, int n, FILE *stream);

void strip_comment (char *line);

int break_line (char *line,
		char *variable,
		char *value);

int upcase (int ch);

void ltrim (char *str);

void rtrim (char *str);

void trim (char *str);

void nice_header (FILE *out,char *fn);

int gmx_strcasecmp(const char *str1, const char *str2);
/* This funny version of strcasecmp, is not only case-insensitive,
 * but does also return 0 when '-' and '_' are compared.
 */

char *gmx_strdup(const char *src);

#ifndef HAVE_STRCASECMP
#define strcasecmp gmx_strcasecmp
#endif

#ifndef HAVE_STRDUP
#define strdup gmx_strdup
#endif

#ifdef CPLUSPLUS
}
#endif

#endif	/* _string2_h */
