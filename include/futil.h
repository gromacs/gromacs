/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef	_futil_h
#define	_futil_h

#ifdef HAVE_IDENT
#ident	"@(#) futil.h 1.1 11/23/92"
#endif /* HAVE_IDENT */
#include <stdio.h>
#include "typedefs.h"

#ifdef CPLUSPLUS
extern "C" { 
#endif

extern bool fexist(char *fname);
/* Return TRUE when fname exists, FALSE otherwise */

extern bool eof(FILE *fp);
/* Return TRUE on end-of-file, FALSE otherwise */

extern char *backup_fn(char *file);
/* Return a backup name for file (name with # before and after) */

extern FILE *ffopen(char *file,char *mode);
/* Return a valid file pointer when succesfull, exits otherwise 
 * If the file is in compressed format, open a pipe which uncompresses
 * the file! Therefore, files must be closed with ffclose (see below)
 */

extern void ffclose(FILE *fp);
/* Close files or pipes */

#define fclose ffclose

extern void frewind(FILE *fp);
/* Does not rewind pipes, but does so for normal files */

#define rewind frewind

extern FILE *uncompress(char *fn);
extern FILE *gunzip(char *fn);
/* Open a pipe to uncompress or unzip files. Must be closed with pclose */

extern char *libfn(char *file);

extern FILE *libopen(char *file);
/* Open a library file for reading. This looks in the current directory
 * first, and then in the library directory. If the file is not found,
 * it terminates with a fatal_error
 */

extern char *low_libfn(char *file,bool bFatal);

extern FILE *low_libopen(char *file,bool bFatal);
/* The same as the above, but does not terminate if (!bFatal) */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _futil_h */
