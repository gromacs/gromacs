/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Grunge ROck MAChoS
 */

#ifndef _futil_h
#define _futil_h

static char *SRCID_futil_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) futil.h 1.1 11/23/92"
#endif /* HAVE_IDENT */
#include <stdio.h>
#include "typedefs.h"

#ifdef CPLUSPLUS
extern "C" { 
#endif
  
extern void no_buffers(void);
/* Turn off buffering of files (which is default) for debugging purposes */

extern bool fexist(char *fname);
/* Return TRUE when fname exists, FALSE otherwise */

extern bool eof(FILE *fp);
/* Return TRUE on end-of-file, FALSE otherwise */

extern bool is_pipe(FILE *fp);
/* Check whether the file (opened by ffopen) is a pipe */

extern char *backup_fn(char *file);
/* Return a backup name for file (name with # before and after) */

/*  Make a backup of file if necessary.  
    Return false if there was a problem.
*/
extern bool make_backup(char * file);

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

bool is_pipe(FILE *fp);

extern FILE *uncompress(char *fn,char *mode);
extern FILE *gunzip(char *fn,char *mode);
/* Open a pipe to uncompress or unzip files. Must be closed with pclose */

extern char *libfn(char *file);

  extern FILE *libopen(char *file);
/* Open a library file for reading. This looks in the current directory
 * first, and then in the library directory. If the file is not found,
 * it terminates with a fatal_error
 */
  
extern bool get_libdir(char *libdir);

extern char *low_libfn(char *file,bool bFatal);

extern FILE *low_libopen(char *file,bool bFatal);
/* The same as the above, but does not terminate if (!bFatal) */

/* Create unique name for temp file (wrapper around mkstemp). 
 * Buf should be at least 7 bytes long 
 */
extern void gmx_tmpnam(char *buf);

#ifdef CPLUSPLUS
}
#endif

#endif	/* _futil_h */
