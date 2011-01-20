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

#ifndef _futil_h
#define _futil_h

#include <stdio.h>
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
}
#endif

/* Native windows uses backslash path separators.
 * Cygwin and everybody else in the world use slash.
 * When reading the PATH environment variable, Unix separates entries
 * with colon, while windows uses semicolon.
 */
#if ((defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64) && !defined __CYGWIN__ && !defined __CYGWIN32__)
#define DIR_SEPARATOR '\\'
#define PATH_SEPARATOR ";"
#else
#define DIR_SEPARATOR '/'
#define PATH_SEPARATOR ":"
#endif


/* Now get the maximum path size. */
#ifdef PATH_MAX
#  define GMX_PATH_MAX PATH_MAX
#elif defined MAX_PATH
#  define GMX_PATH_MAX MAX_PATH
#else
#  define GMX_PATH_MAX 4096
#endif


#ifdef HAVE_FSEEKO
   typedef off_t              gmx_off_t;
#  define SIZEOF_GMX_OFF_T   SIZEOF_OFF_T
#elif defined HAVE__FSEEKI64 
   typedef __int64            gmx_off_t;
#  define SIZEOF_GMX_OFF_T   8
#else
   /* Almost certainly 64 bits, and guaranteed to be available */
   typedef gmx_large_int_t    gmx_off_t;
#  define SIZEOF_GMX_OFF_T   SIZEOF_GMX_LARGE_INT
#endif    


  
void no_buffers(void);
/* Turn off buffering of files (which is default) for debugging purposes */

gmx_bool gmx_fexist(const char *fname);
/* Return TRUE when fname exists, FALSE otherwise */

gmx_bool gmx_fexist_master(const char *fname, t_commrec *cr);
/* Return TRUE when fname exists, FALSE otherwise, bcast from master to others */

gmx_bool gmx_eof(FILE *fp);
/* Return TRUE on end-of-file, FALSE otherwise */

gmx_bool is_pipe(FILE *fp);
/* Check whether the file (opened by ffopen) is a pipe */

/*  Make a backup of file if necessary.  
    Return false if there was a problem.
*/
gmx_bool make_backup(const char * file);

FILE *ffopen(const char *file, const char *mode);
/* Return a valid file pointer when successful, exits otherwise 
 * If the file is in compressed format, open a pipe which uncompresses
 * the file! Therefore, files must be closed with ffclose (see below)
 */

int ffclose(FILE *fp);
/* Close files or pipes */


void frewind(FILE *fp);
/* Does not rewind pipes, but does so for normal files */

#define rewind frewind


int gmx_fseek(FILE *stream, gmx_off_t offset, int whence); 
/* OS-independent fseek. 64-bit when available */

gmx_off_t gmx_ftell(FILE *stream); 
/* OS-independent fseek. 64-bit when available. */


gmx_bool is_pipe(FILE *fp);

char *gmxlibfn(const char *file);
/* allocates and returns a string with the full file name for a library file */

FILE *libopen(const char *file);
/* Open a library file for reading. This looks in the current directory
 * first, and then in the library directory. If the file is not found,
 * it terminates with a fatal_error
 */

/* Opaque data type to list directories */
typedef struct gmx_directory *
gmx_directory_t;

/* Open a directory for reading. The first argument should be a pointer
 * to a declared gmx_directory_t variable. Returns 0 on success.
 */
int
gmx_directory_open(gmx_directory_t *p_gmxdir,const char *dirname);

    
/* Given an initialized gmx_directory_t, if there are more files in
 * the directory this routine returns 0 and write the next name
 * into the USER-PROVIDED buffer name. The last argument is the max
 * number of characters that will be written. Just as strncpy, the
 * string will NOT be terminated it it is longer than maxlength_name.
 */
int
gmx_directory_nextfile(gmx_directory_t gmxdir,char *name,int maxlength_name);
    
/* Release all data for a directory structure */
int 
gmx_directory_close(gmx_directory_t gmxdir);
    

    
gmx_bool get_libdir(char *libdir);

char *low_gmxlibfn(const char *file,gmx_bool bAddCWD,gmx_bool bFatal);

FILE *low_libopen(const char *file,gmx_bool bFatal);
/* The same as the above, but does not terminate if (!bFatal) */

/* Create unique name for temp file (wrapper around mkstemp). 
 * Buf should be at least 7 bytes long 
 */
void gmx_tmpnam(char *buf);

/* truncte the file to the specified length */
int gmx_truncatefile(char *path, gmx_off_t length);

/* rename/move the file (atomically, if the OS makes that available) oldname 
   to newname */
int gmx_file_rename(const char *oldname, const char *newname);

/* copy the file (data only) oldname to newname. if copy_if_empty==FALSE,
   the file won't be copied if it's empty.*/
int gmx_file_copy(const char *oldname, const char *newname, gmx_bool copy_if_empty);

/* do an fsync() on an open file pointer. 
   Only use this during checkpointing! */
int gmx_fsync(FILE *fp);

#ifdef __cplusplus
}
#endif

#endif	/* _futil_h */
