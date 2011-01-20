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

#ifndef _fflibutil_h
#define _fflibutil_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

extern const char *fflib_forcefield_dir_ext();
/* Returns the name of the force field directory extension */

extern const char *fflib_forcefield_itp();
/* Returns the name of the main forcefield itp file */

extern const char *fflib_forcefield_doc();
/* Returns the name of the forcefield documentation file */

extern void fflib_filename_base(const char *filename,char *filebase,int maxlen);
/* Return the base file name of filename in base,
 * i.e. remove path and extension, if present.
 * base should be at least of size maxlen.
 */

extern int fflib_search_file_end(const char *ffdir,
				 const char *file_end,
				 gmx_bool bFatalError,
				 char ***filenames);
/* Search for files ending on file_end in the force field directory fflib.
 * fflib should be in the GROMACS lib.path.
 * Return the number of files and the file names in filenames.
 */

extern int fflib_search_file_in_dirend(const char *filename,const char *dirend,
				       char ***dirnames);
/* Search for files with name filename in subdirectories with names
 * ending on dirend.
 * Return the number of files and the directory names in dirnames.
 */
extern gmx_bool fflib_fexist(const char *file);
/* Check if a file exists in the force field library */

extern FILE *fflib_open(const char *file);
/* Open force field library file "file" for reading.
 * "file" should contain the whole path to the force field library,
 * either absolute or relative to the current dir.
 */

#ifdef __cplusplus
}
#endif

#endif	/* _fflibutil_h */
