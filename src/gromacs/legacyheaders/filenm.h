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

#ifndef _filenm_h
#define _filenm_h

#include "futil.h"

#ifdef __cplusplus
extern "C" {
#endif
  
void set_default_file_name(const char *name);
/* Set the default file name for all file types to name */

const char *ftp2ext(int ftp);
/* Return extension for filetype */

const char *ftp2ext_generic(int ftp);
/* Return extension for filetype, and a generic name for generic types 
   (e.g. trx)*/

const char *ftp2desc(int ftp);
/* Return description for file type */

const char *ftp2defnm(int ftp);
/* Return default file name for file type */

const char *ftp2ftype(int ftp);
/* Return Binary or ASCII depending on file type */

void pr_def(FILE *fp,int ftp);
/* Print definitions for filename ftp */

void pr_defs(FILE *fp);
/* Print definitions for all filename */

void pr_fns(FILE *fp,int nf,const t_filenm tfn[]);
/* Print nf file names and types */

void pr_fopts(FILE *fp,int nf,const t_filenm tfn[], int shell);
/* prints file options in tcsh 'complete' format */

void parse_file_args(int *argc,char *argv[],int nf,t_filenm fnm[],
			    gmx_bool bKeep, gmx_bool bReadNode);
/* Parse command line for file names. When bKeep is set args are 
 * not removed from argv. */

const char *opt2fn(const char *opt,int nfile, const t_filenm fnm[]);
/* Return the filename belonging to cmd-line option opt, or NULL when 
 * no such option. */

const char *opt2fn_master(const char *opt, int nfile, 
                                const t_filenm fnm[], t_commrec *cr);
/* Return the filename belonging to cmd-line option opt, or NULL when 
 * no such option or not running on master */


int opt2fns(char **fns[], const char *opt,int nfile,
                   const t_filenm fnm[]);
/* Return the filenames belonging to cmd-line option opt, or NULL when 
 * no such option. */

#define opt2FILE(opt,nfile,fnm,mode) ffopen(opt2fn(opt,nfile,fnm),mode)
/* Return a file pointer from the filename (see above) */

int fn2ftp(const char *fn);
/* Return the filetype corrsponding to filename */

const char *ftp2fn(int ftp,int nfile,const t_filenm fnm[]);
/* Return the first file name with type ftp, or NULL when none found. */

int ftp2fns(char **fns[], int ftp,int nfile,const t_filenm fnm[]);
/* Return the number of files for the first option with type ftp
   and the files in **fns[] (will be allocated), or NULL when none found. */
 
#if 0
/* This function is not thread-safe and used nowhere: */
char *ftp2filter(int ftp);
/* Return a file extension filter for file type */
#endif

#define ftp2FILE(ftp,nfile,fnm,mode) ffopen(ftp2fn(ftp,nfile,fnm),mode)
/* Return a file pointer from the filename (see above) */

gmx_bool ftp2bSet(int ftp,int nfile,const t_filenm fnm[]);
/* Return TRUE when this file type has been found on the cmd-line */

gmx_bool opt2bSet(const char *opt,int nfile,const t_filenm fnm[]);
/* Return TRUE when this option has been found on the cmd-line */

const char *opt2fn_null(const char *opt,int nfile,const t_filenm fnm[]);
/* Return the filenm belonging top cmd-line option opt, or NULL when 
 * no such option. 
 * Also return NULL when opt is optional and option is not set. 
 */

const char *ftp2fn_null(int ftp,int nfile,const t_filenm fnm[]);
/* Return the first file name with type ftp, or NULL when none found.
 * Also return NULL when ftp is optional and option is not set.
 */

gmx_bool is_optional(const t_filenm *fnm);
/* Return whether or not this filenm is optional */

gmx_bool is_output(const t_filenm *fnm);
/* Return whether or not this filenm is output */

gmx_bool is_set(const t_filenm *fnm);
/* Return whether or not this filenm is set */

/* When we do checkpointing, this routine is called to check for previous
 * output files and append a '.partNNNN' suffix before the (output) file extensions.
 */
int add_suffix_to_output_names(t_filenm *fnm, int nfile, const char *suffix);

/* duplicate the filename list (to make a private copy for each thread, 
   for example) */
t_filenm *dup_tfn(int nf, const t_filenm tfn[]);

/* Free memory allocated for file names by parse_file_args(). */
void done_filenms(int nf, t_filenm fnm[]);
	
#ifdef __cplusplus
}
#endif

#endif	/* _filenm_h */
