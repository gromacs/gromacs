/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _filenm_h
#define _filenm_h

static char *SRCID_filenm_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) filenm.h 1.14 2/2/97"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif

#include "futil.h"

void set_default_file_name(char *name);
/* Set the default file name for all file types to name */

extern char *ftp2ext(int ftp);
/* Return extension for filetype */

extern char *ftp2desc(int ftp);
/* Return description for file type */

extern char *ftp2defnm(int ftp);
/* Return default file name for file type */

extern char *ftp2ftype(int ftp);
/* Return Binary or ASCII depending on file type */

extern void pr_def(FILE *fp,int ftp);
/* Print definitions for filename ftp */

extern void pr_defs(FILE *fp);
/* Print definitions for all filename */

extern void pr_fns(FILE *fp,int nf,t_filenm tfn[]);
/* Print nf file names and types */

extern void pr_fopts(FILE *fp,int nf,t_filenm tfn[]);
/* prints file options in tcsh 'complete' format */

extern void parse_file_args(int *argc,char *argv[],int nf,t_filenm fnm[],
			    bool bKeep);
/* Parse command line for file names. When bKeep is set args are 
 * not removed from argv.
 */

extern char *opt2fn(char *opt,int nfile,t_filenm fnm[]);
/* Return the filenm belonging top cmd-line option opt, or NULL when 
 * no such option. 
 */

#define opt2FILE(opt,nfile,fnm,mode) ffopen(opt2fn(opt,nfile,fnm),mode)
/* Return a file pointer from the filename (see above) */

extern int fn2ftp(char *fn);
/* Return the filetype corrsponding to filename */

extern char *ftp2fn(int ftp,int nfile,t_filenm fnm[]);
/* Return the first file name with type ftp, or NULL when none found. */

extern char *ftp2filter(int ftp);
/* Return a file extension filter for file type */

#define ftp2FILE(ftp,nfile,fnm,mode) ffopen(ftp2fn(ftp,nfile,fnm),mode)
/* Return a file pointer from the filename (see above) */

extern bool ftp2bSet(int ftp,int nfile,t_filenm fnm[]);
/* Return TRUE when this file type has been found on the cmd-line */

extern bool opt2bSet(char *opt,int nfile,t_filenm fnm[]);
/* Return TRUE when this option has been found on the cmd-line */

extern char *opt2fn_null(char *opt,int nfile,t_filenm fnm[]);
/* Return the filenm belonging top cmd-line option opt, or NULL when 
 * no such option. 
 * Also return NULL when opt is optional and option is not set. 
 */

extern char *ftp2fn_null(int ftp,int nfile,t_filenm fnm[]);
/* Return the first file name with type ftp, or NULL when none found.
 * Also return NULL when ftp is optional and option is not set.
 */

extern bool is_optional(t_filenm *fnm);
/* Return whether or not this filenm is optional */

#ifdef CPLUSPLUS
}
#endif

#endif	/* _filenm_h */
