/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _fatal_h
#define _fatal_h

static char *SRCID_fatal_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) fatal.h 1.9 11/23/92"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" {
#endif
  
#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
  
extern void _where(char *file,int line);
#define where() _where(__FILE__,__LINE__)
/* Prints filename and line to stdlog and only on amba memvail */
  
extern void _halt(char *file,int line,char *reason);
#define HALT(reason) _halt(__FILE__,__LINE__,reason)
/* Halts the program with an error message */

extern void _set_fatal_tmp_file(char *fn, char *file, int line);
#define set_fatal_tmp_file(fn) _set_fatal_tmp_file(fn,__FILE__,__LINE__)
/* set filename to be removed when fatal_error is called */

extern void _unset_fatal_tmp_file(char *fn, char *file, int line);
#define unset_fatal_tmp_file(fn) _unset_fatal_tmp_file(fn,__FILE__,__LINE__)
/* unsets filename to be removed */

extern void fatal_error(int fatal_errno,char *fmt,...);
/*
 * Routine fatal_error prints 
 *
 * 	"fatal error file %s line %s \n\t " 
 *
 * followed by the string specified by fmt and supplied parameters. If 
 * errno is 0, only the message and arguments are printed. If errno is 
 * a legal system errno or -1, a perror like message is printed after the
 * first message, if errno is -1, the last system errno will be used.
 * The format of fmt is that like printf etc, only %d, %x, %c, %f and %s
 * are allowed as format specifiers.
 */

/* This include must not be moved upwards, to prevent compilation problems */  
#include "typedefs.h"

extern void init_warning(int maxwarning);
/* Set the max number of warnings */

extern void set_warning_line(char *fn,int line);
/* Set filename and linenumber for the warning */
  
extern char warn_buf[1024];
/* Warning buffer of 1024 bytes, which can be used to print messages to */

extern void warning(char *s);
/* Issue a warning, with the string s. If s == NULL, then warn_buf
 * will be printed instead.
 */
 
extern void print_warn_num(void);
/* Print the total number of warnings, if larger than 0 */
  
extern void _too_few(char *fn,int line);
#define too_few() _too_few(__FILE__,__LINE__)
/* Issue a warning stating 'Too few parameters' */
  
extern void _invalid_case(char *fn,int line);
#define invalid_case() _invalid_case(__FILE__,__LINE__)
/* Issue a warning stating 'Invalid case in switch' */
  
extern void _unexpected_eof(char *fn,int line,char *srcfn,int srcline);
#define unexpected_eof(fn,line) _unexpected_eof(fn,line,__FILE__,__LINE__)
  
/* 
 * Functions can write to this file for debug info
 * Before writing to it, it should be checked whether
 * the file is not 0:
 * if (debug) fprintf(debug,"%s","Hallo");
 */
extern FILE *debug;
  
void init_debug (char *dbgfile);
  
extern bool bDebugMode(void);
/* Return TRUE when the program was started in debug mode */
  
#ifdef USE_SGI_FPE
extern void doexceptions(void);
/* Set exception handlers for debugging */
#endif

extern void check_nprocs_top(char *fn,t_topology *top,int nprocs);
/* Verify whether this tpr file is for nprocs processors, and quit if not */
  
#ifdef CPLUSPLUS
	   }
#endif

#endif	/* _fatal_h */
