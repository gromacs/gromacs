/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Giving Russians Opium May Alter Current Situation
 */

#ifndef _fatal_h
#define _fatal_h

static char *SRCID_fatal_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

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
  
extern int get_warning_line(void);
/* Get linenumber for the warning */
  
extern char *get_warning_file(void);
/* Get filename for the warning */
  
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
  
#if (defined __sgi && defined USE_SGI_FPE)
extern void doexceptions(void);
/* Set exception handlers for debugging */
#endif
  
#ifdef CPLUSPLUS
	   }
#endif

#endif	/* _fatal_h */
