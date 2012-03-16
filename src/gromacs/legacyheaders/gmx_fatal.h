
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

#ifndef _fatal_h
#define _fatal_h

#include <stdio.h>
#include <stdarg.h>
#include <errno.h>
#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif
  
void 
_where(const char *file,int line);
#define where() _where(__FILE__,__LINE__)
/* Prints filename and line to stdlog and only on amba memvail */
  
void 
_set_fatal_tmp_file(const char *fn, const char *file, int line);
#define set_fatal_tmp_file(fn) _set_fatal_tmp_file(fn,__FILE__,__LINE__)
/* set filename to be removed when fatal_error is called */

void 
_unset_fatal_tmp_file(const char *fn, const char *file, int line);
#define unset_fatal_tmp_file(fn) _unset_fatal_tmp_file(fn,__FILE__,__LINE__)
/* unsets filename to be removed */

void 
gmx_fatal(int fatal_errno,const char *file,int line,const char *fmt,...);
#define FARGS 0,__FILE__,__LINE__
/*
 * Routine gmx_fatal prints 
 *
 * 	"fatal error file %s line %s \n\t " 
 *
 * followed by the string specified by fmt and supplied parameters. If 
 * errno is 0, only the message and arguments are printed. If errno is 
 * a legal system errno or -1, a perror like message is printed after the
 * first message, if errno is -1, the last system errno will be used.
 * The format of fmt is that like printf etc, only %d, %x, %c, %f, %g and %s
 * are allowed as format specifiers.
 *
 * Tip of the week:
 * call this function using the FARGS macro:
 * gmx_fatal(FARGS,fmt,...)
 */

void
gmx_fatal_collective(int f_errno,const char *file,int line,
		     t_commrec *cr,gmx_domdec_t *dd,
		     const char *fmt,...);
/* As gmx_fatal, but only the master process prints the error message.
 * This should only be called one of the following two situations:
 * 1) On all nodes in cr->mpi_comm_mysim, with cr!=NULL,dd==NULL.
 * 2) On all nodes in dd->mpi_comm_all,   with cr==NULL,dd!=NULL.
 * This will call MPI_Finalize instead of MPI_Abort when possible,
 * This is useful for handling errors in code that is executed identically
 * for all processes.
 */

void
gmx_fatal_set_log_file(FILE *fp);
/* Set the log file for printing error messages */
  
void 
_invalid_case(const char *fn,int line);
#define invalid_case() _invalid_case(__FILE__,__LINE__)
/* Issue a warning stating 'Invalid case in switch' */
  
void _unexpected_eof(const char *fn,int line,const char *srcfn,int srcline);
#define unexpected_eof(fn,line) _unexpected_eof(fn,line,__FILE__,__LINE__)

/* 
 * Functions can write to this file for debug info
 * Before writing to it, it should be checked whether
 * the file is not NULL:
 * if (debug) fprintf(debug,"%s","Hallo");
 */
extern FILE *debug;
extern gmx_bool gmx_debug_at;

void init_debug (const int dbglevel,const char *dbgfile);
  
gmx_bool bDebugMode(void);
/* Return TRUE when the program was started in debug mode */
  
#if (defined __sgi && defined USE_SGI_FPE)
void doexceptions(void);
/* Set exception handlers for debugging */
#endif

  /* warn_str is allowed to be NULL.
   */
  void _range_check(int n,int n_min,int n_max,const char *warn_str,
			   const char *var,
			   const char *file,int line);

#define range_check_mesg(n,n_min,n_max,str) _range_check(n,n_min,n_max,str,#n,__FILE__,__LINE__)
  /* Range check will terminate with an error message if not
   * n E [ n_min, n_max >
   * That is n_min is inclusive but not n_max.
   */

#define range_check(n,n_min,n_max) _range_check(n,n_min,n_max,NULL,#n,__FILE__,__LINE__)
  /* Range check will terminate with an error message if not
   * n E [ n_min, n_max >
   * That is n_min is inclusive but not n_max.
   */

  char *gmx_strerror(const char *key);
  /* Return error message corresponding to the key.
   * Maybe a multi-line message.
   * The messages are stored in src/gmxlib/fatal.c
   */
  
  void _gmx_error(const char *key,const char *msg,const char *file,int line);
#define gmx_error(key,msg) _gmx_error(key,msg,__FILE__,__LINE__)
  /* Error msg of type key is generated and the program is 
   * terminated unless and error handle is set (see below)
   */

  /* Some common error types */
#define gmx_bug(msg)    gmx_error("bug",msg)
#define gmx_call(msg)   gmx_error("call",msg)
#define gmx_comm(msg)   gmx_error("comm",msg)
#define gmx_file(msg)   gmx_error("file",msg)
#define gmx_cmd(msg)    gmx_error("cmd",msg)
#define gmx_impl(msg)   gmx_error("impl",msg)
#define gmx_incons(msg) gmx_error("incons",msg)
#define gmx_input(msg)  gmx_error("input",msg)
#define gmx_mem(msg)    gmx_error("mem",msg)
#define gmx_open(fn)    gmx_error("open",fn) 
  
void 
set_gmx_error_handler(void (*func)(const char *msg));
/* An error function will be called that terminates the program 
   * with a fatal error, unless you override it with another function.
   * i.e.:
   * set_gmx_error_handler(my_func);
   * where my_func is a function that takes a string as an argument.
   * The string may be a multi-line string.
   */

void gmx_warning(const char *fmt,...);
/* Print a warning message to stderr.
 * The format of fmt is that like printf etc, only %d, %x, %c, %f, %g and %s
 * are allowed as format specifiers.
 * The message string should NOT start with "WARNING"
 * and should NOT end with a newline.
 */

#ifdef __cplusplus
	   }
#endif

#endif	/* _fatal_h */
