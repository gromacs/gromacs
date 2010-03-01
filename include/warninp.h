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

#ifndef _warninp_h
#define _warninp_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

void 
init_warning(int maxwarning);
/* Set the max number of warnings */

void 
set_warning_line(const char *fn,int line);
/* Set filename and linenumber for the warning */
  
int 
get_warning_line(void);
/* Get linenumber for the warning */
  

const char *
get_warning_file(void);
/* Get filename for the warning */
  
extern char 
warn_buf[1024];
/* Warning buffer of 1024 bytes, which can be used to print messages to */

void
warning(const char *s);
/* Issue a warning, with the string s. If s == NULL, then warn_buf
 * will be printed instead. The file and line set by set_warning_line
 * are printed, nwarn_warn (local) is incremented.
 * A fatal error will be generated after processing the input
 * when nwarn_warn is larger than maxwarning passed to init_warning.
 * So warning should only be called for issues that should be resolved,
 * otherwise warning_note should be called.
 */

void 
warning_note(const char *s);
/* Issue a note, with the string s. If s == NULL, then warn_buf
 * will be printed instead. The file and line set by set_warning_line
 * are printed, nwarn_note (local) is incremented.
 * This is for issues which could be a problem for some systems,
 * but 100% ok for other systems.
 */

void 
warning_error(const char *s);
/* Issue an error, with the string s. If s == NULL, then warn_buf
 * will be printed instead. The file and line set by set_warning_line
 * are printed, nwarn_error (local) is incremented.
 */
 
void 
check_warning_error(int f_errno,const char *file,int line);
/* When warning_error has been called at least once gmx_fatal is called,
 * otherwise does nothing.
 */

void 
print_warn_num(bool bFatalError);
/* Print the total number of warnings, if larger than 0.
 * When bFatalError == TRUE generates a fatal error
 * when the number is larger than maxwarn.
 */
  
void 
_too_few(const char *fn,int line);
#define too_few() _too_few(__FILE__,__LINE__)
/* Issue a warning stating 'Too few parameters' */

void 
_incorrect_n_param(const char *fn,int line);
#define incorrect_n_param() _incorrect_n_param(__FILE__,__LINE__)
/* Issue a warning stating 'Incorrect number of parameters' */
  
#ifdef __cplusplus
	   }
#endif

#endif	/* _warninp_h */
