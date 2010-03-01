/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sysstuff.h>
#include <string.h>
#include "smalloc.h"
#include "statutil.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "warninp.h"

#ifdef GMX_THREADS
#include "thread_mpi.h"
#endif

#ifdef GMX_THREADS
static tMPI_Thread_mutex_t warning_mutex=TMPI_THREAD_MUTEX_INITIALIZER;
#endif

static int  nwarn_note  = 0;
static int  nwarn_warn  = 0;
static int  nwarn_error = 0;
static int  maxwarn     = 10;
static int  lineno      = 1;
static char filenm[256] = "";
char   warn_buf[1024]   = "";

void init_warning(int maxwarning)
{
#ifdef GMX_THREADS
  tMPI_Thread_mutex_lock(&warning_mutex);
#endif
  maxwarn     = maxwarning;
  nwarn_note  = 0;
  nwarn_warn  = 0;
  nwarn_error = 0;
#ifdef GMX_THREADS
  tMPI_Thread_mutex_unlock(&warning_mutex);
#endif
}

void set_warning_line(const char *s,int line)
{
    if (s == NULL)
    {
        gmx_incons("Calling set_warning_line with NULL pointer");
    }
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&warning_mutex);
#endif

    strcpy(filenm,s);
    lineno = line;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&warning_mutex);
#endif
}

int get_warning_line()
{
    int ret;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&warning_mutex);
#endif
    ret=lineno;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&warning_mutex);
#endif
    return ret;
}

const char *get_warning_file()
{
    return filenm;
}

static void low_warning(const char *wtype,int n,const char *s)
{
#define indent 2 
  char linenobuf[32], *temp, *temp2;
  int i;

  if (s == NULL)
    s = warn_buf;
  if (lineno != -1)
    sprintf(linenobuf,"%d",lineno);
  else
    strcpy(linenobuf,"unknown");
  snew(temp,strlen(s)+indent+1);
  for(i=0; i<indent; i++)
    temp[i] = ' ';
  temp[indent] = '\0';
  strcat(temp,s);
  temp2 = wrap_lines(temp,78-indent,indent,FALSE);
  if (strlen(filenm) > 0) {
    fprintf(stderr,"\n%s %d [file %s, line %s]:\n%s\n\n",
	    wtype,n,filenm,linenobuf,temp2);
  } else {
    fprintf(stderr,"\n%s %d:\n%s\n\n",wtype,n,temp2);
  }
  sfree(temp);
  sfree(temp2);
}

void warning(const char *s)
{
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&warning_mutex);
#endif
    nwarn_warn++;
    low_warning("WARNING",nwarn_warn,s);
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&warning_mutex);
#endif
}

void warning_note(const char *s)
{
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&warning_mutex);
#endif
    nwarn_note++;
    low_warning("NOTE",nwarn_note,s);
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&warning_mutex);
#endif
}

void warning_error(const char *s)
{
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&warning_mutex);
#endif
    nwarn_error++;
    low_warning("ERROR",nwarn_error,s);
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&warning_mutex);
#endif
}

static void print_warn_count(const char *type,int n)
{
  if (n > 0) {
    fprintf(stderr,"\nThere %s %d %s%s\n",
	    (n==1) ? "was" : "were", n, type, (n==1) ? "" : "s");
  }
}

void print_warn_num(bool bFatalError)
{
    bool cond;
    int nww;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&warning_mutex);
#endif
    print_warn_count("note",nwarn_note);
    print_warn_count("warning",nwarn_warn);
    nww=nwarn_warn;
    cond=bFatalError && nwarn_warn > maxwarn;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&warning_mutex);
#endif

    if (cond) {
        gmx_fatal(FARGS,"Too many warnings (%d), %s terminated.\n"
                  "If you are sure all warnings are harmless, use the -maxwarn option.",
                  nww,Program());
    }
}

void check_warning_error(int f_errno,const char *file,int line)
{
    int nwe;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_lock(&warning_mutex);
#endif
    nwe=nwarn_error;
#ifdef GMX_THREADS
    tMPI_Thread_mutex_unlock(&warning_mutex);
#endif
    if (nwe > 0) {
        print_warn_num(FALSE);
        gmx_fatal(f_errno,file,line,"There %s %d error%s in input file(s)",
                  (nwe==1) ? "was" : "were",nwe,
                  (nwe==1) ? ""    : "s");
    }
}

void _too_few(const char *fn,int line)
{
  sprintf(warn_buf,"Too few parameters on line (source file %s, line %d)",
	  fn,line);
  warning(NULL);
}

void _incorrect_n_param(const char *fn,int line)
{
  sprintf(warn_buf,"Incorrect number of parameters on line (source file %s, line %d)",
	  fn,line);
  warning(NULL);
}
