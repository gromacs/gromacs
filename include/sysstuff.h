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
 * GRoups of Organic Molecules in ACtion for Science
 */

#ifndef _sysstuff_h
#define _sysstuff_h

static char *SRCID_sysstuff_h = "$Id$";
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef HAVE_IDENT
#ident	"@(#) sysstuff.h 1.22 12/16/92"
#endif /* HAVE_IDENT */

#ifdef CPLUSPLUS
extern "C" { 
#endif


#ifndef NO_STDINCLUDE

#ifndef _386_
#include <stdlib.h>
#endif
#include <stdio.h>
#include <errno.h>
#include <signal.h>
#ifndef _860_
#include <unistd.h>
#endif
#include <limits.h>
#include <time.h>

#else

extern	int	kill(int pid,int sig);
extern	int	fork();
extern	int	setpgrp();
extern	int	execv (char *path,char *argv[]);
extern	int	execvp(char *path,char *argv[]);
extern  int     execlp(char *path,...);
extern  int     getpid();
extern  int     getgid();
extern  int     _filbuf(FILE *stream);
extern void	perror(char *s);
extern int	printf(char *format,...);
extern int	fprintf(FILE *stream,char *format,...);
extern void	fflush(FILE *stream);
extern int	fread (char *ptr,int size,int nitems,FILE *stream);
extern int	fwrite(char *ptr,int size,int nitems,FILE *stream);
extern void	ungetc(char c,FILE *stream);
extern int	sscanf(char *s,char *format,...);
extern int	fscanf(FILE *,char *format,...);
extern int	scanf(char *format,...);
extern int	close(int fd);
extern int	fclose(FILE *stream);
extern int	fseek(FILE *stream,long offset,int ptrname);
extern int	read(int fd,char *buf,int nbyte);
extern int	write(int fd,char *buf,int nbyte);
extern int      ioctl(int fildes, int request, void *buf);
extern int	unlink(char *path);

#endif

#ifdef CPLUSPLUS
}
#endif

#endif	/* _sysstuff_h */
