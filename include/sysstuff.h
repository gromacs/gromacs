/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * S  C  A  M  O  R  G
 */

#ifndef	_sysstuff_h
#define	_sysstuff_h

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
