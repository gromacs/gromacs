/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * Great Red Owns Many ACres of Sand 
 */

#include <sysstuff.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include "futil.h"
#include "statutil.h"
#include "main.h"
#include "network.h"
#include "fatal.h"
#include "macros.h"
#include "string2.h"
#include "smalloc.h"


static bool bDebug = FALSE;
static char *fatal_tmp_file = NULL;

bool bDebugMode(void)
{
  return bDebug;
}

void _where(char *file,int line)
{
  static bool bFirst = TRUE;
  static int  nskip  = -1;
  static int  nwhere =  0;
  FILE *fp;
  char *temp; 
  
  if ( bFirst ) {
    if ((temp=getenv("WHERE")) != NULL)
      nskip = atoi(temp);
    bFirst = FALSE;
  } 

  if (nskip >= 0) {
    /* Skip the first n occasions, this allows to see where it goes wrong */
    if (nwhere >= nskip) {
      if (stdlog)
	fp = stdlog;
      else
	fp = stderr;
      fprintf(fp,"WHERE %d, file %s - line %d\n",nwhere,file,line);
    }
    nwhere++;
  }
}

static void bputc(char *msg,int *len,char ch)
{
  msg[(*len)++]=ch;
}

static void bputs(char *msg,int *len,char *s,int fld)
{
  for (fld-=(int)strlen(s); fld>0; fld--) bputc(msg,len,' ');
  while (*s) bputc(msg,len,*(s++));
}

static void bputd(char *msg,int *len,int d)
{
  if (d<10) bputc(msg,len,d+'0'); else bputc(msg,len,d-10+'a');
}

static void bputi(char *msg,int *len,int val,int radix,int fld,bool bNeg)
{
  int fmax=0;
  
  if (bNeg)
    fmax=1;
    
  if (val<radix)
    {
      for (fld--; fld>fmax; fld--) 
	bputc(msg,len,' ');
      if (bNeg)
	bputc(msg,len,'-');
      bputd(msg,len,val);
    }
  else
    {
      if (bNeg)
	bputc(msg,len,'-');
      bputi(msg,len,val/radix,radix,fld-1,FALSE);
      bputd(msg,len,val%radix);
    }
}

static int getfld(char **p)
{
  int fld;

  fld=0;
  while (isdigit(**p)) fld=(fld*10)+((*((*p)++))-'0');
  return fld;
}

void _halt(char *file,int line,char *reason)
{
  fprintf(stderr,"\nHALT in file %s line %d because:\n\t%s\n",
	  file,line,reason);
  exit(1);
}

void quit_gmx(int fatal_errno,char *msg)
{
  if (!fatal_errno) {
    if (stdlog) 
      fprintf(stdlog,"%s\n",msg);
    fprintf(stderr,"%s\n",msg);
  }
  else {
    if (fatal_errno != -1)
      errno=fatal_errno;
    perror(msg);
  }
  
#ifdef USE_MPI
  if (gmx_parallel) {
    int  nnodes;
    int  nodeid;
    
    nnodes = gmx_node_num();
    nodeid    = gmx_node_id();
    
    if (nnodes > 1) 
      fprintf(stderr,"Error on node %d, will try to stop all the nodes\n",nodeid);
    gmx_abort(nodeid,nnodes,-1);
  }
#else
  if (debug)
    fflush(debug);
  if (bDebug) {
    fprintf(stderr,"dump core (y/n):"); 
    fflush(stderr);
    if (toupper(getc(stdin))!='N') 
      (void) abort(); 
  }
#endif 
  exit(-1);
}

void _set_fatal_tmp_file(char *fn, char *file, int line)
{
  if (fatal_tmp_file == NULL)
    fatal_tmp_file = strdup(fn);
  else
    fprintf(stderr,"BUGWARNING: fatal_tmp_file already set at %s:%d",
	    file,line);
}

void _unset_fatal_tmp_file(char *fn, char *file, int line)
{
  if (strcmp(fn,fatal_tmp_file) == 0) {
    sfree(fatal_tmp_file);
    fatal_tmp_file = NULL;
  } else
    fprintf(stderr,"BUGWARNING: file %s not set as fatal_tmp_file at %s:%d",
	    fn,file,line);
}

void fatal_error(int fatal_errno,char *fmt,...)
{
  va_list ap;
  char    *p,cval,*sval,msg[STRLEN];
  char    ibuf[64],ifmt[64];
  int     index,ival,fld,len;
  double  dval;
#ifdef _SPECIAL_VAR_ARG
  int     fatal_errno;
  char    *fmt;
  
  va_start(ap,);
  fatal_errno=va_arg(ap,int);
  fmt=va_arg(ap,char *);
#else
  va_start(ap,fmt);
#endif

  if (fatal_tmp_file) {
    fprintf(stderr,"Cleaning up temporary file %s\n",fatal_tmp_file);
    remove(fatal_tmp_file);
    sfree(fatal_tmp_file);
    fatal_tmp_file = NULL;
  }
  
  len=0;
  bputs(msg,&len,"Fatal error: ",0);
  for (p=fmt; *p; p++) {
    if (*p!='%')
      bputc(msg,&len,*p);
    else {
      p++;
      fld=getfld(&p);
      switch(*p) {
      case 'x':
	ival=va_arg(ap,int);
	sprintf(ifmt,"0x%%%dx",fld);
	sprintf(ibuf,ifmt,(unsigned int)ival);
	for(index=0; (index<(int)strlen(ibuf)); index++)
	  bputc(msg,&len,ibuf[index]);
	break;
      case 'd':
	ival=va_arg(ap,int);
	sprintf(ifmt,"%%%dd",fld);
	sprintf(ibuf,ifmt,ival);
	for(index=0; (index<(int)strlen(ibuf)); index++)
	  bputc(msg,&len,ibuf[index]);
	break;
      case 'f':
	dval=va_arg(ap,double);
	sprintf(ifmt,"%%%df",fld);
	sprintf(ibuf,ifmt,dval);
	for(index=0; (index<(int)strlen(ibuf)); index++)
	  bputc(msg,&len,ibuf[index]);
	break;
      case 'c':
	cval=(char) va_arg(ap,int); /* char is promoted to int */
	bputc(msg,&len,cval);
	break;
      case 's':
	sval=va_arg(ap,char *);
	bputs(msg,&len,sval,fld);
	break;
      default:
	break;
      }
    }
  }
  va_end(ap);
  bputc(msg,&len,'\0');

  quit_gmx(fatal_errno,msg);
}

static int  nwarn       = 0;
static int  maxwarn     = 10;
static int  lineno      = 1;
static char filenm[256] = "";
char   warn_buf[1024];

void init_warning(int maxwarning)
{
  maxwarn = maxwarning;
  nwarn   = 0;
}

void set_warning_line(char *s,int line)
{
  strcpy(filenm,s);
  lineno = line;
}

int get_warning_line()
{
  return lineno;
}

char *get_warning_file()
{
  return filenm;
}

void warning(char *s)
{
#define indent 2 
  char linenobuf[32], *temp, *temp2;
  int i;

  nwarn++;
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
  fprintf(stderr,"WARNING %d [file %s, line %s]:\n%s\n",
	  nwarn,filenm,linenobuf,temp2);
  sfree(temp);
  sfree(temp2);
  if (nwarn >= maxwarn)
    fatal_error(0,"Too many warnings, %s terminated",Program());
}

void print_warn_num(void)
{
  if (nwarn > 0)
    fprintf(stderr,"There %s %d warning%s\n",
	    (nwarn==1) ? "was" : "were", nwarn, (nwarn==1) ? "" : "s");
}

void _too_few(char *fn,int line)
{
  sprintf(warn_buf,"Too few parameters on line (source file %s, line %d)",
	  fn,line);
  warning(NULL);
}

void _invalid_case(char *fn,int line)
{
  fatal_error(0,"Invalid case in switch statement, file %s, line %d",
	      fn,line);
}

void _unexpected_eof(char *fn,int line,char *srcfn,int srcline)
{
  fatal_error(0,"Unexpected end of file in file %s at line %d\n"
	      "(Source file %s, line %d)",fn,line,srcfn,srcline);
}

/* 
 * These files are global variables in the gromacs preprocessor
 * Every routine in a file that includes fatal.h can write to these
 * debug channels. Depending on the debuglevel used
 * 0 to 3 of these filed are redirected to /dev/null
 *
 */
FILE *debug=NULL;

void init_debug (char *dbgfile)
{
  no_buffers();
  debug=ffopen(dbgfile,"w");
  bDebug = TRUE;
}

#if (defined __sgi && defined USE_SGI_FPE)
static void user_routine(unsigned us[5], int ii[2])
{
  fprintf(stderr,"User routine us=(%u,%u,%u,%u,%u) ii=(%d,%d)\n",
	  us[0],us[1],us[2],us[3],us[4],ii[0],ii[1]);
  fprintf(stderr,"Exception encountered! Dumping core\n");
  abort();
}

static void abort_routine(unsigned int **ii)
{
  fprintf(stderr,"Abort routine\n");
  abort();
}

static void handle_signals(int n)
{
  fprintf(stderr,"Handle signals: n = %d\n",n);
  fprintf(stderr,"Dumping core\n");
  abort();
}

void doexceptions(void)
{
#include <sigfpe.h>
#include <signal.h>
  int hs[] = { SIGILL, SIGFPE, SIGTRAP, SIGEMT, SIGSYS };
  
  int onoff,en_mask,abort_action,i;
  
  onoff   = _DEBUG;
  en_mask = _EN_UNDERFL | _EN_OVERFL | _EN_DIVZERO | 
    _EN_INVALID | _EN_INT_OVERFL;
  abort_action = _ABORT_ON_ERROR;
  handle_sigfpes(onoff,en_mask,user_routine,abort_action,abort_routine);
  
  for(i=0; (i<asize(hs)); i++)
    signal(hs[i],handle_signals);
}
#endif /* __sgi and FPE */




