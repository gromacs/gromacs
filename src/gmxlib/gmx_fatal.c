/*
 * $Id$
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
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include "futil.h"
#include "statutil.h"
#include "main.h"
#include "network.h"
#include "gmx_fatal.h"
#include "copyrite.h"
#include "macros.h"
#include "string2.h"
#include "smalloc.h"


static bool bDebug = FALSE;
static char *fatal_tmp_file = NULL;
static FILE *log_file = NULL;

bool bDebugMode(void)
{
  return bDebug;
}

void gmx_fatal_set_log_file(FILE *fp)
{
  log_file = fp;
}

void _where(const char *file,int line)
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
      if (log_file)
	fp = log_file;
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

static void bputs(char *msg,int *len,const char *s,int fld)
{
  for (fld-=(int)strlen(s); fld>0; fld--) 
    bputc(msg,len,' ');
  while (*s) 
    bputc(msg,len,*(s++));
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

static int getfld(const char **p)
{
  int fld;

  fld=0;
  while (isdigit(**p)) fld=(fld*10)+((*((*p)++))-'0');
  return fld;
}

/*static void _halt(char *file,int line,char *reason)
{
  fprintf(stderr,"\nHALT in file %s line %d because:\n\t%s\n",
	  file,line,reason);
  exit(1);
}
*/

static int fatal_errno = 0;

static void quit_gmx(const char *msg)
{
  if (!fatal_errno) {
    if (log_file) 
      fprintf(log_file,"%s\n",msg);
    fprintf(stderr,"%s\n",msg);
  }
  else {
    if (fatal_errno != -1)
      errno=fatal_errno;
    perror(msg);
  }
  
  if (gmx_parallel_env) {
    int  nnodes;
    int  noderank;
    
    nnodes   = gmx_node_num();
    noderank = gmx_node_rank();
    
    if (nnodes > 1) 
      fprintf(stderr,"Error on node %d, will try to stop all the nodes\n",
	      noderank);
    gmx_abort(noderank,nnodes,-1);
  } else {
    if (debug)
      fflush(debug);
    if (bDebug) {
      fprintf(stderr,"dump core (y/n):"); 
      fflush(stderr);
      if (toupper(getc(stdin))!='N') 
	(void) abort(); 
    }
  }

  exit(-1);
}

void _set_fatal_tmp_file(const char *fn, const char *file, int line)
{
  if (fatal_tmp_file == NULL)
    fatal_tmp_file = strdup(fn);
  else
    fprintf(stderr,"BUGWARNING: fatal_tmp_file already set at %s:%d",
	    file,line);
}

void _unset_fatal_tmp_file(const char *fn, const char *file, int line)
{
  if (strcmp(fn,fatal_tmp_file) == 0) {
    sfree(fatal_tmp_file);
    fatal_tmp_file = NULL;
  } else
    fprintf(stderr,"BUGWARNING: file %s not set as fatal_tmp_file at %s:%d",
	    fn,file,line);
}

static void clean_fatal_tmp_file()
{
  if (fatal_tmp_file) {
    fprintf(stderr,"Cleaning up temporary file %s\n",fatal_tmp_file);
    remove(fatal_tmp_file);
    sfree(fatal_tmp_file);
    fatal_tmp_file = NULL;
  }
}

/* Old function do not use */
static void fatal_error(int f_errno,const char *fmt,...)
{
  va_list ap;
  const char    *p;
  char    cval,*sval,msg[STRLEN];
  char    ibuf[64],ifmt[64];
  int     index,ival,fld,len;
  double  dval;
#ifdef _SPECIAL_VAR_ARG
  int     f_errno;
  char    *fmt;
  
  va_start(ap,);
  f_errno=va_arg(ap,int);
  fmt=va_arg(ap,char *);
#else
  va_start(ap,fmt);
#endif

  clean_fatal_tmp_file();
  
  len=0;
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

  fatal_errno = f_errno;
  gmx_error("fatal",msg);
}

void gmx_fatal(int f_errno,const char *file,int line,const char *fmt,...)
{
  va_list ap;
  const char    *p;
  char    cval,*sval,msg[STRLEN];
  char    ibuf[64],ifmt[64];
  int     index,ival,fld,len;
  double  dval;
#ifdef _SPECIAL_VAR_ARG
  int     f_errno,line;
  char    *fmt,*file;
  
  va_start(ap,);
  f_errno = va_arg(ap,int);
  file    = va_arg(ap,char *);
  line    = va_arg(ap,int);
  fmt     = va_arg(ap,char *);
#else
  va_start(ap,fmt);
#endif

  clean_fatal_tmp_file();
  
  len=0;
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
      case 'u':
	ival=va_arg(ap,unsigned);
	sprintf(ifmt,"%%%du",fld);
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
      case 'g':
	dval=va_arg(ap,double);
	sprintf(ifmt,"%%%dg",fld);
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
      case '%':
	bputc(msg,&len,*p);
	break;
      default:
	break;
      }
    }
  }
  va_end(ap);
  bputc(msg,&len,'\0');

  fatal_errno = f_errno;
  
  _gmx_error("fatal",msg,file,line);
}

static int  nwarn_note  = 0;
static int  nwarn_warn  = 0;
static int  nwarn_error = 0;
static int  maxwarn     = 10;
static int  lineno      = 1;
static char filenm[256] = "";
char   warn_buf[1024]   = "";

void init_warning(int maxwarning)
{
  maxwarn     = maxwarning;
  nwarn_note  = 0;
  nwarn_warn  = 0;
  nwarn_error = 0;
}

void set_warning_line(const char *s,int line)
{
  if (s == NULL)
    gmx_incons("Calling set_warning_line with NULL pointer");
  strcpy(filenm,s);
  lineno = line;
}

int get_warning_line()
{
  return lineno;
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
  nwarn_warn++;
  low_warning("WARNING",nwarn_warn,s);
}

void warning_note(const char *s)
{
  nwarn_note++;
  low_warning("NOTE",nwarn_note,s);
}

void warning_error(const char *s)
{
  nwarn_error++;
  low_warning("ERROR",nwarn_error,s);
}

static void print_warn_count(char *type,int n)
{
  if (n > 0) {
    fprintf(stderr,"\nThere %s %d %s%s\n",
	    (n==1) ? "was" : "were", n, type, (n==1) ? "" : "s");
  }
}

void print_warn_num(bool bFatalError)
{
  print_warn_count("note",nwarn_note);

  print_warn_count("warning",nwarn_warn);

  if (bFatalError && nwarn_warn > maxwarn) {
    gmx_fatal(FARGS,"Too many warnings (%d), %s terminated.\n"
	      "If you are sure all warnings are harmless, use the -maxwarn option.",nwarn_warn,Program());
  }
}

void check_warning_error(int f_errno,const char *file,int line)
{
  if (nwarn_error > 0) {
    print_warn_num(FALSE);
    gmx_fatal(f_errno,file,line,"There %s %d error%s in input file(s)",
	      (nwarn_error==1) ? "was" : "were",nwarn_error,
	      (nwarn_error==1) ? ""    : "s");
  }
}

void _too_few(const char *fn,int line)
{
  sprintf(warn_buf,"Too few parameters on line (source file %s, line %d)",
	  fn,line);
  warning(NULL);
}

void _invalid_case(const char *fn,int line)
{
  gmx_fatal(FARGS,"Invalid case in switch statement, file %s, line %d",
	      fn,line);
}

void _unexpected_eof(const char *fn,int line,const char *srcfn,int srcline)
{
  gmx_fatal(FARGS,"Unexpected end of file in file %s at line %d\n"
	      "(Source file %s, line %d)",fn,line,srcfn,srcline);
}

/* 
 * These files are global variables in the gromacs preprocessor
 * Every routine in a file that includes gmx_fatal.h can write to these
 * debug channels. Depending on the debuglevel used
 * 0 to 3 of these filed are redirected to /dev/null
 *
 */
FILE *debug=NULL;
bool gmx_debug_at=FALSE;

void init_debug (const int dbglevel,const char *dbgfile)
{
  no_buffers();
  debug=ffopen(dbgfile,"w");
  bDebug = TRUE;
  if (dbglevel >= 2)
    gmx_debug_at = TRUE;
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

static char *gmxuser = "Please report this to the mailing list (gmx-users@gromacs.org)";

static void (*gmx_error_handler)(const char *msg) = quit_gmx;

void set_gmx_error_handler(void (*func)(const char *msg))
{
  gmx_error_handler = func;
}

char *gmx_strerror(const char *key)
{
  typedef struct {
    char *key,*msg;
  } error_msg_t;
  error_msg_t msg[] = {
    { "bug",    "Possible bug" },
    { "call",   "Routine should not have been called" },
    { "comm",   "Communication (parallel processing) problem" },
    { "fatal",  "Fatal error" },
    { "cmd",    "Invalid command line argument" },
    { "file",   "File input/output error" },
    { "impl",   "Implementation restriction" },
    { "incons", "Software inconsistency error" },
    { "input",  "Input error or input inconsistency" },
    { "mem",    "Memory allocation/freeing error" },
    { "open",   "Can not open file" },
    { "range",  "Range checking error" }
  };
#define NMSG asize(msg)
  char buf[1024];
  int  i;
  
  if (key == NULL)
    return strdup("Empty message");
  else {
    for(i=0; (i<NMSG); i++)
      if (strcmp(key,msg[i].key) == 0) 
	break;
    if (i == NMSG) {
      sprintf(buf,"No error message associated with key %s\n%s",key,gmxuser);
      return strdup(buf);
    }
    else
      return strdup(msg[i].msg);
  }
}

void _gmx_error(const char *key,const char *msg,const char *file,int line)
{
  char buf[10240],tmpbuf[1024];
  int  cqnum;

  /* protect the audience from suggestive discussions */
  char *lines = "-------------------------------------------------------";
  
  cool_quote(tmpbuf,1023,&cqnum);
  sprintf(buf,"\n%s\nProgram %s, %s\n"
	  "Source code file: %s, line: %d\n\n"
	  "%s:\n%s\n%s\n\n%s\n",
	  lines,ShortProgram(),GromacsVersion(),file,line,
	  gmx_strerror(key),msg ? msg : warn_buf,lines,tmpbuf);
  
  gmx_error_handler(buf);
}

void _range_check(int n,int n_min,int n_max,const char *var,const char *file,int line)
{
  char buf[1024];
  
  if ((n < n_min) || (n >= n_max)) {
    if (strlen(warn_buf) > 0) {
      strcpy(buf,warn_buf);
      strcat(buf,"\n");
      warn_buf[0] = '\0';
    }
    else
      buf[0] = '\0';
    
    sprintf(buf+strlen(buf),"Variable %s has value %d. It should have been "
	    "within [ %d .. %d ]\n",var,n,n_min,n_max);
    
    _gmx_error("range",buf,file,line);
  }
}

