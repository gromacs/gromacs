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
 * GROningen MAchine for Chemical Simulation
 */
static char *SRCID_fatal_c = "$Id$";

#include <sysstuff.h>
#include <ctype.h>
#include <errno.h>
#include <stdarg.h>
#include <string.h>
#include "statutil.h"
#include "main.h"
#include "network.h"
#include "fatal.h"
#include "macros.h"

static bool bDebug = FALSE;

bool bDebugMode(void)
{
  return bDebug;
}

void _where(char *file,int line)
{
  static int where =-1;
  static int nw    = 1;
  
  if ( where == -1 ) {
    char *temp; 
    if ((temp=getenv("WHERE")) != NULL)
      where = atoi(temp);
  } 

  if ( where != -1 ) {
    
    /* Skip the first n occasions, this allows to see where it goes wrong */
    if (where > 1) 
      where--;
    else {
      fprintf(stdlog,"WHERE %d, file %s - line %d\n",nw,file,line);
#ifdef _amb_
      fprintf(stdlog,"%512s\n","flushed");
#endif    
    }
    nw++;
  }
}

#define MSGSIZE		256

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
    errno=fatal_errno;
    perror(msg);
  }
  
#ifdef PARALLEL
  {
    int  nprocs;
    int  pid;
    
    nprocs = gmx_cpu_num();
    pid    = gmx_cpu_id();
    
    if (nprocs > 1) 
      fprintf(stderr,"Error on processor %d, will try to stop all the processors\n",pid);
    gmx_abort(pid,nprocs,-1);
  }
#else
  if (bDebug) {
    (void) fprintf(stderr,"dump core (y/n):"); 
    fflush(stderr);
    if (toupper(getc(stdin))=='Y') 
      (void) abort(); 
  }
#endif 
  exit(-1);
}

void fatal_error(int fatal_errno,char *fmt,...)
{
  FILE    *log;
  va_list ap;
  char    *p,cval,*sval,msg[MSGSIZE];
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
	cval=va_arg(ap,char);
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
static char filenm[256] = "\0";
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

void warning(char *s)
{
  char linenobuf[32];
  
  nwarn++;
  if (s == NULL)
    s = warn_buf;
  if (lineno != -1)
    sprintf(linenobuf,"%d",lineno);
  else
    strcpy(linenobuf,"unknown");
  fprintf(stderr,"Warning %d [file %s, line %s]: %s\n",
	  nwarn,filenm,linenobuf,s);
  if (nwarn >= maxwarn) {
    fprintf(stderr,"Too many warnings, %s terminated\n",Program());
    exit(1);
  }
}

void print_warn_num(void)
{
  if (nwarn > 0)
    fprintf(stderr,"There were %d warnings\n",nwarn);
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
  debug=ffopen(dbgfile,"w");
  bDebug = TRUE;
}

#ifdef _SGI_
static void user_routine(unsigned us[5], int ii[2])
{
  fprintf(stderr,"User routine us=(%d,%d,%d,%d,%d) ii=(%d,%d)\n",
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

#ifdef USE_SGI_FPE
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
#endif
#endif
