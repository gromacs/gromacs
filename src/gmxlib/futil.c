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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_futil_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "string2.h"
#include "futil.h"
#include "network.h"
#include "fatal.h"
#include "smalloc.h"

typedef struct t_pstack {
  FILE   *fp;
  struct t_pstack *prev;
} t_pstack;

static t_pstack *pstack=NULL;
static bool     bUnbuffered=FALSE;

void no_buffers(void)
{
  bUnbuffered=TRUE;
}

void push_ps(FILE *fp)
{
  t_pstack *ps;
 
  snew(ps,1);
  ps->fp   = fp;
  ps->prev = pstack;
}

#ifdef fclose
#undef fclose
#endif

void ffclose(FILE *fp)
{
  t_pstack *ps,*tmp;
  
  ps=pstack;
  if (ps == NULL) {
    fclose(fp);
    return;
  }
  if (ps->fp == fp) {
    pclose(fp);
    pstack=pstack->prev;
    sfree(ps);
  }
  else {
    while ((ps->prev != NULL) && (ps->prev->fp != fp))
      ps=ps->prev;
    if (ps->prev->fp == fp) {
      pclose(ps->prev->fp);
      tmp=ps->prev;
      ps->prev=ps->prev->prev;
      sfree(tmp);
    }
    else {
      fclose(fp);
      return;
    }
  }
}

#ifdef rewind
#undef rewind
#endif

void frewind(FILE *fp)
{
  t_pstack *ps;
  
  ps=pstack;
  while (ps != NULL) {
    if (ps->fp == fp) {
      fprintf(stderr,"Cannot rewind compressed file!\n");
      return;
    }
    ps=ps->prev;
  }
  rewind(fp);
}

bool is_pipe(FILE *fp)
{
  t_pstack *ps;
  
  ps=pstack;
  while (ps != NULL) {
    if (ps->fp == fp) {
      return TRUE;
    }
    ps=ps->prev;
  }
  return FALSE;
}

#ifdef NO_PIPE
static FILE *popen(char *nm,char *mode)
{
  fatal_error(0,"Sorry no pipes...");
  
  return NULL;
}

static int pclose(FILE *fp)
{
  fatal_error(0,"Sorry no pipes...");
  
  return 0;
}
#endif


FILE *uncompress(char *fn,char *mode)
{
  FILE *fp;
  char buf[256];
  
  sprintf(buf,"uncompress -c < %s",fn);
  fprintf(stderr,"Going to execute '%s'\n",buf);
  if ((fp=popen(buf,mode)) == NULL)
    fatal_error(0,"Could not open %s",fn);
  push_ps(fp);
  
  return fp;
}

FILE *gunzip(char *fn,char *mode)
{
  FILE *fp;
  char buf[256];
  
  sprintf(buf,"gunzip -c < %s",fn);
  fprintf(stderr,"Going to execute '%s'\n",buf);
  if ((fp=popen(buf,mode)) == NULL)
    fatal_error(0,"Could not open %s",fn);
  push_ps(fp);
  
  return fp;
}

bool fexist(char *fname)
{
  FILE *test;
  
  if (fname == NULL)
    return FALSE;
  test=fopen(fname,"r");
  if (test == NULL) 
    return FALSE;
  else {
    fclose(test);
    return TRUE;
  }
}

bool eof(FILE *fp)
{
  char data[4];
  bool beof;

  if (is_pipe(fp))
    return feof(fp);
  else {
    if ((beof=fread(data,1,1,fp))==1)
      fseek(fp,-1,SEEK_CUR);
    return !beof;
  }
}

char *backup_fn(char *file)
{
  int         i;
  char        *ptr;
  static char buf[256];
  
  for(i=strlen(file)-1; ((i > 0) && (file[i] != '/')); i--);
  if (i > 0) {
    ptr=strdup(file);
    ptr[i]='\0';
    sprintf(buf,"%s/#%s#",ptr,ptr+i+1);
    sfree(ptr);
  }
  else
    sprintf(buf,"#%s#",file);
  
  return buf;
}

FILE *ffopen(char *file,char *mode)
{
  FILE *ff=NULL;
  char buf[256],*bf,*bufsize,*ptr;
  bool bRead;
  int  bs;
  
#ifdef _amb_
  fprintf(stderr,"Going to open %s on NODE %d with mode %s\n",
	  file,gmx_node_id(),mode);
#endif
  if ((mode[0]=='w') && fexist(file)) {
    bf=backup_fn(file);
    if (rename(file,bf) == 0) {
      fprintf(stderr,"\nBack Off! I just backed up %s to %s\n",file,bf);
    }
    else
      fprintf(stderr,"Sorry, I couldn't backup %s to %s\n",file,bf);
  }
  where();
  
  bRead= mode[0]=='r';
  strcpy(buf,file);
  if (fexist(buf) || !bRead) {
    if ((ff=fopen(buf,mode))==NULL)
      fatal_error(0,"Could not open %s",buf);
    where();
    /* Check whether we should be using buffering (default) or not
     * (for debugging)
     */
    if (bUnbuffered || ((bufsize=getenv("LOG_BUFS")) != NULL)) {
      /* Check whether to use completely unbuffered */
      if (bUnbuffered)
	bs = 0;
      else
	bs=atoi(bufsize);
      if (bs <= 0)
	setbuf(ff,NULL); 
      else {
	snew(ptr,bs+8);
	if (setvbuf(ff,ptr,_IOFBF,bs) != 0)
	  fatal_error(0,"Buffering File");
      }
    }
    where();
  }
  else {
    sprintf(buf,"%s.Z",file);
    if (fexist(buf)) {
      ff=uncompress(buf,mode);
    }
    else {
      sprintf(buf,"%s.gz",file);
      if (fexist(buf)) {
	ff=gunzip(buf,mode);
      }
      else 
	fatal_error(0,"%s does not exist",file);
    }
  }
#ifdef _amb_
  fprintf(stderr,"Opened %s on NODE %d with mode %s\n",
	  file,gmx_node_id(),mode);
#endif
  return ff;
}

char *low_libfn(char *file,bool bFatal)
{
  static char *libdir="GMXLIB";
  char *ret=NULL,*lib;
  static char buf[1024];
  
  if (fexist(file))
    ret=file;
  else {
    if ((lib=getenv(libdir)) == NULL) {
      if (bFatal)
	fatal_error(0,"%s environment variable not set and file %s not found",
		    libdir,file);
      else 
	return NULL;
    }
    else {
      sprintf(buf,"%s/%s",lib,file);
      ret=buf;
      if (bFatal && !fexist(ret))
	fatal_error(0,"Library file %s found in current dir nor in libdir %s",
		    ret,lib);
    }
  }
    
  return ret;
}

FILE *low_libopen(char *file,bool bFatal)
{
  FILE *ff;
  char *fn;

  fn=low_libfn(file,bFatal);

  if (fn==NULL) {
    ff=NULL;
  } else {
    if (bFatal)
      fprintf(stderr,"Opening library file %s\n",fn);
    ff=fopen(fn,"r");
  }

  return ff;
}

char *libfn(char *file)
{
  return low_libfn(file,TRUE);
}

FILE *libopen(char *file)
{
  return low_libopen(file,TRUE);
}

