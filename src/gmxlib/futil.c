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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_futil_c = "$Id$";
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "sysstuff.h"
#include "string2.h"
#include "futil.h"
#include "network.h"
#include "fatal.h"
#include "smalloc.h"
#include "statutil.h"

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
  char buf[256],*bf,*bufsize=0,*ptr;
  bool bRead;
  int  bs;
  
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
  return ff;
}



bool search_subdirs(char *parent, char *libdir)
{
  char *ptr;
  bool found;

  /* Search a few common subdirectory names for the gromacs library dir */
  sprintf(libdir,"%s/share/top/FF.dat",parent);
  found=fexist(libdir);
  if(!found) {
    sprintf(libdir,"%s/share/gromacs/top/FF.dat",parent);
    found=fexist(libdir);
  }    
  if(!found) {
    sprintf(libdir,"%s/share/gromacs-%s/top/FF.dat",parent,VERSION);
    found=fexist(libdir);
  }    
  if(!found) {
    sprintf(libdir,"%s/share/gromacs/gromacs-%s/top/FF.dat",parent,VERSION);
    found=fexist(libdir);
  }    
  
  /* Remove the FF.dat part from libdir if we found something */
  if(found) {
    ptr=strrchr(libdir,'/'); /* slash always present, no check necessary */
    *ptr='\0';
  }
  return found;
}

bool get_libdir(char *libdir)
{
  char bin_name[512];
  char buf[512];
  char full_path[512];
  char test_file[512];
  char *system_path;
  char *dir,*ptr,*s;
  bool found=FALSE;
  int i;

  /* First - detect binary name */
  strcpy(bin_name,Program());
  
  /* Only do the smart search part if we got a real name */
  if (bin_name && strcmp(bin_name,"GROMACS")) {
  
    if (!strchr(bin_name,'/')) {
      /* No "/" in name means it must be in the path - search it! */
      system_path=getenv("PATH");
      s=system_path;
      found=FALSE;
      while (!found && (dir=strtok(s,":"))!=NULL) {
	sprintf(full_path,"%s/%s",dir,bin_name);
	found=fexist(full_path);
	s=NULL; /* pointer should be null for subseq. calls to strtok */
      }
      if (!found)
	return FALSE;
    } else if (bin_name[0] != '/') {
      /* name is relative the current dir */
      getcwd(buf,sizeof(buf)-1);
      strcpy(full_path,buf);
      strcat(full_path,bin_name+1);
    } else {
      strcpy(full_path,bin_name);
    }

    /* Now we should have a full path and name in full_path,
     * but it might be a link, or a link to a link to a link...
     */
    while( (i=readlink(full_path,buf,sizeof(buf)-1)) > 0 ) {
      buf[i]='\0';
      /* If it doesn't start with "/" it is relative */
      if (buf[0]!='/') {
	strcpy(strrchr(full_path,'/')+1,buf);
      } else
	strcpy(full_path,buf);
    }
    
    /* Remove the executable name - it always contains at least one slash */
    *(strrchr(full_path,'/')+1)='\0';
    /* Now we have the full path to the gromacs executable.
     * Use it to find the library dir. 
     */
    found=FALSE;
    while(!found && ( (ptr=strrchr(full_path,'/')) != NULL ) ) {
      *ptr='\0';
      found=search_subdirs(full_path,libdir);
    }
  }
  /* End of smart searching. If we didn't find it in our parent tree,
   * or if the program name wasn't set, at least try some standard 
   * locations before giving up, in case we are running from e.g. 
   * a users home directory:
   */
  if(!found) 
    found=search_subdirs("/usr/local",libdir);
  if(!found) 
    found=search_subdirs("/usr",libdir);
  if(!found) 
    found=search_subdirs("/opt",libdir);
  return found;
}


char *low_libfn(char *file, bool bFatal)
{
  char *ret=NULL;
  char *lib,*dir;
  static char buf[1024];
  static char libpath[4096];
  static int  bFirst=1;
  static bool env_is_set;
  char   *s,tmppath[4096];
  bool found;
  
  if (bFirst) {
    /* GMXLIB can be a path now */
    lib=getenv("GMXLIB");
    if (lib != NULL) {
      env_is_set=TRUE;
      strcpy(libpath,lib);
    } 
    else if (!get_libdir(libpath))
      strcpy(libpath,GMXLIBDIR);
    
    bFirst=0;
  }

  if (fexist(file))
    ret=file;
  else {
    found=FALSE;
    strcpy(tmppath,libpath);
    s=tmppath;
    while(!found && (dir=strtok(s,":"))!=NULL) {
      sprintf(buf,"%s/%s",dir,file);
      found=fexist(buf);
      s = NULL;
    }
    ret=buf;
    if (bFatal && !found) {
      if(env_is_set) 
	fatal_error(0,"Library file %s not found in current dir nor in your GMXLIB path.\n",file);
      else
	fatal_error(0,"Library file %s not found in current dir nor in default directories.\n"
		    "(You can set the directories to search with the GMXLIB path variable.)\n",file);
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

void gmx_tmpnam(char *buf)
{
  int i,len,fd;
  
  if ((len = strlen(buf)) < 7)
    fatal_error(0,"Buf passed to gmx_tmpnam must be at least 7 bytes long");
  for(i=len-6; (i<len); i++) {
    buf[i] = 'X';
  }
  fd = mkstemp(buf);
  if (fd == EINVAL)
    fatal_error(0,"Invalid template %s for mkstemp (source %s, line %d)",buf,__FILE__,__LINE__);
  else if (fd == EEXIST)
    fatal_error(0,"mkstemp created existing file %s (source %s, line %d)",buf,__FILE__,__LINE__);
  close(fd);
  /* Buf should now be OK */
}
