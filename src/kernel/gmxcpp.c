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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <limits.h>
#include <ctype.h>

#include "string2.h"
#include "smalloc.h"
#include "futil.h"
#include "macros.h"
#include "gmxcpp.h"

typedef struct {
  char *name;
  char *def;
} t_define;

static int      ndef   = 0;
static t_define *defs  = NULL;
static int      nincl  = 0;
static char     **incl = 0;

/* enum used for handling ifdefs */
enum { eifTRUE, eifFALSE, eifIGNORE, eifNR };

typedef struct gmx_cpp {
  FILE     *fp;
  char     *fn;
  int      line_len;
  char     *line;
  int      line_nr;
  int      nifdef;
  int      *ifdefs;
  struct   gmx_cpp *child,*parent;
} gmx_cpp;

static bool is_blank_end(char c)
{
  return (c == ' ' || c == '\t' || c == '\0' || c == '\n');
}

static char* strstrw(const char *buf,const char *word)
{
  char *ptr;

  ptr = strstr(buf,word);
  if (ptr) {
    /* Check if we did not find part of a longer word */
    if ((ptr > buf && !is_blank_end(ptr[-1])) ||
	!is_blank_end(ptr[strlen(word)]))
      ptr = NULL;
  }

  return ptr;
}

static void add_include(char *include)
{
  int i;
  
  if (include == NULL)
    return;
    
  for(i=0; (i<nincl); i++)
    if (strcmp(incl[i],include) == 0)
      break;
  if (i == nincl) {
    nincl++;
    srenew(incl,nincl);
    incl[nincl-1] = strdup(include);
  }
}

static void add_define(char *define)
{
  int  i;
  char *ptr,name[256];
  
  sscanf(define,"%s",name);
  ptr = define + strlen(name);
  
  while ((*ptr != '\0') && isspace(*ptr))
    ptr++;
    
  for(i=0; (i<ndef); i++) {
    if (strcmp(defs[i].name,name) == 0) {
      break;
    }
  }
  if (i == ndef) {
    ndef++;
    srenew(defs,ndef);
    i = ndef - 1;
    defs[i].name = strdup(name);
  }
  else if (defs[i].def) {
    if (debug)
      fprintf(debug,"Overriding define %s\n",name);
    sfree(defs[i].def);
  }
  if (strlen(ptr) > 0)
    defs[i].def  = strdup(ptr);
  else
    defs[i].def  = NULL;
}

/* Open the file to be processed. The handle variable holds internal
   info for the cpp emulator. Return integer status */
int cpp_open_file(char *filenm,gmx_cpp_t *handle,char **cppopts)
{
  gmx_cpp_t cpp;
  char *buf;
  int i;
  
  /* First process options, they might be necessary for opening files
     (especially include statements). */  
  i  = 0;
  if (cppopts) {
    while(cppopts[i]) {
      if (strstr(cppopts[i],"-I") == cppopts[i])
	add_include(cppopts[i]+2);
      if (strstr(cppopts[i],"-D") == cppopts[i])
	add_define(cppopts[i]+2);
      i++;
    }
  }
  if (debug)
    fprintf(debug,"Added %d command line arguments",i);
  
  snew(cpp,1);
  *handle      = (void *)cpp;
  cpp->fn      = strdup(filenm);
  cpp->line_len= 0;
  cpp->line    = NULL;
  cpp->line_nr = 0;
  cpp->nifdef  = 0;
  cpp->ifdefs  = NULL;
  cpp->child   = NULL;
  cpp->parent  = NULL;
  i = 0;
  while (((cpp->fp = fopen(cpp->fn,"r")) == NULL) && (i<nincl)) {
    snew(buf,strlen(incl[i])+strlen(filenm)+2);
    sprintf(buf,"%s/%s",incl[i],filenm);
    sfree(cpp->fn);
    cpp->fn = strdup(buf);
    sfree(buf);
    i++;
  }
  if (cpp->fp == NULL) {
    sfree(cpp->fn);
    cpp->fn = strdup(filenm);
    cpp->fp = libopen(filenm);
  }
  if (cpp->fp == NULL) {
    switch(errno) {
    case EINVAL:
    default:
      return eCPP_UNKNOWN;
    }
  }
  return eCPP_OK;
}

/* Return one whole line from the file into buf which holds at most n
   characters, for subsequent processing. Returns integer status. This
   routine also does all the "intelligent" work like processing cpp
   directives and so on. Note that often the routine is called
   recursively and no cpp directives are printed. */
int cpp_read_line(gmx_cpp_t *handlep,int n,char buf[])
{
  gmx_cpp_t handle = (gmx_cpp_t)*handlep;
  int  i,i0,nn,len,status;
  char *inc_fn,*ptr,*ptr2,*name;
  int  bIfdef,bIfndef;
  
  if (!handle)
    return eCPP_INVALID_HANDLE;
  if (!handle->fp)
    return eCPP_FILE_NOT_OPEN;
    
  if (feof(handle->fp) || (fgets2(buf,n-1,handle->fp) == NULL)) {
    if (handle->parent == NULL)
      return eCPP_EOF;
    cpp_close_file(handlep);
    *handlep = handle->parent;
    handle->child = NULL;
    return cpp_read_line(handlep,n,buf);
  }
  else {
    if (n > handle->line_len) {
      handle->line_len = n;
      srenew(handle->line,n);
    }
    strcpy(handle->line,buf);
    handle->line_nr++;
  }
  /* Now we've read a line! */
  if (debug) 
    fprintf(debug,"%s : %4d : %s\n",handle->fn,handle->line_nr,buf);
  set_warning_line(handle->fn,handle->line_nr);
    
  /* #ifdef or ifndef statement */
  bIfdef  = (strstrw(buf,"#ifdef") != NULL);
  bIfndef = (strstrw(buf,"#ifndef") != NULL);
  if (bIfdef || bIfndef) {
    if ((handle->nifdef > 0) && (handle->ifdefs[handle->nifdef-1] != eifTRUE)) {
      handle->nifdef++;
      srenew(handle->ifdefs,handle->nifdef);
      handle->ifdefs[handle->nifdef-1] = eifIGNORE;
    }
    else {
      snew(name,strlen(buf));
      status = sscanf(buf,"%*s %s",name);
      for(i=0; (i<ndef); i++) 
	if (strcmp(defs[i].name,name) == 0) 
	  break;
      handle->nifdef++;
      srenew(handle->ifdefs,handle->nifdef);
      if ((bIfdef && (i < ndef)) || (bIfndef && (i == ndef))) 
	handle->ifdefs[handle->nifdef-1] = eifTRUE;
      else
	handle->ifdefs[handle->nifdef-1] = eifFALSE;
      sfree(name);
    }
    /* Don't print lines with ifdef or ifndef, go on to the next */
    return cpp_read_line(handlep,n,buf);
  }
  
  /* #else statement */
  if (strstrw(buf,"#else") != NULL) {
    if (handle->nifdef <= 0)
      return eCPP_SYNTAX;
    if (handle->ifdefs[handle->nifdef-1] == eifTRUE)
      handle->ifdefs[handle->nifdef-1] = eifFALSE;
    else if (handle->ifdefs[handle->nifdef-1] == eifFALSE)
      handle->ifdefs[handle->nifdef-1] = eifTRUE;
    
    /* Don't print lines with else, go on to the next */
    return cpp_read_line(handlep,n,buf);
  }
  
  /* #endif statement */
  if (strstrw(buf,"#endif") != NULL) {
    if (handle->nifdef <= 0)
      return eCPP_SYNTAX;
    handle->nifdef--;
    
    /* Don't print lines with endif, go on to the next */
    return cpp_read_line(handlep,n,buf);
  }

  /* Check whether we're not ifdeffed out. The order of this statement
     is important. It has to come after #ifdef, #else and #endif, but
     anything else should be ignored. */
  if ((handle->nifdef > 0) && (handle->ifdefs[handle->nifdef-1] != eifTRUE)) {
    return cpp_read_line(handlep,n,buf);
  }
  
  /* Check for include statements */
  if (strstrw(buf,"#include") != NULL) {
    len = -1;
    i0  = 0;
    for(i=0; (i<strlen(buf)); i++) {
      if ((buf[i] == '"') || (buf[i] == '<') || (buf[i] == '>'))  {
	if (len == -1) {
	  i0 = i+1;
	  len = 0;
	}
	else
	  break;
      }
      else if (len >= 0)
	len++;
    }
    snew(inc_fn,len+1);
    strncpy(inc_fn,buf+i0,len);
    inc_fn[len] = '\0';
    if (debug)
      fprintf(debug,"Going to open include file '%s' i0 = %d, strlen = %d\n",
	      inc_fn,i0,len);
    /* Open include file and store it as a child in the handle structure */
    status = cpp_open_file(inc_fn,(void *)&(handle->child),NULL);
    sfree(inc_fn);
    if (status != eCPP_OK) {
      handle->child = NULL;
      return status;
    }
    /* Make a linked list of open files and move on to the include file */
    handle->child->parent = handle;
    *handlep = handle->child;
    handle = *handlep;
    /* Don't print lines with include, go on to the next */
    return cpp_read_line(handlep,n,buf);
  }
  
  /* #define statement */
  if (strstrw(buf,"#define") != NULL) {
    add_define(buf+8);
  
    return cpp_read_line(handlep,n,buf);
  }
  
  /* #undef statement */
  if (strstrw(buf,"#undef") != NULL) {
    snew(name,strlen(buf));
    status = sscanf(buf,"%*s %s",name);
    for(i=0; (i<ndef); i++) {
      if (strcmp(defs[i].name,name) == 0) {
	sfree(defs[i].name);
	sfree(defs[i].def);
	break;
      }
    }
    sfree(name);
    for( ; (i<ndef-1); i++) {
      defs[i].name = defs[i+1].name;
      defs[i].def  = defs[i+1].def;
    }
    ndef--;
    
    /* Don't print lines with undef, go on to the next */
    return cpp_read_line(handlep,n,buf);
  }
  
  /* Check whether we have any defines that need to be replaced. Note
     that we have to use a best fit algorithm, rather than first come
     first go. We do this by sorting the defines on length first, and
     then on alphabetical order. */
  for(i=0; (i<ndef); i++) {
    if (defs[i].def) {
      nn  = 0;
      ptr = buf;
      while ((ptr = strstrw(ptr,defs[i].name)) != NULL) {
	nn++;
	ptr += strlen(defs[i].name);
      }
      if (nn > 0) {
	len = strlen(buf) + nn*max(4,4+strlen(defs[i].def)-strlen(defs[i].name));
	snew(name,len);
	ptr = buf;
	while ((ptr2 = strstrw(ptr,defs[i].name)) != NULL) {
	  strncat(name,ptr,(int)(ptr2-ptr));
	  strcat(name,defs[i].def);
	  ptr = ptr2 + strlen(defs[i].name);
	}
	strcat(name,ptr);
	strcpy(buf,name);
	sfree(name);
      }
    }
  }
  
  return eCPP_OK;
}

/* Close the file! Return integer status. */
int cpp_close_file(gmx_cpp_t *handlep)
{
  int i;
  gmx_cpp_t handle = (gmx_cpp_t)*handlep;
  
  if (!handle)
    return eCPP_INVALID_HANDLE;
  if (!handle->fp)
    return eCPP_FILE_NOT_OPEN;
  if (debug)
    fprintf(debug,"Closing file %s\n",handle->fn);
  fclose(handle->fp);
  
  if (0)switch(errno) {
  case 0:
    break;
  case ENOENT:
    return eCPP_FILE_NOT_FOUND;
  case EBADF:
    return eCPP_FILE_NOT_OPEN;
  case EINTR:
    return eCPP_INTERRUPT;
  default:
    if (debug)
      fprintf(debug,"Strange stuff closing file, errno = %d",errno);
    return eCPP_UNKNOWN;
  }
  handle->fp = NULL;
  handle->line_nr = 0;
  if (handle->fn) {
    sfree(handle->fn);
    handle->fn = NULL;
  }
  if (handle->line) {
    sfree(handle->line);
    handle->line = NULL;
  }
  if (handle->ifdefs) 
    sfree(handle->ifdefs);
  handle->nifdef = 0;
  
  return eCPP_OK;
}

/* Return a string containing the error message coresponding to status
   variable */
char *cpp_error(gmx_cpp_t *handlep,int status)
{
  char buf[256];
  char *ecpp[] = {
    "OK", "File not found", "End of file", "Syntax error", "Interrupted",
    "Invalid file handle", 
    "File not open", "Unknown error", "Error status out of range"
  };
  gmx_cpp_t handle = (gmx_cpp_t)*handlep;
  
  if (!handle)
    return ecpp[eCPP_INVALID_HANDLE];
    
  if ((status < 0) || (status >= eCPP_NR))
    status = eCPP_NR;
  
  sprintf(buf,"%s - File %s, line %d\nLast line read:\n'%s'",
	  ecpp[status],
	  (handle && handle->fn) ? handle->fn : "unknown",
	  (handle) ? handle->line_nr : -1,
	  handle->line ? handle->line : "");
  
  return strdup(buf);
}
