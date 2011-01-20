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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "string2.h"
#include "futil.h"
#include "smalloc.h"
#include "gmx_fatal.h"
#include "strdb.h"

gmx_bool get_a_line(FILE *fp,char line[],int n)
{
  char *line0;
  char *dum;
  
  snew(line0,n+1);
 
  do {
    if (!fgets(line0,n+1,fp)) {
      sfree(line0);
      return FALSE;
    }
    dum=strchr(line0,'\n');
    if (dum) 
      dum[0]='\0';
    else if (strlen(line0)==n) {
      fprintf(stderr,"Warning: line length exceeds buffer length (%d), data might be corrupted\n",n);
      line0[n-1] ='\0';
    } else
      fprintf(stderr,"Warning: file does not end with a newline, last line:\n%s\n",
	      line0);
    dum=strchr(line0,';');
    if (dum) 
      dum[0]='\0';
    strncpy(line,line0,n);
    dum=line0;
    ltrim(dum);
  } while (dum[0] == '\0'); 
  
  sfree(line0);
  return TRUE;
}

gmx_bool get_header(char line[],char *header)
{
  char temp[STRLEN],*dum;

  strcpy(temp,line);
  dum=strchr(temp,'[');
  if (dum==NULL)
    return FALSE;
  dum[0]=' ';
  dum=strchr(temp,']');
  if (dum==NULL) {
    gmx_fatal(FARGS,"header is not terminated on line:\n'%s'\n",line); 
    return FALSE;
  }
  dum[0]='\0';
  if (sscanf(temp,"%s%*s",header) != 1)
    return FALSE;

  return TRUE;
}

int get_strings(const char *db,char ***strings)
{
  FILE *in;
  char **ptr;
  char buf[256];
  int  i,nstr;

  in=libopen(db);
  
  if (fscanf(in,"%d",&nstr) != 1) {
    gmx_warning("File %s is empty",db);
    ffclose(in);
    return 0;
  }
  snew(ptr,nstr);
  for(i=0; (i<nstr); i++) {
    if(1 != fscanf(in,"%s",buf))
    { 
      gmx_fatal(FARGS,"Cannot read string from buffer");
    }
#ifdef DEBUG
    fprintf(stderr,"Have read: %s\n",buf);
#endif
    ptr[i] = strdup(buf);
  }
  ffclose(in);

  *strings=ptr;
  
  return nstr;
}

int search_str(int nstr,char **str,char *key)
{
  int i;

  /* Linear search */
  for(i=0; (i<nstr); i++)
    if (gmx_strcasecmp(str[i],key)==0)
      return i;

  return -1;
}

int fget_lines(FILE *in,char ***strings)
{
  char **ptr;
  char buf[256];
  int  i,nstr;
  char *pret;

  pret = fgets(buf,255,in);  
  if ( pret==NULL  || sscanf(buf,"%d",&nstr) != 1) 
  {
    gmx_warning("File is empty");
    ffclose(in);
    
    return 0;
  }
  snew(ptr,nstr);
  for(i=0; (i<nstr); i++) {
    fgets2(buf,255,in);
    ptr[i] = strdup(buf);
  }
  
  (*strings) = ptr;
  
  return nstr;
}

int get_lines(const char *db,char ***strings)
{
  FILE *in;
  int  nstr;
  
  in   = libopen(db);
  nstr = fget_lines(in,strings);
  ffclose(in);

  return nstr;
}

int get_file(const char *db,char ***strings)
{
  FILE *in;
  char **ptr=NULL;
  char buf[STRLEN];
  int  i,nstr,maxi;

  in=libopen(db);
  
  i=maxi=0;
  while (fgets2(buf,STRLEN-1,in)) {
    if (i>=maxi) {
      maxi+=50;
      srenew(ptr,maxi);
    }
    ptr[i] = strdup(buf);
    i++;
  }
  nstr=i;
  ffclose(in);
  srenew(ptr,nstr);
  *strings=ptr;
  
  return nstr;
}

