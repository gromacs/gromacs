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
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include "string2.h"
#include "smalloc.h"
#include "macros.h"
#include "replace.h"

char *replace(const char *string,const char *search,const char *replace)
{
  char *buf=NULL,*ptr=NULL,*bufptr=NULL;
  int  blen,stringlen,slen,rlen;
  int  i,j,tmp;
  
  slen=strlen(search);
  stringlen=strlen(string);
  if ((string == NULL) || (slen == 0) || (stringlen == 0)) {
    if (string)
      buf=gmx_strdup(string);
    return buf;
  }
  rlen=strlen(replace);
  blen=max(stringlen,(rlen*stringlen)/slen);
  snew(buf,blen+1);
  strcpy(buf,string);
  
  bufptr=buf;
  while ((ptr=strstr(bufptr,search)) != NULL) {
    if (rlen <= slen) {
      for(i=0; (i<rlen); i++)
	ptr[i]=replace[i];
      if (rlen < slen) {
	while (ptr[i+slen-rlen] != '\0') {
	  ptr[i]=ptr[i+slen-rlen];
	  i++;
	}
	ptr[i]='\0';
      }
    }
    else {
      tmp=strlen(ptr);
      for(j=tmp; (j>=slen); j--)
	ptr[rlen-slen+j]=ptr[j];
      for(i=0; (i<rlen); i++)
	ptr[i]=replace[i];
    }
    bufptr=ptr+rlen;
  }
  
  return buf;
}

char *replaceww(const char *string,const char *search,const char *replace)
{
  char *buf=NULL,*ptr=NULL,*bufptr=NULL;
  int  buflen,stringlen,searchlen,replacelen;
  int  i,j;
  
  searchlen=strlen(search);
  stringlen=strlen(string);
  if ((string == NULL) || (searchlen == 0) || (stringlen == 0)) {
    if (string)
      buf=gmx_strdup(string);
    return buf;
  }  
  replacelen=strlen(replace);
  buflen=max(stringlen,(replacelen*stringlen)/searchlen);
  snew(buf,buflen+1);
  strcpy(buf,string);
  
  bufptr=buf;
  while ((ptr=strstr(bufptr,search)) != NULL) {
    if (((ptr==bufptr) || !isalnum(ptr[-1])) && !isalnum(ptr[searchlen])) {
      if (replacelen <= searchlen) {
	for(i=0; (i<replacelen); i++)
	  ptr[i]=replace[i];
	if (replacelen < searchlen) {
	  while (ptr[i+searchlen-replacelen] != '\0') {
	    ptr[i]=ptr[i+searchlen-replacelen];
	    i++;
	  }
	  ptr[i]='\0';
	}
      }
      else {
	for(j=strlen(ptr); (j>=searchlen); j--)
	  ptr[replacelen-searchlen+j]=ptr[j];
	for(i=0; (i<replacelen); i++)
	  ptr[i]=replace[i];
      }
      bufptr=ptr+replacelen;
    } else 
      bufptr=ptr+searchlen;
  }
  
  return buf;
}

#ifdef DEBUGREPLACE
void main(int argc,char *argv[])
{
  printf("String was: '%s' Search: '%s' Replace: '%s'\n",
	 argv[1],argv[2],argv[3]);
  printf("String now: '%s'\n\n",replace(argv[1],argv[2],argv[3]));
}
#endif
