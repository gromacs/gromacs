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
 * Grunge ROck MAChoS
 */
static char *SRCID_replace_c = "$Id$";

#include <ctype.h>
#include "string2.h"
#include "smalloc.h"
#include "macros.h"
#include "replace.h"

char *replace(char *string,char *search,char *replace)
{
  char *buf=NULL,*ptr=NULL,*bufptr=NULL;
  int  blen,stringlen,slen,rlen;
  int  i,j,tmp;
  
  slen=strlen(search);
  stringlen=strlen(string);
  if ((string == NULL) || (slen == 0) || (stringlen == 0))
    return string;
  
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

char *replaceww(char *string,char *search,char *replace)
{
  char *buf=NULL,*ptr=NULL,*bufptr=NULL;
  int  buflen,stringlen,searchlen,replacelen;
  int  i,j;
  
  searchlen=strlen(search);
  stringlen=strlen(string);
  if ((string == NULL) || (searchlen == 0) || (stringlen == 0))
    return string;
  
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

#ifdef DEBUG
void main(int argc,char *argv[])
{
  printf("String was: '%s' Search: '%s' Replace: '%s'\n",
	 argv[1],argv[2],argv[3]);
  printf("String now: '%s'\n\n",replace(argv[1],argv[2],argv[3]));
}
#endif
