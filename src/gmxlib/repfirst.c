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
static char *SRCID_repfirst_c = "$Id$";

#include "string2.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "replace.h"
#include "strdb.h"

bool replace1(char *string,char *buf,char *search,char *replace)
{
  char *ptr,*bufptr;
  int  blen,stringlen,slen,rlen;
  int  i,j,tmp;
  
  slen=strlen(search);
  stringlen=strlen(string);
  if ((string == NULL) || (slen == 0) || (stringlen == 0))
    return FALSE;
  
  rlen=strlen(replace);
  sprintf(buf,"%s",string);
  
  bufptr=buf;
  if ((ptr=strstr(bufptr,search)) != NULL) {
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
    
    return TRUE;
  }
  
  return FALSE;
}

void main(int argc,char *argv[])
{
#define MYBUFS 10240
  int  i,nstr;
  char **str;
  char **rep;
  bool *bRep;
  char buffer[MYBUFS];
  char newbuf[MYBUFS];
  
  nstr=get_strings("index.gmx",&str);
  snew(bRep,nstr);
  snew(rep,nstr);
  for(i=0; (i<nstr); i++) {
    sprintf(buffer,"\\myindex{%s}",str[i]);
    rep[i]=strdup(buffer);
  }
  
  while(fgets2(buffer,MYBUFS-1,stdin) != NULL) {
    newbuf[0]='\0';
    for(i=0; (i<nstr); i++)
      if (!bRep[i]) {
	bRep[i]=replace1(buffer,newbuf,str[i],rep[i]);
	if (bRep[i])
	  strcpy(buffer,newbuf);
      }
    printf("%s\n",newbuf);
  }
}
