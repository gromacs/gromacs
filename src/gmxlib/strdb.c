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
static char *SRCID_strdb_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include "string2.h"
#include "futil.h"
#include "smalloc.h"
#include "fatal.h"
#include "strdb.h"

bool get_a_line(FILE *fp,char line[],int n)
{
  static char *line0=NULL;
  static int nalloc=0;
  char *dum;
  
  if (n>nalloc) {
    nalloc=n;
    srenew(line0,nalloc+1);
  }
  do {
    if (!fgets(line0,n+1,fp)) {
      return FALSE;
    }
    dum=strchr(line0,'\n');
    if (dum) 
      dum[0]='\0';
    else if (strlen(line0)==n) {
      fprintf(stderr,"Warning: line length exceeds buffer length (%d), data might be corrupted\n",n);
      line0[n-1] ='\0';
    } else
      fprintf(stderr,"Warning: file does end with a newline, last line:\n%s\n",
	      line0);
    dum=strchr(line0,';');
    if (dum) 
      dum[0]='\0';
    strcpy(line,line0);
    dum=line0;
    ltrim(dum);
  } while (dum[0] == '\0'); 

  return TRUE;
}

bool get_header(char line[],char *header)
{
  char temp[STRLEN],*dum;

  strcpy(temp,line);
  dum=strchr(temp,'[');
  if (dum==NULL)
    return FALSE;
  dum[0]=' ';
  dum=strchr(temp,']');
  if (dum==NULL) {
    fatal_error(0,"header is not terminated on line:\n'%s'\n",line); 
    return FALSE;
  }
  dum[0]='\0';
  if (sscanf(temp,"%s%*s",header) != 1)
    return FALSE;

  return TRUE;
}

int get_strings(char *db,char ***strings)
{
  FILE *in;
  char **ptr;
  char buf[256];
  int  i,nstr;

  in=libopen(db);
  
  set_warning_line(db,1);
  if (fscanf(in,"%d",&nstr) != 1) {
    sprintf(warn_buf,"File %s is empty",db);
    warning(NULL);
    fclose(in);
    return 0;
  }
  snew(ptr,nstr);
  for(i=0; (i<nstr); i++) {
    fscanf(in,"%s",buf);
#ifdef DEBUG
    fprintf(stderr,"Have read: %s\n",buf);
#endif
    ptr[i] = strdup(buf);
  }
  fclose(in);

  *strings=ptr;
  
  return nstr;
}

int search_str(int nstr,char **str,char *key)
{
  int i;

  /* Linear search */
  for(i=0; (i<nstr); i++)
    if (strcasecmp(str[i],key)==0)
      return i;

  return -1;
}

int fget_lines(FILE *in,char ***strings)
{
  char **ptr;
  char buf[256];
  int  i,nstr;
  
  fgets(buf,255,in);  
  if (sscanf(buf,"%d",&nstr) != 1) {
    sprintf(warn_buf,"File is empty");
    warning(NULL);
    fclose(in);
    
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

int get_lines(char *db,char ***strings)
{
  FILE *in;
  int  nstr;
  
  set_warning_line(db,1);
  in   = libopen(db);
  nstr = fget_lines(in,strings);
  fclose(in);

  return nstr;
}

int get_file(char *db,char ***strings)
{
  FILE *in;
  char **ptr=NULL;
  char buf[256];
  int  i,nstr,maxi;

  in=libopen(db);
  
  i=maxi=0;
  while (fgets2(buf,255,in)) {
    if (i>=maxi) {
      maxi+=50;
      srenew(ptr,maxi);
    }
    ptr[i] = strdup(buf);
    i++;
  }
  nstr=i;
  fclose(in);
  srenew(ptr,nstr);
  *strings=ptr;
  
  return nstr;
}

