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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_strdb_c = "$Id$";

#include <stdio.h>
#include <stdlib.h>
#include "string2.h"
#include "futil.h"
#include "smalloc.h"
#include "fatal.h"
#include "strdb.h"

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

