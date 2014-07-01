/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h> 
#include <ctype.h>
#include "gromacs/utility/smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/math/vec.h"
#include "gromacs/commandline/pargs.h"
#include "copyrite.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/strdb.h"

gmx_bool isword(char c)
{
  return (isalnum(c) || (c=='-') || (c=='_'));
}

char *strncasestr(char *line,char *str)
{
  char *dum;

  dum=line;
  if (dum) {
    while (strlen(dum) && strncasecmp(dum,str,strlen(str)))
      dum++;
    if (strlen(dum)==0)
      dum=NULL;
  }

  return dum;
}

char *strstr_href(char *line,gmx_bool *bInHREF,int *i_dat,int n_dat,char **dat)
{
  char *start,*found,*href=NULL;
  gmx_bool bIn;
  int i;

  found=NULL;
  *i_dat=-1;
  bIn=*bInHREF;
  start=line;
  do {
    if (bIn) {
      while (strlen(start) && (strncasecmp(start,"</a",3) != 0))
	start++;
      if (strlen(start)>0) {
	start+=3;
	bIn=FALSE;
      }
    }
    else {
      href=strncasestr(start,"<a href");
      if (href)
	bIn=TRUE;
      i=0;
      while((i<n_dat) && !found) {
	found=strncasestr(start,dat[i]);
	if (found) {
	  if (href && (found>href))
	    found=NULL;
	  else {
	    if (((found!=start) && isword(found[-1])) || 
		isword(found[strlen(dat[i])])) 
	      found=NULL;
	    else
	      *i_dat=i;
	  }
	  i++;
	}
      }
    }
  } while (strlen(start) && !found && href);
  *bInHREF=bIn;

  return found;
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "[TT]hrefify[tt] adds href's for all the words in the input file which are not",
    "already hyperlinked and which appear in the file specified with the",
    "option [TT]-l[tt].[PAR]",
    "If the href's should call a script, text can be added",
    "with the [TT]-t[tt] option."
  };

  int n;
  
  char **text,**str,line[1024],*ptr,*ref,
    start[STRLEN],word[STRLEN],end[STRLEN];
  int n_text,n_str,i_str;
  gmx_bool bInHREF,bIn;
  
  FILE    *fp;
  char    title[STRLEN];
  int     i,l,n_repl;
  t_filenm fnm[] = {
    { efDAT, "-l", "links", ffLIBRD },
  };
#define NFILE asize(fnm)
  static char *in=NULL,*out=NULL,*excl=NULL,*link_text=NULL;
  static gmx_bool peratom=FALSE;
  t_pargs pa[] = {
    { "-f", FALSE, etSTR, { &in } , "HTML input" },
    { "-o", FALSE, etSTR, { &out } , "HTML output" },
    { "-e", FALSE, etSTR, { &excl } , "Exclude a string from HREF's, "
      "when this option is not set, the filename without path and extension "
      "will be excluded from HREF's"},
    { "-t", FALSE, etSTR, { &link_text } , "Insert a string in front of the "
      "href file name, useful for scripts" }
  };
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);

  if (!in || !out)
    gmx_fatal(FARGS,"Input or output filename is not set");

  n_text = get_file(in, &text);
  fprintf(stderr,"Read %d lines from %s\n",n_text,in);

  n_str=get_file(ftp2fn(efDAT,NFILE,fnm),&str);  
  fprintf(stderr,"Read %d strings %s\n",n_str,ftp2fn(efDAT,NFILE,fnm));
  if (!excl) {
    for (i=strlen(in)-1; i>0 && in[i-1]!='/'; i--);
    excl=strdup(in+i);
    for(i=strlen(excl)-1; i>0 && (excl[i]!='.'); i--);
    if (excl[i]=='.')
      excl[i]='\0';
  }
  fprintf(stderr,"Excluding '%s' from references\n",excl);
  for(l=0; l<n_str && strcasecmp(str[l],excl); l++);
  if (l<n_str) {
    for(i=l+1; i<n_str; i++)
      str[i-1]=str[i];
    n_str--;
  }

  if (!link_text)
    link_text=strdup("\0");
  else
    fprintf(stderr,"Adding '%s' to href's\n",link_text);

  fp=gmx_ffopen(out,"w");

  n_repl=0;
  i_str=-1;
  bInHREF=FALSE;
  for(l=0; l<n_text; l++) {
    strcpy(line,text[l]);
    do {
      bIn=bInHREF;
      ptr=strstr_href(line,&bIn,&i_str,n_str,str);
      if (ptr) {
	ref=ptr;
	if ((ref!=line) && (ref[-1]=='.')) {
	  ref--;
	  while((ref>line) && isword(ref[-1]))
	    ref--;
	}
	strcpy(start,line);
	start[ref-line]='\0';
	strcpy(word,ref);
        word[ptr-ref+strlen(str[i_str])]='\0';
	strcpy(end,ptr+strlen(str[i_str]));
	sprintf(line,"%s<a href=\"%s%s.html\">%s</a>%s",
		start,link_text,str[i_str],word,end);
	fprintf(stderr,"line %d: %s\n",l+1,str[i_str]);
	n_repl++;
      }
    } while (ptr);
    bInHREF=bIn;
    fprintf(fp,"%s\n",line);
  }

  gmx_ffclose(fp);
  
  fprintf(stderr,"Added %d HTML references\n",n_repl);

  return 0;
}

