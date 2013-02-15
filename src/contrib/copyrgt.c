/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "stdio.h"
#include "stdlib.h"
#include "macros.h"
#include "string2.h"
#include "futil.h"
#include "copyrite.h"

static char *head1[]= {
  "",
  "               This source code is part of",
  "",
  "                G   R   O   M   A   C   S",
  "",
  "         GROningen MAchine for Chemical Simulations",
  ""
};

static char *head2[]= {
  "This program is free software; you can redistribute it and/or",
  "modify it under the terms of the GNU General Public License",
  "as published by the Free Software Foundation; either version 2",
  "of the License, or (at your option) any later version.",
  "",
  "If you want to redistribute modifications, please consider that",
  "scientific software is very special. Version control is crucial -",
  "bugs must be traceable. We will be happy to consider code for",
  "inclusion in the official distribution, but derived work must not",
  "be called official GROMACS. Details are found in the README & COPYING",
  "files - if they are missing, get the official version at www.gromacs.org.",
  "",
  "To help us fund GROMACS development, we humbly ask that you cite",
  "the papers on the package - you can find them in the top README file.",
  "",
  "For more info, check our website at http://www.gromacs.org",
  "",
  "And Hey:"
};

#define NH1 asize(head1)
#define NCR asize(CopyrightText)
#define NH2 asize(head2)
#define MAXS 10240

void head(FILE *out, char *fn_, gmx_bool bH,
	  char *cstart, char *ccont, char *cend)
{
  char buf[STRLEN];
  int i;

  fprintf(out,"%s\n",cstart);
  /* NOTE: the "" are to mislead CVS so it will not replace by version info */
  fprintf(out,"%s $""Id""$\n",ccont);
  for(i=0; (i<NH1); i++)
    fprintf(out,"%s %s\n",ccont,head1[i]);
  fprintf(out,"%s                        %s\n",ccont,GromacsVersion());
  for(i=0; (i<NCR); i++)
    fprintf(out,"%s %s\n",ccont,CopyrightText[i]);
  for(i=0; (i<NH2); i++)
    fprintf(out,"%s %s\n",ccont,head2[i]);
  bromacs(buf,STRLEN-1);
  fprintf(out,"%s %s\n",ccont,buf);
  fprintf(out,"%s\n",cend);
  if (bH) {
    fprintf(out,"\n");
    fprintf(out,"#ifndef _%s\n",fn_);
    fprintf(out,"#define _%s\n",fn_);
    fprintf(out,"\n");
  }
}

void cr_c(char *fn)
{
  FILE *in,*out;
  char ofn[1024],line[MAXS+1],cwd[1024];
  char *p,*fn_;
  gmx_bool bH;
  
  sprintf(ofn,"%s.bak",fn);
  
  fprintf(stderr,"Processing %s (backed up to %s)\n",
	  fn,ofn);
  
  if (rename(fn,ofn) != 0) {
    perror(ofn);
    exit(1);
  }
  in=ffopen(ofn,"r");
  out=ffopen(fn,"w");
  
  /* Skip over empty lines in the beginning only */
  do { 
    if (fgets2(line,MAXS,in))
      rtrim(line); 
  } while ((strlen(line) == 0) && (!feof(in)));
  
  /* Now we are at end of file, or we have a non-empty string */
  if (strlen(line) != 0) {  
    if (strstr(line,"/*") != NULL) {
      /* File does start with comment, so delete it and add new */
      while ((strstr(line,"*/") == NULL) && (!feof(in)))
	fgets2(line,MAXS,in);
    }
    fn_=strdup(fn);
    p=strchr(fn_,'.');
    if (p)
      p[0]='_';
    bH=FALSE;
    do {
      fgets2(line,MAXS,in);
      if ( (strstr(line,fn_) != NULL) && 
	   (strstr(line,"#define") != NULL) )
	bH=TRUE;
    } while ( ( (strstr(line,fn_) != NULL)  ||
		(strlen(line)==0) ) && (!feof(in) ) );
    getcwd(cwd,STRLEN);
    head(out,fn_,bH,"/*"," *"," */");
    do {
      fprintf(out,"%s\n",line);
    } while (!feof(in) && fgets2(line,MAXS,in));
  }
  ffclose(in);
  ffclose(out);
}

void cr_other(char *fn)
{

  /* Doesnt work right now, so its commented out */
  /*  FILE *in,*out;
  char ofn[1024],line[MAXS+1],line2[MAXS+1],cwd[1024];
  char *p,*fn_,*ptr;
  gmx_bool bH;
  
  sprintf(ofn,"%s.bak",fn);
  
  fprintf(stderr,"Processing %s (backed up to %s)\n",fn,ofn);
  
  if (rename(fn,ofn) != 0) {
    perror(ofn);
    exit(1);
  }
  in=ffopen(ofn,"r");
  out=ffopen(fn,"w");
  */
  /* Skip over empty lines in the beginning only */

}

void cr_tex(char *fn)
{
  FILE *in,*out;
  char ofn[1024],line[MAXS+1];
  char *p;
  
  sprintf(ofn,"%s.bak",fn);
  
  fprintf(stderr,"Processing (as Tex) %s (backed up to %s)\n",
	  fn,ofn);
  
  if (rename(fn,ofn) != 0) {
    perror(ofn);
    exit(1);
  }
  in=ffopen(ofn,"r");
  out=ffopen(fn,"w");
  
  /* Skip over empty lines in the beginning only */
  do
    if (fgets2(line,MAXS,in))
      rtrim(line);
  while ((strlen(line) == 0) && (!feof(in)));
  
  /* Now we are at end of file, or we have a non-empty string */
  if (strlen(line) != 0) {  
    while ((strstr(line,"%") != NULL) && (!feof(in)))
      /* File does start with comment, so delete it and add new */
      fgets2(line,MAXS,in); 
    head(out,"",FALSE,"%","%","%");
    /* Skip over empty lines */
    while ( (strlen(line) == 0) && !feof(in) )
      if (fgets2(line,MAXS,in))
	rtrim(line);
    do
      fprintf(out,"%s\n",line);
    while (!feof(in) && fgets2(line,MAXS,in));
  }
  ffclose(in);
  ffclose(out);
}

int main(int argc,char *argv[])
{
  int i;
  char *fn,*p;
  
  for(i=1; (i<argc); i++) {
    fn=argv[i];
    p=strrchr(fn,'.');
    if ( strcmp(p,".tex")==0 )
      cr_tex(fn);
    else if ((strcmp(p,".c") == 0) || (strcmp(p,".h") == 0))
      cr_c(fn);
    else
      cr_other(fn);
  }
  return 0;
}
