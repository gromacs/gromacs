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
 * Great Red Oystrich Makes All Chemists Sane
 */
static char *SRCID_copyrgt_c = "$Id$";

#include "stdio.h"
#include "stdlib.h"
#include "macros.h"
#include "string2.h"
#include "futil.h"
#include "copyrite.h"

static char *head1[]= {
  "",
  "      This source code is part of",
  "",
  "       G   R   O   M   A   C   S",
  "",
  "GROningen MAchine for Chemical Simulations",
  ""
};
static char *head2[] = {
  "Please refer to:",
  "GROMACS: A message-passing parallel molecular dynamics implementation",
  "H.J.C. Berendsen, D. van der Spoel and R. van Drunen",
  "Comp. Phys. Comm. 91, 43-56 (1995)",
  "",
  "Also check out our WWW page:",
  "http://md.chem.rug.nl/~gmx",
  "or e-mail to:",
  "gromacs@chem.rug.nl",
  "",
  "And Hey:"
};
#define NH1 asize(head1)
#define NCR asize(CopyrightText)
#define NH2 asize(head2)
#define MAXS 10240

void head(FILE *out, char *fn_, bool bH, bool bSRCID,
	  char *cstart, char *ccont, char *cend)
{
  int i;

  fprintf(out,"%s\n",cstart);
  /* NOTE: the "" are to mislead CVS so it will not replace by version info */
  fprintf(out,"%s $""Id""$\n",ccont);
  for(i=0; (i<NH1); i++)
    fprintf(out,"%s %s\n",ccont,head1[i]);
  fprintf(out,"%s               %s\n",ccont,GromacsVersion());
  for(i=0; (i<NCR); i++)
    fprintf(out,"%s %s\n",ccont,CopyrightText[i]);
  for(i=0; (i<NH2); i++)
    fprintf(out,"%s %s\n",ccont,head2[i]);

  fprintf(out,"%s %s\n",ccont,bromacs());
  fprintf(out,"%s\n",cend);
  if (bH) {
    fprintf(out,"\n");
    fprintf(out,"#ifndef _%s\n",fn_);
    fprintf(out,"#define _%s\n",fn_);
    fprintf(out,"\n");
  }
  if (bSRCID)
    fprintf(out,"static char *SRCID_%s = \"$""Id""$\";\n",fn_);
  /* NOTE: the "" are to mislead CVS so it will not replace by version info */
  fprintf(out,"\n");
}

void cr_c(char *fn)
{
  FILE *in,*out;
  char ofn[1024],line[MAXS+1],cwd[1024];
  char *p,*fn_;
  bool bH,bSRCID;
  
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
		(strstr(line,"static char *SRCID") != NULL) ||
		(strlen(line)==0) ) && (!feof(in) ) );
    getcwd(cwd,STRLEN);
    bSRCID = TRUE;
    /* Do not put source id's in include/types since some filenames are
     * be equal to those in include */
    if ((strlen(cwd)>strlen("types")) &&
	(strcmp(cwd+strlen(cwd)-strlen("types"),"types") == NULL))
      bSRCID = FALSE;
    head(out,fn_,bH,bSRCID,"/*"," *"," */");
    do {
      fprintf(out,"%s\n",line);
    } while (!feof(in) && fgets2(line,MAXS,in));
  }
  fclose(in);
  fclose(out);
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
    head(out,"",FALSE,FALSE,"%","%","%");
    /* Skip over empty lines */
    while ( (strlen(line) == 0) && !feof(in) )
      if (fgets2(line,MAXS,in))
	rtrim(line);
    do
      fprintf(out,"%s\n",line);
    while (!feof(in) && fgets2(line,MAXS,in));
  }
  fclose(in);
  fclose(out);
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
    else
      cr_c(fn);
  }
  return 0;
}
