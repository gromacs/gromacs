/*
 *        @(#) copyrgt.c 1.19 10/17/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Good ROcking Metal Altar for Chronical Sinners
 */
#include "stdio.h"
#include "stdlib.h"
#include "macros.h"
#include "string2.h"
#include "futil.h"
#include "copyrite.h"

void head(FILE *out, char *fn)
{
  static char *head1[]= {
    " *",
    " *       This source code is part of",
    " *",
    " *        G   R   O   M   A   C   S",
    " *",
    " * GROningen MAchine for Chemical Simulations",
    " *"
  };
  static char *head2[] = {
    " * Please refer to:",
    " * GROMACS: A message-passing parallel molecular dynamics implementation",
    " * H.J.C. Berendsen, D. van der Spoel and R. van Drunen",
    " * Comp. Phys. Comm. 91, 43-56 (1995)",
    " *",
    " * Also check out our WWW page:",
    " * http://rugmd0.chem.rug.nl/~gmx",
    " * or e-mail to:",
    " * gromacs@chem.rug.nl",
    " *",
    " * And Hey:"
  };
#define NH1 asize(head1)
#define NCR asize(CopyrightText)
#define NH2 asize(head2)
  int i

  fprintf(out,"/*\n");
  fprintf(out," *       $id$\n");
  for(i=0; (i<NH1); i++)
    fprintf(out,"%s\n",head1[i]);
  fprintf(out," *            %s\n",GromacsVersion());
  for(i=0; (i<NCR); i++)
    fprintf(out," * %s\n",CopyrightText[i]);
  for(i=0; (i<NH2); i++)
    fprintf(out,"%s\n",head2[i]);

  fprintf(out," * %s\n */\n",bromacs());
  fprintf(out,"static *char SRCID = \"$Id$\"\n");
}

void cr(char *fn)
{
  FILE *in,*out;
#define MAXS 10240
  char buf[1024],line[MAXS+1];
  bool bSet;
  
  sprintf(buf,"%s.bak",fn);
  
  if (rename(fn,buf) != 0) {
    perror(buf);
    exit(1);
  }
  in=ffopen(buf,"r");
  out=ffopen(fn,"w");
  bSet=FALSE;
  
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
      if ((strstr(line,"*/") != NULL) && (!feof(in)))
	fgets2(line,MAXS,in);
    }
    head(out,fn);
    do {
      fprintf(out,"%s\n",line);
    } while (!feof(in) && fgets2(line,MAXS,in));
  }
  fclose(in);
  fclose(out);
}

int main(int argc,char *argv[])
{
  int i;
  
  for(i=1; (i<argc); i++) {
    fprintf(stderr,"Processing %s (backed up to %s.bak)\n",
	    argv[i],argv[i]);
    cr(argv[i]);
  }
  return 0;
}
