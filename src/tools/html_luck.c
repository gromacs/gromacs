/*
 *       @(#) copyrgt.c 1.12 9/30/97
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
 * GROwing Monsters And Cloning Shrimps
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "copyrite.h"
#include "string2.h"

void catf(char *fn,FILE *fp)
{
  FILE *in;
  char buf[2048];
  
  if ((in=fopen(fn,"r")) == NULL) {
    fprintf(fp,
	    "Welcome to GROMACS. Sorry to say that something went wrong.\n");
  }
  else {
    while(fgets2(buf,2047,in) != NULL)
      fprintf(fp,"%s\n",buf);
  
    fclose(in);
  }
}

int counter(char *fn)
{
  FILE *fp;
  int  cnt;
  
  if ((fp=fopen(fn,"r")) == NULL)
    return 0;
  fscanf(fp,"%d",&cnt);
  fclose(fp);
  cnt++;
  fp=fopen(fn,"w");
  fprintf(fp,"%d\n",cnt);
  fclose(fp);
  
  return cnt;
}

int main () {
  int cnt;
  
  cnt=counter("/home3/gmx/public_html/gmx.cnt");
  printf("Content-type: text/html%c%c",10,10);
  printf("<title>You Are The %dth Reader of the GROMACS Homepage</title>\n",
	 cnt);
  printf("<body>");
  catf("/home3/gmx/public_html/gmx.html",stdout);
  printf("<h3>%s</h3>\n",cool_quote());
  printf("</body>\n");

  return 0;
}

