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
static char *SRCID_prfn_c = "$Id$";

#include "filenm.h"
#include "futil.h"
#include "wman.h"

void pr_texdefs(FILE *fp)
{
  int i;

  fprintf(fp,"\\begin{table}[p]\n");
  fprintf(fp,"\\begin{tabularx}{\\linewidth}{rlccX}\n");
  fprintf(fp,"\\hline\n");
  fprintf(fp,"%8s & %5s & %4s & %5s & %s\n",
	  "Def.","","","Def.","\\\\");
  fprintf(fp,"%8s & %5s & %4s & %5s & %s\n",
	  "Name","Ext.","Type","Opt.","Description \\\\");
  fprintf(fp,"\\hline\n");
  for(i=0; (i<efNR); i++)
    pr_def(fp,i);
  fprintf(fp,"\\hline\n");
  fprintf(fp,"\\end{tabularx}\n");
  fprintf(fp,"\\caption{File types: A = ascii, B = binary, P = XDR portable.}\n");
  fprintf(fp,"\\label{Tab:form}\n");
  fprintf(fp,"\\end{table}\n");
}

void pr_htmldefs(FILE *fp)
{
  int i;

  fprintf(fp,"<title>GROMACS</title>\n");
  fprintf(fp,"<h1>GROMACS Files</h1>\n");
  fprintf(fp,"<b>GRO</b>ningen <b>MA</b>chine for <b>S</b>imulating <b>C</b>hemistry\n");
  fprintf(fp,"<p>\n");
  fprintf(fp,"The following %d filetypes are used by Gromacs:\n",efNR);
  fprintf(fp,"<dl>\n");
  for(i=0; (i<efNR); i++) {
    fprintf(fp,"<dt><a href=\"%s.html\">%s.%s</a> (%s)<dd>%s\n",
	    ftp2ext(i),ftp2defnm(i),ftp2ext(i),ftp2ftype(i),
	    check_html(ftp2desc(i),NULL));
  }
  fprintf(fp,"</dl>\n");
}


void main()
{
  FILE *out;
  
  out=ffopen("files.tex","w");
  pr_texdefs(out);
  fclose(out);
  
  out=ffopen("files.html","w");
  pr_htmldefs(out);
  fclose(out);
}
 
