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
static char *SRCID_prfn_c = "$Id$";

#include "filenm.h"
#include "futil.h"
#include "wman.h"

void pr_texdefs(FILE *fp)
{
  int i;

  fprintf(fp,"\\begin{table}\n");
  fprintf(fp,"\\begin{tabularx}{\\linewidth}{|r@{\\tt.}lccX|}\n");
  fprintf(fp,"\\dline\n");
  fprintf(fp,"\\mc{2}{|c}{%s} & %4s & %7s & %s \\\\[-0.1ex]\n",
	  "Default","","Default","");
  fprintf(fp,"\\mc{1}{|c}{%s} & \\mc{1}{c}{%s} & %4s & %7s & %s "
	  "\\\\[-0.1ex]\n",
	  "Name","Ext.","Type","Option","Description");
  fprintf(fp,"\\hline\n");
  for(i=0; (i<efNR); i++)
    if ( (i!=efGCT) && (i!=efHAT) )
      pr_def(fp,i);
  fprintf(fp,"\\dline\n");
  fprintf(fp,"\\end{tabularx}\n");
  fprintf(fp,"\\caption{The {\\gromacs} file types.}\n");
  fprintf(fp,"\\label{tab:form}\n");
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
 
