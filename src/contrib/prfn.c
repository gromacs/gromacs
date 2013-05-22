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


int main()
{
  FILE *out;
  
  out=ffopen("files.tex","w");
  pr_texdefs(out);
  ffclose(out);
  
  out=ffopen("files.html","w");
  pr_htmldefs(out);
  ffclose(out);

  return 0;
}
 
