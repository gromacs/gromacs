/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * S  C  A  M  O  R  G
 */
static char *SRCID_sheader_c = "$Id$";
#include <string.h>
#include "sheader.h"
#include "txtdump.h"
#include "tpxio.h"

static char *version="@(#) sheader.c 1.5 12/16/92";

void pr_header(FILE *fp,int indent,char *title,t_tpxheader *sh)
{
  if (available(fp,sh,title))
    {
      indent=pr_title(fp,indent,title);
      pr_indent(fp,indent);
      fprintf(fp,"bIr    = %spresent\n",sh->bIr?"":"not ");
      pr_indent(fp,indent);
      fprintf(fp,"bBox   = %spresent\n",sh->bBox?"":"not ");
      pr_indent(fp,indent);
      fprintf(fp,"bTop   = %spresent\n",sh->bTop?"":"not ");
      pr_indent(fp,indent);
      fprintf(fp,"bX     = %spresent\n",sh->bX?"":"not ");
      pr_indent(fp,indent);
      fprintf(fp,"bV     = %spresent\n",sh->bV?"":"not ");
      pr_indent(fp,indent);
      fprintf(fp,"bF     = %spresent\n",sh->bF?"":"not ");
      
      pr_indent(fp,indent);
      fprintf(fp,"natoms = %d\n",sh->natoms);
      pr_indent(fp,indent);
      fprintf(fp,"step   = %d\n",sh->step);
      pr_indent(fp,indent);
      fprintf(fp,"t      = %e\n",sh->t);
      pr_indent(fp,indent);
      fprintf(fp,"lambda = %e\n",sh->lambda);
    }
}
