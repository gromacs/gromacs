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
 * Gyas ROwers Mature At Cryogenic Speed
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
