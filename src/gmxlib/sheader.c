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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#include <string.h>
#include "symtab.h"
#include "magic.h"
#include "fatal.h"
#include "sheader.h"
#include "txtdump.h"
#include "copyrite.h"
#include "binio.h"
#include "recio.h"

#define BUFSIZE		128
#define	SKIP_SCCS	5

static char *version="@(#) sheader.c 1.5 12/16/92";

void wr_header(FILE *fp,t_statheader *sh)
{
  static int magic=GROMACS_MAGIC;
  int len;
  
  len=strlen(version)+1;
  blockwrite(fp,magic);
  blockwrite(fp,len);
  cblockwrite(fp,version,len);
  header_blockio(fp,write,*sh);
}

bool is_status(FILE *fp)
{
  static int magic=0;

  blockread(fp,magic);
  rewind(fp);
  
  if (magic == GROMACS_MAGIC)
    return TRUE;
  else
    return FALSE;
}

char *rd_header(FILE *fp,t_statheader *sh)
{
  static char header_version[BUFSIZE];
  static int magic=0;
  int len;

  blockread(fp,magic);
  blockread(fp,len);
  cblockread(fp,header_version,len);
  header_blockio(fp,read,*sh);
  if (magic!=GROMACS_MAGIC) {
    fprintf(stderr,"ERROR: incorrect tpb or trj file version.\n");
    fprintf(stderr,"       make sure you use GROMACS %s throughout\n",
	    GromacsVersion());
    fprintf(stderr,"       for trj files you can use trjconv -patch to update the version number\n");
    exit(1);
  }
  return header_version;
}

void pr_header(FILE *fp,int indent,char *title,t_statheader *sh)
{
  if (available(fp,sh,title))
    {
      indent=pr_title(fp,indent,title);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"ir_size     = %d\n",sh->ir_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"e_size      = %d\n",sh->e_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"box_size    = %d\n",sh->box_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"vir_size    = %d\n",sh->vir_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"pres_size   = %d\n",sh->pres_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"top_size    = %d\n",sh->top_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"sym_size    = %d\n",sh->sym_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"x_size      = %d\n",sh->x_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"v_size      = %d\n",sh->v_size);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"f_size      = %d\n",sh->f_size);

      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"natoms      = %d\n",sh->natoms);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"step        = %d\n",sh->step);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"nre         = %d\n",sh->nre);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"t           = %e\n",sh->t);
      (void) pr_indent(fp,indent);
      (void) fprintf(fp,"lambda      = %e\n",sh->lambda);
    }
}
