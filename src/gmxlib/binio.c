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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
static char *SRCID_binio_c = "$Id$";

#include <stdio.h>
#include <errno.h>
#include "sysstuff.h"
#include "main.h"
#include "fatal.h"
#include "binio.h"

void _blockwrite(FILE *fp,int nelem,int size,void *data,char *what,char *file,
                 int line)
{
  int len;
  
#ifdef DEBUG
  (void) fprintf(stderr,"blockwrite %s (file %s,line %d,size=%d, nelem=%d)\n",
                 what,file,line,size,nelem);
  fflush(stderr);
#endif
  if ((len=fwrite(data,size,nelem,fp))!=nelem)
    fatal_error(errno,"writing %s (%dx%d, len=%d, pos=%d) from file %s, "
                "line %d",what,nelem,size,len,ftell(fp),file,line);
}

void _blockread(FILE *fp,int nelem,int size,void *data,char *what,char *file,
                int line)
{
  int len;
  
#ifdef DEBUG
  (void) fprintf(stderr,"blockread %s (file %s,line %d,size=%d, nelem=%d)\n",
                 what,file,line,size,nelem);
  fflush(stderr);
#endif
  if ((len=fread(data,size,nelem,fp))!=nelem)
    fatal_error(errno,"reading %s (%dx%d, len=%d, pos=%d) from file %s, "
                "line %d",what,nelem,size,len,ftell(fp),file,line);
}
