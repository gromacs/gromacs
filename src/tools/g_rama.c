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
 * Good ROcking Metal Altar for Chronical Sinners
 */
static char *SRCID_g_rama_c = "$Id$";

#include <math.h>
#include "sysstuff.h"
#include "string.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "xvgr.h"
#include "physics.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "gstat.h"

static void plot_rama(FILE *out,t_xrama *xr)
{
  int i;
  real phi,psi;
  
  for(i=0; (i<xr->npp); i++) {
    phi=xr->dih[xr->pp[i].iphi].ang*RAD2DEG;
    psi=xr->dih[xr->pp[i].ipsi].ang*RAD2DEG;
    fprintf(out,"%g  %g  %s\n",phi,psi,xr->pp[i].label);
  }
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_rama selects the Phi/Psi dihedral combinations from your topology file",
    "and computes these as a function of time.",
    "Using simple Unix tools such as [IT]grep[it] you can select out", 
    "specific residues."
  };

  FILE      *out;
  t_xrama   *xr;
  int       j;
  t_filenm  fnm[] = {
    { efTRX, "-f", NULL,  ffREAD },
    { efTPB, NULL, NULL,  ffREAD },
    { efXVG, NULL, NULL,  ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,0,NULL,asize(desc),desc,0,NULL);
  
  snew(xr,1);
  init_rama(ftp2fn(efTRX,NFILE,fnm),ftp2fn(efTPB,NFILE,fnm),xr);
  
  out=xvgropen(ftp2fn(efXVG,NFILE,fnm),"Ramachandran Plot","Phi","Psi");
  xvgr_line_props(out,0,elNone,ecFrank);
  xvgr_view(out,0.2,0.2,0.8,0.8);
  xvgr_world(out,-180,-180,180,180);
  fprintf(out,"@    xaxis  tick on\n@    xaxis  tick major 60\n@    xaxis  tick minor 30\n");
  fprintf(out,"@    yaxis  tick on\n@    yaxis  tick major 60\n@    yaxis  tick minor 30\n");
  fprintf(out,"@ s0 symbol 2\n@ s0 symbol size 0.4\n@ s0 symbol fill 1\n");
  
  j=0;
  do {
    if ((j%10) == 0)
      fprintf(stderr,"\rFrame: %d",j);
    plot_rama(out,xr);
    j++;
  } while (new_data(xr));
  fprintf(stderr,"\n");
  fclose(out);
  
  xvgr_file(ftp2fn(efXVG,NFILE,fnm),NULL);
  
  thanx(stdout);
  
  return 0;
}
