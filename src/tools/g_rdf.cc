/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * GRoups of Organic Molecules in ACtion for Science
 */
static char *SRCID_g_rdf_cc = "$Id$";

#include <math.h>
#include <ctype.h>
#include "string2.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "pbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "rdgroup.h"
#include "do_rdf.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_rdf calculates radial distribution functions in different ways.",
    "The normal method is spherical around a (set of) particle(s), the",
    "other method is by dividing the sphere into segments (12 x 30 deg.)",
    "and plot the rdf as a function of both distance and angle with",
    "a given axis.[PAR]",
    "The option -c is meant to avoid intramolecular peaks in the rdf plot.",
    "It is however better to make a topology file with a higher number of",
    "exclusions. For eg. benzene a topology with nrexcl set to 5",
    "would eliminate all intramolecular contributions to the rdf." 
  };
  static int bAll=FALSE,bAng=FALSE,bCM=FALSE;
  static real cutoff=0;
  static int  axis=2;
  t_pargs pa[] = {
    { "-a", FALSE, etBOOL, &bAll,
      "calculates rdfs of all combinations of groups" },
    { "-cm",FALSE, etBOOL, &bCM,
      "compute RDF with respect to the center of mass of first group. This may be useful for e.g. micelles." },
    { "-c", FALSE, etREAL, &cutoff,
      "shortest distance that will be considered."}, 
    { "-d",  FALSE, etINT, &axis,
      "calculates angular rdfs using angle relative to axis (0=X,1=Y,2=Z) given on command line." }
  };
#define NPA asize(pa)
  char       *status;
  char       **grpname;
  char       outf1[STRLEN],outf2[STRLEN];
  int        gnx[2];
  atom_id    *index[2];
  int        i,j;
  t_topology *top;
  t_block    *excl,*block;
  t_filenm   fnm[] = {
    { efTRX, "-f", NULL,   ffREAD },
    { efTPB, NULL, NULL,   ffREAD },
    { efNDX, NULL, NULL,   ffREAD },
    { efXVG, NULL, NULL,   ffWRITE },
    { efXVG, "-oa","angular", ffOPTWR }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  bAng = (opt2parg_bSet("-d",NPA,pa));
  if (bCM && bAng)
    fatal_error(0,"Can't do angular and center of mass RDF at the same time!");
  if (bAng) {
    if ((axis < 0) || (axis > 2)) {
      fprintf(stderr,"Invalid axis %d must be 0(X), 1(Y), or 2(Z)\n",axis);
      exit(1);
    }
  }
   
  top=read_top(ftp2fn(efTPB,NFILE,fnm));
  mk_single_top(top);
  excl=&(top->atoms.excl);
  status=ftp2fn(efTRX,NFILE,fnm);

  if (bAll) {
    block=init_index(ftp2fn(efNDX,NFILE,fnm),&grpname);
    for(i=0; (i<block->nr); i++)
      for(j=i; (j<block->nr); j++) {
	sprintf(outf1,"%s-%s.%s",
		grpname[i],grpname[j],ftp2fn(efXVG,NFILE,fnm));
	sprintf(outf2,"%s-%s.%s",
		grpname[i],grpname[j],opt2fn("-oa",NFILE,fnm));
	fprintf(stderr,"RDF of %s-%s output in %s\n",
		grpname[i],grpname[j],outf1);
	if (bAng) {
	  do_rdfang(status,outf1,outf2,cutoff,excl,
		    &(block->a[block->index[i]]),
		    block->index[i+1]-block->index[i],grpname[i],
		    &(block->a[block->index[j]]),
		    block->index[j+1]-block->index[j],grpname[j],
		    axis);
	}
	else {
	  do_rdfn(status,outf1,cutoff,excl,bCM,
		  &(block->a[block->index[i]]),
		  block->index[i+1]-block->index[i],grpname[i],
		  &(block->a[block->index[j]]),
		  block->index[j+1]-block->index[j],grpname[j]);
	}
      }
  }
  else {
    grpname=(char **)calloc(2,sizeof(grpname[0]));
    rd_index(ftp2fn(efNDX,NFILE,fnm),2,gnx,index,grpname);
    
    if (bAng) {
      do_rdfang(status,ftp2fn(efXVG,NFILE,fnm),
		opt2fn("-oa",NFILE,fnm),cutoff,excl,
		index[0],gnx[0],grpname[0],
		index[1],gnx[1],grpname[1],
		axis);
      xvgr_file(opt2fn("-oa",NFILE,fnm),NULL);
    }
    else {
      do_rdfn(status,ftp2fn(efXVG,NFILE,fnm),cutoff,excl,bCM,
	      index[0],gnx[0],grpname[0],
	      index[1],gnx[1],grpname[1]);
    }
    xvgr_file(ftp2fn(efXVG,NFILE,fnm),NULL);
  }
    
  thanx(stdout);
  
  return 0;
}
