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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_g_rdf_c = "$Id$";

#include <math.h>
#include <ctype.h>
#include "string2.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "rmpbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "tpxio.h"
#include "rdgroup.h"
#include "smalloc.h"

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "g_rdf calculates radial distribution functions in different ways.",
    "The normal method is around a (set of) particle(s), the other method",
    "is around the center of mass of a set of particles.[PAR]",
    "If a run input file is supplied ([TT]-s[tt]), exclusions defined",
    "in that file are taken into account when calculating the rdf.",
    "The option [TT]-cut[tt] is meant as an alternative way to avoid",
    "intramolecular peaks in the rdf plot.",
    "It is however better to supply a run input file with a higher number of",
    "exclusions. For eg. benzene a topology with nrexcl set to 5",
    "would eliminate all intramolecular contributions to the rdf.",
    "Note that all atoms in the selected groups are used, also the ones",
    "that don't have Lennard-Jones interactions." 
  };
  static bool bCM=FALSE;
  static real cutoff=0, binwidth=0.005;
  t_pargs pa[] = {
    { "-bin",      FALSE, etREAL, {&binwidth},
      "Binwidth (nm)" },
    { "-com",      FALSE, etBOOL, {&bCM},
      "RDF with respect to the center of mass of first group" },
    { "-cut",      FALSE, etREAL, {&cutoff},
      "Shortest distance (nm) to be considered"},
  };
#define NPA asize(pa)
  FILE       *fp;
  char       *fnTPS,*fnNDX;
  int        status;
  char       **grpname;
  char       outf1[STRLEN],outf2[STRLEN];
  char       title[STRLEN];
  int        natoms,i,j,k,nbin,j0,j1,n;
  int        *count;
  unsigned long int sum;
  real       t,boxmin,hbox,hbox2,cut2,r,r2,invbinw,normfac;
  real       segvol,spherevol,prev_spherevol;
  rvec       *x,xcom,dx;
  real       *inv_segvol;
  bool       *bExcl,bTop;
  matrix     box;
  int        isize[3],*npairs;
  atom_id    *index[3],ix,jx,**pairs;
  t_topology top;
  t_block    *excl;
  t_filenm   fnm[] = {
    { efTRX, "-f", NULL,   ffREAD },
    { efTPS, NULL, NULL,   ffOPTRD },
    { efNDX, NULL, NULL,   ffOPTRD },
    { efXVG, NULL, "rdf",  ffWRITE },
  };
#define NFILE asize(fnm)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_VIEW | PCA_CAN_TIME,TRUE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL);

  fnTPS = ftp2fn_null(efTPS,NFILE,fnm);
  fnNDX = ftp2fn_null(efNDX,NFILE,fnm);
  if (!fnTPS && !fnNDX)
    fatal_error(0,"Neither index file nor topology file specified\n"
		"             Nothing to do!");
  
  excl=NULL;
  if (fnTPS) {
    bTop=read_tps_conf(fnTPS,title,&top,&x,NULL,box,TRUE);
    mk_single_top(&top);
    if (bTop && !bCM)
      /* get exclusions from topology */
      excl=&(top.atoms.excl);
  }
  
  snew(grpname,2);
  if (fnTPS)
    get_index(&top.atoms,fnNDX,2,isize,index,grpname);
  else
    rd_index(fnNDX,2,isize,index,grpname);
  
  natoms=read_first_x(&status,ftp2fn(efTRX,NFILE,fnm),&t,&x,box);
  if ( !natoms )
    fatal_error(0,"Could not read coordinates from statusfile\n");
  if (fnTPS)
    /* check with topology */
    if ( natoms > top.atoms.nr ) 
      fatal_error(0,"Trajectory (%d atoms) does not match topology (%d atoms)",
		  natoms,top.atoms.nr);
  /* check with index groups */
  for (i=0; i<2; i++)
    for (j=0; j<isize[i]; j++)
      if ( index[i][j] >= natoms )
	fatal_error(0,"Atom index (%d) in index group %s (%d atoms) larger "
		    "than number of atoms in trajectory (%d atoms)",
		    index[i][j],grpname[i],isize[i],natoms);
  
  if (bCM) {
    /* move index[0] to index[2] and make 'dummy' index[0] */
    isize[2]=isize[0];
    snew(index[2],isize[2]);
    for(i=0; i<isize[0]; i++)
      index[2][i]=index[0][i];
    isize[0]=1;
    index[0][0]=natoms;
    srenew(index[0],isize[0]);
    /* make space for center of mass */
    srenew(x,natoms+1);
  }
  
  /* initialize some handy things */
  boxmin = min( norm(box[XX]), min( norm(box[YY]), norm(box[ZZ]) ) );
  hbox   = boxmin / 2.0;
  nbin   = (int)(hbox / binwidth) + 1;
  invbinw = 1.0 / binwidth;
  hbox2  = sqr(hbox);
  cut2   = sqr(cutoff);
  
  /* this is THE array */
  snew(count,nbin+1);
  
  /* make pairlist array for groups and exclusions */
  snew(pairs,isize[0]);
  snew(npairs,isize[0]);
  for(i=0; i<isize[0]; i++)
    snew(pairs[i], isize[1]);
  snew(bExcl,natoms);
  for( i = 0; i < isize[0]; i++) {
    ix = index[0][i];
    for( j = 0; j < natoms; j++)
      bExcl[j] = FALSE;
    /* exclusions? */
    if (excl)
      for( j = excl->index[ix]; j < excl->index[ix+1]; j++)
	bExcl[excl->a[j]]=TRUE;
    k = 0;
    for( j = 0; j < isize[1]; j++) {
      jx = index[1][j];
      if (!bExcl[jx])
	pairs[i][k++]=jx;
    }
    npairs[i]=k;
    srenew(pairs[i],npairs[i]);
  }
  sfree(bExcl);
  
  do {
    /* Must init pbc every step because of pressure coupling */
    init_pbc(box,FALSE);
    rm_pbc(&top.idef,natoms,box,x,x);
    
    if (bCM) {
      /* calculate centre of mass */
      clear_rvec(xcom);
      for(i=0; (i < isize[2]); i++) {
	ix = index[2][i];
	rvec_inc(xcom,x[ix]);
      }
      /* store it in the first 'group' */
      for(j=0; (j<DIM); j++)
	x[index[0][0]][j] = xcom[j] / isize[2];
    }
    
    for(i=0; i < isize[0]; i++) {
      ix = index[0][i];
      for(j=0; j < npairs[i]; j++) {
	jx=pairs[i][j];
	pbc_dx(x[ix],x[jx],dx);
	r2=iprod(dx,dx);
	if ( r2 <= hbox2 && r2 > cut2 )
	  count[(int)(sqrt(r2)*invbinw)]++;
      }
    }
  } while (read_next_x(status,&t,natoms,x,box));
  fprintf(stderr,"\n");
  
  close_trj(status);
  
  sfree(x);
  
  /* Calculate volume of sphere segments */
  snew(inv_segvol,nbin);
  prev_spherevol=0;
  for(i=0; (i<nbin); i++) {
    r = (i+1)*binwidth;
    spherevol=(4.0/3.0)*M_PI*r*r*r;
    segvol=spherevol-prev_spherevol;
    inv_segvol[i]=1.0/segvol;
    prev_spherevol=spherevol;
  }
  
  /* Calculate normalization factor and totals */
  sum = 0;
  for(i=0; i<nbin; i++)
    sum += count[i];
  r = nbin*binwidth;
  normfac = (4.0/3.0)*M_PI * r*r*r / sum;
  
  fp=xvgropen(ftp2fn(efXVG,NFILE,fnm),"Radial Distribution","r","");
  fprintf(fp,"@ subtitle \"%s-%s\"\n",grpname[0],grpname[1]);
  sum = 0;
  for(i = 0; i < nbin; i++) {
    sum += count[i];
    fprintf(fp, "%10g %10g %10g\n", (i+0.5)*binwidth, 
	    count[i]*inv_segvol[i]*normfac, (real)sum*normfac);
  }
  ffclose(fp);
  
  xvgr_file(ftp2fn(efXVG,NFILE,fnm),NULL);
  
  thanx(stdout);
  
  return 0;
}
