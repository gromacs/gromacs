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
static char *SRCID_g_confrms_c = "$Id$";

#include "filenm.h"
#include "smalloc.h"
#include "macros.h"
#include "math.h"
#include "typedefs.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "tpxio.h"
#include "string2.h"
#include "vec.h"
#include "rdgroup.h"
#include "pbc.h"
#include "fatal.h"
#include "futil.h"
#include "confio.h"
#include "pdbio.h"
#include "txtdump.h"
#include "do_fit.h"

void rm_gropbc(t_atoms *atoms,rvec x[],matrix box)
{
  real dist;
  int  n,d;
  
  /* check periodic boundary */
  for(d=0; d<DIM; d++)
    for(n=1; n<atoms->nr; n++) {
      dist = x[n][d]-x[n-1][d];
      if ( fabs(dist) > 0.9 * box[d][d]  )
	if ( dist >  0 )
	  x[n][d]-=box[d][d];
	else
	  x[n][d]+=box[d][d];
    }
}

int main (int argc,char *argv[])
{
  static char *desc[] = {
    "g_confrms computes the root mean square deviation (RMSD) of two",
    "structures after LSQ fitting the second structure on the first one.",
    "The two structures do NOT need to have the same number of atoms,",
    "only the two index groups used for the fit need to be identical.",
    "[PAR]",
    "The superimposed structures are written to file. In a [TT].pdb[tt] file",
    "the two structures will be written as separate models",
    "(use [TT]rasmol -nmrpdb[tt]).",
  };
  static bool bOne=FALSE,bRmpbc=FALSE;
  
  t_pargs pa[] = {
    { "-one",FALSE,etBOOL,{&bOne},"Only write the fitted structure to file" },
    { "-pbc",FALSE,etBOOL,{&bRmpbc},"Try to make molecules whole again" }
  };
  t_filenm fnm[] = {
    { efTPS, "-f1",  "conf1.gro", ffREAD  },
    { efSTX, "-f2",  "conf2",     ffREAD  },
    { efSTO, "-o",   "fit.pdb",   ffWRITE },
    { efNDX, "-n1" , "fit1.ndx",  ffOPTRD },
    { efNDX, "-n2" , "fit2.ndx",  ffOPTRD }
  };
#define NFILE asize(fnm)
  
  /* the two structure files */
  FILE    *fp;
  char    title1[STRLEN],title2[STRLEN],*name1,*name2;
  t_topology top;
  t_atoms atoms1,atoms2;
  int     natoms1,natoms2,warn=0;
  matrix  box;
  atom_id at;
  real    *w_rls,mass,totmass;
  rvec    *x1,*v1,*x2,*v2,*fit_x;
  matrix  box1,box2;
  
  /* counters */
  int     i,j,m;
  
  /* center of mass calculation */
  real    tmas1,tmas2;
  rvec    xcm1,xcm2;
  
  /* variables for fit */
  char    *groupnames1,*groupnames2;
  int     isize1,isize2;
  atom_id *index1,*index2;
  real    rms;
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  /* reading reference structure from first structure file */
  fprintf(stderr,"\nReading first structure file\n");
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title1,&top,&x1,&v1,box,TRUE);
  atoms1 = top.atoms;
  fprintf(stderr,"%s\nContaining %d atoms in %d residues\n",
	  title1,atoms1.nr,atoms1.nres);
  srenew(atoms1.resname,atoms1.nres);
  
  if ( bRmpbc ) 
    rm_gropbc(&atoms1,x1,box1);
  
  get_index(&atoms1,opt2fn_null("-n1",NFILE,fnm),
	    1,&isize1,&index1,&groupnames1);
  
  if (isize1 < 3) 
    fatal_error(0,"Need >= 3 points to fit!\n");
  
  /* reading second structure file */
  get_stx_coordnum(opt2fn("-f2",NFILE,fnm),&(atoms2.nr));
  snew(x2,atoms2.nr);
  snew(v2,atoms2.nr);
  snew(atoms2.resname,atoms2.nr);
  snew(atoms2.atom,atoms2.nr);
  snew(atoms2.atomname,atoms2.nr);
  fprintf(stderr,"\nReading second structure file\n");
  read_stx_conf(opt2fn("-f2",NFILE,fnm),title2,&atoms2,x2,v2,box2);
  fprintf(stderr,"%s\nContaining %d atoms in %d residues\n",
	  title2,atoms2.nr,atoms2.nres);
  srenew(atoms2.resname,atoms2.nres);
  
  if ( bRmpbc ) 
    rm_gropbc(&atoms2,x2,box2);
  
  get_index(&atoms2,opt2fn_null("-n2",NFILE,fnm),
	    1,&isize2,&index2,&groupnames2);
  
  /* check isizes, must be equal */
  if ( isize2 != isize1 )
    fatal_error(0,"isize2 != isize2\n");
  
  if (isize2 < 3) 
    fatal_error(0,"Need >= 3 points to fit!\n");
  
  for(i=0; i<isize1; i++) {
    name1=*atoms1.atomname[index1[i]];
    name2=*atoms2.atomname[index2[i]];
    if (strcmp(name1,name2)) {
      if (warn<20)
	fprintf(stderr,
		"Warning: atomnames at index %d don't match: %u %s, %u %s\n",
		i+1,index1[i]+1,name1,index2[i]+1,name2);
      warn++;
    }
  }
  if (warn)
    fprintf(stderr,"%d atomname%s did not match\n",warn,(warn==1) ? "":"s");
  
  /* calculate center of mass of reference structure */
  for(m=0; m<DIM; m++) {
    xcm1[m]=0;
    for(i=0; i<isize1; i++)
      xcm1[m]+=x1[index1[i]][m];
    xcm1[m]/=isize1;
    for(i=0; i<atoms1.nr; i++)
      x1[i][m]-=xcm1[m];
  }
  
  /* calculate center of mass of structure to be fitted */
  for(m=0; m<DIM; m++) {
    xcm2[m]=0;
    for(i=0; i<isize2; i++)
      xcm2[m]+=x2[index2[i]][m];
    xcm2[m]/=isize2;
    for(i=0; i<atoms2.nr; i++)
      x2[i][m]-=xcm2[m];
  }
  
  snew(w_rls,atoms2.nr);
  snew(fit_x,atoms2.nr);
  for(at=0; at<isize1; at++) {
    w_rls[index2[at]] = 1.0;
    copy_rvec(x1[index1[at]],fit_x[index2[at]]);
  }
  
  /* do the least squares fit to the reference structure */
  do_fit(atoms2.nr,w_rls,fit_x,x2);
  
  /* calculate the rms deviation */
  rms     = 0;
  totmass = 0;
  for(at=0; at<isize1; at++) {
    mass = w_rls[index2[at]];
    for(m=0; m<DIM; m++)
      rms += sqr(x1[index1[at]][m] - x2[index2[at]][m])*mass;
    totmass += mass;
  }
  rms = sqrt(rms/totmass);
  
  fprintf(stderr,"Root mean square deviation after lsq fit = %g\n",rms);
  
  /* reset coordinates of reference and fitted structure */
  for(i=0; i<atoms1.nr; i++)
    for(m=0; m<DIM; m++)
      x1[i][m]+=xcm1[m];
  for(i=0; i<atoms2.nr; i++)
    for(m=0; m<DIM; m++)
      x2[i][m]+=xcm1[m];
  
  /* write gromos file of fitted structure(s) */
  fp=ffopen(opt2fn("-o",NFILE,fnm),"w");
  if (fn2ftp(opt2fn("-o",NFILE,fnm))==efGRO) {
    if (!bOne)
      write_hconf_p(fp,title1,&atoms1,3,x1,v1,box1);
    write_hconf_p(fp,title2,&atoms2,3,x2,v2,box2);
  } else {
    if (!bOne)
      write_pdbfile(fp,title1,&atoms1,x1,box1,0,1);
    write_pdbfile(fp,title2,&atoms2,x2,box2,0,bOne ? -1 : 2);
  }
  fclose(fp);
  
  thanx(stderr);
  
  return 0;
}
