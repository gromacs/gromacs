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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
static char *SRCID_g_confrms_c = "$Id$";

/* 
 * fit coordinates of file 1 to coordinates of file 2
 * fit is done on the atoms as indicated by the NDX- file 
 */

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
#include "gstat.h"

#define EPS  1.0e-09

void rm_gropbc(t_atoms *atoms,rvec x[],matrix box)
{
  real dist;
  int  n,d;
  
  /* check periodic boundary */
  for(d=0;(d<DIM);d++) {
    for(n=1;(n<atoms->nr);n++) {
      dist = x[n][d]-x[n-1][d];
      if ( fabs(dist) > 0.9 * box[d][d]  ) {
	if ( dist >  0 )
	  x[n][d]-=box[d][d];
	else
	  x[n][d]+=box[d][d];
      } 	
    }
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
    "the two structures will have chain identifiers 'A' and 'B' respectively.",
    "When the option [TT]-one[tt] is set, only the fitted structure is",
    "written to file and the chain identifiers are not changed."
  };
  static bool bSecond=FALSE,bRmpbc=FALSE;
  
  t_pargs pa[] = {
    { "-one", FALSE, etBOOL, &bSecond, "Only write the fitted structure to file" },
    { "-pbc", FALSE, etBOOL, &bRmpbc, "Try to make molecules whole again" }
  };
  t_filenm fnm[] = {
    { efTPS, "-f1",  "conf1.gro", ffREAD  },
    { efSTX, "-f2",  "conf2", ffREAD  },
    { efSTO, "-o",    "fit.pdb",  ffWRITE },
    { efNDX, "-n1" ,  "fit1.ndx",   ffOPTRD  },
    { efNDX, "-n2" ,  "fit2.ndx",   ffOPTRD  }
  };

#define NFILE asize(fnm)
  
  /* the two gromos files */
  FILE    *fp;
  char    title_1[STRLEN],title_2[STRLEN],*name1,*name2;
  t_topology top;
  t_atoms atoms_1,atoms_2;
  int     natoms_1,natoms_2,warn;
  matrix  box;
  atom_id at;
  real    *w_rls,mass,totmass;
  rvec    *x_1,*v_1,*x_2,*v_2,*xf_1,*xf_2,*fit_x;
  matrix  box_1,box_2;
  
  /* counters */
  int   i,j,m;
  
  /* center of mass calculation */
  real tmas_1,tmas_2;
  rvec xcm_1,xcm_2;

  /* variables for fit */
  char *groupnames_1,*groupnames_2;
  int isize_1,isize_2;
  atom_id *index_1,*index_2;
  real rms;

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  
  /* reading first structure file, nb this is going 
   * to be the reference structure*/
  fprintf(stderr,"\nReading first structure file\n");
  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title_1,&top,&x_1,&v_1,box,TRUE);
  atoms_1 = top.atoms;
  fprintf(stderr,"%s\nContaining %d atoms in %d residues\n",
	  title_1,atoms_1.nr,atoms_1.nres);
  srenew(atoms_1.resname,atoms_1.nres);
    
  if ( bRmpbc ) 
    rm_gropbc(&atoms_1,x_1,box_1);


  /* check if we have an NDX file */
  /*
  if ( opt2bSet("-n1",NFILE,fnm) ) {
    fprintf(stderr,"\nSelect group for root least square fit\n");
    rd_index(opt2fn("-n1",NFILE,fnm),1,&isize_1,&index_1,&groupnames_1);
  } else {
    isize_1 = atoms_1.nr;
    snew(index_1,isize_1);
    for(i=0;(i<isize_1);i++)
      index_1[i]=i;
    groupnames_1 = title_1;
  }
  */
  get_index(&atoms_1,opt2fn_null("-n1",NFILE,fnm),
	    1,&isize_1,&index_1,&groupnames_1);
 
  if (isize_1 < 3) 
    fatal_error(0,"Need >= 3 points to fit!\n");


  /* reading second gromos file */
  get_stx_coordnum(opt2fn("-f2",NFILE,fnm),&(atoms_2.nr));
  snew(x_2,atoms_2.nr);
  snew(v_2,atoms_2.nr);
  snew(atoms_2.resname,atoms_2.nr);
  snew(atoms_2.atom,atoms_2.nr);
  snew(atoms_2.atomname,atoms_2.nr);
  fprintf(stderr,"\nReading second structure file\n");
  read_stx_conf(opt2fn("-f2",NFILE,fnm),title_2,&atoms_2,x_2,v_2,box_2);
  fprintf(stderr,"%s\nContaining %d atoms in %d residues\n",
	  title_2,atoms_2.nr,atoms_2.nres);
  srenew(atoms_2.resname,atoms_2.nres);

  if ( bRmpbc ) 
    rm_gropbc(&atoms_2,x_2,box_2);

  get_index(&atoms_2,opt2fn_null("-n2",NFILE,fnm),
	    1,&isize_2,&index_2,&groupnames_2);

  /* check isizes, must be equal */
  if ( isize_2 != isize_1 )
    fatal_error(0,"isize_2 != isize_2\n");

  if (isize_2 < 3) 
    fatal_error(0,"Need >= 3 points to fit!\n");

  for(i=0; (i<isize_1); i++) {
    name1=*atoms_1.atomname[index_1[i]];
    name2=*atoms_2.atomname[index_2[i]];
    if (strcmp(name1,name2)) {
      if (warn<20)
	fprintf(stderr,
		"Warning: atomnames at index %d don't match: %d %s, %d %s\n",
		i+1,index_1[i]+1,name1,index_2[i]+1,name2);
      warn++;
    }
  }
  if (warn)
    fprintf(stderr,"%d atomname%s did not match\n",warn,(warn==1) ? "":"s");

  /* calculate center of mass of reference structure */
  for(m=0;(m<3);m++)
    xcm_1[m]=0;
  for(i=0;(i<isize_1);i++)
    for(m=0;(m<DIM);m++) 
      xcm_1[m]+=x_1[index_1[i]][m];
  for(m=0;(m<3);m++)
    xcm_1[m]/=isize_1;
  for(i=0;(i<atoms_1.nr);i++)
    for(m=0;(m<DIM);m++)
      x_1[i][m]-=xcm_1[m];
  
  /* calculate center of mass of structure to be fitted */
  for(m=0;(m<3);m++)
    xcm_2[m]=0;
  for(i=0;(i<isize_2);i++)
    for(m=0;(m<DIM);m++) 
      xcm_2[m]+=x_2[index_2[i]][m];
  for(m=0;(m<3);m++)
    xcm_2[m]/=isize_2;
  for(i=0;(i<atoms_2.nr);i++)
    for(m=0;(m<DIM);m++)
      x_2[i][m]-=xcm_2[m];
  
  snew(w_rls,atoms_2.nr);
  snew(fit_x,atoms_2.nr);
  for(at=0; at<isize_1; at++) {
    w_rls[index_2[at]] = 1.0;
    copy_rvec(x_1[index_1[at]],fit_x[index_2[at]]);
  }

  /*do the least squares fit to the reference structure*/
  /*
  rms = my_fit(isize_1,index_1,index_2,x_1,atoms_2.nr,x_2);
  */
  do_fit(atoms_2.nr,w_rls,fit_x,x_2);
 
  rms     = 0;
  totmass = 0;
  for(at=0; at<isize_1; at++) {
    mass     = w_rls[index_2[at]];
    for(m=0; m<3; m++)
      rms     += sqr(x_1[index_1[at]][m] - x_2[index_2[at]][m])*mass;
    totmass += mass;
  }
  rms = sqrt(rms/totmass);

  fprintf(stderr,"Root mean square deviation after lsq fit = %g\n",rms);

  /* reset coordinates of reference and fitted structure */
  for(i=0;(i<atoms_1.nr);i++) {
    for(m=0;(m<3);m++)
      x_1[i][m]+=xcm_1[m];
  }
  for(i=0;(i<atoms_2.nr);i++) {
    for(m=0;(m<3);m++)
      x_2[i][m]+=xcm_1[m];
  }

  /* calculate the rms deviation */
  
  

  /* write gromos file of fitted structure(s) */
  fp=ffopen(opt2fn("-o",NFILE,fnm),"w");
  if (fn2ftp(opt2fn("-o",NFILE,fnm))==efGRO) {
    if (!bSecond)
      write_hconf(fp,title_1,&atoms_1,x_1,v_1,box_1);
    write_hconf(fp,title_2,&atoms_2,x_2,v_2,box_2);
  } else {
    if (bSecond)
      write_pdbfile(fp,title_2,&atoms_2,x_2,box_2,0,TRUE);
    else {
      write_pdbfile(fp,title_1,&atoms_1,x_1,box_1,'A',FALSE);
      write_pdbfile(fp,title_2,&atoms_2,x_2,box_2,'B',TRUE);
    }
  }
  fclose(fp);

  thanx(stdout);
  
  return 0;
}









