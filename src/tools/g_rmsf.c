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
 * Great Red Owns Many ACres of Sand 
 */
#include "smalloc.h"
#include "math.h"
#include "typedefs.h"
#include "xvgr.h"
#include "copyrite.h"
#include "statutil.h"
#include "statusio.h"
#include "string2.h"
#include "vec.h"
#include "rdgroup.h"
#include "pdbio.h"
#include "pbc.h"
#include "fatal.h"
#include "futil.h"
#include "gstat.h"

int find_pdb(int npdb,t_pdbatom pdba[],char *resnm,int resnr,char *atomnm)
{
  int i;
  
  resnm[3]='\0';
  for(i=0; (i<npdb); i++) {
    if ((resnr == pdba[i].resnr) &&
	(strcmp(pdba[i].resnm,resnm) == 0) &&
	(strstr(atomnm,pdba[i].atomnm) != NULL))
      break;
  }
  if (i == npdb) {
    fprintf(stderr,"\rCan not find %s%d-%s in pdbfile\n",resnm,resnr,atomnm);
    return -1;
  }
    
  return i;
}

int main (int argc,char *argv[])
{
  static char *desc[] = {
    "g_rmsf computes the root mean square fluctuations (RMSF)",
    "after first fitting to a reference frame.[PAR]",
    "When the (optional) pdb file is given, the RMSF are be converted",
    "to B-factors and plotted with the experimental data."
  };
  static int r0=1;
  t_pargs pa[] = {
    { "-r0", FALSE,etINT,&r0, 
      "starting residue of trajectory (to match the pdb file)" }
  };
#define MAXFRAME 10000
  int          step,nre,natom,natoms,i,g,m,teller=0;
  real         t,lambda,*w_rls,*w_rms,tmas;
  
  t_statheader header;
  t_inputrec   ir;
  t_topology   top;
  t_pdbatom    *pdba;
  int          npdba;
  bool         bCont;

  matrix       box;
  rvec         *x,*v,*xref;
  int          status;
  char         buf[256];
  
  FILE         *fp;               /* the graphics file */
  int          resnr,pdb_i;
  
  atom_id      *index;
  int          isize;
  char         *grpnames;

  real         bfac;
  real         *time;
  atom_id      aid;
  rvec         *rmsf_xx,*rmsf_x;
  real         *rmsf;
  int          d;
  real         count=0;
  rvec         xcm;

  t_filenm fnm[] = {
    { efTPB, NULL, NULL, ffREAD },
    { efTRX, "-f", NULL, ffREAD },
    { efPDB, "-q", NULL, ffOPTRD },
    { efNDX, NULL, NULL, ffOPTRD },
    { efXVG, NULL, NULL, ffWRITE }
  };
#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_CAN_VIEW,TRUE,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);

  read_status_header(ftp2fn(efTPB,NFILE,fnm),&header);
  snew(x,header.natoms);
  snew(v,header.natoms);
  snew(xref,header.natoms);
  snew(w_rls,header.natoms);
  read_status(ftp2fn(efTPB,NFILE,fnm),
              &step,&t,&lambda,&ir,
              box,NULL,NULL,
              &natom,
              xref,NULL,NULL,&nre,NULL,
              &top);

  /*set box type*/
  init_pbc(box,FALSE);
  
  /* allocate memory for time */
  snew(time,MAXFRAME);
  
  fprintf(stderr,"Select group(s) for root mean square calculation\n");
  get_index(&top.atoms,ftp2fn_null(efNDX,NFILE,fnm),1,&isize,&index,&grpnames);

  /* set the weight */
  for(i=0;(i<isize);i++) 
    w_rls[index[i]]=top.atoms.atom[index[i]].m;

  /* malloc the rmsf arrays */
  snew(rmsf_xx,isize);
  snew(rmsf_x, isize);
  snew(rmsf,isize);

  /* remove pbc */
  rm_pbc(&(top.idef),top.atoms.nr,box,xref,xref);

  /* set center of mass to zero */
  tmas=0;
  for(m=0;(m<3);m++)
    xcm[m]=0;
  for(i=0;(i<top.atoms.nr);i++) {
    tmas+=w_rls[i];
    for(m=0;(m<3);m++) {
      xcm[m]+=xref[i][m]*w_rls[i];
    }
  }
  for(m=0;(m<3);m++)
    xcm[m]/=tmas;
  for(i=0;(i<top.atoms.nr);i++) {
    for(m=0;(m<3);m++)
      xref[i][m]-=xcm[m];
  }


  /* read first frame  */
  if ((natoms=read_first_x(&status,
			   ftp2fn(efTRX,NFILE,fnm),&t,&x,box)) != top.atoms.nr) 
    fatal_error(0,"Topology (%d atoms) does not match trajectory (%d atoms)",
		top.atoms.nr,natoms);

  
  do {
    /* remove periodic boundary */
    rm_pbc(&(top.idef),top.atoms.nr,box,x,x);

    /* set center of mass to zero */
    tmas=0;
    for(m=0;(m<3);m++)
      xcm[m]=0;
    for(i=0;(i<top.atoms.nr);i++) {
      tmas+=w_rls[i];
      for(m=0;(m<3);m++) {
	xcm[m]+=x[i][m]*w_rls[i];
      }
    }
    for(m=0;(m<3);m++)
      xcm[m]/=tmas;
    for(i=0;(i<top.atoms.nr);i++) {
      for(m=0;(m<3);m++)
	x[i][m]-=xcm[m];
    }
    
    /* fit to reference structure */
    do_fit(top.atoms.nr,w_rls,xref,x);
 
    /*print time of frame*/
    if ((teller % 10) == 0)
      fprintf(stderr,"\r %5.2f",t);
      
    /*calculate root_least_squares*/
    for(i=0;(i<isize);i++) {
      aid = index[i];
      for(d=0;(d<DIM);d++) {
	rmsf_xx[i][d]+=x[aid][d]*x[aid][d];
	rmsf_x[i][d] +=x[aid][d];
      }
    }
    time[teller]=t;
    bCont=read_next_x(status,&t,natom,x,box);
    count += 1.0;
    teller++;
  } while (bCont);
  close_trj(status);

  for(i=0;(i<isize);i++) {
    rmsf[i] = (rmsf_xx[i][XX]/count - sqr(rmsf_x[i][XX]/count)+ 
	       rmsf_xx[i][YY]/count - sqr(rmsf_x[i][YY]/count)+ 
	       rmsf_xx[i][ZZ]/count - sqr(rmsf_x[i][ZZ]/count)); 
  } 
  
  if (ftp2bSet(efPDB,NFILE,fnm)) {
    fp    = ftp2FILE(efPDB,NFILE,fnm,"r");
    npdba = read_pdbatoms(fp,&pdba,box,FALSE);
    fclose(fp);
    pdba_trimnames(npdba,pdba);
  }
  else
    npdba=0;
  

  /* write output */
  if (npdba == 0) {
    fp=xvgropen(ftp2fn(efXVG,NFILE,fnm),"RMS fluctuation","Residue","nm");
    for(i=0;(i<isize);i++) 
      fprintf(fp,"%5d %8.4f\n",i,sqrt(rmsf[i]));
  }
  else {
    bfac=8.0*M_PI*M_PI/3.0*100;
    fp=xvgropen(ftp2fn(efXVG,NFILE,fnm),"B-Factors",
		"Residue","A\\b\\S\\So\\N\\S 2");
    for(i=0;(i<isize);i++) {
      resnr=top.atoms.atom[index[i]].resnr;
      pdb_i=find_pdb(npdba,pdba,*(top.atoms.resname[resnr]),resnr+r0,
		     *(top.atoms.atomname[index[i]]));
      
      fprintf(fp,"%5d  %10.5f  %10.5f\n",i,rmsf[i]*bfac,
	      (pdb_i == -1) ? 0.0 : pdba[pdb_i].bfac);
    }
  }
  
  fclose(fp);

  xvgr_file(ftp2fn(efXVG,NFILE,fnm),"-nxy");
    
  thanx(stdout);
  
  return 0;
}
