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
 * Giant Rising Ordinary Mutants for A Clerical Setup
 */
static char *SRCID_gmxcheck_c = "$Id$";

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include "main.h"
#include "macros.h"
#include "math.h"
#include "futil.h"
#include "statutil.h"
#include "copyrite.h"
#include "sysstuff.h"
#include "txtdump.h"
#include "fatal.h"
#include "gmxfio.h"
#include "trnio.h"
#include "xtcio.h"
#include "tpbcmp.h"
#include "vdw.h"
#include "mass.h"
#include "vec.h"
#include "pbc.h"
#include "physics.h"
#include "smalloc.h"
#include "confio.h"
#include "enxio.h"

void chk_trj(char *fn)
{
  t_trnheader  sh,count;
  int          idum,j,natoms,step;
  real         rdum,t,t0,old_t1,old_t2,prec;
  bool         bShowTimestep=TRUE;
  rvec         *x;
  matrix       box;
  size_t       fpos;
  int          xd;
  int          status,ftp;
  
#define BOGUSTIME -1e10

  ftp    = fn2ftp(fn);
  natoms = 0;  
  t      = 0;
  t0     = BOGUSTIME;
  step   = 0;
  count.box_size=0;
  count.vir_size=0;
  count.pres_size=0;
  count.x_size=0;
  count.v_size=0;
  count.f_size=0;

  printf("Checking file %s\n",fn);
  switch (ftp) {
  case efTRJ:
  case efTRR:
    status = open_trn(fn,"r");
    j      =  0;
    t      = -1;
    old_t2 = -2.0;
    old_t1 = -1.0;
    fpos   = fio_ftell(status);
    while (fread_trnheader(status,&sh)) {
      if (j>=2) {
	if ( fabs((sh.t-old_t1)-(old_t1-old_t2)) > 
	     0.1*(fabs(sh.t-old_t1)+fabs(old_t1-old_t2)) ) {
	  bShowTimestep=FALSE;
	  fprintf(stderr,"\nTimesteps at t=%g don't match (%g, %g)\n",
		  old_t1,old_t1-old_t2,sh.t-old_t1);
	}
      }
      old_t2=old_t1;
      old_t1=sh.t;
      if (t0 == BOGUSTIME) t0=sh.t;
      fprintf(stderr,"\rframe: %6d, t: %10.3f bytes: %10u",j,sh.t,fpos);
      if (j == 0)
	fprintf(stderr,"\n");
      fread_htrn(status,&sh,NULL,NULL,NULL,NULL);
      j++;
#define INC(s,n,item) if (s.item  != 0) n.item++
      INC(sh,count,box_size);
      INC(sh,count,vir_size);
      INC(sh,count,pres_size);
      INC(sh,count,x_size);
      INC(sh,count,v_size);
      INC(sh,count,f_size);
    }
    close_trn(status);
    t=sh.t;
    break;
  case efXTC:
    xd = open_xtc(fn,"r");
    if (read_first_xtc(xd,&natoms,&step,&t,box,&x,&prec)) {
      fprintf(stderr,"\nXTC precision %g\n\n",prec);
      j=0;
      old_t2=-2.0;
      old_t1=-1.0;
      do {
	if (j>=2) {
	  if ( fabs((t-old_t1)-(old_t1-old_t2)) > 
	       0.1*(fabs(t-old_t1)+fabs(old_t1-old_t2)) ) {
	    bShowTimestep=FALSE;
	    fprintf(stderr,"\nTimesteps at t=%g don't match (%g, %g)\n",
		    old_t1,old_t1-old_t2,t-old_t1);
	  }
	}
	old_t2=old_t1;
	old_t1=t;
	if (t0 == BOGUSTIME) t0=t;
	fprintf(stderr,"\rframe: %6d, t: %10.3f",j,t);
	if (j == 0)
	  fprintf(stderr,"\n");
	count.x_size++;
	count.box_size++;
	j++;
      } while (read_next_xtc(xd,&natoms,&step,&t,box,x,&prec));
      close_xtc(xd);
    }
    else
      fprintf(stderr,"Empty file %s\n",fn);
    break;
  default:
    fprintf(stderr,"Sorry %s not supported yet\n",fn);
  }

  fprintf(stderr,"\n\n");
  fprintf(stderr,"\n# Atoms     %d\n",natoms);
  
  fprintf(stderr,"\nItem        #frames");
  if (bShowTimestep)
    fprintf(stderr," Timestep (ps)");
  fprintf(stderr,"\n");
#define PRINTITEM(label,item) fprintf(stderr,"%-10s  %6d",label,count.item); if ((bShowTimestep) && (count.item > 1)) fprintf(stderr,"    %g\n",(t-t0)/(count.item-1)); else fprintf(stderr,"\n")
  PRINTITEM ( "Energies",   e_size );
  PRINTITEM ( "Box",        box_size );
  PRINTITEM ( "Virial",     vir_size );
  PRINTITEM ( "Pressure",   pres_size );
  PRINTITEM ( "Coords",     x_size );
  PRINTITEM ( "Velocities", v_size );
  PRINTITEM ( "Forces",     f_size );
}  

void chk_ndx(char *fn)
{
  FILE *in;
  char buf[256];
  int ngrp,nelem,nitem,nat,ng_cnt,na_cnt;
  int nLab;

  fprintf(stderr,"Checking index file %s\n",fn);
  in=ffopen(fn,"r");
  
  if (fscanf(in,"%d%d",&ngrp,&nat) != 2) {
    fprintf(stderr,"Couldn't read NGRP or NAT in file %s\n",fn);
    exit(1);
  }
  nLab=2;
  ng_cnt=0,na_cnt=0;
  fprintf(stderr,"There should be %d groups containing %d atoms in total\n",
	  ngrp,nat);
  while (fscanf(in,"%s",buf) == 1) {
    if (nLab == 2) {
      fprintf(stderr,"Group %4d: %16s",ng_cnt,buf);
      ng_cnt++;
      nLab--;
    }
    else if (nLab == 1) {
      fprintf(stderr,"  Nelem: %16s\n",buf);
      if (sscanf(buf,"%d",&nelem) != 1) {
	fprintf(stderr,"\nERROR in indexfile %s\n",fn);
	fprintf(stderr,"Couldn't find proper NATOMS in buf '%s'\n",buf);
	fprintf(stderr,"May be you have entered too many atoms in the previous group\n\n");
	ffclose(in);
	return;
      }
      nLab--;
    }
    else {
      if (sscanf(buf,"%d",&nitem) != 1) {
	fprintf(stderr,"\nERROR in indexfile %s\n",fn);
	fprintf(stderr,"Couldn't find proper ATOM in buf '%s'\n",buf);
	fprintf(stderr,"You should have entered %d more atom(s) for the previous group\n\n",nelem);
	ffclose(in);
	return;
      }
      nelem--;
      na_cnt++;
      if (nelem == 0)
	nLab=2;
    }
  }
  fprintf(stderr,"\nFound %6d groups, ",ng_cnt);
  if (ng_cnt != ngrp)
    fprintf(stderr,"should be %6d\n",ngrp);
  else
    fprintf(stderr,"as it should be.\n");
  fprintf(stderr,"Found %6d atoms,  ",na_cnt);
  if (na_cnt != nat)
    fprintf(stderr,"should be %6d\n",nat);
  else
    fprintf(stderr,"as it should be.\n");
  
  ffclose(in);
}

void chk_stx(char *fn)
{
  int       natom,i,j,k,nvdw;
  char      title[STRLEN];
  t_atoms   atoms;
  rvec      *x,*v;
  rvec      dx;
  matrix    box;
  bool      bV,bX,bB,bFirst,bOut;
  real      r2,ekin,temp1,temp2;
  t_vdw     *vdw;
  real      *atom_vdw;
  
  fprintf(stderr,"Checking coordinate file %s\n",fn);
  atoms.nr=0;
  atoms.nres=0;
  get_stx_coordnum(fn,&natom);
  fprintf(stderr,"%d atoms in file\n",natom);
  snew(atoms.atomname,natom);
  snew(atoms.resname,natom);
  snew(atoms.atom,natom);
  snew(x,natom);
  snew(v,natom);
  read_stx_conf(fn,title,&atoms,x,v,box);
  
  /* check coordinates and box */
  bV=FALSE;
  bX=FALSE;
  for (i=0; (i<natom) && !(bV && bX); i++)
    for (j=0; (j<DIM) && !(bV && bX); j++) {
      bV=bV || (v[i][j]!=0);
      bX=bX || (x[i][j]!=0);
    }
  bB=FALSE;
  for (i=0; (i<DIM) && !bB; i++)
    for (j=0; (j<DIM) && !bB; j++)
      bB=bB || (box[i][j]!=0);
  
  fprintf(stderr,"%scoordinates found\n",bX?"":"No ");
  fprintf(stderr,"%sbox         found\n",bB?"":"No ");
  fprintf(stderr,"%svelocities  found\n",bV?"":"No ");
  fprintf(stderr,"\n");
  
  /* check velocities */
  if (bV) {
    for (i=0; (i<natom); i++)
      if ((atoms.atom[i].m = get_mass(*atoms.resname[atoms.atom[i].resnr],
				      *atoms.atomname[i])) == 0.0)
	if ( ((*(atoms.atomname[i]))[0]=='H') ||
	     (isdigit((*(atoms.atomname[i]))[0]) && 
	      ((*(atoms.atomname[i]))[1]=='H')) )
	  atoms.atom[i].m=1.008; /* proton mass */
	else
	  atoms.atom[i].m=12.0110; /* carbon mass */
    ekin=0.0;
    for (i=0; (i<natom); i++)
      for (j=0; (j<DIM); j++)
	ekin+=0.5*atoms.atom[i].m*v[i][j]*v[i][j];
    temp1=(2.0*ekin)/(natom*DIM*BOLTZ); 
    temp2=(2.0*ekin)/(natom*(DIM-1)*BOLTZ); 
    fprintf(stderr,"Assuming the number of degrees of freedom to be "
	    "Natoms * 3 or Natoms * 2,\n"
	    "the velocities correspond to a temperature of the system\n"
	    "of %g K or %g K respectively.\n\n",temp1,temp2);
  }
  
  /* check coordinates */
  if (bX) {
    fprintf(stderr,"Checking for atoms closer than 80%% of smallest "
		    "VanderWaals distance:\n");
    nvdw=read_vdw("radii.vdw",&vdw);
    snew(atom_vdw,natom);
    for (i=0; (i<natom); i++)
      if ((atom_vdw[i]=get_vdw(nvdw,vdw,*(atoms.atomname[i])))==0.0)
	if ( ((*atoms.atomname[i])[0]=='H') ||
	     (isdigit((*(atoms.atomname[i]))[0]) && 
	      ((*atoms.atomname[i])[1]=='H')) )
	  atom_vdw[i]=0.1;
	else
	  atom_vdw[i]=0.2;
      
    if (bB) 
      init_pbc(box,FALSE);
      
    bFirst=TRUE;
    for (i=0; (i<natom); i++) {
      if (((i+1)%10)==0)
	fprintf(stderr,"\r%5d",i+1);
      for (j=i+1; (j<natom); j++) {
	if (bB)
	  pbc_dx(box,x[i],x[j],dx);
	else
	  rvec_sub(x[i],x[j],dx);
	r2=iprod(dx,dx);
	if (r2<=sqr(0.8*min(atom_vdw[i],atom_vdw[j]))) {
	  if (bFirst) {
	    fprintf(stderr,"\r%5s %4s %8s %5s  %5s %4s %8s %5s  %6s\n",
		    "atom#","name","residue","r_vdw",
		    "atom#","name","residue","r_vdw","distance");
	    bFirst=FALSE;
	  }
	  fprintf(stderr,
		  "\r%5d %4s %4s%4d %-5.3g  %5d %4s %4s%4d %-5.3g  %-6.4g\n",
		  i+1,*(atoms.atomname[i]),
		  *(atoms.resname[atoms.atom[i].resnr]),atoms.atom[i].resnr+1,
		  atom_vdw[i],
		  j+1,*(atoms.atomname[j]),
		  *(atoms.resname[atoms.atom[j].resnr]),atoms.atom[j].resnr+1,
		  atom_vdw[j],
		  sqrt(r2) );
	}
      }
    }
    if (bFirst) 
      fprintf(stderr,"\rno close atoms found\n");
    fprintf(stderr,"\r      \n");
    
    if (bB) {
      /* check box */
      bFirst=TRUE;
      k=0;
      for (i=0; (i<natom) && (k<10); i++) {
	bOut=FALSE;
	for(j=0; (j<DIM) && !bOut; j++)
	  bOut = bOut || (x[i][j]<0) || (x[i][j]>box[j][j]);
	if (bOut) {
	  k++;
	  if (bFirst) {
	    fprintf(stderr,"Atoms outide box ( ");
	    for (j=0; (j<DIM); j++)
	      fprintf(stderr,"%g ",box[j][j]);
	    fprintf(stderr,"):\n%5s %4s %8s %5s  %s\n",
		    "atom#","name","residue","r_vdw","coordinate");
	    bFirst=FALSE;
	  }
	  fprintf(stderr,
		  "%5d %4s %4s%4d %-5.3g",
		  i,*(atoms.atomname[i]),*(atoms.resname[atoms.atom[i].resnr]),
		  atoms.atom[i].resnr,atom_vdw[i]);
	  for (j=0; (j<DIM); j++)
	    fprintf(stderr," %6.3g",x[i][j]);
	  fprintf(stderr,"\n");
	}
      }
      if (k==10)
	fprintf(stderr,"(maybe more)\n");
      if (bFirst) 
	fprintf(stderr,"no atoms found outside box\n");
      fprintf(stderr,"\n");
    }
  }
}

void chk_enx(char *fn)
{
  int      in,nre,frame,fnr;
  char     **enm=NULL;
  t_energy *ee=NULL;
  bool     bShowTStep;
  real     t,t0,old_t1,old_t2;
  
  fprintf(stderr,"Checking energy file %s\n\n",fn);

  in = open_enx(fn,"r");
  do_enxnms(in,&nre,&enm);
  fprintf(stderr,"%d groups in energy file",nre);
  old_t2=-2.0;
  old_t1=-1.0;
  fnr=0;
  t0=BOGUSTIME;
  bShowTStep=TRUE;
  snew(ee,nre);
  while (do_enx(in,&t,&frame,&nre,ee,NULL)) {
    if (fnr>=2) {
      if ( fabs((t-old_t1)-(old_t1-old_t2)) > 
	   0.1*(fabs(t-old_t1)+fabs(old_t1-old_t2)) ) {
	bShowTStep=FALSE;
	fprintf(stderr,"\nTimesteps at t=%g don't match (%g, %g)\n",
		old_t1,old_t1-old_t2,t-old_t1);
      }
    }
    old_t2=old_t1;
    old_t1=t;
    if (t0 == BOGUSTIME) t0=t;
    fprintf(stderr,"\rframe: %6d (index %6d), t: %10.3f",frame,fnr,t);
    if (fnr == 0)
      fprintf(stderr,"\n");
    fnr++;
  }
  fprintf(stderr,"\n\nFound %d frames",fnr);
  if (bShowTStep)
    fprintf(stderr," with a timestep of %g ps",(t-t0)/(fnr-1));
  fprintf(stderr,".\n",fnr);

}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "gmxcheck reads a trajectory ([TT].trj[tt], [TT].trr[tt] or ",
    "[TT].xtc[tt]), an index file ([TT].ndx[tt]) or an energy file",
    "([TT].ene[tt] or [TT].edr])",
    "and prints out useful information about them.[PAR]",
    "For a coordinate file (generic structure file, e.g. [TT].gro[tt]) ",
    "gmxcheck will check for presence of coordinates, velocities and box",
    "in the file, for close contacts (smaller than 80% of the sum",
    "of both Vanderwaals radii) and atoms outside the box (these may occur",
    "often and are no problem). If velocities are present, an estimated",
    "temperature will be calculated from them.[PAR]",
    "The program will compare run input ([TT].tpr[tt], [TT].tpb[tt] or",
    "[TT].tpa[tt]) files",
    "when both [TT]-s1[tt] and [TT]-s2[tt] are supplied."
  };
  t_filenm fnm[] = {
    { efTRX, "-f", NULL, ffOPTRD },
    { efNDX, "-n", NULL, ffOPTRD },
    { efTPX, "-s1", "top1", ffOPTRD },
    { efTPX, "-s2", "top2", ffOPTRD },
    { efSTX, "-c", NULL, ffOPTRD },
    { efENX, "-e", NULL, ffOPTRD }
  };
#define NFILE asize(fnm)
  char *fn1=NULL,*fn2=NULL;

  CopyRight(stdout,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,0,NULL,
		    asize(desc),desc,0,NULL);
  
  if (ftp2bSet(efTRX,NFILE,fnm))
    chk_trj(ftp2fn(efTRX,NFILE,fnm));
  
  if (ftp2bSet(efNDX,NFILE,fnm))
    chk_ndx(ftp2fn(efNDX,NFILE,fnm));
  
  if (opt2bSet("-s1",NFILE,fnm))
    fn1=opt2fn("-s1",NFILE,fnm);
  if (opt2bSet("-s2",NFILE,fnm))
    fn2=opt2fn("-s2",NFILE,fnm);
  if (fn1 && fn2)
    comp_tpx(fn1,fn2);
  else if (fn1 || fn2)
    fprintf(stderr,"Please give me TWO run input (.tpr/.tpa/.tpb) files!\n");
  
  if (ftp2bSet(efSTX,NFILE,fnm))
    chk_stx(ftp2fn(efSTX,NFILE,fnm));
  
  if (ftp2bSet(efENX, NFILE,fnm))
    chk_enx(ftp2fn(efENX,NFILE,fnm));
  
  thanx(stderr);
  
  return 0;
}
