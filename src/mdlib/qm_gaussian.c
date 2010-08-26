/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_QMMM_GAUSSIAN

#include <math.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "assert.h"
#include "physics.h"
#include "macros.h"
#include "vec.h"
#include "force.h"
#include "invblock.h"
#include "confio.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "copyrite.h"
#include "qmmm.h"
#include <stdio.h>
#include <string.h>
#include "gmx_fatal.h"
#include "typedefs.h"
#include <stdlib.h>


/* TODO: this should be made thread-safe */

/* Gaussian interface routines */

void init_gaussian(t_commrec *cr, t_QMrec *qm, t_MMrec *mm)
{
  FILE    
    *rffile=NULL,*out=NULL;
  ivec
    basissets[eQMbasisNR]={{0,3,0},
			   {0,3,0},/*added for double sto-3g entry in names.c*/
			   {5,0,0},
			   {5,0,1},
			   {5,0,11},
			   {5,6,0},
			   {1,6,0},
			   {1,6,1},
			   {1,6,11},
			   {4,6,0}};
  char
    *buf;
  int
    i;
  
  /* using the ivec above to convert the basis read form the mdp file
   * in a human readable format into some numbers for the gaussian
   * route. This is necessary as we are using non standard routes to
   * do SH.
   */

  /* per layer we make a new subdir for integral file, checkpoint
   * files and such. These dirs are stored in the QMrec for
   * convenience 
   */

  
  if(!qm->nQMcpus){ /* this we do only once per layer 
		     * as we call g01 externally 
		     */

    for(i=0;i<DIM;i++)
      qm->SHbasis[i]=basissets[qm->QMbasis][i];

  /* init gradually switching on of the SA */
    qm->SAstep = 0;
  /* we read the number of cpus and environment from the environment
   * if set.  
   */
    snew(buf,20);
    buf = getenv("NCPUS");
    if (buf)
      sscanf(buf,"%d",&qm->nQMcpus);
    else
      qm->nQMcpus=1;
    fprintf(stderr,"number of CPUs for gaussian = %d\n",qm->nQMcpus);
    snew(buf,50);
    buf = getenv("MEM");
    if (buf)
      sscanf(buf,"%d",&qm->QMmem);
    else
      qm->QMmem=50000000;
    fprintf(stderr,"memory for gaussian = %d\n",qm->QMmem);
    snew(buf,30);
    buf = getenv("ACC");
    if (buf)
      sscanf(buf,"%d",&qm->accuracy);
    else
      qm->accuracy=8;  
    fprintf(stderr,"accuracy in l510 = %d\n",qm->accuracy); 
    snew(buf,30);
    buf = getenv("CPMCSCF");
    if (buf)
	{
		sscanf(buf,"%d",&i);
		qm->cpmcscf = (i!=0);
	}
	else
      qm->cpmcscf=FALSE;
    if (qm->cpmcscf)
      fprintf(stderr,"using cp-mcscf in l1003\n");
    else
      fprintf(stderr,"NOT using cp-mcscf in l1003\n"); 
    snew(buf,50);
    buf = getenv("SASTEP");
    if (buf)
      sscanf(buf,"%d",&qm->SAstep);
    else
      /* init gradually switching on of the SA */
      qm->SAstep = 0;
    /* we read the number of cpus and environment from the environment
     * if set.  
     */
    fprintf(stderr,"Level of SA at start = %d\n",qm->SAstep);
        

    /* punch the LJ C6 and C12 coefficients to be picked up by
     * gaussian and usd to compute the LJ interaction between the
     * MM and QM atoms.
     */
    if(qm->bTS||qm->bOPT){
      out = fopen("LJ.dat","w");
      for(i=0;i<qm->nrQMatoms;i++){

#ifdef GMX_DOUBLE
	fprintf(out,"%3d  %10.7lf  %10.7lf\n",
		qm->atomicnumberQM[i],qm->c6[i],qm->c12[i]);
#else
	fprintf(out,"%3d  %10.7f  %10.7f\n",
		qm->atomicnumberQM[i],qm->c6[i],qm->c12[i]);
#endif
      }
      fclose(out);
    }
    /* gaussian settings on the system */
    snew(buf,200);
    buf = getenv("GAUSS_DIR");
    fprintf(stderr,"%s",buf);

    if (buf){
      snew(qm->gauss_dir,200);
      sscanf(buf,"%s",qm->gauss_dir);
    }
    else
      gmx_fatal(FARGS,"no $GAUSS_DIR, check gaussian manual\n");
    
    snew(buf,200);    
    buf = getenv("GAUSS_EXE");
    if (buf){
      snew(qm->gauss_exe,200);
      sscanf(buf,"%s",qm->gauss_exe);
    }
    else
      gmx_fatal(FARGS,"no $GAUSS_EXE, check gaussian manual\n");
    
    snew(buf,200);
    buf = getenv("DEVEL_DIR");
    if (buf){
      snew(qm->devel_dir,200);
      sscanf(buf,"%s",qm->devel_dir);
    }
    else
      gmx_fatal(FARGS,"no $DEVEL_DIR, this is were the modified links reside.\n");
    
    
    /*  if(fr->bRF){*/
    /* reactionfield, file is needed using gaussian */
    /*    rffile=fopen("rf.dat","w");*/
    /*   fprintf(rffile,"%f %f\n",fr->epsilon_r,fr->rcoulomb/BOHR2NM);*/
    /* fclose(rffile);*/
    /*  }*/
  }
  fprintf(stderr,"gaussian initialised...\n");
}  



void write_gaussian_SH_input(int step,gmx_bool swap,
			     t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
  int
    i;
  gmx_bool
    bSA;
  FILE
    *out;
  t_QMMMrec
    *QMMMrec;
  QMMMrec = fr->qr;
  bSA = (qm->SAstep>0);

  out = fopen("input.com","w");
  /* write the route */
  fprintf(out,"%s","%scr=input\n");
  fprintf(out,"%s","%rwf=input\n");
  fprintf(out,"%s","%int=input\n");
  fprintf(out,"%s","%d2e=input\n");
/*  if(step)
 *   fprintf(out,"%s","%nosave\n");
 */
  fprintf(out,"%s","%chk=input\n");
  fprintf(out,"%s%d\n","%mem=",qm->QMmem);
  fprintf(out,"%s%3d\n","%nprocshare=",qm->nQMcpus);

  /* use the versions of
   * l301 that computes the interaction between MM and QM atoms.
   * l510 that can punch the CI coefficients
   * l701 that can do gradients on MM atoms 
   */

  /* local version */
  fprintf(out,"%s%s%s",
	  "%subst l510 ",
	  qm->devel_dir,
	  "/l510\n");
  fprintf(out,"%s%s%s",
	  "%subst l301 ",
	  qm->devel_dir,
	  "/l301\n");
  fprintf(out,"%s%s%s",
	  "%subst l701 ",
	  qm->devel_dir,
	  "/l701\n");
  
  fprintf(out,"%s%s%s",
	  "%subst l1003 ",
	  qm->devel_dir,
	  "/l1003\n");
  fprintf(out,"%s%s%s",
	  "%subst l9999 ",
	  qm->devel_dir,
	  "/l9999\n");
  /* print the nonstandard route 
   */
  fprintf(out,"%s",
	  "#P nonstd\n 1/18=10,20=1,38=1/1;\n");
  fprintf(out,"%s",
	  " 2/9=110,15=1,17=6,18=5,40=1/2;\n");
  if(mm->nrMMatoms)
    fprintf(out,
	    " 3/5=%d,6=%d,7=%d,25=1,32=1,43=1,94=-2/1,2,3;\n",
	    qm->SHbasis[0],
	    qm->SHbasis[1],
	    qm->SHbasis[2]); /*basisset stuff */
  else
    fprintf(out,
	    " 3/5=%d,6=%d,7=%d,25=1,32=1,43=0,94=-2/1,2,3;\n",
	    qm->SHbasis[0],
	    qm->SHbasis[1],
	    qm->SHbasis[2]); /*basisset stuff */
  /* development */
  if (step+1) /* fetch initial guess from check point file */
    /* hack, to alyays read from chk file!!!!! */
    fprintf(out,"%s%d,%s%d%s"," 4/5=1,7=6,17=",
	    qm->CASelectrons,
	    "18=",qm->CASorbitals,"/1,5;\n");
  else /* generate the first checkpoint file */
    fprintf(out,"%s%d,%s%d%s"," 4/5=0,7=6,17=",
	    qm->CASelectrons,
	    "18=",qm->CASorbitals,"/1,5;\n");
  /* the rest of the input depends on where the system is on the PES 
   */
  if(swap && bSA){ /* make a slide to the other surface */
    if(qm->CASorbitals>6){  /* use direct and no full diag */
      fprintf(out," 5/5=2,16=-2,17=10000000,28=2,32=2,38=6,97=100/10;\n");
    } 
    else {
      if(qm->cpmcscf){
	fprintf(out," 5/5=2,6=%d,17=31000200,28=2,32=2,38=6,97=100/10;\n",
		qm->accuracy);
	if(mm->nrMMatoms>0)
	  fprintf(out," 7/7=1,16=-2,30=1/1;\n");
	fprintf(out," 11/31=1,42=1,45=1/1;\n");
	fprintf(out," 10/6=1,10=700006,28=2,29=1,31=1,97=100/3;\n");
	fprintf(out," 7/30=1/16;\n 99/10=4/99;\n");
      }
      else{
	fprintf(out," 5/5=2,6=%d,17=11000000,28=2,32=2,38=6,97=100/10;\n",
		qm->accuracy);
	fprintf(out," 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
      }
    }
  }
  else if(bSA){ /* do a "state-averaged" CAS calculation */
    if(qm->CASorbitals>6){ /* no full diag */ 
      fprintf(out," 5/5=2,16=-2,17=10000000,28=2,32=2,38=6/10;\n");
    } 
    else {
      if(qm->cpmcscf){
	fprintf(out," 5/5=2,6=%d,17=31000200,28=2,32=2,38=6/10;\n",
		qm->accuracy);
	if(mm->nrMMatoms>0)
	  fprintf(out," 7/7=1,16=-2,30=1/1;\n");
	fprintf(out," 11/31=1,42=1,45=1/1;\n");
	fprintf(out," 10/6=1,10=700006,28=2,29=1,31=1/3;\n");
	fprintf(out," 7/30=1/16;\n 99/10=4/99;\n");
      }
      else{
      	fprintf(out," 5/5=2,6=%d,17=11000000,28=2,32=2,38=6/10;\n",
		qm->accuracy);
	fprintf(out," 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
      }
    }
  }
  else if(swap){/* do a "swapped" CAS calculation */
    if(qm->CASorbitals>6)
      fprintf(out," 5/5=2,16=-2,17=0,28=2,32=2,38=6,97=100/10;\n");
    else
      fprintf(out," 5/5=2,6=%d,17=1000000,28=2,32=2,38=6,97=100/10;\n",
	      qm->accuracy);
    fprintf(out," 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
  }
  else {/* do a "normal" CAS calculation */
    if(qm->CASorbitals>6)
      fprintf(out," 5/5=2,16=-2,17=0,28=2,32=2,38=6/10;\n");
    else
      fprintf(out," 5/5=2,6=%d,17=1000000,28=2,32=2,38=6/10;\n",
	      qm->accuracy);
    fprintf(out," 7/7=1,16=-2,30=1/1,2,3,16;\n 99/10=4/99;\n");
  }
  fprintf(out, "\ninput-file generated by gromacs\n\n");
  fprintf(out,"%2d%2d\n",qm->QMcharge,qm->multiplicity);
  for (i=0;i<qm->nrQMatoms;i++){
#ifdef GMX_DOUBLE
    fprintf(out,"%3d %10.7lf  %10.7lf  %10.7lf\n",
	    qm->atomicnumberQM[i],
	    qm->xQM[i][XX]/BOHR2NM,
	    qm->xQM[i][YY]/BOHR2NM,
	    qm->xQM[i][ZZ]/BOHR2NM);
#else
    fprintf(out,"%3d %10.7f  %10.7f  %10.7f\n",
	    qm->atomicnumberQM[i],
	    qm->xQM[i][XX]/BOHR2NM,
	    qm->xQM[i][YY]/BOHR2NM,
	    qm->xQM[i][ZZ]/BOHR2NM);
#endif
  }
  /* MM point charge data */
  if(QMMMrec->QMMMscheme!=eQMMMschemeoniom && mm->nrMMatoms){
    fprintf(out,"\n");
    for(i=0;i<mm->nrMMatoms;i++){
#ifdef GMX_DOUBLE
      fprintf(out,"%10.7lf  %10.7lf  %10.7lf %8.4lf\n",
	      mm->xMM[i][XX]/BOHR2NM,
	      mm->xMM[i][YY]/BOHR2NM,
	      mm->xMM[i][ZZ]/BOHR2NM,
	      mm->MMcharges[i]);
#else
      fprintf(out,"%10.7f  %10.7f  %10.7f %8.4f\n",
	      mm->xMM[i][XX]/BOHR2NM,
	      mm->xMM[i][YY]/BOHR2NM,
	      mm->xMM[i][ZZ]/BOHR2NM,
	      mm->MMcharges[i]);
#endif
    }
  }
  if(bSA) {/* put the SA coefficients at the end of the file */
#ifdef GMX_DOUBLE
    fprintf(out,"\n%10.8lf %10.8lf\n",
	    qm->SAstep*0.5/qm->SAsteps,
	    1-qm->SAstep*0.5/qm->SAsteps);
#else    
    fprintf(out,"\n%10.8f %10.8f\n",
	    qm->SAstep*0.5/qm->SAsteps,
	    1-qm->SAstep*0.5/qm->SAsteps);
#endif
    fprintf(stderr,"State Averaging level = %d/%d\n",qm->SAstep,qm->SAsteps);
  }
  fprintf(out,"\n");
  fclose(out);
}  /* write_gaussian_SH_input */

void write_gaussian_input(int step ,t_forcerec *fr, t_QMrec *qm, t_MMrec *mm)
{
  int
    i;
  t_QMMMrec
    *QMMMrec;
  FILE
    *out;
  
  QMMMrec = fr->qr;
  out = fopen("input.com","w");
  /* write the route */

  if(qm->QMmethod>=eQMmethodRHF)
    fprintf(out,"%s",
	    "%chk=input\n");
  else
    fprintf(out,"%s",
	    "%chk=se\n");
  if(qm->nQMcpus>1)
    fprintf(out,"%s%3d\n",
	    "%nprocshare=",qm->nQMcpus);
  fprintf(out,"%s%d\n",
	  "%mem=",qm->QMmem);
  /* use the modified links that include the LJ contribution at the QM level */
  if(qm->bTS||qm->bOPT){
    fprintf(out,"%s%s%s",
	    "%subst l701 ",qm->devel_dir,"/l701_LJ\n");
    fprintf(out,"%s%s%s",
	    "%subst l301 ",qm->devel_dir,"/l301_LJ\n");
  }
  else{
    fprintf(out,"%s%s%s",
	    "%subst l701 ",qm->devel_dir,"/l701\n");
    fprintf(out,"%s%s%s",
	    "%subst l301 ",qm->devel_dir,"/l301\n");
  }
  fprintf(out,"%s%s%s",
	  "%subst l9999 ",qm->devel_dir,"/l9999\n");
  if(step){
    fprintf(out,"%s",
	    "#T ");
  }else{
    fprintf(out,"%s",
	    "#P ");
  }
  if(qm->QMmethod==eQMmethodB3LYPLAN){
    fprintf(out," %s", 
	    "B3LYP/GEN Pseudo=Read");
  }
  else{
    fprintf(out," %s", 
	    eQMmethod_names[qm->QMmethod]);
    
    if(qm->QMmethod>=eQMmethodRHF){
      fprintf(out,"/%s",
	      eQMbasis_names[qm->QMbasis]);
      if(qm->QMmethod==eQMmethodCASSCF){
	/* in case of cas, how many electrons and orbitals do we need?
	 */
	fprintf(out,"(%d,%d)",
		qm->CASelectrons,qm->CASorbitals);
      }
    }
  }
  if(QMMMrec->QMMMscheme==eQMMMschemenormal){
    fprintf(out," %s",
	    "Charge ");
  }
  if (step || qm->QMmethod==eQMmethodCASSCF){
    /* fetch guess from checkpoint file, always for CASSCF */
    fprintf(out,"%s"," guess=read");
  }
  fprintf(out,"\nNosymm units=bohr\n");
  
  if(qm->bTS){
    fprintf(out,"OPT=(Redundant,TS,noeigentest,ModRedundant) Punch=(Coord,Derivatives) ");
  }
  else if (qm->bOPT){
    fprintf(out,"OPT=(Redundant,ModRedundant) Punch=(Coord,Derivatives) ");
  }
  else{
    fprintf(out,"FORCE Punch=(Derivatives) ");
  }
  fprintf(out,"iop(3/33=1)\n\n");
  fprintf(out, "input-file generated by gromacs\n\n");
  fprintf(out,"%2d%2d\n",qm->QMcharge,qm->multiplicity);
  for (i=0;i<qm->nrQMatoms;i++){
#ifdef GMX_DOUBLE
    fprintf(out,"%3d %10.7lf  %10.7lf  %10.7lf\n",
	    qm->atomicnumberQM[i],
	    qm->xQM[i][XX]/BOHR2NM,
	    qm->xQM[i][YY]/BOHR2NM,
	    qm->xQM[i][ZZ]/BOHR2NM);
#else
    fprintf(out,"%3d %10.7f  %10.7f  %10.7f\n",
	    qm->atomicnumberQM[i],
	    qm->xQM[i][XX]/BOHR2NM,
	    qm->xQM[i][YY]/BOHR2NM,
	    qm->xQM[i][ZZ]/BOHR2NM);
#endif
  }

  /* Pseudo Potential and ECP are included here if selected (MEthod suffix LAN) */
  if(qm->QMmethod==eQMmethodB3LYPLAN){
    fprintf(out,"\n");
    for(i=0;i<qm->nrQMatoms;i++){
      if(qm->atomicnumberQM[i]<21){
	fprintf(out,"%d ",i+1);
      }
    }
    fprintf(out,"\n%s\n****\n",eQMbasis_names[qm->QMbasis]);
    
    for(i=0;i<qm->nrQMatoms;i++){
      if(qm->atomicnumberQM[i]>21){
	fprintf(out,"%d ",i+1);
      }
    }
    fprintf(out,"\n%s\n****\n\n","lanl2dz");    
    
    for(i=0;i<qm->nrQMatoms;i++){
      if(qm->atomicnumberQM[i]>21){
	fprintf(out,"%d ",i+1);
      }
    }
    fprintf(out,"\n%s\n","lanl2dz");    
  }    
  
    
  
  /* MM point charge data */
  if(QMMMrec->QMMMscheme!=eQMMMschemeoniom && mm->nrMMatoms){
    fprintf(stderr,"nr mm atoms in gaussian.c = %d\n",mm->nrMMatoms);
    fprintf(out,"\n");
    if(qm->bTS||qm->bOPT){
      /* freeze the frontier QM atoms and Link atoms. This is
       * important only if a full QM subsystem optimization is done
       * with a frozen MM environmeent. For dynamics, or gromacs's own
       * optimization routines this is not important.
       */
      for(i=0;i<qm->nrQMatoms;i++){
	if(qm->frontatoms[i]){
	  fprintf(out,"%d F\n",i+1); /* counting from 1 */
	}
      }
      /* MM point charges include LJ parameters in case of QM optimization
       */
      for(i=0;i<mm->nrMMatoms;i++){
#ifdef GMX_DOUBLE
	fprintf(out,"%10.7lf  %10.7lf  %10.7lf %8.4lf 0.0 %10.7lf %10.7lf\n",
		mm->xMM[i][XX]/BOHR2NM,
		mm->xMM[i][YY]/BOHR2NM,
		mm->xMM[i][ZZ]/BOHR2NM,
		mm->MMcharges[i],
		mm->c6[i],mm->c12[i]);
#else
	fprintf(out,"%10.7f  %10.7f  %10.7f %8.4f 0.0 %10.7f %10.7f\n",
		mm->xMM[i][XX]/BOHR2NM,
		mm->xMM[i][YY]/BOHR2NM,
		mm->xMM[i][ZZ]/BOHR2NM,
		mm->MMcharges[i],
		mm->c6[i],mm->c12[i]);
#endif
      }
      fprintf(out,"\n");
    }
    else{
      for(i=0;i<mm->nrMMatoms;i++){
#ifdef GMX_DOUBLE
	fprintf(out,"%10.7lf  %10.7lf  %10.7lf %8.4lf\n",
		mm->xMM[i][XX]/BOHR2NM,
		mm->xMM[i][YY]/BOHR2NM,
		mm->xMM[i][ZZ]/BOHR2NM,
		mm->MMcharges[i]);
#else
	fprintf(out,"%10.7f  %10.7f  %10.7f %8.4f\n",
		mm->xMM[i][XX]/BOHR2NM,
		mm->xMM[i][YY]/BOHR2NM,
		mm->xMM[i][ZZ]/BOHR2NM,
		mm->MMcharges[i]);
#endif
      }
    }
  }
  fprintf(out,"\n");
  

  fclose(out);

}  /* write_gaussian_input */

real read_gaussian_output(rvec QMgrad[],rvec MMgrad[],int step,
			  t_QMrec *qm, t_MMrec *mm)
{
  int
    i,j,atnum;
  char
    buf[300];
  real
    QMener;
  FILE
    *in;
  
  in=fopen("fort.7","r");



  /* in case of an optimization, the coordinates are printed in the
   * fort.7 file first, followed by the energy, coordinates and (if
   * required) the CI eigenvectors.
   */
  if(qm->bTS||qm->bOPT){
    for(i=0;i<qm->nrQMatoms;i++){
      if( NULL == fgets(buf,300,in))
      {
	  gmx_fatal(FARGS,"Error reading Gaussian output - not enough atom lines?");
      }

#ifdef GMX_DOUBLE
      sscanf(buf,"%d %lf %lf %lf\n",
	     &atnum,
	     &qm->xQM[i][XX],
	     &qm->xQM[i][YY],
	     &qm->xQM[i][ZZ]);
#else
      sscanf(buf,"%d %f %f %f\n",
	     &atnum,
	     &qm->xQM[i][XX],
	     &qm->xQM[i][YY],
	     &qm->xQM[i][ZZ]);
#endif     
      for(j=0;j<DIM;j++){
	qm->xQM[i][j]*=BOHR2NM;
      }
    }
  }
  /* the next line is the energy and in the case of CAS, the energy
   * difference between the two states.
   */
  if(NULL == fgets(buf,300,in))
  {
      gmx_fatal(FARGS,"Error reading Gaussian output");
  }

#ifdef GMX_DOUBLE
  sscanf(buf,"%lf\n",&QMener);
#else
  sscanf(buf,"%f\n", &QMener);
#endif
  /* next lines contain the gradients of the QM atoms */
  for(i=0;i<qm->nrQMatoms;i++){
    if(NULL == fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }
#ifdef GMX_DOUBLE
    sscanf(buf,"%lf %lf %lf\n",
	   &QMgrad[i][XX],
	   &QMgrad[i][YY],
	   &QMgrad[i][ZZ]);
#else
    sscanf(buf,"%f %f %f\n",
	   &QMgrad[i][XX],
	   &QMgrad[i][YY],
	   &QMgrad[i][ZZ]);
#endif     
  }
  /* the next lines are the gradients of the MM atoms */
  if(qm->QMmethod>=eQMmethodRHF){  
    for(i=0;i<mm->nrMMatoms;i++){
      if(NULL==fgets(buf,300,in))
      {
          gmx_fatal(FARGS,"Error reading Gaussian output");
      }
#ifdef GMX_DOUBLE
      sscanf(buf,"%lf %lf %lf\n",
	     &MMgrad[i][XX],
	     &MMgrad[i][YY],
	     &MMgrad[i][ZZ]);
#else
      sscanf(buf,"%f %f %f\n",
	     &MMgrad[i][XX],
	     &MMgrad[i][YY],
	     &MMgrad[i][ZZ]);
#endif	
    }
  }
  fclose(in);
  return(QMener);  
}

real read_gaussian_SH_output(rvec QMgrad[],rvec MMgrad[],int step,
			     gmx_bool swapped,t_QMrec *qm, t_MMrec *mm)
{
  int
    i;
  char
    buf[300];
  real
    QMener,DeltaE;
  FILE
    *in;
  
  in=fopen("fort.7","r");
  /* first line is the energy and in the case of CAS, the energy
   * difference between the two states.
   */
  if(NULL == fgets(buf,300,in))
  {
      gmx_fatal(FARGS,"Error reading Gaussian output");
  }

#ifdef GMX_DOUBLE
  sscanf(buf,"%lf %lf\n",&QMener,&DeltaE);
#else
  sscanf(buf,"%f %f\n",  &QMener,&DeltaE);
#endif
  
  /* switch on/off the State Averaging */
  
  if(DeltaE > qm->SAoff){
    if (qm->SAstep > 0){
      qm->SAstep--;
    }
  }
  else if (DeltaE < qm->SAon || (qm->SAstep > 0)){
    if (qm->SAstep < qm->SAsteps){
      qm->SAstep++;
    }
  }
  
  /* for debugging: */
  fprintf(stderr,"Gap = %5f,SA = %3d\n",DeltaE,(qm->SAstep>0));
  /* next lines contain the gradients of the QM atoms */
  for(i=0;i<qm->nrQMatoms;i++){
    if(NULL==fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }

#ifdef GMX_DOUBLE
    sscanf(buf,"%lf %lf %lf\n",
	   &QMgrad[i][XX],
	   &QMgrad[i][YY],
	   &QMgrad[i][ZZ]);
#else
    sscanf(buf,"%f %f %f\n",
	   &QMgrad[i][XX],
	   &QMgrad[i][YY],
	   &QMgrad[i][ZZ]);
#endif     
  }
  /* the next lines, are the gradients of the MM atoms */
  
  for(i=0;i<mm->nrMMatoms;i++){
    if(NULL==fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }
#ifdef GMX_DOUBLE
    sscanf(buf,"%lf %lf %lf\n",
	   &MMgrad[i][XX],
	   &MMgrad[i][YY],
	   &MMgrad[i][ZZ]);
#else
    sscanf(buf,"%f %f %f\n",
	   &MMgrad[i][XX],
	   &MMgrad[i][YY],
	   &MMgrad[i][ZZ]);
#endif	
  }
  
  /* the next line contains the two CI eigenvector elements */
  if(NULL==fgets(buf,300,in))
  {
      gmx_fatal(FARGS,"Error reading Gaussian output");
  }
  if(!step){
    sscanf(buf,"%d",&qm->CIdim);
    snew(qm->CIvec1,qm->CIdim);
    snew(qm->CIvec1old,qm->CIdim);
    snew(qm->CIvec2,qm->CIdim);
    snew(qm->CIvec2old,qm->CIdim);
  } else {
    /* before reading in the new current CI vectors, copy the current
     * CI vector into the old one.
     */
    for(i=0;i<qm->CIdim;i++){
      qm->CIvec1old[i] = qm->CIvec1[i];
      qm->CIvec2old[i] = qm->CIvec2[i];
    }
  }
  /* first vector */
  for(i=0;i<qm->CIdim;i++){
    if(NULL==fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }
#ifdef GMX_DOUBLE
    sscanf(buf,"%lf\n",&qm->CIvec1[i]);
#else
    sscanf(buf,"%f\n", &qm->CIvec1[i]);   
#endif
  }
  /* second vector */
  for(i=0;i<qm->CIdim;i++){
    if(NULL==fgets(buf,300,in))
    {
	gmx_fatal(FARGS,"Error reading Gaussian output");
    }
#ifdef GMX_DOUBLE
    sscanf(buf,"%lf\n",&qm->CIvec2[i]);
#else
    sscanf(buf,"%f\n", &qm->CIvec2[i]);   
#endif
  }
  fclose(in);
  return(QMener);  
}

real inproduct(real *a, real *b, int n)
{
  int
    i;
  real
    dot=0.0;
  
  /* computes the inner product between two vectors (a.b), both of
   * which have length n.
   */  
  for(i=0;i<n;i++){
    dot+=a[i]*b[i];
  }
  return(dot);
}

int hop(int step, t_QMrec *qm)
{
  int
    swap = 0;
  real
    d11=0.0,d12=0.0,d21=0.0,d22=0.0;
  
  /* calculates the inproduct between the current Ci vector and the
   * previous CI vector. A diabatic hop will be made if d12 and d21
   * are much bigger than d11 and d22. In that case hop returns true,
   * otherwise it returns false.
   */  
  if(step){ /* only go on if more than one step has been done */
    d11 = inproduct(qm->CIvec1,qm->CIvec1old,qm->CIdim);
    d12 = inproduct(qm->CIvec1,qm->CIvec2old,qm->CIdim);
    d21 = inproduct(qm->CIvec2,qm->CIvec1old,qm->CIdim);
    d22 = inproduct(qm->CIvec2,qm->CIvec2old,qm->CIdim);
  }
  fprintf(stderr,"-------------------\n");
  fprintf(stderr,"d11 = %13.8f\n",d11);
  fprintf(stderr,"d12 = %13.8f\n",d12);
  fprintf(stderr,"d21 = %13.8f\n",d21);
  fprintf(stderr,"d22 = %13.8f\n",d22);
  fprintf(stderr,"-------------------\n");
  
  if((fabs(d12)>0.5)&&(fabs(d21)>0.5))
    swap = 1;
  
  return(swap);
}

void do_gaussian(int step,char *exe)
{
  char
    buf[100];

  /* make the call to the gaussian binary through system()
   * The location of the binary will be picked up from the 
   * environment using getenv().
   */
  if(step) /* hack to prevent long inputfiles */
    sprintf(buf,"%s < %s > %s",
	    exe,
	    "input.com",
	    "input.log");
  else
    sprintf(buf,"%s < %s > %s",
	    exe,
            "input.com",
	    "input.log");
  fprintf(stderr,"Calling '%s'\n",buf);
#ifdef GMX_NO_SYSTEM
  printf("Warning-- No calls to system(3) supported on this platform.");
  gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#else
  if ( system(buf) != 0 )
    gmx_fatal(FARGS,"Call to '%s' failed\n",buf);
#endif
}

real call_gaussian(t_commrec *cr,  t_forcerec *fr, 
		   t_QMrec *qm, t_MMrec *mm, rvec f[], rvec fshift[])
{
  /* normal gaussian jobs */
  static int
    step=0;
  int
    i,j;
  real
    QMener=0.0;
  rvec
    *QMgrad,*MMgrad;
  char
    *exe;
  
  snew(exe,30);
  sprintf(exe,"%s/%s",qm->gauss_dir,qm->gauss_exe);
  snew(QMgrad,qm->nrQMatoms);
  snew(MMgrad,mm->nrMMatoms);

  write_gaussian_input(step,fr,qm,mm);
  do_gaussian(step,exe);
  QMener = read_gaussian_output(QMgrad,MMgrad,step,qm,mm);
  /* put the QMMM forces in the force array and to the fshift
   */
  for(i=0;i<qm->nrQMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
      fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
    }
  }
  for(i=0;i<mm->nrMMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i+qm->nrQMatoms][j]      = HARTREE_BOHR2MD*MMgrad[i][j];      
      fshift[i+qm->nrQMatoms][j] = HARTREE_BOHR2MD*MMgrad[i][j];
    }
  }
  QMener = QMener*HARTREE2KJ*AVOGADRO;
  step++;
  free(exe);
  return(QMener);

} /* call_gaussian */

real call_gaussian_SH(t_commrec *cr, t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, 
		      rvec f[], rvec fshift[])
{ 
  /* a gaussian call routine intended for doing diabatic surface
   * "sliding". See the manual for the theoretical background of this
   * TSH method.  
   */
  static int
    step=0;
  int
    state,i,j;
  real
    QMener=0.0;
  static  gmx_bool
    swapped=FALSE; /* handle for identifying the current PES */
  gmx_bool
    swap=FALSE; /* the actual swap */
  rvec
    *QMgrad,*MMgrad;
  char
    *buf;
  char
    *exe;
  
  snew(exe,30);
  sprintf(exe,"%s/%s",qm->gauss_dir,qm->gauss_exe);
  /* hack to do ground state simulations */
  if(!step){
    snew(buf,20);
    buf = getenv("STATE");
    if (buf)
      sscanf(buf,"%d",&state);
    else
      state=2;
    if(state==1)
      swapped=TRUE;
  }
  /* end of hack */


  /* copy the QMMMrec pointer */
  snew(QMgrad,qm->nrQMatoms);
  snew(MMgrad,mm->nrMMatoms);
  /* at step 0 there should be no SA */
  /*  if(!step)
   * qr->bSA=FALSE;*/
  /* temporray set to step + 1, since there is a chk start */
  write_gaussian_SH_input(step,swapped,fr,qm,mm);

  do_gaussian(step,exe);
  QMener = read_gaussian_SH_output(QMgrad,MMgrad,step,swapped,qm,mm);

  /* check for a surface hop. Only possible if we were already state
   * averaging.
   */
  if(qm->SAstep>0){
    if(!swapped){
      swap    = (step && hop(step,qm));
      swapped = swap;
    } 
    else { /* already on the other surface, so check if we go back */
      swap    = (step && hop(step,qm));
      swapped =!swap; /* so swapped shoud be false again */
    }
    if (swap){/* change surface, so do another call */
      write_gaussian_SH_input(step,swapped,fr,qm,mm);
      do_gaussian(step,exe);
      QMener = read_gaussian_SH_output(QMgrad,MMgrad,step,swapped,qm,mm);
    }
  }
  /* add the QMMM forces to the gmx force array and fshift
   */
  for(i=0;i<qm->nrQMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i][j]      = HARTREE_BOHR2MD*QMgrad[i][j];
      fshift[i][j] = HARTREE_BOHR2MD*QMgrad[i][j];
    }
  }
  for(i=0;i<mm->nrMMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i+qm->nrQMatoms][j]      = HARTREE_BOHR2MD*MMgrad[i][j];
      fshift[i+qm->nrQMatoms][j] = HARTREE_BOHR2MD*MMgrad[i][j];
    }
  }
  QMener = QMener*HARTREE2KJ*AVOGADRO;
  fprintf(stderr,"step %5d, SA = %5d, swap = %5d\n",
	  step,(qm->SAstep>0),swapped);
  step++;
  free(exe);
  return(QMener);

} /* call_gaussian_SH */
    
/* end of gaussian sub routines */

#else
int
gmx_qmmm_gaussian_empty;
#endif

