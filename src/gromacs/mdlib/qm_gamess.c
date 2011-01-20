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

#ifdef GMX_QMMM_GAMESS

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


/* QMMM sub routines */
/* mopac interface routines */


#ifndef F77_FUNC
#define F77_FUNC(name,NAME) name ## _
#endif


void 
F77_FUNC(inigms,IMIGMS)(void);

void 
F77_FUNC(endgms,ENDGMS)(void);

void 
F77_FUNC(grads,GRADS)(int *nrqmat,real *qmcrd,int *nrmmat, real *mmchrg, 
                      real *mmcrd, real *qmgrad,real *mmgrad, real *energy);



void init_gamess(t_commrec *cr,t_QMrec *qm, t_MMrec *mm){
  /* it works hopelessly complicated :-)
   * first a file is written. Then the standard gamess input/output
   * routine is called (no system()!) to set up all fortran arrays. 
   * this routine writes a punch file, like in a normal gamess run.
   * via this punch file the other games routines, needed for gradient
   * and energy evaluations are called. This setup works fine for 
   * dynamics simulations. 7-6-2002 (London)
   */
  int 
    i,j,rank;
  FILE
    *out;
  char
    periodic_system[37][3]={"XX","H ","He","Li","Be","B ","C ","N ",
			    "O ","F ","Ne","Na","Mg","Al","Si","P ",
			    "S ","Cl","Ar","K ","Ca","Sc","Ti","V ",
			    "Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga",
			    "Ge","As","Se","Br","Kr"};
  
  if (PAR(cr)){

    if MASTER(cr){
      out=fopen("FOR009","w");
      /* of these options I am not completely sure....  the overall
       * preformance on more than 4 cpu's is rather poor at the moment.  
       */
      fprintf(out,"memory 48000000\nPARALLEL IOMODE SCREENED\n");
      fprintf(out,"ELEC %d\nMULT %d\nSUPER ON\nNOSYM\nGEOMETRY ANGSTROM\n",
	      qm->nelectrons,qm->multiplicity);
      for (i=0;i<qm->nrQMatoms;i++){
#ifdef DOUBLE
	fprintf(out,"%10.7lf  %10.7lf  %10.7lf  %5.3lf  %2s\n",
		i/2.,
		i/3.,
		i/4.,  
		qm->atomicnumberQM[i]*1.0,
		periodic_system[qm->atomicnumberQM[i]]);
#else
	fprintf(out,"%10.7f  %10.7f  %10.7f  %5.3f  %2s\n",
		i/2.,
		i/3.,
		i/4.,  
		qm->atomicnumberQM[i]*1.0,
		periodic_system[qm->atomicnumberQM[i]]);
#endif
      }
      if(mm->nrMMatoms){
	for (j=i;j<i+2;j++){
#ifdef DOUBLE
	  fprintf(out,"%10.7lf  %10.7lf  %10.7lf  %5.3lf  BQ\n",
		  j/5.,
		  j/6.,
		  j/7.,
		  1.0);  
#else
	  fprintf(out,"%10.7f  %10.7f  %10.7f  %5.3f  BQ\n",
		  j/5.,
		  j/6.,
		  j/7.,
		  2.0);  
#endif
	}
      }
      if(!qm->bTS)
	fprintf(out,"END\nBASIS %s\nRUNTYPE GRADIENT\nSCFTYPE %s\n",
		eQMbasis_names[qm->QMbasis],
		eQMmethod_names[qm->QMmethod]); /* see enum.h */
      else
	fprintf(out,"END\nBASIS %s\nRUNTYPE SADDLE\nSCFTYPE %s\n",
		eQMbasis_names[qm->QMbasis],
		eQMmethod_names[qm->QMmethod]); /* see enum.h */
      fclose(out);
    }
    gmx_barrier(cr);
    F77_FUNC(inigms,IMIGMS)();
  }
  else{ /* normal serial run */
    
    out=fopen("FOR009","w");
    /* of these options I am not completely sure....  the overall
     * preformance on more than 4 cpu's is rather poor at the moment.  
     */
    fprintf(out,"ELEC %d\nMULT %d\nSUPER ON\nNOSYM\nGEOMETRY ANGSTROM\n",
	    qm->nelectrons,qm->multiplicity);
    for (i=0;i<qm->nrQMatoms;i++){
#ifdef DOUBLE
      fprintf(out,"%10.7lf  %10.7lf  %10.7lf  %5.3lf  %2s\n",
	      i/2.,
	      i/3.,
	      i/4.,  
	      qm->atomicnumberQM[i]*1.0,
	      periodic_system[qm->atomicnumberQM[i]]);
#else
      fprintf(out,"%10.7f  %10.7f  %10.7f  %5.3f  %2s\n",
	      i/2.,
	      i/3.,
	      i/4.,  
	      qm->atomicnumberQM[i]*1.0,
	      periodic_system[qm->atomicnumberQM[i]]);
#endif
    }
    if(mm->nrMMatoms){
      for (j=i;j<i+2;j++){
#ifdef DOUBLE
	fprintf(out,"%10.7lf  %10.7lf  %10.7lf  %5.3lf  BQ\n",
		j/5.,
		j/6.,
		j/7.,
		1.0);  
#else
	fprintf(out,"%10.7f  %10.7f  %10.7f  %5.3f  BQ\n",
		j/5.,
		j/6.,
		j/7.,
		2.0);  
#endif
      }
    }
    if(!qm->bTS)
      fprintf(out,"END\nBASIS %s\nRUNTYPE GRADIENT\nSCFTYPE %s\n",
	      eQMbasis_names[qm->QMbasis],
	      eQMmethod_names[qm->QMmethod]); /* see enum.h */
    else
      fprintf(out,"END\nBASIS %s\nRUNTYPE SADDLE\nSCFTYPE %s\n",
	      eQMbasis_names[qm->QMbasis],
	      eQMmethod_names[qm->QMmethod]); /* see enum.h */
    F77_FUNC(inigms,IMIGMS)();
  }  
}

real call_gamess(t_commrec *cr, t_forcerec *fr, t_QMrec *qm, t_MMrec *mm, 
		 rvec f[], rvec fshift[])
{
  /* do the actual QMMM calculation using GAMESS-UK. In this
   * implementation (3-2001) a system call is made to the GAMESS-UK
   * binary. Now we are working to get the electron integral, SCF, and
   * gradient routines linked directly 
   */
  int 
    i,j,rank;
  real
    QMener=0.0,*qmgrad,*mmgrad,*mmcrd,*qmcrd,energy;
  t_QMMMrec
    *qr;

  /* copy the QMMMrec pointer */
  qr = fr->qr;
  snew(qmcrd, 3*(qm->nrQMatoms));
  snew(mmcrd,3*(mm->nrMMatoms));
  snew(qmgrad,3*(qm->nrQMatoms));
  snew(mmgrad,3*(mm->nrMMatoms));
  
  /* copy the data from qr into the arrays that are going to be used
   * in the fortran routines of gamess
   */
  for(i=0;i<qm->nrQMatoms;i++){
    for (j=0;j<DIM;j++){
      qmcrd[DIM*i+j] = 1/BOHR2NM*qm->xQM[i][j];
    }
  }
  for(i=0;i<mm->nrMMatoms;i++){
    for (j=0;j<DIM;j++){
      mmcrd[DIM*i+j] = 1/BOHR2NM*mm->xMM[i][j];
    }
  }
  for (i=0;i<3*qm->nrQMatoms;i+=3){
    fprintf(stderr,"%8.5f, %8.5f, %8.5f\n",
	    qmcrd[i],
	    qmcrd[i+1],
	    qmcrd[i+2]);
  }

  F77_FUNC(grads,GRADS)(&qm->nrQMatoms,qmcrd,&mm->nrMMatoms,mm->MMcharges,
                        mmcrd,qmgrad,mmgrad,&energy);

  for(i=0;i<qm->nrQMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i][j]      = HARTREE_BOHR2MD*qmgrad[3*i+j];
      fshift[i][j] = HARTREE_BOHR2MD*qmgrad[3*i+j];
    }
  }
  for(i=0;i<mm->nrMMatoms;i++){
    for(j=0;j<DIM;j++){
      f[i][j]      = HARTREE_BOHR2MD*mmgrad[3*i+j];
      fshift[i][j] = HARTREE_BOHR2MD*mmgrad[3*i+j];
    }
  }
  /* convert a.u to kJ/mol */
  QMener=energy*HARTREE2KJ*AVOGADRO;
  return(QMener);
}

#else
int
gmx_qmmm_gamess_empty;
#endif

