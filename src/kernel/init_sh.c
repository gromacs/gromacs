/*
 * $Id$
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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
/* This file is completely threadsafe - keep it that way! */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "xmdrun.h"
#include "smalloc.h"
#include "names.h"
#include "vec.h"
#include "physics.h"
#include "copyrite.h"
	
static void pr_shell(FILE *log,int ns,t_shell s[])
{
  int i;
  
  fprintf(log,"SHELL DATA\n");
  fprintf(log,"%5s  %8s  %5s  %5s  %5s\n",
	  "Shell","Force k","Nucl1","Nucl2","Nucl3");
  for(i=0; (i<ns); i++) {
    fprintf(log,"%5d  %8.3f  %5d",s[i].shell,1.0/s[i].k_1,s[i].nucl1);
    if (s[i].nnucl == 2)
      fprintf(log,"  %5d\n",s[i].nucl2);
    else if (s[i].nnucl == 3)
      fprintf(log,"  %5d  %5d\n",s[i].nucl2,s[i].nucl3);
    else
      fprintf(log,"\n");
  }
}

t_shell *init_shells(FILE *log,int start,int homenr,
		     t_idef *idef,t_mdatoms *md,int *nshell)
{
  t_shell     *shell=NULL;
  int         *shell_index;
  int         n[eptNR],ns,nsi;
  int         i,j,type,ftype,nra;
  real        qS,alpha;
  int         aS,aN=0; /* Shell and nucleus */
  int         bondtypes[] = { F_BONDS, F_CUBICBONDS, F_POLARIZATION, F_WATER_POL };
#define NBT asize(bondtypes)
  bool        bS1,bS2;
  t_iatom     *ia;

  /* Count number of shells, and find their indices */
  for(i=0; (i<eptNR); i++)
    n[i]=0;
  snew(shell_index,homenr);
  nsi = 0;
  for(i=start; (i<start+homenr); i++) {
    n[md->ptype[i]]++;
    if (md->ptype[i] == eptShell)
      shell_index[i-start] = nsi++;
  }
  if (nsi != n[eptShell])
    gmx_fatal(FARGS,"Your number of shells %d is not equal to the number of shells %d",
		nsi,n[eptShell]);

  /* Print the number of each particle type */  
  for(i=0; (i<eptNR); i++)
    if (n[i]!=0)
      fprintf(log,"There are: %d %ss\n",n[i],ptype_str[i]);
  
  ns      = n[eptShell];
  *nshell = ns;
  if (ns > 0) {
    snew(shell,ns);
  
    /* Initiate the shell structures */    
    for(i=0; (i<ns); i++) {
      shell[i].shell=NO_ATID;
      shell[i].nucl1=NO_ATID;
      shell[i].nucl2=NO_ATID;
      shell[i].nucl3=NO_ATID;
      shell[i].nnucl=0;
      shell[i].k_1=0;
      shell[i].k=0;
    }
    
    /* Now fill the structures */
    ns=0;
    for(j=0; (j<NBT); j++) {
      ia=idef->il[bondtypes[j]].iatoms;
      for(i=0; (i<idef->il[bondtypes[j]].nr); ) {
	type  = ia[0];
	ftype = idef->functype[type];
	nra   = interaction_function[ftype].nratoms;
	
	/* Check whether we have a bond with a shell */
	aS = NO_ATID;
	
	switch (bondtypes[j]) {
	case F_BONDS:
	case F_CUBICBONDS:
	case F_POLARIZATION:
	  if (md->ptype[ia[1]] == eptShell) {
	    aS = ia[1];
	    aN = ia[2];
	  }
	  else if (md->ptype[ia[2]] == eptShell) {
	    aS = ia[2];
	    aN = ia[1];
	  }
	  break;
	case F_WATER_POL:
	  aN    = ia[4];  /* Dummy */
	  aS    = ia[5];  /* Shell */
	  break;
	default:
	  gmx_fatal(FARGS,"Death Horror: %s, %d",__FILE__,__LINE__);
	}

	if (aS != NO_ATID) {	  
	  qS = md->chargeA[aS];
	  
	  /* Check whether one of the particles is a shell... */
	  nsi = shell_index[aS-start];
	  if ((nsi < 0) || (nsi >= *nshell))
	    gmx_fatal(FARGS,"nsi is %d should be within 0 - %d. aS = %d",
			nsi,*nshell,aS);
	  if (shell[nsi].shell == NO_ATID) {
	    shell[nsi].shell = aS;
	    ns ++;
	  }
	  else if (shell[nsi].shell != aS)
	    gmx_fatal(FARGS,"Weird stuff in %s, %d",__FILE__,__LINE__);
	  
	  if      (shell[nsi].nucl1 == NO_ATID)
	    shell[nsi].nucl1 = aN;
	  else if (shell[nsi].nucl2 == NO_ATID)
	    shell[nsi].nucl2 = aN;
	  else if (shell[nsi].nucl3 == NO_ATID)
	    shell[nsi].nucl3 = aN;
	  else {
	    pr_shell(log,ns,shell);
	    gmx_fatal(FARGS,"Can not handle more than three bonds per shell\n");
	  }
	  switch (bondtypes[j]) {
	  case F_BONDS:
	    shell[nsi].k    += idef->iparams[type].harmonic.krA;
	    break;
	  case F_CUBICBONDS:
	    shell[nsi].k    += idef->iparams[type].cubic.kb;
	    break;
	  case F_POLARIZATION:
	    if (md->nChargePerturbed && qS != md->chargeB[aS])
	      gmx_fatal(FARGS,"polarize can not be used with qA != qB");
	    shell[nsi].k    += sqr(qS)*ONE_4PI_EPS0/
	      idef->iparams[type].polarize.alpha;
	  break;
	  case F_WATER_POL:
	    if (md->nChargePerturbed && qS != md->chargeB[aS])
	      gmx_fatal(FARGS,"water_pol can not be used with qA != qB");
	    alpha          = (idef->iparams[type].wpol.al_x+
			      idef->iparams[type].wpol.al_y+
			      idef->iparams[type].wpol.al_z)/3.0;
	    shell[nsi].k  += sqr(qS)*ONE_4PI_EPS0/alpha;
	    break;
	  default:
	    gmx_fatal(FARGS,"Death Horror: %s, %d",__FILE__,__LINE__);
	  }
	  shell[nsi].nnucl++;
	}
	ia += nra+1;
	i  += nra+1;
      }
    }

    /* Verify whether it's all correct */
    if (ns != *nshell)
      gmx_fatal(FARGS,"Something weird with shells. They may not be bonded to something");

    for(i=0; (i<ns); i++)
      shell[i].k_1 = 1.0/shell[i].k;

    if (debug)
      pr_shell(debug,ns,shell);
  }
  
  return shell;
}

int count_flexible_constraints(FILE* log,t_forcerec *fr,t_idef *idef)
{
  int nflexcon,i;
  
  nflexcon = 0;
  
  for(i=0; i<idef->il[F_CONSTR].nr; i+=3)
    if (idef->iparams[idef->il[F_CONSTR].iatoms[i]].constr.dA == 0 &&
	idef->iparams[idef->il[F_CONSTR].iatoms[i]].constr.dB == 0)
      nflexcon++;
  
  if (nflexcon > 0) {
    fprintf(log,"There are %d flexible constraints\n",nflexcon);
    if (fr->fc_stepsize == 0) {
      fprintf(log,"WARNING: step size for flexible constraining = 0\n"
	          "         All flexible constraints will be rigid.\n"
	          "         Will try to keep all flexible constraints at their original length,\n"
	          "         but the lengths may exhibit some drift.\n\n");
      nflexcon = 0;
    } else {
      please_cite(log,"Hess2002");
    }
  }
  
  return nflexcon;
}
