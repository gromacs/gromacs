/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */

/* This file is completely threadsafe - keep it that way! */
#include "xmdrun.h"
#include "smalloc.h"
#include "assert.h"
#include "names.h"
	
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
  int         pt1,pt2,a1,a2;
  bool        bS1,bS2;
  t_iatom     *ia;

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
    fatal_error(0,"Your number of shells %d is not equal to the number of shells %d",nsi,n[eptShell]);
  
  for(i=0; (i<eptNR); i++)
    if (n[i]!=0)
      fprintf(log,"There are: %d %s\n",n[i],ptype_str[i]);
  
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
    ia=idef->il[F_BONDS].iatoms;
    for(i=0; (i<idef->il[F_BONDS].nr); ) {
      type  = ia[0];
      ftype = idef->functype[type];
      nra   = interaction_function[ftype].nratoms;
      
      /* Check whether we have a bond */
      
      if (md->ptype[ia[1]] == eptShell) {
	a1 = ia[1];
	a2 = ia[2];
      }
      else if (md->ptype[ia[2]] == eptShell) {
	a1 = ia[2];
	a2 = ia[1];
      }
      else {
	i  += nra+1;
	ia += nra+1;
	continue;
      }
      /* Check whether one of the particles is a shell... */
      nsi = shell_index[a1-start];
      if ((nsi < 0) || (nsi >= *nshell))
	fatal_error(0,"nsi is %d should be within 0 - %d. a1 = %d",
		    nsi,*nshell,a1);
      if (shell[nsi].shell == NO_ATID) {
	shell[nsi].shell = a1;
	ns ++;
      }
      else if (shell[nsi].shell != a1)
	fatal_error(0,"What is this?");
      
      if      (shell[nsi].nucl1 == NO_ATID)
	shell[nsi].nucl1 = a2;
      else if (shell[nsi].nucl2 == NO_ATID)
	shell[nsi].nucl2 = a2;
      else if (shell[nsi].nucl3 == NO_ATID)
	shell[nsi].nucl3 = a2;
      else {
	pr_shell(log,ns,shell);
	fatal_error(0,"Can not handle more than three bonds per shell\n");
      }
      shell[nsi].k    += idef->iparams[type].harmonic.krA;
      shell[nsi].nnucl++;
      
      ia += nra+1;
      i  += nra+1;
    }
    ia=idef->il[F_CUBICBONDS].iatoms;
    for(i=0; (i<idef->il[F_CUBICBONDS].nr); ) {
      type  = ia[0];
      ftype = idef->functype[type];
      nra   = interaction_function[ftype].nratoms;
      
      /* Check whether we have a bond */
      
      if (md->ptype[ia[1]] == eptShell) {
	a1 = ia[1];
	a2 = ia[2];
      }
      else if (md->ptype[ia[2]] == eptShell) {
	a1 = ia[2];
	a2 = ia[1];
      }
      else {
	i  += nra+1;
	ia += nra+1;
	continue;
      }
      /* Check whether one of the particles is a shell... */
      nsi = shell_index[a1-start];
      if ((nsi < 0) || (nsi >= *nshell))
	fatal_error(0,"nsi is %d should be within 0 - %d. a1 = %d",
		    nsi,*nshell,a1);
      if (shell[nsi].shell == NO_ATID) {
	shell[nsi].shell = a1;
	ns ++;
      }
      else if (shell[nsi].shell != a1)
	fatal_error(0,"What is this?");
      
      if      (shell[nsi].nucl1 == NO_ATID)
	shell[nsi].nucl1 = a2;
      else if (shell[nsi].nucl2 == NO_ATID)
	shell[nsi].nucl2 = a2;
      else if (shell[nsi].nucl3 == NO_ATID)
	shell[nsi].nucl3 = a2;
      else {
	pr_shell(log,ns,shell);
	fatal_error(0,"Can not handle more than three bonds per shell\n");
      }
      shell[nsi].k    += idef->iparams[type].cubic.kb;
      shell[nsi].nnucl++;
      
      ia += nra+1;
      i  += nra+1;
    }
    ia = idef->il[F_WPOL].iatoms;
    for(i=0; (i<idef->il[F_WPOL].nr); ) {
      type  = ia[0];
      ftype = idef->functype[type];
      nra   = interaction_function[ftype].nratoms;
      
      a1    = ia[1]+4;  /* Shell */
      a2    = ia[1]+3;  /* Dummy */
      
      /* Check whether one of the particles is a shell... */
      nsi = shell_index[a1-start];
      if ((nsi < 0) || (nsi >= *nshell))
	fatal_error(0,"nsi is %d should be within 0 - %d. a1 = %d",
		    nsi,*nshell,a1);
      if (shell[nsi].shell == NO_ATID) {
	shell[nsi].shell = a1;
	ns ++;
      }
      else if (shell[nsi].shell != a1)
	fatal_error(0,"What is this? shell=%d, a1=%d",shell[nsi].shell,a1);
      
      shell[nsi].nucl1 = a2;
      shell[nsi].k     = (idef->iparams[type].wpol.kx+
			  idef->iparams[type].wpol.ky+
			  idef->iparams[type].wpol.kz)/3.0;
      shell[nsi].nnucl++;
      
      ia += nra+1;
      i  += nra+1;
    }
    /* Verify whether it's all correct */
    if (ns != *nshell)
      fatal_error(0,"Something weird with shells. They may not be bonded to something");

    for(i=0; (i<ns); i++)
      shell[i].k_1 = 1.0/shell[i].k;

    if (debug)
      pr_shell(debug,ns,shell);
  }
  return shell;
}

