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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "smalloc.h"
#include "macros.h"
#include "force.h"
#include "nonbonded.h"
#include "invblock.h"
#include "confio.h"
#include "nsb.h"
#include "names.h"
#include "network.h"
#include "wnblist.h"
#include "pbc.h"
#include "ns.h"
#include "nrnb.h"
#include "bondf.h"
#include "mshift.h"
#include "txtdump.h"
#include "ewald_util.h"
#include "shift_util.h"
#include "pppm.h"
#include "ewald.h"
#include "pme.h"
#include "mdrun.h"

t_forcerec *mk_forcerec(void)
{
  t_forcerec *fr;
  
  snew(fr,1);
  
  return fr;
}

#ifdef DEBUG
static void pr_nbfp(FILE *fp,real *nbfp,bool bBHAM,int atnr)
{
  int i,j;
  
  for(i=0; (i<atnr); i++) {
    for(j=0; (j<atnr); j++) {
      fprintf(fp,"%2d - %2d",i,j);
      if (bBHAM)
	fprintf(fp,"  a=%10g, b=%10g, c=%10g\n",BHAMA(nbfp,atnr,i,j),
		BHAMB(nbfp,atnr,i,j),BHAMC(nbfp,atnr,i,j));
      else
	fprintf(fp,"  c6=%10g, c12=%10g\n",C6(nbfp,atnr,i,j),
		C12(nbfp,atnr,i,j));
    }
  }
}
#endif

static real *mk_nbfp(const t_idef *idef,bool bBHAM)
{
  real *nbfp;
  int  i,j,k,atnr;
  
  atnr=idef->atnr;
  if (bBHAM) {
    snew(nbfp,3*atnr*atnr);
    for(i=k=0; (i<atnr); i++) {
      for(j=0; (j<atnr); j++,k++) {
	BHAMA(nbfp,atnr,i,j) = idef->iparams[k].bham.a;
	BHAMB(nbfp,atnr,i,j) = idef->iparams[k].bham.b;
	BHAMC(nbfp,atnr,i,j) = idef->iparams[k].bham.c;
      }
    }
  }
  else {
    snew(nbfp,2*atnr*atnr);
    for(i=k=0; (i<atnr); i++) {
      for(j=0; (j<atnr); j++,k++) {
	C6(nbfp,atnr,i,j)   = idef->iparams[k].lj.c6;
	C12(nbfp,atnr,i,j)  = idef->iparams[k].lj.c12;
      }
    }
  }
  return nbfp;
}

/* This routine sets fr->solvent_opt to the most common solvent in the 
 * system, e.g. esolSPC or esolTIP4P. It will also mark each charge group in 
 * the fr->solvent_type array with the correct type (or esolNO).
 *
 * Charge groups that fulfill the conditions but are not identical to the
 * most common one will be marked as esolNO in the solvent_type array. 
 *
 * TIP3p is identical to SPC for these purposes, so we call it
 * SPC in the arrays (Apologies to Bill Jorgensen ;-)
 */
static void
check_solvent(FILE *                fp,
              const t_topology *    top,
              t_forcerec *          fr,
			  const t_mdatoms *     md,
              const t_nsborder *    nsb)
{
    const t_block *   cgs;
    const t_block *   excl;
    const t_block *   mols;
    atom_id *         cgid;
    int               i,j,k;
    int               j0,j1,nj;
    int               aj,ak;
    int               n_solvent_parameters;
    bool              is_single_chargegroup;
    bool              excluded[4];
    bool              all_excluded;
	bool              perturbed;
    bool              has_vdw[4];
    bool              match;
    real              tmp_charge[4];
    int               tmp_vdwtype[4];
    int               nexcl;
    int               tjA;
    int               bestsol;
    
    /* We use a list with parameters for each solvent type. 
     * Every time we discover a new molecule that fulfills the basic 
     * conditions for a solvent we compare with the previous entries
     * in these lists. If the parameters are the same we just increment
     * the counter for that type, and otherwise we create a new type
     * based on the current molecule.
     *
     * Once we've finished going through all molecules we check which
     * solvent is most common, and mark all those molecules while we
     * clear the flag on all others.
     */   
    struct 
    {
        int    model;          
        int    count;
        int    vdwtype[4];
        real   charge[4];
    } * solvent_parameters;
    

    cgs  = &(top->blocks[ebCGS]);
    excl = &(top->atoms.excl);
    mols = &(top->blocks[ebMOLS]);
    
    if (debug)
        fprintf(debug,"Going to determine what solvent types we have.\n");
  
    snew(fr->solvent_type,cgs->nr);
   
    /* Overkill to allocate a separate set of parameters for each molecule, 
	 * but it is not a lot of memory, will be released soon, and it is better
     * than calling realloc lots of times.
     */   	 
    snew(solvent_parameters,mols->nr);
    
    
    n_solvent_parameters = 0;
        
    /* Generate charge group number for all atoms */
    cgid = make_invblock(cgs,cgs->nra);

    /* Loop over molecules */
    if (debug)
        fprintf(debug,"There are %d molecules, %d charge groups and %d atoms\n",
                mols->nr,cgs->nr,cgs->nra);


    /* Set a temporary "no solvent" flag everywhere */
    for(i=0; i<cgs->nr; i++) 
    {
        fr->solvent_type[i] = -1;
    }
    
    /* Go through all molecules in system and see if they fulfill the conditions
     * for the water innerloops. We don't care if they really ARE waters or not,
     * but we need to check charges and LJ to make sure all molecules of a 
     * certain solvent type have the same parameters.
     */   
    
    for(i=0; i<mols->nr; i++) {
        
        /* Get the indices of the first atom in this and the next molecule */
        j0     = mols->index[i];
        j1     = mols->index[i+1];

        /* Number of atoms in our molecule */
        nj     = j1-j0;

        /* Check if it could be an SPC (3 atoms) or TIP4p (4) water, otherwise skip it */
        if(nj<3 || nj>4)
        {
            continue;
        }
        
        /* Check that all subsequent atoms in the molecule belong to
         * the same chargegroup as the first.
         */   

        is_single_chargegroup = TRUE;
        
        for( j = j0+1 ; j<j1 && is_single_chargegroup ; j++ ) 
        {
            if(cgid[mols->a[j]] != cgid[mols->a[j-1]])
            {
                is_single_chargegroup = FALSE;
            }
        }
        
        if(is_single_chargegroup == FALSE)
        {
            /* Cannot be a solvent, go to next molecule */
            continue;
        }

        /* It has 3 or 4 atoms in a single chargegroup. Next we check if
         * all intra-molecular interactions are excluded. 
         */
		
        all_excluded = TRUE;

        /* Loop over all atoms in molecule */
        for( j = j0 ; j<j1 && all_excluded ; j++ ) 
        {

            for( k=0 ; k<nj ; k++ )
            {
                excluded[k] = FALSE;
            }

            excluded[j-j0] = TRUE;
            
             /* Count the exclusions for this atom. */
            for(k=excl->index[j]; k<excl->index[j+1]; k++) 
            {
                /* Index of the atom we exclude from j */
                ak = excl->a[k];
               
                /* There cannot be any exclusions to other molecules! */
                if ((ak < j0) || (ak >= j1)) 
                    gmx_fatal(FARGS,"Exclusion outside molecule? ak = %d, j0 = %d, j1 = 5d, mol is %d",ak,j0,j1,i);

                /* Indices are relative to the first atom in the molecule */
                excluded[ak-j0] = TRUE;
            }
            
            /* Make really sure everything is excluded, even if we did something
             * stupid like list one exclusion multiple times in the list above.
             */   
            for(k=0;k<nj && all_excluded;k++)
            {
                all_excluded = all_excluded & excluded[k];
            }
        }
        
        /* Skip to next molecule if any intramolecules weren't excluded */
        if(all_excluded == FALSE)
        {
            continue;
        }
        
        /* Still looks like a solvent, time to check parameters */

        /* If it is perturbed (free energy) we can't use the solvent loops,
         * so then we just skip to the next molecule.
         */   
		perturbed = FALSE; 
		
        for(j=j0 ; j<j1 && !perturbed; j++)
        {
			perturbed =	md->bPerturbed[mols->a[j]];
        }
	
	    if(perturbed)
		{
			continue;
		}

        /* Now it's only a question if the VdW and charge parameters 
         * are OK. Before doing the check we compare and see if they are 
         * identical to a possible previous solvent type.
         * First we assign the current types and charges.    
         */
		for( j=0; j<nj ; j++)
		{
			aj = mols->a[j0+j];
		    tmp_vdwtype[j] = top->atoms.atom[aj].type;
			tmp_charge[j]  = top->atoms.atom[aj].q;
		} 
		 
        /* Does it match any previous solvent type? */
        match = FALSE;

        for(k=0 ; k<n_solvent_parameters && !match; k++)
        {
            match = TRUE;
         
            
            /* We can only match SPC with 3 atoms and TIP4p with 4 atoms */
            if( (solvent_parameters[k].model==esolSPC && nj!=3)  ||
                (solvent_parameters[k].model==esolTIP4P && nj!=4) )
			    match = FALSE;

            /* Check that types & charges match for all atoms in molecule */
            for(j=0 ; j<nj && match==TRUE; j++)
            {			
                if(tmp_vdwtype[j] != solvent_parameters[k].vdwtype[j])
                    match = FALSE;
                
                if(tmp_charge[j] != solvent_parameters[k].charge[j])
                    match = FALSE;
            }
            if(match == TRUE)
            {
                /* Congratulations! We have a matched solvent.
                * Flag it with this type for later processing.
                */
				aj = mols->a[j0];
                fr->solvent_type[cgid[aj]] = k;
                (solvent_parameters[k].count)++;
            }
        }
        
        if(match == TRUE)
        {
            /* Continue to next molecule if we already identified this */
            continue;
        }
        
        /* If we get here, we have a tentative new solvent type.
         * Before we add it we must check that it fulfills the requirements
         * of the solvent optimized loops. First determine which atoms have
         * VdW interactions.   
         */
        for(j=0; j<nj; j++) 
        {
            has_vdw[j] = FALSE;
            tjA        = tmp_vdwtype[j];
            
            /* Go through all other tpes and see if any have non-zero
            * VdW parameters when combined with this one.
            */   
            for(k=0; k<fr->ntype && (has_vdw[j]==FALSE); k++)
            {
                /* We already checked that the atoms weren't perturbed,
                 * so we only need to check state A now.
                 */ 
                if (fr->bBHAM) 
                {
                    has_vdw[j] = (has_vdw[j] || 
                                  (BHAMA(fr->nbfp,fr->ntype,tjA,k) != 0.0) ||
                                  (BHAMB(fr->nbfp,fr->ntype,tjA,k) != 0.0) ||
                                  (BHAMC(fr->nbfp,fr->ntype,tjA,k) != 0.0));
                }
                else
                {
                   /* Standard LJ */
                    has_vdw[j] = (has_vdw[j] || 
                                  (C6(fr->nbfp,fr->ntype,tjA,k)  != 0.0) ||
                                  (C12(fr->nbfp,fr->ntype,tjA,k) != 0.0));
                }
            }
        }
        
        /* Now we know all we need to make the final check and assignment. */
        if(nj==3)
        {
            /* So, is it an SPC?
             * For this we require thatn all atoms have charge, 
             * the charges on atom 2 & 3 should be the same, and only
             * atom 1 should have VdW.
             */
            if(has_vdw[0] == TRUE && 
               has_vdw[1] == FALSE &&
               has_vdw[2] == FALSE &&
               tmp_charge[0]  != 0 &&
               tmp_charge[1]  != 0 &&
               tmp_charge[2]  == tmp_charge[1])
            {
                solvent_parameters[n_solvent_parameters].model=esolSPC;
                solvent_parameters[n_solvent_parameters].count=1;
                for(k=0;k<3;k++)
                {
                    solvent_parameters[n_solvent_parameters].vdwtype[k] = tmp_vdwtype[k];
                    solvent_parameters[n_solvent_parameters].charge[k]  = tmp_charge[k];
                }
            }
			aj = mols->a[j0];
            fr->solvent_type[cgid[aj]] = n_solvent_parameters;
            n_solvent_parameters++;
        }
        else if(nj==4)
        {
            /* Or could it be a TIP4P?
            * For this we require thatn atoms 2,3,4 have charge, but not atom 1. 
            * Only atom 1 should have VdW.
            */
            if(has_vdw[0] == TRUE && 
               has_vdw[1] == FALSE &&
               has_vdw[2] == FALSE &&
               has_vdw[3] == FALSE &&
               tmp_charge[1]  != 0 &&
               tmp_charge[2]  == tmp_charge[1] &&
               tmp_charge[3]  != 0)
            {
                solvent_parameters[n_solvent_parameters].model=esolTIP4P;
                solvent_parameters[n_solvent_parameters].count=1;
                for(k=0;k<4;k++)
                {
                    solvent_parameters[n_solvent_parameters].vdwtype[k] = tmp_vdwtype[k];
                    solvent_parameters[n_solvent_parameters].charge[k]  = tmp_charge[k];
                }
            }
			aj = mols->a[j0];
            fr->solvent_type[cgid[aj]] = n_solvent_parameters;
            n_solvent_parameters++;
        }
    }
    
   
    /* Puh! We finished going through all molecules. Now find the most
     * common solvent model.
     */   
    
    /* Most common solvent this far */
    j = -2;
    for(i=0;i<n_solvent_parameters;i++)
    {
        if(j==-2 || (solvent_parameters[i].count>solvent_parameters[j].count))
            j = i;
    }
    
    if(j>=0)
    {
        bestsol = solvent_parameters[j].model;
    }
    else
    {
        bestsol = esolNO;
    }
    
#ifdef DISABLE_WATER_NLIST
	bestsol = esolNO;
#endif

	fr->nWatMol = 0;
    for(i=0;i<cgs->nr;i++)
    {
        if(fr->solvent_type[i]==j)
        {
            fr->solvent_type[i] = bestsol;
            fr->nWatMol++;
        }
        else
        {
            fr->solvent_type[i] = esolNO;
        }
    }
    
    if(bestsol!=esolNO && fp!=NULL)
    {
        fprintf(fp,"\nEnabling %s water optimization for %d molecules.\n\n",
                esol_names[bestsol],
                solvent_parameters[j].count);
    }
    
    sfree(solvent_parameters);
    fr->solvent_opt = bestsol;
}
 


void set_chargesum(FILE *log,t_forcerec *fr,const t_mdatoms *mdatoms)
{
  double qsum;
  int    i;

  qsum = 0;
  for(i=0; i<mdatoms->nr; i++)
    qsum += mdatoms->chargeA[i];
  fr->qsum[0] = qsum;
  if (fr->efep != efepNO) {
    qsum = 0;
    for(i=0; i<mdatoms->nr; i++)
    qsum += mdatoms->chargeB[i];
    fr->qsum[1] = qsum;
  } else {
    fr->qsum[1] = fr->qsum[0];
  }
  if (log) {
    if (fr->efep == efepNO)
      fprintf(log,"System total charge: %.3f\n",fr->qsum[0]);
    else
      fprintf(log,"System total charge, top. A: %.3f top. B: %.3f\n",
	      fr->qsum[0],fr->qsum[1]);
  }
}

void update_forcerec(FILE *log,t_forcerec *fr,matrix box)
{
  if (fr->epsilon_r != 0)
    fr->epsfac = ONE_4PI_EPS0/fr->epsilon_r;
  else
    /* eps = 0 is infinite dieletric: no coulomb interactions */
    fr->epsfac = 0;
  
  if (EEL_RF(fr->eeltype))
    calc_rffac(log,fr->eeltype,fr->epsilon_r,fr->epsilon_rf,
	       fr->rcoulomb,fr->temp,fr->zsquare,box,
	       &fr->kappa,&fr->k_rf,&fr->c_rf);
}

void set_avcsixtwelve(FILE *log,t_forcerec *fr,
		      const t_mdatoms *mdatoms,const t_block *excl)
{
  int    i,j,tpi,tpj,j1,j2,k,n,nexcl;
  long   npair,npair_ij;
  double csix,ctwelve;
  int    natoms,ntp,*type,*typecount;
  bool   bBHAM;
  real   *nbfp;
  atom_id *AA;

  natoms = mdatoms->nr;
  ntp = fr->ntype;
  type = mdatoms->typeA;
  bBHAM = fr->bBHAM;
  nbfp = fr->nbfp;
  AA = excl->a;

  csix = 0;
  ctwelve = 0;
  npair = 0;
  nexcl = 0;
  if (!fr->bTPI) {
    /* Count the types so we avoid natoms^2 operations */
    snew(typecount,ntp);
    for(i=0; i<natoms; i++)
      typecount[type[i]]++;
    for(tpi=0; tpi<ntp; tpi++) {
      for(tpj=tpi; tpj<ntp; tpj++) {
	if (tpi != tpj)
	  npair_ij = (long)typecount[tpi]*(long)typecount[tpj];
	else
	  npair_ij = (long)typecount[tpi]*((long)typecount[tpi] - 1)/2;
	if (bBHAM) {
	  csix    += npair_ij*BHAMC(nbfp,ntp,tpi,tpj);
	} else {
	  csix    += npair_ij*   C6(nbfp,ntp,tpi,tpj);
	  ctwelve += npair_ij*  C12(nbfp,ntp,tpi,tpj);
	}
	npair += npair_ij;
      }
    }
    sfree(typecount);
    /* Subtract the excluded pairs.
     * The main reason for substracting exclusions is that in some cases some
     * combinations might never occur and the parameters could have any value.
     * These unused values should not influence the dispersion correction.
     */
    for(i=0; (i<natoms); i++) {
      tpi = type[i];
      j1  = excl->index[i];
      j2  = excl->index[i+1];
      for(j=j1; j<j2; j++) {
	k = AA[j];
	if (k > i) {
	  tpj   = type[k];
	  if (bBHAM) {
	    csix -= BHAMC(nbfp,ntp,tpi,tpj);
	  } else {
	    csix    -= C6 (nbfp,ntp,tpi,tpj);
	    ctwelve -= C12(nbfp,ntp,tpi,tpj);
	  }
	  nexcl++;
	}
      }
    }
  } else {
    /* Only correct for the interaction of the test particle
     * with the rest of the system.
     */
    tpi = type[natoms - 1];
    for(j=0; (j<natoms-1); j++) {
      tpj   = type[j];
      if (bBHAM) {
	csix += BHAMC(nbfp,ntp,tpi,tpj);
      } else {
	csix    += C6 (nbfp,ntp,tpi,tpj);
	ctwelve += C12(nbfp,ntp,tpi,tpj);
      }
      npair++;
    }
  }
  csix    /= npair - nexcl;
  ctwelve /= npair - nexcl;
  if (debug) {
    fprintf(debug,"Counted %d exclusions\n",nexcl);
    fprintf(debug,"Average C6 parameter is: %10g\n",csix);
    fprintf(debug,"Average C12 parameter is: %10g\n",ctwelve);
  }
  fr->avcsix = csix;
  fr->avctwelve = ctwelve;
}

static void set_bham_b_max(FILE *log,t_forcerec *fr,const t_mdatoms *mdatoms)
{
  int  i,j,tpi,tpj,ntypes,natoms,*type;
  real b,bmin;
  real *nbfp;

  fprintf(log,"Determining largest Buckingham b parameter for table\n");
  nbfp   = fr->nbfp;
  ntypes = fr->ntype;
  type   = mdatoms->typeA;
  natoms = mdatoms->nr;

  bmin           = -1;
  fr->bham_b_max = 0;
  for(i=0; (i<natoms); i++) {
    tpi = type[i];
    if (tpi >= ntypes)
      gmx_fatal(FARGS,"Atomtype[%d] = %d, maximum = %d",i,tpi,ntypes);
    
    for(j=0; (j<natoms); j++) {
      tpj   = type[j];
      if (tpj >= ntypes)
	gmx_fatal(FARGS,"Atomtype[%d] = %d, maximum = %d",j,tpj,ntypes);
      b = BHAMB(nbfp,ntypes,tpi,tpj);
      if (b > fr->bham_b_max)
	fr->bham_b_max = b;
      if ((b < bmin) || (bmin==-1))
	bmin = b;
    }
  }
  fprintf(log,"Buckingham b parameters, min: %g, max: %g\n",
	  bmin,fr->bham_b_max);
}

static void make_nbf_tables(FILE *fp,t_forcerec *fr,real rtab,
			    const t_commrec *cr,
			    const char *tabfn,char *eg1,char *eg2,
			    t_nblists *nbl)
{
  char buf[STRLEN];
  int i,j;

  sprintf(buf,"%s",tabfn);
  if (eg1 && eg2)
    /* Append the two energy group names */
    sprintf(buf + strlen(tabfn) - strlen(ftp2ext(efXVG)) - 1,"_%s_%s.%s",
	    eg1,eg2,ftp2ext(efXVG));
  nbl->tab = make_tables(fp,fr,MASTER(cr),buf,rtab,FALSE);
  /* Copy the contents of the table to separate coulomb and LJ tables too,
   * to improve cache performance.
   */
  snew(nbl->coultab,4*(nbl->tab.n+1));
  snew(nbl->vdwtab,8*(nbl->tab.n+1));  
  for(i=0; i<=nbl->tab.n; i++) {
    for(j=0; j<4; j++)
      nbl->coultab[4*i+j] = nbl->tab.tab[12*i+j];
    for(j=0; j<8; j++)
      nbl->vdwtab [8*i+j] = nbl->tab.tab[12*i+4+j];
  }
}

void init_forcerec(FILE *fp,
		   t_forcerec *fr,
		   const t_inputrec *ir,
		   const t_topology *top,
		   const t_commrec  *cr,
		   const t_mdatoms  *mdatoms,
		   const t_nsborder *nsb,
		   matrix     box,
		   bool       bMolEpot,
		   const char *tabfn,
		   const char *tabpfn,
		   bool       bNoSolvOpt)
{
  int     i,j,m,natoms,ngrp,negptable,egi,egj;
  real    q,zsq,nrdf,T,rtab;
  rvec    box_size;
  const t_block *mols,*cgs;
  const t_idef  *idef;
  bool    bTab,bSep14tab,bNormalnblists;
  t_nblists *nbl;
  int     *nm_ind,egp_flags;

  if (check_box(box))
    gmx_fatal(FARGS,check_box(box));

  cgs            = &(top->blocks[ebCGS]);
  mols           = &(top->blocks[ebMOLS]);
  idef           = &(top->idef);
  
  natoms         = mdatoms->nr;

  /* Test particle insertion ? */
  fr->bTPI = (ir->eI == eiTPI);

  /* Copy the user determined parameters */
  fr->userint1 = ir->userint1;
  fr->userint2 = ir->userint2;
  fr->userint3 = ir->userint3;
  fr->userint4 = ir->userint4;
  fr->userreal1 = ir->userreal1;
  fr->userreal2 = ir->userreal2;
  fr->userreal3 = ir->userreal3;
  fr->userreal4 = ir->userreal4;

  /* Shell stuff */
  fr->fc_stepsize = ir->fc_stepsize;

  /* Free energy */
  fr->efep       = ir->efep;
  fr->sc_alpha   = ir->sc_alpha;
  fr->sc_sigma6  = pow(ir->sc_sigma,6);

  /* Neighbour searching stuff */
  fr->bGrid      = (ir->ns_type == ensGRID);
  fr->ndelta     = ir->ndelta;
  fr->ePBC       = ir->ePBC;
  fr->rlist      = ir->rlist;
  fr->rlistlong  = max(ir->rlist,max(ir->rcoulomb,ir->rvdw));
  fr->eeltype    = ir->coulombtype;
  fr->vdwtype    = ir->vdwtype;

  fr->bTwinRange = fr->rlistlong > fr->rlist;
  fr->bEwald     = (fr->eeltype==eelPME || fr->eeltype==eelPMEUSER || 
		    fr->eeltype==eelEWALD);
  fr->bvdwtab    = (fr->vdwtype != evdwCUT);
  
  fr->bcoultab   = (fr->eeltype != eelCUT) && !EEL_RF(fr->eeltype);
  
  if (getenv("GMX_FORCE_TABLES")) {
    fr->bvdwtab  = TRUE;
    fr->bcoultab = TRUE;
  }

  if (fp) {
    fprintf(fp,"Table routines are used for coulomb: %s\n",bool_names[fr->bcoultab]);
    fprintf(fp,"Table routines are used for vdw:     %s\n",bool_names[fr->bvdwtab ]);
  }
  
  /* Tables are used for direct ewald sum */
  if(fr->bEwald) {
    fr->ewaldcoeff=calc_ewaldcoeff(ir->rcoulomb, ir->ewald_rtol);
    if (fp)
      fprintf(fp,"Using a Gaussian width (1/beta) of %g nm for Ewald\n",
	      1/fr->ewaldcoeff);
  }

  /* Domain decomposition parallellism... */
  fr->bDomDecomp = ir->bDomDecomp;
  fr->Dimension  = ir->decomp_dir;
  
  /* Electrostatics */
  fr->epsilon_r  = ir->epsilon_r;
  fr->epsilon_rf = ir->epsilon_rf;
  fr->fudgeQQ    = ir->fudgeQQ;
  fr->rcoulomb_switch = ir->rcoulomb_switch;
  fr->rcoulomb        = ir->rcoulomb;

  /* Parameters for generalized RF */
  fr->zsquare = 0.0;
  fr->temp    = 0.0;
  
  if (fr->eeltype == eelGRF) {
    if (ir->efep!=efepNO)
      fprintf(fp,"\nWARNING: the generalized reaction field constants are determined from topology A only\n\n");
    zsq = 0.0;
    for (i=0; (i<cgs->nr); i++) {
      q = 0;
      for(j=cgs->index[i]; (j<cgs->index[i+1]); j++)
	q+=mdatoms->chargeA[cgs->a[j]];
      if (q != 0.0)
	/* Changed from square to fabs 990314 DvdS 
	 * Does not make a difference for monovalent ions, but doe for 
	 * divalent ions (Ca2+!!)
	 */
	zsq += fabs(q);
    }
    fr->zsquare = zsq;
    
    T    = 0.0;
    nrdf = 0.0;
    for(i=0; (i<ir->opts.ngtc); i++) {
      nrdf += ir->opts.nrdf[i];
      T    += (ir->opts.nrdf[i] * ir->opts.ref_t[i]);
    }
    if (nrdf == 0) 
      gmx_fatal(FARGS,"No degrees of freedom!");
    fr->temp   = T/nrdf;
  }
  else if (EEL_FULL(fr->eeltype) || (fr->eeltype == eelSHIFT) || 
	   (fr->eeltype == eelUSER) || (fr->eeltype == eelSWITCH)) {
    /* We must use the long range cut-off for neighboursearching...
     * An extra range of e.g. 0.1 nm (half the size of a charge group)
     * is necessary for neighboursearching. This allows diffusion 
     * into the cut-off range (between neighborlist updates), 
     * and gives more accurate forces because all atoms within the short-range
     * cut-off rc must be taken into account, while the ns criterium takes
     * only those with the center of geometry within the cut-off.
     * (therefore we have to add half the size of a charge group, plus
     * something to account for diffusion if we have nstlist > 1)
     */
    for(m=0; (m<DIM); m++)
      box_size[m]=box[m][m];

    if (fr->phi == NULL)
      snew(fr->phi,mdatoms->nr);
    
    if ((fr->eeltype==eelPPPM) || (fr->eeltype==eelPOISSON) || 
	(fr->eeltype == eelSHIFT && fr->rcoulomb > fr->rcoulomb_switch))
	set_shift_consts(fp,fr->rcoulomb_switch,fr->rcoulomb,box_size,fr);
  }

  /* Initiate arrays */
  if (fr->bTwinRange) {
    snew(fr->f_twin,natoms);
    snew(fr->fshift_twin,SHIFTS);
  }
  
  if (EEL_FULL(fr->eeltype)) {
    if (ir->efep != efepNO) {
      if (fr->eeltype == eelEWALD || fr->eeltype == eelPME || 
	  fr->eeltype == eelPMEUSER) {
	if (fr->sc_alpha != 0)
	  fprintf(fp,
		  "\nWARNING: With %s the soft-core is performed on erfc(beta r)/r instead of on 1/r\n\n",eel_names[fr->eeltype]);
      } else {
	fprintf(fp,"\nWARNING: With %s the reciprocal part only uses the charges from topology A\n\n",eel_names[fr->eeltype]);
      }
    }
    snew(fr->f_el_recip,natoms);
  }
  
  /* Mask that says whether or not this NBF list should be computed */
  /*  if (fr->bMask == NULL) {
      ngrp = ir->opts.ngener*ir->opts.ngener;
      snew(fr->bMask,ngrp);*/
  /* Defaults to always */
  /*    for(i=0; (i<ngrp); i++)
	fr->bMask[i] = TRUE;
	}*/
  
  if (fr->cg_cm == NULL)
    snew(fr->cg_cm,cgs->nr);
  if (fr->shift_vec == NULL)
    snew(fr->shift_vec,SHIFTS);
    
  if (fr->fshift == NULL)
    snew(fr->fshift,SHIFTS);
  
  if (bMolEpot && (fr->nmol==0)) {
    fr->nmol=mols->nr;
    fr->mol_nr=make_invblock(mols,natoms);
    snew(fr->mol_epot,fr->nmol);
    fr->nstcalc=ir->nstenergy;
  }
  
  if (fr->nbfp == NULL) {
    fr->ntype = idef->atnr;
    fr->bBHAM = (idef->functype[0] == F_BHAM);
    fr->nbfp  = mk_nbfp(idef,fr->bBHAM);
  }
  
  /* Copy the energy group exclusions */
  fr->egp_flags = ir->opts.egp_flags;

  /* Van der Waals stuff */
  fr->rvdw        = ir->rvdw;
  fr->rvdw_switch = ir->rvdw_switch;
  if ((fr->vdwtype != evdwCUT) && (fr->vdwtype != evdwUSER) && !fr->bBHAM) {
    if (fr->rvdw_switch >= fr->rvdw)
      gmx_fatal(FARGS,"rvdw_switch (%g) must be < rvdw (%g)",
		  fr->rvdw_switch,fr->rvdw);
    if (fp)
      fprintf(fp,"Using %s Lennard-Jones, switch between %g and %g nm\n",
	      (fr->eeltype==eelSWITCH) ? "switched":"shifted",
	      fr->rvdw_switch,fr->rvdw);
  } 

  if (fr->bBHAM && (fr->vdwtype == evdwSHIFT || fr->vdwtype == evdwSWITCH))
    gmx_fatal(FARGS,"Switch/shift interaction not supported with Buckingham");
  
  if (fp)
    fprintf(fp,"Cut-off's:   NS: %g   Coulomb: %g   %s: %g\n",
	    fr->rlist,fr->rcoulomb,fr->bBHAM ? "BHAM":"LJ",fr->rvdw);
  
  if (ir->eDispCorr != edispcNO)
    set_avcsixtwelve(fp,fr,mdatoms,&top->atoms.excl);

  if (fr->bBHAM)
    set_bham_b_max(fp,fr,mdatoms);

  /* Copy the GBSA data (radius, volume and surftens for each
   * atomtype) from the topology atomtype section to forcerec.
   */
  snew(fr->atype_radius,fr->ntype);
  snew(fr->atype_vol,fr->ntype);
  snew(fr->atype_surftens,fr->ntype);
  if (top->atomtypes.nr > 0) {
    for(i=0;i<fr->ntype;i++)
      fr->atype_radius[i]=top->atomtypes.radius[i];
    for(i=0;i<fr->ntype;i++)
      fr->atype_vol[i]=top->atomtypes.vol[i];
    for(i=0;i<fr->ntype;i++)
      fr->atype_surftens[i]=top->atomtypes.surftens[i];
  }    

  /* Now update the rest of the vars */
  update_forcerec(fp,fr,box);
  
  set_chargesum(fp,fr,mdatoms);

  /* if we are using LR electrostatics, and they are tabulated,
   * the tables will contain modified coulomb interactions.
   * Since we want to use the non-shifted ones for 1-4
   * coulombic interactions, we must have an extra set of tables.
   */

  /* Construct tables.
   * A little unnecessary to make both vdw and coul tables sometimes,
   * but what the heck... */

  bTab = fr->bcoultab || fr->bvdwtab || (fr->efep != efepNO);
  bSep14tab = ((top->idef.il[F_LJ14].multinr[cr->nnodes-1] > 0) &&
	       (!bTab || fr->eeltype!=eelCUT || fr->vdwtype!=evdwCUT));
  
  negptable = 0;
  if (!bTab) {
    bNormalnblists = TRUE;
    fr->nnblists = 1;
  } else {
    bNormalnblists = (ir->eDispCorr != edispcNO);
    for(egi=0; egi<ir->opts.ngener; egi++) {
      for(egj=egi;  egj<ir->opts.ngener; egj++) {
	egp_flags = ir->opts.egp_flags[GID(egi,egj,ir->opts.ngener)];
	if (!(egp_flags & EGP_EXCL)) {
	  if (egp_flags & EGP_TABLE) {
	    negptable++;
	  } else {
	    bNormalnblists = TRUE;
	  }
	}
      }
    }
    if (bNormalnblists) {
      fr->nnblists = negptable + 1;
    } else {
      fr->nnblists = negptable;
    }
    if (fr->nnblists > 1)
      snew(fr->gid2nblists,ir->opts.ngener*ir->opts.ngener);
  }
  snew(fr->nblists,fr->nnblists);

  if (bTab) {
    rtab = fr->rlistlong + ir->tabext;
    if (fr->rlistlong == 0 && fr->efep != efepNO && rtab < 5) {
      rtab = 5;
      if (fp)
	fprintf(fp,
		"\nWARNING: Increasing the free energy table length from\n"
		"         table extension = %f nm to %g nm,\n"
		"         you can set the table extension in the mdp file\n\n",
		ir->tabext,rtab);
    }
    /* make tables for ordinary interactions */
    
    if (bNormalnblists) {
      make_nbf_tables(fp,fr,rtab,cr,tabfn,NULL,NULL,&fr->nblists[0]);
      if (!bSep14tab)
	fr->tab14 = fr->nblists[0].tab;
      m = 1;
    } else {
      m = 0;
    }
    if (negptable > 0) {
      /* Read the special tables for certain energy group pairs */
      nm_ind = top->atoms.grps[egcENER].nm_ind;
      for(egi=0; egi<ir->opts.ngener; egi++) {
	for(egj=egi;  egj<ir->opts.ngener; egj++) {
	  egp_flags = ir->opts.egp_flags[GID(egi,egj,ir->opts.ngener)];
	  if ((egp_flags & EGP_TABLE) && !(egp_flags & EGP_EXCL)) {
	    nbl = &(fr->nblists[m]);
	    fr->gid2nblists[GID(egi,egj,ir->opts.ngener)] = m;
	    /* Read the table file with the two energy groups names appended */
	    make_nbf_tables(fp,fr,rtab,cr,tabfn,
			    *top->atoms.grpname[nm_ind[egi]],
			    *top->atoms.grpname[nm_ind[egj]],
			    &fr->nblists[m]);
	    m++;
	  } else {
	    fr->gid2nblists[GID(egi,egj,ir->opts.ngener)] = 0;
	  }
	}
      }
    }
  }
  if (bSep14tab)
    /* generate extra tables with plain Coulomb for 1-4 interactions only */
    fr->tab14 = make_tables(fp,fr,MASTER(cr),tabpfn,ir->tabext,TRUE);

  check_solvent(fp,top,fr,mdatoms,nsb);

  
  if (getenv("GMX_NO_SOLV_OPT")) {
    if (fp)
      fprintf(fp,"Found environment variable GMX_NO_SOLV_OPT.\n"
  	      "Disabling all solvent optimization\n");
    fr->solvent_opt = esolNO;
  }
  if (bNoSolvOpt)
    fr->solvent_opt = esolNO;
  if (!fr->solvent_opt) {
    for(i=0; i<top->blocks[ebCGS].nr; i++) 
      fr->solvent_type[i] = esolNO;
  }
}
 
#define pr_real(fp,r) fprintf(fp,"%s: %e\n",#r,r)
#define pr_int(fp,i)  fprintf((fp),"%s: %d\n",#i,i)
#define pr_bool(fp,b) fprintf((fp),"%s: %s\n",#b,bool_names[b])

void pr_forcerec(FILE *fp,t_forcerec *fr,t_commrec *cr)
{
  int i;

  pr_real(fp,fr->rlist);
  pr_real(fp,fr->rcoulomb);
  pr_real(fp,fr->fudgeQQ);
  pr_int(fp,fr->ndelta);
  pr_bool(fp,fr->bGrid);
  pr_bool(fp,fr->bTwinRange);
  /*pr_int(fp,fr->cg0);
    pr_int(fp,fr->hcg);*/
  for(i=0; i<fr->nnblists; i++)
    pr_int(fp,fr->nblists[i].tab.n);
  pr_real(fp,fr->rcoulomb_switch);
  pr_real(fp,fr->rcoulomb);
  
  pr_int(fp,fr->nmol);
  pr_int(fp,fr->nstcalc);
  
  fflush(fp);
}

void ns(FILE *fp,
	t_forcerec *fr,
	rvec       x[],
	rvec       f[],
	matrix     box,
	t_groups   *grps,
	t_grpopts  *opts,
	t_topology *top,
	t_mdatoms  *md,
	t_commrec  *cr,
	t_nrnb     *nrnb,
	t_nsborder *nsb,
	int        step,
	real       lambda,
	real       *dvdlambda,
	bool       bFillGrid,
	bool       bDoForces)
{
  static bool bFirst=TRUE;
  static int  nDNL;
  char   *ptr;
  int    nsearch;
  
  if (bFirst) {
    ptr=getenv("DUMPNL");
    if (ptr) {
      nDNL=atoi(ptr);
      fprintf(fp,"nDNL = %d\n",nDNL);  
    } else
      nDNL=0;
    /* Allocate memory for the neighbor lists */
    init_neighbor_list(fp,fr,HOMENR(nsb));
      
    bFirst=FALSE;
  }
    
  if (fr->bTwinRange) 
    fr->nlr=0;

  /* Whether or not we do dynamic load balancing,
   * workload contains the proper numbers of charge groups
   * to be searched.
   */
  if (cr->nodeid == 0)
    fr->cg0=0;
  else
    fr->cg0=nsb->workload[cr->nodeid-1];
  fr->hcg=nsb->workload[cr->nodeid];

  nsearch = search_neighbours(fp,fr,x,box,top,grps,cr,nsb,nrnb,md,
			      lambda,dvdlambda,bFillGrid,bDoForces);
  if (debug)
    fprintf(debug,"nsearch = %d\n",nsearch);
    
  /* Check whether we have to do dynamic load balancing */
  /*if ((nsb->nstDlb > 0) && (mod(step,nsb->nstDlb) == 0))
    count_nb(cr,nsb,&(top->blocks[ebCGS]),nns,fr->nlr,
    &(top->idef),opts->ngener);
  */
  if (nDNL > 0)
    dump_nblist(fp,fr,nDNL);
}

void force(FILE       *fplog,   int        step,
	   t_forcerec *fr,      t_inputrec *ir,
	   t_idef     *idef,    t_nsborder *nsb,
	   t_commrec  *cr,      t_commrec *mcr,
	   t_nrnb     *nrnb,
	   t_groups   *grps,    t_mdatoms  *md,
	   int        ngener,   t_grpopts  *opts,
	   rvec       x[],      rvec       f[],
	   real       epot[],   t_fcdata   *fcd,
	   bool       bVerbose, matrix     box,
	   real       lambda,   t_graph    *graph,
	   t_block    *excl,    
	   bool       bNBFonly, bool bDoForces,
	   rvec       mu_tot[],
	   bool       bGatherOnly)
{
  int     i,nit;
  bool    bDoEpot,bSepDVDL;
  rvec    box_size;
  real    dvdlambda,Vsr,Vlr,Vcorr=0,vdip,vcharge;
  t_pbc   pbc;
#define PRINT_SEPDVDL(s,v,dvdl) if (bSepDVDL) fprintf(fplog,"  %-30s V %12.5e  dVdl %12.5e\n",s,v,dvdl);

  /* Reset box */
  for(i=0; (i<DIM); i++)
    box_size[i]=box[i][i];

  bDoEpot=((fr->nmol > 0) && (fr->nstcalc > 0) && (mod(step,fr->nstcalc)==0));
  bSepDVDL=(fr->bSepDVDL && do_per_step(step,ir->nstlog));
  /* Reset epot... */
  if (bDoEpot) 
    for(i=0; (i<fr->nmol); i++)
      fr->mol_epot[i]=0.0;
  debug_gmx();

  if (bSepDVDL)
    fprintf(fplog,"Step %d: non-bonded V and dVdl for node %d:\n",
	    step,cr->nodeid);
  
  /* Call the short range functions all in one go. */
  dvdlambda = 0;
  where();
  do_nonbonded(fplog,cr,fr,x,f,md,
	  fr->bBHAM ? grps->estat.ee[egBHAMSR] : grps->estat.ee[egLJSR],
	  grps->estat.ee[egCOULSR],box_size,nrnb,
	  lambda,&dvdlambda,FALSE,-1,-1,bDoForces);
  where();

  epot[F_DVDL] += dvdlambda;
  Vsr = 0;
  if (bSepDVDL)
    for(i=0; i<grps->estat.nn; i++)
      Vsr +=
	(fr->bBHAM ? grps->estat.ee[egBHAMSR][i] : grps->estat.ee[egLJSR][i])
	+ grps->estat.ee[egCOULSR][i];
  PRINT_SEPDVDL("VdW and Coulomb SR particle-p.",Vsr,dvdlambda);
  debug_gmx();

  if (debug) 
    pr_rvecs(debug,0,"fshift after SR",fr->fshift,SHIFTS);
  
  /* Shift the coordinates. Must be done before bonded forces and PPPM, 
   * but is also necessary for SHAKE and update, therefore it can NOT 
   * go when no bonded forces have to be evaluated.
   */
  
  /* Check whether we need to do bondeds or correct for exclusions */
  if (!bNBFonly || EEL_RF(fr->eeltype) || EEL_FULL(fr->eeltype)) {
    if (graph) {
      shift_self(graph,box,x);
      if (TRICLINIC(box))
	inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
      else
	inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);
    }
    if (fr->ePBC==epbcFULL || idef->il[F_POSRES].nr>0)
      set_pbc_ss(&pbc,box);
    debug_gmx();
  }
  where();
  if (EEL_FULL(fr->eeltype)) {
    dvdlambda = 0;
    switch (fr->eeltype) {
    case eelPPPM:
      Vlr = do_pppm(fplog,FALSE,x,fr->f_el_recip,md->chargeA,
		    box_size,fr->phi,cr,nsb,nrnb,ir->pme_order);
      break;
    case eelPME:
    case eelPMEUSER:
      Vlr = do_pme(fplog,FALSE,ir,x,fr->f_el_recip,md->chargeA,md->chargeB,
		   box,cr,nsb,nrnb,fr->vir_el_recip,fr->ewaldcoeff,
		   md->bChargePerturbed,lambda,&dvdlambda,bGatherOnly);
      PRINT_SEPDVDL("PME mesh",Vlr,dvdlambda);
      break;
    case eelEWALD:
      Vlr = do_ewald(fplog,FALSE,ir,x,fr->f_el_recip,md->chargeA,md->chargeB,
		     box_size,cr,nsb,fr->vir_el_recip,fr->ewaldcoeff,
		     lambda,&dvdlambda);
      PRINT_SEPDVDL("Ewald long-range",Vlr,dvdlambda);
      break;
    default:
      Vlr = 0;
      gmx_fatal(FARGS,"No such electrostatics method implemented %s",
		eel_names[fr->eeltype]);
    }
    epot[F_DVDL] += dvdlambda;
    if(fr->bEwald) {
      dvdlambda = 0;
      Vcorr = ewald_LRcorrection(fplog,nsb,cr,fr,md->chargeA,md->chargeB,excl,
				 x,box,mu_tot,
				 ir->ewald_geometry,ir->epsilon_surface,
				 lambda,&dvdlambda,&vdip,&vcharge);
      PRINT_SEPDVDL("Ewald excl./charge/dip. corr.",Vcorr,dvdlambda);
      epot[F_DVDL] += dvdlambda;
    } else {
      Vcorr = shift_LRcorrection(fplog,nsb,cr,fr,md->chargeA,excl,x,TRUE,box,
				 fr->vir_el_recip);
    }
    epot[F_COUL_RECIP] = Vlr + Vcorr;
    if (debug)
      fprintf(debug,"Vlr = %g, Vcorr = %g, Vlr_corr = %g\n",
	      Vlr,Vcorr,epot[F_COUL_RECIP]);
    if (debug) {
      pr_rvecs(debug,0,"vir_el_recip after corr",fr->vir_el_recip,DIM);
      pr_rvecs(debug,0,"fshift after LR Corrections",fr->fshift,SHIFTS);
    }
  } else if (EEL_RF(fr->eeltype)) {
    dvdlambda = 0;

      if (fr->eeltype != eelRF_NEC)
      epot[F_RF_EXCL] = RF_excl_correction(fplog,nsb,fr,graph,md,excl,x,f,
					   fr->fshift,&pbc,lambda,&dvdlambda);

    epot[F_DVDL] += dvdlambda;
    PRINT_SEPDVDL("RF exclusion correction",epot[F_RF_EXCL],dvdlambda);
  }
  where();
  debug_gmx();
  
  if (debug)    
    print_nrnb(debug,nrnb); 
  debug_gmx();
  
  if (!bNBFonly) {
    calc_bonds(fplog,cr,mcr,
	       idef,x,f,fr,&pbc,graph,epot,nrnb,lambda,md,
	       opts->ngener,grps->estat.ee[egLJ14],grps->estat.ee[egCOUL14],
	       fcd,step,fr->bSepDVDL && do_per_step(step,ir->nstlog));    
    debug_gmx();
  }

  if (debug) 
    pr_rvecs(debug,0,"fshift after bondeds",fr->fshift,SHIFTS);
}
