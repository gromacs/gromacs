/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#include "physics.h"
#include "force.h"
#include "nonbonded.h"
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
#include "coulomb.h"
#include "pppm.h"
#include "pme.h"
#include "mdrun.h"
#include "domdec.h"
#include "partdec.h"
#include "qmmm.h"
#include "mpelogging.h"
#include "copyrite.h"
#include "mtop_util.h"

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

static real *mk_nbfp(const gmx_ffparams_t *idef,bool bBHAM)
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
 * 
 * NOTE: QM particle should not
 * become an optimized solvent. Not even if there is only one charge
 * group in the Qm 
 */

typedef struct 
{
    int    model;          
    int    count;
    int    vdwtype[4];
    real   charge[4];
} solvent_parameters_t;

static void
check_solvent_cg(const gmx_moltype_t   *molt,
                 int                   cg0,
                 const unsigned char   *qm_grpnr,
                 const t_grps          *qm_grps,
                 t_forcerec *          fr,
                 int                   *n_solvent_parameters,
                 solvent_parameters_t  **solvent_parameters_p,
                 int                   cginfo,
                 int                   *cg_sp)
{
    const t_blocka *  excl;
    t_atom            *atom;
    int               j,k;
    int               j0,j1,nj;
    bool              perturbed;
    bool              has_vdw[4];
    bool              match;
    real              tmp_charge[4];
    int               tmp_vdwtype[4];
    int               tjA;
    bool              qm;
    solvent_parameters_t *solvent_parameters;

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

    solvent_parameters = *solvent_parameters_p;

    /* Mark the cg first as non optimized */
    *cg_sp = -1;
    
    /* Check if this cg has no exclusions with atoms in other charge groups
     * and all atoms inside the charge group excluded.
     * We only have 3 or 4 atom solvent loops.
     */
    if (GET_CGINFO_EXCL_INTER(cginfo) ||
        !GET_CGINFO_EXCL_INTRA(cginfo))
    {
        return;
    }

    /* Get the indices of the first atom in this charge group */
    j0     = molt->cgs.index[cg0];
    j1     = molt->cgs.index[cg0+1];
    
    /* Number of atoms in our molecule */
    nj     = j1 - j0;

    if (debug) {
        fprintf(debug,
                "Moltype '%s': there are %d atoms in this charge group\n",
                *molt->name,nj);
    }
    
    /* Check if it could be an SPC (3 atoms) or TIP4p (4) water,
     * otherwise skip it.
     */
    if (nj<3 || nj>4)
    {
        return;
    }
    
    /* Check if we are doing QM on this group */
    qm = FALSE; 
    if (qm_grpnr != NULL)
    {
        for(j=j0 ; j<j1 && !qm; j++)
        {
            qm = (qm_grpnr[j] < qm_grps->nr - 1);
        }
    }
    /* Cannot use solvent optimization with QM */
    if (qm)
    {
        return;
    }
    
    atom = molt->atoms.atom;

    /* Still looks like a solvent, time to check parameters */
    
    /* If it is perturbed (free energy) we can't use the solvent loops,
     * so then we just skip to the next molecule.
     */   
    perturbed = FALSE; 
    
    for(j=j0; j<j1 && !perturbed; j++)
    {
        perturbed = PERTURBED(atom[j]);
    }
    
    if (perturbed)
    {
        return;
    }
    
    /* Now it's only a question if the VdW and charge parameters 
     * are OK. Before doing the check we compare and see if they are 
     * identical to a possible previous solvent type.
     * First we assign the current types and charges.    
     */
    for(j=0; j<nj; j++)
    {
        tmp_vdwtype[j] = atom[j0+j].type;
        tmp_charge[j]  = atom[j0+j].q;
    } 
    
    /* Does it match any previous solvent type? */
    for(k=0 ; k<*n_solvent_parameters; k++)
    {
        match = TRUE;
        
        
        /* We can only match SPC with 3 atoms and TIP4p with 4 atoms */
        if( (solvent_parameters[k].model==esolSPC   && nj!=3)  ||
            (solvent_parameters[k].model==esolTIP4P && nj!=4) )
            match = FALSE;
        
        /* Check that types & charges match for all atoms in molecule */
        for(j=0 ; j<nj && match==TRUE; j++)
        {			
            if (tmp_vdwtype[j] != solvent_parameters[k].vdwtype[j])
            {
                match = FALSE;
            }
            if(tmp_charge[j] != solvent_parameters[k].charge[j])
            {
                match = FALSE;
            }
        }
        if (match == TRUE)
        {
            /* Congratulations! We have a matched solvent.
             * Flag it with this type for later processing.
             */
            *cg_sp = k;
            (solvent_parameters[k].count)++;

            /* We are done with this charge group */
            return;
        }
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
    if (nj == 3)
    {
        /* So, is it an SPC?
         * For this we require thatn all atoms have charge, 
         * the charges on atom 2 & 3 should be the same, and only
         * atom 1 should have VdW.
         */
        if (has_vdw[0] == TRUE && 
            has_vdw[1] == FALSE &&
            has_vdw[2] == FALSE &&
            tmp_charge[0]  != 0 &&
            tmp_charge[1]  != 0 &&
            tmp_charge[2]  == tmp_charge[1])
        {
            srenew(solvent_parameters,*n_solvent_parameters+1);
            solvent_parameters[*n_solvent_parameters].model = esolSPC;
            solvent_parameters[*n_solvent_parameters].count = 1;
            for(k=0;k<3;k++)
            {
                solvent_parameters[*n_solvent_parameters].vdwtype[k] = tmp_vdwtype[k];
                solvent_parameters[*n_solvent_parameters].charge[k]  = tmp_charge[k];
            }

            *cg_sp = *n_solvent_parameters;
            (*n_solvent_parameters)++;
        }
    }
    else if (nj==4)
    {
        /* Or could it be a TIP4P?
         * For this we require thatn atoms 2,3,4 have charge, but not atom 1. 
         * Only atom 1 should have VdW.
         */
        if(has_vdw[0] == TRUE && 
           has_vdw[1] == FALSE &&
           has_vdw[2] == FALSE &&
           has_vdw[3] == FALSE &&
           tmp_charge[0]  == 0 &&
           tmp_charge[1]  != 0 &&
           tmp_charge[2]  == tmp_charge[1] &&
           tmp_charge[3]  != 0)
        {
            srenew(solvent_parameters,*n_solvent_parameters+1);
            solvent_parameters[*n_solvent_parameters].model=esolTIP4P;
            solvent_parameters[*n_solvent_parameters].count=1;
            for(k=0;k<4;k++)
            {
                solvent_parameters[*n_solvent_parameters].vdwtype[k] = tmp_vdwtype[k];
                solvent_parameters[*n_solvent_parameters].charge[k]  = tmp_charge[k];
            }
            
            *cg_sp = *n_solvent_parameters;
            (*n_solvent_parameters)++;
        }
    }

    *solvent_parameters_p = solvent_parameters;
}

static void
check_solvent(FILE *                fp,
              const gmx_mtop_t *    mtop,
              t_forcerec *          fr,
              int                   *cginfo)
{
    const t_block *   cgs;
    const t_block *   mols;
    const gmx_moltype_t *molt;
    int               mb,mol,cg_mol,cg_offset,at_offset,ncg,i;
    int               n_solvent_parameters;
    solvent_parameters_t *solvent_parameters;
    int               *cg_sp;
    int               bestsp,bestsol;

    if (debug)
    {
        fprintf(debug,"Going to determine what solvent types we have.\n");
    }

    mols = &mtop->mols;

    n_solvent_parameters = 0;
    solvent_parameters = NULL;
    /* Allocate temporary array for solvent type */
    ncg = ncg_mtop(mtop);
    snew(cg_sp,ncg);

    cg_offset = 0;
    at_offset = 0;
    for(mb=0; mb<mtop->nmolblock; mb++) {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        cgs  = &molt->cgs;
        /* Here we have to loop over all individual molecules
         * because we need to check for QMMM particles.
         */
        for(mol=0; mol<mtop->molblock[mb].nmol; mol++)
        {
            for(cg_mol=0; cg_mol<cgs->nr; cg_mol++)
            {
                check_solvent_cg(molt,cg_mol,
                                 mtop->groups.grpnr[egcQMMM] ?
                                 mtop->groups.grpnr[egcQMMM]+at_offset : 0,
                                 &mtop->groups.grps[egcQMMM],
                                 fr,
                                 &n_solvent_parameters,&solvent_parameters,
                                 cginfo[cg_offset+cg_mol],
                                 &cg_sp[cg_offset+cg_mol]);
            }
            cg_offset += cgs->nr;
            at_offset += cgs->index[cgs->nr];
        }
    }

    /* Puh! We finished going through all charge groups.
     * Now find the most common solvent model.
     */   
    
    /* Most common solvent this far */
    bestsp = -2;
    for(i=0;i<n_solvent_parameters;i++)
    {
        if (bestsp == -2 ||
            solvent_parameters[i].count > solvent_parameters[bestsp].count)
        {
            bestsp = i;
        }
    }
    
    if (bestsp >= 0)
    {
        bestsol = solvent_parameters[bestsp].model;
    }
    else
    {
        bestsol = esolNO;
    }
    
#ifdef DISABLE_WATER_NLIST
	bestsol = esolNO;
#endif

    fr->nWatMol = 0;
    for(i=0; i<ncg; i++)
    {
        if (cg_sp[i] == bestsp)
        {
            SET_CGINFO_SOLOPT(cginfo[i],bestsol);
            fr->nWatMol++;
        }
        else
        {
            SET_CGINFO_SOLOPT(cginfo[i],esolNO);
        }
    }
    
    if (bestsol != esolNO && fp!=NULL)
    {
        fprintf(fp,"\nEnabling %s water optimization for %d molecules.\n\n",
                esol_names[bestsol],
                solvent_parameters[bestsp].count);
    }

    sfree(cg_sp);
    sfree(solvent_parameters);
    fr->solvent_opt = bestsol;
}

static int *init_cginfo(FILE *fplog,const gmx_mtop_t *mtop,
                        t_forcerec *fr,bool bNoSolvOpt)
{
    const t_block *cgs;
    const t_blocka *excl;
    const gmx_moltype_t *molt;
    const gmx_molblock_t *molb;
    int  *cginfo;
    int  mb,m,ncg_tot,cg_offset,cg,a_offset,a0,a1,gid,ai,j,aj,excl_nalloc;
    bool *bExcl,bExclIntraAll,bExclInter;

    ncg_tot = ncg_mtop(mtop);
    snew(cginfo,ncg_tot);
    
    excl_nalloc = 10;
    snew(bExcl,excl_nalloc);
    cg_offset = 0;
    a_offset  = 0;
    for(mb=0; mb<mtop->nmolblock; mb++)
    {
        molb = &mtop->molblock[mb];
        molt = &mtop->moltype[molb->type];
        cgs  = &molt->cgs;
        excl = &molt->excls;
        for(m=0; m<molb->nmol; m++)
        {
            for(cg=0; cg<cgs->nr; cg++)
            {
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];

                /* Store the energy group in cginfo */
                gid = ggrpnr(&mtop->groups,egcENER,a_offset+a0);
                SET_CGINFO_GID(cginfo[cg_offset+cg],gid);
                
                /* Check the intra/inter charge group exclusions */
                if (a1-a0 > excl_nalloc) {
                    excl_nalloc = a1 - a0;
                    srenew(bExcl,excl_nalloc);
                }
                /* bExclIntraAll: all intra cg interactions excluded
                 * bExclInter:    any inter cg interactions excluded
                 */
                bExclIntraAll = TRUE;
                bExclInter    = FALSE;
                for(ai=a0; ai<a1; ai++) {
                    /* Clear the exclusion list for atom ai */
                    for(aj=a0; aj<a1; aj++) {
                        bExcl[aj-a0] = FALSE;
                    }
                    /* Loop over all the exclusions of atom ai */
                    for(j=excl->index[ai]; j<excl->index[ai+1]; j++)
                    {
                        aj = excl->a[j];
                        if (aj < a0 || aj >= a1)
                        {
                            bExclInter = TRUE;
                        }
                        else
                        {
                            bExcl[aj-a0] = TRUE;
                        }
                    }
                    /* Check if ai excludes a0 to a1 */
                    for(aj=a0; aj<a1; aj++)
                    {
                        if (!bExcl[aj-a0])
                        {
                            bExclIntraAll = FALSE;
                        }
                    }
                }
                if (bExclIntraAll)
                {
                    SET_CGINFO_EXCL_INTRA(cginfo[cg_offset+cg]);
                }
                if (bExclInter)
                {
                    SET_CGINFO_EXCL_INTER(cginfo[cg_offset+cg]);
                }
            }
            cg_offset += cgs->nr;
            a_offset  += cgs->index[cgs->nr];
        }
    }
    sfree(bExcl);
    
    /* the solvent optimizer is called after the QM is initialized,
     * because we don't want to have the QM subsystemto become an
     * optimized solvent
     */

    check_solvent(fplog,mtop,fr,cginfo);
    
    if (getenv("GMX_NO_SOLV_OPT"))
    {
        if (fplog)
        {
            fprintf(fplog,"Found environment variable GMX_NO_SOLV_OPT.\n"
                    "Disabling all solvent optimization\n");
        }
        fr->solvent_opt = esolNO;
    }
    if (bNoSolvOpt)
    {
        fr->solvent_opt = esolNO;
    }
    if (!fr->solvent_opt)
    {
        for(cg=0; cg<ncg_tot; cg++)
        {
            SET_CGINFO_SOLOPT(cginfo[cg],esolNO);
        }
    }
    
    return cginfo;
}


void set_chargesum(FILE *log,t_forcerec *fr,const gmx_mtop_t *mtop)
{
    double qsum;
    int    mb,nmol,i;
    const t_atoms *atoms;
    
    qsum = 0;
    for(mb=0; mb<mtop->nmolblock; mb++)
    {
        nmol  = mtop->molblock[mb].nmol;
        atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
        for(i=0; i<atoms->nr; i++)
        {
            qsum += nmol*atoms->atom[i].q;
        }
    }
    fr->qsum[0] = qsum;
    if (fr->efep != efepNO)
    {
        qsum = 0;
        for(mb=0; mb<mtop->nmolblock; mb++)
        {
            nmol  = mtop->molblock[mb].nmol;
            atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
            for(i=0; i<atoms->nr; i++)
            {
                qsum += nmol*atoms->atom[i].qB;
            }
            fr->qsum[1] = qsum;
        }
    }
    else
    {
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
  GMX_MPE_LOG(ev_update_fr_start);

  if (fr->eeltype == eelGRF)
    calc_rffac(NULL,fr->eeltype,fr->epsilon_r,fr->epsilon_rf,
	       fr->rcoulomb,fr->temp,fr->zsquare,box,
	       &fr->kappa,&fr->k_rf,&fr->c_rf);


  GMX_MPE_LOG(ev_update_fr_finish);
}

void set_avcsixtwelve(FILE *fplog,t_forcerec *fr,const gmx_mtop_t *mtop)
{
    const t_atoms *atoms;
    const t_blocka *excl;
    int    mb,nmol,nmolc,i,j,tpi,tpj,j1,j2,k,n,nexcl,q;
#if (defined SIZEOF_LONG_LONG_INT) && (SIZEOF_LONG_LONG_INT >= 8)    
    long long int  npair,npair_ij,tmpi,tmpj;
#else
    double npair, npair_ij,tmpi,tmpj;
#endif
    double csix,ctwelve;
    int    ntp,*typecount;
    bool   bBHAM;
    real   *nbfp;

    ntp = fr->ntype;
    bBHAM = fr->bBHAM;
    nbfp = fr->nbfp;
    
    for(q=0; q<(fr->efep==efepNO ? 1 : 2); q++) {
        csix = 0;
        ctwelve = 0;
        npair = 0;
        nexcl = 0;
        if (!fr->n_tpi) {
            /* Count the types so we avoid natoms^2 operations */
            snew(typecount,ntp);
            for(mb=0; mb<mtop->nmolblock; mb++) {
                nmol  = mtop->molblock[mb].nmol;
                atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
                for(i=0; i<atoms->nr; i++) {
                    if (q == 0)
                    {
                        tpi = atoms->atom[i].type;
                    }
                    else
                    {
                        tpi = atoms->atom[i].typeB;
                    }
                    typecount[tpi] += nmol;
                }
            }
            for(tpi=0; tpi<ntp; tpi++) {
                for(tpj=tpi; tpj<ntp; tpj++) {
                    tmpi = typecount[tpi];
                    tmpj = typecount[tpj];
                    if (tpi != tpj)
                    {
                        npair_ij = tmpi*tmpj;
                    }
                    else
                    {
                        npair_ij = tmpi*(tmpi - 1)/2;
                    }
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
             * The main reason for substracting exclusions is that in some cases
             * some combinations might never occur and the parameters could have
             * any value. These unused values should not influence the dispersion
             * correction.
             */
            for(mb=0; mb<mtop->nmolblock; mb++) {
                nmol  = mtop->molblock[mb].nmol;
                atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
                excl  = &mtop->moltype[mtop->molblock[mb].type].excls;
                for(i=0; (i<atoms->nr); i++) {
                    if (q == 0)
                    {
                        tpi = atoms->atom[i].type;
                    }
                    else
                    {
                        tpi = atoms->atom[i].typeB;
                    }
                    j1  = excl->index[i];
                    j2  = excl->index[i+1];
                    for(j=j1; j<j2; j++) {
                        k = excl->a[j];
                        if (k > i)
                        {
                            if (q == 0)
                            {
                                tpj = atoms->atom[k].type;
                            }
                            else
                            {
                                tpj = atoms->atom[k].typeB;
                            }
                            if (bBHAM) {
                                csix -= nmol*BHAMC(nbfp,ntp,tpi,tpj);
                            } else {
                                csix    -= nmol*C6 (nbfp,ntp,tpi,tpj);
                                ctwelve -= nmol*C12(nbfp,ntp,tpi,tpj);
                            }
                            nexcl += nmol;
                        }
                    }
                }
            }
        } else {
            /* Only correct for the interaction of the test particle
             * with the rest of the system.
             */
            atoms = &mtop->moltype[mtop->molblock[mtop->nmolblock-1].type].atoms;
            if (q == 0)
            {
                tpi = atoms->atom[atoms->nr-1].type;
            }
            else
            {
                tpi = atoms->atom[atoms->nr-1].typeB;
            }
            npair = 0;
            for(mb=0; mb<mtop->nmolblock; mb++) {
                nmol  = mtop->molblock[mb].nmol;
                atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
                for(j=0; j<atoms->nr; j++) {
                    nmolc = nmol;
                    /* Remove the interaction of the test charge group
                     * with itself.
                     */
                    if (mb == mtop->nmolblock-1 && j >= atoms->nr - fr->n_tpi)
                    {
                        nmolc--;
                    }
                    if (q == 0)
                    {
                        tpj = atoms->atom[j].type;
                    }
                    else
                    {
                        tpj = atoms->atom[j].typeB;
                    }
                    if (bBHAM)
                    {
                        csix    += nmolc*BHAMC(nbfp,ntp,tpi,tpj);
                    }
                    else
                    {
                        csix    += nmolc*C6 (nbfp,ntp,tpi,tpj);
                        ctwelve += nmolc*C12(nbfp,ntp,tpi,tpj);
                    }
                    npair += nmolc;
                }
            }
        }
        if (npair - nexcl <= 0 && fplog) {
            fprintf(fplog,"\nWARNING: There are no atom pairs for dispersion correction\n\n");
            csix     = 0;
            ctwelve  = 0;
        } else {
            csix    /= npair - nexcl;
            ctwelve /= npair - nexcl;
        }
        if (debug) {
            fprintf(debug,"Counted %d exclusions\n",nexcl);
            fprintf(debug,"Average C6 parameter is: %10g\n",(double)csix);
            fprintf(debug,"Average C12 parameter is: %10g\n",(double)ctwelve);
        }
        fr->avcsix[q]    = csix;
        fr->avctwelve[q] = ctwelve;
    }
}


static void set_bham_b_max(FILE *fplog,t_forcerec *fr,
                           const gmx_mtop_t *mtop)
{
    const t_atoms *at1,*at2;
    int  mt1,mt2,i,j,tpi,tpj,ntypes;
    real b,bmin;
    real *nbfp;

    if (fplog)
    {
        fprintf(fplog,"Determining largest Buckingham b parameter for table\n");
    }
    nbfp   = fr->nbfp;
    ntypes = fr->ntype;
    
    bmin           = -1;
    fr->bham_b_max = 0;
    for(mt1=0; mt1<mtop->nmoltype; mt1++)
    {
        at1 = &mtop->moltype[mt1].atoms;
        for(i=0; (i<at1->nr); i++)
        {
            tpi = at1->atom[i].type;
            if (tpi >= ntypes)
                gmx_fatal(FARGS,"Atomtype[%d] = %d, maximum = %d",i,tpi,ntypes);
            
            for(mt2=mt1; mt2<mtop->nmoltype; mt2++)
            {
                at2 = &mtop->moltype[mt2].atoms;
                for(j=0; (j<at2->nr); j++) {
                    tpj = at2->atom[j].type;
                    if (tpj >= ntypes)
                    {
                        gmx_fatal(FARGS,"Atomtype[%d] = %d, maximum = %d",j,tpj,ntypes);
                    }
                    b = BHAMB(nbfp,ntypes,tpi,tpj);
                    if (b > fr->bham_b_max)
                    {
                        fr->bham_b_max = b;
                    }
                    if ((b < bmin) || (bmin==-1))
                    {
                        bmin = b;
                    }
                }
            }
        }
    }
    if (fplog)
    {
        fprintf(fplog,"Buckingham b parameters, min: %g, max: %g\n",
                bmin,fr->bham_b_max);
    }
}

static void make_nbf_tables(FILE *fp,t_forcerec *fr,real rtab,
			    const t_commrec *cr,
			    const char *tabfn,char *eg1,char *eg2,
			    t_nblists *nbl)
{
  char buf[STRLEN];
  int i,j;
  void *      p_tmp;

  if (tabfn == NULL) {
    if (debug)
      fprintf(debug,"No table file name passed, can not read table, can not do non-bonded interactions\n");
    return;
  }
    
  sprintf(buf,"%s",tabfn);
  if (eg1 && eg2)
    /* Append the two energy group names */
    sprintf(buf + strlen(tabfn) - strlen(ftp2ext(efXVG)) - 1,"_%s_%s.%s",
	    eg1,eg2,ftp2ext(efXVG));
  nbl->tab = make_tables(fp,fr,MASTER(cr),buf,rtab,FALSE,FALSE);
  /* Copy the contents of the table to separate coulomb and LJ tables too,
   * to improve cache performance.
   */

  /* For performance reasons we want
   * the table data to be aligned to 16-byte. This is accomplished
   * by allocating 16 bytes extra to a temporary pointer, and then
   * calculating an aligned pointer. This new pointer must not be
   * used in a free() call, but thankfully we're sloppy enough not
   * to do this...
   */
  
  /* 8 fp entries per vdw table point, n+1 points, and 16 bytes extra to align it. */
  p_tmp = malloc(8*(nbl->tab.n+1)*sizeof(real)+16);
  
  /* align it - size_t has the same same as a pointer */
  nbl->vdwtab = (real *) (((size_t) p_tmp + 16) & (~((size_t) 15)));  

  /* 4 fp entries per coul table point, n+1 points, and 16 bytes extra to align it. */
  p_tmp = malloc(4*(nbl->tab.n+1)*sizeof(real)+16);
  
  /* align it - size_t has the same same as a pointer */
  nbl->coultab = (real *) (((size_t) p_tmp + 16) & (~((size_t) 15)));  

  
  for(i=0; i<=nbl->tab.n; i++) {
    for(j=0; j<4; j++)
      nbl->coultab[4*i+j] = nbl->tab.tab[12*i+j];
    for(j=0; j<8; j++)
      nbl->vdwtab [8*i+j] = nbl->tab.tab[12*i+4+j];
  }
}

static void count_tables(int ftype1,int ftype2,const gmx_mtop_t *mtop,
                         int *ncount,int **count)
{
    const gmx_moltype_t *molt;
    const t_ilist *il;
    int mt,ftype,stride,i,j,tabnr;
    
    for(mt=0; mt<mtop->nmoltype; mt++)
    {
        molt = &mtop->moltype[mt];
        for(ftype=0; ftype<F_NRE; ftype++)
        {
            if (ftype == ftype1 || ftype == ftype2) {
                il = &molt->ilist[ftype];
                stride = 1 + NRAL(ftype);
                for(i=0; i<il->nr; i+=stride) {
                    tabnr = mtop->ffparams.iparams[il->iatoms[i]].tab.table;
                    if (tabnr < 0)
                        gmx_fatal(FARGS,"A bonded table number is smaller than 0: %d\n",tabnr);
                    if (tabnr >= *ncount) {
                        srenew(*count,tabnr+1);
                        for(j=*ncount; j<tabnr+1; j++)
                            (*count)[j] = 0;
                        *ncount = tabnr+1;
                    }
                    (*count)[tabnr]++;
                }
            }
        }
    }
}

static bondedtable_t *make_bonded_tables(FILE *fplog,
                                         int ftype1,int ftype2,
                                         const gmx_mtop_t *mtop,
                                         const char *basefn,const char *tabext)
{
    int  i,ncount,*count;
    char tabfn[STRLEN];
    bondedtable_t *tab;
    
    tab = NULL;
    
    ncount = 0;
    count = NULL;
    count_tables(ftype1,ftype2,mtop,&ncount,&count);
    
    if (ncount > 0) {
        snew(tab,ncount);
        for(i=0; i<ncount; i++) {
            if (count[i] > 0) {
                sprintf(tabfn,"%s",basefn);
                sprintf(tabfn + strlen(basefn) - strlen(ftp2ext(efXVG)) - 1,"_%s%d.%s",
                        tabext,i,ftp2ext(efXVG));
                tab[i] = make_bonded_table(fplog,tabfn,NRAL(ftype1)-2);
            }
        }
        sfree(count);
    }
  
    return tab;
}

void init_forcerec(FILE *fp,
                   t_forcerec *fr,
                   t_fcdata   *fcd,
                   const t_inputrec *ir,
                   const gmx_mtop_t *mtop,
                   const t_commrec  *cr,
                   matrix     box,
                   bool       bMolEpot,
                   const char *tabfn,
                   const char *tabpfn,
                   const char *tabbfn,
                   bool       bNoSolvOpt,
                   real       print_force)
{
    int     i,j,m,natoms,ngrp,negp_pp,negptable,egi,egj;
    real    rtab;
    rvec    box_size;
    const t_block *cgs;
    bool    bTab,bSep14tab,bNormalnblists;
    t_nblists *nbl;
    int     *nm_ind,egp_flags;
    
    fr->bDomDec = DOMAINDECOMP(cr);

    natoms = mtop->natoms;

    if (check_box(ir->ePBC,box))
    {
        gmx_fatal(FARGS,check_box(ir->ePBC,box));
    }
    
    /* Test particle insertion ? */
    if (EI_TPI(ir->eI)) {
        /* Set to the size of the molecule to be inserted (the last one) */
        /* Because of old style topologies, we have to use the last cg
         * instead of the last molecule type.
         */
        cgs = &mtop->moltype[mtop->molblock[mtop->nmolblock-1].type].cgs;
        fr->n_tpi = cgs->index[cgs->nr] - cgs->index[cgs->nr-1];
        if (fr->n_tpi != mtop->mols.index[mtop->mols.nr] - mtop->mols.index[mtop->mols.nr-1]) {
            gmx_fatal(FARGS,"The molecule to insert can not consist of multiple charge groups.\nMake it a single charge group.");
        }
    } else {
        fr->n_tpi = 0;
    }
    
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
    fr->sc_power   = ir->sc_power;
    fr->sc_sigma6  = pow(ir->sc_sigma,6);
    
    /* Neighbour searching stuff */
    fr->bGrid      = (ir->ns_type == ensGRID);
    fr->ePBC       = ir->ePBC;
    fr->bMolPBC    = ir->bPeriodicMols;
    fr->rc_scaling = ir->refcoord_scaling;
    copy_rvec(ir->posres_com,fr->posres_com);
    copy_rvec(ir->posres_comB,fr->posres_comB);
    fr->rlist      = ir->rlist;
    fr->rlistlong  = max(ir->rlist,max(ir->rcoulomb,ir->rvdw));
    fr->eeltype    = ir->coulombtype;
    fr->vdwtype    = ir->vdwtype;
    
    fr->bTwinRange = fr->rlistlong > fr->rlist;
    fr->bEwald     = (EEL_PME(fr->eeltype) || fr->eeltype==eelEWALD);
    
    fr->bvdwtab    = (fr->vdwtype != evdwCUT);
    fr->bcoultab   = (!(fr->eeltype == eelCUT || EEL_RF(fr->eeltype)) ||
                      fr->eeltype == eelRF_ZERO);
    
    if (getenv("GMX_FORCE_TABLES"))
    {
        fr->bvdwtab  = TRUE;
        fr->bcoultab = TRUE;
    }
    
    if (fp) {
        fprintf(fp,"Table routines are used for coulomb: %s\n",bool_names[fr->bcoultab]);
        fprintf(fp,"Table routines are used for vdw:     %s\n",bool_names[fr->bvdwtab ]);
    }
    
    /* Tables are used for direct ewald sum */
    if(fr->bEwald)
    {
        if (EEL_PME(ir->coulombtype))
        {
            if (fp)
                fprintf(fp,"Will do PME sum in reciprocal space.\n");
            please_cite(fp,"Essman95a");
            
            if (ir->ewald_geometry == eewg3DC)
            {
                if (fp)
                {
                    fprintf(fp,"Using the Ewald3DC correction for systems with a slab geometry.\n");
                }
                please_cite(fp,"In-Chul99a");
            }
        }
        fr->ewaldcoeff=calc_ewaldcoeff(ir->rcoulomb, ir->ewald_rtol);
        if (fp)
        {
            fprintf(fp,"Using a Gaussian width (1/beta) of %g nm for Ewald\n",
                    1/fr->ewaldcoeff);
        }
    }
    
    /* Electrostatics */
    fr->epsilon_r  = ir->epsilon_r;
    fr->epsilon_rf = ir->epsilon_rf;
    fr->fudgeQQ    = mtop->ffparams.fudgeQQ;
    fr->rcoulomb_switch = ir->rcoulomb_switch;
    fr->rcoulomb        = ir->rcoulomb;
    
    /* Parameters for generalized RF */
    fr->zsquare = 0.0;
    fr->temp    = 0.0;
    
    if (fr->eeltype == eelGRF)
    {
        init_generalized_rf(fp,mtop,ir,fr);
    }
    else if (EEL_FULL(fr->eeltype) || (fr->eeltype == eelSHIFT) || 
             (fr->eeltype == eelUSER) || (fr->eeltype == eelSWITCH))
    {
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
        
        if (fr->eeltype == eelPPPM && fr->phi == NULL)
            snew(fr->phi,natoms);
        
        if ((fr->eeltype==eelPPPM) || (fr->eeltype==eelPOISSON) || 
            (fr->eeltype == eelSHIFT && fr->rcoulomb > fr->rcoulomb_switch))
            set_shift_consts(fp,fr->rcoulomb_switch,fr->rcoulomb,box_size,fr);
    }
    
    /* Initiate arrays */
    if (fr->bTwinRange) {
        if (!DOMAINDECOMP(cr)) {
            fr->f_twin_n = natoms;
            fr->f_twin_nalloc = fr->f_twin_n;
            snew(fr->f_twin,fr->f_twin_nalloc);
        }
        snew(fr->fshift_twin,SHIFTS);
    }

    fr->bF_NoVirSum = (EEL_FULL(fr->eeltype) ||
                       gmx_mtop_ftype_count(mtop,F_POSRES) > 0);

    if (fr->bF_NoVirSum) {
        if (!DOMAINDECOMP(cr)) {
            fr->f_novirsum_n = natoms;
            fr->f_novirsum_nalloc = fr->f_novirsum_n;
            snew(fr->f_novirsum,fr->f_novirsum_nalloc);
        }
    }
    
    /* Mask that says whether or not this NBF list should be computed */
    /*  if (fr->bMask == NULL) {
        ngrp = ir->opts.ngener*ir->opts.ngener;
        snew(fr->bMask,ngrp);*/
    /* Defaults to always */
    /*    for(i=0; (i<ngrp); i++)
          fr->bMask[i] = TRUE;
          }*/
    
    if (fr->cg_cm == NULL && !DOMAINDECOMP(cr)) {
        /* Count the total number of charge groups */
        fr->cg_nalloc = ncg_mtop(mtop);
        snew(fr->cg_cm,fr->cg_nalloc);
    }
    if (fr->shift_vec == NULL)
        snew(fr->shift_vec,SHIFTS);
    
    if (fr->fshift == NULL)
        snew(fr->fshift,SHIFTS);
    
    if (fr->nbfp == NULL) {
        fr->ntype = mtop->ffparams.atnr;
        fr->bBHAM = (mtop->ffparams.functype[0] == F_BHAM);
        fr->nbfp  = mk_nbfp(&mtop->ffparams,fr->bBHAM);
    }
    
    /* Copy the energy group exclusions */
    fr->egp_flags = ir->opts.egp_flags;
    
    /* Van der Waals stuff */
    fr->rvdw        = ir->rvdw;
    fr->rvdw_switch = ir->rvdw_switch;
    if ((fr->vdwtype != evdwCUT) && (fr->vdwtype != evdwUSER) && !fr->bBHAM) {
        if (fr->rvdw_switch >= fr->rvdw)
            gmx_fatal(FARGS,"rvdw_switch (%f) must be < rvdw (%f)",
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
    {
        set_avcsixtwelve(fp,fr,mtop);
    }
    
    if (fr->bBHAM)
    {
        set_bham_b_max(fp,fr,mtop);
    }
    
    /* Copy the GBSA data (radius, volume and surftens for each
     * atomtype) from the topology atomtype section to forcerec.
     */
    snew(fr->atype_radius,fr->ntype);
    snew(fr->atype_vol,fr->ntype);
    snew(fr->atype_surftens,fr->ntype);
    if (mtop->atomtypes.nr > 0)
    {
        for(i=0;i<fr->ntype;i++)
            fr->atype_radius[i] =mtop->atomtypes.radius[i];
        for(i=0;i<fr->ntype;i++)
            fr->atype_vol[i] = mtop->atomtypes.vol[i];
        for(i=0;i<fr->ntype;i++)
            fr->atype_surftens[i] = mtop->atomtypes.surftens[i];
    }    
    
    /* Set the charge scaling */
    if (fr->epsilon_r != 0)
        fr->epsfac = ONE_4PI_EPS0/fr->epsilon_r;
    else
        /* eps = 0 is infinite dieletric: no coulomb interactions */
        fr->epsfac = 0;
    
    /* Reaction field constants */
    if (EEL_RF(fr->eeltype))
        calc_rffac(fp,fr->eeltype,fr->epsilon_r,fr->epsilon_rf,
                   fr->rcoulomb,fr->temp,fr->zsquare,box,
                   &fr->kappa,&fr->k_rf,&fr->c_rf);
    
    set_chargesum(fp,fr,mtop);
    
    /* if we are using LR electrostatics, and they are tabulated,
     * the tables will contain modified coulomb interactions.
     * Since we want to use the non-shifted ones for 1-4
     * coulombic interactions, we must have an extra set of tables.
     */
    
    /* Construct tables.
     * A little unnecessary to make both vdw and coul tables sometimes,
     * but what the heck... */
    
    bTab = fr->bcoultab || fr->bvdwtab;

    bSep14tab = ((!bTab || fr->eeltype!=eelCUT || fr->vdwtype!=evdwCUT ||
                  fr->bBHAM) &&
                 (gmx_mtop_ftype_count(mtop,F_LJ14) > 0 ||
                  gmx_mtop_ftype_count(mtop,F_LJC14_Q) > 0 ||
                  gmx_mtop_ftype_count(mtop,F_LJC_PAIRS_NB) > 0));

    negp_pp = ir->opts.ngener - ir->nwall;
    negptable = 0;
    if (!bTab) {
        bNormalnblists = TRUE;
        fr->nnblists = 1;
    } else {
        bNormalnblists = (ir->eDispCorr != edispcNO);
        for(egi=0; egi<negp_pp; egi++) {
            for(egj=egi;  egj<negp_pp; egj++) {
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
    
    rtab = fr->rlistlong + ir->tabext;
    if (bTab) {
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
            nm_ind = mtop->groups.grps[egcENER].nm_ind;
            for(egi=0; egi<negp_pp; egi++) {
                for(egj=egi;  egj<negp_pp; egj++) {
                    egp_flags = ir->opts.egp_flags[GID(egi,egj,ir->opts.ngener)];
                    if ((egp_flags & EGP_TABLE) && !(egp_flags & EGP_EXCL)) {
                        nbl = &(fr->nblists[m]);
                        if (fr->nnblists > 1) {
                            fr->gid2nblists[GID(egi,egj,ir->opts.ngener)] = m;
                        }
                        /* Read the table file with the two energy groups names appended */
                        make_nbf_tables(fp,fr,rtab,cr,tabfn,
                                        *mtop->groups.grpname[nm_ind[egi]],
                                        *mtop->groups.grpname[nm_ind[egj]],
                                        &fr->nblists[m]);
                        m++;
                    } else if (fr->nnblists > 1) {
                        fr->gid2nblists[GID(egi,egj,ir->opts.ngener)] = 0;
                    }
                }
            }
        }
    }
    if (bSep14tab)
        /* generate extra tables with plain Coulomb for 1-4 interactions only */
        fr->tab14 = make_tables(fp,fr,MASTER(cr),tabpfn,rtab,FALSE,TRUE);
    
    /* Wall stuff */
    fr->nwall = ir->nwall;
    if (ir->nwall && ir->wall_type==ewtTABLE)
    {
        make_wall_tables(fp,ir,tabfn,&mtop->groups,fr);
    }
    
    if (fcd && tabbfn) {
        fcd->bondtab  = make_bonded_tables(fp,
                                           F_TABBONDS,F_TABBONDSNC,
                                           mtop,tabbfn,"b");
        fcd->angletab = make_bonded_tables(fp,
                                           F_TABANGLES,-1,
                                           mtop,tabbfn,"a");
        fcd->dihtab   = make_bonded_tables(fp,
                                           F_TABDIHS,-1,
                                           mtop,tabbfn,"d");
    } else {
        if (debug)
            fprintf(debug,"No fcdata or table file name passed, can not read table, can not do bonded interactions\n");
    }
    
    /* QM/MM initialization if requested
     */
    if (ir->bQMMM)
    {
        fprintf(stderr,"QM/MM calculation requested.\n");
    }
    
    fr->bQMMM      = ir->bQMMM;   
    fr->qr         = mk_QMMMrec();
    
    /* Set all the static charge group info */
    fr->cginfo_global = init_cginfo(fp,mtop,fr,bNoSolvOpt);
    if (DOMAINDECOMP(cr)) {
        fr->cginfo = NULL;
    } else {
        fr->cginfo = fr->cginfo_global;
    }
    
    /* This is corrected later for domain and particle decomposition */
    fr->cg0 = 0;
    fr->hcg = ncg_mtop(mtop);
    
    fr->print_force = print_force;
    
    /* Initialize neighbor search */
    init_ns(fp,cr,&fr->ns,fr,mtop,box);
    
    if (cr->duty & DUTY_PP)
        gmx_setup_kernels(fp);
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
  pr_bool(fp,fr->bGrid);
  pr_bool(fp,fr->bTwinRange);
  /*pr_int(fp,fr->cg0);
    pr_int(fp,fr->hcg);*/
  for(i=0; i<fr->nnblists; i++)
    pr_int(fp,fr->nblists[i].tab.n);
  pr_real(fp,fr->rcoulomb_switch);
  pr_real(fp,fr->rcoulomb);
  
  fflush(fp);
}

void ns(FILE *fp,
        t_forcerec *fr,
        rvec       x[],
        rvec       f[],
        matrix     box,
        gmx_groups_t *groups,
        t_grpopts  *opts,
        gmx_localtop_t *top,
        t_mdatoms  *md,
        t_commrec  *cr,
        t_nrnb     *nrnb,
        int        step,
        real       lambda,
        real       *dvdlambda,
        gmx_grppairener_t *grppener,
        bool       bFillGrid,
        bool       bDoForces)
{
  static bool bFirst=TRUE;
  static int  nDNL;
  char   *ptr;
  int    nsearch;

  GMX_MPE_LOG(ev_ns_start);

  if (bFirst) {
    ptr=getenv("DUMPNL");
    if (ptr) {
      nDNL = atoi(ptr);
      if (fp) {
	fprintf(fp,"nDNL = %d\n",nDNL);  
      }
    } else
      nDNL=0;
    /* Allocate memory for the neighbor lists */
    init_neighbor_list(fp,fr,md->homenr);
      
    bFirst=FALSE;
  }
    
  if (fr->bTwinRange) 
    fr->nlr=0;

  nsearch = search_neighbours(fp,fr,x,box,top,groups,cr,nrnb,md,
                              lambda,dvdlambda,grppener,
                              bFillGrid,bDoForces);
  if (debug)
    fprintf(debug,"nsearch = %d\n",nsearch);
    
  /* Check whether we have to do dynamic load balancing */
  /*if ((nsb->nstDlb > 0) && (mod(step,nsb->nstDlb) == 0))
    count_nb(cr,nsb,&(top->blocks[ebCGS]),nns,fr->nlr,
    &(top->idef),opts->ngener);
  */
  if (nDNL > 0)
    dump_nblist(fp,cr,fr,nDNL);

  GMX_MPE_LOG(ev_ns_finish);
}

void do_force_lowlevel(FILE       *fplog,   int        step,
                       t_forcerec *fr,      t_inputrec *ir,
                       t_idef     *idef,    t_commrec  *cr,
                       t_nrnb     *nrnb,    gmx_wallcycle_t wcycle,
                       t_mdatoms  *md,
                       t_grpopts  *opts,
                       rvec       x[],      history_t  *hist,
                       rvec       f[],
                       gmx_enerdata_t *enerd,
                       t_fcdata   *fcd,
                       matrix     box,
                       real       lambda,   t_graph    *graph,
                       t_blocka   *excl,    
                       rvec       mu_tot[],
                       int        flags,
                       float      *cycles_force)
{
    int     i,nit,status;
    bool    bDoEpot,bSepDVDL,bSB;
    matrix  boxs;
    rvec    box_size;
    real    dvdlambda,Vsr,Vlr,Vcorr=0,vdip,vcharge;
    t_pbc   pbc;
#ifdef GMX_MPI
    double  t0=0.0,t1,t2,t3; /* time measurement for coarse load balancing */
#endif
    static double  t_fnbf=0.0, t_wait=0.0;
    static int     timesteps=0;
    
#define PRINT_SEPDVDL(s,v,dvdl) if (bSepDVDL) fprintf(fplog,sepdvdlformat,s,v,dvdl);
    
    GMX_MPE_LOG(ev_force_start);
    
    /* Reset box */
    for(i=0; (i<DIM); i++)
    {
        box_size[i]=box[i][i];
    }
    
    bSepDVDL=(fr->bSepDVDL && do_per_step(step,ir->nstlog));
    debug_gmx();
    
    /* do QMMM first if requested */
    if(fr->bQMMM)
    {
        enerd->term[F_EQM] = calculate_QMMM(cr,x,f,fr,md);
    }
    
    if (bSepDVDL)
    {
        fprintf(fplog,"Step %d: non-bonded V and dVdl for node %d:\n",
                step,cr->nodeid);
    }
    
    /* Call the short range functions all in one go. */
    GMX_MPE_LOG(ev_do_fnbf_start);
    
    dvdlambda = 0;
    
#ifdef GMX_MPI
    /*#define TAKETIME ((cr->npmenodes) && (timesteps < 12))*/
#define TAKETIME FALSE
    if (TAKETIME)
    {
        MPI_Barrier(cr->mpi_comm_mygroup);
        t0=MPI_Wtime();
    }
#endif
    
    if (ir->nwall)
    {
        dvdlambda = do_walls(ir,fr,box,md,x,f,lambda,
                             enerd->grpp.ener[egLJSR],nrnb);
        PRINT_SEPDVDL("Walls",0.0,dvdlambda);
        enerd->term[F_DVDL] += dvdlambda;
    }
    
    where();
    do_nonbonded(cr,fr,x,f,md,
                 fr->bBHAM ?
                 enerd->grpp.ener[egBHAMSR] :
                 enerd->grpp.ener[egLJSR],
                 enerd->grpp.ener[egCOULSR],box_size,nrnb,
                 lambda,&dvdlambda,FALSE,-1,-1,
                 flags & GMX_FORCE_FORCES);
    where();
    
#ifdef GMX_MPI
    if (TAKETIME)
    {
        t1=MPI_Wtime();
        t_fnbf += t1-t0;
    }
#endif
    
    enerd->term[F_DVDL] += dvdlambda;
    Vsr = 0;
    if (bSepDVDL)
    {
        for(i=0; i<enerd->grpp.nener; i++)
        {
            Vsr +=
                (fr->bBHAM ?
                 enerd->grpp.ener[egBHAMSR][i] :
                 enerd->grpp.ener[egLJSR][i])
                + enerd->grpp.ener[egCOULSR][i];
        }
    }
    PRINT_SEPDVDL("VdW and Coulomb SR particle-p.",Vsr,dvdlambda);
    debug_gmx();
    
    GMX_MPE_LOG(ev_do_fnbf_finish);
    
    if (debug)
    {
        pr_rvecs(debug,0,"fshift after SR",fr->fshift,SHIFTS);
    }
    
    /* Shift the coordinates. Must be done before bonded forces and PPPM, 
     * but is also necessary for SHAKE and update, therefore it can NOT 
     * go when no bonded forces have to be evaluated.
     */
    
    /* Here sometimes we would not need to shift with NBFonly,
     * but we do so anyhow for consistency of the returned coordinates.
     */
    if (graph)
    {
        shift_self(graph,box,x);
        if (TRICLINIC(box))
        {
            inc_nrnb(nrnb,eNR_SHIFTX,2*graph->nnodes);
        }
        else
        {
            inc_nrnb(nrnb,eNR_SHIFTX,graph->nnodes);
        }
    }
    /* Check whether we need to do bondeds or correct for exclusions */
    if (fr->bMolPBC &&
        ((flags & GMX_FORCE_BONDED)
         || EEL_RF(fr->eeltype) || EEL_FULL(fr->eeltype)))
    {
        /* Since all atoms are in the rectangular or triclinic unit-cell,
         * only single box vector shifts (2 in x) are required.
         */
        set_pbc_dd(&pbc,fr->ePBC,cr->dd,TRUE,box);
    }
    debug_gmx();
    
    if (flags & GMX_FORCE_BONDED)
    {
        GMX_MPE_LOG(ev_calc_bonds_start);
        calc_bonds(fplog,cr->ms,
                   idef,x,hist,f,fr,&pbc,graph,enerd,nrnb,lambda,md,fcd,
                   DOMAINDECOMP(cr) ? cr->dd->gatindex : NULL,
                   fr->bSepDVDL && do_per_step(step,ir->nstlog),step);
        debug_gmx();
        GMX_MPE_LOG(ev_calc_bonds_finish);
    }
    
    where();
    if (EEL_FULL(fr->eeltype))
    {
        bSB = (ir->nwall == 2);
        if (bSB)
        {
            copy_mat(box,boxs);
            svmul(ir->wall_ewald_zfac,boxs[ZZ],boxs[ZZ]);
            box_size[ZZ] *= ir->wall_ewald_zfac;
        }
        
        clear_mat(fr->vir_el_recip);	
        
        if(fr->bEwald)
        {
            dvdlambda = 0;
            Vcorr = ewald_LRcorrection(fplog,md->start,md->start+md->homenr,
                                       cr,fr,
                                       md->chargeA,
                                       md->nChargePerturbed ? md->chargeB : NULL,
                                       excl,x,bSB ? boxs : box,mu_tot,
                                       ir->ewald_geometry,ir->epsilon_surface,
                                       lambda,&dvdlambda,&vdip,&vcharge);
            PRINT_SEPDVDL("Ewald excl./charge/dip. corr.",Vcorr,dvdlambda);
            enerd->term[F_DVDL] += dvdlambda;
        }
        else
        {
            Vcorr = shift_LRcorrection(fplog,md->start,md->homenr,cr,fr,
                                       md->chargeA,excl,x,TRUE,box,
                                       fr->vir_el_recip);
        }
        
        *cycles_force = wallcycle_stop(wcycle,ewcFORCE);
        /* Now we can do communication again */
        
        dvdlambda = 0;
        status = 0;
        switch (fr->eeltype)
        {
        case eelPPPM:
            status = gmx_pppm_do(fplog,fr->pmedata,FALSE,x,fr->f_novirsum,
                                 md->chargeA,
                                 box_size,fr->phi,cr,md->start,md->homenr,
                                 nrnb,ir->pme_order,&Vlr);
            break;
        case eelPME:
        case eelPMESWITCH:
        case eelPMEUSER:
            if (cr->duty & DUTY_PME)
            {
                wallcycle_start(wcycle,ewcPMEMESH);
                status = gmx_pme_do(fr->pmedata,
                                    md->start,md->homenr,
                                    x,fr->f_novirsum,
                                    md->chargeA,md->chargeB,
                                    bSB ? boxs : box,cr,
                                    DOMAINDECOMP(cr) ? dd_pme_maxshift(cr->dd) : 0,
                                    nrnb,
                                    fr->vir_el_recip,fr->ewaldcoeff,
                                    &Vlr,lambda,&dvdlambda,FALSE);
                PRINT_SEPDVDL("PME mesh",Vlr,dvdlambda);
                wallcycle_stop(wcycle,ewcPMEMESH);
            } 
            else
            {
                /* Energies and virial are obtained later from the PME nodes */
                /* but values have to be zeroed out here */
                Vlr=0.0;
            }
            break;
        case eelEWALD:
            Vlr = do_ewald(fplog,FALSE,ir,x,fr->f_novirsum,
                           md->chargeA,md->chargeB,
                           box_size,cr,md->homenr,
                           fr->vir_el_recip,fr->ewaldcoeff,
                           lambda,&dvdlambda);
            PRINT_SEPDVDL("Ewald long-range",Vlr,dvdlambda);
            break;
        default:
            Vlr = 0;
            gmx_fatal(FARGS,"No such electrostatics method implemented %s",
                      eel_names[fr->eeltype]);
        }
        if (status != 0)
        {
            gmx_fatal(FARGS,"Error %d in long range electrostatics routine %s",
                      status,EELTYPE(fr->eeltype));
		}
        enerd->term[F_DVDL] += dvdlambda;
        enerd->term[F_COUL_RECIP] = Vlr + Vcorr;
        if (debug)
        {
            fprintf(debug,"Vlr = %g, Vcorr = %g, Vlr_corr = %g\n",
                    Vlr,Vcorr,enerd->term[F_COUL_RECIP]);
            pr_rvecs(debug,0,"vir_el_recip after corr",fr->vir_el_recip,DIM);
            pr_rvecs(debug,0,"fshift after LR Corrections",fr->fshift,SHIFTS);
        }
    }
    else
    {
        if (EEL_RF(fr->eeltype))
        {
            dvdlambda = 0;
            
            if (fr->eeltype != eelRF_NEC)
            {
                enerd->term[F_RF_EXCL] =
                    RF_excl_correction(fplog,fr,graph,md,excl,x,f,
                                       fr->fshift,&pbc,lambda,&dvdlambda);
            }
            
            enerd->term[F_DVDL] += dvdlambda;
            PRINT_SEPDVDL("RF exclusion correction",
                          enerd->term[F_RF_EXCL],dvdlambda);
        }
        
        *cycles_force = wallcycle_stop(wcycle,ewcFORCE);
        /* Now we can do communication again */
    }
    where();
    debug_gmx();
    
    if (debug)
    {
        print_nrnb(debug,nrnb); 
    }
    debug_gmx();
    
#ifdef GMX_MPI
    if (TAKETIME)
    {
        t2=MPI_Wtime();
        MPI_Barrier(cr->mpi_comm_mygroup);
        t3=MPI_Wtime();
        t_wait += t3-t2;
        if (timesteps == 11)
        {
            fprintf(stderr,"* PP load balancing info: node %d, step %d, rel wait time=%3.0f%% , load string value: %7.2f\n", 
                    cr->nodeid, timesteps, 100*t_wait/(t_wait+t_fnbf), (t_fnbf+t_wait)/t_fnbf);
        }	  
        timesteps++;
    }
#endif
    
    if (debug)
    {
        pr_rvecs(debug,0,"fshift after bondeds",fr->fshift,SHIFTS);
    }
    
    GMX_MPE_LOG(ev_force_finish);
}

void init_enerdata(FILE *log,int ngener,gmx_enerdata_t *enerd)
{
    int i,n2;
    
    for(i=0; i<F_NRE; i++)
    {
        enerd->term[i] = 0;
    }
    
    n2=ngener*ngener;
#ifdef DEBUG
    fprintf(log,"Creating %d sized group matrix for energies\n",n2);
#endif
    enerd->grpp.nener = n2;
    for(i=0; (i<egNR); i++)
    {
        snew(enerd->grpp.ener[i],n2);
    }
}

