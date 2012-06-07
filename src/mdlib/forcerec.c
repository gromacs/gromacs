/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
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
#include "assert.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "maths.h"
#include "macros.h"
#include "smalloc.h"
#include "macros.h"
#include "gmx_fatal.h"
#include "gmx_fatal_collective.h"
#include "physics.h"
#include "force.h"
#include "tables.h"
#include "nonbonded.h"
#include "invblock.h"
#include "names.h"
#include "network.h"
#include "pbc.h"
#include "ns.h"
#include "mshift.h"
#include "txtdump.h"
#include "coulomb.h"
#include "mdrun.h"
#include "domdec.h"
#include "partdec.h"
#include "qmmm.h"
#include "copyrite.h"
#include "mtop_util.h"
#include "nbnxn_search.h"
#include "statutil.h"
#include "gmx_omp_nthreads.h"


#ifdef _MSC_VER
/* MSVC definition for __cpuid() */
#include <intrin.h>
#endif

#include "types/nbnxn_cuda_types_ext.h"
#include "gpu_utils.h"
#include "nbnxn_cuda_data_mgmt.h"
#include "pmalloc_cuda.h"

t_forcerec *mk_forcerec(void)
{
  t_forcerec *fr;
  
  snew(fr,1);
  
  return fr;
}

#ifdef DEBUG
static void pr_nbfp(FILE *fp,real *nbfp,gmx_bool bBHAM,int atnr)
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

static real *mk_nbfp(const gmx_ffparams_t *idef,gmx_bool bBHAM)
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
                 int                   nmol,
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
    gmx_bool              perturbed;
    gmx_bool              has_vdw[4];
    gmx_bool              match;
    real              tmp_charge[4];
    int               tmp_vdwtype[4];
    int               tjA;
    gmx_bool              qm;
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
            solvent_parameters[k].count += nmol;

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
            solvent_parameters[*n_solvent_parameters].count = nmol;
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
            solvent_parameters[*n_solvent_parameters].model = esolTIP4P;
            solvent_parameters[*n_solvent_parameters].count = nmol;
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
              cginfo_mb_t           *cginfo_mb)
{
    const t_block *   cgs;
    const t_block *   mols;
    const gmx_moltype_t *molt;
    int               mb,mol,cg_mol,at_offset,cg_offset,am,cgm,i,nmol_ch,nmol;
    int               n_solvent_parameters;
    solvent_parameters_t *solvent_parameters;
    int               **cg_sp;
    int               bestsp,bestsol;

    if (debug)
    {
        fprintf(debug,"Going to determine what solvent types we have.\n");
    }

    mols = &mtop->mols;

    n_solvent_parameters = 0;
    solvent_parameters = NULL;
    /* Allocate temporary array for solvent type */
    snew(cg_sp,mtop->nmolblock);

    cg_offset = 0;
    at_offset = 0;
    for(mb=0; mb<mtop->nmolblock; mb++)
    {
        molt = &mtop->moltype[mtop->molblock[mb].type];
        cgs  = &molt->cgs;
        /* Here we have to loop over all individual molecules
         * because we need to check for QMMM particles.
         */
        snew(cg_sp[mb],cginfo_mb[mb].cg_mod);
        nmol_ch = cginfo_mb[mb].cg_mod/cgs->nr;
        nmol    = mtop->molblock[mb].nmol/nmol_ch;
        for(mol=0; mol<nmol_ch; mol++)
        {
            cgm = mol*cgs->nr;
            am  = mol*cgs->index[cgs->nr];
            for(cg_mol=0; cg_mol<cgs->nr; cg_mol++)
            {
                check_solvent_cg(molt,cg_mol,nmol,
                                 mtop->groups.grpnr[egcQMMM] ?
                                 mtop->groups.grpnr[egcQMMM]+at_offset+am : 0,
                                 &mtop->groups.grps[egcQMMM],
                                 fr,
                                 &n_solvent_parameters,&solvent_parameters,
                                 cginfo_mb[mb].cginfo[cgm+cg_mol],
                                 &cg_sp[mb][cgm+cg_mol]);
            }
        }
        cg_offset += cgs->nr;
        at_offset += cgs->index[cgs->nr];
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
    for(mb=0; mb<mtop->nmolblock; mb++)
    {
        cgs = &mtop->moltype[mtop->molblock[mb].type].cgs;
        nmol = (mtop->molblock[mb].nmol*cgs->nr)/cginfo_mb[mb].cg_mod;
        for(i=0; i<cginfo_mb[mb].cg_mod; i++)
        {
            if (cg_sp[mb][i] == bestsp)
            {
                SET_CGINFO_SOLOPT(cginfo_mb[mb].cginfo[i],bestsol);
                fr->nWatMol += nmol;
            }
            else
            {
                SET_CGINFO_SOLOPT(cginfo_mb[mb].cginfo[i],esolNO);
            }
        }
        sfree(cg_sp[mb]);
    }
    sfree(cg_sp);
    
    if (bestsol != esolNO && fp!=NULL)
    {
        fprintf(fp,"\nEnabling %s-like water optimization for %d molecules.\n\n",
                esol_names[bestsol],
                solvent_parameters[bestsp].count);
    }

    sfree(solvent_parameters);
    fr->solvent_opt = bestsol;
}

static cginfo_mb_t *init_cginfo_mb(FILE *fplog,const gmx_mtop_t *mtop,
                                   t_forcerec *fr,gmx_bool bNoSolvOpt,
                                   gmx_bool *bExcl_IntraCGAll_InterCGNone)
{
    const t_block *cgs;
    const t_blocka *excl;
    const gmx_moltype_t *molt;
    const gmx_molblock_t *molb;
    cginfo_mb_t *cginfo_mb;
    gmx_bool *type_LJ;
    int  *cginfo;
    int  cg_offset,a_offset,cgm,am;
    int  mb,m,ncg_tot,cg,a0,a1,gid,ai,j,aj,excl_nalloc;
    int  *a_con;
    int  ia;
    gmx_bool bId,*bExcl,bExclIntraAll,bExclInter,bHaveLJ,bHaveQ;

    ncg_tot = ncg_mtop(mtop);
    snew(cginfo_mb,mtop->nmolblock);

    snew(type_LJ,fr->ntype);
    for(ai=0; ai<fr->ntype; ai++)
    {
        type_LJ[ai] = FALSE;
        for(j=0; j<fr->ntype; j++)
        {
            type_LJ[ai] = type_LJ[ai] ||
                fr->bBHAM ||
                C6(fr->nbfp,fr->ntype,ai,j) != 0 ||
                C12(fr->nbfp,fr->ntype,ai,j) != 0;
        }
    }

    *bExcl_IntraCGAll_InterCGNone = TRUE;

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

        /* Check if the cginfo is identical for all molecules in this block.
         * If so, we only need an array of the size of one molecule.
         * Otherwise we make an array of #mol times #cgs per molecule.
         */
        bId = TRUE;
        am = 0;
        for(m=0; m<molb->nmol; m++)
        {
            am = m*cgs->index[cgs->nr];
            for(cg=0; cg<cgs->nr; cg++)
            {
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];
                if (ggrpnr(&mtop->groups,egcENER,a_offset+am+a0) !=
                    ggrpnr(&mtop->groups,egcENER,a_offset   +a0))
                {
                    bId = FALSE;
                }
                if (mtop->groups.grpnr[egcQMMM] != NULL)
                {
                    for(ai=a0; ai<a1; ai++)
                    {
                        if (mtop->groups.grpnr[egcQMMM][a_offset+am+ai] !=
                            mtop->groups.grpnr[egcQMMM][a_offset   +ai])
                        {
                            bId = FALSE;
                        }
                    }
                }
            }
        }

        cginfo_mb[mb].cg_start = cg_offset;
        cginfo_mb[mb].cg_end   = cg_offset + molb->nmol*cgs->nr;
        cginfo_mb[mb].cg_mod   = (bId ? 1 : molb->nmol)*cgs->nr;
        snew(cginfo_mb[mb].cginfo,cginfo_mb[mb].cg_mod);
        cginfo = cginfo_mb[mb].cginfo;

        snew(a_con,molt->atoms.nr);
        for(ia=0; ia<molt->ilist[F_CONSTR].nr; ia+=3)
        {
            a_con[molt->ilist[F_CONSTR].iatoms[ia+1]] = 1;
            a_con[molt->ilist[F_CONSTR].iatoms[ia+2]] = 1;
        }
        for(ia=0; ia<molt->ilist[F_CONSTRNC].nr; ia+=3)
        {
            a_con[molt->ilist[F_CONSTRNC].iatoms[ia+1]] = 1;
            a_con[molt->ilist[F_CONSTRNC].iatoms[ia+2]] = 1;
        }
        for(ia=0; ia<molt->ilist[F_SETTLE].nr; ia+=4)
        {
            a_con[molt->ilist[F_SETTLE].iatoms[ia+1]] = 2;
            a_con[molt->ilist[F_SETTLE].iatoms[ia+2]] = 2;
            a_con[molt->ilist[F_SETTLE].iatoms[ia+3]] = 2;
        }

        for(m=0; m<(bId ? 1 : molb->nmol); m++)
        {
            cgm = m*cgs->nr;
            am  = m*cgs->index[cgs->nr];
            for(cg=0; cg<cgs->nr; cg++)
            {
                a0 = cgs->index[cg];
                a1 = cgs->index[cg+1];

                /* Store the energy group in cginfo */
                gid = ggrpnr(&mtop->groups,egcENER,a_offset+am+a0);
                SET_CGINFO_GID(cginfo[cgm+cg],gid);
                
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
                bHaveLJ       = FALSE;
                bHaveQ        = FALSE;
                for(ai=a0; ai<a1; ai++)
                {
                    /* Check LJ and electrostatic interactions */
                    bHaveLJ = bHaveLJ || (type_LJ[molt->atoms.atom[ai].type] ||
                                          type_LJ[molt->atoms.atom[ai].typeB]);
                    bHaveQ  = bHaveQ  || (molt->atoms.atom[ai].q != 0 ||
                                          molt->atoms.atom[ai].qB != 0);

                    /* Clear the exclusion list for atom ai */
                    for(aj=a0; aj<a1; aj++)
                    {
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

                    if (a_con[ai] == 1)
                    {
                        SET_CGINFO_CONSTR(cginfo[cgm+cg]);
                    }
                    else if (a_con[ai] == 2)
                    {
                        SET_CGINFO_SETTLE(cginfo[cgm+cg]);
                    }
                }
                if (bExclIntraAll)
                {
                    SET_CGINFO_EXCL_INTRA(cginfo[cgm+cg]);
                }
                if (bExclInter)
                {
                    SET_CGINFO_EXCL_INTER(cginfo[cgm+cg]);
                }
                if (a1 - a0 > MAX_CHARGEGROUP_SIZE)
                {
                    /* The size in cginfo is currently only read with DD */
                    gmx_fatal(FARGS,"A charge group has size %d which is larger than the limit of %d atoms",a1-a0,MAX_CHARGEGROUP_SIZE);
                }
                if (bHaveLJ)
                {
                    SET_CGINFO_HAS_LJ(cginfo[cgm+cg]);
                }
                if (bHaveQ)
                {
                    SET_CGINFO_HAS_Q(cginfo[cgm+cg]);
                }
                /* Store the charge group size */
                SET_CGINFO_NATOMS(cginfo[cgm+cg],a1-a0);

                if (!bExclIntraAll || bExclInter)
                {
                    *bExcl_IntraCGAll_InterCGNone = FALSE;
                }
            }
        }

        sfree(a_con);

        cg_offset += molb->nmol*cgs->nr;
        a_offset  += molb->nmol*cgs->index[cgs->nr];
    }
    sfree(bExcl);
    
    /* the solvent optimizer is called after the QM is initialized,
     * because we don't want to have the QM subsystemto become an
     * optimized solvent
     */

    check_solvent(fplog,mtop,fr,cginfo_mb);
    
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
        for(mb=0; mb<mtop->nmolblock; mb++)
        {
            for(cg=0; cg<cginfo_mb[mb].cg_mod; cg++)
            {
                SET_CGINFO_SOLOPT(cginfo_mb[mb].cginfo[cg],esolNO);
            }
        }
    }
    
    return cginfo_mb;
}

static int *cginfo_expand(int nmb,cginfo_mb_t *cgi_mb)
{
    int ncg,mb,cg;
    int *cginfo;

    ncg = cgi_mb[nmb-1].cg_end;
    snew(cginfo,ncg);
    mb = 0;
    for(cg=0; cg<ncg; cg++)
    {
        while (cg >= cgi_mb[mb].cg_end)
        {
            mb++;
        }
        cginfo[cg] =
            cgi_mb[mb].cginfo[(cg - cgi_mb[mb].cg_start) % cgi_mb[mb].cg_mod];
    }

    return cginfo;
}

static void set_chargesum(FILE *log,t_forcerec *fr,const gmx_mtop_t *mtop)
{
    double qsum,q2sum,q;
    int    mb,nmol,i;
    const t_atoms *atoms;
    
    qsum  = 0;
    q2sum = 0;
    for(mb=0; mb<mtop->nmolblock; mb++)
    {
        nmol  = mtop->molblock[mb].nmol;
        atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
        for(i=0; i<atoms->nr; i++)
        {
            q = atoms->atom[i].q;
            qsum  += nmol*q;
            q2sum += nmol*q*q;
        }
    }
    fr->qsum[0]  = qsum;
    fr->q2sum[0] = q2sum;
    if (fr->efep != efepNO)
    {
        qsum  = 0;
        q2sum = 0;
        for(mb=0; mb<mtop->nmolblock; mb++)
        {
            nmol  = mtop->molblock[mb].nmol;
            atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
            for(i=0; i<atoms->nr; i++)
            {
                q = atoms->atom[i].qB;
                qsum  += nmol*q;
                q2sum += nmol*q*q;
            }
            fr->qsum[1]  = qsum;
            fr->q2sum[1] = q2sum;
        }
    }
    else
    {
        fr->qsum[1]  = fr->qsum[0];
        fr->q2sum[1] = fr->q2sum[0];
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
    if (fr->eeltype == eelGRF)
    {
        calc_rffac(NULL,fr->eeltype,fr->epsilon_r,fr->epsilon_rf,
                   fr->rcoulomb,fr->temp,fr->zsquare,box,
                   &fr->kappa,&fr->k_rf,&fr->c_rf);
    }
}

void set_avcsixtwelve(FILE *fplog,t_forcerec *fr,const gmx_mtop_t *mtop)
{
    const t_atoms *atoms,*atoms_tpi;
    const t_blocka *excl;
    int    mb,nmol,nmolc,i,j,tpi,tpj,j1,j2,k,n,nexcl,q;
#if (defined SIZEOF_LONG_LONG_INT) && (SIZEOF_LONG_LONG_INT >= 8)    
    long long int  npair,npair_ij,tmpi,tmpj;
#else
    double npair, npair_ij,tmpi,tmpj;
#endif
    double csix,ctwelve;
    int    ntp,*typecount;
    gmx_bool   bBHAM;
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
            atoms_tpi =
                &mtop->moltype[mtop->molblock[mtop->nmolblock-1].type].atoms;

            npair = 0;
            for(mb=0; mb<mtop->nmolblock; mb++) {
                nmol  = mtop->molblock[mb].nmol;
                atoms = &mtop->moltype[mtop->molblock[mb].type].atoms;
                for(j=0; j<atoms->nr; j++) {
                    nmolc = nmol;
                    /* Remove the interaction of the test charge group
                     * with itself.
                     */
                    if (mb == mtop->nmolblock-1)
                    {
                        nmolc--;
                        
                        if (mb == 0 && nmol == 1)
                        {
                            gmx_fatal(FARGS,"Old format tpr with TPI, please generate a new tpr file");
                        }
                    }
                    if (q == 0)
                    {
                        tpj = atoms->atom[j].type;
                    }
                    else
                    {
                        tpj = atoms->atom[j].typeB;
                    }
                    for(i=0; i<fr->n_tpi; i++)
                    {
                        if (q == 0)
                        {
                            tpi = atoms_tpi->atom[i].type;
                        }
                        else
                        {
                            tpi = atoms_tpi->atom[i].typeB;
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
    if (fplog != NULL)
    {
        if (fr->eDispCorr == edispcAllEner ||
            fr->eDispCorr == edispcAllEnerPres)
        {
            fprintf(fplog,"Long Range LJ corr.: <C6> %10.4e, <C12> %10.4e\n",
                    fr->avcsix[0],fr->avctwelve[0]);
        }
        else
        {
            fprintf(fplog,"Long Range LJ corr.: <C6> %10.4e\n",fr->avcsix[0]);
        }
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

static void make_nbf_tables(FILE *fp,const output_env_t oenv,
                            t_forcerec *fr,real rtab,
			    const t_commrec *cr,
			    const char *tabfn,char *eg1,char *eg2,
			    t_nblists *nbl)
{
  char buf[STRLEN];
  int i,j;

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
  nbl->tab = make_tables(fp,oenv,fr,MASTER(cr),buf,rtab,0);
  /* Copy the contents of the table to separate coulomb and LJ tables too,
   * to improve cache performance.
   */

  /* For performance reasons we want
   * the table data to be aligned to 16-byte. The pointer could be freed
   * but currently isn't.
   */
  snew_aligned(nbl->vdwtab,8*(nbl->tab.n+1),16);
  snew_aligned(nbl->coultab,4*(nbl->tab.n+1),16);
  
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

void forcerec_set_ranges(t_forcerec *fr,
                         int ncg_home,int ncg_force,
                         int natoms_force,
                         int natoms_force_constr,int natoms_f_novirsum)
{
    fr->cg0 = 0;
    fr->hcg = ncg_home;

    /* fr->ncg_force is unused in the standard code,
     * but it can be useful for modified code dealing with charge groups.
     */
    fr->ncg_force           = ncg_force;
    fr->natoms_force        = natoms_force;
    fr->natoms_force_constr = natoms_force_constr;

    if (fr->natoms_force_constr > fr->nalloc_force)
    {
        fr->nalloc_force = over_alloc_dd(fr->natoms_force_constr);

        if (fr->bTwinRange)
        {
            srenew(fr->f_twin,fr->nalloc_force);
        }
    }

    if (fr->bF_NoVirSum)
    {
        fr->f_novirsum_n = natoms_f_novirsum;
        if (fr->f_novirsum_n > fr->f_novirsum_nalloc)
        {
            fr->f_novirsum_nalloc = over_alloc_dd(fr->f_novirsum_n);
            srenew(fr->f_novirsum_alloc,fr->f_novirsum_nalloc);
        }
    }
    else
    {
        fr->f_novirsum_n = 0;
    }
}

static real cutoff_inf(real cutoff)
{
    if (cutoff == 0)
    {
        cutoff = GMX_CUTOFF_INF;
    }

    return cutoff;
}

static void make_adress_tf_tables(FILE *fp,const output_env_t oenv,
                            t_forcerec *fr,const t_inputrec *ir,
			    const char *tabfn, const gmx_mtop_t *mtop,
                            matrix     box)
{
  char buf[STRLEN];
  int i,j;

  if (tabfn == NULL) {
        gmx_fatal(FARGS,"No thermoforce table file given. Use -tabletf to specify a file\n");
    return;
  }

  snew(fr->atf_tabs, ir->adress->n_tf_grps);

  for (i=0; i<ir->adress->n_tf_grps; i++){
    j = ir->adress->tf_table_index[i]; /* get energy group index */
    sprintf(buf + strlen(tabfn) - strlen(ftp2ext(efXVG)) - 1,"tf_%s.%s",
        *(mtop->groups.grpname[mtop->groups.grps[egcENER].nm_ind[j]]) ,ftp2ext(efXVG));
    printf("loading tf table for energygrp index %d from %s\n", ir->adress->tf_table_index[j], buf);
    fr->atf_tabs[i] = make_atf_table(fp,oenv,fr,buf, box);
  }

}

gmx_bool can_use_allvsall(const t_inputrec *ir, const gmx_mtop_t *mtop,
                      gmx_bool bPrintNote,t_commrec *cr,FILE *fp)
{
    gmx_bool bAllvsAll;

    bAllvsAll =
        (
         ir->rlist==0            &&
         ir->rcoulomb==0         &&
         ir->rvdw==0             &&
         ir->ePBC==epbcNONE      &&
         ir->vdwtype==evdwCUT    &&
         ir->coulombtype==eelCUT &&
         ir->efep==efepNO        &&
         (ir->implicit_solvent == eisNO || 
          (ir->implicit_solvent==eisGBSA && (ir->gb_algorithm==egbSTILL || 
                                             ir->gb_algorithm==egbHCT   || 
                                             ir->gb_algorithm==egbOBC))) &&
         getenv("GMX_NO_ALLVSALL") == NULL
            );
    
    if (bAllvsAll && ir->opts.ngener > 1)
    {
        const char *note="NOTE: Can not use all-vs-all force loops, because there are multiple energy monitor groups; you might get significantly higher performance when using only a single energy monitor group.\n";

        if (bPrintNote)
        {
            if (MASTER(cr))
            {
                fprintf(stderr,"\n%s\n",note);
            }
            if (fp != NULL)
            {
                fprintf(fp,"\n%s\n",note);
            }
        }
        bAllvsAll = FALSE;
    }

    if(bAllvsAll && fp && MASTER(cr))
    {
        fprintf(fp,"\nUsing accelerated all-vs-all kernels.\n\n");
    }
    
    return bAllvsAll;
}


/* Return TRUE if the CPU vendor is Intel, or possibly Intel and x86_64 */
static gmx_bool 
check_intel_x86_64()
{
    unsigned int level;
    unsigned int _ebx,_ecx,_edx;
    gmx_bool status;
    int CPUInfo[4];

    level = 0;
#ifdef _MSC_VER
    __cpuid(CPUInfo,0);
	
    _ebx=CPUInfo[1];
    _ecx=CPUInfo[2];
    _edx=CPUInfo[3];

#elif defined(__x86_64__)
    __asm__ ("push %%rbx\n\tcpuid\n\tmov %%rbx, %%rax\n\tpop %%rbx\n"   \
             : "=a" (_ebx), "=c" (_ecx), "=d" (_edx)                    \
             : "1" (level));
#else
    _ebx=_ecx=_edx=0;
#endif

    /* In sequence these ints as chars read "GenuineIntel" */
    status = (_ebx == 0x756e6547 &&
              _edx == 0x49656e69 &&
              _ecx == 0x6c65746e);

    return status;
}

typedef struct
{
    gmx_bool sse2;
    gmx_bool sse4_1;
    gmx_bool avx;
    gmx_bool avx_intel;
} cpuinfo_t;

/* Return 1 if SSE2 support is present, otherwise 0. */
static void 
get_cpuinfo(cpuinfo_t *cpuinfo)
{
#if ( defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) || defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE)|| defined(GMX_SSE2) )
	unsigned int level;
	unsigned int _eax,_ebx,_ecx,_edx;
	int status;
	int CPUInfo[4];
    gmx_bool intel_x86_64;

    /* This code is tricky.
     * Moving this call further down changes the results of the assembly
     * for getting the CPU vendor information.
     */
	intel_x86_64 = check_intel_x86_64();

	level = 1;
#ifdef _MSC_VER
	__cpuid(CPUInfo,1);
	
	_eax=CPUInfo[0];
	_ebx=CPUInfo[1];
	_ecx=CPUInfo[2];
	_edx=CPUInfo[3];
	
#elif defined(__x86_64__)
	/* GCC 64-bit inline asm */
	__asm__ ("push %%rbx\n\tcpuid\n\tpop %%rbx\n"                 \
			 : "=a" (_eax), "=S" (_ebx), "=c" (_ecx), "=d" (_edx) \
			 : "0" (level));
#elif defined(__i386__)
	__asm__ ("push %%ebx\n\tcpuid\n\tpop %%ebx\n"                 \
			 : "=a" (_eax), "=S" (_ebx), "=c" (_ecx), "=d" (_edx) \
			 : "0" (level));
#else
	_eax=_ebx=_ecx=_edx=0;
#endif
    
    /* Features:
     * SSE      Bit 25 of edx should be set
     * SSE2     Bit 26 of edx should be set
     * SSE3     Bit  0 of ecx should be set
     * SSE4.1   Bit 19 of ecx should be set
     * AVX      Bit 28 of ecx should be set
     */
    cpuinfo->sse2      = (_edx & (1 << 26)) != 0;
    cpuinfo->sse4_1    = (_ecx & (1 << 19)) != 0;
    cpuinfo->avx       = (_ecx & (1 << 28)) != 0;
    cpuinfo->avx_intel = (cpuinfo->avx && intel_x86_64);

    if (debug)
    {
        fprintf(debug,"cpuinfo: sse2 %d sse4_1 %d avx %d avx_intel %d\n",
                cpuinfo->sse2,cpuinfo->sse4_1,cpuinfo->avx,cpuinfo->avx_intel);
    }
#else
    cpuinfo->sse2      = FALSE;
    cpuinfo->sse4_1    = FALSE;
    cpuinfo->avx       = FALSE;
    cpuinfo->avx_intel = FALSE;
#endif
}


static void init_forcerec_f_threads(t_forcerec *fr,int grpp_nener)
{
    int t,i;

    fr->nthreads = gmx_omp_nthreads_get(emntBonded);

    snew(fr->f_t,fr->nthreads);
    /* Thread 0 uses the global force and energy arrays */
    for(t=1; t<fr->nthreads; t++)
    {
        fr->f_t[t].f = NULL;
        fr->f_t[t].f_nalloc = 0;
        snew(fr->f_t[t].fshift,SHIFTS);
        /* snew(fr->f_t[t].ener,F_NRE); */
        fr->f_t[t].grpp.nener = grpp_nener;
        for(i=0; i<egNR; i++)
        {
            snew(fr->f_t[t].grpp.ener[i],grpp_nener);
        }
    }
}

/* Parse the GMX_GPU_ID string and return the GPU ID assigned to the current node.
 * The GMX_GPU_ID should connstain a sequence of digits which represent device IDs.
 * The n-th digit in the sequence indicates the ID of the GPU assigned to the n-th
 * process/tMPI thread within a node.
 * */
static int parse_gmx_gpu_id(const t_commrec *cr)
{
    char *env;
    int  intra_id;
    char sbuf[STRLEN];

     intra_id = cr->nc.rank_intra;

    if ((env = getenv("GMX_GPU_ID")) != NULL)
    {
        /* GMX_GPU_ID too short */
        if (strlen(env) < intra_id + 1)
        {
            sprintf(sbuf, "GMX_GPU_ID too short, on node %d GPU ID #%d not specified.",
                    cr->nodeid, cr->nc.rank_intra);
            gmx_fatal(FARGS,sbuf);
        }


        if (env[intra_id] < '0' || env[intra_id] > '9')
        {
            gmx_fatal(FARGS, "Invalid character in GMX_GPU_ID string: '%c'\n",
                      env[intra_id]);
        }

        return env[intra_id] - '0';
    }
    else
    {
        return -1;
    }
}


/* Check the consistency of the GMX_GPU_ID string.

   Make sure that with thread-MPI each thread uses a different, otherwise
   threads mess up each others GPU context. */
static void check_gmx_gpu_id_consistency(const char *gpu_id_str,
                                         const t_commrec *cr)
{
    gmx_bool bAllSame = TRUE, bAllDifferent = TRUE, bSomeSame = FALSE;
    int len, i, j, npppn;
    char sbuf[STRLEN], sbuf1[STRLEN], sbuf2[STRLEN];

    len = strlen(gpu_id_str);
    if (len == 0)
    {
        gmx_fatal(FARGS, "GMX_GPU_ID environment variable set but empty.");
    }

    /* number of PP processes on my node */
    npppn = cr->nnodes_pp_intra;
    sbuf2[0] = '\0';
#ifdef GMX_THREAD_MPI
    sprintf(sbuf1, "thread-MPI thread");
    if (npppn > 1)
    {
        strcat(sbuf1, "s");
    }
#else
#ifdef GMX_MPI
    sprintf(sbuf1, "MPI process");
    if (npppn > 1)
    {
        strcat(sbuf1, "es");
    }
    sprintf(sbuf2, " per node");
#else
    /* neither MPI nor tMPI */
    sprintf(sbuf1, "process");
#endif
#endif

    if (len != npppn)
    {
        gmx_fatal(FARGS, "On node %d incorrect GMX_GPU_ID: %d device IDs passed, but %s is\n"
                  "using %d PP %s%s requiring one GPU for each.",
                  cr->nodeid, len, ShortProgram(), npppn, sbuf1, sbuf2);
    }

    if (len < 2)
    {
        return;
    }

    for (i = 0; i < len - 1; i++)
    {
        for (j = i + 1; j < len; j++)
        {
            bAllSame        &= gpu_id_str[i] == gpu_id_str[j];
            bSomeSame       |= gpu_id_str[i] == gpu_id_str[j];
            bAllDifferent   &= gpu_id_str[i] != gpu_id_str[j];
        }
    }

#ifdef GMX_THREAD_MPI
    if (!bAllDifferent)
    {
        gmx_fatal(FARGS,
                  "Invalid GPU assignment: multiple tMPI threads per GPU (GMX_GPU_ID=%s).\n"
                  "Use MPI if you are sure that you want multiple processes to use the same GPU.",
                  gpu_id_str);
    }
#endif

    if (bSomeSame)
    {
        gmx_warning("GMX_GPU_ID=%s assigns a/some GPU(s) to multiple %s;\n"
                    "         this should be avoided as it generally causes performance loss.",
                    gpu_id_str, sbuf1);
    }
}

/* TODO print some stats to the log! */
static gmx_bool init_cu_nbv(FILE *fp,const t_commrec *cr,gmx_bool forceGPU)
{
    int gpu_device_id;
    char *env;
    gmx_bool GPU_OK;
    int nnodes;
    int *gpus,ngpu,i;

#ifndef GMX_GPU
    if (forceGPU)
    {
        gmx_fatal(FARGS,"GPU requested for non-bonded interactions, "
                  "but %s was compiled without GPU support",ShortProgram());
    }

    return FALSE;
#endif

    /* cppcheck-suppress unreachableCode */
    /* unless GMX_GPU_ID is set, each process will use the GPU with id = intra-node id */
    gpu_device_id = cr->nodeid_group_intra;
    if ((env = getenv("GMX_GPU_ID")) != NULL)
    {
        gpu_device_id = parse_gmx_gpu_id(cr);

        check_gmx_gpu_id_consistency(env, cr);

        /* If you set this env.var, you want to use a GPU */
        forceGPU = TRUE;
    }

    GPU_OK = (init_gpu(fp, gpu_device_id) == 0);

    if (PAR(cr))
    {
        nnodes = cr->dd->nnodes;
    }
    else
    {
        nnodes = 1;
    }
    snew(gpus,cr->nnodes);
    gpus[cr->nodeid] = (GPU_OK ? -1 : 1+gpu_device_id);
    /* In parallel we need to check if every node has a GPU */
    if (PAR(cr))
    {
        gmx_sumi(cr->nnodes,gpus,cr);
    }
    ngpu = 0;
    for(i=0; i<cr->nnodes; i++)
    {
        if (gpus[i] == -1)
        {
            ngpu++;
        }
    }
    if (ngpu < nnodes && (ngpu > 0 || (ngpu == 0 && forceGPU)))
    {
        char *buf,buf2[40];
        
        snew(buf,40+40*cr->nnodes);
        sprintf(buf,"Could not initialize");
        for(i=0; i<cr->nnodes; i++)
        {
            if (gpus[i] > 0)
            {
                sprintf(buf2," on node %d: GPU #%d",gpus[i]-1,i);
                strcat(buf,buf2);
            }
        }
        gmx_fatal_collective(FARGS,cr,NULL,buf);
    }
    sfree(gpus);

    return (ngpu > 0);
}

static void pick_nbnxn_kernel_cpu(FILE *fp,
                                  const t_commrec *cr,
                                  gmx_bool tabulated_force,
                                  int *kernel_type)
{
#if defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) || defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE)|| defined(GMX_SSE2)
    /* by default we'll use the 4xN SSE/AVX 128-bit kernels */
    if (getenv("GMX_NOOPTIMIZEDKERNELS") == NULL)
    {
        cpuinfo_t cpuinfo;
#ifdef GMX_DOUBLE
        gmx_bool double_prec = TRUE;
#else
        gmx_bool double_prec = FALSE;
#endif

        get_cpuinfo(&cpuinfo);

        /* On Intel Sandy-Bridge AVX kernels are always faster, except for
         * tabulated forces in single precision, which are memory limited.
         * On AMD Bulldozer AVX is much slower than SSE.
         */
        if (cpuinfo.avx_intel &&
            (double_prec || !tabulated_force))
        {
#ifdef GMX_AVX
            *kernel_type = nbk4xN_X86_SIMD256;
#else
            md_print_warning(cr,fp,"NOTE: Gromacs was compiled without AVX support,\n      can not use the faster AVX kernels\n");
            *kernel_type = nbk4xN_X86_SIMD128;
#endif
        }
        else
        {
            *kernel_type = nbk4xN_X86_SIMD128;
        }

        if (getenv("GMX_NBNXN_AVX128") != NULL)
        {
            *kernel_type = nbk4xN_X86_SIMD128;
        }
        if (getenv("GMX_NBNXN_AVX256") != NULL)
        {
#ifdef GMX_AVX
            *kernel_type = nbk4xN_X86_SIMD256;
#else
            gmx_fatal(FARGS,"You requested AVX-256 nbnxn kernels, but Gromacs was built without AVX support");
#endif
        }

#ifndef GMX_SSE4_1
        if (*kernel_type == nbk4xN_X86_SIMD128 && cpuinfo.sse4_1)
        {
            md_print_warning(cr,fp,"NOTE: Gromacs was compiled without SSE4.1 support,\n      can not use the faster SSE4.1 kernels\n");
        }
#endif
    }
    else
#endif
    {
        *kernel_type = nbk4x4_PlainC;
    }
}

void pick_nbnxn_kernel(FILE *fp,
                       const t_commrec *cr,
                       gmx_bool tryGPU, gmx_bool *useGPU,
                       gmx_bool forceGPU,
                       gmx_bool tabulated_force,
                       int *kernel_type)
{
    gmx_bool emulateGPU=FALSE;

    *kernel_type = nbkNotSet;

    if (tryGPU)
    {
        /* Try to use a GPU if GMX_GPU is not defined */
        if (getenv("GMX_NO_GPU") != NULL)
        {
            if (forceGPU)
            {
                gmx_fatal(FARGS,"Conflicting command line argument and env.var. for GPU use for non-bonded interactions");
            }

            tryGPU = FALSE;

            if (MASTER(cr))
            {
                gmx_warning("GPU mode turned off by GMX_NO_GPU env var!");
            }
        }
    }

    if (tryGPU)
    {
        /* Run GPU emulation mode if GMX_EMULATE_GPU is defined and also if nobonded 
           calculations are turned off via GMX_NO_NONBONDED. This is the simple way 
       to also turn off GPU/CUDA initializations. */
        emulateGPU = ((getenv("GMX_EMULATE_GPU") != NULL) ||
                      (getenv("GMX_NO_NONBONDED") != NULL));

       *useGPU = FALSE;
    }
    if (emulateGPU)
    {
        *kernel_type = nbk8x8x8_PlainC;

        if (fp != NULL)
        {
            fprintf(fp, "Emulating GPU\n");
        }
    }    
    else if (tryGPU)
    {
        /* Try to use a GPU */
        *useGPU = init_cu_nbv(fp,cr,forceGPU);

        if (*useGPU)
        {
            *kernel_type = nbk8x8x8_CUDA;
        }
    }

    if (*kernel_type == nbkNotSet)
    {
        pick_nbnxn_kernel_cpu(fp,cr,tabulated_force,kernel_type);
    }

    if (fp != NULL)
    {
        if (MASTER(cr))
        {
            fprintf(stderr,"Using %s non-bonded kernels\n",
                    nbk_name[*kernel_type]);
        }
        fprintf(fp,"\nUsing %s non-bonded kernels\n\n",
                nbk_name[*kernel_type]);
    }
}

static void init_ewald_f_table(interaction_const_t *ic,
                               int verlet_kernel_type)
{
#define GPU_REF_EWALD_COULOMB_FORCE_TABLE_SIZE 1536

    if (nbnxn_kernel_pairlist_simple(verlet_kernel_type))
    {
        /* With a spacing of 0.0005 we are at the force summation accuracy
         * for the SSE kernels for "normal" atomistic simulations.
         */
        ic->tabq_scale = ewald_spline3_table_scale(ic->ewaldcoeff,
                                                   ic->rcoulomb);
        ic->tabq_size  = (int)(ic->rcoulomb*ic->tabq_scale) + 2;
#ifndef GMX_DOUBLE
        ic->tabq_format = tableformatFDV0;
#else
        ic->tabq_format = tableformatF;
#endif
    }
    else
    {
        /* Table size identical to the CUDA implementation */
        ic->tabq_size = GPU_REF_EWALD_COULOMB_FORCE_TABLE_SIZE;
        /* Subtract 2 iso 1 to avoid access out of range due to rounding */
        ic->tabq_scale = (ic->tabq_size - 2)/ic->rcoulomb;
        if (verlet_kernel_type == nbk8x8x8_CUDA)
        {
            /* This case is handled in the nbnxn CUDA module */
            ic->tabq_format = tableformatNONE;
        }
        else
        {
            ic->tabq_format = tableformatF;
        }
    }

    switch (ic->tabq_format)
    {
    case tableformatNONE:
        break;
    case tableformatF:
        sfree_aligned(ic->tabq_coul_F);
        sfree_aligned(ic->tabq_coul_V);
        snew_aligned(ic->tabq_coul_F,ic->tabq_size,16);
        snew_aligned(ic->tabq_coul_V,ic->tabq_size,16);
        table_spline3_fill_ewald_lr(ic->tabq_coul_F,ic->tabq_coul_V,
                                    ic->tabq_size,ic->tabq_format,
                                    1/ic->tabq_scale,ic->ewaldcoeff);
        break;
    case tableformatFDV0:
        sfree_aligned(ic->tabq_coul_F);
        snew_aligned(ic->tabq_coul_FDV0,ic->tabq_size*4,16);
        table_spline3_fill_ewald_lr(ic->tabq_coul_FDV0,NULL,
                                    ic->tabq_size,ic->tabq_format,
                                    1/ic->tabq_scale,ic->ewaldcoeff);
        break;
    default:
        gmx_incons("Unknown table format");
    }
}

void init_interaction_const_tables(FILE *fp, 
                                   interaction_const_t *ic,
                                   int verlet_kernel_type)
{
    real spacing;

    if (ic->eeltype == eelEWALD || EEL_PME(ic->eeltype))
    {
        init_ewald_f_table(ic,verlet_kernel_type);

        if (fp != NULL)
        {
            fprintf(fp,"Initialized non-bonded Ewald correction tables, spacing: %.2e size: %d\n\n",
                    1/ic->tabq_scale,ic->tabq_size);
        }
    }
}

void init_interaction_const(FILE *fp, 
                            interaction_const_t **interaction_const,
                            const t_forcerec *fr)
{
    interaction_const_t *ic;
    gmx_bool shLJ,shCoul;

    shLJ   = (getenv("GMX_NO_SHIFT_LJ") == NULL);
    shCoul = (getenv("GMX_NO_SHIFT_COUL") == NULL);

    snew(ic, 1);

    ic->rlist       = fr->rlist;

    /* Lennard-Jones */
    ic->rvdw        = fr->rvdw;
    if (shLJ)
    {
        ic->sh_invrc6 = pow(ic->rvdw,-6.0);
    }
    else
    {
        ic->sh_invrc6 = 0;
    }

    /* Electrostatics */
    ic->eeltype     = fr->eeltype;
    ic->rcoulomb    = fr->rcoulomb;
    ic->epsilon_r   = fr->epsilon_r;
    ic->epsfac      = fr->epsfac;

    /* Ewald */
    ic->ewaldcoeff  = fr->ewaldcoeff;
    if (shCoul)
    {
        ic->sh_ewald = gmx_erfc(ic->ewaldcoeff*ic->rcoulomb);
    }
    else
    {
        ic->sh_ewald = 0;
    }

    /* Reaction-field */
    if (EEL_RF(ic->eeltype))
    {
        ic->epsilon_rf = fr->epsilon_rf;
        ic->k_rf       = fr->k_rf;
        ic->c_rf       = fr->c_rf;
    }
    else
    {
        /* For plain cut-off we might use the reaction-field kernels */
        ic->epsilon_rf = ic->epsilon_r;
        ic->k_rf       = 0;
        if (shCoul)
        {
            ic->c_rf   = 1/ic->rcoulomb;
        }
        else
        {
            ic->c_rf   = 0;
        }
    }

    if (fp != NULL)
    {
        fprintf(fp,"Potential shift: LJ r^-12: %.3f r^-6 %.3f",
                sqr(ic->sh_invrc6),ic->sh_invrc6);
        if (ic->eeltype == eelCUT)
        {
            fprintf(fp,", Coulomb %.3f",ic->c_rf);
        }
        else if (EEL_PME(ic->eeltype))
        {
            fprintf(fp,", Ewald %.3e",ic->sh_ewald);
        }
        fprintf(fp,"\n");
    }

    *interaction_const = ic;

    if (fr->nbv != NULL && fr->nbv->useGPU)
    {
        nbnxn_cuda_init_const(fr->nbv->cu_nbv, ic, fr->nbv);
    }

    if (fr->cutoff_scheme == ecutsVERLET)
    {
        assert(fr->nbv != NULL && fr->nbv->grp != NULL);
        init_interaction_const_tables(fp,ic,fr->nbv->grp[fr->nbv->nloc-1].kernel_type);
    }
}

gmx_bool nb_kernel_pmetune_support(const nonbonded_verlet_t *nbv)
{
    /* PME tuning is currently only supported with only CUDA kernels */
    return (nbv->useGPU &&
            !(nbv->nloc == 2 && nbnxn_kernel_pairlist_simple(nbv->grp[1].kernel_type)));
}

static void init_nb_verlet(FILE *fp, 
                           nonbonded_verlet_t **nb_verlet,
                           const t_inputrec *ir,
                           const t_forcerec *fr,
                           const t_commrec *cr,
                           const char *nbpu_opt,
                           int nbnxn_kernel_preset)
{
    nonbonded_verlet_t *nbv;
    int  i;
    char *env;

    gmx_nbat_alloc_t *nb_alloc;
    gmx_nbat_free_t  *nb_free;

    snew(nbv, 1);

    nbv->nbs                     = NULL;

    nbv->nloc = (DOMAINDECOMP(cr) ? 2 : 1);
    for(i=0; i<nbv->nloc; i++)
    {
        nbv->grp[i].nbl_lists.nnbl = 0;
        nbv->grp[i].nbat           = NULL;
        nbv->grp[i].kernel_type    = nbkNotSet;

        if (i == 0)
        {
            if (nbnxn_kernel_preset < 0)
            {
                pick_nbnxn_kernel(fp, cr,
                                  (nbpu_opt != NULL &&
                                   (strncmp(nbpu_opt,"gpu",3) == 0 ||
                                    strcmp(nbpu_opt,"auto") == 0)),
                                  &nbv->useGPU,
                                  (nbpu_opt != NULL &&
                                   strncmp(nbpu_opt,"gpu",3) == 0),
                                  EEL_FULL(ir->coulombtype),
                                  &nbv->grp[i].kernel_type);
            }
            else
            {
                nbv->useGPU = (nbnxn_kernel_preset == nbk8x8x8_CUDA);
                nbv->grp[i].kernel_type = nbnxn_kernel_preset;
            }
        }
        else
        {
            if (nbpu_opt != NULL && strcmp(nbpu_opt,"gpu_cpu") == 0)
            {
                /* Use GPU for local, select a CPU kernel for non-local */
                pick_nbnxn_kernel(fp, cr, FALSE, NULL, FALSE, FALSE,
                                  &nbv->grp[i].kernel_type);
            }
            else
            {
                /* Use the same kernel for local and non-local interactions */
                nbv->grp[i].kernel_type = nbv->grp[0].kernel_type;
            }
        }
    }

    if (nbv->useGPU)
    {
        /* need to move this  */
        nbnxn_cuda_init(fp, &(nbv->cu_nbv), DOMAINDECOMP(cr));
        env = getenv("GMX_NB_MIN_CI");
        if (env)
        {
            sscanf(env, "%d", &nbv->min_ci_balanced);
            if (debug)
            {
                fprintf(debug, "Neighbor-list balancing parameter: %d (passed as env. var.)\n", 
                        nbv->min_ci_balanced);
            }
        }
        else
        {
            nbv->min_ci_balanced = nbnxn_cuda_min_ci_balanced(nbv->cu_nbv);
            if (debug)
            {
                fprintf(debug, "Neighbor-list balancing parameter: %d (auto-adjusted to the number of GPU multi-processors)\n",
                        nbv->min_ci_balanced);
            }
        }
    }
    else
    {
        nbv->min_ci_balanced = 0;
    }

    *nb_verlet = nbv;

    nbnxn_init_search(&nbv->nbs,
                      DOMAINDECOMP(cr) ? & cr->dd->nc : NULL,
                      DOMAINDECOMP(cr) ? domdec_zones(cr->dd) : NULL,
                      gmx_omp_nthreads_get(emntNonbonded));

    for(i=0; i<nbv->nloc; i++)
    {
        if (nbv->grp[0].kernel_type == nbk8x8x8_CUDA)
        {
            nb_alloc = &pmalloc;
            nb_free  = &pfree;
        }
        else
        {
            nb_alloc = NULL;
            nb_free  = NULL;
        }

        nbnxn_init_pairlist_set(&nbv->grp[i].nbl_lists,
                                nbnxn_kernel_pairlist_simple(nbv->grp[i].kernel_type),
                                /* 8x8x8 "non-simple" lists are ATM always combined */
                                !nbnxn_kernel_pairlist_simple(nbv->grp[i].kernel_type),
                                nb_alloc, nb_free);

        if (i == 0 ||
            nbv->grp[0].kernel_type != nbv->grp[i].kernel_type)
        {
            snew(nbv->grp[i].nbat,1);
            nbnxn_atomdata_init(fp,
                                nbv->grp[i].nbat,
                                nbv->grp[i].kernel_type,
                                fr->ntype,fr->nbfp,
                                ir->opts.ngener,
                                nbnxn_kernel_pairlist_simple(nbv->grp[i].kernel_type) ? gmx_omp_nthreads_get(emntNonbonded) : 1,
                                nb_alloc, nb_free);
        }
        else
        {
            nbv->grp[i].nbat = nbv->grp[0].nbat;
        }
    }
}

void init_forcerec(FILE *fp,
                   const output_env_t oenv,
                   t_forcerec *fr,
                   t_fcdata   *fcd,
                   const t_inputrec *ir,
                   const gmx_mtop_t *mtop,
                   const t_commrec  *cr,
                   matrix     box,
                   gmx_bool       bMolEpot,
                   const char *tabfn,
                   const char *tabafn,
                   const char *tabpfn,
                   const char *tabbfn,
                   const char *nbpu_opt,
                   int        nbnxn_kernel_preset,
                   gmx_bool   bNoSolvOpt,
                   real       print_force)
{
    int     i,j,m,natoms,ngrp,negp_pp,negptable,egi,egj;
    real    rtab;
    char    *env;
    double  dbl;
    rvec    box_size;
    const t_block *cgs;
    gmx_bool    bGenericKernelOnly;
    gmx_bool    bTab,bSep14tab,bNormalnblists;
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
    
    /* Copy AdResS parameters */
    if (ir->bAdress) {
      fr->adress_type     = ir->adress->type;
      fr->adress_const_wf = ir->adress->const_wf;
      fr->adress_ex_width = ir->adress->ex_width;
      fr->adress_hy_width = ir->adress->hy_width;
      fr->adress_icor     = ir->adress->icor;
      fr->adress_site     = ir->adress->site;
      fr->adress_ex_forcecap = ir->adress->ex_forcecap;
      fr->adress_do_hybridpairs = ir->adress->do_hybridpairs;


      snew(fr->adress_group_explicit , ir->adress->n_energy_grps);
      for (i=0; i< ir->adress->n_energy_grps; i++){
          fr->adress_group_explicit[i]= ir->adress->group_explicit[i];
      }

      fr->n_adress_tf_grps = ir->adress->n_tf_grps;
      snew(fr->adress_tf_table_index, fr->n_adress_tf_grps);
      for (i=0; i< fr->n_adress_tf_grps; i++){
          fr->adress_tf_table_index[i]= ir->adress->tf_table_index[i];
      }
      copy_rvec(ir->adress->refs,fr->adress_refs);
    } else {
      fr->adress_type = eAdressOff;
      fr->adress_do_hybridpairs = FALSE;
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
    fr->efep          = ir->efep;
    fr->sc_alpha      = ir->sc_alpha;
    fr->sc_power      = ir->sc_power;
    fr->sc_sigma6_def = pow(ir->sc_sigma,6);
    fr->sc_sigma6_min = pow(ir->sc_sigma_min,6);
    env = getenv("GMX_SCSIGMA_MIN");
    if (env != NULL)
    {
        dbl = 0;
        sscanf(env,"%lf",&dbl);
        fr->sc_sigma6_min = pow(dbl,6);
        if (fp)
        {
            fprintf(fp,"Setting the minimum soft core sigma to %g nm\n",dbl);
        }
    }

    fr->bNonbonded = TRUE;
    if (getenv("GMX_NO_NONBONDED") != NULL)
    {
        /* turn off non-bonded calculations */
        fr->bNonbonded = FALSE;
        if (fp)
        {
            fprintf(fp, "Found environment varialbe GMX_NO_NONBONDED.\n"
                        "Disabling nonbonded calculations.\n\n");
        }
    }

    bGenericKernelOnly = FALSE;
    if (getenv("GMX_NB_GENERIC") != NULL)
    {
        if (fp != NULL)
        {
            fprintf(fp,
                    "Found environment variable GMX_NB_GENERIC.\n"
                    "Disabling interaction-specific nonbonded kernels.\n\n");
        }
        bGenericKernelOnly = TRUE;
        bNoSolvOpt         = TRUE;
    }

    fr->UseOptimizedKernels = (getenv("GMX_NOOPTIMIZEDKERNELS") == NULL);
    if(fp && fr->UseOptimizedKernels==FALSE)
    {
        fprintf(fp,
                "\nFound environment variable GMX_NOOPTIMIZEDKERNELS.\n"
                "Disabling SSE/SSE2/Altivec/ia64/Power6/Bluegene specific kernels.\n\n");
    }    

#if ( defined(GMX_IA32_SSE2) || defined(GMX_X86_64_SSE2) || defined(GMX_IA32_SSE) || defined(GMX_X86_64_SSE)|| defined(GMX_SSE2) )
    {
        cpuinfo_t cpuinfo;

        get_cpuinfo(&cpuinfo);

        fr->UseOptimizedKernels = cpuinfo.sse2;
    }
#endif
    
    /* Check if we can/should do all-vs-all kernels */
    fr->bAllvsAll       = can_use_allvsall(ir,mtop,FALSE,NULL,NULL);
    fr->AllvsAll_work   = NULL;
    fr->AllvsAll_workgb = NULL;


    /* Neighbour searching stuff */
    fr->cutoff_scheme = ir->cutoff_scheme;
    fr->bGrid         = (ir->ns_type == ensGRID);
    fr->ePBC          = ir->ePBC;

    /* Determine if we will do PBC for distances in bonded interactions */
    if (fr->ePBC == epbcNONE)
    {
        fr->bMolPBC = FALSE;
    }
    else
    {
        if (!DOMAINDECOMP(cr))
        {
            if (fr->cutoff_scheme == ecutsGROUP)
            {
                fr->bMolPBC = ir->bPeriodicMols;
            }
            else
            {
                fr->bMolPBC = TRUE;
                if (getenv("GMX_USE_GRAPH") != NULL)
                {
                    fr->bMolPBC = FALSE;
                    if (fp)
                    {
                        fprintf(fp,"\nGMX_MOLPBC is set, using the graph for bonded interactions\n\n");
                    }
                }
            }
        }
        else
        {
            fr->bMolPBC = dd_bonded_molpbc(cr->dd,fr->ePBC);
        }
    }
    fr->rc_scaling = ir->refcoord_scaling;
    copy_rvec(ir->posres_com,fr->posres_com);
    copy_rvec(ir->posres_comB,fr->posres_comB);
    fr->rlist      = cutoff_inf(ir->rlist);
    fr->rlistlong  = cutoff_inf(ir->rlistlong);
    fr->eeltype    = ir->coulombtype;
    fr->vdwtype    = ir->vdwtype;
    
    fr->bTwinRange = fr->rlistlong > fr->rlist;
    fr->bEwald     = (EEL_PME(fr->eeltype) || fr->eeltype==eelEWALD);
    
    fr->reppow     = mtop->ffparams.reppow;

    if (ir->cutoff_scheme == ecutsGROUP)
    {
        fr->bvdwtab    = (fr->vdwtype != evdwCUT ||
                          !gmx_within_tol(fr->reppow,12.0,10*GMX_DOUBLE_EPS));
        fr->bcoultab   = (!(fr->eeltype == eelCUT || EEL_RF(fr->eeltype)) ||
                          fr->eeltype == eelRF_ZERO);

        if (getenv("GMX_FORCE_TABLES"))
        {
            fr->bvdwtab  = TRUE;
            fr->bcoultab = TRUE;
        }

        if (fp)
        {
            fprintf(fp,"Table routines are used for coulomb: %s\n",bool_names[fr->bcoultab]);
            fprintf(fp,"Table routines are used for vdw:     %s\n",bool_names[fr->bvdwtab ]);
        }
    }
    else
    {
        if (!gmx_within_tol(fr->reppow,12.0,10*GMX_DOUBLE_EPS))
        {
            gmx_fatal(FARGS,"Cut-off scheme %S only supports LJ repulsion power 12",ecutscheme_names[ir->cutoff_scheme]);
        }
        fr->bvdwtab  = FALSE;
        fr->bcoultab = FALSE;
    }
    
    /* Tables are used for direct ewald sum */
    if(fr->bEwald)
    {
        if (EEL_PME(ir->coulombtype))
        {
            if (fp)
                fprintf(fp,"Will do PME sum in reciprocal space.\n");
            if (ir->coulombtype == eelP3M_AD)
            {
                please_cite(fp,"Hockney1988");
                please_cite(fp,"Ballenegger2012");
            }
            else
            {
                please_cite(fp,"Essmann95a");
            }
            
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
        init_ewald_tab(&(fr->ewald_table), cr, ir, fp);
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
    fr->rcoulomb        = cutoff_inf(ir->rcoulomb);
    
    /* Parameters for generalized RF */
    fr->zsquare = 0.0;
    fr->temp    = 0.0;
    
    if (fr->eeltype == eelGRF)
    {
        init_generalized_rf(fp,mtop,ir,fr);
    }
    else if (fr->eeltype == eelSHIFT)
    {
        for(m=0; (m<DIM); m++)
            box_size[m]=box[m][m];
        
        if ((fr->eeltype == eelSHIFT && fr->rcoulomb > fr->rcoulomb_switch))
            set_shift_consts(fp,fr->rcoulomb_switch,fr->rcoulomb,box_size,fr);
    }
    
    fr->bF_NoVirSum = (EEL_FULL(fr->eeltype) ||
                       gmx_mtop_ftype_count(mtop,F_POSRES) > 0 ||
                       IR_ELEC_FIELD(*ir) ||
                       (fr->adress_icor != eAdressICOff)
                      );
    
    /* Mask that says whether or not this NBF list should be computed */
    /*  if (fr->bMask == NULL) {
        ngrp = ir->opts.ngener*ir->opts.ngener;
        snew(fr->bMask,ngrp);*/
    /* Defaults to always */
    /*    for(i=0; (i<ngrp); i++)
          fr->bMask[i] = TRUE;
          }*/
    
    if (fr->cutoff_scheme == ecutsGROUP &&
        ncg_mtop(mtop) > fr->cg_nalloc && !DOMAINDECOMP(cr)) {
        /* Count the total number of charge groups */
        fr->cg_nalloc = ncg_mtop(mtop);
        srenew(fr->cg_cm,fr->cg_nalloc);
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
    fr->rvdw        = cutoff_inf(ir->rvdw);
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
    
    fr->eDispCorr = ir->eDispCorr;
    if (ir->eDispCorr != edispcNO)
    {
        set_avcsixtwelve(fp,fr,mtop);
    }
    
    if (fr->bBHAM)
    {
        set_bham_b_max(fp,fr,mtop);
    }

    fr->bGB = (ir->implicit_solvent == eisGBSA);
	fr->gb_epsilon_solvent = ir->gb_epsilon_solvent;

    /* Copy the GBSA data (radius, volume and surftens for each
     * atomtype) from the topology atomtype section to forcerec.
     */
    snew(fr->atype_radius,fr->ntype);
    snew(fr->atype_vol,fr->ntype);
    snew(fr->atype_surftens,fr->ntype);
    snew(fr->atype_gb_radius,fr->ntype);
    snew(fr->atype_S_hct,fr->ntype);

    if (mtop->atomtypes.nr > 0)
    {
        for(i=0;i<fr->ntype;i++)
            fr->atype_radius[i] =mtop->atomtypes.radius[i];
        for(i=0;i<fr->ntype;i++)
            fr->atype_vol[i] = mtop->atomtypes.vol[i];
        for(i=0;i<fr->ntype;i++)
            fr->atype_surftens[i] = mtop->atomtypes.surftens[i];
        for(i=0;i<fr->ntype;i++)
            fr->atype_gb_radius[i] = mtop->atomtypes.gb_radius[i];
        for(i=0;i<fr->ntype;i++)
            fr->atype_S_hct[i] = mtop->atomtypes.S_hct[i];
    }  
	
	/* Generate the GB table if needed */
	if(fr->bGB)
	{
#ifdef GMX_DOUBLE
		fr->gbtabscale=2000;
#else
		fr->gbtabscale=500;
#endif
		
		fr->gbtabr=100;
		fr->gbtab=make_gb_table(fp,oenv,fr,tabpfn,fr->gbtabscale);

        init_gb(&fr->born,cr,fr,ir,mtop,ir->rgbradii,ir->gb_algorithm);

        /* Copy local gb data (for dd, this is done in dd_partition_system) */
        if (!DOMAINDECOMP(cr))
        {
            make_local_gb(cr,fr->born,ir->gb_algorithm);
        }
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
    
    /* This code automatically gives table length tabext without cut-off's,
     * in that case grompp should already have checked that we do not need
     * normal tables and we only generate tables for 1-4 interactions.
     */
    rtab = ir->rlistlong + ir->tabext;

    if (bTab) {
        /* make tables for ordinary interactions */
        if (bNormalnblists) {
            make_nbf_tables(fp,oenv,fr,rtab,cr,tabfn,NULL,NULL,&fr->nblists[0]);
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
                        make_nbf_tables(fp,oenv,fr,rtab,cr,tabfn,
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
    {
        /* generate extra tables with plain Coulomb for 1-4 interactions only */
        fr->tab14 = make_tables(fp,oenv,fr,MASTER(cr),tabpfn,rtab,
                                GMX_MAKETABLES_14ONLY);
    }

    /* Read AdResS Thermo Force table if needed */
    if(fr->adress_icor == eAdressICThermoForce)
    {
        /* old todo replace */ 
        
        if (ir->adress->n_tf_grps > 0){
            make_adress_tf_tables(fp,oenv,fr,ir,tabfn, mtop, box);

        }else{
            /* load the default table */
            snew(fr->atf_tabs, 1);
            fr->atf_tabs[DEFAULT_TF_TABLE] = make_atf_table(fp,oenv,fr,tabafn, box);
        }
    }
    
    /* Wall stuff */
    fr->nwall = ir->nwall;
    if (ir->nwall && ir->wall_type==ewtTABLE)
    {
        make_wall_tables(fp,oenv,ir,tabfn,&mtop->groups,fr);
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
    fr->cginfo_mb = init_cginfo_mb(fp,mtop,fr,bNoSolvOpt,
                                   &fr->bExcl_IntraCGAll_InterCGNone);
    if (DOMAINDECOMP(cr)) {
        fr->cginfo = NULL;
    } else {
        fr->cginfo = cginfo_expand(mtop->nmolblock,fr->cginfo_mb);
    }
    
    if (!DOMAINDECOMP(cr))
    {
        /* When using particle decomposition, the effect of the second argument,
         * which sets fr->hcg, is corrected later in do_md and init_em.
         */
        forcerec_set_ranges(fr,ncg_mtop(mtop),ncg_mtop(mtop),
                            mtop->natoms,mtop->natoms,mtop->natoms);
    }
    
    fr->print_force = print_force;


    /* coarse load balancing vars */
    fr->t_fnbf=0.;
    fr->t_wait=0.;
    fr->timesteps=0;
    
    /* Initialize neighbor search */
    init_ns(fp,cr,&fr->ns,fr,mtop,box);
    
    if (cr->duty & DUTY_PP)
    {
        gmx_setup_kernels(fp,bGenericKernelOnly);

        if (ir->bAdress)
        {
            gmx_setup_adress_kernels(fp,bGenericKernelOnly);
        }
    }

#ifndef GMX_OPENMP
    fr->nthreads = 1;
#else
    init_forcerec_f_threads(fr,mtop->groups.grps[egcENER].nr);
#endif
    
    snew(fr->excl_load,fr->nthreads+1);

    if (fr->cutoff_scheme == ecutsVERLET)
    {
        if (ir->rcoulomb != ir->rvdw)
        {
            gmx_fatal(FARGS,"With Verlet lists rcoulomb and rvdw should be identical");
        }

        init_nb_verlet(fp, &fr->nbv, ir, fr, cr, nbpu_opt, nbnxn_kernel_preset);

        /* initilize interaction constants; 
           TODO should be moved out during modularizzation.
         */
        init_interaction_const(fp, &fr->ic, fr);
    }

    if (ir->eDispCorr != edispcNO)
    {
        calc_enervirdiff(fp,ir->eDispCorr,fr);
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

void forcerec_set_excl_load(t_forcerec *fr,
                            const gmx_localtop_t *top,const t_commrec *cr)
{
    const int *ind,*a;
    int t,i,j,ntot,n,ntarget;

    if (cr != NULL && PARTDECOMP(cr))
    {
        /* No OpenMP with particle decomposition */
        pd_at_range(cr,
                    &fr->excl_load[0],
                    &fr->excl_load[1]);

        return;
    }

    ind = top->excls.index;
    a   = top->excls.a;

    ntot = 0;
    for(i=0; i<top->excls.nr; i++)
    {
        for(j=ind[i]; j<ind[i+1]; j++)
        {
            if (a[j] > i)
            {
                ntot++;
            }
        }
    }

    fr->excl_load[0] = 0;
    n = 0;
    i = 0;
    for(t=1; t<=fr->nthreads; t++)
    {
        ntarget = (ntot*t)/fr->nthreads;
        while(i < top->excls.nr && n < ntarget)
        {
            for(j=ind[i]; j<ind[i+1]; j++)
            {
                if (a[j] > i)
                {
                    n++;
                }
            }
            i++;
        }
        fr->excl_load[t] = i;
    }
}

