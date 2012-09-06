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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef GMX_THREAD_SHM_FDECOMP
#include <thread_mpi.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include "typedefs.h"
#include "txtdump.h"
#include "smalloc.h"
#include "ns.h"
#include "vec.h"
#include "maths.h"
#include "macros.h"
#include "force.h"
#include "names.h"
#include "main.h"
#include "xvgr.h"
#include "gmx_fatal.h"
#include "physics.h"
#include "force.h"
#include "bondf.h"
#include "nrnb.h"
#include "smalloc.h"
#include "nonbonded.h"

#include "nb_kernel_c/nb_kernel_c.h"
#include "nb_kernel_adress_c/nb_kernel_c_adress.h"
#include "nb_free_energy.h"
#include "nb_generic.h"
#include "nb_generic_cg.h"
#include "nb_generic_adress.h"


/* 1,4 interactions uses kernel 330 directly */
#include "nb_kernel_c/nb_kernel330.h" 
#include "nb_kernel_adress_c/nb_kernel330_adress.h"



#if 0 && defined (GMX_X86_SSE2)
#  ifdef GMX_DOUBLE
#    include "nb_kernel_sse2_double/nb_kernel_sse2_double.h"
#  else
#    include "nb_kernel_sse2_single/nb_kernel_sse2_single.h"
#  endif
#endif

#if defined(GMX_FORTRAN)
#  ifdef GMX_DOUBLE
#    include "nb_kernel_f77_double/nb_kernel_f77_double.h"
#  else
#    include "nb_kernel_f77_single/nb_kernel_f77_single.h"
#  endif
#endif


#ifdef GMX_POWER6
#include "nb_kernel_power6/nb_kernel_power6.h"
#endif

#ifdef GMX_BLUEGENE
#include "nb_kernel_bluegene/nb_kernel_bluegene.h"
#endif



enum { TABLE_NONE, TABLE_COMBINED, TABLE_COUL, TABLE_VDW, TABLE_NR };


/* Table version for each kernel.
 */
static const int 
nb_kernel_table[eNR_NBKERNEL_NR] = 
{
  TABLE_NONE,     /* kernel010 */
  TABLE_NONE,     /* kernel020 */
  TABLE_VDW,      /* kernel030 */
  TABLE_NONE,     /* kernel100 */
  TABLE_NONE,     /* kernel101 */
  TABLE_NONE,     /* kernel102 */
  TABLE_NONE,     /* kernel103 */
  TABLE_NONE,     /* kernel104 */
  TABLE_NONE,     /* kernel110 */
  TABLE_NONE,     /* kernel111 */
  TABLE_NONE,     /* kernel112 */
  TABLE_NONE,     /* kernel113 */
  TABLE_NONE,     /* kernel114 */
  TABLE_NONE,     /* kernel120 */
  TABLE_NONE,     /* kernel121 */
  TABLE_NONE,     /* kernel122 */
  TABLE_NONE,     /* kernel123 */
  TABLE_NONE,     /* kernel124 */
  TABLE_VDW,      /* kernel130 */
  TABLE_VDW,      /* kernel131 */
  TABLE_VDW,      /* kernel132 */
  TABLE_VDW,      /* kernel133 */
  TABLE_VDW,      /* kernel134 */
  TABLE_NONE,     /* kernel200 */
  TABLE_NONE,     /* kernel201 */
  TABLE_NONE,     /* kernel202 */
  TABLE_NONE,     /* kernel203 */
  TABLE_NONE,     /* kernel204 */
  TABLE_NONE,     /* kernel210 */
  TABLE_NONE,     /* kernel211 */
  TABLE_NONE,     /* kernel212 */
  TABLE_NONE,     /* kernel213 */
  TABLE_NONE,     /* kernel214 */
  TABLE_NONE,     /* kernel220 */
  TABLE_NONE,     /* kernel221 */
  TABLE_NONE,     /* kernel222 */
  TABLE_NONE,     /* kernel223 */
  TABLE_NONE,     /* kernel224 */
  TABLE_VDW,      /* kernel230 */
  TABLE_VDW,      /* kernel231 */
  TABLE_VDW,      /* kernel232 */
  TABLE_VDW,      /* kernel233 */
  TABLE_VDW,      /* kernel234 */
  TABLE_COUL,     /* kernel300 */
  TABLE_COUL,     /* kernel301 */
  TABLE_COUL,     /* kernel302 */
  TABLE_COUL,     /* kernel303 */
  TABLE_COUL,     /* kernel304 */
  TABLE_COUL,     /* kernel310 */
  TABLE_COUL,     /* kernel311 */
  TABLE_COUL,     /* kernel312 */
  TABLE_COUL,     /* kernel313 */
  TABLE_COUL,     /* kernel314 */
  TABLE_COUL,     /* kernel320 */
  TABLE_COUL,     /* kernel321 */
  TABLE_COUL,     /* kernel322 */
  TABLE_COUL,     /* kernel323 */
  TABLE_COUL,     /* kernel324 */
  TABLE_COMBINED, /* kernel330 */
  TABLE_COMBINED, /* kernel331 */
  TABLE_COMBINED, /* kernel332 */
  TABLE_COMBINED, /* kernel333 */
  TABLE_COMBINED, /* kernel334 */
  TABLE_NONE,     /* kernel400 */
  TABLE_NONE,     /* kernel410 */
  TABLE_VDW       /* kernel430 */
};



static nb_kernel_t **
nb_kernel_list = NULL;

static nb_adress_kernel_t **
nb_kernel_list_adress = NULL;

void
gmx_setup_kernels(FILE *fplog,t_forcerec *fr,gmx_bool bGenericKernelOnly)
{
    int i;
        
    snew(nb_kernel_list,eNR_NBKERNEL_NR);
    
    /* Note that later calls overwrite earlier, so the preferred (fastest)
     * version should be at the end. For instance, we call SSE after 3DNow.
     */
    
    for(i=0; i<eNR_NBKERNEL_NR; i++)
    {
        nb_kernel_list[i] = NULL;
    }
    
    if (bGenericKernelOnly)
    {
        return;
    }
	
	if(fplog)
    {
	    fprintf(fplog,"Configuring nonbonded kernels...\n");
    }
	
    nb_kernel_setup(fplog,nb_kernel_list);
    
    if(fr->use_cpu_acceleration==FALSE)
    {
        return;
    }

    /* Setup kernels. The last called setup routine will overwrite earlier assignments,
	 * so we should e.g. test SSE3 support _after_ SSE2 support,
     * and call e.g. Fortran setup before SSE.
	 */
    
#if defined(GMX_FORTRAN) && defined(GMX_DOUBLE)   
    nb_kernel_setup_f77_double(fplog,nb_kernel_list);
#endif
	
#if defined(GMX_FORTRAN) && !defined(GMX_DOUBLE)   
    nb_kernel_setup_f77_single(fplog,nb_kernel_list);
#endif
	
#ifdef GMX_BLUEGENE
    nb_kernel_setup_bluegene(fplog,nb_kernel_list);
#endif
	
#ifdef GMX_POWER6
    nb_kernel_setup_power6(fplog,nb_kernel_list);
#endif
	
	if(fplog)
    {
	    fprintf(fplog,"\n\n");
    }
}

void
gmx_setup_adress_kernels(FILE *fplog,gmx_bool bGenericKernelOnly) {
    int i;

    snew(nb_kernel_list_adress, eNR_NBKERNEL_NR);

    for (i = 0; i < eNR_NBKERNEL_NR; i++) {
        nb_kernel_list_adress[i] = NULL;
    }

    if (bGenericKernelOnly)
    {
        return;
    }

    nb_kernel_setup_adress(fplog, nb_kernel_list_adress);
}

void do_nonbonded(t_commrec *cr,t_forcerec *fr,
                  rvec x[],rvec f[],t_mdatoms *mdatoms,t_blocka *excl,
                  real egnb[],real egcoul[],real egpol[],rvec box_size,
                  t_nrnb *nrnb,real *lambda, real *dvdl,
                  int nls,int eNL,int flags)
{
    gmx_bool            bLR,bDoForces,bForeignLambda;
	t_nblist *      nlist;
	real *          fshift;
	int             n,n0,n1,i,i0,i1,nrnb_ind,sz;
	t_nblists       *nblists;
	gmx_bool        bWater;
	FILE *          fp;
	int             fac=0;
	int             nthreads = 1;
	int             tabletype;
	int             outeriter,inneriter;
	real *          tabledata = NULL;
	gmx_gbdata_t    gbdata;

        nb_kernel_t         *kernelptr=NULL;
        nb_adress_kernel_t  *adresskernelptr=NULL;
    
        gmx_bool        bCG; /* for AdresS */
        int             k;/* for AdresS */

    bLR            = (flags & GMX_DONB_LR);
    bDoForces      = (flags & GMX_DONB_FORCES);
    bForeignLambda = (flags & GMX_DONB_FOREIGNLAMBDA); 

    bCG = FALSE;  /* for AdresS */

	gbdata.gb_epsilon_solvent = fr->gb_epsilon_solvent;
	gbdata.epsilon_r = fr->epsilon_r;
	gbdata.gpol               = egpol;
    
    if (!fr->adress_type==eAdressOff && !bDoForces){
        gmx_fatal(FARGS,"No force kernels not implemeted for adress");
    }

    if(fr->bAllvsAll) 
    {
        if(fr->bGB)
        {
#if 0 && defined (GMX_X86_SSE2)
# ifdef GMX_DOUBLE
            if(fr->use_cpu_acceleration)
            {
                nb_kernel_allvsallgb_sse2_double(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                                 &outeriter,&inneriter,&fr->AllvsAll_work);
            }
            else
            {
                nb_kernel_allvsallgb(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                     &outeriter,&inneriter,&fr->AllvsAll_work);        
            }
#  else /* not double */
            if(fr->use_cpu_acceleration)
            {
                nb_kernel_allvsallgb_sse2_single(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                                 &outeriter,&inneriter,&fr->AllvsAll_work);
            }
            else
            {
                nb_kernel_allvsallgb(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                     &outeriter,&inneriter,&fr->AllvsAll_work);        
            }
#  endif /* double/single alt. */
#else /* no SSE support compiled in */
            nb_kernel_allvsallgb(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,egpol,
                                 &outeriter,&inneriter,&fr->AllvsAll_work);                    
#endif
            inc_nrnb(nrnb,eNR_NBKERNEL_ALLVSALLGB,inneriter);
        }
        else
        { 
#if 0 && defined (GMX_X86_SSE2)
# ifdef GMX_DOUBLE
            if(fr->use_cpu_acceleration)
            {
                nb_kernel_allvsall_sse2_double(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                                               &outeriter,&inneriter,&fr->AllvsAll_work);
            }
            else 
            {
                nb_kernel_allvsall(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                                   &outeriter,&inneriter,&fr->AllvsAll_work);            
            }
            
#  else /* not double */
            if(fr->use_cpu_acceleration)
            {
                nb_kernel_allvsall_sse2_single(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                                               &outeriter,&inneriter,&fr->AllvsAll_work);
            }
            else 
            {
                nb_kernel_allvsall(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                                   &outeriter,&inneriter,&fr->AllvsAll_work);            
            }

#  endif /* double/single check */
#else /* No SSE2 support compiled in */
            nb_kernel_allvsall(fr,mdatoms,excl,x[0],f[0],egcoul,egnb,
                               &outeriter,&inneriter,&fr->AllvsAll_work);
#endif            
            
            inc_nrnb(nrnb,eNR_NBKERNEL_ALLVSALL,inneriter);
        }
        inc_nrnb(nrnb,eNR_NBKERNEL_OUTER,outeriter);
        return;
    }
	
    if (eNL >= 0) 
    {
		i0 = eNL;
		i1 = i0+1;
    }
    else
    {
		i0 = 0;
		i1 = eNL_NR;
	}
	
	if (nls >= 0) 
	{
		n0 = nls;
		n1 = nls+1;
	}
	else 
	{
		n0 = 0;
		n1 = fr->nnblists;
	}
	
	if(nb_kernel_list == NULL)
    {
		gmx_fatal(FARGS,"gmx_setup_kernels has not been called");
    }
  
    fshift = fr->fshift[0];
  
	for(n=n0; (n<n1); n++) 
	{
		nblists = &fr->nblists[n];
		for(i=i0; (i<i1); i++) 
		{
			outeriter = inneriter = 0;
      
			if (bLR) 
			{
				nlist = &(nblists->nlist_lr[i]);
			}
			else
			{
				nlist = &(nblists->nlist_sr[i]);
			}
			
			if (nlist->nri > 0) 
			{
				nrnb_ind = nlist->il_code;
				
				if(nrnb_ind==eNR_NBKERNEL_FREE_ENERGY)
				{
					/* generic free energy, use combined table */
					tabledata = nblists->tab.tab;
				}
				else
				{
                    if (bForeignLambda)
                    {
                        /* We don't need the non-perturbed interactions */
                        continue;
                    }

					tabletype = nb_kernel_table[nrnb_ind];
					
					/* normal kernels, not free energy */
					if (!bDoForces)
					{
						nrnb_ind += eNR_NBKERNEL_NR/2;
					}
					
					if(tabletype == TABLE_COMBINED)
					{
						tabledata = nblists->tab.tab;
					}
					else if(tabletype == TABLE_COUL)
					{
						tabledata = nblists->coultab;
					}
					else if(tabletype == TABLE_VDW)
					{
						tabledata = nblists->vdwtab;
					}
					else
					{
						tabledata = NULL;
					}
				}
				
				nlist->count = 0;
				
				if(nlist->free_energy)
				{
					if(nlist->ivdw==enbvdwBHAM)
					{
						gmx_fatal(FARGS,"Cannot do free energy Buckingham interactions.");
					}
					
					gmx_nb_free_energy_kernel(nlist->icoul,
											  nlist->ivdw,
											  nlist->nri,
											  nlist->iinr,
											  nlist->jindex,
											  nlist->jjnr,
											  nlist->shift,
											  fr->shift_vec[0],
											  fshift,
											  nlist->gid,
											  x[0],
											  f[0],
											  mdatoms->chargeA,
											  mdatoms->chargeB,
											  fr->epsfac,
											  fr->k_rf,
											  fr->c_rf,
											  fr->ewaldcoeff,
											  egcoul,
											  mdatoms->typeA,
											  mdatoms->typeB,
											  fr->ntype,
											  fr->nbfp,
											  egnb,
											  nblists->tab.scale,
											  tabledata,
											  lambda[efptCOUL],
                                              lambda[efptVDW],
                                              dvdl,
                                              fr->sc_alphacoul,
											  fr->sc_alphavdw,
											  fr->sc_power,
											  fr->sc_r_power,
											  fr->sc_sigma6_def,
                                              fr->sc_sigma6_min,
                                              bDoForces,
											  &outeriter,
											  &inneriter);
                }
                else if (nlist->enlist == enlistCG_CG)
                {
		    if (fr->adress_type==eAdressOff){
                    /* Call the charge group based inner loop */
                       gmx_nb_generic_cg_kernel(nlist,
                                                fr,
                                                mdatoms,
                                                x[0],
                                                f[0],
                                                fshift,
                                                egcoul,
                                                egnb,
                                                nblists->tab.scale,
                                                tabledata,
                                                &outeriter,
                                                &inneriter);
		    }
		    else
		    {
                       /*gmx_nb_generic_adress_kernel(nlist,
                                                fr,
                                                mdatoms,
                                                x[0],
                                                f[0],
                                                fshift,
                                                egcoul,
                                                egnb,
                                                nblists->tab.scale,
                                                tabledata,
                                                &outeriter,
                                                &inneriter);*/
                          gmx_fatal(FARGS,"Death & horror! Adress cgcg kernel not implemented anymore.\n");

		    }
                }
                else
                {
                    /* AdresS*/
                    /* for adress we need to determine for each energy group wether it is explicit or coarse-grained */
                    if (!fr->adress_type == eAdressOff) {                        
                        bCG = FALSE;
                        if ( !fr->adress_group_explicit[ mdatoms->cENER[nlist->iinr[0]] ] ){
                            bCG=TRUE;
                        }
                        /* If this processor has only explicit atoms (w=1)
                          skip the coarse grained force calculation. Same for
                         only coarsegrained atoms and explicit interactions.
                         Last condition is to make sure that generic kernel is not
                         skipped*/
                        if (mdatoms->pureex && bCG && nb_kernel_list[nrnb_ind] != NULL) continue;
                        if (mdatoms->purecg && !bCG && nb_kernel_list[nrnb_ind] != NULL) continue;
                    }

                    if (fr->adress_type == eAdressOff ||
                            mdatoms->pureex ||
                            mdatoms->purecg){
                        /* if we only have to calculate pure cg/ex interactions
                         we can use the faster standard gromacs kernels*/
                        kernelptr = nb_kernel_list[nrnb_ind];
                    }else{
                        /* This processor has hybrid interactions which means
                         * we have to
                         * use our own kernels. We have two kernel types: one that
                         * calculates the forces with the explicit prefactor w1*w2
                         * and one for coarse-grained with (1-w1*w2)
                         * explicit kernels are the second part of the kernel
                         *  list */
                        if (!bCG) nrnb_ind += eNR_NBKERNEL_NR/2;                      
                        adresskernelptr = nb_kernel_list_adress[nrnb_ind];
                    }
                    
                    if (kernelptr == NULL && adresskernelptr == NULL)
                     {
                        /* Call a generic nonbonded kernel */
                        
                        /* If you want to hack/test your own interactions,
                         * do it in this routine and make sure it is called
                         * by setting the environment variable GMX_NB_GENERIC.
                         */
                        if (fr->adress_type==eAdressOff){

                        gmx_nb_generic_kernel(nlist,
                                              fr,
                                              mdatoms,
                                              x[0],
                                              f[0],
                                              fshift,
                                              egcoul,
                                              egnb,
                                              nblists->tab.scale,
                                              tabledata,
                                              &outeriter,
                                              &inneriter);
                        }else /* do generic AdResS kernels (slow)*/
                        {

                            gmx_nb_generic_adress_kernel(nlist,
                                                fr,
                                                mdatoms,
                                                x[0],
                                                f[0],
                                                fshift,
                                                egcoul,
                                                egnb,
                                                nblists->tab.scale,
                                                tabledata,
                                                &outeriter,
                                                &inneriter,
                                                bCG);
                        }


                    }
                    else
                    {
                        /* Call nonbonded kernel from function pointer */
                        if (kernelptr!=NULL){
                        (*kernelptr)( &(nlist->nri),
                                      nlist->iinr,
                                      nlist->jindex,
                                      nlist->jjnr,
                                      nlist->shift,
                                      fr->shift_vec[0],
                                      fshift,
                                      nlist->gid,
                                      x[0],
                                      f[0],
                                      mdatoms->chargeA,
                                      &(fr->epsfac),
                                      &(fr->k_rf),
                                      &(fr->c_rf),
                                      egcoul,
                                      mdatoms->typeA,
                                      &(fr->ntype),
                                      fr->nbfp,
                                      egnb,
                                      &(nblists->tab.scale),
                                      tabledata,
                                      fr->invsqrta,
                                      fr->dvda,
                                      &(fr->gbtabscale),
                                      fr->gbtab.tab,
                                      &nthreads,
                                      &(nlist->count),
                                      nlist->mtx,
                                      &outeriter,
                                      &inneriter,
                                      (real *)&gbdata);
                        }else if (adresskernelptr != NULL)
                        { /* Adress kernels */
                          (*adresskernelptr)( &(nlist->nri),
                                      nlist->iinr,
                                      nlist->jindex,
                                      nlist->jjnr,
                                      nlist->shift,
                                      fr->shift_vec[0],
                                      fshift,
                                      nlist->gid,
                                      x[0],
                                      f[0],
                                      mdatoms->chargeA,
                                      &(fr->epsfac),
                                      &(fr->k_rf),
                                      &(fr->c_rf),
                                      egcoul,
                                      mdatoms->typeA,
                                      &(fr->ntype),
                                      fr->nbfp,
                                      egnb,
                                      &(nblists->tab.scale),
                                      tabledata,
                                      fr->invsqrta,
                                      fr->dvda,
                                      &(fr->gbtabscale),
                                      fr->gbtab.tab,
                                      &nthreads,
                                      &(nlist->count),
                                      nlist->mtx,
                                      &outeriter,
                                      &inneriter,
                                      fr->adress_ex_forcecap,
                                      mdatoms->wf);
                        }
                    }
                }
                
                /* Update flop accounting */
				
				/* Outer loop in kernel */
                switch (nlist->enlist) {
                case enlistATOM_ATOM:   fac =  1; break;
                case enlistSPC_ATOM:    fac =  3; break;
                case enlistSPC_SPC:     fac =  9; break;
                case enlistTIP4P_ATOM:  fac =  4; break;
                case enlistTIP4P_TIP4P: fac = 16; break;
                case enlistCG_CG:       fac =  1; break;
                }
                inc_nrnb(nrnb,eNR_NBKERNEL_OUTER,fac*outeriter);

                /* inner loop in kernel */
                inc_nrnb(nrnb,nrnb_ind,inneriter);
            }
        }
    }
}


real 
do_listed_vdw_q(int ftype,int nbonds,
                const t_iatom iatoms[],const t_iparams iparams[],
                const rvec x[],rvec f[],rvec fshift[],
                const t_pbc *pbc,const t_graph *g,
                real *lambda, real *dvdl,
                const t_mdatoms *md,
                const t_forcerec *fr,gmx_grppairener_t *grppener,
                int *global_atom_index)
{
    static    gmx_bool bWarn=FALSE;
    real      eps,r2,*tab,rtab2=0;
    rvec      dx,x14[2],f14[2];
    int       i,ai,aj,itype;
    int       typeA[2]={0,0},typeB[2]={0,1};
    real      chargeA[2]={0,0},chargeB[2];
    int       gid,shift_vir,shift_f;
    int       j_index[] = { 0, 1 };
    int       i0=0,i1=1,i2=2;
    ivec      dt;
    int       outeriter,inneriter;
    int       nthreads = 1;
    int       count;
    real      krf,crf,tabscale;
    int       ntype=0;
    real      *nbfp=NULL;
    real      *egnb=NULL,*egcoul=NULL;
    t_nblist  tmplist;
    int       icoul,ivdw;
    gmx_bool      bMolPBC,bFreeEnergy;
    
    gmx_bool      bCG; /* AdResS*/
    real      wf14[2]={0,0}; /* AdResS*/
   
#if GMX_THREAD_SHM_FDECOMP
    pthread_mutex_t mtx;
#else
    void *    mtx = NULL;
#endif

    
#if GMX_THREAD_SHM_FDECOMP
    pthread_mutex_initialize(&mtx);
#endif

    bMolPBC = fr->bMolPBC;

    switch (ftype) {
    case F_LJ14:
    case F_LJC14_Q:
        eps = fr->epsfac*fr->fudgeQQ;
        ntype  = 1;
        egnb   = grppener->ener[egLJ14];
        egcoul = grppener->ener[egCOUL14];
        break;
    case F_LJC_PAIRS_NB:
        eps = fr->epsfac;
        ntype  = 1;
        egnb   = grppener->ener[egLJSR];
        egcoul = grppener->ener[egCOULSR];
        break;
    default:
        gmx_fatal(FARGS,"Unknown function type %d in do_nonbonded14",
                  ftype);
    }
    tab = fr->tab14.tab;
    rtab2 = sqr(fr->tab14.r);
    tabscale = fr->tab14.scale;

    krf = fr->k_rf;
    crf = fr->c_rf;

    /* Determine the values for icoul/ivdw. */
    if (fr->bEwald) {
        icoul = enbcoulOOR;
    } 
    else if(fr->bcoultab)
    {
        icoul = enbcoulTAB;
    }
    else if(fr->eeltype == eelRF_NEC)
    {
        icoul = enbcoulRF;
    }
    else 
    {
        icoul = enbcoulOOR;
    }
    
    if(fr->bvdwtab)
    {
        ivdw = enbvdwTAB;
    }
    else
    {
        ivdw = enbvdwLJ;
    }
    
    
    bCG = FALSE; /*Adres*/
    /* We don't do SSE or altivec here, due to large overhead for 4-fold 
     * unrolling on short lists 
     */
    
    bFreeEnergy = FALSE;
    for(i=0; (i<nbonds); ) 
    {
        itype = iatoms[i++];
        ai    = iatoms[i++];
        aj    = iatoms[i++];
        gid   = GID(md->cENER[ai],md->cENER[aj],md->nenergrp);
        
        if (!fr->adress_type == eAdressOff) {
            if (fr->adress_group_explicit[md->cENER[ai]] != fr->adress_group_explicit[md->cENER[aj]]){
                /*exclude cg-ex interaction*/
                continue;
            }           
            bCG = !fr->adress_group_explicit[md->cENER[ai]];
            wf14[0] = md->wf[ai];
            wf14[1] = md->wf[aj];
        }
        switch (ftype) {
        case F_LJ14:
            bFreeEnergy =
                (fr->efep != efepNO &&
                 ((md->nPerturbed && (md->bPerturbed[ai] || md->bPerturbed[aj])) ||
                  iparams[itype].lj14.c6A != iparams[itype].lj14.c6B ||
                  iparams[itype].lj14.c12A != iparams[itype].lj14.c12B));
            chargeA[0] = md->chargeA[ai];
            chargeA[1] = md->chargeA[aj];
            nbfp = (real *)&(iparams[itype].lj14.c6A);
            break;
        case F_LJC14_Q:
            eps = fr->epsfac*iparams[itype].ljc14.fqq;
            chargeA[0] = iparams[itype].ljc14.qi;
            chargeA[1] = iparams[itype].ljc14.qj;
            nbfp = (real *)&(iparams[itype].ljc14.c6);
            break;
        case F_LJC_PAIRS_NB:
            chargeA[0] = iparams[itype].ljcnb.qi;
            chargeA[1] = iparams[itype].ljcnb.qj;
            nbfp = (real *)&(iparams[itype].ljcnb.c6);
            break;
        }
        
        if (!bMolPBC) 
        {
            /* This is a bonded interaction, atoms are in the same box */
            shift_f = CENTRAL;
            r2 = distance2(x[ai],x[aj]);
        }
        else 
        {
            /* Apply full periodic boundary conditions */
            shift_f = pbc_dx_aiuc(pbc,x[ai],x[aj],dx);
            r2 = norm2(dx);
        }

        if (r2 >= rtab2) 
        {
            if (!bWarn) 
            {
                fprintf(stderr,"Warning: 1-4 interaction between %d and %d "
                        "at distance %.3f which is larger than the 1-4 table size %.3f nm\n", 
			glatnr(global_atom_index,ai),
			glatnr(global_atom_index,aj),
			sqrt(r2), sqrt(rtab2));
                fprintf(stderr,"These are ignored for the rest of the simulation\n");
                fprintf(stderr,"This usually means your system is exploding,\n"
                        "if not, you should increase table-extension in your mdp file\n"
                        "or with user tables increase the table size\n");
                bWarn = TRUE;
            }
            if (debug) 
	      fprintf(debug,"%8f %8f %8f\n%8f %8f %8f\n1-4 (%d,%d) interaction not within cut-off! r=%g. Ignored\n",
		      x[ai][XX],x[ai][YY],x[ai][ZZ],
		      x[aj][XX],x[aj][YY],x[aj][ZZ],
		      glatnr(global_atom_index,ai),
		      glatnr(global_atom_index,aj),
		      sqrt(r2));
        }
        else 
        {
            copy_rvec(x[ai],x14[0]);
            copy_rvec(x[aj],x14[1]);
            clear_rvec(f14[0]);
            clear_rvec(f14[1]);
#ifdef DEBUG
            fprintf(debug,"LJ14: grp-i=%2d, grp-j=%2d, ngrp=%2d, GID=%d\n",
                    md->cENER[ai],md->cENER[aj],md->nenergrp,gid);
#endif
            
	    outeriter = inneriter = count = 0;
	    if (bFreeEnergy)
        {
            chargeB[0] = md->chargeB[ai];
            chargeB[1] = md->chargeB[aj];
            /* We pass &(iparams[itype].lj14.c6A) as LJ parameter matrix
             * to the innerloops.
             * Here we use that the LJ-14 parameters are stored in iparams
             * as c6A,c12A,c6B,c12B, which are referenced correctly
             * in the innerloops if we assign type combinations 0-0 and 0-1
             * to atom pair ai-aj in topologies A and B respectively.
             */

            /* need to do a bit of a kludge here -- the way it is set up,
               if the charges change, but the vdw do not, then even though bFreeEnergy is on,
               it won't work, because all the bonds are perturbed.
            */
            if(ivdw==enbvdwBHAM)
            {
                gmx_fatal(FARGS,"Cannot do free energy Buckingham interactions.");
            }
            count = 0;
            gmx_nb_free_energy_kernel(icoul,
                                      ivdw,
                                      i1,
                                      &i0,
                                      j_index,
                                      &i1,
                                      &shift_f,
                                      fr->shift_vec[0],
                                      fshift[0],
                                      &gid,
                                      x14[0],
                                      f14[0],
                                      chargeA,
                                      chargeB,
                                      eps,
                                      krf,
                                      crf,
                                      fr->ewaldcoeff,
                                      egcoul,
                                      typeA,
                                      typeB,
                                      ntype,
                                      nbfp,
                                      egnb,
                                      tabscale,
                                      tab,
                                      lambda[efptCOUL],
                                      lambda[efptVDW],
                                      dvdl,
                                      fr->sc_alphacoul,
                                      fr->sc_alphavdw,
                                      fr->sc_power,
                                      6.0, /* for 1-4's use the 6 power always - 48 power too high because of where they are forced to be */
                                      fr->sc_sigma6_def,
                                      fr->sc_sigma6_min,
                                      TRUE,
                                      &outeriter,
                                      &inneriter);
        }
        else 
        { 
          if (fr->adress_type==eAdressOff || !fr->adress_do_hybridpairs){
            /* Not perturbed - call kernel 330 */
            nb_kernel330
                ( &i1,
                  &i0,
                  j_index,
                  &i1,
                  &shift_f,
                  fr->shift_vec[0],
                  fshift[0],
                  &gid,
                  x14[0],
                  f14[0],
                  chargeA,
                  &eps,
                  &krf,
                  &crf,
                  egcoul,
                  typeA,
                  &ntype,
                  nbfp,
                  egnb,
                  &tabscale,
                  tab,
                  NULL,
                  NULL,
                  NULL,
                  NULL,
                  &nthreads,
                  &count,
                  (void *)&mtx,
                  &outeriter,
                  &inneriter,
                  NULL);                
                } else {
                    if (bCG) {
                        nb_kernel330_adress_cg(&i1,
                                &i0,
                                j_index,
                                &i1,
                                &shift_f,
                                fr->shift_vec[0],
                                fshift[0],
                                &gid,
                                x14[0],
                                f14[0],
                                chargeA,
                                &eps,
                                &krf,
                                &crf,
                                egcoul,
                                typeA,
                                &ntype,
                                nbfp,
                                egnb,
                                &tabscale,
                                tab,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                &nthreads,
                                &count,
                                (void *) &mtx,
                                &outeriter,
                                &inneriter,
                                fr->adress_ex_forcecap,
                                wf14);
                    } else {
                        nb_kernel330_adress_ex(&i1,
                                &i0,
                                j_index,
                                &i1,
                                &shift_f,
                                fr->shift_vec[0],
                                fshift[0],
                                &gid,
                                x14[0],
                                f14[0],
                                chargeA,
                                &eps,
                                &krf,
                                &crf,
                                egcoul,
                                typeA,
                                &ntype,
                                nbfp,
                                egnb,
                                &tabscale,
                                tab,
                                NULL,
                                NULL,
                                NULL,
                                NULL,
                                &nthreads,
                                &count,
                                (void *) &mtx,
                                &outeriter,
                                &inneriter,
                                fr->adress_ex_forcecap,
                                wf14);
                    }

                }
            }
        
        /* Add the forces */
        rvec_inc(f[ai],f14[0]);
        rvec_dec(f[aj],f14[0]);
        
        if (g) 
        {
            /* Correct the shift forces using the graph */
            ivec_sub(SHIFT_IVEC(g,ai),SHIFT_IVEC(g,aj),dt);    
            shift_vir = IVEC2IS(dt);
            rvec_inc(fshift[shift_vir],f14[0]);
            rvec_dec(fshift[CENTRAL],f14[0]);
        }
        
	    /* flops: eNR_KERNEL_OUTER + eNR_KERNEL330 + 12 */
        }
    }
    return 0.0;
}


