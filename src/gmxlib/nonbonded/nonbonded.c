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

#ifdef GMX_THREAD_MPI
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

#include "nb_kernel.h"
#include "nb_kernel_c/nb_kernel_c.h"

#include "nb_free_energy.h"
#include "nb_generic.h"
#include "nb_generic_cg.h"
#include "nb_generic_adress.h"


/* 1,4 interactions uses kernel 330 directly */
#include "nb_kernel_c/nb_kernel330.h" 


#ifdef GMX_THREAD_MPI
static tMPI_Thread_mutex_t nonbonded_setup_mutex = TMPI_THREAD_MUTEX_INITIALIZER;
#endif
static gmx_bool            nonbonded_setup_done  = FALSE;


void
gmx_nonbonded_setup(FILE *         fplog,
                    t_forcerec *   fr,
                    gmx_bool       bGenericKernelOnly)
{
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_lock(&nonbonded_setup_mutex);
#endif
    /* Here we are guaranteed only one thread made it. */
    if(nonbonded_setup_done==FALSE)
    {
        /* Add the generic kernels to the structure stored statically in nb_kernel.c */
        
        /* Add accelerated C kernels */
        nb_kernel_list_add_kernels(kernellist_c,kernellist_c_size);
        
        /* Add accelerated SSE kernels */

        /* Create a hash for faster lookups */
        nb_kernel_list_hash_init();
    
        nonbonded_setup_done=TRUE;
    }
#ifdef GMX_THREAD_MPI
    tMPI_Thread_mutex_unlock(&nonbonded_setup_mutex);
#endif
}



void
gmx_nonbonded_set_kernel_pointers(FILE *log, t_nblist *nl)
{
    const char *     elec;
    const char *     vdw;
    const char *     geom;
    const char *     other;
    const char *     vf;
    
    const char *     archs[] = { "C" };
    int              i;

    
    if(nonbonded_setup_done==FALSE)
    {
        /* This should have been called before, but just in case... */
        gmx_nonbonded_setup(NULL,NULL,FALSE);
    }
    
    /* Not used yet */
    other="";

    if(nl->free_energy==1)
    {
        nl->kernelptr_vf = gmx_nb_free_energy_kernel;
        nl->kernelptr_f  = gmx_nb_free_energy_kernel;
    }
    else
    {
        elec = gmx_nbkernel_elec_names[nl->ielec];
        vdw  = gmx_nbkernel_vdw_names[nl->ivdw];
        geom = gmx_nblist_geometry_names[nl->igeometry];
        
        nl->kernelptr_vf=NULL;
        nl->kernelptr_f=NULL;
        
        for(i=0;i<asize(archs) && nl->kernelptr_vf==NULL ;i++)
        {
                nl->kernelptr_vf = nb_kernel_list_findkernel(log,archs[i],elec,vdw,geom,other,"VF");
        }
        for(i=0;i<asize(archs) && nl->kernelptr_f==NULL ;i++)
        {
            nl->kernelptr_f  = nb_kernel_list_findkernel(log,archs[i],elec,vdw,geom,other,"F");
            /* If there is not force-only optimized kernel, is there a potential & force one? */
            if(nl->kernelptr_f == NULL)
            {
                nl->kernelptr_f  = nb_kernel_list_findkernel(NULL,archs[i],elec,vdw,geom,other,"VF");
            }
        }
    }
    
    if(nl->kernelptr_vf==NULL)
    {
        nl->kernelptr_vf = gmx_nb_generic_kernel;
        nl->kernelptr_f  = gmx_nb_generic_kernel;
    }
    
    return;
}

void do_nonbonded(t_commrec *cr,t_forcerec *fr,
                  rvec x[],rvec f[],t_mdatoms *mdatoms,t_blocka *excl,
                  real egnb[],real egcoul[],real egpol[],rvec box_size,
                  t_nrnb *nrnb,real *lambda, real *dvdl,
                  int nls,int eNL,int flags)
{
    gmx_bool        bLR,bDoForces,bForeignLambda;
	t_nblist *      nlist;
	real *          fshift;
	int             n,n0,n1,i,i0,i1,sz;
	t_nblists       *nblists;
	gmx_bool        bWater;
	FILE *          fp;
	int             fac=0;
	int             nthreads = 1;
	int             tabletype;
    int             count=0; /* dummy, not used */
	int             outeriter,inneriter;
	real *          tabledata = NULL;
	gmx_gbdata_t    gbdata;
    
    nb_kernel_t         *kernelptr=NULL;
        
    bLR            = (flags & GMX_DONB_LR);
    bDoForces      = (flags & GMX_DONB_FORCES);
    bForeignLambda = (flags & GMX_DONB_FOREIGNLAMBDA);
        
	gbdata.gb_epsilon_solvent = fr->gb_epsilon_solvent;
	gbdata.epsilon_r = fr->epsilon_r;
	gbdata.gpol               = egpol;
        
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
				kernelptr = nlist->kernelptr_vf;
				
				if(nlist->free_energy)
				{
					/* generic free energy, use combined table */
					tabledata = nblists->table_elec_vdw.data;
				}
				else
				{
					/* normal kernels, not free energy */

                    if (bForeignLambda)
                    {
                        /* We don't need the non-perturbed interactions */
                        continue;
                    }
					
					if( (nlist->ielec==GMX_NBKERNEL_ELEC_CUBICSPLINETABLE || nlist->ielec==GMX_NBKERNEL_ELEC_EWALD) && nlist->ivdw==GMX_NBKERNEL_VDW_CUBICSPLINETABLE)
					{
						tabledata = nblists->table_elec_vdw.data;
					}
					else if(nlist->ielec==GMX_NBKERNEL_ELEC_CUBICSPLINETABLE  || nlist->ielec==GMX_NBKERNEL_ELEC_EWALD)
					{
						tabledata = nblists->table_elec.data;
					}
					else if(nlist->ivdw==GMX_NBKERNEL_VDW_CUBICSPLINETABLE)
					{
						tabledata = nblists->table_vdw.data;
					}
					else
					{
						tabledata = NULL;
					}
				}
				
				count = 0;
				
				if(nlist->free_energy)
				{
					if(nlist->ivdw==enbvdwBHAM)
					{
						gmx_fatal(FARGS,"Cannot do free energy Buckingham interactions.");
					}
					
					gmx_nb_free_energy_kernel(nlist->ielec,
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
											  nblists->table_elec_vdw.scale,
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
                else if (nlist->igeometry == GMX_NBLIST_GEOMETRY_CG_CG)
                {
                    /* Call the charge group based inner loop */
                    gmx_nb_generic_cg_kernel(nlist,
                                             fr,
                                             mdatoms,
                                             x[0],
                                             f[0],
                                             fshift,
                                             egcoul,
                                             egnb,
                                             nblists->table_elec_vdw.scale,
                                             tabledata,
                                             &outeriter,
                                             &inneriter);
                }
                else
                {
                    if (kernelptr == NULL)
                    {
                        /* Call a generic nonbonded kernel */
                        
                        /* If you want to hack/test your own interactions,
                         * do it in this routine and make sure it is called
                         * by setting the environment variable GMX_NB_GENERIC.
                         */
                        gmx_nb_generic_kernel(nlist,
                                              fr,
                                              mdatoms,
                                              x[0],
                                              f[0],
                                              fshift,
                                              egcoul,
                                              egnb,
                                              nblists->table_elec_vdw.scale,
                                              tabledata,
                                              &outeriter,
                                              &inneriter);
                    }
                    else
                    {
                        /* Call nonbonded kernel from function pointer */
                        if (kernelptr!=NULL)
                        {
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
                                         &(nblists->table_elec_vdw.scale),
                                         tabledata,
                                         fr->invsqrta,
                                         fr->dvda,
                                         &(fr->gbtabscale),
                                         fr->gbtab.data,
                                         &nthreads,
                                         &count,
                                         NULL,
                                         &outeriter,
                                         &inneriter,
                                         (real *)&gbdata);
                        }
                    }
                }
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
    int       ielec,ivdw;
    gmx_bool      bMolPBC,bFreeEnergy;
       
    void *    mtx = NULL;
    
#if GMX_THREAD_SHM_FDECOMP
    pthread_mutex_initialize(&mtx);
#endif
    
    bMolPBC = fr->bMolPBC;
    
    switch (ftype)
    {
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
    tab = fr->tab14.data;
    rtab2 = sqr(fr->tab14.r);
    tabscale = fr->tab14.scale;
    
    krf = fr->k_rf;
    crf = fr->c_rf;
    
    /* Determine the values for ielec/ivdw. */
    if (fr->bEwald) {
        ielec = enbcoulOOR;
    }
    else if(fr->bcoultab)
    {
        ielec = enbcoulTAB;
    }
    else if(fr->eeltype == eelRF_NEC)
    {
        ielec = enbcoulRF;
    }
    else
    {
        ielec = enbcoulOOR;
    }
    
    if(fr->bvdwtab)
    {
        ivdw = enbvdwTAB;
    }
    else
    {
        ivdw = enbvdwLJ;
    }
    
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
                gmx_nb_free_energy_kernel(ielec,
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
                /* Not perturbed - call kernel 330 */
                nb_kernel330( &i1,
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
        }
    }
    return 0.0;
}


