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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <gmx_thread.h>


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
#include "fatal.h"
#include "physics.h"
#include "force.h"
#include "bondf.h"
#include "nrnb.h"
#include "smalloc.h"

#include "nb_kernel/nb_kernel.h"
#include "nb_kernel/nb_kernel330.h"
#include "nb_free_energy.h"


#ifdef GMX_PPC_ALTIVEC   
#include "nb_kernel_ppc_altivec/nb_kernel_ppc_altivec.h"
#endif

#ifdef GMX_IA32_3DNOW  
#include "nb_kernel_ia32_3dnow/nb_kernel_ia32_3dnow.h"
#endif

#ifdef GMX_IA32_SSE   
#include "nb_kernel_ia32_sse/nb_kernel_ia32_sse.h"
#endif

#ifdef GMX_IA32_SSE2    
#include "nb_kernel_ia32_sse2/nb_kernel_ia32_sse2.h"
#endif

#ifdef GMX_X86_64_SSE    
#include "nb_kernel_x86_64_sse/nb_kernel_x86_64_sse.h"
#endif

#ifdef GMX_X86_64_SSE2    
#include "nb_kernel_x86_64_sse2/nb_kernel_x86_64_sse2.h"
#endif

#if (defined GMX_IA64_ASM && defined GMX_DOUBLE) 
#include "nb_kernel_ia64_double/nb_kernel_ia64_double.h"
#endif

#if (defined GMX_IA64_ASM && !defined GMX_DOUBLE)
#include "nb_kernel_ia64_single/nb_kernel_ia64_single.h"
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


static void
setup_kernels(FILE *log)
{
    int i;
    
    snew(nb_kernel_list,eNR_NBKERNEL_NR);
    
    /* Note that later calls overwrite earlier, so the preferred (fastest)
     * version should be at the end. For instance, we call SSE after 3DNow.
     */
    
    for(i=0;i<eNR_NBKERNEL_NR;i++)
        nb_kernel_list[i] = NULL;
    
	if(log)
	    fprintf(log,"Configuring nonbonded kernels...\n");
	
    nb_kernel_setup(log,nb_kernel_list);

    if(getenv("NOASSEMBLYLOOPS") != NULL)
    {
        if(log)
            fprintf(log,
                    "Found environment variable NOASSEMBLYLOOPS.\n"
                    "Disabling all SSE/SSE2/3DNow/Altivec/ia64 asm support.\n");
        return;
    }
    
#ifdef GMX_DOUBLE

    /* Double precision */    

#ifdef GMX_IA32_SSE2    
    nb_kernel_setup_ia32_sse2(log,nb_kernel_list);
#elif defined GMX_X86_64_SSE2 
    nb_kernel_setup_x86_64_sse2(log,nb_kernel_list);
#elif defined GMX_IA64_ASM
    nb_kernel_setup_ia64_double(log,nb_kernel_list);
#endif
    
#else
    /* Single */

    /* Older Athlons only support 3DNow, so check that first, and 
     * then newer AthlonXP/Opteron CPUs with SSE support will override
     * it in the SSE check.
     */   
#ifdef GMX_IA32_3DNOW
    nb_kernel_setup_ia32_3dnow(log,nb_kernel_list);
#endif
#ifdef GMX_IA32_SSE    
    nb_kernel_setup_ia32_sse(log,nb_kernel_list);
#elif defined GMX_X86_64_SSE   
    nb_kernel_setup_x86_64_sse(log,nb_kernel_list);
#elif defined GMX_PPC_ALTIVEC
    nb_kernel_setup_ppc_altivec(log,nb_kernel_list);
#elif defined GMX_IA64_ASM
    nb_kernel_setup_ia64_single(log,nb_kernel_list);
#endif
    
#endif /* precision */

	if(log)
	    fprintf(log,"\n\n");
}


void do_nonbonded(FILE *fplog,t_commrec *cr,t_forcerec *fr,
                  rvec x[],rvec f[],t_mdatoms *mdatoms,
                  real egnb[],real egcoul[],rvec box_size,
                  t_nrnb *nrnb,real lambda,real *dvdlambda,
                  bool bLR,int nls,int eNL,bool bDoForces)
{
    t_nblist *      nlist;
    real *          fshift;
    int             i,i0,i1,nrnb_ind,sz;
    bool            bWater;
    nb_kernel_t *   kernelptr;
    FILE *          fp;
    int             wateratoms;
    int             nthreads = 1;
    int             tabletype;
    int             outeriter,inneriter;
    real *          tabledata = NULL;
    
    if(fr->solvent_opt == esolSPC)
    {
        wateratoms = 3;
    }
    else if(fr->solvent_opt == esolTIP4P)
    {
        wateratoms = 4;
    }
    else
    {
        wateratoms = 1;
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
    
    if(nb_kernel_list == NULL)
    {
        setup_kernels(fplog);    
    }
    
    if (bLR)
    {
        fshift = fr->fshift_twin[0];
    }
    else
    {
        fshift = fr->fshift[0];
    }
    
    for(i=i0; (i<i1); i++) 
    {
        outeriter = inneriter = 0;
        
        if (bLR) 
        {
            nlist  = &(fr->nblists->nlist_lr[i]);
        }
        else
        {
            nlist = &(fr->nblists->nlist_sr[i]);
        }
        
        if (nlist->nri > 0) 
        {
            nrnb_ind = nlist->il_code;
            
            if(nrnb_ind==eNR_NBKERNEL_FREE_ENERGY)
            {
                /* generic free energy */
                if(fr->bvdwtab && fr->bcoultab)
                {
                    tabledata = fr->nblists->tab.tab;
                }
                else if(fr->bvdwtab)
                {
                    tabledata = fr->nblists->vdwtab;
                }
                else if(fr->bcoultab)
                {
                    tabledata = fr->nblists->coultab;
                }
                else
                {
                    tabledata = NULL;
                }                
            }
            else
            {
                /* normal kernels, not free energy */
                tabletype = nb_kernel_table[nrnb_ind];
                
                if(tabletype == TABLE_COMBINED)
                {
                    tabledata = fr->nblists->tab.tab;
                }
                else if(tabletype == TABLE_COUL)
                {
                    tabledata = fr->nblists->coultab;
                }
                else if(tabletype == TABLE_VDW)
                {
                    tabledata = fr->nblists->vdwtab;
                }
                else
                {
                    tabledata = NULL;
                }
            }
            
            nlist->count = 0;
            
            /* Is this a free energy loop? */
            if(nlist->free_energy)
            {
                if(nlist->ivdw==2)
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
                                          egcoul,
                                          mdatoms->typeA,
                                          mdatoms->typeB,
                                          fr->ntype,
                                          fr->nbfp,
                                          egnb,
                                          fr->nblists->tab.scale,
                                          tabledata,
                                          lambda,
                                          dvdlambda,
                                          fr->sc_alpha,
                                          fr->sc_sigma6,
                                          &outeriter,
                                          &inneriter);
            }
            else
            {
                /* Not free energy - call nonbonded kernel from function pointer */
                kernelptr = nb_kernel_list[nrnb_ind];
                
                if(kernelptr == NULL)
                {
                    /* Call generic nonbonded kernel */
                    /*
                    gmx_nb_generic_kernel(nlist,
                                          x,
                                          fr,
                                          f,
                                          fshift,
                                          mdatoms,
                                          lambda,
                                          dvdlambda,
                                          egcoul,
                                          egnb,
                                          &(fr->tab.scale),
                                          tabledata,
                                          &nthreads,
                                          &(nlist->count),
                                          nlist->mtx,
                                          &outeriter,
                                          &inneriter);
                     */
                }
                
                /* Call the appropriate nonbonded kernel function */
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
                              &(fr->nblists->tab.scale),
                              tabledata,
                              NULL,
                              NULL,
                              NULL,
                              NULL,
                              &nthreads,
                              &(nlist->count),
                              nlist->mtx,
                              &outeriter,
                              &inneriter,
                              NULL);
            }
            
            /* Update flop accounting */
            
            /* Outer loop in kernel */
            if (nlist->solvent_opt==enlistWATER)
            {
                inc_nrnb(nrnb,eNR_NBKERNEL_OUTER,wateratoms*outeriter);
            }
            else if (nlist->solvent_opt==enlistWATERWATER)
            {
                inc_nrnb(nrnb,eNR_NBKERNEL_OUTER,wateratoms*wateratoms*outeriter);
            }
            else
            {
                inc_nrnb(nrnb,eNR_NBKERNEL_OUTER,outeriter);
            }
            /* inner loop in kernel */
            inc_nrnb(nrnb,nrnb_ind,inneriter);
        }
    }
}



real 
do_nonbonded14(int nbonds,const t_iatom iatoms[],const t_iparams iparams[],
               const rvec x[],rvec f[],rvec fshift[],
               const t_pbc *pbc,const t_graph *g,
               real lambda,real *dvdlambda,
               const t_mdatoms *md,
               const t_forcerec *fr,int ngrp,real egnb[],real egcoul[])
{
    static    bool bWarn=FALSE;
    bool      bFullPBC;
    real      eps,r2,rtab2;
    rvec      dx,x14[2],f14[2];
    int       i,ai,aj,itype;
    int       typeA[2]={0,0},typeB[2]={0,1};
    real      chargeA[2],chargeB[2];
    int       gid,shift_vir,shift_f;
    int       j_index[] = { 0, 1 };
    int       i0=0,i1=1,i2=2;
    ivec      dt;
    int       outeriter,inneriter;
    int       nthreads = 1;
    int       count;
    real      krf,crf,tabscale;
    real *    nbfp;
    t_nblist  tmplist;
    int       icoul,ivdw;
    
#if GMX_THREADS
    pthread_mutex_t mtx;
#else
    void *    mtx = NULL;
#endif

    
#if GMX_THREADS
    pthread_mutex_initialize(&mtx);
#endif
    
    krf = fr->k_rf;
    crf = fr->c_rf;
    tabscale = fr->tab14.scale;
    
    /* Determine the values for icoul/ivdw. */
    if(fr->bcoultab)
    {
        icoul = 3;
    }
    else if(EEL_RF(fr->eeltype))
    {
        icoul = 2;
    }
    else 
    {
        icoul = 1;
    }
    
    if(fr->bvdwtab)
    {
        ivdw = 3;
    }
    else if(fr->bBHAM)
    {
        ivdw = 2;
    }
    else 
    {
        ivdw = 1;
    }
    
    
    /* We don't do SSE or altivec here, due to large overhead for 4-fold 
     * unrolling on short lists 
     */
    
    bFullPBC = (fr->ePBC == epbcFULL);
    
    /* Reaction field stuff */  
    eps   = fr->epsfac*fr->fudgeQQ;
    
    rtab2 = sqr(fr->tab14.r);
    
    for(i=0; (i<nbonds); ) 
    {
        itype = iatoms[i++];
        ai    = iatoms[i++];
        aj    = iatoms[i++];
        
        if (!bFullPBC) 
        {
            /* This is a bonded interaction, atoms are in the same box */
            shift_f = CENTRAL;
            r2 = distance2(x[ai],x[aj]);
        }
        else 
        {
            /* Apply full periodic boundary conditions */
            shift_f = pbc_dx(pbc,x[ai],x[aj],dx);
            r2 = norm2(dx);
        }

        if (r2 >= rtab2) 
        {
            if (!bWarn) 
            {
                fprintf(stderr,"Warning: 1-4 interaction between %d and %d "
                        "at distance %.3f which is larger than the 1-4 table size %.3f nm\n", ai+1, aj+1, sqrt(r2), sqrt(rtab2));
                fprintf(stderr,"These are ignored for the rest of the simulation\n");
		fprintf(stderr,"This usually means your system is exploding,\n"
			"if not, you should increase table-extension in your mdp file\n");
                bWarn = TRUE;
            }
            if (debug) 
	      fprintf(debug,"%8f %8f %8f\n%8f %8f %8f\n1-4 (%d,%d) interaction not within cut-off! r=%g. Ignored",
		      x[ai][XX],x[ai][YY],x[ai][ZZ],
		      x[aj][XX],x[aj][YY],x[aj][ZZ],
		      (int)ai+1,(int)aj+1,sqrt(r2));
        }
        else 
        {
            chargeA[0] = md->chargeA[ai];
            chargeA[1] = md->chargeA[aj];
            gid  = GID(md->cENER[ai],md->cENER[aj],ngrp);
            copy_rvec(x[ai],x14[0]);
            copy_rvec(x[aj],x14[1]);
            clear_rvec(f14[0]);
            clear_rvec(f14[1]);
#ifdef DEBUG
            fprintf(debug,"LJ14: grp-i=%2d, grp-j=%2d, ngrp=%2d, GID=%d\n",
                    md->cENER[ai],md->cENER[aj],ngrp,gid);
#endif
            
            outeriter = inneriter = count = 0;
            nbfp = (real *)&(iparams[itype].lj14.c6A);

            if (fr->efep != efepNO &&
                (md->bPerturbed[ai] || md->bPerturbed[aj] ||
                 iparams[itype].lj14.c6A != iparams[itype].lj14.c6B ||
                 iparams[itype].lj14.c12A != iparams[itype].lj14.c12B)) 
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
                if(ivdw==2)
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
                                          egcoul,
                                          typeA,
                                          typeB,
                                          i1,
                                          nbfp,
                                          egnb,
                                          tabscale,
                                          fr->tab14.tab,
                                          lambda,
                                          dvdlambda,
                                          fr->sc_alpha,
                                          fr->sc_sigma6,
                                          &outeriter,
                                          &inneriter);
            }
            else 
            { 
                /* Not perturbed - call kernel 330 */
                F77_OR_C_FUNC_(nb_kernel330,NB_KERNEL330)
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
                      &i1,
                      nbfp,
                      egnb,
                      &tabscale,
                      fr->tab14.tab,
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

	    /* flops: eNR_KERNEL_OUTER + eNR_KERNEL330 + 12 */
        }
    }
    return 0.0;
}


