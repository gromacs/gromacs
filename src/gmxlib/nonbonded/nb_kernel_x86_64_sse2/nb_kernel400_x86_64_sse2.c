/*
 * Copyright (c) Erik Lindahl, David van der Spoel 2003
 * 
 * This file is generated automatically at compile time
 * by the program mknb in the Gromacs distribution.
 *
 * Options used when generation this file:
 * Language:         c
 * Precision:        double
 * Threads:          no
 * Software invsqrt: yes
 * Prefetch forces:  no
 * Comments:         no
 */
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#include<math.h>
#include<vec.h>

#include <xmmintrin.h>
#include <emmintrin.h>

#include <gmx_sse2_double.h>

/* get gmx_gbdata_t */
#include "../nb_kerneltype.h"

#include "nb_kernel400_x86_64_sse2.h"


void nb_kernel400_x86_64_sse2(int *           p_nri,
                              int *           iinr,
                              int *           jindex,
                              int *           jjnr,
                              int *           shift,
                              double *         shiftvec,
                              double *         fshift,
                              int *           gid,
                              double *         pos,
                              double *         faction,
                              double *         charge,
                              double *         p_facel,
                              double *         p_krf,
                              double *         p_crf,
                              double *         vc,
                              int *           type,
                              int *           p_ntype,
                              double *         vdwparam,
                              double *         vvdw,
                              double *         p_tabscale,
                              double *         VFtab,
                              double *         invsqrta,
                              double *         dvda,
                              double *         p_gbtabscale,
                              double *         GBtab,
                              int *           p_nthreads,
                              int *           count,
                              void *          mtx,
                              int *           outeriter,
                              int *           inneriter,
                              double *         work)
{
  int           nri,nthreads;
  int           n,ii,is3,ii3,k,nj0,nj1,ggid;
  double        shX,shY,shZ;
  int           jnrA,jnrB;
  int           j3A,j3B;
	gmx_gbdata_t *gbdata;
	double *      gpol;
    
	__m128d  iq,qq,jq,isai;
	__m128d  ix,iy,iz;
	__m128d  jx,jy,jz;
	__m128d  dx,dy,dz;
	__m128d  vctot,vgbtot,dvdasum,gbfactor;
	__m128d  fix,fiy,fiz,tx,ty,tz,rsq;
	__m128d  rinv,isaj,isaprod;
	__m128d  vcoul,fscal,gbscale;
	__m128d  rinvsq,r,rtab;
	__m128d  eps,Y,F,G,H;
	__m128d  vgb,fijGB,dvdatmp;
	__m128d  facel,gbtabscale,dvdaj;
	__m128i  n0, nnn;
	
	const __m128d neg        = _mm_set1_pd(-1.0);
	const __m128d zero       = _mm_set1_pd(0.0);
	const __m128d minushalf  = _mm_set1_pd(-0.5);
	const __m128d two        = _mm_set1_pd(2.0);
	
	gbdata     = (gmx_gbdata_t *)work;
	gpol       = gbdata->gpol;
    
	nri        = *p_nri;
    
  gbfactor   = _mm_set1_pd( - ((1.0/gbdata->epsilon_r) - (1.0/gbdata->gb_epsilon_solvent)));     
  gbtabscale = _mm_load1_pd(p_gbtabscale);  
  facel      = _mm_load1_pd(p_facel);
  
  nj1         = 0;
  jnrA = jnrB = 0;
  j3A = j3B   = 0;
  jx          = _mm_setzero_pd();
  jy          = _mm_setzero_pd();
  jz          = _mm_setzero_pd();
	
	for(n=0;n<nri;n++)
	{
    is3              = 3*shift[n];     
    shX              = shiftvec[is3];  
    shY              = shiftvec[is3+1];
    shZ              = shiftvec[is3+2];
    nj0              = jindex[n];      
    nj1              = jindex[n+1];    
    ii               = iinr[n];        
    ii3              = 3*ii;           
		
		ix               = _mm_set1_pd(shX+pos[ii3+0]);
		iy               = _mm_set1_pd(shY+pos[ii3+1]);
		iz               = _mm_set1_pd(shZ+pos[ii3+2]);
    
		iq               = _mm_load1_pd(charge+ii);
		iq               = _mm_mul_pd(iq,facel);
    
		isai             = _mm_load1_pd(invsqrta+ii);
        		
		vctot            = _mm_setzero_pd();
		vgbtot           = _mm_setzero_pd();
		dvdasum          = _mm_setzero_pd();
		fix              = _mm_setzero_pd();
		fiy              = _mm_setzero_pd();
		fiz              = _mm_setzero_pd();
                
		for(k=nj0;k<nj1-1; k+=2)
		{
			jnrA    = jjnr[k];
			jnrB    = jjnr[k+1];
			
			j3A     = jnrA * 3;
			j3B     = jnrB * 3;
      
      GMX_MM_LOAD_1RVEC_2POINTERS_PD(pos+j3A,pos+j3B,jx,jy,jz);
            
			dx           = _mm_sub_pd(ix,jx);
			dy           = _mm_sub_pd(iy,jy);
			dz           = _mm_sub_pd(iz,jz);
            
      rsq          = gmx_mm_calc_rsq_pd(dx,dy,dz);
      
      rinv         = gmx_mm_invsqrt_pd(rsq);
 			rinvsq       = _mm_mul_pd(rinv,rinv);
      
			/***********************************/
			/* INTERACTION SECTION STARTS HERE */
			/***********************************/
			GMX_MM_LOAD_2VALUES_PD(charge+jnrA,charge+jnrB,jq);
			GMX_MM_LOAD_2VALUES_PD(invsqrta+jnrA,invsqrta+jnrB,isaj);
            
			isaprod      = _mm_mul_pd(isai,isaj);
			qq           = _mm_mul_pd(iq,jq);            
			vcoul        = _mm_mul_pd(qq,rinv);
			fscal        = _mm_mul_pd(vcoul,rinv);                                 
      vctot        = _mm_add_pd(vctot,vcoul);
            
            /* Polarization interaction */
			qq           = _mm_mul_pd(qq,_mm_mul_pd(isaprod,gbfactor));
			gbscale      = _mm_mul_pd(isaprod,gbtabscale);
            
 			/* Calculate GB table index */
			r            = _mm_mul_pd(rsq,rinv);
			rtab         = _mm_mul_pd(r,gbscale);
			
			n0		     = _mm_cvttpd_epi32(rtab);
			eps	     	 = _mm_sub_pd(rtab,_mm_cvtepi32_pd(n0));
			nnn		     = _mm_slli_epi32(n0,2);
			
      /* the tables are 16-byte aligned, so we can use _mm_load_pd */			
      Y            = _mm_load_pd(GBtab+(gmx_mm_extract_epi32(nnn,0))); 
      F            = _mm_load_pd(GBtab+(gmx_mm_extract_epi32(nnn,1)));
      GMX_MM_TRANSPOSE2_PD(Y,F);
      G            = _mm_load_pd(GBtab+(gmx_mm_extract_epi32(nnn,0))+2); 
      H            = _mm_load_pd(GBtab+(gmx_mm_extract_epi32(nnn,1))+2);
      GMX_MM_TRANSPOSE2_PD(G,H);
      
      G       = _mm_mul_pd(G,eps);
      H       = _mm_mul_pd(H, _mm_mul_pd(eps,eps) );
      F       = _mm_add_pd(F, _mm_add_pd( G , H ) );
      Y       = _mm_add_pd(Y, _mm_mul_pd(F, eps));
      F       = _mm_add_pd(F, _mm_add_pd(G , _mm_mul_pd(H,two)));
      vgb     = _mm_mul_pd(Y, qq);           
      fijGB   = _mm_mul_pd(F, _mm_mul_pd(qq,gbscale));
      
      dvdatmp = _mm_mul_pd(_mm_add_pd(vgb, _mm_mul_pd(fijGB,r)) , minushalf);
      
      vgbtot  = _mm_add_pd(vgbtot, vgb);
      
      dvdasum = _mm_add_pd(dvdasum, dvdatmp);
      dvdatmp = _mm_mul_pd(dvdatmp, _mm_mul_pd(isaj,isaj));
      
      GMX_MM_INCREMENT_2VALUES_PD(dvda+jnrA,dvda+jnrB,dvdatmp);
      
      fscal        = _mm_mul_pd( _mm_sub_pd( fscal, fijGB),rinv );
      
      /***********************************/
			/*  INTERACTION SECTION ENDS HERE  */
			/***********************************/
      
      /* Calculate temporary vectorial force */
      tx           = _mm_mul_pd(fscal,dx);
      ty           = _mm_mul_pd(fscal,dy);
      tz           = _mm_mul_pd(fscal,dz);
      
      /* Increment i atom force */
      fix          = _mm_add_pd(fix,tx);
      fiy          = _mm_add_pd(fiy,ty);
      fiz          = _mm_add_pd(fiz,tz);
      
      /* Store j forces back */
			GMX_MM_DECREMENT_1RVEC_2POINTERS_PD(faction+j3A,faction+j3B,tx,ty,tz);
		}
		
		/* In double precision, offset can only be either 0 or 1 */
		if(k<nj1)
		{
			jnrA    = jjnr[k];
			j3A     = jnrA * 3;
      
      GMX_MM_LOAD_1RVEC_1POINTER_PD(pos+j3A,jx,jy,jz);
      
			dx           = _mm_sub_sd(ix,jx);
			dy           = _mm_sub_sd(iy,jy);
			dz           = _mm_sub_sd(iz,jz);
      
      rsq          = gmx_mm_calc_rsq_pd(dx,dy,dz);
      
      rinv         = gmx_mm_invsqrt_pd(rsq);
 			rinvsq       = _mm_mul_sd(rinv,rinv);
      
      /* These reason for zeroing these variables here is for fixing bug 585
       * What happens is that __m128d _mm_add_sd(a,b) gives back r0=a[0]+b[0],
       * and r1=0, but it should be r1=a[1]. 
       * This might be a compiler issue (tested with gcc-4.1.3 and -O3).
       * To work around it, we zero these variables and use _mm_add_pd (**) instead
       * Note that the only variables that get affected are the energies since
       * the total sum needs to be correct 
       */
      vcoul        = _mm_setzero_pd();
      dvdatmp      = _mm_setzero_pd();
      vgb          = _mm_setzero_pd();
      
			/***********************************/
			/* INTERACTION SECTION STARTS HERE */
			/***********************************/
			GMX_MM_LOAD_1VALUE_PD(charge+jnrA,jq);
			GMX_MM_LOAD_1VALUE_PD(invsqrta+jnrA,isaj);
      
			isaprod      = _mm_mul_sd(isai,isaj);
      /* Since we need _mm_add_pd below, the order here og jq,iq becomes important */
			qq           = _mm_mul_sd(jq,iq);  
      vcoul        = _mm_mul_sd(qq,rinv);
      fscal        = _mm_mul_sd(vcoul,rinv);                                 
      vctot        = _mm_add_pd(vctot,vcoul); /* (**) */
      
      /* Polarization interaction */
			qq           = _mm_mul_sd(qq,_mm_mul_sd(isaprod,gbfactor));
			gbscale      = _mm_mul_sd(isaprod,gbtabscale);
      
 			/* Calculate GB table index */
			r            = _mm_mul_sd(rsq,rinv);
			rtab         = _mm_mul_sd(r,gbscale);
      
			n0		     = _mm_cvttpd_epi32(rtab);
			eps	     	 = _mm_sub_sd(rtab,_mm_cvtepi32_pd(n0));
			nnn		     = _mm_slli_epi32(n0,2);
			
      /* the tables are 16-byte aligned, so we can use _mm_load_pd */			
      Y            = _mm_load_pd(GBtab+(gmx_mm_extract_epi32(nnn,0))); 
      F            = _mm_setzero_pd();
      GMX_MM_TRANSPOSE2_PD(Y,F);
      G            = _mm_load_pd(GBtab+(gmx_mm_extract_epi32(nnn,0))+2); 
      H            = _mm_setzero_pd();
      GMX_MM_TRANSPOSE2_PD(G,H);

      G       = _mm_mul_sd(G,eps);
      H       = _mm_mul_sd(H, _mm_mul_sd(eps,eps) );
      F       = _mm_add_sd(F, _mm_add_sd( G , H ) );
      Y       = _mm_add_sd(Y, _mm_mul_sd(F, eps));
      F       = _mm_add_sd(F, _mm_add_sd(G , _mm_mul_sd(H,two)));
      vgb     = _mm_mul_sd(Y, qq);           
      fijGB   = _mm_mul_sd(F, _mm_mul_sd(qq,gbscale));
      dvdatmp = _mm_mul_sd(_mm_add_sd(vgb, _mm_mul_sd(fijGB,r)) , minushalf);
      
      vgbtot  = _mm_add_pd(vgbtot, vgb); /* (**) */
      
      dvdasum = _mm_add_pd(dvdasum, dvdatmp); /* (**) */
      dvdatmp = _mm_mul_sd(dvdatmp, _mm_mul_sd(isaj,isaj));
      
      GMX_MM_INCREMENT_1VALUE_PD(dvda+jnrA,dvdatmp);
			
      fscal        = _mm_mul_sd( _mm_sub_sd( fscal, fijGB),rinv );
      
      /***********************************/
			/*  INTERACTION SECTION ENDS HERE  */
			/***********************************/
      
      /* Calculate temporary vectorial force */
      tx           = _mm_mul_sd(fscal,dx);
      ty           = _mm_mul_sd(fscal,dy);
      tz           = _mm_mul_sd(fscal,dz);
      
      /* Increment i atom force */
      fix          = _mm_add_sd(fix,tx);
      fiy          = _mm_add_sd(fiy,ty);
      fiz          = _mm_add_sd(fiz,tz);
      
      /* Store j forces back */
			GMX_MM_DECREMENT_1RVEC_1POINTER_PD(faction+j3A,tx,ty,tz);
		}
		
    dvdasum = _mm_mul_pd(dvdasum, _mm_mul_pd(isai,isai));
    gmx_mm_update_iforce_1atom_pd(&fix,&fiy,&fiz,faction+ii3,fshift+is3);
    
    ggid     = gid[n];         
   
    gmx_mm_update_1pot_pd(vctot,vc+ggid);
    gmx_mm_update_1pot_pd(vgbtot,gpol+ggid);
    gmx_mm_update_1pot_pd(dvdasum,dvda+ii);
  }
  
	*outeriter   = nri;            
  *inneriter   = nj1; 	
}
