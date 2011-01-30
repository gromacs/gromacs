/*
 * Copyright (c) Erik Lindahl, David van der Spoel 2003
 * 
 * This file is generated automatically at compile time
 * by the program mknb in the Gromacs distribution.
 *
 * Options used when generation this file:
 * Language:         c
 * Precision:        single
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

#include <gmx_sse2_single.h>

/* get gmx_gbdata_t */
#include "../nb_kerneltype.h"

#include "nb_kernel410_x86_64_sse.h"



void nb_kernel410_x86_64_sse(int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    float *         shiftvec,
                    float *         fshift,
                    int *           gid,
                    float *         pos,
                    float *         faction,
                    float *         charge,
                    float *         p_facel,
                    float *         p_krf,
                    float *         p_crf,
                    float *         vc,
                    int *           type,
                    int *           p_ntype,
                    float *         vdwparam,
                    float *         vvdw,
                    float *         p_tabscale,
                    float *         VFtab,
                    float *         invsqrta,
                    float *         dvda,
                    float *         p_gbtabscale,
                    float *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    float *         work)
{
    int           nri,ntype,nthreads;
    int           n,ii,is3,ii3,k,nj0,nj1,ggid;
    float         shX,shY,shZ;
	int			  offset,nti;
    int           jnrA,jnrB,jnrC,jnrD,j3A,j3B,j3C,j3D;
    int           jnrE,jnrF,jnrG,jnrH,j3E,j3F,j3G,j3H;
	int           tjA,tjB,tjC,tjD;
	int           tjE,tjF,tjG,tjH;
	gmx_gbdata_t *gbdata;
	float *        gpol;
	
	__m128   iq,qq,qqB,jq,jqB,isai;
	__m128   ix,iy,iz;
	__m128   jx,jy,jz;
	__m128   jxB,jyB,jzB;
	__m128   dx,dy,dz;
	__m128   dxB,dyB,dzB;
	__m128   vctot,vvdwtot,vgbtot,dvdasum,gbfactor;
	__m128   fix,fiy,fiz,tx,ty,tz,rsq;
	__m128   fixB,fiyB,fizB,txB,tyB,tzB,rsqB;
	__m128   rinv,isaj,isaprod;
	__m128   rinvB,isajB,isaprodB;
	__m128   vcoul,vcoulB,fscal,fscalB,gbscale,gbscaleB,c6,c6B,c12,c12B;
	__m128   rinvsq,r,rtab;
	__m128   rinvsqB,rB,rtabB;
	__m128   eps,Y,F,G,H;
	__m128   epsB,YB,FB,GB,HB;
	__m128   vgb,vgbB,fijGB,fijGBB,dvdatmp,dvdatmpB;
	__m128   rinvsix,vvdw6,vvdw12;
	__m128   rinvsixB,vvdw6B,vvdw12B;
	__m128   facel,gbtabscale,mask,dvdaj;
    
    __m128   mask1 = gmx_mm_castsi128_ps( _mm_set_epi32(0, 0, 0, 0xffffffff) );
	__m128   mask2 = gmx_mm_castsi128_ps( _mm_set_epi32(0, 0, 0xffffffff, 0xffffffff) );
	__m128   mask3 = gmx_mm_castsi128_ps( _mm_set_epi32(0, 0xffffffff, 0xffffffff, 0xffffffff) );
    
	__m128i  n0, nnn;
	__m128i  n0B, nnnB;
	
	const __m128 neg        = _mm_set1_ps(-1.0f);
	const __m128 zero       = _mm_set1_ps(0.0f);
	const __m128 minushalf  = _mm_set1_ps(-0.5f);
	const __m128 two        = _mm_set1_ps(2.0f);
	const __m128 six        = _mm_set1_ps(6.0f);
	const __m128 twelve     = _mm_set1_ps(12.0f);

	gbdata          = (gmx_gbdata_t *)work;
	gpol            = gbdata->gpol;
	
	nri              = *p_nri;         
    ntype            = *p_ntype;       
    
    gbfactor         = _mm_set1_ps( - ((1.0/gbdata->epsilon_r) - (1.0/gbdata->gb_epsilon_solvent)));     
    gbtabscale       = _mm_load1_ps(p_gbtabscale);  
    facel            = _mm_load1_ps(p_facel);

    nj1              = 0;
    jnrA = jnrB = jnrC = jnrD = 0;
    j3A = j3B = j3C = j3D = 0;
    jx               = _mm_setzero_ps();
    jy               = _mm_setzero_ps();
    jz               = _mm_setzero_ps();
    c6               = _mm_setzero_ps();
    c12              = _mm_setzero_ps();

    
    for(n=0; (n<nri); n++)
    {
        is3              = 3*shift[n];     
        shX              = shiftvec[is3];  
        shY              = shiftvec[is3+1];
        shZ              = shiftvec[is3+2];
        nj0              = jindex[n];      
        nj1              = jindex[n+1];    
        ii               = iinr[n];        
        ii3              = 3*ii;           
        
		ix               = _mm_set1_ps(shX+pos[ii3+0]);
		iy               = _mm_set1_ps(shY+pos[ii3+1]);
		iz               = _mm_set1_ps(shZ+pos[ii3+2]);
	
		iq               = _mm_load1_ps(charge+ii);
		iq               = _mm_mul_ps(iq,facel);
			
		isai             = _mm_load1_ps(invsqrta+ii);
								
		nti              = 2*ntype*type[ii];
		
		vctot            = _mm_setzero_ps();
		vvdwtot          = _mm_setzero_ps();
		vgbtot           = _mm_setzero_ps();
		dvdasum          = _mm_setzero_ps();
		fix              = _mm_setzero_ps();
		fiy              = _mm_setzero_ps();
		fiz              = _mm_setzero_ps();
       
        for(k=nj0; k<nj1-7; k+=8)
		{
			jnrA        = jjnr[k];   
			jnrB        = jjnr[k+1];
			jnrC        = jjnr[k+2];
			jnrD        = jjnr[k+3];
			jnrE        = jjnr[k+4];   
			jnrF        = jjnr[k+5];
			jnrG        = jjnr[k+6];
			jnrH        = jjnr[k+7];
            
            j3A         = 3*jnrA;  
			j3B         = 3*jnrB;
			j3C         = 3*jnrC;
			j3D         = 3*jnrD;
            j3E         = 3*jnrE;  
			j3F         = 3*jnrF;
			j3G         = 3*jnrG;
			j3H         = 3*jnrH;
            

            GMX_MM_LOAD_1RVEC_4POINTERS_PS(pos+j3A,pos+j3B,pos+j3C,pos+j3D,jx,jy,jz);
            GMX_MM_LOAD_1RVEC_4POINTERS_PS(pos+j3E,pos+j3F,pos+j3G,pos+j3H,jxB,jyB,jzB);

			dx           = _mm_sub_ps(ix,jx);
			dy           = _mm_sub_ps(iy,jy);
			dz           = _mm_sub_ps(iz,jz);
			dxB          = _mm_sub_ps(ix,jxB);
			dyB          = _mm_sub_ps(iy,jyB);
			dzB          = _mm_sub_ps(iz,jzB);
			
			rsq          = gmx_mm_calc_rsq_ps(dx,dy,dz);
			rsqB         = gmx_mm_calc_rsq_ps(dxB,dyB,dzB);
			
			/* invsqrt */								
			rinv         = gmx_mm_invsqrt_ps(rsq);
			rinvB        = gmx_mm_invsqrt_ps(rsqB);
 			rinvsq       = _mm_mul_ps(rinv,rinv);
 			rinvsqB      = _mm_mul_ps(rinvB,rinvB);

			/***********************************/
			/* INTERACTION SECTION STARTS HERE */
			/***********************************/
			GMX_MM_LOAD_4VALUES_PS(charge+jnrA,charge+jnrB,charge+jnrC,charge+jnrD,jq);
			GMX_MM_LOAD_4VALUES_PS(charge+jnrE,charge+jnrF,charge+jnrG,charge+jnrH,jqB);
			GMX_MM_LOAD_4VALUES_PS(invsqrta+jnrA,invsqrta+jnrB,invsqrta+jnrC,invsqrta+jnrD,isaj);
			GMX_MM_LOAD_4VALUES_PS(invsqrta+jnrE,invsqrta+jnrF,invsqrta+jnrG,invsqrta+jnrH,isajB);

            /* Lennard-Jones */
            tjA          = nti+2*type[jnrA];
			tjB          = nti+2*type[jnrB];
			tjC          = nti+2*type[jnrC];
			tjD          = nti+2*type[jnrD];
            tjE          = nti+2*type[jnrE];
			tjF          = nti+2*type[jnrF];
			tjG          = nti+2*type[jnrG];
			tjH          = nti+2*type[jnrH];
            
            GMX_MM_LOAD_4PAIRS_PS(vdwparam+tjA,vdwparam+tjB,vdwparam+tjC,vdwparam+tjD,c6,c12);
            GMX_MM_LOAD_4PAIRS_PS(vdwparam+tjE,vdwparam+tjF,vdwparam+tjG,vdwparam+tjH,c6B,c12B);
			
			isaprod      = _mm_mul_ps(isai,isaj);
			isaprodB     = _mm_mul_ps(isai,isajB);
			qq           = _mm_mul_ps(iq,jq);            
			qqB          = _mm_mul_ps(iq,jqB);            
			vcoul        = _mm_mul_ps(qq,rinv);
			vcoulB       = _mm_mul_ps(qqB,rinvB);
			fscal        = _mm_mul_ps(vcoul,rinv);                                 
			fscalB       = _mm_mul_ps(vcoulB,rinvB);                                 
            vctot        = _mm_add_ps(vctot,vcoul);
            vctot        = _mm_add_ps(vctot,vcoulB);
            
            /* Polarization interaction */
			qq           = _mm_mul_ps(qq,_mm_mul_ps(isaprod,gbfactor));
			qqB          = _mm_mul_ps(qqB,_mm_mul_ps(isaprodB,gbfactor));
			gbscale      = _mm_mul_ps(isaprod,gbtabscale);
			gbscaleB     = _mm_mul_ps(isaprodB,gbtabscale);

 			/* Calculate GB table index */
			r            = _mm_mul_ps(rsq,rinv);
			rB           = _mm_mul_ps(rsqB,rinvB);
			rtab         = _mm_mul_ps(r,gbscale);
			rtabB        = _mm_mul_ps(rB,gbscaleB);
			
			n0           = _mm_cvttps_epi32(rtab);
			n0B          = _mm_cvttps_epi32(rtabB);
            eps          = _mm_sub_ps(rtab , _mm_cvtepi32_ps(n0) );
            epsB         = _mm_sub_ps(rtabB , _mm_cvtepi32_ps(n0B) );
			nnn          = _mm_slli_epi32(n0,2);
			nnnB         = _mm_slli_epi32(n0B,2);
			
            /* the tables are 16-byte aligned, so we can use _mm_load_ps */			
            Y            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,0)); /* YFGH */
            F            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,1)); /* YFGH */
            G            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,2)); /* YFGH */
            H            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y,F,G,H);
            YB           = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnnB,0)); /* YFGH */
            FB           = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnnB,1)); /* YFGH */
            GB           = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnnB,2)); /* YFGH */
            HB           = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnnB,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(YB,FB,GB,HB);
            
            G       = _mm_mul_ps(G,eps);
            H       = _mm_mul_ps(H, _mm_mul_ps(eps,eps) );
            F       = _mm_add_ps(F, _mm_add_ps( G , H ) );
            Y       = _mm_add_ps(Y, _mm_mul_ps(F, eps));
            F       = _mm_add_ps(F, _mm_add_ps(G , _mm_mul_ps(H,two)));
            vgb     = _mm_mul_ps(Y, qq);           
            fijGB   = _mm_mul_ps(F, _mm_mul_ps(qq,gbscale));

            GB      = _mm_mul_ps(GB,epsB);
            HB      = _mm_mul_ps(HB, _mm_mul_ps(epsB,epsB) );
            FB      = _mm_add_ps(FB, _mm_add_ps( GB , HB ) );
            YB      = _mm_add_ps(YB, _mm_mul_ps(FB, epsB));
            FB      = _mm_add_ps(FB, _mm_add_ps(GB , _mm_mul_ps(HB,two)));
            vgbB    = _mm_mul_ps(YB, qqB);           
            fijGBB  = _mm_mul_ps(FB, _mm_mul_ps(qqB,gbscaleB));
            
            
            dvdatmp = _mm_mul_ps(_mm_add_ps(vgb, _mm_mul_ps(fijGB,r)) , minushalf);
            dvdatmpB = _mm_mul_ps(_mm_add_ps(vgbB, _mm_mul_ps(fijGBB,rB)) , minushalf);

            vgbtot  = _mm_add_ps(vgbtot, vgb);
            vgbtot  = _mm_add_ps(vgbtot, vgbB);
            
            dvdasum = _mm_add_ps(dvdasum, dvdatmp);
            dvdasum = _mm_add_ps(dvdasum, dvdatmpB);
            dvdatmp = _mm_mul_ps(dvdatmp, _mm_mul_ps(isaj,isaj));
            dvdatmpB = _mm_mul_ps(dvdatmpB, _mm_mul_ps(isajB,isajB));
            
            GMX_MM_INCREMENT_4VALUES_PS(dvda+jnrA,dvda+jnrB,dvda+jnrC,dvda+jnrD,dvdatmp);
            GMX_MM_INCREMENT_4VALUES_PS(dvda+jnrE,dvda+jnrF,dvda+jnrG,dvda+jnrH,dvdatmpB);
			
			rinvsix      = _mm_mul_ps(rinvsq,rinvsq);
			rinvsixB     = _mm_mul_ps(rinvsqB,rinvsqB);
			rinvsix      = _mm_mul_ps(rinvsix,rinvsq);
			rinvsixB     = _mm_mul_ps(rinvsixB,rinvsqB);
			
			vvdw6        = _mm_mul_ps(c6,rinvsix);
			vvdw6B       = _mm_mul_ps(c6B,rinvsixB);
			vvdw12       = _mm_mul_ps(c12, _mm_mul_ps(rinvsix,rinvsix));
			vvdw12B      = _mm_mul_ps(c12B, _mm_mul_ps(rinvsixB,rinvsixB));
			vvdwtot      = _mm_add_ps(vvdwtot,_mm_sub_ps(vvdw12,vvdw6));
			vvdwtot      = _mm_add_ps(vvdwtot,_mm_sub_ps(vvdw12B,vvdw6B));

            fscal        = _mm_sub_ps(_mm_mul_ps(rinvsq, 
                                                       _mm_sub_ps(_mm_mul_ps(twelve,vvdw12),
                                                                  _mm_mul_ps(six,vvdw6))),
                                      _mm_mul_ps( _mm_sub_ps( fijGB,fscal),rinv ));

            fscalB       = _mm_sub_ps(_mm_mul_ps(rinvsqB, 
                                                 _mm_sub_ps(_mm_mul_ps(twelve,vvdw12B),
                                                            _mm_mul_ps(six,vvdw6B))),
                                      _mm_mul_ps( _mm_sub_ps( fijGBB,fscalB),rinvB ));
            
            /***********************************/
			/*  INTERACTION SECTION ENDS HERE  */
			/***********************************/

            /* Calculate temporary vectorial force */
            tx           = _mm_mul_ps(fscal,dx);
            ty           = _mm_mul_ps(fscal,dy);
            tz           = _mm_mul_ps(fscal,dz);
            txB          = _mm_mul_ps(fscalB,dxB);
            tyB          = _mm_mul_ps(fscalB,dyB);
            tzB          = _mm_mul_ps(fscalB,dzB);
            
            /* Increment i atom force */
            fix          = _mm_add_ps(fix,tx);
            fiy          = _mm_add_ps(fiy,ty);
            fiz          = _mm_add_ps(fiz,tz);
            fix          = _mm_add_ps(fix,txB);
            fiy          = _mm_add_ps(fiy,tyB);
            fiz          = _mm_add_ps(fiz,tzB);
            
            /* Store j forces back */
			GMX_MM_DECREMENT_1RVEC_4POINTERS_PS(faction+j3A,faction+j3B,faction+j3C,faction+j3D,tx,ty,tz);
			GMX_MM_DECREMENT_1RVEC_4POINTERS_PS(faction+j3E,faction+j3F,faction+j3G,faction+j3H,txB,tyB,tzB);
        }
        for(;k<nj1-3; k+=4)
		{
			jnrA        = jjnr[k];   
			jnrB        = jjnr[k+1];
			jnrC        = jjnr[k+2];
			jnrD        = jjnr[k+3];
            
            j3A         = 3*jnrA;  
			j3B         = 3*jnrB;
			j3C         = 3*jnrC;
			j3D         = 3*jnrD;
            
            
            GMX_MM_LOAD_1RVEC_4POINTERS_PS(pos+j3A,pos+j3B,pos+j3C,pos+j3D,jx,jy,jz);
            
			dx           = _mm_sub_ps(ix,jx);
			dy           = _mm_sub_ps(iy,jy);
			dz           = _mm_sub_ps(iz,jz);
			
			rsq          = gmx_mm_calc_rsq_ps(dx,dy,dz);
			
			/* invsqrt */								
			rinv         = gmx_mm_invsqrt_ps(rsq);
 			rinvsq       = _mm_mul_ps(rinv,rinv);
            
			/***********************************/
			/* INTERACTION SECTION STARTS HERE */
			/***********************************/
			GMX_MM_LOAD_4VALUES_PS(charge+jnrA,charge+jnrB,charge+jnrC,charge+jnrD,jq);
			GMX_MM_LOAD_4VALUES_PS(invsqrta+jnrA,invsqrta+jnrB,invsqrta+jnrC,invsqrta+jnrD,isaj);
            
            /* Lennard-Jones */
            tjA          = nti+2*type[jnrA];
			tjB          = nti+2*type[jnrB];
			tjC          = nti+2*type[jnrC];
			tjD          = nti+2*type[jnrD];
            
            GMX_MM_LOAD_4PAIRS_PS(vdwparam+tjA,vdwparam+tjB,vdwparam+tjC,vdwparam+tjD,c6,c12);
			
			isaprod      = _mm_mul_ps(isai,isaj);
			qq           = _mm_mul_ps(iq,jq);            
			vcoul        = _mm_mul_ps(qq,rinv);
			fscal        = _mm_mul_ps(vcoul,rinv);                                 
            vctot        = _mm_add_ps(vctot,vcoul);
            
            /* Polarization interaction */
			qq           = _mm_mul_ps(qq,_mm_mul_ps(isaprod,gbfactor));
			gbscale      = _mm_mul_ps(isaprod,gbtabscale);
            
 			/* Calculate GB table index */
			r            = _mm_mul_ps(rsq,rinv);
			rtab         = _mm_mul_ps(r,gbscale);
			
			n0           = _mm_cvttps_epi32(rtab);
            eps          = _mm_sub_ps(rtab , _mm_cvtepi32_ps(n0) );
			nnn          = _mm_slli_epi32(n0,2);
			
            /* the tables are 16-byte aligned, so we can use _mm_load_ps */			
            Y            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,0)); /* YFGH */
            F            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,1)); /* YFGH */
            G            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,2)); /* YFGH */
            H            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y,F,G,H);
            
            G       = _mm_mul_ps(G,eps);
            H       = _mm_mul_ps(H, _mm_mul_ps(eps,eps) );
            F       = _mm_add_ps(F, _mm_add_ps( G , H ) );
            Y       = _mm_add_ps(Y, _mm_mul_ps(F, eps));
            F       = _mm_add_ps(F, _mm_add_ps(G , _mm_mul_ps(H,two)));
            vgb     = _mm_mul_ps(Y, qq);           
            fijGB   = _mm_mul_ps(F, _mm_mul_ps(qq,gbscale));
            
            
            dvdatmp = _mm_mul_ps(_mm_add_ps(vgb, _mm_mul_ps(fijGB,r)) , minushalf);
            
            vgbtot  = _mm_add_ps(vgbtot, vgb);
            
            dvdasum = _mm_add_ps(dvdasum, dvdatmp);
            dvdatmp = _mm_mul_ps(dvdatmp, _mm_mul_ps(isaj,isaj));
            
            GMX_MM_INCREMENT_4VALUES_PS(dvda+jnrA,dvda+jnrB,dvda+jnrC,dvda+jnrD,dvdatmp);
            
			
			rinvsix      = _mm_mul_ps(rinvsq,rinvsq);
			rinvsix      = _mm_mul_ps(rinvsix,rinvsq);
			
			vvdw6        = _mm_mul_ps(c6,rinvsix);
			vvdw12       = _mm_mul_ps(c12, _mm_mul_ps(rinvsix,rinvsix));
			vvdwtot      = _mm_add_ps(vvdwtot,_mm_sub_ps(vvdw12,vvdw6));
            
            fscal        = _mm_sub_ps(_mm_mul_ps(rinvsq, 
                                                 _mm_sub_ps(_mm_mul_ps(twelve,vvdw12),
                                                            _mm_mul_ps(six,vvdw6))),
                                      _mm_mul_ps( _mm_sub_ps( fijGB,fscal),rinv ));
            
            /***********************************/
			/*  INTERACTION SECTION ENDS HERE  */
			/***********************************/
            
            /* Calculate temporary vectorial force */
            tx           = _mm_mul_ps(fscal,dx);
            ty           = _mm_mul_ps(fscal,dy);
            tz           = _mm_mul_ps(fscal,dz);
            
            /* Increment i atom force */
            fix          = _mm_add_ps(fix,tx);
            fiy          = _mm_add_ps(fiy,ty);
            fiz          = _mm_add_ps(fiz,tz);
            
            /* Store j forces back */
			GMX_MM_DECREMENT_1RVEC_4POINTERS_PS(faction+j3A,faction+j3B,faction+j3C,faction+j3D,tx,ty,tz);
        }
        
        offset = (nj1-nj0)%4;
        
 		if(offset!=0)
        {
            /* Epilogue loop to take care of non-4-multiple j particles */
            if(offset==1)
            {
                jnrA        = jjnr[k];   
                j3A         = 3*jnrA;  
                tjA          = nti+2*type[jnrA];
                GMX_MM_LOAD_1RVEC_1POINTER_PS(pos+j3A,jx,jy,jz);
                GMX_MM_LOAD_1VALUE_PS(charge+jnrA,jq);
                GMX_MM_LOAD_1VALUE_PS(invsqrta+jnrA,isaj);
                GMX_MM_LOAD_1PAIR_PS(vdwparam+tjA,c6,c12);
                mask        = mask1;
            } 
            else if(offset==2)
            {
                jnrA        = jjnr[k];   
                jnrB        = jjnr[k+1];
                j3A         = 3*jnrA;  
                j3B         = 3*jnrB;
                tjA          = nti+2*type[jnrA];
                tjB          = nti+2*type[jnrB];
                GMX_MM_LOAD_1RVEC_2POINTERS_PS(pos+j3A,pos+j3B,jx,jy,jz);
                GMX_MM_LOAD_2VALUES_PS(charge+jnrA,charge+jnrB,jq);
                GMX_MM_LOAD_2VALUES_PS(invsqrta+jnrA,invsqrta+jnrB,isaj);
                GMX_MM_LOAD_2PAIRS_PS(vdwparam+tjA,vdwparam+tjB,c6,c12);
                mask        = mask2;
            }
            else
            {
                /* offset must be 3 */   
                jnrA        = jjnr[k];   
                jnrB        = jjnr[k+1];
                jnrC        = jjnr[k+2];
                j3A         = 3*jnrA;  
                j3B         = 3*jnrB;
                j3C         = 3*jnrC;
                tjA          = nti+2*type[jnrA];
                tjB          = nti+2*type[jnrB];
                tjC          = nti+2*type[jnrC];
                GMX_MM_LOAD_1RVEC_3POINTERS_PS(pos+j3A,pos+j3B,pos+j3C,jx,jy,jz);
                GMX_MM_LOAD_3VALUES_PS(charge+jnrA,charge+jnrB,charge+jnrC,jq);
                GMX_MM_LOAD_3VALUES_PS(invsqrta+jnrA,invsqrta+jnrB,invsqrta+jnrC,isaj);
                GMX_MM_LOAD_3PAIRS_PS(vdwparam+tjA,vdwparam+tjB,vdwparam+tjC,c6,c12);
                mask        = mask3;
            }
            
            dx           = _mm_sub_ps(ix,jx);
			dy           = _mm_sub_ps(iy,jy);
			dz           = _mm_sub_ps(iz,jz);
			
			rsq          = gmx_mm_calc_rsq_ps(dx,dy,dz);
			
			/* invsqrt */								
			rinv         = gmx_mm_invsqrt_ps(rsq);
            rinv         = _mm_and_ps(rinv,mask);
            
 			rinvsq       = _mm_mul_ps(rinv,rinv);
            
			/***********************************/
			/* INTERACTION SECTION STARTS HERE */
			/***********************************/
            isaprod      = _mm_mul_ps(isai,isaj);
			qq           = _mm_and_ps(_mm_mul_ps(iq,jq),mask);        
			vcoul        = _mm_mul_ps(qq,rinv);
			fscal        = _mm_mul_ps(vcoul,rinv);                                 
            vctot        = _mm_add_ps(vctot,vcoul);
            
            /* Polarization interaction */
			qq           = _mm_mul_ps(qq,_mm_mul_ps(isaprod,gbfactor));
			gbscale      = _mm_mul_ps(isaprod,gbtabscale);
            
 			/* Calculate GB table index */
			r            = _mm_mul_ps(rsq,rinv);
			rtab         = _mm_mul_ps(r,gbscale);
			
			n0           = _mm_cvttps_epi32(rtab);
            eps          = _mm_sub_ps(rtab , _mm_cvtepi32_ps(n0) );
			nnn          = _mm_slli_epi32(n0,2);
			
            /* the tables are 16-byte aligned, so we can use _mm_load_ps */			
            Y            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,0)); /* YFGH */
            F            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,1)); /* YFGH */
            G            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,2)); /* YFGH */
            H            = _mm_load_ps(GBtab + gmx_mm_extract_epi32(nnn,3)); /* YFGH */
            _MM_TRANSPOSE4_PS(Y,F,G,H);
            
            G       = _mm_mul_ps(G,eps);
            H       = _mm_mul_ps(H, _mm_mul_ps(eps,eps) );
            F       = _mm_add_ps(F, _mm_add_ps( G , H ) );
            Y       = _mm_add_ps(Y, _mm_mul_ps(F, eps));
            F       = _mm_add_ps(F, _mm_add_ps(G , _mm_mul_ps(H,two)));
            vgb     = _mm_mul_ps(Y, qq);           
            fijGB   = _mm_mul_ps(F, _mm_mul_ps(qq,gbscale));
            
            dvdatmp = _mm_mul_ps(_mm_add_ps(vgb, _mm_mul_ps(fijGB,r)) , minushalf);            
            vgbtot  = _mm_add_ps(vgbtot, vgb);
            
            dvdasum = _mm_add_ps(dvdasum, dvdatmp);
            dvdatmp = _mm_mul_ps(dvdatmp, _mm_mul_ps(isaj,isaj));
            
            /* Lennard-Jones */
            
			rinvsix      = _mm_mul_ps(rinvsq,rinvsq);
			rinvsix      = _mm_mul_ps(rinvsix,rinvsq);
			
			vvdw6        = _mm_mul_ps(c6,rinvsix);
			vvdw12       = _mm_mul_ps(c12, _mm_mul_ps(rinvsix,rinvsix));
			vvdwtot      = _mm_add_ps(vvdwtot,_mm_sub_ps(vvdw12,vvdw6));
            
            fscal        = _mm_sub_ps(_mm_mul_ps(rinvsq, 
                                                 _mm_sub_ps(_mm_mul_ps(twelve,vvdw12),
                                                            _mm_mul_ps(six,vvdw6))),
                                      _mm_mul_ps( _mm_sub_ps( fijGB,fscal),rinv ));
            
            /***********************************/
			/*  INTERACTION SECTION ENDS HERE  */
			/***********************************/
            
            /* Calculate temporary vectorial force */
            tx           = _mm_mul_ps(fscal,dx);
            ty           = _mm_mul_ps(fscal,dy);
            tz           = _mm_mul_ps(fscal,dz);
            
            /* Increment i atom force */
            fix          = _mm_add_ps(fix,tx);
            fiy          = _mm_add_ps(fiy,ty);
            fiz          = _mm_add_ps(fiz,tz);
            
            if(offset==1)
            {
                GMX_MM_INCREMENT_1VALUE_PS(dvda+jnrA,dvdatmp);
                GMX_MM_DECREMENT_1RVEC_1POINTER_PS(faction+j3A,tx,ty,tz);
            } 
            else if(offset==2)
            {
                GMX_MM_INCREMENT_2VALUES_PS(dvda+jnrA,dvda+jnrB,dvdatmp);
                GMX_MM_DECREMENT_1RVEC_2POINTERS_PS(faction+j3A,faction+j3B,tx,ty,tz);
            }
            else
            {
                /* offset must be 3 */
                GMX_MM_INCREMENT_3VALUES_PS(dvda+jnrA,dvda+jnrB,dvda+jnrC,dvdatmp);
                GMX_MM_DECREMENT_1RVEC_3POINTERS_PS(faction+j3A,faction+j3B,faction+j3C,tx,ty,tz);
            }
        }
        
        dvdasum = _mm_mul_ps(dvdasum, _mm_mul_ps(isai,isai));
		gmx_mm_update_iforce_1atom_ps(&fix,&fiy,&fiz,faction+ii3,fshift+is3);
		
        ggid             = gid[n];         
	GMX_MM_UPDATE_4POT_PS(vctot,vc+ggid,vvdwtot,vvdw+ggid,vgbtot,gpol+ggid,dvdasum,dvda+ii);
    }
	
	*outeriter       = nri;            
    *inneriter       = nj1;            
}



/*
 * Gromacs nonbonded kernel nb_kernel410nf
 * Coulomb interaction:     Generalized-Born
 * VdW interaction:         Lennard-Jones
 * water optimization:      No
 * Calculate forces:        no
 */
void nb_kernel410nf_x86_64_sse(
                    int *           p_nri,
                    int *           iinr,
                    int *           jindex,
                    int *           jjnr,
                    int *           shift,
                    float *         shiftvec,
                    float *         fshift,
                    int *           gid,
                    float *         pos,
                    float *         faction,
                    float *         charge,
                    float *         p_facel,
                    float *         p_krf,
                    float *         p_crf,
                    float *         Vc,
                    int *           type,
                    int *           p_ntype,
                    float *         vdwparam,
                    float *         vvdw,
                    float *         p_tabscale,
                    float *         VFtab,
                    float *         invsqrta,
                    float *         dvda,
                    float *         p_gbtabscale,
                    float *         GBtab,
                    int *           p_nthreads,
                    int *           count,
                    void *          mtx,
                    int *           outeriter,
                    int *           inneriter,
                    float *         work)
{
    int           nri,ntype,nthreads;
    float         facel,krf,crf,tabscale,gbtabscale,vgb,fgb;
    int           n,ii,is3,ii3,k,nj0,nj1,jnr,j3,ggid;
    float         shX,shY,shZ;
    float         rinvsq;
    float         iq;
    float         qq,vcoul,vctot;
    int           nti;
    int           tj;
    float         rinvsix;
    float         vvdw6,vvdwtot;
    float         vvdw12;
    float         r,rt,eps,eps2;
    int           n0,nnn;
    float         Y,F,Geps,Heps2,Fp,VV;
    float         isai,isaj,isaprod,gbscale;
    float         ix1,iy1,iz1;
    float         jx1,jy1,jz1;
    float         dx11,dy11,dz11,rsq11,rinv11;
    float         c6,c12;
    const int     fractshift = 12;
    const int     fractmask = 8388607;
    const int     expshift = 23;
    const int     expmask = 2139095040;
    const int     explsb = 8388608;
    float         lu;
    int           iexp,addr;
    union { unsigned int bval; float fval; } bitpattern,result;

    nri              = *p_nri;         
    ntype            = *p_ntype;       
    nthreads         = *p_nthreads;    
    facel            = *p_facel;       
    krf              = *p_krf;         
    crf              = *p_crf;         
    tabscale         = *p_tabscale;    
    gbtabscale       = *p_gbtabscale;  
    nj1              = 0;              

    for(n=0; (n<nri); n++)
    {
        is3              = 3*shift[n];     
        shX              = shiftvec[is3];  
        shY              = shiftvec[is3+1];
        shZ              = shiftvec[is3+2];
        nj0              = jindex[n];      
        nj1              = jindex[n+1];    
        ii               = iinr[n];        
        ii3              = 3*ii;           
        ix1              = shX + pos[ii3+0];
        iy1              = shY + pos[ii3+1];
        iz1              = shZ + pos[ii3+2];
        iq               = facel*charge[ii];
        isai             = invsqrta[ii];   
        nti              = 2*ntype*type[ii];
        vctot            = 0;              
        vvdwtot          = 0;              
        
        for(k=nj0; (k<nj1); k++)
        {
            jnr              = jjnr[k];        
            j3               = 3*jnr;          
            jx1              = pos[j3+0];      
            jy1              = pos[j3+1];      
            jz1              = pos[j3+2];      
            dx11             = ix1 - jx1;      
            dy11             = iy1 - jy1;      
            dz11             = iz1 - jz1;      
            rsq11            = dx11*dx11+dy11*dy11+dz11*dz11;
            bitpattern.fval  = rsq11;          
            iexp             = (((bitpattern.bval)&expmask)>>expshift);
            addr             = (((bitpattern.bval)&(fractmask|explsb))>>fractshift);
            result.bval      = gmx_invsqrt_exptab[iexp] | gmx_invsqrt_fracttab[addr];
            lu               = result.fval;    
            rinv11           = (0.5*lu*(3.0-((rsq11*lu)*lu)));
            isaj             = invsqrta[jnr];  
            isaprod          = isai*isaj;      
            qq               = iq*charge[jnr]; 
            vcoul            = qq*rinv11;      
            qq               = isaprod*(-qq);  
            gbscale          = isaprod*gbtabscale;
            tj               = nti+2*type[jnr];
            c6               = vdwparam[tj];   
            c12              = vdwparam[tj+1]; 
            rinvsq           = rinv11*rinv11;  
            r                = rsq11*rinv11;   
            rt               = r*gbscale;      
            n0               = rt;             
            eps              = rt-n0;          
            eps2             = eps*eps;        
            nnn              = 4*n0;           
            Y                = GBtab[nnn];     
            F                = GBtab[nnn+1];   
            Geps             = eps*GBtab[nnn+2];
            Heps2            = eps2*GBtab[nnn+3];
            Fp               = F+Geps+Heps2;   
            VV               = Y+eps*Fp;       
            vgb              = qq*VV;          
            vctot            = vctot + vcoul;  
            rinvsix          = rinvsq*rinvsq*rinvsq;
            vvdw6            = c6*rinvsix;     
            vvdw12           = c12*rinvsix*rinvsix;
            vvdwtot          = vvdwtot+vvdw12-vvdw6;
        }
        
        ggid             = gid[n];         
        Vc[ggid]         = Vc[ggid] + vctot;
        vvdw[ggid]       = vvdw[ggid] + vvdwtot;
    }
    
    *outeriter       = nri;            
    *inneriter       = nj1;            
}


