/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*- 
 *
 * 
 * This file is part of Gromacs        Copyright (c) 1991-2004
 * David van der Spoel, Erik Lindahl, University of Groningen.
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * Gnomes, ROck Monsters And Chili Sauce
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* Must come directly after config.h */
#ifdef GMX_THREAD_SHM_FDECOMP
#include <thread_mpi.h>
#endif


#include "ppc_altivec_util.h"
#include "nb_kernel332_ppc_altivec.h"



void 
nb_kernel332_ppc_altivec  (int *             p_nri,
                       int               iinr[],
                       int               jindex[],
                       int               jjnr[],
                       int               shift[],
                       float             shiftvec[],
                       float             fshift[],
                       int               gid[],
                       float             pos[],
                       float             faction[],
                       float             charge[],
                       float *           p_facel,
                       float *           p_krf,
                       float *           p_crf,
                       float             Vc[],
                       int               type[],
                       int *             p_ntype,
                       float             vdwparam[],
                       float             Vvdw[],
                       float *           p_tabscale,
                       float             VFtab[],
                       float             invsqrta[],
                       float             dvda[],
                       float *           p_gbtabscale,
                       float             GBtab[],
                       int *             p_nthreads,
                       int *             count,
                       void *            mtx,
                       int *             outeriter,
                       int *             inneriter,
					   float *           work)
{
	vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
	vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

	vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
	vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
	vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

	vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
	vector float r11,r12,r13,r21,r22,r23,r31,r32,r33;
	vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23;
	vector float rinv31,rinv32,rinv33;

	vector float vfacel,nul;
	vector float fs11,fs12,fs13,fs21,fs22,fs23,fs31,fs32,fs33;
	vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3;
	vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3;
	vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12;
	vector float Vvdwtot,tsc,qqOOt,qqOHt,qqHHt,c6t,c12t;
	vector float VV11c,FF11c,VV12c,FF12c,VV13c,FF13c;
	vector float VV21c,FF21c,VV22c,FF22c,VV23c,FF23c;
	vector float VV31c,FF31c,VV32c,FF32c,VV33c,FF33c;
	vector float VVd,FFd,VVr,FFr;

	int n,k,ii,is3,ii3,nj0,nj1;
	int jnra,jnrb,jnrc,jnrd,tp,tj;
	int j3a,j3b,j3c,j3d;
	int nri, ntype, nouter, ninner;
#ifdef GMX_THREAD_SHM_FDECOMP
	int nn0, nn1;
#endif

    nouter   = 0;
    ninner   = 0;
    nri      = *p_nri;
    ntype    = *p_ntype;
	nul=vec_zero();
	vfacel=load_float_and_splat(p_facel);
	tsc=load_float_and_splat(p_tabscale);
	ii        = iinr[0];
	qO        = load_float_and_splat(charge+ii);
	qH        = load_float_and_splat(charge+ii+1);
	qqOO      = vec_madd(qO,qO,nul);
	qqOH      = vec_madd(qO,qH,nul);
	qqHH      = vec_madd(qH,qH,nul);
	qqOO      = vec_madd(qqOO,vfacel,nul);
	qqOH      = vec_madd(qqOH,vfacel,nul);
	qqHH      = vec_madd(qqHH,vfacel,nul);
	tp        = 2*type[ii];
	tj        = (ntype+1)*tp;
	load_1_pair(vdwparam+tj,&c6,&c12);
	c6        = vec_splat(c6,0);
	c12       = vec_splat(c12,0);

#ifdef GMX_THREAD_SHM_FDECOMP
    nthreads = *p_nthreads;
	do {
		tMPI_Thread_mutex_lock((tMPI_Thread_mutex_t *)mtx);
		nn0              = *count;
		nn1              = nn0+(nri-nn0)/(2*nthreads)+3;
		*count           = nn1;
		tMPI_Thread_mutex_unlock((tMPI_Thread_mutex_t *)mtx);
		if(nn1>nri) nn1=nri;
		for(n=nn0; (n<nn1); n++) {
#if 0
		} /* maintain correct indentation even with conditional left braces */
#endif
#else /* without tMPI_Threads */
		for(n=0;n<nri;n++) {
#endif  
			is3        = 3*shift[n];
			ii         = iinr[n];
			ii3        = 3*ii;
			load_1_3atoms_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
										  &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
			vctot      = nul;
			Vvdwtot     = nul;
			fix1       = nul;
			fiy1       = nul;
			fiz1       = nul;
			fix2       = nul;
			fiy2       = nul;
			fiz2       = nul;
			fix3       = nul;
			fiy3       = nul;
			fiz3       = nul;
			nj0        = jindex[n];
			nj1        = jindex[n+1];
    
			for(k=nj0; k<(nj1-3); k+=4) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				jnrd            = jjnr[k+3];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				j3d             = 3*jnrd;
				load_4_3atoms(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
				dx11            = vec_sub(ix1,jx1);
				dx12            = vec_sub(ix1,jx2);
				dx13            = vec_sub(ix1,jx3);
				dy11            = vec_sub(iy1,jy1);
				dy12            = vec_sub(iy1,jy2);
				dy13            = vec_sub(iy1,jy3);
				dz11            = vec_sub(iz1,jz1);
				dz12            = vec_sub(iz1,jz2);
				dz13            = vec_sub(iz1,jz3);
				dx21            = vec_sub(ix2,jx1);
				dx22            = vec_sub(ix2,jx2);
				dx23            = vec_sub(ix2,jx3);
				dy21            = vec_sub(iy2,jy1);
				dy22            = vec_sub(iy2,jy2);
				dy23            = vec_sub(iy2,jy3);
				dz21            = vec_sub(iz2,jz1);
				dz22            = vec_sub(iz2,jz2);
				dz23            = vec_sub(iz2,jz3);
				dx31            = vec_sub(ix3,jx1);
				dx32            = vec_sub(ix3,jx2);
				dx33            = vec_sub(ix3,jx3);
				dy31            = vec_sub(iy3,jy1);
				dy32            = vec_sub(iy3,jy2);
				dy33            = vec_sub(iy3,jy3);
				dz31            = vec_sub(iz3,jz1);
				dz32            = vec_sub(iz3,jz2);
				dz33            = vec_sub(iz3,jz3);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq12           = vec_madd(dx12,dx12,nul);
				rsq13           = vec_madd(dx13,dx13,nul);
				rsq21           = vec_madd(dx21,dx21,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq31           = vec_madd(dx31,dx31,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq12           = vec_madd(dy12,dy12,rsq12);
				rsq13           = vec_madd(dy13,dy13,rsq13);
				rsq21           = vec_madd(dy21,dy21,rsq21);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq31           = vec_madd(dy31,dy31,rsq31);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq12           = vec_madd(dz12,dz12,rsq12);
				rsq13           = vec_madd(dz13,dz13,rsq13);
				rsq21           = vec_madd(dz21,dz21,rsq21);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq31           = vec_madd(dz31,dz31,rsq31);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				r11             = vec_madd(rsq11,rinv11,nul); 
				r12             = vec_madd(rsq12,rinv12,nul); 
				r13             = vec_madd(rsq13,rinv13,nul); 
				r21             = vec_madd(rsq21,rinv21,nul); 
				r22             = vec_madd(rsq22,rinv22,nul); 
				r23             = vec_madd(rsq23,rinv23,nul); 
				r31             = vec_madd(rsq31,rinv31,nul); 
				r32             = vec_madd(rsq32,rinv32,nul); 
				r33             = vec_madd(rsq33,rinv33,nul); 

				do_4_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
										  &VV11c,&FF11c,&VVd,&FFd,&VVr,&FFr);
				do_4_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
				do_4_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
				do_4_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
				do_4_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
				do_4_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
				do_4_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
				do_4_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
				do_4_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);

				fs11            = vec_nmsub(qqOO,FF11c,nul);
				fs12            = vec_nmsub(qqOH,FF12c,nul);
				fs13            = vec_nmsub(qqOH,FF13c,nul);
				fs21            = vec_nmsub(qqOH,FF21c,nul);
				fs11            = vec_nmsub(c6,FFd,fs11);
				fs22            = vec_nmsub(qqHH,FF22c,nul);
				fs23            = vec_nmsub(qqHH,FF23c,nul);
				fs31            = vec_nmsub(qqOH,FF31c,nul);
				fs32            = vec_nmsub(qqHH,FF32c,nul);
				fs33            = vec_nmsub(qqHH,FF33c,nul);
				fs11            = vec_nmsub(c12,FFr,fs11);

				vctot           = vec_madd(qqOO,VV11c,vctot);
				vctot           = vec_madd(qqOH,VV12c,vctot);
				vctot           = vec_madd(qqOH,VV13c,vctot);
				vctot           = vec_madd(qqOH,VV21c,vctot);
				vctot           = vec_madd(qqHH,VV22c,vctot);
				vctot           = vec_madd(qqHH,VV23c,vctot);
				vctot           = vec_madd(qqOH,VV31c,vctot);
				vctot           = vec_madd(qqHH,VV32c,vctot);
				vctot           = vec_madd(qqHH,VV33c,vctot);

				fs11            = vec_madd(fs11,tsc,nul);
				fs12            = vec_madd(fs12,tsc,nul);
				fs13            = vec_madd(fs13,tsc,nul);
				fs21            = vec_madd(fs21,tsc,nul);
				fs22            = vec_madd(fs22,tsc,nul);
				fs23            = vec_madd(fs23,tsc,nul);
				fs31            = vec_madd(fs31,tsc,nul);
				fs32            = vec_madd(fs32,tsc,nul);
				fs33            = vec_madd(fs33,tsc,nul);
      
				fs11            = vec_madd(fs11,rinv11,nul);
				fs12            = vec_madd(fs12,rinv12,nul);
				fs13            = vec_madd(fs13,rinv13,nul);
				fs21            = vec_madd(fs21,rinv21,nul);
				fs22            = vec_madd(fs22,rinv22,nul);
				fs23            = vec_madd(fs23,rinv23,nul);
				fs31            = vec_madd(fs31,rinv31,nul);
				fs32            = vec_madd(fs32,rinv32,nul);
				fs33            = vec_madd(fs33,rinv33,nul);      

				fix1            = vec_madd(fs11,dx11,fix1);
				fiy1            = vec_madd(fs11,dy11,fiy1);
				fiz1            = vec_madd(fs11,dz11,fiz1);
				fix2            = vec_madd(fs21,dx21,fix2);
				fiy2            = vec_madd(fs21,dy21,fiy2);
				fiz2            = vec_madd(fs21,dz21,fiz2);
				fix3            = vec_madd(fs31,dx31,fix3);
				fiy3            = vec_madd(fs31,dy31,fiy3);
				fiz3            = vec_madd(fs31,dz31,fiz3);

				fix1            = vec_madd(fs12,dx12,fix1);
				fiy1            = vec_madd(fs12,dy12,fiy1);
				fiz1            = vec_madd(fs12,dz12,fiz1);
				fix2            = vec_madd(fs22,dx22,fix2);
				fiy2            = vec_madd(fs22,dy22,fiy2);
				fiz2            = vec_madd(fs22,dz22,fiz2);
				fix3            = vec_madd(fs32,dx32,fix3);
				fiy3            = vec_madd(fs32,dy32,fiy3);
				fiz3            = vec_madd(fs32,dz32,fiz3);

				fix1            = vec_madd(fs13,dx13,fix1);
				fiy1            = vec_madd(fs13,dy13,fiy1);
				fiz1            = vec_madd(fs13,dz13,fiz1);
				fix2            = vec_madd(fs23,dx23,fix2);
				fiy2            = vec_madd(fs23,dy23,fiy2);
				fiz2            = vec_madd(fs23,dz23,fiz2);
				fix3            = vec_madd(fs33,dx33,fix3);
				fiy3            = vec_madd(fs33,dy33,fiy3);
				fiz3            = vec_madd(fs33,dz33,fiz3);

				fjx1            = vec_nmsub(fs11,dx11,nul);
				fjy1            = vec_nmsub(fs11,dy11,nul);
				fjz1            = vec_nmsub(fs11,dz11,nul);
				fjx2            = vec_nmsub(fs12,dx12,nul);
				fjy2            = vec_nmsub(fs12,dy12,nul);
				fjz2            = vec_nmsub(fs12,dz12,nul);
				fjx3            = vec_nmsub(fs13,dx13,nul);
				fjy3            = vec_nmsub(fs13,dy13,nul);
				fjz3            = vec_nmsub(fs13,dz13,nul);

				fjx1            = vec_nmsub(fs21,dx21,fjx1);
				fjy1            = vec_nmsub(fs21,dy21,fjy1);
				fjz1            = vec_nmsub(fs21,dz21,fjz1);
				fjx2            = vec_nmsub(fs22,dx22,fjx2);
				fjy2            = vec_nmsub(fs22,dy22,fjy2);
				fjz2            = vec_nmsub(fs22,dz22,fjz2);
				fjx3            = vec_nmsub(fs23,dx23,fjx3);
				fjy3            = vec_nmsub(fs23,dy23,fjy3);
				fjz3            = vec_nmsub(fs23,dz23,fjz3);

				fjx1            = vec_nmsub(fs31,dx31,fjx1);
				fjy1            = vec_nmsub(fs31,dy31,fjy1);
				fjz1            = vec_nmsub(fs31,dz31,fjz1);
				fjx2            = vec_nmsub(fs32,dx32,fjx2);
				fjy2            = vec_nmsub(fs32,dy32,fjy2);
				fjz2            = vec_nmsub(fs32,dz32,fjz2);
				fjx3            = vec_nmsub(fs33,dx33,fjx3);
				fjy3            = vec_nmsub(fs33,dy33,fjy3);
				fjz3            = vec_nmsub(fs33,dz33,fjz3);
     
				add_force_to_4_3atoms(faction+j3a,faction+j3b,
									  faction+j3c,faction+j3d,
									  fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,
									  fjx3,fjy3,fjz3);
			} 
			if(k<(nj1-2)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				load_3_3atoms(pos+j3a,pos+j3b,pos+j3c,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
				qqOOt           = vec_sld(qqOO,nul,4);
				qqOHt           = vec_sld(qqOH,nul,4);
				qqHHt           = vec_sld(qqHH,nul,4);
				c6t             = vec_sld(c6,nul,4);
				c12t            = vec_sld(c12,nul,4);

				dx11            = vec_sub(ix1,jx1);
				dx12            = vec_sub(ix1,jx2);
				dx13            = vec_sub(ix1,jx3);
				dy11            = vec_sub(iy1,jy1);
				dy12            = vec_sub(iy1,jy2);
				dy13            = vec_sub(iy1,jy3);
				dz11            = vec_sub(iz1,jz1);
				dz12            = vec_sub(iz1,jz2);
				dz13            = vec_sub(iz1,jz3);
				dx21            = vec_sub(ix2,jx1);
				dx22            = vec_sub(ix2,jx2);
				dx23            = vec_sub(ix2,jx3);
				dy21            = vec_sub(iy2,jy1);
				dy22            = vec_sub(iy2,jy2);
				dy23            = vec_sub(iy2,jy3);
				dz21            = vec_sub(iz2,jz1);
				dz22            = vec_sub(iz2,jz2);
				dz23            = vec_sub(iz2,jz3);
				dx31            = vec_sub(ix3,jx1);
				dx32            = vec_sub(ix3,jx2);
				dx33            = vec_sub(ix3,jx3);
				dy31            = vec_sub(iy3,jy1);
				dy32            = vec_sub(iy3,jy2);
				dy33            = vec_sub(iy3,jy3);
				dz31            = vec_sub(iz3,jz1);
				dz32            = vec_sub(iz3,jz2);
				dz33            = vec_sub(iz3,jz3);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq12           = vec_madd(dx12,dx12,nul);
				rsq13           = vec_madd(dx13,dx13,nul);
				rsq21           = vec_madd(dx21,dx21,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq31           = vec_madd(dx31,dx31,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq12           = vec_madd(dy12,dy12,rsq12);
				rsq13           = vec_madd(dy13,dy13,rsq13);
				rsq21           = vec_madd(dy21,dy21,rsq21);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq31           = vec_madd(dy31,dy31,rsq31);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq12           = vec_madd(dz12,dz12,rsq12);
				rsq13           = vec_madd(dz13,dz13,rsq13);
				rsq21           = vec_madd(dz21,dz21,rsq21);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq31           = vec_madd(dz31,dz31,rsq31);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);

				zero_highest_element_in_9_vectors(&rsq11,&rsq12,&rsq13,
												  &rsq21,&rsq22,&rsq23,
												  &rsq31,&rsq32,&rsq33);
      
				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);

				zero_highest_element_in_9_vectors(&rinv11,&rinv12,&rinv13,
												  &rinv21,&rinv22,&rinv23,
												  &rinv31,&rinv32,&rinv33);

				r11             = vec_madd(rsq11,rinv11,nul); 
				r12             = vec_madd(rsq12,rinv12,nul); 
				r13             = vec_madd(rsq13,rinv13,nul); 
				r21             = vec_madd(rsq21,rinv21,nul); 
				r22             = vec_madd(rsq22,rinv22,nul); 
				r23             = vec_madd(rsq23,rinv23,nul); 
				r31             = vec_madd(rsq31,rinv31,nul); 
				r32             = vec_madd(rsq32,rinv32,nul); 
				r33             = vec_madd(rsq33,rinv33,nul); 

				do_3_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
										  &VV11c,&FF11c,&VVd,&FFd,&VVr,&FFr);
				do_3_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
				do_3_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
				do_3_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
				do_3_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
				do_3_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
				do_3_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
				do_3_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
				do_3_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

				Vvdwtot          = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot          = vec_madd(c12t,VVr,Vvdwtot);

				fs11            = vec_nmsub(qqOOt,FF11c,nul);
				fs12            = vec_nmsub(qqOHt,FF12c,nul);
				fs13            = vec_nmsub(qqOHt,FF13c,nul);
				fs21            = vec_nmsub(qqOHt,FF21c,nul);
				fs11            = vec_nmsub(c6t,FFd,fs11);
				fs22            = vec_nmsub(qqHHt,FF22c,nul);
				fs23            = vec_nmsub(qqHHt,FF23c,nul);
				fs31            = vec_nmsub(qqOHt,FF31c,nul);
				fs32            = vec_nmsub(qqHHt,FF32c,nul);
				fs33            = vec_nmsub(qqHHt,FF33c,nul);
				fs11            = vec_nmsub(c12t,FFr,fs11);

				vctot           = vec_madd(qqOOt,VV11c,vctot);
				vctot           = vec_madd(qqOHt,VV12c,vctot);
				vctot           = vec_madd(qqOHt,VV13c,vctot);
				vctot           = vec_madd(qqOHt,VV21c,vctot);
				vctot           = vec_madd(qqHHt,VV22c,vctot);
				vctot           = vec_madd(qqHHt,VV23c,vctot);
				vctot           = vec_madd(qqOHt,VV31c,vctot);
				vctot           = vec_madd(qqHHt,VV32c,vctot);
				vctot           = vec_madd(qqHHt,VV33c,vctot);

				fs11            = vec_madd(fs11,tsc,nul);
				fs12            = vec_madd(fs12,tsc,nul);
				fs13            = vec_madd(fs13,tsc,nul);
				fs21            = vec_madd(fs21,tsc,nul);
				fs22            = vec_madd(fs22,tsc,nul);
				fs23            = vec_madd(fs23,tsc,nul);
				fs31            = vec_madd(fs31,tsc,nul);
				fs32            = vec_madd(fs32,tsc,nul);
				fs33            = vec_madd(fs33,tsc,nul);
      
				fs11            = vec_madd(fs11,rinv11,nul);
				fs12            = vec_madd(fs12,rinv12,nul);
				fs13            = vec_madd(fs13,rinv13,nul);
				fs21            = vec_madd(fs21,rinv21,nul);
				fs22            = vec_madd(fs22,rinv22,nul);
				fs23            = vec_madd(fs23,rinv23,nul);
				fs31            = vec_madd(fs31,rinv31,nul);
				fs32            = vec_madd(fs32,rinv32,nul);
				fs33            = vec_madd(fs33,rinv33,nul);            

				fix1            = vec_madd(fs11,dx11,fix1);
				fiy1            = vec_madd(fs11,dy11,fiy1);
				fiz1            = vec_madd(fs11,dz11,fiz1);
				fix2            = vec_madd(fs21,dx21,fix2);
				fiy2            = vec_madd(fs21,dy21,fiy2);
				fiz2            = vec_madd(fs21,dz21,fiz2);
				fix3            = vec_madd(fs31,dx31,fix3);
				fiy3            = vec_madd(fs31,dy31,fiy3);
				fiz3            = vec_madd(fs31,dz31,fiz3);

				fix1            = vec_madd(fs12,dx12,fix1);
				fiy1            = vec_madd(fs12,dy12,fiy1);
				fiz1            = vec_madd(fs12,dz12,fiz1);
				fix2            = vec_madd(fs22,dx22,fix2);
				fiy2            = vec_madd(fs22,dy22,fiy2);
				fiz2            = vec_madd(fs22,dz22,fiz2);
				fix3            = vec_madd(fs32,dx32,fix3);
				fiy3            = vec_madd(fs32,dy32,fiy3);
				fiz3            = vec_madd(fs32,dz32,fiz3);

				fix1            = vec_madd(fs13,dx13,fix1);
				fiy1            = vec_madd(fs13,dy13,fiy1);
				fiz1            = vec_madd(fs13,dz13,fiz1);
				fix2            = vec_madd(fs23,dx23,fix2);
				fiy2            = vec_madd(fs23,dy23,fiy2);
				fiz2            = vec_madd(fs23,dz23,fiz2);
				fix3            = vec_madd(fs33,dx33,fix3);
				fiy3            = vec_madd(fs33,dy33,fiy3);
				fiz3            = vec_madd(fs33,dz33,fiz3);

				fjx1            = vec_nmsub(fs11,dx11,nul);
				fjy1            = vec_nmsub(fs11,dy11,nul);
				fjz1            = vec_nmsub(fs11,dz11,nul);
				fjx2            = vec_nmsub(fs12,dx12,nul);
				fjy2            = vec_nmsub(fs12,dy12,nul);
				fjz2            = vec_nmsub(fs12,dz12,nul);
				fjx3            = vec_nmsub(fs13,dx13,nul);
				fjy3            = vec_nmsub(fs13,dy13,nul);
				fjz3            = vec_nmsub(fs13,dz13,nul);

				fjx1            = vec_nmsub(fs21,dx21,fjx1);
				fjy1            = vec_nmsub(fs21,dy21,fjy1);
				fjz1            = vec_nmsub(fs21,dz21,fjz1);
				fjx2            = vec_nmsub(fs22,dx22,fjx2);
				fjy2            = vec_nmsub(fs22,dy22,fjy2);
				fjz2            = vec_nmsub(fs22,dz22,fjz2);
				fjx3            = vec_nmsub(fs23,dx23,fjx3);
				fjy3            = vec_nmsub(fs23,dy23,fjy3);
				fjz3            = vec_nmsub(fs23,dz23,fjz3);

				fjx1            = vec_nmsub(fs31,dx31,fjx1);
				fjy1            = vec_nmsub(fs31,dy31,fjy1);
				fjz1            = vec_nmsub(fs31,dz31,fjz1);
				fjx2            = vec_nmsub(fs32,dx32,fjx2);
				fjy2            = vec_nmsub(fs32,dy32,fjy2);
				fjz2            = vec_nmsub(fs32,dz32,fjz2);
				fjx3            = vec_nmsub(fs33,dx33,fjx3);
				fjy3            = vec_nmsub(fs33,dy33,fjy3);
				fjz3            = vec_nmsub(fs33,dz33,fjz3);

				add_force_to_3_3atoms(faction+j3a,faction+j3b,faction+j3c,
									  fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,
									  fjx3,fjy3,fjz3);
			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				load_2_3atoms(pos+j3a,pos+j3b,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
				qqOOt           = vec_sld(qqOO,nul,8);
				qqOHt           = vec_sld(qqOH,nul,8);
				qqHHt           = vec_sld(qqHH,nul,8);
				c6t             = vec_sld(c6,nul,8);
				c12t            = vec_sld(c12,nul,8);

				dx11            = vec_sub(ix1,jx1);
				dx12            = vec_sub(ix1,jx2);
				dx13            = vec_sub(ix1,jx3);
				dy11            = vec_sub(iy1,jy1);
				dy12            = vec_sub(iy1,jy2);
				dy13            = vec_sub(iy1,jy3);
				dz11            = vec_sub(iz1,jz1);
				dz12            = vec_sub(iz1,jz2);
				dz13            = vec_sub(iz1,jz3);
				dx21            = vec_sub(ix2,jx1);
				dx22            = vec_sub(ix2,jx2);
				dx23            = vec_sub(ix2,jx3);
				dy21            = vec_sub(iy2,jy1);
				dy22            = vec_sub(iy2,jy2);
				dy23            = vec_sub(iy2,jy3);
				dz21            = vec_sub(iz2,jz1);
				dz22            = vec_sub(iz2,jz2);
				dz23            = vec_sub(iz2,jz3);
				dx31            = vec_sub(ix3,jx1);
				dx32            = vec_sub(ix3,jx2);
				dx33            = vec_sub(ix3,jx3);
				dy31            = vec_sub(iy3,jy1);
				dy32            = vec_sub(iy3,jy2);
				dy33            = vec_sub(iy3,jy3);
				dz31            = vec_sub(iz3,jz1);
				dz32            = vec_sub(iz3,jz2);
				dz33            = vec_sub(iz3,jz3);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq12           = vec_madd(dx12,dx12,nul);
				rsq13           = vec_madd(dx13,dx13,nul);
				rsq21           = vec_madd(dx21,dx21,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq31           = vec_madd(dx31,dx31,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq12           = vec_madd(dy12,dy12,rsq12);
				rsq13           = vec_madd(dy13,dy13,rsq13);
				rsq21           = vec_madd(dy21,dy21,rsq21);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq31           = vec_madd(dy31,dy31,rsq31);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq12           = vec_madd(dz12,dz12,rsq12);
				rsq13           = vec_madd(dz13,dz13,rsq13);
				rsq21           = vec_madd(dz21,dz21,rsq21);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq31           = vec_madd(dz31,dz31,rsq31);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);

				zero_highest_2_elements_in_9_vectors(&rsq11,&rsq12,&rsq13,
													 &rsq21,&rsq22,&rsq23,
													 &rsq31,&rsq32,&rsq33);

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				zero_highest_2_elements_in_9_vectors(&rinv11,&rinv12,&rinv13,
													 &rinv21,&rinv22,&rinv23,
													 &rinv31,&rinv32,&rinv33);

				r11             = vec_madd(rsq11,rinv11,nul); 
				r12             = vec_madd(rsq12,rinv12,nul); 
				r13             = vec_madd(rsq13,rinv13,nul); 
				r21             = vec_madd(rsq21,rinv21,nul); 
				r22             = vec_madd(rsq22,rinv22,nul); 
				r23             = vec_madd(rsq23,rinv23,nul); 
				r31             = vec_madd(rsq31,rinv31,nul); 
				r32             = vec_madd(rsq32,rinv32,nul); 
				r33             = vec_madd(rsq33,rinv33,nul); 

				do_2_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
										  &VV11c,&FF11c,&VVd,&FFd,&VVr,&FFr);
				do_2_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
				do_2_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
				do_2_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
				do_2_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
				do_2_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
				do_2_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
				do_2_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
				do_2_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

				Vvdwtot          = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot          = vec_madd(c12t,VVr,Vvdwtot);

				fs11            = vec_nmsub(qqOOt,FF11c,nul);
				fs12            = vec_nmsub(qqOHt,FF12c,nul);
				fs13            = vec_nmsub(qqOHt,FF13c,nul);
				fs21            = vec_nmsub(qqOHt,FF21c,nul);
				fs11            = vec_nmsub(c6t,FFd,fs11);
				fs22            = vec_nmsub(qqHHt,FF22c,nul);
				fs23            = vec_nmsub(qqHHt,FF23c,nul);
				fs31            = vec_nmsub(qqOHt,FF31c,nul);
				fs32            = vec_nmsub(qqHHt,FF32c,nul);
				fs33            = vec_nmsub(qqHHt,FF33c,nul);
				fs11            = vec_nmsub(c12t,FFr,fs11);

				vctot           = vec_madd(qqOOt,VV11c,vctot);
				vctot           = vec_madd(qqOHt,VV12c,vctot);
				vctot           = vec_madd(qqOHt,VV13c,vctot);
				vctot           = vec_madd(qqOHt,VV21c,vctot);
				vctot           = vec_madd(qqHHt,VV22c,vctot);
				vctot           = vec_madd(qqHHt,VV23c,vctot);
				vctot           = vec_madd(qqOHt,VV31c,vctot);
				vctot           = vec_madd(qqHHt,VV32c,vctot);
				vctot           = vec_madd(qqHHt,VV33c,vctot);

				fs11            = vec_madd(fs11,tsc,nul);
				fs12            = vec_madd(fs12,tsc,nul);
				fs13            = vec_madd(fs13,tsc,nul);
				fs21            = vec_madd(fs21,tsc,nul);
				fs22            = vec_madd(fs22,tsc,nul);
				fs23            = vec_madd(fs23,tsc,nul);
				fs31            = vec_madd(fs31,tsc,nul);
				fs32            = vec_madd(fs32,tsc,nul);
				fs33            = vec_madd(fs33,tsc,nul);
      
				fs11            = vec_madd(fs11,rinv11,nul);
				fs12            = vec_madd(fs12,rinv12,nul);
				fs13            = vec_madd(fs13,rinv13,nul);
				fs21            = vec_madd(fs21,rinv21,nul);
				fs22            = vec_madd(fs22,rinv22,nul);
				fs23            = vec_madd(fs23,rinv23,nul);
				fs31            = vec_madd(fs31,rinv31,nul);
				fs32            = vec_madd(fs32,rinv32,nul);
				fs33            = vec_madd(fs33,rinv33,nul);      

				fix1            = vec_madd(fs11,dx11,fix1);
				fiy1            = vec_madd(fs11,dy11,fiy1);
				fiz1            = vec_madd(fs11,dz11,fiz1);
				fix2            = vec_madd(fs21,dx21,fix2);
				fiy2            = vec_madd(fs21,dy21,fiy2);
				fiz2            = vec_madd(fs21,dz21,fiz2);
				fix3            = vec_madd(fs31,dx31,fix3);
				fiy3            = vec_madd(fs31,dy31,fiy3);
				fiz3            = vec_madd(fs31,dz31,fiz3);

				fix1            = vec_madd(fs12,dx12,fix1);
				fiy1            = vec_madd(fs12,dy12,fiy1);
				fiz1            = vec_madd(fs12,dz12,fiz1);
				fix2            = vec_madd(fs22,dx22,fix2);
				fiy2            = vec_madd(fs22,dy22,fiy2);
				fiz2            = vec_madd(fs22,dz22,fiz2);
				fix3            = vec_madd(fs32,dx32,fix3);
				fiy3            = vec_madd(fs32,dy32,fiy3);
				fiz3            = vec_madd(fs32,dz32,fiz3);

				fix1            = vec_madd(fs13,dx13,fix1);
				fiy1            = vec_madd(fs13,dy13,fiy1);
				fiz1            = vec_madd(fs13,dz13,fiz1);
				fix2            = vec_madd(fs23,dx23,fix2);
				fiy2            = vec_madd(fs23,dy23,fiy2);
				fiz2            = vec_madd(fs23,dz23,fiz2);
				fix3            = vec_madd(fs33,dx33,fix3);
				fiy3            = vec_madd(fs33,dy33,fiy3);
				fiz3            = vec_madd(fs33,dz33,fiz3);

				fjx1            = vec_nmsub(fs11,dx11,nul);
				fjy1            = vec_nmsub(fs11,dy11,nul);
				fjz1            = vec_nmsub(fs11,dz11,nul);
				fjx2            = vec_nmsub(fs12,dx12,nul);
				fjy2            = vec_nmsub(fs12,dy12,nul);
				fjz2            = vec_nmsub(fs12,dz12,nul);
				fjx3            = vec_nmsub(fs13,dx13,nul);
				fjy3            = vec_nmsub(fs13,dy13,nul);
				fjz3            = vec_nmsub(fs13,dz13,nul);

				fjx1            = vec_nmsub(fs21,dx21,fjx1);
				fjy1            = vec_nmsub(fs21,dy21,fjy1);
				fjz1            = vec_nmsub(fs21,dz21,fjz1);
				fjx2            = vec_nmsub(fs22,dx22,fjx2);
				fjy2            = vec_nmsub(fs22,dy22,fjy2);
				fjz2            = vec_nmsub(fs22,dz22,fjz2);
				fjx3            = vec_nmsub(fs23,dx23,fjx3);
				fjy3            = vec_nmsub(fs23,dy23,fjy3);
				fjz3            = vec_nmsub(fs23,dz23,fjz3);

				fjx1            = vec_nmsub(fs31,dx31,fjx1);
				fjy1            = vec_nmsub(fs31,dy31,fjy1);
				fjz1            = vec_nmsub(fs31,dz31,fjz1);
				fjx2            = vec_nmsub(fs32,dx32,fjx2);
				fjy2            = vec_nmsub(fs32,dy32,fjy2);
				fjz2            = vec_nmsub(fs32,dz32,fjz2);
				fjx3            = vec_nmsub(fs33,dx33,fjx3);
				fjy3            = vec_nmsub(fs33,dy33,fjy3);
				fjz3            = vec_nmsub(fs33,dz33,fjz3);
     
				add_force_to_2_3atoms(faction+j3a,faction+j3b,
									  fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,
									  fjx3,fjy3,fjz3);
			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				load_1_3atoms(pos+j3a,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
				qqOOt           = vec_sld(qqOO,nul,12);
				qqOHt           = vec_sld(qqOH,nul,12);
				qqHHt           = vec_sld(qqHH,nul,12);
				c6t             = vec_sld(c6,nul,12);
				c12t            = vec_sld(c12,nul,12);

				dx11            = vec_sub(ix1,jx1);
				dx12            = vec_sub(ix1,jx2);
				dx13            = vec_sub(ix1,jx3);
				dy11            = vec_sub(iy1,jy1);
				dy12            = vec_sub(iy1,jy2);
				dy13            = vec_sub(iy1,jy3);
				dz11            = vec_sub(iz1,jz1);
				dz12            = vec_sub(iz1,jz2);
				dz13            = vec_sub(iz1,jz3);
				dx21            = vec_sub(ix2,jx1);
				dx22            = vec_sub(ix2,jx2);
				dx23            = vec_sub(ix2,jx3);
				dy21            = vec_sub(iy2,jy1);
				dy22            = vec_sub(iy2,jy2);
				dy23            = vec_sub(iy2,jy3);
				dz21            = vec_sub(iz2,jz1);
				dz22            = vec_sub(iz2,jz2);
				dz23            = vec_sub(iz2,jz3);
				dx31            = vec_sub(ix3,jx1);
				dx32            = vec_sub(ix3,jx2);
				dx33            = vec_sub(ix3,jx3);
				dy31            = vec_sub(iy3,jy1);
				dy32            = vec_sub(iy3,jy2);
				dy33            = vec_sub(iy3,jy3);
				dz31            = vec_sub(iz3,jz1);
				dz32            = vec_sub(iz3,jz2);
				dz33            = vec_sub(iz3,jz3);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq12           = vec_madd(dx12,dx12,nul);
				rsq13           = vec_madd(dx13,dx13,nul);
				rsq21           = vec_madd(dx21,dx21,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq31           = vec_madd(dx31,dx31,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq12           = vec_madd(dy12,dy12,rsq12);
				rsq13           = vec_madd(dy13,dy13,rsq13);
				rsq21           = vec_madd(dy21,dy21,rsq21);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq31           = vec_madd(dy31,dy31,rsq31);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq12           = vec_madd(dz12,dz12,rsq12);
				rsq13           = vec_madd(dz13,dz13,rsq13);
				rsq21           = vec_madd(dz21,dz21,rsq21);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq31           = vec_madd(dz31,dz31,rsq31);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);

				zero_highest_3_elements_in_9_vectors(&rsq11,&rsq12,&rsq13,
													 &rsq21,&rsq22,&rsq23,
													 &rsq31,&rsq32,&rsq33);

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				zero_highest_3_elements_in_9_vectors(&rinv11,&rinv12,&rinv13,
													 &rinv21,&rinv22,&rinv23,
													 &rinv31,&rinv32,&rinv33);

				r11             = vec_madd(rsq11,rinv11,nul); 
				r12             = vec_madd(rsq12,rinv12,nul); 
				r13             = vec_madd(rsq13,rinv13,nul); 
				r21             = vec_madd(rsq21,rinv21,nul); 
				r22             = vec_madd(rsq22,rinv22,nul); 
				r23             = vec_madd(rsq23,rinv23,nul); 
				r31             = vec_madd(rsq31,rinv31,nul); 
				r32             = vec_madd(rsq32,rinv32,nul); 
				r33             = vec_madd(rsq33,rinv33,nul); 

				do_1_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
										  &VV11c,&FF11c,&VVd,&FFd,&VVr,&FFr);
				do_1_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c,&FF12c);
				do_1_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c,&FF13c);
				do_1_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c,&FF21c);
				do_1_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c,&FF22c);
				do_1_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c,&FF23c);
				do_1_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c,&FF31c);
				do_1_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c,&FF32c);
				do_1_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c,&FF33c);

				Vvdwtot          = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot          = vec_madd(c12t,VVr,Vvdwtot);

				fs11            = vec_nmsub(qqOOt,FF11c,nul);
				fs12            = vec_nmsub(qqOHt,FF12c,nul);
				fs13            = vec_nmsub(qqOHt,FF13c,nul);
				fs21            = vec_nmsub(qqOHt,FF21c,nul);
				fs11            = vec_nmsub(c6t,FFd,fs11);
				fs22            = vec_nmsub(qqHHt,FF22c,nul);
				fs23            = vec_nmsub(qqHHt,FF23c,nul);
				fs31            = vec_nmsub(qqOHt,FF31c,nul);
				fs32            = vec_nmsub(qqHHt,FF32c,nul);
				fs33            = vec_nmsub(qqHHt,FF33c,nul);
				fs11            = vec_nmsub(c12t,FFr,fs11);

				vctot           = vec_madd(qqOOt,VV11c,vctot);
				vctot           = vec_madd(qqOHt,VV12c,vctot);
				vctot           = vec_madd(qqOHt,VV13c,vctot);
				vctot           = vec_madd(qqOHt,VV21c,vctot);
				vctot           = vec_madd(qqHHt,VV22c,vctot);
				vctot           = vec_madd(qqHHt,VV23c,vctot);
				vctot           = vec_madd(qqOHt,VV31c,vctot);
				vctot           = vec_madd(qqHHt,VV32c,vctot);
				vctot           = vec_madd(qqHHt,VV33c,vctot);

				fs11            = vec_madd(fs11,tsc,nul);
				fs12            = vec_madd(fs12,tsc,nul);
				fs13            = vec_madd(fs13,tsc,nul);
				fs21            = vec_madd(fs21,tsc,nul);
				fs22            = vec_madd(fs22,tsc,nul);
				fs23            = vec_madd(fs23,tsc,nul);
				fs31            = vec_madd(fs31,tsc,nul);
				fs32            = vec_madd(fs32,tsc,nul);
				fs33            = vec_madd(fs33,tsc,nul);
      
				fs11            = vec_madd(fs11,rinv11,nul);
				fs12            = vec_madd(fs12,rinv12,nul);
				fs13            = vec_madd(fs13,rinv13,nul);
				fs21            = vec_madd(fs21,rinv21,nul);
				fs22            = vec_madd(fs22,rinv22,nul);
				fs23            = vec_madd(fs23,rinv23,nul);
				fs31            = vec_madd(fs31,rinv31,nul);
				fs32            = vec_madd(fs32,rinv32,nul);
				fs33            = vec_madd(fs33,rinv33,nul);      

				fix1            = vec_madd(fs11,dx11,fix1);
				fiy1            = vec_madd(fs11,dy11,fiy1);
				fiz1            = vec_madd(fs11,dz11,fiz1);
				fix2            = vec_madd(fs21,dx21,fix2);
				fiy2            = vec_madd(fs21,dy21,fiy2);
				fiz2            = vec_madd(fs21,dz21,fiz2);
				fix3            = vec_madd(fs31,dx31,fix3);
				fiy3            = vec_madd(fs31,dy31,fiy3);
				fiz3            = vec_madd(fs31,dz31,fiz3);

				fix1            = vec_madd(fs12,dx12,fix1);
				fiy1            = vec_madd(fs12,dy12,fiy1);
				fiz1            = vec_madd(fs12,dz12,fiz1);
				fix2            = vec_madd(fs22,dx22,fix2);
				fiy2            = vec_madd(fs22,dy22,fiy2);
				fiz2            = vec_madd(fs22,dz22,fiz2);
				fix3            = vec_madd(fs32,dx32,fix3);
				fiy3            = vec_madd(fs32,dy32,fiy3);
				fiz3            = vec_madd(fs32,dz32,fiz3);

				fix1            = vec_madd(fs13,dx13,fix1);
				fiy1            = vec_madd(fs13,dy13,fiy1);
				fiz1            = vec_madd(fs13,dz13,fiz1);
				fix2            = vec_madd(fs23,dx23,fix2);
				fiy2            = vec_madd(fs23,dy23,fiy2);
				fiz2            = vec_madd(fs23,dz23,fiz2);
				fix3            = vec_madd(fs33,dx33,fix3);
				fiy3            = vec_madd(fs33,dy33,fiy3);
				fiz3            = vec_madd(fs33,dz33,fiz3);

				fjx1            = vec_nmsub(fs11,dx11,nul);
				fjy1            = vec_nmsub(fs11,dy11,nul);
				fjz1            = vec_nmsub(fs11,dz11,nul);
				fjx2            = vec_nmsub(fs12,dx12,nul);
				fjy2            = vec_nmsub(fs12,dy12,nul);
				fjz2            = vec_nmsub(fs12,dz12,nul);
				fjx3            = vec_nmsub(fs13,dx13,nul);
				fjy3            = vec_nmsub(fs13,dy13,nul);
				fjz3            = vec_nmsub(fs13,dz13,nul);

				fjx1            = vec_nmsub(fs21,dx21,fjx1);
				fjy1            = vec_nmsub(fs21,dy21,fjy1);
				fjz1            = vec_nmsub(fs21,dz21,fjz1);
				fjx2            = vec_nmsub(fs22,dx22,fjx2);
				fjy2            = vec_nmsub(fs22,dy22,fjy2);
				fjz2            = vec_nmsub(fs22,dz22,fjz2);
				fjx3            = vec_nmsub(fs23,dx23,fjx3);
				fjy3            = vec_nmsub(fs23,dy23,fjy3);
				fjz3            = vec_nmsub(fs23,dz23,fjz3);

				fjx1            = vec_nmsub(fs31,dx31,fjx1);
				fjy1            = vec_nmsub(fs31,dy31,fjy1);
				fjz1            = vec_nmsub(fs31,dz31,fjz1);
				fjx2            = vec_nmsub(fs32,dx32,fjx2);
				fjy2            = vec_nmsub(fs32,dy32,fjy2);
				fjz2            = vec_nmsub(fs32,dz32,fjz2);
				fjx3            = vec_nmsub(fs33,dx33,fjx3);
				fjy3            = vec_nmsub(fs33,dy33,fjy3);
				fjz3            = vec_nmsub(fs33,dz33,fjz3);
     
				add_force_to_1_3atoms(faction+j3a,
									  fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,
									  fjx3,fjy3,fjz3);
			}
			/* update outer data */
			update_i_3atoms_forces(faction+ii3,fshift+is3,
								   fix1,fiy1,fiz1,fix2,fiy2,fiz2,
								   fix3,fiy3,fiz3);

			add_vector_to_float(Vc+gid[n],vctot);
			add_vector_to_float(Vvdw+gid[n],Vvdwtot);
			ninner += nj1 - nj0;
		}
#ifdef GMX_THREAD_SHM_FDECOMP
		nouter += nn1 - nn0;
	} while (nn1<nri);
#else
	nouter = nri;
#endif
	*outeriter = nouter;
	*inneriter = ninner;
}







void 
nb_kernel332nf_ppc_altivec(int *             p_nri,
                       int               iinr[],
                       int               jindex[],
                       int               jjnr[],
                       int               shift[],
                       float             shiftvec[],
                       float             fshift[],
                       int               gid[],
                       float             pos[],
                       float             faction[],
                       float             charge[],
                       float *           p_facel,
                       float *           p_krf,
                       float *           p_crf,
                       float             Vc[],
                       int               type[],
                       int *             p_ntype,
                       float             vdwparam[],
                       float             Vvdw[],
                       float *           p_tabscale,
                       float             VFtab[],
                       float             invsqrta[],
                       float             dvda[],
                       float *           p_gbtabscale,
                       float             GBtab[],
                       int *             p_nthreads,
                       int *             count,
                       void *            mtx,
                       int *             outeriter,
                       int *             inneriter,
					   float *           work)
{
	vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3;
	vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3;

	vector float dx11,dy11,dz11,dx12,dy12,dz12,dx13,dy13,dz13;
	vector float dx21,dy21,dz21,dx22,dy22,dz22,dx23,dy23,dz23;
	vector float dx31,dy31,dz31,dx32,dy32,dz32,dx33,dy33,dz33;

	vector float rsq11,rsq12,rsq13,rsq21,rsq22,rsq23,rsq31,rsq32,rsq33;
	vector float r11,r12,r13,r21,r22,r23,r31,r32,r33;
	vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23;
	vector float rinv31,rinv32,rinv33;

	vector float vfacel,nul;
	vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12;
	vector float Vvdwtot,tsc,qqOOt,qqOHt,qqHHt,c6t,c12t;
	vector float VV11c,VV12c,VV13c;
	vector float VV21c,VV22c,VV23c;
	vector float VV31c,VV32c,VV33c;
	vector float VVd,VVr;

	int n,k,ii,is3,ii3,nj0,nj1;
	int jnra,jnrb,jnrc,jnrd,tp,tj;
	int j3a,j3b,j3c,j3d;
	int nri, ntype, nouter, ninner;
#ifdef GMX_THREAD_SHM_FDECOMP
	int nn0, nn1;
#endif

    nouter   = 0;
    ninner   = 0;
    nri      = *p_nri;
    ntype    = *p_ntype;
	nul=vec_zero();
	vfacel=load_float_and_splat(p_facel);
	tsc=load_float_and_splat(p_tabscale);
	ii        = iinr[0];
	qO        = load_float_and_splat(charge+ii);
	qH        = load_float_and_splat(charge+ii+1);
	qqOO      = vec_madd(qO,qO,nul);
	qqOH      = vec_madd(qO,qH,nul);
	qqHH      = vec_madd(qH,qH,nul);
	qqOO      = vec_madd(qqOO,vfacel,nul);
	qqOH      = vec_madd(qqOH,vfacel,nul);
	qqHH      = vec_madd(qqHH,vfacel,nul);
	tp        = 2*type[ii];
	tj        = (ntype+1)*tp;
	load_1_pair(vdwparam+tj,&c6,&c12);
	c6        = vec_splat(c6,0);
	c12       = vec_splat(c12,0);

#ifdef GMX_THREAD_SHM_FDECOMP
    nthreads = *p_nthreads;
	do {
		tMPI_Thread_mutex_lock((tMPI_Thread_mutex_t *)mtx);
		nn0              = *count;
		nn1              = nn0+(nri-nn0)/(2*nthreads)+3;
		*count           = nn1;
		tMPI_Thread_mutex_unlock((tMPI_Thread_mutex_t *)mtx);
		if(nn1>nri) nn1=nri;
		for(n=nn0; (n<nn1); n++) {
#if 0
		} /* maintain correct indentation even with conditional left braces */
#endif
#else /* without tMPI_Threads */
		for(n=0;n<nri;n++) {
#endif  
			is3        = 3*shift[n];
			ii         = iinr[n];
			ii3        = 3*ii;
			load_1_3atoms_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
										  &ix2,&iy2,&iz2,&ix3,&iy3,&iz3);
			vctot      = nul;
			Vvdwtot     = nul;
			nj0        = jindex[n];
			nj1        = jindex[n+1];
    
			for(k=nj0; k<(nj1-3); k+=4) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				jnrd            = jjnr[k+3];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				j3d             = 3*jnrd;
				load_4_3atoms(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
				dx11            = vec_sub(ix1,jx1);
				dx12            = vec_sub(ix1,jx2);
				dx13            = vec_sub(ix1,jx3);
				dy11            = vec_sub(iy1,jy1);
				dy12            = vec_sub(iy1,jy2);
				dy13            = vec_sub(iy1,jy3);
				dz11            = vec_sub(iz1,jz1);
				dz12            = vec_sub(iz1,jz2);
				dz13            = vec_sub(iz1,jz3);
				dx21            = vec_sub(ix2,jx1);
				dx22            = vec_sub(ix2,jx2);
				dx23            = vec_sub(ix2,jx3);
				dy21            = vec_sub(iy2,jy1);
				dy22            = vec_sub(iy2,jy2);
				dy23            = vec_sub(iy2,jy3);
				dz21            = vec_sub(iz2,jz1);
				dz22            = vec_sub(iz2,jz2);
				dz23            = vec_sub(iz2,jz3);
				dx31            = vec_sub(ix3,jx1);
				dx32            = vec_sub(ix3,jx2);
				dx33            = vec_sub(ix3,jx3);
				dy31            = vec_sub(iy3,jy1);
				dy32            = vec_sub(iy3,jy2);
				dy33            = vec_sub(iy3,jy3);
				dz31            = vec_sub(iz3,jz1);
				dz32            = vec_sub(iz3,jz2);
				dz33            = vec_sub(iz3,jz3);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq12           = vec_madd(dx12,dx12,nul);
				rsq13           = vec_madd(dx13,dx13,nul);
				rsq21           = vec_madd(dx21,dx21,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq31           = vec_madd(dx31,dx31,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq12           = vec_madd(dy12,dy12,rsq12);
				rsq13           = vec_madd(dy13,dy13,rsq13);
				rsq21           = vec_madd(dy21,dy21,rsq21);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq31           = vec_madd(dy31,dy31,rsq31);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq12           = vec_madd(dz12,dz12,rsq12);
				rsq13           = vec_madd(dz13,dz13,rsq13);
				rsq21           = vec_madd(dz21,dz21,rsq21);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq31           = vec_madd(dz31,dz31,rsq31);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				r11             = vec_madd(rsq11,rinv11,nul); 
				r12             = vec_madd(rsq12,rinv12,nul); 
				r13             = vec_madd(rsq13,rinv13,nul); 
				r21             = vec_madd(rsq21,rinv21,nul); 
				r22             = vec_madd(rsq22,rinv22,nul); 
				r23             = vec_madd(rsq23,rinv23,nul); 
				r31             = vec_madd(rsq31,rinv31,nul); 
				r32             = vec_madd(rsq32,rinv32,nul); 
				r33             = vec_madd(rsq33,rinv33,nul); 

				do_vonly_4_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
												&VV11c,&VVd,&VVr);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
				vctot           = vec_madd(qqOO,VV11c,vctot);
				vctot           = vec_madd(qqOH,VV12c,vctot);
				vctot           = vec_madd(qqOH,VV13c,vctot);
				vctot           = vec_madd(qqOH,VV21c,vctot);
				vctot           = vec_madd(qqHH,VV22c,vctot);
				vctot           = vec_madd(qqHH,VV23c,vctot);
				vctot           = vec_madd(qqOH,VV31c,vctot);
				vctot           = vec_madd(qqHH,VV32c,vctot);
				vctot           = vec_madd(qqHH,VV33c,vctot);
			} 
			if(k<(nj1-2)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				load_3_3atoms(pos+j3a,pos+j3b,pos+j3c,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
				qqOOt           = vec_sld(qqOO,nul,4);
				qqOHt           = vec_sld(qqOH,nul,4);
				qqHHt           = vec_sld(qqHH,nul,4);
				c6t             = vec_sld(c6,nul,4);
				c12t            = vec_sld(c12,nul,4);

				dx11            = vec_sub(ix1,jx1);
				dx12            = vec_sub(ix1,jx2);
				dx13            = vec_sub(ix1,jx3);
				dy11            = vec_sub(iy1,jy1);
				dy12            = vec_sub(iy1,jy2);
				dy13            = vec_sub(iy1,jy3);
				dz11            = vec_sub(iz1,jz1);
				dz12            = vec_sub(iz1,jz2);
				dz13            = vec_sub(iz1,jz3);
				dx21            = vec_sub(ix2,jx1);
				dx22            = vec_sub(ix2,jx2);
				dx23            = vec_sub(ix2,jx3);
				dy21            = vec_sub(iy2,jy1);
				dy22            = vec_sub(iy2,jy2);
				dy23            = vec_sub(iy2,jy3);
				dz21            = vec_sub(iz2,jz1);
				dz22            = vec_sub(iz2,jz2);
				dz23            = vec_sub(iz2,jz3);
				dx31            = vec_sub(ix3,jx1);
				dx32            = vec_sub(ix3,jx2);
				dx33            = vec_sub(ix3,jx3);
				dy31            = vec_sub(iy3,jy1);
				dy32            = vec_sub(iy3,jy2);
				dy33            = vec_sub(iy3,jy3);
				dz31            = vec_sub(iz3,jz1);
				dz32            = vec_sub(iz3,jz2);
				dz33            = vec_sub(iz3,jz3);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq12           = vec_madd(dx12,dx12,nul);
				rsq13           = vec_madd(dx13,dx13,nul);
				rsq21           = vec_madd(dx21,dx21,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq31           = vec_madd(dx31,dx31,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq12           = vec_madd(dy12,dy12,rsq12);
				rsq13           = vec_madd(dy13,dy13,rsq13);
				rsq21           = vec_madd(dy21,dy21,rsq21);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq31           = vec_madd(dy31,dy31,rsq31);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq12           = vec_madd(dz12,dz12,rsq12);
				rsq13           = vec_madd(dz13,dz13,rsq13);
				rsq21           = vec_madd(dz21,dz21,rsq21);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq31           = vec_madd(dz31,dz31,rsq31);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);

				zero_highest_element_in_9_vectors(&rsq11,&rsq12,&rsq13,
												  &rsq21,&rsq22,&rsq23,
												  &rsq31,&rsq32,&rsq33);

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);

				zero_highest_element_in_9_vectors(&rinv11,&rinv12,&rinv13,
												  &rinv21,&rinv22,&rinv23,
												  &rinv31,&rinv32,&rinv33);

				r11             = vec_madd(rsq11,rinv11,nul); 
				r12             = vec_madd(rsq12,rinv12,nul); 
				r13             = vec_madd(rsq13,rinv13,nul); 
				r21             = vec_madd(rsq21,rinv21,nul); 
				r22             = vec_madd(rsq22,rinv22,nul); 
				r23             = vec_madd(rsq23,rinv23,nul); 
				r31             = vec_madd(rsq31,rinv31,nul); 
				r32             = vec_madd(rsq32,rinv32,nul); 
				r33             = vec_madd(rsq33,rinv33,nul); 

				do_vonly_3_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
												&VV11c,&VVd,&VVr);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

				Vvdwtot          = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot          = vec_madd(c12t,VVr,Vvdwtot);
				vctot           = vec_madd(qqOOt,VV11c,vctot);
				vctot           = vec_madd(qqOHt,VV12c,vctot);
				vctot           = vec_madd(qqOHt,VV13c,vctot);
				vctot           = vec_madd(qqOHt,VV21c,vctot);
				vctot           = vec_madd(qqHHt,VV22c,vctot);
				vctot           = vec_madd(qqHHt,VV23c,vctot);
				vctot           = vec_madd(qqOHt,VV31c,vctot);
				vctot           = vec_madd(qqHHt,VV32c,vctot);
				vctot           = vec_madd(qqHHt,VV33c,vctot);
			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				load_2_3atoms(pos+j3a,pos+j3b,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
				qqOOt           = vec_sld(qqOO,nul,8);
				qqOHt           = vec_sld(qqOH,nul,8);
				qqHHt           = vec_sld(qqHH,nul,8);
				c6t             = vec_sld(c6,nul,8);
				c12t            = vec_sld(c12,nul,8);

				dx11            = vec_sub(ix1,jx1);
				dx12            = vec_sub(ix1,jx2);
				dx13            = vec_sub(ix1,jx3);
				dy11            = vec_sub(iy1,jy1);
				dy12            = vec_sub(iy1,jy2);
				dy13            = vec_sub(iy1,jy3);
				dz11            = vec_sub(iz1,jz1);
				dz12            = vec_sub(iz1,jz2);
				dz13            = vec_sub(iz1,jz3);
				dx21            = vec_sub(ix2,jx1);
				dx22            = vec_sub(ix2,jx2);
				dx23            = vec_sub(ix2,jx3);
				dy21            = vec_sub(iy2,jy1);
				dy22            = vec_sub(iy2,jy2);
				dy23            = vec_sub(iy2,jy3);
				dz21            = vec_sub(iz2,jz1);
				dz22            = vec_sub(iz2,jz2);
				dz23            = vec_sub(iz2,jz3);
				dx31            = vec_sub(ix3,jx1);
				dx32            = vec_sub(ix3,jx2);
				dx33            = vec_sub(ix3,jx3);
				dy31            = vec_sub(iy3,jy1);
				dy32            = vec_sub(iy3,jy2);
				dy33            = vec_sub(iy3,jy3);
				dz31            = vec_sub(iz3,jz1);
				dz32            = vec_sub(iz3,jz2);
				dz33            = vec_sub(iz3,jz3);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq12           = vec_madd(dx12,dx12,nul);
				rsq13           = vec_madd(dx13,dx13,nul);
				rsq21           = vec_madd(dx21,dx21,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq31           = vec_madd(dx31,dx31,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq12           = vec_madd(dy12,dy12,rsq12);
				rsq13           = vec_madd(dy13,dy13,rsq13);
				rsq21           = vec_madd(dy21,dy21,rsq21);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq31           = vec_madd(dy31,dy31,rsq31);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq12           = vec_madd(dz12,dz12,rsq12);
				rsq13           = vec_madd(dz13,dz13,rsq13);
				rsq21           = vec_madd(dz21,dz21,rsq21);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq31           = vec_madd(dz31,dz31,rsq31);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);

				zero_highest_2_elements_in_9_vectors(&rsq11,&rsq12,&rsq13,
													 &rsq21,&rsq22,&rsq23,
													 &rsq31,&rsq32,&rsq33);

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				zero_highest_2_elements_in_9_vectors(&rinv11,&rinv12,&rinv13,
													 &rinv21,&rinv22,&rinv23,
													 &rinv31,&rinv32,&rinv33);

				r11             = vec_madd(rsq11,rinv11,nul); 
				r12             = vec_madd(rsq12,rinv12,nul); 
				r13             = vec_madd(rsq13,rinv13,nul); 
				r21             = vec_madd(rsq21,rinv21,nul); 
				r22             = vec_madd(rsq22,rinv22,nul); 
				r23             = vec_madd(rsq23,rinv23,nul); 
				r31             = vec_madd(rsq31,rinv31,nul); 
				r32             = vec_madd(rsq32,rinv32,nul); 
				r33             = vec_madd(rsq33,rinv33,nul); 

				do_vonly_2_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
												&VV11c,&VVd,&VVr);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

				Vvdwtot          = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot          = vec_madd(c12t,VVr,Vvdwtot);
				vctot           = vec_madd(qqOOt,VV11c,vctot);
				vctot           = vec_madd(qqOHt,VV12c,vctot);
				vctot           = vec_madd(qqOHt,VV13c,vctot);
				vctot           = vec_madd(qqOHt,VV21c,vctot);
				vctot           = vec_madd(qqHHt,VV22c,vctot);
				vctot           = vec_madd(qqHHt,VV23c,vctot);
				vctot           = vec_madd(qqOHt,VV31c,vctot);
				vctot           = vec_madd(qqHHt,VV32c,vctot);
				vctot           = vec_madd(qqHHt,VV33c,vctot);
			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				load_1_3atoms(pos+j3a,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3);
				qqOOt           = vec_sld(qqOO,nul,12);
				qqOHt           = vec_sld(qqOH,nul,12);
				qqHHt           = vec_sld(qqHH,nul,12);
				c6t             = vec_sld(c6,nul,12);
				c12t            = vec_sld(c12,nul,12);

				dx11            = vec_sub(ix1,jx1);
				dx12            = vec_sub(ix1,jx2);
				dx13            = vec_sub(ix1,jx3);
				dy11            = vec_sub(iy1,jy1);
				dy12            = vec_sub(iy1,jy2);
				dy13            = vec_sub(iy1,jy3);
				dz11            = vec_sub(iz1,jz1);
				dz12            = vec_sub(iz1,jz2);
				dz13            = vec_sub(iz1,jz3);
				dx21            = vec_sub(ix2,jx1);
				dx22            = vec_sub(ix2,jx2);
				dx23            = vec_sub(ix2,jx3);
				dy21            = vec_sub(iy2,jy1);
				dy22            = vec_sub(iy2,jy2);
				dy23            = vec_sub(iy2,jy3);
				dz21            = vec_sub(iz2,jz1);
				dz22            = vec_sub(iz2,jz2);
				dz23            = vec_sub(iz2,jz3);
				dx31            = vec_sub(ix3,jx1);
				dx32            = vec_sub(ix3,jx2);
				dx33            = vec_sub(ix3,jx3);
				dy31            = vec_sub(iy3,jy1);
				dy32            = vec_sub(iy3,jy2);
				dy33            = vec_sub(iy3,jy3);
				dz31            = vec_sub(iz3,jz1);
				dz32            = vec_sub(iz3,jz2);
				dz33            = vec_sub(iz3,jz3);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq12           = vec_madd(dx12,dx12,nul);
				rsq13           = vec_madd(dx13,dx13,nul);
				rsq21           = vec_madd(dx21,dx21,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq31           = vec_madd(dx31,dx31,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq12           = vec_madd(dy12,dy12,rsq12);
				rsq13           = vec_madd(dy13,dy13,rsq13);
				rsq21           = vec_madd(dy21,dy21,rsq21);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq31           = vec_madd(dy31,dy31,rsq31);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq12           = vec_madd(dz12,dz12,rsq12);
				rsq13           = vec_madd(dz13,dz13,rsq13);
				rsq21           = vec_madd(dz21,dz21,rsq21);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq31           = vec_madd(dz31,dz31,rsq31);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);

				zero_highest_3_elements_in_9_vectors(&rsq11,&rsq12,&rsq13,
													 &rsq21,&rsq22,&rsq23,
													 &rsq31,&rsq32,&rsq33);

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				zero_highest_3_elements_in_9_vectors(&rinv11,&rinv12,&rinv13,
													 &rinv21,&rinv22,&rinv23,
													 &rinv31,&rinv32,&rinv33);

				r11             = vec_madd(rsq11,rinv11,nul); 
				r12             = vec_madd(rsq12,rinv12,nul); 
				r13             = vec_madd(rsq13,rinv13,nul); 
				r21             = vec_madd(rsq21,rinv21,nul); 
				r22             = vec_madd(rsq22,rinv22,nul); 
				r23             = vec_madd(rsq23,rinv23,nul); 
				r31             = vec_madd(rsq31,rinv31,nul); 
				r32             = vec_madd(rsq32,rinv32,nul); 
				r33             = vec_madd(rsq33,rinv33,nul); 

				do_vonly_1_ljctable_coul_and_lj(VFtab,vec_madd(r11,tsc,nul),
												&VV11c,&VVd,&VVr);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(r12,tsc,nul),&VV12c);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(r13,tsc,nul),&VV13c);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(r21,tsc,nul),&VV21c);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(r22,tsc,nul),&VV22c);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(r23,tsc,nul),&VV23c);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(r31,tsc,nul),&VV31c);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(r32,tsc,nul),&VV32c);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(r33,tsc,nul),&VV33c);

				Vvdwtot          = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot          = vec_madd(c12t,VVr,Vvdwtot);
				vctot           = vec_madd(qqOOt,VV11c,vctot);
				vctot           = vec_madd(qqOHt,VV12c,vctot);
				vctot           = vec_madd(qqOHt,VV13c,vctot);
				vctot           = vec_madd(qqOHt,VV21c,vctot);
				vctot           = vec_madd(qqHHt,VV22c,vctot);
				vctot           = vec_madd(qqHHt,VV23c,vctot);
				vctot           = vec_madd(qqOHt,VV31c,vctot);
				vctot           = vec_madd(qqHHt,VV32c,vctot);
				vctot           = vec_madd(qqHHt,VV33c,vctot);
			}
			/* update outer data */
			add_vector_to_float(Vc+gid[n],vctot);
			add_vector_to_float(Vvdw+gid[n],Vvdwtot);
			ninner += nj1 - nj0;
		}
#ifdef GMX_THREAD_SHM_FDECOMP
		nouter += nn1 - nn0;
	} while (nn1<nri);
#else
	nouter = nri;
#endif
	*outeriter = nouter;
	*inneriter = ninner;
}

