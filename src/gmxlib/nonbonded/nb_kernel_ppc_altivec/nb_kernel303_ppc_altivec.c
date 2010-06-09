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
#include "nb_kernel303_ppc_altivec.h"


void 
nb_kernel303_ppc_altivec  (int *             p_nri,
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
	vector float iMx,iMy,iMz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
	vector float dMx,dMy,dMz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
	vector float vfacel,nul;
	vector float fsM,fsH1,fsH2,tsc,VVcM,FFcM,VVcH1,FFcH1,VVcH2,FFcH2;
	vector float vctot,qqM,qqH,iqM,iqH,jq;
	vector float fiMx,fiMy,fiMz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
	vector float tmp1,tmp2,tmp3,tmp4;
	vector float rinvM,rinvH1,rinvH2,rM,rH1,rH2,rsqM,rsqH1,rsqH2;

	int n,k,ii,is3,ii3,nj0,nj1;
	int jnra,jnrb,jnrc,jnrd;
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
	iqH        = vec_madd(load_float_and_splat(charge+iinr[0]+1),vfacel,nul);
	iqM        = vec_madd(load_float_and_splat(charge+iinr[0]+3),vfacel,nul);
  
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
			load_1_3atoms_shift_and_splat(pos+ii3+3,shiftvec+is3,
										  &iH1x,&iH1y,&iH1z,
										  &iH2x,&iH2y,&iH2z,&iMx,&iMy,&iMz);
			vctot      = nul;
			fiMx       = nul;
			fiMy       = nul;
			fiMz       = nul;
			fiH1x      = nul;
			fiH1y      = nul;
			fiH1z      = nul;
			fiH2x      = nul;
			fiH2y      = nul;
			fiH2z      = nul;
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
				transpose_4_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),
								 load_xyz(pos+j3c),
								 load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
				dMx             = vec_sub(iMx,dH2x);
				dMy             = vec_sub(iMy,dH2y);
				dMz             = vec_sub(iMz,dH2z);
				dH1x            = vec_sub(iH1x,dH2x);
				dH1y            = vec_sub(iH1y,dH2y);
				dH1z            = vec_sub(iH1z,dH2z);
				dH2x            = vec_sub(iH2x,dH2x);
				dH2y            = vec_sub(iH2y,dH2y);
				dH2z            = vec_sub(iH2z,dH2z);
      
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				do_3_invsqrt(rsqM,rsqH1,rsqH2,&rinvM,&rinvH1,&rinvH2);
				rM              = vec_madd(rsqM,rinvM,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);
    
				/* load 4 j charges and multiply by iq */
				jq=load_4_float(charge+jnra,charge+jnrb,
								charge+jnrc,charge+jnrd);
				do_4_ctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM,&FFcM);
				do_4_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
				do_4_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				fsM             = vec_nmsub(qqM,FFcM,nul);
				fsH1            = vec_nmsub(qqH,FFcH1,nul);
				fsH2            = vec_nmsub(qqH,FFcH2,nul);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				fsM             = vec_madd(fsM,tsc,nul);
				fsH1            = vec_madd(fsH1,tsc,nul);
				fsH2            = vec_madd(fsH2,tsc,nul);
				vctot           = vec_madd(qqH,VVcH2,vctot);
				fsM             = vec_madd(fsM,rinvM,nul);
				fsH1            = vec_madd(fsH1,rinvH1,nul);
				fsH2            = vec_madd(fsH2,rinvH2,nul);
      
				fiMx            = vec_madd(fsM,dMx,fiMx); /* +=fx */
				dMx             = vec_nmsub(fsM,dMx,nul); /* -fx */
				fiMy            = vec_madd(fsM,dMy,fiMy); /* +=fy */
				dMy             = vec_nmsub(fsM,dMy,nul); /* -fy */
				fiMz            = vec_madd(fsM,dMz,fiMz); /* +=fz */
				dMz             = vec_nmsub(fsM,dMz,nul); /* -fz */
				fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
				dMx             = vec_nmsub(fsH1,dH1x,dMx); /* -fx */
				fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
				dMy             = vec_nmsub(fsH1,dH1y,dMy); /* -fy */
				fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
				dMz             = vec_nmsub(fsH1,dH1z,dMz); /* -fz */
				fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
				dMx             = vec_nmsub(fsH2,dH2x,dMx); /* -fx */
				fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
				dMy             = vec_nmsub(fsH2,dH2y,dMy); /* -fy */
				fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
				dMz             = vec_nmsub(fsH2,dH2z,dMz); /* -fz */

				transpose_3_to_4(dMx,dMy,dMz,&tmp1,&tmp2,&tmp3,&tmp4);
				add_xyz_to_mem(faction+j3a,tmp1);
				add_xyz_to_mem(faction+j3b,tmp2);
				add_xyz_to_mem(faction+j3c,tmp3);
				add_xyz_to_mem(faction+j3d,tmp4);
			} 
			if(k<(nj1-2)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				transpose_4_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),
								 load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
				dMx             = vec_sub(iMx,dH2x);
				dMy             = vec_sub(iMy,dH2y);
				dMz             = vec_sub(iMz,dH2z);
				dH1x            = vec_sub(iH1x,dH2x);
				dH1y            = vec_sub(iH1y,dH2y);
				dH1z            = vec_sub(iH1z,dH2z);
				dH2x            = vec_sub(iH2x,dH2x);
				dH2y            = vec_sub(iH2y,dH2y);
				dH2z            = vec_sub(iH2z,dH2z);
      
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);

				zero_highest_element_in_3_vectors(&rsqM,&rsqH1,&rsqH2);
				do_3_invsqrt(rsqM,rsqH1,rsqH2,&rinvM,&rinvH1,&rinvH2);
				zero_highest_element_in_3_vectors(&rinvM,&rinvH1,&rinvH2);

				rM              = vec_madd(rsqM,rinvM,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);

				/* load 3 j charges and multiply by iq */
				jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
				do_3_ctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM,&FFcM);
				do_3_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
				do_3_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				fsM             = vec_nmsub(qqM,FFcM,nul);
				fsH1            = vec_nmsub(qqH,FFcH1,nul);
				fsH2            = vec_nmsub(qqH,FFcH2,nul);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				fsM             = vec_madd(fsM,tsc,nul);
				fsH1            = vec_madd(fsH1,tsc,nul);
				fsH2            = vec_madd(fsH2,tsc,nul);
				vctot           = vec_madd(qqH,VVcH2,vctot);
				fsM             = vec_madd(fsM,rinvM,nul);
				fsH1            = vec_madd(fsH1,rinvH1,nul);
				fsH2            = vec_madd(fsH2,rinvH2,nul);
      
				fiMx            = vec_madd(fsM,dMx,fiMx); /* +=fx */
				dMx             = vec_nmsub(fsM,dMx,nul); /* -fx */
				fiMy            = vec_madd(fsM,dMy,fiMy); /* +=fy */
				dMy             = vec_nmsub(fsM,dMy,nul); /* -fy */
				fiMz            = vec_madd(fsM,dMz,fiMz); /* +=fz */
				dMz             = vec_nmsub(fsM,dMz,nul); /* -fz */
				fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
				dMx             = vec_nmsub(fsH1,dH1x,dMx); /* -fx */
				fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
				dMy             = vec_nmsub(fsH1,dH1y,dMy); /* -fy */
				fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
				dMz             = vec_nmsub(fsH1,dH1z,dMz); /* -fz */
				fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
				dMx             = vec_nmsub(fsH2,dH2x,dMx); /* -fx */
				fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
				dMy             = vec_nmsub(fsH2,dH2y,dMy); /* -fy */
				fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
				dMz             = vec_nmsub(fsH2,dH2z,dMz); /* -fz */

				transpose_4_to_3(dMx,dMy,dMz,nul,&tmp1,&tmp2,&tmp3);
				add_xyz_to_mem(faction+j3a,tmp1);
				add_xyz_to_mem(faction+j3b,tmp2);
				add_xyz_to_mem(faction+j3c,tmp3);
			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				transpose_2_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
				dMx             = vec_sub(iMx,dH2x);
				dMy             = vec_sub(iMy,dH2y);
				dMz             = vec_sub(iMz,dH2z);
				dH1x            = vec_sub(iH1x,dH2x);
				dH1y            = vec_sub(iH1y,dH2y);
				dH1z            = vec_sub(iH1z,dH2z);
				dH2x            = vec_sub(iH2x,dH2x);
				dH2y            = vec_sub(iH2y,dH2y);
				dH2z            = vec_sub(iH2z,dH2z);
      
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);

				zero_highest_2_elements_in_3_vectors(&rsqM,&rsqH1,&rsqH2);
				do_3_invsqrt(rsqM,rsqH1,rsqH2,&rinvM,&rinvH1,&rinvH2);
				zero_highest_2_elements_in_3_vectors(&rinvM,&rinvH1,&rinvH2);

				rM              = vec_madd(rsqM,rinvM,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);
    
				/* load 2 j charges and multiply by iq */
				jq=load_2_float(charge+jnra,charge+jnrb);
				do_2_ctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM,&FFcM);
				do_2_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
				do_2_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				fsM             = vec_nmsub(qqM,FFcM,nul);
				fsH1            = vec_nmsub(qqH,FFcH1,nul);
				fsH2            = vec_nmsub(qqH,FFcH2,nul);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				fsM             = vec_madd(fsM,tsc,nul);
				fsH1            = vec_madd(fsH1,tsc,nul);
				fsH2            = vec_madd(fsH2,tsc,nul);
				vctot           = vec_madd(qqH,VVcH2,vctot);
				fsM             = vec_madd(fsM,rinvM,nul);
				fsH1            = vec_madd(fsH1,rinvH1,nul);
				fsH2            = vec_madd(fsH2,rinvH2,nul);
 
				fiMx            = vec_madd(fsM,dMx,fiMx); /* +=fx */
				dMx             = vec_nmsub(fsM,dMx,nul); /* -fx */
				fiMy            = vec_madd(fsM,dMy,fiMy); /* +=fy */
				dMy             = vec_nmsub(fsM,dMy,nul); /* -fy */
				fiMz            = vec_madd(fsM,dMz,fiMz); /* +=fz */
				dMz             = vec_nmsub(fsM,dMz,nul); /* -fz */
				fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
				dMx             = vec_nmsub(fsH1,dH1x,dMx); /* -fx */
				fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
				dMy             = vec_nmsub(fsH1,dH1y,dMy); /* -fy */
				fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
				dMz             = vec_nmsub(fsH1,dH1z,dMz); /* -fz */
				fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
				dMx             = vec_nmsub(fsH2,dH2x,dMx); /* -fx */
				fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
				dMy             = vec_nmsub(fsH2,dH2y,dMy); /* -fy */
				fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
				dMz             = vec_nmsub(fsH2,dH2z,dMz); /* -fz */

				transpose_3_to_2(dMx,dMy,dMz,&tmp1,&tmp2);
				add_xyz_to_mem(faction+j3a,tmp1);
				add_xyz_to_mem(faction+j3b,tmp2);
			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
				dMx             = vec_sub(iMx,dH2x);
				dMy             = vec_sub(iMy,dH2y);
				dMz             = vec_sub(iMz,dH2z);
				dH1x            = vec_sub(iH1x,dH2x);
				dH1y            = vec_sub(iH1y,dH2y);
				dH1z            = vec_sub(iH1z,dH2z);
				dH2x            = vec_sub(iH2x,dH2x);
				dH2y            = vec_sub(iH2y,dH2y);
				dH2z            = vec_sub(iH2z,dH2z);
      
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);

				zero_highest_3_elements_in_3_vectors(&rsqM,&rsqH1,&rsqH2);
				do_3_invsqrt(rsqM,rsqH1,rsqH2,&rinvM,&rinvH1,&rinvH2);
				zero_highest_3_elements_in_3_vectors(&rinvM,&rinvH1,&rinvH2);

				rM              = vec_madd(rsqM,rinvM,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);

				/* load 1 j charges and multiply by iq */
				jq=load_1_float(charge+jnra);
				do_1_ctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM,&FFcM);
				do_1_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
				do_1_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				fsM             = vec_nmsub(qqM,FFcM,nul);
				fsH1            = vec_nmsub(qqH,FFcH1,nul);
				fsH2            = vec_nmsub(qqH,FFcH2,nul);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				fsM             = vec_madd(fsM,tsc,nul);
				fsH1            = vec_madd(fsH1,tsc,nul);
				fsH2            = vec_madd(fsH2,tsc,nul);
				vctot           = vec_madd(qqH,VVcH2,vctot);
				fsM             = vec_madd(fsM,rinvM,nul);
				fsH1            = vec_madd(fsH1,rinvH1,nul);
				fsH2            = vec_madd(fsH2,rinvH2,nul);
      
				fiMx            = vec_madd(fsM,dMx,fiMx); /* +=fx */
				dMx             = vec_nmsub(fsM,dMx,nul); /* -fx */
				fiMy            = vec_madd(fsM,dMy,fiMy); /* +=fy */
				dMy             = vec_nmsub(fsM,dMy,nul); /* -fy */
				fiMz            = vec_madd(fsM,dMz,fiMz); /* +=fz */
				dMz             = vec_nmsub(fsM,dMz,nul); /* -fz */
				fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
				dMx             = vec_nmsub(fsH1,dH1x,dMx); /* -fx */
				fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
				dMy             = vec_nmsub(fsH1,dH1y,dMy); /* -fy */
				fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
				dMz             = vec_nmsub(fsH1,dH1z,dMz); /* -fz */
				fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
				dMx             = vec_nmsub(fsH2,dH2x,dMx); /* -fx */
				fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
				dMy             = vec_nmsub(fsH2,dH2y,dMy); /* -fy */
				fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
				dMz             = vec_nmsub(fsH2,dH2z,dMz); /* -fz */

				transpose_3_to_1(dMx,dMy,dMz,&tmp1);
				add_xyz_to_mem(faction+j3a,tmp1);
			}
			/* update outer data */
			update_i_3atoms_forces(faction+ii3+3,fshift+is3,
								   fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z,
								   fiMx,fiMy,fiMz);

			add_vector_to_float(Vc+gid[n],vctot);
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
nb_kernel303nf_ppc_altivec(int *             p_nri,
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
	vector float iMx,iMy,iMz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z;
	vector float dMx,dMy,dMz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z;
	vector float vfacel,nul;
	vector float tsc,VVcM,VVcH1,VVcH2;
	vector float vctot,qqM,qqH,iqM,iqH,jq;
	vector float rinvM,rinvH1,rinvH2,rM,rH1,rH2,rsqM,rsqH1,rsqH2;  

	int n,k,ii,is3,ii3,nj0,nj1;
	int jnra,jnrb,jnrc,jnrd;
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
	iqH        = vec_madd(load_float_and_splat(charge+iinr[0]+1),vfacel,nul);
	iqM        = vec_madd(load_float_and_splat(charge+iinr[0]+3),vfacel,nul);
  
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
			load_1_3atoms_shift_and_splat(pos+ii3+3,shiftvec+is3,
										  &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z,
										  &iMx,&iMy,&iMz);
			vctot      = nul;
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
				transpose_4_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),
								 load_xyz(pos+j3c),
								 load_xyz(pos+j3d),&dH2x,&dH2y,&dH2z);
				dMx             = vec_sub(iMx,dH2x);
				dMy             = vec_sub(iMy,dH2y);
				dMz             = vec_sub(iMz,dH2z);
				dH1x            = vec_sub(iH1x,dH2x);
				dH1y            = vec_sub(iH1y,dH2y);
				dH1z            = vec_sub(iH1z,dH2z);
				dH2x            = vec_sub(iH2x,dH2x);
				dH2y            = vec_sub(iH2y,dH2y);
				dH2z            = vec_sub(iH2z,dH2z);
      
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				do_3_invsqrt(rsqM,rsqH1,rsqH2,&rinvM,&rinvH1,&rinvH2);
				rM              = vec_madd(rsqM,rinvM,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);
    
				/* load 4 j charges and multiply by iq */
				jq=load_4_float(charge+jnra,charge+jnrb,
								charge+jnrc,charge+jnrd);
				do_vonly_4_ctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM);
				do_vonly_4_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
				do_vonly_4_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				vctot           = vec_madd(qqH,VVcH2,vctot);
			} 
			if(k<(nj1-2)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				transpose_4_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),
								 load_xyz(pos+j3c),nul,&dH2x,&dH2y,&dH2z);
				dMx             = vec_sub(iMx,dH2x);
				dMy             = vec_sub(iMy,dH2y);
				dMz             = vec_sub(iMz,dH2z);
				dH1x            = vec_sub(iH1x,dH2x);
				dH1y            = vec_sub(iH1y,dH2y);
				dH1z            = vec_sub(iH1z,dH2z);
				dH2x            = vec_sub(iH2x,dH2x);
				dH2y            = vec_sub(iH2y,dH2y);
				dH2z            = vec_sub(iH2z,dH2z);
      
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				zero_highest_element_in_3_vectors(&rsqM,&rsqH1,&rsqH2);
				do_3_invsqrt(rsqM,rsqH1,rsqH2,&rinvM,&rinvH1,&rinvH2);
				zero_highest_element_in_3_vectors(&rinvM,&rinvH1,&rinvH2);
				rM              = vec_madd(rsqM,rinvM,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);

				/* load 3 j charges and multiply by iq */
				jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
				do_vonly_3_ctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM);
				do_vonly_3_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
				do_vonly_3_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				vctot           = vec_madd(qqH,VVcH2,vctot);
			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				transpose_2_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),&dH2x,&dH2y,&dH2z);
				dMx             = vec_sub(iMx,dH2x);
				dMy             = vec_sub(iMy,dH2y);
				dMz             = vec_sub(iMz,dH2z);
				dH1x            = vec_sub(iH1x,dH2x);
				dH1y            = vec_sub(iH1y,dH2y);
				dH1z            = vec_sub(iH1z,dH2z);
				dH2x            = vec_sub(iH2x,dH2x);
				dH2y            = vec_sub(iH2y,dH2y);
				dH2z            = vec_sub(iH2z,dH2z);
      
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				zero_highest_2_elements_in_3_vectors(&rsqM,&rsqH1,&rsqH2);
				do_3_invsqrt(rsqM,rsqH1,rsqH2,&rinvM,&rinvH1,&rinvH2);
				zero_highest_2_elements_in_3_vectors(&rinvM,&rinvH1,&rinvH2);
				rM              = vec_madd(rsqM,rinvM,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);
    
				/* load 2 j charges and multiply by iq */
				jq=load_2_float(charge+jnra,charge+jnrb);
				do_vonly_2_ctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM);
				do_vonly_2_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
				do_vonly_2_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				vctot           = vec_madd(qqH,VVcH2,vctot);
			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				transpose_1_to_3(load_xyz(pos+j3a),&dH2x,&dH2y,&dH2z);
				dMx             = vec_sub(iMx,dH2x);
				dMy             = vec_sub(iMy,dH2y);
				dMz             = vec_sub(iMz,dH2z);
				dH1x            = vec_sub(iH1x,dH2x);
				dH1y            = vec_sub(iH1y,dH2y);
				dH1z            = vec_sub(iH1z,dH2z);
				dH2x            = vec_sub(iH2x,dH2x);
				dH2y            = vec_sub(iH2y,dH2y);
				dH2z            = vec_sub(iH2z,dH2z);
      
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				zero_highest_3_elements_in_3_vectors(&rsqM,&rsqH1,&rsqH2);
				do_3_invsqrt(rsqM,rsqH1,rsqH2,&rinvM,&rinvH1,&rinvH2);
				zero_highest_3_elements_in_3_vectors(&rinvM,&rinvH1,&rinvH2);
				rM              = vec_madd(rsqM,rinvM,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);

				/* load 1 j charges and multiply by iq */
				jq=load_1_float(charge+jnra);
				do_vonly_1_ctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM);
				do_vonly_1_ctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
				do_vonly_1_ctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				vctot           = vec_madd(qqH,VVcH2,vctot);
			}
			/* update outer data */
			add_vector_to_float(Vc+gid[n],vctot);
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
