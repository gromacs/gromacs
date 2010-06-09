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
#include "nb_kernel333_ppc_altivec.h"


void 
nb_kernel333_ppc_altivec  (int *             p_nri,
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
	vector float tsc;
	vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z,iMx,iMy,iMz;
	vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z,dMx,dMy,dMz;
	vector float vfacel,nul;
	vector float Vvdwtot,c6,c12;
	vector float fsO,fsH1,fsH2,fsM;
	vector float vctot,qqM,qqH,iqM,iqH,jq;
	vector float fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,fiH2x,fiH2y,fiH2z;
	vector float fiMx,fiMy,fiMz;
	vector float rinvO,rinvH1,rinvH2,rinvM,rsqO,rsqH1,rsqH2,rsqM;
	vector float rO,rH1,rH2,rM,VVcM,FFcM,VVcH1,FFcH1,VVcH2,FFcH2;
	vector float VVd,FFd,VVr,FFr;
	vector float tmp1,tmp2,tmp3,tmp4;

	int n,k,ii,is3,ii3,ntiA,nj0,nj1;
	int jnra,jnrb,jnrc,jnrd;
	int j3a,j3b,j3c,j3d;
	int nri, ntype, nouter, ninner;
	int tja,tjb,tjc,tjd;
#ifdef GMX_THREAD_SHM_FDECOMP
	int nn0, nn1;
#endif

    nouter   = 0;
    ninner   = 0;
    nri      = *p_nri;
    ntype    = *p_ntype;
	nul=vec_zero();
	tsc=load_float_and_splat(p_tabscale);
	vfacel=load_float_and_splat(p_facel);

	ii         = iinr[0];
	iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
	iqM        = vec_madd(load_float_and_splat(charge+ii+3),vfacel,nul);
	ntiA       = 2*ntype*type[ii];
  
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
			load_1_4atoms_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
										  &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z,
										  &iMx,&iMy,&iMz);
			vctot      = nul;
			Vvdwtot     = nul;
			fiOx       = nul;
			fiOy       = nul;
			fiOz       = nul;
			fiH1x      = nul;
			fiH1y      = nul;
			fiH1z      = nul;
			fiH2x      = nul;
			fiH2y      = nul;
			fiH2z      = nul;
			fiMx       = nul;
			fiMy       = nul;
			fiMz       = nul;
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
								 load_xyz(pos+j3d),&dMx,&dMy,&dMz);
				dOx             = vec_sub(iOx,dMx);
				dOy             = vec_sub(iOy,dMy);
				dOz             = vec_sub(iOz,dMz);
				dH1x            = vec_sub(iH1x,dMx);
				dH1y            = vec_sub(iH1y,dMy);
				dH1z            = vec_sub(iH1z,dMz);
				dH2x            = vec_sub(iH2x,dMx);
				dH2y            = vec_sub(iH2y,dMy);
				dH2z            = vec_sub(iH2z,dMz);
				dMx             = vec_sub(iMx,dMx);
				dMy             = vec_sub(iMy,dMy);
				dMz             = vec_sub(iMz,dMz);
      
				rsqO            = vec_madd(dOx,dOx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqO            = vec_madd(dOy,dOy,rsqO);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqO            = vec_madd(dOz,dOz,rsqO);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
				rinvM           = do_invsqrt(rsqM);
				rO              = vec_madd(rsqO,rinvO,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);      
				rM              = vec_madd(rsqM,rinvM,nul);      
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 4 j charges and multiply by iq */
				jq=load_4_float(charge+jnra,charge+jnrb,
								charge+jnrc,charge+jnrd);
				load_4_pair(vdwparam+tja,vdwparam+tjb,vdwparam+tjc,vdwparam+tjd,&c6,&c12);
				do_4_ljctable_lj(VFtab,vec_madd(rO,tsc,nul),&VVd,&FFd,&VVr,&FFr);
				do_4_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
				do_4_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
				do_4_ljctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM,&FFcM);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				fsM             = vec_madd(qqM,FFcM,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				fsO             = vec_madd(c6,FFd,nul);
				fsH1            = vec_madd(qqH,FFcH1,nul);
				fsH2            = vec_madd(qqH,FFcH2,nul);
				fsO             = vec_madd(c12,FFr,fsO);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
				fsO             = vec_nmsub(fsO,tsc,nul);
				vctot           = vec_madd(qqH,VVcH2,vctot);
				fsH1            = vec_nmsub(fsH1,tsc,nul);
				fsH2            = vec_nmsub(fsH2,tsc,nul);
				fsM            = vec_nmsub(fsM,tsc,nul);
				fsO             = vec_madd(fsO,rinvO,nul);
				fsH1            = vec_madd(fsH1,rinvH1,nul);
				fsH2            = vec_madd(fsH2,rinvH2,nul);
				fsM            = vec_madd(fsM,rinvM,nul);

				fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
				dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
				fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
				dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
				fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
				dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
				fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
				dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
				fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
				dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
				fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
				dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
				fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
				dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
				fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
				dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
				fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
				dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */
				fiMx            = vec_madd(fsM,dMx,fiMx); /* +=fx */
				dOx             = vec_nmsub(fsM,dMx,dOx); /* -fx */
				fiMy            = vec_madd(fsM,dMy,fiMy); /* +=fy */
				dOy             = vec_nmsub(fsM,dMy,dOy); /* -fy */
				fiMz            = vec_madd(fsM,dMz,fiMz); /* +=fz */
				dOz             = vec_nmsub(fsM,dMz,dOz); /* -fz */
      
				transpose_3_to_4(dOx,dOy,dOz,&tmp1,&tmp2,&tmp3,&tmp4);
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
								 load_xyz(pos+j3c),
								 nul,&dMx,&dMy,&dMz);
				dOx             = vec_sub(iOx,dMx);
				dOy             = vec_sub(iOy,dMy);
				dOz             = vec_sub(iOz,dMz);
				dH1x            = vec_sub(iH1x,dMx);
				dH1y            = vec_sub(iH1y,dMy);
				dH1z            = vec_sub(iH1z,dMz);
				dH2x            = vec_sub(iH2x,dMx);
				dH2y            = vec_sub(iH2y,dMy);
				dH2z            = vec_sub(iH2z,dMz);
				dMx             = vec_sub(iMx,dMx);
				dMy             = vec_sub(iMy,dMy);
				dMz             = vec_sub(iMz,dMz);
      
				rsqO            = vec_madd(dOx,dOx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqO            = vec_madd(dOy,dOy,rsqO);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqO            = vec_madd(dOz,dOz,rsqO);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
				rinvM           = do_invsqrt(rsqM);

				zero_highest_element_in_3_vectors(&rsqO,&rsqH1,&rsqH2);
				zero_highest_element_in_3_vectors(&rinvO,&rinvH1,&rinvH2);
				zero_highest_element_in_vector(&rsqM);
				zero_highest_element_in_vector(&rinvM);

				rO              = vec_madd(rsqO,rinvO,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);      
				rM              = vec_madd(rsqM,rinvM,nul);      
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 3 j charges and multiply by iq */
				jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
				load_3_pair(vdwparam+tja,vdwparam+tjb,vdwparam+tjc,&c6,&c12);
				do_3_ljctable_lj(VFtab,vec_madd(rO,tsc,nul),
								 &VVd,&FFd,&VVr,&FFr);
				do_3_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
				do_3_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
				do_3_ljctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM,&FFcM);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				fsM             = vec_madd(qqM,FFcM,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				fsO             = vec_madd(c6,FFd,nul);
				fsH1            = vec_madd(qqH,FFcH1,nul);
				fsH2            = vec_madd(qqH,FFcH2,nul);
				fsO             = vec_madd(c12,FFr,fsO);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
				fsO             = vec_nmsub(fsO,tsc,nul);
				vctot           = vec_madd(qqH,VVcH2,vctot);
				fsH1            = vec_nmsub(fsH1,tsc,nul);
				fsH2            = vec_nmsub(fsH2,tsc,nul);
				fsM            = vec_nmsub(fsM,tsc,nul);
				fsO             = vec_madd(fsO,rinvO,nul);
				fsH1            = vec_madd(fsH1,rinvH1,nul);
				fsH2            = vec_madd(fsH2,rinvH2,nul);
				fsM            = vec_madd(fsM,rinvM,nul);

				fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
				dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
				fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
				dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
				fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
				dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
				fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
				dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
				fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
				dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
				fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
				dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
				fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
				dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
				fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
				dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
				fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
				dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */
				fiMx            = vec_madd(fsM,dMx,fiMx); /* +=fx */
				dOx             = vec_nmsub(fsM,dMx,dOx); /* -fx */
				fiMy            = vec_madd(fsM,dMy,fiMy); /* +=fy */
				dOy             = vec_nmsub(fsM,dMy,dOy); /* -fy */
				fiMz            = vec_madd(fsM,dMz,fiMz); /* +=fz */
				dOz             = vec_nmsub(fsM,dMz,dOz); /* -fz */

				transpose_4_to_3(dOx,dOy,dOz,nul,&tmp1,&tmp2,&tmp3);
				add_xyz_to_mem(faction+j3a,tmp1);
				add_xyz_to_mem(faction+j3b,tmp2);
				add_xyz_to_mem(faction+j3c,tmp3);
			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				transpose_2_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),&dMx,&dMy,&dMz);
				dOx             = vec_sub(iOx,dMx);
				dOy             = vec_sub(iOy,dMy);
				dOz             = vec_sub(iOz,dMz);
				dH1x            = vec_sub(iH1x,dMx);
				dH1y            = vec_sub(iH1y,dMy);
				dH1z            = vec_sub(iH1z,dMz);
				dH2x            = vec_sub(iH2x,dMx);
				dH2y            = vec_sub(iH2y,dMy);
				dH2z            = vec_sub(iH2z,dMz);
				dMx             = vec_sub(iMx,dMx);
				dMy             = vec_sub(iMy,dMy);
				dMz             = vec_sub(iMz,dMz);
      
				rsqO            = vec_madd(dOx,dOx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqO            = vec_madd(dOy,dOy,rsqO);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqO            = vec_madd(dOz,dOz,rsqO);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
				rinvM           = do_invsqrt(rsqM);

				zero_highest_2_elements_in_3_vectors(&rsqO,&rsqH1,&rsqH2);
				zero_highest_2_elements_in_3_vectors(&rinvO,&rinvH1,&rinvH2);
				zero_highest_2_elements_in_vector(&rsqM);
				zero_highest_2_elements_in_vector(&rinvM);

				rO              = vec_madd(rsqO,rinvO,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);      
				rM              = vec_madd(rsqM,rinvM,nul);      
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 2 j charges and multiply by iq */
				jq=load_2_float(charge+jnra,charge+jnrb);
				load_2_pair(vdwparam+tja,vdwparam+tjb,&c6,&c12);
				do_2_ljctable_lj(VFtab,vec_madd(rO,tsc,nul),
								 &VVd,&FFd,&VVr,&FFr);
				do_2_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
				do_2_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
				do_2_ljctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM,&FFcM);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				fsM             = vec_madd(qqM,FFcM,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				fsO             = vec_madd(c6,FFd,nul);
				fsH1            = vec_madd(qqH,FFcH1,nul);
				fsH2            = vec_madd(qqH,FFcH2,nul);
				fsO             = vec_madd(c12,FFr,fsO);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
				fsO             = vec_nmsub(fsO,tsc,nul);
				vctot           = vec_madd(qqH,VVcH2,vctot);
				fsH1            = vec_nmsub(fsH1,tsc,nul);
				fsH2            = vec_nmsub(fsH2,tsc,nul);
				fsM            = vec_nmsub(fsM,tsc,nul);
				fsO             = vec_madd(fsO,rinvO,nul);
				fsH1            = vec_madd(fsH1,rinvH1,nul);
				fsH2            = vec_madd(fsH2,rinvH2,nul);
				fsM            = vec_madd(fsM,rinvM,nul);

				fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
				dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
				fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
				dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
				fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
				dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
				fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
				dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
				fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
				dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
				fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
				dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
				fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
				dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
				fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
				dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
				fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
				dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */
				fiMx            = vec_madd(fsM,dMx,fiMx); /* +=fx */
				dOx             = vec_nmsub(fsM,dMx,dOx); /* -fx */
				fiMy            = vec_madd(fsM,dMy,fiMy); /* +=fy */
				dOy             = vec_nmsub(fsM,dMy,dOy); /* -fy */
				fiMz            = vec_madd(fsM,dMz,fiMz); /* +=fz */
				dOz             = vec_nmsub(fsM,dMz,dOz); /* -fz */

				transpose_3_to_2(dOx,dOy,dOz,&tmp1,&tmp2);
				add_xyz_to_mem(faction+j3a,tmp1);
				add_xyz_to_mem(faction+j3b,tmp2);
			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				transpose_1_to_3(load_xyz(pos+j3a),&dMx,&dMy,&dMz);
				dOx             = vec_sub(iOx,dMx);
				dOy             = vec_sub(iOy,dMy);
				dOz             = vec_sub(iOz,dMz);
				dH1x            = vec_sub(iH1x,dMx);
				dH1y            = vec_sub(iH1y,dMy);
				dH1z            = vec_sub(iH1z,dMz);
				dH2x            = vec_sub(iH2x,dMx);
				dH2y            = vec_sub(iH2y,dMy);
				dH2z            = vec_sub(iH2z,dMz);
				dMx             = vec_sub(iMx,dMx);
				dMy             = vec_sub(iMy,dMy);
				dMz             = vec_sub(iMz,dMz);
      
				rsqO            = vec_madd(dOx,dOx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqO            = vec_madd(dOy,dOy,rsqO);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqO            = vec_madd(dOz,dOz,rsqO);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
				rinvM           = do_invsqrt(rsqM);

				zero_highest_3_elements_in_3_vectors(&rsqO,&rsqH1,&rsqH2);
				zero_highest_3_elements_in_3_vectors(&rinvO,&rinvH1,&rinvH2);
				zero_highest_3_elements_in_vector(&rsqM);
				zero_highest_3_elements_in_vector(&rinvM);

				rO              = vec_madd(rsqO,rinvO,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);      
				rM              = vec_madd(rsqM,rinvM,nul);      
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 1 j charge and multiply by iq */
				jq=load_1_float(charge+jnra);
				load_1_pair(vdwparam+tja,&c6,&c12);
				do_1_ljctable_lj(VFtab,vec_madd(rO,tsc,nul),
								 &VVd,&FFd,&VVr,&FFr);
				do_1_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1,&FFcH1);
				do_1_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2,&FFcH2);
				do_1_ljctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM,&FFcM);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				fsM             = vec_madd(qqM,FFcM,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				fsO             = vec_madd(c6,FFd,nul);
				fsH1            = vec_madd(qqH,FFcH1,nul);
				fsH2            = vec_madd(qqH,FFcH2,nul);
				fsO             = vec_madd(c12,FFr,fsO);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
				fsO             = vec_nmsub(fsO,tsc,nul);
				vctot           = vec_madd(qqH,VVcH2,vctot);
				fsH1            = vec_nmsub(fsH1,tsc,nul);
				fsH2            = vec_nmsub(fsH2,tsc,nul);
				fsM            = vec_nmsub(fsM,tsc,nul);
				fsO             = vec_madd(fsO,rinvO,nul);
				fsH1            = vec_madd(fsH1,rinvH1,nul);
				fsH2            = vec_madd(fsH2,rinvH2,nul);
				fsM            = vec_madd(fsM,rinvM,nul);

				fiOx            = vec_madd(fsO,dOx,fiOx); /* +=fx */
				dOx             = vec_nmsub(fsO,dOx,nul); /* -fx */
				fiOy            = vec_madd(fsO,dOy,fiOy); /* +=fy */
				dOy             = vec_nmsub(fsO,dOy,nul); /* -fy */
				fiOz            = vec_madd(fsO,dOz,fiOz); /* +=fz */
				dOz             = vec_nmsub(fsO,dOz,nul); /* -fz */
				fiH1x           = vec_madd(fsH1,dH1x,fiH1x); /* +=fx */
				dOx             = vec_nmsub(fsH1,dH1x,dOx); /* -fx */
				fiH1y           = vec_madd(fsH1,dH1y,fiH1y); /* +=fy */
				dOy             = vec_nmsub(fsH1,dH1y,dOy); /* -fy */
				fiH1z           = vec_madd(fsH1,dH1z,fiH1z); /* +=fz */
				dOz             = vec_nmsub(fsH1,dH1z,dOz); /* -fz */
				fiH2x           = vec_madd(fsH2,dH2x,fiH2x); /* +=fx */
				dOx             = vec_nmsub(fsH2,dH2x,dOx); /* -fx */
				fiH2y           = vec_madd(fsH2,dH2y,fiH2y); /* +=fy */
				dOy             = vec_nmsub(fsH2,dH2y,dOy); /* -fy */
				fiH2z           = vec_madd(fsH2,dH2z,fiH2z); /* +=fz */
				dOz             = vec_nmsub(fsH2,dH2z,dOz); /* -fz */
				fiMx            = vec_madd(fsM,dMx,fiMx); /* +=fx */
				dOx             = vec_nmsub(fsM,dMx,dOx); /* -fx */
				fiMy            = vec_madd(fsM,dMy,fiMy); /* +=fy */
				dOy             = vec_nmsub(fsM,dMy,dOy); /* -fy */
				fiMz            = vec_madd(fsM,dMz,fiMz); /* +=fz */
				dOz             = vec_nmsub(fsM,dMz,dOz); /* -fz */

				transpose_3_to_1(dOx,dOy,dOz,&tmp1);
				add_xyz_to_mem(faction+j3a,tmp1);
			}
			/* update outer data */
			update_i_4atoms_forces(faction+ii3,fshift+is3,
								   fiOx,fiOy,fiOz,fiH1x,fiH1y,fiH1z,
								   fiH2x,fiH2y,fiH2z,fiMx,fiMy,fiMz);

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
nb_kernel333nf_ppc_altivec(int *             p_nri,
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
	vector float tsc;
	vector float iOx,iOy,iOz,iH1x,iH1y,iH1z,iH2x,iH2y,iH2z,iMx,iMy,iMz;
	vector float dOx,dOy,dOz,dH1x,dH1y,dH1z,dH2x,dH2y,dH2z,dMx,dMy,dMz;
	vector float vfacel,nul;
	vector float Vvdwtot,c6,c12;
	vector float vctot,qqM,qqH,iqM,iqH,jq;
	vector float rinvO,rinvH1,rinvH2,rinvM,rsqO,rsqH1,rsqH2,rsqM;
	vector float rO,rH1,rH2,rM,VVcM,VVcH1,VVcH2,VVd,VVr;

	int n,k,ii,is3,ii3,ntiA,nj0,nj1;
	int jnra,jnrb,jnrc,jnrd;
	int j3a,j3b,j3c,j3d;
	int nri, ntype, nouter, ninner;
	int tja,tjb,tjc,tjd;
#ifdef GMX_THREAD_SHM_FDECOMP
	int nn0, nn1;
#endif

    nouter   = 0;
    ninner   = 0;
    nri      = *p_nri;
    ntype    = *p_ntype;
	nul=vec_zero();
	tsc=load_float_and_splat(p_tabscale);
	vfacel=load_float_and_splat(p_facel);

	ii         = iinr[0];
	iqH        = vec_madd(load_float_and_splat(charge+ii+1),vfacel,nul);
	iqM        = vec_madd(load_float_and_splat(charge+ii+3),vfacel,nul);
	ntiA       = 2*ntype*type[ii];
  
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
			load_1_4atoms_shift_and_splat(pos+ii3,shiftvec+is3,&iOx,&iOy,&iOz,
										  &iH1x,&iH1y,&iH1z,&iH2x,&iH2y,&iH2z,
										  &iMx,&iMy,&iMz);
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
				transpose_4_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),
								 load_xyz(pos+j3c),
								 load_xyz(pos+j3d),&dMx,&dMy,&dMz);
				dOx             = vec_sub(iOx,dMx);
				dOy             = vec_sub(iOy,dMy);
				dOz             = vec_sub(iOz,dMz);
				dH1x            = vec_sub(iH1x,dMx);
				dH1y            = vec_sub(iH1y,dMy);
				dH1z            = vec_sub(iH1z,dMz);
				dH2x            = vec_sub(iH2x,dMx);
				dH2y            = vec_sub(iH2y,dMy);
				dH2z            = vec_sub(iH2z,dMz);
				dMx             = vec_sub(iMx,dMx);
				dMy             = vec_sub(iMy,dMy);
				dMz             = vec_sub(iMz,dMz);
      
				rsqO            = vec_madd(dOx,dOx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqO            = vec_madd(dOy,dOy,rsqO);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqO            = vec_madd(dOz,dOz,rsqO);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
				rinvM           = do_invsqrt(rsqM);
				rO              = vec_madd(rsqO,rinvO,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);      
				rM              = vec_madd(rsqM,rinvM,nul);      
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 4 j charges and multiply by iq */
				jq=load_4_float(charge+jnra,charge+jnrb,
								charge+jnrc,charge+jnrd);
				load_4_pair(vdwparam+tja,vdwparam+tjb,vdwparam+tjc,vdwparam+tjd,&c6,&c12);
				do_vonly_4_ljctable_lj(VFtab,vec_madd(rO,tsc,nul),&VVd,&VVr);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
				do_vonly_4_ljctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
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
								 load_xyz(pos+j3c),
								 nul,&dMx,&dMy,&dMz);
				dOx             = vec_sub(iOx,dMx);
				dOy             = vec_sub(iOy,dMy);
				dOz             = vec_sub(iOz,dMz);
				dH1x            = vec_sub(iH1x,dMx);
				dH1y            = vec_sub(iH1y,dMy);
				dH1z            = vec_sub(iH1z,dMz);
				dH2x            = vec_sub(iH2x,dMx);
				dH2y            = vec_sub(iH2y,dMy);
				dH2z            = vec_sub(iH2z,dMz);
				dMx             = vec_sub(iMx,dMx);
				dMy             = vec_sub(iMy,dMy);
				dMz             = vec_sub(iMz,dMz);
      
				rsqO            = vec_madd(dOx,dOx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqO            = vec_madd(dOy,dOy,rsqO);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqO            = vec_madd(dOz,dOz,rsqO);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
				rinvM           = do_invsqrt(rsqM);

				zero_highest_element_in_3_vectors(&rsqO,&rsqH1,&rsqH2);
				zero_highest_element_in_3_vectors(&rinvO,&rinvH1,&rinvH2);
				zero_highest_element_in_vector(&rsqM);
				zero_highest_element_in_vector(&rinvM);

				rO              = vec_madd(rsqO,rinvO,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);      
				rM              = vec_madd(rsqM,rinvM,nul);      
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 3 j charges and multiply by iq */
				jq=load_3_float(charge+jnra,charge+jnrb,charge+jnrc);
				load_3_pair(vdwparam+tja,vdwparam+tjb,vdwparam+tjc,&c6,&c12);
				do_vonly_3_ljctable_lj(VFtab,vec_madd(rO,tsc,nul),&VVd,&VVr);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
				do_vonly_3_ljctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
				vctot           = vec_madd(qqH,VVcH2,vctot);
			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				transpose_2_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),&dMx,&dMy,&dMz);
				dOx             = vec_sub(iOx,dMx);
				dOy             = vec_sub(iOy,dMy);
				dOz             = vec_sub(iOz,dMz);
				dH1x            = vec_sub(iH1x,dMx);
				dH1y            = vec_sub(iH1y,dMy);
				dH1z            = vec_sub(iH1z,dMz);
				dH2x            = vec_sub(iH2x,dMx);
				dH2y            = vec_sub(iH2y,dMy);
				dH2z            = vec_sub(iH2z,dMz);
				dMx             = vec_sub(iMx,dMx);
				dMy             = vec_sub(iMy,dMy);
				dMz             = vec_sub(iMz,dMz);
      
				rsqO            = vec_madd(dOx,dOx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqO            = vec_madd(dOy,dOy,rsqO);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqO            = vec_madd(dOz,dOz,rsqO);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
				rinvM           = do_invsqrt(rsqM);

				zero_highest_2_elements_in_3_vectors(&rsqO,&rsqH1,&rsqH2);
				zero_highest_2_elements_in_3_vectors(&rinvO,&rinvH1,&rinvH2);
				zero_highest_2_elements_in_vector(&rsqM);
				zero_highest_2_elements_in_vector(&rinvM);

				rO              = vec_madd(rsqO,rinvO,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);      
				rM              = vec_madd(rsqM,rinvM,nul);      
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 2 j charges and multiply by iq */
				jq=load_2_float(charge+jnra,charge+jnrb);
				load_2_pair(vdwparam+tja,vdwparam+tjb,&c6,&c12);
				do_vonly_2_ljctable_lj(VFtab,vec_madd(rO,tsc,nul),&VVd,&VVr);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
				do_vonly_2_ljctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
				vctot           = vec_madd(qqH,VVcH2,vctot);
			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				transpose_1_to_3(load_xyz(pos+j3a),&dMx,&dMy,&dMz);
				dOx             = vec_sub(iOx,dMx);
				dOy             = vec_sub(iOy,dMy);
				dOz             = vec_sub(iOz,dMz);
				dH1x            = vec_sub(iH1x,dMx);
				dH1y            = vec_sub(iH1y,dMy);
				dH1z            = vec_sub(iH1z,dMz);
				dH2x            = vec_sub(iH2x,dMx);
				dH2y            = vec_sub(iH2y,dMy);
				dH2z            = vec_sub(iH2z,dMz);
				dMx             = vec_sub(iMx,dMx);
				dMy             = vec_sub(iMy,dMy);
				dMz             = vec_sub(iMz,dMz);
      
				rsqO            = vec_madd(dOx,dOx,nul);
				rsqH1           = vec_madd(dH1x,dH1x,nul);
				rsqH2           = vec_madd(dH2x,dH2x,nul);
				rsqM            = vec_madd(dMx,dMx,nul);
				rsqO            = vec_madd(dOy,dOy,rsqO);
				rsqH1           = vec_madd(dH1y,dH1y,rsqH1);
				rsqH2           = vec_madd(dH2y,dH2y,rsqH2);
				rsqM            = vec_madd(dMy,dMy,rsqM);
				rsqO            = vec_madd(dOz,dOz,rsqO);
				rsqH1           = vec_madd(dH1z,dH1z,rsqH1);
				rsqH2           = vec_madd(dH2z,dH2z,rsqH2);
				rsqM            = vec_madd(dMz,dMz,rsqM);
				do_3_invsqrt(rsqO,rsqH1,rsqH2,&rinvO,&rinvH1,&rinvH2);
				rinvM           = do_invsqrt(rsqM);

				zero_highest_3_elements_in_3_vectors(&rsqO,&rsqH1,&rsqH2);
				zero_highest_3_elements_in_3_vectors(&rinvO,&rinvH1,&rinvH2);
				zero_highest_3_elements_in_vector(&rsqM);
				zero_highest_3_elements_in_vector(&rinvM);

				rO              = vec_madd(rsqO,rinvO,nul);
				rH1             = vec_madd(rsqH1,rinvH1,nul);
				rH2             = vec_madd(rsqH2,rinvH2,nul);      
				rM              = vec_madd(rsqM,rinvM,nul);      
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 1 j charge and multiply by iq */
				jq=load_1_float(charge+jnra);
				load_1_pair(vdwparam+tja,&c6,&c12);
				do_vonly_1_ljctable_lj(VFtab,vec_madd(rO,tsc,nul),&VVd,&VVr);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(rH1,tsc,nul),&VVcH1);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(rH2,tsc,nul),&VVcH2);
				do_vonly_1_ljctable_coul(VFtab,vec_madd(rM,tsc,nul),&VVcM);
				qqM             = vec_madd(iqM,jq,nul);
				qqH             = vec_madd(iqH,jq,nul);
				vctot           = vec_madd(qqM,VVcM,vctot);
				Vvdwtot          = vec_madd(c6,VVd,Vvdwtot);
				vctot           = vec_madd(qqH,VVcH1,vctot);
				Vvdwtot          = vec_madd(c12,VVr,Vvdwtot);
				vctot           = vec_madd(qqH,VVcH2,vctot);
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

