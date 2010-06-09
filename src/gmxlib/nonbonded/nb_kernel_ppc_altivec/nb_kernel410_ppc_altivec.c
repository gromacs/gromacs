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
#include "nb_kernel410_ppc_altivec.h"



void 
nb_kernel410_ppc_altivec  (int *             p_nri,
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
	vector float ix,iy,iz,shvec;
	vector float vfacel,fs,fs2,nul;
	vector float dx,dy,dz;
	vector float Vvdwtot,vctot,qq,iq,c6,c12,VVc,FFc;
	vector float fix,fiy,fiz;
	vector float tmp1,tmp2,tmp3,tmp4;
	vector float rinv,r,rinvsq,rsq,rinvsix,Vvdw6,Vvdw12;
	vector float isai,isaj,isaprod,gbtsc,dvdasum,dvdaj,dvdatmp,gbscale,half;

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
	half=vec_half();
	vfacel=load_float_and_splat(p_facel);
	gbtsc=load_float_and_splat(p_gbtabscale);

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
			shvec      = load_xyz(shiftvec+is3);
			ii         = iinr[n];
			ii3        = 3*ii;
			ix         = load_xyz(pos+ii3);
			Vvdwtot     = nul;
			vctot      = nul;
			dvdasum    = nul;
			fix        = nul;
			fiy        = nul;
			fiz        = nul;
			ix         = vec_add(ix,shvec);
			nj0        = jindex[n];
			nj1        = jindex[n+1];
			splat_xyz_to_vectors(ix,&ix,&iy,&iz);
			ntiA       = 2*ntype*type[ii];
			iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
			isai       = load_float_and_splat(invsqrta+ii);

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
								 load_xyz(pos+j3d),&dx,&dy,&dz);
				dx              = vec_sub(ix,dx);
				dy              = vec_sub(iy,dy);
				dz              = vec_sub(iz,dz);
				rsq             = vec_madd(dx,dx,nul);
				rsq             = vec_madd(dy,dy,rsq);
				rsq             = vec_madd(dz,dz,rsq);
				rinv            = do_invsqrt(rsq);
				rinvsq          = vec_madd(rinv,rinv,nul);
				r               = vec_madd(rinv,rsq,nul);
				rinvsix         = vec_madd(rinvsq,rinvsq,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq,nul);
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 1/sqrt(a2) and multiply with 1/sqrt(a1) */
				isaj    = load_4_float(invsqrta+jnra,invsqrta+jnrb,
									   invsqrta+jnrc,invsqrta+jnrd);
				isaprod = vec_madd(isai,isaj,nul);
				qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
										   charge+jnrc,charge+jnrd),iq,nul);
				qq = vec_madd(isaprod,qq,nul);
				gbscale = vec_madd(isaprod,gbtsc,nul);
				load_4_pair(vdwparam+tja,vdwparam+tjb,vdwparam+tjc,vdwparam+tjd,&c6,&c12);
				do_4_ctable_coul(GBtab,vec_madd(r,gbscale,nul),&VVc,&FFc);
				dvdaj           = load_4_float(dvda+jnra,dvda+jnrb,
											   dvda+jnrc,dvda+jnrd);
				fs              = vec_madd(qq,FFc,nul);
				fs              = vec_madd(fs,gbscale,nul);   /* fijC */
				vctot           = vec_madd(qq,VVc,vctot);
				dvdatmp         = vec_madd(fs,r,nul);
				dvdatmp         = vec_madd(qq,VVc,dvdatmp);
				dvdasum         = vec_sub(dvdasum,dvdatmp);
				dvdaj 	      = vec_sub(dvdaj,dvdatmp);
				Vvdw6            = vec_madd(c6,rinvsix,nul);
				Vvdw12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),
										   nul);
				fs2             = vec_madd(vec_twelve(),Vvdw12,nul);
				fs2             = vec_nmsub(vec_six(),Vvdw6,fs2);
				fs              = vec_nmsub(fs2,rinv,fs);
				fs              = vec_nmsub(fs,rinv,nul);
				Vvdwtot          = vec_add(Vvdwtot,Vvdw12);
				Vvdwtot          = vec_sub(Vvdwtot,Vvdw6);
				store_4_float(dvdaj,dvda+jnra,dvda+jnrb,dvda+jnrc,dvda+jnrd);
				fix             = vec_madd(fs,dx,fix); /* +=fx */
				fiy             = vec_madd(fs,dy,fiy); /* +=fy */
				fiz             = vec_madd(fs,dz,fiz); /* +=fz */
				dx              = vec_nmsub(dx,fs,nul); /* -fx */
				dy              = vec_nmsub(dy,fs,nul); /* -fy */
				dz              = vec_nmsub(dz,fs,nul); /* -fz */
				transpose_3_to_4(dx,dy,dz,&tmp1,&tmp2,&tmp3,&tmp4);
				add_xyz_to_mem(faction+j3a,tmp1);
				add_xyz_to_mem(faction+j3b,tmp2);
				add_xyz_to_mem(faction+j3c,tmp3);
				add_xyz_to_mem(faction+j3d,tmp4);
			}
			if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				transpose_2_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),&dx,&dy,&dz);
				dx              = vec_sub(ix,dx);
				dy              = vec_sub(iy,dy);
				dz              = vec_sub(iz,dz);
				rsq             = vec_madd(dx,dx,nul);
				rsq             = vec_madd(dy,dy,rsq);
				rsq             = vec_madd(dz,dz,rsq);
				zero_highest_2_elements_in_vector(&rsq);
				rinv            = do_invsqrt(rsq);
				zero_highest_2_elements_in_vector(&rinv);
				rinvsq          = vec_madd(rinv,rinv,nul);
				r               = vec_madd(rinv,rsq,nul);
				rinvsix         = vec_madd(rinvsq,rinvsq,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq,nul);
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				/* load 1/sqrt(a2) and multiply with 1/sqrt(a1) */
				isaj    = load_2_float(invsqrta+jnra,invsqrta+jnrb);
				isaprod = vec_madd(isai,isaj,nul);
				qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
				qq = vec_madd(isaprod,qq,nul);
				gbscale = vec_madd(isaprod,gbtsc,nul);
				load_2_pair(vdwparam+tja,vdwparam+tjb,&c6,&c12);
				do_2_ctable_coul(GBtab,vec_madd(r,gbscale,nul),&VVc,&FFc);
				dvdaj           = load_2_float(dvda+jnra,dvda+jnrb);
				fs              = vec_madd(qq,FFc,nul);
				fs              = vec_madd(fs,gbscale,nul);   /* fijC */
				vctot           = vec_madd(qq,VVc,vctot);
				dvdatmp         = vec_madd(fs,r,nul);
				dvdatmp         = vec_madd(qq,VVc,dvdatmp);
				dvdasum         = vec_sub(dvdasum,dvdatmp);
				dvdaj			  = vec_sub(dvdaj,dvdatmp);
				Vvdw6            = vec_madd(c6,rinvsix,nul);
				Vvdw12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),
										   nul);
				fs2             = vec_madd(vec_twelve(),Vvdw12,nul);
				fs2             = vec_nmsub(vec_six(),Vvdw6,fs2);
				fs              = vec_nmsub(fs2,rinv,fs);
				fs              = vec_nmsub(fs,rinv,nul);
				Vvdwtot          = vec_add(Vvdwtot,Vvdw12);
				Vvdwtot          = vec_sub(Vvdwtot,Vvdw6);
				store_2_float(dvdaj,dvda+jnra,dvda+jnrb);
				fix             = vec_madd(fs,dx,fix); /* +=fx */
				fiy             = vec_madd(fs,dy,fiy); /* +=fy */
				fiz             = vec_madd(fs,dz,fiz); /* +=fz */
				dx              = vec_nmsub(dx,fs,nul); /* -fx */
				dy              = vec_nmsub(dy,fs,nul); /* -fy */
				dz              = vec_nmsub(dz,fs,nul); /* -fz */
				transpose_3_to_2(dx,dy,dz,&tmp1,&tmp2);
				add_xyz_to_mem(faction+j3a,tmp1);
				add_xyz_to_mem(faction+j3b,tmp2);
				k              += 2;
			}
			if((nj1-nj0) & 0x1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
				dx              = vec_sub(ix,dx);
				dy              = vec_sub(iy,dy);
				dz              = vec_sub(iz,dz);
				rsq             = vec_madd(dx,dx,nul);
				rsq             = vec_madd(dy,dy,rsq);
				rsq             = vec_madd(dz,dz,rsq);
				zero_highest_3_elements_in_vector(&rsq);
				rinv            = do_invsqrt(rsq);
				zero_highest_3_elements_in_vector(&rinv);
				rinvsq          = vec_madd(rinv,rinv,nul);
				r               = vec_madd(rinv,rsq,nul);
				rinvsix         = vec_madd(rinvsq,rinvsq,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq,nul);
				tja             = ntiA+2*type[jnra];
				/* load 1/sqrt(a2) and multiply with 1/sqrt(a1) */
				isaj    = load_1_float(invsqrta+jnra);
				isaprod = vec_madd(isai,isaj,nul);
				qq = vec_madd(load_1_float(charge+jnra),iq,nul);
				qq = vec_madd(isaprod,qq,nul);
				gbscale = vec_madd(isaprod,gbtsc,nul);
				load_1_pair(vdwparam+tja,&c6,&c12);
				do_1_ctable_coul(GBtab,vec_madd(r,gbscale,nul),&VVc,&FFc);
				dvdaj           = load_1_float(dvda+jnra);
				fs              = vec_madd(qq,FFc,nul);
				fs              = vec_madd(fs,gbscale,nul);   /* fijC */
				vctot           = vec_madd(qq,VVc,vctot);
				dvdatmp         = vec_madd(fs,r,nul);
				dvdatmp         = vec_madd(qq,VVc,dvdatmp);
				dvdasum         = vec_add(dvdasum,dvdatmp);
				dvdaj           = vec_sub(dvdaj,dvdatmp);
				Vvdw6            = vec_madd(c6,rinvsix,nul);
				Vvdw12           = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),
										   nul);
				fs2             = vec_madd(vec_twelve(),Vvdw12,nul);
				fs2             = vec_nmsub(vec_six(),Vvdw6,fs2);
				fs              = vec_nmsub(fs2,rinv,fs);
				fs              = vec_nmsub(fs,rinv,nul);
				Vvdwtot          = vec_add(Vvdwtot,Vvdw12);
				Vvdwtot          = vec_sub(Vvdwtot,Vvdw6);
				store_1_float(dvdaj,dvda+jnra);
				fix             = vec_madd(fs,dx,fix); /* +=fx */
				fiy             = vec_madd(fs,dy,fiy); /* +=fy */
				fiz             = vec_madd(fs,dz,fiz); /* +=fz */
				dx              = vec_nmsub(dx,fs,nul); /* -fx */
				dy              = vec_nmsub(dy,fs,nul); /* -fy */
				dz              = vec_nmsub(dz,fs,nul); /* -fz */
				transpose_3_to_1(dx,dy,dz,&tmp1);
				add_xyz_to_mem(faction+j3a,tmp1);
			}
			/* update outer data */
			transpose_3_to_4(fix,fiy,fiz,&tmp1,&tmp2,&tmp3,&tmp4);
			tmp1 = vec_add(tmp1,tmp3);
			tmp2 = vec_add(tmp2,tmp4);
			tmp1 = vec_add(tmp1,tmp2);

			add_xyz_to_mem(faction+ii3,tmp1);
			add_xyz_to_mem(fshift+is3,tmp1);

			add_vector_to_float(Vc+gid[n],vctot);
			add_vector_to_float(Vvdw+gid[n],Vvdwtot);

			add_vector_to_float(dvda+ii,dvdasum);
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
nb_kernel410nf_ppc_altivec(int *             p_nri,
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
	vector float ix,iy,iz,shvec;
	vector float vfacel,nul;
	vector float dx,dy,dz;
	vector float Vvdwtot,vctot,qq,iq,c6,c12,VVc;
	vector float rinv,r,rinvsq,rsq,rinvsix;
	vector float isai,isaprod,gbtsc,gbscale;

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
	vfacel=load_float_and_splat(p_facel);
	gbtsc=load_float_and_splat(p_gbtabscale);

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
			shvec      = load_xyz(shiftvec+is3);
			ii         = iinr[n];
			ii3        = 3*ii;
			ix         = load_xyz(pos+ii3);
			Vvdwtot     = nul;
			vctot      = nul;
			ix         = vec_add(ix,shvec);
			nj0        = jindex[n];
			nj1        = jindex[n+1];
			splat_xyz_to_vectors(ix,&ix,&iy,&iz);
			ntiA       = 2*ntype*type[ii];
			iq        = vec_madd(load_float_and_splat(charge+ii),vfacel,nul);
			isai       = load_float_and_splat(invsqrta+ii);

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
								 load_xyz(pos+j3d),&dx,&dy,&dz);
				dx              = vec_sub(ix,dx);
				dy              = vec_sub(iy,dy);
				dz              = vec_sub(iz,dz);
				rsq             = vec_madd(dx,dx,nul);
				rsq             = vec_madd(dy,dy,rsq);
				rsq             = vec_madd(dz,dz,rsq);
				rinv            = do_invsqrt(rsq);
				rinvsq          = vec_madd(rinv,rinv,nul);
				r               = vec_madd(rinv,rsq,nul);
				rinvsix         = vec_madd(rinvsq,rinvsq,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq,nul);
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				tjc             = ntiA+2*type[jnrc];
				tjd             = ntiA+2*type[jnrd];
				/* load 1/sqrt(a2) and multiply with 1/sqrt(a1) */
				isaprod = vec_madd(load_4_float(invsqrta+jnra,invsqrta+jnrb,
												invsqrta+jnrc,invsqrta+jnrd),
								   isai,nul);
				qq = vec_madd(load_4_float(charge+jnra,charge+jnrb,
										   charge+jnrc,charge+jnrd),iq,nul);
				qq = vec_madd(isaprod,qq,nul);
				gbscale = vec_madd(isaprod,gbtsc,nul);
				load_4_pair(vdwparam+tja,vdwparam+tjb,vdwparam+tjc,vdwparam+tjd,&c6,&c12);
				do_vonly_4_ctable_coul(GBtab,vec_madd(r,gbscale,nul),&VVc);
				vctot           = vec_madd(qq,VVc,vctot);
				Vvdwtot          = vec_nmsub(c6,rinvsix,Vvdwtot);
				Vvdwtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),
										   Vvdwtot);
			}
			if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				transpose_2_to_3(load_xyz(pos+j3a),
								 load_xyz(pos+j3b),&dx,&dy,&dz);
				dx              = vec_sub(ix,dx);
				dy              = vec_sub(iy,dy);
				dz              = vec_sub(iz,dz);
				rsq             = vec_madd(dx,dx,nul);
				rsq             = vec_madd(dy,dy,rsq);
				rsq             = vec_madd(dz,dz,rsq);
				zero_highest_2_elements_in_vector(&rsq);
				rinv            = do_invsqrt(rsq);
				zero_highest_2_elements_in_vector(&rinv);
				rinvsq          = vec_madd(rinv,rinv,nul);
				r               = vec_madd(rinv,rsq,nul);
				rinvsix         = vec_madd(rinvsq,rinvsq,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq,nul);
				tja             = ntiA+2*type[jnra];
				tjb             = ntiA+2*type[jnrb];
				/* load 1/sqrt(a2) and multiply with 1/sqrt(a1) */
				isaprod = vec_madd(load_2_float(invsqrta+jnra,invsqrta+jnrb),
								   isai,nul);
				qq = vec_madd(load_2_float(charge+jnra,charge+jnrb),iq,nul);
				qq = vec_madd(isaprod,qq,nul);
				gbscale = vec_madd(isaprod,gbtsc,nul);
				load_2_pair(vdwparam+tja,vdwparam+tjb,&c6,&c12);
				do_vonly_2_ctable_coul(GBtab,vec_madd(r,gbscale,nul),&VVc);
				vctot           = vec_madd(qq,VVc,vctot);
				Vvdwtot          = vec_nmsub(c6,rinvsix,Vvdwtot);
				Vvdwtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),
										   Vvdwtot);
				k              += 2;
			}
			if((nj1-nj0) & 0x1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				transpose_1_to_3(load_xyz(pos+j3a),&dx,&dy,&dz);
				dx              = vec_sub(ix,dx);
				dy              = vec_sub(iy,dy);
				dz              = vec_sub(iz,dz);
				rsq             = vec_madd(dx,dx,nul);
				rsq             = vec_madd(dy,dy,rsq);
				rsq             = vec_madd(dz,dz,rsq);
				zero_highest_3_elements_in_vector(&rsq);
				rinv            = do_invsqrt(rsq);
				zero_highest_3_elements_in_vector(&rinv);
				rinvsq          = vec_madd(rinv,rinv,nul);
				r               = vec_madd(rinv,rsq,nul);
				rinvsix         = vec_madd(rinvsq,rinvsq,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq,nul);
				tja             = ntiA+2*type[jnra];
				/* load 1/sqrt(a2) and multiply with 1/sqrt(a1) */
				isaprod = vec_madd(load_1_float(invsqrta+jnra),isai,nul);
				qq = vec_madd(load_1_float(charge+jnra),iq,nul);
				qq = vec_madd(isaprod,qq,nul);
				gbscale = vec_madd(isaprod,gbtsc,nul);
				load_1_pair(vdwparam+tja,&c6,&c12);
				do_vonly_1_ctable_coul(GBtab,vec_madd(r,gbscale,nul),&VVc);
				vctot           = vec_madd(qq,VVc,vctot);
				Vvdwtot          = vec_nmsub(c6,rinvsix,Vvdwtot);
				Vvdwtot          = vec_madd(c12,vec_madd(rinvsix,rinvsix,nul),
										   Vvdwtot);
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
