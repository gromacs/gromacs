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
#include "nb_kernel234_ppc_altivec.h"



void 
nb_kernel234_ppc_altivec  (int *             p_nri,
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
	vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3,ix4,iy4,iz4;
	vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4;

	vector float dx11,dy11,dz11;
	vector float dx22,dy22,dz22,dx23,dy23,dz23,dx24,dy24,dz24;
	vector float dx32,dy32,dz32,dx33,dy33,dz33,dx34,dy34,dz34;
	vector float dx42,dy42,dz42,dx43,dy43,dz43,dx44,dy44,dz44;

	vector float rsq11,rsq22,rsq23,rsq24,rsq32,rsq33,rsq34,rsq42,rsq43,rsq44;
	vector float rinv11,rinv22,rinv23,rinv24,rinv32,rinv33,rinv34,rinv42,rinv43,rinv44;
	vector float rinvsq22,rinvsq23,rinvsq24;
	vector float rinvsq32,rinvsq33,rinvsq34;
	vector float rinvsq42,rinvsq43,rinvsq44;
	vector float vc22,vc23,vc24,vc32,vc33,vc34,vc42,vc43,vc44;

	vector float vfacel,nul;
	vector float fs11,fs22,fs23,fs24,fs32,fs33,fs34,fs42,fs43,fs44;
	vector float fix1,fiy1,fiz1,fix2,fiy2,fiz2,fix3,fiy3,fiz3,fix4,fiy4,fiz4;
	vector float fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,fjx3,fjy3,fjz3,fjx4,fjy4,fjz4;
	vector float vctot,Vvdwtot,qqMM,qqMH,qqHH,qM,qH;
	vector float qqMMt,qqMHt,qqHHt; 
	vector float krsq22,krsq23,krsq24,krsq32,krsq33,krsq34;
	vector float krsq42,krsq43,krsq44;
	vector float vkrf,vcrf;
	vector float c6,c12,c6t,c12t;
	vector float VVd,VVr,FFd,FFr,tsc,r;

	int n,k,ii,is3,ii3,nj0,nj1,tp,tj;
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
	tsc=load_float_and_splat(p_tabscale);
	vfacel=load_float_and_splat(p_facel);
	ii        = iinr[0];
	vkrf=load_float_and_splat(p_krf);
	vcrf=load_float_and_splat(p_crf);
	qH        = load_float_and_splat(charge+iinr[0]+1);
	qM        = load_float_and_splat(charge+iinr[0]+3);
	qqMM      = vec_madd(qM,qM,nul);
	qqMH      = vec_madd(qM,qH,nul);
	qqHH      = vec_madd(qH,qH,nul);
	qqMM      = vec_madd(qqMM,vfacel,nul);
	qqMH      = vec_madd(qqMH,vfacel,nul);
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
			load_1_4atoms_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
										  &ix2,&iy2,&iz2,&ix3,&iy3,&iz3,
										  &ix4,&iy4,&iz4);
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
			fix4       = nul;
			fiy4       = nul;
			fiz4       = nul;
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
				load_4_4atoms(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3,
							  &jx4,&jy4,&jz4);

				dx11            = vec_sub(ix1,jx1);
				dy11            = vec_sub(iy1,jy1);
				dz11            = vec_sub(iz1,jz1);
				dx22            = vec_sub(ix2,jx2);
				dy22            = vec_sub(iy2,jy2);
				dz22            = vec_sub(iz2,jz2);
				dx23            = vec_sub(ix2,jx3);
				dy23            = vec_sub(iy2,jy3);
				dz23            = vec_sub(iz2,jz3);
				dx24            = vec_sub(ix2,jx4);
				dy24            = vec_sub(iy2,jy4);
				dz24            = vec_sub(iz2,jz4);
				dx32            = vec_sub(ix3,jx2);
				dy32            = vec_sub(iy3,jy2);
				dz32            = vec_sub(iz3,jz2);
				dx33            = vec_sub(ix3,jx3);
				dy33            = vec_sub(iy3,jy3);
				dz33            = vec_sub(iz3,jz3);
				dx34            = vec_sub(ix3,jx4);
				dy34            = vec_sub(iy3,jy4);
				dz34            = vec_sub(iz3,jz4);
				dx42            = vec_sub(ix4,jx2);
				dy42            = vec_sub(iy4,jy2);
				dz42            = vec_sub(iz4,jz2);
				dx43            = vec_sub(ix4,jx3);
				dy43            = vec_sub(iy4,jy3);
				dz43            = vec_sub(iz4,jz3);
				dx44            = vec_sub(ix4,jx4);
				dy44            = vec_sub(iy4,jy4);
				dz44            = vec_sub(iz4,jz4);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq24           = vec_madd(dx24,dx24,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq34           = vec_madd(dx34,dx34,nul);
				rsq42           = vec_madd(dx42,dx42,nul);
				rsq43           = vec_madd(dx43,dx43,nul);
				rsq44           = vec_madd(dx44,dx44,nul);

				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq24           = vec_madd(dy24,dy24,rsq24);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq34           = vec_madd(dy34,dy34,rsq34);
				rsq42           = vec_madd(dy42,dy42,rsq42);
				rsq43           = vec_madd(dy43,dy43,rsq43);
				rsq44           = vec_madd(dy44,dy44,rsq44);
     
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq24           = vec_madd(dz24,dz24,rsq24);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);
				rsq34           = vec_madd(dz34,dz34,rsq34);
				rsq42           = vec_madd(dz42,dz42,rsq42);
				rsq43           = vec_madd(dz43,dz43,rsq43);
				rsq44           = vec_madd(dz44,dz44,rsq44);

				rinv11          = do_invsqrt(rsq11);
				do_9_invsqrt(rsq22,rsq23,rsq24,
							 rsq32,rsq33,rsq34,
							 rsq42,rsq43,rsq44,
							 &rinv22,&rinv23,&rinv24,
							 &rinv32,&rinv33,&rinv34,
							 &rinv42,&rinv43,&rinv44);
      
				r               = vec_madd(rinv11,rsq11,nul);
				do_4_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&FFd,&VVr,&FFr);
				Vvdwtot         = vec_madd(c6,VVd,Vvdwtot);
				fs11            = vec_nmsub(c6,FFd,nul);
				Vvdwtot         = vec_madd(c12,VVr,Vvdwtot);
				fs11            = vec_nmsub(c12,FFr,fs11);
				fs11            = vec_madd(fs11,tsc,nul);
				
				krsq22          = vec_madd(vkrf,rsq22,nul);
				krsq23          = vec_madd(vkrf,rsq23,nul);
				krsq24          = vec_madd(vkrf,rsq24,nul);
				krsq32          = vec_madd(vkrf,rsq32,nul);
				krsq33          = vec_madd(vkrf,rsq33,nul);
				krsq34          = vec_madd(vkrf,rsq34,nul);
				krsq42          = vec_madd(vkrf,rsq42,nul);
				krsq43          = vec_madd(vkrf,rsq43,nul);
				krsq44          = vec_madd(vkrf,rsq44,nul);

				rinvsq22        = vec_madd(rinv22,rinv22,nul);
				rinvsq23        = vec_madd(rinv23,rinv23,nul);
				rinvsq24        = vec_madd(rinv24,rinv24,nul);
				rinvsq32        = vec_madd(rinv32,rinv32,nul);
				rinvsq33        = vec_madd(rinv33,rinv33,nul);
				rinvsq34        = vec_madd(rinv34,rinv34,nul);
				rinvsq42        = vec_madd(rinv42,rinv42,nul);
				rinvsq43        = vec_madd(rinv43,rinv43,nul);
				rinvsq44        = vec_madd(rinv44,rinv44,nul);

				vc22            = vec_add(rinv22,krsq22);
				vc23            = vec_add(rinv23,krsq23);
				vc24            = vec_add(rinv24,krsq24);
				vc32            = vec_add(rinv32,krsq32);
				vc33            = vec_add(rinv33,krsq33);
				vc34            = vec_add(rinv34,krsq34);
				vc42            = vec_add(rinv42,krsq42);
				vc43            = vec_add(rinv43,krsq43);
				vc44            = vec_add(rinv44,krsq44);

				vc22            = vec_sub(vc22,vcrf);
				vc23            = vec_sub(vc23,vcrf);
				vc24            = vec_sub(vc24,vcrf);
				vc32            = vec_sub(vc32,vcrf);
				vc33            = vec_sub(vc33,vcrf);
				vc34            = vec_sub(vc34,vcrf);
				vc42            = vec_sub(vc42,vcrf);
				vc43            = vec_sub(vc43,vcrf);
				vc44            = vec_sub(vc44,vcrf);

				fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
				fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
				fs24            = vec_nmsub(vec_two(),krsq24,rinv24);
				fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
				fs33            = vec_nmsub(vec_two(),krsq33,rinv33);
				fs34            = vec_nmsub(vec_two(),krsq34,rinv34);
				fs42            = vec_nmsub(vec_two(),krsq42,rinv42);
				fs43            = vec_nmsub(vec_two(),krsq43,rinv43);
				fs44            = vec_nmsub(vec_two(),krsq44,rinv44);

				fs22            = vec_madd(fs22,qqHH,nul);
				fs23            = vec_madd(fs23,qqHH,nul);
				fs24            = vec_madd(fs24,qqMH,nul);
				fs32            = vec_madd(fs32,qqHH,nul);
				fs33            = vec_madd(fs33,qqHH,nul);
				fs34            = vec_madd(fs34,qqMH,nul);
				fs42            = vec_madd(fs42,qqMH,nul);
				fs43            = vec_madd(fs43,qqMH,nul);
				fs44            = vec_madd(fs44,qqMM,nul);

				fs22            = vec_madd(fs22,rinvsq22,nul);
				fs23            = vec_madd(fs23,rinvsq23,nul);
				fs24            = vec_madd(fs24,rinvsq24,nul);
				fs32            = vec_madd(fs32,rinvsq32,nul);
				fs33            = vec_madd(fs33,rinvsq33,nul);
				fs34            = vec_madd(fs34,rinvsq34,nul);
				fs42            = vec_madd(fs42,rinvsq42,nul);
				fs43            = vec_madd(fs43,rinvsq43,nul);
				fs44            = vec_madd(fs44,rinvsq44,nul);
      
				vctot           = vec_madd(qqHH,vc22,vctot);
				vctot           = vec_madd(qqHH,vc23,vctot);
				vctot           = vec_madd(qqMH,vc24,vctot);
				vctot           = vec_madd(qqHH,vc32,vctot);
				vctot           = vec_madd(qqHH,vc33,vctot);
				vctot           = vec_madd(qqMH,vc34,vctot);
				vctot           = vec_madd(qqMH,vc42,vctot);
				vctot           = vec_madd(qqMH,vc43,vctot);
				vctot           = vec_madd(qqMM,vc44,vctot);

				fs11            = vec_madd(fs11,rinv11,nul);

				fix1            = vec_madd(fs11,dx11,fix1);
				fiy1            = vec_madd(fs11,dy11,fiy1);
				fiz1            = vec_madd(fs11,dz11,fiz1);

				fix2            = vec_madd(fs22,dx22,fix2);
				fiy2            = vec_madd(fs22,dy22,fiy2);
				fiz2            = vec_madd(fs22,dz22,fiz2);
				fix3            = vec_madd(fs32,dx32,fix3);
				fiy3            = vec_madd(fs32,dy32,fiy3);
				fiz3            = vec_madd(fs32,dz32,fiz3);
				fix4            = vec_madd(fs42,dx42,fix4);
				fiy4            = vec_madd(fs42,dy42,fiy4);
				fiz4            = vec_madd(fs42,dz42,fiz4);

				fix2            = vec_madd(fs23,dx23,fix2);
				fiy2            = vec_madd(fs23,dy23,fiy2);
				fiz2            = vec_madd(fs23,dz23,fiz2);
				fix3            = vec_madd(fs33,dx33,fix3);
				fiy3            = vec_madd(fs33,dy33,fiy3);
				fiz3            = vec_madd(fs33,dz33,fiz3);
				fix4            = vec_madd(fs43,dx43,fix4);
				fiy4            = vec_madd(fs43,dy43,fiy4);
				fiz4            = vec_madd(fs43,dz43,fiz4);

				fix2            = vec_madd(fs24,dx24,fix2);
				fiy2            = vec_madd(fs24,dy24,fiy2);
				fiz2            = vec_madd(fs24,dz24,fiz2);
				fix3            = vec_madd(fs34,dx34,fix3);
				fiy3            = vec_madd(fs34,dy34,fiy3);
				fiz3            = vec_madd(fs34,dz34,fiz3);
				fix4            = vec_madd(fs44,dx44,fix4);
				fiy4            = vec_madd(fs44,dy44,fiy4);
				fiz4            = vec_madd(fs44,dz44,fiz4);


				fjx1            = vec_nmsub(fs11,dx11,nul);
				fjy1            = vec_nmsub(fs11,dy11,nul);
				fjz1            = vec_nmsub(fs11,dz11,nul);
				fjx2            = vec_nmsub(fs22,dx22,nul);
				fjy2            = vec_nmsub(fs22,dy22,nul);
				fjz2            = vec_nmsub(fs22,dz22,nul);
				fjx3            = vec_nmsub(fs23,dx23,nul);
				fjy3            = vec_nmsub(fs23,dy23,nul);
				fjz3            = vec_nmsub(fs23,dz23,nul);
				fjx4            = vec_nmsub(fs24,dx24,nul);
				fjy4            = vec_nmsub(fs24,dy24,nul);
				fjz4            = vec_nmsub(fs24,dz24,nul);

				fjx2            = vec_nmsub(fs32,dx32,fjx2);
				fjy2            = vec_nmsub(fs32,dy32,fjy2);
				fjz2            = vec_nmsub(fs32,dz32,fjz2);
				fjx3            = vec_nmsub(fs33,dx33,fjx3);
				fjy3            = vec_nmsub(fs33,dy33,fjy3);
				fjz3            = vec_nmsub(fs33,dz33,fjz3);
				fjx4            = vec_nmsub(fs34,dx34,fjx4);
				fjy4            = vec_nmsub(fs34,dy34,fjy4);
				fjz4            = vec_nmsub(fs34,dz34,fjz4);

				fjx2            = vec_nmsub(fs42,dx42,fjx2);
				fjy2            = vec_nmsub(fs42,dy42,fjy2);
				fjz2            = vec_nmsub(fs42,dz42,fjz2);
				fjx3            = vec_nmsub(fs43,dx43,fjx3);
				fjy3            = vec_nmsub(fs43,dy43,fjy3);
				fjz3            = vec_nmsub(fs43,dz43,fjz3);
				fjx4            = vec_nmsub(fs44,dx44,fjx4);
				fjy4            = vec_nmsub(fs44,dy44,fjy4);
				fjz4            = vec_nmsub(fs44,dz44,fjz4);

				add_force_to_4_4atoms(faction+j3a,faction+j3b,
									  faction+j3c,faction+j3d,
									  fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,
									  fjx3,fjy3,fjz3,fjx4,fjy4,fjz4);
			} 
			if(k<(nj1-2)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				load_3_4atoms(pos+j3a,pos+j3b,pos+j3c,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3,
							  &jx4,&jy4,&jz4);

				qqMMt           = vec_sld(qqMM,nul,4);
				qqMHt           = vec_sld(qqMH,nul,4);
				qqHHt           = vec_sld(qqHH,nul,4);
				c6t             = vec_sld(c6,nul,4);
				c12t            = vec_sld(c12,nul,4);

				dx11            = vec_sub(ix1,jx1);
				dy11            = vec_sub(iy1,jy1);
				dz11            = vec_sub(iz1,jz1);
				dx22            = vec_sub(ix2,jx2);
				dy22            = vec_sub(iy2,jy2);
				dz22            = vec_sub(iz2,jz2);
				dx23            = vec_sub(ix2,jx3);
				dy23            = vec_sub(iy2,jy3);
				dz23            = vec_sub(iz2,jz3);
				dx24            = vec_sub(ix2,jx4);
				dy24            = vec_sub(iy2,jy4);
				dz24            = vec_sub(iz2,jz4);
				dx32            = vec_sub(ix3,jx2);
				dy32            = vec_sub(iy3,jy2);
				dz32            = vec_sub(iz3,jz2);
				dx33            = vec_sub(ix3,jx3);
				dy33            = vec_sub(iy3,jy3);
				dz33            = vec_sub(iz3,jz3);
				dx34            = vec_sub(ix3,jx4);
				dy34            = vec_sub(iy3,jy4);
				dz34            = vec_sub(iz3,jz4);
				dx42            = vec_sub(ix4,jx2);
				dy42            = vec_sub(iy4,jy2);
				dz42            = vec_sub(iz4,jz2);
				dx43            = vec_sub(ix4,jx3);
				dy43            = vec_sub(iy4,jy3);
				dz43            = vec_sub(iz4,jz3);
				dx44            = vec_sub(ix4,jx4);
				dy44            = vec_sub(iy4,jy4);
				dz44            = vec_sub(iz4,jz4);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq24           = vec_madd(dx24,dx24,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq34           = vec_madd(dx34,dx34,nul);
				rsq42           = vec_madd(dx42,dx42,nul);
				rsq43           = vec_madd(dx43,dx43,nul);
				rsq44           = vec_madd(dx44,dx44,nul);

				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq24           = vec_madd(dy24,dy24,rsq24);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq34           = vec_madd(dy34,dy34,rsq34);
				rsq42           = vec_madd(dy42,dy42,rsq42);
				rsq43           = vec_madd(dy43,dy43,rsq43);
				rsq44           = vec_madd(dy44,dy44,rsq44);
     
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq24           = vec_madd(dz24,dz24,rsq24);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);
				rsq34           = vec_madd(dz34,dz34,rsq34);
				rsq42           = vec_madd(dz42,dz42,rsq42);
				rsq43           = vec_madd(dz43,dz43,rsq43);
				rsq44           = vec_madd(dz44,dz44,rsq44);

				rinv11          = do_invsqrt(rsq11);
				do_9_invsqrt(rsq22,rsq23,rsq24,
							 rsq32,rsq33,rsq34,
							 rsq42,rsq43,rsq44,
							 &rinv22,&rinv23,&rinv24,
							 &rinv32,&rinv33,&rinv34,
							 &rinv42,&rinv43,&rinv44);

				zero_highest_element_in_vector(&rinv11);
				zero_highest_element_in_9_vectors(&rinv22,&rinv23,&rinv24,
												  &rinv32,&rinv33,&rinv34,
												  &rinv42,&rinv43,&rinv44);
      
      
				r               = vec_madd(rinv11,rsq11,nul);
				do_3_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&FFd,&VVr,&FFr);
				Vvdwtot         = vec_madd(c6t,VVd,Vvdwtot);
				fs11            = vec_nmsub(c6t,FFd,nul);
				Vvdwtot         = vec_madd(c12t,VVr,Vvdwtot);
				fs11            = vec_nmsub(c12t,FFr,fs11);
				fs11            = vec_madd(fs11,tsc,nul);
				
				krsq22          = vec_madd(vkrf,rsq22,nul);
				krsq23          = vec_madd(vkrf,rsq23,nul);
				krsq24          = vec_madd(vkrf,rsq24,nul);
				krsq32          = vec_madd(vkrf,rsq32,nul);
				krsq33          = vec_madd(vkrf,rsq33,nul);
				krsq34          = vec_madd(vkrf,rsq34,nul);
				krsq42          = vec_madd(vkrf,rsq42,nul);
				krsq43          = vec_madd(vkrf,rsq43,nul);
				krsq44          = vec_madd(vkrf,rsq44,nul);


				rinvsq22        = vec_madd(rinv22,rinv22,nul);
				rinvsq23        = vec_madd(rinv23,rinv23,nul);
				rinvsq24        = vec_madd(rinv24,rinv24,nul);
				rinvsq32        = vec_madd(rinv32,rinv32,nul);
				rinvsq33        = vec_madd(rinv33,rinv33,nul);
				rinvsq34        = vec_madd(rinv34,rinv34,nul);
				rinvsq42        = vec_madd(rinv42,rinv42,nul);
				rinvsq43        = vec_madd(rinv43,rinv43,nul);
				rinvsq44        = vec_madd(rinv44,rinv44,nul);

				vc22            = vec_add(rinv22,krsq22);
				vc23            = vec_add(rinv23,krsq23);
				vc24            = vec_add(rinv24,krsq24);
				vc32            = vec_add(rinv32,krsq32);
				vc33            = vec_add(rinv33,krsq33);
				vc34            = vec_add(rinv34,krsq34);
				vc42            = vec_add(rinv42,krsq42);
				vc43            = vec_add(rinv43,krsq43);
				vc44            = vec_add(rinv44,krsq44);

				vc22            = vec_sub(vc22,vcrf);
				vc23            = vec_sub(vc23,vcrf);
				vc24            = vec_sub(vc24,vcrf);
				vc32            = vec_sub(vc32,vcrf);
				vc33            = vec_sub(vc33,vcrf);
				vc34            = vec_sub(vc34,vcrf);
				vc42            = vec_sub(vc42,vcrf);
				vc43            = vec_sub(vc43,vcrf);
				vc44            = vec_sub(vc44,vcrf);

				fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
				fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
				fs24            = vec_nmsub(vec_two(),krsq24,rinv24);
				fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
				fs33            = vec_nmsub(vec_two(),krsq33,rinv33);
				fs34            = vec_nmsub(vec_two(),krsq34,rinv34);
				fs42            = vec_nmsub(vec_two(),krsq42,rinv42);
				fs43            = vec_nmsub(vec_two(),krsq43,rinv43);
				fs44            = vec_nmsub(vec_two(),krsq44,rinv44);

				fs22            = vec_madd(fs22,qqHHt,nul);
				fs23            = vec_madd(fs23,qqHHt,nul);
				fs24            = vec_madd(fs24,qqMHt,nul);
				fs32            = vec_madd(fs32,qqHHt,nul);
				fs33            = vec_madd(fs33,qqHHt,nul);
				fs34            = vec_madd(fs34,qqMHt,nul);
				fs42            = vec_madd(fs42,qqMHt,nul);
				fs43            = vec_madd(fs43,qqMHt,nul);
				fs44            = vec_madd(fs44,qqMMt,nul);

				fs22            = vec_madd(fs22,rinvsq22,nul);
				fs23            = vec_madd(fs23,rinvsq23,nul);
				fs24            = vec_madd(fs24,rinvsq24,nul);
				fs32            = vec_madd(fs32,rinvsq32,nul);
				fs33            = vec_madd(fs33,rinvsq33,nul);
				fs34            = vec_madd(fs34,rinvsq34,nul);
				fs42            = vec_madd(fs42,rinvsq42,nul);
				fs43            = vec_madd(fs43,rinvsq43,nul);
				fs44            = vec_madd(fs44,rinvsq44,nul);
      
				vctot           = vec_madd(qqHHt,vc22,vctot);
				vctot           = vec_madd(qqHHt,vc23,vctot);
				vctot           = vec_madd(qqMHt,vc24,vctot);
				vctot           = vec_madd(qqHHt,vc32,vctot);
				vctot           = vec_madd(qqHHt,vc33,vctot);
				vctot           = vec_madd(qqMHt,vc34,vctot);
				vctot           = vec_madd(qqMHt,vc42,vctot);
				vctot           = vec_madd(qqMHt,vc43,vctot);
				vctot           = vec_madd(qqMMt,vc44,vctot);

				fs11            = vec_madd(fs11,rinv11,nul);

				fix1            = vec_madd(fs11,dx11,fix1);
				fiy1            = vec_madd(fs11,dy11,fiy1);
				fiz1            = vec_madd(fs11,dz11,fiz1);

				fix2            = vec_madd(fs22,dx22,fix2);
				fiy2            = vec_madd(fs22,dy22,fiy2);
				fiz2            = vec_madd(fs22,dz22,fiz2);
				fix3            = vec_madd(fs32,dx32,fix3);
				fiy3            = vec_madd(fs32,dy32,fiy3);
				fiz3            = vec_madd(fs32,dz32,fiz3);
				fix4            = vec_madd(fs42,dx42,fix4);
				fiy4            = vec_madd(fs42,dy42,fiy4);
				fiz4            = vec_madd(fs42,dz42,fiz4);

				fix2            = vec_madd(fs23,dx23,fix2);
				fiy2            = vec_madd(fs23,dy23,fiy2);
				fiz2            = vec_madd(fs23,dz23,fiz2);
				fix3            = vec_madd(fs33,dx33,fix3);
				fiy3            = vec_madd(fs33,dy33,fiy3);
				fiz3            = vec_madd(fs33,dz33,fiz3);
				fix4            = vec_madd(fs43,dx43,fix4);
				fiy4            = vec_madd(fs43,dy43,fiy4);
				fiz4            = vec_madd(fs43,dz43,fiz4);

				fix2            = vec_madd(fs24,dx24,fix2);
				fiy2            = vec_madd(fs24,dy24,fiy2);
				fiz2            = vec_madd(fs24,dz24,fiz2);
				fix3            = vec_madd(fs34,dx34,fix3);
				fiy3            = vec_madd(fs34,dy34,fiy3);
				fiz3            = vec_madd(fs34,dz34,fiz3);
				fix4            = vec_madd(fs44,dx44,fix4);
				fiy4            = vec_madd(fs44,dy44,fiy4);
				fiz4            = vec_madd(fs44,dz44,fiz4);


				fjx1            = vec_nmsub(fs11,dx11,nul);
				fjy1            = vec_nmsub(fs11,dy11,nul);
				fjz1            = vec_nmsub(fs11,dz11,nul);
				fjx2            = vec_nmsub(fs22,dx22,nul);
				fjy2            = vec_nmsub(fs22,dy22,nul);
				fjz2            = vec_nmsub(fs22,dz22,nul);
				fjx3            = vec_nmsub(fs23,dx23,nul);
				fjy3            = vec_nmsub(fs23,dy23,nul);
				fjz3            = vec_nmsub(fs23,dz23,nul);
				fjx4            = vec_nmsub(fs24,dx24,nul);
				fjy4            = vec_nmsub(fs24,dy24,nul);
				fjz4            = vec_nmsub(fs24,dz24,nul);

				fjx2            = vec_nmsub(fs32,dx32,fjx2);
				fjy2            = vec_nmsub(fs32,dy32,fjy2);
				fjz2            = vec_nmsub(fs32,dz32,fjz2);
				fjx3            = vec_nmsub(fs33,dx33,fjx3);
				fjy3            = vec_nmsub(fs33,dy33,fjy3);
				fjz3            = vec_nmsub(fs33,dz33,fjz3);
				fjx4            = vec_nmsub(fs34,dx34,fjx4);
				fjy4            = vec_nmsub(fs34,dy34,fjy4);
				fjz4            = vec_nmsub(fs34,dz34,fjz4);

				fjx2            = vec_nmsub(fs42,dx42,fjx2);
				fjy2            = vec_nmsub(fs42,dy42,fjy2);
				fjz2            = vec_nmsub(fs42,dz42,fjz2);
				fjx3            = vec_nmsub(fs43,dx43,fjx3);
				fjy3            = vec_nmsub(fs43,dy43,fjy3);
				fjz3            = vec_nmsub(fs43,dz43,fjz3);
				fjx4            = vec_nmsub(fs44,dx44,fjx4);
				fjy4            = vec_nmsub(fs44,dy44,fjy4);
				fjz4            = vec_nmsub(fs44,dz44,fjz4);

				add_force_to_3_4atoms(faction+j3a,faction+j3b,faction+j3c,
									  fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,
									  fjx3,fjy3,fjz3,fjx4,fjy4,fjz4);
			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				load_2_4atoms(pos+j3a,pos+j3b,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3,
							  &jx4,&jy4,&jz4);

				qqMMt           = vec_sld(qqMM,nul,8);
				qqMHt           = vec_sld(qqMH,nul,8);
				qqHHt           = vec_sld(qqHH,nul,8);
				c6t             = vec_sld(c6,nul,8);
				c12t            = vec_sld(c12,nul,8);

				dx11            = vec_sub(ix1,jx1);
				dy11            = vec_sub(iy1,jy1);
				dz11            = vec_sub(iz1,jz1);
				dx22            = vec_sub(ix2,jx2);
				dy22            = vec_sub(iy2,jy2);
				dz22            = vec_sub(iz2,jz2);
				dx23            = vec_sub(ix2,jx3);
				dy23            = vec_sub(iy2,jy3);
				dz23            = vec_sub(iz2,jz3);
				dx24            = vec_sub(ix2,jx4);
				dy24            = vec_sub(iy2,jy4);
				dz24            = vec_sub(iz2,jz4);
				dx32            = vec_sub(ix3,jx2);
				dy32            = vec_sub(iy3,jy2);
				dz32            = vec_sub(iz3,jz2);
				dx33            = vec_sub(ix3,jx3);
				dy33            = vec_sub(iy3,jy3);
				dz33            = vec_sub(iz3,jz3);
				dx34            = vec_sub(ix3,jx4);
				dy34            = vec_sub(iy3,jy4);
				dz34            = vec_sub(iz3,jz4);
				dx42            = vec_sub(ix4,jx2);
				dy42            = vec_sub(iy4,jy2);
				dz42            = vec_sub(iz4,jz2);
				dx43            = vec_sub(ix4,jx3);
				dy43            = vec_sub(iy4,jy3);
				dz43            = vec_sub(iz4,jz3);
				dx44            = vec_sub(ix4,jx4);
				dy44            = vec_sub(iy4,jy4);
				dz44            = vec_sub(iz4,jz4);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq24           = vec_madd(dx24,dx24,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq34           = vec_madd(dx34,dx34,nul);
				rsq42           = vec_madd(dx42,dx42,nul);
				rsq43           = vec_madd(dx43,dx43,nul);
				rsq44           = vec_madd(dx44,dx44,nul);

				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq24           = vec_madd(dy24,dy24,rsq24);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq34           = vec_madd(dy34,dy34,rsq34);
				rsq42           = vec_madd(dy42,dy42,rsq42);
				rsq43           = vec_madd(dy43,dy43,rsq43);
				rsq44           = vec_madd(dy44,dy44,rsq44);
     
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq24           = vec_madd(dz24,dz24,rsq24);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);
				rsq34           = vec_madd(dz34,dz34,rsq34);
				rsq42           = vec_madd(dz42,dz42,rsq42);
				rsq43           = vec_madd(dz43,dz43,rsq43);
				rsq44           = vec_madd(dz44,dz44,rsq44);

				rinv11          = do_invsqrt(rsq11);
				do_9_invsqrt(rsq22,rsq23,rsq24,
							 rsq32,rsq33,rsq34,
							 rsq42,rsq43,rsq44,
							 &rinv22,&rinv23,&rinv24,
							 &rinv32,&rinv33,&rinv34,
							 &rinv42,&rinv43,&rinv44);

				zero_highest_2_elements_in_vector(&rinv11);
				zero_highest_2_elements_in_9_vectors(&rinv22,&rinv23,&rinv24,
													 &rinv32,&rinv33,&rinv34,
													 &rinv42,&rinv43,&rinv44);

      
				r               = vec_madd(rinv11,rsq11,nul);
				do_2_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&FFd,&VVr,&FFr);
				Vvdwtot         = vec_madd(c6t,VVd,Vvdwtot);
				fs11            = vec_nmsub(c6t,FFd,nul);
				Vvdwtot         = vec_madd(c12t,VVr,Vvdwtot);
				fs11            = vec_nmsub(c12t,FFr,fs11);
				fs11            = vec_madd(fs11,tsc,nul);
				
				krsq22          = vec_madd(vkrf,rsq22,nul);
				krsq23          = vec_madd(vkrf,rsq23,nul);
				krsq24          = vec_madd(vkrf,rsq24,nul);
				krsq32          = vec_madd(vkrf,rsq32,nul);
				krsq33          = vec_madd(vkrf,rsq33,nul);
				krsq34          = vec_madd(vkrf,rsq34,nul);
				krsq42          = vec_madd(vkrf,rsq42,nul);
				krsq43          = vec_madd(vkrf,rsq43,nul);
				krsq44          = vec_madd(vkrf,rsq44,nul);

				rinvsq22        = vec_madd(rinv22,rinv22,nul);
				rinvsq23        = vec_madd(rinv23,rinv23,nul);
				rinvsq24        = vec_madd(rinv24,rinv24,nul);
				rinvsq32        = vec_madd(rinv32,rinv32,nul);
				rinvsq33        = vec_madd(rinv33,rinv33,nul);
				rinvsq34        = vec_madd(rinv34,rinv34,nul);
				rinvsq42        = vec_madd(rinv42,rinv42,nul);
				rinvsq43        = vec_madd(rinv43,rinv43,nul);
				rinvsq44        = vec_madd(rinv44,rinv44,nul);

				vc22            = vec_add(rinv22,krsq22);
				vc23            = vec_add(rinv23,krsq23);
				vc24            = vec_add(rinv24,krsq24);
				vc32            = vec_add(rinv32,krsq32);
				vc33            = vec_add(rinv33,krsq33);
				vc34            = vec_add(rinv34,krsq34);
				vc42            = vec_add(rinv42,krsq42);
				vc43            = vec_add(rinv43,krsq43);
				vc44            = vec_add(rinv44,krsq44);

				vc22            = vec_sub(vc22,vcrf);
				vc23            = vec_sub(vc23,vcrf);
				vc24            = vec_sub(vc24,vcrf);
				vc32            = vec_sub(vc32,vcrf);
				vc33            = vec_sub(vc33,vcrf);
				vc34            = vec_sub(vc34,vcrf);
				vc42            = vec_sub(vc42,vcrf);
				vc43            = vec_sub(vc43,vcrf);
				vc44            = vec_sub(vc44,vcrf);

				fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
				fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
				fs24            = vec_nmsub(vec_two(),krsq24,rinv24);
				fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
				fs33            = vec_nmsub(vec_two(),krsq33,rinv33);
				fs34            = vec_nmsub(vec_two(),krsq34,rinv34);
				fs42            = vec_nmsub(vec_two(),krsq42,rinv42);
				fs43            = vec_nmsub(vec_two(),krsq43,rinv43);
				fs44            = vec_nmsub(vec_two(),krsq44,rinv44);

				fs22            = vec_madd(fs22,qqHHt,nul);
				fs23            = vec_madd(fs23,qqHHt,nul);
				fs24            = vec_madd(fs24,qqMHt,nul);
				fs32            = vec_madd(fs32,qqHHt,nul);
				fs33            = vec_madd(fs33,qqHHt,nul);
				fs34            = vec_madd(fs34,qqMHt,nul);
				fs42            = vec_madd(fs42,qqMHt,nul);
				fs43            = vec_madd(fs43,qqMHt,nul);
				fs44            = vec_madd(fs44,qqMMt,nul);

				fs22            = vec_madd(fs22,rinvsq22,nul);
				fs23            = vec_madd(fs23,rinvsq23,nul);
				fs24            = vec_madd(fs24,rinvsq24,nul);
				fs32            = vec_madd(fs32,rinvsq32,nul);
				fs33            = vec_madd(fs33,rinvsq33,nul);
				fs34            = vec_madd(fs34,rinvsq34,nul);
				fs42            = vec_madd(fs42,rinvsq42,nul);
				fs43            = vec_madd(fs43,rinvsq43,nul);
				fs44            = vec_madd(fs44,rinvsq44,nul);
      
				vctot           = vec_madd(qqHHt,vc22,vctot);
				vctot           = vec_madd(qqHHt,vc23,vctot);
				vctot           = vec_madd(qqMHt,vc24,vctot);
				vctot           = vec_madd(qqHHt,vc32,vctot);
				vctot           = vec_madd(qqHHt,vc33,vctot);
				vctot           = vec_madd(qqMHt,vc34,vctot);
				vctot           = vec_madd(qqMHt,vc42,vctot);
				vctot           = vec_madd(qqMHt,vc43,vctot);
				vctot           = vec_madd(qqMMt,vc44,vctot);

				fs11            = vec_madd(fs11,rinv11,nul);

				fix1            = vec_madd(fs11,dx11,fix1);
				fiy1            = vec_madd(fs11,dy11,fiy1);
				fiz1            = vec_madd(fs11,dz11,fiz1);

				fix2            = vec_madd(fs22,dx22,fix2);
				fiy2            = vec_madd(fs22,dy22,fiy2);
				fiz2            = vec_madd(fs22,dz22,fiz2);
				fix3            = vec_madd(fs32,dx32,fix3);
				fiy3            = vec_madd(fs32,dy32,fiy3);
				fiz3            = vec_madd(fs32,dz32,fiz3);
				fix4            = vec_madd(fs42,dx42,fix4);
				fiy4            = vec_madd(fs42,dy42,fiy4);
				fiz4            = vec_madd(fs42,dz42,fiz4);

				fix2            = vec_madd(fs23,dx23,fix2);
				fiy2            = vec_madd(fs23,dy23,fiy2);
				fiz2            = vec_madd(fs23,dz23,fiz2);
				fix3            = vec_madd(fs33,dx33,fix3);
				fiy3            = vec_madd(fs33,dy33,fiy3);
				fiz3            = vec_madd(fs33,dz33,fiz3);
				fix4            = vec_madd(fs43,dx43,fix4);
				fiy4            = vec_madd(fs43,dy43,fiy4);
				fiz4            = vec_madd(fs43,dz43,fiz4);

				fix2            = vec_madd(fs24,dx24,fix2);
				fiy2            = vec_madd(fs24,dy24,fiy2);
				fiz2            = vec_madd(fs24,dz24,fiz2);
				fix3            = vec_madd(fs34,dx34,fix3);
				fiy3            = vec_madd(fs34,dy34,fiy3);
				fiz3            = vec_madd(fs34,dz34,fiz3);
				fix4            = vec_madd(fs44,dx44,fix4);
				fiy4            = vec_madd(fs44,dy44,fiy4);
				fiz4            = vec_madd(fs44,dz44,fiz4);


				fjx1            = vec_nmsub(fs11,dx11,nul);
				fjy1            = vec_nmsub(fs11,dy11,nul);
				fjz1            = vec_nmsub(fs11,dz11,nul);
				fjx2            = vec_nmsub(fs22,dx22,nul);
				fjy2            = vec_nmsub(fs22,dy22,nul);
				fjz2            = vec_nmsub(fs22,dz22,nul);
				fjx3            = vec_nmsub(fs23,dx23,nul);
				fjy3            = vec_nmsub(fs23,dy23,nul);
				fjz3            = vec_nmsub(fs23,dz23,nul);
				fjx4            = vec_nmsub(fs24,dx24,nul);
				fjy4            = vec_nmsub(fs24,dy24,nul);
				fjz4            = vec_nmsub(fs24,dz24,nul);

				fjx2            = vec_nmsub(fs32,dx32,fjx2);
				fjy2            = vec_nmsub(fs32,dy32,fjy2);
				fjz2            = vec_nmsub(fs32,dz32,fjz2);
				fjx3            = vec_nmsub(fs33,dx33,fjx3);
				fjy3            = vec_nmsub(fs33,dy33,fjy3);
				fjz3            = vec_nmsub(fs33,dz33,fjz3);
				fjx4            = vec_nmsub(fs34,dx34,fjx4);
				fjy4            = vec_nmsub(fs34,dy34,fjy4);
				fjz4            = vec_nmsub(fs34,dz34,fjz4);

				fjx2            = vec_nmsub(fs42,dx42,fjx2);
				fjy2            = vec_nmsub(fs42,dy42,fjy2);
				fjz2            = vec_nmsub(fs42,dz42,fjz2);
				fjx3            = vec_nmsub(fs43,dx43,fjx3);
				fjy3            = vec_nmsub(fs43,dy43,fjy3);
				fjz3            = vec_nmsub(fs43,dz43,fjz3);
				fjx4            = vec_nmsub(fs44,dx44,fjx4);
				fjy4            = vec_nmsub(fs44,dy44,fjy4);
				fjz4            = vec_nmsub(fs44,dz44,fjz4);

				add_force_to_2_4atoms(faction+j3a,faction+j3b,
									  fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,
									  fjx3,fjy3,fjz3,fjx4,fjy4,fjz4);
			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				load_1_4atoms(pos+j3a,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3,
							  &jx4,&jy4,&jz4);

				qqMMt           = vec_sld(qqMM,nul,12);
				qqMHt           = vec_sld(qqMH,nul,12);
				qqHHt           = vec_sld(qqHH,nul,12);
				c6t             = vec_sld(c6,nul,12);
				c12t            = vec_sld(c12,nul,12);

				dx11            = vec_sub(ix1,jx1);
				dy11            = vec_sub(iy1,jy1);
				dz11            = vec_sub(iz1,jz1);
				dx22            = vec_sub(ix2,jx2);
				dy22            = vec_sub(iy2,jy2);
				dz22            = vec_sub(iz2,jz2);
				dx23            = vec_sub(ix2,jx3);
				dy23            = vec_sub(iy2,jy3);
				dz23            = vec_sub(iz2,jz3);
				dx24            = vec_sub(ix2,jx4);
				dy24            = vec_sub(iy2,jy4);
				dz24            = vec_sub(iz2,jz4);
				dx32            = vec_sub(ix3,jx2);
				dy32            = vec_sub(iy3,jy2);
				dz32            = vec_sub(iz3,jz2);
				dx33            = vec_sub(ix3,jx3);
				dy33            = vec_sub(iy3,jy3);
				dz33            = vec_sub(iz3,jz3);
				dx34            = vec_sub(ix3,jx4);
				dy34            = vec_sub(iy3,jy4);
				dz34            = vec_sub(iz3,jz4);
				dx42            = vec_sub(ix4,jx2);
				dy42            = vec_sub(iy4,jy2);
				dz42            = vec_sub(iz4,jz2);
				dx43            = vec_sub(ix4,jx3);
				dy43            = vec_sub(iy4,jy3);
				dz43            = vec_sub(iz4,jz3);
				dx44            = vec_sub(ix4,jx4);
				dy44            = vec_sub(iy4,jy4);
				dz44            = vec_sub(iz4,jz4);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq24           = vec_madd(dx24,dx24,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq34           = vec_madd(dx34,dx34,nul);
				rsq42           = vec_madd(dx42,dx42,nul);
				rsq43           = vec_madd(dx43,dx43,nul);
				rsq44           = vec_madd(dx44,dx44,nul);

				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq24           = vec_madd(dy24,dy24,rsq24);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq34           = vec_madd(dy34,dy34,rsq34);
				rsq42           = vec_madd(dy42,dy42,rsq42);
				rsq43           = vec_madd(dy43,dy43,rsq43);
				rsq44           = vec_madd(dy44,dy44,rsq44);
     
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq24           = vec_madd(dz24,dz24,rsq24);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);
				rsq34           = vec_madd(dz34,dz34,rsq34);
				rsq42           = vec_madd(dz42,dz42,rsq42);
				rsq43           = vec_madd(dz43,dz43,rsq43);
				rsq44           = vec_madd(dz44,dz44,rsq44);

				rinv11          = do_invsqrt(rsq11);
				do_9_invsqrt(rsq22,rsq23,rsq24,
							 rsq32,rsq33,rsq34,
							 rsq42,rsq43,rsq44,
							 &rinv22,&rinv23,&rinv24,
							 &rinv32,&rinv33,&rinv34,
							 &rinv42,&rinv43,&rinv44);

				zero_highest_3_elements_in_vector(&rinv11);
				zero_highest_3_elements_in_9_vectors(&rinv22,&rinv23,&rinv24,
													 &rinv32,&rinv33,&rinv34,
													 &rinv42,&rinv43,&rinv44);
      
				r               = vec_madd(rinv11,rsq11,nul);
				do_1_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&FFd,&VVr,&FFr);
				Vvdwtot         = vec_madd(c6t,VVd,Vvdwtot);
				fs11            = vec_nmsub(c6t,FFd,nul);
				Vvdwtot         = vec_madd(c12t,VVr,Vvdwtot);
				fs11            = vec_nmsub(c12t,FFr,fs11);
				fs11            = vec_madd(fs11,tsc,nul);
				
				krsq22          = vec_madd(vkrf,rsq22,nul);
				krsq23          = vec_madd(vkrf,rsq23,nul);
				krsq24          = vec_madd(vkrf,rsq24,nul);
				krsq32          = vec_madd(vkrf,rsq32,nul);
				krsq33          = vec_madd(vkrf,rsq33,nul);
				krsq34          = vec_madd(vkrf,rsq34,nul);
				krsq42          = vec_madd(vkrf,rsq42,nul);
				krsq43          = vec_madd(vkrf,rsq43,nul);
				krsq44          = vec_madd(vkrf,rsq44,nul);

				rinvsq22        = vec_madd(rinv22,rinv22,nul);
				rinvsq23        = vec_madd(rinv23,rinv23,nul);
				rinvsq24        = vec_madd(rinv24,rinv24,nul);
				rinvsq32        = vec_madd(rinv32,rinv32,nul);
				rinvsq33        = vec_madd(rinv33,rinv33,nul);
				rinvsq34        = vec_madd(rinv34,rinv34,nul);
				rinvsq42        = vec_madd(rinv42,rinv42,nul);
				rinvsq43        = vec_madd(rinv43,rinv43,nul);
				rinvsq44        = vec_madd(rinv44,rinv44,nul);

				vc22            = vec_add(rinv22,krsq22);
				vc23            = vec_add(rinv23,krsq23);
				vc24            = vec_add(rinv24,krsq24);
				vc32            = vec_add(rinv32,krsq32);
				vc33            = vec_add(rinv33,krsq33);
				vc34            = vec_add(rinv34,krsq34);
				vc42            = vec_add(rinv42,krsq42);
				vc43            = vec_add(rinv43,krsq43);
				vc44            = vec_add(rinv44,krsq44);

				vc22            = vec_sub(vc22,vcrf);
				vc23            = vec_sub(vc23,vcrf);
				vc24            = vec_sub(vc24,vcrf);
				vc32            = vec_sub(vc32,vcrf);
				vc33            = vec_sub(vc33,vcrf);
				vc34            = vec_sub(vc34,vcrf);
				vc42            = vec_sub(vc42,vcrf);
				vc43            = vec_sub(vc43,vcrf);
				vc44            = vec_sub(vc44,vcrf);

				fs22            = vec_nmsub(vec_two(),krsq22,rinv22);
				fs23            = vec_nmsub(vec_two(),krsq23,rinv23);
				fs24            = vec_nmsub(vec_two(),krsq24,rinv24);
				fs32            = vec_nmsub(vec_two(),krsq32,rinv32);
				fs33            = vec_nmsub(vec_two(),krsq33,rinv33);
				fs34            = vec_nmsub(vec_two(),krsq34,rinv34);
				fs42            = vec_nmsub(vec_two(),krsq42,rinv42);
				fs43            = vec_nmsub(vec_two(),krsq43,rinv43);
				fs44            = vec_nmsub(vec_two(),krsq44,rinv44);

				fs22            = vec_madd(fs22,qqHHt,nul);
				fs23            = vec_madd(fs23,qqHHt,nul);
				fs24            = vec_madd(fs24,qqMHt,nul);
				fs32            = vec_madd(fs32,qqHHt,nul);
				fs33            = vec_madd(fs33,qqHHt,nul);
				fs34            = vec_madd(fs34,qqMHt,nul);
				fs42            = vec_madd(fs42,qqMHt,nul);
				fs43            = vec_madd(fs43,qqMHt,nul);
				fs44            = vec_madd(fs44,qqMMt,nul);

				fs22            = vec_madd(fs22,rinvsq22,nul);
				fs23            = vec_madd(fs23,rinvsq23,nul);
				fs24            = vec_madd(fs24,rinvsq24,nul);
				fs32            = vec_madd(fs32,rinvsq32,nul);
				fs33            = vec_madd(fs33,rinvsq33,nul);
				fs34            = vec_madd(fs34,rinvsq34,nul);
				fs42            = vec_madd(fs42,rinvsq42,nul);
				fs43            = vec_madd(fs43,rinvsq43,nul);
				fs44            = vec_madd(fs44,rinvsq44,nul);
      
				vctot           = vec_madd(qqHHt,vc22,vctot);
				vctot           = vec_madd(qqHHt,vc23,vctot);
				vctot           = vec_madd(qqMHt,vc24,vctot);
				vctot           = vec_madd(qqHHt,vc32,vctot);
				vctot           = vec_madd(qqHHt,vc33,vctot);
				vctot           = vec_madd(qqMHt,vc34,vctot);
				vctot           = vec_madd(qqMHt,vc42,vctot);
				vctot           = vec_madd(qqMHt,vc43,vctot);
				vctot           = vec_madd(qqMMt,vc44,vctot);

				fs11            = vec_madd(fs11,rinv11,nul);

				fix1            = vec_madd(fs11,dx11,fix1);
				fiy1            = vec_madd(fs11,dy11,fiy1);
				fiz1            = vec_madd(fs11,dz11,fiz1);

				fix2            = vec_madd(fs22,dx22,fix2);
				fiy2            = vec_madd(fs22,dy22,fiy2);
				fiz2            = vec_madd(fs22,dz22,fiz2);
				fix3            = vec_madd(fs32,dx32,fix3);
				fiy3            = vec_madd(fs32,dy32,fiy3);
				fiz3            = vec_madd(fs32,dz32,fiz3);
				fix4            = vec_madd(fs42,dx42,fix4);
				fiy4            = vec_madd(fs42,dy42,fiy4);
				fiz4            = vec_madd(fs42,dz42,fiz4);

				fix2            = vec_madd(fs23,dx23,fix2);
				fiy2            = vec_madd(fs23,dy23,fiy2);
				fiz2            = vec_madd(fs23,dz23,fiz2);
				fix3            = vec_madd(fs33,dx33,fix3);
				fiy3            = vec_madd(fs33,dy33,fiy3);
				fiz3            = vec_madd(fs33,dz33,fiz3);
				fix4            = vec_madd(fs43,dx43,fix4);
				fiy4            = vec_madd(fs43,dy43,fiy4);
				fiz4            = vec_madd(fs43,dz43,fiz4);

				fix2            = vec_madd(fs24,dx24,fix2);
				fiy2            = vec_madd(fs24,dy24,fiy2);
				fiz2            = vec_madd(fs24,dz24,fiz2);
				fix3            = vec_madd(fs34,dx34,fix3);
				fiy3            = vec_madd(fs34,dy34,fiy3);
				fiz3            = vec_madd(fs34,dz34,fiz3);
				fix4            = vec_madd(fs44,dx44,fix4);
				fiy4            = vec_madd(fs44,dy44,fiy4);
				fiz4            = vec_madd(fs44,dz44,fiz4);


				fjx1            = vec_nmsub(fs11,dx11,nul);
				fjy1            = vec_nmsub(fs11,dy11,nul);
				fjz1            = vec_nmsub(fs11,dz11,nul);
				fjx2            = vec_nmsub(fs22,dx22,nul);
				fjy2            = vec_nmsub(fs22,dy22,nul);
				fjz2            = vec_nmsub(fs22,dz22,nul);
				fjx3            = vec_nmsub(fs23,dx23,nul);
				fjy3            = vec_nmsub(fs23,dy23,nul);
				fjz3            = vec_nmsub(fs23,dz23,nul);
				fjx4            = vec_nmsub(fs24,dx24,nul);
				fjy4            = vec_nmsub(fs24,dy24,nul);
				fjz4            = vec_nmsub(fs24,dz24,nul);

				fjx2            = vec_nmsub(fs32,dx32,fjx2);
				fjy2            = vec_nmsub(fs32,dy32,fjy2);
				fjz2            = vec_nmsub(fs32,dz32,fjz2);
				fjx3            = vec_nmsub(fs33,dx33,fjx3);
				fjy3            = vec_nmsub(fs33,dy33,fjy3);
				fjz3            = vec_nmsub(fs33,dz33,fjz3);
				fjx4            = vec_nmsub(fs34,dx34,fjx4);
				fjy4            = vec_nmsub(fs34,dy34,fjy4);
				fjz4            = vec_nmsub(fs34,dz34,fjz4);

				fjx2            = vec_nmsub(fs42,dx42,fjx2);
				fjy2            = vec_nmsub(fs42,dy42,fjy2);
				fjz2            = vec_nmsub(fs42,dz42,fjz2);
				fjx3            = vec_nmsub(fs43,dx43,fjx3);
				fjy3            = vec_nmsub(fs43,dy43,fjy3);
				fjz3            = vec_nmsub(fs43,dz43,fjz3);
				fjx4            = vec_nmsub(fs44,dx44,fjx4);
				fjy4            = vec_nmsub(fs44,dy44,fjy4);
				fjz4            = vec_nmsub(fs44,dz44,fjz4);

				add_force_to_1_4atoms(faction+j3a,
									  fjx1,fjy1,fjz1,fjx2,fjy2,fjz2,
									  fjx3,fjy3,fjz3,fjx4,fjy4,fjz4);
			}
			/* update outer data */
			update_i_4atoms_forces(faction+ii3,fshift+is3,
								   fix1,fiy1,fiz1,fix2,fiy2,fiz2,
								   fix3,fiy3,fiz3,fix4,fiy4,fiz4);

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
nb_kernel234nf_ppc_altivec(int *             p_nri,
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
	vector float ix1,iy1,iz1,ix2,iy2,iz2,ix3,iy3,iz3,ix4,iy4,iz4;
	vector float jx1,jy1,jz1,jx2,jy2,jz2,jx3,jy3,jz3,jx4,jy4,jz4;

	vector float dx11,dy11,dz11;
	vector float dx22,dy22,dz22,dx23,dy23,dz23,dx24,dy24,dz24;
	vector float dx32,dy32,dz32,dx33,dy33,dz33,dx34,dy34,dz34;
	vector float dx42,dy42,dz42,dx43,dy43,dz43,dx44,dy44,dz44;

	vector float rsq11,rsq22,rsq23,rsq24,rsq32,rsq33,rsq34,rsq42,rsq43,rsq44;
	vector float rinv11,rinv22,rinv23,rinv24,rinv32,rinv33,rinv34;
	vector float rinv42,rinv43,rinv44;

	vector float vfacel,nul;
	vector float vctot,Vvdwtot,qqMM,qqMH,qqHH,qM,qH,qqMMt,qqMHt,qqHHt; 
	vector float krsq22,krsq23,krsq24,krsq32,krsq33,krsq34;
	vector float krsq42,krsq43,krsq44;
	vector float vkrf,vcrf;
	vector float c6,c12,c6t,c12t;
	vector float VVd,VVr,tsc,r;

	int n,k,ii,is3,ii3,nj0,nj1,tp,tj;
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
	tsc=load_float_and_splat(p_tabscale);
	vfacel=load_float_and_splat(p_facel);
	vkrf=load_float_and_splat(p_krf);
	vcrf=load_float_and_splat(p_crf);
	ii        = iinr[0];
	qH        = load_float_and_splat(charge+iinr[0]+1);
	qM        = load_float_and_splat(charge+iinr[0]+3);
	qqMM      = vec_madd(qM,qM,nul);
	qqMH      = vec_madd(qM,qH,nul);
	qqHH      = vec_madd(qH,qH,nul);
	qqMM      = vec_madd(qqMM,vfacel,nul);
	qqMH      = vec_madd(qqMH,vfacel,nul);
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
			load_1_4atoms_shift_and_splat(pos+ii3,shiftvec+is3,&ix1,&iy1,&iz1,
										  &ix2,&iy2,&iz2,&ix3,&iy3,&iz3,
										  &ix4,&iy4,&iz4);
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
				load_4_4atoms(pos+j3a,pos+j3b,pos+j3c,pos+j3d,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3,
							  &jx4,&jy4,&jz4);

				dx11            = vec_sub(ix1,jx1);
				dy11            = vec_sub(iy1,jy1);
				dz11            = vec_sub(iz1,jz1);
				dx22            = vec_sub(ix2,jx2);
				dy22            = vec_sub(iy2,jy2);
				dz22            = vec_sub(iz2,jz2);
				dx23            = vec_sub(ix2,jx3);
				dy23            = vec_sub(iy2,jy3);
				dz23            = vec_sub(iz2,jz3);
				dx24            = vec_sub(ix2,jx4);
				dy24            = vec_sub(iy2,jy4);
				dz24            = vec_sub(iz2,jz4);
				dx32            = vec_sub(ix3,jx2);
				dy32            = vec_sub(iy3,jy2);
				dz32            = vec_sub(iz3,jz2);
				dx33            = vec_sub(ix3,jx3);
				dy33            = vec_sub(iy3,jy3);
				dz33            = vec_sub(iz3,jz3);
				dx34            = vec_sub(ix3,jx4);
				dy34            = vec_sub(iy3,jy4);
				dz34            = vec_sub(iz3,jz4);
				dx42            = vec_sub(ix4,jx2);
				dy42            = vec_sub(iy4,jy2);
				dz42            = vec_sub(iz4,jz2);
				dx43            = vec_sub(ix4,jx3);
				dy43            = vec_sub(iy4,jy3);
				dz43            = vec_sub(iz4,jz3);
				dx44            = vec_sub(ix4,jx4);
				dy44            = vec_sub(iy4,jy4);
				dz44            = vec_sub(iz4,jz4);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq24           = vec_madd(dx24,dx24,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq34           = vec_madd(dx34,dx34,nul);
				rsq42           = vec_madd(dx42,dx42,nul);
				rsq43           = vec_madd(dx43,dx43,nul);
				rsq44           = vec_madd(dx44,dx44,nul);

				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq24           = vec_madd(dy24,dy24,rsq24);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq34           = vec_madd(dy34,dy34,rsq34);
				rsq42           = vec_madd(dy42,dy42,rsq42);
				rsq43           = vec_madd(dy43,dy43,rsq43);
				rsq44           = vec_madd(dy44,dy44,rsq44);
     
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq24           = vec_madd(dz24,dz24,rsq24);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);
				rsq34           = vec_madd(dz34,dz34,rsq34);
				rsq42           = vec_madd(dz42,dz42,rsq42);
				rsq43           = vec_madd(dz43,dz43,rsq43);
				rsq44           = vec_madd(dz44,dz44,rsq44);

				rinv11          = do_invsqrt(rsq11);
				do_9_invsqrt(rsq22,rsq23,rsq24,
							 rsq32,rsq33,rsq34,
							 rsq42,rsq43,rsq44,
							 &rinv22,&rinv23,&rinv24,
							 &rinv32,&rinv33,&rinv34,
							 &rinv42,&rinv43,&rinv44);
      
				r               = vec_madd(rinv11,rsq11,nul);
				do_vonly_4_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&VVr);
				Vvdwtot         = vec_madd(c6,VVd,Vvdwtot);
				Vvdwtot         = vec_madd(c12,VVr,Vvdwtot);

				krsq22          = vec_madd(vkrf,rsq22,nul);
				krsq23          = vec_madd(vkrf,rsq23,nul);
				krsq24          = vec_madd(vkrf,rsq24,nul);
				krsq32          = vec_madd(vkrf,rsq32,nul);
				krsq33          = vec_madd(vkrf,rsq33,nul);
				krsq34          = vec_madd(vkrf,rsq34,nul);
				krsq42          = vec_madd(vkrf,rsq42,nul);
				krsq43          = vec_madd(vkrf,rsq43,nul);
				krsq44          = vec_madd(vkrf,rsq44,nul);

				rinv22          = vec_add(rinv22,krsq22);
				rinv23          = vec_add(rinv23,krsq23);
				rinv24          = vec_add(rinv24,krsq24);
				rinv32          = vec_add(rinv32,krsq32);
				rinv33          = vec_add(rinv33,krsq33);
				rinv34          = vec_add(rinv34,krsq34);
				rinv42          = vec_add(rinv42,krsq42);
				rinv43          = vec_add(rinv43,krsq43);
				rinv44          = vec_add(rinv44,krsq44);

				rinv22          = vec_sub(rinv22,vcrf);
				rinv23          = vec_sub(rinv23,vcrf);
				rinv24          = vec_sub(rinv24,vcrf);
				rinv32          = vec_sub(rinv32,vcrf);
				rinv33          = vec_sub(rinv33,vcrf);
				rinv34          = vec_sub(rinv34,vcrf);
				rinv42          = vec_sub(rinv42,vcrf);
				rinv43          = vec_sub(rinv43,vcrf);
				rinv44          = vec_sub(rinv44,vcrf);
      
				vctot           = vec_madd(qqHH,rinv22,vctot);
				vctot           = vec_madd(qqHH,rinv23,vctot);
				vctot           = vec_madd(qqMH,rinv24,vctot);
				vctot           = vec_madd(qqHH,rinv32,vctot);
				vctot           = vec_madd(qqHH,rinv33,vctot);
				vctot           = vec_madd(qqMH,rinv34,vctot);
				vctot           = vec_madd(qqMH,rinv42,vctot);
				vctot           = vec_madd(qqMH,rinv43,vctot);
				vctot           = vec_madd(qqMM,rinv44,vctot);
			} 
			if(k<(nj1-2)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				load_3_4atoms(pos+j3a,pos+j3b,pos+j3c,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3,
							  &jx4,&jy4,&jz4);

				qqMMt           = vec_sld(qqMM,nul,4);
				qqMHt           = vec_sld(qqMH,nul,4);
				qqHHt           = vec_sld(qqHH,nul,4);
				c6t             = vec_sld(c6,nul,4);
				c12t            = vec_sld(c12,nul,4);

				dx11            = vec_sub(ix1,jx1);
				dy11            = vec_sub(iy1,jy1);
				dz11            = vec_sub(iz1,jz1);
				dx22            = vec_sub(ix2,jx2);
				dy22            = vec_sub(iy2,jy2);
				dz22            = vec_sub(iz2,jz2);
				dx23            = vec_sub(ix2,jx3);
				dy23            = vec_sub(iy2,jy3);
				dz23            = vec_sub(iz2,jz3);
				dx24            = vec_sub(ix2,jx4);
				dy24            = vec_sub(iy2,jy4);
				dz24            = vec_sub(iz2,jz4);
				dx32            = vec_sub(ix3,jx2);
				dy32            = vec_sub(iy3,jy2);
				dz32            = vec_sub(iz3,jz2);
				dx33            = vec_sub(ix3,jx3);
				dy33            = vec_sub(iy3,jy3);
				dz33            = vec_sub(iz3,jz3);
				dx34            = vec_sub(ix3,jx4);
				dy34            = vec_sub(iy3,jy4);
				dz34            = vec_sub(iz3,jz4);
				dx42            = vec_sub(ix4,jx2);
				dy42            = vec_sub(iy4,jy2);
				dz42            = vec_sub(iz4,jz2);
				dx43            = vec_sub(ix4,jx3);
				dy43            = vec_sub(iy4,jy3);
				dz43            = vec_sub(iz4,jz3);
				dx44            = vec_sub(ix4,jx4);
				dy44            = vec_sub(iy4,jy4);
				dz44            = vec_sub(iz4,jz4);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq24           = vec_madd(dx24,dx24,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq34           = vec_madd(dx34,dx34,nul);
				rsq42           = vec_madd(dx42,dx42,nul);
				rsq43           = vec_madd(dx43,dx43,nul);
				rsq44           = vec_madd(dx44,dx44,nul);

				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq24           = vec_madd(dy24,dy24,rsq24);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq34           = vec_madd(dy34,dy34,rsq34);
				rsq42           = vec_madd(dy42,dy42,rsq42);
				rsq43           = vec_madd(dy43,dy43,rsq43);
				rsq44           = vec_madd(dy44,dy44,rsq44);
     
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq24           = vec_madd(dz24,dz24,rsq24);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);
				rsq34           = vec_madd(dz34,dz34,rsq34);
				rsq42           = vec_madd(dz42,dz42,rsq42);
				rsq43           = vec_madd(dz43,dz43,rsq43);
				rsq44           = vec_madd(dz44,dz44,rsq44);

				rinv11          = do_invsqrt(rsq11);
				do_9_invsqrt(rsq22,rsq23,rsq24,
							 rsq32,rsq33,rsq34,
							 rsq42,rsq43,rsq44,
							 &rinv22,&rinv23,&rinv24,
							 &rinv32,&rinv33,&rinv34,
							 &rinv42,&rinv43,&rinv44);

				zero_highest_element_in_vector(&rinv11);
				zero_highest_element_in_9_vectors(&rinv22,&rinv23,&rinv24,
												  &rinv32,&rinv33,&rinv34,
												  &rinv42,&rinv43,&rinv44);
      
      
				r               = vec_madd(rinv11,rsq11,nul);
				do_vonly_3_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&VVr);
				Vvdwtot         = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot         = vec_madd(c12t,VVr,Vvdwtot);

				krsq22          = vec_madd(vkrf,rsq22,nul);
				krsq23          = vec_madd(vkrf,rsq23,nul);
				krsq24          = vec_madd(vkrf,rsq24,nul);
				krsq32          = vec_madd(vkrf,rsq32,nul);
				krsq33          = vec_madd(vkrf,rsq33,nul);
				krsq34          = vec_madd(vkrf,rsq34,nul);
				krsq42          = vec_madd(vkrf,rsq42,nul);
				krsq43          = vec_madd(vkrf,rsq43,nul);
				krsq44          = vec_madd(vkrf,rsq44,nul);


				rinv22          = vec_add(rinv22,krsq22);
				rinv23          = vec_add(rinv23,krsq23);
				rinv24          = vec_add(rinv24,krsq24);
				rinv32          = vec_add(rinv32,krsq32);
				rinv33          = vec_add(rinv33,krsq33);
				rinv34          = vec_add(rinv34,krsq34);
				rinv42          = vec_add(rinv42,krsq42);
				rinv43          = vec_add(rinv43,krsq43);
				rinv44          = vec_add(rinv44,krsq44);

				rinv22          = vec_sub(rinv22,vcrf);
				rinv23          = vec_sub(rinv23,vcrf);
				rinv24          = vec_sub(rinv24,vcrf);
				rinv32          = vec_sub(rinv32,vcrf);
				rinv33          = vec_sub(rinv33,vcrf);
				rinv34          = vec_sub(rinv34,vcrf);
				rinv42          = vec_sub(rinv42,vcrf);
				rinv43          = vec_sub(rinv43,vcrf);
				rinv44          = vec_sub(rinv44,vcrf);
      
				vctot           = vec_madd(qqHHt,rinv22,vctot);
				vctot           = vec_madd(qqHHt,rinv23,vctot);
				vctot           = vec_madd(qqMHt,rinv24,vctot);
				vctot           = vec_madd(qqHHt,rinv32,vctot);
				vctot           = vec_madd(qqHHt,rinv33,vctot);
				vctot           = vec_madd(qqMHt,rinv34,vctot);
				vctot           = vec_madd(qqMHt,rinv42,vctot);
				vctot           = vec_madd(qqMHt,rinv43,vctot);
				vctot           = vec_madd(qqMMt,rinv44,vctot);
			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				load_2_4atoms(pos+j3a,pos+j3b,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3,
							  &jx4,&jy4,&jz4);

				qqMMt           = vec_sld(qqMM,nul,8);
				qqMHt           = vec_sld(qqMH,nul,8);
				qqHHt           = vec_sld(qqHH,nul,8);
				c6t             = vec_sld(c6,nul,8);
				c12t            = vec_sld(c12,nul,8);

				dx11            = vec_sub(ix1,jx1);
				dy11            = vec_sub(iy1,jy1);
				dz11            = vec_sub(iz1,jz1);
				dx22            = vec_sub(ix2,jx2);
				dy22            = vec_sub(iy2,jy2);
				dz22            = vec_sub(iz2,jz2);
				dx23            = vec_sub(ix2,jx3);
				dy23            = vec_sub(iy2,jy3);
				dz23            = vec_sub(iz2,jz3);
				dx24            = vec_sub(ix2,jx4);
				dy24            = vec_sub(iy2,jy4);
				dz24            = vec_sub(iz2,jz4);
				dx32            = vec_sub(ix3,jx2);
				dy32            = vec_sub(iy3,jy2);
				dz32            = vec_sub(iz3,jz2);
				dx33            = vec_sub(ix3,jx3);
				dy33            = vec_sub(iy3,jy3);
				dz33            = vec_sub(iz3,jz3);
				dx34            = vec_sub(ix3,jx4);
				dy34            = vec_sub(iy3,jy4);
				dz34            = vec_sub(iz3,jz4);
				dx42            = vec_sub(ix4,jx2);
				dy42            = vec_sub(iy4,jy2);
				dz42            = vec_sub(iz4,jz2);
				dx43            = vec_sub(ix4,jx3);
				dy43            = vec_sub(iy4,jy3);
				dz43            = vec_sub(iz4,jz3);
				dx44            = vec_sub(ix4,jx4);
				dy44            = vec_sub(iy4,jy4);
				dz44            = vec_sub(iz4,jz4);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq24           = vec_madd(dx24,dx24,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq34           = vec_madd(dx34,dx34,nul);
				rsq42           = vec_madd(dx42,dx42,nul);
				rsq43           = vec_madd(dx43,dx43,nul);
				rsq44           = vec_madd(dx44,dx44,nul);

				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq24           = vec_madd(dy24,dy24,rsq24);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq34           = vec_madd(dy34,dy34,rsq34);
				rsq42           = vec_madd(dy42,dy42,rsq42);
				rsq43           = vec_madd(dy43,dy43,rsq43);
				rsq44           = vec_madd(dy44,dy44,rsq44);
     
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq24           = vec_madd(dz24,dz24,rsq24);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);
				rsq34           = vec_madd(dz34,dz34,rsq34);
				rsq42           = vec_madd(dz42,dz42,rsq42);
				rsq43           = vec_madd(dz43,dz43,rsq43);
				rsq44           = vec_madd(dz44,dz44,rsq44);

				rinv11          = do_invsqrt(rsq11);
				do_9_invsqrt(rsq22,rsq23,rsq24,
							 rsq32,rsq33,rsq34,
							 rsq42,rsq43,rsq44,
							 &rinv22,&rinv23,&rinv24,
							 &rinv32,&rinv33,&rinv34,
							 &rinv42,&rinv43,&rinv44);

				zero_highest_2_elements_in_vector(&rinv11);
				zero_highest_2_elements_in_9_vectors(&rinv22,&rinv23,&rinv24,
													 &rinv32,&rinv33,&rinv34,
													 &rinv42,&rinv43,&rinv44);

				r               = vec_madd(rinv11,rsq11,nul);
				do_vonly_2_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&VVr);
				Vvdwtot         = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot         = vec_madd(c12t,VVr,Vvdwtot);

				krsq22          = vec_madd(vkrf,rsq22,nul);
				krsq23          = vec_madd(vkrf,rsq23,nul);
				krsq24          = vec_madd(vkrf,rsq24,nul);
				krsq32          = vec_madd(vkrf,rsq32,nul);
				krsq33          = vec_madd(vkrf,rsq33,nul);
				krsq34          = vec_madd(vkrf,rsq34,nul);
				krsq42          = vec_madd(vkrf,rsq42,nul);
				krsq43          = vec_madd(vkrf,rsq43,nul);
				krsq44          = vec_madd(vkrf,rsq44,nul);


				rinv22          = vec_add(rinv22,krsq22);
				rinv23          = vec_add(rinv23,krsq23);
				rinv24          = vec_add(rinv24,krsq24);
				rinv32          = vec_add(rinv32,krsq32);
				rinv33          = vec_add(rinv33,krsq33);
				rinv34          = vec_add(rinv34,krsq34);
				rinv42          = vec_add(rinv42,krsq42);
				rinv43          = vec_add(rinv43,krsq43);
				rinv44          = vec_add(rinv44,krsq44);

				rinv22          = vec_sub(rinv22,vcrf);
				rinv23          = vec_sub(rinv23,vcrf);
				rinv24          = vec_sub(rinv24,vcrf);
				rinv32          = vec_sub(rinv32,vcrf);
				rinv33          = vec_sub(rinv33,vcrf);
				rinv34          = vec_sub(rinv34,vcrf);
				rinv42          = vec_sub(rinv42,vcrf);
				rinv43          = vec_sub(rinv43,vcrf);
				rinv44          = vec_sub(rinv44,vcrf);
      
				vctot           = vec_madd(qqHHt,rinv22,vctot);
				vctot           = vec_madd(qqHHt,rinv23,vctot);
				vctot           = vec_madd(qqMHt,rinv24,vctot);
				vctot           = vec_madd(qqHHt,rinv32,vctot);
				vctot           = vec_madd(qqHHt,rinv33,vctot);
				vctot           = vec_madd(qqMHt,rinv34,vctot);
				vctot           = vec_madd(qqMHt,rinv42,vctot);
				vctot           = vec_madd(qqMHt,rinv43,vctot);
				vctot           = vec_madd(qqMMt,rinv44,vctot);
			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;
				load_1_4atoms(pos+j3a,
							  &jx1,&jy1,&jz1,&jx2,&jy2,&jz2,&jx3,&jy3,&jz3,
							  &jx4,&jy4,&jz4);

				qqMMt           = vec_sld(qqMM,nul,12);
				qqMHt           = vec_sld(qqMH,nul,12);
				qqHHt           = vec_sld(qqHH,nul,12);
				c6t             = vec_sld(c6,nul,12);
				c12t            = vec_sld(c12,nul,12);

				dx11            = vec_sub(ix1,jx1);
				dy11            = vec_sub(iy1,jy1);
				dz11            = vec_sub(iz1,jz1);
				dx22            = vec_sub(ix2,jx2);
				dy22            = vec_sub(iy2,jy2);
				dz22            = vec_sub(iz2,jz2);
				dx23            = vec_sub(ix2,jx3);
				dy23            = vec_sub(iy2,jy3);
				dz23            = vec_sub(iz2,jz3);
				dx24            = vec_sub(ix2,jx4);
				dy24            = vec_sub(iy2,jy4);
				dz24            = vec_sub(iz2,jz4);
				dx32            = vec_sub(ix3,jx2);
				dy32            = vec_sub(iy3,jy2);
				dz32            = vec_sub(iz3,jz2);
				dx33            = vec_sub(ix3,jx3);
				dy33            = vec_sub(iy3,jy3);
				dz33            = vec_sub(iz3,jz3);
				dx34            = vec_sub(ix3,jx4);
				dy34            = vec_sub(iy3,jy4);
				dz34            = vec_sub(iz3,jz4);
				dx42            = vec_sub(ix4,jx2);
				dy42            = vec_sub(iy4,jy2);
				dz42            = vec_sub(iz4,jz2);
				dx43            = vec_sub(ix4,jx3);
				dy43            = vec_sub(iy4,jy3);
				dz43            = vec_sub(iz4,jz3);
				dx44            = vec_sub(ix4,jx4);
				dy44            = vec_sub(iy4,jy4);
				dz44            = vec_sub(iz4,jz4);

				rsq11           = vec_madd(dx11,dx11,nul);
				rsq22           = vec_madd(dx22,dx22,nul);
				rsq23           = vec_madd(dx23,dx23,nul);
				rsq24           = vec_madd(dx24,dx24,nul);
				rsq32           = vec_madd(dx32,dx32,nul);
				rsq33           = vec_madd(dx33,dx33,nul);
				rsq34           = vec_madd(dx34,dx34,nul);
				rsq42           = vec_madd(dx42,dx42,nul);
				rsq43           = vec_madd(dx43,dx43,nul);
				rsq44           = vec_madd(dx44,dx44,nul);

				rsq11           = vec_madd(dy11,dy11,rsq11);
				rsq22           = vec_madd(dy22,dy22,rsq22);
				rsq23           = vec_madd(dy23,dy23,rsq23);
				rsq24           = vec_madd(dy24,dy24,rsq24);
				rsq32           = vec_madd(dy32,dy32,rsq32);
				rsq33           = vec_madd(dy33,dy33,rsq33);
				rsq34           = vec_madd(dy34,dy34,rsq34);
				rsq42           = vec_madd(dy42,dy42,rsq42);
				rsq43           = vec_madd(dy43,dy43,rsq43);
				rsq44           = vec_madd(dy44,dy44,rsq44);
     
				rsq11           = vec_madd(dz11,dz11,rsq11);
				rsq22           = vec_madd(dz22,dz22,rsq22);
				rsq23           = vec_madd(dz23,dz23,rsq23);
				rsq24           = vec_madd(dz24,dz24,rsq24);
				rsq32           = vec_madd(dz32,dz32,rsq32);
				rsq33           = vec_madd(dz33,dz33,rsq33);
				rsq34           = vec_madd(dz34,dz34,rsq34);
				rsq42           = vec_madd(dz42,dz42,rsq42);
				rsq43           = vec_madd(dz43,dz43,rsq43);
				rsq44           = vec_madd(dz44,dz44,rsq44);

				rinv11          = do_invsqrt(rsq11);
				do_9_invsqrt(rsq22,rsq23,rsq24,
							 rsq32,rsq33,rsq34,
							 rsq42,rsq43,rsq44,
							 &rinv22,&rinv23,&rinv24,
							 &rinv32,&rinv33,&rinv34,
							 &rinv42,&rinv43,&rinv44);

				zero_highest_3_elements_in_vector(&rinv11);
				zero_highest_3_elements_in_9_vectors(&rinv22,&rinv23,&rinv24,
													 &rinv32,&rinv33,&rinv34,
													 &rinv42,&rinv43,&rinv44);

      
				r               = vec_madd(rinv11,rsq11,nul);
				do_vonly_1_ljtable_lj(VFtab,vec_madd(r,tsc,nul),&VVd,&VVr);
				Vvdwtot         = vec_madd(c6t,VVd,Vvdwtot);
				Vvdwtot         = vec_madd(c12t,VVr,Vvdwtot);

				krsq22          = vec_madd(vkrf,rsq22,nul);
				krsq23          = vec_madd(vkrf,rsq23,nul);
				krsq24          = vec_madd(vkrf,rsq24,nul);
				krsq32          = vec_madd(vkrf,rsq32,nul);
				krsq33          = vec_madd(vkrf,rsq33,nul);
				krsq34          = vec_madd(vkrf,rsq34,nul);
				krsq42          = vec_madd(vkrf,rsq42,nul);
				krsq43          = vec_madd(vkrf,rsq43,nul);
				krsq44          = vec_madd(vkrf,rsq44,nul);


				rinv22          = vec_add(rinv22,krsq22);
				rinv23          = vec_add(rinv23,krsq23);
				rinv24          = vec_add(rinv24,krsq24);
				rinv32          = vec_add(rinv32,krsq32);
				rinv33          = vec_add(rinv33,krsq33);
				rinv34          = vec_add(rinv34,krsq34);
				rinv42          = vec_add(rinv42,krsq42);
				rinv43          = vec_add(rinv43,krsq43);
				rinv44          = vec_add(rinv44,krsq44);

				rinv22          = vec_sub(rinv22,vcrf);
				rinv23          = vec_sub(rinv23,vcrf);
				rinv24          = vec_sub(rinv24,vcrf);
				rinv32          = vec_sub(rinv32,vcrf);
				rinv33          = vec_sub(rinv33,vcrf);
				rinv34          = vec_sub(rinv34,vcrf);
				rinv42          = vec_sub(rinv42,vcrf);
				rinv43          = vec_sub(rinv43,vcrf);
				rinv44          = vec_sub(rinv44,vcrf);
      
				vctot           = vec_madd(qqHHt,rinv22,vctot);
				vctot           = vec_madd(qqHHt,rinv23,vctot);
				vctot           = vec_madd(qqMHt,rinv24,vctot);
				vctot           = vec_madd(qqHHt,rinv32,vctot);
				vctot           = vec_madd(qqHHt,rinv33,vctot);
				vctot           = vec_madd(qqMHt,rinv34,vctot);
				vctot           = vec_madd(qqMHt,rinv42,vctot);
				vctot           = vec_madd(qqMHt,rinv43,vctot);
				vctot           = vec_madd(qqMMt,rinv44,vctot);
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

