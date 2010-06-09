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
#include "nb_kernel112_ppc_altivec.h"

/* NB: This is one of the most common nonbonded functions called, so I tried
 * to optimize it by doing all register saving/restoring manually. It only gave
 * 5-10% better performance, so I have not implemented it in the other loops.
 * (it makes the code horrible to read). However, since it is slightly faster
 * there is no reason not to keep the optimized code...
 */


void 
nb_kernel112_ppc_altivec  (int *             p_nri,
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
	register vector float v0;
	register vector float v1;
	register vector float v2;
	register vector float v3;
	register vector float v4;
	register vector float v5;
	register vector float v6;
	register vector float v7;
	register vector float v8;
	register vector float v9;
	register vector float v10;
	register vector float v11;
	register vector float v12;
	register vector float v13;
	register vector float v14;
	register vector float v15;
	register vector float v16;
	register vector float v17;
	register vector float v18;
	register vector float v19;
	register vector float v20;
	register vector float v21;
	register vector float v22;
	register vector float v23;
	register vector float v24;
	register vector float v25;
	register vector float v26;
	register vector float v27;
	register vector float v28;
	register vector float v29;
	register vector float v30;
	register vector float v31;

	union vfloat {
		float f[4];
		vector float v;
	} stackdata[52];
  
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

	/* set non java mode */
	v10       = (vector float)vec_mfvscr();
	v11       = (vector float)vec_sl(vec_splat_u32(1),vec_splat_u32(8));

	v12       = (vector float)vec_sl((vector unsigned int)v11,
									 vec_splat_u32(8));

	v10       = (vector float)vec_or((vector unsigned short)v10,
									 (vector unsigned short)v12);

	vec_mtvscr((vector unsigned short)v10);

	v0        = (vector float)vec_splat_u32(0);
	v0        = vec_ctf((vector unsigned int)v0,0);     /* load 0 to v0 */
	v1        = vec_lde(0,p_facel); /* load facel float to a vector */
	v2        = (vector float) vec_lvsl(0,p_facel); 
	v1        = vec_perm(v1,v1,(vector unsigned char) v2); /* move to elem 0 */
	v1        = vec_splat(v1,0); /* splat it to all elem */
  	
	ii        = iinr[0];
  
	v3        = vec_lde(0,charge+ii); /* load qO float to a vector */
	v4        = (vector float) vec_lvsl(0,charge+ii); 
	v3        = vec_perm(v3,v3,(vector unsigned char) v4); /* move to elem 0 */
	v3        = vec_splat(v3,0); /* splat it to all elem */

	v5        = vec_lde(0,charge+ii+1); /* load qH float to a vector */
	v6        = (vector float) vec_lvsl(0,charge+ii+1); 
	v5        = vec_perm(v5,v5,(vector unsigned char) v6); /* move to elem 0 */
	v5        = vec_splat(v5,0); /* splat it to all elem */

	v4        = vec_madd(v3,v5,v0); /* qqOH */
	v3        = vec_madd(v3,v3,v0); /* qqOO */
	v5        = vec_madd(v5,v5,v0); /* qqHH */
	v4        = vec_madd(v4,v1,v0); /* qqOH * facel */
	v3        = vec_madd(v3,v1,v0); /* qqOO * facel */
	v5        = vec_madd(v5,v1,v0); /* qqHH * facel */

	n         = 2*type[ii];
	n        = (ntype+1)*n;
  
	v1        = vec_ld( 0,vdwparam+n);  /* c6a c12a - the vdwparam array is at least
									 * 8-byte aligned and n is even here.
									 */
	v2        = (vector float) vec_lvsl(0,vdwparam+n);
	v1        = vec_perm(v1,v1,(vector unsigned char)v2); /* c6 c12 in 0,1 */
	v2        = vec_splat(v1,1);  /* c12 in all elements */
	v1        = vec_splat(v1,0);  /* c6 in all elements */

	/* store things to stack before starting outer loop */
	vec_st(v3,  0, (float *) stackdata); /* qqOO*facel is in stack pos 0 */
	vec_st(v4, 16, (float *) stackdata); /* qqOH*facel is in stack pos 1 */
	vec_st(v5, 32, (float *) stackdata); /* qqHH*facel is in stack pos 2 */
	vec_st(v1, 48, (float *) stackdata); /* c6 is in stack pos 3  */
	vec_st(v2, 64, (float *) stackdata); /* c12 is in stack pos 4 */
  
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
			/* load shift */
			/* load three consecutive shiftvector floats. 
             * We never access the fourth element,
			 * so this is safe even at the end of an array. 
			 */

			v4         = (vector float)vec_lvsl(0, shiftvec+is3);
			v1         = vec_lde(0, shiftvec+is3);
			v2         = vec_lde(4, shiftvec+is3);
			v3         = vec_lde(8, shiftvec+is3);

			/* Load shX,shY,shZ to elem 0 of v1,v2,v3 */
			v1         = vec_perm(v1,v1,(vector unsigned char)v4); 
			v2         = vec_perm(v2,v2,(vector unsigned char)v4); 
			v3         = vec_perm(v3,v3,(vector unsigned char)v4);
			v2         = vec_sld(v2,v2,4);
			v3         = vec_sld(v3,v3,8);
			v1         = vec_mergeh(v1,v3);
			v1         = vec_mergeh(v1,v2);  /* [ shX shY shZ - ] */
			/* load i coordinates */
			v2         = (vector float)vec_lvsl(0, pos+ii3);
			/* load 3atoms coords into three vectors. 
             * We do not yet know how it is aligned. 
             */
			v3         = vec_ld(0, pos+ii3); 
			v4         = vec_ld(16, pos+ii3);
			v5         = vec_ld(32, pos+ii3);
			v6         = vec_sld(v1,v1,12); /*  - shX shY shZ   */
			v7         = vec_sld(v6,v1,4);  /*  shX shY shZ shX */
			v8         = vec_sld(v6,v1,8);  /*  shY shZ shX shY */
			v9         = vec_sld(v6,v1,12); /*  shZ shX shY shZ */
			/* v3 = Ox  Oy  Oz H1x */
			v3         = vec_perm(v3,v4,(vector unsigned char)v2);
			/* v4 = H1y H1z H2x H2y */
			v4         = vec_perm(v4,v5,(vector unsigned char)v2); 
			/* v5 = H2z   -   -   - */
			v5         = vec_perm(v5,v5,(vector unsigned char)v2); 
			v3         = vec_add(v3,v7);
			v4         = vec_add(v4,v8);
			v5         = vec_add(v5,v9);
			v6         = vec_splat(v3,0);  /* Ox Ox Ox Ox */
			v7         = vec_splat(v3,1);  /* Oy Oy Oy Oy */
			v8         = vec_splat(v3,2);  /* Oz Oz Oz Oz */
			v9         = vec_splat(v3,3);  /* H1x H1x H1x H1x */
			v10        = vec_splat(v4,0);  /* H1y H1y H1y H1y */
			v11        = vec_splat(v4,1);  /* H1z H1z H1z H1z */
			v12        = vec_splat(v4,2);  /* H2x H2x H2x H2x */
			v13        = vec_splat(v4,3);  /* H2y H2y H2y H2y */
			v14        = vec_splat(v5,0);  /* H2z H2z H2z H2z */
			/* Store i 3atoms coordinates to stack */
			vec_st(v6,  80, (float *)stackdata); /* i Ox is in stack pos 5 */
			vec_st(v7,  96, (float *)stackdata); /* i Oy is in stack pos 6 */
			vec_st(v8, 112, (float *)stackdata); /* i Oz is in stack pos 7 */
			vec_st(v9, 128, (float *)stackdata); /* i H1x is in stack pos 8 */
			vec_st(v10,144, (float *)stackdata); /* i H1y is in stack pos 9 */
			vec_st(v11,160, (float *)stackdata); /* i H1z is in stack pos 10 */
			vec_st(v12,176, (float *)stackdata); /* i H2x is in stack pos 11 */
			vec_st(v13,192, (float *)stackdata); /* i H2y is in stack pos 12 */
			vec_st(v14,208, (float *)stackdata); /* i H2z is in stack pos 13 */

			nj0        = jindex[n];
			nj1        = jindex[n+1];
			/*			vec_dst( jjnr + nj1, 0x10010100, 0 ); */
			/* zero vctot, in stack pos 14 */
			vec_st(v0, 224, (float *)stackdata); 
			/* zero vctot, in stack pos 15 */
			vec_st(v0, 240, (float *)stackdata); 
			/* zero fiOx, in stack pos 16 */
			vec_st(v0, 256, (float *)stackdata); 
			/* zero fiOy, in stack pos 17 */
			vec_st(v0, 272, (float *)stackdata); 
			/* zero fiOz, in stack pos 18 */
			vec_st(v0, 288, (float *)stackdata); 
			/* zero fiH1x, in stack pos 19 */
			vec_st(v0, 304, (float *)stackdata); 
			/* zero fiH1y, in stack pos 20 */
			vec_st(v0, 320, (float *)stackdata); 
			/* zero fiH1z, in stack pos 21 */
			vec_st(v0, 336, (float *)stackdata); 
			/* zero fiH2x, in stack pos 22 */
			vec_st(v0, 352, (float *)stackdata); 
			/* zero fiH2y, in stack pos 23 */
			vec_st(v0, 368, (float *)stackdata); 
			/* zero fiH2z, in stack pos 24 */
			vec_st(v0, 384, (float *)stackdata); 

			for(k=nj0; k<(nj1-3); k+=4) { 
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				jnrd            = jjnr[k+3];

				/*				vec_dst( jjnr + k + 4, 0x02020020, 0 ); */

				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;
				j3d             = 3*jnrd;

				/*				vec_dst( pos+j3a, 0x10010100, 1 ); */

				v1              = (vector float)vec_lvsl(0, pos+j3a);
				v8              = (vector float)vec_lvsl(0, pos+j3b);
				v15             = (vector float)vec_lvsl(0, pos+j3c);
				v22             = (vector float)vec_lvsl(0, pos+j3d);
				v2              = vec_ld(0, pos+j3a);
				v9              = vec_ld(0, pos+j3b);
				v16             = vec_ld(0, pos+j3c);
				v23             = vec_ld(0, pos+j3d);

				v3              = vec_ld(16, pos+j3a);
				v10             = vec_ld(16, pos+j3b);
				v17             = vec_ld(16, pos+j3c);
				v24             = vec_ld(16, pos+j3d);
				v4              = vec_ld(32, pos+j3a);
				v11             = vec_ld(32, pos+j3b);
				v18             = vec_ld(32, pos+j3c);
				v25             = vec_ld(32, pos+j3d);
				/* v5  =  Oxa  Oya  Oza H1xa */
				v5              = vec_perm(v2,v3,(vector unsigned char)v1); 
				/* v12 = Oxb  Oyb  Ozb H1xb */
				v12             = vec_perm(v9,v10,(vector unsigned char)v8);  
				/* v19 = Oxc  Oyc  Ozc H1xc */
				v19             = vec_perm(v16,v17,(vector unsigned char)v15);
				/* v26 =  Oxd  Oyd  Ozd H1xd */
				v26             = vec_perm(v23,v24,(vector unsigned char)v22);
				/* H1ya H1za H2xa H2ya */
				v6              = vec_perm(v3,v4,(vector unsigned char)v1); 
				 /* H1yb H1zb H2xb H2yb */
				v13             = vec_perm(v10,v11,(vector unsigned char)v8);
				 /* H1yc H1zc H2xc H2yc */
				v20             = vec_perm(v17,v18,(vector unsigned char)v15);
				 /* H1yd H1zd H2xd H2yd */
				v27             = vec_perm(v24,v25,(vector unsigned char)v22);
				/* H2za   -   -   - */ 
				v7              = vec_perm(v4,v4,(vector unsigned char)v1);    
				 /* H2zb   -   -   - */
				v14             = vec_perm(v11,v11,(vector unsigned char)v8);
				 /* H2zc   -   -   - */
				v21             = vec_perm(v18,v18,(vector unsigned char)v15);
				 /* H2zd   -   -   - */
				v28             = vec_perm(v25,v25,(vector unsigned char)v22);
      
				/* permute 3atoms coordinates */
				/*  Oxa  Oxc  Oya  Oyc */
				v3              = vec_mergeh(v5,v19);  
				/*  Oza  Ozc H1xa H1xc */
				v5              = vec_mergel(v5,v19);  
				/*  Oxb  Oxd  Oyb  Oyd */
				v19             = vec_mergeh(v12,v26); 
				/*  Ozb  Ozd H1xb H1xd */
				v12             = vec_mergel(v12,v26); 
				/* H1ya H1yc H1za H1zc */
				v26             = vec_mergeh(v6,v20);  
				/* H2xa H2xc H2ya H2yc */
				v16              = vec_mergel(v6,v20);  
				/* H1yb H1yd H1zb H1zd */
				v20             = vec_mergeh(v13,v27); 
				/* H2xb H2xd H2yb H2yd */
				v13             = vec_mergel(v13,v27); 

				/* H2za H2zc   -    -  */
				v15             = vec_mergeh(v7,v21);  
				/* H2zb H2zd   -    -  */
				v14             = vec_mergeh(v14,v28); 

				/*  Oxa  Oxb  Oxc  Oxd */
				v1              = vec_mergeh(v3,v19);  
				/* load i H1x */
				v29             = vec_ld(128, (float *) stackdata); 
				/*  Oya  Oyb  Oyc  Oyd */
				v2              = vec_mergel(v3,v19);  
				/* load i H1y */
				v30             = vec_ld(144, (float *) stackdata); 
				/*  Oza  Ozb  Ozc  Ozd */
				v3              = vec_mergeh(v5,v12);  
				/* load i H1z */
				v31             = vec_ld(160, (float *) stackdata); 
				/* H1xa H1xb H1xc H1xd */
				v4              = vec_mergel(v5,v12);  
				/* H1ya H1yb H1yc H1yd */
				v5              = vec_mergeh(v26,v20); 
				/* H1za H1zb H1zc H1zd */
				v6              = vec_mergel(v26,v20); 
				/* H2xa H2xb H2xc H2xd */
				v7              = vec_mergeh(v16,v13); 
				/* H2ya H2yb H2yc H2yd */
				v8              = vec_mergel(v16,v13); 
				/* H2za H2zb H2zc H2zd */
				v9              = vec_mergeh(v15,v14); 

				v10             = vec_sub(v29,v1); /* iH1x - jOx */
				v13             = vec_sub(v29,v4); /* iH1x - jH1x */
				v16             = vec_sub(v29,v7); /* iH1x - jH2x */
				/* load i H2x */     
				v29             = vec_ld(176, (float *) stackdata); 
				v11             = vec_sub(v30,v2); /* iH1y - jOy */
				v14             = vec_sub(v30,v5); /* iH1y - jH1y */
				v17             = vec_sub(v30,v8); /* iH1y - jH2y */
				/* load i H2y */     
				v30             = vec_ld(192, (float *) stackdata); 
				vec_st(v10, 544, (float *)stackdata); /* dx21 */
				vec_st(v13, 592, (float *)stackdata); /* dx22 */
				vec_st(v16, 640, (float *)stackdata); /* dx23 */
				v12             = vec_sub(v31,v3); /* iH1z - jOz */
				v15             = vec_sub(v31,v6); /* iH1z - jH1z */
				v18             = vec_sub(v31,v9); /* iH1z - jH2z */
				/* load i H2z */         
				v31             = vec_ld(208, (float *) stackdata); 
				/* v10-v18 now contains iH1-jO, iH1-jH1 and iJ1-jH2 dist. */
				vec_st(v11, 560, (float *)stackdata); /* dy21 */
				vec_st(v14, 608, (float *)stackdata); /* dy22 */
				vec_st(v17, 656, (float *)stackdata); /* dy23 */
				v19             = vec_sub(v29,v1); /* iH2x - jOx */
				v22             = vec_sub(v29,v4); /* iH2x - jH1x */
				v25             = vec_sub(v29,v7); /* iH2x - jH2x */
				vec_st(v12, 576, (float *)stackdata); /* dz21 */
				vec_st(v15, 624, (float *)stackdata); /* dz22 */
				vec_st(v18, 672, (float *)stackdata); /* dz23 */
				/* load i Ox */     
				v29             = vec_ld(80, (float *) stackdata); 
				v20             = vec_sub(v30,v2); /* iH2y - jOy */
				v23             = vec_sub(v30,v5); /* iH2y - jH1y */
				v26             = vec_sub(v30,v8); /* iH2y - jH2y */
				vec_st(v19, 688, (float *)stackdata); /* dx31 */
				vec_st(v22, 736, (float *)stackdata); /* dx32 */
				vec_st(v25, 784, (float *)stackdata); /* dx33 */
				/* load i Oy */     
				v30             = vec_ld(96, (float *) stackdata); 
				v21             = vec_sub(v31,v3); /* iH2z - jOz */
				v24             = vec_sub(v31,v6); /* iH2z - jH1z */
				v27             = vec_sub(v31,v9); /* iH2z - jH2z */
				/* load i Oz */     
				v31             = vec_ld(112, (float *) stackdata); 
				vec_st(v20, 704, (float *)stackdata); /* dy31 */
				vec_st(v23, 752, (float *)stackdata); /* dy32 */
				vec_st(v26, 800, (float *)stackdata); /* dy33 */

				v1              = vec_sub(v29,v1); /* iOx - jOx */
				v4              = vec_sub(v29,v4); /* iOx - jH1x */
				v7              = vec_sub(v29,v7); /* iOx - jH2x */
				vec_st(v21, 720, (float *)stackdata); /* dz31 */
				vec_st(v24, 768, (float *)stackdata); /* dz32 */
				vec_st(v27, 816, (float *)stackdata); /* dz33 */
				v2              = vec_sub(v30,v2); /* iOy - jOy */
				v5              = vec_sub(v30,v5); /* iOy - jH1y */
				v8              = vec_sub(v30,v8); /* iOy - jH2y */
				vec_st(v1, 400, (float *)stackdata); /* dx11 */
				vec_st(v4, 448, (float *)stackdata); /* dx12 */
				vec_st(v7, 496, (float *)stackdata); /* dx13 */
				v3              = vec_sub(v31,v3); /* iOz - jOz */
				v6              = vec_sub(v31,v6); /* iOz - jH1z */
				v9              = vec_sub(v31,v9); /* iOz - jH2z */
				vec_st(v2, 416, (float *)stackdata); /* dy11 */
				vec_st(v5, 464, (float *)stackdata); /* dy12 */
				vec_st(v8, 512, (float *)stackdata); /* dy13 */

				v1              = vec_madd(v1,v1,v0);
				v4              = vec_madd(v4,v4,v0);
				v7              = vec_madd(v7,v7,v0);
				vec_st(v3, 432, (float *)stackdata); /* dz11 */
				vec_st(v6, 480, (float *)stackdata); /* dz12 */
				vec_st(v9, 528, (float *)stackdata); /* dz13 */
				v10             = vec_madd(v10,v10,v0);
				v13             = vec_madd(v13,v13,v0);
				v16             = vec_madd(v16,v16,v0);
				v19             = vec_madd(v19,v19,v0);
				v22             = vec_madd(v22,v22,v0);
				v25             = vec_madd(v25,v25,v0);
				v1              = vec_madd(v2,v2,v1);
				v4              = vec_madd(v5,v5,v4);
				v7              = vec_madd(v8,v8,v7);
				v10             = vec_madd(v11,v11,v10);
				v13             = vec_madd(v14,v14,v13);
				v16             = vec_madd(v17,v17,v16);
				v19             = vec_madd(v20,v20,v19);
				v22             = vec_madd(v23,v23,v22);
				v25             = vec_madd(v26,v26,v25);
				v1              = vec_madd(v3,v3,v1);
				v2              = vec_madd(v6,v6,v4);
				v3              = vec_madd(v9,v9,v7);
				v4              = vec_madd(v12,v12,v10);
				v5              = vec_madd(v15,v15,v13);
				v6              = vec_madd(v18,v18,v16);
				v7              = vec_madd(v21,v21,v19);
				v8              = vec_madd(v24,v24,v22);
				v9              = vec_madd(v27,v27,v25);
				/* 
				 * v1  = rsq  iO-jO
				 * v2  = rsq  iO-jH1
				 * v3  = rsq  iO-jH2
				 * v4  = rsq  iH1-jO
				 * v5  = rsq  iH1-jH1
				 * v6  = rsq  iH1-jH2
				 * v7  = rsq  iH2-jO
				 * v8  = rsq  iH2-jH1
				 * v9 = rsq  iH2-jH2
				 */

				v10             = vec_rsqrte(v1);
				v11             = vec_rsqrte(v2);
				v12             = vec_rsqrte(v3);
				v13             = vec_rsqrte(v4);
				v14             = vec_rsqrte(v5);
				v15             = vec_rsqrte(v6);
				v16             = vec_rsqrte(v7);
				v17             = vec_rsqrte(v8);
				v18             = vec_rsqrte(v9);
				/* create constant 0.5 */
				v30             = (vector float) vec_splat_u32(1);
				/* 0.5 */
				v31             = vec_ctf((vector unsigned int)v30,1); 
				/* 1.0 */
				v30             = vec_ctf((vector unsigned int)v30,0); 

				v19             = vec_madd(v10,v10,v0); /* lu*lu */
				v20             = vec_madd(v11,v11,v0);
				v21             = vec_madd(v12,v12,v0);
				v22             = vec_madd(v13,v13,v0);
				v23             = vec_madd(v14,v14,v0);
				v24             = vec_madd(v15,v15,v0);
				v25             = vec_madd(v16,v16,v0);
				v26             = vec_madd(v17,v17,v0);
				v27             = vec_madd(v18,v18,v0);

				v19             = vec_nmsub(v1,v19,v30); /* 1.0 - rsq*lu*lu */
				v20             = vec_nmsub(v2,v20,v30);
				v21             = vec_nmsub(v3,v21,v30);
				v22             = vec_nmsub(v4,v22,v30);
				v23             = vec_nmsub(v5,v23,v30);
				v24             = vec_nmsub(v6,v24,v30);
				v25             = vec_nmsub(v7,v25,v30);
				v26             = vec_nmsub(v8,v26,v30);
				v27             = vec_nmsub(v9,v27,v30);

				v1              = vec_madd(v10,v31,v0);/* lu*0.5*/
				v2              = vec_madd(v11,v31,v0);
				v3              = vec_madd(v12,v31,v0);
				v4              = vec_madd(v13,v31,v0);
				v5              = vec_madd(v14,v31,v0);
				v6              = vec_madd(v15,v31,v0);
				v7              = vec_madd(v16,v31,v0);
				v8              = vec_madd(v17,v31,v0);
				v9              = vec_madd(v18,v31,v0);

				/* The rinv values */
				v1              = vec_madd(v1,v19,v10);
				v2              = vec_madd(v2,v20,v11);
				v3              = vec_madd(v3,v21,v12);
				v4              = vec_madd(v4,v22,v13);
				v5              = vec_madd(v5,v23,v14);
				v6              = vec_madd(v6,v24,v15);
				v7              = vec_madd(v7,v25,v16);
				v8              = vec_madd(v8,v26,v17);
				v9              = vec_madd(v9,v27,v18);
      
				/* load qqOO, qqOH and qqHH  to v27,v28,v29 */
				v27             = vec_ld(0, (float *) stackdata);
				v28             = vec_ld(16, (float *) stackdata);
				v29             = vec_ld(32, (float *) stackdata);

				/*				vec_dstst( faction+j3a, 0x10010100, 2 ); */

				/* put rinvsq in v10-v18, rinv6_OO in v30 & rinv12_OO in v31 */
				/* load c6 to v25 and c12 to v26 */
				v25             = vec_ld(48, (float *) stackdata);
				v26             = vec_ld(64, (float *) stackdata);
      
				v10             = vec_madd(v1,v1,v0);
				v1              = vec_madd(v1,v27,v0); /* rinv11*qqOO */
				v11             = vec_madd(v2,v2,v0);
				/* load vctot to v23 and Vvdwtot to v24 */
				v23             = vec_ld(224,(float *) stackdata);
				v24             = vec_ld(240,(float *) stackdata);

				v2              = vec_madd(v2,v28,v0); /* rinv12*qqOH */
				v12             = vec_madd(v3,v3,v0);
				v30             = vec_madd(v10,v10,v0); /* rinv4 */
				v3              = vec_madd(v3,v28,v0); /* rinv13*qqOH */
				v13             = vec_madd(v4,v4,v0);
				v4              = vec_madd(v4,v28,v0); /* rinv21*qqOH */
				v14             = vec_madd(v5,v5,v0);

				v23             = vec_add(v23,v1);

				v30             = vec_madd(v30,v10,v0); /* rinv6 */
				v5              = vec_madd(v5,v29,v0); /* rinv22*qqHH */
				v15             = vec_madd(v6,v6,v0);
				v6              = vec_madd(v6,v29,v0); /* rinv23*qqHH */
				v23             = vec_add(v23,v2);
				v16             = vec_madd(v7,v7,v0);
				v31             = vec_madd(v30,v30,v0); /* rinv12 */
				v25             = vec_madd(v25,v30,v0); /* c6*rinv6 */
				/* load 6.0 to v30 */
				v30             = (vector float)vec_splat_u32(6);
				v30             = vec_ctf((vector unsigned int)v30,0);
				v23             = vec_add(v23,v3);

				v7              = vec_madd(v7,v28,v0); /* rinv31*qqOH */
				v17             = vec_madd(v8,v8,v0);
				v8              = vec_madd(v8,v29,v0); /* rinv32*qqHH */
				v26             = vec_madd(v26,v31,v0); /* c12*rinv12 */
				v23             = vec_add(v23,v4);
				/* load 12.0 to v31 */
				v31             = (vector float)vec_splat_u32(12);
				v31             = vec_ctf((vector unsigned int)v31,0);

				v24             = vec_sub(v24,v25);  /* add Vvdw6 to Vvdwtot */
				v18             = vec_madd(v9,v9,v0);
				v23             = vec_add(v23,v5);
				v9              = vec_madd(v9,v29,v0); /* rinv33*qqHH */

				v24             = vec_add(v24,v26);/* add Vvdw12 to Vvdwtot */
    
				v31             = vec_madd(v31,v26,v0);
				v11             = vec_madd(v11,v2,v0); /* fs12 */
				v23             = vec_add(v23,v6);
				v12             = vec_madd(v12,v3,v0); /* fs13 */
				v13             = vec_madd(v13,v4,v0); /* fs21 */
				v31             = vec_nmsub(v30,v25,v31);

				v14             = vec_madd(v14,v5,v0); /* fs22 */
				v23             = vec_add(v23,v7);
				v15             = vec_madd(v15,v6,v0); /* fs23 */
				v16             = vec_madd(v16,v7,v0); /* fs31 */
				v1              = vec_add(v31,v1);
				v17             = vec_madd(v17,v8,v0); /* fs32 */
				v23             = vec_add(v23,v8);
				v18             = vec_madd(v18,v9,v0); /* fs33 */
				v10             = vec_madd(v10,v1,v0);

				vec_st(v24,240,(float *)stackdata); /* store Vvdwtot */
				/* calculate vectorial forces and accumulate fj. 
				 * v10-v18 has fs11-fs33 now. 
				 * First load iO-* dx,dy,dz vectors to v1-v9 
				 * and load iO forces to v28,v29,v30 
				 * use v19-v27 to accumulate j 3atoms forces 
				 */
				v28             = vec_ld(256, (float *) stackdata);
				v29             = vec_ld(272, (float *) stackdata);
				v30             = vec_ld(288, (float *) stackdata);

				v1              = vec_ld(400, (float *) stackdata);
				v2              = vec_ld(416, (float *) stackdata);
				v23             = vec_add(v23,v9); /* incr. vctot */
				v3              = vec_ld(432, (float *) stackdata);
				v4              = vec_ld(448, (float *) stackdata);
				v5              = vec_ld(464, (float *) stackdata);
				v6              = vec_ld(480, (float *) stackdata);
				vec_st(v23,224,(float *)stackdata); /* store vctot to stack */
				v7              = vec_ld(496, (float *) stackdata);
				v8              = vec_ld(512, (float *) stackdata);
				v9              = vec_ld(528, (float *) stackdata);

				v28             = vec_madd(v10,v1,v28);
				v19             = vec_nmsub(v10,v1,v0);
				v29             = vec_madd(v10,v2,v29);
				v20             = vec_nmsub(v10,v2,v0);
				v30             = vec_madd(v10,v3,v30);
				v21             = vec_nmsub(v10,v3,v0);

				v28             = vec_madd(v11,v4,v28);
				v22             = vec_nmsub(v11,v4,v0);
				v29             = vec_madd(v11,v5,v29);
				v23             = vec_nmsub(v11,v5,v0);
				v30             = vec_madd(v11,v6,v30);
				v24             = vec_nmsub(v11,v6,v0);

				v28             = vec_madd(v12,v7,v28);
				v25             = vec_nmsub(v12,v7,v0);
				v29             = vec_madd(v12,v8,v29);
				v26             = vec_nmsub(v12,v8,v0);
				v30             = vec_madd(v12,v9,v30);
				v27             = vec_nmsub(v12,v9,v0);

				/* store these i forces, and repeat the 
				 * procedure for the iH1-* force 
				 */
				vec_st(v28,256,(float *)stackdata);
				vec_st(v29,272,(float *)stackdata);
				vec_st(v30,288,(float *)stackdata);

				v28             = vec_ld(304,(float *) stackdata);
				v29             = vec_ld(320,(float *) stackdata);
				v30             = vec_ld(336,(float *) stackdata);
				/* load new vectorial distances */
				v1              = vec_ld(544, (float *) stackdata);
				v2              = vec_ld(560, (float *) stackdata);
				v3              = vec_ld(576, (float *) stackdata);
				v4              = vec_ld(592, (float *) stackdata);
				v5              = vec_ld(608, (float *) stackdata);
				v6              = vec_ld(624, (float *) stackdata);
				v7              = vec_ld(640, (float *) stackdata);
				v8              = vec_ld(656, (float *) stackdata);
				v9              = vec_ld(672, (float *) stackdata);
      
				v28             = vec_madd(v13,v1,v28);
				v19             = vec_nmsub(v13,v1,v19);
				v29             = vec_madd(v13,v2,v29);
				v20             = vec_nmsub(v13,v2,v20);
				v30             = vec_madd(v13,v3,v30);
				v21             = vec_nmsub(v13,v3,v21);

				v28             = vec_madd(v14,v4,v28);
				v22             = vec_nmsub(v14,v4,v22);
				v29             = vec_madd(v14,v5,v29);
				v23             = vec_nmsub(v14,v5,v23);
				v30             = vec_madd(v14,v6,v30);
				v24             = vec_nmsub(v14,v6,v24);

				v28             = vec_madd(v15,v7,v28);
				v25             = vec_nmsub(v15,v7,v25);
				v29             = vec_madd(v15,v8,v29);
				v26             = vec_nmsub(v15,v8,v26);
				v30             = vec_madd(v15,v9,v30);
				v27             = vec_nmsub(v15,v9,v27);

				/* store these i forces, and repeat the
				 * procedure for the iH2-* force 
				 */
				vec_st(v28,304,(float *)stackdata);
				vec_st(v29,320,(float *)stackdata);
				vec_st(v30,336,(float *)stackdata);
				v28             = vec_ld(352,(float *) stackdata);
				v29             = vec_ld(368,(float *) stackdata);
				v30             = vec_ld(384,(float *) stackdata);
				/* load new vectorial distances */
				v1              = vec_ld(688, (float *) stackdata);
				v2              = vec_ld(704, (float *) stackdata);
				v3              = vec_ld(720, (float *) stackdata);
				v4              = vec_ld(736, (float *) stackdata);
				v5              = vec_ld(752, (float *) stackdata);
				v6              = vec_ld(768, (float *) stackdata);
				v7              = vec_ld(784, (float *) stackdata);
				v8              = vec_ld(800, (float *) stackdata);
				v9              = vec_ld(816, (float *) stackdata);
      
				v28             = vec_madd(v16,v1,v28);
				v19             = vec_nmsub(v16,v1,v19);
				v29             = vec_madd(v16,v2,v29);
				v20             = vec_nmsub(v16,v2,v20);
				v30             = vec_madd(v16,v3,v30);
				v21             = vec_nmsub(v16,v3,v21);

				v28             = vec_madd(v17,v4,v28);
				v22             = vec_nmsub(v17,v4,v22);
				v29             = vec_madd(v17,v5,v29);
				v23             = vec_nmsub(v17,v5,v23);
				v30             = vec_madd(v17,v6,v30);
				v24             = vec_nmsub(v17,v6,v24);

				v28             = vec_madd(v18,v7,v28);
				v25             = vec_nmsub(v18,v7,v25);
				v29             = vec_madd(v18,v8,v29);
				v26             = vec_nmsub(v18,v8,v26);
				v30             = vec_madd(v18,v9,v30);
				v27             = vec_nmsub(v18,v9,v27);

				/* store these i forces */
				vec_st(v28,352,(float *)stackdata);
				vec_st(v29,368,(float *)stackdata);
				vec_st(v30,384,(float *)stackdata);

				/* j forces present in v19-v27 */    
				
				/*  Oxa  Oza  Oxb  Ozb */
				v1              = vec_mergeh(v19,v21); 
				/*  Oxc  Ozc  Oxd  Ozd */
				v19             = vec_mergel(v19,v21); 
				/*  Oya H1xa  Oyb H1xb */
				v21             = vec_mergeh(v20,v22); 
				/*  Oyc H1xc  Oyd H1xd */
				v20             = vec_mergel(v20,v22); 
				/* H1ya H2xa H1yb H2xb */
				v22             = vec_mergeh(v23,v25); 
				/* H1yc H2xc H1yd H2xd */
				v23             = vec_mergel(v23,v25); 
				/* H1za H2ya H1zb H2yb */
				v25             = vec_mergeh(v24,v26); 
				/* H1zc H2yc H1zd H2yd */
				v24             = vec_mergel(v24,v26); 

				/* H2za   0  H2zb   0  */
				v26             = vec_mergeh(v27,v0);   
				/* H2zc   0  H2zd   0  */
				v27             = vec_mergel(v27,v0);   
      
				/*  Oxa  Oya  Oza H1xa */
				v2              = vec_mergeh(v1,v21);   
				/*  Oxb  Oyb  Ozb H1xb */
				v21             = vec_mergel(v1,v21);   
				/*  Oxc  Oyc  Ozc H1xc */
				v1              = vec_mergeh(v19,v20);    
				/*  Oxd  Oyd  Ozd H1xd */
				v19             = vec_mergel(v19,v20);    
				/* H1ya H1za H2xa H2ya */
				v20             = vec_mergeh(v22,v25);  
				/* H1yb H1zb H2xb H2yb */
				v22             = vec_mergel(v22,v25);  
				/* H1yc H1zc H2xc H2yc */
				v25             = vec_mergeh(v23,v24);  
				/* H1yd H1zd H2xd H2yd */
				v23             = vec_mergel(v23,v24);  
				/* H2za   0    0    0  */
				v24             = vec_mergeh(v26,v0);   
				/* H2zb   0    0    0  */
				v26             = vec_mergel(v26,v0);   
				/* H2zc   0    0    0  */
				v3              = vec_mergeh(v27,v0);   
				/* H2zd   0    0    0  */
				v27             = vec_mergel(v27,v0);   
 
				v29             = (vector float)vec_splat_s32(-1);

				/* move into position, load and add */  
				v30            = (vector float)vec_lvsr(0,(int *)faction+j3a); 
				v4              = vec_ld(  0, faction+j3a);
				v6              = vec_ld( 16, faction+j3a);
				v8              = vec_ld( 32, faction+j3a);
				v10             = vec_perm(v0,v29,(vector unsigned char)v30);
				v12             = vec_perm(v0,v2,(vector unsigned char)v30);
				v12             = vec_add(v12,v4);
				v14             = vec_perm(v2,v20,(vector unsigned char)v30);
				v2              = vec_add(v14,v6);
				v16             = vec_perm(v20,v24,(vector unsigned char)v30);
				v20             = vec_add(v16,v8);
				v12             = vec_sel(v4,v12,(vector unsigned int)v10);
				vec_st(v12,  0, faction+j3a);
				v10             = vec_sld(v0,v10,12);
				vec_st(v2, 16, faction+j3a);
				v12             = vec_sel(v20,v8,(vector unsigned int)v10);
				vec_st(v12, 32, faction+j3a);


				/* Finished 1, now do 2  */
				v30            = (vector float)vec_lvsr(0,(int *)faction+j3b); 
				v4              = vec_ld(  0, faction+j3b);
				v6              = vec_ld( 16, faction+j3b);
				v8              = vec_ld( 32, faction+j3b);
				v10             = vec_perm(v0,v29,(vector unsigned char)v30);
				v12             = vec_perm(v0,v21,(vector unsigned char)v30);
				v24             = vec_add(v12,v4);
				v12             = vec_perm(v21,v22,(vector unsigned char)v30);
				v21             = vec_add(v12,v6);
				v12             = vec_perm(v22,v26,(vector unsigned char)v30);
				v22             = vec_add(v12,v8);
				v12             = vec_sel(v4,v24,(vector unsigned int)v10);
				vec_st(v12,  0, faction+j3b);
				v10             = vec_sld(v0,v10,12);
				vec_st(v21, 16, faction+j3b);
				v12             = vec_sel(v22,v8,(vector unsigned int)v10);
				vec_st(v12, 32, faction+j3b);

				/* water 3 */
				v31            = (vector float)vec_lvsr(0,(int *)faction+j3c);
				v5              = vec_ld(  0, faction+j3c);
				v7              = vec_ld( 16, faction+j3c);
				v9              = vec_ld( 32, faction+j3c);
				v11             = vec_perm(v0,v29,(vector unsigned char)v31);
				v13             = vec_perm(v0,v1,(vector unsigned char)v31);
				v13             = vec_add(v13,v5);
				v15             = vec_perm(v1,v25,(vector unsigned char)v31);
				v1              = vec_add(v15,v7);      
				v17             = vec_perm(v25,v3,(vector unsigned char)v31);
				v25             = vec_add(v17,v9);
				v13             = vec_sel(v5,v13,(vector unsigned int)v11);
				vec_st(v13,  0, faction+j3c);
				v11             = vec_sld(v0,v11,12);
				vec_st(v1, 16, faction+j3c);
				v13             = vec_sel(v25,v9,(vector unsigned int)v11);
				vec_st(v13, 32, faction+j3c);

				/* water 4 */
				v31            = (vector float)vec_lvsr(0,(int *)faction+j3d);    
				v5              = vec_ld(  0, faction+j3d);
				v7              = vec_ld( 16, faction+j3d);
				v9              = vec_ld( 32, faction+j3d);
				v11             = vec_perm(v0,v29,(vector unsigned char)v31);      
				v13             = vec_perm(v0,v19,(vector unsigned char)v31);
				v25             = vec_add(v13,v5);
				v13             = vec_perm(v19,v23,(vector unsigned char)v31);
				v19             = vec_add(v13,v7);
				v13             = vec_perm(v23,v27,(vector unsigned char)v31);
				v23             = vec_add(v13,v9);
				v13             = vec_sel(v5,v25,(vector unsigned int)v11);
				vec_st(v13,  0, faction+j3d);
				v11             = vec_sld(v0,v11,12);
				vec_st(v19, 16, faction+j3d);
				v13             = vec_sel(v23,v9,(vector unsigned int)v11);
				vec_st(v13, 32, faction+j3d);
			} 
			if(k<(nj1-2)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				jnrc            = jjnr[k+2];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;
				j3c             = 3*jnrc;

				v1              = (vector float)vec_lvsl(0, pos+j3a);
				v8              = (vector float)vec_lvsl(0, pos+j3b);
				v15             = (vector float)vec_lvsl(0, pos+j3c);

				v2              = vec_ld(0, pos+j3a);
				v9              = vec_ld(0, pos+j3b);
				v16             = vec_ld(0, pos+j3c);
				v3              = vec_ld(16, pos+j3a);
				v10             = vec_ld(16, pos+j3b);
				v17             = vec_ld(16, pos+j3c);
				v4              = vec_ld(32, pos+j3a);
				v11             = vec_ld(32, pos+j3b);
				v18             = vec_ld(32, pos+j3c);
				/*  Oxa  Oya  Oza H1xa */
				v5              = vec_perm(v2,v3,(vector unsigned char)v1); 
				 /*  Oxb  Oyb  Ozb H1xb */
				v12             = vec_perm(v9,v10,(vector unsigned char)v8); 
				 /*  Oxc  Oyc  Ozc H1xc */
				v19             = vec_perm(v16,v17,(vector unsigned char)v15);

				 /* H1ya H1za H2xa H2ya */
				v6              = vec_perm(v3,v4,(vector unsigned char)v1);
				 /* H1yb H1zb H2xb H2yb */
				v13             = vec_perm(v10,v11,(vector unsigned char)v8);
				 /* H1yc H1zc H2xc H2yc */
				v20             = vec_perm(v17,v18,(vector unsigned char)v15);
				/* H2za   -   -   - */    
				v7              = vec_perm(v4,v4,(vector unsigned char)v1); 
				 /* H2zb   -   -   - */
				v14             = vec_perm(v11,v11,(vector unsigned char)v8);
				/* H2zc   -   -   - */
				v21             = vec_perm(v18,v18,(vector unsigned char)v15); 
      
				/* permute 3atoms coordinates */
				/*  Oxa  Oxc  Oya  Oyc */
				v3              = vec_mergeh(v5,v19);  
				v5              = vec_mergel(v5,v19); /*  Oza  Ozc H1xa H1xc */
				v19             = vec_mergeh(v12,v0); /*  Oxb   -   Oyb   -  */
				v12             = vec_mergel(v12,v0); /*  Ozb   -  H1xb   -  */
      
				v26             = vec_mergeh(v6,v20); /* H1ya H1yc H1za H1zc */
				v16              = vec_mergel(v6,v20);/* H2xa H2xc H2ya H2yc */
				/* H1yb   -  H1zb  -  */
				v20             = vec_mergeh(v13,v0); 
				/* H2xb   -  H2yb  -  */
				v13             = vec_mergel(v13,v0); 
				/* H2za H2zc   -    -  */
				v15             = vec_mergeh(v7,v21);  

				/*  Oxa  Oxb  Oxc  -  */
				v1              = vec_mergeh(v3,v19);  
				/* load i H1x */
				v29             = vec_ld(128, (float *) stackdata); 
				v2              = vec_mergel(v3,v19);  /*  Oya  Oyb  Oyc  -  */
				 /* load i H1y */
				v30             = vec_ld(144, (float *) stackdata);
				v3              = vec_mergeh(v5,v12);  /*  Oza  Ozb  Ozc  - */
				/* load i H1z */
				v31             = vec_ld(160, (float *) stackdata); 
				v4              = vec_mergel(v5,v12);  /* H1xa H1xb H1xc  -  */
				v5              = vec_mergeh(v26,v20); /* H1ya H1yb H1yc  -  */
				v6              = vec_mergel(v26,v20); /* H1za H1zb H1zc  -  */
				v7              = vec_mergeh(v16,v13); /* H2xa H2xb H2xc  -  */
				v8              = vec_mergel(v16,v13); /* H2ya H2yb H2yc  -  */
				v9              = vec_mergeh(v15,v14); /* H2za H2zb H2zc  -  */

				v10             = vec_sub(v29,v1); /* iH1x - jOx */
				v13             = vec_sub(v29,v4); /* iH1x - jH1x */
				v16             = vec_sub(v29,v7); /* iH1x - jH2x */
				/* load i H2x */     
				v29             = vec_ld(176, (float *) stackdata); 
				v11             = vec_sub(v30,v2); /* iH1y - jOy */
				v14             = vec_sub(v30,v5); /* iH1y - jH1y */
				v17             = vec_sub(v30,v8); /* iH1y - jH2y */
				/* load i H2y */     
				v30             = vec_ld(192, (float *) stackdata);
				vec_st(v10, 544, (float *)stackdata); /* dx21 */
				vec_st(v13, 592, (float *)stackdata); /* dx22 */
				vec_st(v16, 640, (float *)stackdata); /* dx23 */
				v12             = vec_sub(v31,v3); /* iH1z - jOz */
				v15             = vec_sub(v31,v6); /* iH1z - jH1z */
				v18             = vec_sub(v31,v9); /* iH1z - jH2z */
				/* load i H2z */         
				v31             = vec_ld(208, (float *) stackdata); 
				/* v10-v18 now contains iH1-jO, iH1-jH1 and iJ1-jH2 dist. */
				vec_st(v11, 560, (float *)stackdata); /* dy21 */
				vec_st(v14, 608, (float *)stackdata); /* dy22 */
				vec_st(v17, 656, (float *)stackdata); /* dy23 */
				v19             = vec_sub(v29,v1); /* iH2x - jOx */
				v22             = vec_sub(v29,v4); /* iH2x - jH1x */
				v25             = vec_sub(v29,v7); /* iH2x - jH2x */
				vec_st(v12, 576, (float *)stackdata); /* dz21 */
				vec_st(v15, 624, (float *)stackdata); /* dz22 */
				vec_st(v18, 672, (float *)stackdata); /* dz23 */
				/* load i Ox */     
				v29             = vec_ld(80, (float *) stackdata); 
				v20             = vec_sub(v30,v2); /* iH2y - jOy */
				v23             = vec_sub(v30,v5); /* iH2y - jH1y */
				v26             = vec_sub(v30,v8); /* iH2y - jH2y */
				vec_st(v19, 688, (float *)stackdata); /* dx31 */
				vec_st(v22, 736, (float *)stackdata); /* dx32 */
				vec_st(v25, 784, (float *)stackdata); /* dx33 */
				/* load i Oy */     
				v30             = vec_ld(96, (float *) stackdata); 
				v21             = vec_sub(v31,v3); /* iH2z - jOz */
				v24             = vec_sub(v31,v6); /* iH2z - jH1z */
				v27             = vec_sub(v31,v9); /* iH2z - jH2z */
				/* load i Oz */     
				v31             = vec_ld(112, (float *) stackdata); 
				vec_st(v20, 704, (float *)stackdata); /* dy31 */
				vec_st(v23, 752, (float *)stackdata); /* dy32 */
				vec_st(v26, 800, (float *)stackdata); /* dy33 */

				v1              = vec_sub(v29,v1); /* iOx - jOx */
				v4              = vec_sub(v29,v4); /* iOx - jH1x */
				v7              = vec_sub(v29,v7); /* iOx - jH2x */
				vec_st(v21, 720, (float *)stackdata); /* dz31 */
				vec_st(v24, 768, (float *)stackdata); /* dz32 */
				vec_st(v27, 816, (float *)stackdata); /* dz33 */
				v2              = vec_sub(v30,v2); /* iOy - jOy */
				v5              = vec_sub(v30,v5); /* iOy - jH1y */
				v8              = vec_sub(v30,v8); /* iOy - jH2y */
				vec_st(v1, 400, (float *)stackdata); /* dx11 */
				vec_st(v4, 448, (float *)stackdata); /* dx12 */
				vec_st(v7, 496, (float *)stackdata); /* dx13 */
				v3              = vec_sub(v31,v3); /* iOz - jOz */
				v6              = vec_sub(v31,v6); /* iOz - jH1z */
				v9              = vec_sub(v31,v9); /* iOz - jH2z */
				vec_st(v2, 416, (float *)stackdata); /* dy11 */
				vec_st(v5, 464, (float *)stackdata); /* dy12 */
				vec_st(v8, 512, (float *)stackdata); /* dy13 */

				v1              = vec_madd(v1,v1,v0);
				v4              = vec_madd(v4,v4,v0);
				v7              = vec_madd(v7,v7,v0);
				vec_st(v3, 432, (float *)stackdata); /* dz11 */
				vec_st(v6, 480, (float *)stackdata); /* dz12 */
				vec_st(v9, 528, (float *)stackdata); /* dz13 */
				v10             = vec_madd(v10,v10,v0);
				v13             = vec_madd(v13,v13,v0);
				v16             = vec_madd(v16,v16,v0);
				v19             = vec_madd(v19,v19,v0);
				v22             = vec_madd(v22,v22,v0);
				v25             = vec_madd(v25,v25,v0);
				v1              = vec_madd(v2,v2,v1);
				v4              = vec_madd(v5,v5,v4);
				v7              = vec_madd(v8,v8,v7);
				v10             = vec_madd(v11,v11,v10);
				v13             = vec_madd(v14,v14,v13);
				v16             = vec_madd(v17,v17,v16);
				v19             = vec_madd(v20,v20,v19);
				v22             = vec_madd(v23,v23,v22);
				v25             = vec_madd(v26,v26,v25);
				v1              = vec_madd(v3,v3,v1);
				v2              = vec_madd(v6,v6,v4);
				v3              = vec_madd(v9,v9,v7);
				v4              = vec_madd(v12,v12,v10);
				v5              = vec_madd(v15,v15,v13);
				v6              = vec_madd(v18,v18,v16);
				v7              = vec_madd(v21,v21,v19);
				v8              = vec_madd(v24,v24,v22);
				v9              = vec_madd(v27,v27,v25);
				/* 
				 * v1  = rsq  iO-jO
				 * v2  = rsq  iO-jH1
				 * v3  = rsq  iO-jH2
				 * v4  = rsq  iH1-jO
				 * v5  = rsq  iH1-jH1
				 * v6  = rsq  iH1-jH2
				 * v7  = rsq  iH2-jO
				 * v8  = rsq  iH2-jH1
				 * v9 = rsq  iH2-jH2
				 */

				v10             = vec_rsqrte(v1);
				v11             = vec_rsqrte(v2);
				v12             = vec_rsqrte(v3);
				v13             = vec_rsqrte(v4);
				v14             = vec_rsqrte(v5);
				v15             = vec_rsqrte(v6);
				v16             = vec_rsqrte(v7);
				v17             = vec_rsqrte(v8);
				v18             = vec_rsqrte(v9);

				/* create constant 0.5 */
				v30             = (vector float) vec_splat_u32(1);
				/* 0.5 */
				v31             = vec_ctf((vector unsigned int)v30,1); 
				/* 1.0 */
				v30             = vec_ctf((vector unsigned int)v30,0); 

				v19             = vec_madd(v10,v10,v0); /* lu*lu */
				v20             = vec_madd(v11,v11,v0);
				v21             = vec_madd(v12,v12,v0);
				v22             = vec_madd(v13,v13,v0);
				v23             = vec_madd(v14,v14,v0);
				v24             = vec_madd(v15,v15,v0);
				v25             = vec_madd(v16,v16,v0);
				v26             = vec_madd(v17,v17,v0);
				v27             = vec_madd(v18,v18,v0);

				v19             = vec_nmsub(v1,v19,v30); /* 1.0 - rsq*lu*lu */
				v20             = vec_nmsub(v2,v20,v30);
				v21             = vec_nmsub(v3,v21,v30);
				v22             = vec_nmsub(v4,v22,v30);
				v23             = vec_nmsub(v5,v23,v30);
				v24             = vec_nmsub(v6,v24,v30);
				v25             = vec_nmsub(v7,v25,v30);
				v26             = vec_nmsub(v8,v26,v30);
				v27             = vec_nmsub(v9,v27,v30);

				v1              = vec_madd(v10,v31,v0);/* lu*0.5*/
				v2              = vec_madd(v11,v31,v0);
				v3              = vec_madd(v12,v31,v0);
				v4              = vec_madd(v13,v31,v0);
				v5              = vec_madd(v14,v31,v0);
				v6              = vec_madd(v15,v31,v0);
				v7              = vec_madd(v16,v31,v0);
				v8              = vec_madd(v17,v31,v0);
				v9              = vec_madd(v18,v31,v0);

				/* The rinv values */
				v1              = vec_madd(v1,v19,v10);
				v2              = vec_madd(v2,v20,v11);
				v3              = vec_madd(v3,v21,v12);
				v4              = vec_madd(v4,v22,v13);
				v5              = vec_madd(v5,v23,v14);
				v6              = vec_madd(v6,v24,v15);
				v7              = vec_madd(v7,v25,v16);
				v8              = vec_madd(v8,v26,v17);
				v9              = vec_madd(v9,v27,v18);
      
				v10             = (vector float)vec_splat_s32(-1);
				v10             = vec_sld(v0,v10,4);

				v1           = (vector float)vec_sel((vector unsigned int)v1,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v2           = (vector float)vec_sel((vector unsigned int)v2,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v3           = (vector float)vec_sel((vector unsigned int)v3,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v4           = (vector float)vec_sel((vector unsigned int)v4,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v5           = (vector float)vec_sel((vector unsigned int)v5,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v6           = (vector float)vec_sel((vector unsigned int)v6,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v7           = (vector float)vec_sel((vector unsigned int)v7,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v8           = (vector float)vec_sel((vector unsigned int)v8,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v9           = (vector float)vec_sel((vector unsigned int)v9,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);

				/* load qqOO, qqOH and qqHH  to v27,v28,v29 */
				v27             = vec_ld(0, (float *) stackdata);
				v28             = vec_ld(16, (float *) stackdata);
				v29             = vec_ld(32, (float *) stackdata);
      

				/*				vec_dstst( faction+j3a, 0x10010100, 2 ); */

				v27             = vec_sld(v27,v0,4);
				v28             = vec_sld(v28,v0,4);
				v29             = vec_sld(v29,v0,4);

				/* put rinvsq in v10-v18, rinv6_OO in v30 
                 *& and rinv12_OO in v31 
				 */
				/* load c6 to v25 and c12 to v26 */
				v25             = vec_ld(48, (float *) stackdata);
				v26             = vec_ld(64, (float *) stackdata);
      
				v10             = vec_madd(v1,v1,v0);
				v1              = vec_madd(v1,v27,v0); /* rinv11*qqOO */
				v11             = vec_madd(v2,v2,v0);
				/* load vctot to v23 and Vvdwtot to v24 */
				v23             = vec_ld(224,(float *) stackdata);
				v24             = vec_ld(240,(float *) stackdata);

				v25             = vec_sld(v25,v0,4);
				v26             = vec_sld(v26,v0,4);

				v2              = vec_madd(v2,v28,v0); /* rinv12*qqOH */
				v12             = vec_madd(v3,v3,v0);
				v30             = vec_madd(v10,v10,v0); /* rinv4 */
				v3              = vec_madd(v3,v28,v0); /* rinv13*qqOH */
				v13             = vec_madd(v4,v4,v0);
				v4              = vec_madd(v4,v28,v0); /* rinv21*qqOH */
				v14             = vec_madd(v5,v5,v0);

				v23             = vec_add(v23,v1);

				v30             = vec_madd(v30,v10,v0); /* rinv6 */
				v5              = vec_madd(v5,v29,v0); /* rinv22*qqHH */
				v15             = vec_madd(v6,v6,v0);
				v6              = vec_madd(v6,v29,v0); /* rinv23*qqHH */
				v23             = vec_add(v23,v2);
				v16             = vec_madd(v7,v7,v0);
				v31             = vec_madd(v30,v30,v0); /* rinv12 */
				v25             = vec_madd(v25,v30,v0); /* c6*rinv6 */
				/* load 6.0 to v30 */
				v30             = (vector float)vec_splat_u32(6);
				v30             = vec_ctf((vector unsigned int)v30,0);
				v23             = vec_add(v23,v3);

				v7              = vec_madd(v7,v28,v0); /* rinv31*qqOH */
				v17             = vec_madd(v8,v8,v0);
				v8              = vec_madd(v8,v29,v0); /* rinv32*qqHH */
				v26             = vec_madd(v26,v31,v0); /* c12*rinv12 */
				v23             = vec_add(v23,v4);
				/* load 12.0 to v31 */
				v31             = (vector float)vec_splat_u32(12);
				v31             = vec_ctf((vector unsigned int)v31,0);


				v24             = vec_sub(v24,v25);  /* add Vvdw6 to Vvdwtot */
				v18             = vec_madd(v9,v9,v0);
				v23             = vec_add(v23,v5);
				v9              = vec_madd(v9,v29,v0); /* rinv33*qqHH */
				v24             = vec_add(v24,v26);/* add Vvdw12 to Vvdwtot */
    
				v31             = vec_madd(v31,v26,v0);
				v11             = vec_madd(v11,v2,v0); /* fs12 */
				v23             = vec_add(v23,v6);
				v12             = vec_madd(v12,v3,v0); /* fs13 */
				v13             = vec_madd(v13,v4,v0); /* fs21 */
				v31             = vec_nmsub(v30,v25,v31);

				v14             = vec_madd(v14,v5,v0); /* fs22 */
				v23             = vec_add(v23,v7);
				v15             = vec_madd(v15,v6,v0); /* fs23 */
				v16             = vec_madd(v16,v7,v0); /* fs31 */
				v1              = vec_add(v31,v1);
				v17             = vec_madd(v17,v8,v0); /* fs32 */
				v23             = vec_add(v23,v8);
				v18             = vec_madd(v18,v9,v0); /* fs33 */
				v10             = vec_madd(v10,v1,v0);

				vec_st(v24,240,(float *)stackdata); /* store Vvdwtot */
				/* calculate vectorial forces and accumulate fj. 
				 * v10-v18 has fs11-fs33 now. 
				 * First load iO-* dx,dy,dz vectors to v1-v9 
				 * and load iO forces to v28,v29,v30 
				 * use v19-v27 to accumulate j 3atoms forces 
				 */
				v28             = vec_ld(256, (float *) stackdata);
				v29             = vec_ld(272, (float *) stackdata);
				v30             = vec_ld(288, (float *) stackdata);

				v1              = vec_ld(400, (float *) stackdata);
				v2              = vec_ld(416, (float *) stackdata);
				v23             = vec_add(v23,v9); /* incr. vctot */
				v3              = vec_ld(432, (float *) stackdata);
				v4              = vec_ld(448, (float *) stackdata);
				v5              = vec_ld(464, (float *) stackdata);
				v6              = vec_ld(480, (float *) stackdata);
				vec_st(v23,224,(float *)stackdata); /* store vctot to stack */
				v7              = vec_ld(496, (float *) stackdata);
				v8              = vec_ld(512, (float *) stackdata);
				v9              = vec_ld(528, (float *) stackdata);

				v28             = vec_madd(v10,v1,v28);
				v19             = vec_nmsub(v10,v1,v0);
				v29             = vec_madd(v10,v2,v29);
				v20             = vec_nmsub(v10,v2,v0);
				v30             = vec_madd(v10,v3,v30);
				v21             = vec_nmsub(v10,v3,v0);

				v28             = vec_madd(v11,v4,v28);
				v22             = vec_nmsub(v11,v4,v0);
				v29             = vec_madd(v11,v5,v29);
				v23             = vec_nmsub(v11,v5,v0);
				v30             = vec_madd(v11,v6,v30);
				v24             = vec_nmsub(v11,v6,v0);

				v28             = vec_madd(v12,v7,v28);
				v25             = vec_nmsub(v12,v7,v0);
				v29             = vec_madd(v12,v8,v29);
				v26             = vec_nmsub(v12,v8,v0);
				v30             = vec_madd(v12,v9,v30);
				v27             = vec_nmsub(v12,v9,v0);

				/* store these i forces, and repeat for the iH1-* force */
				vec_st(v28,256,(float *)stackdata);
				vec_st(v29,272,(float *)stackdata);
				vec_st(v30,288,(float *)stackdata);

				v28             = vec_ld(304,(float *) stackdata);
				v29             = vec_ld(320,(float *) stackdata);
				v30             = vec_ld(336,(float *) stackdata);
				/* load new vectorial distances */
				v1              = vec_ld(544, (float *) stackdata);
				v2              = vec_ld(560, (float *) stackdata);
				v3              = vec_ld(576, (float *) stackdata);
				v4              = vec_ld(592, (float *) stackdata);
				v5              = vec_ld(608, (float *) stackdata);
				v6              = vec_ld(624, (float *) stackdata);
				v7              = vec_ld(640, (float *) stackdata);
				v8              = vec_ld(656, (float *) stackdata);
				v9              = vec_ld(672, (float *) stackdata);
      
				v28             = vec_madd(v13,v1,v28);
				v19             = vec_nmsub(v13,v1,v19);
				v29             = vec_madd(v13,v2,v29);
				v20             = vec_nmsub(v13,v2,v20);
				v30             = vec_madd(v13,v3,v30);
				v21             = vec_nmsub(v13,v3,v21);

				v28             = vec_madd(v14,v4,v28);
				v22             = vec_nmsub(v14,v4,v22);
				v29             = vec_madd(v14,v5,v29);
				v23             = vec_nmsub(v14,v5,v23);
				v30             = vec_madd(v14,v6,v30);
				v24             = vec_nmsub(v14,v6,v24);

				v28             = vec_madd(v15,v7,v28);
				v25             = vec_nmsub(v15,v7,v25);
				v29             = vec_madd(v15,v8,v29);
				v26             = vec_nmsub(v15,v8,v26);
				v30             = vec_madd(v15,v9,v30);
				v27             = vec_nmsub(v15,v9,v27);

				/* store these i forces, and repeat for the iH2-* force */
				vec_st(v28,304,(float *)stackdata);
				vec_st(v29,320,(float *)stackdata);
				vec_st(v30,336,(float *)stackdata);
				v28             = vec_ld(352,(float *) stackdata);
				v29             = vec_ld(368,(float *) stackdata);
				v30             = vec_ld(384,(float *) stackdata);
				/* load new vectorial distances */
				v1              = vec_ld(688, (float *) stackdata);
				v2              = vec_ld(704, (float *) stackdata);
				v3              = vec_ld(720, (float *) stackdata);
				v4              = vec_ld(736, (float *) stackdata);
				v5              = vec_ld(752, (float *) stackdata);
				v6              = vec_ld(768, (float *) stackdata);
				v7              = vec_ld(784, (float *) stackdata);
				v8              = vec_ld(800, (float *) stackdata);
				v9              = vec_ld(816, (float *) stackdata);
      
				v28             = vec_madd(v16,v1,v28);
				v19             = vec_nmsub(v16,v1,v19);
				v29             = vec_madd(v16,v2,v29);
				v20             = vec_nmsub(v16,v2,v20);
				v30             = vec_madd(v16,v3,v30);
				v21             = vec_nmsub(v16,v3,v21);

				v28             = vec_madd(v17,v4,v28);
				v22             = vec_nmsub(v17,v4,v22);
				v29             = vec_madd(v17,v5,v29);
				v23             = vec_nmsub(v17,v5,v23);
				v30             = vec_madd(v17,v6,v30);
				v24             = vec_nmsub(v17,v6,v24);

				v28             = vec_madd(v18,v7,v28);
				v25             = vec_nmsub(v18,v7,v25);
				v29             = vec_madd(v18,v8,v29);
				v26             = vec_nmsub(v18,v8,v26);
				v30             = vec_madd(v18,v9,v30);
				v27             = vec_nmsub(v18,v9,v27);

				/* store these i forces */
				vec_st(v28,352,(float *)stackdata);
				vec_st(v29,368,(float *)stackdata);
				vec_st(v30,384,(float *)stackdata);

				/* j forces present in v19-v27 */      

				v1              = vec_mergeh(v19,v21);/*  Oxa  Oza  Oxb  Ozb */
				v19             = vec_mergel(v19,v21);/*  Oxc  Ozc   -    -  */
				v21             = vec_mergeh(v20,v22);/*  Oya H1xa  Oyb H1xb */
				v20             = vec_mergel(v20,v22);/*  Oyc H1xc   -    -  */
				v22             = vec_mergeh(v23,v25);/* H1ya H2xa H1yb H2xb */
				v23             = vec_mergel(v23,v25);/* H1yc H2xc   -    -  */
				v25             = vec_mergeh(v24,v26);/* H1za H2ya H1zb H2yb */
				v24             = vec_mergel(v24,v26);/* H1zc H2yc   -    -  */

				v26             = vec_mergeh(v27,v0); /* H2za   0  H2zb   0  */
				v27             = vec_mergel(v27,v0); /* H2zc   0   -     0  */
      
				v2              = vec_mergeh(v1,v21); /*  Oxa  Oya  Oza H1xa */
				v21             = vec_mergel(v1,v21); /*  Oxb  Oyb  Ozb H1xb */
				v1              = vec_mergeh(v19,v20);/*  Oxc  Oyc  Ozc H1xc */
				v20             = vec_mergeh(v22,v25);/* H1ya H1za H2xa H2ya */
				v22             = vec_mergel(v22,v25);/* H1yb H1zb H2xb H2yb */
				v25             = vec_mergeh(v23,v24);/* H1yc H1zc H2xc H2yc */
				v24             = vec_mergeh(v26,v0); /* H2za   0    0    0  */
				v26             = vec_mergel(v26,v0); /* H2zb   0    0    0  */
				v3              = vec_mergeh(v27,v0); /* H2zc   0    0    0  */

				v29             = (vector float)vec_splat_s32(-1);

				/* move into position, load and add */  
				v30            = (vector float)vec_lvsr(0,(int *)faction+j3a); 
				v4              = vec_ld(  0, faction+j3a);
				v6              = vec_ld( 16, faction+j3a);
				v8              = vec_ld( 32, faction+j3a);
				v10             = vec_perm(v0,v29,(vector unsigned char)v30);
				v12             = vec_perm(v0,v2,(vector unsigned char)v30);
				v12             = vec_add(v12,v4);
				v14             = vec_perm(v2,v20,(vector unsigned char)v30);
				v2              = vec_add(v14,v6);
				v16             = vec_perm(v20,v24,(vector unsigned char)v30);
				v20             = vec_add(v16,v8);
				v12             = vec_sel(v4,v12,(vector unsigned int)v10);
				vec_st(v12,  0, faction+j3a);
				v10             = vec_sld(v0,v10,12);
				vec_st(v2, 16, faction+j3a);
				v12             = vec_sel(v20,v8,(vector unsigned int)v10);
				vec_st(v12, 32, faction+j3a);


				/* Finished 1, now do 2  */
				v30            = (vector float)vec_lvsr(0,(int *)faction+j3b); 
				v4              = vec_ld(  0, faction+j3b);
				v6              = vec_ld( 16, faction+j3b);
				v8              = vec_ld( 32, faction+j3b);
				v10             = vec_perm(v0,v29,(vector unsigned char)v30);
				v12             = vec_perm(v0,v21,(vector unsigned char)v30);
				v24             = vec_add(v12,v4);
				v12             = vec_perm(v21,v22,(vector unsigned char)v30);
				v21             = vec_add(v12,v6);
				v12             = vec_perm(v22,v26,(vector unsigned char)v30);
				v22             = vec_add(v12,v8);
				v12             = vec_sel(v4,v24,(vector unsigned int)v10);
				vec_st(v12,  0, faction+j3b);
				v10             = vec_sld(v0,v10,12);
				vec_st(v21, 16, faction+j3b);
				v12             = vec_sel(v22,v8,(vector unsigned int)v10);
				vec_st(v12, 32, faction+j3b);

				/* water 3 */
				v31            = (vector float)vec_lvsr(0,(int *)faction+j3c);
				v5              = vec_ld(  0, faction+j3c);
				v7              = vec_ld( 16, faction+j3c);
				v9              = vec_ld( 32, faction+j3c);
				v11             = vec_perm(v0,v29,(vector unsigned char)v31);
				v13             = vec_perm(v0,v1,(vector unsigned char)v31);
				v13             = vec_add(v13,v5);
				v15             = vec_perm(v1,v25,(vector unsigned char)v31);
				v1              = vec_add(v15,v7);      
				v17             = vec_perm(v25,v3,(vector unsigned char)v31);
				v25             = vec_add(v17,v9);
				v13             = vec_sel(v5,v13,(vector unsigned int)v11);
				vec_st(v13,  0, faction+j3c);
				v11             = vec_sld(v0,v11,12);
				vec_st(v1, 16, faction+j3c);
				v13             = vec_sel(v25,v9,(vector unsigned int)v11);
				vec_st(v13, 32, faction+j3c);

			} else if(k<(nj1-1)) {
				jnra            = jjnr[k];
				jnrb            = jjnr[k+1];
				j3a             = 3*jnra;
				j3b             = 3*jnrb;

				v1              = (vector float)vec_lvsl(0, pos+j3a);
				v8              = (vector float)vec_lvsl(0, pos+j3b);

				v2              = vec_ld(0, pos+j3a);
				v9              = vec_ld(0, pos+j3b);
				v3              = vec_ld(16, pos+j3a);
				v10             = vec_ld(16, pos+j3b);
				v4              = vec_ld(32, pos+j3a);
				v11             = vec_ld(32, pos+j3b);
				/*  Oxa  Oya  Oza H1xa */
				v5              = vec_perm(v2,v3,(vector unsigned char)v1); 
				/*  Oxb  Oyb  Ozb H1xb */
				v12             = vec_perm(v9,v10,(vector unsigned char)v8);  

				/* H1ya H1za H2xa H2ya */
				v6              = vec_perm(v3,v4,(vector unsigned char)v1); 
				/* H1yb H1zb H2xb H2yb */
				v13             = vec_perm(v10,v11,(vector unsigned char)v8); 

				/* H2za   -   -   - */    
				v7              = vec_perm(v4,v4,(vector unsigned char)v1); 
				/* H2zb   -   -   - */
				v14             = vec_perm(v11,v11,(vector unsigned char)v8); 
      
				/* permute 3atoms coordinates */
				v1              = vec_mergeh(v5,v12);/*  Oxa  Oxb  Oya  Oyb */
				v3              = vec_mergel(v5,v12);/*  Oza  Ozb H1xa H1xb */
				v5              = vec_mergeh(v6,v13);/* H1ya H1yb H1za H1zb */
				v9              = vec_mergeh(v7,v14);/* H2za H2zb   -    -  */
				v7              = vec_mergel(v6,v13);/* H2xa H2xb H2ya H2yb */

				/* load i H1x */
				v29             = vec_ld(128, (float *) stackdata); 
				v2              = vec_sld(v1,v1,8);  /*  Oya  Oyb   -   -  */
				/* load i H1y */
				v30             = vec_ld(144, (float *) stackdata); 
				v4              = vec_sld(v3,v3,8);  /* H1xa H1xb   -   -  */
				/* load i H1z */
				v31             = vec_ld(160, (float *) stackdata); 
				v6              = vec_sld(v5,v5,8);  /* H1za H1zb   -   -  */
				v8              = vec_sld(v7,v7,8);  /* H2ya H2yb   -   -  */


				v10             = vec_sub(v29,v1); /* iH1x - jOx */
				v13             = vec_sub(v29,v4); /* iH1x - jH1x */
				v16             = vec_sub(v29,v7); /* iH1x - jH2x */
				/* load i H2x */     
				v29             = vec_ld(176, (float *) stackdata); 
				v11             = vec_sub(v30,v2); /* iH1y - jOy */
				v14             = vec_sub(v30,v5); /* iH1y - jH1y */
				v17             = vec_sub(v30,v8); /* iH1y - jH2y */
				/* load i H2y */     
				v30             = vec_ld(192, (float *) stackdata); 
				vec_st(v10, 544, (float *)stackdata); /* dx21 */
				vec_st(v13, 592, (float *)stackdata); /* dx22 */
				vec_st(v16, 640, (float *)stackdata); /* dx23 */
				v12             = vec_sub(v31,v3); /* iH1z - jOz */
				v15             = vec_sub(v31,v6); /* iH1z - jH1z */
				v18             = vec_sub(v31,v9); /* iH1z - jH2z */
				/* load i H2z */         
				v31             = vec_ld(208, (float *) stackdata); 
				/* v10-v18 now contains iH1-jO, iH1-jH1 & iJ1-jH2 distances */
				vec_st(v11, 560, (float *)stackdata); /* dy21 */
				vec_st(v14, 608, (float *)stackdata); /* dy22 */
				vec_st(v17, 656, (float *)stackdata); /* dy23 */
				v19             = vec_sub(v29,v1); /* iH2x - jOx */
				v22             = vec_sub(v29,v4); /* iH2x - jH1x */
				v25             = vec_sub(v29,v7); /* iH2x - jH2x */
				vec_st(v12, 576, (float *)stackdata); /* dz21 */
				vec_st(v15, 624, (float *)stackdata); /* dz22 */
				vec_st(v18, 672, (float *)stackdata); /* dz23 */
				/* load i Ox */     
				v29             = vec_ld(80, (float *) stackdata); 
				v20             = vec_sub(v30,v2); /* iH2y - jOy */
				v23             = vec_sub(v30,v5); /* iH2y - jH1y */
				v26             = vec_sub(v30,v8); /* iH2y - jH2y */
				vec_st(v19, 688, (float *)stackdata); /* dx31 */
				vec_st(v22, 736, (float *)stackdata); /* dx32 */
				vec_st(v25, 784, (float *)stackdata); /* dx33 */
				/* load i Oy */     
				v30             = vec_ld(96, (float *) stackdata); 
				v21             = vec_sub(v31,v3); /* iH2z - jOz */
				v24             = vec_sub(v31,v6); /* iH2z - jH1z */
				v27             = vec_sub(v31,v9); /* iH2z - jH2z */
				/* load i Oz */     
				v31             = vec_ld(112, (float *) stackdata); 
				vec_st(v20, 704, (float *)stackdata); /* dy31 */
				vec_st(v23, 752, (float *)stackdata); /* dy32 */
				vec_st(v26, 800, (float *)stackdata); /* dy33 */

				v1              = vec_sub(v29,v1); /* iOx - jOx */
				v4              = vec_sub(v29,v4); /* iOx - jH1x */
				v7              = vec_sub(v29,v7); /* iOx - jH2x */
				vec_st(v21, 720, (float *)stackdata); /* dz31 */
				vec_st(v24, 768, (float *)stackdata); /* dz32 */
				vec_st(v27, 816, (float *)stackdata); /* dz33 */
				v2              = vec_sub(v30,v2); /* iOy - jOy */
				v5              = vec_sub(v30,v5); /* iOy - jH1y */
				v8              = vec_sub(v30,v8); /* iOy - jH2y */
				vec_st(v1, 400, (float *)stackdata); /* dx11 */
				vec_st(v4, 448, (float *)stackdata); /* dx12 */
				vec_st(v7, 496, (float *)stackdata); /* dx13 */
				v3              = vec_sub(v31,v3); /* iOz - jOz */
				v6              = vec_sub(v31,v6); /* iOz - jH1z */
				v9              = vec_sub(v31,v9); /* iOz - jH2z */
				vec_st(v2, 416, (float *)stackdata); /* dy11 */
				vec_st(v5, 464, (float *)stackdata); /* dy12 */
				vec_st(v8, 512, (float *)stackdata); /* dy13 */

				v1              = vec_madd(v1,v1,v0);
				v4              = vec_madd(v4,v4,v0);
				v7              = vec_madd(v7,v7,v0);
				vec_st(v3, 432, (float *)stackdata); /* dz11 */
				vec_st(v6, 480, (float *)stackdata); /* dz12 */
				vec_st(v9, 528, (float *)stackdata); /* dz13 */
				v10             = vec_madd(v10,v10,v0);
				v13             = vec_madd(v13,v13,v0);
				v16             = vec_madd(v16,v16,v0);
				v19             = vec_madd(v19,v19,v0);
				v22             = vec_madd(v22,v22,v0);
				v25             = vec_madd(v25,v25,v0);
				v1              = vec_madd(v2,v2,v1);
				v4              = vec_madd(v5,v5,v4);
				v7              = vec_madd(v8,v8,v7);
				v10             = vec_madd(v11,v11,v10);
				v13             = vec_madd(v14,v14,v13);
				v16             = vec_madd(v17,v17,v16);
				v19             = vec_madd(v20,v20,v19);
				v22             = vec_madd(v23,v23,v22);
				v25             = vec_madd(v26,v26,v25);
				v1              = vec_madd(v3,v3,v1);
				v2              = vec_madd(v6,v6,v4);
				v3              = vec_madd(v9,v9,v7);
				v4              = vec_madd(v12,v12,v10);
				v5              = vec_madd(v15,v15,v13);
				v6              = vec_madd(v18,v18,v16);
				v7              = vec_madd(v21,v21,v19);
				v8              = vec_madd(v24,v24,v22);
				v9              = vec_madd(v27,v27,v25);
				/* 
				 * v1  = rsq  iO-jO
				 * v2  = rsq  iO-jH1
				 * v3  = rsq  iO-jH2
				 * v4  = rsq  iH1-jO
				 * v5  = rsq  iH1-jH1
				 * v6  = rsq  iH1-jH2
				 * v7  = rsq  iH2-jO
				 * v8  = rsq  iH2-jH1
				 * v9 = rsq  iH2-jH2
				 */

				v10             = vec_rsqrte(v1);
				v11             = vec_rsqrte(v2);
				v12             = vec_rsqrte(v3);
				v13             = vec_rsqrte(v4);
				v14             = vec_rsqrte(v5);
				v15             = vec_rsqrte(v6);
				v16             = vec_rsqrte(v7);
				v17             = vec_rsqrte(v8);
				v18             = vec_rsqrte(v9);
				/* create constant 0.5 */
				v30             = (vector float) vec_splat_u32(1);
				/* 0.5 */
				v31             = vec_ctf((vector unsigned int)v30,1); 
				/* 1.0 */
				v30             = vec_ctf((vector unsigned int)v30,0); 

				v19             = vec_madd(v10,v10,v0); /* lu*lu */
				v20             = vec_madd(v11,v11,v0);
				v21             = vec_madd(v12,v12,v0);
				v22             = vec_madd(v13,v13,v0);
				v23             = vec_madd(v14,v14,v0);
				v24             = vec_madd(v15,v15,v0);
				v25             = vec_madd(v16,v16,v0);
				v26             = vec_madd(v17,v17,v0);
				v27             = vec_madd(v18,v18,v0);

				v19             = vec_nmsub(v1,v19,v30); /* 1.0 - rsq*lu*lu */
				v20             = vec_nmsub(v2,v20,v30);
				v21             = vec_nmsub(v3,v21,v30);
				v22             = vec_nmsub(v4,v22,v30);
				v23             = vec_nmsub(v5,v23,v30);
				v24             = vec_nmsub(v6,v24,v30);
				v25             = vec_nmsub(v7,v25,v30);
				v26             = vec_nmsub(v8,v26,v30);
				v27             = vec_nmsub(v9,v27,v30);

				v1              = vec_madd(v10,v31,v0);/* lu*0.5*/
				v2              = vec_madd(v11,v31,v0);
				v3              = vec_madd(v12,v31,v0);
				v4              = vec_madd(v13,v31,v0);
				v5              = vec_madd(v14,v31,v0);
				v6              = vec_madd(v15,v31,v0);
				v7              = vec_madd(v16,v31,v0);
				v8              = vec_madd(v17,v31,v0);
				v9              = vec_madd(v18,v31,v0);

				/* The rinv values */
				v1              = vec_madd(v1,v19,v10);
				v2              = vec_madd(v2,v20,v11);
				v3              = vec_madd(v3,v21,v12);
				v4              = vec_madd(v4,v22,v13);
				v5              = vec_madd(v5,v23,v14);
				v6              = vec_madd(v6,v24,v15);
				v7              = vec_madd(v7,v25,v16);
				v8              = vec_madd(v8,v26,v17);
				v9              = vec_madd(v9,v27,v18);
      
				v10             = (vector float)vec_splat_s32(-1);
				v10             = vec_sld(v0,v10,8);

				v1           = (vector float)vec_sel((vector unsigned int)v1,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v2           = (vector float)vec_sel((vector unsigned int)v2,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v3           = (vector float)vec_sel((vector unsigned int)v3,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v4           = (vector float)vec_sel((vector unsigned int)v4,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v5           = (vector float)vec_sel((vector unsigned int)v5,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v6           = (vector float)vec_sel((vector unsigned int)v6,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v7           = (vector float)vec_sel((vector unsigned int)v7,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v8           = (vector float)vec_sel((vector unsigned int)v8,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				v9           = (vector float)vec_sel((vector unsigned int)v9,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);

				/* load qqOO, qqOH and qqHH  to v27,v28,v29 */
				v27             = vec_ld(0, (float *) stackdata);
				v28             = vec_ld(16, (float *) stackdata);
				v29             = vec_ld(32, (float *) stackdata);

				/*				vec_dstst( faction+j3a, 0x10010100, 2 ); */
     
				/* put rinvsq in v10-v18, rinv6_OO in v30 & rinv12_OO in v31 */
				/* load c6 to v25 and c12 to v26 */
				v25             = vec_ld(48, (float *) stackdata);
				v26             = vec_ld(64, (float *) stackdata);
      
				v10             = vec_madd(v1,v1,v0);
				v1              = vec_madd(v1,v27,v0); /* rinv11*qqOO */
				v11             = vec_madd(v2,v2,v0);
				/* load vctot to v23 and Vvdwtot to v24 */
				v23             = vec_ld(224,(float *) stackdata);
				v24             = vec_ld(240,(float *) stackdata);

				v2              = vec_madd(v2,v28,v0); /* rinv12*qqOH */
				v12             = vec_madd(v3,v3,v0);
				v30             = vec_madd(v10,v10,v0); /* rinv4 */
				v3              = vec_madd(v3,v28,v0); /* rinv13*qqOH */
				v13             = vec_madd(v4,v4,v0);
				v4              = vec_madd(v4,v28,v0); /* rinv21*qqOH */
				v14             = vec_madd(v5,v5,v0);

				v23             = vec_add(v23,v1);

				v30             = vec_madd(v30,v10,v0); /* rinv6 */
				v5              = vec_madd(v5,v29,v0); /* rinv22*qqHH */
				v15             = vec_madd(v6,v6,v0);
				v6              = vec_madd(v6,v29,v0); /* rinv23*qqHH */
				v23             = vec_add(v23,v2);
				v16             = vec_madd(v7,v7,v0);
				v31             = vec_madd(v30,v30,v0); /* rinv12 */
				v25             = vec_madd(v25,v30,v0); /* c6*rinv6 */
				/* load 6.0 to v30 */
				v30             = (vector float)vec_splat_u32(6);
				v30             = vec_ctf((vector unsigned int)v30,0);
				v23             = vec_add(v23,v3);

				v7              = vec_madd(v7,v28,v0); /* rinv31*qqOH */
				v17             = vec_madd(v8,v8,v0);
				v8              = vec_madd(v8,v29,v0); /* rinv32*qqHH */
				v26             = vec_madd(v26,v31,v0); /* c12*rinv12 */
				v23             = vec_add(v23,v4);
				/* load 12.0 to v31 */
				v31             = (vector float)vec_splat_u32(12);
				v31             = vec_ctf((vector unsigned int)v31,0);

				v24             = vec_sub(v24,v25);  /* add Vvdw6 to Vvdwtot */
				v18             = vec_madd(v9,v9,v0);
				v23             = vec_add(v23,v5);
				v9              = vec_madd(v9,v29,v0); /* rinv33*qqHH */

				v24             = vec_add(v24,v26);/* add Vvdw12 to Vvdwtot */
    
				v31             = vec_madd(v31,v26,v0);
				v11             = vec_madd(v11,v2,v0); /* fs12 */
				v23             = vec_add(v23,v6);
				v12             = vec_madd(v12,v3,v0); /* fs13 */
				v13             = vec_madd(v13,v4,v0); /* fs21 */
				v31             = vec_nmsub(v30,v25,v31);

				v14             = vec_madd(v14,v5,v0); /* fs22 */
				v23             = vec_add(v23,v7);
				v15             = vec_madd(v15,v6,v0); /* fs23 */
				v16             = vec_madd(v16,v7,v0); /* fs31 */
				v1              = vec_add(v31,v1);
				v17             = vec_madd(v17,v8,v0); /* fs32 */
				v23             = vec_add(v23,v8);
				v18             = vec_madd(v18,v9,v0); /* fs33 */
				v10             = vec_madd(v10,v1,v0);

				vec_st(v24,240,(float *)stackdata); /* store Vvdwtot */
				/* calculate vectorial forces and accumulate fj. 
				 * v10-v18 has fs11-fs33 now. 
				 * First load iO-* dx,dy,dz vectors to v1-v9 
				 * and load iO forces to v28,v29,v30 
				 * use v19-v27 to accumulate j 3atoms forces 
				 */
				v28             = vec_ld(256, (float *) stackdata);
				v29             = vec_ld(272, (float *) stackdata);
				v30             = vec_ld(288, (float *) stackdata);

				v1              = vec_ld(400, (float *) stackdata);
				v2              = vec_ld(416, (float *) stackdata);
				v23             = vec_add(v23,v9); /* incr. vctot */
				v3              = vec_ld(432, (float *) stackdata);
				v4              = vec_ld(448, (float *) stackdata);
				v5              = vec_ld(464, (float *) stackdata);
				v6              = vec_ld(480, (float *) stackdata);
				/* store vctot back to stack */
				vec_st(v23,224,(float *)stackdata); 
				v7              = vec_ld(496, (float *) stackdata);
				v8              = vec_ld(512, (float *) stackdata);
				v9              = vec_ld(528, (float *) stackdata);

				v28             = vec_madd(v10,v1,v28);
				v19             = vec_nmsub(v10,v1,v0);
				v29             = vec_madd(v10,v2,v29);
				v20             = vec_nmsub(v10,v2,v0);
				v30             = vec_madd(v10,v3,v30);
				v21             = vec_nmsub(v10,v3,v0);

				v28             = vec_madd(v11,v4,v28);
				v22             = vec_nmsub(v11,v4,v0);
				v29             = vec_madd(v11,v5,v29);
				v23             = vec_nmsub(v11,v5,v0);
				v30             = vec_madd(v11,v6,v30);
				v24             = vec_nmsub(v11,v6,v0);

				v28             = vec_madd(v12,v7,v28);
				v25             = vec_nmsub(v12,v7,v0);
				v29             = vec_madd(v12,v8,v29);
				v26             = vec_nmsub(v12,v8,v0);
				v30             = vec_madd(v12,v9,v30);
				v27             = vec_nmsub(v12,v9,v0);

				/* store these i forces, and repeat for the iH1-* force */
				vec_st(v28,256,(float *)stackdata);
				vec_st(v29,272,(float *)stackdata);
				vec_st(v30,288,(float *)stackdata);

				v28             = vec_ld(304,(float *) stackdata);
				v29             = vec_ld(320,(float *) stackdata);
				v30             = vec_ld(336,(float *) stackdata);
				/* load new vectorial distances */
				v1              = vec_ld(544, (float *) stackdata);
				v2              = vec_ld(560, (float *) stackdata);
				v3              = vec_ld(576, (float *) stackdata);
				v4              = vec_ld(592, (float *) stackdata);
				v5              = vec_ld(608, (float *) stackdata);
				v6              = vec_ld(624, (float *) stackdata);
				v7              = vec_ld(640, (float *) stackdata);
				v8              = vec_ld(656, (float *) stackdata);
				v9              = vec_ld(672, (float *) stackdata);
      
				v28             = vec_madd(v13,v1,v28);
				v19             = vec_nmsub(v13,v1,v19);
				v29             = vec_madd(v13,v2,v29);
				v20             = vec_nmsub(v13,v2,v20);
				v30             = vec_madd(v13,v3,v30);
				v21             = vec_nmsub(v13,v3,v21);

				v28             = vec_madd(v14,v4,v28);
				v22             = vec_nmsub(v14,v4,v22);
				v29             = vec_madd(v14,v5,v29);
				v23             = vec_nmsub(v14,v5,v23);
				v30             = vec_madd(v14,v6,v30);
				v24             = vec_nmsub(v14,v6,v24);

				v28             = vec_madd(v15,v7,v28);
				v25             = vec_nmsub(v15,v7,v25);
				v29             = vec_madd(v15,v8,v29);
				v26             = vec_nmsub(v15,v8,v26);
				v30             = vec_madd(v15,v9,v30);
				v27             = vec_nmsub(v15,v9,v27);

				/* store these i forces, and repeat for the iH2-* force */
				vec_st(v28,304,(float *)stackdata);
				vec_st(v29,320,(float *)stackdata);
				vec_st(v30,336,(float *)stackdata);
				v28             = vec_ld(352,(float *) stackdata);
				v29             = vec_ld(368,(float *) stackdata);
				v30             = vec_ld(384,(float *) stackdata);
				/* load new vectorial distances */
				v1              = vec_ld(688, (float *) stackdata);
				v2              = vec_ld(704, (float *) stackdata);
				v3              = vec_ld(720, (float *) stackdata);
				v4              = vec_ld(736, (float *) stackdata);
				v5              = vec_ld(752, (float *) stackdata);
				v6              = vec_ld(768, (float *) stackdata);
				v7              = vec_ld(784, (float *) stackdata);
				v8              = vec_ld(800, (float *) stackdata);
				v9              = vec_ld(816, (float *) stackdata);
      
				v28             = vec_madd(v16,v1,v28);
				v19             = vec_nmsub(v16,v1,v19);
				v29             = vec_madd(v16,v2,v29);
				v20             = vec_nmsub(v16,v2,v20);
				v30             = vec_madd(v16,v3,v30);
				v21             = vec_nmsub(v16,v3,v21);

				v28             = vec_madd(v17,v4,v28);
				v22             = vec_nmsub(v17,v4,v22);
				v29             = vec_madd(v17,v5,v29);
				v23             = vec_nmsub(v17,v5,v23);
				v30             = vec_madd(v17,v6,v30);
				v24             = vec_nmsub(v17,v6,v24);

				v28             = vec_madd(v18,v7,v28);
				v25             = vec_nmsub(v18,v7,v25);
				v29             = vec_madd(v18,v8,v29);
				v26             = vec_nmsub(v18,v8,v26);
				v30             = vec_madd(v18,v9,v30);
				v27             = vec_nmsub(v18,v9,v27);

				/* store these i forces */
				vec_st(v28,352,(float *)stackdata);
				vec_st(v29,368,(float *)stackdata);
				vec_st(v30,384,(float *)stackdata);

				/* j forces present in v19-v27 */      

				/*  Oxa  Oza  Oxb  Ozb */
				v1              = vec_mergeh(v19,v21); 
				/*  Oya H1xa  Oyb H1xb */
				v21             = vec_mergeh(v20,v22); 
				/* H1ya H2xa H1yb H2xb */
				v22             = vec_mergeh(v23,v25); 
				/* H1za H2ya H1zb H2yb */
				v25             = vec_mergeh(v24,v26); 

				/* H2za   0  H2zb   0  */
				v26             = vec_mergeh(v27,v0);   
      
				/*  Oxa  Oya  Oza H1xa */
				v2              = vec_mergeh(v1,v21);   
				/*  Oxb  Oyb  Ozb H1xb */
				v21             = vec_mergel(v1,v21);   
				/* H1ya H1za H2xa H2ya */
				v20             = vec_mergeh(v22,v25);  
				/* H1yb H1zb H2xb H2yb */
				v22             = vec_mergel(v22,v25);  
				/* H2za   0    0    0  */
				v24             = vec_mergeh(v26,v0);   
				/* H2zb   0    0    0  */
				v26             = vec_mergel(v26,v0);   

				v29             = (vector float)vec_splat_s32(-1);

				/* move into position, load and add */  
				v30            = (vector float)vec_lvsr(0,(int *)faction+j3a); 
				v4              = vec_ld(  0, faction+j3a);
				v6              = vec_ld( 16, faction+j3a);
				v8              = vec_ld( 32, faction+j3a);
				v10             = vec_perm(v0,v29,(vector unsigned char)v30);
				v12             = vec_perm(v0,v2,(vector unsigned char)v30);
				v12             = vec_add(v12,v4);
				v14             = vec_perm(v2,v20,(vector unsigned char)v30);
				v2              = vec_add(v14,v6);
				v16             = vec_perm(v20,v24,(vector unsigned char)v30);
				v20             = vec_add(v16,v8);
				v12             = vec_sel(v4,v12,(vector unsigned int)v10);
				vec_st(v12,  0, faction+j3a);
				v10             = vec_sld(v0,v10,12);
				vec_st(v2, 16, faction+j3a);
				v12             = vec_sel(v20,v8,(vector unsigned int)v10);
				vec_st(v12, 32, faction+j3a);


				/* Finished 1, now do 2  */
				v30            = (vector float)vec_lvsr(0,(int *)faction+j3b); 
				v4              = vec_ld(  0, faction+j3b);
				v6              = vec_ld( 16, faction+j3b);
				v8              = vec_ld( 32, faction+j3b);
				v10             = vec_perm(v0,v29,(vector unsigned char)v30);
				v12             = vec_perm(v0,v21,(vector unsigned char)v30);
				v24             = vec_add(v12,v4);
				v12             = vec_perm(v21,v22,(vector unsigned char)v30);
				v21             = vec_add(v12,v6);
				v12             = vec_perm(v22,v26,(vector unsigned char)v30);
				v22             = vec_add(v12,v8);
				v12             = vec_sel(v4,v24,(vector unsigned int)v10);
				vec_st(v12,  0, faction+j3b);
				v10             = vec_sld(v0,v10,12);
				vec_st(v21, 16, faction+j3b);
				v12             = vec_sel(v22,v8,(vector unsigned int)v10);
				vec_st(v12, 32, faction+j3b);

			} else if(k<nj1) {
				jnra            = jjnr[k];
				j3a             = 3*jnra;

				v10             = (vector float)vec_lvsl(0, pos+j3a);

				v2              = vec_ld(0, pos+j3a);
				v3              = vec_ld(16, pos+j3a);
				v4              = vec_ld(32, pos+j3a);
				/*  Oxa  Oya  Oza H1xa */
				v1              = vec_perm(v2,v3,(vector unsigned char)v10); 
				/* H1ya H1za H2xa H2ya */
				v5              = vec_perm(v3,v4,(vector unsigned char)v10); 
				/* H2za   -   -   - */  
				v9              = vec_perm(v4,v4,(vector unsigned char)v10); 

				/* permute 3atoms coordinates */
				/* just splat things... never mind filling all cells :-) */
				/* load i H1x */
				v29             = vec_ld(128, (float *) stackdata); 
				v2              = vec_splat(v1,1);
				/* load i H1y */
				v30             = vec_ld(144, (float *) stackdata); 
				v3              = vec_splat(v1,2);
				/* load i H1z */
				v31             = vec_ld(160, (float *) stackdata); 
				v4              = vec_splat(v1,3);
				v6              = vec_splat(v5,1);
				v7              = vec_splat(v5,2);
				v8              = vec_splat(v5,3);

				v10             = vec_sub(v29,v1); /* iH1x - jOx */
				v13             = vec_sub(v29,v4); /* iH1x - jH1x */
				v16             = vec_sub(v29,v7); /* iH1x - jH2x */
				/* load i H2x */     
				v29             = vec_ld(176, (float *) stackdata); 
				v11             = vec_sub(v30,v2); /* iH1y - jOy */
				v14             = vec_sub(v30,v5); /* iH1y - jH1y */
				v17             = vec_sub(v30,v8); /* iH1y - jH2y */
				/* load i H2y */     
				v30             = vec_ld(192, (float *) stackdata); 
				vec_st(v10, 544, (float *)stackdata); /* dx21 */
				vec_st(v13, 592, (float *)stackdata); /* dx22 */
				vec_st(v16, 640, (float *)stackdata); /* dx23 */
				v12             = vec_sub(v31,v3); /* iH1z - jOz */
				v15             = vec_sub(v31,v6); /* iH1z - jH1z */
				v18             = vec_sub(v31,v9); /* iH1z - jH2z */
				/* load i H2z */         
				v31             = vec_ld(208, (float *) stackdata); 
				/* v10-v18 now contains iH1-jO, iH1-jH1 and iJ1-jH2 dist. */
				vec_st(v11, 560, (float *)stackdata); /* dy21 */
				vec_st(v14, 608, (float *)stackdata); /* dy22 */
				vec_st(v17, 656, (float *)stackdata); /* dy23 */
				v19             = vec_sub(v29,v1); /* iH2x - jOx */
				v22             = vec_sub(v29,v4); /* iH2x - jH1x */
				v25             = vec_sub(v29,v7); /* iH2x - jH2x */
				vec_st(v12, 576, (float *)stackdata); /* dz21 */
				vec_st(v15, 624, (float *)stackdata); /* dz22 */
				vec_st(v18, 672, (float *)stackdata); /* dz23 */
				/* load i Ox */     
				v29             = vec_ld(80, (float *) stackdata); 
				v20             = vec_sub(v30,v2); /* iH2y - jOy */
				v23             = vec_sub(v30,v5); /* iH2y - jH1y */
				v26             = vec_sub(v30,v8); /* iH2y - jH2y */
				vec_st(v19, 688, (float *)stackdata); /* dx31 */
				vec_st(v22, 736, (float *)stackdata); /* dx32 */
				vec_st(v25, 784, (float *)stackdata); /* dx33 */
				/* load i Oy */     
				v30             = vec_ld(96, (float *) stackdata); 
				v21             = vec_sub(v31,v3); /* iH2z - jOz */
				v24             = vec_sub(v31,v6); /* iH2z - jH1z */
				v27             = vec_sub(v31,v9); /* iH2z - jH2z */
				/* load i Oz */     
				v31             = vec_ld(112, (float *) stackdata); 
				vec_st(v20, 704, (float *)stackdata); /* dy31 */
				vec_st(v23, 752, (float *)stackdata); /* dy32 */
				vec_st(v26, 800, (float *)stackdata); /* dy33 */

				v1              = vec_sub(v29,v1); /* iOx - jOx */
				v4              = vec_sub(v29,v4); /* iOx - jH1x */
				v7              = vec_sub(v29,v7); /* iOx - jH2x */
				vec_st(v21, 720, (float *)stackdata); /* dz31 */
				vec_st(v24, 768, (float *)stackdata); /* dz32 */
				vec_st(v27, 816, (float *)stackdata); /* dz33 */
				v2              = vec_sub(v30,v2); /* iOy - jOy */
				v5              = vec_sub(v30,v5); /* iOy - jH1y */
				v8              = vec_sub(v30,v8); /* iOy - jH2y */
				vec_st(v1, 400, (float *)stackdata); /* dx11 */
				vec_st(v4, 448, (float *)stackdata); /* dx12 */
				vec_st(v7, 496, (float *)stackdata); /* dx13 */
				v3              = vec_sub(v31,v3); /* iOz - jOz */
				v6              = vec_sub(v31,v6); /* iOz - jH1z */
				v9              = vec_sub(v31,v9); /* iOz - jH2z */
				vec_st(v2, 416, (float *)stackdata); /* dy11 */
				vec_st(v5, 464, (float *)stackdata); /* dy12 */
				vec_st(v8, 512, (float *)stackdata); /* dy13 */

				v1              = vec_madd(v1,v1,v0);
				v4              = vec_madd(v4,v4,v0);
				v7              = vec_madd(v7,v7,v0);
				vec_st(v3, 432, (float *)stackdata); /* dz11 */
				vec_st(v6, 480, (float *)stackdata); /* dz12 */
				vec_st(v9, 528, (float *)stackdata); /* dz13 */
				v10             = vec_madd(v10,v10,v0);
				v13             = vec_madd(v13,v13,v0);
				v16             = vec_madd(v16,v16,v0);
				v19             = vec_madd(v19,v19,v0);
				v22             = vec_madd(v22,v22,v0);
				v25             = vec_madd(v25,v25,v0);
				v1              = vec_madd(v2,v2,v1);
				v4              = vec_madd(v5,v5,v4);
				v7              = vec_madd(v8,v8,v7);
				v10             = vec_madd(v11,v11,v10);
				v13             = vec_madd(v14,v14,v13);
				v16             = vec_madd(v17,v17,v16);
				v19             = vec_madd(v20,v20,v19);
				v22             = vec_madd(v23,v23,v22);
				v25             = vec_madd(v26,v26,v25);
				v1              = vec_madd(v3,v3,v1);
				v2              = vec_madd(v6,v6,v4);
				v3              = vec_madd(v9,v9,v7);
				v4              = vec_madd(v12,v12,v10);
				v5              = vec_madd(v15,v15,v13);
				v6              = vec_madd(v18,v18,v16);
				v7              = vec_madd(v21,v21,v19);
				v8              = vec_madd(v24,v24,v22);
				v9              = vec_madd(v27,v27,v25);
				/* 
				 * v1  = rsq  iO-jO
				 * v2  = rsq  iO-jH1
				 * v3  = rsq  iO-jH2
				 * v4  = rsq  iH1-jO
				 * v5  = rsq  iH1-jH1
				 * v6  = rsq  iH1-jH2
				 * v7  = rsq  iH2-jO
				 * v8  = rsq  iH2-jH1
				 * v9 = rsq  iH2-jH2
				 */

				v10             = vec_rsqrte(v1);
				v11             = vec_rsqrte(v2);
				v12             = vec_rsqrte(v3);
				v13             = vec_rsqrte(v4);
				v14             = vec_rsqrte(v5);
				v15             = vec_rsqrte(v6);
				v16             = vec_rsqrte(v7);
				v17             = vec_rsqrte(v8);
				v18             = vec_rsqrte(v9);
				/* create constant 0.5 */
				v30             = (vector float) vec_splat_u32(1);
				/* 0.5 */
				v31             = vec_ctf((vector unsigned int)v30,1); 
				/* 1.0 */
				v30             = vec_ctf((vector unsigned int)v30,0); 

				v19             = vec_madd(v10,v10,v0); /* lu*lu */
				v20             = vec_madd(v11,v11,v0);
				v21             = vec_madd(v12,v12,v0);
				v22             = vec_madd(v13,v13,v0);
				v23             = vec_madd(v14,v14,v0);
				v24             = vec_madd(v15,v15,v0);
				v25             = vec_madd(v16,v16,v0);
				v26             = vec_madd(v17,v17,v0);
				v27             = vec_madd(v18,v18,v0);

				v19             = vec_nmsub(v1,v19,v30); /* 1.0 - rsq*lu*lu */
				v20             = vec_nmsub(v2,v20,v30);
				v21             = vec_nmsub(v3,v21,v30);
				v22             = vec_nmsub(v4,v22,v30);
				v23             = vec_nmsub(v5,v23,v30);
				v24             = vec_nmsub(v6,v24,v30);
				v25             = vec_nmsub(v7,v25,v30);
				v26             = vec_nmsub(v8,v26,v30);
				v27             = vec_nmsub(v9,v27,v30);

				v1              = vec_madd(v10,v31,v0);/* lu*0.5*/
				v2              = vec_madd(v11,v31,v0);
				v3              = vec_madd(v12,v31,v0);
				v4              = vec_madd(v13,v31,v0);
				v5              = vec_madd(v14,v31,v0);
				v6              = vec_madd(v15,v31,v0);
				v7              = vec_madd(v16,v31,v0);
				v8              = vec_madd(v17,v31,v0);
				v9              = vec_madd(v18,v31,v0);

				/* The rinv values */
				v1              = vec_madd(v1,v19,v10);
				v2              = vec_madd(v2,v20,v11);
				v3              = vec_madd(v3,v21,v12);
				v4              = vec_madd(v4,v22,v13);
				v5              = vec_madd(v5,v23,v14);
				v6              = vec_madd(v6,v24,v15);
				v7              = vec_madd(v7,v25,v16);
				v8              = vec_madd(v8,v26,v17);
				v9              = vec_madd(v9,v27,v18);
      
				v10             = (vector float)vec_splat_s32(-1);
				v10             = vec_sld(v0,v10,12);

				v1           = (vector float)vec_sel((vector unsigned int)v1,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);

				v2           = (vector float)vec_sel((vector unsigned int)v2,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				
				v3           = (vector float)vec_sel((vector unsigned int)v3,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				
				v4           = (vector float)vec_sel((vector unsigned int)v4,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				
				v5           = (vector float)vec_sel((vector unsigned int)v5,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				
				v6           = (vector float)vec_sel((vector unsigned int)v6,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);
				
				v7           = (vector float)vec_sel((vector unsigned int)v7,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);

				v8           = (vector float)vec_sel((vector unsigned int)v8,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);

				v9           = (vector float)vec_sel((vector unsigned int)v9,
													 (vector unsigned int)v0,
													 (vector unsigned int)v10);

				/* load qqOO, qqOH and qqHH  to v27,v28,v29 */
				v27             = vec_ld(0, (float *) stackdata);
				v28             = vec_ld(16, (float *) stackdata);
				v29             = vec_ld(32, (float *) stackdata);
				/*				vec_dstst( faction+j3a, 0x10010100, 2 ); */
      
				/* put rinvsq in v10-v18, rinv6_OO in v30 and
				 * rinv12_OO in v31
				 */
				/* load c6 to v25 and c12 to v26 */
				v25             = vec_ld(48, (float *) stackdata);
				v26             = vec_ld(64, (float *) stackdata);
      
				v10             = vec_madd(v1,v1,v0);
				v1              = vec_madd(v1,v27,v0); /* rinv11*qqOO */
				v11             = vec_madd(v2,v2,v0);
				/* load vctot to v23 and Vvdwtot to v24 */
				v23             = vec_ld(224,(float *) stackdata);
				v24             = vec_ld(240,(float *) stackdata);

				v2              = vec_madd(v2,v28,v0); /* rinv12*qqOH */
				v12             = vec_madd(v3,v3,v0);
				v30             = vec_madd(v10,v10,v0); /* rinv4 */
				v3              = vec_madd(v3,v28,v0); /* rinv13*qqOH */
				v13             = vec_madd(v4,v4,v0);
				v4              = vec_madd(v4,v28,v0); /* rinv21*qqOH */
				v14             = vec_madd(v5,v5,v0);

				v23             = vec_add(v23,v1);

				v30             = vec_madd(v30,v10,v0); /* rinv6 */
				v5              = vec_madd(v5,v29,v0); /* rinv22*qqHH */
				v15             = vec_madd(v6,v6,v0);
				v6              = vec_madd(v6,v29,v0); /* rinv23*qqHH */
				v23             = vec_add(v23,v2);
				v16             = vec_madd(v7,v7,v0);
				v31             = vec_madd(v30,v30,v0); /* rinv12 */
				v25             = vec_madd(v25,v30,v0); /* c6*rinv6 */
				/* load 6.0 to v30 */
				v30             = (vector float)vec_splat_u32(6);
				v30             = vec_ctf((vector unsigned int)v30,0);
				v23             = vec_add(v23,v3);

				v7              = vec_madd(v7,v28,v0); /* rinv31*qqOH */
				v17             = vec_madd(v8,v8,v0);
				v8              = vec_madd(v8,v29,v0); /* rinv32*qqHH */
				v26             = vec_madd(v26,v31,v0); /* c12*rinv12 */
				v23             = vec_add(v23,v4);
				/* load 12.0 to v31 */
				v31             = (vector float)vec_splat_u32(12);
				v31             = vec_ctf((vector unsigned int)v31,0);

				v24             = vec_sub(v24,v25);  /* add Vvdw6 to Vvdwtot */
				v18             = vec_madd(v9,v9,v0);
				v23             = vec_add(v23,v5);
				v9              = vec_madd(v9,v29,v0); /* rinv33*qqHH */

				v24             = vec_add(v24,v26);/* add Vvdw12 to Vvdwtot */
    
				v31             = vec_madd(v31,v26,v0);
				v11             = vec_madd(v11,v2,v0); /* fs12 */
				v23             = vec_add(v23,v6);
				v12             = vec_madd(v12,v3,v0); /* fs13 */
				v13             = vec_madd(v13,v4,v0); /* fs21 */
				v31             = vec_nmsub(v30,v25,v31);

				v14             = vec_madd(v14,v5,v0); /* fs22 */
				v23             = vec_add(v23,v7);
				v15             = vec_madd(v15,v6,v0); /* fs23 */
				v16             = vec_madd(v16,v7,v0); /* fs31 */
				v1              = vec_add(v31,v1);
				v17             = vec_madd(v17,v8,v0); /* fs32 */
				v23             = vec_add(v23,v8);
				v18             = vec_madd(v18,v9,v0); /* fs33 */
				v10             = vec_madd(v10,v1,v0);

				vec_st(v24,240,(float *)stackdata); /* store Vvdwtot */
				/* calculate vectorial forces and accumulate fj. 
                 * v10-v18 has fs11-fs33 now. 
				 * First load iO-* dx,dy,dz vectors to v1-v9 
				 * and load iO forces to v28,v29,v30 
				 * use v19-v27 to accumulate j 3atoms forces 
				 */
				v28             = vec_ld(256, (float *) stackdata);
				v29             = vec_ld(272, (float *) stackdata);
				v30             = vec_ld(288, (float *) stackdata);

				v1              = vec_ld(400, (float *) stackdata);
				v2              = vec_ld(416, (float *) stackdata);
				v23             = vec_add(v23,v9); /* incr. vctot */
				v3              = vec_ld(432, (float *) stackdata);
				v4              = vec_ld(448, (float *) stackdata);
				v5              = vec_ld(464, (float *) stackdata);
				v6              = vec_ld(480, (float *) stackdata);
				vec_st(v23,224,(float *)stackdata); /* store vctot to stack */
				v7              = vec_ld(496, (float *) stackdata);
				v8              = vec_ld(512, (float *) stackdata);
				v9              = vec_ld(528, (float *) stackdata);

				v28             = vec_madd(v10,v1,v28);
				v19             = vec_nmsub(v10,v1,v0);
				v29             = vec_madd(v10,v2,v29);
				v20             = vec_nmsub(v10,v2,v0);
				v30             = vec_madd(v10,v3,v30);
				v21             = vec_nmsub(v10,v3,v0);

				v28             = vec_madd(v11,v4,v28);
				v22             = vec_nmsub(v11,v4,v0);
				v29             = vec_madd(v11,v5,v29);
				v23             = vec_nmsub(v11,v5,v0);
				v30             = vec_madd(v11,v6,v30);
				v24             = vec_nmsub(v11,v6,v0);

				v28             = vec_madd(v12,v7,v28);
				v25             = vec_nmsub(v12,v7,v0);
				v29             = vec_madd(v12,v8,v29);
				v26             = vec_nmsub(v12,v8,v0);
				v30             = vec_madd(v12,v9,v30);
				v27             = vec_nmsub(v12,v9,v0);

				/* store these i forces, and repeat the procedure 
				 * for the iH1-* force 
				 */
				vec_st(v28,256,(float *)stackdata);
				vec_st(v29,272,(float *)stackdata);
				vec_st(v30,288,(float *)stackdata);

				v28             = vec_ld(304,(float *) stackdata);
				v29             = vec_ld(320,(float *) stackdata);
				v30             = vec_ld(336,(float *) stackdata);
				/* load new vectorial distances */
				v1              = vec_ld(544, (float *) stackdata);
				v2              = vec_ld(560, (float *) stackdata);
				v3              = vec_ld(576, (float *) stackdata);
				v4              = vec_ld(592, (float *) stackdata);
				v5              = vec_ld(608, (float *) stackdata);
				v6              = vec_ld(624, (float *) stackdata);
				v7              = vec_ld(640, (float *) stackdata);
				v8              = vec_ld(656, (float *) stackdata);
				v9              = vec_ld(672, (float *) stackdata);
      
				v28             = vec_madd(v13,v1,v28);
				v19             = vec_nmsub(v13,v1,v19);
				v29             = vec_madd(v13,v2,v29);
				v20             = vec_nmsub(v13,v2,v20);
				v30             = vec_madd(v13,v3,v30);
				v21             = vec_nmsub(v13,v3,v21);

				v28             = vec_madd(v14,v4,v28);
				v22             = vec_nmsub(v14,v4,v22);
				v29             = vec_madd(v14,v5,v29);
				v23             = vec_nmsub(v14,v5,v23);
				v30             = vec_madd(v14,v6,v30);
				v24             = vec_nmsub(v14,v6,v24);

				v28             = vec_madd(v15,v7,v28);
				v25             = vec_nmsub(v15,v7,v25);
				v29             = vec_madd(v15,v8,v29);
				v26             = vec_nmsub(v15,v8,v26);
				v30             = vec_madd(v15,v9,v30);
				v27             = vec_nmsub(v15,v9,v27);

				/* store these i forces, and repeat the procedure
				 * for the iH2-* force
				 */
				vec_st(v28,304,(float *)stackdata);
				vec_st(v29,320,(float *)stackdata);
				vec_st(v30,336,(float *)stackdata);
				v28             = vec_ld(352,(float *) stackdata);
				v29             = vec_ld(368,(float *) stackdata);
				v30             = vec_ld(384,(float *) stackdata);
				/* load new vectorial distances */
				v1              = vec_ld(688, (float *) stackdata);
				v2              = vec_ld(704, (float *) stackdata);
				v3              = vec_ld(720, (float *) stackdata);
				v4              = vec_ld(736, (float *) stackdata);
				v5              = vec_ld(752, (float *) stackdata);
				v6              = vec_ld(768, (float *) stackdata);
				v7              = vec_ld(784, (float *) stackdata);
				v8              = vec_ld(800, (float *) stackdata);
				v9              = vec_ld(816, (float *) stackdata);
      
				v28             = vec_madd(v16,v1,v28);
				v19             = vec_nmsub(v16,v1,v19);
				v29             = vec_madd(v16,v2,v29);
				v20             = vec_nmsub(v16,v2,v20);
				v30             = vec_madd(v16,v3,v30);
				v21             = vec_nmsub(v16,v3,v21);

				v28             = vec_madd(v17,v4,v28);
				v22             = vec_nmsub(v17,v4,v22);
				v29             = vec_madd(v17,v5,v29);
				v23             = vec_nmsub(v17,v5,v23);
				v30             = vec_madd(v17,v6,v30);
				v24             = vec_nmsub(v17,v6,v24);

				v28             = vec_madd(v18,v7,v28);
				v25             = vec_nmsub(v18,v7,v25);
				v29             = vec_madd(v18,v8,v29);
				v26             = vec_nmsub(v18,v8,v26);
				v30             = vec_madd(v18,v9,v30);
				v27             = vec_nmsub(v18,v9,v27);

				/* store these i forces */
				vec_st(v28,352,(float *)stackdata);
				vec_st(v29,368,(float *)stackdata);
				vec_st(v30,384,(float *)stackdata);

				/* j forces present in v19-v27 */      

				/*  Oxa  Oza   -    -  */
				v1              = vec_mergeh(v19,v21); 
				/*  Oya H1xa   -    -  */
				v21             = vec_mergeh(v20,v22); 
				/* H1ya H2xa   -    -  */
				v22             = vec_mergeh(v23,v25); 
				/* H1za H2ya   -    -  */
				v25             = vec_mergeh(v24,v26); 

				/* H2za   0   -    0  */
				v26             = vec_mergeh(v27,v0);   
      
				/*  Oxa  Oya  Oza H1xa */
				v2              = vec_mergeh(v1,v21);   
				/* H1ya H1za H2xa H2ya */
				v20             = vec_mergeh(v22,v25);  
				/* H2za   0    0    0  */
				v24             = vec_mergeh(v26,v0);   
 

				v29             = (vector float)vec_splat_s32(-1);

				/* move into position, load and add */  
				v30            = (vector float)vec_lvsr(0,(int *)faction+j3a); 
				v4              = vec_ld(  0, faction+j3a);
				v6              = vec_ld( 16, faction+j3a);
				v8              = vec_ld( 32, faction+j3a);
				v10             = vec_perm(v0,v29,(vector unsigned char)v30);
				v12             = vec_perm(v0,v2,(vector unsigned char)v30);
				v12             = vec_add(v12,v4);
				v14             = vec_perm(v2,v20,(vector unsigned char)v30);
				v2              = vec_add(v14,v6);
				v16             = vec_perm(v20,v24,(vector unsigned char)v30);
				v20             = vec_add(v16,v8);
				v12             = vec_sel(v4,v12,(vector unsigned int)v10);
				vec_st(v12,  0, faction+j3a);
				v10             = vec_sld(v0,v10,12);
				vec_st(v2, 16, faction+j3a);
				v12             = vec_sel(v20,v8,(vector unsigned int)v10);
				vec_st(v12, 32, faction+j3a);

			}

			v1          = (vector float)vec_lvsr(0,faction+ii3); 
			v5          = (vector float)vec_splat_s32(-1);
			v2          = vec_ld( 0, faction+ii3);
			v3          = vec_ld(16, faction+ii3);
			v4          = vec_ld(32, faction+ii3);
			v5          = vec_perm(v0, v5,(vector unsigned char)v1); /* mask */
			/* load forces from stack */
			v6          = vec_ld(256, (float *) stackdata); /* Ox */
			v7          = vec_ld(272, (float *) stackdata); /* Oy */
			v8          = vec_ld(288, (float *) stackdata); /* Oz */
			v9          = vec_ld(304, (float *) stackdata); /* H1x */
			v10         = vec_ld(320, (float *) stackdata); /* H1y */
			v11         = vec_ld(336, (float *) stackdata); /* H1z */
			v12         = vec_ld(352, (float *) stackdata); /* H2x */
			v13         = vec_ld(368, (float *) stackdata); /* H2y */
			v14         = vec_ld(384, (float *) stackdata); /* H2z */

			/* accumulate the forces */
			v15         = vec_sld(v6,v6,8);
			v16         = vec_sld(v7,v7,8);
			v17         = vec_sld(v8,v8,8);
			v18         = vec_sld(v9,v9,8);
			v19         = vec_sld(v10,v10,8);
			v20         = vec_sld(v11,v11,8);
			v21         = vec_sld(v12,v12,8);
			v22         = vec_sld(v13,v13,8);
			v23         = vec_sld(v14,v14,8);

			v6          = vec_add(v6,v15);  /*  Ox  Ox' - - */
			v7          = vec_add(v7,v16);  /*  Oy  Oy' - - */
			v8          = vec_add(v8,v17);  /*  Oz  Oz' - - */
			v9          = vec_add(v9,v18);  /* H1x H1x' - - */
			v10         = vec_add(v10,v19); /* H1y H1y' - - */
			v11         = vec_add(v11,v20); /* H1z H1z' - - */
			v12         = vec_add(v12,v21); /* H2x H2x' - - */
			v13         = vec_add(v13,v22); /* H2y H2y' - - */
			v14         = vec_add(v14,v23); /* H2z H2z' - - */
    
			v6          = vec_mergeh(v6,v8);   /*  Ox  Oz  Ox'  Oz' */
			v7          = vec_mergeh(v7,v9);   /*  Oy H1x  Oy' H1x' */
			v10         = vec_mergeh(v10,v12); /* H1y H2x H1y' H2x' */
			v11         = vec_mergeh(v11,v13); /* H1z H2y H1z' H2y' */
			v14         = vec_mergeh(v14,v0);  /* H2z  0  H2z'  0   */

			v15         = vec_sld(v6,v6,8);
			v16         = vec_sld(v7,v7,8);
			v17         = vec_sld(v10,v10,8);
			v18         = vec_sld(v11,v11,8);
			v19         = vec_sld(v14,v14,8);

			v6          = vec_add(v6,v15); /* Ox Oz - - */
			v7          = vec_add(v7,v16); /* Oy H1x - - */
			v10         = vec_add(v10,v17);/* H1y H2x - - */
			v11         = vec_add(v11,v18);/* H1z H2y - - */
			v14         = vec_add(v14,v19);/* H2z 0 - 0 */
    
			v6          = vec_mergeh(v6,v7);   /*  Ox  Oy  Oz H1x */
			v10         = vec_mergeh(v10,v11); /* H1y H1z H2x H2y */
			v14         = vec_mergeh(v14,v0);  /* H2z  0   0   0  */

			v7          = vec_sld(v0,v6,12);   /* 0   Ox  Oy  Oz  */
			v8          = vec_sld(v6,v10,8);   /* -  H1x H1y H1z  */
			v9          = vec_sld(v10,v14,4);  /* -  H2x H2y H2z  */

			/* The part to add to v2 */
			v12         = vec_perm(v0,v6,(vector unsigned char)v1);   
			/* The part to add to v3 */
			v13         = vec_perm(v6,v10,(vector unsigned char)v1);  
			/* The part to add to v4 */
			v14         = vec_perm(v10,v14,(vector unsigned char)v1); 

			v12         = vec_add(v2,v12);
			v13         = vec_add(v3,v13);
			v14         = vec_add(v4,v14);

			v12         = vec_sel(v2,v12,(vector unsigned int)v5);
			v5          = vec_sld(v0,v5,12);
			v14         = vec_sel(v14,v4,(vector unsigned int)v5);

			/* store */
			vec_st(v12, 0, faction+ii3); 
			vec_st(v13,16, faction+ii3);
			vec_st(v14,32, faction+ii3);

			/* accumulate for shift */
			v7          = vec_add(v7,v8);
			v7          = vec_add(v7,v9);
			v7          = vec_sld(v7,v0,4); /* x y z 0 */

			/* add v7 to the memory location fshift+is3 */
			v15         = vec_lde(0, fshift+is3);
			v16         = vec_lde(4, fshift+is3);
			v17         = vec_lde(8, fshift+is3);
			v18         = (vector float)vec_splat(v7,0);
			v19         = (vector float)vec_splat(v7,1);
			v20         = (vector float)vec_splat(v7,2);
			v15         = vec_add(v15,v18);
			v16         = vec_add(v16,v19);
			v17         = vec_add(v17,v20);
			vec_ste(v15,0,fshift+is3);
			vec_ste(v16,4,fshift+is3);
			vec_ste(v17,8,fshift+is3);

			/* update potential energies */
			v1          = vec_ld(224,(float *) stackdata); /* load vctot */
			v2          = vec_ld(240,(float *) stackdata); /* load Vvdwtot */
			v3          = vec_sld(v1,v1,8);
			v4          = vec_sld(v2,v2,8);
			v1          = vec_add(v1,v3);
			v2          = vec_add(v2,v4);
			v3          = vec_sld(v1,v1,4);
			v4          = vec_sld(v2,v2,4);
			v1          = vec_add(v1,v3);
			v2          = vec_add(v2,v4);
			/* all 4 positions in v1, v2 contain the sum now */
			v3          = vec_lde(0, Vc+gid[n]);
			v4          = vec_lde(0, Vvdw+gid[n]);
			v3          = vec_add(v1,v3);
			v4          = vec_add(v2,v4);
			vec_ste(v3,0,Vc+gid[n]);
			vec_ste(v4,0,Vvdw+gid[n]);
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
nb_kernel112nf_ppc_altivec(int *             p_nri,
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
	vector float rinv11,rinv12,rinv13,rinv21,rinv22,rinv23;
	vector float rinv31,rinv32,rinv33;
	vector float rinvsq11;

	vector float vfacel,nul;
	vector float vctot,qqOO,qqOH,qqHH,qO,qH,c6,c12,rinvsix;
	vector float Vvdwtot,qqOOt,qqOHt,qqHHt,c6t,c12t;

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
  
				rinvsq11        = vec_madd(rinv11,rinv11,nul); 
				rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
				Vvdwtot          = vec_nmsub(c6,rinvsix,Vvdwtot);
				rinvsix         = vec_madd(rinvsix,rinvsix,nul);
				vctot           = vec_madd(rinv11,qqOO,vctot);
				vctot           = vec_madd(rinv12,qqOH,vctot);
				vctot           = vec_madd(rinv13,qqOH,vctot);
				Vvdwtot          = vec_madd(c12,rinvsix,Vvdwtot);
				vctot           = vec_madd(rinv21,qqOH,vctot);
				vctot           = vec_madd(rinv22,qqHH,vctot);
				vctot           = vec_madd(rinv23,qqHH,vctot);
				vctot           = vec_madd(rinv31,qqOH,vctot);
				vctot           = vec_madd(rinv32,qqHH,vctot);
				vctot           = vec_madd(rinv33,qqHH,vctot);
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

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				zero_highest_element_in_9_vectors(&rinv11,&rinv12,&rinv13,
												  &rinv21,&rinv22,&rinv23,
												  &rinv31,&rinv32,&rinv33);

				rinvsq11        = vec_madd(rinv11,rinv11,nul); 
				rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
				Vvdwtot          = vec_nmsub(c6t,rinvsix,Vvdwtot);
				rinvsix         = vec_madd(rinvsix,rinvsix,nul);
				vctot           = vec_madd(rinv11,qqOOt,vctot);
				vctot           = vec_madd(rinv12,qqOHt,vctot);
				vctot           = vec_madd(rinv13,qqOHt,vctot);
				Vvdwtot          = vec_madd(c12t,rinvsix,Vvdwtot);
				vctot           = vec_madd(rinv21,qqOHt,vctot);
				vctot           = vec_madd(rinv22,qqHHt,vctot);
				vctot           = vec_madd(rinv23,qqHHt,vctot);
				vctot           = vec_madd(rinv31,qqOHt,vctot);
				vctot           = vec_madd(rinv32,qqHHt,vctot);
				vctot           = vec_madd(rinv33,qqHHt,vctot);
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

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				zero_highest_2_elements_in_9_vectors(&rinv11,&rinv12,&rinv13,
													 &rinv21,&rinv22,&rinv23,
													 &rinv31,&rinv32,&rinv33);

				rinvsq11        = vec_madd(rinv11,rinv11,nul); 
				rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
				Vvdwtot          = vec_nmsub(c6t,rinvsix,Vvdwtot);
				rinvsix         = vec_madd(rinvsix,rinvsix,nul);
				vctot           = vec_madd(rinv11,qqOOt,vctot);
				vctot           = vec_madd(rinv12,qqOHt,vctot);
				vctot           = vec_madd(rinv13,qqOHt,vctot);
				Vvdwtot          = vec_madd(c12t,rinvsix,Vvdwtot);
				vctot           = vec_madd(rinv21,qqOHt,vctot);
				vctot           = vec_madd(rinv22,qqHHt,vctot);
				vctot           = vec_madd(rinv23,qqHHt,vctot);
				vctot           = vec_madd(rinv31,qqOHt,vctot);
				vctot           = vec_madd(rinv32,qqHHt,vctot);
				vctot           = vec_madd(rinv33,qqHHt,vctot);
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

				do_9_invsqrt(rsq11,rsq12,rsq13,
							 rsq21,rsq22,rsq23,
							 rsq31,rsq32,rsq33,
							 &rinv11,&rinv12,&rinv13,
							 &rinv21,&rinv22,&rinv23,
							 &rinv31,&rinv32,&rinv33);
      
				zero_highest_3_elements_in_9_vectors(&rinv11,&rinv12,&rinv13,
													 &rinv21,&rinv22,&rinv23,
													 &rinv31,&rinv32,&rinv33);

				rinvsq11        = vec_madd(rinv11,rinv11,nul); 
				rinvsix         = vec_madd(rinvsq11,rinvsq11,nul);
				rinvsix         = vec_madd(rinvsix,rinvsq11,nul);
				Vvdwtot          = vec_nmsub(c6t,rinvsix,Vvdwtot);
				rinvsix         = vec_madd(rinvsix,rinvsix,nul);
				vctot           = vec_madd(rinv11,qqOOt,vctot);
				vctot           = vec_madd(rinv12,qqOHt,vctot);
				vctot           = vec_madd(rinv13,qqOHt,vctot);
				Vvdwtot          = vec_madd(c12t,rinvsix,Vvdwtot);
				vctot           = vec_madd(rinv21,qqOHt,vctot);
				vctot           = vec_madd(rinv22,qqHHt,vctot);
				vctot           = vec_madd(rinv23,qqHHt,vctot);
				vctot           = vec_madd(rinv31,qqOHt,vctot);
				vctot           = vec_madd(rinv32,qqHHt,vctot);
				vctot           = vec_madd(rinv33,qqHHt,vctot);
			}
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
