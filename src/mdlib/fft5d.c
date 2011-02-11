/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef NOGMX
#define GMX_PARALLEL_ENV_INITIALIZED 1
#else 
#include "main.h"
#define GMX_PARALLEL_ENV_INITIALIZED gmx_parallel_env_initialized()
#endif

#ifdef GMX_LIB_MPI
#include <mpi.h>
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#ifdef FFT5D_THREADS
#include <omp.h>
/* requires fftw compiled with openmp */
#define FFT5D_FFTW_THREADS
#endif

#include "fft5d.h"
#include <float.h>
#include <math.h>
#include <assert.h>

#ifndef __FLT_EPSILON__
#define __FLT_EPSILON__ FLT_EPSILON
#define __DBL_EPSILON__ DBL_EPSILON
#endif

#ifdef NOGMX
FILE* debug=0;
#endif

#include "gmx_fatal.h"


#ifdef GMX_FFT_FFTW3 
#ifdef GMX_THREADS
/* none of the fftw3 calls, except execute(), are thread-safe, so 
   we need to serialize them with this mutex. */
static tMPI_Thread_mutex_t big_fftw_mutex=TMPI_THREAD_MUTEX_INITIALIZER;

#define FFTW_LOCK tMPI_Thread_mutex_lock(&big_fftw_mutex);
#define FFTW_UNLOCK tMPI_Thread_mutex_unlock(&big_fftw_mutex);
#else /* GMX_THREADS */
#define FFTW_LOCK 
#define FFTW_UNLOCK 
#endif /* GMX_THREADS */
#endif /* GMX_FFT_FFTW3 */

static double fft5d_fmax(double a, double b){
	return (a>b)?a:b;
}

/* largest factor smaller than sqrt */
static int lfactor(int z) {  
	int i;
	for (i=sqrt(z);;i--)
		if (z%i==0) return i;
}

/* largest factor */
static int l2factor(int z) {  
	int i;
	if (z==1) return 1;
	for (i=z/2;;i--)
		if (z%i==0) return i;
}

/* largest prime factor: WARNING: slow recursion, only use for small numbers */
static int lpfactor(int z) {
	int f = l2factor(z);
	if (f==1) return z;
	return fft5d_fmax(lpfactor(f),lpfactor(z/f));
}

#ifndef GMX_MPI
#ifdef HAVE_GETTIMEOFDAY
#include <sys/time.h>
double MPI_Wtime() {
    struct timeval tv;
    gettimeofday(&tv,0);
    return tv.tv_sec+tv.tv_usec*1e-6;
}
#else
double MPI_Wtime() {
    return 0.0;
}
#endif
#endif

static int vmax(int* a, int s) {
    int i,max=0;
    for (i=0;i<s;i++) 
    {
        if (a[i]>max) max=a[i];
    }
    return max;
} 

/*
copied here from fftgrid, because:
1. function there not publically available
2. not sure whether we keep fftgrid
3. less dependencies for fft5d

Only used for non-fftw case
*/
static void *
gmx_calloc_aligned(size_t size)
{
    void *p0,*p;
    
    /*We initialize by zero for Valgrind
      For non-divisible case we communicate more than the data.
      If we don't initialize the data we communicate uninitialized data*/
    p0 = calloc(size+32,1);  
    
    if(p0 == NULL)
    {
        gmx_fatal(FARGS,"Failed to allocated %u bytes of aligned memory.",size+32);
    }
    
    p = (void *) (((size_t) p0 + 32) & (~((size_t) 31)));
    
    /* Yeah, yeah, we cannot free this pointer, but who cares... */
    return p;
}


/* NxMxK the size of the data
 * comm communicator to use for fft5d
 * P0 number of processor in 1st axes (can be null for automatic)
 * lin is allocated by fft5d because size of array is only known after planning phase */
fft5d_plan fft5d_plan_3d(int NG, int MG, int KG, MPI_Comm comm[2], int flags, t_complex** rlin, t_complex** rlout)
{

    int P[2],bMaster,prank[2],i;
    int rNG,rMG,rKG;
    int *N0=0, *N1=0, *M0=0, *M1=0, *K0=0, *K1=0, *oN0=0, *oN1=0, *oM0=0, *oM1=0, *oK0=0, *oK1=0;
    int N[3],M[3],K[3],pN[3],pM[3],pK[3],oM[3],oK[3],*iNin[3]={0},*oNin[3]={0},*iNout[3]={0},*oNout[3]={0};
    int C[3],rC[3],nP[2];
    int lsize;
    t_complex *lin=0,*lout=0;
    fft5d_plan plan;
    int s;

    /* comm, prank and P are in the order of the decomposition (plan->cart is in the order of transposes) */
#ifdef GMX_MPI
    if (GMX_PARALLEL_ENV_INITIALIZED && comm[0] != 0)
    {
        MPI_Comm_size(comm[0],&P[0]);
        MPI_Comm_rank(comm[0],&prank[0]);
    }
    else
#endif
    {
        P[0] = 1;
        prank[0] = 0;
    }
#ifdef GMX_MPI
    if (GMX_PARALLEL_ENV_INITIALIZED && comm[1] != 0)
    {
        MPI_Comm_size(comm[1],&P[1]);
        MPI_Comm_rank(comm[1],&prank[1]);
    }
    else
#endif
    {
        P[1] = 1;
        prank[1] = 0;
    }
   
    bMaster=(prank[0]==0&&prank[1]==0);
   
    
    if (debug)
    {
        fprintf(debug,"FFT5D: Using %dx%d processor grid, rank %d,%d\n",
                P[0],P[1],prank[0],prank[1]);
    }
    
    if (bMaster) {
        if (debug) 
            fprintf(debug,"FFT5D: N: %d, M: %d, K: %d, P: %dx%d, real2complex: %d, backward: %d, order yz: %d, debug %d\n",
                NG,MG,KG,P[0],P[1],(flags&FFT5D_REALCOMPLEX)>0,(flags&FFT5D_BACKWARD)>0,(flags&FFT5D_ORDER_YZ)>0,(flags&FFT5D_DEBUG)>0);
        /* The check below is not correct, one prime factor 11 or 13 is ok.
        if (fft5d_fmax(fft5d_fmax(lpfactor(NG),lpfactor(MG)),lpfactor(KG))>7) {
            printf("WARNING: FFT very slow with prime factors larger 7\n");
            printf("Change FFT size or in case you cannot change it look at\n");
            printf("http://www.fftw.org/fftw3_doc/Generating-your-own-code.html\n");
        }
        */
    }
    
    if (NG==0 || MG==0 || KG==0) {
        if (bMaster) printf("FFT5D: FATAL: Datasize cannot be zero in any dimension\n");
        return 0;
    }

    rNG=NG;rMG=MG;rKG=KG;
    
    if (flags&FFT5D_REALCOMPLEX) {
        if (!(flags&FFT5D_BACKWARD)) NG = NG/2+1;
        else {
            if (!(flags&FFT5D_ORDER_YZ)) MG=MG/2+1;
            else KG=KG/2+1;
        }
    }
    
    
    /*for transpose we need to know the size for each processor not only our own size*/

    N0 = (int*)malloc(P[0]*sizeof(int)); N1 = (int*)malloc(P[1]*sizeof(int)); 
    M0 = (int*)malloc(P[0]*sizeof(int)); M1 = (int*)malloc(P[1]*sizeof(int));
    K0 = (int*)malloc(P[0]*sizeof(int)); K1 = (int*)malloc(P[1]*sizeof(int));
    oN0 = (int*)malloc(P[0]*sizeof(int));oN1 = (int*)malloc(P[1]*sizeof(int));
    oM0 = (int*)malloc(P[0]*sizeof(int));oM1 = (int*)malloc(P[1]*sizeof(int));
    oK0 = (int*)malloc(P[0]*sizeof(int));oK1 = (int*)malloc(P[1]*sizeof(int));
    
    for (i=0;i<P[0];i++) 
    {
        #define EVENDIST
        #ifndef EVENDIST
        oN0[i]=i*ceil((double)NG/P[0]);
        oM0[i]=i*ceil((double)MG/P[0]);
        oK0[i]=i*ceil((double)KG/P[0]);
        #else
        oN0[i]=(NG*i)/P[0];
        oM0[i]=(MG*i)/P[0];
        oK0[i]=(KG*i)/P[0];
        #endif
    }
    for (i=0;i<P[1];i++) 
    {
        #ifndef EVENDIST
        oN1[i]=i*ceil((double)NG/P[1]); 
        oM1[i]=i*ceil((double)MG/P[1]); 
        oK1[i]=i*ceil((double)KG/P[1]); 
        #else
        oN1[i]=(NG*i)/P[1]; 
        oM1[i]=(MG*i)/P[1]; 
        oK1[i]=(KG*i)/P[1]; 
        #endif
    }
    for (i=0;i<P[0]-1;i++) 
    {
        N0[i]=oN0[i+1]-oN0[i];
        M0[i]=oM0[i+1]-oM0[i];
        K0[i]=oK0[i+1]-oK0[i];
    }
    N0[P[0]-1]=NG-oN0[P[0]-1];
    M0[P[0]-1]=MG-oM0[P[0]-1];
    K0[P[0]-1]=KG-oK0[P[0]-1];
    for (i=0;i<P[1]-1;i++) 
    {
        N1[i]=oN1[i+1]-oN1[i];
        M1[i]=oM1[i+1]-oM1[i];
        K1[i]=oK1[i+1]-oK1[i];
    }
    N1[P[1]-1]=NG-oN1[P[1]-1];
    M1[P[1]-1]=MG-oM1[P[1]-1];
    K1[P[1]-1]=KG-oK1[P[1]-1];

    /* for step 1-3 the local N,M,K sizes of the transposed system
       C: contiguous dimension, and nP: number of processor in subcommunicator 
       for that step */
    
    
    pM[0] = M0[prank[0]];
    oM[0] = oM0[prank[0]];
    pK[0] = K1[prank[1]];
    oK[0] = oK1[prank[1]];
    C[0] = NG;
    rC[0] = rNG;
    if (!(flags&FFT5D_ORDER_YZ)) {
        N[0] = vmax(N1,P[1]);
        M[0] = M0[prank[0]];
        K[0] = vmax(K1,P[1]);
        pN[0] = N1[prank[1]];
        iNout[0] = N1;
        oNout[0] = oN1;
        nP[0] = P[1];
        C[1] = KG;
        rC[1] =rKG;
        N[1] = vmax(K0,P[0]);
        pN[1] = K0[prank[0]];
        iNin[1] = K1;
        oNin[1] = oK1; 
        iNout[1] = K0;
        oNout[1] = oK0;
        M[1] = vmax(M0,P[0]);
        pM[1] = M0[prank[0]];
        oM[1] = oM0[prank[0]];
        K[1] = N1[prank[1]];
        pK[1] = N1[prank[1]];
        oK[1] = oN1[prank[1]];
        nP[1] = P[0];
        C[2] = MG;
        rC[2] = rMG;
        iNin[2] = M0;
        oNin[2] = oM0;
        M[2] = vmax(K0,P[0]);
        pM[2] = K0[prank[0]];
        oM[2] = oK0[prank[0]];
        K[2] = vmax(N1,P[1]);
        pK[2] = N1[prank[1]];
        oK[2] = oN1[prank[1]];
        free(N0); free(oN0); /*these are not used for this order*/
        free(M1); free(oM1); /*the rest is freed in destroy*/
    } else {
        N[0] = vmax(N0,P[0]);
        M[0] = vmax(M0,P[0]);
        K[0] = K1[prank[1]];
        pN[0] = N0[prank[0]];
        iNout[0] = N0;
        oNout[0] = oN0;
        nP[0] = P[0];
        C[1] = MG;
        rC[1] =rMG;
        N[1] = vmax(M1,P[1]);
        pN[1] = M1[prank[1]];
        iNin[1] = M0;
        oNin[1] = oM0;
        iNout[1] = M1;
        oNout[1] = oM1;
        M[1] = N0[prank[0]];
        pM[1] = N0[prank[0]];
        oM[1] = oN0[prank[0]];
        K[1] = vmax(K1,P[1]);
        pK[1] = K1[prank[1]];
        oK[1] = oK1[prank[1]];
        nP[1] = P[1];
        C[2] = KG;
        rC[2] = rKG;
        iNin[2] = K1;
        oNin[2] = oK1;
        M[2] = vmax(N0,P[0]);
        pM[2] = N0[prank[0]];
        oM[2] = oN0[prank[0]];
        K[2] = vmax(M1,P[1]);
        pK[2] = M1[prank[1]];
        oK[2] = oM1[prank[1]];
        free(N1); free(oN1); /*these are not used for this order*/
        free(K0); free(oK0); /*the rest is freed in destroy*/
    }
    
    /*
      Difference between x-y-z regarding 2d decomposition is whether they are 
      distributed along axis 1, 2 or both 
    */
    
    /* int lsize = fmax(N[0]*M[0]*K[0]*nP[0],N[1]*M[1]*K[1]*nP[1]); */
    lsize = fft5d_fmax(N[0]*M[0]*K[0]*nP[0],fft5d_fmax(N[1]*M[1]*K[1]*nP[1],C[2]*M[2]*K[2])); 
    /* int lsize = fmax(C[0]*M[0]*K[0],fmax(C[1]*M[1]*K[1],C[2]*M[2]*K[2])); */
    if (!(flags&FFT5D_NOMALLOC)) { 
        lin = (t_complex*)gmx_calloc_aligned(sizeof(t_complex) * lsize);   
        lout = (t_complex*)gmx_calloc_aligned(sizeof(t_complex) * lsize); 
    } else {
        lin = *rlin;
        lout = *rlout;
    }

    plan = (fft5d_plan)calloc(1,sizeof(struct fft5d_plan_t));

    
#ifdef FFT5D_THREADS
#ifdef FFT5D_FFTW_THREADS
    FFTW(init_threads)();
    int nthreads;
    #pragma omp parallel
    {
        #pragma omp master
        {
            nthreads = omp_get_num_threads();
        }
    }
    if (prank[0] == 0 && prank[1] == 0)
    {
        printf("Running fftw on %d threads\n",nthreads);        
    }
    FFTW(plan_with_nthreads)(nthreads);
#endif
#endif    

#ifdef GMX_FFT_FFTW3  /*if not FFTW - then we don't do a 3d plan but insead only 1D plans */
    if ((!(flags&FFT5D_INPLACE)) && (!(P[0]>1 || P[1]>1))) {  /*don't do 3d plan in parallel or if in_place requested */  
            int fftwflags=FFTW_DESTROY_INPUT;
            fftw_iodim dims[3];
            int inNG=NG,outMG=MG,outKG=KG;

            FFTW_LOCK;
            if (!(flags&FFT5D_NOMEASURE)) fftwflags|=FFTW_MEASURE;
            if (flags&FFT5D_REALCOMPLEX) {
                if (!(flags&FFT5D_BACKWARD)) {  /*input pointer is not complex*/
                    inNG*=2; 
                } else {                        /*output pointer is not complex*/
                    if (!(flags&FFT5D_ORDER_YZ)) outMG*=2;
                    else outKG*=2;
                }
            }

            if (!(flags&FFT5D_BACKWARD)) {
                dims[0].n  = KG;
                dims[1].n  = MG;
                dims[2].n  = rNG;
                
                dims[0].is = inNG*MG;     /*N M K*/
                dims[1].is = inNG;
                dims[2].is = 1;
                if (!(flags&FFT5D_ORDER_YZ)) {
                    dims[0].os = MG;       /*M K N*/
                    dims[1].os = 1;
                    dims[2].os = MG*KG;
                } else  {
                    dims[0].os = 1;       /*K N M*/
                    dims[1].os = KG*NG;
                    dims[2].os = KG;
                }
            } else {
                if (!(flags&FFT5D_ORDER_YZ)) {
                    dims[0].n  = NG;   
                    dims[1].n  = KG;   
                    dims[2].n  = rMG;  
                    
                    dims[0].is = 1;     
                    dims[1].is = NG*MG;
                    dims[2].is = NG;

                    dims[0].os = outMG*KG;       
                    dims[1].os = outMG;
                    dims[2].os = 1;                  
                } else {
                    dims[0].n  = MG;
                    dims[1].n  = NG;
                    dims[2].n  = rKG;
                    
                    dims[0].is = NG;     
                    dims[1].is = 1;
                    dims[2].is = NG*MG;

                    dims[0].os = outKG*NG;       
                    dims[1].os = outKG;
                    dims[2].os = 1;                  
                }           
            }
            if ((flags&FFT5D_REALCOMPLEX) && !(flags&FFT5D_BACKWARD)) {
                plan->p3d = FFTW(plan_guru_dft_r2c)(/*rank*/ 3, dims,
                                     /*howmany*/ 0, /*howmany_dims*/0 ,
                                     (real*)lin, (FFTW(complex) *)lout,
                                     /*flags*/ fftwflags);              
            } else if ((flags&FFT5D_REALCOMPLEX) && (flags&FFT5D_BACKWARD)) {
                plan->p3d = FFTW(plan_guru_dft_c2r)(/*rank*/ 3, dims,
                                     /*howmany*/ 0, /*howmany_dims*/0 ,
                                     (FFTW(complex) *)lin, (real*)lout,
                                     /*flags*/ fftwflags);              
            } else {
                plan->p3d = FFTW(plan_guru_dft)(/*rank*/ 3, dims,
                                     /*howmany*/ 0, /*howmany_dims*/0 ,
                                     (FFTW(complex) *)lin, (FFTW(complex) *)lout,
                                     /*sign*/ (flags&FFT5D_BACKWARD)?1:-1, /*flags*/ fftwflags);
            }
            FFTW_UNLOCK;
    }
    if (!plan->p3d) {  /* for decomposition and if 3d plan did not work */
#endif /* GMX_FFT_FFTW3 */
        for (s=0;s<3;s++) {
            if (debug)
            {
                fprintf(debug,"FFT5D: Plan s %d rC %d M %d pK %d C %d lsize %d\n",
                        s,rC[s],M[s],pK[s],C[s],lsize);
            }
            if ((flags&FFT5D_REALCOMPLEX) && ((!(flags&FFT5D_BACKWARD) && s==0) || ((flags&FFT5D_BACKWARD) && s==2))) {
                gmx_fft_init_many_1d_real( &plan->p1d[s], rC[s], pM[s]*pK[s], (flags&FFT5D_NOMEASURE)?GMX_FFT_FLAG_CONSERVATIVE:0 );
            } else {
                gmx_fft_init_many_1d     ( &plan->p1d[s],  C[s], pM[s]*pK[s], (flags&FFT5D_NOMEASURE)?GMX_FFT_FLAG_CONSERVATIVE:0 );
            }
        }
#ifdef GMX_FFT_FFTW3 
    }
#endif
    if ((flags&FFT5D_ORDER_YZ)) { /*plan->cart is in the order of transposes */
        plan->cart[0]=comm[0]; plan->cart[1]=comm[1];
    } else {
        plan->cart[1]=comm[0]; plan->cart[0]=comm[1];
    }
#ifdef FFT5D_MPI_TRANSPOSE
    FFTW_LOCK
    for (s=0;s<2;s++) {
        if ((s==0 && !(flags&FFT5D_ORDER_YZ)) || (s==1 && (flags&FFT5D_ORDER_YZ))) 
            plan->mpip[s] = FFTW(mpi_plan_many_transpose)(nP[s], nP[s], N[s]*K[s]*pM[s]*2, 1, 1, (real*)lin, (real*)lout, plan->cart[s], FFTW_PATIENT);
        else
            plan->mpip[s] = FFTW(mpi_plan_many_transpose)(nP[s], nP[s], N[s]*pK[s]*M[s]*2, 1, 1, (real*)lin, (real*)lout, plan->cart[s], FFTW_PATIENT);
    }
    FFTW_UNLOCK
#endif 

    
    plan->lin=lin;
    plan->lout=lout;
    
    plan->NG=NG;plan->MG=MG;plan->KG=KG;
    
    for (s=0;s<3;s++) {
        plan->N[s]=N[s];plan->M[s]=M[s];plan->K[s]=K[s];plan->pN[s]=pN[s];plan->pM[s]=pM[s];plan->pK[s]=pK[s];
        plan->oM[s]=oM[s];plan->oK[s]=oK[s];
        plan->C[s]=C[s];plan->rC[s]=rC[s];
        plan->iNin[s]=iNin[s];plan->oNin[s]=oNin[s];plan->iNout[s]=iNout[s];plan->oNout[s]=oNout[s];
    }
    for (s=0;s<2;s++) {
        plan->P[s]=nP[s];plan->coor[s]=prank[s];
    }
    
/*    plan->fftorder=fftorder;
    plan->direction=direction;    
    plan->realcomplex=realcomplex;
*/
    plan->flags=flags;
    *rlin=lin;
    *rlout=lout;
    return plan;
}


enum order {
    XYZ,
    XZY,
    YXZ,
    YZX,
    ZXY,
    ZYX
};



/*here x,y,z and N,M,K is in rotated coordinate system!!
  x (and N) is mayor (consecutive) dimension, y (M) middle and z (K) major
  maxN,maxM,maxK is max size of local data
  pN, pM, pK is local size specific to current processor (only different to max if not divisible)
  NG, MG, KG is size of global data*/
static void splitaxes(t_complex* lout,const t_complex* lin,
                      int maxN,int maxM,int maxK, int pN, int pM, int pK,
                      int P,int NG,int *N, int* oN)
{
    int x,y,z,i;
    int in_i,out_i,in_z,out_z,in_y,out_y;

#ifdef FFT5D_THREADS
    int zi;

    /* In the thread parallel case we want to loop over z and i
     * in a single for loop to allow for better load balancing.
     */
#pragma omp parallel for private(z,in_z,out_z,i,in_i,out_i,y,in_y,out_y,x) schedule(static)
    for (zi=0; zi<pK*P; zi++)
    {
        z = zi/P;
        i = zi - z*P;
#else
    for (z=0; z<pK; z++) /*3. z l*/ 
    {
#endif
        in_z  = z*maxN*maxM;
        out_z = z*NG*pM;

#ifndef FFT5D_THREADS
        for (i=0; i<P; i++) /*index cube along long axis*/
#endif
        {
            in_i  = in_z  + i*maxN*maxM*maxK;
            out_i = out_z + oN[i];
            for (y=0;y<pM;y++) { /*2. y k*/
                in_y  = in_i  + y*maxN;
                out_y = out_i + y*NG;
                for (x=0;x<N[i];x++) { /*1. x j*/
                    lout[in_y+x] = lin[out_y+x];
                    /*after split important that each processor chunk i has size maxN*maxM*maxK and thus being the same size*/
                    /*before split data contiguos - thus if different processor get different amount oN is different*/
                }
            }
        }
    }
}

/*make axis contiguous again (after AllToAll) and also do local transpose*/
/*transpose mayor and major dimension
  variables see above
  the major, middle, minor order is only correct for x,y,z (N,M,K) for the input
  N,M,K local dimensions
  KG global size*/
static void joinAxesTrans13(t_complex* lin,const t_complex* lout,
                            int maxN,int maxM,int maxK,int pN, int pM, int pK, 
                            int P,int KG, int* K, int* oK)
{
    int i,x,y,z;
    int in_i,out_i,in_x,out_x,in_z,out_z;

#ifdef FFT5D_THREADS
    int xi;

    /* In the thread parallel case we want to loop over x and i
     * in a single for loop to allow for better load balancing.
     */
#pragma omp parallel for private(x,in_x,out_x,i,in_i,out_i,z,in_z,out_z,y) schedule(static)
    for (xi=0; xi<pN*P; xi++)
    {
        x = xi/P;
        i = xi - x*P;
#else
    for (x=0;x<pN;x++) /*1.j*/
    {
#endif
        in_x  = x*KG*pM;
        out_x = x;

#ifndef FFT5D_THREADS
        for (i=0;i<P;i++) /*index cube along long axis*/
#endif
        {
            in_i  = in_x  + oK[i];
            out_i = out_x + i*maxM*maxN*maxK;
            for (z=0;z<K[i];z++) /*3.l*/
            {
                in_z  = in_i  + z;
                out_z = out_i + z*maxM*maxN;
                for (y=0;y<pM;y++) { /*2.k*/
                    lin[in_z+y*KG] = lout[out_z+y*maxN];
                }
            }
        }
    }
}

/*make axis contiguous again (after AllToAll) and also do local transpose
  tranpose mayor and middle dimension
  variables see above
  the minor, middle, major order is only correct for x,y,z (N,M,K) for the input
  N,M,K local size
  MG, global size*/
static void joinAxesTrans12(t_complex* lin,const t_complex* lout,int maxN,int maxM,int maxK,int pN, int pM, int pK,
                int P,int MG, int* M, int* oM) {
    int i,z,y,x;
    int in_i,out_i,in_z,out_z,in_x,out_x;

#ifdef FFT5D_THREADS
    int zi;

    /* In the thread parallel case we want to loop over z and i
     * in a single for loop to allow for better load balancing.
     */
#pragma omp parallel for private(i,in_i,out_i,z,in_z,out_z,in_x,out_x,x,y) schedule(static)
    for (zi=0; zi<pK*P; zi++)
    {
        z = zi/P;
        i = zi - z*P;
#else
    for (z=0; z<pK; z++)
    {
#endif
        in_z  = z*MG*pN;
        out_z = z*maxM*maxN;

#ifndef FFT5D_THREADS
        for (i=0; i<P; i++) /*index cube along long axis*/
#endif
        {
            in_i  = in_z  + oM[i];
            out_i = out_z + i*maxM*maxN*maxK;
            for (x=0;x<pN;x++) {
                in_x  = in_i  + x*MG;
                out_x = out_i + x;
                for (y=0;y<M[i];y++) {
                    lin[in_x+y] = lout[out_x+y*maxN];
                }
            }
        }
    }
}


static void rotate(int x[]) {
    int t=x[0];
/*    x[0]=x[2];
    x[2]=x[1];
    x[1]=t;*/
    x[0]=x[1];
    x[1]=x[2];
    x[2]=t;
}

/*compute the offset to compare or print transposed local data in original input coordinates
  xs matrix dimension size, xl dimension length, xc decomposition offset 
  s: step in computation = number of transposes*/
static void compute_offsets(fft5d_plan plan, int xs[], int xl[], int xc[], int NG[], int s) {
/*    int direction = plan->direction;
    int fftorder = plan->fftorder;*/
    
    int o;
    int pos[3],i;
    int *pM=plan->pM, *pK=plan->pK, *oM=plan->oM, *oK=plan->oK,
        *C=plan->C, *rC=plan->rC;

    NG[0]=plan->NG;NG[1]=plan->MG;NG[2]=plan->KG;

    if (!(plan->flags&FFT5D_ORDER_YZ)) {
        switch (s) {
        case 0: o=XYZ; break;
        case 1: o=ZYX; break;
        case 2: o=YZX; break;
        default: assert(0);
        }
    } else {
        switch (s) {
        case 0: o=XYZ; break;
        case 1: o=YXZ; break;
        case 2: o=ZXY; break;
        default: assert(0);
        }
    }
 
    switch (o) {
        case XYZ:pos[0]=1;pos[1]=2;pos[2]=3;break;
        case XZY:pos[0]=1;pos[1]=3;pos[2]=2;break;
        case YXZ:pos[0]=2;pos[1]=1;pos[2]=3;break;
        case YZX:pos[0]=3;pos[1]=1;pos[2]=2;break;
        case ZXY:pos[0]=2;pos[1]=3;pos[2]=1;break;
        case ZYX:pos[0]=3;pos[1]=2;pos[2]=1;break;
    }
    /*if (debug) printf("pos: %d %d %d\n",pos[0],pos[1],pos[2]);*/
        
    /*xs, xl give dimension size and data length in local transposed coordinate system
      for 0(/1/2): x(/y/z) in original coordinate system*/
    for (i=0;i<3;i++) {
        switch (pos[i]) {
        case 1: xs[i]=1;         xc[i]=0;     xl[i]=C[s];break;
        case 2: xs[i]=C[s];      xc[i]=oM[s]; xl[i]=pM[s];break;
        case 3: xs[i]=C[s]*pM[s];xc[i]=oK[s]; xl[i]=pK[s];break;
        }
    }
    /*input order is different for test program to match FFTW order 
      (important for complex to real)*/
    if (plan->flags&FFT5D_BACKWARD) {
        rotate(xs);
        rotate(xl);
        rotate(xc);
        rotate(NG);
        if (plan->flags&FFT5D_ORDER_YZ) {
            rotate(xs);
            rotate(xl);
            rotate(xc);
            rotate(NG);            
        }
    }
    if (plan->flags&FFT5D_REALCOMPLEX && ((!(plan->flags&FFT5D_BACKWARD) && s==0) || (plan->flags&FFT5D_BACKWARD && s==2))) {
        xl[0] = rC[s];
    }
}

static void print_localdata(const t_complex* lin, const char* txt, int s, fft5d_plan plan) {
    int x,y,z,l;
    int *coor = plan->coor;
    int xs[3],xl[3],xc[3],NG[3];        
    int ll=(plan->flags&FFT5D_REALCOMPLEX)?1:2;
    compute_offsets(plan,xs,xl,xc,NG,s);
    fprintf(debug,txt,coor[0],coor[1],s);
    /*printf("xs: %d %d %d, xl: %d %d %d\n",xs[0],xs[1],xs[2],xl[0],xl[1],xl[2]);*/
    for(z=0;z<xl[2];z++) {
        for(y=0;y<xl[1];y++) {
            fprintf(debug,"%d %d: ",coor[0],coor[1]);
            for (x=0;x<xl[0];x++) {
                for (l=0;l<ll;l++) {
                    fprintf(debug,"%f ",((real*)lin)[(z*xs[2]+y*xs[1])*2+(x*xs[0])*ll+l]);
                }
                fprintf(debug,",");
            }
            fprintf(debug,"\n");
        }
    }
}

void fft5d_execute(fft5d_plan plan,fft5d_time times) {
    t_complex *lin = plan->lin;
    t_complex *lout = plan->lout;

    gmx_fft_t *p1d=plan->p1d;
#ifdef FFT5D_MPI_TRANSPOSE
    FFTW(plan) *mpip=plan->mpip;
#endif
#ifdef GMX_MPI
    MPI_Comm *cart=plan->cart;
#endif

    double time_fft=0,time_local=0,time_mpi[2]={0},time=0;    
    int *N=plan->N,*M=plan->M,*K=plan->K,*pN=plan->pN,*pM=plan->pM,*pK=plan->pK,
        *C=plan->C,*P=plan->P,**iNin=plan->iNin,**oNin=plan->oNin,**iNout=plan->iNout,**oNout=plan->oNout;
    int s=0;
    
    
#ifdef GMX_FFT_FFTW3 
    if (plan->p3d) {
        if (times!=0)
            time=MPI_Wtime();
        FFTW(execute)(plan->p3d); 
        if (times!=0)
            times->fft+=MPI_Wtime()-time;
        return;
    }
#endif

    /*lin: x,y,z*/
    if (plan->flags&FFT5D_DEBUG) print_localdata(lin, "%d %d: copy in lin\n", s, plan);
    for (s=0;s<2;s++) {
        if (times!=0)
            time=MPI_Wtime();
        
        if ((plan->flags&FFT5D_REALCOMPLEX) && !(plan->flags&FFT5D_BACKWARD) && s==0) {
            gmx_fft_many_1d_real(p1d[s],(plan->flags&FFT5D_BACKWARD)?GMX_FFT_COMPLEX_TO_REAL:GMX_FFT_REAL_TO_COMPLEX,lin,lout);
        } else {
            gmx_fft_many_1d(     p1d[s],(plan->flags&FFT5D_BACKWARD)?GMX_FFT_BACKWARD:GMX_FFT_FORWARD,               lin,lout);
        }
        if (times!=0)
            time_fft+=MPI_Wtime()-time;
    
        if (plan->flags&FFT5D_DEBUG) print_localdata(lout, "%d %d: FFT %d\n", s, plan);
        
#ifdef GMX_MPI
        if (GMX_PARALLEL_ENV_INITIALIZED && cart[s] !=0 && P[s]>1 )
        {
            if (times!=0)
                time=MPI_Wtime(); 
            /*prepare for AllToAll
              1. (most outer) axes (x) is split into P[s] parts of size N[s] 
              for sending*/
            splitaxes(lin,lout,N[s],M[s],K[s], pN[s],pM[s],pK[s],P[s],C[s],iNout[s],oNout[s]);

            if (times!=0)
            {
                time_local+=MPI_Wtime()-time;
            
                /*send, recv*/
                time=MPI_Wtime();
            }

#ifdef FFT5D_MPI_TRANSPOSE
            FFTW(execute)(mpip[s]);  
#else
            if ((s==0 && !(plan->flags&FFT5D_ORDER_YZ)) || (s==1 && (plan->flags&FFT5D_ORDER_YZ))) 
                MPI_Alltoall(lin,N[s]*pM[s]*K[s]*sizeof(t_complex)/sizeof(real),GMX_MPI_REAL,lout,N[s]*pM[s]*K[s]*sizeof(t_complex)/sizeof(real),GMX_MPI_REAL,cart[s]);
            else
                MPI_Alltoall(lin,N[s]*M[s]*pK[s]*sizeof(t_complex)/sizeof(real),GMX_MPI_REAL,lout,N[s]*M[s]*pK[s]*sizeof(t_complex)/sizeof(real),GMX_MPI_REAL,cart[s]);
#endif /*FFT5D_MPI_TRANSPOSE*/
            if (times!=0)
                time_mpi[s]=MPI_Wtime()-time;
        }
#endif /*GMX_MPI*/

    
        if (times!=0)
            time=MPI_Wtime();
        /*bring back in matrix form 
          thus make  new 1. axes contiguos
          also local transpose 1 and 2/3 */
        if ((s==0 && !(plan->flags&FFT5D_ORDER_YZ)) || (s==1 && (plan->flags&FFT5D_ORDER_YZ))) 
            joinAxesTrans13(lin,lout,N[s],pM[s],K[s],pN[s],pM[s],pK[s],P[s],C[s+1],iNin[s+1],oNin[s+1]);
        else 
            joinAxesTrans12(lin,lout,N[s],M[s],pK[s],pN[s],pM[s],pK[s],P[s],C[s+1],iNin[s+1],oNin[s+1]);    
        if (times!=0)
            time_local+=MPI_Wtime()-time;
    
        if (plan->flags&FFT5D_DEBUG) print_localdata(lin, "%d %d: tranposed %d\n", s+1, plan);
                
        /*if (debug) print_localdata(lin, "%d %d: transposed x-z\n", N1, M0, K, ZYX, coor);*/
    }    
    
    if (times!=0)
        time=MPI_Wtime();
    if (plan->flags&FFT5D_INPLACE) lout=lin;
    if ((plan->flags&FFT5D_REALCOMPLEX) && (plan->flags&FFT5D_BACKWARD)) {
        gmx_fft_many_1d_real(p1d[s],(plan->flags&FFT5D_BACKWARD)?GMX_FFT_COMPLEX_TO_REAL:GMX_FFT_REAL_TO_COMPLEX,lin,lout);
    } else {
        gmx_fft_many_1d(     p1d[s],(plan->flags&FFT5D_BACKWARD)?GMX_FFT_BACKWARD:GMX_FFT_FORWARD,               lin,lout);
    }

    if (times!=0)
        time_fft+=MPI_Wtime()-time;
    if (plan->flags&FFT5D_DEBUG) print_localdata(lout, "%d %d: FFT %d\n", s, plan);
    /*if (debug) print_localdata(lout, "%d %d: FFT in y\n", N1, M, K0, YZX, coor);*/
    
    if (times!=0) {
        times->fft+=time_fft;
        times->local+=time_local;
        times->mpi2+=time_mpi[1];
        times->mpi1+=time_mpi[0];
    }
}

void fft5d_destroy(fft5d_plan plan) {
    int s;
    for (s=0;s<3;s++) {
        gmx_many_fft_destroy(plan->p1d[s]);
        if (plan->iNin[s]) {
            free(plan->iNin[s]);
            plan->iNin[s]=0;
        }
        if (plan->oNin[s]) {
            free(plan->oNin[s]);
            plan->oNin[s]=0;
        }
        if (plan->iNout[s]) {
            free(plan->iNout[s]);
            plan->iNout[s]=0;
        }
        if (plan->oNout[s]) {
            free(plan->oNout[s]);
            plan->oNout[s]=0;
        }
    }
#ifdef GMX_FFT_FFTW3 
    FFTW_LOCK;
#ifdef FFT5D_MPI_TRANSPOS
    for (s=0;s<2;s++)    
        FFTW(destroy_plan)(plan->mpip[s]);
#endif /* FFT5D_MPI_TRANSPOS */
#endif /* GMX_FFT_FFTW3 */

    /*We can't free lin/lout here - is allocated by gmx_calloc_aligned which can't be freed*/

    
#ifdef FFT5D_THREADS
#ifdef FFT5D_FFTW_THREADS
    FFTW(cleanup_threads)();
#endif
#endif

    free(plan);
}

/*Is this better than direct access of plan? enough data?
  here 0,1 reference divided by which processor grid dimension (not FFT step!)*/
void fft5d_local_size(fft5d_plan plan,int* N1,int* M0,int* K0,int* K1,int** coor) {
    *N1=plan->N[0];
    *M0=plan->M[0];
    *K1=plan->K[0];
    *K0=plan->N[1];
    
    *coor=plan->coor;
}


/*same as fft5d_plan_3d but with cartesian coordinator and automatic splitting 
  of processor dimensions*/
fft5d_plan fft5d_plan_3d_cart(int NG, int MG, int KG, MPI_Comm comm, int P0, int flags, t_complex** rlin, t_complex** rlout) {
    MPI_Comm cart[2]={0};
#ifdef GMX_MPI
    int size=1,prank=0;
    int P[2];
    int coor[2];
    int wrap[]={0,0};
    MPI_Comm gcart;
    int rdim1[] = {0,1}, rdim2[] = {1,0};

    MPI_Comm_size(comm,&size);
    MPI_Comm_rank(comm,&prank);

    if (P0==0) P0 = lfactor(size);
    if (size%P0!=0) {
        if (prank==0) printf("FFT5D: WARNING: Number of processors %d not evenly dividable by %d\n",size,P0);
        P0 = lfactor(size);
    }
        
    P[0] = P0; P[1]=size/P0; /*number of processors in the two dimensions*/
    
    /*Difference between x-y-z regarding 2d decomposition is whether they are 
      distributed along axis 1, 2 or both*/
    
    MPI_Cart_create(comm,2,P,wrap,1,&gcart); /*parameter 4: value 1: reorder*/
    MPI_Cart_get(gcart,2,P,wrap,coor); 
    MPI_Cart_sub(gcart, rdim1 , &cart[0]);
    MPI_Cart_sub(gcart, rdim2 , &cart[1]);
#endif
    return fft5d_plan_3d(NG, MG, KG, cart, flags, rlin, rlout); 
}



/*prints in original coordinate system of data (as the input to FFT)*/
void fft5d_compare_data(const t_complex* lin, const t_complex* in, fft5d_plan plan, int bothLocal, int normalize) {
    int xs[3],xl[3],xc[3],NG[3];
    int x,y,z,l;
    int *coor = plan->coor;
    int ll=2; /*compare ll values per element (has to be 2 for complex)*/
    if (plan->flags&FFT5D_REALCOMPLEX && plan->flags&FFT5D_BACKWARD) 
    {
        ll=1;
    }

    compute_offsets(plan,xs,xl,xc,NG,2);
    if (plan->flags&FFT5D_DEBUG) printf("Compare2\n");
    for (z=0;z<xl[2];z++) {
        for(y=0;y<xl[1];y++) {
            if (plan->flags&FFT5D_DEBUG) printf("%d %d: ",coor[0],coor[1]);
            for (x=0;x<xl[0];x++) {
                for (l=0;l<ll;l++) { /*loop over real/complex parts*/
                    real a,b;
                    a=((real*)lin)[(z*xs[2]+y*xs[1])*2+x*xs[0]*ll+l];
                    if (normalize) a/=plan->rC[0]*plan->rC[1]*plan->rC[2];
                    if (!bothLocal) 
                        b=((real*)in)[((z+xc[2])*NG[0]*NG[1]+(y+xc[1])*NG[0])*2+(x+xc[0])*ll+l];
                    else 
                        b=((real*)in)[(z*xs[2]+y*xs[1])*2+x*xs[0]*ll+l];
                    if (plan->flags&FFT5D_DEBUG) {
                        printf("%f %f, ",a,b);
                    } else {
                        if (fabs(a-b)>2*NG[0]*NG[1]*NG[2]*GMX_REAL_EPS) {
                            printf("result incorrect on %d,%d at %d,%d,%d: FFT5D:%f reference:%f\n",coor[0],coor[1],x,y,z,a,b);
                        }
/*                        assert(fabs(a-b)<2*NG[0]*NG[1]*NG[2]*GMX_REAL_EPS);*/
                    }
                }
                if (plan->flags&FFT5D_DEBUG) printf(",");
            }
            if (plan->flags&FFT5D_DEBUG) printf("\n");
        }
    }
    
}

