#ifndef FFT5D_H_     
#define FFT5D_H_

#include <complex.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef NOGMX
/*#define GMX_MPI*/
/*#define GMX_FFT_FFTW3*/
FILE* debug;
#endif

#ifdef GMX_LIB_MPI
#include <mpi.h>
#elif defined GMX_THREADS
#include "tmpi.h"
#else
#ifndef HAVE_MPICOMM
#ifndef _gmx_parallel_3dfft_h_
typedef void* MPI_Comm;
#endif
#define HAVE_MPICOMM
#endif
double MPI_Wtime();
#endif
#ifdef GMX_THREADS
#include "tmpi.h"
#endif

#include "gmx_fft.h"

#ifdef GMX_FFT_FFTW3
#include <fftw3.h>
#ifdef FFT5D_MPI_TRANSPOSE
#include <fftw3-mpi.h>
#endif
#endif
/* TODO: optional wrapper
#ifdef GMX_FFT_MKL
#include <fftw/fftw3.h>
#ifdef FFT5D_MPI_TRANSPOSE
#include <fftw/fftw3-mpi.h>
#endif
#endif
*/

#ifndef NOGMX
#ifndef GMX_DOUBLE  /*TODO how to not how have to do this GMX specific in here? 
can't be in gmx_parallel_3dfft.h because it has to be also be set when included from fft5d.c */
#define FFT5D_SINGLE
#endif
#endif

typedef t_complex fft5d_type; 
#ifdef FFT5D_SINGLE
#define FFTW(x) fftwf_##x
typedef float fft5d_rtype;  
#define FFT5D_EPS __FLT_EPSILON__
#define FFT5D_MPI_RTYPE MPI_FLOAT 
#else
#define FFTW(x) fftw_##x
typedef double fft5d_rtype; 
#define FFT5D_EPS __DBL_EPSILON__
#define FFT5D_MPI_RTYPE MPI_DOUBLE 
#endif

struct fft5d_time_t {
	double fft,local,mpi1,mpi2;
};
typedef struct fft5d_time_t *fft5d_time;

typedef enum fft5d_flags_t {
	FFT5D_ORDER_YZ=1,
	FFT5D_BACKWARD=2,
	FFT5D_REALCOMPLEX=4,
	FFT5D_DEBUG=8,
	FFT5D_NOMEASURE=16,
	FFT5D_INPLACE=32,
	FFT5D_NOMALLOC=64
} fft5d_flags;

struct fft5d_plan_t {
	fft5d_type *lin;
	fft5d_type *lout;
        gmx_fft_t p1d[3];   /*1D plans*/
#ifdef GMX_FFT_FFTW3 
        FFTW(plan) p2d;  /*2D plan: used for 1D decomposition if FFT supports transposed output*/
        FFTW(plan) p3d;  /*3D plan: used for 0D decomposition if FFT supports transposed output*/
#ifdef FFT5D_MPI_TRANSPOSE
	FFTW(plan) mpip[2];
#endif
#endif
	MPI_Comm cart[2];

    int N[3],M[3],K[3]; /*local length in transposed coordinate system (if not divisisable max)*/
    int pN[3],pM[3], pK[3]; /*local length - not max but length for this processor*/
    int oM[3],oK[3]; /*offset for current processor*/
    int *iNin[3],*oNin[3],*iNout[3],*oNout[3]; /*size for each processor (if divisisable=max) for out(=split) 
						 and in (=join) and offsets in transposed coordinate system*/
    int C[3],rC[3]; /*global length (of the one global axes) */
    /* C!=rC for real<->complex. then C=rC/2 but with potential padding*/
    int P[2]; /*size of processor grid*/
/*	int fftorder;*/
/*	int direction;*/
/*	int realcomplex;*/
	int flags;
    /*int N0,N1,M0,M1,K0,K1;*/
	int NG,MG,KG;
    /*int P[2];*/
	int coor[2];
}; 

typedef struct fft5d_plan_t *fft5d_plan;

void fft5d_execute(fft5d_plan plan,fft5d_time times);
fft5d_plan fft5d_plan_3d(int N, int M, int K, MPI_Comm comm[2], int flags, fft5d_type** lin, fft5d_type** lin2);
void fft5d_local_size(fft5d_plan plan,int* N1,int* M0,int* K0,int* K1,int** coor);
void fft5d_destroy(fft5d_plan plan);
fft5d_plan fft5d_plan_3d_cart(int N, int M, int K, MPI_Comm comm, int P0, int flags, fft5d_type** lin, fft5d_type** lin2);
void fft5d_compare_data(const fft5d_type* lin, const fft5d_type* in, fft5d_plan plan, int bothLocal, int normarlize);
#endif /*FFTLIB_H_*/
