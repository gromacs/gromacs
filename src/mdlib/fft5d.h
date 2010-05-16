#ifndef FFT5D_H_     
#define FFT5D_H_

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef NOGMX
/*#define GMX_MPI*/
/*#define GMX_FFT_FFTW3*/
FILE* debug;
#endif

#include <types/commrec.h>
#include "gmxcomplex.h"
#include "gmx_fft.h"

#ifndef GMX_LIB_MPI
double MPI_Wtime();
#endif

/*currently only special optimization for FFTE*/
#ifdef GMX_FFT_FFTW3
#include <fftw3.h>
#endif

#ifndef GMX_DOUBLE
#define FFTW(x) fftwf_##x
#else
#define FFTW(x) fftw_##x
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
	t_complex *lin;
	t_complex *lout;
        gmx_fft_t p1d[3];   /*1D plans*/
#ifdef GMX_FFT_FFTW3 
        FFTW(plan) p2d;  /*2D plan: used for 1D decomposition if FFT supports transposed output*/
        FFTW(plan) p3d;  /*3D plan: used for 0D decomposition if FFT supports transposed output*/
	FFTW(plan) mpip[2];
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
fft5d_plan fft5d_plan_3d(int N, int M, int K, MPI_Comm comm[2], int flags, t_complex** lin, t_complex** lin2);
void fft5d_local_size(fft5d_plan plan,int* N1,int* M0,int* K0,int* K1,int** coor);
void fft5d_destroy(fft5d_plan plan);
fft5d_plan fft5d_plan_3d_cart(int N, int M, int K, MPI_Comm comm, int P0, int flags, t_complex** lin, t_complex** lin2);
void fft5d_compare_data(const t_complex* lin, const t_complex* in, fft5d_plan plan, int bothLocal, int normarlize);
#endif /*FFTLIB_H_*/
