#ifndef _fftgrid_h
#define _fftgrid_h

#include "typedefs.h"
#include "fftw.h"
#include "complex.h"

#ifdef USE_SGI_FFT

/* Use SGI optimized routines */

#include <fft.h>
typedef real      t_fft_tp;
#define GR_ASSIGN(gr,val)        gr=val
#define GR_INC(gr,val)           gr+=val
#define GR_MULT(gr,val)          gr*=val
#define GR_VALUE(gr)             gr
#define INDEX(i,j,k)             ((i)+(j)*la1+(k)*la12)

#else

/* Use FFTW */

typedef t_complex t_fft_tp;
#define GR_ASSIGN(gr,val)        gr.re=val
#define GR_INC(gr,val)           gr.re+=val
#define GR_MULT(gr,val)          gr.re*=val
#define GR_VALUE(gr)             gr.re
#define INDEX(i,j,k)             ((i)*la12+(j)*la2+k)

#endif

typedef struct {
  t_fft_tp *ptr;
  int      nx,ny,nz,la1,la2,la12;
  int      nptr,nxyz;
} t_fftgrid;


/* Routines ... */
extern t_fftgrid *mk_fftgrid(int nx,int ny,int nz);
/* Create an FFT grid (1 Dimensional), to be indexed by the INDEX macro */

extern void gmxfft3D(FILE *fp,bool bVerbose,t_fftgrid *grid,int dir);
/* Do the FFT, direction may be either 
 * FFTW_FORWARD (sign -1) for real -> complex transform 
 * FFTW_BACKWARD (sign 1) for complex -> real transform
 */
 
extern void clear_fftgrid(t_fftgrid *grid);
/* Set it to zero */

extern void unpack_fftgrid(t_fftgrid *grid,int *nx,int *ny,int *nz,
			   int *la1,int *la2,int *la12,t_fft_tp **ptr);
/* Get the values for the constants into local copies */

extern void print_fftgrid(FILE *out,char *title,t_fftgrid *grid,
			  real factor,char *pdb,rvec box,bool bReal);
/* Print the fftgrid to either the out FILE, or if pdb != NULL, to a file
 * named pdb. All atoms are multiplied by factor. If bReal then the
 * real component is printed, otherwise the imaginary component (pdb only)
 * to the out file, both are printed (if complex at all)
 */

#endif
