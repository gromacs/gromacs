#ifndef _fftgrid_h
#define _fftgrid_h

#include <stdio.h>
#include "typedefs.h"
#include "fftw.h"
#ifdef USE_MPI
#include "mpi.h"
#include "fftwnd_mpi.h"
#endif
#include "complex.h"
#include "network.h"

/* Use FFTW */

typedef t_complex t_fft_tp;
#define GR_ASSIGN(gr,val)        gr.re=val
#define GR_INC(gr,val)           gr.re+=val
#define GR_MULT(gr,val)          gr.re*=val
#define GR_VALUE(gr)             gr.re
#define INDEX(i,j,k)             ((i)*la12+(j)*la2+k)

typedef struct {
  int local_nx,local_x_start,local_ny_after_transpose;
  int local_y_start_after_transpose,total_local_size;
} t_parfft;

typedef struct {
  t_fft_tp *ptr;
  int      nx,ny,nz,la1,la2,la12;
  int      nptr,nxyz;
  fftwnd_plan     plan_fw,plan_bw;         /* fw = FORWARD, bw = BACKWARD */
#ifdef USE_MPI
  fftwnd_mpi_plan plan_mpi_fw,plan_mpi_bw;
  t_parfft        pfft_fw,pfft_bw;
#endif
} t_fftgrid;

extern t_fftgrid *mk_fftgrid(FILE *fp,bool bParallel,int nx,int ny,int nz);
/* Create an FFT grid (1 Dimensional), to be indexed by the INDEX macro 
 * Setup FFTW plans and extract local sizes for the grid.
 * If the file pointer is given, information is printed to it.
 */

extern void done_fftgrid(t_fftgrid *grid);
/* And throw it away again */

extern void gmxfft3D(FILE *fp,bool bVerbose,t_fftgrid *grid,int dir,t_commrec *cr);
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


/************************************************************************
 * 
 * For backward compatibility (for testing the ewald code vs. PPPM etc)
 * some old grid routines are retained here.
 *
 ************************************************************************/
 
extern real ***mk_rgrid(int nx,int ny,int nz);

extern void free_rgrid(real ***grid,int nx,int ny);

extern real print_rgrid(FILE *fp,char *title,int nx,int ny,int nz,
			real ***grid);

extern void print_rgrid_pdb(char *fn,int nx,int ny,int nz,real ***grid);

extern t_complex ***mk_cgrid(int nx,int ny,int nz);

extern void free_cgrid(t_complex ***grid,int nx,int ny);

extern t_complex print_cgrid(FILE *fp,char *title,int nx,int ny,int nz,
			   t_complex ***grid);

extern void clear_cgrid(int nx,int ny,int nz,t_complex ***grid);

extern void clear_rgrid(int nx,int ny,int nz,real ***grid);

#endif
