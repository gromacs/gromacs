/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.3
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

#ifndef _fftgrid_h
#define _fftgrid_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "gmxcomplex.h"
#include "network.h"
#include "gmx_fft.h"

#ifdef GMX_MPI
#include "gmx_parallel_3dfft.h"
#endif

/* Use FFTW */

#define INDEX(i,j,k)             ((i)*la12+(j)*la2+(k))      

typedef struct {
  int local_nx,local_x_start,local_ny_after_transpose;
  int local_y_start_after_transpose;
} t_parfft;

typedef struct {
    real *                 ptr;
    real *                 localptr;
    real *                 workspace;    
    int                    nx,ny,nz,la2r,la2c,la12r,la12c;
    int                    nptr,nxyz;
    gmx_fft_t              fft_setup;
#ifdef GMX_MPI
    gmx_parallel_3dfft_t   mpi_fft_setup;
    t_parfft               pfft;
#endif
} t_fftgrid;

extern t_fftgrid *mk_fftgrid(FILE *       fp,
                             int          nx,
                             int          ny,
                             int          nz,
                             t_commrec *  cr);

/* Create an FFT grid (1 Dimensional), to be indexed by the INDEX macro 
 * Setup FFT plans and extract local sizes for the grid.
 * If the file pointer is given, information is printed to it.
 * If cr is non-NULL and cr->nnodes>1, a parallel grid and FFT will be created.
 */

extern void pr_fftgrid(FILE *fp,char *title,t_fftgrid *grid);
/* Dump a grid to a file */

extern void done_fftgrid(t_fftgrid *grid);
/* And throw it away again */

extern void gmxfft3D(t_fftgrid *grid,int dir,t_commrec *cr);
/* Do the FFT, direction may be either 
 * FFTW_FORWARD (sign -1) for real -> complex transform 
 * FFTW_BACKWARD (sign 1) for complex -> real transform
 */
 
extern void clear_fftgrid(t_fftgrid *grid);
/* Set it to zero */

extern void unpack_fftgrid(t_fftgrid *grid,int *nx,int *ny,int *nz,
			   int *nx2,int *ny2,int *nz2,
			   int *la2, int *la12,bool bReal, real **ptr);

/* Get the values for the constants into local copies */




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






