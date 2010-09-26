/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.0
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <math.h>
#include <string.h>

#include "statutil.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "mshift.h"
#include "xvgr.h"
#include "gstat.h"
#include "txtdump.h"
#include "eigensolver.h"
#include "eigio.h"
#include "mtxio.h"
#include "sparsematrix.h"
#include "physics.h"
#include "main.h"
#include "gmx_ana.h"


static void
nma_full_hessian(real *           hess,
                 int              ndim,
                 gmx_bool             bM,
                 t_topology *     top,
                 int              begin,
                 int              end,
                 real *           eigenvalues,
                 real *           eigenvectors)
{
    int i,j,k,l;
    real mass_fac,rdum;
    int natoms;
    
    natoms = top->atoms.nr;

    /* divide elements hess[i][j] by sqrt(mas[i])*sqrt(mas[j]) when required */

    if (bM)
    {
        for (i=0; (i<natoms); i++) 
        {
            for (j=0; (j<DIM); j++) 
            {
                for (k=0; (k<natoms); k++) 
                {
                    mass_fac=gmx_invsqrt(top->atoms.atom[i].m*top->atoms.atom[k].m);
                    for (l=0; (l<DIM); l++)
                        hess[(i*DIM+j)*ndim+k*DIM+l]*=mass_fac;
                }
            }
        }
    }
    
    /* call diagonalization routine. */
    
    fprintf(stderr,"\nDiagonalizing to find vectors %d through %d...\n",begin,end);
    fflush(stderr);
    
    eigensolver(hess,ndim,begin-1,end-1,eigenvalues,eigenvectors);

    /* And scale the output eigenvectors */
    if (bM && eigenvectors!=NULL)
    {
        for(i=0;i<(end-begin+1);i++)
        {
            for(j=0;j<natoms;j++)
            {
                mass_fac = gmx_invsqrt(top->atoms.atom[j].m);
                for (k=0; (k<DIM); k++) 
                {
                    eigenvectors[i*ndim+j*DIM+k] *= mass_fac;
                }
            }
        }
    }
}



static void
nma_sparse_hessian(gmx_sparsematrix_t *     sparse_hessian,
                   gmx_bool                     bM,
                   t_topology *             top,
                   int                      neig,
                   real *                   eigenvalues,
                   real *                   eigenvectors)
{
    int i,j,k;
    int row,col;
    real mass_fac;
    int iatom,katom;
    int natoms;
    int ndim;
    
    natoms = top->atoms.nr;
    ndim   = DIM*natoms;
    
    /* Cannot check symmetry since we only store half matrix */
    /* divide elements hess[i][j] by sqrt(mas[i])*sqrt(mas[j]) when required */
    
    if (bM)
    {
        for (iatom=0; (iatom<natoms); iatom++) 
        {
            for (j=0; (j<DIM); j++) 
            {
                row = DIM*iatom+j;
                for(k=0;k<sparse_hessian->ndata[row];k++)
                {
                    col = sparse_hessian->data[row][k].col;
                    katom = col/3;
                    mass_fac=gmx_invsqrt(top->atoms.atom[iatom].m*top->atoms.atom[katom].m);
                    sparse_hessian->data[row][k].value *=mass_fac;
                }
            }
        }
    }
    fprintf(stderr,"\nDiagonalizing to find eigenvectors 1 through %d...\n",neig);
    fflush(stderr);
        
    sparse_eigensolver(sparse_hessian,neig,eigenvalues,eigenvectors,10000000);

    /* Scale output eigenvectors */
    if (bM && eigenvectors!=NULL)
    {
        for(i=0;i<neig;i++)
        {
            for(j=0;j<natoms;j++)
            {
                mass_fac = gmx_invsqrt(top->atoms.atom[j].m);
                for (k=0; (k<DIM); k++) 
                {
                    eigenvectors[i*ndim+j*DIM+k] *= mass_fac;
                }
            }
        }
    }
}



int gmx_nmeig(int argc,char *argv[])
{
  const char *desc[] = {
    "g_nmeig calculates the eigenvectors/values of a (Hessian) matrix,",
    "which can be calculated with [TT]mdrun[tt].",
    "The eigenvectors are written to a trajectory file ([TT]-v[tt]).",
    "The structure is written first with t=0. The eigenvectors",
    "are written as frames with the eigenvector number as timestamp.",
    "The eigenvectors can be analyzed with [TT]g_anaeig[tt].",
    "An ensemble of structures can be generated from the eigenvectors with",
    "[TT]g_nmens[tt]. When mass weighting is used, the generated eigenvectors",
    "will be scaled back to plain cartesian coordinates before generating the",
    "output - in this case they will no longer be exactly orthogonal in the",
    "standard cartesian norm (But in the mass weighted norm they would be)."
  };
    
  static gmx_bool bM=TRUE;
  static int  begin=1,end=50;
  t_pargs pa[] = 
  {
    { "-m",  FALSE, etBOOL, {&bM},
      "Divide elements of Hessian by product of sqrt(mass) of involved "
      "atoms prior to diagonalization. This should be used for 'Normal Modes' "
      "analysis" },
    { "-first", FALSE, etINT, {&begin},     
      "First eigenvector to write away" },
    { "-last",  FALSE, etINT, {&end}, 
      "Last eigenvector to write away" }
  };
  FILE       *out;
  int        status,trjout;
  t_topology top;
  int        ePBC;
  rvec       *top_x;
  matrix     box;
  real       *eigenvalues;
  real       *eigenvectors;
  real       rdum,mass_fac;
  int        natoms,ndim,nrow,ncol,count;
  char       *grpname,title[256];
  int        i,j,k,l,d,gnx;
  gmx_bool       bSuck;
  atom_id    *index;
  real       value;
  real       factor_gmx_to_omega2;
  real       factor_omega_to_wavenumber;
  t_commrec  *cr;
  output_env_t oenv;
  
  real *                 full_hessian   = NULL;
  gmx_sparsematrix_t *   sparse_hessian = NULL;

  t_filenm fnm[] = { 
    { efMTX, "-f", "hessian",    ffREAD  }, 
    { efTPS, NULL, NULL,         ffREAD  },
    { efXVG, "-of", "eigenfreq", ffWRITE },
    { efXVG, "-ol", "eigenval",  ffWRITE },
    { efTRN, "-v", "eigenvec",  ffWRITE }
  }; 
#define NFILE asize(fnm) 

	cr = init_par(&argc,&argv);

	if(MASTER(cr))
		CopyRight(stderr,argv[0]); 
	
  parse_common_args(&argc,argv,PCA_BE_NICE | (MASTER(cr) ? 0 : PCA_QUIET),
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv); 

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&top_x,NULL,box,bM);

  natoms = top.atoms.nr;
  ndim = DIM*natoms;

  if(begin<1)
      begin = 1;
  if(end>ndim)
      end = ndim;

  /*open Hessian matrix */
  gmx_mtxio_read(ftp2fn(efMTX,NFILE,fnm),&nrow,&ncol,&full_hessian,&sparse_hessian);
    
  /* Memory for eigenvalues and eigenvectors (begin..end) */
  snew(eigenvalues,nrow);
  snew(eigenvectors,nrow*(end-begin+1));
       
  /* If the Hessian is in sparse format we can calculate max (ndim-1) eigenvectors,
   * and they must start at the first one. If this is not valid we convert to full matrix
   * storage, but warn the user that we might run out of memory...
   */    
  if((sparse_hessian != NULL) && (begin!=1 || end==ndim))
  {
      if(begin!=1)
      {
          fprintf(stderr,"Cannot use sparse Hessian with first eigenvector != 1.\n");
      }
      else if(end==ndim)
      {
          fprintf(stderr,"Cannot use sparse Hessian to calculate all eigenvectors.\n");
      }
      
      fprintf(stderr,"Will try to allocate memory and convert to full matrix representation...\n");
      
      snew(full_hessian,nrow*ncol);
      for(i=0;i<nrow*ncol;i++)
          full_hessian[i] = 0;
      
      for(i=0;i<sparse_hessian->nrow;i++)
      {
          for(j=0;j<sparse_hessian->ndata[i];j++)
          {
              k     = sparse_hessian->data[i][j].col;
              value = sparse_hessian->data[i][j].value;
              full_hessian[i*ndim+k] = value;
              full_hessian[k*ndim+i] = value;
          }
      }
      gmx_sparsematrix_destroy(sparse_hessian);
      sparse_hessian = NULL;
      fprintf(stderr,"Converted sparse to full matrix storage.\n");
  }
  
  if(full_hessian != NULL)
  {
      /* Using full matrix storage */
      nma_full_hessian(full_hessian,nrow,bM,&top,begin,end,eigenvalues,eigenvectors);
  }
  else
  {
      /* Sparse memory storage, allocate memory for eigenvectors */
      snew(eigenvectors,ncol*end);
      nma_sparse_hessian(sparse_hessian,bM,&top,end,eigenvalues,eigenvectors);
  }
  
  
  /* check the output, first 6 eigenvalues should be reasonably small */  
  bSuck=FALSE;
  for (i=begin-1; (i<6); i++) 
  {
      if (fabs(eigenvalues[i]) > 1.0e-3) 
          bSuck=TRUE;
  }
  if (bSuck) 
  {
      fprintf(stderr,"\nOne of the lowest 6 eigenvalues has a non-zero value.\n");
      fprintf(stderr,"This could mean that the reference structure was not\n");
      fprintf(stderr,"properly energy minimized.\n");
  }
                      
                      
  /* now write the output */
  fprintf (stderr,"Writing eigenvalues...\n");
  out=xvgropen(opt2fn("-ol",NFILE,fnm), 
               "Eigenvalues","Eigenvalue index","Eigenvalue [Gromacs units]",
               oenv);
  if (output_env_get_print_xvgr_codes(oenv)) {
    if (bM)
      fprintf(out,"@ subtitle \"mass weighted\"\n");
    else 
      fprintf(out,"@ subtitle \"not mass weighted\"\n");
  }
  
  for (i=0; i<=(end-begin); i++)
      fprintf (out,"%6d %15g\n",begin+i,eigenvalues[i]);
  ffclose(out);
  

  
  fprintf(stderr,"Writing eigenfrequencies - negative eigenvalues will be set to zero.\n");

  out=xvgropen(opt2fn("-of",NFILE,fnm), 
               "Eigenfrequencies","Eigenvector index","Wavenumber [cm\\S-1\\N]",
               oenv);
  if (output_env_get_print_xvgr_codes(oenv)) { 
    if (bM)
      fprintf(out,"@ subtitle \"mass weighted\"\n");
    else 
      fprintf(out,"@ subtitle \"not mass weighted\"\n");
  }
  
  /* Gromacs units are kJ/(mol*nm*nm*amu),
   * where amu is the atomic mass unit.
   *
   * For the eigenfrequencies we want to convert this to spectroscopic absorption
   * wavenumbers given in cm^(-1), which is the frequency divided by the speed of
   * light. Do this by first converting to omega^2 (units 1/s), take the square 
   * root, and finally divide by the speed of light (nm/ps in gromacs).   
   */
  factor_gmx_to_omega2       = 1.0E21/(AVOGADRO*AMU);
  factor_omega_to_wavenumber = 1.0E-5/(2.0*M_PI*SPEED_OF_LIGHT);  
    
  for (i=0; i<=(end-begin); i++)
  {
      value = eigenvalues[i];
      if(value < 0)
          value = 0;
      value=sqrt(value*factor_gmx_to_omega2)*factor_omega_to_wavenumber;
      fprintf (out,"%6d %15g\n",begin+i,value);
  }
  ffclose(out);
  
  /* Writing eigenvectors. Note that if mass scaling was used, the eigenvectors 
   * were scaled back from mass weighted cartesian to plain cartesian in the
   * nma_full_hessian() or nma_sparse_hessian() routines. Mass scaled vectors
   * will not be strictly orthogonal in plain cartesian scalar products.
   */   
  write_eigenvectors(opt2fn("-v",NFILE,fnm),natoms,eigenvectors,FALSE,begin,end,
                     eWXR_NO,NULL,FALSE,top_x,bM,eigenvalues);
  
  thanx(stderr);
  
  return 0;
}

