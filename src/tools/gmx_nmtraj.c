/*
 * $Id$
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
#include "fatal.h"
#include "vec.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "index.h"
#include "pdbio.h"
#include "tpxio.h"
#include "txtdump.h"
#include "physics.h"
#include "random.h"
#include "eigio.h"


int gmx_nmtraj(int argc,char *argv[])
{
    static char *desc[] = 
    {
        "[TT]g_nmtraj[tt] generates an virtual trajectory from an eigenvector, ",
        "corresponding to a harmonic cartesian oscillation around the average ",
        "structure. The eigenvectors should normally be mass-weighted, but you can ",
        "use non-weighted eigenvectors to generate orthogonal motions. ",
        "The output frames are written as a trajectory file covering an entire period, and ",
        "the first frame is the average structure. If you write the trajectory in (or convert to) ",
        "PDB format you can view it directly in PyMol and also render a photorealistic movie. ",
        "Motion amplitudes are calculated from the eigenvalues and a preset temperature, ",
        "assuming equipartition of the energy over all modes. To make the motion clearly visible ",
        "in PyMol you might want to amplify it by setting an unrealistic high temperature. ", 
        "However, be aware that both the linear cartesian displacements and mass weighting will ",
        "lead to serious structure deformation for high amplitudes - this is is simply a limitation ",
        "of the cartesian normal mode model. By default the selected eigenvector is set to 7, since ",
        " the first six normal modes are the translational and rotational degrees of freedom." 
    };

    static int  out_eignr=7;
    static real amplitude=0.25;
    static int  nframes=30;
    static real temp=300.0;
    t_pargs pa[] =
    {
        { "-eignr",     FALSE, etINT,  {&out_eignr}, "Eigenvector to use (first is 1)" },
        { "-temp",      FALSE, etREAL, {&temp},      "Temperature in Kelvin" },
        { "-amplitude", FALSE, etREAL, {&amplitude}, "Amplitude for modes with eigenvalue<=0" },
        { "-nframes",   FALSE, etINT,  {&nframes},   "Number of frames to generate" }
    };
    
#define NPA asize(pa)
  
  int        out;
  t_topology top;
  t_atoms    *atoms;
  rvec       *xtop,*xref,*xav,*xout;
  int        nvec,*eignr=NULL;
  int        *eigvalnr;
  rvec       **eigvec=NULL;
  matrix     box;
  int        natoms;
  int        i,j,d,s,v;
  bool       bDMR,bDMA,bFit;
  char *     indexfile;
  char *     grpname;
  real *     eigval;
  int        neigval;
  int *      dummy;
  real *     invsqrtm;
  char       title[STRLEN];
  real       fraction;
  int        out_eigidx;
  real       out_eigval;
  rvec *     out_eigvec;
  real       omega,Ekin,sum,m,vel;
  bool       found;

  t_filenm fnm[] = 
  { 
      { efTPS, NULL,    NULL,          ffREAD },
      { efTRN, "-v",    "eigenvec",    ffREAD  },
      { efTRX, "-o",    "nmtraj",      ffWRITE }
  }; 
  
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL); 

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&xtop,NULL,box,bDMA);

  read_eigenvectors(opt2fn("-v",NFILE,fnm),&natoms,&bFit,
		    &xref,&bDMR,&xav,&bDMA,&nvec,&eignr,&eigvec,&eigval);
  
  atoms=&top.atoms;

  if(atoms->nr != natoms)
  {
      gmx_fatal(FARGS,"Different number of atoms in topology and eigenvectors.\n");
  }
  
  snew(dummy,natoms);
  for(i=0;i<natoms;i++)
      dummy[i]=i;

  /* Find the eigenvalue/vector to match our select one */
  out_eigidx = 0;
  found      = FALSE;
  
  while(found==FALSE && out_eigidx<nvec)
  {
      if( (out_eignr-1) == eignr[out_eigidx] )
          found = TRUE;
      else
          out_eigidx++;
  }

  if(!found)
  {
      gmx_fatal(FARGS,"Eigenvector #%d not found in input files.\n",out_eignr);
  }
  
  snew(invsqrtm,natoms);
  
  if (bDMA) 
  {
      for(i=0; (i<natoms); i++)
          invsqrtm[i] = invsqrt(atoms->atom[i].m);
  }
  else 
  {
      for(i=0; (i<natoms); i++)
          invsqrtm[i]=1.0;
  }

  snew(xout,natoms);

  out_eigval = eigval[out_eigidx];
  out_eigvec = eigvec[out_eigidx];
  
  if(out_eigval<0)
      out_eigval = 0;

  
  if( (out_eignr >= 7) && (out_eigval > 0))
  {
      /* Derive amplitude from temperature and eigenvalue if we can */
      
      /* Convert eigenvalue to angular frequency, in units s^(-1) */
      omega = sqrt(out_eigval*1.0E21/(AVOGADRO*AMU));
      
      /* Harmonic motion will be x=x0 + A*sin(omega*t)*eigenvec.
      * For t =2*pi*n, all energy will be kinetic, and v=A*omega*eigenvec.
      * The kinetic energy will be sum(0.5*mass*v*v) if we temporarily set A to 1.
      */
      Ekin = 0;
      for(i=0;i<natoms;i++)
      {
          m = atoms->atom[i].m;
          for(d=0;d<DIM;d++)
          {
              vel   = omega*out_eigvec[i][d];
              Ekin += 0.5*m*vel*vel;
          }
      }
      /* Convert Ekin from amu*(nm/s)^2 to J.
       * This will also be proportional to A^2 
       */   
      Ekin *= AMU*1E-18;
      
      /* Set the amplitude so the energy is kT/2 */
      amplitude = sqrt(0.5*BOLTZMANN*temp/Ekin);
  }
      
      
  fprintf(stderr,"Generation normal mode trajectory A=%g...\n",amplitude);

  out=open_trx(ftp2fn(efTRX,NFILE,fnm),"w");
  
  /* Write a sine oscillation around the average structure, 
   * modulated by the eigenvector with selected amplitude.
   */
  
  for(i=0;i<nframes;i++)
  {
      fraction = (real)i/(real)nframes;
      for(j=0;j<natoms;j++)
      {
          for(d=0;d<DIM;d++)
          {
              xout[j][d] = xav[j][d] + amplitude*sin(2*M_PI*fraction)*out_eigvec[j][d];
          }
      }
      write_trx(out,natoms,dummy,atoms,i,(real)i/(real)nframes,box,xout,NULL);      
  }
  
  fprintf(stderr,"\n");
  close_trx(out);
  
  return 0;
}
  
