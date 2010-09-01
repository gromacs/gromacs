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
#include "gmx_fatal.h"
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
#include "gmx_ana.h"


int gmx_nmtraj(int argc,char *argv[])
{
    const char *desc[] = 
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

    static real refamplitude=0.25;
    static int  nframes=30;
    static real temp=300.0;
    static const char *eignrvec = "7";
    static const char *phasevec  = "0.0";
    
    t_pargs pa[] =
    {
        { "-eignr",     FALSE, etSTR,  {&eignrvec}, "String of eigenvectors to use (first is 1)" },
        { "-phases",    FALSE, etSTR,  {&phasevec}, "String of phases (default is 0.0)" },
        { "-temp",      FALSE, etREAL, {&temp},      "Temperature in Kelvin" },
        { "-amplitude", FALSE, etREAL, {&refamplitude}, "Amplitude for modes with eigenvalue<=0" },
        { "-nframes",   FALSE, etINT,  {&nframes},   "Number of frames to generate" }
    };
    
#define NPA asize(pa)
  
  t_trxstatus *out;
  t_topology top;
  int        ePBC;
  t_atoms    *atoms;
  rvec       *xtop,*xref,*xav,*xout;
  int        nvec,*eignr=NULL;
  int        *eigvalnr;
  rvec       **eigvec=NULL;
  matrix     box;
  int        natoms;
  int        i,j,k,kmode,d,s,v;
  gmx_bool       bDMR,bDMA,bFit;
  char *     indexfile;
  
  char *     grpname;
  real *     eigval;
  int        neigval;
  int *      dummy;
  real *     invsqrtm;
  char       title[STRLEN];
  real       fraction;
  int        *out_eigidx;
  real       *out_eigval;
  rvec *     this_eigvec;
  real       omega,Ekin,sum,m,vel;
  gmx_bool       found;
  int        nmodes,nphases;
  int        *imodes;
  real       *amplitude;
  real       *phases;
  real       dum;
  const char       *p;
  char *pe;  
  output_env_t oenv;
    
  t_filenm fnm[] = 
  { 
      { efTPS, NULL,    NULL,          ffREAD },
      { efTRN, "-v",    "eigenvec",    ffREAD  },
      { efTRO, "-o",    "nmtraj",      ffWRITE }
  }; 
  
#define NFILE asize(fnm) 

  CopyRight(stderr,argv[0]); 
  parse_common_args(&argc,argv,PCA_BE_NICE,
		    NFILE,fnm,NPA,pa,asize(desc),desc,0,NULL,&oenv); 

  read_eigenvectors(opt2fn("-v",NFILE,fnm),&natoms,&bFit,
		    &xref,&bDMR,&xav,&bDMA,&nvec,&eignr,&eigvec,&eigval);

  read_tps_conf(ftp2fn(efTPS,NFILE,fnm),title,&top,&ePBC,&xtop,NULL,box,bDMA);
	
  /* Find vectors and phases */
  
  /* first find number of args in string */
  nmodes=0;
  p=eignrvec;
  while(*p!=0)
  {
      dum=strtod(p,&pe);
      p=pe;
      nmodes++;
  }

  snew(imodes,nmodes);
  p=eignrvec;
  for(i=0;i<nmodes;i++)
  {
	  /* C indices start on 0 */
      imodes[i]=strtol(p,&pe,10)-1;
      p = pe;
  }
 
  /* Now read phases */
  nphases=0;
  p=phasevec;
  while(*p!=0)
  {
      dum=strtod(p,&pe);
      p=pe;
      nphases++;
  }
  if(nphases>nmodes)
  {
      gmx_fatal(FARGS,"More phases than eigenvector indices specified.\n");
  }
    
  snew(phases,nmodes);
  p=phasevec;

  for(i=0;i<nphases;i++)
  {
      phases[i]=strtod(p,&pe);
      p = pe;
  }

  if(nmodes>nphases)
  {
      printf("Warning: Setting phase of last %d modes to zero...\n",nmodes-nphases);
  }
    
  for(i=nphases;i<nmodes;i++)
  {
      phases[i]=0;
  }
        
  atoms=&top.atoms;

  if(atoms->nr != natoms)
  {
      gmx_fatal(FARGS,"Different number of atoms in topology and eigenvectors.\n");
  }
  
  snew(dummy,natoms);
  for(i=0;i<natoms;i++)
      dummy[i]=i;

  /* Find the eigenvalue/vector to match our select one */ 
  snew(out_eigidx,nmodes);
  for(i=0;i<nmodes;i++)
    out_eigidx[i]=-1;
  
  for(i=0;i<nvec;i++)
  {
      for(j=0;j<nmodes;j++)
      {
          if(imodes[j]==eignr[i])
              out_eigidx[j]=i;
      }
  }
    for(i=0;i<nmodes;i++)
        if(out_eigidx[i]==-1)
            gmx_fatal(FARGS,"Could not find mode %d in eigenvector file.\n",imodes[i]);
    
  
  snew(invsqrtm,natoms);
  
  if (bDMA) 
  {
      for(i=0; (i<natoms); i++)
          invsqrtm[i] = gmx_invsqrt(atoms->atom[i].m);
  }
  else 
  {
      for(i=0; (i<natoms); i++)
          invsqrtm[i]=1.0;
  }

  snew(xout,natoms);
  snew(amplitude,nmodes);

  	printf("mode phases: %g %g\n",phases[0],phases[1]);
	
  for(i=0;i<nmodes;i++)
  {
      kmode = out_eigidx[i];
      this_eigvec=eigvec[kmode];
   
      if( (kmode >= 6) && (eigval[kmode] > 0))
      {		  		  
          /* Derive amplitude from temperature and eigenvalue if we can */
          
          /* Convert eigenvalue to angular frequency, in units s^(-1) */
          omega = sqrt(eigval[kmode]*1.0E21/(AVOGADRO*AMU));
          /* Harmonic motion will be x=x0 + A*sin(omega*t)*eigenvec.
           * The velocity is thus:
           * 
           * v = A*omega*cos(omega*t)*eigenvec.
           *
           * And the average kinetic energy the integral of mass*v*v/2 over a
           * period:
           *
           * (1/4)*mass*A*omega*eigenvec
           *
           * For t =2*pi*n, all energy will be kinetic, and v=A*omega*eigenvec.
           * The kinetic energy will be sum(0.5*mass*v*v) if we temporarily set A to 1,
		   * and the average over a period half of this.
           */
          
          Ekin = 0;
          for(k=0;k<natoms;k++)
          {
              m = atoms->atom[k].m;
              for(d=0;d<DIM;d++)
              {
                  vel   = omega*this_eigvec[k][d];
                  Ekin += 0.5*0.5*m*vel*vel;
              }
          }
		  
          /* Convert Ekin from amu*(nm/s)^2 to J, i.e., kg*(m/s)^2
           * This will also be proportional to A^2 
           */   
          Ekin *= AMU*1E-18;
          
          /* Set the amplitude so the energy is kT/2 */
          amplitude[i] = sqrt(0.5*BOLTZMANN*temp/Ekin);		  
	  }
	  else
	  {
		  amplitude[i] = refamplitude;
	  }
  }
	
  out=open_trx(ftp2fn(efTRO,NFILE,fnm),"w");
	
    /* Write a sine oscillation around the average structure, 
     * modulated by the eigenvector with selected amplitude.
     */
    
    for(i=0;i<nframes;i++)
    {
        fraction = (real)i/(real)nframes;
		for(j=0;j<natoms;j++)
		{
			copy_rvec(xav[j],xout[j]);
		}
		
        for(k=0;k<nmodes;k++)
        {
            kmode=out_eigidx[k];
            this_eigvec=eigvec[kmode];

            for(j=0;j<natoms;j++)
            {
                for(d=0;d<DIM;d++)
                {
					xout[j][d] += amplitude[k]*sin(2*M_PI*(fraction+phases[k]/360.0))*this_eigvec[j][d];
                }
            }
        }
        write_trx(out,natoms,dummy,atoms,i,(real)i/(real)nframes,box,xout,NULL,NULL);
    }    
    
    fprintf(stderr,"\n");
    close_trx(out);
    
    return 0;
}
  

