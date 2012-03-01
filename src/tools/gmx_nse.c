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

#include <ctype.h>
#include <string.h>
#include "smalloc.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "macros.h"
#include "vec.h"
#include "pbc.h"
#include "xvgr.h"
#include "copyrite.h"
#include "futil.h"
#include "statutil.h"
#include "tpxio.h"
#include "index.h"
#include "gstat.h"
#include "matio.h"
#include "gmx_ana.h"
#include "nsfactor.h"

#ifndef PATH_MAX
#if (defined WIN32 || defined _WIN32 || defined WIN64 || defined _WIN64)
#ifdef MAX_PATH
#define PATH_MAX MAX_PATH
#else
#define PATH_MAX 260
#endif
#endif
#endif

int gmx_nse(int argc,char *argv[])
{
    const char *desc[] = {
        "This is simple tool to compute NSE spectra using Debye formula",
        "It currently uses topology file (since it need to assigne element for each atom)",
        "and trajecory file",
        "[PAR]",
        "[TT]-sqt[TT] Computes NSE intensity curve for each q value",
        "[PAR]",
        "[TT]-startq[TT] Starting q value in nm",
        "[PAR]",
        "[TT]-endq[TT] Ending q value in nm",
        "[PAR]",
        "[TT]-qstep[TT] Stepping in q space",
        "[PAR]",
        "Note: This tools produces large number of sqt files (one file per needed q value)!"
    };
    static gmx_bool bPBC=TRUE;
    static real binwidth=0.2,grid=0.05; /* bins shouldnt be smaller then bond (~0.1nm) length */
    static real start_q=0.01, end_q=2.0, q_step=0.01;
    static gmx_large_int_t  nmc=1048576;
    static unsigned int  seed=0;

    static const char *emode[]= { NULL, "direct", "mc", NULL };
    static const char *emethod[]={ NULL, "debye", "fft", NULL };

    gmx_nentron_atomic_structurefactors_t    *gnsf;
    gmx_nse_t              *gnse;

#define NPA asize(pa)

    t_pargs pa[] = {
        { "-bin", FALSE, etREAL, {&binwidth},
          "Binwidth (nm)" },
        { "-mode", FALSE, etENUM, {emode},
          "Mode for sans spectra calculation" },
        { "-nmc", FALSE, etINT, {&nmc},
          "Number of iterations for Monte-Carlo run"},
        { "-method", FALSE, etENUM, {emethod},
          "[HIDDEN]Method for sans spectra calculation" },
        { "-pbc", FALSE, etBOOL, {&bPBC},
          "Use periodic boundary conditions for computing distances" },
        { "-grid", FALSE, etREAL, {&grid},
          "[HIDDEN]Grid spacing (in nm) for FFTs" },
        {"-startq", FALSE, etREAL, {&start_q},
          "Starting q (1/nm) "},
        {"-endq", FALSE, etREAL, {&end_q},
          "Ending q (1/nm)"},
        { "-qstep", FALSE, etREAL, {&q_step},
          "Stepping in q (1/nm)"},
        { "-seed",     FALSE, etINT,  {&seed},
          "Random seed for Monte-Carlo"},
    };
  FILE      *fp;
  const char *fnTPX,*fnTRX,*fnNDX,*fnDAT=NULL;
  t_trxstatus *status;
  t_topology *top=NULL;
  t_atom    *atom=NULL;
  t_atoms   *atoms=NULL;
  gmx_rmpbc_t  gpbc=NULL;
  gmx_bool  bTPX;
  gmx_bool  bFFT=FALSE, bDEBYE=FALSE;
  gmx_bool  bMC=FALSE, bDIRECT=FALSE;
  int        ePBC=-1;
  matrix     box;
  char       title[STRLEN];
  rvec       *x,*xf;
  int       natoms;
  int       nframes;
  int       nralloc=1;
  real       t;
  char       **grpname=NULL;
  atom_id    *index=NULL;
  int        isize;
  int         i,j;
  char        *sqtf_base=NULL, *sqtf=NULL, *hdr=NULL;
  const char  *sqtf_ext=NULL;
  output_env_t oenv;

#define NFILE asize(fnm)

  t_filenm   fnm[] = {
      { efTPX,  "-s",         NULL,   ffREAD },
      { efTRX,  "-f",         NULL,   ffREAD },
      { efNDX,  NULL,         NULL,   ffOPTRD },
      { efDAT,  "-d",   "nsfactor",   ffOPTRD },
      { efXVG, "-sqt",       "sqt",   ffWRITE }
  };

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_CAN_TIME | PCA_TIME_UNIT | PCA_BE_NICE,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);

  /* Now try to parse opts for modes */
  switch(emethod[0][0]) {
  case 'd':
      bDEBYE=TRUE;
      switch(emode[0][0]) {
      case 'd':
          bDIRECT=TRUE;
          break;
      case 'm':
          bMC=TRUE;
          break;
      default:
          break;
      }
      break;
  case 'f':
      bFFT=TRUE;
      fprintf(stderr,"Using FFT method\n");
      break;
  default:
      break;
  }

  if (!bDEBYE && !bFFT)
      gmx_fatal(FARGS,"Unknown method. Set pr or fft!\n");
  if (!bDIRECT && !bMC)
      gmx_fatal(FARGS,"Unknown mode for g(r) method set to direct or mc!");
  /* Try to read files */
  fnDAT = ftp2fn(efDAT,NFILE,fnm);
  fnTPX = ftp2fn(efTPX,NFILE,fnm);
  fnTRX = ftp2fn(efTRX,NFILE,fnm);

  gnsf = gmx_neutronstructurefactors_init(fnDAT);
  fprintf(stderr,"Read %d atom names from %s with neutron scattering parameters\n\n",gnsf->nratoms,fnDAT);

  snew(top,1);
  snew(gnse,1);
  snew(grpname,1);
  snew(index,1);

  bTPX=read_tps_conf(fnTPX,title,top,&ePBC,&x,NULL,box,TRUE);

  atoms=&(top->atoms);

  printf("\nPlease select group for SANS spectra calculation:\n");
  get_index(&(top->atoms),ftp2fn_null(efNDX,NFILE,fnm),1,&isize,&index,grpname);

  gnse->sans = gmx_sans_init(top,gnsf);

  /* Prepare reference frame */
  if (bPBC) {
      gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
  }

  if (bDEBYE) {
      if (bMC) {
          fprintf(stderr,"Using Monte Carlo Debye method to calculate spectrum\n");
          if (nmc>(gmx_large_int_t)floor(0.5*isize*(isize-1))) {
              fprintf(stderr,"Number of mc iteration larger then number of pairs in index group. Switching to direct method!\n");
              bMC=FALSE;
              bDIRECT=TRUE;
          }
      } else if (bDIRECT) {
          fprintf(stderr,"Using direct Debye method to calculate spectrum\n");
      } else {
          gmx_fatal(FARGS,"Unknown method!\n");
      }
  } else if (bFFT) {
      gmx_fatal(FARGS,"Not implented!\n");
  } else {
      gmx_fatal(FARGS,"Whats this!?\n");
  }

  natoms=read_first_x(oenv,&status,fnTRX,&t,&xf,box);
  if (natoms != atoms->nr)
      fprintf(stderr,"\nWARNING: number of atoms in tpx (%d) and trajectory (%d) do not match\n",natoms,atoms->nr);
  /* copy xf to x in case we dont want to use frame t=0 */
  x = xf;
  /* realy do calc */
  nframes=0;
  snew(gnse->gr,nralloc);
  snew(gnse->sq,nralloc);
  snew(gnse->t,nralloc);
  gnse->gr[0] = calc_radial_distribution_histogram(gnse->sans,x,xf,box,index,isize,binwidth,bMC,nmc,seed);
  gnse->sq[0] = convert_histogram_to_intensity_curve(gnse->gr[nframes],start_q,end_q,q_step);
  gnse->t[0]=t;
  do {
      nframes++;
      if (bPBC) {
          gmx_rmpbc(gpbc,atoms->nr,box,xf);
      }
      gnse->nrframes = nframes;
      if(nralloc<(nframes+1)) {
          nralloc++;
          srenew(gnse->gr,nralloc);
          srenew(gnse->sq,nralloc);
          srenew(gnse->t,nralloc);
      }
      gnse->gr[nframes] = calc_radial_distribution_histogram(gnse->sans,x,xf,box,index,isize,binwidth,bMC,nmc,seed);
      gnse->sq[nframes] = convert_histogram_to_intensity_curve(gnse->gr[nframes],start_q,end_q,q_step);
  } while (read_next_x(oenv,status,&t,natoms,xf,box));
  close_trj(status);

  for(i=1;i<gnse->nrframes;i++) {
      gnse->t[i]= gnse->t[0] + i * (t - gnse->t[0])/(gnse->nrframes - 1);
  }

  snew(gnse->sqt,gnse->sq[0]->qn);

  /* now we will gather s(q(t)) from s(q) spectrums */
  for(i=0;i<gnse->sq[0]->qn;i++) {
      snew(gnse->sqt[i],1);
      gnse->sqt[i]->q = gnse->sq[0]->q[i];
      snew(gnse->sqt[i]->s,gnse->nrframes);
      for(j=0;j<gnse->nrframes;j++) {
          gnse->sqt[i]->s[j] = gnse->sq[j]->s[i];
      }
  }

  /* prepare filenames for sqt output */
  sqtf_ext = strrchr(opt2fn_null("-sqt",NFILE,fnm),'.');
  if (sqtf_ext == NULL)
      gmx_fatal(FARGS,"Output file name '%s' does not contain a '.'",opt2fn_null("-sqt",NFILE,fnm));
  sqtf_base = strdup(opt2fn_null("-sqt",NFILE,fnm));
  sqtf_base[sqtf_ext - opt2fn_null("-sqt",NFILE,fnm)] = '\0';

  /* actualy print data */
  for(i=0;i<gnse->sq[0]->qn;i++) {
      srenew(sqtf,PATH_MAX);
      srenew(hdr,20);
      sprintf(sqtf,"%s_q_%2.2lf%s",sqtf_base,gnse->sqt[i]->q,sqtf_ext);
      sprintf(hdr,"s(q(t)) , q = %2.2lf",gnse->sqt[i]->q);
      fp = xvgropen(sqtf,hdr,"S","q, nm^-1",oenv);
      for(j=0;j<gnse->nrframes;j++) {
          fprintf(fp,"%f\t%f\n",gnse->t[j],gnse->sqt[i]->s[j]);
      }
      fclose(fp);
  }

  thanx(stderr);

  return 0;
}
