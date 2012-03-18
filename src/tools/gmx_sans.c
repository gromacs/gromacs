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


int gmx_sans(int argc,char *argv[])
{
    const char *desc[] = {
        "This is simple tool to compute SANS spectra using Debye formula",
        "It currently uses topology file (since it need to assigne element for each atom)",
        "[PAR]",
        "[TT]-pr[TT] Computes normalized g(r) function",
        "[PAR]",
        "[TT]-sq[TT] Computes SANS intensity curve for needed q diapason",
        "[PAR]",
        "[TT]-startq[TT] Starting q value in nm",
        "[PAR]",
        "[TT]-endq[TT] Ending q value in nm",
        "[PAR]",
        "[TT]-qstep[TT] Stepping in q space",
        "[PAR]",
        "Note: When using Debye direct method computational cost increases as",
        "1/2 * N * (N - 1) where N is atom number in group of interest"
    };
    static gmx_bool bPBC=TRUE;
    static real binwidth=0.2,grid=0.05; /* bins shouldnt be smaller then bond (~0.1nm) length */
    static real start_q=0.0, end_q=2.0, q_step=0.01;
    static gmx_large_int_t  nmc=1048576;
    static unsigned int  seed=0;

    static const char *emode[]= { NULL, "direct", "mc", NULL };
    static const char *emethod[]={ NULL, "debye", "fft", NULL };

    gmx_nentron_atomic_structurefactors_t    *gnsf;
    gmx_sans_t              *gsans;

#define NPA asize(pa)

    t_pargs pa[] = {
        { "-bin", FALSE, etREAL, {&binwidth},
          "[HIDDEN]Binwidth (nm)" },
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
  const char *fnTPX,*fnNDX,*fnDAT=NULL;
  t_trxstatus *status;
  t_topology *top=NULL;
  t_atom    *atom=NULL;
  gmx_rmpbc_t  gpbc=NULL;
  gmx_bool  bTPX;
  gmx_bool  bFFT=FALSE, bDEBYE=FALSE;
  gmx_bool  bMC=FALSE, bDIRECT=FALSE;
  int        ePBC=-1;
  matrix     box;
  char       title[STRLEN];
  rvec       *x;
  int       natoms;
  real       t;
  char       **grpname=NULL;
  atom_id    *index=NULL;
  int        isize;
  int         i,j;
  gmx_radial_distribution_histogram_t  *pr=NULL;
  gmx_static_structurefator_t  *sq=NULL;
  output_env_t oenv;

#define NFILE asize(fnm)

  t_filenm   fnm[] = {
      { efTPX,  "-s",         NULL,   ffREAD },
      { efNDX,  NULL,         NULL,   ffOPTRD },
      { efDAT,  "-d",   "nsfactor",   ffOPTRD },
      { efXVG, "-sq",         "sq",   ffWRITE },
      { efXVG, "-pr",         "pr",   ffWRITE }
  };

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,PCA_BE_NICE,
                    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL,&oenv);

  /* Now try to parse opts for modes */
  switch(emethod[0][0]) {
  case 'd':
      bDEBYE=TRUE;
      switch(emode[0][0]) {
      case 'd':
          bDIRECT=TRUE;
          fprintf(stderr,"Using direct Debye method to calculate spectrum\n");
          break;
      case 'm':
          bMC=TRUE;
          fprintf(stderr,"Using Monte Carlo Debye method to calculate spectrum\n");
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

  gnsf = gmx_neutronstructurefactors_init(fnDAT);
  fprintf(stderr,"Read %d atom names from %s with neutron scattering parameters\n\n",gnsf->nratoms,fnDAT);

  snew(top,1);
  snew(grpname,1);
  snew(index,1);

  bTPX=read_tps_conf(fnTPX,title,top,&ePBC,&x,NULL,box,TRUE);

  printf("\nPlease select group for SANS spectra calculation:\n");
  get_index(&(top->atoms),ftp2fn_null(efNDX,NFILE,fnm),1,&isize,&index,grpname);

  gsans = gmx_sans_init(top,gnsf);

  /* Prepare reference frame */
  if (bPBC) {
      gpbc = gmx_rmpbc_init(&top->idef,ePBC,top->atoms.nr,box);
      gmx_rmpbc(gpbc,top->atoms.nr,box,x);
  }

  natoms=top->atoms.nr;
  if (bDEBYE) {
      if (bDIRECT) {
          /* calc pr */
          pr = calc_radial_distribution_histogram(gsans,x,index,isize,binwidth,bMC,nmc,seed);
      } else if (bMC) {
          if (nmc>(gmx_large_int_t)floor(0.5*isize*(isize-1))) {
              fprintf(stderr,"Number of mc iteration larger then number of pairs in index group. Switching to direct method!\n");
              bMC=FALSE;
              bDIRECT=TRUE;
              pr = calc_radial_distribution_histogram(gsans,x,index,isize,binwidth,bMC,nmc,seed);
          } else {
              pr = calc_radial_distribution_histogram(gsans,x,index,isize,binwidth,bMC,nmc,seed);
          }
      } else {
          gmx_fatal(FARGS,"Unknown method!\n");
      }
  } else if (bFFT) {
      gmx_fatal(FARGS,"Not implented!\n");
  } else {
      gmx_fatal(FARGS,"Whats this!?\n");
  }

  /* prepare pr.xvg */
  fp = xvgropen(opt2fn_null("-pr",NFILE,fnm),"G(r)","Distance (nm)","Probability",oenv);
  for(i=0;i<pr->grn;i++)
      fprintf(fp,"%10.6lf%10.6lf\n",pr->r[i],pr->gr[i]);
  fclose(fp);

  /* prepare sq.xvg */
  sq = convert_histogram_to_intensity_curve(pr,start_q,end_q,q_step);
  fp = xvgropen(opt2fn_null("-sq",NFILE,fnm),"I(q)","q (nm^-1)","s(q)/s(0)",oenv);
  for(i=0;i<sq->qn;i++) {
      fprintf(fp,"%10.6lf%10.6lf\n",sq->q[i],sq->s[i]);
  }
  fclose(fp);

  sfree(pr);

  thanx(stderr);
  
  return 0;
}
