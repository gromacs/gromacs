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
#include <ctype.h>
#include "string2.h"
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
#include "physics.h"
#include "index.h"
#include "smalloc.h"
#include "calcgrid.h"
#include "nrnb.h"
#include "coulomb.h"
#include "gstat.h"
#include "matio.h"
#include "strdb.h"
#include "gmx_ana.h"
#include "names.h"
#include "gmx_random.h"

/*
 * This part will be organized as part of gmxlib in nsfactor.c in near future
 */

typedef struct gmx_nstructurefactors {
    int     nratoms;
    int     *p; /* proton number */
    int     *n; /* neuton number */
    double  *slength; /* scattering length in fm */
    char    **atomnm; /* atom symbol */
} gmx_nstructurefactors;

typedef struct gmx_sans_t {
    t_topology *top; /* topology */
    double *slength; /* scattering length for this topology */
} gmx_sans_t;

typedef struct gmx_gr_t {
    int     grn; /* number of bins */
    double binwidth; /* bin size */
    double *r; /* Distances */
    double *gr; /* Probability */
} gmx_gr_t;

typedef struct gmx_sq_t {
    int     qn; /* number of items */
    double  *s; /* scattering */
    double  *q; /* q vectors */
    double  qstep; /* q increment */
} gmx_sq_t;

void normalize_probability(int n,double *a){
    int i;
    double norm=0.0;
    for (i=0;i<n;i++) norm +=a[i];
    for (i=0;i<n;i++) a[i]/=norm;
}

extern gmx_nstructurefactors *gmx_structurefactors_init(const char *datfn) {
    /* read nsfactor.dat */
    FILE    *fp;
    char    line[STRLEN];
    int     nralloc=10;
    int     n,p;
    int     i, line_no;
    char    atomnm[8];
    double  slength;
    gmx_nstructurefactors   *gnsf;

    fp=libopen(datfn);
    line_no = 0;
    /* allocate memory for structure */
    snew(gnsf,nralloc);
    snew(gnsf->atomnm,nralloc);
    snew(gnsf->p,nralloc);
    snew(gnsf->n,nralloc);
    snew(gnsf->slength,nralloc);

    gnsf->nratoms=line_no;

    while(get_a_line(fp,line,STRLEN)) {
        i=line_no;
        if (sscanf(line,"%s %d %d %lf",atomnm,&p,&n,&slength) == 4) {
            gnsf->atomnm[i]=strdup(atomnm);
            gnsf->n[i]=n;
            gnsf->p[i]=p;
            gnsf->slength[i]=slength;
            line_no++;
            gnsf->nratoms=line_no;
            if (line_no==nralloc){
                nralloc++;
                srenew(gnsf->atomnm,nralloc);
                srenew(gnsf->p,nralloc);
                srenew(gnsf->n,nralloc);
                srenew(gnsf->slength,nralloc);
            }
        } else
            fprintf(stderr,"WARNING: Error in file %s at line %d ignored\n",
                    datfn,line_no);
    }
    srenew(gnsf->atomnm,gnsf->nratoms);
    srenew(gnsf->p,gnsf->nratoms);
    srenew(gnsf->n,gnsf->nratoms);
    srenew(gnsf->slength,gnsf->nratoms);

    fclose(fp);

    return (gmx_nstructurefactors *) gnsf;
}

extern gmx_sans_t *gmx_sans_init (t_topology *top, gmx_nstructurefactors *gnsf) {
    gmx_sans_t    *gsans=NULL;
    int     i,j;
    /* Try to assing scattering length from nsfactor.dat */
    snew(gsans,1);
    snew(gsans->slength,top->atoms.nr);
    /* copy topology data */
    gsans->top = top;
    for(i=0;i<top->atoms.nr;i++) {
        for(j=0;j<gnsf->nratoms;j++) {
            if(top->atoms.atom[i].atomnumber == gnsf->p[j]) {
                /* we need special case for H and D */
                if(top->atoms.atom[i].atomnumber == 1) {
                    if(top->atoms.atom[i].m == 1.008000) {
                        gsans->slength[i] = gnsf->slength[0];
                    } else
                        gsans->slength[i] = gnsf->slength[1];
                } else
                    gsans->slength[i] = gnsf->slength[j];
            }
        }
    }

    return (gmx_sans_t *) gsans;
}

extern gmx_gr_t *calc_pr (gmx_sans_t *gsans, rvec *x, atom_id *index, int isize, double binwidth, gmx_bool bMC, gmx_large_int_t nmc, unsigned int seed) {
    gmx_gr_t    *pr=NULL;
    rvec        xmin, xmax;
    double      rmax;
    int         i,j,d;
    int         mc;
    gmx_rng_t   rng=NULL;

    /* allocate memory for pr */
    snew(pr,1);
    /* set some fields */
    pr->binwidth=binwidth;

    /* Lets try to find min and max distance */
    for(d=0;d<3;d++) {
        xmax[d]=x[index[0]][d];
        xmin[d]=x[index[0]][d];
    }

    for(i=1;i<isize;i++)
        for(d=0;d<3;d++)
            if (xmax[d]<x[index[i]][d]) xmax[d]=x[index[i]][d]; else
                if (xmin[d]>x[index[i]][d]) xmin[d]=x[index[i]][d];

    rmax=sqrt(distance2(xmax,xmin));

    pr->grn=(int)truncf(rmax/pr->binwidth)+1;
    rmax=pr->grn*pr->binwidth;

    snew(pr->gr,pr->grn);

    if(bMC) {
        for(d=0;d<(int)truncf(nmc/524288);d++) {
            rng=gmx_rng_init(seed);
            for(mc=0;mc<524288;mc++) {
                i=(int)truncf(gmx_rng_uniform_real(rng)*isize);
                j=(int)truncf(gmx_rng_uniform_real(rng)*isize);
                if(i!=j)
                    pr->gr[(int)truncf(sqrt(distance2(x[index[i]],x[index[j]]))/binwidth)]+=gsans->slength[index[i]]*gsans->slength[index[j]];
            }
            gmx_rng_destroy(rng);
        }
    } else {
        for(i=0;i<isize;i++)
            for(j=0;j<i;j++)
                pr->gr[(int)truncf(sqrt(distance2(x[index[i]],x[index[j]]))/binwidth)]+=gsans->slength[index[i]]*gsans->slength[index[j]];
    }

    /* normalize */
    normalize_probability(pr->grn,pr->gr);
    snew(pr->r,pr->grn);
    for(i=0;i<pr->grn;i++)
        pr->r[i]=(pr->binwidth*i+pr->binwidth*0.5);

    return (gmx_gr_t *) pr;
}

extern gmx_sq_t *pr2iq (gmx_gr_t *pr, double start_q, double end_q, double q_step) {
    gmx_sq_t    *sq=NULL;
    int         i,j;
    /* init data */
    snew(sq,1);
    sq->qn=(int)truncf((end_q-start_q)/q_step);
    snew(sq->q,sq->qn);
    snew(sq->s,sq->qn);
    for(i=0;i<sq->qn;i++)
        sq->q[i]=start_q+i*q_step;

    if(start_q==0.0) {
        sq->s[0]=1.0;
        for(i=1;i<sq->qn;i++) {
            for(j=0;j<pr->grn;j++)
                sq->s[i]+=(pr->gr[j]/pr->r[j])*sin(sq->q[i]*pr->r[j]);
            sq->s[i] /= sq->q[i];
        }
    } else {
        for(i=0;i<sq->qn;i++) {
            for(j=0;j<pr->grn;j++)
                sq->s[i]+=(pr->gr[j]/pr->r[j])*sin(sq->q[i]*pr->r[j]);
            sq->s[i] /= sq->q[i];
        }
    }

    return (gmx_sq_t *) sq;
}

int gmx_sans(int argc,char *argv[])
{
    const char *desc[] = {
        "This is simple tool to compute SANS spectra using Debay formula",
        "It currently uses topology file (since it need to assigne element for each atom)",
        "[PAR]",
        "[TT]-pr[TT] Computes normalized g(r) function",
        "[PAR]",
        "[TT]-sq[TT] Computes SANS curve for needed q diapason",
        "[PAR]",
        "[TT]-startq[TT] Starting q value in nm",
        "[PAR]",
        "[TT]-endq[TT] Ending q value in nm",
        "[PAR]",
        "[TT]-qstep[TT] Stepping in q space",
        "[PAR]",
        "Note: When using Debay direct method computational cost increases as",
        "1/2 * N * (N - 1) where N is atom number in group of interest"
    };
    static gmx_bool bPBC=TRUE;
    static real binwidth=0.02,grid=0.05;
    static real start_q=0.0, end_q=2.0, q_step=0.01;
    static gmx_large_int_t  nmc=1048576;
    static unsigned int  seed=0;

    static const char *emode[]= { NULL, "direct", "mc", NULL };
    static const char *emethod[]={ NULL, "debay", "fft", NULL };

    gmx_nstructurefactors    *gnsf;
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
  gmx_bool  bFFT=FALSE, bDEBAY=FALSE;
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
  gmx_gr_t  *pr=NULL;
  gmx_sq_t  *sq=NULL;
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
      bDEBAY=TRUE;
      switch(emode[0][0]) {
      case 'd':
          bDIRECT=TRUE;
          printf("Using Debay direct mode method\n");
          break;
      case 'm':
          bMC=TRUE;
          printf("Using Debay Monte Carlo method\n");
          break;
      default:
          break;
      }
      break;
  case 'f':
      bFFT=TRUE;
      printf("Using FFT method\n");
      break;
  default:
      break;
  }

  if (!bDEBAY && !bFFT)
      gmx_fatal(FARGS,"Unknown method. Set pr or fft!\n");
  if (!bDIRECT && !bMC)
      gmx_fatal(FARGS,"Unknown mode for g(r) method set to direct or mc!");
  /* Try to read files */
  fnDAT = ftp2fn(efDAT,NFILE,fnm);
  fnTPX = ftp2fn(efTPX,NFILE,fnm);

  gnsf = gmx_structurefactors_init(fnDAT);
  printf("Read %d atom names from %s with neutron scattering parameters\n\n",gnsf->nratoms,fnDAT);

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
  if (bDEBAY) {
      if (bDIRECT) {
          /* calc pr */
          pr = calc_pr(gsans,x,index,isize,binwidth,bMC,nmc,seed);
          printf("Were on dabay\n");
      } else if (bMC) {
          if (nmc>(gmx_large_int_t)truncf(0.5*isize*(isize-1))) {
              printf("Number of mc iteration larger then number of pairs in index group. Switching to direct method!\n");
              bMC=FALSE;
              bDIRECT=TRUE;
              pr = calc_pr(gsans,x,index,isize,binwidth,bMC,nmc,seed);
          } else {
              pr = calc_pr(gsans,x,index,isize,binwidth,bMC,nmc,seed);
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
  sq = pr2iq(pr,start_q,end_q,q_step);
  fp = xvgropen(opt2fn_null("-sq",NFILE,fnm),"I(q)","q (nm^-1)","s(q)/s(0)",oenv);
  for(i=0;i<sq->qn;i++) {
      fprintf(fp,"%10.6lf%10.6lf\n",sq->q[i],sq->s[i]);
  }
  fclose(fp);

  sfree(pr);

  thanx(stderr);
  
  return 0;
}
