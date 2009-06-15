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

#include "macros.h"
#include "statutil.h"
#include "smalloc.h"
#include "copyrite.h"
#include "gstat.h"
#include "physics.h"
#include "vec.h"
#include "xvgr.h"
#include "pbc.h"
#include "index.h"
#include "tpxio.h"

#include <trajana.h>

typedef struct {
  bool    bMiddle12, bMiddle, bVec23;
  real    invbw, invabw, invrbw;
  real    rcut2, rmin2, rmax2;
  int     *histn;
  int     *hist11, *hist12, *hist21, *hist22;
  int     *hist1a, *hist2a;
  real    *histi11,*histi12,*histi21,*histi22;
  real    *histi1a, *histi2a;
  int     ntot;
  double  sum11,sum12,sum21,sum22;
  double  sum1a, sum2a;
} t_orientdata;

static int
analyze_frame(t_topology *top, t_trxframe *fr, t_pbc *pbc,
              int nr, gmx_ana_selection_t *sel[], void *data)
{
  t_orientdata *d = (t_orientdata *)data;
  int     p, m, n;
  rvec    dxh1, dxh2, dx, outer;
  real    r2, r;
  real    inp, outp, inp2, outp2;
  real    ang1, ang2;

  n = 0;
  for (p = 0; p < sel[0]->p.nr; ++p) {
    for (m = 0; m < sel[1]->p.nr; m += 3) {
      if (d->bMiddle12) {
        rvec_add(sel[1]->p.x[m], sel[1]->p.x[m+1], dxh1);
        svmul(0.5, dxh1, dxh1);
        pbc_dx(pbc, dxh1, sel[0]->p.x[p], dx);
      } else
        pbc_dx(pbc, sel[1]->p.x[m], sel[0]->p.x[p], dx);
      r2  = norm2(dx);
      if (r2 < d->rcut2) {
        r = sqrt(r2);
        rvec_sub(sel[1]->p.x[m+1], sel[1]->p.x[m], dxh1);
        rvec_sub(sel[1]->p.x[m+2], sel[1]->p.x[m], dxh2);
        if (d->bMiddle)
          rvec_inc(dxh1, dxh2);
        svmul(1/r, dx, dx);
        if (!d->bVec23) {
          /* Determine the normal to the plane */
          cprod(dxh1, dxh2, outer);
        } else {
          /* Use the vector between the 2nd and 3rd atom */
          rvec_sub(sel[1]->p.x[m+2], sel[1]->p.x[m+1], outer);
        }
        unitv(dxh1, dxh1);
        inp = iprod(dx, dxh1);
        inp2 = 3*sqr(inp) - 1;
        ang1 = acos(inp) * RAD2DEG;
        /* In some extremely rare cases, rounding errors can cause problems */
        if (isnan(ang1)) { if (inp > 0) ang1 = 0; else ang1 = 180; }
        unitv(outer, outer);
        outp = iprod(dx, outer);
        outp2 = 3*sqr(outp) - 1;
        ang2 = acos(outp) * RAD2DEG;
        if (isnan(ang2)) { if (outp > 0) ang2 = 0; else ang2 = 180; }
        d->histi11[(int)(d->invrbw*r)] += inp;
        d->histi12[(int)(d->invrbw*r)] += inp2;
        d->histi21[(int)(d->invrbw*r)] += outp;
        d->histi22[(int)(d->invrbw*r)] += outp2;
        d->histi1a[(int)(d->invrbw*r)] += ang1;
        d->histi2a[(int)(d->invrbw*r)] += ang2;
        d->histn[(int)(d->invrbw*r)]++;
        if (r2 >= d->rmin2 && r2 < d->rmax2) {
          d->hist11[(int)(d->invbw*(inp + 1))]++;
          d->hist12[(int)(d->invbw*(inp2 + 1))]++;
          d->hist21[(int)(d->invbw*(outp + 1))]++;
          d->hist22[(int)(d->invbw*(outp2 + 1))]++;
          d->hist1a[(int)(d->invabw*ang1)]++;
          d->hist2a[(int)(d->invabw*ang2)]++;
          d->sum11 += inp;
          d->sum12 += inp2;
          d->sum21 += outp;
          d->sum22 += outp2;
          d->sum1a += ang1;
          d->sum2a += ang2;
          n++;
        }
      }
    }
  }
  d->ntot += n;
  return 0;
}

int gmx_sorient(int argc,char *argv[])
{
  FILE    *fp;
  int     i,nbin1,nbin2, nabin, nrbin;
  real    two_pi,nav,normfac,rcut;
  real    c11,c12,c21,c22;
  char    str[STRLEN];
  char    *lege[]  = { "<cos(\\8q\\4\\s1\\N)>",
                       "<cos(\\8q\\4\\s2\\N)>",
                       "<3cos\\S2\\N(\\8q\\4\\s1\\N)-1>",
                       "<3cos\\S2\\N(\\8q\\4\\s2\\N)-1>" };
  char    *legne[] = { "cos(\\8q\\4\\s1\\N)",
                       "cos(\\8q\\4\\s2\\N)",
                       "3cos\\S2\\N(\\8q\\4\\s1\\N)-1",
                       "3cos\\S2\\N(\\8q\\4\\s2\\N)-1" };

  const char *desc[] = {
    "g_sorient analyzes molecule orientation around a group of atoms.",
    "It calculates two angles between the vector from one or more",
    "reference positions to the first atom of each molecule (default) or",
    "to the midpoint of atoms 1 and 2 (with [TT]-m12[tt]):[BR]"
    "theta1: the angle with the vector from the first atom of the",
    "molecule to the midpoint between atoms 2 and 3 (default) or",
    "to atom 2 (with [TT]-nom23[tt]).[BR]",
    "theta2: the angle with the normal of the molecular plane, defined by the",
    "same three atoms, or when the option [TT]-v23[tt] is set",
    "the angle with the vector between atoms 2 and 3.[BR]",
    "The reference can be any set of positions.",
    "The group of molecule atoms should consist of 3 atoms per molecule.",
    "Only molecules between [TT]-rmin[tt] and [TT]-rmax[tt] are",
    "considered for [TT]-o[tt] and [TT]-no[tt] each frame.[PAR]",
    "[TT]-o[tt]: distribtion of cos(theta) and 3cos^2(theta)+1 for",
    "both angles for rmin<=r<=rmax.[PAR]",
    "[TT]-ro[tt]: <cos(theta)> and <3cos^2(theta)-1> for both angles",
    "as a function of the distance.[PAR]",
    "[TT]-co[tt]: the sum over all molecules within distance r",
    "of cos(theta) and 3cos^2(theta)-1 for both angles as a function of r.[PAR]",
    "[TT]-rc[tt]: the distribution of the molecules as a function of r",
  };
  
  static bool bVec23=FALSE,bMiddle12=FALSE,bMiddle=TRUE;
  static real rmin=0.0,rmax=0.5,binwidth=0.02,abinw=1,rbinw=0.02;
  static t_pargs pa[] = {
    { "-m12", FALSE, etBOOL,{&bMiddle12},
      "Use the midpoint of atoms 1 and 2" },
    { "-m23", FALSE, etBOOL,{&bMiddle},
      "Use the midpoint of atoms 2 and 3" },
    { "-v23",  FALSE, etBOOL,  {&bVec23},
      "Use the vector between atoms 2 and 3" },
    { "-rmin",  FALSE, etREAL, {&rmin}, "Minimum distance (nm)" },
    { "-rmax",  FALSE, etREAL, {&rmax}, "Maximum distance (nm)" },
    { "-cbin",  FALSE, etREAL, {&binwidth}, "Binwidth for the cosine" },
    { "-abin",  FALSE, etREAL, {&abinw}, "Binwidth for angles" },
    { "-rbin",  FALSE, etREAL, {&rbinw}, "Binwidth for r (nm)" },
  };
  
  static t_filenm fnm[] = {
    { efXVG, "-o",  "sori.xvg",  ffOPTWR },
    { efXVG, "-ro", "sord.xvg",  ffOPTWR },
    { efXVG, "-oa", "sang.xvg",  ffOPTWR },
    { efXVG, "-roa","sand.xvg",  ffOPTWR },
    { efXVG, "-co", "scum.xvg",  ffOPTWR },
    { efXVG, "-rc", "scount.xvg",  ffOPTWR }
  };
#define NFILE asize(fnm)

  gmx_ana_traj_t       *trj;
  gmx_ana_selection_t **sel;
  t_trxframe           *fr;
  int                   nframes;
  t_orientdata          d;

  CopyRight(stderr,argv[0]);
  gmx_ana_traj_create(&trj, ANA_REQUIRE_WHOLE | ANA_USE_FULLGRPS);
  gmx_ana_set_nrefgrps(trj, 1);
  gmx_ana_set_nanagrps(trj, 1);
  parse_trjana_args(trj, &argc, argv, PCA_CAN_VIEW,
		    NFILE,fnm,asize(pa),pa,asize(desc),desc,0,NULL);
  gmx_ana_get_anagrps(trj, &sel);

  d.bMiddle12 = bMiddle12;
  d.bMiddle   = bMiddle;
  d.bVec23    = bVec23;
  
  two_pi = 2/M_PI;

  if (sel[1]->p.nr % 3 != 0)
    gmx_fatal(FARGS,"The number of atoms (%d) in the analysis group is not a multiple of 3",
              sel[1]->p.nr);

  /* initialize reading trajectory */
  gmx_ana_get_first_frame(trj, &fr);

  d.rmin2 = sqr(rmin);
  d.rmax2 = sqr(rmax);
  rcut  = 0.99*sqrt(max_cutoff2(fr->ePBC,fr->box));
  if (rcut == 0)
    rcut = 10*rmax;
  d.rcut2 = sqr(rcut);

  d.invbw = 1/binwidth;
  nbin1 = (int)(2*d.invbw + 0.5);
  nbin2 = (int)(3*d.invbw + 0.5);

  d.invabw = 1/abinw;
  nabin = (int)(180*d.invabw + 0.5);

  d.invrbw = 1/rbinw;
  nrbin = rcut/rbinw;
  if (nrbin == 0)
    nrbin = 1;

  snew(d.hist11,nbin1+1);
  snew(d.hist12,nbin2+1);
  snew(d.hist21,nbin1+1);
  snew(d.hist22,nbin2+1);
  snew(d.hist1a, nabin+1);
  snew(d.hist2a, nabin+1);

  snew(d.histi11,nrbin);
  snew(d.histi12,nrbin);
  snew(d.histi21,nrbin);
  snew(d.histi22,nrbin);
  snew(d.histi1a,nrbin);
  snew(d.histi2a,nrbin);
  snew(d.histn,nrbin);

  d.ntot = 0;
  d.sum11 = 0;
  d.sum12 = 0;
  d.sum21 = 0;
  d.sum22 = 0;
  d.sum1a = d.sum2a = 0;

  /* Do the actual analysis */
  gmx_ana_do(trj, 0, &analyze_frame, &d);
  gmx_ana_get_nframes(trj, &nframes);

  /* Add the bin for the exact maximum to the previous bin */
  d.hist11[nbin1-1] += d.hist11[nbin1];
  d.hist12[nbin2-1] += d.hist12[nbin2];
  d.hist21[nbin1-1] += d.hist21[nbin1];
  d.hist22[nbin2-1] += d.hist22[nbin2];
  d.hist1a[nabin-1] += d.hist1a[nabin];
  d.hist2a[nabin-1] += d.hist2a[nabin];
  
  nav     = (real)d.ntot/(sel[0]->p.nr*nframes);
  
  fprintf(stderr, "Average nr of molecules between %g and %g nm: %.1f\n",
	  rmin,rmax,nav);
  if (d.ntot > 0) {
    d.sum11 /= d.ntot;
    d.sum12 /= d.ntot;
    d.sum21 /= d.ntot;
    d.sum22 /= d.ntot;
    d.sum1a /= d.ntot; d.sum2a /= d.ntot;
    fprintf(stderr,"Average %-15s between %g and %g nm: %6.3f\n",
            "theta1",rmin,rmax,d.sum1a);
    fprintf(stderr,"Average %-15s between %g and %g nm: %6.3f\n",
            "theta2",rmin,rmax,d.sum2a);
    fprintf(stderr,"Average %-15s between %g and %g nm: %6.3f\n",
            "cos(theta1)",rmin,rmax,d.sum11);
    fprintf(stderr,"Average %-15s between %g and %g nm: %6.3f\n",
            "cos(theta2)",rmin,rmax,d.sum21);
    fprintf(stderr,"Average %-15s between %g and %g nm: %6.3f\n",
            "3cos2(theta1)-1",rmin,rmax,d.sum12);
    fprintf(stderr,"Average %-15s between %g and %g nm: %6.3f\n",
            "3cos2(theta2)-1",rmin,rmax,d.sum22);
  }
  
  if (opt2bSet("-o",NFILE,fnm)) {
    normfac = d.ntot > 0 ? d.invbw / d.ntot : 1;
    sprintf(str,"Molecular orientation between %g and %g nm",rmin,rmax);
    fp=xvgropen(opt2fn("-o",NFILE,fnm),
                str,"","");
    xvgr_selections(fp, trj);
    if (bPrintXvgrCodes())
      fprintf(fp,"@ subtitle \"average shell size %.1f molecules\"\n",nav);
    xvgr_legend(fp,4,legne);
    for(i=0; i<nbin2; i++) {
      fprintf(fp,"%g %g %g %g %g\n",(i+0.5)*binwidth-1,
              i<nbin1?normfac*d.hist11[i]:0.0,
              i<nbin1?normfac*d.hist21[i]:0.0,
              normfac*d.hist12[i],
              normfac*d.hist22[i]);
    }
    fclose(fp);
  }
  
  if (opt2bSet("-ro",NFILE,fnm)) {
    sprintf(str,"Molecular orientation");
    fp=xvgropen(opt2fn("-ro",NFILE,fnm),str,"r (nm)","");
    xvgr_selections(fp, trj);
    if (bPrintXvgrCodes())
      fprintf(fp,"@ subtitle \"as a function of distance\"\n");
    xvgr_legend(fp,4,lege);
    for(i=0; i<nrbin; i++)
      fprintf(fp,"%g %g %g %g %g\n",(i+0.5)*rbinw,
              d.histn[i] ? d.histi11[i]/d.histn[i] : 0,
              d.histn[i] ? d.histi21[i]/d.histn[i] : 0,
              d.histn[i] ? d.histi12[i]/d.histn[i] : 0,
              d.histn[i] ? d.histi22[i]/d.histn[i] : 0);
    fclose(fp);
  }
  
  if (opt2bSet("-oa",NFILE,fnm)) {
    normfac = d.ntot > 0 ? d.invabw / d.ntot : 1;
    sprintf(str,"Molecular orientation between %g and %g nm",rmin,rmax);
    fp=xvgropen(opt2fn("-oa",NFILE,fnm),
                str,"Angle (degrees)","");
    xvgr_selections(fp, trj);
    if (bPrintXvgrCodes())
      fprintf(fp,"@ subtitle \"average shell size %.1f molecules\"\n",nav);
    xvgr_legend(fp,2,legne);
    for(i=0; i<nabin; i++) {
      fprintf(fp,"%g %g %g\n",(i+0.5)*abinw,
              normfac*d.hist1a[i],
              normfac*d.hist2a[i]);
    }
    fclose(fp);
  }

  if (opt2bSet("-roa",NFILE,fnm)) {
    sprintf(str,"Molecular orientation");
    fp=xvgropen(opt2fn("-roa",NFILE,fnm),str,"r (nm)","Average angle (degrees)");
    xvgr_selections(fp, trj);
    if (bPrintXvgrCodes())
      fprintf(fp,"@ subtitle \"as a function of distance\"\n");
    xvgr_legend(fp,2,lege);
    for(i=0; i<nrbin; i++)
      fprintf(fp,"%g %g %g\n",(i+0.5)*rbinw,
              d.histn[i] ? d.histi1a[i]/d.histn[i] : 0,
              d.histn[i] ? d.histi2a[i]/d.histn[i] : 0);
    fclose(fp);
  }

  if (opt2bSet("-co",NFILE,fnm)) {
    sprintf(str,"Cumulative molecular orientation");
    fp=xvgropen(opt2fn("-co",NFILE,fnm),str,"r (nm)","");
    xvgr_selections(fp, trj);
    if (bPrintXvgrCodes())
      fprintf(fp,"@ subtitle \"as a function of distance\"\n");
    xvgr_legend(fp,4,legne);
    normfac = 1.0/(sel[0]->p.nr*nframes);
    c11 = 0;
    c12 = 0;
    c21 = 0;
    c22 = 0;
    fprintf(fp,"%g %g %g %g %g\n",0.0,c11,c21,c12,c22);
    for(i=0; i<nrbin; i++) {
      c11 += d.histi11[i]*normfac;
      c12 += d.histi12[i]*normfac;
      c21 += d.histi21[i]*normfac;
      c22 += d.histi22[i]*normfac;
      fprintf(fp,"%g %g %g %g %g\n",(i+1)*rbinw,c11,c21,c12,c22);
    }
    fclose(fp);
  }

  if (opt2bSet("-rc",NFILE,fnm)) {
    sprintf(str,"Molecular distribution");
    fp=xvgropen(opt2fn("-rc",NFILE,fnm),str,"r (nm)","molecules/nm");
    xvgr_selections(fp, trj);
    if (bPrintXvgrCodes())
      fprintf(fp,"@ subtitle \"as a function of distance\"\n");
    normfac = 1.0/(rbinw*nframes);
    for(i=0; i<nrbin; i++) {
      fprintf(fp,"%g %g\n",(i+0.5)*rbinw,normfac*d.histn[i]);
    }
    fclose(fp);
  }

  do_view(opt2fn("-o",NFILE,fnm),"-nxy");
  do_view(opt2fn("-ro",NFILE,fnm),"-nxy");
  do_view(opt2fn("-co",NFILE,fnm),"-nxy");

  thanx(stderr);
  
  return 0;
}
