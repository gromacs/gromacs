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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "sysstuff.h"
#include "fatal.h"
#include "names.h"
#include "macros.h"
#include "nrnb.h"
#include "main.h"
#include "smalloc.h"
#include "copyrite.h"

typedef struct {
  char *name;
  char *loop;
  int  flop;
} t_nrnb_data;


static const t_nrnb_data nbdata[eNRNB] = {
  { "LJ",                        "inl0100",      31 }, 
  { "LJ(S)",                     "inl0110",      31 },   
  { "Buckingham",                "inl0200",      36 }, 
  { "Buckingham(S)",             "inl0210",      36 },   
  { "LJ(T)",                     "inl0300",      52 }, 
  { "FreeEner LJ(T)",            "inl0301",      65 }, 
  { "Softcore LJ(T)",            "inl0302",     118 }, 
  { "LJ(T)(S)",                  "inl0310",      52 },   
  { "Buckingham(T)",             "inl0400",      57 }, 
  { "FreeEner Bham(T)",          "inl0401",      89 }, 
  { "Softcore Bham(T)",          "inl0402",     128 }, 
  { "Buckingham(T)(S)",          "inl0410",      57 },   
  { "Coulomb",                   "inl1000",      27 }, 
  { "Coulomb(S)",                "inl1010",      27 }, 
  { "Coulomb(W)",                "inl1020",      81 }, 
  { "Coulomb(WW)",               "inl1030",     234 }, 
  { "LJ + Coulomb",              "inl1100",      38 }, 
  { "LJ + Coul(S)",              "inl1110",      38 }, 
  { "LJ + Coul(W)",              "inl1120",      92 }, 
  { "LJ + Coul(WW)",             "inl1130",     245 }, 
  { "Buckingham + Coul",         "inl1200",      40 }, 
  { "Bham + Coul(S)",            "inl1210",      40 }, 
  { "Bham + Coul(W)",            "inl1220",      94 }, 
  { "Bham + Coul(WW)",           "inl1230",     247 }, 
  { "LJ(T) + Coul ",             "inl1300",      58 }, 
  { "LJ(T) + Coul(S)",           "inl1310",      58 }, 
  { "LJ(T) + Coul(W)",           "inl1320",     112 }, 
  { "LJ(T) + Coul(WW)",          "inl1330",     265 }, 
  { "Bham(T) + Coul ",           "inl1400",      63 }, 
  { "Bham(T) + Coul(S)",         "inl1410",      63 }, 
  { "Bham(T) + Coul(W)",         "inl1420",     117 }, 
  { "Bham(T) + Coul(WW)",        "inl1430",     270 }, 
  { "RF Coul",                   "inl2000",      32 }, 
  { "RF Coul(S)",                "inl2010",      32 }, 
  { "RF Coul(W)",                "inl2020",      96 }, 
  { "RF Coul(WW)",               "inl2030",     279 }, 
  { "LJ + RF Coul",              "inl2100",      43 }, 
  { "LJ + RF Coul(S)",           "inl2110",      43 }, 
  { "LJ + RF Coul(W)",           "inl2120",     107 }, 
  { "LJ + RF Coul(WW)",          "inl2130",     290 }, 
  { "Bham + RF Coul",            "inl2200",      45 }, 
  { "Bham + RF Coul(S)",         "inl2210",      45 }, 
  { "Bham + RF Coul(W)",         "inl2220",     109 }, 
  { "Bham + RF Coul(WW)",        "inl2230",     292 }, 
  { "LJ(T) + RF Coul",           "inl2300",      63 }, 
  { "LJ(T) + RF Coul(S)",        "inl2310",      63 }, 
  { "LJ(T) + RF Coul(W)",        "inl2320",     127 }, 
  { "LJ(T) + RF Coul(WW)",       "inl2330",     310 }, 
  { "Bham(T) + RF",              "inl2400",      68 }, 
  { "Bham(T) + RF(S)",           "inl2410",      68 }, 
  { "Bham(T) + RF(W)",           "inl2420",     132 }, 
  { "Bham(T) + RF(WW)",          "inl2430",     310 }, 
  { "Coulomb(T)",                "inl3000",      41 }, 
  { "FreeEner Coul(T)",          "inl3001",      54 }, 
  { "Softcore Coul(T)",          "inl3002",      99 }, 
  { "Coulomb(T)(S)",             "inl3010",      99 }, 
  { "Coulomb(T)(W)",             "inl3020",     123 }, 
  { "Coulomb(T)(WW)",            "inl3030",     360 }, 
  { "LJ + Coulomb(T)",           "inl3100",      54 }, 
  { "LJ + Coulomb(T)(S)",        "inl3110",      54 }, 
  { "LJ + Coulomb(T)(W)",        "inl3120",     136 }, 
  { "LJ + Coulomb(T)(WW)",       "inl3130",     373 }, 
  { "Bham + Coul(T)",            "inl3200",      55 }, 
  { "Bham + Coul(T)(S)",         "inl3210",      55 }, 
  { "Bham + Coul(T)(W)",         "inl3220",     137 }, 
  { "Bham + Coul(T)(WW)",        "inl3230",     374 }, 
  { "LJ(T) + Coul(T)",           "inl3300",      67 }, 
  { "Free LJ(T)+Coul(T)",        "inl3301",      92 }, 
  { "SC LJ(T)+Coul(T)",          "inl3302",     151 }, 
  { "LJ(T) + Coul(T)(S)",        "inl3310",     151 }, 
  { "LJ(T) + Coul(T)(W)",        "inl3320",     149 }, 
  { "LJ(T) + Coul(T)(WW)",       "inl3330",     386 }, 
  { "LJ(T) + Coul(T)",           "inl3400",      71 }, 
  { "Free Bham(T)+Coul(T)",      "inl3401",     116 }, 
  { "SC Bham(T)+Coul(T)",        "inl3402",     161 }, 
  { "Bham(T) + Coul(T)(S)",      "inl3410",     161 }, 
  { "Bham(T) + Coul(T,W)",       "inl3420",     153 }, 
  { "Bham(T) + Coul(T,WW)",      "inl3430",     390 }, 
  { "Innerloop-Iatom",           "-",            10 },
  { "Calc Weights",              "-",            36 },
  { "Spread Q",                  "-",             6 },
  { "Spread Q Bspline",          "-",             2 }, 
  { "Gather F",                  "-",            23 },
  { "Gather F Bspline",          "-",            12 }, 
  { "3D-FFT",                    "-",             8 },
  { "Convolution",               "-",             4 },
  { "Solve PME",                 "-",            64 },
  { "NS-Pairs",                  "-",            21 },
  { "Reset In Box",              "-",             9 },
  { "Shift-X",                   "-",             6 },
  { "CG-CoM",                    "-",            29 },
  { "Sum Forces",                "-",             1 },
  { "Bonds",                     "-",            43 },
  { "G96Bonds",                  "-",            40 },
  { "Angles",                    "-",           163 },
  { "G96Angles",                 "-",           150 },
  { "Propers",                   "-",           229 },
  { "Impropers",                 "-",           208 },
  { "RB-Dihedrals",              "-",           247 },
  { "Four. Dihedrals",           "-",           247 },
  { "Dist. Restr.",              "-",           200 },
  { "Orient. Restr.",            "-",           200 },
  { "Dihedral Restr.",           "-",           200 },
  { "Pos. Restr.",               "-",            50 },
  { "Angle Restr.",              "-",           191 },
  { "Angle Restr. Z",            "-",           164 },
  { "Morse Potent.",             "-",            58 },
  { "Cubic Bonds",               "-",            54 },
  { "Water Pol.",                "-",            62 },
  { "Virial",                    "-",            18 },
  { "Update",                    "-",            31 },
  { "Ext.ens. Update",           "-",            54 },
  { "Stop-CM",                   "-",            10 },
  { "P-Coupling",                "-",             6 },
  { "Calc-Ekin",                 "-",            27 },
  { "Lincs",                     "-",            60 },
  { "Lincs-Mat",                 "-",             4 },
  { "Shake",                     "-",            30 },
  { "Shake-V",                   "-",            15 },
  { "Shake-Init",                "-",            10 },
  { "Shake-Vir",                 "-",            18 },
  { "Settle",                    "-",           323 },
  { "PShake-InitLD",             "-",            59 },    
  { "PShake-InitMD",             "-",            65 },   
  { "PShake",                    "-",             7 },
  { "Dummy2",                    "-",            17 },
  { "Dummy3",                    "-",            28 },
  { "Dummy3fd",                  "-",            95 },
  { "Dummy3fad",                 "-",           176 },
  { "Dummy3out",                 "-",            87 },
  { "Dummy4fd",                  "-",           110 } 
};

void init_nrnb(t_nrnb *nrnb)
{
  int i;

  for(i=0; (i<eNRNB); i++)
    nrnb->n[i]=0.0;
}

void cp_nrnb(t_nrnb *dest, t_nrnb *src)
{
  int i;

  for(i=0; (i<eNRNB); i++)
    dest->n[i]=src->n[i];
}

void add_nrnb(t_nrnb *dest, t_nrnb *s1, t_nrnb *s2)
{
  int i;

  for(i=0; (i<eNRNB); i++)
    dest->n[i]=s1->n[i]+s2->n[i];
}

void print_nrnb(FILE *out, t_nrnb *nrnb)
{
  int i;

  for(i=0; (i<eNRNB); i++)
    fprintf(out,"%14s  %10s  %10.0f.\n",nbdata[i].name,
	    nbdata[i].loop,nrnb->n[i]);
}

void _inc_nrnb(t_nrnb *nrnb,int enr,int inc,char *file,int line)
{
  nrnb->n[enr]+=inc;
#ifdef DEBUG
  printf("nrnb %15s(%2d) incremented with %8d from file %s line %d\n",
	  nbdata[enr].name,enr,inc,file,line);
#endif
}

void print_perf(FILE *out,double nodetime,double realtime,real runtime,
		t_nrnb *nrnb,int nprocs)
{
  int    i;
  double nbfs,mni,frac,tfrac,mflop,tflop;
  char   *myline = "---------------------------------------------------------------------";
  if (nodetime == 0.0) {
    fprintf(out,"nodetime = 0! Infinite Giga flopses!\n");
  }
  
  nbfs=0.0;
  for(i=0; (i<eNR_INLOOP); i++) {
    if (strstr(nbdata[i].name,"(WW)") != NULL)
      nbfs += 9e-6*nrnb->n[i];
    else if (strstr(nbdata[i].name,"(W)") != NULL)
      nbfs += 3e-6*nrnb->n[i];
    else
      nbfs += 1e-6*nrnb->n[i];
  }
  tflop=0;
  for(i=0; (i<eNRNB); i++) 
    tflop+=1e-6*nrnb->n[i]*nbdata[i].flop;
  
  if (tflop == 0) {
    fprintf(out,"No MEGA Flopsen this time\n");
    return;
  }
  fprintf(out,"\n\tM E G A - F L O P S   A C C O U N T I N G\n\n");
  if (nprocs > 1) {
    nodetime = realtime;
    fprintf(out,"\tBased on real time for parallel computer.\n");
  }

  fprintf(out,"   RF=Reaction-field  Free=Free Energy  SC=Softcore\n");
  fprintf(out,"   T=Tabulated        S=Solvent         W=Water     WW=Water-Water\n\n");
  fprintf(out,"%-21s  %8s  %12s  %12s  %8s\n",
	  "Computing:","Function","M-Number","M-Flops","% Flops");
  fprintf(out,"%s\n",myline);
  mflop=0.0;
  tfrac=0.0;
  for(i=0; (i<eNRNB); i++) {
    mni    = 1e-6*nrnb->n[i];
    mflop += mni*nbdata[i].flop;
    frac   = 100.0*mni*nbdata[i].flop/tflop;
    tfrac += frac;
    if (mni != 0)
      fprintf(out,"%-21s  %8s  %12.6f  %12.6f  %6.1f\n",
	      nbdata[i].name,nbdata[i].loop,mni,mni*nbdata[i].flop,frac);
  }
  fprintf(out,"%s\n",myline);
  fprintf(out,"%-21s  %8s  %12s  %12.5f  %6.1f\n",
	  "Total","","",mflop,tfrac);
  fprintf(out,"%s\n\n",myline);
  
  if ((nodetime > 0) && (realtime > 0)) {
    fprintf(out,"%12s %10s %10s %8s\n","","NODE (s)","Real (s)","(%)");
    fprintf(out,"%12s %10.3f %10.3f %8.1f\n","Time:",
	    nodetime, realtime, 100.0*nodetime/realtime);
    if (nodetime > 60) {
      fprintf(out,"%12s %10s","","");
      pr_difftime(out,nodetime);
    }
    if (runtime>0) { /* runtime=0 means calc energies only */
      mflop = mflop/nodetime;
      fprintf(out,"%12s %10s %10s %10s %10s\n",
	      "","(Mnbf/s)",(mflop > 1000) ? "(GFlops)" : "(MFlops)",
	      "(ps/NODE hour)","(NODE hour/ns)");
      fprintf(out,"%12s %10.3f %10.3f %10.3f %10.3f\n","Performance:",
	      nbfs/nodetime,(mflop > 1000) ? (mflop/1000) : mflop,
	      runtime*3600/nodetime,1000*nodetime/(3600*runtime));
    }
  }
}

int cost_nrnb(int enr)
{
  return nbdata[enr].flop;
}

char *nrnb_str(int enr)
{
  return nbdata[enr].name;
}

static const int    force_index[]={ 
  eNR_BONDS,  eNR_ANGLES,  eNR_PROPER, eNR_IMPROPER, 
  eNR_RB,     eNR_DISRES,  eNR_ORIRES, eNR_POSRES,
  eNR_NS,     eNR_INL_IATOM
};
#define NFORCE_INDEX asize(force_index)

static const int    shake_index[]={ 
  eNR_SHAKE,     eNR_SHAKE_RIJ, eNR_SETTLE,       eNR_UPDATE,       eNR_PCOUPL,
  eNR_SHAKE_VIR, eNR_SHAKE_V,   eNR_PSHAKEINITLD, eNR_PSHAKEINITMD, eNR_PSHAKE
};
#define NSHAKE_INDEX asize(shake_index)

static double pr_av(FILE *log,int nprocs,double fav,double ftot[],char *title)
{
  int    i,perc;
  double dperc,unb;
  
  unb=0;
  if (fav > 0) {
    fav/=nprocs;
    fprintf(log,"\n%15s:",title);
    for(i=0; (i<nprocs); i++) {
      dperc=(100.0*ftot[i])/fav;
      unb=max(unb,dperc);
      perc=dperc;
      fprintf(log,"%3d ",perc);
    }
    if (unb > 0) {
      perc=10000.0/unb;
      fprintf(log,"%6d%%\n\n",perc);
    }
    else
      fprintf(log,"\n\n");
  }
  return unb;
}

void pr_load(FILE *log,int nprocs,t_nrnb nrnb[])
{
  int    i,j,perc;
  double dperc,unb,uf,us;
  double *ftot,fav;
  double *stot,sav;
  t_nrnb *av;

  snew(av,1);
  snew(ftot,nprocs);
  snew(stot,nprocs);
  init_nrnb(av);
  for(i=0; (i<nprocs); i++) {
    add_nrnb(av,av,&(nrnb[i]));
    /* Cost due to forces */
    for(j=0; (j<eNR_INLOOP); j++)
      ftot[i]+=nrnb[i].n[j]*cost_nrnb(j);
    for(j=0; (j<NFORCE_INDEX); j++) 
      ftot[i]+=nrnb[i].n[force_index[j]]*cost_nrnb(force_index[j]);
    /* Due to shake */
    for(j=0; (j<NSHAKE_INDEX); j++) {
      stot[i]+=nrnb[i].n[shake_index[j]]*cost_nrnb(shake_index[j]);
    }
  }    
  for(j=0; (j<eNRNB); j++)
    av->n[j]=av->n[j]/(double)nprocs;
  
  fprintf(log,"\nDetailed load balancing info in percentage of average\n");
  
  fprintf(log,"Type        NODE:");
  for(i=0; (i<nprocs); i++)
    fprintf(log,"%3d ",i);
  fprintf(log,"Scaling\n");
  fprintf(log,"----------------");
  for(i=0; (i<nprocs); i++)
    fprintf(log,"----");
  fprintf(log,"-------\n");
  
  for(j=0; (j<eNRNB); j++) {
    unb=100.0;
    if (av->n[j] > 0) {
      fprintf(log,"%15s:",nrnb_str(j));
      for(i=0; (i<nprocs); i++) {
	dperc=(100.0*nrnb[i].n[j])/av->n[j];
	unb=max(unb,dperc);
	perc=dperc;
	fprintf(log,"%3d ",perc);
      }
      if (unb > 0) {
	perc=10000.0/unb;
	fprintf(log,"%6d%%\n",perc);
      }
      else
	fprintf(log,"\n");
    }   
  }
  fav=sav=0;
  for(i=0; (i<nprocs); i++) {
    fav+=ftot[i];
    sav+=stot[i];
  }
  uf=pr_av(log,nprocs,fav,ftot,"Total Force");
  us=pr_av(log,nprocs,sav,stot,"Total Shake");
  
  unb=(uf*fav+us*sav)/(fav+sav);
  if (unb > 0) {
    unb=10000.0/unb;
    fprintf(log,"\nTotal Scaling: %.0f%% of max performance\n\n",unb);
  }
}

