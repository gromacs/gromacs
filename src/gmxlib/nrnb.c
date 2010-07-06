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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <string.h>
#include "sysstuff.h"
#include "gmx_fatal.h"
#include "names.h"
#include "macros.h"
#include "nrnb.h"
#include "main.h"
#include "smalloc.h"
#include "copyrite.h"

typedef struct {
  const char *name;
  int  flop;
} t_nrnb_data;


static const t_nrnb_data nbdata[eNRNB] = {
    { "LJ",                             33 }, /* nb_kernel010 */
    { "Buckingham",                     61 }, /* nb_kernel020 */ 
    { "VdW(T)",                         54 }, /* nb_kernel030 */
    { "Coulomb",                        27 }, /* nb_kernel100 */
    { "Coulomb [W3]",                   80 }, /* nb_kernel101 */
    { "Coulomb [W3-W3]",               234 }, /* nb_kernel102 */
    { "Coulomb [W4]",                   80 }, /* nb_kernel103 */
    { "Coulomb [W4-W4]",               234 }, /* nb_kernel104 */
    { "Coulomb + LJ",                   38 }, /* nb_kernel110 */
    { "Coulomb + LJ [W3]",              91 }, /* nb_kernel111 */
    { "Coulomb + LJ [W3-W3]",          245 }, /* nb_kernel112 */
    { "Coulomb + LJ [W4]",             113 }, /* nb_kernel113 */
    { "Coulomb + LJ [W4-W4]",          267 }, /* nb_kernel114 */
    { "Coulomb + Bham ",                64 }, /* nb_kernel120 */
    { "Coulomb + Bham [W3]",           117 }, /* nb_kernel121 */
    { "Coulomb + Bham [W3-W3]",        271 }, /* nb_kernel122 */
    { "Coulomb + Bham [W4]",           141 }, /* nb_kernel123 */
    { "Coulomb + Bham [W4-W4]",        295 }, /* nb_kernel124 */
    { "Coulomb + VdW(T) ",              59 }, /* nb_kernel130 */
    { "Coulomb + VdW(T) [W3]",         112 }, /* nb_kernel131 */
    { "Coulomb + VdW(T) [W3-W3]",      266 }, /* nb_kernel132 */
    { "Coulomb + VdW(T) [W4]",         134 }, /* nb_kernel133 */
    { "Coulomb + VdW(T) [W4-W4]",      288 }, /* nb_kernel134 */
    { "RF Coul",                        33 }, /* nb_kernel200 */
    { "RF Coul [W3]",                   98 }, /* nb_kernel201 */
    { "RF Coul [W3-W3]",               288 }, /* nb_kernel202 */
    { "RF Coul [W4]",                   98 }, /* nb_kernel203 */
    { "RF Coul [W4-W4]",               288 }, /* nb_kernel204 */
    { "RF Coul + LJ",                   44 }, /* nb_kernel210 */
    { "RF Coul + LJ [W3]",             109 }, /* nb_kernel211 */
    { "RF Coul + LJ [W3-W3]",          299 }, /* nb_kernel212 */
    { "RF Coul + LJ [W4]",             131 }, /* nb_kernel213 */
    { "RF Coul + LJ [W4-W4]",          321 }, /* nb_kernel214 */
    { "RF Coul + Bham ",                70 }, /* nb_kernel220 */
    { "RF Coul + Bham [W3]",           135 }, /* nb_kernel221 */
    { "RF Coul + Bham [W3-W3]",        325 }, /* nb_kernel222 */
    { "RF Coul + Bham [W4]",           159 }, /* nb_kernel223 */
    { "RF Coul + Bham [W4-W4]",        349 }, /* nb_kernel224 */
    { "RF Coul + VdW(T) ",              65 }, /* nb_kernel230 */
    { "RF Coul + VdW(T) [W3]",         130 }, /* nb_kernel231 */
    { "RF Coul + VdW(T) [W3-W3]",      320 }, /* nb_kernel232 */
    { "RF Coul + VdW(T) [W4]",         152 }, /* nb_kernel233 */
    { "RF Coul + VdW(T) [W4-W4]",      342 }, /* nb_kernel234 */
    { "Coul(T)",                        42 }, /* nb_kernel300 */
    { "Coul(T) [W3]",                  125 }, /* nb_kernel301 */
    { "Coul(T) [W3-W3]",               369 }, /* nb_kernel302 */
    { "Coul(T) [W4]",                  125 }, /* nb_kernel303 */
    { "Coul(T) [W4-W4]",               369 }, /* nb_kernel304 */
    { "Coul(T) + LJ",                   55 }, /* nb_kernel310 */
    { "Coul(T) + LJ [W3]",             138 }, /* nb_kernel311 */
    { "Coul(T) + LJ [W3-W3]",          382 }, /* nb_kernel312 */
    { "Coul(T) + LJ [W4]",             158 }, /* nb_kernel313 */
    { "Coul(T) + LJ [W4-W4]",          402 }, /* nb_kernel314 */
    { "Coul(T) + Bham",                 81 }, /* nb_kernel320 */
    { "Coul(T) + Bham [W3]",           164 }, /* nb_kernel321 */
    { "Coul(T) + Bham [W3-W3]",        408 }, /* nb_kernel322 */
    { "Coul(T) + Bham [W4]",           186 }, /* nb_kernel323 */
    { "Coul(T) + Bham [W4-W4]",        430 }, /* nb_kernel324 */
    { "Coul(T) + VdW(T)",               68 }, /* nb_kernel330 */
    { "Coul(T) + VdW(T) [W3]",         151 }, /* nb_kernel331 */
    { "Coul(T) + VdW(T) [W3-W3]",      395 }, /* nb_kernel332 */
    { "Coul(T) + VdW(T) [W4]",         179 }, /* nb_kernel333 */
    { "Coul(T) + VdW(T) [W4-W4]",      423 }, /* nb_kernel334 */
    { "Generalized Born Coulomb",       48 }, /* nb_kernel400 */
    { "GB Coulomb + LJ",                61 }, /* nb_kernel410 */
    { "GB Coulomb + VdW(T)",            79 }, /* nb_kernel430 */
    { "LJ NF",                          19 }, /* nb_kernel010nf */
    { "Buckingham NF",                  48 }, /* nb_kernel020nf */ 
    { "VdW(T) NF",                      33 }, /* nb_kernel030nf */
    { "Coulomb NF",                     16 }, /* nb_kernel100nf */
    { "Coulomb [W3] NF",                47 }, /* nb_kernel101nf */
    { "Coulomb [W3-W3] NF",            135 }, /* nb_kernel102nf */
    { "Coulomb [W4] NF",                47 }, /* nb_kernel103nf */
    { "Coulomb [W4-W4] NF",            135 }, /* nb_kernel104nf */
    { "Coulomb + LJ NF",                24 }, /* nb_kernel110nf */
    { "Coulomb + LJ [W3] NF",           55 }, /* nb_kernel111nf */
    { "Coulomb + LJ [W3-W3] NF",       143 }, /* nb_kernel112nf */
    { "Coulomb + LJ [W4] NF",           66 }, /* nb_kernel113nf */
    { "Coulomb + LJ [W4-W4] NF",       154 }, /* nb_kernel114nf */
    { "Coulomb + Bham  NF",             51 }, /* nb_kernel120nf */
    { "Coulomb + Bham [W3] NF",         82 }, /* nb_kernel121nf */
    { "Coulomb + Bham [W3-W3] NF",     170 }, /* nb_kernel122nf */
    { "Coulomb + Bham [W4] NF",         95 }, /* nb_kernel123nf */
    { "Coulomb + Bham [W4-W4] NF",     183 }, /* nb_kernel124nf */
    { "Coulomb + VdW(T)  NF",           36 }, /* nb_kernel130nf */
    { "Coulomb + VdW(T) [W3] NF",       67 }, /* nb_kernel131nf */
    { "Coulomb + VdW(T) [W3-W3] NF",   155 }, /* nb_kernel132nf */
    { "Coulomb + VdW(T) [W4] NF",       80 }, /* nb_kernel133nf */
    { "Coulomb + VdW(T) [W4-W4] NF",   168 }, /* nb_kernel134nf */
    { "RF Coul NF",                     19 }, /* nb_kernel200nf */
    { "RF Coul [W3] NF",                56 }, /* nb_kernel201nf */
    { "RF Coul [W3-W3] NF",            162 }, /* nb_kernel202nf */
    { "RF Coul [W4] NF",                56 }, /* nb_kernel203nf */
    { "RF Coul [W4-W4] NF",            162 }, /* nb_kernel204nf */
    { "RF Coul + LJ NF",                27 }, /* nb_kernel210nf */
    { "RF Coul + LJ [W3] NF",           64 }, /* nb_kernel211nf */
    { "RF Coul + LJ [W3-W3] NF",       170 }, /* nb_kernel212nf */
    { "RF Coul + LJ [W4] NF",           75 }, /* nb_kernel213nf */
    { "RF Coul + LJ [W4-W4] NF",       181 }, /* nb_kernel214nf */
    { "RF Coul + Bham  NF",             54 }, /* nb_kernel220nf */
    { "RF Coul + Bham [W3] NF",         91 }, /* nb_kernel221nf */
    { "RF Coul + Bham [W3-W3] NF",     197 }, /* nb_kernel222nf */
    { "RF Coul + Bham [W4] NF",        104 }, /* nb_kernel223nf */
    { "RF Coul + Bham [W4-W4] NF",     210 }, /* nb_kernel224nf */
    { "RF Coul + VdW(T)  NF",           39 }, /* nb_kernel230nf */
    { "RF Coul + VdW(T) [W3] NF",       76 }, /* nb_kernel231nf */
    { "RF Coul + VdW(T) [W3-W3] NF",   182 }, /* nb_kernel232nf */
    { "RF Coul + VdW(T) [W4] NF",       89 }, /* nb_kernel233nf */
    { "RF Coul + VdW(T) [W4-W4] NF",   195 }, /* nb_kernel234nf */
    { "Coul(T) NF",                     26 }, /* nb_kernel300nf */
    { "Coul(T) [W3] NF",                77 }, /* nb_kernel301nf */
    { "Coul(T) [W3-W3] NF",            225 }, /* nb_kernel302nf */
    { "Coul(T) [W4] NF",                77 }, /* nb_kernel303nf */
    { "Coul(T) [W4-W4] NF",            225 }, /* nb_kernel304nf */
    { "Coul(T) + LJ NF",                34 }, /* nb_kernel310nf */
    { "Coul(T) + LJ [W3] NF",           85 }, /* nb_kernel311nf */
    { "Coul(T) + LJ [W3-W3] NF",       233 }, /* nb_kernel312nf */
    { "Coul(T) + LJ [W4] NF",           96 }, /* nb_kernel313nf */
    { "Coul(T) + LJ [W4-W4] NF",       244 }, /* nb_kernel314nf */
    { "Coul(T) + Bham NF",              61 }, /* nb_kernel320nf */
    { "Coul(T) + Bham [W3] NF",        112 }, /* nb_kernel321nf */
    { "Coul(T) + Bham [W3-W3] NF",     260 }, /* nb_kernel322nf */
    { "Coul(T) + Bham [W4] NF",        125 }, /* nb_kernel323nf */
    { "Coul(T) + Bham [W4-W4] NF",     273 }, /* nb_kernel324nf */
    { "Coul(T) + VdW(T) NF",            42 }, /* nb_kernel330nf */
    { "Coul(T) + VdW(T) [W3] NF",       93 }, /* nb_kernel331nf */
    { "Coul(T) + VdW(T) [W3-W3] NF",   241 }, /* nb_kernel332nf */
    { "Coul(T) + VdW(T) [W4] NF",      110 }, /* nb_kernel333nf */
    { "Coul(T) + VdW(T) [W4-W4] NF",   258 }, /* nb_kernel334nf */
    { "Generalized Born Coulomb NF",    29 }, /* nb_kernel400nf */
    { "GB Coulomb + LJ NF",             37 }, /* nb_kernel410nf */
    { "GB Coulomb + VdW(T) NF",         49 }, /* nb_kernel430nf */
    { "Free energy innerloop",         150 }, /* free energy, estimate */  
    { "All-vs-All, Coul + LJ",          38 },
    { "All-vs-All, GB + LJ",            61 },
    { "Outer nonbonded loop",           10 },
    { "1,4 nonbonded interactions",     90 },
    { "Born radii (Still)",             47 },
    { "Born radii (HCT/OBC)",          183 },
    { "Born force chain rule",          15 },
    { "All-vs-All Still radii",         47 },
    { "All-vs-All HCT/OBC radii",      183 },
    { "All-vs-All Born chain rule",     15 },
    { "Calc Weights",                   36 },
    { "Spread Q",                        6 },
    { "Spread Q Bspline",                2 }, 
    { "Gather F",                      23  },
    { "Gather F Bspline",              6   }, 
    { "3D-FFT",                        8   },
    { "Convolution",                   4   },
    { "Solve PME",                     64  },
    { "NS-Pairs",                      21  },
    { "Reset In Box",                  3   },
    { "Shift-X",                       6   },
    { "CG-CoM",                        3   },
    { "Sum Forces",                    1   },
    { "Bonds",                         59  },
    { "G96Bonds",                      44  },
    { "FENE Bonds",                    58  },
    { "Tab. Bonds",                    62  },
    { "Restraint Potential",           86  },
    { "Angles",                        168 },
    { "G96Angles",                     150 },
    { "Quartic Angles",                160 },
    { "Tab. Angles",                   169 },
    { "Propers",                       229 },
    { "Impropers",                     208 },
    { "RB-Dihedrals",                  247 },
    { "Four. Dihedrals",               247 },
    { "Tab. Dihedrals",                227 },
    { "Dist. Restr.",                  200 },
    { "Orient. Restr.",                200 },
    { "Dihedral Restr.",               200 },
    { "Pos. Restr.",                   50  },
    { "Angle Restr.",                  191 },
    { "Angle Restr. Z",                164 },
    { "Morse Potent.",                 58  },
    { "Cubic Bonds",                   54  },
    { "Walls",                         31  },
    { "Water Pol.",                    62  },
    { "Thole Pol.",                    296 },
    { "Virial",                        18  },
    { "Update",                        31  },
    { "Ext.ens. Update",               54  },
    { "Stop-CM",                       10  },
    { "P-Coupling",                    6   },
    { "Calc-Ekin",                     27  },
    { "Lincs",                         60  },
    { "Lincs-Mat",                     4   },
    { "Shake",                         30  },
    { "Constraint-V",                   8  },
    { "Shake-Init",                    10  },
    { "Constraint-Vir",                24  },
    { "Settle",                        323 },
    { "Virtual Site 2",                23  },
    { "Virtual Site 3",                37  },
    { "Virtual Site 3fd",              95  },
    { "Virtual Site 3fad",             176 },
    { "Virtual Site 3out",             87  },
    { "Virtual Site 4fd",              110 }, 
    { "Virtual Site 4fdn",             254 }, 
    { "Virtual Site N",                 15 },
    { "Mixed Generalized Born stuff",   10 } 
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
    if (nrnb->n[i] > 0)
      fprintf(out," %-26s %10.0f.\n",nbdata[i].name,nrnb->n[i]);
}

void _inc_nrnb(t_nrnb *nrnb,int enr,int inc,char *file,int line)
{
  nrnb->n[enr]+=inc;
#ifdef DEBUG_NRNB
  printf("nrnb %15s(%2d) incremented with %8d from file %s line %d\n",
	  nbdata[enr].name,enr,inc,file,line);
#endif
}

void print_flop(FILE *out,t_nrnb *nrnb,double *nbfs,double *mflop)
{
  int    i;
  double mni,frac,tfrac,tflop;
  const char   *myline = "-----------------------------------------------------------------------------";
  
  *nbfs = 0.0;
  for(i=0; (i<eNR_NBKERNEL_NR); i++) {
    if (strstr(nbdata[i].name,"W3-W3") != NULL)
      *nbfs += 9e-6*nrnb->n[i];
    else if (strstr(nbdata[i].name,"W3") != NULL)
      *nbfs += 3e-6*nrnb->n[i];
    else if (strstr(nbdata[i].name,"W4-W4") != NULL)
      *nbfs += 10e-6*nrnb->n[i];
    else if (strstr(nbdata[i].name,"W4") != NULL)
      *nbfs += 4e-6*nrnb->n[i];
    else
      *nbfs += 1e-6*nrnb->n[i];
  }
  tflop=0;
  for(i=0; (i<eNRNB); i++) 
    tflop+=1e-6*nrnb->n[i]*nbdata[i].flop;
  
  if (tflop == 0) {
    fprintf(out,"No MEGA Flopsen this time\n");
    return;
  }
  if (out) {
    fprintf(out,"\n\tM E G A - F L O P S   A C C O U N T I N G\n\n");
  }

  if (out) {
    fprintf(out,"   RF=Reaction-Field  FE=Free Energy  SCFE=Soft-Core/Free Energy\n");
    fprintf(out,"   T=Tabulated        W3=SPC/TIP3p    W4=TIP4p (single or pairs)\n");
    fprintf(out,"   NF=No Forces\n\n");
    
    fprintf(out," %-32s %16s %15s  %7s\n",
	    "Computing:","M-Number","M-Flops","% Flops");
    fprintf(out,"%s\n",myline);
  }
  *mflop=0.0;
  tfrac=0.0;
  for(i=0; (i<eNRNB); i++) {
    mni     = 1e-6*nrnb->n[i];
    *mflop += mni*nbdata[i].flop;
    frac    = 100.0*mni*nbdata[i].flop/tflop;
    tfrac  += frac;
    if (out && mni != 0)
      fprintf(out," %-32s %16.6f %15.3f  %6.1f\n",
	      nbdata[i].name,mni,mni*nbdata[i].flop,frac);
  }
  if (out) {
    fprintf(out,"%s\n",myline);
    fprintf(out," %-32s %16s %15.3f  %6.1f\n",
	    "Total","",*mflop,tfrac);
    fprintf(out,"%s\n\n",myline);
  }
}

void print_perf(FILE *out,double nodetime,double realtime,int nprocs,
		gmx_large_int_t nsteps,real delta_t,
		double nbfs,double mflop)
{
  real runtime;

  fprintf(out,"\n");

  if (nodetime == 0.0) {
    fprintf(out,"nodetime = 0! Infinite Giga flopses!\n");
  }
#ifdef GMX_OPENMM
  nodetime = realtime;
  fprintf(out,"\tOpenMM run - timing based on wallclock.\n\n");
#else
  if (nprocs > 1)
  {
      nodetime = realtime;
      fprintf(out,"\tParallel run - timing based on wallclock.\n\n");
  }
#endif

  if ((nodetime > 0) && (realtime > 0)) {
    fprintf(out,"%12s %10s %10s %8s\n","","NODE (s)","Real (s)","(%)");
    fprintf(out,"%12s %10.3f %10.3f %8.1f\n","Time:",
	    nodetime, realtime, 100.0*nodetime/realtime);
    if (nodetime > 60) {
      fprintf(out,"%12s %10s","","");
      pr_difftime(out,nodetime);
    }
    if (delta_t > 0) {
      mflop = mflop/nodetime;
      runtime = nsteps*delta_t;
      fprintf(out,"%12s %10s %10s %10s %10s\n",
	      "","(Mnbf/s)",(mflop > 1000) ? "(GFlops)" : "(MFlops)",
	      "(ns/day)","(hour/ns)");
      fprintf(out,"%12s %10.3f %10.3f %10.3f %10.3f\n","Performance:",
	      nbfs/nodetime,(mflop > 1000) ? (mflop/1000) : mflop,
	      runtime*24*3.6/nodetime,1000*nodetime/(3600*runtime));
    } else {
      fprintf(out,"%12s %10s %10s %14s\n",
	      "","(Mnbf/s)",(mflop > 1000) ? "(GFlops)" : "(MFlops)",
	      "(steps/hour)");
      fprintf(out,"%12s %10.3f %10.3f %14.1f\n","Performance:",
	      nbfs/nodetime,(mflop > 1000) ? (mflop/1000) : mflop,
	      nsteps*3600.0/nodetime);
    }
  }
}

int cost_nrnb(int enr)
{
  return nbdata[enr].flop;
}

const char *nrnb_str(int enr)
{
  return nbdata[enr].name;
}

static const int    force_index[]={ 
  eNR_BONDS,  eNR_ANGLES,  eNR_PROPER, eNR_IMPROPER, 
  eNR_RB,     eNR_DISRES,  eNR_ORIRES, eNR_POSRES,
  eNR_NS,     eNR_NBKERNEL_OUTER
};
#define NFORCE_INDEX asize(force_index)

static const int    constr_index[]={ 
  eNR_SHAKE,     eNR_SHAKE_RIJ, eNR_SETTLE,       eNR_UPDATE,       eNR_PCOUPL,
  eNR_CONSTR_VIR,eNR_CONSTR_V
};
#define NCONSTR_INDEX asize(constr_index)

static double pr_av(FILE *log,t_commrec *cr,
		    double fav,double ftot[],const char *title)
{
  int    i,perc;
  double dperc,unb;
  
  unb=0;
  if (fav > 0) {
    fav /= cr->nnodes - cr->npmenodes;
    fprintf(log,"\n %-26s",title);
    for(i=0; (i<cr->nnodes); i++) {
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

void pr_load(FILE *log,t_commrec *cr,t_nrnb nrnb[])
{
  int    i,j,perc;
  double dperc,unb,uf,us;
  double *ftot,fav;
  double *stot,sav;
  t_nrnb *av;

  snew(av,1);
  snew(ftot,cr->nnodes);
  snew(stot,cr->nnodes);
  init_nrnb(av);
  for(i=0; (i<cr->nnodes); i++) {
      add_nrnb(av,av,&(nrnb[i]));
      /* Cost due to forces */
      for(j=0; (j<eNR_NBKERNEL_NR); j++)
	ftot[i]+=nrnb[i].n[j]*cost_nrnb(j);
      for(j=0; (j<NFORCE_INDEX); j++) 
	ftot[i]+=nrnb[i].n[force_index[j]]*cost_nrnb(force_index[j]);
      /* Due to shake */
      for(j=0; (j<NCONSTR_INDEX); j++) {
	stot[i]+=nrnb[i].n[constr_index[j]]*cost_nrnb(constr_index[j]);
      }
  }   
  for(j=0; (j<eNRNB); j++)
    av->n[j]=av->n[j]/(double)(cr->nnodes - cr->npmenodes);
    
    fprintf(log,"\nDetailed load balancing info in percentage of average\n");
  
  fprintf(log," Type                 NODE:");
  for(i=0; (i<cr->nnodes); i++)
      fprintf(log,"%3d ",i);
  fprintf(log,"Scaling\n");
  fprintf(log,"---------------------------");
  for(i=0; (i<cr->nnodes); i++)
      fprintf(log,"----");
  fprintf(log,"-------\n");
  
  for(j=0; (j<eNRNB); j++) {
    unb=100.0;
    if (av->n[j] > 0) {
      fprintf(log," %-26s",nrnb_str(j));
      for(i=0; (i<cr->nnodes); i++) {
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
  for(i=0; (i<cr->nnodes); i++) {
    fav+=ftot[i];
    sav+=stot[i];
  }
  uf=pr_av(log,cr,fav,ftot,"Total Force");
  us=pr_av(log,cr,sav,stot,"Total Constr.");
  
  unb=(uf*fav+us*sav)/(fav+sav);
  if (unb > 0) {
    unb=10000.0/unb;
    fprintf(log,"\nTotal Scaling: %.0f%% of max performance\n\n",unb);
  }
}

