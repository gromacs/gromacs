/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.3.99_development_20071104
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2006, The GROMACS development team,
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
 
#ifndef _tune_pol_h
#define _tune_pol_h

enum { 
  eelemH,  eelemHe, 
  eelemLi, eelemBe, eelemB,  eelemC,  eelemN, eelemO, eelemF,  eelemNe,
  eelemNa, eelemMg, eelemAl, eelemSi, eelemP, eelemS, eelemCl, eelemAr,
  eelemK,  eelemCa, eelemBr, eelemKr,
  eelemRb, eelemI,  eelemXe, eelemCs,
  eelemNR 
};

enum { 
  emlH,   emlF,    emlCl,   emlBr,  emlI,  
  emlCTE, emlCTR,  emlCBR,  emlCDI, 
  emlNTE, emlNTR2, emlNPI2, emlNDI, 
  emlOTE, emlOTR4, emlOPI2, 
  emlSTE, emlSTR4, emlSPI2,
  emlPTE, emlSi,   emlNR 
};

enum { 
  eqmHF, eqmMP2DZ, eqmMP2TZ, eqmDFT, eqmDFTD, 
  eqmMiller, eqmKang, eqmBosque, eqmSpoel, eqmNR 
};
#define NQUANT (eqmDFTD+1)

/* Global variable! */
extern char *lbasis[eqmNR];
      
#define ecNR  20
#define eatNR 29
#define eatExtra 20

#define nullify(x) memset(&(x),0,sizeof(x))
#define snew(x,n)   { (x) = calloc((n),sizeof(*x));    }
#define srenew(x,n) { if ((x)) { (x) = realloc((x),(n)*sizeof(*x)); } else { snew((x),(n)); } }
#define sfree(x)    { if ((x!=NULL)) { free(x); x = NULL; } }

typedef struct {
  int    nr;
  char   *formula,*molname;
  int    category[ecNR];
  double weight;
  int    nexperiment;
  double *experiment,*error;
  char   **pname,**reference;
  double qm[eqmNR];
  int    frag_comp[eatNR+eatExtra];
  int    elem_comp[eelemNR];
  int    emil_comp[emlNR];
} t_molprop;

typedef struct {
  char   *name,*atom,*miller,*bosque;
  int    nh,maxbond,hybrid;
  double train,all,blen;
} t_spoel;

typedef struct {
  char   *name;
  int    nelec;
  double value,mass;
} t_bosque;

typedef struct {
  char   *name;
  int    nelec;
  double tau_ahc,alpha_ahp;
} t_miller;

/* Global variables! */
extern t_miller miller[emlNR];
extern t_spoel  spoel[eatNR+eatExtra];
extern t_bosque bosque[eelemNR+1];
extern char     *ec_name[ecNR];

extern void write_molprops(char *fn,int npd,t_molprop pd[]);

extern int read_molprops(char *fn,t_molprop **pd,int update_bm);
/* Return number of molprops. If update_bm != 0 then the Bosque and Miller
   compositions will be updated based on the Spoel composition */

extern void gaussj(double **a,int n,double **b,int m);

extern void fatal_error(char *str,char *val);

extern int mp_num_prop(t_molprop *mp,char *prop);
#define mp_num_polar(mp) mp_num_prop(mp,"Polarizability")
#define mp_num_dip(mp)   mp_num_prop(mp,"Dipole")

extern double mp_get_prop(t_molprop *mp,char *prop,int index);
#define mp_get_polar(mp) mp_get_prop(mp,"Polarizability",0)
#define mp_get_dip(mp)   mp_get_prop(mp,"Dipole",0)

extern char *mp_get_ref_prop(t_molprop *mp,char *prop,int index);
#define mp_get_ref_polar(mp) mp_get_ref_prop(mp,"Polarizability",0)
#define mp_get_ref_dip(mp)   mp_get_ref_prop(mp,"Dipole",0)


extern t_molprop *merge_xml(int argc,char *argv[],char *outf,
			    char *sorted,char *doubles,int *nmolprop);

#endif
