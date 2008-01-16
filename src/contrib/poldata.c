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
 
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include "string2.h"
#include "tune_pol.h"

/* Noble gases and cations from G. D. Mahan. Modied Sternheimer
   equation for polarizability. Phys. Rev. A, 22: 1780-1785, 1980. */
   
/* Anions from Christof H{\"a}ttig and Bernd Artur Hess. {TDMP2}
   calculation of dynamic multipole polarizabilities and dispersion
   coefficients for the halogen anions F-, Cl-, Br- and
   I-. J. Chem. Phys., 108:3863-3870, 1998. */
   
t_spoel spoel[eatNR+eatExtra] = {
  { "C33",   "CH3",    "CTE",   "C",  3,  4, 3, 0,   0, 0.15 },
  { "C32",   "CH2",    "CTE",   "C",  2,  4, 3, 0,   0, 0.15 },
  { "C31",   "CH",     "CTE",   "C",  1,  4, 3, 0,   0, 0.15 },
  { "C30",   "C",      "CTE",   "C",  0,  4, 3, 0,   0, 0.15 },
  { "C22",   "CH2",    "CTR",   "C",  2,  3, 2, 0,   0, 0.14 },
  { "C21",   "CH",     "CTR",   "C",  1,  3, 2, 0,   0, 0.14 },
  { "C20",   "C",      "CTR",   "C",  0,  3, 2, 0,   0, 0.14 },
  { "C11",   "CH",     "CDI",   "C",  1,  2, 1, 0,   0, 0.13 },
  { "C10",   "C",      "CDI",   "C",  0,  2, 1, 0,   0, 0.13 },
  { "N32",   "NH2",    "NTE",   "N",  2,  3, 3, 0,   0, 0.13 },
  { "N31",   "NH",     "NTE",   "N",  1,  3, 3, 0,   0, 0.13 },
  { "N30",   "N",      "NTE",   "N",  0,  3, 3, 0,   0, 0.13 },
  { "N22",   "NH2",    "NPI2",  "N",  2,  2, 2, 0,   0, 0.12 },
  { "N21",   "NH",     "NPI2",  "N",  1,  2, 2, 0,   0, 0.12 },
  { "N20",   "N",      "NTR2",  "N",  0,  2, 2, 0,   0, 0.12 },
  { "N10",   "N",      "NDI",   "N",  0,  1, 1, 0,   0, 0.12 },
  { "O31",   "OH",     "OTE",   "O",  1,  2, 3, 0,   0, 0.14 },
  { "O30",   "O",      "OTE",   "O",  0,  2, 3, 0,   0, 0.14 },
  { "O20",   "O",      "OTR4",  "O",  0,  1, 2, 0,   0, 0.13 },
  { "F30",   "F",      "F",     "F",  0,  1, 3, 0,   0, 0.12 },
  { "Si30",  "Si",     "Si",    "Si", 0,  4, 3, 0,   0, 0.16 },
  { "P21",   "PH",     "PTE",   "P",  1,  4, 3, 0,   0, 0.17 },
  { "P20",   "P",      "PTE",   "P",  0,  4, 3, 0,   0, 0.17 },
  { "S31",   "SH",     "STE",   "S",  1,  4, 3, 0,   0, 0.19 },
  { "S30",   "S",      "STE",   "S",  0,  4, 3, 0,   0, 0.19 },
  { "S20",   "S",      "STR4",  "S",  0,  4, 2, 0,   0, 0.18 },
  { "Cl30",  "Cl",     "Cl",    "Cl", 0,  1, 3, 0,   0, 0.20 },
  { "Br30",  "Br",     "Br",    "Br", 0,  1, 3, 0,   0, 0.22 },
  { "I30",   "I",      "I",     "I",  0,  1, 3, 0,   0, 0.24 },
  { "H",     "H",      "H",     "H",  0,  1, 3, 0,   0, 0.108 },
  { "He",    "He",     "He",    "He", 0,  0, 0, 0,   0.24, 0 },
  { "Li+",   "Li+",    "Li",    "Li", 0,  0, 0, 0,   0.032, 0 },
  { "Be2+",  "Be2+",   "Be",    "Be", 0,  0, 0, 0,   0.0083, 0 },
  { "B",     "B",      "B",     "B",  0,  0, 0, 0,   0, 0 },
  { "F-",    "F-",     "F",     "F",  0,  0, 0, 0,   2.467, 0 },
  { "Ne",    "Ne",     "Ne",    "Ne", 0,  0, 0, 0,   0.44, 0 },
  { "Na+",   "Na+",    "Na",    "Na", 0,  0, 0, 0,   0.157, 0 },
  { "Mg2+",  "Mg2+",   "Mg",    "Mg", 0,  0, 0, 0,   0.075, 0 },
  { "Al3+",  "Al3+",   "Al",    "Al", 0,  0, 0, 0,   0, 0 },
  { "Cl-",   "Cl-",    "Cl",    "Cl", 0,  0, 0, 0,   5.482, 0 },
  { "Ar",    "Ar",     "Ar",    "Ar", 0,  0, 0, 0,   1.73, 0 },
  { "K+",    "K+",     "K",     "K",  0,  0, 0, 0,   0.83, 0 },
  { "Ca2+",  "Ca2+",   "Ca",    "Ca", 0,  0, 0, 0,   0.49, 0 },
  { "Br-",   "Br-",    "Br",    "Br", 0,  0, 0, 0,   7.268, 0 },
  { "Kr",    "Kr",     "Kr",    "Kr", 0,  0, 0, 0,   2.58, 0 },
  { "Rb+",   "Rb+",    "Rb",    "Rb", 0,  0, 0, 0,   1.37, 0 },
  { "I-",    "I-",     "I",     "I",  0,  0, 0, 0,   10.275, 0 },
  { "Xe",    "Xe",     "Xe",    "Xe", 0,  0, 0, 0,   4.15, 0 },
  { "Cs+",   "Cs+",    "Cs",    "Cs", 0,  0, 0, 0,   2.36, 0 }
};

t_bosque bosque[eelemNR+1] = {
  { "H",   1, 0.17, 1.0079 }, 
  { "He",  2, 0.0,  4.0026 },
  { "Li",  3, 0.0,  6.941 },
  { "Be",  4, 0.0,  9.01218 },
  { "B",   5, 0.0,  10.81 },
  { "C",   6, 1.51, 12.011 }, 
  { "N",   7, 1.05, 14.00674 }, 
  { "O",   8, 0.57, 15.9994 }, 
  { "F",   9, 0.22, 18.99984 },
  { "Ne", 10, 0.0,  20.179 },
  { "Na", 11, 0.0,  22.98977 },
  { "Mg", 12, 0.0,  24.305 },
  { "Al", 13, 0.0,  26.98154 },
  { "Si", 14, 1.97, 28.086 },
  { "P",  15, 2.48, 30.97376 },
  { "S",  16, 2.99, 32.06 }, 
  { "Cl", 17, 2.16, 35.453 },
  { "Ar", 18, 0.0,  39.948 },
  { "K",  19, 0.0,  39.098 },
  { "Ca", 20, 0.0,  40.08 },
  { "Br", 35, 3.29, 79.904 },
  { "Kr", 36, 0.0,  83.8 }, 
  { "Rb", 37, 0.0,  85.478 },
  { "I",  53, 5.45, 126.9045 },
  { "Xe", 54, 0.0,  131.3 }, 
  { "Cs", 55, 0.0,  132.9054 },
  { "0",   0, 0.32, 0 }
};

t_miller miller[emlNR] = {
  { "H",      1, 0.313, 0.387 },  
  { "F",      9, 1.089, 0.296 },  
  { "Cl",    17, 3.165, 2.315 },  
  { "Br",    35, 5.566, 3.013 },  
  { "I",     53, 8.593, 5.415 },  
  { "CTE",    6, 1.294, 1.061 },  
  { "CTR",    6, 1.433, 1.352 },  
  { "CBR",    6, 1.707, 1.896 },  
  { "CDI",    6, 1.393, 1.283 },  
  { "NTE",    7, 1.373, 0.964 },  
  { "NTR2",   7, 1.262, 1.030 },  
  { "NPI2",   7, 1.220, 1.090 },  
  { "NDI",    7, 1.304, 0.956 },  
  { "OTE",    8, 1.249, 0.637 },  
  { "OTR4",   8, 1.216, 0.569 },  
  { "OPI2",   8, 1.083, 0.274 },  
  { "STE",   16, 3.496, 3.000 },  
  { "STR4",  16, 3.827, 3.729 },  
  { "SPI2",  16, 2.982, 2.700 },  
  { "PTE",   15, 2.485, 1.538 },
  { "Si",    14, 2.700,  1.97  },
};
 #define emlNR asize(miller)
 
char *lbasis[eqmNR] = { 
  "HF/6-31G", "MP2/aug-cc-pVDZ", "MP2/aug-cc-pVTZ",
  "B3LYP/aug-cc-pVTZ","B3LYP/d-aug-cc-pVTZ", 
  "Ahc", "Ahp", "Bosque", "TW" 
};

char *ec_name[ecNR] = {
  "Acid","Alcohol", "Aldehyde", "Alkane", "Alkene", "Alkyne",
  "Amide", "Amine", "Aromatic", 
  "Ester", "Ether", 
  "Halogen", "Heterocyclic",
  "Ketone", 
  "Nitril", "Nitro", "Phospho", "Silicate",
  "Thio",   "Other" 
};

int mp_num_prop(t_molprop *mp,char *prop)
{
  int i,mnp=0;
  
  for(i=0; (i<mp->nexperiment); i++) {
    if (strcasecmp(mp->pname[i],prop) == 0)
      mnp++;
  }
  return mnp;
}

double mp_get_prop(t_molprop *mp,char *prop,int index)
{
  int i,n=0;
  
  for(i=0; (i<mp->nexperiment); i++) {
    if (strcasecmp(mp->pname[i],prop) == 0) {
      if (n == index)
	return mp->experiment[i];
      else
	n++;
    }
  }
  fprintf(stderr,"No %s value for %s (index chosen: %d)\n",
	  prop,mp->molname,index);
  
  return -1;
}

char *mp_get_ref_prop(t_molprop *mp,char *prop,int index)
{
  int i,n=0;
  
  for(i=0; (i<mp->nexperiment); i++) {
    if (strcasecmp(mp->pname[i],prop) == 0) {
      if (n == index)
	return mp->reference[i];
      else
	n++;
    }
  }
  fprintf(stderr,"No %s value for %s (index chosen: %d)\n",
	  prop,mp->molname,index);
  
  return NULL;
}

