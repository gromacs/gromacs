/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.1
 * Copyright (c) 1991-2001, University of Groningen, The Netherlands
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
static char *SRCID_ion_data_h = "$Id$";
	
typedef struct {
  real photo,coh,incoh,incoh_abs;
} t_cross;

/* THIS TABLE HAS ADDED A 12 keV COLUMN TO HYDROGEN, CARBON,  */ 
/* OXYGEN, NITROGEN AND SULPHUR BY FITTING A QUADRATIC TO THE */
/* POINTS 8keV, 10keV and 12keV - now contains 6, 8, 10, 12,  */ 
/* 15 and 20 keV                                              */
/* Units are barn. They are converted to nm^2 by multiplying  */
/* by 1e-10, which is done in Imax (ionize.c)                 */

static t_cross   cross_sec_h[] = {
  { 2.63e-2,     1.01e-1,      5.49e-1,         7.12e-3 },
  { 9.79e-3,     6.18e-2,      5.83e-1,         9.60e-3 },
  { 4.55e-3,     4.16e-2,      5.99e-1,         1.19e-2 },
  { 1.52e-3,     2.79e-2,      6.08e-1,         1.41e-2 },
  { 1.12e-3,     1.96e-2,      6.09e-1,         1.73e-2 },
  { 4.16e-4,     1.13e-2,      6.07e-1,         2.23e-2 }
};
static t_cross   cross_sec_c[] = {
  { 1.99e+2,     5.88e+0,      2.29e+0,         3.06e-2 },
  { 8.01e+1,     4.22e+0,      2.56e+0,         4.38e-2 },
  { 3.92e+1,     3.26e+0,      2.74e+0,         5.72e-2 },
  { 2.19e+1,     2.55e+0,      2.88e+0,         7.07e-2 },
  { 1.06e+1,     1.97e+0,      3.04e+0,         9.15e-2 },
  { 4.15e+0,     1.30e+0,      3.20e+0,         1.24e-1 }
};
static t_cross   cross_sec_n[] = {
  { 3.91e+2,     8.99e+0,      2.49e+0,         3.43e-2 },
  { 1.59e+2,     6.29e+0,      2.86e+0,         5.01e-2 },
  { 7.88e+1,     4.76e+0,      3.10e+0,         6.57e-2 },
  { 4.42e+1,     3.66e+0,      3.28e+0,         8.13e-2 },
  { 2.16e+1,     2.82e+0,      3.46e+0,         1.05e-1 },
  { 8.52e+0,     1.88e+0,      3.65e+0,         1.43e-1 }
};
static t_cross   cross_sec_o[] = {
  { 6.90e+2,     1.33e+1,      2.66e+0,         3.75e-2 },
  { 2.84e+2,     9.21e+0,      3.14e+0,         5.62e-2 },
  { 1.42e+2,     6.85e+0,      3.44e+0,         7.43e-2 },
  { 8.01e+1,     5.18e+0,      3.66e+0,         9.20e-2 },
  { 3.95e+1,     3.97e+0,      3.87e+0,         1.18e-1 },
  { 1.57e+1,     2.64e+0,      4.10e+0,         1.61e-1 }
};
static t_cross   cross_sec_s[] = {
  { 1.10e+4,      5.54e+1,     3.98e+0,         5.42e-2 },
  { 4.91e+3,      4.29e+1,     4.71e+0,         8.38e-2 },
  { 2.58e+3,      3.36e+1,     5.32e+0,         1.16e-1 },
  { 1.52e+3,      2.64e+1,     5.81e+0,         1.48e-1 },
  { 7.82e+2,      1.97e+1,     6.36e+0,         2.00e-1 },
  { 3.29e+2,      1.29e+1,     6.94e+0,         2.80e-1 }
};

typedef struct {
  char *name;
  int  nel;
  t_cross *cross;
} t_element;

static t_element element[] = {
  { "H",   1, cross_sec_h },
  { "C",   6, cross_sec_c },
  { "N",   7, cross_sec_n },
  { "O",   8, cross_sec_o },
  { "S",  16, cross_sec_s },
  { "AR", 20, cross_sec_s },  /* This is not correct! */
  { "EL",  1, cross_sec_h }   /* This is not correct! */
};
#define NELEM asize(element)

/* 
 * In the first column the binding energy of the K-electrons; 
 * THIS IS IN eV,  which matches the photon energies. 
 * In the second column the binding energy of the outer shell electrons
 * The third column describes the photoelectric cross sections, 
 * where this now gives the fraction of photoelectric events 
 * which correspond to K-shell events, I called f_j in my 
 * notes:
 * The final column (a new column) now gives the values for the lifetimes
 * in ps. 
 */
typedef struct {
  real E_K,E_L,Prob_K,tau;
} t_recoil;

t_recoil recoil[] = {
  { 0.0,    0.0,   0.0,   0},                
  { 0.0136, 0.0,   0.0,   0},
  { 0.0246, 0.0,   0.0,   0},
  { 0.055,  0.005, 0.960, 0.012},
  { 0.117,  0.009, 0.956, 0.012},
  { 0.192,  0.008, 0.952, 0.012}, 
  { 0.284,  0.011, 0.950, 0.0113},
  { 0.402,  0.015, 0.950, 0.0083},
  { 0.532,  0.014, 0.936, 0.0066},
  { 0.687,  0.017, 0.928, 0.0045},
  { 0.874,  0.031, 0.922, 0.0033},
  { 1.072,  0.041, 0.933, 0.0028},
  { 1.305,  0.054, 0.927, 0.0022},
  { 1.560,  0.077, 0.922, 0.0019},
  { 1.839,  0.105, 0.918, 0.00165},
  { 2.146,  0.133, 0.921, 0.00145},
  { 2.472,  0.166, 0.908, 0.00130},
  { 2.822,  0.212, 0.902, 0.0012},
  { 3.203,  0.247, 0.902, 0.0010},
  { 3.607,  0.298, 0.894, 0.00095},
  { 4.038,  0.348, 0.890, 0.00085},
  { 4.490,  0.404, 0.886, 0.00078},
  { 4.966,  0.458, 0.882, 0.00073},
  { 5.465,  0.516, 0.885, 0.00062},
  { 5.989,  0.578, 0.883, 0.00055},
  { 6.539,  0.645, 0.880, 0.00049},
  { 7.112,  0.713, 0.877, 0.00044}
};





