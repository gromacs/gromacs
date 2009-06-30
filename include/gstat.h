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
 * Gromacs Runs On Most of All Computer Systems
 */

#ifndef _gstat_h
#define _gstat_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "statutil.h"
#include "mshift.h"
#include "rmpbc.h"

#ifdef CPLUSPLUS
extern "C" {
#endif

/***********************************************
 *
 *     A U T O C O R R E L A T I O N
 *
 ***********************************************/

extern real LegendreP(real x,unsigned long m);

#define eacNormal (1<<0)
#define eacCos    (1<<1)
#define eacVector (1<<2)
#define eacRcross (1<<3  | eacVector)
#define eacP0     (1<<4  | eacVector)
#define eacP1     (1<<5  | eacVector)
#define eacP2     (1<<6  | eacVector)
#define eacP3     (1<<7  | eacVector)
#define eacP4     (1<<8  | eacVector)

enum {
  effnNONE, effnEXP1, effnEXP2, effnEXP3,   effnVAC, 
  effnEXP5, effnEXP7, effnEXP9, effnERREST, effnNR
};

/* must correspond with 'leg' g_chi.c:727 */
enum { edPhi=0, edPsi, edOmega, edChi1, edChi2, edChi3, edChi4, edChi5, edChi6, edMax };

enum { edPrintST=0,edPrintRO } ; 

#define NHISTO 360
#define NONCHI 3
#define MAXCHI edMax-NONCHI
#define NROT 4  /* number of rotamers: 1=g(-), 2=t, 3=g(+), 0=other */ 

typedef struct {
  int minO,minC,H,N,C,O,Cn[MAXCHI+3];
} t_dihatms; /* Cn[0]=N, Cn[1]=Ca, Cn[2]=Cb etc. */

typedef struct {
  char name[12];
  int  resnr;
  int  index;       /* Index for amino acids (histograms) */
  int  j0[edMax];   /* Index in dih array (phi angle is first...) */
  t_dihatms  atm;
  int  b[edMax];
  int  ntr[edMax];
  real S2[edMax];
  real rot_occ[edMax][NROT];

} t_dlist;

extern int  nfp_ffn[effnNR];

extern const char *s_ffn[effnNR+2];

extern const char *longs_ffn[effnNR];

int sffn2effn(const char **sffn);
/* Returns the ffn enum corresponding to the selected enum option in sffn */

extern t_pargs *add_acf_pargs(int *npargs,t_pargs *pa);
/* Add options for autocorr to the current set of options.
 * *npargs must be initialised to the number of elements in pa,
 * it will be incremented appropriately.
 */

extern void cross_corr(int n,real f[],real g[],real corr[]);
/* Simple minded cross correlation algorithm */
  
extern real fit_acf(int ncorr,int fitfn,output_env_t oenv,bool bVerbose,
		    real tbeginfit,real tendfit,real dt,real c1[],real *fit);
  /* Fit an ACF to a given function */
  
extern void do_autocorr(const char *fn,output_env_t oenv,const char *title,
			int nframes,int nitem,real **c1,
			real dt,unsigned long mode,bool bAver);
/* Calls low_do_autocorr (see below). After calling add_acf_pargs */

extern void low_do_autocorr(const char *fn,output_env_t oenv,const char *title,
			    int  nframes,int nitem,int nout,real **c1,
			    real dt,unsigned long mode,int nrestart,
			    bool bAver,bool bNormalize,
			    bool bVerbose,real tbeginfit,real tendfit,
			    int nfitparm,int nskip);
/* 
 * do_autocorr calculates autocorrelation functions for many things.
 * It takes a 2 d array containing nitem arrays of length nframes
 * for each item the ACF is calculated.
 *
 * A number of "modes" exist for computation of the ACF
 *
 * if (mode == eacNormal) {
 *   C(t) = < X (tau) * X (tau+t) >
 * }
 * else if (mode == eacCos) {
 *   C(t) = < cos (X(tau) - X(tau+t)) >
 * }
 * else if (mode == eacVector) {
 *   C(t) = < X(tau) * X(tau+t)
 * }
 * else if (mode == eacP1) {
 *   C(t) = < cos (X(tau) * X(tau+t) >
 * }
 * else if (mode == eacP2) {
 *   C(t) = 1/2 * < 3 cos (X(tau) * X(tau+t) - 1 >
 * }
 * else if (mode == eacRcross) {
 *   C(t) = < ( X(tau) * X(tau+t) )^2 >
 * }
 *
 * For modes eacVector, eacP1, eacP2 and eacRcross the input should be 
 * 3 x nframes long, where each triplet is taken as a 3D vector
 *
 * For mode eacCos inputdata must be in radians, not degrees!
 *
 * Other parameters are:
 *
 * fn is output filename (.xvg) where the correlation function(s) are printed
 * title is the title in the output file
 * nframes is the number of frames in the time series
 * nitem is the number of items
 * c1 		is an array of dimension [ 0 .. nitem-1 ] [ 0 .. nframes-1 ]
 *    		on output, this array is filled with the correlation function
 *    		to reduce storage
 * nrestart     is the number of steps between restarts for direct ACFs 
 *              (i.e. without FFT) When set to 1 all points are used as
 *              time origin for averaging
 * dt 		is the time between frames
 * bAver 	If set, all ndih C(t) functions are averaged into a single 
 *       	C(t)
 * (bFour      	If set, will use fast fourier transform (FFT) for evaluating
 *            	the ACF: removed option, now on the command line only)
 * bNormalize 	If set, all ACFs will be normalized to start at 0
 * nskip        Determines whether steps a re skipped in the output 
 */
 
typedef struct {
  const char *name; /* Description of the J coupling constant */
  real A,B,C;   /* Karplus coefficients */
  real offset;  /* Offset for dihedral angle in histogram (e.g. -M_PI/3) */
  real Jc;      /* Resulting Jcoupling */
  real Jcsig;   /* Standard deviation in Jc */
} t_karplus;

extern void calc_distribution_props(int nh,int histo[],
				    real start,int  nkkk, t_karplus kkk[],
				    real *S2);
/* This routine takes a dihedral distribution and calculates
 * coupling constants and dihedral order parameters of it.
 *
 * nh      is the number of points
 * histo   is the array of datapoints which is assumed to span
 *         2 M_PI radians
 * start   is the starting angle of the histogram, this can be either 0
 *         or -M_PI 
 * nkkk    is the number of karplus sets (multiple coupling constants may be
 *         derived from a single angle)
 * kkk     are the constants for calculating J coupling constants using a
 *         Karplus equation according to
 *
 *                  2
 *         J = A cos theta + B cos theta + C
 *
 *         where theta is phi - offset (phi is the angle in the histogram)
 * offset  is subtracted from phi before substitution in the Karplus
 *         equation
 * S2      is the resulting dihedral order parameter
 *
 */
 
 
/***********************************************
 *
 *     F I T   R O U T I N E S
 *
 ***********************************************/
extern void do_expfit(int ndata,real c1[],real dt,
		      real begintimefit,real endtimefit);

extern void expfit(int n, real x[], real y[], real Dy[], 
		   real *a, real *sa, 
		   real *b, real *sb);
/* This procedure fits y=exp(a+bx) for n (x,y) pairs to determine a and b.
 * The uncertainties in the y values must be in the vector Dy.
 * The standard deviations of a and b, sa and sb, are also calculated.
 *
 * Routine from Computers in physics, 7(3) (1993), p. 280-285.
 */

extern void ana_dih_trans(char *fn_trans,char *fn_histo,
			  real **dih,int nframes,int nangles,
			  char *grpname,real t0,real dt,bool bRb,
                          output_env_t oenv);
/*
 * Analyse dihedral transitions, by counting transitions per dihedral
 * and per frame. The total number of transitions is printed to
 * stderr, as well as the average time between transitions.
 *
 * is wrapper to low_ana_dih_trans, which also passes in and out the 
     number of transitions per dihedral per residue. that uses struc dlist 
     which is not external, so pp2shift.h must be included. 

 * Dihedrals are supposed to be in either of three minima,
 * (trans, gauche+, gauche-)
 * 
 * fn_trans  output file name for #transitions per timeframe
 * fn_histo  output file name for transition time histogram
 * dih       the actual dihedral angles
 * nframes   number of times frames
 * nangles   number of angles
 * grpname   a string for the header of plots
 * t0,dt     starting time resp. time between time frames
 * bRb       determines whether the polymer convention is used
 *           (trans = 0)
 */

extern void low_ana_dih_trans(bool bTrans, char *fn_trans,
			      bool bHisto, char *fn_histo, int maxchi, 
			      real **dih, int nlist, t_dlist dlist[], 
                              int nframes, int nangles, char *grpname, 
                              int xity[], real t0, real dt, bool bRb, 
                              real core_frac, output_env_t oenv); 
  /* as above but passes dlist so can copy occupancies into it, and xity[] 
   *  (1..nangles, corresp to dih[this][], so can have non-3 multiplicity of
   * rotamers. Also production of xvg output files is conditional 
   * and the fractional width of each rotamer can be set ie for a 3 fold 
   * dihedral with core_frac = 0.5 only the central 60 degrees is assigned
   * to each rotamer, the rest goes to rotamer zero */ 



extern void read_ang_dih(char *trj_fn,
			 bool bAngles,bool bSaveAll,bool bRb,bool bPBC,
			 int maxangstat,int angstat[],
			 int *nframes,real **time,
			 int isize,atom_id index[],
			 real **trans_frac,
			 real **aver_angle,
			 real *dih[],
                         output_env_t oenv);
/* 
 * Read a trajectory and calculate angles and dihedrals.
 *
 * trj_fn      file name of trajectory
 * tpb_fn      file name of tpb file
 * bAngles     do we have to read angles or dihedrals
 * bSaveAll    do we have to store all in the dih array
 * bRb         do we have Ryckaert-Bellemans dihedrals (trans = 0)
 * bPBC        compute angles module 2 Pi
 * maxangstat  number of entries in distribution array
 * angstat     angle distribution
 * *nframes    number of frames read
 * time        simulation time at each time frame
 * isize       number of entries in the index, when angles 3*number of angles
 *             else 4*number of angles
 * index       atom numbers that define the angles or dihedrals
 *             (i,j,k) resp (i,j,k,l)
 * trans_frac  number of dihedrals in trans
 * aver_angle  average angle at each time frame
 * dih         all angles at each time frame
 */
 
extern void make_histo(FILE *log,
		       int ndata,real data[],int npoints,int histo[],
		       real minx,real maxx);
/* 
 * Make a histogram from data. The min and max of the data array can
 * be determined (if minx == 0 and maxx == 0)
 * and the index in the histogram is computed from
 * ind = npoints/(max(data) - min(data))
 *
 * log       write error output to this file
 * ndata     number of points in data
 * data      data points
 * npoints   number of points in histogram
 * histo     histogram array. This is NOT set to zero, to allow you
 *           to add multiple histograms
 * minx      start of the histogram
 * maxx      end of the histogram
 *           if both are 0, these values are computed by the routine itself
 */

extern void normalize_histo(int npoints,int histo[],real dx,real normhisto[]);
/* 
 * Normalize a histogram so that the integral over the histo is 1 
 *
 * npoints    number of points in the histo array
 * histo      input histogram
 * dx         distance between points on the X-axis
 * normhisto  normalized output histogram
 */

extern real fit_function(int eFitFn,real *parm,real x);
/* Returns the value of fit function eFitFn at x */

/* Use Levenberg-Marquardt method to fit to a nfitparm parameter exponential */
/* or to a transverse current autocorrelation function */
/* Or: "There is no KILL like OVERKILL", Dr. Ir. D. van der Spoel */
extern real do_lmfit(int ndata,real c1[],real sig[],real dt,real *x,
		     real begintimefit,real endtimefit,output_env_t oenv,
                     bool bVerbose, int eFitFn,real fitparms[],int fix);
/* Returns integral.
 * If x == NULL, the timestep dt will be used to create a time axis.
 * fix fixes fit parameter i at it's starting value, when the i'th bit
 * of fix is set. 
 */

extern real evaluate_integral(int n,real x[],real y[],real dy[],
			      real aver_start,real *stddev);
/* Integrate data in y, and, if given, use dy as weighting 
 * aver_start should be set to a value where the function has 
 * converged to 0.
 */

extern real print_and_integrate(FILE *fp,int n,real dt,
				real c[],real *fit,int nskip);
/* Integrate the data in c[] from 0 to n using trapezium rule.
 * If fp != NULL output is written to it
 * nskip determines whether all elements are written to the output file
 * (written when i % nskip == 0)
 * If fit != NULL the fit is also written.
 */
 
extern int get_acfnout(void);
/* Return the output length for the correlation function 
 * Works only AFTER do_auto_corr has been called!
 */

extern int get_acffitfn(void);
/* Return the fit function type.
 * Works only AFTER do_auto_corr has been called!
 */
 
  /* Routines from pp2shift (anadih.c etc.) */

extern void do_pp2shifts(FILE *fp,int nframes,
			 int nlist,t_dlist dlist[],real **dih);

extern bool has_dihedral(int Dih,t_dlist *dl);

extern t_dlist *mk_dlist(FILE *log, 
			 t_atoms *atoms, int *nlist,
			 bool bPhi, bool bPsi, bool bChi, bool bHChi,
			 int maxchi,int r0,int naa,char **aa);
			 
extern void pr_dlist(FILE *fp,int nl,t_dlist dl[],real dt,  int printtype,
		     bool bPhi, bool bPsi,bool bChi,bool bOmega, int maxchi);

extern int pr_trans(FILE *fp,int nl,t_dlist dl[],real dt,int Xi);

extern void mk_chi_lookup (int **lookup, int maxchi, real **dih, 
			   int nlist, t_dlist dlist[]) ; 

extern void mk_multiplicity_lookup (int *xity, int maxchi, real **dih, 
				    int nlist, t_dlist dlist[],int nangle) ; 

extern void get_chi_product_traj (real **dih,int nframes,int nangles, 
				  int nlist,int maxchi, t_dlist dlist[], 
                                  real time[], int **lookup,int *xity,
                                  bool bRb,bool bNormalize,
				  real core_frac,bool bAll,char *fnall,
                                  output_env_t oenv); 

extern void print_one (output_env_t oenv, const char *base,const char *name,
                       const char *title, const char *ylabel,int nf,
                       real time[],real data[]); 
		      
/* Routines from g_hbond */
extern void analyse_corr(int n,real t[],real ct[],real nt[],real kt[],
                         real sigma_ct[],real sigma_nt[],real sigma_kt[],
                         real fit_start,real temp,real smooth_tail_start,
                         output_env_t oenv);

extern void compute_derivative(int nn,real x[],real y[],real dydx[]);

#ifdef CPLUSPLUS
	     }
#endif
 
#endif


