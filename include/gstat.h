/*
 * $Id$
 * 
 *       This source code is part of
 * 
 *        G   R   O   M   A   C   S
 * 
 * GROningen MAchine for Chemical Simulations
 * 
 *               VERSION 2.0
 * 
 * Copyright (c) 1991-1999
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 * 
 * Also check out our WWW page:
 * http://md.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 * 
 * And Hey:
 * Good ROcking Metal Altar for Chronical Sinners
 */

#ifndef _gstat_h
#define _gstat_h

static char *SRCID_gstat_h = "$Id$";

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
 
extern t_pargs *add_acf_pargs(int *npargs,t_pargs *pa);
/* Add options for autocorr to the current set of options.
 * *npargs must be initialised to the number of elements in pa,
 * it will be incremented appropriately.
 */

extern void do_autocorr(char *fn,char *title,int nframes,int nitem,real **c1,
			real dt,unsigned long mode,bool bAver);
/* Calls low_do_autocorr (see below). After calling add_acf_pargs */

extern void low_do_autocorr(char *fn,char *title,
			    int  nframes,int nitem,int nout,real **c1,
			    real dt,unsigned long mode,int nrestart,
			    bool bAver,bool bFour,bool bNormalize,
			    bool bVerbose,real tbeginfit,real tendfit,
			    int nfitparm);
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
 * bFour      	If set, will use fast fourier transform (FFT) for evaluating
 *            	the ACF
 * bNormalize 	If set, all ACFs will be normalized to start at 0
 */
 
typedef struct {
  char *name;   /* Description of the J coupling constant */
  real A,B,C;   /* Karplus coefficients */
  real offset;  /* Offset for dihedral angle in histogram (e.g. -M_PI/3) */
  real Jc;      /* Resulting Jcoupling */
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
			  char *grpname,real t0,real dt,bool bRb);
/*
 * Analyse dihedral transitions, by counting transitions per dihedral
 * and per frame. The total number of transitions is printed to
 * stderr, as well as the average time between transitions.
 *
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

extern void read_ang_dih(char *trj_fn,char *tpb_fn,
			 bool bAngles,bool bSaveAll,bool bRb,
			 int maxangstat,int angstat[],
			 int *nframes,real **time,
			 int isize,atom_id index[],
			 real **trans_frac,
			 real **aver_angle,
			 real *dih[]);
/* 
 * Read a trajectory and calculate angles and dihedrals.
 *
 * trj_fn      file name of trajectory
 * tpb_fn      file name of tpb file
 * bAngles     do we have to read angles or dihedrals
 * bSaveAll    do we have to store all in the dih array
 * bRb         do we have Ryckaert-Bellemans dihedrals (trans = 0)
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

/* Use Levenberg-Marquardt method to fit to a one parameter exponential */
/* Or: "There is no KILL like OVERKILL", Dr. Ir. D. van der Spoel */
extern real do_lmfit(int ndata,real c1[],real sig[],real dt,real *x,
		     real begintimefit,real endtimefit,
		     bool bVerbose,int nfitparm,
		     real fit[],real fitparms[],char *fix);
/* Returns integral. if fit != NULL, than the fitted function is copied
 * into the original data array.
 * If x == NULL, the timestep dt will be used to create a time axis.
 */

extern real print_and_integrate(FILE *fp,int n,real dt,real c[],int nskip);
/* Integrate the data in c[] from 0 to n using trapezium rule.
 * If fp != NULL output is written to it each nskip steps
 */
 
extern int get_acfnout(void);
/* Return the output length for the correlation function 
 * Works only AFTER do_auto_corr has been called!
 */
 
 
/* Least squares method of computing statistics */
typedef struct {
  double yy,yx,xx,sx,sy;
  int    np;
} t_lsq;

extern void init_lsq(t_lsq *lsq);
extern void done_lsq(t_lsq *lsq);
extern void add_lsq(t_lsq *lsq,real x,real y);
extern void get_lsq_ab(t_lsq *lsq,real *a,real *b);
extern int  npoints_lsq(t_lsq *lsq);
extern real aver_lsq(t_lsq *lsq);
extern real sigma_lsq(t_lsq *lsq);
extern real error_lsq(t_lsq *lsq);

#ifdef CPLUSPLUS
	     }
#endif
 
#endif


