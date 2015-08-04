/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#ifndef GMX_GMXANA_GSTAT_H
#define GMX_GMXANA_GSTAT_H

#include "gromacs/commandline/pargs.h"
#include "gromacs/legacyheaders/oenv.h"
#include "gromacs/legacyheaders/typedefs.h"
#include "gromacs/topology/index.h"

#ifdef __cplusplus
extern "C" {
#endif

struct gmx_residuetype_t;

/* must correspond with 'leg' g_chi.c:727 */
enum {
    edPhi = 0, edPsi, edOmega, edChi1, edChi2, edChi3, edChi4, edChi5, edChi6, edMax
};

enum {
    edPrintST = 0, edPrintRO
};

#define NHISTO 360
#define NONCHI 3
#define MAXCHI edMax-NONCHI
#define NROT 4  /* number of rotamers: 1=g(-), 2=t, 3=g(+), 0=other */

typedef struct {
    int minCalpha, minC, H, N, C, O, Cn[MAXCHI+3];
} t_dihatms; /* Cn[0]=N, Cn[1]=Ca, Cn[2]=Cb etc. */

typedef struct {
    char       name[12];
    int        resnr;
    int        index;     /* Index for amino acids (histograms) */
    int        j0[edMax]; /* Index in dih array (phi angle is first...) */
    t_dihatms  atm;
    int        b[edMax];
    int        ntr[edMax];
    real       S2[edMax];
    real       rot_occ[edMax][NROT];

} t_dlist;

typedef struct {
    const char *name;    /* Description of the J coupling constant */
    real        A, B, C; /* Karplus coefficients */
    real        offset;  /* Offset for dihedral angle in histogram (e.g. -M_PI/3) */
    real        Jc;      /* Resulting Jcoupling */
    real        Jcsig;   /* Standard deviation in Jc */
} t_karplus;

void calc_distribution_props(int nh, int histo[],
                             real start, int  nkkk, t_karplus kkk[],
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

void ana_dih_trans(const char *fn_trans, const char *fn_histo,
                   real **dih, int nframes, int nangles,
                   const char *grpname, real *time, gmx_bool bRb,
                   const output_env_t oenv);
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
 * time      array (size nframes) of times of trajectory frames
 * bRb       determines whether the polymer convention is used
 *           (trans = 0)
 */

void low_ana_dih_trans(gmx_bool bTrans, const char *fn_trans,
                       gmx_bool bHisto, const char *fn_histo, int maxchi,
                       real **dih, int nlist, t_dlist dlist[],
                       int nframes, int nangles, const char *grpname,
                       int multiplicity[], real *time, gmx_bool bRb,
                       real core_frac, const output_env_t oenv);
/* as above but passes dlist so can copy occupancies into it, and multiplicity[]
 *  (1..nangles, corresp to dih[this][], so can have non-3 multiplicity of
 * rotamers. Also production of xvg output files is conditional
 * and the fractional width of each rotamer can be set ie for a 3 fold
 * dihedral with core_frac = 0.5 only the central 60 degrees is assigned
 * to each rotamer, the rest goes to rotamer zero */



void read_ang_dih(const char *trj_fn,
                  gmx_bool bAngles, gmx_bool bSaveAll, gmx_bool bRb, gmx_bool bPBC,
                  int maxangstat, int angstat[],
                  int *nframes, real **time,
                  int isize, atom_id index[],
                  real **trans_frac,
                  real **aver_angle,
                  real *dih[],
                  const output_env_t oenv);
/*
 * Read a trajectory and calculate angles and dihedrals.
 *
 * trj_fn      file name of trajectory
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

void make_histo(FILE *log,
                int ndata, real data[], int npoints, int histo[],
                real minx, real maxx);
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

void normalize_histo(int npoints, int histo[], real dx, real normhisto[]);
/*
 * Normalize a histogram so that the integral over the histo is 1
 *
 * npoints    number of points in the histo array
 * histo      input histogram
 * dx         distance between points on the X-axis
 * normhisto  normalized output histogram
 */

/* Routines from pp2shift (anadih.c etc.) */

void do_pp2shifts(FILE *fp, int nframes,
                  int nlist, t_dlist dlist[], real **dih);

gmx_bool has_dihedral(int Dih, t_dlist *dl);

t_dlist *mk_dlist(FILE *log,
                  t_atoms *atoms, int *nlist,
                  gmx_bool bPhi, gmx_bool bPsi, gmx_bool bChi, gmx_bool bHChi,
                  int maxchi, int r0, struct gmx_residuetype_t *rt);

void pr_dlist(FILE *fp, int nl, t_dlist dl[], real dt,  int printtype,
              gmx_bool bPhi, gmx_bool bPsi, gmx_bool bChi, gmx_bool bOmega, int maxchi);

int pr_trans(FILE *fp, int nl, t_dlist dl[], real dt, int Xi);

void mk_chi_lookup (int **lookup, int maxchi,
                    int nlist, t_dlist dlist[]);

void mk_multiplicity_lookup (int *multiplicity, int maxchi,
                             int nlist, t_dlist dlist[], int nangle);

void get_chi_product_traj (real **dih, int nframes,
                           int nlist, int maxchi, t_dlist dlist[],
                           real time[], int **lookup, int *multiplicity,
                           gmx_bool bRb, gmx_bool bNormalize,
                           real core_frac, gmx_bool bAll, const char *fnall,
                           const output_env_t oenv);

void print_one (const output_env_t oenv, const char *base,
                const char *name,
                const char *title, const char *ylabel, int nf,
                real time[], real data[]);

/* Routines from g_hbond */
void analyse_corr(int n, real t[], real ct[], real nt[], real kt[],
                  real sigma_ct[], real sigma_nt[], real sigma_kt[],
                  real fit_start, real temp);

void compute_derivative(int nn, real x[], real y[], real dydx[]);

#ifdef __cplusplus
}
#endif

#endif
