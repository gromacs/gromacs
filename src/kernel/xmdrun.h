/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef _xmdrun_h
#define _xmdrun_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include "typedefs.h"
#include "network.h"
#include "tgroup.h"
#include "filenm.h"
#include "mshift.h"
#include "force.h"
#include <time.h>
#include "edsam.h"
#include "mdebin.h"
#include "vcm.h"
#include "vsite.h"
#include "gmx_wallcycle.h"

/* This	file contains XMDRUN datatypes and function prototypes, grouped
 * neatly according to parts of the functionalisty
 */

/* GENERAL COUPLING THEORY (GCT) STUFF */
enum {
    eoPres, eoEpot, eoVir, eoDist, eoMu, eoForce, eoFx, eoFy, eoFz,
    eoPx, eoPy, eoPz,
    eoPolarizability, eoDipole, eoObsNR,
    eoMemory = eoObsNR, eoInter, eoUseVirial,  eoCombRule, eoNR
};
extern const char *eoNames[eoNR];

typedef struct {
    int      at_i, at_j;  /* Atom type # for i and j                    */
    int      eObs;        /* Observable to couple to                */
    gmx_bool bPrint;      /* Does this struct have to be printed		*/
    real     c6, c12;     /* Actual value of params			*/
    real     xi_6, xi_12; /* Constants for coupling C6 and C12      */
} t_coupl_LJ;

typedef struct {
    int      at_i, at_j;       /* Atom type # for i and j                       */
    int      eObs;             /* Observable to couple to               */
    gmx_bool bPrint;           /* Does this struct have to be printed		*/
    real     a, b, c;          /* Actual value of params			*/
    real     xi_a, xi_b, xi_c; /* Constants for coupling A, B and C         */
} t_coupl_BU;

typedef struct {
    int      at_i;   /* Atom type					*/
    int      eObs;   /* Observable to couple to                 */
    gmx_bool bPrint; /* Does this struct have to be printed		*/
    real     Q;      /* Actual value of charge			*/
    real     xi_Q;   /* Constant for coupling Q			*/
} t_coupl_Q;

typedef struct {
    int       type; /* Type number in the iparams struct	*/
    int       eObs; /* Observable to couple to              */
    t_iparams xi;   /* Parameters that need to be changed	*/
    t_iparams iprint;
} t_coupl_iparams;

typedef struct {
    real             act_value[eoObsNR];
    real             av_value [eoObsNR];
    real             ref_value[eoObsNR];
    gmx_bool         bObsUsed[eoObsNR];
    int              nLJ, nBU, nQ, nIP;
    t_coupl_LJ      *tcLJ;
    t_coupl_BU      *tcBU;
    t_coupl_Q       *tcQ;
    t_coupl_iparams *tIP;
    int              nmemory;
    gmx_bool         bInter;
    gmx_bool         bVirial;
    int              combrule;
} t_coupl_rec;

extern void write_gct(const char *fn, t_coupl_rec *tcr, t_idef *idef);

extern void read_gct(const char *fn, t_coupl_rec *tcr);

extern void comm_tcr(FILE *log, t_commrec *cr, t_coupl_rec **tcr);

extern void copy_ff(t_coupl_rec *tcr, t_forcerec *fr, t_mdatoms *md,
                    t_idef *idef);

extern t_coupl_rec *init_coupling(FILE *log, int nfile, const t_filenm fnm[],
                                  t_commrec *cr, t_forcerec *fr, t_mdatoms *md,
                                  t_idef *idef);

extern void calc_force(int natom, rvec f[], rvec fff[]);

extern void calc_f_dev(int natoms, real charge[], rvec x[], rvec f[],
                       t_idef *idef, real *xiH, real *xiS);

extern void do_coupling(FILE *log, const output_env_t oenv, int nfile,
                        const t_filenm fnm[],
                        t_coupl_rec *tcr, real t, int step, real ener[],
                        t_forcerec *fr, t_inputrec *ir, gmx_bool bMaster,
                        t_mdatoms *md, t_idef *idef, real mu_aver, int nmols,
                        t_commrec *cr, matrix box, tensor virial,
                        tensor pres, rvec mu_tot,
                        rvec x[], rvec f[], gmx_bool bDoIt);

/* CODE TO ADD SPECIAL 2-DIMENSIONAL LENNARD-JONES CORRECTION TO FORCES AND ENERGY */
extern void do_glas(FILE *log, int start, int homenr, rvec x[], rvec f[],
                    t_forcerec *fr, t_mdatoms *md, int atnr, t_inputrec *ir,
                    real ener[]);

extern real mol_dipole(int k0, int k1, rvec x[], real q[]);
/* Calculate total dipole for group of atoms */

extern real calc_mu_aver(t_commrec *cr, rvec x[], real q[], rvec mu,
                         t_block *mols, t_mdatoms *md, int gnx, atom_id grpindex[]);
/* Compute average dipole */

/********************************************************************/
/* Force field scanning stuff */
typedef struct {
    real     tol, f_max, npow, epot, fac_epot, fac_pres, fac_msf, pres;
    int      molsize, nmol;
    gmx_bool bComb, bVerbose, bLogEps;
} t_ffscan;


extern gmx_bool update_forcefield(FILE *fplog,
                                  int nfile, const t_filenm fnm[], t_forcerec *fr,
                                  int natoms, rvec x[], matrix box);
/* Modify the parameters. Return TRUE when the scan is finished. */

extern gmx_bool print_forcefield(FILE *fp, real ener[], int natoms, rvec f[],
                                 rvec fshake[], rvec x[], t_block *mols, real mass[],
                                 tensor pres);
/* Print results. Return TRUE when the scan is finished. */

extern void set_ffvars(t_ffscan *ff);
/* Set variables for force scanning */

#endif  /* _xmdrun_h */
