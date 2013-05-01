/*
 * $Id: gmx_resp.hpp,v 1.3 2009/01/28 00:04:17 spoel Exp $
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.0.99
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
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

#ifndef _gmx_resp_hpp
#define _gmx_resp_hpp
	
#include <stdio.h>
#include "grompp.h"
#include "poldata.h"
#include "atomprop.h"
#include "gmx_statistics.h"
#include "molprop.hpp"

typedef struct 
{
    int  *row,atomnumber,atype;
    bool bRestrained;
    char *atomtype;
    int  nZeta;
    real *q,*zeta,*zeta_ref;
    int  *iq,*iz;
} gmx_ra;

typedef struct gmx_resp 
{
    int    nesp,nrho,natom,natype,ngridp,iModel;
    double qtot,qsum,watoms;
    double rms,rrms,penalty,pfac,entropy,wtot;
    dvec   origin,space;
    bool   bZatype,bFitZeta,bEntropy;
    bool   bRandZeta,bRandQ;
    bool   bAXpRESP;
    ivec   nxyz;
    real   qfac,b_hyper,zmin,zmax,delta_z,qmin,qmax,rDecrZeta;
    int    nparam; /* Total number of parameters */
    gmx_ra *ra;
    char   **dzatoms;
    const char   *stoichiometry;
    double *pot,*pot_calc,*rho;
    rvec   *x,*esp;
} gmx_resp;

typedef struct gmx_resp *gmx_resp_t;

#ifdef __cplusplus
extern "C" {
#endif
int atomicnumber2row(int elem);

void gmx_ra_init(gmx_ra *ra,int atomnumber,int atype,
                 const char *atomtype,gmx_poldata_t pd,
                 int iModel,char **dzatoms);

gmx_resp_t gmx_resp_init(int eqg_model,
                         bool bAXpRESP,real qfac,real b_hyper,real qtot,
                         real zmin,real zmax,real delta_z,bool bZatyp,
                         real watoms,real rDecrZeta,
                         bool bRandZeta,real penalty_fac,bool bFitZeta,
                         bool bEntropy,const char *dzatoms);
				
void gmx_resp_statistics(gmx_resp_t gr,int len,char buf[]);
  
void gmx_resp_summary(FILE *gp,gmx_resp_t gr,int *symmetric_atoms);
  
void gmx_resp_update_atomtypes(gmx_resp_t gr,t_atoms *atoms);

void gmx_resp_fill_zeta(gmx_resp_t gr,gmx_poldata_t pd);

void gmx_resp_fill_q(gmx_resp_t gr,t_atoms *atoms);

void gmx_resp_add_atom_coords(gmx_resp_t gr,rvec *x);

void gmx_resp_add_atom_info(gmx_resp_t gr,t_atoms *atoms,
                            gmx_poldata_t pd);

void gmx_resp_get_atom_info(gmx_resp_t gr,t_atoms *atoms,
                            t_symtab *symtab,rvec **x);

const char *gmx_resp_get_stoichiometry(gmx_resp_t gr);

void gmx_resp_add_atom_symmetry(gmx_resp_t gr,
                                int *symmetric_atoms);

void gmx_resp_add_point(gmx_resp_t gr,double x,double y,
                        double z,double V);

void gmx_resp_make_grid(gmx_resp_t gr,real spacing,matrix box,rvec x[]);

void gmx_resp_copy_grid(gmx_resp_t dest,gmx_resp_t src);

gmx_resp_t gmx_resp_copy(gmx_resp_t src);

void gmx_resp_calc_rms(gmx_resp_t gr);

double gmx_resp_get_rms(gmx_resp_t gr,real *wtot);

void gmx_resp_pot_lsq(gmx_resp_t gr,gmx_stats_t lsq);

void gmx_resp_calc_rho(gmx_resp_t gr);

void gmx_resp_calc_pot(gmx_resp_t gr);

void gmx_resp_read_cube(gmx_resp_t gr,const char *fn,bool bESPonly);

void gmx_resp_write_cube(gmx_resp_t gr,const char *fn,char *title);

void gmx_resp_write_rho(gmx_resp_t gr,const char *fn,char *title);

void gmx_resp_write_diff_cube(gmx_resp_t grref,gmx_resp_t gr,
                              const char *cube_fn,const char *hist_fn,
                              char *title,output_env_t oenv,int rho);

void gmx_resp_write_histo(gmx_resp_t gr,const char *fn,
                          char *title,output_env_t oenv);

int  gmx_resp_optimize_charges(FILE *fp,gmx_resp_t grt,int maxiter,
                               real toler,real *rms);

void gmx_resp_potcomp(gmx_resp_t gr,const char *potcomp,
                      const char *pdbdiff,output_env_t oenv);

void gmx_resp_destroy(gmx_resp_t grt);

double gmx_resp_get_qtot(gmx_resp_t grt,int atom);

double gmx_resp_get_q(gmx_resp_t grt,int atom,int zz);

double gmx_resp_get_zeta(gmx_resp_t grt,int atom,int zz);

void gmx_resp_set_q(gmx_resp_t grt,int atom,int zz,double q);

void gmx_resp_set_zeta(gmx_resp_t grt,int atom,int zz,double zeta);

#ifdef __cplusplus
}
#endif

#endif
