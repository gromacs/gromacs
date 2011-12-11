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

#ifndef _force_h
#define _force_h


#include "typedefs.h"
#include "pbc.h"
#include "network.h"
#include "tgroup.h"
#include "vsite.h"
#include "genborn.h"

#ifdef __cplusplus
extern "C" {
#endif

static const char *sepdvdlformat="  %-30s V %12.5e  dVdl %12.5e\n";

void calc_vir(FILE *fplog,int nxf,rvec x[],rvec f[],tensor vir,
		     gmx_bool bScrewPBC,matrix box);
/* Calculate virial for nxf atoms, and add it to vir */

void f_calc_vir(FILE *fplog,int i0,int i1,rvec x[],rvec f[],tensor vir,
		       t_graph *g,rvec shift_vec[]);
/* Calculate virial taking periodicity into account */

real RF_excl_correction(FILE *fplog,
			       const t_forcerec *fr,t_graph *g,
			       const t_mdatoms *mdatoms,const t_blocka *excl,
			       rvec x[],rvec f[],rvec *fshift,const t_pbc *pbc,
			       real lambda,real *dvdlambda);
/* Calculate the reaction-field energy correction for this node:
 * epsfac q_i q_j (k_rf r_ij^2 - c_rf)
 * and force correction for all excluded pairs, including self pairs.
 */

void calc_rffac(FILE *fplog,int eel,real eps_r,real eps_rf,
		       real Rc,real Temp,
		       real zsq,matrix box,
		       real *kappa,real *krf,real *crf);
/* Determine the reaction-field constants */

void init_generalized_rf(FILE *fplog,
				const gmx_mtop_t *mtop,const t_inputrec *ir,
				t_forcerec *fr);
/* Initialize the generalized reaction field parameters */


/* In wall.c */
void make_wall_tables(FILE *fplog,const output_env_t oenv,
			     const t_inputrec *ir,const char *tabfn,
			     const gmx_groups_t *groups,
			     t_forcerec *fr);

real do_walls(t_inputrec *ir,t_forcerec *fr,matrix box,t_mdatoms *md,
		     rvec x[],rvec f[],real lambda,real Vlj[],t_nrnb *nrnb);



t_forcerec *mk_forcerec(void);

#define GMX_MAKETABLES_FORCEUSER  (1<<0)
#define GMX_MAKETABLES_14ONLY     (1<<1)

t_forcetable make_tables(FILE *fp,const output_env_t oenv,
                                const t_forcerec *fr, gmx_bool bVerbose,
                                const char *fn, real rtab,int flags);
/* Return tables for inner loops. When bVerbose the tables are printed
 * to .xvg files
 */
 
bondedtable_t make_bonded_table(FILE *fplog,char *fn,int angle);
/* Return a table for bonded interactions,
 * angle should be: bonds 0, angles 1, dihedrals 2
 */

/* Return a table for GB calculations */
t_forcetable make_gb_table(FILE *out,const output_env_t oenv,
                                  const t_forcerec *fr,
                                  const char *fn,
                                  real rtab);

void pr_forcerec(FILE *fplog,t_forcerec *fr,t_commrec *cr);

void
forcerec_set_ranges(t_forcerec *fr,
		    int ncg_home,int ncg_force,
		    int natoms_force,
		    int natoms_force_constr,int natoms_f_novirsum);
/* Set the number of cg's and atoms for the force calculation */

gmx_bool can_use_allvsall(const t_inputrec *ir, const gmx_mtop_t *mtop,
                             gmx_bool bPrintNote,t_commrec *cr,FILE *fp);
/* Returns if we can use all-vs-all loops.
 * If bPrintNote==TRUE, prints a note, if necessary, to stderr
 * and fp (if !=NULL) on the master node.
 */

void init_forcerec(FILE       *fplog,     
                          const output_env_t oenv,
			  t_forcerec *fr,   
			  t_fcdata   *fcd,
			  const t_inputrec *ir,   
			  const gmx_mtop_t *mtop,
			  const t_commrec  *cr,
			  matrix     box,
			  gmx_bool       bMolEpot,
			  const char *tabfn,
			  const char *tabpfn,
			  const char *tabbfn,
			  gmx_bool       bNoSolvOpt,
			  real       print_force);
/* The Force rec struct must be created with mk_forcerec 
 * The gmx_booleans have the following meaning:
 * bSetQ:    Copy the charges [ only necessary when they change ]
 * bMolEpot: Use the free energy stuff per molecule
 * print_force >= 0: print forces for atoms with force >= print_force
 */

void init_enerdata(int ngener,int n_flambda,gmx_enerdata_t *enerd);
/* Intializes the energy storage struct */

void destroy_enerdata(gmx_enerdata_t *enerd);
/* Free all memory associated with enerd */

void reset_enerdata(t_grpopts *opts,
			   t_forcerec *fr,gmx_bool bNS,
			   gmx_enerdata_t *enerd,
			   gmx_bool bMaster);
/* Resets the energy data, if bNS=TRUE also zeros the long-range part */

void sum_epot(t_grpopts *opts,gmx_enerdata_t *enerd);
/* Locally sum the non-bonded potential energy terms */

void sum_dhdl(gmx_enerdata_t *enerd,double lambda,t_inputrec *ir);
/* Sum the free energy contributions */

void update_forcerec(FILE *fplog,t_forcerec *fr,matrix box);
/* Updates parameters in the forcerec that are time dependent */

/* Compute the average C6 and C12 params for LJ corrections */
void set_avcsixtwelve(FILE *fplog,t_forcerec *fr,
			     const gmx_mtop_t *mtop);

/* The state has changed */
#define GMX_FORCE_STATECHANGED (1<<0)
/* The box might have changed */
#define GMX_FORCE_DYNAMICBOX   (1<<1)
/* Do neighbor searching */
#define GMX_FORCE_NS           (1<<2)
/* Calculate bonded energies/forces */
#define GMX_FORCE_DOLR         (1<<3)
/* Calculate long-range energies/forces */
#define GMX_FORCE_BONDED       (1<<4)
/* Store long-range forces in a separate array */
#define GMX_FORCE_SEPLRF       (1<<5)
/* Calculate non-bonded energies/forces */
#define GMX_FORCE_NONBONDED    (1<<6)
/* Calculate forces (not only energies) */
#define GMX_FORCE_FORCES       (1<<7)
/* Calculate the virial */
#define GMX_FORCE_VIRIAL       (1<<8)
/* Calculate dHdl */
#define GMX_FORCE_DHDL         (1<<9)
/* Normally one want all energy terms and forces */
#define GMX_FORCE_ALLFORCES    (GMX_FORCE_BONDED | GMX_FORCE_NONBONDED | GMX_FORCE_FORCES)

void do_force(FILE *log,t_commrec *cr,
		     t_inputrec *inputrec,
		     gmx_large_int_t step,t_nrnb *nrnb,gmx_wallcycle_t wcycle,
		     gmx_localtop_t *top,
		     gmx_mtop_t *mtop,
		     gmx_groups_t *groups,
		     matrix box,rvec x[],history_t *hist,
		     rvec f[],
		     tensor vir_force,
		     t_mdatoms *mdatoms,
		     gmx_enerdata_t *enerd,t_fcdata *fcd,
		     real lambda,t_graph *graph,
		     t_forcerec *fr,gmx_vsite_t *vsite,rvec mu_tot,
		     double t,FILE *field,gmx_edsam_t ed,
		     gmx_bool bBornRadii,
		     int flags);
/* Communicate coordinates (if parallel).
 * Do neighbor searching (if necessary).
 * Calculate forces.
 * Communicate forces (if parallel).
 * Spread forces for vsites (if present).
 *
 * f is always required.
 */

void ns(FILE       *fplog,
	       t_forcerec *fr,
	       rvec       x[],
	       matrix     box,
	       gmx_groups_t *groups,
	       t_grpopts  *opts,
	       gmx_localtop_t *top,
	       t_mdatoms  *md,
	       t_commrec  *cr,
	       t_nrnb     *nrnb,
	       real       lambda,
	       real       *dvdlambda,
	       gmx_grppairener_t *grppener,
	       gmx_bool       bFillGrid,
	       gmx_bool       bDoLongRange,
	       gmx_bool       bDoForces,
	       rvec       *f);
/* Call the neighborsearcher */

void do_force_lowlevel(FILE         *fplog,  
			      gmx_large_int_t   step,
			      t_forcerec   *fr,
			      t_inputrec   *ir,
			      t_idef       *idef,
			      t_commrec    *cr,
			      t_nrnb       *nrnb,
			      gmx_wallcycle_t wcycle,
			      t_mdatoms    *md,
			      t_grpopts    *opts,
			      rvec         x[],
			      history_t    *hist,
			      rvec         f[],    
			      gmx_enerdata_t *enerd,
			      t_fcdata     *fcd,
			      gmx_mtop_t     *mtop,
			      gmx_localtop_t *top,
			      gmx_genborn_t *born,
			      t_atomtypes  *atype,
			      gmx_bool         bBornRadii,
			      matrix       box,
			      real         lambda,
			      t_graph      *graph,
			      t_blocka     *excl,
			      rvec         mu_tot[2],
			      int          flags,
			      float        *cycles_pme);
/* Call all the force routines */

#ifdef __cplusplus
}
#endif

#endif	/* _force_h */
