/*
 * $Id$
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
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif

#include "typedefs.h"

enum { econqCoord,         /* Constrain coordinates (mass weighted)     */ 
       econqVeloc,         /* Constrain velocities (mass weighted)      */
       econqDeriv,         /* Constrain a derivative (mass weighted),   *
			    * for instance velocity or acceleration,    *
			    * constraint virial can not be calculated.  */
       econqDeriv_FlexCon, /* As econqDeriv, but only output flex. con. */
       econqForce          /* Constrain forces (non mass-weighted)      */
};

extern int n_flexible_constraints(struct gmx_constr *constr);
/* Returns the total number of flexible constraints in the system */

extern void too_many_constraint_warnings(int eConstrAlg,int warncount);
/* Generate a fatal error because of too many LINCS/SETTLE warnings */

extern bool bshakef(FILE *log,		/* Log file			*/
		    int natoms,		/* Total number of atoms	*/
		    real invmass[],	/* Atomic masses		*/
		    int nblocks,	/* The number of shake blocks	*/
		    int sblock[],       /* The shake blocks             */
		    t_idef *idef,	/* The interaction def		*/
		    t_inputrec *ir,	/* Input record		        */
		    matrix box,		/* The box			*/
		    rvec x_s[],		/* Coords before update		*/
		    rvec xp[],		/* Output coords		*/
		    t_nrnb *nrnb,       /* Performance measure          */
		    real *lagr,         /* The Lagrange multipliers     */
		    real lambda,        /* FEP lambda                   */
		    real *dvdlambda,    /* FEP force                    */
		    real invdt,         /* 1/delta_t                    */
		    rvec *v,            /* Also constrain v if v!=NULL  */
		    bool bCalcVir,      /* Calculate r x m delta_r      */
		    tensor rmdr,        /* sum r x m delta_r            */
		    bool bDumpOnError); /* Dump debugging stuff on error*/
/* Shake all the atoms blockwise. It is assumed that all the constraints
 * in the idef->shakes field are sorted, to ascending block nr. The
 * sblock array points into the idef->shakes.iatoms field, with block 0 
 * starting
 * at sblock[0] and running to ( < ) sblock[1], block n running from 
 * sblock[n] to sblock[n+1]. Array sblock should be large enough.
 * Return TRUE when OK, FALSE when shake-error
 */
extern void csettle(FILE *log,
		    int nsettle,	/* Number of settles  	        */
		    t_iatom iatoms[],	/* The settle iatom list        */
		    real b4[],		/* Old coordinates		*/
		    real after[],	/* New coords, to be settled	*/
		    real dOH,		/* Constraint length Ox-Hyd	*/
		    real dHH, 		/* Constraint length Hyd-Hyd	*/
		    real mO,  		/* Mass of Oxygen		*/
		    real mH, 		/* Mass of Hydrogen		*/
		    real invdt,         /* 1/delta_t                    */
		    real *v,            /* Also constrain v if v!=NULL  */
		    bool bCalcVir,      /* Calculate r x m delta_r      */
		    tensor rmdr,        /* sum r x m delta_r            */
		    int *xerror);

extern void settle_proj(FILE *fp,int nsettle, t_iatom iatoms[],rvec x[],
			real dOH,real dHH,real invmO,real invmH,
			rvec *der,rvec *derp,
			bool bCalcVir,tensor rmdder);
/* Analytical algorithm to subtract the components of derivatives
 * of coordinates working on settle type constraint.
 */

extern void cshake(atom_id iatom[],int ncon,int *nnit,int maxnit,
		   real dist2[],real xp[],real rij[],real m2[],real omega,
		   real invmass[],real tt[],real lagr[],int *nerror);
/* Regular iterative shake */

extern bool constrain(FILE *log,bool bLog,bool bEner,
		      gmx_constr_t constr,
		      t_topology *top,
		      t_inputrec *ir,
		      t_commrec *cr,
		      int step,int delta_step,
		      t_mdatoms *md,
		      rvec *x,rvec *xprime,rvec *min_proj,matrix box,
		      real lambda,real *dvdlambda,
		      rvec *v,tensor *vir,
		      t_nrnb *nrnb,int econq);
/*
 * When econq=econqCoord constrains coordinates xprime using th
 * directions in x, min_proj is not used.
 *
 * When econq=econqDeriv, calculates the components xprime in
 * the constraint directions and subtracts these components from min_proj.
 * So when min_proj=xprime, the constraint components are projected out.
 *
 * When econq=econqDeriv_FlexCon, the same is done as with econqDeriv,
 * but only the components of the flexible constraints are stored.
 *
 * delta_step is used for determining the constraint reference lengths
 * when lenA != lenB or will the pull code with a pulling rate.
 * step + delta_step is the step at which the final configuration
 * is meant to be; for update delta_step = 1.
 *
 * If v!=NULL also constrain v by adding the constraint corrections / dt.
 *
 * If vir!=NULL calculate the constraint virial.
 *
 * Init_constraints must have be called once, before calling constrain.
 *
 * Return TRUE if OK, FALSE in case of shake error
 *
 */

extern gmx_constr_t init_constraints(FILE *log,
				     gmx_mtop_t *mtop,t_inputrec *ir, 
				     gmx_edsam_t ed,t_state *state,
				     t_commrec *cr);
/* Initialize constraints stuff */

extern void set_constraints(gmx_constr_t constr,
			    t_topology *top,
			    t_inputrec *ir,
			    t_mdatoms *md,
			    gmx_domdec_t *dd);
/* Set up all the local constraints for the node */

extern t_blocka *atom2constraints_moltype(gmx_constr_t constr);
/* Returns the an arry of atom to constraints lists for the moltypes */

extern bool inter_charge_group_constraints(gmx_mtop_t *mtop);
/* Returns if there are inter charge group constraints */

extern real *constr_rmsd_data(gmx_constr_t constr);
/* Return the data for determining constraint RMS relative deviations.
 * Returns NULL when LINCS is not used.
 */

extern real constr_rmsd(gmx_constr_t constr,bool bSD2);
/* Return the RMSD of the constraint, bSD2 selects the second SD step */

extern real *lincs_rmsd_data(gmx_lincsdata_t lincsd);
/* Return the data for determining constraint RMS relative deviations */

extern real lincs_rmsd(gmx_lincsdata_t lincsd,bool bSD2);
/* Return the RMSD of the constraint, bSD2 selects the second SD step */

gmx_lincsdata_t init_lincs(FILE *fplog,gmx_mtop_t *mtop,
			   int nflexcon_global,t_blocka *at2con,
			   bool bPLINCS,int nIter,int nProjOrder);
/* Initializes and returns the lincs data struct */

extern void set_lincs(t_idef *idef,int start,int homenr,
		      t_blocka *at2con,
		      bool bDynamics,gmx_domdec_t *dd,
		      gmx_lincsdata_t li);
/* Initialize lincs stuff */

extern void set_lincs_matrix(gmx_lincsdata_t li,real *invmass,real lambda);
/* Sets the elements of the LINCS constraint coupling matrix */

extern real constr_r_max(FILE *fplog,gmx_mtop_t *mtop,t_inputrec *ir);
/* Returns an estimate of the maximum distance between atoms
 * required for LINCS.
 */

extern bool constrain_lincs(FILE *log,bool bLog,bool bEner,
			    t_inputrec *ir,
			    int step,gmx_lincsdata_t lincsd,t_mdatoms *md,
			    gmx_domdec_t *dd,
			    rvec *x,rvec *xprime,rvec *min_proj,matrix box,
			    real lambda,real *dvdlambda,
			    real invdt,rvec *v,
			    bool bCalcVir,tensor rmdr,
			    int econ,
			    t_nrnb *nrnb,
			    int maxwarn,int *warncount);
/* Returns if the constraining succeeded */
