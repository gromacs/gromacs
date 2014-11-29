/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014, by the GROMACS development team, led by
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

#ifndef _constr_h
#define _constr_h

#include "gromacs/legacyheaders/typedefs.h"

#ifdef __cplusplus
extern "C" {
#endif

struct t_pbc;

enum
{
    econqCoord,         /* Constrain coordinates (mass weighted)           */
    econqVeloc,         /* Constrain velocities (mass weighted)            */
    econqDeriv,         /* Constrain a derivative (mass weighted),         *
                         * for instance velocity or acceleration,          *
                         * constraint virial can not be calculated.        */
    econqDeriv_FlexCon, /* As econqDeriv, but only output flex. con.       */
    econqForce,         /* Constrain forces (non mass-weighted)            */
    econqForceDispl     /* Constrain forces (mass-weighted 1/0 for freeze) */
};

int n_flexible_constraints(struct gmx_constr *constr);
/* Returns the total number of flexible constraints in the system */

void too_many_constraint_warnings(int eConstrAlg, int warncount);
/* Generate a fatal error because of too many LINCS/SETTLE warnings */

gmx_shakedata_t shake_init();
/* Initializes and return the SHAKE data structure */

gmx_bool bshakef(FILE           *log,          /* Log file			*/
                 gmx_shakedata_t shaked,       /* Total number of atoms	*/
                 real            invmass[],    /* Atomic masses		*/
                 int             nblocks,      /* The number of shake blocks	*/
                 int             sblock[],     /* The shake blocks             */
                 t_idef         *idef,         /* The interaction def		*/
                 t_inputrec     *ir,           /* Input record		        */
                 rvec            x_s[],        /* Coords before update		*/
                 rvec            prime[],      /* Output coords		*/
                 t_nrnb         *nrnb,         /* Performance measure          */
                 real           *lagr,         /* The Lagrange multipliers     */
                 real            lambda,       /* FEP lambda                   */
                 real           *dvdlambda,    /* FEP force                    */
                 real            invdt,        /* 1/delta_t                    */
                 rvec           *v,            /* Also constrain v if v!=NULL  */
                 gmx_bool        bCalcVir,     /* Calculate r x m delta_r      */
                 tensor          vir_r_m_dr,   /* sum r x m delta_r            */
                 gmx_bool        bDumpOnError, /* Dump debugging stuff on error*/
                 int             econq);       /* which type of constraint is occurring */
/* Shake all the atoms blockwise. It is assumed that all the constraints
 * in the idef->shakes field are sorted, to ascending block nr. The
 * sblock array points into the idef->shakes.iatoms field, with block 0
 * starting
 * at sblock[0] and running to ( < ) sblock[1], block n running from
 * sblock[n] to sblock[n+1]. Array sblock should be large enough.
 * Return TRUE when OK, FALSE when shake-error
 */

gmx_settledata_t settle_init(real mO, real mH, real invmO, real invmH,
                             real dOH, real dHH);
/* Initializes and returns a structure with SETTLE parameters */

void csettle(gmx_settledata_t    settled,
             int                 nsettle,          /* Number of settles            */
             t_iatom             iatoms[],         /* The settle iatom list        */
             const struct t_pbc *pbc,              /* PBC data pointer, can be NULL */
             real                b4[],             /* Old coordinates              */
             real                after[],          /* New coords, to be settled    */
             real                invdt,            /* 1/delta_t                    */
             real               *v,                /* Also constrain v if v!=NULL  */
             int                 calcvir_atom_end, /* Calculate r x m delta_r up to this atom */
             tensor              vir_r_m_dr,       /* sum r x m delta_r            */
             int                *xerror
             );

void settle_proj(gmx_settledata_t settled, int econq,
                 int nsettle, t_iatom iatoms[],
                 const struct t_pbc *pbc,   /* PBC data pointer, can be NULL  */
                 rvec x[],
                 rvec *der, rvec *derp,
                 int CalcVirAtomEnd, tensor vir_r_m_dder);
/* Analytical algorithm to subtract the components of derivatives
 * of coordinates working on settle type constraint.
 */

void cshake(const atom_id iatom[], int ncon, int *nnit, int maxnit,
            const real dist2[], real xp[], const real rij[], const real m2[], real omega,
            const real invmass[], const real tt[], real lagr[], int *nerror);
/* Regular iterative shake */

void crattle(atom_id iatom[], int ncon, int *nnit, int maxnit,
             real dist2[], real vp[], real rij[], real m2[], real omega,
             real invmass[], real tt[], real lagr[], int *nerror, real invdt);

gmx_bool constrain(FILE *log, gmx_bool bLog, gmx_bool bEner,
                   gmx_constr_t constr,
                   t_idef *idef,
                   t_inputrec *ir,
                   t_commrec *cr,
                   gmx_int64_t step, int delta_step,
                   real step_scaling,
                   t_mdatoms *md,
                   rvec *x, rvec *xprime, rvec *min_proj,
                   gmx_bool bMolPBC, matrix box,
                   real lambda, real *dvdlambda,
                   rvec *v, tensor *vir,
                   t_nrnb *nrnb, int econq);
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
 * When bMolPBC=TRUE, assume that molecules might be broken: correct PBC.
 *
 * delta_step is used for determining the constraint reference lengths
 * when lenA != lenB or will the pull code with a pulling rate.
 * step + delta_step is the step at which the final configuration
 * is meant to be; for update delta_step = 1.
 *
 * step_scaling can be used to update coordinates based on the time
 * step multiplied by this factor. Thus, normally 1.0 is passed. The
 * SD1 integrator uses 0.5 in one of its calls, to correct positions
 * for half a step of changed velocities.
 *
 * If v!=NULL also constrain v by adding the constraint corrections / dt.
 *
 * If vir!=NULL calculate the constraint virial.
 *
 * Return TRUE if OK, FALSE in case of shake error
 *
 */

gmx_constr_t init_constraints(FILE *log,
                              gmx_mtop_t *mtop, t_inputrec *ir,
                              gmx_edsam_t ed, t_state *state,
                              t_commrec *cr);
/* Initialize constraints stuff */

void set_constraints(gmx_constr_t    constr,
                     gmx_localtop_t *top,
                     t_inputrec     *ir,
                     t_mdatoms      *md,
                     t_commrec      *cr);
/* Set up all the local constraints for the node */

/* The at2con t_blocka struct returned by the routines below
 * contains a list of constraints per atom.
 * The F_CONSTRNC constraints in this structure number consecutively
 * after the F_CONSTR constraints.
 */

t_blocka make_at2con(int start, int natoms,
                     t_ilist *ilist, t_iparams *iparams,
                     gmx_bool bDynamics, int *nflexiblecons);
/* Returns a block struct to go from atoms to constraints */

const t_blocka *atom2constraints_moltype(gmx_constr_t constr);
/* Returns the an array of atom to constraints lists for the moltypes */

const int **atom2settle_moltype(gmx_constr_t constr);
/* Returns the an array of atom to settle for the moltypes */

#define constr_iatomptr(nconstr, iatom_constr, iatom_constrnc, con) ((con) < (nconstr) ? (iatom_constr)+(con)*3 : (iatom_constrnc)+(con-nconstr)*3)
/* Macro for getting the constraint iatoms for a constraint number con
 * which comes from a list where F_CONSTR and F_CONSTRNC constraints
 * are concatenated.
 */

gmx_bool inter_charge_group_constraints(const gmx_mtop_t *mtop);
/* Returns if there are inter charge group constraints */

gmx_bool inter_charge_group_settles(const gmx_mtop_t *mtop);
/* Returns if there are inter charge group settles */

real *constr_rmsd_data(gmx_constr_t constr);
/* Return the data for determining constraint RMS relative deviations.
 * Returns NULL when LINCS is not used.
 */

real constr_rmsd(gmx_constr_t constr, gmx_bool bSD2);
/* Return the RMSD of the constraint, bSD2 selects the second SD step */

real *lincs_rmsd_data(gmx_lincsdata_t lincsd);
/* Return the data for determining constraint RMS relative deviations */

real lincs_rmsd(gmx_lincsdata_t lincsd, gmx_bool bSD2);
/* Return the RMSD of the constraint, bSD2 selects the second SD step */

gmx_lincsdata_t init_lincs(FILE *fplog, gmx_mtop_t *mtop,
                           int nflexcon_global, t_blocka *at2con,
                           gmx_bool bPLINCS, int nIter, int nProjOrder);
/* Initializes and returns the lincs data struct */

void set_lincs(t_idef *idef, t_mdatoms *md,
               gmx_bool bDynamics, t_commrec *cr,
               gmx_lincsdata_t li);
/* Initialize lincs stuff */

void set_lincs_matrix(gmx_lincsdata_t li, real *invmass, real lambda);
/* Sets the elements of the LINCS constraint coupling matrix */

real constr_r_max(FILE *fplog, gmx_mtop_t *mtop, t_inputrec *ir);
/* Returns an estimate of the maximum distance between atoms
 * required for LINCS.
 */

gmx_bool
constrain_lincs(FILE *log, gmx_bool bLog, gmx_bool bEner,
                t_inputrec *ir,
                gmx_int64_t step,
                gmx_lincsdata_t lincsd, t_mdatoms *md,
                t_commrec *cr,
                rvec *x, rvec *xprime, rvec *min_proj,
                matrix box, struct t_pbc *pbc,
                real lambda, real *dvdlambda,
                real invdt, rvec *v,
                gmx_bool bCalcVir, tensor vir_r_m_dr,
                int econ,
                t_nrnb *nrnb,
                int maxwarn, int *warncount);
/* Returns if the constraining succeeded */


/* helper functions for andersen temperature control, because the
 * gmx_constr construct is only defined in constr.c. Return the list
 * of blocks (get_sblock) and the number of blocks (get_nblocks).  */

int *get_sblock(struct gmx_constr *constr);

int get_nblocks(struct gmx_constr *constr);

#ifdef __cplusplus
}
#endif
#endif
