/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2004, The GROMACS development team.
 * Copyright (c) 2013,2014,2015,2016,2017,2018, by the GROMACS development team, led by
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

#ifndef GMX_MBLIB_CONSTR_H
#define GMX_MBLIB_CONSTR_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/basedefinitions.h"
#include "gromacs/utility/real.h"

struct gmx_edsam;
struct gmx_localtop_t;
struct gmx_mtop_t;
struct gmx_multisim_t;
struct t_blocka;
struct t_commrec;
struct t_idef;
struct t_ilist;
struct t_inputrec;
struct t_mdatoms;
struct t_nrnb;
struct t_pbc;
union t_iparams;
class t_state;

/* Abstract type for constraints */
typedef struct gmx_constr *gmx_constr_t;

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

int n_flexible_constraints(gmx_constr *constr);
/* Returns the total number of flexible constraints in the system */

void too_many_constraint_warnings(int eConstrAlg, int warncount);
/* Generate a fatal error because of too many LINCS/SETTLE warnings */


gmx_bool constrain(FILE *log, gmx_bool bLog, gmx_bool bEner,
                   gmx_constr_t constr,
                   t_idef *idef,
                   const t_inputrec *ir,
                   const t_commrec *cr,
                   const gmx_multisim_t *ms,
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
                              const gmx_mtop_t *mtop, const t_inputrec *ir,
                              bool doEssentialDynamics,
                              const t_commrec *cr);
/* Initialize constraints stuff */

void saveEdsamPointer(gmx_constr_t      constr,
                      gmx_edsam        *ed);
/* Put a pointer to the essential dynamics constraints into the constr struct */

void set_constraints(gmx_constr_t             constr,
                     gmx_localtop_t          *top,
                     const t_inputrec        *ir,
                     const t_mdatoms         *md,
                     const t_commrec         *cr);
/* Set up all the local constraints for the node */

/* The at2con t_blocka struct returned by the routines below
 * contains a list of constraints per atom.
 * The F_CONSTRNC constraints in this structure number consecutively
 * after the F_CONSTR constraints.
 */

t_blocka make_at2con(int start, int natoms,
                     const t_ilist *ilist, const t_iparams *iparams,
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

real constr_rmsd(gmx_constr_t constr);
/* Return the RMSD of the constraint */

#endif
