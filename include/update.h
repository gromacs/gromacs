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
 * Green Red Orange Magenta Azure Cyan Skyblue
 */

#ifndef _update_h
#define _update_h

static char *SRCID_update_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) update.h 1.33 24 Jun 1996"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "mshift.h"
#include "tgroup.h"
#include "network.h"
#include "force.h"
#include "pull.h"

extern void update(int          natoms,	/* number of atoms in simulation */
		   int      	 start,
		   int          homenr,	/* number of home particles 	*/
		   int          step,
		   real         lambda, /* FEP scaling parameter */
		   real         *dvdlambda, /* FEP stuff */
		   t_inputrec   *ir,    /* input record with constants 	*/
		   t_mdatoms    *md,
		   rvec         x[],	/* coordinates of home particles */
		   t_graph      *graph,
		   rvec         shift_vec[],	
		   rvec         force[],/* forces on home particles 	*/
		   rvec         delta_f[],
		   rvec         vold[],	/* Old velocities		   */
		   rvec         v[], 	/* velocities of home particles */
		   rvec         vt[],  	/* velocity at time t 		*/
		   tensor       pressure,/* instantaneous pressure tensor */
		   tensor       box,  	/* instantaneous box lengths 	*/
		   t_topology   *top,
		   t_groups     *grps,
		   tensor       vir_part,
		   t_commrec    *cr,
		   t_nrnb       *nrnb,
		   bool         bTYZ,
		   bool         bDoUpdate,
		   t_edsamyn    *edyn,
		   t_pull       *pulldata,
		   bool         bNEMD);
     
extern void calc_ke_part(bool bFirstStep,int start,int homenr,
			 rvec vold[],rvec v[],rvec vt[],
			 t_grpopts *opts,t_mdatoms *md,
			 t_groups *grps,t_nrnb *nrnb,
			 real lambda,real *dvdlambda);
/*
 * Compute the partial kinetic energy for home particles;
 * will be accumulated in the calling routine.
 * The tensor is
 *
 * Ekin = SUM(i) 0.5 m[i] v[i] (x) v[i]
 *    
 *     use v[i] = v[i] - u[i] when calculating temperature
 *
 * u must be accumulated already.
 *
 * Now also computes the contribution of the kinetic energy to the
 * free energy
 *
 */

/* Routines from coupling.c to do with Temperature, Pressure and coupling
 * algorithms.
 */
extern real run_aver(real old,real cur,int step,int nmem);

extern void tcoupl(bool bTC,t_grpopts *opts,t_groups *grps,
		   real dt,real SAfactor,int step,int nmem);
/* Compute temperature scaling factors */

extern real calc_temp(real ekin,real nrdf);
/* Calculate the temperature */

extern void calc_pres(int eBox,matrix box,
		      tensor ekin,tensor vir,tensor pres,real Elr);
/* Calculate the pressure. Unit of pressure is bar, If Elr != 0
 * a long range correction based on Ewald/PPPM is made (see c-code)
 */

extern void do_pcoupl(t_inputrec *ir,int step,tensor pres,
		      matrix box,int start,int nr_atoms,
		      rvec x[],unsigned short cFREEZE[],
		      t_nrnb *nrnb,ivec nFreeze[]);
		      
extern void correct_ekin(FILE *log,int start,int end,rvec v[],
			 rvec vcm,real mass[],real tmass,tensor ekin);
/* Correct ekin for vcm */

#endif	/* _update_h */

