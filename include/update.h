/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * Good gRace! Old Maple Actually Chews Slate
 */

#ifndef	_update_h
#define	_update_h

#ifdef HAVE_IDENT
#ident	"@(#) update.h 1.33 24 Jun 1996"
#endif /* HAVE_IDENT */

#include "typedefs.h"
#include "mshift.h"
#include "tgroup.h"
#include "network.h"
#include "force.h"

extern void tcoupl(bool bTC,t_grpopts *opts,t_groups *grps,real dt,real lamb);

extern real calc_temp(real ekin,int nrdf);
/* Calculate the temperature */

extern void calc_pres(matrix box,tensor ekin,tensor vir,tensor pres);
/* Calculate the pressure */

extern void update(int          natoms,	/* number of atoms in simulation */
		   int      	 start,
		   int          homenr,	/* number of home particles 	*/
		   int          step,
		   real         lambda, /* FEP scaling parameter */
		   real         *dvdlambda, /* FEP stuff */
		   t_inputrec   *ir,    /* input record with constants 	*/
		   bool         bFirstStep,   
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
		   t_edsamyn    *edyn);
     
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

/* SHAKE stuff */ 
extern int bshakef(FILE *log,		/* Log file			*/
		   int natoms,		/* Total number of atoms	*/
		   real invmass[],	/* Atomic masses		*/
		   int nblocks,		/* The number of shake blocks	*/
		   int sblock[],        /* The shake blocks             */
		   t_idef *idef,	/* The interaction def		*/
		   t_inputrec *ir,	/* Input record		        */
		   matrix box,		/* The box			*/
		   rvec x_s[],		/* Coords before update		*/
		   rvec xp[],		/* Output coords		*/
		   t_nrnb *nrnb);        /* Performance measure          */
/* Shake all the atoms blockwise. It is assumed that all the constraints
 * in the idef->shakes field are sorted, to ascending block nr. The
 * sblock array points into the idef->shakes.iatoms field, with block 0 
 * starting
 * at sblock[0] and running to ( < ) sblock[1], block n running from 
 * sblock[n] to sblock[n+1]. Array sblock should be large enough.
 * Return 0 when OK, -1 when shake-error
 */
extern void csettle(FILE *log,
		    int nshake,		/* Number of water molecules 	*/
		    int owptr[],	/* pointer to Oxygen in b4 & after */
		    real b4[],		/* Old coordinates		*/
		    real after[],	/* New coords, to be settled	*/
		    real dOH,		/* Constraint length Ox-Hyd	*/
		    real dHH, 		/* Constraint length Hyd-Hyd	*/
		    real mO,  		/* Mass of Oxygen		*/
		    real mH  		/* Mass of Hydrogen		*/
		    );

extern void cshake(atom_id iatom[],int ncon,int *nnit,int maxnit,
		   real dist2[],real xp[],real rij[],real m2[],
		   real invmass[],real tt[],int *nerror);
/* Regular iterative shake */

/* Fortran versions of shake and settle */
DECLAREF77(fsettle) (int *nshake,int owptr[],
		     real b4[],real after[],
		     real *dOH,real *dHH,real *mO,real *mH);
		     
DECLAREF77(fshake)  (atom_id iatom[],int *ncon,int *nit,int *maxnit,
		     real dist2[],real xp[],real rij[],real m2[],
		     real invmass[],real tt[],int *error);

/* C routines for LINCS algorithm */ 
void clincs(rvec *x,rvec *xp,int ncons,int ncm,int cmax,
              int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,
	      real *blc,real *blcc,real *blm,
              int nrec,real *invmass,rvec * r,
	      real *vbo,real *vbn,real *vbt,real wangle,int *warn,
	      real *lambda);

void clincsld(rvec *x,rvec *xp,int ncons,int ncm,int cmax,
              int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,
	      real *blcc,real *blm,
              int nrec,rvec * r,
	      real *rhs1,real *rhs2,real *sol,real wangle,int *warn);

void cconerr(real *max,real *rms,int *imax,rvec *xprime,
	     int ncons,int *bla1,int *bla2,real *bllen);

/* Fortran routines for LINCS algorithm */ 
DECLAREF77(flincs) (real *x,real *xp,int *nc,int *ncm,int *cmax,
		     int *bla1,int *bla2,int *blnr,
		     int *blbnb,real *bllen,real *blc,real *blcc,real *blm,
		     int *nrec,real *invmass,real *r,real *temp1,
		     real *temp2,real *temp3,real *wangle,int *warn,
		     real *lambda);

DECLAREF77(flincsld)  (real *x,real *xp,int *nc,int *ncm,int *cmax,
		     int *bla1,int *bla2,int *blnr,
		     int *blbnb,real *bllen,real *blcc,real *blm,int *nrec,
		     real *r,real *temp1,real *temp2,real *temp3,
		     real *wangle,int *warn);

DECLAREF77(fconerr) (real *max,real *rms,int *imax,
		     rvec *xprime,int *ncons,int *bla1,int *bla2,real *bllen);

#endif	/* _update_h */

