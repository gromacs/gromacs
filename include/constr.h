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
static char *SRCID_constr_h = "$Id$";

#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#include<callf77.h>

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
		   t_nrnb *nrnb,        /* Performance measure          */
		   real lambda,         /* FEP lambda                   */
		   real *dvdlambda);    /* FEP force                    */
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
		    real mH, 		/* Mass of Hydrogen		*/
		    int *xerror);

extern void cshake(atom_id iatom[],int ncon,int *nnit,int maxnit,
		   real dist2[],real xp[],real rij[],real m2[],real omega,
		   real invmass[],real tt[],real lagr[],int *nerror);
/* Regular iterative shake */

extern void constrain(FILE *log,t_topology *top,t_inputrec *ir,int step,
		      t_mdatoms *md,int start,int homenr,
		      rvec *x,rvec *xprime,rvec *min_proj,matrix box,
		      real lambda,real *dvdlambda,t_nrnb *nrnb,
		      bool bCoordinates);
/*
 * When bCoordinates=TRUE constrains coordinates xprime using th
 * directions in x, min_proj is not used.
 *
 * When bCoordinates=FALSE, calculates the components xprime in
 * the constraint directions and subtracts these components from min_proj.
 * So when min_proj=xprime, the constraint components are projected out.
 *
 * Init_constraints must have be called once, before calling constrain.
 */

extern bool init_constraints(FILE *log,t_topology *top,t_inputrec *ir,
			     t_mdatoms *md,int start,int homenr,
			     bool bOnlyCoords);
/* Initialize constraints stuff */

/* C routines for LINCS algorithm */ 
extern void clincsp(rvec *x,rvec *f,rvec *fp,int ncons,
		    int *bla1,int *bla2,int *blnr,int *blbnb,
		    real *blc,real *blcc,real *blm,
		    int nrec,real *invmass,rvec *r,
		    real *vbo,real *vbn,real *vbt);

extern void clincs(rvec *x,rvec *xp,int ncons,
		   int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,
		   real *blc,real *blcc,real *blm,
		   int nit,int nrec,real *invmass,rvec *r,
		   real *vbo,real *vbn,real *vbt,real wangle,int *warn,
		   real *lambda);

extern void cconerr(real *max,real *rms,int *imax,rvec *xprime,
		    int ncons,int *bla1,int *bla2,real *bllen);

void lincs_warning(rvec *x,rvec *xprime,
		   int ncons,int *bla1,int *bla2,real *bllen,real wangle);


#ifdef USE_FORTRAN
extern void F77_FUNC(fsettle,FSETTLE)(int *nshake,int owptr[],
				      real b4[],real after[],
				      real *dOH,real *dHH,
				      real *mO,real *mH,int *error);
extern void F77_FUNC(fshake,FSHAKE)(atom_id iatom[],int *ncon,
				    int *nit, int *maxnit,
				    real dist2[],real xp[],
				    real rij[],real m2[],real *omega,
				    real invmass[],real tt[],
				    real lambda[],int *error);
extern void F77_FUNC(flincsp,FLINCSP)(real *x,real *f,real *fp,
				      int *nc, int *bla1,int *bla2,
				      int *blnr,int *blbnb,
				      real *blc,real *blcc,real *blm,
				      int *nrec,real *invmass,
				      real *r, real *rhs1, real *rhs2,
				      real *sol);
extern void F77_FUNC(flincs,FLINCS)(real *x,real *xp,int *nc,
				    int *bla1,int *bla2,int *blnr,
				    int *blbnb,real *bllen,
				    real *blc,real *blcc,real *blm,
				    int *nit,int *nrec,real *invmass,
				    real *r,real *temp1,real *temp2,
				    real *temp3,real *wangle,
				    int *warn,real *lambda);
#endif 
