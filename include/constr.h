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
		    real mH, 		/* Mass of Hydrogen		*/
		    int *xerror);

extern void cshake(atom_id iatom[],int ncon,int *nnit,int maxnit,
		   real dist2[],real xp[],real rij[],real m2[],
		   real invmass[],real tt[],int *nerror);
/* Regular iterative shake */

extern void constrain(FILE *log,t_topology *top,t_inputrec *ir,int step,
		      t_mdatoms *md,int start,int homenr,
		      rvec *x,rvec *xprime,matrix box,
		      real lambda,real *dvdlambda,t_nrnb *nrnb);
/* Constrain coordinates xprime using the directions in x,
 * init_constraints must have be called once, before calling constrain      
 */

extern bool init_constraints(FILE *log,t_topology *top,t_inputrec *ir,
			     t_mdatoms *md,int start,int homenr);
/* Initialize constraints stuff */

/* C routines for LINCS algorithm */ 
extern void clincs(rvec *x,rvec *xp,int ncons,int ncm,int cmax,
		   int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,
		   real *blc,real *blcc,real *blm,
		   int nit,int nrec,real *invmass,rvec * r,
		   real *vbo,real *vbn,real *vbt,real wangle,int *warn,
		   real *lambda);

extern void clincsld(rvec *x,rvec *xp,int ncons,int ncm,int cmax,
		     int *bla1,int *bla2,int *blnr,int *blbnb,real *bllen,
		     real *blcc,real *blm,int nit,int nrec,rvec * r,
		     real *rhs1,real *rhs2,real *sol,real wangle,int *warn);

extern void cconerr(real *max,real *rms,int *imax,rvec *xprime,
		    int ncons,int *bla1,int *bla2,real *bllen);

void lincs_warning(rvec *x,rvec *xprime,
		   int ncons,int *bla1,int *bla2,real *bllen,real wangle);
