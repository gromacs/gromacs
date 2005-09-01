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
 * GROwing Monsters And Cloning Shrimps
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>
#include <time.h>
#include "typedefs.h"
#include "string2.h"
#include "smalloc.h"
#include "names.h"
#include "confio.h"
#include "mvdata.h"
#include "txtdump.h"
#include "vec.h"
#include "time.h"
#include "nrnb.h"
#include "mshift.h"
#include "mdrun.h"
#include "update.h"
#include "physics.h"
#include "rmpbc.h"
#include "nrjac.h"
#include "edsam.h"



#define EPS  1.0e-9

/*************************** FLOODING ************************************

The flooding ability was added later to edsam. Many of the edsam functionality could be reused for that purpose. 
The flooding covariance matrix, i.e. the selected eigenvectors and their corresponding eigenvalues are 
read as 7th Component Group. The eigenvalues are coded into the stepsize parameter (as used by -linfix or -linacc). 

do_md clls right in the beginning the function init_edsam, which reads the edi file, saves all the necessary information in
the edi structure and calls init_flood, to initialise some extra fields in the edi->flood structure.

since the flooding acts on forces do_flood is called from the function force() (force.c), while the other edsam functionality is hooked
into md via the update() (update.c) function acting as constraint on positions. 
flooding works correctly even if do_edsam() is not called.

do_flood makes a copy of the positions,
fits them, projects them computes flooding_energy, and flooding forces. The forces are computed in the 
space of the eigenvectors and are then blown up to the full cartesian space and rotated back to remove the
fit. Then do_flood adds these forces to the forcefield-forces
(given as parameter) and updates the adaptive flooding parameters Efl and deltaF.

To center the flooding potential at a different location one can use the -ori option in make_edi. The ori
structure is projected to the system of eigenvectors and then this position in the subspace is used as
center of the flooding potential.   If the option is not used, the center will be zero in the subspace,
i.e. the average structure as given in the make_edi file.

To use the flooding potential as restraint, make_edi has the option -restrain, which leads to inverted
signs of alpha2 and Efl, such that the sign in the exponential of Vfl is not inverted but the sign of
Vfl is inverted. Vfl = Efl * exp (- .../Efl/alpha2*x^2...) With tau>0 the negative Efl will grow slowly
so that the restraint is switched off slowly. When Efl==0 and inverted flooding is ON is reached no
 further adaption is applied, Efl will stay constant at zero. 

to use restraints with harmonic potentials switch -restrain and -harmonic. Then the eigenvalues are 
used as spring constants for the harmonic potential. 

to use more than one flooding matrix just concatenate severale .edi files (cat flood1.edi flood2.edi > flood_all.edi )
the routine read_edi_file reads all of theses flooding files.
The structure t_edi is now organized as a list of t_edis  and the function do_flood cycles through the list calling the do_single_flood() routine for every single entry. Since every state variables have been kept in one edi there is no interdependence whatsoever.
The forces are added together. 

  To write energies into the .edr file, call the function 
        get_flood_enx_names(char**, int *nnames) to get the Header (Vfl1 Vfl2... Vfln)
and call
        get_flood_energies(real Vfl[],int nnames); 

  TODO:
- one could program the whole thing such that Efl, Vfl and deltaF is written to the .edr file. -- i don't know how to do that, yet.

  In the moment one can have multiple flooding matrices, but only the first input is used for edsam. Especially with the angular motion
  remover there might be a market for multiple edsam inputs as well. 
 
  Secondly: Maybe one should give a range of atoms for which to remove motion, so that motion is removed with two edsam files from two peptide chains
*/

#define FIT 1
#define BLOWUP 2
#define DOFIT
#define NODEBUG BLOWUP
#ifdef DEBUG
   #define DUMP_FORCES 
   #define DEBUG_PRINT(X) fprintf(stderr,"%s\n",(X))
#else
   #define DEBUG_PRINT(X) {}
#endif

#ifdef DUMP_FORCES
   static FILE* logfile = NULL;
#endif

void do_flood(FILE *log, t_commrec *cr, rvec x[],rvec force[], t_edsamyn *edyn, int step);
/* main hook - called from do_force() in mdrun*/

void finish_edsam(FILE *log,t_topology *top,t_inputrec *ir,
		t_mdatoms *md,int start,int homenr,t_commrec *cr,
	   t_edsamyn *edyn);
/* after calling this routine the pointer to the data-structre t_edpar *edi is invalid, all memory is deallocated. No further calls to any routine of this interface possible */

void init_edsam(FILE *log,t_topology *top,t_inputrec *ir,
		t_mdatoms *md,int start,int homenr,t_commrec *cr,
	   t_edsamyn *edyn);
/* main init_routine - calls init_edi in a loop for every .edi-cycle contained in the input file, 
creates a NULL terminated list of t_edpar structures*/

void init_edi(FILE *log,t_topology *top,t_inputrec *ir,
	      t_mdatoms *md,int start,int homenr,t_commrec *cr,
	      t_edsamyn *edyn,t_edpar *edi);
/* init-routine called for every *.edi-cycle, initialises t_edpar structure */
 
void init_flood(t_edpar *edi, real dt);
/* called by init_edi, configure some flooding related variables and structures, print headers to output files */

void rmrotfit(int ned,rvec *x,matrix rotmat);  /* remove rotation */
void rmtransfit(int ned, rvec *x,rvec *transvec); /* remove fit */
void transfit(int ned, rvec *x, rvec *transvec); /* in fit direction */
void rotate_x(int ned, rvec *x, matrix rotmat); /* apply rotation of rotmat */

void flood_project(rvec *x, t_edpar *edi);
/* project fitted structure onto supbspace -> store in edi->vec.flood.xproj */

real flood_energy(t_edpar *edi);
/* from flood.xproj compute the Vfl(x) at this point -> store in edi->flood.Vfl*/

void flood_forces(t_edpar *edi);
/* from the position and the Vfl compute forces in subspace -> stored in edi->vec.flood.fproj */

void flood_blowup(t_edpar *edi,  rvec *forces_cart);
/* raise forces from subspace into cartesian space */

void update_adaption(t_edpar *edi);
/* update the values of Efl, deltaF depending on tau and Vfl */

void get_flood_enx_names(t_edpar *edi,char**, int *nnames); /* get header of energies */
void get_flood_energies(t_edpar *edi, real Vfl[],int nnames); /*fl has to be big enough to capture nnames-many entries*/
real get_rmsd(t_edpar *edi, rvec *x);


/* print the structure to the stream, used for debug purposes only */
void dump_rotmat(FILE*,matrix rotmat);
void dump_mat(FILE*, int dim, double **mat);
void dump_rvec(FILE *out, int dim, rvec *x);
void dump_edi(FILE*, t_edpar *edi);

/* read the parameter from .edi file and check for its label --> fatal_error if wrong*/
int read_checked_edint(FILE *file,char *label);
real read_checked_edreal(FILE *file,char *label);

 int read_edi_file(t_edsamyn *edyn, t_edpar *edi, int nr_mdatoms);

 int  read_edint(FILE *file,bool *bEOF);
 int  read_edint2(FILE *file);
 real read_edreal(FILE *file);
 void read_edx(FILE *file,int number,int *anrs,rvec *x);

 void read_edvecs(FILE *in,int nr,t_edvecs *vecs);
  /* calls read_edvec for the vector groups, only for flooding there is an extra call */
 void read_edvec(FILE *in,int nr,t_eigvec *tvec);
 void scan_edvec(FILE *in,int nr,rvec *vec);
 /*that is the (very private) scanf routine called by read_edvec */

 void fitit(int nr, rvec *x,t_edpar *edi,rvec *transvec,
		       matrix rmat);
 /* fits x[edi->fitnrs[i]] bzw x[edi->masnrs[i]] which should be equivalent
    X has to contain ALL Atoms */
 void do_edfit(int natoms,rvec *xp,rvec *x,matrix R,t_edpar *edi); 
 void put_in_origin(int nr,rvec *x,int nmass,int *masnrs,
			       real *mass,real tmass);

void project(rvec *x,t_edpar *edi,char *mode);
/* wrapper: project vector x onto all edi->vecs (mon, linfix,...) */

void project_to_eigvectors(rvec *x, t_eigvec *vec, t_edpar *edi,char *mode);
 /* project vector x, uses atoms sav.anrs[i] of x,
  if mode  is "x" it subtract average positions prior to projection
  and add them afterwards to retain the unchanged vector x
  mode = "x","v","f" -> store in xproj, vproj or fproj
  XXX mass-weighting is applied
 */

inline real projectx(t_edpar *edi,rvec *x,rvec *vec);
   /* does not subtract average positions, projection on single eigenvector is returned
   	used by: do_linfix, do_linacc, do_radfix, do_radacc, do_radcon
	here average position is subtracted in ed_cons prior to call to projectx
     */

inline real projectf(t_edpar *edi,rvec *x,rvec *vec);
   /* same as projectx, but mass-weighting is applied differently -> forces  */

void rad_project(t_edpar *edi,rvec *x,t_eigvec *vec);
   /* specialized:   projection is stored in vec->refproj
     ---> used for radacc, radfix, radcon  and center of flooding potential
   
     subtracts average positions, projects vector x, uses atoms sav.anrs[i] of x,
    */

 real calc_radius(t_eigvec *vec);
 void rmfit(int ned,rvec *x,rvec *transvec,matrix rotmat);
 void rotate_vec(int nr,rvec *x,matrix rotmat);
 void ed_cons(rvec *x,t_edpar *edi,int step);
   /* applies all the following  constraints*/

 void do_linfix(rvec *x,t_edpar *edi,int step);
 void do_linacc(rvec *x,t_edpar *edi);
 void do_radfix(rvec *x,t_edpar *edi,int step);
 void do_radacc(rvec *x,t_edpar *edi);
 void do_radcon(rvec *x,t_edpar *edi);

 void write_edo(t_edpar *edi,int step,real rmsd);
 void write_proj(FILE *out,t_edpar *edi,char *mode);
 void do_write_proj(FILE *out,t_eigvec *vec,char *mode);
 void write_edidx(FILE *out,t_edpar *edi);




/* definition of local_buffer structure, 
   the actual structure of the local_buffer of that respective function is defined directly in the function 
*/


struct t_ed_local {
    struct t_fitit *                fitit;
    struct t_do_edfit *             do_edfit;
    struct t_remove_pbc_effect *    remove_pbc_effect;
    struct t_do_edsam *             do_edsam;
    struct t_do_radcon *            do_radcon;
};



static void 
free_t_fitit(struct t_fitit *);

static void 
free_t_do_edfit(struct t_do_edfit *);

static void 
free_t_remove_pbc_effect(struct t_remove_pbc_effect *);

static void 
free_t_do_edsam(struct t_do_edsam *);

static void
free_t_do_radcon(struct t_do_radcon *);




struct t_ed_local *
init_local(void) 
{
    struct t_ed_local *local; 

    snew(local,1); 
    local->fitit=NULL;
    local->do_edfit=NULL;
    local->remove_pbc_effect=NULL;
    local->do_edsam=NULL;
    local->do_radcon=NULL;

    return local;
}



void free_local(struct t_ed_local *local) {
    free_t_fitit(local->fitit);
    free_t_do_edfit(local->do_edfit);
    free_t_remove_pbc_effect(local->remove_pbc_effect);
    free_t_do_edsam(local->do_edsam);
    free_t_do_radcon(local->do_radcon);
    sfree(local);
}


static inline void rvecsub(int dim,rvec *a, rvec *b, rvec *c) {
  /* c=a-b; */
  int i;
  for (i=0;i<dim;i++) {
    rvec_sub(a[i],b[i],c[i]);
  }
}
  
static inline void rvecadd(int dim,rvec *a, rvec *b, rvec *c) {
  /* c=a+b; */
  int i;
  for (i=0;i<dim;i++) {
    rvec_add(a[i],b[i],c[i]);
  }
}
static inline real rvecnorm(int dim,rvec *a) {
  /* c=a+b; */
  int i;
  real sum=0;
  for (i=0;i<dim;i++) {
    sum+=norm2(a[i]);
  }
  return sqrt(sum);
}

static inline void rvecsmul(int dim, real s, rvec *a) {
  int i;
  for (i=0;i<dim;i++) {
    svmul(s,a[i],a[i]);
  }
} 

static inline void rvec_to_one(int dim, rvec *a) {
  rvecsmul(dim,1.0/rvecnorm(dim,a),a);
}


static inline void rveccopy(int dim, rvec *a, rvec *b) {
  /*b=a;*/
  int i;
  for (i=0;i<dim;i++) {
    copy_rvec(a[i],b[i]);
  }
} 



int read_edi(FILE* in, t_edsamyn *edyn,t_edpar *edi,int nr_mdatoms, int edi_nr);





/******************************* IMPLEMENTATION ******************************************/

/* ---------------------- DEBUG Helpers ------------------------*/

void dump_edi(FILE* out, t_edpar *edpars) {
int i;
    fprintf(out,"#NINI\n %d\n#SELMAS\n %d\n#ANALYSIS_MAS\n %d\n",
       edpars->nini,edpars->fitmas,edpars->pcamas);
    fprintf(out,"#OUTFRQ\n %d\n#LOGFRQ\n %d\n#MAXLEN\n %d\n#SLOPECRIT\n %f\n",
        edpars->outfrq,edpars->logfrq,edpars->maxedsteps,edpars->slope);
    fprintf(out,"#PRESTEPS\n %d\n#DELTA_F0\n %f\n#TAU\n %f\n#EFL_NULL\n %f\n#ALPHA2\n %f\n",
        edpars->presteps,edpars->flood.deltaF0,edpars->flood.tau,edpars->flood.constEfl,edpars->flood.alpha2);
    fprintf(out,"REFERENCE: nr %d indices : ",edpars->sref.nr);
for (i=0; i<edpars->sref.nr; i++) 
   fprintf(out," %d",edpars->sref.anrs[i]);
fprintf(out,"\nCoordinates: \n");
dump_rvec(out,edpars->sref.nr,edpars->sref.x);
}

void dump_rotmat(FILE* out,matrix rotmat) {
   fprintf(out,"MATRIX: %f %f %f\n",rotmat[XX][XX],rotmat[XX][YY],rotmat[XX][ZZ]);
   fprintf(out,"ROTMAT: %f %f %f\n",rotmat[YY][XX],rotmat[YY][YY],rotmat[YY][ZZ]);
   fprintf(out,"ROTMAT: %f %f %f\n",rotmat[XX][ZZ],rotmat[XX][ZZ],rotmat[XX][ZZ]);
}

void dump_rvec(FILE *out, int dim, rvec *x) {
int i;
   for (i=0; i<dim;i++)
      fprintf(out,"dim i: %f %f %f\n",x[i][XX],x[i][YY],x[i][ZZ]);
}

void dump_mat(FILE* out, int dim, double** mat) {
int i,j;
fprintf(out,"MATRIX:\n");
for (i=0;i<dim;i++) {
   for (j=0;j<dim;j++)
      fprintf(out,"%f ",mat[i][j]);
   fprintf(out,"\n");
   }
}

/**********************************************************************************
******************** FLOODING *****************************************************
**********************************************************************************/

void write_edo_flood(t_edpar *edi, int step) {
  int i;
  fprintf(edi->edo,"%d.th FL: %d %g %g %g\n",edi->flood.flood_id,step, edi->flood.Efl, edi->flood.Vfl, edi->flood.deltaF);
  fprintf(edi->edo,"FL_FORCES: ");
  for (i=0;i<edi->flood.vecs.neig;i++)
  fprintf(edi->edo," %f",edi->flood.vecs.fproj[i]);
  fprintf(edi->edo,"\n");
  fflush(edi->edo);
}

void flood_project(rvec *x, t_edpar *edi)
{
  /* projects the positions onto the subspace */
  int i;

  /* do projection */
  project_to_eigvectors(x,&edi->flood.vecs,edi,"x");
}  


real flood_energy(t_edpar *edi) {
  /* compute flooding energy Vfl
     Vfl = Efl * exp( - \frac {kT} {2Efl alpha^2} * sum_i { \lambda_i c_i^2 } )
     \lambda_i is the reciproce eigenvalue 1/\sigma_i
         it is already computed by make_edi and stored in stpsz[i]
     bHarmonic:
       Vfl = - Efl * 1/2(sum _i {\frac 1{\lambda_i} c_i^2})
  */
  real summe;
  int i;
  
  summe=0.0;
  /*compute sum which will be the exponent of the exponential */
  if (edi->flood.bHarmonic)
    for (i=0;i<edi->flood.vecs.neig; i++)
      summe+=edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i])*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]);
  else
    for (i=0;i<edi->flood.vecs.neig; i++)
      summe+=edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i])*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]);
#ifdef DUMP_FORCES
  fprintf(logfile, "REFPROJ: ");
  for (i=0;i<edi->flood.vecs.neig; i++)
    fprintf(logfile, "%f ",edi->flood.vecs.refproj[i]);
  fprintf(logfile, "\n");

  fprintf(logfile, "XPROJ: ");
  for (i=0;i<edi->flood.vecs.neig; i++)
    fprintf(logfile, "%f ",edi->flood.vecs.xproj[i]);
  fprintf(logfile, "\n");
  
  fprintf(logfile, "SUMME: %f kT : %f alpha2 : %f Efl %f\n ", summe, edi->flood.kT, edi->flood.alpha2, edi->flood.Efl);
#endif
  
  /*compute the gauss function*/
  if (edi->flood.bHarmonic)
    edi->flood.Vfl=  - 0.5*edi->flood.Efl*summe;  /* minus sign because Efl is negativ, if restrain is on. */
  else
    edi->flood.Vfl= edi->flood.Efl!=0 ? edi->flood.Efl*exp(-edi->flood.kT/2/edi->flood.Efl/edi->flood.alpha2*summe) :0;
  return edi->flood.Vfl;
}

void flood_forces(t_edpar *edi) {
/* compute the forces in the subspace of the flooding eigenvectors
   by the formula F_i= V_{fl}(c) * ( \frac {kT} {E_{fl}} \lambda_i c_i */
  int i;
  real energy=edi->flood.Vfl;
#ifdef DUMP_FORCES
 fprintf(logfile, "Vfl= %f, Efl= %f, xproj= %f, refproj= %f\n",energy, edi->flood.Efl,edi->flood.vecs.xproj[0],edi->flood.vecs.refproj[0]);
#endif
  if (edi->flood.bHarmonic)
    for (i=0; i<edi->flood.vecs.neig; i++) {
      edi->flood.vecs.fproj[i]= edi->flood.Efl* edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]);
#ifdef DUMP_FORCES
    fprintf(logfile, "%f ",edi->flood.vecs.fproj[i]);
#endif
  }

  else
    for (i=0; i<edi->flood.vecs.neig; i++) {
      /* if Efl is zero the forces are zero if not use the formula */
      edi->flood.vecs.fproj[i]= edi->flood.Efl!=0 ? edi->flood.kT/edi->flood.Efl/edi->flood.alpha2*energy*edi->flood.vecs.stpsz[i]*(edi->flood.vecs.xproj[i]-edi->flood.vecs.refproj[i]) : 0;
#ifdef DUMP_FORCES
    fprintf(logfile, "force %f ",edi->flood.vecs.fproj[i]);
#endif
  }
#ifdef DUMP_FORCES
  fprintf(logfile,"\n");
#endif
}

void flood_blowup(t_edpar *edi, rvec *forces_cart) {
  /* this function lifts the forces from the subspace to the cartesian space
     all the values not contained in the subspace are assumed to be zero and then 
     a coordinate transformation from eigenvector to cartesian vectors is performed 
     The nonexistent values don't have to be set to zero explicitly, they would occur 
     as zero valued summands, hence we just stop to compute this part of the sum.
     
     for every atom we add all the contributions to this atom from all the different eigenvectors.
     
     NOTE: one could add directly to the forcefield forces, would mean we wouldn't have to clear the 
     field forces_cart prior the computation, but momentarily we want to compute the forces seperately 
     to have them accessible for diagnostics
  */
  int i,j,eig;
  rvec dum;
  real *forces_sub;
  forces_sub=edi->flood.vecs.fproj;
  /* clear forces first */
  for (j=0; j<edi->ned; j++) 
    clear_rvec(forces_cart[j]);
#if (DEBUG==BLOWUP)
  fprintf(stderr,"cleared cartesian force vector:");
  dump_rvec(stderr, edi->ned, forces_cart);
#endif
  /* now compute atomwise */
  for (j=0; j<edi->sav.nr; j++) {  /* should this be sav.nr ??? */
    /* compute    forces_cart[edi->sav.anrs[j]] */
    for (eig=0; eig<edi->flood.vecs.neig; eig++) {
      /* force vector is force * eigenvector compute only atom j */
      svmul(forces_sub[eig],edi->flood.vecs.vec[eig][j],dum);
      /* add this vector to the cartesian forces */
      rvec_inc(forces_cart[edi->sav.anrs[j]],dum);
    }
#ifdef DUMP_FORCES
    fprintf(logfile,"%d %f %f %f\n",
	    edi->sav.anrs[j], forces_cart[edi->sav.anrs[j]][XX],forces_cart[edi->sav.anrs[j]][YY],forces_cart[edi->sav.anrs[j]][ZZ]); 
  } 
  fprintf(logfile,"\n--------------------------------\n");
#if 0
  {
#endif
#else
}
#endif
}


void update_adaption(t_edpar *edi) {
/* this function updates the parameter Efl and deltaF according to the rules given in 
   'predicting unimolecular chemical reactions: chemical flooding' M Mueller et al, J. chem Phys.
*/

  if ((edi->flood.tau < 0 ? -edi->flood.tau : edi->flood.tau )>0.00000001) {
    edi->flood.Efl=edi->flood.Efl+edi->flood.dt/edi->flood.tau*(edi->flood.deltaF0-edi->flood.deltaF);
    /* check if restrain (inverted flooding) --> don't let EFL become positiv*/
    if (edi->flood.alpha2<0 && edi->flood.Efl>-0.00000001)
      edi->flood.Efl=0;
    
    edi->flood.deltaF=(1-edi->flood.dt/edi->flood.tau)*edi->flood.deltaF+edi->flood.dt/edi->flood.tau*edi->flood.Vfl;
  };
}


void do_single_flood(FILE *log,rvec x_orig[],rvec force[], t_edpar *edi, int step, int nr) {
  
  int i,j,ned=edi->ned,iupdate=500;
   matrix rotmat;
  real mas,rad;
  t_edpar *actual_edi;

 
  
  /*make copy of coordinates */
  
  for (i=0;i<edi->ned;i++) {
    copy_rvec(x_orig[i],edi->flood.loc.x[i]);
    if (!finite(edi->flood.loc.x[i][XX])) {
      fprintf(stderr,"Found invalid coordinate: \n");
      dump_rvec(stderr,edi->ned,edi->flood.loc.x);
      gmx_fatal(FARGS,"Something is wrong with the coordinates: a LINCS error? to much flooding strength? wrong flooding vectors?");
    }
  }



  /* fit the structure */
  DEBUG_PRINT("fit the structure...");
#ifdef DEBUG
  fprintf(stderr,"the structure to be fitted:\n");
  dump_rvec(stderr,edi->ned,edi->flood.loc.x);
#endif
#ifdef DOFIT
  fitit(ned,edi->flood.loc.x,edi,edi->flood.loc.transvec,rotmat);
#endif
#ifdef DEBUG
    dump_rotmat(stderr,rotmat);
#endif
  /* put projected values into edi->flood.vecs.xproj */
  DEBUG_PRINT("put projected values into edi->flood.vecs.xproj");
  flood_project(edi->flood.loc.x,edi);

  DEBUG_PRINT("compute_energy");
  flood_energy(edi);
  update_adaption(edi);

  
  DEBUG_PRINT("compute flooding forces");
  flood_forces(edi);
  
  /* translate them into cartesian coordinates */
  flood_blowup(edi, edi->flood.loc.forces_cartesian);

  /* rotate forces back so that they correspond to the given structure and not to the fitted one */
#ifdef DOFIT
  rmrotfit(ned, edi->flood.loc.forces_cartesian, rotmat);
#endif
  /* and finally: add forces to the master force variable */
  for (i=0; i<edi->ned; i++)
    rvec_inc(force[i],edi->flood.loc.forces_cartesian[i]);
  

  if (do_per_step(step,edi->outfrq)) {
    write_edo_flood(edi,step);
  }
  
}



void do_flood(FILE *log, t_commrec *cr, rvec x_orig[],rvec force[], t_edsamyn *edyn, int step)
{
  /* this is the hook, 
     do_flood is called by mdrun via the force() (force.c) routine.
     the parameters are
     logfile, commrec (to no if we are on the master node), x_orig (positions), 
     force ( forcefield forces , we add the flooding force to them), edi - all the parameters)
  */
  int i,j,ned,iupdate=500;
  int nr;
  matrix rotmat;
  real mas,rad;
  t_edpar *edi;
  t_edpar *actual_edi;
  if (!edyn) 
    return;
  if (!edyn->bEdsam || !edyn->edpar)
    return;

  edi = edyn->edpar;
  if (!edi->flood.vecs.neig) return;
  ned=edi->ned;
  nr=0;
  if(!MASTER(cr))
    return;
  actual_edi=edi;
  while (actual_edi) {
    DEBUG_PRINT("call flooding for one matrix");
    do_single_flood(log,x_orig,force,actual_edi,step,nr++);
    actual_edi=actual_edi->next_edi;
  }
}


void init_flood(t_edpar *edi, real dt) {
  int i;
  matrix rotmat;
  /*  edi->flood.deltaF=0; wird jetzt in edi-Datei initialisiert! */
  edi->flood.Efl=edi->flood.constEfl;
  edi->flood.Vfl=0;
  /*  edi->flood.kT=2.5; */
  edi->flood.dt=dt;
  if (edi->flood.vecs.neig) {
#ifdef DUMP_FORCES
    logfile=ffopen("dump_force2.log","w");
    fprintf(logfile,"Vfl Efl\nForce in Subspace\nForces in cartesian space\n");
#endif
    fprintf(edi->edo,"FL_HEADER: Flooding of matrix %d is switched on! The flooding output will have the following format:\n",edi->flood.flood_id);
    if (edi->flood.flood_id<1)
      fprintf(edi->edo,"FL_HEADER: Step Efl Vfl deltaF \n");
    /* get memory only once, that should speed up computations */

    snew(edi->flood.loc.x,edi->ned);
    snew(edi->flood.loc.transvec,edi->ned);
    snew(edi->flood.loc.forces_cartesian,edi->ned);

    /* set center of flooding potential */
    if (edi->sori.nr > 0) {
      /* ... to the structure given with -ori */
      fitit(edi->ned,edi->sori.x,edi,edi->flood.loc.transvec,rotmat);
      rad_project(edi, edi->sori.x, &edi->flood.vecs);
    } else
      /* ... to the center of the covariance matrix, i.e. the average structure, i.e. zero in the projected system */
      for (i=0;i<edi->flood.vecs.neig;i++)
	edi->flood.vecs.refproj[i]=0.0;
    }
}


/*********** Energy book keeping ******/
void get_flood_enx_names(t_edpar *edi, char** names, int *nnames)  /* get header of energies */
{
  t_edpar *actual;
  int count;
  char buf[STRLEN];
  actual=edi;
  count = 1;
  while (actual) {
    srenew(names,count);
    sprintf(buf,"Vfl_%d",count);
    names[count-1]=strdup(buf);
    actual=actual->next_edi;
    count++;
  };
  *nnames=count-1;
    
}


void get_flood_energies(t_edpar *edi, real Vfl[],int nnames) {
/*fl has to be big enough to capture nnames-many entries*/
  t_edpar *actual;
  int count;
  char buf[STRLEN];
  actual=edi;
  count = 1;
  while (actual) {
    Vfl[count-1]=actual->flood.Vfl;
    actual=actual->next_edi;
    count++;
  };
  if (nnames!=count-1) 
    gmx_fatal(FARGS,"Number of energies is not consistent with t_edi structure");
    
}


/************* END of FLOODING IMPLEMENTATION **************************/

/* ------------------ EDSAM IMPLEMENTATION ---------------------------------- */
/* of course many edsam routines are used or necessary for flooding, too */

void ed_open(int nfile,t_filenm fnm[],t_edsamyn *edyn,t_commrec *cr)
{
  if(!MASTER(cr))
    return;
  fprintf(stderr,"ED sampling will be performed!\n");
  edyn->bEdsam=TRUE;
  edyn->edinam=ftp2fn(efEDI,nfile,fnm);
  edyn->edonam=ftp2fn(efEDO,nfile,fnm); 
}



void init_edi(FILE *log,t_topology *top,t_inputrec *ir,
	      t_mdatoms *md,int start,int homenr,t_commrec *cr,
	      t_edsamyn *edyn,t_edpar *edi) {
  int i,j,ned,*refmasnrs;
  rvec *xdum,*transvec;
 
  matrix rotmat;
  t_edpar *an_edi;


  edi->bNeedDoEdsam=edi->vecs.mon.neig || edi->vecs.linfix.neig || edi->vecs.linacc.neig || edi->vecs.radfix.neig || edi->vecs.radacc.neig || edi->vecs.radcon.neig;
  ned=edi->ned;
  if (start+homenr<ned)
    gmx_fatal(FARGS,"ED sampling currently only works in parallel when all atoms\n"
		"involved in ED constraints are on the master processor - sorry.");
  
  fprintf(log,"Initialising ED sampling: start=%d homenr=%d ned=%d\n\n",start,
	  homenr,ned);
  
  /* evaluate masses */

  edi->tmass=0.0;
  edi->nmass=edi->sref.nr;
  snew(edi->mass,edi->nmass);
  snew(edi->masnrs,edi->nmass);
  snew(refmasnrs,edi->nmass);
  for(i=0; (i < edi->nmass); i++) {
    if (edi->fitmas) {
      edi->mass[i]=top->atoms.atom[edi->sref.anrs[i]].m;
    } else {
      edi->mass[i]=1.0;
    }
    edi->masnrs[i]=edi->sref.anrs[i];
    refmasnrs[i]=i;
    edi->tmass+=edi->mass[i];
  }
  snew(edi->sav.sqrtm,edi->sav.nr);
  for (i=0;i<edi->sav.nr;i++) {
    if (edi->pcamas)
      edi->sav.sqrtm[i]=sqrt(top->atoms.atom[edi->sav.anrs[i]].m);
    else
      edi->sav.sqrtm[i]=1.0;
  }
  /* mark atoms that are to be used for rotational fit */
  edi->nfit = edi->sref.nr;
  snew(edi->fitnrs,edi->nfit);
  for(i=0; (i < edi->nfit); i++) 
    edi->fitnrs[i] = edi->sref.anrs[i];
  
  /* put reference structure in origin */
  put_in_origin(edi->sref.nr,edi->sref.x,edi->nmass,refmasnrs,edi->mass,
		edi->tmass);

  /* init flooding parameters */
  init_flood(edi,ir->delta_t);
  edi->local=init_local();
#ifdef DEBUG
  dump_edi(stderr,edi);
#endif
  sfree(refmasnrs);
}

void init_edsam(FILE *log,t_topology *top,t_inputrec *ir,
                t_mdatoms *md,int start,int homenr,t_commrec *cr,
                t_edsamyn *edyn)
{
    int i,j,ned;
    rvec *xdum,*transvec;
    matrix rotmat;
    t_edpar *an_edi;
    if(!MASTER(cr))
        return;
    if (edyn->bEdsam) {
        fprintf(stderr,"init_edsam\n");
        /* first read the input. All input is stored in edi */
        snew(edyn->edpar,1);
        read_edi_file(edyn,edyn->edpar,top->atoms.nr);
        an_edi=edyn->edpar;
        edyn->edpar->edo=ffopen(edyn->edonam,"w");
        while(an_edi!=NULL) {
            an_edi->edo=edyn->edpar->edo;
            init_edi(log,top,ir,md,start,homenr,cr,edyn,an_edi);
            an_edi=an_edi->next_edi;
        };  
    } else
        edyn->edpar=NULL;
}

void do_first_edsam(FILE *log,t_topology *top,
                    t_mdatoms *md,int start,int homenr,t_commrec *cr,
                    rvec x[],matrix box, t_edsamyn *edyn,bool bHaveConstr) 
{
    
    int i,j,ned;
    rvec *xdum,*transvec;
    matrix rotmat;
    t_edpar *edi;
    
    if(!MASTER(cr) || !edyn->bEdsam)
        return;
    if (!edyn->edpar) 
        return;
    edi=edyn->edpar;
    if (edi->bNeedDoEdsam && ( bHaveConstr || ed_constraints(edyn)) )
        snew(edi->x_unc,edi->ned);
    else
        edi->x_unc=NULL;
    
    fprintf(stderr,"Do first edsam...\n");
    ned=edi->ned;
    
    /* remove pbc */
    snew(xdum,top->atoms.nr);
    rm_pbc(&(top->idef),top->atoms.nr,box,x,xdum);

    /* fit starting positions to reference structure */
    snew(transvec,ned);
    fitit(ned,xdum,edi,transvec,rotmat);
    fprintf(log,"Initial RMSD from reference structure = %10.5f nm\n\n",get_rmsd(edi,xdum));
    sfree(transvec);

    /* calculate initial projections */
    project(xdum,edi,"x");
    fprintf(log,"Initial projections:\n");
    write_edidx(log,edi);

    /* process target structure, if required */
    if (edi->star.nr > 0) {
        snew(transvec,ned);
        fitit(ned,edi->star.x,edi,transvec,rotmat);
        rad_project(edi,edi->star.x,&edi->vecs.radcon);
        sfree(transvec);
    } else {
        rad_project(edi,xdum,&edi->vecs.radcon);
    };
    
    
    /* process structure that will serve as origin of expansion circle */
    if (edi->sori.nr > 0) {
        snew(transvec,ned);
        fitit(ned,edi->sori.x,edi,transvec,rotmat);
        rad_project(edi,edi->sori.x,&edi->vecs.radacc);
        rad_project(edi,edi->sori.x,&edi->vecs.radfix);
        sfree(transvec);
    }
    else {
        rad_project(edi,xdum,&edi->vecs.radacc);
        rad_project(edi,xdum,&edi->vecs.radfix);
    }
    
    /* set starting projections for linsam */
    rad_project(edi,xdum,&edi->vecs.linacc);
    rad_project(edi,xdum,&edi->vecs.linfix);
    sfree(xdum);

    /* calculate initial radii */
    fprintf(log,"Initial fixed increment radius=%f\n",edi->vecs.radfix.radius);
    fprintf(log,"Initial   acceptance    radius=%f\n",edi->vecs.radacc.radius);
    fprintf(log,"Initial   contracting   radius=%f\n",edi->vecs.radcon.radius);
      
    write_edidx(edi->edo,edi);
}

int read_edi_file(t_edsamyn *edyn, t_edpar *edi, int nr_mdatoms) {
    FILE *in;
    t_edpar *actual_edi;
    t_edpar *edi_read;
    int readmagic;
    int edi_nr;
    
    in=ffopen(edyn->edinam,"r");  
    /* now read a sequence of edi files */
    actual_edi=edi;
    edi_nr=0;
    read_edi(in,edyn,actual_edi,nr_mdatoms,edi_nr++);
    if (edi->nini!=nr_mdatoms)
        gmx_fatal(FARGS,"The edsam/flooding input file %s was empty, willste mich verkackeiern?");
    snew(edi_read,1);
    while( read_edi(in, edyn, edi_read, nr_mdatoms, edi_nr++)) {
        actual_edi->next_edi=edi_read;
        actual_edi=edi_read;
        snew(edi_read,1);
    }
    sfree(edi_read);
    actual_edi->next_edi=NULL;
    ffclose(in);
    return 1;
}

int read_edi(FILE* in, t_edsamyn *edyn,t_edpar *edi,int nr_mdatoms, int edi_nr)
{
    int i,j,idum,readmagic;
    static const int magic=668;
    int ignore;
    rvec *xdum;
    bool bEOF;
    
    /* the edi file is not free format, so expect problems if the input
     * is corrupt.
     */
    
    /* check the magic number */
    readmagic=read_edint(in,&bEOF);
    if (bEOF)
        return 0;
    if (readmagic != magic) {
        if (readmagic==666 || readmagic==667)
            gmx_fatal(FARGS,"wrong magic number: Use newest version of make_edi to produce edi file");
        else
            gmx_fatal(FARGS,"Wrong magic number %d in %s",readmagic,edyn->edinam);
    }
    
    
    /* check the number of atoms */
    edi->nini=read_edint(in,&bEOF);
    if (edi->nini != nr_mdatoms)
        gmx_fatal(FARGS,"Nr of atoms in %s (%d) does not match nr of md atoms (%d)",
                  edyn->edinam,edi->nini,nr_mdatoms); 
    
    /* Done checking. For the rest we blindly trust the input */
    edi->fitmas=read_checked_edint(in,"FITMAS");
    edi->pcamas=read_checked_edint(in,"ANALYSIS_MAS");
    edi->outfrq=read_checked_edint(in,"OUTFRQ");
    edi->logfrq=read_checked_edint(in,"LOGFRQ");
    edi->maxedsteps=read_checked_edint(in,"MAXLEN");
    edi->slope=read_checked_edreal(in,"SLOPECRIT");
    
    edi->presteps=read_checked_edint(in,"PRESTEPS");
    edi->flood.deltaF0=read_checked_edreal(in,"DELTA_F0");
    edi->flood.deltaF=read_checked_edreal(in,"INIT_DELTA_F");
    edi->flood.tau=read_checked_edreal(in,"TAU");
    edi->flood.constEfl=read_checked_edreal(in,"EFL_NULL");
    edi->flood.alpha2=read_checked_edreal(in,"ALPHA2");
    edi->flood.kT=read_checked_edreal(in,"KT");
    edi->flood.bHarmonic=read_checked_edint(in,"HARMONIC");
    edi->flood.flood_id=edi_nr;
    edi->sref.nr=read_checked_edint(in,"NREF");
    /*  fprintf(stderr,"read ref\n"); */
    /* allocate space for reference positions and read them */
    snew(edi->sref.anrs,edi->sref.nr);
    snew(edi->sref.x,edi->sref.nr);
    edi->sref.sqrtm=NULL;
    read_edx(in,edi->sref.nr,edi->sref.anrs,edi->sref.x);
    
    /*  fprintf(stderr,"read aver\n"); */
    /* average positions. they define which atoms will be used for ED sampling */
    edi->sav.nr=read_checked_edint(in,"NAV");
    snew(edi->sav.anrs,edi->sav.nr);
    snew(edi->sav.x,edi->sav.nr);
    read_edx(in,edi->sav.nr,edi->sav.anrs,edi->sav.x);
    
    edi->ned=edi->sref.anrs[edi->sref.nr-1]+1;
    if (edi->sav.anrs[edi->sav.nr-1] > edi->ned)
        edi->ned=edi->sav.anrs[edi->sav.nr-1]+1;
    /* fprintf(stderr,"Nr of atoms for ED sampling buffer: %d\n",edi->ned);*/
    
    /* eigenvectors */
    /*  fprintf(stderr,"read edvecs\n"); */
    read_edvecs(in,edi->sav.nr,&edi->vecs);
    read_edvec(in,edi->sav.nr,&edi->flood.vecs);
    /*  fprintf(stderr,"read target_pos\n"); */
    /* target positions */
    edi->star.nr=read_edint(in,&bEOF);
    if (edi->star.nr > 0) {
        snew(edi->star.anrs,edi->star.nr);
        snew(xdum,edi->star.nr);
        edi->star.sqrtm=0;
        read_edx(in,edi->star.nr,edi->star.anrs,xdum);
        snew(edi->star.x,edi->ned);
        
        for(j=0; (j < edi->star.nr); j++) 
            if (edi->star.anrs[j] < 0 || edi->star.anrs[j] > edi->ned)
                gmx_fatal(FARGS,"ED sampling target index out of bounds: %d\n",edi->star.anrs[j]);
        for(i=0; (i < edi->ned); i++) {
            for(j=0; (j < edi->star.nr); j++) {
                if (edi->star.anrs[j] == i) {
                    copy_rvec(xdum[j],edi->star.x[i]);
                }
            }
        }
        sfree(xdum);
    }

    /* positions defining origin of expansion circle */
    edi->sori.nr=read_edint(in,&bEOF);
    if (edi->sori.nr > 0) {
        snew(edi->sori.anrs,edi->sori.nr);
        snew(xdum,edi->sori.nr);
        edi->sori.sqrtm=NULL;
        read_edx(in,edi->sori.nr,edi->sori.anrs,xdum);
        snew(edi->sori.x,edi->ned);
        
        for(j=0; (j < edi->sori.nr); j++) 
            if (edi->sori.anrs[j] < 0 || edi->sori.anrs[j] > edi->ned)
                gmx_fatal(FARGS,"ED sampling origin index out of bounds: %d\n",edi->sori.anrs[j]);
        for(i=0; (i < edi->ned); i++) {
            for(j=0; (j < edi->sori.nr); j++) {
                if (edi->sori.anrs[j] == i) {
                    copy_rvec(xdum[j],edi->sori.x[i]);
                }
            }
        }
        sfree(xdum);
    }
    
    /* all done */
    return 1;
    
}

void check(char *line, char *label) {
  if (!strstr(line,label)) 
        gmx_fatal(FARGS,"Could not find input parameter %s at expected position in edsam input-file (.edi)\nline read instead is %s",label,line);
}

int read_checked_edint(FILE *file,char *label) {
  char line[STRLEN+1];
  int idum;
  
  fgets2 (line,STRLEN,file);
  check(line,label);
  fgets2 (line,STRLEN,file);
  sscanf (line,"%d",&idum);
  return idum;
} 

int read_edint(FILE *file,bool *bEOF)
{
  char line[STRLEN+1];
  int idum;
  char *eof;

  eof=fgets2 (line,STRLEN,file);
  if (eof==NULL) {
    *bEOF = TRUE;
    return -1;
  }
  eof=fgets2 (line,STRLEN,file);  
  if (eof==NULL) {
    *bEOF = TRUE;
    return -1;
  }
  sscanf (line,"%d",&idum);
  *bEOF = FALSE;
  return idum;
}

real read_checked_edreal(FILE *file,char *label)
{
  char line[STRLEN+1];
  double rdum;

  fgets2 (line,STRLEN,file);
  check(line,label);
  fgets2 (line,STRLEN,file);
  sscanf (line,"%lf",&rdum);
  return (real) rdum; /* always read as double and convert to single */
}

real read_edreal(FILE *file)
{
  char line[STRLEN+1];
  double rdum;

  fgets2 (line,STRLEN,file);
  fgets2 (line,STRLEN,file);
  sscanf (line,"%lf",&rdum);
  return (real) rdum; /* always read as double and convert to single */
}

int read_edint2(FILE *file)
{
  char line[STRLEN+1];
  int idum;

  fgets2 (line,STRLEN,file);
  sscanf (line,"%d",&idum);
  return idum;
}

void read_edx(FILE *file,int number,int *anrs,rvec *x)
{
  int i,j;
  char line[STRLEN+1];
  double d[3];

  for(i=0; (i < number); i++) {
    fgets2 (line,STRLEN,file);
    sscanf (line,"%d%lf%lf%lf",&anrs[i],&d[0],&d[1],&d[2]);
    anrs[i]--; /* we are reading FORTRAN indices */
    for(j=0; (j < 3); j++)
      x[i][j]=d[j]; /* always read as double and convert to single */
  }
}

void read_edvecs(FILE *in,int nr,t_edvecs *vecs)
{
  read_edvec(in,nr,&vecs->mon);
  read_edvec(in,nr,&vecs->linfix);
  read_edvec(in,nr,&vecs->linacc);
  read_edvec(in,nr,&vecs->radfix);
  read_edvec(in,nr,&vecs->radacc);
  read_edvec(in,nr,&vecs->radcon);

}

void read_edvec(FILE *in,int nr,t_eigvec *tvec)
{
  int i,idum;
  double rdum;
  char line[STRLEN+1];

  tvec->neig=read_checked_edint(in,"NUMBER OF EIGENVECTORS");
  if (tvec->neig >0) {
    snew(tvec->ieig,tvec->neig);
    snew(tvec->stpsz,tvec->neig);
    snew(tvec->vec,tvec->neig);
    snew(tvec->xproj,tvec->neig);
    snew(tvec->vproj,tvec->neig);
    snew(tvec->fproj,tvec->neig);
    snew(tvec->refproj,tvec->neig);
    for(i=0; (i < tvec->neig); i++) {
      fgets2 (line,STRLEN,in);
      sscanf (line,"%d%lf",&idum,&rdum);
      tvec->ieig[i]=idum;
      tvec->stpsz[i]=rdum;
    }
    for(i=0; (i < tvec->neig); i++) {
      snew(tvec->vec[i],nr);
      scan_edvec(in,nr,tvec->vec[i]);
    }
  }
}

void scan_edvec(FILE *in,int nr,rvec *vec)
{
  char line[STRLEN+1];
  int i;
  double x,y,z;

  for(i=0; (i < nr); i++) {
    fgets2 (line,STRLEN,in);
    sscanf (line,"%le%le%le",&x,&y,&z);
    vec[i][XX]=x;
    vec[i][YY]=y;
    vec[i][ZZ]=z;
  }
}

void rotate_x(int nr,rvec *x,matrix rmat) {
  int i,j,k;
  rvec x_old;
  DEBUG_PRINT(" apply the rotation matrix \n");
  for(i=0;(i<nr);i++) {
    for(j=0;(j<3);j++)
      x_old[j]=x[i][j];
    for(j=0;(j<3);j++) {
      x[i][j]=0;
      for(k=0;(k<3);k++)
        x[i][j]+=rmat[k][j]*x_old[k];
    }
  }
}

struct t_fitit {
  rvec *xdum1;
  int nr1;
  rvec *xdum2;
  int nfit;
};

void free_t_fitit(struct t_fitit *p) { 
  sfree(p->xdum1);
  sfree(p->xdum2);
  sfree(p);
}
	  
void fitit(int nr, rvec *x,t_edpar *edi,rvec *transvec,matrix rmat)
{


  int i,j,k;
  bool bFirst;
  struct t_fitit *loc;
  
  if(edi->local->fitit != NULL)
  {
      bFirst = FALSE;
      loc    = edi->local->fitit;
  }
  else
  {
      bFirst = TRUE;
      snew(edi->local->fitit,1);
  }
  loc = edi->local->fitit;

  if (bFirst) {
    snew(loc->xdum1,nr);
    loc->nr1=nr;
  }
  if (nr>loc->nr1) {
    srenew(loc->xdum1,nr);
    loc->nr1=nr;
  };

  for(i=0; (i<nr); i++)
    copy_rvec(x[i],loc->xdum1[i]); 

  DEBUG_PRINT(" first do translational fit \n");
  put_in_origin(nr,x,edi->nmass,edi->masnrs,edi->mass,edi->tmass);

  /* determine transvec from difference after translational fit */
  for(i=0; (i<nr); i++) 
        rvec_sub(x[i],loc->xdum1[i],transvec[i]); 
 

  DEBUG_PRINT(" now rotational fit \n");
  if (bFirst) {
    loc->nfit = edi->nfit;
    snew(loc->xdum2,edi->nfit);
  }
  if (loc->nfit < edi->nfit) { /* happens in flooding with more than one matrix */
    srenew(loc->xdum2,edi->nfit);
    loc->nfit=edi->nfit;
  }
  for(i=0; (i<edi->nfit); i++)
    copy_rvec(x[edi->fitnrs[i]],loc->xdum2[i]);
  DEBUG_PRINT("do_edfit..");
  do_edfit(edi->nfit,edi->sref.x,loc->xdum2,rmat,edi);
  rotate_x(nr,x,rmat);
}

real get_rmsd(t_edpar *edi, rvec *x) {
  /* fit has to be done previously */
  real rmsd;
  int i,j;

  DEBUG_PRINT(" calculate RMSD \n");
  rmsd=0.0;
  for(i=0; (i<edi->nfit); i++) {
    for(j=0; (j<DIM); j++)
      rmsd+=pow((edi->sref.x[i][j]-x[edi->fitnrs[i]][j]),2);
  }
  rmsd/=(real) edi->nfit;
  rmsd=sqrt(rmsd);
  /*fprintf(stderr,"RMSD to reference structure=%15.10f\n",rmsd);*/
  return rmsd;
}

struct t_do_edfit { 
  double **omega;
  double **om;
};

void free_t_do_edfit(struct t_do_edfit *p) {
  int i;
  for (i=0; i<2*DIM; i++) {
    sfree(p->omega[i]);
    sfree(p->om[i]);
  };
  sfree(p->omega);
  sfree(p->om);
  sfree(p);
}
    

void do_edfit(int natoms,rvec *xp,rvec *x,matrix R,t_edpar *edi)
{
  /* this is a copy of do_fit with some modifications */
  int    c,r,n,j,i,irot;
  double d[6],xnr,xpc;
  matrix vh,vk,u;
  /*
  matrix vh_d,vk_d;
  */
  int    index;
  real   max_d;

  struct t_do_edfit *loc;
  bool bFirst;
  
  if(edi->local->do_edfit != NULL)
  {
      bFirst = FALSE;
      loc    = edi->local->do_edfit;
  }
  else
  {
      bFirst = TRUE;
      snew(edi->local->do_edfit,1);
  }
  loc = edi->local->do_edfit;

  if (bFirst) {
    snew(loc->omega,2*DIM);
    snew(loc->om,2*DIM);
    for(i=0; i<2*DIM; i++) {
      snew(loc->omega[i],2*DIM);
      snew(loc->om[i],2*DIM);
    }
  }

  for(i=0;(i<6);i++) {
    d[i]=0;
    for(j=0;(j<6);j++) {
      loc->omega[i][j]=0;
      loc->om[i][j]=0;
    }
  }

  /*  fprintf(stderr,"calc matrix U natloc->oms=%d",natloc->oms); */

  /*calculate the matrix U*/
  clear_mat(u);
  for(n=0;(n<natoms);n++) {
    /*    fprintf(stderr,"ref + x interlaced: "); */
    for(c=0; (c<DIM); c++) {

      xpc=xp[n][c];
      for(r=0; (r<DIM); r++) {
	xnr=x[n][r];
	/*	fprintf(stderr," %f %f",xp[n][r],xnr) ; */
	u[c][r]+=xnr*xpc;
      }
    }
  }
  /*  fprintf(stderr,"\n"); */
  /*construct loc->omega*/
  /*loc->omega is symmetric -> loc->omega==loc->omega' */
  for(r=0;(r<6);r++)
    for(c=0;(c<=r);c++)
      if ((r>=3) && (c<3)) {
        loc->omega[r][c]=u[r-3][c];
        loc->omega[c][r]=u[r-3][c];
      }
      else {
        loc->omega[r][c]=0;
        loc->omega[c][r]=0;
      }

  /*determine h and k*/
  DEBUG_PRINT("call jacobi");

#ifdef DEBUG
{
  int i;
  dump_mat(stderr,2*DIM,loc->omega);
  for (i=0; i<6; i++)
    fprintf(stderr,"d[%d] = %f\n",i,d[i]);
}
#endif
  jacobi(loc->omega,6,d,loc->om,&irot);

  if (irot==0) {
    fprintf(stderr,"IROT=0\n");
  }

  index=0; /* For the cloc->ompiler only */

  for(j=0;(j<3);j++) {
    max_d=-1000;
    for(i=0;(i<6);i++)
      if (d[i]>max_d) {
        max_d=d[i];
        index=i;
      }
    d[index]=-10000;
    for(i=0;(i<3);i++) {
      vh[j][i]=M_SQRT2*loc->om[i][index];
      vk[j][i]=M_SQRT2*loc->om[i+DIM][index];
    }
  }

  /*determine R*/
  for(c=0;(c<3);c++)
    for(r=0;(r<3);r++)
      R[c][r]=vk[0][r]*vh[0][c]+
              vk[1][r]*vh[1][c]+
              vk[2][r]*vh[2][c];
  if (det(R) < 0)
    for(c=0;(c<3);c++)
      for(r=0;(r<3);r++)
	R[c][r]=vk[0][r]*vh[0][c]+
	        vk[1][r]*vh[1][c]-
	        vk[2][r]*vh[2][c];
}

void put_in_origin(int nr,rvec *x,int nmass,int *masnrs,real *mass,real tmass)
{
  int i,j;
  rvec cm;

  /* calculate CM */
  for (i=0;(i<DIM);i++) 
    cm[i]=0.0;	 
  for (i=0;(i<nmass);i++) {
    for (j=0;(j<DIM);j++)
      cm[j]+=x[masnrs[i]][j]*mass[i];
  }

  for (i=0;(i<DIM);i++) 
    cm[i]/=tmass;

  /* and subtract it */
  for (i=0;(i<nr);i++) 
    rvec_dec(x[i],cm);
}

void project(rvec *x,t_edpar *edi,char *mode)
{
  int i;
  /* make the projections */
  /* it is not more work to subtract the average position in every subroutine again, 
     because these routines are rarely used simultanely */
  project_to_eigvectors(x,&edi->vecs.mon,edi,mode);
  project_to_eigvectors(x,&edi->vecs.linfix,edi,mode);
  project_to_eigvectors(x,&edi->vecs.linacc,edi,mode);
  project_to_eigvectors(x,&edi->vecs.radfix,edi,mode);
  project_to_eigvectors(x,&edi->vecs.radacc,edi,mode);
  project_to_eigvectors(x,&edi->vecs.radcon,edi,mode);
}


void project_to_eigvectors(rvec *x,t_eigvec *vec,t_edpar *edi,char *mode)
{
  int i,j,k;
  real proj;
  if (!vec->neig) return;
  /* subtract average positions */
  if (strcmp(mode,"x") == 0) {
    for (i=0;(i<edi->sav.nr);i++) 
      rvec_dec(x[edi->sav.anrs[i]],edi->sav.x[i]);
  }

  for (i=0;(i<vec->neig);i++) { 
    if (strcmp(mode,"x") == 0)  vec->xproj[i]=projectx(edi,x,vec->vec[i]);
    else if (strcmp(mode,"v") == 0) vec->vproj[i]=projectx(edi,x,vec->vec[i]);
    else if (strcmp(mode,"f") == 0) vec->fproj[i]=projectf(edi,x,vec->vec[i]);
    /* this has no influence on flooding forces, since this routine is called from the edsam branch
       completely disconnected from the do_flood branch */
  };

  /* add average positions */
  if (strcmp(mode,"x") == 0) {
  for (i=0;(i<edi->sav.nr);i++) 
    rvec_inc(x[edi->sav.anrs[i]],edi->sav.x[i]);
  }
}

void rad_project(t_edpar *edi,rvec *x,t_eigvec *vec)
{
  int i,j,k;
  real rad=0.0;

  /* subtract average positions */
  for (i=0;(i<edi->sav.nr);i++)
    rvec_dec(x[edi->sav.anrs[i]],edi->sav.x[i]);

  for (i=0;(i<vec->neig);i++) {
    vec->refproj[i]=projectx(edi,x,vec->vec[i]);
    rad+=pow((vec->refproj[i]-vec->xproj[i]),2);
  }
  vec->radius=sqrt(rad);
  
  /* add average positions */
  for (i=0;(i<edi->sav.nr);i++) 
    rvec_inc(x[edi->sav.anrs[i]],edi->sav.x[i]);
}

inline real projectx(t_edpar *edi,rvec *x,rvec *vec)
{
  int i,j;
  real proj=0.0;

  for (i=0;(i<edi->sav.nr);i++)
    for (j=0;(j<DIM);j++)
      proj+=vec[i][j]*x[edi->sav.anrs[i]][j]*edi->sav.sqrtm[i];

  return proj;
}

inline real projectf(t_edpar *edi,rvec *x,rvec *vec)
{
  int i,j;
  real proj=0.0;

  for (i=0;(i<edi->sav.nr);i++)
    for (j=0;(j<DIM);j++)
      proj+=vec[i][j]*x[edi->sav.anrs[i]][j]/edi->sav.sqrtm[i];

  return proj;
}

	
real calc_radius(t_eigvec *vec)
{
  int i;
  
  vec->radius = 0.0;
  for (i=0;(i<vec->neig);i++) 
    vec->radius+=pow((vec->refproj[i]-vec->xproj[i]),2);
  return vec->radius=sqrt(vec->radius);

}

void correct_vel(int nr,  rvec* vel, rvec *v1, rvec *v2, real fact, rvec* vcorr) {
    int i,m;
    real s;
    for (i=0; i<nr; i++)
        for (m=0; m<DIM; m++) /* x3-2*x2+x1/(dt^2*m) is force */ {
            s=vcorr[i][m]=fact*v2[i][m]-fact*v1[i][m];
            vel[i][m]+=s;
        };
}

void correct_force(int nr, rvec* force, rvec *v1, rvec *v2, real fact, t_topology *top) {
    int i,m;
    real mas;
    for (i=0; i<nr; i++)
        for (m=0; m<DIM; m++)  {
            mas=top->atoms.atom[i].m;
            force[i][m]+=(fact*v2[i][m]-fact*v1[i][m])*mas;
        }
}

struct t_remove_pbc_effect {
    t_pbc pbc;
};

void free_t_remove_pbc_effect(struct t_remove_pbc_effect *p) 
{
    /* there is nothing like free_pbc necessary (->pbc.h) */
    return;
} 


void remove_pbc_effect(int ned,rvec *transvec_compact, matrix box,t_edpar *edi) {
    
    bool bFirst;
    rvec null;
    int i;
    struct t_remove_pbc_effect* loc;
    
    if(edi->local->remove_pbc_effect != NULL)
    {
        bFirst = FALSE;
        loc    = edi->local->remove_pbc_effect;
    }
    else
    {
        bFirst = TRUE;
        snew(edi->local->remove_pbc_effect,1);
    }
    loc = edi->local->remove_pbc_effect;
    
    if (bFirst)
        set_pbc(&(loc->pbc),box);
    for (i=0;i<DIM;i++)
        null[i]=0.0;
    for (i=0;i<ned;i++)
        pbc_dx(&(loc->pbc),null,transvec_compact[i],transvec_compact[i]);
}

void prepare_edsam(int step, int start, int homenr, t_commrec *cr, rvec x[],t_edsamyn *edyn) {
    /* this is to produce the old x_unc, which is no longer provided by update() */
    int n;
    if (!edyn || !MASTER(cr)) 
        return;
    if (!edyn->bEdsam || !edyn->edpar)
        return;
    if (edyn->edpar->x_unc)
        for(n=0; n<edyn->edpar->ned; n++)
            copy_rvec(x[n],edyn->edpar->x_unc[n]);
}

struct t_do_edsam {
    rvec *transvec;
    rvec *vdum;
    rvec *fdum;
    matrix old_rotmat;
    real oldrad;
    rvec *old_transvec,*older_transvec,*transvec_compact;
    rvec *vcorr, *old_vcorr;
    rvec *xdum,*vdum2,*x;
    int cstep;
};

void free_t_do_edsam(struct t_do_edsam *p) {
    sfree(p->transvec);
    sfree(p->vdum);
    sfree(p->fdum);
    sfree(p->old_transvec);
    sfree(p->older_transvec);
    sfree(p->transvec_compact);
    sfree(p->vcorr);
    sfree(p->old_vcorr);
    sfree(p->xdum);
    sfree(p->vdum2);
    sfree(p->x);
    sfree(p);
}


void do_edsam(FILE *log,t_topology *top,t_inputrec *ir,int step,
              t_mdatoms *md,int start,int homenr, t_commrec *cr,
              rvec xs[],rvec xold[],rvec force[],matrix box,
              t_edsamyn *edyn,bool bHave_force)
{
    int i,j,ned,edstep=step,iupdate=500;
    
    matrix rotmat;
    real mas,rad;
    bool bFirst;
    real dt,dt_1,dt_2;
    struct t_do_edsam *loc;
    t_edpar *edi;
    dt = ir->delta_t;
    dt_1 = 1.0/dt;
    dt_2 = 1.0/(dt*dt);
    if (!edyn) 
        return;
    if (!edyn->bEdsam || !edyn->edpar)
        return;
    if(!MASTER(cr))
        return;
    edi=edyn->edpar;  
    if (!edi->bNeedDoEdsam) return;
    ned=edi->ned;
    
    if(edi->local->do_edsam != NULL)
    {
        bFirst = FALSE;
        loc    = edi->local->do_edsam;
    }
    else
    {
        bFirst = TRUE;
        snew(edi->local->do_edsam,1);
    }
    loc = edi->local->do_edsam;
    
    /* initialise radacc radius for slope criterion */
    if (bFirst) {
        loc->oldrad=calc_radius(&edi->vecs.radacc);
        snew(loc->vdum,ned);
        snew(loc->vdum2,ned);
        snew(loc->fdum,ned);
        snew(loc->xdum,ned);
        snew(loc->transvec,ned);
        snew(loc->old_transvec, ned);
        snew(loc->transvec_compact, ned);
        snew(loc->older_transvec, ned);
        snew(loc->vcorr,ned);
        snew(loc->old_vcorr,ned);
        if (!ed_constraints(edyn))
            snew(loc->x,ned); 
    }
    
    if (ed_constraints(edyn))
        loc->x=xs; /* if constraints are applied we want to work on the official copy of the coordinates themselves */
    else
        rveccopy(ned, xs, loc->x);  
    /* copy positions -- fitting and unfitting as previously programmed by Bert 
        introduced small errors in coordinates, this was due to unnecessary matrix inversion
        since this is reprogrammed no errors occur, but it is still slightly slower than copying 
        BUT if we want to use constraints we CANNOT COPY
        */
    
    /* calculate full conservative forces & velocities when required */
    if (bHave_force) {
        for(i=0; (i<ned); i++) {
            for(j=0; (j<DIM); j++)
                loc->vdum[i][j]=(loc->x[i][j]-xold[i][j])*dt_1;
        }
        for(i=0; (i<ned); i++) {
            mas=top->atoms.atom[i].m;
            if (edi->x_unc)
                for(j=0; (j<DIM); j++)
                    loc->fdum[i][j]=force[i][j]+(loc->x[i][j]-edi->x_unc[i][j])*dt_2*mas;
        }
        rveccopy(ned, loc->x, loc->xdum); /* get second copy of positons */
    }
    
    /* fit the structure */
    fitit(ned,loc->x,edi,loc->transvec,rotmat);
    
    rveccopy(ned, loc->transvec, loc->transvec_compact);
    remove_pbc_effect(ned,loc->transvec_compact,box,edi); /*this is to remove virtual jumps in translational velocity due to periodic boundaries*/
    if (bFirst) {
        copy_mat(rotmat, loc->old_rotmat);
        rveccopy(ned, loc->transvec_compact, loc->old_transvec);
        rveccopy(ned, loc->old_transvec,loc->older_transvec);
        dt_1=0; dt_2=0;
    };
    /* correct forces and velocities for the rotational and translational motion */
    if (bHave_force) {
        if (do_per_step(edstep,edi->outfrq)) {
            loc->cstep++;
            /* correct for translational motion */
            correct_vel(ned,loc->vdum,loc->transvec_compact,loc->old_transvec, -1.0*dt_1,loc->vcorr);
            correct_force(ned,loc->fdum,loc->transvec_compact,loc->old_transvec, -1.0*dt_2,top); /* correct forces to incorporate effect of corrected velocities */
            correct_force(ned,loc->fdum,loc->old_transvec,loc->older_transvec, 1.0*dt_2,top); /* opposite sign's here are correct */
            
            /*correct for rotational motion */
            rveccopy(ned, loc->vdum, loc->vdum2);
            rotate_vec(ned, loc->vdum2,loc->old_rotmat);   /* two copies of velocities rotated with different matrices */
            rotate_vec(ned, loc->vdum, rotmat);        /* difference is rotational motion in velocities -- > correct forces */
            
            rotate_vec(ned, loc->fdum, rotmat);       /* rotate forces before you correct them */
            correct_force(ned, loc->fdum, loc->vdum, loc->vdum2, -1.0*dt_1,top); /* correct with rotation in velocities */
        } /* up to here steps need only be done in outfrq steps */
      else
          rotate_vec(ned, loc->vdum, rotmat);        /* rotate velocities even in non - outfrq steps */
 
      rvecadd(ned, loc->xdum, loc->transvec, loc->xdum); /*move copy of actual structure into origin */
      rotate_vec(ned, loc->xdum, loc->old_rotmat); /* rotate positions with old matrix -> rotation in positions can be extracted */
      correct_vel(ned, loc->vdum,loc->x, loc->xdum, -1.0*dt_1,loc->vcorr);
 
      if (loc->cstep > 2 && do_per_step(edstep,edi->outfrq))
          correct_force(ned, loc->fdum, loc->vcorr, loc->old_vcorr, -1.0*dt_1, top);
      rveccopy(ned, loc->vcorr, loc->old_vcorr);
      copy_mat(rotmat, loc->old_rotmat);
      rveccopy(ned,loc->old_transvec,loc->older_transvec);
      rveccopy(ned,loc->transvec_compact,loc->old_transvec);
    } /* end of - correction of forces and velocities for rotational and translational motion */
  
    /* update radsam references, when required */
   if (do_per_step(edstep,edi->maxedsteps) && edstep > edi->presteps) {
       project(loc->x,edi,"x");
       rad_project(edi,loc->x,&edi->vecs.radacc);
       rad_project(edi,loc->x,&edi->vecs.radfix);
       loc->oldrad=-1.e5;
   }

/* update radacc references, when required */
   if (do_per_step(edstep,iupdate) && edstep > edi->presteps) {
       rad=calc_radius(&edi->vecs.radacc);
       if ((rad-loc->oldrad) < edi->slope) {
           project(loc->x,edi,"x");
           rad_project(edi,loc->x,&edi->vecs.radacc);
           loc->oldrad=0.0;
       }
       else
           loc->oldrad=rad;
   }

 
    /* apply the constraints */
       if (edstep> edi->presteps && ed_constraints(edyn))
           ed_cons(loc->x,edi,edstep);

    /* produce output, when required */
    if (do_per_step(edstep,edi->outfrq) && bHave_force) {
        /* rotate forces and velocities */
        project(loc->vdum,edi,"v");
        project(loc->fdum,edi,"f");
        project(loc->x,edi,"x");
        write_edo(edi,edstep,get_rmsd(edi,loc->x));
    }

    /* write to log, when required */
    if ((edstep > 0) && do_per_step(edstep,edi->logfrq)) {
        fprintf(log,"ED sampling information, step %d\n",edstep);
        project(loc->x,edi,"x");
        write_edidx(log,edi);
        fprintf(log,"acceptance radius = %f\n",
                calc_radius(&edi->vecs.radacc));
        fprintf(log,"fixed increment radius = %f\n",
                calc_radius(&edi->vecs.radfix));
        fprintf(log,"contracting radius = %f\n",
                calc_radius(&edi->vecs.radcon));
      /*  fflush(log); */
    }
  
    if (edstep == ir->nsteps) ffclose(edi->edo);
   
   /* remove fitting, if we do not work on throw away copy */
   if (ed_constraints(edyn))
       rmfit(ned,loc->x,loc->transvec,rotmat);
   
   if (bFirst) 
       bFirst = FALSE;
}

void rmrotfit(int ned,rvec *x,matrix rotmat)
{
    int i,j,k;
    matrix r_inv;
    rvec xdum;
#ifdef DEBUG
    dump_rotmat(stderr,rotmat);
#endif
    /* invert the rotation matrix and apply */
    for(i=0;(i<ned);i++) {
        for(j=0;(j<3);j++)
            xdum[j]=x[i][j];
        for(j=0;(j<3);j++) {
            x[i][j]=0;
            for(k=0;(k<3);k++)
                x[i][j]+=rotmat[j][k]*xdum[k];
        }
    }
}

void rmtransfit(int ned, rvec *x, rvec *transvec) { 
    int i;
    /* subtract the translation vector (this vector has only equal elements)*/
    for(i=0;(i<ned);i++)
        rvec_dec(x[i],transvec[0]);
}

void transfit(int ned, rvec *x, rvec *transvec) {
    int i;
    /* add the translation vector (this vector has only equal elements)*/
    for(i=0;(i<ned);i++)
        rvec_inc(x[i],transvec[0]);
}

void rmfit(int ned,rvec *x,rvec *transvec,matrix rotmat) {
    rmrotfit(ned,x,rotmat);
    rmtransfit(ned,x,transvec);
}

void rotate_vec(int nr,rvec *x,matrix rotmat)
{
    int i,j,k;
    rvec xdum;
    
    /* apply the rotation matrix */
    for(i=0;(i<nr);i++) {
        for(j=0;(j<3);j++)
            xdum[j]=x[i][j];
        for(j=0;(j<3);j++) {
            x[i][j]=0;
            for(k=0;(k<3);k++)
                x[i][j]+=rotmat[k][j]*xdum[k];
        }
    }
}

void ed_cons(rvec *x,t_edpar *edi,int step)
{
    int i;
    
    /* subtract the average positions */
    for (i=0;(i<edi->sav.nr);i++) 
        rvec_dec(x[edi->sav.anrs[i]],edi->sav.x[i]);
    
    /* apply the constraints */
    if (step >= 0) do_linfix(x,edi,step);
    do_linacc(x,edi);
    if (step >= 0) do_radfix(x,edi,step);
    do_radacc(x,edi);
    do_radcon(x,edi);
    
    /* add average positions */
    for (i=0;(i<edi->sav.nr);i++) 
        rvec_inc(x[edi->sav.anrs[i]],edi->sav.x[i]);
}

void do_linfix(rvec *x,t_edpar *edi,int step)
{
  int i,j,k;
  real proj,add;

  /* loop over linfix vectors */
  for (i=0;(i < edi->vecs.linfix.neig);i++) {
    /* calculate the projection */
    proj=projectx(edi,x,edi->vecs.linfix.vec[i]);
    /*fprintf(stderr,"Proj[%d]=%f\n",edi->vecs.linfix.ieig[i],proj);*/
    /* calculate the correction */
    add=edi->vecs.linfix.refproj[i]+step*edi->vecs.linfix.stpsz[i]-proj;
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=add*edi->vecs.linfix.vec[i][j][k]/edi->sav.sqrtm[i];
    }
  }
}

void do_linacc(rvec *x,t_edpar *edi)
{
  int i,j,k;
  real proj,add;

  /* loop over linacc vectors */
  for (i=0;(i<edi->vecs.linacc.neig);i++) {
    /* calculate the projection */
    proj=projectx(edi,x,edi->vecs.linacc.vec[i]);
    /* calculate the correction */
    add=0.0;
    if (edi->vecs.linacc.stpsz[i] > 0.0) {
      if ((proj-edi->vecs.linacc.refproj[i]) < 0.0)
	add=edi->vecs.linacc.refproj[i]-proj;
    }
    if (edi->vecs.linacc.stpsz[i] < 0.0) {
      if ((proj-edi->vecs.linacc.refproj[i]) > 0.0)
	add=edi->vecs.linacc.refproj[i]-proj;
    }
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=add*edi->vecs.linacc.vec[i][j][k]/edi->sav.sqrtm[i];
    }
    /* new positions will act as reference */
    edi->vecs.linacc.refproj[i]=proj+add;
  }
}

void do_radfix(rvec *x,t_edpar *edi,int step)
{
  int i,j,k;
  real *proj,rad=0.0,ratio;

  if (edi->vecs.radfix.neig == 0) return;
  snew(proj,edi->vecs.radfix.neig);
  /* loop over radfix vectors */
  for (i=0;(i<edi->vecs.radfix.neig);i++) {
    /* calculate the projections, radius */
    proj[i]=projectx(edi,x,edi->vecs.radfix.vec[i]);
    rad+=pow((proj[i]-edi->vecs.radfix.refproj[i]),2);
  }
  rad=sqrt(rad);
  ratio=(edi->vecs.radfix.stpsz[0]+edi->vecs.radfix.radius)/rad-1.0;
  edi->vecs.radfix.radius+=edi->vecs.radfix.stpsz[0];

  /* loop over radfix vectors */
  for (i=0;(i<edi->vecs.radfix.neig);i++) {
    proj[i]-=edi->vecs.radfix.refproj[i];
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=proj[i]*ratio*edi->vecs.radfix.vec[i][j][k]/edi->sav.sqrtm[i];
    }
  }  
  sfree(proj);
}

void do_radacc(rvec *x,t_edpar *edi)
{
  int i,j,k;
  real *proj,rad=0.0,ratio=0.0;

  if (edi->vecs.radacc.neig == 0) return;
  snew(proj,edi->vecs.radacc.neig);
  /* loop over radacc vectors */
  for (i=0;(i<edi->vecs.radacc.neig);i++) {
    /* calculate the projections, radius */
    proj[i]=projectx(edi,x,edi->vecs.radacc.vec[i]);
    /*    fprintf(stderr,"Proj[%d]=%f\n",edi->vecs.radacc.ieig[i],proj[i]); */
    rad+=pow((proj[i]-edi->vecs.radacc.refproj[i]),2);
  }
  rad=sqrt(rad);
  /*  fprintf(stderr,"Radius=%f old = %f\n",rad,edi->vecs.radacc.radius); */

  /* only correct when radius decreased */
  if (rad < edi->vecs.radacc.radius) {
    ratio=edi->vecs.radacc.radius/rad-1.0;
    rad=edi->vecs.radacc.radius;
  }
  else
    edi->vecs.radacc.radius=rad;

  /* loop over radacc vectors */
  for (i=0;(i<edi->vecs.radacc.neig);i++) {
    proj[i]-=edi->vecs.radacc.refproj[i];
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=proj[i]*ratio*edi->vecs.radacc.vec[i][j][k]/edi->sav.sqrtm[i];
    }
  }  
  sfree(proj);
}

struct t_do_radcon {
  real *proj;
};

void free_t_do_radcon(struct t_do_radcon *p) {
  sfree(p->proj);
  sfree(p);
}

void do_radcon(rvec *x,t_edpar *edi)
{
  int i,j,k;
  real rad=0.0,ratio=0.0;
  struct t_do_radcon *loc;
  bool bFirst;

  if(edi->local->do_radcon != NULL)
  {
      bFirst = FALSE;
      loc    = edi->local->do_radcon;
  }
  else
  {
      bFirst = TRUE;
      snew(edi->local->do_radcon,1);
  }
  loc = edi->local->do_radcon;

  if (edi->vecs.radcon.neig == 0) return;
  if (bFirst)
    snew(loc->proj,edi->vecs.radcon.neig);
  /* loop over radcon vectors */
  for (i=0;(i<edi->vecs.radcon.neig);i++) {
    /* calculate the projections, radius */
    loc->proj[i]=projectx(edi,x,edi->vecs.radcon.vec[i]);
    rad+=pow((loc->proj[i]-edi->vecs.radcon.refproj[i]),2);
  }
  rad=sqrt(rad);
  /*  fprintf(stderr,"Radius=%f old = %f\n",rad,edi->vecs.radcon.radius); */

  /* only correct when radius increased */
  if (rad > edi->vecs.radcon.radius) {
    ratio=edi->vecs.radcon.radius/rad-1.0;
  }
  else
    edi->vecs.radcon.radius=rad;
  /* loop over radcon vectors */
  for (i=0;(i<edi->vecs.radcon.neig);i++) {
    /*    proj[i]-=edi->vecs.radcon.refproj[i]; */
    /* apply the correction */
    for (j=0;(j<edi->sav.nr);j++) {
      for (k=0;(k<DIM);k++)
	x[edi->sav.anrs[j]][k]+=(loc->proj[i]-edi->vecs.radcon.refproj[i])*ratio*edi->vecs.radcon.vec[i][j][k]/edi->sav.sqrtm[i];
    }
  }  
  if (rad!=edi->vecs.radcon.radius) {
    rad=0.0;
    for (i=0;(i<edi->vecs.radcon.neig);i++) {
      /* calculate the projections, radius */
      loc->proj[i]=projectx(edi,x,edi->vecs.radcon.vec[i]);
      rad+=pow((loc->proj[i]-edi->vecs.radcon.refproj[i]),2);
    }
    rad=sqrt(rad);
    /*    fprintf(stderr,"after correction: Radius=%f old = %f\n",rad,edi->vecs.radcon.radius); */
  }
}

void write_edo(t_edpar *edi,int step,real rmsd)
{
  fprintf(edi->edo,"ED: %d\n",step);
  fprintf(edi->edo,"ED: %f %f %f %f \n",rmsd,edi->flood.Efl,edi->flood.deltaF,edi->flood.Vfl);
  write_proj(edi->edo,edi,"x");
  write_proj(edi->edo,edi,"v");
  write_proj(edi->edo,edi,"f");
}

void write_proj(FILE *out,t_edpar *edi,char *mode)
{
  do_write_proj(out,&edi->vecs.mon,mode);
  do_write_proj(out,&edi->vecs.linfix,mode);
  do_write_proj(out,&edi->vecs.linacc,mode);
  do_write_proj(out,&edi->vecs.radfix,mode);
  do_write_proj(out,&edi->vecs.radacc,mode);
  do_write_proj(out,&edi->vecs.radcon,mode);
}

void do_write_proj(FILE *out,t_eigvec *vec,char *mode)
{
  int i;
  if (vec->neig)
    fprintf(out,"ED%s: ",mode);
  for (i=0;(i<vec->neig);i++) {
    if (strcmp(mode,"x") == 0)
      fprintf(out,"%e ",vec->xproj[i]);
    else if (strcmp(mode,"v") == 0)
      fprintf(out,"%e ",vec->vproj[i]);
    else if (strcmp(mode,"f") == 0)
      fprintf(out,"%e ",vec->fproj[i]);
  }
  if (vec->neig > 0) fprintf(out,"\n");
}

void write_edidx(FILE *out,t_edpar *edi)
{
  int i;

  fprintf(out,"monitor eigenvectors");
  for (i=0;(i<edi->vecs.mon.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.mon.ieig[i],edi->vecs.mon.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"linfix  eigenvectors");
  for (i=0;(i<edi->vecs.linfix.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.linfix.ieig[i],edi->vecs.linfix.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"linacc  eigenvectors");
  for (i=0;(i<edi->vecs.linacc.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.linacc.ieig[i],edi->vecs.linacc.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"radfix  eigenvectors");
  for (i=0;(i<edi->vecs.radfix.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.radfix.ieig[i],edi->vecs.radfix.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"radacc  eigenvectors");
  for (i=0;(i<edi->vecs.radacc.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.radacc.ieig[i],edi->vecs.radacc.xproj[i]);
  fprintf(out,"\n");
  fprintf(out,"radcon  eigenvectors");
  for (i=0;(i<edi->vecs.radcon.neig);i++)
    fprintf(out," %d: %f;",edi->vecs.radcon.ieig[i],edi->vecs.radcon.xproj[i]);
  fprintf(out,"\n");
}

extern int ed_constraints(t_edsamyn *edyn){ 
    /* returns if any constraints are switched on */
  t_edpar *edi;
  if (edyn->edpar) {
    edi=edyn->edpar;
    return edi->vecs.linfix.neig || edi->vecs.linacc.neig || edi->vecs.radfix.neig ||  edi->vecs.radacc.neig ||  edi->vecs.radcon.neig;
  } 
  return 0;
}



/****************************************************************************
************************* free_xxx routines *********************************
--------- free_t_xxx does not free the given pointer, just the contents --------
---------- free_xxx does free the given pointer as well -----------------------
******************************************************************************/

void free_t_edvec(t_eigvec *vecs) {
  int i;
  if (vecs->neig>0) {
    sfree(vecs->ieig);
    sfree(vecs->stpsz);
    sfree(vecs->xproj); sfree(vecs->vproj); sfree(vecs->fproj);
    sfree(vecs->refproj);
    for(i=0; i<vecs->neig; i++)
      sfree(vecs->vec[i]);
    sfree(vecs->vec);
  }
}

void free_t_edx(t_edx *edx) {
  if (edx->nr>0) {
    sfree(edx->anrs);
    sfree(edx->x);
    if (edx->sqrtm)
      sfree(edx->sqrtm);
  };
}

void free_t_edlocals(t_edlocals *p) {
  sfree(p->x);  /* is called only when flooding is used */
  sfree(p->transvec);
  sfree(p->forces_cartesian);
}


void free_t_edflood(t_edflood *p) {
  if (p->vecs.neig >0) {
    free_t_edlocals(&(p->loc));
    free_t_edvec(&(p->vecs));
  }
}


void free_single_edpar(t_edpar *edi) {
  if (edi->x_unc) 
    sfree(edi->x_unc);
  free_t_edx(&(edi->sref));
  free_t_edx(&(edi->sav));
  free_t_edx(&(edi->sori));
  free_t_edx(&(edi->star));
  free_t_edvec(&(edi->vecs.mon));
  free_t_edvec(&(edi->vecs.linfix));
  free_t_edvec(&(edi->vecs.linacc));
  free_t_edvec(&(edi->vecs.radfix));
  free_t_edvec(&(edi->vecs.radacc));
  free_t_edvec(&(edi->vecs.radcon));
  free_t_edflood(&(edi->flood));
  sfree(edi->masnrs);
  sfree(edi->mass);
  sfree(edi->fitnrs);
  fclose(edi->edo);
  free_local(edi->local);
  sfree(edi);
}

void finish_edsam(FILE *log,t_topology *top,t_inputrec *ir,
		t_mdatoms *md,int start,int homenr,t_commrec *cr,
	   t_edsamyn *edyn) {
  t_edpar *next,*cur;  
  if(!MASTER(cr))
    return;
  if (!edyn->bEdsam || !edyn->edpar)
    return;
  fprintf(stderr,"finish edsam\n");
  cur=edyn->edpar;
  while (cur) {
    next=cur->next_edi;
    free_single_edpar(cur);
    cur=next;
  };
  fprintf(stderr,"finished finish edsam\n");
}

