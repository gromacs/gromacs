#ifndef _correct_h
#define _correct_h

#include <stdio.h>
#include "typedefs.h"

enum { edcNONBOND, edcBOND, edcDISRE, edcNR }; 

extern char *edc_names[edcNR+1];

typedef struct {
  int  ai,aj,cons_type;
  real lb,ub,wi;
} t_dist;

typedef struct {
  int ai,aj,ak,al,am;
} t_quadruple;

/* Garbage bin with flags and data for the correct routine */
typedef struct {
  int         maxnit,nbcheck,nstprint,nstranlist,ngrow;
  bool        bExplicit,bChiral,bPep,bDump,bLowerOnly;
  real        lodev;
  int         maxdist,ndist;
  int         *ip,*tag;
  t_dist      *d;
  int         npep,nimp;
  t_quadruple *pepbond,*imp;
  real        *omega,*idih;
  bool        *bViol;
} t_correct;

extern t_correct *init_corr(int maxnit,int nstprint,int nbcheck,int nstranlist,
			    int ngrow,bool bExplicit,bool bChiral,bool bPep,
			    bool bDump,real lowdev,bool bLowerOnly);
/* Initiate the data structure and set some of the parameters */

extern void make_tags(t_correct *c,int natom);
/* Make tags (indices) for optimizationof shaking */

extern bool shake_coords(FILE *log,bool bVerbose,int nstruct,
			 int natom,rvec xref[],rvec x[],int *seed,
			 matrix box,t_correct *c,int *niter);
/* Do the actual shaking. Return TRUE when converged */

extern int quick_check(FILE *log,int natom,rvec x[],matrix box,t_correct *c);
/* Check a structure once for the number of violations (return value) */

extern real *read_weights(char *fn,int natom);
/* Read the weights from the occupancy field in the pdb file */  

extern void define_peptide_bonds(FILE *log,t_atoms *atoms,t_correct *c);
/* Fill the peptide bond structure */

extern void define_impropers(FILE *log,t_atoms *atoms,t_correct *c);
/* Fill the improper structure */

extern void define_dist(FILE *log,t_topology *top,t_correct *c,real weight[]);
/* Fill the normal distance restraints */

extern void read_dist(FILE *log,char *fn,int natom,t_correct *c,real weight[]);
/* Read distances from a dist.dat file produced by cdist */

extern void pr_corr(FILE *log,t_correct *c);
/* Print parameters from the corr structure */

extern void pr_distances(FILE *fp,t_correct *c);
/* Print the distances in a cdist compatible file format */

extern void center_in_box(int natom,rvec xin[],matrix box,rvec xout[]);
/* Center the coordinates in the box. xin and xout may be the same array */

extern void measure_dist(FILE *log,int natom,rvec x[],t_correct *c,
			 real weight[],real cutoff);
/* Measure the distances from the structure */

extern void check_dist(FILE *log,t_correct *c);
/* Check internal consistency of the distances */

extern void check_final(FILE *log,t_correct *c,rvec x[]);
/* Check structure for distances */

#endif
