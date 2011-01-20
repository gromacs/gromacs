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
 * GRoups of Organic Molecules in ACtion for Science
 */


#ifndef _idef_h
#define _idef_h

#include "simple.h"

#ifdef __cplusplus
extern "C" {
#endif


/* check kernel/toppush.c when you change these numbers */
#define MAXATOMLIST	6
#define MAXFORCEPARAM	12
#define NR_RBDIHS	6
#define NR_FOURDIHS     4

typedef atom_id t_iatom;

/* this MUST correspond to the 
   t_interaction_function[F_NRE] in gmxlib/ifunc.c */
enum {
  F_BONDS,
  F_G96BONDS,
  F_MORSE,
  F_CUBICBONDS,
  F_CONNBONDS,
  F_HARMONIC,
  F_FENEBONDS,
  F_TABBONDS,
  F_TABBONDSNC,
  F_RESTRBONDS,
  F_ANGLES, 
  F_G96ANGLES,
  F_CROSS_BOND_BONDS,
  F_CROSS_BOND_ANGLES,
  F_UREY_BRADLEY,
  F_QUARTIC_ANGLES,
  F_TABANGLES,
  F_PDIHS,
  F_RBDIHS, 
  F_FOURDIHS,
  F_IDIHS, 
  F_PIDIHS, 
  F_TABDIHS,
  F_CMAP,
  F_GB12,
  F_GB13,
  F_GB14,
  F_GBPOL,
  F_NPSOLVATION,
  F_LJ14,
  F_COUL14,
  F_LJC14_Q,
  F_LJC_PAIRS_NB,
  F_LJ,
  F_BHAM,
  F_LJ_LR,
  F_BHAM_LR,
  F_DISPCORR,
  F_COUL_SR,
  F_COUL_LR,
  F_RF_EXCL,
  F_COUL_RECIP,
  F_DPD,
  F_POLARIZATION,
  F_WATER_POL,
  F_THOLE_POL,
  F_POSRES,
  F_DISRES,
  F_DISRESVIOL,
  F_ORIRES,
  F_ORIRESDEV,
  F_ANGRES,
  F_ANGRESZ,
  F_DIHRES,
  F_DIHRESVIOL,
  F_CONSTR,
  F_CONSTRNC,
  F_SETTLE,
  F_VSITE2,
  F_VSITE3,
  F_VSITE3FD,
  F_VSITE3FAD,
  F_VSITE3OUT,
  F_VSITE4FD,
  F_VSITE4FDN,
  F_VSITEN,
  F_COM_PULL,
  F_EQM,
  F_EPOT,
  F_EKIN,
  F_ETOT,
  F_ECONSERVED,
  F_TEMP,
  F_VTEMP,
  F_PDISPCORR,
  F_PRES,
  F_DVDL,
  F_DKDL,
  F_DHDL_CON,
  F_NRE		/* This number is for the total number of energies	*/
};
  
typedef union
{
  /* Some parameters have A and B values for free energy calculations.
   * The B values are not used for regular simulations of course.
   * Free Energy for nonbondeds can be computed by changing the atom type.
   * The harmonic type is used for all harmonic potentials:
   * bonds, angles and improper dihedrals
   */
  struct {real a,b,c;	                                   } bham;
  struct {real rA,krA,rB,krB;           	           } harmonic;
  struct {real lowA,up1A,up2A,kA,lowB,up1B,up2B,kB;        } restraint;
  /* No free energy supported for cubic bonds, FENE, WPOL or cross terms */ 
  struct {real b0,kb,kcub;                                 } cubic;
  struct {real bm,kb;                                      } fene;
  struct {real r1e,r2e,krr;                                } cross_bb;
  struct {real r1e,r2e,r3e,krt;                            } cross_ba;
  struct {real theta,ktheta,r13,kUB;                       } u_b;
  struct {real theta,c[5];                                 } qangle; 
  struct {real alpha;                                      } polarize;
  struct {real al_x,al_y,al_z,rOH,rHH,rOD;                 } wpol;
  struct {real a,alpha1,alpha2,rfac;                       } thole;
  struct {real c6,c12;				           } lj;
  struct {real c6A,c12A,c6B,c12B;		           } lj14;
  struct {real fqq,qi,qj,c6,c12;	                   } ljc14;
  struct {real qi,qj,c6,c12;		                   } ljcnb;
  /* Proper dihedrals can not have different multiplicity when
   * doing free energy calculations, because the potential would not
   * be periodic anymore.
   */ 
  struct {real phiA,cpA;int mult;real phiB,cpB;            } pdihs;
  struct {real dA,dB;		        	           } constr;
  /* Settle can not be used for Free energy calculations of water bond geometry.
   * Use shake (or lincs) instead if you have to change the water bonds.
   */
  struct {real doh,dhh;                                   } settle;
  /* No free energy supported for morse bonds */ 
  struct {real b0,cb,beta;                        	  } morse;
  struct {real pos0A[DIM],fcA[DIM],pos0B[DIM],fcB[DIM];   } posres;
  struct {real rbcA[NR_RBDIHS], rbcB[NR_RBDIHS];          } rbdihs;
  struct {real a,b,c,d,e,f;                               } vsite;   
  struct {int  n; real a;                                 } vsiten;   
  /* NOTE: npair is only set after reading the tpx file */
  struct {real low,up1,up2,kfac;int type,label,npair;     } disres; 
  struct {real phi,dphi,kfac;int label,power;             } dihres;  
  struct {int  ex,power,label; real c,obs,kfac;           } orires;
  struct {int  table;real kA;real kB;                     } tab;
  struct {real sar,st,pi,gbr,bmlt;                        } gb;
  struct {int cmapA,cmapB;                                } cmap;
  struct {real buf[MAXFORCEPARAM];	  	          } generic; /* Conversion */
} t_iparams;

typedef int t_functype;

/*
 * The nonperturbed/perturbed interactions are now separated (sorted) in the
 * ilist, such that the first 0..(nr_nonperturbed-1) ones are exactly that, and 
 * the remaining ones from nr_nonperturbed..(nr-1) are perturbed bonded 
 * interactions.
 */
typedef struct
{
  int nr;
  int nr_nonperturbed;
  t_iatom *iatoms;
  int nalloc;
} t_ilist;

/*
 * The struct t_ilist defines a list of atoms with their interactions. 
 * General field description:
 *   int nr
 *	the size (nr elements) of the interactions array (iatoms[]).
 *   t_iatom *iatoms
 * 	specifies which atoms are involved in an interaction of a certain 
 *       type. The layout of this array is as follows:
 *
 *	  +-----+---+---+---+-----+---+---+-----+---+---+---+-----+---+---+...
 *	  |type1|at1|at2|at3|type2|at1|at2|type1|at1|at2|at3|type3|at1|at2|
 *	  +-----+---+---+---+-----+---+---+-----+---+---+---+-----+---+---+...
 *
 * 	So for interaction type type1 3 atoms are needed, and for type2 and 
 *      type3 only 2. The type identifier is used to select the function to 
 *	calculate the interaction and its actual parameters. This type 
 *	identifier is an index in a params[] and functype[] array.
 */

typedef struct
{
	real *cmap; /* Has length 4*grid_spacing*grid_spacing, */
	/* there are 4 entries for each cmap type (V,dVdx,dVdy,d2dVdxdy) */
} cmapdata_t;

typedef struct
{
	int ngrid;            /* Number of allocated cmap (cmapdata_t ) grids */
	int grid_spacing;     /* Grid spacing */
	cmapdata_t *cmapdata; /* Pointer to grid with actual, pre-interpolated data */
} gmx_cmap_t;


typedef struct
{
  int        ntypes;
  int        atnr;
  t_functype *functype;
  t_iparams  *iparams;
  double     reppow;     /* The repulsion power for VdW: C12*r^-reppow   */
  real       fudgeQQ;    /* The scaling factor for Coulomb 1-4: f*q1*q2  */
  gmx_cmap_t cmap_grid;  /* The dihedral correction maps                 */
} gmx_ffparams_t;

enum {
  ilsortUNKNOWN, ilsortNO_FE, ilsortFE_UNSORTED, ilsortFE_SORTED
};

typedef struct
{
  int ntypes;
  int atnr;
  t_functype *functype;
  t_iparams  *iparams;
  real fudgeQQ;
  gmx_cmap_t cmap_grid;
  t_iparams  *iparams_posres;
  int iparams_posres_nalloc;

  t_ilist il[F_NRE];
  int ilsort;
} t_idef;

/*
 * The struct t_idef defines all the interactions for the complete
 * simulation. The structure is setup in such a way that the multinode
 * version of the program  can use it as easy as the single node version.
 * General field description:
 *   int ntypes
 *	defines the number of elements in functype[] and param[].
 *   int nodeid
 *      the node id (if parallel machines)
 *   int atnr
 *      the number of atomtypes
 *   t_functype *functype
 *	array of length ntypes, defines for every force type what type of 
 *      function to use. Every "bond" with the same function but different 
 *	force parameters is a different force type. The type identifier in the 
 *	forceatoms[] array is an index in this array.
 *   t_iparams *iparams
 *	array of length ntypes, defines the parameters for every interaction
 *      type. The type identifier in the actual interaction list
 *      (ilist[ftype].iatoms[]) is an index in this array.
 *   gmx_cmap_t cmap_grid
 *      the grid for the dihedral pair correction maps.
 *   t_iparams *iparams_posres
 *	defines the parameters for position restraints only.
 *      Position restraints are the only interactions that have different
 *      parameters (reference positions) for different molecules
 *      of the same type. ilist[F_POSRES].iatoms[] is an index in this array.
 *   t_ilist il[F_NRE]
 *      The list of interactions for each type. Note that some,
 *      such as LJ and COUL will have 0 entries.
 */

typedef struct {
  int  n;         /* n+1 is the number of points */
  real scale;     /* distance between two points */
  real *tab;      /* the actual tables, per point there are  4 numbers */
} bondedtable_t;

#ifdef __cplusplus
}
#endif


#endif
