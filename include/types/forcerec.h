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

enum { eNL_VDWQQ, eNL_VDW, eNL_QQ, 
       eNL_VDWQQ_FREE, eNL_VDW_FREE, eNL_QQ_FREE, 
       eNL_VDWQQ_SOLMNO, eNL_VDW_SOLMNO, eNL_QQ_SOLMNO, 
       eNL_VDWQQ_WATER, eNL_QQ_WATER, 
       eNL_VDWQQ_WATERWATER, eNL_QQ_WATERWATER, 
       eNL_NR };
 
typedef struct {
  /* Cut-Off stuff */
  int  ePBC;
  real rlist,rlistlong;
  
  /* Dielectric constant resp. multiplication factor for charges */
  real zsquare,temp;
  real epsilon_r,epsfac;  
  
  /* Constants for reaction fields */
  bool bRF;
  real kappa,k_rf,c_rf;
  
  /* Constant for long range dispersion correction (average dispersion) */
  real avcsix;
      
  /* Fudge factors */
  real fudgeQQ;

  /* Table stuff */
  bool bcoultab;
  bool bvdwtab;
  real rtab;
  int  ntab;
  real tabscale;
  /* We duplicate tables for cache optimization purposes */
  real *coultab;      /* Coul only */
  real *vdwtab;       /* Vdw only   */
  real *coulvdwtab;   /* Both      */
  real *coulvdw14tab; /* 1,4 table with both */

  /* PPPM & Shifting stuff */
  real rcoulomb_switch,rcoulomb;
  real *phi;

  /* VdW stuff */
  real rvdw_switch,rvdw;
  real bham_b_max;
  real tabscale_exp;

  /* Free energy ? */
  int  efep;
  real sc_alpha;
  real sc_sigma6;

  /* NS Stuff */
  int  eeltype;
  int  vdwtype;
  int  cg0,hcg;
  int  ndelta;
  bool bSolvOpt;
  int  nMNOMol;
  real nMNOav[3];
  int  nWatMol;
  int  Dimension;
  bool bGrid,bDomDecomp;
  int  *solvent_type;
  int  *mno_index;
  rvec *cg_cm;
  rvec *shift_vec;
  
  /* The actual neighbor lists, short and long range, see enum above
   * for definition of neighborlist indices.
   */
  t_nblist nlist_sr[eNL_NR];
  t_nblist nlist_lr[eNL_NR];
  
  /* This mask array of length nn determines whether or not this bit of the
   * neighbourlists should be computed. Usually all these are true of course,
   * but not when shells are used. During minimisation all the forces that 
   * include shells are done, then after minimsation is converged the remaining
   * forces are computed.
   */
  /* bool *bMask; */
    
  /* Twin Range stuff. */
  bool bTwinRange;
  int  nlr;
  rvec *flr;
  rvec *fshift_lr;
  
  /* PME/Ewald stuff */
  bool bEwald;
  real ewaldcoeff;

  /* Virial Stuff */
  rvec *fshift;
  
  /* Free energy stuff */
  int     nmol;
  atom_id *mol_nr;
  real    *mol_epot;
  int     nstcalc;
  
  /* Non bonded Parameter lists */
  int  ntype; /* Number of atom types */
  bool bBHAM;
  real *nbfp;

  /* Energy group exclusions */
  bool *eg_excl;
} t_forcerec;

#define C6(nbfp,ntp,ai,aj)     (nbfp)[2*((ntp)*(ai)+(aj))]
#define C12(nbfp,ntp,ai,aj)    (nbfp)[2*((ntp)*(ai)+(aj))+1]
#define BHAMC(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))]
#define BHAMA(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))+1]
#define BHAMB(nbfp,ntp,ai,aj)  (nbfp)[3*((ntp)*(ai)+(aj))+2]
