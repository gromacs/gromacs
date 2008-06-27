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
 * GROningen Mixture of Alchemy and Childrens' Stories
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

/* This file is completely threadsafe - keep it that way! */
#include <gmx_thread.h>


#include <ctype.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "macros.h"
#include "names.h"
#include "symtab.h"
#include "futil.h"
#include "filenm.h"
#include "gmxfio.h"
#include "topsort.h"
#include "tpxio.h"
#include "txtdump.h"
#include "confio.h"
#include "atomprop.h"
#include "copyrite.h"
#include "vec.h"
#ifdef HAVE_LIBXML2
#include "xmlio.h"
#endif

/* This number should be increased whenever the file format changes! */
static const int tpx_version = 56;

/* This number should only be increased when you edit the TOPOLOGY section
 * of the tpx format. This way we can maintain forward compatibility too
 * for all analysis tools and/or external programs that only need to
 * know the atom/residue names, charges, and bond connectivity.
 *  
 * It first appeared in tpx version 26, when I also moved the inputrecord
 * to the end of the tpx file, so we can just skip it if we only
 * want the topology.
 */
static const int tpx_generation = 16;

/* This number should be the most recent backwards incompatible version 
 * I.e., if this number is 9, we cannot read tpx version 9 with this code.
 */
static const int tpx_incompatible_version = 9;



/* Struct used to maintain tpx compatibility when function types are added */
typedef struct {
  int fvnr; /* file version number in which the function type first appeared */
  int ftype; /* function type */
} t_ftupd;

/* 
 *The entries should be ordered in:
 * 1. ascending file version number
 * 2. ascending function type number
 */
/*static const t_ftupd ftupd[] = {
  { 20, F_CUBICBONDS        },
  { 20, F_CONNBONDS         },
  { 20, F_HARMONIC          },
  { 20, F_EQM,              },
  { 22, F_DISRESVIOL        },
  { 22, F_ORIRES            },
  { 22, F_ORIRESDEV         },
  { 26, F_FOURDIHS          },
  { 26, F_PIDIHS            },
  { 26, F_DIHRES            },
  { 26, F_DIHRESVIOL        },
  { 30, F_CROSS_BOND_BONDS  },
  { 30, F_CROSS_BOND_ANGLES },
  { 30, F_UREY_BRADLEY      },
  { 30, F_POLARIZATION      }
  };*/
  
/* 
 *The entries should be ordered in:
 * 1. ascending function type number
 * 2. ascending file version number
 */
static const t_ftupd ftupd[] = {
  { 20, F_CUBICBONDS        },
  { 20, F_CONNBONDS         },
  { 20, F_HARMONIC          },
  { 34, F_FENEBONDS         },
  { 43, F_TABBONDS          },
  { 43, F_TABBONDSNC        },
  { 30, F_CROSS_BOND_BONDS  },
  { 30, F_CROSS_BOND_ANGLES },
  { 30, F_UREY_BRADLEY      },
  { 34, F_QUARTIC_ANGLES    },
  { 43, F_TABANGLES         },
  { 26, F_FOURDIHS          },
  { 26, F_PIDIHS            },
  { 43, F_TABDIHS           },
  { 41, F_LJC14_Q           },
  { 41, F_LJC_PAIRS_NB      },
  { 32, F_BHAM_LR           },
  { 32, F_RF_EXCL           },
  { 32, F_COUL_RECIP        },
  { 46, F_DPD               },
  { 30, F_POLARIZATION      },
  { 36, F_THOLE_POL         },
  { 22, F_DISRESVIOL        },
  { 22, F_ORIRES            },
  { 22, F_ORIRESDEV         },
  { 26, F_DIHRES            },
  { 26, F_DIHRESVIOL        },
  { 49, F_VSITE4FDN         },
  { 50, F_VSITEN            },
  { 46, F_COM_PULL          },
  { 20, F_EQM               },
  { 46, F_ECONSERVED        },
  { 54, F_DGDL_CON          }
};
#define NFTUPD asize(ftupd)

/* Needed for backward compatibility */
#define MAXNODES 256

void _do_section(int fp,int key,bool bRead,char *src,int line)
{
  char buf[STRLEN];
  bool bDbg;

  if (gmx_fio_getftp(fp) == efTPA) {
    if (!bRead) {
      do_string(itemstr[key]);
      bDbg       = gmx_fio_getdebug(fp);
      gmx_fio_setdebug(fp,FALSE);
      do_string(comment_str[key]);
      gmx_fio_setdebug(fp,bDbg);
    }
    else {
      if (gmx_fio_getdebug(fp))
	fprintf(stderr,"Looking for section %s (%s, %d)",
		itemstr[key],src,line);
      
      do {
	do_string(buf);
      } while ((strcasecmp(buf,itemstr[key]) != 0));
      
      if (strcasecmp(buf,itemstr[key]) != 0) 
	gmx_fatal(FARGS,"\nCould not find section heading %s",itemstr[key]);
      else if (gmx_fio_getdebug(fp))
	fprintf(stderr," and found it\n");
    }
  }
}

#define do_section(key,bRead) _do_section(fp,key,bRead,__FILE__,__LINE__)

/**************************************************************
 *
 * Now the higer level routines that do io of the structures and arrays
 *
 **************************************************************/
static void do_pullgrp(t_pullgrp *pgrp,bool bRead, int file_version)
{
  bool bDum=TRUE;
  int  i;

  do_int(pgrp->nat);
  if (bRead)
    snew(pgrp->ind,pgrp->nat);
  ndo_int(pgrp->ind,pgrp->nat,bDum);
  do_int(pgrp->nweight);
  if (bRead)
    snew(pgrp->weight,pgrp->nweight);
  ndo_real(pgrp->weight,pgrp->nweight,bDum);
  do_int(pgrp->pbcatom);
  do_rvec(pgrp->vec);
  do_rvec(pgrp->init);
  do_real(pgrp->rate);
  do_real(pgrp->k);
  if (file_version >= 56) {
    do_real(pgrp->kB);
  } else {
    pgrp->kB = pgrp->k;
  }
}

static void do_pull(t_pull *pull,bool bRead, int file_version)
{
  int g;

  do_int(pull->ngrp);
  do_int(pull->eGeom);
  do_ivec(pull->dim);
  do_real(pull->cyl_r1);
  do_real(pull->cyl_r0);
  do_real(pull->constr_tol);
  do_int(pull->nstxout);
  do_int(pull->nstfout);
  if (bRead)
    snew(pull->grp,pull->ngrp+1);
  for(g=0; g<pull->ngrp+1; g++)
    do_pullgrp(&pull->grp[g],bRead,file_version);
}

static void do_inputrec(t_inputrec *ir,bool bRead, int file_version,
			real *fudgeQQ)
{
  int  i,j,k,*tmp,idum=0; 
  bool bDum=TRUE;
  real rdum,bd_temp;
  rvec vdum;
  bool bSimAnn;
  real zerotemptime,finish_t,init_temp,finish_temp;
  
  if (file_version != tpx_version) {
    /* Give a warning about features that are not accessible */
    fprintf(stderr,"Note: tpx file_version %d, software version %d\n",
	    file_version,tpx_version);
  }

  if (bRead)
    init_inputrec(ir);

  if (file_version >= 1) {  
    /* Basic inputrec stuff */  
    do_int(ir->eI); 
    do_int(ir->nsteps); 
    if(file_version > 25)
      do_int(ir->init_step);
    else
      ir->init_step=0;

    if (file_version < 53) {
      /* The pbc info has been moved out of do_inputrec,
       * since we always want it, also without reading the inputrec.
       */
      do_int(ir->ePBC);
      if ((file_version <= 15) && (ir->ePBC == 2))
	ir->ePBC = epbcNONE;
      if (file_version >= 45) {
	do_int(ir->bPeriodicMols);
      } else {
	if (ir->ePBC == 2) {
	  ir->ePBC = epbcXYZ;
	  ir->bPeriodicMols = TRUE;
	} else {
	ir->bPeriodicMols = FALSE;
	}
      }
    }
    do_int(ir->ns_type);
    do_int(ir->nstlist);
    do_int(ir->ndelta);
    if (file_version < 41) {
      do_int(idum);
      do_int(idum);
    }
    if (file_version >= 45)
      do_real(ir->rtpi);
    else
      ir->rtpi = 0.05;
    do_int(ir->nstcomm); 
    if (file_version > 34)
      do_int(ir->comm_mode);
    else if (ir->nstcomm < 0) 
      ir->comm_mode = ecmANGULAR;
    else
      ir->comm_mode = ecmLINEAR;
    ir->nstcomm = abs(ir->nstcomm);
    
    if(file_version > 25)
      do_int(ir->nstcheckpoint);
    else
      ir->nstcheckpoint=0;
    
    do_int(ir->nstcgsteep); 

    if(file_version>=30)
      do_int(ir->nbfgscorr); 
    else if (bRead)
      ir->nbfgscorr = 10;

    do_int(ir->nstlog); 
    do_int(ir->nstxout); 
    do_int(ir->nstvout); 
    do_int(ir->nstfout); 
    do_int(ir->nstenergy); 
    do_int(ir->nstxtcout); 
    do_real(ir->init_t); 
    do_real(ir->delta_t); 
    do_real(ir->xtcprec); 
    if (file_version < 19) {
      do_int(idum); 
      do_int(idum);
    }
    if(file_version < 18)
      do_int(idum); 
    do_real(ir->rlist); 
    do_int(ir->coulombtype); 
    if (file_version < 32 && ir->coulombtype == eelRF)
      ir->coulombtype = eelRF_NEC;      
    do_real(ir->rcoulomb_switch); 
    do_real(ir->rcoulomb); 
    do_int(ir->vdwtype);
    do_real(ir->rvdw_switch); 
    do_real(ir->rvdw); 
    do_int(ir->eDispCorr); 
    do_real(ir->epsilon_r);
    if (file_version >= 37) {
      do_real(ir->epsilon_rf);
    } else {
      if (EEL_RF(ir->coulombtype)) {
	ir->epsilon_rf = ir->epsilon_r;
	ir->epsilon_r  = 1.0;
      } else {
	ir->epsilon_rf = 1.0;
      }
    }
    if (file_version >= 29)
      do_real(ir->tabext);
    else
      ir->tabext=1.0;
 
    if(file_version > 25) {
      do_int(ir->gb_algorithm);
      do_int(ir->nstgbradii);
      do_real(ir->rgbradii);
      do_real(ir->gb_saltconc);
      do_int(ir->implicit_solvent);
    } else {
      ir->gb_algorithm=egbSTILL;
      ir->nstgbradii=1;
      ir->rgbradii=1.0;
      ir->gb_saltconc=0;
      ir->implicit_solvent=eisNO;
    }
	if(file_version>=55)
	{
		do_real(ir->gb_epsilon_solvent);
		do_real(ir->gb_obc_alpha);
		do_real(ir->gb_obc_beta);
		do_real(ir->gb_obc_gamma);
		do_real(ir->sa_surface_tension);
	}
	else
	{
		/* Better use sensible values than insane (0.0) ones... */
		ir->gb_epsilon_solvent = 80;
		ir->gb_obc_alpha       = 1.0;
		ir->gb_obc_beta        = 0.8;
		ir->gb_obc_gamma       = 4.85;
		ir->sa_surface_tension = 2.092;
	}

	  
    do_int(ir->nkx); 
    do_int(ir->nky); 
    do_int(ir->nkz);
    do_int(ir->pme_order);
    do_real(ir->ewald_rtol);

    if (file_version >=24) 
      do_int(ir->ewald_geometry);
    else
      ir->ewald_geometry=eewg3D;

    if (file_version <=17) {
      ir->epsilon_surface=0;
      if (file_version==17)
	do_int(idum);
    } 
    else
      do_real(ir->epsilon_surface);
    
    do_int(ir->bOptFFT);

    do_int(ir->bContinuation); 
    do_int(ir->etc);
    /* before version 18, ir->etc was a bool (ir->btc),
     * but the values 0 and 1 still mean no and
     * berendsen temperature coupling, respectively.
     */
    if (file_version <= 15)
      do_int(idum);
    if (file_version <=17) {
      do_int(ir->epct); 
      if (file_version <= 15) {
	if (ir->epct == 5)
	  ir->epct = epctSURFACETENSION;
	do_int(idum);
      }
      ir->epct -= 1;
      /* we have removed the NO alternative at the beginning */
      if(ir->epct==-1) {
	ir->epc=epcNO;
	ir->epct=epctISOTROPIC;
      } 
      else
	ir->epc=epcBERENDSEN;
    } 
    else {
      do_int(ir->epc);
      do_int(ir->epct);
    }
    do_real(ir->tau_p); 
    if (file_version <= 15) {
      do_rvec(vdum);
      clear_mat(ir->ref_p);
      for(i=0; i<DIM; i++)
	ir->ref_p[i][i] = vdum[i];
    } else {
      do_rvec(ir->ref_p[XX]);
      do_rvec(ir->ref_p[YY]);
      do_rvec(ir->ref_p[ZZ]);
    }
    if (file_version <= 15) {
      do_rvec(vdum);
      clear_mat(ir->compress);
      for(i=0; i<DIM; i++)
	ir->compress[i][i] = vdum[i];
    } 
    else {
      do_rvec(ir->compress[XX]);
      do_rvec(ir->compress[YY]);
      do_rvec(ir->compress[ZZ]);
    }
    if (file_version >= 47) {
      do_int(ir->refcoord_scaling);
      do_rvec(ir->posres_com);
      do_rvec(ir->posres_comB);
    } else {
      ir->refcoord_scaling = erscNO;
      clear_rvec(ir->posres_com);
      clear_rvec(ir->posres_comB);
    }
    if(file_version > 25)
      do_int(ir->andersen_seed);
    else
      ir->andersen_seed=0;
    
    if(file_version < 26) {
      do_int(bSimAnn); 
      do_real(zerotemptime);
    }
    
    if (file_version < 37)
      do_real(rdum); 

    do_real(ir->shake_tol);
    if (file_version < 54)
      do_real(*fudgeQQ);
    do_int(ir->efep);
    if (file_version <= 14 && ir->efep > efepNO)
      ir->efep = efepYES;
    do_real(ir->init_lambda); 
    do_real(ir->delta_lambda);
    if (file_version >= 13)
      do_real(ir->sc_alpha);
    else
      ir->sc_alpha = 0;
    if (file_version >= 38)
      do_int(ir->sc_power);
    else
      ir->sc_power = 2;
    if (file_version >= 15)
      do_real(ir->sc_sigma);
    else
      ir->sc_sigma = 0.3;
    do_int(ir->eDisreWeighting); 
    if (file_version < 22) {
      if (ir->eDisreWeighting == 0)
	ir->eDisreWeighting = edrwEqual;
      else
	ir->eDisreWeighting = edrwConservative;
    }
    do_int(ir->bDisreMixed); 
    do_real(ir->dr_fc); 
    do_real(ir->dr_tau); 
    do_int(ir->nstdisreout);
    if (file_version >= 22) {
      do_real(ir->orires_fc);
      do_real(ir->orires_tau);
      do_int(ir->nstorireout);
    } else {
      ir->orires_fc = 0;
      ir->orires_tau = 0;
      ir->nstorireout = 0;
    }
    if(file_version >= 26) {
      do_real(ir->dihre_fc);
      if (file_version < 56) {
	do_real(rdum);
	do_int(idum);
      }
    } else {
      ir->dihre_fc=0;
    }

    do_real(ir->em_stepsize); 
    do_real(ir->em_tol); 
    if (file_version >= 22) 
      do_int(ir->bShakeSOR);
    else if (bRead)
      ir->bShakeSOR = TRUE;
    if (file_version >= 11)
      do_int(ir->niter);
    else if (bRead) {
      ir->niter = 25;
      fprintf(stderr,"Note: niter not in run input file, setting it to %d\n",
	      ir->niter);
    }
    if (file_version >= 21)
      do_real(ir->fc_stepsize);
    else
      ir->fc_stepsize = 0;
    do_int(ir->eConstrAlg);
    do_int(ir->nProjOrder);
    do_real(ir->LincsWarnAngle);
    if (file_version <= 14)
      do_int(idum);
    if (file_version >=26)
      do_int(ir->nLincsIter);
    else if (bRead) {
      ir->nLincsIter = 1;
      fprintf(stderr,"Note: nLincsIter not in run input file, setting it to %d\n",
	      ir->nLincsIter);
    }
    if (file_version < 33)
      do_real(bd_temp);
    do_real(ir->bd_fric);
    do_int(ir->ld_seed);
    if (file_version >= 33) {
      for(i=0; i<DIM; i++)
	do_rvec(ir->deform[i]);
    } else {
      for(i=0; i<DIM; i++)
	clear_rvec(ir->deform[i]);
    }
    if (file_version >= 14)
      do_real(ir->cos_accel);
    else if (bRead)
      ir->cos_accel = 0;
    do_int(ir->userint1); 
    do_int(ir->userint2); 
    do_int(ir->userint3); 
    do_int(ir->userint4); 
    do_real(ir->userreal1); 
    do_real(ir->userreal2); 
    do_real(ir->userreal3); 
    do_real(ir->userreal4); 
    
    /* pull stuff */
    if (file_version >= 48) {
      do_int(ir->ePull);
      if (ir->ePull != epullNO) {
	if (bRead)
	  snew(ir->pull,1);
	do_pull(ir->pull,bRead,file_version);
      }
    } else {
      ir->ePull = epullNO;
    }
    
    /* grpopts stuff */
    do_int(ir->opts.ngtc); 
    do_int(ir->opts.ngacc); 
    do_int(ir->opts.ngfrz); 
    do_int(ir->opts.ngener);

    if (bRead) {
      snew(ir->opts.nrdf,   ir->opts.ngtc); 
      snew(ir->opts.ref_t,  ir->opts.ngtc); 
      snew(ir->opts.annealing, ir->opts.ngtc); 
      snew(ir->opts.anneal_npoints, ir->opts.ngtc); 
      snew(ir->opts.anneal_time, ir->opts.ngtc); 
      snew(ir->opts.anneal_temp, ir->opts.ngtc); 
      snew(ir->opts.tau_t,  ir->opts.ngtc); 
      snew(ir->opts.nFreeze,ir->opts.ngfrz); 
      snew(ir->opts.acc,    ir->opts.ngacc); 
      snew(ir->opts.egp_flags,ir->opts.ngener*ir->opts.ngener);
    } 
    if (ir->opts.ngtc > 0) {
      if (bRead && file_version<13) {
	snew(tmp,ir->opts.ngtc);
	ndo_int (tmp, ir->opts.ngtc,bDum);
	for(i=0; i<ir->opts.ngtc; i++)
	  ir->opts.nrdf[i] = tmp[i];
	sfree(tmp);
      } else {
	ndo_real(ir->opts.nrdf, ir->opts.ngtc,bDum);
      }
      ndo_real(ir->opts.ref_t,ir->opts.ngtc,bDum); 
      ndo_real(ir->opts.tau_t,ir->opts.ngtc,bDum); 
      if (file_version<33 && ir->eI==eiBD) {
	for(i=0; i<ir->opts.ngtc; i++)
	  ir->opts.tau_t[i] = bd_temp;
      }
    }
    if (ir->opts.ngfrz > 0) 
      ndo_ivec(ir->opts.nFreeze,ir->opts.ngfrz,bDum);
    if (ir->opts.ngacc > 0) 
      ndo_rvec(ir->opts.acc,ir->opts.ngacc); 
    if (file_version >= 12)
      ndo_int(ir->opts.egp_flags,ir->opts.ngener*ir->opts.ngener,bDum);

    if(bRead && file_version < 26) {
      for(i=0;i<ir->opts.ngtc;i++) {
	if(bSimAnn) {
	  ir->opts.annealing[i] = eannSINGLE;
	  ir->opts.anneal_npoints[i] = 2;
	  snew(ir->opts.anneal_time[i],2);
	  snew(ir->opts.anneal_temp[i],2);
	  /* calculate the starting/ending temperatures from reft, zerotemptime, and nsteps */
	  finish_t = ir->init_t + ir->nsteps * ir->delta_t;
	  init_temp = ir->opts.ref_t[i]*(1-ir->init_t/zerotemptime);
	  finish_temp = ir->opts.ref_t[i]*(1-finish_t/zerotemptime);
	  ir->opts.anneal_time[i][0] = ir->init_t;
	  ir->opts.anneal_time[i][1] = finish_t;
	  ir->opts.anneal_temp[i][0] = init_temp;
	  ir->opts.anneal_temp[i][1] = finish_temp;
	} else {
	  ir->opts.annealing[i] = eannNO;
	  ir->opts.anneal_npoints[i] = 0;
	}
      }
    } else {
      /* file version 26 or later */
      /* First read the lists with annealing and npoints for each group */
      ndo_int(ir->opts.annealing,ir->opts.ngtc,bDum);
      ndo_int(ir->opts.anneal_npoints,ir->opts.ngtc,bDum);
      for(j=0;j<(ir->opts.ngtc);j++) {
	k=ir->opts.anneal_npoints[j];
	if(bRead) {
	  snew(ir->opts.anneal_time[j],k);
	  snew(ir->opts.anneal_temp[j],k);
	}
	ndo_real(ir->opts.anneal_time[j],k,bDum);
	ndo_real(ir->opts.anneal_temp[j],k,bDum);
      }
    }
    /* Walls */
    if (file_version >= 45) {
      do_int(ir->nwall);
      do_int(ir->wall_type);
      if (file_version >= 50)
	do_real(ir->wall_r_linpot);
      else
	ir->wall_r_linpot = -1;
      do_int(ir->wall_atomtype[0]);
      do_int(ir->wall_atomtype[1]);
      do_real(ir->wall_density[0]);
      do_real(ir->wall_density[1]);
      do_real(ir->wall_ewald_zfac);
    } else {
      ir->nwall = 0;
      ir->wall_type = 0;
      ir->wall_atomtype[0] = -1;
      ir->wall_atomtype[1] = -1;
      ir->wall_density[0] = 0;
      ir->wall_density[1] = 0;
      ir->wall_ewald_zfac = 3;
    }
    /* Cosine stuff for electric fields */
    for(j=0; (j<DIM); j++) {
      do_int  (ir->ex[j].n);
      do_int  (ir->et[j].n);
      if (bRead) {
	snew(ir->ex[j].a,  ir->ex[j].n);
	snew(ir->ex[j].phi,ir->ex[j].n);
	snew(ir->et[j].a,  ir->et[j].n);
	snew(ir->et[j].phi,ir->et[j].n);
      }
      ndo_real(ir->ex[j].a,  ir->ex[j].n,bDum);
      ndo_real(ir->ex[j].phi,ir->ex[j].n,bDum);
      ndo_real(ir->et[j].a,  ir->et[j].n,bDum);
      ndo_real(ir->et[j].phi,ir->et[j].n,bDum);
    }
    
    /* QMMM stuff */
    if(file_version>=39){
      do_int(ir->bQMMM);
      do_int(ir->QMMMscheme);
      do_real(ir->scalefactor);
      do_int(ir->opts.ngQM);
      if (bRead) {
        snew(ir->opts.QMmethod,    ir->opts.ngQM);
        snew(ir->opts.QMbasis,     ir->opts.ngQM);
        snew(ir->opts.QMcharge,    ir->opts.ngQM);
        snew(ir->opts.QMmult,      ir->opts.ngQM);
        snew(ir->opts.bSH,         ir->opts.ngQM);
        snew(ir->opts.CASorbitals, ir->opts.ngQM);
        snew(ir->opts.CASelectrons,ir->opts.ngQM);
        snew(ir->opts.SAon,        ir->opts.ngQM);
        snew(ir->opts.SAoff,       ir->opts.ngQM);
        snew(ir->opts.SAsteps,     ir->opts.ngQM);
        snew(ir->opts.bOPT,        ir->opts.ngQM);
        snew(ir->opts.bTS,         ir->opts.ngQM);
      }
      if (ir->opts.ngQM > 0) {
        ndo_int(ir->opts.QMmethod,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.QMbasis,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.QMcharge,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.QMmult,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.bSH,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.CASorbitals,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.CASelectrons,ir->opts.ngQM,bDum);
        ndo_real(ir->opts.SAon,ir->opts.ngQM,bDum);
        ndo_real(ir->opts.SAoff,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.SAsteps,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.bOPT,ir->opts.ngQM,bDum);
        ndo_int(ir->opts.bTS,ir->opts.ngQM,bDum);
      }
      /* end of QMMM stuff */
    }    
  }
}


static void do_harm(t_iparams *iparams,bool bRead)
{
  do_real(iparams->harmonic.rA);
  do_real(iparams->harmonic.krA);
  do_real(iparams->harmonic.rB);
  do_real(iparams->harmonic.krB);
}

void do_iparams(t_functype ftype,t_iparams *iparams,bool bRead, int file_version)
{
  int i;
  bool bDum;
  real VA[4],VB[4];
  
  if (!bRead)
    set_comment(interaction_function[ftype].name);
  switch (ftype) {
  case F_ANGLES:
  case F_G96ANGLES:
  case F_BONDS:
  case F_G96BONDS:
  case F_HARMONIC:
  case F_IDIHS:
    do_harm(iparams,bRead);
    if ((ftype == F_ANGRES || ftype == F_ANGRESZ) && bRead) {
      /* Correct incorrect storage of parameters */
      iparams->pdihs.phiB = iparams->pdihs.phiA;
      iparams->pdihs.cpB  = iparams->pdihs.cpA;
    }
    break;
  case F_FENEBONDS:
    do_real(iparams->fene.bm);
    do_real(iparams->fene.kb);
    break;
  case F_TABBONDS:
  case F_TABBONDSNC:
  case F_TABANGLES:
  case F_TABDIHS:
    do_real(iparams->tab.kA);
    do_int (iparams->tab.table);
    do_real(iparams->tab.kB);
    break;
  case F_CROSS_BOND_BONDS:
    do_real(iparams->cross_bb.r1e);
    do_real(iparams->cross_bb.r2e);
    do_real(iparams->cross_bb.krr);
    break;
  case F_CROSS_BOND_ANGLES:
    do_real(iparams->cross_ba.r1e);
    do_real(iparams->cross_ba.r2e);
    do_real(iparams->cross_ba.r3e);
    do_real(iparams->cross_ba.krt);
    break;
  case F_UREY_BRADLEY:
    do_real(iparams->u_b.theta);
    do_real(iparams->u_b.ktheta);
    do_real(iparams->u_b.r13);
    do_real(iparams->u_b.kUB);
    break;
  case F_QUARTIC_ANGLES:
    do_real(iparams->qangle.theta);
    ndo_real(iparams->qangle.c,5,bDum);
    break;
  case F_BHAM:
    do_real(iparams->bham.a);
    do_real(iparams->bham.b);
    do_real(iparams->bham.c);
    break;
  case F_MORSE:
    do_real(iparams->morse.b0);
    do_real(iparams->morse.cb);
    do_real(iparams->morse.beta);
    break;
  case F_CUBICBONDS:
    do_real(iparams->cubic.b0);
    do_real(iparams->cubic.kb);
    do_real(iparams->cubic.kcub);
    break;
  case F_CONNBONDS:
    break;
  case F_POLARIZATION:
    do_real(iparams->polarize.alpha);
    break;
  case F_WATER_POL:
    if (file_version < 31) 
      gmx_fatal(FARGS,"Old tpr files with water_polarization not supported. Make a new.");
    do_real(iparams->wpol.al_x);
    do_real(iparams->wpol.al_y);
    do_real(iparams->wpol.al_z);
    do_real(iparams->wpol.rOH);
    do_real(iparams->wpol.rHH);
    do_real(iparams->wpol.rOD);
    break;
  case F_THOLE_POL:
    do_real(iparams->thole.a);
    do_real(iparams->thole.alpha1);
    do_real(iparams->thole.alpha2);
    do_real(iparams->thole.rfac);
    break;
  case F_LJ:
    do_real(iparams->lj.c6);
    do_real(iparams->lj.c12);
    break;
  case F_LJ14:
    do_real(iparams->lj14.c6A);
    do_real(iparams->lj14.c12A);
    do_real(iparams->lj14.c6B);
    do_real(iparams->lj14.c12B);
    break;
  case F_LJC14_Q:
    do_real(iparams->ljc14.fqq);
    do_real(iparams->ljc14.qi);
    do_real(iparams->ljc14.qj);
    do_real(iparams->ljc14.c6);
    do_real(iparams->ljc14.c12);
    break;
  case F_LJC_PAIRS_NB:
    do_real(iparams->ljcnb.qi);
    do_real(iparams->ljcnb.qj);
    do_real(iparams->ljcnb.c6);
    do_real(iparams->ljcnb.c12);
    break;
  case F_PDIHS:
  case F_PIDIHS:
  case F_ANGRES:
  case F_ANGRESZ:
    do_real(iparams->pdihs.phiA);
    do_real(iparams->pdihs.cpA);
    if ((ftype == F_ANGRES || ftype == F_ANGRESZ) && file_version < 42) {
      /* Read the incorrectly stored multiplicity */
      do_real(iparams->harmonic.rB);
      do_real(iparams->harmonic.krB);
      iparams->pdihs.phiB = iparams->pdihs.phiA;
      iparams->pdihs.cpB  = iparams->pdihs.cpA;
    } else {
      do_real(iparams->pdihs.phiB);
      do_real(iparams->pdihs.cpB);
      do_int (iparams->pdihs.mult);
    }
    break;
  case F_DISRES:
    do_int (iparams->disres.label);
    do_int (iparams->disres.type);
    do_real(iparams->disres.low);
    do_real(iparams->disres.up1);
    do_real(iparams->disres.up2);
    do_real(iparams->disres.kfac);
    break;
  case F_ORIRES:
    do_int (iparams->orires.ex);
    do_int (iparams->orires.label);
    do_int (iparams->orires.power);
    do_real(iparams->orires.c);
    do_real(iparams->orires.obs);
    do_real(iparams->orires.kfac);
    break;
  case F_DIHRES:
    do_int (iparams->dihres.power);
    do_int (iparams->dihres.label);
    do_real(iparams->dihres.phi);
    do_real(iparams->dihres.dphi);
    do_real(iparams->dihres.kfac);
    break;
  case F_POSRES:
    do_rvec(iparams->posres.pos0A);
    do_rvec(iparams->posres.fcA);
    if (bRead && file_version < 27) {
      copy_rvec(iparams->posres.pos0A,iparams->posres.pos0B);
      copy_rvec(iparams->posres.fcA,iparams->posres.fcB);
    } else {
      do_rvec(iparams->posres.pos0B);
      do_rvec(iparams->posres.fcB);
    }
    break;
  case F_RBDIHS:
    ndo_real(iparams->rbdihs.rbcA,NR_RBDIHS,bDum);
    if(file_version>=25) 
      ndo_real(iparams->rbdihs.rbcB,NR_RBDIHS,bDum);
    break;
  case F_FOURDIHS:
    /* Fourier dihedrals are internally represented
     * as Ryckaert-Bellemans since those are faster to compute.
     */
    ndo_real(VA,NR_RBDIHS,bDum);
    ndo_real(VB,NR_RBDIHS,bDum);
    break;
  case F_CONSTR:
  case F_CONSTRNC:
    do_real(iparams->constr.dA);
    do_real(iparams->constr.dB);
    break;
  case F_SETTLE:
    do_real(iparams->settle.doh);
    do_real(iparams->settle.dhh);
    break;
  case F_VSITE2:
    do_real(iparams->vsite.a);
    break;
  case F_VSITE3:
  case F_VSITE3FD:
  case F_VSITE3FAD:
    do_real(iparams->vsite.a);
    do_real(iparams->vsite.b);
    break;
  case F_VSITE3OUT:
  case F_VSITE4FD: 
  case F_VSITE4FDN: 
    do_real(iparams->vsite.a);
    do_real(iparams->vsite.b);
    do_real(iparams->vsite.c);
    break;
  case F_VSITEN:
    do_int (iparams->vsiten.n);
    do_real(iparams->vsiten.a);
    break;
  default:
    gmx_fatal(FARGS,"unknown function type %d (%s) in %s line %d",
		
    		ftype,interaction_function[ftype].name,__FILE__,__LINE__);
  }
  if (!bRead)
    unset_comment();
}

static void do_ilist(t_ilist *ilist,bool bRead,char *name,int file_version)
{
  int  i,idum;
  bool bDum=TRUE;
  
  if (!bRead)
    set_comment(name);
  if (file_version < 44)
    for(i=0; i<MAXNODES; i++)
      do_int(idum);
  do_int (ilist->nr);
  if (bRead)
    snew(ilist->iatoms,ilist->nr);
  ndo_int(ilist->iatoms,ilist->nr,bDum);
  if (!bRead)
    unset_comment();
}

static void do_idef(t_idef *idef,bool bRead, int file_version)
{
  int i,j,k,renum[F_NRE];
  bool bDum=TRUE,bClear;
  
  do_int(idef->atnr);
  do_int(idef->nodeid);
  do_int(idef->ntypes);
  if (bRead && debug)
    fprintf(debug,"idef->atnr = %d, nodeid = %d, ntypes = %d\n",
	    idef->atnr,idef->nodeid,idef->ntypes);
  if (bRead) {
    snew(idef->functype,idef->ntypes);
    snew(idef->iparams,idef->ntypes);
  }
  /* Read/write all the function types */
  ndo_int(idef->functype,idef->ntypes,bDum);
  if (bRead && debug)
    pr_ivec(debug,0,"functype",idef->functype,idef->ntypes,TRUE);
    
  /* Check whether all these function types are supported by the code.
   * In practice the code is backwards compatible, which means that the
   * numbering may have to be altered from old numbering to new numbering
   */
  for (i=0; (i<idef->ntypes); i++) {
    if (bRead)
      /* Loop over file versions */
      for (k=0; (k<NFTUPD); k++)
	/* Compare the read file_version to the update table */
	if ((file_version < ftupd[k].fvnr) && 
	    (idef->functype[i] >= ftupd[k].ftype)) {
	  idef->functype[i] += 1;
	  if (debug) {
	    fprintf(debug,"Incrementing function type %d to %d (due to %s)\n",
		    i,idef->functype[i],
		    interaction_function[ftupd[k].ftype].longname);
	    fflush(debug);
	  }
	}
    
    do_iparams(idef->functype[i],&idef->iparams[i],bRead,file_version);
    if (bRead && debug)
      pr_iparams(debug,idef->functype[i],&idef->iparams[i]);
  }
  if (file_version >= 54)
    do_real(idef->fudgeQQ);
  
  for(j=0; (j<F_NRE); j++) {
    bClear = FALSE;
    if (bRead)
      for (k=0; k<NFTUPD; k++)
	if ((file_version < ftupd[k].fvnr) && (j == ftupd[k].ftype))
	  bClear = TRUE;
    if (bClear) {
      idef->il[j].nr = 0;
      idef->il[j].iatoms = NULL;
    } else
      do_ilist(&idef->il[j],bRead,interaction_function[j].name,file_version);
    if (bRead && gmx_debug_at)
      pr_ilist(debug,0,interaction_function[j].longname,
	       idef,&idef->il[j],TRUE);  }
  
  if (bRead) {
    /* We have not checked the sorting */
    idef->ilsort = ilsortUNKNOWN;
  }
}

static void do_block(t_block *block,bool bRead,int file_version)
{
  int  i,idum,dum_nra,*dum_a;
  bool bDum=TRUE;

  if (file_version < 44)
    for(i=0; i<MAXNODES; i++)
      do_int(idum);
  do_int (block->nr);
  if (file_version < 51)
    do_int (dum_nra);
  if (bRead) {
    block->nalloc_index = block->nr+1;
    snew(block->index,block->nalloc_index);
  }
  ndo_int(block->index,block->nr+1,bDum);

  if (file_version < 51 && dum_nra > 0) {
    snew(dum_a,dum_nra);
    ndo_int(dum_a,dum_nra,bDum);
    sfree(dum_a);
  }
}

static void do_blocka(t_blocka *block,bool bRead,int file_version)
{
  int  i,idum;
  bool bDum=TRUE;

  if (file_version < 44)
    for(i=0; i<MAXNODES; i++)
      do_int(idum);
  do_int (block->nr);
  do_int (block->nra);
  if (bRead) {
    block->nalloc_index = block->nr+1;
    snew(block->index,block->nalloc_index);
    block->nalloc_a = block->nra;
    snew(block->a,block->nalloc_a);
  }
  ndo_int(block->index,block->nr+1,bDum);
  ndo_int(block->a,block->nra,bDum);
}

static void do_atom(t_atom *atom,int ngrp,bool bRead, int file_version)
{ 
  int i,myngrp;
  
  do_real (atom->m);
  do_real (atom->q);
  do_real (atom->mB);
  do_real (atom->qB);
  do_ushort(atom->type);
  do_ushort(atom->typeB);
  do_int (atom->ptype);
  do_int (atom->resnr);
  if (file_version >= 52)
    do_int(atom->atomnumber);
  else if (bRead)
    atom->atomnumber = NOTSET;
  if (file_version < 23) 
    myngrp = 8;
  else if (file_version < 39) 
    myngrp = 9;
  else
    myngrp = ngrp;
  
  do_nuchar(atom->grpnr,myngrp);
  for(i=myngrp; (i<ngrp); i++)
    atom->grpnr[i] = 0;
  
}

static void do_grps(int ngrp,t_grps grps[],bool bRead, int file_version)
{
  int i,j,myngrp;
  bool bDum=TRUE;
  
  if (file_version < 23) 
    myngrp = 8;
  else if (file_version < 39) 
    myngrp = 9;
  else
    myngrp = ngrp;

  for(j=0; (j<ngrp); j++) {
    if (j<myngrp) {
      do_int (grps[j].nr);
      if (bRead)
	snew(grps[j].nm_ind,grps[j].nr);
      ndo_int(grps[j].nm_ind,grps[j].nr,bDum);
    }
    else {
      grps[j].nr = 1;
      snew(grps[j].nm_ind,grps[j].nr);
    }
  }
}

static void do_symstr(char ***nm,bool bRead,t_symtab *symtab)
{
  int ls;
  
  if (bRead) {
    do_int(ls);
    *nm = get_symtab_handle(symtab,ls);
  }
  else {
    ls = lookup_symtab(symtab,*nm);
    do_int(ls);
  }
}

static void do_strstr(int nstr,char ***nm,bool bRead,t_symtab *symtab)
{
  int  j;
  
  for (j=0; (j<nstr); j++) 
    do_symstr(&(nm[j]),bRead,symtab);
}

static void do_atoms(t_atoms *atoms,bool bRead,t_symtab *symtab, int file_version)
{
  int i;
  
  do_int(atoms->nr);
  do_int(atoms->nres);
  do_int(atoms->ngrpname);
  if (bRead) {
    snew(atoms->atom,atoms->nr);
    snew(atoms->atomname,atoms->nr);
    snew(atoms->atomtype,atoms->nr);
    snew(atoms->atomtypeB,atoms->nr);
    snew(atoms->resname,atoms->nres);
    snew(atoms->grpname,atoms->ngrpname);
    atoms->pdbinfo = NULL;
  }
  for(i=0; (i<atoms->nr); i++)
    do_atom(&atoms->atom[i],egcNR,bRead, file_version);
  do_strstr(atoms->nr,atoms->atomname,bRead,symtab);
  if (bRead && (file_version <= 20)) {
    for(i=0; i<atoms->nr; i++) {
      atoms->atomtype[i]  = put_symtab(symtab,"?");
      atoms->atomtypeB[i] = put_symtab(symtab,"?");
    }
  } else {
    do_strstr(atoms->nr,atoms->atomtype,bRead,symtab);
    do_strstr(atoms->nr,atoms->atomtypeB,bRead,symtab);
  }
  do_strstr(atoms->nres,atoms->resname,bRead,symtab);
  do_strstr(atoms->ngrpname,atoms->grpname,bRead,symtab);
  
  do_grps(egcNR,atoms->grps,bRead,file_version);
}

static void do_atomtypes(t_atomtypes *atomtypes,bool bRead,
			 t_symtab *symtab,int file_version)
{
  int i,j;
  bool bDum = TRUE;
  
  if (file_version > 25) {
    do_int(atomtypes->nr);
    j=atomtypes->nr;
    if (bRead) {
      snew(atomtypes->radius,j);
      snew(atomtypes->vol,j);
      snew(atomtypes->surftens,j);
      snew(atomtypes->atomnumber,j);
    }
    ndo_real(atomtypes->radius,j,bDum);
    ndo_real(atomtypes->vol,j,bDum);
    ndo_real(atomtypes->surftens,j,bDum);
    if(file_version >= 40)
    {
        ndo_int(atomtypes->atomnumber,j,bDum);
    }
  } else {
    /* File versions prior to 26 cannot do GBSA, 
     * so they dont use this structure 
     */
    atomtypes->nr = 0;
    atomtypes->radius = NULL;
    atomtypes->vol = NULL;
    atomtypes->surftens = NULL;
    atomtypes->atomnumber = NULL;
  }  
}

static void do_symtab(t_symtab *symtab,bool bRead)
{
  int i,nr;
  t_symbuf *symbuf;
  char buf[STRLEN];
  
  do_int(symtab->nr);
  nr     = symtab->nr;
  if (bRead) {
    snew(symtab->symbuf,1);
    symbuf = symtab->symbuf;
    symbuf->bufsize = nr;
    snew(symbuf->buf,nr);
    for (i=0; (i<nr); i++) {
      do_string(buf);
      symbuf->buf[i]=strdup(buf);
    }
  }
  else {
    symbuf = symtab->symbuf;
    while (symbuf!=NULL) {
      for (i=0; (i<symbuf->bufsize) && (i<nr); i++) 
	do_string(symbuf->buf[i]);
      nr-=i;
      symbuf=symbuf->next;
    }
    if (nr != 0)
      gmx_fatal(FARGS,"nr of symtab strings left: %d",nr);
  }
}

static void make_chain_identifiers(t_atoms *atoms,t_block *mols)
{
  int m,a,a0,a1;
  char c,chain;
#define CHAIN_MIN_ATOMS 15

  chain='A';
  for(m=0; m<mols->nr; m++) {
    a0=mols->index[m];
    a1=mols->index[m+1];
    if ((a1-a0 >= CHAIN_MIN_ATOMS) && (chain <= 'Z')) {
      c=chain;
      chain++;
    } else
      c=' ';
    for(a=a0; a<a1; a++)
      atoms->atom[a].chain=c;  
  }
  if (chain == 'B')
    for(a=0; a<atoms->nr; a++)
      atoms->atom[a].chain=' ';
}
  
static void do_top(t_topology *top,bool bRead, int file_version)
{
  int  i;
  t_blocka dumb;

  if (bRead)
    init_top(top);
  do_symtab(&(top->symtab),bRead);
  if (bRead && debug) 
    pr_symtab(debug,0,"symtab",&top->symtab);
  
  do_symstr(&(top->name),bRead,&(top->symtab));
  
  do_atoms (&(top->atoms),bRead,&(top->symtab), file_version);
  if (bRead && gmx_debug_at) 
    pr_atoms(debug,0,"atoms",&top->atoms,TRUE);

  /* This used to be in the atoms struct */
  do_blocka(&top->excls,bRead,file_version);

  do_atomtypes (&(top->atomtypes),bRead,&(top->symtab), file_version);
  if (bRead && debug) 
    pr_atomtypes(debug,0,"atomtypes",&top->atomtypes,TRUE);

  /* Debug statements are inside do_idef */    
  do_idef  (&(top->idef),bRead,file_version);
    
  do_block(&top->cgs,bRead,file_version);
  if (bRead && gmx_debug_at)
    pr_block(debug,0,"cgs",&top->cgs,TRUE);
  do_block(&top->mols,bRead,file_version);
  if (bRead && gmx_debug_at)
    pr_block(debug,0,"mols",&top->mols,TRUE);

  if (file_version < 51) {
    /* Here used to be the shake blocks */
    do_blocka(&dumb,bRead,file_version);
    if (dumb.nr > 0)
      sfree(dumb.index);
    if (dumb.nra > 0)
      sfree(dumb.a);
  }

  if (bRead) {
    close_symtab(&(top->symtab));
    make_chain_identifiers(&(top->atoms),&top->mols);
  }
}

/* If TopOnlyOK is TRUE then we can read even future versions
 * of tpx files, provided the file_generation hasn't changed.
 * If it is FALSE, we need the inputrecord too, and bail out
 * if the file is newer than the program.
 * 
 * The version and generation if the topology (see top of this file)
 * are returned in the two last arguments.
 * 
 * If possible, we will read the inputrec even when TopOnlyOK is TRUE.
 */
static void do_tpxheader(int fp,bool bRead,t_tpxheader *tpx, bool TopOnlyOK, int *file_version, int *file_generation)
{
  char  buf[STRLEN];
  bool  bDouble;
  int   precision;
  int   fver,fgen;
  gmx_fio_select(fp);
  gmx_fio_setdebug(fp,bDebugMode());
  
  /* NEW! XDR tpb file */
  precision = sizeof(real);
  if (bRead) {
    do_string(buf);
    if (strncmp(buf,"VERSION",7))
      gmx_fatal(FARGS,"Can not read file %s,\n"
		  "             this file is from a Gromacs version which is older than 2.0\n"
		  "             Make a new one with grompp or use a gro or pdb file, if possible",
		  gmx_fio_getname(fp));
    do_int(precision);
    bDouble = (precision == sizeof(double));
    if ((precision != sizeof(float)) && !bDouble)
      gmx_fatal(FARGS,"Unknown precision in file %s: real is %d bytes "
		  "instead of %d or %d",
		  gmx_fio_getname(fp),precision,sizeof(float),sizeof(double));
    gmx_fio_setprecision(fp,bDouble);
    fprintf(stderr,"Reading file %s, %s (%s precision)\n",
	    gmx_fio_getname(fp),buf,bDouble ? "double" : "single");
  }
  else {
    do_string(GromacsVersion());
    bDouble = (precision == sizeof(double));
    gmx_fio_setprecision(fp,bDouble);
    do_int(precision);
    fver = tpx_version;
    fgen = tpx_generation;
  }
  
  /* Check versions! */
  do_int(fver);
  
  if(fver>=26)
    do_int(fgen);
  else
    fgen=0;
 
  if(file_version!=NULL)
    *file_version = fver;
  if(file_version!=NULL)
    *file_generation = fgen;
   
  
  if ((fver <= tpx_incompatible_version) ||
      ((fver > tpx_version) && !TopOnlyOK) ||
      (fgen > tpx_generation))
    gmx_fatal(FARGS,"reading tpx file (%s) version %d with version %d program",
		gmx_fio_getname(fp),fver,tpx_version);
  
  do_section(eitemHEADER,bRead);
  do_int (tpx->natoms);
  if (fver >= 28)
    do_int(tpx->ngtc);
  else
    tpx->ngtc = 0;
  do_int (tpx->step);
  do_real(tpx->t);
  do_real(tpx->lambda);
  do_int (tpx->bIr);
  do_int (tpx->bTop);
  do_int (tpx->bX);
  do_int (tpx->bV);
  do_int (tpx->bF);
  do_int (tpx->bBox);

  if((fgen > tpx_generation)) {
    /* This can only happen if TopOnlyOK=TRUE */
    tpx->bIr=FALSE;
  }
}

static int do_tpx(int fp,bool bRead,int *step,real *t,
		  t_inputrec *ir,t_state *state,rvec *f,t_topology *top,
		  bool bXVallocated)
{
  t_tpxheader tpx;
  t_inputrec  dum_ir;
  t_topology  dum_top;
  bool        TopOnlyOK,bDum=TRUE;
  int         file_version,file_generation;
  int         i;
  rvec        *xptr,*vptr;
  int         ePBC;
  bool        bPeriodicMols;

  if (!bRead) {
    tpx.natoms = state->natoms;
    tpx.ngtc   = state->ngtc;
    tpx.step   = *step;
    tpx.t      = *t;
    tpx.lambda = state->lambda;
    tpx.bIr  = (ir       != NULL);
    tpx.bTop = (top      != NULL);
    tpx.bX   = (state->x != NULL);
    tpx.bV   = (state->v != NULL);
    tpx.bF   = (f        != NULL);
    tpx.bBox = TRUE;
  }
  
  TopOnlyOK = (ir==NULL);
  
  do_tpxheader(fp,bRead,&tpx,TopOnlyOK,&file_version,&file_generation);

  if (bRead) {
    *step         = tpx.step;
    *t            = tpx.t;
    state->flags  = 0;
    state->lambda = tpx.lambda;
    /* The init_state calls initialize the Nose-Hoover xi integrals to zero */
    if (bXVallocated) {
      xptr = state->x;
      vptr = state->v;
      init_state(state,0,tpx.ngtc);
      state->natoms = tpx.natoms;
      state->nalloc = tpx.natoms;
      state->x = xptr;
      state->v = vptr;
    } else {
      init_state(state,tpx.natoms,tpx.ngtc);
    }
  }

#define do_test(b,p) if (bRead && (p!=NULL) && !b) gmx_fatal(FARGS,"No %s in %s",#p,gmx_fio_getname(fp)) 

  do_test(tpx.bBox,state->box);
  do_section(eitemBOX,bRead);
  if (tpx.bBox) {
    ndo_rvec(state->box,DIM);
    if (file_version >= 51) {
      ndo_rvec(state->box_rel,DIM);
    } else {
      /* We initialize box_rel after reading the inputrec */
      clear_mat(state->box_rel);
    }
    if (file_version >= 28) {
      ndo_rvec(state->boxv,DIM);
      if (file_version < 56) {
	matrix mdum;
	ndo_rvec(mdum,DIM);
      }
    }
  }
  
  if (state->ngtc > 0 && file_version >= 28) {
    real *dumv;
    ndo_real(state->nosehoover_xi,state->ngtc,bDum);
    snew(dumv,state->ngtc);
    /* These used to be the Berendsen tcoupl_lambda's */
    ndo_real(dumv,state->ngtc,bDum);
    sfree(dumv);
  }

  /* Prior to tpx version 26, the inputrec was here.
   * I moved it to enable partial forward-compatibility
   * for analysis/viewer programs.
   */
  if(file_version<26) {
    do_test(tpx.bIr,ir);
    do_section(eitemIR,bRead);
    if (tpx.bIr) {
      if (ir) {
	do_inputrec(ir,bRead,file_version,top ? &top->idef.fudgeQQ : NULL);
	if (bRead && debug) 
	  pr_inputrec(debug,0,"inputrec",ir,FALSE);
      }
      else {
	do_inputrec(&dum_ir,bRead,file_version,top ? &top->idef.fudgeQQ :NULL);
	if (bRead && debug) 
	  pr_inputrec(debug,0,"inputrec",&dum_ir,FALSE);
	done_inputrec(&dum_ir);
      }
      
    }
  }
  
  do_test(tpx.bTop,top);
  do_section(eitemTOP,bRead);
  if (tpx.bTop) {
    if (top)
      do_top(top,bRead, file_version);
    else {
      do_top(&dum_top,bRead,file_version);
      done_top(&dum_top);
    }
  }
  do_test(tpx.bX,state->x);  
  do_section(eitemX,bRead);
  if (tpx.bX) {
    if (bRead) {
      state->flags |= (1<<estX);
      if (!bXVallocated)
	snew(state->x,state->nalloc);
    }
    ndo_rvec(state->x,state->natoms);
  }
  
  do_test(tpx.bV,state->v);
  do_section(eitemV,bRead);
  if (tpx.bV) {
    if (bRead) {
      state->flags |= (1<<estV);
      if (!bXVallocated)
	snew(state->v,state->nalloc);
    }
    ndo_rvec(state->v,state->natoms);
  }

  do_test(tpx.bF,f);
  do_section(eitemF,bRead);
  if (tpx.bF) ndo_rvec(f,state->natoms);

  /* Starting with tpx version 26, we have the inputrec
   * at the end of the file, so we can ignore it 
   * if the file is never than the software (but still the
   * same generation - see comments at the top of this file.
   *
   * 
   */
  ePBC = -1;
  bPeriodicMols = FALSE;
  if (file_version >= 26) {
    do_test(tpx.bIr,ir);
    do_section(eitemIR,bRead);
    if (tpx.bIr) {
      if (file_version >= 53) {
	/* Removed the pbc info from do_inputrec, since we always want it */
	if (!bRead) {
	  ePBC          = ir->ePBC;
	  bPeriodicMols = ir->bPeriodicMols;
	}
	do_int(ePBC);
	do_int(bPeriodicMols);
      }
      if (file_generation <= tpx_generation && ir) {
	do_inputrec(ir,bRead,file_version,top ? &top->idef.fudgeQQ : NULL);
	if (bRead && debug) 
	  pr_inputrec(debug,0,"inputrec",ir,FALSE);
	if (file_version < 51)
	  set_box_rel(ir,state);
	if (file_version < 53) {
	  ePBC          = ir->ePBC;
	  bPeriodicMols = ir->bPeriodicMols;
	}
      }
      if (bRead && ir && file_version >= 53) {
	/* We need to do this after do_inputrec, since that initializes ir */
	ir->ePBC          = ePBC;
	ir->bPeriodicMols = bPeriodicMols;
      }
    }
  }

  if (bRead && tpx.bIr && ir) {
    /* Check for perturbed bonded interactions */
    gmx_analyze_ilist_fe(&top->idef,ir);

    if (state->ngtc == 0) {
      /* Reading old version without tcoupl state data: set it */
      init_gtc_state(state,ir->opts.ngtc);
    }
  }

  return ePBC;
}

/************************************************************
 *
 *  The following routines are the exported ones
 *
 ************************************************************/

int open_tpx(char *fn,char *mode)
{
  return gmx_fio_open(fn,mode);
}    
 
void close_tpx(int fp)
{
  gmx_fio_close(fp);
}

void read_tpxheader(char *fn,t_tpxheader *tpx, bool TopOnlyOK,int *file_version, int *file_generation)
{
  int fp;

#ifdef HAVE_LIBXML2
  if (fn2ftp(fn) == efXML) {
    gmx_fatal(FARGS,"read_tpxheader called with filename %s",fn);
  }
  else {
#endif
    fp = open_tpx(fn,"r");
    do_tpxheader(fp,TRUE,tpx,TopOnlyOK,file_version,file_generation);
    close_tpx(fp);
#ifdef HAVE_LIBXML2
  }
#endif
}

void write_tpx_state(char *fn,int step,real t,
		     t_inputrec *ir,t_state *state,t_topology *top)
{
  int fp;

#ifdef HAVE_LIBXML2
  if (fn2ftp(fn) == efXML)
    write_xml(fn,*top->name,ir,state->box,state->natoms,
	      state->x,state->v,NULL,1,&top->atoms,&top->idef);
  else {
#endif
    fp = open_tpx(fn,"w");
    do_tpx(fp,FALSE,&step,&t,ir,state,NULL,top,FALSE);
    close_tpx(fp);
#ifdef HAVE_LIBXML2
  }
#endif
}

void read_tpx_state(char *fn,int *step,real *t,
		    t_inputrec *ir,t_state *state,rvec *f,t_topology *top)
{
  int fp;

  fp = open_tpx(fn,"r");
  do_tpx(fp,TRUE,step,t,ir,state,f,top,FALSE);
  close_tpx(fp);
}

int read_tpx(char *fn,int *step,real *t,real *lambda,
	     t_inputrec *ir, matrix box,int *natoms,
	     rvec *x,rvec *v,rvec *f,t_topology *top)
{
  int fp;
  t_state state;
  int ePBC;

#ifdef HAVE_LIBXML2
  if (fn2ftp(fn) == efXML) {
    int  i;
    rvec *xx=NULL,*vv=NULL,*ff=NULL;
    
    read_xml(fn,step,t,lambda,ir,box,natoms,&xx,&vv,&ff,top);
    ePBC = -1;
    for(i=0; (i<*natoms); i++) {
      if (xx) copy_rvec(xx[i],x[i]);
      if (vv) copy_rvec(vv[i],v[i]);
      if (ff) copy_rvec(ff[i],f[i]);
    }
    if (xx) sfree(xx);
    if (vv) sfree(vv);
    if (ff) sfree(ff);
  }
  else {
#endif
    state.x = x;
    state.v = v;
    fp = open_tpx(fn,"r");
    ePBC = do_tpx(fp,TRUE,step,t,ir,&state,f,top,TRUE);
    close_tpx(fp);
    *natoms = state.natoms;
    if (lambda) 
      *lambda = state.lambda;
    if (box) 
      copy_mat(state.box,box);
    state.x = NULL;
    state.v = NULL;
    done_state(&state);
#ifdef HAVE_LIBXML2
  }
#endif

  return ePBC;
}

bool fn2bTPX(char *file)
{
  switch (fn2ftp(file)) {
  case efTPR:
  case efTPB:
  case efTPA:
    return TRUE;
  default:
    return FALSE;
  }
}

bool read_tps_conf(char *infile,char *title,t_topology *top,int *ePBC,
		   rvec **x,rvec **v,matrix box,bool bMass)
{
  t_tpxheader  header;
  real         t,lambda;
  int          natoms,step,i,version,generation;
  bool         bTop,bXNULL;
  gmx_atomprop_t aps;
  
  bTop = fn2bTPX(infile);
  *ePBC = -1;
  if (bTop) {
    read_tpxheader(infile,&header,TRUE,&version,&generation);
    if (x)
      snew(*x,header.natoms);
    if (v)
      snew(*v,header.natoms);
    *ePBC = read_tpx(infile,&step,&t,&lambda,NULL,box,&natoms,
		     (x==NULL) ? NULL : *x,(v==NULL) ? NULL : *v,NULL,top);
    strcpy(title,*top->name);
  }
  else {
    get_stx_coordnum(infile,&natoms);
    init_t_atoms(&top->atoms,natoms,FALSE);
    bXNULL = (x == NULL);
    snew(*x,natoms);
    if (v)
      snew(*v,natoms);
    read_stx_conf(infile,title,&top->atoms,*x,(v==NULL) ? NULL : *v,ePBC,box);
    if (bXNULL) {
      sfree(*x);
      x = NULL;
    }
    if (bMass) {
      aps = gmx_atomprop_init();
      for(i=0; (i<natoms); i++)
	if (!gmx_atomprop_query(aps,epropMass,
				*top->atoms.resname[top->atoms.atom[i].resnr],
				*top->atoms.atomname[i],
				&(top->atoms.atom[i].m))) {
	  if (debug) 
	    fprintf(debug,"Can not find mass for atom %s %d %s, setting to 1\n",
		    *top->atoms.resname[top->atoms.atom[i].resnr],
		    top->atoms.atom[i].resnr+1,
		    *top->atoms.atomname[i]);
	}
      gmx_atomprop_destroy(aps);
    }
    top->idef.ntypes=-1;
  }

  return bTop;
}
