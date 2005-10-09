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
 * Gallium Rubidium Oxygen Manganese Argon Carbon Silicon
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <ctype.h>
#include <stdlib.h>
#include <limits.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "physics.h"
#include "names.h"
#include "fatal.h"
#include "macros.h"
#include "index.h"
#include "symtab.h"
#include "string2.h"
#include "readinp.h"
#include "readir.h" 
#include "toputil.h"
#include "index.h"
#include "network.h"
#include "vec.h"
#include "pbc.h"

#define MAXPTR 254
#define NOGID  255

/* Resource parameters 
 * Do not change any of these until you read the instruction
 * in readinp.h. Some cpp's do not take spaces after the backslash
 * (like the c-shell), which will give you a very weird compiler
 * message.
 */

static char tcgrps[STRLEN],tau_t[STRLEN],ref_t[STRLEN],
  acc[STRLEN],accgrps[STRLEN],freeze[STRLEN],frdim[STRLEN],
  energy[STRLEN],user1[STRLEN],user2[STRLEN],vcm[STRLEN],xtc_grps[STRLEN],
  orirefitgrp[STRLEN],egptable[STRLEN],egpexcl[STRLEN],deform[STRLEN],
  QMMM[STRLEN];
static char anneal[STRLEN],anneal_npoints[STRLEN],
  anneal_time[STRLEN],anneal_temp[STRLEN];
static char QMmethod[STRLEN],QMbasis[STRLEN],QMcharge[STRLEN],QMmult[STRLEN],
  bSH[STRLEN],CASorbitals[STRLEN], CASelectrons[STRLEN],SAon[STRLEN],
  SAoff[STRLEN],SAsteps[STRLEN],bTS[STRLEN],bOPT[STRLEN]; 
static char efield_x[STRLEN],efield_xt[STRLEN],efield_y[STRLEN],
  efield_yt[STRLEN],efield_z[STRLEN],efield_zt[STRLEN];

void init_ir(t_inputrec *ir, t_gromppopts *opts)
{
  snew(opts->title,STRLEN);
  snew(opts->cpp,STRLEN); 
  snew(opts->include,STRLEN); 
  snew(opts->define,STRLEN);
  snew(opts->SolventOpt,STRLEN);
}

static void _low_check(bool b,char *s,int *n)
{
  if (b) {
    fprintf(stderr,"ERROR: %s\n",s);
    (*n)++;
  }
}

void check_ir(t_inputrec *ir, t_gromppopts *opts,int *nerror)
/* Check internal consistency */
{
  /* Strange macro: first one fills the err_buf, and then one can check 
   * the condition, which will print the message and increase the error
   * counter.
   */
#define CHECK(b) _low_check(b,err_buf,nerror)
  char err_buf[256];

  /* TPI STUFF */
  if (ir->eI == eiTPI) {
    sprintf(err_buf,"TPI does not work with pbc = %s",epbc_names[epbcNONE]);
    CHECK(ir->ePBC == epbcNONE);
    sprintf(err_buf,"TPI only works with ns = %s",ens_names[ensGRID]);
    CHECK(ir->ns_type != ensGRID);
    sprintf(err_buf,"with TPI nstlist should be larger than zero");
    CHECK(ir->nstlist <= 0);
    sprintf(err_buf,"TPI does not work with full electrostatics");
    CHECK(EEL_FULL(ir->coulombtype));
  }

  /* SHAKE / LINCS */
  sprintf(err_buf,"constraints with Conjugate Gradients not implemented");
  CHECK((opts->nshake > 0) && (ir->eI == eiCG));

  if ( (opts->nshake > 0) && (opts->bMorse) ) {
    sprintf(warn_buf,
	    "Using morse bond-potentials while constraining bonds is useless");
    warning(NULL);
  }
  
  sprintf(err_buf,"shake_tol must be > 0 instead of %g while using shake",
	  ir->shake_tol);
  CHECK(((ir->shake_tol <= 0.0) && (opts->nshake>0) && 
	 (ir->eConstrAlg == estSHAKE)));
     
     
  /* VACUUM STUFF */
  if (ir->ePBC == epbcNONE) {
    if (ir->epc != epcNO) {
      warning("Turning off pressure coupling for vacuum system");
      ir->epc = epcNO;
    }

    if (ir->ns_type != ensSIMPLE) {
      sprintf(warn_buf,"Can only use nstype=%s with pbc=%s, setting nstype "
	      "to %s\n",
	      ens_names[ensSIMPLE],epbc_names[epbcNONE],ens_names[ensSIMPLE]);
      warning(NULL);
      ir->ns_type = ensSIMPLE;
    }
    if (ir->eDispCorr != edispcNO) {
      warning("Can not have long-range dispersion correction without PBC,"
	      " turned off.");
      ir->eDispCorr = edispcNO;
    }
  }

  if (ir->rlist == 0.0) {
    sprintf(err_buf,"can only have neighborlist cut-off zero (=infinite)\n"
	    "with coulombtype = %s or coulombtype = %s"
	    "and simple neighborsearch\n"
	    "without periodic boundary conditions (pbc = %s) and\n"
	    "rcoulomb and rvdw set to zero",
	    eel_names[eelCUT],eel_names[eelUSER],epbc_names[epbcNONE]);
    CHECK(((ir->coulombtype != eelCUT) && (ir->coulombtype != eelUSER))
	  || (ir->ns_type != ensSIMPLE) || 
	  (ir->ePBC     != epbcNONE) || 
	  (ir->rcoulomb != 0.0)      || (ir->rvdw != 0.0));
  }

  /* COMM STUFF */
  if ((ir->ePBC != epbcNONE) && (ir->nstcomm < 0))
    warning("Removing the rotation around the center of mass in a periodic system (this is not a problem when you have only one molecule).");
    
  if ((ir->eI == eiMD) && (ir->ePBC == epbcNONE)) {
    if ((ir->nstcomm == 0) || (ir->comm_mode != ecmANGULAR))
      warning("Tumbling and or flying ice-cubes: We are not removing rotation around center of mass in a non-periodic system. You should set comm_mode = ANGULAR.");
  }
  
  sprintf(err_buf,"Free-energy not implemented for Ewald and PPPM");
  CHECK((ir->coulombtype==eelEWALD || ir->coulombtype==eelPPPM)
	&& (ir->efep!=efepNO));
  
  sprintf(err_buf,"Domain decomposition can only be used with grid NS");
  CHECK(ir->bDomDecomp && (ir->ns_type == ensSIMPLE));
  sprintf(err_buf,"Twin-range neighbour searching (NS) with simple NS"
	  " algorithm not implemented");
  CHECK(((ir->rcoulomb > ir->rlist) || (ir->rvdw > ir->rlist)) 
	&& (ir->ns_type == ensSIMPLE));
  
  /* PRESSURE COUPLING */
  if(ir->epc == epcISOTROPIC) {
    ir->epc = epcBERENDSEN;
    fprintf(stderr,"Note: Old option for pressure coupling given: "
	           "changing \"Isotropic\" to \"Berendsen\"\n"); 
  }
  
  if (ir->epc != epcNO) {
    sprintf(err_buf,"tau_p must be > 0 instead of %g\n",ir->tau_p);
    CHECK(ir->tau_p <= 0);
       
    sprintf(err_buf,"compressibility must be > 0 when using pressure" 
	    " coupling %s\n",EPCOUPLTYPE(ir->epc));
    CHECK(ir->compress[XX][XX] < 0 || ir->compress[YY][YY] < 0 || 
	  ir->compress[ZZ][ZZ] < 0 || 
	  (trace(ir->compress) == 0 && ir->compress[YY][XX] <= 0 &&
	   ir->compress[ZZ][XX] <= 0 && ir->compress[ZZ][YY] <= 0));
    
    sprintf(err_buf,"pressure coupling with PPPM not implemented, use PME");
    CHECK(ir->coulombtype == eelPPPM);
     
  } else if (ir->coulombtype == eelPPPM) {
    sprintf(warn_buf,"The pressure with PPPM is incorrect, if you need the pressure use PME");
    warning(NULL);
  }
  
  /* TEMPERATURE COUPLING */
  if(ir->etc == etcYES) {
    ir->etc = etcBERENDSEN;
    fprintf(stderr,"Note: Old option for temperature coupling given: "
	           "changing \"yes\" to \"Berendsen\"\n");
  }
  
  if((ir->etc==etcNOSEHOOVER || ir->etc==etcANDERSEN || ir->etc==etcANDERSENINTERVAL ) 
     && ir->epc==epcBERENDSEN) {
    sprintf(warn_buf,"Using Berendsen pressure coupling invalidates the "
	    "true ensemble for the thermostat");
    warning(NULL);
  }

  /* ELECTROSTATICS */
  /* More checks are in triple check (grompp.c) */
  if (EEL_RF(ir->coulombtype) && ir->epsilon_rf==1 && ir->epsilon_r!=1) {
    sprintf(warn_buf,"epsilon_r = %g and epsilon_rf = 1 with reaction field, assuming old format and exchanging epsilon_r and epsilon_rf",ir->epsilon_r);
    warning(NULL);
    ir->epsilon_rf = ir->epsilon_r;
    ir->epsilon_r  = 1.0;
  }
  
  sprintf(err_buf,"epsilon_r must be >= 0 instead of %g\n",ir->epsilon_r);
  CHECK(ir->epsilon_r < 0);
  
  if (EEL_RF(ir->coulombtype)) {
    /* reaction field (at the cut-off) */
    
    sprintf(err_buf,"epsilon_rf must be >= epsilon_r");
    CHECK((ir->epsilon_rf < ir->epsilon_r && ir->epsilon_rf != 0) ||
	  (ir->epsilon_r == 0));
    if (ir->epsilon_rf == ir->epsilon_r) {
      sprintf(warn_buf,"Using epsilon_rf = epsilon_r with %s does not make sense",
	      eel_names[ir->coulombtype]);
      warning(NULL);
    }
  }
  /* Allow rlist>rcoulomb for tabulated long range stuff. This just
   * means the interaction is zero outside rcoulomb, but it helps to
   * provide accurate energy conservation.
   */
  if ((ir->coulombtype == eelSHIFT) || (ir->coulombtype == eelSWITCH) || (ir->coulombtype == eelENCADSHIFT) )
  {
    sprintf(err_buf,"With coulombtype = %s rcoulomb_switch must be < rcoulomb",
	    eel_names[ir->coulombtype]);
    CHECK(ir->rcoulomb_switch >= ir->rcoulomb);
  } else if (EEL_RF(ir->coulombtype)) {
    sprintf(err_buf,"With coulombtype = %s, rcoulomb must be >= rlist",eel_names[ir->coulombtype]);
    CHECK(ir->rlist > ir->rcoulomb);
  }

  if ((ir->coulombtype == eelPME) ||
      (ir->coulombtype == eelPMEUSER)) {
    if ((ir->pme_order < 4) || ((ir->pme_order % 2) == 1)) {
      if (ir->pme_order < 4)
	ir->pme_order = 4;
      else if ((ir->pme_order % 2) == 1)
	ir->pme_order++;
      sprintf(err_buf,"pme_order should be even and at least 4, modified to %d",
	      ir->pme_order);
      warning(NULL);
    }
  }

  if ((ir->vdwtype == evdwSWITCH) || (ir->vdwtype == evdwSHIFT) || (ir->vdwtype == evdwENCADSHIFT) )
  {
    sprintf(err_buf,"With vdwtype = %s rvdw_switch must be < rvdw",
	    evdw_names[ir->vdwtype]);
    CHECK(ir->rvdw_switch >= ir->rvdw);
  } else {
    sprintf(err_buf,"With vdwtype = %s,rvdw must be >= rlist",evdw_names[ir->vdwtype]);
    CHECK(ir->rlist > ir->rvdw);
  }
  if ((((ir->coulombtype == eelSHIFT) || (ir->coulombtype == eelSWITCH))
       && (ir->rlist == ir->rcoulomb)) ||
      (((ir->vdwtype == evdwSWITCH) || (ir->vdwtype == evdwSHIFT))
       && (ir->rlist == ir->rvdw))) {
    sprintf(warn_buf,"For energy conservation with switch/shift potentials, rlist should be 0.1 to 0.3 nm larger than rcoulomb/rvdw.");
    warning(NULL);
  }
  
  if(ir->eI == eiLBFGS && (ir->coulombtype==eelCUT || ir->vdwtype==evdwCUT)) {
    sprintf(warn_buf,"For efficient BFGS minimization, use switch/shift/pme instead of cut-off.");
    warning(NULL);
  }

  if(ir->eI == eiLBFGS && ir->nbfgscorr<=0) {
    sprintf(warn_buf,"Using L-BFGS with nbfgscorr<=0 just gets you steepest descent.");
    warning(NULL);
  }
}

static int str_nelem(char *str,int maxptr,char *ptr[])
{
  int  np=0;
  char *copy0,*copy;
  
  copy0=strdup(str); 
  copy=copy0;
  ltrim(copy);
  while (*copy != '\0') {
    if (np >= maxptr)
      gmx_fatal(FARGS,"Too many groups on line: '%s' (max is %d)",
		  str,maxptr);
    if (ptr) 
      ptr[np]=copy;
    np++;
    while ((*copy != '\0') && !isspace(*copy))
      copy++;
    if (*copy != '\0') {
      *copy='\0';
      copy++;
    }
    ltrim(copy);
  }
  if (ptr == NULL)
    sfree(copy0);

  return np;
}


void get_ir(char *mdparin,char *mdparout,
	    t_inputrec *ir,t_gromppopts *opts,int *nerror)
{
  char      *dumstr[2];
  double    dumdub[2][6];
  t_inpfile *inp;
  char      *tmp;
  int       i,j,m,ninp;
  char      dummy[STRLEN];
  
  inp=read_inpfile(mdparin,&ninp);

  snew(dumstr[0],STRLEN);
  snew(dumstr[1],STRLEN);

  CCTYPE ("VARIOUS PREPROCESSING OPTIONS");
  STYPE ("title",	opts->title,	NULL);
  STYPE ("cpp",		opts->cpp,	"/usr/bin/cpp");
  STYPE ("include",	opts->include,	NULL);
  STYPE ("define",	opts->define,	NULL);
    
  CCTYPE ("RUN CONTROL PARAMETERS");
  EETYPE("integrator",  ir->eI,         ei_names, nerror, TRUE);
  CTYPE ("Start time and timestep in ps");
  RTYPE ("tinit",	ir->init_t,	0.0);
  RTYPE ("dt",		ir->delta_t,	0.001);
  ITYPE ("nsteps",      ir->nsteps,     0);
  CTYPE ("For exact run continuation or redoing part of a run");
  ITYPE ("init_step",   ir->init_step,  0);
  CTYPE ("mode for center of mass motion removal");
  EETYPE("comm-mode",   ir->comm_mode,  ecm_names, nerror, TRUE);
  CTYPE ("number of steps for center of mass motion removal");
  ITYPE ("nstcomm",	ir->nstcomm,	1);
  CTYPE ("group(s) for center of mass motion removal");
  STYPE ("comm-grps",   vcm,            NULL);
  
  CCTYPE ("LANGEVIN DYNAMICS OPTIONS");
  CTYPE ("Friction coefficient (amu/ps) and random seed");
  RTYPE ("bd-fric",     ir->bd_fric,    0.0);
  ITYPE ("ld-seed",     ir->ld_seed,    1993);
  
  /* Em stuff */
  CCTYPE ("ENERGY MINIMIZATION OPTIONS");
  CTYPE ("Force tolerance and initial step-size");
  RTYPE ("emtol",       ir->em_tol,     10.0);
  RTYPE ("emstep",      ir->em_stepsize,0.01);
  CTYPE ("Max number of iterations in relax_shells");
  ITYPE ("niter",       ir->niter,      20);
  CTYPE ("Step size (ps^2) for minimization of flexible constraints");
  RTYPE ("fcstep",      ir->fc_stepsize, 0);
  CTYPE ("Frequency of steepest descents steps when doing CG");
  ITYPE ("nstcgsteep",	ir->nstcgsteep,	1000);
  ITYPE ("nbfgscorr",   ir->nbfgscorr,  10); 

  /* Output options */
  CCTYPE ("OUTPUT CONTROL OPTIONS");
  CTYPE ("Output frequency for coords (x), velocities (v) and forces (f)");
  ITYPE ("nstxout",	ir->nstxout,	100);
  ITYPE ("nstvout",	ir->nstvout,	100);
  ITYPE ("nstfout",	ir->nstfout,	0);
  CTYPE ("Checkpointing helps you continue after crashes");
  ITYPE ("nstcheckpoint",  ir->nstcheckpoint,	1000);
  CTYPE ("Output frequency for energies to log file and energy file");
  ITYPE ("nstlog",	ir->nstlog,	100);
  ITYPE ("nstenergy",   ir->nstenergy,  100);
  CTYPE ("Output frequency and precision for xtc file");
  ITYPE ("nstxtcout",   ir->nstxtcout,  0);
  RTYPE ("xtc-precision",ir->xtcprec,   1000.0);
  CTYPE ("This selects the subset of atoms for the xtc file. You can");
  CTYPE ("select multiple groups. By default all atoms will be written.");
  STYPE ("xtc-grps",    xtc_grps,       NULL);
  CTYPE ("Selection of energy groups");
  STYPE ("energygrps",  energy,         NULL);

  /* Neighbor searching */  
  CCTYPE ("NEIGHBORSEARCHING PARAMETERS");
  CTYPE ("nblist update frequency");
  ITYPE ("nstlist",	ir->nstlist,	10);
  CTYPE ("ns algorithm (simple or grid)");
  EETYPE("ns-type",     ir->ns_type,    ens_names, nerror, TRUE);
  /* set ndelta to the optimal value of 2 */
  ir->ndelta = 2;
  CTYPE ("Periodic boundary conditions: xyz (default), no (vacuum)");
  CTYPE ("or full (infinite systems only)");
  EETYPE("pbc",         ir->ePBC,       epbc_names, nerror, TRUE);
  CTYPE ("nblist cut-off");
  RTYPE ("rlist",	ir->rlist,	1.0);
  EETYPE("domain-decomposition",ir->bDomDecomp, yesno_names, nerror, TRUE);
  
  /* Electrostatics */
  CCTYPE ("OPTIONS FOR ELECTROSTATICS AND VDW");
  CTYPE ("Method for doing electrostatics");
  EETYPE("coulombtype",	ir->coulombtype,    eel_names, nerror, TRUE);
  CTYPE ("cut-off lengths");
  RTYPE ("rcoulomb-switch",	ir->rcoulomb_switch,	0.0);
  RTYPE ("rcoulomb",	ir->rcoulomb,	1.0);
  CTYPE ("Relative dielectric constant for the medium and the reaction field");
  RTYPE ("epsilon_r",   ir->epsilon_r,  1.0);
  RTYPE ("epsilon_rf",  ir->epsilon_rf, 1.0);
  CTYPE ("Method for doing Van der Waals");
  EETYPE("vdw-type",	ir->vdwtype,    evdw_names, nerror, TRUE);
  CTYPE ("cut-off lengths");
  RTYPE ("rvdw-switch",	ir->rvdw_switch,	0.0);
  RTYPE ("rvdw",	ir->rvdw,	1.0);
  CTYPE ("Apply long range dispersion corrections for Energy and Pressure");
  EETYPE("DispCorr",    ir->eDispCorr,  edispc_names, nerror, TRUE);
  CTYPE ("Extension of the potential lookup tables beyond the cut-off");
  RTYPE ("table-extension", ir->tabext, 1.0);
  CTYPE ("Seperate tables between energy group pairs");
  STYPE ("energygrp_table", egptable,   NULL);
  CTYPE ("Spacing for the PME/PPPM FFT grid");
  RTYPE ("fourierspacing", opts->fourierspacing,0.12);
  CTYPE ("FFT grid size, when a value is 0 fourierspacing will be used");
  ITYPE ("fourier_nx",  ir->nkx,         0);
  ITYPE ("fourier_ny",  ir->nky,         0);
  ITYPE ("fourier_nz",  ir->nkz,         0);
  CTYPE ("EWALD/PME/PPPM parameters");
  ITYPE ("pme_order",   ir->pme_order,   4);
  RTYPE ("ewald_rtol",  ir->ewald_rtol, 0.00001);
  EETYPE("ewald_geometry", ir->ewald_geometry, eewg_names, nerror,  TRUE);
  RTYPE ("epsilon_surface", ir->epsilon_surface, 0.0);
  EETYPE("optimize_fft",ir->bOptFFT,  yesno_names, nerror, TRUE);

  CCTYPE ("GENERALIZED BORN ELECTROSTATICS"); 
  CTYPE ("Algorithm for calculating Born radii");
  EETYPE("gb_algorithm", ir->gb_algorithm, egb_names, nerror, TRUE);
  CTYPE ("Frequency of calculating the Born radii inside rlist");
  ITYPE ("nstgbradii", ir->nstgbradii, 1);
  CTYPE ("Cutoff for Born radii calculation; the contribution from atoms");
  CTYPE ("between rlist and rgbradii is updated every nstlist steps");
  RTYPE ("rgbradii",  ir->rgbradii, 2.0);
  CTYPE ("Salt concentration in M for Generalized Born models");
  RTYPE ("gb_saltconc",  ir->gb_saltconc, 0.0); 

  CCTYPE("IMPLICIT SOLVENT (for use with Generalized Born electrostatics)");
  EETYPE("implicit_solvent", ir->implicit_solvent, eis_names, nerror, TRUE);
  
  /* Coupling stuff */
  CCTYPE ("OPTIONS FOR WEAK COUPLING ALGORITHMS");
  CTYPE ("Temperature coupling");
  EETYPE("tcoupl",	ir->etc,        etcoupl_names, nerror, TRUE);
  CTYPE ("Groups to couple separately");
  STYPE ("tc-grps",     tcgrps,         NULL);
  CTYPE ("Time constant (ps) and reference temperature (K)");
  STYPE ("tau-t",	tau_t,		NULL);
  STYPE ("ref-t",	ref_t,		NULL);
  CTYPE ("Pressure coupling");
  EETYPE("Pcoupl",	ir->epc,        epcoupl_names, nerror, TRUE);
  EETYPE("Pcoupltype",	ir->epct,       epcoupltype_names, nerror, TRUE);
  CTYPE ("Time constant (ps), compressibility (1/bar) and reference P (bar)");
  RTYPE ("tau-p",	ir->tau_p,	1.0);
  STYPE ("compressibility",	dumstr[0],	NULL);
  STYPE ("ref-p",       dumstr[1],      NULL);
  
  CTYPE ("Random seed for Andersen thermostat");
  ITYPE ("andersen_seed", ir->andersen_seed, 815131);

  /* QMMM */
  CCTYPE ("OPTIONS FOR QMMM calculations");
  EETYPE("QMMM", ir->bQMMM, yesno_names, nerror, TRUE);
  CTYPE ("Groups treated Quantum Mechanically");
  STYPE ("QMMM-grps",  QMMM,          NULL);
  CTYPE ("QM method");
  STYPE("QMmethod",     QMmethod, NULL);
  CTYPE ("QMMM scheme");
  EETYPE("QMMMscheme",  ir->QMMMscheme,    eQMMMscheme_names, nerror, TRUE);
  CTYPE ("QM basisset");
  STYPE("QMbasis",      QMbasis, NULL);
  CTYPE ("QM charge");
  STYPE ("QMcharge",    QMcharge,NULL);
  CTYPE ("QM multiplicity");
  STYPE ("QMmult",      QMmult,NULL);
  CTYPE ("Surface Hopping");
  STYPE ("SH",          bSH, NULL);
  CTYPE ("CAS space options");
  STYPE ("CASorbitals",      CASorbitals,   NULL);
  STYPE ("CASelectrons",     CASelectrons,  NULL);
  STYPE ("SAon", SAon, NULL);
  STYPE ("SAoff",SAoff,NULL);
  STYPE ("SAsteps",  SAsteps, NULL);
  CTYPE ("Scale factor for MM charges");
  RTYPE ("MMChargeScaleFactor", ir->scalefactor, 1.0);
  CTYPE ("Optimization of QM subsystem");
  STYPE ("bOPT",          bOPT, NULL);
  STYPE ("bTS",          bTS, NULL);

  /* Simulated annealing */
  CCTYPE("SIMULATED ANNEALING");
  CTYPE ("Type of annealing for each temperature group (no/single/periodic)");
  STYPE ("annealing",   anneal,      NULL);
  CTYPE ("Number of time points to use for specifying annealing in each group");
  STYPE ("annealing_npoints", anneal_npoints, NULL);
  CTYPE ("List of times at the annealing points for each group");
  STYPE ("annealing_time",       anneal_time,       NULL);
  CTYPE ("Temp. at each annealing point, for each group.");
  STYPE ("annealing_temp",  anneal_temp,  NULL);
  
  /* Startup run */
  CCTYPE ("GENERATE VELOCITIES FOR STARTUP RUN");
  EETYPE("gen-vel",     opts->bGenVel,  yesno_names, nerror, TRUE);
  RTYPE ("gen-temp",    opts->tempi,    300.0);
  ITYPE ("gen-seed",    opts->seed,     173529);
  
  /* Shake stuff */
  CCTYPE ("OPTIONS FOR BONDS");
  EETYPE("constraints",	opts->nshake,	constraints, nerror, TRUE);
  CTYPE ("Type of constraint algorithm");
  EETYPE("constraint-algorithm",  ir->eConstrAlg, eshake_names, nerror, TRUE);
  CTYPE ("Do not constrain the start configuration");
  EETYPE("unconstrained-start", ir->bUncStart, yesno_names, nerror, TRUE);
  CTYPE ("Use successive overrelaxation to reduce the number of shake iterations");
  EETYPE("Shake-SOR", ir->bShakeSOR, yesno_names, nerror, TRUE);
  CTYPE ("Relative tolerance of shake");
  RTYPE ("shake-tol", ir->shake_tol, 0.0001);
  CTYPE ("Highest order in the expansion of the constraint coupling matrix");
  ITYPE ("lincs-order", ir->nProjOrder, 4);
  CTYPE ("Number of iterations in the final step of LINCS. 1 is fine for");
  CTYPE ("normal simulations, but use 2 to conserve energy in NVE runs.");
  CTYPE ("For energy minimization with constraints it should be 4 to 8.");
  ITYPE ("lincs-iter", ir->nLincsIter, 1);
  CTYPE ("Lincs will write a warning to the stderr if in one step a bond"); 
  CTYPE ("rotates over more degrees than");
  RTYPE ("lincs-warnangle", ir->LincsWarnAngle, 30.0);
  CTYPE ("Convert harmonic bonds to morse potentials");
  EETYPE("morse",       opts->bMorse,yesno_names, nerror, TRUE);

  /* Energy group exclusions */
  CCTYPE ("ENERGY GROUP EXCLUSIONS");
  CTYPE ("Pairs of energy groups for which all non-bonded interactions are excluded");
  STYPE ("energygrp_excl", egpexcl,     NULL);
  
  /* Refinement */
  CCTYPE("NMR refinement stuff");
  CTYPE ("Distance restraints type: No, Simple or Ensemble");
  EETYPE("disre",       opts->eDisre,   edisre_names, nerror, TRUE);
  CTYPE ("Force weighting of pairs in one distance restraint: Conservative or Equal");
  EETYPE("disre-weighting", ir->eDisreWeighting, edisreweighting_names, nerror, TRUE);
  CTYPE ("Use sqrt of the time averaged times the instantaneous violation");
  EETYPE("disre-mixed", ir->bDisreMixed, yesno_names, nerror, TRUE);
  RTYPE ("disre-fc",	ir->dr_fc,	1000.0);
  RTYPE ("disre-tau",	ir->dr_tau,	0.0);
  CTYPE ("Output frequency for pair distances to energy file");
  ITYPE ("nstdisreout", ir->nstdisreout, 100);
  CTYPE ("Orientation restraints: No or Yes");
  EETYPE("orire",       opts->bOrire,   yesno_names, nerror, TRUE);
  CTYPE ("Orientation restraints force constant and tau for time averaging");
  RTYPE ("orire-fc",	ir->orires_fc,	0.0);
  RTYPE ("orire-tau",	ir->orires_tau,	0.0);
  STYPE ("orire-fitgrp",orirefitgrp,    NULL);
  CTYPE ("Output frequency for trace(SD) and S to energy file");
  ITYPE ("nstorireout", ir->nstorireout, 100);
  CTYPE ("Dihedral angle restraints: No, Simple or Ensemble");
  EETYPE("dihre",       opts->eDihre,   edisre_names, nerror, TRUE);
  RTYPE ("dihre-fc",	ir->dihre_fc,	1000.0);
  RTYPE ("dihre-tau",	ir->dihre_tau,	0.0);
  CTYPE ("Output frequency for dihedral values to energy file");
  ITYPE ("nstdihreout", ir->nstdihreout, 100);

  /* Free energy stuff */
  CCTYPE ("Free energy control stuff");
  EETYPE("free-energy",	ir->efep, efep_names, nerror, TRUE);
  RTYPE ("init-lambda",	ir->init_lambda,0.0);
  RTYPE ("delta-lambda",ir->delta_lambda,0.0);
  RTYPE ("sc-alpha",ir->sc_alpha,0.0);
  ITYPE ("sc-power",ir->sc_power,1);
  RTYPE ("sc-sigma",ir->sc_sigma,0.3);

  /* Non-equilibrium MD stuff */  
  CCTYPE("Non-equilibrium MD stuff");
  STYPE ("acc-grps",    accgrps,        NULL);
  STYPE ("accelerate",  acc,            NULL);
  STYPE ("freezegrps",  freeze,         NULL);
  STYPE ("freezedim",   frdim,          NULL);
  RTYPE ("cos-acceleration", ir->cos_accel, 0);
  STYPE ("deform",      deform,         NULL);

  /* Electric fields */
  CCTYPE("Electric fields");
  CTYPE ("Format is number of terms (int) and for all terms an amplitude (real)");
  CTYPE ("and a phase angle (real)");
  STYPE ("E-x",   	efield_x,	NULL);
  STYPE ("E-xt",	efield_xt,	NULL);
  STYPE ("E-y",   	efield_y,	NULL);
  STYPE ("E-yt",	efield_yt,	NULL);
  STYPE ("E-z",   	efield_z,	NULL);
  STYPE ("E-zt",	efield_zt,	NULL);
  
  /* User defined thingies */
  CCTYPE ("User defined thingies");
  STYPE ("user1-grps",  user1,          NULL);
  STYPE ("user2-grps",  user2,          NULL);
  ITYPE ("userint1",    ir->userint1,   0);
  ITYPE ("userint2",    ir->userint2,   0);
  ITYPE ("userint3",    ir->userint3,   0);
  ITYPE ("userint4",    ir->userint4,   0);
  RTYPE ("userreal1",   ir->userreal1,  0);
  RTYPE ("userreal2",   ir->userreal2,  0);
  RTYPE ("userreal3",   ir->userreal3,  0);
  RTYPE ("userreal4",   ir->userreal4,  0);
#undef CTYPE

  write_inpfile(mdparout,ninp,inp,FALSE);
  for (i=0; (i<ninp); i++) {
    sfree(inp[i].name);
    sfree(inp[i].value);
  }
  sfree(inp);

  /* Process options if necessary */
  for(m=0; m<2; m++) {
    for(i=0; i<2*DIM; i++)
      dumdub[m][i]=0.0;
    if(ir->epc) {
      switch (ir->epct) {
      case epctISOTROPIC:
	if (sscanf(dumstr[m],"%lf",&(dumdub[m][XX]))!=1) {
	  fprintf(stderr,
		  "ERROR: pressure coupling not enough values (I need 1)\n");
	  (*nerror)++;
	}
	dumdub[m][YY]=dumdub[m][ZZ]=dumdub[m][XX];
	break;
      case epctSEMIISOTROPIC:
      case epctSURFACETENSION:
	if (sscanf(dumstr[m],"%lf%lf",
		   &(dumdub[m][XX]),&(dumdub[m][ZZ]))!=2) {
	  fprintf(stderr,
		  "ERROR: pressure coupling not enough values (I need 2)\n");
	  (*nerror)++;
	}
	dumdub[m][YY]=dumdub[m][XX];
	break;
      case epctANISOTROPIC:
	if (sscanf(dumstr[m],"%lf%lf%lf%lf%lf%lf",
		   &(dumdub[m][XX]),&(dumdub[m][YY]),&(dumdub[m][ZZ]),
		   &(dumdub[m][3]),&(dumdub[m][4]),&(dumdub[m][5]))!=6) {
	  fprintf(stderr,
		  "ERROR: pressure coupling not enough values (I need 6)\n");
	  (*nerror)++;
	}
	break;
      default:
	gmx_fatal(FARGS,"Pressure coupling type %s not implemented yet",
		    epcoupltype_names[ir->epct]);
      }
    }
  }
  clear_mat(ir->ref_p);
  clear_mat(ir->compress);
  for(i=0; i<DIM; i++) {
    ir->ref_p[i][i]    = dumdub[1][i];
    ir->compress[i][i] = dumdub[0][i];
  }
  if (ir->epct == epctANISOTROPIC) {
    ir->ref_p[XX][YY] = dumdub[1][3];
    ir->ref_p[XX][ZZ] = dumdub[1][4];
    ir->ref_p[YY][ZZ] = dumdub[1][5];
    ir->compress[XX][YY] = dumdub[0][3];
    ir->compress[XX][ZZ] = dumdub[0][4];
    ir->compress[YY][ZZ] = dumdub[0][5];
    for(i=0; i<DIM; i++)
      for(m=0; m<i; m++) {
	ir->ref_p[i][m] = ir->ref_p[m][i];
	ir->compress[i][m] = ir->compress[m][i];
      }
  } 
  
  if (ir->comm_mode == ecmNO)
    ir->nstcomm = 0;
  
  if (opts->bOrire && str_nelem(orirefitgrp,MAXPTR,NULL)!=1) {
    fprintf(stderr,"ERROR: Need one orientation restraint fit group\n");
    (*nerror)++;
  }

  clear_mat(ir->deform);
  for(i=0; i<6; i++)
    dumdub[0][i] = 0;
  m = sscanf(deform,"%lf %lf %lf %lf %lf %lf",
	     &(dumdub[0][0]),&(dumdub[0][1]),&(dumdub[0][2]),
	     &(dumdub[0][3]),&(dumdub[0][4]),&(dumdub[0][5]));
  for(i=0; i<3; i++)
    ir->deform[i][i] = dumdub[0][i];
  ir->deform[YY][XX] = dumdub[0][3];
  ir->deform[ZZ][XX] = dumdub[0][4];
  ir->deform[ZZ][YY] = dumdub[0][5];
  if (ir->epc != epcNO) {
    for(i=0; i<3; i++)
      for(j=0; j<=i; j++)
	if (ir->deform[i][j]!=0 && ir->compress[i][j]!=0) {
	  fprintf(stderr,"ERROR: A box element has deform set and compressibility > 0\n");
	  (*nerror)++;
	}
    for(i=0; i<3; i++)
      for(j=0; j<i; j++)
	if (ir->deform[i][j]!=0) {
	  for(m=j; m<DIM; m++)
	    if (ir->compress[m][j]!=0) {
	      sprintf(warn_buf,"An off-diagonal box element has deform set while compressibility > 0 for the same component of another box vector, this might lead to spurious periodicity effects.");
	      warning(NULL);
	    }
	}
  }

  sfree(dumstr[0]);
  sfree(dumstr[1]);
}

static int search_QMstring(char *s,int ng,const char *gn[])
{
  /* same as normal search_string, but this one searches QM strings */
  int i;

  for(i=0; (i<ng); i++)
    if (strcasecmp(s,gn[i]) == 0)
      return i;

  gmx_fatal(FARGS,"this QM method or basisset (%s) is not implemented\n!",s);

  return -1;

} /* search_QMstring */


static int search_string(char *s,int ng,char *gn[])
{
  int i;
  
  for(i=0; (i<ng); i++)
    if (strcasecmp(s,gn[i]) == 0)
      return i;
      
  gmx_fatal(FARGS,"Group %s not found in indexfile.\nMaybe you have non-default goups in your mdp file, while not using the '-n' option of grompp.\nIn that case use the '-n' option.\n",s);
  
  return -1;
}

static void do_numbering(t_atoms *atoms,int ng,char *ptrs[],
			 t_block *block,char *gnames[],
			 int gtype,int restnm,
			 bool bOneGroup,bool bVerbose)
{
  unsigned short *cbuf;
  t_grps *groups=&(atoms->grps[gtype]);
  int    i,j,gid,aj,ognr,ntot=0;
  const char *title;

  if (debug)
    fprintf(debug,"Starting numbering %d groups of type %d\n",ng,gtype);
  
  title = gtypes[gtype];
    
  snew(cbuf,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    cbuf[i]=NOGID;
  
  snew(groups->nm_ind,ng+1); /* +1 for possible rest group */
  for(i=0; (i<ng); i++) {
    /* Lookup the group name in the block structure */
    gid = search_string(ptrs[i],block->nr,gnames);
    if ((i==0) || !bOneGroup)
      groups->nm_ind[groups->nr++]=gid;
    if (debug) 
      fprintf(debug,"Found gid %d for group %s\n",gid,ptrs[i]);
    
    /* Now go over the atoms in the group */
    for(j=block->index[gid]; (j<block->index[gid+1]); j++) {
      aj=block->a[j];
      
      /* Range checking */
      if ((aj < 0) || (aj >= atoms->nr)) 
	gmx_fatal(FARGS,"Invalid atom number %d in indexfile",aj);
	
      /* Lookup up the old group number */
      ognr = cbuf[aj];
      if (ognr != NOGID) 
	gmx_fatal(FARGS,"Atom %d in multiple %s groups (%d and %d)",
		    aj+1,title,ognr+1,i+1);
      else {
	/* Store the group number in buffer */
	if (bOneGroup)
	  cbuf[aj] = 0;
	else
	  cbuf[aj] = i;
	ntot++;
      }
    }
  }
  
  /* Now check whether we have done all atoms */
  if (ntot != atoms->nr) {
    if (bVerbose)
      fprintf(stderr,"Making dummy/rest group for %s containing %d elements\n",
	      title,atoms->nr-ntot);
    /* Add group name "rest" */ 
    groups->nm_ind[groups->nr] = restnm;
    
    /* Assign the rest name to all atoms not currently assigned to a group */
    for(j=0; (j<atoms->nr); j++) {
      if (cbuf[j] == NOGID)
	cbuf[j] = groups->nr;
    }
    groups->nr++;
  }
  for(j=0; (j<atoms->nr); j++) 
    atoms->atom[j].grpnr[gtype]=cbuf[j];
  
  sfree(cbuf);
}

static void calc_nrdf(t_atoms *atoms,t_idef *idef,t_grpopts *opts,
		      char **gnames,int nstcomm,int comm_mode)
{
  int     ai,aj,i,j,d,g,imin,jmin;
  t_iatom *ia;
  int     *nrdf,*na_vcm,na_tot;
  double  *nrdf_vcm,nrdf_uc,n_sub;

  /* Calculate nrdf. 
   * First calc 3xnr-atoms for each group
   * then subtract half a degree of freedom for each constraint
   *
   * Only atoms and nuclei contribute to the degrees of freedom...
   */

  snew(nrdf_vcm,atoms->grps[egcVCM].nr);
  snew(na_vcm,atoms->grps[egcVCM].nr);
  
  for(i=0; i<atoms->grps[egcTC].nr; i++)
    opts->nrdf[i] = 0;
  for(i=0; i<atoms->grps[egcVCM].nr; i++)
    nrdf_vcm[i] = 0;

  snew(nrdf,atoms->nr);
  for(i=0; i<atoms->nr; i++) {
    nrdf[i] = 0;
    if ((atoms->atom[i].ptype == eptAtom) ||
	(atoms->atom[i].ptype == eptNucleus)) {
      g = atoms->atom[i].grpnr[egcFREEZE];
      /* Double count nrdf for particle i */
      for(d=0; d<DIM; d++)
	if (opts->nFreeze[g][d] == 0)
	  nrdf[i] += 2; 
      opts->nrdf[atoms->atom[i].grpnr[egcTC]] += 0.5*nrdf[i];
      nrdf_vcm[atoms->atom[i].grpnr[egcVCM]] += 0.5*nrdf[i];
    }
  }
  ia=idef->il[F_SHAKE].iatoms;
  for(i=0; (i<idef->il[F_SHAKE].nr); ) {
    /* Subtract degrees of freedom for the constraints,
     * if the particles still have degrees of freedom left.
     * If one of the particles is a vsite or a shell, then all
     * constraint motion will go there, but since they do not
     * contribute to the constraints the degrees of freedom do not
     * change.
     */
    ai=ia[1];
    aj=ia[2];
    if (((atoms->atom[ai].ptype == eptNucleus) ||
	 (atoms->atom[ai].ptype == eptAtom)) &&
	((atoms->atom[aj].ptype == eptNucleus) ||
	 (atoms->atom[aj].ptype == eptAtom))) {
      if (nrdf[ai] > 0) 
	jmin = 1;
      else
	jmin = 2;
      if (nrdf[aj] > 0)
	imin = 1;
      else
	imin = 2;
      imin = min(imin,nrdf[ai]);
      jmin = min(jmin,nrdf[aj]);
      nrdf[ai] -= imin;
      nrdf[aj] -= jmin;
      opts->nrdf[atoms->atom[ai].grpnr[egcTC]] -= 0.5*imin;
      opts->nrdf[atoms->atom[aj].grpnr[egcTC]] -= 0.5*jmin;
      nrdf_vcm[atoms->atom[ai].grpnr[egcVCM]] -= 0.5*imin;
      nrdf_vcm[atoms->atom[aj].grpnr[egcVCM]] -= 0.5*jmin;
    }
    ia += interaction_function[F_SHAKE].nratoms+1;
    i  += interaction_function[F_SHAKE].nratoms+1;
  }
  ia=idef->il[F_SETTLE].iatoms;
  for(i=0; i<idef->il[F_SETTLE].nr; ) {
    for(ai=ia[1]; ai<ia[1]+3; ai++) {
      opts->nrdf[atoms->atom[ai].grpnr[egcTC]] -= 1;
      nrdf_vcm[atoms->atom[ai].grpnr[egcVCM]] -= 1;
    }
    ia+=2;
    i+=2;
  }

  if (nstcomm != 0) {
    /* Subtract 3 from the number of degrees of freedom in each vcm group
     * when com translation is removed and 6 when rotation is removed
     * as well.
     */
    switch (comm_mode) {
    case ecmLINEAR:
      n_sub = 3;
      break;
    case ecmANGULAR:
      n_sub = 6;
      break;
    default:
      n_sub = 0;
      gmx_incons("Checking comm_mode");
    }
    
    for(i=0; i<atoms->grps[egcTC].nr; i++) {
      /* Count the number of atoms of TC group i for every VCM group */
      for(j=0; j<atoms->grps[egcVCM].nr; j++)
	na_vcm[j] = 0;
      na_tot = 0;
      for(ai=0; ai<atoms->nr; ai++)
	if (atoms->atom[ai].grpnr[egcTC] == i) {
	  na_vcm[atoms->atom[ai].grpnr[egcVCM]]++;
	  na_tot++;
	}
      /* Correct for VCM removal according to the fraction of each VCM
       * group present in this TC group.
       */
      nrdf_uc = opts->nrdf[i];
      opts->nrdf[i] = 0;
      for(j=0; j<atoms->grps[egcVCM].nr; j++)
	if (nrdf_vcm[j] > n_sub) {
	  opts->nrdf[i] += nrdf_uc*((double)na_vcm[j]/(double)na_tot)*
	    (nrdf_vcm[j] - n_sub)/nrdf_vcm[j];
	}
    }
  }
  for(i=0; (i<atoms->grps[egcTC].nr); i++) {
    if (opts->nrdf[i] < 0)
      opts->nrdf[i] = 0;
    fprintf(stderr,
	    "Number of degrees of freedom in T-Coupling group %s is %.2f\n",
	    gnames[atoms->grps[egcTC].nm_ind[i]],opts->nrdf[i]);
  }
  
  sfree(nrdf);
  sfree(nrdf_vcm);
  sfree(na_vcm);
}

static void decode_cos(char *s,t_cosines *cosine,bool bTime)
{
  char   *t;
  char   format[STRLEN],f1[STRLEN];
  double a,phi;
  int    i;
  
  t=strdup(s);
  trim(t);
  
  cosine->n=0;
  cosine->a=NULL;
  cosine->phi=NULL;
  if (strlen(t)) {
    sscanf(t,"%d",&(cosine->n));
    if (cosine->n <= 0) {
      cosine->n=0;
    } else {
      snew(cosine->a,cosine->n);
      snew(cosine->phi,cosine->n);
      
      sprintf(format,"%%*d");
      for(i=0; (i<cosine->n); i++) {
	strcpy(f1,format);
	strcat(f1,"%lf%lf");
	if (sscanf(t,f1,&a,&phi) < 2)
	  gmx_fatal(FARGS,"Invalid input for electric field shift: '%s'",t);
	cosine->a[i]=a;
	cosine->phi[i]=phi;
	strcat(format,"%*lf%*lf");
      }
    }
  }
  sfree(t);
}

static bool do_egp_flag(t_inputrec *ir,t_atoms *atoms,char **gnames,
			char *option,char *val,int flag)
{
  int  nelem,i,j,k,nr;
  char *names[MAXPTR];
  bool bSet;

  nelem = str_nelem(val,MAXPTR,names);
  if (nelem % 2 != 0)
    gmx_fatal(FARGS,"The number of groups for %s is odd",option);
  nr=atoms->grps[egcENER].nr;
  bSet = FALSE;
  for(i=0; i<nelem/2; i++) {
    j = 0;
    while ((j < nr) &&
	   strcasecmp(names[2*i],gnames[atoms->grps[egcENER].nm_ind[j]]))
      j++;
    if (j == nr)
      gmx_fatal(FARGS,"%s in %s is not an energy group\n",
		  names[2*i],option);
    k = 0;
    while ((k < nr) &&
	   strcasecmp(names[2*i+1],gnames[atoms->grps[egcENER].nm_ind[k]]))
      k++;
    if (k==nr)
      gmx_fatal(FARGS,"%s in %s is not an energy group\n",
	      names[2*i+1],option);
    if ((j < nr) && (k < nr)) {
      ir->opts.egp_flags[nr*j+k] |= flag;
      ir->opts.egp_flags[nr*k+j] |= flag;
      bSet = TRUE;
    }
  }

  return bSet;
}

void do_index(char *ndx,
	      t_symtab   *symtab,
	      t_atoms    *atoms,bool bVerbose,
	      t_inputrec *ir,t_idef *idef,rvec *v)
{
  t_block *grps;
  char    warnbuf[STRLEN],**gnames;
  int     nr,ntcg,ntau_t,nref_t,nacc,nofg,nSA,nSA_points,nSA_time,nSA_temp;
  int     nacg,nfreeze,nfrdim,nenergy,nuser;
  char    *ptr1[MAXPTR],*ptr2[MAXPTR],*ptr3[MAXPTR];
  int     i,j,k,restnm;
  real    SAtime;
  bool    bExcl,bTable,bSetTCpar,bAnneal;
  int     nQMmethod,nQMbasis,nQMcharge,nQMmult,nbSH,nCASorb,nCASelec,
    nSAon,nSAoff,nSAsteps,nQMg,nbOPT,nbTS;


  if (bVerbose)
    fprintf(stderr,"processing index file...\n");
  debug_gmx();
  if (ndx == NULL) {
    snew(grps,1);
    snew(grps->index,1);
    snew(gnames,1);
    analyse(atoms,grps,&gnames,FALSE,TRUE);
  } else
    grps = init_index(ndx,&gnames);
  
  snew(atoms->grpname,grps->nr+1);
  
  for(i=0; (i<grps->nr); i++)
    atoms->grpname[i]=put_symtab(symtab,gnames[i]);
  atoms->grpname[i]=put_symtab(symtab,"rest");
  restnm=i;
  srenew(gnames,grps->nr+1);
  gnames[restnm]=*(atoms->grpname[i]);
  atoms->ngrpname=grps->nr+1;
  
  for(i=0; (i<atoms->nr); i++)
    for(j=0; (j<egcNR); j++)
      atoms->atom[i].grpnr[j]=NOGID;

  ntau_t = str_nelem(tau_t,MAXPTR,ptr1);
  nref_t = str_nelem(ref_t,MAXPTR,ptr2);
  ntcg   = str_nelem(tcgrps,MAXPTR,ptr3);
  if ((ntau_t != ntcg) || (nref_t != ntcg)) 
    gmx_fatal(FARGS,"Invalid T coupling input: %d groups, %d ref_t values and "
		"%d tau_t values",ntcg,nref_t,ntau_t);  

  if (ir->eI != eiMD)
    ir->etc = etcNO;
  do_numbering(atoms,ntcg,ptr3,grps,gnames,egcTC,
	       restnm,FALSE,bVerbose);
  nr=atoms->grps[egcTC].nr;
  ir->opts.ngtc=nr;
  snew(ir->opts.nrdf,nr);
  snew(ir->opts.tau_t,nr);
  snew(ir->opts.ref_t,nr);
  bSetTCpar = ir->etc || ir->eI==eiSD || ir->eI==eiBD || ir->eI==eiTPI;
  if (ir->eI==eiBD && ir->bd_fric==0) {
    fprintf(stderr,"bd_fric=0, so tau_t will be used as the inverse friction constant(s)\n"); 
  }
  if (bSetTCpar) {
    if (nr != nref_t)
      gmx_fatal(FARGS,"Not enough ref_t and tau_t values!");
    for(i=0; (i<nr); i++) {
      ir->opts.tau_t[i]=atof(ptr1[i]);
      if (ir->opts.tau_t[i] < 0)
	gmx_fatal(FARGS,"tau_t for group %d negative",i);
    }
    for(i=0; (i<nr); i++) {
      ir->opts.ref_t[i]=atof(ptr2[i]);
      if (ir->opts.ref_t[i] < 0)
	gmx_fatal(FARGS,"ref_t for group %d negative",i);
    }
  }
  /* Simulated annealing for each group. There are nr groups */
  nSA = str_nelem(anneal,MAXPTR,ptr1);
  if (nSA == 1 && (ptr1[0][0]=='n' || ptr1[0][0]=='N'))
     nSA = 0;
  if(nSA>0 && nSA != nr) 
    gmx_fatal(FARGS,"Not enough annealing values: %d (for %d groups)\n",nSA,nr);
  else {
    snew(ir->opts.annealing,nr);
    snew(ir->opts.anneal_npoints,nr);
    snew(ir->opts.anneal_time,nr);
    snew(ir->opts.anneal_temp,nr);
    for(i=0;i<nr;i++) {
      ir->opts.annealing[i]=eannNO;
      ir->opts.anneal_npoints[i]=0;
      ir->opts.anneal_time[i]=NULL;
      ir->opts.anneal_temp[i]=NULL;
    }
    if (nSA > 0) {
      bAnneal=FALSE;
      for(i=0;i<nr;i++) { 
	if(ptr1[i][0]=='n' || ptr1[i][0]=='N') {
	  ir->opts.annealing[i]=eannNO;
	} else if(ptr1[i][0]=='s'|| ptr1[i][0]=='S') {
	  ir->opts.annealing[i]=eannSINGLE;
	  bAnneal=TRUE;
	} else if(ptr1[i][0]=='p'|| ptr1[i][0]=='P') {
	  ir->opts.annealing[i]=eannPERIODIC;
	  bAnneal=TRUE;
	} 
      } 
      if(bAnneal) {
	/* Read the other fields too */
	nSA_points = str_nelem(anneal_npoints,MAXPTR,ptr1);
	if(nSA_points!=nSA) 
	  gmx_fatal(FARGS,"Found %d annealing_npoints values for %d groups\n",nSA_points,nSA);
	for(k=0,i=0;i<nr;i++) {
	  ir->opts.anneal_npoints[i]=strtol(ptr1[i],NULL,10);
	  if(ir->opts.anneal_npoints[i]==1)
	    gmx_fatal(FARGS,"Please specify at least a start and an end point for annealing\n");
	  snew(ir->opts.anneal_time[i],ir->opts.anneal_npoints[i]);
	  snew(ir->opts.anneal_temp[i],ir->opts.anneal_npoints[i]);
	  k += ir->opts.anneal_npoints[i];
	}

	nSA_time = str_nelem(anneal_time,MAXPTR,ptr1);
	if(nSA_time!=k) 
	  gmx_fatal(FARGS,"Found %d annealing_time values, wanter %d\n",nSA_time,k);
	nSA_temp = str_nelem(anneal_temp,MAXPTR,ptr2);
	if(nSA_temp!=k) 
	  gmx_fatal(FARGS,"Found %d annealing_temp values, wanted %d\n",nSA_temp,k);

	for(i=0,k=0;i<nr;i++) {
	  
	  for(j=0;j<ir->opts.anneal_npoints[i];j++) {
	    ir->opts.anneal_time[i][j]=atof(ptr1[k]);
	    ir->opts.anneal_temp[i][j]=atof(ptr2[k]);
	    if(j==0) {
	      if(ir->opts.anneal_time[i][0] > (ir->init_t+GMX_REAL_EPS))
		gmx_fatal(FARGS,"First time point for annealing > init_t.\n");      
	    } else { 
	      /* j>0 */
	      if(ir->opts.anneal_time[i][j]<ir->opts.anneal_time[i][j-1])
		gmx_fatal(FARGS,"Annealing timepoints out of order: t=%f comes after t=%f\n",
			    ir->opts.anneal_time[i][j],ir->opts.anneal_time[i][j-1]);
	    }
	    if(ir->opts.anneal_temp[i][j]<0) 
	      gmx_fatal(FARGS,"Found negative temperature in annealing: %f\n",ir->opts.anneal_temp[i][j]);    
	    k++;
	  }
	}
	/* Print out some summary information, to make sure we got it right */
	for(i=0,k=0;i<nr;i++) {
	  if(ir->opts.annealing[i]!=eannNO) {
	    j=atoms->grps[egcTC].nm_ind[i];
	    fprintf(stderr,"Simulated annealing for group %s: %s, %d timepoints\n",
		    *(atoms->grpname[j]),eann_names[ir->opts.annealing[i]],
		    ir->opts.anneal_npoints[i]);
	    fprintf(stderr,"Time (ps)   Temperature (K)\n");
	    /* All terms except the last one */
	    for(j=0;j<(ir->opts.anneal_npoints[i]-1);j++) 
		fprintf(stderr,"%9.1f      %5.1f\n",ir->opts.anneal_time[i][j],ir->opts.anneal_temp[i][j]);
	    
	    /* Finally the last one */
	    j = ir->opts.anneal_npoints[i]-1;
	    if(ir->opts.annealing[i]==eannSINGLE)
	      fprintf(stderr,"%9.1f-     %5.1f\n",ir->opts.anneal_time[i][j],ir->opts.anneal_temp[i][j]);
	    else {
	      fprintf(stderr,"%9.1f      %5.1f\n",ir->opts.anneal_time[i][j],ir->opts.anneal_temp[i][j]);
	      if(fabs(ir->opts.anneal_temp[i][j]-ir->opts.anneal_temp[i][0])>GMX_REAL_EPS)
		fprintf(stderr,"Note: There is a temperature jump when your annealing loops back.\n");
	    }
	  }
	} 
      }
    }
  }	

  nacc = str_nelem(acc,MAXPTR,ptr1);
  nacg = str_nelem(accgrps,MAXPTR,ptr2);
  if (nacg*DIM != nacc)
    gmx_fatal(FARGS,"Invalid Acceleration input: %d groups and %d acc. values",
		nacg,nacc);
  do_numbering(atoms,nacg,ptr2,grps,gnames,egcACC,
	       restnm,FALSE,bVerbose);  nr=atoms->grps[egcACC].nr;
  snew(ir->opts.acc,nr);
  ir->opts.ngacc=nr;
  
  for(i=k=0; (i<nacg); i++)
    for(j=0; (j<DIM); j++,k++)
      ir->opts.acc[i][j]=atof(ptr1[k]);
  for( ;(i<nr); i++)
    for(j=0; (j<DIM); j++)
      ir->opts.acc[i][j]=0;
  
  nfrdim  = str_nelem(frdim,MAXPTR,ptr1);
  nfreeze = str_nelem(freeze,MAXPTR,ptr2);
  if (nfrdim != DIM*nfreeze)
    gmx_fatal(FARGS,"Invalid Freezing input: %d groups and %d freeze values",
		nfreeze,nfrdim);
  do_numbering(atoms,nfreeze,ptr2,grps,gnames,egcFREEZE,
	       restnm,FALSE,bVerbose);
  nr=atoms->grps[egcFREEZE].nr;
  ir->opts.ngfrz=nr;
  snew(ir->opts.nFreeze,nr);
  for(i=k=0; (i<nfreeze); i++)
    for(j=0; (j<DIM); j++,k++) {
      ir->opts.nFreeze[i][j]=(strncasecmp(ptr1[k],"Y",1)==0);
      if (!ir->opts.nFreeze[i][j]) {
	if (strncasecmp(ptr1[k],"N",1) != 0) {
	  sprintf(warnbuf,"Please use Y(ES) or N(O) for freezedim only "
		  "(not %s)", ptr1[k]);
	  warning(NULL);
	}
      }
    }
  for( ; (i<nr); i++)
    for(j=0; (j<DIM); j++)
      ir->opts.nFreeze[i][j]=0;
  
  nenergy=str_nelem(energy,MAXPTR,ptr1);
  do_numbering(atoms,nenergy,ptr1,grps,gnames,egcENER,
	       restnm,FALSE,bVerbose);
  ir->opts.ngener=atoms->grps[egcENER].nr;
  nuser=str_nelem(vcm,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcVCM,
	       restnm,FALSE,bVerbose);

  /* Now we have filled the freeze struct, so we can calculate NRDF */ 
  calc_nrdf(atoms,idef,&(ir->opts),gnames,ir->nstcomm,ir->comm_mode);
  if (v && NULL) {
    real fac,ntot=0;
    
    /* Must check per group! */
    for(i=0; (i<ir->opts.ngtc); i++) 
      ntot += ir->opts.nrdf[i];
    if (ntot != (DIM*atoms->nr)) {
      fac = sqrt(ntot/(DIM*atoms->nr));
      if (bVerbose)
	fprintf(stderr,"Scaling velocities by a factor of %.3f to account for constraints\n"
		"and removal of center of mass motion\n",fac);
      for(i=0; (i<atoms->nr); i++)
	svmul(fac,v[i],v[i]);
    }
  }
  
  nuser=str_nelem(user1,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcUser1,
	       restnm,FALSE,bVerbose);
  nuser=str_nelem(user2,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcUser2,
	       restnm,FALSE,bVerbose);
  nuser=str_nelem(xtc_grps,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcXTC,
	       restnm,TRUE,bVerbose);
  nofg = str_nelem(orirefitgrp,MAXPTR,ptr1);
  do_numbering(atoms,nofg,ptr1,grps,gnames,egcORFIT,
	       restnm,FALSE,bVerbose);

  /* QMMM input processing */
  nQMg          = str_nelem(QMMM,MAXPTR,ptr1);
  nQMmethod     = str_nelem(QMmethod,MAXPTR,ptr2);
  nQMbasis      = str_nelem(QMbasis,MAXPTR,ptr3);
  if((nQMmethod != nQMg)||(nQMbasis != nQMg)){
    gmx_fatal(FARGS,"Invalid QMMM input: %d groups %d basissets"
	      " and %d methods\n",nQMg,nQMbasis,nQMmethod);
  }
  /* group rest, if any, is always MM! */
  do_numbering(atoms,nQMg,ptr1,grps,gnames,egcQMMM,
               restnm,FALSE,bVerbose);
  nr = nQMg; /*atoms->grps[egcQMMM].nr;*/
  ir->opts.ngQM = nQMg;
  snew(ir->opts.QMmethod,nr);
  snew(ir->opts.QMbasis,nr);
  for(i=0;i<nr;i++){
    /* input consists of strings: RHF CASSCF PM3 .. These need to be
     * converted to the corresponding enum in names.c
     */
    ir->opts.QMmethod[i] = search_QMstring(ptr2[i],eQMmethodNR,
                                           eQMmethod_names);
    ir->opts.QMbasis[i]  = search_QMstring(ptr3[i],eQMbasisNR,
                                           eQMbasis_names);

  }
  nQMmult   = str_nelem(QMmult,MAXPTR,ptr1);
  nQMcharge = str_nelem(QMcharge,MAXPTR,ptr2);
  nbSH      = str_nelem(bSH,MAXPTR,ptr3);
  snew(ir->opts.QMmult,nr);
  snew(ir->opts.QMcharge,nr);
  snew(ir->opts.bSH,nr);

  for(i=0;i<nr;i++){
    ir->opts.QMmult[i]   = atoi(ptr1[i]);
    ir->opts.QMcharge[i] = atoi(ptr2[i]);
    ir->opts.bSH[i]      = (strncasecmp(ptr3[i],"Y",1)==0);
  }

  nCASelec  = str_nelem(CASelectrons,MAXPTR,ptr1);
  nCASorb   = str_nelem(CASorbitals,MAXPTR,ptr2);
  snew(ir->opts.CASelectrons,nr);
  snew(ir->opts.CASorbitals,nr);
  for(i=0;i<nr;i++){
    ir->opts.CASelectrons[i]= atoi(ptr1[i]);
    ir->opts.CASorbitals[i] = atoi(ptr2[i]);
  }
  /* special optimization options */

  nbOPT = str_nelem(bOPT,MAXPTR,ptr1);
  nbTS = str_nelem(bTS,MAXPTR,ptr2);
  snew(ir->opts.bOPT,nr);
  snew(ir->opts.bTS,nr);
  for(i=0;i<nr;i++){
    ir->opts.bOPT[i] = (strncasecmp(ptr1[i],"Y",1)==0);
    ir->opts.bTS[i]  = (strncasecmp(ptr2[i],"Y",1)==0);
  }
  nSAon     = str_nelem(SAon,MAXPTR,ptr1);
  nSAoff    = str_nelem(SAoff,MAXPTR,ptr2);
  nSAsteps  = str_nelem(SAsteps,MAXPTR,ptr3);
  snew(ir->opts.SAon,nr);
  snew(ir->opts.SAoff,nr);
  snew(ir->opts.SAsteps,nr);

  for(i=0;i<nr;i++){
    ir->opts.SAon[i]    = atof(ptr1[i]);
    ir->opts.SAoff[i]   = atof(ptr2[i]);
    ir->opts.SAsteps[i] = atoi(ptr3[i]);
  }
  /* end of QMMM input */

  if (bVerbose)
    for(i=0; (i<egcNR); i++) {
      fprintf(stderr,"%-16s has %d element(s):",gtypes[i],atoms->grps[i].nr); 
      for(j=0; (j<atoms->grps[i].nr); j++)
	fprintf(stderr," %s",*(atoms->grpname[atoms->grps[i].nm_ind[j]]));
      fprintf(stderr,"\n");
    }

  nr=atoms->grps[egcENER].nr;
  snew(ir->opts.egp_flags,nr*nr);

  bExcl = do_egp_flag(ir,atoms,gnames,"energygrp_excl",egpexcl,EGP_EXCL);
  if (bExcl && EEL_FULL(ir->coulombtype))
    warning("Can not exclude the lattice Coulomb energy between energy groups");

  bTable = do_egp_flag(ir,atoms,gnames,"energygrp_table",egptable,EGP_TABLE);
  if (bTable && !(ir->vdwtype == evdwUSER) && 
      !(ir->coulombtype == eelUSER) &&!(ir->coulombtype == eelPMEUSER))
    gmx_fatal(FARGS,"Can only have energy group pair tables in combination with user tables for VdW and/or Coulomb");

  decode_cos(efield_x,&(ir->ex[XX]),FALSE);
  decode_cos(efield_xt,&(ir->et[XX]),TRUE);
  decode_cos(efield_y,&(ir->ex[YY]),FALSE);
  decode_cos(efield_yt,&(ir->et[YY]),TRUE);
  decode_cos(efield_z,&(ir->ex[ZZ]),FALSE);
  decode_cos(efield_zt,&(ir->et[ZZ]),TRUE);
  
  for(i=0; (i<grps->nr); i++)
    sfree(gnames[i]);
  sfree(gnames);
  done_block(grps);
  sfree(grps);

}



static void check_disre(t_topology *sys)
{
  t_functype *functype;
  t_iparams  *ip;
  int i,ndouble,ftype;
  int label,old_label;
  
  if (sys->idef.il[F_DISRES].nr) {
    functype  = sys->idef.functype;
    ip        = sys->idef.iparams;
    ndouble   = 0;
    old_label = -1;
    for(i=0; i<sys->idef.ntypes; i++) {
      ftype = functype[i];
      if (ftype == F_DISRES) {
	label = ip[i].disres.label;
	if (label == old_label) {
	  fprintf(stderr,"Distance restraint index %d occurs twice\n",label);
	  ndouble++;
	}
	old_label = label;
      }
    }
    if (ndouble>0)
      gmx_fatal(FARGS,"Found %d double distance restraint indices,\n"
		  "probably the parameters for multiple pairs in one restraint "
		  "are not identical\n",ndouble);
  }
}

void triple_check(char *mdparin,t_inputrec *ir,t_topology *sys,int *nerror)
{
  char err_buf[256];
  int  i,m,npct;
  real *mgrp,mt;
  rvec acc;

  /* Generalized reaction field */  
  if (ir->opts.ngtc == 0) {
    sprintf(err_buf,"No temperature coupling while using coulombtype %s",
	    eel_names[eelGRF]);
    CHECK(ir->coulombtype == eelGRF);
  }
  else {
    sprintf(err_buf,"When using coulombtype = %s"
	    " ref_t for temperature coupling should be > 0",
	    eel_names[eelGRF]);
    CHECK((ir->coulombtype == eelGRF) && (ir->opts.ref_t[0] <= 0));
  }

  clear_rvec(acc);
  snew(mgrp,sys->atoms.grps[egcACC].nr);
  for(i=0; (i<sys->atoms.nr); i++) 
    mgrp[sys->atoms.atom[i].grpnr[egcACC]] += sys->atoms.atom[i].m;
  mt=0.0;
  for(i=0; (i<sys->atoms.grps[egcACC].nr); i++) {
    for(m=0; (m<DIM); m++)
      acc[m]+=ir->opts.acc[i][m]*mgrp[i];
    mt+=mgrp[i];
  }
  for(m=0; (m<DIM); m++) {
    if (fabs(acc[m]) > 1e-6) {
      char *dim[DIM] = { "X", "Y", "Z" };
      fprintf(stderr,
	      "Net Acceleration in %s direction, will %s be corrected\n",
	      dim[m],ir->nstcomm != 0 ? "" : "not");
      if (ir->nstcomm != 0) {
	acc[m]/=mt;
	for (i=0; (i<sys->atoms.grps[egcACC].nr); i++)
	  ir->opts.acc[i][m]-=acc[m];
      }
    }
  }
  sfree(mgrp);

  check_disre(sys);
}

void double_check(t_inputrec *ir,matrix box,t_molinfo *mol,int *nerror)
{
  real min_size,rlong;
  bool bTWIN;
  char *ptr;
  
  ptr = check_box(box);
  if (ptr) {
    fprintf(stderr,
	    "ERROR: %s\n",ptr);
    (*nerror)++;
  }  

  if( (ir->eConstrAlg==estSHAKE) && 
      (mol->plist[F_SHAKE].nr > 0) && 
      (ir->shake_tol <= 0.0) ) {
    fprintf(stderr,"ERROR: shake_tol must be > 0 instead of %g\n",
	    ir->shake_tol);
    (*nerror)++;
  }

  if( (ir->eConstrAlg==estLINCS) && mol->plist[F_SHAKE].nr>0) {
    /* If we have Lincs constraints: */
    if(ir->eI==eiMD && ir->etc==etcNO &&
       ir->eConstrAlg==estLINCS && ir->nLincsIter==1) {
      sprintf(warn_buf,"For energy conservation with LINCS, lincs_iter should be 2 or larger.\n"
	      "You can safely ignore this if your system doesn't have any LINCS-constrained bonds;\n"
	      "for water molecules we normally use the analytical SETTLE algorithm instead."); 
      warning(NULL);
    }
    
    if((ir->eI == eiSteep) && (ir->nLincsIter<4)) {
      sprintf(warn_buf,"For minimization with LINCS constraints, lincs_iter should be 4 to 8.");
      warning(NULL);
    }
  }

  if (ir->ePBC != epbcNONE) {
    if (ir->nstlist == 0) {
      warning("With nstlist=0 atoms are only put into the box at step 0, therefore drifting atoms might cause the simulation to crash.");
    }
    rlong = max(ir->rlist,max(ir->rcoulomb,ir->rvdw));
    bTWIN = (rlong > ir->rlist);
    if (ir->ns_type==ensGRID) {
      if (sqr(rlong) >= max_cutoff2(box)) {
	fprintf(stderr,"ERROR: The cut-off length is longer than half the shortest box vector or longer than the smallest box diagonal element. Increase the box size or decrease %s.\n",
		bTWIN ? (ir->rcoulomb==rlong ? "rcoulomb" : "rvdw"):"rlist");
	(*nerror)++;
      }
    } else {
      min_size = min(box[XX][XX],min(box[YY][YY],box[ZZ][ZZ]));
      if (2*rlong >= min_size) {
	fprintf(stderr,"ERROR: One of the box lengths is smaller than twice the cut-off length. Increase the box size or decrease rlist.");
	(*nerror)++;
	if (TRICLINIC(box))
	  fprintf(stderr,"Grid search might allow larger cut-off's than simple search with triclinic boxes.");
      }
    }
  }
}

