/*
 * $Id$
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.0
 * 
 * Copyright (c) 1991-2001
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
 * University of Groningen, The Netherlands
 * 
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
 * Do check out http://www.gromacs.org , or mail us at gromacs@gromacs.org .
 * 
 * And Hey:
 * Gravel Rubs Often Many Awfully Cauterized Sores
 */
static char *SRCID_readir_c = "$Id$";
#include <ctype.h>
#include <stdlib.h>
#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "physics.h"
#include "names.h"
#include "fatal.h"
#include "macros.h"
#include "rdgroup.h"
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
  egexcl[STRLEN];
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
    ir->epc = epcNO;

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
    if (ir->nstcomm > 0)
      warning("Tumbling ice-cubes: We are not removing rotation around center of mass in a non-periodic system. You should set nstcomm = -1.");
    if (ir->nstcomm == 0)
      warning("Flying ice-cubes: We are not removing center of mass motion in a non-periodic system. You should set nstcomm = -1 (will also stop rotation).");
  }
  
  if ((EEL_LR(ir->coulombtype)) && (ir->efep!=efepNO)) {
    warning("You are using lattice sum electrostatics with free energy integration. "
	    "This might give wrong results, since the lattice contribution "
	    "to the free energy not calculated.");
  }
  
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
	  ir->compress[ZZ][ZZ] < 0 || trace(ir->compress) == 0);
    
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
  
  if(ir->etc==etcNOSEHOOVER && ir->epc==epcBERENDSEN) {
    sprintf(warn_buf,"Using Berendsen pressure coupling invalidates the "
	    "true ensemble for the Nose-Hoover thermostat");
    warning(NULL);
  }
  /* ELECTROSTATICS */
  /* More checks are in triple check (grompp.c) */
  sprintf(err_buf,"epsilon_r must be >= 0 instead of %g\n",ir->epsilon_r);
  CHECK(ir->epsilon_r < 0);
  
  if (ir->epsilon_r == 0) {
    sprintf(err_buf,"epsilon_r can only be %f (=infinity) with (generalized)"
	    " reaction field",ir->epsilon_r);
    CHECK((ir->coulombtype != eelRF) && (ir->coulombtype != eelGRF));
  }
  
  if ((ir->coulombtype == eelCUT) || 
      (ir->coulombtype == eelRF) || (ir->coulombtype == eelGRF)) {
    /* Cut-off electrostatics (may be tabulated) */
    sprintf(err_buf,"With coulombtype = %s rcoulomb must be >= rlist",
	    eel_names[ir->coulombtype]);
    CHECK(ir->rlist > ir->rcoulomb);

    if ((ir->coulombtype == eelRF) || (ir->coulombtype == eelGRF)) {
      /* reaction field (at the cut-off) */
      if (ir->epsilon_r == 1.0) {
	sprintf(warn_buf,"Using epsilon_r = 1.0 with %s does not make sense",
		eel_names[ir->coulombtype]);
	warning(NULL);
      }
    }
  } else if (ir->coulombtype != eelUSER) {
    /* Tabulated electrostatics (no cut-off) */
    sprintf(err_buf,"With coulombtype = %s rcoulomb must be <= rlist",
	    eel_names[ir->coulombtype]);
    CHECK(ir->rcoulomb > ir->rlist);
    if(ir->coulombtype != eelEWALD && ir->coulombtype != eelPME) {
      sprintf(err_buf,"With coulombtype = %s rcoulomb_switch must be < rcoulomb",
		eel_names[ir->coulombtype]);
	CHECK(ir->rcoulomb_switch >= ir->rcoulomb);
	if (ir->rcoulomb_switch > ir->rcoulomb-0.0999) { 
      sprintf(warn_buf,"rcoulomb should be 0.1 to 0.3 nm larger than rcoulomb_switch to account for diffusion and the size of charge groups"); 
      warning(NULL);
    }
    }
  }

  
  if (ir->vdwtype == evdwCUT) {
    sprintf(err_buf,"With vdwtype = %s rvdw must be >= rlist",
	    eel_names[ir->vdwtype]);
    CHECK(ir->rlist > ir->rvdw);
  } else if (ir->vdwtype != evdwUSER) {
    sprintf(err_buf,"With vdwtype = %s rvdw_switch must be < rvdw",
	    evdw_names[ir->vdwtype]);
    CHECK(ir->rvdw_switch >= ir->rvdw);
    if (ir->rvdw_switch > ir->rvdw-0.0999) { 
      sprintf(warn_buf,"rvdw should be 0.1 to 0.3 nm larger than rvdw_switch to account for diffusion and the size of charge groups"); 
      warning(NULL);
    }
  }
}

void get_ir(char *mdparin,char *mdparout,
	    t_inputrec *ir,t_gromppopts *opts,int *nerror)
{
  char      *dumstr[2];
  double    dumdub[2][6];
  char      epsbuf[STRLEN];
  t_inpfile *inp;
  char      *tmp;
  int       i,m,ninp;
  char      dummy[STRLEN];
  double    epsje;
  
  inp=read_inpfile(mdparin,&ninp);

  snew(dumstr[0],STRLEN);
  snew(dumstr[1],STRLEN);
  CCTYPE ("VARIOUS PREPROCESSING OPTIONS");
  STYPE ("title",	opts->title,	NULL);
  STYPE ("cpp",		opts->cpp,	"/lib/cpp");
  STYPE ("include",	opts->include,	NULL);
  STYPE ("define",	opts->define,	NULL);
    
  CCTYPE ("RUN CONTROL PARAMETERS");
  EETYPE("integrator",  ir->eI,         ei_names, nerror, TRUE);
  CTYPE ("start time and timestep in ps");
  RTYPE ("tinit",	ir->init_t,	0.0);
  RTYPE ("dt",		ir->delta_t,	0.001);
  ITYPE ("nsteps",      ir->nsteps,     1);
  CTYPE ("number of steps for center of mass motion removal");
  ITYPE ("nstcomm",	ir->nstcomm,	1);
  STYPE ("comm-grps",   vcm,            NULL);
  
  CCTYPE ("LANGEVIN DYNAMICS OPTIONS");
  CTYPE ("Temperature, friction coefficient (amu/ps) and random seed");
  RTYPE ("bd-temp",     ir->bd_temp,    300.0);
  RTYPE ("bd-fric",     ir->bd_fric,    0.0);
  ITYPE ("ld-seed",     ir->ld_seed,    1993);
  
  /* Em stuff */
  CCTYPE ("ENERGY MINIMIZATION OPTIONS");
  CTYPE ("Force tolerance and initial step-size");
  RTYPE ("emtol",       ir->em_tol,     100.0);
  RTYPE ("emstep",      ir->em_stepsize,0.01);
  CTYPE ("Max number of iterations in relax_shells");
  ITYPE ("niter",       ir->niter,      20);
  CTYPE ("Frequency of steepest descents steps when doing CG");
  ITYPE ("nstcgsteep",	ir->nstcgsteep,	1000);
  
  /* Output options */
  CCTYPE ("OUTPUT CONTROL OPTIONS");
  CTYPE ("Output frequency for coords (x), velocities (v) and forces (f)");
  ITYPE ("nstxout",	ir->nstxout,	100);
  ITYPE ("nstvout",	ir->nstvout,	100);
  ITYPE ("nstfout",	ir->nstfout,	0);
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
  CTYPE ("Periodic boundary conditions: xyz or no");
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
  CTYPE ("Dielectric constant (DC) for cut-off or DC of reaction field");
  STYPE ("epsilon-r",   epsbuf,         "1");
  CTYPE ("Method for doing Van der Waals");
  EETYPE("vdw-type",	ir->vdwtype,    evdw_names, nerror, TRUE);
  CTYPE ("cut-off lengths");
  RTYPE ("rvdw-switch",	ir->rvdw_switch,	0.0);
  RTYPE ("rvdw",	ir->rvdw,	1.0);
    
  CTYPE ("Apply long range dispersion corrections for Energy and Pressure");
  EETYPE("DispCorr",    ir->eDispCorr,  edispc_names, nerror, TRUE);
  CTYPE ("Spacing for the PME/PPPM FFT grid");
  RTYPE ("fourierspacing", opts->fourierspacing,0.12);
  CTYPE ("FFT grid size, when a value is 0 fourierspacing will be used");
  ITYPE ("fourier_nx",  ir->nkx,         0);
  ITYPE ("fourier_ny",  ir->nky,         0);
  ITYPE ("fourier_nz",  ir->nkz,         0);
  CTYPE ("EWALD/PME/PPPM parameters");
  ITYPE ("pme_order",   ir->pme_order,   4);
  RTYPE ("ewald_rtol",  ir->ewald_rtol, 0.00001);
  RTYPE ("epsilon_surface", ir->epsilon_surface, 0.0);
  EETYPE("optimize_fft",ir->bOptFFT,  yesno_names, nerror, FALSE);

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
  
  /* Simulated annealing */
  CCTYPE ("SIMULATED ANNEALING CONTROL");
  EETYPE("annealing",	ir->bSimAnn,    yesno_names, nerror, TRUE);
  CTYPE ("Time at which temperature should be zero (ps)");
  RTYPE ("zero-temp_time",ir->zero_temp_time,0.0);
  
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
  CTYPE ("Relative tolerance of shake");
  RTYPE ("shake-tol", ir->shake_tol, 0.0001);
  CTYPE ("Highest order in the expansion of the constraint coupling matrix");
  ITYPE ("lincs-order", ir->nProjOrder, 4);
  CTYPE ("Lincs will write a warning to the stderr if in one step a bond"); 
  CTYPE ("rotates over more degrees than");
  RTYPE ("lincs-warnangle", ir->LincsWarnAngle, 30.0);
  CTYPE ("Convert harmonic bonds to morse potentials");
  EETYPE("morse",       opts->bMorse,yesno_names, nerror, TRUE);
  
  /* Refinement */
  CCTYPE("NMR refinement stuff");
  CTYPE ("Distance restraints type: No, Simple or Ensemble");
  EETYPE("disre",       opts->eDisre,   edisre_names, nerror, TRUE);
  CTYPE ("Force weighting of pairs in one distance restraint: Equal or Conservative");
  EETYPE("disre-weighting", ir->eDisreWeighting, edisreweighting_names, nerror, TRUE);
  CTYPE ("Use sqrt of the time averaged times the instantaneous violation");
  EETYPE("disre-mixed", ir->bDisreMixed, yesno_names, nerror, TRUE);
  RTYPE ("disre-fc",	ir->dr_fc,	1000.0);
  RTYPE ("disre-tau",	ir->dr_tau,	0.0);
  CTYPE ("Output frequency for pair distances to energy file");
  ITYPE ("nstdisreout", ir->nstdisreout, 100);
  
  /* Free energy stuff */
  CCTYPE ("Free energy control stuff");
  EETYPE("free-energy",	ir->efep, efep_names, nerror, TRUE);
  RTYPE ("init-lambda",	ir->init_lambda,0.0);
  RTYPE ("delta-lambda",ir->delta_lambda,0.0);
  RTYPE ("sc-alpha",ir->sc_alpha,0.0);
  RTYPE ("sc-sigma",ir->sc_sigma,0.3);

  /* Non-equilibrium MD stuff */  
  CCTYPE("Non-equilibrium MD stuff");
  STYPE ("acc-grps",    accgrps,        NULL);
  STYPE ("accelerate",  acc,            NULL);
  STYPE ("freezegrps",  freeze,         NULL);
  STYPE ("freezedim",   frdim,          NULL);
  RTYPE ("cos-acceleration", ir->cos_accel, 0);
  STYPE ("energygrp_excl", egexcl,      NULL);
  
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

  write_inpfile(mdparout,ninp,inp);
  for (i=0; (i<ninp); i++) {
    sfree(inp[i].name);
    sfree(inp[i].value);
  }
  sfree(inp);

  /* Process options if necessary */
  /* Convert to uppercase and trim the epsilon_r string buffer */
  upstring(epsbuf);
  trim(epsbuf);
  if (strlen(epsbuf) == 0) 
    ir->epsilon_r = 1;
  else {
    if (strstr(epsbuf,"INF") != NULL)
      ir->epsilon_r = 0;
    else if (sscanf(epsbuf,"%lf",&epsje) == 1)
      ir->epsilon_r = epsje;
    else {
      sprintf(warn_buf,"Invalid value for epsilon_r: %s, setting to 1",epsbuf);
      warning(NULL);
      ir->epsilon_r = 1;
    }
  }
  for(m=0; m<2; m++) {
    for(i=0; i<2*DIM; i++)
      dumdub[m][i]=0.0;
    if(ir->epc) {
      switch (ir->epct) {
      case epctISOTROPIC:
	if (sscanf(dumstr[m],"%lf",&(dumdub[m][XX]))==1)
	  dumdub[m][YY]=dumdub[m][ZZ]=dumdub[m][XX];
	else
	  fprintf(stderr,"pressure coupling not enough values (I need 1)\n");
	break;
      case epctSEMIISOTROPIC:
      case epctSURFACETENSION:
	if (sscanf(dumstr[m],"%lf%lf",
		   &(dumdub[m][XX]),&(dumdub[m][ZZ]))==2) 
	  dumdub[m][YY]=dumdub[m][XX];
	else
	  fprintf(stderr,"pressure coupling not enough values (I need 2)\n");
	break;
      case epctANISOTROPIC:
	if (sscanf(dumstr[m],"%lf%lf%lf%lf%lf%lf",
		   &(dumdub[m][XX]),&(dumdub[m][YY]),&(dumdub[m][ZZ]),
		   &(dumdub[m][3]),&(dumdub[m][4]),&(dumdub[m][5]))!=6) 
	  fprintf(stderr,"pressure coupling not enough values (I need 6)\n");
	break;
      default:
	fatal_error(0,"Pressure coupling type %s not implemented yet",
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
  fprintf(stderr,"Warning: as of GMX v 2.0 unit of compressibility is truly 1/bar\n");

  sfree(dumstr[0]);
  sfree(dumstr[1]);
}

static int str_nelem(char *str,int maxptr,char *ptr[])
{
  int  np=0;
  char *copy;
  
  copy=strdup(str);
  ltrim(copy);
  while (*copy != '\0') {
    if (np >= maxptr)
      fatal_error(0,"Too many groups on line: '%s' (max is %d)",
		  str,maxptr);
    ptr[np++]=copy;
    while ((*copy != '\0') && !isspace(*copy))
      copy++;
    if (*copy != '\0') {
      *copy='\0';
      copy++;
    }
    ltrim(copy);
  }

  return np;
}

static int search_string(char *s,int ng,char *gn[])
{
  int i;
  
  for(i=0; (i<ng); i++)
    if (strcasecmp(s,gn[i]) == 0)
      return i;
      
  fatal_error(0,"Group %s not found in indexfile.\nMaybe you have non-default goups in your mdp file, while not using the '-n' option of grompp.\nIn that case use the '-n' option.\n",s);
  
  return -1;
}

static void do_numbering(t_atoms *atoms,int ng,char *ptrs[],
			 t_block *block,char *gnames[],
			 int gtype,char *title,int restnm,
			 int *forward,bool bOneGroup,bool bVerbose)
{
  unsigned short *cbuf;
  t_grps *groups=&(atoms->grps[gtype]);
  int i,j,gid,aj,ognr,ntot=0;

  if (debug)
    fprintf(debug,"Starting numbering %d groups of type %d\n",ng,gtype);
    
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
	fatal_error(0,"Invalid atom number %d in indexfile",aj);
	
      /* Lookup up the old group number */
      ognr = cbuf[aj];
      if (ognr != NOGID) 
	fatal_error(0,"Atom %d in multiple %s groups (%d and %d)",
		    aj,title,gid,ognr);
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
  if (forward != NULL) {
    for(j=0; (j<atoms->nr); j++) 
      atoms->atom[j].grpnr[gtype]=cbuf[forward[j]];
  }
  else {
    for(j=0; (j<atoms->nr); j++) 
      atoms->atom[j].grpnr[gtype]=cbuf[j];
  }
  sfree(cbuf);
}

static void calc_nrdf(t_atoms *atoms,t_idef *idef,t_grpopts *opts,
		      char **gnames,int nstcomm)
{
  int     ai,aj,i,j,d,g,imin,jmin;
  t_iatom *ia;
  int     *nrdf,*na_vcm,na_tot;
  double  *nrdf_vcm,nrdf_uc,n_sub;

  if (nstcomm<0 && atoms->grps[egcVCM].nr>1)
    fatal_error(0,"Can not remove rotation on more than one group");

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
     * If one of the particles is a dummy or a shell, then all
     * constraints motion will go there, but since they do not
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
    if (nstcomm > 0)
      n_sub = 3;
    else {
      if (atoms->nr == 2)
	n_sub = 5;
      else
	n_sub = 6;
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
	opts->nrdf[i] += nrdf_uc*((double)na_vcm[j]/(double)na_tot)*
	  (nrdf_vcm[j] - n_sub)/nrdf_vcm[j];
    }
  }
  for(i=0; i<atoms->grps[egcTC].nr; i++) {
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

static void decode_cos(char *s,t_cosines *cosine)
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
	if (sscanf(t,f1,&a,&phi) != 2)
	  fatal_error(0,"Invalid input for electric field shift: '%s'",t);
	cosine->a[i]=a;
	cosine->phi[i]=phi;
	strcat(format,"%*lf%*lf");
      }
    }
  }
  sfree(t);
}

void do_index(char *ndx,
	      t_symtab   *symtab,
	      t_atoms    *atoms,bool bVerbose,
	      t_inputrec *ir,t_idef *idef,int *forward)
{
  t_block *grps;
  char    warnbuf[STRLEN],**gnames;
  int     nr,ntcg,ntau_t,nref_t,nacc,nacg,nfreeze,nfrdim,nenergy,nuser,negexcl;
  char    *ptr1[MAXPTR],*ptr2[MAXPTR],*ptr3[MAXPTR];
  int     i,j,k,restnm;
  bool    bExcl,bSetTCpar;
  
  if (bVerbose)
    fprintf(stderr,"processing index file...\n");
  debug_gmx();
  if (ndx == NULL) {
    snew(grps,1);
    snew(grps->index,1);
    snew(gnames,1);
    analyse(atoms,grps,&gnames,FALSE,TRUE);
    /* Do not shuffle the index when it is based on atoms */
    forward = NULL;
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
  
  if (ir->bSimAnn) {
    if (ir->eI==eiMD && ir->etc==etcNO)
      fatal_error(0,"You must select a temperature coupling algorithm "
		  "for simulated annealing");
    if (ir->zero_temp_time == 0)
      fatal_error(0,"Cannot anneal to zero temp at t=0");
  }  

  ntau_t = str_nelem(tau_t,MAXPTR,ptr1);
  nref_t = str_nelem(ref_t,MAXPTR,ptr2);
  ntcg   = str_nelem(tcgrps,MAXPTR,ptr3);
  if ((ntau_t != ntcg) || (nref_t != ntcg)) 
    fatal_error(0,"Invalid T coupling input: %d groups, %d ref_t values and "
		"%d tau_t values",ntcg,nref_t,ntau_t);  

  if (ir->eI != eiMD)
    ir->etc = etcNO;
  do_numbering(atoms,ntcg,ptr3,grps,gnames,egcTC,"T-Coupling",
	       restnm,forward,FALSE,bVerbose);
  nr=atoms->grps[egcTC].nr;
  ir->opts.ngtc=nr;
  snew(ir->opts.nrdf,nr);
  snew(ir->opts.tau_t,nr);
  snew(ir->opts.ref_t,nr);
  bSetTCpar = ir->etc || ir->eI==eiSD;
  if (ir->eI==eiBD && ir->bd_fric==0) {
    fprintf(stderr,"bd_fric=0, so tau_t will be used as the inverse friction constant(s)\n"); 
    bSetTCpar = TRUE;
  }
  if (bSetTCpar) {
    if (nr != nref_t)
      fatal_error(0,"Not enough ref_t and tau_t values!");
    for(i=0; (i<nr); i++) {
      ir->opts.tau_t[i]=atof(ptr1[i]);
      if (ir->opts.tau_t[i] < 0)
	fatal_error(0,"tau_t for group %d negative",i);
    }
    for(i=0; (i<nr); i++) {
      ir->opts.ref_t[i]=atof(ptr2[i]);
      if (ir->opts.ref_t[i] < 0)
	fatal_error(0,"ref_t for group %d negative",i);
    }
  }
  
  nacc = str_nelem(acc,MAXPTR,ptr1);
  nacg = str_nelem(accgrps,MAXPTR,ptr2);
  if (nacg*DIM != nacc)
    fatal_error(0,"Invalid Acceleration input: %d groups and %d acc. values",
		nacg,nacc);
  do_numbering(atoms,nacg,ptr2,grps,gnames,egcACC,"Acceleration",
	       restnm,forward,FALSE,bVerbose);
  nr=atoms->grps[egcACC].nr;
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
    fatal_error(0,"Invalid Freezing input: %d groups and %d freeze values",
		nfreeze,nfrdim);
  do_numbering(atoms,nfreeze,ptr2,grps,gnames,egcFREEZE,"Freeze",
	       restnm,forward,FALSE,bVerbose);
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
  do_numbering(atoms,nenergy,ptr1,grps,gnames,egcENER,"Energy",
	       restnm,forward,FALSE,bVerbose);
  ir->opts.ngener=atoms->grps[egcENER].nr;
  nuser=str_nelem(vcm,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcVCM,"VCM",
	       restnm,forward,FALSE,bVerbose);

  /* Now we have filled the freeze struct, so we can calculate NRDF */ 
  calc_nrdf(atoms,idef,&(ir->opts),gnames,ir->nstcomm);
  
  nuser=str_nelem(user1,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcUser1,"User1",
	       restnm,forward,FALSE,bVerbose);
  nuser=str_nelem(user2,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcUser2,"User2",
	       restnm,forward,FALSE,bVerbose);
  nuser=str_nelem(xtc_grps,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcXTC,"xtc_grps",
	       restnm,forward,TRUE,bVerbose);
  if (bVerbose)
    for(i=0; (i<egcNR); i++) {
      fprintf(stderr,"%-16s has %d element(s):",gtypes[i],atoms->grps[i].nr); 
      for(j=0; (j<atoms->grps[i].nr); j++)
	fprintf(stderr," %s",*(atoms->grpname[atoms->grps[i].nm_ind[j]]));
      fprintf(stderr,"\n");
    }

  negexcl=str_nelem(egexcl,MAXPTR,ptr1);
  if (negexcl % 2 != 0)
    fatal_error(0,"The number of groups for energygrp_excl is odd");
  nr=atoms->grps[egcENER].nr;
  snew(ir->opts.eg_excl,nr*nr);
  bExcl=FALSE;
  for(i=0; i<negexcl/2; i++) {
    j=0;
    while ((j < nr) &&
	   strcasecmp(ptr1[2*i],gnames[atoms->grps[egcENER].nm_ind[j]]))
      j++;
    if (j==nr)
      fatal_error(0,"%s in energygrp_excl is not an energy group\n",
		  ptr1[2*i]);
    k=0;
    while ((k < nr) &&
	   strcasecmp(ptr1[2*i+1],gnames[atoms->grps[egcENER].nm_ind[k]]))
      k++;
    if (k==nr)
      fatal_error(0,"%s in energygrp_excl is not an energy group\n",
	      ptr1[2*i+1]);
    if ((j < nr) && (k < nr)) {
      ir->opts.eg_excl[nr*j+k] = TRUE;
      ir->opts.eg_excl[nr*k+j] = TRUE;
      bExcl = TRUE;
    }
  }
  if (bExcl && EEL_LR(ir->coulombtype))
    warning("Can not exclude the lattice Coulomb energy between energy groups");
  
  decode_cos(efield_x,&(ir->ex[XX]));
  decode_cos(efield_xt,&(ir->et[XX]));
  decode_cos(efield_y,&(ir->ex[YY]));
  decode_cos(efield_yt,&(ir->et[YY]));
  decode_cos(efield_z,&(ir->ex[ZZ]));
  decode_cos(efield_zt,&(ir->et[ZZ]));
  
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
  int index,old_index;
  
  if (sys->idef.il[F_DISRES].nr) {
    functype  = sys->idef.functype;
    ip        = sys->idef.iparams;
    ndouble   = 0;
    old_index = -1;
    for(i=0; i<sys->idef.ntypes; i++) {
      ftype = functype[i];
      if (ftype == F_DISRES) {
	index = ip[i].disres.index;
	if (index == old_index) {
	  fprintf(stderr,"Distance restraint index %d occurs twice\n",index);
	  ndouble++;
	}
	old_index = index;
      }
    }
    if (ndouble>0)
      fatal_error(0,"Found %d double distance restraint indices,\n"
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
      fprintf(stderr,"Net Acceleration in %s direction, will be corrected\n",
	      dim[m]);
      acc[m]/=mt;
      for (i=0; (i<sys->atoms.grps[egcACC].nr); i++)
	ir->opts.acc[i][m]-=acc[m];
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

  if (ir->ePBC != epbcNONE) {
    rlong = max(ir->rlist,max(ir->rcoulomb,ir->rvdw));
    bTWIN = (rlong > ir->rlist);
    if (ir->ns_type==ensGRID) {
      min_size = min(norm2(box[XX]),min(norm2(box[YY]),norm2(box[ZZ])));
      if (sqr(2*rlong) >= min_size) {
	fprintf(stderr,"ERROR: One of the box vectors is shorter than twice the cut-off length. Increase the box size or decrease %s.\n",
		bTWIN ? "rcoulomb":"rlist");
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

