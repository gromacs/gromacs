/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * GRoups of Organic Molecules in ACtion for Science
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

#define MAXPTR 254
#define NOGID  255

/* Resource parameters 
 * Do not change any of these until you read the instruction
 * in readinp.h. Some cpp's do not take spaces after the backslash
 * (like the c-shell), which will give you a very weird compiler
 * message.
 */


static char tcgrps[STRLEN],tau_t[STRLEN],ref_t[STRLEN],
  acc[STRLEN],accgrps[STRLEN],freeze[STRLEN],fdim[STRLEN],
  energy[STRLEN],user1[STRLEN],user2[STRLEN],user3[STRLEN],xtc_grps[STRLEN];
static char efield_x[STRLEN],efield_xt[STRLEN],efield_y[STRLEN],
  efield_yt[STRLEN],efield_z[STRLEN],efield_zt[STRLEN];

void init_ir(t_inputrec *ir, t_gromppopts *opts)
{
  snew(opts->title,STRLEN);
  snew(opts->cpp,STRLEN); 
  snew(opts->include,STRLEN); 
  snew(opts->define,STRLEN);
  snew(opts->watertype,STRLEN);
}

void check_ir(t_inputrec *ir, t_gromppopts *opts)
/* Check internal consistency */
{
  bool bStop=FALSE;

#define BS(b,s,val) if (b) fprintf(stderr,s,val), bStop=TRUE
  BS(((ir->tol <= 0.0) && (opts->nshake > 0)),
     "tol must be > 0 instead of %10.5e while using shake\n",ir->tol);
  BS(((ir->ndelta < 1) && (ir->ns_type==ensGRID)),
     "when you use ndelta=%d the neighbour search will be done\n"
     "on one grid cell, and take much longer than the simple method\n",
     ir->ndelta);
  BS((ir->epsilon_r <= 0),"Epsilon-R must be > 0 instead of %e\n",
     ir->epsilon_r);
  BS((ir->epc && (ir->tau_p <= 0)),
     "tau_p must be > 0 instead of %10.5e\n",ir->tau_p);
  if ((ir->rshort == 0.0) && (ir->rlong == 0.0)) {
    if (ir->eBox != ebtNONE) {
      fprintf(stderr,"Can not have cut-off to zero (=infinite) with periodic\n"
	      "boundary conditions. Either set the box type to %s, or increase\n"
	      "the cut-off radii\n",eboxtype_names[ebtNONE]);
      bStop = TRUE;
    }
  }
  if ((ir->epc == epcTRICLINIC) && (ir->eBox != ebtTRICLINIC)) {
    sprintf(warn_buf,
	    "%s pressure coupling does not make sense without %s box\n"
	    "resetting presure coupling to %s\n",
	    epcoupl_names[epcTRICLINIC],eboxtype_names[ebtTRICLINIC],
	    epcoupl_names[epcISOTROPIC]);
    warning(NULL);
    ir->epc = epcISOTROPIC;
  }
  if ((ir->eeltype == eelPPPM) && (ir->epc != epcNO)) {
    sprintf(warn_buf,"Can not (yet) use pressure coupling with PPPM: pressure"
	    " coupling turned off");
    warning(NULL);
    ir->epc = epcNO;
  }

  if ((ir->eeltype == eelSHIFT) || (ir->eeltype == eelSWITCH)) {
    if (ir->rshort >= ir->rlong) {
      fprintf(stderr,"rlong (%g) must be longer than rshort (%g) "
	      "when using %s\n",ir->rlong,ir->rshort,eel_names[ir->eeltype]);
      bStop=TRUE;
    }
  } 
  else if ((ir->eeltype == eelRF) || (ir->eeltype == eelGRF)) {
    if (ir->rlong != ir->rshort) {
      sprintf(warn_buf,
	      "specifying different rlong (%g) and rshort (%g) "
	      "when using %s is useless, since only rshort is used.\n",
	      ir->rlong,ir->rshort,eel_names[ir->eeltype]);
      warning(NULL);
      ir->rlong=ir->rshort;
    }
    if (ir->epsilon_r == 1.0) {
      sprintf(warn_buf,"Using epsilon_r = 1.0 with %s does not make sense",
	      eel_names[eelRF]);
      warning(NULL);
      ir->eeltype = eelTWIN; 
    }
  } else
    BS((ir->rshort > ir->rlong),
       "rshort (%g) must be <= rlong\n",ir->rshort);
       
  if (ir->delta_t > 0.005) {
    sprintf(warn_buf,"time step > 0.005! (%10.5e)\n",ir->delta_t);
    warning(NULL);
  }
  if (bStop) {
    fprintf(stderr,"program terminated\n");
    exit(1);
  }
}

void get_ir(char *mdparin,char *mdparout,
	    t_inputrec *ir,t_gromppopts *opts)
{
  static char *shake[eshNR+1]    = { 
    "none", "h-bonds", "all-bonds", "h-angles", "all-angles", NULL 
  };
  
  char      *dumstr[2];
  double    dumdub[2][DIM];
  t_inpfile *inp;
  char      *tmp;
  int       i,m,ninp;
  char      dummy[STRLEN];

  inp=read_inpfile(mdparin,&ninp);

  snew(dumstr[0],STRLEN);
  snew(dumstr[1],STRLEN);
  CCTYPE ("VARIOUS PREPROCESSING OPTIONS");
  STYPE ("title",	opts->title,	NULL);
  STYPE ("cpp",		opts->cpp,	"/lib/cpp");
  STYPE ("include",	opts->include,	NULL);
  STYPE ("define",	opts->define,	NULL);
    
  CCTYPE ("RUN CONTROL PARAMETERS");
  ETYPE ("integrator",  ir->eI,         ei_names);
  CTYPE ("start time and timestep in ps");
  RTYPE ("tinit",	ir->init_t,	0.0);
  RTYPE ("dt",		ir->delta_t,	0.001);
  CTYPE ("number of steps and, number of steps for center of mass motion removal");
  ITYPE ("nsteps",      ir->nsteps,     1);
  ITYPE ("nstcomm",	ir->nstcomm,	1);
  
  CCTYPE ("LANGEVIN DYNAMICS OPTIONS");
  CTYPE ("Temparature, friction coefficient (amu/ps) and random seed");
  RTYPE ("ld_temp",     ir->ld_temp,    300.0);
  RTYPE ("ld_fric",     ir->ld_fric,    0.0);
  ITYPE ("ld_seed",     ir->ld_seed,    1993);
  
  /* Em stuff */
  CCTYPE ("ENERGY MINIMIZATION OPTIONS");
  RTYPE ("emtol",       ir->em_tol,     0.001);
  RTYPE ("emstep",      ir->em_stepsize,0.1);
  
  /* Output options */
  CCTYPE ("OUTPUT CONTROL OPTIONS");
  CTYPE ("Output frequency for coords (x), velocities (v) and forces (f)");
  ITYPE ("nstxout",	ir->nstxout,	1);
  ITYPE ("nstvout",	ir->nstvout,	1);
  ITYPE ("nstfout",	ir->nstfout,	1);
  CTYPE ("Output frequency for group stuff and for energies (nstprint)");
  ITYPE ("nstgrp",	ir->nstgrp,	1);
  ITYPE ("nstprint",	ir->nstprint,	1);
  CTYPE ("Output frequency for xtc files, and associated precision");
  ITYPE ("nstxtcout",   ir->nstxtcout,  0);
  RTYPE ("xtc_precision",ir->xtcprec,   1000.0);
  CTYPE ("This selects the subset of atoms for the XTC file.");
  CTYPE ("Only the first group gets written out, it does not make sense to have");
  CTYPE ("multiple groups. By  default all atoms will be written");
  STYPE ("xtc_grps",    xtc_grps,       NULL);
  CTYPE ("Selection of energy groups");
  STYPE ("energygrps",  energy,         NULL);

  /* Neighbor searching */  
  CCTYPE ("NEIGHBORSEARCHING PARAMETERS");
  CTYPE ("nblist update frequency");
  ITYPE ("nstlist",	ir->nstlist,	10);
  CTYPE ("ns algorithm (simple or grid)");
  ETYPE ("ns_type",     ir->ns_type,    ens_names);
  ITYPE ("deltagrid",	ir->ndelta,	2);
  CTYPE ("Box type, rectangular, triclinic, none");
  ETYPE ("box",         ir->eBox,       eboxtype_names);
  
  /* Electrostatics */
  CCTYPE ("OPTIONS FOR ELECTROSTATICS");
  CTYPE ("Method for doing electrostatics");
  ETYPE ("eel_type",	ir->eeltype,    eel_names);
  CTYPE ("cut-off lengths");
  RTYPE ("rshort",	ir->rshort,	1.0);
  RTYPE ("rlong",	ir->rlong,	1.0);
  CTYPE ("Dielectric constant (DC) for twin-range or DC of reaction field");
  RTYPE ("epsilon_r",   ir->epsilon_r,  1.0);
  CTYPE ("Apply long range dispersion corrections for Energy and Pressure");
  ETYPE ("bLJcorr",     ir->bLJcorr,    yesno_names);
  CTYPE ("Some thingies for future use");
  ITYPE ("niter",       ir->niter,      100);
  RTYPE ("gauss_width", ir->gausswidth, 0.1);
  ITYPE ("fourier_nx",  ir->nkx,        10);
  ITYPE ("fourier_ny",  ir->nky,        10);
  ITYPE ("fourier_nz",  ir->nkz,        10);

  /* 1-4 Stuff */  
  CCTYPE ("LENNARD JONES 1-4 INTERACTION THINGIES");
  CTYPE ("Compute LJ 1-4 parameters by scaling normal LJ parameters by this");
  RTYPE ("fudgeLJ",     ir->fudgeLJ,    1.0);
  CTYPE ("Scale charge in 1-4 interaction by this");
  RTYPE ("fudgeQQ",     ir->fudgeQQ,    1.0);
  CTYPE ("Generate 1-4 parameters (sometimes for non GROMOS force fields)");
  ETYPE ("gen_pairs",   opts->bGenPairs,yesno_names);

  /* Coupling stuff */
  CCTYPE ("OPTIONS FOR WEAK COUPLING ALGORITHMS");
  CTYPE ("Temperature coupling");
  ETYPE ("tcoupl",	ir->btc,        yesno_names);
  CTYPE ("Memory for running average (steps)");
  ITYPE ("ntcmemory",   ir->ntcmemory,  1);
  CTYPE ("Groups to couple separately");
  STYPE ("tc_grps",     tcgrps,         NULL);
  CTYPE ("Time constant (ps) and reference temperature (K)");
  STYPE ("tau_t",	tau_t,		NULL);
  STYPE ("ref_t",	ref_t,		NULL);
  CTYPE ("Pressure coupling");
  ETYPE ("Pcoupl",	ir->epc,        epcoupl_names);
  CTYPE ("Memory for running average (steps)");
  ITYPE ("npcmemory",   ir->npcmemory,  1);
  CTYPE ("Time constant (ps), compressibility (1/bar) and reference P (bar)");
  RTYPE ("tau_p",	ir->tau_p,	1.0);
  STYPE ("compress",	dumstr[0],	NULL);
  STYPE ("ref_p",       dumstr[1],      NULL);
  
  /* Simulated annealing */
  CCTYPE ("SIMULATED ANNEALING CONTROL");
  ETYPE ("annealing",	ir->bSimAnn,    yesno_names);
  CTYPE ("Time at which temperature should be zero (ps)");
  RTYPE ("zero_temp_time",ir->zero_temp_time,0.0);
  
  /* Startup run */
  CCTYPE ("GENERATE VELOCITIES FOR STARTUP RUN");
  ETYPE ("gen_vel",     opts->bGenVel,  yesno_names);
  RTYPE ("gen_temp",    opts->tempi,    300.0);
  ITYPE ("gen_seed",    opts->seed,     173529);
  
  /* Optimization */
  CCTYPE ("OPTIMIZATIONS FOR WATER MODELS");
  CTYPE ("Name of water in the index file");
  STYPE ("watertype",   opts->watertype,NULL);
  CTYPE ("Number of atoms in water model. (Not implemented for non-three atom models)");
  ITYPE ("nwatoms",     ir->nwatoms,    3);

  /* Shake stuff */
  CCTYPE ("OPTIONS FOR CONSTRAINTS");
  ETYPE ("constraints",	opts->nshake,	shake);
  CTYPE ("Type of constraint solver");
  ETYPE ("shake_type",  ir->eShakeType, eshake_names);
  CTYPE ("Shake the start configuration (using either shake of lincs)");
  ETYPE ("shake_first", ir->bShakeFirst,yesno_names);
  CTYPE ("Relative tolerance of shake");
  RTYPE ("tol",		ir->tol,	0.0001);
  CTYPE ("Highest order in the expansion of the constraint coupling matrix");
  ITYPE ("lincs_order", ir->nProjOrder, 4);
  CTYPE ("Lincs will write a warning to the stderr if in one step a bond rotates"); 
  CTYPE ("over more degrees than");
  RTYPE ("lincs_warnangle", ir->LincsWarnAngle, 30.0);
  CTYPE ("Output frequency of the constraint accuracy");
  ITYPE ("nstLincsout",	ir->nstLincsout,100);
  
  /* Refinement */
  CCTYPE ("NMR refinement stuff");
  RTYPE ("dihre_fc",	ir->dihr_fc,	1000.0);
  RTYPE ("disre_fc",	ir->dr_fc,	1000.0);
  RTYPE ("disre_tau",	ir->dr_tau,	1.25);
  
  /* Free energy stuff */
  CCTYPE ("Free energy control stuff");
  ETYPE ("free_energy", ir->bPert,      yesno_names);
  RTYPE ("init_lambda",	ir->init_lambda,0.0);
  RTYPE ("delta_lambda",ir->delta_lambda,0.0);

  /* Non-equilibrium MD stuff */  
  CCTYPE ("Non-equilibrium MD stuff");
  STYPE ("acc_grps",    accgrps,        NULL);
  STYPE ("accelerate",  acc,            NULL);
  STYPE ("freezegrps",  freeze,         NULL);
  STYPE ("freezedim",   fdim,           NULL);
  
  /* Electric fields */
  CCTYPE ("Electric fields");
  STYPE ("E_x",   	efield_x,	NULL);
  STYPE ("E_xt",	efield_xt,	NULL);
  STYPE ("E_y",   	efield_y,	NULL);
  STYPE ("E_yt",	efield_yt,	NULL);
  STYPE ("E_z",   	efield_z,	NULL);
  STYPE ("E_zt",	efield_zt,	NULL);
  
  /* User defined thnigies */
  CCTYPE ("User defined thingies");
  STYPE ("user1_grps",  user1,          NULL);
  STYPE ("user2_grps",  user2,          NULL);
  STYPE ("user3_grps",  user3,          NULL);
  ITYPE ("userint1",    ir->userint1,   0);
  ITYPE ("userint2",    ir->userint2,   0);
  ITYPE ("userint3",    ir->userint3,   0);
  ITYPE ("userint4",    ir->userint4,   0);
  RTYPE ("userreal1",   ir->userreal1,  0);
  RTYPE ("userreal2",   ir->userreal2,  0);
  RTYPE ("userreal3",   ir->userreal3,  0);
  RTYPE ("userreal4",   ir->userreal4,  0);
  
  write_inpfile(mdparout,ninp,inp);
  for (i=0; (i<ninp); i++) {
    sfree(inp[i].name);
    sfree(inp[i].value);
  }
  sfree(inp);
#undef CTYPE
  for(m=0; (m<2); m++) {
    for(i=0; (i<DIM); i++)
      dumdub[m][i]=0.0;
    switch (ir->epc) {
    case epcNO:
      break;
    case epcISOTROPIC:
      if (sscanf(dumstr[m],"%lf",&(dumdub[m][XX]))==1)
	dumdub[m][YY]=dumdub[m][ZZ]=dumdub[m][XX];
      else
	fprintf(stderr,"Warning: pressure coupling not enough vals\n");
      break;
    case epcANISOTROPIC:
      if (sscanf(dumstr[m],"%lf%lf%lf",
		 &(dumdub[m][XX]),&(dumdub[m][YY]),&(dumdub[m][ZZ]))!=3) 
	fprintf(stderr,"Warning: pressure coupling not enough vals\n");
      break;
    default:
      fprintf(stderr,"Pressure coupling type %s not implemented yet\n",
	      epcoupl_names[ir->epc]);
      exit(1);
    }
  }
  for(i=0; (i<DIM); i++) {
    ir->ref_p[i]    = dumdub[1][i];
    ir->compress[i] = dumdub[0][i];
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
			 int *forward,bool bVerbose)
{
  ushort *cbuf;
  t_grps *groups=&(atoms->grps[gtype]);
  int i,j,gid,aj,ognr,ntot=0;
  
  snew(cbuf,atoms->nr);
  for(i=0; (i<atoms->nr); i++)
    cbuf[i]=NOGID;
  
  snew(groups->nm_ind,ng+1); /* +1 for possible rest group */
  for(i=0; (i<ng); i++) {
    gid=search_string(ptrs[i],block->nr,gnames);
    groups->nm_ind[groups->nr++]=gid;
    
    for(j=block->index[gid]; (j<block->index[gid+1]); j++) {
      aj=block->a[j];
      if ((aj < 0) || (aj >= atoms->nr)) 
	fatal_error(0,"Invalid atom number %d in indexfile",aj);
      ognr=cbuf[aj];/*atoms->atom[aj].grpnr[gtype];*/
      if (ognr != NOGID) 
	fatal_error(0,"Atom %d in multiple %s groups (%d and %d)",
		    aj,title,gid,ognr);
      else {
	cbuf[aj]=i;
	/*atoms->atom[aj].grpnr[gtype] = i;*/
	ntot++;
      }
    }
  }
  if (ntot != atoms->nr) {
    if (bVerbose)
      fprintf(stderr,"Making dummy/rest group for %s containing %d elements\n",
	      title,atoms->nr-ntot);
    i=groups->nr;
    groups->nm_ind[groups->nr++] = restnm;
    
    for(j=0; (j<atoms->nr); j++) {
      /*if (atoms->atom[j].grpnr[gtype] == NOGID)
	atoms->atom[j].grpnr[gtype]=i;*/
      if (cbuf[j] == NOGID)
	cbuf[j]=i;
    }
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

static void calc_nrdf(t_atoms *atoms,t_idef *idef,t_grpopts *opts)
{
  int     ai,aj,i,g,gi,gj;
  t_iatom *ia;
  real    *ndf;

  /* Calculate nrdf. 
   * First calc 3xnr-atoms for each group
   * then subtract half a degree of freedom for each constraint
   * This is done in integer by calcing 6xnr-atoms and subtracting
   * one degree of freedom and finally division by two.
   *
   * Only atoms and nuclei contribute to the degrees of freedom...
   */
  snew(ndf,atoms->grps[egcTC].nr);
  for(i=0; (i<atoms->nr); i++) {
    if ((atoms->atom[i].ptype == eptAtom) ||
	(atoms->atom[i].ptype == eptNucleus)) {
      g=atoms->atom[i].grpnr[egcTC];
      ndf[g]+=3;
    }
  }
  ia=idef->il[F_SHAKE].iatoms;
  for(i=0; (i<idef->il[F_SHAKE].nr); ) {
    ai=ia[1];
    aj=ia[2];
    gi=atoms->atom[ai].grpnr[egcTC];
    gj=atoms->atom[aj].grpnr[egcTC];
    ndf[gi]-=0.5;
    ndf[gj]-=0.5;
    ia+=3;
    i+=3;
  }
  ia=idef->il[F_SETTLE].iatoms;
  for(i=0; (i<idef->il[F_SETTLE].nr); ) {
    ai=ia[1];
    gi=atoms->atom[ai].grpnr[egcTC];
    ndf[gi]-=3;
    ia+=2;
    i+=2;
  }
  for(i=0; (i<atoms->grps[egcTC].nr); i++)
    opts->nrdf[i]=ndf[i];
    
  sfree(ndf);
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
  if (strlen(t) == 0) 
    return;
  sscanf(t,"%d",&(cosine->n));
  if (cosine->n <= 0) {
    cosine->n=0;
    return;
  }
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

void do_index(char *ndx,
	      t_symtab   *symtab,
	      t_atoms    *atoms,bool bVerbose,
	      t_inputrec *ir,t_idef *idef,int *forward)
{
  static char *gtypes[egcNR] = {
    "T-Coupling", "Energy Mon.", "Acceleration", "Freeze",
    "User1", "User2", "User3", "XTC"
  };
  t_block *grps;
  char    **gnames;
  int     nr,ntcg,ntau_t,nref_t,nacc,nacg,nfreeze,nfdim,nenergy,nuser;
  char    *ptr1[MAXPTR],*ptr2[MAXPTR],*ptr3[MAXPTR];
  int     i,j,k,restnm;
  
  if (bVerbose)
    fprintf(stderr,"processing index file...\n");
  if (ndx == NULL) {
    snew(grps,1);
    snew(grps->index,1);
    snew(gnames,1);
    analyse(atoms,grps,&gnames,FALSE,TRUE);
  }
  else {
    grps=init_index(ndx,&gnames);
  }
  snew(atoms->grpname,grps->nr+1);

  for(i=0; (i<grps->nr); i++)
    atoms->grpname[i]=put_symtab(symtab,gnames[i]);
  atoms->grpname[i]=put_symtab(symtab,"rest");
  restnm=i;
  atoms->ngrpname=grps->nr+1;

  for(i=0; (i<atoms->nr); i++)
    for(j=0; (j<egcNR); j++)
      atoms->atom[i].grpnr[j]=NOGID;
  
  if (ir->bSimAnn) {
    ir->btc=TRUE;
    if (ir->zero_temp_time == 0)
      fatal_error(0,"Cannot anneal to zero temp at t=0");
  }
  
  ntau_t=str_nelem(tau_t,MAXPTR,ptr1);
  nref_t=str_nelem(ref_t,MAXPTR,ptr2);
  ntcg=str_nelem(tcgrps,MAXPTR,ptr3);
  if ((ntau_t != ntcg) || (nref_t != ntcg)) 
    fatal_error(0,"Invalid T coupling input: %d groups, %d ref_t values and"
		"%d tau_t values",ntcg,nref_t,ntau_t);
 
  do_numbering(atoms,ntcg,ptr3,grps,gnames,egcTC,"T-Coupling",
	       restnm,forward,bVerbose);
  nr=atoms->grps[egcTC].nr;
  ir->opts.ngtc=nr;
  snew(ir->opts.nrdf,nr);
  snew(ir->opts.tau_t,nr);
  snew(ir->opts.ref_t,nr);
  if (ir->btc) {
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
  calc_nrdf(atoms,idef,&(ir->opts));
  
  nacc=str_nelem(acc,MAXPTR,ptr1);
  nacg=str_nelem(accgrps,MAXPTR,ptr2);
  if (nacg*DIM != nacc)
    fatal_error(0,"Invalid Acceleration input: %d groups and %d acc. values",
		nacg,nacc);
  do_numbering(atoms,nacg,ptr2,grps,gnames,egcACC,"Acceleration",
	       restnm,forward,bVerbose);
  nr=atoms->grps[egcACC].nr;
  snew(ir->opts.acc,nr);
  ir->opts.ngacc=nr;
  
  for(i=k=0; (i<nacg); i++)
    for(j=0; (j<DIM); j++,k++)
      ir->opts.acc[i][j]=atof(ptr1[k]);
  for( ;(i<nr); i++)
    for(j=0; (j<DIM); j++)
      ir->opts.acc[i][j]=0;
  
  nfdim=str_nelem(fdim,MAXPTR,ptr1);
  nfreeze=str_nelem(freeze,MAXPTR,ptr2);
  if (nfdim != DIM*nfreeze)
    fatal_error(0,"Invalid Freezing input: %d groups and %d freeze values",
		nfreeze,nfdim);
  do_numbering(atoms,nfreeze,ptr2,grps,gnames,egcFREEZE,"Freeze",
	       restnm,forward,bVerbose);
  nr=atoms->grps[egcFREEZE].nr;
  ir->opts.ngfrz=nr;
  snew(ir->opts.nFreeze,nr);
  for(i=k=0; (i<nfreeze); i++)
    for(j=0; (j<DIM); j++,k++) {
      ir->opts.nFreeze[i][j]=(strcasecmp(ptr1[k],"YES")==0);
      if (!ir->opts.nFreeze[i][j]) {
	if (strcasecmp(ptr1[k],"NO") != 0)
	  fprintf(stderr,"Please use YES or NO for freezedim only (not %s)\n",
		  ptr1[k]);
      }
    }
  for( ; (i<nr); i++)
    for(j=0; (j<DIM); j++)
      ir->opts.nFreeze[i][j]=0;
  
  nenergy=str_nelem(energy,MAXPTR,ptr1);
  do_numbering(atoms,nenergy,ptr1,grps,gnames,egcENER,"Energy",
	       restnm,forward,bVerbose);
  nr=atoms->grps[egcENER].nr;
  ir->opts.ngener=nr;
  
  nuser=str_nelem(user1,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcUser1,"User1",
	       restnm,forward,bVerbose);
  nuser=str_nelem(user2,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcUser2,"User2",
	       restnm,forward,bVerbose);
  nuser=str_nelem(user3,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcUser3,"User3",
	       restnm,forward,bVerbose);
  nuser=str_nelem(xtc_grps,MAXPTR,ptr1);
  do_numbering(atoms,nuser,ptr1,grps,gnames,egcXTC,"xtc_grps",
	       restnm,forward,bVerbose);
  if (bVerbose)
    for(i=0; (i<egcNR); i++) {
      fprintf(stderr,"%-16s has %d element(s):",gtypes[i],atoms->grps[i].nr); 
      for(j=0; (j<atoms->grps[i].nr); j++)
	fprintf(stderr," %s",*(atoms->grpname[atoms->grps[i].nm_ind[j]]));
      fprintf(stderr,"\n");
    }
  
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

