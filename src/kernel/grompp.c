/*
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 3.2.03
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

#include <sys/types.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <limits.h>

#include "sysstuff.h"
#include "smalloc.h"
#include "macros.h"
#include "string2.h"
#include "readir.h"
#include "toputil.h"
#include "topio.h"
#include "confio.h"
#include "copyrite.h"
#include "readir.h"
#include "symtab.h"
#include "names.h"
#include "grompp.h"
#include "random.h"
#include "vec.h"
#include "futil.h"
#include "statutil.h"
#include "splitter.h"
#include "sortwater.h"
#include "convparm.h"
#include "gmx_fatal.h"
#include "index.h"
#include "gmxfio.h"
#include "trnio.h"
#include "tpxio.h"
#include "vsite_parm.h"
#include "txtdump.h"
#include "calcgrid.h"
#include "add_par.h"
#include "enxio.h"
#include "perf_est.h"
#include "compute_io.h"
#include "gpp_atomtype.h"
#include "gpp_tomorse.h"
#include "mtop_util.h"
#include "genborn.h"

static int rm_interactions(int ifunc,int nrmols,t_molinfo mols[])
{
  int  i,n;
  
  n=0;
  /* For all the molecule types */
  for(i=0; i<nrmols; i++) {
    n += mols[i].plist[ifunc].nr;
    mols[i].plist[ifunc].nr=0;
  }
  return n;
}

static int check_atom_names(const char *fn1, const char *fn2, 
			    gmx_mtop_t *mtop, t_atoms *at)
{
  int mb,m,i,j,nmismatch;
  t_atoms *tat;
#define MAXMISMATCH 20

  if (mtop->natoms != at->nr)
    gmx_incons("comparing atom names");
  
  nmismatch=0;
  i = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    tat = &mtop->moltype[mtop->molblock[mb].type].atoms;
    for(m=0; m<mtop->molblock[mb].nmol; m++) {
      for(j=0; j < tat->nr; j++) {
	if (strcmp( *(tat->atomname[j]) , *(at->atomname[i]) ) != 0) {
	  if (nmismatch < MAXMISMATCH) {
	    fprintf(stderr,
		    "Warning: atom name %d in %s and %s does not match (%s - %s)\n",
		    i+1, fn1, fn2, *(tat->atomname[j]), *(at->atomname[i]));
	  } else if (nmismatch == MAXMISMATCH) {
	    fprintf(stderr,"(more than %d non-matching atom names)\n",MAXMISMATCH);
	  }
	  nmismatch++;
	}
	i++;
      }
    }
  }

  return nmismatch;
}

static void check_eg_vs_cg(gmx_mtop_t *mtop)
{
  int astart,mb,m,cg,j,firstj;
  unsigned char firsteg,eg;
  gmx_moltype_t *molt;
  
  /* Go through all the charge groups and make sure all their
   * atoms are in the same energy group.
   */
  
  astart = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    molt = &mtop->moltype[mtop->molblock[mb].type];
    for(m=0; m<mtop->molblock[mb].nmol; m++) {
      for(cg=0; cg<molt->cgs.nr;cg++) {
	/* Get the energy group of the first atom in this charge group */
	firstj = astart + molt->cgs.index[cg];
	firsteg = ggrpnr(&mtop->groups,egcENER,firstj);
	for(j=molt->cgs.index[cg]+1;j<molt->cgs.index[cg+1];j++) {
	  eg = ggrpnr(&mtop->groups,egcENER,astart+j);
	  if(eg != firsteg) {
	    gmx_fatal(FARGS,"atoms %d and %d in charge group %d of molecule type '%s' are in different energy groups",
		      firstj+1,astart+j+1,cg+1,*molt->name);
	  }
	}
      }
      astart += molt->atoms.nr;
    }
  }  
}

static void check_cg_sizes(const char *topfn,t_block *cgs)
{
  int maxsize,cg;

  maxsize = 0;
  for(cg=0; cg<cgs->nr; cg++)
    maxsize = max(maxsize,cgs->index[cg+1]-cgs->index[cg]);
 
  if (maxsize > 10) {
    set_warning_line(topfn,-1);
    sprintf(warn_buf,
	    "The largest charge group contains %d atoms.\n"
	    "Since atoms only see each other when the centers of geometry of the charge groups they belong to are within the cut-off distance, too large charge groups can lead to serious cut-off artifacts.\n"
	    "For efficiency and accuracy, charge group should consist of a few atoms.\n"
	    "For all-atom force fields use: CH3, CH2, CH, NH2, NH, OH, CO2, CO, etc.",
	    maxsize);
    warning_note(warn_buf);
  }
}

static void check_vel(gmx_mtop_t *mtop,rvec v[])
{
  gmx_mtop_atomloop_all_t aloop;
  t_atom *atom;
  int a;

  aloop = gmx_mtop_atomloop_all_init(mtop);
  while (gmx_mtop_atomloop_all_next(aloop,&a,&atom)) {
    if (atom->ptype == eptShell ||
	atom->ptype == eptBond  ||
	atom->ptype == eptVSite) {
      clear_rvec(v[a]);
    }
  }
}

static bool nint_ftype(gmx_mtop_t *mtop,t_molinfo *mi,int ftype)
{
  int nint,mb;

  nint = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    nint += mtop->molblock[mb].nmol*mi[mtop->molblock[mb].type].plist[ftype].nr;
  }

  return nint;
}

/* This routine reorders the molecule type array
 * in the order of use in the molblocks,
 * unused molecule types are deleted.
 */
static void renumber_moltypes(gmx_mtop_t *sys,
			      int *nmolinfo,t_molinfo **molinfo)
{
  int *order,norder,i;
  int mb,mi;
  t_molinfo *minew;

  snew(order,*nmolinfo);
  norder = 0;
  for(mb=0; mb<sys->nmolblock; mb++) {
    for(i=0; i<norder; i++) {
      if (order[i] == sys->molblock[mb].type) {
	break;
      }
    }
    if (i == norder) {
      /* This type did not occur yet, add it */
      order[norder] = sys->molblock[mb].type;
      /* Renumber the moltype in the topology */
      norder++;
    }
    sys->molblock[mb].type = i;
  }
  
  /* We still need to reorder the molinfo structs */
  snew(minew,norder);
  for(mi=0; mi<*nmolinfo; mi++) {
    for(i=0; i<norder; i++) {
      if (order[i] == mi) {
	break;
      }
    }
    if (i == norder) {
      done_mi(&(*molinfo)[mi]);
    } else {
      minew[i] = (*molinfo)[mi];
    }
  }
  sfree(*molinfo);

  *nmolinfo = norder;
  *molinfo  = minew;
}

static void molinfo2mtop(int nmi,t_molinfo *mi,gmx_mtop_t *mtop)
{
  int m;
  gmx_moltype_t *molt;

  mtop->nmoltype = nmi;
  snew(mtop->moltype,nmi);
  for(m=0; m<nmi; m++) {
    molt = &mtop->moltype[m];
    molt->name  = mi[m].name;
    molt->atoms = mi[m].atoms;
    /* ilists are copied later */
    molt->cgs   = mi[m].cgs;
    molt->excls = mi[m].excls;
  }
}

static void
new_status(const char *topfile,const char *topppfile,const char *confin,
	   t_gromppopts *opts,t_inputrec *ir,bool bZero,
	   bool bGenVel,bool bVerbose,t_state *state,
	   gpp_atomtype_t atype,gmx_mtop_t *sys,
	   int *nmi,t_molinfo **mi,t_params plist[],
	   int *comb,double *reppow,real *fudgeQQ,
	   bool bMorse,
	   int *nerror)
{
  t_molinfo   *molinfo=NULL;
  int         nmolblock;
  gmx_molblock_t *molblock,*molbs;
  t_atoms     *confat;
  int         mb,mbs,i,nrmols,nmismatch;
  char        buf[STRLEN];
  bool        bGB=FALSE;

  init_mtop(sys);

  /* Set boolean for GB */
  if(ir->implicit_solvent)
    bGB=TRUE;
  
  /* TOPOLOGY processing */
  sys->name = do_top(bVerbose,topfile,topppfile,opts,bZero,&(sys->symtab),
		     plist,comb,reppow,fudgeQQ,
		     atype,&nrmols,&molinfo,ir,
		     &nmolblock,&molblock,bGB);
  
  sys->nmolblock = 0;
  snew(sys->molblock,nmolblock);
  mbs;
  sys->natoms = 0;
  for(mb=0; mb<nmolblock; mb++) {
    if (sys->nmolblock > 0 &&
	molblock[mb].type == sys->molblock[sys->nmolblock-1].type) {
      /* Merge consecutive blocks with the same molecule type */
      sys->molblock[sys->nmolblock-1].nmol += molblock[mb].nmol;
      sys->natoms += molblock[mb].nmol*sys->molblock[sys->nmolblock-1].natoms_mol;
    } else if (molblock[mb].nmol > 0) {
      /* Add a new molblock to the topology */
      molbs = &sys->molblock[sys->nmolblock];
      *molbs = molblock[mb];
      molbs->natoms_mol = molinfo[molbs->type].atoms.nr;
      molbs->nposres_xA = 0;
      molbs->nposres_xB = 0;
      sys->natoms += molbs->nmol*molbs->natoms_mol;
      sys->nmolblock++;
    }
  }
  if (sys->nmolblock == 0) {
    gmx_fatal(FARGS,"No molecules were defined in the system");
  }

  renumber_moltypes(sys,&nrmols,&molinfo);

  if (bMorse)
    convert_harmonics(nrmols,molinfo,atype);

  if (ir->eDisre == edrNone) {
    i = rm_interactions(F_DISRES,nrmols,molinfo);
    if (i > 0) {
      set_warning_line("unknown",-1);
      sprintf(warn_buf,"disre = no, removed %d distance restraints",i);
      warning_note(NULL);
    }
  }
  if (opts->bOrire == FALSE) {
    i = rm_interactions(F_ORIRES,nrmols,molinfo);
    if (i > 0) {
      set_warning_line("unknown",-1);
      sprintf(warn_buf,"orire = no, removed %d orientation restraints",i);
      warning_note(NULL);
    }
  }
  if (opts->bDihre == FALSE) {
    i = rm_interactions(F_DIHRES,nrmols,molinfo);
    if (i > 0) {
      set_warning_line("unknown",-1);
      sprintf(warn_buf,"dihre = no, removed %d dihedral restraints",i);
      warning_note(NULL);
    }
  }
  
  /* Copy structures from msys to sys */
  molinfo2mtop(nrmols,molinfo,sys);
  
  /* COORDINATE file processing */
  if (bVerbose) 
    fprintf(stderr,"processing coordinates...\n");

  get_stx_coordnum(confin,&state->natoms);
  if (state->natoms != sys->natoms)
    gmx_fatal(FARGS,"number of coordinates in coordinate file (%s, %d)\n"
		"             does not match topology (%s, %d)",
	      confin,state->natoms,topfile,sys->natoms);
  else {
    /* make space for coordinates and velocities */
    char title[STRLEN];
    snew(confat,1);
    init_t_atoms(confat,state->natoms,FALSE);
    init_state(state,state->natoms,0);
    read_stx_conf(confin,title,confat,state->x,state->v,NULL,state->box);
    /* This call fixes the box shape for runs with pressure scaling */
    set_box_rel(ir,state);

    nmismatch = check_atom_names(topfile, confin, sys, confat);
    free_t_atoms(confat,TRUE);
    sfree(confat);
    
    if (nmismatch) {
      sprintf(buf,"%d non-matching atom name%s\n"
	      "atom names from %s will be used\n"
	      "atom names from %s will be ignored\n",
	      nmismatch,(nmismatch == 1) ? "" : "s",topfile,confin);
      warning(buf);
    }    
    if (bVerbose) 
      fprintf(stderr,"double-checking input for internal consistency...\n");
    double_check(ir,state->box,nint_ftype(sys,molinfo,F_CONSTR),nerror);
  }

  if (bGenVel) {
    real *mass;
    gmx_mtop_atomloop_all_t aloop;
    t_atom *atom;

    snew(mass,state->natoms);
    aloop = gmx_mtop_atomloop_all_init(sys);
    while (gmx_mtop_atomloop_all_next(aloop,&i,&atom)) {
      mass[i] = atom->m;
    }

    if (opts->seed == -1) {
      opts->seed = make_seed();
      fprintf(stderr,"Setting gen_seed to %d\n",opts->seed);
    }
    maxwell_speed(opts->tempi,opts->seed,sys,state->v);

    stop_cm(stdout,state->natoms,mass,state->x,state->v);
    sfree(mass);
  }

  *nmi = nrmols;
  *mi  = molinfo;
}

static void cont_status(const char *slog,const char *ener,
			bool bNeedVel,bool bGenVel, real fr_time,
			t_inputrec *ir,t_state *state,
			gmx_mtop_t *sys,
                        const output_env_t oenv)
     /* If fr_time == -1 read the last frame available which is complete */
{
  t_trxframe  fr;
  int         fp;

  fprintf(stderr,
	  "Reading Coordinates%s and Box size from old trajectory\n",
	  (!bNeedVel || bGenVel) ? "" : ", Velocities");
  if (fr_time == -1)
    fprintf(stderr,"Will read whole trajectory\n");
  else
    fprintf(stderr,"Will read till time %g\n",fr_time);
  if (!bNeedVel || bGenVel) {
    if (bGenVel)
      fprintf(stderr,"Velocities generated: "
	      "ignoring velocities in input trajectory\n");
    read_first_frame(oenv,&fp,slog,&fr,TRX_NEED_X);
  } else
    read_first_frame(oenv,&fp,slog,&fr,TRX_NEED_X | TRX_NEED_V);
  
  state->natoms = fr.natoms;

  if (sys->natoms != state->natoms)
    gmx_fatal(FARGS,"Number of atoms in Topology "
		"is not the same as in Trajectory");

  /* Find the appropriate frame */
  while ((fr_time == -1 || fr.time < fr_time) && read_next_frame(oenv,fp,&fr));
  
  close_trj(fp);

  if (fr.not_ok & FRAME_NOT_OK)
    gmx_fatal(FARGS,"Can not start from an incomplete frame");

  state->x = fr.x;
  if (bNeedVel && !bGenVel)
    state->v = fr.v;
  copy_mat(fr.box,state->box);
  /* Set the relative box lengths for preserving the box shape.
   * Note that this call can lead to differences in the last bit
   * with respect to using tpbconv to create a tpx file.
   */
  set_box_rel(ir,state);

  fprintf(stderr,"Using frame at t = %g ps\n",fr.time);
  fprintf(stderr,"Starting time for run is %g ps\n",ir->init_t); 
  
  if ((ir->epc != epcNO  || ir->etc ==etcNOSEHOOVER) && ener) {
    get_enx_state(ener,fr.time,&sys->groups,ir,state);
    preserve_box_shape(ir,state->box_rel,state->boxv);
  }
}

static void read_posres(gmx_mtop_t *mtop,t_molinfo *molinfo,bool bTopB,
			char *fn,
			int rc_scaling, int ePBC, 
			rvec com)
{
  bool   bFirst = TRUE;
  rvec   *x,*v,*xp;
  dvec   sum;
  double totmass;
  t_atoms dumat;
  matrix box,invbox;
  int    natoms,npbcdim=0;
  char   title[STRLEN];
  int    a,i,ai,j,k,mb,nat_molb;
  gmx_molblock_t *molb;
  t_params *pr;
  t_atom *atom;

  get_stx_coordnum(fn,&natoms);
  if (natoms != mtop->natoms) {
    sprintf(warn_buf,"The number of atoms in %s (%d) does not match the number of atoms in the topology (%d). Will assume that the first %d atoms in the topology and %s match.",fn,natoms,mtop->natoms,min(mtop->natoms,natoms),fn);
    warning(NULL);
  }
  snew(x,natoms);
  snew(v,natoms);
  init_t_atoms(&dumat,natoms,FALSE);
  read_stx_conf(fn,title,&dumat,x,v,NULL,box);
  
  npbcdim = ePBC2npbcdim(ePBC);
  clear_rvec(com);
  if (rc_scaling != erscNO) {
    copy_mat(box,invbox);
    for(j=npbcdim; j<DIM; j++) {
      clear_rvec(invbox[j]);
      invbox[j][j] = 1;
    }
    m_inv_ur0(invbox,invbox);
  }

  /* Copy the reference coordinates to mtop */
  clear_dvec(sum);
  totmass = 0;
  a = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    molb = &mtop->molblock[mb];
    nat_molb = molb->nmol*mtop->moltype[molb->type].atoms.nr;
    pr = &(molinfo[molb->type].plist[F_POSRES]);
    if (pr->nr > 0) {
      atom = mtop->moltype[molb->type].atoms.atom;
      for(i=0; (i<pr->nr); i++) {
	ai=pr->param[i].AI;
	if (ai >= natoms) {
	  gmx_fatal(FARGS,"Position restraint atom index (%d) in moltype '%s' is larger than number of atoms in %s (%d).\n",
		    ai+1,*molinfo[molb->type].name,fn,natoms);
	}
	if (rc_scaling == erscCOM) {
	  /* Determine the center of mass of the posres reference coordinates */
	  for(j=0; j<npbcdim; j++) {
	    sum[j] += atom[ai].m*x[a+ai][j];
	  }
	  totmass  += atom[ai].m;
	}
      }
      if (!bTopB) {
	molb->nposres_xA = nat_molb;
	snew(molb->posres_xA,molb->nposres_xA);
	for(i=0; i<nat_molb; i++) {
	  copy_rvec(x[a+i],molb->posres_xA[i]);
	}
      } else {
	molb->nposres_xB = nat_molb;
	snew(molb->posres_xB,molb->nposres_xB);
	for(i=0; i<nat_molb; i++) {
	  copy_rvec(x[a+i],molb->posres_xB[i]);
	}
      }
    }
    a += nat_molb;
  }
  if (rc_scaling == erscCOM) {
    if (totmass == 0)
      gmx_fatal(FARGS,"The total mass of the position restraint atoms is 0");
    for(j=0; j<npbcdim; j++)
      com[j] = sum[j]/totmass;
    fprintf(stderr,"The center of mass of the position restraint coord's is %6.3f %6.3f %6.3f\n",com[XX],com[YY],com[ZZ]);
  }

  if (rc_scaling != erscNO) {
    for(mb=0; mb<mtop->nmolblock; mb++) {
      molb = &mtop->molblock[mb];
      nat_molb = molb->nmol*mtop->moltype[molb->type].atoms.nr;
      if (molb->nposres_xA > 0 || molb->nposres_xB > 0) {
	xp = (!bTopB ? molb->posres_xA : molb->posres_xB);
	for(i=0; i<nat_molb; i++) {
	  for(j=0; j<npbcdim; j++) {
	    if (rc_scaling == erscALL) {
	      /* Convert from Cartesian to crystal coordinates */
	      xp[i][j] *= invbox[j][j];
	      for(k=j+1; k<npbcdim; k++) {
		xp[i][j] += invbox[k][j]*xp[i][k];
	      }
	    } else if (rc_scaling == erscCOM) {
	      /* Subtract the center of mass */
	      xp[i][j] -= com[j];
	    }
	  }
	}
      }
    }

    if (rc_scaling == erscCOM) {
      /* Convert the COM from Cartesian to crystal coordinates */
      for(j=0; j<npbcdim; j++) {
	com[j] *= invbox[j][j];
	for(k=j+1; k<npbcdim; k++) {
	  com[j] += invbox[k][j]*com[k];
	}
      }
    }
  }
  
  free_t_atoms(&dumat,TRUE);
  sfree(x);
  sfree(v);
}

static void gen_posres(gmx_mtop_t *mtop,t_molinfo *mi,
		       char *fnA, char *fnB,
		       int rc_scaling, int ePBC,
		       rvec com, rvec comB)
{
  int i,j;

  read_posres  (mtop,mi,FALSE,fnA,rc_scaling,ePBC,com);
  if (strcmp(fnA,fnB) != 0) {
    read_posres(mtop,mi,TRUE ,fnB,rc_scaling,ePBC,comB);
  }
}

static void set_wall_atomtype(gpp_atomtype_t at,t_gromppopts *opts,
			      t_inputrec *ir)
{
  int i;

  if (ir->nwall > 0)
    fprintf(stderr,"Searching the wall atom type(s)\n");
  for(i=0; i<ir->nwall; i++)
    ir->wall_atomtype[i] = get_atomtype_type(opts->wall_atomtype[i],at);
}

static int nrdf_internal(t_atoms *atoms)
{
  int i,nmass,nrdf;

  nmass = 0;
  for(i=0; i<atoms->nr; i++) {
    /* Vsite ptype might not be set here yet, so also check the mass */
    if ((atoms->atom[i].ptype == eptAtom ||
	 atoms->atom[i].ptype == eptNucleus)
	&& atoms->atom[i].m > 0) {
      nmass++;
    }
  }
  switch (nmass) {
  case 0:  nrdf = 0; break;
  case 1:  nrdf = 0; break;
  case 2:  nrdf = 1; break;
  default: nrdf = nmass*3 - 6; break;
  }
  
  return nrdf;
}

void
spline1d( double        dx,
		 double *      y,
		 int           n,
		 double *      u,
		 double *      y2 )
{
    int i;
    double p,q;
	
    y2[0] = 0.0;
    u[0]  = 0.0;
	
    for(i=1;i<n-1;i++)
    {
		p = 0.5*y2[i-1]+2.0;
        y2[i] = -0.5/p;
        q = (y[i+1]-2.0*y[i]+y[i-1])/dx;
		u[i] = (3.0*q/dx-0.5*u[i-1])/p;
    }
	
    y2[n-1] = 0.0;
	
    for(i=n-2;i>=0;i--)
    {
        y2[i] = y2[i]*y2[i+1]+u[i];
    }
}


void
interpolate1d( double     xmin,
			  double     dx,
			  double *   ya,
			  double *   y2a,
			  double     x,
			  double *   y,
			  double *   y1)
{
    int ix;
    double a,b;
	
    ix = (x-xmin)/dx;
	
    a = (xmin+(ix+1)*dx-x)/dx;
    b = (x-xmin-ix*dx)/dx;
	
    *y  = a*ya[ix]+b*ya[ix+1]+((a*a*a-a)*y2a[ix]+(b*b*b-b)*y2a[ix+1])*(dx*dx)/6.0;
    *y1 = (ya[ix+1]-ya[ix])/dx-(3.0*a*a-1.0)/6.0*dx*y2a[ix]+(3.0*b*b-1.0)/6.0*dx*y2a[ix+1];
}


void
setup_cmap (int              grid_spacing,
			int              nc,
			real *           grid ,
			gmx_cmap_t *     cmap_grid)
{
 	double *tmp_u,*tmp_u2,*tmp_yy,*tmp_y1,*tmp_t2,*tmp_grid;
	
    int    i,j,k,ii,jj,kk,idx;
	int    offset;
    double dx,xmin,v,v1,v2,v12;
    double phi,psi;
	
	snew(tmp_u,2*grid_spacing);
	snew(tmp_u2,2*grid_spacing);
	snew(tmp_yy,2*grid_spacing);
	snew(tmp_y1,2*grid_spacing);
	snew(tmp_t2,2*grid_spacing*2*grid_spacing);
	snew(tmp_grid,2*grid_spacing*2*grid_spacing);
	
    dx = 360.0/grid_spacing;
    xmin = -180.0-dx*grid_spacing/2;
	
	for(kk=0;kk<nc;kk++)
	{
		/* Compute an offset depending on which cmap we are using                                 
		 * Offset will be the map number multiplied with the grid_spacing * grid_spacing * 2      
		 */
		offset = kk * grid_spacing * grid_spacing * 2;
		
		for(i=0;i<2*grid_spacing;i++)
		{
			ii=(i+grid_spacing-grid_spacing/2)%grid_spacing;
			
			for(j=0;j<2*grid_spacing;j++)
			{
				jj=(j+grid_spacing-grid_spacing/2)%grid_spacing;
				tmp_grid[i*grid_spacing*2+j] = grid[offset+ii*grid_spacing+jj];
			}
		}
		
		for(i=0;i<2*grid_spacing;i++)
		{
			spline1d(dx,&(tmp_grid[2*grid_spacing*i]),2*grid_spacing,tmp_u,&(tmp_t2[2*grid_spacing*i]));
		}
		
		for(i=grid_spacing/2;i<grid_spacing+grid_spacing/2;i++)
		{
			ii = i-grid_spacing/2;
			phi = ii*dx-180.0;
			
			for(j=grid_spacing/2;j<grid_spacing+grid_spacing/2;j++)
			{
				jj = j-grid_spacing/2;
				psi = jj*dx-180.0;
				
				for(k=0;k<2*grid_spacing;k++)
				{
					interpolate1d(xmin,dx,&(tmp_grid[2*grid_spacing*k]),
								  &(tmp_t2[2*grid_spacing*k]),psi,&tmp_yy[k],&tmp_y1[k]);
				}
				
				spline1d(dx,tmp_yy,2*grid_spacing,tmp_u,tmp_u2);
				interpolate1d(xmin,dx,tmp_yy,tmp_u2,phi,&v,&v1);
				spline1d(dx,tmp_y1,2*grid_spacing,tmp_u,tmp_u2);
				interpolate1d(xmin,dx,tmp_y1,tmp_u2,phi,&v2,&v12);
				
				idx = ii*grid_spacing+jj;
				cmap_grid->cmapdata[kk].cmap[idx*4] = grid[offset+ii*grid_spacing+jj];
				cmap_grid->cmapdata[kk].cmap[idx*4+1] = v1;
				cmap_grid->cmapdata[kk].cmap[idx*4+2] = v2;
				cmap_grid->cmapdata[kk].cmap[idx*4+3] = v12;
			}
		}
	}
}				
				
void init_cmap_grid(gmx_cmap_t *cmap_grid, int ngrid, int grid_spacing)
{
	int i,k,nelem;
	
	cmap_grid->ngrid        = ngrid;
	cmap_grid->grid_spacing = grid_spacing;
	nelem                   = cmap_grid->grid_spacing*cmap_grid->grid_spacing;
	
	snew(cmap_grid->cmapdata,ngrid);
	
	for(i=0;i<cmap_grid->ngrid;i++)
	{
		snew(cmap_grid->cmapdata[i].cmap,4*nelem);
	}
}


static int count_constraints(gmx_mtop_t *mtop,t_molinfo *mi)
{
  int count,count_mol,i,mb;
  gmx_molblock_t *molb;
  t_params *plist;
  char buf[STRLEN];

  count = 0;
  for(mb=0; mb<mtop->nmolblock; mb++) {
    count_mol = 0;
    molb  = &mtop->molblock[mb];
    plist = mi[molb->type].plist;
      
    for(i=0; i<F_NRE; i++) {
      if (i == F_SETTLE)
	count_mol += 3*plist[i].nr;
      else if (interaction_function[i].flags & IF_CONSTRAINT)
	count_mol += plist[i].nr;
    }
      
    if (count_mol > nrdf_internal(&mi[molb->type].atoms)) {
      sprintf(buf,
	      "Molecule type '%s' has %d constraints.\n"
	      "For stability and efficiency there should not be more constraints than internal number of degrees of freedom: %d.\n",
	      *mi[molb->type].name,count_mol,
	      nrdf_internal(&mi[molb->type].atoms));
      warning(buf);
    }
    count += molb->nmol*count_mol;
  }

  return count;
}

int main (int argc, char *argv[])
{
  static const char *desc[] = {
    "The gromacs preprocessor",
    "reads a molecular topology file, checks the validity of the",
    "file, expands the topology from a molecular description to an atomic",
    "description. The topology file contains information about",
    "molecule types and the number of molecules, the preprocessor",
    "copies each molecule as needed. ",
    "There is no limitation on the number of molecule types. ",
    "Bonds and bond-angles can be converted into constraints, separately",
    "for hydrogens and heavy atoms.",
    "Then a coordinate file is read and velocities can be generated",
    "from a Maxwellian distribution if requested.",
    "grompp also reads parameters for the mdrun ",
    "(eg. number of MD steps, time step, cut-off), and others such as",
    "NEMD parameters, which are corrected so that the net acceleration",
    "is zero.",
    "Eventually a binary file is produced that can serve as the sole input",
    "file for the MD program.[PAR]",
    
    "grompp uses the atom names from the topology file. The atom names",
    "in the coordinate file (option [TT]-c[tt]) are only read to generate",
    "warnings when they do not match the atom names in the topology.",
    "Note that the atom names are irrelevant for the simulation as",
    "only the atom types are used for generating interaction parameters.[PAR]",

    "grompp calls a preprocessor to resolve includes, macros ",
    "etcetera. By default we use the cpp in your path. To specify a "
    "different macro-preprocessor (e.g. m4) or alternative location",

    "you can put a line in your parameter file specifying the path",
    "to that program. Specifying [TT]-pp[tt] will get the pre-processed",
    "topology file written out.[PAR]",
    
    "If your system does not have a c-preprocessor, you can still",
    "use grompp, but you do not have access to the features ",
    "from the cpp. Command line options to the c-preprocessor can be given",
    "in the [TT].mdp[tt] file. See your local manual (man cpp).[PAR]",
    
    "When using position restraints a file with restraint coordinates",
    "can be supplied with [TT]-r[tt], otherwise restraining will be done",
    "with respect to the conformation from the [TT]-c[tt] option.",
    "For free energy calculation the the coordinates for the B topology",
    "can be supplied with [TT]-rb[tt], otherwise they will be equal to",
    "those of the A topology.[PAR]",
    
    "Starting coordinates can be read from trajectory with [TT]-t[tt].",
    "The last frame with coordinates and velocities will be read,",
    "unless the [TT]-time[tt] option is used.",
    "Note that these velocities will not be used when [TT]gen_vel = yes[tt]",
    "in your [TT].mdp[tt] file. An energy file can be supplied with",
    "[TT]-e[tt] to have exact restarts when using pressure and/or",
    "Nose-Hoover temperature coupling. For an exact restart do not forget",
    "to turn off velocity generation and turn on unconstrained starting",
    "when constraints are present in the system.",
    "If you want to continue a crashed run, it is",
    "easier to use [TT]tpbconv[tt].[PAR]",

    "Using the [TT]-morse[tt] option grompp can convert the harmonic bonds",
    "in your topology to morse potentials. This makes it possible to break",
    "bonds. For this option to work you need an extra file in your $GMXLIB",
    "with dissociation energy. Use the -debug option to get more information",
    "on the workings of this option (look for MORSE in the grompp.log file",
    "using less or something like that).[PAR]",
    
    "By default all bonded interactions which have constant energy due to",
    "virtual site constructions will be removed. If this constant energy is",
    "not zero, this will result in a shift in the total energy. All bonded",
    "interactions can be kept by turning off [TT]-rmvsbds[tt]. Additionally,",
    "all constraints for distances which will be constant anyway because",
    "of virtual site constructions will be removed. If any constraints remain",
    "which involve virtual sites, a fatal error will result.[PAR]"
    
    "To verify your run input file, please make notice of all warnings",
    "on the screen, and correct where necessary. Do also look at the contents",
    "of the [TT]mdout.mdp[tt] file, this contains comment lines, as well as",
    "the input that [TT]grompp[tt] has read. If in doubt you can start grompp",
    "with the [TT]-debug[tt] option which will give you more information",
    "in a file called grompp.log (along with real debug info). Finally, you",
    "can see the contents of the run input file with the [TT]gmxdump[tt]",
    "program."
  };
  t_gromppopts *opts;
  gmx_mtop_t   *sys;
  int          nmi;
  t_molinfo    *mi;
  gpp_atomtype_t atype;
  t_inputrec   *ir;
  int          natoms,nvsite,comb,mt;
  t_params     *plist;
  t_state      state;
  matrix       box;
  real         max_spacing,fudgeQQ;
  double       reppow;
  char         fn[STRLEN],fnB[STRLEN];
  const char   *mdparin;
  int          nerror,ntype;
  bool         bNeedVel,bGenVel;
  bool         have_radius,have_vol,have_surftens,have_gb_radius,have_S_hct;
  bool         have_atomnumber;
  int		   n12,n13,n14;
  t_params     *gb_plist = NULL;
  gmx_genborn_t *born = NULL;
  output_env_t oenv;

  t_filenm fnm[] = {
    { efMDP, NULL,  NULL,        ffOPTRD },
    { efMDP, "-po", "mdout",     ffWRITE },
    { efSTX, "-c",  NULL,        ffREAD  },
    { efSTX, "-r",  NULL,        ffOPTRD },
    { efSTX, "-rb", NULL,        ffOPTRD },
    { efNDX, NULL,  NULL,        ffOPTRD },
    { efTOP, NULL,  NULL,        ffREAD  },
    { efTOP, "-pp", "processed", ffOPTWR },
    { efTPX, "-o",  NULL,        ffWRITE },
    { efTRN, "-t",  NULL,        ffOPTRD },
    { efEDR, "-e",  NULL,        ffOPTRD }
  };
#define NFILE asize(fnm)

  /* Command line options */
  static bool bVerbose=TRUE,bRenum=TRUE;
  static bool bRmVSBds=TRUE,bZero=FALSE;
  static int  i,maxwarn=0;
  static real fr_time=-1;
  t_pargs pa[] = {
    { "-v",       FALSE, etBOOL, {&bVerbose},
      "Be loud and noisy" },
    { "-time",    FALSE, etREAL, {&fr_time},
      "Take frame at or first after this time." },
    { "-rmvsbds",FALSE, etBOOL, {&bRmVSBds},
      "Remove constant bonded interactions with virtual sites" },
    { "-maxwarn", FALSE, etINT,  {&maxwarn},
      "Number of allowed warnings during input processing" },
    { "-zero",    FALSE, etBOOL, {&bZero},
      "Set parameters for bonded interactions without defaults to zero instead of generating an error" },
    { "-renum",   FALSE, etBOOL, {&bRenum},
      "Renumber atomtypes and minimize number of atomtypes" }
  };
  
  CopyRight(stdout,argv[0]);
  
  /* Initiate some variables */
  nerror=0;
  snew(ir,1);
  snew(opts,1);
  init_ir(ir,opts);
  
  /* Parse the command line */
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
                    asize(desc),desc,0,NULL,&oenv);
  
  init_warning(maxwarn);
  
  /* PARAMETER file processing */
  mdparin = opt2fn("-f",NFILE,fnm);
  set_warning_line(mdparin,-1);    
  get_ir(mdparin,opt2fn("-po",NFILE,fnm),ir,opts,&nerror);
  
  if (bVerbose) 
    fprintf(stderr,"checking input for internal consistency...\n");
  check_ir(mdparin,ir,opts,&nerror);

  if (ir->ld_seed == -1) {
    ir->ld_seed = make_seed();
    fprintf(stderr,"Setting the LD random seed to %d\n",ir->ld_seed);
  }

  bNeedVel = EI_STATE_VELOCITY(ir->eI);
  bGenVel  = (bNeedVel && opts->bGenVel);

  snew(plist,F_NRE);
  init_plist(plist);
  snew(sys,1);
  atype = init_atomtype();
  if (debug)
    pr_symtab(debug,0,"Just opened",&sys->symtab);
    
  strcpy(fn,ftp2fn(efTOP,NFILE,fnm));
  if (!gmx_fexist(fn)) 
    gmx_fatal(FARGS,"%s does not exist",fn);
  new_status(fn,opt2fn_null("-pp",NFILE,fnm),opt2fn("-c",NFILE,fnm),
	     opts,ir,bZero,bGenVel,bVerbose,&state,
	     atype,sys,&nmi,&mi,plist,&comb,&reppow,&fudgeQQ,
	     opts->bMorse,
	     &nerror);
  
  if (debug)
    pr_symtab(debug,0,"After new_status",&sys->symtab);
  
  if (count_constraints(sys,mi) && (ir->eConstrAlg == econtSHAKE)) {
    if (ir->eI == eiCG || ir->eI == eiLBFGS) {
      fprintf(stderr,
	      "ERROR: Can not do %s with %s, use %s\n",
	      EI(ir->eI),econstr_names[econtSHAKE],econstr_names[econtLINCS]);
      nerror++;
    }
    if (ir->bPeriodicMols) {
      fprintf(stderr,
	      "ERROR: can not do periodic molecules with %s, use %s\n",
	      econstr_names[econtSHAKE],econstr_names[econtLINCS]);
      nerror++;
    }
  }

  /* If we are doing GBSA, check that we got the parameters we need                                                            
   * This checking is to see if there are GBSA paratmeters for all                                                             
   * atoms in the force field. To go around this for testing purposes                                                          
   * comment out the nerror++ counter temporarliy                                                                              
   */
  have_radius=have_vol=have_surftens=have_gb_radius=have_S_hct=TRUE;
  for(i=0;i<get_atomtype_ntypes(atype);i++) {
    have_radius=have_radius       && (get_atomtype_radius(i,atype) > 0);
    have_vol=have_vol             && (get_atomtype_vol(i,atype) > 0);
    have_surftens=have_surftens   && (get_atomtype_surftens(i,atype) > 0);
    have_gb_radius=have_gb_radius && (get_atomtype_gb_radius(i,atype) > 0);
    have_S_hct=have_S_hct         && (get_atomtype_S_hct(i,atype) > 0);
  }
  if(!have_radius && ir->implicit_solvent==eisGBSA) {
    fprintf(stderr,"Can't do GB electrostatics; the forcefield is missing values for\n"
	    "atomtype radii, or they might be zero\n.");
    /* nerror++; */
  }
  /*
  if(!have_surftens && ir->implicit_solvent!=eisNO) {
    fprintf(stderr,"Can't do implicit solvent; the forcefield is missing values\n"
	    " for atomtype surface tension\n.");
    nerror++;                                                                                                                
  }
  */
  
  /* If we are doing QM/MM, check that we got the atom numbers */
  have_atomnumber = TRUE;
  for (i=0; i<get_atomtype_ntypes(atype); i++) {
    have_atomnumber = have_atomnumber && (get_atomtype_atomnumber(i,atype) >= 0);
  }
  if (!have_atomnumber && ir->bQMMM)
  {
    fprintf(stderr,"\n"
            "It appears as if you are trying to run a QM/MM calculation, but the force\n"
            "field you are using does not contain atom numbers fields. This is an\n"
            "optional field (introduced in Gromacs 3.3) for general runs, but mandatory\n"
            "for QM/MM. The good news is that it is easy to add - put the atom number as\n"
            "an integer just before the mass column in ffXXXnb.itp.\n"
            "NB: United atoms have the same atom numbers as normal ones.\n\n"); 
    nerror++;
  }

  if (nerror) {
    print_warn_num(FALSE);
    
    gmx_fatal(FARGS,"There were %d error(s) processing your input",nerror);
  }
  if (opt2bSet("-r",NFILE,fnm))
    sprintf(fn,"%s",opt2fn("-r",NFILE,fnm));
  else
    sprintf(fn,"%s",opt2fn("-c",NFILE,fnm));
  if (opt2bSet("-rb",NFILE,fnm))
    sprintf(fnB,"%s",opt2fn("-rb",NFILE,fnm));
  else
    strcpy(fnB,fn);

  if (nint_ftype(sys,mi,F_POSRES) > 0) {
    if (bVerbose) {
      fprintf(stderr,"Reading position restraint coords from %s",fn);
      if (strcmp(fn,fnB) == 0) {
	fprintf(stderr,"\n");
      } else {
	fprintf(stderr," and %s\n",fnB);
	if (ir->efep != efepNO && ir->n_flambda > 0) {
	  fprintf(stderr,"ERROR: can not change the position restraint reference coordinates with lambda togther with foreign lambda calculation.\n");
	  nerror++;
	}
      }
    }
    gen_posres(sys,mi,fn,fnB,
	       ir->refcoord_scaling,ir->ePBC,
	       ir->posres_com,ir->posres_comB);
  }
		
  nvsite = 0;
  /* set parameters for virtual site construction (not for vsiten) */
  for(mt=0; mt<sys->nmoltype; mt++) {
    nvsite +=
      set_vsites(bVerbose, &sys->moltype[mt].atoms, atype, mi[mt].plist);
  }
  /* now throw away all obsolete bonds, angles and dihedrals: */
  /* note: constraints are ALWAYS removed */
  if (nvsite) {
    for(mt=0; mt<sys->nmoltype; mt++) {
      clean_vsite_bondeds(mi[mt].plist,sys->moltype[mt].atoms.nr,bRmVSBds);
    }
  }
  
	/* If we are using CMAP, setup the pre-interpolation grid */
	if(plist->ncmap>0)
	{
		init_cmap_grid(&sys->cmap_grid, plist->nc, plist->grid_spacing);
		setup_cmap(plist->grid_spacing, plist->nc, plist->cmap,&sys->cmap_grid);
	}
	
  set_wall_atomtype(atype,opts,ir);
  if (bRenum) {
    renum_atype(plist, sys, ir->wall_atomtype, atype, bVerbose);
    ntype = get_atomtype_ntypes(atype);
  }
  
	/* PELA: Copy the atomtype data to the topology atomtype list */
	copy_atomtype_atomtypes(atype,&(sys->atomtypes));

	if (debug)
    pr_symtab(debug,0,"After renum_atype",&sys->symtab);

  if (bVerbose) 
    fprintf(stderr,"converting bonded parameters...\n");
	
  ntype = get_atomtype_ntypes(atype);
  convert_params(ntype, plist, mi, comb, reppow, fudgeQQ, sys);
  	
	if(ir->implicit_solvent)
	{
		printf("Constructing Generalized Born topology...\n");

		/* Check for -normvsbds switch to grompp, necessary for gb together with vsites */
		if(bRmVSBds && nvsite)
		{
			fprintf(stderr, "ERROR: Must use -normvsbds switch to grompp when doing Generalized Born\n"
					"together with virtual sites\n");
			nerror++;
		}
		
		if (nerror)
		{
			print_warn_num(FALSE);
			gmx_fatal(FARGS,"There were %d error(s) processing your input",nerror);
		}
		
		generate_gb_topology(sys,mi);
	}
	
  if (debug)
    pr_symtab(debug,0,"After convert_params",&sys->symtab);

  /* set ptype to VSite for virtual sites */
  for(mt=0; mt<sys->nmoltype; mt++) {
    set_vsites_ptype(FALSE,&sys->moltype[mt]);
  }
  if (debug) {
    pr_symtab(debug,0,"After virtual sites",&sys->symtab);
  }
  /* Check velocity for virtual sites and shells */
  if (bGenVel) {
    check_vel(sys,state.v);
  }
    
  /* check masses */
  check_mol(sys);
  
  for(i=0; i<sys->nmoltype; i++) {
    check_cg_sizes(ftp2fn(efTOP,NFILE,fnm),&sys->moltype[i].cgs);
  }

  check_warning_error(FARGS);
	
  if (bVerbose) 
    fprintf(stderr,"initialising group options...\n");
  do_index(mdparin,ftp2fn_null(efNDX,NFILE,fnm),
	   sys,bVerbose,ir,
	   bGenVel ? state.v : NULL);
	
  /* Init the temperature coupling state */
  init_gtc_state(&state,ir->opts.ngtc);

  if (bVerbose)
    fprintf(stderr,"Checking consistency between energy and charge groups...\n");
  check_eg_vs_cg(sys);
  
  if (debug)
    pr_symtab(debug,0,"After index",&sys->symtab);
  triple_check(mdparin,ir,sys,&nerror);
  close_symtab(&sys->symtab);
  if (debug)
    pr_symtab(debug,0,"After close",&sys->symtab);

  /* make exclusions between QM atoms */
  if (ir->bQMMM) {
    generate_qmexcl(sys,ir);
  }

  if (ftp2bSet(efTRN,NFILE,fnm)) {
    if (bVerbose)
      fprintf(stderr,"getting data from old trajectory ...\n");
    cont_status(ftp2fn(efTRN,NFILE,fnm),ftp2fn_null(efEDR,NFILE,fnm),
		bNeedVel,bGenVel,fr_time,ir,&state,sys,oenv);
  }

  if (ir->ePBC==epbcXY && ir->nwall!=2)
    clear_rvec(state.box[ZZ]);
  
  if (EEL_FULL(ir->coulombtype)) {
    /* Calculate the optimal grid dimensions */
    copy_mat(state.box,box);
    if (ir->ePBC==epbcXY && ir->nwall==2)
      svmul(ir->wall_ewald_zfac,box[ZZ],box[ZZ]);
    max_spacing = calc_grid(stdout,box,opts->fourierspacing,
			    &(ir->nkx),&(ir->nky),&(ir->nkz),1);
    if ((ir->coulombtype == eelPPPM) && (max_spacing > 0.1)) {
      set_warning_line(mdparin,-1);
      sprintf(warn_buf,"Grid spacing larger then 0.1 while using PPPM.");
      warning_note(NULL);
    }
  }

  if (ir->ePull != epullNO)
    set_pull_init(ir,sys,state.x,state.box,oenv,opts->pull_start);

  /*  reset_multinr(sys); */
  
  if (EEL_PME(ir->coulombtype)) {
    float ratio = pme_load_estimate(sys,ir,state.box);
    fprintf(stderr,"Estimate for the relative computational load of the PME mesh part: %.2f\n",ratio);
    if (ratio > 0.5)
      warning_note("The optimal PME mesh load for parallel simulations is below 0.5\n"
		   "and for highly parallel simulations between 0.25 and 0.33,\n"
		   "for higher performance, increase the cut-off and the PME grid spacing");
  }

  {
    double cio = compute_io(ir,sys->natoms,&sys->groups,F_NRE,1);
    sprintf(warn_buf,"This run will generate roughly %.0f Mb of data",cio);
    if (cio > 2000) {
      set_warning_line(mdparin,-1);
      warning_note(NULL);
    } else {
      printf("%s\n",warn_buf);
    }
  }
	
  if (bVerbose) 
    fprintf(stderr,"writing run input file...\n");

  print_warn_num(TRUE);
  state.lambda = ir->init_lambda;
  write_tpx_state(ftp2fn(efTPX,NFILE,fnm),ir,&state,sys);
  
  thanx(stderr);
  
  return 0;
}
