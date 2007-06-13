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
#include "topcat.h"
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
#include "compute_io.h"

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

static int check_atom_names(char *fn1, char *fn2, t_atoms *at1, t_atoms *at2)
{
  int i,nmismatch,idx;
#define MAXMISMATCH 20

  if (at1->nr != at2->nr)
    gmx_incons("comparing atom names");
  
  nmismatch=0;
  for(i=0; i < at1->nr; i++) {
    idx=i;
    if (strcmp( *(at1->atomname[i]) , *(at2->atomname[idx]) ) != 0) {
      if (nmismatch < MAXMISMATCH)
	fprintf(stderr,
		"Warning: atom names in %s and %s don't match (%s - %s)\n",
		fn1, fn2, *(at1->atomname[i]), *(at2->atomname[idx]));
      else if (nmismatch == MAXMISMATCH)
	fprintf(stderr,"(more than %d non-matching atom names)\n",MAXMISMATCH);
      nmismatch++;
    }
  }
  return nmismatch;
}

static void check_eg_vs_cg(t_atoms *atoms,t_block *cgblock)
{
  int i,j,firstj;
  unsigned char firsteg,eg;
  
  /* Go through all the charge groups and make sure all their
   * atoms are in the same energy group.
   */
  
  for(i=0;i<cgblock->nr;i++) {
    /* Get the energy group of the first atom in this charge group */
    firstj=cgblock->index[i];
    firsteg=atoms->atom[firstj].grpnr[egcENER];
    for(j=cgblock->index[i]+1;j<cgblock->index[i+1];j++) {
      eg=atoms->atom[j].grpnr[egcENER];
      if(eg!=firsteg) {
	gmx_fatal(FARGS,"atoms %d and %d in charge group %d are in different energy groups",
		    firstj+1,j+1,i+1);
      }
    }
  }  
}


static void check_pairs(int nrmols,t_molinfo mi[],int ntype,t_param *nb)
{
  int      i,j,jj,k,ai,aj,taiA,tajA,taiB,tajB,bLJ;
  t_params *p;
  
  for(i=0; (i<nrmols); i++) {
    p = &(mi[i].plist[F_LJ14]);
    for(j=jj=0; (j<p->nr); j++) {
      /* Extract atom types and hence the index in the nbf matrix */
      ai = p->param[j].a[0];
      aj = p->param[j].a[1];
      taiA = mi[i].atoms.atom[ai].type;
      tajA = mi[i].atoms.atom[aj].type;
      taiB = mi[i].atoms.atom[ai].typeB;
      tajB = mi[i].atoms.atom[aj].typeB;
      
      bLJ  = FALSE;
      for(k=0; (k<MAXFORCEPARAM); k++)
	bLJ = bLJ || (((nb[taiA].c[k]*nb[tajA].c[k]) != 0) || 
		      ((nb[taiB].c[k]*nb[tajB].c[k]) != 0));
      if (bLJ) {
	cp_param(&(p->param[jj]),&(p->param[j]));
	jj++;
      }
      else if (debug) {
	fprintf(debug,"Removed 1-4 interaction between atoms %d and %d (within mol %s)\n",
		ai+1,aj+1,*(mi[i].name));
      }
    }
    fprintf(stderr,"Removed %d 1-4 interactions for molecule %s\n",
	    p->nr-jj,*(mi[i].name));
    p->nr = jj;
  }
}

static void check_vel(t_atoms *atoms,rvec v[])
{
  int i;
  
  for(i=0; (i<atoms->nr); i++) {
    if ((atoms->atom[i].ptype == eptShell) ||
	(atoms->atom[i].ptype == eptBond)  ||
	(atoms->atom[i].ptype == eptVSite))
      clear_rvec(v[i]);
  }
}

static int num_real_atoms(t_topology *sys)
{
  int i, nr, n;
  t_atom* atoms;
  atoms = sys->atoms.atom;
  nr = sys->atoms.nr;
  n = 0;
  for (i = 0; i < nr; i++) {
    if (atoms[i].ptype == eptAtom || atoms[i].ptype == eptNucleus) {
      n++;
    }
  }
  return n;
}


static void
new_status(char *topfile,char *topppfile,char *confin,
	   t_gromppopts *opts,t_inputrec *ir,
	   bool bGenVel,bool bVerbose,t_state *state,
	   t_atomtype *atype,t_topology *sys,
	   t_molinfo *msys,t_params plist[],int *comb,real *reppow,
	   bool bEnsemble,bool bMorse,
	   bool bCheckPairs,int *nerror)
{
  t_molinfo   *molinfo=NULL;
  t_simsystem *Sims=NULL;
  t_atoms     *confat;
  int         i,nrmols,Nsim,nmismatch;
  char        buf[STRLEN];

  init_top(sys);
  init_molinfo(msys);
  
  /* TOPOLOGY processing */
  msys->name=do_top(bVerbose,topfile,topppfile,opts,&(sys->symtab),
		    plist,comb,reppow,atype,&nrmols,&molinfo,ir,&Nsim,&Sims);
  
  if (bCheckPairs)
    check_pairs(nrmols,molinfo,atype->nr,atype->nb);
  
  if (bMorse)
    convert_harmonics(nrmols,molinfo,atype);

  if (opts->eDisre == edrNone) {
    i = rm_interactions(F_DISRES,nrmols,molinfo);
    if (bVerbose && i)
      fprintf(stderr,"removed %d distance restraints\n",i);
  }
  if (opts->bOrire == FALSE) {
    i = rm_interactions(F_ORIRES,nrmols,molinfo);
    if (bVerbose && i)
      fprintf(stderr,"removed %d orientation restraints\n",i);
  }
  if (opts->eDihre == edrNone) {
    i = rm_interactions(F_DIHRES,nrmols,molinfo);
    if (bVerbose && i)
      fprintf(stderr,"removed %d dihedral restraints\n",i);
  }
  
  topcat(msys,nrmols,molinfo,Nsim,Sims,bEnsemble);
  
  /* Copy structures from msys to sys */
  mi2top(sys,msys);
  
  /* COORDINATE file processing */
  if (bVerbose) 
    fprintf(stderr,"processing coordinates...\n");
  
  get_stx_coordnum(confin,&state->natoms);
  if (state->natoms != sys->atoms.nr)
    gmx_fatal(FARGS,"number of coordinates in coordinate file (%s, %d)\n"
		"             does not match topology (%s, %d)",
		confin,state->natoms,topfile,sys->atoms.nr);
  else {
    /* make space for coordinates and velocities */
    snew(confat,1);
    init_t_atoms(confat,state->natoms,FALSE);
    init_state(state,state->natoms,0);
    read_stx_conf(confin,opts->title,confat,state->x,state->v,state->box);

    nmismatch=check_atom_names(topfile, confin, &(sys->atoms), confat);
    free_t_atoms(confat);
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
    double_check(ir,state->box,msys,nerror);
  }

  if (bGenVel) {
    real *mass;
    
    snew(mass,msys->atoms.nr);
    for(i=0; (i<msys->atoms.nr); i++)
      mass[i]=msys->atoms.atom[i].m;
    
    if (opts->seed == -1) {
      opts->seed = make_seed();
      fprintf(stderr,"Setting gen_seed to %d\n",opts->seed);
    }
    maxwell_speed(opts->tempi,num_real_atoms(sys)*DIM,
		  opts->seed,&(sys->atoms),state->v);
    stop_cm(stdout,sys->atoms.nr,mass,state->x,state->v);
    sfree(mass);
  }
  for(i=0; (i<nrmols); i++)
    done_mi(&(molinfo[i]));
  sfree(molinfo);
  sfree(Sims);

}

static void cont_status(char *slog,char *ener,
			bool bNeedVel,bool bGenVel, real fr_time,
			t_inputrec *ir,t_state *state,
			t_topology *sys)
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
    read_first_frame(&fp,slog,&fr,TRX_NEED_X);
  } else
    read_first_frame(&fp,slog,&fr,TRX_NEED_X | TRX_NEED_V);
  
  state->natoms = fr.natoms;

  if(sys->atoms.nr != state->natoms)
    gmx_fatal(FARGS,"Number of atoms in Topology "
		"is not the same as in Trajectory");

  /* Find the appropriate frame */
  while ((fr_time == -1 || fr.time < fr_time) && read_next_frame(fp,&fr));
  
  close_trj(fp);

  if (fr.not_ok & FRAME_NOT_OK)
    gmx_fatal(FARGS,"Can not start from an incomplete frame");

  state->x = fr.x;
  if (bNeedVel && !bGenVel)
    state->v = fr.v;
  copy_mat(fr.box,state->box);

  fprintf(stderr,"Using frame at t = %g ps\n",fr.time);
  fprintf(stderr,"Starting time for run is %g ps\n",ir->init_t); 
  
  if ((ir->epc != epcNO  || ir->etc ==etcNOSEHOOVER) && ener) {
    get_enx_state(ener,fr.time,&sys->atoms,ir,state);
  }
}

static void read_posres(t_params *pr, char *fn, int offset, 
			int rc_scaling, int ePBC, 
			int natoms_top, t_atom *atom, rvec com)
{
  bool   bFirst = TRUE;
  rvec   *x,*v;
  dvec   sum;
  double totmass;
  t_atoms dumat;
  matrix box,invbox;
  int    natoms,npbcdim=0;
  char   title[256];
  int    i,ai,j,k;
  
  get_stx_coordnum(fn,&natoms);
  if (natoms != natoms_top) {
    sprintf(warn_buf,"The number of atoms in %s (%d) does not match the number of atoms in the topology (%d). Will assume that the first %d atoms in the topology and %s match.",fn,natoms,natoms_top,min(natoms_top,natoms),fn);
    warning(NULL);
  }
  snew(x,natoms);
  snew(v,natoms);
  init_t_atoms(&dumat,natoms,FALSE);
  read_stx_conf(fn,title,&dumat,x,v,box);
  
  npbcdim = ePBC2npbcdim(ePBC);
  clear_rvec(com);
  if (rc_scaling != erscNO) {
    copy_mat(box,invbox);
    for(j=npbcdim; j<DIM; j++) {
      clear_rvec(invbox[j]);
      invbox[j][j] = 1;
    }
    m_inv_lowerleft0(invbox,invbox);
  }
  if (rc_scaling == erscCOM) {
    /* Determine the center of mass of the posres reference coordinates */
    clear_dvec(sum);
    totmass = 0;
    for(i=0; (i<pr->nr); i++) {
      ai=pr->param[i].AI;
      if (ai >= natoms)
	gmx_fatal(FARGS,"Position restraint atom index (%d) is larger than number of atoms in %s (%d).\n",ai+1,fn,natoms);
      for(j=0; j<npbcdim; j++)
	sum[j] += atom[ai].m*x[ai][j];
      totmass += atom[ai].m;
    }
    if (totmass == 0)
      gmx_fatal(FARGS,"The total mass of the position restraint atoms is 0");
    for(j=0; j<npbcdim; j++)
      com[j] = sum[j]/totmass;
    fprintf(stderr,"The center of mass of the position restraint coord's is %6.3f %6.3f %6.3f\n",com[XX],com[YY],com[ZZ]);
  }

  for(i=0; (i<pr->nr); i++) {
    ai=pr->param[i].AI;
    if (ai >= natoms)
      gmx_fatal(FARGS,"Position restraint atom index (%d) is larger than number of atoms in %s (%d).\n",ai+1,fn,natoms);
    for(j=0; j<DIM; j++) {
      if (j < npbcdim) {
	if (rc_scaling == erscALL) {
	  /* Convert from Cartesian to crystal coordinates */
	  x[ai][j] *= invbox[j][j];
	  for(k=j+1; k<npbcdim; k++)
	    x[ai][j] += invbox[k][j]*x[ai][k];
	} else if (rc_scaling == erscCOM) {
	  x[ai][j] -= com[j];
	}
      }
      pr->param[i].c[offset + j] = x[ai][j];
    }
  }

  if (rc_scaling == erscCOM) {
    /* Convert the COM from Cartesian to crystal coordinates */
    for(j=0; j<npbcdim; j++) {
      com[j] *= invbox[j][j];
      for(k=j+1; k<npbcdim; k++)
	com[j] += invbox[k][j]*com[k];
    }
  }
  
  free_t_atoms(&dumat);
  sfree(x);
  sfree(v);
}

static void gen_posres(t_params *pr, char *fnA, char *fnB,
		       int rc_scaling, int ePBC, t_atoms *atoms,
		       rvec com, rvec comB)
{
  int i,j;

  read_posres(pr,fnA,2*DIM,rc_scaling,ePBC,atoms->nr,atoms->atom,com);
  if (strcmp(fnA,fnB) == 0) {
    for(i=0; (i<pr->nr); i++)
      for(j=0; (j<DIM); j++)
	pr->param[i].c[3*DIM + j] = pr->param[i].c[2*DIM + j];
    copy_rvec(com,comB);
  } else {
    read_posres(pr,fnB,3*DIM,rc_scaling,ePBC,atoms->nr,atoms->atom,comB);
  }
}

static int search_atomtypes(t_atomtype *at,int *n,int typelist[],int thistype,
			    t_param param[],int ftype)
{
  int i,nn,nrfp,j,k,found;

  nn    = *n;
  nrfp  = NRFP(ftype);

  for(i=0; (i<nn); i++) 
  {
    if (typelist[i] == thistype)
    {
      /* This type number has already been added */
      break;
    }

    /* Otherwise, check if the parameters are identical to any previously added type */
    
    found=1;
    for(j=0;j<at->nr && found;j++) 
    {
      /* Check nonbonded parameters */
      for(k=0;k<nrfp && found;k++) 
      {
        found=(param[at->nr*typelist[i]+j].c[k]==param[at->nr*thistype+j].c[k]);
      }

      /* Check radius, volume, surftens */
      found = found && 
          ((at->radius[typelist[i]] == at->radius[thistype]) &&
           (at->vol[typelist[i]] == at->vol[thistype]) &&
           (at->surftens[typelist[i]] == at->surftens[thistype]) &&	
           (at->atomnumber[typelist[i]] ==	at->atomnumber[thistype]));
    }
    if (found)
    {
      break;
    }
  }
  
  if (i == nn) {
    if (debug)
      fprintf(debug,"Renumbering atomtype %d to %d\n",thistype,nn);
    if (nn == at->nr)
      gmx_fatal(FARGS,"Atomtype horror n = %d, %s, %d",nn,__FILE__,__LINE__);
    typelist[nn]=thistype;
    nn++;
  }
  *n = nn;
  
  return i;
}

static int renum_atype(t_params plist[],t_topology *top,int *wall_atomtype,
		       t_atomtype *at,bool bVerbose)
{
  int      i,j,k,l,mi,mj,nat,nrfp,ftype;
  t_param  *nbsnew;
  int      *typelist;
  real     *new_radius;
  real     *new_vol;
  real     *new_surftens;
  int      *new_atomnumber;
  
  snew(typelist,at->nr);

  if (bVerbose)
    fprintf(stderr,"renumbering atomtypes...\n");

  /* Since the bonded interactions have been assigned now,
   * we want to reduce the number of atom types by merging 
   * ones with identical nonbonded interactions, in addition
   * to removing unused ones.
   *
   * With Generalized-Born electrostatics, or implicit solvent
   * we also check that the atomtype radius, effective_volume
   * and surface tension match.
   *
   * With QM/MM we also check that the atom numbers match
   */
  
  /* Get nonbonded interaction type */
  if (plist[F_LJ].nr > 0)
    ftype=F_LJ;
  else
    ftype=F_BHAM;
   
  /* Renumber atomtypes by first making a list of which ones are actually used.
   * We provide the list of nonbonded parameters so search_atomtypes
   * can determine if two types should be merged. 
   */    
  nat=0;
  for(i=0; (i<top->atoms.nr); i++) {
    top->atoms.atom[i].type=
      search_atomtypes(at,&nat,typelist,top->atoms.atom[i].type,
		       plist[ftype].param,ftype);
    top->atoms.atom[i].typeB=
      search_atomtypes(at,&nat,typelist,top->atoms.atom[i].typeB,
		       plist[ftype].param,ftype);
  }

  for(i=0; i<2; i++) {
    if (wall_atomtype[i] >= 0)
      wall_atomtype[i] = search_atomtypes(at,&nat,typelist,wall_atomtype[i],
					  plist[ftype].param,ftype);
  }

  snew(new_radius,nat);
  snew(new_vol,nat);
  snew(new_surftens,nat);
  snew(new_atomnumber,nat);  

  /* We now have a list of unique atomtypes in typelist */

  if (debug)
    pr_ivec(debug,0,"typelist",typelist,nat,TRUE);
    
  /* Renumber nlist */ 
  nbsnew = NULL;
  snew(nbsnew,plist[ftype].nr);

  nrfp  = NRFP(ftype);
  
  for(i=k=0; (i<nat); i++)
  {
    mi=typelist[i];
    for(j=0; (j<nat); j++,k++) 
    {
      mj=typelist[j];
      for(l=0; (l<nrfp); l++)
      {
        nbsnew[k].c[l]=plist[ftype].param[at->nr*mi+mj].c[l];
      }
    }  
    new_radius[i]     = at->radius[mi];
    new_vol[i]        = at->vol[mi];
    new_surftens[i]   = at->surftens[mi];
    new_atomnumber[i] = at->atomnumber[mi];
  }
  
  for(i=0; (i<nat*nat); i++) {
    for(l=0; (l<nrfp); l++)
      plist[ftype].param[i].c[l]=nbsnew[i].c[l];
  }
  plist[ftype].nr=i;
  
  sfree(at->radius);
  sfree(at->vol);
  sfree(at->surftens);
  sfree(at->atomnumber);
  
  at->radius     = new_radius;
  at->vol        = new_vol;
  at->surftens   = new_surftens;
  at->atomnumber = new_atomnumber;
  
  at->nr=nat;

  sfree(nbsnew);
  sfree(typelist);
  
  return nat;
}

static void set_wall_atomtype(t_atomtype *at,t_gromppopts *opts,
			      t_inputrec *ir)
{
  int i;

  if (ir->nwall > 0)
    fprintf(stderr,"Searching the wall atom type(s)\n");
  for(i=0; i<ir->nwall; i++)
    ir->wall_atomtype[i] = at2type(opts->wall_atomtype[i],at);
}

static int count_constraints(t_params plist[])
{
  int count,i;

  count = 0;
  for(i=0; i<F_NRE; i++)
    if (i == F_SETTLE)
      count += 3*plist[i].nr;
    else if (interaction_function[i].flags & IF_CONSTRAINT)
      count += plist[i].nr;
  
  return count;
}

int main (int argc, char *argv[])
{
  static char *desc[] = {
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
  t_topology   *sys;
  t_molinfo    msys;
  t_atomtype   *atype;
  t_inputrec   *ir;
  int          natoms,nvsite,nc,comb;
  t_params     *plist;
  t_state      state;
  matrix       box;
  real         max_spacing,reppow;
  char         fn[STRLEN],fnB[STRLEN],*mdparin;
  int          nerror;
  bool         bNeedVel,bGenVel;
  bool         have_radius,have_vol,have_surftens;
  bool         have_atomnumber;
  
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
    { efENX, "-e",  NULL,        ffOPTRD }
  };
#define NFILE asize(fnm)

  /* Command line options */
  static bool bVerbose=TRUE,bRenum=TRUE;
  static bool bRmVSBds=TRUE,bCheckPairs=FALSE;
  static int  i,maxwarn=10;
  static real fr_time=-1;
  t_pargs pa[] = {
    { "-v",       FALSE, etBOOL, {&bVerbose},
      "Be loud and noisy" },
    { "-time",    FALSE, etREAL, {&fr_time},
      "Take frame at or first after this time." },
    { "-rmvsbds",FALSE, etBOOL, {&bRmVSBds},
      "Remove constant bonded interactions with virtual sites" },
    { "-maxwarn", FALSE, etINT,  {&maxwarn},
      "Number of warnings after which input processing stops" },
    { "-check14", FALSE, etBOOL, {&bCheckPairs},
      "Remove 1-4 interactions without Van der Waals" },
    { "-renum",   FALSE, etBOOL, {&bRenum},
      "Renumber atomtypes and minimize number of atomtypes" },
  };
  
  CopyRight(stdout,argv[0]);
  
  /* Initiate some variables */
  nerror=0;
  snew(ir,1);
  snew(opts,1);
  init_ir(ir,opts);
  
  /* Parse the command line */
  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  init_warning(maxwarn);
  
  /* PARAMETER file processing */
  mdparin = opt2fn("-f",NFILE,fnm);
  set_warning_line(mdparin,-1);    
  get_ir(mdparin,opt2fn("-po",NFILE,fnm),ir,opts,&nerror);
  
  if (bVerbose) 
    fprintf(stderr,"checking input for internal consistency...\n");
  check_ir(ir,opts,&nerror);

  if (ir->ld_seed == -1) {
    ir->ld_seed = make_seed();
    fprintf(stderr,"Setting the LD random seed to %d\n",ir->ld_seed);
  }

  bNeedVel = EI_STATE_VELOCITY(ir->eI);
  bGenVel  = (bNeedVel && opts->bGenVel);

  snew(plist,F_NRE);
  init_plist(plist);
  snew(sys,1);
  snew(atype,1);
  if (debug)
    pr_symtab(debug,0,"Just opened",&sys->symtab);
    
  strcpy(fn,ftp2fn(efTOP,NFILE,fnm));
  if (!fexist(fn)) 
    gmx_fatal(FARGS,"%s does not exist",fn);
  new_status(fn,opt2fn_null("-pp",NFILE,fnm),opt2fn("-c",NFILE,fnm),
	     opts,ir,bGenVel,bVerbose,&state,
	     atype,sys,&msys,plist,&comb,&reppow,
	     (opts->eDisre==edrEnsemble),opts->bMorse,
	     bCheckPairs,&nerror);
  
  if (debug)
    pr_symtab(debug,0,"After new_status",&sys->symtab);
  
  nc = count_constraints(msys.plist);
  if ((ir->eI == eiCG) && nc) {
    fprintf(stderr,
	    "ERROR: can not do Conjugate Gradients with constraints (%d)\n",
	    nc);
    nerror++;
  }
  if (ir->bPeriodicMols && (ir->eConstrAlg == estSHAKE) && nc) {
    fprintf(stderr,
	    "ERROR: can not do periodic molecules with %s, use %s\n",
	    eshake_names[estSHAKE],eshake_names[estLINCS]);
    nerror++;
  }

  /* If we are doing GBSA, check that we got the parameters we need */
  have_radius=have_vol=have_surftens=TRUE;
  for(i=0;i<atype->nr;i++) {
    have_radius=have_radius && (atype->radius[i]>0);
    have_vol=have_vol && (atype->vol[i]>0);
    have_surftens=have_surftens && (atype->surftens[i]>=0);
  }
  if(!have_radius && ir->coulombtype==eelGB) {
    fprintf(stderr,"Can't do GB electrostatics; the forcefield is missing values for\n"
		"atomtype radii, or they might be zero.");
    nerror++;
  }
  if(!have_vol && ir->gb_algorithm==egbKARPLUS) {
    fprintf(stderr,"Can't calculate Karplus Born radii; the forcefield is missing values\n"
		" for atomtype effective volumes, or they might be zero.");
    nerror++;
  }
  if(!have_surftens && ir->implicit_solvent!=eisNO) {
    fprintf(stderr,"Can't do implicit solvent; the forcefield is missing values\n"
		" for atomtype surface tension.");
    nerror++;
  }
  
  /* If we are doing QM/MM, check that we got the atom numbers */
  have_atomnumber = TRUE;
  for (i=0;i<atype->nr;i++) {
    have_atomnumber = have_atomnumber && (atype->atomnumber[i]>=0);
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
    print_warn_num();
    
    gmx_fatal(FARGS,"There were %d error(s) processing your input",nerror);
  }
  if (opt2bSet("-r",NFILE,fnm))
    sprintf(fn,opt2fn("-r",NFILE,fnm));
  else
    sprintf(fn,opt2fn("-c",NFILE,fnm));
  if (opt2bSet("-rb",NFILE,fnm))
    sprintf(fnB,opt2fn("-rb",NFILE,fnm));
  else
    strcpy(fnB,fn);

  if (msys.plist[F_POSRES].nr > 0) {
    if (bVerbose) {
      fprintf(stderr,"Reading position restraint coords from %s",fn);
      if (strcmp(fn,fnB) ==0) {
	fprintf(stderr,"\n");
      } else {
	fprintf(stderr," and %s\n",fnB);
      }
    }
    gen_posres(&(msys.plist[F_POSRES]),fn,fnB,
	       ir->refcoord_scaling,ir->ePBC,&sys->atoms,
	       ir->posres_com,ir->posres_comB);
  }
  
  /* set parameters for virtual site construction */
  nvsite=set_vsites(bVerbose, &sys->atoms, atype, msys.plist);
  /* now throw away all obsolete bonds, angles and dihedrals: */
  /* note: constraints are ALWAYS removed */
  if (nvsite)
    clean_vsite_bondeds(msys.plist,sys->atoms.nr,bRmVSBds);
  
  set_wall_atomtype(atype,opts,ir);
  if (bRenum) 
    atype->nr = renum_atype(plist, sys, ir->wall_atomtype, atype, bVerbose);
  
  /* Copy the atomtype data to the topology atomtype list */
  sys->atomtypes.nr=atype->nr;
  snew(sys->atomtypes.radius,atype->nr);
  snew(sys->atomtypes.vol,atype->nr);
  snew(sys->atomtypes.surftens,atype->nr);
  snew(sys->atomtypes.atomnumber,atype->nr);


  for(i=0;i<atype->nr;i++) {
    sys->atomtypes.radius[i]=atype->radius[i];
    sys->atomtypes.vol[i]=atype->vol[i];
    sys->atomtypes.surftens[i]=atype->surftens[i];
    sys->atomtypes.atomnumber[i] = atype->atomnumber[i];
  }

  if (debug)
    pr_symtab(debug,0,"After renum_atype",&sys->symtab);

  if (bVerbose) 
    fprintf(stderr,"converting bonded parameters...\n");
  convert_params(atype->nr, plist, msys.plist, comb, reppow, &sys->idef);
  
  if (debug)
    pr_symtab(debug,0,"After convert_params",&sys->symtab);

  /* set ptype to VSite for virtual sites */
  if (nvsite) {
    set_vsites_ptype(bVerbose,&sys->idef,&sys->atoms);
    if (debug)
      pr_symtab(debug,0,"After virtual sites",&sys->symtab);
  }
  /* Check velocity for virtual sites and shells */
  if (bGenVel) 
    check_vel(&sys->atoms,state.v);
    
  /* check masses */
  check_mol(&(sys->atoms));
  
  check_warning_error(FARGS);

  /* Now build the shakeblocks from the shakes */
  gen_sblocks(stdout,sys->atoms.nr,&(sys->idef),
	      &(sys->blocks[ebSBLOCKS]),FALSE);
  if (debug)
    pr_symtab(debug,0,"After gen_sblocks",&sys->symtab);
   
  if (bVerbose) 
    fprintf(stderr,"initialising group options...\n");
  do_index(ftp2fn_null(efNDX,NFILE,fnm),
	   &sys->symtab,&(sys->atoms),bVerbose,ir,&sys->idef,
	   bGenVel ? state.v : NULL);

  /* Init the temperature coupling state */
  init_gtc_state(&state,ir->opts.ngtc);

  if (bVerbose)
    fprintf(stderr,"Checking consistency between energy and charge groups...\n");
  check_eg_vs_cg(&(sys->atoms),&(sys->blocks[ebCGS]));
  
  if (debug)
    pr_symtab(debug,0,"After index",&sys->symtab);
  triple_check(mdparin,ir,sys,&nerror);
  close_symtab(&sys->symtab);
  if (debug)
    pr_symtab(debug,0,"After close",&sys->symtab);

  /* make exclusions between QM atoms */
  if(ir->bQMMM)
    generate_qmexcl(sys,ir);

  if (ftp2bSet(efTRN,NFILE,fnm)) {
    if (bVerbose)
      fprintf(stderr,"getting data from old trajectory ...\n");
    cont_status(ftp2fn(efTRN,NFILE,fnm),ftp2fn_null(efENX,NFILE,fnm),
		bNeedVel,bGenVel,fr_time,ir,&state,sys);
  }

  if (ir->ePBC==epbcXY && ir->nwall!=2)
    clear_rvec(state.box[ZZ]);
  
  if ((ir->coulombtype == eelPPPM) || (ir->coulombtype == eelPME) || 
      (ir->coulombtype == eelPMEUSER)|| (ir->coulombtype == eelEWALD)) {
    /* Calculate the optimal grid dimensions */
    copy_mat(state.box,box);
    if (ir->ePBC==epbcXY && ir->nwall==2)
      svmul(ir->wall_ewald_zfac,box[ZZ],box[ZZ]);
    max_spacing = calc_grid(stdout,box,opts->fourierspacing,
			    &(ir->nkx),&(ir->nky),&(ir->nkz),1);
    if ((ir->coulombtype == eelPPPM) && (max_spacing > 0.1)) {
      set_warning_line(mdparin,-1);
      sprintf(warn_buf,"Grid spacing larger then 0.1 while using PPPM.");
      warning(NULL);
    }
  }
  /*  reset_multinr(sys); */
  
  {
    double cio = compute_io(ir,&sys->atoms,F_NRE,1);
    sprintf(warn_buf,"This run will generate roughly %.0f Mb of data",cio);
    if (cio > 2000)
      warning(NULL);
    else
      printf("%s\n",warn_buf);
  }
    
  if (bVerbose) 
    fprintf(stderr,"writing run input file...\n");

  state.lambda = ir->init_lambda;
  write_tpx_state(ftp2fn(efTPX,NFILE,fnm),0,ir->init_t,ir,&state,sys);
  print_warn_num();
  
  thanx(stderr);
  
  return 0;
}
