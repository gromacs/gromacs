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
 * GRowing Old MAkes el Chrono Sweat
 */
static char *SRCID_grompp_c = "$Id$";

#include <sys/types.h>
#include <math.h>
#include <string.h>
#include <assert.h>
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
#include "fatal.h"
#include "rdgroup.h"
#include "gmxfio.h"
#include "trnio.h"
#include "tpxio.h"
#include "dum_parm.h"
#include "txtdump.h"
#include "calcgrid.h"

void check_solvent(bool bVerbose,t_molinfo msys[],
		   int Nsim,t_simsystem Sims[],t_inputrec *ir,char *SolventOpt)
{
  int  i,wmol,nwt;
  char buf[128];
  
  ir->solvent_opt=-1;
  if (!SolventOpt || strlen(SolventOpt)==0) {
    if (bVerbose)
      fprintf(stderr,"no solvent optimizations...\n");
  }
  else {
    if (bVerbose)
      fprintf(stderr,"checking for solvent...\n");
    for(i=0; (i<Nsim); i++) {
      wmol = Sims[i].whichmol;
      if ((strcmp(SolventOpt,*(msys[wmol].name)) == 0) && 
	  (Sims[i].nrcopies > 0)) {
	nwt = msys[wmol].atoms.atom[0].type;
	if (ir->solvent_opt == -1) {
	  if (msys[wmol].atoms.nr == 3)
	    ir->solvent_opt=nwt;
	  else {
	    sprintf(buf,"Sorry, can only do solvent optimization with SPC-like models\n");
	    warning(buf);
	  }
	}
	else if (ir->solvent_opt == nwt) {
	  if (debug)
	    fprintf(debug,"Remark: Multiple topology entries for %s\n",
		    SolventOpt);
	} else
	  fatal_error(0,"Multiple non-matching topology entries for %s",
		      SolventOpt);
      }
    }
    if (bVerbose) {
      if (ir->solvent_opt != -1)
	fprintf(stderr,"...using solvent optimization for atomtype %d\n",
		ir->solvent_opt);
      else
	fprintf(stderr,"...no solvent\n");
    }
  }
}

static int *shuffle_xv(char *ndx,bool bSort,bool bVerbose,
		       int ntab,int *tab,int nmol,t_molinfo *mol,
		       int natoms,rvec *x,rvec *v,
		       int Nsim,t_simsystem Sims[])
{
  FILE *out;
  rvec *xbuf,*vbuf;
  int  *nindex,**index,xind,*done,*forward,*backward;
  int  i,j,j0,k,n,mi,nnat;

  fprintf(stderr,"Entering shuffle_xv\n");
  /* Determine order in old x array! 
   * the index array holds for each molecule type
   * a pointer to the position at which the coordinates start
   */
  snew(index,nmol);
  snew(nindex,nmol);
  snew(done,nmol);
  
  /* Compute the number of copies for each molecule type */
  for(i=0; (i<Nsim); i++) {
    mi=Sims[i].whichmol;
    nindex[mi]+=Sims[i].nrcopies;
  }
  /* Allocate space */
  for(i=0; (i<nmol); i++) {
    snew(index[i],nindex[i]);
    nindex[i]=0;
  }
  xind=0;
  for(i=0; (i<Nsim); i++) {
    /* Mol-index */
    mi=Sims[i].whichmol;
    
    /* Current number of molecules processed for this mol-type */
    k=nindex[mi];
    
    /* Number of atoms in this mol-type */
    nnat = mol[mi].atoms.nr;
    
    for(j=0; (j<Sims[i].nrcopies); j++,k++) {
      index[mi][k]=xind;
      xind += nnat;
    }
    nindex[mi]=k;
  }
  assert(xind == natoms);
  
  /* Buffers for x and v */  
  snew(xbuf,natoms);
  snew(vbuf,natoms);
  
  /* Sort the coordinates if necessary */
  if (bSort) {
    for(i=0; (i<nmol); i++) {
      if (bVerbose)
	fprintf(stderr,"Sorting coordinates for %5d copies of molecule %s\n",
		nindex[i],*mol[i].name);
      nnat = mol[i].atoms.nr;
      /* Copy molecules into buffer arrays */
      for(j=n=0; (j<nindex[i]); j++) {
	for(k=0; (k<nnat); k++,n++) {
	  copy_rvec(x[index[i][j]+k],xbuf[n]);
	  copy_rvec(v[index[i][j]+k],vbuf[n]);
	}
      }
      /* Sort the coords */
      sortwater(0,j,nnat,xbuf,vbuf);
      /* Copy molecules back from buffer arrays */
      for(j=n=0; (j<nindex[i]); j++) {
	for(k=0; (k<nnat); k++,n++) {
	  copy_rvec(xbuf[n],x[index[i][j]+k]);
	  copy_rvec(vbuf[n],v[index[i][j]+k]);
	}
      }
    }
  }
  
  /* Make a forward shuffle array, i.e. the old numbers of
   * the current order: this makes a shuffled order from the
   * original.
   * Simultaneously copy the coordinates..
   */
  snew(forward,natoms);
  for(i=k=0; (i<ntab); i++) {
    /* Get the molecule type */
    mi   = tab[i];
    
    /* Determine number of atoms in thsi mol-type */
    nnat = mol[mi].atoms.nr;
    
    /* Find the starting index in the x & v arrays */
    j0   = index[mi][done[mi]];
    
    /* Copy the coordinates */
    for(j=j0; (j<j0+nnat); j++,k++) {
      copy_rvec(x[j],xbuf[k]);
      copy_rvec(v[j],vbuf[k]);
      /* Store the old index of the new one */
      forward[k]=j;
    }
    /* Increment the number of molecules processed for this type */
    done[mi]++;
  }

  /* Now copy the buffers back to the original x and v */
  for(i=0; (i<natoms); i++) {
    copy_rvec(xbuf[i],x[i]);
    copy_rvec(vbuf[i],v[i]);
  }
  
  /* Delete buffers */
  sfree(xbuf);
  sfree(vbuf);
  
  /* Now make an inverse shuffle index:
   * this transforms the new order into the original one...
   */
  snew(backward,natoms);
  for(i=0; (i<natoms); i++)
    backward[forward[i]] = i;
    
  /* Make an index file for deshuffling the atoms */
  out=ffopen(ndx,"w");
  fprintf(out,"1  %d\nDeShuffle  %d\n",natoms,natoms);
  for(i=0; (i<natoms); i++)
    fprintf(out,"  %d",backward[i]);
  fprintf(out,"\n");
  ffclose(out);
  
  sfree(backward);
  for(i=0; (i<nmol); i++)
    sfree(index[i]);
  sfree(index);
  sfree(done);
  
  return forward;
}

int rm_disre(int nrmols,t_molinfo mols[])
{
  int  i,n;
  
  n=0;
  /* For all the molecule types */
  for(i=0; (i<nrmols); i++) {
    n+=mols[i].plist[F_DISRES].nr;
    mols[i].plist[F_DISRES].nr=0;
  }
  return n;
}

static int check_atom_names(char *fn1, char *fn2, t_atoms *at1, t_atoms *at2,
			    int *forward)
{
  int i,nmismatch,idx;
#define MAXMISMATCH 20

  assert(at1->nr==at2->nr);
  
  nmismatch=0;
  for(i=0; i < at1->nr; i++) {
    if(forward)
      idx=forward[i];
    else
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

static int *new_status(char *topfile,char *topppfile,char *confin,
		       char *ndxout,
		       t_gromppopts *opts,t_inputrec *ir,
		       bool bGenVel,bool bVerbose,
		       bool bSort,int *natoms,
		       rvec **x,rvec **v,matrix box,
		       t_atomtype *atype,t_topology *sys,
		       t_molinfo *msys,t_params plist[],
		       int nprocs,bool bEnsemble,bool bMorse,int *nerror)
{
  t_molinfo   *molinfo=NULL;
  t_simsystem *Sims=NULL;
  t_atoms     *confat;
  int         *forward=NULL;
  int         i,nrmols,Nsim,nmismatch;
  int         ntab,*tab;
  char        buf[STRLEN];

  init_top(sys);
  init_molinfo(msys);
  
  /* TOPOLOGY processing */
  msys->name=do_top(bVerbose,topfile,topppfile,opts,&(sys->symtab),
		    plist,atype,&nrmols,&molinfo,ir,&Nsim,&Sims);
  
  check_solvent(bVerbose,molinfo,Nsim,Sims,ir,opts->SolventOpt);
  
  ntab = 0;
  tab  = NULL;
  if (nprocs > 1) {
    tab=mk_shuffle_tab(nrmols,molinfo,nprocs,&ntab,Nsim,Sims,bVerbose);
    if (debug) {
      for(i=0; (i<ntab); i++)
	fprintf(debug,"Mol[%5d] = %s\n",i,*molinfo[tab[i]].name);
      fflush(debug);
    }
    fprintf(stderr,"Made a shuffling table with %d entries [molecules]\n",
	    ntab);
  }
  if (bMorse)
    convert_harmonics(nrmols,molinfo,atype);
  
  if (opts->eDisre==edrNone) {
    i=rm_disre(nrmols,molinfo);
    if (bVerbose && i)
      fprintf(stderr,"removed %d distance restraints\n",i);
  }
  
  topcat(msys,nrmols,molinfo,ntab,tab,Nsim,Sims,bEnsemble);
  
  /* Copy structures from msys to sys */
  mi2top(sys,msys);
  
  /* COORDINATE file processing */
  if (bVerbose) 
    fprintf(stderr,"processing coordinates...\n");
  
  get_stx_coordnum(confin,natoms);
  if (*natoms != sys->atoms.nr)
    fatal_error(0,"number of coordinates in coordinate file (%s, %d)\n"
		"             does not match topology (%s, %d)",
		confin,*natoms,topfile,sys->atoms.nr);
  else {
    /* make space for coordinates and velocities */
    snew(confat,1);
    init_t_atoms(confat,*natoms,FALSE);
    snew(*x,*natoms);
    snew(*v,*natoms);
    read_stx_conf(confin,opts->title,confat,*x,*v,box);

    if (ntab > 0) {
      if (bVerbose)
	fprintf(stderr,"Shuffling coordinates...\n");
      forward=shuffle_xv(ndxout,bSort,bVerbose,
			 ntab,tab,nrmols,molinfo,
			 *natoms,*x,*v,Nsim,Sims);
    }
    
    nmismatch=check_atom_names(topfile, confin, &(sys->atoms), confat,forward);
    free_t_atoms(confat);
    sfree(confat);
    
    if (nmismatch) {
      sprintf(buf,"%d non-matching atom name%s\n",nmismatch,
	      (nmismatch == 1) ? "" : "s");
      warning(buf);
    }    
    if (bVerbose) 
      fprintf(stderr,"double-checking input for internal consistency...\n");
    double_check(ir,box,msys,nerror);
  }
  
  if (bGenVel) {
    real *mass;
    
    snew(mass,msys->atoms.nr);
    for(i=0; (i<msys->atoms.nr); i++)
      mass[i]=msys->atoms.atom[i].m;
    
    maxwell_speed(opts->tempi,sys->atoms.nr*DIM,
		  opts->seed,&(sys->atoms),*v);
    stop_cm(stdout,sys->atoms.nr,mass,*x,*v);
    sfree(mass);
  }
  for(i=0; (i<nrmols); i++)
    done_mi(&(molinfo[i]));
  sfree(molinfo);
  sfree(Sims);
  
  return forward;
}

static void cont_status(char *slog,bool bNeedVel,bool bGenVel, real fr_time,
			t_inputrec *ir,int *natoms,
			rvec **x,rvec **v,matrix box,
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
  
  *natoms = fr.natoms;

  if(sys->atoms.nr != *natoms)
    fatal_error(0,"Number of atoms in Topology "
		"is not the same as in Trajectory");

  /* Find the appropriate frame */
  while ((fr_time == -1 || fr.time < fr_time) && read_next_frame(fp,&fr));
  
  close_trj(fp);

  if (fr.not_ok & FRAME_NOT_OK)
    fatal_error(0,"Can not start from an incomplete frame");

  *x = fr.x;
  if (bNeedVel && !bGenVel)
    *v = fr.v;
  copy_mat(fr.box,box);

  fprintf(stderr,"Using frame at t = %g ps\n",fr.time);
  fprintf(stderr,"Starting time for run is %g ps\n",ir->init_t); 
}

static void gen_posres(t_params *pr,char *fn)
{
  rvec   *x,*v;
  t_atoms dumat;
  matrix box;
  int    natoms;
  char   title[256];
  int    i,ai,j;
  
  get_stx_coordnum(fn,&natoms);
  snew(x,natoms);
  snew(v,natoms);
  init_t_atoms(&dumat,natoms,FALSE);
  read_stx_conf(fn,title,&dumat,x,v,box);
  
  for(i=0; (i<pr->nr); i++) {
    ai=pr->param[i].AI;
    if (ai >= natoms)
      fatal_error(0,"Position restraint atom index (%d) is larger than natoms (%d)\n",
		  ai+1,natoms);
    for(j=0; (j<DIM); j++)
      pr->param[i].c[j+DIM]=x[ai][j];
  }
  /*pr->nrfp+=DIM;*/
  
  free_t_atoms(&dumat);
  sfree(x);
  sfree(v);
}

static int search_array(int atnr,int *n,int map[],int key)
{
  int i,nn;
  
  nn = *n;
  for(i=0; (i<nn); i++)
    if (map[i] == key)
      break;
  
  if (i == nn) {
    if (debug)
      fprintf(debug,"Renumbering atomtype %d to %d\n",key,nn);
    if (nn == atnr)
      fatal_error(0,"Atomtype horror n = %d, %s, %d",nn,__FILE__,__LINE__);
    map[nn]=key;
    nn++;
  }
  *n = nn;
  
  return i;
}

static int renum_atype(t_params plist[],t_topology *top,
		       int atnr,t_inputrec *ir,bool bVerbose)
{
  int      i,j,k,l,mi,mj,nat,nrfp,ftype;
  t_param  *nbsnew;
  int      *map;

  snew(map,atnr);
  if (bVerbose)
    fprintf(stderr,"renumbering atomtypes...\n");
  /* Renumber atomtypes and meanwhile make a list of atomtypes */    
  nat=0;
  for(i=0; (i<top->atoms.nr); i++) {
    top->atoms.atom[i].type=
      search_array(atnr,&nat,map,top->atoms.atom[i].type);
    top->atoms.atom[i].typeB=
      search_array(atnr,&nat,map,top->atoms.atom[i].typeB);
  }
  
  if (ir->solvent_opt != -1) {
    if (debug)
      fprintf(debug,"Solvent_type before: %d\n",ir->solvent_opt);
    for(j=0; (j<nat); j++)
      if (map[j] == ir->solvent_opt) {
	ir->solvent_opt=j;
	break;
      }
    if (debug)
      fprintf(debug,"Renumbering solvent_opt (atomtype for OW) to %d\n",
	      ir->solvent_opt);
  }
  if (debug)
    pr_ivec(debug,0,"map",map,nat);
    
  /* Renumber nlist */
  if (plist[F_LJ].nr > 0)
    ftype=F_LJ;
  else
    ftype=F_BHAM;
    
  nbsnew = NULL;
  snew(nbsnew,plist[ftype].nr);
  nrfp  = NRFP(ftype);
  
  for(i=k=0; (i<nat); i++) {
    mi=map[i];
    for(j=0; (j<nat); j++,k++) {
      mj=map[j];
      for(l=0; (l<nrfp); l++)
	nbsnew[k].c[l]=plist[ftype].param[atnr*mi+mj].c[l];
    }
  }
  for(i=0; (i<nat*nat); i++) {
    for(l=0; (l<nrfp); l++)
      plist[ftype].param[i].c[l]=nbsnew[i].c[l];
  }
  plist[ftype].nr=i;
  
  sfree(nbsnew);
  sfree(map);
  
  return nat;
}

int count_constraints(t_params plist[])
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
    
    "grompp calls the c-preprocessor to resolve includes, macros ",
    "etcetera. To specify a macro-preprocessor other than /lib/cpp ",
    "(such as m4)",
    "you can put a line in your parameter file specifying the path",
    "to that cpp. Specifying [TT]-pp[tt] will get the pre-processed",
    "topology file written out.[PAR]",
    
    "If your system does not have a c-preprocessor, you can still",
    "use grompp, but you do not have access to the features ",
    "from the cpp. Command line options to the c-preprocessor can be given",
    "in the [TT].mdp[tt] file. See your local manual (man cpp).[PAR]",
    
    "When using position restraints a file with restraint coordinates",
    "can be supplied with [TT]-r[tt], otherwise constraining will be done",
    "relative to the conformation from the [TT]-c[tt] option.[PAR]",
    
    "Starting coordinates can be read from trajectory with [TT]-t[tt].",
    "The last frame with coordinates and velocities will be read,",
    "unless the [TT]-time[tt] option is used.",
    "Note that these velocities will not be used when [TT]gen_vel = yes[tt]",
    "in your [TT].mdp[tt] file. If you want to continue a crashed run, it is",
    "easier to use [TT]tpbconv[tt].[PAR]",

    "When preparing an input file for parallel [TT]mdrun[tt] it may",
    "be advantageous to partition the simulation system over the",
    "processors in a way in which each processor ahs a similar amount of",
    "work. The -shuffle option does just that. For a single protein",
    "in water this does not make a difference, however for a system where",
    "you have many copies of different molecules  (e.g. liquid mixture",
    "or membrane/water system) the option is definitely a must.[PAR]",
    
    "A further optimization for parallel systems is the [TT]-sort[tt]",
    "option which sorts molecules according to coordinates. This must",
    "always be used in conjunction with [TT]-shuffle[tt], however",
    "sorting also works when you have only one molecule type.[PAR]", 
        
    "Using the [TT]-morse[tt] option grompp can convert the harmonic bonds",
    "in your topology to morse potentials. This makes it possible to break",
    "bonds. For this option to work you need an extra file in your $GMXLIB",
    "with dissociation energy. Use the -debug option to get more information",
    "on the workings of this option (look for MORSE in the grompp.log file",
    "using less or something like that).[PAR]",
    
    "By default all bonded interactions which have constant energy due to",
    "dummy atom constructions will be removed. If this constant energy is",
    "not zero, this will result in a shift in the total energy. All bonded",
    "interactions can be kept by turning off [TT]-rmdumbds[tt]. Additionally,",
    "all constraints for distances which will be constant anyway because",
    "of dummy atom constructions will be removed. If any constraints remain",
    "which involve dummy atoms, a fatal error will result.[PAR]"
    
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
  t_atomtype   atype;
  t_inputrec   *ir;
  int          natoms,ndum;
  int          *forward=NULL;
  t_params     *plist;
  rvec         *x=NULL,*v=NULL;
  matrix       box;
  real         max_spacing;
  char         fn[STRLEN],*mdparin;
  int          nerror;
  bool         bNeedVel,bGenVel;

  t_filenm fnm[] = {
    { efMDP, NULL,  NULL,        ffREAD  },
    { efMDP, "-po", "mdout",     ffWRITE },
    { efSTX, "-c",  NULL,        ffREAD  },
    { efSTX, "-r",  NULL,        ffOPTRD },
    { efNDX, NULL,  NULL,        ffOPTRD },
    { efNDX, "-deshuf", "deshuf",ffOPTWR },
    { efTOP, NULL,  NULL,        ffREAD  },
    { efTOP, "-pp", "processed", ffOPTWR },
    { efTPX, "-o",  NULL,        ffWRITE },
    { efTRN, "-t",  NULL,        ffOPTRD }
  };
#define NFILE asize(fnm)

  /* Command line options */
  static bool bVerbose=TRUE,bRenum=TRUE,bShuffle=TRUE;
  static bool bRmDumBds=TRUE,bSort=FALSE;
  static int  nprocs=1,maxwarn=10;
  static real fr_time=-1;
  t_pargs pa[] = {
    { "-v",       FALSE, etBOOL, {&bVerbose},
      "Be loud and noisy" },
    { "-time",    FALSE, etREAL, {&fr_time},
      "Take frame at or first after this time." },
    { "-np",      FALSE, etINT,  {&nprocs},
      "Generate statusfile for # processors" },
    { "-shuffle", FALSE, etBOOL, {&bShuffle},
      "Shuffle molecules over processors" },
    { "-sort",    FALSE, etBOOL, {&bSort},
      "Sort molecules according to X coordinate" },
    { "-rmdumbds",FALSE, etBOOL, {&bRmDumBds},
      "Remove constant bonded interactions with dummies" },
    { "-maxwarn", FALSE, etINT,  {&maxwarn},
      "Number of warnings after which input processing stops" },
    { "-renum",   FALSE, etBOOL, {&bRenum},
      "HIDDENRenumber atomtypes and minimize number of atomtypes" },
  };
  
  CopyRight(stdout,argv[0]);
  
  /* Initiate some variables */
  nerror=0;
  snew(ir,1);
  snew(opts,1);
  init_ir(ir,opts);
  
  /* Parse the command line */
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,0,NULL);
  
  if ((nprocs > 0) && (nprocs <= MAXPROC))
    printf("creating statusfile for %d processor%s...\n",
	   nprocs,nprocs==1?"":"s");
  else 
    fatal_error(0,"invalid number of processors %d\n",nprocs);
    
  if (bShuffle && opt2bSet("-r",NFILE,fnm)) {
    fprintf(stderr,"Can not shuffle and do position restraints, "
	    "turning off shuffle\n");
    bShuffle=FALSE;
  }
	       
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

  if (bShuffle && (opts->eDisre==edrEnsemble)) {
    fprintf(stderr,"Can not shuffle and do ensemble averaging, "
	    "turning off shuffle\n");
    bShuffle=FALSE;
  }
  if (bSort && !bShuffle) 
    warning("Can not do sorting without shuffling. Sorting turned off.");
    
  if ( (nprocs > 1) && (ir->nstcomm < 0) ) {
    /* this check is also in md.c, if it becomes obsolete here,
       also remove it there */
    fprintf(stderr,
	    "ERROR: removing rotation around center of mass (nstcomm=%d)"
	    "in a parallel run (np=%d) not implemented\n",ir->nstcomm,nprocs);
    nerror++;
  }
    
  bNeedVel = (ir->eI == eiMD || ir->eI == eiSD);
  bGenVel  = (bNeedVel && opts->bGenVel);

  snew(plist,F_NRE);
  init_plist(plist);
  snew(sys,1);
  if (debug)
    pr_symtab(debug,0,"Just opened",&sys->symtab);
    
  strcpy(fn,ftp2fn(efTOP,NFILE,fnm));
  if (!fexist(fn)) 
    fatal_error(0,"%s does not exist",fn);
  forward=new_status(fn,opt2fn_null("-pp",NFILE,fnm),opt2fn("-c",NFILE,fnm),
		     opt2fn("-deshuf",NFILE,fnm),
		     opts,ir,bGenVel,bVerbose,bSort,&natoms,&x,&v,box,
		     &atype,sys,&msys,plist,bShuffle ? nprocs : 1,
		     (opts->eDisre==edrEnsemble),opts->bMorse,&nerror);
  
  if (debug)
    pr_symtab(debug,0,"After new_status",&sys->symtab);
  
  if (ir->eI == eiCG) {
    int nc;
    nc = count_constraints(msys.plist);
    if (nc) {
      fprintf(stderr,
	      "ERROR: can not do Conjugate Gradients with constraints (%d)\n",
	      nc);
      nerror++;
    }
  }

  if (nerror) {
    print_warn_num();
    if (nerror==1)
      fatal_error(0,"There was %d error",nerror);
    else
      fatal_error(0,"There were %d errors",nerror);
  }
  if (opt2bSet("-r",NFILE,fnm))
    sprintf(fn,opt2fn("-r",NFILE,fnm));
  else
    sprintf(fn,opt2fn("-c",NFILE,fnm));
  
  if (msys.plist[F_POSRES].nr > 0) {
    if (bVerbose)
      fprintf(stderr,"Reading position restraint coords from %s\n",fn);
    gen_posres(&(msys.plist[F_POSRES]),fn);
  }
  
  /* set parameters for Dummy construction */
  ndum=set_dummies(bVerbose, &sys->atoms, atype, msys.plist);
  /* now throw away all obsolete bonds, angles and dihedrals: */
  /* note: constraints are ALWAYS removed */
  if (ndum)
    clean_dum_bondeds(msys.plist,sys->atoms.nr,bRmDumBds);
  
  if (bRenum) 
    atype.nr=renum_atype(plist, sys, atype.nr, ir, bVerbose);
  
  if (debug)
    pr_symtab(debug,0,"After renum_atype",&sys->symtab);

  if (bVerbose) 
    fprintf(stderr,"converting bonded parameters...\n");
  convert_params(atype.nr, plist, msys.plist, &sys->idef);
  
  if (debug)
    pr_symtab(debug,0,"After convert_params",&sys->symtab);

  /* set ptype to Dummy for dummy atoms */
  if (ndum) {
    set_dummies_ptype(bVerbose,&sys->idef,&sys->atoms);
    if (debug)
      pr_symtab(debug,0,"After dummy",&sys->symtab);
  }
  
  /* check masses */
  check_mol(&(sys->atoms));
  
  /* Now build the shakeblocks from the shakes */
  gen_sblocks(bVerbose,sys->atoms.nr,&(sys->idef),
	      &(sys->blocks[ebSBLOCKS]),FALSE);
  if (debug)
    pr_symtab(debug,0,"After gen_sblocks",&sys->symtab);
   
  if (bVerbose) 
    fprintf(stderr,"initialising group options...\n");
  do_index(ftp2fn_null(efNDX,NFILE,fnm),
	   &sys->symtab,&(sys->atoms),bVerbose,ir,&sys->idef,
	   forward);
  if (debug)
    pr_symtab(debug,0,"After index",&sys->symtab);
  triple_check(mdparin,ir,sys,&nerror);
  close_symtab(&sys->symtab);
  if (debug)
    pr_symtab(debug,0,"After close",&sys->symtab);

  if (ftp2bSet(efTRN,NFILE,fnm)) {
    if (bVerbose)
      fprintf(stderr,"getting data from old trajectory ...\n");
    cont_status(ftp2fn(efTRN,NFILE,fnm),bNeedVel,bGenVel,fr_time,ir,&natoms,
		&x,&v,box,sys);
  }
  
  if ((ir->coulombtype == eelPPPM) || (ir->coulombtype == eelPME)) {
    /* Calculate the optimal grid dimensions */
    max_spacing = calc_grid(box,opts->fourierspacing,
			    &(ir->nkx),&(ir->nky),&(ir->nkz),nprocs);
    if ((ir->coulombtype == eelPPPM) && (max_spacing > 0.1)) {
      set_warning_line(mdparin,-1);
      sprintf(warn_buf,"Grid spacing larger then 0.1 while using PPPM.");
      warning(NULL);
    }
  }
  
  /* This is also necessary for setting the multinr arrays */
  split_top(bVerbose,nprocs,sys);

  if (bVerbose) 
    fprintf(stderr,"writing run input file...\n");

  write_tpx(ftp2fn(efTPX,NFILE,fnm),
	    0,ir->init_t,ir->init_lambda,ir,box,
	    natoms,x,v,NULL,sys);
  
  print_warn_num();
  
  thanx(stdout);
  
  return 0;
}
