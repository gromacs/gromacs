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
 * Giant Rising Ordinary Mutants for A Clerical Setup
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
#include "convparm.h"
#include "fatal.h"
#include "rdgroup.h"
#include "gmxfio.h"
#include "trnio.h"
#include "tpxio.h"
#include "dum_parm.h"

void check_solvent(bool bVerbose,int nmol,t_molinfo msys[],
		   int Nsim,t_simsystem Sims[],t_inputrec *ir,char *SolventOpt)
{
  int i,wmol,nwt;
  
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
	if (ir->solvent_opt == -1) 
	  ir->solvent_opt=nwt;
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

static int *shuffle_xv(char *ndx,
		       int ntab,int *tab,int nmol,t_molinfo *mol,
		       int natoms,rvec *x,rvec *v,
		       int Nsim,t_simsystem Sims[])
{
  FILE *out;
  rvec *xbuf,*vbuf;
  int  *nindex,**index,xind,*done,*forward,*backward;
  int  i,j,j0,k,mi,nnat;
  
  /* Determine order in old x array! 
   * the index array holds for each molecule type
   * a pointer to the position at which the coordinates start
   */
  snew(index,nmol);
  snew(nindex,nmol);
  snew(done,nmol);
  
  /* COmpute the number of copies for each molecule type */
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

static int check_atom_names(char *fn1, char *fn2, t_atoms *at1, t_atoms *at2)
{
  int i,nmismatch;
#define MAXMISMATCH 20

  assert(at1->nr==at2->nr);
  
  nmismatch=0;
  for(i=0; i < at1->nr; i++)
    if (strcmp( *(at1->atomname[i]) , *(at2->atomname[i]) ) != 0) {
      if (nmismatch < MAXMISMATCH)
	fprintf(stderr,
		"Warning: atom names in %s and %s don't match (%s - %s)\n",
		fn1, fn2, *(at1->atomname[i]), *(at2->atomname[i]));
      else if (nmismatch == MAXMISMATCH)
	fprintf(stderr,"(more than %d non-matching atom names)\n",MAXMISMATCH);
      nmismatch++;
    }
  
  return nmismatch;
}

static int *new_status(char *topfile,char *topppfile,char *confin,
		       t_gromppopts *opts,t_inputrec *ir,
		       bool bGenVel,bool bVerbose,int *natoms,
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

  check_solvent(bVerbose,nrmols,molinfo,Nsim,Sims,ir,opts->SolventOpt);
  
  ntab = 0;
  tab  = NULL;
  if (nprocs > 1) {
    tab=mk_shuffle_tab(nrmols,molinfo,nprocs,&ntab,Nsim,Sims,bVerbose);
    if (debug)
      for(i=0; (i<ntab); i++)
	fprintf(debug,"Mol[%5d] = %s\n",i,*molinfo[tab[i]].name);
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
  if (*natoms != sys->atoms.nr) {
    fprintf(stderr,
	    "ERROR: number of coordinates in confin (%d) "
	    "does not match topology (%d)\n",*natoms,sys->atoms.nr);
    (*nerror)++;
  } else {
    /* make space for coordinates and velocities */
    snew(confat,1);
    init_t_atoms(confat,*natoms,FALSE);
    snew(*x,*natoms);
    snew(*v,*natoms);
    read_stx_conf(confin,opts->title,confat,*x,*v,box);
    /* We can not dispose of atoms, as that might fuck up the stack,
     * at least in gcc/linux it does.
     *
     * free_t_atoms(dumat);
     * sfree(dumat);
     */
    nmismatch=check_atom_names(topfile, confin, &(sys->atoms), confat);
    if (nmismatch) {
      sprintf(buf,"%d non-matching atom name%s\n",nmismatch,
	      (nmismatch == 1) ? "" : "s");
      warning(buf);
    }
    if (ntab > 0) {
      if (bVerbose)
	fprintf(stderr,"Shuffling coordinates...\n");
      forward=shuffle_xv("deshuf.ndx",ntab,tab,nrmols,molinfo,
			 *natoms,*x,*v,Nsim,Sims);
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

static void cont_status(char *slog,bool bNeedVel,bool bGenVel, real time,
			t_inputrec *ir,int *natoms,
			rvec **x,rvec **v,matrix box,
			int *nre,
			t_energy **e,
			t_topology *sys)
     /* If time == -1 read the last frame available which is complete */
{
  int         fp;
  real        tt;

  tt      = ir->init_t;
  fprintf(stderr,
	  "Reading Coordinates%s and Box size from old trajectory\n",
	  (!bNeedVel || bGenVel) ? "" : ", Velocities");
  if (time == -1)
    fprintf(stderr,"Will read whole trajectory\n");
  else
    fprintf(stderr,"Will read till time %g\n",time);
  if (!bNeedVel || bGenVel) {
    if (bGenVel)
      fprintf(stderr,"Velocities generated: "
	      "ignoring velocities in input trajectory\n");
    *natoms = read_first_x(&fp,slog,&tt,x,box);
  } else
    *natoms = read_first_x_v(&fp,slog,&tt,x,v,box);
  
  if(sys->atoms.nr != *natoms)
    fatal_error(0,"Number of atoms in Topology "
		"is not the same as in Trajectory");

  /* Now scan until the last set of box, x and v (step == 0)
   * or the ones at step step.
   * Or only until box and x if gen_vel is set.
   */
  while ((!bNeedVel || bGenVel) ?
	 read_next_x(fp,&tt,*natoms,*x,box):
	 read_next_x_v(fp,&tt,*natoms,*x,*v,box)) {
    if ( (time != -1) && (tt >= time) )
      break;
  }
  close_trj(fp);
  
  /*change the input record to the actual data*/
  ir->init_t = tt;
  
  fprintf(stderr,"Read frame from t = %g\n",tt);
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
  dumat.nr=natoms;
  snew(dumat.resname,natoms);
  snew(dumat.atom,natoms);
  snew(dumat.atomname,natoms);
  read_stx_conf(fn,title,&dumat,x,v,box);
  
  for(i=0; (i<pr->nr); i++) {
    ai=pr->param[i].AI;
    assert(ai < natoms);
    for(j=0; (j<DIM); j++)
      pr->param[i].c[j+DIM]=x[ai][j];
  }
  /*pr->nrfp+=DIM;*/
  
  sfree(x);
  sfree(v);
}

static int search_array(int *n,int map[],int key)
{
  int i;
  
  for(i=0; (i<*n); i++)
    if (map[i] == key)
      break;
  
  if (i == *n) {
    if (debug)
      fprintf(debug,"Renumbering atomtype %d to %d",key,*n);
    map[(*n)++]=key;
  }
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
      search_array(&nat,map,top->atoms.atom[i].type);
    top->atoms.atom[i].typeB=
      search_array(&nat,map,top->atoms.atom[i].typeB);
  }
  
  if (ir->solvent_opt != -1) {
    for(j=0; (j<nat); j++)
      if (map[j] == ir->solvent_opt)
	ir->solvent_opt=j;
    if (debug)
      fprintf(debug,"Renumbering solvent_opt (atomtype for OW) to %d\n",
	      ir->solvent_opt);
  }
    
  /* Renumber nlist */
  if (plist[F_LJ].nr > 0)
    ftype=F_LJ;
  else
    ftype=F_BHAM;
    
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
    "Using the [TT]-morse[tt] option grompp can convert the harmonic bonds",
    "in your topology to morse potentials. This makes it possible to break",
    "bonds. For this option to work you need an extra file in your $GMXLIB",
    "with dissociation energy. Use the -debug option to get more information",
    "on the workings of this option (look for MORSE in the grompp.log file",
    "using less or something like that).[PAR]",
    "To verify your run input file, please make notice of all warnings",
    "on the screen, and correct where necessary. Do also look at the contents",
    "of the [TT]mdout.mdp[tt] file, this contains comment lines, as well as",
    "the input that [TT]grompp[tt] has read. If in doubt you can start grompp",
    "with the [TT]-debug[tt] option which will give you more information",
    "in a file called grompp.log (along with real debug info). Finally, you",
    "can see the contents of the run input file with the [TT]gmxdump[tt]",
    "program."
  };
  static char *bugs[] = {
    "shuffling is sometimes buggy when used on systems when the number of "
    "molecules of a certain type is smaller than the number of processors."
  };
  t_gromppopts *opts;
  t_topology   sys;
  t_molinfo    msys;
  t_atomtype   atype;
  t_inputrec   *ir;
  t_energy     *e=NULL;
  int          nre=0,natoms;
  int          *forward=NULL;
  t_params     *plist;
  rvec         *x=NULL,*v=NULL;
  matrix       box;
  char         fn[STRLEN];
  int          nerror;
  bool         bNeedVel,bGenVel;

  t_filenm fnm[] = {
    { efMDP, NULL,  NULL,        ffREAD  },
    { efMDP, "-po", "mdout",     ffWRITE },
    { efSTX, "-c",  NULL,        ffREAD  },
    { efSTX, "-r",  NULL,        ffOPTRD },
    { efNDX, NULL,  NULL,        ffOPTRD },
    { efTOP, NULL,  NULL,        ffREAD  },
    { efTOP, "-pp", "processed", ffOPTWR },
    { efTPX, "-o",  NULL,        ffWRITE },
    { efTRN, "-t",  NULL,        ffOPTRD }
  };
#define NFILE asize(fnm)

  /* Command line options */
  static bool bVerbose=TRUE,bRenum=TRUE,bShuffle=FALSE;
  static int  nprocs=1,maxwarn=10;
  static real time=-1;
  t_pargs pa[] = {
    { "-np",      FALSE, etINT,  &nprocs,
      "Generate statusfile for # processors" },
    { "-time",    FALSE, etREAL, &time,
      "Take frame at or first after this time." },
    { "-v",       FALSE, etBOOL, &bVerbose,
      "Be loud and noisy" },
    { "-renum",   FALSE, etBOOL, &bRenum,
      "HIDDENRenumber atomtypes and minimize number of atomtypes" },
    { "-shuffle", FALSE, etBOOL, &bShuffle,
      "Shuffle molecules over processors (only with N > 1)" },
    { "-maxwarn", FALSE, etINT,  &maxwarn,
      "Number of warnings after which input processing stops" }
  };
  
  CopyRight(stdout,argv[0]);
  
  /* Initiate some variables */
  nerror=0;
  snew(ir,1);
  snew(opts,1);
  init_ir(ir,opts);
  
  /* Parse the command line */
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,asize(bugs),bugs);
  
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
  set_warning_line(ftp2fn(efMDP,NFILE,fnm),-1);    
  get_ir(ftp2fn(efMDP,NFILE,fnm),opt2fn("-po",NFILE,fnm),ir,opts,&nerror);

  if (bVerbose) 
    fprintf(stderr,"checking input for internal consistency...\n");
  check_ir(ir,opts,&nerror);

  if (bShuffle && (opts->eDisre==edrEnsemble)) {
    fprintf(stderr,"Can not shuffle and do ensemble averaging, "
	    "turning off shuffle\n");
    bShuffle=FALSE;
  }
  
  bNeedVel = (ir->eI == eiMD);
  bGenVel  = (bNeedVel && opts->bGenVel);

  snew(plist,F_NRE);
  init_plist(plist);
  open_symtab(&sys.symtab);
  strcpy(fn,ftp2fn(efTOP,NFILE,fnm));
  if (!fexist(fn)) 
    fatal_error(0,"%s does not exist",fn);
  forward=new_status(fn,opt2fn_null("-pp",NFILE,fnm),opt2fn("-c",NFILE,fnm),
		     opts,ir,bGenVel,bVerbose,&natoms,&x,&v,box,
		     &atype,&sys,&msys,plist,bShuffle ? nprocs : 1,
		     (opts->eDisre==edrEnsemble),opts->bMorse,&nerror);
  
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
  set_dummies(bVerbose, &sys.atoms, atype, plist, msys.plist);
  
  if (bRenum) 
    atype.nr=renum_atype(plist, &sys, atype.nr, ir, bVerbose);
  
  if (bVerbose) 
    fprintf(stderr,"converting bonded parameters...\n");
  convert_params(atype.nr, plist, msys.plist, &sys.idef);
  
  /* set ptype to Dummy for dummy atoms */
  set_dummies_ptype(bVerbose,&sys.idef,&sys.atoms);
  
  /* check masses */
  check_mol(&(sys.atoms));
  
  /* Now build the shakeblocks from the shakes */
  gen_sblocks(bVerbose,sys.atoms.nr,&(sys.idef),&(sys.blocks[ebSBLOCKS]));
   
  if (bVerbose) 
    fprintf(stderr,"initialising group options...\n");
  do_index(ftp2fn_null(efNDX,NFILE,fnm),
	   &sys.symtab,&(sys.atoms),bVerbose,ir,&sys.idef,
	   forward);
  triple_check(ir,&sys,&nerror);
  close_symtab(&sys.symtab);
  
  if (ftp2bSet(efTRN,NFILE,fnm)) {
    if (bVerbose)
      fprintf(stderr,"getting data from old trajectory ...\n");
    cont_status(ftp2fn(efTRN,NFILE,fnm),bNeedVel,bGenVel,time,ir,&natoms,
		&x,&v,box,&nre,&e,&sys);
  }
  /* This is also necessary for setting the multinr arrays */
  split_top(bVerbose,nprocs,&sys);

  if (nerror) {
    fprintf(stderr,"%d error%s, program terminated\n",
	    nerror,nerror==1 ? "" : "s");
    exit(1);
  }
  
  if (bVerbose) 
    fprintf(stderr,"writing run input file...\n");
  
  write_tpx(ftp2fn(efTPX,NFILE,fnm),
	    0,ir->init_t,ir->init_lambda,ir,box,
	    natoms,x,v,NULL,&sys);
  
  print_warn_num();
  
  thanx(stdout);
  
  return 0;
}
