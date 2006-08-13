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

#include "maths.h"
#include "macros.h"
#include "copyrite.h"
#include "bondf.h"
#include "string2.h"
#include "smalloc.h"
#include "strdb.h"
#include "sysstuff.h"
#include "confio.h"
#include "physics.h"
#include "statutil.h"
#include "vec.h"
#include "random.h"
#include "3dview.h"
#include "txtdump.h"
#include "readinp.h"
#include "names.h"
#include "toppush.h"
#include "pdb2top.h"
#include "gen_ad.h"
#include "topexcl.h"
#include "vec.h"
#include "x2top.h"
#include "atomprop.h"

char atp[6] = "HCNOSX";
#define NATP asize(atp)

real blen[NATP][NATP] = { 
  {  0.00,  0.108, 0.105, 0.10, 0.10, 0.10 },
  {  0.108, 0.15,  0.14,  0.14, 0.16, 0.14 },
  {  0.105, 0.14,  0.14,  0.14, 0.16, 0.14 },
  {  0.10,  0.14,  0.14,  0.14, 0.17, 0.14 },
  {  0.10,  0.16,  0.16,  0.17, 0.20, 0.17 },
  {  0.10,  0.14,  0.14,  0.14, 0.17, 0.17 }
};

#define MARGIN_FAC 1.1

static real set_x_blen(real scale)
{
  real maxlen;
  int  i,j;

  for(i=0; i<NATP-1; i++) {
    blen[NATP-1][i] *= scale;
    blen[i][NATP-1] *= scale;
  }
  blen[NATP-1][NATP-1] *= scale;
  
  maxlen = 0;
  for(i=0; i<NATP; i++)
    for(j=0; j<NATP; j++)
      if (blen[i][j] > maxlen)
	maxlen = blen[i][j];
  
  return maxlen*MARGIN_FAC;
}

static bool is_bond(int aai,int aaj,real len2)
{
  bool bIsBond;
  real bl;
    
  if (len2 == 0.0)
    bIsBond = FALSE;
  else {
    /* There is a bond when the deviation from ideal length is less than
     * 10 %
     */
    bl = blen[aai][aaj]*MARGIN_FAC;
    
    bIsBond = (len2 < bl*bl);
  }
#ifdef DEBUG
  if (debug)
    fprintf(debug,"ai: %5s  aj: %5s  len: %8.3f  bond: %s\n",
	    ai,aj,len,BOOL(bIsBond));
#endif
  return bIsBond;
}

static int get_atype(char *nm)
{
  int i,aai=NATP-1;
  
  for(i=0; (i<NATP-1); i++) {
    if (nm[0] == atp[i]) {
      aai=i;
      break;
    }
  }
  return aai;
}

void mk_bonds(t_atoms *atoms,rvec x[],t_params *bond,int nbond[],char *ff,
	      real cutoff,bool bPBC,matrix box,void *atomprop)
{
  t_param b;
  int     i,j,*elem;
  t_pbc   pbc;
  rvec    dx;
  real    dx2,c2;
  
  for(i=0; (i<MAXATOMLIST); i++)
    b.a[i] = -1;
  for(i=0; (i<MAXFORCEPARAM); i++)
    b.c[i] = 0.0;
    
  snew(elem,atoms->nr);
  for(i=0; (i<atoms->nr); i++) 
    elem[i] = get_atype(*atoms->atomname[i]);
  
  c2 = sqr(cutoff);
  if (bPBC)
    set_pbc(&pbc,box);
  for(i=0; (i<atoms->nr); i++) {
    if ((i % 10) == 0)
      fprintf(stderr,"\ratom %d",i);
    for(j=i+1; (j<atoms->nr); j++) {
      if (bPBC)
	pbc_dx(&pbc,x[i],x[j],dx);
      else
	rvec_sub(x[i],x[j],dx);
      
      dx2 = iprod(dx,dx);
      if ((dx2 < c2) && (is_bond(elem[i],elem[j],dx2))) {
	b.AI = i;
	b.AJ = j;
	b.C0 = sqrt(dx2);
	push_bondnow (bond,&b);
	nbond[i]++;
	nbond[j]++;
	if (debug) 
	  fprintf(debug,"Bonding atoms %s-%d and %s-%d\n",
		  *atoms->atomname[i],i+1,*atoms->atomname[j],j+1);
      }
    }
  }
  fprintf(stderr,"\ratom %d\n",i);
}

int *set_cgnr(t_atoms *atoms,bool bUsePDBcharge,real *qtot,real *mtot)
{
  int    i,n=1;
  int    *cgnr;
  double qt=0,mt=0;
  
  *qtot = *mtot = 0;
  snew(cgnr,atoms->nr);
  for(i=0; (i<atoms->nr); i++) {
    if (atoms->pdbinfo && bUsePDBcharge)
      atoms->atom[i].q = atoms->pdbinfo[i].bfac;
    qt += atoms->atom[i].q;
    *qtot += atoms->atom[i].q;
    *mtot += atoms->atom[i].m;
    cgnr[i] = n;
    if (is_int(qt)) {
      n++;
      qt=0;
    }
  }
  return cgnr;
}

t_atomtype *set_atom_type(t_symtab *tab,t_atoms *atoms,t_params *bonds,
			  int *nbonds,int nnm,t_nm2type nm2t[])
{
  t_atomtype *atype;
  int nresolved;
  
  snew(atype,1);
  snew(atoms->atomtype,atoms->nr);
  nresolved = nm2type(nnm,nm2t,tab,atoms,atype,nbonds,bonds);
  if (nresolved != atoms->nr)
    gmx_fatal(FARGS,"Could only find a forcefield type for %d out of %d atoms",
	      nresolved,atoms->nr);
  
  fprintf(stderr,"There are %d different atom types in your sample\n",
	  atype->nr);
    
  return atype;
}

void lo_set_force_const(t_params *plist,real c[],int nrfp,bool bRound,
			bool bDih,bool bParam)
{
  int    i,j;
  double cc;
  char   buf[32];
  
  for(i=0; (i<plist->nr); i++) {
    if (!bParam)
      for(j=0; j<nrfp; j++)
	c[j] = NOTSET;
    else {
      if (bRound) {
	sprintf(buf,"%.2e",plist->param[i].c[0]);
	sscanf(buf,"%lf",&cc);
	c[0] = cc;
      }
      else 
	c[0] = plist->param[i].c[0];
      if (bDih) {
	c[0] *= c[2];
	c[0] = ((int)(c[0] + 3600)) % 360;
	if (c[0] > 180)
	  c[0] -= 360;
	/* To put the minimum at the current angle rather than the maximum */
	c[0] += 180; 
      }
    }
    for(j=0; (j<nrfp); j++) {
      plist->param[i].c[j]      = c[j];
      plist->param[i].c[nrfp+j] = c[j];
    }
    set_p_string(&(plist->param[i]),"");
  }
}

void set_force_const(t_params plist[],real kb,real kt,real kp,bool bRound,
		     bool bParam)
{
  int i;
  real c[MAXFORCEPARAM];
  
  c[0] = 0;
  c[1] = kb;
  lo_set_force_const(&plist[F_BONDS],c,2,bRound,FALSE,bParam);
  c[1] = kt;
  lo_set_force_const(&plist[F_ANGLES],c,2,bRound,FALSE,bParam);
  c[1] = kp;
  c[2] = 3;
  lo_set_force_const(&plist[F_PDIHS],c,3,bRound,TRUE,bParam);
}

void calc_angles_dihs(t_params *ang,t_params *dih,rvec x[],bool bPBC,
		      matrix box)
{
  int    i,ai,aj,ak,al,t1,t2,t3;
  rvec   r_ij,r_kj,r_kl,m,n;
  real   sign,th,costh,ph,cosph;
  t_pbc  pbc;

  if (bPBC)
    set_pbc(&pbc,box);
  if (debug)
    pr_rvecs(debug,0,"X2TOP",box,DIM);
  for(i=0; (i<ang->nr); i++) {
    ai = ang->param[i].AI;
    aj = ang->param[i].AJ;
    ak = ang->param[i].AK;
    th = RAD2DEG*bond_angle(x[ai],x[aj],x[ak],bPBC ? &pbc : NULL,
			    r_ij,r_kj,&costh,&t1,&t2);
    if (debug)
      fprintf(debug,"X2TOP: ai=%3d aj=%3d ak=%3d r_ij=%8.3f r_kj=%8.3f th=%8.3f\n",
	      ai,aj,ak,norm(r_ij),norm(r_kj),th);
    ang->param[i].C0 = th;
  }
  for(i=0; (i<dih->nr); i++) {
    ai = dih->param[i].AI;
    aj = dih->param[i].AJ;
    ak = dih->param[i].AK;
    al = dih->param[i].AL;
    ph = RAD2DEG*dih_angle(x[ai],x[aj],x[ak],x[al],bPBC ? & pbc : NULL,
			   r_ij,r_kj,r_kl,m,n,&cosph,&sign,&t1,&t2,&t3);
    if (debug)
      fprintf(debug,"X2TOP: ai=%3d aj=%3d ak=%3d al=%3d r_ij=%8.3f r_kj=%8.3f r_kl=%8.3f ph=%8.3f\n",
	      ai,aj,ak,al,norm(r_ij),norm(r_kj),norm(r_kl),ph);
    dih->param[i].C0 = ph;
  }
}

static void dump_hybridization(FILE *fp,t_atoms *atoms,int nbonds[])
{
  int i;
  
  for(i=0; (i<atoms->nr); i++) {
    fprintf(fp,"Atom %5s has %1d bonds\n",*atoms->atomname[i],nbonds[i]);
  }
}

static void print_pl(FILE *fp,t_params plist[],int ftp,char *name,
		     char ***atomname)
{ 
  int i,j,nral,nrfp;

  if (plist[ftp].nr > 0) {
    fprintf(fp,"\n");
    fprintf(fp,"[ %s ]\n",name);
    nral = interaction_function[ftp].nratoms;
    nrfp = interaction_function[ftp].nrfpA;
    for(i=0; (i<plist[ftp].nr); i++) {
      for(j=0; (j<nral); j++) 
	fprintf(fp,"  %5s",*atomname[plist[ftp].param[i].a[j]]);
      for(j=0; (j<nrfp); j++) 
	fprintf(fp,"  %10.3e",plist[ftp].param[i].c[j]);
      fprintf(fp,"\n");
    }
  }
}

static void print_rtp(char *filenm,char *title,t_atoms *atoms,
		      t_params plist[],int cgnr[],int nbts,int bts[])
{
  FILE *fp;
  int i;
  
  fp = ffopen(filenm,"w");
  fprintf(fp,"; %s\n",title);
  fprintf(fp,"\n");
  fprintf(fp,"[ %s ]\n",*atoms->resname[0]);
  fprintf(fp,"\n");
  fprintf(fp,"[ atoms ]\n");
  for(i=0; (i<atoms->nr); i++) {
    fprintf(fp,"%-8s  %12s  %8.4f  %5d\n",
	    *atoms->atomname[i],*atoms->atomtype[i],
	    atoms->atom[i].q,cgnr[i]);
  }
  print_pl(fp,plist,F_BONDS,"bonds",atoms->atomname);
  print_pl(fp,plist,F_ANGLES,"angles",atoms->atomname);
  print_pl(fp,plist,F_PDIHS,"dihedrals",atoms->atomname);
  print_pl(fp,plist,F_IDIHS,"impropers",atoms->atomname);
  
  fclose(fp);
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "x2top generates a primitive topology from a coordinate file.",
    "The program assumes all hydrogens are present when defining",
    "the hybridization from the atom name and the number of bonds.",
    "The program can also make an rtp entry, which you can then add",
    "to the rtp database.[PAR]",
    "When [TT]-param[tt] is set, equilibrium distances and angles",
    "and force constants will be printed in the topology for all",
    "interactions. The equilibrium distances and angles are taken",
    "from the input coordinates, the force constant are set with",
    "command line options."
    "The force fields supported currently are:[PAR]",
    "G43a1  GROMOS96 43a1 Forcefield (official distribution)[PAR]",
    "oplsaa OPLS-AA/L all-atom force field (2001 aminoacid dihedrals)[PAR]",
    "G43b1  GROMOS96 43b1 Vacuum Forcefield (official distribution)[PAR]",
    "gmx    Gromacs Forcefield (a modified GROMOS87, see manual)[PAR]",
    "G43a2  GROMOS96 43a2 Forcefield (development) (improved alkane dihedrals)[PAR]",
    "The corresponding data files can be found in the library directory",
    "with names like ffXXXX.YYY. Check chapter 5 of the manual for more",
    "information about file formats. By default the forcefield selection",
    "is interactive, but you can use the [TT]-ff[tt] option to specify",
    "one of the short names above on the command line instead. In that",
    "case pdb2gmx just looks for the corresponding file.[PAR]",
  };
  static char *bugs[] = {
    "The atom type selection is primitive. Virtually no chemical knowledge is used",
    "Periodic boundary conditions screw up the bonding",
    "No improper dihedrals are generated",
    "The atoms to atomtype translation table is incomplete (ffG43a1.n2t file in the $GMXLIB directory). Please extend it and send the results back to the GROMACS crew."
  };
  FILE       *fp;
  t_params   plist[F_NRE];
  t_excls    *excls;
  t_atoms    *atoms;       /* list with all atoms */
  t_atomtype *atype;
  t_nextnb   nnb;
  t_nm2type  *nm2t;
  t_mols     mymol;
  void       *atomprop;
  int        nnm;
  char       title[STRLEN],forcefield[32];
  rvec       *x;        /* coordinates? */
  int        *nbonds,*cgnr;
  int        bts[] = { 1,1,1,2 };
  matrix     box;          /* box length matrix */
  int        natoms;       /* number of atoms in one molecule  */
  int        nres;         /* number of molecules? */
  int        i,j,k,l,m,ndih;
  bool       bRTP,bTOP,bOPLS;
  t_symtab   symtab;
  real       cutoff,qtot,mtot;
  
  t_filenm fnm[] = {
    { efSTX, "-f", "conf", ffREAD  },
    { efTOP, "-o", "out",  ffOPTWR },
    { efRTP, "-r", "out",  ffOPTWR }
  };
#define NFILE asize(fnm)
  static real scale = 1.1, kb = 4e5,kt = 400,kp = 5;
  static int  nexcl = 3;
  static bool bRemoveDih = FALSE;
  static bool bParam = TRUE, bH14 = TRUE,bAllDih = FALSE,bRound = TRUE;
  static bool bPairs = TRUE, bPBC = TRUE;
  static bool bUsePDBcharge = FALSE;
  static char *molnm = "ICE";
  static char *ff = "select";
  t_pargs pa[] = {
    { "-scale", FALSE, etREAL, {&scale},
      "Scaling factor for bonds with unknown atom types relative to atom type O" },
    { "-ff",     FALSE, etSTR, {&ff},
      "Select the force field for your simulation." },
    { "-nexcl", FALSE, etINT,  {&nexcl},
      "Number of exclusions" },
    { "-H14",    FALSE, etBOOL, {&bH14}, 
      "Use 3rd neighbour interactions for hydrogen atoms" },
    { "-alldih", FALSE, etBOOL, {&bAllDih}, 
      "Generate all proper dihedrals" },
    { "-remdih", FALSE, etBOOL, {&bRemoveDih}, 
      "Remove dihedrals on the same bond as an improper" },
    { "-pairs",  FALSE, etBOOL, {&bPairs},
      "Output 1-4 interactions (pairs) in topology file" },
    { "-name",   FALSE, etSTR,  {&molnm},
      "Name of your molecule" },
    { "-pbc",    FALSE, etBOOL, {&bPBC},
      "Use periodic boundary conditions." },
    { "-pdbq",  FALSE, etBOOL, {&bUsePDBcharge},
      "Use the B-factor supplied in a pdb file for the atomic charges" },
    { "-param", FALSE, etBOOL, {&bParam},
      "Print parameters in the output" },
    { "-round",  FALSE, etBOOL, {&bRound},
      "Round off measured values" },
    { "-kb",    FALSE, etREAL, {&kb},
      "Bonded force constant (kJ/mol/nm^2)" },
    { "-kt",    FALSE, etREAL, {&kt},
      "Angle force constant (kJ/mol/rad^2)" },
    { "-kp",    FALSE, etREAL, {&kp},
      "Dihedral angle force constant (kJ/mol/rad^2)" }
  };
  
  CopyRight(stdout,argv[0]);

  parse_common_args(&argc,argv,0,NFILE,fnm,asize(pa),pa,
		    asize(desc),desc,asize(bugs),bugs);
  bRTP = opt2bSet("-r",NFILE,fnm);
  bTOP = opt2bSet("-o",NFILE,fnm);
  
  if (!bRTP && !bTOP)
    gmx_fatal(FARGS,"Specify at least one output file");

  atomprop = get_atomprop();
    
  cutoff = set_x_blen(scale);

  if(!strncmp(ff,"select",6)) {
    /* Interactive forcefield selection */
    choose_ff(forcefield,sizeof(forcefield));
  } else {
    sprintf(forcefield,"ff%s",ff);
  }
  {
    char rtp[STRLEN];

    sprintf(rtp,"%s.rtp",forcefield);
    printf("Looking whether force field file %s exists\n",rtp);
    fclose(libopen(rtp));
  }

  bOPLS = (strcmp(forcefield,"ffoplsaa") == 0);
  
    
  mymol.name = strdup(molnm);
  mymol.nr   = 1;
	
  /* Init parameter lists */
  init_plist(plist);
  
  /* Read coordinates */
  get_stx_coordnum(opt2fn("-f",NFILE,fnm),&natoms); 
  snew(atoms,1);
  
  /* make space for all the atoms */
  init_t_atoms(atoms,natoms,TRUE);
  snew(x,natoms);              

  read_stx_conf(opt2fn("-f",NFILE,fnm),title,atoms,x,NULL,box);

  snew(nbonds,atoms->nr);
  
  printf("Generating bonds from distances...\n");
  mk_bonds(atoms,x,&(plist[F_BONDS]),nbonds,forcefield,cutoff,
	   bPBC,box,atomprop);

  nm2t = rd_nm2type(forcefield,&nnm);
  printf("There are %d name to type translations\n",nnm);
  if (debug)
    dump_nm2type(debug,nnm,nm2t);
  
  open_symtab(&symtab);
  atype = set_atom_type(&symtab,atoms,&(plist[F_BONDS]),nbonds,nnm,nm2t);
  
  /* Make Angles and Dihedrals */
  snew(excls,atoms->nr);
  printf("Generating angles and dihedrals from bonds...\n");
  init_nnb(&nnb,atoms->nr,4);
  gen_nnb(&nnb,plist);
  print_nnb(&nnb,"NNB");
  gen_pad(&nnb,atoms,bH14,nexcl,plist,excls,NULL,bAllDih,bRemoveDih,TRUE);
  done_nnb(&nnb);

  if (!bPairs)
    plist[F_LJ14].nr = 0;
  fprintf(stderr,
	  "There are %4d %s dihedrals, %4d impropers, %4d angles\n"
	  "          %4d pairs,     %4d bonds and  %4d atoms\n",
	  plist[F_PDIHS].nr, 
	  bOPLS ? "Ryckaert-Bellemans" : "proper",
	  plist[F_IDIHS].nr, plist[F_ANGLES].nr,
	  plist[F_LJ14].nr, plist[F_BONDS].nr,atoms->nr);

  calc_angles_dihs(&plist[F_ANGLES],&plist[F_PDIHS],x,bPBC,box);
  
  set_force_const(plist,kb,kt,kp,bRound,bParam);

  cgnr = set_cgnr(atoms,bUsePDBcharge,&qtot,&mtot);
  printf("Total charge is %g, total mass is %g\n",qtot,mtot);
  
  if (bTOP) {    
    fp = ftp2FILE(efTOP,NFILE,fnm,"w");
    print_top_header(fp,ftp2fn(efTOP,NFILE,fnm),
		     "Generated by x2top",TRUE, forcefield,1.0);
    
    write_top(fp,NULL,mymol.name,atoms,bts,plist,excls,atype,cgnr,nexcl);
    print_top_mols(fp,mymol.name,NULL,0,NULL,1,&mymol);
    
    fclose(fp);
  }
  if (bRTP)
    print_rtp(ftp2fn(efRTP,NFILE,fnm),"Generated by x2top",
	      atoms,plist,cgnr,asize(bts),bts);
  
  if (debug) {
    dump_hybridization(debug,atoms,nbonds);
  }
  close_symtab(&symtab);
    
  thanx(stderr);
  
  return 0;
}
