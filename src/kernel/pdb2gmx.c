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
static char *SRCID_pdb2gmx_c = "$Id$";

#include <time.h>
#include <ctype.h>
#include "assert.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "copyrite.h"
#include "string2.h"
#include "confio.h"
#include "symtab.h"
#include "vec.h"
#include "statutil.h"
#include "futil.h"
#include "fatal.h"
#include "pdbio.h"
#include "toputil.h"
#include "h_db.h"
#include "physics.h"
#include "pgutil.h"
#include "calch.h"
#include "resall.h"
#include "pdb2top.h"
#include "ter_db.h"
#include "strdb.h"
#include "gbutil.h"
#include "genhydro.h"
#include "readinp.h"
#include "xlate.h"
#include "specbond.h"
#include "index.h"
#include "hizzie.h"

#define NREXCL 3

static char *select_res(int nr,int resnr,char *name[],char *expl[],char *title)
{
  int sel=0;

  printf("Which %s type do you want for residue %d\n",title,resnr+1);
  for(sel=0; (sel < nr); sel++)
    printf("%d. %s (%s)\n",sel,expl[sel],name[sel]);
  printf("\nType a number:"); fflush(stdout);

  if (scanf("%d",&sel) != 1)
    fatal_error(0,"Answer me for res %s %d!",title,resnr+1);
  
  return name[sel];
}

static char *get_asptp(int resnr)
{
  enum { easp, easpH, easpNR };
  static char *lh[easpNR] = { "ASP", "ASPH" };
  static char *expl[easpNR] = {
    "Not protonated (charge -1)",
    "Protonated (charge 0)"
  };

  return select_res(easpNR,resnr,lh,expl,"ASPARTIC ACID");
}

static char *get_glutp(int resnr)
{
  enum { eglu, egluH, egluNR };
  static char *lh[egluNR] = { "GLU", "GLUH" };
  static char *expl[egluNR] = {
    "Not protonated (charge -1)",
    "Protonated (charge 0)"
  };

  return select_res(egluNR,resnr,lh,expl,"GLUTAMIC ACID");
}

static char *get_lystp(int resnr)
{
  enum { elys, elysH, elysNR };
  static char *lh[elysNR] = { "LYS", "LYSH" };
  static char *expl[elysNR] = {
    "Not protonated (charge 0)",
    "Protonated (charge +1)"
  };

  return select_res(elysNR,resnr,lh,expl,"LYSINE");
}

static char *get_cystp(int resnr)
{
  enum { ecys, ecysH, ecysNR };
  static char *lh[ecysNR] = { "CYS", "CYSH" };
  static char *expl[ecysNR] = {
    "Cysteine in disulfide bridge",
    "Protonated"
  };

  return select_res(ecysNR,resnr,lh,expl,"CYSTEINE");

}

static char *get_histp(int resnr)
{
  static char *expl[ehisNR] = {
    "H on ND1 only",
    "H on NE2 only",
    "H on ND1 and NE2",
    "Coupled to Heme"
  };
  
  return select_res(ehisNR,resnr,hh,expl,"HISTIDINE");
}

static void rename_pdbres(t_atoms *pdba,char *oldnm,char *newnm,
			  bool bFullCompare)
{
  char *resnm;
  int i;
  
  for(i=0; (i<pdba->nres); i++) {
    resnm=*pdba->resname[i];
    if ((bFullCompare && (strcasecmp(resnm,oldnm) == 0)) ||
	(!bFullCompare && strstr(resnm,oldnm) != NULL)) {
      sfree(*pdba->resname[i]);
      *pdba->resname[i]=strdup(newnm);
    }
  }
}

static void rename_pdbresint(t_atoms *pdba,char *oldnm,
			     char *gettp(int),bool bFullCompare)
{
  int  i;
  char *ptr,*resnm;
  
  for(i=0; i<pdba->nres; i++) {
    resnm=*pdba->resname[i];
    if ((bFullCompare && (strcmp(resnm,oldnm) == 0)) ||
	(!bFullCompare && strstr(resnm,oldnm) != NULL)) {
      ptr=gettp(i);
      sfree(*pdba->resname[i]);
      *pdba->resname[i]=strdup(ptr);
    }
  }
}

static void check_occupancy(t_atoms *atoms,char *filename)
{
  int i,ftp;
  int nzero=0;
  int nnotone=0;
  
  ftp = fn2ftp(filename);
  if (!atoms->pdbinfo || ((ftp != efPDB) && (ftp != efBRK) && (ftp != efENT)))
    fprintf(stderr,"No occupancies in %s\n",filename);
  else {
    for(i=0; (i<atoms->nr); i++) {
      if (atoms->pdbinfo[i].occup == 0)
	nzero++;
      else if (atoms->pdbinfo[i].occup != 1)
	nnotone++;
    }
    if (nzero == atoms->nr)
      fprintf(stderr,"All occupancy fields zero. This is probably not an X-Ray structure\n");
    else if ((nzero > 0) || (nnotone > 0))
      fprintf(stderr,
	      "WARNING: there were %d atoms with zero occupancy and %d atoms"
	      " with\n         occupancy unequal to one (out of %d atoms)."
	      " Check your pdb file.\n",nzero,nnotone,atoms->nr);
    else
      fprintf(stderr,"All occupancies are one\n");
  }
}

void write_posres(char *fn,t_atoms *pdba)
{
  FILE *fp;
  int  i;
  
  fp=ffopen(fn,"w");
  fprintf(fp,
	  "; In this topology include file, you will find position restraint\n"
	  "; entries for all the heavy atoms in your original pdb file.\n"
	  "; This means that all the protons which were added by pdb2gmx are\n"
	  "; not restrained.\n"
	  "\n"
	  "[ position_restraints ]\n"
	  "; %4s%6s%8s%8s%8s\n","atom","type","fx","fy","fz"
	  );
  for(i=0; (i<pdba->nr); i++) {
    if (!is_hydrogen(*pdba->atomname[i]) && !is_dummymass(*pdba->atomname[i]))
      fprintf(fp,"%6d%6d%8.1f%8.1f%8.1f\n",i+1,1,1000.0,1000.0,1000.0);
  }
  ffclose(fp);
}

int read_pdball(char *inf, char *outf,char *title,
		t_atoms *atoms, rvec **x,matrix box, bool bRetainH)
/* Read a pdb file. (containing proteins) */
{
  int       natom,new_natom,i;
  
  /* READ IT */
  printf("Reading %s...\n",inf);
  get_stx_coordnum(inf,&natom);
  init_t_atoms(atoms,natom,TRUE);
  snew(*x,natom);
  read_stx_conf(inf,title,atoms,*x,NULL,box);
  if (!bRetainH) {
    new_natom=0;
    for(i=0; i<atoms->nr; i++)
      if (!is_hydrogen(*atoms->atomname[i])) {
	atoms->atom[new_natom]=atoms->atom[i];
	atoms->atomname[new_natom]=atoms->atomname[i];
	atoms->pdbinfo[new_natom]=atoms->pdbinfo[i];
	copy_rvec((*x)[i],(*x)[new_natom]);
	new_natom++;
      }
    atoms->nr=new_natom;
    natom=new_natom;
  }
    
  printf("Read");
  if (title && title[0])
    printf(" '%s',",title);
  printf(" %d atoms\n",natom);
  
  /* Rename residues */
  rename_pdbres(atoms,"SOL","HOH",FALSE);
  rename_pdbres(atoms,"WAT","HOH",FALSE);
  rename_pdbres(atoms,"HEM","HEME",FALSE);

  rename_atoms(atoms);
  
  if (natom == 0)
    return 0;

  if (outf)
    write_sto_conf(outf,title,atoms,*x,NULL,box);
 
  return natom;
}

void process_chain(t_atoms *pdba, rvec *x, 
		   bool bTrpU,bool bPheU,bool bTyrU,
		   bool bLysMan,bool bAspMan,bool bGluMan,
		   bool bHisMan,bool bCysMan,
		   int *nssbonds,t_ssbond **ssbonds,
		   real angle,real distance)
{
  /* Rename aromatics, lys, asp and histidine */
  if (bTyrU) rename_pdbres(pdba,"TYR","TYRU",FALSE);
  if (bTrpU) rename_pdbres(pdba,"TRP","TRPU",FALSE);
  if (bPheU) rename_pdbres(pdba,"PHE","PHEU",FALSE);
  if (bLysMan) 
    rename_pdbresint(pdba,"LYS",get_lystp,FALSE);
  else
    rename_pdbres(pdba,"LYS","LYSH",FALSE);
  if (bAspMan) 
    rename_pdbresint(pdba,"ASP",get_asptp,FALSE);
  else
    rename_pdbres(pdba,"ASPH","ASP",FALSE);
  if (bGluMan) 
    rename_pdbresint(pdba,"GLU",get_glutp,FALSE);
  else
    rename_pdbres(pdba,"GLUH","GLU",FALSE);

  /* Make sure we don't have things like CYS? */ 
  rename_pdbres(pdba,"CYS","CYS",FALSE);
  *nssbonds=mk_specbonds(pdba,x,bCysMan,ssbonds);
  rename_pdbres(pdba,"CYS","CYSH",TRUE);

  if (!bHisMan)
    set_histp(pdba,x,angle,distance);
  else
    rename_pdbresint(pdba,"HIS",get_histp,TRUE);
}

/* struct for sorting the atoms from the pdb file */
typedef struct {
  int  resnr;  /* residue number               */
  int  j;      /* database order index         */
  int  index;  /* original atom number         */
  char anm1;   /* second letter of atom name   */
  char altloc; /* alternate location indicator */
} t_pdbindex;
  
int pdbicomp(const void *a,const void *b)
{
  t_pdbindex *pa,*pb;
  int d;

  pa=(t_pdbindex *)a;
  pb=(t_pdbindex *)b;

  d = (pa->resnr - pb->resnr);
  if (d==0) {
    d = (pa->j - pb->j);
    if (d==0) {
      d = (pa->anm1 - pb->anm1);
      if (d==0)
	d = (pa->altloc - pb->altloc);
    }
  }

  return d;
}

static void sort_pdbatoms(int nrtp,t_restp restp[],
			  int natoms,t_atoms **pdbaptr,rvec **x,
			  t_block *block,char ***gnames)
{
  t_atoms *pdba,*pdbnew;
  rvec **xnew;
  int     i,j;
  t_restp *rptr;
  t_pdbindex *pdbi;
  atom_id *a;
  char *atomnm,*resnm;
  
  pdba=*pdbaptr;
  natoms=pdba->nr;
  pdbnew=NULL;
  snew(xnew,1);
  snew(pdbi, natoms);
  
  for(i=0; i<natoms; i++) {
    atomnm=*pdba->atomname[i];
    resnm=*pdba->resname[pdba->atom[i].resnr];
    if ((rptr=search_rtp(resnm,nrtp,restp)) == NULL)
      fatal_error(0,"Residue type %s not found",resnm);
    for(j=0; (j<rptr->natom); j++)
      if (strcasecmp(atomnm,*(rptr->atomname[j])) == 0)
	break;
    if (j==rptr->natom) {
      if ( ( ( pdba->atom[i].resnr == 0) && (atomnm[0] == 'H') &&
	     ( (atomnm[1] == '1') || (atomnm[1] == '2') || 
	       (atomnm[1] == '3') ) ) )
	j=1;
      else {
	char buf[STRLEN];
	
	sprintf(buf,"Atom %s in residue %s %d not found in database\n"
		"             while sorting atoms%n",atomnm,
		rptr->resname,pdba->atom[i].resnr+1,&i);
	if ( is_hydrogen(atomnm) )
	  sprintf(buf+i,". Maybe different protonation state.\n"
		  "             Remove this hydrogen or choose a different "
		  "protonation state.");
	fatal_error(0,buf);
      }
    }
    /* make shadow array to be sorted into indexgroup */
    pdbi[i].resnr  = pdba->atom[i].resnr;
    pdbi[i].j      = j;
    pdbi[i].index  = i;
    pdbi[i].anm1   = atomnm[1];
    pdbi[i].altloc = pdba->pdbinfo[i].altloc;
  }
  qsort(pdbi,natoms,(size_t)sizeof(pdbi[0]),pdbicomp);
  
  /* pdba is sorted in pdbnew using the pdbi index */ 
  snew(a,natoms);
  snew(pdbnew,1);
  init_t_atoms(pdbnew,natoms,TRUE);
  snew(*xnew,natoms);
  pdbnew->nr=pdba->nr;
  pdbnew->nres=pdba->nres;
  sfree(pdbnew->resname);
  pdbnew->resname=pdba->resname;
  for (i=0; i<natoms; i++) {
    pdbnew->atom[i]     = pdba->atom[pdbi[i].index];
    pdbnew->atomname[i] = pdba->atomname[pdbi[i].index];
    pdbnew->pdbinfo[i]  = pdba->pdbinfo[pdbi[i].index];
    copy_rvec((*x)[pdbi[i].index],(*xnew)[i]);
     /* make indexgroup in block */
    a[i]=pdbi[i].index;
  }
  /* clean up */
  sfree(pdba->atomname);
  sfree(pdba->atom);
  sfree(pdba->pdbinfo);
  done_block(&pdba->excl);
  sfree(pdba);
  sfree(*x);
  /* copy the sorted pdbnew back to pdba */
  *pdbaptr=pdbnew;
  *x=*xnew;
  add_grp(block, gnames, natoms, a, "prot_sort");
  sfree(xnew);
  sfree(a);
  sfree(pdbi);
}

static int remove_duplicate_atoms(t_atoms *pdba,rvec x[])
{
  int     i,j,oldnatoms;
  
  printf("Checking for duplicate atoms....\n");
  oldnatoms    = pdba->nr;
  
  /* NOTE: pdba->nr is modified inside the loop */
  for(i=1; (i < pdba->nr); i++) {
    /* compare 'i' and 'i-1', throw away 'i' if they are identical 
       this is a 'while' because multiple alternate locations can be present */
    while ( (i < pdba->nr) &&
	    (pdba->atom[i-1].resnr == pdba->atom[i].resnr) &&
	    (strcmp(*pdba->atomname[i-1],*pdba->atomname[i])==0) ) {
      printf("deleting duplicate atom %4s  %s%4d",
	     *pdba->atomname[i], *pdba->resname[pdba->atom[i].resnr], 
	     pdba->atom[i].resnr+1);
      if (pdba->atom[i].chain && (pdba->atom[i].chain!=' '))
	printf(" ch %c", pdba->atom[i].chain);
      if (pdba->pdbinfo) {
	if (pdba->pdbinfo[i].atomnr)
	  printf("  pdb nr %4d",pdba->pdbinfo[i].atomnr);
	if (pdba->pdbinfo[i].altloc && (pdba->pdbinfo[i].altloc!=' '))
	  printf("  altloc %c",pdba->pdbinfo[i].altloc);
      }
      printf("\n");
      pdba->nr--;
      sfree(pdba->atomname[i]);
      for (j=i; j < pdba->nr; j++) {
	pdba->atom[j]     = pdba->atom[j+1];
	pdba->atomname[j] = pdba->atomname[j+1];
	pdba->pdbinfo[j]  = pdba->pdbinfo[j+1];
	copy_rvec(x[j+1],x[j]);
      }
      srenew(pdba->atom,     pdba->nr);
      srenew(pdba->atomname, pdba->nr);
      srenew(pdba->pdbinfo,  pdba->nr);
    }
  }
  if (pdba->nr != oldnatoms)
    printf("Now there are %d atoms\n",pdba->nr);
  
  return pdba->nr;
}

void find_nc_ter(t_atoms *pdba,int *rn,int *rc)
{
  int rnr;
  
  *rn=-1;
  *rc=-1;
  for(rnr=0; (rnr<pdba->nres); rnr++) {
    if ((*rn == -1) && (is_protein(*pdba->resname[rnr])))
	*rn=rnr;
    if ((*rc != rnr) && (is_protein(*pdba->resname[rnr])))
      *rc=rnr;
  }
  if (debug) fprintf(debug,"nres: %d, rN: %d, rC: %d\n",pdba->nres,*rn,*rc);
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "This program reads a pdb file, lets you choose a forcefield, reads",
    "some database files, adds hydrogens to the molecules and generates",
    "coordinates in Gromacs (Gromos) format and a topology in Gromacs format.",
    "These files can subsequently be processed to generate a run input file.",
    "[PAR]",
    
    "Note that a pdb file is nothing more than a file format, and it",
    "need not necessarily contain a protein structure. Every kind of",
    "molecule for which there is support in the database can be converted.",
    "If there is no support in the database, you can add it yourself.[PAR]",
    
    "The program has limited intelligence, it reads a number of database",
    "files, that allow it to make special bonds (Cys-Cys, Heme-His, etc.),",
    "if necessary this can be done manually. The program can prompt the",
    "user to select which kind of LYS, ASP, GLU, CYS or HIS residue she",
    "wants. For LYS the choice is between LYS (two protons on NZ) or LYSH",
    "(three protons, default), for ASP and GLU unprotonated (default) or",
    "protonated, for HIS the proton can be either on ND1 (HISA), on NE2",
    "(HISB) or on both (HISH). By default these selections are done",
    "automatically. For His, this is based on an optimal hydrogen bonding",
    "conformation. Hydrogen bonds are defined based on a simple geometric",
    "criterium, specified by the maximum hydrogen-donor-acceptor angle",
    "and donor-acceptor distance, which are set by [TT]-angle[tt] and",
    "[TT]-dist[tt] respectively.[PAR]",
    
    "pdb2gmx will also check the occupancy field of the pdb file.",
    "If any of the occupanccies are not one, indicating that the atom is",
    "not resolved well in the structure, a warning message is issued.",
    "When a pdb file does not originate from an X-Ray structure determination",
    "all occupancy fields may be zero. Either way, it is up to the user",
    "to verify the correctness of the input data (read the article!).[PAR]", 
    
    "During processing the atoms will be reordered according to Gromacs",
    "conventions. With [TT]-n[tt] an index file can be generated that",
    "contains one group reordered in the same way. This allows you to",
    "convert a Gromos trajectory and coordinate file to Gromos. There is",
    "one limitation: reordering is done after the hydrogens are stripped",
    "from the input and before new hydrogens are added. This means that",
    "you should not turn off [TT]-reth[tt].[PAR]",

    "The [TT].gro[tt] and [TT].g96[tt] file formats do not support chain",
    "identifiers. Therefore it is useful to enter a pdb file name at",
    "the [TT]-o[tt] option when you want to convert a multichain pdb file.",
    "[PAR]",
    
    "When using [TT]-reth[tt] to keep all hydrogens from the [TT].pdb[tt]",
    "file, the names of the hydrogens in the [TT].pdb[tt] file [IT]must[it]",
    "match the names in the database.[PAR]", 
    
    "[TT]-sort[tt] will sort all residues according to the order in the",
    "database, sometimes this is necessary to get charge groups",
    "together.[PAR]",
    
    "[TT]-alldih[tt] will generate all proper dihedrals instead of only",
    "those with as few hydrogens as possible, this is useful for use with",
    "the Charmm forcefield.[PAR]",
    
    "The option [TT]-dummy[tt] removes hydrogen and fast improper dihedral",
    "motions. Angular and out-of-plane motions can be removed by changing",
    "hydrogens into dummy atoms and fixing angles, which fixes their",
    "position relative to neighboring atoms. Additionally, all atoms in the",
    "aromatic rings of the standard amino acids (i.e. PHE, TRP, TYR and HIS)",
    "can be converted into dummy atoms, elminating the fast improper dihedral",
    "fluctuations in these rings. Note that in this case all other hydrogen",
    "atoms are also converted to dummy atoms. The mass of all atoms that are",
    "converted into dummy atoms, is added to the heavy atoms.[PAR]",
    "Also slowing down of dihedral motion can be done with [TT]-heavyh[tt]",
    "done by increasing the hydrogen-mass by a factor of 4. This is also",
    "done for water hydrogens to slow down the rotational motion of water.",
    "The increase in mass of the hydrogens is subtracted from the bonded",
    "(heavy) atom so that the total mass of the system remains the same."
  };

  typedef struct {
    char chain;
    int  start;
    int  natom;
    bool bAllWat;
  } t_pdbchain;

  typedef struct {
    char chain;
    bool bAllWat;
    t_atoms *pdba;
    rvec *x;
  } t_chain;
  
  FILE       *fp,*top_file,*top_file2,*itp_file=NULL;
  int        natom,nres;
  t_atoms    pdba_all,*pdba;
  t_atoms    *atoms;
  t_block    *block;
  int        chain,nchain,nwaterchain;
  t_pdbchain *pdbchains;
  t_chain    *chains;
  char       pchain;
  int        nincl,nmol;
  char       **incls;
  t_mols     *mols;
  char       **gnames;
  matrix     box;
  rvec       box_space;
  char       *ff;
  int        i,j,k,l,nrtp,rN,rC;
  int        *swap_index,si;
  int        bts[ebtsNR];
  t_restp    *restp;
  t_hackblock *ah;
  t_symtab   tab;
  t_atomtype *atype;
  char       fn[256],*top_fn,itp_fn[STRLEN],posre_fn[STRLEN];
  char       molname[STRLEN],title[STRLEN];
  char       *c;
  int        nah,nNtdb,nCtdb;
  t_hackblock *ntdb,*ctdb,*sel_ntdb=NULL,*sel_ctdb=NULL;
  int        nssbonds;
  t_ssbond   *ssbonds;
  rvec       *pdbx,*x;
  bool       bUsed,bDummies=FALSE,bWat,bPrevWat=FALSE,bITP,bDummyAromatics=FALSE;
  real       mHmult=0;
  
  t_filenm   fnm[] = { 
    { efSTX, "-f", "eiwit.pdb", ffREAD  },
    { efSTO, "-o", "conf",      ffWRITE },
    { efTOP, NULL, NULL,        ffWRITE },
    { efITP, "-i", "posre",     ffWRITE },
    { efNDX, "-n", "clean",     ffOPTWR },
    { efSTO, "-q", "clean.pdb", ffOPTWR }
  };
#define NFILE asize(fnm)
  
  /* Command line arguments msut be static */
  static bool bNewRTP=FALSE;
  static bool bInter=FALSE, bCysMan=FALSE; 
  static bool bLysMan=FALSE, bAspMan=FALSE, bGluMan=FALSE, bHisMan=FALSE;
  static bool bTerMan=FALSE, bUnA=FALSE, bHeavyH;
  static bool bH14=FALSE, bSort=TRUE, bRetainH=TRUE;
  static bool bAlldih=FALSE;
  static real angle=135.0, distance=0.3;
  static real long_bond_dist=0.25, short_bond_dist=0.05;
  static char *dumstr[] = { NULL, "none", "hydrogens", "aromatics", NULL };
  t_pargs pa[] = {
    { "-newrtp", FALSE, etBOOL, {&bNewRTP},
      "HIDDENWrite the residue database in new format to 'new.rtp'"},
    { "-lb",     FALSE, etREAL, {&long_bond_dist},
      "HIDDENLong bond warning distance" },
    { "-sb",     FALSE, etREAL, {&short_bond_dist},
      "HIDDENShort bond warning distance" },
    { "-inter",  FALSE, etBOOL, {&bInter},
      "Set the next 6 options to interactive"},
    { "-ss",     FALSE, etBOOL, {&bCysMan}, 
      "Interactive SS bridge selection" },
    { "-ter",    FALSE, etBOOL, {&bTerMan}, 
      "Interactive termini selection, iso charged" },
    { "-lys",    FALSE, etBOOL, {&bLysMan}, 
      "Interactive Lysine selection, iso charged" },
    { "-asp",    FALSE, etBOOL, {&bAspMan}, 
      "Interactive Aspartic Acid selection, iso charged" },
    { "-glu",    FALSE, etBOOL, {&bGluMan}, 
      "Interactive Glutamic Acid selection, iso charged" },
    { "-his",    FALSE, etBOOL, {&bHisMan},
      "Interactive Histidine selection, iso checking H-bonds" },
    { "-angle",  FALSE, etREAL, {&angle}, 
      "Minimum hydrogen-donor-acceptor angle for a H-bond (degrees)" },
    { "-dist",   FALSE, etREAL, {&distance},
      "Maximum donor-acceptor distance for a H-bond (nm)" },
    { "-una",    FALSE, etBOOL, {&bUnA}, 
      "Select aromatic rings with united CH atoms on Phenylalanine, "
      "Tryptophane and Tyrosine" },
    { "-sort",   FALSE, etBOOL, {&bSort}, 
      "Sort the residues according to database" },
    { "-H14",    FALSE, etBOOL, {&bH14}, 
      "Use 1-4 interactions for hydrogen atoms" },
    { "-reth",   FALSE, etBOOL, {&bRetainH}, 
      "Retain hydrogen atoms that are in the pdb file" },
    { "-alldih", FALSE, etBOOL, {&bAlldih}, 
      "Generate all proper dihedrals" },
    { "-dummy",  FALSE, etENUM, {dumstr}, 
      "Convert atoms to dummy atoms" },
    { "-heavyh", FALSE, etBOOL, {&bHeavyH},
      "Make hydrogen atoms heavy" }
  };
#define NPARGS asize(pa)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    0,NULL);
  if (bInter) {
    /* if anything changes here, also change description of -inter */
    bCysMan = TRUE;
    bTerMan = TRUE;
    bLysMan = TRUE;
    bAspMan = TRUE;
    bGluMan = TRUE;
    bHisMan = TRUE;
  }
  
  if (bHeavyH)
    mHmult=4.0;
  else
    mHmult=1.0;
  
  switch(dumstr[0][0]) {
  case 'n': /* none */
    bDummies=FALSE;
    bDummyAromatics=FALSE;
    break;
  case 'h': /* hydrogens */
    bDummies=TRUE;
    bDummyAromatics=FALSE;
    break;
  case 'a': /* aromatics */
    bDummies=TRUE;
    bDummyAromatics=TRUE;
    break;
  default:
    fatal_error(0,"DEATH HORROR in $s (%d): dumstr[0]='%s'",
		__FILE__,__LINE__,dumstr[0]);
  }/* end switch */
  
  clear_mat(box);
  natom=read_pdball(opt2fn("-f",NFILE,fnm),opt2fn_null("-q",NFILE,fnm),title,
		    &pdba_all,&pdbx,box,bRetainH);
  
  if (natom==0)
    fatal_error(0,"No atoms found in pdb file %s\n",opt2fn("-f",NFILE,fnm));

  printf("Analyzing pdb file\n");
  nchain=0;
  nwaterchain=0;
  /* keep the compiler happy */
  pchain='?';
  pdbchains=NULL;
  for (i=0; (i<natom); i++) {
    bWat = strcasecmp(*pdba_all.resname[pdba_all.atom[i].resnr],"HOH") == 0;
    if ((i==0) || (pdba_all.atom[i].chain!=pchain) || (bWat != bPrevWat)) {
      pchain=pdba_all.atom[i].chain;
      /* set natom for previous chain */
      if (nchain > 0)
	pdbchains[nchain-1].natom=i-pdbchains[nchain-1].start;
      if (bWat) {
	nwaterchain++;
	pdba_all.atom[i].chain='\0';
      }
      /* check if chain identifier was used before */
      for (j=0; (j<nchain); j++)
	if ((pdbchains[j].chain != '\0') && (pdbchains[j].chain != ' ') &&
	    (pdbchains[j].chain == pdba_all.atom[i].chain))
	  fatal_error(0,"Chain identifier '%c' was used "
		      "in two non-sequential blocks (residue %d, atom %d)",
		      pdba_all.atom[i].chain,pdba_all.atom[i].resnr+1,i+1);
      srenew(pdbchains,nchain+1);
      pdbchains[nchain].chain=pdba_all.atom[i].chain;
      pdbchains[nchain].start=i;
      pdbchains[nchain].bAllWat=bWat;
      nchain++;
    }
    bPrevWat=bWat;
  }
  pdbchains[nchain-1].natom=natom-pdbchains[nchain-1].start;
  
  /* set all the water blocks at the end of the chain */
  snew(swap_index,nchain);
  j=0;
  for(i=0; i<nchain; i++)
    if (!pdbchains[i].bAllWat) {
      swap_index[j]=i;
      j++;
    }
  for(i=0; i<nchain; i++)
    if (pdbchains[i].bAllWat) {
      swap_index[j]=i;
      j++;
    }
  if (nwaterchain>1)
    printf("Moved all the water blocks to the end\n");

  snew(chains,nchain);
  /* copy pdb data and x for all chains */
  for (i=0; (i<nchain); i++) {
    si=swap_index[i];
    chains[i].chain   = pdbchains[si].chain;
    chains[i].bAllWat = pdbchains[si].bAllWat; 
    /* check for empty chain identifiers */
    if ((nchain-nwaterchain>1) && !pdbchains[si].bAllWat && 
	((chains[i].chain=='\0') || (chains[i].chain==' '))) {
      bUsed=TRUE;
      for(k='A'; (k<='Z') && bUsed; k++) {
	bUsed=FALSE;
	for(j=0; j<nchain; j++)
	  bUsed=bUsed || (!pdbchains[j].bAllWat && (chains[j].chain==k));
	if (!bUsed) {
	  printf("Gave chain %d chain identifier '%c'\n",i+1,k);
	  chains[i].chain=k;
	}
      }
    }
    snew(chains[i].pdba,1);
    init_t_atoms(chains[i].pdba,pdbchains[si].natom,TRUE);
    snew(chains[i].x,chains[i].pdba->nr);
    for (j=0; j<chains[i].pdba->nr; j++) {
      chains[i].pdba->atom[j]=pdba_all.atom[pdbchains[si].start+j];
      snew(chains[i].pdba->atomname[j],1);
      *chains[i].pdba->atomname[j] = 
	strdup(*pdba_all.atomname[pdbchains[si].start+j]);
      /* make all chain identifiers equal to that off the chain */
      chains[i].pdba->atom[j].chain=pdbchains[si].chain;
      chains[i].pdba->pdbinfo[j]=pdba_all.pdbinfo[pdbchains[si].start+j];
      copy_rvec(pdbx[pdbchains[si].start+j],chains[i].x[j]);
    }
    /* Renumber the residues assuming that the numbers are continuous */
    k    = chains[i].pdba->atom[0].resnr;
    nres = chains[i].pdba->atom[chains[i].pdba->nr-1].resnr - k + 1;
    chains[i].pdba->nres = nres;
    for(j=0; j < chains[i].pdba->nr; j++)
      chains[i].pdba->atom[j].resnr -= k;
    srenew(chains[i].pdba->resname,nres);
    for(j=0; j<nres; j++) {
      snew(chains[i].pdba->resname[j],1);
      *chains[i].pdba->resname[j] = strdup(*pdba_all.resname[k+j]);
    }
  }

  printf("There are %d chains and %d blocks of water and "
	 "%d residues with %d atoms\n",
	 nchain-nwaterchain,nwaterchain,
	 pdba_all.atom[natom-1].resnr+1,natom);
	  
  printf("\n  %5s  %4s %6s\n","chain","#res","#atoms");
  for (i=0; (i<nchain); i++)
    printf("  %d '%c'%s  %4d %6d  %s\n",
	   i+1, chains[i].chain, chains[i].chain?"":" ",
	   chains[i].pdba->nres, chains[i].pdba->nr,
	   chains[i].bAllWat ? "(only water)":"");
  printf("\n");
  
  check_occupancy(&pdba_all,opt2fn("-f",NFILE,fnm));
  
  ff=choose_ff();
  printf("Using %s force field\n",ff);
  
  /* Read atomtypes... */
  open_symtab(&tab);
  atype=read_atype(ff,&tab);
    
  /* read residue database */
  printf("Reading residue database... (%s)\n",ff);
  nrtp=read_resall(ff,bts,&restp,atype,&tab);
  if (bNewRTP) {
    fp=ffopen("new.rtp","w");
    print_resall(fp,bts,nrtp,restp,atype);
    fclose(fp);
  }
    
  /* read hydrogen database */
  nah=read_h_db(ff,&ah);
  
  /* Read Termini database... */
  sprintf(fn,"%s-n.tdb",ff);
  nNtdb=read_ter_db(fn,&ntdb,atype);
  sprintf(fn,"%s-c.tdb",ff);
  nCtdb=read_ter_db(fn,&ctdb,atype);
  
  top_fn=ftp2fn(efTOP,NFILE,fnm);
  top_file=ffopen(top_fn,"w");
  print_top_header(top_file,title,FALSE,ff,mHmult);

  nincl=0;
  nmol=0;
  incls=NULL;
  mols=NULL;
  nres=0;
  for(chain=0; (chain<nchain); chain++) {
    /* set pdba, natom and nres to the current chain */
    pdba =chains[chain].pdba;
    x    =chains[chain].x;
    natom=chains[chain].pdba->nr;
    nres =chains[chain].pdba->nres;
    
    if (chains[chain].chain && ( chains[chain].chain != ' ' ) )
      printf("Processing chain %d '%c' (%d atoms, %d residues)\n",
	      chain+1,chains[chain].chain,natom,nres);
    else
      printf("Processing chain %d (%d atoms, %d residues)\n",
	      chain+1,natom,nres);

    process_chain(pdba,x,bUnA,bUnA,bUnA,bLysMan,bAspMan,bGluMan,
		  bHisMan,bCysMan,&nssbonds,&ssbonds,angle,distance);
		  
    if (bSort) {
      block = new_block();
      snew(gnames,1);
      sort_pdbatoms(nrtp,restp,natom,&pdba,&x,block,&gnames);
      natom = remove_duplicate_atoms(pdba,x);
      if (ftp2bSet(efNDX,NFILE,fnm)) {
	if (!bRetainH)
	  fprintf(stderr,"WARNING: without the -reth option the generated "
		  "index file (%s) might be useless\n"
		  "(the index file is generated before hydrogens are added)",
		  ftp2fn(efNDX,NFILE,fnm));
	write_index(ftp2fn(efNDX,NFILE,fnm),block,gnames);
      }
      for(i=0; i < block->nr; i++)
	sfree(gnames[i]);
      sfree(gnames);
      done_block(block);
    } else 
      fprintf(stderr,"WARNING: "
	      "without sorting no check for duplicate atoms can be done\n");
    
    if (debug) {
      if ( chains[chain].chain == '\0' || chains[chain].chain == ' ')
	sprintf(fn,"chain.pdb");
      else
	sprintf(fn,"chain_%c.pdb",chains[chain].chain);
      write_sto_conf(fn,title,pdba,x,NULL,box);
    }
    
    find_nc_ter(pdba,&rN,&rC);
    
    if ( (rN<0) || (rC<0) ) {
      printf("No N- or C-terminus found: "
	     "this chain appears to contain no protein\n");
    } else {
      /* set termini */
      if ( (rN>=0) && (bTerMan || (nNtdb<4)) )
	sel_ntdb=choose_ter(nNtdb,ntdb,"Select N-terminus type (start)");
      else
	if (strncmp(*pdba->resname[pdba->atom[rN].resnr],"PRO",3))
	  sel_ntdb=&(ntdb[1]);
	else
	  sel_ntdb=&(ntdb[3]);
      printf("N-terminus: %s\n",sel_ntdb->name);
      
      if ( (rC>=0) && (bTerMan || (nCtdb<2)) )
	sel_ctdb=choose_ter(nCtdb,ctdb,"Select C-terminus type (end)");
      else
	sel_ctdb=&(ctdb[1]);
      printf("C-terminus: %s\n",sel_ctdb->name);
    }
    
    /* Generate Hydrogen atoms (and termini) in the sequence */
    natom=add_h(&pdba,&x,nah,ah,sel_ntdb,sel_ctdb,rN,rC);
    printf("Now there are %d residues with %d atoms\n",
	   pdba->nres,pdba->nr);
    if (debug) write_pdbfile(debug,title,pdba,x,box,0,TRUE);

    if (debug)
      for(i=0; (i<natom); i++)
	fprintf(debug,"Res %s%d atom %d %s\n",
		*(pdba->resname[pdba->atom[i].resnr]),
		pdba->atom[i].resnr+1,i+1,*pdba->atomname[i]);
    
    strcpy(posre_fn,ftp2fn(efITP,NFILE,fnm));
    
    /* make up molecule name(s) */
    if (chains[chain].bAllWat) 
      sprintf(molname,"Water");
    else if ( chains[chain].chain == '\0' || chains[chain].chain == ' ' )
      sprintf(molname,"Protein");
    else
      sprintf(molname,"Protein_%c",chains[chain].chain);
    
    if ((nchain-nwaterchain>1) && !chains[chain].bAllWat) {
      bITP=TRUE;
      strcpy(itp_fn,top_fn);
      printf("Chain time...\n");
      c=strrchr(itp_fn,'.');
      if ( chains[chain].chain == '\0' || chains[chain].chain == ' ' )
	sprintf(c,".itp");
      else
	sprintf(c,"_%c.itp",chains[chain].chain);
      c=strrchr(posre_fn,'.');
      if ( chains[chain].chain == '\0' || chains[chain].chain == ' ' )
	sprintf(c,".itp");
      else
	sprintf(c,"_%c.itp",chains[chain].chain);
      
      nincl++;
      srenew(incls,nincl);
      incls[nincl-1]=strdup(itp_fn);
      itp_file=ffopen(itp_fn,"w");
    } else
      bITP=FALSE;

    srenew(mols,nmol+1);
    if (chains[chain].bAllWat) {
      mols[nmol].name = strdup("SOL");
      mols[nmol].nr   = pdba->nres;
    } else {
      mols[nmol].name = strdup(molname);
      mols[nmol].nr   = 1;
    }
    nmol++;

    if (bITP)
      print_top_comment(itp_file,title,TRUE);

    if (chains[chain].bAllWat)
      top_file2=NULL;
    else
      if (bITP)
	top_file2=itp_file;
      else
	top_file2=top_file;
    
    pdb2top(top_file2,posre_fn,molname,pdba,&x,atype,&tab,bts,nrtp,restp,
	    sel_ntdb,sel_ctdb,bH14,rN,rC,bAlldih,
	    bDummies,bDummyAromatics,mHmult,nssbonds,ssbonds,NREXCL, 
	    long_bond_dist, short_bond_dist);
    
    if (!chains[chain].bAllWat)
      write_posres(posre_fn,pdba);

    if (bITP)
      fclose(itp_file);

    /* pdba and natom have been reassigned somewhere so: */
    chains[chain].pdba = pdba;
    chains[chain].x = x;
    
    if (debug) {
      if ( chains[chain].chain == '\0' || chains[chain].chain == ' ' )
	sprintf(fn,"chain.pdb");
      else
	sprintf(fn,"chain_%c.pdb",chains[chain].chain);
      write_sto_conf(fn,cool_quote(),pdba,x,NULL,box);
    }
  }
  
  print_top_mols(top_file,title,nincl,incls,nmol,mols);
  fclose(top_file);
  
  /* now merge all chains back together */
  natom=0;
  nres=0;
  for (i=0; (i<nchain); i++) {
    natom+=chains[i].pdba->nr;
    nres+=chains[i].pdba->nres;
  }
  snew(atoms,1);
  init_t_atoms(atoms,natom,FALSE);
  for(i=0; i < atoms->nres; i++)
    sfree(atoms->resname[i]);
  sfree(atoms->resname);
  atoms->nres=nres;
  snew(atoms->resname,nres);
  snew(x,natom);
  k=0;
  l=0;
  for (i=0; (i<nchain); i++) {
    if (nchain>1)
      printf("Including chain %d in system: %d atoms %d residues\n",
	     i+1,chains[i].pdba->nr,chains[i].pdba->nres);
    for (j=0; (j<chains[i].pdba->nr); j++) {
      atoms->atom[k]=chains[i].pdba->atom[j];
      atoms->atom[k].resnr+=l; /* l is processed nr of residues */
      atoms->atomname[k]=chains[i].pdba->atomname[j];
      atoms->atom[k].chain=chains[i].chain;
      copy_rvec(chains[i].x[j],x[k]);
      k++;
    }
    for (j=0; (j<chains[i].pdba->nres); j++) {
      atoms->resname[l]=chains[i].pdba->resname[j];
      l++;
    }
  }
  
  if (nchain>1) {
    fprintf(stderr,"Now there are %d atoms and %d residues\n",k,l);
    print_sums(atoms, TRUE);
  }
  
  fprintf(stderr,"\nWriting coordinate file...\n");
  clear_rvec(box_space);
  if (box[0][0] == 0) 
    gen_box(0,atoms->nr,x,box,box_space,FALSE);
  write_sto_conf(ftp2fn(efSTO,NFILE,fnm),title,atoms,x,NULL,box);

  thanx(stdout);
  
  return 0;
}
