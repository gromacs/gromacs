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
static char *SRCID_pdb2gmx_c = "$Id$";

#include <time.h>
#include <ctype.h>
#include "assert.h"
#include "sysstuff.h"
#include "typedefs.h"
#include "smalloc.h"
#include "copyrite.h"
#include "string2.h"
#include "pdb2gmx.h"
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
#include "gstat.h"
#include "genhydro.h"
#include "readinp.h"
#include "xlate.h"
#include "specbond.h"
#include "index.h"

#define NREXCL 3
char *hh[ehisNR]   = { "HISA", "HISB", "HISH", "HIS1" };

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

static char *get_lystp(int resnr)
{
  enum { elys, elysH, elysNR };
  static char *lh[elysNR] = { "LYS", "LYSH" };
  static char *expl[elysNR] = {
    "Not protonated",
    "Protonated"
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

void write_posres(char *fn,t_atoms *pdba)
{
  FILE *fp;
  int  i;
  
  fp=ffopen(fn,"w");
  fprintf(fp,
	  "; In this topology include file, you will find position restraint\n"
	  "; entries for all the heavy atoms in your original pdb file.\n"
	  "; This means that all the protons which were added by pdb2gmx are\n"
	  "; not restraint. This is especially useful for crystal waters.\n"
	  "\n"
	  "[ position_restraints ]\n"
	  "; %6s%6s%8s%8s%8s\n","atom","type","fx","fy","fz"
	  );
  for(i=0; (i<pdba->nr); i++) {
    if (!is_hydrogen(*pdba->atomname[i])) 
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
    
  printf("'%s': %d atoms\n",title,natom);
  
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
		   bool bLysH,bool bHisMan,bool bCysMan,
		   int *nssbonds,t_ssbond **ssbonds,
		   real angle,real distance)
{
  int i;

  /* Rename aromatics, lys and histidine */
  if (bTyrU) rename_pdbres(pdba,"TYR","TYRU",FALSE);
  if (bTrpU) rename_pdbres(pdba,"TRP","TRPU",FALSE);
  if (bPheU) rename_pdbres(pdba,"PHE","PHEU",FALSE);
  if (bLysH) 
    rename_pdbres(pdba,"LYS","LYSH",FALSE);
  else
    rename_pdbresint(pdba,"LYS",get_lystp,FALSE);

  *nssbonds=mk_specbonds(pdba,x,bCysMan,ssbonds);
  rename_pdbres(pdba,"CYS","CYSH",FALSE);
  for(i=0; i<*nssbonds; i++) {
    if (strcmp(*pdba->resname[(*ssbonds)[i].res1],"CYSH")==0) {
      sfree(*pdba->resname[(*ssbonds)[i].res1]);
      *pdba->resname[(*ssbonds)[i].res1]=strdup("CYS");
    }
    if (strcmp(*pdba->resname[(*ssbonds)[i].res2],"CYSH")==0) {
      sfree(*pdba->resname[(*ssbonds)[i].res2]);
      *pdba->resname[(*ssbonds)[i].res2]=strdup("CYS");
    }
  }
  
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
  
  for(i=0; (i<natoms); i++) {
    atomnm=*pdba->atomname[i];
    resnm=*pdba->resname[pdba->atom[i].resnr];
    if ((rptr=search_rtp(resnm,nrtp,restp)) == NULL)
      fatal_error(0,"Residue type %s not found",resnm);
    for(j=0; (j<rptr->natom); j++)
      if (strcasecmp(atomnm,*(rptr->atomname[j])) == 0)
	break;
    if (j==rptr->natom)
      if ( ( ( pdba->atom[i].resnr == 0) && (atomnm[0] == 'H') &&
	     ( (atomnm[1] == '1') || (atomnm[1] == '2') || 
	       (atomnm[1] == '3') ) ) ||
	   (strcasecmp(atomnm,"OXT")==0) )
	j=1;
      else 
	fatal_error(0,"Atom %s in residue %s %d not found in database"
		    " while sorting atoms",atomnm,
		    rptr->resname,pdba->atom[i].resnr+1);
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
  pdbnew->resname=pdba->resname;
  for (i=0; (i<natoms); i++) {
    pdbnew->atom[i]     = pdba->atom[pdbi[i].index];
    pdbnew->atomname[i] = pdba->atomname[pdbi[i].index];
    pdbnew->pdbinfo[i]  = pdba->pdbinfo[pdbi[i].index];
    copy_rvec((*x)[pdbi[i].index],(*xnew)[i]);
     /* make indexgroup in block */
    a[i]=pdbi[i].index;
  }
  sfree(pdba);
  sfree(*x);
  /* copy the sorted pdbnew back to pdba */
  *pdbaptr=pdbnew;
  *x=*xnew;
  add_grp(block, gnames, natoms, a, "prot_sort");
  sfree(a);
  sfree(pdbi);
}

static int remove_double_atoms(t_atoms *pdba,rvec x[])
{
  int     i,j,natoms,oldnatoms;
  
  printf("Checking for double atoms....\n");
  natoms    = pdba->nr;
  oldnatoms = natoms;
  
  /* NOTE: natoms is modified inside the loop */
  for(i=1; (i<natoms); i++) {
    if ( (pdba->atom[i-1].resnr == pdba->atom[i].resnr) &&
	 (strcmp(*pdba->atomname[i-1],*pdba->atomname[i])==0) 
	 ) {
      printf("deleting double atom %4s  %s%4d",
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
      natoms--;
      for (j=i; (j<natoms); j++) {
	pdba->atom[j]     = pdba->atom[j+1];
	pdba->atomname[j] = pdba->atomname[j+1];
	pdba->pdbinfo[j]  = pdba->pdbinfo[j+1];
	copy_rvec(x[j+1],x[j]);
      }
    }
  }
  pdba->nr=natoms;
  if (natoms != oldnatoms)
    printf("Now there are %d atoms\n",natoms);
  
  return natoms;
}

static char *choose_ff(bool bFFMan)
{
  typedef struct { char *desc,*fn; } t_fff;
  static  char *fnsel;
  FILE    *in;
  t_fff   *fff;
  int     i,nff,sel;
  char    *c,buf[256],fn[32];
  
  in=libopen("FF.dat");
  fgets2(buf,255,in);
  sscanf(buf,"%d",&nff);
  snew(fff,nff);
  for(i=0; (i<nff); i++) {
    fgets2(buf,255,in);
    sscanf(buf,"%s",fn);
    fff[i].fn=strdup(fn);
    /* Search for next non-space character, there starts description */
    c=&(buf[strlen(fn)+1]);
    while (isspace(*c)) c++;
    fff[i].desc=strdup(c);
  }
  fclose(in);

  if (bFFMan && (nff > 1)) {
    printf("\nSelect the Force Field:\n");
    for(i=0; (i<nff); i++)
      printf("%2d: %s\n",i,fff[i].desc);
    do {
      scanf("%d",&sel);
    } while ((sel < 0) || (sel >= nff));
  }
  else
    sel=0;

  fnsel=strdup(fff[sel].fn);

  for(i=0; (i<nff); i++) {
    sfree(fff[i].desc);
    sfree(fff[i].fn);
  }
  sfree(fff);
  
  return fnsel;
}

void find_nc_ter(int natom,t_atoms *pdba,int *rn,int *rc/*,int ter_type[]*/)
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
    "This program reads a pdb file, reads up some database files, adds",
    "hydrogens to the molecules if requested and generates coordinates",
    "in Gromacs (Gromos) format and a topology in Gromacs format. These",
    "files can subsequently be processed to generate a status file.[PAR]",
    
    "Note that a pdb file is nothing more than a file format, and it",
    "need not necessarily contain a protein structure. Every kind of",
    "molecule for which there is support in the database can be converted.",
    "If there is no support in the database, you can add it yourself.[PAR]",
    
    "The program has limited intelligence, it reads a number of database",
    "files, that allow it to make special bonds (Cys-Cys, Heme-His, etc.),",
    "if necessary this can be done manually. The program can prompt the",
    "user to select which kind of LYS, CYS or HIS residue she wants. For",
    "LYS the choice is between LYS (two protons on NZ) or LYSH (three",
    "protons), for HIS the proton can be either on ND1 (HISA), on NE2",
    "(HISB) or on both (HISH). By default these selections are done",
    "automatically. For His, this is based on an optimal hydrogen bonding",
    "conformation. Hydrogen bonds are defined based on a simple geometric",
    "criterium, specified by the maximum hydrogen-donor-acceptor angle",
    "and donor-acceptor distance, which are set by [TT]-angle[tt] and",
    "[TT]-dist[tt] respectively.[PAR]",
    
    "During processing the atoms will be reordered according to Gromacs",
    "conventions. With [TT]-n[tt] an index file can be generated that",
    "contains one group reordered in the same way. This allows you to",
    "convert a Gromos trajectory and coordinate file to Gromos. There is",
    "one limitation: reordering is done after the hydrogens are stripped",
    "from the input and before new hydrogens are added. This means that",
    "if you have hydrogens in your input file, you [BB]must[bb] select",
    "the [TT]-reth[tt] option to obtain a useful index file.[PAR]",
    
    "When using [TT]-reth[tt] to keep all hydrogens from the [TT].pdb[tt]",
    "file, the names of the hydrogens in the [TT].pdb[tt] file [IT]must[it]",
    "match names in the database files used by pdb2gmx. Except for residues",
    "Tyr, Trp, Phe, Lys and His, no additional hydrogen atoms will be",
    "added.[PAR]",
    
    "[TT]-sort[tt] will sort all residues according to the order in the",
    "database, sometimes this is necessary to get charge groups",
    "together.[PAR]",
    
    "[TT]-alldih[tt] will generate all proper dihedrals instead of only",
    "those with as few hydrogens as possible, this is useful for use with",
    "the Charmm forcefield.[PAR]",
    
    "The option -convert removes or slows down hydrogen motions. Angular",
    "and out-of-plane motions can be removed by changing hydrogens into",
    "dummy atoms and fixing angles, which fixes their position relative",
    "to neighboring atoms. Slowing down of dihedral motion is done by",
    "increasing the hydrogen-mass by a factor of 4. This is also done",
    "for water hydrogens to slow down the rotational motion of water.",
    "The increase in mass of the hydrogens is subtracted from the bonded",
    "(heavy) atom so that the total mass of the system remains the same.",
    "Three options are available: normal topology; dummy hydrogens and",
    "fixed angles; increase mass of hydrogens (heavy); both dummies and",
    "increased mass."
  };
  static char *bugs[] = {
    "Generation of N-terminal hydrogen atoms on OPLS files does not work.",
    "Deuterium (D) is not recognized as a hydrogen and will crash the "
    "program.",
    "It is assumed that atomic coordinates in pdb files are in Angstrom.",
    "The program should be able to select the protonation on pH, and also "
    "should allow for selection of ASPH instead of ASP.",
  };
  
  typedef struct {
    char chain;
    int  start;
    int  natom;
    int  nres;
    bool bAllWat;
    t_atoms *pdba;
    rvec *x;
  } t_chain;
  
  FILE       *fp,*top_file,*top_file2,*itp_file;
  int        natom,nres;
  t_atoms    pdba_all,*pdba;
  t_atoms    *atoms;
  t_block    *block;
  int        chain,nchain,nwaterchain;
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
  t_restp    *restp;
  t_resbond  *rb;
  t_resang   *ra;
  t_resdih   *rd;
  t_idihres  *idih;
  t_addh     *ah;
  t_symtab   tab;
  t_atomtype *atype;
  char       fn[256],*top_fn,itp_fn[STRLEN],posre_fn[STRLEN];
  char       molname[STRLEN],title[STRLEN];
  char       *c;
  int        nah,nNtdb,nCtdb;
  t_hackblock *ntdb,*ctdb,*sel_ntdb,*sel_ctdb;
  int        nddb;
  t_dumblock *ddb;
  int        nssbonds;
  t_ssbond   *ssbonds;
  rvec       *pdbx,*x;
  bool       bUsed,bDummies,bWat,bPrevWat,bITP;
  real       mHmult;
  
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
  static bool bInter=FALSE, bLysH=TRUE, bFFMan=FALSE, bCysMan=FALSE; 
  static bool bTerMan=FALSE, bUnA=FALSE;
  static bool bH14= FALSE,bSort=TRUE, bRetainH=FALSE;
  static bool bAlldih=FALSE,bHisMan = FALSE;
  static real angle=135.0,distance=0.3;
  static char *dumstr[] = { NULL, "normal", "dummy", "heavy", "both", NULL };
  t_pargs pa[] = {
    { "-newrtp", FALSE, etBOOL, &bNewRTP,
      "HIDDENWrite the residue database in new format to 'new.rtp'"},
    { "-inter",  FALSE, etBOOL, &bInter,
      "Set the next 5 options to interactive"},
    { "-ff",     FALSE, etBOOL, &bFFMan, 
      "Interactive Force Field selection, instead of the first one" },
    { "-ss",     FALSE, etBOOL, &bCysMan, 
      "Interactive SS bridge selection" },
    { "-ter",    FALSE, etBOOL, &bTerMan, 
      "Interactive termini selection, instead of charged" },
    { "-lysh",   FALSE, etBOOL, &bLysH, 
      "Select the LysH (charge +1) residue type, "
      "instead of interactive selection" },
    { "-his",    FALSE, etBOOL, &bHisMan,
      "Interactive Histidine selection, instead of checking H-bonds" },
    { "-angle",  FALSE, etREAL, &angle, 
      "Minimum hydrogen-donor-acceptor angle for a H-bond (degrees)" },
    { "-dist",   FALSE, etREAL, &distance,
      "Maximum donor-acceptor distance for a H-bond (nm)" },
    { "-una",    FALSE, etBOOL, &bUnA, 
      "Select aromatic rings with united CH atoms on Phenylalanine, "
      "Tryptophane and Tyrosine" },
    { "-sort",   FALSE, etBOOL, &bSort, 
      "Sort the residues according to database" },
    { "-H14",    FALSE, etBOOL, &bH14, 
      "Use 3rd neighbour interactions for hydrogen atoms" },
    { "-reth",   FALSE, etBOOL, &bRetainH, 
      "Retain hydrogen atoms that are in the pdb file" },
    { "-alldih", FALSE, etBOOL, &bAlldih, 
      "Generate all proper dihedrals" },
    { "-convert",FALSE, etENUM, &dumstr, 
      "Convert atoms" }
  };
#define NPARGS asize(pa)
  
  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    asize(bugs),bugs);
  if (bInter) {
    /* if anything changes here, also change description of -inter */
    bFFMan=TRUE;
    bCysMan=TRUE;
    bTerMan=TRUE;
    bLysH=FALSE;
    bHisMan=TRUE;
  }
  
  switch(dumstr[0][0]) {
  case 'n': /* normal */
    bDummies=FALSE;
    mHmult=1.0;
    break;
  case 'd': /* dummy */
    bDummies=TRUE;
    mHmult=1.0;
    break;
  case 'h': /* heavy */
    bDummies=FALSE;
    mHmult=4.0;
    break;
  case 'b': /* both */
    bDummies=TRUE;
    mHmult=4.0;
    break;
  default:
    fatal_error(0,"DEATH HORROR in pdb2gmx: dumstr[0]='%s'",dumstr[0]);
  }/* end switch */
  
  clear_mat(box);
  natom=read_pdball(opt2fn("-f",NFILE,fnm),opt2fn_null("-q",NFILE,fnm),title,
		    &pdba_all,&pdbx,box,bRetainH);
  
  if (natom==0)
    fatal_error(0,"No atoms found in pdb file %s.\n",opt2fn("-f",NFILE,fnm));

  printf("Analyzing pdb file\n");
  nchain=0;
  nwaterchain=0;
  /* keep the compiler happy */
  pchain='?';
  chains=NULL;
  for (i=0; (i<natom); i++) {
    bWat = strcasecmp(*pdba_all.resname[pdba_all.atom[i].resnr],"HOH") == 0;
    if ((i==0) || (pdba_all.atom[i].chain!=pchain) || (bWat != bPrevWat)) {
      pchain=pdba_all.atom[i].chain;
      /* set natom for previous chain */
      if (nchain > 0)
	chains[nchain-1].natom=i-chains[nchain-1].start;
      /* check if chain identifier was used before */
      for (j=0; (j<nchain); j++)
	if ((chains[j].chain != '\0') && (chains[j].chain != ' ') &&
	    (chains[j].chain == pdba_all.atom[i].chain))
	  fatal_error(0,"Chain identifier '%c' was used "
		      "in two non-sequential blocks (residue %d, atom %d)",
		      pdba_all.atom[i].chain,pdba_all.atom[i].resnr+1,i+1);
      srenew(chains,nchain+1);
      chains[nchain].chain=pdba_all.atom[i].chain;
      chains[nchain].start=i;
      chains[nchain].bAllWat=bWat;
      if (bWat)
	nwaterchain++;
      nchain++;
    }
    bPrevWat=bWat;
  }
  chains[nchain-1].natom=natom-chains[nchain-1].start;

  /* copy pdb data and x for all chains */
  for (i=0; (i<nchain); i++) {
    /* check for empty chain identifiers */
    if ((nchain-nwaterchain>1) && !chains[i].bAllWat && 
	((chains[i].chain=='\0') || (chains[i].chain==' '))) {
      bUsed=TRUE;
      for(k='A'; (k<='Z') && bUsed; k++) {
	bUsed=FALSE;
	for(j=0; j<nchain; j++)
	  bUsed=bUsed || (!chains[j].bAllWat && (chains[j].chain==k));
	if (!bUsed) {
	  printf("Gave chain %d chain identifier '%c'\n",i+1,k);
	  chains[i].chain=k;
	}
      }
    }
    snew(chains[i].pdba,1);
    init_t_atoms(chains[i].pdba,chains[i].natom,TRUE);
    snew(chains[i].x,chains[i].natom);
    for (j=0; j<chains[i].natom; j++) {
      chains[i].pdba->atom[j]=pdba_all.atom[chains[i].start+j];
      snew(chains[i].pdba->atomname[j],1);
      *chains[i].pdba->atomname[j] = 
	strdup(*pdba_all.atomname[chains[i].start+j]);
      /* make all chain identifiers equal to that off the chain */
      chains[i].pdba->atom[j].chain=chains[i].chain;
      chains[i].pdba->pdbinfo[j]=pdba_all.pdbinfo[chains[i].start+j];
      copy_rvec(pdbx[chains[i].start+j],chains[i].x[j]);
    }
    /* Renumber the residues assuming that the numbers are continuous */
    k=chains[i].pdba->atom[0].resnr;
    nres=chains[i].pdba->atom[chains[i].natom-1].resnr - k + 1;
    chains[i].pdba->nres=nres;
    for(j=0; j<chains[i].natom; j++)
      chains[i].pdba->atom[j].resnr-=k;
    snew(chains[i].pdba->resname,nres);
    for(j=0; j<nres; j++) {
      snew(chains[i].pdba->resname[j],1);
      *chains[i].pdba->resname[j] = strdup(*pdba_all.resname[k+j]);
    }
  }

  /* chains[nchain].start will simply contain the total natom */
  srenew(chains,nchain+1);
  chains[nchain].chain='\0';
  chains[nchain].start=natom;
  
  for (i=0; (i<nchain); i++)
    chains[i].nres = (pdba_all.atom[chains[i+1].start-1].resnr+1 -
		      pdba_all.atom[chains[i  ].start  ].resnr);
  
  printf("There are %d chains and %d blocks of water and %d residues with %d atoms\n",
	  nchain-nwaterchain,nwaterchain,pdba_all.atom[natom-1].resnr+1,natom);
	  
  printf("\n  %5s  %4s %6s\n","chain","#res","#atoms");
  for (i=0; (i<nchain); i++)
    printf("  %d '%c'  %4d %6d  %s\n",
	   i+1,chains[i].chain,chains[i].nres,
	   chains[i+1].start-chains[i].start,
	   chains[i].bAllWat ? "(only water)":"");
  
  ff=choose_ff(bFFMan);
  printf("\nUsing %s force field\n",ff);
  
  /* Read atomtypes... */
  open_symtab(&tab);
  atype=read_atype(ff,&tab);
    
  /* read residue database */
  printf("Reading residue database... (%s)\n",ff);
  nrtp=read_resall(ff,&restp,&rb,&ra,&rd,&idih,atype,&tab);
  if (debug) {
    fprintf(debug,"\nResidue database with %d residues:\n\n\n",nrtp);
    print_resall(debug,nrtp,restp,rb,ra,rd,idih,atype);
    fprintf(debug,"\n");
  }
  if (bNewRTP) {
    fp=ffopen("new.rtp","w");
    print_resall(fp,nrtp,restp,rb,ra,rd,idih,atype);
    fclose(fp);
  }
    
  /* read hydrogen database */
  if (bRetainH) {
    nah=0;
    ah=NULL;
  } else {
    nah=read_h_db(ff,&ah);
  }
  if (debug) {
    fprintf(debug,"Hydrogen database:\n");
    print_h_db(debug,nah,ah);
  }  
    
  /* Read Termini database... */
  sprintf(fn,"%s-n.tdb",ff);
  nNtdb=read_ter_db(fn,&ntdb,atype);
  sprintf(fn,"%s-c.tdb",ff);
  nCtdb=read_ter_db(fn,&ctdb,atype);
  
  /* Read dummies database */
  nddb=0;
  ddb=NULL;
  if (bDummies)
    nddb=read_dum_db(ff,&ddb);
  if (debug) print_dum_db(stderr,nddb,ddb);

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
    natom=chains[chain].natom;
    nres =chains[chain].nres;
    
    if (chains[chain].chain && ( chains[chain].chain != ' ' ) )
      printf("Processing chain %d '%c' (%d atoms, %d residues)\n",
	      chain+1,chains[chain].chain,natom,nres);
    else
      printf("Processing chain %d (%d atoms, %d residues)\n",
	      chain+1,natom,nres);

    process_chain(pdba,x,bUnA,bUnA,bUnA,bLysH,bHisMan,bCysMan,
		  &nssbonds,&ssbonds,angle,distance);
		  
    if (bSort) {
      block = new_block();
      snew(gnames,1);
      sort_pdbatoms(nrtp,restp,natom,&pdba,&x,block,&gnames);
      natom = remove_double_atoms(pdba,x);
      if (ftp2bSet(efNDX,NFILE,fnm)) {
	if (!bRetainH)
	  fprintf(stderr,"WARNING: without the -reth option the generated "
		  "index file (%s) might be useless\n"
		  "(the index file is generated before hydrogens are added)",
		  ftp2fn(efNDX,NFILE,fnm));
	write_index(ftp2fn(efNDX,NFILE,fnm),block,gnames);
      }
      done_block(block);
      sfree(gnames);
    } else 
      fprintf(stderr,"WARNING: "
	      "without sorting no check for double atoms can be done\n");
    
    if (debug) {
      if ( chains[chain].chain == '\0' || chains[chain].chain == ' ')
	sprintf(fn,"chain.pdb");
      else
	sprintf(fn,"chain_%c.pdb",chains[chain].chain);
      write_sto_conf(fn,title,pdba,x,NULL,box);
    }
    
    find_nc_ter(natom,pdba,&rN,&rC);
    
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
      printf("N-terminus: %s\n",sel_ntdb->bname);
      
      if ( (rC>=0) && (bTerMan || (nCtdb<2)) )
	sel_ctdb=choose_ter(nCtdb,ctdb,"Select C-terminus type (end)");
      else
	sel_ctdb=&(ctdb[1]);
      printf("C-terminus: %s\n",sel_ctdb->bname);
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
    
    /* make filenames for topol.top/.itp and for posre.itp */
    /*
    if ( ! ( (nchain==1) || 
	     ( (chain==nchain-1) && chains[chain].bAllWat ) ) ) {
	     */
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

    write_posres(posre_fn,pdba);

    if (bITP)
      print_top_comment(itp_file,title,TRUE);

    if (chains[chain].bAllWat)
      top_file2=NULL;
    else
      if (bITP)
	top_file2=itp_file;
      else
	top_file2=top_file;

    pdb2top(ff,top_file2,posre_fn,molname,nincl,incls,nmol,mols,pdba,nah,ah,
	    &x,atype,&tab,nrtp,rb,nrtp,restp,nrtp,ra,nrtp,rd,nrtp,idih,
	    sel_ntdb,sel_ctdb,bH14,rN,rC,bAlldih,
	    nddb,ddb,bDummies,mHmult,nssbonds,ssbonds,NREXCL);
    
    if (bITP)
      fclose(itp_file);

    /* pdba and natom have been reassigned somewhere so: */
    chains[chain].pdba = pdba;
    chains[chain].natom= pdba->nr;
    chains[chain].x = x;
    
    if (debug) {
      if ( chains[chain].chain == '\0' || chains[chain].chain == ' ' )
	sprintf(fn,"chain.pdb");
      else
	sprintf(fn,"chain_%c.pdb",chains[chain].chain);
      write_sto_conf(fn,cool_quote(),pdba,x,NULL,box);
    }
  }
  /* check if .top file was already written */
  /*
  if (!bTopWritten)
    write_top(ff,ftp2fn(efTOP,NFILE,fnm),NULL,title,NULL,
	      nincl,incls,nmol,mols,NULL,NULL,NULL,NULL,NULL,NREXCL,mHmult);
	      */
  print_top_mols(top_file,nincl,incls,nmol,mols);
  fclose(top_file);
  
  /* now merge all chains back together */
  natom=0;
  nres=0;
  for (i=0; (i<nchain); i++) {
    natom+=chains[i].natom;
    nres+=chains[i].pdba->nres;
  }
  snew(atoms,1);
  init_t_atoms(atoms,natom,FALSE);
  atoms->nres=nres;
  snew(atoms->resname,nres);
  snew(x,natom);
  k=0;
  l=0;
  for (i=0; (i<nchain); i++) {
    if (nchain>1)
      printf("Including chain %d in system: %d atoms %d residues\n",
	     i+1,chains[i].natom,chains[i].pdba->nres);
    for (j=0; (j<chains[i].natom); j++) {
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
