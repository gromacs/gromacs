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

static void rename_pdbres(int natom,t_pdbatom pdba[],char *oldnm,char *newnm,
			  bool bFullCompare)
{
  int i;
  
  for(i=0; (i<natom); i++) {
    if ((bFullCompare && (strcasecmp(pdba[i].resnm,oldnm) == 0)) ||
	(!bFullCompare && strstr(pdba[i].resnm,oldnm) != NULL)) 
      strcpy(pdba[i].resnm,newnm);
  }
}

static void rename_pdbresint(int natom,t_pdbatom pdba[],char *oldnm,
			     char *gettp(int),bool bFullCompare)
{
  int  i,rnr;
  char *ptr;
  
  for(i=0; (i<natom); ) {
    if ((bFullCompare && (strcmp(pdba[i].resnm,oldnm) == 0)) ||
	(!bFullCompare && strstr(pdba[i].resnm,oldnm) != NULL)) {
      rnr=pdba[i].resnr;
      ptr=gettp(rnr);
      while((pdba[i].resnr == rnr) && (i < natom)) {
	strcpy(pdba[i].resnm,ptr);
	i++;
      }
    }
    else
      i++;
  }
}

void write_posres(char *fn,int natom,t_pdbatom pdba[])
{
  FILE *fp;
  int  i;
  
  fp=ffopen(fn,"w");
  fprintf(fp,
	  "; In this topology include file, you will find position\n"
	  "; restraint entries for all the heavy atoms in your original pdb file.\n"
	  "; This means that all the protons which were added by pdb2gmx\n"
	  "; are not restraint. This is especially useful for crystal waters.\n\n"
	  "[ position_restraints ]\n"
	  "; %6s%6s%8s%8s%8s\n","atom","type","fx","fy","fz"
	  );
  for(i=0; (i<natom); i++) {
    if (pdba[i].atomnm[0] != 'H') 
      fprintf(fp,"%6d%6d%8.1f%8.1f%8.1f\n",i+1,1,1000.0,1000.0,1000.0);
  }
  ffclose(fp);
}

void pr_seqres(FILE *fp,int natom,t_pdbatom pdba[])
{
  int i,resnr,nres;
  int line;
  
  resnr = -1;
  nres  = pdba[natom-1].resnr+1;
  line  = 0;
  for(i=0; (i<natom); i++) {
    if (resnr != pdba[i].resnr) {
      if ((line % 13) == 0) 
	fprintf(fp,"\nSEQRES%4d%7d ",1+line/13,nres);
      resnr=pdba[i].resnr;
      fprintf(fp,"%4s",pdba[i].resnm);
      line++;
    }
  }
  fprintf(fp,"\n");
}

int read_pdball(char *inf, char *outf,char *title,
		t_pdbatom **pdbaptr, matrix box, bool bRetainH)
/* Read a pdb file. (containing proteins) */
{
  FILE      *in,*out;
  t_pdbatom *pdba;
  int       natom;
  
  /* READ IT */
  printf("Reading pdb file...\n");
  in=ffopen(inf,"r");
  natom=read_pdbatoms(in,title,&pdba,box,!bRetainH);
  ffclose(in);
  *pdbaptr=pdba;
  printf("'%s': %d atoms\n",title,natom);
  
  /* renumber and sort atoms */
  renumber_pdb(natom,pdba);
  pdba_trimnames(natom,pdba);

  /* Rename residues */
  rename_pdbres(natom,pdba,"SOL","HOH",FALSE);
  rename_pdbres(natom,pdba,"WAT","HOH",FALSE);
  rename_pdbres(natom,pdba,"HEM","HEME",FALSE);

  rename_atoms(natom,pdba);
  
  if (natom == 0)
    return 0;

  /* Dump a 'clean' PDB file */
  out=ffopen(outf,"w");
  /* pr_seqres(out,natom,pdba); */
  print_pdbatoms(out,title,natom,pdba,box);
  ffclose(out);
  
  return natom;
}

void process_chain(int natom, t_pdbatom *pdba, 
		   bool bTrpU,bool bPheU,bool bTyrU,
		   bool bLysH,bool bHisMan,bool bCysMan,
		   int *nssbonds,t_ssbond **ssbonds,
		   real angle,real distance)
{
  /* did this for the whole pdb file already, repeat for every chain */
  renumber_pdb(natom,pdba);
  
  /* Rename aromatics, lys and histidine */
  if (bTyrU) rename_pdbres(natom,pdba,"TYR","TYRU",FALSE);
  if (bTrpU) rename_pdbres(natom,pdba,"TRP","TRPU",FALSE);
  if (bPheU) rename_pdbres(natom,pdba,"PHE","PHEU",FALSE);
  if (bLysH) 
    rename_pdbres(natom,pdba,"LYS","LYSH",FALSE);
  else
    rename_pdbresint(natom,pdba,"LYS",get_lystp,FALSE);
  
  *nssbonds=mk_specbonds(natom,pdba,bCysMan,ssbonds);
  
  if (!bHisMan)
    set_histp(natom,pdba,angle,distance);
  else
    rename_pdbresint(natom,pdba,"HIS",get_histp,TRUE);
}

int pdbcompare(const void *a,const void *b)
{
  t_pdbatom *pa,*pb;
  
  pa=(t_pdbatom *)a;
  pb=(t_pdbatom *)b;
  
  if (pa->resnr == pb->resnr)
    if (pa->atomnr == pb->atomnr)
      return (pa->atomnm[1] - pb->atomnm[1]);
    else
      return (pa->atomnr - pb->atomnr);
  else
    return (pa->resnr - pb->resnr);
}

typedef struct {
  int resnr;
  int j;
  int index;
  char NtH;
} t_pdbindex;
  
int pdbicomp(const void *a,const void *b)
{
  t_pdbindex *pa,*pb;
  
  pa=(t_pdbindex *)a;
  pb=(t_pdbindex *)b;
  
  if (pa->resnr == pb->resnr)
    if (pa->j == pb->j)
      return (pa->NtH - pb->NtH);
    else
      return (pa->j - pb->j);
  else
    return (pa->resnr - pb->resnr);
}

static void sort_pdbatoms(int nrtp,t_restp restp[],
			  int natoms,t_pdbatom *pdba,
			  t_block *block,char ***gnames)
{
  int     i,j;
  t_restp *rptr;
  t_pdbindex *pdbi;
  atom_id *a;
  
  snew(pdbi, natoms);
  
  for(i=0; (i<natoms); i++) {
    if ((rptr=search_rtp(pdba[i].resnm,nrtp,restp)) == NULL)
      fatal_error(0,"Residue type %s not found",pdba[i].resnm);
    for(j=0; (j<rptr->natom); j++)
      if (strcasecmp(pdba[i].atomnm,*(rptr->atomname[j])) == 0)
	break;
    if (j==rptr->natom)
      if ( ( ( pdba[i].resnr == 0) && (pdba[i].atomnm[0] == 'H') &&
	     ( (pdba[i].atomnm[1] == '1') || (pdba[i].atomnm[1] == '2') || 
	       (pdba[i].atomnm[1] == '3') ) ) ||
	   (strcasecmp(pdba[i].atomnm,"OXT")==0) )
	j=1;
      else 
	fatal_error(0,"Atom %s not found in residue %s (looking for '%s %d')"
		    " while sorting atoms",pdba[i].atomnm,
		    rptr->resname,pdba[i].resnm,pdba[i].resnr+1);
    /* make shadow array to be sorted into indexgroup */
    pdbi[i].resnr=pdba[i].resnr;
    pdbi[i].j    =j;
    pdbi[i].index=i;
    pdbi[i].NtH  =pdba[i].atomnm[1];
    
    pdba[i].atomnr=j;
  }
  qsort(pdba,natoms,(size_t)sizeof(pdba[0]),pdbcompare);
  qsort(pdbi,natoms,(size_t)sizeof(pdbi[0]),pdbicomp);
  
  /* make indexgroup in block */
  snew(a,natoms);
  for (i=0; (i<natoms); i++)
    a[i]=pdbi[i].index;
  add_grp(block, gnames, natoms, a, "prot_sort");
  sfree(a);
  sfree(pdbi);
}

static int remove_double_atoms(int natoms,t_pdbatom **pdbaptr)
{
  int     i,j,nres,oldnatoms;
  t_pdbatom *pdba;
  
  pdba=*pdbaptr;
  
  oldnatoms=natoms;
  nres=pdba[natoms-1].resnr;
  printf("Checking for double atoms....\n");
  /* NOTE: natoms is modified inside the loop */
  for(i=1; (i<natoms); i++) {
    if ( (pdba[i-1].resnr     == pdba[i].resnr) &&
	 (pdba[i-1].atomnr    == pdba[i].atomnr) &&
	 (strcmp(pdba[i-1].atomnm,pdba[i].atomnm)==0) && 
	 !( (pdba[i].resnr==nres) && (strcasecmp(pdba[i-1].atomnm,"O")==0) ) ) {
      printf("deleting double atom (%d %s %s %c %d)\n",
	     i+1, pdba[i].atomnm, pdba[i].resnm, 
	     pdba[i].chain, pdba[i].resnr+1);
      natoms--;
      for (j=i; (j<natoms); j++)
	pdba[j]=pdba[j+1];
    }
  }
  if (natoms!=oldnatoms)
    printf("Now there are %d atoms\n",natoms);
  *pdbaptr=pdba;
  
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

void find_nc_ter(int natom,t_pdbatom pdba[],int *rn,int *rc/*,int ter_type[]*/)
{
  int i,rnr;
  
  *rn=-1;
  *rc=-1;
  for(i=0; (i<natom); i++) {
    rnr=pdba[i].resnr;
    /*ter_type[i]=eterNormal;*/
    if (*rn == -1) {
      if (is_protein(pdba[i].resnm)) {
	*rn=rnr;
	/*ter_type[i]=eterNterm;*/
      }
    }
    if (*rc != rnr) {
      if (is_protein(pdba[i].resnm))
	*rc=rnr;
    }
  }
  if (debug)
    fprintf(debug,"nres: %d, rN: %d, rC: %d\n",pdba[natom-1].resnr,*rn,*rc);
}

int main(int argc, char *argv[])
{
  static char *desc[] = {
    "This program reads a pdb file, reads up some database files,",
    "adds hydrogens to the molecules if requested",
    "and generates coordinates in Gromacs (Gromos) format and a topology",
    "in Gromacs format. These files can subsequently be processed to generate",
    "a status file.[PAR]",
    "Note that a pdb file is nothing more than a file format, and it need",
    "not necessarily contain a protein structure. Every kind of molecule",
    "for which there is support in the database can be converted. If there",
    "is no support in the database, you can add it yourself.[PAR]",
    "The program has limited intelligence, it reads a number of",
    "database files, that allow it to make special bonds",
    "(Cys-Cys, Heme-His, etc.), if necessary this can be done manually.",
    "The program can prompt the user to select which kind of LYS, CYS or HIS", 
    "residue she wants. For LYS the choice is between LYS (two protons on",
    "NZ) or LYSH (three protons), for HIS the proton can be either on ND1",
    "(HISA), on NE2 (HISB) or on both (HISH). By default these selections",
    "are done automatically.[PAR]",
    "During processing the atoms will be reorderd according to Gromacs",
    "conventions. With -n an index file can be generated that contains",
    "one group reordered in the same way. This allows you to convert a",
    "GROMOS trajectory and coordinate file to GROMACS. There is one",
    "limitation: reordering is done after the hydrogens are stripped",
    "from the input and before new hydrogens are added. This means that",
    "if you have hydrogens in your input file, you [BB]must[bb] select",
    "the -reth option to obtain a usefull index file."
  };
  static char *bugs[] = {
    "Generation of N-terminal hydrogen atoms on OPLS files does not work",
    "Deuterium (D) is not recognized as a hydrogen and will crash the program.",
    "It is assumed that atomic coordinates pdb files are allways in "
    "Angstrom.",
    "The program should be able to select the protonation on pH, and also "
    "should allow for selection of ASPH instead of ASP.",
  };

  typedef struct {
    char chain;
    int  start;
    int  natom;
    bool bAllWat;
    t_pdbatom *pdba;
    rvec *x;
    t_atoms *atoms;
  } t_chain;
  
  FILE       *fp;
  int        natom,nres;
  t_pdbatom  *pdba_all,*pdba;
  t_atoms    *atoms;
  t_block    *block;
  int        chain,nchain;
  t_chain    *chains;
  char       pchain;
  int        nincl,nmol;
  char       **incls,**mols;
  char       **gnames;
  matrix     box;
  rvec       box_space;
  char       *ff;
  int        i,j,k,l,nrtp,rN,rC;
  t_restp    *restp;
  t_resbond  *rb;
  t_resang   *ra;
  t_idihres  *idih;
  t_addh     *ah;
  t_symtab   tab;
  t_atomtype *atype;
  char       fn[256],top_fn[STRLEN],itp_fn[STRLEN];
  char       molname[STRLEN],title[STRLEN];
  char       *c;
  int        nah,nNtdb,nCtdb;
  t_terblock *ntdb,*ctdb,*sel_ntdb,*sel_ctdb;
  int        nssbonds;
  t_ssbond   *ssbonds;
  rvec       *x,*dummy;
  bool       bTopWritten;
  t_filenm   fnm[] = { 
    { efPDB, "-f", NULL,    ffREAD  },
    { efPDB, "-q", "clean", ffWRITE },
    { efPDB, "-op", "pdbout",ffWRITE },
    { efGRO, "-o", NULL,    ffWRITE },
    { efNDX, "-n", "clean", ffOPTWR },
    { efTOP, NULL, NULL,    ffWRITE },
    { efITP, "-i", "posre", ffWRITE }
  };
#define NFILE asize(fnm)

  /* Command line arguments msut be static */
  static bool bNewRTP=FALSE;
  static bool bInter=FALSE, bLysH=TRUE, bFFMan=FALSE, bCysMan=FALSE; 
  static bool bTerMan=FALSE, bUnA=FALSE;
  static bool bH14= FALSE,bSort=TRUE, bRetainH=FALSE;
  static bool bAlldih=FALSE,bHisMan = FALSE; 
  static real angle=135.0,distance=0.3;
  t_pargs pa[] = {
    { "-newrtp", FALSE,   etBOOL, &bNewRTP,
        "HIDDENWrite the residue database in new format to 'new.rtp'"},
    { "-inter", FALSE,    etBOOL, &bInter,
        "Overrides the next 5 options and makes their selections interactive"},
    { "-ff", FALSE,    etBOOL, &bFFMan, 
	"Interactive Force Field selection, instead of the first one" },
    { "-ss", FALSE,    etBOOL, &bCysMan, 
	"Interactive SS bridge selection" },
    { "-ter", FALSE,    etBOOL, &bTerMan, 
	"Interactive termini selection, instead of charged" },
    { "-lysh", FALSE,  etBOOL, &bLysH,
        "Selects the LysH (charge +1) residue type, instead of interactive "
        "selection" },
    { "-his", FALSE, etBOOL, &bHisMan,
	"Interactive Histidine selection, instead of checking H-bonds" },
    { "-angle", FALSE, etREAL, &angle,
	"Minimum angle for a hydrogen bond (180 = ideal)" },
    { "-dist", FALSE, etREAL,  &distance,
	"Maximum distance for a hydrogen bond  (in nm)" },
    { "-una", FALSE,  etBOOL, &bUnA, 
	"Selects aromatic rings with united CH atoms on Phenylalanine, "
	"Tryptophane and Tyrosine. " },
    { "-sort", FALSE,  etBOOL, &bSort,  
	"Sort the residues according to database, sometimes this is necessary "
	"to get charge groups together" },
    { "-H14",  FALSE,  etBOOL, &bH14, 
      "Use 3rd neighbour interactions for hydrogen atoms" },
    { "-reth", FALSE,  etBOOL, &bRetainH, 
        "Retain hydrogen atoms that are in the pdb file. Their names *must* "
        "match names in the database files used by pdb2gmx. Except for "
        "residues Tyr, Trp, Phe, Lys and His, no additional "
        "hydrogen atoms will be added." },
    { "-alldih",FALSE, etBOOL, &bAlldih, 
        "Generate all proper dihedrals instead of only those with as few "
	"hydrogens as possible (useful for use with Charmm)" }
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
  
  clear_mat(box);
  natom=read_pdball(opt2fn("-f",NFILE,fnm),opt2fn("-q",NFILE,fnm),title,
		    &pdba_all,box,bRetainH);
  
  if (natom==0)
    fatal_error(0,"No atoms found in pdb file %s\n",opt2fn("-f",NFILE,fnm));

  printf("Analyzing pdb file\n");
  pchain='\0';
  nchain=0;
  chains=NULL;
  for (i=0; (i<natom); i++)
    if (pdba_all[i].chain!=pchain) {
      pchain=pdba_all[i].chain;
      /* set natom for previous chain */
      if (nchain > 0)
	chains[nchain-1].natom=i-chains[nchain-1].start;
      /* check if chain identifier was used before */
      for (j=0; (j<nchain); j++)
	if (chains[j].chain == pdba_all[i].chain)
	  fatal_error(0,"Chain identifier '%c' was used "
		      "in two non-sequential blocks (residue %d, atom %d)",
		      pdba_all[i].chain,pdba_all[i].resnr+1,i+1);
      nchain++;
      srenew(chains,nchain);
      chains[nchain-1].chain=pdba_all[i].chain;
      chains[nchain-1].start=i;
      chains[nchain-1].bAllWat=TRUE;
    }
  chains[nchain-1].natom=natom-chains[nchain-1].start;
  /* watch out: dirty loop! 
   * (both 'i' and 'nchain' are modified within the loopbody)         */
  for (i=0; (i<nchain); i++) {
    for (j=0; (j<chains[i].natom) && chains[i].bAllWat; j++)
      chains[i].bAllWat = chains[i].bAllWat && 
	(strcasecmp(pdba_all[chains[i].start+j].resnm,"HOH") == 0);
    if (chains[i].bAllWat && (i>0) && chains[i-1].bAllWat) {
      /* this chain and previous chain contain only water: merge them */
      printf("Merging chains %c and %c\n",chains[i-1].chain,chains[i].chain);
      chains[i-1].natom += chains[i].natom;
      for (j=i+1; (j<nchain); j++) {
	chains[j-1].chain = chains[j].chain;
	chains[j-1].start = chains[j].start;
	chains[j-1].natom = chains[j].natom;
      }
      nchain--;
      srenew(chains,nchain);
      /* this is dirty but necessary: */
      i--;
    }
  }
  /* copy pdb data for all chains */
  for (i=0; (i<nchain); i++) {
    snew(chains[i].pdba, chains[i].natom);
    for (j=0; (j<chains[i].natom); j++)
      chains[i].pdba[j]=pdba_all[chains[i].start+j];
  }
  
  if ((nchain==2) && ( (chains[0].chain==' ') || (chains[1].chain==' ') ) ){
    nchain=1;
    srenew(chains,nchain);
    chains[0].chain=' ';
    chains[0].start=0;
    chains[0].natom=natom;
  }
  /* chains[nchain].start will simply contain the total natom */
  srenew(chains,nchain+1);
  chains[nchain].chain='\0';
  chains[nchain].start=natom;
  
  j=nchain;
  for (i=0; (i<nchain); i++)
    if (chains[i].chain==' ')
      j--;
  if (j==0) j=1;
  
  printf("There are %d chains and %d residues with %d atoms\n",
	  j,pdba_all[natom-1].resnr+1,natom);
	  
  printf("%5s %5s %4s %6s\n","chain","start","#res","#atoms");
  for (i=0; (i<nchain); i++)
    printf("%d '%c' %5d %4d %6d %s\n",
	   i+1,chains[i].chain,chains[i].start+1,
	   pdba_all[chains[i+1].start-1].resnr+1 -
	   pdba_all[chains[i  ].start  ].resnr,
	   chains[i+1].start-chains[i].start,
	   chains[i].bAllWat ? "(only water)":"");
  
  ff=choose_ff(bFFMan);
  printf("Using %s force field\n",ff);
  
  /* Read atomtypes... */
  open_symtab(&tab);
  atype=read_atype(ff,&tab);
    
  /* read residue database */
  printf("Reading residue database... (%s)\n",ff);
  nrtp=read_resall(ff,&restp,&rb,&ra,&idih,atype,&tab);
  if (debug) {
    fprintf(debug,"\nResidue database with %d residues:\n\n\n",nrtp);
    print_resall(debug,nrtp,restp,rb,ra,idih,atype);
    fprintf(debug,"\n");
  }
  if (bNewRTP) {
    fp=ffopen("new.rtp","w");
    print_resall(fp,nrtp,restp,rb,ra,idih,atype);
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
  
  bTopWritten=FALSE;
  nincl=0;
  nmol=0;
  incls=NULL;
  mols=NULL;
  nres=0;
  for(chain=0; (chain<nchain); chain++) {
    /* set pdba and natom to the current chain */
    pdba =chains[chain].pdba;
    natom=chains[chain].natom;
    nres =pdba[natom-1].resnr + 1 - nres;
    
    if (chains[chain].chain && ( chains[chain].chain != ' ' ) )
      printf("Processing chain %d '%c' (%d atoms, %d residues)\n",
	      chain+1,chains[chain].chain,natom,nres);
    else
      printf("Processing chain %d (%d atoms, %d residues)\n",
	      chain+1,natom,nres);

    process_chain(natom,pdba,bUnA,bUnA,bUnA,bLysH,bHisMan,bCysMan,
		  &nssbonds,&ssbonds,angle,distance);
		  
    if (bSort) {
      block = new_block();
      snew(gnames,1);
      sort_pdbatoms(nrtp,restp,natom,pdba,block,&gnames);
      natom = remove_double_atoms(natom, &pdba);
      if (ftp2bSet(efNDX,NFILE,fnm)) {
	if (!bRetainH)
	  fprintf(stderr,"WARNING: without the -reth option the generated "
		  "index file (%s) might be useless\n"
		  "(the index file is generated before hydrogens are added)",
		  ftp2fn(efNDX,NFILE,fnm));
	write_index(ftp2fn(efNDX,NFILE,fnm),block,gnames);
      }
      sfree(block[0].index);
      sfree(block);
      sfree(gnames);
    } else 
      fprintf(stderr,"WARNING: "
	      "without sorting no check for double atoms can be done\n");
    
    if (debug) {
      sprintf(fn,"%c.pdb",chains[chain].chain);
      fp=ffopen(fn,"w");
      print_pdbatoms(fp,title,natom,pdba,box);
      fclose(fp);
    }

    find_nc_ter(natom,pdba,&rN,&rC);

    if ( (rN<0) || (rC<0) )
      printf("No N- or C-terminus found: "
	     "assuming this chain contains no protein\n");
    
    /* set termini */
    if ( (rN>=0) && (bTerMan || (nNtdb<4)) )
      sel_ntdb=choose_ter(nNtdb,ntdb,"Select N-terminus type (start)");
    else {
      if (strncmp(pdba[rN].resnm,"PRO",3))
	sel_ntdb=&(ntdb[1]);
      else
	sel_ntdb=&(ntdb[3]);
      printf("N-terminus: %s\n",sel_ntdb->bname);
    }
    
    if ( (rC>=0) && (bTerMan || (nCtdb<2)) )
      sel_ctdb=choose_ter(nCtdb,ctdb,"Select C-terminus type (end)");
    else {
      sel_ctdb=&(ctdb[1]);
      printf("C-terminus: %s\n",ctdb[1].bname);
    }
    
    /* Generate Hydrogen atoms (and termini) in the sequence */
    natom=add_h(natom,&pdba,nah,ah,&x,sel_ntdb,sel_ctdb,rN,rC);
    printf("Now there are %d residues with %d atoms\n",
	    pdba[natom-1].resnr+1,natom);
    if (debug)
      print_pdbatoms(debug,title,natom,pdba,box);
    
    snew(atoms,1);
    pdb2atoms(natom,pdba,atoms,&dummy,&tab);
    sfree(dummy);
    
    if (debug)
      for(i=0; (i<natom); i++)
	fprintf(debug,"Res %d atom %d %s\n",
		atoms->atom[i].resnr,i,*atoms->atomname[i]);
    
    strcpy(top_fn,ftp2fn(efTOP,NFILE,fnm));
    strcpy(itp_fn,ftp2fn(efITP,NFILE,fnm));
    /* make up molecule name(s) */
    if (chains[chain].bAllWat) 
      sprintf(molname,"Water");
    else if (chains[chain].chain==' ')
      sprintf(molname,"Protein");
    else
      sprintf(molname,"Protein_%c",chains[chain].chain);
    /* make filenames for topol.top/.itp and for posre.itp */
    if ( ! ( (nchain==1) || 
	     ( (chain==nchain-1) && chains[chain].bAllWat ) ) ) {
      c=strrchr(top_fn,'.');
      if ( chains[chain].chain != ' ' )
	sprintf(c,"_%c.itp",chains[chain].chain);
      else 
	sprintf(c,".itp");
      c=strrchr(itp_fn,'.');
      if ( chains[chain].chain != ' ' )
	sprintf(c,"_%c.itp",chains[chain].chain);
      else 
	sprintf(c,".itp");
      
      nincl++;
      srenew(incls,nincl);
      incls[nincl-1]=strdup(top_fn);
    } else  
      bTopWritten=TRUE;
      
    nmol++;
    srenew(mols,nmol);
    mols[nmol-1]=strdup(molname);

    write_posres(itp_fn,natom,pdba);
    
    pdb2top(ff,top_fn,itp_fn,title,molname,
	    nincl,incls,nmol,mols,atoms,nah,ah,x,
	    atype,&tab,nrtp,rb,nrtp,restp,nrtp,ra,nrtp,idih,
	    sel_ntdb,sel_ctdb,bH14,rN,rC,bAlldih,nssbonds,ssbonds);
    
    /* pdba and natom have been reassigned somewhere so: */
    chains[chain].pdba = pdba;
    chains[chain].natom=natom;
    
    /* also save x and atoms */
    chains[chain].x = x;
    chains[chain].atoms=atoms;
    
    if (debug) {
      sprintf(fn,"%c.gro",chains[chain].chain);
      write_conf(fn,cool_quote(),atoms,x,NULL,box);
    }
  }
  /* check if .top file was already written */
  if (!bTopWritten)
    write_top(ff,ftp2fn(efTOP,NFILE,fnm),NULL,title,NULL,
	      nincl,incls,nmol,mols,NULL,NULL,0,NULL,NULL,NULL);
	      
  /* now merge all chains back together */
  natom=0;
  nres=0;
  for (i=0; (i<nchain); i++) {
    natom+=chains[i].natom;
    nres+=chains[i].atoms->nres;
  }
  snew(atoms,1);
  atoms->nr=natom;
  snew(atoms->atom,natom);
  snew(atoms->atomname,natom);
  snew(atoms->chain,natom);
  atoms->nres=nres;
  snew(atoms->resname,nres);
  snew(x,natom);
  k=0;
  l=0;
  for (i=0; (i<nchain); i++) {
    printf("Including chain %d in system: %d atoms %d residues\n",
	   i+1,chains[i].natom,chains[i].atoms->nres);
    for (j=0; (j<chains[i].natom); j++) {
      atoms->atom[k]=chains[i].atoms->atom[j];
      atoms->atom[k].resnr+=l; /* l is processed nr of residues */
      atoms->atomname[k]=chains[i].atoms->atomname[j];
      atoms->chain[k]=chains[i].chain;
      copy_rvec(chains[i].x[j],x[k]);
      k++;
    }
    for (j=0; (j<chains[i].atoms->nres); j++) {
      atoms->resname[l]=chains[i].atoms->resname[j];
      l++;
    }
  }
  printf("Now there are %d atoms and %d residues\n",k,l);
  
  printf("Writing coordinate file...\n");
  clear_rvec(box_space);
  if (box[0][0] == 0) 
    gen_box(0,atoms->nr,x,box,box_space,FALSE);
  write_conf(ftp2fn(efGRO,NFILE,fnm),title,atoms,x,NULL,box);
  write_pdb_conf(opt2fn("-op",NFILE,fnm),title,atoms,x,box,FALSE);

  thanx(stdout);
  
  return 0;
}
