/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
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
 * Great Red Oystrich Makes All Chemists Sane
 */
#include <pwd.h>
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

  printf("Which %s type do you want for residue %d\n",title,resnr);
  for(sel=0; (sel < nr); sel++)
    printf("%d. %s (%s)\n",sel,expl[sel],name[sel]);
  printf("\nType a number:"); fflush(stdout);

  if (scanf("%d",&sel) != 1)
    fatal_error(0,"Answer me for res %s %d!",title,resnr);
  
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
  static char *message =
    "; In this topology include file, you will find position\n"
    "; restraint entries for all the heavy atoms in your original pdb file.\n"
    "; This means that all the protons which where added by pdb2gmx\n"
    "; are not restraint. This is especially useful for crystal waters.\n\n"
    "[ position_restraints ]\n"
    "; atom  type      fx      fy      fz\n";
  FILE *fp;
  int  i;
  
  fp=ffopen(fn,"w");
  fprintf(fp,message);
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

int read_pdball(char *inf,char *outf,
		t_pdbatom **pdbaptr, matrix box,
		bool bTrpU,bool bPheU,bool bTyrU,
		bool bLysH,bool bHisMan,
		bool bRetainH,bool bCysMan,
		int *nssbonds,t_ssbond **ssbonds,
		real angle,real distance)
/* Read a pdb file. (containing proteins) */
{
  FILE      *in,*out;
  t_pdbatom *pdba;
  int       natom;
  
  /* READ IT */
  fprintf(stderr,"Reading pdb file...\n");
  in=ffopen(inf,"r");
  natom=read_pdbatoms(in,&pdba,box,!bRetainH);
  ffclose(in);
  
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
  print_pdbatoms(out,natom,pdba,box);
  ffclose(out);
  
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

  /* *nssbonds=ss_bonds(natom,pdba,bCysMan,ssbonds); */
  
  *pdbaptr=pdba;
  
  return natom;
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

static void renumber_res(int nrtp,t_restp restp[],
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
		    " while renumbering residues",
		    pdba[i].atomnm,rptr->resname,pdba[i].resnm,pdba[i].resnr);
    /* make shadow array to be sorted into indexgroup */
    pdbi[i].resnr=pdba[i].resnr;
    pdbi[i].j    =j;
    pdbi[i].index=i;
    pdbi[i].NtH  =pdba[i].atomnm[1];
    pdba[i].atomnr=j;
  }
  qsort(pdba,natoms,sizeof(pdba[0]),pdbcompare);
  qsort(pdbi,natoms,sizeof(pdbi[0]),pdbicomp);
  
  /* make indexgroup in block */
  snew(a,natoms);
  for (i=0; (i<natoms); i++)
    a[i]=pdbi[i].index;
  add_grp(block, gnames, natoms, a, "prot_sort");
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
    for(c=&(buf[strlen(fn)+1]); isspace(*c); c++)
      ;
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

bool is_prot(char *key)
{
  static char **prot=NULL;
  static int  nprot;
  int    i;
  
  if (prot == NULL) {
    nprot=get_strings("aminoacids.dat",&prot);
  }
  for(i=0; (i<nprot); i++) 
    if (strcasecmp(prot[i],key) == 0)
      return TRUE;
  return FALSE;
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
      if (is_prot(pdba[i].resnm)) {
	*rn=rnr;
	/*ter_type[i]=eterNterm;*/
      }
    }
    if (*rc != rnr) {
      if (is_prot(pdba[i].resnm))
	*rc=rnr;
    }
  }
  fprintf(stderr,"nres: %d, rN: %d, rC: %d\n",pdba[natom-1].resnr,*rn,*rc);
}

static void analyse_pdba(int natom,t_pdbatom pdba[])
{
  int  i,nres;
  
  if (natom == 0) {
    fprintf(stderr,"No atoms in pdb file");
    exit(1);
  }
  
  nres = pdba[natom-1].resnr+1;
  fprintf(stderr,"There are %d residues with %d atoms\n",nres,natom);
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
    "residue she wants. For LYS the choice is between LYS (two protons on NZ)", 
    "or LYSH (three protons), ",
    "for HIS the proton can be either on ND1 (HISA), on NE2 (HISB)",
    "or on both (HISH)."
  };
  static char *bugs[] = {
    "Generation of N-terminal hydrogen atoms on OPLS files does not work",
    "Deuterium (D) is not recognized as a hydrogen and will  crash the program.",
    "The HETATM format in pdbfiles is not yet supported (this is used for ions, solvent etc.).",
    "It is assumed that atomic coordinates pdb files are allways in Angstrom.",
    "The program should be able to select the protonation on pH, and also should allow for selection of ASPH instead of ASP.",
  };

  FILE       *fp;
  int        natom;
  t_pdbatom  *pdba;
  t_atoms    *atoms;
  t_block    *block;
  char       **gnames;
  matrix     box;
  rvec       box_space;
  char       *ff;
  int        i,nrtp,rN,rC;
  int        *molnr;
  t_restp    *restp;
  t_resbond  *rb;
  t_resang   *ra;
  t_idihres  *idih;
  t_addh     *ah;
  t_symtab   tab;
  t_atomtype *atype;
  char       fn[256];
  int        nah,nNtdb,nCtdb;
  t_terblock *ntdb,*ctdb,*sel_ntdb,*sel_ctdb;
  int        nssbonds;
  t_ssbond   *ssbonds;
  rvec       *x,*v,*dummy;
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
        "Overrides the next 7 options and makes their selections interactive"},
    { "-ff", FALSE,    etBOOL, &bFFMan, 
	"Interactive Force Field selection, instead of the first one" },
    { "-ss", FALSE,    etBOOL, &bCysMan, 
	"Interactive SS bridge selection" },
    { "-ter", FALSE,    etBOOL, &bTerMan, 
	"Interactive termini selection, instead of charged" },
    { "-his", FALSE, etBOOL, &bHisMan,
	"Interactive Histidine selection, instead of checking H-bonds" },
    { "-angle", FALSE, etREAL, &angle,
	"Minimum angle for hydrogen bonds (180 = ideal)" },
    { "-dist", FALSE, etREAL,  &distance,
	"Maximum distance for hydrogen bonds (in nm)" },
    { "-lysh", FALSE,  etBOOL, &bLysH,
        "Selects the LysH (charge +1) residue type, instead of interactive selection" },
    { "-una", FALSE,  etBOOL, &bUnA, 
	"Selects aromatic rings with united CH atoms on Phenylalanine, Tryptophane and Tyrosine. " },
    { "-sort", FALSE,  etBOOL, &bSort,  
	"Sort the residues according to database, sometimes this is necessary to get charge groups together" },
    { "-H14",  FALSE,  etBOOL, &bH14, 
      "Use 3rd neighbour interactions for hydrogen atoms" },
    { "-reth", FALSE,  etBOOL, &bRetainH, 
      "Retain hydrogen atoms that are in the pdb file. Their names *must* "
      "match names in the database files used by pdb2gmx. Except for "
      "residues Tyr, Trp, Phe, Lys and His, no additional "
      "hydrogen atoms will be added." },
    { "-alldih",FALSE, etBOOL, &bAlldih, 
      "Generate all proper dihedrals instead of only those with as few hydrogens as possible"
    "(useful for use with Charmm)" }
  };
#define NPARGS asize(pa)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,asize(desc),desc,
		    asize(bugs),bugs);
  if (bInter) {
    bFFMan=TRUE;
    bCysMan=TRUE;
    bTerMan=TRUE;
    bHisMan=TRUE;
    bLysH=FALSE;
  }

  natom=read_pdball(opt2fn("-f",NFILE,fnm),
		    opt2fn("-q",NFILE,fnm),&pdba,box,
		    bUnA,bUnA,bUnA,bLysH,bHisMan,
		    bRetainH,bCysMan,&nssbonds,&ssbonds,
		    angle,distance);
  analyse_pdba(natom,pdba);
  
  ff=choose_ff(bFFMan);
  fprintf(stderr,"Using %s force field\n",ff);
  
  /* Read atomtypes... */
  open_symtab(&tab);
  atype=read_atype(ff,&tab);
  
  fprintf(stderr,"Reading residue database... (%s)\n",ff);
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

  if (bSort) {
    block = new_block();
    snew(gnames,1);
    renumber_res(nrtp,restp,natom,pdba,block,&gnames);
    if (ftp2bSet(efNDX,NFILE,fnm)) {
      if (!opt2bSet("-reth",NFILE,fnm))
	fprintf(stderr,"WARNING: without the -reth option the generated .ndx file might be useless (the index file is generated before hydrogens are added)");
      write_index(ftp2fn(efNDX,NFILE,fnm),block,gnames);
    }
  }

  find_nc_ter(natom,pdba,&rN,&rC);

  /* Read Termini database... */
  sprintf(fn,"%s-n.tdb",ff);
  nNtdb=read_ter_db(fn,&ntdb,atype);
  if (bTerMan || (nNtdb<4))
    sel_ntdb=choose_ter(nNtdb,ntdb,"Select N-terminus type (start)");
  else {
    if (strncmp(pdba[rN].resnm,"PRO",3))
      sel_ntdb=&(ntdb[1]);
    else
      sel_ntdb=&(ntdb[3]);
    fprintf(stderr,"N-terminus: %s\n",sel_ntdb->bname);
  }

  sprintf(fn,"%s-c.tdb",ff);
  nCtdb=read_ter_db(fn,&ctdb,atype);
  if (bTerMan || (nCtdb<2))
    sel_ctdb=choose_ter(nCtdb,ctdb,"Select C-terminus type (end)");
  else {
    sel_ctdb=&(ctdb[1]);
    fprintf(stderr,"C-terminus: %s\n",ctdb[1].bname);
  }
  
  if (bRetainH) {
    nah=0;
    ah=NULL;
  } else {
    nah=read_h_db(ff,&ah);
  }
  /* Generate Hydrogen atoms in the sequence */
  if (debug) {
    fprintf(debug,"Hydrogen database:\n");
    print_h_db(debug,nah,ah);
  }  
  natom=add_h(natom,&pdba,nah,ah,&x,sel_ntdb,sel_ctdb,rN,rC);
  fprintf(stderr,"Now there are %d residues with %d atoms\n",
	  pdba[natom-1].resnr+1,natom);
  snew(v,natom);
  if (debug)
    print_pdbatoms(debug,natom,pdba,box);

  write_posres(ftp2fn(efITP,NFILE,fnm),natom,pdba);
  
  snew(atoms,1);
  pdb2atoms(natom,pdba,atoms,&dummy,&tab);
  sfree(dummy);
  
  if (debug)
    for(i=0; (i<natom); i++)
      fprintf(debug,"Res %d atom %d %s\n",
	      atoms->atom[i].resnr,i,*atoms->atomname[i]);

  /* Computing molecule numbers */
  snew(molnr,natom);
  for (i=1; (i<natom); i++) {
    if (pdba[i].chain == pdba[i-1].chain) 
      molnr[i]=molnr[i-1];
    else
      molnr[i]=molnr[i-1]+1;
  }
  
  pdb2top(ff,ftp2fn(efTOP,NFILE,fnm),ftp2fn(efITP,NFILE,fnm),
	  atoms,nah,ah,x,
	  atype,&tab,nrtp,rb,nrtp,restp,nrtp,ra,nrtp,idih,sel_ntdb,sel_ctdb,
	  molnr,bH14,rN,rC,bAlldih,nssbonds,ssbonds);

  fprintf(stderr,"Writing coordinate file...\n");
  clear_rvec(box_space);
  if (box[0][0] == 0) 
    gen_box(0,atoms->nr,x,box,box_space,FALSE);
  write_conf(ftp2fn(efGRO,NFILE,fnm),cool_quote(),atoms,x,v,box);
  write_pdb_conf(opt2fn("-op",NFILE,fnm),atoms,x,box,FALSE);

  thanx(stdout);
  
  return 0;
}
