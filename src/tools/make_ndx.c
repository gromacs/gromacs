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
 * GROwing Monsters And Cloning Shrimps
 */
#include <ctype.h>
#include "sysstuff.h"
#include "strdb.h"
#include "futil.h"
#include "macros.h"
#include "string2.h"
#include "statutil.h"
#include "confio.h"
#include "assert.h"
#include "copyrite.h"
#include "typedefs.h"
#include "rdgroup.h"
#include "smalloc.h"
#include "index.h"

typedef enum { etOther, etProt, etDNA, erestNR } eRestp;
static  char *ResTP[erestNR] = { "OTHER", "PROTEIN", "DNA" };

static char **AminoAcids;   
static int NAA;

static char   *Sugars[]     = { "A", "T", "G", "C", "U" };
#define  NDNA asize(Sugars)

static void do_select(char *ndxfile,char *simpndx)
{
  FILE    *out;
  int     ngrps;
  int     i,j,nra;
  int     *isize;
  atom_id **index;
  char    **grpnames;

  do {
    printf("How many groups ? ");
    fflush(stdout);
  } while (scanf("%d",&ngrps) != 1);

  if (ngrps <= 0) 
    return;

  snew(isize,ngrps);
  snew(index,ngrps);
  snew(grpnames,ngrps);

  rd_index(ndxfile,ngrps,isize,index,grpnames);

  nra=0;
  for(i=0; (i<ngrps); i++)
    nra+=isize[i];
  out=ffopen(simpndx,"w");
  fprintf(out,"%8d  %8d\n",ngrps,nra);
  for(i=0; (i<ngrps); i++) {
    fprintf(out,"%10s  %8d",grpnames[i],isize[i]);
    for(j=0; (j<isize[i]); j++) {
      fprintf(out,"  %d",index[i][j]);
      if ((j % 15) == 0)
	fprintf(out,"\n");
    }
    fprintf(out,"\n");
  }
  fclose(out);
}

static bool bASK=TRUE;

static bool yn(void)
{
  char c;

  if (bASK) {
    do {
      c=toupper(fgetc(stdin));
    } while ((c != 'Y') && (c != 'N'));

    return (c == 'Y');
  }
  else
    return FALSE;
}

static t_atoms *get_atoms(char *infile)
{
  int     natoms;
  rvec    *x,*v;
  matrix  box;
  char    title[256];
  t_atoms *atoms;

  (void) get_coordnum(infile,&natoms);
  snew(x,natoms);
  snew(v,natoms);
  snew(atoms,1);
  snew(atoms->atom,natoms);
  snew(atoms->atomname,natoms);
  snew(atoms->resname,natoms);
  atoms->nr=0;
  atoms->nres=0;
  read_whole_conf(infile,title,atoms,x,v,box);
  printf("Read coordinate file. Title was:\n");
  printf("%s\n",title);

  sfree(x);
  sfree(v);

  return atoms;
}

static void p_status(int nres,eRestp restp[],int natres[])
{
  int i,j,ntp[erestNR];

  for(i=0; (i<erestNR); i++)
    ntp[i]=0;
  for(j=0; (j<nres); j++)
    ntp[restp[j]]++;
  
  for(i=0; (i<erestNR); i++) 
    printf("There are: %5d %10s residues\n",ntp[i],ResTP[i]);
}

int main(int argc,char *argv[])
{
  static char *desc[] = {
    "Indexfiles are necessary for almost every gromacs program. The",
    "gromacs preprocessor grompp needs it to determine the",
    "[IT]group number[it] of each particle, analysis tools use it to",
    "select groups of atoms to perform an analysis on,",
    "eg. to calculate a diffusion constant of a number of particles",
    "in a larger system. By using an indexfile to list the particles",
    "programs do not need to interpret the data that they get, ie.",
    "programs do not have to select which particles are the oxygen atoms",
    "and which are the hydrogens. This setup makes the programs more",
    "generally applicable and faster.[PAR]",
    "Because it is general not easy, and no fun to type in numbers,",
    "and more important, the probability of typing an incorrect number",
    "is proportional to the number of digits,",
    "the gromacs software contains a program to generate an indexfile",
    "called make_ndx. The program reads a coordinate file",
    "and uses atom names and numbers, as well as residue names and numbers",
    "to build up groups of atoms. The make_ndx program [IT]does[it]",
    "interpret the names, it does know things, especially about proteins.",
    "Proteins are split in groups Ca, Backbone, MainChain, MainChain+H,",
    "MainChain+Cb and Protein (containing all atoms). At that each of these",
    "categories can be split in residues. Non protein residues are grouped",
    "together, eg. SOL, Na, Cl, DMSO. Each of these groups can be split into",
    "atoms to be able to perfrom analysis on subgroups.[PAR]",
    "Normally make_ndx creates one indexfile in which each atom may",
    "be incorporated many times (eg. in sidechain, protein, system).",
    "make_ndx has an option",
    "to make a selection out of all the groups that are generated",
    "(eg. Protein and SOL) and write that to a separate indexfile",
    "(option -i indexfile)."
  };
  int      i;
  char     title[STRLEN];
  t_atoms  atoms;
  rvec     *x,*v;
  matrix   box;
  t_block  *block;
  char     **gnames;
  t_filenm fnm[] = {
    { efSTX, "-f", NULL,     ffREAD  },
    { efNDX, "-o", NULL,     ffWRITE },
    { efNDX, "-i", "simple", ffOPTWR }
  };
#define NFILE asize(fnm)
  t_pargs pa[] = {
    { "-ask", FALSE, etBOOL, &bASK, "ask about splitting groups." }
  };
  
  CopyRight(stdout,argv[0]);
  
  parse_common_args(&argc,argv,0,FALSE,NFILE,fnm,asize(pa),pa,asize(desc),
		    desc,0,NULL);
  
  NAA = get_strings("aminoacids.dat",&AminoAcids);

  block = new_block();
  snew(gnames,1);

  get_stx_coordnum(ftp2fn(efSTX,NFILE,fnm),&(atoms.nr));
  snew(x,atoms.nr);
  snew(v,atoms.nr);
  snew(atoms.resname,atoms.nr);
  snew(atoms.atom,atoms.nr);
  snew(atoms.atomname,atoms.nr);
  fprintf(stderr,"\nReading structure file\n");
  read_stx_conf(ftp2fn(efSTX,NFILE,fnm),title,&atoms,x,v,box);

  analyse(&atoms,block,&gnames,bASK,TRUE);

  printf("The following groups have been found:\n");
  for(i=0; (i<block->nr); i++)
    printf("%-20s: %5d atoms\n",gnames[i],
	   block->index[i+1]-block->index[i]);
  write_index(opt2fn("-o",NFILE,fnm),block,gnames);

  if (opt2bSet("-i",NFILE,fnm))
    do_select(opt2fn("-o",NFILE,fnm),opt2fn("-i",NFILE,fnm));
  
  thanx(stdout);
    
  return 0;
}
