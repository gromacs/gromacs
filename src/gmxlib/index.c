/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 1.6
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
 * Grunge ROck MAChoS
 */
static char *SRCID_index_c = "$Id$";

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
#include "smalloc.h"

typedef enum { etOther, etProt, etDNA, erestNR } eRestp;
static  char *ResTP[erestNR] = { "OTHER", "PROTEIN", "DNA" };

static char **AminoAcids;   
static int NAA;

static char   *Sugars[]     = { "A", "T", "G", "C", "U" };
#define  NDNA asize(Sugars)

static bool yn(bool bASK)
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

t_block *new_block(void)
{
  t_block *block;

  snew(block,1);
  snew(block->index,1);

  return block;
}

void write_index(char *outf, t_block *b,char **gnames)
{
  FILE *out;
  int  i,j,k;

  out=ffopen(outf,"w");
  fprintf(out,"%5d  %5d\n",b->nr,b->nra);
  for(i=0; (i<b->nr); i++) {
    fprintf(out,"%s  %5d\n",gnames[i],b->index[i+1]-b->index[i]);
    for(k=0,j=b->index[i]; j<b->index[i+1]; j++,k++) {
      fprintf(out,"%5d  ",b->a[j]);
      if ((k % 15) == 14)
	fprintf(out,"\n");
    }
    fprintf(out,"\n");
  }
  fclose(out);
}

void add_grp(t_block *b,char ***gnames,int nra,atom_id a[],char *name)
{
  int i;


  srenew(b->index,b->nr+2);
  srenew(*gnames,b->nr+1);
  (*gnames)[b->nr]=strdup(name);

  srenew(b->a,b->nra+nra);
  for(i=0; (i<nra); i++)
    b->a[b->nra++]=a[i];
  b->nr++;
  b->index[b->nr]=b->nra;
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

static void p_status(int nres,eRestp restp[],int natres[],bool bVerb)
{
  int i,j,ntp[erestNR];

  for(i=0; (i<erestNR); i++)
    ntp[i]=0;
  for(j=0; (j<nres); j++)
    ntp[restp[j]]++;
  
  if (bVerb)
    for(i=0; (i<erestNR); i++) 
      printf("There are: %5d %10s residues\n",ntp[i],ResTP[i]);
}

atom_id *mk_aid(t_atoms *atoms,eRestp restp[],eRestp res,int *nra,
		bool bTrue)
/* Make an array of atom_ids for all atoms with:
 * (restp[i] == res) == bTrue
 */
{
  atom_id *a;
  int     i;

  snew(a,atoms->nr);
  *nra=0;
  for(i=0; (i<atoms->nr); i++) 
    if ((restp[atoms->atom[i].resnr] == res) == bTrue)
      a[(*nra)++]=i;
  
  return a;
}

static void analyse_other(int nres,eRestp Restp[],int natres[],t_atoms *atoms,
			  t_block *gb,char ***gn,bool bASK,bool bVerb)
{
  char **restp=NULL;
  char **attp=NULL;
  char *rname,*aname;
  atom_id *other_ndx,*aid,*aaid;
  int  i,j,k,l,resnr,naid,naaid,natp,nrestp=0;
  
  for(i=0; (i<nres); i++)
    if (Restp[i] == etOther)
      break;
  if (i==gb->nr)
    return; /* No Others */
  
  if (bVerb)
    printf("Analysing Other...\n");
  snew(other_ndx,atoms->nr);
  for(k=0; (k<atoms->nr); k++) {
    resnr=atoms->atom[k].resnr;
    rname=*atoms->resname[resnr];
    if (Restp[resnr] ==  etOther) {
      for(l=0; (l<nrestp); l++)
        if (strcmp(restp[l],rname) == 0)
          break;
      if (l==nrestp) {
        srenew(restp,++nrestp);
        restp[nrestp-1]=strdup(rname);
      }
    }
  }
  for(i=0; (i<nrestp); i++) {
    snew(aid,atoms->nr);
    naid=0;
    for(j=0; (j<atoms->nr); j++) {
      rname=*atoms->resname[atoms->atom[j].resnr];
      if (strcmp(restp[i],rname) == 0) 
        aid[naid++] = j;
    }
    add_grp(gb,gn,naid,aid,restp[i]);
    if (bASK) {
      printf("split %s into atoms (y/n) ? ",restp[i]);
      fflush(stdout);
      if (yn(bASK)) {
	natp=0;
	for(k=0; (k<naid); k++) {
	  aname=*atoms->atomname[aid[k]];
	  for(l=0; (l<natp); l++)
	    if (strcmp(aname,attp[l]) == 0)
	      break;
	  if (l == natp) {
	    srenew(attp,++natp);
	    attp[natp-1]=aname;
	  }
	}
	if (natp > 1) {
	  for(l=0; (l<natp); l++) {
	    snew(aaid,naid);
	    naaid=0;
	    for(k=0; (k<naid); k++) {
	      aname=*atoms->atomname[aid[k]];
	      if (strcmp(aname,attp[l])==0) 
		aaid[naaid++]=aid[k];
	    }
	    add_grp(gb,gn,naaid,aaid,attp[l]);
	    sfree(aaid);
	  }
	}
	sfree(attp);
	attp=NULL;
      }
      sfree(aid);
    }
  }
  sfree(other_ndx);
}

static void analyse_prot(int nres,eRestp restp[],int natres[],t_atoms *atoms,
			 t_block *gb,char ***gn,bool bASK,bool bVerb)
{
  static char *pnoh[]    = { "H" };
  static char *calpha[]  = { "CA" };
  static char *bb[]      = { "N", "CA", "C" };
  static char *mc[]      = { "N", "CA", "C", "O" };
  static char *mcb[]     = { "N", "CA", "CB", "C", "O" };
  static char *mch[]     = { "N", "CA", "C", "O" , "H" };
  static char *sc[]      = { "N", "CA", "C", "O" , "H" };
  static char *scnoh[]   = { "N", "CA", "C", "O" , "H" };
  static char **chains[] = { NULL, pnoh, calpha, bb, mc, mcb ,mch, sc, scnoh};

#define NCH asize(chains)

  static int       sizes[NCH] = { 0, asize(pnoh), asize(calpha), asize(bb), asize(mc), asize(mcb), asize(mch), asize(sc), asize(scnoh)};
  static char   *ch_name[NCH] = { "Protein", "Protein-H", "C-Alpha", "Backbone", "MainChain", "MainChain+Cb", "MainChain+H", "SideChain", "SideChain-H" };
  static bool complement[NCH] = { 1, 1, 0, 0, 0, 0, 0, 1, 1};
  static int  wholename[NCH]  = {-1, 0,-1,-1,-1,-1,-1,-1, 4};
  /* the index in wholename gives the first item in the arrays of 
   * atomtypes that should be tested with 'strncmp' in stead of
   * strcasecmp, or -1 if all items should be tested with strcasecmp
   */

  int     i,n,j;
  atom_id *aid;
  int     nra,nnpres,npres;
  bool    match;
  char    ndx_name[STRLEN],*atnm;

  if (bVerb)
    printf("Analysing Protein...\n");
  snew(aid,atoms->nr);

  /* calculate the number of protein residues */
  npres=0;
  for(i=0; (i<atoms->nres); i++)
    if (restp[i] == etProt)
      npres++;

  /* find matching or complement atoms */
  for(i=0; (i<NCH); i++) {
    nra=0;
    for(n=0; (n<atoms->nr); n++) {
      if (restp[atoms->atom[n].resnr] == etProt) {
	for(j=0,match=0;(j<sizes[i]); j++) {
	  /* skip digits at beginning of atomname, e.g. 1H */
	  atnm=*atoms->atomname[n];
	  while (isdigit(atnm[0]))
	    atnm++;
	  if ( (wholename[i]==-1) || (j<wholename[i]) ) {
	    if (strcasecmp(chains[i][j],atnm) == 0)
	      match=1;
	  } else {
	    if (strncasecmp(chains[i][j],atnm,strlen(chains[i][j])) == 0)
	      match=1;
	  }
	}
	if (match != complement[i])
	  aid[nra++]=n;
      }
    }
    add_grp(gb,gn,nra,aid,ch_name[i]); 
  }
  
  /* find matching or complement atoms*/
  if (bASK) {
    for(i=0; (i<NCH); i++) {
      printf("Split %12s into %5d residues (y/n) ? ",ch_name[i],npres);
      if (yn(bASK)) {
	int resnr;
	nra = 0;
	for(n=0;((atoms->atom[n].resnr<npres) && (n<atoms->nr));) {
	  resnr = atoms->atom[n].resnr;
	  for(;((atoms->atom[n].resnr==resnr) && (n<atoms->nr));n++) {
	    for(j=0,match=0;(j<sizes[i]); j++) 
	      if (strcasecmp(chains[i][j],*atoms->atomname[n]) == 0)
		match=1;
	    if (match!=complement[i])
	      aid[nra++]=n;
	  }
	  /* copy the residuename to the tail of the groupname */
	  if (nra > 0) {
	    sprintf(ndx_name,"%s_%s%d",
		    ch_name[i],*atoms->resname[resnr],resnr+1);
	    add_grp(gb,gn,nra,aid,ndx_name);
	    nra = 0;
	  }
	}
      } 
    }
    printf("Make group with sidechain and C=O swapped (y/n) ? ");
    if (yn(bASK)) {
      /* Make swap sidechain C=O index */
      {
	int resnr,hold;
	nra = 0;
	for(n=0;((atoms->atom[n].resnr<npres) && (n<atoms->nr));) {
	  resnr = atoms->atom[n].resnr;
	  hold  = -1;
	  for(;((atoms->atom[n].resnr==resnr) && (n<atoms->nr));n++) {
	    if (strcmp("CA",*atoms->atomname[n]) == 0) {
	      aid[nra++]=n;
	      hold=nra;
	      nra+=2;
	    }
	    else if (strcmp("C",*atoms->atomname[n]) == 0) {
	      assert(hold != -1);
	      aid[hold]=n;
	    }
	    else if (strcmp("O",*atoms->atomname[n]) == 0) {
	      assert(hold != -1);
	      aid[hold+1]=n;
	    }
	    else if (strcmp("O1",*atoms->atomname[n]) == 0) {
	      assert(hold != -1);
	      aid[hold+1]=n;
	    }
	    else 
	      aid[nra++]=n;
	  }
	}
	/* copy the residuename to the tail of the groupname */
	if (nra > 0) {
	  add_grp(gb,gn,nra,aid,"SwapSC-CO");
	  nra = 0;
	}
      } 
    }
  }
  sfree(aid);
}

static void analyse_dna(int nres,eRestp restp[],int natres[],t_atoms *atoms,
			t_block *gb,char ***gn,bool bASK,bool bVerb)
{
  if (bVerb)
    printf("Analysing DNA... (not really)\n");
}

void analyse(t_atoms *atoms,t_block *gb,char ***gn,bool bASK,bool bVerb)
{
  eRestp  *restp;
  int     *natres;
  char    *resnm;
  atom_id *aid;
  int     nra;
  int     i,j;

  NAA = get_strings("aminoacids.dat",&AminoAcids);

  if (bVerb)
    printf("Analysing residue names:\n");
  snew(restp,atoms->nres);
  aid=mk_aid(atoms,restp,etOther,&nra,TRUE);
  add_grp(gb,gn,nra,aid,"System"); 
  sfree(aid);

  for(i=0; (i<atoms->nres); i++) {
    resnm=*atoms->resname[i];
    if (restp[i] == etOther)
      for(j=0; (j<NAA); j++) {
	if (strcasecmp(AminoAcids[j],resnm) == 0)
	  restp[i] = etProt;
      }
    if (restp[i] == etOther)
      for(j=0; (j<NDNA);  j++) {
	if (strcasecmp(Sugars[j],resnm) == 0)
	  restp[i] = etDNA;
      }
  }
  snew(natres,atoms->nres);
  for(i=0; (i<atoms->nr); i++)
    natres[atoms->atom[i].resnr]++;
  p_status(atoms->nres,restp,natres,bVerb);

  /* Protein */
  aid=mk_aid(atoms,restp,etProt,&nra,TRUE);
  if (nra > 0) 
    analyse_prot(atoms->nres,restp,natres,atoms,gb,gn,bASK,bVerb);
  
  sfree(aid);

  /* Non-Protein */
  aid=mk_aid(atoms,restp,etProt,&nra,FALSE);
  if ((nra > 0) && (nra < atoms->nr))
    add_grp(gb,gn,nra,aid,"Non-Protein"); 
  sfree(aid);

  /* DNA */
  aid=mk_aid(atoms,restp,etDNA,&nra,TRUE);
  if (nra > 0) {
    add_grp(gb,gn,nra,aid,"DNA"); 
    analyse_dna(atoms->nres,restp,natres,atoms,gb,gn,bASK,bVerb);
  }
  sfree(aid);

  /* Other */
  analyse_other(atoms->nres,restp,natres,atoms,gb,gn,bASK,bVerb);
  aid=mk_aid(atoms,restp,etOther,&nra,TRUE);
  if ((nra > 0) && (nra < atoms->nr))
    add_grp(gb,gn,nra,aid,"Other"); 
  sfree(aid);
}






