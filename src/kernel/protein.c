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
#include <math.h>
#include "typedefs.h"
#include "smalloc.h"
#include "macros.h"
#include "sysstuff.h"
#include "string2.h"
#include "physics.h"
#include "symtab.h"
#include "vec.h"
#include "fatal.h"
#include "protein.h"

static const real dH=0.1;
static const real mH=1.008;
static const real Alfa=(DEG2RAD*109.5);

static void calc_h_res(char *restp,int nah,t_addh ah[],
		       int natoms,t_atom atom[],char **aname[],
		       int res_start,rvec x[])
{
  int          i,m;
  int          nh[3],na[3];
  t_addh      *ah0,key;
  t_add_block *ab;
  
  key.resname=restp;
  if ((ah0=(t_addh *)bsearch(&key,ah,nah,sizeof(key),compaddh))==NULL)
    return;

  for(i=0; (i<ah0->n_add); i++) {
    /* Nasty shortcut again */
    ab=&(ah0->ab[i]);

    for(m=0; (m<3); m++)
      na[m]=search_atom(ab->na[m],res_start,natoms,aname);
    for(m=0; (m<ab->nh); m++)
      nh[m]=na[0]+1+m;
    for(   ; (m<3); m++)
      nh[m]=-1;
    calc_h_pos(ab->tp,nh,na,dH,Alfa,x);
  }
}

void calc_h(int natom,t_atom atom[],char **ai[],char **res[],rvec x[],
	    int nah,t_addh ah[])
{
  int i,resnr;
  
  fprintf(stderr,"Calculating hydrogen positions...\n");
  if (natom > 0) {
    resnr=atom[0].resnr-1;
    for(i=0; (i<natom); i++) {
      if (atom[i].resnr != resnr) {
	resnr=atom[i].resnr;
	calc_h_res(*(res[resnr]),nah,ah,natom,atom,ai,i,x);
      }
    }
  }
}

t_seq *search_res(t_seq *seq,int i)
{
  t_seq *loop;

  for (loop=seq; (loop != NULL); loop=loop->next)
    if (loop->pdba.resnr == i)
      return loop;

  return NULL;
}


void gen_h(t_seq *seq,int nah,t_addh ah[])
{
  int     i,j,k,l;
  rvec    xH={0.0, 0.0, 0.0};
  t_addh  *ahptr;
  t_seq   *ptr,*rstart,*nal;
  t_add_block *ab;
  char    *after,resnm[20];
  int     answer,nres;
  char    **name;
  char    buf[12];
  t_atom  at;
  char    yes[STRLEN];

  fprintf(stderr,"Adding Hydrogens...\n");
  
  ptr=seq;
  for(nres=0; (ptr != NULL); ptr=ptr->next)
    nres=ptr->pdba.resnr+1;
  
  /* Loop over residues */
  for(i=j=0; (i<nres); i++) {

    if ((rstart=search_res(seq,i)) == NULL) {
      fprintf(stderr,"Residue number %d not found, skipping rest of pdb\n",i);
      return;
    }
    /* Search entry in the database */
    strcpy(resnm,seq->pdba.resnm);
    if ((ahptr=search_h_db(nah,ah,i+1,seq->pdba.resnm)) != NULL) {
      
      /* Loop over add blocks */
      for(k=0; (k<ahptr->n_add); k++) {
	/* Nasty shortcut */
	ab=&(ahptr->ab[k]);

	/* Insert atoms after 'after' */
	after=ab->na[0];

	/* Search this atom */
	ptr=rstart;
	if (strcasecmp(*ptr->pdba.anm,after)!=0) {
	  while((ptr->next != NULL) && 
		(strcasecmp(*ptr->next->pdba.anm,after)!=0)) {
	    ptr=ptr->next;
	  }
	  ptr=ptr->next;
	  if (ptr == NULL)
	    fatal_error(0,"Atom %s not found in residue %s-%d\n",
			after,rstart->pdba.resnm,i+1);
	}

	for(l=0; (l<ab->nh); l++) {
	  char       buf[12];
	  int        blen;
	  t_seq      *nal;
	  t_atom     at;
	  char       **name;
	  
	  set_at(&at,mH,0.0,0,i);
	  strcpy(buf,after);
	  buf[0]='H';
	  
	  if (ab->nh > 1) {
	    blen=strlen(buf);
	    buf[blen]='0'+ab->nh-l;
	    buf[blen+1]='\0';
	  }
	  name=put_symtab(symtab,buf);
	  nal=init_al(&at,xH,name);
	  nal->next=ptr->next;
	  ptr->next=nal;
	  j++;
	  seq->natom++;
	}
      }
    }
#ifdef DEBUG
    else 
      fprintf(stderr,"No H's on residue type: %s\n",*seq->resname[i]);
#endif
  }

  fprintf(stderr,"Add hydrogens to the N-terminus (Y/N) ");
  scanf("%s",yes);
  if ( (yes[0]=='Y') || (yes[0]=='y')) { 
    /* now add hydrogens at the N terminal part */
    printf("\nWhich N-Terminal  type do you want for residue 1\n"); 
    printf("0. Not Protonated ( No Charge , 2 Hydrogens ) \n");
    printf("1. Protonated ( Positive charge , 3 Hydrogens )\n");
    printf("Type a number:");
    fflush(stdout);
    scanf("%d",&answer);
    switch ( answer ) {
    case 1:
      /* add Hydrogen */
      set_at(&at,mH,0.0,0,0);
      strcpy(buf,"H");
      name = put_symtab(&(seq->symtab),buf);   
      nal=init_al(&at,xH,name);
      nal->next=seq->al;
      seq->al = nal;
      seq->natom++;
      
    case 0:
      /* add Hydrogen  */
      set_at(&at,mH,0.0,0,0);
      strcpy(buf,"H");
      name = put_symtab(&(seq->symtab),buf);   
      nal=init_al(&at,xH,name);
      nal->next=seq->al;
      seq->al = nal;
      seq->natom++;
      
      break;
    default:
      exit(0);
    }
  }    

  
  fprintf(stderr,"Add hydrogens to the C-terminus (Y/N) ");
  scanf("%s",yes);
  if ( (yes[0]=='Y') || (yes[0]=='y')) { 
    /* now add hydrogens at the C terminal part */ 
    printf("\nWhich C-Terminal  type do you want for residue %d\n",seq->nres); 
    printf("0. Protonated ( One Hydrogen atom, No charge )\n");
    printf("1. Not Protonated ( No hydrogen atoms, Negative charge ) \n");
    printf("Type a number:");
    scanf("%d",&answer);
    fflush(stdout);
    switch ( answer ) {
    case 0:
      set_at(&at,mH,0.0,0,seq->nres-1);
      strcpy(buf,"H");
      name = put_symtab(&(seq->symtab),buf);   
      append_al(seq->al,&at,xH,name);
      seq->natom++;
      break;
    case 1:
      break;
    default:
      exit(0);
    }
  }
}


