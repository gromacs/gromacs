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
 * Great Red Owns Many ACres of Sand 
 */
static char *SRCID_gro2top_c = "$Id$";

#include <math.h>
#include "smalloc.h"
#include "typedefs.h"
#include "macros.h"
#include "string2.h"
#include "confio.h"
#include "vec.h"
#include "statutil.h"
#include "copyrite.h"

#define  TOLERANCE  0.02

typedef struct t_nb {
  int     naid;
  atom_id aid[10];
  int     nnb;
  struct t_nb **nb;
} t_nb;

typedef struct {
  char  atomname[2];
  real  dist;
} t_bdist;

typedef struct {
  atom_id one;
  atom_id two;
} t_bonds;

typedef struct {
  atom_id one;
  atom_id two;
  atom_id three;
} t_angles;

typedef struct {
  atom_id one;
  atom_id two;
  atom_id three;
  atom_id four;
} t_impropers;

typedef struct {
  atom_id one;
  atom_id two;
  atom_id three;
  atom_id four;
} t_dihedrals;


typedef struct {
  char resname[10];
  char atomname[10];
  char type[10];
} t_types;

typedef struct {
  char name[10];
} t_attype;

typedef struct {
  atom_id one;
  atom_id two;
} t_pairs;

real distance(rvec a,rvec b)
{
  return sqrt(sqr(a[XX]-b[XX])+sqr(a[YY]-b[YY])+sqr(a[ZZ]-b[ZZ]));
}

find_sons(int nbonds,t_bonds *bonds,t_nb *nb,int level,int *color)
{
  int b,i;
  atom_id this,parent;
  /* print message */
  /* initialise variables */
  this = nb->aid[nb->naid-1];
  if (nb->naid==1)
    parent=nb->aid[0];
  else
    parent=nb->aid[nb->naid-2];
  nb->nnb=0;
  snew((nb->nb),7);
  
  if (color[this]==1)
    return;

  if ( level != 0 ) {
    for(b=0;(b<nbonds);b++) {
      bool bFound=FALSE;
      atom_id new_atom;
      if ((bonds[b].one==this)&&
	  (bonds[b].two!=parent)) {
	new_atom = bonds[b].two;
	bFound = TRUE;
      } else if ((bonds[b].two==this)&&
		 (bonds[b].one!=parent)) {
	new_atom = bonds[b].one;
	bFound=TRUE;
      }

      if (bFound ==TRUE ) {
	int index;

	/* add new son to the list */
	nb->nnb++;
	index = nb->nnb-1;
	nb->nb[index]=(struct t_nb *)malloc(sizeof(struct t_nb));
	

	for(i=0;(i<nb->naid);i++)
	  nb->nb[index]->aid[i]=nb->aid[i];
	nb->nb[index]->aid[nb->naid]=new_atom;
	nb->nb[index]->naid=nb->naid+1;
		
      }  
    }

    
    for(b=0;(b<nb->nnb);b++) {
      find_sons(nbonds,bonds,nb->nb[b],level-1,color);
    }
  } else {

    /* haal begin en eind effecten weg */

    if (nb->aid[0]<nb->aid[nb->naid-1]) {
      for(i=0;(i<nb->naid);i++)
	printf(" %5d ",nb->aid[i]+1);
      printf("\n");
    }
  }

}

int gendihedrals(t_atoms *atoms,int nbonds,t_bonds *bonds,t_dihedrals **dihedrals)
{
  int ndihedrals=0;
  int n,i,b;
  
  int *color;
  int *coordno;

  fprintf(stderr,"\n\n\n *** Generetaing Dihedrals ***\n");
  
  dihedrals = NULL;

  snew(coordno,atoms->nr);

  snew(color,atoms->nr);
  for(n=0;(n<atoms->nr);n++)
    color[n]=0;
  
  for(n=0;(n<atoms->nr);n++) {
    coordno[n]=0;
    for(b=0;(b<nbonds);b++) {
      if (n==bonds[b].one)
	coordno[n]++;
      else if (n==bonds[b].two)
	coordno[n]++;
    }
  }

  for(n=0;(n<atoms->nr);n++) {
    int b;
    int nlist = 0;
    atom_id list[10];

    for(b=0;(b<nbonds);b++) {
      if (bonds[b].one==n)
	list[nlist++]=bonds[b].two;
      else if (bonds[b].two==n)
	list[nlist++]=bonds[b].one;
    }

    if ( nlist > 2 ) {
      int nlistold;
      nlistold = nlist;
      /* remove neighbours that have less than 2 coordnumms */
      for(i=nlistold-1;(i>0);i--) 
	if ((coordno[list[i]]==1)&&(nlist>2)) {
	  color[list[i]]=1;
	  nlist--;
	}
      
    }  
  }
  

  for(n=0;(n<atoms->nr);n++) {
    t_nb nb;
    nb.naid=1;
    nb.aid[0]=n;
    find_sons(nbonds,bonds,&nb,3,color);
  }
  fprintf(stderr,"\n\n\n");
  return ndihedrals;
}

int genangles(t_atoms *atoms,int nbonds,t_bonds *bonds,t_angles **angles)
{
  int i,b,n,nn;
  int nangles=0;

  *angles=NULL;

  fprintf(stderr,"\n\n\n *** Generating Angles ***\n");
  for(i=0;(i<atoms->nr);i++) {
    /* maak een neighbour list van atoom i */
    int     nnlist = 0;
    atom_id nlist[10];  
    for(b=0;(b<nbonds);b++) {
      if (bonds[b].one==i)
	nlist[nnlist++]=bonds[b].two;
      else if (bonds[b].two==i)
	nlist[nnlist++]=bonds[b].one;
    }
    for(n=0;(n<nnlist-1);n++) {
      for(nn=n+1;(nn<nnlist);nn++) {
	fprintf(stderr,"\rNew angle ; Total is %d ",++nangles);
	srenew(*angles,nangles);
	(*angles)[nangles-1].one=nlist[n];
	(*angles)[nangles-1].two=i;
	(*angles)[nangles-1].three=nlist[nn];
	
      }
    }
  }
  return nangles;
}

int genimpropers(t_atoms *atoms,int nbonds,t_bonds *bonds,t_impropers **impropers)
{
  int i,b,n,nn;
  int nimpropers=0;

  *impropers=NULL;

  fprintf(stderr,"\n\n\n *** Generating Impropers ***\n");
  for(i=0;(i<atoms->nr);i++) {
    /* maak een neighbour list van atoom i */
    int     nnlist = 0;
    atom_id nlist[10];  
    for(b=0;(b<nbonds);b++) {
      if (bonds[b].one==i)
	nlist[nnlist++]=bonds[b].two;
      else if (bonds[b].two==i)
	nlist[nnlist++]=bonds[b].one;
    }
    /* if we have 3 neighbours we have an improper */
    if ( nnlist == 3) {
      fprintf(stderr,"\rNew Improper Total is %d ",++nimpropers);
      srenew(*impropers,nimpropers);
      (*impropers)[nimpropers-1].one=i;
      (*impropers)[nimpropers-1].two=nlist[0];
      (*impropers)[nimpropers-1].three=nlist[1];
      (*impropers)[nimpropers-1].four=nlist[2];
    }
  }
  return nimpropers;
}

void gentypes(char *typedata,t_atoms *atoms,t_attype *attype)
{
  FILE *fp;
  char line[STRLEN];
  char name[5];
  int i,n;
  t_types *types;
  int     ntypes;
  bool    bFound;
  char thistype[10];
  static bSave = FALSE;
  
  fprintf(stderr,"\n\n\n *** Assigning Types ***\n");
  fprintf(stderr,"I assume that a certain atomname occurs only once\n");
  fprintf(stderr,"in each residue\n");

  /* read type database */
  fp = ffopen(typedata,"r");
  fgets(line,STRLEN,fp);
  sscanf(line,"%d",&ntypes);
  fprintf(stderr,"Number of types is %d\n",ntypes);
  snew(types,ntypes);
  for(i=0;(i<ntypes);i++) {
    fgets(line,STRLEN,fp);
    memcpy(name,line,5);
    sscanf(name,"%s",types[i].atomname);
    memcpy(name,line+5,5);
    sscanf(name,"%s",types[i].resname);
    memcpy(name,line+10,5);
    sscanf(name,"%s",types[i].type);
  }
  fclose(fp);
  fprintf(stderr,"Atomtypes read\n");

  for(i=0;(i<ntypes);i++)
    fprintf(stderr,"%5d %5s %5s %5s \n",i+1,types[i].atomname,types[i].resname,types[i].type);

  /* scan along atoms and prompt if a certain type could not
     be assigned */
  for(n=0;(n<atoms->nr);n++) {
    bFound=FALSE;
    for(i=0;((i<ntypes)&&(bFound==FALSE));i++) {
      if ((strcmp(*atoms->atomname[n],types[i].atomname)==0)&&
	  (strcmp(*atoms->resname[atoms->atom[i].resnr],types[i].resname)==0)) {
	strcpy(thistype,types[i].type);
	bFound=TRUE;
      }
      
      
    }
    if ( bFound==FALSE ) {
      bSave=TRUE;
      /* this atom type in this residue has not been found before */
      /* ask to set the type */
      fprintf(stderr,"This attype was not found in the residue type database\n");
      fprintf(stderr,"atomnumber    : %5d\n",n+1);
      fprintf(stderr,"atomname      : %5s\n",*atoms->atomname[n]);
      fprintf(stderr,"residuenumber : %5d\n",atoms->atom[n].resnr+1);
      fprintf(stderr,"residuename   : %5s\n",*atoms->resname[atoms->atom[n].resnr]);
      fprintf(stderr,"Type the type :");
      scanf("%s",&thistype);
      ntypes++;
      srenew(types,ntypes);
      strcpy(types[ntypes-1].atomname,*atoms->atomname[n]);
      strcpy(types[ntypes-1].resname,*atoms->resname[atoms->atom[n].resnr]);
      strcpy(types[ntypes-1].type,thistype);
    }
    strcpy(attype[n].name,thistype);
  }

  /* save the type database */
  if (bSave==TRUE) {
    fp = ffopen(typedata,"w");
    fprintf(fp,"%5d\n",ntypes);
    for(i=0;(i<ntypes);i++) {
      fprintf(fp,"%-5s",types[i].atomname);
      fprintf(fp,"%-5s",types[i].resname);
      fprintf(fp,"%-5s\n",types[i].type);
    }
    fclose(fp);
  }
}

int genpairs(t_atoms *atoms,int nbonds,t_bonds *bonds,t_pairs **pairs,t_attype *attype)
{
  static int npairs=0;
  int n,m,i;
  int *color;

#define BLACK  0
#define WHITE  1
#define RED    2
#define YELLOW 3
#define BLUE   4

  snew(color,atoms->nr);
  
  *pairs=NULL;
  /* */
  fprintf(stderr,"\n\n\n *** Searching Pairs ***\n");
  
  /* walk over all atoms */
  for(n=0;(n<atoms->nr-1);n++) { 

    /* color all atoms white */
    for(m=0;(m<atoms->nr);m++)
      color[m]=WHITE;
    color[n]=BLACK;
    
    /* find the next neighbours of the BLACK atoms and color them RED */
    for(i=0;(i<nbonds);i++) {
      if (n==bonds[i].one) {
	if (color[bonds[i].two]==WHITE)
	  color[bonds[i].two]=RED;
      } else if (n==bonds[i].two) {
	if (color[bonds[i].one]==WHITE)
	  color[bonds[i].one]=RED;
      }
    }
	
    /* find the neighbours of the RED atoms and color them BLUE */
    for(i=0;(i<nbonds);i++) {
      if (color[bonds[i].one]==RED) {
	if (color[bonds[i].two]==WHITE)
	  color[bonds[i].two]=BLUE;
      } else if (color[bonds[i].two]==RED) {
	if (color[bonds[i].one]==WHITE)
	  color[bonds[i].one]=BLUE;
      }
    }
    
    /* find the neighbours of the BLUE atoms and color them YELLOW */
    for(i=0;(i<nbonds);i++) {
      if (color[bonds[i].one]==BLUE) {
	if (color[bonds[i].two]==WHITE)
	  color[bonds[i].two]=YELLOW;
      } else if (color[bonds[i].two]==BLUE) {
	if (color[bonds[i].one]==WHITE)
	  color[bonds[i].one]=YELLOW;
      }
    }

    /* now assign all YELLOW and BLACK combinations */
    for(m=n+1;(m<atoms->nr);m++)
      if (color[m]==YELLOW) {
	fprintf(stderr,"\r Number of pairs = %d",++npairs);
	srenew(*pairs,npairs);
	(*pairs)[npairs-1].one=n;
	(*pairs)[npairs-1].two=m;
      }
  }
  
  
  return npairs;
}


int genbonds(char *bdistdata,t_atoms *atoms,rvec *x,t_bonds **bonds)
{
  FILE *fp;
  t_bdist *bdist;
  int i;
  char line[STRLEN];
  real dd;
  static int nbdist=0;
  int     nbonds=0;

  /* initialise */
  fprintf(stderr,"\n\n\n *** Searching Bonds ***\n");
  /* read bond distance list */
  fp=ffopen(bdistdata,"r");
  fgets(line,STRLEN,fp);
  sscanf(line,"%d",&nbdist);
  snew(bdist,nbdist);

  for(i=0;(i<nbdist);i++) {
    fgets(line,STRLEN,fp);
    bdist[i].atomname[0]=line[0];
    bdist[i].atomname[1]=line[1];
    sscanf(line+5,"%f",&dd);
    bdist[i].dist=dd;
  }
  fclose(fp);
  /* now the bond distance database is read */
  

  /* this is a check if everything is read correctly */
  fprintf(stderr,"Bond distance database read\n");
  
  *bonds=NULL;

  {
    int i,j,n;

    for(i=0;(i<atoms->nr-1);i++)
      for(j=i+1;(j<atoms->nr);j++) {
	real dist,compdist=-1.0;
	dist = distance(x[i],x[j]);
	
	/* look up the compare distance */
	for(n=0;(n<nbdist);n++) {
	  if ((*atoms->atomname[i][0]==bdist[n].atomname[0])&&
	      (*atoms->atomname[j][0]==bdist[n].atomname[1])) {
	    compdist=bdist[n].dist;
	    n = nbdist;   
	  } else if ((*atoms->atomname[i][0]==bdist[n].atomname[1])&&
		     (*atoms->atomname[j][0]==bdist[n].atomname[0])) {
	    compdist=bdist[n].dist;
	    n = nbdist;   
	  } 
	}

	/* if the compare distance is still -1 this means it is not found */
	/* now switch to interactive mode to get this distance */
	if ( compdist==-1.0) {
	  fprintf(stderr,"The Distance between the following atoms is not defined\n");
	  fprintf(stderr,"%5d%-5s  %5d%-5s\n",
		  i+1,*atoms->atomname[i],
		  j+1,*atoms->atomname[j]);
	  fprintf(stderr,"Enter the distance\n");
	  while ((compdist<0)||(compdist>0.5))
	    scanf("%f",&compdist);
	  
	  /* add this distance to the entry */
	  nbdist++;
	  srenew(bdist,nbdist);
	  bdist[nbdist-1].dist=compdist;
	  bdist[nbdist-1].atomname[0]=*atoms->atomname[i][0];
	  bdist[nbdist-1].atomname[1]=*atoms->atomname[j][0];
  	  
	  /* save the new bond list */
	  {
	    FILE *fp;
	    int i;
	    fp=ffopen(bdistdata,"w");
	    fprintf(fp,"%d\n",nbdist);
	    for(i=0;(i<nbdist);i++) {
	      fprintf(fp,"%c%c   %8.3f\n",
		      bdist[i].atomname[0],
		      bdist[i].atomname[1],
		      bdist[i].dist);
	    }
	    fclose(fp);
	  }
	}
	
	/* now there is a compare distance and a distance */

	/* compare them, and if they are within tolerance add a bond */
	if (fabs(compdist-dist)<TOLERANCE) {
	  fprintf(stderr,"tolerance %f",fabs(compdist-dist));
	  fprintf(stderr,"New bond ; Total is %d \n",++nbonds);
	  srenew(*bonds,nbonds);
	  (*bonds)[nbonds-1].one=i;
	  (*bonds)[nbonds-1].two=j;
	}
	

      }
  }
  fprintf(stderr,"\n");
  
}
  

int main(int argc, char *argv[])
{
  char    title[STRLEN];
  int     natoms;
  rvec    *x,*v;
  matrix  box;
  t_atoms atoms;

  t_attype    *attype;
  t_bonds     *bonds;
  t_angles    *angles;
  t_pairs     *pairs;
  t_impropers *impropers;
  t_dihedrals *dihedrals;
  
  int     ndihedrals=0;
  int     npairs=0;
  int     nangles=0;
  int     nbonds=0;
  int     nimpropers=0;

  t_filenm fnm[] = {
    efGRO, "-f", NULL, ffREAD,
  };

#define NFILE asize(fnm)

  CopyRight(stderr,argv[0]);
  parse_common_args(&argc,argv,0,NFILE,fnm,FALSE,NULL);

  /* read configuration */
  atoms.nr=0;
  atoms.nres=0;
  get_coordnum(fnm[0].fn,&natoms);
  snew(atoms.atomname,natoms);
  snew(atoms.resname,natoms);
  snew(atoms.atom,natoms);
  snew(x,natoms);
  snew(v,natoms);
  read_whole_conf(fnm[0].fn,title,&atoms,x,v,box);

  /* generate bonds */
  nbonds = genbonds("bdist.dat",&atoms,x,&bonds);

  /* define types for atoms */
  snew(attype,atoms.nr);
  gentypes("typedata.dat",&atoms,attype);

  /* generate pairs */
  npairs = genpairs(&atoms,nbonds,bonds,&pairs,attype);

  /* generate angles */
  nangles = genangles(&atoms,nbonds,bonds,&angles);

  /* generate impropers */
  nimpropers = genimpropers(&atoms,nbonds,bonds,&impropers);

  /* generate dihedrals */
  ndihedrals = gendihedrals(&atoms,nbonds,bonds,&dihedrals);


  thanx(stdout);
  
  return 0;
}


