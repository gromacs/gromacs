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
static char *SRCID_ter_db_c = "$Id$";

#include "sysstuff.h"
#include "smalloc.h"
#include "typedefs.h"
#include "symtab.h"
#include "futil.h"
#include "resall.h"
#include "h_db.h"
#include "string2.h"
#include "strdb.h"
#include "fatal.h"
#include "ter_db.h"
#include "toputil.h"

/* use bonded types definitions in hackblock.h */
#define ekwRepl ebtsNR+1
#define ekwAdd  ebtsNR+2
#define ekwDel  ebtsNR+3
#define ekwNR   3
char *kw_names[ekwNR] = {
  "replace", "add", "delete" 
};

int find_kw(char *keyw)
{
  int i;
  
  for(i=0; i<ebtsNR; i++)
    if (strcasecmp(btsNames[i],keyw) == 0)
      return i;
  for(i=0; i<ekwNR; i++)
    if (strcasecmp(kw_names[i],keyw) == 0)
      return ebtsNR + 1 + i;
  
  return NOTSET;
}

#define FATAL() fatal_error(0,"Reading Termini Database: not enough items on line\n%s",line)

static void read_atom(char *line, t_atom *a, t_atomtype *atype, 
		      char nnew[], int *cgnr)
{
  int    i, n;
  char   type[12];
  double m, q;
  
  if ( (i=sscanf(line,"%s%s%lf%lf%n", nnew, type, &m, &q, &n)) != 4 ) 
    fatal_error(0,"Reading Termini Database: expected %d items of atom data in stead of %d on line\n%s", 4, i, line);
  a->m=m;
  a->q=q;
  a->type=at2type(type,atype);
  if ( sscanf(line+n,"%d", cgnr) != 1 )
    *cgnr = NOTSET;
}

int read_ter_db(char *FF,char ter,t_hackblock **tbptr,t_atomtype *atype)
{
  FILE       *in;
  char       inf[STRLEN],header[STRLEN],buf[STRLEN],line[STRLEN];
 t_hackblock *tb;
  int        i,j,n,ni,kwnr,nb,maxnb,nh;
  
  sprintf(inf,"%s-%c.tdb",FF,ter);
  in=libopen(inf);
  if (debug)
    fprintf(debug,"Opened %s\n",inf);
  
  tb=NULL;
  nb=-1;
  maxnb=0;
  kwnr=NOTSET;
  get_a_line(in,line,STRLEN);
  while (!feof(in)) {
    if (get_header(line,header)) {
      /* this is a new block, or a new keyword */
      kwnr=find_kw(header);
      
      if (kwnr == NOTSET) {
	nb++;
	/* here starts a new block */
	if ( nb >= maxnb ) {
	  maxnb+=100;
	  srenew(tb,maxnb);
	}
	clear_t_hackblock(&tb[nb]);
	tb[nb].name=strdup(header);
      }
    } else {
      if (nb < 0)
	fatal_error(0,"reading termini database: "
		    "directive expected before line:\n%s\n"
		    "This might be a file in an old format.",line);
      /* this is not a header, so it must be data */
      if (kwnr >= ebtsNR) {
	/* this is a hack: add/rename/delete atoms */
	/* make space for hacks */
	if (tb[nb].nhack >= tb[nb].maxhack) {
	  tb[nb].maxhack+=10;
	  srenew(tb[nb].hack, tb[nb].maxhack);
	}
	nh=tb[nb].nhack;
	clear_t_hack(&(tb[nb].hack[nh]));
	for(i=0; i<4; i++)
	  tb[nb].hack[nh].a[i]=NULL;
	tb[nb].nhack++;
	
	/* get data */
	n=0;
	if ( kwnr==ekwRepl || kwnr==ekwDel ) {
	  if (sscanf(line, "%s%n", buf, &n) != 1) 
	    fatal_error(0,"Reading Termini Database: "
			"expected atom name on line\n%s",line);
	  tb[nb].hack[nh].oname = strdup(buf);
	  /* we only replace or delete one atom at a time */
	  tb[nb].hack[nh].nr = 1;
	} else if ( kwnr==ekwAdd ) {
	  read_ab(line, inf, &(tb[nb].hack[nh]));
	  get_a_line(in, line, STRLEN);
	} else
	  fatal_error(0,"unimplemented keyword number %d (%s:%d)",
		      kwnr,__FILE__,__LINE__);
	if ( kwnr==ekwRepl || kwnr==ekwAdd ) {
	  snew(tb[nb].hack[nh].atom, 1);
	  read_atom(line+n, tb[nb].hack[nh].atom, atype, 
		    buf, &tb[nb].hack[nh].cgnr);
	  tb[nb].hack[nh].nname = strdup(buf);
	}
      } else if (kwnr >= 0 && kwnr < ebtsNR) {
	/* this is bonded data: bonds, angles, dihedrals or impropers */
	srenew(tb[nb].rb[kwnr].b,tb[nb].rb[kwnr].nb+1);
	n=0;
	for(j=0; j<btsNiatoms[kwnr]; j++) {
	  if ( sscanf(line+n, "%s%n", buf, &ni) == 1 )
	    tb[nb].rb[kwnr].b[tb[nb].rb[kwnr].nb].a[j] = strdup(buf);
	  else
	    fatal_error(0,"Reading Termini Database: expected %d atom names (found %d) on line\n%s", btsNiatoms[kwnr], j-1, line);
	  n+=ni;
	}
	for(   ; j<MAXATOMLIST; j++)
	  tb[nb].rb[kwnr].b[tb[nb].rb[kwnr].nb].a[j] = NULL;
	strcpy(buf, "");
	sscanf(line+n, "%s", buf);
	tb[nb].rb[kwnr].b[tb[nb].rb[kwnr].nb].s = strdup(buf);
	tb[nb].rb[kwnr].nb++;
      } else
	fatal_error(0,"Reading Termini Database: Expecting a header at line\n"
		    "%s",line);
    }
    get_a_line(in,line,STRLEN);
  }
  nb++;
  srenew(tb,nb);
  
  fclose(in);
  
  if (debug) print_ter_db(debug, nb, tb, atype);
  
  *tbptr=tb;
  return nb;
}

t_hackblock *choose_ter(int nb,t_hackblock tb[],char *title)
{
  int i,sel;
  
  if (nb == 1)
    return &(tb[0]);
    
  printf("%s\n",title);
  for(i=0; (i<nb); i++) 
    printf("%2d: %s\n",i,tb[i].name);
  do {
    fscanf(stdin,"%d",&sel);
  } while ((sel < 0) || (sel >= nb));
  
  return &(tb[sel]);
}

static void print_atom(FILE *out,t_atom *a,t_atomtype *atype,char *newnm)
{
  fprintf(out,"%s\t%s\t%g\t%g\n",
	  newnm,type2nm(a->type,atype),a->m,a->q);
}

void print_ter_db(FILE *out,int nb,t_hackblock tb[],t_atomtype *atype) 
{
  int i,j,k,bt,nrepl,nadd,ndel;
  
  for(i=0; (i<nb); i++) {
    fprintf(out,"[ %s ]\n",tb[i].name);
    
    /* first count: */
    nrepl=0;
    nadd=0;
    ndel=0;
    for(j=0; j<tb[i].nhack; j++) 
      if ( tb[i].hack[j].oname!=NULL && tb[i].hack[j].nname!=NULL )
	nrepl++;
      else if ( tb[i].hack[j].oname==NULL && tb[i].hack[j].nname!=NULL )
	nadd++;
      else if ( tb[i].hack[j].oname!=NULL && tb[i].hack[j].nname==NULL )
	ndel++;
      else if ( tb[i].hack[j].oname==NULL && tb[i].hack[j].nname==NULL )
	fatal_error(0,"invalid hack (%s) in termini database",tb[i].name);
    if (nrepl) {
      fprintf(out,"[ %s ]\n",kw_names[ekwRepl-ebtsNR-1]);
      for(j=0; j<tb[i].nhack; j++) 
	if ( tb[i].hack[j].oname!=NULL && tb[i].hack[j].nname!=NULL ) {
	  fprintf(out,"%s\t",tb[i].hack[j].oname);
	  print_atom(out,tb[i].hack[j].atom,atype,tb[i].hack[j].nname);
	}
    }
    if (nadd) {
      fprintf(out,"[ %s ]\n",kw_names[ekwAdd-ebtsNR-1]);
      for(j=0; j<tb[i].nhack; j++) 
	if ( tb[i].hack[j].oname==NULL && tb[i].hack[j].nname!=NULL ) {
	  print_ab(out,&(tb[i].hack[j]));
	  fprintf(out,"\t");
	  print_atom(out,tb[i].hack[j].atom,atype,tb[i].hack[j].nname);
	}
    }
    if (ndel) {
      fprintf(out,"[ %s ]\n",kw_names[ekwDel-ebtsNR-1]);
      for(j=0; j<tb[i].nhack; j++)
	if ( tb[i].hack[j].oname!=NULL && tb[i].hack[j].nname==NULL )
	  fprintf(out,"%s\n",tb[i].hack[j].oname);
    }
    for(bt=0; bt<ebtsNR; bt++)
      if (tb[i].rb[bt].nb) {
	fprintf(out,"[ %s ]\n", btsNames[bt]);
	for(j=0; j<tb[i].rb[bt].nb; j++) {
	  for(k=0; k<btsNiatoms[bt]; k++) 
	    fprintf(out,"%s%s",k?"\t":"",tb[i].rb[bt].b[j].a[k]);
	  if ( tb[i].rb[bt].b[j].s )
	    fprintf(out,"\t%s",tb[i].rb[bt].b[j].s);
	  fprintf(out,"\n");
	}
      }
    fprintf(out,"\n");
  }
}

#ifdef DBTDB
int main(int argc,char *argv[])
{
  t_symtab   tab;
  t_atomtype *atype;
  t_hackblock *tb,*seltb;
  int nb;
  
  open_symtab(&tab);
  atype=read_atype("ffgmx",&tab);
  nb=read_ter_db("ffgmx-c.tdb",&tb,atype);
  seltb=choose_ter(nb,tb,"What do you want ?");
  print_ter_db(stdout,1,seltb,atype);
  
  return 0; 
}
#endif
