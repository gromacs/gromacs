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
 * GROtesk MACabre and Sinister
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
#include "fatal.h"
#include "ter_db.h"
#include "toputil.h"

void check_kw(char *inf,char *x,char *keyw,char *file,int line)
{
  if (strcasecmp(x,keyw) != 0)
    fatal_error(0,"Reading Termini Database %s (source %s, line %d) looking for keyword %s, found %s",inf,file,line,keyw,x);
}

#define FATAL() fatal_error(0,"Reading Termini Database (source code line: %d)",__LINE__)
#define CHECK_KW(x,kw) check_kw(inf,x,kw,__FILE__,__LINE__)

static void read_atom(FILE *in,t_atom *a,t_atomtype *atype,char nnew[])
{
  char   type[12];
  double m,q;
  
  if (fscanf(in,"%s%s%lf%lf",nnew,type,&m,&q) != 4) FATAL();
  a->m=m;
  a->q=q;
  a->type=at2type(type,atype);
}

int read_ter_db(char *inf,t_hackblock **tbptr,t_atomtype *atype)
{
  FILE       *in;
  char       bname[124],nnew[24],oldnm[24],buf[24],keyw[128];
  t_hackblock *tb=NULL;
  int        i,j,nb=0;

  in=libopen(inf);
  if (debug)
    fprintf(debug,"Opened %s\n",inf);
    
  while (fscanf(in,"%s%s",keyw,bname) == 2) {
    if (debug) 
      fprintf(debug,"block name %s\n",bname);
    CHECK_KW(keyw,"BLOCK_NAME");
    srenew(tb,nb+1);
    
    /* Name of block */
    tb[nb].bname=strdup(bname);
    
    /* Number of replacements */
    if (fscanf(in,"%s%d",keyw,&(tb[nb].nreplace)) != 2) FATAL();
    CHECK_KW(keyw,"N_REPLACE");

    if (debug) 
      fprintf(debug,"N_REPLACE %d\n",tb[nb].nreplace);
    snew(tb[nb].nm_repl,tb[nb].nreplace);
    snew(tb[nb].new_nm,tb[nb].nreplace);
    snew(tb[nb].repl_by,tb[nb].nreplace);
    for(i=0; (i<tb[nb].nreplace); i++) {
      if (fscanf(in,"%s",oldnm) != 1) FATAL();
      tb[nb].nm_repl[i]=strdup(oldnm);
      read_atom(in,&(tb[nb].repl_by[i]),atype,nnew);
      tb[nb].new_nm[i]=strdup(nnew);
    }
    /* Number of additions */
    if (fscanf(in,"%s%d",keyw,&(tb[nb].nadd)) != 2) FATAL();
    CHECK_KW(keyw,"N_ADD");
    if (debug) 
      fprintf(debug,"N_ADD %d\n",tb[nb].nadd);
    snew(tb[nb].ab,tb[nb].nadd);
    snew(tb[nb].adder,tb[nb].nadd);
    snew(tb[nb].add_nm,tb[nb].nadd);
    for(i=0; (i<tb[nb].nadd); i++) {
      read_ab(in,inf,&(tb[nb].ab[i]));
      read_atom(in,&(tb[nb].adder[i]),atype,nnew);
      tb[nb].add_nm[i]=strdup(nnew);
    }
    /* Number of impropers */
    if (fscanf(in,"%s%d",keyw,&(tb[nb].nidih)) != 2) FATAL();
    CHECK_KW(keyw,"N_IDIH");
    if (debug) 
      fprintf(debug,"N_IDIH %d\n",tb[nb].nidih);
    snew(tb[nb].idih,tb[nb].nidih);
    for(i=0; (i<tb[nb].nidih); i++) {
      for(j=0; (j<4); j++) {
	if (fscanf(in,"%s",buf) != 1) FATAL();
	tb[nb].idih[i].ai[j]=strdup(buf);
      }
    }
    /* Number of propers */
    if (fscanf(in,"%s%d",keyw,&(tb[nb].ndih)) != 2) FATAL();
    CHECK_KW(keyw,"N_DIH");
    if (debug) 
      fprintf(debug,"N_DIH %d\n",tb[nb].ndih);
    snew(tb[nb].dih,tb[nb].ndih);
    for(i=0; (i<tb[nb].ndih); i++) {
      for(j=0; (j<4); j++) {
	if (fscanf(in,"%s",buf) != 1) FATAL();
	tb[nb].dih[i].ai[j]=strdup(buf);
      }
    }
    /* Number of delete atoms! */
    if (fscanf(in,"%s%d",keyw,&(tb[nb].ndel)) != 2) FATAL();
    CHECK_KW(keyw,"N_DELETE");
    if (debug) 
      fprintf(debug,"N_DELETE %d\n",tb[nb].ndel);
    snew(tb[nb].nm_del,tb[nb].ndel);
    for(i=0; (i<tb[nb].ndel); i++) {
      if (fscanf(in,"%s",nnew) != 1) FATAL();
      tb[nb].nm_del[i]=strdup(nnew);
    }
    nb++;
  } 
  
  fclose(in);

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
    printf("%2d: %s\n",i,tb[i].bname);
  do {
    fscanf(stdin,"%d",&sel);
  } while ((sel < 0) || (sel >= nb));
  
  return &(tb[sel]);
}

static void print_atom(FILE *out,t_atom *a,t_atomtype *atype,char *newnm)
{
  fprintf(out,"%6s  %6s  %10e  %10e\n",
	  newnm,type2nm(a->type,atype),a->m,a->q);
}

void print_ter_db(FILE *out,int nb,t_hackblock tb[],t_atomtype *atype) 
{
  int i,j,k;
  
  for(i=0; (i<nb); i++) {
    fprintf(out,"%s\n",tb[i].bname);
    fprintf(out,"%d\n",tb[i].nreplace);
    for(j=0; (j<tb[i].nreplace); j++) {
      fprintf(out,"%6s  ",tb[i].nm_repl[j]);
      print_atom(out,&(tb[i].repl_by[j]),atype,tb[i].new_nm[j]);
    }
    fprintf(out,"%d\n",tb[i].nadd); 
    for(j=0; (j<tb[i].nadd); j++) {
      print_ab(out,&(tb[i].ab[j]));
      print_atom(out,&(tb[i].adder[j]),atype,tb[i].add_nm[j]);
    }
    fprintf(out,"%d\n",tb[i].nidih);
    for(j=0; (j<tb[i].nidih); j++) {
      for(k=0; (k<4); k++) 
	fprintf(out,"%6s  ",tb[i].idih[j].ai[k]);
      fprintf(out,"\n");
    }
    fprintf(out,"%d\n",tb[i].ndih);
    for(j=0; (j<tb[i].ndih); j++) {
      for(k=0; (k<4); k++) 
	fprintf(out,"%6s  ",tb[i].dih[j].ai[k]);
      fprintf(out,"\n");
    }
    fprintf(out,"%d\n",tb[i].ndel);
    for(j=0; (j<tb[i].ndel); j++)
      fprintf(out,"%6s  \n",tb[i].nm_del[j]);
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
