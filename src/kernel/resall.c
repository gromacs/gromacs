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
static char *SRCID_resall_c = "$Id$";

#include "sysstuff.h"
#include "assert.h"
#include "string2.h"
#include "strdb.h"
#include "futil.h"
#include "smalloc.h"
#include "fatal.h"
#include "symtab.h"
#include "macros.h"
#include "resall.h"
#include "pgutil.h"

#define PREM(str) fatal_error(0,"Premature EOF on %s",str)

t_atomtype *read_atype(char *adb,t_symtab *tab)
{
#define MAXAT 1000
  FILE       *in;
  char       aadb[STRLEN];
  char       buf[STRLEN],name[STRLEN];
  double     m;
  int        nratt;
  t_atomtype *at;

  sprintf(aadb,"%s.atp",adb);
  in=libopen(aadb);
  
  snew(at,1);
  snew(at->atom,MAXAT);
  snew(at->atomname,MAXAT);
  
  for(nratt=0; ; nratt++) {
    if (nratt >= MAXAT)
      fatal_error(0,"nratt >= MAXAT(%d). Increase the latter",MAXAT);
    if (feof(in))
      break;
    if (fgets2(buf,STRLEN,in) == NULL)
      break;
    if (sscanf(buf,"%s%lf",name,&m) != 2)
      break;
    set_at(&(at->atom[nratt]),m,0.0,nratt,0);
    at->atomname[nratt]=put_symtab(tab,name);
    fprintf(stderr,"\rAtomtype %d",nratt+1);
  }
  fclose(in);
  fprintf(stderr,"\n");
  at->nr=nratt;
  
  return at;
}

static void print_resatoms(FILE *out,t_atomtype *atype,t_restp *rtp)
{
  int j,tp;
  
  /* fprintf(out,"%5s\n",rtp->resname);
     fprintf(out,"%5d\n",rtp->natom); */
  fprintf(out,"[ %s ]\n",rtp->resname);
  fprintf(out," [ atoms ]\n");
  
  for(j=0; (j<rtp->natom); j++) {
    tp=rtp->atom[j].type;
    assert (tp >= 0);
    assert (tp < atype->nr);
    fprintf(out,"%6s%6s%8.3f%6d\n",
	    *(rtp->atomname[j]),*(atype->atomname[tp]),
	    rtp->atom[j].q,rtp->cgnr[j]);
  }
}

static bool read_atoms(FILE *in,char *line,
		       t_restp *r0,t_symtab *tab,t_atomtype *atype)
{
  int    i,j,cg,maxentries;
  char   buf[256],buf1[256];
  double q;

  /* Read Atoms */
  maxentries=0;
  r0->atom=     NULL;
  r0->atomname= NULL;
  r0->cgnr=     NULL;
  i=0;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) { 
    if (sscanf(line,"%s%s%lf%d",buf,buf1,&q,&cg) != 4)
      return FALSE;
    if (i>=maxentries) {
      maxentries+=100;
      srenew(r0->atom,     maxentries);
      srenew(r0->atomname, maxentries);
      srenew(r0->cgnr,     maxentries);
    }
    r0->atomname[i]=put_symtab(tab,buf);
    r0->atom[i].q=q;
    r0->cgnr[i]=cg;
    for(j=0; (j<atype->nr); j++)
      if (strcasecmp(buf1,*(atype->atomname[j])) == 0)
	break;
    if (j == atype->nr)
      fatal_error(0,"Atom type %s (residue %s) not found in atomtype "
		  "database",buf1,r0->resname);
    r0->atom[i].type=j;
    r0->atom[i].m=atype->atom[j].m;
    i++;
  }
  r0->natom=i;
  srenew(r0->atom,i);
  srenew(r0->atomname,i);
  srenew(r0->cgnr,i);

  return TRUE;
}
	
static bool read_atoms_old(FILE *in,char *line,
			   t_restp *r0,t_symtab *tab,t_atomtype *atype)
{
  int    i,j,cg,nat;
  char   buf[256],buf1[256];
  double q;

  if (sscanf(line,"%s",buf) != 1)
    return FALSE;
  r0->resname=strdup(buf);
  
  get_a_line(in,line,STRLEN);
  /* Read Atoms */
  if (sscanf(line,"%d",&nat) != 1)
    return FALSE;
  get_a_line(in,line,STRLEN);
  r0->natom=nat;
  snew(r0->atom,nat);
  snew(r0->atomname,nat);
  snew(r0->cgnr,nat);
  for(i=0; (i<nat); i++) {
    if (sscanf(line,"%s%s%lf%d",buf,buf1,&q,&cg) != 4)
      return FALSE;
    r0->atomname[i]=put_symtab(tab,buf);
    r0->atom[i].q=q;
    r0->cgnr[i]=cg;
    for(j=0; (j<atype->nr); j++)
      if (strcasecmp(buf1,*(atype->atomname[j])) == 0)
	break;
    if (j == atype->nr)
      fatal_error(0,"Atom type %s (residue %s) not found in atomtype "
		  "database",buf1,r0->resname);
    r0->atom[i].type=j;
    r0->atom[i].m=atype->atom[j].m;
    get_a_line(in,line,STRLEN);
  }
  return TRUE;
}

static void print_resbonds(FILE *out,t_resbond *rb)
{
  int j;

  /* fprintf(out,"%5d\n",rb->nb); */
  fprintf(out," [ bonds ]\n");

  for(j=0; (j<rb->nb); j++)
    fprintf(out,"%6s%6s\n",rb->rbond[j].ai,rb->rbond[j].aj);
}

static bool read_bonds(FILE *in,char *line,t_resbond *rb)
{
  char ai[12],aj[12];
  double c[MAXFORCEPARAM];
  int  i,j,n,maxentries;

  i=0;
  maxentries=0;
  rb->rbond=NULL;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) {
    if ((n=sscanf(line,"%s%s%lf%lf%lf%lf",ai,aj,&c[0],&c[1],&c[2],&c[3])) < 2)
      return FALSE;
    if (i>=maxentries) {
      maxentries+=100;
      srenew(rb->rbond,maxentries);
    }
    rb->rbond[i].ai=strdup(ai);
    rb->rbond[i].aj=strdup(aj);
    for (j=0; j<MAXFORCEPARAM; j++)
      if (j<n-2)
	rb->rbond[i].c[j]=c[j];
      else
	rb->rbond[i].c[j]=NOTSET;
    i++;
  }
  rb->nb=i;
  srenew(rb->rbond,i);

  return TRUE;
}

static bool read_bonds_old(FILE *in,char *line,t_resbond *rb)
{
  char ai[12],aj[12];
  double c[MAXFORCEPARAM];
  int  i,j,n,nb;

  if (sscanf(line,"%d",&nb) != 1)
    return FALSE;
  get_a_line(in,line,STRLEN);
  rb->nb=nb;
  snew(rb->rbond,nb);
  for(i=0; (i<nb); i++) {
    if ((n=sscanf(line,"%s%s%lf%lf%lf%lf",ai,aj,&c[0],&c[1],&c[2],&c[3])) < 2)
      return FALSE;
    rb->rbond[i].ai=strdup(ai);
    rb->rbond[i].aj=strdup(aj);
    for (j=0; j<MAXFORCEPARAM; j++)
      if (j<n-2)
	rb->rbond[i].c[j]=c[j];
      else
	rb->rbond[i].c[j]=NOTSET;
    get_a_line(in,line,STRLEN);
  }
  return TRUE;
}

static void print_resangs(FILE *out,t_resang *ra)
{
  int i,j;

  /* fprintf(out,"%5d\n",ra->na); */
  fprintf(out," [ angles ]\n");
  fprintf(out,";   ai    aj    ak            c0            c1\n");

  for(j=0; (j<ra->na); j++) {
    fprintf(out,"%6s%6s%6s",ra->rang[j].ai,ra->rang[j].aj,ra->rang[j].ak);
    for(i=0; (i<MAXFORCEPARAM && (ra->rang[j].c[i] != NOTSET)); i++)
      fprintf(out,"%14g",ra->rang[j].c[i]);
    fprintf(out,"\n");
  }
  
}

static bool read_angles(FILE *in,char *line,t_resang *ra)
{
  char ai[12],aj[12],ak[12];
  double c[MAXFORCEPARAM];
  int  i,j,n,maxentries;

  maxentries=0;
  ra->rang=NULL;
  i=0;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) {
    if ((n=sscanf(line,"%s%s%s%lf%lf%lf%lf",
		  ai,aj,ak,&c[0],&c[1],&c[2],&c[3])) < 3) 
      return FALSE;
    if (i>=maxentries) {
      maxentries+=100;
      srenew(ra->rang,maxentries);
    }
    ra->rang[i].ai=strdup(ai);
    ra->rang[i].aj=strdup(aj);
    ra->rang[i].ak=strdup(ak);
    for (j=0; j<MAXFORCEPARAM; j++)
      if (j<n-3)
	ra->rang[i].c[j]=c[j];
      else
	ra->rang[i].c[j]=NOTSET;
    i++;
  }
  ra->na=i;
  srenew(ra->rang,i);

  return TRUE ;
}

static void print_resdihs(FILE *out,t_resdih *rd)
{
  int i,j;

  fprintf(out," [ dihedrals ]\n");
  fprintf(out,";   ai    aj    ak    al            c0            c1            c2\n");
  for(j=0; (j<rd->nd); j++) {
    fprintf(out,"%6s%6s%6s%6s",rd->rdih[j].ai,rd->rdih[j].aj,
	                       rd->rdih[j].ak,rd->rdih[j].al);
    for(i=0; (i<MAXFORCEPARAM && (rd->rdih[j].c[i] != NOTSET)); i++)
      fprintf(out,"%14g",rd->rdih[j].c[i]);
    fprintf(out,"\n");
  }
  
}

static bool read_dihedrals(FILE *in,char *line,t_resdih *rd)
{
  int       i,j,n,maxentries;
  double    c[MAXFORCEPARAM];
  char      ai[4][12];

  maxentries=0;
  rd->rdih=NULL;
  i=0;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) {
    if ((n=sscanf(line,"%s%s%s%s%lf%lf%lf%lf",
		 ai[0],ai[1],ai[2],ai[3],&c[0],&c[1],&c[2],&c[3])) < 4)
      return FALSE;     
    if (i>=maxentries) {
      maxentries+=100;
      srenew(rd->rdih,maxentries);
    }
    rd->rdih[i].ai=strdup(ai[0]);
    rd->rdih[i].aj=strdup(ai[1]);
    rd->rdih[i].ak=strdup(ai[2]);
    rd->rdih[i].al=strdup(ai[3]);
    for (j=0; j<MAXFORCEPARAM; j++)
      if (j<n-4)
	rd->rdih[i].c[j]=c[j];
      else
	rd->rdih[i].c[j]=NOTSET;
    i++;
  }
  rd->nd=i;
  srenew(rd->rdih,i);

  return TRUE;
}

static void print_idihs(FILE *out,t_idihres *ires)
{
  int  j,k;
  
  /* fprintf(out,"%5d\n",ires->nidih); */
  fprintf(out," [ impropers ]\n");

  for(j=0; (j<ires->nidih); j++) {
    for(k=0; (k<4); k++)
      fprintf(out,"%6s",ires->idih[j].ai[k]);
    fprintf(out,"\n");
  }
}

static bool read_idihs(FILE *in,char *line,t_idihres *ires)
{
  int       i,j,n,maxentries;
  double    c[MAXFORCEPARAM];
  char      ai[4][12];

  maxentries=0;
  ires->idih=NULL;
  i=0;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) {
    if ((n=sscanf(line,"%s%s%s%s%lf%lf%lf%lf",
		 ai[0],ai[1],ai[2],ai[3],&c[0],&c[1],&c[2],&c[3])) < 4)
      return FALSE;     
    if (i>=maxentries) {
      maxentries+=100;
      srenew(ires->idih,maxentries);
    }
    for(j=0; (j<4); j++)
      ires->idih[i].ai[j]=strdup(ai[j]);
    for (j=0; j<MAXFORCEPARAM; j++)
      if (j<n-4)
	ires->idih[i].c[j]=c[j];
      else
	ires->idih[i].c[j]=NOTSET;
    i++;
  }
  ires->nidih=i;
  srenew(ires->idih,i);

  return TRUE;
}

static bool read_idihs_old(FILE *in,char *line,t_idihres *ires)
{
  int       i,j,n,nidih;
  double    c[MAXFORCEPARAM];
  char      ai[4][12];

  
  if (sscanf(line,"%d",&nidih) != 1)
    return FALSE;
  get_a_line(in,line,STRLEN);
  ires->nidih=nidih;
  snew(ires->idih,nidih);
  for(i=0; (i<nidih); i++) {
    if ((n=sscanf(line,"%s%s%s%s%lf%lf%lf%lf",
		 ai[0],ai[1],ai[2],ai[3],&c[0],&c[1],&c[2],&c[3])) < 4)
      return FALSE;     
    for(j=0; (j<4); j++)
      ires->idih[i].ai[j]=strdup(ai[j]);
    for (j=0; j<MAXFORCEPARAM; j++)
      if (j<n-4)
	ires->idih[i].c[j]=c[j];
      else
	ires->idih[i].c[j]=NOTSET;
    get_a_line(in,line,STRLEN);
  }
  return TRUE;
}

static void check_rtp(int nrtp,t_restp rtp[],char *libfn)
{
  int i;
  
  for(i=1; (i<nrtp); i++) {
    if (strcasecmp(rtp[i-1].resname,rtp[i].resname) == 0)
      fprintf(stderr,"WARNING double entry %s in file %s\n",rtp[i].resname,libfn);
  }
}

int read_resall(char       *ff,
		t_restp    **rtp,
		t_resbond  **rb,
		t_resang   **ra,
		t_resdih   **rd,
		t_idihres  **ires,
		t_atomtype *atype,
		t_symtab   *tab)
{
#define MAXRTP 1000
  FILE      *in;
  char      rrdb[STRLEN],line[STRLEN],*dum,header[STRLEN];
  int       nrtp,i;
  t_restp   *rrtp;
  t_resbond *rrbd;
  t_resang  *rran;
  t_resdih  *rrdi;
  t_idihres *rrid;
  bool      bNextResidue,bError;
  
  sprintf(rrdb,"%s.rtp",ff);
  in=libopen(rrdb);
  /* fscanf(in,"%d",&nrtp); */
  snew(rrtp,MAXRTP);
  snew(rrbd,MAXRTP);
  snew(rran,MAXRTP);
  snew(rrdi,MAXRTP);
  snew(rrid,MAXRTP);

  get_a_line(in,line,STRLEN);
  if (strchr(line,'[') == 0) {
    fprintf(stderr,"\n\n\tREADING .rtp FILE WITH OLD FORMAT\n\n"
	    "\tTO CONVERT TO NEW FORMAT USE THE HIDDEN OPTION -newrtp\n"
	    "\tWHICH WILL PRODUCE A FILE new.rtp\n\n\n");
    for(nrtp=0; ; nrtp++) {
      if (nrtp >= MAXRTP)
	fatal_error(0,"nrtp >= MAXRTP(%d). Increase the latter",MAXRTP);
      if (feof(in))
	break;
      if (!read_atoms_old(in,line,&(rrtp[nrtp]),tab,atype))
	fatal_error(0,"in .rtp file in atoms of residue %s:\n%s\n",
		    rrtp[nrtp].resname,line);
      if (!read_bonds_old(in,line,&(rrbd[nrtp])))
	fatal_error(0,"in .rtp file in bonds of residue %s:\n%s\n",
		    rrtp[nrtp].resname,line);
      rrbd[nrtp].resname=rrtp[nrtp].resname;
      rran[nrtp].resname=rrtp[nrtp].resname;
      rrdi[nrtp].resname=rrtp[nrtp].resname;
      if (!read_idihs_old(in,line,&(rrid[nrtp])))
	fatal_error(0,"in .rtp file in impropers of residue %s:\n%s\n",
		    rrtp[nrtp].resname,line);
      rrid[nrtp].resname=rrtp[nrtp].resname;
    
      if (debug) {
	fprintf(debug,"Residue %d",nrtp+1);
	fprintf(debug,"(%s): %d atoms,",rrtp[nrtp].resname,rrtp[nrtp].natom);
	fprintf(debug," %d bonds and",rrbd[nrtp].nb);
	fprintf(debug," %d angles and",rran[nrtp].na);
	fprintf(debug," %d dihedrals and",rrdi[nrtp].nd);
	fprintf(debug," %d impropers\n",rrid[nrtp].nidih);
      }
      fprintf(stderr,"\rResidue %d",nrtp+1);
    }
  }
  else {
    nrtp=0;
    while (!feof(in)) {
      if (nrtp >= MAXRTP)
	fatal_error(0,"nrtp >= MAXRTP(%d). Increase the latter",MAXRTP);
      if (!get_header(line,header))
	fatal_error(0,"in .rtp file at line:\n%s\n",line);
      rrtp[nrtp].resname=strdup(header);
  
      get_a_line(in,line,STRLEN);
      rrbd[nrtp].resname=rrtp[nrtp].resname;
      rran[nrtp].resname=rrtp[nrtp].resname;
      rrdi[nrtp].resname=rrtp[nrtp].resname;
      rrid[nrtp].resname=rrtp[nrtp].resname;
      bError=FALSE;
      bNextResidue=FALSE;
      while (!bNextResidue) {
	if (!get_header(line,header))
	  if (feof(in))
	    bNextResidue=TRUE;
	  else
	    bError=TRUE;
	else {
	  if (strncasecmp("atoms",header,5)==0) 
	    bError=!read_atoms(in,line,&(rrtp[nrtp]),tab,atype);
	  else if (strncasecmp("bonds",header,5)==0)
	    bError=!read_bonds(in,line,&(rrbd[nrtp]));
	  else if (strncasecmp("angles",header,6)==0) 
	    bError=!read_angles(in,line,&(rran[nrtp]));
	  else if (strncasecmp("dihedrals",header,6)==0) 
	    bError=!read_dihedrals(in,line,&(rrdi[nrtp]));
	  else if (strncasecmp("impropers",header,9)==0)
	    bError=!read_idihs(in,line,&(rrid[nrtp]));
	  else {
	    if (!feof(in) && !get_header(line,header))
	      bError=TRUE;
	    bNextResidue=TRUE;
	  }
	}
	if (bError)
	  fatal_error(0,"in .rtp file in residue %s at line:\n%s\n",
		      rrtp[nrtp].resname,line);
      }
      if (rrtp[nrtp].natom == 0)
	fatal_error(0,"No atoms found in .rtp file in residue %s\n",
		    rrtp[nrtp].resname);

      if (debug) {
	fprintf(debug,"Residue %d(%s): %d atoms, %d bonds, %d angles "
		"%d dihedrals and %d impropers\n",nrtp+1,rrtp[nrtp].resname,
		rrtp[nrtp].natom,rrbd[nrtp].nb,rran[nrtp].na,rrdi[nrtp].nd,
		rrid[nrtp].nidih);
      }
      nrtp++;
      fprintf(stderr,"\rResidue %d",nrtp);
    }
  }
  fclose(in);

  fprintf(stderr,"\nSorting it all out...\n");
  qsort(rrtp,nrtp,(size_t)sizeof(rrtp[0]),comprtp);
  qsort(rrbd,nrtp,(size_t)sizeof(rrbd[0]),comprb);
  qsort(rran,nrtp,(size_t)sizeof(rran[0]),comprang);
  qsort(rrdi,nrtp,(size_t)sizeof(rrdi[0]),comprdih);
  qsort(rrid,nrtp,(size_t)sizeof(rrid[0]),icomp);
  
  /* Check for consistency and doubles */
  check_rtp(nrtp,rrtp,rrdb);
  
  *rtp  = rrtp;
  *rb   = rrbd;
  *ra   = rran;
  *rd   = rrdi;
  *ires = rrid;
  
  return nrtp;
}

void print_resall(FILE *out,
		  int nrtp,
		  t_restp rtp[],
		  t_resbond rb[],
		  t_resang ra[],
		  t_resdih rd[],
		  t_idihres ires[],
		  t_atomtype *atype)
{
  int i;
  
  for(i=0; (i<nrtp); i++) {
    if (rtp[i].natom > 0) {
      print_resatoms(out,atype,&rtp[i]);
      if (rb[i].nb)
	print_resbonds(out,&rb[i]);
      if (ra[i].na)
	print_resangs(out,&ra[i]);
      if (rd[i].nd)
	print_resdihs(out,&rd[i]);
      if (ires[i].nidih)
	print_idihs(out,&(ires[i]));
      fprintf(out,"\n");
    }
  }
}

/************************************************************
 *
 *                  SEARCH   ROUTINES
 * 
 ***********************************************************/
int comprb(const void *a,const void *b)
{
  t_resbond *r1,*r2;

  r1=(t_resbond *)a;
  r2=(t_resbond *)b;

  return strcasecmp(r1->resname,r2->resname);
}

t_resbond *search_rb(char *key,int nrb,t_resbond rb[])
{
  t_resbond rbkey;
  
  rbkey.resname=key;
  return (t_resbond *)bsearch((void *)&rbkey,rb,nrb,
			      (size_t)sizeof(rbkey),comprb);
}

int comprang(const void *a,const void *b)
{
  t_resang *r1,*r2;

  r1=(t_resang *)a;
  r2=(t_resang *)b;

  return strcasecmp(r1->resname,r2->resname);
}

t_resang *search_rang(char *key,int nrang,t_resang rang[])
{
  t_resang rangkey;
  
  rangkey.resname=key;
  return (t_resang *)bsearch((void *)&rangkey,rang,nrang,
			      (size_t)sizeof(rangkey),comprang);
}

int comprdih(const void *a,const void *b)
{
  t_resdih *r1,*r2;

  r1=(t_resdih *)a;
  r2=(t_resdih *)b;

  return strcasecmp(r1->resname,r2->resname);
}

t_resdih *search_rdih(char *key,int ndih,t_resdih rdih[])
{
  t_resdih rdihkey;
  
  rdihkey.resname=key;
  return (t_resdih *)bsearch((void *)&rdihkey,rdih,ndih,
			      (size_t)sizeof(rdihkey),comprdih);
}

int comprtp(const void *a,const void *b)
{
  t_restp *ra,*rb;

  ra=(t_restp *)a;
  rb=(t_restp *)b;

  return strcasecmp(ra->resname,rb->resname);
}

int neq_str(char *a1,char *a2)
{
  int j,l;
  
  l=min((int)strlen(a1),(int)strlen(a2));
  for(j=0; (j<l) && (a1[j] == a2[j]); j++)
    ;
  
  return j;
}

t_restp *search_rtp(char *key,int nrtp,t_restp rtp[])
{
  int i,n,best,besti;

  besti=-1;
  best=1;
  for(i=0; (i<nrtp); i++) {
    n=neq_str(key,rtp[i].resname);
    if (n > best) {
      besti=i;
      best=n;
    }
  }
  if (besti == -1)
    fatal_error(0,"Residue '%s' not found in residue topology database\n",key);
  if (strlen(rtp[besti].resname) != strlen(key))
    fprintf(stderr,"Warning: '%s' not found in residue topology database, trying to use %s\n",
	    key,rtp[besti].resname);
	    
  return &rtp[besti];
}

int icomp(const void *a,const void *b)
{
  t_idihres *ra,*rb;

  ra=(t_idihres *)a;
  rb=(t_idihres *)b;
  
  return strcasecmp(ra->resname,rb->resname);
}

t_idihres *search_idih(char *key,int nrdh,t_idihres ires[])
{
  t_idihres rkey;

  rkey.resname=key;

  return (t_idihres *)bsearch((void *)&rkey,ires,nrdh,
			      (size_t)sizeof(rkey),icomp);
}

/************************************************************
 *
 *                  DUMMY DATABASE ROUTINES
 * 
 ***********************************************************/
static bool read_mass(FILE *in,char *line,t_dumblock *db)
{
  char name[12],type[12],ai[12];
  real c[MAXFORCEPARAM];
  int  i,j,n,maxentries,tp,nm;
  
  maxentries=0;
  db->mass=NULL;
  i=0;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) {
    if (sscanf(line,"%s%s%d%d",name,type,&nm,&tp) < 4)
      return FALSE;
    if (i>=maxentries) {
      maxentries+=100;
      srenew(db->mass,maxentries);
    }
    n=sscanf(line,"%*s%*s%*d%*d%g%g%g%g%g%g",
	     &c[0],&c[1],&c[2],&c[3],&c[4],&c[5]);
    for(j=n; (j<MAXFORCEPARAM); j++)
      c[j]=NOTSET;
    switch(tp) {
    case 1:
      if (n < 3)
	return FALSE;
      else
	if (n > 4) {
	  fprintf(stderr,"WARNING: Ignoring arguments 9..%d in .ddb database "
		  "file on line:\n%s\n",4+n,line);
	  for(j=4; (j<n); j++)
	    c[j]=NOTSET;
	}
      break;
    default:
      fatal_error(0,"Invalid dummy mass type (%d) in .ddb database file "
		  "(resall)\n",tp);
    }
    db->mass[i].mname=strdup(name);
    db->mass[i].mtype=strdup(type);
    db->mass[i].nm=nm;
    db->mass[i].tp=tp;
    for(j=0; (j<MAXFORCEPARAM); j++)
      db->mass[i].c[j]=c[j];
    i++;
  }
  db->nmass=i;
  srenew(db->mass,i);
  
  return TRUE ;
}

static bool read_dum(FILE *in,char *line,t_dumblock *db)
{
  char ai[12],aj[12],ak[12],al[12];
  real a,b,c;
  int  i,n,maxentries,tp;
  
  maxentries=0;
  db->dum=NULL;
  i=0;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) {
    if ((n=sscanf(line,"%s%s%s%s%d%g%g%g",ai,aj,ak,al,&tp,&a,&b,&c)) < 5) 
      return FALSE;
    if (i>=maxentries) {
      maxentries+=100;
      srenew(db->dum,maxentries);
    }
    db->dum[i].na[0]=strdup(ai);
    db->dum[i].na[1]=strdup(aj);
    db->dum[i].na[2]=strdup(ak);
    db->dum[i].na[3]=strdup(al);
    db->dum[i].tp=tp;
    if ( ( ( (tp==2) || (tp==5) ) && (n<7) ) || ( (tp==3) && (n<8) ) )
      return FALSE;
    if (n>=7) {
      db->dum[i].a = a;
      db->dum[i].b = b;
    } else {
      db->dum[i].a = NOTSET;
      db->dum[i].b = NOTSET;
    }
    if (n>=8) {
      db->dum[i].c = c;
    } else {
      db->dum[i].c = NOTSET;
    }
    i++;
  }
  db->ndum=i;
  srenew(db->dum,i);
  
  return TRUE ;
}

static bool read_angcon(FILE *in,char *line,t_dumblock *db)
{
  char ai[12],aj[12],ak[12];
  real c0,c1;
  int  i,k,maxentries;
  
  maxentries=0;
  db->ang=NULL;
  i=0;
  while (get_a_line(in,line,STRLEN) && (strchr(line,'[')==NULL)) {
    if (sscanf(line,"%s%s%s%g%g",ai,aj,ak,&c0,&c1) < 5) 
      return FALSE;
    if (i>=maxentries) {
      maxentries+=100;
      srenew(db->ang,maxentries);
    }
    db->ang[i].na[0]=strdup(ai);
    db->ang[i].na[1]=strdup(aj);
    db->ang[i].na[2]=strdup(ak);
    db->ang[i].c[0]=c0;
    db->ang[i].c[1]=c1;
    for (k=2; (k<MAXFORCEPARAM); k++)
      db->ang[i].c[k]=NOTSET;
    i++;
  }
  db->nang=i;
  srenew(db->ang,i);
  
  return TRUE ;
}

int read_dum_db(char *ff,t_dumblock **dbptr)
{
  FILE       *in;
  char       fn[STRLEN],line[STRLEN],*dum,header[STRLEN];
  int        nddb,maxddb,i,j;
  t_dumblock *ddb;
  bool       bNextBlock,bError;
  
  sprintf(fn,"%s.ddb",ff);
  in=libopen(fn);
  
  get_a_line(in,line,STRLEN);
  nddb=0;
  maxddb=100;
  snew(ddb,100);
  while (!feof(in)) {
    if (nddb >= maxddb) {
      fatal_error(0,"too much in .ddb file");
    }
    if (!get_header(line,header))
      fatal_error(0,"in .ddb file at line:\n%s",line);
    ddb[nddb].bname=strdup(header);
    get_a_line(in,line,STRLEN);
    bError=FALSE;
    bNextBlock=FALSE;
    while (!bNextBlock) {
      if (!get_header(line,header))
	if (feof(in))
	  bNextBlock=TRUE;
	else
	  bError=TRUE;
      else {
	if (strncasecmp("masses",header,6)==0) 
	  bError=!read_mass  (in,line,&(ddb[nddb]));
	else if (strncasecmp("dummies",header,7)==0)
	  bError=!read_dum   (in,line,&(ddb[nddb]));
	else if (strncasecmp("angles",header,6)==0) 
	  bError=!read_angcon(in,line,&(ddb[nddb]));
	else {
	  if (!feof(in) && !get_header(line,header))
	    bError=TRUE;
	  bNextBlock=TRUE;
	}
      }
    }
    if (bError)
      fatal_error(0,"in .ddb file in block %s at line:\n%s\n",
		  ddb[nddb].bname,line);
    if (ddb[nddb].nmass + ddb[nddb].ndum + ddb[nddb].nang == 0)
      fatal_error(0,"No data found in .ddb file in block %s\n",
		  ddb[nddb].bname);
    if (debug) {
      fprintf(debug,"Block %d (%s): %d masses, %d dummies and %d angles\n",
	      nddb+1, ddb[nddb].bname,
	      ddb[nddb].nmass, ddb[nddb].ndum, ddb[nddb].nang);
    }
    nddb++;
  }
  fclose(in);
  
  *dbptr=ddb;
  return nddb;
}

t_dumblock *search_dum_db(int nddb, t_dumblock *ddb, char resnm[])
{
  int i;
  
  for (i=0; (i<nddb); i++)
    if (strcmp(resnm,ddb[i].bname)==0)
      break;
  
  return &ddb[i];
}

void print_dum_db(FILE *out,int nddb,t_dumblock *ddb)
{
  int i,j,k,l;
  
  fprintf(out,"; Dummies Database\n");
  for (i=0; (i<nddb); i++) {
    fprintf(out,"[ %s ]\n",ddb[i].bname);
    if (ddb[i].nmass) {
      fprintf(out," [ masses ]\n");
      fprintf(out,"; %4s %5s %5s %5s %5s\n","name","type","n","type","from");
    }
    for(j=0; (j<ddb[i].nmass); j++) {
      fprintf(out," %5s %5s %5d %5d",ddb[i].mass[j].mname,
	      ddb[i].mass[j].mtype,ddb[i].mass[j].nm,ddb[i].mass[j].tp);
      fprintf(out,"\n");
    }
    if (ddb[i].ndum) {
      fprintf(out," [ dummies ]\n");
      fprintf(out,"; %4s %5s %5s %5s %5s %9s %9s %9s\n",
	      "atom","from","","","type","a","b","c");
    }
    for(j=0; (j<ddb[i].ndum); j++) {
      for(k=0; (k<4); k++)
	fprintf(out," %5s",ddb[i].dum[j].na[k]);
      fprintf(out," %5d",ddb[i].dum[j].tp);
      if (ddb[i].dum[j].a != NOTSET)
	fprintf(out," %9g %9g",ddb[i].dum[j].a,ddb[i].dum[j].b);
      if (ddb[i].dum[j].c != NOTSET)
	fprintf(out," %9g",ddb[i].dum[j].c);
      fprintf(out,"\n");
    }
    if (ddb[i].nang) fprintf(out," [ angles ]\n");
    for(j=0; (j<ddb[i].nang); j++) {
      for(k=0; (k<3); k++)
	fprintf(out," %5s",ddb[i].ang[j].na[k]);
      fprintf(out,"\n");
    }
  }
}

#ifdef DEBUG
int main(int argc,char *argv[])
{
  int        nrtp;
  t_atomtype *atype;
  t_restp    *rtp;
  t_resbond  *rb;
  t_resang   *ra;
  t_idihres  *id;
  t_symtab   tab;
    
  open_symtab(&tab);
  atype=read_atype("atomtype",&tab);
  nrtp=read_resall("residue",&rtp,&rb,&ra,&id,atype,&tab);
  close_symtab(&tab);
  
  print_resall(stdout,nrtp,rtp,rb,ra,id,atype);
  
  return 0;
}
#endif
