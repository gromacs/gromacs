#include <stdio.h>
#include "smalloc.h"
#include "strdb.h"
#include "futil.h"
#include "pdbio.h"
#include "vec.h"
#include "names.h"
#include "disco.h"
#include "pp2shift.h"

char *edc_names[edcNR+1] = { "NONBND", "BOND", "DISRE", NULL };

bool newres(int i,t_atom atom[])
{
  return ((i == 0) || (atom[i].resnr != atom[i-1].resnr));
}

t_correct *init_corr(int maxnit,int nstprint,int nbcheck,int nstranlist,
		     int ngrow,bool bExplicit,bool bChiral,bool bPep,
		     bool bDump,real lowdev)
{
  t_correct *c;
  
  snew(c,1);
  c->maxnit     = maxnit;
  c->nstprint   = nstprint;
  c->nbcheck    = nbcheck;
  c->nstranlist = nstranlist;
  c->bExplicit  = bExplicit;
  c->bChiral    = bChiral;
  c->bPep       = bPep;
  c->bDump      = bDump;
  c->lodev      = lowdev;
  c->maxdist    = 0;
  c->ndist      = 0;
  c->ngrow      = ngrow;
  
  return c;
}

void make_tags(t_correct *c,int natom)
{
  int i,n,ai;
  
  if (c->ndist == 0)
    fatal_error(0,"Can not make tags without distances. %s, %d",
		__FILE__,__LINE__);
		
  /* First sort the distances */
  
  /* Now make the tags */
  snew(c->tag,natom+1);
  ai = c->d[0].ai;
  n  = 0;
  for(i=0; (i<c->ndist); i++) {
    if (c->d[i].ai != ai) {
      /* New tag starts here, but skip over possible missing atoms */
      while ((n < natom) && (n<c->d[i].ai)) {
	n++;
	c->tag[n] = i;
	ai = c->d[i].ai;
      }
      /*     else
	     fatal_error(0,"Too many atoms, or distances not sorted");*/
    }
  }
  if (n < natom) {
    while (n < natom) {
      n++;
      c->tag[n] = i;
    }
  }
  else
    fatal_error(0,"Too many atoms, or distances not sorted");
  fprintf(stderr,"There are %d tags for %d atoms\n",n,natom);
  
  if (debug) 
    for(i=0; (i<=natom); i++)
      fprintf(debug,"tag[%5d] = %d\n",i,c->tag[i]);
}

void center_in_box(int natom,rvec xin[],matrix box,rvec xout[])
{
  int  i,m;
  rvec xcm,dx;
  
  /* Dump the optimization trajectory to an xtc file */
  /* Put the whole thing in the center of the box 1st */
  clear_rvec(xcm);
  for(i=0; (i<natom); i++) {
    copy_rvec(xin[i],xout[i]);
    rvec_inc(xcm,xin[i]);
  }
  for(m=0; (m<DIM); m++)
    dx[m] = 0.5*box[m][m]-xcm[m]/natom;
  for(i=0; (i<natom); i++) 
    rvec_inc(xout[i],dx);
}

void define_peptide_bonds(FILE *log,t_atoms *atoms,t_correct *c)
{
  int    i,npep,naa,nlist;
  char   **aa;
  t_dlist *dlist;
  
  naa   = get_strings("aminoacids.dat",&aa);
  dlist = mk_dlist(log,atoms,&nlist,TRUE,TRUE,FALSE,0,1,naa,aa);
  for(i=0; (i<naa); i++)
    sfree(aa[i]);
  sfree(aa);

  npep  = 0;
  snew(c->pepbond,nlist);
  snew(c->omega,nlist);
  for(i=0; (i<nlist); i++) {
    if (has_dihedral(edOmega,&dlist[i])) {
      c->pepbond[npep].ai = dlist[i].atm.minC;
      c->pepbond[npep].aj = dlist[i].atm.minO;
      c->pepbond[npep].ak = dlist[i].atm.N;
      c->pepbond[npep].al = dlist[i].atm.H;
      c->pepbond[npep].am = dlist[i].atm.Cn[1];
      npep++;
    }
  }
  c->npep = npep;
  if (debug)
    pr_dlist(debug,nlist,dlist,1.0);
  sfree(dlist);
  fprintf(log,"There are %d peptide bonds\n",npep);
}

static char *aname(t_atoms *atoms,int i)
{
  static char buf[32];
  
  sprintf(buf,"%4s%3d:%4s", 
	  *(atoms->resname[atoms->atom[i].resnr]),
	  atoms->atom[i].resnr+1,
	  *(atoms->atomname[i]));
  
  return buf;
}

int find_atom(t_atoms *atoms,int resnr,int j0,char *nm)
{
  int j,aa=-1;
  
  for(j=j0; (j < atoms->nr) && (atoms->atom[j].resnr == resnr); j++)
    if (strcmp(*atoms->atomname[j],nm) == 0) {
      aa = j;
      break;
    }
  return aa;
}

void define_impropers(FILE *log,t_atoms *atoms,t_correct *c)
{
  typedef struct {
    char *res,*aa[4];
  } t_impdef;
  
  t_impdef id[] = {
    { NULL,   "CA", "N", "C", "CB" },
    { NULL,   "N",  "CA", "H1",  "H3" },
    { "LYSH", "NZ", "CE", "HZ2", "HZ1" },
    { "LYSH", "NZ", "CE", "HZ1", "HZ3" },
    { "LEU",  "CG", "CD2", "CD1", "CB" },
    { "VAL",  "CB", "CG2", "CG1", "CA" },
    { "ILE",  "CB", "CG1", "CG2", "CA" },
    { "THR",  "CB", "OG1", "CG2", "CA" }
  };
#define NID asize(id)
  int i,j,k,l,aa[4],nimp;
  
  nimp = 0;
  for(i=0; (i<atoms->nres); i++) {
    /* Skip until residue number */
    for(j=0; (j<atoms->nr) && (atoms->atom[j].resnr != i); j++) 
      ;
    for(k=0; (k<NID); k++) {
      if ((id[k].res == NULL) || 
	  (strcmp(id[k].res,*atoms->resname[i]) == 0)) {
	/* This (i) is the right residue to look for this improper (k) */
	for(l=0; (l<4); l++)
	  aa[l] = find_atom(atoms,i,j,id[k].aa[l]);
	if ((aa[0] != -1) && (aa[1] != -1) && (aa[2] != -1) && (aa[3] != -1)) {
	  srenew(c->imp,nimp+1);
	  srenew(c->idih,nimp+1);
	  c->imp[nimp].ai = aa[0];
	  c->imp[nimp].aj = aa[1];
	  c->imp[nimp].ak = aa[2];
	  c->imp[nimp].al = aa[3];
	  nimp++;
	}
      }
    }
  }
  c->nimp = nimp;

  fprintf(log,"There are %d impropers\n",c->nimp);
  if (debug) {
    fprintf(debug,"Overview of improper dihedrals\n");
    for(i=0; (i<c->nimp); i++) { 
      fprintf(debug,"  %s",aname(atoms,c->imp[i].ai));
      fprintf(debug,"  %s",aname(atoms,c->imp[i].aj));
      fprintf(debug,"  %s",aname(atoms,c->imp[i].ak));
      fprintf(debug,"  %s\n",aname(atoms,c->imp[i].al));
    }
  }
}

void pr_dist(FILE *fp,bool bHeader,t_correct *c,int i)
{
  real ideal=0;
  
  if (bHeader)
    fprintf(fp,"#%4s%5s%10s%10s%10s\n","ai","aj","ideal","lb","ub");
  switch (c->d[i].cons_type) {
  case edcBOND:
    ideal = 5*(c->d[i].lb+c->d[i].ub);
    break;
  case edcNONBOND:
    ideal = 0;
    break;
  case edcDISRE:
    ideal = -1;
    break;
  default:
    fatal_error(0,"cons_type for distance %d = %d\n",i,c->d[i].cons_type);
  }
  fprintf(fp,"%5d%5d%10.5f%10.5f%10.5f\n",1+c->d[i].ai,1+c->d[i].aj,
	  ideal,10*c->d[i].lb,10*c->d[i].ub);
}

void pr_distances(FILE *fp,t_correct *c)
{
  int i;
  
  for(i=0; (i<c->ndist); i++) 
    pr_dist(fp,(i==0),c,i);
}

void add_dist(t_correct *c,int ai,int aj,real lb,real ideal,real ub,real w[])
{
  int n = c->ndist;
  
  if ((w[ai] != 0) || (w[aj] != 0)) {
    if (n == c->maxdist) {
      c->maxdist += 100;
      srenew(c->d,c->maxdist);
      srenew(c->ip,c->maxdist);
    }
    c->d[n].ai = ai;
    c->d[n].aj = aj;
    if (ideal > 0)
      c->d[n].cons_type = edcBOND;
    else if (ideal == 0.0)
      c->d[n].cons_type = edcNONBOND;
    else
      c->d[n].cons_type = edcDISRE;
    c->d[n].lb = lb;
    c->d[n].ub = ub;
    c->d[n].wi = w[ai]/(w[ai]+w[aj]);
    c->ip[n] = n;
    c->ndist++;
  }
}

void define_dist(FILE *log,t_topology *top,t_correct *c,real weight[])
{ 
  int  i,k,nra,ndist;
  int  ai,aj,tp,ftype;
  real lb=0,ub=0,b0;
  
  c->ndist = 0;
  c->d     = NULL;
  c->ip    = NULL;
  snew(c->bViol,top->atoms.nr);
  for(ftype=0; (ftype<F_NRE); ftype++) {
    if ((interaction_function[ftype].flags & IF_BOND) ||
	(ftype == F_DISRES)) {
      nra   = interaction_function[ftype].nratoms+1;
      
      ndist = top->idef.il[ftype].nr/nra;
  
      for(i=k=0; (i<ndist); i++,k+=nra) {
	tp = top->idef.il[ftype].iatoms[k];
	ai = top->idef.il[ftype].iatoms[k+1];
	aj = top->idef.il[ftype].iatoms[k+2];
	b0 = 0;
	switch (ftype) {
	case F_DISRES:
	  lb = top->idef.iparams[tp].disres.low;
	  ub = top->idef.iparams[tp].disres.up1;
	  break;
	case F_SHAKE:
	  b0 = top->idef.iparams[tp].shake.dA;
	  break;
	case F_BONDS:
	case F_G96BONDS:
	  b0 = top->idef.iparams[tp].harmonic.rA;
	  break;
	case F_MORSE:
	  b0 = top->idef.iparams[tp].morse.b0;
	  break;
	}
	if (b0 != 0.0) {
	  lb = b0*0.98;
	  ub = b0*1.02;
	}
	/* CHECK THIS! */
	b0 = (lb + ub)*0.5;
	add_dist(c,ai,aj,lb,b0,ub,weight);
      }
    }
  }
}

void read_dist(FILE *log,char *fn,int natom,t_correct *c,real weight[])
{
  FILE   *fp;
  char   buf[256];
  int    ai,aj,i,nline=0;
  double ideal,lb,ub;
  int    nnn[edcNR];
  
  for(i=0; (i<edcNR); i++)
    nnn[i] = 0;
    
  snew(c->bViol,natom);
  
  fp = ffopen(fn,"r");
  while (fgets(buf,255,fp)) {
    nline ++;
    if (buf[0] != '#') {
      if (sscanf(buf,"%d%d%lf%lf%lf",&ai,&aj,&ideal,&lb,&ub) != 5)
	fprintf(stderr,"Error in %s, line %d\n",fn,nline);
      else {
	add_dist(c,ai-1,aj-1,lb*0.1,ideal*0.1,ub*0.1,weight);
	nnn[c->d[c->ndist-1].cons_type]++;
      }
    }
  }
  ffclose(fp);
  
  fprintf(log,"Read distances from file %s\n",fn);
  for(i=0; (i<edcNR); i++)
    fprintf(log,"There were %d %s distances\n",
	    nnn[i],edc_names[i]);
}

real *read_weights(char *fn,int natom)
{
  t_atoms newatoms;
  int     i;
  rvec    *x;
  matrix  box;
  char    title[256];
  real    *w;

  /* Read the weights from the occupancy field in the pdb file */  
  snew(x,natom);
  init_t_atoms(&newatoms,natom,TRUE);
  read_pdb_conf(fn,title,&newatoms,x,box,FALSE);
  snew(w,newatoms.nr);
  for(i=0; (i<newatoms.nr); i++)
    w[i] = newatoms.pdbinfo[i].occup;
  free_t_atoms(&newatoms);
  sfree(x);
  
  return w;
}

void pr_corr(FILE *log,t_correct *c)
{
  fprintf(log,"Parameters for this disco run\n");
  fprintf(log,"maxnit     = %d\n",c->maxnit);
  fprintf(log,"nbcheck    = %d\n",c->nbcheck);
  fprintf(log,"nstprint   = %d\n",c->nstprint);
  fprintf(log,"nstranlist = %d\n",c->nstranlist);
  fprintf(log,"ngrow      = %d\n",c->ngrow);
  fprintf(log,"bExplicit  = %s\n",BOOL(c->bExplicit));
  fprintf(log,"bChiral    = %s\n",BOOL(c->bChiral));
  fprintf(log,"bPep       = %s\n",BOOL(c->bPep));
  fprintf(log,"bDump      = %s\n",BOOL(c->bDump));
  fprintf(log,"ndist      = %d\n",c->ndist);
  fprintf(log,"npep       = %d\n",c->npep);
  fprintf(log,"nimp       = %d\n",c->nimp);
  fprintf(log,"lowdev     = %g\n",c->lodev);
  fflush(log);
}

void measure_dist(FILE *log,int natom,rvec x[],t_correct *c,real weight[],
		  real cutoff)
{
  int    ai,aj;
  rvec   dx;
  real   ideal;
  
  snew(c->bViol,natom);

  for(ai=0; (ai<natom); ai++)
    for(aj=ai+1; (aj<natom); aj++) {
      rvec_sub(x[ai],x[aj],dx);
      ideal = 10.0*norm(dx);
      if ((ideal < cutoff)  || (cutoff == 0))
	add_dist(c,ai,aj,ideal*0.98,ideal,ideal*1.02,weight);
    }
}

void check_dist(FILE *log,t_correct *c)
{
  int  i;
  real tol=0.001;
  
  fprintf(log,"Checking distances for internal consistency\n");
  for(i=0; (i<c->ndist); i++) {
    if ((c->d[i].ub - c->d[i].lb) < tol) {
      pr_dist(log,TRUE,c,i);
    }
  }
}

void check_final(FILE *log,t_correct *c,rvec x[],bool bVerbose)
{
  int  i,ai,aj,nviol;
  rvec dx;
  real len;
  
  fprintf(log,"Structure check after optimisation\n");
  nviol = 0;
  for(i=0; (i<c->ndist); i++) {
    ai = c->d[i].ai;
    aj = c->d[i].aj;
    rvec_sub(x[ai],x[aj],dx);
    len = norm(dx);
    if ((len < c->d[i].lb) || (len > c->d[i].ub)) {
      if (bVerbose) {
	if ((nviol % 20) == 0)
	  fprintf(log,"%5s  %5s  %8s  %8s  %8s\n","ai","aj","len","lb","ub");
	fprintf(log,"%5d  %5d  %8.3f  %8.3f  %8.3f\n",ai,aj,len,
		c->d[i].lb,c->d[i].ub);
      }
      nviol++;
    }
  }
  fprintf(log,"Totally %d violations\n",nviol);
}
