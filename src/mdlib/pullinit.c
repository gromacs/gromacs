#include <ctype.h>
#include "princ.h"
#include <stdlib.h>
#include "readinp.h"
#include "sysstuff.h"
#include "futil.h"
#include "statutil.h"
#include "vec.h"
#include "smalloc.h"
#include "typedefs.h"
#include "names.h"
#include "fatal.h"
#include "macros.h"
#include "rdgroup.h"
#include "symtab.h"
#include "index.h"
#include "confio.h"
#include "pull.h"

static void get_pullmemory(t_pullgrps *pull, int ngrps) 
{
  snew(pull->ngx,ngrps);   
  snew(pull->x_con,ngrps);     
  snew(pull->xprev,ngrps);
  snew(pull->tmass,ngrps); 
  snew(pull->idx,ngrps);      
  snew(pull->f,ngrps);        
  snew(pull->spring,ngrps);        
  snew(pull->dir,ngrps);
  snew(pull->x_unc,ngrps);   
  snew(pull->x_ref,ngrps);
  snew(pull->d_ref,ngrps);
  snew(pull->x0,ngrps);
  snew(pull->xp,ngrps);
  snew(pull->comhist,ngrps);
}

static void print_info(FILE *log,t_pull *pull) 
{

  fprintf(log,"\n**************************************************\n"
	  "                         PULL INFO                      \n"
	  "**************************************************\n");

  switch (pull->runtype) {
  case eStart:
    fprintf(log,"RUN TYPE: Generating Starting structures\n");
    break;
  case eAfm:
    fprintf(log,"RUN TYPE: Afm\n");
    break;
  case eConstraint:
    fprintf(log,"RUN TYPE: Constraint\n");
    break;
  case eUmbrella:
    fprintf(log,"RUN TYPE: Umbrella sampling\n");
    break;
  case eTest:
    fprintf(log,"RUN TYPE: Test run\n");
    break;
  default:
    fprintf(log,"RUN TYPE: WARNING! pullinit does not know this runtype\n");
  }
  
  if (pull->runtype == eConstraint || pull->runtype == eStart) 
    switch (pull->reftype) {
    case eCom: 
      fprintf(log,"REFERENCE TYPE: center of mass of reference group\n");
      break;
    case eComT0:
      fprintf(log,
	      "REFERENCE TYPE: center of mass of reference group at t=0\n");
      break;
    case eDyn:
      fprintf(log,
	      "REFERENCE TYPE: center of mass of dynamically made groups\n");
      fprintf(log,"Using dynamic reference groups: r=%8.3f, rc=%8.3f\n",
	      pull->r,pull->rc);
      break;
    case eDynT0:
      fprintf(log,
	      "REFERENCE TYPE: center of mass of dynamically made groups,\n"
	      "based on the positions of its atoms at t=0\n");
      fprintf(log,"Using dynamic reference groups: r=%8.3f, rc=%8.3f\n",
	      pull->r,pull->rc);
      break;
    default:
      fprintf(log,"REFERENCE TYPE: no clue! What did you do wrong?\n");
    }
}

static void get_named_indexgroup(FILE *log,atom_id **target,int *isize,
				 char *name,atom_id **index,int ngx[],
				 char **grpnames,int ngrps) 
{
  int i,j;
  bool bFound = FALSE;
  atom_id *tmp;

  fprintf(log,"Looking for group %s: ",name);
  for (i=0;i<ngrps;i++) {
    if (!strcmp(name,grpnames[i])) {
      /* found the group we're looking for */
      snew(tmp,ngx[i]);
      for (j=0;j<ngx[i];j++) 
	tmp[j]=index[i][j];
      *isize=ngx[i];
      bFound = TRUE;
      fprintf(log,"found group %s: %d elements. First: %d\n",name,ngx[i],
	      tmp[0]);
    }
  }
  
  *target=tmp;
  if (!bFound) 
    fatal_error(0,"Can't find group %s in the index file",name);
}

static void read_whole_index(char *indexfile,char ***grpnames,atom_id ***index,
			      int **ngx,int *totalgrps)
{
  t_block *grps;
  char    **gnames;
  int i,j;

  if (!indexfile)
    fatal_error(0,"No index file specified");

  grps = init_index(indexfile,&gnames);
  if (grps->nr==0) 
    fatal_error(0,"No groups in indexfile\n");
  
  *totalgrps = grps->nr;
  snew(*index,grps->nr);
  snew(*grpnames,grps->nr);
  snew(*ngx,grps->nr);

  for(i=0; (i<*totalgrps); i++) {
    (*grpnames)[i]=strdup(gnames[i]);
    (*ngx)[i]=grps->index[i+1]-grps->index[i];
    snew((*index)[i],(*ngx)[i]);
    for(j=0; (j<(*ngx)[i]); j++)
      (*index)[i][j]=grps->a[grps->index[i]+j];
  }
}

static void print_whole_index(char **grpnames, atom_id **index, int *ngx, int
			      ngrps) 
{
  int i,j;
  FILE *tmp;
  
  tmp = ffopen("indexcheck","w");
  for (i=0;i<ngrps;i++) {
    fprintf(tmp,"\nGrp %d: %s. %d elements\n",i,grpnames[i],ngx[i]);
    for (j=0;j<ngx[i];j++)
      fprintf(tmp," %d ",index[i][j]);
  }
  fflush(tmp);
}

void init_pull(FILE *log,int nfile,t_filenm fnm[],t_pull *pull,rvec *x,
	       t_mdatoms *md,rvec boxsize) 
{
  int i,j,m,ii;
  int ngrps;
  real tm;
  rvec tmp;
  matrix box;
  char **grpnames;
  atom_id **index;
  int *ngx;
  int totalgrps;    /* total number of groups in the index file */

  clear_mat(box);
  for (m=0;m<DIM;m++) 
    box[m][m]=boxsize[m];
  
  /* do we have to do any pulling at all? If not return */
  pull->bPull = opt2bSet("-pi",nfile,fnm);
  if (!pull->bPull) return;

  pull->out = ffopen(opt2fn("-po",nfile,fnm),"w");
  read_pullparams(pull, opt2fn("-pi",nfile,fnm));
  ngrps = pull->pull.n;

  if (pull->reftype == eDyn || pull->reftype == eDynT0) 
    pull->bCyl = TRUE;
  else
    pull->bCyl = FALSE;
  
  if (pull->bCyl && (pull->rc < 0.01 || pull->r < 0.01)) 
    fatal_error(0,"rc or r is too small or zero.");

  print_info(log,pull);

  get_pullmemory(&pull->pull,ngrps);
  get_pullmemory(&pull->dyna,ngrps);
  get_pullmemory(&pull->ref,1);
  
  /* read the whole index file */
  read_whole_index(opt2fn("-n",nfile,fnm),&grpnames,&index,&ngx,&totalgrps);
  
  /* DEBUG */
  if (pull->bVerbose) {
    fprintf(stderr,"read_whole_index: %d groups total\n",totalgrps);
    for(i=0;i<totalgrps;i++) 
      fprintf(stderr,"group %i (%s) %d elements\n",
	      i,grpnames[i],ngx[i]);
    /*    print_whole_index(grpnames,index,ngx,totalgrps); */
  } /* END DEBUG */
  
  /* grab the groups that are specified in the param file */
  for (i=0;i<pull->pull.n;i++) 
    get_named_indexgroup(log,&pull->pull.idx[i],&pull->pull.ngx[i],
			 pull->pull.grps[i],index,ngx,grpnames,totalgrps) ;
  get_named_indexgroup(log,&pull->ref.idx[0],&pull->ref.ngx[0],
		       pull->ref.grps[0],index,ngx,grpnames,totalgrps);
  
  /* get more memory! Don't we love C? */
  snew(pull->ref.x0[0],pull->ref.ngx[0]);
  snew(pull->ref.xp[0],pull->ref.ngx[0]);
  snew(pull->ref.comhist[0],pull->reflag);
  for (i=0;i<pull->pull.n;i++) 
    snew(pull->dyna.comhist[i],pull->reflag);
  
  for (i=0;i<ngrps;i++) {
    tm = calc_com(x,pull->pull.ngx[i],pull->pull.idx[i],
 		  md,tmp,box);
    copy_rvec(tmp,pull->pull.x_con[i]);
    copy_rvec(tmp,pull->pull.x_unc[i]);
    copy_rvec(tmp,pull->pull.x_ref[i]);
    copy_rvec(tmp,pull->pull.spring[i]);
    fprintf(log,"Initializing pull groups. Mass of group %d: %8.3f\n"
	    "Initial coordinates center of mass: %8.3f %8.3f %8.3f\n",
	    i,tm,tmp[XX],tmp[YY],tmp[ZZ]);
    pull->pull.tmass[i] = tm;
  }
  
  /* initialize the reference group, in all cases */
  tm = calc_com(x,pull->ref.ngx[0],pull->ref.idx[0],md,
		tmp,box);
  
  copy_rvec(tmp,pull->ref.x_unc[0]);
  copy_rvec(tmp,pull->ref.x_con[0]);
  copy_rvec(tmp,pull->ref.x_ref[0]);
  
  for (j=0;j<pull->reflag;j++) 
    copy_rvec(pull->ref.x_unc[0],pull->ref.comhist[0][j]);
  
  fprintf(log,"Initializing reference group. Mass: %8.3f\n"
	  "Initial coordinates center of mass: %8.3f %8.3f %8.3f\n",
	  tm,tmp[XX],tmp[YY],tmp[ZZ]);
  /* keep the initial coordinates for center of mass at t0 */
  for (j=0;j<pull->ref.ngx[0];j++) {
    copy_rvec(x[pull->ref.idx[0][j]],pull->ref.x0[0][j]);
    copy_rvec(x[pull->ref.idx[0][j]],pull->ref.xp[0][j]);
  }
  pull->ref.tmass[0] = tm;

  /* if we use dynamic reference groups, do some initialising for them */
  if (pull->bCyl) {
    make_refgrps(pull,x,md);
    for (i=0;i<ngrps;i++) {
      copy_rvec(pull->dyna.x_unc[i],pull->dyna.x_con[i]);
      copy_rvec(pull->dyna.x_unc[i],pull->dyna.x_ref[i]);
      
      /* initialize comhist values for running averages */
      for (j=0;j<pull->reflag;j++) 
	copy_rvec(pull->dyna.x_unc[i],pull->dyna.comhist[i][j]);
      
      fprintf(log,"Initializing dynamic groups. %d atoms. Weighted mass"
	      "for %d:%8.3f\n"
	      "Initial coordinates center of mass: %8.3f %8.3f %8.3f\n",
	      pull->dyna.ngx[i],i,pull->dyna.tmass[i],pull->dyna.x_ref[i][XX],
	      pull->dyna.x_ref[i][YY],pull->dyna.x_ref[i][ZZ]);
    }
  } 
  
  /* set the reference distances and directions, taking into account pbc */
  for (i=0;i<ngrps;i++) {
    if (pull->bCyl) {
      rvec_sub(pull->pull.x_ref[i],pull->dyna.x_ref[i],tmp);
      for (m=0;m<DIM;m++) {
	if (tmp[m] < -0.5*boxsize[m]) tmp[m] += boxsize[m];
	if (tmp[m] >  0.5*boxsize[m]) tmp[m] -= boxsize[m];
      }
    } else { 
      rvec_sub(pull->pull.x_ref[i],pull->ref.x_ref[0],tmp);
      for (m=0;m<DIM;m++) {
	if (tmp[m] < -0.5*boxsize[m]) tmp[m] += boxsize[m];
	if (tmp[m] >  0.5*boxsize[m]) tmp[m] -= boxsize[m];
      }
    }
    
    /* reference distance for constraint run */
    pull->pull.d_ref[i] = norm(tmp);
    
    /* select elements of direction vector to use for Afm and Start runs */
    for (m=0;m<DIM;m++)  
      tmp[m] = tmp[m] * pull->dims[m];
    svmul(1/norm(tmp),tmp,pull->pull.dir[i]);
    if (pull->bReverse) 
      svmul(-1.0,pull->pull.dir[i],pull->pull.dir[i]);

    if (pull->runtype == eAfm || pull->runtype == eStart)  
      fprintf(log,"\nPull direction: %8.3f %8.3f %8.3f\n",
	      pull->pull.dir[i][XX],pull->pull.dir[i][YY],
	      pull->pull.dir[i][ZZ]);
    
  }
  fprintf(log,"**************************************************\n"
	  "                      END   PULL INFO                    \n"
	  "**************************************************\n\n");
}


