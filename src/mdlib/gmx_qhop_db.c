#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "string2.h"
#include "futil.h"
#include "symtab.h"
#include "hackblock.h"
#include "resall.h"
#include "gmx_qhop_parm.h"
#include "gmx_qhop_xml.h"
#include "gmx_qhop_db.h"

#if 0
/* These routines link to invalid dependencies in kernel/resall.c! */

typedef struct {
  int id,charge,natom;
  int ndonor,*donor,nacceptor,*acceptor;
} qhop_resinfo_t;
  
typedef struct {
  int               nrtp;
  t_restp           *rtp;
  qhop_resinfo_t    *resinfo;
  int               bts[4];
  int               nrexcl;
  bool              bAllDih,bHH14,bRemoveDih;
  t_atomtype        atype;
  t_symtab          tab;
  int               ngqh;
  gmx_qhop          *gqh;
  t_qhop_parameters *qhop_param;
} gmx_qhop_db_t;

static void fill_resinfo(t_restp *rtp,qhop_resinfo_t *ri)
{
  int    j,k,m;
  bool   bDonor;
  double qtot;
  
  qtot = 0;
  ri->natom = rtp->natom;
  ri->ndonor = 0;
  snew(ri->donor,ri->natom);
  ri->nacceptor = 0;
  snew(ri->acceptor,ri->natom);
  for(j=0; (j<rtp->natom); j++) {
    qtot += rtp->atom[j].q;
    /* First check atoms usually involved in Hbonds */
    if ((rtp->atom[j].atomnumber == 7) ||
	(rtp->atom[j].atomnumber == 8) ||
	(rtp->atom[j].atomnumber == 9) ||
	(rtp->atom[j].atomnumber == 16)) {
      ri->acceptor[ri->nacceptor++] = j;
      /* Now check whether this atom has a proton, i.e. is a donor */
      bDonor = FALSE;
      for(k=0; (k<rtp->rb[ebtsBONDS].nb); k++) {
	if (strcasecmp(*(rtp->atomname[j]),rtp->rb[ebtsBONDS].b[k].a[0]) == 0)
	  for(m=0; (m<rtp->natom); m++) 
	    if ((strcasecmp(*(rtp->atomname[m]),rtp->rb[ebtsBONDS].b[k].a[1]) == 0) &&
		(rtp->atom[m].atomnumber == 1) &&
		(rtp->atom[m].q > 0))
	      bDonor = TRUE;
      }
      if (bDonor) 
	ri->donor[ri->ndonor++] = j;
    }
  }
  ri->charge = qtot;
}

static void lo_fill_qhp(gmx_qhop gqh,char *name,real *xx)
{
  double x;
  int    test;
  
  test = gmx_qhop_get_value(gqh,name,&x);
  if (test == 1) 
    *xx = x;
  else {
    fprintf(stderr,"WARNING: Could not extract %s from qhop data for %s-%s\n",
	    name,gmx_qhop_get_donor(gqh),gmx_qhop_get_acceptor(gqh));
    *xx = 0;
  }
}

static void fill_qhp(t_qhop_parameters *qhp,gmx_qhop gqh)
{
  lo_fill_qhp(gqh,"alpha",&qhp->alpha);
  lo_fill_qhp(gqh,"beta",&qhp->beta);
  lo_fill_qhp(gqh,"gamma",&qhp->gamma);
  lo_fill_qhp(gqh,"k_1",&qhp->k_1);
  lo_fill_qhp(gqh,"k_2",&qhp->k_2);
  lo_fill_qhp(gqh,"k_3",&qhp->k_3);
  lo_fill_qhp(gqh,"m_1",&qhp->m_1);
  lo_fill_qhp(gqh,"m_2",&qhp->m_2);
  lo_fill_qhp(gqh,"m_3",&qhp->m_3);
  lo_fill_qhp(gqh,"s_A",&qhp->s_A);
  lo_fill_qhp(gqh,"t_A",&qhp->t_A);
  lo_fill_qhp(gqh,"v_A",&qhp->v_A);
  lo_fill_qhp(gqh,"s_B",&qhp->s_B);
  lo_fill_qhp(gqh,"s_C",&qhp->s_C);
  lo_fill_qhp(gqh,"t_C",&qhp->t_C);
  lo_fill_qhp(gqh,"v_C",&qhp->v_C);
  lo_fill_qhp(gqh,"f",&qhp->f);
  lo_fill_qhp(gqh,"g",&qhp->g);
  lo_fill_qhp(gqh,"h",&qhp->h);
  lo_fill_qhp(gqh,"p_1",&qhp->p_1);
  lo_fill_qhp(gqh,"q_1",&qhp->q_1);
  lo_fill_qhp(gqh,"q_2",&qhp->q_2);
  lo_fill_qhp(gqh,"q_3",&qhp->q_3);
  lo_fill_qhp(gqh,"r_1",&qhp->r_1);
  lo_fill_qhp(gqh,"r_2",&qhp->r_2);
  lo_fill_qhp(gqh,"r_3",&qhp->r_3);
}

static void dump_qhp(FILE *fp,t_qhop_parameters *qhp)
{
  fprintf(fp,"alpha = %g\n",qhp->alpha);
  fprintf(fp,"beta  = %g\n",qhp->beta);
  fprintf(fp,"gamma = %g\n",qhp->gamma);
  fprintf(fp,"k_1   = %g\n",qhp->k_1);
  fprintf(fp,"k_2   = %g\n",qhp->k_2);
  fprintf(fp,"k_3   = %g\n",qhp->k_3);
  fprintf(fp,"m_1   = %g\n",qhp->m_1);
  fprintf(fp,"m_2   = %g\n",qhp->m_2);
  fprintf(fp,"m_3   = %g\n",qhp->m_3);
  fprintf(fp,"s_A   = %g\n",qhp->s_A);
  fprintf(fp,"t_A   = %g\n",qhp->t_A);
  fprintf(fp,"v_A   = %g\n",qhp->v_A);
  fprintf(fp,"s_B   = %g\n",qhp->s_B);
  fprintf(fp,"s_C   = %g\n",qhp->s_C);
  fprintf(fp,"t_C   = %g\n",qhp->t_C);
  fprintf(fp,"v_C   = %g\n",qhp->v_C);
  fprintf(fp,"f     = %g\n",qhp->f);
  fprintf(fp,"g     = %g\n",qhp->g);
  fprintf(fp,"h     = %g\n",qhp->h);
  fprintf(fp,"p_1   = %g\n",qhp->p_1);
  fprintf(fp,"q_1   = %g\n",qhp->q_1);
  fprintf(fp,"q_2   = %g\n",qhp->q_2);
  fprintf(fp,"q_3   = %g\n",qhp->q_3);
  fprintf(fp,"r_1   = %g\n",qhp->r_1);
  fprintf(fp,"r_2   = %g\n",qhp->r_2);
  fprintf(fp,"r_3   = %g\n",qhp->r_3);
}

gmx_qhop_db gmx_qhop_db_read(char *forcefield)
{
  gmx_qhop_db_t *qdb;
  char buf[256];
  char *fn;
  int i,j;
  double qtot;
  
  snew(qdb,1);
  open_symtab(&(qdb->tab));
  sprintf(buf,"%s-qhop",forcefield);
  qdb->atype = read_atype(forcefield,&(qdb->tab));
  qdb->nrtp = read_resall(buf,qdb->bts,&(qdb->rtp),qdb->atype,
			  &(qdb->tab),&(qdb->bAllDih),
			  &(qdb->nrexcl),&(qdb->bHH14),&(qdb->bRemoveDih));
  snew(qdb->resinfo,qdb->nrtp);
  for(i=0; (i<qdb->nrtp); i++) 
    fill_resinfo(&qdb->rtp[i],&qdb->resinfo[i]);
    
  sprintf(buf,"%s-qhop.dat",forcefield);
  fn = (char *)libfn(buf);
  /* Read the xml data file */
  qdb->gqh = gmx_qhops_read(fn,&qdb->ngqh);
  sprintf(buf,"%s-qhop-debug.dat",forcefield);
  gmx_qhops_write(buf,qdb->ngqh,qdb->gqh);
  snew(qdb->qhop_param,qdb->ngqh);
  for(i=0; (i<qdb->ngqh); i++) {
    fill_qhp(&(qdb->qhop_param[i]),qdb->gqh[i]);
    if (debug) {
      fprintf(debug,"Donor: %s  Acceptor: %s\n",
	      gmx_qhop_get_donor(qdb->gqh[i]),
	      gmx_qhop_get_acceptor(qdb->gqh[i]));
      dump_qhp(debug,&(qdb->qhop_param[i]));
    }
  }
  sfree(fn);
  return (gmx_qhop_db) qdb;
}

int gmx_qhop_db_write(char *fn,gmx_qhop_db qdb)
{
  gmx_qhop_db_t *db = (gmx_qhop_db_t *) qdb;
  FILE *fp;
  
  fp=ffopen(fn,"w");
  print_resall(fp,db->bts,db->nrtp,db->rtp,db->atype,db->bAllDih,
	       db->nrexcl,db->bHH14,db->bRemoveDih);
  fclose(fp);
  
  return 1;
}

int gmx_qhop_db_done(gmx_qhop_db qdb)
{
  fprintf(stderr,"gmx_qhop_db_done not implemented yet.\n");

  return 1;
}

int gmx_qhop_db_get_nstates(gmx_qhop_db qdb,char *resname)
{
  gmx_qhop_db_t *db = (gmx_qhop_db_t *) qdb;
  int i,nstate=0;
  
  for(i=0; (i<db->nrtp); i++) {
    if (strncmp(db->rtp[i].resname,resname,3) == 0) {
      nstate++;
    }
  }
  return nstate;
}

int gmx_qhop_db_get_qstate(gmx_qhop_db qdb,char *resname,int istate)
{
  gmx_qhop_db_t *db = (gmx_qhop_db_t *) qdb;
  int i,nstate=0;
  
  for(i=0; (i<db->nrtp); i++) {
    if (strncmp(db->rtp[i].resname,resname,3) == 0) {
      nstate++;
      if (nstate == istate)
	return db->resinfo[i].charge;
    }
  }
  return nstate;
}

char **gmx_qhop_db_get_donors(gmx_qhop_db qdb,char *resname,int state)
{
  return NULL;
}
 
char **gmx_qhop_db_get_acceptors(gmx_qhop_db qdb,char *resname,int state)
{
  return NULL;
}

int gmx_qhop_db_set_charges(gmx_qhop_db qdb,char *resname,int state,
			    int natoms,real q[])
{
  return 0;
}

int gmx_qhop_db_get_parameters(gmx_qhop_db qdb,
			       char *donor,char *acceptor,
			       t_qhop_parameters *qp)
{
  gmx_qhop_db_t *db = (gmx_qhop_db_t *) qdb;
  char *aa,*dd;
  int i;
  
  for(i=0; (i<db->ngqh); i++) {
    aa = gmx_qhop_get_acceptor(db->gqh[i]);
    dd = gmx_qhop_get_donor(db->gqh[i]);
    if (strncmp(donor,dd,3) == 0) {
      if (strncmp(acceptor,aa,3) == 0) {
	memcpy(qp,&(db->qhop_param[i]),sizeof(*qp));
	
	return 1;
      }
    }
  }
  return 0;
}

#else
int gmx_qhop_db_dummy=0;
#endif
