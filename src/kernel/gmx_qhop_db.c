#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "typedefs.h"
#include "smalloc.h"
#include "string2.h"
#include "gmx_fatal.h"
#include "futil.h"
#include "symtab.h"
#include "hackblock.h"
#include "resall.h"
#include "mdatoms.h"
#include "types/gmx_qhop_types.h"
#include "gmx_qhop_parm.h"
#include "gmx_qhop_xml.h"
#include "gmx_qhop_db.h"

#define _ERIKS_SPARKLING_NEW_CODE
#define BUFSIZE 1000
/* First of all, a safe strdup, more or less copied from hackblock.c. */
char *safe_strdup(const char *s)
{
  if (s == NULL)
    return NULL;
  else
    return strdup(s);
}

/* This one was copied from trim_string in symtab.c and modified. */
char *trim_strndup(const char *s, int maxlen)
     /*
      * Returns a pointer to a newly allocated area which contains a copy 
      * of s without leading or trailing spaces. Strings are
      * truncated to BUFSIZE positions.
      */
{
  int len,i;
  if(strlen(s)>(size_t)(maxlen-1))
    gmx_fatal(FARGS,"Character buffer size too small\n");
  
  for (; (*s)&&((*s)==' '); s++);
  for (len=strlen(s); (len>0); len--) if (s[len-1]!=' ') break;
  if (len>=BUFSIZE) len=BUFSIZE-1;

  if (s == NULL)
    return NULL;
  else
    return strndup(s, len);
}


/* It seems we need afew special copy routines here, since
 * copy_t_restp() doesn't seem to work for our purpouses. */

static void qhop_copy_t_rbondeds(t_rbondeds *s, t_rbondeds *d, t_symtab *tab)
{
  int i,j;
  d->nb = s->nb;
  if (d->b == NULL)
    snew(d->b, d->nb);
  
  for (i=0; i<s->nb; i++)
    {
      for(j=0; j<MAXATOMLIST; j++)
	(d->b[i]).a[j] = *(put_symtab(tab,(s->b[i]).a[j]));
      (d->b[i]).s = *(put_symtab(tab, (s->b[i]).s));
    }
}

static void qhop_copy_t_restp(t_restp *src, t_restp *dest, t_symtab *tab)
{
  int i, j;
  dest->resname = *(put_symtab(tab,src->resname));
  dest->natom = src->natom;
  snew(dest->atom, dest->natom);
  snew(dest->atomname, dest->natom);
    
  for (i=0; i<dest->natom; i++)
    {
      dest->atom[i].m           = src->atom[i].m;
      dest->atom[i].q           = src->atom[i].q;
      dest->atom[i].type        = src->atom[i].type;
      dest->atom[i].typeB       = src->atom[i].typeB;
      dest->atom[i].ptype       = src->atom[i].ptype;
      dest->atom[i].resind      = src->atom[i].resind;
      dest->atom[i].atomnumber  = src->atom[i].atomnumber;
      strncpy(dest->atom[i].elem, src->atom[i].elem,4);

      snew(dest->atomname[i],1);
      dest->atomname[i] = put_symtab(tab, *(src->atomname[i]));
    }
  snew(dest->cgnr, dest->natom);
  for (i=0; i<dest->natom; i++)
    dest->cgnr[i] = src->cgnr[i];
  qhop_copy_t_rbondeds(src->rb, dest->rb, tab);
}

/* Checks to see if name refers to an entry in the symtab,
 * not checking the content but the pointer itself.
 * Then moves the entry to the symtab, frees the old string data
 * and has name point to the entry in the symtab.
 * NOT safe if name points to another symtab, cuz that symtab may be messed up.
 * USE WITH CARE! */
/* ****************<  Redundant now. >*************** */
/* static void safe_move_to_symtab(t_symtab *tab, char *name) */
/* { */
/*   t_symbuf *sb=tab->symbuf; */
/*   bool hit = FALSE; */
/*   char **newname; */
/*   if (name[0] != '\0') */
/*     { */
/*       fprintf(stderr, "   '%s' encountered\n", name); */
/*       while (!hit && sb != NULL) */
/* 	{ */
/* 	  hit = (*sb->buf == name); */
/* 	  sb = sb->next; */
/* 	} */
/*       /\* Hmm. What if it points to ANOTHER symtab? Make sure that is not the case first. *\/ */
/*       if (!hit) */
/* 	{ */
/* 	  newname = put_symtab(tab, name); */
/* 	  sfree(name); */
/* 	  name = *newname; */
/* 	} */
/*     } */
/* } */

/* Moves strings from atomnames and resnames and such to
 * the symtab if they're not there already.
 * Will NOT work if some names point to the symtab from the beginning. */
/* *******************< This is redundant nowadays. >****************** */
/* static void move_strings_to_symtab(qhop_db_t qdb) */
/* { */
/*   t_restp *r; */
/*   t_symtab *tab = &(qdb->tab); */
/*   t_symbuf *sb; */
/*   qhop_resblocks    *rb = &(qdb->rb); */
/*   qhop_reactant     *qr; */
/*   int i,j,k,l; */
/*   bool points_to_tab; */
/*   pr_symtab(stderr, 2, "Before the move\n", tab); */
/*   fprintf(stderr, " -- RESIDUES --\n"); */
/*   /\* Residues *\/ */
/*   for (i=0; i<qdb->nrtp; i++) */
/*     { */
/*       r = &(qdb->rtp[i]); */
/*       safe_move_to_symtab(tab, r->resname); */
/*       for (j=0; j<r->natom; j++) */
/* 	safe_move_to_symtab(tab, *(r->atomname[j])); */
/*     } */

/*   fprintf(stderr, " -- INTERACTIONS --\n"); */
/*   /\* Interactions *\/ */

/*   for (i=0; i<rb->nf; i++) */
/*     safe_move_to_symtab(tab, rb->files[i]); */
  
/*   /\* Resblocks *\/ */

/*   for (i=0; i<rb->nrestypes; i++) */
/*     { */
/*       safe_move_to_symtab(tab, rb->restype[i]); */
/*       for (j=0; j<rb->nres[i]; j++) */
/* 	{ */
/* 	  /\* Qhop_res *\/ */
/* 	  safe_move_to_symtab(tab, rb->res[i][j].name); */

/* 	  /\* Acceptors *\/ */
/* 	  for (k=0; k<rb->res[i][j].na; k++) */
/* 	    { */
/* 	      qr = &(rb->res[i][j].acc[k]); */

/* 	      /\*safe_move_to_symtab(tab, qr[k].name);*\/ */
/* 	      safe_move_to_symtab(tab, qr[k].product); */
/* 	      for (l=0; l<qr[k].nH; l++) */
/* 		{ */
/* 		  /\* Acceptor hydrogens *\/ */
/* 		  safe_move_to_symtab(tab, qr[k].H[l]); */
/* 		} */
/* 	    } */

/* 	  /\* Donors *\/ */
/* 	  for (k=0; k<rb->res[i][j].nd; k++) */
/* 	    { */
/* 	      qr = &(rb->res[i][j].don[k]); */
	  
/* 	      /\*safe_move_to_symtab(tab, qr[k].name);*\/ */
/* 	      safe_move_to_symtab(tab, qr[k].product); */
/* 	      for (l=0; l<qr[k].nH; l++) */
/* 		{ */
/* 		  /\* Donor hydrogens *\/ */
/* 		  safe_move_to_symtab(tab, qr[k].H[l]); */
/* 		} */
/* 	    } */
	    
/* 	} */
/*     } */
/*   pr_symtab(stderr, 2, "After the move\n", tab); */
/* } */

void set_reactant_products(qhop_db *qdb)
{
  int rb, r, p, reac, nreac[2], i;
  qhop_reactant *da;
  for (rb=0; rb<qdb->rb.nrestypes; rb++)
    {
      for (r=0; r<qdb->rb.nres[rb]; r++)
	{
	  nreac[0] = qdb->rb.res[rb][r].na;
	  nreac[1] = qdb->rb.res[rb][r].nd;

	  for (i=0; i<2; i++) /* Acceptors and donors */
	    {
	      da = (i==0) ? qdb->rb.res[rb][r].acc : qdb->rb.res[rb][r].don;
	      for (reac=0; reac<nreac[i]; reac++)
		for (p=0; p<qdb->rb.nres[rb]; p++)
		  if (p!=r)
		    if (strcmp(da[reac].product, qdb->rb.res[rb][p].name) == 0)
		      { /* This is the product. Set the pointer and proceed
			 * with the other acceptors/donors */
			da[reac].productdata = &(qdb->rb.res[rb][p]);
			break;
		      }
	    }
	}
    }
}
  


/* Takes a rtp database and picks out the residues found in qdb->rb->res[][].
 * The resulting slice is stored in qdb->rtp, and its size in qdb->nrtp.
 * Move strings to symtab in qdb. */
static void strip_rtp(char *ff, qhop_db *qdb, t_restp *bigrtp, int nbigrtp)
{
  printf("trim the rtp\n");
  int i,rt,r,a;
  bool match;
  qdb->nrtp = 0;
  for (i=0; i<nbigrtp; i++)
    {
      fprintf(stderr, "rtp entry no %d", i);
      match = FALSE;
      for (rt=0; rt<qdb->rb.nrestypes && !match; rt++) /* resblock */
	{
	  for (r=0; r<qdb->rb.nres[rt] && !match; r++) /* res */
	    {
	      if (match = ((strcmp(bigrtp[i].resname, qdb->rb.res[rt][r].name) == 0)))
		{
		  fprintf(stderr, " FOUND!");
		  /* Keep this entry */
		  if (qdb->nrtp == 0)
		    snew(qdb->rtp,1);
		  else
		    srenew(qdb->rtp,(qdb->nrtp)+1);
		  memset((void*)(&(qdb->rtp[qdb->nrtp])), 0, sizeof(t_restp));
		  /* Copy content, don't trust copy t_restp. */
		  /* copy_t_restp(&bigrtp[i], &(qdb->rtp[qdb->nrtp++])); */
		  
		  qhop_copy_t_restp(&(bigrtp[i]), &(qdb->rtp[(qdb->nrtp)++]), &(qdb->tab));
		}
	    }
	}
      fprintf(stderr, "\n");
    }
}

static int init_qhop_H_exist(gmx_mtop_t *top, t_mdatoms *mda, qhop_H_exist *Hext)
{
  
  int a, nH = 0, n;
  char *H;
  atom_id *H2atomid;
  int *atomid2H;
  t_atom *atom;
  n = mda->nr;
  snew(atomid2H, n);
  /* Find the hydrogens */  
  for (a=0; a<n; a++)
    {
      atomid2H[a] = -1;
      gmx_mtop_atomnr_to_atom(top,a,&atom);
      
      if (atom->atomnumber == 1 || (atom->atomnumber == NOTSET && atom->elem[0] == 'H' ))
	{ /* It's a H allright. */
	  srenew(H, nH+1);
	  srenew(H2atomid, nH+1);
	  H[nH] = (char)0;
	  H2atomid[nH] = a;
	  atomid2H[a] = nH;
	  nH++;
	}
    }
  Hext->H = H;
  Hext->H2atomid = H2atomid;
  Hext->atomid2H = atomid2H;
  return nH;
}

static void clear_qhop_H_exist(qhop_H_exist Hext)
{
  sfree(Hext.H);
  sfree(Hext.H2atomid);
  sfree(Hext.atomid2H);
}

/* These routines link to invalid dependencies in kernel/resall.c! */
/* Not anymore they don't! */
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

static void lo_fill_qhp(qhop_t gqh,char *name,real *xx)
{
  double x;
  int    test;
  
  test = qhop_get_value(gqh,name,&x);
  if (test == 1) 
    *xx = x;
  else {
    fprintf(stderr,"WARNING: Could not extract %s from qhop data for %s-%s\n",
	    name,qhop_get_donor(gqh),qhop_get_acceptor(gqh));
    *xx = 0;
  }
}

static void fill_qhp(qhop_parameters *qhp,qhop_t gqh)
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

static void dump_qhp(FILE *fp,qhop_parameters *qhp)
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

void qhop_db_print (qhop_parameters *qhp)
{
  dump_qhp(stderr, qhp);
}

/* static void fill_resblocks(qhop_parameters *qhp, qhop_resblocks_t rb) */
/* { */
/*   int i; */
/*   i = 0; */
/* } */

qhop_db_t qhop_db_read(char *forcefield, gmx_mtop_t *top, t_mdatoms *mda)
{
  qhop_db_t qdb;
  char buf[256];
  char *fn;
  int i,j,nrtp=0, nah;
  double qtot;
  t_atomtype atype;
  t_restp *bigrtp=NULL;
  t_symtab *stab=NULL;
  t_hackblock *ah=NULL;

  snew(stab,1);
  open_symtab(stab);

  snew(qdb,1);
  open_symtab(&(qdb->tab));
  /*  sprintf(buf,"%s-qhop",
      forcefield);*/
  atype = read_atype(forcefield, TRUE ,stab);
  /* read_resall has changed a lot! Rework this part. */
  /* nrtp = read_resall(forcefield,qdb->bts,&(bigrtp),atype, */
/* 		     stab,&(qdb->bAllDih), */
/* 		     &(qdb->nrexcl),&(qdb->bHH14),&(qdb->bRemoveDih)); */
  read_resall(forcefield, &nrtp, &bigrtp, atype, stab, FALSE);

  nah=read_h_db(forcefield,&ah);

  /* Read termini database too.
   * Not implemented yet, which is why amber is the
   * first forcefield to play with. */

  

  /* snew(qdb->resinfo,qdb->nrtp); */

  /* I'm redoing this part from scratch in a new manner.
   * No longer reading a special file; rather, the usual
   * rtp-file will be read along with extra itp-files if needed,
   * to get the interaction parameters for the slice of residues
   * listed in ffXXX-qhop.dat and the ones in the tpr-file.
   * / Erik Marklund */


  /* Process rtp-info AFTER reading the ffXXX-qhop.dat */
    
  sprintf(buf,"%s-qhop.dat",forcefield);
  fn = (char *)gmxlibfn(buf);
  /* Read the xml data file */
  qhops_read(fn, qdb);
  sprintf(buf,"%s-qhop-debug.dat",forcefield);
  qhops_write(buf,qdb->ngqh,qdb->gqh);
  snew(qdb->qhop_param,qdb->ngqh);
  for(i=0; (i<qdb->ngqh); i++) {
    fill_qhp(&(qdb->qhop_param[i]),qdb->gqh[i]);
    if (debug) {
      fprintf(debug,"Donor: %s  Acceptor: %s\n",
	      qhop_get_donor(qdb->gqh[i]),
	      qhop_get_acceptor(qdb->gqh[i]));
      dump_qhp(debug,&(qdb->qhop_param[i]));
    }
  }

  /* Must process the resblocks stuff */
  strip_rtp(forcefield, qdb, bigrtp, nrtp); /* Take away stuff not found in ffXXX-qhop.dat */
  done_symtab(stab);

  set_reactant_products(qdb);
  /* qhop_build_interaction_lib(forcefield, qdb, &atype, &stab); */

  /* This should be coded into the xml read instead. */
  //  move_strings_to_symtab(qdb); /* Copy names to reduced symtab in qdb->tab*/
  printf("free the redundant rtp parts.\n");
  for(i=0; i<nrtp; i++) /* Can't use free_t_restp. Not entirely sure why, though.*/
    {
      sfree(bigrtp[i].resname);
      sfree(bigrtp[i].atom);
      /* Will NOT free the actual atomnames, as thy may be used elsewhere. */  
      sfree(bigrtp[i].atomname);
      sfree(bigrtp[i].cgnr);
    }

  /* init the hydrogen existence array. */
  qdb->H_map.nH = init_qhop_H_exist(top, mda, &(qdb->H_map));
  
  done_atomtype(atype);
  sfree(fn);
  return qdb;
}

int qhop_db_write(char *fn,qhop_db *qdb)
{
  /*  qhop_db *db = (qhop_db *) qdb;*/
  int i;
  FILE *fp;
  
  fp=ffopen(fn,"w");
/*   for (i ) */
  
  /*print_resall(fp,db->bts,db->nrtp,db->rtp,db->atype,db->bAllDih,
    db->nrexcl,db->bHH14,db->bRemoveDih);*/
  fclose(fp);
  
  return 1;
}

static void clear_qhop_res(qhop_res res)
{
  if (res.acc != NULL) sfree(res.acc);
  if (res.don != NULL) sfree(res.don);
/*   if (res.ft  != NULL) sfree(res.ft); */
}

static void clear_qhop_rb(qhop_resblocks rb)
{
  /* All stings point to a symtab, so don't bother sfreeing them here */
  int rt, r, f;
  for (rt=0; rt<rb.nrestypes; rt++)
    {
      for (r=0; r<rb.nres[rt]; r++)
	{
	  clear_qhop_res(rb.res[rt][r]);
	}
    }
  if (rb.restype != NULL) sfree(rb.restype);
  if (rb.files != NULL) sfree(rb.files);
  if (rb.ilib != NULL) sfree(rb.ilib);
}

int qhop_db_done(qhop_db_t qdb)
{
  int i;
  /* Free a ton of structures. */
  free_t_restp(qdb->nrtp, &(qdb->rtp));
  free_symtab(&qdb->tab);
  qhop_done(*qdb->gqh);
  sfree(qdb->qhop_param);
  clear_qhop_rb(qdb->rb);;
  clear_qhop_H_exist(qdb->H_map);

  return 1;
}

int qhop_db_get_nstates(qhop_db_t qdb,char *resname)
{
  qhop_db_t db = qdb;
  int i,nstate=0;
  
  for(i=0; (i<db->nrtp); i++) {
    if (strncmp(db->rtp[i].resname,resname,3) == 0) {
      nstate++;
    }
  }
  return nstate;
}

int qhop_db_get_qstate(qhop_db_t qdb,char *resname,int istate)
{
  qhop_db_t db = qdb;
  int i,nstate=0;
  printf("qhop_db_get_qstate is called. Thta's no good since it's crippled\n"
	 "and will be replaced by other things.\n");
  for(i=0; (i<db->nrtp); i++) {
    if (strncmp(db->rtp[i].resname,resname,3) == 0) {
      nstate++;
      /*if (nstate == istate)
	return db->resinfo[i].charge;*/
    }
  }
  return nstate;
}

char **qhop_db_get_donors(qhop_db_t qdb,char *resname,int state)
{
  return NULL;
}
 
char **qhop_db_get_acceptors(qhop_db_t qdb,char *resname,int state)
{
  return NULL;
}

/* Obsolete. For now. */
int qhop_db_set_charges(qhop_db_t qdb,char *resname,int state,
			    int natoms,real q[])
{
  return 0;
}

int qhop_db_get_parameters(qhop_db_t qdb,
			       char *donor,char *acceptor,
			       qhop_parameters *qp)
{
  qhop_db_t db = qdb;
  char *aa,*dd;
  int i;
  
  for(i=0; (i<db->ngqh); i++) {
    aa = qhop_get_acceptor(db->gqh[i]);
    dd = qhop_get_donor(db->gqh[i]);
    if (strncmp(donor,dd,3) == 0) {
      if (strncmp(acceptor,aa,3) == 0) {
	memcpy(qp,&(db->qhop_param[i]),sizeof(*qp));
	
	return 1;
      }
    }
  }
  return 0;
}
