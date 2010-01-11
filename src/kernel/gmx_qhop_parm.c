#include <stdlib.h>
#include <string.h>
#include "smalloc.h"
#include "string2.h"

#include "hackblock.h"
#include "gpp_atomtype.h"
#include "types/gmx_qhop_types.h"
#include "gmx_qhop_parm.h"

#define assign_str(dst,src)  if (NULL != src) { if (NULL != dst) *dst = strdup(src); } else { *dst = NULL; }
#define assign_scal(dst,src) if (NULL != dst) *dst = src

/* Return a new qhop_resblocks structure */
qhop_resblocks_t qhop_resblock_init()
{
  qhop_resblocks *rb;
  snew(rb,1);
  rb->nrestypes = 0;
  return rb;
}

/* Scan files for residues/moleculetypes
 * ff - forcefield
 * ffiles - list of extra itp-files to be scanned
 * nfiles - number of extra files
 * res    - list of residues to look for
 * nres   - size of first dimension of res
 */
qhop_interactions* qhop_scan_rtpitp(char *ff, char **files, int nfiles, char **res, int nres)
{
  int i, ni;
  bool *resfound; /* bolean record of residues found in rtp + extra files. */
  char *dot;
  FILE fp;
  t_symtab symtab;
  qhop_interactions *qi;

  /* rtp-stuff */
  bool       bAlldih,HH14,bRemoveDih;
  int        nrexcl;
  gpp_atomtype_t atype;
  int        bts[ebtsNR];
  t_restp    *restp;
  t_aa_names *aan;

  /* itp-stuff */


  printf("Scanning rtp and itp-files...\n"
	 "  Will store interaction parameters for protonation events\n"
	 "  Can not do this for termini just yet.\n");

  snew(qi,nres);

  open_symtab(&symtab);
  aan = get_aa_names();

  snew(resfound, nres);
  for (i=0; i<nres; i++)
    resfound = 0;

  /* -----------------------------------------.
   * First scan the rtp.                       \
   */
  read_resall(ff, bts, restp, atype, &symtab, &bAlldih, &nrexcl, &HH14, &bRemoveDih);

  /* We'll also need some atomtype data. */
  /*  -  read ffXXXnb.itp */
  /*  -  read ffXXXbon.itp */

  /* rtp scanned.                              /
   * -----------------------------------------'
   */

  /* -----------------------------------------.
   * Now scan extra itp-files                  \
   */
  for (i=0; i<nfiles; i++){
    /* determine filetype */
    if (dot = rindex(files[i],'.'))
      if (!strcmp(".itp", dot))
	{
	}
      else
	dot = NULL;

    if (!dot)
      gmx_fatal(FARGS, "QHOP: %s appears not to be an itp-file, by judging from its name.", files[i]);
  }
  /* itp-files scanned.                        /
   * -----------------------------------------'
   */
  
  return qi;
}

/* Add a res to a restype. Requires a set of residues collected in a qhop_res_t. */
void qhop_add_res(qhop_resblocks_t rb, int resblocknr, qhop_res_t res, int nres)
{
  int i,n;
  if (!rb->res)
    snew(rb->res, 1);
  else
    /* Let's not assume that we're adding stuff to a new resblock. */
    if (resblocknr > rb->nrestypes+1)
      srenew(rb->res,rb->nrestypes+1);
  
  /* The following allows for addition of residues to a non-empty resblock. */
  n = rb->nres[resblocknr];
  if (n > 0)
    n--;

  for (i=0; i<nres; i++) {
    rb->res[resblocknr][i + n] = res[i];
    rb->nres[resblocknr]++;
  }
}

/* Add a restype. Requires a set of residues collected in a qhop_res_t. */
void qhop_add_restype(qhop_resblocks_t rb, char *name, int nres, qhop_res_t res)
{
  int i, j;
  if (!rb)
    gmx_fatal(FARGS,"The resblock is not initialized!\n");

  /* (re)allocate restype */
  if (!rb->restype) {
    snew(rb->restype,1);
    i = 1;
  } else {
    i = rb->nrestypes+1;
    srenew(rb->restype, i);
  }

  /* (re)allocate nres */
  if (!rb->nres)
    snew(rb->nres, 1);
  else
    srenew(rb->nres, i);

  /*rb->nres[rb->nrestypes] = nres;*/

  qhop_add_res(rb, rb->nrestypes, res, nres);

  rb->nrestypes++;
}


void qhop_set_protonation(qhop_resblocks_t rb, qhop_res_t res)
{
  int i;
}

/* Return a new qhop structure */
qhop_t qhop_init()
{
  struct qhop *qht;
  
  snew(qht,1);

  return qht;
}

void qhop_set_donor(qhop_t gqh,char *donor)
{
  if (donor != NULL) {
    printf(" donor %s,", donor);
    gqh->donor = strdup(donor);
  }
}

void qhop_set_acceptor(qhop_t gqh,char *acceptor)
{
  if (acceptor != NULL) {
    printf(" acceptor %s\n", acceptor);
    gqh->acceptor = strdup(acceptor);
  }
}

void qhop_set_don_atom(qhop_t gqh,char *don_atom)
{
  if (don_atom != NULL)
    {
      printf("Setting don_atom %s\n", don_atom);
      gqh->don_atom = strdup(don_atom);
    }
}

void qhop_set_acc_atom(qhop_t gqh,char *acc_atom)
{
  if (acc_atom != NULL)
    {
      printf("Setting acc_atom %s\n", acc_atom);
      gqh->acc_atom = strdup(acc_atom);
    }
}


char *qhop_get_donor(qhop_t gqh)
{
  return gqh->donor;
}

char *qhop_get_acceptor(qhop_t gqh)
{
  return gqh->acceptor;
}

/* Add parameter to gqh, return 1 if OK, 0 if not OK */
int qhop_add_param(qhop_t gqh,char *name,char *value,char *unit)
{
  srenew(gqh->name,gqh->nparam+1);
  srenew(gqh->value,gqh->nparam+1);
  srenew(gqh->unit,gqh->nparam+1);
  gqh->name[gqh->nparam]  = strdup(name);
  gqh->value[gqh->nparam] = strdup(value);
  gqh->unit[gqh->nparam]  = strdup(unit);
  gqh->nparam++;
  
  return 1;
}

/* Lists the parameters, one by one on repeatedly calling the
   function. Returns 1 if OK, 0 if not OK */
int qhop_get_param(qhop_t gqh,char **name,char **value,char **unit)
{
  if (gqh->nparam_c < gqh->nparam) {
    assign_str(name,gqh->name[gqh->nparam_c]);
    assign_str(value,gqh->value[gqh->nparam_c]);
    assign_str(unit,gqh->unit[gqh->nparam_c]);
    gqh->nparam_c++;
    
    return 1;
  }
  else
    gqh->nparam_c = 0;
    
  return 0;
}

/* Return a value corresponding to name */
int qhop_get_value(qhop_t gqh,char *name,double *x)
{
  int i;
  
  for(i=0; (i<gqh->nparam); i++) 
    if (gmx_strcasecmp(gqh->name[i],name) == 0) {
      *x = strtod(gqh->value[i],NULL);
      return 1;
    }
    
  return 0;
}

/* Liberate memory */
void qhop_done(qhop_t gqh)
{
  int i;
  
  for(i=0; (i<gqh->nparam); i++) {
    sfree(gqh->name[i]);
    sfree(gqh->value[i]);
    sfree(gqh->unit[i]);
  }
  if (gqh->nparam > 0) {
    sfree(gqh->name);
    sfree(gqh->value);
    sfree(gqh->unit);
  }
  if (gqh->donor)
    sfree(gqh->donor);
  if (gqh->acceptor)
    sfree(gqh->acceptor);
}

