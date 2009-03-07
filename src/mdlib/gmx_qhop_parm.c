#include <stdlib.h>
#include <string.h>
#include "smalloc.h"
#include "string2.h"

#include "gmx_qhop_parm.h"	
	
typedef struct {
  char *donor,*acceptor;
  int nparam,nparam_c;
  char **value,**unit,**name;
} gmx_qhop_t;

#define assign_str(dst,src)  if (NULL != src) { if (NULL != dst) *dst = strdup(src); } else { *dst = NULL; }
#define assign_scal(dst,src) if (NULL != dst) *dst = src

/* Return a new gmx_qhop structure */
gmx_qhop gmx_qhop_init()
{
  gmx_qhop_t *qht;
  
  snew(qht,1);

  return (gmx_qhop) qht;
}

void gmx_qhop_set_donor(gmx_qhop gqh,char *donor)
{
  gmx_qhop_t *qht = (gmx_qhop_t *)gqh;

  qht->donor = strdup(donor);
}

void gmx_qhop_set_acceptor(gmx_qhop gqh,char *acceptor)
{
  gmx_qhop_t *qht = (gmx_qhop_t *)gqh;

  qht->acceptor = strdup(acceptor);
}

char *gmx_qhop_get_donor(gmx_qhop gqh)
{
  gmx_qhop_t *qht = (gmx_qhop_t *)gqh;

  return qht->donor;
}

char *gmx_qhop_get_acceptor(gmx_qhop gqh)
{
  gmx_qhop_t *qht = (gmx_qhop_t *)gqh;

  return qht->acceptor;
}

/* Add parameter to gqh, return 1 if OK, 0 if not OK */
int gmx_qhop_add_param(gmx_qhop gqh,char *name,char *value,char *unit)
{
  gmx_qhop_t *qht = (gmx_qhop_t *)gqh;
  
  srenew(qht->name,qht->nparam+1);
  srenew(qht->value,qht->nparam+1);
  srenew(qht->unit,qht->nparam+1);
  qht->name[qht->nparam]  = strdup(name);
  qht->value[qht->nparam] = strdup(value);
  qht->unit[qht->nparam]  = strdup(unit);
  qht->nparam++;
  
  return 1;
}

/* Lists the parameters, one by one on repeatedly calling the
   function. Returns 1 if OK, 0 if not OK */
int gmx_qhop_get_param(gmx_qhop gqh,char **name,char **value,char **unit)
{
  gmx_qhop_t *qht = (gmx_qhop_t *)gqh;

  if (qht->nparam_c < qht->nparam) {
    assign_str(name,qht->name[qht->nparam_c]);
    assign_str(value,qht->value[qht->nparam_c]);
    assign_str(unit,qht->unit[qht->nparam_c]);
    qht->nparam_c++;
    
    return 1;
  }
  else
    qht->nparam_c = 0;
    
  return 0;
}

/* Return a value corresponding to name */
int gmx_qhop_get_value(gmx_qhop gqh,char *name,double *x)
{
  gmx_qhop_t *qht = (gmx_qhop_t *)gqh;
  int i;
  
  for(i=0; (i<qht->nparam); i++) 
    if (gmx_strcasecmp(qht->name[i],name) == 0) {
      *x = atof(qht->value[i]);
      return 1;
    }
    
  return 0;
}

/* Liberate memory */
void gmx_qhop_done(gmx_qhop gqh)
{
  gmx_qhop_t *qht = (gmx_qhop_t *)gqh;
  int i;
  
  for(i=0; (i<qht->nparam); i++) {
    sfree(qht->name[i]);
    sfree(qht->value[i]);
    sfree(qht->unit[i]);
  }
  if (qht->nparam > 0) {
    sfree(qht->name);
    sfree(qht->value);
    sfree(qht->unit);
  }
  if (qht->donor)
    sfree(qht->donor);
  if (qht->acceptor)
    sfree(qht->acceptor);
}

