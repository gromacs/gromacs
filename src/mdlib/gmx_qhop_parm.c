#include <stdlib.h>
#include <string.h>
#include "smalloc.h"
#include "string2.h"

#include "gmx_qhop_parm.h"	
	
typedef struct gmx_qhop {
  char *donor,*acceptor;
  int nparam,nparam_c;
  char **value,**unit,**name;
} t_gmx_qhop;

#define assign_str(dst,src)  if (NULL != src) { if (NULL != dst) *dst = strdup(src); } else { *dst = NULL; }
#define assign_scal(dst,src) if (NULL != dst) *dst = src

/* Return a new gmx_qhop structure */
gmx_qhop_t gmx_qhop_init()
{
  struct gmx_qhop *qht;
  
  snew(qht,1);

  return qht;
}

void gmx_qhop_set_donor(gmx_qhop_t gqh,char *donor)
{
  gqh->donor = strdup(donor);
}

void gmx_qhop_set_acceptor(gmx_qhop_t gqh,char *acceptor)
{
  gqh->acceptor = strdup(acceptor);
}

char *gmx_qhop_get_donor(gmx_qhop_t gqh)
{
  return gqh->donor;
}

char *gmx_qhop_get_acceptor(gmx_qhop_t gqh)
{
  return gqh->acceptor;
}

/* Add parameter to gqh, return 1 if OK, 0 if not OK */
int gmx_qhop_add_param(gmx_qhop_t gqh,char *name,char *value,char *unit)
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
int gmx_qhop_get_param(gmx_qhop_t gqh,char **name,char **value,char **unit)
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
int gmx_qhop_get_value(gmx_qhop_t gqh,char *name,double *x)
{
  int i;
  
  for(i=0; (i<gqh->nparam); i++) 
    if (gmx_strcasecmp(gqh->name[i],name) == 0) {
      *x = atof(gqh->value[i]);
      return 1;
    }
    
  return 0;
}

/* Liberate memory */
void gmx_qhop_done(gmx_qhop_t gqh)
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

