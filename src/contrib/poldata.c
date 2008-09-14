#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "smalloc.h"
#include "gmx_fatal.h"
#include "poldata.h"

typedef struct {
  char   *name,*elem,*miller_equiv;
  int    nhydrogen,hybridization,charge;
  double polarizability,blength;
} t_spoel;

typedef struct {
  char   *elem;
  double polarizability;
} t_bosque;

typedef struct {
  char   *name,*spoel_equiv;
  int    atomnumber;
  double tau_ahc,alpha_ahp;
} t_miller;

typedef struct {
  char *base;
  char *atom1,*atom2,*atom3;
} t_symcharges;

typedef struct {
  int nspoel,nspoel_c;
  t_spoel *spoel;
  char *spoel_polar_unit,*spoel_blength_unit;
  int nmiller,nmiller_c;
  t_miller *miller;
  char *miller_tau_unit,*miller_ahp_unit;
  int nbosque,nbosque_c;
  t_bosque *bosque;
  char *bosque_polar_unit;
  int nsymcharges,nsymcharges_c;
  t_symcharges *symcharges;
} gmx_poldata;

#define assign_str(dst,src)  if (dst) { if (src) *dst = strdup(src); else *dst = NULL; }
#define assign_scal(dst,src) if (dst) *dst = src

gmx_poldata_t gmx_poldata_init()
{
  gmx_poldata *gpd;
  
  snew(gpd,1);
  
  return (gmx_poldata_t) gpd;
}

char *gmx_poldata_add_spoel(gmx_poldata_t pd,char *elem,
			    char *miller_equiv,
			    int nhydrogen,int charge,
			    int hybridization,
			    double polarizability,double blength)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;
  t_spoel *sp;
  char    buf[32];
  
  gpd->nspoel++;
  srenew(gpd->spoel,gpd->nspoel);
  sp = &(gpd->spoel[gpd->nspoel-1]);
  sp->elem           = strdup(elem);
  sp->miller_equiv   = strdup(miller_equiv);
  sp->nhydrogen      = nhydrogen;
  sp->charge         = charge;
  sp->hybridization  = hybridization;
  sp->polarizability = polarizability;
  sp->blength        = blength;
  if (hybridization > 0) 
    sprintf(buf,"%s%d%d",elem,hybridization,nhydrogen);
  else {
    strcpy(buf,elem);
    if ((charge < -1) || (charge > 1)) 
      sprintf(buf,"%s%d",elem,abs(charge));
    if (charge < 0)
      strcat(buf,"-");
    else if (charge > 0)
      strcat(buf,"+");
  }
  sp->name           = strdup(buf);
  
  return sp->name;
}

void gmx_poldata_set_spoel(gmx_poldata_t pd,char *name,
			   double polarizability)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;
  int i;
  
  for(i=0; (i<gpd->nspoel); i++) 
    if (strcmp(name,gpd->spoel[i].name) == 0) 
      break;
  
  if (i<gpd->nspoel) 
    gpd->spoel[i].polarizability = polarizability;
  else
    gmx_fatal(FARGS,"Could not find atom %s in poldata series",name);
}		

void gmx_poldata_set_spoel_units(gmx_poldata_t pd,char *polar_unit,char *blength_unit)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;

  gpd->spoel_polar_unit   = strdup(polar_unit);
  gpd->spoel_blength_unit = strdup(blength_unit);
}

char *gmx_poldata_get_spoel(gmx_poldata_t pd,char *name,char **elem,
			    char **miller_equiv,int *nhydrogen,
			    int *charge,int *hybridization,
			    double *polarizability,double *blength)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;
  t_spoel *sp;
  int i;
  
  if (name) {
    for(i=0; (i<gpd->nspoel); i++) 
      if (strcmp(name,gpd->spoel[i].name) == 0) 
	break;
  }
  else
    i = gpd->nspoel_c;
      
  if (i<gpd->nspoel) {
    sp = &(gpd->spoel[i]);
    assign_scal(nhydrogen,sp->nhydrogen);
    assign_scal(charge,sp->charge);
    assign_scal(hybridization,sp->hybridization);
    assign_scal(polarizability,sp->polarizability);
    assign_scal(blength,sp->blength);
    assign_str(elem,sp->elem);
    assign_str(miller_equiv,sp->miller_equiv);
    if (!name)
      gpd->nspoel_c++;
    
    return sp->name;
  }
  else
    gpd->nspoel_c = 0;
    
  return NULL;
}
				  
void gmx_poldata_add_miller(gmx_poldata_t pd,char *name,
			    int atomnumber,
			    double tau_ahc,double alpha_ahp,
			    char *spoel_equiv)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;
  t_miller *mil;
  
  gpd->nmiller++;
  srenew(gpd->miller,gpd->nmiller);
  mil = &(gpd->miller[gpd->nmiller-1]);
  mil->name       = strdup(name);
  mil->atomnumber = atomnumber;
  mil->tau_ahc    = tau_ahc;
  mil->alpha_ahp  = alpha_ahp;
  if (spoel_equiv)
    mil->spoel_equiv = strdup(spoel_equiv);
  else
    mil->spoel_equiv = NULL;
}
				  
void gmx_poldata_set_miller_units(gmx_poldata_t pd,char *tau_unit,char *ahp_unit)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;

  gpd->miller_tau_unit = strdup(tau_unit);
  gpd->miller_ahp_unit = strdup(ahp_unit);
}

char *gmx_poldata_get_miller(gmx_poldata_t pd,char *name,
			     int *atomnumber,double *tau_ahc,
			     double *alpha_ahp,char **spoel_equiv)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;
  t_miller *mil;
  int i;
  
  if (name) {
    for(i=0; (i<gpd->nmiller); i++) 
      if (strcmp(name,gpd->miller[i].name) == 0) 
	break;
  }
  else
    i = gpd->nmiller_c;
    
  if (i < gpd->nmiller) {
    mil = &(gpd->miller[i]);
    assign_scal(atomnumber,mil->atomnumber);
    assign_scal(tau_ahc,mil->tau_ahc);
    assign_scal(alpha_ahp,mil->alpha_ahp);
    assign_str(spoel_equiv,mil->spoel_equiv);
    if (!name)
      gpd->nmiller_c++;
    
    return mil->name;
  }
  else
    gpd->nmiller_c = 0;
    
  return NULL;
}

void gmx_poldata_add_bosque(gmx_poldata_t pd,char *elem,
			    double polarizability)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;
  t_bosque *bs;
  
  gpd->nbosque++;
  srenew(gpd->bosque,gpd->nbosque);
  bs = &(gpd->bosque[gpd->nbosque-1]);
  bs->elem           = strdup(elem);
  bs->polarizability = polarizability;
}

char *gmx_poldata_get_bosque(gmx_poldata_t pd,char *elem,
			     double *polarizability)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;
  int i;
  
  if (elem) {
    for(i=0; (i<gpd->nbosque); i++) 
      if (strcmp(elem,gpd->bosque[i].elem) == 0) 
	break;
  }
  else 
    i = gpd->nbosque_c;
    
  if (i < gpd->nbosque) {
    assign_scal(polarizability,gpd->bosque[i].polarizability);
    if (!elem)
      gpd->nbosque_c++;
    
    return gpd->bosque[i].elem;
  }
  else
    gpd->nbosque_c = 0;
    
  return NULL;
}
				  
void gmx_poldata_set_bosque_units(gmx_poldata_t pd,char *polar_unit)
{
  gmx_poldata *gpd = (gmx_poldata *) pd;

  gpd->bosque_polar_unit   = strdup(polar_unit);
}
				  

