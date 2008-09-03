#include "smalloc.h"
#include "molprop.h"
	
typedef struct {
  int    eMP;
  char   *pname;
  double pvalue,perror;
  char   *pmethod;
  char   *preference;
} t_property;

typedef struct {
  char *cname;
  int  cnumber;
} t_catom;

typedef struct {
  char    *compname;
  int     ncatom;
  t_catom *catom;
} t_composition;

typedef struct {
  char          *formula,*molname,*reference;
  double        weight;
  int           nproperty,nprop_c;
  t_property    *property;
  int           ncomposition,ncomp_c;
  t_composition *composition;
  int           ncategory,ncat_c;
  char          **category;
} gmx_molprop;

gmx_molprop_t gmx_molprop_init()
{
  gmx_molprop *mp;
  
  snew(mp,1);
  
  return mp;
}

void gmx_molprop_set_weight(gmx_molprop_t mpt,double weight)
{
  gmx_molprop *mp = mpt;
  
  mp->weight=weight;
}

double gmx_molprop_get_weight(gmx_molprop_t mpt)
{
  gmx_molprop *mp = mpt;
  
  return mp->weight;
}

void gmx_molprop_set_formula(gmx_molprop_t mpt,char *formula)
{
  gmx_molprop *mp = mpt;
  
  if (mp->formula != NULL)
    sfree(mp->formula);
  mp->formula = strdup(formula);
}

char *gmx_molprop_get_formula(gmx_molprop_t mpt)
{
  gmx_molprop *mp = mpt;
  
  return mp->formula;
}

void gmx_molprop_set_molname(gmx_molprop_t mpt,char *molname)
{
  gmx_molprop *mp = mpt;
  
  if (mp->molname != NULL)
    sfree(mp->molname);
  mp->molname = strdup(molname);
}

char *gmx_molprop_get_molname(gmx_molprop_t mpt)
{
  gmx_molprop *mp = mpt;
  
  return mp->molname;
}

void gmx_molprop_set_reference(gmx_molprop_t mpt,char *reference)
{
  gmx_molprop *mp = mpt;
  
  if (mp->reference != NULL)
    sfree(mp->reference);
  mp->reference = strdup(reference);
}

char *gmx_molprop_get_reference(gmx_molprop_t mpt)
{
  gmx_molprop *mp = mpt;
  
  return mp->reference;
}

int gmx_molprop_add_composition_atom(gmx_molprop_t mpt,char *compname,
				     char *atomname,int natom)
{
  gmx_molprop *mp = mpt;
  int i,j;
  
  for(i=0; (i<mp->ncomposition); i++)
    if(strcmp(compname,mp->composition[i].compname) == 0)
      break;
  if (i == mp->ncomposition) {
    mp->ncomposition++;
    srenew(mp->composition,mp->ncomposition);
    mp->composition[i].compname = strdup(compname);
    mp->composition[i].ncatom = 0;
    mp->composition[i].catom = NULL;
  }
  for(j=0; (j<mp->composition[i].ncatom); j++) {
    if (strdup(mp->composition[i].catom[j].cname,atomname) == 0) {
      mp->composition[i].catom[j].cnumber += natom;
      break;
    }
  }
  if (j == mp->composition[i].ncatom) {
    mp->composition[i].ncatom++;
    srenew(mp->composition[i].catom,mp->composition[i].ncatom);
    mp->composition[i].catom[j].cnumber = 0;
    mp->composition[i].catom[j].cname = strdup(atomname);
  }
  mp->composition[i].catom[j].cnumber += natom;
}

void gmx_molprop_add_property(gmx_molprop_t mpt,int eMP,
			      char *prop_name,double value,double error,
			      char *prop_method,char *prop_reference)
{
  gmx_molprop *mp = mpt;
  t_property  *pp;
  
  mp->nproperty++;
  srenew(mp->property,mp->nproperty);
  pp = &(mp->property[mp->nproperty-1]);
  pp->eMP       = eMP;
  pp->name      = strdup(prop_name);
  pp->value     = value;
  pp->error     = error;
  pp->method    = strdup(prop_method);
  pp->reference = strdup(prop_reference);
}

void gmx_molprop_add_category(gmx_molprop_t,char *category)
{
  gmx_molprop *mp = mpt;
  int i;
  
  for(i=0; (i<mp->ncategory); i++) 
    if (strcasecmp(mp->categories[i],category) == 0)
      break;
  if (i == mp->ncategory) {
    mp->ncategory++;
    mp->categories[i] = strdup(category);
  }
} 

char *gmx_molprop_get_category(gmx_molprop_t)
{
  gmx_molprop *mp = mpt;
  int i;
  
  if (mp->ncat_c < mp->ncategory) {
    mp->ncat_c++;
    return mp->category[mp->ncat_c-1];
  }
  return NULL;
}

void gmx_molprop_reset_category(gmx_molprop_t mpt)
{
  gmx_molprop *mp = mpt;

  mp->ncat_c = 0;
}

