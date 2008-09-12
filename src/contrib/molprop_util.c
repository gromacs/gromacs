#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "smalloc.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"

void genenerate_composition(int nmol,gmx_molprop_t mp[],gmx_poldata_t pd)
{
  int i,j,cnumber,nh,hybrid,charge;
  char *compname,*catom,*miller_equiv,*elem;
  double pol,blength;
  
  for(j=0; (j<nmol); j++) {
    gmx_molprop_delete_composition(mp[j],"miller");
    gmx_molprop_add_composition(mp[j],"miller");
    gmx_molprop_delete_composition(mp[j],"bosque");
    gmx_molprop_add_composition(mp[j],"bosque");
    
    while ((compname = gmx_molprop_get_composition(mp[j])) != NULL) {
      if (strcasecmp(compname,"spoel") == 0) {
	while (gmx_molprop_get_composition_atom(mp[j],&catom,&cnumber) == 1) {
	  if (gmx_poldata_get_spoel(pd,catom,&elem,&miller_equiv,
				    &nh,&charge,&hybrid,&pol,&blength) != NULL) {
	    gmx_molprop_add_composition_atom(mp[j],"miller",miller_equiv,cnumber);
	    gmx_molprop_add_composition_atom(mp[j],"bosque",elem,cnumber);
	    if (nh > 0)
	      gmx_molprop_add_composition_atom(mp[j],"bosque","H",nh*cnumber);
	    sfree(miller_equiv);
	    sfree(elem);
	  }
	  else {
	    fprintf(stderr,"Atom %s not found in poldata\n",catom);
	  }
	  sfree(catom);
	}
      }
      sfree(compname);
    }
  }
}

void generate_formula(int nmol,gmx_molprop_t mp[],gmx_atomprop_t ap)
{
  int  i,j,jj,cnumber,an;
  char formula[1280],number[32];
  char *mform,*compname,*catom;
  int  *ncomp;
  real value;

  for(i=0; (i<nmol); i++) {
    snew(ncomp,110);  
    formula[0] = '\0';
    while ((compname = gmx_molprop_get_composition(mp[i])) != NULL) {
      if (strcasecmp(compname,"bosque") == 0) {
	while (gmx_molprop_get_composition_atom(mp[i],&catom,&cnumber) == 1) {
	  if (gmx_atomprop_query(ap,epropElement,"???",catom,&value)) {
	    an = gmx_nint(value);
	    range_check(an,1,110);
	    ncomp[an]++;
	  }
	  sfree(catom);
	}
      }
    }
    for(jj=0; (jj<2); jj++) {
      for(j=1; (j<110); j++) {
	if ((((jj == 0) && (j == 6)) || ((jj == 1) && (j != 6))) &&
	    (ncomp[j] > 0)) {
	  strcat(formula,gmx_atomprop_element(ap,j));
	  if (ncomp[j] > 1) {
	    sprintf(number,"%d",ncomp[j]);
	    strcat(formula,number);
	  }
	}
      }
    }
    if (debug) {
      mform = gmx_molprop_get_formula(mp[i]);
      if (strcasecmp(formula,mform) != 0) 
	fprintf(debug,"Formula '%s' does match '%s' based on composition for %s.\n",
		mform,formula,gmx_molprop_get_molname(mp[i]));
    }
    gmx_molprop_set_formula(mp[i],formula);
    sfree(compname);
    sfree(ncomp);
  }
}
  
gmx_molprop_t atoms_2_molprop(char *molname,t_atoms*atoms,t_params *bonds,
			      gmx_atomprop_t ap,gmx_poldata_t pd)
{
  gmx_molprop_t mp;
  int ai,aj,i,j,*nh;
  char buf[32],*ptr;
  
  mp = gmx_molprop_init();
  gmx_molprop_set_molname(mp,molname);
  gmx_molprop_add_composition(mp,"spoel");
  
  snew(nh,atoms->nr);
  for(i=0; (i<bonds->nr); i++) {
    ai = bonds->param[i].AI;
    aj = bonds->param[i].AJ;
    if ((*atoms->atomtype[ai])[0] == 'H') {
      nh[aj]++;
    }
    else if ((*atoms->atomtype[aj])[0] == 'H') {
      nh[ai]++;
    }
  }
  for(i=0; (i<atoms->nr); i++) {
    ptr = *(atoms->atomtype[i]);
    if (ptr[0] != 'H') {
      sprintf(buf,"%s%d",ptr,nh[i]);
      gmx_molprop_add_composition_atom(mp,"spoel",buf,1);
    }
  }
  sfree(nh);
  
  genenerate_composition(1,&mp,pd);
  generate_formula(1,&mp,ap);
  
  return mp;  
}

static int comp_mp(const void *a,const void *b)
{
  gmx_molprop_t ma = (gmx_molprop_t)a;
  gmx_molprop_t mb = (gmx_molprop_t)b;
  char *mma = gmx_molprop_get_molname(ma);
  char *mmb = gmx_molprop_get_molname(mb);
  
  if (mma && mmb)
    return strcasecmp(mma,mmb);
  else
    return 0;
}

static int comp_mp_formula(const void *a,const void *b)
{
  int r;
  gmx_molprop_t ma = (gmx_molprop_t)a;
  gmx_molprop_t mb = (gmx_molprop_t)b;
  char *fma = gmx_molprop_get_formula(ma);
  char *fmb = gmx_molprop_get_formula(mb);
  
  r = strcasecmp(fma,fmb);
  
  if (r == 0) 
    return comp_mp(a,b);
  else 
    return r;
}

static void merge_doubles(int *np,gmx_molprop_t mp[],char *doubles)
{
  int i,j,ndouble=0;
  FILE *fp;
  
  fp = fopen(doubles,"w");
  for(i=1; (i<*np); i++) {
    if (strcasecmp(gmx_molprop_get_molname(mp[i]),
		   gmx_molprop_get_molname(mp[i-1])) == 0) {
      if (strcasecmp(gmx_molprop_get_formula(mp[i]),
		     gmx_molprop_get_formula(mp[i-1])) == 0) {
	fprintf(fp,"%5d  %s\n",ndouble+1,
		gmx_molprop_get_molname(mp[i-1]));
	gmx_molprop_merge(mp[i-1],mp[i]);
	for(j=i+1; (j<*np); j++) {
	  gmx_molprop_delete(mp[j-1]);
	  mp[j-1] = gmx_molprop_copy(mp[j]);
	}
	ndouble++;
	(*np)--;
      }
      else {
	printf("Molecules %s, %s have formulae %s resp. %s\n",
	       gmx_molprop_get_molname(mp[i]),
	       gmx_molprop_get_formula(mp[i]));
      }
      
    }
  }
  fclose(fp);
  printf("There were %d double entries\n",ndouble);
}

static void dump_mp(int np,gmx_molprop_t mp[])
{
  FILE *fp;
  int  i,j,k;
  
  fp = fopen("dump_mp.dat","w");
  
  for(i=0; (i<np); ) {
    for(j=i; (j<np-1) && (strcasecmp(gmx_molprop_get_formula(mp[i]),
				     gmx_molprop_get_formula(mp[j+1])) == 0); j++)
      ;
    if (j > i) {
      for(k=i; (k<=j); k++)
	fprintf(fp,"%-20s  %s\n",
		gmx_molprop_get_formula(mp[k]),
		gmx_molprop_get_molname(mp[k]));
      fprintf(fp,"\n");
    }
    i=j+1;
  }
  
  fclose(fp);
}

gmx_molprop_t *merge_xml(int argc,char *argv[],char *outf,
			 char *sorted,char *doubles,int *nmolprop)
{
  gmx_molprop_t *mp=NULL,*mpout=NULL;
  int       i,j,np,npout=0;
  char      buf[100];
  
  for(i=1; (i<argc); i++) {
    mp = gmx_molprops_read(argv[i],&np);
    sprintf(buf,"%s-test",argv[i]);
    gmx_molprops_write(buf,np,mp);
    srenew(mpout,npout+np);
    for(j=0; (j<np); j++)
      mpout[npout+j] = gmx_molprop_copy(mp[j]);
    npout += np;
  }
  gmx_molprops_write("allmols.xml",npout,mpout);
  
  gmx_molprop_sort(npout,mpout);
  
  //qsort(mpout,npout,sizeof(mpout[0]),comp_mp);
  merge_doubles(&npout,mpout,doubles);
      
  if (outf) {
    printf("There are %d entries to store in output file %s\n",npout,outf);
    gmx_molprops_write(outf,npout,mpout);
  }
  if (sorted) {
    qsort(mpout,npout,sizeof(mp[0]),comp_mp_formula);
    gmx_molprops_write(sorted,npout,mpout);
    dump_mp(npout,mpout);
  }
  
  *nmolprop = npout; 
  
  return mpout;
}
