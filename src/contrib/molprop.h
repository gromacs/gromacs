#ifndef _molprop_h
#define _molprop_h

#include "atomprop.h"
	
typedef struct gmx_molprop *gmx_molprop_t;

enum { eMOLPROP_Exp, eMOLPROP_Calc, eMOLPROP_Any, eMOLPROP_NR };
	
extern gmx_molprop_t gmx_molprop_init();

extern void gmx_molprop_delete(gmx_molprop_t mpt);

extern void gmx_molprop_set_weight(gmx_molprop_t mpt,double weight);

extern double gmx_molprop_get_weight(gmx_molprop_t mpt);

extern void gmx_molprop_set_formula(gmx_molprop_t mpt,char *formula);

extern char *gmx_molprop_get_formula(gmx_molprop_t mpt);

extern void gmx_molprop_set_molname(gmx_molprop_t mpt,char *molname);

extern char *gmx_molprop_get_molname(gmx_molprop_t mpt);

extern void gmx_molprop_set_reference(gmx_molprop_t mpt,char *reference);

extern char *gmx_molprop_get_reference(gmx_molprop_t mpt);

extern void gmx_molprop_add_composition(gmx_molprop_t mpt,char *compname);

extern void gmx_molprop_reset_composition(gmx_molprop_t mpt,char *compname);

extern char *gmx_molprop_get_composition(gmx_molprop_t mpt);

extern void gmx_molprop_add_composition_atom(gmx_molprop_t mpt,char *compname,
					     char *atomname,int natom);
					     
extern int gmx_molprop_get_composition_atom(gmx_molprop_t mpt,char *compname,
					    char **atomname,int *natom);

extern int gmx_molprop_count_composition_atoms(gmx_molprop_t mpt,
					       char *compname,char *atom);
					       
extern void gmx_molprop_add_property(gmx_molprop_t mpt,int eMP,
				     char *prop_name,double value,double error,
				     char *prop_method,char *prop_reference);

extern int gmx_molprop_get_property(gmx_molprop_t mpt,int *eMP,
				    char **prop_name,double *value,double *error,
				    char **prop_method,char **prop_reference);
			
extern int gmx_molprop_search_property(gmx_molprop_t mpt,int eMP,
				       char *prop_name,double *value,double *error,
				       char *prop_method,char **prop_reference);
	    
extern void gmx_molprop_add_category(gmx_molprop_t mpt,char *category);

/* Returns one category at a time. If NULL, you got them all. */
extern char *gmx_molprop_get_category(gmx_molprop_t mpt);

extern int gmx_molprop_search_category(gmx_molprop_t mpt,char *catname);

extern gmx_molprop_t gmx_molprop_copy(gmx_molprop_t mpt);

extern void gmx_molprop_merge(gmx_molprop_t dst,gmx_molprop_t src);

#endif
