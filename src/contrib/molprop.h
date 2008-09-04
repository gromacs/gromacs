typedef struct gmx_molprop *gmx_molprop_t;

enum { eMOLPROP_Exp, eMOLPROP_Calc, eMOLPROP_NR };
	
extern gmx_molprop_t gmx_molprop_init();

extern void gmx_molprop_set_weight(gmx_molprop_t mpt,double weight);

extern double gmx_molprop_get_weight(gmx_molprop_t mpt);

extern void gmx_molprop_set_formula(gmx_molprop_t mpt,char *formula);

extern char *gmx_molprop_get_formula(gmx_molprop_t mpt);

extern void gmx_molprop_set_molname(gmx_molprop_t mpt,char *molname);

extern char *gmx_molprop_get_molname(gmx_molprop_t mpt);

extern void gmx_molprop_set_reference(gmx_molprop_t mpt,char *reference);

extern char *gmx_molprop_get_reference(gmx_molprop_t mpt);

extern int gmx_molprop_add_composition(gmx_molprop_t mpt,char *compname);

extern void gmx_molprop_add_composition_atom(gmx_molprop_t mpt,char *compname,
					     char *atomname,int natom);

extern void gmx_molprop_add_property(gmx_molprop_t mpt,int eMP,
				     char *prop_name,double value,double error,
				     char *prop_method,char *prop_reference);

extern int gmx_molprop_get_property(gmx_molprop_t mpt,int *eMP,
				    char **prop_name,double *value,double *error,
				    char **prop_method,char **prop_reference);
				    
extern void gmx_molprop_add_category(gmx_molprop_t mpt,char *category);

/* Returns one category at a time. If NULL */
extern char *gmx_molprop_get_category(gmx_molprop_t mpt);

extern void gmx_molprop_reset_category(gmx_molprop_t mpt);

