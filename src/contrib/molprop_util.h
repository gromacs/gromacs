#ifndef _molprop_util_h
#define _molprop_util_h
	
#include "grompp.h"
#include "atomprop.h"
#include "molprop.h"
#include "poldata.h"

extern void generate_composition(int nmol,gmx_molprop_t mp[],gmx_poldata_t pd,
				 gmx_atomprop_t ap);

extern void check_formula(int np,gmx_molprop_t pd[],gmx_atomprop_t ap);

extern gmx_molprop_t atoms_2_molprop(char *molname,t_atoms*atoms,t_params *bonds,
				     gmx_atomprop_t ap,gmx_poldata_t pd);
				     
extern gmx_molprop_t *merge_xml(int argc,char *argv[],char *outf,
				char *sorted,char *doubles,int *nmolprop);

#endif
