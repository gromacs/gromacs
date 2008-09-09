#ifndef _molprop_xml_h
#define _molprop_xml_h
	
#include "molprop.h"

extern void gmx_molprops_write(char *fn,int nmolprop,gmx_molprop_t mpt[]);

extern gmx_molprop_t *gmx_molprops_read(char *fn,int *nmolprop);

#endif
