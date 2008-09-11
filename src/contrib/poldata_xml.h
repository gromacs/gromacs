#ifndef _poldata_xml_h
#define _poldata_xml_h

#include "poldata.h"

extern void gmx_poldata_write(char *fn,gmx_poldata_t pd);

extern gmx_poldata_t gmx_poldata_read(char *fn);
			  
#endif
