#ifndef _ionize_h
#define _ionize_h

#include <stdio.h>
#include "typedefs.h"
	
extern void ionize(FILE *log,t_mdatoms *md,char **atomname[],
		   real t,t_inputrec *ir);

#endif
