#ifndef _interf_h
#define _interf_h

#include "typedefs.h"

typedef struct {
	real Z; /* Interface height-coordinate */
	real t; /*Interface thickness*/
} t_interf;

static void init_interf(t_interf *surf)
{
	surf->Z=0;
	surf->t=0;
}

#endif
