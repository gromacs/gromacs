#ifndef _sortwater_h
#define _sortwater_h

#include "typedefs.h"

extern void randwater(int astart,int nwater,int nwatom,
		      rvec x[],rvec v[],int *seed);
/* Randomize the order of nwater molecules of length nwatom, the
 * first atom of which is at astart.
 * If v is not NULL it will be shuffled along
 */


extern void sortwater(int nwater,int nwatom,rvec x[]);
/* Sort the order of nwater molecules of length nwatom on X coordinate
 * v not supported yet.
 */

#endif
