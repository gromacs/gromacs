#ifndef _sortwater_h
#define _sortwater_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"

extern void randwater(int astart,int nwater,int nwatom,
		      rvec x[],rvec v[],int *seed);
/* Randomize the order of nwater molecules of length nwatom, the
 * first atom of which is at astart.
 * If v is not NULL it will be shuffled along
 */


extern void sortwater(int astart,int nwater,int nwatom,rvec x[],rvec v[]);
/* Sort the order of nwater molecules of length nwatom on X coordinate
 * If v is not NULL it will be shuffled along
 */

extern void mkcompact(int astart,int nwater,int nwatom,rvec x[],rvec v[],
		      int nnode,matrix box);
/* Make compact subboxes */

#endif
