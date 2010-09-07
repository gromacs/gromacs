#ifndef CUTYPEDEFS_H
#define CUTYPEDEFS_H

#include "cutypedefs_ext.h"

struct cudata 
{
    int         natoms;
    float3 *    f;  /* forces, size natoms */
    float3 *    x;  /* coordinates, size ntypes */

    int     ntypes;
    int *   atom_types; /* atom type indices, size natoms */
    
    float * charges; /* size ntypes */
    float * masses; /* size ntypes */
    float * exclusions;

    float   eps_r;
    float   eps_rf;
    float * nbfp; /* nonbonded parameters C12, C6 */    
};

#endif
