#include "gmxpre.h"

#include "gromacs/legacyheaders/constr.h"
#include "gromacs/legacyheaders/types/simple.h"

TEST_F(Shake, ConstrainsOneBond)
{
    const int numAtoms = 2;
    atom_id
    cshake(iatom, ncon, nnit, maxnit, constrained_distance_squared,
           xp, rij, half_of_reduced_mass, omega, invmass, twice_the_tolerance, g, nerror);
}
