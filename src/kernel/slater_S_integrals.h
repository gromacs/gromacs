/* slater_S_integrals.h (c) 2008 Paul J. van Maaren and David van der Spoel
 * Use the functions by indexing with the Slater index - 1
 */

typedef double t_slater_SS_func(double rij,double xii,double xij);

typedef double t_slater_NS_func(double rij,double xii);

#define SLATER_MAX 6

extern t_slater_SS_func (*Slater_SS[SLATER_MAX][SLATER_MAX]);

extern t_slater_SS_func (*DSlater_SS[SLATER_MAX][SLATER_MAX]);

extern t_slater_NS_func (*Slater_NS[SLATER_MAX]);

extern t_slater_NS_func (*DSlater_NS[SLATER_MAX]);

