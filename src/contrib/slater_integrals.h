/* slater_S_int2.h (c) 2008 Paul J. van Maaren and David van der Spoel
 * Use the functions by indexing with the Slater index - 1
 */

#define SLATER_MAX 5

#ifdef __cplusplus
extern "C"
#endif
double Coulomb_SS(double r,int i,int j,double xi,double xj);

#ifdef __cplusplus
extern "C"
#endif
double Nuclear_SS(double r,int i,double xi);

#ifdef __cplusplus
extern "C"
#endif
double DCoulomb_SS(double r,int i,int j,double xi,double xj);

#ifdef __cplusplus
extern "C"
#endif
double DNuclear_SS(double r,int i,double xi);

