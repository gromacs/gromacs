/* slater_integrals.cpp (c) 2008 Paul J. van Maaren and David van der Spoel */
#define SLATER_MAX 6

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

