typedef struct {
  int  np;          /* Number of points           */
  int  atype;       /* Atom type                  */
  int  ptype;       /* Parameter type             */
  real rmin,rmax;   /* Minimum and maximum value  */
  real dr;
  real rval;        /* Current value              */
} t_range;

typedef struct {
  int  NP,D;
  int  strategy;
  int  seed;
  int  ipop,gen;        /* Current population member and generation */
  int  imin;            /* Member with lowest energy */
  real CR,FF;
  real **pold,**pnew;   /* Old and new populations */
  real *best,*bestit,*cost,*tmp,*rmsf,*energy;
} t_genalg;

enum { eseSIGMA, eseEPSILON, eseBHAMA, eseBHAMB, eseBHAMC, eseNR };

extern real value_rand(t_range *r,int *seed);

extern t_genalg *init_ga(char *infile,int D,t_range range[]);

extern void update_ga(FILE *fpout_ptr,t_range range[],t_genalg *ga);

extern bool print_ga(FILE *fp,t_genalg *ga,real rmsf,real energy,t_range range[],
		     real tol);

extern real cost(real rmsf,real energy);
