#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

typedef struct {
  int           nrQMatoms;      /* total nr of QM atoms              */
  rvec          *xQM;           /* shifted to center of box          */  
  int           *indexQM;       /* atom i = atom indexQM[i] in mdrun */
  int           *atomicnumberQM;/* atomic numbers of QM atoms        */  
  real          *QMcharges;     /* atomic charges of QM atoms(ONIOM) */
  int           *shiftQM;
  int           QMcharge;       /* charge of the QM system           */
  int           multiplicity;   /* multipicity (no of unpaired eln)  */
  int           QMmethod;       /* see enums.h for all methods       */
  int           QMbasis;        /* see enums.h for all bases         */
  int           nelectrons;     /* total number of elecs in QM region*/
  bool          bTS;            /* Optimize a TS, only steep, no md  */
  bool          bOPT;          /* Optimize QM subsys, only steep, no md  */
  bool          *frontatoms;   /* qm atoms on the QM side of a QM-MM bond */
  /* Gaussian specific stuff */
  int           nQMcpus;        /* no. of CPUs used for the QM calc. */
  int           QMmem;          /* memory for the gaussian calc.     */
  int           accuracy;       /* convergence criterium (E(-x))     */
  bool          cpmcscf;        /* using cpmcscf(l1003)*/
  char          *gauss_dir;
  char          *gauss_exe;
  char          *devel_dir;
  real          *c6;
  real          *c12;
  /* Surface hopping stuff */
  bool          bSH;            /* surface hopping (diabatic only)   */
  real          SAon;           /* at which energy gap the SA starts */
  real          SAoff;          /* at which energy gap the SA stops  */
  int           SAsteps;        /* stepwise switchinng on the SA     */
  int           SAstep;         /* current state of SA               */
  int           CIdim;
  real          *CIvec1;
  real          *CIvec2;
  real          *CIvec1old;
  real          *CIvec2old;
  ivec          SHbasis;
  int           CASelectrons;
  int           CASorbitals;
} t_QMrec;

typedef struct {
  int           nrMMatoms;      /* nr of MM atoms, updated every step*/
  rvec          *xMM;           /* shifted to center of box          */
  int           *indexMM;       /* atom i = atom indexMM[I] in mdrun */
  real          *MMcharges;     /* MM point charges in std QMMM calc.*/
  int           *shiftMM;
  int           *MMatomtype;    /* only important for semi-emp.      */
  real          scalefactor;
  /* gaussian specific stuff */
  real          *c6;
  real          *c12;
} t_MMrec;


typedef struct {
  int           QMMMscheme; /* ONIOM (multi-layer) or normal          */
  int           nrQMlayers; /* number of QM layers (total layers +1 (MM)) */
  t_QMrec       **qm;        /* atoms and run params for each QM group */
  t_MMrec       *mm;        /* there can only be one MM subsystem !   */
} t_QMMMrec;


