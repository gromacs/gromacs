#ifndef _mkinl_h
#define _mkinl_h

static char *SRCID_mkinl_h = "";
#include <config.h>
#include <types/simple.h>
#include <metacode.h>

#define STRINGSIZE 1024

#define WATERATOMS 3

/* Most options are controlled by these defines.
 * There are two major reasons for this:
 * First, it makes the code somewhat less cluttered.
 * Second, it makes it easier to introduce new interactions and/or
 * change the code without changing all if statements - just set
 * the options in mkinl.c and most things should follow automatically.
 */

#define DO_RF        (loop.coul==COUL_RF)
#define DO_COULTAB     (loop.coul==COUL_TAB)
#define DO_VDWTAB       ((loop.vdw==VDW_TAB) || (loop.vdw==VDW_BHAMTAB))
#define DO_TAB         (DO_COULTAB || DO_VDWTAB)
#define DO_BHAM        ((loop.vdw==VDW_BHAM) || (loop.vdw==VDW_BHAMTAB))

#define DO_PREFETCH_X  ((loop.sol==SOL_NO && (opt.prefetch_x & TYPE_NORMAL)) || \
                        (loop.sol==SOL_MNO && (opt.prefetch_x & TYPE_SOLVENT)) || \
                        (loop.sol==SOL_WATER && (opt.prefetch_x & TYPE_WATER)) || \
                        (loop.sol==SOL_WATERWATER && (opt.prefetch_x & TYPE_WATERWATER)))

#define DO_PREFETCH_F  ((loop.sol==SOL_NO && (opt.prefetch_f & TYPE_NORMAL)) || \
                        (loop.sol==SOL_MNO && (opt.prefetch_f & TYPE_SOLVENT)) || \
                        (loop.sol==SOL_WATER && (opt.prefetch_f & TYPE_WATER)) || \
                        (loop.sol==SOL_WATERWATER && (opt.prefetch_f & TYPE_WATERWATER)))
  
#define DO_SOL         (loop.sol==SOL_MNO) /* non-water solvent optimzations */
#define DO_WATER       (loop.sol==SOL_WATER || loop.sol==SOL_WATERWATER)
#define DO_SOFTCORE    (loop.free==FREE_SOFTCORE)

#define DO_VECTORIZE    (loop.vectorize_invsqrt || loop.vectorize_recip)

/* For pure LJ-only loops we can save some cycles by just calculating the
 * reciprocal. In all other cases we need to do the invsqrt
 */
#define DO_INLINE_INVSQRT   (loop.invsqrt && !loop.vectorize_invsqrt && opt.inline_gmxcode)
#define DO_INLINE_RECIP     (loop.recip && !loop.vectorize_recip && opt.inline_gmxcode)
#define OVERWRITE_RSQ  (!arch.vector && !loop.vdw_needs_r && !loop.vdw_needs_rsq && \
                        !loop.coul_needs_r && !loop.coul_needs_rsq)
			/* Can we overwrite the vectorization array
			 * holding rsq with the resulting rinv or rinvsq,
			 * or do we need to keep rsq a while?
			 */
#define N_VDWPARAM       (DO_BHAM ? 3 : 2)
			
/* Enumerations for different types of innerloops. These might be similar,
 * but not identical, to some of the definitions in include/types/enums.h.
 * Here we also need to take into account whether we use tables or not,
 * and whether a certain innerloop should calculate free energy or not.
 */


/* Electrostatics */
typedef enum {
  COUL_NO,
  COUL_NORMAL,          /* ordinary coulomb, constant epsilon */
  COUL_RF,              /* reaction field */
  COUL_TAB,             /* tabulated interactions */
  COUL_NR
} coul_t;

/* Alternatives for nonbonded interactions */
typedef enum {
  VDW_NO,
  VDW_LJ,              /* Lennard-Jones 6-12 interactions */
  VDW_BHAM,            /* Buckingham */
  VDW_TAB,             /* tabulated interactions */
  VDW_BHAMTAB,         /* Exponential tab, for buckingham interactions */
  VDW_NR
} vdw_t;

typedef enum {
  SOL_NO,             /* no water or solvent optimization */
  SOL_MNO,            /* interactions with solvent consisting of N
                       * lj+coul atoms and M atoms with coul only */ 
  SOL_WATER,          /* interactions with a 3-atom water molecule */
  SOL_WATERWATER,     /* interactions between two water molecules */
  SOL_NR
} sol_t;

typedef enum {
  FREE_NO,
  FREE_LAMBDA,      /* A separate free energy innerloop (lambda)      */
  FREE_SOFTCORE,    /* is called for alpha=0 or lambda=0 or lambda=1, */ 
  FREE_NR           /* to avoid if statements or extra table lookups. */
} free_t;


typedef enum {
  THREAD_NO,
  THREAD_MUTEX,
  THREAD_STRIDE
} thread_t;


/* This is essentially a boolean telling for which loop types
 * some features are enabled.
 * Bit 0 concerns normal loops,
 * bit 1 solvent loops and
 * bit 2 water loops and
 * bit 3 water-water loops.
 */
typedef int looptype_t;

#define TYPE_NORMAL     1
#define TYPE_SOLVENT    2
#define TYPE_WATER      4
#define TYPE_WATERWATER 8


/* Global structure determining architectural options */
typedef struct {
  bool     gmx_invsqrt;    /* Use gmx software routines? */
  bool     gmx_recip;    
  thread_t threads;
  bool     simplewater;
  bool     vector;          
  bool     cray_pragma;
} arch_t;


/* Global structure determining optimization options */
typedef struct {
  bool         inline_gmxcode;            /* Inline gmx software routines */
  looptype_t   prefetch_x;       
  looptype_t   prefetch_f;
  bool         decrease_square_latency;   /* Try to hide latencies */ 
  bool         decrease_lookup_latency;   
  bool         delay_invsqrt;
  /* Dont try to hide latency for (hardware) invsqrt or reciprocal
   * by calculating the before we need them. 
   */
  
  looptype_t   vectorize_invsqrt;         /* vectorize distance calculation when possible */
  bool         vectorize_recip;           /* not meaningful on water loops */
} opt_t;


/* Global structure determining options for the individual
 * loop currently being written
 */ 

typedef struct {
    coul_t   coul;         /* What kind of loop is this? */
    vdw_t     vdw;
    sol_t    sol;
    free_t   free;
    bool     coul_needs_rinv;    /* Which power of r are needed for the */
    bool     coul_needs_rinvsq;  /* coulombic interactions?             */
    bool     coul_needs_rsq;
    bool     coul_needs_r;
    bool     vdw_needs_rinv;    /* And for the nonbonded? */
    bool     vdw_needs_rinvsq;
    bool     vdw_needs_rsq;
    bool     vdw_needs_r;
    bool     invsqrt;           /* tells what we are calculating */
    bool     recip;
    bool     vectorize_invsqrt; 
    bool     vectorize_recip;
    int      ni;
    int      nj;
} loop_t;


extern arch_t arch;
extern opt_t  opt;
extern loop_t loop;
extern bool writeback;

extern char loopname[STRINGSIZE];

extern int table_element_size;

extern char forcebuf[255];
extern char tabforcebuf[255];

/* Exported routines from each file */
 
/* mkinl_declarations.c */
void file_header(void);
/* Adds include files for c, and calls fortran functions to be written
 * at the top of the source file
 */
void func_header(char *appendname);
/* Writes the header with some info and assigns the function name
 * to each loop
 */
void func_args(void);
/* Declares function arguments to the variable buffer */
void func_localvars(void);
/* Declares function variables to the buffer */
void func_init_vars(void);
/* Initiates some local data before the outermost loop */


/* Outer loop routines found in mkinl_outerloop.c
 */
void outer_loop(void);
/* Called to construct the loop over i particles, and in some
 * cases also loops over e.g. threads.
 */


/* Innerloop routines found in mkinl_innerloop.c
 */
void inner_loop(bool calcdist, bool calcforce);
/* Constructs the loop over j particles in a neighbour list */

/* Invsqrt calculation, mkinl_invsqrt.c
 */
int calc_invsqrt(void);
void invsqrt_vars(void);
void fortran_invsqrt(void);
int calc_recip(void);
void recip_vars(void);
void fortran_recip(void);
void fortran_vecinvsqrt(void);
void fortran_vecsqrt(void);


/* The interactions, mkinl_interactions.c
 */

void start_stripmine_loop(void);
void close_stripmine_loop(void);

void call_vectorized_routines(void);

int calc_interactions(bool prefetchx);

int calc_dist(void);
int calc_rinv_and_rinvsq(void);

int calc_rsquare(char *atom1, char *atom2);

void prefetch_forces(void);

void unpack_inner_data(bool calcdist, bool calcforce);

void fetch_coord(bool bPrefetch);

int update_inner_forces(int i,int j);

#endif
