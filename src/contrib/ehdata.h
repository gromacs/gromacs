/*******************************************************************
 *
 * Functions to read the tables
 * The functions must the very first time be called with the file
 * name containing the table. After that, the routines may be called
 * with NULL for a filename. The init_tables routine basically
 * does all the initiation.
 *
 *******************************************************************/	
	
extern real get_omega(real ekin,int *seed,FILE *fp,char *fn);

extern real get_q_inel(real ekin,real omega,int *seed,FILE *fp,char *fn);

extern real get_theta_el(real ekin,int *seed,FILE *fp,char *fn);

extern real cross_inel(real ekin,real rho,char *fn);

extern real cross_el(real ekin,real rho,char *fn);

extern real band_ener(int *seed,FILE *fp,char *fn);

extern void init_tables(int nfile,t_filenm fnm[]);
/* Must be called before any of the table lookup thingies */

extern void test_tables(int *seed,char *fn,real rho);

/*******************************************************************
 *
 * Functions to make histograms
 *
 *******************************************************************/	
	
enum { enormNO, enormFAC, enormNP, enormNR };

typedef struct {
  int np;
  real minx,maxx,dx,dx_1;
  real *y;
  int  *nh;
} t_histo;

extern t_histo *init_histo(int np,real minx,real maxx);

extern void done_histo(t_histo *h);

extern void add_histo(t_histo *h,real x,real y);

extern void dump_histo(t_histo *h,char *fn,char *title,char *xaxis,char *yaxis,
		       int enorm,real norm_fac);


/*******************************************************************
 *
 * Functions to analyse and monitor scattering
 *
 *******************************************************************/	
	
typedef struct {
  int  np,maxp;
  real *time;
  real *ekin;
  bool *bInel;
  rvec *pos;
} t_ana_scat;

extern void add_scatter_event(t_ana_scat *scatter,rvec pos,bool bInel,
			      real t,real ekin);
			      
extern void reset_ana_scat(t_ana_scat *scatter);

extern void done_scatter(t_ana_scat *scatter);

extern void analyse_scatter(t_ana_scat *scatter,t_histo *hmfp);

/*******************************************************************
 *
 * Functions to analyse structural changes
 *
 *******************************************************************/	

typedef struct {
  int  nanal,index;
  real dt;
  real *t;
  real *maxdist;
  real *averdist,*ad2;
  int  *nion;
} t_ana_struct;

extern t_ana_struct *init_ana_struct(int nstep,int nsave,real timestep);

extern void done_ana_struct(t_ana_struct *anal);

extern void reset_ana_struct(t_ana_struct *anal);

extern void add_ana_struct(t_ana_struct *total,t_ana_struct *add);

extern void analyse_structure(t_ana_struct *anal,real t,rvec center,
			      rvec x[],int nparticle,real charge[]);

extern void dump_ana_struct(char *rmax,char *nion,char *gyr,
			    t_ana_struct *anal,int nsim);

/*******************************************************************
 *
 * Functions to analyse energies
 *
 *******************************************************************/	
			    
enum { eCOUL, eREPULS, ePOT, eHOLE, eELECTRON, eLATTICE, eKIN, eTOT, eNR };
extern char *enms[eNR];

typedef real evec[eNR];

typedef struct {
  int  nx,maxx;
  evec *e;
} t_ana_ener;

extern void add_ana_ener(t_ana_ener *ae,int nn,real e[]);

extern void dump_ana_ener(t_ana_ener *ae,int nsim,real dt,char *edump,
			  t_ana_struct *total);
