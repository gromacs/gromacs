#ifndef _poldata_h
#define _poldata_h
	
typedef struct gmx_poldata *gmx_poldata_t;

extern gmx_poldata_t gmx_poldata_init();

extern char *gmx_poldata_add_spoel(gmx_poldata_t pd,char *elem,
				   char *miller_equiv,int nhydrogen,
				   int charge,int hybridization,
				   double polarizability,double blength);
				   
extern void gmx_poldata_set_spoel(gmx_poldata_t pd,char *name,
				  double polarizability);
		
extern void gmx_poldata_set_spoel_units(gmx_poldata_t pd,char *polar_unit,char *blength_unit);

/* Returns name or NULL if last or not found */
extern char *gmx_poldata_get_spoel(gmx_poldata_t pd,char *name,char **elem,
				   char **miller_equiv,int *nhydrogen,
				   int *charge,int *hybridization,
				   double *polarizability,double *blength);
				  
extern void gmx_poldata_add_miller(gmx_poldata_t pd,char *name,
				   int atomnumber,
				   double tau_ahc,double alpha_ahp,
				   char *spoel_equiv);

extern void gmx_poldata_set_miller_units(gmx_poldata_t pd,char *tau_unit,char *ahp_unit);

/* Returns name or NULL if last or not found */
extern char *gmx_poldata_get_miller(gmx_poldata_t pd,char *name,
				    int *atomnumber,
				    double *tau_ahc,double *alpha_ahp,
				    char **spoel_equiv);

extern void gmx_poldata_add_bosque(gmx_poldata_t pd,char *elem,
				   double polarizability);
				  
				  
extern void gmx_poldata_set_bosque_units(gmx_poldata_t pd,char *polar_unit);

/* Returns elem or NULL if last or not found */
extern char *gmx_poldata_get_bosque(gmx_poldata_t pd,char *elem,
				    double *polarizability);
				  
#endif
