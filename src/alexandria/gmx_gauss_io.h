
typedef struct gau_atomprop *gau_atomprop_t;

#ifdef __cplusplus
extern "C"
#endif
/* read composite method atom data */
gau_atomprop_t read_gauss_data(void);

#ifdef __cplusplus
extern "C"
#endif

void done_gauss_data(gau_atomprop_t gaps);

#ifdef __cplusplus
extern "C"
#endif
gmx_molprop_t gmx_molprop_read_gauss(const char *g98,gmx_bool bBabel,
                                     gmx_atomprop_t aps,gmx_poldata_t pd,
                                     char *molnm,char *iupac,char *conformation,
                                     char *basisset,gau_atomprop_t gaps,
                                     real th_toler,real ph_toler,
                                     int maxpot,gmx_bool bVerbose);

