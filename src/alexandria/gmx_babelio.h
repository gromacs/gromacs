
typedef struct gap {
    char *element, *method,*desc;
    real temp;
    real value;
} gap_t;

typedef struct gau_atomprop
{
    int   ngap;
    gap_t *gap;
} gau_atomprop;

typedef struct gau_atomprop *gau_atomprop_t;

#ifdef __cplusplus
extern "C"
#endif

int gau_atomprop_get_value(gau_atomprop_t gaps,char *element,char *method,
                           char *desc,double temp,double *value);

/* Read a line from a G03/G09 composite method (G3, G4, etc) record */
#ifdef __cplusplus
extern "C"
#endif

int gau_comp_meth_read_line(char *line,real *temp,real *pres);

#ifdef __cplusplus
extern "C"
#endif

gmx_bool read_polar(char *str,tensor T);

#ifdef __cplusplus
extern "C"
#endif

gmx_bool read_dipole(char *str,rvec mu);

#ifdef __cplusplus
extern "C"
#endif

gmx_bool read_quad(char *str1,char *str2,tensor Q);

#ifdef __cplusplus
extern "C"
#endif

gmx_molprop_t gmx_molprop_read_gauss(const char *g98,
                                     gmx_atomprop_t aps,gmx_poldata_t pd,
                                     char *molnm,char *iupac,char *conformation,
                                     gau_atomprop_t gaps,
                                     real th_toler,real ph_toler,
                                     int maxpot,gmx_bool bVerbose);
