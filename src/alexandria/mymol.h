
#ifndef _MYMOL_H
#define _MYMOL_H

#include "typedefs.h"
#include "vsite.h"
#include "atomprop.h"
#include "gmx_resp.h"
#include "gentop_qgen.h"
#include "molprop.h"
#include "molselect.h"
#include "poldata.h"

enum 
{ 
    immOK, immZeroDip, immNoQuad, immCharged, immError, 
    immAtomTypes, immAtomNumber, immMolpropConv,
    immQMInconsistency, immTest, immNR 
};

enum { eSupportNo, eSupportLocal, eSupportRemote, eSupportNR };

typedef struct {
    char           *molname,*lot,*ref;
    int            eSupport;
    int            qtotal,mult,natom,nalloc,nshell;
    real           dip_exp,mu_exp2,dip_err,dip_weight,dip_calc,chieq;
    real           *qESP;
    gmx_mtop_t     mtop;
    gmx_localtop_t *ltop;
    t_symtab       symtab;
    t_inputrec     ir;
    gmx_shellfc_t  shell;
    t_mdatoms      *md;
    t_atoms        *atoms;
    t_forcerec     *fr;
    gmx_vsite_t    *vs;
    gmx_resp_t     gr;
    int            *symmetric_charges;
    rvec           *x,*f,*buf,mu_exp,mu_calc,mu_esp,coq;
    tensor         Q_exp,Q_calc,Q_esp;
    t_state        state;
    matrix         box;
    gentop_qgen_t  qgen;
} t_mymol; 

typedef struct {
    int  n,nopt,nconst,nopt_c;
    char **name;
    int  *tot_count,*count;
    gmx_bool *bConst;
} t_index_count;

enum { ermsBOUNDS, ermsMU, ermsQUAD, ermsCHARGE, ermsESP, ermsTOT, ermsNR };

typedef struct {
    gmx_bool    bDone,bFinal,bGaussianBug,bFitZeta;
    int     nmol,nmol_support,iModel;
    t_mymol *mymol;
    t_index_count *ic;
    real    J0_0,Chi0_0,w_0,J0_1,Chi0_1,w_1,hfac,hfac0,decrzeta,epsr;
    real    ener[ermsNR],fc[ermsNR];
    gmx_bool    bOptHfac,bPol,bQM;
    char    *fixchi;
    gmx_poldata_t  pd;
    gmx_atomprop_t atomprop;
    t_commrec *cr;
} t_moldip;

extern const char *immsg(int imm);

extern int init_mymol(t_mymol *mymol,gmx_molprop_t mp,
                      gmx_bool bQM,char *lot,gmx_bool bZero,
                      gmx_poldata_t pd,gmx_atomprop_t aps,
                      int  iModel,t_commrec *cr,int *nwarn,
                      gmx_bool bCharged,const output_env_t oenv,
                      real th_toler,real ph_toler,
                      real dip_toler,real hfac,gmx_bool bESP,
                      real watoms,real rDecrZeta,gmx_bool bPol,gmx_bool bFitZeta);

extern t_moldip *init_moldip(t_commrec *cr,gmx_bool bQM,gmx_bool bGaussianBug,
                             int  iModel,real rDecrZeta,real epsr,
                             real J0_0,real Chi0_0,real w_0,
                             real J0_1,real Chi0_1,real w_1,
                             real fc_bound,real fc_mu,real fc_quad,real fc_charge,
                             real fc_esp,char *fixchi,
                             gmx_bool bOptHfac,real hfac,
                             gmx_bool bPol,gmx_bool bFitZeta);
                             
extern void read_moldip(t_moldip *md,
                        FILE *fp,const char *fn,const char *pd_fn,
                        int minimum_data,
                        gmx_bool bZero,gmx_bool bWeighted,
                        char *opt_elem,char *const_elem,
                        char *lot,gmx_bool bCharged,
                        output_env_t oenv,gmx_molselect_t gms,
                        real th_toler,real ph_toler,real dip_toler,
                        real watoms,gmx_bool bCheckSupport);

extern char *opt_index_count(t_index_count *ic);

#define gmx_assert(n,m) if (n != m) { gmx_fatal(FARGS,"Variable %s = %d, should have been %d",#n,n,m); }

#endif
