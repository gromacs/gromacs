
#ifndef _MYMOL_H
#define _MYMOL_H

#include "typedefs.h"
#include "vsite.h"
#include "gpp_atomtype.h"
#include "atomprop.h"
#include "gmx_resp.hpp"
#include "gentop_qgen.h"
#include "molprop.hpp"
#include "molselect.h"
#include "poldata.h"

enum immStatus { 
    immOK, immZeroDip, immNoQuad, immCharged, immError, 
    immAtomTypes, immAtomNumber, immMolpropConv, immBondOrder,
    immQMInconsistency, immTest, immNoData, immNR 
};

enum eSupport { eSupportNo, eSupportLocal, eSupportRemote, eSupportNR };

typedef struct {
  /* Primary properties */
    char           *molname,*lot,*ref;
    eSupport       eSupp;
    int            qtotal,mult,natom,nalloc,nshell;
    real           dip_exp,mu_exp2,dip_err,dip_weight,dip_calc,chieq,Hform,Emol,Ecalc,Force2;
    real           *qESP;
    gmx_bool       *bRing;
    int            nbond;
    double         *bondorder;
    gmx_mtop_t     mtop;
    gmx_localtop_t *ltop;
    gpp_atomtype_t atype;
    t_symtab       symtab;
    t_inputrec     ir;
    gmx_shellfc_t  shell;
    gmx_enerdata_t enerd;
    t_mdatoms      *md;
    t_topology     *topology;
    t_excls        *excls;
    char           **smnames;
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

enum { ermsBOUNDS, ermsMU, ermsQUAD, ermsCHARGE, ermsESP, 
       ermsEPOT, ermsForce2, ermsTOT, ermsNR };

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

extern void mv_plists(gmx_poldata_t pd,t_params plist[],gmx_bool bForward);

extern int mk_bonds(gmx_poldata_t pd,t_atoms *atoms,rvec x[],
                    gmx_conect gc,t_params plist[],int nbond[],
                    gmx_bool bRing[],
                    gmx_bool bH14,gmx_bool bAllDihedrals,
                    gmx_bool bRemoveDoubleDihedrals,
                    int nexcl,t_excls **excls,
                    gmx_bool bPBC,matrix box,gmx_atomprop_t aps,real tol,
                    gmx_bool bMovePlists);
		     
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
                        gmx_bool bH14,gmx_bool bAllDihedrals,gmx_bool bRemoveDoubleDihedrals,
                        real watoms,gmx_bool bCheckSupport);

extern char *opt_index_count(t_index_count *ic);

#define gmx_assert(n,m) if (n != m) { gmx_fatal(FARGS,"Variable %s = %d, should have been %d",#n,n,m); }

#endif
