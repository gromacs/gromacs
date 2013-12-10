#ifndef MOLDIP_HPP
#define MOLDIP_HPP

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include "typedefs.h"
#include "mymol.hpp"

typedef struct {
    int  n,nopt,nconst,nopt_c;
    char **name;
    int  *tot_count,*count;
    gmx_bool *bConst;
} t_index_count;

extern char *opt_index_count(t_index_count *ic);

enum { ermsBOUNDS, ermsMU, ermsQUAD, ermsCHARGE, ermsESP, 
       ermsEPOT, ermsForce2, ermsTOT, ermsNR };

namespace alexandria {

class MolDip {
private:
public:
    bool       _bDone,_bFinal,_bGaussianBug;
    bool       _bFitZeta;
    std::vector<alexandria::MyMol> _mymol;
    int            _nmol_support;
    ChargeGenerationModel _iModel;
    t_index_count *_ic;
    real           _J0_0,_Chi0_0,_w_0,_J0_1,_Chi0_1,_w_1;
    real           _hfac,_hfac0,_decrzeta,_epsr;
    real           _ener[ermsNR],_fc[ermsNR];
    gmx_bool       _bOptHfac,_bPol,_bQM;
    char          *_fixchi;
    gmx_poldata_t  _pd;
    gmx_atomprop_t _atomprop;
    t_commrec     *_cr;
    
    //! Constructor
    MolDip();
    
    //! Destructor
    ~MolDip() {};
    
    void Init(t_commrec *cr,gmx_bool bQM,gmx_bool bGaussianBug,
              ChargeGenerationModel iModel,real rDecrZeta,real epsr,
              real J0_0,real Chi0_0,real w_0,
              real J0_1,real Chi0_1,real w_1,
              real fc_bound,real fc_mu,real fc_quad,real fc_charge,
              real fc_esp,real fc_epot,real fc_force,char *fixchi,
              gmx_bool bOptHfac,real hfac,
              gmx_bool bPol,gmx_bool bFitZeta);
    void Read(FILE *fp,const char *fn,const char *pd_fn,
              int minimum_data,
              gmx_bool bZero,
              char *opt_elem,char *const_elem,
              char *lot,
              output_env_t oenv,gmx_molselect_t gms,
              real watoms,gmx_bool bCheckSupport, unsigned int seed);
              
    void CalcDeviation();
};

}

#endif
