/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef MOLDIP_H
#define MOLDIP_H

#include "gromacs/utility/real.h"

#include "mymol.h"

typedef struct {
    int       n, nopt, nconst, nopt_c;
    char    **name;
    int      *tot_count, *count;
    gmx_bool *bConst;
} t_index_count;

extern char *opt_index_count(t_index_count *ic);

enum {
    ermsBOUNDS, ermsMU, ermsQUAD, ermsCHARGE, ermsESP,
    ermsEPOT, ermsForce2, ermsTOT, ermsNR
};

namespace alexandria
{

class MolDip
{
    private:
    public:
        bool                           _bDone, _bFinal, _bGaussianBug;
        bool                           _bFitZeta;
        std::vector<alexandria::MyMol> _mymol;
        int                            _nmol_support;
        ChargeDistributionModel        _iChargeDistributionModel;
        ChargeGenerationAlgorithm      _iChargeGenerationAlgorithm;
        t_index_count                 *_ic;
        real                           _J0_0, _Chi0_0, _w_0, _J0_1, _Chi0_1, _w_1;
        real                           _hfac, _hfac0, _decrzeta;
        real                           _ener[ermsNR], _fc[ermsNR];
        gmx_bool                       _bOptHfac, _bPol, _bQM;
        char                          *_fixchi;
        Poldata                        pd_;
        gmx_atomprop_t                 _atomprop;
        t_commrec                     *_cr;

        //! Constructor
        MolDip();

        //! Destructor
        ~MolDip() {};

        void Init(t_commrec *cr, gmx_bool bQM, gmx_bool bGaussianBug,
                  ChargeDistributionModel iChargeDistributionModel,
                  ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
                  real rDecrZeta,
                  real J0_0, real Chi0_0, real w_0,
                  real J0_1, real Chi0_1, real w_1,
                  real fc_bound, real fc_mu, real fc_quad, real fc_charge,
                  real fc_esp, real fc_epot, real fc_force, char *fixchi,
                  gmx_bool bOptHfac, real hfac,
                  gmx_bool bPol, gmx_bool bFitZeta);

        void Read(FILE *fp, const char *fn, const char *pd_fn,
                  int minimum_data, gmx_bool bZero,
                  char *opt_elem, char *const_elem,
                  char *lot, const MolSelect &gms,
                  real watoms, gmx_bool bCheckSupport,
                  bool bPairs, bool bDihedral, 
		  bool bPolar, const char *tabfn);

        void CalcDeviation();
};

}

#endif
