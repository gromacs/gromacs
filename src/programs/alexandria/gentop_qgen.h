/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef GENTOP_QGEN_H
#define GENTOP_QGEN_H

#include <stdio.h>

#include "gmx_resp.h"
#include "poldata.h"

enum {
    eQGEN_OK, eQGEN_NOTCONVERGED, eQGEN_NOSUPPORT, eQGEN_ERROR, eQGEN_NR
};

enum ChargeGenerationAlgorithm {
    eqgNONE, eqgEEM, eqgESP, eqgRESP, eqgNR
};

namespace alexandria
{

class GentopQgen
{

    public:

        ~GentopQgen();


        GentopQgen(const Poldata &pd,
                   t_atoms *atoms,
                   gmx_atomprop_t aps,
                   rvec *x,
                   ChargeDistributionModel   iChargeDistributionModel,
                   ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
                   real hfac, int qtotal,
                   real epsr);

        // void done();

        int generateChargesSm(FILE *fp,
                              const Poldata &pd,
                              t_atoms *atoms,
                              real tol, int maxiter, gmx_atomprop_t aps,
                              real *chieq);

        int generateCharges(FILE *fp,
                            Resp * gr, const std::string molname,
                            const Poldata &pd,
                            t_atoms *atoms,
                            real tol, int maxiter, int maxcycle,
                            gmx_atomprop_t aps);

        void message(int len, char buf[], Resp * gr);
        gmx_bool SplitQ(ChargeDistributionModel iDistributionModel);

        /* The routines below return NOTSET if something is out of the ordinary */
        int getNzeta(int atom);

        int getRow(int atom, int z);

        double getQ(int atom, int z);

        void checkSupport(const Poldata &pd, gmx_atomprop_t aps);

        void saveParams( Resp * gr);

        void getParams( Resp *  gr);
        double getZeta(int atom, int z);


        void print(FILE *fp, t_atoms *atoms);

        void debugFun(FILE *fp);

    private:
        gmx_bool                                           _bWarned;
        ChargeDistributionModel                            _iChargeDistributionModel;
        ChargeGenerationAlgorithm                          _iChargeGenerationAlgorithm;
        int                                                _natom, _eQGEN;
        real                                               _qtotal, _chieq, _hfac, _epsr;
        /* For each atom i there is an elem, atomnr, chi0, rhs, j00 and x */
        std::vector<std::string>                           _elem;
        std::vector<int>                                   _atomnr;
        std::vector<real>                                  _chi0, _rhs, _j00;
        rvec                         *                     _x;
        /* Jab is a matrix over atom pairs */
        std::vector<std::vector<real> >                    _Jab;
        /* For each atom i there are nZeta[i] row, q and zeta entries */
        std::vector<int>                                   _nZeta;
        std::vector<std::vector<int> >                     _row;
        gmx_bool                                           _bAllocSave;
        std::vector<std::vector<real> >                    _q, _zeta, _qsave, _zetasave;


        real calcJab(ChargeDistributionModel iChargeDistributionModel,
                     rvec xi, rvec xj,
                     int nZi, int nZj,
                     std::vector<real> zetaI,
                     std::vector<real> zetaJ,
                     std::vector<int> rowI,
                     std::vector<int> rowJ);

        void calcJab();

        void solveQEem(FILE *fp,  real hardsnesFactor);

        void updateJ00();

        real calcSij(int i, int j);

        void updateFromPoldata(t_atoms *atoms,
                               const Poldata &pd);

        int generateChargesBultinck(FILE *fp,
                                    const Poldata &pd,
                                    t_atoms *atoms,
                                    gmx_atomprop_t aps);

        void calcRhs();
};
}
#endif
