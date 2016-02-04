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

class QgenEem
{

    public:

        ~QgenEem();


        QgenEem(const Poldata &pd,
                t_atoms *atoms,
                rvec *x,
                ChargeDistributionModel   iChargeDistributionModel,
                double hfac, int qtotal,
                double epsr);

        int generateChargesSm(FILE *fp,
                              const Poldata &pd,
                              t_atoms *atoms,
                              double tol, int maxiter, 
                              double *chieq);

        int generateCharges(FILE *fp,
                            const std::string molname,
                            const Poldata &pd,
                            t_atoms *atoms,
                            double tol, int maxiter);

        void message(int len, char buf[],
                     Resp &gr);
        gmx_bool SplitQ(ChargeDistributionModel iDistributionModel);

        /* The routines below return NOTSET if something is out of the ordinary */
        int getNzeta(int atom);

        int getRow(int atom, int z);

        double getQ(int atom, int z);

        void checkSupport(const Poldata &pd);

        double getZeta(int atom, int z);


        void print(FILE *fp, t_atoms *atoms);

        void debugFun(FILE *fp);

    private:
        gmx_bool                                           _bWarned;
        ChargeDistributionModel                            _iChargeDistributionModel;
        int                                                _natom, _eQGEN;
        double                                             _qtotal, _chieq, _hfac, _epsr;
        /* For each atom i there is an elem, atomnr, chi0, rhs, j00 and x */
        std::vector<std::string>                           _elem;
        std::vector<int>                                   _atomnr;
        std::vector<double>                                _chi0, _rhs, _j00;
        rvec                         *                     _x;
        /* Jab is a matrix over atom pairs */
        std::vector<std::vector<double> >                  _Jab;
        /* For each atom i there are nZeta[i] row, q and zeta entries */
        std::vector<int>                                   _nZeta;
        std::vector<std::vector<int> >                     _row;
        bool                                               _bAllocSave;
        std::vector<std::vector<double> >                  _q, _zeta, _qsave, _zetasave;


        double calcJab(ChargeDistributionModel iChargeDistributionModel,
                       rvec xi, rvec xj,
                       int nZi, int nZj,
                       std::vector<double> zetaI,
                       std::vector<double> zetaJ,
                       std::vector<int> rowI,
                       std::vector<int> rowJ);

        void copyChargesToAtoms(t_atoms *atoms);
        
        void calcJab();

        void solveQEem(FILE *fp, double hardnessFactor);

        void updateJ00();

        double calcSij(int i, int j);

        void updateFromPoldata(t_atoms *atoms,
                               const Poldata &pd);

        int generateChargesBultinck(FILE *fp,
                                    const Poldata &pd,
                                    t_atoms *atoms);
        void saveParams(Resp &gr);

        void restoreParams(Resp &gr);

        void calcRhs();
};
}
#endif
