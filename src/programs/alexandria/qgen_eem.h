/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#ifndef QGEN_EEM_H
#define QGEN_EEM_H

#include <cstdio>

#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/state.h"

#include "poldata.h"

struct t_atoms;

enum {
    eQGEN_OK, eQGEN_NOTCONVERGED, eQGEN_NOSUPPORT, eQGEN_ERROR, eQGEN_NR
};

enum ChargeGenerationAlgorithm {
    eqgNONE, eqgEEM, eqgESP
};

namespace alexandria
{

class QgenEem
{
    public:
        QgenEem(const Poldata &pd,
                t_atoms *atoms,
                PaddedRVecVector x,
                ChargeDistributionModel   iChargeDistributionModel,
                double hfac, int qtotal);

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

        const char *message() const;
        
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
        gmx_bool                                           bWarned_;
        ChargeDistributionModel                            _iChargeDistributionModel;
        int                                                _natom, _eQGEN;
        double                                             _qtotal, _chieq, _hfac;
        /* For each atom i there is an elem, atomnr, chi0, rhs, j00 and x */
        std::vector<std::string>                           _elem;
        std::vector<int>                                   _atomnr;
        std::vector<double>                                _chi0, _rhs, _j00;
        std::vector<gmx::RVec>                             _x;
        /* Jab is a matrix over atom pairs */
        std::vector<std::vector<double>>                  _Jab;
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
        void calcRhs();
};
}
#endif
