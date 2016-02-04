/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "gentop_qgen.h"

#include <ctype.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/random/random.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/txtdump.h"

#include "gmx_resp.h"
#include "molprop.h"
#include "poldata.h"
#include "coulombintegrals/coulombintegrals.h"

namespace alexandria
{

QgenEem::QgenEem(const Poldata &pd, 
                       t_atoms *atoms, 
                       rvec *x,
                       ChargeDistributionModel   iChargeDistributionModel,
                       double hfac, 
                       int  qtotal, 
                       double epsr)
{
    _bWarned    = false;
    _bAllocSave = false;
    _eQGEN      = eQGEN_OK;
    _chieq      = 0;
    std::string  atp;
    bool         bSupport = true;
    int          i, j, k, atm;

    _iChargeDistributionModel   = iChargeDistributionModel;
    _hfac                       = hfac;
    _qtotal                     = qtotal;
    _epsr                       = std::max(1.0, epsr);

    _natom      = 0;
    for (i = j = 0; (i < atoms->nr); i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            _natom++;
        }
    }

    _chi0.resize(_natom, 0);
    _rhs.resize(_natom+1, 0);
    _elem.resize(_natom);
    _atomnr.resize(_natom, 0);
    _row.resize(_natom);
    _Jab.resize(_natom+1);
    _zeta.resize(_natom);
    _j00.resize(_natom, 0);
    _q.resize(_natom+1);
    _nZeta.resize(_natom+1, 0);

    snew(_x, _natom);

    /* Special case for chi_eq */
    _nZeta[_natom] = 1;
    _q[_natom].resize(_nZeta[_natom], 0);

    for (i = j = 0; (i < atoms->nr) && bSupport; i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            _Jab[j].resize(_natom+1, 0);
            atm = atoms->atom[i].atomnumber;
            if (atm == 0)
            {
                gmx_fatal(FARGS, "Don't know atomic number for %s %s",
                          *(atoms->resinfo[i].name),
                          *(atoms->atomname[j]));
            }
            atp.assign(*atoms->atomtype[j]);
            if (!pd.haveEemSupport(_iChargeDistributionModel, atp, TRUE))
            {
                fprintf(stderr, "No charge distribution support for atom %s, model %s\n",
                        *atoms->atomtype[j], 
                        getEemtypeName(_iChargeDistributionModel));
                bSupport = false;
            }
            if (bSupport)
            {
                _elem[j]   = atp;
                _atomnr[j] = atm;
                int nz     = pd.getNzeta(_iChargeDistributionModel, atp);
                _nZeta[j]  = nz;
                _q[j].resize(nz, 0);
                _zeta[j].resize(nz, 0);
                _row[j].resize(nz, 0);
                for (k = 0; (k < nz); k++)
                {
                    _q[j][k]    = pd.getQ(_iChargeDistributionModel, atp, k);
                    _zeta[j][k] = pd.getZeta(_iChargeDistributionModel, atp, k);
                    _row[j][k]  = pd.getRow(_iChargeDistributionModel, atp, k);
                    char buf[256];
                    snprintf(buf, sizeof(buf), "Row (in the periodic table) should be at least 1. Here: atype = %s q = %g zeta = %g row = %d model = %s",
                             atp.c_str(), _q[j][k], _zeta[j][k], _row[j][k], 
                             getEemtypeName(_iChargeDistributionModel));
                    GMX_RELEASE_ASSERT(iChargeDistributionModel == eqdAXp ||
                                       _row[j][k] != 0, buf);
                    if (_row[j][k] > SLATER_MAX)
                    {
                        if (debug)
                        {
                            fprintf(debug, "Can not handle higher slaters than %d for atom %s %s\n",
                                    SLATER_MAX,
                                    *(atoms->resinfo[i].name),
                                    *(atoms->atomname[j]));
                        }
                        _row[j][k] = SLATER_MAX;
                    }
                }
                _chi0[j]  = 0;
                _j00[j]   = 0;
                copy_rvec(x[i], _x[j]);
                j++;
            }
        }
    }
}

QgenEem::~QgenEem()
{
    sfree(_x);
}

void QgenEem::saveParams(Resp &gr)
{
    int i, j;
    if (!_bAllocSave)
    {
        _qsave.resize(_natom);
        _zetasave.resize(_natom);
    }
    for (i = 0; (i < _natom); i++)
    {
        if (!_bAllocSave)
        {
            _qsave[i].resize(_nZeta[i], 0);
            _zetasave[i].resize(_nZeta[i], 0);
        }
        for (j = 0; (j < _nZeta[i]); j++)
        {
            if (gr.nAtom() > 0)
            {
                _q[i][j]    = (double)gr.getCharge(i, j);
                _zeta[i][j] = gr.getZeta(i, j);
                if (j == _nZeta[i]-1)
                {
                    _q[i][j] = gr.getAtomCharge(i);
                }
            }
            _qsave[i][j]    = _q[i][j];
            _zetasave[i][j] = _zeta[i][j];
        }
    }
    _bAllocSave = true;
}

void QgenEem::restoreParams(Resp &gr)
{
    int i, j;

    if (_bAllocSave)
    {
        for (i = 0; (i < _natom); i++)
        {
            for (j = 0; (j < _nZeta[i]); j++)
            {
                _q[i][j]    = _qsave[i][j];
                _zeta[i][j] = _zetasave[i][j];
                if (static_cast<int>(gr.nAtom()) == _natom)
                {
                    gr.setCharge(i, j, _q[i][j]);
                    gr.setZeta(i, j, _zeta[i][j]);
                }
            }
        }
    }
    else
    {
        fprintf(stderr, "WARNING: no ESP charges generated.\n");
    }
}

int QgenEem::getNzeta( int atom)
{
    if ((0 <= atom) && (atom < _natom))
    {
        return _nZeta[atom];
    }
    return 0;

}

int QgenEem::getRow( int atom, int z)
{
    if ((0 <= atom) && (atom < _natom) &&
        (0 <= z) && (z <= _nZeta[atom]))
    {
        return _row[atom][z];
    }
    return 0;

}
double QgenEem::getQ(int atom, int z)
{
    if ((0 <= atom) && (atom < _natom) &&
        (0 <= z) && (z <= _nZeta[atom]))
    {
        return _q[atom][z];
    }
    return 0;

}


double QgenEem::getZeta(int atom, int z)
{
    if ((0 <= atom) && (atom < _natom) &&
        (0 <= z) && (z <= _nZeta[atom]))
    {
        return _zeta[atom][z];
    }
    return 0;
}

double CoulombNN(double r)
{
    return 1/r;
}

double QgenEem::calcJab(ChargeDistributionModel iChargeDistributionModel,
                         rvec xi, rvec xj,
                         int nZi, int nZj,
                         std::vector<double> zetaI, std::vector<double> zetaJ,
                         std::vector<int> rowI, std::vector<int> rowJ)
{
    int  i, j;
    rvec dx;
    double r;
    double eTot = 0;

    rvec_sub(xi, xj, dx);
    r = norm(dx);
    if (r == 0)
    {
        gmx_fatal(FARGS, "Zero distance between atoms!\n");
    }
    if ((zetaI[0] <= 0) || (zetaJ[0] <= 0))
    {
        iChargeDistributionModel = eqdAXp;
    }
    switch (iChargeDistributionModel)
    {
        case eqdAXp:
            eTot = CoulombNN(r);
            break;
        case eqdAXs:
        case eqdRappe:
        case eqdYang:
            eTot = 0;
            for (i = nZi-1; (i < nZi); i++)
            {
                for (j = nZj-1; (j < nZj); j++)
                {
                    eTot += Coulomb_SS(r, rowI[i], rowJ[j], zetaI[i], zetaJ[j]);
                }
            }
            break;
        case eqdAXg:
            eTot = 0;
            for (i = nZi-1; (i < nZi); i++)
            {
                for (j = nZj-1; (j < nZj); j++)
                {
                    eTot += Coulomb_GG(r, zetaI[i], zetaJ[j]);
                }
            }
            break;
        default:
            gmx_fatal(FARGS, "Unsupported model %d in calc_jab", iChargeDistributionModel);
    }

    return ONE_4PI_EPS0*(eTot)/ELECTRONVOLT;
}

void QgenEem::solveQEem(FILE *fp,  double hardnessFactor)
{
    double **a, qtot, q;
    int      i, j, n;

    n = _natom+1;
    a = alloc_matrix(n, n);
    for (i = 0; (i < n-1); i++)
    {
        for (j = 0; (j < n-1); j++)
        {
            a[i][j] = _Jab[i][j];
        }
        a[i][i] = hardnessFactor*_Jab[i][i];
    }
    for (j = 0; (j < n-1); j++)
    {
        a[n-1][j] = 1;
    }
    for (i = 0; (i < n-1); i++)
    {
        a[i][n-1] = -1;
    }
    a[n-1][n-1] = 0;

    if (matrix_invert(fp, n, a) == 0)
    {
        for (i = 0; (i < n); i++)
        {
            q = 0;
            for (j = 0; (j < n); j++)
            {
                q += a[i][j]*_rhs[j];
            }
            this->_q[i][_nZeta[i]-1] = q;
            if (fp)
            {
                fprintf(fp, "%2d _RHS = %10g Charge= %10g\n", i, _rhs[i], q);
            }
        }
    }
    else
    {
        for (i = 0; (i < n); i++)
        {
            this->_q[i][_nZeta[i]] = 1;
        }
    }
    _chieq      = this->_q[n-1][_nZeta[n-1]-1];
    qtot        = 0;
    for (i = 0; (i < n-1); i++)
    {
        for (j = 0; (j < _nZeta[i]); j++)
        {
            qtot += this->_q[i][j];
        }
    }

    if (fp && (fabs(qtot - _qtotal) > 1e-2))
    {
        fprintf(fp, "qtot = %g, it should be %g\n", qtot, _qtotal);
    }
    free_matrix(a);
}

void QgenEem::updateJ00()
{
    int    i;
    double j0, qq;
    double zetaH = 1.0698;

    for (i = 0; (i < _natom); i++)
    {
        j0 = _j00[i]/_epsr;
        if (((_iChargeDistributionModel == eqdYang) ||
             (_iChargeDistributionModel == eqdRappe)) &&
            (_atomnr[i] == 1))
        {
            qq = _q[i][_nZeta[i]-1];
            j0 = (1+qq/zetaH)*j0;

            if (debug && (j0 < 0) && !_bWarned)
            {
                fprintf(debug, "WARNING: _J00 = %g for atom %d. The equations will be instable.\n", j0, i+1);
                _bWarned = TRUE;
            }
        }
        _Jab[i][i] = (j0 > 0) ? j0 : 0;
    }
}

void QgenEem::debugFun(FILE *fp)
{
    int i, j;

    for (i = 0; (i < _natom); i++)
    {
        fprintf(fp, "THIS: i: %2d _chi0: %8g J0: %8g q:",
                i+1, _chi0[i], _Jab[i][i]);
        for (j = 0; (j < _nZeta[i]); j++)
        {
            fprintf(fp, " %8g", _q[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "qgen _Jab matrix:\n");
    for (i = 0; (i < _natom); i++)
    {
        for (j = 0; (j <= i); j++)
        {
            fprintf(fp, "  %6.2f", _Jab[i][j]);
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
}

double QgenEem::calcSij(int i, int j)
{
    double dist, dism, Sij = 1.0;
    rvec dx;
    int  l, m, tag;

    rvec_sub(_x[i], _x[j], dx);
    dist = norm(dx);
    if ((dist < 0.118) && (_atomnr[i] != 1) && (_atomnr[j] != 1))
    {
        Sij = Sij*1.64;
    }
    else if ((dist < 0.122) && (_atomnr[i] != 1) && (_atomnr[j] != 1))
    {
        if ((_atomnr[i] != 8) && (_atomnr[j] != 8))
        {
            Sij = Sij*2.23;
        }
        else
        {
            Sij = Sij*2.23;
        }
    }
    else if (dist < 0.125)
    {
        tag = 0;
        if ((_atomnr[i] == 6) && (_atomnr[j] == 8))
        {
            tag = i;
        }
        else if ((_atomnr[i] == 8) && (_atomnr[j] == 6))
        {
            tag = j;
        }
        if (tag != 0)
        {
            printf("found CO\n");
            for (l = 0; (l < _natom); l++)
            {
                if (_atomnr[l] == 1)
                {
                    printf("found H\n");
                    dism = 0.0;
                    for (m = 0; (m < DIM); m++)
                    {
                        dism = dism+gmx::square(_x[tag][m]-_x[l][m]);
                    }

                    printf("dist: %8.3f\n", sqrt(dism));
                    if (sqrt(dism) < 0.105)
                    {
                        printf("dist %5d %5d %5s  %5s %8.3f\n",
                               i, l, _elem[tag].c_str(), _elem[l].c_str(), sqrt(dism));
                        Sij = Sij*1.605;
                    }
                }
            }
        }
    }
    else if ((_atomnr[i] == 6) && (_atomnr[j] == 8))
    {
        Sij = Sij*1.03;
    }
    else if (((_atomnr[j] == 6) && (_atomnr[i] == 7) && (dist < 0.15)) ||
             ((_atomnr[i] == 6) && (_atomnr[j] == 7) && (dist < 0.15)))
    {
        if (_atomnr[i] == 6)
        {
            tag = i;
        }
        else
        {
            tag = j;
        }
        for (l = 0; (l < _natom); l++)
        {
            if (_atomnr[l] == 8)
            {
                printf("found Oxy\n");
                dism = 0.0;
                for (m = 0; (m < DIM); m++)
                {
                    dism = dism+gmx::square(_x[tag][m]-_x[l][m]);
                }
                if (sqrt(dism) < 0.130)
                {
                    printf("found peptide bond\n");
                    Sij = Sij*0.66;
                }
                else
                {
                    Sij = Sij*1.1;
                }
            }
        }
    }
    return Sij;
}

void QgenEem::calcJab()
{
    int    i, j;
    double Jab;

    for (i = 0; (i < _natom); i++)
    {
        for (j = i+1; (j < _natom); j++)
        {
            Jab = calcJab(_iChargeDistributionModel,
                          _x[i], _x[j],
                          _nZeta[i], _nZeta[j],
                          _zeta[i], _zeta[j],
                          _row[i], _row[j]);
            if (_iChargeDistributionModel == eqdYang)
            {
                Jab = Jab*calcSij(i, j);
            }
            _Jab[j][i] = _Jab[i][j] = Jab/_epsr;
        }
    }
}

void QgenEem::calcRhs()
{
    int    i, j, k, l;
    rvec   dx;
    double   r, j1, j1q, qcore;

    /* This right hand side is for all models */
    for (i = 0; (i < _natom); i++)
    {
        _rhs[i] = -_chi0[i];
    }
    _rhs[_natom] = _qtotal;

    /* In case the charge is split in nuclear charge and electronic charge
     * we need to add some more stuff. See paper for details.
     */
    for (i = 0; (i < _natom); i++)
    {
        j1q   = 0;
        qcore = 0;
        for (k = 0; (k < _nZeta[i]-1); k++)
        {
            j1q   += _j00[i]*_q[i][k];
            qcore += _q[i][k];
        }
        j1 = 0;
        /* This assignment for k is superfluous because of the previous loop,
         * but if I take it out it will at some stage break the loop below where
         * exactly this value of k is needed.
         */
        k  = _nZeta[i]-1;
        for (j = 0; (j < _natom); j++)
        {
            if (i != j)
            {
                rvec_sub(_x[i], _x[j], dx);
                r = norm(dx);
                switch (_iChargeDistributionModel)
                {
                    case eqdAXs:
                    case eqdRappe:
                    case eqdYang:
                        for (l = 0; (l < _nZeta[j]-1); l++)
                        {
                            j1 += _q[j][l]*Coulomb_SS(r, k, l, _zeta[i][k], _zeta[j][l]);
                        }
                        break;
                    case eqdAXg:
                        for (l = 0; (l < _nZeta[j]-1); l++)
                        {
                            j1 += _q[j][l]*Coulomb_GG(r, _zeta[i][k], _zeta[j][l]);
                        }
                        break;
                    default:
                        break;
                }
            }
        }
        _rhs[i]           -= j1q + ONE_4PI_EPS0*j1/ELECTRONVOLT;
        _rhs[_natom]      -= qcore;
    }
}

int atomicnumber2rowXX(int elem)
{
    int row;

    /* Compute which row in the periodic table is this element */
    if (elem <= 2)
    {
        row = 1;
    }
    else if (elem <= 10)
    {
        row = 2;
    }
    else if (elem <= 18)
    {
        row = 3;
    }
    else if (elem <= 36)
    {
        row = 4;
    }
    else if (elem <= 54)
    {
        row = 5;
    }
    else if (elem <= 86)
    {
        row = 6;
    }
    else
    {
        row = 7;
    }

    return row;
}

void QgenEem::copyChargesToAtoms(t_atoms *atoms)
{
    int j;
    for (int i = j= 0; (i < atoms->nr); i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            double qq = 0;
            for (int k = 0; (k < _nZeta[j]); k++)
            {
                qq += _q[j][k];
            }
            atoms->atom[i].q = atoms->atom[i].qB = qq;
            j++;
        }
    }
} 

void QgenEem::print(FILE *fp, t_atoms *atoms)
{
    int  i, j;
    rvec mu = { 0, 0, 0 };

    if (_eQGEN == eQGEN_OK)
    {
        if (fp)
        {
            fprintf(fp, "Res  Atom   Nr       J0     _chi0 row        q zeta (1/nm)\n");
        }
        for (i = j = 0; (i < atoms->nr); i++)
        {
            if (atoms->atom[i].ptype == eptAtom)
            {
                for (int m = 0; (m < DIM); m++)
                {
                    mu[m] += atoms->atom[i].q * _x[i][m] * ENM2DEBYE;
                }
                if (fp)
                {
                    fprintf(fp, "%4s %4s%5d %8g %8g",
                            *(atoms->resinfo[atoms->atom[i].resind].name),
                            *(atoms->atomname[i]), i+1, _j00[j], _chi0[j]);
                    for (int k = 0; (k < _nZeta[j]); k++)
                    {
                        fprintf(fp, " %3d %8.5f %8.4f", _row[j][k], _q[j][k],
                                _zeta[j][k]);
                    }
                    fprintf(fp, "\n");
                }
                j++;
            }
        }
        if (fp)
        {
            fprintf(fp, "<chieq> = %10g\n|mu| = %8.3f ( %8.3f  %8.3f  %8.3f )\n",
                    _chieq, norm(mu), mu[XX], mu[YY], mu[ZZ]);
        }
    }
}

void QgenEem::message( int len, char buf[], Resp &gr)
{
    switch (_eQGEN)
    {
        case eQGEN_OK:
            if (gr.nAtom() > 0)
            {
                gr.calcPot();
                gr.calcRms();
                gr.statistics( len, buf);
            }
            else
            {
                sprintf(buf, "Charge generation finished correctly.\n");
            }
            break;
        case eQGEN_NOTCONVERGED:
            sprintf(buf, "Charge generation did not converge.\n");
            break;
        case eQGEN_NOSUPPORT:
            sprintf(buf, "No charge generation support for (some of) the atomtypes.\n");
            break;
        case eQGEN_ERROR:
        default:
            sprintf(buf, "Unknown status %d in charge generation.\n", _eQGEN);
    }
}

void QgenEem::checkSupport(const Poldata &pd)
{
    int  i;
    bool bSupport = true;

    for (i = 0; (i < _natom); i++)
    {
        if (!pd.haveEemSupport(_iChargeDistributionModel, _elem[i].c_str(), TRUE))
        {
            fprintf(stderr, "No charge generation support for atom %s, model %s\n",
                    _elem[i].c_str(), getEemtypeName(_iChargeDistributionModel));
            bSupport = false;
        }
    }
    if (bSupport)
    {
        _eQGEN = eQGEN_OK;
    }
    else
    {
        _eQGEN = eQGEN_NOSUPPORT;
    }
}

void QgenEem::updateFromPoldata(t_atoms *atoms, const Poldata &pd)
{
    int i, j, n, nz;

    for (i = j = 0; (i < atoms->nr); i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            _chi0[j]       = pd.getChi0(_iChargeDistributionModel, _elem[j].c_str());
            _j00[j]        = pd.getJ00(_iChargeDistributionModel, _elem[j].c_str());
            nz             = pd.getNzeta(_iChargeDistributionModel, _elem[j].c_str());
            for (n = 0; (n < nz); n++)
            {
                _zeta[j][n] = pd.getZeta(_iChargeDistributionModel, _elem[j].c_str(), n);
                _q[j][n]    = pd.getQ(_iChargeDistributionModel, _elem[j].c_str(), n);
                _row[j][n]  = pd.getRow(_iChargeDistributionModel, _elem[j].c_str(), n);
            }
            j++;
        }
    }
}

int QgenEem::generateChargesSm(FILE *fp,
                                  const Poldata &pd,
                                  t_atoms *atoms,
                                  double tol, 
                                  int maxiter, 
                                  double *chieq)
{
    std::vector<double>       qq;
    int                     i, j, iter;
    double                    rms;

    checkSupport(pd);
    if (eQGEN_OK == _eQGEN)
    {

        updateFromPoldata(atoms, pd);

        qq.resize(atoms->nr+1);
        for (i = j = 0; (i < atoms->nr); i++)
        {
            if (atoms->atom[i].ptype != eptShell)
            {
                qq[j] = _q[j][_nZeta[j]-1];
                j++;
            }
        }
        iter = 0;
        calcJab();
        calcRhs();
        do
        {
            updateJ00();
            if (debug)
            {
                debugFun(debug);
            }
            solveQEem(debug, 1.0);
            rms = 0;
            for (i = j = 0; (i < atoms->nr); i++)
            {
                if (atoms->atom[i].ptype != eptShell)
                {
                    rms  += gmx::square(qq[j] - _q[j][_nZeta[j]-1]);
                    qq[j] = _q[j][_nZeta[j]-1];
                    j++;
                }
            }
            rms = sqrt(rms/atoms->nr);
            iter++;
        }
        while ((rms > tol) && (iter < maxiter));

        if (iter < maxiter)
        {
            _eQGEN = eQGEN_OK;
        }
        else
        {
            _eQGEN = eQGEN_NOTCONVERGED;
        }

        if (fp)
        {
            if (_eQGEN == eQGEN_OK)
            {
                fprintf(fp, "Converged to tolerance %g after %d iterations\n",
                        tol, iter);
            }
            else
            {
                fprintf(fp, "Did not converge within %d iterations. RMS = %g\n",
                        maxiter, rms);
            }
        }
        *chieq = _chieq;
    }

    if (eQGEN_OK == _eQGEN)
    {
        copyChargesToAtoms(atoms);
        print(fp, atoms);
    }

    return _eQGEN;
}

int QgenEem::generateChargesBultinck(FILE *fp,
                                     const Poldata &pd,
                                     t_atoms *atoms)
{
    checkSupport(pd);
    if (eQGEN_OK == _eQGEN)
    {
        updateFromPoldata(atoms, pd);

        calcJab();
        calcRhs();
        updateJ00();
        solveQEem(debug, 2.0);
        copyChargesToAtoms(atoms);
        print(fp, atoms);
    }

    return _eQGEN;
}

int QgenEem::generateCharges(FILE              *fp,
                             const std::string  molname, 
                             const Poldata     &pd,
                             t_atoms           *atoms,
                             double             tol,
                             int                maxiter)
{
    double chieq;
    /* Generate charges using empirical algorithms */
    if (fp)
    {
        fprintf(fp, "Generating charges for %s using %s algorithm\n",
                molname.c_str(), 
                getEemtypeName(_iChargeDistributionModel));
    }
    if (_iChargeDistributionModel == eqdBultinck)
    {
        (void) generateChargesBultinck(fp, pd, atoms);
    }
    else
    {
        (void) generateChargesSm(fp, pd, atoms, tol, maxiter, &chieq);
    }
    copyChargesToAtoms(atoms);
    
    return _eQGEN;
}

}
