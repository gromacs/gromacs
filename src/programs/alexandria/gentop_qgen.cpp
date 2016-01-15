/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "gentop_qgen.h"

#include <ctype.h>

#include "gromacs/fileio/confio.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/txtdump.h"
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
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/smalloc.h"

#include "gmx_resp.h"
#include "molprop.h"
#include "poldata.h"
#include "coulombintegrals/coulombintegrals.h"

namespace alexandria
{


GentopQgen::GentopQgen(Poldata * pd, t_atoms *atoms, gmx_atomprop_t aps,
                       rvec *x,
                       ChargeDistributionModel   iChargeDistributionModel,
                       ChargeGenerationAlgorithm iChargeGenerationAlgorithm,
                       real hfac, int qtotal, real epsr)
{
    _bWarned    = false;
    _bAllocSave = false;
    _natom      = 0;
    _eQGEN      = 0;
    _qtotal     = 0;
    _chieq      = 0;
    _hfac       = 0;
    _epsr       = 0;
    std::string  atp;
    gmx_bool     bSup = TRUE;
    int          i, j, k, atm, nz;


    _iChargeDistributionModel   = iChargeDistributionModel;
    _iChargeGenerationAlgorithm = iChargeGenerationAlgorithm;
    _hfac                       = hfac;
    _qtotal                     = qtotal;
    if (epsr <= 1)
    {
        epsr = 1;
    }
    _epsr   = epsr;
    for (i = j = 0; (i < atoms->nr); i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            _natom++;
        }
    }

    _chi0.resize(_natom);
    _rhs.resize(_natom+1);
    _elem.resize(_natom);
    _atomnr.resize(_natom);
    _row.resize(_natom);
    _Jab.resize(_natom+1);
    _zeta.resize(_natom);
    _j00.resize(_natom);
    _q.resize(_natom+1);

    _nZeta.resize(_natom+1);

    snew(_x, _natom);

    _bAllocSave = FALSE;


    /* Special case for chi_eq */
    _nZeta[_natom] = 1;
    _q[_natom].resize(_nZeta[_natom]);


    for (i = j = 0; (i < atoms->nr) && bSup; i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {

            _Jab[j].resize(_natom+1);
            atm = atoms->atom[i].atomnumber;
            if (atm == 0)
            {
                gmx_fatal(FARGS, "Don't know atomic number for %s %s",
                          *(atoms->resinfo[i].name),
                          *(atoms->atomname[j]));
            }
            atp = *atoms->atomtype[j];
            if (pd->haveEemSupport(_iChargeDistributionModel, atp, TRUE) == 0)
            {
                atp = gmx_atomprop_element(aps, atm);
                if (pd->haveEemSupport(_iChargeDistributionModel, atp, TRUE) == 0)
                {
                    fprintf(stderr, "No charge distribution support for atom %s (element %s), model %s\n",
                            *atoms->atomtype[j], atp.c_str(), Poldata::getEemtypeName(_iChargeDistributionModel).c_str());
                    bSup = FALSE;
                }
            }
            if (bSup)
            {
                _elem[j].assign(atp);
                _atomnr[j]      = atm;
                nz              = pd->getNzeta(_iChargeDistributionModel, atp);
                _nZeta[j]       = nz;

                _q[j].resize(nz);

                _zeta[j].resize(nz);

                _row[j].resize(nz);
                for (k = 0; (k < nz); k++)
                {
                    _q[j][k]    = pd->getQ(_iChargeDistributionModel, *atoms->atomtype[j], k);
                    _zeta[j][k] = pd->getZeta(_iChargeDistributionModel, *atoms->atomtype[j], k);
                    _row[j][k]  = pd->getRow(_iChargeDistributionModel, *atoms->atomtype[j], k);
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




GentopQgen::~GentopQgen()
{
    sfree(_x);
}


void GentopQgen::saveParams( Resp * gr)
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
            _qsave[i].resize(_nZeta[i]);
            _zetasave[i].resize(_nZeta[i]);
        }
        for (j = 0; (j < _nZeta[i]); j++)
        {
            if (NULL != gr)
            {
                _q[i][j]    = (real)gr->getQ( i, j);
                _zeta[i][j] = gr->getZeta( i, j);
            }
            _qsave[i][j]    = _q[i][j];
            _zetasave[i][j] = _zeta[i][j];
        }
    }
    _bAllocSave = TRUE;
}

void GentopQgen::getParams( Resp * gr)
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
                if (NULL != gr)
                {
                    gr->setQ( i, j, _q[i][j]);
                    gr->setZeta( i, j, _zeta[i][j]);
                }
            }
        }
    }
    else
    {
        fprintf(stderr, "WARNING: no ESP charges generated.\n");
    }
}

int GentopQgen::getNzeta( int atom)
{
    if ((0 <= atom) && (atom < _natom))
    {
        return _nZeta[atom];
    }
    return 0;

}

int GentopQgen::getRow( int atom, int z)
{
    if ((0 <= atom) && (atom < _natom) &&
        (0 <= z) && (z <= _nZeta[atom]))
    {
        return _row[atom][z];
    }
    return 0;

}
double GentopQgen::getQ(int atom, int z)
{
    if ((0 <= atom) && (atom < _natom) &&
        (0 <= z) && (z <= _nZeta[atom]))
    {
        return _q[atom][z];
    }
    return 0;

}


double GentopQgen::getZeta(int atom, int z)
{
    if ((0 <= atom) && (atom < _natom) &&
        (0 <= z) && (z <= _nZeta[atom]))
    {
        return _zeta[atom][z];
    }
    return 0;
}

real CoulombNN(real r)
{
    return 1/r;
}

real GentopQgen::calcJab(ChargeDistributionModel iChargeDistributionModel,
                         rvec xi, rvec xj,
                         int nZi, int nZj,
                         std::vector<real> zetaI, std::vector<real> zetaJ,
                         std::vector<int> rowI, std::vector<int> rowJ)
{
    int  i, j;
    rvec dx;
    real r;
    real eTot = 0;

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

void GentopQgen::solveQEem(FILE *fp,  real hardnessFactor)
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

void GentopQgen::updateJ00()
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

void GentopQgen::debugFun(FILE *fp)
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

real GentopQgen::calcSij(int i, int j)
{
    real dist, dism, Sij = 1.0;
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

void GentopQgen::calcJab()
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

void GentopQgen::calcRhs()
{
    int    i, j, k, l;
    rvec   dx;
    real   r, j1, j1q, qcore;

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


void GentopQgen::print(FILE *fp, t_atoms *atoms)
{
    int  i, j, k, m;
    rvec mu = { 0, 0, 0 };
    real qq;

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
                qq = 0;
                for (k = 0; (k < _nZeta[j]); k++)
                {
                    qq += _q[j][k];
                }

                atoms->atom[i].q = qq;
                for (m = 0; (m < DIM); m++)
                {
                    mu[m] += qq* _x[i][m] * ENM2DEBYE;
                }
                if (fp)
                {
                    fprintf(fp, "%4s %4s%5d %8g %8g",
                            *(atoms->resinfo[atoms->atom[i].resind].name),
                            *(atoms->atomname[i]), i+1, _j00[j], _chi0[j]);
                    for (k = 0; (k < _nZeta[j]); k++)
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

void GentopQgen::message( int len, char buf[], Resp * gr)
{
    switch (_eQGEN)
    {
        case eQGEN_OK:
            if (NULL != gr)
            {
                gr->calcPot();
                gr->calcRms();
                gr->statistics( len, buf);
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

void GentopQgen::checkSupport(Poldata * pd, gmx_atomprop_t aps)
{
    int      i;
    gmx_bool bSup = TRUE;

    for (i = 0; (i < _natom); i++)
    {
        if (pd->haveEemSupport(_iChargeDistributionModel, _elem[i].c_str(), TRUE) == 0)
        {
            /*sfree(elem[i]);*/
            _elem[i].assign(gmx_atomprop_element(aps, _atomnr[i]));
            if (pd->haveEemSupport(_iChargeDistributionModel, _elem[i].c_str(), TRUE) == 0)
            {
                fprintf(stderr, "No charge generation support for atom %s, model %s\n",
                        _elem[i].c_str(), Poldata::getEemtypeName(_iChargeDistributionModel).c_str());
                bSup = FALSE;
            }
        }
    }
    if (bSup)
    {
        _eQGEN = eQGEN_OK;
    }
    else
    {
        _eQGEN = eQGEN_NOSUPPORT;
    }
}

void GentopQgen::updatePd(t_atoms *atoms, Poldata * pd)
{
    int i, j, n, nz;

    for (i = j = 0; (i < atoms->nr); i++)
    {
        if (atoms->atom[i].ptype == eptAtom)
        {
            _chi0[j]       = pd->getChi0(_iChargeDistributionModel, _elem[j].c_str());
            _j00[j]        = pd->getJ00(_iChargeDistributionModel, _elem[j].c_str());
            nz             = pd->getNzeta(_iChargeDistributionModel, _elem[j].c_str());
            for (n = 0; (n < nz); n++)
            {
                _zeta[j][n] = pd->getZeta(_iChargeDistributionModel, _elem[j].c_str(), n);
                _q[j][n]    = pd->getQ(_iChargeDistributionModel, _elem[j].c_str(), n);
                _row[j][n]  = pd->getRow(_iChargeDistributionModel, _elem[j].c_str(), n);
            }
            j++;
        }
    }
}

int GentopQgen::generateChargesSm(FILE *fp,
                                  Poldata * pd,
                                  t_atoms *atoms,
                                  real tol, int maxiter, gmx_atomprop_t aps,
                                  real *chieq)
{
    std::vector<real>       qq;
    int                     i, j, iter;
    real                    rms;

    checkSupport(pd, aps);
    if (eQGEN_OK == _eQGEN)
    {

        updatePd(atoms, pd);

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
        print(fp, atoms);
    }

    return _eQGEN;
}

int GentopQgen::generateChargesBultinck(FILE *fp,
                                        Poldata * pd, t_atoms *atoms,
                                        gmx_atomprop_t aps)
{
    checkSupport(pd, aps);
    if (eQGEN_OK == _eQGEN)
    {
        updatePd(atoms, pd);

        calcJab();
        calcRhs();
        updateJ00();
        solveQEem(debug, 2.0);

        print(fp, atoms);
    }

    return _eQGEN;
}

int GentopQgen::generateCharges(FILE *fp,
                                Resp * gr,
                                const std::string molname, Poldata * pd,
                                t_atoms *atoms,
                                real tol, int maxiter, int maxcycle,
                                gmx_atomprop_t aps)
{
    int  cc, eQGEN_min = eQGEN_NOTCONVERGED;
    real chieq, chi2, chi2min = GMX_REAL_MAX;

    /* Generate charges */
    switch (_iChargeGenerationAlgorithm)
    {
        case eqgRESP:
            if (NULL == gr)
            {
                gmx_incons("No RESP data structure");
            }
            if (fp)
            {
                fprintf(fp, "Generating %s charges for %s using RESP algorithm\n",
                        Poldata::getEemtypeName(_iChargeDistributionModel).c_str(), molname.c_str());
            }
            for (cc = 0; (cc < maxcycle); cc++)
            {
                if (fp)
                {
                    fprintf(fp, "Cycle %d/%d\n", cc+1, maxcycle);
                }
                /* Fit charges to electrostatic potential */
                _eQGEN = gr->optimizeCharges(fp, maxiter, tol, &chi2);
                if (_eQGEN == eQGEN_OK)
                {
                    eQGEN_min = _eQGEN;
                    if (chi2 <= chi2min)
                    {
                        saveParams(gr);
                        chi2min = chi2;
                    }

                    if (NULL != fp)
                    {
                        fprintf(fp, "chi2 = %g kJ/mol e\n", chi2);
                    }
                    print(fp, atoms);
                }
            }
            if (maxcycle > 1)
            {
                if (fp)
                {
                    fprintf(fp, "---------------------------------\nchi2 at minimum is %g\n", chi2min);
                }
                getParams(gr);
                print(fp, atoms);
            }
            _eQGEN = eQGEN_min;
            break;
        default:
            /* Use empirical algorithms */
            if (fp)
            {
                fprintf(fp, "Generating charges for %s using %s algorithm\n",
                        molname.c_str(), Poldata::getEemtypeName(_iChargeDistributionModel).c_str());
            }
            if (_iChargeDistributionModel == eqdBultinck)
            {
                (void) generateChargesBultinck(fp, pd, atoms, aps);
            }
            else
            {
                (void) generateChargesSm(fp, pd, atoms, tol, maxiter, aps, &chieq);
            }
            saveParams(gr);
    }
    return _eQGEN;
}

}
