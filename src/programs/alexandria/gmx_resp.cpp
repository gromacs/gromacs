/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include "gmx_resp.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>

#include "gromacs/commandline/filenm.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/pdbio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/toputil.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/random/random.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/atoms.h"
#include "gromacs/topology/symtab.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"
#include "gromacs/utility/textreader.h"
#include "gromacs/utility/txtdump.h"

#include "gentop_qgen.h"
#include "nmsimplex.h"
#include "poldata.h"
#include "stringutil.h"
#include "coulombintegrals/coulombintegrals.h"

namespace alexandria
{
Resp::Resp() 
{
    rnd_                = nullptr; 
    setOptions(eqdAXp, 0, false, 5, 100, -1, false,
               0, -2, 2, true, 0);
    _bAXpRESP           = false;
    _qfac               = 1e-3;
    _bHyper             = 0.1;
    _wtot               = 0;
    _pfac               = 1;
    _bEntropy           = false;
    _rDecrZeta          = true;
}

void Resp::setOptions(ChargeDistributionModel c,
                      unsigned int            seed,
                      bool                    fitZeta, 
                      real                    zetaMin,
                      real                    zetaMax,
                      real                    deltaZeta,
                      bool                    randomZeta,
                      real                    qtot,
                      real                    qmin,
                      real                    qmax,
                      bool                    randomQ,
                      real                    watoms)
{
    _iDistributionModel = c;
    _bFitZeta           = fitZeta;
    _zmin               = zetaMin;
    _zmax               = zetaMax;
    _deltaZ             = deltaZeta;
    _bRandZeta          = randomZeta;
    _rDecrZeta          = true;
    
    _qtot               = qtot;
    _qmin               = qmin;
    _qmax               = qmax; /* e */
    _bRandQ             = randomQ;
    _watoms             = watoms;
    if (seed <= 0)
    {
        seed = gmx_rng_make_seed();
    }
    if (nullptr != rnd_)
    {
        gmx_rng_destroy(rnd_);
    }   
    rnd_ = gmx_rng_init(seed);
}

Resp::~Resp()
{
    gmx_rng_destroy(rnd_);
}

void Resp::setAtomInfo(t_atoms                   *atoms, 
                       const alexandria::Poldata &pd,
                       const rvec                 x[])
{
    for (int i = 0; (i < atoms->nr); i++)
    {
        ra_.push_back(RespAtom(atoms->atom[i].atomnumber,
                               atoms->atom[i].type, 
                               0,
                               x[i]));
        if (findRAT(atoms->atom[i].type) == endRAT())
        {
            ratype_.push_back(RespAtomType(atoms->atom[i].type,
                                           *(atoms->atomtype[i]), pd, 
                                           _iDistributionModel, _dzatoms));
        }
    }
}

void Resp::summary(FILE             *fp,
                   std::vector<int> &symmetricAtoms)
{
    if (NULL != fp)
    {
        fprintf(fp, "There are %d atoms, %d atomtypes %d parameters for (R)ESP fitting.\n",
                static_cast<int>(nAtom()), 
                static_cast<int>(nAtomType()),
                static_cast<int>(nParam()));
        for (size_t i = 0; (i < nAtom()); i++)
        {
            fprintf(fp, " %d", symmetricAtoms[i]);
        }
        fprintf(fp, "\n");
    }
}

void Resp::addParam(int aindex, eParm eparm, size_t zz)
{
    if (eparm == eparmQ)
    {
        range_check(aindex, 0, nAtom());
        raparam_.push_back(RespParam(eparm, aindex, zz));
        if (debug)
        {
            fprintf(debug, "GRESP: Adding parameter %d for atom %d\n", 
                    eparm, aindex);
        }
    }
    else if (_bFitZeta && (zz < ratype_[aindex].getNZeta()))
    {
        range_check(aindex, 0, nAtomType());
        raparam_.push_back(RespParam(eparm, aindex, zz));
        FILE * debug = stdout;
        if (debug)
        {
            fprintf(debug, "GRESP: Adding parameter %d for atom type %d zz %d\n", 
                    eparm, aindex, static_cast<int>(zz));
        }
    }
}

void Resp::setAtomSymmetry(const std::vector<int> &symmetricAtoms)
{
    GMX_RELEASE_ASSERT(!ra_.empty(), "RespAtom vector not initialized");
    GMX_RELEASE_ASSERT(!ratype_.empty(), "RespAtomType vector not initialized");
    GMX_RELEASE_ASSERT(nParam() == 0, "There are parameters already in the Resp structure");

    /* Map the symmetric atoms */
    for (size_t i = 0; (i < nAtom()); i++)
    {
        int atype = ra_[i].atype();
        RespAtomTypeIterator rai = findRAT(atype);
        if (0 == i)
        {
            /* The first charge is not a free variable, it follows from the 
             * total charge. Only add the zeta values here.
             */
            for (size_t zz = 0; (zz < rai->getNZeta()); zz++)
            {
                addParam(atype, eparmZ, zz);
                (rai->beginRZ()+zz)->setZindex(nParam()-1);
            }
        }
        else if (symmetricAtoms[i] == static_cast<int>(i))
        {
            // We optimize at most 1 charge per atom, so use index 0
            addParam(i, eparmQ, 0);
            ra_[i].setQindex(nParam()-1);
            // Check if we have this atype covered already
            if (rai == ratype_.end())
            {
                // New atom type
                for (size_t zz = 0; (zz < rai->getNZeta()); zz++)
                {
                    addParam(atype, eparmZ, zz);
                }
            }
        }
        else if (symmetricAtoms[i] > static_cast<int>(i))
        {
            gmx_fatal(FARGS, "The symmetricAtoms array can not point to larger atom numbers");
        }
        else
        {
            // Symmetric atom 
            ra_[i].setQindex(ra_[symmetricAtoms[i]].qIndex());
            if (debug)
            {
                fprintf(debug, "Atom %d is a copy of atom %d\n",
                        static_cast<int>(i+1), symmetricAtoms[i]+1);
            }
        }
    }
    if (debug)
    {
        size_t maxz = 0;
        for (size_t i = 0; (i < nAtomType()); i++)
        {
            if (ratype_[i].getNZeta() > maxz)
            {
                maxz = ratype_[i].getNZeta();
            }
        }

        fprintf(debug, "GRQ: %3s %5s", "nr", "type");
        fprintf(debug, " %8s %8s\n", "q", "zeta");
        for (size_t i = 0; (i < nAtom()); i++)
        {
            int atype = ra_[i].atype();
            RespAtomTypeIterator rai = findRAT(atype);
            fprintf(debug, "GRQ: %3d %5s", static_cast<int>(i+1),
                    rai->getAtomtype().c_str());
            for (RespZetaConstIterator zz = rai->beginRZ(); zz <  rai->endRZ(); ++zz)
            {
                fprintf(debug, " %8.4f %8.4f\n",
                        ra_[i].charge(), zz->zeta());
            }
            fprintf(debug, "\n");
        }
    }
    printf("There are %d variables to optimize for %d atoms and %d ESP points.\n",
           static_cast<int>(nParam()), static_cast<int>(nAtom()),
           static_cast<int>(nEsp()));
}

void Resp::writeHisto(const std::string &fn,
                      const std::string &title,
                      const gmx_output_env_t *oenv)
{
    FILE       *fp;
    gmx_stats_t gs;
    real       *x, *y;
    int         nbin = 100;

    if (0 == fn.size())
    {
        return;
    }
    gs = gmx_stats_init();
    for (size_t i = 0; (i < nEsp()); i++)
    {
        gmx_stats_add_point(gs, i, gmx2convert(_potCalc[i], eg2cHartree_e), 0, 0);
    }

    gmx_stats_make_histogram(gs, 0, &nbin, ehistoY, 1, &x, &y);

    fp = xvgropen(fn.c_str(), title.c_str(), "Pot (1/a.u.)", "()", oenv);
    for (int i = 0; (i < nbin); i++)
    {
        fprintf(fp, "%10g  %10g\n", x[i], y[i]);
    }
    sfree(x);
    sfree(y);
    fclose(fp);
    gmx_stats_free(gs);
}

void Resp::writeDiffCube(Resp              &src, 
                         const std::string &cubeFn,
                         const std::string &histFn, 
                         const std::string &title, 
                         const gmx_output_env_t *oenv,
                         int rho)
{
    FILE       *fp;
    int         i, m, ix, iy, iz;
    real        pp, q, r, rmin;
    gmx_stats_t gst = NULL, ppcorr = NULL;

    if (0 != histFn.size())
    {
        gst    = gmx_stats_init();
        ppcorr = gmx_stats_init();
    }
    if (0 != cubeFn.size())
    {
        fp = gmx_ffopen(cubeFn.c_str(), "w");
        fprintf(fp, "%s\n", title.c_str());
        fprintf(fp, "POTENTIAL\n");
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n",
                static_cast<int>(nAtom()),
                gmx2convert(_origin[XX], eg2cBohr),
                gmx2convert(_origin[YY], eg2cBohr),
                gmx2convert(_origin[ZZ], eg2cBohr));
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", _nxyz[XX],
                gmx2convert(_space[XX], eg2cBohr), 0.0, 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", _nxyz[YY],
                0.0, gmx2convert(_space[YY], eg2cBohr), 0.0);
        fprintf(fp, "%5d%12.6f%12.6f%12.6f\n", _nxyz[ZZ],
                0.0, 0.0, gmx2convert(_space[ZZ], eg2cBohr));

        for (size_t m = 0; (m < nAtom()); m++)
        {
            q = ra_[m].charge();
            fprintf(fp, "%5d%12.6f%12.6f%12.6f%12.6f\n",
                    ra_[m].atomnumber(), q,
                    gmx2convert(ra_[m].x()[XX], eg2cBohr),
                    gmx2convert(ra_[m].x()[YY], eg2cBohr),
                    gmx2convert(ra_[m].x()[ZZ], eg2cBohr));
        }

        for (ix = m = 0; ix < _nxyz[XX]; ix++)
        {
            for (iy = 0; iy < _nxyz[YY]; iy++)
            {
                for (iz = 0; iz < _nxyz[ZZ]; iz++, m++)
                {
                    if (src.nEsp() > 0)
                    {
                        pp = _potCalc[m] - src._pot[m];
                        if (NULL != ppcorr)
                        {
                            gmx_stats_add_point(ppcorr,
                                                gmx2convert(src._pot[m], eg2cHartree_e),
                                                gmx2convert(_potCalc[m], eg2cHartree_e), 0, 0);
                        }
                    }
                    else
                    {
                        if (rho == 0)
                        {
                            pp = gmx2convert(_potCalc[m], eg2cHartree_e);
                        }
                        else
                        {
                            pp = _rho[m]*pow(BOHR2NM, 3);
                        }
                    }
                    fprintf(fp, "%13.5e", pp);
                    if (iz % 6 == 5)
                    {
                        fprintf(fp, "\n");
                    }
                    if (NULL != gst)
                    {
                        rmin = 1000;
                        /* Add point to histogram! */
                        for (RespAtomIterator i = ra_.begin(); i < ra_.end(); ++i)
                        {
                            gmx::RVec dx;
                            rvec_sub(i->x(), _esp[m], dx);
                            r = norm(dx);
                            if (r < rmin)
                            {
                                rmin = r;
                            }
                        }
                        gmx_stats_add_point(gst, rmin, pp, 0, 0);
                    }
                }
                if ((iz % 6) != 0)
                {
                    fprintf(fp, "\n");
                }
            }
        }
        fclose(fp);
    }
    if (NULL != gst)
    {
        int   nb = 0;
        real *x  = NULL, *y = NULL;

        fp = xvgropen(histFn.c_str(), "Absolute deviation from QM", "Distance (nm)",
                      "Potential", oenv);
        gmx_stats_dump_xy(gst, fp);
        if (0)
        {
            gmx_stats_make_histogram(gst, 0.01, &nb, ehistoX, 0, &x, &y);
            gmx_stats_free(gst);
            for (i = 0; (i < nb); i++)
            {
                fprintf(fp, "%10g  %10g\n", x[i], y[i]);
            }
            sfree(x);
            sfree(y);
        }
        fclose(fp);
        fp = xvgropen("diff-pot.xvg", "Correlation between QM and Calc", "Pot (QM)",
                      "Pot (Calc)", oenv);
        gmx_stats_dump_xy(ppcorr, fp);
        fclose(fp);
    }
}

void Resp::writeCube(const std::string &fn, const std::string &title)
{
    Resp dummy;
    writeDiffCube(dummy,  fn, NULL, title, NULL, 0);
}

void Resp::writeRho(const std::string &fn, const std::string &title)
{
    Resp dummy;
    writeDiffCube(dummy,  fn, NULL, title, NULL, 1);
}

void Resp::readCube(const std::string &fn, bool bESPonly)
{
    int    natom, nxyz[DIM] = { 0, 0, 0 };
    double space[DIM] = { 0, 0, 0 };
    
    gmx::TextReader tr(fn);
    std::string     tmp;
    int             line = 0;
    bool            bOK  = true;
    while(bOK && tr.readLine(&tmp)) 
    {
        while (!tmp.empty() && tmp[tmp.length()-1] == '\n') 
        {
            tmp.erase(tmp.length()-1);
        }
        if (0 == line)
        {
            printf("%s\n", tmp.c_str());
        }
        else if (1 == line && tmp.compare("POTENTIAL") != 0)
        {
            bOK = false;
        }
        else if (2 == line)
        {
            double origin[DIM];
            bOK = (4 == sscanf(tmp.c_str(), "%d%lf%lf%lf",
                               &natom, &origin[XX], &origin[YY], &origin[ZZ]));
            if (bOK && !bESPonly)
            {
                _origin[XX] = origin[XX];
                _origin[YY] = origin[YY];
                _origin[ZZ] = origin[ZZ];
            }
        }
        else if (3 == line)
        {
            bOK = (2 == sscanf(tmp.c_str(), "%d%lf",
                               &nxyz[XX], &space[XX]));
        }
        else if (4 == line)
        {
            bOK = (2 == sscanf(tmp.c_str(), "%d%*s%lf",
                               &nxyz[YY], &space[YY]));
        }
        else if (5 == line)
        {
            bOK = (2 == sscanf(tmp.c_str(), "%d%*s%*s%lf",
                               &nxyz[ZZ], &space[ZZ]));
            if (bOK)
            {
                for (int m = 0; (m < DIM); m++)
                {
                    _nxyz[m]  = nxyz[m];
                    _space[m] = space[m];
                }
                for (int m = 0; (m < DIM); m++)
                {
                    _origin[m] = convert2gmx(_origin[m], eg2cBohr);
                    _space[m]  = convert2gmx(_space[m], eg2cBohr);
                }
            }
            _pot.clear();
        }
        else if (line >= 6 && line < 6+natom)
        {
            double lx, ly, lz, qq;
            int    anr, m = line - 6;
            bOK = (5 == sscanf(tmp.c_str(), "%d%lf%lf%lf%lf",
                               &anr, &qq, &lx, &ly, &lz));
            if (bOK)
            {
                if (!bESPonly)
                {
                    ra_[m].setAtomnumber(anr);
                    ra_[m].setCharge(qq);
                }
                RVec xx;
                xx[XX] = convert2gmx(lx, eg2cBohr);
                xx[YY] = convert2gmx(ly, eg2cBohr);
                xx[ZZ] = convert2gmx(lz, eg2cBohr);
                ra_[m].setX(xx);
            }
        }
        else if (line >= 6+natom)
        {
            std::vector<std::string> ss = gmx::splitString(tmp);
            for(const auto &s : ss)
            {
                _pot.push_back(convert2gmx(atof(s.c_str()), eg2cHartree_e));
            }
        }
        
        line++;
    }
    if (bOK)
    {
        _esp.clear();
        for(int ix = 0; ix < _nxyz[XX]; ix++)
        {
            for(int iy = 0; iy < _nxyz[YY]; iy++)
            {
                for(int iz = 0; iz < _nxyz[ZZ]; iz++)
                {
                    rvec e;
                    e[XX] = _origin[XX] + ix*_space[XX];
                    e[YY] = _origin[YY] + iy*_space[YY];
                    e[ZZ] = _origin[ZZ] + iz*_space[ZZ];
                    
                    _esp.push_back(e);
                }
            }
        }
    }
    bOK = bOK && (_esp.size() == _pot.size());
    if (!bOK)
    {
        gmx_fatal(FARGS, "Error reading %s. Found %d potential values, %d coordinates and %d atoms",
                  fn.c_str(), static_cast<int>(_pot.size()), static_cast<int>(_esp.size()),
                  static_cast<int>(ra_.size()));
    }
}

void Resp::copyGrid(Resp &src)
{
    int m;

    for (m = 0; (m < DIM); m++)
    {
        _origin[m] = src._origin[m];
        _space[m]  = src._space[m];
        _nxyz[m]   = src._nxyz[m];
    }
    int nesp = src.nEsp();
    _esp.clear();
    _pot.resize(nesp, 0);
    _potCalc.resize(nesp, 0);
    for (m = 0; (m < nesp); m++)
    {
        _esp.push_back(src._esp[m]);
    }
}

//Resp * Resp::copy()
//{
//   Resp * dest = new Resp(_iDistributionModel, _qtot);

//  memcpy(dest, this, sizeof(*this));

//    return dest;
//}

void Resp::makeGrid(real spacing, matrix box, rvec x[])
{
    if (0 != nEsp())
    {
        fprintf(stderr, "Overwriting existing ESP grid\n");
    }
    if (spacing <= 0)
    {
        spacing = 0.1;
        fprintf(stderr, "spacing too small, setting it to %g\n", spacing);
    }
    for (size_t i = 0; (i < nAtom()); i++)
    {
        ra_[i].setX(x[i]);
    }
    for (int m = 0; (m < DIM); m++)
    {
        _nxyz[m]  = 1+(int) (box[m][m]/spacing);
        _space[m] = box[m][m]/_nxyz[m];
    }
    _esp.clear();
    _potCalc.clear();
    for (int i = 0; (i < _nxyz[XX]); i++)
    {
        rvec xyz;
        xyz[XX] = (i-0.5*_nxyz[XX])*_space[XX];
        for (int j = 0; (j < _nxyz[YY]); j++)
        {
            xyz[YY] = (j-0.5*_nxyz[YY])*_space[YY];
            for (int k = 0; (k < _nxyz[ZZ]); k++)
            {
                xyz[ZZ] = (k-0.5*_nxyz[ZZ])*_space[ZZ];
                _esp.push_back(xyz);
                _potCalc.push_back(0);
            }
        }
    }
}

void Resp::calcRho()
{
    double pi32 = pow(M_PI, -1.5);
    if (_rho.size() < nEsp())
    {
        _rho.resize(nEsp(), 0);
    }
    for (size_t i = 0; (i < _rho.size()); i++)
    {
        double V = 0;
        for (const auto &ra : ra_)
        {
            double vv = 0;
            gmx::RVec dx;
            rvec_sub(_esp[i], ra.x(), dx);
            double r = norm(dx);
            int atype = ra.atype();
            RespAtomTypeIterator rat = findRAT(atype);
            GMX_RELEASE_ASSERT(rat == endRAT(), "Can not find atomtype");
            switch (_iDistributionModel)
            {
                case eqdYang:
                case eqdRappe:
                    vv = ra.charge()*Nuclear_SS(r, 
                                                rat->beginRZ()->row(),
                                                rat->beginRZ()->zeta());
                    break;
                case eqdAXg:
                    vv = 0;
                    for (auto k = rat->beginRZ(); k < rat->endRZ(); ++k)
                    {
                        real z = k->zeta();
                        real q = k->q();
                        if (k == rat->endRZ()-1)
                        {
                            q = ra.charge();
                        }
                        if (z > 0 && q != 0)
                        {
                            vv -= (q*pi32*exp(-gmx::square(r*z))*
                                   pow(z, 3));
                        }
                    }
                    break;
                case eqdBultinck:
                case eqdAXp:
                case eqdAXs:
                default:
                    gmx_fatal(FARGS, "Krijg nou wat, iDistributionModel = %d!",
                              _iDistributionModel);
            }
            V  += vv;
        }
        _rho[i] = V;
    }
}

void Resp::calcPot()
{
    std::fill(_potCalc.begin(), _potCalc.end(), 0);
    for (size_t i = 0; (i < nEsp()); i++)
    {
        double V = 0;
        for (auto ra : ra_)
        {
            gmx::RVec            dx;
            rvec_sub(_esp[i], ra.x(), dx);
            double               r     = norm(dx);
            int                  atype = ra.atype();
            RespAtomTypeIterator rat   = findRAT(atype);
            double vv                  = 0;
            switch (_iDistributionModel)
            {
                case eqdBultinck:
                case eqdAXp:
                    if (r > 0.01)
                    {
                        vv = ra.charge()/r;
                    }
                    break;
                case eqdAXs:
                    vv = 0;
                    for(auto k = rat->beginRZ(); k < rat->endRZ(); ++k)
                    {
                        real q = k->q();
                        if (k == rat->endRZ()-1)
                        {
                            q = ra.charge();
                        }
                        vv += q*Nuclear_SS(r, 
                                           k->row(),
                                           k->zeta());
                    }
                    break;
                case eqdYang:
                case eqdRappe:
                    vv = ra.charge()*Nuclear_SS(r, 
                                                rat->beginRZ()->row(),
                                                rat->beginRZ()->zeta());
                    break;
                case eqdAXg:
                    vv = 0;
                    for(auto k = rat->beginRZ(); k < rat->endRZ(); ++k)
                    {
                        real q = k->q();
                        if (k == rat->endRZ()-1)
                        {
                            q = ra.charge();
                        }
                        vv += q*Nuclear_GG(r, k->zeta());
                    }
                    break;
                default:
                    gmx_fatal(FARGS, "Krijg nou wat, iDistributionModel = %s!",
                              getEemtypeName(_iDistributionModel));
            }
            V  += vv;
        }
        _potCalc[i] = V*ONE_4PI_EPS0;
    }
}

void Resp::warning(const std::string fn, int line)
{
    fprintf(stderr, "WARNING: It seems like you have two sets of ESP data in your file\n         %s\n", fn.c_str());
    fprintf(stderr, "         using the second set, starting at line %d\n", line);
}

void Resp::setVector(double *params)
{
    size_t n = 0;

    for (const auto &rp : raparam_)
    {
        if (rp.eParam() == eparmQ)
        {
            /* First do charges */
            int atom = rp.aIndex();
            if (_bRandQ)
            {
                params[n] = 0.2*(gmx_rng_uniform_real(rnd_)-0.5);
            }
            else
            {
                params[n] = ra_[atom].charge();
            }
        }
        else 
        {
            int atype = rp.aIndex();
            int zz    = rp.zIndex();
            RespAtomTypeIterator rai = findRAT(atype);
            RespZetaIterator rzi = rai->beginRZ() + zz;
            if (_bRandZeta)
            {
                real zmin = _zmin;
                real zmax = _zmax;
                
                if ((_deltaZ > 0) && (rai->getBRestrained()))
                {
                    zmin = rzi->zetaRef()-_deltaZ;
                    zmax = rzi->zetaRef()+_deltaZ;
                }
                /*if ((zz > 1) && (_rDecrZeta >= 0))
                    {
                    zmax = ra_[i].getZeta(zz-1)-_rDecrZeta;
                    if (zmax < zmin)
                    {
                    zmax = (zmin+ra_[i].getZeta(zz-1))/2;
                    }
                    }*/
                params[n] = zmin + (zmax-zmin)*gmx_rng_uniform_real(rnd_);
            }
            else
            {
                params[n] = rzi->zeta();
            }
        }
        n++;
    }
}

void Resp::getVector(double *params)
{
    double qtot = 0;
    for(auto &ra : ra_)
    {
        /* First do charges */
        int qi = ra.qIndex();
        if (qi >= 0)
        {
            ra.setCharge(params[qi]);
            qtot += params[qi];
        }
    }
    ra_[0].setCharge(_qtot-qtot);
    
    for(auto &rat : ratype_)
    {
        for(auto rz = rat.beginRZ(); rz < rat.endRZ(); ++rz)
        {
            int zi = rz->zIndex();
            if (zi >= 0)
            {
                rz->setZeta(params[zi]);
            }
        }
    }
}

void Resp::addEspPoint(double x, double y,
                       double z, double V)
{
    rvec e;
    e[XX] = x;
    e[YY] = y;
    e[ZZ] = z;
    
    _esp.push_back(e);
    _pot.push_back(V);
    _potCalc.push_back(0);
}

real Resp::myWeight(size_t iatom) const
{
    if (iatom < nAtom())
    {
        return _watoms;
    }
    else
    {
        return 1.0;
    }
}

void Resp::potLsq( gmx_stats_t lsq)
{
    for(size_t i = 0; (i < nEsp()); i++)
    {
        double w = myWeight( i);
        if (w > 0)
        {
            gmx_stats_add_point(lsq,
                                gmx2convert(_pot[i], eg2cHartree_e),
                                gmx2convert(_potCalc[i], eg2cHartree_e), 0, 0);
        }
    }
}

void Resp::calcRms()
{
    double pot2, s2, sum2, w, wtot, entropy;
    char   buf[STRLEN];

    pot2 = sum2 = wtot = entropy = 0;
    sprintf(buf, " - weight %g in fit", _watoms);
    for (size_t i = 0; (i < nEsp()); i++)
    {
        w = myWeight(i);
        if ((NULL != debug) && (i < 4*nAtom()))
        {
            fprintf(debug, "ESP %d QM: %g EEM: %g DIFF: %g%s\n",
                    static_cast<int>(i), _pot[i], _potCalc[i],
                    _pot[i]-_potCalc[i],
                    (i < nAtom())  ? buf : "");
        }
        s2    = w*gmx::square(_pot[i]-_potCalc[i]);
        if ((s2 > 0) && (_bEntropy))
        {
            entropy += s2*log(s2);
        }
        sum2 += s2;
        pot2 += w*gmx::square(_pot[i]);
        wtot += w;
    }
    _wtot = wtot;
    if (wtot > 0)
    {
        _rms     = gmx2convert(sqrt(sum2/wtot), eg2cHartree_e);
        _entropy = gmx2convert(entropy/wtot, eg2cHartree_e);
    }
    else
    {
        _rms     = 0;
        _entropy = 0;
    }
    _rrms = sqrt(sum2/pot2);
}

real Resp::getRms(real *wtot, real *rrms)
{
    calcRms();
    *wtot = _wtot;
    *rrms = _rrms;
    if (_bEntropy)
    {
        return _entropy;
    }
    else
    {
        return _rms;
    }
}

double Resp::calcPenalty()
{
    double p, b2;

    p = 0;
    /* Check for excessive charges */
    for (auto &ra : ra_)
    {
        real qi    = ra.charge();
        int  atype = ra.atype();
        RespAtomTypeIterator rat = findRAT(atype);
        for (auto z = rat->beginRZ(); z < rat->endRZ(); ++z)
        {
            qi += z->q();
        }
        if (qi < _qmin)
        {
            p += gmx::square(_qmin-qi);
        }
        else if (qi > _qmax)
        {
            p += gmx::square(_qmax-qi);
        }
        else if ((qi < -0.02) && (ra.atomnumber() == 1))
        {
            p += qi*qi;
        }
    }
    p *= _pfac;
    if (_bAXpRESP && (_iDistributionModel == eqdAXp))
    {
        b2 = gmx::square(_bHyper);
        for (size_t i = 0; (i < nAtom()); i++)
        {
            p += sqrt(gmx::square(ra_[i].charge()) + b2) - _bHyper;
        }
        p = (_qfac * p);
    }
    _penalty = p;
    
    return _penalty;
}

// Writen in C style, needed as an function argument
double chargeFunction(void *gr, double v[])
{
    Resp *resp = (Resp *)gr;
    real  rrms, rms  = 0;
    real  wtot;

    resp->getVector(v);
    resp->calcPot();
    double penalty = resp->calcPenalty();
    rms = resp->getRms(&wtot, &rrms);
    
    return rms + penalty;
}

void Resp::statistics( int len, char buf[])
{
    if (len >= 100)
    {
        sprintf(buf, "RMS: %10e [Hartree/e] RRMS: %10e Entropy: %10e Penalty: %10e",
                _rms, _rrms, _entropy, _penalty);
    }
    else
    {
        fprintf(stderr, "buflen too small (%d) in gmx_resp_statistics\n", len);
    }
}

int Resp::optimizeCharges(FILE *fp,  int maxiter,
                          real toler, real *rms)
{
    std::vector<double> param;
    double  ccc;
    int     bConv;
    char    buf[STRLEN];

    param.resize(nParam(), 0);
    setVector(param.data());

    bConv = nmsimplex(fp, (void *)this, chargeFunction, 
                      param.data(), nParam(),
                      toler, 1, maxiter, &ccc);
    if (bConv)
    {
        statistics(STRLEN-1, buf);
    }
    else
    {
        printf("NM Simplex did not converge\n\n");
    }

    if (_bEntropy)
    {
        *rms = _entropy;
    }
    else
    {
        *rms = _rms;
    }

    getVector(param.data());

    if (bConv)
    {
        return eQGEN_OK;
    }
    else
    {
        return eQGEN_NOTCONVERGED;
    }
}


void Resp::potcomp(const std::string &potcomp,
                   const std::string &pdbdiff,
                   const gmx_output_env_t *oenv)
{
    double  pp, exp, eem;
    FILE   *fp;
    int     unit = eg2cHartree_e;

    if (0 != potcomp.size())
    {
        const char *pcleg[2] = { "Atoms", "ESP points" };
        fp = xvgropen(potcomp.c_str(), "Electrostatic potential", unit2string(unit), unit2string(unit), oenv);
        xvgr_legend(fp, 2, pcleg, oenv);
        fprintf(fp, "@type xy\n");
        for (size_t i = 0; (i < nEsp()); i++)
        {
            /* Conversion may or may not be in vain depending on unit */
            exp = gmx2convert(_pot[i], unit);
            eem = gmx2convert(_potCalc[i], unit);
            if (i == nAtom())
            {
                fprintf(fp, "&\n");
                fprintf(fp, "@type xy\n");
            }
            fprintf(fp, "%10.5e  %10.5e\n", exp, eem);
        }
        fprintf(fp, "&\n");
        fclose(fp);
    }
    if (0 != pdbdiff.c_str())
    {
        fp = fopen(pdbdiff.c_str(), "w");
        fprintf(fp, "REMARK All distances are scaled by a factor of two.\n");
        for (size_t i = 0; (i < nEsp()); i++)
        {
            exp = gmx2convert(_pot[i], eg2cHartree_e);
            eem = gmx2convert(_potCalc[i], eg2cHartree_e);
            pp  = _pot[i]-_potCalc[i];
            fprintf(fp, "%-6s%5u  %-4.4s%3.3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f\n",
                    "ATOM", 1, "HE", "HE", ' ', static_cast<int>(i+1),
                    ' ', 20*_esp[i][XX], 20*_esp[i][YY], 20*_esp[i][ZZ], 0.0, pp);
        }
        fclose(fp);
    }
}

double Resp::getAtomCharge(int atom) const
{
    range_check(atom, 0, nAtom());
    double                    q     = ra_[atom].charge();
    int                       atype = ra_[atom].atype();
    RespAtomTypeConstIterator rat   = findRAT(atype);
    for (auto z = rat->beginRZ(); z < rat->endRZ()-1; ++z)
    {
        q += z->q();
    }
    return q;
}

double Resp::getCharge(int atom, size_t zz) const
{
    range_check(atom, 0, nAtom());
    double                    q     = ra_[atom].charge();
    int                       atype = ra_[atom].atype();
    RespAtomTypeConstIterator rat   = findRAT(atype);
    if (zz < rat->getNZeta())
    {
        q = (rat->beginRZ()+zz)->q();
    }
    return q;
}

double Resp::getZeta(int atom, int zz) const
{
    range_check(atom, 0, nAtom());
    int atype = ra_[atom].atype();
    RespAtomTypeConstIterator rat = findRAT(atype);
    range_check(zz, 0, rat->getNZeta());
    
    return (rat->beginRZ()+zz)->zeta();
}

void Resp::setCharge(int atom, int zz, double q)
{
    range_check(atom, 0, nAtom());
    int atype = ra_[atom].atype();
    RespAtomTypeIterator rat = findRAT(atype);
    range_check(zz, 0, rat->getNZeta());
    (rat->beginRZ()+zz)->setQ(q);
}

void Resp::setZeta(int atom, int zz, double zeta)
{
    range_check(atom, 0, nAtom());
    int atype = ra_[atom].atype();
    RespAtomTypeIterator rat = findRAT(atype);
    range_check(zz, 0, rat->getNZeta());
    (rat->beginRZ()+zz)->setZeta(zeta);
}

} // namespace
