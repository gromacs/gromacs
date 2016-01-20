/*! \internal \brief
 * Implements part of the alexandria program.
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */
#include "gmxpre.h"

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/gmxlib/nrnb.h"
#include "gromacs/gmxpreprocess/convparm.h"
#include "gromacs/gmxpreprocess/gpp_atomtype.h"
#include "gromacs/gmxpreprocess/grompp.h"
#include "gromacs/gmxpreprocess/pdb2top.h"
#include "gromacs/linearalgebra/matrix.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdlib/force.h"
#include "gromacs/mdlib/mdatoms.h"
#include "gromacs/mdlib/shellfc.h"
#include "gromacs/mdlib/vsite.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/state.h"
#include "gromacs/random/random.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/timing/wallcycle.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/init.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/stringutil.h"

#include "nmsimplex.h"

// alexandria stuff
#include "gentop_core.h"
#include "gentop_qgen.h"
#include "gmx_simple_comm.h"
#include "moldip.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "mymol.h"
#include "poldata.h"
#include "poldata-low.h"
#include "poldata_xml.h"
#include "stringutil.h"

typedef struct OptMask
{
    int         ngtb;
    std::string cgt; 
};
typedef std::vector<OptMask>::iterator OptMaskIterator;

class ForceConstants
{
private:
    std::vector<opt_mask> om_[ebtsNR];
public:
    ForceConstants() {}

    void addForceConstant(int bt, OptMask om) { om_[bt].push_back(om); }

    void analyzeIdef(FILE                      *fp, 
                     t_idef                    *idef,
                     const alexandria::Mymol   &mymol,
                     const alexandria::Poldata &pd);

    OptMaskIterator beginOM(int bt) { return om_[bt].begin(); }

    OptMaskIterator endOM(int bt) { return om_[bt].end(); }
};


void ForceConstants::analyzeIdef(FILE                          *fp,
                                 std::vector<alexandria::MyMol> mm,
                                 const Poldata                  pd,
                                 bool                           bOpt[])
{
    int          gt, i, ak, al, ft;
    int          ntot[ebtsNR];
    std::string  aai, aaj, aak, aal, params;
    const char  *btsnames[ebtsNR] =  { "bond", "angle", "proper", "improper", NULL, NULL };

    for (int bt = 0; (bt <= ebtsIDIHS); bt++)
    {
        if (!bOpt[bt])
        {
            continue;
        }
        switch (bt)
        {
        case ebtsBONDS:
            omt->nb[bt] = pd->getNgtBond();
            ft          = pd->getBondFtype();
            break;
        case ebtsANGLES:
            omt->nb[bt] = pd->getNgtAngle();
            ft          = pd->getAngleFtype();
            break;
        case ebtsPDIHS:
            omt->nb[bt] = pd->getNgtDihedral( egdPDIHS);
            ft          = pd->getDihedralFtype( egdPDIHS);
            break;
        case ebtsIDIHS:
            omt->nb[bt] = pd->getNgtDihedral( egdIDIHS);
            ft          = pd->getDihedralFtype( egdIDIHS);
            break;
        default:
            gmx_fatal(FARGS, "Boe");
        }
        
        for (std::vector<alexandria::MyMol>::iterator mymol = mm.begin();
             (mymol < mm.end()); mymol++)
        {
            for (i = 0; (i < mymol->ltop_->idef.il[ft].nr); i += interaction_function[ft].nratoms+1)
            {
                double value, sigma, bondorder;
                int ntrain;
                bool found = false;
                int ai  = mymol->ltop_->idef.il[ft].iatoms[i+1];
                int aj  = mymol->ltop_->idef.il[ft].iatoms[i+2];
                if (pd->atypeToBtype( *mymol->topology_->atoms.atomtype[ai], aai) &&
                    pd->atypeToBtype( *mymol->topology_->atoms.atomtype[aj], aaj))
                { 
                    char buf[STRLEN];
                    switch (bt)
                    {
                    case ebtsBONDS:
                        GtBondConstIterator gtb = pd->findBond(aai, aaj, 0);
                        if (pd->getBondEnd() != gtb)
                        {
                            sprintf(buf, "%s-%s", 
                                    gtb->atom1().c_str(), 
                                    gtb->atom2().c_str());
                            found = true;
                        }
                        break;
                    case ebtsANGLES:
                        ak  = mymol->ltop_->idef.il[ft].iatoms[i+3];
                        if (pd->atypeToBtype( *mymol->topology_->atoms.atomtype[ak], aak))
                        {
                            if (pd->searchAngle(aai, aaj, aak, &value,
                                                &sigma, &ntrain, params))
                            {
                                sprintf(buf, "%s-%s-%s", aai.c_str(), aaj.c_str(), aak.c_str());
                                found = true;
                            }
                        }
                        break;
                    case ebtsPDIHS:
                    case ebtsIDIHS:
                        ak  = mymol->ltop_->idef.il[ft].iatoms[i+3];
                        al  = mymol->ltop_->idef.il[ft].iatoms[i+4];
                        if (pd->atypeToBtype( *mymol->topology_->atoms.atomtype[ak], aak) &&
                            pd->atypeToBtype( *mymol->topology_->atoms.atomtype[al], aal))
                        {
                            if (pd->searchDihedral( (bt == ebtsPDIHS) ? egdPDIHS : egdIDIHS,
                                                    aai, aaj, aak, aal,
                                                    &value, &sigma, &ntrain, params))
                            {
                                sprintf(buf, "%s-%s-%s-%s", aai.c_str(), aaj.c_str(), aak.c_str(), aal.c_str());
                                    found = true;
                            }
                        }
                        break;
                    }
                    if (found)
                    {
                        omt->ngtb[bt][gt-1]++;
                        if (NULL == omt->cgt[bt][gt-1])
                        {
                            omt->cgt[bt][gt-1] = strdup(buf);
                        }
                    }
                }
            }
        }
        for (i = 0; (i < omt->nb[bt]); i++)
        {
            if (omt->ngtb[bt][i] > 0)
            {
                fprintf(fp, "%-8s  %6d  %-20s  %d\n", btsnames[bt], i,
                        omt->cgt[bt][i], omt->ngtb[bt][i]);
                sfree(omt->cgt[bt][i]);
                ntot[bt]++;
            }
        }
        sfree(omt->cgt[bt]);
    }
    for (bt = 0; (bt < ebtsNR); bt++)
    {
        if (bOpt[bt])
        {
            fprintf(fp, "%-8s %d of %4d types\n", btsnames[bt], ntot[bt], omt->nb[bt]);
        }
    }

    return omt;
}

static void check_support(FILE                           *fp,
                          std::vector<alexandria::MyMol> &mm,
                          const Poldata                  &pd,
                          t_commrec                      *cr,
                          bool                            bOpt[])
{
    int ntotal  = (int) mm.size();
    int nlocal  = 0;

    for (std::vector<alexandria::MyMol>::iterator mymol = mm.begin(); (mymol < mm.end()); )
    {
        if (mymol->eSupp != eSupportLocal)
        {
            continue;
        }
        bool bSupport = true;

        for (int bt = 0; bSupport && (bt <= ebtsIDIHS); bt++)
        {
            int  ft;
            if (bOpt[bt])
            {
                switch (bt)
                {
                    case ebtsBONDS:
                        ft = pd.getBondFtype();
                        break;
                    case ebtsANGLES:
                        ft = pd.getAngleFtype();
                        break;
                    case ebtsPDIHS:
                        ft = pd.getDihedralFtype(egdPDIHS);
                        break;
                    case ebtsIDIHS:
                        ft = pd.getDihedralFtype( egdIDIHS);
                        break;
                    default:
                        gmx_fatal(FARGS, "Boe");
                }

                for (int i = 0; bSupport && (i < mymol->ltop_->idef.il[ft].nr); i += interaction_function[ft].nratoms+1)
                {
                    int         ai, aj, ak, al;
                    std::string aai, aaj, aak, aal;

                    ai  = mymol->ltop_->idef.il[ft].iatoms[i+1];
                    aj  = mymol->ltop_->idef.il[ft].iatoms[i+2];
                    if (!(pd.atypeToBtype(*mymol->topology_->atoms.atomtype[ai], aai) &&
                          pd.atypeToBtype(*mymol->topology_->atoms.atomtype[aj], aaj)))
                    {
                        bSupport = false;
                    }
                    switch (bt)
                    {
                    case ebtsBONDS:
                        {
                        double value, sigma, bondorder;
                        int ntrain;
                        std::string params;
                        bSupport = (bSupport &&
                                    pd.searchBond(aai, aaj, &value, &sigma,
                                                   &ntrain, &bondorder, params));
                        break;
                        }
                    case ebtsANGLES:
                        {
                        ak  = mymol->ltop_->idef.il[ft].iatoms[i+3];
                        if (!pd.atypeToBtype( *mymol->topology_->atoms.atomtype[ak], aak))
                        {
                            bSupport = false;
                        }
                        else
                        {
                                bSupport = (bSupport &&
                                            pd.findAngle(aai, aaj, aak) != pd.getAngleEnd());
                        }
                        break;
                        {
                    case ebtsPDIHS:
                    case ebtsIDIHS:
                        {
                        ak  = mymol->ltop_->idef.il[ft].iatoms[i+3];
                        al  = mymol->ltop_->idef.il[ft].iatoms[i+4];
                        if (!(pd.atypeToBtype( *mymol->topology_->atoms.atomtype[ak], aak) &&
                              pd.atypeToBtype( *mymol->topology_->atoms.atomtype[al], aal)))
                        {
                            bSupport = false;
                        }
                        else
                        {
                            int egd = (bt == ebtsPDIHS) ? egdPDIHS : egdIDIHS;
                            bSupport = (bSupport &&
                                        pd.findDihedral(egd, aai, aaj, aak, aal) != pd.getDihedralEnd(egd));
                        }
                        break;
                        }
                    }
                }
            }
        }
        if (!bSupport)
        {
            fprintf(stderr, "No force field support for %s\n",
                    mymol->molProp()->getMolname().c_str());
            mymol = mm.erase(mymol);
        }
        else
        {
            mymol++;
            nlocal++;
        }
    }
    if (PAR(cr))
    {
        gmx_sumi(1, &nlocal, cr);
    }
    if (NULL != fp)
    {
        fprintf(fp, "%d out of %d molecules have support in the force field.\n",
                nlocal, ntotal);
    }
}

static void done_opt_mask(opt_mask_t **omt)
{
    int i;

    for (i = 0; (i < ebtsNR); i++)
    {
        sfree((*omt)->ngtb[i]);
    }
    sfree(*omt);
    *omt = NULL;
}

typedef struct {
    int     gt_index;
    int     natom;
    double  bondorder;
    char  **ai;
    int     nparam;
    double *param;
} opt_bad_t;

static void add_obt(opt_bad_t *obt, int natom, std::string *ai,
                    int nparam, double *param, double bondorder, int gt)
{
    int i;

    obt->gt_index = gt-1;
    obt->natom    = natom;
    snew(obt->ai, natom);
    for (i = 0; (i < natom); i++)
    {
        obt->ai[i] = (char *)ai[i].c_str();
    }
    obt->nparam = nparam;
    snew(obt->param, nparam);
    for (i = 0; (i < nparam); i++)
    {
        obt->param[i] = param[i];
    }
    obt->bondorder = bondorder;
}

namespace alexandria
{

class OptParam : public MolDip
{
    public:
        bool       _bOpt[ebtsNR];
        int        _nbad[ebtsNR];
        opt_bad_t *_bad[ebtsNR];
        int       *_inv_gt[ebtsNR];
        int        _nparam;
        double    *_param, *_orig, *_best, *_lower, *_upper, *_psigma;

        OptParam();
        ~OptParam() {};

        void InitOpt(FILE *fplog, int *nparam,
                     bool bOpt[ebtsNR],
                     real D0, real beta0,
                     opt_mask_t *omt, real factor);

        void MkInvGt();
        void getDissociationEnergy(FILE *fplog);
        void Opt2List();
        void List2Opt();
        void Add(int otype, int natom, std::string *ai,
                 int nparam, double *param, double bondorder, int gt);
        void Print(FILE *fp);
        double CalcDeviation();
        double EnergyFunction(double v[]);
        void GuessAll(int iter, real stepsize,
                      bool bRandom, gmx_rng_t rng);
        void Optimize(FILE *fp, FILE *fplog,
                      int maxiter, //real tol,
                      int nrun, real stepsize, int seed,
                      bool bRandom, const gmx_output_env_t *oenv,
                      int nprint,
                      const char *xvgconv, const char *xvgepot,
                      real temperature);
        void Bayes(FILE *fplog, const char *xvgconv, const char *xvgepot,
                   double start[], double sig[],
                   int n,
                   int nprint,
                   double step,
                   unsigned int seed,
                   real temperature,
                   int    maxiter,
                   double *chi2,
                   const gmx_output_env_t *oenv);
        void PrintSpecs(FILE *fp, char *title,
                        const char *xvg, const gmx_output_env_t *oenv,
                        bool bCheckOutliers);
};

OptParam::OptParam()
{
    int i;

    for (i = 0; (i < ebtsNR); i++)
    {
        _bOpt[i]   = false;
        _nbad[i]   = 0;
        _bad[i]    = NULL;
        _inv_gt[i] = NULL;
    }
    _nparam = 0;
    _param  = NULL;
    _orig   = NULL;
    _best   = NULL;
    _lower  = NULL;
    _upper  = NULL;
    _psigma = NULL;
}

void OptParam::Opt2List()
{
    int i, j, p, n, nparam;

    for (i = n = 0; (i < ebtsNR); i++)
    {
        for (j = 0; (j < _nbad[i]); j++)
        {
            nparam = n+_bad[i][j].nparam;
            if (nparam > _nparam)
            {
                srenew(_param, nparam);
                _nparam = nparam;
            }
            for (p = 0; (p < _bad[i][j].nparam); p++)
            {
                _param[n++] = _bad[i][j].param[p];
            }
        }
    }
}

void OptParam::List2Opt()
{
    int  i, j, p, n;
    char buf[STRLEN];

    for (i = n = 0; (i < ebtsNR); i++)
    {
        for (j = 0; (j < _nbad[i]); j++)
        {
            buf[0] = '\0';
            for (p = 0; (p < _bad[i][j].nparam); p++)
            {
                _bad[i][j].param[p] = _param[n++];
                strcat(buf, " ");
                const char *ptr = gmx_ftoa(_bad[i][j].param[p]).c_str();
                strcat(buf, ptr);
            }
            switch (i)
            {
                case ebtsBONDS:
                    _pd->setBondParams(_bad[i][j].ai[0],
                                       _bad[i][j].ai[1],
                                       0, 0, 0,
                                       _bad[i][j].bondorder, buf);
                    break;
                case ebtsANGLES:
                    _pd->setAngleParams(_bad[i][j].ai[0],
                                        _bad[i][j].ai[1],
                                        _bad[i][j].ai[2],
                                        0, 0, 0, buf);
                    break;
                case ebtsPDIHS:
                    _pd->setDihedralParams(egdPDIHS,
                                           _bad[i][j].ai[0],
                                           _bad[i][j].ai[1],
                                           _bad[i][j].ai[2],
                                           _bad[i][j].ai[3],
                                           0, 0, 0, buf);
                    break;
                case ebtsIDIHS:
                    _pd->setDihedralParams(egdIDIHS,
                                           _bad[i][j].ai[0],
                                           _bad[i][j].ai[1],
                                           _bad[i][j].ai[2],
                                           _bad[i][j].ai[3],
                                           0, 0, 0, buf);
                    break;
                default:
                    gmx_fatal(FARGS, "Unsupported bts %d", i);
            }
        }
    }
}

void OptParam::Add(int otype, int natom, std::string *ai,
                   int nparam, double *param, double bondorder, int gt)
{
    assert(otype < ebtsNR);
    assert(otype >= 0);
    srenew(_bad[otype], ++_nbad[otype]);
    add_obt(&_bad[otype][_nbad[otype]-1], natom, ai, nparam, param, bondorder, gt);
}


void OptParam::MkInvGt()
{
    int i, j, gt_max;

    /* Make in index to look up gt */
    for (i = 0; (i < ebtsNR); i++)
    {
        gt_max = -1;
        for (j = 0; (j < _nbad[i]); j++)
        {
            if (_bad[i][j].gt_index > gt_max)
            {
                gt_max = _bad[i][j].gt_index;
            }
        }
        snew(_inv_gt[i], gt_max+1);
        for (j = 0; (j <= gt_max); j++)
        {
            _inv_gt[i][j] = -666;
        }
        for (j = 0; (j < _nbad[i]); j++)
        {
            _inv_gt[i][_bad[i][j].gt_index] = j;
        }
    }
}

static void dump_csv(int                             nD,
                     char                           *ctest[],
                     std::vector<alexandria::MyMol> &_mymol,
                     double                        **a,
                     double                         *x)
{
    FILE *csv = gmx_ffopen("out.csv", "w");
    fprintf(csv, "\"\",");
    for (int j = 0; (j < nD); j++)
    {
        fprintf(csv, "\"%s\",", ctest[j]);
    }
    fprintf(csv, "\n");
    int i = 0;
    for (std::vector<alexandria::MyMol>::iterator mymol = _mymol.begin();
         (mymol < _mymol.end()); mymol++)
    {
        fprintf(csv, "\"%s\",", mymol->molProp()->getMolname().c_str());
        for (int j = 0; (j < nD); j++)
        {
            fprintf(csv, "%g,", a[i][j]);
        }
        fprintf(csv, "%.3f\n", x[i]);
        i++;
    }
    fclose(csv);

}
void OptParam::getDissociationEnergy(FILE *fplog)
{
    int        nD, nMol, ftb, gt, gti, ai, aj, niter, row;
    char     **ctest;
    char       buf[STRLEN];
    double   **a, **at, **ata;
    double    *x, *atx, *fpp;
    double     a0, da0, ax, chi2;
    int       *test;
    gmx_bool   bZero = FALSE;

    nD     = _nbad[ebtsBONDS];
    nMol   = _mymol.size();
    if ((0 == nD) || (0 == nMol))
    {
        gmx_fatal(FARGS, "Number of variables is %d and number of molecules is %d",
                  nD, nMol);
    }
    a      = alloc_matrix(nMol, nD);
    at     = alloc_matrix(nD, nMol);
    ata    = alloc_matrix(nD, nD);
    snew(x, nMol+1);
    snew(atx, nD+1);
    snew(fpp, nD+1);
    snew(test, nD+1);
    snew(ctest, nD+1);

    fprintf(fplog, "There are %d different bondtypes to optimize the heat of formation\n",
            nD);
    fprintf(fplog, "There are %d (experimental) reference heat of formation.\n", nMol);

    ftb = _pd->getBondFtype();
    int j   = 0;
    for (std::vector<alexandria::MyMol>::iterator mymol = _mymol.begin();
         (mymol < _mymol.end()); mymol++, j++)
    {
        for (int i = 0; (i < mymol->ltop_->idef.il[ftb].nr); i += interaction_function[ftb].nratoms+1)
        {
            ai = mymol->ltop_->idef.il[ftb].iatoms[i+1];
            aj = mymol->ltop_->idef.il[ftb].iatoms[i+2];
            std::string aai = _pd->atypeToBtype(*mymol->topology_->atoms.atomtype[ai]);
            std::string aaj = _pd->atypeToBtype(*mymol->topology_->atoms.atomtype[aj]);
            if ((gt = _pd->searchBond(aai, aaj,
                                      NULL, NULL, NULL, NULL, NULL)) != 0)
            {
                gti = _inv_gt[ebtsBONDS][gt-1];
                at[gti][j]++;
                a[j][gti]++;
                test[gti]++;
                if (NULL == ctest[gti])
                {
                    sprintf(buf, "%s-%s", aai.c_str(), aaj.c_str());
                    ctest[gti] = strdup(buf);
                }
            }
            else
            {
                gmx_fatal(FARGS, "No parameters for bond %s-%s in the force field, atoms %s-%s mol %s",
                          aai.c_str(), aaj.c_str(),
                          *mymol->topology_->atoms.atomtype[ai],
                          *mymol->topology_->atoms.atomtype[aj],
                          mymol->molProp()->getIupac().c_str());
            }
        }
        x[j] = mymol->Emol;
    }

    //double chi22 = multi_regression(stderr, nD, x, nMol, a, fpp);
    //fprintf(stderr, "chi22 = %g\n", chi22);

    matrix_multiply(debug, nMol, nD, a, at, ata);
    dump_csv(nD, ctest, _mymol, a, x);
    if ((row = matrix_invert(debug, nD, ata)) != 0)
    {
        int k = row - 1;
        for (int m = 0; (m < nD); m++)
        {
            if (m == k)
            {
                continue;
            }
            bool   bSame = true;
            double bfac1 = 0, bfac2 = 0;
            for (int l = 0; bSame && (l < nMol); l++)
            {
                if ((a[m][l] != 0) || (a[k][l] != 0))
                {
                    if (a[m][l] != 0)
                    {
                        bfac2 = (1.0*a[k][l])/a[m][l];
                        if ((bfac1 == 0) && (bfac2 != 0))
                        {
                            bfac1 = bfac2;
                        }
                        else if (bfac1 != 0)
                        {
                            bSame = (bfac1 == bfac2);
                        }
                    }
                }
            }
            if (bSame)
            {
                gmx_fatal(FARGS, "Colums %d and %d are identical bfac1 = %g",
                          k + 1, m + 1, bfac1);
            }
        }
        gmx_fatal(FARGS, "Matrix inversion failed. Incorrect column = %d (%s), nD = %d.\nThis probably indicates that you do not have sufficient data points, or that some parameters are linearly dependent.",
                  row, ctest[row-1], nD);
    }
    a0    = 0;
    niter = 0;
    do
    {
        for (int i = 0; (i < nD); i++)
        {
            atx[i] = 0;
            for (int j = 0; (j < nMol); j++)
            {
                atx[i] += at[i][j]*(x[j]-a0);
            }
        }
        for (int i = 0; (i < nD); i++)
        {
            fpp[i] = 0;
            for (int j = 0; (j < nD); j++)
            {
                fpp[i] += ata[i][j]*atx[j];
            }
        }
        da0  = 0;
        chi2 = 0;
        if (bZero)
        {
            for (j = 0; (j < nMol); j++)
            {
                ax = a0;
                for (int i = 0; (i < nD); i++)
                {
                    ax += fpp[i]*a[j][i];
                }
                da0  += (x[j]-ax);
                chi2 += gmx::square(x[j]-ax);
            }
            da0 = da0 / nMol;
            a0 += da0;
            niter++;
            printf("iter: %d, a0 = %g, chi2 = %g\n",
                   niter, a0, chi2/nMol);
        }
    }
    while ((fabs(da0) > 1e-5) && (niter < 1000));

    for (int i = 0; (i < nD); i++)
    {
        _bad[ebtsBONDS][i].param[0] = -fpp[i];
        if (fplog)
        {
            fprintf(fplog, "Optimized dissociation energy for %8s with %4d copies to %g\n",
                    ctest[i], test[i], fpp[i]);
        }
        sfree(ctest[i]);
    }
    sfree(ctest);
    sfree(test);
    sfree(fpp);
    sfree(atx);
    sfree(x);
    free_matrix(a);
    free_matrix(at);
    free_matrix(ata);
}

void OptParam::InitOpt(FILE *fplog, int *nparam,
                       bool bOpt[ebtsNR],
                       real D0, real beta0,
                       opt_mask_t *omt, real factor)
{
    std::string ai[4];
    int         gt, i, n = 0, maxfc = 0;
    double     *fc = NULL;

    for (i = 0; (i < ebtsNR); i++)
    {
        _bOpt[i] = bOpt[i];
    }
    *nparam = 0;
    if (bOpt[ebtsBONDS])
    {
        GtBondIterator bond;
        for (bond = _pd->getBondBegin(), gt = 1;
             bond != _pd->getBondEnd(); bond++, gt++)
        {
            ai[0] = bond->getAtom1();
            ai[1] = bond->getAtom2();
            if (omt->ngtb[ebtsBONDS][gt-1] > 0)
            {
                std::vector<std::string> ptr = splitString(bond->getParams());
                for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                {
                    if (pi->length() > 0)
                    {
                        if (n >= maxfc)
                        {
                            srenew(fc, ++maxfc);
                        }
                        fc[n] = atof(pi->c_str());
                    }
                }
                if (D0 > 0)
                {
                    fc[0] = D0;
                }
                if (beta0 > 0)
                {
                    fc[1] = beta0;
                }
                Add(ebtsBONDS, 2, ai, n, fc, bond->getBondorder(), gt);
                *nparam += n;
            }
        }
    }
    if (bOpt[ebtsANGLES])
    {
        GtAngleIterator angle;
        for (angle = _pd->getAngleBegin(), gt = 1;
             angle != _pd->getAngleEnd(); angle++, gt++)
        {
            ai[0] = angle->getAtom1();
            ai[1] = angle->getAtom2();
            ai[2] = angle->getAtom3();

            if (omt->ngtb[ebtsANGLES][gt-1] > 0)
            {
                std::vector<std::string> ptr = splitString(angle->getParams());
                for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                {
                    if (pi->length() > 0)
                    {
                        if (n > maxfc)
                        {
                            srenew(fc, ++maxfc);
                        }
                        fc[n] = atof(pi->c_str());
                    }
                }
                Add(ebtsANGLES, 3, ai, n, fc, 0, gt);
                *nparam += n;
            }
        }
    }
    if (bOpt[ebtsPDIHS])
    {
        DihedralIterator dihydral;
        for (dihydral = _pd->getDihedralBegin(egdPDIHS), gt = 1;
             dihydral != _pd->getDihedralEnd(egdPDIHS); dihydral++, gt++)
        {
            ai[0]  = dihydral->getAtom1();
            ai[1]  = dihydral->getAtom2();
            ai[2]  = dihydral->getAtom3();
            ai[3]  = dihydral->getAtom4();
            if (omt->ngtb[ebtsPDIHS][gt-1] > 0)
            {
                std::vector<std::string> ptr = splitString(dihydral->getParams());
                for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                {
                    if (pi->length() > 0)
                    {
                        if (n > maxfc)
                        {
                            srenew(fc, ++maxfc);
                        }
                        fc[n] = atof(pi->c_str());
                    }
                }
                Add(ebtsPDIHS, 4, ai, n, fc, 0, gt);
                *nparam += n;
            }
        }
    }
    if (bOpt[ebtsIDIHS])
    {
        DihedralIterator dihydral;
        for (dihydral = _pd->getDihedralBegin(egdIDIHS), gt = 1;
             dihydral != _pd->getDihedralEnd(egdIDIHS); dihydral++, gt++)
        {
            ai[0]  = dihydral->getAtom1();
            ai[1]  = dihydral->getAtom2();
            ai[2]  = dihydral->getAtom3();
            ai[3]  = dihydral->getAtom4();

            if (omt->ngtb[ebtsIDIHS][gt-1] > 0)
            {
                std::vector<std::string> ptr = splitString(dihydral->getParams());
                for (std::vector<std::string>::iterator pi = ptr.begin(); (pi < ptr.end()); ++pi)
                {
                    if (pi->length() > 0)
                    {
                        if (n > maxfc)
                        {
                            srenew(fc, ++maxfc);
                        }
                        fc[n] = atof(pi->c_str());
                    }
                }
                Add(ebtsIDIHS, 4, ai, n, fc, 0, gt);
                *nparam += n;
            }
        }
    }
    sfree(fc);

    MkInvGt();
    Opt2List();
    getDissociationEnergy(fplog);
    Opt2List();
    List2Opt();
    snew(_best, _nparam);
    snew(_orig, _nparam);
    snew(_lower, _nparam);
    snew(_upper, _nparam);
    snew(_psigma, _nparam);
    if (factor < 1)
    {
        factor = 1/factor;
    }
    for (i = 0; (i < _nparam); i++)
    {
        _best[i]  = _orig[i] = _param[i];
        _lower[i] = _orig[i]/factor;
        _upper[i] = _orig[i]*factor;
    }
}

void OptParam::Print(FILE *fp)
{
    int k;

    fprintf(fp, "Param        Orig        Best\n");
    for (k = 0; (k < _nparam); k++)
    {
        fprintf(fp, "%-5d  %10g  %10g\n", k, _orig[k], _best[k]);
    }
}

static void print_stats(FILE *fp, const char *prop, gmx_stats_t lsq, gmx_bool bHeader,
                        char *xaxis, char *yaxis)
{
    real a, da, b, db, chi2, rmsd, Rfit;
    int  n;

    if (bHeader)
    {
        fprintf(fp, "Fitting data to y = ax+b, where x = %s and y = %s\n",
                xaxis, yaxis);
        fprintf(fp, "%-12s %5s %13s %13s %8s %8s\n",
                "Property", "N", "a", "b", "R", "RMSD");
        fprintf(fp, "---------------------------------------------------------------\n");
    }
    gmx_stats_get_ab(lsq, elsqWEIGHT_NONE, &a, &b, &da, &db, &chi2, &Rfit);
    gmx_stats_get_rmsd(lsq, &rmsd);
    gmx_stats_get_npoints(lsq, &n);
    fprintf(fp, "%-12s %5d %6.3f(%5.3f) %6.3f(%5.3f) %7.2f%% %8.4f\n",
            prop, n, a, da, b, db, Rfit*100, rmsd);
}

static void print_lsq_set(FILE *fp, gmx_stats_t lsq)
{
    real   x, y;

    fprintf(fp, "@type xy\n");
    while (gmx_stats_get_point(lsq, &x, &y, NULL, NULL, 0) == estatsOK)
    {
        fprintf(fp, "%10g  %10g\n", x, y);
    }
    fprintf(fp, "&\n");
}

static void xvgr_symbolize(FILE *xvgf, int nsym, const char *leg[],
                           const gmx_output_env_t *oenv)
{
    int i;

    xvgr_legend(xvgf, nsym, leg, oenv);
    for (i = 0; (i < nsym); i++)
    {
        xvgr_line_props(xvgf, i, elNone, ecBlack+i, oenv);
        fprintf(xvgf, "@ s%d symbol %d\n", i, i+1);
    }
}

double OptParam::CalcDeviation()
{
    int             j;
    int             flags;
    double          ener;
    real            t         = 0;
    rvec            mu_tot    = {0, 0, 0};
    tensor          force_vir = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    t_nrnb          my_nrnb;
    gmx_wallcycle_t wcycle;
    FILE           *dbcopy;

    if (PAR(_cr))
    {
        gmx_bcast(sizeof(_bDone), &_bDone, _cr);
        gmx_bcast(sizeof(_bFinal), &_bFinal, _cr);
    }
    if (_bDone)
    {
        return _ener[ermsTOT];
    }
    if (NULL == debug)
    {
        fprintf(debug, "Begin communicating force parameters\n");
        fflush(debug);
    }
    if (PAR(_cr))
    {
        _pd->commForceParameters( _cr);
    }
    if (NULL == debug)
    {
        fprintf(debug, "Done communicating force parameters\n");
        fflush(debug);
    }
    init_nrnb(&my_nrnb);

    wcycle  = wallcycle_init(stdout, 0, _cr);
    for (j = 0; (j < ermsNR); j++)
    {
        _ener[j] = 0;
    }
    flags = GMX_FORCE_NS | GMX_FORCE_LISTED | GMX_FORCE_NONBONDED | GMX_FORCE_FORCES | GMX_FORCE_ENERGY | GMX_FORCE_STATECHANGED;
    flags = ~0;
    for (std::vector<alexandria::MyMol>::iterator mymol = _mymol.begin(); (mymol < _mymol.end()); mymol++)
    {
        if ((mymol->eSupp == eSupportLocal) ||
            (_bFinal && (mymol->eSupp == eSupportRemote)))
        {
            /* Update topology for this molecule */
            mymol->UpdateIdef(_pd, _bOpt);

            /* Now compute energy */
            atoms2md(mymol->mtop_, mymol->inputrec_, 0, NULL, 0,
                     mymol->md_);

            for (j = 0; (j < mymol->molProp()->NAtom()); j++)
            {
                clear_rvec(mymol->f_[j]);
            }

            /* Now optimize the shell positions */
            dbcopy = debug;
            debug  = NULL;
            if (mymol->shellfc_)
            {
                (void)
                    relax_shell_flexcon(debug, _cr, FALSE, 0,
                                    mymol->inputrec_, TRUE, flags,
                                    mymol->ltop_, NULL, mymol->enerd_,
                                    NULL, mymol->state_,
                                    mymol->f_, force_vir, mymol->md_,
                                    &my_nrnb, wcycle, NULL,
                                    &(mymol->mtop_->groups),
                                    mymol->shellfc_, mymol->fr_, FALSE, t, mu_tot,
                                    NULL, NULL);
            }
            else
            {
                do_force(debug, _cr, mymol->inputrec_, 0,
                         &my_nrnb, wcycle, mymol->ltop_,
                         &(mymol->mtop_->groups),
                         mymol->box, mymol->x_, NULL,
                         mymol->f_, force_vir, mymol->md_,
                         mymol->enerd_, NULL,
                         mymol->state_->lambda, NULL,
                         mymol->fr_,
                         NULL, mu_tot, t, NULL, NULL, FALSE,
                         flags);
            }
            debug         = dbcopy;
            mymol->Force2 = 0;
            for (j = 0; (j < mymol->molProp()->NAtom()); j++)
            {
                mymol->Force2 += iprod(mymol->f_[j], mymol->f_[j]);
            }
            mymol->Force2     /= mymol->molProp()->NAtom();
            _ener[ermsForce2] += _fc[ermsForce2]*mymol->Force2;
            mymol->Ecalc       = mymol->enerd_->term[F_EPOT];
            ener               = gmx::square(mymol->Ecalc-mymol->Emol);
            _ener[ermsEPOT]   += _fc[ermsEPOT]*ener/_nmol_support;

            if (NULL != debug)
            {
                fprintf(debug, "%s ener %g Epot %g Force2 %g\n",
                        mymol->molProp()->getMolname().c_str(), ener,
                        mymol->Ecalc, mymol->Force2);
            }
        }
    }
    /* Compute E-bounds */
    for (j = 0; (j < _nparam); j++)
    {
        if (_param[j] < _lower[j])
        {
            _ener[ermsBOUNDS] += _fc[ermsBOUNDS]*gmx::square(_param[j]-_lower[j]);
        }
        else if (_param[j] > _upper[j])
        {
            _ener[ermsBOUNDS] += _fc[ermsBOUNDS]*gmx::square(_param[j]-_upper[j]);
        }
    }

    for (j = 0; (j < ermsTOT); j++)
    {
        _ener[ermsTOT] += _ener[j];
    }

    if (debug)
    {
        fprintf(debug, "ENER:");
        for (j = 0; (j < ermsNR); j++)
        {
            fprintf(debug, "  %8.3f", _ener[j]);
        }
        fprintf(debug, "\n");
    }
    /* Global sum energies */
    if (PAR(_cr))
    {
#if GMX_DOUBLE
        gmx_sumd(ermsNR, _ener, _cr);
#else
        gmx_sumf(ermsNR, _ener, _cr);
#endif
    }
    return _ener[ermsTOT];
}

double OptParam::EnergyFunction(double v[])
{
    int      i;

    /* Copy parameters to topologies */
    for (i = 0; (i < _nparam); i++)
    {
        _param[i] = v[i];
    }
    List2Opt();

    return CalcDeviation();
}

static real guess_new_param(real x, real step, real x0, real x1, gmx_rng_t rng,
                            gmx_bool bRandom)
{
    real r = gmx_rng_uniform_real(rng);

    if (bRandom)
    {
        x = x0+(x1-x0)*r;
    }
    else
    {
        x = x*(1-step+2*step*r);
    }

    if (x < x0)
    {
        return x0;
    }
    else if (x > x1)
    {
        return x1;
    }
    else
    {
        return x;
    }
}

void OptParam::GuessAll(int iter, real stepsize,
                        bool bRandom, gmx_rng_t rng)
{
    double   ppp, xxx;
    gmx_bool bStart = (iter == 0);
    gmx_bool bRand  = bRandom && (iter == 0);
    int      n;

    for (n = 0; (n < _nparam); n++)
    {
        if (bStart)
        {
            ppp = _param[n];
            xxx = guess_new_param(ppp, stepsize, _lower[n], _upper[n], rng, bRand);
            if (bRand)
            {
                _orig[n] = xxx;
            }
            else
            {
                _orig[n] = ppp;
            }
            ppp = xxx;
        }
        else
        {
            ppp = guess_new_param(_orig[n], stepsize, _lower[n],
                                  _upper[n], rng, bRand);
        }
        _param[n] = ppp;
    }
}

void OptParam::Bayes(FILE *fplog, const char *xvgconv, const char *xvgepot,
                     double start[], double sig[],
                     int n,
                     int nprint,
                     double step,
                     unsigned int seed,
                     real temperature,
                     int    maxiter,
                     double *chi2,
                     const gmx_output_env_t *oenv)
{
    int       iter, j, k, nsum, cur = 0;
    double    ds, sorig, DE, E[2] = {0, 0}, beta;
    double   *ssum, *s2sum;
#define prev (1-cur)
    gmx_rng_t rng;
    FILE     *fpc = NULL, *fpe = NULL;

    beta = 1/(BOLTZ*temperature);
    if (NULL != xvgconv)
    {
        fpc = xvgropen(xvgconv, "Parameter convergence", "iteration", "", oenv);
    }
    if (NULL != xvgepot)
    {
        fpe = xvgropen(xvgepot, "Parameter energy", "iteration", "kT", oenv);
    }
    rng = gmx_rng_init(seed);

    E[prev] = EnergyFunction(start);
    *chi2   = E[prev];
    snew(ssum, n);
    snew(s2sum, n);
    nsum = 0;
    for (j = iter = 0; (iter < maxiter); iter++)
    {
        if ((NULL != fpc) && ((j % nprint) == 0))
        {
            fprintf(fpc, "%5d", iter);
            for (k = 0; (k < n); k++)
            {
                fprintf(fpc, "  %10g", start[k]);
            }
            fprintf(fpc, "\n");
        }
        if ((NULL != fpe) && ((j % nprint) == 0))
        {
            fprintf(fpe, "%5d  %10g\n", iter, E[prev]);
        }
        ds        = (2*gmx_rng_uniform_real(rng)-1)*step*fabs(start[j]);
        sorig     = start[j];
        start[j] += ds;
        E[cur]    = EnergyFunction(start);
        DE        = E[cur]-E[prev];
        if (NULL != debug)
        {
            fprintf(debug, "DE = %g ds = %g\n", DE, ds);
        }
        if ((DE < 0) || (exp(-beta*DE) > gmx_rng_uniform_real(rng)))
        {
            if (NULL != debug)
            {
                fprintf(debug, "Changing parameter %3d from %.3f to %.3f. DE = %.3f 'kT'\n",
                        j, sorig, start[j], beta*DE);
            }
            *chi2 = E[cur];
            cur   = prev;
        }
        else
        {
            start[j] = sorig;
        }
        if (iter >= maxiter/2)
        {
            for (k = 0; (k < n); k++)
            {
                ssum[k]  += start[k];
                s2sum[k] += gmx::square(start[k]);
            }
            nsum++;
        }
        j = (j+1) % n;
    }
    gmx_rng_destroy(rng);
    if (NULL != fpc)
    {
        xvgrclose(fpc);
    }
    if (NULL != fpe)
    {
        xvgrclose(fpe);
    }
    if (nsum > 0)
    {
        for (k = 0; (k < n); k++)
        {
            ssum[k]  /= nsum;
            s2sum[k] /= nsum;
        }
    }
    if (NULL != fplog)
    {
        fprintf(fplog, "Average and standard deviation of parameters\n");
        for (k = 0; (k < n); k++)
        {
            sig[k] = sqrt(s2sum[k]-gmx::square(ssum[k]));
            fprintf(fplog, "%5d  %10g  %10g\n",
                    k, ssum[k], sig[k]);
            start[k] = ssum[k];
        }
    }
    sfree(ssum);
    sfree(s2sum);
}

void OptParam::Optimize(FILE *fp, FILE *fplog,
                        int maxiter,
                        int nrun, real stepsize, int seed,
                        bool bRandom, const gmx_output_env_t *oenv,
                        int nprint,
                        const char *xvgconv, const char *xvgepot,
                        real temperature)
{
    double    chi2, chi2_min;
    int       k, n;
    gmx_bool  bMinimum = FALSE;
    gmx_rng_t rng;

    if (MASTER(_cr))
    {
        rng = gmx_rng_init(seed);

        chi2 = chi2_min = GMX_REAL_MAX;
        for (n = 0; (n < nrun); n++)
        {
            if ((NULL != fp) && (0 == n))
            {
                fprintf(fp, "\nStarting run %d out of %d\n", n+1, nrun);
            }

            GuessAll(n, stepsize, bRandom, rng);

            Bayes(fplog, xvgconv, xvgepot, _param, _psigma, _nparam, nprint,
                  stepsize, seed, temperature, maxiter, &chi2, oenv);

            if (chi2 < chi2_min)
            {
                bMinimum = TRUE;
                /* Print convergence if needed */
                for (k = 0; (k < _nparam); k++)
                {
                    _best[k] = _param[k];
                }
                chi2_min   = chi2;
            }

            if (NULL != fp)
            {
                fprintf(fp, "%5d  %8.3f  %8.3f  %8.3f\n", n, chi2, _ener[ermsTOT], _ener[ermsBOUNDS]);
            }
            if (NULL != fplog)
            {
                fprintf(fplog, "%5d  %8.3f  %8.3f  %8.3f\n", n, chi2, _ener[ermsTOT], _ener[ermsBOUNDS]);
                fflush(fplog);
            }
        }

        if (bMinimum)
        {
            for (k = 0; (k < _nparam); k++)
            {
                _param[k] = _best[k];
            }

            double emin = EnergyFunction(_best);
            if (fplog)
            {
                fprintf(fplog, "\nMinimum chi^2 value during optimization: %.3f.\n",
                        chi2_min);
                fprintf(fplog, "\nMinimum RMSD value during optimization: %.3f (kJ/mol).\n",
                        emin);
                //print_opt(fplog,opt);
            }
        }
        CalcDeviation();
        _bDone = TRUE;
        gmx_rng_destroy(rng);
    }
    else
    {
        /* Slave calculators */
        do
        {
            CalcDeviation();
        }
        while (!_bDone);
    }
    CalcDeviation();
}

static void print_moldip_mols(FILE *fp, std::vector<alexandria::MyMol> mol,
                              gmx_bool bForce, gmx_bool bMtop)
{
    int j, k;

    for (std::vector<alexandria::MyMol>::iterator mi = mol.begin(); (mi < mol.end()); mi++)
    {
        fprintf(fp, "%-30s  %d\n", mi->molProp()->getMolname().c_str(), mi->molProp()->NAtom());
        for (j = 0; (j < mi->molProp()->NAtom()); j++)
        {
            fprintf(fp, "  %-5s  %-5s  q = %10g", *(mi->topology_->atoms.atomname[j]),
                    *(mi->topology_->atoms.atomtype[j]), mi->topology_->atoms.atom[j].q);
            if (bForce)
            {
                fprintf(fp, "  %8.3f  %8.3f  %8.3f",
                        mi->f_[j][XX],
                        mi->f_[j][YY],
                        mi->f_[j][ZZ]);
            }
            fprintf(fp, "\n");
        }
        if (bForce)
        {
            for (k = 0; (k < F_NRE); k++)
            {
                if ((mi->enerd_->term[k] != 0) ||
                    (mi->mtop_->moltype[0].ilist[k].nr > 0))
                {
                    fprintf(fp, "%s %d %g\n", interaction_function[k].name,
                            mi->mtop_->moltype[0].ilist[k].nr,
                            mi->enerd_->term[k]);
                }
            }
        }
        if (bMtop)
        {
            pr_mtop(fp, 0, mi->molProp()->getMolname().c_str(), mi->mtop_, TRUE);
        }
    }
}

void OptParam::PrintSpecs(FILE *fp, char *title,
                          const char *xvg, const gmx_output_env_t *oenv,
                          bool bCheckOutliers)
{
    FILE       *xfp;
    int         i;
    double      msd;
    gmx_stats_t gst;

    gst = gmx_stats_init();
    if (NULL != xvg)
    {
        xfp = xvgropen(xvg, "Entalpy of Formation", "Experiment (kJ/mol)", "Calculated (kJ/mol)",
                       oenv);
    }
    fprintf(fp, "%s\n", title);
    fprintf(fp, "Nr.   %-30s %10s %10s %10s %10s %10s\n",
            "Molecule", "DHf@298K", "Emol@0K", "Calc-Exp", "rms F", "Outlier?");
    msd = 0;
    i   = 0;
    for (std::vector<alexandria::MyMol>::iterator mi = _mymol.begin(); (mi < _mymol.end()); mi++, i++)
    {
        real DeltaE = mi->Ecalc - mi->Emol;
        fprintf(fp, "%-5d %-30s %10g %10g %10g %10g %-10s\n",
                i,
                mi->molProp()->getMolname().c_str(),
                mi->Hform, mi->Emol, DeltaE,
                sqrt(mi->Force2),
                (bCheckOutliers && (fabs(DeltaE) > 1000)) ? "XXX" : "");
        msd += gmx::square(mi->Emol-mi->Ecalc);
        gmx_stats_add_point(gst, mi->Hform, mi->Hform + DeltaE, 0, 0);
        if (NULL != xvg)
        {
            fprintf(xfp, "%10g  %10g\n", mi->Hform, mi->Hform + DeltaE);
        }
    }
    fprintf(fp, "\n");
    fprintf(fp, "RMSD is %g kJ/mol for %d molecules.\n\n",
            sqrt(msd/_mymol.size()), (int)_mymol.size());
    fflush(fp);
    if (NULL != xvg)
    {
        xvgrclose(xfp);
        do_view(oenv, xvg, NULL);
    }
    //! Do statistics
    real a, b, da, db, chi2, Rfit;
    int  N;
    gmx_stats_get_ab(gst, 1, &a, &b, &da, &db, &chi2, &Rfit);
    gmx_stats_get_npoints(gst, &N);
    fprintf(fp, "Regression analysis fit to y = ax + b:\n");
    fprintf(fp, "a = %.3f  b = %3f  R2 = %.1f%%  chi2 = %.1f N = %d\n",
            a, b, Rfit*100, chi2, N);
    gmx_stats_free(gst);
    fflush(fp);
}
}

int alex_tune_fc(int argc, char *argv[])
{
    static const char    *desc[] = {
        "tune_fc read a series of molecules and corresponding experimental",
        "heats of formation from a file, and tunes parameters in an algorithm",
        "until the experimental energies are reproduced by the force field.[PAR]",
        "Minima and maxima for the parameters can be set, these are however",
        "not strictly enforced, but rather they are penalized with a harmonic",
        "function, for which the force constant can be set explicitly.[PAR]",
        "At every reinit step parameters are changed by a random amount within",
        "the fraction set by step size, and within the boundaries given",
        "by the minima and maxima. If the [TT]-random[tt] flag is",
        "given a completely random set of parameters is generated at the start",
        "of each run. At reinit steps however, the parameters are only changed",
        "slightly, in order to speed-up local search but not global search."
        "In other words, complete random starts are done only at the beginning of each",
        "run, and only when explicitly requested.[PAR]",
        "The absolut dipole moment of a molecule remains unchanged if all the",
        "atoms swap the sign of the charge. To prevent this kind of mirror",
        "effects a penalty is added to the square deviation ",
        "if hydrogen atoms have a negative charge. Similarly a penalty is",
        "added if atoms from row VI or VII in the periodic table have a positive",
        "charge. The penalty is equal to the force constant given on the command line",
        "time the square of the charge.[PAR]",
        "One of the electronegativities (chi) is redundant in the optimization,",
        "only the relative values are meaningful.",
        "Therefore by default we fix the value for hydrogen to what is written",
        "in the eemprops.dat file (or whatever is given with the [tt]-d[TT] flag).",
        "A suitable value would be 2.3, the original, value due to Pauling,",
        "this can by overridden by setting the [tt]-fixchi[TT] flag to something else (e.g. a non-existing atom).[PAR]",
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };

    t_filenm              fnm[] = {
        { efDAT, "-f", "allmols",    ffREAD  },
        { efDAT, "-d", "gentop",     ffOPTRD },
        { efDAT, "-o", "tune_fc",    ffWRITE },
        { efDAT, "-sel", "molselect", ffREAD },
        { efLOG, "-g", "tune_fc",    ffWRITE },
        { efXVG, "-x", "hform-corr", ffWRITE },
        { efXVG, "-conv", "param-conv", ffWRITE },
        { efXVG, "-epot", "param-epot", ffWRITE }
    };
#define NFILE sizeof(fnm)/sizeof(fnm[0])
    static int            nrun          = 1, maxiter = 100, reinit = 0, seed = 0;
    static int            minimum_data  = 3, compress = 0;
    static real           tol           = 1e-3, stol = 1e-6, watoms = 1;
    static gmx_bool       bRandom       = FALSE, bZero = TRUE, bWeighted = TRUE, bOptHfac = FALSE, bQM = FALSE, bGaussianBug = TRUE, bPol = FALSE, bFitZeta = TRUE;
    static real           J0_0          = 5, Chi0_0 = 1, w_0 = 5, step = 0.01, hfac = 0, rDecrZeta = -1;
    static real           J0_1          = 30, Chi0_1 = 30, w_1 = 50, epsr = 1;
    static real           fc_mu         = 1, fc_bound = 1, fc_quad = 1, fc_charge = 0, fc_esp = 0, fc_epot = 1, fc_force = 0.001;
    static real           factor        = 0.8;
    static char          *opt_elem      = NULL, *const_elem = NULL, *fixchi = (char *)"H";
    static char          *lot           = (char *)"B3LYP/aug-cc-pVTZ";
    static const char    *cqdist[]      = {
        NULL, "AXp", "AXg", "AXs",
        "Yang", "Bultinck", "Rappe", NULL
    };
    static const char    *cqgen[]      = {
        NULL, "None", "EEM", "ESP", "RESP", NULL
    };
    static bool           bOpt[ebtsNR] = { true, false, false, false, false, false };
    static real           beta0        = 0, D0 = 0, beta_min = 10, D0_min = 50, temperature;
    static int            nprint       = 10;
    t_pargs               pa[]         = {
        { "-tol",   FALSE, etREAL, {&tol},
          "Tolerance for convergence in optimization" },
        { "-maxiter", FALSE, etINT, {&maxiter},
          "Max number of iterations for optimization" },
        { "-nprint", FALSE, etINT, {&nprint},
          "How often to print the parameters during the simulation" },
        { "-temp",  FALSE, etREAL, {&temperature},
          "'Temperature' for the Monte Carlo simulation" },
        { "-reinit", FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A vlue of 0 means this is never done at all." },
        { "-stol",   FALSE, etREAL, {&stol},
          "If reinit is -1 then a reinit will be done as soon as the simplex size is below this treshold." },
        { "-nrun",   FALSE, etINT,  {&nrun},
          "This many runs will be done, before each run a complete randomization will be done" },
        { "-qdist",   FALSE, etENUM, {cqdist},
          "Model used for charge distribution" },
        { "-qgen",   FALSE, etENUM, {cqgen},
          "Algorithm used for charge generation" },
        { "-bonds",  FALSE, etBOOL, {&bOpt[ebtsBONDS]},
          "Optimize bond parameters" },
        { "-angles",  FALSE, etBOOL, {&bOpt[ebtsANGLES]},
          "Optimize angle parameters" },
        { "-dihedrals",  FALSE, etBOOL, {&bOpt[ebtsPDIHS]},
          "Optimize proper dihedral parameters" },
        { "-impropers",  FALSE, etBOOL, {&bOpt[ebtsIDIHS]},
          "Optimize improper dihedral parameters" },
        { "-beta0", FALSE, etREAL, {&beta0},
          "Reset the initial beta for Morse potentials to this value, independent of gentop.dat. If value is <= 0 gentop.dat value is used." },
        { "-D0", FALSE, etREAL, {&D0},
          "Reset the initial D for Morse potentials to this value, independent of gentop.dat. If value is <= 0 gentop.dat value is used." },
        { "-beta_min", FALSE, etREAL, {&beta_min},
          "Minimum value for beta in Morse potential" },
        { "-DO_min", FALSE, etREAL, {&D0_min},
          "Minimum value for D0 in Morse potential" },
        { "-qm",     FALSE, etBOOL, {&bQM},
          "Use only quantum chemistry results (from the levels of theory below) in order to fit the parameters. If not set, experimental values will be used as reference with optional quantum chemistry results, in case no experimental results are available" },
        { "-lot",    FALSE, etSTR,  {&lot},
          "Use this method and level of theory when selecting coordinates and charges. Multiple levels can be specified which will be used in the order given, e.g.  B3LYP/aug-cc-pVTZ:HF/6-311G**" },
        { "-fc_bound",    FALSE, etREAL, {&fc_bound},
          "Force constant in the penalty function for going outside the borders given with the above six options." },
        { "-fc_epot",    FALSE, etREAL, {&fc_epot},
          "Force constant in the penalty function for the energy term" },
        { "-fc_force",    FALSE, etREAL, {&fc_force},
          "Force constant in the penalty function for the force term" },
        { "-step",  FALSE, etREAL, {&step},
          "Step size in parameter optimization. Is used as a fraction of the starting value, should be less than 10%. At each reinit step the step size is updated." },
        { "-min_data",  FALSE, etINT, {&minimum_data},
          "Minimum number of data points in order to be able to optimize the parameters for a given atomtype" },
        { "-opt_elem",  FALSE, etSTR, {&opt_elem},
          "Space-separated list of elements to optimize, e.g. \"H C Br\". The other available elements in gentop.dat are left unmodified. If this variable is not set, all elements will be optimized." },
        { "-const_elem",  FALSE, etSTR, {&const_elem},
          "Space-separated list of elements to include but keep constant, e.g. \"O N\". These elements from gentop.dat are left unmodified" },
        { "-seed", FALSE, etINT, {&seed},
          "Random number seed for reinit" },
        { "-factor", FALSE, etREAL, {&factor},
          "Factor for generating random parameters. Parameters will be taken within the limit factor*x - x/factor" },
        { "-random", FALSE, etBOOL, {&bRandom},
          "Generate completely random starting parameters within the limits set by the options. This will be done at the very first step and before each subsequent run." },
        { "-weight", FALSE, etBOOL, {&bWeighted},
          "Perform a weighted fit, by using the errors in the dipoles presented in the input file. This may or may not improve convergence." },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" }
    };
    FILE                 *fp;
    t_commrec            *cr;
    gmx_output_env_t     *oenv;
    gmx_molselect *      gms;
    time_t                my_t;
    opt_mask_t           *omt = NULL;

    cr = init_commrec();
    if (MASTER(cr))
    {
        printf("There are %d threads/processes.\n", cr->nnodes);
    }
    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW,
                           NFILE, fnm,
                           sizeof(pa)/sizeof(pa[0]), pa,
                           sizeof(desc)/sizeof(desc[0]), desc,
                           0, NULL, &oenv))
    {
        return 0;
    }

    if (MASTER(cr))
    {
        fp = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

        time(&my_t);
        fprintf(fp, "# This file was created %s", ctime(&my_t));
        fprintf(fp, "# alexandria is part of G R O M A C S:\n#\n");
        fprintf(fp, "# %s\n#\n", gmx::bromacs().c_str());
    }
    else
    {
        fp = NULL;
    }

    if (MASTER(cr))
    {
        gms = gmx_molselect_init(opt2fn_null("-sel", NFILE, fnm));
    }
    else
    {
        gms = NULL;
    }

    alexandria::OptParam      opt;
    int                       nparam;
    ChargeDistributionModel   iDistributionModel         = Poldata::name2eemtype(cqdist[0]);
    ChargeGenerationAlgorithm iChargeGenerationAlgorithm = (ChargeGenerationAlgorithm) get_option(cqgen);

    opt.Init(cr, bQM, bGaussianBug, iDistributionModel,
             iChargeGenerationAlgorithm,
             rDecrZeta, epsr,
             J0_0, Chi0_0, w_0, J0_1, Chi0_1, w_1,
             fc_bound, fc_mu, fc_quad, fc_charge,
             fc_esp, fc_epot, fc_force, fixchi, bOptHfac, hfac, bPol, bFitZeta);
    if (0 == seed)
    {
        seed = gmx_rng_make_seed();
    }
    opt.Read(fp ? fp : (debug ? debug : NULL),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             minimum_data, bZero,
             opt_elem, const_elem,
             lot, gms, watoms, FALSE, seed);

    check_support(fp, opt._mymol, opt._pd, opt._cr, bOpt);

    fprintf(fp, "In the total data set of %d molecules we have:\n", 
            opt_.mymol.size());
    ForceConstants fc;
    fc.analyzeIdef(fp, opt._mymol, opt._pd, bOpt);

    opt.InitOpt(fp, &nparam, bOpt, D0, beta0,
                omt, factor);

    print_moldip_mols(fp, opt._mymol, FALSE, FALSE);
    omt = analyze_idef(fp, opt._mymol, opt._pd, bOpt);
    if (MASTER(cr))
    {
        opt.PrintSpecs(fp, (char *)"Before optimization", NULL, oenv, false);
    }
    opt.Optimize(MASTER(cr) ? stderr : NULL, fp,
                 maxiter, nrun, step, seed,
                 bRandom, oenv, nprint,
                 opt2fn("-conv", NFILE, fnm),
                 opt2fn("-epot", NFILE, fnm),
                 temperature);
    done_opt_mask(&omt);
    if (MASTER(cr))
    {
        print_moldip_mols(fp, opt._mymol, TRUE, FALSE);
        opt.PrintSpecs(fp, (char *)"After optimization",
                       opt2fn("-x", NFILE, fnm), oenv, true);

        alexandria::PoldataXml::write(opt2fn("-o", NFILE, fnm), opt._pd, compress);

        gmx_molselect_done(gms);
        done_filenms(NFILE, fnm);

        gmx_ffclose(fp);
    }

    return 0;
}
