/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2013,2014,2015,2016, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Implements gmx::analysismodules::Dos.
 *
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 * \ingroup module_trajectoryanalysis
 */
#include "gmxpre.h"

#include "dos.h"

#include <algorithm>
#include <string>
#include <vector>

#include "gromacs/analysisdata/analysisdata.h"
#include "gromacs/analysisdata/modules/average.h"
#include "gromacs/analysisdata/modules/plot.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/correlationfunctions/autocorr.h"
#include "gromacs/correlationfunctions/integrate.h"
#include "gromacs/correlationfunctions/manyautocorrelation.h"
#include "gromacs/fft/fft.h"
#include "gromacs/fileio/confio.h"
#include "gromacs/fileio/oenv.h"
#include "gromacs/fileio/trxio.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/gmxana/gmx_ana.h"
#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vecdump.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/options/basicoptions.h"
#include "gromacs/options/filenameoption.h"
#include "gromacs/options/ioptionscontainer.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/pbcutil/rmpbc.h"
#include "gromacs/random/random.h"
#include "gromacs/topology/atomprop.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/trajectory/trajectoryframe.h"
#include "gromacs/trajectoryanalysis/analysissettings.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/exceptions.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/gmxomp.h"
#include "gromacs/utility/pleasecite.h"
#include "gromacs/utility/programcontext.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/snprintf.h"
#include "gromacs/utility/strdb.h"
#include "gromacs/utility/stringutil.h"

#include "dosutils.h"

namespace gmx
{

namespace analysismodules
{

namespace
{

static const real toler = 1e-5;

struct velName
{
    const char *shortName;
    const char *longName;
};

static velName velocityName[V_NR] = {
    { "Trans",  "Translation" },
    { "Angul",  "Angular motion" },
    { "Vibr",   "Vibration" },
    { "Rotat",  "Rotational motion" },
    { "Atomic", "Atomic motion" },
    { "Sum",    "Sum of it all" }
};

enum Component {
    C_SOLID, C_GAS, C_TOTAL, C_NR
};

static const char *cName[C_NR] = {
    "Solid", "Gas", "Total"
};

static const char *tmName[TM_NR] = {
    "1PT", "2PT", "MD"
};

struct DosNameUnit
{
    const char *shortName;
    const char *longName;
    const char *unitName;
};

static DosNameUnit dosNameUnit[DOS_NR] = {
    { "cV", "Heat capacity",    "(J/mol K)" },
    { "S",  "Entropy",          "(J/mol K)" },
    { "A",  "Helmholtz energy", "(kJ/mol)" },
    { "E",  "Internal energy",  "(kJ/mol)" }
};

//! Helper class for storing data
class DosData
{
    public:
        DosData();

        /*! \brief Resize the internal array
         *
         * \param[in] Ntot the length of the array
         */
        void resizeDosData(int Ntot);
        double            X[DOS_NR][C_NR];
        std::vector<real> dos[C_NR];
        std::vector<real> wdos[DOS_NR][C_NR];
};

DosData::DosData()
{
    for(int i = 0; i < DOS_NR; i++)
    {
        for(int j = 0; j < C_NR; j++)
        {
            X[i][j] = 0;
        }
    }
}

void DosData::resizeDosData(int Ntot)
{
    for (int i = 0; i < C_NR; i++)
    {
        // TODO check whether we need alignment!
        dos[i].resize(Ntot, 0);
        for (int j = 0; (j < DOS_NR); j++)
        {
            X[j][i] = 0;
            wdos[j][i].resize(Ntot, 0);
        }
    }
}

//! Helper class for computing properties for each of components
class DosProps
{
    private:
        int                 velocity_;
        int                 numParticles_;
        std::vector<real>   alldos_;
        double              wDiff_[DOS_NR];
        double              Emd_, E0_;
        double              f_, y_, fy_, ry_, z_, nrho_, Delta_, beta_, bfac_, partMass_;
        double              hsdf_, ttdf_;
        double              sigHS_, Shs_, Sig_;
        double              DiffACF_, DiffDoS_;
        double              dostot_, dof_, DoS0_;
        double              Temperature_;
        /* Velocity components. Dimension DIM*npart x nFrames_ */
        std::vector<real *> vec_;
    public:
        DosData             dd_[LOT_NR][TM_NR];
        
    DosProps();
        ~DosProps();

        int velocity() const { return velocity_; }

        void setVelocity(int velocity) { velocity_ = velocity; }

        void resizeDosProps(int Ntot);

        void analyze(FILE                    *fplog,
                     const std::vector<real> &nu,
                     const std::vector<real> &tt,
                     real                    Volume,
                     const std::vector<real> &dos_vacf,
                     const std::vector<RVec> vecI,
                     const std::string      &fnDos,
                     bool                    bRecip,
                     real                    rot_symm,
                     real                    massSystem,
                     const gmx_output_env_t *oenv);

        void print(FILE *fplog_);

        int numParticles() const { return numParticles_; }

        void setParticles(int numParticles)
        {
            numParticles_ = numParticles;
            vec_.resize(DIM*numParticles_, nullptr);
        }
        double temperature() const { return Temperature_; }

        void setTemperature(double T) { Temperature_ = T; }

        double particleMass() const { return partMass_; }

        void setParticleMass(double mass) { partMass_ = mass; }

        double dosTot() const { return dostot_; }

        void setDostot(double dt) { dostot_ = dt; }

        void setAllDos(int i, double val) { alldos_[i] = val; }

        void incrementAllDos(int i, double val) { alldos_[i] += val; }

        double dos0() const { return DoS0_; }

        void setDos0(double d0) { DoS0_ = d0; }

        double fluidicity() const { return f_; }

        void setFluidicity(double f) { f_ = f; }

        double eMd() const { return Emd_; }

        void setEmd(double emd) { Emd_ = emd; }

        double e0() const { return E0_; }

        void setE0(double e0) { E0_ = e0; }

        void resizeVec(int index, int size) { srenew(vec_[index], size); }

        void setVec(int i, int j, double value) { vec_[i][j] = value; }

        double vec(int i, int j) const { return vec_[i][j]; }

        real **rawVector() { return vec_.data(); }

        double numDoF() const { return dof_; }

        void setDoF(double dof) { dof_ = dof; }

        double beta() const { return beta_; }

        void setBeta(double beta) { beta_ = beta; }

        double bfac() const { return bfac_; }

        void setBfac(double bfac) { bfac_ = bfac; }

        double allDoS(int j) const
        {
            range_check(j, 0, alldos_.size());
            return alldos_[j];
        }

        void printComponent(const char              *fn,
                            bool                     bRecip,
                            const std::vector<real> &nu,
                            const gmx_output_env_t  *oenv);

};

DosProps::DosProps()
{ 
    velocity_ = 0;
    numParticles_ = 0;
    Emd_ = 0;
    E0_ = 0;
    f_ = y_ = fy_ = ry_ = z_ = nrho_ = Delta_ = beta_ = bfac_ = partMass_ = 0;
    hsdf_ = ttdf_ = 0;
    sigHS_ = Shs_ = Sig_ = 0;
    DiffACF_ = DiffDoS_ = 0;
    dostot_ = dof_ = DoS0_ = 0;
    Temperature_ = 0;
    
    for(int d = 0; d<DOS_NR; d++) 
    {
        wDiff_[d] = 0;
    } 
}

DosProps::~DosProps()
{
    for(auto &i : vec_)
    {
        sfree(i);
    }
}

void DosProps::resizeDosProps(int Ntot)
{
    alldos_.resize(Ntot, 0);
    for (int k = 0; (k < LOT_NR); k++)
    {
        for (int l = 0; (l < TM_NR); l++)
        {
            dd_[k][l].resizeDosData(Ntot);
        }
    }
}

void DosProps::analyze(FILE                   *fplog,
                       const std::vector<real> &nu,
                       const std::vector<real> &tt,
                       real                    Volume,
                       const std::vector<real> &dos_vacf,
                       const std::vector<RVec> momentOfInertia,
                       const std::string      &fnDos,
                       bool                    bRecip,
                       real                    rot_symm,
                       real                    massSystem,
                       const gmx_output_env_t *oenv)
{
    setDos0(alldos_[0]);
    if (V_V != velocity_)
    {
        printf("DosProps::analyze numParticles = %d particleMass = %g Volume = %g, T = %g dos0 = %g\n",
               numParticles(), particleMass(), Volume, Temperature_, dos0());
        /* eq 13 JPC B, 114, 8191 (2010) */
        nrho_  = numParticles()/Volume;
        Delta_ = ((2*dos0()/(9*numParticles()))*
                  sqrt(M_PI*BOLTZ*Temperature_/particleMass())*
                  pow(nrho_, 1.0/3.0)*pow(6/M_PI, 2.0/3.0));
        f_     = calc_fluidicity(Delta_, toler);
        ry_    = calc_y(f_, Delta_, toler);

        for (size_t j = 0; (j < nu.size()); j++)
        {
            real dj = (dos0()/(1+gmx::square(dos0()*M_PI*nu[j]/
                                            (6*f_*numParticles()))));
            dd_[LOT_CLASSICAL][TM_2PT].dos[C_GAS][j] = dj;
            dd_[LOT_QUANTUM][TM_2PT].dos[C_GAS][j]   = dj;
            /* Difference between the two */
            dd_[LOT_QUANTUM_CORRECTION][TM_2PT].dos[C_GAS][j] = 0;
        }

        /* Determine hard sphere degrees of freedom hsdf_.
         * Since all the 2PT variants have the same diffusivity, we just take the first.
         */
        real stddev = 0;
        hsdf_  = evaluate_integral(nu.size(), &nu[0],
                                   &dd_[LOT_QUANTUM][TM_2PT].dos[C_GAS][0],
                                   NULL, 0, &stddev);
        /* Determine total degrees of freedom ttdf_ */
        stddev = 0;
        ttdf_  = evaluate_integral(nu.size(), nu.data(), alldos_.data(), 
                                   NULL, 0, &stddev);

        y_     = ry_*(hsdf_/ttdf_);
        fy_    = f_ * y_;
        z_     = calc_compress(y_);
        /* Lin 2010, eq. 16, Sackur-Tetrode equation */

        double xlog = (pow(2*M_PI*particleMass()*BOLTZ*Temperature_/
                           (gmx::square(PLANCK)), 1.5)*Volume/(f_ * numParticles()));
        Sig_   = BOLTZ*(2.5+log(xlog));
        fprintf(fplog, "m %g k %g T %g h %g V %g N %d xlog %g\n",
                particleMass(), BOLTZ, Temperature_, PLANCK, Volume, numParticles(), xlog);

        Shs_   = Sig_ + calc_Shs(y_);
        sigHS_ = pow(6*ry_*Volume/(M_PI*numParticles()), 1.0/3.0);
    }
    else        /* velocity_ == V_V */
    {
        /* All other parameters are 0 */
        y_     = 1;
        for (size_t j = 0; j < nu.size(); j++)
        {
            for (int lot = 0; (lot < LOT_NR); lot++)
            {
                dd_[lot][TM_2PT].dos[C_GAS][j] = 0;
            }
        }
    }
    /* Now compute solid (2) component, for 2PT only.
     * Fill the total arrays for the other methods.
     */
    for (size_t j = 0; (j < nu.size()); j++)
    {
        for (int lot = 0; (lot < LOT_NR); lot++)
        {
            for (int tm = 0; (tm < TM_NR); tm++)
            {
                dd_[lot][tm].dos[C_TOTAL][j] = alldos_[j];
                if (TM_2PT == tm)
                {
                    dd_[lot][tm].dos[C_SOLID][j] = (dd_[lot][tm].dos[C_TOTAL][j]-
                                                    dd_[lot][tm].dos[C_GAS][j]);
                }
            }
        }
    }
    real stddev = 0;
    dostot_ = evaluate_integral(nu.size(), &nu[0], &alldos_[0], NULL, 0, &stddev);

    if (fnDos.size() > 0)
    {
        char buf[256];
        snprintf(buf, sizeof(buf), "%s-%s",
                 velocityName[velocity_].shortName, fnDos.c_str());
        printComponent(buf, bRecip, nu, oenv);
    }
    /* Finally analyze the results! */
    for (int i = 0; (i < DOS_NR); i++)
    {
        wDiff_[i] = 0;
    }
    if (V_A == velocity_)
    {
        const double Q = PLANCK*PLANCK/(8.0*M_PI*M_PI*BOLTZ);
        rvec         I, rT;
        clear_rvec(rT);
        clear_rvec(I);

        for (const auto &mOI : momentOfInertia) 
        {
            for (int k = 0; (k < DIM); k++)
            {
                if (mOI[k] > toler)
                {
                    rT[k] += Q/mOI[k];
                    I[k]  += mOI[k];
                }
            }
        }
        svmul(1.0/numParticles(), I,  I);
        svmul(1.0/numParticles(), rT, rT);

        fprintf(fplog, "Momenta of Inertia Ix %g Iy %g Iz %g amu nm^2\n",
                I[XX], I[YY], I[ZZ]);
        fprintf(fplog, "Characteristic Rotational T x %g y %g z %g K\n",
                rT[XX], rT[YY], rT[ZZ]);

        double crT  = 1.0; /* characteristic rotational temperature */
        int    rdof = 0;   /* rotational degrees of freedom */

        for (int i = 0; i < DIM; i++)
        {
            if (rT[i] > 0.0)
            {
                crT *= rT[i];
                rdof++;
            }
        }
        if (rdof == 0)
        {
            crT = -1.0;
        }
        real SR  = (3.0/2.0 + log(sqrt(M_PI*pow(Temperature_, rdof)/crT)/rot_symm))/3;
        wDiff_[DOS_S]  = SR;
        fprintf(fplog, "SR = %g, rdof = %d, crT = %g rot_symm = %g\n",
                SR, rdof, crT, rot_symm);
    }
    else
    {
        wDiff_[DOS_S] = Shs_/(3*BOLTZ);
    }
    if (V_V != velocity_)
    {
        wDiff_[DOS_CV] = 0.5;
        wDiff_[DOS_E]  = 0.5;
        wDiff_[DOS_A]  = wDiff_[DOS_E] - wDiff_[DOS_S];
    }

    real   fff[DOS_NR];
    int    Nmol = momentOfInertia.size();
    fff[DOS_CV] = fff[DOS_S] = (BOLTZ*KILO/Nmol);
    fff[DOS_A]  = fff[DOS_E] = (BOLTZ*Temperature_/Nmol);
    for (int dos = 0; (dos < DOS_NR); dos++)
    {
        for (size_t j = 0; j < nu.size(); j++)
        {
            for (int lot = 0; lot < LOT_QUANTUM_CORRECTION; lot++)
            {
                for (int tm = 0; (tm < TM_NR); tm++)
                {
                    if (TM_2PT == tm)
                    {
                        dd_[lot][tm].wdos[dos][C_SOLID][j] =
                            (dd_[lot][tm].dos[C_SOLID][j] *
                             weight(dos, lot, tm, C_SOLID, nu[j], beta_));
                        dd_[lot][tm].wdos[dos][C_GAS][j]   =
                            (dd_[lot][tm].dos[C_GAS][j] *
                             wDiff_[dos]);
                        dd_[lot][tm].wdos[dos][C_TOTAL][j] =
                            (dd_[lot][tm].wdos[dos][C_SOLID][j] +
                             dd_[lot][tm].wdos[dos][C_GAS][j]);
                    }
                    else
                    {
                        dd_[lot][tm].wdos[dos][C_TOTAL][j] =
                            (dd_[lot][tm].dos[C_TOTAL][j] *
                             weight(dos, lot, tm, C_TOTAL, nu[j], beta_));
                    }
                }
                for (int tm = 0; (tm < TM_NR); tm++)
                {
                    for (int cc = 0; (cc < C_NR); cc++)
                    {
                        dd_[LOT_QUANTUM_CORRECTION][tm].wdos[dos][cc][j] =
                            (dd_[LOT_QUANTUM][tm].wdos[dos][cc][j] -
                             dd_[LOT_CLASSICAL][tm].wdos[dos][cc][j]);
                    }
                }
            }
        }
        for (int lot = 0; (lot < LOT_NR); lot++)
        {
            for (int tm = 0; (tm < TM_NR); tm++)
            {
                for (int cc = 0; (cc < C_NR); cc++)
                {
                    dd_[lot][tm].X[dos][cc] =
                        fff[dos]*evaluate_integral(nu.size(), nu.data(),
                                                   dd_[lot][tm].wdos[dos][cc].data(),
                                                   NULL, 0, &stddev);
                }
            }
        }

        if (V_V != velocity_)
        {
            real stddev = 0;
            DiffACF_ = (1000/3.0)*evaluate_integral(tt.size(), tt.data(),
                                                    dos_vacf.data(), NULL,
                                                    0, &stddev);
            DiffDoS_ = 1000*dos0()/(12*massSystem*beta_);
        }
    }
}

void DosProps::print(FILE *fplog)
{
    if (numParticles_ == 0)
    {
        fprintf(fplog, "\n No particles in %s\n",
                velocityName[velocity_].longName);
        return;
    }
    fprintf(fplog, "\n+++ Analyzing %s +++\n",
            velocityName[velocity_].longName);
    fprintf(fplog, "Npart                           = %d\n", numParticles_);
    fprintf(fplog, "T                               = %g\n", Temperature_);
    fprintf(fplog, "Number density                  = %g particles/nm^3\n", nrho_);
    fprintf(fplog, "Delta                           = %g\n", Delta_);
    fprintf(fplog, "fluidicity f                    = %g\n", f_);
    fprintf(fplog, "ry                              = %g\n", ry_);
    fprintf(fplog, "hsdf                            = %g\n", hsdf_);
    fprintf(fplog, "ttdf                            = %g\n", ttdf_);
    fprintf(fplog, "fy                              = %g nm\n", fy_);
    fprintf(fplog, "hard sphere packing fraction y  = %g\n", y_);
    fprintf(fplog, "hard sphere compressibility z   = %g\n", z_);
    fprintf(fplog, "Particle mass                   = %g amu\n", partMass_);
    fprintf(fplog, "ideal gas entropy Sig           = %g 1/K\n", Sig_ / BOLTZ);
    fprintf(fplog, "hard sphere entropy Shs         = %g 1/K\n", Shs_ / BOLTZ);
    fprintf(fplog, "sigma_HS                        = %g nm\n", sigHS_);
    fprintf(fplog, "DoS0                            = %g\n", dos0());
    fprintf(fplog, "DoSTot                          = %g\n", dostot_);
    if (V_V != velocity_)
    {
        fprintf(fplog, "Diffusion coefficient from dos[VACF] %g 10^-5 cm^2/s\n",
                DiffACF_);
        fprintf(fplog, "Diffusion coefficient from DoS0 %g 10^-5 cm^2/s\n",
                DiffDoS_);
    }
}

static double reciprocalFactor(bool recip)
{
    if (recip)
    {
        return (1e7/SPEED_OF_LIGHT);
    }
    return 1;
}

void DosProps::printComponent(const char              *fn,
                              bool                     bRecip,
                              const std::vector<real> &nu,
                              const gmx_output_env_t  *oenv)
{
    FILE       *fp;
    int         k;

    fp = xvgropen(fn, "Density of states",
                  bRecip ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)",
                  "\\f{4}S(\\f{12}n\\f{4})", oenv);
    double recip_fac = reciprocalFactor(bRecip);
    xvgr_legend(fp, C_NR, cName, oenv);
    for (size_t j = 0; (j < nu.size()); j++)
    {
        fprintf(fp, "%10g", recip_fac*nu[j]);
        for (k = 0; (k < C_NR); k++)
        {
            fprintf(fp, "  %10g", dd_[LOT_QUANTUM][TM_2PT].dos[k][j]/recip_fac);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}



/*! \brief
 * Class used to perform density of states analysis
 *
 * Inherits TrajectoryAnalysisModule and all functions from there.
 * Does not implement any new functionality but rather uses helper
 * class DosProps to do the real work.
 *
 * \ingroup module_trajectoryanalysis
 */
class Dos : public TrajectoryAnalysisModule
{
    public:
        Dos();
        virtual ~Dos();

        virtual void initOptions(IOptionsContainer          *options,
                                 TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                                  const TopologyInformation        &top);
        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                                  TrajectoryAnalysisModuleData *pdata);
        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        void sumProperties();

        void printResume();

        void printDosPropsLegend(FILE *fp);

        void printAcfs(const std::vector<real> &tt,
                       const std::vector<real> &dos_vacf);

        void printDosGraph(const std::vector<real> &nu);

        int nAtoms() { return atom_.size(); }
        FILE                *fplog_;
        std::string          fnDos_;
        std::string          fnLog_;
        std::string          fnVacf_;
        std::string          fnMvacf_;
        std::string          fnEnergy_;
        bool                 bRecip_;
        bool                 dumpW_;
        bool                 bAtomic_;
        DosProps             dp_[V_NR];
        t_topology          *top_;
        gmx_output_env_t    *oenv_;
        rvec                *r1_;
        rvec                *v1_;
        real                 Emd_;
        real                 Volume_;
        real                 Vsum_;
        real                 V2sum_;
        real                 VolError_;
        int                  nV_;
        int                  nFrames_;
        int                  nAlloc_;
        real                 massSystem_;
        real                 cPclassical_;
        real                 cVclassical_;
        real                 rot_symm_;
        int                  nthreads_;
        int                  NatmXmol_;
        bool                 t0set_;
        double               t0_;
        double               t1_;
        double               dt_;
        std::vector<double **> intensity_;
        std::vector<double **> eigenVector_;
        gmx_rmpbc_t          gRmPbc_;
        std::vector<t_atom>  atom_;
        std::vector<int>     particleIndex_;
        std::vector<int>     dummyIndex_;
        std::vector<real>    massMol_;
        //! Moments of inertia for each molecule
        std::vector<RVec>    MomentOfInertia_;
        // Copy and assign disallowed by base.
};

// Constructor. Here it is important to initialize the pointer to
// subclasses that are elements of the main class. Here we have only
// one. The type of this depends on what kind of tool you need.
// Here we only have simple value/time kind of data.
Dos::Dos() : top_(NULL), oenv_(NULL), r1_(NULL), v1_(NULL)
{
    fplog_ = NULL;
    for (int t = 0; t < V_NR; t++)
    {
        dp_[t].setVelocity(t);
    }
    bRecip_  = false;
    dumpW_   = false;
    bAtomic_ = false;
    Emd_     = 0;
    Volume_  = 0;
    Vsum_    = 0;
    V2sum_   = 0;
    VolError_ = 0;
    nV_       = 0;
    nFrames_ = 0;
    nAlloc_  = 0;
    massSystem_ = 0;
    cPclassical_ = 0;
    cVclassical_ = 0;
    rot_symm_  = 0;
    nthreads_ = 0;
    NatmXmol_ = 0;
    t0set_ = false;
    t0_ = 0;
    t1_ = 0;
    dt_ = 0;
    gRmPbc_ = FALSE;
}

Dos::~Dos()
{
    sfree(v1_);
}

void Dos::printDosPropsLegend(FILE *fp)
{
    const char *leg[V_NR];
    int         nleg = 0;

    for (int t = 0; (t < V_SUM); t++)
    {
        if (dp_[t].numParticles() > 0)
        {
            leg[nleg++] = velocityName[t].shortName;
        }
    }
    xvgr_legend(fp, nleg, leg, oenv_);
}

void Dos::printAcfs(const std::vector<real> &tt,
                    const std::vector<real> &dos_vacf)
{
    if (fnVacf_.size() > 0 && dos_vacf.size() > 0)
    {
        FILE *fp = xvgropen(fnVacf_.c_str(), "Velocity ACF",
                            "Time (ps)", "C(t)", oenv_);
        for (size_t j = 0; (j < dos_vacf.size()); j++)
        {
            fprintf(fp, "%10g %10g\n", tt[j], dos_vacf[j]);
        }
        xvgrclose(fp);
    }
    if (fnMvacf_.size() > 0)
    {
        FILE *fp = xvgropen(fnMvacf_.c_str(), "Mass-weighted velocity ACF",
                            "Time (ps)", "C(t)", oenv_);
        printDosPropsLegend(fp);
        for (size_t j = 0; (j < tt.size()); j++)
        {
            fprintf(fp, "%10g", tt[j]);
            for (int t = 0; (t < V_SUM); t++)
            {
                if ((dp_[t].numParticles() > 0) && (dp_[t].numDoF() > 0))
                {
                    fprintf(fp, "  %10g", dp_[t].allDoS(j));
                }
            }
            fprintf(fp, "\n");
        }
        xvgrclose(fp);
    }
}

void Dos::printDosGraph(const std::vector<real> &nu)
{
    if (fnDos_.size() == 0)
    {
        return;
    }
    FILE *fp = xvgropen(fnDos_.c_str(), "DoS",
                        bRecip_ ? "E (cm\\S-1\\N)" : "\\f{12}n\\f{4} (1/ps)", "DoS(v)",
                        oenv_);

    double recip_fac = reciprocalFactor(bRecip_);
    printDosPropsLegend(fp);
    for (size_t j = 0; (j < nu.size()); j++)
    {
        fprintf(fp, "%10g", recip_fac*nu[j]);
        for (int t = 0; (t < V_SUM); t++)
        {
            if ((dp_[t].numParticles() > 0) && (dp_[t].numDoF() > 0))
            {
                fprintf(fp, " %10g", dp_[t].allDoS(j)/recip_fac);
            }
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

void Dos::sumProperties()
{
    double Ttot = 0;
    for (int t = 0; (t <= V_V); t++)
    {
        for (int dos = 0; (dos < DOS_NR); dos++)
        {
            for (int lot = 0; (lot < LOT_NR); lot++)
            {
                for (int tm = 0; (tm < TM_NR); tm++)
                {
                    for (int cc = 0; (cc < C_NR); cc++)
                    {
                        dp_[V_SUM].dd_[lot][tm].X[dos][cc] +=
                            dp_[t].dd_[lot][tm].X[dos][cc];
                    }
                }
            }
        }
        dp_[V_SUM].setDostot(dp_[V_SUM].dosTot() + dp_[t].dosTot());;
        dp_[V_SUM].setDos0(dp_[V_SUM].dos0() + dp_[t].dos0());
        dp_[V_SUM].setEmd(dp_[V_SUM].eMd() + dp_[t].eMd());
        dp_[V_SUM].setE0(dp_[V_SUM].e0() + dp_[t].e0());
        dp_[V_SUM].setDoF(dp_[V_SUM].numDoF() + dp_[t].numDoF());
        Ttot += dp_[t].temperature()*dp_[t].numDoF();
    }
    dp_[V_SUM].setTemperature(Ttot / dp_[V_SUM].numDoF());
}

void Dos::initOptions(IOptionsContainer          *options,
                      TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "[THISMODULE] computes the Density of States from a constant volume simulation.",
        "From the Density of States as a function of frequency, DoS(v), the",
        "system entropy can be computed and quantum corrections to the heat",
        "capacity at constant volume. A number of methods is used to compute",
        "these properties. Properties based on the DoS(v) are printed in the",
        "log file and recommended values given at the bottom of",
        "this file.[PAR]",
        "In order for this to be work accurately,",
        "the coordinates and velocities must be saved",
        "in the trajectory with sufficiently high frequency such as to cover",
        "all vibrations, typically at most 5 fs.",
        "between saving. [PAR]",
        "It is highly recommended you provide an energy file corresponding to",
        "the simulation performed rather than relying on the command-line options",
        "to enter the MD energy and/or heat capacities."
    };


    // Control input settings
    settings->setFlags(TrajectoryAnalysisSettings::efRequireTop |
                       TrajectoryAnalysisSettings::efUseTopX |
                       TrajectoryAnalysisSettings::efNoUserRmPBC |
                       TrajectoryAnalysisSettings::efNoUserPBC);
    settings->setPBC(true);
    settings->setRmPBC(true);
    settings->setHelpText(desc);
    settings->setFrameFlags(TRX_NEED_X | TRX_NEED_V);

    // Add option for output file
    options->addOption(FileNameOption("enx").filetype(eftEnergy).inputFile()
                           .store(&fnEnergy_).defaultBasename("ener")
                           .description("Energy file corresponding to the simulation trajectory"));
    options->addOption(FileNameOption("g").filetype(eftLog).outputFile()
                           .store(&fnLog_).defaultBasename("dos")
                           .description("Details of DoS calculations")
                           .required());
    options->addOption(FileNameOption("dos").filetype(eftPlot).outputFile()
                           .store(&fnDos_).defaultBasename("dos")
                           .description("Density of states components"));
    options->addOption(FileNameOption("vacf").filetype(eftPlot).outputFile()
                           .store(&fnVacf_).defaultBasename("vacf")
                           .description("Velocity autocorrelation function"));
    options->addOption(FileNameOption("mvacf").filetype(eftPlot).outputFile()
                           .store(&fnMvacf_).defaultBasename("vacf")
                           .description("Mass-weighted velocity autocorrelation function"));

    // Options for controlling behavior of the program
    options->addOption(BooleanOption("atomic").store(&bAtomic_).
                           description("Do the Density of State calculation treating the system as independent atoms as well. This requires a lot of extra memory and computer time and is therefore by default turned off."));
    options->addOption(BooleanOption("recip").store(&bRecip_)
                           .description("Use cm^-1 on X-axis instead of 1/ps for DoS plots."));
    options->addOption(BooleanOption("dumpw").store(&dumpW_)
                           .description("Dump plots from papers and quit."));
    options->addOption(RealOption("Emd").store(&Emd_)
                           .description("The average total energy in the MD simulation"));
    options->addOption(RealOption("cp").store(&cPclassical_)
                           .description("Reference heat capacity at constant pressure. This will be corrected using the quantum corrections to cV due to Berens et al."));
    options->addOption(RealOption("cv").store(&cVclassical_)
                           .description("Reference heat capacity at constant volume. This will be used instead of the one computed from the energy file."));
    options->addOption(RealOption("rot_symm").store(&rot_symm_)
                           .description("Rotational symmetry in the molecule"));
    options->addOption(IntegerOption("nt").store(&nthreads_)
                           .description("Number of threads to run the analysis (<= 0 means determine automatically)"));

}


void
Dos::initAnalysis(const TrajectoryAnalysisSettings &settings,
                  const TopologyInformation        &top)
{
    time_unit_t  time_unit
        = static_cast<time_unit_t>(settings.timeUnit() + 1);
    xvg_format_t xvg = static_cast<xvg_format_t>(settings.plotSettings().plotFormat());
    output_env_init(&oenv_, getProgramContext(), time_unit, FALSE, xvg, 0);

    if (dumpW_)
    {
        double beta = BOLTZ*300;
        dump_w(oenv_, beta);
        exit(0);
    }
    // Extracts number of molecules
    top_ = top.topology();
    GMX_RELEASE_ASSERT(top_->mols.nr > 0, "You need at least one molecule");

    // Open log file etc.
    fplog_ = gmx_ffopen(fnLog_.c_str(), "w");
    fprintf(fplog_, "Doing density of states analysis based on trajectory.\n");
    please_cite(fplog_, "Berens1983a");
    please_cite(fplog_, "Pascal2011a");
    please_cite(fplog_, "Caleman2011b");

    /* Read topology, allocate memory and such. */
    matrix box;
    top.getTopologyConf(&r1_, box);
    snew(v1_, top_->atoms.nr);

    Volume_    = det(box);

    // Check whether all molecules are the same size
    {
        int molsize = top_->mols.index[1]-top_->mols.index[0];
        for (int j = 1; (j < top_->mols.nr); j++)
        {
            if (molsize != (top_->mols.index[j+1]-
                            top_->mols.index[j]))
            {
                gmx_fatal(FARGS, "The system contains more than 1 molecule type.");
            }
        }
    }
    if (fnEnergy_.size() > 0)
    {
        /* Get the energy from the energy file */
        t_energy_term *et;

        snew(et, 3);
        et[0].ftype       = F_PRES;
        et[1].ftype       = F_TEMP;
        et[2].ftype       = F_ETOT;
        et[1].bCheckDrift = et[2].bCheckDrift = TRUE;
        get_energy_terms2(fnEnergy_.c_str(), 3, et);

        if (cVclassical_ == 0)
        {
            cVclassical_ = KILO*((pow(et[2].stddev, 2))/
                                 (top_->mols.nr*BOLTZ*pow(et[1].energy, 2)));
        }
        fprintf(fplog_, "Read energy file %s.\n", fnEnergy_.c_str());
        fprintf(fplog_, "Etot        = %10g +/- %10g\n", et[2].energy, et[2].stddev);
        fprintf(fplog_, "T           = %10g +/- %10g\n", et[1].energy, et[1].stddev);
        fprintf(fplog_, "P           = %10g +/- %10g\n", et[0].energy, et[0].stddev);
        fprintf(fplog_, "cVclassical = %10g\n", cVclassical_);
        if (Emd_ != 0.0)
        {
            fprintf(stderr, "WARNING: You set both -Emd %g and -enx %s.\n", Emd_, fnEnergy_.c_str());
            fprintf(stderr, "Using energy = %g from energy file.\n", et[2].energy);
        }
        Emd_ = et[2].energy;
        sfree(et);
    }
    else if (Emd_ == 0)
    {
        gmx_fatal(FARGS, "Neither the -Emd nor the -enx (preferred) option were set.\n"
                  "You need to supply at least one of these to get correct energies.");
    }
    else
    {
        fprintf(stderr, "WARNING: without energy file no correct cV can be computed.\n");
    }
    /* allocation of memory */

    massSystem_ = 0;

    massMol_.resize(top_->mols.nr);
    for (int k = 0; (k < top_->mols.nr); k++)
    {
        for (int i = top_->mols.index[k]; (i < top_->mols.index[k+1]); i++)
        {
            massMol_[k] += top_->atoms.atom[i].m;
        }
    }
    for (int j = 0; (j < top_->atoms.nr); j++)
    {
        if (eptAtom == top_->atoms.atom[j].ptype)
        {
            /* we need to store the propertie of each atom */
            particleIndex_.push_back(j);
            dummyIndex_.push_back(nAtoms());
            atom_.push_back(top_->atoms.atom[j]);
            massSystem_ += atom_[nAtoms()-1].m;
        }
    }
    NatmXmol_ = nAtoms()/top_->mols.nr;

    // Compute degrees of freedom.
    // TODO: Should be molecule specific!
    int nconstr = ((top_->idef.il[F_CONSTR].nr/(interaction_function[F_CONSTR].nratoms+1)) +
                   3*top_->idef.il[F_SETTLE].nr/(interaction_function[F_SETTLE].nratoms+1))/top_->mols.nr;
    fprintf(fplog_, "Number of constraints in molecule type 0 = %d\n", nconstr);
    dp_[V_ATOM].setDoF(std::max(0, 3*NatmXmol_ - nconstr));
    dp_[V_T].setDoF(std::min(3.0, dp_[V_ATOM].numDoF()));
    dp_[V_A].setDoF(std::min(3.0, dp_[V_ATOM].numDoF() - dp_[V_T].numDoF()));
    if (NatmXmol_ == 1)
    {
        dp_[V_A].setDoF(0);
    }
    else if (NatmXmol_ == 2)
    {
        dp_[V_A].setDoF(std::min(2.0, dp_[V_A].numDoF()));
    }
    dp_[V_R].setDoF(dp_[V_A].numDoF());
    dp_[V_V].setDoF(dp_[V_ATOM].numDoF() - dp_[V_T].numDoF() -
                    dp_[V_A].numDoF());

    /* allocation of memory */
    dp_[V_T].setParticles(top_->mols.nr);
    dp_[V_A].setParticles(top_->mols.nr);
    dp_[V_V].setParticles(nAtoms());
    /* To avoid allocating memory */
    dp_[V_SUM].setParticles(0);
    if (bAtomic_)
    {
        dp_[V_R].setParticles(nAtoms());
        dp_[V_ATOM].setParticles(nAtoms());
    }

    /* Pbc corrections.
     * Note that this routine breaks the mtop structure so it can not be
     * used afterwards!
     */
    int ePBC = guess_ePBC(box);
    gRmPbc_ = gmx_rmpbc_init(&top_->idef, ePBC, top_->atoms.nr);

    fprintf(fplog_, "Natom  %d, Nmol %d, NatmXmol_ %d, NatomTotal %d.\n",
            nAtoms(), top_->mols.nr, NatmXmol_, top_->atoms.nr);

    nAlloc_  = 0;
    nFrames_ = 0;
    Vsum_    = V2sum_ = VolError_ = 0;
    nV_      = 0;

    if (nthreads_ > 1)
    {
        gmx_omp_set_num_threads(nthreads_);
    }
    else
    {
        nthreads_  = gmx_omp_get_max_threads();
    }
    fprintf(fplog_, "Running the analysis on %d threads\n", nthreads_);
    /* Allocate temp variables for principal component calc */
    intensity_.resize(nthreads_);
    eigenVector_.resize(nthreads_);
    for(int i =0; i < nthreads_; i++)
    {
        snew(intensity_[i], NDIM);
        snew(eigenVector_[i], NDIM);
        for(int j = 0; (j<NDIM); j++)
        {
            snew(intensity_[i][j], NDIM);
            snew(eigenVector_[i][j], NDIM);
        }
    }
    MomentOfInertia_.resize(top_->mols.nr);
}

void
Dos::analyzeFrame(int /*frnr*/, const t_trxframe &fr, t_pbc *pbc,
                  TrajectoryAnalysisModuleData */*pdata*/)
{
    if (!fr.bX || !fr.bV)
    {
        return;
    }
    nFrames_++;
    if (!t0set_)
    {
        t0_    = fr.time;
        t0set_ = true;
    }
    t1_ = fr.time;

    GMX_RELEASE_ASSERT(NULL != pbc,
                       "You have no periodic boundary conditions");

    if (fr.bBox)
    {
        Volume_ = det(fr.box);
        V2sum_ += Volume_*Volume_;
        Vsum_  += Volume_;
        nV_++;
    }

    if (nFrames_ >= nAlloc_)
    {
        nAlloc_ += 12800;
        for (int t = 0; (t < V_NR); t++)
        {
#pragma omp parallel
            {
                int i0, i1, thread_id, NN;

                NN        = dp_[t].numParticles()*DIM;
                thread_id = gmx_omp_get_thread_num();
                i0        = thread_id*NN/nthreads_;
                i1        = std::min(NN, (thread_id+1)*NN/nthreads_);

                for (int i = i0; i < i1; i++)
                {
                    dp_[t].resizeVec(i, nAlloc_);
                }
            }
            // End parallel section.
        }
        /* Remove periodicity */
        matrix box;
        copy_mat(fr.box, box);
        gmx_rmpbc_copy(gRmPbc_, top_->atoms.nr, box, fr.x, fr.x);

        /* Read velocities and coordinates of all atoms:
         * each atoms is pointed to by j
         */
#pragma omp parallel
        {
            int thread_id = gmx_omp_get_thread_num();
            int j0        = thread_id*nAtoms()/nthreads_;
            int j1        = std::min(nAtoms(), (thread_id+1)*nAtoms()/nthreads_);

            for (int j = j0; j < j1; j++)
            {
                int  aj  = particleIndex_[j];

                /* load atomic positions and velocities */
                copy_rvec(fr.v[aj], v1_[j]);
                copy_rvec(fr.x[aj], r1_[j]);

                if (bAtomic_)
                {
                    for (int k = 0; (k < DIM); k++)
                    {
                        /* save total velocities */
                        dp_[V_ATOM].setVec(DIM*j+k, nFrames_, v1_[j][k]);
                    }
                }
            }
        }
        // End parallel section
        /* Dividing atoms in molecules, id. Do it in parallel using OpenMP */
#pragma omp parallel
        {
            int    i0, i1, thread_id;
            matrix trans;
            rvec   tmp, I, angmom, pomega, omega, vcm, xcm, vr;

            thread_id = gmx_omp_get_thread_num();
            i0        = (thread_id*top_->mols.nr)/nthreads_;
            i1        = std::min(top_->mols.nr,
                                 ((thread_id+1)*top_->mols.nr)/nthreads_);

            for (int id = i0; id < i1; id++)
            {
                int  k, lid  = id*DIM;
                int *d_index = &dummyIndex_[id*NatmXmol_];

                /* Find center of mass and shift all position to this new origin:
                 * it must be done for each atom in the molecule.
                 * Returns total mass of molecule
                 */
                dsub_xcm(r1_, NatmXmol_, d_index, atom_.data(), xcm);
                /* find the velocity of center of mass and shift all position to
                 * this new origin: must be done for each molecule
                 */
                dsub_xcm(v1_, NatmXmol_, d_index, atom_.data(), vcm);
                /* save translational velocities for the molecule id */
                for (k = 0; (k < DIM); k++)
                {
                    dp_[V_T].setVec(lid+k, nFrames_, vcm[k]);
                }

                clear_mat(trans);
                clear_rvec(I);
                if (NatmXmol_ > 1.0)
                {
                    /* compute principal moment of inertia
                     * return trans eigenvector tensor and I the eigenvalues
                     * (principal moment of inertia)
                     */
                    principal(NatmXmol_, d_index, atom_.data(), r1_, trans, I,
                              intensity_[thread_id], eigenVector_[thread_id]);
                    /* Save the principal moment: it will be used to compute the
                     * rotational temperatures.
                     */
                    for(int m = 0; m < DIM; m++)
                    {
                        MomentOfInertia_[id][m] += I[m];
                    }
                    /* we run through the j-th atom of id-th molecule and we
                     * save velocities for i-th atoms from d_index
                     * reset the needed variables
                     */
                    clear_rvec(angmom);
                    clear_rvec(pomega);
                    clear_rvec(omega);

                    for (int j = 0; j < NatmXmol_; j++)
                    {
                        int ai = d_index[j];
                        /* compute the angular momentum angmom */
                        clear_rvec(tmp);
                        cprod(r1_[ai], v1_[ai], tmp);
                        /* mass weigth it for each atom in the molecule */
                        angmom[XX] += atom_[ai].m*tmp[XX];
                        angmom[YY] += atom_[ai].m*tmp[YY];
                        angmom[ZZ] += atom_[ai].m*tmp[ZZ];
                    }
                    /* compute angular velocity along principle axis, first step */
                    for (int j = 0; (j < DIM); j++)
                    {
                        if (I[j] > 0.0)
                        {
                            for (k = 0; k < 3; k++)
                            {
                                pomega[j] += angmom[k]*trans[k][j];
                            }
                            pomega[j] /= I[j];
                        }
                    }
                    /* calculate angular velocities. Here we use the transpose of the trans
                     * matrix by swapping the indexing
                     */
                    for (int j = 0; j < DIM; j++)
                    {
                        for (int k = 0; k < DIM; k++)
                        {
                            omega[j] += pomega[k]*trans[j][k];
                        }
                    }
                    /* calculate inertia weighted angular velocities and save */
                    for (int j = 0; j < DIM; j++)
                    {
                        real vvv = 0;
                        for (int k = 0; k < DIM; k++)
                        {
                            if (I[k] > 0)
                            {
                                vvv += pomega[k]*trans[j][k]*sqrt(I[k]);
                            }
                        }
                        /* save angular velocities of molecule id,
                         * weighted by sqrt(I)
                         */
                        dp_[V_A].setVec(id*DIM+j, nFrames_, vvv);
                    }

                    for (int j = 0; j < NatmXmol_; j++)
                    {
                        int ai = d_index[j];
                        /* calculate velocity due to rotation vr = w x r */
                        cprod(omega, r1_[ai], vr);
                        for (int k = 0; k < DIM; k++)
                        {
                            if (bAtomic_)
                            {
                                /* save rotational velocities */
                                dp_[V_R].setVec(ai*DIM+k, nFrames_, vr[k]);
                            }
                            /* calculate vibrational velocities and save
                             * v1_ is the relative velocity versus center of mass
                             */
                            dp_[V_V].setVec(ai*DIM+k, nFrames_, v1_[ai][k] - vr[k]);
                        }
                    }
                }
            }
        }
        // End parallel section
    }
}

void Dos::printResume()
{
    fprintf(fplog_, "\n\t\t Final resume \n");
    fprintf(fplog_, "System  = %s\n", *top_->name);
    fprintf(fplog_, "Nmol    = %d\n", top_->mols.nr);
    fprintf(fplog_, "Natom   = %d\n", nAtoms());
    fprintf(fplog_, "dt      = %g ps\n", dt_);
    fprintf(fplog_, "tmass   = %g amu\n", massSystem_);
    if (VolError_/Volume_ > toler)
    {
        fprintf(fplog_, "Volume  = %g +/- %g nm^3\n", Volume_, VolError_);
        fprintf(fplog_, "WARNING: The two phase thermodynamics method may not work well\n");
        fprintf(fplog_, "in constant pressure simulations.\n");
    }
    fprintf(fplog_, "Volume  = %g nm^3\n", Volume_);
    double rho = (massSystem_*AMU)/(Volume_*NANO*NANO*NANO);
    fprintf(fplog_, "Density = %g g/l\n", rho);
    fprintf(fplog_, "Emd     = %g kJ/mol/system\n", Emd_);
    fprintf(fplog_, "bRecip  = %s\n", boolToString(bRecip_));

#define NITEM 8
    const char *items[NITEM] = {
        "", "Temperature", "DOF", "Dos integral", "DoS0",
        "E0", "Emd", "fluidicity",
    };

    for (int k = 0; (k < NITEM); k++)
    {
        fprintf(fplog_, "%-29s", items[k]);
        for (int t = 0; (t <= V_SUM); t++)
        {
            if ((t != V_SUM) &&
                ((dp_[t].numDoF() == 0) || (dp_[t].numParticles() == 0)))
            {
                continue;
            }
            switch (k)
            {
                case 0:
                    fprintf(fplog_, "  %10s", velocityName[t].shortName);
                    break;
                case 1:
                    fprintf(fplog_, "  %10g", dp_[t].temperature());
                    break;
                case 2:
                    fprintf(fplog_, "  %10g", dp_[t].numDoF());
                    break;
                case 3:
                    fprintf(fplog_, "  %10g", dp_[t].dosTot());
                    break;
                case 4:
                    fprintf(fplog_, "  %10g", dp_[t].dos0());
                    break;
                case 5:
                    fprintf(fplog_, "  %10g", dp_[t].e0());
                    break;
                case 6:
                    fprintf(fplog_, "  %10g", dp_[t].eMd());
                    break;
                case 7:
                    fprintf(fplog_, "  %10g", dp_[t].fluidicity());
                    break;
            }
        }
        fprintf(fplog_, "\n");
    }
    for (int dos = 0; (dos < DOS_NR); dos++)
    {
        for (int tm = 0; (tm < TM_NR); tm++)
        {
            for (int lot = 0; (lot < LOT_NR); lot++)
            {
                int j0 = C_TOTAL;
                if (TM_2PT == tm)
                {
                    j0 = 0;
                }
                for (int j = j0; (j < C_NR); j++)
                {
                    fprintf(fplog_, "%-2s  %-10s  %-6s %-6s",
                            dosNameUnit[dos].shortName,
                            lotName(static_cast<LevelOfTheory>(lot)),
                            tmName[tm], cName[j]);
                    for (int t = 0; (t <= V_SUM); t++)
                    {
                        if ((t != V_SUM) &&
                            ((dp_[t].numDoF() == 0) || (dp_[t].numParticles() == 0)))
                        {
                            continue;
                        }
                        fprintf(fplog_, "  %10g", dp_[t].dd_[lot][tm].X[dos][j]);
                    }
                    fprintf(fplog_, " %s\n", dosNameUnit[dos].unitName);
                }
            }
        }
    }

    fprintf(fplog_, "\nRecommended thermodynamics values:\n");
    fprintf(fplog_, "----------------------------------\n");
    fprintf(fplog_, "%-16s %2s = %10.3f %s (%s/%s method)\n",
            dosNameUnit[DOS_S].longName,
            dosNameUnit[DOS_S].shortName,
            dp_[V_SUM].dd_[LOT_QUANTUM][TM_2PT].X[DOS_S][C_TOTAL],
            dosNameUnit[DOS_S].unitName,
            tmName[TM_2PT], lotName(LOT_QUANTUM));

    if (cVclassical_ > 0)
    {
        fprintf(fplog_, "%-16s %2s = %10.3f %s (%s/%s method)\n",
                dosNameUnit[DOS_CV].longName,
                dosNameUnit[DOS_CV].shortName,
                cVclassical_+dp_[V_SUM].dd_[LOT_QUANTUM_CORRECTION][TM_1PT].X[DOS_CV][C_TOTAL],
                dosNameUnit[DOS_CV].unitName,
                tmName[TM_MD], lotName(LOT_QUANTUM_CORRECTION));
        fprintf(fplog_, "Note: that heat capacities in MD simulations converge poorly.\n");
        fprintf(fplog_, "Please check the convergence using gmx fluctprops.\n");
    }
    else
    {
        fprintf(fplog_, "No cV based on the MD simulation. Please rerun with the -enx option.\n");
    }
    if (cPclassical_ > 0)
    {
        fprintf(fplog_, "%-16s %2s = %10.3f %s (%s/%s method)\n",
                "Heat capacity", "cP",
                cPclassical_+dp_[V_SUM].dd_[LOT_QUANTUM_CORRECTION][TM_1PT].X[DOS_CV][C_TOTAL],
                dosNameUnit[DOS_CV].unitName,
                tmName[TM_MD], lotName(LOT_QUANTUM_CORRECTION));
    }


    fprintf(fplog_, "\nNote that the Helmholtz energy and the internal energy\n");
    fprintf(fplog_, "include the MD energy which usually is on an unknown scale\n");
    fprintf(fplog_, "meaning the energies can not be compared to experiment.\n");
    fprintf(fplog_, "%-16s %2s = %10.3f %s (%s/%s method)\n",
            dosNameUnit[DOS_A].longName,
            dosNameUnit[DOS_A].shortName,
            dp_[V_SUM].dd_[LOT_QUANTUM][TM_2PT].X[DOS_A][C_TOTAL],
            dosNameUnit[DOS_A].unitName,
            tmName[TM_2PT], lotName(LOT_QUANTUM));

    fprintf(fplog_, "%-16s %2s = %10.3f %s (%s/%s method)\n",
            dosNameUnit[DOS_E].longName,
            dosNameUnit[DOS_E].shortName,
            dp_[V_SUM].eMd() + dp_[V_SUM].dd_[LOT_QUANTUM_CORRECTION][TM_1PT].X[DOS_E][C_TOTAL],
            dosNameUnit[DOS_E].unitName,
            tmName[TM_MD], lotName(LOT_QUANTUM_CORRECTION));
    fprintf(fplog_, "\nArrivederci!\n");
}

void
Dos::finishAnalysis(int /* nframes */)
{
    if (nFrames_ <= 4)
    {
        fprintf(stderr, "Not enough frames (%d) for analysis.\n", nFrames_);
        return;
    }

    std::vector<real> nu, tt, dos_vacf;
    double            fudge = 1.5;
    int               N_tot = fudge*nFrames_;
    int               nfreq = 1+N_tot/2;
    printf("N_tot = %d nAlloc_ = %d\n", N_tot, nAlloc_);

    gmx_rmpbc_done(gRmPbc_);

    /* Normalize moment of inertia */
    real ddd = 1.0/nFrames_;
    rvec sumVec;
    clear_rvec(sumVec);
    for (const auto &v : MomentOfInertia_)
    {
        rvec v2;
        svmul(ddd, v, v2);
        rvec_inc(sumVec, v2);
    }
    if (NULL != debug)
    {
        pr_rvecs(debug, 0, "Average MomentOfInertia_", as_rvec_array(MomentOfInertia_.data()), MomentOfInertia_.size());
        svmul(1.0/MomentOfInertia_.size(), sumVec, sumVec);
        pr_rvec(debug, 0, "Overall Average momentOfInertia", sumVec, DIM, FALSE);
    }

    dt_ = (t1_-t0_)/(nFrames_-1);
    if (nV_ > 0)
    {
        double V2aver = V2sum_/nV_;
        double Vol2   = Volume_*Volume_;
        Volume_   = Vsum_/nV_;
        VolError_ = 0;
        if (V2aver > Vol2)
        {
            VolError_ = sqrt(V2aver - Vol2);
        }
        fprintf(fplog_, "Vsum = %g V2sum = %g Volume = %g VolError = %g\n",
                Vsum_, V2sum_, Volume_, VolError_);
    }

    dos_vacf.resize(N_tot);
    int nfft = 0;
    for (int t = 0; (t < V_NR); t++)
    {
        /* Set the particle mass */
        if (dp_[t].numParticles() > 0)
        {
            dp_[t].setParticleMass(massSystem_/dp_[t].numParticles());
        }

        /* Padding the vec array with zeroes */
#pragma omp parallel
        {
            int i0, i1, thread_id, NN;

            NN        = DIM*dp_[t].numParticles();
            thread_id = gmx_omp_get_thread_num();
            i0        = thread_id*NN/nthreads_;
            i1        = std::min(NN, (thread_id+1)*NN/nthreads_);
            for (int i = i0; (i < i1); i++)
            {
                int j;
                if (N_tot > nAlloc_)
                {
                    dp_[t].resizeVec(i, N_tot);
                }
                for (j = nFrames_; (j < N_tot); j++)
                {
                    dp_[t].setVec(i, j, 0);
                }
            }
        }
        // End parallel section
        nfft += DIM*dp_[t].numParticles();
    }
    printf("Going to determine %d correlation functions of length %d. Hang on.\n",
           nfft, N_tot);
    printf("Going to merge the ACFs into the mass-weighted and plain ACF\n");

    /* Allocate arrays */
    for (int t = 0; (t < V_NR); t++)
    {
        dp_[t].resizeDosProps(N_tot);
    }

    for (int t = 0; (t < V_NR); t++)
    {
        if ((dp_[t].numDoF() > 0) && (dp_[t].numParticles() > 0))
        {
            int fftcode;
            fftcode = many_auto_correl(DIM*dp_[t].numParticles(),
                                       nFrames_, N_tot, dp_[t].rawVector());
            if (0 != fftcode)
            {
                GMX_THROW(APIError("FFT failed"));
            }
#pragma omp parallel
            //TODO: Swap i and j loops for parallellization.
            {
                int j0, j1, thread_id;

                thread_id = gmx_omp_get_thread_num();
                j0        = thread_id*N_tot/nthreads_;
                j1        = std::min(N_tot, (thread_id+1)*N_tot/nthreads_);

                for (int j = j0; (j < j1); j++)
                {
                    for (int i = 0; (i < dp_[t].numParticles()); i++)
                    {
                        int  ai    = i*DIM;
                        real v1j, mass  = atom_[i].m;
                        if (V_T == t)
                        {
                            /* NOTE: should be molecule specific */
                            mass = massMol_[0];
                        }
                        else if (V_A == t)
                        {
                            mass = 1;
                        }

                        v1j = (dp_[t].vec(ai+XX, j) +
                               dp_[t].vec(ai+YY, j) +
                               dp_[t].vec(ai+ZZ, j));
                        dp_[t].incrementAllDos(j, v1j*mass);
                        if (V_ATOM == t)
                        {
                            dos_vacf[j] += v1j/dp_[t].numParticles();
                        }
                    }
                }
            }
            // End parallel section
        }
    }

    /* Compute the temperatures. The number of degrees of freedom
     * are per molecule, therefore we divide the DoS by the number
     * of molecules.
     * From here, the temperature is taken from simulation.
     */
    for (int t = 0; (t < V_NR); t++)
    {
        if ((dp_[t].numDoF() > 0) && (dp_[t].numParticles() > 0))
        {
            dp_[t].setTemperature(dp_[t].allDoS(0)/(dp_[t].numDoF()*BOLTZ*top_->mols.nr));
            dp_[t].setBeta(1.0/(dp_[t].temperature()*BOLTZ));
            dp_[t].setBfac(2.0*dt_*dp_[t].beta());
        }
    }
    /* Make time- and frequency arrays */
    tt.resize(N_tot);
    nu.resize(N_tot);
    {
        double pwrfreq = 1/(dt_*(N_tot+1));
        for (int j = 0; (j < N_tot); j++)
        {
            tt[j] = j*dt_;
            nu[j] = j*pwrfreq;
        }
    }
    printAcfs(tt, dos_vacf);

    // Time to prepare for ffts.
    int       fftcode;
    gmx_fft_t fft2;
    if ((fftcode = gmx_fft_init_1d(&fft2, N_tot,
                                   GMX_FFT_FLAG_CONSERVATIVE)) != 0)
    {
        gmx_fatal(FARGS, "gmx_fft_init_1d returned %d", fftcode);
    }

    /* compute density of state for each kind of velocities,
     * return in same vector
     */
    for (int t = 0; t < V_NR; t++)
    {
        if ((dp_[t].numDoF() > 0) && (dp_[t].numParticles() > 0))
        {
            typedef real complex[2];
            complex *in, *out;
            int      j;
            snew(in, N_tot);
            snew(out, N_tot);
            for (j = 0; (j < N_tot); j++)
            {
                in[j][0] = dp_[t].bfac()*dp_[t].allDoS(j);
                in[j][1] = 0;
            }
            for (; (j < N_tot); j++)
            {
                in[j][0] = in[j][1] = 0;
            }
            /* The 2pt code uses a backword transform, but the articles
             * write a forward transform. However for real data it
             * does not affect the real component of the output.
             */
            if ((fftcode = gmx_fft_1d(fft2, GMX_FFT_BACKWARD,
                                      (void *)in, (void *)out)) != 0)
            {
                gmx_fatal(FARGS, "gmx_fft_1d_real returned %d", fftcode);
            }
            for (j = 0; (j < nfreq); j++)
            {
                dp_[t].setAllDos(j, sqrt(gmx::square(out[j][0])+gmx::square(out[j][1])));
            }
            sfree(in);
            sfree(out);
        }
    }
    gmx_fft_destroy(fft2);

    printDosGraph(nu);

    /* Analyze the DoS[t] */
    for (int t = 0; (t < V_NR); t++)
    {
        if ((dp_[t].numDoF() == 0) || (dp_[t].numParticles() == 0))
        {
            continue;
        }
        dp_[t].analyze(fplog_, nu, tt, Volume_,
                       dos_vacf, MomentOfInertia_, fnDos_, bRecip_,
                       rot_symm_, massSystem_, oenv_);

    }
    for (int t = 0; (t <= V_ATOM); t++)
    {
        double E0[C_NR];
        double dof_frac;

        dof_frac  = (dp_[t].numDoF()/dp_[V_ATOM].numDoF());
        dp_[t].setEmd(Emd_*dof_frac/top_->mols.nr);
        dp_[t].setE0(dp_[t].eMd() - dp_[t].dd_[LOT_CLASSICAL][TM_1PT].X[DOS_E][C_TOTAL]);

        E0[C_SOLID] = (1-dp_[t].fluidicity())*dp_[t].e0();
        E0[C_GAS]   = dp_[t].fluidicity()*dp_[t].e0();
        E0[C_TOTAL] = dp_[t].e0();
        /* Since we add E0 to both QM and Classical
         * the difference to be added to the quantum corrections is 0.
         */
        for (int lot = 0; (lot < LOT_QUANTUM_CORRECTION); lot++)
        {
            for (int tm = 0; (tm < TM_NR); tm++)
            {
                for (int cc = 0; (cc < C_NR); cc++)
                {
                    dp_[t].dd_[lot][tm].X[DOS_A][cc] += E0[cc];
                    dp_[t].dd_[lot][tm].X[DOS_E][cc] += E0[cc];
                }
            }
        }
        dp_[t].print(fplog_);
    }

    /* Sum the energies etc. */
    sumProperties();
}

void
Dos::writeOutput()
{
    printResume();
    fclose(fplog_);
    do_view(oenv_, fnDos_.c_str(), "-nxy");
}

}       // namespace

const char DosInfo::name[]             = "dos";
const char DosInfo::shortDescription[] =
    "Calculate thermodynamics based on density of states";

TrajectoryAnalysisModulePointer DosInfo::create()
{
    return TrajectoryAnalysisModulePointer(new Dos);
}

} // namespace analysismodules

} // namespace gmx
