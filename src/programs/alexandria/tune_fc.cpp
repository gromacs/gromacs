/*
 * This source file is part of the Alexandria program.
 *
 * Copyright (C) 2014-2019
 *
 * Developers:
 *             Mohammad Mehdi Ghahremanpour,
 *             Paul J. van Maaren,
 *             David van der Spoel (Project leader)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA  02110-1301, USA.
 */

/*! \internal \brief
 * Implements part of the alexandria program.
 * \author Mohammad Mehdi Ghahremanpour <mohammad.ghahremanpour@icm.uu.se>
 * \author David van der Spoel <david.vanderspoel@icm.uu.se>
 */

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <map>
#include <random>

#include "gromacs/commandline/pargs.h"
#include "gromacs/commandline/viewit.h"
#include "gromacs/fileio/xvgr.h"
#include "gromacs/hardware/detecthardware.h"
#include "gromacs/listed-forces/bonded.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/statistics/statistics.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/utility/arraysize.h"
#include "gromacs/utility/coolstuff.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/utility/smalloc.h"
#include "gromacs/utility/strconvert.h"

#include "alex_modules.h"
#include "communication.h"
#include "gentop_core.h"
#include "gmx_simple_comm.h"
#include "molgen.h"
#include "molprop.h"
#include "molprop_util.h"
#include "molprop_xml.h"
#include "molselect.h"
#include "mymol.h"
#include "optparam.h"
#include "poldata.h"
#include "poldata_xml.h"
#include "regression.h"
#include "tune_fc_utils.h"

/*! \brief Write a csv file containing molecule names and bond energy
 *
 * Writes the whole bond energy matrix.
 */
static void dump_csv(const std::vector<std::string>        &ctest,
                     const std::vector<alexandria::MyMol>  &mm,
                     const std::vector<int>                &ntest,
                     const std::vector<double>             &Edissoc,
                     const MatrixWrapper                   &a,
                     const double                           x[])
{
    FILE *csv = gmx_ffopen("tune_fc.csv", "w");
    fprintf(csv, ",");
    for (auto j : ctest)
    {
        fprintf(csv, "%s,", j.c_str());
    }
    fprintf(csv, "\n");
    int i = 0;
    for (auto &mymol : mm)
    {
        fprintf(csv, "%s,", mymol.molProp()->getMolname().c_str());
        for (size_t j = 0; (j < ctest.size()); j++)
        {
            fprintf(csv, "%g,", a.get(j, i));
        }
        fprintf(csv, "%.3f\n", x[i]);
        i++;
    }
    fprintf(csv, "Total,");
    for (auto j : ntest)
    {
        fprintf(csv, "%d,", j);
    }
    fprintf(csv, "\n");
    fprintf(csv, "Edissoc,");
    for (auto j : Edissoc)
    {
        fprintf(csv, "%.3f,", j);
    }
    fprintf(csv, "\n");
    fclose(csv);
}

namespace alexandria
{

typedef struct MolEnergyEntry
{
    // Job type
    int    jobType;
    // Reference energy
    double EReference;
    // Alexandria energy
    double EAlexandria;
} MolEnergyEntry;

class MolEnergy
{
    public:
        MolEnergy() : Force2_(0.0) {}

        void setForce2(double Force2) { Force2_ = Force2; }

        double force2() const { return Force2_; }

        void setTerms(real term[F_NRE]);

        real term(int ftype) const { return term_[ftype]; }

        void addE(int jobType, double Ereference, double Ealexandria);

        void clear();

        const std::vector<MolEnergyEntry> &entries() { return MolEnergyEntry_; }
    private:
        // Sum of forces squared at the OPT geometry
        double                      Force2_;
        // Array of energies at all geometries
        std::vector<MolEnergyEntry> MolEnergyEntry_;
        // Array of energy components at the OPT geometry
        real                        term_[F_NRE] = { 0 };
};

void MolEnergy::addE(int jobType, double Ereference, double Ealexandria)
{
    MolEnergyEntry_.push_back({ jobType, Ereference, Ealexandria });
}

void MolEnergy::clear()
{
    Force2_ = 0;
    for (int i = 0; i < F_NRE; i++)
    {
        term_[i] = 0;
    }
    MolEnergyEntry_.clear();
}

void MolEnergy::setTerms(real term[F_NRE])
{
    for (int i = 0; i < F_NRE; i++)
    {
        term_[i] = term[i];
    }
}

class Optimization : public MolGen, Bayes
{
    using param_type = std::vector<double>;

    private:
        std::vector<ForceConstants> ForceConstants_;
        NonBondParams               NonBondParams_;
        int                         iOpt_[eitNR]      = { 0 };
        bool                        optimizeGeometry_ = false;
        std::vector<PoldataUpdate>  poldataUpdates_;
        real                        factor_           = 1;
        bool                        bDissoc_          = true;
        real                        w_dhf_            = 1;
        const char                 *lot_              = nullptr;
        bool                        calcAll_          = false;
        // Map from molid to MolEnergy
        std::map<int, MolEnergy>    MolEnergyMap_;
    public:


        /*! \brief
         *
         * Constructor
         */
        Optimization()
        {
            iOpt_[eitBONDS] = 1;
        };

        /*! \brief
         *
         * Destructor
         */
        ~Optimization() {};

        /*! \brief
         * Add command line options.
         * \param[inout] pargs Command line options
         */
        void add_pargs(std::vector<t_pargs> *pargs);

        //! \brief To be called after processing options
        void optionsFinished()
        {
            MolGen::optionsFinished();
            setBounds(weight(ermsBOUNDS) > 0);
        }

        int iOpt(int itype) { return iOpt_[itype]; }
    
        void setCalcAll(bool calcAll) { calcAll_ = calcAll; }

        //! \brief Return the level of theory used as reference
        const char *lot() const { return lot_; }

        /*! \brief
         *
         * Check whether molecules are supported by the force field.
         * Check whether all the atomtypes in the molecules are present
         * in the force field file. If not the molecules are removed.
         *
         * Routine also divides molecules over processors.
         */
        void checkSupport(FILE *fp);

        /*! \brief
         *
         * Fill parameter vector using ForceConstants which
         * is built based on Poldata.
         * \param[in] factor Scaling factor for parameters
         */
        void polData2TuneFc(real factor);

        /*! \brief
         *
         * Copy the optimized parameters back to Poldata
         *
         * \param[in] changed Boolean list stating whether a parameter has changed
         */
        void toPolData(const std::vector<bool> &changed);

        /*! \brief
         * Broadcast changes in Poldata to the
         * slaves when in parallel.
         */
        void broadcastPoldataUpdate();

        /*! \brief
         *
         * Compute the dissociation energies for all the bonds.
         * Given all the bonds and the enthalpies of formation of all
         * molecules, we can approximate the dissociation enthalpy (D0 in the
         * Morse potential by least squares fitting the D0 to reproduce the
         * molecular energy (Delta H formation of molecule - Delta H formation of
         * the atoms). This is a crude approximation since all other energy
         * terms in the force field are ignored, however the dissociation
         * energy is the largest contribution to the molecular energy.
         */
        void getDissociationEnergy(FILE *fplog);

        /*! \brief
         * Initialize the optimization algorithm.
         * \param[in] fplog Log file for dumping information.
         */
        void InitOpt(FILE *fplog);

        /*! \brief
         * Do the actual optimization.
         * \param[in] fplog  FILE pointer for logging
         * \param[in] oenv   Output environment for managing xvg files etc.
         * \param[in] xvgconv Output file monitoring parameters
         * \param[in] xvgepot Output file monitoring penalty function
         * \return true if better parameters were found.
         */
        bool optRun(FILE                   *fplog,
                    const gmx_output_env_t *oenv,
                    const char             *xvgconv,
                    const char             *xvgepot);

        /*! \brief
         * Print the results of the optimization per molecule
         * as well as statistics over the minimum energy
         * structures.
         * \param[in] fp     File pointer to write to
         * \param[in] title  Text string to write
         * \param[in] xvg    Filename for correlation plot for DHformation at
         *                   minimum energy structures.
         * \param[in] HF_xvg Filename for energy correlation per compound
         * \param[in] oenv   GROMACS file management structure.
         */
        void printResults(FILE                   *fp,
                          const char             *title,
                          const char             *xvg,
                          const char             *HF_xvg,
                          const gmx_output_env_t *oenv);

        /*! \brief
         * Print data for all the molecules.
         * \param[in] fp     File pointer to write to
         * \param[in] bForce Print forces and energy terms at optimal geometry
         * \param[in] bMtop  Print the molecular topolog7
         */
        void printMolecules(FILE *fp,
                            bool  bForce,
                            bool  bMtop);
        /*! \brief
         * Compute the deviation from the reference energies etc.
         */
        virtual double calcDeviation();

        /*! \brief
         * Send the information from this processor to destination
         * \param[in] cr   Communication data structure
         * \param[in] dest Destination processor
         */
        CommunicationStatus Send(t_commrec *cr, int dest);

        /*! \brief
         * Receive information from another processor to this
         * \param[in] cr   Communication data structure
         * \param[in] src Source processor
         */
        CommunicationStatus Receive(t_commrec *cr, int src);

        //! Distribute all information over all processors
        void broadcast();
};

CommunicationStatus Optimization::Send(t_commrec *cr, int dest)
{
    CommunicationStatus cs;
    cs = gmx_send_data(cr, dest);
    if (CS_OK == cs)
    {
        gmx_send_int(cr, dest, ForceConstants_.size());
        for (auto &fc : ForceConstants_)
        {
            fc.Send(cr, dest);
        }
        NonBondParams_.Send(cr, dest);
    }
    return cs;
}

CommunicationStatus Optimization::Receive(t_commrec *cr, int src)
{
    CommunicationStatus cs;
    cs = gmx_recv_data(cr, src);
    if (CS_OK == cs)
    {
        int nfc = gmx_recv_int(cr, src);

        for (int n = 0; (CS_OK == cs) && (n < nfc); n++)
        {
            alexandria::ForceConstants fc;
            cs = fc.Receive(cr, src);
            if (CS_OK == cs)
            {
                ForceConstants_.push_back(fc);
            }
        }
        cs = NonBondParams_.Receive(cr, src);
    }
    return cs;
}

void Optimization::add_pargs(std::vector<t_pargs> *pargs)
{
    t_pargs pa[] =
    {
        { "-dissoc",  FALSE, etBOOL, {&bDissoc_},
          "Derive dissociation energy from the enthalpy of formation. If not chosen, the dissociation energy will be read from the gentop.dat file." },
        { "-opt_geom", FALSE, etBOOL, {&optimizeGeometry_},
          "Optimize bond lengths, angles and dihedral angles" },
        { "-weight_dhf", FALSE, etREAL, {&w_dhf_},
          "Fitting weight of the minimum energy structure, representing the enthalpy of formation relative to high energy structures." },
        { "-angles",  FALSE, etINT, {&(iOpt_[eitANGLES])},
          "Optimize angle parameters" },
        { "-langles", FALSE, etINT, {&(iOpt_[eitLINEAR_ANGLES])},
          "Optimize linear angle parameters" },
        { "-dihedrals", FALSE, etINT, {&(iOpt_[eitPROPER_DIHEDRALS])},
          "Optimize proper dihedral parameters" },
        { "-impropers", FALSE, etINT, {&(iOpt_[eitIMPROPER_DIHEDRALS])},
          "Optimize improper dihedral parameters" },
        { "-vdw", FALSE, etINT, {&(iOpt_[eitVDW])},
          "Optimize van der Waals parameters" },
        { "-factor", FALSE, etREAL, {&factor_},
          "Parameters will be taken within the limit factor*x - x/factor." }
    };
    for (size_t i = 0; i < sizeof(pa)/sizeof(pa[0]); i++)
    {
        pargs->push_back(pa[i]);
    }
    addOptions(pargs, etuneFC);
    Bayes::add_pargs(pargs);
}

void Optimization::broadcastPoldataUpdate()
{
    t_commrec *cr = commrec();
    if (!PAR(cr))
    {
        return;
    }

    if (MASTER(cr))
    {
        for (int i = 1; i < cr->nnodes; i++)
        {
            gmx_send_int(cr, i, static_cast<int>(poldataUpdates_.size()));
            for (auto &p : poldataUpdates_)
            {
                p.Send(cr, i);
            }
        }
    }
    else
    {
        int n = gmx_recv_int(cr, 0);
        poldataUpdates_.resize(n);
        for (int i = 0; i < n; i++)
        {
            poldataUpdates_[i].Receive(cr, 0);
        }
        if (debug)
        {
            fprintf(debug, "broadcastPoldataUpdate: received %d updates\n",
                    n);
        }
    }
}

void Optimization::broadcast()
{
    t_commrec *cr = commrec();
    if (!PAR(cr))
    {
        return;
    }
    const int src = 0;
    if (MASTER(cr))
    {
        for (int dest = 1; dest < cr->nnodes; dest++)
        {
            Send(cr, dest);
        }
    }
    else
    {
        if (nullptr != debug)
        {
            fprintf(debug, "Updating ForceConstants and NonBondParams on node %d\n",
                    cr->nodeid);
        }
        Receive(cr, src);
    }
    poldata()->broadcast(cr);
}

void Optimization::checkSupport(FILE *fp)
{
    auto ntotal  = mymols().size();
    auto nlocal  = 0;

    for (auto mymol =  mymols().begin(); mymol <  mymols().end(); )
    {
        if (mymol->eSupp_ != eSupportLocal)
        {
            mymol++;
            continue;
        }

        bool bSupport = true;
        for (int bt = 0; bSupport && (bt < eitNR); bt++)
        {
            if (iOpt_[bt])
            {
                auto iType = static_cast<InteractionType>(bt);
                auto fs    = poldata()->findForces(iType);
                if (fs == poldata()->forcesEnd())
                {
                    continue;
                }
                int ft     = fs->fType();
                bSupport   = (mymol->ltop_ != nullptr);
                for (int i = 0; (i < mymol->ltop_->idef.il[ft].nr) && bSupport;
                     i += interaction_function[ft].nratoms+1)
                {
                    // Loop starts from 1 because the first value is the function type
                    std::vector<std::string> atoms;
                    for(int j = 1; j <= interaction_function[ft].nratoms && bSupport; j++)
                    {
                        std::string aa;
                        int         ai = mymol->ltop_->idef.il[ft].iatoms[i+j];
                        if (!poldata()->atypeToBtype(*mymol->atoms_->atomtype[ai], aa))
                        {
                            bSupport = false;
                        }
                        else
                        {
                            atoms.push_back(aa);
                        }
                    }
                    if (bSupport)
                    {
                        bSupport = fs->findForce(atoms) != fs->forceEnd();
                    }
                }
            }
        }
        if (!bSupport)
        {
            fprintf(stderr, "No force field support for %s\n",
                    mymol->molProp()->getMolname().c_str());
            mymol =  mymols().erase(mymol);
        }
        else
        {
            mymol++;
            nlocal++;
        }
    }
    if (PAR(commrec()))
    {
        gmx_sumi(1, &nlocal, commrec());
    }
    if (nullptr != fp)
    {
        fprintf(fp, "%d out of %zu molecules have support in the force field.\n",
                nlocal, ntotal);
    }
}

void Optimization::polData2TuneFc(real factor)
{
    for (auto &fc : ForceConstants_)
    {
        for (auto b = fc.beginBN(); b  < fc.endBN(); ++b)
        {
            if (optimizeGeometry_)
            {
                Bayes::addParam(b->geometry(), factor);
            }
            for (const auto &p : b->paramValues())
            {
                if (fc.ftype() == F_FOURDIHS)
                {
                    Bayes::addParam(p, -10, 10);
                }
                else
                {
                    Bayes::addParam(p, factor);
                }
            }
        }
    }
    for (auto at = NonBondParams_.beginAT(); at < NonBondParams_.endAT(); ++at)
    {
        for (const auto &p : at->paramValues())
        {
            Bayes::addParam(p, factor);
        }
    }
    if (debug)
    {
        fprintf(debug, "poldata2TuneFc: Copied %zu parameters.\n",
                Bayes::nParam());
    }
}

void Optimization::toPolData(const std::vector<bool> &changed)
{
    size_t n   = 0;
    poldataUpdates_.clear();
    auto param = Bayes::getParam();
    for (auto &fc : ForceConstants_)
    {
        const auto iType = fc.interactionType();
        for (auto b = fc.beginBN(); b  < fc.endBN(); ++b)
        {
            std::string paramString;
            // Check for geometry changes
            bool        bondChanged = false;
            double      geometry    = 0;
            if (optimizeGeometry_)
            {
                bondChanged = changed[n];
                geometry    = param[n++]; 
            }
            for (size_t p = 0; p < b->nParams(); p++)
            {
                bondChanged = bondChanged || changed[n];
                paramString.append(gmx::formatString(" %g", param[n++]));
            }
            if (bondChanged)
            {
                b->setParamString(paramString);
                poldataUpdates_.push_back(PoldataUpdate(iType, 
                                                        b->poldataIndex(),
                                                        geometry,
                                                        paramString));
            }
        }
    }
    for (auto at = NonBondParams_.beginAT(); at < NonBondParams_.endAT(); ++at)
    {
        bool        bondChanged = false;
        double      geometry    = 0;
        std::string paramString;
        for (size_t p = 0; p < at->nParams(); p++)
        {
            bondChanged = bondChanged || changed[n];
            paramString.append(gmx::formatString(" %g", param[n++]));
        }
        if (bondChanged)
        {
            at->setParamString(paramString);
            poldataUpdates_.push_back(PoldataUpdate(eitVDW,
                                                    at->poldataIndex(),
                                                    geometry,
                                                    paramString));
        }
    }
    GMX_RELEASE_ASSERT(n == param.size(), "Number of parameters set should be equal to the length of the parameter array");
}

void Optimization::getDissociationEnergy(FILE *fplog)
{
    std::vector<double>         rhs;
    std::vector<int>            ntest;
    std::vector<std::string>    ctest;

    int nD   = ForceConstants_[eitBONDS].nbad();
    int nMol = mymols().size();

    if ((0 == nD) || (0 == nMol))
    {
        gmx_fatal(FARGS, "Number of variables is %d and number of molecules is %d",
                  nD, nMol);
    }

    MatrixWrapper a(nD, nMol);
    MatrixWrapper a_copy(nD, nMol);
    ntest.resize(nD, 0);
    ctest.resize(nD);

    fprintf(fplog, "There are %d different bondtypes to optimize the heat of formation\n", nD);
    fprintf(fplog, "There are %d (experimental) reference heat of formation.\n", nMol);

    auto fs  = poldata()->findForces(eitBONDS);
    auto ftb = fs->fType();
    auto j   = 0;

    for (auto mymol =  mymols().begin(); mymol <  mymols().end(); mymol++, j++)
    {
        for (auto i = 0; i < mymol->ltop_->idef.il[ftb].nr;
             i += interaction_function[ftb].nratoms+1)
        {
            auto                     ai = mymol->ltop_->idef.il[ftb].iatoms[i+1];
            auto                     aj = mymol->ltop_->idef.il[ftb].iatoms[i+2];
            std::string              aai, aaj;
            std::vector<std::string> atoms;
            if (poldata()->atypeToBtype(*mymol->atoms_->atomtype[ai], aai) &&
                poldata()->atypeToBtype(*mymol->atoms_->atomtype[aj], aaj))
            {
                atoms  = {aai, aaj};
                auto f = fs->findForce(atoms);
                if (fs->forceEnd() != f)
                {
                    auto  gt  = f - fs->forceBegin();
                    auto  gti = ForceConstants_[eitBONDS].reverseIndex(gt);
                    a.set(gti, j, a.get(gti, j) + 1);
                    a_copy.set(gti, j, a.get(gti, j));
                    ntest[gti]++;
                    if (ctest[gti].empty())
                    {
                        char buf[STRLEN];
                        snprintf(buf, sizeof(buf), "%s-%s", aai.c_str(), aaj.c_str());
                        ctest[gti].assign(buf);
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "No parameters for bond %s-%s in the force field, atoms %s-%s mol %s",
                          aai.c_str(), aaj.c_str(),
                          *mymol->atoms_->atomtype[ai],
                          *mymol->atoms_->atomtype[aj],
                          mymol->molProp()->getIupac().c_str());
            }
        }
        rhs.push_back(-mymol->Emol_);
    }

    char buf[STRLEN];
    snprintf(buf, sizeof(buf), "Inconsistency in number of energies nMol %d != #rhs %zu", nMol, rhs.size());
    GMX_RELEASE_ASSERT(static_cast<int>(rhs.size()) == nMol, buf);

    auto nzero = std::count_if(ntest.begin(), ntest.end(), [](const int n)
                               {
                                   return n == 0;
                               });

    GMX_RELEASE_ASSERT(nzero == 0, "Inconsistency in the number of bonds in poldata and ForceConstants_");

    std::vector<double> Edissoc(nD);
    a.solve(rhs, &Edissoc);
    if (debug)
    {
        dump_csv(ctest,  mymols(), ntest, Edissoc, a_copy, rhs.data());
    }
    for (size_t i = 0; i < ctest.size(); i++)
    {
        if (fplog)
        {
            fprintf(fplog, "Optimized dissociation energy for %8s with %4d copies to %g\n",
                    ctest[i].c_str(), ntest[i], Edissoc[i]);
        }
    }

    auto i = 0;
    for (auto b = ForceConstants_[eitBONDS].beginBN();
         b < ForceConstants_[eitBONDS].endBN(); ++b)
    {
        const auto       atoms = gmx::splitString(b->name());
        auto             fs    = poldata()->findForces(eitBONDS);
        auto             f     = fs->findForce(atoms);
        GMX_RELEASE_ASSERT(fs->forceEnd() != f, "Cannot find my bonds");
        const auto       pp    = gmx::splitString(b->paramString());
        char             buf[256];

        if (!f->fixed())
        {
            /* TODO: De is assumed to be the first parameter. Please fix. */
            snprintf(buf, sizeof(buf), "%.2f  %s", std::max(100.0, Edissoc[i]), pp[1].c_str());
            f->setParams(buf);
            b->setParamString(buf);
            i++;
        }
    }
}

void Optimization::InitOpt(FILE *fplog)
{
    for (auto fs = poldata()->forcesBegin();
         fs != poldata()->forcesEnd(); ++fs)
    {
        int bt = static_cast<int>(fs->iType());
        if (iOpt_[bt])
        {
            ForceConstants fc(fs->fType(), fs->iType(), iOpt_[bt]);
            fc.analyzeIdef(mymols(), poldata(), optimizeGeometry_);
            fc.makeReverseIndex();
            fc.dump(fplog);
            ForceConstants_.push_back(std::move(fc));
        }
    }
    if (bDissoc_)
    {
        if (ForceConstants_[eitBONDS].nbad() <= mymols().size())
        {
            getDissociationEnergy(fplog);
        }
        else
        {
            printf("\n"
                   "WARNING: %zu molecule(s) is (are) not enough to calculate dissociation\n"
                   "         energy for %zu bond type(s) using linear regression. Default\n"
                   "         values from gentop.dat will be used as the initial guess.\n"
                   "         Recomendation is to add more molecules having the same bond types.\n\n",
                   mymols().size(), ForceConstants_[eitBONDS].nbad());
        }
    }
    
    if (iOpt_[eitVDW])
    {
        NonBondParams_.setOpt(true);
        NonBondParams_.analyzeIdef(mymols(), poldata());
        NonBondParams_.makeReverseIndex();
        NonBondParams_.dump(fplog);
    }
    polData2TuneFc(factor_);
}

double Optimization::calcDeviation()
{
    if (!calcAll_)
    {
        if (PAR(commrec()))
        {
            broadcastPoldataUpdate();
        }
    }
    for (auto &pUpd : poldataUpdates_)
    {
        pUpd.execute(poldata());
    }
    poldataUpdates_.clear();
    if (calcAll_)
    {
        Bayes::dumpParam(debug);
    }
    for (auto &mem : MolEnergyMap_)
    {
        mem.second.clear();
    }
    // First compute all the energies and store them
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupportLocal) ||
            (calcAll_ && (mymol.eSupp_ == eSupportRemote)))
        {
            int      nSP    = 0, nOpt = 0;
            int      natoms = mymol.atoms_->nr;
            gmx_bool bpolar = (mymol.shellfc_ != nullptr);
            double   optHF;
            if (mymol.molProp()->getOptHF(&optHF))
            {
                for (const auto fc : ForceConstants_)
                {
                    if (fc.nbad() > 0)
                    {
                        mymol.UpdateIdef(poldata(), fc.interactionType());
                    }
                }

                if (NonBondParams_.nAT() > 0)
                {
                    mymol.UpdateIdef(poldata(),
                                     NonBondParams_.interactionType());
                }
                mymol.f_.resizeWithPadding(natoms);
                mymol.optf_.resizeWithPadding(natoms);
                int  molid          = mymol.molProp()->getIndex();
                auto molEnergyEntry = MolEnergyMap_.find(molid);
                if (molEnergyEntry == MolEnergyMap_.end())
                {
                    // Create new (empty) MolEnergy and insert it into the map
                    MolEnergy me;
                    molEnergyEntry = 
                        MolEnergyMap_.insert(MolEnergyMap_.begin(),
                                             std::pair<int, MolEnergy>(molid, std::move(me)));
                }
                molEnergyEntry->second.clear();
                // Now loop over experiments!
                for (auto ei = mymol.molProp()->BeginExperiment();
                     ei < mymol.molProp()->EndExperiment(); ++ei)
                {
                    auto jtype = ei->getJobtype();
                    // Exclude other experimental data points!
                    if (jtype == JOB_OPT || jtype == JOB_SP)
                    {
                        double spHF;
                        if (!ei->getHF(&spHF))
                        {
                            continue;
                        }

                        FILE *dbcopy = debug;
                        debug  = nullptr;
                        mymol.changeCoordinate(ei, bpolar);
                        mymol.computeForces(debug, commrec());
                        debug  = dbcopy;

                        double deltaEref  = spHF - optHF;
                        if (deltaEref < 0)
                        {
                            gmx_fatal(FARGS, "Energy mismatch for %s. spHF = %g (%s) optHF = %g",
                                      mymol.molProp()->getMolname().c_str(),
                                      spHF,
                                      ei->getDatafile().c_str(),
                                      optHF);
                        }
                        double deltaEalex = mymol.enerd_->term[F_EPOT] - mymol.Emol_;
                        molEnergyEntry->second.addE(jtype, deltaEref, deltaEalex);

                        // Force is added only for the opt geometry
                        double OptForce2 = 0.0;
                        if (jtype == JOB_OPT)
                        {
                            for (int j = 0; j < natoms; j++)
                            {
                                OptForce2 += iprod(mymol.f_[j], mymol.f_[j]);
                                copy_rvec(mymol.f_[j], mymol.optf_[j]);
                            }
                            OptForce2 /= natoms;
                            molEnergyEntry->second.setForce2(OptForce2);
                            molEnergyEntry->second.setTerms(mymol.enerd_->term);
                            nOpt += 1;
                        }
                        else
                        {
                            nSP += 1;
                        }
                        if (nullptr != debug)
                        {
                            int angleType = poldata()->findForces(eitANGLES)->fType();
                            int pdihType  = poldata()->findForces(eitPROPER_DIHEDRALS)->fType();
                            int idihType  = poldata()->findForces(eitIMPROPER_DIHEDRALS)->fType();
                            int vdwType   = poldata()->getVdwFtype();
                            fprintf(debug, "spHF: %g  optHF: %g  deltaRef: %g  deltaEalex: %g\n",
                                    spHF, optHF, deltaEref, deltaEalex);
                            fprintf(debug, "%s Chi2 %g Morse %g  "
                                    "%s %g Langle %g %s %g %s %g Coul %g VdW %g POL %g  Force %g\n",
                                    mymol.molProp()->getMolname().c_str(),
                                    gmx::square(deltaEalex-deltaEref),
                                    mymol.enerd_->term[F_MORSE],
                                    interaction_function[angleType].name,
                                    mymol.enerd_->term[angleType],
                                    mymol.enerd_->term[F_LINEAR_ANGLES],
                                    interaction_function[pdihType].name,
                                    mymol.enerd_->term[pdihType],
                                    interaction_function[idihType].name,
                                    mymol.enerd_->term[idihType],
                                    mymol.enerd_->term[F_COUL_SR],
                                    mymol.enerd_->term[vdwType],
                                    mymol.enerd_->term[F_POLARIZATION],
                                    std::sqrt(OptForce2));
                        }
                    }
                }
            }
            else
            {
                gmx_fatal(FARGS, "There is no optimized structure for %s\n",
                          mymol.molProp()->getMolname().c_str());
            }
            if (debug)
            {
                fprintf(debug, "%d: %s nOpt %d nSP %d\n",
                        commrec()->nodeid,
                        mymol.molProp()->getMolname().c_str(),
                        nOpt, nSP);
            }
        }
    }
    // Now compute the deviation for the fitting or otherwise
    resetEnergies();
    double nCalc = 0;
    double ePot2 = 0;
    for (auto &mymol : mymols())
    {
        if ((mymol.eSupp_ == eSupportLocal) ||
            (calcAll_ && (mymol.eSupp_ == eSupportRemote)))
        {
            int  molid          = mymol.molProp()->getIndex();
            auto molEnergyEntry = MolEnergyMap_.find(molid);
            if (molEnergyEntry != MolEnergyMap_.end())
            {
                increaseEnergy(ermsForce2, molEnergyEntry->second.force2());
                if (debug)
                {
                    fprintf(debug, "%d: %s molid: %d nEntries: %zu\n",
                            commrec()->nodeid,
                            mymol.molProp()->getMolname().c_str(),
                            molid,
                            molEnergyEntry->second.entries().size());
                }
                for (const auto &d : molEnergyEntry->second.entries())
                {
                    ePot2 += gmx::square(d.EReference-d.EAlexandria);
                    nCalc += 1;
                }
            }
        }
    }
    if (debug)
    {
        fprintf(debug, "%d: ePot2 = %g nCalc = %g\n",
                commrec()->nodeid, ePot2, nCalc);
    }
    if (!calcAll_ && PAR(commrec()))
    {
        gmx_sumd(1, &ePot2, commrec());
        gmx_sumd(1, &nCalc, commrec());
    }
    double chi2 = ePot2/nCalc;
    if (debug)
    {
        fprintf(debug, "%d: ePot2 = %g nCalc = %g chi2 = %g\n",
                commrec()->nodeid, ePot2, nCalc, chi2);
    }
    setEnergy(ermsEPOT, chi2);
    setEnergy(ermsTOT, chi2);
    printEnergies(debug);
    
    return energy(ermsTOT);
}

bool Optimization::optRun(FILE                   *fplog,
                          const gmx_output_env_t *oenv,
                          const char             *xvgconv,
                          const char             *xvgepot)
{
    bool bMinimum = false;
    if (MASTER(commrec()))
    {
        if (PAR(commrec()))
        {
            // Tell the slave nodes how many times they have
            // to run calcDeviation.
            int niter = 1+Bayes::numberObjectiveFunctionCalls();
            for (int dest = 1; dest < commrec()->nnodes; dest++)
            {
                gmx_send_int(commrec(), dest, niter);
            }
        }
        double chi2_min = Bayes::objFunction(Bayes::getParam().data());
        double chi2     = chi2_min;
        if (fplog)
        {
            fprintf(fplog, "Initial chi2 %g\n", chi2_min);
        }
        {
            auto func = [&] (const double v[]) {
                            return Bayes::objFunction(v);
                        };
            Bayes::setFunc(func, &chi2);
            Bayes::setOutputFiles(xvgconv, xvgepot, oenv);
        }
        Bayes::MCMC();
        if (fplog)
        {
            fprintf(fplog, "Final chi2 %g\n", chi2_min);
        }
        if (chi2 < chi2_min)
        {
            chi2_min = chi2;
            bMinimum = true;

            setCalcAll(true);
            if (fplog)
            {
                auto pmean  = Bayes::getPmean();
                auto psigma = Bayes::getPsigma();
                auto best   = Bayes::getBestParam();
                // This call copies data to poldata as well.
                double chi2 = Bayes::objFunction(best.data());
                fprintf(fplog, "\nLowest RMSD value during optimization: %g.\n",
                        std::sqrt(chi2));
                fprintf(fplog, "Parameters after the optimization:\n");
                fprintf(fplog, "%-5s  %10s  %10s  %10s\n", "Index",
                        "Average", "Std. Dev.", "Optimum");
                for (size_t k = 0; k < Bayes::nParam(); k++)
                {
                    fprintf(fplog, "%5zu  %10g  %10g  %10g\n",
                            k, pmean[k], psigma[k], best[k]);
                }
            }
            printEnergies(fplog);
        }
    }
    else
    {
        /* S L A V E   N O D E S */
        int niter = gmx_recv_int(commrec(), 0);
        for (int n = 0; n < niter; n++)
        {
            calcAll_ = false;
            (void) calcDeviation();
        }
    }
    return bMinimum;
}

void Optimization::printMolecules(FILE *fp,
                                  bool  bForce,
                                  bool  bMtop)
{
    int j, k;
    for (const auto &mi : mymols())
    {
        int nSP = 0, nOpt = 0;
        for (auto ei = mi.molProp()->BeginExperiment(); ei < mi.molProp()->EndExperiment(); ++ei)
        {
            auto jtype = ei->getJobtype();
            if (jtype == JOB_SP)
            {
                nSP += 1;
            }
            else if (jtype == JOB_OPT)
            {
                nOpt += 1;
            }
        }
        if (nOpt != 1)
        {
            fprintf(stderr, "Number of optimized conformations is %d. Check your input.\n", nOpt);
        }
        fprintf(fp, "%s natoms: %d Opt conformations: %d SP conformations: %d\n",
                mi.molProp()->getMolname().c_str(),
                mi.atoms_->nr,
                nOpt,
                nSP);
        for (j = 0; j < mi.atoms_->nr; j++)
        {
            fprintf(fp, "  %-5s  %-5s  q = %10g",
                    *(mi.atoms_->atomname[j]),
                    *(mi.atoms_->atomtype[j]),
                    mi.atoms_->atom[j].q);
            if (bForce)
            {
                fprintf(fp, "   f = %8.3f  %8.3f  %8.3f",
                        mi.optf_[j][XX], mi.optf_[j][YY], mi.optf_[j][ZZ]);
            }
            fprintf(fp, "\n");
        }

        if (bMtop)
        {
            pr_mtop(fp, 0, mi.molProp()->getMolname().c_str(), mi.mtop_, true, false);
        }
    }
    if (bForce)
    {
        for (const auto &mi : mymols())
        {
            int  molid     = mi.molProp()->getIndex();
            auto molEnergy = MolEnergyMap_.find(molid);
            if (molEnergy != MolEnergyMap_.end())
            {
                fprintf(fp, "%-20s", mi.molProp()->getMolname().c_str());
                for (k = 0; k < F_NRE; k++)
                {
                    real term = molEnergy->second.term(k);
                    if (term != 0 ||
                        (mi.mtop_->moltype[0].ilist[k].size() > 0))
                    {
                        fprintf(fp, " %s: %.2f", interaction_function[k].name,
                                mi.enerd_->term[k]);
                    }
                }
                fprintf(fp, "\n");
            }
        }
    }
}

void Optimization::printResults(FILE                   *fp,
                                const char             *title,
                                const char             *hform_xvg,
                                const char             *HF_xvg,
                                const gmx_output_env_t *oenv)
{
    FILE       *xfp = nullptr, *hfp = nullptr;
    gmx_stats_t gst;

    gst = gmx_stats_init();
    {
        const char *title = "Enthalpy of Formation";
        const char *yaxis = "Alexandria (kJ/mol)";
        if (nullptr != hform_xvg)
        {
            xfp = xvgropen(hform_xvg, title, "Experiment (kJ/mol)", yaxis, oenv);
        }
        if (nullptr != HF_xvg)
        {
            hfp = xvgropen(HF_xvg, title,
                           "Experiment + \\f{12}DD\\f{4}E(B3LYP/aug-cc-pVTZ) (kJ/mol)",
                           yaxis, oenv);
            xvgr_view(hfp, 0.15, 0.15, 0.75, 0.85, oenv);
        }
    }
    fprintf(fp, "%s\n", title);
    fprintf(fp, "Fit of energy at different conformations to y = ax\n");
    fprintf(fp, "Nr.   %-30s %10s %10s %7s %7s %7s %7s %7s %4s\n",
            "Molecule", "DHf@298K", "Emol@0K", "rmsF@0K",
            "rms E", "MSE E", "a", "r2", "N");

    int    imol          = 0;
    int    nconformation = 0;
    for (auto mi = mymols().begin(); mi < mymols().end(); mi++, imol++)
    {
        int  molid      = mi->molProp()->getIndex();
        auto molEnergy  = MolEnergyMap_.find(molid);
        if (molEnergy != MolEnergyMap_.end())
        {
            if (nullptr != hfp)
            {
                fprintf(hfp, "@ s%d legend \"%s\"\n", imol,
                        mi->molProp()->getMolname().c_str());
                fprintf(hfp, "@type xy\n");
            }
            gmx_stats_t gmol = gmx_stats_init();
            for (const auto &entry : molEnergy->second.entries())
            {
                real deltaE = entry.EAlexandria - entry.EReference;
                if (nullptr != xfp && entry.jobType == JOB_OPT)
                {
                    fprintf(xfp, "%10g  %10g\n", mi->Hform_, mi->Hform_ + deltaE);
                }
                real Hexper = mi->Hform_ + entry.EReference;
                real Halex  = mi->Hform_ + entry.EAlexandria;
                if (nullptr != hfp)
                {
                    fprintf(hfp, "%10g  %10g\n", Hexper, Halex);
                }
                gmx_stats_add_point(gmol, Hexper, Halex, 0, 0);
                if (entry.jobType == JOB_OPT)
                {
                    gmx_stats_add_point(gst, Hexper, Halex, 0, 0);
                }
            }
            if (nullptr != hfp)
            {
                fprintf(hfp, "&\n");
            }

            real a, chi2, da, rmsd, Rfit, Rdata, mse, mae;
            int  N;
            // Note we are ignoring the return value for these functions
            (void) gmx_stats_get_a(gmol, 0, &a, &da, &chi2, &Rfit);
            (void) gmx_stats_get_rmsd(gmol, &rmsd);
            (void) gmx_stats_get_mse_mae(gmol, &mse, &mae);
            (void) gmx_stats_get_npoints(gmol, &N);
            nconformation += N;
            gmx_stats_get_corr_coeff(gmol, &Rdata);
            fprintf(fp, "%-5d %-30s %10g %10g %7.1f %7.3f %7.3f %7.3f %6.1f%% %4d\n",
                    imol,
                    mi->molProp()->getMolname().c_str(),
                    mi->Hform_,
                    mi->Emol_,
                    std::sqrt(molEnergy->second.force2()), rmsd,
                    mse, a, Rdata*100, N);
            gmx_stats_free(gmol);
        }
    }
    fprintf(fp, "RMSD from target energies for %zu compounds and %d conformation is %g.\n",
            mymols().size(), nconformation, std::sqrt(energy(ermsTOT)));
    fprintf(fp, "\n");
    if (nullptr != hform_xvg)
    {
        xvgrclose(xfp);
        do_view(oenv, hform_xvg, nullptr);
    }
    if (nullptr != HF_xvg)
    {
        xvgrclose(hfp);
        do_view(oenv, HF_xvg, nullptr);
    }

    real a, da, chi2, rmsd, Rfit, Rdata;
    int  N;
    gmx_stats_get_a(gst, 1, &a, &da, &chi2, &Rfit);
    gmx_stats_get_npoints(gst, &N);
    gmx_stats_get_corr_coeff(gst, &Rdata);
    gmx_stats_get_rmsd(gst, &rmsd);
    fprintf(fp, "Regression analysis fit to y = ax of %d minimum energy structures:\n", N);
    fprintf(fp, "a = %.3f  r2 = %.1f%%  rmsd = %.2f kJ/mol\n",
            a, Rdata*100, rmsd);
    gmx_stats_free(gst);
    fflush(fp);
}

} // namespace alexandria

int alex_tune_fc(int argc, char *argv[])
{
    const char           *desc[] = {
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
        "A selection of molecules into a training set and a test set (or ignore set)",
        "can be made using option [TT]-sel[tt]. The format of this file is:[BR]",
        "iupac|Train[BR]",
        "iupac|Test[BR]",
        "iupac|Ignore[BR]",
        "and you should ideally have a line for each molecule in the molecule database",
        "([TT]-f[tt] option). Missing molecules will be ignored."
    };

    t_filenm              fnm[] = {
        { efDAT, "-f",     "allmols",    ffREAD  },
        { efDAT, "-d",     "gentop",     ffOPTRD },
        { efDAT, "-o",     "tune_fc",    ffWRITE },
        { efDAT, "-sel",   "molselect",  ffREAD  },
        { efXVG, "-table", "table",      ffOPTRD },
        { efLOG, "-g",     "tune_fc",    ffWRITE },
        { efXVG, "-x",     "hform-corr", ffWRITE },
        { efXVG, "-hf",    "hf-corr",    ffWRITE },
        { efXVG, "-conv",  "param-conv", ffWRITE },
        { efXVG, "-epot",  "param-epot", ffWRITE }
    };

    const int             NFILE         = asize(fnm);

    int                   reinit        = 0;
    int                   compress      = 0;
    int                   nmultisim     = 0;
    gmx_bool              bZPE          = false;
    gmx_bool              bZero         = true;
    gmx_bool              bTestPar      = false;
    gmx_bool              bForceOutput  = true;
    t_pargs               pa[]          = {
        { "-multi",   FALSE, etINT, {&nmultisim},
          "Do optimization in multiple simulation" },
        { "-reinit",  FALSE, etINT, {&reinit},
          "After this many iterations the search vectors are randomized again. A value of 0 means this is never done at all." },
        { "-zpe",     FALSE, etBOOL, {&bZPE},
          "Consider zero-point energy from thermochemistry calculations in order to calculate the reference enthalpy of the molecule. If set, the zero point energy will be subtracted from the target energy when optimizing the force field model." },
        { "-testpar", FALSE, etBOOL, {&bTestPar},
          "Test the parallel execution gives the same result." },
        { "-compress", FALSE, etBOOL, {&compress},
          "Compress output XML file" },
        { "-force_output", FALSE, etBOOL, {&bForceOutput},
          "Write output even if no new minimum is found" }
    };

    FILE                 *fplog;
    gmx_output_env_t     *oenv;
    time_t                my_t;
    MolSelect             gms;

    std::vector<t_pargs>  pargs;
    for (size_t i = 0; i < asize(pa); i++)
    {
        pargs.push_back(pa[i]);
    }
    alexandria::Optimization opt;
    opt.add_pargs(&pargs);

    if (!parse_common_args(&argc, argv, PCA_CAN_VIEW, NFILE, fnm,
                           pargs.size(), pargs.data(),
                           asize(desc), desc, 0, nullptr, &oenv))
    {
        return 0;
    }
    opt.optionsFinished();

    if (MASTER(opt.commrec()))
    {
        fplog = gmx_ffopen(opt2fn("-g", NFILE, fnm), "w");

        time(&my_t);
        fprintf(fplog, "# This file was created %s", ctime(&my_t));
        fprintf(fplog, "# alexandria is part of G R O M A C S:\n#\n");
        fprintf(fplog, "# %s\n#\n", gmx::bromacs().c_str());
    }
    else
    {
        fplog = nullptr;
    }

    if (MASTER(opt.commrec()))
    {
        gms.read(opt2fn_null("-sel", NFILE, fnm));
    }

    const char *tabfn = opt2fn_null("-table", NFILE, fnm);

    opt.Read(fplog ? fplog : (debug ? debug : nullptr),
             opt2fn("-f", NFILE, fnm),
             opt2fn_null("-d", NFILE, fnm),
             bZero,
             nullptr,
             gms,
             false,
             false,
             opt.iOpt(eitPROPER_DIHEDRALS),
             bZPE,
             false,
             true,
             tabfn);

    opt.checkSupport(fplog);

    if (MASTER(opt.commrec()))
    {
        opt.InitOpt(fplog);
        opt.printMolecules(fplog, false, false);
    }

    opt.broadcast();

    if (bTestPar)
    {
        opt.setCalcAll(false);
        auto chi2 = opt.calcDeviation();
        if (MASTER(opt.commrec()))
        {
            fprintf(fplog, "chi2 = %g\n", chi2);
            opt.printResults(fplog, (char *)"Before optimization - parallel test",
                             nullptr, nullptr, oenv);
        }
    }
    if (MASTER(opt.commrec()))
    {
        opt.setCalcAll(true);
        auto chi2 = opt.calcDeviation();
        fprintf(fplog, "chi2 = %g\n", chi2);
        opt.printResults(fplog, (char *)"Before optimization",
                         nullptr, nullptr, oenv);
        opt.setCalcAll(false);
    }
    if (bTestPar)
    {
        gmx_ffclose(fplog);
        return 0;
    }
    // Optimize the parameters to minimize the penalty function.
    bool bMinimum = opt.optRun(fplog,
                               oenv,
                               opt2fn("-conv", NFILE, fnm),
                               opt2fn("-epot", NFILE, fnm));

    if (MASTER(opt.commrec()))
    {
        if (bMinimum || bForceOutput)
        {
            // Now print the output.
            opt.printMolecules(fplog, true, false);
            opt.printResults(fplog, (char *)"After optimization",
                             opt2fn("-x", NFILE, fnm),
                             opt2fn("-hf", NFILE, fnm),
                             oenv);
            writePoldata(opt2fn("-o", NFILE, fnm), opt.poldata(), compress);
        }
        else if (!bMinimum)
        {
            printf("No improved parameters found. Please try again with more iterations.\n");
        }
        gmx_ffclose(fplog);
    }

    return 0;
}
