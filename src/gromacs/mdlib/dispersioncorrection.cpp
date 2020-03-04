/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019,2020, by the GROMACS development team, led by
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
#include "gmxpre.h"

#include "dispersioncorrection.h"

#include <cstdio>

#include "gromacs/math/units.h"
#include "gromacs/math/utilities.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/forcerec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/interaction_const.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/tables/forcetable.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/logger.h"

/* Implementation here to avoid other files needing to include the file that defines t_nblists */
DispersionCorrection::InteractionParams::~InteractionParams() = default;

/* Returns a matrix, as flat list, of combination rule combined LJ parameters */
static std::vector<real> mk_nbfp_combination_rule(const gmx_ffparams_t& ffparams, const int comb_rule)
{
    const int atnr = ffparams.atnr;

    std::vector<real> nbfp(atnr * atnr * 2);

    for (int i = 0; i < atnr; ++i)
    {
        for (int j = 0; j < atnr; ++j)
        {
            const real c6i  = ffparams.iparams[i * (atnr + 1)].lj.c6;
            const real c12i = ffparams.iparams[i * (atnr + 1)].lj.c12;
            const real c6j  = ffparams.iparams[j * (atnr + 1)].lj.c6;
            const real c12j = ffparams.iparams[j * (atnr + 1)].lj.c12;
            real       c6   = std::sqrt(c6i * c6j);
            real       c12  = std::sqrt(c12i * c12j);
            if (comb_rule == eCOMB_ARITHMETIC && !gmx_numzero(c6) && !gmx_numzero(c12))
            {
                const real sigmai = gmx::sixthroot(c12i / c6i);
                const real sigmaj = gmx::sixthroot(c12j / c6j);
                const real epsi   = c6i * c6i / c12i;
                const real epsj   = c6j * c6j / c12j;
                const real sigma  = 0.5 * (sigmai + sigmaj);
                const real eps    = std::sqrt(epsi * epsj);
                c6                = eps * gmx::power6(sigma);
                c12               = eps * gmx::power12(sigma);
            }
            C6(nbfp, atnr, i, j)  = c6 * 6.0;
            C12(nbfp, atnr, i, j) = c12 * 12.0;
        }
    }

    return nbfp;
}

/* Returns the A-topology atom type when aOrB=0, the B-topology atom type when aOrB=1 */
static int atomtypeAOrB(const t_atom& atom, int aOrB)
{
    if (aOrB == 0)
    {
        return atom.type;
    }
    else
    {
        return atom.typeB;
    }
}

DispersionCorrection::TopologyParams::TopologyParams(const gmx_mtop_t&         mtop,
                                                     const t_inputrec&         inputrec,
                                                     bool                      useBuckingham,
                                                     int                       numAtomTypes,
                                                     gmx::ArrayRef<const real> nonbondedForceParameters)
{
    const int      ntp   = numAtomTypes;
    const gmx_bool bBHAM = useBuckingham;

    gmx::ArrayRef<const real> nbfp = nonbondedForceParameters;
    std::vector<real>         nbfp_comb;
    /* For LJ-PME, we want to correct for the difference between the
     * actual C6 values and the C6 values used by the LJ-PME based on
     * combination rules. */
    if (EVDW_PME(inputrec.vdwtype))
    {
        nbfp_comb = mk_nbfp_combination_rule(
                mtop.ffparams,
                (inputrec.ljpme_combination_rule == eljpmeLB) ? eCOMB_ARITHMETIC : eCOMB_GEOMETRIC);
        for (int tpi = 0; tpi < ntp; ++tpi)
        {
            for (int tpj = 0; tpj < ntp; ++tpj)
            {
                C6(nbfp_comb, ntp, tpi, tpj) = C6(nbfp, ntp, tpi, tpj) - C6(nbfp_comb, ntp, tpi, tpj);
                C12(nbfp_comb, ntp, tpi, tpj) = C12(nbfp, ntp, tpi, tpj);
            }
        }
        nbfp = nbfp_comb;
    }

    for (int q = 0; q < (inputrec.efep == efepNO ? 1 : 2); q++)
    {
        double  csix    = 0;
        double  ctwelve = 0;
        int64_t npair   = 0;
        int64_t nexcl   = 0;
        if (!EI_TPI(inputrec.eI))
        {
            numAtomsForDensity_ = mtop.natoms;
            numCorrections_     = 0.5 * mtop.natoms;

            /* Count the types so we avoid natoms^2 operations */
            std::vector<int> typecount(ntp);
            gmx_mtop_count_atomtypes(&mtop, q, typecount.data());

            for (int tpi = 0; tpi < ntp; tpi++)
            {
                for (int tpj = tpi; tpj < ntp; tpj++)
                {
                    const int64_t iCount = typecount[tpi];
                    const int64_t jCount = typecount[tpj];
                    int64_t       npair_ij;
                    if (tpi != tpj)
                    {
                        npair_ij = iCount * jCount;
                    }
                    else
                    {
                        npair_ij = iCount * (iCount - 1) / 2;
                    }
                    if (bBHAM)
                    {
                        /* nbfp now includes the 6.0 derivative prefactor */
                        csix += npair_ij * BHAMC(nbfp, ntp, tpi, tpj) / 6.0;
                    }
                    else
                    {
                        /* nbfp now includes the 6.0/12.0 derivative prefactors */
                        csix += npair_ij * C6(nbfp, ntp, tpi, tpj) / 6.0;
                        ctwelve += npair_ij * C12(nbfp, ntp, tpi, tpj) / 12.0;
                    }
                    npair += npair_ij;
                }
            }
            /* Subtract the excluded pairs.
             * The main reason for substracting exclusions is that in some cases
             * some combinations might never occur and the parameters could have
             * any value. These unused values should not influence the dispersion
             * correction.
             */
            for (const gmx_molblock_t& molb : mtop.molblock)
            {
                const int      nmol  = molb.nmol;
                const t_atoms* atoms = &mtop.moltype[molb.type].atoms;
                const auto&    excl  = mtop.moltype[molb.type].excls;
                for (int i = 0; (i < atoms->nr); i++)
                {
                    const int tpi = atomtypeAOrB(atoms->atom[i], q);
                    for (const int k : excl[i])
                    {
                        if (k > i)
                        {
                            const int tpj = atomtypeAOrB(atoms->atom[k], q);
                            if (bBHAM)
                            {
                                /* nbfp now includes the 6.0 derivative prefactor */
                                csix -= nmol * BHAMC(nbfp, ntp, tpi, tpj) / 6.0;
                            }
                            else
                            {
                                /* nbfp now includes the 6.0/12.0 derivative prefactors */
                                csix -= nmol * C6(nbfp, ntp, tpi, tpj) / 6.0;
                                ctwelve -= nmol * C12(nbfp, ntp, tpi, tpj) / 12.0;
                            }
                            nexcl += molb.nmol;
                        }
                    }
                }
            }
        }
        else
        {
            const t_atoms& atoms_tpi = mtop.moltype[mtop.molblock.back().type].atoms;

            /* Only correct for the interaction of the test particle
             * with the rest of the system.
             */
            numAtomsForDensity_ = mtop.natoms - atoms_tpi.nr;
            numCorrections_     = atoms_tpi.nr;

            npair = 0;
            for (size_t mb = 0; mb < mtop.molblock.size(); mb++)
            {
                const gmx_molblock_t& molb  = mtop.molblock[mb];
                const t_atoms&        atoms = mtop.moltype[molb.type].atoms;
                for (int j = 0; j < atoms.nr; j++)
                {
                    int nmolc = molb.nmol;
                    /* Remove the interaction of the test charge group
                     * with itself.
                     */
                    if (mb == mtop.molblock.size() - 1)
                    {
                        nmolc--;

                        if (mb == 0 && molb.nmol == 1)
                        {
                            gmx_fatal(FARGS,
                                      "Old format tpr with TPI, please generate a new tpr file");
                        }
                    }
                    const int tpj = atomtypeAOrB(atoms.atom[j], q);
                    for (int i = 0; i < atoms_tpi.nr; i++)
                    {
                        const int tpi = atomtypeAOrB(atoms_tpi.atom[i], q);
                        if (bBHAM)
                        {
                            /* nbfp now includes the 6.0 derivative prefactor */
                            csix += nmolc * BHAMC(nbfp, ntp, tpi, tpj) / 6.0;
                        }
                        else
                        {
                            /* nbfp now includes the 6.0/12.0 derivative prefactors */
                            csix += nmolc * C6(nbfp, ntp, tpi, tpj) / 6.0;
                            ctwelve += nmolc * C12(nbfp, ntp, tpi, tpj) / 12.0;
                        }
                        npair += nmolc;
                    }
                }
            }
        }
        if (npair - nexcl <= 0)
        {
            csix    = 0;
            ctwelve = 0;
        }
        else
        {
            csix /= npair - nexcl;
            ctwelve /= npair - nexcl;
        }
        if (debug)
        {
            fprintf(debug, "Counted %" PRId64 " exclusions\n", nexcl);
            fprintf(debug, "Average C6 parameter is: %10g\n", csix);
            fprintf(debug, "Average C12 parameter is: %10g\n", ctwelve);
        }
        avcsix_[q]    = csix;
        avctwelve_[q] = ctwelve;
    }
}

static void integrate_table(const real vdwtab[],
                            const real scale,
                            const int  offstart,
                            const int  rstart,
                            const int  rend,
                            double*    enerout,
                            double*    virout)
{
    const double invscale  = 1.0 / scale;
    const double invscale2 = invscale * invscale;
    const double invscale3 = invscale * invscale2;

    /* Following summation derived from cubic spline definition,
     * Numerical Recipies in C, second edition, p. 113-116.  Exact for
     * the cubic spline.  We first calculate the negative of the
     * energy from rvdw to rvdw_switch, assuming that g(r)=1, and then
     * add the more standard, abrupt cutoff correction to that result,
     * yielding the long-range correction for a switched function.  We
     * perform both the pressure and energy loops at the same time for
     * simplicity, as the computational cost is low. */

    /* Since the dispersion table has been scaled down a factor
     * 6.0 and the repulsion a factor 12.0 to compensate for the
     * c6/c12 parameters inside nbfp[] being scaled up (to save
     * flops in kernels), we need to correct for this.
     */
    const double tabfactor = (offstart == 0 ? 6.0 : 12.0);

    double enersum = 0.0;
    double virsum  = 0.0;
    for (int ri = rstart; ri < rend; ++ri)
    {
        const double r  = ri * invscale;
        const double ea = invscale3;
        const double eb = 2.0 * invscale2 * r;
        const double ec = invscale * r * r;

        const double pa = invscale3;
        const double pb = 3.0 * invscale2 * r;
        const double pc = 3.0 * invscale * r * r;
        const double pd = r * r * r;

        /* this "8" is from the packing in the vdwtab array - perhaps
           should be defined? */

        const int    offset = 8 * ri + offstart;
        const double y0     = vdwtab[offset];
        const double f      = vdwtab[offset + 1];
        const double g      = vdwtab[offset + 2];
        const double h      = vdwtab[offset + 3];

        enersum += y0 * (ea / 3 + eb / 2 + ec) + f * (ea / 4 + eb / 3 + ec / 2)
                   + g * (ea / 5 + eb / 4 + ec / 3) + h * (ea / 6 + eb / 5 + ec / 4);
        virsum += f * (pa / 4 + pb / 3 + pc / 2 + pd) + 2 * g * (pa / 5 + pb / 4 + pc / 3 + pd / 2)
                  + 3 * h * (pa / 6 + pb / 5 + pc / 4 + pd / 3);
    }
    *enerout = 4.0 * M_PI * enersum * tabfactor;
    *virout  = 4.0 * M_PI * virsum * tabfactor;
}

/* Struct for storing and passing energy or virial corrections */
struct InteractionCorrection
{
    real dispersion = 0;
    real repulsion  = 0;
};

/* Adds the energy and virial corrections beyond the cut-off */
static void addCorrectionBeyondCutoff(InteractionCorrection* energy,
                                      InteractionCorrection* virial,
                                      const double           cutoffDistance)
{
    const double rc3 = cutoffDistance * cutoffDistance * cutoffDistance;
    const double rc9 = rc3 * rc3 * rc3;

    energy->dispersion += -4.0 * M_PI / (3.0 * rc3);
    energy->repulsion += 4.0 * M_PI / (9.0 * rc9);
    virial->dispersion += 8.0 * M_PI / rc3;
    virial->repulsion += -16.0 * M_PI / (3.0 * rc9);
}

void DispersionCorrection::setInteractionParameters(InteractionParams*         iParams,
                                                    const interaction_const_t& ic,
                                                    const char*                tableFileName)
{
    /* We only need to set the tables at first call, i.e. tableFileName!=nullptr
     * or when we changed the cut-off with LJ-PME tuning.
     */
    if (tableFileName || EVDW_PME(ic.vdwtype))
    {
        iParams->dispersionCorrectionTable_ =
                makeDispersionCorrectionTable(nullptr, &ic, ic.rvdw, tableFileName);
    }

    InteractionCorrection energy;
    InteractionCorrection virial;

    if ((ic.vdw_modifier == eintmodPOTSHIFT) || (ic.vdw_modifier == eintmodPOTSWITCH)
        || (ic.vdw_modifier == eintmodFORCESWITCH) || (ic.vdwtype == evdwSHIFT)
        || (ic.vdwtype == evdwSWITCH))
    {
        if (((ic.vdw_modifier == eintmodPOTSWITCH) || (ic.vdw_modifier == eintmodFORCESWITCH)
             || (ic.vdwtype == evdwSWITCH))
            && ic.rvdw_switch == 0)
        {
            gmx_fatal(FARGS,
                      "With dispersion correction rvdw-switch can not be zero "
                      "for vdw-type = %s",
                      evdw_names[ic.vdwtype]);
        }

        GMX_ASSERT(iParams->dispersionCorrectionTable_, "We need an initialized table");

        /* TODO This code depends on the logic in tables.c that
           constructs the table layout, which should be made
           explicit in future cleanup. */
        GMX_ASSERT(iParams->dispersionCorrectionTable_->interaction == GMX_TABLE_INTERACTION_VDWREP_VDWDISP,
                   "Dispersion-correction code needs a table with both repulsion and dispersion "
                   "terms");
        const real  scale  = iParams->dispersionCorrectionTable_->scale;
        const real* vdwtab = iParams->dispersionCorrectionTable_->data;

        /* Round the cut-offs to exact table values for precision */
        int ri0 = static_cast<int>(std::floor(ic.rvdw_switch * scale));
        int ri1 = static_cast<int>(std::ceil(ic.rvdw * scale));

        /* The code below has some support for handling force-switching, i.e.
         * when the force (instead of potential) is switched over a limited
         * region. This leads to a constant shift in the potential inside the
         * switching region, which we can handle by adding a constant energy
         * term in the force-switch case just like when we do potential-shift.
         *
         * For now this is not enabled, but to keep the functionality in the
         * code we check separately for switch and shift. When we do force-switch
         * the shifting point is rvdw_switch, while it is the cutoff when we
         * have a classical potential-shift.
         *
         * For a pure potential-shift the potential has a constant shift
         * all the way out to the cutoff, and that is it. For other forms
         * we need to calculate the constant shift up to the point where we
         * start modifying the potential.
         */
        ri0 = (ic.vdw_modifier == eintmodPOTSHIFT) ? ri1 : ri0;

        const double r0  = ri0 / scale;
        const double rc3 = r0 * r0 * r0;
        const double rc9 = rc3 * rc3 * rc3;

        if ((ic.vdw_modifier == eintmodFORCESWITCH) || (ic.vdwtype == evdwSHIFT))
        {
            /* Determine the constant energy shift below rvdw_switch.
             * Table has a scale factor since we have scaled it down to compensate
             * for scaling-up c6/c12 with the derivative factors to save flops in analytical kernels.
             */
            iParams->enershiftsix_ = static_cast<real>(-1.0 / (rc3 * rc3)) - 6.0 * vdwtab[8 * ri0];
            iParams->enershifttwelve_ = static_cast<real>(1.0 / (rc9 * rc3)) - 12.0 * vdwtab[8 * ri0 + 4];
        }
        else if (ic.vdw_modifier == eintmodPOTSHIFT)
        {
            iParams->enershiftsix_    = static_cast<real>(-1.0 / (rc3 * rc3));
            iParams->enershifttwelve_ = static_cast<real>(1.0 / (rc9 * rc3));
        }

        /* Add the constant part from 0 to rvdw_switch.
         * This integration from 0 to rvdw_switch overcounts the number
         * of interactions by 1, as it also counts the self interaction.
         * We will correct for this later.
         */
        energy.dispersion += 4.0 * M_PI * iParams->enershiftsix_ * rc3 / 3.0;
        energy.repulsion += 4.0 * M_PI * iParams->enershifttwelve_ * rc3 / 3.0;

        /* Calculate the contribution in the range [r0,r1] where we
         * modify the potential. For a pure potential-shift modifier we will
         * have ri0==ri1, and there will not be any contribution here.
         */
        double enersum = 0;
        double virsum  = 0;
        integrate_table(vdwtab, scale, 0, ri0, ri1, &enersum, &virsum);
        energy.dispersion -= enersum;
        virial.dispersion -= virsum;
        integrate_table(vdwtab, scale, 4, ri0, ri1, &enersum, &virsum);
        energy.repulsion -= enersum;
        virial.repulsion -= virsum;


        /* Alright: Above we compensated by REMOVING the parts outside r0
         * corresponding to the ideal VdW 1/r6 and /r12 potentials.
         *
         * Regardless of whether r0 is the point where we start switching,
         * or the cutoff where we calculated the constant shift, we include
         * all the parts we are missing out to infinity from r0 by
         * calculating the analytical dispersion correction.
         */
        addCorrectionBeyondCutoff(&energy, &virial, r0);
    }
    else if (ic.vdwtype == evdwCUT || EVDW_PME(ic.vdwtype) || ic.vdwtype == evdwUSER)
    {
        /* Note that with LJ-PME, the dispersion correction is multiplied
         * by the difference between the actual C6 and the value of C6
         * that would produce the combination rule.
         * This means the normal energy and virial difference formulas
         * can be used here.
         */

        const double rc3 = ic.rvdw * ic.rvdw * ic.rvdw;
        const double rc9 = rc3 * rc3 * rc3;
        if (ic.vdw_modifier == eintmodPOTSHIFT)
        {
            /* Contribution within the cut-off */
            energy.dispersion += -4.0 * M_PI / (3.0 * rc3);
            energy.repulsion += 4.0 * M_PI / (3.0 * rc9);
        }
        /* Contribution beyond the cut-off */
        addCorrectionBeyondCutoff(&energy, &virial, ic.rvdw);
    }
    else
    {
        gmx_fatal(FARGS, "Dispersion correction is not implemented for vdw-type = %s",
                  evdw_names[ic.vdwtype]);
    }

    iParams->enerdiffsix_    = energy.dispersion;
    iParams->enerdifftwelve_ = energy.repulsion;
    /* The 0.5 is due to the Gromacs definition of the virial */
    iParams->virdiffsix_    = 0.5 * virial.dispersion;
    iParams->virdifftwelve_ = 0.5 * virial.repulsion;
}

DispersionCorrection::DispersionCorrection(const gmx_mtop_t&          mtop,
                                           const t_inputrec&          inputrec,
                                           bool                       useBuckingham,
                                           int                        numAtomTypes,
                                           gmx::ArrayRef<const real>  nonbondedForceParameters,
                                           const interaction_const_t& ic,
                                           const char*                tableFileName) :
    eDispCorr_(inputrec.eDispCorr),
    vdwType_(inputrec.vdwtype),
    eFep_(inputrec.efep),
    topParams_(mtop, inputrec, useBuckingham, numAtomTypes, nonbondedForceParameters)
{
    if (eDispCorr_ != edispcNO)
    {
        GMX_RELEASE_ASSERT(tableFileName, "Need a table file name");

        setInteractionParameters(&iParams_, ic, tableFileName);
    }
}

bool DispersionCorrection::correctFullInteraction() const
{
    return (eDispCorr_ == edispcAllEner || eDispCorr_ == edispcAllEnerPres);
}

void DispersionCorrection::print(const gmx::MDLogger& mdlog) const
{
    if (topParams_.avcsix_[0] == 0 && topParams_.avctwelve_[0] == 0)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText("WARNING: There are no atom pairs for dispersion correction");
    }
    else if (vdwType_ == evdwUSER)
    {
        GMX_LOG(mdlog.warning)
                .asParagraph()
                .appendText("WARNING: using dispersion correction with user tables\n");
    }

    std::string text = gmx::formatString("Long Range LJ corr.: <C6> %10.4e", topParams_.avcsix_[0]);
    if (correctFullInteraction())
    {
        text += gmx::formatString(" <C12> %10.4e", topParams_.avctwelve_[0]);
    }
    GMX_LOG(mdlog.info).appendText(text);
}

void DispersionCorrection::setParameters(const interaction_const_t& ic)
{
    if (eDispCorr_ != edispcNO)
    {
        setInteractionParameters(&iParams_, ic, nullptr);
    }
}

DispersionCorrection::Correction DispersionCorrection::calculate(const matrix box, const real lambda) const

{
    Correction corr;

    if (eDispCorr_ == edispcNO)
    {
        return corr;
    }

    const bool bCorrAll  = correctFullInteraction();
    const bool bCorrPres = (eDispCorr_ == edispcEnerPres || eDispCorr_ == edispcAllEnerPres);

    const real invvol  = 1 / det(box);
    const real density = topParams_.numAtomsForDensity_ * invvol;
    const real numCorr = topParams_.numCorrections_;

    real avcsix;
    real avctwelve;
    if (eFep_ == efepNO)
    {
        avcsix    = topParams_.avcsix_[0];
        avctwelve = topParams_.avctwelve_[0];
    }
    else
    {
        avcsix    = (1 - lambda) * topParams_.avcsix_[0] + lambda * topParams_.avcsix_[1];
        avctwelve = (1 - lambda) * topParams_.avctwelve_[0] + lambda * topParams_.avctwelve_[1];
    }

    const real enerdiff = numCorr * (density * iParams_.enerdiffsix_ - iParams_.enershiftsix_);
    corr.energy += avcsix * enerdiff;
    real dvdlambda = 0;
    if (eFep_ != efepNO)
    {
        dvdlambda += (topParams_.avcsix_[1] - topParams_.avcsix_[0]) * enerdiff;
    }
    if (bCorrAll)
    {
        const real enerdiff = numCorr * (density * iParams_.enerdifftwelve_ - iParams_.enershifttwelve_);
        corr.energy += avctwelve * enerdiff;
        if (eFep_ != efepNO)
        {
            dvdlambda += (topParams_.avctwelve_[1] - topParams_.avctwelve_[0]) * enerdiff;
        }
    }

    if (bCorrPres)
    {
        corr.virial = numCorr * density * avcsix * iParams_.virdiffsix_ / 3.0;
        if (eDispCorr_ == edispcAllEnerPres)
        {
            corr.virial += numCorr * density * avctwelve * iParams_.virdifftwelve_ / 3.0;
        }
        /* The factor 2 is because of the Gromacs virial definition */
        corr.pressure = -2.0 * invvol * corr.virial * PRESFAC;
    }

    if (eFep_ != efepNO)
    {
        corr.dvdl += dvdlambda;
    }

    return corr;
}
