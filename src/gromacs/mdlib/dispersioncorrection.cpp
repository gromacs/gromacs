/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2019, by the GROMACS development team, led by
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
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/mdtypes/nblist.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/arrayref.h"
#include "gromacs/utility/fatalerror.h"

/* Returns a matrix, as flat list, of combination rule combined LJ parameters */
static std::vector<real> mk_nbfp_combination_rule(const gmx_ffparams_t &ffparams,
                                                  const int             comb_rule)
{
    const int         atnr = ffparams.atnr;

    std::vector<real> nbfp(atnr*atnr*2);

    for (int i = 0; i < atnr; ++i)
    {
        for (int j = 0; j < atnr; ++j)
        {
            const real c6i  = ffparams.iparams[i*(atnr + 1)].lj.c6;
            const real c12i = ffparams.iparams[i*(atnr + 1)].lj.c12;
            const real c6j  = ffparams.iparams[j*(atnr + 1)].lj.c6;
            const real c12j = ffparams.iparams[j*(atnr + 1)].lj.c12;
            real       c6   = std::sqrt(c6i  * c6j);
            real       c12  = std::sqrt(c12i * c12j);
            if (comb_rule == eCOMB_ARITHMETIC
                && !gmx_numzero(c6) && !gmx_numzero(c12))
            {
                const real sigmai = gmx::sixthroot(c12i / c6i);
                const real sigmaj = gmx::sixthroot(c12j / c6j);
                const real epsi   = c6i * c6i / c12i;
                const real epsj   = c6j * c6j / c12j;
                const real sigma  = 0.5*(sigmai + sigmaj);
                const real eps    = std::sqrt(epsi * epsj);
                c6                = eps * gmx::power6(sigma);
                c12               = eps * gmx::power12(sigma);
            }
            C6(nbfp, atnr, i, j)  = c6*6.0;
            C12(nbfp, atnr, i, j) = c12*12.0;
        }
    }

    return nbfp;
}

/* Returns the A-topology atom type when aOrB=0, the B-topology atom type when aOrB=1 */
static int atomtypeAOrB(const t_atom &atom,
                        int           aOrB)
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

void set_avcsixtwelve(FILE             *fplog,
                      t_forcerec       *fr,
                      const gmx_mtop_t *mtop)
{
    const int                 ntp   = fr->ntype;
    const gmx_bool            bBHAM = fr->bBHAM;

    gmx::ArrayRef<const real> nbfp  = gmx::arrayRefFromArray(fr->nbfp, ntp*ntp*2);
    std::vector<real>         nbfp_comb;
    /* For LJ-PME, we want to correct for the difference between the
     * actual C6 values and the C6 values used by the LJ-PME based on
     * combination rules. */
    if (EVDW_PME(fr->ic->vdwtype))
    {
        nbfp_comb = mk_nbfp_combination_rule(mtop->ffparams,
                                             (fr->ljpme_combination_rule == eljpmeLB) ? eCOMB_ARITHMETIC : eCOMB_GEOMETRIC);
        for (int tpi = 0; tpi < ntp; ++tpi)
        {
            for (int tpj = 0; tpj < ntp; ++tpj)
            {
                C6(nbfp_comb, ntp, tpi, tpj) =
                    C6(nbfp, ntp, tpi, tpj) - C6(nbfp_comb, ntp, tpi, tpj);
                C12(nbfp_comb, ntp, tpi, tpj) = C12(nbfp, ntp, tpi, tpj);
            }
        }
        nbfp = nbfp_comb;
    }

    for (int q = 0; q < (fr->efep == efepNO ? 1 : 2); q++)
    {
        double csix    = 0;
        double ctwelve = 0;
        int    npair   = 0;
        int    nexcl   = 0;
        if (!fr->n_tpi)
        {
            /* Count the types so we avoid natoms^2 operations */
            std::vector<int> typecount(ntp);
            gmx_mtop_count_atomtypes(mtop, q, typecount.data());

            for (int tpi = 0; tpi < ntp; tpi++)
            {
                for (int tpj = tpi; tpj < ntp; tpj++)
                {
                    const int iCount = typecount[tpi];
                    const int jCount = typecount[tpj];
                    int       npair_ij;
                    if (tpi != tpj)
                    {
                        npair_ij = iCount*jCount;
                    }
                    else
                    {
                        npair_ij = iCount*(iCount - 1)/2;
                    }
                    if (bBHAM)
                    {
                        /* nbfp now includes the 6.0 derivative prefactor */
                        csix    += npair_ij*BHAMC(nbfp, ntp, tpi, tpj)/6.0;
                    }
                    else
                    {
                        /* nbfp now includes the 6.0/12.0 derivative prefactors */
                        csix    += npair_ij*   C6(nbfp, ntp, tpi, tpj)/6.0;
                        ctwelve += npair_ij*  C12(nbfp, ntp, tpi, tpj)/12.0;
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
            for (const gmx_molblock_t &molb : mtop->molblock)
            {
                const int       nmol  = molb.nmol;
                const t_atoms  *atoms = &mtop->moltype[molb.type].atoms;
                const t_blocka *excl  = &mtop->moltype[molb.type].excls;
                for (int i = 0; (i < atoms->nr); i++)
                {
                    const int tpi = atomtypeAOrB( atoms->atom[i], q);
                    const int j1  = excl->index[i];
                    const int j2  = excl->index[i+1];
                    for (int j = j1; j < j2; j++)
                    {
                        const int k = excl->a[j];
                        if (k > i)
                        {
                            const int tpj = atomtypeAOrB(atoms->atom[k], q);
                            if (bBHAM)
                            {
                                /* nbfp now includes the 6.0 derivative prefactor */
                                csix -= nmol*BHAMC(nbfp, ntp, tpi, tpj)/6.0;
                            }
                            else
                            {
                                /* nbfp now includes the 6.0/12.0 derivative prefactors */
                                csix    -= nmol*C6 (nbfp, ntp, tpi, tpj)/6.0;
                                ctwelve -= nmol*C12(nbfp, ntp, tpi, tpj)/12.0;
                            }
                            nexcl += molb.nmol;
                        }
                    }
                }
            }
        }
        else
        {
            /* Only correct for the interaction of the test particle
             * with the rest of the system.
             */
            const t_atoms *atoms_tpi =
                &mtop->moltype[mtop->molblock.back().type].atoms;

            npair = 0;
            for (size_t mb = 0; mb < mtop->molblock.size(); mb++)
            {
                const gmx_molblock_t &molb  = mtop->molblock[mb];
                const t_atoms        *atoms = &mtop->moltype[molb.type].atoms;
                for (int j = 0; j < atoms->nr; j++)
                {
                    int nmolc = molb.nmol;
                    /* Remove the interaction of the test charge group
                     * with itself.
                     */
                    if (mb == mtop->molblock.size() - 1)
                    {
                        nmolc--;

                        if (mb == 0 && molb.nmol == 1)
                        {
                            gmx_fatal(FARGS, "Old format tpr with TPI, please generate a new tpr file");
                        }
                    }
                    const int tpj = atomtypeAOrB(atoms->atom[j], q);
                    for (int i = 0; i < fr->n_tpi; i++)
                    {
                        const int tpi = atomtypeAOrB(atoms_tpi->atom[i], q);
                        if (bBHAM)
                        {
                            /* nbfp now includes the 6.0 derivative prefactor */
                            csix    += nmolc*BHAMC(nbfp, ntp, tpi, tpj)/6.0;
                        }
                        else
                        {
                            /* nbfp now includes the 6.0/12.0 derivative prefactors */
                            csix    += nmolc*C6 (nbfp, ntp, tpi, tpj)/6.0;
                            ctwelve += nmolc*C12(nbfp, ntp, tpi, tpj)/12.0;
                        }
                        npair += nmolc;
                    }
                }
            }
        }
        if (npair - nexcl <= 0 && fplog)
        {
            fprintf(fplog, "\nWARNING: There are no atom pairs for dispersion correction\n\n");
            csix     = 0;
            ctwelve  = 0;
        }
        else
        {
            csix    /= npair - nexcl;
            ctwelve /= npair - nexcl;
        }
        if (debug)
        {
            fprintf(debug, "Counted %d exclusions\n", nexcl);
            fprintf(debug, "Average C6 parameter is: %10g\n", csix);
            fprintf(debug, "Average C12 parameter is: %10g\n", ctwelve);
        }
        fr->avcsix[q]    = csix;
        fr->avctwelve[q] = ctwelve;
    }

    if (fplog != nullptr)
    {
        if (fr->eDispCorr == edispcAllEner ||
            fr->eDispCorr == edispcAllEnerPres)
        {
            fprintf(fplog, "Long Range LJ corr.: <C6> %10.4e, <C12> %10.4e\n",
                    fr->avcsix[0], fr->avctwelve[0]);
        }
        else
        {
            fprintf(fplog, "Long Range LJ corr.: <C6> %10.4e\n", fr->avcsix[0]);
        }
    }
}

static void
integrate_table(const real  vdwtab[],
                const real  scale,
                const int   offstart,
                const int   rstart,
                const int   rend,
                double     *enerout,
                double     *virout)
{
    const double invscale  = 1.0/scale;
    const double invscale2 = invscale*invscale;
    const double invscale3 = invscale*invscale2;

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

    double       enersum = 0.0;
    double       virsum  = 0.0;
    for (int ri = rstart; ri < rend; ++ri)
    {
        const double r  = ri*invscale;
        const double ea = invscale3;
        const double eb = 2.0*invscale2*r;
        const double ec = invscale*r*r;

        const double pa = invscale3;
        const double pb = 3.0*invscale2*r;
        const double pc = 3.0*invscale*r*r;
        const double pd = r*r*r;

        /* this "8" is from the packing in the vdwtab array - perhaps
           should be defined? */

        const int    offset = 8*ri + offstart;
        const double y0     = vdwtab[offset];
        const double f      = vdwtab[offset+1];
        const double g      = vdwtab[offset+2];
        const double h      = vdwtab[offset+3];

        enersum += y0*(ea/3 + eb/2 + ec) + f*(ea/4 + eb/3 + ec/2) + g*(ea/5 + eb/4 + ec/3) + h*(ea/6 + eb/5 + ec/4);
        virsum  +=  f*(pa/4 + pb/3 + pc/2 + pd) + 2*g*(pa/5 + pb/4 + pc/3 + pd/2) + 3*h*(pa/6 + pb/5 + pc/4 + pd/3);
    }
    *enerout = 4.0*M_PI*enersum*tabfactor;
    *virout  = 4.0*M_PI*virsum*tabfactor;
}

void calc_enervirdiff(FILE       *fplog,
                      const int   eDispCorr,
                      t_forcerec *fr)
{
    fr->enershiftsix    = 0;
    fr->enershifttwelve = 0;
    fr->enerdiffsix     = 0;
    fr->enerdifftwelve  = 0;
    fr->virdiffsix      = 0;
    fr->virdifftwelve   = 0;

    const interaction_const_t *ic = fr->ic;

    if (eDispCorr != edispcNO)
    {
        double eners[2] = { 0 };
        double virs[2]  = { 0 };

        if ((ic->vdw_modifier == eintmodPOTSHIFT) ||
            (ic->vdw_modifier == eintmodPOTSWITCH) ||
            (ic->vdw_modifier == eintmodFORCESWITCH) ||
            (ic->vdwtype == evdwSHIFT) ||
            (ic->vdwtype == evdwSWITCH))
        {
            if (((ic->vdw_modifier == eintmodPOTSWITCH) ||
                 (ic->vdw_modifier == eintmodFORCESWITCH) ||
                 (ic->vdwtype == evdwSWITCH)) && ic->rvdw_switch == 0)
            {
                gmx_fatal(FARGS,
                          "With dispersion correction rvdw-switch can not be zero "
                          "for vdw-type = %s", evdw_names[ic->vdwtype]);
            }

            /* TODO This code depends on the logic in tables.c that
               constructs the table layout, which should be made
               explicit in future cleanup. */
            GMX_ASSERT(fr->dispersionCorrectionTable->interaction == GMX_TABLE_INTERACTION_VDWREP_VDWDISP,
                       "Dispersion-correction code needs a table with both repulsion and dispersion terms");
            const real  scale  = fr->dispersionCorrectionTable->scale;
            const real *vdwtab = fr->dispersionCorrectionTable->data;

            /* Round the cut-offs to exact table values for precision */
            int ri0  = static_cast<int>(std::floor(ic->rvdw_switch*scale));
            int ri1  = static_cast<int>(std::ceil(ic->rvdw*scale));

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
            ri0               = (ic->vdw_modifier == eintmodPOTSHIFT) ? ri1 : ri0;

            const double r0   = ri0/scale;
            const double rc3  = r0*r0*r0;
            const double rc9  = rc3*rc3*rc3;

            if ((ic->vdw_modifier == eintmodFORCESWITCH) ||
                (ic->vdwtype == evdwSHIFT))
            {
                /* Determine the constant energy shift below rvdw_switch.
                 * Table has a scale factor since we have scaled it down to compensate
                 * for scaling-up c6/c12 with the derivative factors to save flops in analytical kernels.
                 */
                fr->enershiftsix    = static_cast<real>(-1.0/(rc3*rc3)) - 6.0*vdwtab[8*ri0];
                fr->enershifttwelve = static_cast<real>( 1.0/(rc9*rc3)) - 12.0*vdwtab[8*ri0 + 4];
            }
            else if (ic->vdw_modifier == eintmodPOTSHIFT)
            {
                fr->enershiftsix    = static_cast<real>(-1.0/(rc3*rc3));
                fr->enershifttwelve = static_cast<real>( 1.0/(rc9*rc3));
            }

            /* Add the constant part from 0 to rvdw_switch.
             * This integration from 0 to rvdw_switch overcounts the number
             * of interactions by 1, as it also counts the self interaction.
             * We will correct for this later.
             */
            eners[0] += 4.0*M_PI*fr->enershiftsix*rc3/3.0;
            eners[1] += 4.0*M_PI*fr->enershifttwelve*rc3/3.0;

            /* Calculate the contribution in the range [r0,r1] where we
             * modify the potential. For a pure potential-shift modifier we will
             * have ri0==ri1, and there will not be any contribution here.
             */
            for (int i = 0; i < 2; i++)
            {
                double enersum = 0;
                double virsum  = 0;
                integrate_table(vdwtab, scale, (i == 0 ? 0 : 4), ri0, ri1, &enersum, &virsum);
                eners[i] -= enersum;
                virs[i]  -= virsum;
            }

            /* Alright: Above we compensated by REMOVING the parts outside r0
             * corresponding to the ideal VdW 1/r6 and /r12 potentials.
             *
             * Regardless of whether r0 is the point where we start switching,
             * or the cutoff where we calculated the constant shift, we include
             * all the parts we are missing out to infinity from r0 by
             * calculating the analytical dispersion correction.
             */
            eners[0] += -4.0*M_PI/(3.0*rc3);
            eners[1] +=  4.0*M_PI/(9.0*rc9);
            virs[0]  +=  8.0*M_PI/rc3;
            virs[1]  += -16.0*M_PI/(3.0*rc9);
        }
        else if (ic->vdwtype == evdwCUT ||
                 EVDW_PME(ic->vdwtype) ||
                 ic->vdwtype == evdwUSER)
        {
            if (ic->vdwtype == evdwUSER && fplog)
            {
                fprintf(fplog,
                        "WARNING: using dispersion correction with user tables\n");
            }

            /* Note that with LJ-PME, the dispersion correction is multiplied
             * by the difference between the actual C6 and the value of C6
             * that would produce the combination rule.
             * This means the normal energy and virial difference formulas
             * can be used here.
             */

            const double rc3 = ic->rvdw*ic->rvdw*ic->rvdw;
            const double rc9 = rc3*rc3*rc3;
            /* Contribution beyond the cut-off */
            eners[0] += -4.0*M_PI/(3.0*rc3);
            eners[1] +=  4.0*M_PI/(9.0*rc9);
            if (ic->vdw_modifier == eintmodPOTSHIFT)
            {
                /* Contribution within the cut-off */
                eners[0] += -4.0*M_PI/(3.0*rc3);
                eners[1] +=  4.0*M_PI/(3.0*rc9);
            }
            /* Contribution beyond the cut-off */
            virs[0]  +=  8.0*M_PI/rc3;
            virs[1]  += -16.0*M_PI/(3.0*rc9);
        }
        else
        {
            gmx_fatal(FARGS,
                      "Dispersion correction is not implemented for vdw-type = %s",
                      evdw_names[ic->vdwtype]);
        }

        /* When we deprecate the group kernels the code below can go too */
        if (ic->vdwtype == evdwPME && fr->cutoff_scheme == ecutsGROUP)
        {
            /* Calculate self-interaction coefficient (assuming that
             * the reciprocal-space contribution is constant in the
             * region that contributes to the self-interaction).
             */
            fr->enershiftsix = gmx::power6(ic->ewaldcoeff_lj) / 6.0;

            eners[0] += -gmx::power3(std::sqrt(M_PI)*ic->ewaldcoeff_lj)/3.0;
            virs[0]  +=  gmx::power3(std::sqrt(M_PI)*ic->ewaldcoeff_lj);
        }

        fr->enerdiffsix    = eners[0];
        fr->enerdifftwelve = eners[1];
        /* The 0.5 is due to the Gromacs definition of the virial */
        fr->virdiffsix     = 0.5*virs[0];
        fr->virdifftwelve  = 0.5*virs[1];
    }
}

void calc_dispcorr(const t_inputrec *ir, const t_forcerec *fr,
                   const matrix box, real lambda, tensor pres, tensor virial,
                   real *prescorr, real *enercorr, real *dvdlcorr)
{
    *prescorr = 0;
    *enercorr = 0;
    *dvdlcorr = 0;

    clear_mat(virial);
    clear_mat(pres);

    if (ir->eDispCorr != edispcNO)
    {
        const bool bCorrAll  = (ir->eDispCorr == edispcAllEner ||
                                ir->eDispCorr == edispcAllEnerPres);
        const bool bCorrPres = (ir->eDispCorr == edispcEnerPres ||
                                ir->eDispCorr == edispcAllEnerPres);

        const real invvol = 1/det(box);
        real       dens;
        real       ninter;
        if (fr->n_tpi)
        {
            /* Only correct for the interactions with the inserted molecule */
            dens   = (fr->numAtomsForDispersionCorrection - fr->n_tpi)*invvol;
            ninter = fr->n_tpi;
        }
        else
        {
            dens   = fr->numAtomsForDispersionCorrection*invvol;
            ninter = 0.5*fr->numAtomsForDispersionCorrection;
        }

        real avcsix;
        real avctwelve;
        if (ir->efep == efepNO)
        {
            avcsix    = fr->avcsix[0];
            avctwelve = fr->avctwelve[0];
        }
        else
        {
            avcsix    = (1 - lambda)*fr->avcsix[0]    + lambda*fr->avcsix[1];
            avctwelve = (1 - lambda)*fr->avctwelve[0] + lambda*fr->avctwelve[1];
        }

        const real enerdiff   = ninter*(dens*fr->enerdiffsix - fr->enershiftsix);
        *enercorr            += avcsix*enerdiff;
        real       dvdlambda  = 0;
        if (ir->efep != efepNO)
        {
            dvdlambda += (fr->avcsix[1] - fr->avcsix[0])*enerdiff;
        }
        if (bCorrAll)
        {
            const real enerdiff  = ninter*(dens*fr->enerdifftwelve - fr->enershifttwelve);
            *enercorr           += avctwelve*enerdiff;
            if (fr->efep != efepNO)
            {
                dvdlambda += (fr->avctwelve[1] - fr->avctwelve[0])*enerdiff;
            }
        }

        real svir  = 0;
        real spres = 0;
        if (bCorrPres)
        {
            svir = ninter*dens*avcsix*fr->virdiffsix/3.0;
            if (ir->eDispCorr == edispcAllEnerPres)
            {
                svir += ninter*dens*avctwelve*fr->virdifftwelve/3.0;
            }
            /* The factor 2 is because of the Gromacs virial definition */
            spres = -2.0*invvol*svir*PRESFAC;

            for (int m = 0; m < DIM; m++)
            {
                virial[m][m] += svir;
                pres[m][m]   += spres;
            }
            *prescorr += spres;
        }

        /* Can't currently control when it prints, for now, just print when degugging */
        if (debug)
        {
            if (bCorrAll)
            {
                fprintf(debug, "Long Range LJ corr.: <C6> %10.4e, <C12> %10.4e\n",
                        avcsix, avctwelve);
            }
            if (bCorrPres)
            {
                fprintf(debug,
                        "Long Range LJ corr.: Epot %10g, Pres: %10g, Vir: %10g\n",
                        *enercorr, spres, svir);
            }
            else
            {
                fprintf(debug, "Long Range LJ corr.: Epot %10g\n", *enercorr);
            }
        }

        if (fr->efep != efepNO)
        {
            *dvdlcorr += dvdlambda;
        }
    }
}
