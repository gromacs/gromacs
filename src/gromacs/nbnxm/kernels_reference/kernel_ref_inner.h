/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 2012- The GROMACS Authors
 * and the project initiators Erik Lindahl, Berk Hess and David van der Spoel.
 * Consult the AUTHORS/COPYING files and https://www.gromacs.org for details.
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
 * https://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at https://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out https://www.gromacs.org.
 */

/* When calculating RF or Ewald interactions we calculate the electrostatic
 * forces and energies on excluded atom pairs here in the non-bonded loops.
 */
#if defined CHECK_EXCLS && (defined CALC_COULOMB || defined LJ_EWALD)
#    define EXCL_FORCES
#endif

{
    const int cj = l_cj[cjind].cj;

    for (int i = 0; i < UNROLLI; i++)
    {
        const int ai = ci * UNROLLI + i;

        const int type_i_off = type[ai] * ntype2;

        for (int j = 0; j < UNROLLJ; j++)
        {
            real FrLJ6 = 0, FrLJ12 = 0, frLJ = 0;

            /* A multiply mask used to zero an interaction
             * when either the distance cutoff is exceeded, or
             * (if appropriate) the i and j indices are
             * unsuitable for this kind of inner loop. */
#ifdef CHECK_EXCLS
            /* A multiply mask used to zero an interaction
             * when that interaction should be excluded
             * (e.g. because of bonding). */
            const real interact = static_cast<real>((l_cj[cjind].excl >> (i * UNROLLI + j)) & 1);
#    ifndef EXCL_FORCES
            real skipmask = interact;
#    else
            real skipmask = (cj == ci_sh && j <= i) ? 0.0 : 1.0;
#    endif
#else
            constexpr real interact = 1.0;
            real           skipmask = interact;
#endif

            real gmx_unused VLJ = 0;

            const int aj = cj * UNROLLJ + j;

            const real dx = xi[i * XI_STRIDE + XX] - x[aj * X_STRIDE + XX];
            const real dy = xi[i * XI_STRIDE + YY] - x[aj * X_STRIDE + YY];
            const real dz = xi[i * XI_STRIDE + ZZ] - x[aj * X_STRIDE + ZZ];

            real rsq = dx * dx + dy * dy + dz * dz;

            /* Prepare to enforce the cut-off. */
            skipmask = (rsq >= rcut2) ? 0 : skipmask;
            /* 9 flops for r^2 + cut-off check */

            // Ensure the distances do not fall below the limit where r^-12 overflows.
            // This should never happen for normal interactions.
            rsq = std::max(rsq, c_nbnxnMinDistanceSquared);

#ifdef COUNT_PAIRS
            npair++;
#endif

            real rinv = gmx::invsqrt(rsq);
            /* 5 flops for invsqrt */

            /* Partially enforce the cut-off (and perhaps
             * exclusions) to avoid possible overflow of
             * rinvsix when computing LJ, and/or overflowing
             * the Coulomb table during lookup. */
            rinv = rinv * skipmask;

            const real rinvsq = rinv * rinv;

#ifdef ENERGY_GROUPS
            const int egpJ = nbatParams.energyGroupsPerCluster->getEnergyGroup(cj, j);
#endif

#ifdef HALF_LJ
            if (i < UNROLLI / 2)
#endif
            {
                const real c6  = nbfp[type_i_off + type[aj] * 2];
                const real c12 = nbfp[type_i_off + type[aj] * 2 + 1];

#if defined LJ_CUT || defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
                real rinvsix = interact * rinvsq * rinvsq * rinvsq;
                FrLJ6        = c6 * rinvsix;
                FrLJ12       = c12 * rinvsix * rinvsix;
                frLJ         = FrLJ12 - FrLJ6;
                /* 7 flops for r^-2 + LJ force */
#    if defined CALC_ENERGIES || defined LJ_POT_SWITCH
                VLJ = (FrLJ12 + c12 * ic->repulsion_shift.cpot) / 12
                      - (FrLJ6 + c6 * ic->dispersion_shift.cpot) / 6;
                /* 7 flops for LJ energy */
#    endif
#endif

#if defined LJ_FORCE_SWITCH || defined LJ_POT_SWITCH
                /* Force or potential switching from ic->rvdw_switch */
                real r   = rsq * rinv;
                real rsw = r - ic->rvdw_switch;
                rsw      = (rsw >= 0.0 ? rsw : 0.0);
#endif
#ifdef LJ_FORCE_SWITCH
                frLJ += -c6 * (ic->dispersion_shift.c2 + ic->dispersion_shift.c3 * rsw) * rsw * rsw * r
                        + c12 * (ic->repulsion_shift.c2 + ic->repulsion_shift.c3 * rsw) * rsw * rsw * r;
#    if defined CALC_ENERGIES
                VLJ += -c6 * (-ic->dispersion_shift.c2 / 3 - ic->dispersion_shift.c3 / 4 * rsw)
                               * rsw * rsw * rsw
                       + c12 * (-ic->repulsion_shift.c2 / 3 - ic->repulsion_shift.c3 / 4 * rsw)
                                 * rsw * rsw * rsw;
#    endif
#endif

#if defined CALC_ENERGIES || defined LJ_POT_SWITCH
                /* Masking should be done after force switching,
                 * but before potential switching.
                 */
                /* Need to zero the interaction if there should be exclusion. */
                VLJ = VLJ * interact;
#endif

#ifdef LJ_POT_SWITCH
                {
                    const real sw  = 1.0 + (swV3 + (swV4 + swV5 * rsw) * rsw) * rsw * rsw * rsw;
                    const real dsw = (swF2 + (swF3 + swF4 * rsw) * rsw) * rsw * rsw;

                    frLJ = frLJ * sw - r * VLJ * dsw;
                    VLJ *= sw;
                }
#endif

#ifdef LJ_EWALD
                {
#    ifdef LJ_EWALD_COMB_GEOM
                    const real c6grid = ljc[type[ai] * 2] * ljc[type[aj] * 2];
#    elif defined LJ_EWALD_COMB_LB
                    real c6grid = NAN;
                    {
                        /* These sigma and epsilon are scaled to give 6*C6 */
                        const real sigma   = ljc[type[ai] * 2] + ljc[type[aj] * 2];
                        const real epsilon = ljc[type[ai] * 2 + 1] * ljc[type[aj] * 2 + 1];

                        const real sigma2 = sigma * sigma;
                        c6grid            = epsilon * sigma2 * sigma2 * sigma2;
                    }
#    else
#        error "No LJ Ewald combination rule defined"
#    endif

#    ifdef CHECK_EXCLS
                    /* Recalculate rinvsix without exclusion mask */
                    const real rinvsix_nm = rinvsq * rinvsq * rinvsq;
#    else
                    const real rinvsix_nm = rinvsix;
#    endif
                    const real cr2 = lje_coeff2 * rsq;
#    if GMX_DOUBLE
                    const real expmcr2 = exp(-cr2);
#    else
                    const real expmcr2    = expf(-cr2);
#    endif
                    const real poly = 1 + cr2 + 0.5 * cr2 * cr2;

                    /* Subtract the grid force from the total LJ force */
                    frLJ += c6grid * (rinvsix_nm - expmcr2 * (rinvsix_nm * poly + lje_coeff6_6));
#    ifdef CALC_ENERGIES
                    /* Shift should only be applied to real LJ pairs */
                    const real sh_mask = lje_vc * interact;

                    VLJ += c6grid / 6 * (rinvsix_nm * (1 - expmcr2 * poly) + sh_mask);
#    endif
                }
#endif /* LJ_EWALD */

#ifdef VDW_CUTOFF_CHECK
                /* Mask for VdW cut-off shorter than Coulomb cut-off */
                {
                    real skipmask_rvdw = (rsq < rvdw2) ? 1.0 : 0.0;
                    frLJ *= skipmask_rvdw;
#    ifdef CALC_ENERGIES
                    VLJ *= skipmask_rvdw;
#    endif
                }
#else
#    if defined CALC_ENERGIES
                /* Need to zero the interaction if r >= rcut */
                VLJ = VLJ * skipmask;
                /* 1 more flop for LJ energy */
#    endif
#endif /* VDW_CUTOFF_CHECK */


#ifdef CALC_ENERGIES
#    ifdef ENERGY_GROUPS
                Vvdw[egp_sh_i[i] + egpJ] += VLJ;
#    else
                Vvdw_ci += VLJ;
                /* 1 flop for LJ energy addition */
#    endif
#endif
            }

#ifdef CALC_COULOMB
            /* Enforce the cut-off and perhaps exclusions. In
             * those cases, rinv is zero because of skipmask,
             * but fcoul and vcoul will later be non-zero (in
             * both RF and table cases) because of the
             * contributions that do not depend on rinv. These
             * contributions cannot be allowed to accumulate
             * to the force and potential, and the easiest way
             * to do this is to zero the charges in
             * advance. */
            const real qq = skipmask * qi[i] * q[aj];

#    ifdef CALC_COUL_RF
            real fcoul = qq * (interact * rinv * rinvsq - k_rf2);
            /* 4 flops for RF force */
#        ifdef CALC_ENERGIES
            real vcoul = qq * (interact * rinv + reactionFieldCoefficient * rsq - reactionFieldShift);
            /* 4 flops for RF energy */
#        endif
#    endif

#    ifdef CALC_COUL_TAB
            const real rs   = rsq * rinv * tab_coul_scale;
            const int  ri   = int(rs);
            const real frac = rs - static_cast<real>(ri);
#        if !GMX_DOUBLE
            /* fexcl = F_i + frac * (F_(i+1)-F_i) */
            const real fexcl = tab_coul_FDV0[ri * 4] + frac * tab_coul_FDV0[ri * 4 + 1];
#        else
            /* fexcl = (1-frac) * F_i + frac * F_(i+1) */
            const real fexcl = (1 - frac) * tab_coul_F[ri] + frac * tab_coul_F[ri + 1];
#        endif
            real fcoul = interact * rinvsq - fexcl;
            /* 7 flops for float 1/r-table force */
#        ifdef CALC_ENERGIES
#            if !GMX_DOUBLE
            real vcoul =
                    qq
                    * (interact * (rinv - ic->sh_ewald)
                       - (tab_coul_FDV0[ri * 4 + 2] - halfsp * frac * (tab_coul_FDV0[ri * 4] + fexcl)));
            /* 7 flops for float 1/r-table energy (8 with excls) */
#            else
            real vcoul = qq
                         * (interact * (rinv - ic->sh_ewald)
                            - (tab_coul_V[ri] - halfsp * frac * (tab_coul_F[ri] + fexcl)));
#            endif
#        endif
            fcoul *= qq * rinv;
#    endif

#    ifdef CALC_ENERGIES
#        ifdef ENERGY_GROUPS
            Vc[egp_sh_i[i] + egpJ] += vcoul;
#        else
            Vc_ci += vcoul;
            /* 1 flop for Coulomb energy addition */
#        endif
#    endif
#endif

#ifdef CALC_COULOMB
            /* 2 flops for scalar LJ+Coulomb force if !HALF_LJ || (i < UNROLLI / 2) */
#    ifdef HALF_LJ
            const real fscal = (i < UNROLLI / 2) ? frLJ * rinvsq + fcoul : fcoul;
#    else
            const real fscal = frLJ * rinvsq + fcoul;
#    endif
#else
            const real fscal = frLJ * rinvsq;
#endif
            const real fx = fscal * dx;
            const real fy = fscal * dy;
            const real fz = fscal * dz;

            /* Increment i-atom force */
            fi[i * FI_STRIDE + XX] += fx;
            fi[i * FI_STRIDE + YY] += fy;
            fi[i * FI_STRIDE + ZZ] += fz;
            /* Decrement j-atom force */
            f[aj * F_STRIDE + XX] -= fx;
            f[aj * F_STRIDE + YY] -= fy;
            f[aj * F_STRIDE + ZZ] -= fz;
            /* 9 flops for force addition */
        }
    }
}

#undef interact
#undef EXCL_FORCES
