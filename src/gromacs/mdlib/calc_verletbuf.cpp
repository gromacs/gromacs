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
#include "gmxpre.h"

#include "calc_verletbuf.h"

#include <cmath>
#include <cstdlib>

#include <algorithm>

#include "gromacs/ewald/ewald_utils.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/units.h"
#include "gromacs/math/vec.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/nbnxm/nbnxm.h"
#include "gromacs/nbnxm/nbnxm_geometry.h"
#include "gromacs/nbnxm/nbnxm_simd.h"
#include "gromacs/topology/block.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/topology.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/strconvert.h"

/* The code in this file estimates a pairlist buffer length
 * given a target energy drift per atom per picosecond.
 * This is done by estimating the drift given a buffer length.
 * Ideally we would like to have a tight overestimate of the drift,
 * but that can be difficult to achieve.
 *
 * Significant approximations used:
 *
 * A uniform effective particle density is used to determine what the density
 * of particles around each particle is. When there are regions with higher
 * and lower, but non-zero, particle density in the system, this approximation
 * can slightly OVER/UNDERESTIMATE the drift, depending on the differences
 * in properties of the particles between those regions.
 *
 * Interactions don't affect particle motion. OVERESTIMATES the drift on longer
 * time scales. This approximation probably introduces the largest errors.
 *
 * Only take one constraint per particle into account: OVERESTIMATES the drift.
 *
 * For rotating constraints assume the same functional shape for time scales
 * where the constraints rotate significantly as the exact expression for
 * short time scales. OVERESTIMATES the drift on long time scales.
 *
 * For non-linear virtual sites use the mass of the lightest constructing atom
 * to determine the displacement. OVER/UNDERESTIMATES the drift, depending on
 * the geometry and masses of constructing atoms.
 *
 * Note that the formulas for normal atoms and linear virtual sites are exact,
 * apart from the first two approximations.
 *
 * Note that apart from the effect of the above approximations, the actual
 * drift of the total energy of a system can be orders of magnitude smaller
 * due to cancellation of positive and negative drift for different pairs.
 */


/* Struct for unique atom type for calculating the energy drift.
 * The atom displacement depends on mass and constraints.
 * The energy jump for given distance depend on LJ type and q.
 */
struct VerletbufAtomtype
{
    atom_nonbonded_kinetic_prop_t prop; /* non-bonded and kinetic atom prop. */
    int                           n;    /* #atoms of this type in the system */
};

// Struct for derivatives of a non-bonded interaction potential
struct pot_derivatives_t
{
    real md1; // -V' at the cutoff
    real d2;  //  V'' at the cutoff
    real md3; // -V''' at the cutoff
};

VerletbufListSetup verletbufGetListSetup(Nbnxm::KernelType nbnxnKernelType)
{
    /* Note that the current buffer estimation code only handles clusters
     * of size 1, 2 or 4, so for 4x8 or 8x8 we use the estimate for 4x4.
     */
    VerletbufListSetup listSetup;

    listSetup.cluster_size_i = Nbnxm::IClusterSizePerKernelType[nbnxnKernelType];
    listSetup.cluster_size_j = Nbnxm::JClusterSizePerKernelType[nbnxnKernelType];

    if (!Nbnxm::kernelTypeUsesSimplePairlist(nbnxnKernelType))
    {
        /* The GPU kernels (except for OpenCL) split the j-clusters in two halves */
        listSetup.cluster_size_j /= 2;
    }

    return listSetup;
}

VerletbufListSetup verletbufGetSafeListSetup(ListSetupType listType)
{
    /* When calling this function we often don't know which kernel type we
     * are going to use. We choose the kernel type with the smallest possible
     * i- and j-cluster sizes, so we potentially overestimate, but never
     * underestimate, the buffer drift.
     */
    Nbnxm::KernelType nbnxnKernelType;

    if (listType == ListSetupType::Gpu)
    {
        nbnxnKernelType = Nbnxm::KernelType::Gpu8x8x8;
    }
    else if (GMX_SIMD && listType == ListSetupType::CpuSimdWhenSupported)
    {
#ifdef GMX_NBNXN_SIMD_2XNN
        /* We use the smallest cluster size to be on the safe side */
        nbnxnKernelType = Nbnxm::KernelType::Cpu4xN_Simd_2xNN;
#else
        nbnxnKernelType = Nbnxm::KernelType::Cpu4xN_Simd_4xN;
#endif
    }
    else
    {
        nbnxnKernelType = Nbnxm::KernelType::Cpu4x4_PlainC;
    }

    return verletbufGetListSetup(nbnxnKernelType);
}

// Returns whether prop1 and prop2 are identical
static bool atom_nonbonded_kinetic_prop_equal(const atom_nonbonded_kinetic_prop_t& prop1,
                                              const atom_nonbonded_kinetic_prop_t& prop2)
{
    return (prop1.mass == prop2.mass && prop1.type == prop2.type && prop1.q == prop2.q
            && prop1.bConstr == prop2.bConstr && prop1.con_mass == prop2.con_mass
            && prop1.con_len == prop2.con_len);
}

static void addAtomtype(std::vector<VerletbufAtomtype>* att, const atom_nonbonded_kinetic_prop_t& prop, int nmol)
{
    if (prop.mass == 0)
    {
        /* Ignore massless particles */
        return;
    }

    size_t i = 0;
    while (i < att->size() && !atom_nonbonded_kinetic_prop_equal(prop, (*att)[i].prop))
    {
        i++;
    }

    if (i < att->size())
    {
        (*att)[i].n += nmol;
    }
    else
    {
        att->push_back({ prop, nmol });
    }
}

/* Returns the mass of atom atomIndex or 1 when setMassesToOne=true */
static real getMass(const t_atoms& atoms, int atomIndex, bool setMassesToOne)
{
    if (!setMassesToOne)
    {
        return atoms.atom[atomIndex].m;
    }
    else
    {
        return 1;
    }
}

// Set the masses of a vsites in vsite_m and the non-linear vsite count in n_nonlin_vsite
static void get_vsite_masses(const gmx_moltype_t&  moltype,
                             const gmx_ffparams_t& ffparams,
                             bool                  setMassesToOne,
                             gmx::ArrayRef<real>   vsite_m)
{
    int numNonlinearVsites = 0;

    /* Check for virtual sites, determine mass from constructing atoms */
    for (const auto& ilist : extractILists(moltype.ilist, IF_VSITE))
    {
        for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
        {
            const t_iparams& ip = ffparams.iparams[ilist.iatoms[i]];
            const int        a1 = ilist.iatoms[i + 1];

            if (ilist.functionType != F_VSITEN)
            {
                /* Only vsiten can have more than four
                   constructing atoms, so NRAL(ft) <= 5 */
                const int         maxj = NRAL(ilist.functionType);
                std::vector<real> cam(maxj, 0);
                GMX_ASSERT(maxj <= 5, "This code expect at most 5 atoms in a vsite");
                for (int j = 1; j < maxj; j++)
                {
                    const int aj = ilist.iatoms[i + 1 + j];
                    cam[j]       = getMass(moltype.atoms, aj, setMassesToOne);
                    if (cam[j] == 0)
                    {
                        cam[j] = vsite_m[aj];
                    }
                    /* A vsite should be constructed from normal atoms or
                     * vsites of lower complexity, which we have processed
                     * in a previous iteration.
                     */
                    GMX_ASSERT(cam[j] != 0, "We should have a non-zero mass");
                }

                switch (ilist.functionType)
                {
                    case F_VSITE2:
                        /* Exact */
                        vsite_m[a1] = (cam[1] * cam[2])
                                      / (cam[2] * gmx::square(1 - ip.vsite.a)
                                         + cam[1] * gmx::square(ip.vsite.a));
                        break;
                    case F_VSITE3:
                        /* Exact */
                        vsite_m[a1] = (cam[1] * cam[2] * cam[3])
                                      / (cam[2] * cam[3] * gmx::square(1 - ip.vsite.a - ip.vsite.b)
                                         + cam[1] * cam[3] * gmx::square(ip.vsite.a)
                                         + cam[1] * cam[2] * gmx::square(ip.vsite.b));
                        break;
                    case F_VSITEN:
                        GMX_RELEASE_ASSERT(false, "VsiteN should not end up in this code path");
                        break;
                    default:
                        /* Use the mass of the lightest constructing atom.
                         * This is an approximation.
                         * If the distance of the virtual site to the
                         * constructing atom is less than all distances
                         * between constructing atoms, this is a safe
                         * over-estimate of the displacement of the vsite.
                         * This condition holds for all H mass replacement
                         * vsite constructions, except for SP2/3 groups.
                         * In SP3 groups one H will have a F_VSITE3
                         * construction, so even there the total drift
                         * estimate shouldn't be far off.
                         */
                        vsite_m[a1] = cam[1];
                        for (int j = 2; j < maxj; j++)
                        {
                            vsite_m[a1] = std::min(vsite_m[a1], cam[j]);
                        }
                        numNonlinearVsites++;
                        break;
                }
            }
            else
            {
                /* Exact */
                real inv_mass             = 0;
                int  numConstructingAtoms = ffparams.iparams[ilist.iatoms[i]].vsiten.n;
                for (int j = 0; j < 3 * numConstructingAtoms; j += 3)
                {
                    int  aj    = ilist.iatoms[i + j + 2];
                    real coeff = ffparams.iparams[ilist.iatoms[i + j]].vsiten.a;
                    real m_aj;
                    if (moltype.atoms.atom[aj].ptype == ParticleType::VSite)
                    {
                        m_aj = vsite_m[aj];
                    }
                    else
                    {
                        m_aj = moltype.atoms.atom[aj].m;
                    }
                    if (m_aj <= 0)
                    {
                        gmx_incons("The mass of a vsiten constructing atom is <= 0");
                    }
                    inv_mass += coeff * coeff / m_aj;
                }
                vsite_m[a1] = 1 / inv_mass;
                /* Correct the loop increment of i for processes more than 1 entry */
                i += (numConstructingAtoms - 1) * ilistStride(ilist);
            }
            if (gmx_debug_at)
            {
                fprintf(debug,
                        "atom %4d %-20s mass %6.3f\n",
                        a1,
                        interaction_function[ilist.functionType].longname,
                        vsite_m[a1]);
            }
        }
    }

    if (debug && numNonlinearVsites > 0)
    {
        fprintf(debug, "The molecule type has %d non-linear virtual constructions\n", numNonlinearVsites);
    }
}

static std::vector<VerletbufAtomtype> getVerletBufferAtomtypes(const gmx_mtop_t& mtop,
                                                               const bool        setMassesToOne,
                                                               const bool        useFep)
{
    std::vector<VerletbufAtomtype> att;
    int                            ft, i, a1, a2, a3, a;

    for (const gmx_molblock_t& molblock : mtop.molblock)
    {
        int                  nmol    = molblock.nmol;
        const gmx_moltype_t& moltype = mtop.moltype[molblock.type];
        const t_atoms*       atoms   = &moltype.atoms;

        /* Check for constraints, as they affect the kinetic energy.
         * For virtual sites we need the masses and geometry of
         * the constructing atoms to determine their velocity distribution.
         * Thus we need a list of properties for all atoms which
         * we partially fill when looping over constraints.
         */
        std::vector<atom_nonbonded_kinetic_prop_t> prop(atoms->nr);

        for (ft = F_CONSTR; ft <= F_CONSTRNC; ft++)
        {
            const InteractionList& il = moltype.ilist[ft];

            for (i = 0; i < il.size(); i += 1 + NRAL(ft))
            {
                const t_iparams& ip = mtop.ffparams.iparams[il.iatoms[i]];
                /* When using free-energy perturbation constraint can be perturbed.
                 * As we can have a dynamic lambda, we might not know the constraint length.
                 * And even with fixed lambda we would here need to have the constraint lambda
                 * value. So we skip the optimization for a perturbed constraint, this results in a
                 * more conservative buffer estimate.
                 */
                if (useFep && ip.constr.dB != ip.constr.dA)
                {
                    continue;
                }
                /* Check for flexible constraints, indicated by length=0.
                 * As flexible constraints have varying length, we will not take
                 * them into account, which gives a more conservative estimate.
                 */
                if (ip.constr.dA == 0)
                {
                    continue;
                }
                GMX_RELEASE_ASSERT(ip.constr.dA > 0,
                                   "We should only have positive constraint lengths here");

                a1         = il.iatoms[i + 1];
                a2         = il.iatoms[i + 2];
                real mass1 = getMass(*atoms, a1, setMassesToOne);
                real mass2 = getMass(*atoms, a2, setMassesToOne);
                if (mass2 > prop[a1].con_mass)
                {
                    prop[a1].con_mass = mass2;
                    prop[a1].con_len  = ip.constr.dA;
                }
                if (mass1 > prop[a2].con_mass)
                {
                    prop[a2].con_mass = mass1;
                    prop[a2].con_len  = ip.constr.dA;
                }
            }
        }

        const InteractionList& il = moltype.ilist[F_SETTLE];

        for (i = 0; i < il.size(); i += 1 + NRAL(F_SETTLE))
        {
            const t_iparams* ip = &mtop.ffparams.iparams[il.iatoms[i]];

            a1 = il.iatoms[i + 1];
            a2 = il.iatoms[i + 2];
            a3 = il.iatoms[i + 3];
            /* Usually the mass of a1 (usually oxygen) is larger than a2/a3.
             * If this is not the case, we overestimate the displacement,
             * which leads to a larger buffer (ok since this is an exotic case).
             */
            prop[a1].con_mass = getMass(*atoms, a2, setMassesToOne);
            prop[a1].con_len  = ip->settle.doh;

            prop[a2].con_mass = getMass(*atoms, a1, setMassesToOne);
            prop[a2].con_len  = ip->settle.doh;

            prop[a3].con_mass = getMass(*atoms, a1, setMassesToOne);
            prop[a3].con_len  = ip->settle.doh;
        }

        std::vector<real> vsite_m(atoms->nr);
        get_vsite_masses(moltype, mtop.ffparams, setMassesToOne, vsite_m);

        for (a = 0; a < atoms->nr; a++)
        {
            if (atoms->atom[a].ptype == ParticleType::VSite)
            {
                prop[a].mass = vsite_m[a];
            }
            else
            {
                prop[a].mass = getMass(*atoms, a, setMassesToOne);
            }
            prop[a].type = atoms->atom[a].type;
            prop[a].q    = atoms->atom[a].q;
            /* We consider an atom constrained, #DOF=2, when it is
             * connected with constraints to (at least one) atom with
             * a mass of more than 0.4x its own mass. This is not a critical
             * parameter, since with roughly equal masses the unconstrained
             * and constrained displacement will not differ much (and both
             * overestimate the displacement).
             */
            prop[a].bConstr = (prop[a].con_mass > 0.4 * prop[a].mass);

            addAtomtype(&att, prop[a], nmol);
        }
    }

    if (gmx_debug_at)
    {
        for (size_t a = 0; a < att.size(); a++)
        {
            fprintf(debug,
                    "type %zu: m %5.2f t %d q %6.3f con %s con_m %5.3f con_l %5.3f n %d\n",
                    a,
                    att[a].prop.mass,
                    att[a].prop.type,
                    att[a].prop.q,
                    gmx::boolToString(att[a].prop.bConstr),
                    att[a].prop.con_mass,
                    att[a].prop.con_len,
                    att[a].n);
        }
    }

    return att;
}

/* This function computes two components of the estimate of the variance
 * in the displacement of one atom in a system of two constrained atoms.
 * Returns in sigma2_2d the variance due to rotation of the constrained
 * atom around the atom to which it constrained.
 * Returns in sigma2_3d the variance due to displacement of the COM
 * of the whole system of the two constrained atoms.
 *
 * Note that we only take a single constraint (the one to the heaviest atom)
 * into account. If an atom has multiple constraints, this will result in
 * an overestimate of the displacement, which gives a larger drift and buffer.
 */
void constrained_atom_sigma2(real kT_fac, const atom_nonbonded_kinetic_prop_t* prop, real* sigma2_2d, real* sigma2_3d)
{
    /* Here we decompose the motion of a constrained atom into two
     * components: rotation around the COM and translation of the COM.
     */

    /* Determine the variance of the arc length for the two rotational DOFs */
    real massFraction = prop->con_mass / (prop->mass + prop->con_mass);
    real sigma2_rot   = kT_fac * massFraction / prop->mass;

    /* The distance from the atom to the COM, i.e. the rotational arm */
    real comDistance = prop->con_len * massFraction;

    /* The variance relative to the arm */
    real sigma2_rel = sigma2_rot / gmx::square(comDistance);

    /* For sigma2_rel << 1 we don't notice the rotational effect and
     * we have a normal, Gaussian displacement distribution.
     * For larger sigma2_rel the displacement is much less, in fact it can
     * not exceed 2*comDistance. We can calculate MSD/arm^2 as:
     *   integral_x=0-inf distance2(x) x/sigma2_rel exp(-x^2/(2 sigma2_rel)) dx
     * where x is angular displacement and distance2(x) is the distance^2
     * between points at angle 0 and x:
     *   distance2(x) = (sin(x) - sin(0))^2 + (cos(x) - cos(0))^2
     * The limiting value of this MSD is 2, which is also the value for
     * a uniform rotation distribution that would be reached at long time.
     * The maximum is 2.5695 at sigma2_rel = 4.5119.
     * We approximate this integral with a rational polynomial with
     * coefficients from a Taylor expansion. This approximation is an
     * overestimate for all values of sigma2_rel. Its maximum value
     * of 2.6491 is reached at sigma2_rel = sqrt(45/2) = 4.7434.
     * We keep the approximation constant after that.
     * We use this approximate MSD as the variance for a Gaussian distribution.
     *
     * NOTE: For any sensible buffer tolerance this will result in a (large)
     * overestimate of the buffer size, since the Gaussian has a long tail,
     * whereas the actual distribution can not reach values larger than 2.
     */
    /* Coeffients obtained from a Taylor expansion */
    const real a = 1.0 / 3.0;
    const real b = 2.0 / 45.0;

    /* Our approximation is constant after sigma2_rel = 1/sqrt(b) */
    sigma2_rel = std::min(sigma2_rel, 1 / std::sqrt(b));

    /* Compute the approximate sigma^2 for 2D motion due to the rotation */
    *sigma2_2d =
            gmx::square(comDistance) * sigma2_rel / (1 + a * sigma2_rel + b * gmx::square(sigma2_rel));

    /* The constrained atom also moves (in 3D) with the COM of both atoms */
    *sigma2_3d = kT_fac / (prop->mass + prop->con_mass);
}

static void get_atom_sigma2(real kT_fac, const atom_nonbonded_kinetic_prop_t* prop, real* sigma2_2d, real* sigma2_3d)
{
    if (prop->bConstr)
    {
        /* Complicated constraint calculation in a separate function */
        constrained_atom_sigma2(kT_fac, prop, sigma2_2d, sigma2_3d);
    }
    else
    {
        /* Unconstrained atom: trivial */
        *sigma2_2d = 0;
        *sigma2_3d = kT_fac / prop->mass;
    }
}

static void approx_2dof(real s2, real x, real* shift, real* scale)
{
    /* A particle with 1 DOF constrained has 2 DOFs instead of 3.
     * This code is also used for particles with multiple constraints,
     * in which case we overestimate the displacement.
     * The 2DOF distribution is sqrt(pi/2)*erfc(r/(sqrt(2)*s))/(2*s).
     * We approximate this with scale*Gaussian(s,r+shift),
     * by matching the distribution value and derivative at x.
     * This is a tight overestimate for all r>=0 at any s and x.
     */
    real ex, er;

    ex = std::exp(-x * x / (2 * s2));
    er = std::erfc(x / std::sqrt(2 * s2));

    *shift = -x + std::sqrt(2 * s2 / M_PI) * ex / er;
    *scale = 0.5 * M_PI * std::exp(ex * ex / (M_PI * er * er)) * er;
}

// Returns an (over)estimate of the energy drift for a single atom pair,
// given the kinetic properties, displacement variances and list buffer.
static real energyDriftAtomPair(bool                     isConstrained_i,
                                bool                     isConstrained_j,
                                real                     s2,
                                real                     s2i_2d,
                                real                     s2j_2d,
                                real                     r_buffer,
                                const pot_derivatives_t* der)
{
    // For relatively small arguments erfc() is so small that if will be 0.0
    // when stored in a float. We set an argument limit of 8 (Erfc(8)=1e-29),
    // such that we can divide by erfc and have some space left for arithmetic.
    const real erfc_arg_max = 8.0;

    real rsh    = r_buffer;
    real sc_fac = 1.0;

    real c_exp, c_erfc;

    if (rsh * rsh > 2 * s2 * erfc_arg_max * erfc_arg_max)
    {
        // Below we calculate c_erfc = 0.5*erfc(rsh/sqrt(2*s2))
        // When rsh/sqrt(2*s2) increases, this erfc will be the first
        // result that underflows and becomes 0.0. To avoid this,
        // we set c_exp=0 and c_erfc=0 for large arguments.
        // This also avoids NaN in approx_2dof().
        // In any relevant case this has no effect on the results,
        // since c_exp < 6e-29, so the displacement is completely
        // negligible for such atom pairs (and an overestimate).
        // In nearly all use cases, there will be other atom pairs
        // that contribute much more to the total, so zeroing
        // this particular contribution has no effect at all.
        c_exp  = 0;
        c_erfc = 0;
    }
    else
    {
        /* For constraints: adapt r and scaling for the Gaussian */
        if (isConstrained_i)
        {
            real sh, sc;

            approx_2dof(s2i_2d, r_buffer * s2i_2d / s2, &sh, &sc);
            rsh += sh;
            sc_fac *= sc;
        }
        if (isConstrained_j)
        {
            real sh, sc;

            approx_2dof(s2j_2d, r_buffer * s2j_2d / s2, &sh, &sc);
            rsh += sh;
            sc_fac *= sc;
        }

        /* Exact contribution of an atom pair with Gaussian displacement
         * with sigma s to the energy drift for a potential with
         * derivative -md and second derivative dd at the cut-off.
         * The only catch is that for potentials that change sign
         * near the cut-off there could be an unlucky compensation
         * of positive and negative energy drift.
         * Such potentials are extremely rare though.
         *
         * Note that pot has unit energy*length, as the linear
         * atom density still needs to be put in.
         */
        c_exp  = std::exp(-rsh * rsh / (2 * s2)) / std::sqrt(2 * M_PI);
        c_erfc = 0.5 * std::erfc(rsh / (std::sqrt(2 * s2)));
    }
    real s    = std::sqrt(s2);
    real rsh2 = rsh * rsh;

    real pot1 = sc_fac * der->md1 / 2 * ((rsh2 + s2) * c_erfc - rsh * s * c_exp);
    real pot2 = sc_fac * der->d2 / 6 * (s * (rsh2 + 2 * s2) * c_exp - rsh * (rsh2 + 3 * s2) * c_erfc);
    real pot3 = sc_fac * der->md3 / 24
                * ((rsh2 * rsh2 + 6 * rsh2 * s2 + 3 * s2 * s2) * c_erfc - rsh * s * (rsh2 + 5 * s2) * c_exp);

    return pot1 + pot2 + pot3;
}

// Computes and returns an estimate of the energy drift for the whole system
static real energyDrift(gmx::ArrayRef<const VerletbufAtomtype> att,
                        const gmx_ffparams_t*                  ffp,
                        real                                   kT_fac,
                        const pot_derivatives_t&               ljDisp,
                        const pot_derivatives_t&               ljRep,
                        const pot_derivatives_t&               elec,
                        real                                   rlj,
                        real                                   rcoulomb,
                        real                                   rlist,
                        const int                              totNumAtoms,
                        const real                             effectiveAtomDensity)
{
    double drift_tot = 0;

    if (kT_fac == 0)
    {
        /* No atom displacements: no drift, avoid division by 0 */
        return drift_tot;
    }

    // Here add up the contribution of all atom pairs in the system to
    // (estimated) energy drift by looping over all atom type pairs.
    for (gmx::index i = 0; i < att.ssize(); i++)
    {
        // Get the thermal displacement variance for the i-atom type
        const atom_nonbonded_kinetic_prop_t* prop_i = &att[i].prop;
        real                                 s2i_2d, s2i_3d;
        get_atom_sigma2(kT_fac, prop_i, &s2i_2d, &s2i_3d);

        for (gmx::index j = i; j < att.ssize(); j++)
        {
            // Get the thermal displacement variance for the j-atom type
            const atom_nonbonded_kinetic_prop_t* prop_j = &att[j].prop;
            real                                 s2j_2d, s2j_3d;
            get_atom_sigma2(kT_fac, prop_j, &s2j_2d, &s2j_3d);

            /* Add up the up to four independent variances */
            real s2 = s2i_2d + s2i_3d + s2j_2d + s2j_3d;

            // Set -V', V'' and -V''' at the cut-off for LJ */
            real              c6  = ffp->iparams[prop_i->type * ffp->atnr + prop_j->type].lj.c6;
            real              c12 = ffp->iparams[prop_i->type * ffp->atnr + prop_j->type].lj.c12;
            pot_derivatives_t lj;
            lj.md1 = c6 * ljDisp.md1 + c12 * ljRep.md1;
            lj.d2  = c6 * ljDisp.d2 + c12 * ljRep.d2;
            lj.md3 = c6 * ljDisp.md3 + c12 * ljRep.md3;

            real pot_lj = energyDriftAtomPair(
                    prop_i->bConstr, prop_j->bConstr, s2, s2i_2d, s2j_2d, rlist - rlj, &lj);

            // Set -V' and V'' at the cut-off for Coulomb
            pot_derivatives_t elec_qq;
            elec_qq.md1 = elec.md1 * prop_i->q * prop_j->q;
            elec_qq.d2  = elec.d2 * prop_i->q * prop_j->q;
            elec_qq.md3 = 0;

            real pot_q = energyDriftAtomPair(
                    prop_i->bConstr, prop_j->bConstr, s2, s2i_2d, s2j_2d, rlist - rcoulomb, &elec_qq);

            // Note that attractive and repulsive potentials for individual
            // pairs can partially cancel.
            real pot = pot_lj + pot_q;

            /* Multiply by the number of atom pairs */
            if (j == i)
            {
                pot *= static_cast<double>(att[i].n) * (att[i].n - 1) / 2;
            }
            else
            {
                pot *= static_cast<double>(att[i].n) * att[j].n;
            }
            /* We need the line density to get the energy drift of the system.
             * The effective average r^2 is close to (rlist+sigma)^2.
             */
            pot *= 4 * M_PI * gmx::square(rlist + std::sqrt(s2)) * effectiveAtomDensity / totNumAtoms;

            /* Add the unsigned drift to avoid cancellation of errors */
            drift_tot += std::abs(pot);
        }
    }

    return drift_tot;
}

// Returns the chance that a particle in a cluster is at distance rlist
// when the cluster is at distance rlist
static real surface_frac(int cluster_size, real particle_distance, real rlist)
{
    real d, area_rel;

    if (rlist < 0.5 * particle_distance)
    {
        /* We have non overlapping spheres */
        return 1.0;
    }

    /* Half the inter-particle distance relative to rlist */
    d = 0.5 * particle_distance / rlist;

    /* Determine the area of the surface at distance rlist to the closest
     * particle, relative to surface of a sphere of radius rlist.
     * The formulas below assume close to cubic cells for the pair search grid,
     * which the pair search code tries to achieve.
     * Note that in practice particle distances will not be delta distributed,
     * but have some spread, often involving shorter distances,
     * as e.g. O-H bonds in a water molecule. Thus the estimates below will
     * usually be slightly too high and thus conservative.
     */
    switch (cluster_size)
    {
        case 1:
            /* One particle: trivial */
            area_rel = 1.0;
            break;
        case 2:
            /* Two particles: two spheres at fractional distance 2*a */
            area_rel = 1.0 + d;
            break;
        case 4:
            /* We assume a perfect, symmetric tetrahedron geometry.
             * The surface around a tetrahedron is too complex for a full
             * analytical solution, so we use a Taylor expansion.
             */
            area_rel = (1.0
                        + 1 / M_PI
                                  * (6 * std::acos(1 / std::sqrt(3)) * d
                                     + std::sqrt(3) * d * d
                                               * (1.0 + 5.0 / 18.0 * d * d + 7.0 / 45.0 * d * d * d * d
                                                  + 83.0 / 756.0 * d * d * d * d * d * d)));
            break;
        default: gmx_incons("surface_frac called with unsupported cluster_size");
    }

    return area_rel / cluster_size;
}

/* Returns the negative of the third derivative of a potential r^-p
 * with a force-switch function, evaluated at the cut-off rc.
 */
static real md3_force_switch(real p, real rswitch, real rc)
{
    /* The switched force function is:
     * p*r^-(p+1) + a*(r - rswitch)^2 + b*(r - rswitch)^3
     */
    real a, b;
    real md3_pot, md3_sw;

    a = -((p + 4) * rc - (p + 1) * rswitch) / (pow(rc, p + 2) * gmx::square(rc - rswitch));
    b = ((p + 3) * rc - (p + 1) * rswitch) / (pow(rc, p + 2) * gmx::power3(rc - rswitch));

    md3_pot = (p + 2) * (p + 1) * p * pow(rc, p + 3);
    md3_sw  = 2 * a + 6 * b * (rc - rswitch);

    return md3_pot + md3_sw;
}

// Returns the derivatives of the Van der Waals dispersion and repulsion
static std::pair<pot_derivatives_t, pot_derivatives_t> getVdwDerivatives(const t_inputrec& ir,
                                                                         const real        repPow)
{
    pot_derivatives_t ljDisp = { 0, 0, 0 };
    pot_derivatives_t ljRep  = { 0, 0, 0 };

    if (ir.vdwtype == VanDerWaalsType::Cut)
    {
        real sw_range, md3_pswf;

        switch (ir.vdw_modifier)
        {
            case InteractionModifiers::None:
            case InteractionModifiers::PotShift:
                /* Derivatives of -r^-6 and r^-reppow */
                ljDisp.md1 = -6 * std::pow(ir.rvdw, -7.0);
                ljDisp.d2  = 7 * ljDisp.md1 / ir.rvdw;
                ljDisp.md3 = 8 * ljDisp.d2 / ir.rvdw;
                ljRep.md1  = repPow * std::pow(ir.rvdw, -(repPow + 1));
                ljRep.d2   = (repPow + 1) * ljRep.md1 / ir.rvdw;
                ljRep.md3  = (repPow + 2) * ljRep.d2 / ir.rvdw;
                break;
            case InteractionModifiers::ForceSwitch:
                /* At the cut-off: V=V'=V''=0, so we use only V''' */
                ljDisp.md3 = -md3_force_switch(6.0, ir.rvdw_switch, ir.rvdw);
                ljRep.md3  = md3_force_switch(repPow, ir.rvdw_switch, ir.rvdw);
                break;
            case InteractionModifiers::PotSwitch:
                /* At the cut-off: V=V'=V''=0.
                 * V''' is given by the original potential times
                 * the third derivative of the switch function.
                 */
                sw_range = ir.rvdw - ir.rvdw_switch;
                md3_pswf = 60.0 / gmx::power3(sw_range);

                ljDisp.md3 = -std::pow(ir.rvdw, -6.0) * md3_pswf;
                ljRep.md3  = std::pow(ir.rvdw, -repPow) * md3_pswf;
                break;
            default: gmx_incons("Unimplemented VdW modifier");
        }
    }
    else if (usingLJPme(ir.vdwtype))
    {
        real b   = calc_ewaldcoeff_lj(ir.rvdw, ir.ewald_rtol_lj);
        real r   = ir.rvdw;
        real br  = b * r;
        real br2 = br * br;
        real br4 = br2 * br2;
        real br6 = br4 * br2;
        // -dV/dr of g(br)*r^-6 [where g(x) = exp(-x^2)(1+x^2+x^4/2),
        // see LJ-PME equations in manual] and r^-reppow
        ljDisp.md1 = -std::exp(-br2) * (br6 + 3.0 * br4 + 6.0 * br2 + 6.0) * std::pow(r, -7.0);
        ljRep.md1  = repPow * pow(r, -(repPow + 1));
        // The contribution of the higher derivatives is negligible
    }
    else
    {
        gmx_fatal(FARGS,
                  "Energy drift calculation is only implemented for plain cut-off Lennard-Jones "
                  "interactions");
    }

    return { ljDisp, ljRep };
}

// Returns the derivatives of the Electrostatics interaction function, including 4 pi eps0 / eps_r
static pot_derivatives_t getElecDerivatives(const t_inputrec& ir)
{
    const real elfac = gmx::c_one4PiEps0 / ir.epsilon_r;

    pot_derivatives_t elec = { 0, 0, 0 };

    if (ir.coulombtype == CoulombInteractionType::Cut || usingRF(ir.coulombtype))
    {
        real eps_rf, k_rf;

        if (ir.coulombtype == CoulombInteractionType::Cut)
        {
            eps_rf = 1;
            k_rf   = 0;
        }
        else
        {
            eps_rf = ir.epsilon_rf / ir.epsilon_r;
            if (eps_rf != 0)
            {
                k_rf = (eps_rf - ir.epsilon_r) / (gmx::power3(ir.rcoulomb) * (2 * eps_rf + ir.epsilon_r));
            }
            else
            {
                /* reactionFieldPermitivity = infinity */
                k_rf = 0.5 / gmx::power3(ir.rcoulomb);
            }
        }

        if (eps_rf > 0)
        {
            elec.md1 = elfac * (1.0 / gmx::square(ir.rcoulomb) - 2 * k_rf * ir.rcoulomb);
        }
        elec.d2 = elfac * (2.0 / gmx::power3(ir.rcoulomb) + 2 * k_rf);
    }
    else if (usingPme(ir.coulombtype) || ir.coulombtype == CoulombInteractionType::Ewald)
    {
        real b, rc, br;

        b        = calc_ewaldcoeff_q(ir.rcoulomb, ir.ewald_rtol);
        rc       = ir.rcoulomb;
        br       = b * rc;
        elec.md1 = elfac * (b * std::exp(-br * br) * M_2_SQRTPI / rc + std::erfc(br) / (rc * rc));
        elec.d2  = elfac / (rc * rc)
                  * (2 * b * (1 + br * br) * std::exp(-br * br) * M_2_SQRTPI + 2 * std::erfc(br) / rc);
    }
    else
    {
        gmx_fatal(FARGS,
                  "Energy drift calculation is only implemented for Reaction-Field and Ewald "
                  "electrostatics");
    }

    return elec;
}

/* Returns the variance of the atomic displacement over timePeriod.
 *
 * Note: When not using BD with a non-mass dependendent friction coefficient,
 *       the return value still needs to be divided by the particle mass.
 */
static real displacementVariance(const t_inputrec& ir, real temperature, real timePeriod)
{
    real kT_fac;

    if (ir.eI == IntegrationAlgorithm::BD)
    {
        /* Get the displacement distribution from the random component only.
         * With accurate integration the systematic (force) displacement
         * should be negligible (unless nstlist is extremely large, which
         * you wouldn't do anyhow).
         */
        kT_fac = 2 * gmx::c_boltz * temperature * timePeriod;
        if (ir.bd_fric > 0)
        {
            /* This is directly sigma^2 of the displacement */
            kT_fac /= ir.bd_fric;
        }
        else
        {
            /* Per group tau_t is not implemented yet, use the maximum */
            real tau_t = ir.opts.tau_t[0];
            for (int i = 1; i < ir.opts.ngtc; i++)
            {
                tau_t = std::max(tau_t, ir.opts.tau_t[i]);
            }

            kT_fac *= tau_t;
            /* This kT_fac needs to be divided by the mass to get sigma^2 */
        }
    }
    else
    {
        kT_fac = gmx::c_boltz * temperature * gmx::square(timePeriod);
    }

    return kT_fac;
}

/* Returns the largest sigma of the Gaussian displacement over all particle
 * types. This ignores constraints, so is an overestimate.
 */
static real maxSigma(real kT_fac, gmx::ArrayRef<const VerletbufAtomtype> att)
{
    GMX_ASSERT(!att.empty(), "We should have at least one type");
    real smallestMass = att[0].prop.mass;
    for (const auto& atomType : att)
    {
        smallestMass = std::min(smallestMass, atomType.prop.mass);
    }

    return 2 * std::sqrt(kT_fac / smallestMass);
}

/* Returns the density weighted over cells weighted by the atom count per cell */
static real computeEffectiveAtomDensity(gmx::ArrayRef<const gmx::RVec> coordinates,
                                        const matrix                   box,
                                        const real                     cutoff)
{
    GMX_RELEASE_ASSERT(!coordinates.empty(), "Need coordinates to compute a density");

    gmx::IVec numCells;
    gmx::RVec invCellSize;

    for (int d = 0; d < DIM; d++)
    {
        GMX_RELEASE_ASSERT(cutoff < box[d][d], "The cutoff should be smaller than the boxsize");
        numCells[d]    = int(lround(box[d][d] / cutoff));
        invCellSize[d] = numCells[d] / box[d][d];
    }

    std::vector<int> atomCount(numCells[XX] * numCells[YY] * numCells[ZZ], 0);

    for (const gmx::RVec& coord : coordinates)
    {
        gmx::IVec indices;
        for (int d = 0; d < DIM; d++)
        {
            indices[d] = (int(coord[d] * invCellSize[d]) + numCells[d]) % numCells[d];
        }
        int index = (indices[XX] * numCells[YY] + indices[YY]) * numCells[ZZ] + indices[ZZ];
        atomCount[index]++;
    }

    int64_t sumSquares = 0;
    for (int count : atomCount)
    {
        sumSquares += gmx::square(int64_t(count));
    }

    return (double(sumSquares) / coordinates.size()) * invCellSize[XX] * invCellSize[YY] * invCellSize[ZZ];
}

real computeEffectiveAtomDensity(gmx::ArrayRef<const gmx::RVec> coordinates,
                                 const matrix                   box,
                                 const real                     cutoff,
                                 MPI_Comm                       communicator)
{
    int ourMpiRank = 0;
#if GMX_MPI
    if (communicator != MPI_COMM_NULL)
    {
        MPI_Comm_rank(communicator, &ourMpiRank);
    }
#endif

    real effectiveAtomDensity;

    if (ourMpiRank == 0)
    {
        effectiveAtomDensity = computeEffectiveAtomDensity(coordinates, box, cutoff);
    }
    if (communicator != MPI_COMM_NULL)
    {
        gmx_bcast(sizeof(effectiveAtomDensity), &effectiveAtomDensity, communicator);
    }

    return effectiveAtomDensity;
}

real calcVerletBufferSize(const gmx_mtop_t&         mtop,
                          const real                effectiveAtomDensity,
                          const t_inputrec&         ir,
                          const int                 nstlist,
                          const int                 listLifetime,
                          real                      ensembleTemperature,
                          const VerletbufListSetup& listSetup)
{
    double resolution;
    char*  env;

    real particle_distance;
    real nb_clust_frac_pairs_not_in_list_at_cutoff;

    int  ib0, ib1, ib;
    real rb, rl;
    real drift;

    if (!EI_DYNAMICS(ir.eI))
    {
        gmx_incons(
                "Can only determine the Verlet buffer size for integrators that perform dynamics");
    }
    if (ir.verletbuf_tol <= 0)
    {
        gmx_incons("The Verlet buffer tolerance needs to be larger than zero");
    }

    if (ensembleTemperature < 0)
    {
        /* We use the maximum temperature with multiple T-coupl groups.
         * We could use a per particle temperature, but since particles
         * interact, this might underestimate the buffer size.
         */
        ensembleTemperature = maxReferenceTemperature(ir);

        GMX_RELEASE_ASSERT(ensembleTemperature >= 0, "Without T-coupling we should not end up here");
    }

    /* Resolution of the buffer size */
    resolution = 0.001;

    env = getenv("GMX_VERLET_BUFFER_RES");
    if (env != nullptr)
    {
        sscanf(env, "%lf", &resolution);
    }

    /* In an atom wise pair-list there would be no pairs in the list
     * beyond the pair-list cut-off.
     * However, we use a pair-list of groups vs groups of atoms.
     * For groups of 4 atoms, the parallelism of SSE instructions, only
     * 10% of the atoms pairs are not in the list just beyond the cut-off.
     * As this percentage increases slowly compared to the decrease of the
     * Gaussian displacement distribution over this range, we can simply
     * reduce the drift by this fraction.
     * For larger groups, e.g. of 8 atoms, this fraction will be lower,
     * so then buffer size will be on the conservative (large) side.
     *
     * Note that the formulas used here do not take into account
     * cancellation of errors which could occur by missing both
     * attractive and repulsive interactions.
     *
     * The only major assumption is homogeneous particle distribution.
     * For an inhomogeneous system, such as a liquid-vapor system,
     * the buffer will be underestimated. The actual energy drift
     * will be higher by the factor: local/homogeneous particle density.
     *
     * The results of this estimate have been checked againt simulations.
     * In most cases the real drift differs by less than a factor 2.
     */

    /* Worst case assumption: HCP packing of particles gives largest distance */
    particle_distance = std::cbrt(std::sqrt(2) / effectiveAtomDensity);

    /* TODO: Obtain masses through (future) integrator functionality
     *       to avoid scattering the code with (or forgetting) checks.
     */
    const bool setMassesToOne = (ir.eI == IntegrationAlgorithm::BD && ir.bd_fric > 0);
    const auto att =
            getVerletBufferAtomtypes(mtop, setMassesToOne, ir.efep != FreeEnergyPerturbationType::No);
    GMX_ASSERT(!att.empty(), "We expect at least one type");

    if (debug)
    {
        fprintf(debug, "Using an effective atom density of: %f atoms/nm^3\n", effectiveAtomDensity);
        fprintf(debug, "particle distance assuming HCP packing: %f nm\n", particle_distance);
        fprintf(debug, "energy drift atom types: %zu\n", att.size());
    }

    pot_derivatives_t ljDisp;
    pot_derivatives_t ljRep;
    std::tie(ljDisp, ljRep) = getVdwDerivatives(ir, mtop.ffparams.reppow);

    // Determine the 1st and 2nd derivative for the electostatics
    const pot_derivatives_t elec = getElecDerivatives(ir);

    /* Determine the variance of the atomic displacement
     * over list_lifetime steps: kT_fac
     * For inertial dynamics (not Brownian dynamics) the mass factor
     * is not included in kT_fac, it is added later.
     */
    const real kT_fac = displacementVariance(ir, ensembleTemperature, listLifetime * ir.delta_t);

    if (debug)
    {
        fprintf(debug, "Derivatives of non-bonded potentials at the cut-off:\n");
        fprintf(debug, "LJ disp. -V' %9.2e V'' %9.2e -V''' %9.2e\n", ljDisp.md1, ljDisp.d2, ljDisp.md3);
        fprintf(debug, "LJ rep.  -V' %9.2e V'' %9.2e -V''' %9.2e\n", ljRep.md1, ljRep.d2, ljRep.md3);
        fprintf(debug, "Electro. -V' %9.2e V'' %9.2e\n", elec.md1, elec.d2);
        fprintf(debug, "sqrt(kT_fac) %f\n", std::sqrt(kT_fac));
    }

    /* Search using bisection */
    ib0 = -1;
    /* The drift will be neglible at 5 times the max sigma */
    ib1 = static_cast<int>(5 * maxSigma(kT_fac, att) / resolution) + 1;
    while (ib1 - ib0 > 1)
    {
        ib = (ib0 + ib1) / 2;
        rb = ib * resolution;
        rl = std::max(ir.rvdw, ir.rcoulomb) + rb;

        /* Calculate the average energy drift at the last step
         * of the nstlist steps at which the pair-list is used.
         */
        drift = energyDrift(
                att, &mtop.ffparams, kT_fac, ljDisp, ljRep, elec, ir.rvdw, ir.rcoulomb, rl, mtop.natoms, effectiveAtomDensity);

        /* Correct for the fact that we are using a Ni x Nj particle pair list
         * and not a 1 x 1 particle pair list. This reduces the drift.
         */
        /* We don't have a formula for 8 (yet), use 4 which is conservative */
        nb_clust_frac_pairs_not_in_list_at_cutoff =
                surface_frac(std::min(listSetup.cluster_size_i, 4), particle_distance, rl)
                * surface_frac(std::min(listSetup.cluster_size_j, 4), particle_distance, rl);
        drift *= nb_clust_frac_pairs_not_in_list_at_cutoff;

        /* Convert the drift to drift per unit time per atom */
        drift /= nstlist * ir.delta_t * mtop.natoms;

        if (debug)
        {
            fprintf(debug,
                    "ib %3d %3d %3d rb %.3f %dx%d fac %.3f drift %.1e\n",
                    ib0,
                    ib,
                    ib1,
                    rb,
                    listSetup.cluster_size_i,
                    listSetup.cluster_size_j,
                    nb_clust_frac_pairs_not_in_list_at_cutoff,
                    drift);
        }

        if (std::abs(drift) > ir.verletbuf_tol)
        {
            ib0 = ib;
        }
        else
        {
            ib1 = ib;
        }
    }

    return std::max(ir.rvdw, ir.rcoulomb) + ib1 * resolution;
}

/* Returns the pairlist buffer size for use as a minimum buffer size
 *
 * Note that this is a rather crude estimate. It is ok for a buffer
 * set for good energy conservation or RF electrostatics. But it is
 * too small with PME and the buffer set with the default tolerance.
 */
static real minCellSizeFromPairlistBuffer(const t_inputrec& ir)
{
    return ir.rlist - std::max(ir.rvdw, ir.rcoulomb);
}

static real chanceOfAtomCrossingCell(gmx::ArrayRef<const VerletbufAtomtype> atomtypes, real kT_fac, real cellSize)
{
    /* We assume atoms are distributed uniformly over the cell width.
     * Once an atom has moved by more than the cellSize (as passed
     * as the buffer argument to energyDriftAtomPair() below),
     * the chance of crossing the boundary of the neighbor cell
     * thus increases as 1/cellSize with the additional displacement
     * on top of cellSize. We thus create a linear interaction with
     * derivative = -1/cellSize. Using this in the energyDriftAtomPair
     * function will return the chance of crossing the next boundary.
     */
    const pot_derivatives_t boundaryInteraction = { 1 / cellSize, 0, 0 };

    real chance = 0;
    for (const VerletbufAtomtype& att : atomtypes)
    {
        const atom_nonbonded_kinetic_prop_t& propAtom = att.prop;
        real                                 s2_2d;
        real                                 s2_3d;
        get_atom_sigma2(kT_fac, &propAtom, &s2_2d, &s2_3d);

        real chancePerAtom = energyDriftAtomPair(
                propAtom.bConstr, false, s2_2d + s2_3d, s2_2d, 0, cellSize, &boundaryInteraction);

        if (propAtom.bConstr)
        {
            /* energyDriftAtomPair() uses an unlimited Gaussian displacement
             * distribution for constrained atoms, whereas they can
             * actually not move more than the COM of the two constrained
             * atoms plus twice the distance from the COM.
             * Use this maximum, limited displacement when this results in
             * a smaller chance (note that this is still an overestimate).
             */
            real massFraction = propAtom.con_mass / (propAtom.mass + propAtom.con_mass);
            real comDistance  = propAtom.con_len * massFraction;

            real chanceWithMaxDistance = energyDriftAtomPair(
                    false, false, s2_3d, 0, 0, cellSize - 2 * comDistance, &boundaryInteraction);
            chancePerAtom = std::min(chancePerAtom, chanceWithMaxDistance);
        }

        /* Take into account the line density of the boundary */
        chancePerAtom /= cellSize;

        chance += att.n * chancePerAtom;
    }

    return chance;
}

/* Struct for storing constraint properties of atoms */
struct AtomConstraintProps
{
    void addConstraint(real length)
    {
        numConstraints += 1;
        sumLengths += length;
    }

    int  numConstraints = 0; /* The number of constraints of an atom */
    real sumLengths     = 0; /* The sum of constraint length over the constraints */
};

/* Constructs and returns a list of constraint properties per atom */
static std::vector<AtomConstraintProps> getAtomConstraintProps(const gmx_moltype_t&  moltype,
                                                               const gmx_ffparams_t& ffparams)
{
    const t_atoms&                   atoms = moltype.atoms;
    std::vector<AtomConstraintProps> props(atoms.nr);

    for (const auto& ilist : extractILists(moltype.ilist, IF_CONSTRAINT))
    {
        // Settles are handled separately
        if (ilist.functionType == F_SETTLE)
        {
            continue;
        }

        for (size_t i = 0; i < ilist.iatoms.size(); i += ilistStride(ilist))
        {
            int  type   = ilist.iatoms[i];
            int  a1     = ilist.iatoms[i + 1];
            int  a2     = ilist.iatoms[i + 2];
            real length = ffparams.iparams[type].constr.dA;
            props[a1].addConstraint(length);
            props[a2].addConstraint(length);
        }
    }

    return props;
}

/* Return the chance of at least one update group in a molecule crossing a cell of size cellSize */
static real chanceOfUpdateGroupCrossingCell(const gmx_moltype_t&          moltype,
                                            const gmx_ffparams_t&         ffparams,
                                            const gmx::RangePartitioning& updateGrouping,
                                            real                          kT_fac,
                                            real                          cellSize)
{
    const t_atoms& atoms = moltype.atoms;
    GMX_ASSERT(updateGrouping.fullRange().end() == atoms.nr,
               "The update groups should match the molecule type");

    const pot_derivatives_t boundaryInteraction = { 1 / cellSize, 0, 0 };

    const auto atomConstraintProps = getAtomConstraintProps(moltype, ffparams);

    real chance = 0;
    for (int group = 0; group < updateGrouping.numBlocks(); group++)
    {
        const auto& block = updateGrouping.block(group);
        /* Determine the number of atoms with constraints and the mass of the COG */
        int  numAtomsWithConstraints = 0;
        real massSum                 = 0;
        for (const int atom : block)
        {
            if (atomConstraintProps[atom].numConstraints > 0)
            {
                numAtomsWithConstraints++;
            }
            massSum += moltype.atoms.atom[atom].m;
        }
        /* Determine the maximum possible distance between the center of mass
         * and the center of geometry of the update group
         */
        real maxComCogDistance = 0;
        if (numAtomsWithConstraints == 2)
        {
            for (const int atom : block)
            {
                if (atomConstraintProps[atom].numConstraints > 0)
                {
                    GMX_ASSERT(atomConstraintProps[atom].numConstraints == 1,
                               "Two atoms should be connected by one constraint");
                    maxComCogDistance = std::abs(atoms.atom[atom].m / massSum - 0.5)
                                        * atomConstraintProps[atom].sumLengths;
                    break;
                }
            }
        }
        else if (numAtomsWithConstraints > 2)
        {
            for (const int atom : block)
            {
                if (atomConstraintProps[atom].numConstraints == numAtomsWithConstraints - 1)
                {
                    real comCogDistance = atomConstraintProps[atom].sumLengths / numAtomsWithConstraints;
                    maxComCogDistance = std::max(maxComCogDistance, comCogDistance);
                }
            }
        }
        else if (block.size() > 1)
        {
            // All normal atoms must be connected by SETTLE
            for (const int atom : block)
            {
                const auto& ilist = moltype.ilist[F_SETTLE];
                GMX_RELEASE_ASSERT(!ilist.empty(),
                                   "There should be at least one settle in this moltype");
                for (int i = 0; i < ilist.size(); i += 1 + NRAL(F_SETTLE))
                {
                    if (atom == ilist.iatoms[i + 1])
                    {
                        const t_iparams& iparams = ffparams.iparams[ilist.iatoms[i]];
                        real             dOH     = iparams.settle.doh;
                        real             dHH     = iparams.settle.dhh;
                        real             dOMidH  = std::sqrt(dOH * dOH - 0.25_real * dHH * dHH);
                        maxComCogDistance =
                                std::abs(atoms.atom[atom].m / massSum - 1.0_real / 3.0_real) * dOMidH;
                    }
                }
            }
        }
        real s2_3d = kT_fac / massSum;
        chance += energyDriftAtomPair(
                false, false, s2_3d, 0, 0, cellSize - 2 * maxComCogDistance, &boundaryInteraction);
    }

    return chance;
}

/* Return the chance of at least one update group in the system crossing a cell of size cellSize */
static real chanceOfUpdateGroupCrossingCell(const gmx_mtop_t&      mtop,
                                            PartitioningPerMoltype updateGrouping,
                                            real                   kT_fac,
                                            real                   cellSize)
{
    GMX_RELEASE_ASSERT(updateGrouping.size() == mtop.moltype.size(),
                       "The update groups should match the topology");

    real chance = 0;
    for (const gmx_molblock_t& molblock : mtop.molblock)
    {
        const gmx_moltype_t& moltype = mtop.moltype[molblock.type];
        chance += molblock.nmol
                  * chanceOfUpdateGroupCrossingCell(
                          moltype, mtop.ffparams, updateGrouping[molblock.type], kT_fac, cellSize);
    }

    return chance;
}

real minCellSizeForAtomDisplacement(const gmx_mtop_t&      mtop,
                                    const t_inputrec&      ir,
                                    PartitioningPerMoltype updateGrouping,
                                    real                   chanceRequested,
                                    const ChanceTarget     chanceTarget)
{
    if (!EI_DYNAMICS(ir.eI) || (EI_MD(ir.eI) && ir.etc == TemperatureCoupling::No))
    {
        return minCellSizeFromPairlistBuffer(ir);
    }

    /* We will compute the chance for the whole system, so if the requested
     * chance is per atom, we need to multiply the target chance
     * by the number of atoms (since chance << 1).
     */
    switch (chanceTarget)
    {
        case ChanceTarget::System: break;
        case ChanceTarget::Atom: chanceRequested *= mtop.natoms; break;
        default: GMX_RELEASE_ASSERT(false, "Unhandled ChanceTarget");
    }

    /* We use the maximum temperature with multiple T-coupl groups.
     * We could use a per particle temperature, but since particles
     * interact, this might underestimate the displacements.
     */
    const real temperature = maxReferenceTemperature(ir);

    const bool setMassesToOne = (ir.eI == IntegrationAlgorithm::BD && ir.bd_fric > 0);

    const auto atomtypes =
            getVerletBufferAtomtypes(mtop, setMassesToOne, ir.efep != FreeEnergyPerturbationType::No);

    const real kT_fac = displacementVariance(ir, temperature, ir.nstlist * ir.delta_t);

    /* Resolution of the cell size */
    real resolution = 0.001;

    /* Search using bisection, avoid 0 and start at 1 */
    int ib0 = 0;
    /* The chance will be neglible at 10 times the max sigma */
    int  ib1      = int(10 * maxSigma(kT_fac, atomtypes) / resolution) + 1;
    real cellSize = 0;
    while (ib1 - ib0 > 1)
    {
        int ib   = (ib0 + ib1) / 2;
        cellSize = ib * resolution;

        real chance;
        if (updateGrouping.empty())
        {
            chance = chanceOfAtomCrossingCell(atomtypes, kT_fac, cellSize);
        }
        else
        {
            chance = chanceOfUpdateGroupCrossingCell(mtop, updateGrouping, kT_fac, cellSize);
        }

        /* Note: chance is for every nstlist steps */
        if (chance > chanceRequested * ir.nstlist)
        {
            ib0 = ib;
        }
        else
        {
            ib1 = ib;
        }
    }

    return cellSize;
}
