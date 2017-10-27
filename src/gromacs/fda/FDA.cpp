/*
 * FDA.cpp
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#include <limits>
#include <sstream>
#include "FDA.h"
#include "gromacs/fileio/readinp.h"
#include "gromacs/math/vec.h"
#include "gromacs/math/vectypes.h"
#include "gromacs/utility/cstringutil.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/real.h"
#include "gromacs/utility/smalloc.h"
#include "Utilities.h"

using namespace fda;

static const real HALF    = 1.0 / 2.0;
static const real THIRD   = 1.0 / 3.0;
static const real QUARTER = 0.25;

static const int FDA_GROUP_IDX_FDA1 = 0;
static const int FDA_GROUP_IDX_FDA2 = 1;
static const int FDA_GROUP_IDX_FDA12 = 2;
static const int FDA_GROUP_IDX_REST = 3;
static const int FDA_GROUP_DIM = 4;
static const int FDA_GROUP_DIM2 = FDA_GROUP_DIM * FDA_GROUP_DIM;

static int FDA_egp_flags[FDA_GROUP_DIM2] = {
    1, 0, 0, 1,
    0, 1, 0, 1,
    0, 0, 0, 1,
    1, 1, 1, 1
};

FDA::FDA(fda::FDASettings const& fda_settings)
 : fda_settings(fda_settings),
   atom_based(fda_settings.atom_based_result_type,
              fda_settings.syslen_atoms,
              fda_settings.atom_based_result_filename,
              fda_settings),
   residue_based(fda_settings.residue_based_result_type,
                 fda_settings.syslen_residues,
                 fda_settings.residue_based_result_filename,
                 fda_settings),
   time_averaging_steps(0),
   time_averaging_com(nullptr),
   nsteps(0)
{
    if (fda_settings.time_averaging_period != 1) {
        if (residue_based.PF_or_PS_mode()) {
            snew(time_averaging_com, fda_settings.syslen_residues);
            clear_rvecs(fda_settings.syslen_residues, time_averaging_com);
        }
    }
}

FDA::~FDA()
{
    atom_based.write_compat_header(nsteps);
    residue_based.write_compat_header(nsteps);
}

void FDA::add_bonded_nocheck(int i, int j, fda::InteractionType type, rvec force)
{
    // the calling functions will not have i == j, but there is not such guarantee for ri and rj;
    // and it makes no sense to look at the interaction of a residue to itself
    if (residue_based.PF_or_PS_mode()) {
        int ri = fda_settings.get_atom2residue(i);
        int rj = fda_settings.get_atom2residue(j);
        rvec force_residue;
        if (ri != rj) {
            switch(fda_settings.one_pair) {
                case fda::OnePair::DETAILED:
                    if (ri > rj) {
                        int_swap(&ri, &rj);
                        clear_rvec(force_residue);
                        rvec_dec(force_residue, force);
                        residue_based.distributed_forces.add_detailed(ri, rj, force_residue, to_pure(type));
                    } else {
                        residue_based.distributed_forces.add_detailed(ri, rj, force, to_pure(type));
                    }
                    break;
                case fda::OnePair::SUMMED:
                    if (ri > rj) {
                        int_swap(&ri, &rj);
                        clear_rvec(force_residue);
                        rvec_dec(force_residue, force);
                        residue_based.distributed_forces.add_summed(ri, rj, force_residue, type);
                    } else {
                        residue_based.distributed_forces.add_summed(ri, rj, force, type);
                    }
                    break;
            }
        }
    }

    if (atom_based.PF_or_PS_mode()) {
        if (i > j) {
            int_swap(&i, &j);
            rvec_opp(force);
        }
        switch(fda_settings.one_pair) {
            case fda::OnePair::DETAILED:
                atom_based.distributed_forces.add_detailed(i, j, force, to_pure(type));
                break;
            case fda::OnePair::SUMMED:
                atom_based.distributed_forces.add_summed(i, j, force, type);
                break;
        }
    }
}

void FDA::add_bonded(int i, int j, fda::InteractionType type, rvec force)
{
    // leave early if the interaction is not interesting
    if (!(fda_settings.type & type)) return;
    if (!fda_settings.atoms_in_groups(i, j)) return;

    add_bonded_nocheck(i, j, type, force);
}

void FDA::add_nonbonded_single(int i, int j, fda::InteractionType type, real force, real dx, real dy, real dz)
{
    // leave early if the interaction is not interesting
    if (!(fda_settings.type & type)) return;
    if (!fda_settings.atoms_in_groups(i, j)) return;

    rvec force_v;
    force_v[0] = force * dx;
    force_v[1] = force * dy;
    force_v[2] = force * dz;
    add_bonded_nocheck(i, j, type, force_v);
}

void FDA::add_nonbonded(int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz)
{
    real pf_lj_residue, pf_coul_residue, pf_lj_coul;
    rvec pf_lj_atom_v, pf_lj_residue_v, pf_coul_atom_v, pf_coul_residue_v;

    /* first check that the interaction is interesting before doing expensive calculations and atom lookup*/
    if (!(fda_settings.type & fda::InteractionType_COULOMB))
        if (!(fda_settings.type & fda::InteractionType_LJ))
            return;
        else {
            add_nonbonded_single(i, j, fda::InteractionType_LJ, pf_lj, dx, dy, dz);
            return;
        }
    else
        if (!(fda_settings.type & fda::InteractionType_LJ)) {
            add_nonbonded_single(i, j, fda::InteractionType_COULOMB, pf_coul, dx, dy, dz);
            return;
        }

    if (!fda_settings.atoms_in_groups(i, j)) return;

    /* checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2
     * however, if only residue_based_result_type is non-zero, atoms won't be initialized... so the conversion to residue numebers needs to be done here already;
     * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
     * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
     * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
     */
    if (residue_based.PF_or_PS_mode()) {
        /* the calling functions will not have i == j, but there is not such guarantee for ri and rj;
         * and it makes no sense to look at the interaction of a residue to itself
         */
        int ri = fda_settings.get_atom2residue(i);
        int rj = fda_settings.get_atom2residue(j);
        if (ri != rj) {
            if (ri > rj) {
                int tmp = rj; rj = ri; ri = tmp; // swap
                pf_lj_residue = -pf_lj;
                pf_coul_residue = -pf_coul;
            } else {
                pf_lj_residue = pf_lj;
                pf_coul_residue = pf_coul;
            }

            switch(fda_settings.one_pair) {
                case fda::OnePair::DETAILED:
                    pf_coul_residue_v[0] = pf_coul_residue * dx;
                    pf_coul_residue_v[1] = pf_coul_residue * dy;
                    pf_coul_residue_v[2] = pf_coul_residue * dz;
                    residue_based.distributed_forces.add_detailed(ri, rj, pf_coul_residue_v, fda::PureInteractionType::COULOMB);
                    pf_lj_residue_v[0] = pf_lj_residue * dx;
                    pf_lj_residue_v[1] = pf_lj_residue * dy;
                    pf_lj_residue_v[2] = pf_lj_residue * dz;
                    residue_based.distributed_forces.add_detailed(ri, rj, pf_lj_residue_v, fda::PureInteractionType::LJ);
                    break;
                case fda::OnePair::SUMMED:
                    pf_lj_coul = pf_lj_residue + pf_coul_residue;
                    pf_coul_residue_v[0] = pf_lj_coul * dx;
                    pf_coul_residue_v[1] = pf_lj_coul * dy;
                    pf_coul_residue_v[2] = pf_lj_coul * dz;
                    residue_based.distributed_forces.add_summed(ri, rj, pf_coul_residue_v, fda::InteractionType_COULOMB | fda::InteractionType_LJ);
                    break;
            }
        }
    }

    /* i & j as well as pf_lj & pf_coul are not used after this point, so it's safe to operate on their values directly */
    if (atom_based.PF_or_PS_mode()) {
        if (i > j) {
            int tmp = j; j = i; i = tmp; // swap
            pf_lj = -pf_lj;
            pf_coul = -pf_coul;
        }
        switch(fda_settings.one_pair) {
            case fda::OnePair::DETAILED:
                pf_coul_atom_v[0] = pf_coul * dx;
                pf_coul_atom_v[1] = pf_coul * dy;
                pf_coul_atom_v[2] = pf_coul * dz;
                atom_based.distributed_forces.add_detailed(i, j, pf_coul_atom_v, fda::PureInteractionType::COULOMB);
                pf_lj_atom_v[0] = pf_lj * dx;
                pf_lj_atom_v[1] = pf_lj * dy;
                pf_lj_atom_v[2] = pf_lj * dz;
                atom_based.distributed_forces.add_detailed(i, j, pf_lj_atom_v, fda::PureInteractionType::LJ);
                break;
            case fda::OnePair::SUMMED:
                pf_lj_coul = pf_lj + pf_coul;
                pf_coul_atom_v[0] = pf_lj_coul * dx;
                pf_coul_atom_v[1] = pf_lj_coul * dy;
                pf_coul_atom_v[2] = pf_lj_coul * dz;
                atom_based.distributed_forces.add_summed(i, j, pf_coul_atom_v, fda::InteractionType_COULOMB | fda::InteractionType_LJ);
                break;
        }
    }
}

void FDA::add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k)
{
    rvec uf_i, uf_j, uf_k, f_j_i, f_j_k, f_i_k;
    real nf_j_i, nf_j_k;

    // below computation can sometimes return before finishing to avoid division with very small numbers;
    // this situation can occur f.e. when all f_i, f_j, f_k and f_l are (almost) zero;
    // in this case there is no call to fda_add_bonded, no pairwise forces are recorded (which is different from recording zero forces!)
    if (norm(f_i) + norm(f_j) + norm(f_k) == 0.0) return;

    unitv(f_i, uf_i);
    unitv(f_j, uf_j);
    unitv(f_k, uf_k);
    nf_j_i = - norm(f_i) * iprod(uf_i, uf_j);
    nf_j_k = - norm(f_k) * iprod(uf_k, uf_j);
    svmul(nf_j_i, uf_j, f_j_i);
    svmul(nf_j_k, uf_j, f_j_k);
    rvec_add(f_j_i, f_i, f_i_k);
    add_bonded(aj, ai, fda::InteractionType_ANGLE, f_j_i);
    add_bonded(ai, ak, fda::InteractionType_ANGLE, f_i_k);
    add_bonded(aj, ak, fda::InteractionType_ANGLE, f_j_k);
}

void FDA::add_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_j, rvec f_k, rvec f_l)
{
    rvec f_mj, f_mk;
    rvec f_ipl, f_jpk, f_jpk_i, f_jpk_l, f_j_ipl, f_jpk_ipl;
    rvec f_j_i, f_k_i, f_j_l, f_k_l, f_j_k, f_l_i;
    rvec f_jpk_c_f_j, f_jpk_c_f_k;
    rvec uf_jpk, uf_j, uf_k;
    real cos_a, sin_a, cos_b, sin_b, sinacosbpsinbcosa;
    real nf_ipl, nf_jpk, nf_j, nf_k, nf_j_i, nf_j_l, nf_k_i, nf_k_l, nf_jpkxnf_j, nf_jpkxnf_k, nf_jpk_i, nf_jpk_l;

    // below computation can sometimes return before finishing to avoid division with very small numbers;
    // this situation can occur f.e. when all f_i, f_j, f_k and f_l are (almost) zero;
    // in this case there is no call to fda_add_bonded, no pairwise forces are recorded (which is different from recording zero forces!)
    if (norm(f_i) + norm(f_j) + norm(f_k) + norm(f_l) == 0.0) return;

    // below computation needs -f_j and -f_k
    clear_rvec(f_mj);
    rvec_sub(f_mj, f_j, f_mj);
    clear_rvec(f_mk);
    rvec_sub(f_mk, f_k, f_mk);
    rvec_add(f_i, f_l, f_ipl);
    rvec_add(f_mj, f_mk, f_jpk);
    unitv(f_jpk, uf_jpk);
    unitv(f_mj, uf_j);
    unitv(f_mk, uf_k);

    nf_ipl = norm(f_ipl);
    if (nf_ipl < GMX_FLOAT_EPS) return;

    svmul(iprod(f_i, f_ipl) / nf_ipl, uf_jpk, f_jpk_i);
    svmul(iprod(f_l, f_ipl) / nf_ipl, uf_jpk, f_jpk_l);
    rvec_add(f_jpk_i, f_jpk_l, f_jpk_ipl);

    // decompose f_jpk_i in 2 forces in directions of f_j and f_k; decompose f_jpk_l in 2 forces in directions of f_j and f_k
    nf_jpk = norm(f_jpk);
    nf_j = norm(f_mj);
    nf_k = norm(f_mk);

    // a = angle between f_jpk and -f_j, b = angle between f_jpk and -f_k
    // obtain cos from dot product and sin from cross product
    cprod(f_jpk, f_mj, f_jpk_c_f_j);
    cprod(f_jpk, f_mk, f_jpk_c_f_k);
    nf_jpkxnf_j = nf_jpk * nf_j;
    nf_jpkxnf_k = nf_jpk * nf_k;

    if ((nf_jpkxnf_j < GMX_FLOAT_EPS) || (nf_jpkxnf_k < GMX_FLOAT_EPS)) return;

    cos_a = iprod(f_jpk, f_mj) / nf_jpkxnf_j;
    sin_a = norm(f_jpk_c_f_j) / nf_jpkxnf_j;
    cos_b = iprod(f_jpk, f_mk) / nf_jpkxnf_k;
    sin_b = norm(f_jpk_c_f_k) / nf_jpkxnf_k;

    // in a triangle, known: length of one side and 2 angles; unknown: lengths of the 2 other sides
    sinacosbpsinbcosa = sin_a * cos_b + sin_b * cos_a;
    if (sinacosbpsinbcosa < GMX_FLOAT_EPS) return;

    nf_jpk_i = norm(f_jpk_i);
    nf_jpk_l = norm(f_jpk_l);
    nf_j_i = nf_jpk_i * sin_b / sinacosbpsinbcosa;
    nf_k_i = nf_jpk_i * sin_a / sinacosbpsinbcosa;
    nf_j_l = nf_jpk_l * sin_b / sinacosbpsinbcosa;
    nf_k_l = nf_jpk_l * sin_a / sinacosbpsinbcosa;

    // make vectors from lengths
    // f_j_i and f_j_l are in the direction of f_j, f_k_i and f_k_l are in the direction of f_k
    svmul(nf_j_i, uf_j, f_j_i);
    svmul(nf_j_l, uf_j, f_j_l);
    svmul(nf_k_i, uf_k, f_k_i);
    svmul(nf_k_l, uf_k, f_k_l);

    // get f_j_k from difference
    rvec_add(f_j_i, f_j_l, f_j_ipl);
    rvec_sub(f_mj, f_j_ipl, f_j_k);

    // do the same for f_k_j, just to check by comparing with f_j_k
    // f_i is minus (f_j_i + f_k_i + f_l_i) because these are forces from i on these atoms, in the opposite direction from f_i
    rvec_add(f_i, f_jpk_i, f_l_i);
    rvec_opp(f_l_i);

    add_bonded(j, i, fda::InteractionType_DIHEDRAL, f_j_i);
    add_bonded(k, i, fda::InteractionType_DIHEDRAL, f_k_i);
    add_bonded(l, i, fda::InteractionType_DIHEDRAL, f_l_i);
    add_bonded(j, k, fda::InteractionType_DIHEDRAL, f_j_k);
    add_bonded(j, l, fda::InteractionType_DIHEDRAL, f_j_l);
    add_bonded(k, l, fda::InteractionType_DIHEDRAL, f_k_l);
}

void FDA::add_virial(int ai, tensor v, real s)
{
    // Only symmetric tensor is used, therefore full multiplication is not as efficient
    // atom_vir[ai] += s * v;

    atom_based.virial_stress[ai](XX, XX) += s * v[XX][XX];
    atom_based.virial_stress[ai](YY, YY) += s * v[YY][YY];
    atom_based.virial_stress[ai](ZZ, ZZ) += s * v[ZZ][ZZ];
    atom_based.virial_stress[ai](XX, YY) += s * v[XX][YY];
    atom_based.virial_stress[ai](XX, ZZ) += s * v[XX][ZZ];
    atom_based.virial_stress[ai](YY, ZZ) += s * v[YY][ZZ];
}

void FDA::add_virial_bond(int ai, int aj, real f, real dx, real dy, real dz)
{
    if (!atom_based.VS_mode()) return;

    tensor v;
    v[XX][XX] = dx * dx * f;
    v[YY][YY] = dy * dy * f;
    v[ZZ][ZZ] = dz * dz * f;
    v[XX][YY] = dx * dy * f;
    v[XX][ZZ] = dx * dz * f;
    v[YY][ZZ] = dy * dz * f;
    add_virial(ai, v, HALF);
    add_virial(aj, v, HALF);
}

void FDA::add_virial_angle(int ai, int aj, int ak,
    rvec r_ij, rvec r_kj, rvec f_i, rvec f_k)
{
    if (!atom_based.VS_mode()) return;

    tensor v;
    v[XX][XX] = r_ij[XX] * f_i[XX] + r_kj[XX] * f_k[XX];
    v[YY][YY] = r_ij[YY] * f_i[YY] + r_kj[YY] * f_k[YY];
    v[ZZ][ZZ] = r_ij[ZZ] * f_i[ZZ] + r_kj[ZZ] * f_k[ZZ];
    v[XX][YY] = r_ij[XX] * f_i[YY] + r_kj[XX] * f_k[YY];
    v[XX][ZZ] = r_ij[XX] * f_i[ZZ] + r_kj[XX] * f_k[ZZ];
    v[YY][ZZ] = r_ij[YY] * f_i[ZZ] + r_kj[YY] * f_k[ZZ];
    add_virial(ai, v, THIRD);
    add_virial(aj, v, THIRD);
    add_virial(ak, v, THIRD);
}

void FDA::add_virial_dihedral(int i, int j, int k, int l,
    rvec f_i, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl)
{
    if (!atom_based.VS_mode()) return;

    rvec r_lj;
    tensor v;
    rvec_sub(r_kj, r_kl, r_lj);
    v[XX][XX] = r_ij[XX] * f_i[XX] + r_kj[XX] * f_k[XX] + r_lj[XX] * f_l[XX];
    v[YY][YY] = r_ij[YY] * f_i[YY] + r_kj[YY] * f_k[YY] + r_lj[YY] * f_l[YY];
    v[ZZ][ZZ] = r_ij[ZZ] * f_i[ZZ] + r_kj[ZZ] * f_k[ZZ] + r_lj[ZZ] * f_l[ZZ];
    v[XX][YY] = r_ij[XX] * f_i[YY] + r_kj[XX] * f_k[YY] + r_lj[XX] * f_l[YY];
    v[XX][ZZ] = r_ij[XX] * f_i[ZZ] + r_kj[XX] * f_k[ZZ] + r_lj[XX] * f_l[ZZ];
    v[YY][ZZ] = r_ij[YY] * f_i[ZZ] + r_kj[YY] * f_k[ZZ] + r_lj[YY] * f_l[ZZ];
    add_virial(i, v, QUARTER);
    add_virial(j, v, QUARTER);
    add_virial(k, v, QUARTER);
    add_virial(l, v, QUARTER);
}

void FDA::save_and_write_scalar_time_averages(PaddedRVecVector const& x, gmx_mtop_t *mtop)
{
    if (fda_settings.time_averaging_period != 1) {
        if (atom_based.PF_or_PS_mode())
            atom_based.distributed_forces.summed_merge_to_scalar(x);
        if (residue_based.PF_or_PS_mode()) {
        	PaddedRVecVector com = get_residues_com(x, mtop);
            residue_based.distributed_forces.summed_merge_to_scalar(com);
            for (int i = 0; i != fda_settings.syslen_residues; ++i) {
                rvec_inc(time_averaging_com[i], com[i]);
            }
        }
        ++time_averaging_steps;
        if (fda_settings.time_averaging_period != 0 and time_averaging_steps >= fda_settings.time_averaging_period)
            write_scalar_time_averages();
    } else {
        write_frame(x, mtop);
    }
    // Clear arrays for next frame
    atom_based.distributed_forces.clear();
    residue_based.distributed_forces.clear();
}

void FDA::write_scalar_time_averages()
{
    if (time_averaging_steps == 0) return;

    if (atom_based.PF_or_PS_mode()) {
        atom_based.distributed_forces.scalar_real_divide(time_averaging_steps);
        if (atom_based.compatibility_mode())
            atom_based.write_frame_scalar_compat(nsteps);
        else
            atom_based.write_frame_scalar(nsteps);
        atom_based.distributed_forces.clear_scalar();
    }

    if (residue_based.PF_or_PS_mode()) {
        residue_based.distributed_forces.scalar_real_divide(time_averaging_steps);
        //pf_x_real_div(time_averaging_com, fda_settings.syslen_residues, time_averaging_steps);
        for (size_t i = 0; i != residue_based.distributed_forces.scalar.size(); ++i)
            svdiv(time_averaging_steps, time_averaging_com[i]);
        if (residue_based.compatibility_mode())
            residue_based.write_frame_scalar_compat(nsteps);
        else
            residue_based.write_frame_scalar(nsteps);
        residue_based.distributed_forces.clear_scalar();
        clear_rvecs(fda_settings.syslen_residues, time_averaging_com);
    }

    time_averaging_steps = 0;
}

void FDA::write_frame(PaddedRVecVector const& x, gmx_mtop_t *mtop)
{
    atom_based.write_frame(x, nsteps);
    residue_based.write_frame(get_residues_com(x, mtop), nsteps);
    ++nsteps;
}

void FDA::modify_energy_group_exclusions(gmx_mtop_t *mtop, t_inputrec *inputrec) const
{
    if (!fda_settings.nonbonded_exclusion_on) {
        printf("WARNING: FDA energy group exclusion is set off.\n");
        return;
    }

    // If no nonbonded interactions are needed we simply exclude all energy groups
    if (!(fda_settings.type & fda::InteractionType_NONBONDED)) {

        int i;
        for (i = 0; i < inputrec->opts.ngener; ++i) inputrec->opts.egp_flags[i] = 1;

    // If only one energy group is defined, e.g. the automatic generated group 'rest',
    // the FDA groups can easily added to the energy groups and energy group exclusions.
    // There must only be checked that no charge group will be split.
    } else if (inputrec->opts.ngener == 1) {

        #ifdef FDA_PRINT_DEBUG_ON
            printf("FDA: No energy groups are defined.\n");
        #endif

        mtop->groups.ngrpnr[egcENER] = mtop->natoms;
        snew(mtop->groups.grpnr[egcENER],mtop->natoms);

        // First set all atoms to energy group "rest" ...
        for (int i = 0; i < mtop->natoms; ++i) {
            mtop->groups.grpnr[egcENER][i] = FDA_GROUP_IDX_REST;
        }

        // ... and then overwrite the defined atoms with FDA group 1
        for (int i = fda_settings.groups->index[fda_settings.index_group1]; i < fda_settings.groups->index[fda_settings.index_group1 + 1]; ++i) {
            mtop->groups.grpnr[egcENER][fda_settings.groups->a[i]] = FDA_GROUP_IDX_FDA1;
        }

        // ... and FDA group 2
        for (int i = fda_settings.groups->index[fda_settings.index_group2]; i < fda_settings.groups->index[fda_settings.index_group2 + 1]; ++i) {
            if (mtop->groups.grpnr[egcENER][fda_settings.groups->a[i]] == FDA_GROUP_IDX_FDA1)
                mtop->groups.grpnr[egcENER][fda_settings.groups->a[i]] = FDA_GROUP_IDX_FDA12;
            else
                mtop->groups.grpnr[egcENER][fda_settings.groups->a[i]] = FDA_GROUP_IDX_FDA2;
        }

        respect_charge_groups(mtop->groups.grpnr[egcENER],mtop);

        // Search FDA group names in mtop->groups (tpr-file)
        int mtop_g1idx = add_name_to_energygrp("FDA1", &mtop->groups);
        int mtop_g2idx = add_name_to_energygrp("FDA2", &mtop->groups);
        int mtop_g3idx = add_name_to_energygrp("FDA12", &mtop->groups);
        int mtop_rest_idx = get_index_in_energygrp("rest", &mtop->groups);

        #ifdef FDA_PRINT_DEBUG_ON
            printf("=== DEBUG === mtop_g1idx = %i\n", mtop_g1idx);
            printf("=== DEBUG === mtop_g2idx = %i\n", mtop_g2idx);
            printf("=== DEBUG === mtop_g3idx = %i\n", mtop_g3idx);
            printf("=== DEBUG === mtop_rest_idx = %i\n", mtop_rest_idx);
        #endif

        // Lookup table energy group index to group index
        mtop->groups.grps[egcENER].nr = FDA_GROUP_DIM;
        snew(mtop->groups.grps[egcENER].nm_ind, FDA_GROUP_DIM);
        mtop->groups.grps[egcENER].nm_ind[0] = mtop_g1idx;
        mtop->groups.grps[egcENER].nm_ind[1] = mtop_g2idx;
        mtop->groups.grps[egcENER].nm_ind[1] = mtop_g3idx;
        mtop->groups.grps[egcENER].nm_ind[2] = mtop_rest_idx;

        // Write egp_flags table
        inputrec->opts.ngener = FDA_GROUP_DIM;
        snew(inputrec->opts.egp_flags,FDA_GROUP_DIM2);
        for (int i = 0; i < FDA_GROUP_DIM2; ++i) inputrec->opts.egp_flags[i] = FDA_egp_flags[i];

    } else {

        #ifdef FDA_PRINT_DEBUG_ON
            printf("FDA: %i energy groups are defined.\n", inputrec->opts.ngener);
        #endif

        unsigned char* FDA_eg;
        snew(FDA_eg,mtop->natoms);

        // First set all atoms to energy group "rest" ...
        for (int i = 0; i < mtop->natoms; ++i) {
            FDA_eg[i] = FDA_GROUP_IDX_REST;
        }

        // ... and then overwrite the defined atoms with FDA group 1
        for (int i = fda_settings.groups->index[fda_settings.index_group1]; i < fda_settings.groups->index[fda_settings.index_group1 + 1]; ++i) {
            FDA_eg[fda_settings.groups->a[i]] = FDA_GROUP_IDX_FDA1;
        }

        // ... and FDA group 2
        for (int i = fda_settings.groups->index[fda_settings.index_group2]; i < fda_settings.groups->index[fda_settings.index_group2 + 1]; ++i) {
            FDA_eg[fda_settings.groups->a[i]] = FDA_GROUP_IDX_FDA2;
        }

        respect_charge_groups(FDA_eg, mtop);

        // new values with prefix c_ for combined groups
        int c_ngener = FDA_GROUP_DIM * inputrec->opts.ngener;
        unsigned char* c_eg;
        snew(c_eg, mtop->natoms);

        // Generate new energy group indices as combination of original energy groups and FDA groups
        for (int i = 0; i < mtop->natoms; ++i) {
            c_eg[i] = mtop->groups.grpnr[egcENER][i] * FDA_GROUP_DIM + FDA_eg[i];
        }

        int* c_egp_flags;
        snew(c_egp_flags,c_ngener * c_ngener);

        #ifdef FDA_PRINT_DEBUG_ON
            printf("=== DEBUG === inputrec->opts.egp_flags\n");
            print_exclusion_table(inputrec->opts.egp_flags, inputrec->opts.ngener);
        #endif

        #ifdef FDA_PRINT_DEBUG_ON
            printf("=== DEBUG === FDA_egp_flags\n");
            print_exclusion_table(FDA_egp_flags, FDA_GROUP_DIM);
        #endif

        // Build new exclusion table
        for (int i = 0; i < inputrec->opts.ngener; ++i) {
            for (int j = 0; j < inputrec->opts.ngener; ++j) {
                for (int k = 0; k < FDA_GROUP_DIM; ++k) {
                    for (int l = 0; l < FDA_GROUP_DIM; ++l) {
                        c_egp_flags[(i*FDA_GROUP_DIM+k)*c_ngener + (j*FDA_GROUP_DIM+l)] =
                            inputrec->opts.egp_flags[i*inputrec->opts.ngener+j] || FDA_egp_flags[k*FDA_GROUP_DIM+l];
                    }
                }
            }
        }

        // Determine start index of energy groups
        int startIdx = std::numeric_limits<int>::max();
        for (int i = 0; i < mtop->groups.grps[egcENER].nr; ++i) {
            startIdx = std::min(startIdx, mtop->groups.grps[egcENER].nm_ind[i]);
        }

        #ifdef FDA_PRINT_DEBUG_ON
            printf("=== DEBUG === startIdx = %i\n", startIdx);
        #endif

        // Rename existing groups
        char buffer[15];
        for (int i = 0; i < inputrec->opts.ngener; ++i) {
            sprintf(buffer, "FDA%d", i);
            char **tmp;
            snew(tmp,1);
            tmp[0] = strdup(&buffer[0]);
            mtop->groups.grpname[i+startIdx] = tmp;
        }

        // Add additional group names
        for (int i = 0; i < inputrec->opts.ngener * (FDA_GROUP_DIM - 1); ++i) {
            sprintf(buffer, "FDA%d", i + inputrec->opts.ngener);
            add_name_to_energygrp(buffer, &mtop->groups);
        }

        // Update lookup table energy group index to group index
        mtop->groups.grps[egcENER].nr = c_ngener;
        srenew(mtop->groups.grps[egcENER].nm_ind, c_ngener);
        for (int i = 0; i < c_ngener; ++i) {
            mtop->groups.grps[egcENER].nm_ind[i] = i + startIdx;
        }

        // Set new energy group array
        sfree(mtop->groups.grpnr[egcENER]);
        mtop->groups.grpnr[egcENER] = c_eg;

        // Set new egp_flags table
        inputrec->opts.ngener = c_ngener;
        sfree(inputrec->opts.egp_flags);
        inputrec->opts.egp_flags = c_egp_flags;

        sfree(FDA_eg);
    }

#ifdef FDA_PRINT_DEBUG_ON

    printf("=== DEBUG === Print final arrays:\n");

    printf("=== DEBUG === mtop->groups.ngrpname = %i\n", mtop->groups.ngrpname);

    int i;
    for (i = 0; i < mtop->groups.ngrpname; ++i) {
        printf("=== DEBUG === mtop->groups.grpname[%i] = %s\n", i, *mtop->groups.grpname[i]);
    }

    printf("=== DEBUG === mtop->groups.ngrpnr[egcENER] = %i\n", mtop->groups.ngrpnr[egcENER]);

    for (i = 0; i < mtop->natoms; ++i) {
        printf("=== DEBUG === ggrpnr(&mtop->groups, egcENER, %i) = %i\n", i, (unsigned int)ggrpnr(&mtop->groups, egcENER, i));
    }

    printf("=== DEBUG === inputrec->opts->ngener = %i\n", inputrec->opts.ngener);
    print_exclusion_table(inputrec->opts.egp_flags, inputrec->opts.ngener);

    printf("=== DEBUG === mtop->nmoltype = %i\n", mtop->nmoltype);

    int j;
    for (i = 0; i < mtop->nmoltype; ++i) {
        printf("=== DEBUG === mtop->moltype[%i].cgs.nr = %i\n", i, mtop->moltype[i].cgs.nr);

        for (j = 0; j < mtop->moltype[i].cgs.nr + 1; ++j) {
            printf("=== DEBUG === mtop->moltype[%i].cgs.index[%i] = %i\n",
                i, j, mtop->moltype[i].cgs.index[j]);
        }
    }

    printf("=== DEBUG === mtop->nmolblock = %i\n", mtop->nmolblock);

    for (i = 0; i < mtop->nmolblock; ++i) {
        printf("=== DEBUG === mtop->molblock[%i].type = %i\n", i, mtop->molblock[i].type);
        printf("=== DEBUG === mtop->molblock[%i].nmol = %i\n", i, mtop->molblock[i].nmol);
    }

#endif

}

PaddedRVecVector FDA::get_residues_com(PaddedRVecVector const& x, gmx_mtop_t *mtop) const
{
    int moltype_index, mol_index, d;
    int i, atom_index, atom_global_index, residue_global_index;
    t_atoms *atoms;
    t_atom *atom_info;
    gmx_molblock_t *mb;
    rvec r;

    std::vector<real> mass(fda_settings.syslen_residues);
    PaddedRVecVector com(fda_settings.syslen_residues);

    for (i = 0; i < fda_settings.syslen_residues; i++) {
        mass[i] = 0.0;
        for (d = 0; d < DIM; d++)
            com[i][d] = 0.0;
    }

    atom_global_index = 0;
    for (moltype_index = 0; moltype_index < mtop->nmolblock; moltype_index++) {
        mb = &mtop->molblock[moltype_index];
        for (mol_index = 0; mol_index < mb->nmol; mol_index++) {
        atoms = &mtop->moltype[mb->type].atoms;
            for(atom_index = 0; atom_index < atoms->nr; atom_index++) {
                if (fda_settings.atom_in_groups(atom_global_index)) {
                    residue_global_index = fda_settings.atom_2_residue[atom_global_index];
                    atom_info=&atoms->atom[atom_index];
                    mass[residue_global_index] += atom_info->m;
                    svmul(atom_info->m, x[atom_global_index], r);
                    rvec_inc(com[residue_global_index], r);
                }
                atom_global_index++;
            }
        }
    }

    // There might be residues with no interesting atoms, so their mass would be set to the initial 0.0;
    // float/double comparison can be tricky in general, but here it's used only to prevent division by 0
    for (i = 0; i < fda_settings.syslen_residues; ++i) {
        if (mass[i] != 0.0)
            for (d = 0; d < DIM; d++)
                com[i][d] /= mass[i];
    }

    return com;
}

int FDA::add_name_to_energygrp(char const* name, gmx_groups_t* groups) const
{
    int index = groups->ngrpname;
    if (index == 255)
        gmx_fatal(FARGS, "FDA error: Limit of energy groups (256) exceeded.");
    groups->ngrpname += 1;
    srenew(groups->grpname, groups->ngrpname);
    char **tmp1;
    snew(tmp1,1);
    tmp1[0] = strdup(&name[0]);
    groups->grpname[index] = tmp1;
    return index;
}

void FDA::respect_charge_groups(unsigned char* array, gmx_mtop_t const* mtop) const
{
    bool found_FDA1, found_FDA2, found_REST;
    int i, j, k, l, length;
    int nbMolecules = 0;
    int molTypeIndex = 0;
    int nbChargeGroups = 0;
    int *chargeGroupsIndex = NULL;
    unsigned char* p = array;

    for (i = 0; i < mtop->nmolblock; ++i) {

        nbMolecules = mtop->molblock[i].nmol;
        molTypeIndex = mtop->molblock[i].type;
        nbChargeGroups = mtop->moltype[molTypeIndex].cgs.nr;
        chargeGroupsIndex = mtop->moltype[molTypeIndex].cgs.index;

        for (j = 0; j < nbMolecules; ++j) {
            for (k = 0; k < nbChargeGroups; ++k) {
               length = chargeGroupsIndex[k+1] - chargeGroupsIndex[k];
               found_FDA1 = false; found_FDA2 = false; found_REST = false;
               for (l = 0; l < length; ++l) {
                   if (p[l] == FDA_GROUP_IDX_FDA1) found_FDA1 = true;
                   if (p[l] == FDA_GROUP_IDX_FDA2) found_FDA2 = true;
                   if (p[l] == FDA_GROUP_IDX_REST) found_REST = true;
               }
               if (found_FDA1 && found_FDA2) for (l = 0; l < length; ++l) p[l] = FDA_GROUP_IDX_FDA12;
               else if (found_FDA1 && found_REST) for (l = 0; l < length; ++l) p[l] = FDA_GROUP_IDX_FDA1;
               else if (found_FDA2 && found_REST) for (l = 0; l < length; ++l) p[l] = FDA_GROUP_IDX_FDA2;
               p += length;
            }
        }
    }
}

int FDA::get_index_in_energygrp(char const* name, gmx_groups_t const* groups) const
{
    int i, index = -1;
    for (i = 0; i < groups->ngrpname; ++i) {
        if (gmx_strcasecmp(*groups->grpname[i], name) == 0) {
            index = i;
            break;
        }
    }

    if (index == -1)
        gmx_fatal(FARGS, "FDA error: group \"rest\" not found.");
    return index;
}

void FDA::print_exclusion_table(int* egp_flags, int dim) const
{
    int i, j;
    for (i = 0; i < dim; ++i) {
        for (j = 0; j < dim; ++j) {
            printf("%i ", egp_flags[i*dim + j]);
        }
        printf("\n");
    }
}

void fda_add_nonbonded(FDA *fda, int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz)
{
    fda->add_nonbonded(i, j, pf_coul, pf_lj, dx, dy, dz);
}

void fda_add_nonbonded_coulomb(FDA *fda, int i, int j, real force, real dx, real dy, real dz)
{
    fda->add_nonbonded_single(i, j, fda::InteractionType_COULOMB, force, dx, dy, dz);
}

void fda_add_nonbonded_lj(FDA *fda, int i, int j, real force, real dx, real dy, real dz)
{
    fda->add_nonbonded_single(i, j, fda::InteractionType_LJ, force, dx, dy, dz);
}

void fda_virial_bond(FDA *fda, int ai, int aj, real f, real dx, real dy, real dz)
{
    fda->add_virial_bond(ai, aj, f, dx, dy, dz);
}
