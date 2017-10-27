/*
 * FDA.h
 *
 *  Created on: Oct 31, 2016
 *      Author: Bernd Doser, HITS gGmbH <bernd.doser@h-its.org>
 */

#ifndef SRC_GROMACS_FDA_FDA_H
#define SRC_GROMACS_FDA_FDA_H

#ifdef BUILD_WITH_FDA

#include "gromacs/math/vectypes.h"
#include "gromacs/utility/real.h"

#ifdef __cplusplus
#include <cstdio>
#include <vector>
#include "FDABase.h"
#include "FDASettings.h"
#include "gromacs/math/functions.h"
#include "gromacs/math/paddedvector.h"
#include "gromacs/mdtypes/inputrec.h"
#include "InteractionType.h"
#include "PureInteractionType.h"

class FDA {
public:

    /// Default constructor
    FDA(fda::FDASettings const& fda_settings = fda::FDASettings());

    /// Destructor
    /// Write compat footer
    ~FDA();

    /**
     * Checking is symmetrical for atoms i and j; one of them has to be from g1, the other one from g2;
     * the check below makes the atoms equivalent, make them always have the same order (i,j) and not (j,i) where i < j;
     * force is the force atom j exerts on atom i; if i and j are switched, the force goes in opposite direction
     * it's possible that i > j, but ri < rj, so the force has to be handled separately for each of them
     *
     * the logic is a bit complicated by the fact that atom_based_result_type and residue_based_result_type are independent;
     * if residue_based_result_type part is done first, the atom_based_result_type part can use i/j/force directly, without saving them
     * first in intermediate variables, as the initial values of i/j/force are no longer needed; if atom_based_result_type
     * is done first (the original code), i/j/force are needed for the later atom->residue mapping
     * and saving in intermediate variables is needed
     */
    void add_bonded_nocheck(int i, int j, fda::InteractionType type, rvec force);

    void add_bonded(int i, int j, fda::InteractionType type, rvec force);

    /**
     * Add a particular type of nonbonded interaction for the kernels where only one type of interaction is computed;
     * force is passed as scalar along with the distance vector (as dx, dy, dz) from which the vector force is
     * computed, the same way it's done in the nonbonded kernels
     */
    void add_nonbonded_single(int i, int j, fda::InteractionType type, real force, real dx, real dy, real dz);

    /**
     * Add a nonbonded interaction for kernels where both Coulomb and LJ are computed;
     * this is more efficient than calling the previous one twice because some of the tests are made only once;
     * forces are passed as scalars along with the distance vector (as dx, dy, dz) from which the vector forces are
     * computed, the same way it's done in the nonbonded kernels
     */
    void add_nonbonded(int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz);

    void add_angle(int ai, int aj, int ak, rvec f_i, rvec f_j, rvec f_k);

    void add_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_j, rvec f_k, rvec f_l);

    /**
     * The atom virial can be expressed as a 6-real tensor, as it's symmetric.
     * To avoid defining a new tensor type, the 9-real tensor is used instead.
     */
    void add_virial(int ai, tensor v, real s);

    /**
     * Origin on j, but for 2 atoms it doesn't matter.
     */
    void add_virial_bond(int ai, int aj, real fbond, real dx, real dy, real dz);

    /**
     * Translate to origin on the middle (j) atom:
     * vir = ri*Fi + rj*Fj + rk*Fk
     *         = (ri-rj)*Fi + (rk-rj)*Fk
     *         = r_ij[dim1]*f_i[dim2] + r_kj[dim1]*f_k[dim2]
     */
    void add_virial_angle(int ai, int aj, int ak, rvec r_ij, rvec r_kj, rvec f_i, rvec f_k);

    /**
     * Translate to origin on the second (j) atom:
     * vir = ri*Fi + rj*Fj + rk*Fk + rl*Fl
     *         = (ri-rj)*Fi + (rk-rj)*Fk + (rl-rj)*Fl
     *         = (ri-rj)*Fi + (rk-rj)*Fk + ((rl-rk) + (rk-rj))*Fl
     *         = r_ij[dim1]*f_i[dim2] + r_kj[dim1]*f_k[dim2] + (r_kj-r_kl)[dim1]*f_l[dim2]
     */
    void add_virial_dihedral(int i, int j, int k, int l, rvec f_i, rvec f_k, rvec f_l, rvec r_ij, rvec r_kj, rvec r_kl);

    /**
     * Main function for scalar time averages; saves data and decides when to write it out
     *
     * dealing with residues is more complicated because COMs have to be averaged over time;
     * it's not possible to average the atom positions and calculate COMs only once before
     * writing because for this fda->atoms would need to be initialized (to get the atom
     * number or to get the sys2ps mapping) which only happens when AtomBased is non-zero
     */
    void save_and_write_scalar_time_averages(PaddedRVecVector const& x, gmx_mtop_t *mtop);

    /**
     * Write scalar time averages; this is similar to pf_write_frame, except that time averages are used
     *
     * There are several cases:
     * steps == 0 - there was no frame saved at all or there are no frames saved since the last writing, so no writing needed
     * steps > 0, period == 0 - no frames written so far, but frames saved, so writing one frame
     * steps > 0, period != 0 - frames saved, so writing the last frame
     *
     * if this is called from pf_write_and_save... then steps is certainly > 0
     * period is certainly != 1, otherwise the averages functions are not called
     */
    void write_scalar_time_averages();

    void write_frame(PaddedRVecVector const& x, gmx_mtop_t *mtop);

    /// Main routine for FDA exclusions
    void modify_energy_group_exclusions(gmx_mtop_t *mtop, t_inputrec *inputrec) const;

    fda::FDASettings get_settings() const { return fda_settings; }

private:

    /**
     * Computes the COM for residues in system;
     * only the atoms for which sys_in_g is non-zero are considered, such that the COM might
     * not express the COM of the whole residue but the COM of the atoms of the residue which
     * are interesting for PF
     */
    PaddedRVecVector get_residues_com(PaddedRVecVector const& x, gmx_mtop_t *mtop) const;

    /// Append group to energy groups, returns the position index
    int add_name_to_energygrp(char const* name, gmx_groups_t* groups) const;

    /// FDA groups must not be defined over complete charge groups.
    /// This group redefine the energy group array with respect to the charge groups.
    void respect_charge_groups(unsigned char* array, gmx_mtop_t const* mtop) const;

    /// Return position of group in energy groups array, if not found exit with fatal error
    int get_index_in_energygrp(char const* name, gmx_groups_t const* groups) const;

    /// Print exclusion table as matrix
    void print_exclusion_table(int* egp_flags, int dim) const;

    /// Settings
    fda::FDASettings const& fda_settings;

    /// Atom-based operation
    fda::FDABase<fda::Atom> atom_based;

    /// Residue-based operation
    fda::FDABase<fda::Residue> residue_based;

    /// Counter for current step, incremented for every call of pf_save_and_write_scalar_averages()
    /// When it reaches time_averages_steps, data is written
    int time_averaging_steps;

    /// Averaged residue COM coordinates over steps, needed for COM calculations
    /// Only initialized when residue_based_forces is non-zero
    rvec *time_averaging_com;

    /// Current number of steps
    /// Incremented for each step written during run, also used to write the total number of steps at the end
    int nsteps;

};

#else

/// Dummy class for nonbonded kernels
struct FDA {};

#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Add a nonbonded interaction for kernels where both Coulomb and LJ are computed;
 * this is more efficient than calling the previous one twice because some of the tests are made only once;
 * forces are passed as scalars along with the distance vector (as dx, dy, dz) from which the vector forces are
 * computed, the same way it's done in the nonbonded kernels
 */
void fda_add_nonbonded(struct FDA *fda, int i, int j, real pf_coul, real pf_lj, real dx, real dy, real dz);

/**
 * Add coulomb part of nonbonded interaction
 */
void fda_add_nonbonded_coulomb(struct FDA *fda, int i, int j, real force, real dx, real dy, real dz);

/**
 * Add lennard-jones part of nonbonded interaction
 */
void fda_add_nonbonded_lj(struct FDA *fda, int i, int j, real force, real dx, real dy, real dz);

/**
 * Origin on j, but for 2 atoms it doesn't matter.
 */
void fda_virial_bond(struct FDA *fda, int ai, int aj, real fbond, real dx, real dy, real dz);

#ifdef __cplusplus
}
#endif

#endif // BUILD_WITH_FDA
#endif // SRC_GROMACS_FDA_FDA_H
