/// -*- c++ -*-
/* Based on Jeff Comer's old code */
#ifndef GMX_COLVARS_COLVARPROXY_GROMACS_H
#define GMX_COLVARS_COLVARPROXY_GROMACS_H

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarproxy.h"
#include "gromacs/random/tabulatednormaldistribution.h"
#include "gromacs/random/threefry.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdlib/groupcoord.h"
#include "gromacs/mdtypes/iforceprovider.h"
#include "gromacs/mdtypes/observableshistory.h"
#include "gromacs/topology/atoms.h"
#include "colvarproxy_gromacs_version.h"

/// \brief Communication between colvars and Gromacs (implementation of
/// \link colvarproxy \endlink)
class colvarproxy_gromacs : public colvarproxy, public gmx::IForceProvider {
public:
  // GROMACS structures.
  //PBC struct
  t_pbc gmx_pbc;
  //Box
  const real (*gmx_box)[3];
  //
  t_atoms gmx_atoms;
protected:
  bool first_timestep;

  bool colvars_restart;
  std::string config_file;
  size_t restart_frequency_s;
  int64_t previous_gmx_step;
  double bias_energy;

  bool gmx_bNS; // Is this a neighbor-search step? Eventually will become unnecessary

  // GROMACS random number generation.
  gmx::DefaultRandomEngine           rng;   // gromacs random number generator
  gmx::TabulatedNormalDistribution<> normal_distribution;


  // Node-local bookkepping data
  //! Total number of Colvars atoms
  int        n_colvars_atoms = 0;
  //! Part of the atoms that are local.
  int        nat_loc = 0;
  //! Global indices of the Colvars atoms.
  int       *ind = nullptr;
  //! Local indices of the Colvars atoms.
  int       *ind_loc = nullptr;
  //! Allocation size for ind_loc.
  int        nalloc_loc = 0;
  //! Unwrapped positions for all Colvars atoms, communicated to all nodes.
  rvec      *x_colvars_unwrapped = nullptr;
  //! Shifts for all Colvars atoms, to make molecule(s) whole.
  ivec      *xa_shifts = nullptr;
  //! Extra shifts since last DD step.
  ivec      *xa_eshifts = nullptr;
  //! Old positions for all Colvars atoms on master.
  rvec      *xa_old_whole = nullptr;
  //! Position of each local atom in the collective array.
  int       *xa_ind = nullptr;
  //! Bias forces on all Colvars atoms
  rvec      *f_colvars = nullptr;
public:
  friend class cvm::atom;
  colvarproxy_gromacs();
  ~colvarproxy_gromacs();

  // Initialize colvars.
  void init(t_inputrec *gmx_inp, int64_t step, const gmx_mtop_t &mtop, ObservablesHistory* oh,
            const std::string &prefix, gmx::ArrayRef<const std::string> filenames_config,
            const std::string &filename_restart, const t_commrec *cr,
            const rvec x[], rvec **xa_old_whole_colvars_state_p, int *n_colvars_atoms_state_p);

  void dd_make_local_atoms(const t_commrec *cr);
  // Called each step before evaluating the force provider
  // Should eventually be replaced by the MDmodule interface?
  void update_data(const t_commrec *cr, int64_t const step, t_pbc const &pbc, const matrix box, bool bNS);
  /*! \brief
    * Computes forces.
    *
    * \param[in]    forceProviderInput    struct that collects input data for the force providers
    * \param[in,out] forceProviderOutput   struct that collects output data of the force providers
    */
  virtual void calculateForces(const gmx::ForceProviderInput &forceProviderInput,
                                gmx::ForceProviderOutput      *forceProviderOutput);

  // Compute virial tensor for position r and force f, and add to matrix vir
  void add_virial_term(matrix vir, rvec const f, gmx::RVec const r);

  void add_energy (cvm::real energy);
  void finish(const t_commrec *cr);

  // **************** SYSTEM-WIDE PHYSICAL QUANTITIES ****************
  cvm::real rand_gaussian();
  // **************** SIMULATION PARAMETERS ****************
  size_t restart_frequency();
  std::string restart_output_prefix();
  std::string output_prefix();
  // **************** PERIODIC BOUNDARY CONDITIONS ****************
  cvm::rvector position_distance (cvm::atom_pos const &pos1,
                                  cvm::atom_pos const &pos2) const;

  // **************** INPUT/OUTPUT ****************
  /// Print a message to the main log
  void log (std::string const &message);
  /// Print a message to the main log and let the rest of the program handle the error
  void error (std::string const &message);
  /// Print a message to the main log and exit with error code
  void fatal_error (std::string const &message);
  /// Print a message to the main log and exit normally
  void exit (std::string const &message);
  /// Request to set the units used internally by Colvars
  int set_unit_system(std::string const &units_in, bool colvars_defined);
  int backup_file (char const *filename);
  /// Read atom identifiers from a file \param filename name of
  /// the file (usually a PDB) \param atoms array to which atoms read
  /// from "filename" will be appended \param pdb_field (optiona) if
  /// "filename" is a PDB file, use this field to determine which are
  /// the atoms to be set
  int load_atoms (char const *filename,
                           std::vector<cvm::atom> &atoms,
                           std::string const &pdb_field,
                           double const pdb_field_value = 0.0);
  /// Load the coordinates for a group of atoms from a file
  /// (usually a PDB); if "pos" is already allocated, the number of its
  /// elements must match the number of atoms in "filename"
  int load_coords (char const *filename,
                            std::vector<cvm::atom_pos> &pos,
                            const std::vector<int> &indices,
                            std::string const &pdb_field,
                            double const pdb_field_value = 0.0);

  int init_atom(int atom_number);

  int check_atom_id(int atom_number);
  void update_atom_properties(int index);
};

#endif
