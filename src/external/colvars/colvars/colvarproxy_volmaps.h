// -*- c++ -*-

#ifndef COLVARPROXY_VOLMAPS_H
#define COLVARPROXY_VOLMAPS_H

#include "colvarmodule.h"


/// \brief Container of grid-based objects
class colvarproxy_volmaps {

public:

  /// Contructor
  colvarproxy_volmaps();

  /// Destructor
  virtual ~colvarproxy_volmaps();

  /// Clear volumetric map data
  int reset();

  /// Test whether this implementation can use volumetric maps as CVs
  virtual int check_volmaps_available();

  /// Create a slot for a volumetric map not requested yet
  int add_volmap_slot(int volmap_id);

  /// Request this map for computation by the MD engine
  /// \param volmap_id Numeric ID used by the MD engine
  /// \returns Index of the map in the colvarproxy arrays
  virtual int request_engine_volmap_by_id(int volmap_id);

  /// Request this map for computation by the MD engine
  /// \param volmap_name Name used by the MD engine
  /// \returns Index of the map in the colvarproxy arrays
  virtual int request_engine_volmap_by_name(std::string const &volmap_name);

  /// Select a map from the MD engine for internal computation
  /// \param volmap_id Numeric ID used by the MD engine
  /// \returns Index of the map in the colvarproxy arrays
  virtual int init_internal_volmap_by_id(int volmap_id);

  /// Select a map from the MD engine for internal computation
  /// \param volmap_name Name used by the MD engine
  /// \returns Index of the map in the colvarproxy arrays
  virtual int init_internal_volmap_by_name(std::string const &filename);

  /// Load a map internally in Colvars, and select it for internal computation
  /// \param filename Name of file containing the map
  /// \returns Index of the map in the colvarproxy arrays
  virtual int load_internal_volmap_from_file(std::string const &filename);

  /// Used by the CVC destructors
  virtual int clear_volmap(int index);

  /// Get the numeric ID of the given volumetric map (for the MD program)
  inline int get_volmap_id(int index) const
  {
    return volmaps_ids[index];
  }

  /// Read the current value of the volumetric map
  inline cvm::real get_engine_volmap_value(int index) const
  {
    return volmaps_values[index];
  }

  /// Request that this force is applied to the given volumetric map by the engine
  inline void apply_engine_volmap_force(int index, cvm::real const &new_force)
  {
    volmaps_new_colvar_forces[index] += new_force;
  }

  /// Re-weigh an atomic field (e.g. a colvar) by the value of a volumetric map

  /// \param flags Combination of flags
  /// \param volmap_index Index of the map in the proxy arrays
  /// \param ag Pointer to the SOA atom group
  /// \param value Pointer to location of total to increment
  /// \param atom_field Array of atomic field values (if NULL, ones are used)
  virtual int compute_volmap(int flags,
                             int volmap_index,
                             cvm::atom_group* ag,
                             cvm::real *value,
                             cvm::real *atom_field);

  /// Flags controlling what computation is done on the map
  enum {
    volmap_flag_null = 0,
    volmap_flag_gradients = 1,
    volmap_flag_use_atom_field = (1<<8)
  };

  /// Compute the root-mean-square of the applied forces
  void compute_rms_volmaps_applied_force();

  /// Compute the maximum norm among all applied forces
  void compute_max_volmaps_applied_force();

protected:

  /// Array of numeric IDs of volumetric maps (-1 if loaded internally)
  std::vector<int>          volmaps_ids;

  /// Keep track of how many times each vol map is used by a separate colvar object
  std::vector<size_t>       volmaps_refcount;

  /// Current total values of the volmaps (when computed by the MD engine)
  std::vector<cvm::real>    volmaps_values;

  /// Forces applied from colvars, to be communicated to the MD engine
  std::vector<cvm::real>    volmaps_new_colvar_forces;

  /// Names of files containing maps (empty when maps are loaded by MD engine)
  std::vector<std::string>  volmaps_filenames;

  /// Root-mean-square of the the applied forces
  cvm::real volmaps_rms_applied_force_;

  /// Maximum norm among all applied forces
  cvm::real volmaps_max_applied_force_;
};


#endif
