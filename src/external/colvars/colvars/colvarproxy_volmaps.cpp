// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvarproxy_volmaps.h"
#include "colvarmodule_utils.h"


colvarproxy_volmaps::colvarproxy_volmaps()
{
  volmaps_rms_applied_force_ = volmaps_max_applied_force_ = 0.0;
}


colvarproxy_volmaps::~colvarproxy_volmaps() {}


int colvarproxy_volmaps::check_volmaps_available()
{
  return COLVARS_NOT_IMPLEMENTED;
}


int colvarproxy_volmaps::reset()
{
  int error_code = COLVARS_OK;
  for (size_t i = 0; i < volmaps_ids.size(); i++) {
    error_code |= clear_volmap(i);
  }
  volmaps_ids.clear();
  volmaps_refcount.clear();
  volmaps_values.clear();
  volmaps_new_colvar_forces.clear();
  volmaps_filenames.clear();
  return error_code;
}


int colvarproxy_volmaps::add_volmap_slot(int volmap_id)
{
  volmaps_ids.push_back(volmap_id);
  volmaps_refcount.push_back(1);
  volmaps_values.push_back(0.0);
  volmaps_new_colvar_forces.push_back(0.0);
  volmaps_filenames.push_back(std::string(""));
  return (volmaps_ids.size() - 1);
}


int colvarproxy_volmaps::request_engine_volmap_by_id(int /* volmap_id */)
{
  return -1;
}


int colvarproxy_volmaps::request_engine_volmap_by_name(std::string const & /* volmap_name */)
{
  return -1;
}


int colvarproxy_volmaps::init_internal_volmap_by_id(int /* volmap_id */)
{
  return -1;
}


int colvarproxy_volmaps::init_internal_volmap_by_name(std::string const & /* volmap_name */)
{
  return -1;
}


int colvarproxy_volmaps::load_internal_volmap_from_file(std::string const & /* volmap_filename */)
{
  return -1;
}


int colvarproxy_volmaps::clear_volmap(int index)
{
  if (index < 0 || ((size_t) index) >= volmaps_ids.size()) {
    return COLVARS_INPUT_ERROR;
  }
  volmaps_filenames[index].clear();

  if (volmaps_refcount[index] > 0) {
    volmaps_refcount[index] -= 1;
  }

  return COLVARS_OK;
}


int colvarproxy_volmaps::compute_volmap(int /* flags */,
                                        int /* volmap_id */,
                                        cvm::atom_group* ag,
                                        cvm::real * /* value */,
                                        cvm::real * /* atom_field */)
{
  return COLVARS_NOT_IMPLEMENTED;
}


void colvarproxy_volmaps::compute_rms_volmaps_applied_force()
{
  volmaps_rms_applied_force_ =
    compute_norm2_stats<decltype(volmaps_new_colvar_forces), 0, false>(volmaps_new_colvar_forces);
}


void colvarproxy_volmaps::compute_max_volmaps_applied_force()
{
  volmaps_max_applied_force_ =
    compute_norm2_stats<decltype(volmaps_new_colvar_forces), 1, false>(volmaps_new_colvar_forces);
}
