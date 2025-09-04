// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#ifndef COLVARMODULE_H
#define COLVARMODULE_H

#include <cstdint>

#include "colvars_version.h"

#ifndef COLVARS_DEBUG
#define COLVARS_DEBUG false
#endif

/*! \mainpage Main page
This is the Developer's documentation for the Collective Variables module (Colvars).

You can browse the class hierarchy or the list of source files.

Please note that this documentation is only supported for the master branch, and its features may differ from those in a given release of a simulation package.
 */

/// \file colvarmodule.h
/// \brief Collective variables main module
///
/// This file declares the main class for defining and manipulating
/// collective variables: there should be only one instance of this
/// class, because several variables are made static (i.e. they are
/// shared between all object instances) to be accessed from other
/// objects.

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>

class colvarparse;
class colvar;
class colvarbias;
class colvarproxy;
class colvarvalue;


/// \brief Collective variables module (main class)
///
/// Class to control the collective variables calculation.  An object
/// (usually one) of this class is spawned from the MD program,
/// containing all i/o routines and general interface.
///
/// At initialization, the colvarmodule object creates a proxy object
/// to provide a transparent interface between the MD program and the
/// child objects
class colvarmodule {

public:

  /// Get the version string (YYYY-MM-DD format)
  std::string version() const
  {
    return std::string(COLVARS_VERSION);
  }

  /// Get the version number (higher = more recent)
  int version_number() const
  {
    return version_int;
  }

  /// Get the patch version number (non-zero in patch releases of other packages)
  int patch_version_number() const
  {
    return patch_version_int;
  }

private:

  /// Integer representing the version string (allows comparisons)
  int version_int = 0;

  /// Patch version number (non-zero in patch releases of other packages)
  int patch_version_int = 4;

public:

  /// Use a 64-bit integer to store the step number
  typedef long long step_number;

  /// Defining an abstract real number allows to switch precision
  typedef  double    real;


  // Math functions

  /// Override the STL pow() with a product for n integer
  static inline real integer_power(real const &x, int const n)
  {
    // Original code: math_special.h in LAMMPS
    double yy, ww;
    if (x == 0.0) return 0.0;
    int nn = (n > 0) ? n : -n;
    ww = x;
    for (yy = 1.0; nn != 0; nn >>= 1, ww *=ww) {
      if (nn & 1) yy *= ww;
    }
    return (n > 0) ? yy : 1.0/yy;
  }

  /// Reimplemented to work around MS compiler issues
  static inline real pow(real const &x, real const &y)
  {
    return ::pow(static_cast<double>(x), static_cast<double>(y));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real floor(real const &x)
  {
    return ::floor(static_cast<double>(x));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real fabs(real const &x)
  {
    return ::fabs(static_cast<double>(x));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real sqrt(real const &x)
  {
    return ::sqrt(static_cast<double>(x));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real sin(real const &x)
  {
    return ::sin(static_cast<double>(x));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real cos(real const &x)
  {
    return ::cos(static_cast<double>(x));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real asin(real const &x)
  {
    return ::asin(static_cast<double>(x));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real acos(real const &x)
  {
    return ::acos(static_cast<double>(x));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real atan2(real const &x, real const &y)
  {
    return ::atan2(static_cast<double>(x), static_cast<double>(y));
  }

  /// Reimplemented to work around MS compiler issues
  static inline real exp(real const &x)
  {
    return ::exp(static_cast<double>(x));
  }

  /// Reimplemented to work around MS compiler issues.  Note: log() is
  /// currently defined as the text logging function, but this can be changed
  /// at a later time
  static inline real logn(real const &x)
  {
    return ::log(static_cast<double>(x));
  }

  // Forward declarations
  class rvector;
  template <class T> class vector1d;
  template <class T> class matrix2d;
  class quaternion;
  class rotation;

  class usage;
  class memory_stream;

  /// Residue identifier
  typedef int residue_id;

  /// \brief Atom position (different type name from rvector, to make
  /// possible future PBC-transparent implementations)
  typedef rvector atom_pos;

  /// \brief 3x3 matrix of real numbers
  class rmatrix;

  // allow these classes to access protected data
  class atom;
  class atom_group;
  typedef std::vector<atom>::iterator       atom_iter;
  typedef std::vector<atom>::const_iterator atom_const_iter;

  /// Module-wide error state
  /// see constants at the top of this file
private:

  static int errorCode;

public:

  static void set_error_bits(int code);

  static bool get_error_bit(int code);

  static inline int get_error()
  {
    return errorCode;
  }

  static void clear_error();

  /// Current step number
  static step_number it;
  /// Starting step number for this run
  static step_number it_restart;

  /// Return the current step number from the beginning of this run
  static inline step_number step_relative()
  {
    return it - it_restart;
  }

  /// Return the current step number from the beginning of the whole
  /// calculation
  static inline step_number step_absolute()
  {
    return it;
  }

  bool binary_restart;

  /// \brief Finite difference step size (if there is no dynamics, or
  /// if gradients need to be tested independently from the size of
  /// dt)
  static real debug_gradients_step_size;

private:

  /// Prefix for all output files for this run
  std::string cvm_output_prefix;

public:
  /// Accessor for the above
  static inline std::string &output_prefix()
  {
    colvarmodule *cv = colvarmodule::main();
    return cv->cvm_output_prefix;
  }

private:

  /// Array of collective variables
  std::vector<colvar *> colvars;

  /// Array of collective variables
  std::vector<colvar *> colvars_active;

  /// Collective variables to be calculated on different threads;
  /// colvars with multple items (e.g. multiple active CVCs) are duplicated
  std::vector<colvar *> colvars_smp;
  /// Indexes of the items to calculate for each colvar
  std::vector<int> colvars_smp_items;

  /// Array of named atom groups
  std::vector<atom_group *> named_atom_groups;
public:
  /// Register a named atom group into named_atom_groups
  void register_named_atom_group(atom_group *ag);

  /// Remove a named atom group from named_atom_groups
  void unregister_named_atom_group(atom_group *ag);

  /// Array of collective variables
  std::vector<colvar *> *variables();

  /* TODO: implement named CVCs
  /// Array of named (reusable) collective variable components
  static std::vector<cvc *>     cvcs;
  /// Named cvcs register themselves at initialization time
  inline void register_cvc(cvc *p) {
    cvcs.push_back(p);
  }
  */

  /// Collective variables with the active flag on
  std::vector<colvar *> *variables_active();

  /// Collective variables to be calculated on different threads;
  /// colvars with multple items (e.g. multiple active CVCs) are duplicated
  std::vector<colvar *> *variables_active_smp();

  /// Indexes of the items to calculate for each colvar
  std::vector<int> *variables_active_smp_items();

  /// Array of collective variable biases
  std::vector<colvarbias *> biases;

  /// Energy of built-in and scripted biases, summed per time-step
  real total_bias_energy;

private:

  /// Pointer to a map counting how many biases of each type were used
  void *num_biases_types_used_;

  /// Array of active collective variable biases
  std::vector<colvarbias *> biases_active_;

public:

  /// Array of active collective variable biases
  std::vector<colvarbias *> *biases_active();

  /// \brief Whether debug output should be enabled (compile-time option)
  static inline bool debug()
  {
    return COLVARS_DEBUG;
  }

  /// How many objects (variables and biases) are configured yet?
  size_t size() const;

  /// Constructor
  /// \param Pointer to instance of the proxy class (communicate with engine)
  colvarmodule(colvarproxy *proxy);

private:

  /// Cannot initialize the main object without a proxy
  colvarmodule();

public:

  /// Destructor
  ~colvarmodule();

  /// Set the initial step number (it is 0 otherwise); may be overridden when reading a state
  void set_initial_step(step_number it);

  /// Actual function called by the destructor
  int reset();

  /// Open a config file, load its contents, and pass it to config_string()
  /// \param config_file_name Configuration file name
  int read_config_file(char const *config_file_name);

  /// \brief Parse a config string assuming it is a complete configuration
  /// (i.e. calling all parse functions)
  int read_config_string(std::string const &conf);

  /// \brief Parse a "clean" config string (no comments)
  int parse_config(std::string &conf);

  /// Get the configuration string read so far (includes comments)
  std::string const & get_config() const;

  // Parse functions (setup internal data based on a string)

  /// Allow reading from Windows text files using using std::getline
  /// (which can still be used when the text is produced by Colvars itself)
  static std::istream & getline(std::istream &is, std::string &line);

  /// Parse the few module's global parameters
  int parse_global_params(std::string const &conf);

  /// Parse and initialize collective variables
  int parse_colvars(std::string const &conf);

  /// Run provided Tcl script
  int run_tcl_script(std::string const &filename);

  /// Parse and initialize collective variable biases
  int parse_biases(std::string const &conf);

  /// \brief Add new configuration during parsing (e.g. to implement
  /// back-compatibility); cannot be nested, i.e. conf should not contain
  /// anything that triggers another call
  int append_new_config(std::string const &conf);

  /// Signals to the module object that the configuration has changed
  void config_changed();

private:

  /// Configuration string read so far by the module (includes comments)
  std::string config_string;

  /// Auto-generated configuration during parsing (e.g. to implement
  /// back-compatibility)
  std::string extra_conf;

  /// Parse and initialize collective variable biases of a specific type
  template <class bias_type>
  int parse_biases_type(std::string const &conf, char const *keyword);

  /// Test error condition and keyword parsing
  /// on error, delete new bias
  bool check_new_bias(std::string &conf, char const *key);

  /// Initialization Tcl script, user-provided
  std::string source_Tcl_script;

public:

  /// Return how many variables are defined
  size_t num_variables() const;

  /// Return how many variables have this feature enabled
  size_t num_variables_feature(int feature_id) const;

  /// Return how many biases are defined
  size_t num_biases() const;

  /// Return how many biases have this feature enabled
  size_t num_biases_feature(int feature_id) const;

  /// Return how many biases of this type are defined
  size_t num_biases_type(std::string const &type) const;

  /// Return the names of time-dependent biases with forces enabled (ABF,
  /// metadynamics, etc)
  std::vector<std::string> const time_dependent_biases() const;

private:
  /// Useful wrapper to interrupt parsing if any error occurs
  int catch_input_errors(int result);

public:

  // "Setup" functions (change internal data based on related data
  // from the proxy that may change during program execution)
  // No additional parsing is done within these functions

  /// (Re)initialize any internal data affected by changes in the engine
  /// Also calls setup() member functions of colvars and biases
  int update_engine_parameters();

  /// (Re)initialize and (re)read the input state file calling read_restart()
  int setup_input();

  /// (Re)initialize the output trajectory and state file (does not write it yet)
  int setup_output();

private:

  template <typename IST> IST & read_state_template_(IST &is);

  /// Default input state file; if given, it is read unless the MD engine provides it
  std::string default_input_state_file_;

  /// Internal state buffer, to be read as an unformatted stream
  std::vector<unsigned char> input_state_buffer_;

public:

  /// Read all objects' state fron a formatted (text) stream
  std::istream & read_state(std::istream &is);

  /// Read all objects' state fron an unformatted (binary) stream
  memory_stream & read_state(memory_stream &is);

  /// Set an internal state buffer, to be read later as an unformatted stream when ready
  int set_input_state_buffer(size_t n, unsigned char *buf);

  /// Same as set_input_state_buffer() for C array, but uses std::move
  int set_input_state_buffer(std::vector<unsigned char> &buf);

  /// Read the states of individual objects; allows for changes
  std::istream & read_objects_state(std::istream &is);

  /// Read the states of individual objects; allows for changes
  memory_stream & read_objects_state(memory_stream &is);

  /// If needed (old restart file), print the warning that cannot be ignored
  int print_total_forces_errning(bool warn_total_forces);

private:
  template <typename OST> OST &write_state_template_(OST &os);

public:

  /// Write the state of the module to a formatted (text) file
  std::ostream & write_state(std::ostream &os);

  /// Write the state of the module to an unformatted (binary) file
  memory_stream & write_state(memory_stream &os);

  /// Write the state of the module to an array of bytes (wrapped as a memory_stream object)
  int write_state_buffer(std::vector<unsigned char> &buffer);

  /// Strips .colvars.state from filename and checks that it is not empty
  static std::string state_file_prefix(char const *filename);

  /// Open a trajectory file if requested (and leave it open)
  int open_traj_file(std::string const &file_name);
  /// Close it (note: currently unused)
  int close_traj_file();
  /// Write in the trajectory file
  std::ostream & write_traj(std::ostream &os);
  /// Write explanatory labels in the trajectory file
  std::ostream & write_traj_label(std::ostream &os);

  /// Write all trajectory files
  int write_traj_files();
  /// Write a state file useful to resume the simulation
  int write_restart_file(std::string const &out_name);
  /// Write all other output files
  int write_output_files();
  /// Backup a file before writing it
  static int backup_file(char const *filename);

  /// Write the state into a string
  int write_restart_string(std::string &output);

  /// Look up a bias by name; returns NULL if not found
  static colvarbias * bias_by_name(std::string const &name);

  /// Look up a colvar by name; returns NULL if not found
  static colvar * colvar_by_name(std::string const &name);

  /// Look up a named atom group by name; returns NULL if not found
  static atom_group * atom_group_by_name(std::string const &name);

  /// Load new configuration for the given bias -
  /// currently works for harmonic (force constant and/or centers)
  int change_configuration(std::string const &bias_name, std::string const &conf);

  /// Read a colvar value
  std::string read_colvar(std::string const &name);

  /// Calculate change in energy from using alt. config. for the given bias -
  /// currently works for harmonic (force constant and/or centers)
  real energy_difference(std::string const &bias_name, std::string const &conf);

  /// Main worker function
  int calc();

  /// Calculate collective variables
  int calc_colvars();

  /// Calculate biases
  int calc_biases();

  /// Integrate bias and restraint forces, send colvar forces to atoms
  int update_colvar_forces();

  /// Perform analysis
  int analyze();

  /// Carry out operations needed before next step is run
  int end_of_step();

  /// \brief Read a collective variable trajectory (post-processing
  /// only, not called at runtime)
  int read_traj(char const *traj_filename,
                long        traj_read_begin,
                long        traj_read_end);

  /// Convert to string for output purposes
  static std::string to_str(char const *s);

  /// Convert to string for output purposes
  static std::string to_str(std::string const &s);

  /// Convert to string for output purposes
  static std::string to_str(bool x);

  /// Convert to string for output purposes
  static std::string to_str(int const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(size_t const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(long int const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(step_number const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(real const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(rvector const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(quaternion const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(colvarvalue const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(vector1d<real> const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(matrix2d<real> const &x,
                            size_t width = 0, size_t prec = 0);


  /// Convert to string for output purposes
  static std::string to_str(std::vector<int> const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(std::vector<size_t> const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(std::vector<long int> const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(std::vector<real> const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(std::vector<rvector> const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(std::vector<quaternion> const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(std::vector<colvarvalue> const &x,
                            size_t width = 0, size_t prec = 0);

  /// Convert to string for output purposes
  static std::string to_str(std::vector<std::string> const &x,
                            size_t width = 0, size_t prec = 0);


  /// Reduce the number of characters in a string
  static std::string wrap_string(std::string const &s,
                                 size_t nchars);

  /// Number of characters to represent a time step
  static size_t const it_width;
  /// Number of digits to represent a collective variables value(s)
  static size_t const cv_prec;
  /// Number of characters to represent a collective variables value(s)
  static size_t const cv_width;
  /// Number of digits to represent the collective variables energy
  static size_t const en_prec;
  /// Number of characters to represent the collective variables energy
  static size_t const en_width;
  /// Line separator in the log output
  static const char * const line_marker;


  // proxy functions

  /// \brief Time step of MD integrator (fs)
  static real dt();

  /// Request calculation of total force from MD engine
  static void request_total_force();

  /// Track usage of the given Colvars feature
  int cite_feature(std::string const &feature);

  /// Report usage of the Colvars features
  std::string feature_report(int flag = 0);

  /// Print a message to the main log
  /// \param message Message to print
  /// \param min_log_level Only print if cvm::log_level() >= min_log_level
  static void log(std::string const &message, int min_log_level = 10);

  /// Print a message to the main log and set global error code
  static int error(std::string const &message, int code = -1);

private:

  /// Level of logging requested by the user
  static int log_level_;

public:

  /// Level of logging requested by the user
  static inline int log_level()
  {
    return log_level_;
  }

  /// Level at which initialization messages are logged
  static inline int log_init_messages()
  {
    return 1;
  }

  /// Level at which a keyword's user-provided value is logged
  static inline int log_user_params()
  {
    return 2;
  }

  /// Level at which a keyword's default value is logged
  static inline int log_default_params()
  {
    return 3;
  }

  /// Level at which output-file operations are logged
  static inline int log_output_files()
  {
    return 4;
  }

  /// Level at which input-file operations (configuration, state) are logged
  static inline int log_input_files()
  {
    return 5;
  }

  /// \brief Get the distance between two atomic positions with pbcs handled
  /// correctly
  static rvector position_distance(atom_pos const &pos1,
                                   atom_pos const &pos2);

  /// \brief Names of .ndx files that have been loaded
  std::vector<std::string> index_file_names;

  /// \brief Names of groups from one or more Gromacs .ndx files
  std::vector<std::string> index_group_names;

  /// \brief Groups from one or more Gromacs .ndx files
  std::vector<std::vector<int> *> index_groups;

  /// \brief Read a Gromacs .ndx file
  int read_index_file(char const *filename);

  /// Clear the index groups loaded so far
  int reset_index_groups();

  /// \brief Load coordinates for a group of atoms from a file (PDB or XYZ);
  /// if "pos" is already allocated, the number of its elements must match the
  /// number of entries in "filename" \param filename name of the file \param
  /// pos array of coordinates \param atoms group containing the atoms (used
  /// to obtain internal IDs) \param pdb_field (optional) if the file is a PDB
  /// and this string is non-empty, select atoms for which this field is
  /// non-zero \param pdb_field_value (optional) if non-zero, select only
  /// atoms whose pdb_field equals this
  static int load_coords(char const *filename,
                         std::vector<rvector> *pos,
                         atom_group *atoms,
                         std::string const &pdb_field,
                         double pdb_field_value = 0.0);

  /// Load coordinates into an atom group from an XYZ file (assumes Angstroms)
  int load_coords_xyz(char const *filename,
                      std::vector<rvector> *pos,
                      atom_group *atoms,
                      bool keep_open = false);

  /// Frequency for collective variables trajectory output
  static size_t cv_traj_freq;

  /// Frequency for saving output restarts
  static size_t restart_out_freq;
  /// Output restart file name
  std::string   restart_out_name;

  /// Pseudo-random number with Gaussian distribution
  static real rand_gaussian();

protected:

  /// Configuration file parser object
  colvarparse *parse;

  /// Name of the trajectory file
  std::string cv_traj_name;

  /// Write labels at the next iteration
  bool cv_traj_write_labels;

  /// Version of the most recent state file read
  std::string restart_version_str;

  /// Integer version of the most recent state file read
  int restart_version_int;

  /// Counter for the current depth in the object hierarchy (useg e.g. in output)
  size_t depth_s;

  /// Thread-specific depth
  std::vector<size_t> depth_v;

  /// Track how many times the XYZ reader has been used
  int xyz_reader_use_count;

  /// Track usage of Colvars features
  usage *usage_;

public:

  /// Version of the most recent state file read
  inline std::string restart_version() const
  {
    return restart_version_str;
  }

  /// Integer version of the most recent state file read
  inline int restart_version_number() const
  {
    return restart_version_int;
  }

  /// Get the current object depth in the hierarchy
  static size_t & depth();

  /// Increase the depth (number of indentations in the output)
  static void increase_depth();

  /// Decrease the depth (number of indentations in the output)
  static void decrease_depth();

  static inline bool scripted_forces()
  {
    return use_scripted_forces;
  }

  /// Use scripted colvars forces?
  static bool use_scripted_forces;

  /// Wait for all biases before calculating scripted forces?
  static bool scripting_after_biases;

  /// Calculate the energy and forces of scripted biases
  int calc_scripted_forces();

  /// \brief Pointer to the proxy object, used to retrieve atomic data
  /// from the hosting program; it is static in order to be accessible
  /// from static functions in the colvarmodule class
  static colvarproxy *proxy;

  /// \brief Access the one instance of the Colvars module
  static colvarmodule *main();

};


/// Shorthand for the frequently used type prefix
typedef colvarmodule cvm;


std::ostream & operator << (std::ostream &os, cvm::rvector const &v);
std::istream & operator >> (std::istream &is, cvm::rvector &v);


namespace {
  constexpr int32_t COLVARS_OK = 0;
  constexpr int32_t COLVARS_ERROR = 1;
  constexpr int32_t COLVARS_NOT_IMPLEMENTED = (1<<1);
  constexpr int32_t COLVARS_INPUT_ERROR     = (1<<2); // out of bounds or inconsistent input
  constexpr int32_t COLVARS_BUG_ERROR       = (1<<3); // Inconsistent state indicating bug
  constexpr int32_t COLVARS_FILE_ERROR      = (1<<4);
  constexpr int32_t COLVARS_MEMORY_ERROR    = (1<<5);
  constexpr int32_t COLVARS_NO_SUCH_FRAME   = (1<<6); // Cannot load the requested frame
}


#endif
