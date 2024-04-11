/// -*- c++ -*-

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cerrno>


#include "gromacs/math/units.h"
#include "gromacs/mdtypes/inputrec.h"
#include "gromacs/mdtypes/forceoutput.h"
#include "gromacs/mdtypes/enerdata.h"
#include "gromacs/pbcutil/pbc.h"
#include "gromacs/utility/fatalerror.h"
#include "gromacs/utility/futil.h"
#include "gromacs/mdtypes/commrec.h"
#include "gromacs/domdec/domdec_struct.h"
#include "gromacs/gmxlib/network.h"
#include "gromacs/domdec/ga2la.h"
#include "gromacs/mdtypes/colvarshistory.h"
#include "gromacs/topology/ifunc.h"
#include "gromacs/topology/mtop_util.h"
#include "gromacs/mdlib/broadcaststructs.h"


#include "colvarproxy_gromacs.h"


//************************************************************
// colvarproxy_gromacs
colvarproxy_gromacs::colvarproxy_gromacs() : colvarproxy() {}

// Colvars Initialization
void colvarproxy_gromacs::init(t_inputrec *ir, int64_t step, const gmx_mtop_t &mtop,
                               ObservablesHistory* oh,
                               const std::string &prefix,
                               gmx::ArrayRef<const std::string> filenames_config,
                               const std::string &filename_restart,
                               const t_commrec *cr,
                               const rvec x[],
                               rvec **xa_old_whole_colvars_state_p,
                               int *n_colvars_atoms_state_p) {


  // Initialize colvars.
  first_timestep = true;
  restart_frequency_s = 0;

  // User-scripted forces are not available in GROMACS
  have_scripts = false;

  angstrom_value_ = 0.1;

  boltzmann_ = gmx::c_boltz;

  // Get the thermostat temperature.
  // NOTE: Considers only the first temperature coupling group!
  set_target_temperature(ir->opts.ref_t[0]);

  // GROMACS random number generation.
  // Seed with the mdp parameter ld_seed, the Langevin dynamics seed.
  rng.seed(ir->ld_seed);

  /// Handle input filenames and prefix/suffix for colvars files.
  ///
  /// filename_config is the colvars configuration file collected from "-colvars" option.
  /// The output prefix will be the prefix of Gromacs log filename.
  /// or "output" otherwise.
  ///
  /// For restart, 'filename_restart' is the colvars input file for restart,
  /// set by the "-cv_restart" option. It will be NULL otherwise.
  ///

  if(!prefix.empty())
  {
    output_prefix_str = prefix;
  }
  else {
    output_prefix_str = "output";
  }

  restart_output_prefix_str = prefix + ".restart";

  colvars_restart = false;

  if(!filename_restart.empty())
  {
    colvars_restart = true;
    input_prefix_str = filename_restart;
    // Don't strip the input_prefix_str because colvarmodule.cpp doesn't know that restart file from GROMACS needs the .dat extension.
  }

  // Retrieve masses and charges from input file
  updated_masses_ = updated_charges_ = true;

  // Get GROMACS timestep (picosecond to femtosecond)
  set_integration_timestep(ir->delta_t * 1000.0);
  // Retrieve the topology of all atoms
  gmx_atoms = gmx_mtop_global_atoms(mtop);

  // Read configuration file and set up the proxy only on the master node.
  if (MASTER(cr))
  {

    // initiate module: this object will be the communication proxy
    // colvarmodule pointer is only defined on the Master due to the static pointer to colvarproxy.
    colvars = new colvarmodule(this);

    version_int = get_version_from_string(COLVARPROXY_VERSION);

    colvars->cite_feature("GROMACS engine");
    colvars->cite_feature("Colvars-GROMACS interface");

    if (cvm::debug()) {
      log("Initializing the colvars proxy object.\n");
    }

    cvm::log("Using GROMACS interface, version "+
      cvm::to_str(COLVARPROXY_VERSION)+".\n");

    auto i = filenames_config.begin();
    for(; i != filenames_config.end(); ++i) {
        add_config("configfile", i->c_str());
    }

    colvarproxy::parse_module_config();
    colvars->update_engine_parameters();
    colvars->setup_input();

    // Citation Reporter
    cvm::log(std::string("\n")+colvars->feature_report(0)+std::string("\n"));

    colvars->setup_output();

    if (step != 0) {
      cvm::log("Initializing step number to "+cvm::to_str(step)+".\n");
    }

    colvars->it = colvars->it_restart = step;


  } // end master


  // MPI initialisation

  // Initialise attributs for the MPI communication
  if(MASTER(cr)) {
    // Retrieve the number of colvar atoms
    n_colvars_atoms = atoms_ids.size();
    // Copy their global indices
    ind = atoms_ids.data(); // This has to be updated if the vector is reallocated
  }


  if(PAR(cr)) {
    // Let the other nodes know the number of colvar atoms.
    block_bc(cr->mpi_comm_mygroup, n_colvars_atoms);

    // Initialise atoms_new_colvar_forces on non-master nodes
    if(!MASTER(cr)) {
      atoms_new_colvar_forces.reserve(n_colvars_atoms);
    }
  }

  snew(x_colvars_unwrapped,         n_colvars_atoms);
  snew(xa_ind,       n_colvars_atoms);
  snew(xa_shifts,    n_colvars_atoms);
  snew(xa_eshifts,   n_colvars_atoms);
  snew(xa_old_whole, n_colvars_atoms);
  snew(f_colvars,    n_colvars_atoms);

  // Prepare data

  // Manage restart with .cpt
  if (MASTER(cr))
  {
      /* colvarsHistory is the struct holding the data saved in the cpt

       If we dont start with from a .cpt, prepare the colvarsHistory struct for proper .cpt writing,
       If we did start from .cpt, we copy over the last whole structures from .cpt,
       In any case, for subsequent checkpoint writing, we set the pointers (xa_old_whole_p) in
       the xa_old_whole arrays, which contain the correct PBC representation of
       colvars atoms at the last time step.
      */

      if (oh->colvarsHistory == nullptr)
      {
          oh->colvarsHistory = std::make_unique<colvarshistory_t>(colvarshistory_t{});
      }
      colvarshistory_t *colvarshist = oh->colvarsHistory.get();


      snew(colvarshist->xa_old_whole_p, n_colvars_atoms);

      /* We always need the last whole positions such that
      * in the next time step we can make the colvars atoms whole again in PBC */
      if (colvarshist->bFromCpt)
      {
          for (int i = 0; i < n_colvars_atoms; i++)
            {
                copy_rvec(colvarshist->xa_old_whole[i], xa_old_whole[i]);
            }
      }
      else
      {
          colvarshist->n_atoms = n_colvars_atoms;
          for (int i = 0; i < n_colvars_atoms; i++)
          {
              int ii = ind[i];
              copy_rvec(x[ii], xa_old_whole[i]);
          }
      }

      /* For subsequent checkpoint writing, set the pointers (xa_old_whole_p) to the xa_old_whole
      * arrays that get updated at every NS step */
      colvarshist->xa_old_whole_p = xa_old_whole;
      //Initialize number of colvars atoms from the global state
      *n_colvars_atoms_state_p = n_colvars_atoms;
      // Point the shifts array from the  global state to the local shifts array
      *xa_old_whole_colvars_state_p = xa_old_whole;
  }


  // Communicate initial coordinates and global indices to all processes
  if (PAR(cr))
  {
    nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, xa_old_whole);
    snew_bc(MASTER(cr), ind, n_colvars_atoms);
    nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, ind);
  }

  // Serial Run
  if (!PAR(cr))
  {
    nat_loc = n_colvars_atoms;
    nalloc_loc = n_colvars_atoms;
    ind_loc = ind;

    // xa_ind[i] needs to be set to i for serial runs
    for (int i = 0; i < n_colvars_atoms; i++)
    {
        xa_ind[i] = i;
    }
  }

  if (MASTER(cr) && cvm::debug()) {
    cvm::log ("atoms_ids = "+cvm::to_str (atoms_ids)+"\n");
    cvm::log ("atoms_refcount = "+cvm::to_str (atoms_refcount)+"\n");
    cvm::log ("positions = "+cvm::to_str (atoms_positions)+"\n");
    cvm::log ("atoms_new_colvar_forces = "+cvm::to_str (atoms_new_colvar_forces)+"\n");
    cvm::log (cvm::line_marker);
    log("done initializing the colvars proxy object.\n");
  }


} // End colvars initialization.


colvarproxy_gromacs::~colvarproxy_gromacs()
{}

void colvarproxy_gromacs::finish(const t_commrec *cr)
{
  if(MASTER(cr)) {
    colvars->write_restart_file(output_prefix_str+".colvars.state");
    colvars->write_output_files();
  }
}

cvm::real colvarproxy_gromacs::rand_gaussian()
{
  return  normal_distribution(rng);
}

size_t colvarproxy_gromacs::restart_frequency()
{
  return restart_frequency_s;
}

// **************** PERIODIC BOUNDARY CONDITIONS ****************
//  Get the PBC-aware distance vector between two positions
cvm::rvector colvarproxy_gromacs::position_distance (cvm::atom_pos const &pos1,
                                                     cvm::atom_pos const &pos2) const
{
  rvec r1, r2, dr;
  r1[0] = pos1.x;
  r1[1] = pos1.y;
  r1[2] = pos1.z;
  r2[0] = pos2.x;
  r2[1] = pos2.y;
  r2[2] = pos2.z;

  pbc_dx(&gmx_pbc, r2, r1, dr);
  return cvm::atom_pos( dr[0], dr[1], dr[2] );
}


void colvarproxy_gromacs::log (std::string const &message)
{
  std::istringstream is(message);
  std::string line;
  while (std::getline(is, line))
    // Gromacs prints messages on the stderr FILE
    fprintf(stderr, "colvars: %s\n", line.c_str());
}

void colvarproxy_gromacs::error (std::string const &message)
{
  // In GROMACS, all errors are fatal.
  fatal_error (message);
}

void colvarproxy_gromacs::fatal_error (std::string const &message)
{
  log(message);
  if (!cvm::debug())
    log("If this error message is unclear, "
	"try recompiling with -DCOLVARS_DEBUG.\n");
  gmx_fatal(FARGS,"Error in collective variables module.\n");
}

void colvarproxy_gromacs::exit (std::string const gmx_unused &message)
{
  gmx_fatal(FARGS,"SUCCESS: %s\n", message.c_str());
}

int colvarproxy_gromacs::load_atoms (char const gmx_unused *filename, std::vector<cvm::atom> gmx_unused &atoms,
                                     std::string const gmx_unused &pdb_field, double const gmx_unused pdb_field_value)
{
  cvm::error("Selecting collective variable atoms "
		   "from a PDB file is currently not supported.\n");
  return COLVARS_NOT_IMPLEMENTED;
}

int colvarproxy_gromacs::load_coords (char const gmx_unused *filename, std::vector<cvm::atom_pos> gmx_unused &pos,
                                      const std::vector<int> gmx_unused &indices, std::string const gmx_unused &pdb_field_str,
                                      double const gmx_unused pdb_field_value)
{
  cvm::error("Loading atoms coordinates from a PDB or GRO file is currently not supported."
             "Please use an XYZ file.\n");
  return COLVARS_NOT_IMPLEMENTED;
}

int colvarproxy_gromacs::set_unit_system(std::string const &units_in, bool /*colvars_defined*/)
{
  if (units_in != "gromacs") {
    cvm::error("Specified unit system \"" + units_in + "\" is unsupported in Gromacs. Supported units are \"gromacs\" (nm, kJ/mol).\n");
    return COLVARS_ERROR;
  }
  return COLVARS_OK;
}

int colvarproxy_gromacs::backup_file (char const *filename)
{
  // Incremental gromacs backup system will be use only for those file
  if (std::string(filename).rfind(std::string(".colvars.traj")) != std::string::npos) {

    // GROMACS function
    make_backup(filename);

  // Otherwise, just keep one backup.
  } else {

    //Handle filename of the backup file
    const char *extension = ".old";
    char *backup = new char[strlen(filename)+strlen(extension)+1];
    strcpy(backup, filename);
    strcat(backup, extension);

    gmx_file_copy(filename, backup, FALSE);

    delete [] backup;

  }
  return COLVARS_OK;
}


void colvarproxy_gromacs::update_data(const t_commrec *cr, int64_t const step, t_pbc const &pbc, const matrix box, bool bNS)
{

  if (MASTER(cr)) {

    if(cvm::debug()) {
      cvm::log(cvm::line_marker);
      cvm::log("colvarproxy_gromacs, step no. "+cvm::to_str(colvars->it)+"\n"+
              "Updating internal data.\n");
    }

    // step update on master only due to the call of colvars pointer.
    if (first_timestep) {
      first_timestep = false;
    } else {
      if ( step - previous_gmx_step == 1 )
        colvars->it++;
      // Other cases?
    }
  } // end master

  gmx_pbc = pbc;
  gmx_box = box;
  gmx_bNS = bNS;

  previous_gmx_step = step;

  // Prepare data for MPI communication
  if(PAR(cr) && bNS) {
    dd_make_local_group_indices(cr->dd->ga2la, n_colvars_atoms, ind, &nat_loc, &ind_loc, &nalloc_loc, xa_ind);
  }
}


void colvarproxy_gromacs::calculateForces(
                    const gmx::ForceProviderInput &forceProviderInput,
                    gmx::ForceProviderOutput      *forceProviderOutput)
{

  const t_commrec *cr           = &(forceProviderInput.cr_);
  // Local atom coords
  const gmx::ArrayRef<const gmx::RVec> x  = forceProviderInput.x_;
  // Local atom coords (coerced into into old gmx type)
  const rvec *x_pointer          = &(x.data()->as_vec());


  // Eventually there needs to be an interface to update local data upon neighbor search
  // We could check if by chance all atoms are in one node, and skip communication
  communicate_group_positions(cr, x_colvars_unwrapped, xa_shifts, xa_eshifts,
                              gmx_bNS, x_pointer, n_colvars_atoms, nat_loc,
                              ind_loc, xa_ind, xa_old_whole, gmx_box);

  // Communicate_group_positions takes care of removing shifts (unwrapping)
  // in single node jobs, communicate_group_positions() is efficient and adds no overhead

  if (MASTER(cr))
  {
    // On non-master nodes, jump directly to applying the forces

    // Zero the forces on the atoms, so that they can be accumulated by the colvars.
    for (size_t i = 0; i < atoms_new_colvar_forces.size(); i++) {
      atoms_new_colvar_forces[i].x = atoms_new_colvar_forces[i].y = atoms_new_colvar_forces[i].z = 0.0;
    }

    // Get the atom positions from the Gromacs array.
    for (size_t i = 0; i < atoms_ids.size(); i++) {
      atoms_positions[i] = cvm::rvector(x_colvars_unwrapped[i][0], x_colvars_unwrapped[i][1], x_colvars_unwrapped[i][2]);
    }

    bias_energy = 0.0;
    // Call the collective variable module to fill atoms_new_colvar_forces
    if (colvars->calc() != COLVARS_OK) {
      cvm::error("Error calling colvars->calc()\n");
    }

    // Copy the forces to C array for broadcasting
    for (int i = 0; i < n_colvars_atoms; i++)
    {
      f_colvars[i][0] = atoms_new_colvar_forces[i].x;
      f_colvars[i][1] = atoms_new_colvar_forces[i].y;
      f_colvars[i][2] = atoms_new_colvar_forces[i].z;
    }

    forceProviderOutput->enerd_.term[F_COM_PULL] += bias_energy;
  } // master node

  //Broadcast the forces to all the nodes
  if (PAR(cr))
  {
    nblock_bc(cr->mpi_comm_mygroup, n_colvars_atoms, f_colvars);
  }

  const gmx::ArrayRef<gmx::RVec> &f_out = forceProviderOutput->forceWithVirial_.force_;
  matrix local_colvars_virial = { { 0 } };
  const bool computeVirial = forceProviderOutput->forceWithVirial_.computeVirial_;

  // Pass the applied forces back to GROMACS
  for (int i = 0; i < n_colvars_atoms; i++)
  {
    int i_global = ind[i];

    // check if this is a local atom and find out locndx
    if (PAR(cr)) {
      const int *locndx = cr->dd->ga2la->findHome(i_global);
      if (locndx) {
        f_out[*locndx] += f_colvars[i];
        if (computeVirial) {
          add_virial_term(local_colvars_virial, f_colvars[i], x_colvars_unwrapped[i]);
        }
      }
      // Do nothing if atom is not local
    } else { // Non MPI-parallel
      f_out[i_global] += f_colvars[i];
      if (computeVirial) {
        add_virial_term(local_colvars_virial, f_colvars[i], x_colvars_unwrapped[i]);
      }
    }
  }

  if (computeVirial) {
    forceProviderOutput->forceWithVirial_.addVirialContribution(local_colvars_virial);
  }
  return;
}


void colvarproxy_gromacs::add_virial_term(matrix vir, rvec const f, gmx::RVec const x)
{
  for (int j = 0; j < DIM; j++) {
    for (int m = 0; m < DIM; m++) {
      vir[j][m] -= 0.5 * f[j] * x[m];
    }
  }
}


// Pass restraint energy value for current timestep to MD engine
void colvarproxy_gromacs::add_energy (cvm::real energy)
{
  bias_energy += energy;
}

// **************** ATOMS ****************

int colvarproxy_gromacs::check_atom_id(int atom_number)
{
  // GROMACS uses zero-based arrays.
  int const aid = (atom_number-1);

  if (cvm::debug())
    log("Adding atom "+cvm::to_str(atom_number)+
        " for collective variables calculation.\n");

  if ( (aid < 0) || (aid >= gmx_atoms.nr) ) {
    cvm::error("Error: invalid atom number specified, "+
               cvm::to_str(atom_number)+"\n", COLVARS_INPUT_ERROR);
    return COLVARS_INPUT_ERROR;
  }

  return aid;
}


int colvarproxy_gromacs::init_atom(int atom_number)
{
  // GROMACS uses zero-based arrays.
  int aid = atom_number-1;

  for (size_t i = 0; i < atoms_ids.size(); i++) {
    if (atoms_ids[i] == aid) {
      // this atom id was already recorded
      atoms_refcount[i] += 1;
      return i;
    }
  }

  aid = check_atom_id(atom_number);

  if(aid < 0) {
    return COLVARS_INPUT_ERROR;
  }

  int const index = add_atom_slot(aid);
  update_atom_properties(index);
  return index;
}

void colvarproxy_gromacs::update_atom_properties(int index)
{

  // update mass
  double const mass = gmx_atoms.atom[atoms_ids[index]].m;
  if (mass <= 0.001) {
    this->log("Warning: near-zero mass for atom "+
              cvm::to_str(atoms_ids[index]+1)+
              "; expect unstable dynamics if you apply forces to it.\n");
  }
  atoms_masses[index] = mass;
  // update charge
  atoms_charges[index] = gmx_atoms.atom[atoms_ids[index]].q;
}
