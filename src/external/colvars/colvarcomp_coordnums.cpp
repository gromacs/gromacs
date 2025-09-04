// -*- c++ -*-

// This file is part of the Collective Variables module (Colvars).
// The original version of Colvars and its updates are located at:
// https://github.com/Colvars/colvars
// Please update all Colvars source files before making any changes.
// If you wish to distribute your changes, please submit them to the
// Colvars repository at GitHub.

#include "colvarmodule.h"
#include "colvaratoms.h"
#include "colvarvalue.h"
#include "colvar.h"
#include "colvarcomp.h"



template<int flags>
cvm::real colvar::coordnum::switching_function(cvm::real const &r0,
                                               cvm::rvector const &r0_vec,
                                               int en,
                                               int ed,
                                               cvm::atom &A1,
                                               cvm::atom &A2,
                                               bool **pairlist_elem,
                                               cvm::real pairlist_tol)
{
  if ((flags & ef_use_pairlist) && !(flags & ef_rebuild_pairlist)) {
    bool const within = **pairlist_elem;
    (*pairlist_elem)++;
    if (!within) {
      return 0.0;
    }
  }

  cvm::rvector const r0sq_vec(r0_vec.x*r0_vec.x,
                              r0_vec.y*r0_vec.y,
                              r0_vec.z*r0_vec.z);

  cvm::rvector const diff = cvm::position_distance(A1.pos, A2.pos);

  cvm::rvector const scal_diff(diff.x/((flags & ef_anisotropic) ?
                                       r0_vec.x : r0),
                               diff.y/((flags & ef_anisotropic) ?
                                       r0_vec.y : r0),
                               diff.z/((flags & ef_anisotropic) ?
                                       r0_vec.z : r0));
  cvm::real const l2 = scal_diff.norm2();

  // Assume en and ed are even integers, and avoid sqrt in the following
  int const en2 = en/2;
  int const ed2 = ed/2;

  cvm::real const xn = cvm::integer_power(l2, en2);
  cvm::real const xd = cvm::integer_power(l2, ed2);
  //The subtraction and division stretches the function back to the range of [0,1] from [pairlist_tol,1]
  cvm::real const func_no_pairlist = (1.0-xn)/(1.0-xd);
  cvm::real func, inv_one_pairlist_tol;
  if (flags & ef_use_pairlist) {
    inv_one_pairlist_tol = 1 / (1.0-pairlist_tol);
    func = (func_no_pairlist - pairlist_tol) * inv_one_pairlist_tol;
  } else {
    func = func_no_pairlist;
  }

  if (flags & ef_rebuild_pairlist) {
    //Particles just outside of the cutoff also are considered if they come near.
    **pairlist_elem = (func > (-pairlist_tol * 0.5)) ? true : false;
    (*pairlist_elem)++;
  }
  //If the value is too small, we need to exclude it, rather than let it contribute to the sum or the gradients.
  if (func < 0)
    return 0;

  if (flags & ef_gradients) {
    //This is the old, completely correct expression for dFdl2:
    //cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2) -
    //                                        func*ed2*(xd/l2))*(-1.0);
    //This can become:
    //cvm::real const dFdl2 = (1.0/(1.0-xd))*(en2*(xn/l2)*(1.0-xn)/(1.0-xn) -
    //                                        func*ed2*(xd/l2))*(-1.0);
    //Recognizing that func = (1.0-xn)/(1.0-xd), we can group together the "func" and get a version of dFdl2 that is 0
    //when func=0, which lets us skip this gradient calculation when func=0.
    cvm::real dFdl2;
    if (flags & ef_use_pairlist) {
      dFdl2 = func_no_pairlist * inv_one_pairlist_tol * ((ed2*xd/((1.0-xd)*l2)) - (en2*xn/((1.0-xn)*l2)));
    } else {
      dFdl2 = func * ((ed2*xd/((1.0-xd)*l2)) - (en2*xn/((1.0-xn)*l2)));
    }
    cvm::rvector const dl2dx((2.0/((flags & ef_anisotropic) ? r0sq_vec.x :
                                   r0*r0)) * diff.x,
                             (2.0/((flags & ef_anisotropic) ? r0sq_vec.y :
                                   r0*r0)) * diff.y,
                             (2.0/((flags & ef_anisotropic) ? r0sq_vec.z :
                                   r0*r0)) * diff.z);
    A1.grad += (-1.0)*dFdl2*dl2dx;
    A2.grad +=        dFdl2*dl2dx;
  }

  return func;
}


colvar::coordnum::coordnum()
{
  set_function_type("coordNum");
  x.type(colvarvalue::type_scalar);
  colvarproxy *proxy = cvm::main()->proxy;
  r0 = proxy->angstrom_to_internal(4.0);
  r0_vec = cvm::rvector(proxy->angstrom_to_internal(4.0),
                        proxy->angstrom_to_internal(4.0),
                        proxy->angstrom_to_internal(4.0));
}


int colvar::coordnum::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  group1 = parse_group(conf, "group1");
  group2 = parse_group(conf, "group2");

  if (!group1 || !group2) {
    return error_code | COLVARS_INPUT_ERROR;
  }

  if (int atom_number = cvm::atom_group::overlap(*group1, *group2)) {
    error_code |= cvm::error(
        "Error: group1 and group2 share a common atom (number: " + cvm::to_str(atom_number) + ")\n",
        COLVARS_INPUT_ERROR);
  }

  if (group1->b_dummy) {
    error_code |=
        cvm::error("Error: only group2 is allowed to be a dummy atom\n", COLVARS_INPUT_ERROR);
  }

  bool const b_isotropic = get_keyval(conf, "cutoff", r0, r0);

  if (get_keyval(conf, "cutoff3", r0_vec, r0_vec)) {
    if (b_isotropic) {
      error_code |= cvm::error("Error: cannot specify \"cutoff\" and \"cutoff3\" "
                               "at the same time.\n",
                               COLVARS_INPUT_ERROR);
    }

    b_anisotropic = true;
    // remove meaningless negative signs
    if (r0_vec.x < 0.0) r0_vec.x *= -1.0;
    if (r0_vec.y < 0.0) r0_vec.y *= -1.0;
    if (r0_vec.z < 0.0) r0_vec.z *= -1.0;
  }

  get_keyval(conf, "expNumer", en, en);
  get_keyval(conf, "expDenom", ed, ed);

  if ( (en%2) || (ed%2) ) {
    error_code |= cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
                             COLVARS_INPUT_ERROR);
  }

  if ( (en <= 0) || (ed <= 0) ) {
    error_code |= cvm::error("Error: negative exponent(s) provided.\n",
                             COLVARS_INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

  get_keyval(conf, "group2CenterOnly", b_group2_center_only, group2->b_dummy);

  get_keyval(conf, "tolerance", tolerance, tolerance);
  if (tolerance > 0) {
    cvm::main()->cite_feature("coordNum pairlist");
    get_keyval(conf, "pairListFrequency", pairlist_freq, pairlist_freq);
    if ( ! (pairlist_freq > 0) ) {
      return cvm::error("Error: non-positive pairlistfrequency provided.\n",
                        COLVARS_INPUT_ERROR);
      // return and do not allocate the pairlists below
    }
    size_t n;
    if (b_group2_center_only) {
      n = group1->size();
    } else {
      n = group1->size() * group2->size();
    }
    pairlist = new bool[n];
    for (size_t i = 0; i < n; i++) pairlist[i] = true;
  }

  init_scalar_boundaries(0.0, b_group2_center_only ?
                         static_cast<cvm::real>(group1->size()) :
                         static_cast<cvm::real>(group1->size() *
                                                group2->size()));

  return error_code;
}


colvar::coordnum::~coordnum()
{
  if (pairlist) {
    delete [] pairlist;
  }
}


template<int flags> void colvar::coordnum::main_loop(bool **pairlist_elem)
{
  if (b_group2_center_only) {
    cvm::atom group2_com_atom;
    group2_com_atom.pos = group2->center_of_mass();
    for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++) {
      x.real_value += switching_function<flags>(r0, r0_vec, en, ed,
                                                *ai1, group2_com_atom,
                                                pairlist_elem,
                                                tolerance);
    }
    if (b_group2_center_only) {
      group2->set_weighted_gradient(group2_com_atom.grad);
    }
  } else {
    for (cvm::atom_iter ai1 = group1->begin(); ai1 != group1->end(); ai1++) {
      for (cvm::atom_iter ai2 = group2->begin(); ai2 != group2->end(); ai2++) {
        x.real_value += switching_function<flags>(r0, r0_vec, en, ed,
                                                  *ai1, *ai2,
                                                  pairlist_elem,
                                                  tolerance);
      }
    }
  }
}


template<int compute_flags> int colvar::coordnum::compute_coordnum()
{
  bool const use_pairlist = (pairlist != NULL);
  bool const rebuild_pairlist = (pairlist != NULL) &&
    (cvm::step_relative() % pairlist_freq == 0);

  bool *pairlist_elem = use_pairlist ? pairlist : NULL;

  if (b_anisotropic) {

    if (use_pairlist) {
      if (rebuild_pairlist) {
        int const flags = compute_flags | ef_anisotropic | ef_use_pairlist |
          ef_rebuild_pairlist;
        main_loop<flags>(&pairlist_elem);
      } else {
        int const flags = compute_flags | ef_anisotropic | ef_use_pairlist;
        main_loop<flags>(&pairlist_elem);
      }

    } else {

      int const flags = compute_flags | ef_anisotropic;
      main_loop<flags>(NULL);
    }

  } else {

    if (use_pairlist) {

      if (rebuild_pairlist) {
        int const flags = compute_flags | ef_use_pairlist | ef_rebuild_pairlist;
        main_loop<flags>(&pairlist_elem);
      } else {
        int const flags = compute_flags | ef_use_pairlist;
        main_loop<flags>(&pairlist_elem);
      }

    } else {

      int const flags = compute_flags;
      main_loop<flags>(NULL);
    }
  }

  return COLVARS_OK;
}


void colvar::coordnum::calc_value()
{
  x.real_value = 0.0;
  if (is_enabled(f_cvc_gradient)) {
    compute_coordnum<ef_gradients>();
  } else {
    compute_coordnum<ef_null>();
  }
}


void colvar::coordnum::calc_gradients()
{
  // Gradients are computed by calc_value() if f_cvc_gradients is enabled
}



// h_bond member functions

colvar::h_bond::h_bond()
{
  colvarproxy *proxy = cvm::main()->proxy;
  r0 = proxy->angstrom_to_internal(3.3);
  set_function_type("hBond");
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);
}


int colvar::h_bond::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  if (cvm::debug())
    cvm::log("Initializing h_bond object.\n");

  set_function_type("hBond");
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);

  int a_num = -1, d_num = -1;
  get_keyval(conf, "acceptor", a_num, a_num);
  get_keyval(conf, "donor",    d_num, a_num);

  if ( (a_num == -1) || (d_num == -1) ) {
    error_code |= cvm::error("Error: either acceptor or donor undefined.\n", COLVARS_INPUT_ERROR);
  }

  cvm::atom acceptor = cvm::atom(a_num);
  cvm::atom donor    = cvm::atom(d_num);
  register_atom_group(new cvm::atom_group);
  atom_groups[0]->add_atom(acceptor);
  atom_groups[0]->add_atom(donor);

  get_keyval(conf, "cutoff",   r0, r0);
  get_keyval(conf, "expNumer", en, en);
  get_keyval(conf, "expDenom", ed, ed);

  if ((en % 2) || (ed % 2)) {
    error_code |= cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
                             COLVARS_INPUT_ERROR);
  }

  if ((en <= 0) || (ed <= 0)) {
    error_code |= cvm::error("Error: negative exponent(s) provided.\n", COLVARS_INPUT_ERROR);
  }

  if (cvm::debug())
    cvm::log("Done initializing h_bond object.\n");

  return error_code;
}


colvar::h_bond::h_bond(cvm::atom const &acceptor,
                       cvm::atom const &donor,
                       cvm::real r0_i, int en_i, int ed_i)
  : h_bond()
{
  r0 = r0_i;
  en = en_i;
  ed = ed_i;
  register_atom_group(new cvm::atom_group);
  atom_groups[0]->add_atom(acceptor);
  atom_groups[0]->add_atom(donor);
}


void colvar::h_bond::calc_value()
{
  int const flags = coordnum::ef_null;
  cvm::rvector const r0_vec(0.0); // TODO enable the flag?
  x.real_value =
    coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                        (*atom_groups[0])[0],
                                        (*atom_groups[0])[1],
                                        NULL, 0.0);
}


void colvar::h_bond::calc_gradients()
{
  int const flags = coordnum::ef_gradients;
  cvm::rvector const r0_vec(0.0); // TODO enable the flag?
  coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                      (*atom_groups[0])[0],
                                      (*atom_groups[0])[1],
                                      NULL, 0.0);
}



colvar::selfcoordnum::selfcoordnum()
{
  set_function_type("selfCoordNum");
  x.type(colvarvalue::type_scalar);
  r0 = cvm::main()->proxy->angstrom_to_internal(4.0);
}


int colvar::selfcoordnum::init(std::string const &conf)
{
  int error_code = cvc::init(conf);

  group1 = parse_group(conf, "group1");

  if (!group1 || group1->size() == 0) {
    return error_code | COLVARS_INPUT_ERROR;
  }

  get_keyval(conf, "cutoff", r0, r0);
  get_keyval(conf, "expNumer", en, en);
  get_keyval(conf, "expDenom", ed, ed);


  if ((en % 2) || (ed % 2)) {
    error_code |= cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
                             COLVARS_INPUT_ERROR);
  }

  if ((en <= 0) || (ed <= 0)) {
    error_code |= cvm::error("Error: negative exponent(s) provided.\n", COLVARS_INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

  get_keyval(conf, "tolerance", tolerance, tolerance);
  if (tolerance > 0) {
    get_keyval(conf, "pairListFrequency", pairlist_freq, pairlist_freq);
    if ( ! (pairlist_freq > 0) ) {
      error_code |= cvm::error("Error: non-positive pairlistfrequency provided.\n",
                               COLVARS_INPUT_ERROR);
    }
    size_t n = (group1->size()-1) * (group1->size()-1);
    pairlist = new bool[n];
    for (size_t i = 0; i < n; i++) pairlist[i] = true;
  }

  init_scalar_boundaries(0.0, static_cast<cvm::real>((group1->size()-1) *
                                                     (group1->size()-1)));

  return error_code;
}


colvar::selfcoordnum::~selfcoordnum()
{
  if (pairlist) {
    delete [] pairlist;
  }
}


template<int compute_flags> int colvar::selfcoordnum::compute_selfcoordnum()
{
  cvm::rvector const r0_vec(0.0); // TODO enable the flag?

  bool const use_pairlist = (pairlist != NULL);
  bool const rebuild_pairlist = (pairlist != NULL) &&
    (cvm::step_relative() % pairlist_freq == 0);

  bool *pairlist_elem = use_pairlist ? pairlist : NULL;
  size_t i = 0, j = 0;
  size_t const n = group1->size();

  // Always isotropic (TODO: enable the ellipsoid?)

  if (use_pairlist) {

    if (rebuild_pairlist) {
      int const flags = compute_flags | coordnum::ef_use_pairlist |
        coordnum::ef_rebuild_pairlist;
      for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
          x.real_value +=
            coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                                (*group1)[i],
                                                (*group1)[j],
                                                &pairlist_elem,
                                                tolerance);
        }
      }
    } else {
      int const flags = compute_flags | coordnum::ef_use_pairlist;
      for (i = 0; i < n - 1; i++) {
        for (j = i + 1; j < n; j++) {
          x.real_value +=
            coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                                (*group1)[i],
                                                (*group1)[j],
                                                &pairlist_elem,
                                                tolerance);
        }
      }
    }

  } else { // if (use_pairlist) {

    int const flags = compute_flags | coordnum::ef_null;
    for (i = 0; i < n - 1; i++) {
      for (j = i + 1; j < n; j++) {
        x.real_value +=
          coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                              (*group1)[i],
                                              (*group1)[j],
                                              &pairlist_elem,
                                              tolerance);
      }
    }
  }

  return COLVARS_OK;
}


void colvar::selfcoordnum::calc_value()
{
  x.real_value = 0.0;
  if (is_enabled(f_cvc_gradient)) {
    compute_selfcoordnum<coordnum::ef_gradients>();
  } else {
    compute_selfcoordnum<coordnum::ef_null>();
  }
}


void colvar::selfcoordnum::calc_gradients()
{
  // Gradients are computed by calc_value() if f_cvc_gradients is enabled
}



colvar::groupcoordnum::groupcoordnum()
{
  set_function_type("groupCoord");
  x.type(colvarvalue::type_scalar);
  init_scalar_boundaries(0.0, 1.0);
  colvarproxy *proxy = cvm::main()->proxy;
  r0 = proxy->angstrom_to_internal(4.0);
  r0_vec = cvm::rvector(proxy->angstrom_to_internal(4.0),
                        proxy->angstrom_to_internal(4.0),
                        proxy->angstrom_to_internal(4.0));
}


int colvar::groupcoordnum::init(std::string const &conf)
{
  int error_code = distance::init(conf);

  // group1 and group2 are already initialized by distance()
  if (group1->b_dummy || group2->b_dummy) {
    return cvm::error("Error: neither group can be a dummy atom\n", COLVARS_INPUT_ERROR);
  }

  bool const b_scale = get_keyval(conf, "cutoff", r0, r0);

  if (get_keyval(conf, "cutoff3", r0_vec, r0_vec)) {
    if (b_scale) {
      error_code |=
          cvm::error("Error: cannot specify \"cutoff\" and \"cutoff3\" at the same time.\n",
                     COLVARS_INPUT_ERROR);
    }
    b_anisotropic = true;
    // remove meaningless negative signs
    if (r0_vec.x < 0.0) r0_vec.x *= -1.0;
    if (r0_vec.y < 0.0) r0_vec.y *= -1.0;
    if (r0_vec.z < 0.0) r0_vec.z *= -1.0;
  }

  get_keyval(conf, "expNumer", en, en);
  get_keyval(conf, "expDenom", ed, ed);

  if ((en % 2) || (ed % 2)) {
    error_code |= cvm::error("Error: odd exponent(s) provided, can only use even ones.\n",
                             COLVARS_INPUT_ERROR);
  }

  if ((en <= 0) || (ed <= 0)) {
    error_code |= cvm::error("Error: negative exponent(s) provided.\n", COLVARS_INPUT_ERROR);
  }

  if (!is_enabled(f_cvc_pbc_minimum_image)) {
    cvm::log("Warning: only minimum-image distances are used by this variable.\n");
  }

  return error_code;
}


void colvar::groupcoordnum::calc_value()
{
  // create fake atoms to hold the com coordinates
  cvm::atom group1_com_atom;
  cvm::atom group2_com_atom;
  group1_com_atom.pos = group1->center_of_mass();
  group2_com_atom.pos = group2->center_of_mass();
  if (b_anisotropic) {
    int const flags = coordnum::ef_anisotropic;
    x.real_value = coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                                       group1_com_atom,
                                                       group2_com_atom,
                                                       NULL, 0.0);
  } else {
    int const flags = coordnum::ef_null;
    x.real_value = coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                                       group1_com_atom,
                                                       group2_com_atom,
                                                       NULL, 0.0);
  }
}


void colvar::groupcoordnum::calc_gradients()
{
  cvm::atom group1_com_atom;
  cvm::atom group2_com_atom;
  group1_com_atom.pos = group1->center_of_mass();
  group2_com_atom.pos = group2->center_of_mass();

  if (b_anisotropic) {
    int const flags = coordnum::ef_gradients | coordnum::ef_anisotropic;
    coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                        group1_com_atom,
                                        group2_com_atom,
                                        NULL, 0.0);
  } else {
    int const flags = coordnum::ef_gradients;
    coordnum::switching_function<flags>(r0, r0_vec, en, ed,
                                        group1_com_atom,
                                        group2_com_atom,
                                        NULL, 0.0);
  }

  group1->set_weighted_gradient(group1_com_atom.grad);
  group2->set_weighted_gradient(group2_com_atom.grad);
}
