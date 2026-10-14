#include "colvargrid_integrate.h"
#include "colvarproxy.h"

#include <cstdio>
#include <iomanip>
#include <iostream>

colvargrid_integrate::colvargrid_integrate(std::vector<colvar *> &colvars,
    std::shared_ptr<colvar_grid_gradient> gradients_in,
    bool weighted_in)
  : colvar_grid_scalar(colvars, gradients_in, true),
    b_smoothed(false),
    weighted(weighted_in),
    gradients(gradients_in)

{
  // parent class colvar_grid_scalar is constructed with add_extra_bin option set to true
  // hence PMF grid is wider than gradient grid if non-PBC

  if (nd > 1) {
    cvmodule->cite_feature("Poisson integration of 2D/3D free energy surfaces");
    divergence.resize(nt);
  }
  init_Poisson_computation();
}


colvargrid_integrate::colvargrid_integrate(std::shared_ptr<colvar_grid_gradient> gradients_in,
    bool weighted_in)
  : b_smoothed(false),
    weighted(weighted_in),
    gradients(gradients_in)
{
  nd = gradients->num_variables();
  nx = gradients->number_of_points_vec();
  widths = gradients->widths;
  periodic = gradients->periodic;


  // Expand grid by 1 bin in non-periodic dimensions
  for (size_t i = 0; i < nd; i++) {
    if (!periodic[i]) nx[i]++;
    // Shift the grid by half the bin width (values at edges instead of center of bins)
    lower_boundaries.push_back(gradients->lower_boundaries[i].real_value - 0.5 * widths[i]);
  }
  setup(nx);

  init_Poisson_computation();
}


int colvargrid_integrate::integrate(const int itmax, const cvm::real &tol, cvm::real &err,
                                    bool verbose)
{
  int iter = 0;
  if (nd == 1 && !weighted) {
    cvm::real sum = 0.0;
    cvm::real corr;
    if (periodic[0]) {
      corr = gradients->average(); // Enforce PBC by subtracting average gradient
    } else {
      corr = 0.0;
    }
    std::vector<int> ix;
    // Iterate over valid indices in gradient grid
    for (ix = new_index(); gradients->index_ok(ix); incr(ix)) {
      set_value(ix, sum);
      cvm::real val = gradients->value_output_smoothed(ix, b_smoothed);
      sum += (val - corr) * widths[0];
    }
    if (index_ok(ix)) {
      // This will happen if non-periodic: then PMF grid has one extra bin wrt gradient grid
      // If not, sum should be zero
      set_value(ix, sum);
    }
  } else {
    divergence.clear();
    divergence.resize(computation_nt);
    if (weighted || nd > 3) {
      if (weighted && !gradients->samples)
        cvmodule->error("Error: Trying to perform weighted Poisson integration without a samples grid.", COLVARS_BUG_ERROR);
      prepare_calculations();
      // extrapolate_data(); // Potential enhancement, needs testing
      set_weighted_div();
      std::vector<cvm::real> temp = divergence;
      std::vector<cvm::real> temp2(computation_grid->data);
      linewise_laplacian_weighted =
        precompute ? &colvargrid_integrate::linewise_laplacian_weighted_precomputed<true>
        : &colvargrid_integrate::linewise_laplacian_weighted_otf<true>;
      laplacian_weighted<true>(divergence, temp2);
      temp.clear();
      temp2.clear();
      for (size_t i = 0; i < computation_nt; i++) {
        divergence[i] += div_border_supplement[i];
      }
      div_border_supplement.clear();
    } else {
      set_unweighted_div();
    }

    cvm::real average = 0;
    for (size_t i = 0; i < computation_nt; i++) {
      average += divergence[i];
    }
    average = average / static_cast<cvm::real>(computation_nt);
    for (size_t i = 0; i < computation_nt; i++) {
      divergence[i] = divergence[i] - average;
    }

    if (weighted || nd > 3) {
      if (nd > 3) {
        cvmodule->log("WARNING: Integration of potential of dimension higher than 3 requires a lot "
                      "of memory.");
      }
      nr_linbcg_sym(weighted, divergence, computation_grid->data, tol, itmax, iter, err);
    } else
      nr_linbcg_sym(weighted, divergence, computation_grid->data, tol, itmax, iter, err);
    if (verbose)
      cvmodule->log("Integrated in " + cvm::to_str(iter) + " steps, error: " + cvm::to_str(err));
    if (need_to_extrapolate_solution) {
      extrapolate_potential();
    }
  }

  return iter;
}


void colvargrid_integrate::set_unweighted_div()
{
  if (nd == 1)
    return;
  for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
       computation_grid->incr(ix)) {
    update_div_local(ix);
  }
}

void colvargrid_integrate::set_weighted_div()
{
  for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
       computation_grid->incr(ix)) {
    update_weighted_div_local(ix);
  }
}

void colvargrid_integrate::update_div_neighbors(const std::vector<int> &ix0)
{
  std::vector<int> ix(ix0);
  int i, j, k;
  bool edge;
  // If not periodic, expanded grid ensures that upper neighbors of ix0 are valid grid points
  if (nd == 1) {
    return;
  } else if (nd == 2) {
    edge = computation_grid->wrap_detect_edge(ix);
    if (!edge) update_div_local(ix);
    ix[0]++;
    edge = computation_grid->wrap_detect_edge(ix);
    if (!edge) update_div_local(ix);
    ix[1]++;
    edge = computation_grid->wrap_detect_edge(ix);
    if (!edge) update_div_local(ix);
    ix[0]--;
    edge = computation_grid->wrap_detect_edge(ix);
    if (!edge) update_div_local(ix);
  } else if (nd == 3) {
    for (i = 0; i < 2; i++) {
      ix[1] = ix0[1];
      for (j = 0; j < 2; j++) {
        ix[2] = ix0[2];
        for (k = 0; k < 2; k++) {
          edge = computation_grid->wrap_detect_edge(ix);
          if (!edge) update_div_local(ix);
          ix[2]++;
        }
        ix[1]++;
      }
      ix[0]++;
    }
  }
}

void colvargrid_integrate::get_grad(std::vector<cvm::real> &g, std::vector<int> ix)
{
  gradients->vector_value(ix, g);
}


void colvargrid_integrate::update_div_local(const std::vector<int> &ix0)
{
  const size_t linear_index = computation_grid->address(ix0);
  std::vector<int> ix = ix0;
  if (nd == 2) {
    // gradients at grid points surrounding the current scalar grid point
    std::vector<cvm::real> g00(2, 0), g01(2, 0), g10(2, 0), g11(2, 0);
    get_grad(g00, ix);

    ix[0] = ix0[0] + 1;
    wrap_detect_edge(ix);
    get_grad(g10, ix);

    ix[0] = ix0[0] + 1;
    ix[1] = ix0[1] + 1;
    wrap_detect_edge(ix);
    get_grad(g11, ix);

    ix[0] = ix0[0];
    ix[1] = ix0[1] + 1;
    wrap_detect_edge(ix);
    get_grad(g01, ix);

    ix = ix0;

    divergence[linear_index] = ((g10[0] - g00[0] + g11[0] - g01[0]) / widths[0] +
                                (g01[1] - g00[1] + g11[1] - g10[1]) / widths[1]) *
                               0.5;
    // handle boundary conditions
    cvm::real supplement = 0;
    cvm::real supplement_coefficient = 1;
    if ((ix0[0] == 0 || ix0[0] == computation_nx[0] - 1) && !periodic[0]) {
      cvm::real coefficient = (ix[0] == 0) ? 1.0 : -1.0;
      supplement_coefficient = 0.5 * supplement_coefficient;
      supplement += coefficient * (g10[0] + g00[0] + g11[0] + g01[0]) * 0.25 * 2 / widths[0];
    }
    if ((ix0[1] == 0 || ix0[1] == computation_nx[1] - 1) && !periodic[1]) {
      cvm::real coefficient = (ix[1] == 0) ? 1.0 : -1.0;
      supplement_coefficient = 0.5 * supplement_coefficient;
      supplement += coefficient * (g01[1] + g00[1] + g11[1] + g10[1]) * 0.25 * 2 / widths[1];
    }
    divergence[linear_index] += supplement;
    // See Long Chen's paper: multiply boundary line to conserve symmetry
    divergence[linear_index] *= supplement_coefficient;

  } else if (nd == 3) {
    // 8 gradient vectors, each of size 3
    std::vector<cvm::real> g000(3, 0.0), g001(3, 0.0), g010(3, 0.0), g011(3, 0.0);
    std::vector<cvm::real> g100(3, 0.0), g101(3, 0.0), g110(3, 0.0), g111(3, 0.0);

    // Fill the 8 corners, shifting with +1
    get_grad(g000, ix);            // (0,0,0)
    ix[2] = ix0[2] + 1;
    wrap_detect_edge(ix);
    get_grad(g001, ix);            // (0,0,1)

    ix[1] = ix0[1] + 1;
    wrap_detect_edge(ix);
    get_grad(g011, ix);            // (0,1,1)

    ix[2] = ix0[2];
    wrap_detect_edge(ix);
    get_grad(g010, ix);            // (0,1,0)

    ix[0] = ix0[0] + 1;
    ix[1] = ix0[1];
    ix[2] = ix0[2];
    wrap_detect_edge(ix);
    get_grad(g100, ix);            // (1,0,0)

    ix[2] = ix0[2] + 1;
    wrap_detect_edge(ix);
    get_grad(g101, ix);            // (1,0,1)

    ix[1] = ix0[1] + 1;
    wrap_detect_edge(ix);
    get_grad(g111, ix);            // (1,1,1)

    ix[2] = ix0[2];
    wrap_detect_edge(ix);
    get_grad(g110, ix);            // (1,1,0)
    divergence[linear_index] =
      ((g100[0] - g000[0] + g101[0] - g001[0] + g110[0] - g010[0] + g111[0] - g011[0]) /
       widths[0] +
       (g010[1] - g000[1] + g011[1] - g001[1] + g110[1] - g100[1] + g111[1] - g101[1]) /
       widths[1] +
       (g001[2] - g000[2] + g011[2] - g010[2] + g101[2] - g100[2] + g111[2] - g110[2]) /
       widths[2]) *
      0.25;
    cvm::real supplement = 0;
    cvm::real supplement_coefficient = 1;

    ix = ix0;
    if ((ix0[0] == 0 || ix0[0] == computation_nx[0] - 1) && !periodic[0]) {
      cvm::real coefficient = (ix[0] == 0) ? 1.0 : -1.0;
      supplement_coefficient *= 0.5;
      supplement += coefficient * (g100[0] + g000[0] + g101[0] + g001[0] + g110[0] + g010[0] + g111[0] + g011[0]) *
        0.125 * 2 / widths[0];
    }
    if ((ix0[1] == 0 || ix0[1] == computation_nx[1] - 1) && !periodic[1]) {
      cvm::real coefficient = (ix[1] == 0) ? 1.0 : -1.0;
      supplement_coefficient *= 0.5;
      supplement +=
        coefficient *
        (g010[1] + g000[1] + g011[1] + g001[1] + g110[1] + g100[1] + g111[1] + g101[1]) * 2 * 0.125 / widths[1] ;
    }
    if ((ix0[2] == 0 || ix0[2] == computation_nx[2] - 1) && !periodic[2]) {
      cvm::real coefficient = (ix[2] == 0) ? 1.0 : -1.0;
      supplement_coefficient *= 0.5;
      supplement +=
        coefficient *
        (g001[2] + g000[2] + g011[2] + g010[2] + g101[2] + g100[2] + g111[2] + g110[2]) * 0.125 * 2 / widths[2];
    }
    divergence[linear_index] += supplement;
    divergence[linear_index] *= supplement_coefficient;
  }
}


void colvargrid_integrate::prepare_divergence_stencils()
{
  surrounding_points_relative_positions.clear();
  size_t n_combinations = cvm::integer_power(2, nd);
  for (size_t i = 0; i < n_combinations; i++) {
    surrounding_points_relative_positions.push_back(convert_base_two(i, nd));
  }
}

/// Updates the divergence at the point ix0
void colvargrid_integrate::update_weighted_div_local(const std::vector<int> &ix0)
{
  const size_t linear_index = computation_grid->address(ix0);
  cvm::real div_at_point = 0;
  for (const std::vector<int>& surrounding_point_relative_position :
       surrounding_points_relative_positions) {
    std::vector<int> surrounding_point_coordinates = ix0;
    std::vector<cvm::real> gradient_at_surrounding_point(nd, 0);
    for (size_t i = 0; i < nd; i++) {
      surrounding_point_coordinates[i] += surrounding_point_relative_position[i];
    }

    gradients->wrap_detect_edge(surrounding_point_coordinates);
    get_regularized_grad(gradient_at_surrounding_point, surrounding_point_coordinates);
    cvm::real weight = regularized_weights[gradients->address(surrounding_point_coordinates)];

    for (size_t i = 0; i < nd; i++) {
      div_at_point +=  (((surrounding_point_relative_position[i] & 1) << 1) - 1) *
                       gradient_at_surrounding_point[i] * weight / widths[i];
    }
  }
  divergence[linear_index] = div_at_point / cvm::integer_power(2, nd - 1);
}

/// Multiplication by sparse matrix representing Laplacian
/// NOTE: Laplacian must be symmetric for solving with CG
void colvargrid_integrate::laplacian(const std::vector<cvm::real> &A, std::vector<cvm::real> &LA)
{
  if (nd == 2) {
    // DIMENSION 2

    size_t li, li2;
    int i, j;
    cvm::real fact;
    const cvm::real ffx = 1.0 / (widths[0] * widths[0]);
    const cvm::real ffy = 1.0 / (widths[1] * widths[1]);
    const int h = computation_nx[1];
    const int w = computation_nx[0];
    // offsets for 4 reference points of the Laplacian stencil
    int xm = -h;
    int xp =  h;
    int ym = -1;
    int yp =  1;

    // NOTE on performance: this version is slightly sub-optimal because
    // it contains two double loops on the core of the array (for x and y terms)
    // The slightly faster version is in commit 0254cb5a2958cb2e135f268371c4b45fad34866b
    // yet it is much uglier, and probably horrible to extend to dimension 3
    // All terms in the matrix are assigned (=) during the x loops, then updated (+=)
    // with the y (and z) contributions


    // All x components except on x edges
    li = h; // Skip first column

    // Halve the term on y edges (if any) to preserve symmetry of the Laplacian matrix
    // (Long Chen, Finite Difference Methods, UCI, 2017)
    fact = periodic[1] ? 1.0 : 0.5;

    for (i = 1; i < w - 1; i++) {
      // Full range of j, but factor may change on y edges (j == 0 and j == h-1)
      LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
      li++;
      for (j = 1; j < h - 1; j++) {
        LA[li] = ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        li++;
      }
      LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
      li++;
    }
    // Edges along x (x components only)
    li = 0L; // Follows left edge
    li2 = h * static_cast<size_t>(w - 1); // Follows right edge
    if (periodic[0]) {
      xm =  h * (w - 1);
      xp =  h;
      fact = periodic[1] ? 1.0 : 0.5;
      LA[li]  = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
      LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
      li++;
      li2++;
      for (j = 1; j < h - 1; j++) {
        LA[li]  = ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        LA[li2] = ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
        li++;
        li2++;
      }
      LA[li]  = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
      LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
    } else {
      xm = -h;
      xp =  h;
      fact = periodic[1] ? 1.0 : 0.5; // Halve in corners in full PBC only
      // lower corner, "j == 0"
      LA[li]  = fact * ffx * (A[li + xp] - A[li]);
      LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
      li++;
      li2++;
      for (j = 1; j < h - 1; j++) {
        // x gradient (+ y term of laplacian, calculated below)
        LA[li]  = ffx * (A[li + xp] - A[li]);
        LA[li2] = ffx * (A[li2 + xm] - A[li2]);
        li++;
        li2++;
      }
      // upper corner, j == h-1
      LA[li]  = fact * ffx * (A[li + xp] - A[li]);
      LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
    }

    // Now adding all y components
    // All y components except on y edges
    li = 1; // Skip first element (in first row)

    fact = periodic[0] ? 1.0 : 0.5; // for i == 0
    for (i = 0; i < w; i++) {
      // Factor of 1/2 on x edges if non-periodic
      if (i == 1)
        fact = 1.0;
      if (i == w - 1)
        fact = periodic[0] ? 1.0 : 0.5;
      for (j = 1; j < h - 1; j++) {
        LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        li++;
      }
      li += 2; // skip the edges and move to next column
    }
    // Edges along y (y components only)
    li = 0L; // Follows bottom edge
    li2 = h - 1; // Follows top edge
    if (periodic[1]) {
      fact = periodic[0] ? 1.0 : 0.5;
      ym = h - 1;
      yp = 1;
      LA[li]  += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
      LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
      li  += h;
      li2 += h;
      for (i = 1; i < w - 1; i++) {
        LA[li]  += ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        LA[li2] += ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
        li  += h;
        li2 += h;
      }
      LA[li]  += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
      LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
    } else {
      ym = -1;
      yp = 1;
      fact = periodic[0] ? 1.0 : 0.5; // Halve in corners in full PBC only
      // Left corner
      LA[li]  += fact * ffy * (A[li + yp] - A[li]);
      LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
      li  += h;
      li2 += h;
      for (i = 1; i < w - 1; i++) {
        // y gradient (+ x term of laplacian, calculated above)
        LA[li]  += ffy * (A[li + yp] - A[li]);
        LA[li2] += ffy * (A[li2 + ym] - A[li2]);
        li  += h;
        li2 += h;
      }
      // Right corner
      LA[li]  += fact * ffy * (A[li + yp] - A[li]);
      LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
    }
  } else if (nd == 3) {
    // DIMENSION 3

    int i, j, k;
    size_t li, li2;
    cvm::real fact = 1.0;
    const cvm::real ffx = 1.0 / (widths[0] * widths[0]);
    const cvm::real ffy = 1.0 / (widths[1] * widths[1]);
    const cvm::real ffz = 1.0 / (widths[2] * widths[2]);
    const int h = computation_nx[2]; // height
    const int d = computation_nx[1]; // depth
    const int w = computation_nx[0]; // width
    // offsets for 6 reference points of the Laplacian stencil
    int xm = -d * h;
    int xp =  d * h;
    int ym = -h;
    int yp =  h;
    int zm = -1;
    int zp =  1;

    cvm::real factx = periodic[0] ? 1 : 0.5; // factor to be applied on x edges
    cvm::real facty = periodic[1] ? 1 : 0.5; // same for y
    cvm::real factz = periodic[2] ? 1 : 0.5; // same for z
    cvm::real ifactx = 1 / factx;
    cvm::real ifacty = 1 / facty;
    cvm::real ifactz = 1 / factz;

    // All x components except on x edges
    li = d * static_cast<size_t>(h); // Skip left slab
    fact = facty * factz;
    for (i = 1; i < w - 1; i++) {
      for (j = 0; j < d; j++) {
        // full range of y
        if (j == 1)
          fact *= ifacty;
        if (j == d - 1)
          fact *= facty;
        LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        li++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          // full range of z
          LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
          li++;
        }
        fact *= factz;
        LA[li] = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        li++;
      }
    }
    // Edges along x (x components only)
    li = 0L; // Follows left slab
    li2 = static_cast<size_t>(d) * h * (w - 1); // Follows right slab
    if (periodic[0]) {
      xm =  d * h * (w - 1);
      xp =  d * h;
      fact = facty * factz;
      for (j = 0; j < d; j++) {
        if (j == 1)
          fact *= ifacty;
        if (j == d - 1)
          fact *= facty;
        LA[li]  = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
        li++;
        li2++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          LA[li]  = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
          LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
          li++;
          li2++;
        }
        fact *= factz;
        LA[li]  = fact * ffx * (A[li + xm] + A[li + xp] - 2.0 * A[li]);
        LA[li2] = fact * ffx * (A[li2 - xp] + A[li2 - xm] - 2.0 * A[li2]);
        li++;
        li2++;
      }
    } else {
      xm = -d * h;
      xp =  d * h;
      fact = facty * factz;
      for (j = 0; j < d; j++) {
        if (j == 1)
          fact *= ifacty;
        if (j == d - 1)
          fact *= facty;
        LA[li]  = fact * ffx * (A[li + xp] - A[li]);
        LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
        li++;
        li2++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          // x gradient (+ y, z terms of laplacian, calculated below)
          LA[li]  = fact * ffx * (A[li + xp] - A[li]);
          LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
          li++;
          li2++;
        }
        fact *= factz;
        LA[li]  = fact * ffx * (A[li + xp] - A[li]);
        LA[li2] = fact * ffx * (A[li2 + xm] - A[li2]);
        li++;
        li2++;
      }
    }

    // Now adding all y components
    // All y components except on y edges
    li = h; // Skip first column (in front slab)
    fact = factx * factz;
    for (i = 0; i < w; i++) {
      // full range of x
      if (i == 1)
        fact *= ifactx;
      if (i == w - 1)
        fact *= factx;
      for (j = 1; j < d - 1; j++) {
        LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        li++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
          li++;
        }
        fact *= factz;
        LA[li] += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        li++;
      }
      li += 2 * h; // skip columns in front and back slabs
    }
    // Edges along y (y components only)
    li = 0L; // Follows front slab
    li2 = h * static_cast<size_t>(d - 1); // Follows back slab
    if (periodic[1]) {
      ym = h * (d - 1);
      yp = h;
      fact = factx * factz;
      for (i = 0; i < w; i++) {
        if (i == 1)
          fact *= ifactx;
        if (i == w - 1)
          fact *= factx;
        LA[li]  += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
        li++;
        li2++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          LA[li]  += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
          LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
          li++;
          li2++;
        }
        fact *= factz;
        LA[li]  += fact * ffy * (A[li + ym] + A[li + yp] - 2.0 * A[li]);
        LA[li2] += fact * ffy * (A[li2 - yp] + A[li2 - ym] - 2.0 * A[li2]);
        li++;
        li2++;
        li  += h * static_cast<size_t>(d - 1);
        li2 += h * static_cast<size_t>(d - 1);
      }
    } else {
      ym = -h;
      yp =  h;
      fact = factx * factz;
      for (i = 0; i < w; i++) {
        if (i == 1)
          fact *= ifactx;
        if (i == w - 1)
          fact *= factx;
        LA[li]  += fact * ffy * (A[li + yp] - A[li]);
        LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
        li++;
        li2++;
        fact *= ifactz;
        for (k = 1; k < h - 1; k++) {
          // y gradient (+ x, z terms of laplacian, calculated above and below)
          LA[li]  += fact * ffy * (A[li + yp] - A[li]);
          LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
          li++;
          li2++;
        }
        fact *= factz;
        LA[li]  += fact * ffy * (A[li + yp] - A[li]);
        LA[li2] += fact * ffy * (A[li2 + ym] - A[li2]);
        li++;
        li2++;
        li  += h * static_cast<size_t>(d - 1);
        li2 += h * static_cast<size_t>(d - 1);
      }
    }

    // Now adding all z components
    // All z components except on z edges
    li = 1; // Skip first element (in bottom slab)
    fact = factx * facty;
    for (i = 0; i < w; i++) {
      // full range of x
      if (i == 1)
        fact *= ifactx;
      if (i == w - 1)
        fact *= factx;
      for (k = 1; k < h - 1; k++) {
        LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
        li++;
      }
      fact *= ifacty;
      li += 2; // skip edge slabs
      for (j = 1; j < d - 1; j++) {
        // full range of y
        for (k = 1; k < h - 1; k++) {
          LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
          li++;
        }
        li += 2; // skip edge slabs
      }
      fact *= facty;
      for (k = 1; k < h - 1; k++) {
        LA[li] += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
        li++;
      }
      li += 2; // skip edge slabs
    }
    // Edges along z (z components onlz)
    li = 0; // Follows bottom slab
    li2 = h - 1; // Follows top slab
    if (periodic[2]) {
      zm = h - 1;
      zp = 1;
      fact = factx * facty;
      for (i = 0; i < w; i++) {
        if (i == 1)
          fact *= ifactx;
        if (i == w - 1)
          fact *= factx;
        LA[li]  += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
        LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
        li  += h;
        li2 += h;
        fact *= ifacty;
        for (j = 1; j < d - 1; j++) {
          LA[li]  += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
          LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
          li  += h;
          li2 += h;
        }
        fact *= facty;
        LA[li]  += fact * ffz * (A[li + zm] + A[li + zp] - 2.0 * A[li]);
        LA[li2] += fact * ffz * (A[li2 - zp] + A[li2 - zm] - 2.0 * A[li2]);
        li  += h;
        li2 += h;
      }
    } else {
      zm = -1;
      zp = 1;
      fact = factx * facty;
      for (i = 0; i < w; i++) {
        if (i == 1)
          fact *= ifactx;
        if (i == w - 1)
          fact *= factx;
        LA[li]  += fact * ffz * (A[li + zp] - A[li]);
        LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
        li  += h;
        li2 += h;
        fact *= ifacty;
        for (j = 1; j < d - 1; j++) {
          // z gradient (+ x, y terms of laplacian, calculated above)
          LA[li]  += fact * ffz * (A[li + zp] - A[li]);
          LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
          li  += h;
          li2 += h;
        }
        fact *= facty;
        LA[li]  += fact * ffz * (A[li + zp] - A[li]);
        LA[li2] += fact * ffz * (A[li2 + zm] - A[li2]);
        li  += h;
        li2 += h;
      }
    }
  }
}


template <bool initialize_div_supplement>
void colvargrid_integrate::linewise_laplacian_weighted_precomputed(const std::vector<cvm::real> &A,
    std::vector<cvm::real> &LA,
    size_t grid_address)
{
  std::vector<int> neighbor_coordinate(nd);
  std::vector<int> reference_point_coordinates(nd);
  std::vector<cvm::real> averaged_normal_vector(nd);
  cvm::real multiplicity = laplacian_stencil.size();
  std::vector<int> ix(nd);
  computation_grid->index(grid_address, ix);
  for (size_t i = 0; i < laplacian_stencil.size(); i++) {
    std::vector<int> neighbor_relative_position = laplacian_stencil[i];
    for (size_t dim = 0; dim < nd; dim++) {
      neighbor_coordinate[dim] = ix[dim] + neighbor_relative_position[dim];
    }
    bool is_ghost_point = computation_grid->wrap_detect_edge(neighbor_coordinate);
    cvm::real coefficient = laplacian_coefficients[grid_address * multiplicity + i];
    if (!is_ghost_point) {
      LA[grid_address] += coefficient * A[computation_grid->address(neighbor_coordinate)];
      // laplacian_matrix_test[computation_grid->address(ix) * computation_nt +
      // computation_grid->address(neighbor_coordinate)] += coefficient;
    } else {
      computation_grid->wrap_to_edge(neighbor_coordinate, reference_point_coordinates);
      LA[grid_address] += coefficient * A[computation_grid->address(reference_point_coordinates)];
      // laplacian_matrix_test[computation_grid->address(ix) * computation_nt +
      // computation_grid->address(reference_point_coordinates)] += coefficient;

      if (initialize_div_supplement) {
        cvm::real div_supplement_term = 0;
        averaged_normal_vector = compute_averaged_border_normal_gradient(neighbor_coordinate);
        for (size_t dim = 0; dim < nd; dim++) {
          div_supplement_term +=
            averaged_normal_vector[dim] * neighbor_relative_position[dim] * widths[dim];
        }
        div_border_supplement[grid_address] -= div_supplement_term * coefficient;
      }
    }
  }
}
/// Calculate the multiplication by the laplacian matrix at the line grid_address
/// I.e. Calculates the discrete laplacian of a function at a given point of the grid whose address is grid_address
/// Note: To store all the coefficients of the laplacian in a test matrix, uncomment the lines with
/// "laplacian_matrix_test"
template <bool initialize_div_supplement>
void colvargrid_integrate::linewise_laplacian_weighted_otf(const std::vector<cvm::real> &A,
    std::vector<cvm::real> &LA,
    size_t grid_address)
{
  std::vector<int> neighbor_coordinate(nd);
  std::vector<int> reference_point_coordinates(nd);
  std::vector<cvm::real> averaged_normal_vector(nd);
  std::vector<int> ix(nd);
  computation_grid->index(grid_address, ix);
  cvm::real diagonal_coeff = 0;
  for (size_t i = 0; i < laplacian_stencil.size() - 1; i++) {

    std::vector<int> neighbor_relative_position = laplacian_stencil[i];
    for (size_t dim = 0; dim < nd; dim++) {
      neighbor_coordinate[dim] = ix[dim] + neighbor_relative_position[dim];
    }

    cvm::real coefficient = 0;
    // sum all the weights that will be averaged
    for (std::vector<int> direction : weight_stencil[i]) {
      std::vector<int> weight_coordinate = ix;
      // compute weights' positions in the data grid
      for (size_t n = 0; n < nd && n < direction.size(); n++) {
        weight_coordinate[n] += direction[n];
      }
      gradients->wrap_detect_edge(weight_coordinate);
      coefficient += regularized_weights[gradients->address(weight_coordinate)];
    }
    coefficient *= weight_counts[i];
    diagonal_coeff += coefficient;
    bool is_ghost_point = computation_grid->wrap_detect_edge(neighbor_coordinate);
    if (!is_ghost_point) {
      LA[grid_address] += coefficient * A[computation_grid->address(neighbor_coordinate)];
      // laplacian_matrix_test[computation_grid->address(ix) * computation_nt +
      // computation_grid->address(neighbor_coordinate)] += coefficient;
    } else {
      computation_grid->wrap_to_edge(neighbor_coordinate, reference_point_coordinates);
      LA[grid_address] += coefficient * A[computation_grid->address(reference_point_coordinates)];
      // laplacian_matrix_test[computation_grid->address(ix) * computation_nt +
      // computation_grid->address(reference_point_coordinates)] += coefficient;

      if (initialize_div_supplement) {
        cvm::real div_supplement_term = 0;
        averaged_normal_vector = compute_averaged_border_normal_gradient(neighbor_coordinate);
        for (size_t dim = 0; dim < nd; dim++) {
          div_supplement_term +=
            averaged_normal_vector[dim] * neighbor_relative_position[dim] * widths[dim];
        }
        div_border_supplement[grid_address] -= div_supplement_term * coefficient;
      }
    }
  }
  LA[grid_address] += -diagonal_coeff * A[grid_address];
}

/// Multiplication by sparse matrix representing Laplacian
/// NOTE: Laplacian must be symmetric for solving with CG
template <bool initialize_div_supplement>
void colvargrid_integrate::laplacian_weighted(const std::vector<cvm::real> &A,
    std::vector<cvm::real> &LA)
{

  if (initialize_div_supplement) {
    div_border_supplement.resize(divergence.size());
  }
  for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
       computation_grid->incr(ix)) {
    LA[computation_grid->address(ix)] = 0;
    if (initialize_div_supplement)
      div_border_supplement[computation_grid->address(ix)] = 0;
  }
  // for testing, uncomment lines involving laplacian_matrix_test, also in
  // linewise_laplacian_weighted laplacian_matrix_test =
  // std::vector<cvm::real>(computation_nt*computation_nt, 0);
  if (m_num_threads == 1) {
    for (size_t grid_address = 0; grid_address < computation_grid->nt; grid_address++) {
      (this->*linewise_laplacian_weighted)(A, LA, grid_address);
    }
  } else {
#ifdef _OPENMP
    #pragma omp parallel for // schedule(static, 128)

    for (int grid_address_int = 0; grid_address_int < static_cast<int>(computation_grid->nt);
         grid_address_int++) {
      (this->*linewise_laplacian_weighted)(A, LA, static_cast<size_t>(grid_address_int));
    }
#else
    cvmodule->error("multiple threads required in weighted poisson integration, but this binary "
                    "is not linked with a supported threading library.\n");
#endif
  }
}

void colvargrid_integrate::prepare_laplacian_stencils()
{
  laplacian_stencil.resize(2 * nd + 1);
  weight_stencil.resize(2 * nd);
  weight_counts.resize(2 * nd);
  size_t const stencil_size = cvm::integer_power(2, nd - 1);

  for (size_t dim = 0; dim < nd; dim++) {
    std::vector<int> relative_position_left(nd, 0);
    relative_position_left[dim] = -1;
    std::vector<int> relative_position_right(nd, 0);
    relative_position_right[dim] = 1;
    laplacian_stencil[2 * dim] = relative_position_left;
    laplacian_stencil[2 * dim + 1] = relative_position_right;
    weight_counts[2 * dim] = 1. / (widths[dim] * widths[dim]) * (1. / stencil_size);
    weight_counts[2 * dim + 1] = 1. / (widths[dim] * widths[dim]) * (1. / stencil_size);
    weight_stencil[2 * dim].resize(stencil_size);
    weight_stencil[2 * dim + 1].resize(stencil_size);

    for (int weights_to_average_relative_pos = 0;
         weights_to_average_relative_pos < static_cast<int>(stencil_size); weights_to_average_relative_pos++) {
      std::vector<int> binary = convert_base_two(weights_to_average_relative_pos, nd - 1);
      weight_stencil[2 * dim][weights_to_average_relative_pos].resize(nd);
      weight_stencil[2 * dim + 1][weights_to_average_relative_pos].resize(nd);
      int offset = 0;
      for (size_t k = 0; k < nd; k++) {
        if (k != dim) {
          weight_stencil[2 * dim][weights_to_average_relative_pos][k] = binary[k + offset];
          weight_stencil[2 * dim + 1][weights_to_average_relative_pos][k] =
            binary[k + offset];
        } else {
          offset = -1;
          weight_stencil[2 * dim][weights_to_average_relative_pos][dim] = 0;
          weight_stencil[2 * dim + 1][weights_to_average_relative_pos][dim] = 1;
        }
      }
    }
  }
  laplacian_stencil[2 * nd] = std::vector<int>(nd, 0);
  // coefficient for the point where we compute the laplacian
  // is the sum of the stencil's point coefficient. we don't count it in the stencil
}


void colvargrid_integrate::prepare_calculations()
{
  if (weighted) {
    std::vector<int> max_position;
    std::vector<int> min_position;
    sorted_counts = {};

    for (std::vector<int> ix = gradients->new_index(); gradients->index_ok(ix);
         gradients->incr(ix)) {
      size_t count = gradients->samples->value(ix);
      if (count > 0) {
        insert_into_sorted_list<size_t>(sorted_counts, count);
      }
    }

    if (sorted_counts.size()) {
      lower_threshold_count = sorted_counts[static_cast<int>(sorted_counts.size() * lambda_min)];
      upper_threshold_count = sorted_counts[static_cast<int>(sorted_counts.size() * (1 - lambda_max))];
    } else {
      lower_threshold_count = 1;
      upper_threshold_count = 1;
    }
    sorted_counts.clear();
    regularized_weights.resize(gradients->nt);
    for (std::vector<int> ix = gradients->new_index(); gradients->index_ok(ix);
         gradients->incr(ix)) {
      regularized_weights[gradients->address(ix)] = get_regularized_weight(ix);
    }
    std::size_t required = (laplacian_stencil.size() + 2) * computation_nt * sizeof(cvm::real);
    double gigabytes = required / (1024.0 * 1024.0 * 1024.0);
    precompute = gigabytes < 2;
    std::string print_precompute = precompute ? "precomputed; precomputing will use "
                                   : "computed on the fly, as precomputing would use ";
    cvmodule->log("Laplacian computation will be " + print_precompute + cvm::to_str(gigabytes)
                  + " GB of memory.");
    if (precompute) {
      laplacian_coefficients.clear();
      laplacian_coefficients.resize(computation_grid->data.size() * laplacian_stencil.size(), 0);
      // laplacian_matrix_test = std::vector<cvm::real>(computation_nt*computation_nt, 0);

      std::vector<int> neighbor_coordinate(nd);
      std::vector<int> reference_point_coordinates(nd);
      std::vector<cvm::real> averaged_normal_vector(nd);
      cvm::real multiplicity = laplacian_stencil.size();
      // Precalculation of the laplacian coefficient
      cvm::real diagonal_coeff;
      size_t coeff_address;
      for (std::vector<int> ix = computation_grid->new_index(); computation_grid->index_ok(ix);
           computation_grid->incr(ix)) {
        coeff_address = computation_grid->address(ix) * multiplicity;
        diagonal_coeff = 0;
        for (size_t i = 0; i < laplacian_stencil.size() - 1; i++) {
          std::vector<int> neighbor_relative_position = laplacian_stencil[i];
          for (size_t dim = 0; dim < nd; dim++) {
            neighbor_coordinate[dim] = ix[dim] + neighbor_relative_position[dim];
          }
          cvm::real coefficient = 0;
          //  --> Calling the function makes the code simpler but takes longer to compute
          // coefficient = calculate_weight_sum(ix,weight_stencil[i]);
          for (std::vector<int> direction : weight_stencil[i]) {
            std::vector<int> weight_coordinate = ix;
            // Initialize with stencil_point instead of size
            for (size_t n = 0; n < nd && n < direction.size(); n++) {
              weight_coordinate[n] += direction[n];
            }
            gradients->wrap_detect_edge(weight_coordinate);
            coefficient += regularized_weights[gradients->address(weight_coordinate)];
          }
          coefficient *= weight_counts[i];
          laplacian_coefficients[coeff_address + i] += coefficient;
          diagonal_coeff += coefficient;
        }
        laplacian_coefficients[coeff_address + (2 * nd)] +=
          -diagonal_coeff;
      }
    }
  }
  is_calculations_prepared = true;
}


cvm::real colvargrid_integrate::get_regularized_weight(std::vector<int> &ix)
{
  cvm::real regularized_weight;
  size_t count = gradients->samples->value(ix);
  if (count < lower_threshold_count) {
    regularized_weight = lower_threshold_count;
  } else if (count > upper_threshold_count) {
    regularized_weight = upper_threshold_count;
  } else {
    regularized_weight = count;
  }
  return regularized_weight;
}


void colvargrid_integrate::get_regularized_grad(std::vector<cvm::real> &F, std::vector<int> &ix)
{
  size_t count = gradients->samples->value(ix);
  gradients->vector_value(ix, F);
  cvm::real multiplier = 1.;
  if (count < min_count_F) {
    multiplier = 0;
  } else if (count < max_count_F) {
    multiplier = static_cast<cvm::real>(count - min_count_F) / static_cast<cvm::real>(max_count_F - min_count_F);
  }
  for (size_t i = 0; i < nd; i++) {
    F[i] = multiplier * F[i];
  }
}


std::vector<cvm::real> colvargrid_integrate::compute_averaged_border_normal_gradient(
  std::vector<int> ghost_point_coordinates)
{
  std::vector<int> reference_point_coordinates(nd, 0);
  gradients->wrap_to_edge(ghost_point_coordinates, reference_point_coordinates);

  std::vector<int> directions_to_average_along;
  directions_to_average_along.reserve(nd);           // at most nd entries
  std::vector<bool> normal_directions(nd);

  for (size_t i = 0; i < nd; i++) {
    if ((0 <= ghost_point_coordinates[i] && ghost_point_coordinates[i] < computation_nx[i]) ||
        periodic[i]) {
      directions_to_average_along.push_back(i);
      normal_directions[i] = false;
    } else {
      normal_directions[i] = true;
    }
  }

  const int nb_averaged_gradients =
    static_cast<int>(cvm::integer_power(2, static_cast<int>(directions_to_average_along.size())));

  std::vector<std::vector<int>> gradients_to_average_relative_positions;
  gradients_to_average_relative_positions.reserve(nb_averaged_gradients); // known size

  // The i=0 case (all bits zero) produces an all-zero vector regardless of
  // directions_to_average_along
  const size_t n_dims_to_average = directions_to_average_along.size();
  for (int i = 0; i < nb_averaged_gradients; i++) {
    std::vector<int> pos(nd, 0);
    for (size_t bit = 0; bit < n_dims_to_average; bit++) {
      pos[directions_to_average_along[bit]] = (i >> bit) & 1;
    }
    gradients_to_average_relative_positions.push_back(std::move(pos));
  }
  const size_t nb_averaged = gradients_to_average_relative_positions.size();
  const cvm::real inv_nb_averaged = cvm::real(1) / static_cast<cvm::real>(nb_averaged);

  std::vector<cvm::real> averaged_border_normal_gradient(nd, 0);
  std::vector<int>       gradient_position(nd);
  std::vector<cvm::real> gradient(nd);

  for (size_t i = 0; i < nb_averaged; i++) {
    const auto& relative_pos = gradients_to_average_relative_positions[i];
    for (size_t j = 0; j < nd; j++) {
      gradient_position[j] = reference_point_coordinates[j] + relative_pos[j];
    }
    wrap_detect_edge(gradient_position);

    if (weighted) {
      get_regularized_grad(gradient, gradient_position);
    } else {
      get_grad(gradient, gradient_position);
    }

    for (size_t j = 0; j < nd; j++) {
      averaged_border_normal_gradient[j] += gradient[j];
    }
  }

  for (size_t j = 0; j < nd; j++) {
    averaged_border_normal_gradient[j] = normal_directions[j]
      ? averaged_border_normal_gradient[j] * inv_nb_averaged
      : cvm::real(0);
  }
  return averaged_border_normal_gradient;
}


std::vector<int> colvargrid_integrate::convert_base_two(int n, size_t length)
{
  unsigned int un = static_cast<unsigned int>(n);
  std::vector<int> result(length, 0);
  size_t i = length;
  while (un > 0) {
    i--;
    result[i] = un & 1;
    un >>= 1;
  }

  return result;
}


void colvargrid_integrate::nr_linbcg_sym(const bool weighted, const std::vector<cvm::real> &b,
    std::vector<cvm::real> &x, const cvm::real &tol,
    const int itmax, int &iter, cvm::real &err)
{
  cvm::real ak, akden, bk, bkden, bknum, bnrm;
  const cvm::real EPS = 1.0e-14;
  int j;
  std::vector<cvm::real> p(computation_nt), r(computation_nt), z(computation_nt);
  typedef void (colvargrid_integrate::*func_pointer)(const std::vector<double> &,
      std::vector<double> &);
  func_pointer atimes = weighted ? &colvargrid_integrate::laplacian_weighted<false>
                        : &colvargrid_integrate::laplacian;
  linewise_laplacian_weighted =
    precompute ? &colvargrid_integrate::linewise_laplacian_weighted_precomputed<false>
    : &colvargrid_integrate::linewise_laplacian_weighted_otf<false>;
  iter = 0;
  (this->*atimes)(x, r);
  for (j = 0; j < int(computation_nt); j++) {
    r[j] = b[j] - r[j];
  }
  bnrm = l2norm(b);
  if (bnrm < EPS) {
    return; // Target is zero, will break relative error calc
  }
  //   asolve(r,z); // precon
  bkden = 1.0;
  while (iter < itmax) {
    ++iter;
    for (bknum = 0.0, j = 0; j < int(computation_nt); j++) {
      bknum += r[j] * r[j]; // precon: z[j]*r[j]
    }
    if (iter == 1) {
      for (j = 0; j < int(computation_nt); j++) {
        p[j] = r[j];  // precon: p[j] = z[j]
      }
    } else {
      bk = bknum / bkden;
      for (j = 0; j < int(computation_nt); j++) {
        p[j] = bk * p[j] + r[j]; // precon:  bk*p[j] + z[j]
      }
    }
    bkden = bknum;
    (this->*atimes)(p, z);
    for (akden = 0.0, j = 0; j < int(computation_nt); j++) {
      akden += z[j] * p[j];
    }
    ak = bknum / akden;
    for (j = 0; j < int(computation_nt); j++) {
      x[j] += ak * p[j];
      r[j] -= ak * z[j];
    }
    //     asolve(r,z);  // precon
    err = l2norm(r) / bnrm;
    if (cvm::debug())
      std::cout << "iter=" << std::setw(4) << iter + 1 << std::setw(12) << err << std::endl;
    if (err <= tol)
      break;
  }
}


void colvargrid_integrate::extrapolate_potential()
{
  for (std::vector<int> ix = new_index(); index_ok(ix); incr(ix)) {
    std::vector<int> corresponding_index_in_small_grid = ix;
    for (size_t i = 0; i < nd; i++) {
      if (!periodic[i]) {
        corresponding_index_in_small_grid[i] = ix[i] - 1;
      }
    }
    std::vector<int> reference_index_in_small_grid(nd, 0);
    bool need_to_extrapolate = computation_grid->wrap_to_edge(corresponding_index_in_small_grid,
                               reference_index_in_small_grid);
    cvm::real potential_value =
      computation_grid->data[computation_grid->address(reference_index_in_small_grid)];
    std::vector<cvm::real> relative_position(nd, 0);
    if (need_to_extrapolate) {
      for (size_t i = 0; i < nd; i++) {
        relative_position[i] =
          std::min(corresponding_index_in_small_grid[i] - reference_index_in_small_grid[i], 1);
      }
      std::vector<cvm::real> averaged_normal_vector =
        compute_averaged_border_normal_gradient(corresponding_index_in_small_grid);
      for (size_t i = 0; i < nd; i++) {
        potential_value += averaged_normal_vector[i] * relative_position[i] * widths[i];
      }
    }
    data[address(ix)] = potential_value;
  }
};


template <typename T>
typename std::vector<T>::iterator
colvargrid_integrate::insert_into_sorted_list(std::vector<T> &sortedList, const T &value)
{
  // Find the first position where the element is not less than value
  auto it = std::lower_bound(sortedList.begin(), sortedList.end(), value);

  // Insert the value at the found position and return iterator to inserted element
  return sortedList.insert(it, value);
}

cvm::real colvargrid_integrate::l2norm(const std::vector<cvm::real> &x)
{
  size_t i;
  cvm::real sum = 0.0;
  for (i = 0; i < x.size(); i++)
    sum += x[i] * x[i];
  return sqrt(sum);
}


int colvargrid_integrate::init_Poisson_computation()
{
  if (nd == 1 && !weighted) {
    return COLVARS_OK;
  }

  cvm::main()->cite_feature("Poisson integration of 2D/3D free energy surfaces");
  computation_nx.resize(nd);
  computation_nt = 1;
  computation_nxc.resize(nd);
  for (size_t i = 0; i < nd; i++) {
    if (periodic[i]) {
      computation_nx[i] = nx[i];
    } else {
      computation_nx[i] = nx[i] - 2; // One point less than data size for non-periodic dimensions
    }
    computation_nt *= computation_nx[i];
    computation_nxc[i] = computation_nt;
  }
  divergence.resize(computation_nt);

  if (weighted) {
    div_border_supplement.resize(computation_nt);
    prepare_divergence_stencils();
    prepare_laplacian_stencils();
  }
  need_to_extrapolate_solution = false;
  for (size_t i = 0; i < nd; i++) {
    if (!periodic[i])
      need_to_extrapolate_solution = true;
  }
  if (!need_to_extrapolate_solution) {
    computation_grid = this;
  } else {
    computation_grid = new colvar_grid_scalar();
    computation_grid->periodic = periodic;
    computation_grid->setup(computation_nx);
  }
  cvmodule->log(cvm::to_str(nx[0]) + " " + "computation : " + cvm::to_str(computation_nx[0]));

#ifdef _OPENMP
  m_num_threads = cvmodule->proxy->smp_num_threads();
#else
  if (m_num_threads > 1) {
    return cvmodule->error("Multi-threading requested in weighted integrator, which is not supported "
                           "by this build.\n");
  }
#endif
  if (weighted) {
    cvmodule->log("Will perform weighted Poisson integrator using " + cvm::to_str(m_num_threads) + " threads.");
  }
  return COLVARS_OK;
}

// Potential enhancement, needs testing
// void colvargrid_integrate::extrapolate_data()
// {
//   // TODO: add the count thresholds' calculation ?
//   bool converged = false;
//   size_t step = 1;
//   std::vector<cvm::real> interpolated(gradients->data);
//   std::vector<size_t> known(gradients->nt, nt);
//   for (std::vector<int> ix = new_index(); gradients->index_ok(ix); incr(ix)) {
//     if (gradients->samples->value(ix) > lower_threshold_count) {
//       known[gradients->address(ix)] = step;
//       interpolated[gradients->address(ix)] = gradients->data[address(ix)];
//     }
//   }
//   std::vector<cvm::real> neighbor_grad(nd, 0);
//   cvm::real weight_total = 0.;
//   cvm::real weight_neighbor = 0.;
//   size_t known_neighbors = 0;
//   while (!converged) {
//     step++;
//     converged = true;
//     for (std::vector<int> ix = new_index(); gradients->index_ok(ix); incr(ix)) {
//       if (known[gradients->address(ix)] > step) {
//         known_neighbors = 0;
//         std::vector<cvm::real> interpolated_gradient(nd, 0);
//         weight_total = 0;
//         for (size_t neighbor = 0; neighbor < laplacian_stencil.size(); neighbor++) {
//           std::vector<int> neighbor_coordinates(ix);
//           for (size_t d = 0; d < nd; d++) {
//             neighbor_coordinates[d] += laplacian_stencil[neighbor][d];
//           }
//           if (known[gradients->address(neighbor_coordinates)] < step) {
//             known_neighbors++;
//             get_grad(neighbor_grad, neighbor_coordinates);
//             weight_neighbor = get_regularized_weight(neighbor_coordinates);
//             weight_total += weight_neighbor;
//             for (size_t d = 0; d < nd; d++) {
//               interpolated_gradient[d] += neighbor_grad[d] * weight_neighbor;
//             }
//           }
//         }
//         if (known_neighbors > 0) {
//           for (size_t d = 0; d < nd; d++) {
//             interpolated[address(ix) * nd + d] = interpolated_gradient[d] / weight_total;
//           }
//           known[gradients->address(ix)] = step;
//         } else
//           converged = false;
//       }
//     }
//   }
//   gradients->data = interpolated;
//   interpolated.clear();
//   known.clear();
// }
