/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright 1991- The GROMACS Authors
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

#ifndef GMX_MDLIB_MDEBIN_BAR_H
#define GMX_MDLIB_MDEBIN_BAR_H

#include <cstdint>

#include <array>
#include <vector>

#include "gromacs/utility/real.h"

/* The functions & data structures here describe writing
   energy differences (or their histogram )for use with gmx bar */

class delta_h_history_t;
struct t_enxframe;
class energyhistory_t;
struct t_inputrec;

namespace gmx
{
template<typename>
class ArrayRef;
}

/* Data for one foreign lambda, or derivative. */
typedef struct
{
    std::vector<real>  dh;     /* the raw energy data. */
    std::vector<float> dhf;    /* raw difference data -- in floats, for storage. */
    unsigned int       ndh;    /* number of data points */
    unsigned int       ndhmax; /* the maximum number of points */

    int nhist;                           /* the number of histograms. There can either be
                                            0 (for no histograms)
                                            1 (for 'foreign lambda' histograms)
                                            2 (for derivative histograms: there's
                                               a 'forward' and 'backward' histogram
                                               containing the minimum and maximum
                                               values, respectively). */
    std::array<std::vector<int>, 2> bin; /* the histogram(s) */
    double                          dx;  /* the histogram spacing in kJ/mol. This is the
                                            same for the two histograms? */
    unsigned int           nbins;        /* the number of bins in the histograms*/
    std::array<int64_t, 2> x0;           /* the starting point in units of spacing
                                        of the histogram */
    std::array<unsigned int, 2> maxbin;  /* highest bin number with data */

    int type;                    /* the block type according to dhbtDH, etc. */
    int derivative;              /* The derivative direction (as an index in the lambda
                                    vector) if this delta_h contains derivatives */
    std::vector<double> lambda;  /* lambda vector (or NULL if not applicable) */
    int                 nlambda; /* length of the lambda vector */
    bool                written; /* whether this data has already been written out */

    std::array<int64_t, 5> subblock_meta_l; /* metadata for an mdebin subblock for
                                           I/O: for histogram counts, etc.*/
    std::vector<double> subblock_meta_d;    /* metadata subblock for I/O, used for
                                   communicating doubles (i.e. the lambda
                                   vector) */
    std::array<int, 4> subblock_meta_i;     /* metadata subblock for I/O, used for
                                   communicating ints (i.e. derivative indices,
                                   etc.) */
} t_mde_delta_h;

/* the type definition is in mdebin_bar.h */
struct t_mde_delta_h_coll
{
    t_mde_delta_h_coll(const t_inputrec& inputrec);

    std::vector<t_mde_delta_h> dh;  /* the delta h data */
    int                        ndh; /* the number of delta_h structures */

    int nlambda;     /* number of bar dU delta_h structures */
    int dh_du_index; /* the delta h data (index into dh) */

    int ndhdl;         /* number of bar dU delta_h structures */
    int dh_dhdl_index; /* the dhdl data (index into dh) */

    int dh_energy_index;   /* energy output block (index into dh) */
    int dh_pv_index;       /* pV output block (index into dh) */
    int dh_expanded_index; /* expanded ensemble output block (index into dh) */

    double start_time;     /* start time of the current dh collection */
    double delta_time;     /* time difference between samples */
    bool   start_time_set; /* whether the start time has been set */
    double start_lambda;   /* starting lambda for continuous motion of state*/
    double delta_lambda;   /* delta lambda, for continuous motion of state */
    double temperature;    /* the temperature of the samples*/

    std::vector<double> native_lambda_vec;     /* The lambda vector describing the current
                                      lambda state if it is set (NULL otherwise) */
    int              n_lambda_vec;             /* the size of the native lambda vector */
    std::vector<int> native_lambda_components; /* the native lambda (and by extension,
                                      foreign lambda) components in terms
                                      of efptFEP, efptMASS, etc. */
    int lambda_index;                          /* the lambda_fep_state */

    std::vector<double> subblock_d; /* for writing a metadata mdebin subblock for I/O */
    std::vector<int>    subblock_i; /* for writing a metadata mdebin subblock for I/O */
};

/* add a bunch of samples to the delta_h collection
    dhc = the collection
    dhdl = the hamiltonian derivatives
    U = the array with energies: from enerd->enerpart_lambda.
    time = the current simulation time.
    current_lambda = current lambda values : primarily useful for continuous processes
    fep_state = current fep_state
 */

/* add a bunch of samples - note fep_state is double to allow for better data storage */
void mde_delta_h_coll_add_dh(t_mde_delta_h_coll*   dhc,
                             double                fep_state,
                             double                energy,
                             double                pV,
                             gmx::ArrayRef<double> dhdl,
                             double*               foreign_dU,
                             double                time);

/* write the data associated with the du blocks collection as a collection
    of mdebin blocks.
    dhc = the collection
    fr = the enxio frame
    nblock = the current number of blocks */
void mde_delta_h_coll_handle_block(t_mde_delta_h_coll* dhc, t_enxframe* fr, int nblock);


/* reset the collection of delta_h buffers for a new round of
   data gathering */
void mde_delta_h_coll_reset(t_mde_delta_h_coll* dhc);


/* set the energyhistory variables to save state */
void mde_delta_h_coll_update_energyhistory(const t_mde_delta_h_coll* dhc, energyhistory_t* enerhist);

/* restore the variables from an energyhistory */
void mde_delta_h_coll_restore_energyhistory(t_mde_delta_h_coll* dhc, const delta_h_history_t* deltaH);

#endif /* GMX_MDLIB_MDEBIN_BAR_H */
