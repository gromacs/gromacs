/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2017, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
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
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
/*! \internal \file
 * \brief
 * Declares nbnxn cuda cache and texture helper functions
 */
#ifndef GMX_MDLIB_NBNXN_CUDA_NBNXN_CUDA_H
#define GMX_MDLIB_NBNXN_CUDA_NBNXN_CUDA_H

#include "nbnxn_cuda_types.h"

//! Set up the cache configuration for the non-bonded kernels.
void nbnxn_cuda_set_cacheconfig(const gmx_device_info_t *devinfo);
/*! \brief Return the reference to the nbfp texture.
 *
 *  Note: it can return junk when c_disableCudaTextures==false, but we don't
 *  assert on that condition because the data_mgmt module ends up calling this
 *  function even if texture references are not used.
 */
const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_nbfp_texref();
/*! \brief Return the reference to the nbfp_comb texture.
 *
 *  Note: it can return junk when c_disableCudaTextures==false, but we don't
 *  assert on that condition because the data_mgmt module ends up calling this
 *  function even if texture references are not used.
 */
const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_nbfp_comb_texref();
/*! \brief Return the reference to the coulomb_tab texture.
 *
 *  Note: it can return junk when c_disableCudaTextures==false, but we don't
 *  assert on that condition because the data_mgmt module ends up calling this
 *  function even if texture references are not used.
 */
const struct texture<float, 1, cudaReadModeElementType> &nbnxn_cuda_get_coulomb_tab_texref();

#endif
