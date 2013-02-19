/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2012,2013, by the GROMACS development team, led by
 * David van der Spoel, Berk Hess, Erik Lindahl, and including many
 * others, as listed in the AUTHORS file in the top-level source
 * directory and at http://www.gromacs.org.
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

#ifndef GMX_HARDWARE_DETECT_H
#define GMX_HARDWARE_DETECT_H
#include "visibility.h"
#include "types/hw_info.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

GMX_LIBGMX_EXPORT
void gmx_detect_hardware(FILE *fplog, gmx_hw_info_t *hwinfo,
                         const t_commrec *cr,
                         gmx_bool bForceUseGPU, gmx_bool bTryUseGPU,
                         const char *gpu_id);

GMX_LIBGMX_EXPORT
void gmx_hardware_info_free(gmx_hw_info_t *hwinfo);

GMX_LIBGMX_EXPORT
void gmx_check_hw_runconf_consistency(FILE *fplog, gmx_hw_info_t *hwinfo,
                                      const t_commrec *cr, int ntmpi_requsted,
                                      gmx_bool bUseGPU);

#ifdef __cplusplus
}
#endif


#endif /* GMX_HARDWARE_DETECT_H */
