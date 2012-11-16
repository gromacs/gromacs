/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 *
 * 
 * This file is part of GROMACS.
 * Copyright (c) 2012-  
 *
 * Written by the Gromacs development team under coordination of
 * David van der Spoel, Berk Hess, and Erik Lindahl.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org
 * 
 * And Hey:
 * GROup of MAchos and Cynical Suckers
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
