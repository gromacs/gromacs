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

#include "types/hwinfo.h"

#ifdef __cplusplus
extern "C" {
#endif
#if 0
} /* fixes auto-indentation problems */
#endif

void gmx_hw_detect(FILE *fplog, gmx_hwinfo_t *hwinfo,
                   const t_commrec *cr,
                   int cutoff_scheme, const char *nbpu_opt, const char *gpu_id);

void gmx_hw_info_free(gmx_hwinfo_t *hwinfo);

void gmx_check_hw_runconf_consistency(FILE *fplog, gmx_hwinfo_t *hwinfo,
                                      const t_commrec *cr, int ntmpi_requsted,
                                      const char *nbpu_opt);

#ifdef __cplusplus
}
#endif


#endif /* GMX_HARDWARE_DETECT_H */
