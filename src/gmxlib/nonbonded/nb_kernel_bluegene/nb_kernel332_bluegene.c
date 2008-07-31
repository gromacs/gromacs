/*
 * Copyright (c) 2006, International Business Machines (IBM) Inc.
 *
 * Author: Mathias Puetz
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 4. Neither the name of IBM nor the names of its contributors may be used
 * to endorse or promote products derived from this software without
 * specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY IBM AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
 * IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL IBM OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 * OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 * OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include<config.h>
#endif
#include <stdio.h>
#include <math.h>

#define COULOMB_NONE     0
#define COULOMB_CUTOFF   1
#define REACTION_FIELD   2
#define COULOMB_TAB      3
#define GENERALIZED_BORN 4

#define VDW_NONE      0
#define LENNARD_JONES 1
#define BUCKINGHAM    2
#define VDW_TAB       3

#define COULOMB   COULOMB_TAB
#define VDW       VDW_TAB
#define CONST_LJ  1

#include "interaction.h"

#undef NO_FORCE
#undef NB_KERNEL
#define NB_KERNEL nb_kernel332_bluegene

#include "nb_kernel_w3w3_bluegene.h"

#define NO_FORCE 1
#undef NB_KERNEL
#define NB_KERNEL nb_kernel332nf_bluegene

#include "nb_kernel_w3w3_bluegene.h"
