/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/*
 * Copyright (c) 1997-1999 Massachusetts Institute of Technology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *
 */

/* rconfig.c -- this file contains all the real-complex codelets
 * the system knows about */

#include <fftw-int.h>
#include <rfftw.h>

#define NOTW_CODELET(x) \
	 &fftw_real2hc_##x##_desc
#define NOTWI_CODELET(x) \
	 &fftw_hc2real_##x##_desc

#define TWIDDLE_CODELET(x) \
	 &fftw_hc2hc_forward_##x##_desc
#define TWIDDLEI_CODELET(x) \
	 &fftw_hc2hc_backward_##x##_desc

/* automatically-generated list of codelets */

extern fftw_codelet_desc fftw_real2hc_1_desc;
extern fftw_codelet_desc fftw_hc2real_1_desc;
extern fftw_codelet_desc fftw_real2hc_2_desc;
extern fftw_codelet_desc fftw_hc2real_2_desc;
extern fftw_codelet_desc fftw_real2hc_3_desc;
extern fftw_codelet_desc fftw_hc2real_3_desc;
extern fftw_codelet_desc fftw_real2hc_4_desc;
extern fftw_codelet_desc fftw_hc2real_4_desc;
extern fftw_codelet_desc fftw_real2hc_5_desc;
extern fftw_codelet_desc fftw_hc2real_5_desc;
extern fftw_codelet_desc fftw_real2hc_6_desc;
extern fftw_codelet_desc fftw_hc2real_6_desc;
extern fftw_codelet_desc fftw_real2hc_7_desc;
extern fftw_codelet_desc fftw_hc2real_7_desc;
extern fftw_codelet_desc fftw_real2hc_8_desc;
extern fftw_codelet_desc fftw_hc2real_8_desc;
extern fftw_codelet_desc fftw_real2hc_9_desc;
extern fftw_codelet_desc fftw_hc2real_9_desc;
extern fftw_codelet_desc fftw_real2hc_10_desc;
extern fftw_codelet_desc fftw_hc2real_10_desc;
extern fftw_codelet_desc fftw_real2hc_11_desc;
extern fftw_codelet_desc fftw_hc2real_11_desc;
extern fftw_codelet_desc fftw_real2hc_12_desc;
extern fftw_codelet_desc fftw_hc2real_12_desc;
extern fftw_codelet_desc fftw_real2hc_13_desc;
extern fftw_codelet_desc fftw_hc2real_13_desc;
extern fftw_codelet_desc fftw_real2hc_14_desc;
extern fftw_codelet_desc fftw_hc2real_14_desc;
extern fftw_codelet_desc fftw_real2hc_15_desc;
extern fftw_codelet_desc fftw_hc2real_15_desc;
extern fftw_codelet_desc fftw_real2hc_16_desc;
extern fftw_codelet_desc fftw_hc2real_16_desc;
extern fftw_codelet_desc fftw_real2hc_32_desc;
extern fftw_codelet_desc fftw_hc2real_32_desc;
extern fftw_codelet_desc fftw_real2hc_64_desc;
extern fftw_codelet_desc fftw_hc2real_64_desc;
extern fftw_codelet_desc fftw_real2hc_128_desc;
extern fftw_codelet_desc fftw_hc2real_128_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_2_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_2_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_3_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_3_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_4_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_4_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_5_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_5_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_6_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_6_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_7_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_7_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_8_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_8_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_9_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_9_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_10_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_10_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_16_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_16_desc;
extern fftw_codelet_desc fftw_hc2hc_forward_32_desc;
extern fftw_codelet_desc fftw_hc2hc_backward_32_desc;

fftw_codelet_desc *rfftw_config[] =
{
     NOTW_CODELET(1),
     NOTWI_CODELET(1),
     NOTW_CODELET(2),
     NOTWI_CODELET(2),
     NOTW_CODELET(3),
     NOTWI_CODELET(3),
     NOTW_CODELET(4),
     NOTWI_CODELET(4),
     NOTW_CODELET(5),
     NOTWI_CODELET(5),
     NOTW_CODELET(6),
     NOTWI_CODELET(6),
     NOTW_CODELET(7),
     NOTWI_CODELET(7),
     NOTW_CODELET(8),
     NOTWI_CODELET(8),
     NOTW_CODELET(9),
     NOTWI_CODELET(9),
     NOTW_CODELET(10),
     NOTWI_CODELET(10),
     NOTW_CODELET(11),
     NOTWI_CODELET(11),
     NOTW_CODELET(12),
     NOTWI_CODELET(12),
     NOTW_CODELET(13),
     NOTWI_CODELET(13),
     NOTW_CODELET(14),
     NOTWI_CODELET(14),
     NOTW_CODELET(15),
     NOTWI_CODELET(15),
     NOTW_CODELET(16),
     NOTWI_CODELET(16),
     NOTW_CODELET(32),
     NOTWI_CODELET(32),
     NOTW_CODELET(64),
     NOTWI_CODELET(64),
     NOTW_CODELET(128),
     NOTWI_CODELET(128),
     TWIDDLE_CODELET(2),
     TWIDDLEI_CODELET(2),
     TWIDDLE_CODELET(3),
     TWIDDLEI_CODELET(3),
     TWIDDLE_CODELET(4),
     TWIDDLEI_CODELET(4),
     TWIDDLE_CODELET(5),
     TWIDDLEI_CODELET(5),
     TWIDDLE_CODELET(6),
     TWIDDLEI_CODELET(6),
     TWIDDLE_CODELET(7),
     TWIDDLEI_CODELET(7),
     TWIDDLE_CODELET(8),
     TWIDDLEI_CODELET(8),
     TWIDDLE_CODELET(9),
     TWIDDLEI_CODELET(9),
     TWIDDLE_CODELET(10),
     TWIDDLEI_CODELET(10),
     TWIDDLE_CODELET(16),
     TWIDDLEI_CODELET(16),
     TWIDDLE_CODELET(32),
     TWIDDLEI_CODELET(32),
     (fftw_codelet_desc *) 0
};
