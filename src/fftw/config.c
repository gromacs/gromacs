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

/* config.c -- this file contains all the codelets the system knows about */

/* $Id$ */

#include <fftw-int.h>

/* the signature is the same as the size, for now */
#define NOTW_CODELET(x) \
	 &fftw_no_twiddle_##x##_desc
#define NOTWI_CODELET(x) \
	 &fftwi_no_twiddle_##x##_desc

#define TWIDDLE_CODELET(x) \
	 &fftw_twiddle_##x##_desc

#define TWIDDLEI_CODELET(x) \
	 &fftwi_twiddle_##x##_desc

/* automatically-generated list of codelets */

extern fftw_codelet_desc fftw_no_twiddle_1_desc;
extern fftw_codelet_desc fftwi_no_twiddle_1_desc;
extern fftw_codelet_desc fftw_no_twiddle_2_desc;
extern fftw_codelet_desc fftwi_no_twiddle_2_desc;
extern fftw_codelet_desc fftw_no_twiddle_3_desc;
extern fftw_codelet_desc fftwi_no_twiddle_3_desc;
extern fftw_codelet_desc fftw_no_twiddle_4_desc;
extern fftw_codelet_desc fftwi_no_twiddle_4_desc;
extern fftw_codelet_desc fftw_no_twiddle_5_desc;
extern fftw_codelet_desc fftwi_no_twiddle_5_desc;
extern fftw_codelet_desc fftw_no_twiddle_6_desc;
extern fftw_codelet_desc fftwi_no_twiddle_6_desc;
extern fftw_codelet_desc fftw_no_twiddle_7_desc;
extern fftw_codelet_desc fftwi_no_twiddle_7_desc;
extern fftw_codelet_desc fftw_no_twiddle_8_desc;
extern fftw_codelet_desc fftwi_no_twiddle_8_desc;
extern fftw_codelet_desc fftw_no_twiddle_9_desc;
extern fftw_codelet_desc fftwi_no_twiddle_9_desc;
extern fftw_codelet_desc fftw_no_twiddle_10_desc;
extern fftw_codelet_desc fftwi_no_twiddle_10_desc;
extern fftw_codelet_desc fftw_no_twiddle_11_desc;
extern fftw_codelet_desc fftwi_no_twiddle_11_desc;
extern fftw_codelet_desc fftw_no_twiddle_12_desc;
extern fftw_codelet_desc fftwi_no_twiddle_12_desc;
extern fftw_codelet_desc fftw_no_twiddle_13_desc;
extern fftw_codelet_desc fftwi_no_twiddle_13_desc;
extern fftw_codelet_desc fftw_no_twiddle_14_desc;
extern fftw_codelet_desc fftwi_no_twiddle_14_desc;
extern fftw_codelet_desc fftw_no_twiddle_15_desc;
extern fftw_codelet_desc fftwi_no_twiddle_15_desc;
extern fftw_codelet_desc fftw_no_twiddle_16_desc;
extern fftw_codelet_desc fftwi_no_twiddle_16_desc;
extern fftw_codelet_desc fftw_no_twiddle_32_desc;
extern fftw_codelet_desc fftwi_no_twiddle_32_desc;
extern fftw_codelet_desc fftw_no_twiddle_64_desc;
extern fftw_codelet_desc fftwi_no_twiddle_64_desc;
extern fftw_codelet_desc fftw_twiddle_2_desc;
extern fftw_codelet_desc fftwi_twiddle_2_desc;
extern fftw_codelet_desc fftw_twiddle_3_desc;
extern fftw_codelet_desc fftwi_twiddle_3_desc;
extern fftw_codelet_desc fftw_twiddle_4_desc;
extern fftw_codelet_desc fftwi_twiddle_4_desc;
extern fftw_codelet_desc fftw_twiddle_5_desc;
extern fftw_codelet_desc fftwi_twiddle_5_desc;
extern fftw_codelet_desc fftw_twiddle_6_desc;
extern fftw_codelet_desc fftwi_twiddle_6_desc;
extern fftw_codelet_desc fftw_twiddle_7_desc;
extern fftw_codelet_desc fftwi_twiddle_7_desc;
extern fftw_codelet_desc fftw_twiddle_8_desc;
extern fftw_codelet_desc fftwi_twiddle_8_desc;
extern fftw_codelet_desc fftw_twiddle_9_desc;
extern fftw_codelet_desc fftwi_twiddle_9_desc;
extern fftw_codelet_desc fftw_twiddle_10_desc;
extern fftw_codelet_desc fftwi_twiddle_10_desc;
extern fftw_codelet_desc fftw_twiddle_16_desc;
extern fftw_codelet_desc fftwi_twiddle_16_desc;
extern fftw_codelet_desc fftw_twiddle_32_desc;
extern fftw_codelet_desc fftwi_twiddle_32_desc;
extern fftw_codelet_desc fftw_twiddle_64_desc;
extern fftw_codelet_desc fftwi_twiddle_64_desc;

fftw_codelet_desc *fftw_config[] =
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
     TWIDDLE_CODELET(64),
     TWIDDLEI_CODELET(64),
     (fftw_codelet_desc *) 0
};
