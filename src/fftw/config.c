/*
 * Copyright (c) 1997 Massachusetts Institute of Technology
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to use, copy, modify, and distribute the Software without
 * restriction, provided the Software, including any modified copies made
 * under this license, is not distributed for a fee, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE MASSACHUSETTS INSTITUTE OF TECHNOLOGY BE LIABLE
 * FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
 * CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * 
 * Except as contained in this notice, the name of the Massachusetts
 * Institute of Technology shall not be used in advertising or otherwise
 * to promote the sale, use or other dealings in this Software without
 * prior written authorization from the Massachusetts Institute of
 * Technology.
 *  
 */

/* config.c -- this file contains all the codelets the system knows about */

/* $Id$ */

#ifdef FFTW_USING_CILK
#include <cilk.h>
#include <cilk-compat.h>
#endif

#include <fftw.h>

/* the signature is the same as the size, for now */
#define NOTW_CODELET(x)  { x, x, fftw_no_twiddle_##x }
#define NOTWI_CODELET(x)  { x, x, fftwi_no_twiddle_##x }

extern notw_codelet fftw_no_twiddle_1;
extern notw_codelet fftw_no_twiddle_2;
extern notw_codelet fftw_no_twiddle_3;
extern notw_codelet fftw_no_twiddle_4;
extern notw_codelet fftw_no_twiddle_5;
extern notw_codelet fftw_no_twiddle_6;
extern notw_codelet fftw_no_twiddle_7;
extern notw_codelet fftw_no_twiddle_8;
extern notw_codelet fftw_no_twiddle_9;
extern notw_codelet fftw_no_twiddle_10;
extern notw_codelet fftw_no_twiddle_11;
extern notw_codelet fftw_no_twiddle_12;
extern notw_codelet fftw_no_twiddle_13;
extern notw_codelet fftw_no_twiddle_14;
extern notw_codelet fftw_no_twiddle_15;
extern notw_codelet fftw_no_twiddle_16;
extern notw_codelet fftw_no_twiddle_32;
extern notw_codelet fftw_no_twiddle_64;

extern notw_codelet fftwi_no_twiddle_1;
extern notw_codelet fftwi_no_twiddle_2;
extern notw_codelet fftwi_no_twiddle_3;
extern notw_codelet fftwi_no_twiddle_4;
extern notw_codelet fftwi_no_twiddle_5;
extern notw_codelet fftwi_no_twiddle_6;
extern notw_codelet fftwi_no_twiddle_7;
extern notw_codelet fftwi_no_twiddle_8;
extern notw_codelet fftwi_no_twiddle_9;
extern notw_codelet fftwi_no_twiddle_10;
extern notw_codelet fftwi_no_twiddle_11;
extern notw_codelet fftwi_no_twiddle_12;
extern notw_codelet fftwi_no_twiddle_13;
extern notw_codelet fftwi_no_twiddle_14;
extern notw_codelet fftwi_no_twiddle_15;
extern notw_codelet fftwi_no_twiddle_16;
extern notw_codelet fftwi_no_twiddle_32;
extern notw_codelet fftwi_no_twiddle_64;

config_notw fftw_config_notw[] =
{
     NOTW_CODELET(1),
     NOTW_CODELET(2),
     NOTW_CODELET(3),
     NOTW_CODELET(4),
     NOTW_CODELET(5),
     NOTW_CODELET(6),
     NOTW_CODELET(7),
     NOTW_CODELET(8),
     NOTW_CODELET(9),
     NOTW_CODELET(10),
     NOTW_CODELET(11),
     NOTW_CODELET(12),
     NOTW_CODELET(13),
     NOTW_CODELET(14),
     NOTW_CODELET(15),
     NOTW_CODELET(16),
     NOTW_CODELET(32),
     NOTW_CODELET(64),
     {0, 0, (notw_codelet *) 0}
};

config_notw fftwi_config_notw[] =
{
     NOTWI_CODELET(1),
     NOTWI_CODELET(2),
     NOTWI_CODELET(3),
     NOTWI_CODELET(4),
     NOTWI_CODELET(5),
     NOTWI_CODELET(6),
     NOTWI_CODELET(7),
     NOTWI_CODELET(8),
     NOTWI_CODELET(9),
     NOTWI_CODELET(10),
     NOTWI_CODELET(11),
     NOTWI_CODELET(12),
     NOTWI_CODELET(13),
     NOTWI_CODELET(14),
     NOTWI_CODELET(15),
     NOTWI_CODELET(16),
     NOTWI_CODELET(32),
     NOTWI_CODELET(64),
     {0, 0, (notw_codelet *) 0}
};

/* the signature is the same as the size, for now */
#define TWIDDLE_CODELET(x)  { x, x, fftw_twiddle_##x }
#define TWIDDLEI_CODELET(x)  { x, x, fftwi_twiddle_##x }

extern twiddle_codelet fftw_twiddle_2;
extern twiddle_codelet fftw_twiddle_3;
extern twiddle_codelet fftw_twiddle_4;
extern twiddle_codelet fftw_twiddle_5;
extern twiddle_codelet fftw_twiddle_6;
extern twiddle_codelet fftw_twiddle_7;
extern twiddle_codelet fftw_twiddle_8;
extern twiddle_codelet fftw_twiddle_9;
extern twiddle_codelet fftw_twiddle_10;
extern twiddle_codelet fftw_twiddle_16;
extern twiddle_codelet fftw_twiddle_32;
extern twiddle_codelet fftw_twiddle_64;

extern twiddle_codelet fftwi_twiddle_2;
extern twiddle_codelet fftwi_twiddle_3;
extern twiddle_codelet fftwi_twiddle_4;
extern twiddle_codelet fftwi_twiddle_5;
extern twiddle_codelet fftwi_twiddle_6;
extern twiddle_codelet fftwi_twiddle_7;
extern twiddle_codelet fftwi_twiddle_8;
extern twiddle_codelet fftwi_twiddle_9;
extern twiddle_codelet fftwi_twiddle_10;
extern twiddle_codelet fftwi_twiddle_16;
extern twiddle_codelet fftwi_twiddle_32;
extern twiddle_codelet fftwi_twiddle_64;

config_twiddle fftw_config_twiddle[] =
{
     TWIDDLE_CODELET(2),
     TWIDDLE_CODELET(3),
     TWIDDLE_CODELET(4),
     TWIDDLE_CODELET(5),
     TWIDDLE_CODELET(6),
     TWIDDLE_CODELET(7),
     TWIDDLE_CODELET(8),
     TWIDDLE_CODELET(9),
     TWIDDLE_CODELET(10),
     TWIDDLE_CODELET(16),
     TWIDDLE_CODELET(32),
     TWIDDLE_CODELET(64),
     {0, 0, (twiddle_codelet *) 0}
};

config_twiddle fftwi_config_twiddle[] =
{
     TWIDDLEI_CODELET(2),
     TWIDDLEI_CODELET(3),
     TWIDDLEI_CODELET(4),
     TWIDDLEI_CODELET(5),
     TWIDDLEI_CODELET(6),
     TWIDDLEI_CODELET(7),
     TWIDDLEI_CODELET(8),
     TWIDDLEI_CODELET(9),
     TWIDDLEI_CODELET(10),
     TWIDDLEI_CODELET(16),
     TWIDDLEI_CODELET(32),
     TWIDDLEI_CODELET(64),
     {0, 0, (twiddle_codelet *) 0}
};
