/*
 * This code is part of the tng binary trajectory format.
 *
 * Copyright (c) 2010,2013, The GROMACS development team.
 * Copyright (c) 2020, by the GROMACS development team.
 * TNG was orginally written by Magnus Lundborg, Daniel Spångberg and
 * Rossen Apostolov. The API is implemented mainly by Magnus Lundborg,
 * Daniel Spångberg and Anders Gärdenäs.
 *
 * Please see the AUTHORS file for more information.
 *
 * The TNG library is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 *
 * To help us fund future development, we humbly ask that you cite
 * the research papers on the package.
 *
 * Check out http://www.gromacs.org for more information.
 */

/* This code is part of the tng compression routines
 * Written by Daniel Spangberg
 */

#ifndef FIXPOINT_H
#define FIXPOINT_H

#include "../compression/my64bit.h"

/* There are at least 32 bits available in a long. */
typedef unsigned long fix_t;

/* Positive double to 32 bit fixed point value */
fix_t Ptngc_ud_to_fix_t(double d, double max);

/* double to signed 32 bit fixed point value */
fix_t Ptngc_d_to_fix_t(double d, double max);

/* 32 bit fixed point value to positive double */
double Ptngc_fix_t_to_ud(fix_t f, double max);

/* signed 32 bit fixed point value to double */
double Ptngc_fix_t_to_d(fix_t f, double max);

/* Convert a floating point variable to two 32 bit integers with range
   -2.1e9 to 2.1e9 and precision to somewhere around 1e-9. */
void Ptngc_d_to_i32x2(double d, fix_t* hi, fix_t* lo);

/* Convert two 32 bit integers to a floating point variable
   -2.1e9 to 2.1e9 and precision to somewhere around 1e-9. */
double Ptngc_i32x2_to_d(fix_t hi, fix_t lo);

#endif
