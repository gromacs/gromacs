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

#ifndef WIDEMULDIV_H
#define WIDEMULDIV_H

/* Add a unsigned int to a largeint. */
void Ptngc_largeint_add(unsigned int v1, unsigned int* largeint, int n);

/* Multiply v1 with largeint_in and return result in largeint_out */
void Ptngc_largeint_mul(unsigned int v1, unsigned int* largeint_in, unsigned int* largeint_out, int n);

/* Return the remainder from dividing largeint_in with v1. Result of the division is returned in largeint_out */
unsigned int Ptngc_largeint_div(unsigned int v1, unsigned int* largeint_in, unsigned int* largeint_out, int n);

#endif
