/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */


#ifndef WIDEMULDIV_H
#define WIDEMULDIV_H

/* Add a unsigned int to a largeint. */
void Ptngc_largeint_add(const unsigned int v1, unsigned int *largeint, const int n);

/* Multiply v1 with largeint_in and return result in largeint_out */
void Ptngc_largeint_mul(const unsigned int v1, unsigned int *largeint_in, unsigned int *largeint_out, const int n);

/* Return the remainder from dividing largeint_in with v1. Result of the division is returned in largeint_out */
unsigned int Ptngc_largeint_div(const unsigned int v1, unsigned int *largeint_in, unsigned int *largeint_out, const int n);

#endif
