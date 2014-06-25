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

/* Multiply two 32 bit unsigned integers returning a 64 bit unsigned value (in two integers) */
void Ptngc_widemul(unsigned int i1, unsigned int i2, unsigned int *ohi, unsigned int *olo);

/* Divide a 64 bit unsigned value in hi:lo with the 32 bit value i and
   return the result in result and the remainder in remainder */
void Ptngc_widediv(unsigned int hi, unsigned int lo, unsigned int i, unsigned int *result, unsigned int *remainder);

/* Add a unsigned int to a largeint. */
void Ptngc_largeint_add(unsigned int v1, unsigned int *largeint, int n);

/* Multiply v1 with largeint_in and return result in largeint_out */
void Ptngc_largeint_mul(unsigned int v1, unsigned int *largeint_in, unsigned int *largeint_out, int n);

/* Return the remainder from dividing largeint_in with v1. Result of the division is returned in largeint_out */
unsigned int Ptngc_largeint_div(unsigned int v1, unsigned int *largeint_in, unsigned int *largeint_out, int n);

#endif
