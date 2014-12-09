/* This code is part of the tng compression routines.
 *
 * Written by Daniel Spangberg
 * Copyright (c) 2010, 2013, The GROMACS development team.
 *
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the Revised BSD License.
 */

#include <stdio.h>
#include <math.h>
#include "../../include/compression/fixpoint.h"

#define MAX32BIT 4294967295UL
#define MAX31BIT 2147483647UL
#define SIGN32BIT 2147483648UL

/* Conversion routines from / to double precision */

/* Positive double to 32 bit fixed point value */
fix_t Ptngc_ud_to_fix_t(double d, const double max)
{
  fix_t val;
  if (d<0.)
    d=0.;
  if (d>max)
    d=max;
  val=(fix_t)(MAX32BIT*(d/max));
  if (val>MAX32BIT)
    val=MAX32BIT;
  return val;
}

/* double to signed 32 bit fixed point value */
fix_t Ptngc_d_to_fix_t(double d, const double max)
{
  fix_t val;
  int sign=0;
  if (d<0.)
    {
      sign=1;
      d=-d;
    }
  if (d>max)
    d=max;
  val=(fix_t)(MAX31BIT*(d/max));
  if (val>MAX31BIT)
    val=MAX31BIT;
  if (sign)
    val|=SIGN32BIT;
  return val;
}


/* 32 bit fixed point value to positive double */
double Ptngc_fix_t_to_ud(fix_t f, const double max)
{
  return (double)f*(max/MAX32BIT);
}

/* signed 32 bit fixed point value to double */
double Ptngc_fix_t_to_d(fix_t f, const double max)
{
  int sign=0;
  double d;
  if (f&SIGN32BIT)
    {
      sign=1;
      f&=MAX31BIT;
    }
  d=(double)f*(max/MAX31BIT);
  if (sign)
    d=-d;
  return d;
}


/* Convert a floating point variable to two 32 bit integers with range
   -2.1e9 to 2.1e9 and precision to somewhere around 1e-9. */
void Ptngc_d_to_i32x2(double d, fix_t *hi, fix_t *lo)
{
  int sign=0;
  double frac;
  double ent;
  fix_t val,vallo;
  if (d<0.)
    {
      sign=1;
      d=-d;
    }
  /* First the integer part */
  ent=floor(d);
  /* Then the fractional part */
  frac=d-ent;

  val=(fix_t)ent;
  if (sign)
    val|=SIGN32BIT;

  vallo=Ptngc_ud_to_fix_t(frac,1.);

  *hi=val;
  *lo=vallo;
}

/* Convert two 32 bit integers to a floating point variable
   -2.1e9 to 2.1e9 and precision to somewhere around 1e-9. */
double Ptngc_i32x2_to_d(fix_t hi, fix_t lo)
{
  double ent,frac=0.;
  double val=0.;
  int sign=0;
  if (hi&SIGN32BIT)
    {
      sign=1;
      hi&=MAX31BIT;
    }
  ent=(double)hi;
  frac=Ptngc_fix_t_to_ud(lo,1.);
  val=ent+frac;
  if (sign)
    val=-val;
  return val;
}

