/* -*- mode: c; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4; c-file-style: "stroustrup"; -*-
 * 
 *                This source code is part of
 * 
 *                 G   R   O   M   A   C   S
 * 
 *          GROningen MAchine for Chemical Simulations
 * 
 *                        VERSION 4.5
 * Written by David van der Spoel, Erik Lindahl, Berk Hess, and others.
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2008, The GROMACS development team,
 * check out http://www.gromacs.org for more information.
 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * If you want to redistribute modifications, please consider that
 * scientific software is very special. Version control is crucial -
 * bugs must be traceable. We will be happy to consider code for
 * inclusion in the official distribution, but derived work must not
 * be called official GROMACS. Details are found in the README & COPYING
 * files - if they are missing, get the official version at www.gromacs.org.
 * 
 * To help us fund GROMACS development, we humbly ask that you cite
 * the papers on the package - you can find them in the top README file.
 * 
 * For more info, check our website at http://www.gromacs.org
 * 
 * And Hey:
 * Groningen Machine for Chemical Simulation
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


/* slater_S_integrals.c (c) 2008 Paul J. van Maaren and David van der Spoel */
#include <stdio.h>
#include <math.h>
#include "slater_S_integrals.h"

#define Power(x, y)     (pow((double)(x), (double)(y)))
#define Pi              3.14159265358979323846264
#define E               2.71828182845904523536029

static double Slater_1S_1S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-24.0 + 24.0*Power(E,2.0*rij*xii) - 33.0*rij*xii - 18.0*Power(rij,2.0)*Power(xii,2.0) - 
        4.0*Power(rij,3.0)*Power(xii,3.0))/(24.*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),3.0) + 
        Power(E,2.0*rij*xij)*Power(xij,4.0)*
         (-3.0*Power(xii,2.0) - rij*Power(xii,3.0) + Power(xij,2.0) + 
           rij*xii*Power(xij,2.0)) - 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (Power(xii,2.0)*(1.0 + rij*xij) - Power(xij,2.0)*(3.0 + rij*xij)))/
      (Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),3.0))
    ;
  }
  return S;
}

static double Slater_1S_2S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-240.0 + 240.0*Power(E,2.0*rij*xii) - 375.0*rij*xii - 
        270.0*Power(rij,2.0)*Power(xii,2.0) - 115.0*Power(rij,3.0)*Power(xii,3.0) - 
        30.0*Power(rij,4.0)*Power(xii,4.0) - 4.0*Power(rij,5.0)*Power(xii,5.0))/
      (240.*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (6.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),5.0) + 
        6.0*Power(E,2.0*rij*xij)*Power(xij,6.0)*
         (-4.0*Power(xii,4.0) - rij*Power(xii,5.0) - 5.0*Power(xii,2.0)*Power(xij,2.0) + 
           Power(xij,4.0) + rij*xii*Power(xij,4.0)) - 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (Power(xii,6.0)*(6.0 + 9.0*rij*xij + 6.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,3.0)*Power(xij,3.0)) - 
           3.0*Power(xii,4.0)*Power(xij,2.0)*
            (10.0 + 15.0*rij*xij + 10.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,3.0)*Power(xij,3.0)) + 
           3.0*Power(xii,2.0)*Power(xij,4.0)*
            (20.0 + 33.0*rij*xij + 14.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,3.0)*Power(xij,3.0)) - 
           Power(xij,6.0)*(84.0 + 63.0*rij*xij + 18.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,3.0)*Power(xij,3.0))))/
      (6.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),5.0))
    ;
  }
  return S;
}

static double Slater_1S_3S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-120960.0 + 120960.0*Power(E,2.0*rij*xii) - 203175.0*rij*xii - 
        164430.0*Power(rij,2.0)*Power(xii,2.0) - 84420.0*Power(rij,3.0)*Power(xii,3.0) - 
        30240.0*Power(rij,4.0)*Power(xii,4.0) - 7728.0*Power(rij,5.0)*Power(xii,5.0) - 
        1344.0*Power(rij,6.0)*Power(xii,6.0) - 128.0*Power(rij,7.0)*Power(xii,7.0))/
      (120960.*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (45.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),7.0) + 
        15.0*Power(E,2.0*rij*xij)*Power(xij,8.0)*
         (-15.0*Power(xii,6.0) - 3.0*rij*Power(xii,7.0) - 
           63.0*Power(xii,4.0)*Power(xij,2.0) - 7.0*rij*Power(xii,5.0)*Power(xij,2.0) - 
           21.0*Power(xii,2.0)*Power(xij,4.0) + 7.0*rij*Power(xii,3.0)*Power(xij,4.0) + 
           3.0*Power(xij,6.0) + 3.0*rij*xii*Power(xij,6.0)) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-10.0*Power(xii,2.0)*Power(xij,8.0)*
            (135.0 + 333.0*rij*xij + 228.0*Power(rij,2.0)*Power(xij,2.0) + 
              75.0*Power(rij,3.0)*Power(xij,3.0) + 13.0*Power(rij,4.0)*Power(xij,4.0) + 
              Power(rij,5.0)*Power(xij,5.0)) + 
           2.0*Power(xij,10.0)*(945.0 + 945.0*rij*xij + 420.0*Power(rij,2.0)*Power(xij,2.0) + 
              105.0*Power(rij,3.0)*Power(xij,3.0) + 15.0*Power(rij,4.0)*Power(xij,4.0) + 
              Power(rij,5.0)*Power(xij,5.0)) - 
           Power(xii,10.0)*(45.0 + 75.0*rij*xij + 60.0*Power(rij,2.0)*Power(xij,2.0) + 
              30.0*Power(rij,3.0)*Power(xij,3.0) + 10.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) + 
           5.0*Power(xii,8.0)*Power(xij,2.0)*
            (63.0 + 105.0*rij*xij + 84.0*Power(rij,2.0)*Power(xij,2.0) + 
              42.0*Power(rij,3.0)*Power(xij,3.0) + 14.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) - 
           5.0*Power(xii,6.0)*Power(xij,4.0)*
            (189.0 + 315.0*rij*xij + 252.0*Power(rij,2.0)*Power(xij,2.0) + 
              132.0*Power(rij,3.0)*Power(xij,3.0) + 36.0*Power(rij,4.0)*Power(xij,4.0) + 
              4.0*Power(rij,5.0)*Power(xij,5.0)) + 
           5.0*Power(xii,4.0)*Power(xij,6.0)*
            (315.0 + 513.0*rij*xij + 468.0*Power(rij,2.0)*Power(xij,2.0) + 
              204.0*Power(rij,3.0)*Power(xij,3.0) + 44.0*Power(rij,4.0)*Power(xij,4.0) + 
              4.0*Power(rij,5.0)*Power(xij,5.0))))/
      (45.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),7.0))
    ;
  }
  return S;
}

static double Slater_1S_4S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-2903040.0 + 2903040.0*Power(E,2.0*rij*xii) - 5088825.0*rij*xii - 
        4371570.0*Power(rij,2.0)*Power(xii,2.0) - 2439990.0*Power(rij,3.0)*Power(xii,3.0) - 
        986580.0*Power(rij,4.0)*Power(xii,4.0) - 303912.0*Power(rij,5.0)*Power(xii,5.0) - 
        72576.0*Power(rij,6.0)*Power(xii,6.0) - 13248.0*Power(rij,7.0)*Power(xii,7.0) - 
        1728.0*Power(rij,8.0)*Power(xii,8.0) - 128.0*Power(rij,9.0)*Power(xii,9.0))/
      (2.90304e6*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (1260.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),9.0) + 
        1260.0*Power(E,2.0*rij*xij)*Power(xij,10.0)*
         (-6.0*Power(xii,8.0) - rij*Power(xii,9.0) - 51.0*Power(xii,6.0)*Power(xij,2.0) - 
           6.0*rij*Power(xii,7.0)*Power(xij,2.0) - 63.0*Power(xii,4.0)*Power(xij,4.0) - 
           9.0*Power(xii,2.0)*Power(xij,6.0) + 6.0*rij*Power(xii,3.0)*Power(xij,6.0) + 
           Power(xij,8.0) + rij*xii*Power(xij,8.0)) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-42.0*Power(xii,10.0)*Power(xij,4.0)*
            (1080.0 + 1890.0*rij*xij + 1620.0*Power(rij,2.0)*Power(xij,2.0) + 
              900.0*Power(rij,3.0)*Power(xij,3.0) + 360.0*Power(rij,4.0)*Power(xij,4.0) + 
              111.0*Power(rij,5.0)*Power(xij,5.0) + 22.0*Power(rij,6.0)*Power(xij,6.0) + 
              2.0*Power(rij,7.0)*Power(xij,7.0)) + 
           70.0*Power(xii,8.0)*Power(xij,6.0)*
            (1512.0 + 2646.0*rij*xij + 2268.0*Power(rij,2.0)*Power(xij,2.0) + 
              1248.0*Power(rij,3.0)*Power(xij,3.0) + 528.0*Power(rij,4.0)*Power(xij,4.0) + 
              153.0*Power(rij,5.0)*Power(xij,5.0) + 26.0*Power(rij,6.0)*Power(xij,6.0) + 
              2.0*Power(rij,7.0)*Power(xij,7.0)) - 
           14.0*Power(xii,2.0)*Power(xij,12.0)*
            (2970.0 + 16335.0*rij*xij + 15390.0*Power(rij,2.0)*Power(xij,2.0) + 
              7110.0*Power(rij,3.0)*Power(xij,3.0) + 1980.0*Power(rij,4.0)*Power(xij,4.0) + 
              351.0*Power(rij,5.0)*Power(xij,5.0) + 38.0*Power(rij,6.0)*Power(xij,6.0) + 
              2.0*Power(rij,7.0)*Power(xij,7.0)) + 
           2.0*Power(xij,14.0)*(62370.0 + 72765.0*rij*xij + 
              39690.0*Power(rij,2.0)*Power(xij,2.0) + 
              13230.0*Power(rij,3.0)*Power(xij,3.0) + 
              2940.0*Power(rij,4.0)*Power(xij,4.0) + 441.0*Power(rij,5.0)*Power(xij,5.0) + 
              42.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,7.0)*Power(xij,7.0)) - 
           Power(xii,14.0)*(1260.0 + 2205.0*rij*xij + 
              1890.0*Power(rij,2.0)*Power(xij,2.0) + 1050.0*Power(rij,3.0)*Power(xij,3.0) + 
              420.0*Power(rij,4.0)*Power(xij,4.0) + 126.0*Power(rij,5.0)*Power(xij,5.0) + 
              28.0*Power(rij,6.0)*Power(xij,6.0) + 4.0*Power(rij,7.0)*Power(xij,7.0)) + 
           7.0*Power(xii,12.0)*Power(xij,2.0)*
            (1620.0 + 2835.0*rij*xij + 2430.0*Power(rij,2.0)*Power(xij,2.0) + 
              1350.0*Power(rij,3.0)*Power(xij,3.0) + 540.0*Power(rij,4.0)*Power(xij,4.0) + 
              162.0*Power(rij,5.0)*Power(xij,5.0) + 36.0*Power(rij,6.0)*Power(xij,6.0) + 
              4.0*Power(rij,7.0)*Power(xij,7.0)) - 
           35.0*Power(xii,6.0)*Power(xij,8.0)*
            (4536.0 + 7983.0*rij*xij + 6534.0*Power(rij,2.0)*Power(xij,2.0) + 
              4014.0*Power(rij,3.0)*Power(xij,3.0) + 1644.0*Power(rij,4.0)*Power(xij,4.0) + 
              414.0*Power(rij,5.0)*Power(xij,5.0) + 60.0*Power(rij,6.0)*Power(xij,6.0) + 
              4.0*Power(rij,7.0)*Power(xij,7.0)) + 
           21.0*Power(xii,4.0)*Power(xij,10.0)*
            (7920.0 + 11385.0*rij*xij + 12330.0*Power(rij,2.0)*Power(xij,2.0) + 
              7410.0*Power(rij,3.0)*Power(xij,3.0) + 2580.0*Power(rij,4.0)*Power(xij,4.0) + 
              546.0*Power(rij,5.0)*Power(xij,5.0) + 68.0*Power(rij,6.0)*Power(xij,6.0) + 
              4.0*Power(rij,7.0)*Power(xij,7.0))))/
      (1260.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),9.0))
    ;
  }
  return S;
}

static double Slater_1S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-1596672000.0 + 1596672000.0*Power(E,2.0*rij*xii) - 2875101075.0*rij*xii - 
        2556858150.0*Power(rij,2.0)*Power(xii,2.0) - 
        1492929900.0*Power(rij,3.0)*Power(xii,3.0) - 
        641163600.0*Power(rij,4.0)*Power(xii,4.0) - 
        214719120.0*Power(rij,5.0)*Power(xii,5.0) - 
        57879360.0*Power(rij,6.0)*Power(xii,6.0) - 
        12735360.0*Power(rij,7.0)*Power(xii,7.0) - 2280960.0*Power(rij,8.0)*Power(xii,8.0) - 
        323840.0*Power(rij,9.0)*Power(xii,9.0) - 33792.0*Power(rij,10.0)*Power(xii,10.0) - 
        2048.0*Power(rij,11.0)*Power(xii,11.0))/(1.596672e9*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (14175.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),11.0) + 
        2835.0*Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-35.0*Power(xii,10.0) - 5.0*rij*Power(xii,11.0) - 
           495.0*Power(xii,8.0)*Power(xij,2.0) - 55.0*rij*Power(xii,9.0)*Power(xij,2.0) - 
           1254.0*Power(xii,6.0)*Power(xij,4.0) - 66.0*rij*Power(xii,7.0)*Power(xij,4.0) - 
           726.0*Power(xii,4.0)*Power(xij,6.0) + 66.0*rij*Power(xii,5.0)*Power(xij,6.0) - 
           55.0*Power(xii,2.0)*Power(xij,8.0) + 55.0*rij*Power(xii,3.0)*Power(xij,8.0) + 
           5.0*Power(xij,10.0) + 5.0*rij*xii*Power(xij,10.0)) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-(Power(xii,18.0)*(14175.0 + 25515.0*rij*xij + 
                22680.0*Power(rij,2.0)*Power(xij,2.0) + 
                13230.0*Power(rij,3.0)*Power(xij,3.0) + 
                5670.0*Power(rij,4.0)*Power(xij,4.0) + 
                1890.0*Power(rij,5.0)*Power(xij,5.0) + 
                504.0*Power(rij,6.0)*Power(xij,6.0) + 108.0*Power(rij,7.0)*Power(xij,7.0) + 
                18.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0))) + 
           9.0*Power(xii,16.0)*Power(xij,2.0)*
            (17325.0 + 31185.0*rij*xij + 27720.0*Power(rij,2.0)*Power(xij,2.0) + 
              16170.0*Power(rij,3.0)*Power(xij,3.0) + 
              6930.0*Power(rij,4.0)*Power(xij,4.0) + 2310.0*Power(rij,5.0)*Power(xij,5.0) + 
              616.0*Power(rij,6.0)*Power(xij,6.0) + 132.0*Power(rij,7.0)*Power(xij,7.0) + 
              22.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) - 
           126.0*Power(xii,10.0)*Power(xij,8.0)*
            (37125.0 + 66825.0*rij*xij + 59400.0*Power(rij,2.0)*Power(xij,2.0) + 
              34725.0*Power(rij,3.0)*Power(xij,3.0) + 
              14625.0*Power(rij,4.0)*Power(xij,4.0) + 
              5043.0*Power(rij,5.0)*Power(xij,5.0) + 1396.0*Power(rij,6.0)*Power(xij,6.0) + 
              276.0*Power(rij,7.0)*Power(xij,7.0) + 34.0*Power(rij,8.0)*Power(xij,8.0) + 
              2.0*Power(rij,9.0)*Power(xij,9.0)) + 
           126.0*Power(xii,8.0)*Power(xij,10.0)*
            (51975.0 + 93420.0*rij*xij + 84240.0*Power(rij,2.0)*Power(xij,2.0) + 
              46815.0*Power(rij,3.0)*Power(xij,3.0) + 
              20835.0*Power(rij,4.0)*Power(xij,4.0) + 
              7485.0*Power(rij,5.0)*Power(xij,5.0) + 1964.0*Power(rij,6.0)*Power(xij,6.0) + 
              348.0*Power(rij,7.0)*Power(xij,7.0) + 38.0*Power(rij,8.0)*Power(xij,8.0) + 
              2.0*Power(rij,9.0)*Power(xij,9.0)) - 
           9.0*Power(xii,2.0)*Power(xij,16.0)*
            (-135135.0 + 405405.0*rij*xij + 582120.0*Power(rij,2.0)*Power(xij,2.0) + 
              346500.0*Power(rij,3.0)*Power(xij,3.0) + 
              124740.0*Power(rij,4.0)*Power(xij,4.0) + 
              30492.0*Power(rij,5.0)*Power(xij,5.0) + 
              5264.0*Power(rij,6.0)*Power(xij,6.0) + 636.0*Power(rij,7.0)*Power(xij,7.0) + 
              50.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) + 
           Power(xij,18.0)*(2837835.0 + 3648645.0*rij*xij + 
              2245320.0*Power(rij,2.0)*Power(xij,2.0) + 
              873180.0*Power(rij,3.0)*Power(xij,3.0) + 
              238140.0*Power(rij,4.0)*Power(xij,4.0) + 
              47628.0*Power(rij,5.0)*Power(xij,5.0) + 
              7056.0*Power(rij,6.0)*Power(xij,6.0) + 756.0*Power(rij,7.0)*Power(xij,7.0) + 
              54.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) - 
           9.0*Power(xii,14.0)*Power(xij,4.0)*
            (86625.0 + 155925.0*rij*xij + 138600.0*Power(rij,2.0)*Power(xij,2.0) + 
              80850.0*Power(rij,3.0)*Power(xij,3.0) + 
              34650.0*Power(rij,4.0)*Power(xij,4.0) + 
              11550.0*Power(rij,5.0)*Power(xij,5.0) + 
              3080.0*Power(rij,6.0)*Power(xij,6.0) + 672.0*Power(rij,7.0)*Power(xij,7.0) + 
              104.0*Power(rij,8.0)*Power(xij,8.0) + 8.0*Power(rij,9.0)*Power(xij,9.0)) + 
           21.0*Power(xii,12.0)*Power(xij,6.0)*
            (111375.0 + 200475.0*rij*xij + 178200.0*Power(rij,2.0)*Power(xij,2.0) + 
              103950.0*Power(rij,3.0)*Power(xij,3.0) + 
              44550.0*Power(rij,4.0)*Power(xij,4.0) + 
              14778.0*Power(rij,5.0)*Power(xij,5.0) + 
              4056.0*Power(rij,6.0)*Power(xij,6.0) + 864.0*Power(rij,7.0)*Power(xij,7.0) + 
              120.0*Power(rij,8.0)*Power(xij,8.0) + 8.0*Power(rij,9.0)*Power(xij,9.0)) - 
           21.0*Power(xii,6.0)*Power(xij,12.0)*
            (307125.0 + 594945.0*rij*xij + 456840.0*Power(rij,2.0)*Power(xij,2.0) + 
              281790.0*Power(rij,3.0)*Power(xij,3.0) + 
              137430.0*Power(rij,4.0)*Power(xij,4.0) + 
              47250.0*Power(rij,5.0)*Power(xij,5.0) + 
              11064.0*Power(rij,6.0)*Power(xij,6.0) + 
              1728.0*Power(rij,7.0)*Power(xij,7.0) + 168.0*Power(rij,8.0)*Power(xij,8.0) + 
              8.0*Power(rij,9.0)*Power(xij,9.0)) + 
           9.0*Power(xii,4.0)*Power(xij,14.0)*
            (675675.0 + 675675.0*rij*xij + 748440.0*Power(rij,2.0)*Power(xij,2.0) + 
              561330.0*Power(rij,3.0)*Power(xij,3.0) + 
              256410.0*Power(rij,4.0)*Power(xij,4.0) + 
              76230.0*Power(rij,5.0)*Power(xij,5.0) + 
              15400.0*Power(rij,6.0)*Power(xij,6.0) + 2112.0*Power(rij,7.0)*Power(xij,7.0) + 
              184.0*Power(rij,8.0)*Power(xij,8.0) + 8.0*Power(rij,9.0)*Power(xij,9.0))))/
      (14175.*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),11.0))
    ;
  }
  return S;
}

static double Slater_1S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-74724249600.0 + 74724249600.0*Power(E,2.0*rij*xii) - 137006619750.0*rij*xii - 
        124564740300.0*Power(rij,2.0)*Power(xii,2.0) - 
        74754654975.0*Power(rij,3.0)*Power(xii,3.0) - 
        33239155950.0*Power(rij,4.0)*Power(xii,4.0) - 
        11644853220.0*Power(rij,5.0)*Power(xii,5.0) - 
        3334050720.0*Power(rij,6.0)*Power(xii,6.0) - 
        797528160.0*Power(rij,7.0)*Power(xii,7.0) - 
        161235360.0*Power(rij,8.0)*Power(xii,8.0) - 
        27593280.0*Power(rij,9.0)*Power(xii,9.0) - 
        3953664.0*Power(rij,10.0)*Power(xii,10.0) - 
        459264.0*Power(rij,11.0)*Power(xii,11.0) - 39936.0*Power(rij,12.0)*Power(xii,12.0) - 
        2048.0*Power(rij,13.0)*Power(xii,13.0))/(7.47242496e10*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (935550.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),13.0) + 
        311850.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-24.0*Power(xii,12.0) - 3.0*rij*Power(xii,13.0) - 
           507.0*Power(xii,10.0)*Power(xij,2.0) - 52.0*rij*Power(xii,11.0)*Power(xij,2.0) - 
           2145.0*Power(xii,8.0)*Power(xij,4.0) - 143.0*rij*Power(xii,9.0)*Power(xij,4.0) - 
           2574.0*Power(xii,6.0)*Power(xij,6.0) - 858.0*Power(xii,4.0)*Power(xij,8.0) + 
           143.0*rij*Power(xii,5.0)*Power(xij,8.0) - 39.0*Power(xii,2.0)*Power(xij,10.0) + 
           52.0*rij*Power(xii,3.0)*Power(xij,10.0) + 3.0*Power(xij,12.0) + 
           3.0*rij*xii*Power(xij,12.0)) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-110.0*Power(xii,18.0)*Power(xij,4.0)*
            (663390.0 + 1216215.0*rij*xij + 1105650.0*Power(rij,2.0)*Power(xij,2.0) + 
              663390.0*Power(rij,3.0)*Power(xij,3.0) + 
              294840.0*Power(rij,4.0)*Power(xij,4.0) + 
              103194.0*Power(rij,5.0)*Power(xij,5.0) + 
              29484.0*Power(rij,6.0)*Power(xij,6.0) + 
              7020.0*Power(rij,7.0)*Power(xij,7.0) + 1404.0*Power(rij,8.0)*Power(xij,8.0) + 
              237.0*Power(rij,9.0)*Power(xij,9.0) + 30.0*Power(rij,10.0)*Power(xij,10.0) + 
              2.0*Power(rij,11.0)*Power(xij,11.0)) + 
           330.0*Power(xii,16.0)*Power(xij,6.0)*
            (810810.0 + 1486485.0*rij*xij + 1351350.0*Power(rij,2.0)*Power(xij,2.0) + 
              810810.0*Power(rij,3.0)*Power(xij,3.0) + 
              360360.0*Power(rij,4.0)*Power(xij,4.0) + 
              126126.0*Power(rij,5.0)*Power(xij,5.0) + 
              36036.0*Power(rij,6.0)*Power(xij,6.0) + 
              8556.0*Power(rij,7.0)*Power(xij,7.0) + 1740.0*Power(rij,8.0)*Power(xij,8.0) + 
              291.0*Power(rij,9.0)*Power(xij,9.0) + 34.0*Power(rij,10.0)*Power(xij,10.0) + 
              2.0*Power(rij,11.0)*Power(xij,11.0)) - 
           330.0*Power(xii,6.0)*Power(xij,16.0)*
            (3169530.0 + 7960680.0*rij*xij + 5798520.0*Power(rij,2.0)*Power(xij,2.0) + 
              3144960.0*Power(rij,3.0)*Power(xij,3.0) + 
              1572480.0*Power(rij,4.0)*Power(xij,4.0) + 
              638001.0*Power(rij,5.0)*Power(xij,5.0) + 
              191646.0*Power(rij,6.0)*Power(xij,6.0) + 
              41886.0*Power(rij,7.0)*Power(xij,7.0) + 
              6630.0*Power(rij,8.0)*Power(xij,8.0) + 741.0*Power(rij,9.0)*Power(xij,9.0) + 
              54.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) + 
           110.0*Power(xii,4.0)*Power(xij,18.0)*
            (12162150.0 + 8108100.0*rij*xij + 6486480.0*Power(rij,2.0)*Power(xij,2.0) + 
              5675670.0*Power(rij,3.0)*Power(xij,3.0) + 
              3243240.0*Power(rij,4.0)*Power(xij,4.0) + 
              1216215.0*Power(rij,5.0)*Power(xij,5.0) + 
              319410.0*Power(rij,6.0)*Power(xij,6.0) + 
              61074.0*Power(rij,7.0)*Power(xij,7.0) + 
              8586.0*Power(rij,8.0)*Power(xij,8.0) + 867.0*Power(rij,9.0)*Power(xij,9.0) + 
              58.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) - 
           Power(xii,22.0)*(935550.0 + 1715175.0*rij*xij + 
              1559250.0*Power(rij,2.0)*Power(xij,2.0) + 
              935550.0*Power(rij,3.0)*Power(xij,3.0) + 
              415800.0*Power(rij,4.0)*Power(xij,4.0) + 
              145530.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              9900.0*Power(rij,7.0)*Power(xij,7.0) + 1980.0*Power(rij,8.0)*Power(xij,8.0) + 
              330.0*Power(rij,9.0)*Power(xij,9.0) + 44.0*Power(rij,10.0)*Power(xij,10.0) + 
              4.0*Power(rij,11.0)*Power(xij,11.0)) + 
           11.0*Power(xii,20.0)*Power(xij,2.0)*
            (1105650.0 + 2027025.0*rij*xij + 1842750.0*Power(rij,2.0)*Power(xij,2.0) + 
              1105650.0*Power(rij,3.0)*Power(xij,3.0) + 
              491400.0*Power(rij,4.0)*Power(xij,4.0) + 
              171990.0*Power(rij,5.0)*Power(xij,5.0) + 
              49140.0*Power(rij,6.0)*Power(xij,6.0) + 
              11700.0*Power(rij,7.0)*Power(xij,7.0) + 
              2340.0*Power(rij,8.0)*Power(xij,8.0) + 390.0*Power(rij,9.0)*Power(xij,9.0) + 
              52.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           11.0*Power(xii,2.0)*Power(xij,20.0)*
            (-48648600.0 + 2027025.0*rij*xij + 44594550.0*Power(rij,2.0)*Power(xij,2.0) + 
              36486450.0*Power(rij,3.0)*Power(xij,3.0) + 
              16216200.0*Power(rij,4.0)*Power(xij,4.0) + 
              4864860.0*Power(rij,5.0)*Power(xij,5.0) + 
              1065960.0*Power(rij,6.0)*Power(xij,6.0) + 
              176040.0*Power(rij,7.0)*Power(xij,7.0) + 
              21960.0*Power(rij,8.0)*Power(xij,8.0) + 
              2010.0*Power(rij,9.0)*Power(xij,9.0) + 
              124.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) + 
           Power(xij,22.0)*(340540200.0 + 468242775.0*rij*xij + 
              312161850.0*Power(rij,2.0)*Power(xij,2.0) + 
              133783650.0*Power(rij,3.0)*Power(xij,3.0) + 
              41164200.0*Power(rij,4.0)*Power(xij,4.0) + 
              9604980.0*Power(rij,5.0)*Power(xij,5.0) + 
              1746360.0*Power(rij,6.0)*Power(xij,6.0) + 
              249480.0*Power(rij,7.0)*Power(xij,7.0) + 
              27720.0*Power(rij,8.0)*Power(xij,8.0) + 
              2310.0*Power(rij,9.0)*Power(xij,9.0) + 
              132.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           165.0*Power(xii,14.0)*Power(xij,8.0)*
            (4054050.0 + 7432425.0*rij*xij + 6756750.0*Power(rij,2.0)*Power(xij,2.0) + 
              4054050.0*Power(rij,3.0)*Power(xij,3.0) + 
              1801800.0*Power(rij,4.0)*Power(xij,4.0) + 
              631260.0*Power(rij,5.0)*Power(xij,5.0) + 
              178920.0*Power(rij,6.0)*Power(xij,6.0) + 
              43176.0*Power(rij,7.0)*Power(xij,7.0) + 
              8904.0*Power(rij,8.0)*Power(xij,8.0) + 1428.0*Power(rij,9.0)*Power(xij,9.0) + 
              152.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           231.0*Power(xii,12.0)*Power(xij,10.0)*
            (5212350.0 + 9555975.0*rij*xij + 8687250.0*Power(rij,2.0)*Power(xij,2.0) + 
              5209650.0*Power(rij,3.0)*Power(xij,3.0) + 
              2327400.0*Power(rij,4.0)*Power(xij,4.0) + 
              801540.0*Power(rij,5.0)*Power(xij,5.0) + 
              230040.0*Power(rij,6.0)*Power(xij,6.0) + 
              57240.0*Power(rij,7.0)*Power(xij,7.0) + 
              11640.0*Power(rij,8.0)*Power(xij,8.0) + 
              1740.0*Power(rij,9.0)*Power(xij,9.0) + 
              168.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) - 
           231.0*Power(xii,10.0)*Power(xij,12.0)*
            (6949800.0 + 12746025.0*rij*xij + 11535750.0*Power(rij,2.0)*Power(xij,2.0) + 
              7056450.0*Power(rij,3.0)*Power(xij,3.0) + 
              3040200.0*Power(rij,4.0)*Power(xij,4.0) + 
              1051920.0*Power(rij,5.0)*Power(xij,5.0) + 
              316800.0*Power(rij,6.0)*Power(xij,6.0) + 
              79680.0*Power(rij,7.0)*Power(xij,7.0) + 
              15360.0*Power(rij,8.0)*Power(xij,8.0) + 
              2100.0*Power(rij,9.0)*Power(xij,9.0) + 
              184.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           165.0*Power(xii,8.0)*Power(xij,14.0)*
            (9775080.0 + 17424855.0*rij*xij + 17019450.0*Power(rij,2.0)*Power(xij,2.0) + 
              9519930.0*Power(rij,3.0)*Power(xij,3.0) + 
              4059720.0*Power(rij,4.0)*Power(xij,4.0) + 
              1519056.0*Power(rij,5.0)*Power(xij,5.0) + 
              475776.0*Power(rij,6.0)*Power(xij,6.0) + 
              114720.0*Power(rij,7.0)*Power(xij,7.0) + 
              20256.0*Power(rij,8.0)*Power(xij,8.0) + 2508.0*Power(rij,9.0)*Power(xij,9.0) + 
              200.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0))))/
      (935550.*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),13.0))
    ;
  }
  return S;
}

static double Slater_2S_2S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-80640.0 + 80640.0*Power(E,2.0*rij*xii) - 131985.0*rij*xii - 
        102690.0*Power(rij,2.0)*Power(xii,2.0) - 49980.0*Power(rij,3.0)*Power(xii,3.0) - 
        16800.0*Power(rij,4.0)*Power(xii,4.0) - 4032.0*Power(rij,5.0)*Power(xii,5.0) - 
        672.0*Power(rij,6.0)*Power(xii,6.0) - 64.0*Power(rij,7.0)*Power(xii,7.0))/
      (80640.*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (6.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),7.0) - 
        Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (21.0*Power(xii,4.0)*Power(xij,4.0)*
            (6.0 + 11.0*rij*xij + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           2.0*Power(xij,8.0)*(90.0 + 54.0*rij*xij + 12.0*Power(rij,2.0)*Power(xij,2.0) + 
              Power(rij,3.0)*Power(xij,3.0)) + 
           Power(xii,8.0)*(6.0 + 9.0*rij*xij + 6.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,3.0)*Power(xij,3.0)) + 
           Power(xii,2.0)*Power(xij,6.0)*
            (-390.0 - 69.0*rij*xij + 18.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,3.0)*Power(xij,3.0)) - 
           Power(xii,6.0)*Power(xij,2.0)*
            (42.0 + 63.0*rij*xij + 42.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,3.0)*Power(xij,3.0))) + 
        Power(E,2.0*rij*xij)*Power(xij,6.0)*
         (-24.0*Power(rij,2.0)*Power(xii,10.0) - 2.0*Power(rij,3.0)*Power(xii,11.0) - 
           69.0*rij*Power(xii,7.0)*Power(xij,2.0) + 6.0*Power(xij,8.0) + 
           9.0*rij*xii*Power(xij,8.0) + 
           4.0*rij*Power(xii,9.0)*(-27.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           18.0*Power(xii,8.0)*(-10.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*Power(xii,2.0)*Power(xij,6.0)*(-7.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           42.0*Power(xii,4.0)*Power(xij,4.0)*(-3.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           rij*Power(xii,3.0)*Power(xij,6.0)*(-63.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*Power(xii,6.0)*Power(xij,2.0)*(-65.0 + 7.0*Power(rij,2.0)*Power(xij,2.0)) + 
           Power(xii,5.0)*(231.0*rij*Power(xij,4.0) - 4.0*Power(rij,3.0)*Power(xij,6.0))))/
      (6.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),7.0))
    ;
  }
  return S;
}

static double Slater_2S_3S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-4354560.0 + 4354560.0*Power(E,2.0*rij*xii) - 7430535.0*rij*xii - 
        6151950.0*Power(rij,2.0)*Power(xii,2.0) - 3275370.0*Power(rij,3.0)*Power(xii,3.0) - 
        1251180.0*Power(rij,4.0)*Power(xii,4.0) - 361368.0*Power(rij,5.0)*Power(xii,5.0) - 
        80640.0*Power(rij,6.0)*Power(xii,6.0) - 13824.0*Power(rij,7.0)*Power(xii,7.0) - 
        1728.0*Power(rij,8.0)*Power(xii,8.0) - 128.0*Power(rij,9.0)*Power(xii,9.0))/
      (4.35456e6*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (90.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),9.0) + 
        5.0*Power(E,2.0*rij*xij)*Power(xij,8.0)*
         (-90.0*Power(rij,2.0)*Power(xii,12.0) - 6.0*Power(rij,3.0)*Power(xii,13.0) + 
           18.0*Power(xij,10.0) + 27.0*rij*xii*Power(xij,10.0) + 
           18.0*Power(xii,2.0)*Power(xij,8.0)*(-9.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           162.0*Power(xii,4.0)*Power(xij,6.0)*(-4.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           198.0*Power(xii,10.0)*(5.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           108.0*Power(xii,6.0)*Power(xij,4.0)*(36.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           2.0*rij*Power(xii,5.0)*Power(xij,6.0)*(675.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           18.0*rij*Power(xii,7.0)*Power(xij,4.0)*
            (-81.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*rij*Power(xii,3.0)*Power(xij,8.0)*
            (-81.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           rij*Power(xii,11.0)*(495.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           9.0*rij*Power(xii,9.0)*Power(xij,2.0)*
            (-233.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*Power(xii,8.0)*Power(xij,2.0)*(-1063.0 + 90.0*Power(rij,2.0)*Power(xij,2.0))) - 
        2.0*Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (-90.0*Power(xii,6.0)*Power(xij,6.0)*
            (42.0 + 65.0*rij*xij + 76.0*Power(rij,2.0)*Power(xij,2.0) + 
              22.0*Power(rij,3.0)*Power(xij,3.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) - 
           2.0*Power(xij,12.0)*(2970.0 + 2475.0*rij*xij + 
              900.0*Power(rij,2.0)*Power(xij,2.0) + 180.0*Power(rij,3.0)*Power(xij,3.0) + 
              20.0*Power(rij,4.0)*Power(xij,4.0) + Power(rij,5.0)*Power(xij,5.0)) + 
           10.0*Power(xii,8.0)*Power(xij,4.0)*
            (162.0 + 270.0*rij*xij + 216.0*Power(rij,2.0)*Power(xij,2.0) + 
              122.0*Power(rij,3.0)*Power(xij,3.0) + 22.0*Power(rij,4.0)*Power(xij,4.0) + 
              Power(rij,5.0)*Power(xij,5.0)) - 
           5.0*Power(xii,4.0)*Power(xij,8.0)*
            (-639.0 - 3555.0*rij*xij - 1452.0*Power(rij,2.0)*Power(xij,2.0) - 
              174.0*Power(rij,3.0)*Power(xij,3.0) + 6.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) + 
           Power(xii,12.0)*(45.0 + 75.0*rij*xij + 60.0*Power(rij,2.0)*Power(xij,2.0) + 
              30.0*Power(rij,3.0)*Power(xij,3.0) + 10.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) - 
           Power(xii,10.0)*Power(xij,2.0)*
            (405.0 + 675.0*rij*xij + 540.0*Power(rij,2.0)*Power(xij,2.0) + 
              270.0*Power(rij,3.0)*Power(xij,3.0) + 90.0*Power(rij,4.0)*Power(xij,4.0) + 
              8.0*Power(rij,5.0)*Power(xij,5.0)) + 
           Power(xii,2.0)*Power(xij,10.0)*
            (-21615.0 - 9075.0*rij*xij - 300.0*Power(rij,2.0)*Power(xij,2.0) + 
              490.0*Power(rij,3.0)*Power(xij,3.0) + 110.0*Power(rij,4.0)*Power(xij,4.0) + 
              8.0*Power(rij,5.0)*Power(xij,5.0))))/
      (90.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),9.0))
    ;
  }
  return S;
}

static double Slater_2S_4S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-638668800.0 + 638668800.0*Power(E,2.0*rij*xii) - 1125310725.0*rij*xii - 
        973283850.0*Power(rij,2.0)*Power(xii,2.0) - 
        549063900.0*Power(rij,3.0)*Power(xii,3.0) - 
        226195200.0*Power(rij,4.0)*Power(xii,4.0) - 
        72099720.0*Power(rij,5.0)*Power(xii,5.0) - 
        18350640.0*Power(rij,6.0)*Power(xii,6.0) - 3785760.0*Power(rij,7.0)*Power(xii,7.0) - 
        633600.0*Power(rij,8.0)*Power(xii,8.0) - 84480.0*Power(rij,9.0)*Power(xii,9.0) - 
        8448.0*Power(rij,10.0)*Power(xii,10.0) - 512.0*Power(rij,11.0)*Power(xii,11.0))/
      (6.386688e8*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (1260.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),11.0) + 
        210.0*Power(E,2.0*rij*xij)*Power(xij,10.0)*
         (-36.0*Power(rij,2.0)*Power(xii,14.0) - 2.0*Power(rij,3.0)*Power(xii,15.0) - 
           1287.0*rij*Power(xii,9.0)*Power(xij,4.0) + 6.0*Power(xij,12.0) + 
           9.0*rij*xii*Power(xij,12.0) - 
           22.0*rij*Power(xii,7.0)*Power(xij,6.0)*
            (-135.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*Power(xii,2.0)*Power(xij,10.0)*(-11.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           66.0*Power(xii,4.0)*Power(xij,8.0)*(-5.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           8.0*rij*Power(xii,5.0)*Power(xij,8.0)*(99.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           rij*Power(xii,3.0)*Power(xij,10.0)*(-99.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           132.0*Power(xii,6.0)*Power(xij,6.0)*(27.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           78.0*Power(xii,12.0)*(7.0 + 3.0*Power(rij,2.0)*Power(xij,2.0)) - 
           2.0*rij*Power(xii,13.0)*(117.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) + 
           66.0*Power(xii,8.0)*Power(xij,4.0)*(-191.0 + 6.0*Power(rij,2.0)*Power(xij,2.0)) + 
           rij*Power(xii,11.0)*Power(xij,2.0)*
            (-2151.0 + 22.0*Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*Power(xii,10.0)*Power(xij,2.0)*(-1099.0 + 33.0*Power(rij,2.0)*Power(xij,2.0))) \
    + Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (-385.0*Power(xii,8.0)*Power(xij,8.0)*
            (1080.0 + 1935.0*rij*xij + 1350.0*Power(rij,2.0)*Power(xij,2.0) + 
              1170.0*Power(rij,3.0)*Power(xij,3.0) + 420.0*Power(rij,4.0)*Power(xij,4.0) + 
              66.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           7.0*Power(xii,6.0)*Power(xij,10.0)*
            (99540.0 + 58095.0*rij*xij + 190710.0*Power(rij,2.0)*Power(xij,2.0) + 
              100950.0*Power(rij,3.0)*Power(xij,3.0) + 
              21660.0*Power(rij,4.0)*Power(xij,4.0) + 
              1938.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0) - 
              8.0*Power(rij,7.0)*Power(xij,7.0)) + 
           4.0*Power(xij,16.0)*(135135.0 + 135135.0*rij*xij + 
              62370.0*Power(rij,2.0)*Power(xij,2.0) + 
              17325.0*Power(rij,3.0)*Power(xij,3.0) + 
              3150.0*Power(rij,4.0)*Power(xij,4.0) + 378.0*Power(rij,5.0)*Power(xij,5.0) + 
              28.0*Power(rij,6.0)*Power(xij,6.0) + Power(rij,7.0)*Power(xij,7.0)) - 
           Power(xii,16.0)*(1260.0 + 2205.0*rij*xij + 
              1890.0*Power(rij,2.0)*Power(xij,2.0) + 1050.0*Power(rij,3.0)*Power(xij,3.0) + 
              420.0*Power(rij,4.0)*Power(xij,4.0) + 126.0*Power(rij,5.0)*Power(xij,5.0) + 
              28.0*Power(rij,6.0)*Power(xij,6.0) + 4.0*Power(rij,7.0)*Power(xij,7.0)) + 
           7.0*Power(xii,4.0)*Power(xij,12.0)*
            (114660.0 - 343395.0*rij*xij - 242910.0*Power(rij,2.0)*Power(xij,2.0) - 
              61950.0*Power(rij,3.0)*Power(xij,3.0) - 
              6060.0*Power(rij,4.0)*Power(xij,4.0) + 282.0*Power(rij,5.0)*Power(xij,5.0) + 
              116.0*Power(rij,6.0)*Power(xij,6.0) + 8.0*Power(rij,7.0)*Power(xij,7.0)) - 
           7.0*Power(xii,12.0)*Power(xij,4.0)*
            (9900.0 + 17325.0*rij*xij + 14850.0*Power(rij,2.0)*Power(xij,2.0) + 
              8250.0*Power(rij,3.0)*Power(xij,3.0) + 3300.0*Power(rij,4.0)*Power(xij,4.0) + 
              1074.0*Power(rij,5.0)*Power(xij,5.0) + 164.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0)) + 
           7.0*Power(xii,10.0)*Power(xij,6.0)*
            (29700.0 + 51975.0*rij*xij + 44550.0*Power(rij,2.0)*Power(xij,2.0) + 
              23850.0*Power(rij,3.0)*Power(xij,3.0) + 
              11700.0*Power(rij,4.0)*Power(xij,4.0) + 
              2814.0*Power(rij,5.0)*Power(xij,5.0) + 284.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0)) + 
           Power(xii,14.0)*Power(xij,2.0)*
            (13860.0 + 24255.0*rij*xij + 20790.0*Power(rij,2.0)*Power(xij,2.0) + 
              11550.0*Power(rij,3.0)*Power(xij,3.0) + 
              4620.0*Power(rij,4.0)*Power(xij,4.0) + 1386.0*Power(rij,5.0)*Power(xij,5.0) + 
              308.0*Power(rij,6.0)*Power(xij,6.0) + 24.0*Power(rij,7.0)*Power(xij,7.0)) - 
           Power(xii,2.0)*Power(xij,14.0)*
            (-3063060.0 - 1936935.0*rij*xij - 408870.0*Power(rij,2.0)*Power(xij,2.0) + 
              11550.0*Power(rij,3.0)*Power(xij,3.0) + 
              23100.0*Power(rij,4.0)*Power(xij,4.0) + 5082.0*Power(rij,5.0)*Power(xij,5.0) + 
              532.0*Power(rij,6.0)*Power(xij,6.0) + 24.0*Power(rij,7.0)*Power(xij,7.0))))/
      (1260.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),11.0))
    ;
  }
  return S;
}

static double Slater_2S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-124540416000.0 + 124540416000.0*Power(E,2.0*rij*xii) - 224622748350.0*rij*xii - 
        200164664700.0*Power(rij,2.0)*Power(xii,2.0) - 
        117249207075.0*Power(rij,3.0)*Power(xii,3.0) - 
        50639138550.0*Power(rij,4.0)*Power(xii,4.0) - 
        17132415300.0*Power(rij,5.0)*Power(xii,5.0) - 
        4704860160.0*Power(rij,6.0)*Power(xii,6.0) - 
        1071195840.0*Power(rij,7.0)*Power(xii,7.0) - 
        204478560.0*Power(rij,8.0)*Power(xii,8.0) - 
        32809920.0*Power(rij,9.0)*Power(xii,9.0) - 
        4392960.0*Power(rij,10.0)*Power(xii,10.0) - 
        479232.0*Power(rij,11.0)*Power(xii,11.0) - 39936.0*Power(rij,12.0)*Power(xii,12.0) - 
        2048.0*Power(rij,13.0)*Power(xii,13.0))/(1.24540416e11*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (28350.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),13.0) + 
        945.0*Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-210.0*Power(rij,2.0)*Power(xii,16.0) - 10.0*Power(rij,3.0)*Power(xii,17.0) + 
           30.0*Power(xij,14.0) + 45.0*rij*xii*Power(xij,14.0) + 
           39.0*rij*Power(xii,7.0)*Power(xij,8.0)*
            (1309.0 - 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           858.0*Power(xii,8.0)*Power(xij,6.0)*(-305.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           30.0*Power(xii,2.0)*Power(xij,12.0)*(-13.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           390.0*Power(xii,4.0)*Power(xij,10.0)*(-6.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           143.0*rij*Power(xii,9.0)*Power(xij,6.0)*
            (-153.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           5.0*rij*Power(xii,3.0)*Power(xij,12.0)*
            (-117.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           45.0*rij*Power(xii,15.0)*(35.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           138.0*Power(xii,12.0)*Power(xij,2.0)*
            (580.0 + 13.0*Power(rij,2.0)*Power(xij,2.0)) - 
           150.0*Power(xii,14.0)*(28.0 + 17.0*Power(rij,2.0)*Power(xij,2.0)) + 
           13.0*rij*Power(xii,11.0)*Power(xij,4.0)*
            (-4071.0 + 22.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*rij*Power(xii,13.0)*Power(xij,2.0)*
            (-8135.0 + 26.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*rij*Power(xii,5.0)*Power(xij,10.0)*
            (2171.0 + 30.0*Power(rij,2.0)*Power(xij,2.0)) + 
           234.0*Power(xii,10.0)*Power(xij,4.0)*
            (-1235.0 + 33.0*Power(rij,2.0)*Power(xij,2.0)) - 
           78.0*Power(xii,6.0)*Power(xij,8.0)*(550.0 + 47.0*Power(rij,2.0)*Power(xij,2.0))) - 
        2.0*Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (-819.0*Power(xii,10.0)*Power(xij,10.0)*
            (22275.0 + 39780.0*rij*xij + 38160.0*Power(rij,2.0)*Power(xij,2.0) + 
              16560.0*Power(rij,3.0)*Power(xij,3.0) + 
              9840.0*Power(rij,4.0)*Power(xij,4.0) + 3900.0*Power(rij,5.0)*Power(xij,5.0) + 
              816.0*Power(rij,6.0)*Power(xij,6.0) + 88.0*Power(rij,7.0)*Power(xij,7.0) + 
              4.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xii,20.0)*(14175.0 + 25515.0*rij*xij + 
              22680.0*Power(rij,2.0)*Power(xij,2.0) + 
              13230.0*Power(rij,3.0)*Power(xij,3.0) + 
              5670.0*Power(rij,4.0)*Power(xij,4.0) + 1890.0*Power(rij,5.0)*Power(xij,5.0) + 
              504.0*Power(rij,6.0)*Power(xij,6.0) + 108.0*Power(rij,7.0)*Power(xij,7.0) + 
              18.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) - 
           Power(xij,20.0)*(16216200.0 + 18243225.0*rij*xij + 
              9729720.0*Power(rij,2.0)*Power(xij,2.0) + 
              3243240.0*Power(rij,3.0)*Power(xij,3.0) + 
              748440.0*Power(rij,4.0)*Power(xij,4.0) + 
              124740.0*Power(rij,5.0)*Power(xij,5.0) + 
              15120.0*Power(rij,6.0)*Power(xij,6.0) + 
              1296.0*Power(rij,7.0)*Power(xij,7.0) + 72.0*Power(rij,8.0)*Power(xij,8.0) + 
              2.0*Power(rij,9.0)*Power(xij,9.0)) + 
           18.0*Power(xii,16.0)*Power(xij,4.0)*
            (61425.0 + 110565.0*rij*xij + 98280.0*Power(rij,2.0)*Power(xij,2.0) + 
              57330.0*Power(rij,3.0)*Power(xij,3.0) + 
              24570.0*Power(rij,4.0)*Power(xij,4.0) + 
              8190.0*Power(rij,5.0)*Power(xij,5.0) + 2184.0*Power(rij,6.0)*Power(xij,6.0) + 
              496.0*Power(rij,7.0)*Power(xij,7.0) + 64.0*Power(rij,8.0)*Power(xij,8.0) + 
              3.0*Power(rij,9.0)*Power(xij,9.0)) - 
           18.0*Power(xii,4.0)*Power(xij,16.0)*
            (6572475.0 - 3161340.0*rij*xij - 4782960.0*Power(rij,2.0)*Power(xij,2.0) - 
              1912365.0*Power(rij,3.0)*Power(xij,3.0) - 
              378105.0*Power(rij,4.0)*Power(xij,4.0) - 
              34125.0*Power(rij,5.0)*Power(xij,5.0) + 
              1092.0*Power(rij,6.0)*Power(xij,6.0) + 650.0*Power(rij,7.0)*Power(xij,7.0) + 
              71.0*Power(rij,8.0)*Power(xij,8.0) + 3.0*Power(rij,9.0)*Power(xij,9.0)) - 
           21.0*Power(xii,8.0)*Power(xij,12.0)*
            (-1063800.0 - 2775735.0*rij*xij - 862920.0*Power(rij,2.0)*Power(xij,2.0) - 
              1132020.0*Power(rij,3.0)*Power(xij,3.0) - 
              698580.0*Power(rij,4.0)*Power(xij,4.0) - 
              196920.0*Power(rij,5.0)*Power(xij,5.0) - 
              28992.0*Power(rij,6.0)*Power(xij,6.0) - 
              2064.0*Power(rij,7.0)*Power(xij,7.0) - 24.0*Power(rij,8.0)*Power(xij,8.0) + 
              4.0*Power(rij,9.0)*Power(xij,9.0)) + 
           21.0*Power(xii,12.0)*Power(xij,8.0)*
            (482625.0 + 868725.0*rij*xij + 772200.0*Power(rij,2.0)*Power(xij,2.0) + 
              455400.0*Power(rij,3.0)*Power(xij,3.0) + 
              178200.0*Power(rij,4.0)*Power(xij,4.0) + 
              72180.0*Power(rij,5.0)*Power(xij,5.0) + 
              19920.0*Power(rij,6.0)*Power(xij,6.0) + 
              2952.0*Power(rij,7.0)*Power(xij,7.0) + 204.0*Power(rij,8.0)*Power(xij,8.0) + 
              4.0*Power(rij,9.0)*Power(xij,9.0)) + 
           6.0*Power(xii,6.0)*Power(xij,14.0)*
            (-10357200.0 + 5071815.0*rij*xij - 6463800.0*Power(rij,2.0)*Power(xij,2.0) - 
              7151130.0*Power(rij,3.0)*Power(xij,3.0) - 
              2572290.0*Power(rij,4.0)*Power(xij,4.0) - 
              468720.0*Power(rij,5.0)*Power(xij,5.0) - 
              42672.0*Power(rij,6.0)*Power(xij,6.0) - 648.0*Power(rij,7.0)*Power(xij,7.0) + 
              228.0*Power(rij,8.0)*Power(xij,8.0) + 16.0*Power(rij,9.0)*Power(xij,9.0)) - 
           Power(xii,18.0)*Power(xij,2.0)*
            (184275.0 + 331695.0*rij*xij + 294840.0*Power(rij,2.0)*Power(xij,2.0) + 
              171990.0*Power(rij,3.0)*Power(xij,3.0) + 
              73710.0*Power(rij,4.0)*Power(xij,4.0) + 
              24570.0*Power(rij,5.0)*Power(xij,5.0) + 
              6552.0*Power(rij,6.0)*Power(xij,6.0) + 1404.0*Power(rij,7.0)*Power(xij,7.0) + 
              234.0*Power(rij,8.0)*Power(xij,8.0) + 16.0*Power(rij,9.0)*Power(xij,9.0)) + 
           Power(xii,2.0)*Power(xij,18.0)*
            (-133783650.0 - 107432325.0*rij*xij - 
              35675640.0*Power(rij,2.0)*Power(xij,2.0) - 
              5135130.0*Power(rij,3.0)*Power(xij,3.0) + 
              270270.0*Power(rij,4.0)*Power(xij,4.0) + 
              270270.0*Power(rij,5.0)*Power(xij,5.0) + 
              57960.0*Power(rij,6.0)*Power(xij,6.0) + 
              6948.0*Power(rij,7.0)*Power(xij,7.0) + 486.0*Power(rij,8.0)*Power(xij,8.0) + 
              16.0*Power(rij,9.0)*Power(xij,9.0)) - 
           6.0*Power(xii,14.0)*Power(xij,6.0)*
            (675675.0 + 1216215.0*rij*xij + 1081080.0*Power(rij,2.0)*Power(xij,2.0) + 
              630630.0*Power(rij,3.0)*Power(xij,3.0) + 
              270270.0*Power(rij,4.0)*Power(xij,4.0) + 
              88200.0*Power(rij,5.0)*Power(xij,5.0) + 
              26544.0*Power(rij,6.0)*Power(xij,6.0) + 5160.0*Power(rij,7.0)*Power(xij,7.0) + 
              492.0*Power(rij,8.0)*Power(xij,8.0) + 16.0*Power(rij,9.0)*Power(xij,9.0))))/
      (28350.*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),13.0))
    ;
  }
  return S;
}

static double Slater_2S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-125536739328000.0 + 125536739328000.0*Power(E,2.0*rij*xii) - 
        230286692010375.0*rij*xii - 209499905364750.0*Power(rij,2.0)*Power(xii,2.0) - 
        125847482260500.0*Power(rij,3.0)*Power(xii,3.0) - 
        56052916920000.0*Power(rij,4.0)*Power(xii,4.0) - 
        19698207328800.0*Power(rij,5.0)*Power(xii,5.0) - 
        5671583517600.0*Power(rij,6.0)*Power(xii,6.0) - 
        1370593224000.0*Power(rij,7.0)*Power(xii,7.0) - 
        282291609600.0*Power(rij,8.0)*Power(xii,8.0) - 
        49989139200.0*Power(rij,9.0)*Power(xii,9.0) - 
        7633866240.0*Power(rij,10.0)*Power(xii,10.0) - 
        1002193920.0*Power(rij,11.0)*Power(xii,11.0) - 
        111820800.0*Power(rij,12.0)*Power(xii,12.0) - 
        10321920.0*Power(rij,13.0)*Power(xii,13.0) - 
        737280.0*Power(rij,14.0)*Power(xii,14.0) - 32768.0*Power(rij,15.0)*Power(xii,15.0))/
      (1.25536739328e14*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (935550.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),15.0) + 
        51975.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-144.0*Power(rij,2.0)*Power(xii,18.0) - 6.0*Power(rij,3.0)*Power(xii,19.0) - 
           63999.0*rij*Power(xii,11.0)*Power(xij,6.0) + 18.0*Power(xij,16.0) + 
           27.0*rij*xii*Power(xij,16.0) + 
           18.0*Power(xii,2.0)*Power(xij,14.0)*(-15.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           270.0*Power(xii,4.0)*Power(xij,12.0)*(-7.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*rij*Power(xii,3.0)*Power(xij,14.0)*
            (-135.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           918.0*Power(xii,16.0)*(4.0 + 3.0*Power(rij,2.0)*Power(xij,2.0)) - 
           117.0*rij*Power(xii,9.0)*Power(xij,8.0)*
            (-1045.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) - 
           4.0*rij*Power(xii,17.0)*(306.0 + 23.0*Power(rij,2.0)*Power(xij,2.0)) - 
           3.0*rij*Power(xii,15.0)*Power(xij,2.0)*
            (9441.0 + 28.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*rij*Power(xii,7.0)*Power(xij,10.0)*
            (27261.0 + 28.0*Power(rij,2.0)*Power(xij,2.0)) + 
           9.0*rij*Power(xii,13.0)*Power(xij,4.0)*
            (-12915.0 + 52.0*Power(rij,2.0)*Power(xij,2.0)) + 
           234.0*Power(xii,10.0)*Power(xij,6.0)*
            (-4209.0 + 55.0*Power(rij,2.0)*Power(xij,2.0)) - 
           78.0*Power(xii,8.0)*Power(xij,8.0)*(6655.0 + 69.0*Power(rij,2.0)*Power(xij,2.0)) - 
           90.0*Power(xii,14.0)*Power(xij,2.0)*
            (1117.0 + 77.0*Power(rij,2.0)*Power(xij,2.0)) + 
           rij*Power(xii,5.0)*Power(xij,12.0)*
            (6111.0 + 92.0*Power(rij,2.0)*Power(xij,2.0)) - 
           18.0*Power(xii,6.0)*Power(xij,10.0)*
            (3107.0 + 259.0*Power(rij,2.0)*Power(xij,2.0)) + 
           18.0*Power(xii,12.0)*Power(xij,4.0)*
            (-31885.0 + 403.0*Power(rij,2.0)*Power(xij,2.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (-3465.0*Power(xii,12.0)*Power(xij,12.0)*
            (1351350.0 + 2483775.0*rij*xij + 2189250.0*Power(rij,2.0)*Power(xij,2.0) + 
              1499400.0*Power(rij,3.0)*Power(xij,3.0) + 
              512400.0*Power(rij,4.0)*Power(xij,4.0) + 
              191940.0*Power(rij,5.0)*Power(xij,5.0) + 
              73080.0*Power(rij,6.0)*Power(xij,6.0) + 
              18200.0*Power(rij,7.0)*Power(xij,7.0) + 
              2680.0*Power(rij,8.0)*Power(xij,8.0) + 220.0*Power(rij,9.0)*Power(xij,9.0) + 
              8.0*Power(rij,10.0)*Power(xij,10.0)) + 
           330.0*Power(xii,8.0)*Power(xij,16.0)*
            (-2409750.0 - 79762725.0*rij*xij - 9440550.0*Power(rij,2.0)*Power(xij,2.0) - 
              6036975.0*Power(rij,3.0)*Power(xij,3.0) - 
              10098900.0*Power(rij,4.0)*Power(xij,4.0) - 
              4800285.0*Power(rij,5.0)*Power(xij,5.0) - 
              1163190.0*Power(rij,6.0)*Power(xij,6.0) - 
              164670.0*Power(rij,7.0)*Power(xij,7.0) - 
              13110.0*Power(rij,8.0)*Power(xij,8.0) - 365.0*Power(rij,9.0)*Power(xij,9.0) + 
              26.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) + 
           2.0*Power(xij,24.0)*(1240539300.0 + 1516214700.0*rij*xij + 
              891891000.0*Power(rij,2.0)*Power(xij,2.0) + 
              334459125.0*Power(rij,3.0)*Power(xij,3.0) + 
              89189100.0*Power(rij,4.0)*Power(xij,4.0) + 
              17837820.0*Power(rij,5.0)*Power(xij,5.0) + 
              2744280.0*Power(rij,6.0)*Power(xij,6.0) + 
              326700.0*Power(rij,7.0)*Power(xij,7.0) + 
              29700.0*Power(rij,8.0)*Power(xij,8.0) + 
              1980.0*Power(rij,9.0)*Power(xij,9.0) + 88.0*Power(rij,10.0)*Power(xij,10.0) + 
              2.0*Power(rij,11.0)*Power(xij,11.0)) - 
           Power(xii,24.0)*(935550.0 + 1715175.0*rij*xij + 
              1559250.0*Power(rij,2.0)*Power(xij,2.0) + 
              935550.0*Power(rij,3.0)*Power(xij,3.0) + 
              415800.0*Power(rij,4.0)*Power(xij,4.0) + 
              145530.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              9900.0*Power(rij,7.0)*Power(xij,7.0) + 1980.0*Power(rij,8.0)*Power(xij,8.0) + 
              330.0*Power(rij,9.0)*Power(xij,9.0) + 44.0*Power(rij,10.0)*Power(xij,10.0) + 
              4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           110.0*Power(xii,6.0)*Power(xij,18.0)*
            (-313749450.0 + 140006475.0*rij*xij + 
              40682250.0*Power(rij,2.0)*Power(xij,2.0) - 
              63603225.0*Power(rij,3.0)*Power(xij,3.0) - 
              41107500.0*Power(rij,4.0)*Power(xij,4.0) - 
              11688705.0*Power(rij,5.0)*Power(xij,5.0) - 
              1918350.0*Power(rij,6.0)*Power(xij,6.0) - 
              179550.0*Power(rij,7.0)*Power(xij,7.0) - 
              5670.0*Power(rij,8.0)*Power(xij,8.0) + 735.0*Power(rij,9.0)*Power(xij,9.0) + 
              98.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           10.0*Power(xii,2.0)*Power(xij,22.0)*
            (-2825672850.0 - 2653375725.0*rij*xij - 
              1114863750.0*Power(rij,2.0)*Power(xij,2.0) - 
              260134875.0*Power(rij,3.0)*Power(xij,3.0) - 
              29729700.0*Power(rij,4.0)*Power(xij,4.0) + 
              1486485.0*Power(rij,5.0)*Power(xij,5.0) + 
              1295910.0*Power(rij,6.0)*Power(xij,6.0) + 
              272250.0*Power(rij,7.0)*Power(xij,7.0) + 
              34650.0*Power(rij,8.0)*Power(xij,8.0) + 
              2915.0*Power(rij,9.0)*Power(xij,9.0) + 
              154.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           165.0*Power(xii,16.0)*Power(xij,8.0)*
            (7739550.0 + 14189175.0*rij*xij + 12899250.0*Power(rij,2.0)*Power(xij,2.0) + 
              7739550.0*Power(rij,3.0)*Power(xij,3.0) + 
              3439800.0*Power(rij,4.0)*Power(xij,4.0) + 
              1210860.0*Power(rij,5.0)*Power(xij,5.0) + 
              330120.0*Power(rij,6.0)*Power(xij,6.0) + 
              86400.0*Power(rij,7.0)*Power(xij,7.0) + 
              18480.0*Power(rij,8.0)*Power(xij,8.0) + 
              2460.0*Power(rij,9.0)*Power(xij,9.0) + 
              168.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) + 
           5.0*Power(xii,22.0)*Power(xij,2.0)*
            (2806650.0 + 5145525.0*rij*xij + 4677750.0*Power(rij,2.0)*Power(xij,2.0) + 
              2806650.0*Power(rij,3.0)*Power(xij,3.0) + 
              1247400.0*Power(rij,4.0)*Power(xij,4.0) + 
              436590.0*Power(rij,5.0)*Power(xij,5.0) + 
              124740.0*Power(rij,6.0)*Power(xij,6.0) + 
              29700.0*Power(rij,7.0)*Power(xij,7.0) + 
              5940.0*Power(rij,8.0)*Power(xij,8.0) + 990.0*Power(rij,9.0)*Power(xij,9.0) + 
              132.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           55.0*Power(xii,18.0)*Power(xij,6.0)*
            (7739550.0 + 14189175.0*rij*xij + 12899250.0*Power(rij,2.0)*Power(xij,2.0) + 
              7739550.0*Power(rij,3.0)*Power(xij,3.0) + 
              3439800.0*Power(rij,4.0)*Power(xij,4.0) + 
              1203930.0*Power(rij,5.0)*Power(xij,5.0) + 
              343980.0*Power(rij,6.0)*Power(xij,6.0) + 
              80820.0*Power(rij,7.0)*Power(xij,7.0) + 
              17460.0*Power(rij,8.0)*Power(xij,8.0) + 
              2790.0*Power(rij,9.0)*Power(xij,9.0) + 
              244.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           22.0*Power(xii,4.0)*Power(xij,20.0)*
            (2199137850.0 + 366522975.0*rij*xij - 
              665232750.0*Power(rij,2.0)*Power(xij,2.0) - 
              422542575.0*Power(rij,3.0)*Power(xij,3.0) - 
              123095700.0*Power(rij,4.0)*Power(xij,4.0) - 
              20724795.0*Power(rij,5.0)*Power(xij,5.0) - 
              1838970.0*Power(rij,6.0)*Power(xij,6.0) + 
              12150.0*Power(rij,7.0)*Power(xij,7.0) + 
              26910.0*Power(rij,8.0)*Power(xij,8.0) + 
              3735.0*Power(rij,9.0)*Power(xij,9.0) + 
              258.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) - 
           33.0*Power(xii,10.0)*Power(xij,14.0)*
            (-188215650.0 - 280764225.0*rij*xij - 
              416886750.0*Power(rij,2.0)*Power(xij,2.0) - 
              131922000.0*Power(rij,3.0)*Power(xij,3.0) - 
              59043600.0*Power(rij,4.0)*Power(xij,4.0) - 
              34671420.0*Power(rij,5.0)*Power(xij,5.0) - 
              11740680.0*Power(rij,6.0)*Power(xij,6.0) - 
              2266200.0*Power(rij,7.0)*Power(xij,7.0) - 
              255000.0*Power(rij,8.0)*Power(xij,8.0) - 
              15060.0*Power(rij,9.0)*Power(xij,9.0) - 
              216.0*Power(rij,10.0)*Power(xij,10.0) + 16.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 11.0*Power(xii,20.0)*Power(xij,4.0)*
            (8930250.0 + 16372125.0*rij*xij + 14883750.0*Power(rij,2.0)*Power(xij,2.0) + 
              8930250.0*Power(rij,3.0)*Power(xij,3.0) + 
              3969000.0*Power(rij,4.0)*Power(xij,4.0) + 
              1389150.0*Power(rij,5.0)*Power(xij,5.0) + 
              396900.0*Power(rij,6.0)*Power(xij,6.0) + 
              94500.0*Power(rij,7.0)*Power(xij,7.0) + 
              18900.0*Power(rij,8.0)*Power(xij,8.0) + 
              3290.0*Power(rij,9.0)*Power(xij,9.0) + 
              364.0*Power(rij,10.0)*Power(xij,10.0) + 16.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 33.0*Power(xii,14.0)*Power(xij,10.0)*
            (85135050.0 + 156080925.0*rij*xij + 
              141891750.0*Power(rij,2.0)*Power(xij,2.0) + 
              84848400.0*Power(rij,3.0)*Power(xij,3.0) + 
              38984400.0*Power(rij,4.0)*Power(xij,4.0) + 
              12157740.0*Power(rij,5.0)*Power(xij,5.0) + 
              3814440.0*Power(rij,6.0)*Power(xij,6.0) + 
              1072200.0*Power(rij,7.0)*Power(xij,7.0) + 
              198120.0*Power(rij,8.0)*Power(xij,8.0) + 
              21020.0*Power(rij,9.0)*Power(xij,9.0) + 
              1096.0*Power(rij,10.0)*Power(xij,10.0) + 16.0*Power(rij,11.0)*Power(xij,11.0)))\
    )/(935550.*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),15.0))
    ;
  }
  return S;
}

double Slater_2S_1S(double rij,double xii,double xij)
{
  return Slater_1S_2S(rij,xij,xii);
}

static double Slater_3S_3S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-1437004800.0 + 1437004800.0*Power(E,2.0*rij*xii) - 2503064025.0*rij*xii - 
        2132118450.0*Power(rij,2.0)*Power(xii,2.0) - 
        1180664100.0*Power(rij,3.0)*Power(xii,3.0) - 
        476506800.0*Power(rij,4.0)*Power(xii,4.0) - 
        148856400.0*Power(rij,5.0)*Power(xii,5.0) - 
        37255680.0*Power(rij,6.0)*Power(xii,6.0) - 7603200.0*Power(rij,7.0)*Power(xii,7.0) - 
        1267200.0*Power(rij,8.0)*Power(xii,8.0) - 168960.0*Power(rij,9.0)*Power(xii,9.0) - 
        16896.0*Power(rij,10.0)*Power(xii,10.0) - 1024.0*Power(rij,11.0)*Power(xii,11.0))/
      (1.4370048e9*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (135.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),11.0) + 
        Power(E,2.0*rij*xij)*Power(xij,8.0)*
         (-150.0*Power(rij,4.0)*Power(xii,18.0) - 6.0*Power(rij,5.0)*Power(xii,19.0) + 
           135.0*Power(xij,14.0) + 225.0*rij*xii*Power(xij,14.0) + 
           10.0*Power(rij,3.0)*Power(xii,17.0)*(-165.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           30.0*Power(rij,2.0)*Power(xii,16.0)*(330.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           45.0*rij*Power(xii,3.0)*Power(xij,12.0)*
            (-55.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           45.0*Power(xii,2.0)*Power(xij,12.0)*(-33.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) + 
           rij*Power(xii,9.0)*Power(xij,6.0)*
            (234135.0 - 4950.0*Power(rij,2.0)*Power(xij,2.0) - 
              34.0*Power(rij,4.0)*Power(xij,4.0)) - 
           5.0*rij*Power(xii,7.0)*Power(xij,8.0)*
            (6237.0 - 1242.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           3.0*rij*Power(xii,5.0)*Power(xij,10.0)*
            (4125.0 - 330.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           15.0*Power(xii,4.0)*Power(xij,10.0)*
            (495.0 - 132.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    - 165.0*Power(xii,6.0)*Power(xij,8.0)*
            (135.0 - 60.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    - 5.0*rij*Power(xii,13.0)*Power(xij,2.0)*
            (43875.0 - 3438.0*Power(rij,2.0)*Power(xij,2.0) + 
              22.0*Power(rij,4.0)*Power(xij,4.0)) + 
           5.0*rij*Power(xii,11.0)*Power(xij,4.0)*
            (7695.0 - 2442.0*Power(rij,2.0)*Power(xij,2.0) + 
              22.0*Power(rij,4.0)*Power(xij,4.0)) + 
           15.0*Power(xii,8.0)*Power(xij,6.0)*
            (-33.0 - 3564.0*Power(rij,2.0)*Power(xij,2.0) + 
              26.0*Power(rij,4.0)*Power(xij,4.0)) + 
           rij*Power(xii,15.0)*(-32175.0 - 3690.0*Power(rij,2.0)*Power(xij,2.0) + 
              34.0*Power(rij,4.0)*Power(xij,4.0)) + 
           15.0*Power(xii,10.0)*Power(xij,4.0)*
            (-32277.0 + 1364.0*Power(rij,2.0)*Power(xij,2.0) + 
              66.0*Power(rij,4.0)*Power(xij,4.0)) + 
           15.0*Power(xii,14.0)*(-3003.0 - 2932.0*Power(rij,2.0)*Power(xij,2.0) + 
              94.0*Power(rij,4.0)*Power(xij,4.0)) - 
           15.0*Power(xii,12.0)*Power(xij,2.0)*
            (28119.0 - 5252.0*Power(rij,2.0)*Power(xij,2.0) + 
              154.0*Power(rij,4.0)*Power(xij,4.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,8.0)*
         (-5.0*Power(xii,2.0)*Power(xij,12.0)*
            (-84357.0 - 43875.0*rij*xij - 8796.0*Power(rij,2.0)*Power(xij,2.0) - 
              738.0*Power(rij,3.0)*Power(xij,3.0) - 6.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) - 
           3.0*Power(xii,14.0)*(45.0 + 75.0*rij*xij + 60.0*Power(rij,2.0)*Power(xij,2.0) + 
              30.0*Power(rij,3.0)*Power(xij,3.0) + 10.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) - 
           55.0*Power(xii,8.0)*Power(xij,6.0)*
            (-405.0 - 567.0*rij*xij - 972.0*Power(rij,2.0)*Power(xij,2.0) - 
              90.0*Power(rij,3.0)*Power(xij,3.0) + 18.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) + 
           55.0*Power(xii,6.0)*Power(xij,8.0)*
            (9.0 - 4257.0*rij*xij - 372.0*Power(rij,2.0)*Power(xij,2.0) + 
              222.0*Power(rij,3.0)*Power(xij,3.0) + 42.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) + 
           3.0*Power(xij,14.0)*(15015.0 + 10725.0*rij*xij + 
              3300.0*Power(rij,2.0)*Power(xij,2.0) + 550.0*Power(rij,3.0)*Power(xij,3.0) + 
              50.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,5.0)*Power(xij,5.0)) + 
           5.0*Power(xii,12.0)*Power(xij,2.0)*
            (297.0 + 495.0*rij*xij + 396.0*Power(rij,2.0)*Power(xij,2.0) + 
              198.0*Power(rij,3.0)*Power(xij,3.0) + 66.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) + 
           Power(xii,10.0)*Power(xij,4.0)*
            (-7425.0 - 12375.0*rij*xij - 9900.0*Power(rij,2.0)*Power(xij,2.0) - 
              6210.0*Power(rij,3.0)*Power(xij,3.0) - 390.0*Power(rij,4.0)*Power(xij,4.0) + 
              34.0*Power(rij,5.0)*Power(xij,5.0)) - 
           Power(xii,4.0)*Power(xij,10.0)*
            (-484155.0 + 38475.0*rij*xij + 78780.0*Power(rij,2.0)*Power(xij,2.0) + 
              17190.0*Power(rij,3.0)*Power(xij,3.0) + 1410.0*Power(rij,4.0)*Power(xij,4.0) + 
              34.0*Power(rij,5.0)*Power(xij,5.0))))/
      (135.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),11.0))
    ;
  }
  return S;
}

static double Slater_3S_4S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-74724249600.0 + 74724249600.0*Power(E,2.0*rij*xii) - 132871488750.0*rij*xii - 
        116294478300.0*Power(rij,2.0)*Power(xii,2.0) - 
        66678987375.0*Power(rij,3.0)*Power(xii,3.0) - 
        28114836750.0*Power(rij,4.0)*Power(xii,4.0) - 
        9274044780.0*Power(rij,5.0)*Power(xii,5.0) - 
        2484321840.0*Power(rij,6.0)*Power(xii,6.0) - 
        553204080.0*Power(rij,7.0)*Power(xii,7.0) - 
        103783680.0*Power(rij,8.0)*Power(xii,8.0) - 
        16473600.0*Power(rij,9.0)*Power(xii,9.0) - 
        2196480.0*Power(rij,10.0)*Power(xii,10.0) - 
        239616.0*Power(rij,11.0)*Power(xii,11.0) - 19968.0*Power(rij,12.0)*Power(xii,12.0) - 
        1024.0*Power(rij,13.0)*Power(xii,13.0))/(7.47242496e10*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (3780.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),13.0) + 
        84.0*Power(E,2.0*rij*xij)*Power(xij,10.0)*
         (-60.0*Power(rij,4.0)*Power(xii,20.0) - 2.0*Power(rij,5.0)*Power(xii,21.0) + 
           45.0*Power(xij,16.0) + 75.0*rij*xii*Power(xij,16.0) - 
           4.0*Power(rij,3.0)*Power(xii,19.0)*(195.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           15.0*rij*Power(xii,3.0)*Power(xij,14.0)*
            (-65.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           15.0*Power(xii,2.0)*Power(xij,14.0)*(-39.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) - 
           30.0*Power(rij,2.0)*Power(xii,18.0)*(182.0 + 9.0*Power(rij,2.0)*Power(xij,2.0)) + 
           30.0*rij*Power(xii,13.0)*Power(xij,4.0)*
            (-13047.0 + 377.0*Power(rij,2.0)*Power(xij,2.0)) + 
           2.0*rij*Power(xii,5.0)*Power(xij,12.0)*
            (2925.0 - 195.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) \
    + 10.0*Power(xii,4.0)*Power(xij,12.0)*
            (351.0 - 78.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) - 
           130.0*Power(xii,6.0)*Power(xij,10.0)*
            (99.0 - 36.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           13.0*rij*Power(xii,11.0)*Power(xij,6.0)*
            (30735.0 - 1650.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,4.0)*Power(xij,4.0)) + 
           rij*Power(xii,7.0)*Power(xij,10.0)*
            (-15015.0 + 3330.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,4.0)*Power(xij,4.0)) + 
           210.0*Power(xii,16.0)*(-156.0 - 262.0*Power(rij,2.0)*Power(xij,2.0) + 
              5.0*Power(rij,4.0)*Power(xij,4.0)) - 
           6.0*rij*Power(xii,9.0)*Power(xij,8.0)*
            (-48620.0 - 715.0*Power(rij,2.0)*Power(xij,2.0) + 
              6.0*Power(rij,4.0)*Power(xij,4.0)) + 
           3.0*rij*Power(xii,17.0)*(-6825.0 - 1870.0*Power(rij,2.0)*Power(xij,2.0) + 
              12.0*Power(rij,4.0)*Power(xij,4.0)) - 
           30.0*Power(xii,14.0)*Power(xij,2.0)*
            (17934.0 - 12.0*Power(rij,2.0)*Power(xij,2.0) + 
              13.0*Power(rij,4.0)*Power(xij,4.0)) - 
           15.0*Power(xii,8.0)*Power(xij,8.0)*
            (2145.0 + 2860.0*Power(rij,2.0)*Power(xij,2.0) + 
              14.0*Power(rij,4.0)*Power(xij,4.0)) + 
           65.0*Power(xii,10.0)*Power(xij,6.0)*
            (-13725.0 - 792.0*Power(rij,2.0)*Power(xij,2.0) + 
              22.0*Power(rij,4.0)*Power(xij,4.0)) - 
           10.0*Power(xii,12.0)*Power(xij,4.0)*
            (153630.0 - 15054.0*Power(rij,2.0)*Power(xij,2.0) + 
              143.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,15.0)*(-269325.0*rij*Power(xij,2.0) + 
              9270.0*Power(rij,3.0)*Power(xij,4.0) - 52.0*Power(rij,5.0)*Power(xij,6.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,8.0)*
         (Power(xii,2.0)*Power(xij,16.0)*
            (70073640.0 + 47669895.0*rij*xij + 13931190.0*Power(rij,2.0)*Power(xij,2.0) + 
              2170350.0*Power(rij,3.0)*Power(xij,3.0) + 
              169260.0*Power(rij,4.0)*Power(xij,4.0) + 
              1638.0*Power(rij,5.0)*Power(xij,5.0) - 756.0*Power(rij,6.0)*Power(xij,6.0) - 
              44.0*Power(rij,7.0)*Power(xij,7.0)) + 
           364.0*Power(xii,10.0)*Power(xij,8.0)*
            (-7425.0 - 13860.0*rij*xij - 5940.0*Power(rij,2.0)*Power(xij,2.0) - 
              11880.0*Power(rij,3.0)*Power(xij,3.0) - 
              2640.0*Power(rij,4.0)*Power(xij,4.0) - 45.0*Power(rij,5.0)*Power(xij,5.0) + 
              30.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,7.0)*Power(xij,7.0)) - 
           364.0*Power(xii,8.0)*Power(xij,10.0)*
            (-20925.0 + 18270.0*rij*xij - 58320.0*Power(rij,2.0)*Power(xij,2.0) - 
              17730.0*Power(rij,3.0)*Power(xij,3.0) - 300.0*Power(rij,4.0)*Power(xij,4.0) + 
              423.0*Power(rij,5.0)*Power(xij,5.0) + 54.0*Power(rij,6.0)*Power(xij,6.0) + 
              2.0*Power(rij,7.0)*Power(xij,7.0)) - 
           3.0*Power(xii,18.0)*(1260.0 + 2205.0*rij*xij + 
              1890.0*Power(rij,2.0)*Power(xij,2.0) + 1050.0*Power(rij,3.0)*Power(xij,3.0) + 
              420.0*Power(rij,4.0)*Power(xij,4.0) + 126.0*Power(rij,5.0)*Power(xij,5.0) + 
              28.0*Power(rij,6.0)*Power(xij,6.0) + 4.0*Power(rij,7.0)*Power(xij,7.0)) + 
           3.0*Power(xij,18.0)*(1801800.0 + 1576575.0*rij*xij + 
              630630.0*Power(rij,2.0)*Power(xij,2.0) + 
              150150.0*Power(rij,3.0)*Power(xij,3.0) + 
              23100.0*Power(rij,4.0)*Power(xij,4.0) + 
              2310.0*Power(rij,5.0)*Power(xij,5.0) + 140.0*Power(rij,6.0)*Power(xij,6.0) + 
              4.0*Power(rij,7.0)*Power(xij,7.0)) + 
           2.0*Power(xii,14.0)*Power(xij,4.0)*
            (-147420.0 - 257985.0*rij*xij - 221130.0*Power(rij,2.0)*Power(xij,2.0) - 
              122850.0*Power(rij,3.0)*Power(xij,3.0) - 
              49140.0*Power(rij,4.0)*Power(xij,4.0) - 
              17388.0*Power(rij,5.0)*Power(xij,5.0) - 
              1512.0*Power(rij,6.0)*Power(xij,6.0) + 8.0*Power(rij,7.0)*Power(xij,7.0)) - 
           42.0*Power(xii,12.0)*Power(xij,6.0)*
            (-25740.0 - 45045.0*rij*xij - 38610.0*Power(rij,2.0)*Power(xij,2.0) - 
              19470.0*Power(rij,3.0)*Power(xij,3.0) - 
              12540.0*Power(rij,4.0)*Power(xij,4.0) - 
              1836.0*Power(rij,5.0)*Power(xij,5.0) - 8.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0)) + 
           42.0*Power(xii,6.0)*Power(xij,12.0)*
            (921600.0 - 1640835.0*rij*xij - 546030.0*Power(rij,2.0)*Power(xij,2.0) + 
              20730.0*Power(rij,3.0)*Power(xij,3.0) + 
              30180.0*Power(rij,4.0)*Power(xij,4.0) + 
              5028.0*Power(rij,5.0)*Power(xij,5.0) + 344.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0)) - 
           2.0*Power(xii,4.0)*Power(xij,14.0)*
            (-67767840.0 - 13377735.0*rij*xij + 6601770.0*Power(rij,2.0)*Power(xij,2.0) + 
              3115350.0*Power(rij,3.0)*Power(xij,3.0) + 
              548940.0*Power(rij,4.0)*Power(xij,4.0) + 
              48132.0*Power(rij,5.0)*Power(xij,5.0) + 
              1848.0*Power(rij,6.0)*Power(xij,6.0) + 8.0*Power(rij,7.0)*Power(xij,7.0)) + 
           Power(xii,16.0)*Power(xij,2.0)*
            (49140.0 + 85995.0*rij*xij + 73710.0*Power(rij,2.0)*Power(xij,2.0) + 
              40950.0*Power(rij,3.0)*Power(xij,3.0) + 
              16380.0*Power(rij,4.0)*Power(xij,4.0) + 4914.0*Power(rij,5.0)*Power(xij,5.0) + 
              1092.0*Power(rij,6.0)*Power(xij,6.0) + 44.0*Power(rij,7.0)*Power(xij,7.0))))/
      (3780.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),13.0))
    ;
  }
  return S;
}

static double Slater_3S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-313841848320000.0 + 313841848320000.0*Power(E,2.0*rij*xii) - 
        568188982486125.0*rij*xii - 508694268332250.0*Power(rij,2.0)*Power(xii,2.0) - 
        299892470377500.0*Power(rij,3.0)*Power(xii,3.0) - 
        130753815192000.0*Power(rij,4.0)*Power(xii,4.0) - 
        44881155118800.0*Power(rij,5.0)*Power(xii,5.0) - 
        12601803614400.0*Power(rij,6.0)*Power(xii,6.0) - 
        2967953788800.0*Power(rij,7.0)*Power(xii,7.0) - 
        596237241600.0*Power(rij,8.0)*Power(xii,8.0) - 
        103264761600.0*Power(rij,9.0)*Power(xii,9.0) - 
        15498362880.0*Power(rij,10.0)*Power(xii,10.0) - 
        2012774400.0*Power(rij,11.0)*Power(xii,11.0) - 
        223641600.0*Power(rij,12.0)*Power(xii,12.0) - 
        20643840.0*Power(rij,13.0)*Power(xii,13.0) - 
        1474560.0*Power(rij,14.0)*Power(xii,14.0) - 65536.0*Power(rij,15.0)*Power(xii,15.0))/
      (3.1384184832e14*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (42525.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),15.0) + 
        189.0*Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-350.0*Power(rij,4.0)*Power(xii,22.0) - 10.0*Power(rij,5.0)*Power(xii,23.0) + 
           225.0*Power(xij,18.0) + 375.0*rij*xii*Power(xij,18.0) - 
           70.0*Power(rij,3.0)*Power(xii,21.0)*(75.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           75.0*rij*Power(xii,3.0)*Power(xij,16.0)*
            (-75.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           75.0*Power(xii,2.0)*Power(xij,16.0)*(-45.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) - 
           50.0*Power(rij,2.0)*Power(xii,20.0)*(840.0 + 71.0*Power(rij,2.0)*Power(xij,2.0)) + 
           rij*Power(xii,9.0)*Power(xij,10.0)*
            (4694625.0 + 124800.0*Power(rij,2.0)*Power(xij,2.0) - 
              248.0*Power(rij,4.0)*Power(xij,4.0)) + 
           20.0*rij*Power(xii,17.0)*Power(xij,2.0)*
            (-185895.0 - 948.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           5.0*rij*Power(xii,5.0)*Power(xij,14.0)*
            (7875.0 - 450.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           25.0*Power(xii,4.0)*Power(xij,14.0)*
            (945.0 - 180.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    - 375.0*Power(xii,6.0)*Power(xij,12.0)*
            (273.0 - 84.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    - 5.0*rij*Power(xii,11.0)*Power(xij,8.0)*
            (-2803125.0 + 49140.0*Power(rij,2.0)*Power(xij,2.0) + 
              8.0*Power(rij,4.0)*Power(xij,4.0)) + 
           5.0*rij*Power(xii,7.0)*Power(xij,12.0)*
            (-16965.0 + 5152.0*Power(rij,2.0)*Power(xij,2.0) + 
              14.0*Power(rij,4.0)*Power(xij,4.0)) + 
           325.0*Power(xii,10.0)*Power(xij,8.0)*
            (-60117.0 - 5340.0*Power(rij,2.0)*Power(xij,2.0) + 
              40.0*Power(rij,4.0)*Power(xij,4.0)) - 
           15.0*rij*Power(xii,15.0)*Power(xij,4.0)*
            (845085.0 - 22960.0*Power(rij,2.0)*Power(xij,2.0) + 
              52.0*Power(rij,4.0)*Power(xij,4.0)) + 
           15.0*rij*Power(xii,13.0)*Power(xij,6.0)*
            (-139125.0 - 10140.0*Power(rij,2.0)*Power(xij,2.0) + 
              52.0*Power(rij,4.0)*Power(xij,4.0)) + 
           75.0*Power(xii,12.0)*Power(xij,6.0)*
            (-729687.0 + 25532.0*Power(rij,2.0)*Power(xij,2.0) + 
              52.0*Power(rij,4.0)*Power(xij,4.0)) + 
           60.0*Power(xii,18.0)*(-5355.0 - 11940.0*Power(rij,2.0)*Power(xij,2.0) + 
              86.0*Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*rij*Power(xii,19.0)*(-89250.0 - 35425.0*Power(rij,2.0)*Power(xij,2.0) + 
              124.0*Power(rij,4.0)*Power(xij,4.0)) + 
           100.0*Power(xii,16.0)*Power(xij,2.0)*
            (-79713.0 - 13311.0*Power(rij,2.0)*Power(xij,2.0) + 
              146.0*Power(rij,4.0)*Power(xij,4.0)) - 
           5.0*Power(xii,8.0)*Power(xij,10.0)*
            (157365.0 + 95940.0*Power(rij,2.0)*Power(xij,2.0) + 
              952.0*Power(rij,4.0)*Power(xij,4.0)) - 
           15.0*Power(xii,14.0)*Power(xij,4.0)*
            (2638467.0 - 157500.0*Power(rij,2.0)*Power(xij,2.0) + 
              1820.0*Power(rij,4.0)*Power(xij,4.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,8.0)*
         (2.0*Power(xii,2.0)*Power(xij,20.0)*
            (1782492075.0 + 1449175455.0*rij*xij + 
              533365560.0*Power(rij,2.0)*Power(xij,2.0) + 
              114631335.0*Power(rij,3.0)*Power(xij,3.0) + 
              15221115.0*Power(rij,4.0)*Power(xij,4.0) + 
              1142505.0*Power(rij,5.0)*Power(xij,5.0) + 
              18396.0*Power(rij,6.0)*Power(xij,6.0) - 
              5238.0*Power(rij,7.0)*Power(xij,7.0) - 513.0*Power(rij,8.0)*Power(xij,8.0) - 
              17.0*Power(rij,9.0)*Power(xij,9.0)) + 
           42.0*Power(xii,4.0)*Power(xij,18.0)*
            (251336925.0 + 104824125.0*rij*xij + 340200.0*Power(rij,2.0)*Power(xij,2.0) - 
              9122085.0*Power(rij,3.0)*Power(xij,3.0) - 
              2798145.0*Power(rij,4.0)*Power(xij,4.0) - 
              433755.0*Power(rij,5.0)*Power(xij,5.0) - 
              39060.0*Power(rij,6.0)*Power(xij,6.0) - 
              1890.0*Power(rij,7.0)*Power(xij,7.0) - 27.0*Power(rij,8.0)*Power(xij,8.0) + 
              Power(rij,9.0)*Power(xij,9.0)) + 
           6.0*Power(xij,22.0)*(34459425.0 + 34459425.0*rij*xij + 
              16216200.0*Power(rij,2.0)*Power(xij,2.0) + 
              4729725.0*Power(rij,3.0)*Power(xij,3.0) + 
              945945.0*Power(rij,4.0)*Power(xij,4.0) + 
              135135.0*Power(rij,5.0)*Power(xij,5.0) + 
              13860.0*Power(rij,6.0)*Power(xij,6.0) + 990.0*Power(rij,7.0)*Power(xij,7.0) + 
              45.0*Power(rij,8.0)*Power(xij,8.0) + Power(rij,9.0)*Power(xij,9.0)) - 
           3.0*Power(xii,22.0)*(14175.0 + 25515.0*rij*xij + 
              22680.0*Power(rij,2.0)*Power(xij,2.0) + 
              13230.0*Power(rij,3.0)*Power(xij,3.0) + 
              5670.0*Power(rij,4.0)*Power(xij,4.0) + 1890.0*Power(rij,5.0)*Power(xij,5.0) + 
              504.0*Power(rij,6.0)*Power(xij,6.0) + 108.0*Power(rij,7.0)*Power(xij,7.0) + 
              18.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) - 
           21.0*Power(xii,18.0)*Power(xij,4.0)*
            (212625.0 + 382725.0*rij*xij + 340200.0*Power(rij,2.0)*Power(xij,2.0) + 
              198450.0*Power(rij,3.0)*Power(xij,3.0) + 
              85050.0*Power(rij,4.0)*Power(xij,4.0) + 
              28350.0*Power(rij,5.0)*Power(xij,5.0) + 
              7560.0*Power(rij,6.0)*Power(xij,6.0) + 1836.0*Power(rij,7.0)*Power(xij,7.0) + 
              162.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) + 
           54.0*Power(xii,6.0)*Power(xij,16.0)*
            (133451955.0 - 73700865.0*rij*xij - 
              54096840.0*Power(rij,2.0)*Power(xij,2.0) - 
              8306235.0*Power(rij,3.0)*Power(xij,3.0) + 
              966945.0*Power(rij,4.0)*Power(xij,4.0) + 
              516747.0*Power(rij,5.0)*Power(xij,5.0) + 
              80724.0*Power(rij,6.0)*Power(xij,6.0) + 
              6434.0*Power(rij,7.0)*Power(xij,7.0) + 251.0*Power(rij,8.0)*Power(xij,8.0) + 
              3.0*Power(rij,9.0)*Power(xij,9.0)) - 
           315.0*Power(xii,12.0)*Power(xij,10.0)*
            (-405405.0 - 710073.0*rij*xij - 805896.0*Power(rij,2.0)*Power(xij,2.0) - 
              101556.0*Power(rij,3.0)*Power(xij,3.0) - 
              258804.0*Power(rij,4.0)*Power(xij,4.0) - 
              90972.0*Power(rij,5.0)*Power(xij,5.0) - 
              9744.0*Power(rij,6.0)*Power(xij,6.0) + 120.0*Power(rij,7.0)*Power(xij,7.0) + 
              84.0*Power(rij,8.0)*Power(xij,8.0) + 4.0*Power(rij,9.0)*Power(xij,9.0)) + 
           315.0*Power(xii,10.0)*Power(xij,12.0)*
            (-482895.0 - 2656395.0*rij*xij + 1186920.0*Power(rij,2.0)*Power(xij,2.0) - 
              1155420.0*Power(rij,3.0)*Power(xij,3.0) - 
              643356.0*Power(rij,4.0)*Power(xij,4.0) - 
              93492.0*Power(rij,5.0)*Power(xij,5.0) + 336.0*Power(rij,6.0)*Power(xij,6.0) + 
              1368.0*Power(rij,7.0)*Power(xij,7.0) + 132.0*Power(rij,8.0)*Power(xij,8.0) + 
              4.0*Power(rij,9.0)*Power(xij,9.0)) - 
           27.0*Power(xii,16.0)*Power(xij,6.0)*
            (-716625.0 - 1289925.0*rij*xij - 1146600.0*Power(rij,2.0)*Power(xij,2.0) - 
              668850.0*Power(rij,3.0)*Power(xij,3.0) - 
              286650.0*Power(rij,4.0)*Power(xij,4.0) - 
              90006.0*Power(rij,5.0)*Power(xij,5.0) - 
              32872.0*Power(rij,6.0)*Power(xij,6.0) - 
              4812.0*Power(rij,7.0)*Power(xij,7.0) - 178.0*Power(rij,8.0)*Power(xij,8.0) + 
              6.0*Power(rij,9.0)*Power(xij,9.0)) + 
           Power(xii,20.0)*Power(xij,2.0)*
            (637875.0 + 1148175.0*rij*xij + 1020600.0*Power(rij,2.0)*Power(xij,2.0) + 
              595350.0*Power(rij,3.0)*Power(xij,3.0) + 
              255150.0*Power(rij,4.0)*Power(xij,4.0) + 
              85050.0*Power(rij,5.0)*Power(xij,5.0) + 
              22680.0*Power(rij,6.0)*Power(xij,6.0) + 
              4860.0*Power(rij,7.0)*Power(xij,7.0) + 810.0*Power(rij,8.0)*Power(xij,8.0) + 
              34.0*Power(rij,9.0)*Power(xij,9.0)) + 
           3.0*Power(xii,14.0)*Power(xij,8.0)*
            (-19348875.0 - 34827975.0*rij*xij - 
              30958200.0*Power(rij,2.0)*Power(xij,2.0) - 
              18689580.0*Power(rij,3.0)*Power(xij,3.0) - 
              5847660.0*Power(rij,4.0)*Power(xij,4.0) - 
              3723300.0*Power(rij,5.0)*Power(xij,5.0) - 
              845040.0*Power(rij,6.0)*Power(xij,6.0) - 
              58680.0*Power(rij,7.0)*Power(xij,7.0) + 
              1548.0*Power(rij,8.0)*Power(xij,8.0) + 236.0*Power(rij,9.0)*Power(xij,9.0)) - 
           3.0*Power(xii,8.0)*Power(xij,14.0)*
            (-593408025.0 + 946053675.0*rij*xij - 
              394427880.0*Power(rij,2.0)*Power(xij,2.0) - 
              315870660.0*Power(rij,3.0)*Power(xij,3.0) - 
              53891460.0*Power(rij,4.0)*Power(xij,4.0) + 
              910980.0*Power(rij,5.0)*Power(xij,5.0) + 
              1409520.0*Power(rij,6.0)*Power(xij,6.0) + 
              192168.0*Power(rij,7.0)*Power(xij,7.0) + 
              11196.0*Power(rij,8.0)*Power(xij,8.0) + 236.0*Power(rij,9.0)*Power(xij,9.0))))/
      (42525.*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),15.0))
    ;
  }
  return S;
}

static double Slater_3S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-12804747411456000.0 + 12804747411456000.0*Power(E,2.0*rij*xii) - 
        23523793155237375.0*rij*xii - 
        21438091487562750.0*Power(rij,2.0)*Power(xii,2.0) - 
        12909495448599750.0*Power(rij,3.0)*Power(xii,3.0) - 
        5771367188086500.0*Power(rij,4.0)*Power(xii,4.0) - 
        2040067705876200.0*Power(rij,5.0)*Power(xii,5.0) - 
        592812380160000.0*Power(rij,6.0)*Power(xii,6.0) - 
        145331660073600.0*Power(rij,7.0)*Power(xii,7.0) - 
        30604380206400.0*Power(rij,8.0)*Power(xii,8.0) - 
        5606134934400.0*Power(rij,9.0)*Power(xii,9.0) - 
        900980720640.0*Power(rij,10.0)*Power(xii,10.0) - 
        127672796160.0*Power(rij,11.0)*Power(xii,11.0) - 
        15968010240.0*Power(rij,12.0)*Power(xii,12.0) - 
        1754726400.0*Power(rij,13.0)*Power(xii,13.0) - 
        167116800.0*Power(rij,14.0)*Power(xii,14.0) - 
        13369344.0*Power(rij,15.0)*Power(xii,15.0) - 
        835584.0*Power(rij,16.0)*Power(xii,16.0) - 32768.0*Power(rij,17.0)*Power(xii,17.0))/
      (1.2804747411456e16*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (2806650.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),17.0) + 
        20790.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-240.0*Power(rij,4.0)*Power(xii,24.0) - 6.0*Power(rij,5.0)*Power(xii,25.0) + 
           135.0*Power(xij,20.0) + 225.0*rij*xii*Power(xij,20.0) - 
           80.0*Power(rij,3.0)*Power(xii,23.0)*(51.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           45.0*rij*Power(xii,3.0)*Power(xij,18.0)*
            (-85.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           45.0*Power(xii,2.0)*Power(xij,18.0)*(-51.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) - 
           30.0*Power(rij,2.0)*Power(xii,22.0)*
            (1224.0 + 137.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3060.0*rij*Power(xii,15.0)*Power(xij,6.0)*
            (-11875.0 + 146.0*Power(rij,2.0)*Power(xij,2.0)) + 
           2.0*rij*Power(xii,9.0)*Power(xij,12.0)*
            (3977235.0 + 115260.0*Power(rij,2.0)*Power(xij,2.0) - 
              47.0*Power(rij,4.0)*Power(xij,4.0)) + 
           1020.0*rij*Power(xii,13.0)*Power(xij,8.0)*
            (23775.0 - 741.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) \
    + 6.0*rij*Power(xii,5.0)*Power(xij,16.0)*
            (5100.0 - 255.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) \
    + 30.0*Power(xii,4.0)*Power(xij,16.0)*
            (612.0 - 102.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) - 
           510.0*Power(xii,6.0)*Power(xij,14.0)*
            (180.0 - 48.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           20.0*rij*Power(xii,7.0)*Power(xij,14.0)*
            (-1683.0 + 1158.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,4.0)*Power(xij,4.0)) + 
           510.0*Power(xii,10.0)*Power(xij,10.0)*
            (-83889.0 - 7948.0*Power(rij,2.0)*Power(xij,2.0) + 
              12.0*Power(rij,4.0)*Power(xij,4.0)) - 
           34.0*rij*Power(xii,11.0)*Power(xij,10.0)*
            (-1158885.0 + 3450.0*Power(rij,2.0)*Power(xij,2.0) + 
              16.0*Power(rij,4.0)*Power(xij,4.0)) - 
           90.0*Power(xii,20.0)*(3876.0 + 10354.0*Power(rij,2.0)*Power(xij,2.0) + 
              29.0*Power(rij,4.0)*Power(xij,4.0)) + 
           1020.0*Power(xii,12.0)*Power(xij,8.0)*
            (-172098.0 - 26.0*Power(rij,2.0)*Power(xij,2.0) + 
              31.0*Power(rij,4.0)*Power(xij,4.0)) - 
           1020.0*Power(xii,14.0)*Power(xij,6.0)*
            (210168.0 - 8596.0*Power(rij,2.0)*Power(xij,2.0) + 
              39.0*Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*rij*Power(xii,21.0)*(-87210.0 - 43125.0*Power(rij,2.0)*Power(xij,2.0) + 
              47.0*Power(rij,4.0)*Power(xij,4.0)) - 
           15.0*rij*Power(xii,17.0)*Power(xij,4.0)*
            (1992273.0 - 31144.0*Power(rij,2.0)*Power(xij,2.0) + 
              68.0*Power(rij,4.0)*Power(xij,4.0)) - 
           90.0*Power(xii,8.0)*Power(xij,12.0)*
            (17425.0 + 6664.0*Power(rij,2.0)*Power(xij,2.0) + 
              76.0*Power(rij,4.0)*Power(xij,4.0)) + 
           rij*Power(xii,19.0)*Power(xij,2.0)*
            (-5204385.0 - 202710.0*Power(rij,2.0)*Power(xij,2.0) + 
              544.0*Power(rij,4.0)*Power(xij,4.0)) + 
           45.0*Power(xii,18.0)*Power(xij,2.0)*
            (-267615.0 - 83676.0*Power(rij,2.0)*Power(xij,2.0) + 
              680.0*Power(rij,4.0)*Power(xij,4.0)) - 
           15.0*Power(xii,16.0)*Power(xij,4.0)*
            (6000651.0 - 41616.0*Power(rij,2.0)*Power(xij,2.0) + 
              952.0*Power(rij,4.0)*Power(xij,4.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,8.0)*
         (2.0*Power(xii,2.0)*Power(xij,24.0)*
            (436049563950.0 + 402658381125.0*rij*xij + 
              173330907750.0*Power(rij,2.0)*Power(xij,2.0) + 
              45555359850.0*Power(rij,3.0)*Power(xij,3.0) + 
              7994586600.0*Power(rij,4.0)*Power(xij,4.0) + 
              948782835.0*Power(rij,5.0)*Power(xij,5.0) + 
              69999930.0*Power(rij,6.0)*Power(xij,6.0) + 
              1737450.0*Power(rij,7.0)*Power(xij,7.0) - 
              254430.0*Power(rij,8.0)*Power(xij,8.0) - 
              34155.0*Power(rij,9.0)*Power(xij,9.0) - 
              1914.0*Power(rij,10.0)*Power(xij,10.0) - 46.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 44.0*Power(xii,20.0)*Power(xij,6.0)*
            (-43375500.0 - 79521750.0*rij*xij - 
              72292500.0*Power(rij,2.0)*Power(xij,2.0) - 
              43375500.0*Power(rij,3.0)*Power(xij,3.0) - 
              19278000.0*Power(rij,4.0)*Power(xij,4.0) - 
              6747300.0*Power(rij,5.0)*Power(xij,5.0) - 
              1927800.0*Power(rij,6.0)*Power(xij,6.0) - 
              441180.0*Power(rij,7.0)*Power(xij,7.0) - 
              109620.0*Power(rij,8.0)*Power(xij,8.0) - 
              14715.0*Power(rij,9.0)*Power(xij,9.0) - 
              690.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) + 
           6.0*Power(xij,26.0)*(6547290750.0 + 7202019825.0*rij*xij + 
              3790536750.0*Power(rij,2.0)*Power(xij,2.0) + 
              1263512250.0*Power(rij,3.0)*Power(xij,3.0) + 
              297297000.0*Power(rij,4.0)*Power(xij,4.0) + 
              52026975.0*Power(rij,5.0)*Power(xij,5.0) + 
              6936930.0*Power(rij,6.0)*Power(xij,6.0) + 
              707850.0*Power(rij,7.0)*Power(xij,7.0) + 
              54450.0*Power(rij,8.0)*Power(xij,8.0) + 
              3025.0*Power(rij,9.0)*Power(xij,9.0) + 
              110.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) + 
           44.0*Power(xii,6.0)*Power(xij,20.0)*
            (100049928300.0 - 5205782925.0*rij*xij - 
              25852279950.0*Power(rij,2.0)*Power(xij,2.0) - 
              8238935250.0*Power(rij,3.0)*Power(xij,3.0) - 
              784614600.0*Power(rij,4.0)*Power(xij,4.0) + 
              136745280.0*Power(rij,5.0)*Power(xij,5.0) + 
              52950240.0*Power(rij,6.0)*Power(xij,6.0) + 
              7931520.0*Power(rij,7.0)*Power(xij,7.0) + 
              685440.0*Power(rij,8.0)*Power(xij,8.0) + 
              34425.0*Power(rij,9.0)*Power(xij,9.0) + 
              822.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) - 
           3.0*Power(xii,26.0)*(935550.0 + 1715175.0*rij*xij + 
              1559250.0*Power(rij,2.0)*Power(xij,2.0) + 
              935550.0*Power(rij,3.0)*Power(xij,3.0) + 
              415800.0*Power(rij,4.0)*Power(xij,4.0) + 
              145530.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              9900.0*Power(rij,7.0)*Power(xij,7.0) + 1980.0*Power(rij,8.0)*Power(xij,8.0) + 
              330.0*Power(rij,9.0)*Power(xij,9.0) + 44.0*Power(rij,10.0)*Power(xij,10.0) + 
              4.0*Power(rij,11.0)*Power(xij,11.0)) + 
           2244.0*Power(xii,14.0)*Power(xij,12.0)*
            (-15479100.0 - 28676025.0*rij*xij - 
              22821750.0*Power(rij,2.0)*Power(xij,2.0) - 
              22689450.0*Power(rij,3.0)*Power(xij,3.0) - 
              1852200.0*Power(rij,4.0)*Power(xij,4.0) - 
              2372580.0*Power(rij,5.0)*Power(xij,5.0) - 
              1252440.0*Power(rij,6.0)*Power(xij,6.0) - 
              228600.0*Power(rij,7.0)*Power(xij,7.0) - 
              15000.0*Power(rij,8.0)*Power(xij,8.0) + 450.0*Power(rij,9.0)*Power(xij,9.0) + 
              108.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           2244.0*Power(xii,12.0)*Power(xij,14.0)*
            (-27556200.0 - 14104125.0*rij*xij - 
              108438750.0*Power(rij,2.0)*Power(xij,2.0) + 
              15375150.0*Power(rij,3.0)*Power(xij,3.0) - 
              5632200.0*Power(rij,4.0)*Power(xij,4.0) - 
              8370180.0*Power(rij,5.0)*Power(xij,5.0) - 
              2119320.0*Power(rij,6.0)*Power(xij,6.0) - 
              198000.0*Power(rij,7.0)*Power(xij,7.0) + 
              2400.0*Power(rij,8.0)*Power(xij,8.0) + 2010.0*Power(rij,9.0)*Power(xij,9.0) + 
              156.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) + 
           330.0*Power(xii,18.0)*Power(xij,8.0)*
            (-20241900.0 - 37110150.0*rij*xij - 
              33736500.0*Power(rij,2.0)*Power(xij,2.0) - 
              20241900.0*Power(rij,3.0)*Power(xij,3.0) - 
              8996400.0*Power(rij,4.0)*Power(xij,4.0) - 
              3211803.0*Power(rij,5.0)*Power(xij,5.0) - 
              773514.0*Power(rij,6.0)*Power(xij,6.0) - 
              263898.0*Power(rij,7.0)*Power(xij,7.0) - 
              53202.0*Power(rij,8.0)*Power(xij,8.0) - 
              4393.0*Power(rij,9.0)*Power(xij,9.0) - 62.0*Power(rij,10.0)*Power(xij,10.0) + 
              6.0*Power(rij,11.0)*Power(xij,11.0)) - 
           165.0*Power(xii,8.0)*Power(xij,18.0)*
            (-11754743490.0 + 11330341155.0*rij*xij + 
              1384290810.0*Power(rij,2.0)*Power(xij,2.0) - 
              2116476810.0*Power(rij,3.0)*Power(xij,3.0) - 
              782225640.0*Power(rij,4.0)*Power(xij,4.0) - 
              97437186.0*Power(rij,5.0)*Power(xij,5.0) + 
              2679012.0*Power(rij,6.0)*Power(xij,6.0) + 
              2436804.0*Power(rij,7.0)*Power(xij,7.0) + 
              347316.0*Power(rij,8.0)*Power(xij,8.0) + 
              25014.0*Power(rij,9.0)*Power(xij,9.0) + 
              916.0*Power(rij,10.0)*Power(xij,10.0) + 12.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 4.0*Power(xii,4.0)*Power(xij,22.0)*
            (921052717200.0 + 543777678675.0*rij*xij + 
              99905825250.0*Power(rij,2.0)*Power(xij,2.0) - 
              10883876850.0*Power(rij,3.0)*Power(xij,3.0) - 
              9266934600.0*Power(rij,4.0)*Power(xij,4.0) - 
              2236505040.0*Power(rij,5.0)*Power(xij,5.0) - 
              316673280.0*Power(rij,6.0)*Power(xij,6.0) - 
              28779300.0*Power(rij,7.0)*Power(xij,7.0) - 
              1601820.0*Power(rij,8.0)*Power(xij,8.0) - 
              40095.0*Power(rij,9.0)*Power(xij,9.0) + 
              726.0*Power(rij,10.0)*Power(xij,10.0) + 58.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 4.0*Power(xii,22.0)*Power(xij,4.0)*
            (95426100.0 + 174947850.0*rij*xij + 
              159043500.0*Power(rij,2.0)*Power(xij,2.0) + 
              95426100.0*Power(rij,3.0)*Power(xij,3.0) + 
              42411600.0*Power(rij,4.0)*Power(xij,4.0) + 
              14844060.0*Power(rij,5.0)*Power(xij,5.0) + 
              4241160.0*Power(rij,6.0)*Power(xij,6.0) + 
              1009800.0*Power(rij,7.0)*Power(xij,7.0) + 
              201960.0*Power(rij,8.0)*Power(xij,8.0) + 
              37125.0*Power(rij,9.0)*Power(xij,9.0) + 
              3102.0*Power(rij,10.0)*Power(xij,10.0) + 58.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 66.0*Power(xii,16.0)*Power(xij,10.0)*
            (-263144700.0 - 482431950.0*rij*xij - 
              438574500.0*Power(rij,2.0)*Power(xij,2.0) - 
              259704900.0*Power(rij,3.0)*Power(xij,3.0) - 
              130712400.0*Power(rij,4.0)*Power(xij,4.0) - 
              27031095.0*Power(rij,5.0)*Power(xij,5.0) - 
              13816530.0*Power(rij,6.0)*Power(xij,6.0) - 
              4240170.0*Power(rij,7.0)*Power(xij,7.0) - 
              537330.0*Power(rij,8.0)*Power(xij,8.0) - 
              20565.0*Power(rij,9.0)*Power(xij,9.0) + 
              1146.0*Power(rij,10.0)*Power(xij,10.0) + 86.0*Power(rij,11.0)*Power(xij,11.0)) \
    + Power(xii,24.0)*Power(xij,2.0)*(47713050.0 + 87473925.0*rij*xij + 
              79521750.0*Power(rij,2.0)*Power(xij,2.0) + 
              47713050.0*Power(rij,3.0)*Power(xij,3.0) + 
              21205800.0*Power(rij,4.0)*Power(xij,4.0) + 
              7422030.0*Power(rij,5.0)*Power(xij,5.0) + 
              2120580.0*Power(rij,6.0)*Power(xij,6.0) + 
              504900.0*Power(rij,7.0)*Power(xij,7.0) + 
              100980.0*Power(rij,8.0)*Power(xij,8.0) + 
              16830.0*Power(rij,9.0)*Power(xij,9.0) + 
              2244.0*Power(rij,10.0)*Power(xij,10.0) + 92.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 33.0*Power(xii,10.0)*Power(xij,16.0)*
            (5519319750.0 - 27722883825.0*rij*xij + 
              11646151650.0*Power(rij,2.0)*Power(xij,2.0) + 
              955234350.0*Power(rij,3.0)*Power(xij,3.0) - 
              2729953800.0*Power(rij,4.0)*Power(xij,4.0) - 
              902572650.0*Power(rij,5.0)*Power(xij,5.0) - 
              105286860.0*Power(rij,6.0)*Power(xij,6.0) + 
              622260.0*Power(rij,7.0)*Power(xij,7.0) + 
              1538340.0*Power(rij,8.0)*Power(xij,8.0) + 
              178830.0*Power(rij,9.0)*Power(xij,9.0) + 
              9060.0*Power(rij,10.0)*Power(xij,10.0) + 172.0*Power(rij,11.0)*Power(xij,11.0))\
    ))/(2.80665e6*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),17.0))
    ;
  }
  return S;
}

double Slater_3S_1S(double rij,double xii,double xij)
{
  return Slater_1S_3S(rij,xij,xii);
}

double Slater_3S_2S(double rij,double xii,double xij)
{
  return Slater_2S_3S(rij,xij,xii);
}

static double Slater_4S_4S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-83691159552000.0 + 83691159552000.0*Power(E,2.0*rij*xii) - 
        150568359566625.0*rij*xii - 133754400029250.0*Power(rij,2.0)*Power(xii,2.0) - 
        78142908343500.0*Power(rij,3.0)*Power(xii,3.0) - 
        33740723016000.0*Power(rij,4.0)*Power(xii,4.0) - 
        11470756096800.0*Power(rij,5.0)*Power(xii,5.0) - 
        3193358968800.0*Power(rij,6.0)*Power(xii,6.0) - 
        747112766400.0*Power(rij,7.0)*Power(xii,7.0) - 
        149448499200.0*Power(rij,8.0)*Power(xii,8.0) - 
        25830604800.0*Power(rij,9.0)*Power(xii,9.0) - 
        3874590720.0*Power(rij,10.0)*Power(xii,10.0) - 
        503193600.0*Power(rij,11.0)*Power(xii,11.0) - 
        55910400.0*Power(rij,12.0)*Power(xii,12.0) - 
        5160960.0*Power(rij,13.0)*Power(xii,13.0) - 
        368640.0*Power(rij,14.0)*Power(xii,14.0) - 16384.0*Power(rij,15.0)*Power(xii,15.0))/
      (8.3691159552e13*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (1260.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),15.0) + 
        Power(E,2.0*rij*xij)*Power(xij,10.0)*
         (-3276.0*Power(rij,5.0)*Power(xii,25.0) - 168.0*Power(rij,6.0)*Power(xii,26.0) - 
           4.0*Power(rij,7.0)*Power(xii,27.0) + 1260.0*Power(xij,20.0) + 
           2205.0*rij*xii*Power(xij,20.0) + 
           1890.0*Power(xii,2.0)*Power(xij,18.0)*(-10.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           420.0*Power(rij,4.0)*Power(xii,24.0)*(91.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           525.0*rij*Power(xii,3.0)*Power(xij,18.0)*
            (-63.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           42.0*Power(rij,3.0)*Power(xii,23.0)*
            (-6825.0 - 405.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           63.0*rij*Power(xii,5.0)*Power(xij,16.0)*
            (3675.0 - 250.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           210.0*Power(xii,4.0)*Power(xij,16.0)*
            (630.0 - 135.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 252.0*Power(rij,2.0)*Power(xii,22.0)*
            (-5460.0 - 1225.0*Power(rij,2.0)*Power(xij,2.0) + 
              17.0*Power(rij,4.0)*Power(xij,4.0)) - 
           1260.0*rij*Power(xii,17.0)*Power(xij,4.0)*
            (141729.0 - 10145.0*Power(rij,2.0)*Power(xij,2.0) + 
              116.0*Power(rij,4.0)*Power(xij,4.0)) + 
           21.0*rij*Power(xii,9.0)*Power(xij,12.0)*
            (164775.0 - 18460.0*Power(rij,2.0)*Power(xij,2.0) + 
              828.0*Power(rij,4.0)*Power(xij,4.0)) + 
           14.0*Power(xii,6.0)*Power(xij,14.0)*
            (-40950.0 + 14175.0*Power(rij,2.0)*Power(xij,2.0) - 
              450.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) - 
           210.0*Power(xii,8.0)*Power(xij,12.0)*
            (-8190.0 + 4095.0*Power(rij,2.0)*Power(xij,2.0) - 
              210.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           42.0*Power(xii,10.0)*Power(xij,10.0)*
            (-209430.0 - 2925.0*Power(rij,2.0)*Power(xij,2.0) - 
              8840.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           rij*Power(xii,7.0)*Power(xij,14.0)*
            (-1003275.0 + 110250.0*Power(rij,2.0)*Power(xij,2.0) - 
              1890.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           21.0*rij*Power(xii,11.0)*Power(xij,10.0)*
            (-1033695.0 - 218400.0*Power(rij,2.0)*Power(xij,2.0) + 
              552.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           280.0*Power(xii,18.0)*Power(xij,2.0)*
            (-385560.0 - 73953.0*Power(rij,2.0)*Power(xij,2.0) + 
              2370.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           35.0*rij*Power(xii,15.0)*Power(xij,6.0)*
            (-1565613.0 + 359520.0*Power(rij,2.0)*Power(xij,2.0) - 
              7020.0*Power(rij,4.0)*Power(xij,4.0) + 8.0*Power(rij,6.0)*Power(xij,6.0)) + 
           14.0*rij*Power(xii,19.0)*Power(xij,2.0)*
            (-4980150.0 + 126765.0*Power(rij,2.0)*Power(xij,2.0) - 
              3852.0*Power(rij,4.0)*Power(xij,4.0) + 20.0*Power(rij,6.0)*Power(xij,6.0)) - 
           630.0*Power(xii,14.0)*Power(xij,6.0)*
            (708714.0 - 14385.0*Power(rij,2.0)*Power(xij,2.0) - 
              2340.0*Power(rij,4.0)*Power(xij,4.0) + 20.0*Power(rij,6.0)*Power(xij,6.0)) + 
           210.0*Power(xii,16.0)*Power(xij,4.0)*
            (-2087532.0 + 328491.0*Power(rij,2.0)*Power(xij,2.0) - 
              11740.0*Power(rij,4.0)*Power(xij,4.0) + 52.0*Power(rij,6.0)*Power(xij,6.0)) - 
           84.0*Power(xii,20.0)*(59670.0 + 236250.0*Power(rij,2.0)*Power(xij,2.0) - 
              8745.0*Power(rij,4.0)*Power(xij,4.0) + 92.0*Power(rij,6.0)*Power(xij,6.0)) - 
           2.0*rij*Power(xii,21.0)*(1949220.0 + 1598625.0*Power(rij,2.0)*Power(xij,2.0) - 
              41391.0*Power(rij,4.0)*Power(xij,4.0) + 128.0*Power(rij,6.0)*Power(xij,6.0)) \
    + rij*Power(xii,13.0)*Power(xij,8.0)*
            (173037375.0 - 2784600.0*Power(rij,2.0)*Power(xij,2.0) - 
              112140.0*Power(rij,4.0)*Power(xij,4.0) + 256.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 14.0*Power(xii,12.0)*Power(xij,8.0)*
            (-7260750.0 - 2521935.0*Power(rij,2.0)*Power(xij,2.0) + 
              19500.0*Power(rij,4.0)*Power(xij,4.0) + 344.0*Power(rij,6.0)*Power(xij,6.0))) \
    + Power(E,2.0*rij*xii)*Power(xii,10.0)*
         (210.0*Power(xii,2.0)*Power(xij,18.0)*
            (514080.0 + 332010.0*rij*xij + 94500.0*Power(rij,2.0)*Power(xij,2.0) + 
              15225.0*Power(rij,3.0)*Power(xij,3.0) + 
              1470.0*Power(rij,4.0)*Power(xij,4.0) + 81.0*Power(rij,5.0)*Power(xij,5.0) + 
              2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           105.0*Power(xii,18.0)*Power(xij,2.0)*
            (180.0 + 315.0*rij*xij + 270.0*Power(rij,2.0)*Power(xij,2.0) + 
              150.0*Power(rij,3.0)*Power(xij,3.0) + 60.0*Power(rij,4.0)*Power(xij,4.0) + 
              18.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           1365.0*Power(xii,10.0)*Power(xij,10.0)*
            (-6444.0 + 15903.0*rij*xij - 25866.0*Power(rij,2.0)*Power(xij,2.0) - 
              2040.0*Power(rij,3.0)*Power(xij,3.0) + 1080.0*Power(rij,4.0)*Power(xij,4.0) + 
              180.0*Power(rij,5.0)*Power(xij,5.0) + 8.0*Power(rij,6.0)*Power(xij,6.0)) + 
           Power(xii,14.0)*Power(xij,6.0)*
            (573300.0 + 1003275.0*rij*xij + 859950.0*Power(rij,2.0)*Power(xij,2.0) + 
              387660.0*Power(rij,3.0)*Power(xij,3.0) + 
              371280.0*Power(rij,4.0)*Power(xij,4.0) + 
              11592.0*Power(rij,5.0)*Power(xij,5.0) - 
              4816.0*Power(rij,6.0)*Power(xij,6.0) - 256.0*Power(rij,7.0)*Power(xij,7.0)) + 
           2.0*Power(xij,20.0)*(2506140.0 + 1949220.0*rij*xij + 
              687960.0*Power(rij,2.0)*Power(xij,2.0) + 
              143325.0*Power(rij,3.0)*Power(xij,3.0) + 
              19110.0*Power(rij,4.0)*Power(xij,4.0) + 
              1638.0*Power(rij,5.0)*Power(xij,5.0) + 84.0*Power(rij,6.0)*Power(xij,6.0) + 
              2.0*Power(rij,7.0)*Power(xij,7.0)) - 
           42.0*Power(xii,4.0)*Power(xij,16.0)*
            (-10437660.0 - 4251870.0*rij*xij - 493020.0*Power(rij,2.0)*Power(xij,2.0) + 
              42255.0*Power(rij,3.0)*Power(xij,3.0) + 
              17490.0*Power(rij,4.0)*Power(xij,4.0) + 
              1971.0*Power(rij,5.0)*Power(xij,5.0) + 102.0*Power(rij,6.0)*Power(xij,6.0) + 
              2.0*Power(rij,7.0)*Power(xij,7.0)) + 
           21.0*Power(xii,16.0)*Power(xij,4.0)*
            (-6300.0 - 11025.0*rij*xij - 9450.0*Power(rij,2.0)*Power(xij,2.0) - 
              5250.0*Power(rij,3.0)*Power(xij,3.0) - 2100.0*Power(rij,4.0)*Power(xij,4.0) - 
              828.0*Power(rij,5.0)*Power(xij,5.0) - 8.0*Power(rij,6.0)*Power(xij,6.0) + 
              4.0*Power(rij,7.0)*Power(xij,7.0)) - 
           Power(xii,20.0)*(1260.0 + 2205.0*rij*xij + 
              1890.0*Power(rij,2.0)*Power(xij,2.0) + 1050.0*Power(rij,3.0)*Power(xij,3.0) + 
              420.0*Power(rij,4.0)*Power(xij,4.0) + 126.0*Power(rij,5.0)*Power(xij,5.0) + 
              28.0*Power(rij,6.0)*Power(xij,6.0) + 4.0*Power(rij,7.0)*Power(xij,7.0)) - 
           35.0*Power(xii,8.0)*Power(xij,12.0)*
            (-2904300.0 + 4943925.0*rij*xij + 258930.0*Power(rij,2.0)*Power(xij,2.0) - 
              359520.0*Power(rij,3.0)*Power(xij,3.0) - 
              70440.0*Power(rij,4.0)*Power(xij,4.0) - 
              4176.0*Power(rij,5.0)*Power(xij,5.0) + 32.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0)) + 
           35.0*Power(xii,12.0)*Power(xij,8.0)*
            (-49140.0 - 98865.0*rij*xij + 3510.0*Power(rij,2.0)*Power(xij,2.0) - 
              131040.0*Power(rij,3.0)*Power(xij,3.0) - 
              7800.0*Power(rij,4.0)*Power(xij,4.0) + 3204.0*Power(rij,5.0)*Power(xij,5.0) + 
              360.0*Power(rij,6.0)*Power(xij,6.0) + 8.0*Power(rij,7.0)*Power(xij,7.0)) + 
           Power(xii,6.0)*Power(xij,14.0)*
            (446489820.0 - 54796455.0*rij*xij - 68983110.0*Power(rij,2.0)*Power(xij,2.0) - 
              12782700.0*Power(rij,3.0)*Power(xij,3.0) - 
              663600.0*Power(rij,4.0)*Power(xij,4.0) + 
              53928.0*Power(rij,5.0)*Power(xij,5.0) + 7728.0*Power(rij,6.0)*Power(xij,6.0) + 
              256.0*Power(rij,7.0)*Power(xij,7.0))))/
      (1260.*Power(E,2.0*rij*(xii + xij))*rij*Power(Power(xii,2.0) - Power(xij,2.0),15.0))
    ;
  }
  return S;
}

static double Slater_4S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-14227497123840000.0 + 14227497123840000.0*Power(E,2.0*rij*xii) - 
        25913502934444125.0*rij*xii - 
        23372011621208250.0*Power(rij,2.0)*Power(xii,2.0) - 
        13907709869303250.0*Power(rij,3.0)*Power(xii,3.0) - 
        6137735659555500.0*Power(rij,4.0)*Power(xii,4.0) - 
        2140857388870200.0*Power(rij,5.0)*Power(xii,5.0) - 
        614116575072000.0*Power(rij,6.0)*Power(xii,6.0) - 
        148809580920000.0*Power(rij,7.0)*Power(xii,7.0) - 
        31036639233600.0*Power(rij,8.0)*Power(xii,8.0) - 
        5645342102400.0*Power(rij,9.0)*Power(xii,9.0) - 
        903333150720.0*Power(rij,10.0)*Power(xii,10.0) - 
        127744081920.0*Power(rij,11.0)*Power(xii,11.0) - 
        15968010240.0*Power(rij,12.0)*Power(xii,12.0) - 
        1754726400.0*Power(rij,13.0)*Power(xii,13.0) - 
        167116800.0*Power(rij,14.0)*Power(xii,14.0) - 
        13369344.0*Power(rij,15.0)*Power(xii,15.0) - 
        835584.0*Power(rij,16.0)*Power(xii,16.0) - 32768.0*Power(rij,17.0)*Power(xii,17.0))/
      (1.422749712384e16*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (56700.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),17.0) + 
        9.0*Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-980.0*Power(rij,6.0)*Power(xii,28.0) - 20.0*Power(rij,7.0)*Power(xii,29.0) + 
           6300.0*Power(xij,22.0) + 11025.0*rij*xii*Power(xij,22.0) - 
           50.0*Power(rij,5.0)*Power(xii,27.0)*(441.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3150.0*Power(xii,2.0)*Power(xij,20.0)*
            (-34.0 + 3.0*Power(rij,2.0)*Power(xij,2.0)) + 
           525.0*rij*Power(xii,3.0)*Power(xij,20.0)*
            (-357.0 + 10.0*Power(rij,2.0)*Power(xij,2.0)) - 
           420.0*Power(rij,4.0)*Power(xii,26.0)*
            (700.0 + 19.0*Power(rij,2.0)*Power(xij,2.0)) + 
           1050.0*Power(xii,4.0)*Power(xij,18.0)*
            (816.0 - 153.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 210.0*rij*Power(xii,5.0)*Power(xij,18.0)*
            (7140.0 - 425.0*Power(rij,2.0)*Power(xij,2.0) + 
              3.0*Power(rij,4.0)*Power(xij,4.0)) + 
           42.0*Power(rij,3.0)*Power(xii,25.0)*
            (-59500.0 - 6035.0*Power(rij,2.0)*Power(xij,2.0) + 
              18.0*Power(rij,4.0)*Power(xij,4.0)) + 
           84.0*Power(rij,2.0)*Power(xii,24.0)*
            (-160650.0 - 52700.0*Power(rij,2.0)*Power(xij,2.0) + 
              397.0*Power(rij,4.0)*Power(xij,4.0)) - 
           28.0*Power(xii,12.0)*Power(xij,10.0)*
            (100849950.0 + 27100125.0*Power(rij,2.0)*Power(xij,2.0) + 
              186150.0*Power(rij,4.0)*Power(xij,4.0) - 2177.0*Power(rij,6.0)*Power(xij,6.0)\
    ) + 140.0*Power(xii,6.0)*Power(xij,16.0)*
            (-30600.0 + 9180.0*Power(rij,2.0)*Power(xij,2.0) - 
              255.0*Power(rij,4.0)*Power(xij,4.0) + Power(rij,6.0)*Power(xij,6.0)) - 
           2380.0*Power(xii,8.0)*Power(xij,14.0)*
            (-6300.0 + 2700.0*Power(rij,2.0)*Power(xij,2.0) - 
              120.0*Power(rij,4.0)*Power(xij,4.0) + Power(rij,6.0)*Power(xij,6.0)) + 
           10.0*rij*Power(xii,7.0)*Power(xij,16.0)*
            (-749700.0 + 71400.0*Power(rij,2.0)*Power(xij,2.0) - 
              1071.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           204.0*rij*Power(xii,15.0)*Power(xij,8.0)*
            (28962255.0 - 1744750.0*Power(rij,2.0)*Power(xij,2.0) + 
              9555.0*Power(rij,4.0)*Power(xij,4.0) + 6.0*Power(rij,6.0)*Power(xij,6.0)) - 
           42.0*rij*Power(xii,11.0)*Power(xij,12.0)*
            (-12911925.0 - 1634550.0*Power(rij,2.0)*Power(xij,2.0) - 
              7103.0*Power(rij,4.0)*Power(xij,4.0) + 18.0*Power(rij,6.0)*Power(xij,6.0)) + 
           2.0*rij*Power(xii,9.0)*Power(xij,14.0)*
            (16948575.0 - 1184400.0*Power(rij,2.0)*Power(xij,2.0) + 
              63861.0*Power(rij,4.0)*Power(xij,4.0) + 50.0*Power(rij,6.0)*Power(xij,6.0)) + 
           28.0*Power(xii,22.0)*(-2180250.0 - 10993050.0*Power(rij,2.0)*Power(xij,2.0) + 
              14925.0*Power(rij,4.0)*Power(xij,4.0) + 73.0*Power(rij,6.0)*Power(xij,6.0)) - 
           952.0*Power(xii,14.0)*Power(xij,8.0)*
            (16966215.0 + 725175.0*Power(rij,2.0)*Power(xij,2.0) - 
              36075.0*Power(rij,4.0)*Power(xij,4.0) + 79.0*Power(rij,6.0)*Power(xij,6.0)) - 
           84.0*Power(xii,10.0)*Power(xij,12.0)*
            (1723800.0 + 279225.0*Power(rij,2.0)*Power(xij,2.0) + 
              45600.0*Power(rij,4.0)*Power(xij,4.0) + 107.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 35.0*rij*Power(xii,17.0)*Power(xij,6.0)*
            (132637869.0 - 2205240.0*Power(rij,2.0)*Power(xij,2.0) - 
              48348.0*Power(rij,4.0)*Power(xij,4.0) + 136.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 6.0*rij*Power(xii,21.0)*Power(xij,2.0)*
            (192298050.0 + 12644275.0*Power(rij,2.0)*Power(xij,2.0) - 
              218029.0*Power(rij,4.0)*Power(xij,4.0) + 204.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 4.0*rij*Power(xii,13.0)*Power(xij,10.0)*
            (1259522775.0 + 15895425.0*Power(rij,2.0)*Power(xij,2.0) - 
              493017.0*Power(rij,4.0)*Power(xij,4.0) + 263.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 140.0*Power(xii,16.0)*Power(xij,6.0)*
            (180826281.0 - 15101406.0*Power(rij,2.0)*Power(xij,2.0) + 
              160140.0*Power(rij,4.0)*Power(xij,4.0) + 442.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 2.0*rij*Power(xii,23.0)*(21366450.0 + 23526300.0*Power(rij,2.0)*Power(xij,2.0) - 
              246729.0*Power(rij,4.0)*Power(xij,4.0) + 526.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 7.0*rij*Power(xii,19.0)*Power(xij,4.0)*
            (-811081215.0 + 39095550.0*Power(rij,2.0)*Power(xij,2.0) - 
              515916.0*Power(rij,4.0)*Power(xij,4.0) + 680.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 70.0*Power(xii,18.0)*Power(xij,4.0)*
            (-180554454.0 + 9873711.0*Power(rij,2.0)*Power(xij,2.0) - 
              414120.0*Power(rij,4.0)*Power(xij,4.0) + 2924.0*Power(rij,6.0)*Power(xij,6.0)\
    ) - 14.0*Power(xii,20.0)*Power(xij,2.0)*
            (136919700.0 + 71867115.0*Power(rij,2.0)*Power(xij,2.0) - 
              2154150.0*Power(rij,4.0)*Power(xij,4.0) + 
              10268.0*Power(rij,6.0)*Power(xij,6.0))) - 
        4.0*Power(E,2.0*rij*xii)*Power(xii,10.0)*
         (-10710.0*Power(xii,12.0)*Power(xij,12.0)*
            (-3555.0 - 127008.0*rij*xij + 138384.0*Power(rij,2.0)*Power(xij,2.0) - 
              74556.0*Power(rij,3.0)*Power(xij,3.0) - 
              22284.0*Power(rij,4.0)*Power(xij,4.0) + 408.0*Power(rij,5.0)*Power(xij,5.0) + 
              576.0*Power(rij,6.0)*Power(xij,6.0) + 60.0*Power(rij,7.0)*Power(xij,7.0) + 
              2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           2.0*Power(xii,20.0)*Power(xij,4.0)*
            (963900.0 + 1735020.0*rij*xij + 1542240.0*Power(rij,2.0)*Power(xij,2.0) + 
              899640.0*Power(rij,3.0)*Power(xij,3.0) + 
              385560.0*Power(rij,4.0)*Power(xij,4.0) + 
              128520.0*Power(rij,5.0)*Power(xij,5.0) + 
              34272.0*Power(rij,6.0)*Power(xij,6.0) + 
              9126.0*Power(rij,7.0)*Power(xij,7.0) + 333.0*Power(rij,8.0)*Power(xij,8.0) - 
              20.0*Power(rij,9.0)*Power(xij,9.0)) - 
           2.0*Power(xij,24.0)*(119041650.0 + 107137485.0*rij*xij + 
              45110520.0*Power(rij,2.0)*Power(xij,2.0) + 
              11695320.0*Power(rij,3.0)*Power(xij,3.0) + 
              2063880.0*Power(rij,4.0)*Power(xij,4.0) + 
              257985.0*Power(rij,5.0)*Power(xij,5.0) + 
              22932.0*Power(rij,6.0)*Power(xij,6.0) + 
              1404.0*Power(rij,7.0)*Power(xij,7.0) + 54.0*Power(rij,8.0)*Power(xij,8.0) + 
              Power(rij,9.0)*Power(xij,9.0)) + 
           2.0*Power(xii,2.0)*Power(xij,22.0)*
            (-3264488325.0 - 2505368880.0*rij*xij - 
              881390160.0*Power(rij,2.0)*Power(xij,2.0) - 
              185775660.0*Power(rij,3.0)*Power(xij,3.0) - 
              25639740.0*Power(rij,4.0)*Power(xij,4.0) - 
              2361555.0*Power(rij,5.0)*Power(xij,5.0) - 
              139356.0*Power(rij,6.0)*Power(xij,6.0) - 
              4482.0*Power(rij,7.0)*Power(xij,7.0) - 27.0*Power(rij,8.0)*Power(xij,8.0) + 
              2.0*Power(rij,9.0)*Power(xij,9.0)) + 
           Power(xii,24.0)*(14175.0 + 25515.0*rij*xij + 
              22680.0*Power(rij,2.0)*Power(xij,2.0) + 
              13230.0*Power(rij,3.0)*Power(xij,3.0) + 
              5670.0*Power(rij,4.0)*Power(xij,4.0) + 1890.0*Power(rij,5.0)*Power(xij,5.0) + 
              504.0*Power(rij,6.0)*Power(xij,6.0) + 108.0*Power(rij,7.0)*Power(xij,7.0) + 
              18.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) - 
           102.0*Power(xii,10.0)*Power(xij,14.0)*
            (44986725.0 - 97433280.0*rij*xij + 44467920.0*Power(rij,2.0)*Power(xij,2.0) + 
              15857100.0*Power(rij,3.0)*Power(xij,3.0) - 
              457380.0*Power(rij,4.0)*Power(xij,4.0) - 
              620550.0*Power(rij,5.0)*Power(xij,5.0) - 
              83160.0*Power(rij,6.0)*Power(xij,6.0) - 
              4068.0*Power(rij,7.0)*Power(xij,7.0) - 6.0*Power(rij,8.0)*Power(xij,8.0) + 
              4.0*Power(rij,9.0)*Power(xij,9.0)) + 
           102.0*Power(xii,14.0)*Power(xij,10.0)*
            (-859950.0 - 1437345.0*rij*xij - 2260440.0*Power(rij,2.0)*Power(xij,2.0) + 
              810810.0*Power(rij,3.0)*Power(xij,3.0) - 
              1056510.0*Power(rij,4.0)*Power(xij,4.0) - 
              217854.0*Power(rij,5.0)*Power(xij,5.0) + 
              6552.0*Power(rij,6.0)*Power(xij,6.0) + 3852.0*Power(rij,7.0)*Power(xij,7.0) + 
              258.0*Power(rij,8.0)*Power(xij,8.0) + 4.0*Power(rij,9.0)*Power(xij,9.0)) - 
           Power(xii,22.0)*Power(xij,2.0)*
            (240975.0 + 433755.0*rij*xij + 385560.0*Power(rij,2.0)*Power(xij,2.0) + 
              224910.0*Power(rij,3.0)*Power(xij,3.0) + 
              96390.0*Power(rij,4.0)*Power(xij,4.0) + 
              32130.0*Power(rij,5.0)*Power(xij,5.0) + 
              8568.0*Power(rij,6.0)*Power(xij,6.0) + 1836.0*Power(rij,7.0)*Power(xij,7.0) + 
              306.0*Power(rij,8.0)*Power(xij,8.0) + 4.0*Power(rij,9.0)*Power(xij,9.0)) + 
           2.0*Power(xii,4.0)*Power(xij,20.0)*
            (-18032978565.0 - 9823683240.0*rij*xij - 
              2047323600.0*Power(rij,2.0)*Power(xij,2.0) - 
              129098340.0*Power(rij,3.0)*Power(xij,3.0) + 
              26410860.0*Power(rij,4.0)*Power(xij,4.0) + 
              7094304.0*Power(rij,5.0)*Power(xij,5.0) + 
              788256.0*Power(rij,6.0)*Power(xij,6.0) + 
              48654.0*Power(rij,7.0)*Power(xij,7.0) + 
              1593.0*Power(rij,8.0)*Power(xij,8.0) + 20.0*Power(rij,9.0)*Power(xij,9.0)) - 
           6.0*Power(xii,16.0)*Power(xij,8.0)*
            (-5622750.0 - 10120950.0*rij*xij - 8996400.0*Power(rij,2.0)*Power(xij,2.0) - 
              5698350.0*Power(rij,3.0)*Power(xij,3.0) - 
              897750.0*Power(rij,4.0)*Power(xij,4.0) - 
              1641591.0*Power(rij,5.0)*Power(xij,5.0) - 
              211932.0*Power(rij,6.0)*Power(xij,6.0) + 
              10224.0*Power(rij,7.0)*Power(xij,7.0) + 
              2364.0*Power(rij,8.0)*Power(xij,8.0) + 73.0*Power(rij,9.0)*Power(xij,9.0)) + 
           2.0*Power(xii,18.0)*Power(xij,6.0)*
            (-4819500.0 - 8675100.0*rij*xij - 7711200.0*Power(rij,2.0)*Power(xij,2.0) - 
              4498200.0*Power(rij,3.0)*Power(xij,3.0) - 
              1927800.0*Power(rij,4.0)*Power(xij,4.0) - 
              561519.0*Power(rij,5.0)*Power(xij,5.0) - 
              279468.0*Power(rij,6.0)*Power(xij,6.0) - 
              20682.0*Power(rij,7.0)*Power(xij,7.0) + 
              1305.0*Power(rij,8.0)*Power(xij,8.0) + 106.0*Power(rij,9.0)*Power(xij,9.0)) + 
           3.0*Power(xii,8.0)*Power(xij,16.0)*
            (-9364244085.0 + 6940428705.0*rij*xij + 
              2117684520.0*Power(rij,2.0)*Power(xij,2.0) - 
              230268150.0*Power(rij,3.0)*Power(xij,3.0) - 
              149610510.0*Power(rij,4.0)*Power(xij,4.0) - 
              21824334.0*Power(rij,5.0)*Power(xij,5.0) - 
              1223208.0*Power(rij,6.0)*Power(xij,6.0) + 
              12708.0*Power(rij,7.0)*Power(xij,7.0) + 
              4470.0*Power(rij,8.0)*Power(xij,8.0) + 146.0*Power(rij,9.0)*Power(xij,9.0)) - 
           Power(xii,6.0)*Power(xij,18.0)*
            (57304872765.0 + 7147185255.0*rij*xij - 
              5801702760.0*Power(rij,2.0)*Power(xij,2.0) - 
              2053388610.0*Power(rij,3.0)*Power(xij,3.0) - 
              271655370.0*Power(rij,4.0)*Power(xij,4.0) - 
              10864854.0*Power(rij,5.0)*Power(xij,5.0) + 
              1337112.0*Power(rij,6.0)*Power(xij,6.0) + 
              202716.0*Power(rij,7.0)*Power(xij,7.0) + 
              10746.0*Power(rij,8.0)*Power(xij,8.0) + 212.0*Power(rij,9.0)*Power(xij,9.0))))/
      (56700.*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),17.0))
    ;
  }
  return S;
}

static double Slater_4S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-2919482409811968000.0 + 2919482409811968000.0*Power(E,2.0*rij*xii) - 
        5378825373422626125.0*rij*xii - 
        4918685927221316250.0*Power(rij,2.0)*Power(xii,2.0) - 
        2974825584766035000.0*Power(rij,3.0)*Power(xii,3.0) - 
        1337724873111627000.0*Power(rij,4.0)*Power(xii,4.0) - 
        476688322649038500.0*Power(rij,5.0)*Power(xii,5.0) - 
        140080945989184200.0*Power(rij,6.0)*Power(xii,6.0) - 
        34878402537778800.0*Power(rij,7.0)*Power(xii,7.0) - 
        7501749557702400.0*Power(rij,8.0)*Power(xii,8.0) - 
        1413711970070400.0*Power(rij,9.0)*Power(xii,9.0) - 
        235878458175360.0*Power(rij,10.0)*Power(xii,10.0) - 
        35103763618560.0*Power(rij,11.0)*Power(xii,11.0) - 
        4680908144640.0*Power(rij,12.0)*Power(xii,12.0) - 
        560108666880.0*Power(rij,13.0)*Power(xii,13.0) - 
        60011642880.0*Power(rij,14.0)*Power(xii,14.0) - 
        5715394560.0*Power(rij,15.0)*Power(xii,15.0) - 
        476282880.0*Power(rij,16.0)*Power(xii,16.0) - 
        33619968.0*Power(rij,17.0)*Power(xii,17.0) - 
        1867776.0*Power(rij,18.0)*Power(xii,18.0) - 65536.0*Power(rij,19.0)*Power(xii,19.0))/
      (2.919482409811968e18*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (1871100.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),19.0) + 
        495.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-672.0*Power(rij,6.0)*Power(xii,30.0) - 12.0*Power(rij,7.0)*Power(xii,31.0) + 
           3780.0*Power(xij,24.0) + 6615.0*rij*xii*Power(xij,24.0) - 
           136.0*Power(rij,5.0)*Power(xii,29.0)*(126.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           1890.0*Power(xii,2.0)*Power(xij,22.0)*
            (-38.0 + 3.0*Power(rij,2.0)*Power(xij,2.0)) + 
           315.0*rij*Power(xii,3.0)*Power(xij,22.0)*
            (-399.0 + 10.0*Power(rij,2.0)*Power(xij,2.0)) - 
           84.0*Power(rij,4.0)*Power(xii,28.0)*
            (3060.0 + 121.0*Power(rij,2.0)*Power(xij,2.0)) + 
           630.0*Power(xii,4.0)*Power(xij,20.0)*
            (1026.0 - 171.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           63.0*rij*Power(xii,5.0)*Power(xij,20.0)*
            (17955.0 - 950.0*Power(rij,2.0)*Power(xij,2.0) + 
              6.0*Power(rij,4.0)*Power(xij,4.0)) + 
           84.0*Power(rij,2.0)*Power(xii,26.0)*
            (-174420.0 - 71535.0*Power(rij,2.0)*Power(xij,2.0) + 
              179.0*Power(rij,4.0)*Power(xij,4.0)) - 
           63.0*rij*Power(xii,19.0)*Power(xij,6.0)*
            (468377895.0 - 14898090.0*Power(rij,2.0)*Power(xij,2.0) + 
              78812.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,27.0)*(-2441880.0*Power(rij,3.0) - 
              327978.0*Power(rij,5.0)*Power(xij,2.0) + 496.0*Power(rij,7.0)*Power(xij,4.0)) \
    + 2.0*rij*Power(xii,11.0)*Power(xij,14.0)*
            (613624095.0 + 56366730.0*Power(rij,2.0)*Power(xij,2.0) + 
              383607.0*Power(rij,4.0)*Power(xij,4.0) - 248.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 42.0*Power(xii,6.0)*Power(xij,18.0)*
            (-87210.0 + 23085.0*Power(rij,2.0)*Power(xij,2.0) - 
              570.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) - 
           798.0*Power(xii,8.0)*Power(xij,16.0)*
            (-18360.0 + 6885.0*Power(rij,2.0)*Power(xij,2.0) - 
              270.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           3.0*rij*Power(xii,7.0)*Power(xij,18.0)*
            (-2136645.0 + 179550.0*Power(rij,2.0)*Power(xij,2.0) - 
              2394.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           1596.0*Power(xii,14.0)*Power(xij,10.0)*
            (-34484670.0 - 2408985.0*Power(rij,2.0)*Power(xij,2.0) + 
              32810.0*Power(rij,4.0)*Power(xij,4.0) + 22.0*Power(rij,6.0)*Power(xij,6.0)) - 
           7980.0*Power(xii,16.0)*Power(xij,8.0)*
            (15696909.0 - 494343.0*Power(rij,2.0)*Power(xij,2.0) - 
              4182.0*Power(rij,4.0)*Power(xij,4.0) + 34.0*Power(rij,6.0)*Power(xij,6.0)) + 
           2.0*rij*Power(xii,9.0)*Power(xij,16.0)*
            (19433295.0 - 690795.0*Power(rij,2.0)*Power(xij,2.0) + 
              55251.0*Power(rij,4.0)*Power(xij,4.0) + 68.0*Power(rij,6.0)*Power(xij,6.0)) + 
           6.0*rij*Power(xii,25.0)*(-8546580.0 - 
              11329605.0*Power(rij,2.0)*Power(xij,2.0) - 
              24003.0*Power(rij,4.0)*Power(xij,4.0) + 92.0*Power(rij,6.0)*Power(xij,6.0)) - 
           6.0*rij*Power(xii,13.0)*Power(xij,12.0)*
            (-2361196215.0 - 54738810.0*Power(rij,2.0)*Power(xij,2.0) + 
              388626.0*Power(rij,4.0)*Power(xij,4.0) + 92.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 38.0*rij*Power(xii,15.0)*Power(xij,10.0)*
            (808181955.0 - 17168130.0*Power(rij,2.0)*Power(xij,2.0) - 
              32130.0*Power(rij,4.0)*Power(xij,4.0) + 106.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 84.0*Power(xii,10.0)*Power(xij,14.0)*
            (3168630.0 + 683145.0*Power(rij,2.0)*Power(xij,2.0) + 
              54315.0*Power(rij,4.0)*Power(xij,4.0) + 193.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 19.0*rij*Power(xii,17.0)*Power(xij,8.0)*
            (-2525985.0 + 33479460.0*Power(rij,2.0)*Power(xij,2.0) - 
              406980.0*Power(rij,4.0)*Power(xij,4.0) + 272.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 84.0*Power(xii,12.0)*Power(xij,12.0)*
            (-88925130.0 - 19869345.0*Power(rij,2.0)*Power(xij,2.0) - 
              235790.0*Power(rij,4.0)*Power(xij,4.0) + 643.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 210.0*Power(xii,18.0)*Power(xij,6.0)*
            (-496605582.0 + 32638599.0*Power(rij,2.0)*Power(xij,2.0) - 
              564604.0*Power(rij,4.0)*Power(xij,4.0) + 1292.0*Power(rij,6.0)*Power(xij,6.0)\
    ) + 42.0*Power(xii,20.0)*Power(xij,4.0)*
            (-777723210.0 - 46394505.0*Power(rij,2.0)*Power(xij,2.0) + 
              625670.0*Power(rij,4.0)*Power(xij,4.0) + 1292.0*Power(rij,6.0)*Power(xij,6.0)\
    ) + 42.0*Power(xii,24.0)*(-1918620.0 - 11344995.0*Power(rij,2.0)*Power(xij,2.0) - 
              323070.0*Power(rij,4.0)*Power(xij,4.0) + 2114.0*Power(rij,6.0)*Power(xij,6.0)\
    ) - rij*Power(xii,23.0)*Power(xij,2.0)*
            (1919335635.0 + 275096430.0*Power(rij,2.0)*Power(xij,2.0) - 
              3302586.0*Power(rij,4.0)*Power(xij,4.0) + 
              4028.0*Power(rij,6.0)*Power(xij,6.0)) + 
           rij*Power(xii,21.0)*Power(xij,4.0)*
            (-14708379735.0 + 255168270.0*Power(rij,2.0)*Power(xij,2.0) - 
              2899134.0*Power(rij,4.0)*Power(xij,4.0) + 
              5168.0*Power(rij,6.0)*Power(xij,6.0)) - 
           42.0*Power(xii,22.0)*Power(xij,2.0)*
            (81654210.0 + 66273255.0*Power(rij,2.0)*Power(xij,2.0) - 
              1203870.0*Power(rij,4.0)*Power(xij,4.0) + 5206.0*Power(rij,6.0)*Power(xij,6.0)\
    )) - 2.0*Power(E,2.0*rij*xii)*Power(xii,10.0)*
         (21318.0*Power(xii,14.0)*Power(xij,14.0)*
            (-3146850.0 + 4890375.0*rij*xij - 24522750.0*Power(rij,2.0)*Power(xij,2.0) + 
              12162150.0*Power(rij,3.0)*Power(xij,3.0) - 
              1549800.0*Power(rij,4.0)*Power(xij,4.0) - 
              1615950.0*Power(rij,5.0)*Power(xij,5.0) - 
              185220.0*Power(rij,6.0)*Power(xij,6.0) + 
              12240.0*Power(rij,7.0)*Power(xij,7.0) + 
              3960.0*Power(rij,8.0)*Power(xij,8.0) + 300.0*Power(rij,9.0)*Power(xij,9.0) + 
              8.0*Power(rij,10.0)*Power(xij,10.0)) + 
           3.0*Power(xii,24.0)*Power(xij,4.0)*
            (53326350.0 + 97764975.0*rij*xij + 88877250.0*Power(rij,2.0)*Power(xij,2.0) + 
              53326350.0*Power(rij,3.0)*Power(xij,3.0) + 
              23700600.0*Power(rij,4.0)*Power(xij,4.0) + 
              8295210.0*Power(rij,5.0)*Power(xij,5.0) + 
              2370060.0*Power(rij,6.0)*Power(xij,6.0) + 
              564300.0*Power(rij,7.0)*Power(xij,7.0) + 
              112860.0*Power(rij,8.0)*Power(xij,8.0) + 
              22440.0*Power(rij,9.0)*Power(xij,9.0) + 
              1056.0*Power(rij,10.0)*Power(xij,10.0) - 20.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 4.0*Power(xij,28.0)*(13749310575.0 + 13749310575.0*rij*xij + 
              6547290750.0*Power(rij,2.0)*Power(xij,2.0) + 
              1964187225.0*Power(rij,3.0)*Power(xij,3.0) + 
              413513100.0*Power(rij,4.0)*Power(xij,4.0) + 
              64324260.0*Power(rij,5.0)*Power(xij,5.0) + 
              7567560.0*Power(rij,6.0)*Power(xij,6.0) + 
              675675.0*Power(rij,7.0)*Power(xij,7.0) + 
              45045.0*Power(rij,8.0)*Power(xij,8.0) + 
              2145.0*Power(rij,9.0)*Power(xij,9.0) + 66.0*Power(rij,10.0)*Power(xij,10.0) + 
              Power(rij,11.0)*Power(xij,11.0)) - 
           1254.0*Power(xii,16.0)*Power(xij,12.0)*
            (-20241900.0 - 38315025.0*rij*xij - 
              21687750.0*Power(rij,2.0)*Power(xij,2.0) - 
              50122800.0*Power(rij,3.0)*Power(xij,3.0) + 
              14137200.0*Power(rij,4.0)*Power(xij,4.0) - 
              5853330.0*Power(rij,5.0)*Power(xij,5.0) - 
              2687580.0*Power(rij,6.0)*Power(xij,6.0) - 
              208530.0*Power(rij,7.0)*Power(xij,7.0) + 
              19530.0*Power(rij,8.0)*Power(xij,8.0) + 
              3630.0*Power(rij,9.0)*Power(xij,9.0) + 
              172.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) + 
           627.0*Power(xii,12.0)*Power(xij,16.0)*
            (-1240964550.0 + 4740389325.0*rij*xij - 
              3311818650.0*Power(rij,2.0)*Power(xij,2.0) + 
              134804250.0*Power(rij,3.0)*Power(xij,3.0) + 
              407673000.0*Power(rij,4.0)*Power(xij,4.0) + 
              58641030.0*Power(rij,5.0)*Power(xij,5.0) - 
              3549420.0*Power(rij,6.0)*Power(xij,6.0) - 
              1641060.0*Power(rij,7.0)*Power(xij,7.0) - 
              167940.0*Power(rij,8.0)*Power(xij,8.0) - 
              6990.0*Power(rij,9.0)*Power(xij,9.0) - 36.0*Power(rij,10.0)*Power(xij,10.0) + 
              4.0*Power(rij,11.0)*Power(xij,11.0)) + 
           Power(xii,28.0)*(935550.0 + 1715175.0*rij*xij + 
              1559250.0*Power(rij,2.0)*Power(xij,2.0) + 
              935550.0*Power(rij,3.0)*Power(xij,3.0) + 
              415800.0*Power(rij,4.0)*Power(xij,4.0) + 
              145530.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              9900.0*Power(rij,7.0)*Power(xij,7.0) + 1980.0*Power(rij,8.0)*Power(xij,8.0) + 
              330.0*Power(rij,9.0)*Power(xij,9.0) + 44.0*Power(rij,10.0)*Power(xij,10.0) + 
              4.0*Power(rij,11.0)*Power(xij,11.0)) + 
           2.0*Power(xii,2.0)*Power(xij,26.0)*
            (-937068397650.0 - 815439881025.0*rij*xij - 
              332904552750.0*Power(rij,2.0)*Power(xij,2.0) - 
              84006776700.0*Power(rij,3.0)*Power(xij,3.0) - 
              14504767200.0*Power(rij,4.0)*Power(xij,4.0) - 
              1786235220.0*Power(rij,5.0)*Power(xij,5.0) - 
              157754520.0*Power(rij,6.0)*Power(xij,6.0) - 
              9667350.0*Power(rij,7.0)*Power(xij,7.0) - 
              367290.0*Power(rij,8.0)*Power(xij,8.0) - 
              5115.0*Power(rij,9.0)*Power(xij,9.0) + 
              198.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           6.0*Power(xii,4.0)*Power(xij,24.0)*
            (-2262441500550.0 - 1503711230175.0*rij*xij - 
              426178264050.0*Power(rij,2.0)*Power(xij,2.0) - 
              60134347350.0*Power(rij,3.0)*Power(xij,3.0) - 
              2014551000.0*Power(rij,4.0)*Power(xij,4.0) + 
              846111420.0*Power(rij,5.0)*Power(xij,5.0) + 
              184864680.0*Power(rij,6.0)*Power(xij,6.0) + 
              20183130.0*Power(rij,7.0)*Power(xij,7.0) + 
              1367190.0*Power(rij,8.0)*Power(xij,8.0) + 
              57255.0*Power(rij,9.0)*Power(xij,9.0) + 
              1298.0*Power(rij,10.0)*Power(xij,10.0) + 10.0*Power(rij,11.0)*Power(xij,11.0)) \
    - Power(xii,26.0)*Power(xij,2.0)*(17775450.0 + 32588325.0*rij*xij + 
              29625750.0*Power(rij,2.0)*Power(xij,2.0) + 
              17775450.0*Power(rij,3.0)*Power(xij,3.0) + 
              7900200.0*Power(rij,4.0)*Power(xij,4.0) + 
              2765070.0*Power(rij,5.0)*Power(xij,5.0) + 
              790020.0*Power(rij,6.0)*Power(xij,6.0) + 
              188100.0*Power(rij,7.0)*Power(xij,7.0) + 
              37620.0*Power(rij,8.0)*Power(xij,8.0) + 
              6270.0*Power(rij,9.0)*Power(xij,9.0) + 
              836.0*Power(rij,10.0)*Power(xij,10.0) + 16.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 9.0*Power(xii,22.0)*Power(xij,6.0)*
            (-100727550.0 - 184667175.0*rij*xij - 
              167879250.0*Power(rij,2.0)*Power(xij,2.0) - 
              100727550.0*Power(rij,3.0)*Power(xij,3.0) - 
              44767800.0*Power(rij,4.0)*Power(xij,4.0) - 
              15668730.0*Power(rij,5.0)*Power(xij,5.0) - 
              4476780.0*Power(rij,6.0)*Power(xij,6.0) - 
              971520.0*Power(rij,7.0)*Power(xij,7.0) - 
              307560.0*Power(rij,8.0)*Power(xij,8.0) - 
              27060.0*Power(rij,9.0)*Power(xij,9.0) + 
              264.0*Power(rij,10.0)*Power(xij,10.0) + 64.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 9.0*Power(xii,6.0)*Power(xij,22.0)*
            (3452543428950.0 + 1097992509075.0*rij*xij - 
              101420792550.0*Power(rij,2.0)*Power(xij,2.0) - 
              110557373850.0*Power(rij,3.0)*Power(xij,3.0) - 
              24909330600.0*Power(rij,4.0)*Power(xij,4.0) - 
              2686726350.0*Power(rij,5.0)*Power(xij,5.0) - 
              93485700.0*Power(rij,6.0)*Power(xij,6.0) + 
              12941280.0*Power(rij,7.0)*Power(xij,7.0) + 
              2081640.0*Power(rij,8.0)*Power(xij,8.0) + 
              137940.0*Power(rij,9.0)*Power(xij,9.0) + 
              4664.0*Power(rij,10.0)*Power(xij,10.0) + 64.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 22.0*Power(xii,20.0)*Power(xij,8.0)*
            (-164826900.0 - 302182650.0*rij*xij - 
              274711500.0*Power(rij,2.0)*Power(xij,2.0) - 
              164826900.0*Power(rij,3.0)*Power(xij,3.0) - 
              73256400.0*Power(rij,4.0)*Power(xij,4.0) - 
              26991090.0*Power(rij,5.0)*Power(xij,5.0) - 
              4622940.0*Power(rij,6.0)*Power(xij,6.0) - 
              2941110.0*Power(rij,7.0)*Power(xij,7.0) - 
              438930.0*Power(rij,8.0)*Power(xij,8.0) - 
              5505.0*Power(rij,9.0)*Power(xij,9.0) + 
              2082.0*Power(rij,10.0)*Power(xij,10.0) + 82.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 22.0*Power(xii,18.0)*Power(xij,10.0)*
            (-494480700.0 - 906547950.0*rij*xij - 
              824134500.0*Power(rij,2.0)*Power(xij,2.0) - 
              475684650.0*Power(rij,3.0)*Power(xij,3.0) - 
              294953400.0*Power(rij,4.0)*Power(xij,4.0) + 
              2663010.0*Power(rij,5.0)*Power(xij,5.0) - 
              40797540.0*Power(rij,6.0)*Power(xij,6.0) - 
              10248390.0*Power(rij,7.0)*Power(xij,7.0) - 
              434610.0*Power(rij,8.0)*Power(xij,8.0) + 
              65865.0*Power(rij,9.0)*Power(xij,9.0) + 
              6366.0*Power(rij,10.0)*Power(xij,10.0) + 136.0*Power(rij,11.0)*Power(xij,11.0)\
    ) + 11.0*Power(xii,8.0)*Power(xij,20.0)*
            (-2338604626050.0 + 656001834075.0*rij*xij + 
              504510561450.0*Power(rij,2.0)*Power(xij,2.0) + 
              51560967150.0*Power(rij,3.0)*Power(xij,3.0) - 
              15574998600.0*Power(rij,4.0)*Power(xij,4.0) - 
              5055778350.0*Power(rij,5.0)*Power(xij,5.0) - 
              626213700.0*Power(rij,6.0)*Power(xij,6.0) - 
              34768620.0*Power(rij,7.0)*Power(xij,7.0) + 
              207540.0*Power(rij,8.0)*Power(xij,8.0) + 
              150240.0*Power(rij,9.0)*Power(xij,9.0) + 
              8464.0*Power(rij,10.0)*Power(xij,10.0) + 164.0*Power(rij,11.0)*Power(xij,11.0)\
    ) - 11.0*Power(xii,10.0)*Power(xij,18.0)*
            (742805182350.0 - 933111659025.0*rij*xij + 
              57080542050.0*Power(rij,2.0)*Power(xij,2.0) + 
              129505209750.0*Power(rij,3.0)*Power(xij,3.0) + 
              19066887000.0*Power(rij,4.0)*Power(xij,4.0) - 
              1817573310.0*Power(rij,5.0)*Power(xij,5.0) - 
              810647460.0*Power(rij,6.0)*Power(xij,6.0) - 
              97669980.0*Power(rij,7.0)*Power(xij,7.0) - 
              5173020.0*Power(rij,8.0)*Power(xij,8.0) - 
              37770.0*Power(rij,9.0)*Power(xij,9.0) + 
              8212.0*Power(rij,10.0)*Power(xij,10.0) + 272.0*Power(rij,11.0)*Power(xij,11.0))\
    ))/(1.8711e6*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),19.0))
    ;
  }
  return S;
}

double Slater_4S_1S(double rij,double xii,double xij)
{
  return Slater_1S_4S(rij,xij,xii);
}

double Slater_4S_2S(double rij,double xii,double xij)
{
  return Slater_2S_4S(rij,xij,xii);
}

double Slater_4S_3S(double rij,double xii,double xij)
{
  return Slater_3S_4S(rij,xij,xii);
}

static double Slater_5S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-12164510040883200000.0 + 12164510040883200000.0*Power(E,2.0*rij*xii) - 
        22324788235240115625.0*rij*xii - 
        20320556388713831250.0*Power(rij,2.0)*Power(xii,2.0) - 
        12225924086428552500.0*Power(rij,3.0)*Power(xii,3.0) - 
        5467446348494130000.0*Power(rij,4.0)*Power(xii,4.0) - 
        1937619942864606000.0*Power(rij,5.0)*Power(xii,5.0) - 
        566528792821992000.0*Power(rij,6.0)*Power(xii,6.0) - 
        140462831126217600.0*Power(rij,7.0)*Power(xii,7.0) - 
        30115609927603200.0*Power(rij,8.0)*Power(xii,8.0) - 
        5663731244371200.0*Power(rij,9.0)*Power(xii,9.0) - 
        943983142502400.0*Power(rij,10.0)*Power(xii,10.0) - 
        140427244339200.0*Power(rij,11.0)*Power(xii,11.0) - 
        18723632578560.0*Power(rij,12.0)*Power(xii,12.0) - 
        2240434667520.0*Power(rij,13.0)*Power(xii,13.0) - 
        240046571520.0*Power(rij,14.0)*Power(xii,14.0) - 
        22861578240.0*Power(rij,15.0)*Power(xii,15.0) - 
        1905131520.0*Power(rij,16.0)*Power(xii,16.0) - 
        134479872.0*Power(rij,17.0)*Power(xii,17.0) - 
        7471104.0*Power(rij,18.0)*Power(xii,18.0) - 262144.0*Power(rij,19.0)*Power(xii,19.0))/
      (1.21645100408832e19*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (70875.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),19.0) + 
        Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-630.0*Power(rij,8.0)*Power(xii,34.0) - 10.0*Power(rij,9.0)*Power(xii,35.0) + 
           70875.0*Power(xij,26.0) + 127575.0*rij*xii*Power(xij,26.0) - 
           30.0*Power(rij,7.0)*Power(xii,33.0)*(630.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           14175.0*Power(xii,2.0)*Power(xij,24.0)*
            (-95.0 + 8.0*Power(rij,2.0)*Power(xij,2.0)) + 
           4725.0*rij*Power(xii,3.0)*Power(xij,24.0)*
            (-513.0 + 14.0*Power(rij,2.0)*Power(xij,2.0)) - 
           90.0*Power(rij,6.0)*Power(xii,32.0)*
            (3920.0 + 43.0*Power(rij,2.0)*Power(xij,2.0)) + 
           4725.0*rij*Power(xii,5.0)*Power(xij,22.0)*
            (4617.0 - 266.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           14175.0*Power(xii,4.0)*Power(xij,22.0)*
            (855.0 - 152.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 36.0*Power(rij,5.0)*Power(xii,31.0)*
            (-124950.0 - 4985.0*Power(rij,2.0)*Power(xij,2.0) + 
              13.0*Power(rij,4.0)*Power(xij,4.0)) + 
           36.0*Power(rij,4.0)*Power(xii,30.0)*
            (-1124550.0 - 127960.0*Power(rij,2.0)*Power(xij,2.0) + 
              863.0*Power(rij,4.0)*Power(xij,4.0)) + 
           135.0*rij*Power(xii,7.0)*Power(xij,20.0)*
            (-915705.0 + 83790.0*Power(rij,2.0)*Power(xij,2.0) - 
              1330.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           315.0*Power(xii,6.0)*Power(xij,20.0)*
            (-218025.0 + 61560.0*Power(rij,2.0)*Power(xij,2.0) - 
              1710.0*Power(rij,4.0)*Power(xij,4.0) + 8.0*Power(rij,6.0)*Power(xij,6.0)) - 
           36.0*Power(rij,3.0)*Power(xii,29.0)*
            (7122150.0 + 2102730.0*Power(rij,2.0)*Power(xij,2.0) - 
              23294.0*Power(rij,4.0)*Power(xij,4.0) + 37.0*Power(rij,6.0)*Power(xij,6.0)) - 
           36.0*Power(rij,2.0)*Power(xii,28.0)*
            (30523500.0 + 23401350.0*Power(rij,2.0)*Power(xij,2.0) - 
              299250.0*Power(rij,4.0)*Power(xij,4.0) + 1297.0*Power(rij,6.0)*Power(xij,6.0)\
    ) + rij*Power(xii,17.0)*Power(xij,10.0)*
            (1073961177975.0 - 21753487980.0*Power(rij,2.0)*Power(xij,2.0) - 
              745994340.0*Power(rij,4.0)*Power(xij,4.0) + 
              5307156.0*Power(rij,6.0)*Power(xij,6.0) - 818.0*Power(rij,8.0)*Power(xij,8.0)\
    ) + 10.0*rij*Power(xii,9.0)*Power(xij,18.0)*
            (49448070.0 - 6409935.0*Power(rij,2.0)*Power(xij,2.0) + 
              161595.0*Power(rij,4.0)*Power(xij,4.0) - 
              1026.0*Power(rij,6.0)*Power(xij,6.0) + Power(rij,8.0)*Power(xij,8.0)) + 
           90.0*Power(xii,8.0)*Power(xij,18.0)*
            (3052350.0 - 1220940.0*Power(rij,2.0)*Power(xij,2.0) + 
              53865.0*Power(rij,4.0)*Power(xij,4.0) - 
              532.0*Power(rij,6.0)*Power(xij,6.0) + Power(rij,8.0)*Power(xij,8.0)) - 
           1710.0*Power(xii,10.0)*Power(xij,16.0)*
            (481950.0 - 257040.0*Power(rij,2.0)*Power(xij,2.0) + 
              16065.0*Power(rij,4.0)*Power(xij,4.0) - 
              252.0*Power(rij,6.0)*Power(xij,6.0) + Power(rij,8.0)*Power(xij,8.0)) + 
           6.0*rij*Power(xii,11.0)*Power(xij,16.0)*
            (-207559800.0 + 50390550.0*Power(rij,2.0)*Power(xij,2.0) - 
              1165815.0*Power(rij,4.0)*Power(xij,4.0) + 
              21396.0*Power(rij,6.0)*Power(xij,6.0) + 5.0*Power(rij,8.0)*Power(xij,8.0)) - 
           18.0*rij*Power(xii,13.0)*Power(xij,14.0)*
            (-1703720025.0 - 155669850.0*Power(rij,2.0)*Power(xij,2.0) - 
              7410270.0*Power(rij,4.0)*Power(xij,4.0) - 
              1532.0*Power(rij,6.0)*Power(xij,6.0) + 26.0*Power(rij,8.0)*Power(xij,8.0)) + 
           18.0*rij*Power(xii,15.0)*Power(xij,12.0)*
            (19380896325.0 + 1329128850.0*Power(rij,2.0)*Power(xij,2.0) - 
              7608930.0*Power(rij,4.0)*Power(xij,4.0) - 
              116238.0*Power(rij,6.0)*Power(xij,6.0) + 74.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 18.0*Power(xii,12.0)*Power(xij,14.0)*
            (89026875.0 + 179071200.0*Power(rij,2.0)*Power(xij,2.0) + 
              1552950.0*Power(rij,4.0)*Power(xij,4.0) + 
              295820.0*Power(rij,6.0)*Power(xij,6.0) + 146.0*Power(rij,8.0)*Power(xij,8.0)) \
    + 18.0*rij*Power(xii,25.0)*Power(xij,2.0)*
            (-5449970925.0 - 1137574935.0*Power(rij,2.0)*Power(xij,2.0) + 
              37834755.0*Power(rij,4.0)*Power(xij,4.0) - 
              273062.0*Power(rij,6.0)*Power(xij,6.0) + 171.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 9.0*rij*Power(xii,19.0)*Power(xij,8.0)*
            (-37914907275.0 + 7613889570.0*Power(rij,2.0)*Power(xij,2.0) - 
              170524620.0*Power(rij,4.0)*Power(xij,4.0) + 
              397936.0*Power(rij,6.0)*Power(xij,6.0) + 342.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 3.0*rij*Power(xii,23.0)*Power(xij,4.0)*
            (219130630425.0 - 11118046590.0*Power(rij,2.0)*Power(xij,2.0) + 
              327611970.0*Power(rij,4.0)*Power(xij,4.0) - 
              2920908.0*Power(rij,6.0)*Power(xij,6.0) + 
              2584.0*Power(rij,8.0)*Power(xij,8.0)) + 
           3.0*rij*Power(xii,21.0)*Power(xij,6.0)*
            (-345162539925.0 + 19030764690.0*Power(rij,2.0)*Power(xij,2.0) - 
              141976170.0*Power(rij,4.0)*Power(xij,4.0) - 
              1441872.0*Power(rij,6.0)*Power(xij,6.0) + 
              2584.0*Power(rij,8.0)*Power(xij,8.0)) + 
           63.0*Power(xii,20.0)*Power(xij,6.0)*
            (-50980542525.0 + 6240202920.0*Power(rij,2.0)*Power(xij,2.0) - 
              201314310.0*Power(rij,4.0)*Power(xij,4.0) + 
              956080.0*Power(rij,6.0)*Power(xij,6.0) + 2584.0*Power(rij,8.0)*Power(xij,8.0)\
    ) + 18.0*Power(xii,14.0)*Power(xij,12.0)*
            (-7803332775.0 - 2519206200.0*Power(rij,2.0)*Power(xij,2.0) - 
              119719950.0*Power(rij,4.0)*Power(xij,4.0) + 
              182280.0*Power(rij,6.0)*Power(xij,6.0) + 2734.0*Power(rij,8.0)*Power(xij,8.0)\
    ) - 18.0*Power(xii,26.0)*(195859125.0 + 1794781800.0*Power(rij,2.0)*Power(xij,2.0) + 
              67337235.0*Power(rij,4.0)*Power(xij,4.0) - 
              1659700.0*Power(rij,6.0)*Power(xij,6.0) + 
              4089.0*Power(rij,8.0)*Power(xij,8.0)) + 
           9.0*Power(xii,18.0)*Power(xij,8.0)*
            (-357591274425.0 + 8328390840.0*Power(rij,2.0)*Power(xij,2.0) + 
              912042180.0*Power(rij,4.0)*Power(xij,4.0) - 
              12842480.0*Power(rij,6.0)*Power(xij,6.0) + 
              10678.0*Power(rij,8.0)*Power(xij,8.0)) - 
           9.0*Power(xii,16.0)*Power(xij,10.0)*
            (128599724925.0 + 21298077360.0*Power(rij,2.0)*Power(xij,2.0) - 
              267928500.0*Power(rij,4.0)*Power(xij,4.0) - 
              5458320.0*Power(rij,6.0)*Power(xij,6.0) + 
              14722.0*Power(rij,8.0)*Power(xij,8.0)) + 
           18.0*Power(xii,24.0)*Power(xij,2.0)*
            (-7604930025.0 - 8866107180.0*Power(rij,2.0)*Power(xij,2.0) + 
              399272265.0*Power(rij,4.0)*Power(xij,4.0) - 
              5925780.0*Power(rij,6.0)*Power(xij,6.0) + 
              17651.0*Power(rij,8.0)*Power(xij,8.0)) - 
           9.0*Power(xii,22.0)*Power(xij,4.0)*
            (129194933175.0 + 3909863160.0*Power(rij,2.0)*Power(xij,2.0) + 
              91420770.0*Power(rij,4.0)*Power(xij,4.0) - 
              8762040.0*Power(rij,6.0)*Power(xij,6.0) + 
              43928.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xii,27.0)*(-2884470750.0*rij - 
              6409935000.0*Power(rij,3.0)*Power(xij,2.0) + 
              28332990.0*Power(rij,5.0)*Power(xij,4.0) + 
              58104.0*Power(rij,7.0)*Power(xij,6.0) + 818.0*Power(rij,9.0)*Power(xij,8.0))) \
    + Power(E,2.0*rij*xii)*Power(xii,12.0)*
         (Power(xii,8.0)*Power(xij,18.0)*
            (3218321469825.0 - 341234165475.0*rij*xij - 
              393132783960.0*Power(rij,2.0)*Power(xij,2.0) - 
              57092294070.0*Power(rij,3.0)*Power(xij,3.0) + 
              822786930.0*Power(rij,4.0)*Power(xij,4.0) + 
              982835910.0*Power(rij,5.0)*Power(xij,5.0) + 
              106664040.0*Power(rij,6.0)*Power(xij,6.0) + 
              4915116.0*Power(rij,7.0)*Power(xij,7.0) + 
              73602.0*Power(rij,8.0)*Power(xij,8.0) - 818.0*Power(rij,9.0)*Power(xij,9.0)) + 
           10.0*Power(xij,26.0)*(352546425.0 + 288447075.0*rij*xij + 
              109884600.0*Power(rij,2.0)*Power(xij,2.0) + 
              25639740.0*Power(rij,3.0)*Power(xij,3.0) + 
              4048380.0*Power(rij,4.0)*Power(xij,4.0) + 
              449820.0*Power(rij,5.0)*Power(xij,5.0) + 
              35280.0*Power(rij,6.0)*Power(xij,6.0) + 
              1890.0*Power(rij,7.0)*Power(xij,7.0) + 63.0*Power(rij,8.0)*Power(xij,8.0) + 
              Power(rij,9.0)*Power(xij,9.0)) + 
           30.0*Power(xii,2.0)*Power(xij,24.0)*
            (4562958015.0 + 3269982555.0*rij*xij + 
              1076869080.0*Power(rij,2.0)*Power(xij,2.0) + 
              213664500.0*Power(rij,3.0)*Power(xij,3.0) + 
              28081620.0*Power(rij,4.0)*Power(xij,4.0) + 
              2523276.0*Power(rij,5.0)*Power(xij,5.0) + 
              153552.0*Power(rij,6.0)*Power(xij,6.0) + 
              5982.0*Power(rij,7.0)*Power(xij,7.0) + 129.0*Power(rij,8.0)*Power(xij,8.0) + 
              Power(rij,9.0)*Power(xij,9.0)) - 
           15.0*Power(xii,24.0)*Power(xij,2.0)*
            (-89775.0 - 161595.0*rij*xij - 143640.0*Power(rij,2.0)*Power(xij,2.0) - 
              83790.0*Power(rij,3.0)*Power(xij,3.0) - 
              35910.0*Power(rij,4.0)*Power(xij,4.0) - 
              11970.0*Power(rij,5.0)*Power(xij,5.0) - 
              3192.0*Power(rij,6.0)*Power(xij,6.0) - 684.0*Power(rij,7.0)*Power(xij,7.0) - 
              114.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) - 
           5.0*Power(xii,26.0)*(14175.0 + 25515.0*rij*xij + 
              22680.0*Power(rij,2.0)*Power(xij,2.0) + 
              13230.0*Power(rij,3.0)*Power(xij,3.0) + 
              5670.0*Power(rij,4.0)*Power(xij,4.0) + 1890.0*Power(rij,5.0)*Power(xij,5.0) + 
              504.0*Power(rij,6.0)*Power(xij,6.0) + 108.0*Power(rij,7.0)*Power(xij,7.0) + 
              18.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,9.0)*Power(xij,9.0)) - 
           1938.0*Power(xii,14.0)*Power(xij,12.0)*
            (-826875.0 + 15824025.0*rij*xij - 23398200.0*Power(rij,2.0)*Power(xij,2.0) + 
              12344850.0*Power(rij,3.0)*Power(xij,3.0) + 
              1244250.0*Power(rij,4.0)*Power(xij,4.0) - 
              384930.0*Power(rij,5.0)*Power(xij,5.0) - 
              59640.0*Power(rij,6.0)*Power(xij,6.0) - 
              1848.0*Power(rij,7.0)*Power(xij,7.0) + 84.0*Power(rij,8.0)*Power(xij,8.0) + 
              4.0*Power(rij,9.0)*Power(xij,9.0)) + 
           1938.0*Power(xii,12.0)*Power(xij,14.0)*
            (72476775.0 - 180008325.0*rij*xij + 
              98907480.0*Power(rij,2.0)*Power(xij,2.0) + 
              11224710.0*Power(rij,3.0)*Power(xij,3.0) - 
              4235490.0*Power(rij,4.0)*Power(xij,4.0) - 
              791910.0*Power(rij,5.0)*Power(xij,5.0) - 
              31080.0*Power(rij,6.0)*Power(xij,6.0) + 
              2232.0*Power(rij,7.0)*Power(xij,7.0) + 204.0*Power(rij,8.0)*Power(xij,8.0) + 
              4.0*Power(rij,9.0)*Power(xij,9.0)) + 
           342.0*Power(xii,16.0)*Power(xij,10.0)*
            (2409750.0 + 3641400.0*rij*xij + 9424800.0*Power(rij,2.0)*Power(xij,2.0) - 
              8193150.0*Power(rij,3.0)*Power(xij,3.0) + 
              6301050.0*Power(rij,4.0)*Power(xij,4.0) + 
              400470.0*Power(rij,5.0)*Power(xij,5.0) - 
              143640.0*Power(rij,6.0)*Power(xij,6.0) - 
              15518.0*Power(rij,7.0)*Power(xij,7.0) - 281.0*Power(rij,8.0)*Power(xij,8.0) + 
              9.0*Power(rij,9.0)*Power(xij,9.0)) - 
           171.0*Power(xii,10.0)*Power(xij,16.0)*
            (-6768406575.0 + 6280474725.0*rij*xij + 
              438336360.0*Power(rij,2.0)*Power(xij,2.0) - 
              400731030.0*Power(rij,3.0)*Power(xij,3.0) - 
              74168430.0*Power(rij,4.0)*Power(xij,4.0) - 
              2490810.0*Power(rij,5.0)*Power(xij,5.0) + 
              461160.0*Power(rij,6.0)*Power(xij,6.0) + 
              51244.0*Power(rij,7.0)*Power(xij,7.0) + 
              1858.0*Power(rij,8.0)*Power(xij,8.0) + 18.0*Power(rij,9.0)*Power(xij,9.0)) + 
           9.0*Power(xii,22.0)*Power(xij,4.0)*
            (-1346625.0 - 2423925.0*rij*xij - 2154600.0*Power(rij,2.0)*Power(xij,2.0) - 
              1256850.0*Power(rij,3.0)*Power(xij,3.0) - 
              538650.0*Power(rij,4.0)*Power(xij,4.0) - 
              179550.0*Power(rij,5.0)*Power(xij,5.0) - 
              47880.0*Power(rij,6.0)*Power(xij,6.0) - 
              14264.0*Power(rij,7.0)*Power(xij,7.0) + 292.0*Power(rij,8.0)*Power(xij,8.0) + 
              52.0*Power(rij,9.0)*Power(xij,9.0)) - 
           9.0*Power(xii,4.0)*Power(xij,22.0)*
            (-129194933175.0 - 73043543475.0*rij*xij - 
              17732214360.0*Power(rij,2.0)*Power(xij,2.0) - 
              2275149870.0*Power(rij,3.0)*Power(xij,3.0) - 
              134674470.0*Power(rij,4.0)*Power(xij,4.0) + 
              3148110.0*Power(rij,5.0)*Power(xij,5.0) + 
              1197000.0*Power(rij,6.0)*Power(xij,6.0) + 
              93176.0*Power(rij,7.0)*Power(xij,7.0) + 
              3452.0*Power(rij,8.0)*Power(xij,8.0) + 52.0*Power(rij,9.0)*Power(xij,9.0)) + 
           9.0*Power(xii,6.0)*Power(xij,20.0)*
            (356863797675.0 + 115054179975.0*rij*xij + 
              3909863160.0*Power(rij,2.0)*Power(xij,2.0) - 
              3706015530.0*Power(rij,3.0)*Power(xij,3.0) - 
              798544530.0*Power(rij,4.0)*Power(xij,4.0) - 
              75669510.0*Power(rij,5.0)*Power(xij,5.0) - 
              3319400.0*Power(rij,6.0)*Power(xij,6.0) - 
              6456.0*Power(rij,7.0)*Power(xij,7.0) + 5188.0*Power(rij,8.0)*Power(xij,8.0) + 
              148.0*Power(rij,9.0)*Power(xij,9.0)) - 
           9.0*Power(xii,20.0)*Power(xij,6.0)*
            (-7630875.0 - 13735575.0*rij*xij - 12209400.0*Power(rij,2.0)*Power(xij,2.0) - 
              7122150.0*Power(rij,3.0)*Power(xij,3.0) - 
              3052350.0*Power(rij,4.0)*Power(xij,4.0) - 
              777210.0*Power(rij,5.0)*Power(xij,5.0) - 
              591640.0*Power(rij,6.0)*Power(xij,6.0) + 
              3064.0*Power(rij,7.0)*Power(xij,7.0) + 5468.0*Power(rij,8.0)*Power(xij,8.0) + 
              148.0*Power(rij,9.0)*Power(xij,9.0)) + 
           2.0*Power(xii,18.0)*Power(xij,8.0)*
            (-137355750.0 - 247240350.0*rij*xij - 
              219769200.0*Power(rij,2.0)*Power(xij,2.0) - 
              151171650.0*Power(rij,3.0)*Power(xij,3.0) + 
              13976550.0*Power(rij,4.0)*Power(xij,4.0) - 
              66692430.0*Power(rij,5.0)*Power(xij,5.0) - 
              1640520.0*Power(rij,6.0)*Power(xij,6.0) + 
              1046142.0*Power(rij,7.0)*Power(xij,7.0) + 
              66249.0*Power(rij,8.0)*Power(xij,8.0) + 409.0*Power(rij,9.0)*Power(xij,9.0))))/
      (70875.*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),19.0))
    ;
  }
  return S;
}

static double Slater_5S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-1532728265151283200000.0 + 1532728265151283200000.0*Power(E,2.0*rij*xii) - 
        2836013677800293615625.0*rij*xii - 
        2606570825298020831250.0*Power(rij,2.0)*Power(xii,2.0) - 
        1586115932378174071875.0*Power(rij,3.0)*Power(xii,3.0) - 
        718622941126509168750.0*Power(rij,4.0)*Power(xii,4.0) - 
        258482050835109601500.0*Power(rij,5.0)*Power(xii,5.0) - 
        76853380678272198000.0*Power(rij,6.0)*Power(xii,6.0) - 
        19417985233400754000.0*Power(rij,7.0)*Power(xii,7.0) - 
        4253183134704504000.0*Power(rij,8.0)*Power(xii,8.0) - 
        819670099680432000.0*Power(rij,9.0)*Power(xii,9.0) - 
        140553592289510400.0*Power(rij,10.0)*Power(xii,10.0) - 
        21625475644281600.0*Power(rij,11.0)*Power(xii,11.0) - 
        3003582726144000.0*Power(rij,12.0)*Power(xii,12.0) - 
        378073350144000.0*Power(rij,13.0)*Power(xii,13.0) - 
        43208382873600.0*Power(rij,14.0)*Power(xii,14.0) - 
        4480869335040.0*Power(rij,15.0)*Power(xii,15.0) - 
        420081500160.0*Power(rij,16.0)*Power(xii,16.0) - 
        35300966400.0*Power(rij,17.0)*Power(xii,17.0) - 
        2614886400.0*Power(rij,18.0)*Power(xii,18.0) - 
        165150720.0*Power(rij,19.0)*Power(xii,19.0) - 
        8257536.0*Power(rij,20.0)*Power(xii,20.0) - 262144.0*Power(rij,21.0)*Power(xii,21.0))/
      (1.5327282651512832e21*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (4677750.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),21.0) + 
        110.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-432.0*Power(rij,8.0)*Power(xii,36.0) - 6.0*Power(rij,9.0)*Power(xii,37.0) + 
           42525.0*Power(xij,28.0) + 76545.0*rij*xii*Power(xij,28.0) + 
           19845.0*rij*Power(xii,3.0)*Power(xij,26.0)*
            (-81.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           1134.0*Power(rij,6.0)*Power(xii,34.0)*
            (272.0 + 5.0*Power(rij,2.0)*Power(xij,2.0)) - 
           8.0*Power(rij,7.0)*Power(xii,35.0)*(1836.0 + 7.0*Power(rij,2.0)*Power(xij,2.0)) + 
           8505.0*Power(xii,2.0)*Power(xij,26.0)*
            (-105.0 + 8.0*Power(rij,2.0)*Power(xij,2.0)) + 
           378.0*Power(rij,5.0)*Power(xii,33.0)*
            (-11628.0 - 666.0*Power(rij,2.0)*Power(xij,2.0) + 
              Power(rij,4.0)*Power(xij,4.0)) + 
           5670.0*rij*Power(xii,5.0)*Power(xij,24.0)*
            (2835.0 - 147.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) \
    + 17010.0*Power(xii,4.0)*Power(xij,24.0)*
            (525.0 - 84.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           378.0*Power(rij,4.0)*Power(xii,32.0)*
            (-116280.0 - 17444.0*Power(rij,2.0)*Power(xij,2.0) + 
              59.0*Power(rij,4.0)*Power(xij,4.0)) + 
           162.0*rij*Power(xii,7.0)*Power(xij,22.0)*
            (-628425.0 + 51450.0*Power(rij,2.0)*Power(xij,2.0) - 
              735.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           378.0*Power(xii,6.0)*Power(xij,22.0)*
            (-149625.0 + 37800.0*Power(rij,2.0)*Power(xij,2.0) - 
              945.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           18.0*Power(rij,3.0)*Power(xii,31.0)*
            (17093160.0 + 6309387.0*Power(rij,2.0)*Power(xij,2.0) - 
              23562.0*Power(rij,4.0)*Power(xij,4.0) + 16.0*Power(rij,6.0)*Power(xij,6.0)) + 
           54.0*Power(rij,2.0)*Power(xii,30.0)*
            (-26860680.0 - 24843735.0*Power(rij,2.0)*Power(xij,2.0) - 
              40180.0*Power(rij,4.0)*Power(xij,4.0) + 578.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 378.0*rij*Power(xii,23.0)*Power(xij,6.0)*
            (-14625683325.0 + 704051250.0*Power(rij,2.0)*Power(xij,2.0) - 
              10752861.0*Power(rij,4.0)*Power(xij,4.0) + 
              33478.0*Power(rij,6.0)*Power(xij,6.0)) + 
           3.0*rij*Power(xii,9.0)*Power(xij,20.0)*
            (152707275.0 - 17595900.0*Power(rij,2.0)*Power(xij,2.0) + 
              396900.0*Power(rij,4.0)*Power(xij,4.0) - 
              2268.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           27.0*Power(xii,8.0)*Power(xij,20.0)*
            (9426375.0 - 3351600.0*Power(rij,2.0)*Power(xij,2.0) + 
              132300.0*Power(rij,4.0)*Power(xij,4.0) - 
              1176.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           567.0*Power(xii,10.0)*Power(xij,18.0)*
            (1526175.0 - 718200.0*Power(rij,2.0)*Power(xij,2.0) + 
              39900.0*Power(rij,4.0)*Power(xij,4.0) - 
              560.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           54.0*rij*Power(xii,13.0)*Power(xij,16.0)*
            (-1356769575.0 - 127011675.0*Power(rij,2.0)*Power(xij,2.0) - 
              3867843.0*Power(rij,4.0)*Power(xij,4.0) - 
              8556.0*Power(rij,6.0)*Power(xij,6.0) + 7.0*Power(rij,8.0)*Power(xij,8.0)) + 
           7.0*rij*Power(xii,11.0)*Power(xij,18.0)*
            (-151091325.0 + 45272250.0*Power(rij,2.0)*Power(xij,2.0) - 
              647676.0*Power(rij,4.0)*Power(xij,4.0) + 
              15336.0*Power(rij,6.0)*Power(xij,6.0) + 8.0*Power(rij,8.0)*Power(xij,8.0)) + 
           18.0*rij*Power(xii,15.0)*Power(xij,14.0)*
            (63046289250.0 + 3917182500.0*Power(rij,2.0)*Power(xij,2.0) + 
              10158435.0*Power(rij,4.0)*Power(xij,4.0) - 
              178842.0*Power(rij,6.0)*Power(xij,6.0) + 16.0*Power(rij,8.0)*Power(xij,8.0)) \
    + 378.0*rij*Power(xii,21.0)*Power(xij,8.0)*
            (-8559820125.0 + 17573325.0*Power(rij,2.0)*Power(xij,2.0) + 
              7421001.0*Power(rij,4.0)*Power(xij,4.0) - 
              49096.0*Power(rij,6.0)*Power(xij,6.0) + 19.0*Power(rij,8.0)*Power(xij,8.0)) - 
           378.0*Power(xii,12.0)*Power(xij,16.0)*
            (17296650.0 + 14244300.0*Power(rij,2.0)*Power(xij,2.0) + 
              360525.0*Power(rij,4.0)*Power(xij,4.0) + 
              15928.0*Power(rij,6.0)*Power(xij,6.0) + 22.0*Power(rij,8.0)*Power(xij,8.0)) - 
           189.0*rij*Power(xii,25.0)*Power(xij,4.0)*
            (9994948425.0 + 63821700.0*Power(rij,2.0)*Power(xij,2.0) - 
              1458540.0*Power(rij,4.0)*Power(xij,4.0) - 
              18756.0*Power(rij,6.0)*Power(xij,6.0) + 38.0*Power(rij,8.0)*Power(xij,8.0)) - 
           189.0*Power(xii,24.0)*Power(xij,4.0)*
            (17962854525.0 + 4036942800.0*Power(rij,2.0)*Power(xij,2.0) - 
              126472500.0*Power(rij,4.0)*Power(xij,4.0) + 
              765464.0*Power(rij,6.0)*Power(xij,6.0) + 190.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 21.0*rij*Power(xii,19.0)*Power(xij,10.0)*
            (-228066210225.0 + 13487616450.0*Power(rij,2.0)*Power(xij,2.0) - 
              85465800.0*Power(rij,4.0)*Power(xij,4.0) - 
              320112.0*Power(rij,6.0)*Power(xij,6.0) + 328.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 189.0*Power(xii,18.0)*Power(xij,10.0)*
            (86069971575.0 + 2157712200.0*Power(rij,2.0)*Power(xij,2.0) - 
              158179560.0*Power(rij,4.0)*Power(xij,4.0) + 
              578816.0*Power(rij,6.0)*Power(xij,6.0) + 978.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 2.0*rij*Power(xii,29.0)*(2085060285.0 + 5450330025.0*Power(rij,2.0)*Power(xij,2.0) + 
              127424745.0*Power(rij,4.0)*Power(xij,4.0) - 
              1398276.0*Power(rij,6.0)*Power(xij,6.0) + 
              1159.0*Power(rij,8.0)*Power(xij,8.0)) - 
           378.0*Power(xii,22.0)*Power(xij,6.0)*
            (37244490525.0 - 2411839800.0*Power(rij,2.0)*Power(xij,2.0) + 
              92951775.0*Power(rij,4.0)*Power(xij,4.0) - 
              942172.0*Power(rij,6.0)*Power(xij,6.0) + 1292.0*Power(rij,8.0)*Power(xij,8.0)\
    ) - 27.0*Power(xii,16.0)*Power(xij,12.0)*
            (164245367475.0 + 26909517600.0*Power(rij,2.0)*Power(xij,2.0) + 
              62674920.0*Power(rij,4.0)*Power(xij,4.0) - 
              3885112.0*Power(rij,6.0)*Power(xij,6.0) + 
              2122.0*Power(rij,8.0)*Power(xij,8.0)) + 
           3.0*rij*Power(xii,27.0)*Power(xij,2.0)*
            (-63819198135.0 - 21841975890.0*Power(rij,2.0)*Power(xij,2.0) + 
              442430100.0*Power(rij,4.0)*Power(xij,4.0) - 
              2756664.0*Power(rij,6.0)*Power(xij,6.0) + 
              2296.0*Power(rij,8.0)*Power(xij,8.0)) + 
           rij*Power(xii,17.0)*Power(xij,12.0)*
            (4851990871875.0 + 21622847400.0*Power(rij,2.0)*Power(xij,2.0) - 
              2153738160.0*Power(rij,4.0)*Power(xij,4.0) + 
              3608388.0*Power(rij,6.0)*Power(xij,6.0) + 
              2318.0*Power(rij,8.0)*Power(xij,8.0)) + 
           18.0*Power(xii,14.0)*Power(xij,14.0)*
            (-23418646650.0 - 6922729800.0*Power(rij,2.0)*Power(xij,2.0) - 
              259958475.0*Power(rij,4.0)*Power(xij,4.0) - 
              697732.0*Power(rij,6.0)*Power(xij,6.0) + 3030.0*Power(rij,8.0)*Power(xij,8.0)\
    ) + 126.0*Power(xii,20.0)*Power(xij,8.0)*
            (-186637212225.0 + 13028280300.0*Power(rij,2.0)*Power(xij,2.0) - 
              116198775.0*Power(rij,4.0)*Power(xij,4.0) - 
              1266160.0*Power(rij,6.0)*Power(xij,6.0) + 
              4332.0*Power(rij,8.0)*Power(xij,8.0)) - 
           54.0*Power(xii,28.0)*(102965940.0 + 
              1089437580.0*Power(rij,2.0)*Power(xij,2.0) + 
              102508245.0*Power(rij,4.0)*Power(xij,4.0) - 
              1593144.0*Power(rij,6.0)*Power(xij,6.0) + 
              4538.0*Power(rij,8.0)*Power(xij,8.0)) + 
           63.0*Power(xii,26.0)*Power(xij,2.0)*
            (-4544129205.0 - 7396000920.0*Power(rij,2.0)*Power(xij,2.0) + 
              149614020.0*Power(rij,4.0)*Power(xij,4.0) - 
              1684112.0*Power(rij,6.0)*Power(xij,6.0) + 5922.0*Power(rij,8.0)*Power(xij,8.0)\
    )) + Power(E,2.0*rij*xii)*Power(xii,12.0)*
         (6.0*Power(xii,24.0)*Power(xij,6.0)*
            (1036901250.0 + 1900985625.0*rij*xij + 
              1728168750.0*Power(rij,2.0)*Power(xij,2.0) + 
              1036901250.0*Power(rij,3.0)*Power(xij,3.0) + 
              460845000.0*Power(rij,4.0)*Power(xij,4.0) + 
              161295750.0*Power(rij,5.0)*Power(xij,5.0) + 
              46084500.0*Power(rij,6.0)*Power(xij,6.0) + 
              9084900.0*Power(rij,7.0)*Power(xij,7.0) + 
              4082100.0*Power(rij,8.0)*Power(xij,8.0) + 
              121935.0*Power(rij,9.0)*Power(xij,9.0) - 
              21494.0*Power(rij,10.0)*Power(xij,10.0) - 
              766.0*Power(rij,11.0)*Power(xij,11.0)) + 
           5.0*Power(xii,28.0)*Power(xij,2.0)*
            (19646550.0 + 36018675.0*rij*xij + 32744250.0*Power(rij,2.0)*Power(xij,2.0) + 
              19646550.0*Power(rij,3.0)*Power(xij,3.0) + 
              8731800.0*Power(rij,4.0)*Power(xij,4.0) + 
              3056130.0*Power(rij,5.0)*Power(xij,5.0) + 
              873180.0*Power(rij,6.0)*Power(xij,6.0) + 
              207900.0*Power(rij,7.0)*Power(xij,7.0) + 
              41580.0*Power(rij,8.0)*Power(xij,8.0) + 
              6930.0*Power(rij,9.0)*Power(xij,9.0) + 
              924.0*Power(rij,10.0)*Power(xij,10.0) - 4.0*Power(rij,11.0)*Power(xij,11.0)) + 
           26334.0*Power(xii,16.0)*Power(xij,14.0)*
            (43880400.0 - 186686775.0*rij*xij + 
              576771750.0*Power(rij,2.0)*Power(xij,2.0) - 
              398603250.0*Power(rij,3.0)*Power(xij,3.0) + 
              72552600.0*Power(rij,4.0)*Power(xij,4.0) + 
              27903120.0*Power(rij,5.0)*Power(xij,5.0) - 
              342720.0*Power(rij,6.0)*Power(xij,6.0) - 
              574800.0*Power(rij,7.0)*Power(xij,7.0) - 
              50800.0*Power(rij,8.0)*Power(xij,8.0) - 945.0*Power(rij,9.0)*Power(xij,9.0) + 
              58.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) + 
           10.0*Power(xij,30.0)*(97302813300.0 + 89194245525.0*rij*xij + 
              38780106750.0*Power(rij,2.0)*Power(xij,2.0) + 
              10576392750.0*Power(rij,3.0)*Power(xij,3.0) + 
              2014551000.0*Power(rij,4.0)*Power(xij,4.0) + 
              282037140.0*Power(rij,5.0)*Power(xij,5.0) + 
              29688120.0*Power(rij,6.0)*Power(xij,6.0) + 
              2356200.0*Power(rij,7.0)*Power(xij,7.0) + 
              138600.0*Power(rij,8.0)*Power(xij,8.0) + 
              5775.0*Power(rij,9.0)*Power(xij,9.0) + 
              154.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) + 
           10.0*Power(xii,2.0)*Power(xij,28.0)*
            (4582499159700.0 + 3733416276975.0*rij*xij + 
              1428215931450.0*Power(rij,2.0)*Power(xij,2.0) + 
              338545295550.0*Power(rij,3.0)*Power(xij,3.0) + 
              55198697400.0*Power(rij,4.0)*Power(xij,4.0) + 
              6486854220.0*Power(rij,5.0)*Power(xij,5.0) + 
              558419400.0*Power(rij,6.0)*Power(xij,6.0) + 
              34939080.0*Power(rij,7.0)*Power(xij,7.0) + 
              1532520.0*Power(rij,8.0)*Power(xij,8.0) + 
              43285.0*Power(rij,9.0)*Power(xij,9.0) + 
              638.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) - 
           110.0*Power(xii,10.0)*Power(xij,20.0)*
            (-14063418170550.0 + 6795156458475.0*rij*xij + 
              2067471236250.0*Power(rij,2.0)*Power(xij,2.0) - 
              214664924250.0*Power(rij,3.0)*Power(xij,3.0) - 
              124416469800.0*Power(rij,4.0)*Power(xij,4.0) - 
              14935545450.0*Power(rij,5.0)*Power(xij,5.0) - 
              256688460.0*Power(rij,6.0)*Power(xij,6.0) + 
              105750900.0*Power(rij,7.0)*Power(xij,7.0) + 
              11502180.0*Power(rij,8.0)*Power(xij,8.0) + 
              518085.0*Power(rij,9.0)*Power(xij,9.0) + 
              9294.0*Power(rij,10.0)*Power(xij,10.0) + 2.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 55.0*Power(xii,20.0)*Power(xij,10.0)*
            (1730682450.0 + 3172917825.0*rij*xij + 
              2884470750.0*Power(rij,2.0)*Power(xij,2.0) + 
              1571960250.0*Power(rij,3.0)*Power(xij,3.0) + 
              1404081000.0*Power(rij,4.0)*Power(xij,4.0) - 
              426654270.0*Power(rij,5.0)*Power(xij,5.0) + 
              283536540.0*Power(rij,6.0)*Power(xij,6.0) + 
              39116700.0*Power(rij,7.0)*Power(xij,7.0) - 
              2659860.0*Power(rij,8.0)*Power(xij,8.0) - 
              528850.0*Power(rij,9.0)*Power(xij,9.0) - 
              18236.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 5.0*Power(xii,30.0)*(935550.0 + 1715175.0*rij*xij + 
              1559250.0*Power(rij,2.0)*Power(xij,2.0) + 
              935550.0*Power(rij,3.0)*Power(xij,3.0) + 
              415800.0*Power(rij,4.0)*Power(xij,4.0) + 
              145530.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              9900.0*Power(rij,7.0)*Power(xij,7.0) + 1980.0*Power(rij,8.0)*Power(xij,8.0) + 
              330.0*Power(rij,9.0)*Power(xij,9.0) + 44.0*Power(rij,10.0)*Power(xij,10.0) + 
              4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           13167.0*Power(xii,14.0)*Power(xij,16.0)*
            (-2319354450.0 + 8540029575.0*rij*xij - 
              7335672750.0*Power(rij,2.0)*Power(xij,2.0) + 
              1133154750.0*Power(rij,3.0)*Power(xij,3.0) + 
              575014200.0*Power(rij,4.0)*Power(xij,4.0) - 
              913710.0*Power(rij,5.0)*Power(xij,5.0) - 
              14863940.0*Power(rij,6.0)*Power(xij,6.0) - 
              1687300.0*Power(rij,7.0)*Power(xij,7.0) - 
              46900.0*Power(rij,8.0)*Power(xij,8.0) + 
              3210.0*Power(rij,9.0)*Power(xij,9.0) + 
              236.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           770.0*Power(xii,18.0)*Power(xij,12.0)*
            (329653800.0 + 654729075.0*rij*xij + 
              45785250.0*Power(rij,2.0)*Power(xij,2.0) + 
              1602483750.0*Power(rij,3.0)*Power(xij,3.0) - 
              915705000.0*Power(rij,4.0)*Power(xij,4.0) + 
              266036400.0*Power(rij,5.0)*Power(xij,5.0) + 
              63745920.0*Power(rij,6.0)*Power(xij,6.0) - 
              2304000.0*Power(rij,7.0)*Power(xij,7.0) - 
              1074240.0*Power(rij,8.0)*Power(xij,8.0) - 
              64635.0*Power(rij,9.0)*Power(xij,9.0) - 
              514.0*Power(rij,10.0)*Power(xij,10.0) + 34.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 385.0*Power(xii,12.0)*Power(xij,18.0)*
            (973565393850.0 - 1429122323475.0*rij*xij + 
              298281831750.0*Power(rij,2.0)*Power(xij,2.0) + 
              138841148250.0*Power(rij,3.0)*Power(xij,3.0) - 
              2454240600.0*Power(rij,4.0)*Power(xij,4.0) - 
              4925394810.0*Power(rij,5.0)*Power(xij,5.0) - 
              623832300.0*Power(rij,6.0)*Power(xij,6.0) - 
              19098540.0*Power(rij,7.0)*Power(xij,7.0) + 
              2083140.0*Power(rij,8.0)*Power(xij,8.0) + 
              212430.0*Power(rij,9.0)*Power(xij,9.0) + 
              7012.0*Power(rij,10.0)*Power(xij,10.0) + 68.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 14.0*Power(xii,26.0)*Power(xij,4.0)*
            (-70166250.0 - 128638125.0*rij*xij - 
              116943750.0*Power(rij,2.0)*Power(xij,2.0) - 
              70166250.0*Power(rij,3.0)*Power(xij,3.0) - 
              31185000.0*Power(rij,4.0)*Power(xij,4.0) - 
              10914750.0*Power(rij,5.0)*Power(xij,5.0) - 
              3118500.0*Power(rij,6.0)*Power(xij,6.0) - 
              742500.0*Power(rij,7.0)*Power(xij,7.0) - 
              148500.0*Power(rij,8.0)*Power(xij,8.0) - 
              32615.0*Power(rij,9.0)*Power(xij,9.0) - 
              154.0*Power(rij,10.0)*Power(xij,10.0) + 74.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 7.0*Power(xii,4.0)*Power(xij,26.0)*
            (-69822945249750.0 - 46669577290875.0*rij*xij - 
              14025037430250.0*Power(rij,2.0)*Power(xij,2.0) - 
              2430881664750.0*Power(rij,3.0)*Power(xij,3.0) - 
              251629270200.0*Power(rij,4.0)*Power(xij,4.0) - 
              12434519790.0*Power(rij,5.0)*Power(xij,5.0) + 
              452930940.0*Power(rij,6.0)*Power(xij,6.0) + 
              131125500.0*Power(rij,7.0)*Power(xij,7.0) + 
              11018700.0*Power(rij,8.0)*Power(xij,8.0) + 
              514470.0*Power(rij,9.0)*Power(xij,9.0) + 
              13332.0*Power(rij,10.0)*Power(xij,10.0) + 
              148.0*Power(rij,11.0)*Power(xij,11.0)) - 
           50.0*Power(xii,8.0)*Power(xij,22.0)*
            (-51768833574150.0 - 5003280391725.0*rij*xij + 
              4493439477450.0*Power(rij,2.0)*Power(xij,2.0) + 
              1286866176750.0*Power(rij,3.0)*Power(xij,3.0) + 
              111437476920.0*Power(rij,4.0)*Power(xij,4.0) - 
              6620313546.0*Power(rij,5.0)*Power(xij,5.0) - 
              2406603276.0*Power(rij,6.0)*Power(xij,6.0) - 
              242686620.0*Power(rij,7.0)*Power(xij,7.0) - 
              12228876.0*Power(rij,8.0)*Power(xij,8.0) - 
              256223.0*Power(rij,9.0)*Power(xij,9.0) + 
              2486.0*Power(rij,10.0)*Power(xij,10.0) + 158.0*Power(rij,11.0)*Power(xij,11.0)\
    ) + 25.0*Power(xii,22.0)*Power(xij,8.0)*
            (-1119853350.0 - 2053064475.0*rij*xij - 
              1866422250.0*Power(rij,2.0)*Power(xij,2.0) - 
              1119853350.0*Power(rij,3.0)*Power(xij,3.0) - 
              497712600.0*Power(rij,4.0)*Power(xij,4.0) - 
              194415606.0*Power(rij,5.0)*Power(xij,5.0) - 
              9338868.0*Power(rij,6.0)*Power(xij,6.0) - 
              31217076.0*Power(rij,7.0)*Power(xij,7.0) - 
              2256804.0*Power(rij,8.0)*Power(xij,8.0) + 
              246774.0*Power(rij,9.0)*Power(xij,9.0) + 
              22836.0*Power(rij,10.0)*Power(xij,10.0) + 
              316.0*Power(rij,11.0)*Power(xij,11.0)) + 
           3.0*Power(xii,6.0)*Power(xij,24.0)*
            (596006592662250.0 + 266778699697125.0*rij*xij + 
              37515651153750.0*Power(rij,2.0)*Power(xij,2.0) - 
              2214626163750.0*Power(rij,3.0)*Power(xij,3.0) - 
              1538075107800.0*Power(rij,4.0)*Power(xij,4.0) - 
              248955308910.0*Power(rij,5.0)*Power(xij,5.0) - 
              21434337540.0*Power(rij,6.0)*Power(xij,6.0) - 
              957980100.0*Power(rij,7.0)*Power(xij,7.0) - 
              4874100.0*Power(rij,8.0)*Power(xij,8.0) + 
              1831830.0*Power(rij,9.0)*Power(xij,9.0) + 
              91828.0*Power(rij,10.0)*Power(xij,10.0) + 
              1532.0*Power(rij,11.0)*Power(xij,11.0))))/
      (4.67775e6*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),21.0))
    ;
  }
  return S;
}

double Slater_5S_1S(double rij,double xii,double xij)
{
  return Slater_1S_5S(rij,xij,xii);
}

double Slater_5S_2S(double rij,double xii,double xij)
{
  return Slater_2S_5S(rij,xij,xii);
}

double Slater_5S_3S(double rij,double xii,double xij)
{
  return Slater_3S_5S(rij,xij,xii);
}

double Slater_5S_4S(double rij,double xii,double xij)
{
  return Slater_4S_5S(rij,xij,xii);
}

static double Slater_6S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     -(3722690410399436636160000.0 - 3722690410399436636160000.0*Power(E,2.0*rij*xii) + 
         6924936452406883646360625.0*rij*xii + 
         6404492084014894020401250.0*Power(rij,2.0)*Power(xii,2.0) + 
         3925597144715015967697500.0*Power(rij,3.0)*Power(xii,3.0) + 
         1793665117676464332300000.0*Power(rij,4.0)*Power(xii,4.0) + 
         651524259419605812240000.0*Power(rij,5.0)*Power(xii,5.0) + 
         195930326813816174580000.0*Power(rij,6.0)*Power(xii,6.0) + 
         50160444229615663944000.0*Power(rij,7.0)*Power(xii,7.0) + 
         11155494661051156416000.0*Power(rij,8.0)*Power(xii,8.0) + 
         2188143143401479264000.0*Power(rij,9.0)*Power(xii,9.0) + 
         382976811299821939200.0*Power(rij,10.0)*Power(xii,10.0) + 
         60350063176103500800.0*Power(rij,11.0)*Power(xii,11.0) + 
         8621483857123737600.0*Power(rij,12.0)*Power(xii,12.0) + 
         1122323342347468800.0*Power(rij,13.0)*Power(xii,13.0) + 
         133609921708032000.0*Power(rij,14.0)*Power(xii,14.0) + 
         14575627822694400.0*Power(rij,15.0)*Power(xii,15.0) + 
         1457562782269440.0*Power(rij,16.0)*Power(xii,16.0) + 
         133371757854720.0*Power(rij,17.0)*Power(xii,17.0) + 
         11114313154560.0*Power(rij,18.0)*Power(xii,18.0) + 
         835662643200.0*Power(rij,19.0)*Power(xii,19.0) + 
         55710842880.0*Power(rij,20.0)*Power(xii,20.0) + 
         3183476736.0*Power(rij,21.0)*Power(xii,21.0) + 
         144703488.0*Power(rij,22.0)*Power(xii,22.0) + 
         4194304.0*Power(rij,23.0)*Power(xii,23.0))/
      (3.7226904103994365e24*Power(E,2.0*rij*xii)*rij)
    ;
  }
  else {
    S =     (2806650.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),23.0) + 
        Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-1056.0*Power(rij,10.0)*Power(xii,42.0) - 12.0*Power(rij,11.0)*Power(xii,43.0) + 
           2806650.0*Power(xij,32.0) + 5145525.0*rij*xii*Power(xij,32.0) - 
           88.0*Power(rij,9.0)*Power(xii,41.0)*(510.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           935550.0*Power(xii,2.0)*Power(xij,30.0)*
            (-69.0 + 5.0*Power(rij,2.0)*Power(xij,2.0)) + 
           467775.0*rij*Power(xii,3.0)*Power(xij,30.0)*
            (-253.0 + 6.0*Power(rij,2.0)*Power(xij,2.0)) - 
           132.0*Power(rij,8.0)*Power(xii,40.0)*
            (9180.0 + 89.0*Power(rij,2.0)*Power(xij,2.0)) + 
           311850.0*Power(xii,4.0)*Power(xij,28.0)*
            (2277.0 - 345.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,4.0)*Power(xij,4.0)) + 
           31185.0*rij*Power(xii,5.0)*Power(xij,28.0)*
            (41745.0 - 2070.0*Power(rij,2.0)*Power(xij,2.0) + 
              14.0*Power(rij,4.0)*Power(xij,4.0)) + 
           1980.0*Power(rij,6.0)*Power(xii,38.0)*
            (-162792.0 - 11859.0*Power(rij,2.0)*Power(xij,2.0) + 
              41.0*Power(rij,4.0)*Power(xij,4.0)) + 
           22.0*Power(rij,7.0)*Power(xii,39.0)*
            (-1046520.0 - 30885.0*Power(rij,2.0)*Power(xij,2.0) + 
              44.0*Power(rij,4.0)*Power(xij,4.0)) + 
           62370.0*Power(xii,6.0)*Power(xij,26.0)*
            (-79695.0 + 18975.0*Power(rij,2.0)*Power(xij,2.0) - 
              460.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) - 
           110.0*Power(rij,5.0)*Power(xii,37.0)*
            (30767688.0 + 4989438.0*Power(rij,2.0)*Power(xij,2.0) - 
              25359.0*Power(rij,4.0)*Power(xij,4.0) + 20.0*Power(rij,6.0)*Power(xij,6.0)) + 
           1485.0*rij*Power(xii,7.0)*Power(xij,26.0)*
            (-6136515.0 + 478170.0*Power(rij,2.0)*Power(xij,2.0) - 
              6762.0*Power(rij,4.0)*Power(xij,4.0) + 20.0*Power(rij,6.0)*Power(xij,6.0)) - 
           132.0*Power(rij,4.0)*Power(xii,36.0)*
            (201455100.0 + 69647445.0*Power(rij,2.0)*Power(xij,2.0) - 
              318735.0*Power(rij,4.0)*Power(xij,4.0) + 353.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 495.0*rij*Power(xii,9.0)*Power(xij,24.0)*
            (92047725.0 - 10041570.0*Power(rij,2.0)*Power(xij,2.0) + 
              223146.0*Power(rij,4.0)*Power(xij,4.0) - 
              1380.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           2970.0*Power(xii,8.0)*Power(xij,24.0)*
            (8367975.0 - 2789325.0*Power(rij,2.0)*Power(xij,2.0) + 
              106260.0*Power(rij,4.0)*Power(xij,4.0) - 
              966.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           22.0*Power(rij,3.0)*Power(xii,35.0)*
            (6950200950.0 + 5142653145.0*Power(rij,2.0)*Power(xij,2.0) + 
              7644510.0*Power(rij,4.0)*Power(xij,4.0) - 
              235635.0*Power(rij,6.0)*Power(xij,6.0) + 124.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 132.0*Power(rij,2.0)*Power(xii,34.0)*
            (4633467300.0 + 7767871650.0*Power(rij,2.0)*Power(xij,2.0) + 
              160904205.0*Power(rij,4.0)*Power(xij,4.0) - 
              2493315.0*Power(rij,6.0)*Power(xij,6.0) + 
              5281.0*Power(rij,8.0)*Power(xij,8.0)) - 
           495.0*rij*Power(xii,27.0)*Power(xij,6.0)*
            (8395934795325.0 - 439434024750.0*Power(rij,2.0)*Power(xij,2.0) + 
              11948496210.0*Power(rij,4.0)*Power(xij,4.0) - 
              118623972.0*Power(rij,6.0)*Power(xij,6.0) + 
              248906.0*Power(rij,8.0)*Power(xij,8.0)) + 
           11.0*rij*Power(xii,15.0)*Power(xij,18.0)*
            (1488922594425.0 + 252796524750.0*Power(rij,2.0)*Power(xij,2.0) + 
              6172031250.0*Power(rij,4.0)*Power(xij,4.0) + 
              104343660.0*Power(rij,6.0)*Power(xij,6.0) + 
              66810.0*Power(rij,8.0)*Power(xij,8.0) - 88.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 66.0*Power(xii,10.0)*Power(xij,22.0)*
            (-1430923725.0 + 627598125.0*Power(rij,2.0)*Power(xij,2.0) - 
              33471900.0*Power(rij,4.0)*Power(xij,4.0) + 
              478170.0*Power(rij,6.0)*Power(xij,6.0) - 
              2070.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,10.0)*Power(xij,10.0)) - 
           1518.0*Power(xii,12.0)*Power(xij,20.0)*
            (-186642225.0 + 103690125.0*Power(rij,2.0)*Power(xij,2.0) - 
              7276500.0*Power(rij,4.0)*Power(xij,4.0) + 
              145530.0*Power(rij,6.0)*Power(xij,6.0) - 
              990.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,10.0)*Power(xij,10.0)) + 
           3.0*rij*Power(xii,11.0)*Power(xij,22.0)*
            (-57713923575.0 + 8284295250.0*Power(rij,2.0)*Power(xij,2.0) - 
              257733630.0*Power(rij,4.0)*Power(xij,4.0) + 
              2504700.0*Power(rij,6.0)*Power(xij,6.0) - 
              7590.0*Power(rij,8.0)*Power(xij,8.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) + 
           11.0*rij*Power(xii,13.0)*Power(xij,20.0)*
            (56066193225.0 - 6918959250.0*Power(rij,2.0)*Power(xij,2.0) + 
              430816050.0*Power(rij,4.0)*Power(xij,4.0) - 
              3349620.0*Power(rij,6.0)*Power(xij,6.0) + 
              33690.0*Power(rij,8.0)*Power(xij,8.0) + 8.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 55.0*rij*Power(xii,17.0)*Power(xij,16.0)*
            (7416068831325.0 + 658162968750.0*Power(rij,2.0)*Power(xij,2.0) + 
              11421785970.0*Power(rij,4.0)*Power(xij,4.0) - 
              22800852.0*Power(rij,6.0)*Power(xij,6.0) - 
              224214.0*Power(rij,8.0)*Power(xij,8.0) + 40.0*Power(rij,10.0)*Power(xij,10.0)\
    ) - 198.0*Power(xii,14.0)*Power(xij,18.0)*
            (12601626975.0 + 2529410625.0*Power(rij,2.0)*Power(xij,2.0) + 
              582340500.0*Power(rij,4.0)*Power(xij,4.0) + 
              3239250.0*Power(rij,6.0)*Power(xij,6.0) + 
              132690.0*Power(rij,8.0)*Power(xij,8.0) + 74.0*Power(rij,10.0)*Power(xij,10.0)\
    ) - 231.0*rij*Power(xii,25.0)*Power(xij,8.0)*
            (21444497452125.0 - 909858116250.0*Power(rij,2.0)*Power(xij,2.0) + 
              1447333650.0*Power(rij,4.0)*Power(xij,4.0) + 
              178686540.0*Power(rij,6.0)*Power(xij,6.0) - 
              747270.0*Power(rij,8.0)*Power(xij,8.0) + 
              184.0*Power(rij,10.0)*Power(xij,10.0)) - 
           198.0*Power(xii,20.0)*Power(xij,12.0)*
            (42449899182075.0 + 4344172457625.0*Power(rij,2.0)*Power(xij,2.0) - 
              85249741500.0*Power(rij,4.0)*Power(xij,4.0) - 
              1059301110.0*Power(rij,6.0)*Power(xij,6.0) + 
              6582370.0*Power(rij,8.0)*Power(xij,8.0) + 
              194.0*Power(rij,10.0)*Power(xij,10.0)) + 
           11.0*rij*Power(xii,19.0)*Power(xij,14.0)*
            (239338679943825.0 + 8851966719750.0*Power(rij,2.0)*Power(xij,2.0) - 
              112537092150.0*Power(rij,4.0)*Power(xij,4.0) - 
              1100275380.0*Power(rij,6.0)*Power(xij,6.0) + 
              2919090.0*Power(rij,8.0)*Power(xij,8.0) + 
              248.0*Power(rij,10.0)*Power(xij,10.0)) - 
           330.0*Power(xii,28.0)*Power(xij,4.0)*
            (4860066085875.0 + 2524912849305.0*Power(rij,2.0)*Power(xij,2.0) - 
              109538431380.0*Power(rij,4.0)*Power(xij,4.0) + 
              1633704282.0*Power(rij,6.0)*Power(xij,6.0) - 
              6421278.0*Power(rij,8.0)*Power(xij,8.0) + 
              322.0*Power(rij,10.0)*Power(xij,10.0)) + 
           33.0*rij*Power(xii,29.0)*Power(xij,4.0)*
            (-31641507079875.0 - 2157639318450.0*Power(rij,2.0)*Power(xij,2.0) + 
              74910015810.0*Power(rij,4.0)*Power(xij,4.0) - 
              522003060.0*Power(rij,6.0)*Power(xij,6.0) - 
              250470.0*Power(rij,8.0)*Power(xij,8.0) + 
              1288.0*Power(rij,10.0)*Power(xij,10.0)) - 
           330.0*Power(xii,18.0)*Power(xij,14.0)*
            (4867016286825.0 + 1199363925375.0*Power(rij,2.0)*Power(xij,2.0) + 
              26817947100.0*Power(rij,4.0)*Power(xij,4.0) - 
              167333418.0*Power(rij,6.0)*Power(xij,6.0) - 
              1476138.0*Power(rij,8.0)*Power(xij,8.0) + 
              1294.0*Power(rij,10.0)*Power(xij,10.0)) + 
           66.0*Power(xii,16.0)*Power(xij,16.0)*
            (-1657759205025.0 - 682207855875.0*Power(rij,2.0)*Power(xij,2.0) - 
              31509229500.0*Power(rij,4.0)*Power(xij,4.0) - 
              492146550.0*Power(rij,6.0)*Power(xij,6.0) - 
              11910.0*Power(rij,8.0)*Power(xij,8.0) + 
              2594.0*Power(rij,10.0)*Power(xij,10.0)) + 
           1386.0*Power(xii,26.0)*Power(xij,6.0)*
            (-6066588045375.0 + 98854491375.0*Power(rij,2.0)*Power(xij,2.0) - 
              12496954500.0*Power(rij,4.0)*Power(xij,4.0) + 
              420813750.0*Power(rij,6.0)*Power(xij,6.0) - 
              2881210.0*Power(rij,8.0)*Power(xij,8.0) + 
              2622.0*Power(rij,10.0)*Power(xij,10.0)) + 
           11.0*rij*Power(xii,23.0)*Power(xij,10.0)*
            (149900659402725.0 - 26541339882750.0*Power(rij,2.0)*Power(xij,2.0) + 
              594745455150.0*Power(rij,4.0)*Power(xij,4.0) - 
              1399125420.0*Power(rij,6.0)*Power(xij,6.0) - 
              7887390.0*Power(rij,8.0)*Power(xij,8.0) + 
              4232.0*Power(rij,10.0)*Power(xij,10.0)) - 
           11.0*rij*Power(xii,31.0)*Power(xij,2.0)*
            (7685082491625.0 + 5034333946950.0*Power(rij,2.0)*Power(xij,2.0) - 
              108088893990.0*Power(rij,4.0)*Power(xij,4.0) + 
              1254174300.0*Power(rij,6.0)*Power(xij,6.0) - 
              6355950.0*Power(rij,8.0)*Power(xij,8.0) + 
              4232.0*Power(rij,10.0)*Power(xij,10.0)) - 
           462.0*Power(xii,24.0)*Power(xij,8.0)*
            (40495013164125.0 - 3973079865375.0*Power(rij,2.0)*Power(xij,2.0) + 
              110288047500.0*Power(rij,4.0)*Power(xij,4.0) - 
              381623130.0*Power(rij,6.0)*Power(xij,6.0) - 
              4811370.0*Power(rij,8.0)*Power(xij,8.0) + 
              9338.0*Power(rij,10.0)*Power(xij,10.0)) + 
           198.0*Power(xii,32.0)*(-9126526500.0 - 
              152866565775.0*Power(rij,2.0)*Power(xij,2.0) - 
              32383266300.0*Power(rij,4.0)*Power(xij,4.0) + 
              709444890.0*Power(rij,6.0)*Power(xij,6.0) - 
              5562070.0*Power(rij,8.0)*Power(xij,8.0) + 
              11042.0*Power(rij,10.0)*Power(xij,10.0)) + 
           2.0*rij*Power(xii,33.0)*(-764522104500.0 - 
              3357151476525.0*Power(rij,2.0)*Power(xij,2.0) - 
              242177564475.0*Power(rij,4.0)*Power(xij,4.0) + 
              4513719870.0*Power(rij,6.0)*Power(xij,6.0) - 
              20531775.0*Power(rij,8.0)*Power(xij,8.0) + 
              11236.0*Power(rij,10.0)*Power(xij,10.0)) - 
           rij*Power(xii,21.0)*Power(xij,12.0)*
            (-5533525427435775.0 + 138591131159250.0*Power(rij,2.0)*Power(xij,2.0) + 
              2815739907750.0*Power(rij,4.0)*Power(xij,4.0) - 
              32922004500.0*Power(rij,6.0)*Power(xij,6.0) + 
              11347050.0*Power(rij,8.0)*Power(xij,8.0) + 
              22472.0*Power(rij,10.0)*Power(xij,10.0)) + 
           66.0*Power(xii,22.0)*Power(xij,10.0)*
            (-283522589265825.0 + 7639225988625.0*Power(rij,2.0)*Power(xij,2.0) + 
              480728209500.0*Power(rij,4.0)*Power(xij,4.0) - 
              8458349130.0*Power(rij,6.0)*Power(xij,6.0) + 
              9771090.0*Power(rij,8.0)*Power(xij,8.0) + 
              31786.0*Power(rij,10.0)*Power(xij,10.0)) - 
           66.0*Power(xii,30.0)*Power(xij,2.0)*
            (1678609807875.0 + 4713298976925.0*Power(rij,2.0)*Power(xij,2.0) - 
              30578971500.0*Power(rij,4.0)*Power(xij,4.0) + 
              53723250.0*Power(rij,6.0)*Power(xij,6.0) - 
              9140190.0*Power(rij,8.0)*Power(xij,8.0) + 
              38042.0*Power(rij,10.0)*Power(xij,10.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,14.0)*
         (-302841.0*Power(xii,16.0)*Power(xij,16.0)*
            (-361285650.0 + 1346857875.0*rij*xij - 
              1306923750.0*Power(rij,2.0)*Power(xij,2.0) + 
              321527250.0*Power(rij,3.0)*Power(xij,3.0) + 
              55737000.0*Power(rij,4.0)*Power(xij,4.0) - 
              9297750.0*Power(rij,5.0)*Power(xij,5.0) - 
              1843380.0*Power(rij,6.0)*Power(xij,6.0) - 
              50820.0*Power(rij,7.0)*Power(xij,7.0) + 
              7340.0*Power(rij,8.0)*Power(xij,8.0) + 570.0*Power(rij,9.0)*Power(xij,9.0) + 
              12.0*Power(rij,10.0)*Power(xij,10.0)) + 
           12.0*Power(xij,32.0)*(150587687250.0 + 127420350750.0*rij*xij + 
              50968140300.0*Power(rij,2.0)*Power(xij,2.0) + 
              12742035075.0*Power(rij,3.0)*Power(xij,3.0) + 
              2216006100.0*Power(rij,4.0)*Power(xij,4.0) + 
              282037140.0*Power(rij,5.0)*Power(xij,5.0) + 
              26860680.0*Power(rij,6.0)*Power(xij,6.0) + 
              1918620.0*Power(rij,7.0)*Power(xij,7.0) + 
              100980.0*Power(rij,8.0)*Power(xij,8.0) + 
              3740.0*Power(rij,9.0)*Power(xij,9.0) + 88.0*Power(rij,10.0)*Power(xij,10.0) + 
              Power(rij,11.0)*Power(xij,11.0)) - 
           3.0*Power(xii,32.0)*(935550.0 + 1715175.0*rij*xij + 
              1559250.0*Power(rij,2.0)*Power(xij,2.0) + 
              935550.0*Power(rij,3.0)*Power(xij,3.0) + 
              415800.0*Power(rij,4.0)*Power(xij,4.0) + 
              145530.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              9900.0*Power(rij,7.0)*Power(xij,7.0) + 1980.0*Power(rij,8.0)*Power(xij,8.0) + 
              330.0*Power(rij,9.0)*Power(xij,9.0) + 44.0*Power(rij,10.0)*Power(xij,10.0) + 
              4.0*Power(rij,11.0)*Power(xij,11.0)) - 
           11.0*Power(xii,30.0)*Power(xij,2.0)*
            (-5868450.0 - 10758825.0*rij*xij - 9780750.0*Power(rij,2.0)*Power(xij,2.0) - 
              5868450.0*Power(rij,3.0)*Power(xij,3.0) - 
              2608200.0*Power(rij,4.0)*Power(xij,4.0) - 
              912870.0*Power(rij,5.0)*Power(xij,5.0) - 
              260820.0*Power(rij,6.0)*Power(xij,6.0) - 
              62100.0*Power(rij,7.0)*Power(xij,7.0) - 
              12420.0*Power(rij,8.0)*Power(xij,8.0) - 
              2070.0*Power(rij,9.0)*Power(xij,9.0) - 
              276.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) - 
           5313.0*Power(xii,14.0)*Power(xij,18.0)*
            (-302299148250.0 + 495525217275.0*rij*xij - 
              161894625750.0*Power(rij,2.0)*Power(xij,2.0) - 
              26085287250.0*Power(rij,3.0)*Power(xij,3.0) + 
              5971779000.0*Power(rij,4.0)*Power(xij,4.0) + 
              1231357050.0*Power(rij,5.0)*Power(xij,5.0) + 
              33184620.0*Power(rij,6.0)*Power(xij,6.0) - 
              7768980.0*Power(rij,7.0)*Power(xij,7.0) - 
              751620.0*Power(rij,8.0)*Power(xij,8.0) - 
              23190.0*Power(rij,9.0)*Power(xij,9.0) - 
              20.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           5313.0*Power(xii,18.0)*Power(xij,14.0)*
            (469625850.0 - 3082655475.0*rij*xij + 
              8474631750.0*Power(rij,2.0)*Power(xij,2.0) - 
              6813281250.0*Power(rij,3.0)*Power(xij,3.0) + 
              1665711000.0*Power(rij,4.0)*Power(xij,4.0) + 
              232996050.0*Power(rij,5.0)*Power(xij,5.0) - 
              39477060.0*Power(rij,6.0)*Power(xij,6.0) - 
              6196500.0*Power(rij,7.0)*Power(xij,7.0) - 
              121380.0*Power(rij,8.0)*Power(xij,8.0) + 
              16330.0*Power(rij,9.0)*Power(xij,9.0) + 
              812.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           11.0*Power(xii,2.0)*Power(xij,30.0)*
            (10071658847250.0 + 7685082491625.0*rij*xij + 
              2751598183950.0*Power(rij,2.0)*Power(xij,2.0) + 
              610391177550.0*Power(rij,3.0)*Power(xij,3.0) + 
              93214459800.0*Power(rij,4.0)*Power(xij,4.0) + 
              10285306290.0*Power(rij,5.0)*Power(xij,5.0) + 
              835769340.0*Power(rij,6.0)*Power(xij,6.0) + 
              49894380.0*Power(rij,7.0)*Power(xij,7.0) + 
              2134620.0*Power(rij,8.0)*Power(xij,8.0) + 
              61770.0*Power(rij,9.0)*Power(xij,9.0) + 
              1068.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) \
    + 11.0*Power(xii,28.0)*Power(xij,4.0)*
            (-64552950.0 - 118347075.0*rij*xij - 
              107588250.0*Power(rij,2.0)*Power(xij,2.0) - 
              64552950.0*Power(rij,3.0)*Power(xij,3.0) - 
              28690200.0*Power(rij,4.0)*Power(xij,4.0) - 
              10041570.0*Power(rij,5.0)*Power(xij,5.0) - 
              2869020.0*Power(rij,6.0)*Power(xij,6.0) - 
              683100.0*Power(rij,7.0)*Power(xij,7.0) - 
              136620.0*Power(rij,8.0)*Power(xij,8.0) - 
              33690.0*Power(rij,9.0)*Power(xij,9.0) + 
              1332.0*Power(rij,10.0)*Power(xij,10.0) + 88.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 11.0*Power(xii,4.0)*Power(xij,28.0)*
            (-145801982576250.0 - 94924521239625.0*rij*xij - 
              28279793861550.0*Power(rij,2.0)*Power(xij,2.0) - 
              5034333946950.0*Power(rij,3.0)*Power(xij,3.0) - 
              582898793400.0*Power(rij,4.0)*Power(xij,4.0) - 
              44032284450.0*Power(rij,5.0)*Power(xij,5.0) - 
              1930850460.0*Power(rij,6.0)*Power(xij,6.0) - 
              15289020.0*Power(rij,7.0)*Power(xij,7.0) + 
              3824820.0*Power(rij,8.0)*Power(xij,8.0) + 
              253590.0*Power(rij,9.0)*Power(xij,9.0) + 
              7380.0*Power(rij,10.0)*Power(xij,10.0) + 88.0*Power(rij,11.0)*Power(xij,11.0)) \
    - 253.0*Power(xii,20.0)*Power(xij,12.0)*
            (1119853350.0 + 2437660575.0*rij*xij - 
              1979538750.0*Power(rij,2.0)*Power(xij,2.0) + 
              10991153250.0*Power(rij,3.0)*Power(xij,3.0) - 
              8219799000.0*Power(rij,4.0)*Power(xij,4.0) + 
              2482996950.0*Power(rij,5.0)*Power(xij,5.0) + 
              218260980.0*Power(rij,6.0)*Power(xij,6.0) - 
              47838060.0*Power(rij,7.0)*Power(xij,7.0) - 
              5151420.0*Power(rij,8.0)*Power(xij,8.0) - 
              44850.0*Power(rij,9.0)*Power(xij,9.0) + 
              8292.0*Power(rij,10.0)*Power(xij,10.0) + 184.0*Power(rij,11.0)*Power(xij,11.0)\
    ) + 253.0*Power(xii,12.0)*Power(xij,20.0)*
            (33221660229450.0 - 21871642005675.0*rij*xij - 
              1992841562250.0*Power(rij,2.0)*Power(xij,2.0) + 
              1153971299250.0*Power(rij,3.0)*Power(xij,3.0) + 
              201395565000.0*Power(rij,4.0)*Power(xij,4.0) + 
              1321478550.0*Power(rij,5.0)*Power(xij,5.0) - 
              2305327500.0*Power(rij,6.0)*Power(xij,6.0) - 
              232090380.0*Power(rij,7.0)*Power(xij,7.0) - 
              8375580.0*Power(rij,8.0)*Power(xij,8.0) + 
              32670.0*Power(rij,9.0)*Power(xij,9.0) + 
              9924.0*Power(rij,10.0)*Power(xij,10.0) + 184.0*Power(rij,11.0)*Power(xij,11.0)\
    ) + 11.0*Power(xii,6.0)*Power(xij,26.0)*
            (764390093717250.0 + 377817065789625.0*rij*xij + 
              75747385479150.0*Power(rij,2.0)*Power(xij,2.0) + 
              6472917955350.0*Power(rij,3.0)*Power(xij,3.0) - 
              183473829000.0*Power(rij,4.0)*Power(xij,4.0) - 
              108088893990.0*Power(rij,5.0)*Power(xij,5.0) - 
              12770008020.0*Power(rij,6.0)*Power(xij,6.0) - 
              820676340.0*Power(rij,7.0)*Power(xij,7.0) - 
              29919780.0*Power(rij,8.0)*Power(xij,8.0) - 
              471270.0*Power(rij,9.0)*Power(xij,9.0) + 
              4236.0*Power(rij,10.0)*Power(xij,10.0) + 200.0*Power(rij,11.0)*Power(xij,11.0)\
    ) - 11.0*Power(xii,26.0)*Power(xij,6.0)*
            (-451870650.0 - 828429525.0*rij*xij - 
              753117750.0*Power(rij,2.0)*Power(xij,2.0) - 
              451870650.0*Power(rij,3.0)*Power(xij,3.0) - 
              200831400.0*Power(rij,4.0)*Power(xij,4.0) - 
              70290990.0*Power(rij,5.0)*Power(xij,5.0) - 
              20083140.0*Power(rij,6.0)*Power(xij,6.0) - 
              3349620.0*Power(rij,7.0)*Power(xij,7.0) - 
              2388420.0*Power(rij,8.0)*Power(xij,8.0) + 
              66810.0*Power(rij,9.0)*Power(xij,9.0) + 
              15564.0*Power(rij,10.0)*Power(xij,10.0) + 
              200.0*Power(rij,11.0)*Power(xij,11.0)) - 
           11.0*Power(xii,24.0)*Power(xij,8.0)*
            (2259353250.0 + 4142147625.0*rij*xij + 
              3765588750.0*Power(rij,2.0)*Power(xij,2.0) + 
              2259353250.0*Power(rij,3.0)*Power(xij,3.0) + 
              1004157000.0*Power(rij,4.0)*Power(xij,4.0) + 
              430816050.0*Power(rij,5.0)*Power(xij,5.0) - 
              58306500.0*Power(rij,6.0)*Power(xij,6.0) + 
              104343660.0*Power(rij,7.0)*Power(xij,7.0) - 
              71460.0*Power(rij,8.0)*Power(xij,8.0) - 
              1121070.0*Power(rij,9.0)*Power(xij,9.0) - 
              38820.0*Power(rij,10.0)*Power(xij,10.0) + 
              248.0*Power(rij,11.0)*Power(xij,11.0)) + 
           11.0*Power(xii,8.0)*Power(xij,24.0)*
            (1700790552893250.0 + 450334446494625.0*rij*xij - 
              12455665913250.0*Power(rij,2.0)*Power(xij,2.0) - 
              19774531113750.0*Power(rij,3.0)*Power(xij,3.0) - 
              3286152941400.0*Power(rij,4.0)*Power(xij,4.0) - 
              224730047430.0*Power(rij,5.0)*Power(xij,5.0) + 
              322339500.0*Power(rij,6.0)*Power(xij,6.0) + 
              1254174300.0*Power(rij,7.0)*Power(xij,7.0) + 
              100117260.0*Power(rij,8.0)*Power(xij,8.0) + 
              3733050.0*Power(rij,9.0)*Power(xij,9.0) + 
              63372.0*Power(rij,10.0)*Power(xij,10.0) + 
              248.0*Power(rij,11.0)*Power(xij,11.0)) + 
           Power(xii,22.0)*Power(xij,10.0)*
            (94440965850.0 + 173141770725.0*rij*xij + 
              157401609750.0*Power(rij,2.0)*Power(xij,2.0) + 
              76108551750.0*Power(rij,3.0)*Power(xij,3.0) + 
              115303419000.0*Power(rij,4.0)*Power(xij,4.0) - 
              67892343750.0*Power(rij,5.0)*Power(xij,5.0) + 
              32481672300.0*Power(rij,6.0)*Power(xij,6.0) + 
              1254046860.0*Power(rij,7.0)*Power(xij,7.0) - 
              487125540.0*Power(rij,8.0)*Power(xij,8.0) - 
              32109990.0*Power(rij,9.0)*Power(xij,9.0) + 
              38412.0*Power(rij,10.0)*Power(xij,10.0) + 
              22472.0*Power(rij,11.0)*Power(xij,11.0)) - 
           Power(xii,10.0)*Power(xij,22.0)*
            (-18712490891544450.0 + 1648907253429975.0*rij*xij + 
              1835562897803250.0*Power(rij,2.0)*Power(xij,2.0) + 
              210177224853750.0*Power(rij,3.0)*Power(xij,3.0) - 
              17320778937000.0*Power(rij,4.0)*Power(xij,4.0) - 
              5914505623950.0*Power(rij,5.0)*Power(xij,5.0) - 
              539122413060.0*Power(rij,6.0)*Power(xij,6.0) - 
              17226100980.0*Power(rij,7.0)*Power(xij,7.0) + 
              603252540.0*Power(rij,8.0)*Power(xij,8.0) + 
              69915450.0*Power(rij,9.0)*Power(xij,9.0) + 
              2186316.0*Power(rij,10.0)*Power(xij,10.0) + 
              22472.0*Power(rij,11.0)*Power(xij,11.0))))/
      (2.80665e6*Power(E,2.0*rij*(xii + xij))*rij*
        Power(Power(xii,2.0) - Power(xij,2.0),23.0))
    ;
  }
  return S;
}

double Slater_6S_1S(double rij,double xii,double xij)
{
  return Slater_1S_6S(rij,xij,xii);
}

double Slater_6S_2S(double rij,double xii,double xij)
{
  return Slater_2S_6S(rij,xij,xii);
}

double Slater_6S_3S(double rij,double xii,double xij)
{
  return Slater_3S_6S(rij,xij,xii);
}

double Slater_6S_4S(double rij,double xii,double xij)
{
  return Slater_4S_6S(rij,xij,xii);
}

double Slater_6S_5S(double rij,double xii,double xij)
{
  return Slater_5S_6S(rij,xij,xii);
}

static double DSlater_1S_1S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-6.0 + 6.0*Power(E,2.0*rij*xii) - 12.0*rij*xii - 12.0*Power(rij,2.0)*Power(xii,2.0) - 
        7.0*Power(rij,3.0)*Power(xii,3.0) - 2.0*Power(rij,4.0)*Power(xii,4.0))/
      (6.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),3.0) + 
        Power(E,2.0*rij*xij)*Power(xij,4.0)*
         (-6.0*rij*Power(xii,3.0) - 2.0*Power(rij,2.0)*Power(xii,4.0) + Power(xij,2.0) + 
           2.0*rij*xii*Power(xij,2.0) + 
           Power(xii,2.0)*(-3.0 + 2.0*Power(rij,2.0)*Power(xij,2.0))) - 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (Power(xii,2.0)*(1.0 + 2.0*rij*xij + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           Power(xij,2.0)*(3.0 + 6.0*rij*xij + 2.0*Power(rij,2.0)*Power(xij,2.0))))/
      (Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),3.0))
    ;
  }
  return S;
}

static double DSlater_1S_2S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-120.0 + 120.0*Power(E,2.0*rij*xii) - 240.0*rij*xii - 
        240.0*Power(rij,2.0)*Power(xii,2.0) - 155.0*Power(rij,3.0)*Power(xii,3.0) - 
        70.0*Power(rij,4.0)*Power(xii,4.0) - 22.0*Power(rij,5.0)*Power(xii,5.0) - 
        4.0*Power(rij,6.0)*Power(xii,6.0))/(120.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (3.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),5.0) + 
        3.0*Power(E,2.0*rij*xij)*Power(xij,6.0)*
         (-4.0*Power(xii,4.0) - 8.0*rij*Power(xii,5.0) - 2.0*Power(rij,2.0)*Power(xii,6.0) - 
           10.0*rij*Power(xii,3.0)*Power(xij,2.0) + Power(xij,4.0) + 
           2.0*rij*xii*Power(xij,4.0) + 
           Power(xii,2.0)*Power(xij,2.0)*(-5.0 + 2.0*Power(rij,2.0)*Power(xij,2.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-6.0*Power(xii,2.0)*Power(xij,4.0)*
            (5.0 + 10.0*rij*xij + 13.0*Power(rij,2.0)*Power(xij,2.0) + 
              6.0*Power(rij,3.0)*Power(xij,3.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*Power(xij,6.0)*(21.0 + 42.0*rij*xij + 27.0*Power(rij,2.0)*Power(xij,2.0) + 
              8.0*Power(rij,3.0)*Power(xij,3.0) + Power(rij,4.0)*Power(xij,4.0)) - 
           Power(xii,6.0)*(3.0 + 6.0*rij*xij + 6.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,3.0)*Power(xij,3.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           3.0*Power(xii,4.0)*Power(xij,2.0)*
            (5.0 + 10.0*rij*xij + 10.0*Power(rij,2.0)*Power(xij,2.0) + 
              8.0*Power(rij,3.0)*Power(xij,3.0) + 2.0*Power(rij,4.0)*Power(xij,4.0))))/
      (3.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),5.0))
    ;
  }
  return S;
}

static double DSlater_1S_3S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-30240.0 + 30240.0*Power(E,2.0*rij*xii) - 60480.0*rij*xii - 
        60480.0*Power(rij,2.0)*Power(xii,2.0) - 40005.0*Power(rij,3.0)*Power(xii,3.0) - 
        19530.0*Power(rij,4.0)*Power(xii,4.0) - 7392.0*Power(rij,5.0)*Power(xii,5.0) - 
        2184.0*Power(rij,6.0)*Power(xii,6.0) - 480.0*Power(rij,7.0)*Power(xii,7.0) - 
        64.0*Power(rij,8.0)*Power(xii,8.0))/(30240.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (45.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),7.0) + 
        15.0*Power(E,2.0*rij*xij)*Power(xij,8.0)*
         (-30.0*rij*Power(xii,7.0) - 6.0*Power(rij,2.0)*Power(xii,8.0) - 
           126.0*rij*Power(xii,5.0)*Power(xij,2.0) - 
           42.0*rij*Power(xii,3.0)*Power(xij,4.0) + 3.0*Power(xij,6.0) + 
           6.0*rij*xii*Power(xij,6.0) + 
           7.0*Power(xii,4.0)*Power(xij,2.0)*(-9.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*Power(xii,2.0)*Power(xij,4.0)*(-7.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           Power(xii,6.0)*(15.0 + 14.0*Power(rij,2.0)*Power(xij,2.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-10.0*Power(xii,2.0)*Power(xij,8.0)*
            (135.0 + 270.0*rij*xij + 438.0*Power(rij,2.0)*Power(xij,2.0) + 
              306.0*Power(rij,3.0)*Power(xij,3.0) + 111.0*Power(rij,4.0)*Power(xij,4.0) + 
              22.0*Power(rij,5.0)*Power(xij,5.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           2.0*Power(xij,10.0)*(945.0 + 1890.0*rij*xij + 
              1470.0*Power(rij,2.0)*Power(xij,2.0) + 630.0*Power(rij,3.0)*Power(xij,3.0) + 
              165.0*Power(rij,4.0)*Power(xij,4.0) + 26.0*Power(rij,5.0)*Power(xij,5.0) + 
              2.0*Power(rij,6.0)*Power(xij,6.0)) - 
           Power(xii,10.0)*(45.0 + 90.0*rij*xij + 90.0*Power(rij,2.0)*Power(xij,2.0) + 
              60.0*Power(rij,3.0)*Power(xij,3.0) + 30.0*Power(rij,4.0)*Power(xij,4.0) + 
              12.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           5.0*Power(xii,8.0)*Power(xij,2.0)*
            (63.0 + 126.0*rij*xij + 126.0*Power(rij,2.0)*Power(xij,2.0) + 
              84.0*Power(rij,3.0)*Power(xij,3.0) + 42.0*Power(rij,4.0)*Power(xij,4.0) + 
              20.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           5.0*Power(xii,6.0)*Power(xij,4.0)*
            (189.0 + 378.0*rij*xij + 378.0*Power(rij,2.0)*Power(xij,2.0) + 
              240.0*Power(rij,3.0)*Power(xij,3.0) + 156.0*Power(rij,4.0)*Power(xij,4.0) + 
              56.0*Power(rij,5.0)*Power(xij,5.0) + 8.0*Power(rij,6.0)*Power(xij,6.0)) + 
           5.0*Power(xii,4.0)*Power(xij,6.0)*
            (315.0 + 630.0*rij*xij + 558.0*Power(rij,2.0)*Power(xij,2.0) + 
              528.0*Power(rij,3.0)*Power(xij,3.0) + 276.0*Power(rij,4.0)*Power(xij,4.0) + 
              72.0*Power(rij,5.0)*Power(xij,5.0) + 8.0*Power(rij,6.0)*Power(xij,6.0))))/
      (45.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),7.0))
    ;
  }
  return S;
}

static double DSlater_1S_4S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-362880.0 + 362880.0*Power(E,2.0*rij*xii) - 725760.0*rij*xii - 
        725760.0*Power(rij,2.0)*Power(xii,2.0) - 482895.0*Power(rij,3.0)*Power(xii,3.0) - 
        240030.0*Power(rij,4.0)*Power(xii,4.0) - 94689.0*Power(rij,5.0)*Power(xii,5.0) - 
        30618.0*Power(rij,6.0)*Power(xii,6.0) - 8208.0*Power(rij,7.0)*Power(xii,7.0) - 
        1800.0*Power(rij,8.0)*Power(xii,8.0) - 304.0*Power(rij,9.0)*Power(xii,9.0) - 
        32.0*Power(rij,10.0)*Power(xii,10.0))/(362880.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (315.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),9.0) + 
        315.0*Power(E,2.0*rij*xij)*Power(xij,10.0)*
         (-12.0*rij*Power(xii,9.0) - 2.0*Power(rij,2.0)*Power(xii,10.0) - 
           51.0*Power(xii,6.0)*Power(xij,2.0) - 102.0*rij*Power(xii,7.0)*Power(xij,2.0) - 
           126.0*rij*Power(xii,5.0)*Power(xij,4.0) - 
           18.0*rij*Power(xii,3.0)*Power(xij,6.0) + Power(xij,8.0) + 
           2.0*rij*xii*Power(xij,8.0) + 
           Power(xii,2.0)*Power(xij,6.0)*(-9.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           6.0*Power(xii,8.0)*(1.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*Power(xii,4.0)*Power(xij,4.0)*(-21.0 + 4.0*Power(rij,2.0)*Power(xij,2.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-42.0*Power(xii,10.0)*Power(xij,4.0)*
            (270.0 + 540.0*rij*xij + 540.0*Power(rij,2.0)*Power(xij,2.0) + 
              360.0*Power(rij,3.0)*Power(xij,3.0) + 180.0*Power(rij,4.0)*Power(xij,4.0) + 
              69.0*Power(rij,5.0)*Power(xij,5.0) + 28.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0) + Power(rij,8.0)*Power(xij,8.0)) + 
           70.0*Power(xii,8.0)*Power(xij,6.0)*
            (378.0 + 756.0*rij*xij + 756.0*Power(rij,2.0)*Power(xij,2.0) + 
              510.0*Power(rij,3.0)*Power(xij,3.0) + 228.0*Power(rij,4.0)*Power(xij,4.0) + 
              111.0*Power(rij,5.0)*Power(xij,5.0) + 44.0*Power(rij,6.0)*Power(xij,6.0) + 
              10.0*Power(rij,7.0)*Power(xij,7.0) + Power(rij,8.0)*Power(xij,8.0)) - 
           70.0*Power(xii,6.0)*Power(xij,8.0)*
            (567.0 + 1134.0*rij*xij + 1179.0*Power(rij,2.0)*Power(xij,2.0) + 
              630.0*Power(rij,3.0)*Power(xij,3.0) + 387.0*Power(rij,4.0)*Power(xij,4.0) + 
              204.0*Power(rij,5.0)*Power(xij,5.0) + 66.0*Power(rij,6.0)*Power(xij,6.0) + 
              12.0*Power(rij,7.0)*Power(xij,7.0) + Power(rij,8.0)*Power(xij,8.0)) + 
           42.0*Power(xii,4.0)*Power(xij,10.0)*
            (990.0 + 1980.0*rij*xij + 1305.0*Power(rij,2.0)*Power(xij,2.0) + 
              1230.0*Power(rij,3.0)*Power(xij,3.0) + 885.0*Power(rij,4.0)*Power(xij,4.0) + 
              372.0*Power(rij,5.0)*Power(xij,5.0) + 94.0*Power(rij,6.0)*Power(xij,6.0) + 
              14.0*Power(rij,7.0)*Power(xij,7.0) + Power(rij,8.0)*Power(xij,8.0)) - 
           Power(xii,14.0)*(315.0 + 630.0*rij*xij + 630.0*Power(rij,2.0)*Power(xij,2.0) + 
              420.0*Power(rij,3.0)*Power(xij,3.0) + 210.0*Power(rij,4.0)*Power(xij,4.0) + 
              84.0*Power(rij,5.0)*Power(xij,5.0) + 28.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           7.0*Power(xii,12.0)*Power(xij,2.0)*
            (405.0 + 810.0*rij*xij + 810.0*Power(rij,2.0)*Power(xij,2.0) + 
              540.0*Power(rij,3.0)*Power(xij,3.0) + 270.0*Power(rij,4.0)*Power(xij,4.0) + 
              108.0*Power(rij,5.0)*Power(xij,5.0) + 36.0*Power(rij,6.0)*Power(xij,6.0) + 
              12.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           7.0*Power(xii,2.0)*Power(xij,12.0)*
            (1485.0 + 2970.0*rij*xij + 8640.0*Power(rij,2.0)*Power(xij,2.0) + 
              8280.0*Power(rij,3.0)*Power(xij,3.0) + 4140.0*Power(rij,4.0)*Power(xij,4.0) + 
              1278.0*Power(rij,5.0)*Power(xij,5.0) + 256.0*Power(rij,6.0)*Power(xij,6.0) + 
              32.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xij,14.0)*(31185.0 + 62370.0*rij*xij + 
              52920.0*Power(rij,2.0)*Power(xij,2.0) + 
              26460.0*Power(rij,3.0)*Power(xij,3.0) + 8820.0*Power(rij,4.0)*Power(xij,4.0) + 
              2058.0*Power(rij,5.0)*Power(xij,5.0) + 336.0*Power(rij,6.0)*Power(xij,6.0) + 
              36.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0))))/
      (315.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),9.0))
    ;
  }
  return S;
}

static double DSlater_1S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-399168000.0 + 399168000.0*Power(E,2.0*rij*xii) - 798336000.0*rij*xii - 
        798336000.0*Power(rij,2.0)*Power(xii,2.0) - 
        531964125.0*Power(rij,3.0)*Power(xii,3.0) - 
        265592250.0*Power(rij,4.0)*Power(xii,4.0) - 
        105862680.0*Power(rij,5.0)*Power(xii,5.0) - 
        35010360.0*Power(rij,6.0)*Power(xii,6.0) - 9836640.0*Power(rij,7.0)*Power(xii,7.0) - 
        2376000.0*Power(rij,8.0)*Power(xii,8.0) - 492800.0*Power(rij,9.0)*Power(xii,9.0) - 
        85888.0*Power(rij,10.0)*Power(xii,10.0) - 11776.0*Power(rij,11.0)*Power(xii,11.0) - 
        1024.0*Power(rij,12.0)*Power(xii,12.0))/
      (3.99168e8*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (14175.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),11.0) + 
        2835.0*Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-70.0*rij*Power(xii,11.0) - 10.0*Power(rij,2.0)*Power(xii,12.0) - 
           990.0*rij*Power(xii,9.0)*Power(xij,2.0) - 
           2508.0*rij*Power(xii,7.0)*Power(xij,4.0) - 
           1452.0*rij*Power(xii,5.0)*Power(xij,6.0) - 
           110.0*rij*Power(xii,3.0)*Power(xij,8.0) + 5.0*Power(xij,10.0) + 
           10.0*rij*xii*Power(xij,10.0) + 
           66.0*Power(xii,6.0)*Power(xij,4.0)*(-19.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           5.0*Power(xii,2.0)*Power(xij,8.0)*(-11.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           33.0*Power(xii,8.0)*Power(xij,2.0)*(15.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) + 
           22.0*Power(xii,4.0)*Power(xij,6.0)*(-33.0 + 5.0*Power(rij,2.0)*Power(xij,2.0)) - 
           5.0*Power(xii,10.0)*(7.0 + 22.0*Power(rij,2.0)*Power(xij,2.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-(Power(xii,18.0)*(14175.0 + 28350.0*rij*xij + 
                28350.0*Power(rij,2.0)*Power(xij,2.0) + 
                18900.0*Power(rij,3.0)*Power(xij,3.0) + 
                9450.0*Power(rij,4.0)*Power(xij,4.0) + 
                3780.0*Power(rij,5.0)*Power(xij,5.0) + 
                1260.0*Power(rij,6.0)*Power(xij,6.0) + 
                360.0*Power(rij,7.0)*Power(xij,7.0) + 90.0*Power(rij,8.0)*Power(xij,8.0) + 
                20.0*Power(rij,9.0)*Power(xij,9.0) + 4.0*Power(rij,10.0)*Power(xij,10.0))) + 
           9.0*Power(xii,16.0)*Power(xij,2.0)*
            (17325.0 + 34650.0*rij*xij + 34650.0*Power(rij,2.0)*Power(xij,2.0) + 
              23100.0*Power(rij,3.0)*Power(xij,3.0) + 
              11550.0*Power(rij,4.0)*Power(xij,4.0) + 
              4620.0*Power(rij,5.0)*Power(xij,5.0) + 1540.0*Power(rij,6.0)*Power(xij,6.0) + 
              440.0*Power(rij,7.0)*Power(xij,7.0) + 110.0*Power(rij,8.0)*Power(xij,8.0) + 
              28.0*Power(rij,9.0)*Power(xij,9.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           126.0*Power(xii,10.0)*Power(xij,8.0)*
            (37125.0 + 74250.0*rij*xij + 74250.0*Power(rij,2.0)*Power(xij,2.0) + 
              49350.0*Power(rij,3.0)*Power(xij,3.0) + 
              25575.0*Power(rij,4.0)*Power(xij,4.0) + 
              9078.0*Power(rij,5.0)*Power(xij,5.0) + 3106.0*Power(rij,6.0)*Power(xij,6.0) + 
              1136.0*Power(rij,7.0)*Power(xij,7.0) + 314.0*Power(rij,8.0)*Power(xij,8.0) + 
              52.0*Power(rij,9.0)*Power(xij,9.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) + 
           126.0*Power(xii,8.0)*Power(xij,10.0)*
            (51975.0 + 103950.0*rij*xij + 102600.0*Power(rij,2.0)*Power(xij,2.0) + 
              74850.0*Power(rij,3.0)*Power(xij,3.0) + 
              31125.0*Power(rij,4.0)*Power(xij,4.0) + 
              11730.0*Power(rij,5.0)*Power(xij,5.0) + 
              5150.0*Power(rij,6.0)*Power(xij,6.0) + 1840.0*Power(rij,7.0)*Power(xij,7.0) + 
              430.0*Power(rij,8.0)*Power(xij,8.0) + 60.0*Power(rij,9.0)*Power(xij,9.0) + 
              4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           9.0*Power(xii,2.0)*Power(xij,16.0)*
            (-135135.0 - 270270.0*rij*xij + 228690.0*Power(rij,2.0)*Power(xij,2.0) + 
              471240.0*Power(rij,3.0)*Power(xij,3.0) + 
              318780.0*Power(rij,4.0)*Power(xij,4.0) + 
              127512.0*Power(rij,5.0)*Power(xij,5.0) + 
              34664.0*Power(rij,6.0)*Power(xij,6.0) + 
              6712.0*Power(rij,7.0)*Power(xij,7.0) + 922.0*Power(rij,8.0)*Power(xij,8.0) + 
              84.0*Power(rij,9.0)*Power(xij,9.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) + 
           Power(xij,18.0)*(2837835.0 + 5675670.0*rij*xij + 
              5051970.0*Power(rij,2.0)*Power(xij,2.0) + 
              2744280.0*Power(rij,3.0)*Power(xij,3.0) + 
              1031940.0*Power(rij,4.0)*Power(xij,4.0) + 
              285768.0*Power(rij,5.0)*Power(xij,5.0) + 
              59976.0*Power(rij,6.0)*Power(xij,6.0) + 
              9576.0*Power(rij,7.0)*Power(xij,7.0) + 1134.0*Power(rij,8.0)*Power(xij,8.0) + 
              92.0*Power(rij,9.0)*Power(xij,9.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           9.0*Power(xii,14.0)*Power(xij,4.0)*
            (86625.0 + 173250.0*rij*xij + 173250.0*Power(rij,2.0)*Power(xij,2.0) + 
              115500.0*Power(rij,3.0)*Power(xij,3.0) + 
              57750.0*Power(rij,4.0)*Power(xij,4.0) + 
              23100.0*Power(rij,5.0)*Power(xij,5.0) + 
              7700.0*Power(rij,6.0)*Power(xij,6.0) + 2128.0*Power(rij,7.0)*Power(xij,7.0) + 
              616.0*Power(rij,8.0)*Power(xij,8.0) + 144.0*Power(rij,9.0)*Power(xij,9.0) + 
              16.0*Power(rij,10.0)*Power(xij,10.0)) + 
           21.0*Power(xii,12.0)*Power(xij,6.0)*
            (111375.0 + 222750.0*rij*xij + 222750.0*Power(rij,2.0)*Power(xij,2.0) + 
              148500.0*Power(rij,3.0)*Power(xij,3.0) + 
              74250.0*Power(rij,4.0)*Power(xij,4.0) + 
              29988.0*Power(rij,5.0)*Power(xij,5.0) + 
              9276.0*Power(rij,6.0)*Power(xij,6.0) + 2928.0*Power(rij,7.0)*Power(xij,7.0) + 
              888.0*Power(rij,8.0)*Power(xij,8.0) + 176.0*Power(rij,9.0)*Power(xij,9.0) + 
              16.0*Power(rij,10.0)*Power(xij,10.0)) - 
           21.0*Power(xii,6.0)*Power(xij,12.0)*
            (307125.0 + 614250.0*rij*xij + 733050.0*Power(rij,2.0)*Power(xij,2.0) + 
              350100.0*Power(rij,3.0)*Power(xij,3.0) + 
              151290.0*Power(rij,4.0)*Power(xij,4.0) + 
              85860.0*Power(rij,5.0)*Power(xij,5.0) + 
              39180.0*Power(rij,6.0)*Power(xij,6.0) + 
              11760.0*Power(rij,7.0)*Power(xij,7.0) + 
              2280.0*Power(rij,8.0)*Power(xij,8.0) + 272.0*Power(rij,9.0)*Power(xij,9.0) + 
              16.0*Power(rij,10.0)*Power(xij,10.0)) + 
           9.0*Power(xii,4.0)*Power(xij,14.0)*
            (675675.0 + 1351350.0*rij*xij + 602910.0*Power(rij,2.0)*Power(xij,2.0) + 
              374220.0*Power(rij,3.0)*Power(xij,3.0) + 
              353430.0*Power(rij,4.0)*Power(xij,4.0) + 
              207900.0*Power(rij,5.0)*Power(xij,5.0) + 
              75460.0*Power(rij,6.0)*Power(xij,6.0) + 
              18128.0*Power(rij,7.0)*Power(xij,7.0) + 2936.0*Power(rij,8.0)*Power(xij,8.0) + 
              304.0*Power(rij,9.0)*Power(xij,9.0) + 16.0*Power(rij,10.0)*Power(xij,10.0))))/
      (14175.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),11.0))
    ;
  }
  return S;
}

static double DSlater_1S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-37362124800.0 + 37362124800.0*Power(E,2.0*rij*xii) - 74724249600.0*rij*xii - 
        74724249600.0*Power(rij,2.0)*Power(xii,2.0) - 
        49810085325.0*Power(rij,3.0)*Power(xii,3.0) - 
        24895921050.0*Power(rij,4.0)*Power(xii,4.0) - 
        9949449510.0*Power(rij,5.0)*Power(xii,5.0) - 
        3309726420.0*Power(rij,6.0)*Power(xii,6.0) - 
        941466240.0*Power(rij,7.0)*Power(xii,7.0) - 
        233204400.0*Power(rij,8.0)*Power(xii,8.0) - 
        50862240.0*Power(rij,9.0)*Power(xii,9.0) - 
        9801792.0*Power(rij,10.0)*Power(xii,10.0) - 
        1657344.0*Power(rij,11.0)*Power(xii,11.0) - 
        239616.0*Power(rij,12.0)*Power(xii,12.0) - 27648.0*Power(rij,13.0)*Power(xii,13.0) - 
        2048.0*Power(rij,14.0)*Power(xii,14.0))/
      (3.73621248e10*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (467775.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),13.0) + 
        155925.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-48.0*rij*Power(xii,13.0) - 6.0*Power(rij,2.0)*Power(xii,14.0) - 
           1014.0*rij*Power(xii,11.0)*Power(xij,2.0) - 
           2145.0*Power(xii,8.0)*Power(xij,4.0) - 
           4290.0*rij*Power(xii,9.0)*Power(xij,4.0) - 
           5148.0*rij*Power(xii,7.0)*Power(xij,6.0) - 
           1716.0*rij*Power(xii,5.0)*Power(xij,8.0) - 
           78.0*rij*Power(xii,3.0)*Power(xij,10.0) + 3.0*Power(xij,12.0) + 
           6.0*rij*xii*Power(xij,12.0) + 
           286.0*Power(xii,6.0)*Power(xij,6.0)*(-9.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*Power(xii,2.0)*Power(xij,10.0)*(-13.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           26.0*Power(xii,4.0)*Power(xij,8.0)*(-33.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) - 
           8.0*Power(xii,12.0)*(3.0 + 13.0*Power(rij,2.0)*Power(xij,2.0)) - 
           13.0*Power(xii,10.0)*Power(xij,2.0)*(39.0 + 22.0*Power(rij,2.0)*Power(xij,2.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,4.0)*
         (-110.0*Power(xii,18.0)*Power(xij,4.0)*
            (331695.0 + 663390.0*rij*xij + 663390.0*Power(rij,2.0)*Power(xij,2.0) + 
              442260.0*Power(rij,3.0)*Power(xij,3.0) + 
              221130.0*Power(rij,4.0)*Power(xij,4.0) + 
              88452.0*Power(rij,5.0)*Power(xij,5.0) + 
              29484.0*Power(rij,6.0)*Power(xij,6.0) + 
              8424.0*Power(rij,7.0)*Power(xij,7.0) + 2106.0*Power(rij,8.0)*Power(xij,8.0) + 
              456.0*Power(rij,9.0)*Power(xij,9.0) + 102.0*Power(rij,10.0)*Power(xij,10.0) + 
              20.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) + 
           330.0*Power(xii,16.0)*Power(xij,6.0)*
            (405405.0 + 810810.0*rij*xij + 810810.0*Power(rij,2.0)*Power(xij,2.0) + 
              540540.0*Power(rij,3.0)*Power(xij,3.0) + 
              270270.0*Power(rij,4.0)*Power(xij,4.0) + 
              108108.0*Power(rij,5.0)*Power(xij,5.0) + 
              36036.0*Power(rij,6.0)*Power(xij,6.0) + 
              10368.0*Power(rij,7.0)*Power(xij,7.0) + 
              2466.0*Power(rij,8.0)*Power(xij,8.0) + 576.0*Power(rij,9.0)*Power(xij,9.0) + 
              138.0*Power(rij,10.0)*Power(xij,10.0) + 
              24.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) - 
           330.0*Power(xii,6.0)*Power(xij,16.0)*
            (1584765.0 + 3169530.0*rij*xij + 5061420.0*Power(rij,2.0)*Power(xij,2.0) + 
              2653560.0*Power(rij,3.0)*Power(xij,3.0) + 
              786240.0*Power(rij,4.0)*Power(xij,4.0) + 
              296478.0*Power(rij,5.0)*Power(xij,5.0) + 
              158886.0*Power(rij,6.0)*Power(xij,6.0) + 
              65988.0*Power(rij,7.0)*Power(xij,7.0) + 
              18681.0*Power(rij,8.0)*Power(xij,8.0) + 
              3666.0*Power(rij,9.0)*Power(xij,9.0) + 
              498.0*Power(rij,10.0)*Power(xij,10.0) + 
              44.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) + 
           110.0*Power(xii,4.0)*Power(xij,18.0)*
            (6081075.0 + 12162150.0*rij*xij + 4864860.0*Power(rij,2.0)*Power(xij,2.0) + 
              810810.0*Power(rij,3.0)*Power(xij,3.0) + 
              810810.0*Power(rij,4.0)*Power(xij,4.0) + 
              810810.0*Power(rij,5.0)*Power(xij,5.0) + 
              417690.0*Power(rij,6.0)*Power(xij,6.0) + 
              136188.0*Power(rij,7.0)*Power(xij,7.0) + 
              31023.0*Power(rij,8.0)*Power(xij,8.0) + 
              5118.0*Power(rij,9.0)*Power(xij,9.0) + 
              606.0*Power(rij,10.0)*Power(xij,10.0) + 
              48.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) - 
           22.0*Power(xii,2.0)*Power(xij,20.0)*
            (-12162150.0 - 24324300.0*rij*xij - 
              10135125.0*Power(rij,2.0)*Power(xij,2.0) + 
              4054050.0*Power(rij,3.0)*Power(xij,3.0) + 
              6081075.0*Power(rij,4.0)*Power(xij,4.0) + 
              3243240.0*Power(rij,5.0)*Power(xij,5.0) + 
              1099980.0*Power(rij,6.0)*Power(xij,6.0) + 
              268920.0*Power(rij,7.0)*Power(xij,7.0) + 
              49590.0*Power(rij,8.0)*Power(xij,8.0) + 
              6960.0*Power(rij,9.0)*Power(xij,9.0) + 
              726.0*Power(rij,10.0)*Power(xij,10.0) + 
              52.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) + 
           2.0*Power(xij,22.0)*(85135050.0 + 170270100.0*rij*xij + 
              156080925.0*Power(rij,2.0)*Power(xij,2.0) + 
              89189100.0*Power(rij,3.0)*Power(xij,3.0) + 
              36018675.0*Power(rij,4.0)*Power(xij,4.0) + 
              10977120.0*Power(rij,5.0)*Power(xij,5.0) + 
              2619540.0*Power(rij,6.0)*Power(xij,6.0) + 
              498960.0*Power(rij,7.0)*Power(xij,7.0) + 
              76230.0*Power(rij,8.0)*Power(xij,8.0) + 
              9240.0*Power(rij,9.0)*Power(xij,9.0) + 
              858.0*Power(rij,10.0)*Power(xij,10.0) + 
              56.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) - 
           Power(xii,22.0)*(467775.0 + 935550.0*rij*xij + 
              935550.0*Power(rij,2.0)*Power(xij,2.0) + 
              623700.0*Power(rij,3.0)*Power(xij,3.0) + 
              311850.0*Power(rij,4.0)*Power(xij,4.0) + 
              124740.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              11880.0*Power(rij,7.0)*Power(xij,7.0) + 
              2970.0*Power(rij,8.0)*Power(xij,8.0) + 660.0*Power(rij,9.0)*Power(xij,9.0) + 
              132.0*Power(rij,10.0)*Power(xij,10.0) + 
              24.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) + 
           11.0*Power(xii,20.0)*Power(xij,2.0)*
            (552825.0 + 1105650.0*rij*xij + 1105650.0*Power(rij,2.0)*Power(xij,2.0) + 
              737100.0*Power(rij,3.0)*Power(xij,3.0) + 
              368550.0*Power(rij,4.0)*Power(xij,4.0) + 
              147420.0*Power(rij,5.0)*Power(xij,5.0) + 
              49140.0*Power(rij,6.0)*Power(xij,6.0) + 
              14040.0*Power(rij,7.0)*Power(xij,7.0) + 
              3510.0*Power(rij,8.0)*Power(xij,8.0) + 780.0*Power(rij,9.0)*Power(xij,9.0) + 
              156.0*Power(rij,10.0)*Power(xij,10.0) + 
              32.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           462.0*Power(xii,10.0)*Power(xij,12.0)*
            (1737450.0 + 3474900.0*rij*xij + 3489075.0*Power(rij,2.0)*Power(xij,2.0) + 
              2239650.0*Power(rij,3.0)*Power(xij,3.0) + 
              1248075.0*Power(rij,4.0)*Power(xij,4.0) + 
              468180.0*Power(rij,5.0)*Power(xij,5.0) + 
              129960.0*Power(rij,6.0)*Power(xij,6.0) + 
              38880.0*Power(rij,7.0)*Power(xij,7.0) + 
              12960.0*Power(rij,8.0)*Power(xij,8.0) + 
              3480.0*Power(rij,9.0)*Power(xij,9.0) + 
              636.0*Power(rij,10.0)*Power(xij,10.0) + 
              72.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) + 
           330.0*Power(xii,8.0)*Power(xij,14.0)*
            (2443770.0 + 4887540.0*rij*xij + 4457565.0*Power(rij,2.0)*Power(xij,2.0) + 
              3749760.0*Power(rij,3.0)*Power(xij,3.0) + 
              1715175.0*Power(rij,4.0)*Power(xij,4.0) + 
              510804.0*Power(rij,5.0)*Power(xij,5.0) + 
              164808.0*Power(rij,6.0)*Power(xij,6.0) + 
              65808.0*Power(rij,7.0)*Power(xij,7.0) + 
              21912.0*Power(rij,8.0)*Power(xij,8.0) + 
              5112.0*Power(rij,9.0)*Power(xij,9.0) + 
              804.0*Power(rij,10.0)*Power(xij,10.0) + 
              80.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           165.0*Power(xii,14.0)*Power(xij,8.0)*
            (2027025.0 + 4054050.0*rij*xij + 4054050.0*Power(rij,2.0)*Power(xij,2.0) + 
              2702700.0*Power(rij,3.0)*Power(xij,3.0) + 
              1351350.0*Power(rij,4.0)*Power(xij,4.0) + 
              539280.0*Power(rij,5.0)*Power(xij,5.0) + 
              183960.0*Power(rij,6.0)*Power(xij,6.0) + 
              49392.0*Power(rij,7.0)*Power(xij,7.0) + 
              12012.0*Power(rij,8.0)*Power(xij,8.0) + 
              3192.0*Power(rij,9.0)*Power(xij,9.0) + 
              744.0*Power(rij,10.0)*Power(xij,10.0) + 
              112.0*Power(rij,11.0)*Power(xij,11.0) + 8.0*Power(rij,12.0)*Power(xij,12.0)) + 
           231.0*Power(xii,12.0)*Power(xij,10.0)*
            (2606175.0 + 5212350.0*rij*xij + 5212350.0*Power(rij,2.0)*Power(xij,2.0) + 
              3477600.0*Power(rij,3.0)*Power(xij,3.0) + 
              1718550.0*Power(rij,4.0)*Power(xij,4.0) + 
              724320.0*Power(rij,5.0)*Power(xij,5.0) + 
              226440.0*Power(rij,6.0)*Power(xij,6.0) + 
              58320.0*Power(rij,7.0)*Power(xij,7.0) + 
              16500.0*Power(rij,8.0)*Power(xij,8.0) + 4680.0*Power(rij,9.0)*Power(xij,9.0) + 
              984.0*Power(rij,10.0)*Power(xij,10.0) + 
              128.0*Power(rij,11.0)*Power(xij,11.0) + 8.0*Power(rij,12.0)*Power(xij,12.0))))/
      (467775.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),13.0))
    ;
  }
  return S;
}

static double DSlater_2S_2S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-20160.0 + 20160.0*Power(E,2.0*rij*xii) - 40320.0*rij*xii - 
        40320.0*Power(rij,2.0)*Power(xii,2.0) - 26355.0*Power(rij,3.0)*Power(xii,3.0) - 
        12390.0*Power(rij,4.0)*Power(xii,4.0) - 4368.0*Power(rij,5.0)*Power(xii,5.0) - 
        1176.0*Power(rij,6.0)*Power(xii,6.0) - 240.0*Power(rij,7.0)*Power(xii,7.0) - 
        32.0*Power(rij,8.0)*Power(xii,8.0))/(20160.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (3.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),7.0) + 
        Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (-21.0*Power(xii,4.0)*Power(xij,4.0)*
            (3.0 + 6.0*rij*xij + 10.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,3.0)*Power(xij,3.0)) + 
           Power(xii,2.0)*Power(xij,6.0)*
            (195.0 + 390.0*rij*xij + 78.0*Power(rij,2.0)*Power(xij,2.0) - 
              14.0*Power(rij,3.0)*Power(xij,3.0) - 4.0*Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*Power(xij,8.0)*(45.0 + 90.0*rij*xij + 48.0*Power(rij,2.0)*Power(xij,2.0) + 
              11.0*Power(rij,3.0)*Power(xij,3.0) + Power(rij,4.0)*Power(xij,4.0)) - 
           Power(xii,8.0)*(3.0 + 6.0*rij*xij + 6.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,3.0)*Power(xij,3.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,6.0)*Power(xij,2.0)*
            (21.0 + 42.0*rij*xij + 42.0*Power(rij,2.0)*Power(xij,2.0) + 
              38.0*Power(rij,3.0)*Power(xij,3.0) + 4.0*Power(rij,4.0)*Power(xij,4.0))) + 
        Power(E,2.0*rij*xij)*Power(xij,6.0)*
         (-22.0*Power(rij,3.0)*Power(xii,11.0) - 2.0*Power(rij,4.0)*Power(xii,12.0) + 
           3.0*Power(xij,8.0) + 6.0*rij*xii*Power(xij,8.0) + 
           4.0*Power(rij,2.0)*Power(xii,10.0)*(-24.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*Power(xii,2.0)*Power(xij,6.0)*(-7.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           2.0*rij*Power(xii,9.0)*(-90.0 + 7.0*Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*rij*Power(xii,7.0)*Power(xij,2.0)*(-65.0 + 7.0*Power(rij,2.0)*Power(xij,2.0)) - 
           6.0*Power(xii,8.0)*(15.0 + 13.0*Power(rij,2.0)*Power(xij,2.0)) + 
           Power(xii,4.0)*Power(xij,4.0)*
            (63.0 - 42.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,5.0)*(126.0*rij*Power(xij,4.0) - 38.0*Power(rij,3.0)*Power(xij,6.0)) + 
           Power(xii,6.0)*(-195.0*Power(xij,2.0) + 210.0*Power(rij,2.0)*Power(xij,4.0) - 
              4.0*Power(rij,4.0)*Power(xij,6.0)) + 
           Power(xii,3.0)*(-42.0*rij*Power(xij,6.0) + 4.0*Power(rij,3.0)*Power(xij,8.0))))/
      (3.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),7.0))
    ;
  }
  return S;
}

static double DSlater_2S_3S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-544320.0 + 544320.0*Power(E,2.0*rij*xii) - 1088640.0*rij*xii - 
        1088640.0*Power(rij,2.0)*Power(xii,2.0) - 719145.0*Power(rij,3.0)*Power(xii,3.0) - 
        349650.0*Power(rij,4.0)*Power(xii,4.0) - 132111.0*Power(rij,5.0)*Power(xii,5.0) - 
        39942.0*Power(rij,6.0)*Power(xii,6.0) - 9792.0*Power(rij,7.0)*Power(xii,7.0) - 
        1944.0*Power(rij,8.0)*Power(xii,8.0) - 304.0*Power(rij,9.0)*Power(xii,9.0) - 
        32.0*Power(rij,10.0)*Power(xii,10.0))/(544320.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (45.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),9.0) + 
        Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (180.0*Power(xii,6.0)*Power(xij,6.0)*
            (21.0 + 42.0*rij*xij + 27.0*Power(rij,2.0)*Power(xij,2.0) + 
              54.0*Power(rij,3.0)*Power(xij,3.0) + 19.0*Power(rij,4.0)*Power(xij,4.0) + 
              2.0*Power(rij,5.0)*Power(xij,5.0)) + 
           Power(xii,2.0)*Power(xij,10.0)*
            (21615.0 + 43230.0*rij*xij + 17850.0*Power(rij,2.0)*Power(xij,2.0) + 
              1580.0*Power(rij,3.0)*Power(xij,3.0) - 650.0*Power(rij,4.0)*Power(xij,4.0) - 
              188.0*Power(rij,5.0)*Power(xij,5.0) - 16.0*Power(rij,6.0)*Power(xij,6.0)) + 
           4.0*Power(xij,12.0)*(1485.0 + 2970.0*rij*xij + 
              2025.0*Power(rij,2.0)*Power(xij,2.0) + 720.0*Power(rij,3.0)*Power(xij,3.0) + 
              150.0*Power(rij,4.0)*Power(xij,4.0) + 18.0*Power(rij,5.0)*Power(xij,5.0) + 
              Power(rij,6.0)*Power(xij,6.0)) - 
           20.0*Power(xii,8.0)*Power(xij,4.0)*
            (81.0 + 162.0*rij*xij + 162.0*Power(rij,2.0)*Power(xij,2.0) + 
              94.0*Power(rij,3.0)*Power(xij,3.0) + 89.0*Power(rij,4.0)*Power(xij,4.0) + 
              20.0*Power(rij,5.0)*Power(xij,5.0) + Power(rij,6.0)*Power(xij,6.0)) + 
           5.0*Power(xii,4.0)*Power(xij,8.0)*
            (-639.0 - 1278.0*rij*xij - 5658.0*Power(rij,2.0)*Power(xij,2.0) - 
              2556.0*Power(rij,3.0)*Power(xij,3.0) - 366.0*Power(rij,4.0)*Power(xij,4.0) + 
              4.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           Power(xii,12.0)*(45.0 + 90.0*rij*xij + 90.0*Power(rij,2.0)*Power(xij,2.0) + 
              60.0*Power(rij,3.0)*Power(xij,3.0) + 30.0*Power(rij,4.0)*Power(xij,4.0) + 
              12.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           Power(xii,10.0)*Power(xij,2.0)*
            (405.0 + 810.0*rij*xij + 810.0*Power(rij,2.0)*Power(xij,2.0) + 
              540.0*Power(rij,3.0)*Power(xij,3.0) + 270.0*Power(rij,4.0)*Power(xij,4.0) + 
              148.0*Power(rij,5.0)*Power(xij,5.0) + 16.0*Power(rij,6.0)*Power(xij,6.0))) + 
        5.0*Power(E,2.0*rij*xij)*Power(xij,8.0)*
         (-84.0*Power(rij,3.0)*Power(xii,13.0) - 6.0*Power(rij,4.0)*Power(xii,14.0) + 
           9.0*Power(xij,10.0) + 18.0*rij*xii*Power(xij,10.0) - 
           72.0*rij*Power(xii,7.0)*Power(xij,4.0)*(54.0 + Power(rij,2.0)*Power(xij,2.0)) - 
           2.0*Power(rij,2.0)*Power(xii,12.0)*(225.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*rij*Power(xii,3.0)*Power(xij,8.0)*(-27.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           9.0*Power(xii,2.0)*Power(xij,8.0)*(-9.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*rij*Power(xii,9.0)*Power(xij,2.0)*
            (-1063.0 + 84.0*Power(rij,2.0)*Power(xij,2.0)) - 
           2.0*rij*Power(xii,11.0)*(495.0 + 98.0*Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*Power(xii,4.0)*Power(xij,6.0)*
            (54.0 - 27.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*Power(xii,6.0)*Power(xij,4.0)*
            (-972.0 + 702.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           9.0*Power(xii,10.0)*(-55.0 - 222.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,4.0)*Power(xij,4.0)) - 
           3.0*Power(xii,8.0)*Power(xij,2.0)*
            (1063.0 - 396.0*Power(rij,2.0)*Power(xij,2.0) + 
              12.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,5.0)*(648.0*rij*Power(xij,6.0) - 164.0*Power(rij,3.0)*Power(xij,8.0))))/
      (45.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),9.0))
    ;
  }
  return S;
}

static double DSlater_2S_4S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-159667200.0 + 159667200.0*Power(E,2.0*rij*xii) - 319334400.0*rij*xii - 
        319334400.0*Power(rij,2.0)*Power(xii,2.0) - 
        212109975.0*Power(rij,3.0)*Power(xii,3.0) - 
        104885550.0*Power(rij,4.0)*Power(xii,4.0) - 
        40997880.0*Power(rij,5.0)*Power(xii,5.0) - 
        13111560.0*Power(rij,6.0)*Power(xii,6.0) - 3496680.0*Power(rij,7.0)*Power(xii,7.0) - 
        784080.0*Power(rij,8.0)*Power(xii,8.0) - 147840.0*Power(rij,9.0)*Power(xii,9.0) - 
        23232.0*Power(rij,10.0)*Power(xii,10.0) - 2944.0*Power(rij,11.0)*Power(xii,11.0) - 
        256.0*Power(rij,12.0)*Power(xii,12.0))/
      (1.596672e8*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (315.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),11.0) + 
        Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (-770.0*Power(xii,8.0)*Power(xij,8.0)*
            (135.0 + 270.0*rij*xij + 315.0*Power(rij,2.0)*Power(xij,2.0) + 
              45.0*Power(rij,3.0)*Power(xij,3.0) + 135.0*Power(rij,4.0)*Power(xij,4.0) + 
              72.0*Power(rij,5.0)*Power(xij,5.0) + 14.0*Power(rij,6.0)*Power(xij,6.0) + 
              Power(rij,7.0)*Power(xij,7.0)) + 
           Power(xii,2.0)*Power(xij,14.0)*
            (765765.0 + 1531530.0*rij*xij + 866250.0*Power(rij,2.0)*Power(xij,2.0) + 
              210210.0*Power(rij,3.0)*Power(xij,3.0) + 
              11550.0*Power(rij,4.0)*Power(xij,4.0) - 
              6468.0*Power(rij,5.0)*Power(xij,5.0) - 
              1876.0*Power(rij,6.0)*Power(xij,6.0) - 230.0*Power(rij,7.0)*Power(xij,7.0) - 
              12.0*Power(rij,8.0)*Power(xij,8.0)) - 
           Power(xii,16.0)*(315.0 + 630.0*rij*xij + 630.0*Power(rij,2.0)*Power(xij,2.0) + 
              420.0*Power(rij,3.0)*Power(xij,3.0) + 210.0*Power(rij,4.0)*Power(xij,4.0) + 
              84.0*Power(rij,5.0)*Power(xij,5.0) + 28.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xij,16.0)*(135135.0 + 270270.0*rij*xij + 
              207900.0*Power(rij,2.0)*Power(xij,2.0) + 
              90090.0*Power(rij,3.0)*Power(xij,3.0) + 
              25200.0*Power(rij,4.0)*Power(xij,4.0) + 
              4788.0*Power(rij,5.0)*Power(xij,5.0) + 616.0*Power(rij,6.0)*Power(xij,6.0) + 
              50.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           7.0*Power(xii,6.0)*Power(xij,10.0)*
            (-24885.0 - 49770.0*rij*xij + 18630.0*Power(rij,2.0)*Power(xij,2.0) - 
              44880.0*Power(rij,3.0)*Power(xij,3.0) - 
              34230.0*Power(rij,4.0)*Power(xij,4.0) - 
              8892.0*Power(rij,5.0)*Power(xij,5.0) - 964.0*Power(rij,6.0)*Power(xij,6.0) - 
              14.0*Power(rij,7.0)*Power(xij,7.0) + 4.0*Power(rij,8.0)*Power(xij,8.0)) + 
           7.0*Power(xii,4.0)*Power(xij,12.0)*
            (28665.0 + 57330.0*rij*xij - 110970.0*Power(rij,2.0)*Power(xij,2.0) - 
              90480.0*Power(rij,3.0)*Power(xij,3.0) - 
              26430.0*Power(rij,4.0)*Power(xij,4.0) - 
              3312.0*Power(rij,5.0)*Power(xij,5.0) - 4.0*Power(rij,6.0)*Power(xij,6.0) + 
              46.0*Power(rij,7.0)*Power(xij,7.0) + 4.0*Power(rij,8.0)*Power(xij,8.0)) - 
           7.0*Power(xii,12.0)*Power(xij,4.0)*
            (2475.0 + 4950.0*rij*xij + 4950.0*Power(rij,2.0)*Power(xij,2.0) + 
              3300.0*Power(rij,3.0)*Power(xij,3.0) + 
              1650.0*Power(rij,4.0)*Power(xij,4.0) + 576.0*Power(rij,5.0)*Power(xij,5.0) + 
              332.0*Power(rij,6.0)*Power(xij,6.0) + 70.0*Power(rij,7.0)*Power(xij,7.0) + 
              4.0*Power(rij,8.0)*Power(xij,8.0)) + 
           7.0*Power(xii,10.0)*Power(xij,6.0)*
            (7425.0 + 14850.0*rij*xij + 14850.0*Power(rij,2.0)*Power(xij,2.0) + 
              10350.0*Power(rij,3.0)*Power(xij,3.0) + 
              3150.0*Power(rij,4.0)*Power(xij,4.0) + 
              3036.0*Power(rij,5.0)*Power(xij,5.0) + 
              1052.0*Power(rij,6.0)*Power(xij,6.0) + 130.0*Power(rij,7.0)*Power(xij,7.0) + 
              4.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xii,14.0)*Power(xij,2.0)*
            (3465.0 + 6930.0*rij*xij + 6930.0*Power(rij,2.0)*Power(xij,2.0) + 
              4620.0*Power(rij,3.0)*Power(xij,3.0) + 2310.0*Power(rij,4.0)*Power(xij,4.0) + 
              924.0*Power(rij,5.0)*Power(xij,5.0) + 308.0*Power(rij,6.0)*Power(xij,6.0) + 
              118.0*Power(rij,7.0)*Power(xij,7.0) + 12.0*Power(rij,8.0)*Power(xij,8.0))) + 
        105.0*Power(E,2.0*rij*xij)*Power(xij,10.0)*
         (-34.0*Power(rij,3.0)*Power(xii,15.0) - 2.0*Power(rij,4.0)*Power(xii,16.0) + 
           3.0*Power(xij,12.0) + 6.0*rij*xii*Power(xij,12.0) - 
           8.0*Power(rij,2.0)*Power(xii,14.0)*(27.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           2.0*rij*Power(xii,3.0)*Power(xij,10.0)*
            (-33.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*Power(xii,2.0)*Power(xij,10.0)*(-11.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           66.0*rij*Power(xii,9.0)*Power(xij,4.0)*
            (-191.0 + 6.0*Power(rij,2.0)*Power(xij,2.0)) - 
           22.0*rij*Power(xii,7.0)*Power(xij,6.0)*
            (162.0 + 11.0*Power(rij,2.0)*Power(xij,2.0)) - 
           21.0*Power(xii,10.0)*Power(xij,2.0)*(157.0 + 66.0*Power(rij,2.0)*Power(xij,2.0)) + 
           2.0*rij*Power(xii,11.0)*Power(xij,2.0)*
            (-3297.0 + 88.0*Power(rij,2.0)*Power(xij,2.0)) - 
           2.0*rij*Power(xii,13.0)*(273.0 + 113.0*Power(rij,2.0)*Power(xij,2.0)) - 
           11.0*Power(xii,8.0)*Power(xij,4.0)*
            (573.0 - 252.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + Power(xii,4.0)*Power(xij,8.0)*(165.0 - 66.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*Power(xii,6.0)*Power(xij,6.0)*
            (-891.0 + 462.0*Power(rij,2.0)*Power(xij,2.0) + 4.0*Power(rij,4.0)*Power(xij,4.0)) \
    + Power(xii,12.0)*(-273.0 - 2034.0*Power(rij,2.0)*Power(xij,2.0) + 
              22.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,5.0)*(330.0*rij*Power(xij,8.0) - 74.0*Power(rij,3.0)*Power(xij,10.0))))/
      (315.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),11.0))
    ;
  }
  return S;
}

static double DSlater_2S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-62270208000.0 + 62270208000.0*Power(E,2.0*rij*xii) - 124540416000.0*rij*xii - 
        124540416000.0*Power(rij,2.0)*Power(xii,2.0) - 
        82915457625.0*Power(rij,3.0)*Power(xii,3.0) - 
        41290499250.0*Power(rij,4.0)*Power(xii,4.0) - 
        16374307950.0*Power(rij,5.0)*Power(xii,5.0) - 
        5370264900.0*Power(rij,6.0)*Power(xii,6.0) - 
        1491272640.0*Power(rij,7.0)*Power(xii,7.0) - 
        355520880.0*Power(rij,8.0)*Power(xii,8.0) - 
        73238880.0*Power(rij,9.0)*Power(xii,9.0) - 
        13041600.0*Power(rij,10.0)*Power(xii,10.0) - 
        1996800.0*Power(rij,11.0)*Power(xii,11.0) - 
        259584.0*Power(rij,12.0)*Power(xii,12.0) - 27648.0*Power(rij,13.0)*Power(xii,13.0) - 
        2048.0*Power(rij,14.0)*Power(xii,14.0))/
      (6.2270208e10*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (14175.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),13.0) + 
        945.0*Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-200.0*Power(rij,3.0)*Power(xii,17.0) - 10.0*Power(rij,4.0)*Power(xii,18.0) + 
           15.0*Power(xij,14.0) + 30.0*rij*xii*Power(xij,14.0) + 
           10.0*rij*Power(xii,3.0)*Power(xij,12.0)*
            (-39.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           15.0*Power(xii,2.0)*Power(xij,12.0)*(-13.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           30.0*Power(rij,2.0)*Power(xii,16.0)*(49.0 + 3.0*Power(rij,2.0)*Power(xij,2.0)) + 
           286.0*rij*Power(xii,9.0)*Power(xij,6.0)*
            (-915.0 + 4.0*Power(rij,2.0)*Power(xij,2.0)) - 
           60.0*rij*Power(xii,5.0)*Power(xij,10.0)*
            (-39.0 + 8.0*Power(rij,2.0)*Power(xij,2.0)) + 
           338.0*rij*Power(xii,11.0)*Power(xij,4.0)*
            (-855.0 + 22.0*Power(rij,2.0)*Power(xij,2.0)) - 
           156.0*rij*Power(xii,7.0)*Power(xij,8.0)*
            (275.0 + 23.0*Power(rij,2.0)*Power(xij,2.0)) - 
           60.0*rij*Power(xii,15.0)*(70.0 + 41.0*Power(rij,2.0)*Power(xij,2.0)) - 
           24.0*rij*Power(xii,13.0)*Power(xij,2.0)*
            (3335.0 + 78.0*Power(rij,2.0)*Power(xij,2.0)) + 
           10.0*Power(xii,4.0)*Power(xij,10.0)*
            (117.0 - 39.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) - 
           39.0*Power(xii,8.0)*Power(xij,6.0)*
            (3355.0 - 1298.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           6.0*Power(xii,14.0)*(-350.0 - 3855.0*Power(rij,2.0)*Power(xij,2.0) + 
              13.0*Power(rij,4.0)*Power(xij,4.0)) + 
           6.0*Power(xii,6.0)*Power(xij,8.0)*
            (-3575.0 + 1391.0*Power(rij,2.0)*Power(xij,2.0) + 
              15.0*Power(rij,4.0)*Power(xij,4.0)) - 
           13.0*Power(xii,10.0)*Power(xij,4.0)*
            (11115.0 - 1386.0*Power(rij,2.0)*Power(xij,2.0) + 
              22.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,12.0)*(-40020.0*Power(xij,2.0) - 
              52026.0*Power(rij,2.0)*Power(xij,4.0) + 286.0*Power(rij,4.0)*Power(xij,6.0))) \
    + Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (819.0*Power(xii,10.0)*Power(xij,10.0)*
            (22275.0 + 44550.0*rij*xij + 41400.0*Power(rij,2.0)*Power(xij,2.0) + 
              43200.0*Power(rij,3.0)*Power(xij,3.0) + 
              3600.0*Power(rij,4.0)*Power(xij,4.0) + 4080.0*Power(rij,5.0)*Power(xij,5.0) + 
              3720.0*Power(rij,6.0)*Power(xij,6.0) + 1104.0*Power(rij,7.0)*Power(xij,7.0) + 
              148.0*Power(rij,8.0)*Power(xij,8.0) + 8.0*Power(rij,9.0)*Power(xij,9.0)) + 
           2.0*Power(xij,20.0)*(8108100.0 + 16216200.0*rij*xij + 
              13378365.0*Power(rij,2.0)*Power(xij,2.0) + 
              6486480.0*Power(rij,3.0)*Power(xij,3.0) + 
              2120580.0*Power(rij,4.0)*Power(xij,4.0) + 
              498960.0*Power(rij,5.0)*Power(xij,5.0) + 
              86940.0*Power(rij,6.0)*Power(xij,6.0) + 
              11232.0*Power(rij,7.0)*Power(xij,7.0) + 
              1044.0*Power(rij,8.0)*Power(xij,8.0) + 64.0*Power(rij,9.0)*Power(xij,9.0) + 
              2.0*Power(rij,10.0)*Power(xij,10.0)) + 
           42.0*Power(xii,8.0)*Power(xij,12.0)*
            (-531900.0 - 1063800.0*rij*xij - 2344275.0*Power(rij,2.0)*Power(xij,2.0) + 
              269100.0*Power(rij,3.0)*Power(xij,3.0) - 
              84150.0*Power(rij,4.0)*Power(xij,4.0) - 
              304740.0*Power(rij,5.0)*Power(xij,5.0) - 
              124440.0*Power(rij,6.0)*Power(xij,6.0) - 
              22800.0*Power(rij,7.0)*Power(xij,7.0) - 
              1980.0*Power(rij,8.0)*Power(xij,8.0) - 40.0*Power(rij,9.0)*Power(xij,9.0) + 
              4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           Power(xii,20.0)*(14175.0 + 28350.0*rij*xij + 
              28350.0*Power(rij,2.0)*Power(xij,2.0) + 
              18900.0*Power(rij,3.0)*Power(xij,3.0) + 
              9450.0*Power(rij,4.0)*Power(xij,4.0) + 3780.0*Power(rij,5.0)*Power(xij,5.0) + 
              1260.0*Power(rij,6.0)*Power(xij,6.0) + 360.0*Power(rij,7.0)*Power(xij,7.0) + 
              90.0*Power(rij,8.0)*Power(xij,8.0) + 20.0*Power(rij,9.0)*Power(xij,9.0) + 
              4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           18.0*Power(xii,16.0)*Power(xij,4.0)*
            (61425.0 + 122850.0*rij*xij + 122850.0*Power(rij,2.0)*Power(xij,2.0) + 
              81900.0*Power(rij,3.0)*Power(xij,3.0) + 
              40950.0*Power(rij,4.0)*Power(xij,4.0) + 
              16380.0*Power(rij,5.0)*Power(xij,5.0) + 
              5460.0*Power(rij,6.0)*Power(xij,6.0) + 1392.0*Power(rij,7.0)*Power(xij,7.0) + 
              544.0*Power(rij,8.0)*Power(xij,8.0) + 104.0*Power(rij,9.0)*Power(xij,9.0) + 
              6.0*Power(rij,10.0)*Power(xij,10.0)) + 
           18.0*Power(xii,4.0)*Power(xij,16.0)*
            (6572475.0 + 13144950.0*rij*xij - 1539720.0*Power(rij,2.0)*Power(xij,2.0) - 
              5741190.0*Power(rij,3.0)*Power(xij,3.0) - 
              2690415.0*Power(rij,4.0)*Power(xij,4.0) - 
              619710.0*Power(rij,5.0)*Power(xij,5.0) - 
              73710.0*Power(rij,6.0)*Power(xij,6.0) - 
              1716.0*Power(rij,7.0)*Power(xij,7.0) + 803.0*Power(rij,8.0)*Power(xij,8.0) + 
              118.0*Power(rij,9.0)*Power(xij,9.0) + 6.0*Power(rij,10.0)*Power(xij,10.0)) - 
           21.0*Power(xii,12.0)*Power(xij,8.0)*
            (482625.0 + 965250.0*rij*xij + 965250.0*Power(rij,2.0)*Power(xij,2.0) + 
              633600.0*Power(rij,3.0)*Power(xij,3.0) + 
              376200.0*Power(rij,4.0)*Power(xij,4.0) + 
              67680.0*Power(rij,5.0)*Power(xij,5.0) + 
              44760.0*Power(rij,6.0)*Power(xij,6.0) + 
              22128.0*Power(rij,7.0)*Power(xij,7.0) + 
              4476.0*Power(rij,8.0)*Power(xij,8.0) + 376.0*Power(rij,9.0)*Power(xij,9.0) + 
              8.0*Power(rij,10.0)*Power(xij,10.0)) - 
           12.0*Power(xii,6.0)*Power(xij,14.0)*
            (-5178600.0 - 10357200.0*rij*xij + 8303715.0*Power(rij,2.0)*Power(xij,2.0) + 
              687330.0*Power(rij,3.0)*Power(xij,3.0) - 
              3292695.0*Power(rij,4.0)*Power(xij,4.0) - 
              1634850.0*Power(rij,5.0)*Power(xij,5.0) - 
              362040.0*Power(rij,6.0)*Power(xij,6.0) - 
              40728.0*Power(rij,7.0)*Power(xij,7.0) - 
              1446.0*Power(rij,8.0)*Power(xij,8.0) + 164.0*Power(rij,9.0)*Power(xij,9.0) + 
              16.0*Power(rij,10.0)*Power(xij,10.0)) - 
           2.0*Power(xii,2.0)*Power(xij,18.0)*
            (-66891825.0 - 133783650.0*rij*xij - 
              89594505.0*Power(rij,2.0)*Power(xij,2.0) - 
              30540510.0*Power(rij,3.0)*Power(xij,3.0) - 
              5540535.0*Power(rij,4.0)*Power(xij,4.0) - 
              270270.0*Power(rij,5.0)*Power(xij,5.0) + 
              125370.0*Power(rij,6.0)*Power(xij,6.0) + 
              37116.0*Power(rij,7.0)*Power(xij,7.0) + 
              5247.0*Power(rij,8.0)*Power(xij,8.0) + 422.0*Power(rij,9.0)*Power(xij,9.0) + 
              16.0*Power(rij,10.0)*Power(xij,10.0)) + 
           Power(xii,18.0)*Power(xij,2.0)*
            (184275.0 + 368550.0*rij*xij + 368550.0*Power(rij,2.0)*Power(xij,2.0) + 
              245700.0*Power(rij,3.0)*Power(xij,3.0) + 
              122850.0*Power(rij,4.0)*Power(xij,4.0) + 
              49140.0*Power(rij,5.0)*Power(xij,5.0) + 
              16380.0*Power(rij,6.0)*Power(xij,6.0) + 
              4680.0*Power(rij,7.0)*Power(xij,7.0) + 1170.0*Power(rij,8.0)*Power(xij,8.0) + 
              340.0*Power(rij,9.0)*Power(xij,9.0) + 32.0*Power(rij,10.0)*Power(xij,10.0)) + 
           6.0*Power(xii,14.0)*Power(xij,6.0)*
            (675675.0 + 1351350.0*rij*xij + 1351350.0*Power(rij,2.0)*Power(xij,2.0) + 
              900900.0*Power(rij,3.0)*Power(xij,3.0) + 
              450450.0*Power(rij,4.0)*Power(xij,4.0) + 
              187740.0*Power(rij,5.0)*Power(xij,5.0) + 
              43680.0*Power(rij,6.0)*Power(xij,6.0) + 
              22128.0*Power(rij,7.0)*Power(xij,7.0) + 6876.0*Power(rij,8.0)*Power(xij,8.0) + 
              856.0*Power(rij,9.0)*Power(xij,9.0) + 32.0*Power(rij,10.0)*Power(xij,10.0))))/
      (14175.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),13.0))
    ;
  }
  return S;
}

static double DSlater_2S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-31384184832000.0 + 31384184832000.0*Power(E,2.0*rij*xii) - 
        62768369664000.0*rij*xii - 62768369664000.0*Power(rij,2.0)*Power(xii,2.0) - 
        41826211552125.0*Power(rij,3.0)*Power(xii,3.0) - 
        20884053440250.0*Power(rij,4.0)*Power(xii,4.0) - 
        8328251131200.0*Power(rij,5.0)*Power(xii,5.0) - 
        2759624267400.0*Power(rij,6.0)*Power(xii,6.0) - 
        779901922800.0*Power(rij,7.0)*Power(xii,7.0) - 
        191286295200.0*Power(rij,8.0)*Power(xii,8.0) - 
        41167526400.0*Power(rij,9.0)*Power(xii,9.0) - 
        7818370560.0*Power(rij,10.0)*Power(xii,10.0) - 
        1311448320.0*Power(rij,11.0)*Power(xii,11.0) - 
        193589760.0*Power(rij,12.0)*Power(xii,12.0) - 
        24944640.0*Power(rij,13.0)*Power(xii,13.0) - 
        2764800.0*Power(rij,14.0)*Power(xii,14.0) - 
        253952.0*Power(rij,15.0)*Power(xii,15.0) - 16384.0*Power(rij,16.0)*Power(xii,16.0))/
      (3.1384184832e13*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (467775.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),15.0) + 
        51975.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-138.0*Power(rij,3.0)*Power(xii,19.0) - 6.0*Power(rij,4.0)*Power(xii,20.0) + 
           9.0*Power(xij,16.0) + 18.0*rij*xii*Power(xij,16.0) + 
           2.0*rij*Power(xii,5.0)*Power(xij,12.0)*
            (945.0 - 181.0*Power(rij,2.0)*Power(xij,2.0)) + 
           6.0*rij*Power(xii,3.0)*Power(xij,14.0)*
            (-45.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           9.0*Power(xii,2.0)*Power(xij,14.0)*(-15.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           4.0*Power(rij,2.0)*Power(xii,18.0)*(288.0 + 23.0*Power(rij,2.0)*Power(xij,2.0)) + 
           234.0*rij*Power(xii,11.0)*Power(xij,6.0)*
            (-4209.0 + 55.0*Power(rij,2.0)*Power(xij,2.0)) - 
           78.0*rij*Power(xii,9.0)*Power(xij,8.0)*
            (6655.0 + 63.0*Power(rij,2.0)*Power(xij,2.0)) + 
           18.0*rij*Power(xii,13.0)*Power(xij,4.0)*
            (-31885.0 + 377.0*Power(rij,2.0)*Power(xij,2.0)) - 
           6.0*rij*Power(xii,7.0)*Power(xij,10.0)*
            (9321.0 + 791.0*Power(rij,2.0)*Power(xij,2.0)) - 
           6.0*rij*Power(xii,15.0)*Power(xij,2.0)*
            (16755.0 + 1141.0*Power(rij,2.0)*Power(xij,2.0)) - 
           2.0*rij*Power(xii,17.0)*(1836.0 + 1331.0*Power(rij,2.0)*Power(xij,2.0)) - 
           9.0*Power(xii,12.0)*Power(xij,4.0)*
            (31885.0 + 7514.0*Power(rij,2.0)*Power(xij,2.0)) + 
           3.0*Power(xii,4.0)*Power(xij,12.0)*
            (315.0 - 90.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    - 117.0*Power(xii,10.0)*Power(xij,6.0)*
            (4209.0 - 990.0*Power(rij,2.0)*Power(xij,2.0) + 
              4.0*Power(rij,4.0)*Power(xij,4.0)) - 
           6.0*Power(xii,16.0)*(306.0 + 4491.0*Power(rij,2.0)*Power(xij,2.0) + 
              14.0*Power(rij,4.0)*Power(xij,4.0)) + 
           3.0*Power(xii,8.0)*Power(xij,8.0)*
            (-86515.0 + 28158.0*Power(rij,2.0)*Power(xij,2.0) + 
              28.0*Power(rij,4.0)*Power(xij,4.0)) + 
           9.0*Power(xii,14.0)*Power(xij,2.0)*
            (-5585.0 - 12530.0*Power(rij,2.0)*Power(xij,2.0) + 
              52.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,6.0)*Power(xij,10.0)*
            (-27963.0 + 8442.0*Power(rij,2.0)*Power(xij,2.0) + 
              92.0*Power(rij,4.0)*Power(xij,4.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,6.0)*
         (-3465.0*Power(xii,12.0)*Power(xij,12.0)*
            (675675.0 + 1351350.0*rij*xij + 1389150.0*Power(rij,2.0)*Power(xij,2.0) + 
              689850.0*Power(rij,3.0)*Power(xij,3.0) + 
              730800.0*Power(rij,4.0)*Power(xij,4.0) + 
              128520.0*Power(rij,5.0)*Power(xij,5.0) + 
              9240.0*Power(rij,6.0)*Power(xij,6.0) + 
              18480.0*Power(rij,7.0)*Power(xij,7.0) + 
              8820.0*Power(rij,8.0)*Power(xij,8.0) + 1800.0*Power(rij,9.0)*Power(xij,9.0) + 
              184.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           330.0*Power(xii,8.0)*Power(xij,16.0)*
            (-1204875.0 - 2409750.0*rij*xij - 75042450.0*Power(rij,2.0)*Power(xij,2.0) - 
              3403575.0*Power(rij,3.0)*Power(xij,3.0) + 
              9111375.0*Power(rij,4.0)*Power(xij,4.0) - 
              498330.0*Power(rij,5.0)*Power(xij,5.0) - 
              1892310.0*Power(rij,6.0)*Power(xij,6.0) - 
              669180.0*Power(rij,7.0)*Power(xij,7.0) - 
              118785.0*Power(rij,8.0)*Power(xij,8.0) - 
              11650.0*Power(rij,9.0)*Power(xij,9.0) - 
              482.0*Power(rij,10.0)*Power(xij,10.0) + 
              16.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) + 
           2.0*Power(xij,24.0)*(620269650.0 + 1240539300.0*rij*xij + 
              1070269200.0*Power(rij,2.0)*Power(xij,2.0) + 
              557431875.0*Power(rij,3.0)*Power(xij,3.0) + 
              200675475.0*Power(rij,4.0)*Power(xij,4.0) + 
              53513460.0*Power(rij,5.0)*Power(xij,5.0) + 
              10977120.0*Power(rij,6.0)*Power(xij,6.0) + 
              1764180.0*Power(rij,7.0)*Power(xij,7.0) + 
              222750.0*Power(rij,8.0)*Power(xij,8.0) + 
              21780.0*Power(rij,9.0)*Power(xij,9.0) + 
              1584.0*Power(rij,10.0)*Power(xij,10.0) + 
              78.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) - 
           Power(xii,24.0)*(467775.0 + 935550.0*rij*xij + 
              935550.0*Power(rij,2.0)*Power(xij,2.0) + 
              623700.0*Power(rij,3.0)*Power(xij,3.0) + 
              311850.0*Power(rij,4.0)*Power(xij,4.0) + 
              124740.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              11880.0*Power(rij,7.0)*Power(xij,7.0) + 
              2970.0*Power(rij,8.0)*Power(xij,8.0) + 660.0*Power(rij,9.0)*Power(xij,9.0) + 
              132.0*Power(rij,10.0)*Power(xij,10.0) + 
              24.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           110.0*Power(xii,6.0)*Power(xij,18.0)*
            (-156874725.0 - 313749450.0*rij*xij + 
              119665350.0*Power(rij,2.0)*Power(xij,2.0) + 
              104285475.0*Power(rij,3.0)*Power(xij,3.0) - 
              1941975.0*Power(rij,4.0)*Power(xij,4.0) - 
              17730090.0*Power(rij,5.0)*Power(xij,5.0) - 
              6892830.0*Power(rij,6.0)*Power(xij,6.0) - 
              1379700.0*Power(rij,7.0)*Power(xij,7.0) - 
              159705.0*Power(rij,8.0)*Power(xij,8.0) - 
              8610.0*Power(rij,9.0)*Power(xij,9.0) + 
              294.0*Power(rij,10.0)*Power(xij,10.0) + 
              78.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           10.0*Power(xii,2.0)*Power(xij,22.0)*
            (-1412836425.0 - 2825672850.0*rij*xij - 
              2095943850.0*Power(rij,2.0)*Power(xij,2.0) - 
              854728875.0*Power(rij,3.0)*Power(xij,3.0) - 
              215540325.0*Power(rij,4.0)*Power(xij,4.0) - 
              32702670.0*Power(rij,5.0)*Power(xij,5.0) - 
              1753290.0*Power(rij,6.0)*Power(xij,6.0) + 
              479160.0*Power(rij,7.0)*Power(xij,7.0) + 
              150975.0*Power(rij,8.0)*Power(xij,8.0) + 
              22990.0*Power(rij,9.0)*Power(xij,9.0) + 
              2222.0*Power(rij,10.0)*Power(xij,10.0) + 
              134.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           165.0*Power(xii,16.0)*Power(xij,8.0)*
            (3869775.0 + 7739550.0*rij*xij + 7739550.0*Power(rij,2.0)*Power(xij,2.0) + 
              5159700.0*Power(rij,3.0)*Power(xij,3.0) + 
              2579850.0*Power(rij,4.0)*Power(xij,4.0) + 
              1018080.0*Power(rij,5.0)*Power(xij,5.0) + 
              385560.0*Power(rij,6.0)*Power(xij,6.0) + 
              70920.0*Power(rij,7.0)*Power(xij,7.0) + 
              21720.0*Power(rij,8.0)*Power(xij,8.0) + 
              8640.0*Power(rij,9.0)*Power(xij,9.0) + 
              1704.0*Power(rij,10.0)*Power(xij,10.0) + 
              148.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) + 
           5.0*Power(xii,22.0)*Power(xij,2.0)*
            (1403325.0 + 2806650.0*rij*xij + 2806650.0*Power(rij,2.0)*Power(xij,2.0) + 
              1871100.0*Power(rij,3.0)*Power(xij,3.0) + 
              935550.0*Power(rij,4.0)*Power(xij,4.0) + 
              374220.0*Power(rij,5.0)*Power(xij,5.0) + 
              124740.0*Power(rij,6.0)*Power(xij,6.0) + 
              35640.0*Power(rij,7.0)*Power(xij,7.0) + 
              8910.0*Power(rij,8.0)*Power(xij,8.0) + 1980.0*Power(rij,9.0)*Power(xij,9.0) + 
              396.0*Power(rij,10.0)*Power(xij,10.0) + 
              92.0*Power(rij,11.0)*Power(xij,11.0) + 8.0*Power(rij,12.0)*Power(xij,12.0)) + 
           55.0*Power(xii,18.0)*Power(xij,6.0)*
            (3869775.0 + 7739550.0*rij*xij + 7739550.0*Power(rij,2.0)*Power(xij,2.0) + 
              5159700.0*Power(rij,3.0)*Power(xij,3.0) + 
              2579850.0*Power(rij,4.0)*Power(xij,4.0) + 
              1031940.0*Power(rij,5.0)*Power(xij,5.0) + 
              343980.0*Power(rij,6.0)*Power(xij,6.0) + 
              101520.0*Power(rij,7.0)*Power(xij,7.0) + 
              19710.0*Power(rij,8.0)*Power(xij,8.0) + 
              6300.0*Power(rij,9.0)*Power(xij,9.0) + 
              1692.0*Power(rij,10.0)*Power(xij,10.0) + 
              204.0*Power(rij,11.0)*Power(xij,11.0) + 8.0*Power(rij,12.0)*Power(xij,12.0)) + 
           22.0*Power(xii,4.0)*Power(xij,20.0)*
            (1099568925.0 + 2199137850.0*rij*xij + 
              699139350.0*Power(rij,2.0)*Power(xij,2.0) - 
              242690175.0*Power(rij,3.0)*Power(xij,3.0) - 
              237899025.0*Power(rij,4.0)*Power(xij,4.0) - 
              81646110.0*Power(rij,5.0)*Power(xij,5.0) - 
              16127370.0*Power(rij,6.0)*Power(xij,6.0) - 
              1875420.0*Power(rij,7.0)*Power(xij,7.0) - 
              82035.0*Power(rij,8.0)*Power(xij,8.0) + 
              11970.0*Power(rij,9.0)*Power(xij,9.0) + 
              2574.0*Power(rij,10.0)*Power(xij,10.0) + 
              218.0*Power(rij,11.0)*Power(xij,11.0) + 8.0*Power(rij,12.0)*Power(xij,12.0)) - 
           33.0*Power(xii,10.0)*Power(xij,14.0)*
            (-94107825.0 - 188215650.0*rij*xij - 
              72320850.0*Power(rij,2.0)*Power(xij,2.0) - 
              284964750.0*Power(rij,3.0)*Power(xij,3.0) - 
              43356600.0*Power(rij,4.0)*Power(xij,4.0) + 
              10299240.0*Power(rij,5.0)*Power(xij,5.0) - 
              5319720.0*Power(rij,6.0)*Power(xij,6.0) - 
              4942080.0*Power(rij,7.0)*Power(xij,7.0) - 
              1373700.0*Power(rij,8.0)*Power(xij,8.0) - 
              194760.0*Power(rij,9.0)*Power(xij,9.0) - 
              14088.0*Power(rij,10.0)*Power(xij,10.0) - 
              296.0*Power(rij,11.0)*Power(xij,11.0) + 16.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 11.0*Power(xii,20.0)*Power(xij,4.0)*
            (4465125.0 + 8930250.0*rij*xij + 8930250.0*Power(rij,2.0)*Power(xij,2.0) + 
              5953500.0*Power(rij,3.0)*Power(xij,3.0) + 
              2976750.0*Power(rij,4.0)*Power(xij,4.0) + 
              1190700.0*Power(rij,5.0)*Power(xij,5.0) + 
              396900.0*Power(rij,6.0)*Power(xij,6.0) + 
              113400.0*Power(rij,7.0)*Power(xij,7.0) + 
              28350.0*Power(rij,8.0)*Power(xij,8.0) + 
              5740.0*Power(rij,9.0)*Power(xij,9.0) + 
              1652.0*Power(rij,10.0)*Power(xij,10.0) + 
              284.0*Power(rij,11.0)*Power(xij,11.0) + 16.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 33.0*Power(xii,14.0)*Power(xij,10.0)*
            (42567525.0 + 85135050.0*rij*xij + 85135050.0*Power(rij,2.0)*Power(xij,2.0) + 
              57043350.0*Power(rij,3.0)*Power(xij,3.0) + 
              26371800.0*Power(rij,4.0)*Power(xij,4.0) + 
              14668920.0*Power(rij,5.0)*Power(xij,5.0) + 
              2621640.0*Power(rij,6.0)*Power(xij,6.0) + 
              597840.0*Power(rij,7.0)*Power(xij,7.0) + 
              378780.0*Power(rij,8.0)*Power(xij,8.0) + 
              114040.0*Power(rij,9.0)*Power(xij,9.0) + 
              16088.0*Power(rij,10.0)*Power(xij,10.0) + 
              1016.0*Power(rij,11.0)*Power(xij,11.0) + 16.0*Power(rij,12.0)*Power(xij,12.0)))\
    )/(467775.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),15.0))
    ;
  }
  return S;
}

double DSlater_2S_1S(double rij,double xii,double xij)
{
  return DSlater_1S_2S(rij,xij,xii);
}

static double DSlater_3S_3S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-359251200.0 + 359251200.0*Power(E,2.0*rij*xii) - 718502400.0*rij*xii - 
        718502400.0*Power(rij,2.0)*Power(xii,2.0) - 
        475727175.0*Power(rij,3.0)*Power(xii,3.0) - 
        232951950.0*Power(rij,4.0)*Power(xii,4.0) - 
        89397000.0*Power(rij,5.0)*Power(xii,5.0) - 
        27858600.0*Power(rij,6.0)*Power(xii,6.0) - 7223040.0*Power(rij,7.0)*Power(xii,7.0) - 
        1584000.0*Power(rij,8.0)*Power(xii,8.0) - 295680.0*Power(rij,9.0)*Power(xii,9.0) - 
        46464.0*Power(rij,10.0)*Power(xii,10.0) - 5888.0*Power(rij,11.0)*Power(xii,11.0) - 
        512.0*Power(rij,12.0)*Power(xii,12.0))/
      (3.592512e8*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (135.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),11.0) + 
        Power(E,2.0*rij*xij)*Power(xij,8.0)*
         (-276.0*Power(rij,5.0)*Power(xii,19.0) - 12.0*Power(rij,6.0)*Power(xii,20.0) + 
           135.0*Power(xij,14.0) + 270.0*rij*xii*Power(xij,14.0) - 
           100.0*Power(rij,3.0)*Power(xii,17.0)*(165.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           10.0*Power(rij,4.0)*Power(xii,18.0)*(-285.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           90.0*rij*Power(xii,3.0)*Power(xij,12.0)*
            (-33.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           135.0*Power(xii,2.0)*Power(xij,12.0)*(-11.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           18.0*rij*Power(xii,5.0)*Power(xij,10.0)*
            (825.0 - 110.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 45.0*Power(xii,4.0)*Power(xij,10.0)*
            (165.0 - 66.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    - 10.0*rij*Power(xii,7.0)*Power(xij,8.0)*
            (4455.0 - 738.0*Power(rij,2.0)*Power(xij,2.0) + 
              62.0*Power(rij,4.0)*Power(xij,4.0)) + 
           10.0*rij*Power(xii,11.0)*Power(xij,4.0)*
            (-96831.0 + 6534.0*Power(rij,2.0)*Power(xij,2.0) + 
              154.0*Power(rij,4.0)*Power(xij,4.0)) - 
           10.0*rij*Power(xii,13.0)*Power(xij,2.0)*
            (84357.0 - 12318.0*Power(rij,2.0)*Power(xij,2.0) + 
              418.0*Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*rij*Power(xii,9.0)*Power(xij,6.0)*
            (-495.0 - 48510.0*Power(rij,2.0)*Power(xij,2.0) + 
              458.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,15.0)*(-90090.0*rij - 80580.0*Power(rij,3.0)*Power(xij,2.0) + 
              2684.0*Power(rij,5.0)*Power(xij,4.0)) + 
           Power(xii,16.0)*(-54450.0*Power(rij,2.0) - 
              7290.0*Power(rij,4.0)*Power(xij,2.0) + 68.0*Power(rij,6.0)*Power(xij,4.0)) - 
           5.0*Power(xii,8.0)*Power(xij,6.0)*
            (99.0 + 1782.0*Power(rij,2.0)*Power(xij,2.0) - 
              2250.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           3.0*Power(xii,6.0)*Power(xij,8.0)*
            (-7425.0 + 4950.0*Power(rij,2.0)*Power(xij,2.0) - 
              330.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           5.0*Power(xii,14.0)*(9009.0 + 78954.0*Power(rij,2.0)*Power(xij,2.0) - 
              6030.0*Power(rij,4.0)*Power(xij,4.0) + 44.0*Power(rij,6.0)*Power(xij,6.0)) + 
           5.0*Power(xii,12.0)*Power(xij,2.0)*
            (-84357.0 - 366.0*Power(rij,2.0)*Power(xij,2.0) - 
              3498.0*Power(rij,4.0)*Power(xij,4.0) + 44.0*Power(rij,6.0)*Power(xij,6.0)) - 
           Power(xii,10.0)*Power(xij,4.0)*
            (484155.0 - 447810.0*Power(rij,2.0)*Power(xij,2.0) + 
              12870.0*Power(rij,4.0)*Power(xij,4.0) + 68.0*Power(rij,6.0)*Power(xij,6.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,8.0)*
         (Power(xii,4.0)*Power(xij,10.0)*
            (484155.0 + 968310.0*rij*xij + 1830.0*Power(rij,2.0)*Power(xij,2.0) - 
              123180.0*Power(rij,3.0)*Power(xij,3.0) - 
              30150.0*Power(rij,4.0)*Power(xij,4.0) - 
              2684.0*Power(rij,5.0)*Power(xij,5.0) - 68.0*Power(rij,6.0)*Power(xij,6.0)) + 
           5.0*Power(xii,2.0)*Power(xij,12.0)*
            (84357.0 + 168714.0*rij*xij + 78954.0*Power(rij,2.0)*Power(xij,2.0) + 
              16116.0*Power(rij,3.0)*Power(xij,3.0) + 
              1458.0*Power(rij,4.0)*Power(xij,4.0) + 20.0*Power(rij,5.0)*Power(xij,5.0) - 
              4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           3.0*Power(xii,14.0)*(45.0 + 90.0*rij*xij + 90.0*Power(rij,2.0)*Power(xij,2.0) + 
              60.0*Power(rij,3.0)*Power(xij,3.0) + 30.0*Power(rij,4.0)*Power(xij,4.0) + 
              12.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           55.0*Power(xii,8.0)*Power(xij,6.0)*
            (-405.0 - 810.0*rij*xij - 162.0*Power(rij,2.0)*Power(xij,2.0) - 
              1764.0*Power(rij,3.0)*Power(xij,3.0) - 234.0*Power(rij,4.0)*Power(xij,4.0) + 
              28.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           55.0*Power(xii,6.0)*Power(xij,8.0)*
            (9.0 + 18.0*rij*xij - 8142.0*Power(rij,2.0)*Power(xij,2.0) - 
              1188.0*Power(rij,3.0)*Power(xij,3.0) + 318.0*Power(rij,4.0)*Power(xij,4.0) + 
              76.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           3.0*Power(xij,14.0)*(15015.0 + 30030.0*rij*xij + 
              18150.0*Power(rij,2.0)*Power(xij,2.0) + 
              5500.0*Power(rij,3.0)*Power(xij,3.0) + 950.0*Power(rij,4.0)*Power(xij,4.0) + 
              92.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           5.0*Power(xii,12.0)*Power(xij,2.0)*
            (297.0 + 594.0*rij*xij + 594.0*Power(rij,2.0)*Power(xij,2.0) + 
              396.0*Power(rij,3.0)*Power(xij,3.0) + 198.0*Power(rij,4.0)*Power(xij,4.0) + 
              124.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           Power(xii,10.0)*Power(xij,4.0)*
            (-7425.0 - 14850.0*rij*xij - 14850.0*Power(rij,2.0)*Power(xij,2.0) - 
              7380.0*Power(rij,3.0)*Power(xij,3.0) - 11250.0*Power(rij,4.0)*Power(xij,4.0) - 
              916.0*Power(rij,5.0)*Power(xij,5.0) + 68.0*Power(rij,6.0)*Power(xij,6.0))))/
      (135.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),11.0))
    ;
  }
  return S;
}

static double DSlater_3S_4S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-37362124800.0 + 37362124800.0*Power(E,2.0*rij*xii) - 74724249600.0*rij*xii - 
        74724249600.0*Power(rij,2.0)*Power(xii,2.0) - 
        49615490925.0*Power(rij,3.0)*Power(xii,3.0) - 
        24506732250.0*Power(rij,4.0)*Power(xii,4.0) - 
        9566747190.0*Power(rij,5.0)*Power(xii,5.0) - 
        3063240180.0*Power(rij,6.0)*Power(xii,6.0) - 
        824709600.0*Power(rij,7.0)*Power(xii,7.0) - 
        189961200.0*Power(rij,8.0)*Power(xii,8.0) - 
        37889280.0*Power(rij,9.0)*Power(xii,9.0) - 
        6589440.0*Power(rij,10.0)*Power(xii,10.0) - 
        998400.0*Power(rij,11.0)*Power(xii,11.0) - 
        129792.0*Power(rij,12.0)*Power(xii,12.0) - 13824.0*Power(rij,13.0)*Power(xii,13.0) - 
        1024.0*Power(rij,14.0)*Power(xii,14.0))/
      (3.73621248e10*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (945.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),13.0) + 
        21.0*Power(E,2.0*rij*xij)*Power(xij,10.0)*
         (-112.0*Power(rij,5.0)*Power(xii,21.0) - 4.0*Power(rij,6.0)*Power(xii,22.0) + 
           45.0*Power(xij,16.0) + 90.0*rij*xii*Power(xij,16.0) + 
           30.0*rij*Power(xii,3.0)*Power(xij,14.0)*
            (-39.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           45.0*Power(xii,2.0)*Power(xij,14.0)*(-13.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           4.0*Power(rij,4.0)*Power(xii,20.0)*(345.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           4.0*Power(rij,3.0)*Power(xii,19.0)*
            (2340.0 + 131.0*Power(rij,2.0)*Power(xij,2.0)) + 
           12.0*rij*Power(xii,5.0)*Power(xij,12.0)*
            (585.0 - 65.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           30.0*Power(xii,4.0)*Power(xij,12.0)*
            (117.0 - 39.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           6.0*Power(rij,2.0)*Power(xii,18.0)*
            (-5915.0 - 1735.0*Power(rij,2.0)*Power(xij,2.0) + 
              12.0*Power(rij,4.0)*Power(xij,4.0)) - 
           12.0*rij*Power(xii,7.0)*Power(xij,10.0)*
            (2145.0 - 225.0*Power(rij,2.0)*Power(xij,2.0) + 
              23.0*Power(rij,4.0)*Power(xij,4.0)) + 
           78.0*rij*Power(xii,11.0)*Power(xij,6.0)*
            (-22875.0 - 770.0*Power(rij,2.0)*Power(xij,2.0) + 
              34.0*Power(rij,4.0)*Power(xij,4.0)) - 
           6.0*rij*Power(xii,9.0)*Power(xij,8.0)*
            (10725.0 + 15730.0*Power(rij,2.0)*Power(xij,2.0) + 
              46.0*Power(rij,4.0)*Power(xij,4.0)) - 
           20.0*rij*Power(xii,13.0)*Power(xij,4.0)*
            (153630.0 - 13923.0*Power(rij,2.0)*Power(xij,2.0) + 
              143.0*Power(rij,4.0)*Power(xij,4.0)) - 
           4.0*rij*Power(xii,15.0)*Power(xij,2.0)*
            (269010.0 + 4455.0*Power(rij,2.0)*Power(xij,2.0) + 
              143.0*Power(rij,4.0)*Power(xij,4.0)) + 
           12.0*rij*Power(xii,17.0)*(-5460.0 - 8235.0*Power(rij,2.0)*Power(xij,2.0) + 
              163.0*Power(rij,4.0)*Power(xij,4.0)) + 
           30.0*Power(xii,14.0)*Power(xij,2.0)*
            (-17934.0 - 26106.0*Power(rij,2.0)*Power(xij,2.0) + 
              793.0*Power(rij,4.0)*Power(xij,4.0)) + 
           3.0*Power(xii,10.0)*Power(xij,6.0)*
            (-297375.0 + 211640.0*Power(rij,2.0)*Power(xij,2.0) + 
              1430.0*Power(rij,4.0)*Power(xij,4.0) - 24.0*Power(rij,6.0)*Power(xij,6.0)) + 
           2.0*Power(xii,6.0)*Power(xij,10.0)*
            (-6435.0 + 3510.0*Power(rij,2.0)*Power(xij,2.0) - 
              195.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           Power(xii,8.0)*Power(xij,8.0)*
            (-32175.0 + 12870.0*Power(rij,2.0)*Power(xij,2.0) + 
              7290.0*Power(rij,4.0)*Power(xij,4.0) + 8.0*Power(rij,6.0)*Power(xij,6.0)) + 
           2.0*Power(xii,12.0)*Power(xij,4.0)*
            (-768150.0 + 324285.0*Power(rij,2.0)*Power(xij,2.0) - 
              19305.0*Power(rij,4.0)*Power(xij,4.0) + 52.0*Power(rij,6.0)*Power(xij,6.0)) - 
           2.0*Power(xii,16.0)*(16380.0 + 241815.0*Power(rij,2.0)*Power(xij,2.0) - 
              7695.0*Power(rij,4.0)*Power(xij,4.0) + 52.0*Power(rij,6.0)*Power(xij,6.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,8.0)*
         (2.0*Power(xii,2.0)*Power(xij,16.0)*
            (8759205.0 + 17518410.0*rij*xij + 10176075.0*Power(rij,2.0)*Power(xij,2.0) + 
              2940210.0*Power(rij,3.0)*Power(xij,3.0) + 
              479115.0*Power(rij,4.0)*Power(xij,4.0) + 
              41496.0*Power(rij,5.0)*Power(xij,5.0) + 882.0*Power(rij,6.0)*Power(xij,6.0) - 
              156.0*Power(rij,7.0)*Power(xij,7.0) - 11.0*Power(rij,8.0)*Power(xij,8.0)) + 
           6.0*Power(xij,18.0)*(225225.0 + 450450.0*rij*xij + 
              315315.0*Power(rij,2.0)*Power(xij,2.0) + 
              120120.0*Power(rij,3.0)*Power(xij,3.0) + 
              28875.0*Power(rij,4.0)*Power(xij,4.0) + 
              4620.0*Power(rij,5.0)*Power(xij,5.0) + 490.0*Power(rij,6.0)*Power(xij,6.0) + 
              32.0*Power(rij,7.0)*Power(xij,7.0) + Power(rij,8.0)*Power(xij,8.0)) - 
           3.0*Power(xii,18.0)*(315.0 + 630.0*rij*xij + 630.0*Power(rij,2.0)*Power(xij,2.0) + 
              420.0*Power(rij,3.0)*Power(xij,3.0) + 210.0*Power(rij,4.0)*Power(xij,4.0) + 
              84.0*Power(rij,5.0)*Power(xij,5.0) + 28.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           84.0*Power(xii,6.0)*Power(xij,12.0)*
            (115200.0 + 230400.0*rij*xij - 341955.0*Power(rij,2.0)*Power(xij,2.0) - 
              141690.0*Power(rij,3.0)*Power(xij,3.0) - 
              6135.0*Power(rij,4.0)*Power(xij,4.0) + 5031.0*Power(rij,5.0)*Power(xij,5.0) + 
              1042.0*Power(rij,6.0)*Power(xij,6.0) + 80.0*Power(rij,7.0)*Power(xij,7.0) + 
              2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           4.0*Power(xii,4.0)*Power(xij,14.0)*
            (-8470980.0 - 16941960.0*rij*xij - 4169655.0*Power(rij,2.0)*Power(xij,2.0) + 
              871605.0*Power(rij,3.0)*Power(xij,3.0) + 
              572985.0*Power(rij,4.0)*Power(xij,4.0) + 
              113169.0*Power(rij,5.0)*Power(xij,5.0) + 
              10878.0*Power(rij,6.0)*Power(xij,6.0) + 456.0*Power(rij,7.0)*Power(xij,7.0) + 
              2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           2.0*Power(xii,14.0)*Power(xij,4.0)*
            (-36855.0 - 73710.0*rij*xij - 73710.0*Power(rij,2.0)*Power(xij,2.0) - 
              49140.0*Power(rij,3.0)*Power(xij,3.0) - 
              24570.0*Power(rij,4.0)*Power(xij,4.0) - 
              7182.0*Power(rij,5.0)*Power(xij,5.0) - 6804.0*Power(rij,6.0)*Power(xij,6.0) - 
              768.0*Power(rij,7.0)*Power(xij,7.0) + 4.0*Power(rij,8.0)*Power(xij,8.0)) - 
           42.0*Power(xii,12.0)*Power(xij,6.0)*
            (-6435.0 - 12870.0*rij*xij - 12870.0*Power(rij,2.0)*Power(xij,2.0) - 
              9570.0*Power(rij,3.0)*Power(xij,3.0) - 330.0*Power(rij,4.0)*Power(xij,4.0) - 
              4434.0*Power(rij,5.0)*Power(xij,5.0) - 908.0*Power(rij,6.0)*Power(xij,6.0) - 
              16.0*Power(rij,7.0)*Power(xij,7.0) + 4.0*Power(rij,8.0)*Power(xij,8.0)) + 
           91.0*Power(xii,10.0)*Power(xij,8.0)*
            (-7425.0 - 14850.0*rij*xij - 21780.0*Power(rij,2.0)*Power(xij,2.0) + 
              11880.0*Power(rij,3.0)*Power(xij,3.0) - 
              15840.0*Power(rij,4.0)*Power(xij,4.0) - 
              5100.0*Power(rij,5.0)*Power(xij,5.0) - 240.0*Power(rij,6.0)*Power(xij,6.0) + 
              48.0*Power(rij,7.0)*Power(xij,7.0) + 4.0*Power(rij,8.0)*Power(xij,8.0)) - 
           91.0*Power(xii,8.0)*Power(xij,10.0)*
            (-20925.0 - 41850.0*rij*xij + 94860.0*Power(rij,2.0)*Power(xij,2.0) - 
              81180.0*Power(rij,3.0)*Power(xij,3.0) - 
              34560.0*Power(rij,4.0)*Power(xij,4.0) - 
              2292.0*Power(rij,5.0)*Power(xij,5.0) + 576.0*Power(rij,6.0)*Power(xij,6.0) + 
              96.0*Power(rij,7.0)*Power(xij,7.0) + 4.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xii,16.0)*Power(xij,2.0)*
            (12285.0 + 24570.0*rij*xij + 24570.0*Power(rij,2.0)*Power(xij,2.0) + 
              16380.0*Power(rij,3.0)*Power(xij,3.0) + 8190.0*Power(rij,4.0)*Power(xij,4.0) + 
              3276.0*Power(rij,5.0)*Power(xij,5.0) + 1092.0*Power(rij,6.0)*Power(xij,6.0) + 
              480.0*Power(rij,7.0)*Power(xij,7.0) + 22.0*Power(rij,8.0)*Power(xij,8.0))))/
      (945.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),13.0))
    ;
  }
  return S;
}

static double DSlater_3S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-78460462080000.0 + 78460462080000.0*Power(E,2.0*rij*xii) - 
        156920924160000.0*rij*xii - 156920924160000.0*Power(rij,2.0)*Power(xii,2.0) - 
        104400898977375.0*Power(rij,3.0)*Power(xii,3.0) - 
        51880873794750.0*Power(rij,4.0)*Power(xii,4.0) - 
        20495752477200.0*Power(rij,5.0)*Power(xii,5.0) - 
        6688323041400.0*Power(rij,6.0)*Power(xii,6.0) - 
        1848971124000.0*Power(rij,7.0)*Power(xii,7.0) - 
        440561721600.0*Power(rij,8.0)*Power(xii,8.0) - 
        91589097600.0*Power(rij,9.0)*Power(xii,9.0) - 
        16761064320.0*Power(rij,10.0)*Power(xii,10.0) - 
        2717245440.0*Power(rij,11.0)*Power(xii,11.0) - 
        391372800.0*Power(rij,12.0)*Power(xii,12.0) - 
        49889280.0*Power(rij,13.0)*Power(xii,13.0) - 
        5529600.0*Power(rij,14.0)*Power(xii,14.0) - 
        507904.0*Power(rij,15.0)*Power(xii,15.0) - 32768.0*Power(rij,16.0)*Power(xii,16.0))/
      (7.846046208e13*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (42525.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),15.0) + 
        189.0*Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-660.0*Power(rij,5.0)*Power(xii,23.0) - 20.0*Power(rij,6.0)*Power(xii,24.0) + 
           225.0*Power(xij,18.0) + 450.0*rij*xii*Power(xij,18.0) + 
           150.0*rij*Power(xii,3.0)*Power(xij,16.0)*
            (-45.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           225.0*Power(xii,2.0)*Power(xij,16.0)*(-15.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           70.0*Power(rij,4.0)*Power(xii,22.0)*(135.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           20.0*Power(rij,3.0)*Power(xii,21.0)*
            (3675.0 + 341.0*Power(rij,2.0)*Power(xij,2.0)) + 
           30.0*rij*Power(xii,5.0)*Power(xij,14.0)*
            (1575.0 - 150.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           75.0*Power(xii,4.0)*Power(xij,14.0)*
            (315.0 - 90.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 30.0*rij*Power(xii,13.0)*Power(xij,6.0)*
            (-3648435.0 + 137800.0*Power(rij,2.0)*Power(xij,2.0) + 
              156.0*Power(rij,4.0)*Power(xij,4.0)) - 
           10.0*rij*Power(xii,7.0)*Power(xij,12.0)*
            (20475.0 - 1148.0*Power(rij,2.0)*Power(xij,2.0) + 
              178.0*Power(rij,4.0)*Power(xij,4.0)) + 
           120.0*rij*Power(xii,17.0)*Power(xij,2.0)*
            (-132855.0 - 21869.0*Power(rij,2.0)*Power(xij,2.0) + 
              242.0*Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*Power(rij,2.0)*Power(xii,20.0)*
            (-157500.0 - 65525.0*Power(rij,2.0)*Power(xij,2.0) + 
              248.0*Power(rij,4.0)*Power(xij,4.0)) - 
           26.0*rij*Power(xii,9.0)*Power(xij,10.0)*
            (60525.0 + 46500.0*Power(rij,2.0)*Power(xij,2.0) + 
              328.0*Power(rij,4.0)*Power(xij,4.0)) + 
           30.0*rij*Power(xii,11.0)*Power(xij,8.0)*
            (-1302535.0 - 99320.0*Power(rij,2.0)*Power(xij,2.0) + 
              872.0*Power(rij,4.0)*Power(xij,4.0)) - 
           30.0*rij*Power(xii,15.0)*Power(xij,4.0)*
            (2638467.0 - 134540.0*Power(rij,2.0)*Power(xij,2.0) + 
              1716.0*Power(rij,4.0)*Power(xij,4.0)) + 
           4.0*rij*Power(xii,19.0)*(-160650.0 - 322775.0*Power(rij,2.0)*Power(xij,2.0) + 
              2332.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,10.0)*Power(xij,8.0)*
            (-19538025.0 + 11124750.0*Power(rij,2.0)*Power(xij,2.0) + 
              210600.0*Power(rij,4.0)*Power(xij,4.0) - 496.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 20.0*Power(xii,18.0)*(-16065.0 - 335970.0*Power(rij,2.0)*Power(xij,2.0) - 
              2670.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           5.0*Power(xii,6.0)*Power(xij,12.0)*
            (-20475.0 + 9450.0*Power(rij,2.0)*Power(xij,2.0) - 
              450.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           5.0*Power(xii,12.0)*Power(xij,6.0)*
            (10945305.0 - 5223270.0*Power(rij,2.0)*Power(xij,2.0) + 
              100620.0*Power(rij,4.0)*Power(xij,4.0) + 16.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 5.0*Power(xii,8.0)*Power(xij,10.0)*
            (-157365.0 + 62010.0*Power(rij,2.0)*Power(xij,2.0) + 
              13160.0*Power(rij,4.0)*Power(xij,4.0) + 28.0*Power(rij,6.0)*Power(xij,6.0)) - 
           30.0*Power(xii,16.0)*Power(xij,2.0)*
            (265710.0 + 800715.0*Power(rij,2.0)*Power(xij,2.0) - 
              21500.0*Power(rij,4.0)*Power(xij,4.0) + 52.0*Power(rij,6.0)*Power(xij,6.0)) + 
           15.0*Power(xii,14.0)*Power(xij,4.0)*
            (-2638467.0 - 435750.0*Power(rij,2.0)*Power(xij,2.0) - 
              14820.0*Power(rij,4.0)*Power(xij,4.0) + 104.0*Power(rij,6.0)*Power(xij,6.0))) \
    + Power(E,2.0*rij*xii)*Power(xii,8.0)*
         (2.0*Power(xii,2.0)*Power(xij,20.0)*
            (1782492075.0 + 3564984150.0*rij*xij + 
              2364985350.0*Power(rij,2.0)*Power(xij,2.0) + 
              837468450.0*Power(rij,3.0)*Power(xij,3.0) + 
              183599325.0*Power(rij,4.0)*Power(xij,4.0) + 
              25872210.0*Power(rij,5.0)*Power(xij,5.0) + 
              2193030.0*Power(rij,6.0)*Power(xij,6.0) + 
              68220.0*Power(rij,7.0)*Power(xij,7.0) - 
              6885.0*Power(rij,8.0)*Power(xij,8.0) - 890.0*Power(rij,9.0)*Power(xij,9.0) - 
              34.0*Power(rij,10.0)*Power(xij,10.0)) + 
           42.0*Power(xii,4.0)*Power(xij,18.0)*
            (251336925.0 + 502673850.0*rij*xij + 
              209308050.0*Power(rij,2.0)*Power(xij,2.0) + 
              18924570.0*Power(rij,3.0)*Power(xij,3.0) - 
              9849735.0*Power(rij,4.0)*Power(xij,4.0) - 
              3861270.0*Power(rij,5.0)*Power(xij,5.0) - 
              672210.0*Power(rij,6.0)*Power(xij,6.0) - 
              66780.0*Power(rij,7.0)*Power(xij,7.0) - 
              3591.0*Power(rij,8.0)*Power(xij,8.0) - 62.0*Power(rij,9.0)*Power(xij,9.0) + 
              2.0*Power(rij,10.0)*Power(xij,10.0)) + 
           6.0*Power(xij,22.0)*(34459425.0 + 68918850.0*rij*xij + 
              52702650.0*Power(rij,2.0)*Power(xij,2.0) + 
              22972950.0*Power(rij,3.0)*Power(xij,3.0) + 
              6621615.0*Power(rij,4.0)*Power(xij,4.0) + 
              1351350.0*Power(rij,5.0)*Power(xij,5.0) + 
              200970.0*Power(rij,6.0)*Power(xij,6.0) + 
              21780.0*Power(rij,7.0)*Power(xij,7.0) + 
              1665.0*Power(rij,8.0)*Power(xij,8.0) + 82.0*Power(rij,9.0)*Power(xij,9.0) + 
              2.0*Power(rij,10.0)*Power(xij,10.0)) - 
           3.0*Power(xii,22.0)*(14175.0 + 28350.0*rij*xij + 
              28350.0*Power(rij,2.0)*Power(xij,2.0) + 
              18900.0*Power(rij,3.0)*Power(xij,3.0) + 
              9450.0*Power(rij,4.0)*Power(xij,4.0) + 3780.0*Power(rij,5.0)*Power(xij,5.0) + 
              1260.0*Power(rij,6.0)*Power(xij,6.0) + 360.0*Power(rij,7.0)*Power(xij,7.0) + 
              90.0*Power(rij,8.0)*Power(xij,8.0) + 20.0*Power(rij,9.0)*Power(xij,9.0) + 
              4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           21.0*Power(xii,18.0)*Power(xij,4.0)*
            (212625.0 + 425250.0*rij*xij + 425250.0*Power(rij,2.0)*Power(xij,2.0) + 
              283500.0*Power(rij,3.0)*Power(xij,3.0) + 
              141750.0*Power(rij,4.0)*Power(xij,4.0) + 
              56700.0*Power(rij,5.0)*Power(xij,5.0) + 
              18900.0*Power(rij,6.0)*Power(xij,6.0) + 
              4104.0*Power(rij,7.0)*Power(xij,7.0) + 2538.0*Power(rij,8.0)*Power(xij,8.0) + 
              308.0*Power(rij,9.0)*Power(xij,9.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) + 
           54.0*Power(xii,6.0)*Power(xij,16.0)*
            (133451955.0 + 266903910.0*rij*xij - 
              93304890.0*Power(rij,2.0)*Power(xij,2.0) - 
              91581210.0*Power(rij,3.0)*Power(xij,3.0) - 
              19513305.0*Power(rij,4.0)*Power(xij,4.0) - 
              133098.0*Power(rij,5.0)*Power(xij,5.0) + 
              629874.0*Power(rij,6.0)*Power(xij,6.0) + 
              122844.0*Power(rij,7.0)*Power(xij,7.0) + 
              11111.0*Power(rij,8.0)*Power(xij,8.0) + 478.0*Power(rij,9.0)*Power(xij,9.0) + 
              6.0*Power(rij,10.0)*Power(xij,10.0)) - 
           315.0*Power(xii,12.0)*Power(xij,10.0)*
            (-405405.0 - 810810.0*rij*xij - 614250.0*Power(rij,2.0)*Power(xij,2.0) - 
              1408680.0*Power(rij,3.0)*Power(xij,3.0) + 
              573300.0*Power(rij,4.0)*Power(xij,4.0) - 
              153720.0*Power(rij,5.0)*Power(xij,5.0) - 
              133224.0*Power(rij,6.0)*Power(xij,6.0) - 
              20208.0*Power(rij,7.0)*Power(xij,7.0) - 348.0*Power(rij,8.0)*Power(xij,8.0) + 
              136.0*Power(rij,9.0)*Power(xij,9.0) + 8.0*Power(rij,10.0)*Power(xij,10.0)) + 
           315.0*Power(xii,10.0)*Power(xij,12.0)*
            (-482895.0 - 965790.0*rij*xij - 6499710.0*Power(rij,2.0)*Power(xij,2.0) + 
              4684680.0*Power(rij,3.0)*Power(xij,3.0) - 
              380772.0*Power(rij,4.0)*Power(xij,4.0) - 
              912744.0*Power(rij,5.0)*Power(xij,5.0) - 
              188664.0*Power(rij,6.0)*Power(xij,6.0) - 
              7536.0*Power(rij,7.0)*Power(xij,7.0) + 1812.0*Power(rij,8.0)*Power(xij,8.0) + 
              232.0*Power(rij,9.0)*Power(xij,9.0) + 8.0*Power(rij,10.0)*Power(xij,10.0)) - 
           27.0*Power(xii,16.0)*Power(xij,6.0)*
            (-716625.0 - 1433250.0*rij*xij - 1433250.0*Power(rij,2.0)*Power(xij,2.0) - 
              955500.0*Power(rij,3.0)*Power(xij,3.0) - 
              477750.0*Power(rij,4.0)*Power(xij,4.0) - 
              213276.0*Power(rij,5.0)*Power(xij,5.0) - 
              15652.0*Power(rij,6.0)*Power(xij,6.0) - 
              36872.0*Power(rij,7.0)*Power(xij,7.0) - 
              8378.0*Power(rij,8.0)*Power(xij,8.0) - 404.0*Power(rij,9.0)*Power(xij,9.0) + 
              12.0*Power(rij,10.0)*Power(xij,10.0)) + 
           Power(xii,20.0)*Power(xij,2.0)*
            (637875.0 + 1275750.0*rij*xij + 1275750.0*Power(rij,2.0)*Power(xij,2.0) + 
              850500.0*Power(rij,3.0)*Power(xij,3.0) + 
              425250.0*Power(rij,4.0)*Power(xij,4.0) + 
              170100.0*Power(rij,5.0)*Power(xij,5.0) + 
              56700.0*Power(rij,6.0)*Power(xij,6.0) + 
              16200.0*Power(rij,7.0)*Power(xij,7.0) + 
              4050.0*Power(rij,8.0)*Power(xij,8.0) + 1348.0*Power(rij,9.0)*Power(xij,9.0) + 
              68.0*Power(rij,10.0)*Power(xij,10.0)) + 
           3.0*Power(xii,14.0)*Power(xij,8.0)*
            (-19348875.0 - 38697750.0*rij*xij - 
              38697750.0*Power(rij,2.0)*Power(xij,2.0) - 
              24537240.0*Power(rij,3.0)*Power(xij,3.0) - 
              19836180.0*Power(rij,4.0)*Power(xij,4.0) + 
              3197880.0*Power(rij,5.0)*Power(xij,5.0) - 
              3221400.0*Power(rij,6.0)*Power(xij,6.0) - 
              1338000.0*Power(rij,7.0)*Power(xij,7.0) - 
              128196.0*Power(rij,8.0)*Power(xij,8.0) + 
              1208.0*Power(rij,9.0)*Power(xij,9.0) + 472.0*Power(rij,10.0)*Power(xij,10.0)) \
    - 3.0*Power(xii,8.0)*Power(xij,14.0)*
            (-593408025.0 - 1186816050.0*rij*xij + 
              2286535230.0*Power(rij,2.0)*Power(xij,2.0) - 
              157114440.0*Power(rij,3.0)*Power(xij,3.0) - 
              470066940.0*Power(rij,4.0)*Power(xij,4.0) - 
              111426840.0*Power(rij,5.0)*Power(xij,5.0) - 
              5225640.0*Power(rij,6.0)*Power(xij,6.0) + 
              1666032.0*Power(rij,7.0)*Power(xij,7.0) + 
              305964.0*Power(rij,8.0)*Power(xij,8.0) + 
              20504.0*Power(rij,9.0)*Power(xij,9.0) + 472.0*Power(rij,10.0)*Power(xij,10.0)))\
    )/(42525.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),15.0))
    ;
  }
  return S;
}

static double DSlater_3S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-1600593426432000.0 + 1600593426432000.0*Power(E,2.0*rij*xii) - 
        3201186852864000.0*rij*xii - 3201186852864000.0*Power(rij,2.0)*Power(xii,2.0) - 
        2132149009740750.0*Power(rij,3.0)*Power(xii,3.0) - 
        1063111166617500.0*Power(rij,4.0)*Power(xii,4.0) - 
        422807944083525.0*Power(rij,5.0)*Power(xii,5.0) - 
        139509188869050.0*Power(rij,6.0)*Power(xii,6.0) - 
        39204349984800.0*Power(rij,7.0)*Power(xii,7.0) - 
        9554082337800.0*Power(rij,8.0)*Power(xii,8.0) - 
        2044960117200.0*Power(rij,9.0)*Power(xii,9.0) - 
        387930422880.0*Power(rij,10.0)*Power(xii,10.0) - 
        65654184960.0*Power(rij,11.0)*Power(xii,11.0) - 
        9962184960.0*Power(rij,12.0)*Power(xii,12.0) - 
        1359912960.0*Power(rij,13.0)*Power(xii,13.0) - 
        167116800.0*Power(rij,14.0)*Power(xii,14.0) - 
        18382848.0*Power(rij,15.0)*Power(xii,15.0) - 
        1775616.0*Power(rij,16.0)*Power(xii,16.0) - 
        143360.0*Power(rij,17.0)*Power(xii,17.0) - 8192.0*Power(rij,18.0)*Power(xii,18.0))/
      (1.600593426432e15*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (1403325.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),17.0) + 
        10395.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-456.0*Power(rij,5.0)*Power(xii,25.0) - 12.0*Power(rij,6.0)*Power(xii,26.0) + 
           135.0*Power(xij,20.0) + 270.0*rij*xii*Power(xij,20.0) + 
           90.0*rij*Power(xii,3.0)*Power(xij,18.0)*
            (-51.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           135.0*Power(xii,2.0)*Power(xij,18.0)*(-17.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           80.0*Power(rij,4.0)*Power(xii,24.0)*(93.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           20.0*Power(rij,3.0)*Power(xii,23.0)*
            (3264.0 + 395.0*Power(rij,2.0)*Power(xij,2.0)) + 
           36.0*rij*Power(xii,5.0)*Power(xij,16.0)*
            (1020.0 - 85.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           90.0*Power(xii,4.0)*Power(xij,16.0)*
            (204.0 - 51.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           2040.0*rij*Power(xii,13.0)*Power(xij,8.0)*
            (-172098.0 + 715.0*Power(rij,2.0)*Power(xij,2.0) + 
              29.0*Power(rij,4.0)*Power(xij,4.0)) - 
           2040.0*rij*Power(xii,15.0)*Power(xij,6.0)*
            (210168.0 - 8158.0*Power(rij,2.0)*Power(xij,2.0) + 
              39.0*Power(rij,4.0)*Power(xij,4.0)) - 
           20.0*rij*Power(xii,7.0)*Power(xij,14.0)*
            (9180.0 - 132.0*Power(rij,2.0)*Power(xij,2.0) + 
              67.0*Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*Power(rij,2.0)*Power(xii,22.0)*
            (-156060.0 - 80085.0*Power(rij,2.0)*Power(xij,2.0) + 
              94.0*Power(rij,4.0)*Power(xij,4.0)) + 
           68.0*rij*Power(xii,11.0)*Power(xij,10.0)*
            (-1258335.0 - 115770.0*Power(rij,2.0)*Power(xij,2.0) + 
              212.0*Power(rij,4.0)*Power(xij,4.0)) - 
           30.0*rij*Power(xii,17.0)*Power(xij,4.0)*
            (6000651.0 - 10472.0*Power(rij,2.0)*Power(xij,2.0) + 
              816.0*Power(rij,4.0)*Power(xij,4.0)) - 
           4.0*rij*Power(xii,21.0)*(174420.0 + 422805.0*Power(rij,2.0)*Power(xij,2.0) + 
              1399.0*Power(rij,4.0)*Power(xij,4.0)) - 
           4.0*rij*Power(xii,9.0)*Power(xij,12.0)*
            (784125.0 + 415140.0*Power(rij,2.0)*Power(xij,2.0) + 
              3326.0*Power(rij,4.0)*Power(xij,4.0)) + 
           45.0*Power(xii,16.0)*Power(xij,4.0)*
            (-2000217.0 - 1628872.0*Power(rij,2.0)*Power(xij,2.0) + 
              20808.0*Power(rij,4.0)*Power(xij,4.0)) + 
           Power(xii,19.0)*(-24085350.0*rij*Power(xij,2.0) - 
              7125420.0*Power(rij,3.0)*Power(xij,4.0) + 
              59024.0*Power(rij,5.0)*Power(xij,6.0)) + 
           1020.0*Power(xii,14.0)*Power(xij,6.0)*
            (-210168.0 + 38954.0*Power(rij,2.0)*Power(xij,2.0) - 
              1365.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           6.0*Power(xii,6.0)*Power(xij,14.0)*
            (-15300.0 + 6120.0*Power(rij,2.0)*Power(xij,2.0) - 
              255.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) - 
           68.0*Power(xii,12.0)*Power(xij,8.0)*
            (2581470.0 - 1159275.0*Power(rij,2.0)*Power(xij,2.0) + 
              4845.0*Power(rij,4.0)*Power(xij,4.0) + 16.0*Power(rij,6.0)*Power(xij,6.0)) + 
           10.0*Power(xii,8.0)*Power(xij,12.0)*
            (-156825.0 + 53244.0*Power(rij,2.0)*Power(xij,2.0) + 
              6684.0*Power(rij,4.0)*Power(xij,4.0) + 16.0*Power(rij,6.0)*Power(xij,6.0)) - 
           2.0*Power(xii,10.0)*Power(xij,10.0)*
            (21391695.0 - 9981210.0*Power(rij,2.0)*Power(xij,2.0) - 
              221340.0*Power(rij,4.0)*Power(xij,4.0) + 94.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 15.0*Power(xii,18.0)*Power(xij,2.0)*
            (802845.0 + 3733518.0*Power(rij,2.0)*Power(xij,2.0) - 
              56168.0*Power(rij,4.0)*Power(xij,4.0) + 136.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 2.0*Power(xii,20.0)*(-174420.0 - 4738455.0*Power(rij,2.0)*Power(xij,2.0) - 
              198795.0*Power(rij,4.0)*Power(xij,4.0) + 544.0*Power(rij,6.0)*Power(xij,6.0))) \
    + Power(E,2.0*rij*xii)*Power(xii,8.0)*
         (2.0*Power(xii,2.0)*Power(xij,24.0)*
            (218024781975.0 + 436049563950.0*rij*xij + 
              315992927250.0*Power(rij,2.0)*Power(xij,2.0) + 
              127775547900.0*Power(rij,3.0)*Power(xij,3.0) + 
              33563479950.0*Power(rij,4.0)*Power(xij,4.0) + 
              6097020930.0*Power(rij,5.0)*Power(xij,5.0) + 
              773783010.0*Power(rij,6.0)*Power(xij,6.0) + 
              64787580.0*Power(rij,7.0)*Power(xij,7.0) + 
              2627955.0*Power(rij,8.0)*Power(xij,8.0) - 
              117810.0*Power(rij,9.0)*Power(xij,9.0) - 
              25542.0*Power(rij,10.0)*Power(xij,10.0) - 
              1684.0*Power(rij,11.0)*Power(xij,11.0) - 46.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 88.0*Power(xii,20.0)*Power(xij,6.0)*
            (-10843875.0 - 21687750.0*rij*xij - 
              21687750.0*Power(rij,2.0)*Power(xij,2.0) - 
              14458500.0*Power(rij,3.0)*Power(xij,3.0) - 
              7229250.0*Power(rij,4.0)*Power(xij,4.0) - 
              2891700.0*Power(rij,5.0)*Power(xij,5.0) - 
              963900.0*Power(rij,6.0)*Power(xij,6.0) - 
              302130.0*Power(rij,7.0)*Power(xij,7.0) - 
              28755.0*Power(rij,8.0)*Power(xij,8.0) - 
              25380.0*Power(rij,9.0)*Power(xij,9.0) - 
              5805.0*Power(rij,10.0)*Power(xij,10.0) - 
              350.0*Power(rij,11.0)*Power(xij,11.0) + Power(rij,12.0)*Power(xij,12.0)) + 
           88.0*Power(xii,6.0)*Power(xij,20.0)*
            (25012482075.0 + 50024964150.0*rij*xij + 
              3860178525.0*Power(rij,2.0)*Power(xij,2.0) - 
              8806672350.0*Power(rij,3.0)*Power(xij,3.0) - 
              3531006675.0*Power(rij,4.0)*Power(xij,4.0) - 
              529052580.0*Power(rij,5.0)*Power(xij,5.0) + 
              2184840.0*Power(rij,6.0)*Power(xij,6.0) + 
              14577840.0*Power(rij,7.0)*Power(xij,7.0) + 
              2766240.0*Power(rij,8.0)*Power(xij,8.0) + 
              273870.0*Power(rij,9.0)*Power(xij,9.0) + 
              15363.0*Power(rij,10.0)*Power(xij,10.0) + 
              406.0*Power(rij,11.0)*Power(xij,11.0) + Power(rij,12.0)*Power(xij,12.0)) + 
           4488.0*Power(xii,14.0)*Power(xij,12.0)*
            (-3869775.0 - 7739550.0*rij*xij - 8632575.0*Power(rij,2.0)*Power(xij,2.0) - 
              66150.0*Power(rij,3.0)*Power(xij,3.0) - 
              9955575.0*Power(rij,4.0)*Power(xij,4.0) + 
              1446480.0*Power(rij,5.0)*Power(xij,5.0) + 
              379260.0*Power(rij,6.0)*Power(xij,6.0) - 
              283320.0*Power(rij,7.0)*Power(xij,7.0) - 
              88050.0*Power(rij,8.0)*Power(xij,8.0) - 
              8400.0*Power(rij,9.0)*Power(xij,9.0) - 18.0*Power(rij,10.0)*Power(xij,10.0) + 
              44.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) - 
           4488.0*Power(xii,12.0)*Power(xij,14.0)*
            (-6889050.0 - 13778100.0*rij*xij + 20057625.0*Power(rij,2.0)*Power(xij,2.0) - 
              61906950.0*Power(rij,3.0)*Power(xij,3.0) + 
              11911725.0*Power(rij,4.0)*Power(xij,4.0) + 
              5554080.0*Power(rij,5.0)*Power(xij,5.0) - 
              1535940.0*Power(rij,6.0)*Power(xij,6.0) - 
              762660.0*Power(rij,7.0)*Power(xij,7.0) - 
              103200.0*Power(rij,8.0)*Power(xij,8.0) - 
              2820.0*Power(rij,9.0)*Power(xij,9.0) + 
              654.0*Power(rij,10.0)*Power(xij,10.0) + 
              68.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) + 
           6.0*Power(xij,26.0)*(3273645375.0 + 6547290750.0*rij*xij + 
              5306751450.0*Power(rij,2.0)*Power(xij,2.0) + 
              2527024500.0*Power(rij,3.0)*Power(xij,3.0) + 
              817566750.0*Power(rij,4.0)*Power(xij,4.0) + 
              193243050.0*Power(rij,5.0)*Power(xij,5.0) + 
              34684650.0*Power(rij,6.0)*Power(xij,6.0) + 
              4813380.0*Power(rij,7.0)*Power(xij,7.0) + 
              517275.0*Power(rij,8.0)*Power(xij,8.0) + 
              42350.0*Power(rij,9.0)*Power(xij,9.0) + 
              2530.0*Power(rij,10.0)*Power(xij,10.0) + 
              100.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) - 
           3.0*Power(xii,26.0)*(467775.0 + 935550.0*rij*xij + 
              935550.0*Power(rij,2.0)*Power(xij,2.0) + 
              623700.0*Power(rij,3.0)*Power(xij,3.0) + 
              311850.0*Power(rij,4.0)*Power(xij,4.0) + 
              124740.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              11880.0*Power(rij,7.0)*Power(xij,7.0) + 
              2970.0*Power(rij,8.0)*Power(xij,8.0) + 660.0*Power(rij,9.0)*Power(xij,9.0) + 
              132.0*Power(rij,10.0)*Power(xij,10.0) + 
              24.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) + 
           330.0*Power(xii,18.0)*Power(xij,8.0)*
            (-10120950.0 - 20241900.0*rij*xij - 
              20241900.0*Power(rij,2.0)*Power(xij,2.0) - 
              13494600.0*Power(rij,3.0)*Power(xij,3.0) - 
              6747300.0*Power(rij,4.0)*Power(xij,4.0) - 
              2572794.0*Power(rij,5.0)*Power(xij,5.0) - 
              1278018.0*Power(rij,6.0)*Power(xij,6.0) + 
              18180.0*Power(rij,7.0)*Power(xij,7.0) - 
              77691.0*Power(rij,8.0)*Power(xij,8.0) - 
              35630.0*Power(rij,9.0)*Power(xij,9.0) - 
              4114.0*Power(rij,10.0)*Power(xij,10.0) - 
              92.0*Power(rij,11.0)*Power(xij,11.0) + 6.0*Power(rij,12.0)*Power(xij,12.0)) - 
           165.0*Power(xii,8.0)*Power(xij,18.0)*
            (-5877371745.0 - 11754743490.0*rij*xij + 
              10638195750.0*Power(rij,2.0)*Power(xij,2.0) + 
              3500767620.0*Power(rij,3.0)*Power(xij,3.0) - 
              943138350.0*Power(rij,4.0)*Power(xij,4.0) - 
              587351268.0*Power(rij,5.0)*Power(xij,5.0) - 
              104134716.0*Power(rij,6.0)*Power(xij,6.0) - 
              4631400.0*Power(rij,7.0)*Power(xij,7.0) + 
              1221198.0*Power(rij,8.0)*Power(xij,8.0) + 
              247260.0*Power(rij,9.0)*Power(xij,9.0) + 
              20892.0*Power(rij,10.0)*Power(xij,10.0) + 
              856.0*Power(rij,11.0)*Power(xij,11.0) + 12.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 8.0*Power(xii,4.0)*Power(xij,22.0)*
            (230263179300.0 + 460526358600.0*rij*xij + 
              246912383025.0*Power(rij,2.0)*Power(xij,2.0) + 
              55394851050.0*Power(rij,3.0)*Power(xij,3.0) + 
              1508262525.0*Power(rij,4.0)*Power(xij,4.0) - 
              2396962260.0*Power(rij,5.0)*Power(xij,5.0) - 
              722410920.0*Power(rij,6.0)*Power(xij,6.0) - 
              115167690.0*Power(rij,7.0)*Power(xij,7.0) - 
              11586465.0*Power(rij,8.0)*Power(xij,8.0) - 
              720720.0*Power(rij,9.0)*Power(xij,9.0) - 
              21681.0*Power(rij,10.0)*Power(xij,10.0) + 
              218.0*Power(rij,11.0)*Power(xij,11.0) + 29.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 8.0*Power(xii,22.0)*Power(xij,4.0)*
            (23856525.0 + 47713050.0*rij*xij + 47713050.0*Power(rij,2.0)*Power(xij,2.0) + 
              31808700.0*Power(rij,3.0)*Power(xij,3.0) + 
              15904350.0*Power(rij,4.0)*Power(xij,4.0) + 
              6361740.0*Power(rij,5.0)*Power(xij,5.0) + 
              2120580.0*Power(rij,6.0)*Power(xij,6.0) + 
              605880.0*Power(rij,7.0)*Power(xij,7.0) + 
              151470.0*Power(rij,8.0)*Power(xij,8.0) + 
              26730.0*Power(rij,9.0)*Power(xij,9.0) + 
              11583.0*Power(rij,10.0)*Power(xij,10.0) + 
              1406.0*Power(rij,11.0)*Power(xij,11.0) + 29.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 66.0*Power(xii,16.0)*Power(xij,10.0)*
            (-131572350.0 - 263144700.0*rij*xij - 
              263144700.0*Power(rij,2.0)*Power(xij,2.0) - 
              178869600.0*Power(rij,3.0)*Power(xij,3.0) - 
              63636300.0*Power(rij,4.0)*Power(xij,4.0) - 
              76650210.0*Power(rij,5.0)*Power(xij,5.0) + 
              7510230.0*Power(rij,6.0)*Power(xij,6.0) - 
              1096020.0*Power(rij,7.0)*Power(xij,7.0) - 
              2359515.0*Power(rij,8.0)*Power(xij,8.0) - 
              455070.0*Power(rij,9.0)*Power(xij,9.0) - 
              25722.0*Power(rij,10.0)*Power(xij,10.0) + 
              716.0*Power(rij,11.0)*Power(xij,11.0) + 86.0*Power(rij,12.0)*Power(xij,12.0)) \
    + Power(xii,24.0)*Power(xij,2.0)*(23856525.0 + 47713050.0*rij*xij + 
              47713050.0*Power(rij,2.0)*Power(xij,2.0) + 
              31808700.0*Power(rij,3.0)*Power(xij,3.0) + 
              15904350.0*Power(rij,4.0)*Power(xij,4.0) + 
              6361740.0*Power(rij,5.0)*Power(xij,5.0) + 
              2120580.0*Power(rij,6.0)*Power(xij,6.0) + 
              605880.0*Power(rij,7.0)*Power(xij,7.0) + 
              151470.0*Power(rij,8.0)*Power(xij,8.0) + 
              33660.0*Power(rij,9.0)*Power(xij,9.0) + 
              6732.0*Power(rij,10.0)*Power(xij,10.0) + 
              1784.0*Power(rij,11.0)*Power(xij,11.0) + 92.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 33.0*Power(xii,10.0)*Power(xij,16.0)*
            (2759659875.0 + 5519319750.0*rij*xij - 
              33545959650.0*Power(rij,2.0)*Power(xij,2.0) + 
              10690917300.0*Power(rij,3.0)*Power(xij,3.0) + 
              5050165050.0*Power(rij,4.0)*Power(xij,4.0) - 
              924808500.0*Power(rij,5.0)*Power(xij,5.0) - 
              639355500.0*Power(rij,6.0)*Power(xij,6.0) - 
              107153640.0*Power(rij,7.0)*Power(xij,7.0) - 
              4761930.0*Power(rij,8.0)*Power(xij,8.0) + 
              823020.0*Power(rij,9.0)*Power(xij,9.0) + 
              138060.0*Power(rij,10.0)*Power(xij,10.0) + 
              8200.0*Power(rij,11.0)*Power(xij,11.0) + 172.0*Power(rij,12.0)*Power(xij,12.0))\
    ))/(1.403325e6*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),17.0))
    ;
  }
  return S;
}

double DSlater_3S_1S(double rij,double xii,double xij)
{
  return DSlater_1S_3S(rij,xij,xii);
}

double DSlater_3S_2S(double rij,double xii,double xij)
{
  return DSlater_2S_3S(rij,xij,xii);
}

static double DSlater_4S_4S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-20922789888000.0 + 20922789888000.0*Power(E,2.0*rij*xii) - 
        41845579776000.0*rij*xii - 41845579776000.0*Power(rij,2.0)*Power(xii,2.0) - 
        27805745842875.0*Power(rij,3.0)*Power(xii,3.0) - 
        13765911909750.0*Power(rij,4.0)*Power(xii,4.0) - 
        5399605411200.0*Power(rij,5.0)*Power(xii,5.0) - 
        1743679337400.0*Power(rij,6.0)*Power(xii,6.0) - 
        476010334800.0*Power(rij,7.0)*Power(xii,7.0) - 
        112021509600.0*Power(rij,8.0)*Power(xii,8.0) - 
        23063040000.0*Power(rij,9.0)*Power(xii,9.0) - 
        4197473280.0*Power(rij,10.0)*Power(xii,10.0) - 
        679311360.0*Power(rij,11.0)*Power(xii,11.0) - 
        97843200.0*Power(rij,12.0)*Power(xii,12.0) - 
        12472320.0*Power(rij,13.0)*Power(xii,13.0) - 
        1382400.0*Power(rij,14.0)*Power(xii,14.0) - 
        126976.0*Power(rij,15.0)*Power(xii,15.0) - 8192.0*Power(rij,16.0)*Power(xii,16.0))/
      (2.0922789888e13*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (315.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),15.0) + 
        Power(E,2.0*rij*xij)*Power(xij,10.0)*
         (-1428.0*Power(rij,6.0)*Power(xii,26.0) - 78.0*Power(rij,7.0)*Power(xii,27.0) - 
           2.0*Power(rij,8.0)*Power(xii,28.0) + 315.0*Power(xij,20.0) + 
           630.0*rij*xii*Power(xij,20.0) + 
           210.0*rij*Power(xii,3.0)*Power(xij,18.0)*
            (-45.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           315.0*Power(xii,2.0)*Power(xij,18.0)*(-15.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           42.0*Power(rij,5.0)*Power(xii,25.0)*(377.0 + 5.0*Power(rij,2.0)*Power(xij,2.0)) + 
           42.0*Power(rij,4.0)*Power(xii,24.0)*
            (-2730.0 - 190.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) \
    + 42.0*rij*Power(xii,5.0)*Power(xij,16.0)*
            (1575.0 - 150.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           105.0*Power(xii,4.0)*Power(xij,16.0)*
            (315.0 - 90.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 63.0*Power(rij,3.0)*Power(xii,23.0)*
            (-8645.0 - 2180.0*Power(rij,2.0)*Power(xij,2.0) + 
              32.0*Power(rij,4.0)*Power(xij,4.0)) + 
           2.0*rij*Power(xii,7.0)*Power(xij,14.0)*
            (-143325.0 + 22050.0*Power(rij,2.0)*Power(xij,2.0) - 
              630.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           7.0*Power(xii,6.0)*Power(xij,14.0)*
            (-20475.0 + 9450.0*Power(rij,2.0)*Power(xij,2.0) - 
              450.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           21.0*rij*Power(xii,11.0)*Power(xij,10.0)*
            (-209430.0 - 112125.0*Power(rij,2.0)*Power(xij,2.0) - 
              8288.0*Power(rij,4.0)*Power(xij,4.0) + 10.0*Power(rij,6.0)*Power(xij,6.0)) - 
           21.0*rij*Power(xii,9.0)*Power(xij,12.0)*
            (-40950.0 + 11245.0*Power(rij,2.0)*Power(xij,2.0) - 
              222.0*Power(rij,4.0)*Power(xij,4.0) + 10.0*Power(rij,6.0)*Power(xij,6.0)) + 
           7.0*rij*Power(xii,19.0)*Power(xij,2.0)*
            (-7711200.0 - 1605825.0*Power(rij,2.0)*Power(xij,2.0) + 
              55104.0*Power(rij,4.0)*Power(xij,4.0) + 20.0*Power(rij,6.0)*Power(xij,6.0)) - 
           4.0*Power(rij,2.0)*Power(xii,22.0)*
            (401310.0 + 341775.0*Power(rij,2.0)*Power(xij,2.0) - 
              9009.0*Power(rij,4.0)*Power(xij,4.0) + 32.0*Power(rij,6.0)*Power(xij,6.0)) + 
           105.0*rij*Power(xii,17.0)*Power(xij,4.0)*
            (-2087532.0 + 267621.0*Power(rij,2.0)*Power(xij,2.0) - 
              10348.0*Power(rij,4.0)*Power(xij,4.0) + 52.0*Power(rij,6.0)*Power(xij,6.0)) - 
           105.0*rij*Power(xii,15.0)*Power(xij,6.0)*
            (2126142.0 - 103075.0*Power(rij,2.0)*Power(xij,2.0) - 
              4680.0*Power(rij,4.0)*Power(xij,4.0) + 56.0*Power(rij,6.0)*Power(xij,6.0)) + 
           21.0*Power(xii,10.0)*Power(xij,10.0)*
            (-104715.0 + 83850.0*Power(rij,2.0)*Power(xij,2.0) + 
              4030.0*Power(rij,4.0)*Power(xij,4.0) + 404.0*Power(rij,6.0)*Power(xij,6.0)) - 
           70.0*Power(xii,18.0)*Power(xij,2.0)*
            (385560.0 + 1201608.0*Power(rij,2.0)*Power(xij,2.0) - 
              84195.0*Power(rij,4.0)*Power(xij,4.0) + 1064.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 3.0*rij*Power(xii,21.0)*(835380.0 + 2774625.0*Power(rij,2.0)*Power(xij,2.0) - 
              94836.0*Power(rij,4.0)*Power(xij,4.0) + 1160.0*Power(rij,6.0)*Power(xij,6.0)) \
    + rij*Power(xii,13.0)*Power(xij,8.0)*
            (-50825250.0 - 16261245.0*Power(rij,2.0)*Power(xij,2.0) + 
              248640.0*Power(rij,4.0)*Power(xij,4.0) + 2024.0*Power(rij,6.0)*Power(xij,6.0)\
    ) - 70.0*Power(xii,16.0)*Power(xij,4.0)*
            (1565649.0 - 145035.0*Power(rij,2.0)*Power(xij,2.0) + 
              63465.0*Power(rij,4.0)*Power(xij,4.0) - 
              1560.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xii,8.0)*Power(xij,12.0)*
            (429975.0 - 286650.0*Power(rij,2.0)*Power(xij,2.0) + 
              22050.0*Power(rij,4.0)*Power(xij,4.0) - 
              420.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           7.0*Power(xii,12.0)*Power(xij,8.0)*
            (3630375.0 - 2811510.0*Power(rij,2.0)*Power(xij,2.0) - 
              298350.0*Power(rij,4.0)*Power(xij,4.0) + 
              1688.0*Power(rij,6.0)*Power(xij,6.0) + 6.0*Power(rij,8.0)*Power(xij,8.0)) + 
           14.0*Power(xii,20.0)*(-89505.0 - 2135700.0*Power(rij,2.0)*Power(xij,2.0) + 
              24030.0*Power(rij,4.0)*Power(xij,4.0) - 
              1236.0*Power(rij,6.0)*Power(xij,6.0) + 10.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xii,14.0)*Power(xij,6.0)*
            (-111622455.0 + 84253050.0*Power(rij,2.0)*Power(xij,2.0) - 
              2497950.0*Power(rij,4.0)*Power(xij,4.0) - 
              40320.0*Power(rij,6.0)*Power(xij,6.0) + 128.0*Power(rij,8.0)*Power(xij,8.0))) \
    + Power(E,2.0*rij*xii)*Power(xii,10.0)*
         (105.0*Power(xii,18.0)*Power(xij,2.0)*
            (45.0 + 90.0*rij*xij + 90.0*Power(rij,2.0)*Power(xij,2.0) + 
              60.0*Power(rij,3.0)*Power(xij,3.0) + 30.0*Power(rij,4.0)*Power(xij,4.0) + 
              12.0*Power(rij,5.0)*Power(xij,5.0) + 4.0*Power(rij,6.0)*Power(xij,6.0) + 
              2.0*Power(rij,7.0)*Power(xij,7.0)) + 
           105.0*Power(xii,2.0)*Power(xij,18.0)*
            (257040.0 + 514080.0*rij*xij + 284760.0*Power(rij,2.0)*Power(xij,2.0) + 
              79275.0*Power(rij,3.0)*Power(xij,3.0) + 
              13020.0*Power(rij,4.0)*Power(xij,4.0) + 
              1308.0*Power(rij,5.0)*Power(xij,5.0) + 76.0*Power(rij,6.0)*Power(xij,6.0) + 
              2.0*Power(rij,7.0)*Power(xij,7.0)) - 
           1365.0*Power(xii,10.0)*Power(xij,10.0)*
            (-1611.0 - 3222.0*rij*xij + 14418.0*Power(rij,2.0)*Power(xij,2.0) - 
              11913.0*Power(rij,3.0)*Power(xij,3.0) - 
              1830.0*Power(rij,4.0)*Power(xij,4.0) + 360.0*Power(rij,5.0)*Power(xij,5.0) + 
              80.0*Power(rij,6.0)*Power(xij,6.0) + 4.0*Power(rij,7.0)*Power(xij,7.0)) + 
           Power(xii,14.0)*Power(xij,6.0)*
            (143325.0 + 286650.0*rij*xij + 286650.0*Power(rij,2.0)*Power(xij,2.0) + 
              236145.0*Power(rij,3.0)*Power(xij,3.0) - 
              84630.0*Power(rij,4.0)*Power(xij,4.0) + 
              174048.0*Power(rij,5.0)*Power(xij,5.0) + 
              11816.0*Power(rij,6.0)*Power(xij,6.0) - 
              2024.0*Power(rij,7.0)*Power(xij,7.0) - 128.0*Power(rij,8.0)*Power(xij,8.0)) + 
           21.0*Power(xii,16.0)*Power(xij,4.0)*
            (-1575.0 - 3150.0*rij*xij - 3150.0*Power(rij,2.0)*Power(xij,2.0) - 
              2100.0*Power(rij,3.0)*Power(xij,3.0) - 1050.0*Power(rij,4.0)*Power(xij,4.0) - 
              222.0*Power(rij,5.0)*Power(xij,5.0) - 404.0*Power(rij,6.0)*Power(xij,6.0) - 
              10.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           Power(xii,20.0)*(315.0 + 630.0*rij*xij + 630.0*Power(rij,2.0)*Power(xij,2.0) + 
              420.0*Power(rij,3.0)*Power(xij,3.0) + 210.0*Power(rij,4.0)*Power(xij,4.0) + 
              84.0*Power(rij,5.0)*Power(xij,5.0) + 28.0*Power(rij,6.0)*Power(xij,6.0) + 
              8.0*Power(rij,7.0)*Power(xij,7.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xij,20.0)*(1253070.0 + 2506140.0*rij*xij + 
              1605240.0*Power(rij,2.0)*Power(xij,2.0) + 
              544635.0*Power(rij,3.0)*Power(xij,3.0) + 
              114660.0*Power(rij,4.0)*Power(xij,4.0) + 
              15834.0*Power(rij,5.0)*Power(xij,5.0) + 
              1428.0*Power(rij,6.0)*Power(xij,6.0) + 78.0*Power(rij,7.0)*Power(xij,7.0) + 
              2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           21.0*Power(xii,4.0)*Power(xij,16.0)*
            (-5218830.0 - 10437660.0*rij*xij - 4005360.0*Power(rij,2.0)*Power(xij,2.0) - 
              535275.0*Power(rij,3.0)*Power(xij,3.0) + 
              16020.0*Power(rij,4.0)*Power(xij,4.0) + 
              13548.0*Power(rij,5.0)*Power(xij,5.0) + 
              1716.0*Power(rij,6.0)*Power(xij,6.0) + 96.0*Power(rij,7.0)*Power(xij,7.0) + 
              2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           35.0*Power(xii,8.0)*Power(xij,12.0)*
            (-726075.0 - 1452150.0*rij*xij + 2407230.0*Power(rij,2.0)*Power(xij,2.0) + 
              309225.0*Power(rij,3.0)*Power(xij,3.0) - 
              126930.0*Power(rij,4.0)*Power(xij,4.0) - 
              31044.0*Power(rij,5.0)*Power(xij,5.0) - 
              2128.0*Power(rij,6.0)*Power(xij,6.0) + 4.0*Power(rij,7.0)*Power(xij,7.0) + 
              4.0*Power(rij,8.0)*Power(xij,8.0)) + 
           35.0*Power(xii,12.0)*Power(xij,8.0)*
            (-12285.0 - 24570.0*rij*xij - 50310.0*Power(rij,2.0)*Power(xij,2.0) + 
              67275.0*Power(rij,3.0)*Power(xij,3.0) - 
              59670.0*Power(rij,4.0)*Power(xij,4.0) - 
              7104.0*Power(rij,5.0)*Power(xij,5.0) + 1152.0*Power(rij,6.0)*Power(xij,6.0) + 
              168.0*Power(rij,7.0)*Power(xij,7.0) + 4.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xii,6.0)*Power(xij,14.0)*
            (111622455.0 + 223244910.0*rij*xij - 
              10152450.0*Power(rij,2.0)*Power(xij,2.0) - 
              28100205.0*Power(rij,3.0)*Power(xij,3.0) - 
              5893650.0*Power(rij,4.0)*Power(xij,4.0) - 
              385728.0*Power(rij,5.0)*Power(xij,5.0) + 
              17304.0*Power(rij,6.0)*Power(xij,6.0) + 3480.0*Power(rij,7.0)*Power(xij,7.0) + 
              128.0*Power(rij,8.0)*Power(xij,8.0))))/
      (315.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),15.0))
    ;
  }
  return S;
}

static double DSlater_4S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-1778437140480000.0 + 1778437140480000.0*Power(E,2.0*rij*xii) - 
        3556874280960000.0*rij*xii - 3556874280960000.0*Power(rij,2.0)*Power(xii,2.0) - 
        2366075437976250.0*Power(rij,3.0)*Power(xii,3.0) - 
        1175276594992500.0*Power(rij,4.0)*Power(xii,4.0) - 
        464005220453775.0*Power(rij,5.0)*Power(xii,5.0) - 
        151391487797550.0*Power(rij,6.0)*Power(xii,6.0) - 
        41921958078000.0*Power(rij,7.0)*Power(xii,7.0) - 
        10045335900600.0*Power(rij,8.0)*Power(xii,8.0) - 
        2113817706000.0*Power(rij,9.0)*Power(xii,9.0) - 
        395085731040.0*Power(rij,10.0)*Power(xii,10.0) - 
        66153185280.0*Power(rij,11.0)*Power(xii,11.0) - 
        9980006400.0*Power(rij,12.0)*Power(xii,12.0) - 
        1359912960.0*Power(rij,13.0)*Power(xii,13.0) - 
        167116800.0*Power(rij,14.0)*Power(xii,14.0) - 
        18382848.0*Power(rij,15.0)*Power(xii,15.0) - 
        1775616.0*Power(rij,16.0)*Power(xii,16.0) - 
        143360.0*Power(rij,17.0)*Power(xii,17.0) - 8192.0*Power(rij,18.0)*Power(xii,18.0))/
      (1.77843714048e15*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (14175.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),17.0) + 
        9.0*Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-460.0*Power(rij,7.0)*Power(xii,29.0) - 10.0*Power(rij,8.0)*Power(xii,30.0) + 
           1575.0*Power(xij,22.0) + 3150.0*rij*xii*Power(xij,22.0) - 
           50.0*Power(rij,6.0)*Power(xii,28.0)*(196.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           1050.0*rij*Power(xii,3.0)*Power(xij,20.0)*
            (-51.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           1575.0*Power(xii,2.0)*Power(xij,20.0)*
            (-17.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           30.0*Power(rij,5.0)*Power(xii,27.0)*
            (4165.0 + 128.0*Power(rij,2.0)*Power(xij,2.0)) + 
           420.0*rij*Power(xii,5.0)*Power(xij,18.0)*
            (1020.0 - 85.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           1050.0*Power(xii,4.0)*Power(xij,18.0)*
            (204.0 - 51.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           42.0*Power(rij,4.0)*Power(xii,26.0)*
            (-24500.0 - 2780.0*Power(rij,2.0)*Power(xij,2.0) + 
              9.0*Power(rij,4.0)*Power(xij,4.0)) + 
           210.0*Power(rij,3.0)*Power(xii,25.0)*
            (-26180.0 - 9333.0*Power(rij,2.0)*Power(xij,2.0) + 
              74.0*Power(rij,4.0)*Power(xij,4.0)) + 
           20.0*rij*Power(xii,7.0)*Power(xij,16.0)*
            (-107100.0 + 14280.0*Power(rij,2.0)*Power(xij,2.0) - 
              357.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           70.0*Power(xii,6.0)*Power(xij,16.0)*
            (-15300.0 + 6120.0*Power(rij,2.0)*Power(xij,2.0) - 
              255.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) - 
           42.0*rij*Power(xii,11.0)*Power(xij,12.0)*
            (1723800.0 + 1096500.0*Power(rij,2.0)*Power(xij,2.0) + 
              52703.0*Power(rij,4.0)*Power(xij,4.0) + 80.0*Power(rij,6.0)*Power(xij,6.0)) - 
           340.0*rij*Power(xii,15.0)*Power(xij,8.0)*
            (23752701.0 + 491820.0*Power(rij,2.0)*Power(xij,2.0) - 
              44772.0*Power(rij,4.0)*Power(xij,4.0) + 116.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 2.0*Power(rij,2.0)*Power(xii,24.0)*
            (8996400.0 + 10103100.0*Power(rij,2.0)*Power(xij,2.0) - 
              102522.0*Power(rij,4.0)*Power(xij,4.0) + 263.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 70.0*rij*Power(xii,17.0)*Power(xij,6.0)*
            (180826281.0 - 14550096.0*Power(rij,2.0)*Power(xij,2.0) + 
              184314.0*Power(rij,4.0)*Power(xij,4.0) + 340.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 68.0*rij*Power(xii,13.0)*Power(xij,10.0)*
            (-20763225.0 - 6046950.0*Power(rij,2.0)*Power(xij,2.0) - 
              9324.0*Power(rij,4.0)*Power(xij,4.0) + 425.0*Power(rij,6.0)*Power(xij,6.0)) + 
           4.0*rij*Power(xii,23.0)*(-7630875.0 - 
              32594100.0*Power(rij,2.0)*Power(xij,2.0) - 
              71127.0*Power(rij,4.0)*Power(xij,4.0) + 650.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 2.0*rij*Power(xii,9.0)*Power(xij,14.0)*
            (-3748500.0 + 1014300.0*Power(rij,2.0)*Power(xij,2.0) - 
              7539.0*Power(rij,4.0)*Power(xij,4.0) + 670.0*Power(rij,6.0)*Power(xij,6.0)) + 
           14.0*rij*Power(xii,19.0)*Power(xij,4.0)*
            (-451386135.0 + 14910390.0*Power(rij,2.0)*Power(xij,2.0) - 
              777342.0*Power(rij,4.0)*Power(xij,4.0) + 6800.0*Power(rij,6.0)*Power(xij,6.0)\
    ) - 4.0*rij*Power(xii,21.0)*Power(xij,2.0)*
            (239609475.0 + 116284245.0*Power(rij,2.0)*Power(xij,2.0) - 
              3442719.0*Power(rij,4.0)*Power(xij,4.0) + 
              17510.0*Power(rij,6.0)*Power(xij,6.0)) + 
           10.0*Power(xii,8.0)*Power(xij,14.0)*
            (374850.0 - 214200.0*Power(rij,2.0)*Power(xij,2.0) + 
              14280.0*Power(rij,4.0)*Power(xij,4.0) - 
              238.0*Power(rij,6.0)*Power(xij,6.0) + Power(rij,8.0)*Power(xij,8.0)) + 
           2.0*Power(xii,10.0)*Power(xij,12.0)*
            (-18099900.0 + 11406150.0*Power(rij,2.0)*Power(xij,2.0) + 
              844200.0*Power(rij,4.0)*Power(xij,4.0) + 
              37548.0*Power(rij,6.0)*Power(xij,6.0) + 25.0*Power(rij,8.0)*Power(xij,8.0)) - 
           14.0*Power(xii,12.0)*Power(xij,10.0)*
            (50424975.0 - 32917950.0*Power(rij,2.0)*Power(xij,2.0) - 
              2731050.0*Power(rij,4.0)*Power(xij,4.0) - 
              5212.0*Power(rij,6.0)*Power(xij,6.0) + 27.0*Power(rij,8.0)*Power(xij,8.0)) - 
           35.0*Power(xii,18.0)*Power(xij,4.0)*
            (90277227.0 + 71255790.0*Power(rij,2.0)*Power(xij,2.0) - 
              1723800.0*Power(rij,4.0)*Power(xij,4.0) - 
              16864.0*Power(rij,6.0)*Power(xij,6.0) + 68.0*Power(rij,8.0)*Power(xij,8.0)) + 
           14.0*Power(xii,20.0)*Power(xij,2.0)*
            (-34229925.0 - 184803525.0*Power(rij,2.0)*Power(xij,2.0) + 
              8158275.0*Power(rij,4.0)*Power(xij,4.0) - 
              116144.0*Power(rij,6.0)*Power(xij,6.0) + 170.0*Power(rij,8.0)*Power(xij,8.0)) \
    + 2.0*Power(xii,14.0)*Power(xij,8.0)*
            (-2018979585.0 + 1345818600.0*Power(rij,2.0)*Power(xij,2.0) + 
              3016650.0*Power(rij,4.0)*Power(xij,4.0) - 
              446012.0*Power(rij,6.0)*Power(xij,6.0) + 263.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 2.0*Power(xii,22.0)*(7630875.0 + 249971400.0*Power(rij,2.0)*Power(xij,2.0) + 
              19123125.0*Power(rij,4.0)*Power(xij,4.0) - 
              325766.0*Power(rij,6.0)*Power(xij,6.0) + 306.0*Power(rij,8.0)*Power(xij,8.0)) \
    + Power(xii,16.0)*Power(xij,6.0)*(-6328919835.0 + 
              2425600800.0*Power(rij,2.0)*Power(xij,2.0) - 
              161149800.0*Power(rij,4.0)*Power(xij,4.0) + 
              1051960.0*Power(rij,6.0)*Power(xij,6.0) + 612.0*Power(rij,8.0)*Power(xij,8.0))\
    ) + Power(E,2.0*rij*xii)*Power(xii,10.0)*
         (10710.0*Power(xii,12.0)*Power(xij,12.0)*
            (-3555.0 - 7110.0*rij*xij - 392400.0*Power(rij,2.0)*Power(xij,2.0) + 
              425880.0*Power(rij,3.0)*Power(xij,3.0) - 
              82260.0*Power(rij,4.0)*Power(xij,4.0) - 
              46200.0*Power(rij,5.0)*Power(xij,5.0) - 
              2064.0*Power(rij,6.0)*Power(xij,6.0) + 792.0*Power(rij,7.0)*Power(xij,7.0) + 
              106.0*Power(rij,8.0)*Power(xij,8.0) + 4.0*Power(rij,9.0)*Power(xij,9.0)) + 
           2.0*Power(xii,18.0)*Power(xij,6.0)*
            (4819500.0 + 9639000.0*rij*xij + 9639000.0*Power(rij,2.0)*Power(xij,2.0) + 
              6426000.0*Power(rij,3.0)*Power(xij,3.0) + 
              3213000.0*Power(rij,4.0)*Power(xij,4.0) + 
              1609524.0*Power(rij,5.0)*Power(xij,5.0) - 
              274302.0*Power(rij,6.0)*Power(xij,6.0) + 
              434844.0*Power(rij,7.0)*Power(xij,7.0) + 
              50499.0*Power(rij,8.0)*Power(xij,8.0) - 
              1762.0*Power(rij,9.0)*Power(xij,9.0) - 212.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 2.0*Power(xii,2.0)*Power(xij,22.0)*
            (3264488325.0 + 6528976650.0*rij*xij + 
              4129347600.0*Power(rij,2.0)*Power(xij,2.0) + 
              1391229000.0*Power(rij,3.0)*Power(xij,3.0) + 
              294632100.0*Power(rij,4.0)*Power(xij,4.0) + 
              41833260.0*Power(rij,5.0)*Power(xij,5.0) + 
              4026330.0*Power(rij,6.0)*Power(xij,6.0) + 
              251820.0*Power(rij,7.0)*Power(xij,7.0) + 
              8775.0*Power(rij,8.0)*Power(xij,8.0) + 70.0*Power(rij,9.0)*Power(xij,9.0) - 
              4.0*Power(rij,10.0)*Power(xij,10.0)) + 
           4.0*Power(xij,24.0)*(59520825.0 + 119041650.0*rij*xij + 
              84582225.0*Power(rij,2.0)*Power(xij,2.0) + 
              33415200.0*Power(rij,3.0)*Power(xij,3.0) + 
              8599500.0*Power(rij,4.0)*Power(xij,4.0) + 
              1547910.0*Power(rij,5.0)*Power(xij,5.0) + 
              200655.0*Power(rij,6.0)*Power(xij,6.0) + 
              18720.0*Power(rij,7.0)*Power(xij,7.0) + 
              1215.0*Power(rij,8.0)*Power(xij,8.0) + 50.0*Power(rij,9.0)*Power(xij,9.0) + 
              Power(rij,10.0)*Power(xij,10.0)) - 
           Power(xii,24.0)*(14175.0 + 28350.0*rij*xij + 
              28350.0*Power(rij,2.0)*Power(xij,2.0) + 
              18900.0*Power(rij,3.0)*Power(xij,3.0) + 
              9450.0*Power(rij,4.0)*Power(xij,4.0) + 3780.0*Power(rij,5.0)*Power(xij,5.0) + 
              1260.0*Power(rij,6.0)*Power(xij,6.0) + 360.0*Power(rij,7.0)*Power(xij,7.0) + 
              90.0*Power(rij,8.0)*Power(xij,8.0) + 20.0*Power(rij,9.0)*Power(xij,9.0) + 
              4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           204.0*Power(xii,14.0)*Power(xij,10.0)*
            (-429975.0 - 859950.0*rij*xij - 307125.0*Power(rij,2.0)*Power(xij,2.0) - 
              3071250.0*Power(rij,3.0)*Power(xij,3.0) + 
              2395575.0*Power(rij,4.0)*Power(xij,4.0) - 
              620802.0*Power(rij,5.0)*Power(xij,5.0) - 
              234234.0*Power(rij,6.0)*Power(xij,6.0) - 
              5004.0*Power(rij,7.0)*Power(xij,7.0) + 2949.0*Power(rij,8.0)*Power(xij,8.0) + 
              242.0*Power(rij,9.0)*Power(xij,9.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) + 
           102.0*Power(xii,10.0)*Power(xij,14.0)*
            (44986725.0 + 89973450.0*rij*xij - 
              239334480.0*Power(rij,2.0)*Power(xij,2.0) + 
              57221640.0*Power(rij,3.0)*Power(xij,3.0) + 
              33086340.0*Power(rij,4.0)*Power(xij,4.0) + 
              1567440.0*Power(rij,5.0)*Power(xij,5.0) - 
              825300.0*Power(rij,6.0)*Power(xij,6.0) - 
              141912.0*Power(rij,7.0)*Power(xij,7.0) - 
              8094.0*Power(rij,8.0)*Power(xij,8.0) - 44.0*Power(rij,9.0)*Power(xij,9.0) + 
              8.0*Power(rij,10.0)*Power(xij,10.0)) + 
           Power(xii,22.0)*Power(xij,2.0)*
            (240975.0 + 481950.0*rij*xij + 481950.0*Power(rij,2.0)*Power(xij,2.0) + 
              321300.0*Power(rij,3.0)*Power(xij,3.0) + 
              160650.0*Power(rij,4.0)*Power(xij,4.0) + 
              64260.0*Power(rij,5.0)*Power(xij,5.0) + 
              21420.0*Power(rij,6.0)*Power(xij,6.0) + 
              6120.0*Power(rij,7.0)*Power(xij,7.0) + 1530.0*Power(rij,8.0)*Power(xij,8.0) + 
              580.0*Power(rij,9.0)*Power(xij,9.0) + 8.0*Power(rij,10.0)*Power(xij,10.0)) + 
           2.0*Power(xii,20.0)*Power(xij,4.0)*
            (-963900.0 - 1927800.0*rij*xij - 1927800.0*Power(rij,2.0)*Power(xij,2.0) - 
              1285200.0*Power(rij,3.0)*Power(xij,3.0) - 
              642600.0*Power(rij,4.0)*Power(xij,4.0) - 
              257040.0*Power(rij,5.0)*Power(xij,5.0) - 
              85680.0*Power(rij,6.0)*Power(xij,6.0) - 
              13788.0*Power(rij,7.0)*Power(xij,7.0) - 
              15921.0*Power(rij,8.0)*Power(xij,8.0) - 826.0*Power(rij,9.0)*Power(xij,9.0) + 
              40.0*Power(rij,10.0)*Power(xij,10.0)) - 
           2.0*Power(xii,4.0)*Power(xij,20.0)*
            (-18032978565.0 - 36065957130.0*rij*xij - 
              17600042880.0*Power(rij,2.0)*Power(xij,2.0) - 
              3836450520.0*Power(rij,3.0)*Power(xij,3.0) - 
              337429260.0*Power(rij,4.0)*Power(xij,4.0) + 
              24444504.0*Power(rij,5.0)*Power(xij,5.0) + 
              10247328.0*Power(rij,6.0)*Power(xij,6.0) + 
              1284588.0*Power(rij,7.0)*Power(xij,7.0) + 
              86157.0*Power(rij,8.0)*Power(xij,8.0) + 
              3026.0*Power(rij,9.0)*Power(xij,9.0) + 40.0*Power(rij,10.0)*Power(xij,10.0)) + 
           12.0*Power(xii,16.0)*Power(xij,8.0)*
            (-2811375.0 - 5622750.0*rij*xij - 5622750.0*Power(rij,2.0)*Power(xij,2.0) - 
              3298050.0*Power(rij,3.0)*Power(xij,3.0) - 
              4351725.0*Power(rij,4.0)*Power(xij,4.0) + 
              2385432.0*Power(rij,5.0)*Power(xij,5.0) - 
              1111761.0*Power(rij,6.0)*Power(xij,6.0) - 
              242604.0*Power(rij,7.0)*Power(xij,7.0) + 
              1950.0*Power(rij,8.0)*Power(xij,8.0) + 2072.0*Power(rij,9.0)*Power(xij,9.0) + 
              73.0*Power(rij,10.0)*Power(xij,10.0)) - 
           3.0*Power(xii,8.0)*Power(xij,16.0)*
            (-9364244085.0 - 18728488170.0*rij*xij + 
              11763172890.0*Power(rij,2.0)*Power(xij,2.0) + 
              4695905340.0*Power(rij,3.0)*Power(xij,3.0) - 
              11704770.0*Power(rij,4.0)*Power(xij,4.0) - 
              211923684.0*Power(rij,5.0)*Power(xij,5.0) - 
              37532628.0*Power(rij,6.0)*Power(xij,6.0) - 
              2522664.0*Power(rij,7.0)*Power(xij,7.0) - 
              5874.0*Power(rij,8.0)*Power(xij,8.0) + 7772.0*Power(rij,9.0)*Power(xij,9.0) + 
              292.0*Power(rij,10.0)*Power(xij,10.0)) + 
           Power(xii,6.0)*Power(xij,18.0)*
            (57304872765.0 + 114609745530.0*rij*xij + 
              20096073270.0*Power(rij,2.0)*Power(xij,2.0) - 
              7496628300.0*Power(rij,3.0)*Power(xij,3.0) - 
              3291811110.0*Power(rij,4.0)*Power(xij,4.0) - 
              499851324.0*Power(rij,5.0)*Power(xij,5.0) - 
              28415268.0*Power(rij,6.0)*Power(xij,6.0) + 
              1457928.0*Power(rij,7.0)*Power(xij,7.0) + 
              330210.0*Power(rij,8.0)*Power(xij,8.0) + 
              19796.0*Power(rij,9.0)*Power(xij,9.0) + 424.0*Power(rij,10.0)*Power(xij,10.0)))\
    )/(14175.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),17.0))
    ;
  }
  return S;
}

static double DSlater_4S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-729870602452992000.0 + 729870602452992000.0*Power(E,2.0*rij*xii) - 
        1459741204905984000.0*rij*xii - 
        1459741204905984000.0*Power(rij,2.0)*Power(xii,2.0) - 
        971930171227640625.0*Power(rij,3.0)*Power(xii,3.0) - 
        484119137549297250.0*Power(rij,4.0)*Power(xii,4.0) - 
        192174113906775000.0*Power(rij,5.0)*Power(xii,5.0) - 
        63242978838039000.0*Power(rij,6.0)*Power(xii,6.0) - 
        17722869187923900.0*Power(rij,7.0)*Power(xii,7.0) - 
        4311139542910200.0*Power(rij,8.0)*Power(xii,8.0) - 
        923450838710400.0*Power(rij,9.0)*Power(xii,9.0) - 
        176129454140640.0*Power(rij,10.0)*Power(xii,10.0) - 
        30179820041280.0*Power(rij,11.0)*Power(xii,11.0) - 
        4679384411520.0*Power(rij,12.0)*Power(xii,12.0) - 
        660128071680.0*Power(rij,13.0)*Power(xii,13.0) - 
        85016494080.0*Power(rij,14.0)*Power(xii,14.0) - 
        10001940480.0*Power(rij,15.0)*Power(xii,15.0) - 
        1071636480.0*Power(rij,16.0)*Power(xii,16.0) - 
        103661568.0*Power(rij,17.0)*Power(xii,17.0) - 
        8871936.0*Power(rij,18.0)*Power(xii,18.0) - 
        638976.0*Power(rij,19.0)*Power(xii,19.0) - 32768.0*Power(rij,20.0)*Power(xii,20.0))/
      (7.29870602452992e17*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (467775.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),19.0) + 
        495.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-318.0*Power(rij,7.0)*Power(xii,31.0) - 6.0*Power(rij,8.0)*Power(xii,32.0) + 
           945.0*Power(xij,24.0) + 1890.0*rij*xii*Power(xij,24.0) + 
           630.0*rij*Power(xii,3.0)*Power(xij,22.0)*
            (-57.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           945.0*Power(xii,2.0)*Power(xij,22.0)*(-19.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           4.0*Power(rij,6.0)*Power(xii,30.0)*(1932.0 + 17.0*Power(rij,2.0)*Power(xij,2.0)) - 
           18.0*Power(rij,5.0)*Power(xii,29.0)*
            (6188.0 + 271.0*Power(rij,2.0)*Power(xij,2.0)) + 
           126.0*rij*Power(xii,5.0)*Power(xij,20.0)*
            (2565.0 - 190.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           315.0*Power(xii,4.0)*Power(xij,20.0)*
            (513.0 - 114.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 4.0*Power(rij,4.0)*Power(xii,28.0)*
            (-257040.0 - 37821.0*Power(rij,2.0)*Power(xij,2.0) + 
              62.0*Power(rij,4.0)*Power(xij,4.0)) + 
           6.0*Power(rij,3.0)*Power(xii,27.0)*
            (-1017450.0 - 446082.0*Power(rij,2.0)*Power(xij,2.0) + 
              1129.0*Power(rij,4.0)*Power(xij,4.0)) - 
           6.0*rij*Power(xii,25.0)*(6715170.0 + 34042680.0*Power(rij,2.0)*Power(xij,2.0) + 
              1106742.0*Power(rij,4.0)*Power(xij,4.0) - 
              7261.0*Power(rij,6.0)*Power(xij,6.0)) - 
           6.0*rij*Power(xii,13.0)*Power(xij,12.0)*
            (622475910.0 + 166454820.0*Power(rij,2.0)*Power(xij,2.0) + 
              1261904.0*Power(rij,4.0)*Power(xij,4.0) - 
              4639.0*Power(rij,6.0)*Power(xij,6.0)) - 
           6.0*Power(rij,2.0)*Power(xii,26.0)*
            (3662820.0 + 4913685.0*Power(rij,2.0)*Power(xij,2.0) + 
              15134.0*Power(rij,4.0)*Power(xij,4.0) - 46.0*Power(rij,6.0)*Power(xij,6.0)) + 
           6.0*rij*Power(xii,7.0)*Power(xij,18.0)*
            (-305235.0 + 35910.0*Power(rij,2.0)*Power(xij,2.0) - 
              798.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           21.0*Power(xii,6.0)*Power(xij,18.0)*
            (-43605.0 + 15390.0*Power(rij,2.0)*Power(xij,2.0) - 
              570.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           114.0*rij*Power(xii,15.0)*Power(xij,10.0)*
            (-241392690.0 - 14001540.0*Power(rij,2.0)*Power(xij,2.0) + 
              240380.0*Power(rij,4.0)*Power(xij,4.0) + 101.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 6.0*rij*Power(xii,9.0)*Power(xij,16.0)*
            (-1220940.0 + 342720.0*Power(rij,2.0)*Power(xij,2.0) + 
              462.0*Power(rij,4.0)*Power(xij,4.0) + 167.0*Power(rij,6.0)*Power(xij,6.0)) - 
           342.0*rij*Power(xii,17.0)*Power(xij,8.0)*
            (183130605.0 - 6697320.0*Power(rij,2.0)*Power(xij,2.0) - 
              26180.0*Power(rij,4.0)*Power(xij,4.0) + 374.0*Power(rij,6.0)*Power(xij,6.0)) \
    - 18.0*rij*Power(xii,11.0)*Power(xij,14.0)*
            (7393470.0 + 4725490.0*Power(rij,2.0)*Power(xij,2.0) + 
              169358.0*Power(rij,4.0)*Power(xij,4.0) + 409.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 42.0*rij*Power(xii,19.0)*Power(xij,6.0)*
            (-1241513955.0 + 70422930.0*Power(rij,2.0)*Power(xij,2.0) - 
              1293292.0*Power(rij,4.0)*Power(xij,4.0) + 
              3230.0*Power(rij,6.0)*Power(xij,6.0)) + 
           6.0*rij*Power(xii,21.0)*Power(xij,4.0)*
            (-2722031235.0 - 183644790.0*Power(rij,2.0)*Power(xij,2.0) + 
              2673034.0*Power(rij,4.0)*Power(xij,4.0) + 
              3230.0*Power(rij,6.0)*Power(xij,6.0)) - 
           18.0*rij*Power(xii,23.0)*Power(xij,2.0)*
            (95263245.0 + 69677230.0*Power(rij,2.0)*Power(xij,2.0) - 
              1221038.0*Power(rij,4.0)*Power(xij,4.0) + 
              5738.0*Power(rij,6.0)*Power(xij,6.0)) - 
           21.0*Power(xii,20.0)*Power(xij,4.0)*
            (388861605.0 + 679369590.0*Power(rij,2.0)*Power(xij,2.0) - 
              21408630.0*Power(rij,4.0)*Power(xij,4.0) + 
              121448.0*Power(rij,6.0)*Power(xij,6.0)) + 
           6.0*Power(xii,8.0)*Power(xij,16.0)*
            (610470.0 - 305235.0*Power(rij,2.0)*Power(xij,2.0) + 
              17955.0*Power(rij,4.0)*Power(xij,4.0) - 
              266.0*Power(rij,6.0)*Power(xij,6.0) + Power(rij,8.0)*Power(xij,8.0)) + 
           2.0*Power(xii,10.0)*Power(xij,14.0)*
            (-33270615.0 + 16889670.0*Power(rij,2.0)*Power(xij,2.0) + 
              1365525.0*Power(rij,4.0)*Power(xij,4.0) + 
              37758.0*Power(rij,6.0)*Power(xij,6.0) + 34.0*Power(rij,8.0)*Power(xij,8.0)) - 
           6.0*Power(xii,14.0)*Power(xij,10.0)*
            (2293230555.0 - 1340795610.0*Power(rij,2.0)*Power(xij,2.0) - 
              20823810.0*Power(rij,4.0)*Power(xij,4.0) + 
              201628.0*Power(rij,6.0)*Power(xij,6.0) + 46.0*Power(rij,8.0)*Power(xij,8.0)) \
    + 19.0*Power(xii,16.0)*Power(xij,8.0)*
            (-1648175445.0 + 756275940.0*Power(rij,2.0)*Power(xij,2.0) - 
              18485460.0*Power(rij,4.0)*Power(xij,4.0) - 
              14280.0*Power(rij,6.0)*Power(xij,6.0) + 106.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 2.0*Power(xii,12.0)*Power(xij,12.0)*
            (933713865.0 - 515440170.0*Power(rij,2.0)*Power(xij,2.0) - 
              35610750.0*Power(rij,4.0)*Power(xij,4.0) - 
              158046.0*Power(rij,6.0)*Power(xij,6.0) + 124.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 2.0*Power(xii,24.0)*(10072755.0 + 420272685.0*Power(rij,2.0)*Power(xij,2.0) + 
              63685755.0*Power(rij,4.0)*Power(xij,4.0) - 
              770154.0*Power(rij,6.0)*Power(xij,6.0) + 1007.0*Power(rij,8.0)*Power(xij,8.0)\
    ) - Power(xii,18.0)*Power(xij,6.0)*
            (26071793055.0 + 1689529590.0*Power(rij,2.0)*Power(xij,2.0) + 
              229129740.0*Power(rij,4.0)*Power(xij,4.0) - 
              3527160.0*Power(rij,6.0)*Power(xij,6.0) + 
              2584.0*Power(rij,8.0)*Power(xij,8.0)) + 
           Power(xii,22.0)*Power(xij,2.0)*
            (-857369205.0 - 6658320690.0*Power(rij,2.0)*Power(xij,2.0) + 
              89662230.0*Power(rij,4.0)*Power(xij,4.0) - 
              1176252.0*Power(rij,6.0)*Power(xij,6.0) + 2584.0*Power(rij,8.0)*Power(xij,8.0)\
    )) + Power(E,2.0*rij*xii)*Power(xii,10.0)*
         (-21318.0*Power(xii,14.0)*Power(xij,14.0)*
            (-1573425.0 - 3146850.0*rij*xij + 17151750.0*Power(rij,2.0)*Power(xij,2.0) - 
              36684900.0*Power(rij,3.0)*Power(xij,3.0) + 
              14486850.0*Power(rij,4.0)*Power(xij,4.0) + 
              1682100.0*Power(rij,5.0)*Power(xij,5.0) - 
              1152900.0*Power(rij,6.0)*Power(xij,6.0) - 
              221940.0*Power(rij,7.0)*Power(xij,7.0) - 
              1620.0*Power(rij,8.0)*Power(xij,8.0) + 2760.0*Power(rij,9.0)*Power(xij,9.0) + 
              264.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,11.0)*Power(xij,11.0)) + 
           9.0*Power(xii,22.0)*Power(xij,6.0)*
            (50363775.0 + 100727550.0*rij*xij + 
              100727550.0*Power(rij,2.0)*Power(xij,2.0) + 
              67151700.0*Power(rij,3.0)*Power(xij,3.0) + 
              33575850.0*Power(rij,4.0)*Power(xij,4.0) + 
              13430340.0*Power(rij,5.0)*Power(xij,5.0) + 
              4476780.0*Power(rij,6.0)*Power(xij,6.0) + 
              1562220.0*Power(rij,7.0)*Power(xij,7.0) - 
              104940.0*Power(rij,8.0)*Power(xij,8.0) + 
              199320.0*Power(rij,9.0)*Power(xij,9.0) + 
              28248.0*Power(rij,10.0)*Power(xij,10.0) + 
              56.0*Power(rij,11.0)*Power(xij,11.0) - 64.0*Power(rij,12.0)*Power(xij,12.0)) + 
           2.0*Power(xii,2.0)*Power(xij,26.0)*
            (468534198825.0 + 937068397650.0*rij*xij + 
              648987604650.0*Power(rij,2.0)*Power(xij,2.0) + 
              248897776050.0*Power(rij,3.0)*Power(xij,3.0) + 
              62249625900.0*Power(rij,4.0)*Power(xij,4.0) + 
              10932296760.0*Power(rij,5.0)*Power(xij,5.0) + 
              1391848920.0*Power(rij,6.0)*Power(xij,6.0) + 
              128752470.0*Power(rij,7.0)*Power(xij,7.0) + 
              8381835.0*Power(rij,8.0)*Power(xij,8.0) + 
              346830.0*Power(rij,9.0)*Power(xij,9.0) + 
              6006.0*Power(rij,10.0)*Power(xij,10.0) - 
              158.0*Power(rij,11.0)*Power(xij,11.0) - 8.0*Power(rij,12.0)*Power(xij,12.0)) + 
           2.0*Power(xij,28.0)*(13749310575.0 + 27498621150.0*rij*xij + 
              20951330400.0*Power(rij,2.0)*Power(xij,2.0) + 
              9166207050.0*Power(rij,3.0)*Power(xij,3.0) + 
              2687835150.0*Power(rij,4.0)*Power(xij,4.0) + 
              569729160.0*Power(rij,5.0)*Power(xij,5.0) + 
              90810720.0*Power(rij,6.0)*Power(xij,6.0) + 
              11081070.0*Power(rij,7.0)*Power(xij,7.0) + 
              1036035.0*Power(rij,8.0)*Power(xij,8.0) + 
              72930.0*Power(rij,9.0)*Power(xij,9.0) + 
              3696.0*Power(rij,10.0)*Power(xij,10.0) + 
              122.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) + 
           1254.0*Power(xii,16.0)*Power(xij,12.0)*
            (-10120950.0 - 20241900.0*rij*xij - 
              27471150.0*Power(rij,2.0)*Power(xij,2.0) + 
              28435050.0*Power(rij,3.0)*Power(xij,3.0) - 
              71328600.0*Power(rij,4.0)*Power(xij,4.0) + 
              25843860.0*Power(rij,5.0)*Power(xij,5.0) + 
              865620.0*Power(rij,6.0)*Power(xij,6.0) - 
              2061990.0*Power(rij,7.0)*Power(xij,7.0) - 
              276885.0*Power(rij,8.0)*Power(xij,8.0) + 
              5010.0*Power(rij,9.0)*Power(xij,9.0) + 
              2856.0*Power(rij,10.0)*Power(xij,10.0) + 
              162.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) - 
           627.0*Power(xii,12.0)*Power(xij,16.0)*
            (-620482275.0 - 1240964550.0*rij*xij + 
              6396298650.0*Power(rij,2.0)*Power(xij,2.0) - 
              3446622900.0*Power(rij,3.0)*Power(xij,3.0) - 
              476705250.0*Power(rij,4.0)*Power(xij,4.0) + 
              290390940.0*Power(rij,5.0)*Power(xij,5.0) + 
              67514580.0*Power(rij,6.0)*Power(xij,6.0) + 
              1373760.0*Power(rij,7.0)*Power(xij,7.0) - 
              1053270.0*Power(rij,8.0)*Power(xij,8.0) - 
              139980.0*Power(rij,9.0)*Power(xij,9.0) - 
              6828.0*Power(rij,10.0)*Power(xij,10.0) - 
              56.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           Power(xii,28.0)*(467775.0 + 935550.0*rij*xij + 
              935550.0*Power(rij,2.0)*Power(xij,2.0) + 
              623700.0*Power(rij,3.0)*Power(xij,3.0) + 
              311850.0*Power(rij,4.0)*Power(xij,4.0) + 
              124740.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              11880.0*Power(rij,7.0)*Power(xij,7.0) + 
              2970.0*Power(rij,8.0)*Power(xij,8.0) + 660.0*Power(rij,9.0)*Power(xij,9.0) + 
              132.0*Power(rij,10.0)*Power(xij,10.0) + 
              24.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           6.0*Power(xii,4.0)*Power(xij,24.0)*
            (-1131220750275.0 - 2262441500550.0*rij*xij - 
              1290622098150.0*Power(rij,2.0)*Power(xij,2.0) - 
              366043916700.0*Power(rij,3.0)*Power(xij,3.0) - 
              57112520850.0*Power(rij,4.0)*Power(xij,4.0) - 
              3706773840.0*Power(rij,5.0)*Power(xij,5.0) + 
              383949720.0*Power(rij,6.0)*Power(xij,6.0) + 
              124315290.0*Power(rij,7.0)*Power(xij,7.0) + 
              15397965.0*Power(rij,8.0)*Power(xij,8.0) + 
              1138170.0*Power(rij,9.0)*Power(xij,9.0) + 
              51414.0*Power(rij,10.0)*Power(xij,10.0) + 
              1248.0*Power(rij,11.0)*Power(xij,11.0) + 10.0*Power(rij,12.0)*Power(xij,12.0)) \
    + Power(xii,26.0)*Power(xij,2.0)*(8887725.0 + 17775450.0*rij*xij + 
              17775450.0*Power(rij,2.0)*Power(xij,2.0) + 
              11850300.0*Power(rij,3.0)*Power(xij,3.0) + 
              5925150.0*Power(rij,4.0)*Power(xij,4.0) + 
              2370060.0*Power(rij,5.0)*Power(xij,5.0) + 
              790020.0*Power(rij,6.0)*Power(xij,6.0) + 
              225720.0*Power(rij,7.0)*Power(xij,7.0) + 
              56430.0*Power(rij,8.0)*Power(xij,8.0) + 
              12540.0*Power(rij,9.0)*Power(xij,9.0) + 
              2508.0*Power(rij,10.0)*Power(xij,10.0) + 
              756.0*Power(rij,11.0)*Power(xij,11.0) + 16.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 3.0*Power(xii,24.0)*Power(xij,4.0)*
            (-26663175.0 - 53326350.0*rij*xij - 
              53326350.0*Power(rij,2.0)*Power(xij,2.0) - 
              35550900.0*Power(rij,3.0)*Power(xij,3.0) - 
              17775450.0*Power(rij,4.0)*Power(xij,4.0) - 
              7110180.0*Power(rij,5.0)*Power(xij,5.0) - 
              2370060.0*Power(rij,6.0)*Power(xij,6.0) - 
              677160.0*Power(rij,7.0)*Power(xij,7.0) - 
              169290.0*Power(rij,8.0)*Power(xij,8.0) - 
              23100.0*Power(rij,9.0)*Power(xij,9.0) - 
              17688.0*Power(rij,10.0)*Power(xij,10.0) - 
              1156.0*Power(rij,11.0)*Power(xij,11.0) + 20.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 9.0*Power(xii,6.0)*Power(xij,22.0)*
            (1726271714475.0 + 3452543428950.0*rij*xij + 
              1148702905350.0*Power(rij,2.0)*Power(xij,2.0) + 
              9136581300.0*Power(rij,3.0)*Power(xij,3.0) - 
              73193377950.0*Power(rij,4.0)*Power(xij,4.0) - 
              19535877900.0*Power(rij,5.0)*Power(xij,5.0) - 
              2453012100.0*Power(rij,6.0)*Power(xij,6.0) - 
              132309540.0*Power(rij,7.0)*Power(xij,7.0) + 
              5655540.0*Power(rij,8.0)*Power(xij,8.0) + 
              1529880.0*Power(rij,9.0)*Power(xij,9.0) + 
              116952.0*Power(rij,10.0)*Power(xij,10.0) + 
              4344.0*Power(rij,11.0)*Power(xij,11.0) + 64.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 22.0*Power(xii,20.0)*Power(xij,8.0)*
            (-82413450.0 - 164826900.0*rij*xij - 
              164826900.0*Power(rij,2.0)*Power(xij,2.0) - 
              109884600.0*Power(rij,3.0)*Power(xij,3.0) - 
              54942300.0*Power(rij,4.0)*Power(xij,4.0) - 
              19274220.0*Power(rij,5.0)*Power(xij,5.0) - 
              15433740.0*Power(rij,6.0)*Power(xij,6.0) + 
              4200390.0*Power(rij,7.0)*Power(xij,7.0) - 
              1404855.0*Power(rij,8.0)*Power(xij,8.0) - 
              416910.0*Power(rij,9.0)*Power(xij,9.0) - 
              14874.0*Power(rij,10.0)*Power(xij,10.0) + 
              1672.0*Power(rij,11.0)*Power(xij,11.0) + 82.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 22.0*Power(xii,18.0)*Power(xij,10.0)*
            (-247240350.0 - 494480700.0*rij*xij - 
              494480700.0*Power(rij,2.0)*Power(xij,2.0) - 
              348449850.0*Power(rij,3.0)*Power(xij,3.0) - 
              33254550.0*Power(rij,4.0)*Power(xij,4.0) - 
              300279420.0*Power(rij,5.0)*Power(xij,5.0) + 
              104656860.0*Power(rij,6.0)*Power(xij,6.0) - 
              10052370.0*Power(rij,7.0)*Power(xij,7.0) - 
              8727255.0*Power(rij,8.0)*Power(xij,8.0) - 
              698070.0*Power(rij,9.0)*Power(xij,9.0) + 
              37218.0*Power(rij,10.0)*Power(xij,10.0) + 
              5686.0*Power(rij,11.0)*Power(xij,11.0) + 136.0*Power(rij,12.0)*Power(xij,12.0)\
    ) - 11.0*Power(xii,8.0)*Power(xij,20.0)*
            (-1169302313025.0 - 2338604626050.0*rij*xij + 
              403746553350.0*Power(rij,2.0)*Power(xij,2.0) + 
              452949594300.0*Power(rij,3.0)*Power(xij,3.0) + 
              74923465050.0*Power(rij,4.0)*Power(xij,4.0) - 
              5463441900.0*Power(rij,5.0)*Power(xij,5.0) - 
              3490244100.0*Power(rij,6.0)*Power(xij,6.0) - 
              521907840.0*Power(rij,7.0)*Power(xij,7.0) - 
              35495010.0*Power(rij,8.0)*Power(xij,8.0) - 
              393420.0*Power(rij,9.0)*Power(xij,9.0) + 
              112152.0*Power(rij,10.0)*Power(xij,10.0) + 
              7644.0*Power(rij,11.0)*Power(xij,11.0) + 164.0*Power(rij,12.0)*Power(xij,12.0)\
    ) + 11.0*Power(xii,10.0)*Power(xij,18.0)*
            (371402591175.0 + 742805182350.0*rij*xij - 
              961651930050.0*Power(rij,2.0)*Power(xij,2.0) - 
              72424667700.0*Power(rij,3.0)*Power(xij,3.0) + 
              100904879250.0*Power(rij,4.0)*Power(xij,4.0) + 
              22702033620.0*Power(rij,5.0)*Power(xij,5.0) + 
              209045340.0*Power(rij,6.0)*Power(xij,6.0) - 
              517637520.0*Power(rij,7.0)*Power(xij,7.0) - 
              79564410.0*Power(rij,8.0)*Power(xij,8.0) - 
              5021940.0*Power(rij,9.0)*Power(xij,9.0) - 
              74724.0*Power(rij,10.0)*Power(xij,10.0) + 
              6852.0*Power(rij,11.0)*Power(xij,11.0) + 272.0*Power(rij,12.0)*Power(xij,12.0))\
    ))/(467775.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),19.0))
    ;
  }
  return S;
}

double DSlater_4S_1S(double rij,double xii,double xij)
{
  return DSlater_1S_4S(rij,xij,xii);
}

double DSlater_4S_2S(double rij,double xii,double xij)
{
  return DSlater_2S_4S(rij,xij,xii);
}

double DSlater_4S_3S(double rij,double xii,double xij)
{
  return DSlater_3S_4S(rij,xij,xii);
}

static double DSlater_5S_5S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-3041127510220800000.0 + 3041127510220800000.0*Power(E,2.0*rij*xii) - 
        6082255020441600000.0*rij*xii - 
        6082255020441600000.0*Power(rij,2.0)*Power(xii,2.0) - 
        4047316151142639375.0*Power(rij,3.0)*Power(xii,3.0) - 
        2012377281843678750.0*Power(rij,4.0)*Power(xii,4.0) - 
        796103231382459000.0*Power(rij,5.0)*Power(xii,5.0) - 
        260648980404813000.0*Power(rij,6.0)*Power(xii,6.0) - 
        72570149721669600.0*Power(rij,7.0)*Power(xii,7.0) - 
        17529098189803200.0*Power(rij,8.0)*Power(xii,8.0) - 
        3730342475059200.0*Power(rij,9.0)*Power(xii,9.0) - 
        707903551555200.0*Power(rij,10.0)*Power(xii,10.0) - 
        120923460403200.0*Power(rij,11.0)*Power(xii,11.0) - 
        18723632578560.0*Power(rij,12.0)*Power(xii,12.0) - 
        2640512286720.0*Power(rij,13.0)*Power(xii,13.0) - 
        340065976320.0*Power(rij,14.0)*Power(xii,14.0) - 
        40007761920.0*Power(rij,15.0)*Power(xii,15.0) - 
        4286545920.0*Power(rij,16.0)*Power(xii,16.0) - 
        414646272.0*Power(rij,17.0)*Power(xii,17.0) - 
        35487744.0*Power(rij,18.0)*Power(xii,18.0) - 
        2555904.0*Power(rij,19.0)*Power(xii,19.0) - 131072.0*Power(rij,20.0)*Power(xii,20.0))/
      (3.0411275102208e18*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (70875.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),19.0) + 
        Power(E,2.0*rij*xii)*Power(xii,12.0)*
         (Power(xii,8.0)*Power(xij,18.0)*
            (3218321469825.0 + 6436642939650.0*rij*xij - 
              289335546990.0*Power(rij,2.0)*Power(xij,2.0) - 
              672080979780.0*Power(rij,3.0)*Power(xij,3.0) - 
              116652948930.0*Power(rij,4.0)*Power(xij,4.0) - 
              2285769780.0*Power(rij,5.0)*Power(xij,5.0) + 
              1432351620.0*Power(rij,6.0)*Power(xij,6.0) + 
              183837384.0*Power(rij,7.0)*Power(xij,7.0) + 
              9315018.0*Power(rij,8.0)*Power(xij,8.0) + 
              153748.0*Power(rij,9.0)*Power(xij,9.0) - 
              1636.0*Power(rij,10.0)*Power(xij,10.0)) + 
           10.0*Power(xij,26.0)*(352546425.0 + 705092850.0*rij*xij + 
              467009550.0*Power(rij,2.0)*Power(xij,2.0) + 
              168489720.0*Power(rij,3.0)*Power(xij,3.0) + 
              39134340.0*Power(rij,4.0)*Power(xij,4.0) + 
              6297480.0*Power(rij,5.0)*Power(xij,5.0) + 
              723240.0*Power(rij,6.0)*Power(xij,6.0) + 
              59220.0*Power(rij,7.0)*Power(xij,7.0) + 
              3339.0*Power(rij,8.0)*Power(xij,8.0) + 118.0*Power(rij,9.0)*Power(xij,9.0) + 
              2.0*Power(rij,10.0)*Power(xij,10.0)) + 
           30.0*Power(xii,2.0)*Power(xij,24.0)*
            (4562958015.0 + 9125916030.0*rij*xij + 
              5463096030.0*Power(rij,2.0)*Power(xij,2.0) + 
              1726409160.0*Power(rij,3.0)*Power(xij,3.0) + 
              343084140.0*Power(rij,4.0)*Power(xij,4.0) + 
              46070136.0*Power(rij,5.0)*Power(xij,5.0) + 
              4278792.0*Power(rij,6.0)*Power(xij,6.0) + 
              271212.0*Power(rij,7.0)*Power(xij,7.0) + 
              11061.0*Power(rij,8.0)*Power(xij,8.0) + 
              250.0*Power(rij,9.0)*Power(xij,9.0) + 2.0*Power(rij,10.0)*Power(xij,10.0)) - 
           15.0*Power(xii,24.0)*Power(xij,2.0)*
            (-89775.0 - 179550.0*rij*xij - 179550.0*Power(rij,2.0)*Power(xij,2.0) - 
              119700.0*Power(rij,3.0)*Power(xij,3.0) - 
              59850.0*Power(rij,4.0)*Power(xij,4.0) - 
              23940.0*Power(rij,5.0)*Power(xij,5.0) - 
              7980.0*Power(rij,6.0)*Power(xij,6.0) - 
              2280.0*Power(rij,7.0)*Power(xij,7.0) - 570.0*Power(rij,8.0)*Power(xij,8.0) - 
              244.0*Power(rij,9.0)*Power(xij,9.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           5.0*Power(xii,26.0)*(14175.0 + 28350.0*rij*xij + 
              28350.0*Power(rij,2.0)*Power(xij,2.0) + 
              18900.0*Power(rij,3.0)*Power(xij,3.0) + 
              9450.0*Power(rij,4.0)*Power(xij,4.0) + 
              3780.0*Power(rij,5.0)*Power(xij,5.0) + 
              1260.0*Power(rij,6.0)*Power(xij,6.0) + 360.0*Power(rij,7.0)*Power(xij,7.0) + 
              90.0*Power(rij,8.0)*Power(xij,8.0) + 20.0*Power(rij,9.0)*Power(xij,9.0) + 
              4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           1938.0*Power(xii,14.0)*Power(xij,12.0)*
            (-826875.0 - 1653750.0*rij*xij + 55046250.0*Power(rij,2.0)*Power(xij,2.0) - 
              71486100.0*Power(rij,3.0)*Power(xij,3.0) + 
              20956950.0*Power(rij,4.0)*Power(xij,4.0) + 
              4028220.0*Power(rij,5.0)*Power(xij,5.0) - 
              471660.0*Power(rij,6.0)*Power(xij,6.0) - 
              108192.0*Power(rij,7.0)*Power(xij,7.0) - 
              4284.0*Power(rij,8.0)*Power(xij,8.0) + 136.0*Power(rij,9.0)*Power(xij,9.0) + 
              8.0*Power(rij,10.0)*Power(xij,10.0)) + 
           1938.0*Power(xii,12.0)*Power(xij,14.0)*
            (72476775.0 + 144953550.0*rij*xij - 
              458924130.0*Power(rij,2.0)*Power(xij,2.0) + 
              175365540.0*Power(rij,3.0)*Power(xij,3.0) + 
              35155890.0*Power(rij,4.0)*Power(xij,4.0) - 
              5303340.0*Power(rij,5.0)*Power(xij,5.0) - 
              1428420.0*Power(rij,6.0)*Power(xij,6.0) - 
              75552.0*Power(rij,7.0)*Power(xij,7.0) + 
              3036.0*Power(rij,8.0)*Power(xij,8.0) + 376.0*Power(rij,9.0)*Power(xij,9.0) + 
              8.0*Power(rij,10.0)*Power(xij,10.0)) + 
           342.0*Power(xii,16.0)*Power(xij,10.0)*
            (2409750.0 + 4819500.0*rij*xij - 2142000.0*Power(rij,2.0)*Power(xij,2.0) + 
              35235900.0*Power(rij,3.0)*Power(xij,3.0) - 
              35289450.0*Power(rij,4.0)*Power(xij,4.0) + 
              11000220.0*Power(rij,5.0)*Power(xij,5.0) + 
              1519140.0*Power(rij,6.0)*Power(xij,6.0) - 
              194172.0*Power(rij,7.0)*Power(xij,7.0) - 
              29069.0*Power(rij,8.0)*Power(xij,8.0) - 
              634.0*Power(rij,9.0)*Power(xij,9.0) + 18.0*Power(rij,10.0)*Power(xij,10.0)) - 
           171.0*Power(xii,10.0)*Power(xij,16.0)*
            (-6768406575.0 - 13536813150.0*rij*xij + 
              12122613090.0*Power(rij,2.0)*Power(xij,2.0) + 
              1678134780.0*Power(rij,3.0)*Power(xij,3.0) - 
              578956770.0*Power(rij,4.0)*Power(xij,4.0) - 
              138373620.0*Power(rij,5.0)*Power(xij,5.0) - 
              7287420.0*Power(rij,6.0)*Power(xij,6.0) + 
              614856.0*Power(rij,7.0)*Power(xij,7.0) + 
              89482.0*Power(rij,8.0)*Power(xij,8.0) + 
              3572.0*Power(rij,9.0)*Power(xij,9.0) + 36.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 9.0*Power(xii,22.0)*Power(xij,4.0)*
            (-1346625.0 - 2693250.0*rij*xij - 2693250.0*Power(rij,2.0)*Power(xij,2.0) - 
              1795500.0*Power(rij,3.0)*Power(xij,3.0) - 
              897750.0*Power(rij,4.0)*Power(xij,4.0) - 
              359100.0*Power(rij,5.0)*Power(xij,5.0) - 
              119700.0*Power(rij,6.0)*Power(xij,6.0) - 
              10176.0*Power(rij,7.0)*Power(xij,7.0) - 
              30572.0*Power(rij,8.0)*Power(xij,8.0) + 
              168.0*Power(rij,9.0)*Power(xij,9.0) + 104.0*Power(rij,10.0)*Power(xij,10.0)) \
    - 9.0*Power(xii,4.0)*Power(xij,22.0)*
            (-129194933175.0 - 258389866350.0*rij*xij - 
              128354872590.0*Power(rij,2.0)*Power(xij,2.0) - 
              30914128980.0*Power(rij,3.0)*Power(xij,3.0) - 
              4146276330.0*Power(rij,4.0)*Power(xij,4.0) - 
              281941380.0*Power(rij,5.0)*Power(xij,5.0) + 
              311220.0*Power(rij,6.0)*Power(xij,6.0) + 
              1834944.0*Power(rij,7.0)*Power(xij,7.0) + 
              162188.0*Power(rij,8.0)*Power(xij,8.0) + 
              6488.0*Power(rij,9.0)*Power(xij,9.0) + 104.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 9.0*Power(xii,6.0)*Power(xij,20.0)*
            (356863797675.0 + 713727595350.0*rij*xij + 
              226198496790.0*Power(rij,2.0)*Power(xij,2.0) + 
              15231757380.0*Power(rij,3.0)*Power(xij,3.0) - 
              5016397470.0*Power(rij,4.0)*Power(xij,4.0) - 
              1294411020.0*Power(rij,5.0)*Power(xij,5.0) - 
              134742020.0*Power(rij,6.0)*Power(xij,6.0) - 
              6600064.0*Power(rij,7.0)*Power(xij,7.0) - 
              49228.0*Power(rij,8.0)*Power(xij,8.0) + 
              9192.0*Power(rij,9.0)*Power(xij,9.0) + 296.0*Power(rij,10.0)*Power(xij,10.0)) \
    - 9.0*Power(xii,20.0)*Power(xij,6.0)*
            (-7630875.0 - 15261750.0*rij*xij - 
              15261750.0*Power(rij,2.0)*Power(xij,2.0) - 
              10174500.0*Power(rij,3.0)*Power(xij,3.0) - 
              5087250.0*Power(rij,4.0)*Power(xij,4.0) - 
              2995860.0*Power(rij,5.0)*Power(xij,5.0) + 
              1403780.0*Power(rij,6.0)*Power(xij,6.0) - 
              1201664.0*Power(rij,7.0)*Power(xij,7.0) - 
              32148.0*Power(rij,8.0)*Power(xij,8.0) + 
              9752.0*Power(rij,9.0)*Power(xij,9.0) + 296.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 2.0*Power(xii,18.0)*Power(xij,8.0)*
            (-137355750.0 - 274711500.0*rij*xij - 
              274711500.0*Power(rij,2.0)*Power(xij,2.0) - 
              137195100.0*Power(rij,3.0)*Power(xij,3.0) - 
              344272950.0*Power(rij,4.0)*Power(xij,4.0) + 
              294722820.0*Power(rij,5.0)*Power(xij,5.0) - 
              125182260.0*Power(rij,6.0)*Power(xij,6.0) - 
              9557892.0*Power(rij,7.0)*Power(xij,7.0) + 
              1628541.0*Power(rij,8.0)*Power(xij,8.0) + 
              129226.0*Power(rij,9.0)*Power(xij,9.0) + 818.0*Power(rij,10.0)*Power(xij,10.0)\
    )) + Power(E,2.0*rij*xij)*Power(xij,12.0)*
         (-1180.0*Power(rij,9.0)*Power(xii,35.0) - 20.0*Power(rij,10.0)*Power(xii,36.0) + 
           70875.0*Power(xij,26.0) + 141750.0*rij*xii*Power(xij,26.0) + 
           47250.0*rij*Power(xii,3.0)*Power(xij,24.0)*
            (-57.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           70875.0*Power(xii,2.0)*Power(xij,24.0)*
            (-19.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           30.0*Power(rij,8.0)*Power(xii,34.0)*(1113.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           300.0*Power(rij,7.0)*Power(xii,33.0)*
            (1974.0 + 25.0*Power(rij,2.0)*Power(xij,2.0)) + 
           9450.0*rij*Power(xii,5.0)*Power(xij,22.0)*
            (2565.0 - 190.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 23625.0*Power(xii,4.0)*Power(xij,22.0)*
            (513.0 - 114.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 18.0*Power(rij,6.0)*Power(xii,32.0)*
            (-401800.0 - 18435.0*Power(rij,2.0)*Power(xij,2.0) + 
              52.0*Power(rij,4.0)*Power(xij,4.0)) + 
           72.0*Power(rij,5.0)*Power(xii,31.0)*
            (-874650.0 - 113005.0*Power(rij,2.0)*Power(xij,2.0) + 
              811.0*Power(rij,4.0)*Power(xij,4.0)) + 
           450.0*rij*Power(xii,7.0)*Power(xij,20.0)*
            (-305235.0 + 35910.0*Power(rij,2.0)*Power(xij,2.0) - 
              798.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           1575.0*Power(xii,6.0)*Power(xij,20.0)*
            (-43605.0 + 15390.0*Power(rij,2.0)*Power(xij,2.0) - 
              570.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           36.0*Power(rij,4.0)*Power(xii,30.0)*
            (10870650.0 + 3565660.0*Power(rij,2.0)*Power(xij,2.0) - 
              40547.0*Power(rij,4.0)*Power(xij,4.0) + 74.0*Power(rij,6.0)*Power(xij,6.0)) - 
           216.0*Power(rij,3.0)*Power(xii,29.0)*
            (7800450.0 + 6398630.0*Power(rij,2.0)*Power(xij,2.0) - 
              76456.0*Power(rij,4.0)*Power(xij,4.0) + 383.0*Power(rij,6.0)*Power(xij,6.0)) + 
           100.0*rij*Power(xii,9.0)*Power(xij,18.0)*
            (5494230.0 - 915705.0*Power(rij,2.0)*Power(xij,2.0) + 
              32319.0*Power(rij,4.0)*Power(xij,4.0) - 342.0*Power(rij,6.0)*Power(xij,6.0) + 
              Power(rij,8.0)*Power(xij,8.0)) + 
           450.0*Power(xii,8.0)*Power(xij,18.0)*
            (610470.0 - 305235.0*Power(rij,2.0)*Power(xij,2.0) + 
              17955.0*Power(rij,4.0)*Power(xij,4.0) - 266.0*Power(rij,6.0)*Power(xij,6.0) + 
              Power(rij,8.0)*Power(xij,8.0)) - 
           36.0*rij*Power(xii,13.0)*Power(xij,14.0)*
            (89026875.0 + 334741050.0*Power(rij,2.0)*Power(xij,2.0) + 
              16373490.0*Power(rij,4.0)*Power(xij,4.0) + 
              300416.0*Power(rij,6.0)*Power(xij,6.0) + 42.0*Power(rij,8.0)*Power(xij,8.0)) - 
           12.0*rij*Power(xii,11.0)*Power(xij,16.0)*
            (137355750.0 - 22865850.0*Power(rij,2.0)*Power(xij,2.0) + 
              2246895.0*Power(rij,4.0)*Power(xij,4.0) - 
              7632.0*Power(rij,6.0)*Power(xij,6.0) + 305.0*Power(rij,8.0)*Power(xij,8.0)) + 
           4.0*Power(rij,2.0)*Power(xii,28.0)*
            (-1167523875.0 - 2573131050.0*Power(rij,2.0)*Power(xij,2.0) + 
              700245.0*Power(rij,4.0)*Power(xij,4.0) + 
              110763.0*Power(rij,6.0)*Power(xij,6.0) + 409.0*Power(rij,8.0)*Power(xij,8.0)) \
    + 36.0*rij*Power(xii,15.0)*Power(xij,12.0)*
            (-7803332775.0 - 3848335050.0*Power(rij,2.0)*Power(xij,2.0) - 
              104502090.0*Power(rij,4.0)*Power(xij,4.0) + 
              530994.0*Power(rij,6.0)*Power(xij,6.0) + 2438.0*Power(rij,8.0)*Power(xij,8.0)) \
    + 18.0*rij*Power(xii,19.0)*Power(xij,8.0)*
            (-357591274425.0 + 15942280410.0*Power(rij,2.0)*Power(xij,2.0) + 
              570992940.0*Power(rij,4.0)*Power(xij,4.0) - 
              11648672.0*Power(rij,6.0)*Power(xij,6.0) + 
              12046.0*Power(rij,8.0)*Power(xij,8.0)) + 
           36.0*rij*Power(xii,25.0)*Power(xij,2.0)*
            (-7604930025.0 - 7728532245.0*Power(rij,2.0)*Power(xij,2.0) + 
              323602755.0*Power(rij,4.0)*Power(xij,4.0) - 
              5106594.0*Power(rij,6.0)*Power(xij,6.0) + 
              16967.0*Power(rij,8.0)*Power(xij,8.0)) - 
           4.0*rij*Power(xii,27.0)*(1762732125.0 + 
              12948068700.0*Power(rij,2.0)*Power(xij,2.0) + 
              634368105.0*Power(rij,4.0)*Power(xij,4.0) - 
              14850144.0*Power(rij,6.0)*Power(xij,6.0) + 
              38437.0*Power(rij,8.0)*Power(xij,8.0)) + 
           6.0*rij*Power(xii,21.0)*Power(xij,6.0)*
            (-1070591393025.0 + 112013496630.0*Power(rij,2.0)*Power(xij,2.0) - 
              3943648170.0*Power(rij,4.0)*Power(xij,4.0) + 
              24403296.0*Power(rij,6.0)*Power(xij,6.0) + 
              43928.0*Power(rij,8.0)*Power(xij,8.0)) - 
           6.0*rij*Power(xii,23.0)*Power(xij,4.0)*
            (387584799525.0 + 22847636070.0*Power(rij,2.0)*Power(xij,2.0) - 
              380961630.0*Power(rij,4.0)*Power(xij,4.0) - 
              17523396.0*Power(rij,6.0)*Power(xij,6.0) + 
              121448.0*Power(rij,8.0)*Power(xij,8.0)) - 
           2.0*rij*Power(xii,17.0)*Power(xij,10.0)*
            (1157397524325.0 + 169929208260.0*Power(rij,2.0)*Power(xij,2.0) - 
              3903345180.0*Power(rij,4.0)*Power(xij,4.0) - 
              33203412.0*Power(rij,6.0)*Power(xij,6.0) + 
              129226.0*Power(rij,8.0)*Power(xij,8.0)) + 
           10.0*Power(xii,10.0)*Power(xij,16.0)*
            (-82413450.0 + 54942300.0*Power(rij,2.0)*Power(xij,2.0) - 
              4578525.0*Power(rij,4.0)*Power(xij,4.0) + 
              107730.0*Power(rij,6.0)*Power(xij,6.0) - 
              855.0*Power(rij,8.0)*Power(xij,8.0) + 2.0*Power(rij,10.0)*Power(xij,10.0)) + 
           6.0*Power(xii,12.0)*Power(xij,14.0)*
            (-267080625.0 + 122094000.0*Power(rij,2.0)*Power(xij,2.0) + 
              114757650.0*Power(rij,4.0)*Power(xij,4.0) + 
              2105670.0*Power(rij,6.0)*Power(xij,6.0) + 
              45858.0*Power(rij,8.0)*Power(xij,8.0) + 10.0*Power(rij,10.0)*Power(xij,10.0)) \
    - 18.0*Power(xii,14.0)*Power(xij,12.0)*
            (7803332775.0 - 5926646250.0*Power(rij,2.0)*Power(xij,2.0) - 
              670499550.0*Power(rij,4.0)*Power(xij,4.0) - 
              13909140.0*Power(rij,6.0)*Power(xij,6.0) + 
              16074.0*Power(rij,8.0)*Power(xij,8.0) + 52.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 9.0*Power(xii,16.0)*Power(xij,10.0)*
            (-128599724925.0 + 98821662660.0*Power(rij,2.0)*Power(xij,2.0) + 
              4512729900.0*Power(rij,4.0)*Power(xij,4.0) - 
              57727320.0*Power(rij,6.0)*Power(xij,6.0) - 
              361898.0*Power(rij,8.0)*Power(xij,8.0) + 296.0*Power(rij,10.0)*Power(xij,10.0)\
    ) + 18.0*Power(xii,26.0)*(-195859125.0 - 9105160050.0*Power(rij,2.0)*Power(xij,2.0) - 
              2073138165.0*Power(rij,4.0)*Power(xij,4.0) + 
              67371010.0*Power(rij,6.0)*Power(xij,6.0) - 
              517501.0*Power(rij,8.0)*Power(xij,8.0) + 342.0*Power(rij,10.0)*Power(xij,10.0)\
    ) - 9.0*Power(xii,20.0)*Power(xij,6.0)*
            (356863797675.0 - 32148394110.0*Power(rij,2.0)*Power(xij,2.0) + 
              11000178630.0*Power(rij,4.0)*Power(xij,4.0) - 
              307586440.0*Power(rij,6.0)*Power(xij,6.0) + 
              922488.0*Power(rij,8.0)*Power(xij,8.0) + 684.0*Power(rij,10.0)*Power(xij,10.0)\
    ) - Power(xii,18.0)*Power(xij,8.0)*
            (3218321469825.0 - 2072966838390.0*Power(rij,2.0)*Power(xij,2.0) + 
              68132114820.0*Power(rij,4.0)*Power(xij,4.0) + 
              914077080.0*Power(rij,6.0)*Power(xij,6.0) - 
              9941598.0*Power(rij,8.0)*Power(xij,8.0) + 
              1636.0*Power(rij,10.0)*Power(xij,10.0)) - 
           6.0*Power(xii,24.0)*Power(xij,2.0)*
            (22814790075.0 + 192532308885.0*Power(rij,2.0)*Power(xij,2.0) - 
              7524596205.0*Power(rij,4.0)*Power(xij,4.0) + 
              238725270.0*Power(rij,6.0)*Power(xij,6.0) - 
              2550237.0*Power(rij,8.0)*Power(xij,8.0) + 
              2584.0*Power(rij,10.0)*Power(xij,10.0)) + 
           3.0*Power(xii,22.0)*Power(xij,4.0)*
            (-387584799525.0 - 678595490370.0*Power(rij,2.0)*Power(xij,2.0) + 
              38884316310.0*Power(rij,4.0)*Power(xij,4.0) - 
              415382940.0*Power(rij,6.0)*Power(xij,6.0) - 
              1961256.0*Power(rij,8.0)*Power(xij,8.0) + 
              5168.0*Power(rij,10.0)*Power(xij,10.0))))/
      (70875.*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),19.0))
    ;
  }
  return S;
}

static double DSlater_5S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     (-766364132575641600000.0 + 766364132575641600000.0*Power(E,2.0*rij*xii) - 
        1532728265151283200000.0*rij*xii - 
        1532728265151283200000.0*Power(rij,2.0)*Power(xii,2.0) - 
        1020454892919846759375.0*Power(rij,3.0)*Power(xii,3.0) - 
        508181520688410318750.0*Power(rij,4.0)*Power(xii,4.0) - 
        201658839456289965750.0*Power(rij,5.0)*Power(xii,5.0) - 
        66348599139429106500.0*Power(rij,6.0)*Power(xii,6.0) - 
        18599424978069936000.0*Power(rij,7.0)*Power(xii,7.0) - 
        4531844261934990000.0*Power(rij,8.0)*Power(xii,8.0) - 
        974502735982776000.0*Power(rij,9.0)*Power(xii,9.0) - 
        187178934377635200.0*Power(rij,10.0)*Power(xii,10.0) - 
        32426214068102400.0*Power(rij,11.0)*Power(xii,11.0) - 
        5105770650489600.0*Power(rij,12.0)*Power(xii,12.0) - 
        735142625280000.0*Power(rij,13.0)*Power(xii,13.0) - 
        97218861465600.0*Power(rij,14.0)*Power(xii,14.0) - 
        11842297528320.0*Power(rij,15.0)*Power(xii,15.0) - 
        1330258083840.0*Power(rij,16.0)*Power(xii,16.0) - 
        137673768960.0*Power(rij,17.0)*Power(xii,17.0) - 
        13074432000.0*Power(rij,18.0)*Power(xii,18.0) - 
        1128529920.0*Power(rij,19.0)*Power(xii,19.0) - 
        86704128.0*Power(rij,20.0)*Power(xii,20.0) - 
        5636096.0*Power(rij,21.0)*Power(xii,21.0) - 262144.0*Power(rij,22.0)*Power(xii,22.0))/
      (7.663641325756416e20*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (2338875.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),21.0) + 
        55.0*Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-816.0*Power(rij,9.0)*Power(xii,37.0) - 12.0*Power(rij,10.0)*Power(xii,38.0) + 
           42525.0*Power(xij,28.0) + 85050.0*rij*xii*Power(xij,28.0) + 
           28350.0*rij*Power(xii,3.0)*Power(xij,26.0)*
            (-63.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           42525.0*Power(xii,2.0)*Power(xij,26.0)*
            (-21.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           16.0*Power(rij,8.0)*Power(xii,36.0)*(1647.0 + 7.0*Power(rij,2.0)*Power(xij,2.0)) - 
           4.0*Power(rij,7.0)*Power(xii,35.0)*
            (132192.0 + 2723.0*Power(rij,2.0)*Power(xij,2.0)) + 
           11340.0*rij*Power(xii,5.0)*Power(xij,24.0)*
            (1575.0 - 105.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) \
    + 28350.0*Power(xii,4.0)*Power(xij,24.0)*
            (315.0 - 63.0*Power(rij,2.0)*Power(xij,2.0) + Power(rij,4.0)*Power(xij,4.0)) + 
           378.0*Power(rij,6.0)*Power(xii,34.0)*
            (-19176.0 - 1227.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           756.0*Power(rij,5.0)*Power(xii,33.0)*
            (-93024.0 - 15446.0*Power(rij,2.0)*Power(xij,2.0) + 
              55.0*Power(rij,4.0)*Power(xij,4.0)) + 
           540.0*rij*Power(xii,7.0)*Power(xij,22.0)*
            (-209475.0 + 22050.0*Power(rij,2.0)*Power(xij,2.0) - 
              441.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) + 
           1890.0*Power(xii,6.0)*Power(xij,22.0)*
            (-29925.0 + 9450.0*Power(rij,2.0)*Power(xij,2.0) - 
              315.0*Power(rij,4.0)*Power(xij,4.0) + 2.0*Power(rij,6.0)*Power(xij,6.0)) - 
           18.0*Power(rij,4.0)*Power(xii,32.0)*
            (26860680.0 + 10787154.0*Power(rij,2.0)*Power(xij,2.0) - 
              38451.0*Power(rij,4.0)*Power(xij,4.0) + 32.0*Power(rij,6.0)*Power(xij,6.0)) + 
           36.0*Power(rij,3.0)*Power(xii,31.0)*
            (-63488880.0 - 61912431.0*Power(rij,2.0)*Power(xij,2.0) - 
              191226.0*Power(rij,4.0)*Power(xij,4.0) + 1798.0*Power(rij,6.0)*Power(xij,6.0)\
    ) + 30.0*rij*Power(xii,9.0)*Power(xij,20.0)*
            (16967475.0 - 2513700.0*Power(rij,2.0)*Power(xij,2.0) + 
              79380.0*Power(rij,4.0)*Power(xij,4.0) - 
              756.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           135.0*Power(xii,8.0)*Power(xij,20.0)*
            (1885275.0 - 837900.0*Power(rij,2.0)*Power(xij,2.0) + 
              44100.0*Power(rij,4.0)*Power(xij,4.0) - 
              588.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           378.0*rij*Power(xii,25.0)*Power(xij,4.0)*
            (17962854525.0 + 3973121100.0*Power(rij,2.0)*Power(xij,2.0) - 
              123555420.0*Power(rij,4.0)*Power(xij,4.0) + 
              821732.0*Power(rij,6.0)*Power(xij,6.0) + 38.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 108.0*rij*Power(xii,13.0)*Power(xij,16.0)*
            (121076550.0 + 226721775.0*Power(rij,2.0)*Power(xij,2.0) + 
              10259361.0*Power(rij,4.0)*Power(xij,4.0) + 
              137164.0*Power(rij,6.0)*Power(xij,6.0) + 126.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 14.0*rij*Power(xii,11.0)*Power(xij,18.0)*
            (123620175.0 - 12901950.0*Power(rij,2.0)*Power(xij,2.0) + 
              1936548.0*Power(rij,4.0)*Power(xij,4.0) + 
              648.0*Power(rij,6.0)*Power(xij,6.0) + 194.0*Power(rij,8.0)*Power(xij,8.0)) - 
           756.0*rij*Power(xii,23.0)*Power(xij,6.0)*
            (37244490525.0 - 1707788550.0*Power(rij,2.0)*Power(xij,2.0) + 
              71446053.0*Power(rij,4.0)*Power(xij,4.0) - 
              841738.0*Power(rij,6.0)*Power(xij,6.0) + 1292.0*Power(rij,8.0)*Power(xij,8.0)\
    ) - 210.0*rij*Power(xii,19.0)*Power(xij,10.0)*
            (154925948835.0 + 1186358670.0*Power(rij,2.0)*Power(xij,2.0) - 
              250536888.0*Power(rij,4.0)*Power(xij,4.0) + 
              1233936.0*Power(rij,6.0)*Power(xij,6.0) + 
              1498.0*Power(rij,8.0)*Power(xij,8.0)) - 
           2.0*Power(rij,2.0)*Power(xii,30.0)*
            (3444882210.0 + 8888317515.0*Power(rij,2.0)*Power(xij,2.0) + 
              249425190.0*Power(rij,4.0)*Power(xij,4.0) - 
              2687310.0*Power(rij,6.0)*Power(xij,6.0) + 
              2318.0*Power(rij,8.0)*Power(xij,8.0)) + 
           36.0*rij*Power(xii,15.0)*Power(xij,14.0)*
            (-23418646650.0 - 10839912300.0*Power(rij,2.0)*Power(xij,2.0) - 
              280275345.0*Power(rij,4.0)*Power(xij,4.0) - 
              161206.0*Power(rij,6.0)*Power(xij,6.0) + 2966.0*Power(rij,8.0)*Power(xij,8.0)\
    ) + 252.0*rij*Power(xii,21.0)*Power(xij,8.0)*
            (-186637212225.0 + 12975560325.0*Power(rij,2.0)*Power(xij,2.0) - 
              160724781.0*Power(rij,4.0)*Power(xij,4.0) - 
              824296.0*Power(rij,6.0)*Power(xij,6.0) + 4104.0*Power(rij,8.0)*Power(xij,8.0)\
    ) - 20.0*rij*Power(xii,29.0)*(556016076.0 + 
              4792896927.0*Power(rij,2.0)*Power(xij,2.0) + 
              502574625.0*Power(rij,4.0)*Power(xij,4.0) - 
              7764012.0*Power(rij,6.0)*Power(xij,6.0) + 
              23578.0*Power(rij,8.0)*Power(xij,8.0)) - 
           2.0*rij*Power(xii,17.0)*Power(xij,12.0)*
            (4434624921825.0 + 748179822600.0*Power(rij,2.0)*Power(xij,2.0) - 
              2615253480.0*Power(rij,4.0)*Power(xij,4.0) - 
              94072860.0*Power(rij,6.0)*Power(xij,6.0) + 
              66566.0*Power(rij,8.0)*Power(xij,8.0)) + 
           6.0*rij*Power(xii,27.0)*Power(xij,2.0)*
            (-95426713305.0 - 133474043430.0*Power(rij,2.0)*Power(xij,2.0) + 
              2257034220.0*Power(rij,4.0)*Power(xij,4.0) - 
              27096360.0*Power(rij,6.0)*Power(xij,6.0) + 
              115178.0*Power(rij,8.0)*Power(xij,8.0)) + 
           189.0*Power(xii,24.0)*Power(xij,4.0)*
            (-17962854525.0 - 54465790500.0*Power(rij,2.0)*Power(xij,2.0) + 
              2436787500.0*Power(rij,4.0)*Power(xij,4.0) - 
              39184124.0*Power(rij,6.0)*Power(xij,6.0) + 
              135242.0*Power(rij,8.0)*Power(xij,8.0)) + 
           3.0*Power(xii,10.0)*Power(xij,18.0)*
            (-288447075.0 + 169674750.0*Power(rij,2.0)*Power(xij,2.0) - 
              12568500.0*Power(rij,4.0)*Power(xij,4.0) + 
              264600.0*Power(rij,6.0)*Power(xij,6.0) - 
              1890.0*Power(rij,8.0)*Power(xij,8.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) + 
           14.0*Power(xii,12.0)*Power(xij,16.0)*
            (-467009550.0 + 233504775.0*Power(rij,2.0)*Power(xij,2.0) + 
              74474775.0*Power(rij,4.0)*Power(xij,4.0) + 
              1502604.0*Power(rij,6.0)*Power(xij,6.0) + 
              19494.0*Power(rij,8.0)*Power(xij,8.0) + 8.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 378.0*Power(xii,22.0)*Power(xij,6.0)*
            (-37244490525.0 - 19531480050.0*Power(rij,2.0)*Power(xij,2.0) + 
              314001975.0*Power(rij,4.0)*Power(xij,4.0) + 
              10131142.0*Power(rij,6.0)*Power(xij,6.0) - 
              89148.0*Power(rij,8.0)*Power(xij,8.0) + 38.0*Power(rij,10.0)*Power(xij,10.0)) \
    - 18.0*Power(xii,14.0)*Power(xij,14.0)*
            (23418646650.0 - 15063347250.0*Power(rij,2.0)*Power(xij,2.0) - 
              1541945475.0*Power(rij,4.0)*Power(xij,4.0) - 
              26695718.0*Power(rij,6.0)*Power(xij,6.0) - 
              30126.0*Power(rij,8.0)*Power(xij,8.0) + 42.0*Power(rij,10.0)*Power(xij,10.0)) \
    + 9.0*Power(xii,16.0)*Power(xij,12.0)*
            (-492736102425.0 + 332913709800.0*Power(rij,2.0)*Power(xij,2.0) + 
              16232804280.0*Power(rij,4.0)*Power(xij,4.0) - 
              17642940.0*Power(rij,6.0)*Power(xij,6.0) - 
              670806.0*Power(rij,8.0)*Power(xij,8.0) + 64.0*Power(rij,10.0)*Power(xij,10.0)\
    ) - 63.0*Power(xii,26.0)*Power(xij,2.0)*
            (4544129205.0 + 52573689630.0*Power(rij,2.0)*Power(xij,2.0) + 
              831772260.0*Power(rij,4.0)*Power(xij,4.0) - 
              17171800.0*Power(rij,6.0)*Power(xij,6.0) - 
              71082.0*Power(rij,8.0)*Power(xij,8.0) + 228.0*Power(rij,10.0)*Power(xij,10.0)\
    ) - 42.0*Power(xii,20.0)*Power(xij,8.0)*
            (559911636675.0 - 188981369325.0*Power(rij,2.0)*Power(xij,2.0) + 
              12441827475.0*Power(rij,4.0)*Power(xij,4.0) - 
              104458200.0*Power(rij,6.0)*Power(xij,6.0) - 
              229140.0*Power(rij,8.0)*Power(xij,8.0) + 
              328.0*Power(rij,10.0)*Power(xij,10.0)) + 
           6.0*Power(xii,28.0)*(-926693460.0 - 
              54014259915.0*Power(rij,2.0)*Power(xij,2.0) - 
              19074253275.0*Power(rij,4.0)*Power(xij,4.0) + 
              370738620.0*Power(rij,6.0)*Power(xij,6.0) - 
              2470770.0*Power(rij,8.0)*Power(xij,8.0) + 
              2296.0*Power(rij,10.0)*Power(xij,10.0)) + 
           Power(xii,18.0)*Power(xij,10.0)*
            (-16267224627675.0 + 10111789349550.0*Power(rij,2.0)*Power(xij,2.0) - 
              46442115720.0*Power(rij,4.0)*Power(xij,4.0) - 
              3760495200.0*Power(rij,6.0)*Power(xij,6.0) + 
              8510670.0*Power(rij,8.0)*Power(xij,8.0) + 
              4636.0*Power(rij,10.0)*Power(xij,10.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,12.0)*
         (5.0*Power(xii,28.0)*Power(xij,2.0)*
            (9823275.0 + 19646550.0*rij*xij + 19646550.0*Power(rij,2.0)*Power(xij,2.0) + 
              13097700.0*Power(rij,3.0)*Power(xij,3.0) + 
              6548850.0*Power(rij,4.0)*Power(xij,4.0) + 
              2619540.0*Power(rij,5.0)*Power(xij,5.0) + 
              873180.0*Power(rij,6.0)*Power(xij,6.0) + 
              249480.0*Power(rij,7.0)*Power(xij,7.0) + 
              62370.0*Power(rij,8.0)*Power(xij,8.0) + 
              13860.0*Power(rij,9.0)*Power(xij,9.0) + 
              2772.0*Power(rij,10.0)*Power(xij,10.0) + 
              944.0*Power(rij,11.0)*Power(xij,11.0) - 4.0*Power(rij,12.0)*Power(xij,12.0)) + 
           52668.0*Power(xii,16.0)*Power(xij,14.0)*
            (10970100.0 + 21940200.0*rij*xij - 
              237536325.0*Power(rij,2.0)*Power(xij,2.0) + 
              487687500.0*Power(rij,3.0)*Power(xij,3.0) - 
              253716075.0*Power(rij,4.0)*Power(xij,4.0) + 
              8373180.0*Power(rij,5.0)*Power(xij,5.0) + 
              14379960.0*Power(rij,6.0)*Power(xij,6.0) + 
              690840.0*Power(rij,7.0)*Power(xij,7.0) - 
              198500.0*Power(rij,8.0)*Power(xij,8.0) - 
              23510.0*Power(rij,9.0)*Power(xij,9.0) - 
              603.0*Power(rij,10.0)*Power(xij,10.0) + 
              24.0*Power(rij,11.0)*Power(xij,11.0) + Power(rij,12.0)*Power(xij,12.0)) + 
           20.0*Power(xij,30.0)*(24325703325.0 + 48651406650.0*rij*xij + 
              34902096075.0*Power(rij,2.0)*Power(xij,2.0) + 
              14101857000.0*Power(rij,3.0)*Power(xij,3.0) + 
              3777283125.0*Power(rij,4.0)*Power(xij,4.0) + 
              725238360.0*Power(rij,5.0)*Power(xij,5.0) + 
              103908420.0*Power(rij,6.0)*Power(xij,6.0) + 
              11309760.0*Power(rij,7.0)*Power(xij,7.0) + 
              935550.0*Power(rij,8.0)*Power(xij,8.0) + 
              57750.0*Power(rij,9.0)*Power(xij,9.0) + 
              2541.0*Power(rij,10.0)*Power(xij,10.0) + 
              72.0*Power(rij,11.0)*Power(xij,11.0) + Power(rij,12.0)*Power(xij,12.0)) + 
           20.0*Power(xii,2.0)*Power(xij,28.0)*
            (1145624789925.0 + 2291249579850.0*rij*xij + 
              1509654155625.0*Power(rij,2.0)*Power(xij,2.0) + 
              544835317950.0*Power(rij,3.0)*Power(xij,3.0) + 
              127873624725.0*Power(rij,4.0)*Power(xij,4.0) + 
              21112494480.0*Power(rij,5.0)*Power(xij,5.0) + 
              2545402860.0*Power(rij,6.0)*Power(xij,6.0) + 
              226801080.0*Power(rij,7.0)*Power(xij,7.0) + 
              14787630.0*Power(rij,8.0)*Power(xij,8.0) + 
              679690.0*Power(rij,9.0)*Power(xij,9.0) + 
              20207.0*Power(rij,10.0)*Power(xij,10.0) + 
              314.0*Power(rij,11.0)*Power(xij,11.0) + Power(rij,12.0)*Power(xij,12.0)) - 
           110.0*Power(xii,10.0)*Power(xij,20.0)*
            (-7031709085275.0 - 14063418170550.0*rij*xij + 
              5761420840350.0*Power(rij,2.0)*Power(xij,2.0) + 
              2282136160500.0*Power(rij,3.0)*Power(xij,3.0) - 
              28040219550.0*Power(rij,4.0)*Power(xij,4.0) - 
              94545378900.0*Power(rij,5.0)*Power(xij,5.0) - 
              14293824300.0*Power(rij,6.0)*Power(xij,6.0) - 
              573941160.0*Power(rij,7.0)*Power(xij,7.0) + 
              65493270.0*Power(rij,8.0)*Power(xij,8.0) + 
              9429840.0*Power(rij,9.0)*Power(xij,9.0) + 
              476262.0*Power(rij,10.0)*Power(xij,10.0) + 
              9284.0*Power(rij,11.0)*Power(xij,11.0) + 2.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 55.0*Power(xii,20.0)*Power(xij,10.0)*
            (865341225.0 + 1730682450.0*rij*xij + 
              1730682450.0*Power(rij,2.0)*Power(xij,2.0) + 
              1312510500.0*Power(rij,3.0)*Power(xij,3.0) - 
              534161250.0*Power(rij,4.0)*Power(xij,4.0) + 
              2257389540.0*Power(rij,5.0)*Power(xij,5.0) - 
              1135495620.0*Power(rij,6.0)*Power(xij,6.0) + 
              166186440.0*Power(rij,7.0)*Power(xij,7.0) + 
              48426210.0*Power(rij,8.0)*Power(xij,8.0) - 
              544460.0*Power(rij,9.0)*Power(xij,9.0) - 
              446788.0*Power(rij,10.0)*Power(xij,10.0) - 
              18256.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 5.0*Power(xii,30.0)*(467775.0 + 935550.0*rij*xij + 
              935550.0*Power(rij,2.0)*Power(xij,2.0) + 
              623700.0*Power(rij,3.0)*Power(xij,3.0) + 
              311850.0*Power(rij,4.0)*Power(xij,4.0) + 
              124740.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              11880.0*Power(rij,7.0)*Power(xij,7.0) + 
              2970.0*Power(rij,8.0)*Power(xij,8.0) + 660.0*Power(rij,9.0)*Power(xij,9.0) + 
              132.0*Power(rij,10.0)*Power(xij,10.0) + 
              24.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           13167.0*Power(xii,14.0)*Power(xij,16.0)*
            (-1159677225.0 - 2319354450.0*rij*xij + 
              12207865950.0*Power(rij,2.0)*Power(xij,2.0) - 
              8468827500.0*Power(rij,3.0)*Power(xij,3.0) + 
              270633450.0*Power(rij,4.0)*Power(xij,4.0) + 
              576841620.0*Power(rij,5.0)*Power(xij,5.0) + 
              36246140.0*Power(rij,6.0)*Power(xij,6.0) - 
              9802040.0*Power(rij,7.0)*Power(xij,7.0) - 
              1523150.0*Power(rij,8.0)*Power(xij,8.0) - 
              59740.0*Power(rij,9.0)*Power(xij,9.0) + 
              2148.0*Power(rij,10.0)*Power(xij,10.0) + 
              216.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           1540.0*Power(xii,18.0)*Power(xij,12.0)*
            (82413450.0 + 164826900.0*rij*xij + 
              315918225.0*Power(rij,2.0)*Power(xij,2.0) - 
              778349250.0*Power(rij,3.0)*Power(xij,3.0) + 
              1488020625.0*Power(rij,4.0)*Power(xij,4.0) - 
              723888900.0*Power(rij,5.0)*Power(xij,5.0) + 
              53335800.0*Power(rij,6.0)*Power(xij,6.0) + 
              35328960.0*Power(rij,7.0)*Power(xij,7.0) + 
              727920.0*Power(rij,8.0)*Power(xij,8.0) - 
              407850.0*Power(rij,9.0)*Power(xij,9.0) - 
              31161.0*Power(rij,10.0)*Power(xij,10.0) - 
              342.0*Power(rij,11.0)*Power(xij,11.0) + 17.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 385.0*Power(xii,12.0)*Power(xij,18.0)*
            (486782696925.0 + 973565393850.0*rij*xij - 
              1578263239350.0*Power(rij,2.0)*Power(xij,2.0) + 
              159440683500.0*Power(rij,3.0)*Power(xij,3.0) + 
              142522509150.0*Power(rij,4.0)*Power(xij,4.0) + 
              7396549020.0*Power(rij,5.0)*Power(xij,5.0) - 
              3365814060.0*Power(rij,6.0)*Power(xij,6.0) - 
              566536680.0*Power(rij,7.0)*Power(xij,7.0) - 
              26389530.0*Power(rij,8.0)*Power(xij,8.0) + 
              1233420.0*Power(rij,9.0)*Power(xij,9.0) + 
              180876.0*Power(rij,10.0)*Power(xij,10.0) + 
              6672.0*Power(rij,11.0)*Power(xij,11.0) + 68.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 14.0*Power(xii,26.0)*Power(xij,4.0)*
            (-35083125.0 - 70166250.0*rij*xij - 
              70166250.0*Power(rij,2.0)*Power(xij,2.0) - 
              46777500.0*Power(rij,3.0)*Power(xij,3.0) - 
              23388750.0*Power(rij,4.0)*Power(xij,4.0) - 
              9355500.0*Power(rij,5.0)*Power(xij,5.0) - 
              3118500.0*Power(rij,6.0)*Power(xij,6.0) - 
              891000.0*Power(rij,7.0)*Power(xij,7.0) - 
              222750.0*Power(rij,8.0)*Power(xij,8.0) - 
              18040.0*Power(rij,9.0)*Power(xij,9.0) - 
              31922.0*Power(rij,10.0)*Power(xij,10.0) - 
              524.0*Power(rij,11.0)*Power(xij,11.0) + 74.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 7.0*Power(xii,4.0)*Power(xij,26.0)*
            (-34911472624875.0 - 69822945249750.0*rij*xij - 
              39657058575750.0*Power(rij,2.0)*Power(xij,2.0) - 
              11594155765500.0*Power(rij,3.0)*Power(xij,3.0) - 
              2053437759450.0*Power(rij,4.0)*Power(xij,4.0) - 
              226760230620.0*Power(rij,5.0)*Power(xij,5.0) - 
              13566847140.0*Power(rij,6.0)*Power(xij,6.0) + 
              59554440.0*Power(rij,7.0)*Power(xij,7.0) + 
              92560050.0*Power(rij,8.0)*Power(xij,8.0) + 
              8960820.0*Power(rij,9.0)*Power(xij,9.0) + 
              454476.0*Power(rij,10.0)*Power(xij,10.0) + 
              12592.0*Power(rij,11.0)*Power(xij,11.0) + 
              148.0*Power(rij,12.0)*Power(xij,12.0)) - 
           50.0*Power(xii,8.0)*Power(xij,22.0)*
            (-25884416787075.0 - 51768833574150.0*rij*xij - 
              7250000130450.0*Power(rij,2.0)*Power(xij,2.0) + 
              3206573300700.0*Power(rij,3.0)*Power(xij,3.0) + 
              1119709961370.0*Power(rij,4.0)*Power(xij,4.0) + 
              124678104012.0*Power(rij,5.0)*Power(xij,5.0) - 
              603805356.0*Power(rij,6.0)*Power(xij,6.0) - 
              1678543416.0*Power(rij,7.0)*Power(xij,7.0) - 
              199885554.0*Power(rij,8.0)*Power(xij,8.0) - 
              11203984.0*Power(rij,9.0)*Power(xij,9.0) - 
              267410.0*Power(rij,10.0)*Power(xij,10.0) + 
              1696.0*Power(rij,11.0)*Power(xij,11.0) + 158.0*Power(rij,12.0)*Power(xij,12.0)\
    ) + 25.0*Power(xii,22.0)*Power(xij,8.0)*
            (-559926675.0 - 1119853350.0*rij*xij - 
              1119853350.0*Power(rij,2.0)*Power(xij,2.0) - 
              746568900.0*Power(rij,3.0)*Power(xij,3.0) - 
              373284450.0*Power(rij,4.0)*Power(xij,4.0) - 
              108881388.0*Power(rij,5.0)*Power(xij,5.0) - 
              171068436.0*Power(rij,6.0)*Power(xij,6.0) + 
              84312360.0*Power(rij,7.0)*Power(xij,7.0) - 
              23318262.0*Power(rij,8.0)*Power(xij,8.0) - 
              3243900.0*Power(rij,9.0)*Power(xij,9.0) + 
              144012.0*Power(rij,10.0)*Power(xij,10.0) + 
              21256.0*Power(rij,11.0)*Power(xij,11.0) + 
              316.0*Power(rij,12.0)*Power(xij,12.0)) - 
           6.0*Power(xii,24.0)*Power(xij,6.0)*
            (-518450625.0 - 1036901250.0*rij*xij - 
              1036901250.0*Power(rij,2.0)*Power(xij,2.0) - 
              691267500.0*Power(rij,3.0)*Power(xij,3.0) - 
              345633750.0*Power(rij,4.0)*Power(xij,4.0) - 
              138253500.0*Power(rij,5.0)*Power(xij,5.0) - 
              46084500.0*Power(rij,6.0)*Power(xij,6.0) - 
              18829800.0*Power(rij,7.0)*Power(xij,7.0) + 
              5202450.0*Power(rij,8.0)*Power(xij,8.0) - 
              3594360.0*Power(rij,9.0)*Power(xij,9.0) - 
              218658.0*Power(rij,10.0)*Power(xij,10.0) + 
              17664.0*Power(rij,11.0)*Power(xij,11.0) + 
              766.0*Power(rij,12.0)*Power(xij,12.0)) + 
           3.0*Power(xii,6.0)*Power(xij,24.0)*
            (298003296331125.0 + 596006592662250.0*rij*xij + 
              248020874120250.0*Power(rij,2.0)*Power(xij,2.0) + 
              39730277317500.0*Power(rij,3.0)*Power(xij,3.0) + 
              92486497950.0*Power(rij,4.0)*Power(xij,4.0) - 
              1040164489980.0*Power(rij,5.0)*Power(xij,5.0) - 
              195369465060.0*Power(rij,6.0)*Power(xij,6.0) - 
              18560397240.0*Power(rij,7.0)*Power(xij,7.0) - 
              940920750.0*Power(rij,8.0)*Power(xij,8.0) - 
              12201420.0*Power(rij,9.0)*Power(xij,9.0) + 
              1418604.0*Power(rij,10.0)*Power(xij,10.0) + 
              84168.0*Power(rij,11.0)*Power(xij,11.0) + 
              1532.0*Power(rij,12.0)*Power(xij,12.0))))/
      (2.338875e6*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),21.0))
    ;
  }
  return S;
}

double DSlater_5S_1S(double rij,double xii,double xij)
{
  return DSlater_1S_5S(rij,xij,xii);
}

double DSlater_5S_2S(double rij,double xii,double xij)
{
  return DSlater_2S_5S(rij,xij,xii);
}

double DSlater_5S_3S(double rij,double xii,double xij)
{
  return DSlater_3S_5S(rij,xij,xii);
}

double DSlater_5S_4S(double rij,double xii,double xij)
{
  return DSlater_4S_5S(rij,xij,xii);
}

static double DSlater_6S_6S(double rij,double xij,double xii)
{
  double S;

  if (xii == xij) {
    S =     -(930672602599859159040000.0 - 930672602599859159040000.0*Power(E,2.0*rij*xii) + 
         1861345205199718318080000.0*rij*xii + 
         1861345205199718318080000.0*Power(rij,2.0)*Power(xii,2.0) + 
         1239447469649939026351875.0*Power(rij,3.0)*Power(xii,3.0) + 
         617549734100159734623750.0*Power(rij,4.0)*Power(xii,4.0) + 
         245308299418626353910000.0*Power(rij,5.0)*Power(xii,5.0) + 
         80849221192532687895000.0*Power(rij,6.0)*Power(xii,6.0) + 
         22724497062484591374000.0*Power(rij,7.0)*Power(xii,7.0) + 
         5558106457968308244000.0*Power(rij,8.0)*Power(xii,8.0) + 
         1201461043722619680000.0*Power(rij,9.0)*Power(xii,9.0) + 
         232373746276140268800.0*Power(rij,10.0)*Power(xii,10.0) + 
         40613247709652217600.0*Power(rij,11.0)*Power(xii,11.0) + 
         6465950980961472000.0*Power(rij,12.0)*Power(xii,12.0) + 
         943771901519462400.0*Power(rij,13.0)*Power(xii,13.0) + 
         126929425622630400.0*Power(rij,14.0)*Power(xii,14.0) + 
         15790263474585600.0*Power(rij,15.0)*Power(xii,15.0) + 
         1821953477836800.0*Power(rij,16.0)*Power(xii,16.0) + 
         195294359715840.0*Power(rij,17.0)*Power(xii,17.0) + 
         19450048020480.0*Power(rij,18.0)*Power(xii,18.0) + 
         1796674682880.0*Power(rij,19.0)*Power(xii,19.0) + 
         153204817920.0*Power(rij,20.0)*Power(xii,20.0) + 
         11938037760.0*Power(rij,21.0)*Power(xii,21.0) + 
         832045056.0*Power(rij,22.0)*Power(xii,22.0) + 
         49283072.0*Power(rij,23.0)*Power(xii,23.0) + 
         2097152.0*Power(rij,24.0)*Power(xii,24.0))/
      (9.306726025998591e23*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  }
  else {
    S =     (1403325.0*Power(E,2.0*rij*(xii + xij))*Power(Power(xii,2.0) - Power(xij,2.0),23.0) + 
        Power(E,2.0*rij*xij)*Power(xij,14.0)*
         (-996.0*Power(rij,11.0)*Power(xii,43.0) - 12.0*Power(rij,12.0)*Power(xii,44.0) + 
           1403325.0*Power(xij,32.0) + 2806650.0*rij*xii*Power(xij,32.0) - 
           88.0*Power(rij,10.0)*Power(xii,42.0)*(456.0 + Power(rij,2.0)*Power(xij,2.0)) + 
           935550.0*rij*Power(xii,3.0)*Power(xij,30.0)*
            (-69.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) + 
           1403325.0*Power(xii,2.0)*Power(xij,30.0)*
            (-23.0 + 2.0*Power(rij,2.0)*Power(xij,2.0)) - 
           44.0*Power(rij,9.0)*Power(xii,41.0)*
            (23460.0 + 257.0*Power(rij,2.0)*Power(xij,2.0)) + 
           187110.0*rij*Power(xii,5.0)*Power(xij,28.0)*
            (3795.0 - 230.0*Power(rij,2.0)*Power(xij,2.0) + 
              2.0*Power(rij,4.0)*Power(xij,4.0)) + 
           467775.0*Power(xii,4.0)*Power(xij,28.0)*
            (759.0 - 138.0*Power(rij,2.0)*Power(xij,2.0) + 2.0*Power(rij,4.0)*Power(xij,4.0)) \
    + 44.0*Power(rij,8.0)*Power(xii,40.0)*
            (-426870.0 - 14241.0*Power(rij,2.0)*Power(xij,2.0) + 
              22.0*Power(rij,4.0)*Power(xij,4.0)) + 
           220.0*Power(rij,7.0)*Power(xii,39.0)*
            (-1151172.0 - 94377.0*Power(rij,2.0)*Power(xij,2.0) + 
              347.0*Power(rij,4.0)*Power(xij,4.0)) + 
           8910.0*rij*Power(xii,7.0)*Power(xij,26.0)*
            (-557865.0 + 53130.0*Power(rij,2.0)*Power(xij,2.0) - 
              966.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) + 
           31185.0*Power(xii,6.0)*Power(xij,26.0)*
            (-79695.0 + 22770.0*Power(rij,2.0)*Power(xij,2.0) - 
              690.0*Power(rij,4.0)*Power(xij,4.0) + 4.0*Power(rij,6.0)*Power(xij,6.0)) - 
           110.0*Power(rij,6.0)*Power(xii,38.0)*
            (23442048.0 + 4242321.0*Power(rij,2.0)*Power(xij,2.0) - 
              22038.0*Power(rij,4.0)*Power(xij,4.0) + 20.0*Power(rij,6.0)*Power(xij,6.0)) - 
           44.0*Power(rij,5.0)*Power(xii,37.0)*
            (450526860.0 + 171521550.0*Power(rij,2.0)*Power(xij,2.0) - 
              702615.0*Power(rij,4.0)*Power(xij,4.0) + 809.0*Power(rij,6.0)*Power(xij,6.0)) \
    + 990.0*rij*Power(xii,9.0)*Power(xij,24.0)*
            (25103925.0 - 3347190.0*Power(rij,2.0)*Power(xij,2.0) + 
              95634.0*Power(rij,4.0)*Power(xij,4.0) - 
              828.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) + 
           4455.0*Power(xii,8.0)*Power(xij,24.0)*
            (2789325.0 - 1115730.0*Power(rij,2.0)*Power(xij,2.0) + 
              53130.0*Power(rij,4.0)*Power(xij,4.0) - 
              644.0*Power(rij,6.0)*Power(xij,6.0) + 2.0*Power(rij,8.0)*Power(xij,8.0)) - 
           22.0*Power(rij,4.0)*Power(xii,36.0)*
            (5137105050.0 + 4097941470.0*Power(rij,2.0)*Power(xij,2.0) + 
              14337945.0*Power(rij,4.0)*Power(xij,4.0) - 
              245166.0*Power(rij,6.0)*Power(xij,6.0) + 124.0*Power(rij,8.0)*Power(xij,8.0)) \
    - 44.0*Power(rij,3.0)*Power(xii,35.0)*
            (10425301425.0 + 18160961805.0*Power(rij,2.0)*Power(xij,2.0) + 
              471245850.0*Power(rij,4.0)*Power(xij,4.0) - 
              7008675.0*Power(rij,6.0)*Power(xij,6.0) + 
              15533.0*Power(rij,8.0)*Power(xij,8.0)) - 
           2.0*Power(rij,2.0)*Power(xii,34.0)*
            (611617683600.0 + 2588132183175.0*Power(rij,2.0)*Power(xij,2.0) + 
              215628370650.0*Power(rij,4.0)*Power(xij,4.0) - 
              3937764105.0*Power(rij,6.0)*Power(xij,6.0) + 
              18963318.0*Power(rij,8.0)*Power(xij,8.0) - 
              11236.0*Power(rij,10.0)*Power(xij,10.0)) + 
           18.0*rij*Power(xii,11.0)*Power(xij,22.0)*
            (-5246720325.0 + 920477250.0*Power(rij,2.0)*Power(xij,2.0) - 
              36819090.0*Power(rij,4.0)*Power(xij,4.0) + 
              500940.0*Power(rij,6.0)*Power(xij,6.0) - 
              2530.0*Power(rij,8.0)*Power(xij,8.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) + 
           99.0*Power(xii,10.0)*Power(xij,22.0)*
            (-476974575.0 + 251039250.0*Power(rij,2.0)*Power(xij,2.0) - 
              16735950.0*Power(rij,4.0)*Power(xij,4.0) + 
              318780.0*Power(rij,6.0)*Power(xij,6.0) - 
              2070.0*Power(rij,8.0)*Power(xij,8.0) + 4.0*Power(rij,10.0)*Power(xij,10.0)) - 
           22.0*rij*Power(xii,13.0)*Power(xij,20.0)*
            (-12878313525.0 + 3695139000.0*Power(rij,2.0)*Power(xij,2.0) - 
              71262450.0*Power(rij,4.0)*Power(xij,4.0) + 
              5017140.0*Power(rij,6.0)*Power(xij,6.0) - 
              930.0*Power(rij,8.0)*Power(xij,8.0) + 158.0*Power(rij,10.0)*Power(xij,10.0)) \
    - 990.0*rij*Power(xii,29.0)*Power(xij,4.0)*
            (1620022028625.0 + 769716305820.0*Power(rij,2.0)*Power(xij,2.0) - 
              31518809406.0*Power(rij,4.0)*Power(xij,4.0) + 
              492367788.0*Power(rij,6.0)*Power(xij,6.0) - 
              2173822.0*Power(rij,8.0)*Power(xij,8.0) + 
              322.0*Power(rij,10.0)*Power(xij,10.0)) - 
           22.0*rij*Power(xii,15.0)*Power(xij,18.0)*
            (113414642775.0 + 149162958000.0*Power(rij,2.0)*Power(xij,2.0) + 
              11413095750.0*Power(rij,4.0)*Power(xij,4.0) + 
              185668740.0*Power(rij,6.0)*Power(xij,6.0) + 
              1327830.0*Power(rij,8.0)*Power(xij,8.0) + 
              446.0*Power(rij,10.0)*Power(xij,10.0)) - 
           110.0*rij*Power(xii,19.0)*Power(xij,14.0)*
            (14601048860475.0 + 4483288448100.0*Power(rij,2.0)*Power(xij,2.0) + 
              57946422870.0*Power(rij,4.0)*Power(xij,4.0) - 
              832082868.0*Power(rij,6.0)*Power(xij,6.0) - 
              3260778.0*Power(rij,8.0)*Power(xij,8.0) + 
              4006.0*Power(rij,10.0)*Power(xij,10.0)) + 
           22.0*rij*Power(xii,17.0)*Power(xij,16.0)*
            (-4973277615075.0 - 3692030989500.0*Power(rij,2.0)*Power(xij,2.0) - 
              151636618350.0*Power(rij,4.0)*Power(xij,4.0) - 
              1305433260.0*Power(rij,6.0)*Power(xij,6.0) + 
              2206410.0*Power(rij,8.0)*Power(xij,8.0) + 
              7282.0*Power(rij,10.0)*Power(xij,10.0)) - 
           462.0*rij*Power(xii,25.0)*Power(xij,8.0)*
            (40495013164125.0 - 3518150807250.0*Power(rij,2.0)*Power(xij,2.0) + 
              108840713850.0*Power(rij,4.0)*Power(xij,4.0) - 
              649652940.0*Power(rij,6.0)*Power(xij,6.0) - 
              3316830.0*Power(rij,8.0)*Power(xij,8.0) + 
              8878.0*Power(rij,10.0)*Power(xij,10.0)) + 
           198.0*rij*Power(xii,27.0)*Power(xij,6.0)*
            (-42466116317625.0 - 406603622250.0*Power(rij,2.0)*Power(xij,2.0) - 
              27736200450.0*Power(rij,4.0)*Power(xij,4.0) + 
              2056016460.0*Power(rij,6.0)*Power(xij,6.0) - 
              17679410.0*Power(rij,8.0)*Power(xij,8.0) + 
              18354.0*Power(rij,10.0)*Power(xij,10.0)) + 
           2.0*rij*Power(xii,21.0)*Power(xij,12.0)*
            (-4202540019025425.0 - 360777507725250.0*Power(rij,2.0)*Power(xij,2.0) + 
              11255464316250.0*Power(rij,4.0)*Power(xij,4.0) + 
              55487803140.0*Power(rij,6.0)*Power(xij,6.0) - 
              628960530.0*Power(rij,8.0)*Power(xij,8.0) + 
              36974.0*Power(rij,10.0)*Power(xij,10.0)) - 
           2475.0*Power(xii,28.0)*Power(xij,4.0)*
            (324004405725.0 + 1510859435778.0*Power(rij,2.0)*Power(xij,2.0) - 
              65979118674.0*Power(rij,4.0)*Power(xij,4.0) + 
              1845131148.0*Power(rij,6.0)*Power(xij,6.0) - 
              20728198.0*Power(rij,8.0)*Power(xij,8.0) + 
              49588.0*Power(rij,10.0)*Power(xij,10.0)) + 
           22.0*rij*Power(xii,23.0)*Power(xij,10.0)*
            (-850567767797475.0 + 36188347907250.0*Power(rij,2.0)*Power(xij,2.0) + 
              847439173350.0*Power(rij,4.0)*Power(xij,4.0) - 
              23276359260.0*Power(rij,6.0)*Power(xij,6.0) + 
              45088050.0*Power(rij,8.0)*Power(xij,8.0) + 
              84778.0*Power(rij,10.0)*Power(xij,10.0)) - 
           22.0*rij*Power(xii,31.0)*Power(xij,2.0)*
            (5035829423625.0 + 11622729957300.0*Power(rij,2.0)*Power(xij,2.0) + 
              16351979490.0*Power(rij,4.0)*Power(xij,4.0) - 
              1720091700.0*Power(rij,6.0)*Power(xij,6.0) - 
              14708670.0*Power(rij,8.0)*Power(xij,8.0) + 
              103546.0*Power(rij,10.0)*Power(xij,10.0)) + 
           4.0*rij*Power(xii,33.0)*(-451763061750.0 - 
              5888319267600.0*Power(rij,2.0)*Power(xij,2.0) - 
              1360794117375.0*Power(rij,4.0)*Power(xij,4.0) + 
              28346942250.0*Power(rij,6.0)*Power(xij,6.0) - 
              234258915.0*Power(rij,8.0)*Power(xij,8.0) + 
              518489.0*Power(rij,10.0)*Power(xij,10.0)) + 
           3.0*Power(xii,12.0)*Power(xij,20.0)*
            (47220482925.0 - 31480321950.0*Power(rij,2.0)*Power(xij,2.0) + 
              2761431750.0*Power(rij,4.0)*Power(xij,4.0) - 
              73638180.0*Power(rij,6.0)*Power(xij,6.0) + 
              751410.0*Power(rij,8.0)*Power(xij,8.0) - 
              3036.0*Power(rij,10.0)*Power(xij,10.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 11.0*Power(xii,14.0)*Power(xij,18.0)*
            (-113414642775.0 + 78830888850.0*Power(rij,2.0)*Power(xij,2.0) + 
              8804234250.0*Power(rij,4.0)*Power(xij,4.0) + 
              576582300.0*Power(rij,6.0)*Power(xij,6.0) + 
              5009850.0*Power(rij,8.0)*Power(xij,8.0) + 
              39684.0*Power(rij,10.0)*Power(xij,10.0) + 8.0*Power(rij,12.0)*Power(xij,12.0)\
    ) + 55.0*Power(xii,18.0)*Power(xij,14.0)*
            (-14601048860475.0 + 11014160607450.0*Power(rij,2.0)*Power(xij,2.0) + 
              899524492650.0*Power(rij,4.0)*Power(xij,4.0) + 
              8911784700.0*Power(rij,6.0)*Power(xij,6.0) - 
              53799750.0*Power(rij,8.0)*Power(xij,8.0) - 
              189276.0*Power(rij,10.0)*Power(xij,10.0) + 
              40.0*Power(rij,12.0)*Power(xij,12.0)) - 
           11.0*Power(xii,16.0)*Power(xij,16.0)*
            (4973277615075.0 - 3535546162050.0*Power(rij,2.0)*Power(xij,2.0) - 
              536379590250.0*Power(rij,4.0)*Power(xij,4.0) - 
              13554229500.0*Power(rij,6.0)*Power(xij,6.0) - 
              104593770.0*Power(rij,8.0)*Power(xij,8.0) + 
              3228.0*Power(rij,10.0)*Power(xij,10.0) + 88.0*Power(rij,12.0)*Power(xij,12.0)\
    ) - 231.0*Power(xii,26.0)*Power(xij,6.0)*
            (18199764136125.0 + 21741060926250.0*Power(rij,2.0)*Power(xij,2.0) - 
              1022330706750.0*Power(rij,4.0)*Power(xij,4.0) + 
              7759539900.0*Power(rij,6.0)*Power(xij,6.0) + 
              118181130.0*Power(rij,8.0)*Power(xij,8.0) - 
              676476.0*Power(rij,10.0)*Power(xij,10.0) + 
              184.0*Power(rij,12.0)*Power(xij,12.0)) + 
           11.0*Power(xii,20.0)*Power(xij,12.0)*
            (-382049092638675.0 + 278436232062450.0*Power(rij,2.0)*Power(xij,2.0) + 
              6550223699250.0*Power(rij,4.0)*Power(xij,4.0) - 
              160205642100.0*Power(rij,6.0)*Power(xij,6.0) - 
              685586070.0*Power(rij,8.0)*Power(xij,8.0) + 
              2934804.0*Power(rij,10.0)*Power(xij,10.0) + 
              248.0*Power(rij,12.0)*Power(xij,12.0)) + 
           33.0*Power(xii,30.0)*Power(xij,2.0)*
            (-1678609807875.0 - 26928208102950.0*Power(rij,2.0)*Power(xij,2.0) - 
              2249376232950.0*Power(rij,4.0)*Power(xij,4.0) + 
              75178632060.0*Power(rij,6.0)*Power(xij,6.0) - 
              585984390.0*Power(rij,8.0)*Power(xij,8.0) + 
              91908.0*Power(rij,10.0)*Power(xij,10.0) + 
              1288.0*Power(rij,12.0)*Power(xij,12.0)) - 
           22.0*Power(xii,32.0)*(41069369250.0 + 
              3154641699825.0*Power(rij,2.0)*Power(xij,2.0) + 
              2079992878425.0*Power(rij,4.0)*Power(xij,4.0) - 
              38081936970.0*Power(rij,6.0)*Power(xij,6.0) + 
              451881945.0*Power(rij,8.0)*Power(xij,8.0) - 
              2730774.0*Power(rij,10.0)*Power(xij,10.0) + 
              2116.0*Power(rij,12.0)*Power(xij,12.0)) + 
           11.0*Power(xii,24.0)*Power(xij,8.0)*
            (-850395276446625.0 + 66465982229850.0*Power(rij,2.0)*Power(xij,2.0) - 
              19593192890250.0*Power(rij,4.0)*Power(xij,4.0) + 
              554675026500.0*Power(rij,6.0)*Power(xij,6.0) - 
              2106396810.0*Power(rij,8.0)*Power(xij,8.0) - 
              6122508.0*Power(rij,10.0)*Power(xij,10.0) + 
              4232.0*Power(rij,12.0)*Power(xij,12.0)) - 
           Power(xii,22.0)*Power(xij,10.0)*
            (9356245445772225.0 - 5281430969811150.0*Power(rij,2.0)*Power(xij,2.0) + 
              186183223899750.0*Power(rij,4.0)*Power(xij,4.0) + 
              1420112301300.0*Power(rij,6.0)*Power(xij,6.0) - 
              30664882710.0*Power(rij,8.0)*Power(xij,8.0) + 
              20787492.0*Power(rij,10.0)*Power(xij,10.0) + 
              22472.0*Power(rij,12.0)*Power(xij,12.0))) + 
        Power(E,2.0*rij*xii)*Power(xii,14.0)*
         (-302841.0*Power(xii,16.0)*Power(xij,16.0)*
            (-180642825.0 - 361285650.0*rij*xij + 
              2000319750.0*Power(rij,2.0)*Power(xij,2.0) - 
              1628451000.0*Power(rij,3.0)*Power(xij,3.0) + 
              237921750.0*Power(rij,4.0)*Power(xij,4.0) + 
              74332500.0*Power(rij,5.0)*Power(xij,5.0) - 
              4689300.0*Power(rij,6.0)*Power(xij,6.0) - 
              1690920.0*Power(rij,7.0)*Power(xij,7.0) - 
              76510.0*Power(rij,8.0)*Power(xij,8.0) + 
              5060.0*Power(rij,9.0)*Power(xij,9.0) + 
              516.0*Power(rij,10.0)*Power(xij,10.0) + 12.0*Power(rij,11.0)*Power(xij,11.0)) \
    + Power(xii,10.0)*Power(xij,22.0)*(9356245445772225.0 + 
              18712490891544450.0*rij*xij - 
              731125804528350.0*Power(rij,2.0)*Power(xij,2.0) - 
              1625385672949500.0*Power(rij,3.0)*Power(xij,3.0) - 
              236158393259250.0*Power(rij,4.0)*Power(xij,4.0) + 
              5491767689100.0*Power(rij,5.0)*Power(xij,5.0) + 
              4566699591300.0*Power(rij,6.0)*Power(xij,6.0) + 
              487444110120.0*Power(rij,7.0)*Power(xij,7.0) + 
              19337484870.0*Power(rij,8.0)*Power(xij,8.0) - 
              323590740.0*Power(rij,9.0)*Power(xij,9.0) - 
              60077028.0*Power(rij,10.0)*Power(xij,10.0) - 
              2073956.0*Power(rij,11.0)*Power(xij,11.0) - 
              22472.0*Power(rij,12.0)*Power(xij,12.0)) + 
           11.0*Power(xii,30.0)*Power(xij,2.0)*
            (2934225.0 + 5868450.0*rij*xij + 5868450.0*Power(rij,2.0)*Power(xij,2.0) + 
              3912300.0*Power(rij,3.0)*Power(xij,3.0) + 
              1956150.0*Power(rij,4.0)*Power(xij,4.0) + 
              782460.0*Power(rij,5.0)*Power(xij,5.0) + 
              260820.0*Power(rij,6.0)*Power(xij,6.0) + 
              74520.0*Power(rij,7.0)*Power(xij,7.0) + 
              18630.0*Power(rij,8.0)*Power(xij,8.0) + 
              4140.0*Power(rij,9.0)*Power(xij,9.0) + 
              828.0*Power(rij,10.0)*Power(xij,10.0) + 
              316.0*Power(rij,11.0)*Power(xij,11.0) - 8.0*Power(rij,12.0)*Power(xij,12.0)) + 
           12.0*Power(xij,32.0)*(75293843625.0 + 150587687250.0*rij*xij + 
              101936280600.0*Power(rij,2.0)*Power(xij,2.0) + 
              38226105225.0*Power(rij,3.0)*Power(xij,3.0) + 
              9418025925.0*Power(rij,4.0)*Power(xij,4.0) + 
              1651931820.0*Power(rij,5.0)*Power(xij,5.0) + 
              214885440.0*Power(rij,6.0)*Power(xij,6.0) + 
              21104820.0*Power(rij,7.0)*Power(xij,7.0) + 
              1565190.0*Power(rij,8.0)*Power(xij,8.0) + 
              86020.0*Power(rij,9.0)*Power(xij,9.0) + 
              3344.0*Power(rij,10.0)*Power(xij,10.0) + 
              83.0*Power(rij,11.0)*Power(xij,11.0) + Power(rij,12.0)*Power(xij,12.0)) - 
           3.0*Power(xii,32.0)*(467775.0 + 935550.0*rij*xij + 
              935550.0*Power(rij,2.0)*Power(xij,2.0) + 
              623700.0*Power(rij,3.0)*Power(xij,3.0) + 
              311850.0*Power(rij,4.0)*Power(xij,4.0) + 
              124740.0*Power(rij,5.0)*Power(xij,5.0) + 
              41580.0*Power(rij,6.0)*Power(xij,6.0) + 
              11880.0*Power(rij,7.0)*Power(xij,7.0) + 
              2970.0*Power(rij,8.0)*Power(xij,8.0) + 660.0*Power(rij,9.0)*Power(xij,9.0) + 
              132.0*Power(rij,10.0)*Power(xij,10.0) + 
              24.0*Power(rij,11.0)*Power(xij,11.0) + 4.0*Power(rij,12.0)*Power(xij,12.0)) - 
           5313.0*Power(xii,14.0)*Power(xij,18.0)*
            (-151149574125.0 - 302299148250.0*rij*xij + 
              576472530150.0*Power(rij,2.0)*Power(xij,2.0) - 
              135809338500.0*Power(rij,3.0)*Power(xij,3.0) - 
              35042955750.0*Power(rij,4.0)*Power(xij,4.0) + 
              3509064900.0*Power(rij,5.0)*Power(xij,5.0) + 
              1148395500.0*Power(rij,6.0)*Power(xij,6.0) + 
              56491560.0*Power(rij,7.0)*Power(xij,7.0) - 
              5138310.0*Power(rij,8.0)*Power(xij,8.0) - 
              658860.0*Power(rij,9.0)*Power(xij,9.0) - 
              23100.0*Power(rij,10.0)*Power(xij,10.0) - 
              60.0*Power(rij,11.0)*Power(xij,11.0) + 8.0*Power(rij,12.0)*Power(xij,12.0)) + 
           5313.0*Power(xii,18.0)*Power(xij,14.0)*
            (234812925.0 + 469625850.0*rij*xij - 
              7319971350.0*Power(rij,2.0)*Power(xij,2.0) + 
              15287913000.0*Power(rij,3.0)*Power(xij,3.0) - 
              9311847750.0*Power(rij,4.0)*Power(xij,4.0) + 
              1199718900.0*Power(rij,5.0)*Power(xij,5.0) + 
              331688700.0*Power(rij,6.0)*Power(xij,6.0) - 
              20887560.0*Power(rij,7.0)*Power(xij,7.0) - 
              5771670.0*Power(rij,8.0)*Power(xij,8.0) - 
              186700.0*Power(rij,9.0)*Power(xij,9.0) + 
              12676.0*Power(rij,10.0)*Power(xij,10.0) + 
              772.0*Power(rij,11.0)*Power(xij,11.0) + 8.0*Power(rij,12.0)*Power(xij,12.0)) + 
           11.0*Power(xii,2.0)*Power(xij,30.0)*
            (5035829423625.0 + 10071658847250.0*rij*xij + 
              6309283399650.0*Power(rij,2.0)*Power(xij,2.0) + 
              2141207006400.0*Power(rij,3.0)*Power(xij,3.0) + 
              470569487850.0*Power(rij,4.0)*Power(xij,4.0) + 
              72643847220.0*Power(rij,5.0)*Power(xij,5.0) + 
              8195882940.0*Power(rij,6.0)*Power(xij,6.0) + 
              686086200.0*Power(rij,7.0)*Power(xij,7.0) + 
              42423210.0*Power(rij,8.0)*Power(xij,8.0) + 
              1887540.0*Power(rij,9.0)*Power(xij,9.0) + 
              56964.0*Power(rij,10.0)*Power(xij,10.0) + 
              1028.0*Power(rij,11.0)*Power(xij,11.0) + 8.0*Power(rij,12.0)*Power(xij,12.0)) \
    + 11.0*Power(xii,28.0)*Power(xij,4.0)*
            (-32276475.0 - 64552950.0*rij*xij - 
              64552950.0*Power(rij,2.0)*Power(xij,2.0) - 
              43035300.0*Power(rij,3.0)*Power(xij,3.0) - 
              21517650.0*Power(rij,4.0)*Power(xij,4.0) - 
              8607060.0*Power(rij,5.0)*Power(xij,5.0) - 
              2869020.0*Power(rij,6.0)*Power(xij,6.0) - 
              819720.0*Power(rij,7.0)*Power(xij,7.0) - 
              204930.0*Power(rij,8.0)*Power(xij,8.0) - 
              1860.0*Power(rij,9.0)*Power(xij,9.0) - 
              39684.0*Power(rij,10.0)*Power(xij,10.0) + 
              892.0*Power(rij,11.0)*Power(xij,11.0) + 88.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 11.0*Power(xii,4.0)*Power(xij,28.0)*
            (-72900991288125.0 - 145801982576250.0*rij*xij - 
              80784624308850.0*Power(rij,2.0)*Power(xij,2.0) - 
              23245459914600.0*Power(rij,3.0)*Power(xij,3.0) - 
              4159985756850.0*Power(rij,4.0)*Power(xij,4.0) - 
              494834224500.0*Power(rij,5.0)*Power(xij,5.0) - 
              39205158300.0*Power(rij,6.0)*Power(xij,6.0) - 
              1884983400.0*Power(rij,7.0)*Power(xij,7.0) - 
              28675890.0*Power(rij,8.0)*Power(xij,8.0) + 
              2810460.0*Power(rij,9.0)*Power(xij,9.0) + 
              220380.0*Power(rij,10.0)*Power(xij,10.0) + 
              6940.0*Power(rij,11.0)*Power(xij,11.0) + 88.0*Power(rij,12.0)*Power(xij,12.0)) \
    - 253.0*Power(xii,20.0)*Power(xij,12.0)*
            (559926675.0 + 1119853350.0*rij*xij + 
              3427429950.0*Power(rij,2.0)*Power(xij,2.0) - 
              12970692000.0*Power(rij,3.0)*Power(xij,3.0) + 
              23320851750.0*Power(rij,4.0)*Power(xij,4.0) - 
              13185792900.0*Power(rij,5.0)*Power(xij,5.0) + 
              1937344500.0*Power(rij,6.0)*Power(xij,6.0) + 
              361775160.0*Power(rij,7.0)*Power(xij,7.0) - 
              29808090.0*Power(rij,8.0)*Power(xij,8.0) - 
              4972020.0*Power(rij,9.0)*Power(xij,9.0) - 
              82164.0*Power(rij,10.0)*Power(xij,10.0) + 
              7372.0*Power(rij,11.0)*Power(xij,11.0) + 184.0*Power(rij,12.0)*Power(xij,12.0)\
    ) + 253.0*Power(xii,12.0)*Power(xij,20.0)*
            (16610830114725.0 + 33221660229450.0*rij*xij - 
              20875221224550.0*Power(rij,2.0)*Power(xij,2.0) - 
              3146812861500.0*Power(rij,3.0)*Power(xij,3.0) + 
              851877951750.0*Power(rij,4.0)*Power(xij,4.0) + 
              198752607900.0*Power(rij,5.0)*Power(xij,5.0) + 
              7084797300.0*Power(rij,6.0)*Power(xij,6.0) - 
              1609056360.0*Power(rij,7.0)*Power(xij,7.0) - 
              202775850.0*Power(rij,8.0)*Power(xij,8.0) - 
              8506260.0*Power(rij,9.0)*Power(xij,9.0) - 
              11988.0*Power(rij,10.0)*Power(xij,10.0) + 
              9004.0*Power(rij,11.0)*Power(xij,11.0) + 184.0*Power(rij,12.0)*Power(xij,12.0)\
    ) + 11.0*Power(xii,6.0)*Power(xij,26.0)*
            (382195046858625.0 + 764390093717250.0*rij*xij + 
              339943373050050.0*Power(rij,2.0)*Power(xij,2.0) + 
              69274467523800.0*Power(rij,3.0)*Power(xij,3.0) + 
              6748128698850.0*Power(rij,4.0)*Power(xij,4.0) + 
              32703958980.0*Power(rij,5.0)*Power(xij,5.0) - 
              76163873940.0*Power(rij,6.0)*Power(xij,6.0) - 
              10307979000.0*Power(rij,7.0)*Power(xij,7.0) - 
              715957110.0*Power(rij,8.0)*Power(xij,8.0) - 
              28034700.0*Power(rij,9.0)*Power(xij,9.0) - 
              490332.0*Power(rij,10.0)*Power(xij,10.0) + 
              3236.0*Power(rij,11.0)*Power(xij,11.0) + 200.0*Power(rij,12.0)*Power(xij,12.0)\
    ) - 11.0*Power(xii,26.0)*Power(xij,6.0)*
            (-225935325.0 - 451870650.0*rij*xij - 
              451870650.0*Power(rij,2.0)*Power(xij,2.0) - 
              301247100.0*Power(rij,3.0)*Power(xij,3.0) - 
              150623550.0*Power(rij,4.0)*Power(xij,4.0) - 
              60249420.0*Power(rij,5.0)*Power(xij,5.0) - 
              20083140.0*Power(rij,6.0)*Power(xij,6.0) - 
              10034280.0*Power(rij,7.0)*Power(xij,7.0) + 
              5009850.0*Power(rij,8.0)*Power(xij,8.0) - 
              2655660.0*Power(rij,9.0)*Power(xij,9.0) - 
              3228.0*Power(rij,10.0)*Power(xij,10.0) + 
              14564.0*Power(rij,11.0)*Power(xij,11.0) + 
              200.0*Power(rij,12.0)*Power(xij,12.0)) - 
           11.0*Power(xii,24.0)*Power(xij,8.0)*
            (1129676625.0 + 2259353250.0*rij*xij + 
              2259353250.0*Power(rij,2.0)*Power(xij,2.0) + 
              1506235500.0*Power(rij,3.0)*Power(xij,3.0) + 
              753117750.0*Power(rij,4.0)*Power(xij,4.0) + 
              142524900.0*Power(rij,5.0)*Power(xij,5.0) + 
              576582300.0*Power(rij,6.0)*Power(xij,6.0) - 
              371337480.0*Power(rij,7.0)*Power(xij,7.0) + 
              104593770.0*Power(rij,8.0)*Power(xij,8.0) + 
              4412820.0*Power(rij,9.0)*Power(xij,9.0) - 
              946380.0*Power(rij,10.0)*Power(xij,10.0) - 
              40060.0*Power(rij,11.0)*Power(xij,11.0) + 
              248.0*Power(rij,12.0)*Power(xij,12.0)) + 
           11.0*Power(xii,8.0)*Power(xij,24.0)*
            (850395276446625.0 + 1700790552893250.0*rij*xij + 
              456562279451250.0*Power(rij,2.0)*Power(xij,2.0) + 
              7318865200500.0*Power(rij,3.0)*Power(xij,3.0) - 
              14845301701650.0*Power(rij,4.0)*Power(xij,4.0) - 
              2836692846540.0*Power(rij,5.0)*Power(xij,5.0) - 
              225535896180.0*Power(rij,6.0)*Power(xij,6.0) - 
              3440183400.0*Power(rij,7.0)*Power(xij,7.0) + 
              903763890.0*Power(rij,8.0)*Power(xij,8.0) + 
              85185060.0*Power(rij,9.0)*Power(xij,9.0) + 
              3447876.0*Power(rij,10.0)*Power(xij,10.0) + 
              62132.0*Power(rij,11.0)*Power(xij,11.0) + 
              248.0*Power(rij,12.0)*Power(xij,12.0)) + 
           Power(xii,22.0)*Power(xij,10.0)*
            (47220482925.0 + 94440965850.0*rij*xij + 
              94440965850.0*Power(rij,2.0)*Power(xij,2.0) + 
              81293058000.0*Power(rij,3.0)*Power(xij,3.0) - 
              96846576750.0*Power(rij,4.0)*Power(xij,4.0) + 
              251088106500.0*Power(rij,5.0)*Power(xij,5.0) - 
              149096524500.0*Power(rij,6.0)*Power(xij,6.0) + 
              28719531720.0*Power(rij,7.0)*Power(xij,7.0) + 
              2958986250.0*Power(rij,8.0)*Power(xij,8.0) - 
              358685580.0*Power(rij,9.0)*Power(xij,9.0) - 
              32282844.0*Power(rij,10.0)*Power(xij,10.0) - 
              73948.0*Power(rij,11.0)*Power(xij,11.0) + 
              22472.0*Power(rij,12.0)*Power(xij,12.0))))/
      (1.403325e6*Power(E,2.0*rij*(xii + xij))*Power(rij,2.0)*
        Power(Power(xii,2.0) - Power(xij,2.0),23.0))
    ;
  }
  return S;
}

double DSlater_6S_1S(double rij,double xii,double xij)
{
  return DSlater_1S_6S(rij,xij,xii);
}

double DSlater_6S_2S(double rij,double xii,double xij)
{
  return DSlater_2S_6S(rij,xij,xii);
}

double DSlater_6S_3S(double rij,double xii,double xij)
{
  return DSlater_3S_6S(rij,xij,xii);
}

double DSlater_6S_4S(double rij,double xii,double xij)
{
  return DSlater_4S_6S(rij,xij,xii);
}

double DSlater_6S_5S(double rij,double xii,double xij)
{
  return DSlater_5S_6S(rij,xij,xii);
}

double Nuclear_1S(double rij,double xii)
{
  double S;

  S = 
  1.0/rij - (1.0 + rij*xii)/(Power(E,2.0*rij*xii)*rij)
    ;
  return S;
}

double Nuclear_2S(double rij,double xii)
{
  double S;

  S = 
  1.0/rij - (6.0 + 9.0*rij*xii + 6.0*Power(rij,2.0)*Power(xii,2.0) + 
       2.0*Power(rij,3.0)*Power(xii,3.0))/(6.*Power(E,2.0*rij*xii)*rij)
    ;
  return S;
}

double Nuclear_3S(double rij,double xii)
{
  double S;

  S = 
  1.0/rij - (45.0 + 75.0*rij*xii + 60.0*Power(rij,2.0)*Power(xii,2.0) + 
       30.0*Power(rij,3.0)*Power(xii,3.0) + 10.0*Power(rij,4.0)*Power(xii,4.0) + 
       2.0*Power(rij,5.0)*Power(xii,5.0))/(45.*Power(E,2.0*rij*xii)*rij)
    ;
  return S;
}

double Nuclear_4S(double rij,double xii)
{
  double S;

  S = 
  1.0/rij - (1260.0 + 2205.0*rij*xii + 1890.0*Power(rij,2.0)*Power(xii,2.0) + 
       1050.0*Power(rij,3.0)*Power(xii,3.0) + 420.0*Power(rij,4.0)*Power(xii,4.0) + 
       126.0*Power(rij,5.0)*Power(xii,5.0) + 28.0*Power(rij,6.0)*Power(xii,6.0) + 
       4.0*Power(rij,7.0)*Power(xii,7.0))/(1260.*Power(E,2.0*rij*xii)*rij)
    ;
  return S;
}

double Nuclear_5S(double rij,double xii)
{
  double S;

  S = 
  1.0/rij - (14175.0 + 25515.0*rij*xii + 22680.0*Power(rij,2.0)*Power(xii,2.0) + 
       13230.0*Power(rij,3.0)*Power(xii,3.0) + 5670.0*Power(rij,4.0)*Power(xii,4.0) + 
       1890.0*Power(rij,5.0)*Power(xii,5.0) + 504.0*Power(rij,6.0)*Power(xii,6.0) + 
       108.0*Power(rij,7.0)*Power(xii,7.0) + 18.0*Power(rij,8.0)*Power(xii,8.0) + 
       2.0*Power(rij,9.0)*Power(xii,9.0))/(14175.*Power(E,2.0*rij*xii)*rij)
    ;
  return S;
}

double Nuclear_6S(double rij,double xii)
{
  double S;

  S = 
  1.0/rij - (935550.0 + 1715175.0*rij*xii + 1559250.0*Power(rij,2.0)*Power(xii,2.0) + 
       935550.0*Power(rij,3.0)*Power(xii,3.0) + 415800.0*Power(rij,4.0)*Power(xii,4.0) + 
       145530.0*Power(rij,5.0)*Power(xii,5.0) + 41580.0*Power(rij,6.0)*Power(xii,6.0) + 
       9900.0*Power(rij,7.0)*Power(xii,7.0) + 1980.0*Power(rij,8.0)*Power(xii,8.0) + 
       330.0*Power(rij,9.0)*Power(xii,9.0) + 44.0*Power(rij,10.0)*Power(xii,10.0) + 
       4.0*Power(rij,11.0)*Power(xii,11.0))/(935550.*Power(E,2.0*rij*xii)*rij)
    ;
  return S;
}

double DNuclear_1S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2.0) - (1.0 + 2.0*rij*xii + 2.0*Power(rij,2.0)*Power(xii,2.0))/
     (Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  return S;
}

double DNuclear_2S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2.0) - (3.0 + 6.0*rij*xii + 6.0*Power(rij,2.0)*Power(xii,2.0) + 
       4.0*Power(rij,3.0)*Power(xii,3.0) + 2.0*Power(rij,4.0)*Power(xii,4.0))/
     (3.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  return S;
}

double DNuclear_3S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2.0) - (45.0 + 90.0*rij*xii + 90.0*Power(rij,2.0)*Power(xii,2.0) + 
       60.0*Power(rij,3.0)*Power(xii,3.0) + 30.0*Power(rij,4.0)*Power(xii,4.0) + 
       12.0*Power(rij,5.0)*Power(xii,5.0) + 4.0*Power(rij,6.0)*Power(xii,6.0))/
     (45.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  return S;
}

double DNuclear_4S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2.0) - (315.0 + 630.0*rij*xii + 630.0*Power(rij,2.0)*Power(xii,2.0) + 
       420.0*Power(rij,3.0)*Power(xii,3.0) + 210.0*Power(rij,4.0)*Power(xii,4.0) + 
       84.0*Power(rij,5.0)*Power(xii,5.0) + 28.0*Power(rij,6.0)*Power(xii,6.0) + 
       8.0*Power(rij,7.0)*Power(xii,7.0) + 2.0*Power(rij,8.0)*Power(xii,8.0))/
     (315.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  return S;
}

double DNuclear_5S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2.0) - (14175.0 + 28350.0*rij*xii + 28350.0*Power(rij,2.0)*Power(xii,2.0) + 
       18900.0*Power(rij,3.0)*Power(xii,3.0) + 9450.0*Power(rij,4.0)*Power(xii,4.0) + 
       3780.0*Power(rij,5.0)*Power(xii,5.0) + 1260.0*Power(rij,6.0)*Power(xii,6.0) + 
       360.0*Power(rij,7.0)*Power(xii,7.0) + 90.0*Power(rij,8.0)*Power(xii,8.0) + 
       20.0*Power(rij,9.0)*Power(xii,9.0) + 4.0*Power(rij,10.0)*Power(xii,10.0))/
     (14175.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  return S;
}

double DNuclear_6S(double rij,double xii)
{
  double S;

  S = 
  Power(rij,-2.0) - (467775.0 + 935550.0*rij*xii + 
       935550.0*Power(rij,2.0)*Power(xii,2.0) + 623700.0*Power(rij,3.0)*Power(xii,3.0) + 
       311850.0*Power(rij,4.0)*Power(xii,4.0) + 124740.0*Power(rij,5.0)*Power(xii,5.0) + 
       41580.0*Power(rij,6.0)*Power(xii,6.0) + 11880.0*Power(rij,7.0)*Power(xii,7.0) + 
       2970.0*Power(rij,8.0)*Power(xii,8.0) + 660.0*Power(rij,9.0)*Power(xii,9.0) + 
       132.0*Power(rij,10.0)*Power(xii,10.0) + 24.0*Power(rij,11.0)*Power(xii,11.0) + 
       4.0*Power(rij,12.0)*Power(xii,12.0))/(467775.*Power(E,2.0*rij*xii)*Power(rij,2.0))
    ;
  return S;
}

t_slater_SS_func (*Slater_SS[SLATER_MAX][SLATER_MAX]) = {
  {  Slater_1S_1S,  Slater_2S_1S,  Slater_3S_1S,  Slater_4S_1S,  Slater_5S_1S,  Slater_6S_1S},
  {  Slater_1S_2S,  Slater_2S_2S,  Slater_3S_2S,  Slater_4S_2S,  Slater_5S_2S,  Slater_6S_2S},
  {  Slater_1S_3S,  Slater_2S_3S,  Slater_3S_3S,  Slater_4S_3S,  Slater_5S_3S,  Slater_6S_3S},
  {  Slater_1S_4S,  Slater_2S_4S,  Slater_3S_4S,  Slater_4S_4S,  Slater_5S_4S,  Slater_6S_4S},
  {  Slater_1S_5S,  Slater_2S_5S,  Slater_3S_5S,  Slater_4S_5S,  Slater_5S_5S,  Slater_6S_5S},
  {  Slater_1S_6S,  Slater_2S_6S,  Slater_3S_6S,  Slater_4S_6S,  Slater_5S_6S,  Slater_6S_6S}
};

t_slater_SS_func (*DSlater_SS[SLATER_MAX][SLATER_MAX]) = {
  {  DSlater_1S_1S,  DSlater_2S_1S,  DSlater_3S_1S,  DSlater_4S_1S,  DSlater_5S_1S,  DSlater_6S_1S},
  {  DSlater_1S_2S,  DSlater_2S_2S,  DSlater_3S_2S,  DSlater_4S_2S,  DSlater_5S_2S,  DSlater_6S_2S},
  {  DSlater_1S_3S,  DSlater_2S_3S,  DSlater_3S_3S,  DSlater_4S_3S,  DSlater_5S_3S,  DSlater_6S_3S},
  {  DSlater_1S_4S,  DSlater_2S_4S,  DSlater_3S_4S,  DSlater_4S_4S,  DSlater_5S_4S,  DSlater_6S_4S},
  {  DSlater_1S_5S,  DSlater_2S_5S,  DSlater_3S_5S,  DSlater_4S_5S,  DSlater_5S_5S,  DSlater_6S_5S},
  {  DSlater_1S_6S,  DSlater_2S_6S,  DSlater_3S_6S,  DSlater_4S_6S,  DSlater_5S_6S,  DSlater_6S_6S}
};

t_slater_NS_func (*Slater_NS[SLATER_MAX]) = {
  Nuclear_1S,  Nuclear_2S,  Nuclear_3S,  Nuclear_4S,  Nuclear_5S,  Nuclear_6S
};

t_slater_NS_func (*DSlater_NS[SLATER_MAX]) = {
  DNuclear_1S,  DNuclear_2S,  DNuclear_3S,  DNuclear_4S,  DNuclear_5S,  DNuclear_6S
};

