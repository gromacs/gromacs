/*
 *       @(#) copyrgt.c 1.12 9/30/97
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0b
 * 
 * Copyright (c) 1990-1997,
 * BIOSON Research Institute, Dept. of Biophysical Chemistry,
 * University of Groningen, The Netherlands
 *
 * Please refer to:
 * GROMACS: A message-passing parallel molecular dynamics implementation
 * H.J.C. Berendsen, D. van der Spoel and R. van Drunen
 * Comp. Phys. Comm. 91, 43-56 (1995)
 *
 * Also check out our WWW page:
 * http://rugmd0.chem.rug.nl/~gmx
 * or e-mail to:
 * gromacs@chem.rug.nl
 *
 * And Hey:
 * GRoups of Organic Molecules in ACtion for Science
 */

#include "main.h"
#include "fatal.h"
#include "led.h"
#include "typedefs.h"
#include "ns.h"

/* dummylib, to be used on single processor systems, other than i860 */

/* Led routines */
void delay(int ms)					{}
void put_port(int value)				{}
int  get_port(void)					{ return 0; }
void put_leds(int value)				{}
int  get_leds(void)					{ return 0; }
void put_led(int nr,int value)				{}
int  get_led(int nr)					{ return 0; }
void set_led(int nr)					{}
void clr_led(int nr)					{}
void inv_led(int nr)					{}
void flash_led(int nr,int count)			{}
void flash_leds(int mask,int count)			{}

