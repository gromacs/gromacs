/*
 *       $Id$
 *
 *       This source code is part of
 *
 *        G   R   O   M   A   C   S
 *
 * GROningen MAchine for Chemical Simulations
 *
 *            VERSION 2.0
 * 
 * Copyright (c) 1991-1997
 * BIOSON Research Institute, Dept. of Biophysical Chemistry
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
 * GROup of MAchos and Cynical Suckers
 */
static char *SRCID_ambaled_c = "$Id$";

#include <stdio.h>
#include "led.h"
#include "delay.h"
#include "ambadma.h"

typedef unsigned char	byte;

static byte ledport = INIT_LED;

void put_port (int value)
{
  int i;
  byte l ;

  l = 0x01;
  ledport = value;
  for (i=0; i<4; i++)
    { 
      if ((ledport & l) != 0) ambaLed (i, ON);
      else ambaLed(i, OFF);
      l = l << 1;
    }
}


int get_port ()
{
  return ledport;
}

void put_leds (int value)
{
  put_port ( value & MASK_LEDS);
}

int get_leds()
{
  return (get_port()&MASK_LEDS);
}

void put_led(int nr,int value)
{
  if (value!=0) 
    put_leds(ledport|led_mask(nr)); 
  else 
    put_leds(ledport&~led_mask(nr));
}

int get_led(int nr)
{
  return (ledport&led_mask(nr));
}

void set_led(int nr)
{
  put_led(nr,1);
}

void clr_led(int nr)
{
  put_led(nr,0);
}

void inv_led(int nr)
{
  put_led(nr,get_led(nr)==0);
}

void flash_led(int nr,int count)
{
  for (;count!=0 ;count--) { inv_led(nr); delay(200); }
}

void flash_leds(int mask,int count)
{
  int i;
  for (;count!=0 ;count--) 
    { 
      for (i=0; i<=MAXLED; i++) if ((mask>>i)&1) inv_led(i); 
      delay(200); 
    }
}

