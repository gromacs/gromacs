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
 * Great Red Owns Many ACres of Sand 
 */

#ifndef _led_h
#define _led_h

static char *SRCID_led_h = "$Id$";

#ifdef HAVE_IDENT
#ident	"@(#) led.h 1.8 11/23/92"
#endif /* HAVE_IDENT */

/*
 * The ledport routines are used to control the leds of the SPC-860 board
 *
 * The ledport routines must be initialised by putting a value in the
 * ledport. This can be done by calling put_port. The argument of put_port
 * will be the initial value. Take care that bit 7 has the right value
 * because it spefies whether the boot rom is enabled. In most cases
 * the value of INIT_LED will do: boot rom disabled and all leds dark.
 *
 * Layout of the ledport:
 *
 * +----+----+----+----+----+----+----+----+
 * | D7 | D6 | D5 | D4 | D3 | D2 | D1 | D0 |
 * +----+----+----+----+----+----+----+----+
 *   |    |    |    |    |    |    |    |
 *   |    |    |    |    |    |    |    +-- Led 0 (panel name D2)
 *   |    |    |    |    |    |    +------- Led 1 (panel name D3)
 *   |    |    |    |    |    +------------ Led 2 (panel name D4)
 *   |    |    |    |    +----------------- Led 3 (panel name D5)
 *   |    |    |    +---------------------- Not used
 *   |    |    +--------------------------- Not used
 *   |    +-------------------------------- Not used
 *   +------------------------------------- Boot rom: 0=enabled, 1=disabled
 *
 * All bits: 0=enabled, 1=disabled (inverted!), see remark below.
 *
 * The panel name refers to the name at the leds name on the I860 board.
 * Contrary to the enable/disable state of the output latch, the led
 * routines return and accept values according to the state of the
 * leds where a '1' means enabled and a '0' means disabled. This is a more
 * convenient and logical approach. So get_leds and put_leds invert the
 * actual hardware state in the more convenient state and v.v..
 *
 */

#define MAXLED		3	/* leds 0..MAXLED are available          */
#define INIT_LED	0x00	/* boot rom disabled and all leds dark   */
#define MASK_LEDS	(~(-1<<((MAXLED)+1)))	/* masks out led bits    */
#define led_mask(nr)	((1<<nr)&(MASK_LEDS))	/* generate led bit      */

typedef enum {FORCE_LED,NS_LED,UPDATE_LED,COMM_LED} t_ledid;

extern void put_port(int value);
     /*
      * Lowlevel write to port, converts from software state into hardware
      * state.
      */

extern int get_port(void);
     /*
      * Hardware is write only, so a copy is actually read from. The correct
      * state is only returned when all led addressing is done via this
      * module. Eventually led addressing in the monitor will corrupt this
      * state, but any write to the led port (via this module) will correct
      * this.
      */

extern void put_leds(int value);
     /*
      * Write value into the led port, does not affect the boot rom bit.
      */

extern int get_leds(void);
     /*
      * Returns the status of the leds only.
      */

extern void put_led(int nr,int value);
     /*
      * Puts the led nr (0..MAXLED) in the specified state.
      */

extern int get_led(int nr);
     /*
      * Returns the the state of led nr (0..MAXLED).
      */

extern void set_led(int nr);
     /*
      * Puts the led nr (0..MAXLED) in on ('1') state.
      */

extern void clr_led(int nr);
     /*
      * Puts the led nr (0..MAXLED) in off ('0') state.
      */

extern void inv_led(int nr);
     /*
      * Inverts the state of led nr (0..MAXLED).
      */

extern void flash_led(int nr,int count);
     /*
      * Flashes led nr (0..MAXLED) for count times, count<0 defines virtually
      * forever. This function can be used for example to signal an error
      * state .
      */

extern void flash_leds(int mask,int count);
     /*
      * Flashes the leds specified by mask (every bit set flashes the
      * according led) for count times, count<0 defines virtually forever.
      * This function can be used for example to signal an error state.
      */

#endif	/* _led_h */
