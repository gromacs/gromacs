//
//	@(#) delay.s 1.5 11/5/92
//
//
//	GROMACS - Groningen Machine for Chemical Simulation
//	Copyright (c) 1990, 1991, 1992, Groningen University
//
// Delay routine for the SPC-i860 by A. Sijbers, (c) 1991,1992
//
// General description:
//
// The routine delay(ms) delays for the number of milli seconds, specified
// in the argument. This routine is written in assembler to be independent 
// of compiler options and optimalisations. Caching effects are ignored.
//
// C-interface:
//
//   extern void delay(int ms)
//      Delays for ms milli seconds.
//
.globl  _delay
.globl  _delay01
//
//   extern void delay(int ms)
//      delays for ms milli seconds.
//   extern void delay01(int ms01)
//      delays for ms 0.01 milli seconds.
//    in:
//       r16 = ms
//    modifies:
//       r16,r17,r18
//
MS_LOOP		=	6666	// Experimental value for loop count,
				// accuracy: 1 percent
//
_delay:
	adds	1,r0,r18	// decrement	
	subs	0,r16,r0
	bnc	delay_done
	subs	r16,r18,r16
	mov	MS_LOOP,r17
loop:
	subs	0,r17,r0
	bnc	_delay
	subs	r17,r18,r17
	br	loop
	nop
delay_done:
	bri r1
	nop
//
MS01_LOOP	=	66	// Experimental value for loop count,
				// accuracy: 1 percent
//
_delay01:
	adds	1,r0,r18	// decrement	
	subs	0,r16,r0
	bnc	delay01_done
	subs	r16,r18,r16
	mov	MS01_LOOP,r17
loop01:
	subs	0,r17,r0
	bnc	_delay01
	subs	r17,r18,r17
	br	loop01
	nop
delay01_done:
	bri r1
	nop
//
