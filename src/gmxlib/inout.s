//
//	@(#) inout.s 1.5 11/5/92
//
//
//	GROMACS - Groningen Machine for Chemical Simulation
//	Copyright (c) 1990, 1991, 1992, Groningen University
//
// Memory mapped io for the SPC-i860 by A. Sijbers, (c) 1991,1992
//
// General description:
//
// The i860 io buffers (both dual ported rams, idt7130 and hd63310) are 
// organised as bytes located at the lowest byte address in a 64 bit word.
// So to address consecutive bytes of the io buffers, it is necessary to
// increment the i860 byte address by 8. 
//
// This module implements a general method of accessing the io buffers by
// specifying the base address of the buffer and an offset. The base 
// address is the physical address of the buffer, the offset is byte count 
// in the buffers local address space (to address the first byte of an io 
// buffer, specify 0, for the second byte 1 etc.). 
//
// Although it is possible to use any address and buffer size combination
// for the put_io_buf and get_io_buf routines, it is strongly recommended 
// to use only long word aligned addresses and sizes which are multiple 
// of 4 for maximum speed.
//
// C-interface:
//
//   extern void poke_io_buf(int iobase,int offset,int byte);
//      Puts byte at io buffer, pointed to by iobase, location offset.
//   extern int peek_io_buf(int iobase,int offset);
//      Returns the byte at offset from io buffer pointed to by iobase.
//   extern void put_io_buf(int iobase,int offset,void *buf,int bufsize);
//      Puts bufsize bytes of buffer, pointed to by buf, into io buffer, 
//      pointed to by iobase, starting at offset.
//   extern void get_io_buf(int iobase,int offset,void *buf,int bufsize);
//      Puts bufsize bytes from io buffer, pointed to by iobase, location 
//      offset into the buffer pointed to by buf.
//
.globl  _poke_io_buf
.globl  _peek_io_buf
.globl  _put_io_buf
.globl  _get_io_buf
//
//   extern void poke_io_buf(int iobase,int offset,int byte);
//      puts byte at io buffer, pointed to by iobase, location offset.
//    in:
//       r16 = iobase
//       r17 = offset
//       r18 = byte
//    modifies:
//       r31
//
_poke_io_buf:
        shl     3,r17,r31
        addu    r16,r31,r31
        bri     r1
        st.b    r18,0(r31)
//
//   extern int peek_io_buf(int iobase,int offset);
//      returns the byte at offset from io buffer pointed to by iobase.
//    in:
//       r16 = iobase
//       r17 = offset
//    out:
//       r16 = function result
//    modifies:
//       r31
//
_peek_io_buf:
        shl     3,r17,r31
        addu    r16,r31,r31
        ld.b    0(r31),r16
        bri     r1
        and     0xff,r16,r16
//
//   extern void put_io_buf(int iobase,int offset,void *buf,int bufsize);
//      puts bufsize bytes of buffer, pointed to by buf, into io buffer, 
//      pointed to by iobase, starting at offset.
//
//    in:
//       r16 = iobase
//       r17 = offset
//       r18 = bufptr
//       r19 = bufsize
//    modifies:
//       r18,r19,r20,r21,r22,r31
//
_put_io_buf:
	subs	0,r19,r0		// bufsize > 0 ?
	bnc	put_done		// jump if bufsize <= 0
        shl     3,r17,r31
        addu    r16,r31,r31             // calculated base address in r31
	and	3,r18,r20
	bc	put_prologue_done	// don't do prologue if aligned
put_prologue:				// avoid this: very expensive
	ld.b	0(r18),r21
	st.b	r21,0(r31)
	adds	-1,r19,r19
	bte	r0,r19,put_done		// jump if all bytes done
        addu	1,r18,r18
        addu	8,r31,r31
	and	3,r18,r0
	bnc	put_prologue		// until adress long word aligned
put_prologue_done:
	and	3,r19,r22		// save epilogue part in r22
        shr     2,r19,r19               // in steps of 4 bytes (1 long word)
	bte	r0,r19,put_epilogue	// jump if all bytes done
        adds    -1,r19,r19              // loop count
        adds    -1,r0,r20               // loop increment
        bla     r20,r19,put_loop
        nop
put_loop:
        ld.l    0(r18),r21		// fetch long word from memory
        addu    4,r18,r18		// point to next long word
        st.b    r21,0(r31)		// 1st byte to io buffer
        shr     8,r21,r21		// 2nd byte into right position
        st.b    r21,8(r31)		// 2nd byte to io buffer
        shr     8,r21,r21		// 3rd byte into right position 
        st.b    r21,16(r31)             // 3rd byte to io buffer          
        shr     8,r21,r21		// 4th byte into right position   
        st.b    r21,24(r31)             // 4th byte to io buffer          
        bla     r20,r19,put_loop	// decrement loop count
        addu    32,r31,r31		// next io buffer address
put_epilogue:
	bte	r0,r22,put_done		// jump if no epilogue part
put_epilogue1:				// avoid this: very expensive
	ld.b	0(r18),r21
	st.b	r21,0(r31)
	addu	1,r18,r18
        addu	8,r31,r31
	adds	-1,r22,r22
	btne	r0,r22,put_epilogue1	// until epilogue done (r22=0)
put_done:
        bri     r1
        nop
//
//   extern void get_io_buf(int iobase,int offset,void *buf,int bufsize);
//      puts bufsize bytes from io buffer, pointed to by iobase, location 
//      offset into the buffer pointed to by buf.
//    in:
//       r16 = iobase
//       r17 = offset
//       r18 = bufptr
//       r19 = bufsize
//    modifies:
//       r16,r17,r18,r19,r20,r21,r22,r23,r24,r31
//
_get_io_buf:
	subs	0,r19,r0		// bufsize > 0 ?
	bnc	get_done		// jump if bufsize <= 0
        shl     3,r17,r31
        addu    r16,r31,r31             // calculated base address in r31
	and	3,r18,r17
	bc	get_prologue_done	// don't do prologue if aligned
get_prologue:				// avoid this: very expensive
	ld.b	0(r31),r21
	st.b	r21,0(r18)
	adds	-1,r19,r19
	bte	r0,r19,get_done		// jump if all bytes done
        addu	1,r18,r18
        addu	8,r31,r31
	and	3,r18,r0
	bnc	get_prologue		// until adress long word aligned
get_prologue_done:
	and	3,r19,r24		// save epilogue part in r24
        shr     2,r19,r19               // in steps of 4 bytes (1 long word)
	bte	r0,r19,get_epilogue	// jump if all bytes done

	adds    -8,r31,r31              // start at 8(-8+offset) = 0 (offset)
        addu    8,r0,r23                // increment for PFLD load instr.
                                        // now ready for pflds
        pfld.l  r23(r31)++,f0		// Initialise pipe (pipe 1)
        pfld.l  r23(r31)++,f0		// pipe 2
        adds    -1,r19,r19		// loop count
        pfld.l  r23(r31)++,f0		// pipe 3
        adds    -1,r0,r17		// loop increment
        pfld.l  r23(r31)++,f16          // save first result in f16
        shr     8,r0,r0                 // set SC field in psr for SHRD
        pfld.l  r23(r31)++,f17          // next result in f17
        fxfr    f16,r21                 // start in pipeline with defined r21
        bla     r17,r19,get_loop
        nop
get_loop:                               // result of load in r20
        pfld.l  r23(r31)++,f16
        fxfr    f17,r22
        shrd    r21,r0,r20
        pfld.l  r23(r31)++,f17
        fxfr    f16,r21
        shrd    r22,r20,r20
        pfld.l  r23(r31)++,f16
        fxfr    f17,r22
        shrd    r21,r20,r20
        pfld.l  r23(r31)++,f17          // the last are dummy reads at the end
        fxfr    f16,r21
        shrd    r22,r20,r20             // r20 needs to settle
        addu    4,r18,r18
        bla     r17,r19,get_loop
        st.l    r20,-4(r18)
	adds	-32,r31,r31		// correct for last dummy reads
get_epilogue:
	bte	r0,r24,get_done		// jump if no epilogue part
get_epilogue1:
	ld.b	0(r31),r20
	st.b	r20,0(r18)
	addu	1,r18,r18
        addu	8,r31,r31
	adds	-1,r24,r24
	btne	r0,r24,get_epilogue1	// until epilogue done (r24=0)
get_done:
        bri     r1
        nop
//
