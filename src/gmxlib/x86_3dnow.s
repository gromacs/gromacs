 ;#
;#  $Id$
;#
;#                This source code is part of
;#
;#                 G   R   O   M   A   C   S
;#
;#          GROningen MAchine for Chemical Simulations
;#
;#                        VERSION 3.1
;# Copyright (c) 1991-2001, University of Groningen, The Netherlands
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;# 
;# If you want to redistribute modifications, please consider that
;# scientific software is very special. Version control is crucial -
;# bugs must be traceable. We will be happy to consider code for
;# inclusion in the official distribution, but derived work must not
;# be called official GROMACS. Details are found in the README & COPYING
;# files - if they are missing, get the official version at www.gromacs.org.
;# 
;# To help us fund GROMACS development, we humbly ask that you cite
;# the papers on the package - you can find them in the top README file.
;# 
;# For more info, check our website at http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
;#

;# This file contains a subset of the gromacs innerloops
;# manually written in assembly to optimize performance
;# on AMD extended 3DNow-enabled processors like Athlon 
;# and later generations. 
;# Erik Lindahl, 2000-2001, erik@theophys.kth.se
;#
;#
;# These files require GNU binutils 2.10 or later, since we
;# use intel syntax for portability, or a recent version 
;# of NASM that understands Extended 3DNow and SSE2 instructions.
;# (NASM is normally only used with MS Visual C++).

;# Since NASM and gnu as disagree on some definitions and use 
;# completely different preprocessing options I have to introduce a
;# trick: NASM uses ';' for comments, while gnu as uses '#' on x86.
;# Gnu as treats ';' as a line break, i.e. ignores it. This is the
;# reason why all comments need both symbols...
;# The source is written for GNU as, with intel syntax. When you use
;# NASM we redefine a couple of things. The false if-statement around 
;# the following code is seen by GNU as (NASM doesn't understant this
;# if syntax), but NASM doesn't see it, so the code inside is only 
;# read by NASM (NASM doesn't understand .if):

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM wants 'dword' only, not 'dword ptr'. 
%define ptr
%macro .equiv 2
   %1 equ %2
%endmacro
; .endif  # End of NASM-specific block

; .intel_syntax noprefix   # Line only read by gnu as


.section .text
	
mm_two:	
	.long 0x40000000
	.long 0x40000000
mm_six:	
	.long 0x40c00000
	.long 0x40c00000
mm_twelve:	
	.long 0x41400000
	.long 0x41400000

	.align 4

.globl check3dnow  ;# try to issue an Extended 3DNow instruction 
.globl _check3dnow
check3dnow:	
_check3dnow:	
	femms
	pswapd mm0,mm0
	femms
	ret

			
.globl vecrecip_3dnow
.globl _vecrecip_3dnow
vecrecip_3dnow:	
_vecrecip_3dnow:	
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx

	mov eax, [ebp + 8]
	mov ebx, [ebp + 12]	
	mov ecx, [ebp + 16]
    mov edx, ecx
    shr ecx, 2 
    jecxz .vecrecip_tail
    emms	
.vecrecip_mainloop:	
    movq mm0,[eax]
	add eax,  8
    pfrcp mm1,mm0
	movq mm4,[eax]
	pswapd mm0,mm0
	add eax,  8 
    pfrcp mm2,mm0
	pswapd mm0,mm0
    pfrcp mm5,mm4
	pswapd mm4,mm4	
	punpckldq mm1,mm2
	pfrcp mm6,mm4
	pswapd mm4,mm4
	pfrcpit1 mm0,mm1
	punpckldq mm5,mm6	
	pfrcpit2 mm0,mm1
    movq [ebx],mm0
	pfrcpit1 mm4,mm5
	add ebx,  8
	pfrcpit2 mm4,mm5	
    movq [ebx],mm4
	add ebx,  8	
    dec ecx
    jecxz .vecrecip_tail
    jmp short .vecrecip_mainloop
.vecrecip_tail:
    mov ecx,edx
    and ecx,3
    jecxz .vecrecip_end
.vecrecip_tailloop:	
    movd mm0,[eax]
	add eax,  4
    pfrcp mm1,mm0
    pfrcpit1 mm0,mm1
    pfrcpit2 mm0,mm1
    movd [ebx],mm0	
	add ebx,  4
	dec ecx
	jecxz .vecrecip_end
	jmp short .vecrecip_tailloop
.vecrecip_end:	
	emms
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	

.globl vecinvsqrt_3dnow
.globl _vecinvsqrt_3dnow
vecinvsqrt_3dnow:	
_vecinvsqrt_3dnow:	
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx

	mov eax, [ebp + 8]
	mov ebx, [ebp + 12]	
	mov ecx, [ebp + 16]
    mov edx, ecx
    shr ecx, 2 
    jecxz .vecinvsqrt_tail
    emms	
.vecinvsqrt_mainloop:	
    movq mm0,[eax]
	add eax,  8
    pfrsqrt mm1,mm0
	movq mm4,[eax]
	pswapd mm0,mm0
	add eax,  8
    pfrsqrt mm2,mm0
	pswapd mm0,mm0
    pfrsqrt mm5,mm4
	pswapd mm4,mm4	
	punpckldq mm1,mm2
	pfrsqrt mm6,mm4
	movq mm3,mm1
	pswapd mm4,mm4
	pfmul mm1,mm1
	punpckldq mm5,mm6	
	pfrsqit1 mm1,mm0
	movq mm7,mm5	
	pfrcpit2 mm1,mm3
	pfmul mm5,mm5
    movq [ebx],mm1
	pfrsqit1 mm5,mm4
	add ebx,  8
	pfrcpit2 mm5,mm7	
    movq [ebx],mm5
	add ebx,  8	
    dec ecx
    jecxz .vecinvsqrt_tail
    jmp short .vecinvsqrt_mainloop
.vecinvsqrt_tail:
    mov ecx,edx
    and ecx,3
    jecxz .vecinvsqrt_end
.vecinvsqrt_tailloop:	
    movd mm0,[eax]
	add eax,  4
    pfrsqrt mm1,mm0
    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0
    pfrcpit2 mm1,mm2
    movd [ebx],mm1		
	add ebx,  4
	dec ecx
	jecxz .vecinvsqrt_end
	jmp short .vecinvsqrt_tailloop
.vecinvsqrt_end:	
	emms
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
	

.globl inl0100_3dnow
.globl _inl0100_3dnow
inl0100_3dnow:	
_inl0100_3dnow:	
.equiv		i0100_nri,		8
.equiv		i0100_iinr,		12
.equiv		i0100_jindex,		16
.equiv		i0100_jjnr,		20
.equiv		i0100_shift,		24
.equiv		i0100_shiftvec,		28
.equiv		i0100_fshift,		32
.equiv		i0100_gid,		36
.equiv		i0100_pos,		40
.equiv		i0100_faction,		44
.equiv		i0100_type,		48
.equiv		i0100_ntype,		52
.equiv		i0100_nbfp,		56
.equiv		i0100_Vnb,		60
	;# stack offsets for local variables 
.equiv		i0100_is3,		0
.equiv		i0100_ii3,		4
.equiv		i0100_ix,		8
.equiv		i0100_iy,		12
.equiv		i0100_iz,		16
.equiv		i0100_vnbtot,		20  
.equiv		i0100_c6,		28  
.equiv		i0100_c12,		36  
.equiv		i0100_six,		44  
.equiv		i0100_twelve,		52  
.equiv		i0100_ntia,		60
.equiv		i0100_innerjjnr,	64
.equiv		i0100_innerk,		68		
.equiv		i0100_fix,		72
.equiv		i0100_fiy,		76
.equiv		i0100_fiz,		80
.equiv		i0100_dx1,		84
.equiv		i0100_dy1,		88
.equiv		i0100_dz1,		92
.equiv		i0100_dx2,		96
.equiv		i0100_dy2,		100
.equiv		i0100_dz2,		104						
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 108						;# local stack space 
	femms
	;# move data to local stack  
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + i0100_six ], mm0
	movq  [esp + i0100_twelve ], mm1
	;# assume we have at least one i particle - start directly 	
.i0100_outer:
	mov   eax, [ebp + i0100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]					;# ebx=shift[n] 
	add   dword ptr [ebp + i0100_shift], 4		;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]			;# ebx=3*is 
	mov   [esp + i0100_is3],ebx    		;# store is3 

	mov   eax, [ebp + i0100_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]			;# move shX/shY to mm0 and shZ to mm1. 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + i0100_iinr]       ;# ecx = pointer into iinr[] 	
	add   dword ptr [ebp + i0100_iinr],  4		;# advance pointer 
	mov   ebx, [ecx]					;# ebx =ii 

	mov   edx, [ebp + i0100_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i0100_ntype]
	shl   edx, 1
	mov   [esp + i0100_ntia], edx

	lea   ebx, [ebx + ebx*2]			;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i0100_pos]		;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]			;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]		;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i0100_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + i0100_ix], mm0	
	movd  [esp + i0100_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + i0100_vnbtot], mm7
	movq  [esp + i0100_fix],    mm7
	movd  [esp + i0100_fiz],    mm7

	mov   eax, [ebp + i0100_jindex]
	mov   ecx, [eax]					;# jindex[n] 
	mov   edx, [eax + 4]				;# jindex[n+1] 
	add   dword ptr [ebp + i0100_jindex],  4
	sub   edx, ecx						;# number of innerloop atoms 

	mov   esi, [ebp + i0100_pos]
	mov   edi, [ebp + i0100_faction]	
	mov   eax, [ebp + i0100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i0100_innerjjnr], eax	;#  pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + i0100_innerk], edx		;# number of innerloop atoms 
	jge   .i0100_unroll_loop
	jmp   .i0100_finish_inner
.i0100_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + i0100_innerjjnr]	;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]				;# eax/ebx=jnr 
	add   dword ptr [esp + i0100_innerjjnr],  8	;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]					;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + i0100_type]
	mov edx, [ecx + eax*4]        		;# type [jnr1] 
	mov ecx, [ecx + ebx*4]				;# type [jnr2] 

	mov esi, [ebp + i0100_nbfp]			;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i0100_ntia]			;# tja = ntia + 2*type 
	add ecx, [esp + i0100_ntia]

	movq mm5, [esi + edx*4]				;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]				;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7					;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7					;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i0100_c6], mm5
	movq [esp + i0100_c12], mm6

	lea   eax, [eax + eax*2]			;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i0100_pos]

	movq  mm0, [esp + i0100_ix]
	movd  mm1, [esp + i0100_iz]	 	
	movq  mm4, [esi + eax*4]			;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0						;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i0100_dx1], mm4	    ;# store dr 
	movd  [esp + i0100_dz1], mm5
	pfmul mm4,mm4						;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5						;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5						;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]			;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0						;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i0100_dx2], mm6	    ;# store dr 
	movd  [esp + i0100_dz2], mm7
	pfmul mm6,mm6						;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7						;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7						;# second rsq in lower mm6 

    pfrcp mm0, mm4						;# lookup reciprocal seed  
    pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        			;# now 4 has rsq and 0 the seed for both pairs. 
                  	        			;# amd 3dnow N-R iteration to get full precision. 
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	;# do potential and fscal
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             		;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5				;# mm5=rinvtwelve 

	pfmul mm5, [esp + i0100_c12]
	pfmul mm4, [esp + i0100_c6]	
	movq mm6, mm5				;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i0100_six]

	pfmul mm5, [esp + i0100_twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5				;# mm0 is total fscal now 	

	prefetchw [esp + i0100_dx1]		;# prefetch i forces to cache 

						;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]		;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i0100_dx1]	;# fetch dr 
	movd mm3,  [esp + i0100_dz1]

	;# update vnbtot  
	pfadd mm6, [esp + i0100_vnbtot]     ;# add the earlier value 
	movq [esp + i0100_vnbtot], mm6      ;# store the sum 

	prefetchw [edi + ebx*4]		;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0				;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i0100_dx2] 		;# fetch dr 
	movd mm5,  [esp + i0100_dz2]
	pfmul mm4, mm1   			;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i0100_fix]
	movd mm1,  [esp + i0100_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i0100_fix], mm0
	movd [esp + i0100_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub  dword ptr [esp + i0100_innerk],  2
	jl    .i0100_finish_inner
	jmp   .i0100_unroll_loop
.i0100_finish_inner:	
	and dword ptr [esp + i0100_innerk],  1
	jnz  .i0100_single_inner
	jmp  .i0100_updateouterdata		
.i0100_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i0100_innerjjnr]
	mov   eax, [eax]					;# eax=jnr offset 

	mov esi, [ebp + i0100_nbfp]
	mov ecx, [ebp + i0100_type]
	mov edx, [ecx + eax*4]        		;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i0100_ntia]			;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]				;# mm5 = 1st c6 		
	movq [esp + i0100_c6], mm5
	movd mm5, [esi + edx*4 + 4]			;# mm5 = 1st c12 		
	movq [esp + i0100_c12], mm5

	mov   esi, [ebp + i0100_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i0100_ix]
	movd  mm1, [esp + i0100_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i0100_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i0100_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5						;# mm4=rsq 
	
    pfrcp mm0,mm4
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0					;# mm4=invsq 
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i0100_c12]
	pfmul mm4, [esp + i0100_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6 
	pfsub mm6, mm4

	pfmul mm4, [esp + i0100_six]

	pfmul mm5, [esp + i0100_twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5    ;# mm0 is total fscal now 

	;# update vnbtot 
	pfadd mm6, [esp + i0100_vnbtot]      ;# add the earlier value 
	movq [esp + i0100_vnbtot], mm6       ;# store the sum   

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache 
	movq mm2,  [esp + i0100_dx1]
	movd mm3,  [esp + i0100_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i0100_fix]
	movd mm1,  [esp + i0100_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i0100_fix], mm0
	movd [esp + i0100_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i0100_updateouterdata:	
	mov   ecx, [esp + i0100_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i0100_fix]
	pfadd mm7, [esp + i0100_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i0100_fshift]    ;# increment fshift force 
	mov   edx, [esp + i0100_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i0100_fix]
	pfadd mm7, [esp + i0100_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + i0100_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i0100_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i0100_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i0100_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + i0100_nri]
	dec ecx
	jecxz .i0100_end
	;# not last, iterate once more! 
	mov [ebp + i0100_nri], ecx
	jmp .i0100_outer
.i0100_end:
	femms
	add esp, 108
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret
	



		
		
.globl inl0110_3dnow
.globl _inl0110_3dnow
inl0110_3dnow:	 
_inl0110_3dnow:	 
.equiv		i0110_nri,		8
.equiv		i0110_iinr,		12
.equiv		i0110_jindex,		16
.equiv		i0110_jjnr,		20
.equiv		i0110_shift,		24
.equiv		i0110_shiftvec,		28
.equiv		i0110_fshift,		32
.equiv		i0110_gid,		36
.equiv		i0110_pos,		40		
.equiv		i0110_faction,		44
.equiv		i0110_type,		48
.equiv		i0110_ntype,		52
.equiv		i0110_nbfp,		56	
.equiv		i0110_Vnb,		60				
.equiv		i0110_nsatoms,		64		
	;# stack offsets for local variables 
.equiv		i0110_is3,		0
.equiv		i0110_ii3,		4
.equiv		i0110_shX,		8
.equiv		i0110_shY,		12 
.equiv		i0110_shZ,		16	
.equiv		i0110_ix,		20
.equiv		i0110_iy,		24
.equiv		i0110_iz,		28	
.equiv		i0110_vnbtot,		32 
.equiv		i0110_c6,		40 
.equiv		i0110_c12,		48 
.equiv		i0110_six,		56 
.equiv		i0110_twelve,		64 
.equiv		i0110_ntia,		72	
.equiv		i0110_innerjjnr0,	76
.equiv		i0110_innerk0,		80		
.equiv		i0110_innerjjnr,	84
.equiv		i0110_innerk,		88	
.equiv		i0110_fix,		92
.equiv		i0110_fiy,		96
.equiv		i0110_fiz,		100
.equiv		i0110_dx1,		104
.equiv		i0110_dy1,		108
.equiv		i0110_dz1,		112
.equiv		i0110_dx2,		116
.equiv		i0110_dy2,		120
.equiv		i0110_dz2,		124					
.equiv		i0110_nsvdwc,		128
.equiv		i0110_nscoul,		132
.equiv		i0110_nsvdw,		136
.equiv		i0110_solnr,		140		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 144		;# local stack space 
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + i0110_six],    mm0
	movq  [esp + i0110_twelve], mm1
	;# assume we have at least one i particle - start directly 		
.i0110_outer:
	mov   eax, [ebp + i0110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i0110_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i0110_is3],ebx    	;# store is3 

	mov   eax, [ebp + i0110_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + i0110_shX], mm0
	movd  [esp + i0110_shZ], mm1

	mov   ecx, [ebp + i0110_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i0110_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + i0110_nsatoms]
	add dword ptr [ebp + i0110_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + i0110_nsvdwc], edx
	mov   [esp + i0110_nscoul], eax
	mov   [esp + i0110_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + i0110_vnbtot], mm7
	mov   [esp + i0110_solnr],  ebx
	
	mov   eax, [ebp + i0110_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i0110_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + i0110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i0110_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + i0110_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + i0110_pos]
	mov   edi, [ebp + i0110_faction]
	
	mov   ecx, [esp + i0110_nsvdwc]
	cmp   ecx,  0
	jnz   .i0110_mno_vdwc
	jmp   .i0110_testvdw
.i0110_mno_vdwc:
	mov   ebx, [esp + i0110_solnr]
	inc   dword ptr [esp + i0110_solnr]

	mov   edx, [ebp + i0110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i0110_ntype]
	shl   edx, 1
	mov   [esp + i0110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i0110_pos]    ;# eax = base of pos[] 
	mov   [esp + i0110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i0110_shX]
	pfadd mm1, [esp + i0110_shZ]
	movq  [esp + i0110_ix], mm0	
	movd  [esp + i0110_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i0110_fix],   mm7
	movd  [esp + i0110_fiz],   mm7

	mov   ecx, [esp + i0110_innerjjnr0]
	mov   [esp + i0110_innerjjnr], ecx
	mov   edx, [esp + i0110_innerk0]
    sub   edx,  2
    mov   [esp + i0110_innerk], edx    ;# number of innerloop atoms 
	jge   .i0110_unroll_vdwc_loop
	jmp   .i0110_finish_vdwc_inner
.i0110_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i0110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i0110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + i0110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i0110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i0110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i0110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 	
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i0110_c6], mm5
	movq [esp + i0110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i0110_pos]

	movq  mm0, [esp + i0110_ix]
	movd  mm1, [esp + i0110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i0110_dx1], mm4	     ;# store dr 
	movd  [esp + i0110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 	         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i0110_dx2], mm6	     ;# store dr 
	movd  [esp + i0110_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrcp mm0, mm4	             ;# lookup reciprocal seed  
    pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
                  	        	;# amd 3dnow N-R iteration to get full precision 
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	 ;# do potential and fscal
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i0110_c12]
	pfmul mm4, [esp + i0110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i0110_six]

	pfmul mm5, [esp + i0110_twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5    ;# mm0 is total fscal now 	

	prefetchw [esp + i0110_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i0110_dx1]	;# fetch dr 
	movd mm3,  [esp + i0110_dz1]

	;# update vnbtot 
	pfadd mm6, [esp + i0110_vnbtot]      ;# add the earlier value 
	movq [esp + i0110_vnbtot], mm6       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i0110_dx2] 	;# fetch dr 
	movd mm5,  [esp + i0110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i0110_fix]
	movd mm1,  [esp + i0110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i0110_fix], mm0
	movd [esp + i0110_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr  [esp + i0110_innerk],  2
	jl    .i0110_finish_vdwc_inner
	jmp   .i0110_unroll_vdwc_loop
.i0110_finish_vdwc_inner:	
	and dword ptr [esp + i0110_innerk],  1
	jnz  .i0110_single_vdwc_inner
	jmp  .i0110_updateouterdata_vdwc		
.i0110_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i0110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + i0110_nbfp]
	mov ecx, [ebp + i0110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i0110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i0110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i0110_c12], mm5

	mov   esi, [ebp + i0110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i0110_ix]
	movd  mm1, [esp + i0110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i0110_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i0110_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrcp mm0,mm4
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	;# mm4=invsq  
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i0110_c12]
	pfmul mm4, [esp + i0110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i0110_six]

	pfmul mm5, [esp + i0110_twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5    ;# mm0 is total fscal now 

	;# update vnbtot 
	pfadd mm6, [esp + i0110_vnbtot]      ;# add the earlier value 
	movq [esp + i0110_vnbtot], mm6       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i0110_dx1]
	movd mm3,  [esp + i0110_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i0110_fix]
	movd mm1,  [esp + i0110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i0110_fix], mm0
	movd [esp + i0110_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i0110_updateouterdata_vdwc:	
	mov   ecx, [esp + i0110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i0110_fix]
	pfadd mm7, [esp + i0110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i0110_fshift]    ;# increment fshift force 
	mov   edx, [esp + i0110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i0110_fix]
	pfadd mm7, [esp + i0110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i0110_nsvdwc]
	jz  .i0110_testvdw
	jmp .i0110_mno_vdwc
.i0110_testvdw:	
	mov  ebx,  [esp + i0110_nscoul]
	add  [esp + i0110_solnr],  ebx

	mov  ecx, [esp + i0110_nsvdw]
	cmp  ecx,  0
	jnz  .i0110_mno_vdw
	jmp  .i0110_last_mno
.i0110_mno_vdw:
	mov   ebx,  [esp + i0110_solnr]
	inc   dword ptr [esp + i0110_solnr]

	mov   edx, [ebp + i0110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i0110_ntype]
	shl   edx, 1
	mov   [esp + i0110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i0110_pos]    ;# eax = base of pos[] 
	mov   [esp + i0110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i0110_shX]
	pfadd mm1, [esp + i0110_shZ]
	movq  [esp + i0110_ix], mm0	
	movd  [esp + i0110_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i0110_fix],   mm7
	movd  [esp + i0110_fiz],   mm7

	mov   ecx, [esp + i0110_innerjjnr0]
	mov   [esp + i0110_innerjjnr], ecx
	mov   edx, [esp + i0110_innerk0]
    sub   edx,  2
    mov   [esp + i0110_innerk], edx    ;# number of innerloop atoms 
	jge   .i0110_unroll_vdw_loop
	jmp   .i0110_finish_vdw_inner
.i0110_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i0110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i0110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + i0110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i0110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i0110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i0110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i0110_c6], mm5
	movq [esp + i0110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i0110_pos]

	movq  mm0, [esp + i0110_ix]
	movd  mm1, [esp + i0110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i0110_dx1], mm4	     ;# store dr 
	movd  [esp + i0110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i0110_dx2], mm6	     ;# store dr 
	movd  [esp + i0110_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrcp mm0, mm4	             ;# lookup reciprocal seed  
    pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
                  	        	;# amd 3dnow N-R iteration to get full precision 
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	 ;# do potential and fscal
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i0110_c12]
	pfmul mm4, [esp + i0110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i0110_six]

	pfmul mm5, [esp + i0110_twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5    ;# mm0 is total fscal now 	

	prefetchw [esp + i0110_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i0110_dx1]	;# fetch dr 
	movd mm3,  [esp + i0110_dz1]

	;# update vnbtot 
	pfadd mm6, [esp + i0110_vnbtot]      ;# add the earlier value 
	movq [esp + i0110_vnbtot], mm6       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i0110_dx2] 	;# fetch dr 
	movd mm5,  [esp + i0110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i0110_fix]
	movd mm1,  [esp + i0110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i0110_fix], mm0
	movd [esp + i0110_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	;# should we do one more iteration? 
	sub dword ptr  [esp + i0110_innerk],  2
	jl    .i0110_finish_vdw_inner
	jmp   .i0110_unroll_vdw_loop
.i0110_finish_vdw_inner:	
	and dword ptr [esp + i0110_innerk],  1
	jnz  .i0110_single_vdw_inner
	jmp  .i0110_updateouterdata_vdw		
.i0110_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i0110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + i0110_nbfp]
	mov ecx, [ebp + i0110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i0110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i0110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i0110_c12], mm5

	mov   esi, [ebp + i0110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i0110_ix]
	movd  mm1, [esp + i0110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i0110_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i0110_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrcp mm0,mm4
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	;# mm4=invsq  
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i0110_c12]
	pfmul mm4, [esp + i0110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i0110_six]

	pfmul mm5, [esp + i0110_twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5    ;# mm0 is total fscal now 

	;# update vnbtot 
	pfadd mm6, [esp + i0110_vnbtot]      ;# add the earlier value 
	movq [esp + i0110_vnbtot], mm6       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i0110_dx1]
	movd mm3,  [esp + i0110_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i0110_fix]
	movd mm1,  [esp + i0110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i0110_fix], mm0
	movd [esp + i0110_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i0110_updateouterdata_vdw:	
	mov   ecx, [esp + i0110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i0110_fix]
	pfadd mm7, [esp + i0110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i0110_fshift]    ;# increment fshift force 
	mov   edx, [esp + i0110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i0110_fix]
	pfadd mm7, [esp + i0110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i0110_nsvdw]
	jz  .i0110_last_mno
	jmp .i0110_mno_vdw
	
.i0110_last_mno:	
	mov   edx, [ebp + i0110_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i0110_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i0110_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i0110_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + i0110_nri]
	dec ecx
	jecxz .i0110_end
	;# not last, iterate once more! 
	mov [ebp + i0110_nri], ecx
	jmp .i0110_outer
.i0110_end:
	femms
	add esp, 144
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret



.globl inl0300_3dnow
.globl _inl0300_3dnow
inl0300_3dnow:	
_inl0300_3dnow:	
.equiv		i0300_nri,		8
.equiv		i0300_iinr,		12
.equiv		i0300_jindex,		16
.equiv		i0300_jjnr,		20
.equiv		i0300_shift,		24
.equiv		i0300_shiftvec,		28
.equiv		i0300_fshift,		32
.equiv		i0300_gid,		36
.equiv		i0300_pos,		40		
.equiv		i0300_faction,		44
.equiv		i0300_type,		48
.equiv		i0300_ntype,		52
.equiv		i0300_nbfp,		56	
.equiv		i0300_Vnb,		60
.equiv		i0300_tabscale,		64
.equiv		i0300_VFtab,		68
	;# stack offsets for local variables 
.equiv		i0300_is3,		0
.equiv		i0300_ii3,		4
.equiv		i0300_ix,		8
.equiv		i0300_iy,		12
.equiv		i0300_iz,		16
.equiv		i0300_vnbtot,		20 
.equiv		i0300_c6,		28 
.equiv		i0300_c12,		36 
.equiv		i0300_two,		44 
.equiv		i0300_n1,		52 
.equiv		i0300_tsc,		60 
.equiv		i0300_ntia,		68
.equiv		i0300_innerjjnr,	72
.equiv		i0300_innerk,		76		
.equiv		i0300_fix,		80
.equiv		i0300_fiy,		84
.equiv		i0300_fiz,		88
.equiv		i0300_dx1,		92
.equiv		i0300_dy1,		96
.equiv		i0300_dz1,		100
.equiv		i0300_dx2,		104
.equiv		i0300_dy2,		108
.equiv		i0300_dz2,		112
    push ebp
    mov ebp,esp
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 116		;# local stack space 
	femms
	;# move data to local stack  
	movq  mm0, [mm_two]
	movd  mm3, [ebp + i0300_tabscale]
	movq  [esp + i0300_two],    mm0
	punpckldq mm3,mm3
	movq  [esp + i0300_tsc], mm3	
	;# assume we have at least one i particle - start directly 	
.i0300_outer:
	mov   eax, [ebp + i0300_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i0300_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i0300_is3],ebx    	;# store is3 

	mov   eax, [ebp + i0300_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + i0300_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i0300_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i0300_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i0300_ntype]
	shl   edx, 1
	mov   [esp + i0300_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i0300_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i0300_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + i0300_ix], mm0	
	movd  [esp + i0300_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + i0300_vnbtot], mm7
	movq  [esp + i0300_fix],    mm7
	movd  [esp + i0300_fiz],    mm7

	mov   eax, [ebp + i0300_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i0300_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + i0300_pos]
	mov   edi, [ebp + i0300_faction]	
	mov   eax, [ebp + i0300_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i0300_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + i0300_innerk], edx    ;# number of innerloop atoms 
	jge   .i0300_unroll_loop
	jmp   .i0300_finish_inner
.i0300_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + i0300_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i0300_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best  
	
	mov ecx, [ebp + i0300_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i0300_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i0300_ntia]	     ;# tja = ntia + 2*type  
	add ecx, [esp + i0300_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i0300_c6], mm5
	movq [esp + i0300_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i0300_pos]

	movq  mm0, [esp + i0300_ix]
	movd  mm1, [esp + i0300_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i0300_dx1], mm4	     ;# store dr 
	movd  [esp + i0300_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i0300_dx2], mm6	     ;# store dr 
	movd  [esp + i0300_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 
	;# do potential and fscal 
	pfmul mm1, [esp + i0300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i0300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i0300_VFtab]
	;# dispersion table 
	mov ecx, [esp + i0300_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i0300_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i0300_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	movq mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i0300_vnbtot]      ;# add the earlier value 
	movq [esp + i0300_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + i0300_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + i0300_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i0300_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijD+ fijR 

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + i0300_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i0300_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i0300_dx1]	;# fetch dr 
	movd mm3,  [esp + i0300_dz1]

	;# update vnbtot 
	pfadd mm5, [esp + i0300_vnbtot]      ;# add the earlier value 
	movq [esp + i0300_vnbtot], mm5       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i0300_dx2] 	;# fetch dr 
	movd mm5,  [esp + i0300_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i0300_fix]
	movd mm1,  [esp + i0300_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i0300_fix], mm0
	movd [esp + i0300_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i0300_innerk],  2
	jl    .i0300_finish_inner
	jmp   .i0300_unroll_loop
.i0300_finish_inner:	
	and dword ptr [esp + i0300_innerk],  1
	jnz  .i0300_single_inner
	jmp  .i0300_updateouterdata		
.i0300_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i0300_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + i0300_nbfp]
	mov ecx, [ebp + i0300_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i0300_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i0300_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i0300_c12], mm5

	mov   esi, [ebp + i0300_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i0300_ix]
	movd  mm1, [esp + i0300_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i0300_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i0300_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i0300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i0300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i0300_VFtab]
	mov ecx, [esp + i0300_n1]
	shl ecx, 3
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i0300_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	movq mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i0300_vnbtot]      ;# add the earlier value 
	movq [esp + i0300_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i0300_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i0300_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# update vnbtot 
	pfadd mm5, [esp + i0300_vnbtot]      ;# add the earlier value 
	movq [esp + i0300_vnbtot], mm5       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i0300_dx1]
	movd mm3,  [esp + i0300_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i0300_fix]
	movd mm1,  [esp + i0300_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i0300_fix], mm0
	movd [esp + i0300_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i0300_updateouterdata:	
	mov   ecx, [esp + i0300_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i0300_fix]
	pfadd mm7, [esp + i0300_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i0300_fshift]    ;# increment fshift force 
	mov   edx, [esp + i0300_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i0300_fix]
	pfadd mm7, [esp + i0300_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + i0300_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i0300_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i0300_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i0300_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + i0300_nri]
	dec ecx
	jecxz .i0300_end
	;# not last, iterate once more! 
	mov [ebp + i0300_nri], ecx
	jmp .i0300_outer
.i0300_end:
	femms
	add esp, 116
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


			
	
.globl inl0310_3dnow
.globl _inl0310_3dnow
inl0310_3dnow:	
_inl0310_3dnow:	
.equiv		i0310_nri,		8
.equiv		i0310_iinr,		12
.equiv		i0310_jindex,		16
.equiv		i0310_jjnr,		20
.equiv		i0310_shift,		24
.equiv		i0310_shiftvec,		28
.equiv		i0310_fshift,		32
.equiv		i0310_gid,		36
.equiv		i0310_pos,		40		
.equiv		i0310_faction,		44
.equiv		i0310_type,		48
.equiv		i0310_ntype,		52
.equiv		i0310_nbfp,		56	
.equiv		i0310_Vnb,		60
.equiv		i0310_tabscale,		64
.equiv		i0310_VFtab,		68
.equiv		i0310_nsatoms,		72		
	;# stack offsets for local variables 
.equiv		i0310_is3,		0
.equiv		i0310_ii3,		4
.equiv		i0310_shX,		8
.equiv		i0310_shY,		12 
.equiv		i0310_shZ,		16	
.equiv		i0310_ix,		20
.equiv		i0310_iy,		24
.equiv		i0310_iz,		28	
.equiv		i0310_vnbtot,		32 
.equiv		i0310_c6,		40	 
.equiv		i0310_c12,		48 
.equiv		i0310_two,		56 
.equiv		i0310_n1,		64 
.equiv		i0310_tsc,		72 
.equiv		i0310_ntia,		80	
.equiv		i0310_innerjjnr0,	84
.equiv		i0310_innerk0,		88		
.equiv		i0310_innerjjnr,	92
.equiv		i0310_innerk,		96	
.equiv		i0310_fix,		100
.equiv		i0310_fiy,		104
.equiv		i0310_fiz,		108
.equiv		i0310_dx1,		112
.equiv		i0310_dy1,		116
.equiv		i0310_dz1,		120
.equiv		i0310_dx2,		124
.equiv		i0310_dy2,		128
.equiv		i0310_dz2,		132
.equiv		i0310_nsvdwc,		136
.equiv		i0310_nscoul,		140
.equiv		i0310_nsvdw,		144
.equiv		i0310_solnr,		148		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 152		;# local stack space 
	femms
	movq  mm0, [mm_two]
	movd  mm3, [ebp + i0310_tabscale]
	movq  [esp + i0310_two],    mm0
	punpckldq mm3,mm3
	movq  [esp + i0310_tsc], mm3	
	
	;# assume we have at least one i particle - start directly 		
.i0310_outer:
	mov   eax, [ebp + i0310_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i0310_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i0310_is3],ebx    	;# store is3 

	mov   eax, [ebp + i0310_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + i0310_shX], mm0
	movd  [esp + i0310_shZ], mm1

	mov   ecx, [ebp + i0310_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i0310_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + i0310_nsatoms]
	add dword ptr [ebp + i0310_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + i0310_nsvdwc], edx
	mov   [esp + i0310_nscoul], eax
	mov   [esp + i0310_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + i0310_vnbtot], mm7
	mov   [esp + i0310_solnr],  ebx
	
	mov   eax, [ebp + i0310_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i0310_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + i0310_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i0310_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + i0310_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + i0310_pos]
	mov   edi, [ebp + i0310_faction]
	
	mov   ecx, [esp + i0310_nsvdwc]
	cmp   ecx,  0
	jnz   .i0310_mno_vdwc
	jmp   .i0310_testvdw
.i0310_mno_vdwc:
	mov   ebx,  [esp + i0310_solnr]
	inc   dword ptr [esp + i0310_solnr]

	mov   edx, [ebp + i0310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i0310_ntype]
	shl   edx, 1
	mov   [esp + i0310_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i0310_pos]    ;# eax = base of pos[] 
	mov   [esp + i0310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i0310_shX]
	pfadd mm1, [esp + i0310_shZ]
	movq  [esp + i0310_ix], mm0	
	movd  [esp + i0310_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i0310_fix],   mm7
	movd  [esp + i0310_fiz],   mm7

	mov   ecx, [esp + i0310_innerjjnr0]
	mov   [esp + i0310_innerjjnr], ecx
	mov   edx, [esp + i0310_innerk0]
    sub   edx,  2
    mov   [esp + i0310_innerk], edx    ;# number of innerloop atoms 
	jge   .i0310_unroll_vdwc_loop
	jmp   .i0310_finish_vdwc_inner
.i0310_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i0310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i0310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i0310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i0310_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i0310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i0310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i0310_c6], mm5
	movq [esp + i0310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i0310_pos]

	movq  mm0, [esp + i0310_ix]
	movd  mm1, [esp + i0310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i0310_dx1], mm4	     ;# store dr 
	movd  [esp + i0310_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i0310_dx2], mm6	     ;# store dr 
	movd  [esp + i0310_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 
	;# do potential and fscal 
	pfmul mm1, [esp + i0310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i0310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i0310_VFtab]
	;# dispersion table 
	mov ecx, [esp + i0310_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i0310_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i0310_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	movq mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i0310_vnbtot]      ;# add the earlier value 
	movq [esp + i0310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + i0310_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + i0310_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i0310_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijD+ fijR 

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + i0310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i0310_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i0310_dx1]	;# fetch dr 
	movd mm3,  [esp + i0310_dz1]

	;# update vnbtot 
	pfadd mm5, [esp + i0310_vnbtot]      ;# add the earlier value 
	movq [esp + i0310_vnbtot], mm5       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i0310_dx2] 	;# fetch dr 
	movd mm5,  [esp + i0310_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i0310_fix]
	movd mm1,  [esp + i0310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i0310_fix], mm0
	movd [esp + i0310_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
		
	;# should we do one more iteration? 
	sub dword ptr [esp + i0310_innerk],  2
	jl    .i0310_finish_vdwc_inner
	jmp   .i0310_unroll_vdwc_loop
.i0310_finish_vdwc_inner:	
	and dword ptr [esp + i0310_innerk],  1
	jnz  .i0310_single_vdwc_inner
	jmp  .i0310_updateouterdata_vdwc		
.i0310_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i0310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + i0310_nbfp]
	mov ecx, [ebp + i0310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i0310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i0310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i0310_c12], mm5

	mov   esi, [ebp + i0310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i0310_ix]
	movd  mm1, [esp + i0310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i0310_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i0310_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i0310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i0310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i0310_VFtab]
	mov ecx, [esp + i0310_n1]
	shl ecx, 3
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i0310_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	movq mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i0310_vnbtot]      ;# add the earlier value 
	movq [esp + i0310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i0310_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i0310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# update vnbtot 
	pfadd mm5, [esp + i0310_vnbtot]      ;# add the earlier value 
	movq [esp + i0310_vnbtot], mm5       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i0310_dx1]
	movd mm3,  [esp + i0310_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i0310_fix]
	movd mm1,  [esp + i0310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i0310_fix], mm0
	movd [esp + i0310_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i0310_updateouterdata_vdwc:	
	mov   ecx, [esp + i0310_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i0310_fix]
	pfadd mm7, [esp + i0310_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i0310_fshift]    ;# increment fshift force 
	mov   edx, [esp + i0310_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i0310_fix]
	pfadd mm7, [esp + i0310_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i0310_nsvdwc]
	jz  .i0310_testvdw
	jmp .i0310_mno_vdwc
.i0310_testvdw:	
	mov  ebx,  [esp + i0310_nscoul]
	add  [esp + i0310_solnr],  ebx

	mov  ecx, [esp + i0310_nsvdw]
	cmp  ecx,  0
	jnz  .i0310_mno_vdw
	jmp  .i0310_last_mno
.i0310_mno_vdw:
	mov   ebx,  [esp + i0310_solnr]
	inc   dword ptr [esp + i0310_solnr]

	mov   edx, [ebp + i0310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i0310_ntype]
	shl   edx, 1
	mov   [esp + i0310_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i0310_pos]    ;# eax = base of pos[] 
	mov   [esp + i0310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i0310_shX]
	pfadd mm1, [esp + i0310_shZ]
	movq  [esp + i0310_ix], mm0	
	movd  [esp + i0310_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i0310_fix],   mm7
	movd  [esp + i0310_fiz],   mm7

	mov   ecx, [esp + i0310_innerjjnr0]
	mov   [esp + i0310_innerjjnr], ecx
	mov   edx, [esp + i0310_innerk0]
    sub   edx,  2
    mov   [esp + i0310_innerk], edx    ;# number of innerloop atoms 
	jge   .i0310_unroll_vdw_loop
	jmp   .i0310_finish_vdw_inner
.i0310_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i0310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i0310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i0310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i0310_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i0310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i0310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i0310_c6], mm5
	movq [esp + i0310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i0310_pos]

	movq  mm0, [esp + i0310_ix]
	movd  mm1, [esp + i0310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i0310_dx1], mm4	     ;# store dr 
	movd  [esp + i0310_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i0310_dx2], mm6	     ;# store dr 
	movd  [esp + i0310_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 
	;# do potential and fscal 
	pfmul mm1, [esp + i0310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i0310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i0310_VFtab]
	;# dispersion table 
	mov ecx, [esp + i0310_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i0310_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i0310_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	movq mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i0310_vnbtot]      ;# add the earlier value 
	movq [esp + i0310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + i0310_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + i0310_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i0310_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijD+ fijR 

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + i0310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i0310_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i0310_dx1]	;# fetch dr 
	movd mm3,  [esp + i0310_dz1]

	;# update vnbtot 
	pfadd mm5, [esp + i0310_vnbtot]      ;# add the earlier value 
	movq [esp + i0310_vnbtot], mm5       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i0310_dx2] 	;# fetch dr 
	movd mm5,  [esp + i0310_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i0310_fix]
	movd mm1,  [esp + i0310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i0310_fix], mm0
	movd [esp + i0310_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i0310_innerk],  2
	jl    .i0310_finish_vdw_inner
	jmp   .i0310_unroll_vdw_loop
.i0310_finish_vdw_inner:	
	and dword ptr [esp + i0310_innerk],  1
	jnz  .i0310_single_vdw_inner
	jmp  .i0310_updateouterdata_vdw		
.i0310_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i0310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + i0310_nbfp]
	mov ecx, [ebp + i0310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i0310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i0310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i0310_c12], mm5

	mov   esi, [ebp + i0310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i0310_ix]
	movd  mm1, [esp + i0310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i0310_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i0310_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i0310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i0310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i0310_VFtab]
	mov ecx, [esp + i0310_n1]
	shl ecx, 3
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i0310_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	movq mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i0310_vnbtot]      ;# add the earlier value 
	movq [esp + i0310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i0310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i0310_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i0310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# update vnbtot 
	pfadd mm5, [esp + i0310_vnbtot]      ;# add the earlier value 
	movq [esp + i0310_vnbtot], mm5       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i0310_dx1]
	movd mm3,  [esp + i0310_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i0310_fix]
	movd mm1,  [esp + i0310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i0310_fix], mm0
	movd [esp + i0310_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i0310_updateouterdata_vdw:	
	mov   ecx, [esp + i0310_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i0310_fix]
	pfadd mm7, [esp + i0310_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i0310_fshift]    ;# increment fshift force 
	mov   edx, [esp + i0310_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i0310_fix]
	pfadd mm7, [esp + i0310_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i0310_nsvdw]
	jz  .i0310_last_mno
	jmp .i0310_mno_vdw
	
.i0310_last_mno:	
	mov   edx, [ebp + i0310_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i0310_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i0310_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i0310_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + i0310_nri]
	dec ecx
	jecxz .i0310_end
	;# not last, iterate once more! 
	mov [ebp + i0310_nri], ecx
	jmp .i0310_outer
.i0310_end:
	femms
	add esp, 152
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl inl1000_3dnow
.globl _inl1000_3dnow
inl1000_3dnow:	
_inl1000_3dnow:	
.equiv		i1000_nri,		8
.equiv		i1000_iinr,		12
.equiv		i1000_jindex,		16
.equiv		i1000_jjnr,		20
.equiv		i1000_shift,		24
.equiv		i1000_shiftvec,		28
.equiv		i1000_fshift,		32
.equiv		i1000_gid,		36
.equiv		i1000_pos,		40		
.equiv		i1000_faction,		44
.equiv		i1000_charge,		48
.equiv		i1000_facel,		52
.equiv		i1000_Vc,		56			
	;# stack offsets for local variables 
.equiv		i1000_is3,		0
.equiv		i1000_ii3,		4
.equiv		i1000_ix,		8
.equiv		i1000_iy,		12
.equiv		i1000_iz,		16
.equiv		i1000_iq,		20		
.equiv		i1000_vctot,		28 
.equiv		i1000_innerjjnr,	36
.equiv		i1000_innerk,		40		
.equiv		i1000_fix,		44
.equiv		i1000_fiy,		48
.equiv		i1000_fiz,		52
.equiv		i1000_dx1,		56
.equiv		i1000_dy1,		60
.equiv		i1000_dz1,		64
.equiv		i1000_dx2,		68
.equiv		i1000_dy2,		72
.equiv		i1000_dz2,		76			
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 80		;# 80 bytes local stack space 
	femms
	;# assume we have at least one i particle - start directly 	
.i1000_outer:
	mov   eax, [ebp + i1000_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i1000_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i1000_is3],ebx    	;# store is3 

	mov   eax, [ebp + i1000_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + i1000_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i1000_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i1000_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i1000_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i1000_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1000_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i1000_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + i1000_ix], mm0	
	movd  [esp + i1000_iz], mm1	
				
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i1000_vctot], mm7
	movq  [esp + i1000_fix],   mm7
	movd  [esp + i1000_fiz],   mm7

	mov   eax, [ebp + i1000_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i1000_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + i1000_pos]
	mov   edi, [ebp + i1000_faction]	
	mov   eax, [ebp + i1000_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i1000_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + i1000_innerk], edx    ;# number of innerloop atoms 
	jge   .i1000_unroll_loop
	jmp   .i1000_finish_inner
.i1000_unroll_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i1000_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i1000_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i1000_charge]    ;# base of charge[] 
	movq mm5, [esp + i1000_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
	movd mm7, [ecx + ebx*4]  	 ;# charge[jnr2] 
	punpckldq mm3,mm7	     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + i1000_ix]
	movd  mm1, [esp + i1000_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i1000_dx1], mm4	     ;# store dr 
	movd  [esp + i1000_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i1000_dx2], mm6	     ;# store dr 
	movd  [esp + i1000_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt
	 ;# do potential and fscal
	 
	prefetchw [esp + i1000_dx1]	;# prefetch i forces to cache 
	
	pfmul mm3,mm1		;# 3 has both vcoul 
	pfmul mm0,mm3		;# 0 has both fscal 

	;# update vctot 

	pfadd mm3, [esp + i1000_vctot]      ;# add the earlier value  
	movq [esp + i1000_vctot], mm3       ;# store the sum 
	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i1000_dx1]	;# fetch dr 
	movd mm3,  [esp + i1000_dz1]
	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i1000_dx2] 	;# fetch dr 
	movd mm5,  [esp + i1000_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i1000_fix]
	movd mm1,  [esp + i1000_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i1000_fix], mm0
	movd [esp + i1000_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i1000_innerk],  2
	jl    .i1000_finish_inner
	jmp   .i1000_unroll_loop
.i1000_finish_inner:	
	and dword ptr [esp + i1000_innerk],  1
	jnz  .i1000_single_inner
	jmp  .i1000_updateouterdata		
.i1000_single_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i1000_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i1000_charge]
	movd mm6, [esp + i1000_iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;# mm6=qq 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + i1000_ix]
	movd  mm1, [esp + i1000_iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + i1000_dx1], mm0
	pfmul mm0,mm0
	movd  [esp + i1000_dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;# mm0=rsq 
	
	pfrsqrt mm1,mm0
	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  
	;# update vctot 
	pfadd mm6, [esp + i1000_vctot]
	movq [esp + i1000_vctot], mm6
	;# spread fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + i1000_dx1]
	movd mm1,  [esp + i1000_dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	;# update i particle force 
	movq mm2,  [esp + i1000_fix]
	movd mm3,  [esp + i1000_fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1000_fix], mm2
	movd [esp + i1000_fiz], mm3
	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	;# done! 
.i1000_updateouterdata:	
	mov   ecx, [esp + i1000_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1000_fix]
	pfadd mm7, [esp + i1000_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i1000_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1000_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1000_fix]
	pfadd mm7, [esp + i1000_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + i1000_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i1000_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i1000_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i1000_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + i1000_nri]
	dec ecx
	jecxz .i1000_end
	;# not last, iterate once more! 
	mov [ebp + i1000_nri], ecx
	jmp .i1000_outer
.i1000_end:
	femms
	add esp, 80
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl inl1010_3dnow
.globl _inl1010_3dnow
inl1010_3dnow:	
_inl1010_3dnow:	
.equiv		i1010_nri,		8
.equiv		i1010_iinr,		12
.equiv		i1010_jindex,		16
.equiv		i1010_jjnr,		20
.equiv		i1010_shift,		24
.equiv		i1010_shiftvec,		28
.equiv		i1010_fshift,		32
.equiv		i1010_gid,		36
.equiv		i1010_pos,		40		
.equiv		i1010_faction,		44
.equiv		i1010_charge,		48
.equiv		i1010_facel,		52
.equiv		i1010_Vc,		56
.equiv		i1010_nsatoms,		60		
	;# stack offsets for local variables 
.equiv		i1010_is3,		0
.equiv		i1010_ii3,		4
.equiv		i1010_shX,		8
.equiv		i1010_shY,		12 
.equiv		i1010_shZ,		16	
.equiv		i1010_ix,		20
.equiv		i1010_iy,		24
.equiv		i1010_iz,		28	
.equiv		i1010_iq,		32 
.equiv		i1010_vctot,		40 
.equiv		i1010_innerjjnr0,	48
.equiv		i1010_innerk0,		52		
.equiv		i1010_innerjjnr,	56
.equiv		i1010_innerk,		60		
.equiv		i1010_fix,		64
.equiv		i1010_fiy,		68
.equiv		i1010_fiz,		72
.equiv		i1010_dx1,		76
.equiv		i1010_dy1,		80
.equiv		i1010_dz1,		84
.equiv		i1010_dx2,		88
.equiv		i1010_dy2,		92
.equiv		i1010_dz2,		96
.equiv		i1010_nscoul,		100
.equiv		i1010_solnr,		104		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 108		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	
	add dword ptr [ebp + i1010_nsatoms],  8

.i1010_outer:
	mov   eax, [ebp + i1010_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i1010_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i1010_is3],ebx    	;# store is3 

	mov   eax, [ebp + i1010_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + i1010_shX], mm0
	movd  [esp + i1010_shZ], mm1

	mov   ecx, [ebp + i1010_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i1010_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + i1010_nsatoms]
	mov   ecx, [eax]
	add dword ptr [ebp + i1010_nsatoms],  12
	mov   [esp + i1010_nscoul], ecx

	;# clear potential
	pxor  mm7,mm7
	movq  [esp + i1010_vctot], mm7
	mov   [esp + i1010_solnr], ebx
	
	mov   eax, [ebp + i1010_jindex]	;# current pointer to jindex list
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i1010_jindex],  4 ;# advance pointer
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + i1010_jjnr]
	shl   ecx, 2
	add   eax, ecx          ;# pointer to index of the first j atom
	mov   [esp + i1010_innerjjnr0], eax     ;# save pointer to jjnr[nj0] 

	mov   [esp + i1010_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + i1010_pos]
	mov   edi, [ebp + i1010_faction]

	mov   ecx, [esp + i1010_nscoul]
	cmp   ecx,  0
	jnz   .i1010_mno_coul
	jmp   .i1010_last_mno
.i1010_mno_coul:				
	mov   ebx,  [esp + i1010_solnr]
	inc   dword ptr [esp + i1010_solnr]
	mov   edx, [ebp + i1010_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i1010_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i1010_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1010_pos]    ;# eax = base pointer of pos[] 
	mov   [esp + i1010_ii3], ebx    ;# store ii3 
	
	movq  mm0, [eax + ebx*4]        ;# load x and y coords to mm0 
	movd  mm1, [eax + ebx*4 + 8]    ;# load z coord to mm1
	pfadd mm0, [esp + i1010_shX]    ;# add shift vector
	pfadd mm1, [esp + i1010_shZ]
	movq  [esp + i1010_ix], mm0     ;# store shifted coords
	movd  [esp + i1010_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i1010_fix],   mm7
	movd  [esp + i1010_fiz],   mm7

	mov   ecx, [esp + i1010_innerjjnr0]
	mov   [esp + i1010_innerjjnr], ecx
	mov   edx, [esp + i1010_innerk0]
        sub   edx,  2
        mov   [esp + i1010_innerk], edx    ;# number of innerloop atoms 
	jge   .i1010_unroll_coul_loop
	jmp   .i1010_finish_coul_inner
.i1010_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i1010_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i1010_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i1010_charge]    ;# base of charge[] 
	movq mm5, [esp + i1010_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
	movd mm7, [ecx + ebx*4]  	 ;# charge[jnr2] 
	punpckldq mm3,mm7	     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + i1010_ix]
	movd  mm1, [esp + i1010_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i1010_dx1], mm4	     ;# store dr 
	movd  [esp + i1010_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i1010_dx2], mm6	     ;# store dr 
	movd  [esp + i1010_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	prefetchw [esp + i1010_dx1]	;# prefetch i forces to cache 
	
	pfmul mm3,mm1		;# 3 has both vcoul 
	pfmul mm0,mm3		;# 0 has both fscal 

	;# update vctot 

	pfadd mm3, [esp + i1010_vctot]      ;# add the earlier value  
	movq [esp + i1010_vctot], mm3       ;# store the sum 
	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i1010_dx1]	;# fetch dr 
	movd mm3,  [esp + i1010_dz1]
	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i1010_dx2] 	;# fetch dr 
	movd mm5,  [esp + i1010_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i1010_fix]
	movd mm1,  [esp + i1010_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i1010_fix], mm0
	movd [esp + i1010_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i1010_innerk],  2
	jl    .i1010_finish_coul_inner
	jmp   .i1010_unroll_coul_loop
.i1010_finish_coul_inner:	
	and dword ptr [esp + i1010_innerk],  1
	jnz  .i1010_single_coul_inner
	jmp  .i1010_updateouterdata_coul		
.i1010_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i1010_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i1010_charge]
	movd mm6, [esp + i1010_iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;# mm6=qq 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + i1010_ix]
	movd  mm1, [esp + i1010_iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + i1010_dx1], mm0
	pfmul mm0,mm0
	movd  [esp + i1010_dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;# mm0=rsq 
	
    pfrsqrt mm1,mm0
    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  
	;# update vctot 
	pfadd mm6, [esp + i1010_vctot]
	movq [esp + i1010_vctot], mm6
	;# spread fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + i1010_dx1]
	movd mm1,  [esp + i1010_dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	;# update i particle force 
	movq mm2,  [esp + i1010_fix]
	movd mm3,  [esp + i1010_fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1010_fix], mm2
	movd [esp + i1010_fiz], mm3
	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	;# done! 
.i1010_updateouterdata_coul:	
	mov   ecx, [esp + i1010_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1010_fix]
	pfadd mm7, [esp + i1010_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i1010_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1010_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1010_fix]
	pfadd mm7, [esp + i1010_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i1010_nscoul]
	jz  .i1010_last_mno
	jmp .i1010_mno_coul
.i1010_last_mno:	
	mov   edx, [ebp + i1010_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i1010_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i1010_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i1010_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + i1010_nri]
	dec ecx
	jecxz .i1010_end
	;# not last, iterate once more! 
	mov [ebp + i1010_nri], ecx
	jmp .i1010_outer
.i1010_end:
	femms
	add esp, 108
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret
	

			
.globl inl1020_3dnow
.globl _inl1020_3dnow
inl1020_3dnow:	
_inl1020_3dnow:	
.equiv		i1020_nri,		8
.equiv		i1020_iinr,		12
.equiv		i1020_jindex,		16
.equiv		i1020_jjnr,		20
.equiv		i1020_shift,		24
.equiv		i1020_shiftvec,		28
.equiv		i1020_fshift,		32
.equiv		i1020_gid,		36
.equiv		i1020_pos,		40		
.equiv		i1020_faction,		44
.equiv		i1020_charge,		48
.equiv		i1020_facel,		52
.equiv		i1020_Vc,		56			
			;# stack offsets for local variables 
.equiv		i1020_is3,		0
.equiv		i1020_ii3,		4
.equiv		i1020_ixO,		8
.equiv		i1020_iyO,		12
.equiv		i1020_izO,		16	
.equiv		i1020_ixH,		20 
.equiv		i1020_iyH,		28 
.equiv		i1020_izH,		36 
.equiv		i1020_iqO,		44	
.equiv		i1020_iqH,		52		
.equiv		i1020_vctot,		60 
.equiv		i1020_innerjjnr,	68
.equiv		i1020_innerk,		72		
.equiv		i1020_fixO,		76 
.equiv		i1020_fiyO,		80
.equiv		i1020_fizO,		84
.equiv		i1020_fixH,		88
.equiv		i1020_fiyH,		96
.equiv		i1020_fizH,		104         
.equiv		i1020_dxO,		112
.equiv		i1020_dyO,		116
.equiv		i1020_dzO,		120
.equiv		i1020_dxH,		124         
.equiv		i1020_dyH,		132         
.equiv		i1020_dzH,		140         
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 148		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + i1020_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i1020_charge]
	movd  mm1, [ebp + i1020_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + i1020_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i1020_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 
.i1020_outer:
	mov   eax, [ebp + i1020_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i1020_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i1020_is3],ebx    	;# store is3 

	mov   eax, [ebp + i1020_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i1020_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i1020_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1020_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i1020_ii3], ebx	    ;# (use mm7 as temp storage for iz) 
	pfadd mm6, mm7
	movq  [esp + i1020_ixO], mm5	
	movq  [esp + i1020_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i1020_ixH], mm0	
	movq [esp + i1020_iyH], mm1	
	movq [esp + i1020_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i1020_vctot], mm7
	movq  [esp + i1020_fixO],   mm7
	movd  [esp + i1020_fizO],   mm7
	movq  [esp + i1020_fixH],   mm7
	movq  [esp + i1020_fiyH],   mm7
	movq  [esp + i1020_fizH],   mm7

	mov   eax, [ebp + i1020_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i1020_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i1020_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + i1020_pos]
	mov   edi, [ebp + i1020_faction]	
	mov   eax, [ebp + i1020_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i1020_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i1020_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i1020_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i1020_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + i1020_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + i1020_iqO]
	pfmul mm7, [esp + i1020_iqH]	;# mm6=qqO, mm7=qqH 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i1020_ixO]
	pfsubr mm1, [esp + i1020_izO]
		
	movq  [esp + i1020_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i1020_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i1020_ixH]
	pfsubr mm3, [esp + i1020_iyH]
	pfsubr mm4, [esp + i1020_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i1020_dxH], mm2
	movq [esp + i1020_dyH], mm3
	movq [esp + i1020_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

	pfrsqrt mm1,mm0

	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + i1020_vctot]
	movq [esp + i1020_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + i1020_dxO]
	movd mm1,  [esp + i1020_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + i1020_dxH]
	movq mm6, [esp + i1020_dyH]
	movq mm7, [esp + i1020_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + i1020_fixO]
	movd mm3,  [esp + i1020_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1020_fixO], mm2
	movd [esp + i1020_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i1020_fixH]
	movq mm3, [esp + i1020_fiyH]
	movq mm4, [esp + i1020_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i1020_fixH], mm2
	movq [esp + i1020_fiyH], mm3
	movq [esp + i1020_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i1020_innerk]
	jz  .i1020_updateouterdata
	jmp .i1020_inner_loop
.i1020_updateouterdata:	
	mov   ecx, [esp + i1020_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1020_fixO]
	pfadd mm7, [esp + i1020_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i1020_fixH]
	movq  mm3, [esp + i1020_fiyH]
	movq  mm1, [esp + i1020_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 
	
	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i1020_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1020_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1020_fixO]
	pfadd mm7, [esp + i1020_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i1020_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i1020_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i1020_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i1020_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	;# finish if last 
	dec dword ptr [ebp + i1020_nri]
	jz  .i1020_end
	;# not last, iterate once more! 
	jmp .i1020_outer
.i1020_end:
	femms
	add esp, 148
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl inl1030_3dnow
.globl _inl1030_3dnow
inl1030_3dnow:	
_inl1030_3dnow:	
.equiv		i1030_nri,		8
.equiv		i1030_iinr,		12
.equiv		i1030_jindex,		16
.equiv		i1030_jjnr,		20
.equiv		i1030_shift,		24
.equiv		i1030_shiftvec,		28
.equiv		i1030_fshift,		32
.equiv		i1030_gid,		36
.equiv		i1030_pos,		40		
.equiv		i1030_faction,		44
.equiv		i1030_charge,		48
.equiv		i1030_facel,		52
.equiv		i1030_Vc,		56
			;# stack offsets for local variables 
.equiv		i1030_is3,		0
.equiv		i1030_ii3,		4
.equiv		i1030_ixO,		8
.equiv		i1030_iyO,		12
.equiv		i1030_izO,		16	
.equiv		i1030_ixH,		20
.equiv		i1030_iyH,		28
.equiv		i1030_izH,		36
.equiv		i1030_qqOO,		44		
.equiv		i1030_qqOH,		52		
.equiv		i1030_qqHH,		60     	
.equiv		i1030_vctot,		68
.equiv		i1030_innerjjnr,	76
.equiv		i1030_innerk,		80		
.equiv		i1030_fixO,		84 
.equiv		i1030_fiyO,		88
.equiv		i1030_fizO,		92
.equiv		i1030_fixH,		96
.equiv		i1030_fiyH,		104         
.equiv		i1030_fizH,		112         
.equiv		i1030_dxO,		120
.equiv		i1030_dyO,		124
.equiv		i1030_dzO,		128
.equiv		i1030_dxH,		132         
.equiv		i1030_dyH,		140         
.equiv		i1030_dzH,		148         
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 156		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + i1030_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i1030_charge]
	movd  mm1, [ebp + i1030_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + i1030_qqOO], mm4
	movq  [esp + i1030_qqOH], mm5
	movq  [esp + i1030_qqHH], mm6
.i1030_outer:
	mov   eax, [ebp + i1030_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i1030_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i1030_is3],ebx    	;# store is3 

	mov   eax, [ebp + i1030_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i1030_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i1030_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1030_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i1030_ii3], ebx	    ;# (use mm7 as temp storage for iz) 
	pfadd mm6, mm7
	movq  [esp + i1030_ixO], mm5	
	movq  [esp + i1030_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i1030_ixH], mm0	
	movq [esp + i1030_iyH], mm1	
	movq [esp + i1030_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i1030_vctot], mm7
	movq  [esp + i1030_fixO],  mm7
	movq  [esp + i1030_fizO],  mm7
	movq  [esp + i1030_fixH],  mm7
	movq  [esp + i1030_fiyH],  mm7
	movq  [esp + i1030_fizH],  mm7

	mov   eax, [ebp + i1030_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i1030_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i1030_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + i1030_pos]
	mov   edi, [ebp + i1030_faction]	
	mov   eax, [ebp + i1030_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i1030_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i1030_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i1030_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i1030_innerjjnr],  4 ;# advance pointer 

	movd  mm6, [esp + i1030_qqOO]
	movq  mm7, [esp + i1030_qqOH]

	lea   eax, [eax + eax*2]
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i1030_ixO]
	pfsubr mm1, [esp + i1030_izO]
		
	movq  [esp + i1030_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i1030_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i1030_ixH]
	pfsubr mm3, [esp + i1030_iyH]
	pfsubr mm4, [esp + i1030_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i1030_dxH], mm2
	movq [esp + i1030_dyH], mm3
	movq [esp + i1030_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    pfrsqit1 mm5,mm3				
    pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + i1030_vctot]
	movq [esp + i1030_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + i1030_dxO]
	movd mm1,  [esp + i1030_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + i1030_dxH]
	movq mm6, [esp + i1030_dyH]
	movq mm7, [esp + i1030_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + i1030_fixO]
	movd mm3,  [esp + i1030_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1030_fixO], mm2
	movd [esp + i1030_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i1030_fixH]
	movq mm3, [esp + i1030_fiyH]
	movq mm4, [esp + i1030_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i1030_fixH], mm2
	movq [esp + i1030_fiyH], mm3
	movq [esp + i1030_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	;# interactions with j H1 
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	movd mm6, [esp + i1030_qqOH]
	movq mm7, [esp + i1030_qqHH]
	
	pfsubr mm0, [esp + i1030_ixO]
	pfsubr mm1, [esp + i1030_izO]
		
	movq  [esp + i1030_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i1030_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i1030_ixH]
	pfsubr mm3, [esp + i1030_iyH]
	pfsubr mm4, [esp + i1030_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i1030_dxH], mm2
	movq [esp + i1030_dyH], mm3
	movq [esp + i1030_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

	pfrsqrt mm1,mm0

	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + i1030_vctot]
	movq [esp + i1030_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + i1030_dxO]
	movd mm1,  [esp + i1030_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + i1030_dxH]
	movq mm6, [esp + i1030_dyH]
	movq mm7, [esp + i1030_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + i1030_fixO]
	movd mm3,  [esp + i1030_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1030_fixO], mm2
	movd [esp + i1030_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i1030_fixH]
	movq mm3, [esp + i1030_fiyH]
	movq mm4, [esp + i1030_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i1030_fixH], mm2
	movq [esp + i1030_fiyH], mm3
	movq [esp + i1030_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	movd mm6, [esp + i1030_qqOH]
	movq mm7, [esp + i1030_qqHH]

	pfsubr mm0, [esp + i1030_ixO]
	pfsubr mm1, [esp + i1030_izO]
		
	movq  [esp + i1030_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i1030_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i1030_ixH]
	pfsubr mm3, [esp + i1030_iyH]
	pfsubr mm4, [esp + i1030_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i1030_dxH], mm2
	movq [esp + i1030_dyH], mm3
	movq [esp + i1030_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

	pfrsqrt mm1,mm0

	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + i1030_vctot]
	movq [esp + i1030_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + i1030_dxO]
	movd mm1,  [esp + i1030_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + i1030_dxH]
	movq mm6, [esp + i1030_dyH]
	movq mm7, [esp + i1030_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + i1030_fixO]
	movd mm3,  [esp + i1030_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1030_fixO], mm2
	movd [esp + i1030_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i1030_fixH]
	movq mm3, [esp + i1030_fiyH]
	movq mm4, [esp + i1030_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i1030_fixH], mm2
	movq [esp + i1030_fiyH], mm3
	movq [esp + i1030_fizH], mm4	

	;# pack j forces from H in the same form as the oxygen force 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i1030_innerk]
	jz  .i1030_updateouterdata
	jmp .i1030_inner_loop	
.i1030_updateouterdata:	
	mov   ecx, [esp + i1030_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1030_fixO]
	pfadd mm7, [esp + i1030_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i1030_fixH]
	movq  mm3, [esp + i1030_fiyH]
	movq  mm1, [esp + i1030_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 

	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i1030_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1030_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1030_fixO]
	pfadd mm7, [esp + i1030_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i1030_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i1030_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i1030_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + i1030_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	dec dword ptr [ebp + i1030_nri]
	jz  .i1030_end
	;# not last, iterate once more! 
	jmp .i1030_outer
.i1030_end:
	femms
	add esp, 156
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl inl1100_3dnow
.globl _inl1100_3dnow
inl1100_3dnow:	
_inl1100_3dnow:	
.equiv		i1100_nri,		8
.equiv		i1100_iinr,		12
.equiv		i1100_jindex,		16
.equiv		i1100_jjnr,		20
.equiv		i1100_shift,		24
.equiv		i1100_shiftvec,		28
.equiv		i1100_fshift,		32
.equiv		i1100_gid,		36
.equiv		i1100_pos,		40		
.equiv		i1100_faction,		44
.equiv		i1100_charge,		48
.equiv		i1100_facel,		52
.equiv		i1100_Vc,		56			
.equiv		i1100_type,		60
.equiv		i1100_ntype,		64
.equiv		i1100_nbfp,		68	
.equiv		i1100_Vnb,		72	
	;# stack offsets for local variables 
.equiv		i1100_is3,		0
.equiv		i1100_ii3,		4
.equiv		i1100_ix,		8
.equiv		i1100_iy,		12
.equiv		i1100_iz,		16
.equiv		i1100_iq,		20 
.equiv		i1100_vctot,		28 
.equiv		i1100_vnbtot,		36 
.equiv		i1100_c6,		44 
.equiv		i1100_c12,		52 
.equiv		i1100_six,		60 
.equiv		i1100_twelve,		68 
.equiv		i1100_ntia,		76
.equiv		i1100_innerjjnr,	80
.equiv		i1100_innerk,		84		
.equiv		i1100_fix,		88
.equiv		i1100_fiy,		92
.equiv		i1100_fiz,		96
.equiv		i1100_dx1,		100
.equiv		i1100_dy1,		104
.equiv		i1100_dz1,		108
.equiv		i1100_dx2,		112
.equiv		i1100_dy2,		116
.equiv		i1100_dz2,		120			
	push ebp
	mov ebp,esp
	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 124		;# local stack space 
	femms
	;# move data to local stack  
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + i1100_six],    mm0
	movq  [esp + i1100_twelve], mm1
	;# assume we have at least one i particle - start directly 	
.i1100_outer:
	mov   eax, [ebp + i1100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i1030_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i1030_is3],ebx    	;# store is3 

	mov   eax, [ebp + i1100_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + i1100_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i1100_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i1100_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i1100_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i1100_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + i1100_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i1100_ntype]
	shl   edx, 1
	mov   [esp + i1100_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1100_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i1100_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + i1100_ix], mm0	
	movd  [esp + i1100_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + i1100_vctot],  mm7
	movq  [esp + i1100_vnbtot], mm7
	movq  [esp + i1100_fix],    mm7
	movd  [esp + i1100_fiz],    mm7

	mov   eax, [ebp + i1100_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i1100_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + i1100_pos]
	mov   edi, [ebp + i1100_faction]	
	mov   eax, [ebp + i1100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i1100_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + i1100_innerk], edx    ;# number of innerloop atoms 
	jge   .i1100_unroll_loop
	jmp   .i1100_finish_inner
.i1100_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + i1100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i1100_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i1100_charge]    ;# base of charge[] 
	movq mm5, [esp + i1100_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + i1100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i1100_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i1100_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i1100_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i1100_c6], mm5
	movq [esp + i1100_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i1100_pos]

	movq  mm0, [esp + i1100_ix]
	movd  mm1, [esp + i1100_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i1100_dx1], mm4	     ;# store dr 
	movd  [esp + i1100_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i1100_dx2], mm6	     ;# store dr 
	movd  [esp + i1100_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	movq  mm7, mm3	    ;# use mm7 for sum to make fscal  

	pfmul mm5, [esp + i1100_c12]
	pfmul mm4, [esp + i1100_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i1100_six]

	pfmul mm5, [esp + i1100_twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7    ;# mm0 is total fscal now 	

	prefetchw [esp + i1100_dx1]	;# prefetch i forces to cache 

	;# update vctot 
	pfadd mm3, [esp + i1100_vctot]      ;# add the earlier value 
	movq [esp + i1100_vctot], mm3       ;# store the sum       

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i1100_dx1]	;# fetch dr 
	movd mm3,  [esp + i1100_dz1]

	;# update vnbtot 
	pfadd mm6, [esp + i1100_vnbtot]      ;# add the earlier value 
	movq [esp + i1100_vnbtot], mm6       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i1100_dx2] 	;# fetch dr 
	movd mm5,  [esp + i1100_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i1100_fix]
	movd mm1,  [esp + i1100_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i1100_fix], mm0
	movd [esp + i1100_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i1100_innerk],  2
	jl    .i1100_finish_inner
	jmp   .i1100_unroll_loop
.i1100_finish_inner:	
	and dword ptr [esp + i1100_innerk],  1
	jnz  .i1100_single_inner
	jmp  .i1100_updateouterdata		
.i1100_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i1100_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i1100_charge]
	movd mm5, [esp + i1100_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + i1100_nbfp]
	mov ecx, [ebp + i1100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i1100_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i1100_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i1100_c12], mm5


	mov   esi, [ebp + i1100_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i1100_ix]
	movd  mm1, [esp + i1100_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i1100_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i1100_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	movq  mm1, mm0
	pfmul mm0, mm0		;# mm0=invsq 
	;# calculate potentials and scalar force 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	movq  mm7, mm3	    ;# use mm7 for sum to make fscal  

	pfmul mm5, [esp + i1100_c12]
	pfmul mm4, [esp + i1100_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i1100_six]

	pfmul mm5, [esp + i1100_twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7    ;# mm0 is total fscal now 

	;# update vctot 
	pfadd mm3, [esp + i1100_vctot]
	movq [esp + i1100_vctot], mm3

	;# update vnbtot 
	pfadd mm6, [esp + i1100_vnbtot]      ;# add the earlier value 
	movq [esp + i1100_vnbtot], mm6       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i1100_dx1]
	movd mm3,  [esp + i1100_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i1100_fix]
	movd mm1,  [esp + i1100_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i1100_fix], mm0
	movd [esp + i1100_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i1100_updateouterdata:	
	mov   ecx, [esp + i1100_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1100_fix]
	pfadd mm7, [esp + i1100_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i1100_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1100_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1100_fix]
	pfadd mm7, [esp + i1100_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + i1100_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i1100_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i1100_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i1100_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i1100_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i1100_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + i1100_nri]
	dec ecx
	jecxz .i1100_end
	;# not last, iterate once more! 
	mov [ebp + i1100_nri], ecx
	jmp .i1100_outer
.i1100_end:
	femms
	add esp, 124
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	



.globl inl1110_3dnow
.globl _inl1110_3dnow
inl1110_3dnow:	
_inl1110_3dnow:	
.equiv		i1110_nri,		8
.equiv		i1110_iinr,		12
.equiv		i1110_jindex,		16
.equiv		i1110_jjnr,		20
.equiv		i1110_shift,		24
.equiv		i1110_shiftvec,		28
.equiv		i1110_fshift,		32
.equiv		i1110_gid,		36
.equiv		i1110_pos,		40		
.equiv		i1110_faction,		44
.equiv		i1110_charge,		48
.equiv		i1110_facel,		52
.equiv		i1110_Vc,		56	
.equiv		i1110_type,		60
.equiv		i1110_ntype,		64
.equiv		i1110_nbfp,		68	
.equiv		i1110_Vnb,		72				
.equiv		i1110_nsatoms,		76		
	;# stack offsets for local variables 
.equiv		i1110_is3,		0
.equiv		i1110_ii3,		4
.equiv		i1110_shX,		8
.equiv		i1110_shY,		12 
.equiv		i1110_shZ,		16	
.equiv		i1110_ix,		20
.equiv		i1110_iy,		24
.equiv		i1110_iz,		28	
.equiv		i1110_iq,		32		 
.equiv		i1110_vctot,		40 
.equiv		i1110_vnbtot,		48 
.equiv		i1110_c6,		56 
.equiv		i1110_c12,		64 
.equiv		i1110_six,		72 
.equiv		i1110_twelve,		80 
.equiv		i1110_ntia,		88	
.equiv		i1110_innerjjnr0,	92
.equiv		i1110_innerk0,		96		
.equiv		i1110_innerjjnr,	100
.equiv		i1110_innerk,		104	
.equiv		i1110_fix,		108
.equiv		i1110_fiy,		112
.equiv		i1110_fiz,		116
.equiv		i1110_dx1,		120
.equiv		i1110_dy1,		124
.equiv		i1110_dz1,		128
.equiv		i1110_dx2,		132
.equiv		i1110_dy2,		136
.equiv		i1110_dz2,		140
.equiv		i1110_nsvdwc,		144
.equiv		i1110_nscoul,		148
.equiv		i1110_nsvdw,		152
.equiv		i1110_solnr,		156		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 160		;# local stack space 
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + i1110_six],    mm0
	movq  [esp + i1110_twelve], mm1
	;# assume we have at least one i particle - start directly 		
.i1110_outer:
	mov   eax, [ebp + i1110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i1110_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i1110_is3],ebx    	;# store is3 

	mov   eax, [ebp + i1110_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + i1110_shX], mm0
	movd  [esp + i1110_shZ], mm1

	mov   ecx, [ebp + i1110_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i1110_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + i1110_nsatoms]
	add dword ptr [ebp + i1110_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + i1110_nsvdwc], edx
	mov   [esp + i1110_nscoul], eax
	mov   [esp + i1110_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + i1110_vctot],  mm7
	movq  [esp + i1110_vnbtot], mm7
	mov   [esp + i1110_solnr],  ebx
	
	mov   eax, [ebp + i1110_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i1110_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + i1110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i1110_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + i1110_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + i1110_pos]
	mov   edi, [ebp + i1110_faction]
	
	mov   ecx, [esp + i1110_nsvdwc]
	cmp   ecx,  0
	jnz   .i1110_mno_vdwc
	jmp   .i1110_testcoul
.i1110_mno_vdwc:
	mov   ebx,  [esp + i1110_solnr]
	inc   dword ptr [esp + i1110_solnr]
	mov   edx, [ebp + i1110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i1110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i1110_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + i1110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i1110_ntype]
	shl   edx, 1
	mov   [esp + i1110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1110_pos]    ;# eax = base of pos[] 
	mov   [esp + i1110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i1110_shX]
	pfadd mm1, [esp + i1110_shZ]
	movq  [esp + i1110_ix], mm0	
	movd  [esp + i1110_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i1110_fix],   mm7
	movd  [esp + i1110_fiz],   mm7

	mov   ecx, [esp + i1110_innerjjnr0]
	mov   [esp + i1110_innerjjnr], ecx
	mov   edx, [esp + i1110_innerk0]
    sub   edx,  2
    mov   [esp + i1110_innerk], edx    ;# number of innerloop atoms 
	jge   .i1110_unroll_vdwc_loop
	jmp   .i1110_finish_vdwc_inner
.i1110_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i1110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i1110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i1110_charge]    ;# base of charge[] 
	movq mm5, [esp + i1110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + i1110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i1110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i1110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i1110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i1110_c6], mm5
	movq [esp + i1110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i1110_pos]

	movq  mm0, [esp + i1110_ix]
	movd  mm1, [esp + i1110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i1110_dx1], mm4	     ;# store dr 
	movd  [esp + i1110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i1110_dx2], mm6	     ;# store dr 
	movd  [esp + i1110_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	movq  mm7, mm3	    ;# use mm7 for sum to make fscal  

	pfmul mm5, [esp + i1110_c12]
	pfmul mm4, [esp + i1110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i1110_six]

	pfmul mm5, [esp + i1110_twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7    ;# mm0 is total fscal now 	

	prefetchw [esp + i1110_dx1]	;# prefetch i forces to cache 

	;# update vctot 
	pfadd mm3, [esp + i1110_vctot]      ;# add the earlier value 
	movq [esp + i1110_vctot], mm3       ;# store the sum       

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i1110_dx1]	;# fetch dr 
	movd mm3,  [esp + i1110_dz1]

	;# update vnbtot 
	pfadd mm6, [esp + i1110_vnbtot]      ;# add the earlier value 
	movq [esp + i1110_vnbtot], mm6       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i1110_dx2] 	;# fetch dr 
	movd mm5,  [esp + i1110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i1110_fix]
	movd mm1,  [esp + i1110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i1110_fix], mm0
	movd [esp + i1110_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i1110_innerk],  2
	jl    .i1110_finish_vdwc_inner
	jmp   .i1110_unroll_vdwc_loop
.i1110_finish_vdwc_inner:	
	and dword ptr [esp + i1110_innerk],  1
	jnz  .i1110_single_vdwc_inner
	jmp  .i1110_updateouterdata_vdwc		
.i1110_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i1110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i1110_charge]
	movd mm5, [esp + i1110_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + i1110_nbfp]
	mov ecx, [ebp + i1110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i1110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i1110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i1110_c12], mm5


	mov   esi, [ebp + i1110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i1110_ix]
	movd  mm1, [esp + i1110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i1110_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i1110_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
	pfrsqrt mm0,mm4
	movq mm2,mm0
	pfmul mm0,mm0
	pfrsqit1 mm0,mm4				
	pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	movq  mm1, mm0
	pfmul mm0, mm0		;# mm0=invsq 
	;# calculate potentials and scalar force 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	movq  mm7, mm3	    ;# use mm7 for sum to make fscal  

	pfmul mm5, [esp + i1110_c12]
	pfmul mm4, [esp + i1110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i1110_six]

	pfmul mm5, [esp + i1110_twelve]
	pfsub mm7,mm4
	pfadd mm7, mm5
 	pfmul mm0, mm7    ;# mm0 is total fscal now 

	;# update vctot 
	pfadd mm3, [esp + i1110_vctot]
	movq [esp + i1110_vctot], mm3

	;# update vnbtot 
	pfadd mm6, [esp + i1110_vnbtot]      ;# add the earlier value 
	movq [esp + i1110_vnbtot], mm6       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i1110_dx1]
	movd mm3,  [esp + i1110_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i1110_fix]
	movd mm1,  [esp + i1110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i1110_fix], mm0
	movd [esp + i1110_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i1110_updateouterdata_vdwc:	
	mov   ecx, [esp + i1110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1110_fix]
	pfadd mm7, [esp + i1110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i1110_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1110_fix]
	pfadd mm7, [esp + i1110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec  dword ptr [esp + i1110_nsvdwc]
	jz  .i1110_testcoul
	jmp .i1110_mno_vdwc
.i1110_testcoul:	
	mov  ecx, [esp + i1110_nscoul]
	cmp  ecx,  0
	jnz  .i1110_mno_coul
	jmp  .i1110_testvdw
.i1110_mno_coul:
	mov   ebx,  [esp + i1110_solnr]
	inc   dword ptr [esp + i1110_solnr]
	mov   edx, [ebp + i1110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i1110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i1110_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1110_pos]    ;# eax = base of pos[] 
	mov   [esp + i1110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i1110_shX]
	pfadd mm1, [esp + i1110_shZ]
	movq  [esp + i1110_ix], mm0	
	movd  [esp + i1110_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i1110_fix],   mm7
	movd  [esp + i1110_fiz],   mm7

	mov   ecx, [esp + i1110_innerjjnr0]
	mov   [esp + i1110_innerjjnr], ecx
	mov   edx, [esp + i1110_innerk0]
    sub   edx,  2
    mov   [esp + i1110_innerk], edx    ;# number of innerloop atoms 
	jge   .i1110_unroll_coul_loop
	jmp   .i1110_finish_coul_inner
.i1110_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i1110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i1110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i1110_charge]    ;# base of charge[] 
	movq mm5, [esp + i1110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
	movd mm7, [ecx + ebx*4]  	 ;# charge[jnr2] 
	punpckldq mm3,mm7	     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + i1110_ix]
	movd  mm1, [esp + i1110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i1110_dx1], mm4	     ;# store dr 
	movd  [esp + i1110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i1110_dx2], mm6	     ;# store dr 
	movd  [esp + i1110_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

	pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
	pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
	movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
	pfrsqit1 mm0,mm4				
	pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	prefetchw [esp + i1110_dx1]	;# prefetch i forces to cache 
	
	pfmul mm3,mm1		;# 3 has both vcoul 
	pfmul mm0,mm3		;# 0 has both fscal 

	;# update vctot

	pfadd mm3, [esp + i1110_vctot]      ;# add the earlier value  
	movq [esp + i1110_vctot], mm3       ;# store the sum 
	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1
	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i1110_dx1]	;# fetch dr 
	movd mm3,  [esp + i1110_dz1]
	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i1110_dx2] 	;# fetch dr 
	movd mm5,  [esp + i1110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i1110_fix]
	movd mm1,  [esp + i1110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i1110_fix], mm0
	movd [esp + i1110_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i1110_innerk],  2
	jl    .i1110_finish_coul_inner
	jmp   .i1110_unroll_coul_loop
.i1110_finish_coul_inner:	
	and dword ptr [esp + i1110_innerk],  1
	jnz  .i1110_single_coul_inner
	jmp  .i1110_updateouterdata_coul		
.i1110_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i1110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i1110_charge]
	movd mm6, [esp + i1110_iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;# mm6=qq 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + i1110_ix]
	movd  mm1, [esp + i1110_iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq  [esp + i1110_dx1], mm0
	pfmul mm0,mm0
	movd  [esp + i1110_dz1], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;# mm0=rsq 
	
    pfrsqrt mm1,mm0 
    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  
	;# update vctot 
	pfadd mm6, [esp + i1110_vctot]
	movq [esp + i1110_vctot], mm6
	;# spread fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + i1110_dx1]
	movd mm1,  [esp + i1110_dz1]
	pfmul mm0, mm4
	pfmul mm1, mm4
	;# update i particle force 
	movq mm2,  [esp + i1110_fix]
	movd mm3,  [esp + i1110_fiz]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1110_fix], mm2
	movd [esp + i1110_fiz], mm3
	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax *4+ 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	;# done! 
.i1110_updateouterdata_coul:	
	mov   ecx, [esp + i1110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1110_fix]
	pfadd mm7, [esp + i1110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i1110_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1110_fix]
	pfadd mm7, [esp + i1110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i1110_nscoul]
	jz  .i1110_testvdw
	jmp .i1110_mno_coul
.i1110_testvdw:	
	mov  ecx, [esp + i1110_nsvdw]
	cmp  ecx,  0
	jnz  .i1110_mno_vdw
	jmp  .i1110_last_mno
.i1110_mno_vdw:
	mov   ebx,  [esp + i1110_solnr]
	inc   dword ptr [esp + i1110_solnr]

	mov   edx, [ebp + i1110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i1110_ntype]
	shl   edx, 1
	mov   [esp + i1110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1110_pos]    ;# eax = base of pos[] 
	mov   [esp + i1110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i1110_shX]
	pfadd mm1, [esp + i1110_shZ]
	movq  [esp + i1110_ix], mm0	
	movd  [esp + i1110_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i1110_fix],   mm7
	movd  [esp + i1110_fiz],   mm7

	mov   ecx, [esp + i1110_innerjjnr0]
	mov   [esp + i1110_innerjjnr], ecx
	mov   edx, [esp + i1110_innerk0]
    sub   edx,  2
    mov   [esp + i1110_innerk], edx    ;# number of innerloop atoms 
	jge   .i1110_unroll_vdw_loop
	jmp   .i1110_finish_vdw_inner
.i1110_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i1110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i1110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + i1110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i1110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i1110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i1110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i1110_c6], mm5
	movq [esp + i1110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i1110_pos]

	movq  mm0, [esp + i1110_ix]
	movd  mm1, [esp + i1110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i1110_dx1], mm4	     ;# store dr 
	movd  [esp + i1110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i1110_dx2], mm6	     ;# store dr 
	movd  [esp + i1110_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i1110_c12]
	pfmul mm4, [esp + i1110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i1110_six]

	pfmul mm5, [esp + i1110_twelve]
	movq  mm7, mm5
	pfsub mm7,mm4
 	pfmul mm0, mm7    ;# mm0 is total fscal now 	

	prefetchw [esp + i1110_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i1110_dx1]	;# fetch dr 
	movd mm3,  [esp + i1110_dz1]

	;# update vnbtot 
	pfadd mm6, [esp + i1110_vnbtot]      ;# add the earlier value 
	movq [esp + i1110_vnbtot], mm6       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i1110_dx2] 	;# fetch dr 
	movd mm5,  [esp + i1110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i1110_fix]
	movd mm1,  [esp + i1110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i1110_fix], mm0
	movd [esp + i1110_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i1110_innerk],  2
	jl    .i1110_finish_vdw_inner
	jmp   .i1110_unroll_vdw_loop
.i1110_finish_vdw_inner:	
	and dword ptr [esp + i1110_innerk],  1
	jnz  .i1110_single_vdw_inner
	jmp  .i1110_updateouterdata_vdw		
.i1110_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + i1110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + i1110_nbfp]
	mov ecx, [ebp + i1110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i1110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i1110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i1110_c12], mm5


	mov   esi, [ebp + i1110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i1110_ix]
	movd  mm1, [esp + i1110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i1110_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i1110_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
	pfrsqrt mm0,mm4
	movq mm2,mm0
	pfmul mm0,mm0
	pfrsqit1 mm0,mm4				
	pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	movq  mm1, mm0
	pfmul mm0, mm0		;# mm0=invsq 
	;# calculate potentials and scalar force 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i1110_c12]
	pfmul mm4, [esp + i1110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i1110_six]

	pfmul mm5, [esp + i1110_twelve]
	movq  mm7, mm5
	pfsub mm7,mm4
 	pfmul mm0, mm7    ;# mm0 is total fscal now 

	;# update vnbtot 
	pfadd mm6, [esp + i1110_vnbtot]      ;# add the earlier value 
	movq [esp + i1110_vnbtot], mm6       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i1110_dx1]
	movd mm3,  [esp + i1110_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i1110_fix]
	movd mm1,  [esp + i1110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i1110_fix], mm0
	movd [esp + i1110_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i1110_updateouterdata_vdw:	
	mov   ecx, [esp + i1110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1110_fix]
	pfadd mm7, [esp + i1110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i1110_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1110_fix]
	pfadd mm7, [esp + i1110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i1110_nsvdw]
	jz  .i1110_last_mno
	jmp .i1110_mno_vdw
	
.i1110_last_mno:	
	mov   edx, [ebp + i1110_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i1110_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i1110_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i1110_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i1110_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i1110_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + i1110_nri]
	dec ecx
	jecxz .i1110_end
	;# not last, iterate once more! 
	mov [ebp + i1110_nri], ecx
	jmp .i1110_outer
.i1110_end:
	femms
	add esp, 160
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret



.globl inl1120_3dnow
.globl _inl1120_3dnow
inl1120_3dnow:	
_inl1120_3dnow:	
.equiv		i1120_nri,		8
.equiv		i1120_iinr,		12
.equiv		i1120_jindex,		16
.equiv		i1120_jjnr,		20
.equiv		i1120_shift,		24
.equiv		i1120_shiftvec,		28
.equiv		i1120_fshift,		32
.equiv		i1120_gid,		36
.equiv		i1120_pos,		40		
.equiv		i1120_faction,		44
.equiv		i1120_charge,		48
.equiv		i1120_facel,		52
.equiv		i1120_Vc,		56	
.equiv		i1120_type,		60
.equiv		i1120_ntype,		64
.equiv		i1120_nbfp,		68	
.equiv		i1120_Vnb,		72				
			;# stack offsets for local variables 
.equiv		i1120_is3,		0
.equiv		i1120_ii3,		4
.equiv		i1120_ixO,		8
.equiv		i1120_iyO,		12
.equiv		i1120_izO,		16	
.equiv		i1120_ixH,		20  
.equiv		i1120_iyH,		28  
.equiv		i1120_izH,		36  
.equiv		i1120_iqO,		44  
.equiv		i1120_iqH,		52  
.equiv		i1120_vctot,		60  
.equiv		i1120_vnbtot,		68  
.equiv		i1120_c6,		76  
.equiv		i1120_c12,		84  
.equiv		i1120_six,		92  
.equiv		i1120_twelve,		100 
.equiv		i1120_ntia,		108 
.equiv		i1120_innerjjnr,	116
.equiv		i1120_innerk,		120	
.equiv		i1120_fixO,		124
.equiv		i1120_fiyO,		128
.equiv		i1120_fizO,		132
.equiv		i1120_fixH,		136  
.equiv		i1120_fiyH,		144  
.equiv		i1120_fizH,		152  
.equiv		i1120_dxO,		160
.equiv		i1120_dyO,		164
.equiv		i1120_dzO,		168
.equiv		i1120_dxH,		172  
.equiv		i1120_dyH,		180  
.equiv		i1120_dzH,		188  
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 196		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + i1120_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i1120_charge]
	movd  mm1, [ebp + i1120_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + i1120_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i1120_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + i1120_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [ebp + i1120_ntype]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + i1120_ntia], ecx
	
	movq  mm3, [mm_six]
	movq  mm4, [mm_twelve]
	movq  [esp + i1120_six],    mm3
	movq  [esp + i1120_twelve], mm4  
.i1120_outer:
	mov   eax, [ebp + i1120_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i1120_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i1120_is3],ebx    	;# store is3 

	mov   eax, [ebp + i1120_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i1120_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i1120_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1120_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i1120_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + i1120_ixO], mm5	
	movq  [esp + i1120_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i1120_ixH], mm0	
	movq [esp + i1120_iyH], mm1	
	movq [esp + i1120_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i1120_vctot], mm7
	movq  [esp + i1120_vnbtot], mm7
	movq  [esp + i1120_fixO],   mm7
	movd  [esp + i1120_fizO],   mm7
	movq  [esp + i1120_fixH],   mm7
	movq  [esp + i1120_fiyH],   mm7
	movq  [esp + i1120_fizH],   mm7

	mov   eax, [ebp + i1120_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i1120_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i1120_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + i1120_pos]
	mov   edi, [ebp + i1120_faction]	
	mov   eax, [ebp + i1120_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i1120_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i1120_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i1120_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i1120_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + i1120_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + i1120_iqO]
	pfmul mm7, [esp + i1120_iqH]	;# mm6=qqO, mm7=qqH 

	mov ecx, [ebp + i1120_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + i1120_nbfp]
	shl edx, 1
	add edx, [esp + i1120_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i1120_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i1120_c12], mm5	
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i1120_ixO]
	pfsubr mm1, [esp + i1120_izO]
		
	movq  [esp + i1120_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i1120_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i1120_ixH]
	pfsubr mm3, [esp + i1120_iyH]
	pfsubr mm4, [esp + i1120_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i1120_dxH], mm2
	movq [esp + i1120_dyH], mm3
	movq [esp + i1120_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 

	movq  mm0, mm4
	pfmul mm0, mm4
	pfmul mm0, mm4		;# mm0=rinvsix 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm2=rintwelve 
	
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	movq  mm1, mm6		;# use mm1 for fscal sum 

	;# LJ for the oxygen 
	pfmul mm0, [esp + i1120_c6]	 
	pfmul mm2, [esp + i1120_c12]	 

	;# calc nb potential 
	movq mm5, mm2
	pfsub mm5, mm0

	;# calc nb force 
	pfmul mm0, [esp + i1120_six]
	pfmul mm2, [esp + i1120_twelve]
	
	;# increment scalar force 
	pfsub mm1, mm0
	pfadd mm1, mm2
	pfmul mm4, mm1		;# total scalar force on oxygen. 
	
	;# update nb potential 
	pfadd mm5, [esp + i1120_vnbtot]
	movq [esp + i1120_vnbtot], mm5
	
	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3. 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's. 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + i1120_vctot]
	movq [esp + i1120_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + i1120_dxO]
	movd mm1,  [esp + i1120_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + i1120_dxH]
	movq mm6, [esp + i1120_dyH]
	movq mm7, [esp + i1120_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + i1120_fixO]
	movd mm3,  [esp + i1120_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1120_fixO], mm2
	movd [esp + i1120_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i1120_fixH]
	movq mm3, [esp + i1120_fiyH]
	movq mm4, [esp + i1120_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i1120_fixH], mm2
	movq [esp + i1120_fiyH], mm3
	movq [esp + i1120_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i1120_innerk]
	jz  .i1120_updateouterdata
	jmp .i1120_inner_loop
.i1120_updateouterdata:	
	mov   ecx, [esp + i1120_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1120_fixO]
	pfadd mm7, [esp + i1120_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i1120_fixH]
	movq  mm3, [esp + i1120_fiyH]
	movq  mm1, [esp + i1120_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 
	
	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i1120_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1120_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1120_fixO]
	pfadd mm7, [esp + i1120_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i1120_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i1120_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i1120_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i1120_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i1120_vnbtot]     
	pfacc mm7,mm7	          ;# same for Vnb 
	
	mov   eax, [ebp + i1120_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 
	;# finish if last 
	dec dword ptr [ebp + i1120_nri]
	jz  .i1120_end
	;# not last, iterate once more! 
	jmp .i1120_outer
.i1120_end:
	femms
	add esp, 196
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	

.globl inl1130_3dnow
.globl _inl1130_3dnow
inl1130_3dnow:	
_inl1130_3dnow:	
.equiv		i1130_nri,		8
.equiv		i1130_iinr,		12
.equiv		i1130_jindex,		16
.equiv		i1130_jjnr,		20
.equiv		i1130_shift,		24
.equiv		i1130_shiftvec,		28
.equiv		i1130_fshift,		32
.equiv		i1130_gid,		36
.equiv		i1130_pos,		40		
.equiv		i1130_faction,		44
.equiv		i1130_charge,		48
.equiv		i1130_facel,		52
.equiv		i1130_Vc,		56						
.equiv		i1130_type,		60
.equiv		i1130_ntype,		64
.equiv		i1130_nbfp,		68	
.equiv		i1130_Vnb,		72			
			;# stack offsets for local variables 
.equiv		i1130_is3,		0
.equiv		i1130_ii3,		4
.equiv		i1130_ixO,		8
.equiv		i1130_iyO,		12
.equiv		i1130_izO,		16	
.equiv		i1130_ixH,		20  
.equiv		i1130_iyH,		28  
.equiv		i1130_izH,		36  
.equiv		i1130_qqOO,		44  
.equiv		i1130_qqOH,		52  
.equiv		i1130_qqHH,		60  
.equiv		i1130_c6,		68  
.equiv		i1130_c12,		76  
.equiv		i1130_six,		84  
.equiv		i1130_twelve,		92  
.equiv		i1130_vctot,		100 
.equiv		i1130_vnbtot,		108 
.equiv		i1130_innerjjnr,	116
.equiv		i1130_innerk,		120	
.equiv		i1130_fixO,		124
.equiv		i1130_fiyO,		128
.equiv		i1130_fizO,		132
.equiv		i1130_fixH,		136 
.equiv		i1130_fiyH,		144 
.equiv		i1130_fizH,		152 
.equiv		i1130_dxO,		160
.equiv		i1130_dyO,		164
.equiv		i1130_dzO,		168
.equiv		i1130_dxH,		172 
.equiv		i1130_dyH,		180 
.equiv		i1130_dzH,		188 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 196		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + i1130_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i1130_charge]
	movd  mm1, [ebp + i1130_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + i1130_qqOO], mm4
	movq  [esp + i1130_qqOH], mm5
	movq  [esp + i1130_qqHH], mm6
	mov   edx, [ebp + i1130_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + i1130_ntype]
	add   edx, ecx
	mov   eax, [ebp + i1130_nbfp]
	movd  mm0, [eax + edx*4]          
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + i1130_c6], mm0
	movq  [esp + i1130_c12], mm1
	movq  mm2, [mm_six]
	movq  mm3, [mm_twelve]
	movq  [esp + i1130_six], mm2
	movq  [esp + i1130_twelve], mm3
.i1130_outer:
	mov   eax, [ebp + i1130_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i1130_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i1130_is3],ebx    	;# store is3 

	mov   eax, [ebp + i1130_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i1130_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i1130_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i1130_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i1130_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + i1130_ixO], mm5	
	movq  [esp + i1130_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i1130_ixH], mm0	
	movq [esp + i1130_iyH], mm1	
	movq [esp + i1130_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i1130_vctot], mm7
	movq  [esp + i1130_vnbtot], mm7
	movq  [esp + i1130_fixO],  mm7
	movq  [esp + i1130_fizO],  mm7
	movq  [esp + i1130_fixH],  mm7
	movq  [esp + i1130_fiyH],  mm7
	movq  [esp + i1130_fizH],  mm7

	mov   eax, [ebp + i1130_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i1130_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i1130_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + i1130_pos]
	mov   edi, [ebp + i1130_faction]	
	mov   eax, [ebp + i1130_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i1130_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i1130_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i1130_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i1130_innerjjnr],  4 ;# advance pointer 

	movd  mm6, [esp + i1130_qqOO]
	movq  mm7, [esp + i1130_qqOH]

	lea   eax, [eax + eax*2]
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i1130_ixO]
	pfsubr mm1, [esp + i1130_izO]
		
	movq  [esp + i1130_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i1130_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i1130_ixH]
	pfsubr mm3, [esp + i1130_iyH]
	pfsubr mm4, [esp + i1130_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i1130_dxH], mm2
	movq [esp + i1130_dyH], mm3
	movq [esp + i1130_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

        pfrsqrt mm1,mm0

        movq mm2,mm1
        pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq  

	movq mm2, mm4
	pfmul mm2, mm4
	pfmul mm2, mm4
	movq mm0, mm2
	pfmul mm0,mm0
	pfmul mm2, [esp + i1130_c6]
	pfmul mm0, [esp + i1130_c12]
	movq mm5, mm0
	pfsub mm5, mm2		;# vnb 

	pfmul mm2, [esp + i1130_six]
	pfmul mm0, [esp + i1130_twelve]

	pfsub mm0, mm2
	
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfadd mm0, mm6
	pfmul mm4, mm0		;# mm4=fscalar  

	;# update nb potential 
	pfadd mm5, [esp + i1130_vnbtot]
	movq [esp + i1130_vnbtot], mm5

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's. 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + i1130_vctot]
	movq [esp + i1130_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + i1130_dxO]
	movd mm1,  [esp + i1130_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + i1130_dxH]
	movq mm6, [esp + i1130_dyH]
	movq mm7, [esp + i1130_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + i1130_fixO]
	movd mm3,  [esp + i1130_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1130_fixO], mm2
	movd [esp + i1130_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i1130_fixH]
	movq mm3, [esp + i1130_fiyH]
	movq mm4, [esp + i1130_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i1130_fixH], mm2
	movq [esp + i1130_fiyH], mm3
	movq [esp + i1130_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	;# interactions with j H1 
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	movd mm6, [esp + i1130_qqOH]
	movq mm7, [esp + i1130_qqHH]
	
	pfsubr mm0, [esp + i1130_ixO]
	pfsubr mm1, [esp + i1130_izO]
		
	movq  [esp + i1130_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i1130_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i1130_ixH]
	pfsubr mm3, [esp + i1130_iyH]
	pfsubr mm4, [esp + i1130_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i1130_dxH], mm2
	movq [esp + i1130_dyH], mm3
	movq [esp + i1130_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

	pfrsqrt mm1,mm0

	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's. 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + i1130_vctot]
	movq [esp + i1130_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + i1130_dxO]
	movd mm1,  [esp + i1130_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + i1130_dxH]
	movq mm6, [esp + i1130_dyH]
	movq mm7, [esp + i1130_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + i1130_fixO]
	movd mm3,  [esp + i1130_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1130_fixO], mm2
	movd [esp + i1130_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i1130_fixH]
	movq mm3, [esp + i1130_fiyH]
	movq mm4, [esp + i1130_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i1130_fixH], mm2
	movq [esp + i1130_fiyH], mm3
	movq [esp + i1130_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	movd mm6, [esp + i1130_qqOH]
	movq mm7, [esp + i1130_qqHH]

	pfsubr mm0, [esp + i1130_ixO]
	pfsubr mm1, [esp + i1130_izO]
		
	movq  [esp + i1130_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i1130_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i1130_ixH]
	pfsubr mm3, [esp + i1130_iyH]
	pfsubr mm4, [esp + i1130_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i1130_dxH], mm2
	movq [esp + i1130_dyH], mm3
	movq [esp + i1130_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

	pfrsqrt mm1,mm0

	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfmul mm4, mm6		;# mm4=fscalar  

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3. 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	movq mm3,mm5
	pfmul mm3,mm3		;# mm3=invsq 
	pfmul mm7, mm5		;# mm7=vcoul 
	pfmul mm3, mm7		;# mm3=fscal for the two H's. 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + i1130_vctot]
	movq [esp + i1130_vctot], mm7
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm4,mm4
	;# calc vectorial force for O 
	movq mm0,  [esp + i1130_dxO]
	movd mm1,  [esp + i1130_dzO]
	pfmul mm0, mm4
	pfmul mm1, mm4

	;# calc vectorial force for H's 
	movq mm5, [esp + i1130_dxH]
	movq mm6, [esp + i1130_dyH]
	movq mm7, [esp + i1130_dzH]
	pfmul mm5, mm3
	pfmul mm6, mm3
	pfmul mm7, mm3
	
	;# update iO particle force 
	movq mm2,  [esp + i1130_fixO]
	movd mm3,  [esp + i1130_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i1130_fixO], mm2
	movd [esp + i1130_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i1130_fixH]
	movq mm3, [esp + i1130_fiyH]
	movq mm4, [esp + i1130_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i1130_fixH], mm2
	movq [esp + i1130_fiyH], mm3
	movq [esp + i1130_fizH], mm4	

	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i1130_innerk]
	jz  .i1130_updateouterdata
	jmp .i1130_inner_loop	
.i1130_updateouterdata:	
	mov   ecx, [esp + i1130_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i1130_fixO]
	pfadd mm7, [esp + i1130_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i1130_fixH]
	movq  mm3, [esp + i1130_fiyH]
	movq  mm1, [esp + i1130_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 

	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i1130_fshift]    ;# increment fshift force 
	mov   edx, [esp + i1130_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i1130_fixO]
	pfadd mm7, [esp + i1130_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i1130_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i1130_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i1130_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + i1130_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i1130_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + i1130_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnbtot[gid] 
	;# finish if last 
	dec dword ptr [ebp + i1130_nri]
	jz  .i1130_end
	;# not last, iterate once more! 
	jmp .i1130_outer
.i1130_end:
	femms
	add esp, 196
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret




.globl inl3000_3dnow
.globl _inl3000_3dnow
inl3000_3dnow:	
_inl3000_3dnow:	
.equiv		i3000_nri,		8
.equiv		i3000_iinr,		12
.equiv		i3000_jindex,		16
.equiv		i3000_jjnr,		20
.equiv		i3000_shift,		24
.equiv		i3000_shiftvec,		28
.equiv		i3000_fshift,		32
.equiv		i3000_gid,		36
.equiv		i3000_pos,		40		
.equiv		i3000_faction,		44
.equiv		i3000_charge,		48
.equiv		i3000_facel,		52
.equiv		i3000_Vc,		56			
.equiv		i3000_tabscale,		60
.equiv		i3000_VFtab,		64
	;# stack offsets for local variables 
.equiv		i3000_is3,		0
.equiv		i3000_ii3,		4
.equiv		i3000_ix,		8
.equiv		i3000_iy,		12
.equiv		i3000_iz,		16
.equiv		i3000_iq,		20 
.equiv		i3000_vctot,		28 
.equiv		i3000_two,		36 
.equiv		i3000_n1,		44 
.equiv		i3000_tsc,		52 
.equiv		i3000_ntia,		60
.equiv		i3000_innerjjnr,	64
.equiv		i3000_innerk,		68		
.equiv		i3000_fix,		72
.equiv		i3000_fiy,		76
.equiv		i3000_fiz,		80
.equiv		i3000_dx1,		84
.equiv		i3000_dy1,		88
.equiv		i3000_dz1,		92
.equiv		i3000_dx2,		96
.equiv		i3000_dy2,		100
.equiv		i3000_dz2,		104
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 108		;# local stack space 
	femms
	;# move data to local stack  
	movq  mm0, [mm_two]
	movd  mm3, [ebp + i3000_tabscale]
	movq  [esp + i3000_two],    mm0
	punpckldq mm3,mm3
	movq  [esp + i3000_tsc], mm3	
	;# assume we have at least one i particle - start directly 	
.i3000_outer:
	mov   eax, [ebp + i3000_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3000_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3000_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3000_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + i3000_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3000_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3000_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i3000_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3000_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3000_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3000_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + i3000_ix], mm0	
	movd  [esp + i3000_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + i3000_vctot],  mm7
	movq  [esp + i3000_fix],    mm7
	movd  [esp + i3000_fiz],    mm7

	mov   eax, [ebp + i3000_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3000_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + i3000_pos]
	mov   edi, [ebp + i3000_faction]	
	mov   eax, [ebp + i3000_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3000_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + i3000_innerk], edx    ;# number of innerloop atoms 
	jge   .i3000_unroll_loop
	jmp   .i3000_finish_inner
.i3000_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + i3000_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3000_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3000_charge]    ;# base of charge[] 
	movq mm5, [esp + i3000_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3000_pos]

	movq  mm0, [esp + i3000_ix]
	movd  mm1, [esp + i3000_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3000_dx1], mm4	     ;# store dr 
	movd  [esp + i3000_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3000_dx2], mm6	     ;# store dr 
	movd  [esp + i3000_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3000_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i3000_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3000_VFtab]
	mov ecx, [esp + i3000_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3000_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3000_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC. 
	;# increment vcoul - then we can get rid of mm5. 
	;# update vctot 
	pfadd mm5, [esp + i3000_vctot]      ;# add the earlier value 
	movq [esp + i3000_vctot], mm5       ;# store the sum       

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + i3000_tsc]	
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3000_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3000_dx1]	;# fetch dr 
	movd mm3,  [esp + i3000_dz1]

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3000_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3000_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3000_fix]
	movd mm1,  [esp + i3000_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3000_fix], mm0
	movd [esp + i3000_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3000_innerk],  2
	jl    .i3000_finish_inner
	jmp   .i3000_unroll_loop
.i3000_finish_inner:	
	and dword ptr [esp + i3000_innerk],  1
	jnz  .i3000_single_inner
	jmp  .i3000_updateouterdata		
.i3000_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3000_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i3000_charge]
	movd mm5, [esp + i3000_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + i3000_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3000_ix]
	movd  mm1, [esp + i3000_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3000_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3000_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3000_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3000_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3000_VFtab]
	mov ecx, [esp + i3000_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3000_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3000_vctot]      ;# add the earlier value 
	movq [esp + i3000_vctot], mm5       ;# store the sum       
	
	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3000_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3000_dx1]
	movd mm3,  [esp + i3000_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3000_fix]
	movd mm1,  [esp + i3000_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3000_fix], mm0
	movd [esp + i3000_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3000_updateouterdata:	
	mov   ecx, [esp + i3000_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3000_fix]
	pfadd mm7, [esp + i3000_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3000_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3000_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3000_fix]
	pfadd mm7, [esp + i3000_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + i3000_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3000_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3000_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3000_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	;# finish if last 
	mov   ecx, [ebp + i3000_nri]
	dec ecx
	jecxz .i3000_end
	;# not last, iterate once more! 
	mov [ebp + i3000_nri], ecx
	jmp .i3000_outer
.i3000_end:
	femms
	add esp, 108
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret



	
.globl inl3010_3dnow
.globl _inl3010_3dnow
inl3010_3dnow:	
_inl3010_3dnow:	
.equiv		i3010_nri,		8
.equiv		i3010_iinr,		12
.equiv		i3010_jindex,		16
.equiv		i3010_jjnr,		20
.equiv		i3010_shift,		24
.equiv		i3010_shiftvec,		28
.equiv		i3010_fshift,		32
.equiv		i3010_gid,		36
.equiv		i3010_pos,		40		
.equiv		i3010_faction,		44
.equiv		i3010_charge,		48
.equiv		i3010_facel,		52
.equiv		i3010_Vc,		56
.equiv		i3010_tabscale,		60		
.equiv		i3010_VFtab,		64
.equiv		i3010_nsatoms,		68		
	;# stack offsets for local variables 
.equiv		i3010_is3,		0
.equiv		i3010_ii3,		4
.equiv		i3010_shX,		8
.equiv		i3010_shY,		12 
.equiv		i3010_shZ,		16	
.equiv		i3010_ix,		20
.equiv		i3010_iy,		24
.equiv		i3010_iz,		28	
.equiv		i3010_iq,		32 
.equiv		i3010_vctot,		40 
.equiv		i3010_two,		48 
.equiv		i3010_n1,		56 
.equiv		i3010_tsc,		64 			
.equiv		i3010_innerjjnr0,	72
.equiv		i3010_innerk0,		76		
.equiv		i3010_innerjjnr,	80
.equiv		i3010_innerk,		84		
.equiv		i3010_fix,		88
.equiv		i3010_fiy,		92
.equiv		i3010_fiz,		96
.equiv		i3010_dx1,		100
.equiv		i3010_dy1,		104
.equiv		i3010_dz1,		108
.equiv		i3010_dx2,		112
.equiv		i3010_dy2,		116
.equiv		i3010_dz2,		120
.equiv		i3010_nscoul,		124
.equiv		i3010_solnr,		128		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 132		;# local stack space 
	femms
	
	add dword ptr [ebp + i3010_nsatoms],  8
	movq  mm2, [mm_two]
	movq  [esp + i3010_two], mm2
	movd  mm3, [ebp + i3010_tabscale]
	punpckldq mm3,mm3
	movq  [esp + i3010_tsc], mm3
	
	;# assume we have at least one i particle - start directly 		
.i3010_outer:
	mov   eax, [ebp + i3010_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3010_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3010_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3010_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + i3010_shX], mm0
	movd  [esp + i3010_shZ], mm1

	mov   ecx, [ebp + i3010_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3010_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + i3010_nsatoms]
	mov   ecx, [eax]
	add dword ptr [ebp + i3010_nsatoms],  12
	mov   [esp + i3010_nscoul], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + i3010_vctot], mm7
	mov   [esp + i3010_solnr], ebx
	
	mov   eax, [ebp + i3010_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3010_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + i3010_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3010_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + i3010_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + i3010_pos]
	mov   edi, [ebp + i3010_faction]
	mov   ecx, [esp + i3010_nscoul]
	cmp   ecx,  0
	jnz  .i3010_mno_coul
	jmp  .i3010_last_mno
.i3010_mno_coul:				
	mov   ebx,  [esp + i3010_solnr]
	inc   dword ptr [esp + i3010_solnr]
	mov   edx, [ebp + i3010_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i3010_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3010_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3010_pos]    ;# eax = base of pos[] 
	mov   [esp + i3010_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i3010_shX]
	pfadd mm1, [esp + i3010_shZ]
	movq  [esp + i3010_ix], mm0	
	movd  [esp + i3010_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i3010_fix],   mm7
	movd  [esp + i3010_fiz],   mm7

	mov   ecx, [esp + i3010_innerjjnr0]
	mov   [esp + i3010_innerjjnr], ecx
	mov   edx, [esp + i3010_innerk0]
    sub   edx,  2
    mov   [esp + i3010_innerk], edx    ;# number of innerloop atoms 
	jge   .i3010_unroll_coul_loop
	jmp   .i3010_finish_coul_inner
.i3010_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i3010_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3010_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3010_charge]    ;# base of charge[] 
	movq mm5, [esp + i3010_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3010_pos]

	movq  mm0, [esp + i3010_ix]
	movd  mm1, [esp + i3010_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3010_dx1], mm4	     ;# store dr 
	movd  [esp + i3010_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3010_dx2], mm6	     ;# store dr 
	movd  [esp + i3010_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3010_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i3010_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3010_VFtab]
	mov ecx, [esp + i3010_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3010_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3010_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3010_vctot]      ;# add the earlier value 
	movq [esp + i3010_vctot], mm5       ;# store the sum       

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + i3010_tsc]	
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3010_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3010_dx1]	;# fetch dr 
	movd mm3,  [esp + i3010_dz1]

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3010_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3010_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3010_fix]
	movd mm1,  [esp + i3010_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3010_fix], mm0
	movd [esp + i3010_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3010_innerk],  2
	jl    .i3010_finish_coul_inner
	jmp   .i3010_unroll_coul_loop
.i3010_finish_coul_inner:	
	and dword ptr [esp + i3010_innerk],  1
	jnz  .i3010_single_coul_inner
	jmp  .i3010_updateouterdata_coul		
.i3010_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3010_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i3010_charge]
	movd mm5, [esp + i3010_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + i3010_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3010_ix]
	movd  mm1, [esp + i3010_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3010_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3010_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3010_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3010_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3010_VFtab]
	mov ecx, [esp + i3010_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3010_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3010_vctot]      ;# add the earlier value 
	movq [esp + i3010_vctot], mm5       ;# store the sum       
	
	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3010_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3010_dx1]
	movd mm3,  [esp + i3010_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3010_fix]
	movd mm1,  [esp + i3010_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3010_fix], mm0
	movd [esp + i3010_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3010_updateouterdata_coul:	
	mov   ecx, [esp + i3010_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3010_fix]
	pfadd mm7, [esp + i3010_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3010_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3010_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3010_fix]
	pfadd mm7, [esp + i3010_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i3010_nscoul]
	jz  .i3010_last_mno
	jmp .i3010_mno_coul
.i3010_last_mno:	
	mov   edx, [ebp + i3010_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3010_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3010_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3010_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + i3010_nri]
	dec ecx
	jecxz .i3010_end
	;# not last, iterate once more! 
	mov [ebp + i3010_nri], ecx
	jmp .i3010_outer
.i3010_end:
	femms
	add esp, 132
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret

	


.globl inl3020_3dnow
.globl _inl3020_3dnow
inl3020_3dnow:	
_inl3020_3dnow:	
.equiv		i3020_nri,		8
.equiv		i3020_iinr,		12
.equiv		i3020_jindex,		16
.equiv		i3020_jjnr,		20
.equiv		i3020_shift,		24
.equiv		i3020_shiftvec,		28
.equiv		i3020_fshift,		32
.equiv		i3020_gid,		36
.equiv		i3020_pos,		40		
.equiv		i3020_faction,		44
.equiv		i3020_charge,		48
.equiv		i3020_facel,		52
.equiv		i3020_Vc,		56			
.equiv		i3020_tabscale,		60
.equiv		i3020_VFtab,		64
			;# stack offsets for local variables 
.equiv		i3020_is3,		0
.equiv		i3020_ii3,		4
.equiv		i3020_ixO,		8
.equiv		i3020_iyO,		12
.equiv		i3020_izO,		16	
.equiv		i3020_ixH,		20  
.equiv		i3020_iyH,		28  
.equiv		i3020_izH,		36  
.equiv		i3020_iqO,		44  
.equiv		i3020_iqH,		52  
.equiv		i3020_qqO,		60  
.equiv		i3020_qqH,		68  
.equiv		i3020_vctot,		76  
.equiv		i3020_two,		84  
.equiv		i3020_n1,		92  
.equiv		i3020_tsc,		100 
.equiv		i3020_innerjjnr,	108
.equiv		i3020_innerk,		112	
.equiv		i3020_fixO,		116
.equiv		i3020_fiyO,		120
.equiv		i3020_fizO,		124
.equiv		i3020_fixH,		128 
.equiv		i3020_fiyH,		136 
.equiv		i3020_fizH,		144 
.equiv		i3020_dxO,		152
.equiv		i3020_dyO,		156
.equiv		i3020_dzO,		160
.equiv		i3020_dxH,		164 
.equiv		i3020_dyH,		172 
.equiv		i3020_dzH,		180 
.equiv		i3020_tmprsqH,		188 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 196		;# local stack space 
	femms

	mov   ecx, [ebp + i3020_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3020_charge]
	movd  mm1, [ebp + i3020_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + i3020_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3020_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	movq  mm3, [mm_two]
	movd  mm4, [ebp + i3020_tabscale]
	punpckldq mm4,mm4	    ;# spread to both halves 
	movq  [esp + i3020_two],    mm3
	movq  [esp + i3020_tsc], mm4	      
	;# assume we have at least one i particle - start directly 	 
.i3020_outer:
	mov   eax, [ebp + i3020_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3020_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3020_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3020_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i3020_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3020_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3020_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3020_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + i3020_ixO], mm5	
	movq  [esp + i3020_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i3020_ixH], mm0	
	movq [esp + i3020_iyH], mm1	
	movq [esp + i3020_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i3020_vctot], mm7
	movq  [esp + i3020_fixO],   mm7
	movd  [esp + i3020_fizO],   mm7
	movq  [esp + i3020_fixH],   mm7
	movq  [esp + i3020_fiyH],   mm7
	movq  [esp + i3020_fizH],   mm7

	mov   eax, [ebp + i3020_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3020_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i3020_innerk], edx        

	mov   esi, [ebp + i3020_pos]
	mov   edi, [ebp + i3020_faction]	
	mov   eax, [ebp + i3020_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3020_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i3020_inner_loop:	
	;# a single j particle iteration 
	mov   eax, [esp + i3020_innerjjnr]
	mov   eax, [eax]	 ;# eax=jnr offset 
    add dword ptr [esp + i3020_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + i3020_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + i3020_iqO]
	pfmul mm7, [esp + i3020_iqH]	 ;# mm6=qqO, mm7=qqH 
	movd [esp + i3020_qqO], mm6
	movq [esp + i3020_qqH], mm7
		
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3020_ixO]
	pfsubr mm1, [esp + i3020_izO]
		
	movq  [esp + i3020_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3020_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3020_ixH]
	pfsubr mm3, [esp + i3020_iyH]
	pfsubr mm4, [esp + i3020_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3020_dxH], mm2
	movq [esp + i3020_dyH], mm3
	movq [esp + i3020_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3020_tmprsqH], mm3
	
    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + i3020_tsc]
	pf2iw mm4, mm0
	movd [esp + i3020_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3020_VFtab]
	mov ecx, [esp + i3020_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3020_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3020_qqO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3020_qqO]	;# fijC=qq*FF 
	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + i3020_vctot]
	movq [esp + i3020_vctot], mm5
	movq mm3, mm7	

	;# change sign of fscal and multiply with rinv  
	pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + i3020_tsc]
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 	
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3. 
	;# now do the two hydrogens. 
	 
	movq mm0, [esp + i3020_tmprsqH] ;# mm0=rsqH 

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3020_tsc]
	pf2iw mm4, mm0
	movq [esp + i3020_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3020_VFtab]
	mov ecx, [esp + i3020_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3020_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3020_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3020_qqH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3020_qqH]	;# fijC=qq*FF 

	;# update vctot 
	pfadd mm5, [esp + i3020_vctot]
	movq [esp + i3020_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3020_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + i3020_dxO]
	movd mm1,  [esp + i3020_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3020_dxH]
	movq mm6, [esp + i3020_dyH]
	movq mm7, [esp + i3020_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3020_fixO]
	movd mm3,  [esp + i3020_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3020_fixO], mm2
	movd [esp + i3020_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3020_fixH]
	movq mm3, [esp + i3020_fiyH]
	movq mm4, [esp + i3020_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3020_fixH], mm2
	movq [esp + i3020_fiyH], mm3
	movq [esp + i3020_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 + 8], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i3020_innerk]
	jz  .i3020_updateouterdata
	jmp .i3020_inner_loop
.i3020_updateouterdata:	
	mov   ecx, [esp + i3020_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3020_fixO]
	pfadd mm7, [esp + i3020_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i3020_fixH]
	movq  mm3, [esp + i3020_fiyH]
	movq  mm1, [esp + i3020_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3, mm3	
	;# mm1 is fzH1 
	;# mm3 is fzH2 
	
	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i3020_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3020_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3020_fixO]
	pfadd mm7, [esp + i3020_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i3020_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3020_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3020_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3020_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	;# finish if last 
	dec dword ptr [ebp + i3020_nri]
	jz  .i3020_end
	;# not last, iterate once more! 
	jmp .i3020_outer
.i3020_end:
	femms
	add esp, 196
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	

.globl inl3030_3dnow
.globl _inl3030_3dnow
inl3030_3dnow:	
_inl3030_3dnow:	
.equiv		i3030_nri,		8
.equiv		i3030_iinr,		12
.equiv		i3030_jindex,		16
.equiv		i3030_jjnr,		20
.equiv		i3030_shift,		24
.equiv		i3030_shiftvec,		28
.equiv		i3030_fshift,		32
.equiv		i3030_gid,		36
.equiv		i3030_pos,		40		
.equiv		i3030_faction,		44
.equiv		i3030_charge,		48
.equiv		i3030_facel,		52
.equiv		i3030_Vc,		56			
.equiv		i3030_tabscale,		60
.equiv		i3030_VFtab,		64
			;# stack offsets for local variables 
.equiv		i3030_is3,		0
.equiv		i3030_ii3,		4
.equiv		i3030_ixO,		8
.equiv		i3030_iyO,		12
.equiv		i3030_izO,		16	
.equiv		i3030_ixH,		20  
.equiv		i3030_iyH,		28  
.equiv		i3030_izH,		36  
.equiv		i3030_qqOO,		44  
.equiv		i3030_qqOH,		52  
.equiv		i3030_qqHH,		60  
.equiv		i3030_two,		68  
.equiv		i3030_n1,		76  
.equiv		i3030_tsc,		84  
.equiv		i3030_vctot,		92  
.equiv		i3030_innerjjnr,	100
.equiv		i3030_innerk,		104	
.equiv		i3030_fixO,		108
.equiv		i3030_fiyO,		112
.equiv		i3030_fizO,		116
.equiv		i3030_fixH,		120 
.equiv		i3030_fiyH,		128 
.equiv		i3030_fizH,		136 
.equiv		i3030_dxO,		144
.equiv		i3030_dyO,		148
.equiv		i3030_dzO,		152
.equiv		i3030_dxH,		156 
.equiv		i3030_dyH,		164 
.equiv		i3030_dzH,		172 
.equiv		i3030_tmprsqH,		180 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 188		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + i3030_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3030_charge]
	movd  mm1, [ebp + i3030_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + i3030_qqOO], mm4
	movq  [esp + i3030_qqOH], mm5
	movq  [esp + i3030_qqHH], mm6
	movq  mm2, [mm_two]
	movq  [esp + i3030_two], mm2
	movd  mm3, [ebp + i3030_tabscale]
	punpckldq mm3,mm3
	movq  [esp + i3030_tsc], mm3
.i3030_outer:
	mov   eax, [ebp + i3030_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3030_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3030_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3030_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i3030_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3030_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3030_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3030_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + i3030_ixO], mm5	
	movq  [esp + i3030_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i3030_ixH], mm0	
	movq [esp + i3030_iyH], mm1	
	movq [esp + i3030_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i3030_vctot], mm7
	movq  [esp + i3030_fixO],  mm7
	movq  [esp + i3030_fizO],  mm7
	movq  [esp + i3030_fixH],  mm7
	movq  [esp + i3030_fiyH],  mm7
	movq  [esp + i3030_fizH],  mm7

	mov   eax, [ebp + i3030_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3030_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i3030_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + i3030_pos]
	mov   edi, [ebp + i3030_faction]	
	mov   eax, [ebp + i3030_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3030_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i3030_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3030_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i3030_innerjjnr],  4 ;# advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3030_ixO]
	pfsubr mm1, [esp + i3030_izO]
		
	movq  [esp + i3030_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3030_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3030_ixH]
	pfsubr mm3, [esp + i3030_iyH]
	pfsubr mm4, [esp + i3030_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3030_dxH], mm2
	movq [esp + i3030_dyH], mm3
	movq [esp + i3030_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3030_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	pfmul mm0, mm1		;# mm0=rsq  

	pfmul mm0, [esp + i3030_tsc]
	pf2iw mm4, mm0
	movd [esp + i3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3030_VFtab]
	mov ecx, [esp + i3030_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3030_qqOO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3030_qqOO]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + i3030_vctot]
	movq [esp + i3030_vctot], mm5
	movq mm3, mm7

	;# change sign of fscal and multiply with rinv  
	pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + i3030_tsc]
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3. 
	;# time for hydrogens! 

	movq mm0, [esp + i3030_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3030_tsc]
	pf2iw mm4, mm0
	movq [esp + i3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3030_VFtab]
	mov ecx, [esp + i3030_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3030_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3030_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3030_qqOH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3030_vctot]
	movq [esp + i3030_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3030_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3030_dxO]
	movd mm1,  [esp + i3030_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3030_dxH]
	movq mm6, [esp + i3030_dyH]
	movq mm7, [esp + i3030_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3030_fixO]
	movd mm3,  [esp + i3030_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3030_fixO], mm2
	movd [esp + i3030_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3030_fixH]
	movq mm3, [esp + i3030_fiyH]
	movq mm4, [esp + i3030_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3030_fixH], mm2
	movq [esp + i3030_fiyH], mm3
	movq [esp + i3030_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	;# interactions with j H1 

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3030_ixO]
	pfsubr mm1, [esp + i3030_izO]
		
	movq  [esp + i3030_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3030_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3030_ixH]
	pfsubr mm3, [esp + i3030_iyH]
	pfsubr mm4, [esp + i3030_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3030_dxH], mm2
	movq [esp + i3030_dyH], mm3
	movq [esp + i3030_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3030_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1		;# mm0=rsq  
	
	pfmul mm0, [esp + i3030_tsc]
	pf2iw mm4, mm0
	movd [esp + i3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3030_VFtab]
	mov ecx, [esp + i3030_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3030_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3030_qqOH]	;# fijC=qq*FF 

	;# update vctot  directly, force is moved to mm3 
	pfadd mm5, [esp + i3030_vctot]
	movq [esp + i3030_vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + i3030_tsc]
	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 

	movq mm0, [esp + i3030_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3030_tsc]
	pf2iw mm4, mm0
	movq [esp + i3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3030_VFtab]
	mov ecx, [esp + i3030_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3030_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3030_qqHH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3030_qqHH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3030_vctot]
	movq [esp + i3030_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3030_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 		

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3030_dxO]
	movd mm1,  [esp + i3030_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3030_dxH]
	movq mm6, [esp + i3030_dyH]
	movq mm7, [esp + i3030_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3030_fixO]
	movd mm3,  [esp + i3030_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3030_fixO], mm2
	movd [esp + i3030_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3030_fixH]
	movq mm3, [esp + i3030_fiyH]
	movq mm4, [esp + i3030_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3030_fixH], mm2
	movq [esp + i3030_fiyH], mm3
	movq [esp + i3030_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + i3030_ixO]
	pfsubr mm1, [esp + i3030_izO]
		
	movq  [esp + i3030_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3030_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3030_ixH]
	pfsubr mm3, [esp + i3030_iyH]
	pfsubr mm4, [esp + i3030_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3030_dxH], mm2
	movq [esp + i3030_dyH], mm3
	movq [esp + i3030_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3030_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1

	pfmul mm0, [esp + i3030_tsc]
	pf2iw mm4, mm0
	movd [esp + i3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3030_VFtab]
	mov ecx, [esp + i3030_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3030_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3030_qqOH]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + i3030_vctot]
	movq [esp + i3030_vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + i3030_tsc]
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 

	movq mm0, [esp + i3030_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3030_tsc]
	pf2iw mm4, mm0
	movq [esp + i3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3030_VFtab]
	mov ecx, [esp + i3030_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3030_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3030_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3030_qqHH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3030_qqHH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3030_vctot]
	movq [esp + i3030_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3030_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3030_dxO]
	movd mm1,  [esp + i3030_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3030_dxH]
	movq mm6, [esp + i3030_dyH]
	movq mm7, [esp + i3030_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3030_fixO]
	movd mm3,  [esp + i3030_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3030_fixO], mm2
	movd [esp + i3030_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3030_fixH]
	movq mm3, [esp + i3030_fiyH]
	movq mm4, [esp + i3030_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3030_fixH], mm2
	movq [esp + i3030_fiyH], mm3
	movq [esp + i3030_fizH], mm4	

	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i3030_innerk]
	jz  .i3030_updateouterdata
	jmp .i3030_inner_loop	
.i3030_updateouterdata:	
	mov   ecx, [esp + i3030_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3030_fixO]
	pfadd mm7, [esp + i3030_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i3030_fixH]
	movq  mm3, [esp + i3030_fiyH]
	movq  mm1, [esp + i3030_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 

	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i3030_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3030_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3030_fixO]
	pfadd mm7, [esp + i3030_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i3030_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3030_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3030_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + i3030_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	;# finish if last 
	dec dword ptr [ebp + i3030_nri]
	jz  .i3030_end
	;# not last, iterate once more! 
	jmp .i3030_outer
.i3030_end:
	femms
	add esp, 188
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret




.globl inl3100_3dnow
.globl _inl3100_3dnow
inl3100_3dnow:	
_inl3100_3dnow:	
.equiv		i3100_nri,		8
.equiv		i3100_iinr,		12
.equiv		i3100_jindex,		16
.equiv		i3100_jjnr,		20
.equiv		i3100_shift,		24
.equiv		i3100_shiftvec,		28
.equiv		i3100_fshift,		32
.equiv		i3100_gid,		36
.equiv		i3100_pos,		40		
.equiv		i3100_faction,		44
.equiv		i3100_charge,		48
.equiv		i3100_facel,		52
.equiv		i3100_Vc,		56			
.equiv		i3100_type,		60
.equiv		i3100_ntype,		64
.equiv		i3100_nbfp,		68	
.equiv		i3100_Vnb,		72
.equiv		i3100_tabscale,		76
.equiv		i3100_VFtab,		80
	;# stack offsets for local variables 
.equiv		i3100_is3,		0 
.equiv		i3100_ii3,		4
.equiv		i3100_ix,		8
.equiv		i3100_iy,		12
.equiv		i3100_iz,		16
.equiv		i3100_iq,		20 
.equiv		i3100_vctot,		28 
.equiv		i3100_vnbtot,		36 
.equiv		i3100_c6,		44 
.equiv		i3100_c12,		52 
.equiv		i3100_six,		60 
.equiv		i3100_twelve,		68 
.equiv		i3100_two,		76 
.equiv		i3100_n1,		84 
.equiv		i3100_tsc,		92 
.equiv		i3100_ntia,		100
.equiv		i3100_innerjjnr,	104
.equiv		i3100_innerk,		108	
.equiv		i3100_fix,		112
.equiv		i3100_fiy,		116
.equiv		i3100_fiz,		120
.equiv		i3100_dx1,		124
.equiv		i3100_dy1,		128
.equiv		i3100_dz1,		132
.equiv		i3100_dx2,		136
.equiv		i3100_dy2,		140
.equiv		i3100_dz2,		144
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 148		;# local stack space 
	femms
	;# move data to local stack  
	movq  mm0, [mm_two]
	movq  mm1, [mm_six]
	movq  mm2, [mm_twelve]
	movd  mm3, [ebp + i3100_tabscale]
	movq  [esp + i3100_two],    mm0
	movq  [esp + i3100_six],    mm1
	movq  [esp + i3100_twelve],    mm2
	punpckldq mm3,mm3
	movq  [esp + i3100_tsc], mm3	
	;# assume we have at least one i particle - start directly 	
.i3100_outer:
	mov   eax, [ebp + i3100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3100_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3100_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3100_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + i3100_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3100_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3100_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i3100_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3100_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + i3100_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i3100_ntype]
	shl   edx, 1
	mov   [esp + i3100_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3100_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3100_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + i3100_ix], mm0	
	movd  [esp + i3100_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + i3100_vctot],  mm7
	movq  [esp + i3100_vnbtot], mm7
	movq  [esp + i3100_fix],    mm7
	movd  [esp + i3100_fiz],    mm7

	mov   eax, [ebp + i3100_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3100_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + i3100_pos]
	mov   edi, [ebp + i3100_faction]	
	mov   eax, [ebp + i3100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3100_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + i3100_innerk], edx    ;# number of innerloop atoms 
	jge   .i3100_unroll_loop
	jmp   .i3100_finish_inner
.i3100_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + i3100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3100_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3100_charge]    ;# base of charge[] 
	movq mm5, [esp + i3100_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + i3100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i3100_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i3100_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i3100_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i3100_c6], mm5
	movq [esp + i3100_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3100_pos]

	movq  mm0, [esp + i3100_ix]
	movd  mm1, [esp + i3100_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3100_dx1], mm4	     ;# store dr 
	movd  [esp + i3100_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3100_dx2], mm6	     ;# store dr 
	movd  [esp + i3100_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3100_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i3100_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3100_VFtab]
	mov ecx, [esp + i3100_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3100_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3100_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + i3100_tsc]
	
	pfmul mm1, [esp + i3100_c12]

	pfmul mm2, [esp + i3100_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = vnb12-vnb6 

	pfmul mm2, [esp + i3100_six]
	pfmul mm1, [esp + i3100_twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	;# mm1=	(12*vnb12-6*vnb6)*rinv11 

	pfsub mm1, mm3

	;# update vctot 
	pfadd mm5, [esp + i3100_vctot]      ;# add the earlier value 
	movq [esp + i3100_vctot], mm5       ;# store the sum       

 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3100_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3100_dx1]	;# fetch dr 
	movd mm3,  [esp + i3100_dz1]

	;# update vnbtot 
	pfadd mm4, [esp + i3100_vnbtot]      ;# add the earlier value 
	movq [esp + i3100_vnbtot], mm4       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3100_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3100_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3100_fix]
	movd mm1,  [esp + i3100_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3100_fix], mm0
	movd [esp + i3100_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3100_innerk],  2
	jl    .i3100_finish_inner
	jmp   .i3100_unroll_loop
.i3100_finish_inner:	
	and dword ptr [esp + i3100_innerk],  1
	jnz  .i3100_single_inner
	jmp  .i3100_updateouterdata		
.i3100_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3100_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i3100_charge]
	movd mm5, [esp + i3100_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + i3100_nbfp]
	mov ecx, [ebp + i3100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i3100_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i3100_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i3100_c12], mm5


	mov   esi, [ebp + i3100_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3100_ix]
	movd  mm1, [esp + i3100_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3100_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3100_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3100_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3100_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3100_VFtab]
	mov ecx, [esp + i3100_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3100_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  
	
	;# at this point mm5 contains vcoul and mm3 fijC 

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + i3100_tsc]
	
	pfmul mm1, [esp + i3100_c12]

	pfmul mm2, [esp + i3100_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = vnb12-vnb6 

	pfmul mm2, [esp + i3100_six]
	pfmul mm1, [esp + i3100_twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	;# mm1=	(12*vnb12-6*vnb6)*rinv11 

	pfsub mm1, mm3

	;# update vctot 
	pfadd mm5, [esp + i3100_vctot]      ;# add the earlier value 
	movq [esp + i3100_vctot], mm5       ;# store the sum       

 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3100_dx1]
	movd mm3,  [esp + i3100_dz1]

	;# update vnbtot 
	pfadd mm4, [esp + i3100_vnbtot]      ;# add the earlier value 
	movq [esp + i3100_vnbtot], mm4       ;# store the sum       

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3100_fix]
	movd mm1,  [esp + i3100_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3100_fix], mm0
	movd [esp + i3100_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3100_updateouterdata:	
	mov   ecx, [esp + i3100_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3100_fix]
	pfadd mm7, [esp + i3100_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3100_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3100_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3100_fix] 
	pfadd mm7, [esp + i3100_fiz] 
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + i3100_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3100_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3100_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3100_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
 
	movq  mm7, [esp + i3100_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3100_Vnb] 
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + i3100_nri]
	dec ecx
	jecxz .i3100_end
	;# not last, iterate once more! 
	mov [ebp + i3100_nri], ecx
	jmp .i3100_outer
.i3100_end:
	femms
	add esp, 148
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret







.globl inl3110_3dnow
.globl _inl3110_3dnow
inl3110_3dnow:	
_inl3110_3dnow:	
.equiv		i3110_nri,		8
.equiv		i3110_iinr,		12
.equiv		i3110_jindex,		16
.equiv		i3110_jjnr,		20
.equiv		i3110_shift,		24
.equiv		i3110_shiftvec,		28
.equiv		i3110_fshift,		32
.equiv		i3110_gid,		36
.equiv		i3110_pos,		40		
.equiv		i3110_faction,		44
.equiv		i3110_charge,		48
.equiv		i3110_facel,		52
.equiv		i3110_Vc,		56			
.equiv		i3110_type,		60
.equiv		i3110_ntype,		64
.equiv		i3110_nbfp,		68	
.equiv		i3110_Vnb,		72
.equiv		i3110_tabscale,		76
.equiv		i3110_VFtab,		80
.equiv		i3110_nsatoms,		84	
	;# stack offsets for local variables 
.equiv		i3110_is3,		0
.equiv		i3110_ii3,		4
.equiv		i3110_shX,		8
.equiv		i3110_shY,		12 
.equiv		i3110_shZ,		16	
.equiv		i3110_ix,		20
.equiv		i3110_iy,		24
.equiv		i3110_iz,		28	
.equiv		i3110_iq,		32  
.equiv		i3110_vctot,		40  
.equiv		i3110_vnbtot,		48  
.equiv		i3110_c6,		56  
.equiv		i3110_c12,		64  
.equiv		i3110_six,		72  
.equiv		i3110_twelve,		80  
.equiv		i3110_two,		88  
.equiv		i3110_n1,		96  
.equiv		i3110_tsc,		104 
.equiv		i3110_ntia,		112
.equiv		i3110_innerjjnr0,	116
.equiv		i3110_innerk0,		120	
.equiv		i3110_innerjjnr,	124
.equiv		i3110_innerk,		128	
.equiv		i3110_fix,		132
.equiv		i3110_fiy,		136
.equiv		i3110_fiz,		140
.equiv		i3110_dx1,		144
.equiv		i3110_dy1,		148
.equiv		i3110_dz1,		152
.equiv		i3110_dx2,		156
.equiv		i3110_dy2,		160
.equiv		i3110_dz2,		164
.equiv		i3110_nsvdwc,		168
.equiv		i3110_nscoul,		172
.equiv		i3110_nsvdw,		176
.equiv		i3110_solnr,		180		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 184		;# local stack space 
	femms
	movq  mm0, [mm_six]
	movq  mm1, [mm_twelve]
	movq  [esp + i3110_six],    mm0
	movq  [esp + i3110_twelve], mm1
	movq  mm2, [mm_two]
	movd  mm3, [ebp + i3110_tabscale]
	movq  [esp + i3110_two],    mm2
	punpckldq mm3,mm3
	movq  [esp + i3110_tsc], mm3	
	;# assume we have at least one i particle - start directly 		
.i3110_outer:
	mov   eax, [ebp + i3110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3110_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3110_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3110_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + i3110_shX], mm0
	movd  [esp + i3110_shZ], mm1

	mov   ecx, [ebp + i3110_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3110_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + i3110_nsatoms]
	add dword ptr [ebp + i3110_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + i3110_nsvdwc], edx
	mov   [esp + i3110_nscoul], eax
	mov   [esp + i3110_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + i3110_vctot],  mm7
	movq  [esp + i3110_vnbtot], mm7
	mov   [esp + i3110_solnr],  ebx
	
	mov   eax, [ebp + i3110_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3110_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + i3110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3110_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + i3110_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + i3110_pos]
	mov   edi, [ebp + i3110_faction]
	
	mov   ecx, [esp + i3110_nsvdwc]
	cmp   ecx,  0
	jnz   .i3110_mno_vdwc
	jmp   .i3110_testcoul
.i3110_mno_vdwc:
	mov   ebx,  [esp + i3110_solnr]
	inc   dword ptr [esp + i3110_solnr]
	mov   edx, [ebp + i3110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i3110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3110_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + i3110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i3110_ntype]
	shl   edx, 1
	mov   [esp + i3110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3110_pos]    ;# eax = base of pos[] 
	mov   [esp + i3110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i3110_shX]
	pfadd mm1, [esp + i3110_shZ]
	movq  [esp + i3110_ix], mm0	
	movd  [esp + i3110_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i3110_fix],   mm7
	movd  [esp + i3110_fiz],   mm7

	mov   ecx, [esp + i3110_innerjjnr0]
	mov   [esp + i3110_innerjjnr], ecx
	mov   edx, [esp + i3110_innerk0]
    sub   edx,  2
    mov   [esp + i3110_innerk], edx    ;# number of innerloop atoms 
	jge   .i3110_unroll_vdwc_loop
	jmp   .i3110_finish_vdwc_inner
.i3110_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i3110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3110_charge]    ;# base of charge[] 
	movq mm5, [esp + i3110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + i3110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i3110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i3110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i3110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i3110_c6], mm5
	movq [esp + i3110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3110_pos]

	movq  mm0, [esp + i3110_ix]
	movd  mm1, [esp + i3110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3110_dx1], mm4	     ;# store dr 
	movd  [esp + i3110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3110_dx2], mm6	     ;# store dr 
	movd  [esp + i3110_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3110_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i3110_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3110_VFtab]
	mov ecx, [esp + i3110_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3110_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3110_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + i3110_tsc]
	
	pfmul mm1, [esp + i3110_c12]

	pfmul mm2, [esp + i3110_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = vnb12-vnb6 

	pfmul mm2, [esp + i3110_six]
	pfmul mm1, [esp + i3110_twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	;# mm1=	(12*vnb12-6*vnb6)*rinv11 

	pfsub mm1, mm3

	;# update vctot 
	pfadd mm5, [esp + i3110_vctot]      ;# add the earlier value 
	movq [esp + i3110_vctot], mm5       ;# store the sum       

 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3110_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3110_dx1]	;# fetch dr 
	movd mm3,  [esp + i3110_dz1]

	;# update vnbtot 
	pfadd mm4, [esp + i3110_vnbtot]      ;# add the earlier value 
	movq [esp + i3110_vnbtot], mm4       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3110_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3110_fix]
	movd mm1,  [esp + i3110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3110_fix], mm0
	movd [esp + i3110_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7	
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3110_innerk],  2
	jl    .i3110_finish_vdwc_inner
	jmp   .i3110_unroll_vdwc_loop
.i3110_finish_vdwc_inner:	
	and dword ptr [esp + i3110_innerk],  1
	jnz  .i3110_single_vdwc_inner
	jmp  .i3110_updateouterdata_vdwc		
.i3110_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i3110_charge]
	movd mm5, [esp + i3110_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + i3110_nbfp]
	mov ecx, [ebp + i3110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i3110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i3110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i3110_c12], mm5


	mov   esi, [ebp + i3110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3110_ix]
	movd  mm1, [esp + i3110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3110_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3110_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3110_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3110_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3110_VFtab]
	mov ecx, [esp + i3110_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3110_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + i3110_tsc]
	
	pfmul mm1, [esp + i3110_c12]

	pfmul mm2, [esp + i3110_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = vnb12-vnb6 

	pfmul mm2, [esp + i3110_six]
	pfmul mm1, [esp + i3110_twelve]

	pfsub mm1, mm2
	pfmul mm1, mm0	;# mm1=	(12*vnb12-6*vnb6)*rinv11 

	pfsub mm1, mm3

	;# update vctot 
	pfadd mm5, [esp + i3110_vctot]      ;# add the earlier value 
	movq [esp + i3110_vctot], mm5       ;# store the sum       

 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3110_dx1]
	movd mm3,  [esp + i3110_dz1]

	;# update vnbtot 
	pfadd mm4, [esp + i3110_vnbtot]      ;# add the earlier value 
	movq [esp + i3110_vnbtot], mm4       ;# store the sum       

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3110_fix]
	movd mm1,  [esp + i3110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3110_fix], mm0
	movd [esp + i3110_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3110_updateouterdata_vdwc:	
	mov   ecx, [esp + i3110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3110_fix]
	pfadd mm7, [esp + i3110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3110_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3110_fix]
	pfadd mm7, [esp + i3110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i3110_nsvdwc]
	jz  .i3110_testcoul
	jmp .i3110_mno_vdwc
.i3110_testcoul:	
	mov  ecx, [esp + i3110_nscoul]
	cmp  ecx,  0
	jnz  .i3110_mno_coul
	jmp  .i3110_testvdw
.i3110_mno_coul:
	mov   ebx,  [esp + i3110_solnr]
	inc   dword ptr [esp + i3110_solnr]
	mov   edx, [ebp + i3110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i3110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3110_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3110_pos]    ;# eax = base of pos[] 
	mov   [esp + i3110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i3110_shX]
	pfadd mm1, [esp + i3110_shZ]
	movq  [esp + i3110_ix], mm0	
	movd  [esp + i3110_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i3110_fix],   mm7
	movd  [esp + i3110_fiz],   mm7

	mov   ecx, [esp + i3110_innerjjnr0]
	mov   [esp + i3110_innerjjnr], ecx
	mov   edx, [esp + i3110_innerk0]
    sub   edx,  2
    mov   [esp + i3110_innerk], edx    ;# number of innerloop atoms 
	jge   .i3110_unroll_coul_loop
	jmp   .i3110_finish_coul_inner
.i3110_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i3110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3110_charge]    ;# base of charge[] 
	movq mm5, [esp + i3110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3110_pos]

	movq  mm0, [esp + i3110_ix]
	movd  mm1, [esp + i3110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3110_dx1], mm4	     ;# store dr 
	movd  [esp + i3110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3110_dx2], mm6	     ;# store dr 
	movd  [esp + i3110_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3110_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i3110_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3110_VFtab]
	mov ecx, [esp + i3110_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3110_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3110_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3110_vctot]      ;# add the earlier value 
	movq [esp + i3110_vctot], mm5       ;# store the sum       

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + i3110_tsc]	
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3110_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3110_dx1]	;# fetch dr 
	movd mm3,  [esp + i3110_dz1]

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3110_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3110_fix]
	movd mm1,  [esp + i3110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3110_fix], mm0
	movd [esp + i3110_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3110_innerk],  2
	jl    .i3110_finish_coul_inner
	jmp   .i3110_unroll_coul_loop
.i3110_finish_coul_inner:	
	and dword ptr [esp + i3110_innerk],  1
	jnz  .i3110_single_coul_inner
	jmp  .i3110_updateouterdata_coul		
.i3110_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i3110_charge]
	movd mm5, [esp + i3110_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + i3110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3110_ix]
	movd  mm1, [esp + i3110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3110_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3110_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3110_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3110_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3110_VFtab]
	mov ecx, [esp + i3110_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3110_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3110_vctot]      ;# add the earlier value 
	movq [esp + i3110_vctot], mm5       ;# store the sum       
	
	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3110_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3110_dx1]
	movd mm3,  [esp + i3110_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3110_fix]
	movd mm1,  [esp + i3110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3110_fix], mm0
	movd [esp + i3110_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3110_updateouterdata_coul:	
	mov   ecx, [esp + i3110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3110_fix]
	pfadd mm7, [esp + i3110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3110_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3110_fix]
	pfadd mm7, [esp + i3110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i3110_nscoul]
	jz  .i3110_testvdw
	jmp .i3110_mno_coul
.i3110_testvdw:	
	mov  ecx, [esp + i3110_nsvdw]
	cmp  ecx,  0
	jnz  .i3110_mno_vdw
	jmp  .i3110_last_mno
.i3110_mno_vdw:
	mov   ebx,  [esp + i3110_solnr]
	inc   dword ptr [esp + i3110_solnr]

	mov   edx, [ebp + i3110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i3110_ntype]
	shl   edx, 1
	mov   [esp + i3110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3110_pos]    ;# eax = base of pos[] 
	mov   [esp + i3110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i3110_shX]
	pfadd mm1, [esp + i3110_shZ]
	movq  [esp + i3110_ix], mm0	
	movd  [esp + i3110_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i3110_fix],   mm7
	movd  [esp + i3110_fiz],   mm7

	mov   ecx, [esp + i3110_innerjjnr0]
	mov   [esp + i3110_innerjjnr], ecx
	mov   edx, [esp + i3110_innerk0]
    sub   edx,  2
    mov   [esp + i3110_innerk], edx    ;# number of innerloop atoms 
	jge   .i3110_unroll_vdw_loop
	jmp   .i3110_finish_vdw_inner
.i3110_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i3110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + i3110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i3110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i3110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i3110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i3110_c6], mm5
	movq [esp + i3110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3110_pos]

	movq  mm0, [esp + i3110_ix]
	movd  mm1, [esp + i3110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3110_dx1], mm4	     ;# store dr 
	movd  [esp + i3110_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3110_dx2], mm6	     ;# store dr 
	movd  [esp + i3110_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrcp mm0, mm4	             ;# lookup reciprocal seed  
    pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
                  	        	;# amd 3dnow N-R iteration to get full precision. 
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	 ;# do potential and fscal 
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i3110_c12]
	pfmul mm4, [esp + i3110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i3110_six]

	pfmul mm5, [esp + i3110_twelve]
	pfsub mm5,mm4
 	pfmul mm0, mm5    ;# mm0 is total fscal now 	

	prefetchw [esp + i3110_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3110_dx1]	;# fetch dr 
	movd mm3,  [esp + i3110_dz1]

	;# update vnbtot 
	pfadd mm6, [esp + i3110_vnbtot]      ;# add the earlier value 
	movq [esp + i3110_vnbtot], mm6       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3110_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3110_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3110_fix]
	movd mm1,  [esp + i3110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3110_fix], mm0
	movd [esp + i3110_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3110_innerk],  2
	jl    .i3110_finish_vdw_inner
	jmp   .i3110_unroll_vdw_loop
.i3110_finish_vdw_inner:	
	and dword ptr [esp + i3110_innerk],  1
	jnz  .i3110_single_vdw_inner
	jmp  .i3110_updateouterdata_vdw		
.i3110_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + i3110_nbfp]
	mov ecx, [ebp + i3110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i3110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i3110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i3110_c12], mm5

	mov   esi, [ebp + i3110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3110_ix]
	movd  mm1, [esp + i3110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3110_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3110_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
	pfrcp mm0,mm4
	pfrcpit1 mm4,mm0				
	pfrcpit2 mm4,mm0	;# mm4=invsq  
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + i3110_c12]
	pfmul mm4, [esp + i3110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4

	pfmul mm4, [esp + i3110_six]

	pfmul mm5, [esp + i3110_twelve]
	pfsub mm5, mm4
 	pfmul mm0, mm5    ;# mm0 is total fscal now 

	;# update vnbtot 
	pfadd mm6, [esp + i3110_vnbtot]      ;# add the earlier value 
	movq [esp + i3110_vnbtot], mm6       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3110_dx1]
	movd mm3,  [esp + i3110_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3110_fix]
	movd mm1,  [esp + i3110_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3110_fix], mm0
	movd [esp + i3110_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3110_updateouterdata_vdw:	
	mov   ecx, [esp + i3110_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3110_fix]
	pfadd mm7, [esp + i3110_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3110_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3110_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3110_fix]
	pfadd mm7, [esp + i3110_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i3110_nsvdw]
	jz  .i3110_last_mno
	jmp .i3110_mno_vdw
	
.i3110_last_mno:	
	mov   edx, [ebp + i3110_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3110_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3110_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3110_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i3110_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3110_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + i3110_nri]
	dec ecx
	jecxz .i3110_end
	;# not last, iterate once more! 
	mov [ebp + i3110_nri], ecx
	jmp .i3110_outer
.i3110_end:
	femms
	add esp, 184
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret


	

.globl inl3120_3dnow
.globl _inl3120_3dnow
inl3120_3dnow:	
_inl3120_3dnow:	
.equiv		i3120_nri,		8
.equiv		i3120_iinr,		12
.equiv		i3120_jindex,		16
.equiv		i3120_jjnr,		20
.equiv		i3120_shift,		24
.equiv		i3120_shiftvec,		28
.equiv		i3120_fshift,		32
.equiv		i3120_gid,		36
.equiv		i3120_pos,		40		
.equiv		i3120_faction,		44
.equiv		i3120_charge,		48
.equiv		i3120_facel,		52
.equiv		i3120_Vc,		56			
.equiv		i3120_type,		60
.equiv		i3120_ntype,		64
.equiv		i3120_nbfp,		68	
.equiv		i3120_Vnb,		72
.equiv		i3120_tabscale,		76
.equiv		i3120_VFtab,		80
			;# stack offsets for local variables 
.equiv		i3120_is3,		0
.equiv		i3120_ii3,		4
.equiv		i3120_ixO,		8
.equiv		i3120_iyO,		12
.equiv		i3120_izO,		16	
.equiv		i3120_ixH,		20  
.equiv		i3120_iyH,		28  
.equiv		i3120_izH,		36  
.equiv		i3120_iqO,		44  
.equiv		i3120_iqH,		52  
.equiv		i3120_qqO,		60  
.equiv		i3120_qqH,		68  
.equiv		i3120_vctot,		76  
.equiv		i3120_vnbtot,		84  
.equiv		i3120_c6,		92  
.equiv		i3120_c12,		100 
.equiv		i3120_six,		108 
.equiv		i3120_twelve,		116 
.equiv		i3120_two,		124 
.equiv		i3120_n1,		132 
.equiv		i3120_tsc,		140 
.equiv		i3120_ntia,		148 
.equiv		i3120_innerjjnr,	156
.equiv		i3120_innerk,		160	
.equiv		i3120_fixO,		164
.equiv		i3120_fiyO,		168
.equiv		i3120_fizO,		172
.equiv		i3120_fixH,		176 
.equiv		i3120_fiyH,		184 
.equiv		i3120_fizH,		192 
.equiv		i3120_dxO,		200
.equiv		i3120_dyO,		204
.equiv		i3120_dzO,		208
.equiv		i3120_dxH,		212 
.equiv		i3120_dyH,		220 
.equiv		i3120_dzH,		228 
.equiv		i3120_tmprsqH,		236 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 244		;# local stack space 
	femms

	mov   ecx, [ebp + i3120_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3120_charge]
	movd  mm1, [ebp + i3120_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1
	movq  [esp + i3120_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3120_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + i3120_type] 	
	mov   edx, [edx + ebx*4]
	shl   edx, 1
	mov   ecx, edx		        
	imul  ecx, [ebp + i3120_ntype]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + i3120_ntia], ecx
	 	
	movq  mm3, [mm_two]
	movq  mm4, [mm_six]
	movq  mm5, [mm_twelve]
	movq  mm6, [ebp + i3120_tabscale]
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + i3120_two], mm3
	movq  [esp + i3120_six], mm4
	movq  [esp + i3120_twelve], mm5
	movq  [esp + i3120_tsc], mm6	      
 	;# assume we have at least one i particle - start directly 	
.i3120_outer:
	mov   eax, [ebp + i3120_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3120_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3120_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3120_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i3120_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3120_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3120_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3120_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + i3120_ixO], mm5	
	movq  [esp + i3120_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i3120_ixH], mm0	
	movq [esp + i3120_iyH], mm1	
	movq [esp + i3120_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i3120_vctot], mm7
	movq  [esp + i3120_vnbtot], mm7
	movq  [esp + i3120_fixO],   mm7
	movd  [esp + i3120_fizO],   mm7
	movq  [esp + i3120_fixH],   mm7
	movq  [esp + i3120_fiyH],   mm7
	movq  [esp + i3120_fizH],   mm7

	mov   eax, [ebp + i3120_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3120_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i3120_innerk], edx        

	mov   esi, [ebp + i3120_pos]
	mov   edi, [ebp + i3120_faction]	
	mov   eax, [ebp + i3120_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3120_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i3120_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3120_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i3120_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + i3120_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + i3120_iqO]
	pfmul mm7, [esp + i3120_iqH]	;# mm6=qqO, mm7=qqH 
	movd [esp + i3120_qqO], mm6
	movq [esp + i3120_qqH], mm7

	mov ecx, [ebp + i3120_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + i3120_nbfp]
	shl edx, 1
	add edx, [esp + i3120_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i3120_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i3120_c12], mm5	
			
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3120_ixO]
	pfsubr mm1, [esp + i3120_izO]
		
	movq  [esp + i3120_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3120_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3120_ixH]
	pfsubr mm3, [esp + i3120_iyH]
	pfsubr mm4, [esp + i3120_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3120_dxH], mm2
	movq [esp + i3120_dyH], mm3
	movq [esp + i3120_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3120_tmprsqH], mm3
	
    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + i3120_tsc]
	pf2iw mm4, mm0
	movd [esp + i3120_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3120_VFtab]
	mov ecx, [esp + i3120_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3120_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3120_qqO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3120_qqO]	;# fijC=qq*FF 
	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + i3120_vctot]
	movq [esp + i3120_vctot], mm5

	movq mm3, mm7
	pfmul mm3, [esp + i3120_tsc]
	
	;# nontabulated LJ - mm1 is invsqrt. - keep mm1! 
	movq mm0, mm1
	pfmul mm0, mm0		;# mm0 is invsq 
	movq mm2, mm0
	pfmul mm2, mm0
	pfmul mm2, mm0		;# mm2 = rinvsix 
	movq mm4, mm2
	pfmul mm4, mm4		;# mm4=rinvtwelve 

	pfmul mm4, [esp + i3120_c12]
	pfmul mm2, [esp + i3120_c6]
	movq mm5, mm4
	pfsub mm5, mm2		;# mm5=vnb12-vnb6 

	pfmul mm2, [esp + i3120_six]
	pfmul mm4, [esp + i3120_twelve]
	pfsub mm4, mm2
	pfmul mm4, mm1    ;# mm4=(12*vnb12-6*vnb6)*rinv11 

	pfsubr mm3, mm4 
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 
	
	;# update vnbtot  
	pfadd mm5, [esp + i3120_vnbtot]      ;# add the earlier value 
	movq [esp + i3120_vnbtot], mm5       ;# store the sum       
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3. 
	;# now do the two hydrogens. 
	movq mm0, [esp + i3120_tmprsqH] ;# mm0=rsqH 

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3120_tsc]
	pf2iw mm4, mm0
	movq [esp + i3120_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3120_VFtab]
	mov ecx, [esp + i3120_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3120_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3120_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3120_qqH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3120_qqH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3120_vctot]
	movq [esp + i3120_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7
	pfmul mm4, [esp + i3120_tsc]	
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + i3120_dxO]
	movd mm1,  [esp + i3120_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3120_dxH]
	movq mm6, [esp + i3120_dyH]
	movq mm7, [esp + i3120_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3120_fixO]
	movd mm3,  [esp + i3120_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3120_fixO], mm2
	movd [esp + i3120_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3120_fixH]
	movq mm3, [esp + i3120_fiyH]
	movq mm4, [esp + i3120_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3120_fixH], mm2
	movq [esp + i3120_fiyH], mm3
	movq [esp + i3120_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i3120_innerk]
	jz  .i3120_updateouterdata
	jmp .i3120_inner_loop
.i3120_updateouterdata:	
	mov   ecx, [esp + i3120_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3120_fixO]
	pfadd mm7, [esp + i3120_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i3120_fixH]
	movq  mm3, [esp + i3120_fiyH]
	movq  mm1, [esp + i3120_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 
	
	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i3120_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3120_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3120_fixO]
	pfadd mm7, [esp + i3120_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i3120_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3120_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3120_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3120_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i3120_vnbtot]     
	pfacc mm7,mm7	          ;# same for Vnb 
	
	mov   eax, [ebp + i3120_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 
	;# finish if last 
	dec dword ptr [ebp + i3120_nri]
	jz  .i3120_end
	;# not last, iterate once more! 
	jmp .i3120_outer
.i3120_end:
	femms
	add esp, 244
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret





.globl inl3130_3dnow
.globl _inl3130_3dnow
inl3130_3dnow:	
_inl3130_3dnow:	
.equiv		i3130_nri,		8
.equiv		i3130_iinr,		12
.equiv		i3130_jindex,		16
.equiv		i3130_jjnr,		20
.equiv		i3130_shift,		24
.equiv		i3130_shiftvec,		28
.equiv		i3130_fshift,		32
.equiv		i3130_gid,		36
.equiv		i3130_pos,		40		
.equiv		i3130_faction,		44
.equiv		i3130_charge,		48
.equiv		i3130_facel,		52
.equiv		i3130_Vc,		56			
.equiv		i3130_type,		60
.equiv		i3130_ntype,		64
.equiv		i3130_nbfp,		68	
.equiv		i3130_Vnb,		72
.equiv		i3130_tabscale,		76
.equiv		i3130_VFtab,		80
			;# stack offsets for local variables 
.equiv		i3130_is3,		0
.equiv		i3130_ii3,		4
.equiv		i3130_ixO,		8
.equiv		i3130_iyO,		12
.equiv		i3130_izO,		16	
.equiv		i3130_ixH,		20  
.equiv		i3130_iyH,		28  
.equiv		i3130_izH,		36  
.equiv		i3130_qqOO,		44  
.equiv		i3130_qqOH,		52  
.equiv		i3130_qqHH,		60  
.equiv		i3130_c6,		68  
.equiv		i3130_c12,		76  
.equiv		i3130_six,		84  
.equiv		i3130_twelve,		92  
.equiv		i3130_two,		100 
.equiv		i3130_n1,		108 
.equiv		i3130_tsc,		116 
.equiv		i3130_vctot,		124 
.equiv		i3130_vnbtot,		132 
.equiv		i3130_innerjjnr,	140
.equiv		i3130_innerk,		144	
.equiv		i3130_fixO,		148
.equiv		i3130_fiyO,		152
.equiv		i3130_fizO,		156
.equiv		i3130_fixH,		160 
.equiv		i3130_fiyH,		168 
.equiv		i3130_fizH,		176 
.equiv		i3130_dxO,		184
.equiv		i3130_dyO,		188
.equiv		i3130_dzO,		192
.equiv		i3130_dxH,		200 
.equiv		i3130_dyH,		208 
.equiv		i3130_dzH,		216 
.equiv		i3130_tmprsqH,		224 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 232		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + i3130_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3130_charge]
	movd  mm1, [ebp + i3130_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + i3130_qqOO], mm4
	movq  [esp + i3130_qqOH], mm5
	movq  [esp + i3130_qqHH], mm6
	mov   edx, [ebp + i3130_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + i3130_ntype]
	add   edx, ecx
	mov   eax, [ebp + i3130_nbfp]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + i3130_c6], mm0
	movq  [esp + i3130_c12], mm1
	movq  mm2, [mm_two]
	movq  mm3, [mm_six]
	movq  mm4, [mm_twelve]
	movq  [esp + i3130_two], mm2
	movq  [esp + i3130_six], mm3
	movq  [esp + i3130_twelve], mm4
	movd  mm5, [ebp + i3130_tabscale]
	punpckldq mm5,mm5
	movq  [esp + i3130_tsc], mm5
.i3130_outer:
	mov   eax, [ebp + i3130_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3130_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3130_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3130_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i3130_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3130_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3130_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3130_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + i3130_ixO], mm5	
	movq  [esp + i3130_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i3130_ixH], mm0	
	movq [esp + i3130_iyH], mm1	
	movq [esp + i3130_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i3130_vctot], mm7
	movq  [esp + i3130_vnbtot], mm7
	movq  [esp + i3130_fixO],  mm7
	movq  [esp + i3130_fizO],  mm7
	movq  [esp + i3130_fixH],  mm7
	movq  [esp + i3130_fiyH],  mm7
	movq  [esp + i3130_fizH],  mm7

	mov   eax, [ebp + i3130_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3130_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i3130_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + i3130_pos]
	mov   edi, [ebp + i3130_faction]	
	mov   eax, [ebp + i3130_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3130_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i3130_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3130_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i3130_innerjjnr],  4 ;# advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3130_ixO]
	pfsubr mm1, [esp + i3130_izO]
		
	movq  [esp + i3130_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3130_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3130_ixH]
	pfsubr mm3, [esp + i3130_iyH]
	pfsubr mm4, [esp + i3130_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3130_dxH], mm2
	movq [esp + i3130_dyH], mm3
	movq [esp + i3130_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3130_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	pfmul mm0, mm1		;# mm0=rsq  

	pfmul mm0, [esp + i3130_tsc]
	pf2iw mm4, mm0
	movd [esp + i3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3130_VFtab]
	mov ecx, [esp + i3130_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3130_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3130_qqOO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3130_qqOO]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + i3130_vctot]
	movq [esp + i3130_vctot], mm5
	movq mm3, mm7
	pfmul mm3, [esp + i3130_tsc]
	
	movq mm5, mm1
	pfmul mm5,mm5
	movq mm4, mm5
	pfmul mm4,mm5
	pfmul mm4,mm5
	movq mm5, mm4
	pfmul mm5,mm5	;# mm4=rinvsix, mm5=rinvtwelve 

	pfmul mm4, [esp + i3130_c6]
	pfmul mm5, [esp + i3130_c12]
	movq mm6,mm5
	pfsub mm6,mm4

	pfmul mm4, [esp + i3130_six]
	pfmul mm5, [esp + i3130_twelve]
	pfsub mm5,mm4
	pfmul mm5, mm1
	pfsubr mm3, mm5

 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 
	
	;# update vnbtot  
	pfadd mm6, [esp + i3130_vnbtot]      ;# add the earlier value 
	movq [esp + i3130_vnbtot], mm6       ;# store the sum       
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3. 
	;# time for hydrogens! 

	movq mm0, [esp + i3130_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3130_tsc]
	pf2iw mm4, mm0
	movq [esp + i3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3130_VFtab]
	mov ecx, [esp + i3130_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3130_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3130_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3130_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3130_qqOH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3130_vctot]
	movq [esp + i3130_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3130_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3130_dxO]
	movd mm1,  [esp + i3130_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3130_dxH]
	movq mm6, [esp + i3130_dyH]
	movq mm7, [esp + i3130_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3130_fixO]
	movd mm3,  [esp + i3130_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3130_fixO], mm2
	movd [esp + i3130_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3130_fixH]
	movq mm3, [esp + i3130_fiyH]
	movq mm4, [esp + i3130_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3130_fixH], mm2
	movq [esp + i3130_fiyH], mm3
	movq [esp + i3130_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	;# interactions with j H1 

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3130_ixO]
	pfsubr mm1, [esp + i3130_izO]
		
	movq  [esp + i3130_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3130_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3130_ixH]
	pfsubr mm3, [esp + i3130_iyH]
	pfsubr mm4, [esp + i3130_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3130_dxH], mm2
	movq [esp + i3130_dyH], mm3
	movq [esp + i3130_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3130_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1		;# mm0=rsq  
	
	pfmul mm0, [esp + i3130_tsc]
	pf2iw mm4, mm0
	movd [esp + i3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3130_VFtab]
	mov ecx, [esp + i3130_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3130_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3130_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3130_qqOH]	;# fijC=qq*FF 

	;# update vctot  directly, force is moved to mm3 
	pfadd mm5, [esp + i3130_vctot]
	movq [esp + i3130_vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + i3130_tsc]
	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 

	movq mm0, [esp + i3130_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3130_tsc]
	pf2iw mm4, mm0
	movq [esp + i3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3130_VFtab]
	mov ecx, [esp + i3130_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3130_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3130_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3130_qqHH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3130_qqHH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3130_vctot]
	movq [esp + i3130_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3130_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 		

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3130_dxO]
	movd mm1,  [esp + i3130_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3130_dxH]
	movq mm6, [esp + i3130_dyH]
	movq mm7, [esp + i3130_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3130_fixO]
	movd mm3,  [esp + i3130_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3130_fixO], mm2
	movd [esp + i3130_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3130_fixH]
	movq mm3, [esp + i3130_fiyH]
	movq mm4, [esp + i3130_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3130_fixH], mm2
	movq [esp + i3130_fiyH], mm3
	movq [esp + i3130_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + i3130_ixO]
	pfsubr mm1, [esp + i3130_izO]
		
	movq  [esp + i3130_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3130_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3130_ixH]
	pfsubr mm3, [esp + i3130_iyH]
	pfsubr mm4, [esp + i3130_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3130_dxH], mm2
	movq [esp + i3130_dyH], mm3
	movq [esp + i3130_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3130_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1

	pfmul mm0, [esp + i3130_tsc]
	pf2iw mm4, mm0
	movd [esp + i3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3130_VFtab]
	mov ecx, [esp + i3130_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3130_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3130_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3130_qqOH]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + i3130_vctot]
	movq [esp + i3130_vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + i3130_tsc]
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 

	movq mm0, [esp + i3130_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3130_tsc]
	pf2iw mm4, mm0
	movq [esp + i3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3130_VFtab]
	mov ecx, [esp + i3130_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3130_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3130_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3130_qqHH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3130_qqHH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3130_vctot]
	movq [esp + i3130_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3130_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3130_dxO]
	movd mm1,  [esp + i3130_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3130_dxH]
	movq mm6, [esp + i3130_dyH]
	movq mm7, [esp + i3130_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3130_fixO]
	movd mm3,  [esp + i3130_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3130_fixO], mm2
	movd [esp + i3130_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3130_fixH]
	movq mm3, [esp + i3130_fiyH]
	movq mm4, [esp + i3130_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3130_fixH], mm2
	movq [esp + i3130_fiyH], mm3
	movq [esp + i3130_fizH], mm4	

	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i3130_innerk]
	jz  .i3130_updateouterdata
	jmp .i3130_inner_loop	
.i3130_updateouterdata:	
	mov   ecx, [esp + i3130_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3130_fixO]
	pfadd mm7, [esp + i3130_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i3130_fixH]
	movq  mm3, [esp + i3130_fiyH]
	movq  mm1, [esp + i3130_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 

	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i3130_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3130_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3130_fixO]
	pfadd mm7, [esp + i3130_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i3130_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3130_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3130_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + i3130_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i3130_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + i3130_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnbtot[gid] 
	;# finish if last 
	dec dword ptr [ebp + i3130_nri]
	jz  .i3130_end
	;# not last, iterate once more! 
	jmp .i3130_outer
.i3130_end:
	femms
	add esp, 232
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret


.globl inl3300_3dnow
.globl _inl3300_3dnow
inl3300_3dnow:	
_inl3300_3dnow:	
.equiv		i3300_nri,		8
.equiv		i3300_iinr,		12
.equiv		i3300_jindex,		16
.equiv		i3300_jjnr,		20
.equiv		i3300_shift,		24
.equiv		i3300_shiftvec,		28
.equiv		i3300_fshift,		32
.equiv		i3300_gid,		36
.equiv		i3300_pos,		40		
.equiv		i3300_faction,		44
.equiv		i3300_charge,		48
.equiv		i3300_facel,		52
.equiv		i3300_Vc,		56			
.equiv		i3300_type,		60
.equiv		i3300_ntype,		64
.equiv		i3300_nbfp,		68	
.equiv		i3300_Vnb,		72
.equiv		i3300_tabscale,		76
.equiv		i3300_VFtab,		80
	;# stack offsets for local variables 
.equiv		i3300_is3,		0
.equiv		i3300_ii3,		4
.equiv		i3300_ix,		8
.equiv		i3300_iy,		12
.equiv		i3300_iz,		16
.equiv		i3300_iq,		20  
.equiv		i3300_vctot,		28  
.equiv		i3300_vnbtot,		36  
.equiv		i3300_c6,		44  
.equiv		i3300_c12,		52  
.equiv		i3300_two,		60  
.equiv		i3300_n1,		68  
.equiv		i3300_tsc,		76  
.equiv		i3300_ntia,		84
.equiv		i3300_innerjjnr,	88
.equiv		i3300_innerk,		92		
.equiv		i3300_fix,		96
.equiv		i3300_fiy,		100
.equiv		i3300_fiz,		104
.equiv		i3300_dx1,		108
.equiv		i3300_dy1,		112
.equiv		i3300_dz1,		116
.equiv		i3300_dx2,		120
.equiv		i3300_dy2,		124
.equiv		i3300_dz2,		128						
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 132		;# local stack space 
	femms
	;# move data to local stack  
	movq  mm0, [mm_two]
	movd  mm3, [ebp + i3300_tabscale]
	movq  [esp + i3300_two],    mm0
	punpckldq mm3,mm3
	movq  [esp + i3300_tsc], mm3	
	;# assume we have at least one i particle - start directly 	
.i3300_outer:
	mov   eax, [ebp + i3300_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3300_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3300_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3300_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + i3300_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3300_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3300_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i3300_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3300_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + i3300_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i3300_ntype]
	shl   edx, 1
	mov   [esp + i3300_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3300_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3300_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + i3300_ix], mm0	
	movd  [esp + i3300_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + i3300_vctot],  mm7
	movq  [esp + i3300_vnbtot], mm7
	movq  [esp + i3300_fix],    mm7
	movd  [esp + i3300_fiz],    mm7

	mov   eax, [ebp + i3300_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3300_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + i3300_pos]
	mov   edi, [ebp + i3300_faction]	
	mov   eax, [ebp + i3300_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3300_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + i3300_innerk], edx    ;# number of innerloop atoms 
	jge   .i3300_unroll_loop
	jmp   .i3300_finish_inner
.i3300_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + i3300_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3300_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3300_charge]    ;# base of charge[] 
	movq mm5, [esp + i3300_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + i3300_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i3300_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i3300_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i3300_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i3300_c6], mm5
	movq [esp + i3300_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3300_pos]

	movq  mm0, [esp + i3300_ix]
	movd  mm1, [esp + i3300_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3300_dx1], mm4	     ;# store dr 
	movd  [esp + i3300_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3300_dx2], mm6	     ;# store dr 
	movd  [esp + i3300_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i3300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3300_VFtab]
	mov ecx, [esp + i3300_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3300_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3300_vctot]      ;# add the earlier value 
	movq [esp + i3300_vctot], mm5       ;# store the sum       

	;# dispersion table 
	mov ecx, [esp + i3300_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + i3300_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i3300_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	pfadd mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i3300_vnbtot]      ;# add the earlier value 
	movq [esp + i3300_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + i3300_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + i3300_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i3300_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3300_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3300_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3300_dx1]	;# fetch dr 
	movd mm3,  [esp + i3300_dz1]

	;# update vnbtot 
	pfadd mm5, [esp + i3300_vnbtot]      ;# add the earlier value 
	movq [esp + i3300_vnbtot], mm5       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3300_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3300_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3300_fix]
	movd mm1,  [esp + i3300_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3300_fix], mm0
	movd [esp + i3300_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3300_innerk],  2
	jl    .i3300_finish_inner
	jmp   .i3300_unroll_loop
.i3300_finish_inner:	
	and dword ptr [esp + i3300_innerk],  1
	jnz  .i3300_single_inner
	jmp  .i3300_updateouterdata		
.i3300_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3300_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i3300_charge]
	movd mm5, [esp + i3300_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + i3300_nbfp]
	mov ecx, [ebp + i3300_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i3300_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i3300_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i3300_c12], mm5

	mov   esi, [ebp + i3300_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3300_ix]
	movd  mm1, [esp + i3300_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3300_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3300_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3300_VFtab]
	mov ecx, [esp + i3300_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3300_vctot]      ;# add the earlier value 
	movq [esp + i3300_vctot], mm5       ;# store the sum       
	
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i3300_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	pfadd mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i3300_vnbtot]      ;# add the earlier value 
	movq [esp + i3300_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3300_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i3300_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3300_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# update vnbtot 
	pfadd mm5, [esp + i3300_vnbtot]      ;# add the earlier value 
	movq [esp + i3300_vnbtot], mm5       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3300_dx1]
	movd mm3,  [esp + i3300_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3300_fix]
	movd mm1,  [esp + i3300_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3300_fix], mm0
	movd [esp + i3300_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3300_updateouterdata:	
	mov   ecx, [esp + i3300_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3300_fix]
	pfadd mm7, [esp + i3300_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3300_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3300_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3300_fix]
	pfadd mm7, [esp + i3300_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	mov   edx, [ebp + i3300_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3300_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3300_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3300_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i3300_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3300_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + i3300_nri]
	dec ecx
	jecxz .i3300_end
	;# not last, iterate once more! 
	mov [ebp + i3300_nri], ecx
	jmp .i3300_outer
.i3300_end:
	femms
	add esp, 132
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret





.globl inl3310_3dnow
.globl _inl3310_3dnow
inl3310_3dnow:	
_inl3310_3dnow:	
.equiv		i3310_nri,		8
.equiv		i3310_iinr,		12
.equiv		i3310_jindex,		16
.equiv		i3310_jjnr,		20
.equiv		i3310_shift,		24
.equiv		i3310_shiftvec,		28
.equiv		i3310_fshift,		32
.equiv		i3310_gid,		36
.equiv		i3310_pos,		40		
.equiv		i3310_faction,		44
.equiv		i3310_charge,		48
.equiv		i3310_facel,		52
.equiv		i3310_Vc,		56			
.equiv		i3310_type,		60
.equiv		i3310_ntype,		64
.equiv		i3310_nbfp,		68	
.equiv		i3310_Vnb,		72
.equiv		i3310_tabscale,		76
.equiv		i3310_VFtab,		80
.equiv		i3310_nsatoms,		84		
	;# stack offsets for local variables 
.equiv		i3310_is3,		0
.equiv		i3310_ii3,		4
.equiv		i3310_shX,		8
.equiv		i3310_shY,		12 
.equiv		i3310_shZ,		16	
.equiv		i3310_ix,		20
.equiv		i3310_iy,		24
.equiv		i3310_iz,		28	
.equiv		i3310_iq,		32  
.equiv		i3310_vctot,		40  
.equiv		i3310_vnbtot,		48  
.equiv		i3310_c6,		56  
.equiv		i3310_c12,		64  
.equiv		i3310_two,		72  
.equiv		i3310_n1,		80  
.equiv		i3310_tsc,		88  
.equiv		i3310_ntia,		96	
.equiv		i3310_innerjjnr0,	100
.equiv		i3310_innerk0,		104	
.equiv		i3310_innerjjnr,	108
.equiv		i3310_innerk,		112	
.equiv		i3310_fix,		116
.equiv		i3310_fiy,		120
.equiv		i3310_fiz,		124
.equiv		i3310_dx1,		128
.equiv		i3310_dy1,		132
.equiv		i3310_dz1,		136
.equiv		i3310_dx2,		140
.equiv		i3310_dy2,		144
.equiv		i3310_dz2,		148
.equiv		i3310_nsvdwc,		152
.equiv		i3310_nscoul,		156
.equiv		i3310_nsvdw,		160
.equiv		i3310_solnr,		164		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 168		;# local stack space 
	femms
	movq  mm0, [mm_two]
	movd  mm3, [ebp + i3310_tabscale]
	movq  [esp + i3310_two],    mm0
	punpckldq mm3,mm3
	movq  [esp + i3310_tsc], mm3	
	;# assume we have at least one i particle - start directly 		
.i3310_outer:
	mov   eax, [ebp + i3310_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3310_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3310_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3310_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + i3310_shX], mm0
	movd  [esp + i3310_shZ], mm1

	mov   ecx, [ebp + i3310_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3310_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + i3310_nsatoms]
	add dword ptr [ebp + i3310_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + i3310_nsvdwc], edx
	mov   [esp + i3310_nscoul], eax
	mov   [esp + i3310_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + i3310_vctot],  mm7
	movq  [esp + i3310_vnbtot], mm7
	mov   [esp + i3310_solnr],  ebx
	
	mov   eax, [ebp + i3310_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3310_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + i3310_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3310_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + i3310_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + i3310_pos]
	mov   edi, [ebp + i3310_faction]
	
	mov   ecx, [esp + i3310_nsvdwc]
	cmp   ecx,  0
	jnz   .i3310_mno_vdwc
	jmp   .i3310_testcoul
.i3310_mno_vdwc:
	mov   ebx,  [esp + i3310_solnr]
	inc   dword ptr [esp + i3310_solnr]
	mov   edx, [ebp + i3310_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i3310_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3310_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + i3310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i3310_ntype]
	shl   edx, 1
	mov   [esp + i3310_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3310_pos]    ;# eax = base of pos[] 
	mov   [esp + i3310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i3310_shX]
	pfadd mm1, [esp + i3310_shZ]
	movq  [esp + i3310_ix], mm0	
	movd  [esp + i3310_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i3310_fix],   mm7
	movd  [esp + i3310_fiz],   mm7

	mov   ecx, [esp + i3310_innerjjnr0]
	mov   [esp + i3310_innerjjnr], ecx
	mov   edx, [esp + i3310_innerk0]
    sub   edx,  2
    mov   [esp + i3310_innerk], edx    ;# number of innerloop atoms 
	jge   .i3310_unroll_vdwc_loop
	jmp   .i3310_finish_vdwc_inner
.i3310_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i3310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3310_charge]    ;# base of charge[] 
	movq mm5, [esp + i3310_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + i3310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i3310_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i3310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i3310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i3310_c6], mm5
	movq [esp + i3310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3310_pos]

	movq  mm0, [esp + i3310_ix]
	movd  mm1, [esp + i3310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3310_dx1], mm4	     ;# store dr 
	movd  [esp + i3310_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3310_dx2], mm6	     ;# store dr 
	movd  [esp + i3310_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3310_VFtab]
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3310_vctot]      ;# add the earlier value 
	movq [esp + i3310_vctot], mm5       ;# store the sum       

	;# dispersion table 
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + i3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i3310_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	pfadd mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i3310_vnbtot]      ;# add the earlier value 
	movq [esp + i3310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + i3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i3310_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3310_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3310_dx1]	;# fetch dr 
	movd mm3,  [esp + i3310_dz1]

	;# update vnbtot 
	pfadd mm5, [esp + i3310_vnbtot]      ;# add the earlier value 
	movq [esp + i3310_vnbtot], mm5       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3310_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3310_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3310_fix]
	movd mm1,  [esp + i3310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3310_fix], mm0
	movd [esp + i3310_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3310_innerk],  2
	jl    .i3310_finish_vdwc_inner
	jmp   .i3310_unroll_vdwc_loop
.i3310_finish_vdwc_inner:	
	and dword ptr [esp + i3310_innerk],  1
	jnz  .i3310_single_vdwc_inner
	jmp  .i3310_updateouterdata_vdwc		
.i3310_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i3310_charge]
	movd mm5, [esp + i3310_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + i3310_nbfp]
	mov ecx, [ebp + i3310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i3310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i3310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i3310_c12], mm5

	mov   esi, [ebp + i3310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3310_ix]
	movd  mm1, [esp + i3310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3310_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3310_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3310_VFtab]
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3310_vctot]      ;# add the earlier value 
	movq [esp + i3310_vctot], mm5       ;# store the sum       
	
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i3310_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	pfadd mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i3310_vnbtot]      ;# add the earlier value 
	movq [esp + i3310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i3310_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# update vnbtot 
	pfadd mm5, [esp + i3310_vnbtot]      ;# add the earlier value 
	movq [esp + i3310_vnbtot], mm5       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3310_dx1]
	movd mm3,  [esp + i3310_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3310_fix]
	movd mm1,  [esp + i3310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3310_fix], mm0
	movd [esp + i3310_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3310_updateouterdata_vdwc:	
	mov   ecx, [esp + i3310_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3310_fix]
	pfadd mm7, [esp + i3310_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3310_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3310_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3310_fix]
	pfadd mm7, [esp + i3310_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i3310_nsvdwc]
	jz  .i3310_testcoul
	jmp .i3310_mno_vdwc
.i3310_testcoul:	
	mov  ecx, [esp + i3310_nscoul]
	cmp  ecx,  0
	jnz  .i3310_mno_coul
	jmp  .i3310_testvdw
.i3310_mno_coul:
	mov   ebx,  [esp + i3310_solnr]
	inc   dword ptr [esp + i3310_solnr]
	mov   edx, [ebp + i3310_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + i3310_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3310_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3310_pos]    ;# eax = base of pos[] 
	mov   [esp + i3310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i3310_shX]
	pfadd mm1, [esp + i3310_shZ]
	movq  [esp + i3310_ix], mm0	
	movd  [esp + i3310_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i3310_fix],   mm7
	movd  [esp + i3310_fiz],   mm7

	mov   ecx, [esp + i3310_innerjjnr0]
	mov   [esp + i3310_innerjjnr], ecx
	mov   edx, [esp + i3310_innerk0]
    sub   edx,  2
    mov   [esp + i3310_innerk], edx    ;# number of innerloop atoms 
	jge   .i3310_unroll_coul_loop
	jmp   .i3310_finish_coul_inner
.i3310_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i3310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3310_charge]    ;# base of charge[] 
	movq mm5, [esp + i3310_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3310_pos]

	movq  mm0, [esp + i3310_ix]
	movd  mm1, [esp + i3310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3310_dx1], mm4	     ;# store dr 
	movd  [esp + i3310_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	movq  [esp + i3310_dx2], mm6	     ;# store dr 
	movd  [esp + i3310_dz2], mm7
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

	pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
	pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
	movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
	pfrsqit1 mm0,mm4				
	pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + i3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3310_VFtab]
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3310_vctot]      ;# add the earlier value 
	movq [esp + i3310_vctot], mm5       ;# store the sum       

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3
	pfmul mm1, [esp + i3310_tsc]	
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3310_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3310_dx1]	;# fetch dr 
	movd mm3,  [esp + i3310_dz1]

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3310_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3310_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 

	movq mm0,  [esp + i3310_fix]
	movd mm1,  [esp + i3310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3310_fix], mm0
	movd [esp + i3310_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3310_innerk],  2
	jl    .i3310_finish_coul_inner
	jmp   .i3310_unroll_coul_loop
.i3310_finish_coul_inner:	
	and dword ptr [esp + i3310_innerk],  1
	jnz  .i3310_single_coul_inner
	jmp  .i3310_updateouterdata_coul		
.i3310_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + i3310_charge]
	movd mm5, [esp + i3310_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + i3310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3310_ix]
	movd  mm1, [esp + i3310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3310_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3310_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3310_VFtab]
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	pfmul mm3, mm7	;# fijC=FF*qq  

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + i3310_vctot]      ;# add the earlier value 
	movq [esp + i3310_vctot], mm5       ;# store the sum       
	
	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3310_dx1]
	movd mm3,  [esp + i3310_dz1]


	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3310_fix]
	movd mm1,  [esp + i3310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3310_fix], mm0
	movd [esp + i3310_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3310_updateouterdata_coul:	
	mov   ecx, [esp + i3310_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3310_fix]
	pfadd mm7, [esp + i3310_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3310_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3310_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3310_fix]
	pfadd mm7, [esp + i3310_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i3310_nscoul]
	jz  .i3310_testvdw
	jmp .i3310_mno_coul
.i3310_testvdw:	
	mov  ecx, [esp + i3310_nsvdw]
	cmp  ecx,  0
	jnz  .i3310_mno_vdw
	jmp  .i3310_last_mno
.i3310_mno_vdw:
	mov   ebx,  [esp + i3310_solnr]
	inc   dword ptr [esp + i3310_solnr]

	mov   edx, [ebp + i3310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + i3310_ntype]
	shl   edx, 1
	mov   [esp + i3310_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3310_pos]    ;# eax = base of pos[] 
	mov   [esp + i3310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + i3310_shX]
	pfadd mm1, [esp + i3310_shZ]
	movq  [esp + i3310_ix], mm0	
	movd  [esp + i3310_iz], mm1	

	;# clear forces 
	pxor  mm7,mm7
	movq  [esp + i3310_fix],   mm7
	movd  [esp + i3310_fiz],   mm7

	mov   ecx, [esp + i3310_innerjjnr0]
	mov   [esp + i3310_innerjjnr], ecx
	mov   edx, [esp + i3310_innerk0]
    sub   edx,  2
    mov   [esp + i3310_innerk], edx    ;# number of innerloop atoms 
	jge   .i3310_unroll_vdw_loop
	jmp   .i3310_finish_vdw_inner
.i3310_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + i3310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + i3310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + i3310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + i3310_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + i3310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + i3310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6, mm5			
	punpckldq mm5, mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6, mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + i3310_c6], mm5
	movq [esp + i3310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + i3310_pos]

	movq  mm0, [esp + i3310_ix]
	movd  mm1, [esp + i3310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	movq  [esp + i3310_dx1], mm4	     ;# store dr 
	movd  [esp + i3310_dz1], mm5
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6, mm0	             ;# dr = ir - jr  
	pfsubr mm7, mm1
	movq  [esp + i3310_dx2], mm6	     ;# store dr 
	movd  [esp + i3310_dz2], mm7
	pfmul mm6, mm6	             ;# square dx,dy,dz 
	pfmul mm7, mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0, mm1
	punpckldq mm4, mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2, mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0, mm0
    pfrsqit1 mm0, mm4				
    pfrcpit2 mm0, mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + i3310_tsc]	;# mm1=rt 
	pf2iw mm4, mm1
	movq [esp + i3310_n1], mm4
	pi2fd mm4, mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2, mm1
	pfmul mm2, mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3310_VFtab]
	;# dispersion table 
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]	
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	mov ecx, [esp + i3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i3310_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	movq mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i3310_vnbtot]      ;# add the earlier value 
	movq [esp + i3310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + i3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i3310_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijD+ fijR 

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm1, [esp + i3310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	prefetchw [esp + i3310_dx1]	;# prefetch i forces to cache 

	;# spread fscalar to both positions 
	movq mm1,mm0
	punpckldq mm0,mm0
	punpckhdq mm1,mm1

	;# calc vector force 
	prefetchw [edi + eax*4]	;# prefetch the 1st faction to cache 
	movq mm2,  [esp + i3310_dx1]	;# fetch dr 
	movd mm3,  [esp + i3310_dz1]

	;# update vnbtot 
	pfadd mm5, [esp + i3310_vnbtot]      ;# add the earlier value 
	movq [esp + i3310_vnbtot], mm5       ;# store the sum       

	prefetchw [edi + ebx*4]	;# prefetch the 2nd faction to cache 
	pfmul mm2, mm0		;# mult by fs  
	pfmul mm3, mm0

	movq mm4,  [esp + i3310_dx2] 	;# fetch dr 
	movd mm5,  [esp + i3310_dz2]
	pfmul mm4, mm1   	;# mult by fs  
	pfmul mm5, mm1
	;# update i forces 
	
	movq mm0,  [esp + i3310_fix]
	movd mm1,  [esp + i3310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3

	pfadd mm0, mm4
	pfadd mm1, mm5
	movq [esp + i3310_fix], mm0
	movd [esp + i3310_fiz], mm1
	;# update j forces 

	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax*4 + 8]
	movq mm6,  [edi + ebx*4]
	movd mm7,  [edi + ebx*4 + 8]
	
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfsub mm6, mm4
	pfsub mm7, mm5
	
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	movq [edi + ebx*4], mm6
	movd [edi + ebx*4 + 8], mm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + i3310_innerk],  2
	jl    .i3310_finish_vdw_inner
	jmp   .i3310_unroll_vdw_loop
.i3310_finish_vdw_inner:	
	and dword ptr [esp + i3310_innerk],  1
	jnz  .i3310_single_vdw_inner
	jmp  .i3310_updateouterdata_vdw		
.i3310_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + i3310_nbfp]
	mov ecx, [ebp + i3310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + i3310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i3310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i3310_c12], mm5

	mov   esi, [ebp + i3310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + i3310_ix]
	movd  mm1, [esp + i3310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	movq  [esp + i3310_dx1], mm4
	pfmul mm4,mm4
	movd  [esp + i3310_dz1], mm5	
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + i3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + i3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 n0. 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + i3310_VFtab]
	mov ecx, [esp + i3310_n1]
	lea ecx, [ecx + ecx*2]	
	shl ecx, 2
	;# dispersion table
	 ;# load all the table values we need
	 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i3310_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	movq mm3, mm7	;# add to fscal 

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i3310_vnbtot]      ;# add the earlier value 
	movq [esp + i3310_vnbtot], mm5       ;# store the sum       

	;# repulsion table
	 ;# load all the table values we need
	 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3310_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i3310_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of mm3 
	pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + i3310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# update vnbtot 
	pfadd mm5, [esp + i3310_vnbtot]      ;# add the earlier value 
	movq [esp + i3310_vnbtot], mm5       ;# store the sum       

	;# spread fscalar to both positions 
	punpckldq mm0,mm0
	;# calc vectorial force 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm2,  [esp + i3310_dx1]
	movd mm3,  [esp + i3310_dz1]

	pfmul mm2, mm0
	pfmul mm3, mm0

	;# update i particle force 
	movq mm0,  [esp + i3310_fix]
	movd mm1,  [esp + i3310_fiz]
	pfadd mm0, mm2
	pfadd mm1, mm3
	movq [esp + i3310_fix], mm0
	movd [esp + i3310_fiz], mm1
	;# update j particle force 
	movq mm0,  [edi + eax*4]
	movd mm1,  [edi + eax *4+ 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	movq [edi + eax*4], mm0
	movd [edi + eax*4 +8], mm1
	;# done! 
.i3310_updateouterdata_vdw:	
	mov   ecx, [esp + i3310_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment i force 
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3310_fix]
	pfadd mm7, [esp + i3310_fiz]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	mov   ebx, [ebp + i3310_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3310_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3310_fix]
	pfadd mm7, [esp + i3310_fiz]
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7

	;# loop back to mno 
	dec dword ptr [esp + i3310_nsvdw]
	jz  .i3310_last_mno
	jmp .i3310_mno_vdw
	
.i3310_last_mno:	
	mov   edx, [ebp + i3310_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3310_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3310_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3310_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i3310_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3310_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + i3310_nri]
	dec ecx
	jecxz .i3310_end
	;# not last, iterate once more! 
	mov [ebp + i3310_nri], ecx
	jmp .i3310_outer
.i3310_end:
	femms
	add esp, 168
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret


.globl inl3320_3dnow
.globl _inl3320_3dnow
inl3320_3dnow:	
_inl3320_3dnow:	
.equiv		i3320_nri,		8
.equiv		i3320_iinr,		12
.equiv		i3320_jindex,		16
.equiv		i3320_jjnr,		20
.equiv		i3320_shift,		24
.equiv		i3320_shiftvec,		28
.equiv		i3320_fshift,		32
.equiv		i3320_gid,		36
.equiv		i3320_pos,		40		
.equiv		i3320_faction,		44
.equiv		i3320_charge,		48
.equiv		i3320_facel,		52
.equiv		i3320_Vc,		56			
.equiv		i3320_type,		60
.equiv		i3320_ntype,		64
.equiv		i3320_nbfp,		68	
.equiv		i3320_Vnb,		72
.equiv		i3320_tabscale,		76
.equiv		i3320_VFtab,		80
			;# stack offsets for local variables 
.equiv		i3320_is3,		0
.equiv		i3320_ii3,		4
.equiv		i3320_ixO,		8
.equiv		i3320_iyO,		12
.equiv		i3320_izO,		16	
.equiv		i3320_ixH,		20  
.equiv		i3320_iyH,		28  
.equiv		i3320_izH,		36  
.equiv		i3320_iqO,		44  
.equiv		i3320_iqH,		52  
.equiv		i3320_qqO,		60  
.equiv		i3320_qqH,		68  
.equiv		i3320_vctot,		76  
.equiv		i3320_vnbtot,		84  
.equiv		i3320_c6,		92  
.equiv		i3320_c12,		100 
.equiv		i3320_two,		108 
.equiv		i3320_n1,		116 
.equiv		i3320_tsc,		124 
.equiv		i3320_ntia,		132 
.equiv		i3320_innerjjnr,	140
.equiv		i3320_innerk,		144	
.equiv		i3320_fixO,		148
.equiv		i3320_fiyO,		152
.equiv		i3320_fizO,		156
.equiv		i3320_fixH,		160 
.equiv		i3320_fiyH,		168 
.equiv		i3320_fizH,		176 
.equiv		i3320_dxO,		184
.equiv		i3320_dyO,		188
.equiv		i3320_dzO,		192
.equiv		i3320_dxH,		196 
.equiv		i3320_dyH,		204 
.equiv		i3320_dzH,		212 
.equiv		i3320_tmprsqH,		220 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 228		;# local stack space 
	femms

	mov   ecx, [ebp + i3320_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3320_charge]
	movd  mm1, [ebp + i3320_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + i3320_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + i3320_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + i3320_type] 	
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1		        
	imul  ecx, [ebp + i3320_ntype]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + i3320_ntia], ecx
	 	
	movq  mm3, [mm_two]
	movq  mm4, [ebp + i3320_tabscale]
	punpckldq mm4,mm4	    ;# spread to both halves 
	movq  [esp + i3320_two],    mm3
	movq  [esp + i3320_tsc], mm4	      
	;# assume we have at least one i particle - start directly 	 
.i3320_outer:
	mov   eax, [ebp + i3320_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3320_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3320_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3320_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i3320_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3320_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3320_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3320_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + i3320_ixO], mm5	
	movq  [esp + i3320_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i3320_ixH], mm0	
	movq [esp + i3320_iyH], mm1	
	movq [esp + i3320_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i3320_vctot], mm7
	movq  [esp + i3320_vnbtot], mm7
	movq  [esp + i3320_fixO],   mm7
	movd  [esp + i3320_fizO],   mm7
	movq  [esp + i3320_fixH],   mm7
	movq  [esp + i3320_fiyH],   mm7
	movq  [esp + i3320_fizH],   mm7

	mov   eax, [ebp + i3320_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3320_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i3320_innerk], edx        

	mov   esi, [ebp + i3320_pos]
	mov   edi, [ebp + i3320_faction]	
	mov   eax, [ebp + i3320_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3320_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i3320_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3320_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i3320_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + i3320_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + i3320_iqO]
	pfmul mm7, [esp + i3320_iqH]	;# mm6=qqO, mm7=qqH 
	movd [esp + i3320_qqO], mm6
	movq [esp + i3320_qqH], mm7

	mov ecx, [ebp + i3320_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + i3320_nbfp]
	shl edx, 1
	add edx, [esp + i3320_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + i3320_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + i3320_c12], mm5	
			
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3320_ixO]
	pfsubr mm1, [esp + i3320_izO]
		
	movq  [esp + i3320_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3320_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3320_ixH]
	pfsubr mm3, [esp + i3320_iyH]
	pfsubr mm4, [esp + i3320_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3320_dxH], mm2
	movq [esp + i3320_dyH], mm3
	movq [esp + i3320_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3320_tmprsqH], mm3
	
    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + i3320_tsc]
	pf2iw mm4, mm0
	movd [esp + i3320_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3320_VFtab]
	mov ecx, [esp + i3320_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3320_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3320_qqO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3320_qqO]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + i3320_vctot]
	movq [esp + i3320_vctot], mm5
	movq mm3, mm7
	
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3320_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i3320_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	pfadd mm3, mm7	;# add to fscal  

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i3320_vnbtot]      ;# add the earlier value 
	movq [esp + i3320_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3320_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i3320_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of fscal and multiply with rinv  
	pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + i3320_tsc]
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 
	
	;# update vnbtot  
	pfadd mm5, [esp + i3320_vnbtot]      ;# add the earlier value 
	movq [esp + i3320_vnbtot], mm5       ;# store the sum       
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3. 
	;# now do the two hydrogens. 
	movq mm0, [esp + i3320_tmprsqH] ;# mm0=rsqH 

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3320_tsc]
	pf2iw mm4, mm0
	movq [esp + i3320_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3320_VFtab]
	mov ecx, [esp + i3320_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3320_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3320_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3320_qqH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3320_qqH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3320_vctot]
	movq [esp + i3320_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
	pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3320_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	prefetchw [edi + eax*4]	;# prefetch faction to cache  
	movq mm0,  [esp + i3320_dxO]
	movd mm1,  [esp + i3320_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3320_dxH]
	movq mm6, [esp + i3320_dyH]
	movq mm7, [esp + i3320_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3320_fixO]
	movd mm3,  [esp + i3320_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3320_fixO], mm2
	movd [esp + i3320_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3320_fixH]
	movq mm3, [esp + i3320_fiyH]
	movq mm4, [esp + i3320_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3320_fixH], mm2
	movq [esp + i3320_fiyH], mm3
	movq [esp + i3320_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i3320_innerk]
	jz  .i3320_updateouterdata
	jmp .i3320_inner_loop
.i3320_updateouterdata:	
	mov   ecx, [esp + i3320_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3320_fixO]
	pfadd mm7, [esp + i3320_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i3320_fixH]
	movq  mm3, [esp + i3320_fiyH]
	movq  mm1, [esp + i3320_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 
	
	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i3320_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3320_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3320_fixO]
	pfadd mm7, [esp + i3320_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i3320_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3320_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3320_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + i3320_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i3320_vnbtot]     
	pfacc mm7,mm7	          ;# same for Vnb 
	
	mov   eax, [ebp + i3320_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 
	;# finish if last 
	dec dword ptr [ebp + i3320_nri]
	jz  .i3320_end
	;# not last, iterate once more! 
	jmp .i3320_outer
.i3320_end:
	femms
	add esp, 228
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret

	

.globl inl3330_3dnow
.globl _inl3330_3dnow
inl3330_3dnow:	
_inl3330_3dnow:	
.equiv		i3330_nri,		8
.equiv		i3330_iinr,		12
.equiv		i3330_jindex,		16
.equiv		i3330_jjnr,		20
.equiv		i3330_shift,		24
.equiv		i3330_shiftvec,		28
.equiv		i3330_fshift,		32
.equiv		i3330_gid,		36
.equiv		i3330_pos,		40		
.equiv		i3330_faction,		44
.equiv		i3330_charge,		48
.equiv		i3330_facel,		52
.equiv		i3330_Vc,		56			
.equiv		i3330_type,		60
.equiv		i3330_ntype,		64
.equiv		i3330_nbfp,		68	
.equiv		i3330_Vnb,		72
.equiv		i3330_tabscale,		76
.equiv		i3330_VFtab,		80
			;# stack offsets for local variables 
.equiv		i3330_is3,		0
.equiv		i3330_ii3,		4
.equiv		i3330_ixO,		8
.equiv		i3330_iyO,		12
.equiv		i3330_izO,		16	
.equiv		i3330_ixH,		20  
.equiv		i3330_iyH,		28  
.equiv		i3330_izH,		36  
.equiv		i3330_qqOO,		44  
.equiv		i3330_qqOH,		52  
.equiv		i3330_qqHH,		60  
.equiv		i3330_c6,		68  
.equiv		i3330_c12,		76  
.equiv		i3330_two,		84  
.equiv		i3330_n1,		92  
.equiv		i3330_tsc,		100 
.equiv		i3330_vctot,		108 
.equiv		i3330_vnbtot,		116 
.equiv		i3330_innerjjnr,	124
.equiv		i3330_innerk,		128	
.equiv		i3330_fixO,		132
.equiv		i3330_fiyO,		136
.equiv		i3330_fizO,		140
.equiv		i3330_fixH,		144 
.equiv		i3330_fiyH,		152 
.equiv		i3330_fizH,		160 
.equiv		i3330_dxO,		168
.equiv		i3330_dyO,		172
.equiv		i3330_dzO,		176
.equiv		i3330_dxH,		180  
.equiv		i3330_dyH,		188  
.equiv		i3330_dzH,		196  
.equiv		i3330_tmprsqH,		204  
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 212		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + i3330_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + i3330_charge]
	movd  mm1, [ebp + i3330_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + i3330_qqOO], mm4
	movq  [esp + i3330_qqOH], mm5
	movq  [esp + i3330_qqHH], mm6
	mov   edx, [ebp + i3330_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + i3330_ntype]
	add   edx, ecx
	mov   eax, [ebp + i3330_nbfp]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + i3330_c6], mm0
	movq  [esp + i3330_c12], mm1
	movq  mm2, [mm_two]
	movq  [esp + i3330_two], mm2
	movd  mm3, [ebp + i3330_tabscale]
	punpckldq mm3,mm3
	movq  [esp + i3330_tsc], mm3
.i3330_outer:
	mov   eax, [ebp + i3330_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + i3330_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + i3330_is3],ebx    	;# store is3 

	mov   eax, [ebp + i3330_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + i3330_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + i3330_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + i3330_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + i3330_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + i3330_ixO], mm5	
	movq  [esp + i3330_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + i3330_ixH], mm0	
	movq [esp + i3330_iyH], mm1	
	movq [esp + i3330_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + i3330_vctot], mm7
	movq  [esp + i3330_vnbtot], mm7
	movq  [esp + i3330_fixO],  mm7
	movq  [esp + i3330_fizO],  mm7
	movq  [esp + i3330_fixH],  mm7
	movq  [esp + i3330_fiyH],  mm7
	movq  [esp + i3330_fizH],  mm7

	mov   eax, [ebp + i3330_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + i3330_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + i3330_innerk], edx        

	mov   esi, [ebp + i3330_pos]
	mov   edi, [ebp + i3330_faction]	
	mov   eax, [ebp + i3330_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + i3330_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.i3330_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + i3330_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + i3330_innerjjnr],  4 ;# advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3330_ixO]
	pfsubr mm1, [esp + i3330_izO]
		
	movq  [esp + i3330_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3330_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3330_ixH]
	pfsubr mm3, [esp + i3330_iyH]
	pfsubr mm4, [esp + i3330_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3330_dxH], mm2
	movq [esp + i3330_dyH], mm3
	movq [esp + i3330_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3330_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	pfmul mm0, mm1		;# mm0=rsq  

	pfmul mm0, [esp + i3330_tsc]
	pf2iw mm4, mm0
	movd [esp + i3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3330_VFtab]
	mov ecx, [esp + i3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3330_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3330_qqOO]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3330_qqOO]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + i3330_vctot]
	movq [esp + i3330_vctot], mm5
	movq mm3, mm7

	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3330_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + i3330_c6]
	pfmul mm7, mm4	;# fijD 
	pfmul mm5, mm4	;# vnb6            
	pfadd mm3, mm7	;# add to fscal  

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + i3330_vnbtot]      ;# add the earlier value 
	movq [esp + i3330_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm7, [esp + i3330_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + i3330_c12]
	pfmul mm7, mm6	;# fijR 
	pfmul mm5, mm6	;# vnb12 
	pfadd mm3, mm7	;# total fscal fijC+ fijD+ fijR 

	;# change sign of fscal and multiply with rinv  
    pxor mm0,mm0
	pfsubr mm3, mm0	
	pfmul mm3, [esp + i3330_tsc]
 	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 
	
	;# update vnbtot  
	pfadd mm5, [esp + i3330_vnbtot]      ;# add the earlier value 
	movq [esp + i3330_vnbtot], mm5       ;# store the sum       
	
	;# Ready with the oxygen - potential is updated, fscal is in mm3.
	 ;# time for hydrogens!
         
	
	movq mm0, [esp + i3330_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3330_tsc]
	pf2iw mm4, mm0
	movq [esp + i3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3330_VFtab]
	mov ecx, [esp + i3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3330_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3330_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3330_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3330_qqOH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3330_vctot]
	movq [esp + i3330_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
    pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3330_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	
	
	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3330_dxO]
	movd mm1,  [esp + i3330_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3330_dxH]
	movq mm6, [esp + i3330_dyH]
	movq mm7, [esp + i3330_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3330_fixO]
	movd mm3,  [esp + i3330_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3330_fixO], mm2
	movd [esp + i3330_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3330_fixH]
	movq mm3, [esp + i3330_fiyH]
	movq mm4, [esp + i3330_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3330_fixH], mm2
	movq [esp + i3330_fiyH], mm3
	movq [esp + i3330_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4]
	movd mm3,  [edi + eax*4 + 8]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4], mm2
	movd [edi + eax*4 +8], mm3

	;# interactions with j H1 

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + i3330_ixO]
	pfsubr mm1, [esp + i3330_izO]
		
	movq  [esp + i3330_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3330_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3330_ixH]
	pfsubr mm3, [esp + i3330_iyH]
	pfsubr mm4, [esp + i3330_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3330_dxH], mm2
	movq [esp + i3330_dyH], mm3
	movq [esp + i3330_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3330_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1		;# mm0=rsq  
	
	pfmul mm0, [esp + i3330_tsc]
	pf2iw mm4, mm0
	movd [esp + i3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3330_VFtab]
	mov ecx, [esp + i3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3330_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3330_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3330_qqOH]	;# fijC=qq*FF 

	;# update vctot directly, force is moved to mm3. 
	pfadd mm5, [esp + i3330_vctot]
	movq [esp + i3330_vctot], mm5
	pxor mm3, mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + i3330_tsc]
	pfmul mm3, mm1    ;# mm3 is total fscal (for the oxygen) now 

	movq mm0, [esp + i3330_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3330_tsc]
	pf2iw mm4, mm0
	movq [esp + i3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3330_VFtab]
	mov ecx, [esp + i3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3330_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3330_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3330_qqHH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3330_qqHH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3330_vctot]
	movq [esp + i3330_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
    pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3330_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 		

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3330_dxO]
	movd mm1,  [esp + i3330_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3330_dxH]
	movq mm6, [esp + i3330_dyH]
	movq mm7, [esp + i3330_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3330_fixO]
	movd mm3,  [esp + i3330_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3330_fixO], mm2
	movd [esp + i3330_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3330_fixH]
	movq mm3, [esp + i3330_fiyH]
	movq mm4, [esp + i3330_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3330_fixH], mm2
	movq [esp + i3330_fiyH], mm3
	movq [esp + i3330_fizH], mm4
	
	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 12]
	movd mm3,  [edi + eax*4 + 20]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 12], mm2
	movd [edi + eax*4 + 20], mm3

	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + i3330_ixO]
	pfsubr mm1, [esp + i3330_izO]
		
	movq  [esp + i3330_dxO], mm0
	pfmul mm0,mm0
	movd  [esp + i3330_dzO], mm1	
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + i3330_ixH]
	pfsubr mm3, [esp + i3330_iyH]
	pfsubr mm4, [esp + i3330_izH] ;# mm2-mm4 is dxH-dzH 
	
	movq [esp + i3330_dxH], mm2
	movq [esp + i3330_dyH], mm3
	movq [esp + i3330_dzH], mm4
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + i3330_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1

	pfmul mm0, [esp + i3330_tsc]
	pf2iw mm4, mm0
	movd [esp + i3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + i3330_VFtab]
	mov ecx, [esp + i3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3330_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3330_qqOH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3330_qqOH]	;# fijC=qq*FF 

	;# update vctot directly, use mm3 for fscal sum 
	pfadd mm5, [esp + i3330_vctot]
	movq [esp + i3330_vctot], mm5
	pxor mm3,mm3
	pfsub mm3, mm7
 	pfmul mm3, [esp + i3330_tsc]
 	pfmul mm3, mm1     ;# mm3 is total fscal (for the oxygen) now 	

	movq mm0, [esp + i3330_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + i3330_tsc]
	pf2iw mm4, mm0
	movq [esp + i3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + i3330_VFtab]
	mov ecx, [esp + i3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + i3330_n1 + 4];# mm5 = Fp 
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + i3330_two]	;# two*Heps2 
	pfadd mm7, mm6
	pfadd mm7, mm5	;# mm7=FF 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + i3330_qqHH]	;# vcoul=qq*VV 
	pfmul mm7, [esp + i3330_qqHH]	;# fijC=qq*FF 
	;# update vctot 
	pfadd mm5, [esp + i3330_vctot]
	movq [esp + i3330_vctot], mm5
	
	;# change sign of fijC and multiply by rinv 
    pxor mm4,mm4
	pfsub mm4, mm7	
	pfmul mm4, [esp + i3330_tsc]
 	pfmul mm4, mm1    ;# mm4 is total fscal (for the hydrogens) now 	

	;# spread oxygen fscalar to both positions 
	punpckldq mm3,mm3
	;# calc vectorial force for O 
	movq mm0,  [esp + i3330_dxO]
	movd mm1,  [esp + i3330_dzO]
	pfmul mm0, mm3
	pfmul mm1, mm3

	;# calc vectorial force for H's 
	movq mm5, [esp + i3330_dxH]
	movq mm6, [esp + i3330_dyH]
	movq mm7, [esp + i3330_dzH]
	pfmul mm5, mm4
	pfmul mm6, mm4
	pfmul mm7, mm4
	
	;# update iO particle force 
	movq mm2,  [esp + i3330_fixO]
	movd mm3,  [esp + i3330_fizO]
	pfadd mm2, mm0
	pfadd mm3, mm1
	movq [esp + i3330_fixO], mm2
	movd [esp + i3330_fizO], mm3

	;# update iH forces 
	movq mm2, [esp + i3330_fixH]
	movq mm3, [esp + i3330_fiyH]
	movq mm4, [esp + i3330_fizH]
	pfadd mm2, mm5
	pfadd mm3, mm6
	pfadd mm4, mm7
	movq [esp + i3330_fixH], mm2
	movq [esp + i3330_fiyH], mm3
	movq [esp + i3330_fizH], mm4	

	;# pack j forces from H in the same form as the oxygen force. 
	pfacc mm5, mm6		;# mm5(l)=fjx(H1+ h2) mm5(h)=fjy(H1+ h2) 
	pfacc mm7, mm7		;# mm7(l)=fjz(H1+ h2) 
	
	pfadd mm0, mm5		;# add up total force on j particle.  
	pfadd mm1, mm7

	;# update j particle force 
	movq mm2,  [edi + eax*4 + 24]
	movd mm3,  [edi + eax*4 + 32]
	pfsub mm2, mm0
	pfsub mm3, mm1
	movq [edi + eax*4 + 24], mm2
	movd [edi + eax*4 + 32], mm3
	
	;#  done  - one more? 
	dec dword ptr [esp + i3330_innerk]
	jz  .i3330_updateouterdata
	jmp .i3330_inner_loop	
.i3330_updateouterdata:	
	mov   ecx, [esp + i3330_ii3]

	movq  mm6, [edi + ecx*4]       ;# increment iO force  
	movd  mm7, [edi + ecx*4 + 8]	
	pfadd mm6, [esp + i3330_fixO]
	pfadd mm7, [esp + i3330_fizO]
	movq  [edi + ecx*4],    mm6
	movd  [edi + ecx*4 +8], mm7

	movq  mm0, [esp + i3330_fixH]
	movq  mm3, [esp + i3330_fiyH]
	movq  mm1, [esp + i3330_fizH]
	movq  mm2, mm0
	punpckldq mm0, mm3	;# mm0(l)=fxH1, mm0(h)=fyH1 
	punpckhdq mm2, mm3	;# mm2(l)=fxH2, mm2(h)=fyH2 
	movq mm3, mm1
	pswapd mm3,mm3		
	;# mm1 is fzH1 
	;# mm3 is fzH2 

	movq  mm6, [edi + ecx*4 + 12]       ;# increment iH1 force  
	movd  mm7, [edi + ecx*4 + 20] 	
	pfadd mm6, mm0
	pfadd mm7, mm1
	movq  [edi + ecx*4 + 12],  mm6
	movd  [edi + ecx*4 + 20],  mm7
	
	movq  mm6, [edi + ecx*4 + 24]       ;# increment iH2 force 
	movd  mm7, [edi + ecx*4 + 32] 	
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [edi + ecx*4 + 24],  mm6
	movd  [edi + ecx*4 + 32],  mm7

	
	mov   ebx, [ebp + i3330_fshift]    ;# increment fshift force 
	mov   edx, [esp + i3330_is3]

	movq  mm6, [ebx + edx*4]	
	movd  mm7, [ebx + edx*4 + 8]	
	pfadd mm6, [esp + i3330_fixO]
	pfadd mm7, [esp + i3330_fizO]
	pfadd mm6, mm0
	pfadd mm7, mm1
	pfadd mm6, mm2
	pfadd mm7, mm3
	movq  [ebx + edx*4],     mm6
	movd  [ebx + edx*4 + 8], mm7
	
	mov   edx, [ebp + i3330_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + i3330_gid],  4  ;# advance pointer 

	movq  mm7, [esp + i3330_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + i3330_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + i3330_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + i3330_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnbtot[gid] 
	;# finish if last 
	dec dword ptr [ebp + i3330_nri]
	jz  .i3330_end
	;# not last, iterate once more! 
	jmp .i3330_outer
.i3330_end:
	femms
	add esp, 212
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	
 
	
	

.globl mcinl0100_3dnow
.globl _mcinl0100_3dnow
mcinl0100_3dnow:	
_mcinl0100_3dnow:	
.equiv		mci0100_nri,		8
.equiv		mci0100_iinr,		12
.equiv		mci0100_jindex,		16
.equiv		mci0100_jjnr,		20
.equiv		mci0100_shift,		24
.equiv		mci0100_shiftvec,	28
.equiv		mci0100_gid,		32
.equiv		mci0100_pos,		36
.equiv		mci0100_type,		40
.equiv		mci0100_ntype,		44
.equiv		mci0100_nbfp,		48
.equiv		mci0100_Vnb,		52
	;# stack offsets for local variables 
.equiv		mci0100_is3,		0
.equiv		mci0100_ii3,		4
.equiv		mci0100_ix,		8
.equiv		mci0100_iy,		12
.equiv		mci0100_iz,		16
.equiv		mci0100_vnbtot,		20  
.equiv		mci0100_c6,		28  
.equiv		mci0100_c12,		36  
.equiv		mci0100_ntia,		44
.equiv		mci0100_innerjjnr,	48
.equiv		mci0100_innerk,		52	
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 56		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	
.mci0100_outer:
	mov   eax, [ebp + mci0100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci0100_shift], 4		;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci0100_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci0100_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1. 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + mci0100_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci0100_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx =ii 

	mov   edx, [ebp + mci0100_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci0100_ntype]
	shl   edx, 1
	mov   [esp + mci0100_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci0100_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci0100_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + mci0100_ix], mm0	
	movd  [esp + mci0100_iz], mm1	
				
	;# clear total potential 
	pxor  mm7,mm7
	movq  [esp + mci0100_vnbtot], mm7

	mov   eax, [ebp + mci0100_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci0100_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + mci0100_pos]	
	mov   eax, [ebp + mci0100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci0100_innerjjnr], eax     ;#  pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + mci0100_innerk], edx    ;# number of innerloop atoms 
	jge   .mci0100_unroll_loop
	jmp   .mci0100_finish_inner
.mci0100_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + mci0100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci0100_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + mci0100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci0100_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci0100_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci0100_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci0100_c6], mm5
	movq [esp + mci0100_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci0100_pos]

	movq  mm0, [esp + mci0100_ix]
	movd  mm1, [esp + mci0100_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrcp mm0, mm4	             ;# lookup reciprocal seed  
    pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
                  	        	;# amd 3dnow N-R iteration to get full precision. 
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	 ;# do potential and fscal
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci0100_c12]
	pfmul mm4, [esp + mci0100_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot  
	pfadd mm6, [esp + mci0100_vnbtot]      ;# add the earlier value 
	movq [esp + mci0100_vnbtot], mm6       ;# store the sum 
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci0100_innerk],  2
	jl    .mci0100_finish_inner
	jmp   .mci0100_unroll_loop
.mci0100_finish_inner:	
	and dword ptr [esp + mci0100_innerk],  1
	jnz  .mci0100_single_inner
	jmp  .mci0100_updateouterdata
.mci0100_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci0100_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci0100_nbfp]
	mov ecx, [ebp + mci0100_type]
	mov edx, [ecx + eax*4]        	;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci0100_ntia]	    ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6 		
	movq [esp + mci0100_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12 		
	movq [esp + mci0100_c12], mm5

	mov   esi, [ebp + mci0100_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci0100_ix]
	movd  mm1, [esp + mci0100_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrcp mm0,mm4
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	;# mm4=invsq 
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci0100_c12]
	pfmul mm4, [esp + mci0100_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6 
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci0100_vnbtot]      ;# add the earlier value 
	movq [esp + mci0100_vnbtot], mm6       ;# store the sum   

.mci0100_updateouterdata:	
	mov   edx, [ebp + mci0100_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci0100_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci0100_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci0100_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + mci0100_nri]
	dec ecx
	jecxz .mci0100_end
	;# not last, iterate once more! 
	mov [ebp + mci0100_nri], ecx
	jmp .mci0100_outer
.mci0100_end:
	femms
	add esp, 56
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret
	



		
		
.globl mcinl0110_3dnow
.globl _mcinl0110_3dnow
mcinl0110_3dnow:	 
_mcinl0110_3dnow:	 
.equiv		mci0110_nri,		8
.equiv		mci0110_iinr,		12
.equiv		mci0110_jindex,		16
.equiv		mci0110_jjnr,		20
.equiv		mci0110_shift,		24
.equiv		mci0110_shiftvec,	28
.equiv		mci0110_gid,		32
.equiv		mci0110_pos,		36		
.equiv		mci0110_type,		40
.equiv		mci0110_ntype,		44
.equiv		mci0110_nbfp,		48	
.equiv		mci0110_Vnb,		52				
.equiv		mci0110_nsatoms,	56		
	;# stack offsets for local variables 
.equiv		mci0110_is3,		0
.equiv		mci0110_ii3,		4
.equiv		mci0110_shX,		8
.equiv		mci0110_shY,		12 
.equiv		mci0110_shZ,		16	
.equiv		mci0110_ix,			20
.equiv		mci0110_iy,			24
.equiv		mci0110_iz,			28	
.equiv		mci0110_vnbtot,		32 
.equiv		mci0110_c6,			40 
.equiv		mci0110_c12,		48 
.equiv		mci0110_ntia,		56
.equiv		mci0110_innerjjnr0,	60
.equiv		mci0110_innerk0,	64		
.equiv		mci0110_innerjjnr,	68
.equiv		mci0110_innerk,		72	
.equiv		mci0110_nsvdwc,		76
.equiv		mci0110_nscoul,		80
.equiv		mci0110_nsvdw,		84
.equiv		mci0110_solnr,		88		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 92		;# local stack space 
	femms
	
	;# assume we have at least one i particle - start directly 		
.mci0110_outer:
	mov   eax, [ebp + mci0110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci0110_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci0110_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci0110_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + mci0110_shX], mm0
	movd  [esp + mci0110_shZ], mm1

	mov   ecx, [ebp + mci0110_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci0110_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + mci0110_nsatoms]
	add dword ptr [ebp + mci0110_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + mci0110_nsvdwc], edx
	mov   [esp + mci0110_nscoul], eax
	mov   [esp + mci0110_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + mci0110_vnbtot], mm7
	mov   [esp + mci0110_solnr],  ebx
	
	mov   eax, [ebp + mci0110_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci0110_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + mci0110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci0110_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + mci0110_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + mci0110_pos]
	
	mov   ecx, [esp + mci0110_nsvdwc]
	cmp   ecx,  0
	jnz   .mci0110_mno_vdwc
	jmp   .mci0110_testvdw
.mci0110_mno_vdwc:
	mov   ebx, [esp + mci0110_solnr]
	inc   dword ptr [esp + mci0110_solnr]

	mov   edx, [ebp + mci0110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci0110_ntype]
	shl   edx, 1
	mov   [esp + mci0110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci0110_pos]    ;# eax = base of pos[] 
	mov   [esp + mci0110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci0110_shX]
	pfadd mm1, [esp + mci0110_shZ]
	movq  [esp + mci0110_ix], mm0	
	movd  [esp + mci0110_iz], mm1	

	mov   ecx, [esp + mci0110_innerjjnr0]
	mov   [esp + mci0110_innerjjnr], ecx
	mov   edx, [esp + mci0110_innerk0]
    sub   edx,  2
    mov   [esp + mci0110_innerk], edx    ;# number of innerloop atoms 
	jge   .mci0110_unroll_vdwc_loop
	jmp   .mci0110_finish_vdwc_inner
.mci0110_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci0110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci0110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + mci0110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci0110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci0110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci0110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci0110_c6], mm5
	movq [esp + mci0110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci0110_pos]

	movq  mm0, [esp + mci0110_ix]
	movd  mm1, [esp + mci0110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrcp mm0, mm4	             ;# lookup reciprocal seed  
    pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
                  	        	;# amd 3dnow N-R iteration to get full precision 
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	 ;# do potential and fscal
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci0110_c12]
	pfmul mm4, [esp + mci0110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci0110_vnbtot]      ;# add the earlier value 
	movq [esp + mci0110_vnbtot], mm6       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci0110_innerk],  2
	jl    .mci0110_finish_vdwc_inner
	jmp   .mci0110_unroll_vdwc_loop
.mci0110_finish_vdwc_inner:	
	and dword ptr [esp + mci0110_innerk],  1
	jnz  .mci0110_single_vdwc_inner
	jmp  .mci0110_updateouterdata_vdwc		
.mci0110_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci0110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci0110_nbfp]
	mov ecx, [ebp + mci0110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci0110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci0110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci0110_c12], mm5

	mov   esi, [ebp + mci0110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci0110_ix]
	movd  mm1, [esp + mci0110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrcp mm0,mm4
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	;# mm4=invsq  
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci0110_c12]
	pfmul mm4, [esp + mci0110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci0110_vnbtot]      ;# add the earlier value 
	movq [esp + mci0110_vnbtot], mm6       ;# store the sum       

.mci0110_updateouterdata_vdwc:	
	;# loop back to mno 
	dec dword ptr [esp + mci0110_nsvdwc]
	jz  .mci0110_testvdw
	jmp .mci0110_mno_vdwc
.mci0110_testvdw:	
	mov  ebx,  [esp + mci0110_nscoul]
	add  [esp + mci0110_solnr],  ebx

	mov  ecx, [esp + mci0110_nsvdw]
	cmp  ecx,  0
	jnz  .mci0110_mno_vdw
	jmp  .mci0110_last_mno
.mci0110_mno_vdw:
	mov   ebx,  [esp + mci0110_solnr]
	inc   dword ptr [esp + mci0110_solnr]

	mov   edx, [ebp + mci0110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci0110_ntype]
	shl   edx, 1
	mov   [esp + mci0110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci0110_pos]    ;# eax = base of pos[] 
	mov   [esp + mci0110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci0110_shX]
	pfadd mm1, [esp + mci0110_shZ]
	movq  [esp + mci0110_ix], mm0	
	movd  [esp + mci0110_iz], mm1	

	mov   ecx, [esp + mci0110_innerjjnr0]
	mov   [esp + mci0110_innerjjnr], ecx
	mov   edx, [esp + mci0110_innerk0]
    sub   edx,  2
    mov   [esp + mci0110_innerk], edx    ;# number of innerloop atoms 
	jge   .mci0110_unroll_vdw_loop
	jmp   .mci0110_finish_vdw_inner
.mci0110_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci0110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci0110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + mci0110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci0110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci0110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci0110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci0110_c6], mm5
	movq [esp + mci0110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci0110_pos]

	movq  mm0, [esp + mci0110_ix]
	movd  mm1, [esp + mci0110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrcp mm0, mm4	             ;# lookup reciprocal seed  
    pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
                  	        	;# amd 3dnow N-R iteration to get full precision 
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	 ;# do potential and fscal
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci0110_c12]
	pfmul mm4, [esp + mci0110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci0110_vnbtot]      ;# add the earlier value 
	movq [esp + mci0110_vnbtot], mm6       ;# store the sum       

	;# should we do one more iteration? 
	sub dword ptr [esp + mci0110_innerk],  2
	jl    .mci0110_finish_vdw_inner
	jmp   .mci0110_unroll_vdw_loop
.mci0110_finish_vdw_inner:	
	and dword ptr [esp + mci0110_innerk],  1
	jnz  .mci0110_single_vdw_inner
	jmp  .mci0110_updateouterdata_vdw		
.mci0110_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci0110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci0110_nbfp]
	mov ecx, [ebp + mci0110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci0110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci0110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci0110_c12], mm5

	mov   esi, [ebp + mci0110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci0110_ix]
	movd  mm1, [esp + mci0110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrcp mm0,mm4
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	;# mm4=invsq  
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci0110_c12]
	pfmul mm4, [esp + mci0110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci0110_vnbtot]      ;# add the earlier value 
	movq [esp + mci0110_vnbtot], mm6       ;# store the sum       

.mci0110_updateouterdata_vdw:	
	;# loop back to mno 
	dec dword ptr [esp + mci0110_nsvdw]
	jz  .mci0110_last_mno
	jmp .mci0110_mno_vdw
	
.mci0110_last_mno:	
	mov   edx, [ebp + mci0110_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci0110_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci0110_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci0110_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + mci0110_nri]
	dec ecx
	jecxz .mci0110_end
	;# not last, iterate once more! 
	mov [ebp + mci0110_nri], ecx
	jmp .mci0110_outer
.mci0110_end:
	femms
	add esp, 92
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret



.globl mcinl0300_3dnow
.globl _mcinl0300_3dnow
mcinl0300_3dnow:	
_mcinl0300_3dnow:	
.equiv		mci0300_nri,		8
.equiv		mci0300_iinr,		12
.equiv		mci0300_jindex,		16
.equiv		mci0300_jjnr,		20
.equiv		mci0300_shift,		24
.equiv		mci0300_shiftvec,	28
.equiv		mci0300_gid,		32
.equiv		mci0300_pos,		36		
.equiv		mci0300_type,		40
.equiv		mci0300_ntype,		44
.equiv		mci0300_nbfp,		48	
.equiv		mci0300_Vnb,		52
.equiv		mci0300_tabscale,	56
.equiv		mci0300_VFtab,		60
	;# stack offsets for local variables 
.equiv		mci0300_is3,		0
.equiv		mci0300_ii3,		4
.equiv		mci0300_ix,			8
.equiv		mci0300_iy,			12
.equiv		mci0300_iz,			16
.equiv		mci0300_vnbtot,		20 
.equiv		mci0300_c6,			28 
.equiv		mci0300_c12,		36
.equiv		mci0300_n1,			44 
.equiv		mci0300_tsc,		52 
.equiv		mci0300_ntia,		60
.equiv		mci0300_innerjjnr,	64
.equiv		mci0300_innerk,		68
    push ebp
    mov ebp,esp
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 72		;# local stack space 
	femms
	;# move data to local stack  
	movd  mm3, [ebp + mci0300_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci0300_tsc], mm3	
	;# assume we have at least one i particle - start directly 	
.mci0300_outer:
	mov   eax, [ebp + mci0300_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci0300_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci0300_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci0300_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + mci0300_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci0300_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci0300_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci0300_ntype]
	shl   edx, 1
	mov   [esp + mci0300_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci0300_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci0300_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + mci0300_ix], mm0	
	movd  [esp + mci0300_iz], mm1	
				
	;# clear total potential 
	pxor  mm7,mm7
	movq  [esp + mci0300_vnbtot], mm7

	mov   eax, [ebp + mci0300_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci0300_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + mci0300_pos]
	mov   eax, [ebp + mci0300_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci0300_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + mci0300_innerk], edx    ;# number of innerloop atoms 
	jge   .mci0300_unroll_loop
	jmp   .mci0300_finish_inner
.mci0300_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + mci0300_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci0300_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best  
	
	mov ecx, [ebp + mci0300_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci0300_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci0300_ntia]	     ;# tja = ntia + 2*type  
	add ecx, [esp + mci0300_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci0300_c6], mm5
	movq [esp + mci0300_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci0300_pos]

	movq  mm0, [esp + mci0300_ix]
	movd  mm1, [esp + mci0300_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 
	;# do potential and fscal 
	pfmul mm1, [esp + mci0300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci0300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci0300_VFtab]
	;# dispersion table 
	mov ecx, [esp + mci0300_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci0300_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci0300_c6]
	pfmul mm5, mm4	;# vnb6        
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci0300_vnbtot]      ;# add the earlier value 
	movq [esp + mci0300_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + mci0300_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + mci0300_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci0300_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci0300_vnbtot]      ;# add the earlier value 
	movq [esp + mci0300_vnbtot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci0300_innerk],  2
	jl    .mci0300_finish_inner
	jmp   .mci0300_unroll_loop
.mci0300_finish_inner:	
	and dword ptr [esp + mci0300_innerk],  1
	jnz  .mci0300_single_inner
	jmp  .mci0300_updateouterdata		
.mci0300_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci0300_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci0300_nbfp]
	mov ecx, [ebp + mci0300_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci0300_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci0300_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci0300_c12], mm5

	mov   esi, [ebp + mci0300_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci0300_ix]
	movd  mm1, [esp + mci0300_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci0300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci0300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci0300_VFtab]
	mov ecx, [esp + mci0300_n1]
	shl ecx, 3
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci0300_c6]
	pfmul mm5, mm4	;# vnb6            
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci0300_vnbtot]      ;# add the earlier value 
	movq [esp + mci0300_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci0300_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci0300_vnbtot]      ;# add the earlier value 
	movq [esp + mci0300_vnbtot], mm5       ;# store the sum       

.mci0300_updateouterdata:	
	mov   edx, [ebp + mci0300_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci0300_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci0300_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci0300_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + mci0300_nri]
	dec ecx
	jecxz .mci0300_end
	;# not last, iterate once more! 
	mov [ebp + mci0300_nri], ecx
	jmp .mci0300_outer
.mci0300_end:
	femms
	add esp, 72
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


			
	
.globl mcinl0310_3dnow
.globl _mcinl0310_3dnow
mcinl0310_3dnow:	
_mcinl0310_3dnow:	
.equiv		mci0310_nri,		8
.equiv		mci0310_iinr,		12
.equiv		mci0310_jindex,		16
.equiv		mci0310_jjnr,		20
.equiv		mci0310_shift,		24
.equiv		mci0310_shiftvec,	28
.equiv		mci0310_gid,		32
.equiv		mci0310_pos,		36		
.equiv		mci0310_type,		40
.equiv		mci0310_ntype,		44
.equiv		mci0310_nbfp,		48	
.equiv		mci0310_Vnb,		52
.equiv		mci0310_tabscale,	56
.equiv		mci0310_VFtab,		60
.equiv		mci0310_nsatoms,	64		
	;# stack offsets for local variables 
.equiv		mci0310_is3,		0
.equiv		mci0310_ii3,		4
.equiv		mci0310_shX,		8
.equiv		mci0310_shY,		12 
.equiv		mci0310_shZ,		16	
.equiv		mci0310_ix,		20
.equiv		mci0310_iy,		24
.equiv		mci0310_iz,		28	
.equiv		mci0310_vnbtot,		32 
.equiv		mci0310_c6,		40 
.equiv		mci0310_c12,		48
.equiv		mci0310_n1,		56 
.equiv		mci0310_tsc,		64 
.equiv		mci0310_ntia,		72	
.equiv		mci0310_innerjjnr0,	76
.equiv		mci0310_innerk0,	80		
.equiv		mci0310_innerjjnr,	84
.equiv		mci0310_innerk,		88					
.equiv		mci0310_nsvdwc,		92
.equiv		mci0310_nscoul,		96
.equiv		mci0310_nsvdw,		100
.equiv		mci0310_solnr,		104		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 108		;# local stack space 
	femms
	movd  mm3, [ebp + mci0310_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci0310_tsc], mm3	
	
	;# assume we have at least one i particle - start directly 		
.mci0310_outer:
	mov   eax, [ebp + mci0310_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci0310_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci0310_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci0310_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + mci0310_shX], mm0
	movd  [esp + mci0310_shZ], mm1

	mov   ecx, [ebp + mci0310_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci0310_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + mci0310_nsatoms]
	add dword ptr [ebp + mci0310_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + mci0310_nsvdwc], edx
	mov   [esp + mci0310_nscoul], eax
	mov   [esp + mci0310_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + mci0310_vnbtot], mm7
	mov   [esp + mci0310_solnr],  ebx
	
	mov   eax, [ebp + mci0310_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci0310_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + mci0310_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci0310_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + mci0310_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + mci0310_pos]
	
	mov   ecx, [esp + mci0310_nsvdwc]
	cmp   ecx,  0
	jnz   .mci0310_mno_vdwc
	jmp   .mci0310_testvdw
.mci0310_mno_vdwc:
	mov   ebx,  [esp + mci0310_solnr]
	inc   dword ptr [esp + mci0310_solnr]

	mov   edx, [ebp + mci0310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci0310_ntype]
	shl   edx, 1
	mov   [esp + mci0310_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci0310_pos]    ;# eax = base of pos[] 
	mov   [esp + mci0310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci0310_shX]
	pfadd mm1, [esp + mci0310_shZ]
	movq  [esp + mci0310_ix], mm0	
	movd  [esp + mci0310_iz], mm1	

	mov   ecx, [esp + mci0310_innerjjnr0]
	mov   [esp + mci0310_innerjjnr], ecx
	mov   edx, [esp + mci0310_innerk0]
    sub   edx,  2
    mov   [esp + mci0310_innerk], edx    ;# number of innerloop atoms 
	jge   .mci0310_unroll_vdwc_loop
	jmp   .mci0310_finish_vdwc_inner
.mci0310_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci0310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci0310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci0310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci0310_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci0310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci0310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci0310_c6], mm5
	movq [esp + mci0310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci0310_pos]

	movq  mm0, [esp + mci0310_ix]
	movd  mm1, [esp + mci0310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 
	;# do potential and fscal 
	pfmul mm1, [esp + mci0310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci0310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci0310_VFtab]
	;# dispersion table 
	mov ecx, [esp + mci0310_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci0310_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci0310_c6]
	pfmul mm5, mm4	;# vnb6            
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci0310_vnbtot]      ;# add the earlier value 
	movq [esp + mci0310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + mci0310_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + mci0310_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci0310_c12]
	pfmul mm5, mm6	;# vnb12 
	pfadd mm5, [esp + mci0310_vnbtot]      ;# add the earlier value 
	movq [esp + mci0310_vnbtot], mm5       ;# store the sum       
		
	;# should we do one more iteration? 
	sub dword ptr [esp + mci0310_innerk],  2
	jl    .mci0310_finish_vdwc_inner
	jmp   .mci0310_unroll_vdwc_loop
.mci0310_finish_vdwc_inner:	
	and dword ptr [esp + mci0310_innerk],  1
	jnz  .mci0310_single_vdwc_inner
	jmp  .mci0310_updateouterdata_vdwc		
.mci0310_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci0310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci0310_nbfp]
	mov ecx, [ebp + mci0310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci0310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci0310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci0310_c12], mm5

	mov   esi, [ebp + mci0310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci0310_ix]
	movd  mm1, [esp + mci0310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci0310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci0310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci0310_VFtab]
	mov ecx, [esp + mci0310_n1]
	shl ecx, 3
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci0310_c6]
	pfmul mm5, mm4	;# vnb6            
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci0310_vnbtot]      ;# add the earlier value 
	movq [esp + mci0310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci0310_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci0310_vnbtot]      ;# add the earlier value 
	movq [esp + mci0310_vnbtot], mm5       ;# store the sum       

.mci0310_updateouterdata_vdwc:	
	;# loop back to mno 
	dec dword ptr [esp + mci0310_nsvdwc]
	jz  .mci0310_testvdw
	jmp .mci0310_mno_vdwc
.mci0310_testvdw:	
	mov  ebx,  [esp + mci0310_nscoul]
	add  [esp + mci0310_solnr],  ebx

	mov  ecx, [esp + mci0310_nsvdw]
	cmp  ecx,  0
	jnz  .mci0310_mno_vdw
	jmp  .mci0310_last_mno
.mci0310_mno_vdw:
	mov   ebx,  [esp + mci0310_solnr]
	inc   dword ptr [esp + mci0310_solnr]

	mov   edx, [ebp + mci0310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci0310_ntype]
	shl   edx, 1
	mov   [esp + mci0310_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci0310_pos]    ;# eax = base of pos[] 
	mov   [esp + mci0310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci0310_shX]
	pfadd mm1, [esp + mci0310_shZ]
	movq  [esp + mci0310_ix], mm0	
	movd  [esp + mci0310_iz], mm1	

	mov   ecx, [esp + mci0310_innerjjnr0]
	mov   [esp + mci0310_innerjjnr], ecx
	mov   edx, [esp + mci0310_innerk0]
    sub   edx,  2
    mov   [esp + mci0310_innerk], edx    ;# number of innerloop atoms 
	jge   .mci0310_unroll_vdw_loop
	jmp   .mci0310_finish_vdw_inner
.mci0310_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci0310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci0310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci0310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci0310_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci0310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci0310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci0310_c6], mm5
	movq [esp + mci0310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci0310_pos]

	movq  mm0, [esp + mci0310_ix]
	movd  mm1, [esp + mci0310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 
	;# do potential and fscal 
	pfmul mm1, [esp + mci0310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci0310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci0310_VFtab]
	;# dispersion table 
	mov ecx, [esp + mci0310_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci0310_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci0310_c6]
	pfmul mm5, mm4	;# vnb6            
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci0310_vnbtot]      ;# add the earlier value 
	movq [esp + mci0310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + mci0310_n1]
	shl ecx, 3
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + mci0310_n1 + 4]
	shl ecx, 3
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci0310_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci0310_vnbtot]      ;# add the earlier value 
	movq [esp + mci0310_vnbtot], mm5       ;# store the sum 
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci0310_innerk],  2
	jl    .mci0310_finish_vdw_inner
	jmp   .mci0310_unroll_vdw_loop
.mci0310_finish_vdw_inner:	
	and dword ptr [esp + mci0310_innerk],  1
	jnz  .mci0310_single_vdw_inner
	jmp  .mci0310_updateouterdata_vdw		
.mci0310_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci0310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci0310_nbfp]
	mov ecx, [ebp + mci0310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci0310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci0310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci0310_c12], mm5

	mov   esi, [ebp + mci0310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci0310_ix]
	movd  mm1, [esp + mci0310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci0310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci0310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci0310_VFtab]
	mov ecx, [esp + mci0310_n1]
	shl ecx, 3
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci0310_c6]
	pfmul mm5, mm4	;# vnb6  
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci0310_vnbtot]      ;# add the earlier value 
	movq [esp + mci0310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci0310_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci0310_vnbtot]      ;# add the earlier value 
	movq [esp + mci0310_vnbtot], mm5       ;# store the sum       

.mci0310_updateouterdata_vdw:	
	;# loop back to mno 
	dec dword ptr [esp + mci0310_nsvdw]
	jz  .mci0310_last_mno
	jmp .mci0310_mno_vdw
	
.mci0310_last_mno:	
	mov   edx, [ebp + mci0310_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci0310_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci0310_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci0310_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + mci0310_nri]
	dec ecx
	jecxz .mci0310_end
	;# not last, iterate once more! 
	mov [ebp + mci0310_nri], ecx
	jmp .mci0310_outer
.mci0310_end:
	femms
	add esp, 108
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl mcinl1000_3dnow
.globl _mcinl1000_3dnow
mcinl1000_3dnow:	
_mcinl1000_3dnow:	
.equiv		mci1000_nri,		8
.equiv		mci1000_iinr,		12
.equiv		mci1000_jindex,		16
.equiv		mci1000_jjnr,		20
.equiv		mci1000_shift,		24
.equiv		mci1000_shiftvec,	28
.equiv		mci1000_gid,		32
.equiv		mci1000_pos,		36		
.equiv		mci1000_charge,		40
.equiv		mci1000_facel,		44
.equiv		mci1000_Vc,		48			
	;# stack offsets for local variables 
.equiv		mci1000_is3,		0
.equiv		mci1000_ii3,		4
.equiv		mci1000_ix,		8
.equiv		mci1000_iy,		12
.equiv		mci1000_iz,		16
.equiv		mci1000_iq,		20		
.equiv		mci1000_vctot,		28 
.equiv		mci1000_innerjjnr,	36
.equiv		mci1000_innerk,		40
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 44		;# 80 bytes local stack space 
	femms
	;# assume we have at least one i particle - start directly 	
.mci1000_outer:
	mov   eax, [ebp + mci1000_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci1000_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci1000_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci1000_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + mci1000_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci1000_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci1000_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci1000_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci1000_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1000_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci1000_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + mci1000_ix], mm0	
	movd  [esp + mci1000_iz], mm1	
				
	;# clear vctot 
	pxor  mm7,mm7
	movq  [esp + mci1000_vctot], mm7

	mov   eax, [ebp + mci1000_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci1000_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + mci1000_pos]
	mov   eax, [ebp + mci1000_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci1000_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + mci1000_innerk], edx    ;# number of innerloop atoms 
	jge   .mci1000_unroll_loop
	jmp   .mci1000_finish_inner
.mci1000_unroll_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci1000_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci1000_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci1000_charge]    ;# base of charge[] 
	movq mm5, [esp + mci1000_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
	movd mm7, [ecx + ebx*4]  	 ;# charge[jnr2] 
	punpckldq mm3,mm7	     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + mci1000_ix]
	movd  mm1, [esp + mci1000_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	;# mm1=invsqrt
	 ;# do potential and fscal
	 
	
	pfmul mm3,mm1		;# 3 has both vcoul 
	pfadd mm3, [esp + mci1000_vctot]      ;# add the earlier value  
	movq [esp + mci1000_vctot], mm3       ;# store the sum 
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci1000_innerk],  2
	jl    .mci1000_finish_inner
	jmp   .mci1000_unroll_loop
.mci1000_finish_inner:	
	and dword ptr [esp + mci1000_innerk],  1
	jnz  .mci1000_single_inner
	jmp  .mci1000_updateouterdata		
.mci1000_single_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci1000_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci1000_charge]
	movd mm6, [esp + mci1000_iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;# mm6=qq 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + mci1000_ix]
	movd  mm1, [esp + mci1000_iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;# mm0=rsq 
	
        pfrsqrt mm1,mm0
        movq mm2,mm1
	pfmul mm1,mm1
        pfrsqit1 mm1,mm0				
        pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfadd mm6, [esp + mci1000_vctot]
	movq [esp + mci1000_vctot], mm6
	
.mci1000_updateouterdata:	
	mov   edx, [ebp + mci1000_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci1000_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci1000_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci1000_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + mci1000_nri]
	dec ecx
	jecxz .mci1000_end
	;# not last, iterate once more! 
	mov [ebp + mci1000_nri], ecx
	jmp .mci1000_outer
.mci1000_end:
	femms
	add esp, 44
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl mcinl1010_3dnow
.globl _mcinl1010_3dnow
mcinl1010_3dnow:	
_mcinl1010_3dnow:	
.equiv		mci1010_nri,		8
.equiv		mci1010_iinr,		12
.equiv		mci1010_jindex,		16
.equiv		mci1010_jjnr,		20
.equiv		mci1010_shift,		24
.equiv		mci1010_shiftvec,	28
.equiv		mci1010_gid,		32
.equiv		mci1010_pos,		36		
.equiv		mci1010_charge,		40
.equiv		mci1010_facel,		44
.equiv		mci1010_Vc,		48
.equiv		mci1010_nsatoms,	52		
	;# stack offsets for local variables 
.equiv		mci1010_is3,		0
.equiv		mci1010_ii3,		4
.equiv		mci1010_shX,		8
.equiv		mci1010_shY,		12 
.equiv		mci1010_shZ,		16	
.equiv		mci1010_ix,		20
.equiv		mci1010_iy,		24
.equiv		mci1010_iz,		28	
.equiv		mci1010_iq,		32 
.equiv		mci1010_vctot,		40 
.equiv		mci1010_innerjjnr0,	48
.equiv		mci1010_innerk0,	52		
.equiv		mci1010_innerjjnr,	56
.equiv		mci1010_innerk,		60
.equiv		mci1010_nscoul,		64
.equiv		mci1010_solnr,		68		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 72		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	
	add dword ptr [ebp + mci1010_nsatoms],  8

.mci1010_outer:
	mov   eax, [ebp + mci1010_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci1010_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci1010_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci1010_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + mci1010_shX], mm0
	movd  [esp + mci1010_shZ], mm1

	mov   ecx, [ebp + mci1010_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci1010_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + mci1010_nsatoms]
	mov   ecx, [eax]
	add dword ptr [ebp + mci1010_nsatoms],  12
	mov   [esp + mci1010_nscoul], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + mci1010_vctot], mm7
	mov   [esp + mci1010_solnr], ebx
	
	mov   eax, [ebp + mci1010_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci1010_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + mci1010_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci1010_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + mci1010_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + mci1010_pos]

	mov   ecx, [esp + mci1010_nscoul]
	cmp   ecx,  0
	jnz   .mci1010_mno_coul
	jmp   .mci1010_last_mno
.mci1010_mno_coul:				
	mov   ebx,  [esp + mci1010_solnr]
	inc   dword ptr [esp + mci1010_solnr]
	mov   edx, [ebp + mci1010_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci1010_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci1010_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1010_pos]    ;# eax = base of pos[] 
	mov   [esp + mci1010_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci1010_shX]
	pfadd mm1, [esp + mci1010_shZ]
	movq  [esp + mci1010_ix], mm0	
	movd  [esp + mci1010_iz], mm1	

	mov   ecx, [esp + mci1010_innerjjnr0]
	mov   [esp + mci1010_innerjjnr], ecx
	mov   edx, [esp + mci1010_innerk0]
    sub   edx,  2
    mov   [esp + mci1010_innerk], edx    ;# number of innerloop atoms 
	jge   .mci1010_unroll_coul_loop
	jmp   .mci1010_finish_coul_inner
.mci1010_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci1010_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci1010_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci1010_charge]    ;# base of charge[] 
	movq mm5, [esp + mci1010_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
	movd mm7, [ecx + ebx*4]  	 ;# charge[jnr2] 
	punpckldq mm3,mm7	     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + mci1010_ix]
	movd  mm1, [esp + mci1010_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	;# mm1=invsqrt 
	;# do potential 
	
	pfmul mm3,mm1			;# 3 has both vcoul 
	pfadd mm3, [esp + mci1010_vctot]      ;# add the earlier value  
	movq [esp + mci1010_vctot], mm3       ;# store the sum 
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci1010_innerk],  2
	jl    .mci1010_finish_coul_inner
	jmp   .mci1010_unroll_coul_loop
.mci1010_finish_coul_inner:	
	and dword ptr [esp + mci1010_innerk],  1
	jnz  .mci1010_single_coul_inner
	jmp  .mci1010_updateouterdata_coul		
.mci1010_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci1010_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci1010_charge]
	movd mm6, [esp + mci1010_iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;# mm6=qq 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + mci1010_ix]
	movd  mm1, [esp + mci1010_iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;# mm0=rsq 
	
    pfrsqrt mm1,mm0
    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	pfadd mm6, [esp + mci1010_vctot] 
	movq [esp + mci1010_vctot], mm6
	
.mci1010_updateouterdata_coul:	
	;# loop back to mno 
	dec dword ptr [esp + mci1010_nscoul]
	jz  .mci1010_last_mno
	jmp .mci1010_mno_coul
.mci1010_last_mno:	
	mov   edx, [ebp + mci1010_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci1010_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci1010_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci1010_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + mci1010_nri]
	dec ecx
	jecxz .mci1010_end
	;# not last, iterate once more! 
	mov [ebp + mci1010_nri], ecx
	jmp .mci1010_outer
.mci1010_end:
	femms
	add esp, 72
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

			
.globl mcinl1020_3dnow
.globl _mcinl1020_3dnow
mcinl1020_3dnow:	
_mcinl1020_3dnow:	
.equiv		mci1020_nri,		8
.equiv		mci1020_iinr,		12
.equiv		mci1020_jindex,		16
.equiv		mci1020_jjnr,		20
.equiv		mci1020_shift,		24
.equiv		mci1020_shiftvec,	28
.equiv		mci1020_gid,		32
.equiv		mci1020_pos,		36		
.equiv		mci1020_charge,		40
.equiv		mci1020_facel,		44
.equiv		mci1020_Vc,			48			
			;# stack offsets for local variables 
.equiv		mci1020_is3,		0
.equiv		mci1020_ii3,		4
.equiv		mci1020_ixO,		8
.equiv		mci1020_iyO,		12
.equiv		mci1020_izO,		16	 
.equiv		mci1020_ixH,		20 
.equiv		mci1020_iyH,		28 
.equiv		mci1020_izH,		36 
.equiv		mci1020_iqO,		44	
.equiv		mci1020_iqH,		52		
.equiv		mci1020_vctot,		60 
.equiv		mci1020_innerjjnr,	68
.equiv		mci1020_innerk,		72      
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 76		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + mci1020_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci1020_charge]
	movd  mm1, [ebp + mci1020_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, mm1		
	movq  [esp + mci1020_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci1020_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 
.mci1020_outer:
	mov   eax, [ebp + mci1020_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci1020_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci1020_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci1020_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci1020_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci1020_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1020_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci1020_ii3], ebx	    ;# (use mm7 as temp storage for iz) 
	pfadd mm6, mm7
	movq  [esp + mci1020_ixO], mm5	
	movq  [esp + mci1020_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci1020_ixH], mm0	
	movq [esp + mci1020_iyH], mm1	
	movq [esp + mci1020_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci1020_vctot], mm7

	mov   eax, [ebp + mci1020_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci1020_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci1020_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + mci1020_pos]
	mov   eax, [ebp + mci1020_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci1020_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci1020_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci1020_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci1020_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + mci1020_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + mci1020_iqO]
	pfmul mm7, [esp + mci1020_iqH]	;# mm6=qqO, mm7=qqH 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci1020_ixO]
	pfsubr mm1, [esp + mci1020_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci1020_ixH]
	pfsubr mm3, [esp + mci1020_iyH]
	pfsubr mm4, [esp + mci1020_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    pfrsqit1 mm5,mm3				
    pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + mci1020_vctot]
	movq [esp + mci1020_vctot], mm7
	
	;#  done  - one more? 
	dec dword ptr [esp + mci1020_innerk]
	jz  .mci1020_updateouterdata
	jmp .mci1020_inner_loop
.mci1020_updateouterdata:	
	mov   edx, [ebp + mci1020_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci1020_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci1020_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci1020_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	;# finish if last 
	dec dword ptr [ebp + mci1020_nri]
	jz  .mci1020_end
	;# not last, iterate once more! 
	jmp .mci1020_outer
.mci1020_end:
	femms
	add esp, 76
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl mcinl1030_3dnow
.globl _mcinl1030_3dnow
mcinl1030_3dnow:	
_mcinl1030_3dnow:	
.equiv		mci1030_nri,		8
.equiv		mci1030_iinr,		12
.equiv		mci1030_jindex,		16
.equiv		mci1030_jjnr,		20
.equiv		mci1030_shift,		24
.equiv		mci1030_shiftvec,	28
.equiv		mci1030_gid,		32
.equiv		mci1030_pos,		36		
.equiv		mci1030_charge,		40
.equiv		mci1030_facel,		44
.equiv		mci1030_Vc,			48
			;# stack offsets for local variables 
.equiv		mci1030_is3,		0
.equiv		mci1030_ii3,		4
.equiv		mci1030_ixO,		8
.equiv		mci1030_iyO,		12
.equiv		mci1030_izO,		16	
.equiv		mci1030_ixH,		20
.equiv		mci1030_iyH,		28
.equiv		mci1030_izH,		36
.equiv		mci1030_qqOO,		44		
.equiv		mci1030_qqOH,		52		
.equiv		mci1030_qqHH,		60     	
.equiv		mci1030_vctot,		68
.equiv		mci1030_innerjjnr,	76
.equiv		mci1030_innerk,		80   
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 84		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + mci1030_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci1030_charge]
	movd  mm1, [ebp + mci1030_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + mci1030_qqOO], mm4
	movq  [esp + mci1030_qqOH], mm5
	movq  [esp + mci1030_qqHH], mm6
.mci1030_outer:
	mov   eax, [ebp + mci1030_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci1030_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci1030_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci1030_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci1030_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci1030_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1030_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci1030_ii3], ebx	    ;# (use mm7 as temp storage for iz) 
	pfadd mm6, mm7
	movq  [esp + mci1030_ixO], mm5	
	movq  [esp + mci1030_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci1030_ixH], mm0	
	movq [esp + mci1030_iyH], mm1	
	movq [esp + mci1030_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci1030_vctot], mm7
	
	mov   eax, [ebp + mci1030_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci1030_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci1030_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + mci1030_pos]
	mov   eax, [ebp + mci1030_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci1030_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci1030_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci1030_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci1030_innerjjnr],  4 ;# advance pointer 

	movd  mm6, [esp + mci1030_qqOO]
	movq  mm7, [esp + mci1030_qqOH]

	lea   eax, [eax + eax*2]
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci1030_ixO]
	pfsubr mm1, [esp + mci1030_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci1030_ixH]
	pfsubr mm3, [esp + mci1030_iyH]
	pfsubr mm4, [esp + mci1030_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    pfrsqit1 mm5,mm3				
    pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + mci1030_vctot]
	movq [esp + mci1030_vctot], mm7
	
	;# interactions with j H1 
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	movd mm6, [esp + mci1030_qqOH]
	movq mm7, [esp + mci1030_qqHH]
	
	pfsubr mm0, [esp + mci1030_ixO]
	pfsubr mm1, [esp + mci1030_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci1030_ixH]
	pfsubr mm3, [esp + mci1030_iyH]
	pfsubr mm4, [esp + mci1030_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    pfrsqit1 mm5,mm3				
    pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + mci1030_vctot]
	movq [esp + mci1030_vctot], mm7
	
	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	movd mm6, [esp + mci1030_qqOH]
	movq mm7, [esp + mci1030_qqHH]

	pfsubr mm0, [esp + mci1030_ixO]
	pfsubr mm1, [esp + mci1030_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci1030_ixH]
	pfsubr mm3, [esp + mci1030_iyH]
	pfsubr mm4, [esp + mci1030_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    pfrsqit1 mm5,mm3				
    pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + mci1030_vctot]
	movq [esp + mci1030_vctot], mm7
	
	;#  done  - one more? 
	dec dword ptr [esp + mci1030_innerk]
	jz  .mci1030_updateouterdata
	jmp .mci1030_inner_loop	
.mci1030_updateouterdata:	
	mov   edx, [ebp + mci1030_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci1030_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci1030_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + mci1030_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	dec dword ptr [ebp + mci1030_nri]
	jz  .mci1030_end
	;# not last, iterate once more! 
	jmp .mci1030_outer
.mci1030_end:
	femms
	add esp, 84
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl mcinl1100_3dnow
.globl _mcinl1100_3dnow
mcinl1100_3dnow:	
_mcinl1100_3dnow:	
.equiv		mci1100_nri,		8
.equiv		mci1100_iinr,		12
.equiv		mci1100_jindex,		16
.equiv		mci1100_jjnr,		20
.equiv		mci1100_shift,		24
.equiv		mci1100_shiftvec,	28
.equiv		mci1100_gid,		32
.equiv		mci1100_pos,		36		
.equiv		mci1100_charge,		40
.equiv		mci1100_facel,		44
.equiv		mci1100_Vc,			48			
.equiv		mci1100_type,		52
.equiv		mci1100_ntype,		56
.equiv		mci1100_nbfp,		60	
.equiv		mci1100_Vnb,		64	
	;# stack offsets for local variables 
.equiv		mci1100_is3,		0
.equiv		mci1100_ii3,		4
.equiv		mci1100_ix,			8
.equiv		mci1100_iy,			12
.equiv		mci1100_iz,			16
.equiv		mci1100_iq,			20 
.equiv		mci1100_vctot,		28 
.equiv		mci1100_vnbtot,		36 
.equiv		mci1100_c6,			44 
.equiv		mci1100_c12,		52
.equiv		mci1100_ntia,		60
.equiv		mci1100_innerjjnr,	64
.equiv		mci1100_innerk,		68					
	push ebp
	mov ebp,esp
	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 72		;# local stack space 
	femms
	;# move data to local stack  
	;# assume we have at least one i particle - start directly 	
.mci1100_outer:
	mov   eax, [ebp + mci1100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci1100_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci1100_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci1100_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + mci1100_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci1100_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci1100_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci1100_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci1100_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + mci1100_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci1100_ntype]
	shl   edx, 1
	mov   [esp + mci1100_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1100_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci1100_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + mci1100_ix], mm0	
	movd  [esp + mci1100_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + mci1100_vctot],  mm7
	movq  [esp + mci1100_vnbtot], mm7

	mov   eax, [ebp + mci1100_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci1100_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + mci1100_pos]
	mov   eax, [ebp + mci1100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci1100_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + mci1100_innerk], edx    ;# number of innerloop atoms 
	jge   .mci1100_unroll_loop
	jmp   .mci1100_finish_inner
.mci1100_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + mci1100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci1100_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci1100_charge]    ;# base of charge[] 
	movq mm5, [esp + mci1100_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + mci1100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci1100_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci1100_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci1100_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci1100_c6], mm5
	movq [esp + mci1100_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci1100_pos]

	movq  mm0, [esp + mci1100_ix]
	movd  mm1, [esp + mci1100_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	pfmul mm5, [esp + mci1100_c12]
	pfmul mm4, [esp + mci1100_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vctot 
	pfadd mm3, [esp + mci1100_vctot]      ;# add the earlier value 
	movq [esp + mci1100_vctot], mm3       ;# store the sum       
	;# update vnbtot 
	pfadd mm6, [esp + mci1100_vnbtot]      ;# add the earlier value 
	movq [esp + mci1100_vnbtot], mm6       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci1100_innerk],  2
	jl    .mci1100_finish_inner
	jmp   .mci1100_unroll_loop
.mci1100_finish_inner:	
	and dword ptr [esp + mci1100_innerk],  1
	jnz  .mci1100_single_inner
	jmp  .mci1100_updateouterdata		
.mci1100_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci1100_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci1100_charge]
	movd mm5, [esp + mci1100_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + mci1100_nbfp]
	mov ecx, [ebp + mci1100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci1100_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci1100_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci1100_c12], mm5


	mov   esi, [ebp + mci1100_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci1100_ix]
	movd  mm1, [esp + mci1100_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	movq  mm1, mm0
	pfmul mm0, mm0		;# mm0=invsq 
	;# calculate potentials and scalar force 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	pfmul mm5, [esp + mci1100_c12]
	pfmul mm4, [esp + mci1100_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vctot 
	pfadd mm3, [esp + mci1100_vctot]
	movq [esp + mci1100_vctot], mm3
	;# update vnbtot 
	pfadd mm6, [esp + mci1100_vnbtot]      ;# add the earlier value 
	movq [esp + mci1100_vnbtot], mm6       ;# store the sum       

.mci1100_updateouterdata:	
	mov   edx, [ebp + mci1100_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci1100_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci1100_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci1100_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci1100_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci1100_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + mci1100_nri]
	dec ecx
	jecxz .mci1100_end
	;# not last, iterate once more! 
	mov [ebp + mci1100_nri], ecx
	jmp .mci1100_outer
.mci1100_end:
	femms
	add esp, 72
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	



.globl mcinl1110_3dnow
.globl _mcinl1110_3dnow
mcinl1110_3dnow:	
_mcinl1110_3dnow:	
.equiv		mci1110_nri,		8
.equiv		mci1110_iinr,		12
.equiv		mci1110_jindex,		16
.equiv		mci1110_jjnr,		20
.equiv		mci1110_shift,		24
.equiv		mci1110_shiftvec,	28
.equiv		mci1110_gid,		32
.equiv		mci1110_pos,		36		
.equiv		mci1110_charge,		40
.equiv		mci1110_facel,		44
.equiv		mci1110_Vc,			48	
.equiv		mci1110_type,		52
.equiv		mci1110_ntype,		56
.equiv		mci1110_nbfp,		60	
.equiv		mci1110_Vnb,		64				
.equiv		mci1110_nsatoms,	68		
	;# stack offsets for local variables 
.equiv		mci1110_is3,		0
.equiv		mci1110_ii3,		4
.equiv		mci1110_shX,		8
.equiv		mci1110_shY,		12 
.equiv		mci1110_shZ,		16	
.equiv		mci1110_ix,			20
.equiv		mci1110_iy,			24
.equiv		mci1110_iz,			28	
.equiv		mci1110_iq,			32		 
.equiv		mci1110_vctot,		40 
.equiv		mci1110_vnbtot,		48 
.equiv		mci1110_c6,			56 
.equiv		mci1110_c12,		64
.equiv		mci1110_ntia,		72	
.equiv		mci1110_innerjjnr0,	76
.equiv		mci1110_innerk0,	80		
.equiv		mci1110_innerjjnr,	84
.equiv		mci1110_innerk,		88				
.equiv		mci1110_nsvdwc,		92
.equiv		mci1110_nscoul,		96
.equiv		mci1110_nsvdw,		100
.equiv		mci1110_solnr,		104		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 108		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 		
.mci1110_outer:
	mov   eax, [ebp + mci1110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci1110_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci1110_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci1110_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + mci1110_shX], mm0
	movd  [esp + mci1110_shZ], mm1

	mov   ecx, [ebp + mci1110_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci1110_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + mci1110_nsatoms]
	add dword ptr [ebp + mci1110_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + mci1110_nsvdwc], edx
	mov   [esp + mci1110_nscoul], eax
	mov   [esp + mci1110_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + mci1110_vctot],  mm7
	movq  [esp + mci1110_vnbtot], mm7
	mov   [esp + mci1110_solnr],  ebx
	
	mov   eax, [ebp + mci1110_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci1110_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + mci1110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci1110_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + mci1110_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + mci1110_pos]
	
	mov   ecx, [esp + mci1110_nsvdwc]
	cmp   ecx,  0
	jnz   .mci1110_mno_vdwc
	jmp   .mci1110_testcoul
.mci1110_mno_vdwc:
	mov   ebx,  [esp + mci1110_solnr]
	inc   dword ptr [esp + mci1110_solnr]
	mov   edx, [ebp + mci1110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci1110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci1110_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + mci1110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci1110_ntype]
	shl   edx, 1
	mov   [esp + mci1110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1110_pos]    ;# eax = base of pos[] 
	mov   [esp + mci1110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci1110_shX]
	pfadd mm1, [esp + mci1110_shZ]
	movq  [esp + mci1110_ix], mm0	
	movd  [esp + mci1110_iz], mm1	

	mov   ecx, [esp + mci1110_innerjjnr0]
	mov   [esp + mci1110_innerjjnr], ecx
	mov   edx, [esp + mci1110_innerk0]
    sub   edx,  2
    mov   [esp + mci1110_innerk], edx    ;# number of innerloop atoms 
	jge   .mci1110_unroll_vdwc_loop
	jmp   .mci1110_finish_vdwc_inner
.mci1110_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci1110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci1110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci1110_charge]    ;# base of charge[] 
	movq mm5, [esp + mci1110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + mci1110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci1110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci1110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci1110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci1110_c6], mm5
	movq [esp + mci1110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci1110_pos]

	movq  mm0, [esp + mci1110_ix]
	movd  mm1, [esp + mci1110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 
	pfmul mm5, [esp + mci1110_c12]
	pfmul mm4, [esp + mci1110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vctot 
	pfadd mm3, [esp + mci1110_vctot]      ;# add the earlier value 
	movq [esp + mci1110_vctot], mm3       ;# store the sum       
	;# update vnbtot 
	pfadd mm6, [esp + mci1110_vnbtot]      ;# add the earlier value 
	movq [esp + mci1110_vnbtot], mm6       ;# store the sum       

	;# should we do one more iteration? 
	sub dword ptr [esp + mci1110_innerk],  2
	jl    .mci1110_finish_vdwc_inner
	jmp   .mci1110_unroll_vdwc_loop
.mci1110_finish_vdwc_inner:	
	and dword ptr [esp + mci1110_innerk],  1
	jnz  .mci1110_single_vdwc_inner
	jmp  .mci1110_updateouterdata_vdwc		
.mci1110_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci1110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci1110_charge]
	movd mm5, [esp + mci1110_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + mci1110_nbfp]
	mov ecx, [ebp + mci1110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci1110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci1110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci1110_c12], mm5


	mov   esi, [ebp + mci1110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci1110_ix]
	movd  mm1, [esp + mci1110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	movq  mm1, mm0
	pfmul mm0, mm0		;# mm0=invsq 
	;# calculate potentials and scalar force 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm3, mm1		;# mm3 has vcoul for both interactions 

	pfmul mm5, [esp + mci1110_c12]
	pfmul mm4, [esp + mci1110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vctot 
	pfadd mm3, [esp + mci1110_vctot]
	movq [esp + mci1110_vctot], mm3
	;# update vnbtot 
	pfadd mm6, [esp + mci1110_vnbtot]      ;# add the earlier value 
	movq [esp + mci1110_vnbtot], mm6       ;# store the sum       
.mci1110_updateouterdata_vdwc:	
	;# loop back to mno 
	dec  dword ptr [esp + mci1110_nsvdwc]
	jz  .mci1110_testcoul
	jmp .mci1110_mno_vdwc
.mci1110_testcoul:	
	mov  ecx, [esp + mci1110_nscoul]
	cmp  ecx,  0
	jnz  .mci1110_mno_coul
	jmp  .mci1110_testvdw
.mci1110_mno_coul:
	mov   ebx,  [esp + mci1110_solnr]
	inc   dword ptr [esp + mci1110_solnr]
	mov   edx, [ebp + mci1110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci1110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci1110_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1110_pos]    ;# eax = base of pos[] 
	mov   [esp + mci1110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci1110_shX]
	pfadd mm1, [esp + mci1110_shZ]
	movq  [esp + mci1110_ix], mm0	
	movd  [esp + mci1110_iz], mm1	

	mov   ecx, [esp + mci1110_innerjjnr0]
	mov   [esp + mci1110_innerjjnr], ecx
	mov   edx, [esp + mci1110_innerk0]
    sub   edx,  2
    mov   [esp + mci1110_innerk], edx    ;# number of innerloop atoms 
	jge   .mci1110_unroll_coul_loop
	jmp   .mci1110_finish_coul_inner
.mci1110_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci1110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci1110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci1110_charge]    ;# base of charge[] 
	movq mm5, [esp + mci1110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
	movd mm7, [ecx + ebx*4]  	 ;# charge[jnr2] 
	punpckldq mm3,mm7	     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	movq  mm0, [esp + mci1110_ix]
	movd  mm1, [esp + mci1110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	;# mm1 is invsqrt 
	;# do potential and fscal 
	pfmul mm3,mm1		;# 3 has both vcoul 
	pfadd mm3, [esp + mci1110_vctot]      ;# add the earlier value  
	movq [esp + mci1110_vctot], mm3       ;# store the sum 
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci1110_innerk],  2
	jl    .mci1110_finish_coul_inner
	jmp   .mci1110_unroll_coul_loop
.mci1110_finish_coul_inner:	
	and dword ptr [esp + mci1110_innerk],  1
	jnz  .mci1110_single_coul_inner
	jmp  .mci1110_updateouterdata_coul		
.mci1110_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci1110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci1110_charge]
	movd mm6, [esp + mci1110_iq]
	movd mm7, [ecx + eax*4]
	pfmul mm6, mm7	  	;# mm6=qq 
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esp + mci1110_ix]
	movd  mm1, [esp + mci1110_iz]
	movq  mm2, [esi + eax*4]
	movd  mm3, [esi + eax*4 + 8]
	pfsub mm0, mm2
	pfsub mm1, mm3
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfacc mm0, mm1		;# mm0=rsq 
	
    pfrsqrt mm1,mm0
    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	;# update vctot 
	pfadd mm6, [esp + mci1110_vctot]
	movq [esp + mci1110_vctot], mm6
	
.mci1110_updateouterdata_coul:	
	;# loop back to mno 
	dec dword ptr [esp + mci1110_nscoul]
	jz  .mci1110_testvdw
	jmp .mci1110_mno_coul
.mci1110_testvdw:	
	mov  ecx, [esp + mci1110_nsvdw]
	cmp  ecx,  0
	jnz  .mci1110_mno_vdw
	jmp  .mci1110_last_mno
.mci1110_mno_vdw:
	mov   ebx,  [esp + mci1110_solnr]
	inc   dword ptr [esp + mci1110_solnr]

	mov   edx, [ebp + mci1110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci1110_ntype]
	shl   edx, 1
	mov   [esp + mci1110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1110_pos]    ;# eax = base of pos[] 
	mov   [esp + mci1110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci1110_shX]
	pfadd mm1, [esp + mci1110_shZ]
	movq  [esp + mci1110_ix], mm0	
	movd  [esp + mci1110_iz], mm1

	mov   ecx, [esp + mci1110_innerjjnr0]
	mov   [esp + mci1110_innerjjnr], ecx
	mov   edx, [esp + mci1110_innerk0]
    sub   edx,  2
    mov   [esp + mci1110_innerk], edx    ;# number of innerloop atoms 
	jge   .mci1110_unroll_vdw_loop
	jmp   .mci1110_finish_vdw_inner
.mci1110_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci1110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci1110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + mci1110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci1110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci1110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci1110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci1110_c6], mm5
	movq [esp + mci1110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci1110_pos]

	movq  mm0, [esp + mci1110_ix]
	movd  mm1, [esp + mci1110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	movq mm1,mm0
	pfmul mm0,mm0
	;# mm0 now contains invsq, and mm1 invsqrt 
	;# do potential and fscal 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci1110_c12]
	pfmul mm4, [esp + mci1110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci1110_vnbtot]      ;# add the earlier value 
	movq [esp + mci1110_vnbtot], mm6       ;# store the sum       

	;# should we do one more iteration? 
	sub dword ptr [esp + mci1110_innerk],  2
	jl    .mci1110_finish_vdw_inner
	jmp   .mci1110_unroll_vdw_loop
.mci1110_finish_vdw_inner:	
	and dword ptr [esp + mci1110_innerk],  1
	jnz  .mci1110_single_vdw_inner
	jmp  .mci1110_updateouterdata_vdw		
.mci1110_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments 
	mov   eax, [esp + mci1110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci1110_nbfp]
	mov ecx, [ebp + mci1110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci1110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci1110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci1110_c12], mm5


	mov   esi, [ebp + mci1110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci1110_ix]
	movd  mm1, [esp + mci1110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	movq  mm1, mm0
	pfmul mm0, mm0		;# mm0=invsq 
	;# calculate potentials and scalar force 
	movq mm4, mm0
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci1110_c12]
	pfmul mm4, [esp + mci1110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci1110_vnbtot]      ;# add the earlier value 
	movq [esp + mci1110_vnbtot], mm6       ;# store the sum       

.mci1110_updateouterdata_vdw:	
	;# loop back to mno 
	dec dword ptr [esp + mci1110_nsvdw]
	jz  .mci1110_last_mno
	jmp .mci1110_mno_vdw
	
.mci1110_last_mno:	
	mov   edx, [ebp + mci1110_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci1110_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci1110_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci1110_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci1110_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci1110_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + mci1110_nri]
	dec ecx
	jecxz .mci1110_end
	;# not last, iterate once more! 
	mov [ebp + mci1110_nri], ecx
	jmp .mci1110_outer
.mci1110_end:
	femms
	add esp, 108
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret



.globl mcinl1120_3dnow
.globl _mcinl1120_3dnow
mcinl1120_3dnow:	
_mcinl1120_3dnow:	
.equiv		mci1120_nri,		8
.equiv		mci1120_iinr,		12
.equiv		mci1120_jindex,		16
.equiv		mci1120_jjnr,		20
.equiv		mci1120_shift,		24
.equiv		mci1120_shiftvec,	28
.equiv		mci1120_gid,		32
.equiv		mci1120_pos,		36		
.equiv		mci1120_charge,		40
.equiv		mci1120_facel,		44
.equiv		mci1120_Vc,		48	
.equiv		mci1120_type,		52
.equiv		mci1120_ntype,		56
.equiv		mci1120_nbfp,		60	
.equiv		mci1120_Vnb,		64				
			;# stack offsets for local variables 
.equiv		mci1120_is3,		0
.equiv		mci1120_ii3,		4
.equiv		mci1120_ixO,		8
.equiv		mci1120_iyO,		12
.equiv		mci1120_izO,		16	
.equiv		mci1120_ixH,		20  
.equiv		mci1120_iyH,		28  
.equiv		mci1120_izH,		36  
.equiv		mci1120_iqO,		44  
.equiv		mci1120_iqH,		52  
.equiv		mci1120_vctot,		60  
.equiv		mci1120_vnbtot,		68  
.equiv		mci1120_c6,		76  
.equiv		mci1120_c12,		84  
.equiv		mci1120_ntia,		92
.equiv		mci1120_innerjjnr,	96
.equiv		mci1120_innerk,		100
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 104		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + mci1120_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci1120_charge]
	movd  mm1, [ebp + mci1120_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + mci1120_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci1120_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + mci1120_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	imul  ecx, [ebp + mci1120_ntype]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + mci1120_ntia], ecx
	
.mci1120_outer:
	mov   eax, [ebp + mci1120_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci1120_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci1120_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci1120_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci1120_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci1120_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1120_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci1120_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + mci1120_ixO], mm5	
	movq  [esp + mci1120_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci1120_ixH], mm0	
	movq [esp + mci1120_iyH], mm1	
	movq [esp + mci1120_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci1120_vctot], mm7
	movq  [esp + mci1120_vnbtot], mm7

	mov   eax, [ebp + mci1120_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci1120_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci1120_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + mci1120_pos]
	mov   eax, [ebp + mci1120_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci1120_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci1120_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci1120_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci1120_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + mci1120_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + mci1120_iqO]
	pfmul mm7, [esp + mci1120_iqH]	;# mm6=qqO, mm7=qqH 

	mov ecx, [ebp + mci1120_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + mci1120_nbfp]
	shl edx, 1
	add edx, [esp + mci1120_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci1120_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci1120_c12], mm5	
	
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci1120_ixO]
	pfsubr mm1, [esp + mci1120_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci1120_ixH]
	pfsubr mm3, [esp + mci1120_iyH]
	pfsubr mm4, [esp + mci1120_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq 

	movq  mm0, mm4
	pfmul mm0, mm4
	pfmul mm0, mm4		;# mm0=rinvsix 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm2=rintwelve 
	
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	movq  mm1, mm6		;# use mm1 for fscal sum 

	;# LJ for the oxygen 
	pfmul mm0, [esp + mci1120_c6]	 
	pfmul mm2, [esp + mci1120_c12]	 

	;# calc nb potential 
	pfsub mm2, mm0
	;# update nb potential 
	pfadd mm2, [esp + mci1120_vnbtot]
	movq [esp + mci1120_vnbtot], mm2
	
	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3. 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + mci1120_vctot]
	movq [esp + mci1120_vctot], mm7
		
	;#  done  - one more? 
	dec dword ptr [esp + mci1120_innerk]
	jz  .mci1120_updateouterdata
	jmp .mci1120_inner_loop
.mci1120_updateouterdata:	
	mov   edx, [ebp + mci1120_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci1120_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci1120_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci1120_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci1120_vnbtot]     
	pfacc mm7,mm7	          ;# same for Vnb 
	
	mov   eax, [ebp + mci1120_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 
	;# finish if last 
	dec dword ptr [ebp + mci1120_nri]
	jz  .mci1120_end
	;# not last, iterate once more! 
	jmp .mci1120_outer
.mci1120_end:
	femms
	add esp, 104
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	

.globl mcinl1130_3dnow
.globl _mcinl1130_3dnow
mcinl1130_3dnow:	
_mcinl1130_3dnow:	
.equiv		mci1130_nri,		8
.equiv		mci1130_iinr,		12
.equiv		mci1130_jindex,		16
.equiv		mci1130_jjnr,		20
.equiv		mci1130_shift,		24
.equiv		mci1130_shiftvec,	28
.equiv		mci1130_gid,		32
.equiv		mci1130_pos,		36		
.equiv		mci1130_charge,		40
.equiv		mci1130_facel,		44
.equiv		mci1130_Vc,		48					
.equiv		mci1130_type,		52
.equiv		mci1130_ntype,		56
.equiv		mci1130_nbfp,		60	
.equiv		mci1130_Vnb,		64
			;# stack offsets for local variables 
.equiv		mci1130_is3,		0
.equiv		mci1130_ii3,		4
.equiv		mci1130_ixO,		8
.equiv		mci1130_iyO,		12
.equiv		mci1130_izO,		16	
.equiv		mci1130_ixH,		20  
.equiv		mci1130_iyH,		28  
.equiv		mci1130_izH,		36  
.equiv		mci1130_qqOO,		44  
.equiv		mci1130_qqOH,		52  
.equiv		mci1130_qqHH,		60  
.equiv		mci1130_c6,		68  
.equiv		mci1130_c12,		76 
.equiv		mci1130_vctot,		84 
.equiv		mci1130_vnbtot,		92 
.equiv		mci1130_innerjjnr,	100
.equiv		mci1130_innerk,		104
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 108		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + mci1130_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci1130_charge]
	movd  mm1, [ebp + mci1130_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + mci1130_qqOO], mm4
	movq  [esp + mci1130_qqOH], mm5
	movq  [esp + mci1130_qqHH], mm6
	mov   edx, [ebp + mci1130_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + mci1130_ntype]
	add   edx, ecx
	mov   eax, [ebp + mci1130_nbfp]
	movd  mm0, [eax + edx*4]          
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + mci1130_c6], mm0
	movq  [esp + mci1130_c12], mm1
	
.mci1130_outer:
	mov   eax, [ebp + mci1130_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci1130_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci1130_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci1130_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci1130_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci1130_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci1130_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci1130_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + mci1130_ixO], mm5	
	movq  [esp + mci1130_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci1130_ixH], mm0	
	movq [esp + mci1130_iyH], mm1	
	movq [esp + mci1130_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci1130_vctot], mm7
	movq  [esp + mci1130_vnbtot], mm7

	mov   eax, [ebp + mci1130_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci1130_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci1130_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + mci1130_pos]
	mov   eax, [ebp + mci1130_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci1130_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci1130_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci1130_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci1130_innerjjnr],  4 ;# advance pointer 

	movd  mm6, [esp + mci1130_qqOO]
	movq  mm7, [esp + mci1130_qqOH]

	lea   eax, [eax + eax*2]
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci1130_ixO]
	pfsubr mm1, [esp + mci1130_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci1130_ixH]
	pfsubr mm3, [esp + mci1130_iyH]
	pfsubr mm4, [esp + mci1130_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	movq  mm4, mm1
	pfmul mm4, mm4		;# mm4=invsq  

	movq mm2, mm4
	pfmul mm2, mm4
	pfmul mm2, mm4
	movq mm0, mm2
	pfmul mm0,mm0
	pfmul mm2, [esp + mci1130_c6]
	pfmul mm0, [esp + mci1130_c12]
	movq mm5, mm0
	pfsub mm5, mm2		;# vnb 

	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 
	;# update nb potential 
	pfadd mm5, [esp + mci1130_vnbtot]
	movq [esp + mci1130_vnbtot], mm5

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
    pfrsqit1 mm5,mm3				
    pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + mci1130_vctot]
	movq [esp + mci1130_vctot], mm7
	
	;# interactions with j H1 
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	movd mm6, [esp + mci1130_qqOH]
	movq mm7, [esp + mci1130_qqHH]
	
	pfsubr mm0, [esp + mci1130_ixO]
	pfsubr mm1, [esp + mci1130_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci1130_ixH]
	pfsubr mm3, [esp + mci1130_iyH]
	pfsubr mm4, [esp + mci1130_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

	pfrsqrt mm1,mm0

	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 
	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + mci1130_vctot]
	movq [esp + mci1130_vctot], mm7
	
	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	movd mm6, [esp + mci1130_qqOH]
	movq mm7, [esp + mci1130_qqHH]

	pfsubr mm0, [esp + mci1130_ixO]
	pfsubr mm1, [esp + mci1130_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci1130_ixH]
	pfsubr mm3, [esp + mci1130_iyH]
	pfsubr mm4, [esp + mci1130_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 

	pfrsqrt mm1,mm0

	movq mm2,mm1
	pfmul mm1,mm1
	pfrsqit1 mm1,mm0				
	pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	;# calculate potential and scalar force 
	pfmul mm6, mm1		;# mm6=vcoul 

	pfrsqrt mm5, mm3
	pswapd mm3,mm3
	pfrsqrt mm2, mm3
	pswapd mm3,mm3
	punpckldq mm5,mm2	;# seeds are in mm5 now, and rsq in mm3. 

	movq mm2, mm5
	pfmul mm5,mm5
	pfrsqit1 mm5,mm3				
	pfrcpit2 mm5,mm2	;# mm5=invsqrt 
	pfmul mm7, mm5		;# mm7=vcoul 

	;# update vctot 
	pfadd mm7, mm6
	pfadd mm7, [esp + mci1130_vctot]
	movq [esp + mci1130_vctot], mm7
		
	;#  done  - one more? 
	dec dword ptr [esp + mci1130_innerk]
	jz  .mci1130_updateouterdata
	jmp .mci1130_inner_loop	
.mci1130_updateouterdata:	
	mov   edx, [ebp + mci1130_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci1130_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci1130_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + mci1130_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci1130_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + mci1130_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnbtot[gid] 
	;# finish if last 
	dec dword ptr [ebp + mci1130_nri]
	jz  .mci1130_end
	;# not last, iterate once more! 
	jmp .mci1130_outer
.mci1130_end:
	femms
	add esp, 108
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret




.globl mcinl3000_3dnow
.globl _mcinl3000_3dnow
mcinl3000_3dnow:	
_mcinl3000_3dnow:	
.equiv		mci3000_nri,		8
.equiv		mci3000_iinr,		12
.equiv		mci3000_jindex,		16
.equiv		mci3000_jjnr,		20
.equiv		mci3000_shift,		24
.equiv		mci3000_shiftvec,	28
.equiv		mci3000_gid,		32
.equiv		mci3000_pos,		36		
.equiv		mci3000_charge,		40
.equiv		mci3000_facel,		44
.equiv		mci3000_Vc,		48			
.equiv		mci3000_tabscale,	52
.equiv		mci3000_VFtab,		56
	;# stack offsets for local variables 
.equiv		mci3000_is3,		0
.equiv		mci3000_ii3,		4
.equiv		mci3000_ix,		8
.equiv		mci3000_iy,		12
.equiv		mci3000_iz,		16
.equiv		mci3000_iq,		20 
.equiv		mci3000_vctot,		28 
.equiv		mci3000_n1,		36
.equiv		mci3000_tsc,		44 
.equiv		mci3000_ntia,		52
.equiv		mci3000_innerjjnr,	56
.equiv		mci3000_innerk,		60						
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 64		;# local stack space 
	femms
	;# move data to local stack  
	movd  mm3, [ebp + mci3000_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci3000_tsc], mm3	
	;# assume we have at least one i particle - start directly 	
.mci3000_outer:
	mov   eax, [ebp + mci3000_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3000_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3000_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3000_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + mci3000_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3000_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3000_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci3000_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3000_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3000_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3000_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + mci3000_ix], mm0	
	movd  [esp + mci3000_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3000_vctot],  mm7

	mov   eax, [ebp + mci3000_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3000_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + mci3000_pos]	
	mov   eax, [ebp + mci3000_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3000_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + mci3000_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3000_unroll_loop
	jmp   .mci3000_finish_inner
.mci3000_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3000_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3000_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3000_charge]    ;# base of charge[] 
	movq mm5, [esp + mci3000_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3000_pos]

	movq  mm0, [esp + mci3000_ix]
	movd  mm1, [esp + mci3000_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3000_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci3000_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3000_VFtab]
	mov ecx, [esp + mci3000_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3000_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5. 
	;# update vctot 
	pfadd mm5, [esp + mci3000_vctot]      ;# add the earlier value 
	movq [esp + mci3000_vctot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3000_innerk],  2
	jl    .mci3000_finish_inner
	jmp   .mci3000_unroll_loop
.mci3000_finish_inner:	
	and dword ptr [esp + mci3000_innerk],  1
	jnz  .mci3000_single_inner
	jmp  .mci3000_updateouterdata		
.mci3000_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3000_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci3000_charge]
	movd mm5, [esp + mci3000_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + mci3000_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3000_ix]
	movd  mm1, [esp + mci3000_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3000_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3000_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3000_VFtab]
	mov ecx, [esp + mci3000_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3000_vctot]      ;# add the earlier value 
	movq [esp + mci3000_vctot], mm5       ;# store the sum       
	
.mci3000_updateouterdata:	
	mov   edx, [ebp + mci3000_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3000_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3000_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3000_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	;# finish if last 
	mov   ecx, [ebp + mci3000_nri]
	dec ecx
	jecxz .mci3000_end
	;# not last, iterate once more! 
	mov [ebp + mci3000_nri], ecx
	jmp .mci3000_outer
.mci3000_end:
	femms
	add esp, 64
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret



	
.globl mcinl3010_3dnow
.globl _mcinl3010_3dnow
mcinl3010_3dnow:	
_mcinl3010_3dnow:	
.equiv		mci3010_nri,		8
.equiv		mci3010_iinr,		12
.equiv		mci3010_jindex,		16
.equiv		mci3010_jjnr,		20
.equiv		mci3010_shift,		24
.equiv		mci3010_shiftvec,	28
.equiv		mci3010_gid,		32
.equiv		mci3010_pos,		36	
.equiv		mci3010_charge,		40
.equiv		mci3010_facel,		44
.equiv		mci3010_Vc,		48
.equiv		mci3010_tabscale,	52		
.equiv		mci3010_VFtab,		56
.equiv		mci3010_nsatoms,	60		
	;# stack offsets for local variables 
.equiv		mci3010_is3,		0
.equiv		mci3010_ii3,		4
.equiv		mci3010_shX,		8
.equiv		mci3010_shY,		12 
.equiv		mci3010_shZ,		16	
.equiv		mci3010_ix,		20
.equiv		mci3010_iy,		24
.equiv		mci3010_iz,		28	
.equiv		mci3010_iq,		32 
.equiv		mci3010_vctot,		40 
.equiv		mci3010_n1,		48
.equiv		mci3010_tsc,		56 			
.equiv		mci3010_innerjjnr0,	64
.equiv		mci3010_innerk0,	68		
.equiv		mci3010_innerjjnr,	72
.equiv		mci3010_innerk,		76				
.equiv		mci3010_nscoul,		80
.equiv		mci3010_solnr,		84		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 88		;# local stack space 
	femms
	
	add dword ptr [ebp + mci3010_nsatoms],  8
	movd  mm3, [ebp + mci3010_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci3010_tsc], mm3
	
	;# assume we have at least one i particle - start directly 		
.mci3010_outer:
	mov   eax, [ebp + mci3010_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3010_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3010_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3010_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + mci3010_shX], mm0
	movd  [esp + mci3010_shZ], mm1

	mov   ecx, [ebp + mci3010_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3010_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + mci3010_nsatoms]
	mov   ecx, [eax]
	add dword ptr [ebp + mci3010_nsatoms],  12
	mov   [esp + mci3010_nscoul], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + mci3010_vctot], mm7
	mov   [esp + mci3010_solnr], ebx
	
	mov   eax, [ebp + mci3010_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3010_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + mci3010_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3010_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + mci3010_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + mci3010_pos]
	mov   ecx, [esp + mci3010_nscoul]
	cmp   ecx,  0
	jnz  .mci3010_mno_coul
	jmp  .mci3010_last_mno
.mci3010_mno_coul:				
	mov   ebx,  [esp + mci3010_solnr]
	inc   dword ptr [esp + mci3010_solnr]
	mov   edx, [ebp + mci3010_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci3010_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3010_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3010_pos]    ;# eax = base of pos[] 
	mov   [esp + mci3010_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci3010_shX]
	pfadd mm1, [esp + mci3010_shZ]
	movq  [esp + mci3010_ix], mm0	
	movd  [esp + mci3010_iz], mm1	

	mov   ecx, [esp + mci3010_innerjjnr0]
	mov   [esp + mci3010_innerjjnr], ecx
	mov   edx, [esp + mci3010_innerk0]
    sub   edx,  2
    mov   [esp + mci3010_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3010_unroll_coul_loop
	jmp   .mci3010_finish_coul_inner
.mci3010_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3010_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3010_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3010_charge]    ;# base of charge[] 
	movq mm5, [esp + mci3010_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3010_pos]

	movq  mm0, [esp + mci3010_ix]
	movd  mm1, [esp + mci3010_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3010_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci3010_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3010_VFtab]
	mov ecx, [esp + mci3010_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3010_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3010_vctot]      ;# add the earlier value 
	movq [esp + mci3010_vctot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3010_innerk],  2
	jl    .mci3010_finish_coul_inner
	jmp   .mci3010_unroll_coul_loop
.mci3010_finish_coul_inner:	
	and dword ptr [esp + mci3010_innerk],  1
	jnz  .mci3010_single_coul_inner
	jmp  .mci3010_updateouterdata_coul		
.mci3010_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3010_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci3010_charge]
	movd mm5, [esp + mci3010_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + mci3010_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3010_ix]
	movd  mm1, [esp + mci3010_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3010_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3010_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3010_VFtab]
	mov ecx, [esp + mci3010_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3010_vctot]      ;# add the earlier value 
	movq [esp + mci3010_vctot], mm5       ;# store the sum       
	
.mci3010_updateouterdata_coul:	
	;# loop back to mno 
	dec dword ptr [esp + mci3010_nscoul]
	jz  .mci3010_last_mno
	jmp .mci3010_mno_coul
.mci3010_last_mno:	
	mov   edx, [ebp + mci3010_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3010_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3010_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3010_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + mci3010_nri]
	dec ecx
	jecxz .mci3010_end
	;# not last, iterate once more! 
	mov [ebp + mci3010_nri], ecx
	jmp .mci3010_outer
.mci3010_end:
	femms
	add esp, 88
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	


.globl mcinl3020_3dnow
.globl _mcinl3020_3dnow
mcinl3020_3dnow:	
_mcinl3020_3dnow:	
.equiv		mci3020_nri,		8
.equiv		mci3020_iinr,		12
.equiv		mci3020_jindex,		16
.equiv		mci3020_jjnr,		20
.equiv		mci3020_shift,		24
.equiv		mci3020_shiftvec,	28
.equiv		mci3020_gid,		32
.equiv		mci3020_pos,		36		
.equiv		mci3020_charge,		40
.equiv		mci3020_facel,		44
.equiv		mci3020_Vc,		48			
.equiv		mci3020_tabscale,	52
.equiv		mci3020_VFtab,		56
			;# stack offsets for local variables 
.equiv		mci3020_is3,		0
.equiv		mci3020_ii3,		4
.equiv		mci3020_ixO,		8
.equiv		mci3020_iyO,		12
.equiv		mci3020_izO,		16	
.equiv		mci3020_ixH,		20  
.equiv		mci3020_iyH,		28  
.equiv		mci3020_izH,		36  
.equiv		mci3020_iqO,		44  
.equiv		mci3020_iqH,		52  
.equiv		mci3020_qqO,		60  
.equiv		mci3020_qqH,		68  
.equiv		mci3020_vctot,		76  
.equiv		mci3020_n1,		84  
.equiv		mci3020_tsc,		92 
.equiv		mci3020_innerjjnr,	100
.equiv		mci3020_innerk,		104
.equiv		mci3020_tmprsqH,	108 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 116		;# local stack space 
	femms

	mov   ecx, [ebp + mci3020_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3020_charge]
	movd  mm1, [ebp + mci3020_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + mci3020_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3020_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	movd  mm4, [ebp + mci3020_tabscale]
	punpckldq mm4,mm4	    ;# spread to both halves 
	movq  [esp + mci3020_tsc], mm4	      
	;# assume we have at least one i particle - start directly 	 
.mci3020_outer:
	mov   eax, [ebp + mci3020_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3020_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3020_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3020_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci3020_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3020_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3020_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3020_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + mci3020_ixO], mm5	
	movq  [esp + mci3020_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci3020_ixH], mm0	
	movq [esp + mci3020_iyH], mm1	
	movq [esp + mci3020_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3020_vctot], mm7

	mov   eax, [ebp + mci3020_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3020_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci3020_innerk], edx        

	mov   esi, [ebp + mci3020_pos]	
	mov   eax, [ebp + mci3020_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3020_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci3020_inner_loop:	
	;# a single j particle iteration 
	mov   eax, [esp + mci3020_innerjjnr]
	mov   eax, [eax]	 ;# eax=jnr offset 
    add dword ptr [esp + mci3020_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + mci3020_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + mci3020_iqO]
	pfmul mm7, [esp + mci3020_iqH]	 ;# mm6=qqO, mm7=qqH 
	movd [esp + mci3020_qqO], mm6
	movq [esp + mci3020_qqH], mm7
		
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3020_ixO]
	pfsubr mm1, [esp + mci3020_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3020_ixH]
	pfsubr mm3, [esp + mci3020_iyH]
	pfsubr mm4, [esp + mci3020_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3020_tmprsqH], mm3
	
    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + mci3020_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3020_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3020_VFtab]
	mov ecx, [esp + mci3020_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3020_qqO]	;# vcoul=qq*VV 
	;# update vctot directly 
	pfadd mm5, [esp + mci3020_vctot]
	movq [esp + mci3020_vctot], mm5
	
	;# now do the two hydrogens. 
	movq mm0, [esp + mci3020_tmprsqH] ;# mm0=rsqH 

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3020_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3020_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3020_VFtab]
	mov ecx, [esp + mci3020_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3020_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3020_qqH]	;# vcoul=qq*VV 

	;# update vctot 
	pfadd mm5, [esp + mci3020_vctot]
	movq [esp + mci3020_vctot], mm5
		
	;#  done  - one more? 
	dec dword ptr [esp + mci3020_innerk]
	jz  .mci3020_updateouterdata
	jmp .mci3020_inner_loop
.mci3020_updateouterdata:	
	mov   edx, [ebp + mci3020_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3020_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3020_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3020_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	;# finish if last 
	dec dword ptr [ebp + mci3020_nri]
	jz  .mci3020_end
	;# not last, iterate once more! 
	jmp .mci3020_outer
.mci3020_end:
	femms
	add esp, 116
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

	

.globl mcinl3030_3dnow
.globl _mcinl3030_3dnow
mcinl3030_3dnow:	
_mcinl3030_3dnow:	
.equiv		mci3030_nri,		8
.equiv		mci3030_iinr,		12
.equiv		mci3030_jindex,		16
.equiv		mci3030_jjnr,		20
.equiv		mci3030_shift,		24
.equiv		mci3030_shiftvec,	28
.equiv		mci3030_gid,		32
.equiv		mci3030_pos,		36		
.equiv		mci3030_charge,		40
.equiv		mci3030_facel,		44
.equiv		mci3030_Vc,			48			
.equiv		mci3030_tabscale,	52
.equiv		mci3030_VFtab,		56
			;# stack offsets for local variables 
.equiv		mci3030_is3,		0
.equiv		mci3030_ii3,		4
.equiv		mci3030_ixO,		8
.equiv		mci3030_iyO,		12
.equiv		mci3030_izO,		16	
.equiv		mci3030_ixH,		20  
.equiv		mci3030_iyH,		28  
.equiv		mci3030_izH,		36  
.equiv		mci3030_qqOO,		44  
.equiv		mci3030_qqOH,		52  
.equiv		mci3030_qqHH,		60 
.equiv		mci3030_n1,			68  
.equiv		mci3030_tsc,		76  
.equiv		mci3030_vctot,		84  
.equiv		mci3030_innerjjnr,	92
.equiv		mci3030_innerk,		96
.equiv		mci3030_tmprsqH,	100 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 108		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + mci3030_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3030_charge]
	movd  mm1, [ebp + mci3030_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + mci3030_qqOO], mm4
	movq  [esp + mci3030_qqOH], mm5
	movq  [esp + mci3030_qqHH], mm6
	movd  mm3, [ebp + mci3030_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci3030_tsc], mm3
.mci3030_outer:
	mov   eax, [ebp + mci3030_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3030_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3030_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3030_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci3030_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3030_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3030_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3030_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + mci3030_ixO], mm5	
	movq  [esp + mci3030_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci3030_ixH], mm0	
	movq [esp + mci3030_iyH], mm1	
	movq [esp + mci3030_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3030_vctot], mm7

	mov   eax, [ebp + mci3030_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3030_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci3030_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + mci3030_pos]	
	mov   eax, [ebp + mci3030_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3030_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci3030_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3030_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci3030_innerjjnr],  4 ;# advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3030_ixO]
	pfsubr mm1, [esp + mci3030_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3030_ixH]
	pfsubr mm3, [esp + mci3030_iyH]
	pfsubr mm4, [esp + mci3030_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3030_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	pfmul mm0, mm1		;# mm0=rsq  

	pfmul mm0, [esp + mci3030_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3030_VFtab]
	mov ecx, [esp + mci3030_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3030_qqOO]	;# vcoul=qq*VV 
	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + mci3030_vctot]
	movq [esp + mci3030_vctot], mm5
	
	;# time for hydrogens! 

	movq mm0, [esp + mci3030_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3030_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3030_VFtab]
	mov ecx, [esp + mci3030_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3030_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3030_qqOH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3030_vctot]
	movq [esp + mci3030_vctot], mm5
	
	;# interactions with j H1 

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3030_ixO]
	pfsubr mm1, [esp + mci3030_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3030_ixH]
	pfsubr mm3, [esp + mci3030_iyH]
	pfsubr mm4, [esp + mci3030_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3030_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1		;# mm0=rsq  
	
	pfmul mm0, [esp + mci3030_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3030_VFtab]
	mov ecx, [esp + mci3030_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3030_qqOH]	;# vcoul=qq*VV 

	;# update vctot  directly, force is moved to mm3 
	pfadd mm5, [esp + mci3030_vctot]
	movq [esp + mci3030_vctot], mm5
	
	movq mm0, [esp + mci3030_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3030_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3030_VFtab]
	mov ecx, [esp + mci3030_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3030_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3030_qqHH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3030_vctot]
	movq [esp + mci3030_vctot], mm5
	
	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + mci3030_ixO]
	pfsubr mm1, [esp + mci3030_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3030_ixH]
	pfsubr mm3, [esp + mci3030_iyH]
	pfsubr mm4, [esp + mci3030_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3030_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1

	pfmul mm0, [esp + mci3030_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3030_VFtab]
	mov ecx, [esp + mci3030_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3030_qqOH]	;# vcoul=qq*VV 

	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + mci3030_vctot]
	movq [esp + mci3030_vctot], mm5

	movq mm0, [esp + mci3030_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3030_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3030_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3030_VFtab]
	mov ecx, [esp + mci3030_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3030_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3030_qqHH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3030_vctot]
	movq [esp + mci3030_vctot], mm5
	
	;#  done  - one more? 
	dec dword ptr [esp + mci3030_innerk]
	jz  .mci3030_updateouterdata
	jmp .mci3030_inner_loop	
.mci3030_updateouterdata:	
	mov   edx, [ebp + mci3030_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3030_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3030_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + mci3030_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	;# finish if last 
	dec dword ptr [ebp + mci3030_nri]
	jz  .mci3030_end
	;# not last, iterate once more! 
	jmp .mci3030_outer
.mci3030_end:
	femms
	add esp, 108
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret




.globl mcinl3100_3dnow
.globl _mcinl3100_3dnow
mcinl3100_3dnow:	
_mcinl3100_3dnow:	
.equiv		mci3100_nri,		8
.equiv		mci3100_iinr,		12
.equiv		mci3100_jindex,		16
.equiv		mci3100_jjnr,		20
.equiv		mci3100_shift,		24
.equiv		mci3100_shiftvec,	28
.equiv		mci3100_gid,		32
.equiv		mci3100_pos,		36
.equiv		mci3100_charge,		40
.equiv		mci3100_facel,		44
.equiv		mci3100_Vc,			48			
.equiv		mci3100_type,		52
.equiv		mci3100_ntype,		56
.equiv		mci3100_nbfp,		60	
.equiv		mci3100_Vnb,		64
.equiv		mci3100_tabscale,	68
.equiv		mci3100_VFtab,		72
	;# stack offsets for local variables 
.equiv		mci3100_is3,		0 
.equiv		mci3100_ii3,		4
.equiv		mci3100_ix,			8
.equiv		mci3100_iy,			12
.equiv		mci3100_iz,			16
.equiv		mci3100_iq,			20 
.equiv		mci3100_vctot,		28 
.equiv		mci3100_vnbtot,		36 
.equiv		mci3100_c6,			44 
.equiv		mci3100_c12,		52
.equiv		mci3100_n1,			60 
.equiv		mci3100_tsc,		68 
.equiv		mci3100_ntia,		76
.equiv		mci3100_innerjjnr,	80
.equiv		mci3100_innerk,		84						
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 88		;# local stack space 
	femms
	;# move data to local stack  
	movd  mm3, [ebp + mci3100_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci3100_tsc], mm3	
	;# assume we have at least one i particle - start directly 	
.mci3100_outer:
	mov   eax, [ebp + mci3100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3100_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3100_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3100_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + mci3100_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3100_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3100_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci3100_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3100_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + mci3100_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci3100_ntype]
	shl   edx, 1
	mov   [esp + mci3100_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3100_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3100_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + mci3100_ix], mm0	
	movd  [esp + mci3100_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3100_vctot],  mm7
	movq  [esp + mci3100_vnbtot], mm7

	mov   eax, [ebp + mci3100_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3100_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + mci3100_pos]
	mov   eax, [ebp + mci3100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3100_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + mci3100_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3100_unroll_loop
	jmp   .mci3100_finish_inner
.mci3100_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3100_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3100_charge]    ;# base of charge[] 
	movq mm5, [esp + mci3100_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + mci3100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci3100_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci3100_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci3100_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci3100_c6], mm5
	movq [esp + mci3100_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3100_pos]

	movq  mm0, [esp + mci3100_ix]
	movd  mm1, [esp + mci3100_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3100_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci3100_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3100_VFtab]
	mov ecx, [esp + mci3100_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3100_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + mci3100_tsc]
	
	pfmul mm1, [esp + mci3100_c12]

	pfmul mm2, [esp + mci3100_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = vnb12-vnb6 
	;# update vctot 
	pfadd mm5, [esp + mci3100_vctot]      ;# add the earlier value 
	movq [esp + mci3100_vctot], mm5       ;# store the sum       
	;# update vnbtot 
	pfadd mm4, [esp + mci3100_vnbtot]      ;# add the earlier value 
	movq [esp + mci3100_vnbtot], mm4       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3100_innerk],  2
	jl    .mci3100_finish_inner
	jmp   .mci3100_unroll_loop
.mci3100_finish_inner:	
	and dword ptr [esp + mci3100_innerk],  1
	jnz  .mci3100_single_inner
	jmp  .mci3100_updateouterdata		
.mci3100_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3100_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci3100_charge]
	movd mm5, [esp + mci3100_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + mci3100_nbfp]
	mov ecx, [ebp + mci3100_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci3100_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci3100_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci3100_c12], mm5


	mov   esi, [ebp + mci3100_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3100_ix]
	movd  mm1, [esp + mci3100_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3100_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3100_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3100_VFtab]
	mov ecx, [esp + mci3100_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	;# at this point mm5 contains vcoul 

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + mci3100_tsc]
	
	pfmul mm1, [esp + mci3100_c12]

	pfmul mm2, [esp + mci3100_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = vnb12-vnb6 
	;# update vctot 
	pfadd mm5, [esp + mci3100_vctot]      ;# add the earlier value 
	movq [esp + mci3100_vctot], mm5       ;# store the sum       
	;# update vnbtot 
	pfadd mm4, [esp + mci3100_vnbtot]      ;# add the earlier value 
	movq [esp + mci3100_vnbtot], mm4       ;# store the sum       

.mci3100_updateouterdata:	
	mov   edx, [ebp + mci3100_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3100_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3100_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3100_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
 
	movq  mm7, [esp + mci3100_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3100_Vnb] 
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + mci3100_nri]
	dec ecx
	jecxz .mci3100_end
	;# not last, iterate once more! 
	mov [ebp + mci3100_nri], ecx
	jmp .mci3100_outer
.mci3100_end:
	femms
	add esp, 88
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret




.globl mcinl3110_3dnow
.globl _mcinl3110_3dnow
mcinl3110_3dnow:	
_mcinl3110_3dnow:	
.equiv		mci3110_nri,		8
.equiv		mci3110_iinr,		12
.equiv		mci3110_jindex,		16
.equiv		mci3110_jjnr,		20
.equiv		mci3110_shift,		24
.equiv		mci3110_shiftvec,	28
.equiv		mci3110_gid,		32
.equiv		mci3110_pos,		36		
.equiv		mci3110_charge,		40
.equiv		mci3110_facel,		44
.equiv		mci3110_Vc,			48			
.equiv		mci3110_type,		52
.equiv		mci3110_ntype,		56
.equiv		mci3110_nbfp,		60	
.equiv		mci3110_Vnb,		64
.equiv		mci3110_tabscale,	68
.equiv		mci3110_VFtab,		72
.equiv		mci3110_nsatoms,	76	
	;# stack offsets for local variables 
.equiv		mci3110_is3,		0
.equiv		mci3110_ii3,		4
.equiv		mci3110_shX,		8
.equiv		mci3110_shY,		12 
.equiv		mci3110_shZ,		16	
.equiv		mci3110_ix,			20
.equiv		mci3110_iy,			24
.equiv		mci3110_iz,			28	
.equiv		mci3110_iq,			32  
.equiv		mci3110_vctot,		40  
.equiv		mci3110_vnbtot,		48  
.equiv		mci3110_c6,			56  
.equiv		mci3110_c12,		64 
.equiv		mci3110_two,		72 
.equiv		mci3110_n1,			80  
.equiv		mci3110_tsc,		88 
.equiv		mci3110_ntia,		96
.equiv		mci3110_innerjjnr0,	104
.equiv		mci3110_innerk0,	108
.equiv		mci3110_innerjjnr,	112
.equiv		mci3110_innerk,		116
.equiv		mci3110_nsvdwc,		120
.equiv		mci3110_nscoul,		124
.equiv		mci3110_nsvdw,		128
.equiv		mci3110_solnr,		132		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 136		;# local stack space 
	femms
	movq  mm2, [mm_two]
	movd  mm3, [ebp + mci3110_tabscale]
	movq  [esp + mci3110_two],    mm2
	punpckldq mm3,mm3
	movq  [esp + mci3110_tsc], mm3	
	;# assume we have at least one i particle - start directly 		
.mci3110_outer:
	mov   eax, [ebp + mci3110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3110_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3110_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3110_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + mci3110_shX], mm0
	movd  [esp + mci3110_shZ], mm1

	mov   ecx, [ebp + mci3110_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3110_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + mci3110_nsatoms]
	add dword ptr [ebp + mci3110_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + mci3110_nsvdwc], edx
	mov   [esp + mci3110_nscoul], eax
	mov   [esp + mci3110_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + mci3110_vctot],  mm7
	movq  [esp + mci3110_vnbtot], mm7
	mov   [esp + mci3110_solnr],  ebx
	
	mov   eax, [ebp + mci3110_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3110_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + mci3110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3110_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + mci3110_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + mci3110_pos]
	
	mov   ecx, [esp + mci3110_nsvdwc]
	cmp   ecx,  0
	jnz   .mci3110_mno_vdwc
	jmp   .mci3110_testcoul
.mci3110_mno_vdwc:
	mov   ebx,  [esp + mci3110_solnr]
	inc   dword ptr [esp + mci3110_solnr]
	mov   edx, [ebp + mci3110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci3110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3110_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + mci3110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci3110_ntype]
	shl   edx, 1
	mov   [esp + mci3110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3110_pos]    ;# eax = base of pos[] 
	mov   [esp + mci3110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci3110_shX]
	pfadd mm1, [esp + mci3110_shZ]
	movq  [esp + mci3110_ix], mm0	
	movd  [esp + mci3110_iz], mm1	

	mov   ecx, [esp + mci3110_innerjjnr0]
	mov   [esp + mci3110_innerjjnr], ecx
	mov   edx, [esp + mci3110_innerk0]
    sub   edx,  2
    mov   [esp + mci3110_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3110_unroll_vdwc_loop
	jmp   .mci3110_finish_vdwc_inner
.mci3110_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3110_charge]    ;# base of charge[] 
	movq mm5, [esp + mci3110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + mci3110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci3110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci3110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci3110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci3110_c6], mm5
	movq [esp + mci3110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3110_pos]

	movq  mm0, [esp + mci3110_ix]
	movd  mm1, [esp + mci3110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3110_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci3110_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3110_VFtab]
	mov ecx, [esp + mci3110_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3110_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + mci3110_two]	;# two*Heps2 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + mci3110_tsc]
	
	pfmul mm1, [esp + mci3110_c12]

	pfmul mm2, [esp + mci3110_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = vnb12-vnb6 
	;# update vctot 
	pfadd mm5, [esp + mci3110_vctot]      ;# add the earlier value 
	movq [esp + mci3110_vctot], mm5       ;# store the sum       
	;# update vnbtot 
	pfadd mm4, [esp + mci3110_vnbtot]      ;# add the earlier value 
	movq [esp + mci3110_vnbtot], mm4       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3110_innerk],  2
	jl    .mci3110_finish_vdwc_inner
	jmp   .mci3110_unroll_vdwc_loop
.mci3110_finish_vdwc_inner:	
	and dword ptr [esp + mci3110_innerk],  1
	jnz  .mci3110_single_vdwc_inner
	jmp  .mci3110_updateouterdata_vdwc		
.mci3110_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci3110_charge]
	movd mm5, [esp + mci3110_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + mci3110_nbfp]
	mov ecx, [ebp + mci3110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci3110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci3110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci3110_c12], mm5


	mov   esi, [ebp + mci3110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3110_ix]
	movd  mm1, [esp + mci3110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3110_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3110_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3110_VFtab]
	mov ecx, [esp + mci3110_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + mci3110_two]	;# two*Heps2 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 
	
	movq mm1, mm0
	pfmul mm1,mm1 	;# mm1=invsq 
	movq mm2, mm1
	pfmul mm2,mm1
	pfmul mm2,mm1	;# mm2=rinvsix 
	movq  mm1,mm2
	pfmul mm1,mm1	;# mm1=rinvtwelve 
	
	pfmul mm3, [esp + mci3110_tsc]
	
	pfmul mm1, [esp + mci3110_c12]

	pfmul mm2, [esp + mci3110_c6]

	movq mm4, mm1
	pfsub mm4, mm2	;# mm4 = vnb12-vnb6 
	;# update vctot 
	pfadd mm5, [esp + mci3110_vctot]      ;# add the earlier value 
	movq [esp + mci3110_vctot], mm5       ;# store the sum       
	;# update vnbtot 
	pfadd mm4, [esp + mci3110_vnbtot]      ;# add the earlier value 
	movq [esp + mci3110_vnbtot], mm4       ;# store the sum       

.mci3110_updateouterdata_vdwc:	
	;# loop back to mno 
	dec dword ptr [esp + mci3110_nsvdwc]
	jz  .mci3110_testcoul
	jmp .mci3110_mno_vdwc
.mci3110_testcoul:	
	mov  ecx, [esp + mci3110_nscoul]
	cmp  ecx,  0
	jnz  .mci3110_mno_coul
	jmp  .mci3110_testvdw
.mci3110_mno_coul:
	mov   ebx,  [esp + mci3110_solnr]
	inc   dword ptr [esp + mci3110_solnr]
	mov   edx, [ebp + mci3110_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci3110_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3110_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3110_pos]    ;# eax = base of pos[] 
	mov   [esp + mci3110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci3110_shX]
	pfadd mm1, [esp + mci3110_shZ]
	movq  [esp + mci3110_ix], mm0	
	movd  [esp + mci3110_iz], mm1	

	mov   ecx, [esp + mci3110_innerjjnr0]
	mov   [esp + mci3110_innerjjnr], ecx
	mov   edx, [esp + mci3110_innerk0]
    sub   edx,  2
    mov   [esp + mci3110_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3110_unroll_coul_loop
	jmp   .mci3110_finish_coul_inner
.mci3110_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3110_charge]    ;# base of charge[] 
	movq mm5, [esp + mci3110_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3110_pos]

	movq  mm0, [esp + mci3110_ix]
	movd  mm1, [esp + mci3110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3110_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci3110_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3110_VFtab]
	mov ecx, [esp + mci3110_n1]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3110_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + mci3110_two]	;# two*Heps2 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3110_vctot]      ;# add the earlier value 
	movq [esp + mci3110_vctot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3110_innerk],  2
	jl    .mci3110_finish_coul_inner
	jmp   .mci3110_unroll_coul_loop
.mci3110_finish_coul_inner:	
	and dword ptr [esp + mci3110_innerk],  1
	jnz  .mci3110_single_coul_inner
	jmp  .mci3110_updateouterdata_coul		
.mci3110_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci3110_charge]
	movd mm5, [esp + mci3110_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + mci3110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3110_ix]
	movd  mm1, [esp + mci3110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3110_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3110_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3110_VFtab]
	mov ecx, [esp + mci3110_n1]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm7, [esp + mci3110_two]	;# two*Heps2 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3110_vctot]      ;# add the earlier value 
	movq [esp + mci3110_vctot], mm5       ;# store the sum       
	
.mci3110_updateouterdata_coul:	
	;# loop back to mno 
	dec dword ptr [esp + mci3110_nscoul]
	jz  .mci3110_testvdw
	jmp .mci3110_mno_coul
.mci3110_testvdw:	
	mov  ecx, [esp + mci3110_nsvdw]
	cmp  ecx,  0
	jnz  .mci3110_mno_vdw
	jmp  .mci3110_last_mno
.mci3110_mno_vdw:
	mov   ebx,  [esp + mci3110_solnr]
	inc   dword ptr [esp + mci3110_solnr]

	mov   edx, [ebp + mci3110_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci3110_ntype]
	shl   edx, 1
	mov   [esp + mci3110_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3110_pos]    ;# eax = base of pos[] 
	mov   [esp + mci3110_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci3110_shX]
	pfadd mm1, [esp + mci3110_shZ]
	movq  [esp + mci3110_ix], mm0	
	movd  [esp + mci3110_iz], mm1	

	mov   ecx, [esp + mci3110_innerjjnr0]
	mov   [esp + mci3110_innerjjnr], ecx
	mov   edx, [esp + mci3110_innerk0]
    sub   edx,  2
    mov   [esp + mci3110_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3110_unroll_vdw_loop
	jmp   .mci3110_finish_vdw_inner
.mci3110_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3110_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 	

	mov ecx, [ebp + mci3110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci3110_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci3110_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci3110_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci3110_c6], mm5
	movq [esp + mci3110_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3110_pos]

	movq  mm0, [esp + mci3110_ix]
	movd  mm1, [esp + mci3110_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrcp mm0, mm4	             ;# lookup reciprocal seed  
    pfrcp mm1, mm6
 
	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
                  	        	;# amd 3dnow N-R iteration to get full precision. 
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	
	;# mm4 now contains invsq,
	 ;# do potential and fscal 
	 
	movq  mm0, mm4
	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci3110_c12]
	pfmul mm4, [esp + mci3110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci3110_vnbtot]      ;# add the earlier value 
	movq [esp + mci3110_vnbtot], mm6       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3110_innerk],  2
	jl    .mci3110_finish_vdw_inner
	jmp   .mci3110_unroll_vdw_loop
.mci3110_finish_vdw_inner:	
	and dword ptr [esp + mci3110_innerk],  1
	jnz  .mci3110_single_vdw_inner
	jmp  .mci3110_updateouterdata_vdw		
.mci3110_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3110_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci3110_nbfp]
	mov ecx, [ebp + mci3110_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci3110_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci3110_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci3110_c12], mm5

	mov   esi, [ebp + mci3110_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3110_ix]
	movd  mm1, [esp + mci3110_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm4=rsq 
	
    pfrcp mm0,mm4
    pfrcpit1 mm4,mm0				
    pfrcpit2 mm4,mm0	;# mm4=invsq  
	;# calculate potentials and scalar force 
	movq  mm0, mm4

	pfmul mm4, mm0
	pfmul mm4, mm0             	;# mm4=rinvsix 
	movq  mm5, mm4	
	pfmul mm5, mm5	            ;# mm5=rinvtwelve 

	pfmul mm5, [esp + mci3110_c12]
	pfmul mm4, [esp + mci3110_c6]	
	movq mm6, mm5	;# mm6 is vnb12-vnb6  
	pfsub mm6, mm4
	;# update vnbtot 
	pfadd mm6, [esp + mci3110_vnbtot]      ;# add the earlier value 
	movq [esp + mci3110_vnbtot], mm6       ;# store the sum       

.mci3110_updateouterdata_vdw:	
	;# loop back to mno 
	dec dword ptr [esp + mci3110_nsvdw]
	jz  .mci3110_last_mno
	jmp .mci3110_mno_vdw
	
.mci3110_last_mno:	
	mov   edx, [ebp + mci3110_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3110_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3110_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3110_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci3110_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3110_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + mci3110_nri]
	dec ecx
	jecxz .mci3110_end
	;# not last, iterate once more! 
	mov [ebp + mci3110_nri], ecx
	jmp .mci3110_outer
.mci3110_end:
	femms
	add esp, 136
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


	

.globl mcinl3120_3dnow
.globl _mcinl3120_3dnow
mcinl3120_3dnow:	
_mcinl3120_3dnow:	
.equiv		mci3120_nri,		8
.equiv		mci3120_iinr,		12
.equiv		mci3120_jindex,		16
.equiv		mci3120_jjnr,		20
.equiv		mci3120_shift,		24
.equiv		mci3120_shiftvec,	28
.equiv		mci3120_gid,		32
.equiv		mci3120_pos,		36		
.equiv		mci3120_charge,		40
.equiv		mci3120_facel,		44
.equiv		mci3120_Vc,			48			
.equiv		mci3120_type,		52
.equiv		mci3120_ntype,		56
.equiv		mci3120_nbfp,		60	
.equiv		mci3120_Vnb,		64
.equiv		mci3120_tabscale,	68
.equiv		mci3120_VFtab,		72
			;# stack offsets for local variables 
.equiv		mci3120_is3,		0
.equiv		mci3120_ii3,		4
.equiv		mci3120_ixO,		8
.equiv		mci3120_iyO,		12
.equiv		mci3120_izO,		16	
.equiv		mci3120_ixH,		20  
.equiv		mci3120_iyH,		28  
.equiv		mci3120_izH,		36  
.equiv		mci3120_iqO,		44  
.equiv		mci3120_iqH,		52  
.equiv		mci3120_qqO,		60  
.equiv		mci3120_qqH,		68  
.equiv		mci3120_vctot,		76  
.equiv		mci3120_vnbtot,		84  
.equiv		mci3120_c6,			92  
.equiv		mci3120_c12,		100 
.equiv		mci3120_n1,			108
.equiv		mci3120_tsc,		116 
.equiv		mci3120_ntia,		124 
.equiv		mci3120_innerjjnr,	128
.equiv		mci3120_innerk,		132
.equiv		mci3120_tmprsqH,	136 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 144		;# local stack space 
	femms

	mov   ecx, [ebp + mci3120_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3120_charge]
	movd  mm1, [ebp + mci3120_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1
	movq  [esp + mci3120_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3120_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + mci3120_type] 	
	mov   edx, [edx + ebx*4]
	shl   edx, 1
	mov   ecx, edx		        
	imul  ecx, [ebp + mci3120_ntype]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + mci3120_ntia], ecx
	 	
	movq  mm6, [ebp + mci3120_tabscale]
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + mci3120_tsc], mm6	      
 	;# assume we have at least one i particle - start directly 	
.mci3120_outer:
	mov   eax, [ebp + mci3120_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3120_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3120_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3120_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci3120_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3120_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3120_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3120_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + mci3120_ixO], mm5	
	movq  [esp + mci3120_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci3120_ixH], mm0	
	movq [esp + mci3120_iyH], mm1	
	movq [esp + mci3120_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3120_vctot], mm7
	movq  [esp + mci3120_vnbtot], mm7

	mov   eax, [ebp + mci3120_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3120_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci3120_innerk], edx        

	mov   esi, [ebp + mci3120_pos]	
	mov   eax, [ebp + mci3120_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3120_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci3120_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3120_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci3120_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + mci3120_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + mci3120_iqO]
	pfmul mm7, [esp + mci3120_iqH]	;# mm6=qqO, mm7=qqH 
	movd [esp + mci3120_qqO], mm6
	movq [esp + mci3120_qqH], mm7

	mov ecx, [ebp + mci3120_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + mci3120_nbfp]
	shl edx, 1
	add edx, [esp + mci3120_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci3120_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci3120_c12], mm5	
			
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3120_ixO]
	pfsubr mm1, [esp + mci3120_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3120_ixH]
	pfsubr mm3, [esp + mci3120_iyH]
	pfsubr mm4, [esp + mci3120_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3120_tmprsqH], mm3
	
    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + mci3120_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3120_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3120_VFtab]
	mov ecx, [esp + mci3120_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3120_qqO]	;# vcoul=qq*VV 
	;# update vctot directly 
	pfadd mm5, [esp + mci3120_vctot]
	movq [esp + mci3120_vctot], mm5
	
	;# nontabulated LJ - mm1 is invsqrt. - keep mm1! 
	movq mm0, mm1
	pfmul mm0, mm0		;# mm0 is invsq 
	movq mm2, mm0
	pfmul mm2, mm0
	pfmul mm2, mm0		;# mm2 = rinvsix 
	movq mm4, mm2
	pfmul mm4, mm4		;# mm4=rinvtwelve 

	pfmul mm4, [esp + mci3120_c12]
	pfmul mm2, [esp + mci3120_c6]
	pfsub mm4, mm2		;# mm4=vnb12-vnb6 

	;# update vnbtot  
	pfadd mm4, [esp + mci3120_vnbtot]      ;# add the earlier value 
	movq [esp + mci3120_vnbtot], mm4       ;# store the sum       

		;# now do the two hydrogens. 
	movq mm0, [esp + mci3120_tmprsqH] ;# mm0=rsqH 

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3120_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3120_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3120_VFtab]
	mov ecx, [esp + mci3120_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3120_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 


	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3120_qqH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3120_vctot]
	movq [esp + mci3120_vctot], mm5

	;#  done  - one more? 
	dec dword ptr [esp + mci3120_innerk]
	jz  .mci3120_updateouterdata
	jmp .mci3120_inner_loop
.mci3120_updateouterdata:	

	mov   edx, [ebp + mci3120_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3120_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3120_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3120_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci3120_vnbtot]     
	pfacc mm7,mm7	          ;# same for Vnb 
	
	mov   eax, [ebp + mci3120_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 
	;# finish if last 
	dec dword ptr [ebp + mci3120_nri]
	jz  .mci3120_end
	;# not last, iterate once more! 
	jmp .mci3120_outer
.mci3120_end:
	femms
	add esp, 144
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret





.globl mcinl3130_3dnow
.globl _mcinl3130_3dnow
mcinl3130_3dnow:	
_mcinl3130_3dnow:	
.equiv		mci3130_nri,		8
.equiv		mci3130_iinr,		12
.equiv		mci3130_jindex,		16
.equiv		mci3130_jjnr,		20
.equiv		mci3130_shift,		24
.equiv		mci3130_shiftvec,	28
.equiv		mci3130_gid,		32
.equiv		mci3130_pos,		36		
.equiv		mci3130_charge,		40
.equiv		mci3130_facel,		44
.equiv		mci3130_Vc,			48			
.equiv		mci3130_type,		52
.equiv		mci3130_ntype,		56
.equiv		mci3130_nbfp,		60	
.equiv		mci3130_Vnb,		64
.equiv		mci3130_tabscale,	68
.equiv		mci3130_VFtab,		72
	;# stack offsets for local variables 
.equiv		mci3130_is3,		0
.equiv		mci3130_ii3,		4
.equiv		mci3130_ixO,		8
.equiv		mci3130_iyO,		12
.equiv		mci3130_izO,		16	
.equiv		mci3130_ixH,		20  
.equiv		mci3130_iyH,		28  
.equiv		mci3130_izH,		36  
.equiv		mci3130_qqOO,		44  
.equiv		mci3130_qqOH,		52  
.equiv		mci3130_qqHH,		60  
.equiv		mci3130_c6,			68  
.equiv		mci3130_c12,		76 
.equiv		mci3130_n1,			84
.equiv		mci3130_tsc,		92 
.equiv		mci3130_vctot,		100 
.equiv		mci3130_vnbtot,		108 
.equiv		mci3130_innerjjnr,	116
.equiv		mci3130_innerk,		120
.equiv		mci3130_tmprsqH,	124 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 132		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + mci3130_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3130_charge]
	movd  mm1, [ebp + mci3130_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + mci3130_qqOO], mm4
	movq  [esp + mci3130_qqOH], mm5
	movq  [esp + mci3130_qqHH], mm6
	mov   edx, [ebp + mci3130_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + mci3130_ntype]
	add   edx, ecx
	mov   eax, [ebp + mci3130_nbfp]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + mci3130_c6], mm0
	movq  [esp + mci3130_c12], mm1
	movd  mm5, [ebp + mci3130_tabscale]
	punpckldq mm5,mm5
	movq  [esp + mci3130_tsc], mm5
.mci3130_outer:
	mov   eax, [ebp + mci3130_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3130_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3130_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3130_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci3130_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3130_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3130_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3130_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + mci3130_ixO], mm5	
	movq  [esp + mci3130_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci3130_ixH], mm0	
	movq [esp + mci3130_iyH], mm1	
	movq [esp + mci3130_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3130_vctot], mm7
	movq  [esp + mci3130_vnbtot], mm7

	mov   eax, [ebp + mci3130_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3130_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci3130_innerk], edx    ;# number of innerloop atoms 

	mov   esi, [ebp + mci3130_pos]
	mov   eax, [ebp + mci3130_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3130_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci3130_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3130_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci3130_innerjjnr],  4 ;# advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3130_ixO]
	pfsubr mm1, [esp + mci3130_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3130_ixH]
	pfsubr mm3, [esp + mci3130_iyH]
	pfsubr mm4, [esp + mci3130_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3130_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	pfmul mm0, mm1		;# mm0=rsq  

	pfmul mm0, [esp + mci3130_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3130_VFtab]
	mov ecx, [esp + mci3130_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3130_qqOO]	;# vcoul=qq*VV 

	;# update vctot directly 
	pfadd mm5, [esp + mci3130_vctot]
	movq [esp + mci3130_vctot], mm5
	
	movq mm5, mm1
	pfmul mm5,mm5
	movq mm4, mm5
	pfmul mm4,mm5
	pfmul mm4,mm5
	movq mm5, mm4
	pfmul mm5,mm5	;# mm4=rinvsix, mm5=rinvtwelve 

	pfmul mm4, [esp + mci3130_c6]
	pfmul mm5, [esp + mci3130_c12]
	movq mm6,mm5
	pfsub mm6,mm4

	;# update vnbtot  
	pfadd mm6, [esp + mci3130_vnbtot]      ;# add the earlier value 
	movq [esp + mci3130_vnbtot], mm6       ;# store the sum       
	
	;# time for hydrogens! 

	movq mm0, [esp + mci3130_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3130_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3130_VFtab]
	mov ecx, [esp + mci3130_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3130_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3130_qqOH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3130_vctot]
	movq [esp + mci3130_vctot], mm5
	
	;# interactions with j H1 

	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3130_ixO]
	pfsubr mm1, [esp + mci3130_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3130_ixH]
	pfsubr mm3, [esp + mci3130_iyH]
	pfsubr mm4, [esp + mci3130_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3130_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1		;# mm0=rsq  
	
	pfmul mm0, [esp + mci3130_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3130_VFtab]
	mov ecx, [esp + mci3130_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3130_qqOH]	;# vcoul=qq*VV 

	;# update vctot  directly, force is moved to mm3 
	pfadd mm5, [esp + mci3130_vctot]
	movq [esp + mci3130_vctot], mm5
	
	movq mm0, [esp + mci3130_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3130_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3130_VFtab]
	mov ecx, [esp + mci3130_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3130_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3130_qqHH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3130_vctot]
	movq [esp + mci3130_vctot], mm5
	
	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + mci3130_ixO]
	pfsubr mm1, [esp + mci3130_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3130_ixH]
	pfsubr mm3, [esp + mci3130_iyH]
	pfsubr mm4, [esp + mci3130_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3130_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1

	pfmul mm0, [esp + mci3130_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3130_VFtab]
	mov ecx, [esp + mci3130_n1]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3130_qqOH]	;# vcoul=qq*VV 

	;# update vctot directly 
	pfadd mm5, [esp + mci3130_vctot]
	movq [esp + mci3130_vctot], mm5
	
	movq mm0, [esp + mci3130_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3130_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3130_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3130_VFtab]
	mov ecx, [esp + mci3130_n1]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3130_n1 + 4]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3130_qqHH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3130_vctot]
	movq [esp + mci3130_vctot], mm5
		
	;#  done  - one more? 
	dec dword ptr [esp + mci3130_innerk]
	jz  .mci3130_updateouterdata
	jmp .mci3130_inner_loop	
.mci3130_updateouterdata:	
	mov   edx, [ebp + mci3130_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3130_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3130_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + mci3130_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci3130_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + mci3130_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnbtot[gid] 
	;# finish if last 
	dec dword ptr [ebp + mci3130_nri]
	jz  .mci3130_end
	;# not last, iterate once more! 
	jmp .mci3130_outer
.mci3130_end:
	femms
	add esp, 132
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl mcinl3300_3dnow
.globl _mcinl3300_3dnow
mcinl3300_3dnow:	
_mcinl3300_3dnow:	
.equiv		mci3300_nri,		8
.equiv		mci3300_iinr,		12
.equiv		mci3300_jindex,		16
.equiv		mci3300_jjnr,		20
.equiv		mci3300_shift,		24
.equiv		mci3300_shiftvec,	28
.equiv		mci3300_gid,		32
.equiv		mci3300_pos,		36		
.equiv		mci3300_charge,		40
.equiv		mci3300_facel,		44
.equiv		mci3300_Vc,			48			
.equiv		mci3300_type,		52
.equiv		mci3300_ntype,		56
.equiv		mci3300_nbfp,		60	
.equiv		mci3300_Vnb,		64
.equiv		mci3300_tabscale,	68
.equiv		mci3300_VFtab,		72
	;# stack offsets for local variables 
.equiv		mci3300_is3,		0
.equiv		mci3300_ii3,		4
.equiv		mci3300_ix,			8
.equiv		mci3300_iy,			12
.equiv		mci3300_iz,			16
.equiv		mci3300_iq,			20  
.equiv		mci3300_vctot,		28  
.equiv		mci3300_vnbtot,		36  
.equiv		mci3300_c6,			44  
.equiv		mci3300_c12,		52  
.equiv		mci3300_n1,			60
.equiv		mci3300_tsc,		68  
.equiv		mci3300_ntia,		76
.equiv		mci3300_innerjjnr,	80
.equiv		mci3300_innerk,		84					
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 88		;# local stack space 
	femms
	;# move data to local stack  
	movd  mm3, [ebp + mci3300_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci3300_tsc], mm3	
	;# assume we have at least one i particle - start directly 	
.mci3300_outer:
	mov   eax, [ebp + mci3300_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3300_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3300_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3300_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]

	mov   ecx, [ebp + mci3300_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3300_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3300_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci3300_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3300_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + mci3300_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci3300_ntype]
	shl   edx, 1
	mov   [esp + mci3300_ntia], edx

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3300_pos]    ;# eax = base of pos[] 
	
	pfadd mm0, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm3, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3300_ii3], ebx
	pfadd mm1, mm3
	movq  [esp + mci3300_ix], mm0	
	movd  [esp + mci3300_iz], mm1	
				
	;# clear total potential and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3300_vctot],  mm7
	movq  [esp + mci3300_vnbtot], mm7

	mov   eax, [ebp + mci3300_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3300_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + mci3300_pos]
	mov   eax, [ebp + mci3300_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3300_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	sub   edx,  2
	mov   [esp + mci3300_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3300_unroll_loop
	jmp   .mci3300_finish_inner
.mci3300_unroll_loop:
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3300_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3300_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3300_charge]    ;# base of charge[] 
	movq mm5, [esp + mci3300_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + mci3300_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci3300_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci3300_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci3300_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci3300_c6], mm5
	movq [esp + mci3300_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3300_pos]

	movq  mm0, [esp + mci3300_ix]
	movd  mm1, [esp + mci3300_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci3300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3300_VFtab]
	mov ecx, [esp + mci3300_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3300_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3300_vctot]      ;# add the earlier value 
	movq [esp + mci3300_vctot], mm5       ;# store the sum       

	;# dispersion table 
	mov ecx, [esp + mci3300_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + mci3300_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci3300_c6]
	pfmul mm5, mm4	;# vnb6    
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci3300_vnbtot]      ;# add the earlier value 
	movq [esp + mci3300_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + mci3300_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + mci3300_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci3300_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci3300_vnbtot]      ;# add the earlier value 
	movq [esp + mci3300_vnbtot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3300_innerk],  2
	jl    .mci3300_finish_inner
	jmp   .mci3300_unroll_loop
.mci3300_finish_inner:	
	and dword ptr [esp + mci3300_innerk],  1
	jnz  .mci3300_single_inner
	jmp  .mci3300_updateouterdata		
.mci3300_single_inner:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3300_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci3300_charge]
	movd mm5, [esp + mci3300_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + mci3300_nbfp]
	mov ecx, [ebp + mci3300_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci3300_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci3300_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci3300_c12], mm5

	mov   esi, [ebp + mci3300_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3300_ix]
	movd  mm1, [esp + mci3300_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3300_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3300_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3300_VFtab]
	mov ecx, [esp + mci3300_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3300_vctot]      ;# add the earlier value 
	movq [esp + mci3300_vctot], mm5       ;# store the sum       
	
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci3300_c6]
	pfmul mm5, mm4	;# vnb6  

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci3300_vnbtot]      ;# add the earlier value 
	movq [esp + mci3300_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci3300_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci3300_vnbtot]      ;# add the earlier value 
	movq [esp + mci3300_vnbtot], mm5       ;# store the sum       

.mci3300_updateouterdata:	
	mov   edx, [ebp + mci3300_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3300_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3300_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3300_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci3300_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3300_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 

	;# finish if last 
	mov   ecx, [ebp + mci3300_nri]
	dec ecx
	jecxz .mci3300_end
	;# not last, iterate once more! 
	mov [ebp + mci3300_nri], ecx
	jmp .mci3300_outer
.mci3300_end:
	femms
	add esp, 88
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret





.globl mcinl3310_3dnow
.globl _mcinl3310_3dnow
mcinl3310_3dnow:	
_mcinl3310_3dnow:	
.equiv		mci3310_nri,		8
.equiv		mci3310_iinr,		12
.equiv		mci3310_jindex,		16
.equiv		mci3310_jjnr,		20
.equiv		mci3310_shift,		24
.equiv		mci3310_shiftvec,	28
.equiv		mci3310_gid,		32
.equiv		mci3310_pos,		36		
.equiv		mci3310_charge,		40
.equiv		mci3310_facel,		44
.equiv		mci3310_Vc,			48			
.equiv		mci3310_type,		52
.equiv		mci3310_ntype,		56
.equiv		mci3310_nbfp,		60	
.equiv		mci3310_Vnb,		64
.equiv		mci3310_tabscale,	68
.equiv		mci3310_VFtab,		72
.equiv		mci3310_nsatoms,	76		
	;# stack offsets for local variables 
.equiv		mci3310_is3,		0
.equiv		mci3310_ii3,		4
.equiv		mci3310_shX,		8
.equiv		mci3310_shY,		12 
.equiv		mci3310_shZ,		16	
.equiv		mci3310_ix,			20
.equiv		mci3310_iy,			24
.equiv		mci3310_iz,			28	
.equiv		mci3310_iq,			32  
.equiv		mci3310_vctot,		40  
.equiv		mci3310_vnbtot,		48  
.equiv		mci3310_c6,			56  
.equiv		mci3310_c12,		64 
.equiv		mci3310_n1,			72  
.equiv		mci3310_tsc,		80  
.equiv		mci3310_ntia,		88	
.equiv		mci3310_innerjjnr0,	92
.equiv		mci3310_innerk0,	96	
.equiv		mci3310_innerjjnr,	100
.equiv		mci3310_innerk,		104
.equiv		mci3310_nsvdwc,		108
.equiv		mci3310_nscoul,		112
.equiv		mci3310_nsvdw,		116
.equiv		mci3310_solnr,		120		
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 124		;# local stack space 
	femms
	movd  mm3, [ebp + mci3310_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci3310_tsc], mm3	
	;# assume we have at least one i particle - start directly 		
.mci3310_outer:
	mov   eax, [ebp + mci3310_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3310_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3310_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3310_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm0, [eax + ebx*4]	;# move shX/shY to mm0 and shZ to mm1 
	movd  mm1, [eax + ebx*4 + 8]
	movq  [esp + mci3310_shX], mm0
	movd  [esp + mci3310_shZ], mm1

	mov   ecx, [ebp + mci3310_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3310_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   eax, [ebp + mci3310_nsatoms]
	add dword ptr [ebp + mci3310_nsatoms],  12
	mov   ecx, [eax]	
	mov   edx, [eax + 4]
	mov   eax, [eax + 8]	
	sub   ecx, eax
	sub   eax, edx
	
	mov   [esp + mci3310_nsvdwc], edx
	mov   [esp + mci3310_nscoul], eax
	mov   [esp + mci3310_nsvdw], ecx
		
	;# clear potential 
	pxor  mm7,mm7
	movq  [esp + mci3310_vctot],  mm7
	movq  [esp + mci3310_vnbtot], mm7
	mov   [esp + mci3310_solnr],  ebx
	
	mov   eax, [ebp + mci3310_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3310_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   eax, [ebp + mci3310_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3310_innerjjnr0], eax     ;# pointer to jjnr[nj0] 

	mov   [esp + mci3310_innerk0], edx    ;# number of innerloop atoms 
	mov   esi, [ebp + mci3310_pos]
	
	mov   ecx, [esp + mci3310_nsvdwc]
	cmp   ecx,  0
	jnz   .mci3310_mno_vdwc
	jmp   .mci3310_testcoul
.mci3310_mno_vdwc:
	mov   ebx,  [esp + mci3310_solnr]
	inc   dword ptr [esp + mci3310_solnr]
	mov   edx, [ebp + mci3310_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci3310_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3310_iq], mm2	    ;# iq =facel*charge[ii] 

	mov   edx, [ebp + mci3310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci3310_ntype]
	shl   edx, 1
	mov   [esp + mci3310_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3310_pos]    ;# eax = base of pos[] 
	mov   [esp + mci3310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci3310_shX]
	pfadd mm1, [esp + mci3310_shZ]
	movq  [esp + mci3310_ix], mm0	
	movd  [esp + mci3310_iz], mm1	

	mov   ecx, [esp + mci3310_innerjjnr0]
	mov   [esp + mci3310_innerjjnr], ecx
	mov   edx, [esp + mci3310_innerk0]
    sub   edx,  2
    mov   [esp + mci3310_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3310_unroll_vdwc_loop
	jmp   .mci3310_finish_vdwc_inner
.mci3310_unroll_vdwc_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3310_charge]    ;# base of charge[] 
	movq mm5, [esp + mci3310_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	mov ecx, [ebp + mci3310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci3310_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci3310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci3310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6,mm5			
	punpckldq mm5,mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6,mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci3310_c6], mm5
	movq [esp + mci3310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3310_pos]

	movq  mm0, [esp + mci3310_ix]
	movd  mm1, [esp + mci3310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3310_VFtab]
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3310_vctot]      ;# add the earlier value 
	movq [esp + mci3310_vctot], mm5       ;# store the sum       

	;# dispersion table 
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + mci3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci3310_c6]
	pfmul mm5, mm4	;# vnb6   

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci3310_vnbtot]      ;# add the earlier value 
	movq [esp + mci3310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]
	mov ecx, [esp + mci3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 32]
	punpckldq mm5, [edx + ecx*4 + 36]
	punpckldq mm6, [edx + ecx*4 + 40]
	punpckldq mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci3310_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci3310_vnbtot]      ;# add the earlier value 
	movq [esp + mci3310_vnbtot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3310_innerk],  2
	jl    .mci3310_finish_vdwc_inner
	jmp   .mci3310_unroll_vdwc_loop
.mci3310_finish_vdwc_inner:	
	and dword ptr [esp + mci3310_innerk],  1
	jnz  .mci3310_single_vdwc_inner
	jmp  .mci3310_updateouterdata_vdwc		
.mci3310_single_vdwc_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci3310_charge]
	movd mm5, [esp + mci3310_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov esi, [ebp + mci3310_nbfp]
	mov ecx, [ebp + mci3310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci3310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci3310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci3310_c12], mm5

	mov   esi, [ebp + mci3310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3310_ix]
	movd  mm1, [esp + mci3310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3310_VFtab]
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3310_vctot]      ;# add the earlier value 
	movq [esp + mci3310_vctot], mm5       ;# store the sum       
	
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci3310_c6]
	pfmul mm5, mm4	;# vnb6            
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci3310_vnbtot]      ;# add the earlier value 
	movq [esp + mci3310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci3310_c12]
	pfmul mm5, mm6	;# vnb12 

	;# change sign of mm3 
    pxor mm1,mm1
	pfsub mm1, mm3	
	pfmul mm0, [esp + mci3310_tsc]
 	pfmul mm0, mm1    ;# mm0 is total fscal now 	

	;# update vnbtot 
	pfadd mm5, [esp + mci3310_vnbtot]      ;# add the earlier value 
	movq [esp + mci3310_vnbtot], mm5       ;# store the sum       

.mci3310_updateouterdata_vdwc:	
	;# loop back to mno 
	dec dword ptr [esp + mci3310_nsvdwc]
	jz  .mci3310_testcoul
	jmp .mci3310_mno_vdwc
.mci3310_testcoul:	
	mov  ecx, [esp + mci3310_nscoul]
	cmp  ecx,  0
	jnz  .mci3310_mno_coul
	jmp  .mci3310_testvdw
.mci3310_mno_coul:
	mov   ebx,  [esp + mci3310_solnr]
	inc   dword ptr [esp + mci3310_solnr]
	mov   edx, [ebp + mci3310_charge]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii] 
	pfmul mm2, [ebp + mci3310_facel]
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3310_iq], mm2	    ;# iq =facel*charge[ii] 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3310_pos]    ;# eax = base of pos[] 
	mov   [esp + mci3310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci3310_shX]
	pfadd mm1, [esp + mci3310_shZ]
	movq  [esp + mci3310_ix], mm0	
	movd  [esp + mci3310_iz], mm1	

	mov   ecx, [esp + mci3310_innerjjnr0]
	mov   [esp + mci3310_innerjjnr], ecx
	mov   edx, [esp + mci3310_innerk0]
    sub   edx,  2
    mov   [esp + mci3310_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3310_unroll_coul_loop
	jmp   .mci3310_finish_coul_inner
.mci3310_unroll_coul_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3310_charge]    ;# base of charge[] 
	movq mm5, [esp + mci3310_iq]
	movd mm3, [ecx + eax*4]	     ;# charge[jnr1] 
    punpckldq mm3, [ecx + ebx*4]     ;# move charge 2 to high part of mm3 
	pfmul mm3,mm5		     ;# mm3 now has qq for both particles 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3310_pos]

	movq  mm0, [esp + mci3310_ix]
	movd  mm1, [esp + mci3310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6,mm0	             ;# dr = ir - jr  
	pfsubr mm7,mm1
	pfmul mm6,mm6	             ;# square dx,dy,dz 
	pfmul mm7,mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0,mm1
	punpckldq mm4,mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2,mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movq [esp + mci3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3310_VFtab]
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# coulomb table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3310_vctot]      ;# add the earlier value 
	movq [esp + mci3310_vctot], mm5       ;# store the sum       

	;# should we do one more iteration? 
	sub dword ptr [esp + mci3310_innerk],  2
	jl    .mci3310_finish_coul_inner
	jmp   .mci3310_unroll_coul_loop
.mci3310_finish_coul_inner:	
	and dword ptr [esp + mci3310_innerk],  1
	jnz  .mci3310_single_coul_inner
	jmp  .mci3310_updateouterdata_coul		
.mci3310_single_coul_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov ecx, [ebp + mci3310_charge]
	movd mm5, [esp + mci3310_iq]
	movd mm3, [ecx + eax*4]
	pfmul mm3, mm5	  	;# mm3=qq 

	mov   esi, [ebp + mci3310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3310_ix]
	movd  mm1, [esp + mci3310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3310_VFtab]
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 

	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, mm3	;# vcoul=qq*VV 

	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	pfadd mm5, [esp + mci3310_vctot]      ;# add the earlier value 
	movq [esp + mci3310_vctot], mm5       ;# store the sum       
	
.mci3310_updateouterdata_coul:	
	;# loop back to mno 
	dec dword ptr [esp + mci3310_nscoul]
	jz  .mci3310_testvdw
	jmp .mci3310_mno_coul
.mci3310_testvdw:	
	mov  ecx, [esp + mci3310_nsvdw]
	cmp  ecx,  0
	jnz  .mci3310_mno_vdw
	jmp  .mci3310_last_mno
.mci3310_mno_vdw:
	mov   ebx,  [esp + mci3310_solnr]
	inc   dword ptr [esp + mci3310_solnr]

	mov   edx, [ebp + mci3310_type] 	
	mov   edx, [edx + ebx*4]	
	imul  edx, [ebp + mci3310_ntype]
	shl   edx, 1
	mov   [esp + mci3310_ntia], edx	

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3310_pos]    ;# eax = base of pos[] 
	mov   [esp + mci3310_ii3], ebx
	
	movq  mm0, [eax + ebx*4]
	movd  mm1, [eax + ebx*4 + 8]
	pfadd mm0, [esp + mci3310_shX]
	pfadd mm1, [esp + mci3310_shZ]
	movq  [esp + mci3310_ix], mm0	
	movd  [esp + mci3310_iz], mm1	

	mov   ecx, [esp + mci3310_innerjjnr0]
	mov   [esp + mci3310_innerjjnr], ecx
	mov   edx, [esp + mci3310_innerk0]
    sub   edx,  2
    mov   [esp + mci3310_innerk], edx    ;# number of innerloop atoms 
	jge   .mci3310_unroll_vdw_loop
	jmp   .mci3310_finish_vdw_inner
.mci3310_unroll_vdw_loop:	
	;# paired innerloop starts here 
	mov   ecx, [esp + mci3310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]         ;# eax/ebx=jnr 
	add dword ptr [esp + mci3310_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	prefetch [ecx + 16]	     ;# prefetch data - trial and error says 16 is best 
	
	mov ecx, [ebp + mci3310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	mov ecx, [ecx + ebx*4]       ;# type [jnr2] 

	mov esi, [ebp + mci3310_nbfp]		;# base of nbfp  
	shl edx, 1
	shl ecx, 1
	add edx, [esp + mci3310_ntia]	     ;# tja = ntia + 2*type 
	add ecx, [esp + mci3310_ntia]

	movq mm5, [esi + edx*4]		;# mm5 = 1st c6 / c12 		
	movq mm7, [esi + ecx*4]		;# mm7 = 2nd c6 / c12 	
	movq mm6, mm5			
	punpckldq mm5, mm7		;# mm5 = 1st c6 / 2nd c6 
	punpckhdq mm6, mm7		;# mm6 = 1st c12 / 2nd c12 
	movq [esp + mci3310_c6], mm5
	movq [esp + mci3310_c12], mm6

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]		

	mov   esi, [ebp + mci3310_pos]

	movq  mm0, [esp + mci3310_ix]
	movd  mm1, [esp + mci3310_iz]	 	
	movq  mm4, [esi + eax*4]     ;# fetch first j coordinates 
	movd  mm5, [esi + eax*4 + 8]		
	pfsubr mm4,mm0		     ;# dr = ir - jr  
	pfsubr mm5,mm1
	pfmul mm4,mm4	             ;# square dx,dy,dz 		         
	pfmul mm5,mm5		
	pfacc mm4, mm5               ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm4, mm5		     ;# first rsq in lower mm4 

	movq  mm6, [esi + ebx*4]     ;# fetch second j coordinates  
	movd  mm7, [esi + ebx*4 + 8]
	
	pfsubr mm6, mm0	             ;# dr = ir - jr  
	pfsubr mm7, mm1
	pfmul mm6, mm6	             ;# square dx,dy,dz 
	pfmul mm7, mm7
	pfacc mm6, mm7		     ;# accumulate to get dx*dx+ dy*dy+ dz*dz 
	pfacc mm6, mm7	             ;# second rsq in lower mm6 

    pfrsqrt mm0, mm4	     ;# lookup inverse square root seed 
    pfrsqrt mm1, mm6
 

	punpckldq mm0, mm1
	punpckldq mm4, mm6        	;# now 4 has rsq and 0 the seed for both pairs. 
    movq mm2, mm0	        	;# amd 3dnow N-R iteration to get full precision. 
	pfmul mm0, mm0
    pfrsqit1 mm0, mm4				
    pfrcpit2 mm0, mm2	
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 
	;# do potential and fscal 
	pfmul mm1, [esp + mci3310_tsc]	;# mm1=rt 
	pf2iw mm4, mm1
	movq [esp + mci3310_n1], mm4
	pi2fd mm4, mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 is n0 
	
	movq mm2, mm1
	pfmul mm2, mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3310_VFtab]
	;# dispersion table 
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]	
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci3310_c6]
	pfmul mm5, mm4	;# vnb6            

	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci3310_vnbtot]      ;# add the earlier value 
	movq [esp + mci3310_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	mov ecx, [esp + mci3310_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4 + 16]
	punpckldq mm5, [edx + ecx*4 + 20]
	punpckldq mm6, [edx + ecx*4 + 24]
	punpckldq mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci3310_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci3310_vnbtot]      ;# add the earlier value 
	movq [esp + mci3310_vnbtot], mm5       ;# store the sum       
	
	;# should we do one more iteration? 
	sub dword ptr [esp + mci3310_innerk],  2
	jl    .mci3310_finish_vdw_inner
	jmp   .mci3310_unroll_vdw_loop
.mci3310_finish_vdw_inner:	
	and dword ptr [esp + mci3310_innerk],  1
	jnz  .mci3310_single_vdw_inner
	jmp  .mci3310_updateouterdata_vdw		
.mci3310_single_vdw_inner:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3310_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 

	mov esi, [ebp + mci3310_nbfp]
	mov ecx, [ebp + mci3310_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr1] 
	shl edx, 1
	add edx, [esp + mci3310_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [esi + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci3310_c6], mm5
	movd mm5, [esi + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci3310_c12], mm5

	mov   esi, [ebp + mci3310_pos]
	lea   eax, [eax + eax*2]

	movq  mm0, [esp + mci3310_ix]
	movd  mm1, [esp + mci3310_iz]
	movq  mm4, [esi + eax*4]
	movd  mm5, [esi + eax*4 + 8]
	pfsubr mm4, mm0
	pfsubr mm5, mm1
	pfmul mm4,mm4
	pfmul mm5,mm5
	pfacc mm4, mm5
	pfacc mm4, mm5		;# mm0=rsq 
	
    pfrsqrt mm0,mm4
    movq mm2,mm0
    pfmul mm0,mm0
    pfrsqit1 mm0,mm4				
    pfrcpit2 mm0,mm2	;# mm1=invsqrt 
	pfmul mm4, mm0
	movq mm1, mm4
	;# mm0 is invsqrt, and mm1 r. 

	;# calculate potentials and scalar force 
	pfmul mm1, [esp + mci3310_tsc]	;# mm1=rt 
	pf2iw mm4,mm1
	movd [esp + mci3310_n1], mm4
	pi2fd mm4,mm4
	pfsub mm1, mm4               ;# now mm1 is eps and mm4 n0. 

	movq mm2,mm1
	pfmul mm2,mm2	;# mm1 is eps, mm2 is eps2 
	
	mov edx, [ebp + mci3310_VFtab]
	mov ecx, [esp + mci3310_n1]
	lea ecx, [ecx + ecx*2]	
	shl ecx, 2
	;# dispersion table
    ;# load all the table values we need
	 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci3310_c6]
	pfmul mm5, mm4	;# vnb6            
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci3310_vnbtot]      ;# add the earlier value 
	movq [esp + mci3310_vnbtot], mm5       ;# store the sum       

	;# repulsion table
    ;# load all the table values we need
	 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]

	pfmul mm6, mm1  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm1  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci3310_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot 
	pfadd mm5, [esp + mci3310_vnbtot]      ;# add the earlier value 
	movq [esp + mci3310_vnbtot], mm5       ;# store the sum       

.mci3310_updateouterdata_vdw:	
	;# loop back to mno 
	dec dword ptr [esp + mci3310_nsvdw]
	jz  .mci3310_last_mno
	jmp .mci3310_mno_vdw
	
.mci3310_last_mno:	
	mov   edx, [ebp + mci3310_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3310_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3310_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3310_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci3310_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3310_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 
	;# finish if last 
	mov   ecx, [ebp + mci3310_nri]
	dec ecx
	jecxz .mci3310_end
	;# not last, iterate once more! 
	mov [ebp + mci3310_nri], ecx
	jmp .mci3310_outer
.mci3310_end:
	femms
	add esp, 124
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


.globl mcinl3320_3dnow
.globl _mcinl3320_3dnow
mcinl3320_3dnow:	
_mcinl3320_3dnow:	
.equiv		mci3320_nri,		8
.equiv		mci3320_iinr,		12
.equiv		mci3320_jindex,		16
.equiv		mci3320_jjnr,		20
.equiv		mci3320_shift,		24
.equiv		mci3320_shiftvec,	28
.equiv		mci3320_gid,		32
.equiv		mci3320_pos,		36		
.equiv		mci3320_charge,		40
.equiv		mci3320_facel,		44
.equiv		mci3320_Vc,			48			
.equiv		mci3320_type,		52
.equiv		mci3320_ntype,		56
.equiv		mci3320_nbfp,		60	
.equiv		mci3320_Vnb,		64
.equiv		mci3320_tabscale,	68
.equiv		mci3320_VFtab,		72
			;# stack offsets for local variables 
.equiv		mci3320_is3,		0
.equiv		mci3320_ii3,		4
.equiv		mci3320_ixO,		8
.equiv		mci3320_iyO,		12
.equiv		mci3320_izO,		16	
.equiv		mci3320_ixH,		20  
.equiv		mci3320_iyH,		28  
.equiv		mci3320_izH,		36  
.equiv		mci3320_iqO,		44  
.equiv		mci3320_iqH,		52  
.equiv		mci3320_qqO,		60  
.equiv		mci3320_qqH,		68  
.equiv		mci3320_vctot,		76  
.equiv		mci3320_vnbtot,		84  
.equiv		mci3320_c6,			92  
.equiv		mci3320_c12,		100
.equiv		mci3320_n1,			108 
.equiv		mci3320_tsc,		116 
.equiv		mci3320_ntia,		124 
.equiv		mci3320_innerjjnr,	128
.equiv		mci3320_innerk,		132
.equiv		mci3320_tmprsqH,	136 
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 144		;# local stack space 
	femms

	mov   ecx, [ebp + mci3320_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3320_charge]
	movd  mm1, [ebp + mci3320_facel]
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] 
	pfmul mm2, mm1		
	movq  [esp + mci3320_iqO], mm2	    ;# iqO = facel*charge[ii] 
	
	movd  mm2, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] 
	pfmul mm2, mm1
	punpckldq mm2,mm2	    ;# spread to both halves 
	movq  [esp + mci3320_iqH], mm2	    ;# iqH = facel*charge[ii0+1] 

	mov   edx, [ebp + mci3320_type] 	
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1		        
	imul  ecx, [ebp + mci3320_ntype]      ;# ecx = ntia = 2*ntype*type[ii0]  
	mov   [esp + mci3320_ntia], ecx
	 	
	movq  mm4, [ebp + mci3320_tabscale]
	punpckldq mm4,mm4	    ;# spread to both halves 
	movq  [esp + mci3320_tsc], mm4	      
	;# assume we have at least one i particle - start directly 	 
.mci3320_outer:
	mov   eax, [ebp + mci3320_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3320_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3320_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3320_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci3320_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3320_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3320_pos]    ;# eax = base of pos[] 
	
	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3320_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + mci3320_ixO], mm5	
	movq  [esp + mci3320_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci3320_ixH], mm0	
	movq [esp + mci3320_iyH], mm1	
	movq [esp + mci3320_izH], mm2	
					
	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3320_vctot], mm7
	movq  [esp + mci3320_vnbtot], mm7

	mov   eax, [ebp + mci3320_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3320_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci3320_innerk], edx        

	mov   esi, [ebp + mci3320_pos]
	mov   eax, [ebp + mci3320_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3320_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci3320_inner_loop:	
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3320_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci3320_innerjjnr],  4 ;# advance pointer 
	prefetch [ecx + 16]	   ;# prefetch data - trial and error says 16 is best 

	mov ecx, [ebp + mci3320_charge]
	movd mm7, [ecx + eax*4]
	punpckldq mm7,mm7
	movq mm6,mm7
	pfmul mm6, [esp + mci3320_iqO]
	pfmul mm7, [esp + mci3320_iqH]	;# mm6=qqO, mm7=qqH 
	movd [esp + mci3320_qqO], mm6
	movq [esp + mci3320_qqH], mm7

	mov ecx, [ebp + mci3320_type]
	mov edx, [ecx + eax*4]        	 ;# type [jnr] 
	mov ecx, [ebp + mci3320_nbfp]
	shl edx, 1
	add edx, [esp + mci3320_ntia]	     ;# tja = ntia + 2*type 
	movd mm5, [ecx + edx*4]		;# mm5 = 1st c6  		
	movq [esp + mci3320_c6], mm5
	movd mm5, [ecx + edx*4 + 4]	;# mm5 = 1st c12  		
	movq [esp + mci3320_c12], mm5	
			
	lea   eax, [eax + eax*2]
	
	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3320_ixO]
	pfsubr mm1, [esp + mci3320_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3320_ixH]
	pfsubr mm3, [esp + mci3320_iyH]
	pfsubr mm4, [esp + mci3320_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3320_tmprsqH], mm3
	
    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 

	pfmul mm0, mm1		;# mm0=r 

	pfmul mm0, [esp + mci3320_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3320_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3320_VFtab]
	mov ecx, [esp + mci3320_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3320_qqO]	;# vcoul=qq*VV 

	;# update vctot directly 
	pfadd mm5, [esp + mci3320_vctot]
	movq [esp + mci3320_vctot], mm5
	
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci3320_c6]
	pfmul mm5, mm4	;# vnb6            
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci3320_vnbtot]      ;# add the earlier value 
	movq [esp + mci3320_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci3320_c12]
	pfmul mm5, mm6	;# vnb12 
	;# update vnbtot  
	pfadd mm5, [esp + mci3320_vnbtot]      ;# add the earlier value 
	movq [esp + mci3320_vnbtot], mm5       ;# store the sum       
	
	;# now do the two hydrogens. 
	movq mm0, [esp + mci3320_tmprsqH] ;# mm0=rsqH 

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3320_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3320_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3320_VFtab]
	mov ecx, [esp + mci3320_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3320_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3320_qqH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3320_vctot]
	movq [esp + mci3320_vctot], mm5
	
	;#  done  - one more? 
	dec dword ptr [esp + mci3320_innerk]
	jz  .mci3320_updateouterdata
	jmp .mci3320_inner_loop
.mci3320_updateouterdata:	
	mov   edx, [ebp + mci3320_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3320_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3320_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 
	
	mov   eax, [ebp + mci3320_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci3320_vnbtot]     
	pfacc mm7,mm7	          ;# same for Vnb 
	
	mov   eax, [ebp + mci3320_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnb[gid] 
	;# finish if last 
	dec dword ptr [ebp + mci3320_nri]
	jz  .mci3320_end
	;# not last, iterate once more! 
	jmp .mci3320_outer
.mci3320_end:
	femms
	add esp, 144
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret
 
	

.globl mcinl3330_3dnow
.globl _mcinl3330_3dnow
mcinl3330_3dnow:	
_mcinl3330_3dnow:	
.equiv		mci3330_nri,		8
.equiv		mci3330_iinr,		12
.equiv		mci3330_jindex,		16
.equiv		mci3330_jjnr,		20
.equiv		mci3330_shift,		24
.equiv		mci3330_shiftvec,	28
.equiv		mci3330_gid,		32
.equiv		mci3330_pos,		36		
.equiv		mci3330_charge,		40
.equiv		mci3330_facel,		44
.equiv		mci3330_Vc,			48			
.equiv		mci3330_type,		52
.equiv		mci3330_ntype,		56
.equiv		mci3330_nbfp,		60	
.equiv		mci3330_Vnb,		64
.equiv		mci3330_tabscale,	68
.equiv		mci3330_VFtab,		72
			;# stack offsets for local variables 
.equiv		mci3330_is3,		0
.equiv		mci3330_ii3,		4
.equiv		mci3330_ixO,		8
.equiv		mci3330_iyO,		12
.equiv		mci3330_izO,		16	
.equiv		mci3330_ixH,		20  
.equiv		mci3330_iyH,		28  
.equiv		mci3330_izH,		36  
.equiv		mci3330_qqOO,		44  
.equiv		mci3330_qqOH,		52  
.equiv		mci3330_qqHH,		60  
.equiv		mci3330_c6,			68  
.equiv		mci3330_c12,		76 
.equiv		mci3330_n1,			84  
.equiv		mci3330_tsc,		92 
.equiv		mci3330_vctot,		100 
.equiv		mci3330_vnbtot,		108 
.equiv		mci3330_innerjjnr,	116
.equiv		mci3330_innerk,		120
.equiv		mci3330_tmprsqH,	124  
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 132		;# local stack space 
	femms
	;# assume we have at least one i particle - start directly 	

	mov   ecx, [ebp + mci3330_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx]	    ;# ebx=ii 

	mov   edx, [ebp + mci3330_charge]
	movd  mm1, [ebp + mci3330_facel]	;# mm1=facel 
	movd  mm2, [edx + ebx*4]    ;# mm2=charge[ii0] (O) 
	movd  mm3, [edx + ebx*4 + 4]    ;# mm2=charge[ii0+1] (H)  
	movq  mm4, mm2	
	pfmul mm4, mm1
	movq  mm6, mm3
	pfmul mm6, mm1
	movq  mm5, mm4
	pfmul mm4, mm2			;# mm4=qqOO*facel 
	pfmul mm5, mm3			;# mm5=qqOH*facel 
	pfmul mm6, mm3			;# mm6=qqHH*facel 
	punpckldq mm5,mm5	    ;# spread to both halves 
	punpckldq mm6,mm6	    ;# spread to both halves 
	movq  [esp + mci3330_qqOO], mm4
	movq  [esp + mci3330_qqOH], mm5
	movq  [esp + mci3330_qqHH], mm6
	mov   edx, [ebp + mci3330_type]
	mov   ecx, [edx + ebx*4]
	shl   ecx, 1
	mov   edx, ecx
	imul  ecx, [ebp + mci3330_ntype]
	add   edx, ecx
	mov   eax, [ebp + mci3330_nbfp]
	movd  mm0, [eax + edx*4]
	movd  mm1, [eax + edx*4 + 4]
	movq  [esp + mci3330_c6], mm0
	movq  [esp + mci3330_c12], mm1
	movd  mm3, [ebp + mci3330_tabscale]
	punpckldq mm3,mm3
	movq  [esp + mci3330_tsc], mm3
.mci3330_outer:
	mov   eax, [ebp + mci3330_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax]		;# ebx=shift[n] 
	add dword ptr [ebp + mci3330_shift],  4  ;# advance pointer one step 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + mci3330_is3],ebx    	;# store is3 

	mov   eax, [ebp + mci3330_shiftvec]   ;# eax = base of shiftvec[] 
	
	movq  mm5, [eax + ebx*4]	;# move shX/shY to mm5 and shZ to mm6. 
	movd  mm6, [eax + ebx*4 + 8]
	movq  mm0, mm5
	movq  mm1, mm5
	movq  mm2, mm6
	punpckldq mm0,mm0	    ;# also expand shX,Y,Z in mm0--mm2. 
	punpckhdq mm1,mm1
	punpckldq mm2,mm2		
	
	mov   ecx, [ebp + mci3330_iinr]       ;# ecx = pointer into iinr[] 	
	add dword ptr [ebp + mci3330_iinr],  4   ;# advance pointer 
	mov   ebx, [ecx]	    ;# ebx=ii 

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + mci3330_pos]    ;# eax = base of pos[] 

	pfadd mm5, [eax + ebx*4]    ;# ix = shX + posX (and iy too) 
	movd  mm7, [eax + ebx*4 + 8]    ;# cant use direct memory add for 4 bytes (iz) 
	mov   [esp + mci3330_ii3], ebx	    ;# (use mm7 as temp. storage for iz.) 
	pfadd mm6, mm7
	movq  [esp + mci3330_ixO], mm5	
	movq  [esp + mci3330_izO], mm6

	movd  mm3, [eax + ebx*4 + 12]
	movd  mm4, [eax + ebx*4 + 16]
	movd  mm5, [eax + ebx*4 + 20]
	punpckldq  mm3, [eax + ebx*4 + 24]
	punpckldq  mm4, [eax + ebx*4 + 28]
	punpckldq  mm5, [eax + ebx*4 + 32] ;# coords of H1 in low mm3-mm5, H2 in high 
	
	pfadd mm0, mm3
	pfadd mm1, mm4
	pfadd mm2, mm5		
	movq [esp + mci3330_ixH], mm0	
	movq [esp + mci3330_iyH], mm1	
	movq [esp + mci3330_izH], mm2	

	;# clear vctot and i forces 
	pxor  mm7,mm7
	movq  [esp + mci3330_vctot], mm7
	movq  [esp + mci3330_vnbtot], mm7

	mov   eax, [ebp + mci3330_jindex]
	mov   ecx, [eax]	     ;# jindex[n] 
	mov   edx, [eax + 4]	     ;# jindex[n+1] 
	add dword ptr [ebp + mci3330_jindex],  4
	sub   edx, ecx               ;# number of innerloop atoms 
	mov   [esp + mci3330_innerk], edx        

	mov   esi, [ebp + mci3330_pos]	
	mov   eax, [ebp + mci3330_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + mci3330_innerjjnr], eax     ;# pointer to jjnr[nj0] 
.mci3330_inner_loop:
	;# a single j particle iteration here - compare with the unrolled code for comments. 
	mov   eax, [esp + mci3330_innerjjnr]
	mov   eax, [eax]	;# eax=jnr offset 
    add dword ptr [esp + mci3330_innerjjnr],  4 ;# advance pointer 

	lea   eax, [eax + eax*2]

	movq  mm0, [esi + eax*4]
	movd  mm1, [esi + eax*4 + 8]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3330_ixO]
	pfsubr mm1, [esp + mci3330_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm0
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3330_ixH]
	pfsubr mm3, [esp + mci3330_iyH]
	pfsubr mm4, [esp + mci3330_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3330_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt  
	pfmul mm0, mm1		;# mm0=rsq  

	pfmul mm0, [esp + mci3330_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3330_VFtab]
	mov ecx, [esp + mci3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3330_qqOO]	;# vcoul=qq*VV 
	;# update vctot directly, use mm3 for fscal sum. 
	pfadd mm5, [esp + mci3330_vctot]
	movq [esp + mci3330_vctot], mm5
	
	;# dispersion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 16]
	movd mm5, [edx + ecx*4 + 20]
	movd mm6, [edx + ecx*4 + 24]
	movd mm7, [edx + ecx*4 + 28]
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 	

	movq mm4, [esp + mci3330_c6]
	pfmul mm5, mm4	;# vnb6            
	;# update vnbtot to release mm5! 
	pfadd mm5, [esp + mci3330_vnbtot]      ;# add the earlier value 
	movq [esp + mci3330_vnbtot], mm5       ;# store the sum       

	;# repulsion table 
	;# load all the table values we need 
	movd mm4, [edx + ecx*4 + 32]
	movd mm5, [edx + ecx*4 + 36]
	movd mm6, [edx + ecx*4 + 40]
	movd mm7, [edx + ecx*4 + 44]

	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 
	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	movq mm6, [esp + mci3330_c12]
	pfmul mm5, mm6	;# vnb12 
	;# change sign of fscal and multiply with rinv  
 	;# update vnbtot  
	pfadd mm5, [esp + mci3330_vnbtot]      ;# add the earlier value 
	movq [esp + mci3330_vnbtot], mm5       ;# store the sum       
	
	;# Ready with the oxygen - time for hydrogens 
	
	movq mm0, [esp + mci3330_tmprsqH]

	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3330_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3330_VFtab]
	mov ecx, [esp + mci3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3330_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3330_qqOH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3330_vctot]
	movq [esp + mci3330_vctot], mm5

	;# interactions with j H1 
	movq  mm0, [esi + eax*4 + 12]
	movd  mm1, [esi + eax*4 + 20]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4
	
	pfsubr mm0, [esp + mci3330_ixO]
	pfsubr mm1, [esp + mci3330_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3330_ixH]
	pfsubr mm3, [esp + mci3330_iyH]
	pfsubr mm4, [esp + mci3330_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3330_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1		;# mm0=rsq  
	
	pfmul mm0, [esp + mci3330_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3330_VFtab]
	mov ecx, [esp + mci3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3330_qqOH]	;# vcoul=qq*VV 
	;# update vctot directly 
	pfadd mm5, [esp + mci3330_vctot]
	movq [esp + mci3330_vctot], mm5
	
	movq mm0, [esp + mci3330_tmprsqH]
	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3330_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3330_VFtab]
	mov ecx, [esp + mci3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3330_n1 + 4]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3330_qqHH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3330_vctot]
	movq [esp + mci3330_vctot], mm5
	
	;# interactions with j H2 
	movq  mm0, [esi + eax*4 + 24]
	movd  mm1, [esi + eax*4 + 32]
	;# copy & expand to mm2-mm4 for the H interactions 
	movq  mm2, mm0
	movq  mm3, mm0
	movq  mm4, mm1
	punpckldq mm2,mm2
	punpckhdq mm3,mm3
	punpckldq mm4,mm4

	pfsubr mm0, [esp + mci3330_ixO]
	pfsubr mm1, [esp + mci3330_izO]
		
	pfmul mm0,mm0
	pfmul mm1,mm1
	pfacc mm0, mm1
	pfadd mm0, mm1		;# mm0=rsqO 
	
	punpckldq mm2, mm2
	punpckldq mm3, mm3
	punpckldq mm4, mm4  ;# mm2-mm4 is jx-jz 
	pfsubr mm2, [esp + mci3330_ixH]
	pfsubr mm3, [esp + mci3330_iyH]
	pfsubr mm4, [esp + mci3330_izH] ;# mm2-mm4 is dxH-dzH 
	
	pfmul mm2,mm2
	pfmul mm3,mm3
	pfmul mm4,mm4

	pfadd mm3,mm2
	pfadd mm3,mm4		;# mm3=rsqH 
	movq [esp + mci3330_tmprsqH], mm3

    pfrsqrt mm1,mm0

    movq mm2,mm1
    pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	pfmul mm0, mm1

	pfmul mm0, [esp + mci3330_tsc]
	pf2iw mm4, mm0
	movd [esp + mci3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 

	;# coulomb table 
	mov edx, [ebp + mci3330_VFtab]
	mov ecx, [esp + mci3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2

	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3330_qqOH]	;# vcoul=qq*VV 
	;# update vctot directly 
	pfadd mm5, [esp + mci3330_vctot]
	movq [esp + mci3330_vctot], mm5
	
	movq mm0, [esp + mci3330_tmprsqH]
	pfrsqrt mm1, mm0
	pswapd mm0,mm0
	pfrsqrt mm2, mm0
	pswapd mm0,mm0
	punpckldq mm1,mm2	;# seeds are in mm1 now, and rsq in mm0. 

	movq mm2, mm1
	pfmul mm1,mm1
    pfrsqit1 mm1,mm0				
    pfrcpit2 mm1,mm2	;# mm1=invsqrt 
	
	pfmul mm0,mm1		;# mm0=r 
	pfmul mm0, [esp + mci3330_tsc]
	pf2iw mm4, mm0
	movq [esp + mci3330_n1], mm4
	pi2fd mm4,mm4
	pfsub mm0, mm4               ;# now mm0 is eps and mm4 n0 
	movq  mm2, mm0
	pfmul mm2, mm2		;# mm0 is eps, mm2 eps2 
	
	;# coulomb table 
	mov edx, [ebp + mci3330_VFtab]
	mov ecx, [esp + mci3330_n1]
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	;# load all values we need 
	movd mm4, [edx + ecx*4]
	movd mm5, [edx + ecx*4 + 4]
	movd mm6, [edx + ecx*4 + 8]
	movd mm7, [edx + ecx*4 + 12]
	mov ecx, [esp + mci3330_n1 + 4] ;# mm5 = Fp 
	lea ecx, [ecx + ecx*2]
	shl ecx, 2
	punpckldq mm4, [edx + ecx*4]
	punpckldq mm5, [edx + ecx*4 + 4]
	punpckldq mm6, [edx + ecx*4 + 8]
	punpckldq mm7, [edx + ecx*4 + 12]

	
	pfmul mm6, mm0  ;# mm6 = Geps 		
	pfmul mm7, mm2	;# mm7 = Heps2 
	
	pfadd mm5, mm6
	pfadd mm5, mm7	;# mm5 = Fp 

	pfmul mm5, mm0  ;# mm5=eps*Fp 
	pfadd mm5, mm4	;#  mm5= VV 

	pfmul mm5, [esp + mci3330_qqHH]	;# vcoul=qq*VV 
	;# update vctot 
	pfadd mm5, [esp + mci3330_vctot]
	movq [esp + mci3330_vctot], mm5
		
	;#  done  - one more? 
	dec dword ptr [esp + mci3330_innerk]
	jz  .mci3330_updateouterdata
	jmp .mci3330_inner_loop	
.mci3330_updateouterdata:	
	mov   edx, [ebp + mci3330_gid]      ;# get group index for this i particle 
	mov   edx, [edx]
	add dword ptr [ebp + mci3330_gid],  4  ;# advance pointer 

	movq  mm7, [esp + mci3330_vctot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + mci3330_Vc]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vc[gid] 

	movq  mm7, [esp + mci3330_vnbtot]     
	pfacc mm7,mm7	          ;# get and sum the two parts of total potential 

	mov   eax, [ebp + mci3330_Vnb]
	movd  mm6, [eax + edx*4] 
	pfadd mm6, mm7
	movd  [eax + edx*4], mm6          ;# increment vnbtot[gid] 
	;# finish if last 
	dec dword ptr [ebp + mci3330_nri]
	jz  .mci3330_end
	;# not last, iterate once more! 
	jmp .mci3330_outer
.mci3330_end:
	femms
	add esp, 132
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret
 

