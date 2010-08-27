;#
;#
;# Gromacs 4.0                         Copyright (c) 1991-2003 
;# David van der Spoel, Erik Lindahl
;#
;# This program is free software; you can redistribute it and/or
;# modify it under the terms of the GNU General Public License
;# as published by the Free Software Foundation; either version 2
;# of the License, or (at your option) any later version.
;#
;# To help us fund GROMACS development, we humbly ask that you cite
;# the research papers on the package. Check out http://www.gromacs.org
;# 
;# And Hey:
;# Gnomes, ROck Monsters And Chili Sauce
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
;# the following code is seen by GNU as, but NASM doesn't see it, so 
;# the code inside is read by NASM but not gcc.

; .if 0    # block below only read by NASM
%define .section	section
%define .long		dd
%define .align		align
%define .globl		global
;# NASM only wants 'dword', not 'dword ptr'.
%define ptr
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as

.section .text


 

.globl nb_kernel100_ia32_sse
.globl _nb_kernel100_ia32_sse
nb_kernel100_ia32_sse:	
_nb_kernel100_ia32_sse:	
.equiv          nb100_p_nri,            8
.equiv          nb100_iinr,             12
.equiv          nb100_jindex,           16
.equiv          nb100_jjnr,             20
.equiv          nb100_shift,            24
.equiv          nb100_shiftvec,         28
.equiv          nb100_fshift,           32
.equiv          nb100_gid,              36
.equiv          nb100_pos,              40
.equiv          nb100_faction,          44
.equiv          nb100_charge,           48
.equiv          nb100_p_facel,          52
.equiv          nb100_p_krf,            56
.equiv          nb100_p_crf,            60
.equiv          nb100_Vc,               64
.equiv          nb100_type,             68
.equiv          nb100_p_ntype,          72
.equiv          nb100_vdwparam,         76
.equiv          nb100_Vvdw,             80
.equiv          nb100_p_tabscale,       84
.equiv          nb100_VFtab,            88
.equiv          nb100_invsqrta,         92
.equiv          nb100_dvda,             96
.equiv          nb100_p_gbtabscale,     100
.equiv          nb100_GBtab,            104
.equiv          nb100_p_nthreads,       108
.equiv          nb100_count,            112
.equiv          nb100_mtx,              116
.equiv          nb100_outeriter,        120
.equiv          nb100_inneriter,        124
.equiv          nb100_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb100_ix,               0
.equiv          nb100_iy,               16
.equiv          nb100_iz,               32
.equiv          nb100_iq,               48
.equiv          nb100_dx,               64
.equiv          nb100_dy,               80
.equiv          nb100_dz,               96
.equiv          nb100_vctot,            112
.equiv          nb100_fix,              128
.equiv          nb100_fiy,              144
.equiv          nb100_fiz,              160
.equiv          nb100_half,             176
.equiv          nb100_three,            192
.equiv          nb100_is3,              208
.equiv          nb100_ii3,              212
.equiv          nb100_innerjjnr,        216
.equiv          nb100_innerk,           220
.equiv          nb100_n,                224
.equiv          nb100_nn1,              228
.equiv          nb100_nri,              232
.equiv          nb100_facel,            236
.equiv          nb100_nouter,           240
.equiv          nb100_ninner,           244
.equiv          nb100_salign,           248
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 252		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb100_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb100_p_nri]
	mov esi, [ebp + nb100_p_facel]
	mov ecx, [ecx]
	mov esi, [esi]
	mov [esp + nb100_nri], ecx
	mov [esp + nb100_facel], esi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb100_nouter], eax
	mov [esp + nb100_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [esp + nb100_half], eax
	movss xmm1, [esp + nb100_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [esp + nb100_half],  xmm1
	movaps [esp + nb100_three],  xmm3

.nb100_threadloop:
        mov   esi, [ebp + nb100_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb100_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                            ;# ebx=nn1=nn0+10
        lock 
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb100_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb100_n], eax
        mov [esp + nb100_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb100_outerstart
        jmp .nb100_end

.nb100_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb100_nouter]
	mov [esp + nb100_nouter], ebx

.nb100_outer:
	mov   eax, [ebp + nb100_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb100_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb100_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb100_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb100_charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [esp + nb100_facel]
	shufps xmm3, xmm3, 0
	
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb100_pos]    ;# eax = base of pos[]  

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + nb100_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [esp + nb100_ix], xmm0
	movaps [esp + nb100_iy], xmm1
	movaps [esp + nb100_iz], xmm2

	mov   [esp + nb100_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb100_vctot], xmm4
	movaps [esp + nb100_fix], xmm4
	movaps [esp + nb100_fiy], xmm4
	movaps [esp + nb100_fiz], xmm4
	
	mov   eax, [ebp + nb100_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb100_pos]
	mov   edi, [ebp + nb100_faction]	
	mov   eax, [ebp + nb100_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb100_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb100_ninner]
	mov   [esp + nb100_ninner], ecx
	add   edx, 0
	mov   [esp + nb100_innerk], edx    ;# number of innerloop atoms 
	jge   .nb100_unroll_loop
	jmp   .nb100_finish_inner
.nb100_unroll_loop:	
	;# quad-unrolled innerloop here 
	mov   edx, [esp + nb100_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	add dword ptr [esp + nb100_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb100_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm5, [esp + nb100_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# constant 10001000	      
	mov esi, [ebp + nb100_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm5
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	

	movlps xmm4, [esi + eax*4]	;# x1 y1 - - 
	movlps xmm5, [esi + ecx*4]	;# x3 y3 - - 
	movss xmm2, [esi + eax*4 + 8]	;# z1 -  - - 
	movss xmm6, [esi + ecx*4 + 8]   ;# z3 -  - - 

	movhps xmm4, [esi + ebx*4]	;# x1 y1 x2 y2 
	movhps xmm5, [esi + edx*4]	;# x3 y3 x4 y4 

	movss xmm0, [esi + ebx*4 + 8]	;# z2 - - - 
	movss xmm1, [esi + edx*4 + 8]	;# z4 - - - 

	shufps xmm2, xmm0, 0		;# z1 z1 z2 z2 
	shufps xmm6, xmm1, 0		;# z3 z3 z4 z4 
	
	movaps xmm0, xmm4		;# x1 y1 x2 y2 	
	movaps xmm1, xmm4		;# x1 y1 x2 y2 

	shufps xmm2, xmm6, 136  ;# constant 10001000	;# z1 z2 z3 z4 
	
	shufps xmm0, xmm5, 136  ;# constant 10001000	;# x1 x2 x3 x4 
	shufps xmm1, xmm5, 221  ;# constant 11011101	;# y1 y2 y3 y4 		

	mov    edi, [ebp + nb100_faction]

	;# move nb100_ix-iz to xmm4-xmm6 
	movaps xmm4, [esp + nb100_ix]
	movaps xmm5, [esp + nb100_iy]
	movaps xmm6, [esp + nb100_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb100_dx], xmm4
	movaps [esp + nb100_dy], xmm5
	movaps [esp + nb100_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb100_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb100_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 

	movaps xmm5, [esp + nb100_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm4, xmm3	;# xmm4=fscal 
	addps  xmm5, xmm3

	movaps xmm0, [esp + nb100_dx]
	movaps xmm1, [esp + nb100_dy]
	movaps xmm2, [esp + nb100_dz]

	movaps [esp + nb100_vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb100_fix]
	movaps xmm4, [esp + nb100_fiy]
	movaps xmm5, [esp + nb100_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb100_fix], xmm3
	movaps [esp + nb100_fiy], xmm4
	movaps [esp + nb100_fiz], xmm5
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm4, [edi + eax*4]
	movlps xmm6, [edi + ecx*4]
	movhps xmm4, [edi + ebx*4]
	movhps xmm6, [edi + edx*4]

	movaps xmm3, xmm4
	shufps xmm3, xmm6, 136  ;# constant 10001000
	shufps xmm4, xmm6, 221  ;# constant 11011101			      

	;# now xmm3-xmm5 contains fjx, fjy, fjz 
	subps  xmm3, xmm0
	subps  xmm4, xmm1
	
	;# unpack them back so we can store them - first x & y in xmm3/xmm4 

	movaps xmm6, xmm3
	unpcklps xmm6, xmm4
	unpckhps xmm3, xmm4	
	;# xmm6(l)=x & y for j1, (h) for j2 
	;# xmm3(l)=x & y for j3, (h) for j4 
	movlps [edi + eax*4], xmm6
	movlps [edi + ecx*4], xmm3
	
	movhps [edi + ebx*4], xmm6
	movhps [edi + edx*4], xmm3

	;# and the z forces 
	movss  xmm4, [edi + eax*4 + 8]
	movss  xmm5, [edi + ebx*4 + 8]
	movss  xmm6, [edi + ecx*4 + 8]
	movss  xmm7, [edi + edx*4 + 8]
	subss  xmm4, xmm2
	shufps xmm2, xmm2, 229  ;# constant 11100101
	subss  xmm5, xmm2
	shufps xmm2, xmm2, 234  ;# constant 11101010
	subss  xmm6, xmm2
	shufps xmm2, xmm2, 255  ;# constant 11111111
	subss  xmm7, xmm2
	movss  [edi + eax*4 + 8], xmm4
	movss  [edi + ebx*4 + 8], xmm5
	movss  [edi + ecx*4 + 8], xmm6
	movss  [edi + edx*4 + 8], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb100_innerk],  4
	jl    .nb100_finish_inner
	jmp   .nb100_unroll_loop
.nb100_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [esp + nb100_innerk],  4
	mov   edx, [esp + nb100_innerk]
	and   edx, 2
	jnz   .nb100_dopair
	jmp   .nb100_checksingle
.nb100_dopair:	
	mov esi, [ebp + nb100_charge]
	mov edi, [ebp + nb100_pos]
    	mov   ecx, [esp + nb100_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add dword ptr [esp + nb100_innerjjnr],  8

	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 0 
	shufps xmm3, xmm3, 8 ;# constant 00001000 ;# xmm3(0,1) has the charges 

	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + nb100_iq]
	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	mov    edi, [ebp + nb100_faction]
	;# move nb100_ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb100_ix]
	movaps xmm5, [esp + nb100_iy]
	movaps xmm6, [esp + nb100_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb100_dx], xmm4
	movaps [esp + nb100_dy], xmm5
	movaps [esp + nb100_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb100_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb100_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 

	movaps xmm5, [esp + nb100_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm4, xmm3	;# xmm4=fscal 
	addps  xmm5, xmm3

	movaps xmm0, [esp + nb100_dx]
	movaps xmm1, [esp + nb100_dy]
	movaps xmm2, [esp + nb100_dz]

	movaps [esp + nb100_vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb100_fix]
	movaps xmm4, [esp + nb100_fiy]
	movaps xmm5, [esp + nb100_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb100_fix], xmm3
	movaps [esp + nb100_fiy], xmm4
	movaps [esp + nb100_fiz], xmm5
	;# update the fj's  
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5

	shufps  xmm0, xmm0, 225  ;# constant 11100001
	shufps  xmm1, xmm1, 225  ;# constant 11100001
	shufps  xmm2, xmm2, 225  ;# constant 11100001
	
	movss   xmm3, [edi + ebx*4]
	movss   xmm4, [edi + ebx*4 + 4]
	movss   xmm5, [edi + ebx*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movss   [edi + ebx*4], xmm3
	movss   [edi + ebx*4 + 4], xmm4
	movss   [edi + ebx*4 + 8], xmm5	
.nb100_checksingle:				
	mov   edx, [esp + nb100_innerk]
	and   edx, 1
	jnz    .nb100_dosingle
	jmp    .nb100_updateouterdata
.nb100_dosingle:			
	mov esi, [ebp + nb100_charge]
	mov edi, [ebp + nb100_pos]
	mov   ecx, [esp + nb100_innerjjnr]
	mov   eax, [ecx]	
	movss xmm3, [esi + eax*4]	;# xmm3(0) has the charge 	
	
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + nb100_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb100_ix]
	movaps xmm5, [esp + nb100_iy]
	movaps xmm6, [esp + nb100_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb100_dx], xmm4
	movaps [esp + nb100_dy], xmm5
	movaps [esp + nb100_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb100_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb100_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	mov    edi, [ebp + nb100_faction]
	movaps xmm5, [esp + nb100_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm4, xmm3	;# xmm4=fscal 
	addss  xmm5, xmm3

	movaps xmm0, [esp + nb100_dx]
	movaps xmm1, [esp + nb100_dy]
	movaps xmm2, [esp + nb100_dz]

	movaps [esp + nb100_vctot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb100_fix]
	movaps xmm4, [esp + nb100_fiy]
	movaps xmm5, [esp + nb100_fiz]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movaps [esp + nb100_fix], xmm3
	movaps [esp + nb100_fiy], xmm4
	movaps [esp + nb100_fiz], xmm5
	;# update fj 	
	movss   xmm3, [edi + eax*4]
	movss   xmm4, [edi + eax*4 + 4]
	movss   xmm5, [edi + eax*4 + 8]
	subss   xmm3, xmm0
	subss   xmm4, xmm1
	subss   xmm5, xmm2
	movss   [edi + eax*4], xmm3
	movss   [edi + eax*4 + 4], xmm4
	movss   [edi + eax*4 + 8], xmm5	
.nb100_updateouterdata:
	mov   ecx, [esp + nb100_ii3]
	mov   edi, [ebp + nb100_faction]
	mov   esi, [ebp + nb100_fshift]
	mov   edx, [esp + nb100_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb100_fix]
	movaps xmm1, [esp + nb100_fiy]
	movaps xmm2, [esp + nb100_fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addps  xmm0, xmm3
	addps  xmm1, xmm4
	addps  xmm2, xmm5 ;# sum is in 1/2 in xmm0-xmm2 

	movaps xmm3, xmm0	
	movaps xmm4, xmm1	
	movaps xmm5, xmm2	

	shufps xmm3, xmm3, 1
	shufps xmm4, xmm4, 1
	shufps xmm5, xmm5, 1
	addss  xmm0, xmm3
	addss  xmm1, xmm4
	addss  xmm2, xmm5	;# xmm0-xmm2 has single force in pos0 

	;# increment i force 
	movss  xmm3, [edi + ecx*4]
	movss  xmm4, [edi + ecx*4 + 4]
	movss  xmm5, [edi + ecx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [edi + ecx*4],     xmm3
	movss  [edi + ecx*4 + 4], xmm4
	movss  [edi + ecx*4 + 8], xmm5

	;# increment fshift force  
	movss  xmm3, [esi + edx*4]
	movss  xmm4, [esi + edx*4 + 4]
	movss  xmm5, [esi + edx*4 + 8]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [esi + edx*4],     xmm3
	movss  [esi + edx*4 + 4], xmm4
	movss  [esi + edx*4 + 8], xmm5

	;# get n from stack
	mov esi, [esp + nb100_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb100_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb100_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb100_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb100_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb100_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb100_n], esi
        jmp .nb100_outer
.nb100_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb100_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb100_end
        ;# non-zero, do one more workunit
        jmp   .nb100_threadloop
.nb100_end:
	emms

	mov eax, [esp + nb100_nouter]
	mov ebx, [esp + nb100_ninner]
	mov ecx, [ebp + nb100_outeriter]
	mov edx, [ebp + nb100_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb100_salign]
	add esp, eax
	add esp, 252
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




.globl nb_kernel100nf_ia32_sse
.globl _nb_kernel100nf_ia32_sse
nb_kernel100nf_ia32_sse:	
_nb_kernel100nf_ia32_sse:	
.equiv          nb100nf_p_nri,          8
.equiv          nb100nf_iinr,           12
.equiv          nb100nf_jindex,         16
.equiv          nb100nf_jjnr,           20
.equiv          nb100nf_shift,          24
.equiv          nb100nf_shiftvec,       28
.equiv          nb100nf_fshift,         32
.equiv          nb100nf_gid,            36
.equiv          nb100nf_pos,            40
.equiv          nb100nf_faction,        44
.equiv          nb100nf_charge,         48
.equiv          nb100nf_p_facel,        52
.equiv          nb100nf_p_krf,          56
.equiv          nb100nf_p_crf,          60
.equiv          nb100nf_Vc,             64
.equiv          nb100nf_type,           68
.equiv          nb100nf_p_ntype,        72
.equiv          nb100nf_vdwparam,       76
.equiv          nb100nf_Vvdw,           80
.equiv          nb100nf_p_tabscale,     84
.equiv          nb100nf_VFtab,          88
.equiv          nb100nf_invsqrta,       92
.equiv          nb100nf_dvda,           96
.equiv          nb100nf_p_gbtabscale,   100
.equiv          nb100nf_GBtab,          104
.equiv          nb100nf_p_nthreads,     108
.equiv          nb100nf_count,          112
.equiv          nb100nf_mtx,            116
.equiv          nb100nf_outeriter,      120
.equiv          nb100nf_inneriter,      124
.equiv          nb100nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb100nf_ix,             0
.equiv          nb100nf_iy,             16
.equiv          nb100nf_iz,             32
.equiv          nb100nf_iq,             48
.equiv          nb100nf_vctot,          64
.equiv          nb100nf_half,           80
.equiv          nb100nf_three,          96
.equiv          nb100nf_is3,            112
.equiv          nb100nf_ii3,            116
.equiv          nb100nf_innerjjnr,      120
.equiv          nb100nf_innerk,         124
.equiv          nb100nf_n,              128
.equiv          nb100nf_nn1,            132
.equiv          nb100nf_nri,            136
.equiv          nb100nf_facel,          140
.equiv          nb100nf_nouter,         144
.equiv          nb100nf_ninner,         148
.equiv          nb100nf_salign,         152
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx 
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 156		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb100nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb100nf_p_nri]
	mov esi, [ebp + nb100nf_p_facel]
	mov ecx, [ecx]
	mov esi, [esi]
	mov [esp + nb100nf_nri], ecx
	mov [esp + nb100nf_facel], esi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb100nf_nouter], eax
	mov [esp + nb100nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [esp + nb100nf_half], eax
	movss xmm1, [esp + nb100nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [esp + nb100nf_half],  xmm1
	movaps [esp + nb100nf_three],  xmm3

.nb100nf_threadloop:
        mov   esi, [ebp + nb100nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb100nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                            ;# ebx=nn1=nn0+10
        lock 
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb100nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb100nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb100nf_n], eax
        mov [esp + nb100nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb100nf_outerstart
        jmp .nb100nf_end

.nb100nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb100nf_nouter]
	mov [esp + nb100nf_nouter], ebx

.nb100nf_outer:
	mov   eax, [ebp + nb100nf_shift]      	;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]			;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb100nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb100nf_shiftvec]   	;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb100nf_iinr]       	;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    		;# ebx =ii 

	mov   edx, [ebp + nb100nf_charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [esp + nb100nf_facel]
	shufps xmm3, xmm3, 0
	
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb100nf_pos]    ;# eax = base of pos[]  

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + nb100nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [esp + nb100nf_ix], xmm0
	movaps [esp + nb100nf_iy], xmm1
	movaps [esp + nb100nf_iz], xmm2

	mov   [esp + nb100nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb100nf_vctot], xmm4

	mov   eax, [ebp + nb100nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb100nf_pos]

	mov   eax, [ebp + nb100nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb100nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb100nf_ninner]
	mov   [esp + nb100nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb100nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb100nf_unroll_loop
	jmp   .nb100nf_finish_inner
.nb100nf_unroll_loop:	
	;# quad-unrolled innerloop here 
	mov   edx, [esp + nb100nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	add dword ptr [esp + nb100nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb100nf_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm5, [esp + nb100nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# constant 10001000	      
	mov esi, [ebp + nb100nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm5
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	

	movlps xmm4, [esi + eax*4]	;# x1 y1 - - 
	movlps xmm5, [esi + ecx*4]	;# x3 y3 - - 
	movss xmm2, [esi + eax*4 + 8]	;# z1 -  - - 
	movss xmm6, [esi + ecx*4 + 8]   ;# z3 -  - - 

	movhps xmm4, [esi + ebx*4]	;# x1 y1 x2 y2 
	movhps xmm5, [esi + edx*4]	;# x3 y3 x4 y4 

	movss xmm0, [esi + ebx*4 + 8]	;# z2 - - - 
	movss xmm1, [esi + edx*4 + 8]	;# z4 - - - 

	shufps xmm2, xmm0, 0		;# z1 z1 z2 z2 
	shufps xmm6, xmm1, 0		;# z3 z3 z4 z4 
	
	movaps xmm0, xmm4		;# x1 y1 x2 y2 	
	movaps xmm1, xmm4		;# x1 y1 x2 y2 

	shufps xmm2, xmm6, 136  ;# constant 10001000	;# z1 z2 z3 z4 
	
	shufps xmm0, xmm5, 136  ;# constant 10001000	;# x1 x2 x3 x4 
	shufps xmm1, xmm5, 221  ;# constant 11011101	;# y1 y2 y3 y4 		

	;# move nb100nf_ix-iz to xmm4-xmm6 
	movaps xmm4, [esp + nb100nf_ix]
	movaps xmm5, [esp + nb100nf_iy]
	movaps xmm6, [esp + nb100nf_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb100nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb100nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	movaps xmm5, [esp + nb100nf_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	addps  xmm5, xmm3
	movaps [esp + nb100nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb100nf_innerk],  4
	jl    .nb100nf_finish_inner
	jmp   .nb100nf_unroll_loop
.nb100nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [esp + nb100nf_innerk],  4
	mov   edx, [esp + nb100nf_innerk]
	and   edx, 2
	jnz   .nb100nf_dopair
	jmp   .nb100nf_checksingle
.nb100nf_dopair:	
	mov esi, [ebp + nb100nf_charge]
	mov edi, [ebp + nb100nf_pos]
    	mov   ecx, [esp + nb100nf_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add dword ptr [esp + nb100nf_innerjjnr],  8

	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 0 
	shufps xmm3, xmm3, 8 ;# constant 00001000 ;# xmm3(0,1) has the charges 

	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + nb100nf_iq]
	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	;# move nb100nf_ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb100nf_ix]
	movaps xmm5, [esp + nb100nf_iy]
	movaps xmm6, [esp + nb100nf_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb100nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb100nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 

	movaps xmm5, [esp + nb100nf_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	addps  xmm5, xmm3
	movaps [esp + nb100nf_vctot], xmm5
.nb100nf_checksingle:				
	mov   edx, [esp + nb100nf_innerk]
	and   edx, 1
	jnz    .nb100nf_dosingle
	jmp    .nb100nf_updateouterdata
.nb100nf_dosingle:			
	mov esi, [ebp + nb100nf_charge]
	mov edi, [ebp + nb100nf_pos]
	mov   ecx, [esp + nb100nf_innerjjnr]
	mov   eax, [ecx]	
	movss xmm3, [esi + eax*4]	;# xmm3(0) has the charge 	
	
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + nb100nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb100nf_ix]
	movaps xmm5, [esp + nb100nf_iy]
	movaps xmm6, [esp + nb100nf_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb100nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb100nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm5, [esp + nb100nf_vctot]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	addss  xmm5, xmm3
	movaps [esp + nb100nf_vctot], xmm5

.nb100nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb100nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb100nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb100nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb100nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb100nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb100nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb100nf_n], esi
        jmp .nb100nf_outer
.nb100nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb100nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb100nf_end
        ;# non-zero, do one more workunit
        jmp   .nb100nf_threadloop
.nb100nf_end:
	emms

	mov eax, [esp + nb100nf_nouter]
	mov ebx, [esp + nb100nf_ninner]
	mov ecx, [ebp + nb100nf_outeriter]
	mov edx, [ebp + nb100nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb100nf_salign]
	add esp, eax
	add esp, 156
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
