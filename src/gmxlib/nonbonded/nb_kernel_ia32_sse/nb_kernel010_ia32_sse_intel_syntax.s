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

;# nb010 - forces are calculated
.globl nb_kernel010_ia32_sse
.globl _nb_kernel010_ia32_sse
nb_kernel010_ia32_sse:	
_nb_kernel010_ia32_sse:	
.equiv          nb010_p_nri,            8
.equiv          nb010_iinr,             12
.equiv          nb010_jindex,           16
.equiv          nb010_jjnr,             20
.equiv          nb010_shift,            24
.equiv          nb010_shiftvec,         28
.equiv          nb010_fshift,           32
.equiv          nb010_gid,              36
.equiv          nb010_pos,              40
.equiv          nb010_faction,          44
.equiv          nb010_charge,           48
.equiv          nb010_p_facel,          52
.equiv          nb010_p_krf,            56
.equiv          nb010_p_crf,            60
.equiv          nb010_Vc,               64
.equiv          nb010_type,             68
.equiv          nb010_p_ntype,          72
.equiv          nb010_vdwparam,         76
.equiv          nb010_Vvdw,             80
.equiv          nb010_p_tabscale,       84
.equiv          nb010_VFtab,            88
.equiv          nb010_invsqrta,         92
.equiv          nb010_dvda,             96
.equiv          nb010_p_gbtabscale,     100
.equiv          nb010_GBtab,            104
.equiv          nb010_p_nthreads,       108
.equiv          nb010_count,            112
.equiv          nb010_mtx,              116
.equiv          nb010_outeriter,        120
.equiv          nb010_inneriter,        124
.equiv          nb010_work,             128
        ;# The mutex (last arg) is not used in assembly.
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb010_ix,               0
.equiv          nb010_iy,               16
.equiv          nb010_iz,               32
.equiv          nb010_dx,               48
.equiv          nb010_dy,               64
.equiv          nb010_dz,               80
.equiv          nb010_two,              96
.equiv          nb010_c6,               112
.equiv          nb010_c12,              128
.equiv          nb010_six,              144
.equiv          nb010_twelve,           160
.equiv          nb010_Vvdwtot,          176
.equiv          nb010_fix,              192
.equiv          nb010_fiy,              208
.equiv          nb010_fiz,              224
.equiv          nb010_half,             240
.equiv          nb010_three,            256
.equiv          nb010_is3,              272
.equiv          nb010_ii3,              276
.equiv          nb010_ntia,             280
.equiv          nb010_innerjjnr,        284
.equiv          nb010_innerk,           288
.equiv          nb010_n,                292
.equiv          nb010_nn1,              296
.equiv          nb010_nri,              300
.equiv          nb010_ntype,            304
.equiv          nb010_nouter,           308
.equiv          nb010_ninner,           312
.equiv          nb010_salign,           316
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 320		;# local stack space 
	mov  eax, esp
	and  eax, 0xf   	;# constant 16-byte align bottom of stack 
	sub esp, eax
	mov [esp + nb010_salign], eax

	emms
	
	;# Move args passed by reference to stack
	mov ecx, [ebp + nb010_p_nri]
	mov edi, [ebp + nb010_p_ntype]
	mov ecx, [ecx]
	mov edi, [edi]
	mov [esp + nb010_nri], ecx
	mov [esp + nb010_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb010_nouter], eax
	mov [esp + nb010_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x40000000     ;# constant 2.0 in IEEE (hex)
	mov [esp + nb010_two], eax
	movss xmm1, [esp + nb010_two]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1
	addps  xmm2, xmm1       ;# 4.0
	addps  xmm2, xmm1       ;# 6.0
	movaps xmm3, xmm2
	addps  xmm3, xmm3       ;# constant 12.0
	movaps [esp + nb010_two], xmm1
	movaps [esp + nb010_six],  xmm2
	movaps [esp + nb010_twelve], xmm3

.nb010_threadloop:
        mov   esi, [ebp + nb010_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb010_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                            ;# ebx=nn1=nn0+10
        lock 
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb010_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb010_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb010_n], eax
        mov [esp + nb010_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb010_outerstart
        jmp .nb010_end

.nb010_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb010_nouter]
	mov [esp + nb010_nouter], ebx

.nb010_outer:
        mov   eax, [ebp + nb010_shift]      	;# eax = base of shift[] 
        mov   ebx, [eax + esi*4]            	;# ebx=shift[n] 

        lea   ebx, [ebx + ebx*2]    		;# ebx=3*is 
        mov   [esp + nb010_is3],ebx     	;# store is3 

        mov   eax, [ebp + nb010_shiftvec]   	;# eax = base of shiftvec[] 

        movss xmm0, [eax + ebx*4]
        movss xmm1, [eax + ebx*4 + 4]
        movss xmm2, [eax + ebx*4 + 8] 

        mov   ecx, [ebp + nb010_iinr]       	;# ecx = base of iinr[] 
        mov   ebx, [ecx + esi*4]            	;# ebx =ii 

    	mov  edx, [ebp + nb010_type] 
    	mov  edx, [edx + ebx*4]
    	imul edx, [esp + nb010_ntype]
    	shl  edx, 1
    	mov  [esp + nb010_ntia], edx

        lea   ebx, [ebx + ebx*2]        	;# ebx = 3*ii=ii3 
        mov   eax, [ebp + nb010_pos]    	;# eax = base of pos[]  

        addss xmm0, [eax + ebx*4]
        addss xmm1, [eax + ebx*4 + 4]
        addss xmm2, [eax + ebx*4 + 8]

        shufps xmm0, xmm0, 0
        shufps xmm1, xmm1, 0
        shufps xmm2, xmm2, 0

        movaps [esp + nb010_ix], xmm0
        movaps [esp + nb010_iy], xmm1
        movaps [esp + nb010_iz], xmm2

        mov   [esp + nb010_ii3], ebx

        ;# clear Vvdwtot and i forces 
        xorps xmm4, xmm4
        movaps [esp + nb010_Vvdwtot], xmm4
        movaps [esp + nb010_fix], xmm4
        movaps [esp + nb010_fiy], xmm4
        movaps [esp + nb010_fiz], xmm4

        mov   eax, [ebp + nb010_jindex]
        mov   ecx, [eax  + esi*4]    		;# jindex[n] 
        mov   edx, [eax + esi*4 + 4]         	;# jindex[n+1] 
        sub   edx, ecx               		;# number of innerloop atoms 

        mov   esi, [ebp + nb010_pos]
        mov   edi, [ebp + nb010_faction]
        mov   eax, [ebp + nb010_jjnr]
        shl   ecx, 2
        add   eax, ecx
        mov   [esp + nb010_innerjjnr], eax      ;# pointer to jjnr[nj0] 
	mov   ecx, edx
        sub   edx,  4
	add   ecx, [esp + nb010_ninner]
	mov   [esp + nb010_ninner], ecx
	add   edx, 0
        mov   [esp + nb010_innerk], edx         ;# number of innerloop atoms 

        jge   .nb010_unroll_loop
        jmp   .nb010_finish_inner
.nb010_unroll_loop:
	;# quad-unrolled innerloop starts here 
	mov   edx, [esp + nb010_innerjjnr]	;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]			;# eax-edx=jnr1-4 
	;# advance pointer (unrolled 4) 
	add   dword ptr [esp + nb010_innerjjnr],  16 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + nb010_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb010_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb010_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + nb010_c6], xmm4
	movaps [esp + nb010_c12], xmm6
	
	mov esi, [ebp + nb010_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ix-iz to xmm4-xmm6 
	movaps xmm4, [esp + nb010_ix]
	movaps xmm5, [esp + nb010_iy]
	movaps xmm6, [esp + nb010_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb010_dx], xmm4 
	movaps [esp + nb010_dy], xmm5
	movaps [esp + nb010_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6

	;# rsq in xmm4 
	rcpps xmm5, xmm4
	;# constant 1/x lookup seed in xmm5 
	movaps xmm0, [esp + nb010_two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;# xmm0=rinvsq
	
	movaps xmm4, xmm0

	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm1, [esp + nb010_c6]
	mulps  xmm2, [esp + nb010_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb010_Vvdwtot]
	mulps  xmm1, [esp + nb010_six]
	mulps  xmm2, [esp + nb010_twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	;# xmm4=total fscal 

	movaps xmm0, [esp + nb010_dx]
	movaps xmm1, [esp + nb010_dy]
	movaps xmm2, [esp + nb010_dz]

	movaps [esp + nb010_Vvdwtot], xmm5

	mov    edi, [ebp + nb010_faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb010_fix]
	movaps xmm4, [esp + nb010_fiy]
	movaps xmm5, [esp + nb010_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb010_fix], xmm3
	movaps [esp + nb010_fiy], xmm4
	movaps [esp + nb010_fiz], xmm5
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
	sub   dword ptr [esp + nb010_innerk],  4
	jl    .nb010_finish_inner
	jmp   .nb010_unroll_loop
.nb010_finish_inner:
	;# check if at least two particles remain 
	add   dword ptr [esp + nb010_innerk],  4
	mov   edx, [esp + nb010_innerk]
	and   edx, 2
	jnz   .nb010_dopair
	jmp   .nb010_checksingle
.nb010_dopair:	
	mov   ecx, [esp + nb010_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   dword ptr [esp + nb010_innerjjnr],  8

	mov esi, [ebp + nb010_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + nb010_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb010_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + nb010_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# constant 00001000 	
	shufps xmm6, xmm6, 13 ;# constant 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + nb010_c6], xmm4
	movaps [esp + nb010_c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	


	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	mov    edi, [ebp + nb010_faction]
	;# move nb010_ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb010_ix]
	movaps xmm5, [esp + nb010_iy]
	movaps xmm6, [esp + nb010_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb010_dx], xmm4
	movaps [esp + nb010_dy], xmm5
	movaps [esp + nb010_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 


	rcpps xmm5, xmm4
	;# constant 1/x lookup seed in xmm5 
	movaps xmm0, [esp + nb010_two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;# xmm0=rinvsq 
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm1, [esp + nb010_c6]
	mulps  xmm2, [esp + nb010_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb010_Vvdwtot]
	mulps  xmm1, [esp + nb010_six]
	mulps  xmm2, [esp + nb010_twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	;# xmm4=total fscal 

	movaps xmm0, [esp + nb010_dx]
	movaps xmm1, [esp + nb010_dy]
	movaps xmm2, [esp + nb010_dz]

	movaps [esp + nb010_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb010_fix]
	movaps xmm4, [esp + nb010_fiy]
	movaps xmm5, [esp + nb010_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb010_fix], xmm3
	movaps [esp + nb010_fiy], xmm4
	movaps [esp + nb010_fiz], xmm5
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

.nb010_checksingle:				
	mov   edx, [esp + nb010_innerk]
	and   edx, 1
	jnz    .nb010_dosingle
	jmp    .nb010_updateouterdata
.nb010_dosingle:
	mov edi, [ebp + nb010_pos]
	mov   ecx, [esp + nb010_innerjjnr]
	mov   eax, [ecx]		

	mov esi, [ebp + nb010_type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + nb010_vdwparam]
	shl ecx, 1
	add ecx, [esp + nb010_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# constant 11111100
	shufps xmm6, xmm6, 253  ;# constant 11111101	
			
	movaps [esp + nb010_c6], xmm4
	movaps [esp + nb010_c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb010_ix]
	movaps xmm5, [esp + nb010_iy]
	movaps xmm6, [esp + nb010_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb010_dx], xmm4
	movaps [esp + nb010_dy], xmm5
	movaps [esp + nb010_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	rcpps xmm5, xmm4
	;# constant 1/x lookup seed in xmm5 
	movaps xmm0, [esp + nb010_two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;# xmm0=rinvsq 
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm1, [esp + nb010_c6]
	mulps  xmm2, [esp + nb010_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addss  xmm5, [esp + nb010_Vvdwtot]
	mulps  xmm1, [esp + nb010_six]
	mulps  xmm2, [esp + nb010_twelve]
	subps  xmm2, xmm1
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	
	mov    edi, [ebp + nb010_faction]

	movaps xmm0, [esp + nb010_dx]
	movaps xmm1, [esp + nb010_dy]
	movaps xmm2, [esp + nb010_dz]

	movss [esp + nb010_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb010_fix]
	movaps xmm4, [esp + nb010_fiy]
	movaps xmm5, [esp + nb010_fiz]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movaps [esp + nb010_fix], xmm3
	movaps [esp + nb010_fiy], xmm4
	movaps [esp + nb010_fiz], xmm5
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
	
.nb010_updateouterdata:
	mov   ecx, [esp + nb010_ii3]
	mov   edi, [ebp + nb010_faction]
	mov   esi, [ebp + nb010_fshift]
	mov   edx, [esp + nb010_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb010_fix]
	movaps xmm1, [esp + nb010_fiy]
	movaps xmm2, [esp + nb010_fiz]

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
	mov esi, [esp + nb010_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb010_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

         ;# accumulate total lj energy and update it 
        movaps xmm7, [esp + nb010_Vvdwtot]
        ;# accumulate 
        movhlps xmm6, xmm7
        addps  xmm7, xmm6       ;# pos 0-1 in xmm7 have the sum now 
        movaps xmm6, xmm7
        shufps xmm6, xmm6, 1
        addss  xmm7, xmm6

        ;# add earlier value from mem 
        mov   eax, [ebp + nb010_Vvdw]
        addss xmm7, [eax + edx*4] 
        ;# move back to mem 
        movss [eax + edx*4], xmm7 

        ;# finish if last 
        mov ecx, [esp + nb010_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb010_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb010_n], esi
        jmp .nb010_outer
.nb010_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb010_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb010_end
        ;# non-zero, do one more workunit
        jmp   .nb010_threadloop
.nb010_end:
	emms

	mov eax, [esp + nb010_nouter]
	mov ebx, [esp + nb010_ninner]
	mov ecx, [ebp + nb010_outeriter]
	mov edx, [ebp + nb010_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb010_salign]
	add esp, eax 		;# account for stack alignment 
	add esp, 320 		
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret


	
.globl nb_kernel010nf_ia32_sse
.globl _nb_kernel010nf_ia32_sse
nb_kernel010nf_ia32_sse:	
_nb_kernel010nf_ia32_sse:	
.equiv          nb010nf_p_nri,          8
.equiv          nb010nf_iinr,           12
.equiv          nb010nf_jindex,         16
.equiv          nb010nf_jjnr,           20
.equiv          nb010nf_shift,          24
.equiv          nb010nf_shiftvec,       28
.equiv          nb010nf_fshift,         32
.equiv          nb010nf_gid,            36
.equiv          nb010nf_pos,            40
.equiv          nb010nf_faction,        44
.equiv          nb010nf_charge,         48
.equiv          nb010nf_p_facel,        52
.equiv          nb010nf_p_krf,          56
.equiv          nb010nf_p_crf,          60
.equiv          nb010nf_Vc,             64
.equiv          nb010nf_type,           68
.equiv          nb010nf_p_ntype,        72
.equiv          nb010nf_vdwparam,       76
.equiv          nb010nf_Vvdw,           80
.equiv          nb010nf_p_tabscale,     84
.equiv          nb010nf_VFtab,          88
.equiv          nb010nf_invsqrta,       92
.equiv          nb010nf_dvda,           96
.equiv          nb010nf_p_gbtabscale,   100
.equiv          nb010nf_GBtab,          104
.equiv          nb010nf_p_nthreads,     108
.equiv          nb010nf_count,          112
.equiv          nb010nf_mtx,            116
.equiv          nb010nf_outeriter,      120
.equiv          nb010nf_inneriter,      124
.equiv          nb010nf_work,           128
        ;# The mutex (last arg) is not used in assembly.
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb010nf_ix,             0
.equiv          nb010nf_iy,             16
.equiv          nb010nf_iz,             32
.equiv          nb010nf_two,            48
.equiv          nb010nf_c6,             64
.equiv          nb010nf_c12,            80
.equiv          nb010nf_Vvdwtot,        96
.equiv          nb010nf_half,           112
.equiv          nb010nf_three,          128
.equiv          nb010nf_is3,            144
.equiv          nb010nf_ii3,            148
.equiv          nb010nf_ntia,           152
.equiv          nb010nf_innerjjnr,      156
.equiv          nb010nf_innerk,         160
.equiv          nb010nf_n,              164
.equiv          nb010nf_nn1,            168
.equiv          nb010nf_nri,            172
.equiv          nb010nf_ntype,          176
.equiv          nb010nf_nouter,         180
.equiv          nb010nf_ninner,         184
.equiv          nb010nf_salign,         188
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi 
	sub esp, 192		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb010nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb010nf_p_nri]
	mov edi, [ebp + nb010nf_p_ntype]
	mov ecx, [ecx]
	mov edi, [edi]
	mov [esp + nb010nf_nri], ecx
	mov [esp + nb010nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb010nf_nouter], eax
	mov [esp + nb010nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x40000000     ;# constant 2.0 in IEEE (hex)
	mov [esp + nb010nf_two], eax
	movss xmm1, [esp + nb010nf_two]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps [esp + nb010nf_two], xmm1

.nb010nf_threadloop:
        mov   esi, [ebp + nb010nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb010nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock 
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb010nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb010nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb010nf_n], eax
        mov [esp + nb010nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb010nf_outerstart
        jmp .nb010nf_end

.nb010nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb010nf_nouter]
	mov [esp + nb010nf_nouter], ebx

.nb010nf_outer:
        mov   eax, [ebp + nb010nf_shift]      	;# eax = base of shift[] 
        mov   ebx, [eax +esi*4]                	;# ebx=shift[n] 
        
        lea   ebx, [ebx + ebx*2]    		;# ebx=3*is 
        mov   [esp + nb010nf_is3],ebx           ;# store is3 

        mov   eax, [ebp + nb010nf_shiftvec]   	;# eax = base of shiftvec[] 

        movss xmm0, [eax + ebx*4]
        movss xmm1, [eax + ebx*4 + 4]
        movss xmm2, [eax + ebx*4 + 8] 

        mov   ecx, [ebp + nb010nf_iinr]       	;# ecx = base of iinr[] 
        mov   ebx, [ecx + esi*4]            	;# ebx =ii 

        mov   edx, [ebp + nb010nf_type] 
        mov   edx, [edx + ebx*4]
        imul  edx, [esp + nb010nf_ntype]
        shl   edx, 1
        mov   [esp + nb010nf_ntia], edx

        lea   ebx, [ebx + ebx*2]        ;# ebx = 3*ii=ii3 
        mov   eax, [ebp + nb010nf_pos]    ;# eax = base of pos[]  

        addss xmm0, [eax + ebx*4]
        addss xmm1, [eax + ebx*4 + 4]
        addss xmm2, [eax + ebx*4 + 8]

        shufps xmm0, xmm0, 0
        shufps xmm1, xmm1, 0
        shufps xmm2, xmm2, 0

	movaps [esp + nb010nf_ix], xmm0
	movaps [esp + nb010nf_iy], xmm1
	movaps [esp + nb010nf_iz], xmm2

	mov   [esp + nb010nf_ii3], ebx
	
	;# clear Vvdwtot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb010nf_Vvdwtot], xmm4
	
        mov   eax, [ebp + nb010nf_jindex]
        mov   ecx, [eax + esi*4]             ;# jindex[n] 
        mov   edx, [eax + esi*4 + 4]         ;# jindex[n+1] 
        sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb010nf_pos]
	mov   eax, [ebp + nb010nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb010nf_innerjjnr], eax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb010nf_ninner]
	mov   [esp + nb010nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb010nf_innerk], edx    	;# number of innerloop atoms 
	
	jge   .nb010nf_unroll_loop
	jmp   .nb010nf_finish_inner
.nb010nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb010nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	;# advance pointer (unrolled 4) 
	add   dword ptr [esp + nb010nf_innerjjnr],  16 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + nb010nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb010nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb010nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [esi + eax*4]
	movlps xmm7, [esi + ecx*4]
	movhps xmm6, [esi + ebx*4]
	movhps xmm7, [esi + edx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [esp + nb010nf_c6], xmm4
	movaps [esp + nb010nf_c12], xmm6
	
	mov esi, [ebp + nb010nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	mulps xmm3, xmm2
	lea   ecx, [ecx + ecx*2]     ;# replace jnr with j3 
	lea   edx, [edx + edx*2]	

	;# move four coordinates to xmm0-xmm2 	

	movlps xmm4, [esi + eax*4]
	movlps xmm5, [esi + ecx*4]
	movss xmm2, [esi + eax*4 + 8]
	movss xmm6, [esi + ecx*4 + 8]

	movhps xmm4, [esi + ebx*4]
	movhps xmm5, [esi + edx*4]

	movss xmm0, [esi + ebx*4 + 8]
	movss xmm1, [esi + edx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ix-iz to xmm4-xmm6 
	movaps xmm4, [esp + nb010nf_ix]
	movaps xmm5, [esp + nb010nf_iy]
	movaps xmm6, [esp + nb010nf_iz]

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
	rcpps xmm5, xmm4
	;# constant 1/x lookup seed in xmm5 
	movaps xmm0, [esp + nb010nf_two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;# xmm0=rinvsq 
	movaps xmm4, xmm0

	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm1, [esp + nb010nf_c6]
	mulps  xmm2, [esp + nb010nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb010nf_Vvdwtot]
	movaps [esp + nb010nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub   dword ptr [esp + nb010nf_innerk],  4
	jl    .nb010nf_finish_inner
	jmp   .nb010nf_unroll_loop
.nb010nf_finish_inner:
	;# check if at least two particles remain 
	add   dword ptr [esp + nb010nf_innerk],  4
	mov   edx, [esp + nb010nf_innerk]
	and   edx, 2
	jnz   .nb010nf_dopair
	jmp   .nb010nf_checksingle
.nb010nf_dopair:	
	mov   ecx, [esp + nb010nf_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add   dword ptr [esp + nb010nf_innerjjnr],  8

	mov esi, [ebp + nb010nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + nb010nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb010nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + nb010nf_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# constant 00001000 	
	shufps xmm6, xmm6, 13 ;# constant 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + nb010nf_c6], xmm4
	movaps [esp + nb010nf_c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	;# move nb010nf_ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb010nf_ix]
	movaps xmm5, [esp + nb010nf_iy]
	movaps xmm6, [esp + nb010nf_iz]

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


	rcpps xmm5, xmm4
	;# constant 1/x lookup seed in xmm5 
	movaps xmm0, [esp + nb010nf_two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;# xmm0=rinvsq 
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm1, [esp + nb010nf_c6]
	mulps  xmm2, [esp + nb010nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb010nf_Vvdwtot]
	movaps [esp + nb010nf_Vvdwtot], xmm5

.nb010nf_checksingle:				
	mov   edx, [esp + nb010nf_innerk]
	and   edx, 1
	jnz    .nb010nf_dosingle
	jmp    .nb010nf_updateouterdata
.nb010nf_dosingle:
	mov edi, [ebp + nb010nf_pos]
	mov   ecx, [esp + nb010nf_innerjjnr]
	mov   eax, [ecx]		

	mov esi, [ebp + nb010nf_type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + nb010nf_vdwparam]
	shl ecx, 1
	add ecx, [esp + nb010nf_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# constant 11111100
	shufps xmm6, xmm6, 253  ;# constant 11111101	
			
	movaps [esp + nb010nf_c6], xmm4
	movaps [esp + nb010nf_c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb010nf_ix]
	movaps xmm5, [esp + nb010nf_iy]
	movaps xmm6, [esp + nb010nf_iz]

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

	rcpps xmm5, xmm4
	;# constant 1/x lookup seed in xmm5 
	movaps xmm0, [esp + nb010nf_two]
	mulps xmm4, xmm5
	subps xmm0, xmm4
	mulps xmm0, xmm5	;# xmm0=rinvsq 
	movaps xmm4, xmm0
	
	movaps xmm1, xmm0
	mulps  xmm1, xmm0
	mulps  xmm1, xmm0	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm1, [esp + nb010nf_c6]
	mulps  xmm2, [esp + nb010nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addss  xmm5, [esp + nb010nf_Vvdwtot]
	movss [esp + nb010nf_Vvdwtot], xmm5
	
.nb010nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb010nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb010nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

        ;# accumulate total lj energy and update it 
        movaps xmm7, [esp + nb010nf_Vvdwtot]
        ;# accumulate 
        movhlps xmm6, xmm7
        addps  xmm7, xmm6       ;# pos 0-1 in xmm7 have the sum now 
        movaps xmm6, xmm7
        shufps xmm6, xmm6, 1
        addss  xmm7, xmm6

        ;# add earlier value from mem 
        mov   eax, [ebp + nb010nf_Vvdw]
        addss xmm7, [eax + edx*4] 
        ;# move back to mem 
        movss [eax + edx*4], xmm7 

        ;# finish if last 
        mov ecx, [esp + nb010nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb010nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb010nf_n], esi
        jmp .nb010nf_outer
.nb010nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb010nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb010nf_end
        ;# non-zero, do one more workunit
        jmp   .nb010nf_threadloop
.nb010nf_end:
	emms

	mov eax, [esp + nb010nf_nouter]
	mov ebx, [esp + nb010nf_ninner]
	mov ecx, [ebp + nb010nf_outeriter]
	mov edx, [ebp + nb010nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb010nf_salign]
	add esp, eax
	add esp, 192
	pop edi
	pop esi
	pop edx
	pop ecx
	pop ebx
	pop eax
	leave
	ret
