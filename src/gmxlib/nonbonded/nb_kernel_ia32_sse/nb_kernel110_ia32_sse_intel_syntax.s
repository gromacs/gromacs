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




.globl nb_kernel110_ia32_sse
.globl _nb_kernel110_ia32_sse
nb_kernel110_ia32_sse:	
_nb_kernel110_ia32_sse:	
.equiv          nb110_p_nri,            8
.equiv          nb110_iinr,             12
.equiv          nb110_jindex,           16
.equiv          nb110_jjnr,             20
.equiv          nb110_shift,            24
.equiv          nb110_shiftvec,         28
.equiv          nb110_fshift,           32
.equiv          nb110_gid,              36
.equiv          nb110_pos,              40
.equiv          nb110_faction,          44
.equiv          nb110_charge,           48
.equiv          nb110_p_facel,          52
.equiv          nb110_p_krf,            56
.equiv          nb110_p_crf,            60
.equiv          nb110_Vc,               64
.equiv          nb110_type,             68
.equiv          nb110_p_ntype,          72
.equiv          nb110_vdwparam,         76
.equiv          nb110_Vvdw,             80
.equiv          nb110_p_tabscale,       84
.equiv          nb110_VFtab,            88
.equiv          nb110_invsqrta,         92
.equiv          nb110_dvda,             96
.equiv          nb110_p_gbtabscale,     100
.equiv          nb110_GBtab,            104
.equiv          nb110_p_nthreads,       108
.equiv          nb110_count,            112
.equiv          nb110_mtx,              116
.equiv          nb110_outeriter,        120
.equiv          nb110_inneriter,        124
.equiv          nb110_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb110_ix,               0
.equiv          nb110_iy,               16
.equiv          nb110_iz,               32
.equiv          nb110_iq,               48
.equiv          nb110_dx,               64
.equiv          nb110_dy,               80
.equiv          nb110_dz,               96
.equiv          nb110_c6,               112
.equiv          nb110_c12,              128
.equiv          nb110_six,              144
.equiv          nb110_twelve,           160
.equiv          nb110_vctot,            176
.equiv          nb110_Vvdwtot,          192
.equiv          nb110_fix,              208
.equiv          nb110_fiy,              224
.equiv          nb110_fiz,              240
.equiv          nb110_half,             256
.equiv          nb110_three,            272
.equiv          nb110_is3,              288
.equiv          nb110_ii3,              292
.equiv          nb110_ntia,             296
.equiv          nb110_innerjjnr,        300
.equiv          nb110_innerk,           304
.equiv          nb110_n,                308
.equiv          nb110_nn1,              312
.equiv          nb110_nri,              316
.equiv          nb110_facel,            320
.equiv          nb110_ntype,            324
.equiv          nb110_nouter,           328
.equiv          nb110_ninner,           332
.equiv          nb110_salign,           336
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  340		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb110_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb110_p_nri]
	mov esi, [ebp + nb110_p_facel]
	mov edi, [ebp + nb110_p_ntype]
	mov ecx, [ecx]
	mov esi, [esi]
	mov edi, [edi]
	mov [esp + nb110_nri], ecx
	mov [esp + nb110_facel], esi
	mov [esp + nb110_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb110_nouter], eax
	mov [esp + nb110_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb110_half], eax
	movss xmm1, [esp + nb110_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# 6.0
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# constant 12.0
	movaps [esp + nb110_half],  xmm1
	movaps [esp + nb110_three],  xmm3
	movaps [esp + nb110_six],  xmm4
	movaps [esp + nb110_twelve],  xmm5

.nb110_threadloop:
        mov   esi, [ebp + nb110_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb110_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb110_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb110_n], eax
        mov [esp + nb110_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110_outerstart
        jmp .nb110_end

.nb110_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb110_nouter]
	mov [esp + nb110_nouter], ebx

.nb110_outer:
	mov   eax, [ebp + nb110_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb110_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb110_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb110_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx +esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb110_charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [esp + nb110_facel]
	shufps xmm3, xmm3, 0

    	mov   edx, [ebp + nb110_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb110_ntype]
    	shl   edx, 1
    	mov   [esp + nb110_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb110_pos]    ;# eax = base of pos[]  

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + nb110_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [esp + nb110_ix], xmm0
	movaps [esp + nb110_iy], xmm1
	movaps [esp + nb110_iz], xmm2

	mov   [esp + nb110_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb110_vctot], xmm4
	movaps [esp + nb110_Vvdwtot], xmm4
	movaps [esp + nb110_fix], xmm4
	movaps [esp + nb110_fiy], xmm4
	movaps [esp + nb110_fiz], xmm4
	
	mov   eax, [ebp + nb110_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb110_pos]
	mov   edi, [ebp + nb110_faction]	
	mov   eax, [ebp + nb110_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb110_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb110_ninner]
	mov   [esp + nb110_ninner], ecx
	add   edx, 0
	mov   [esp + nb110_innerk], edx    ;# number of innerloop atoms 
	jge   .nb110_unroll_loop
	jmp   .nb110_finish_inner
.nb110_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb110_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	add dword ptr [esp + nb110_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb110_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + nb110_iq]
	shufps xmm3, xmm6, 0
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + nb110_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb110_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb110_ntia]
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

	movaps [esp + nb110_c6], xmm4
	movaps [esp + nb110_c12], xmm6
	
	mov esi, [ebp + nb110_pos]       ;# base of pos[] 

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
	movaps xmm4, [esp + nb110_ix]
	movaps xmm5, [esp + nb110_iy]
	movaps xmm6, [esp + nb110_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb110_dx], xmm4
	movaps [esp + nb110_dy], xmm5
	movaps [esp + nb110_dz], xmm6
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
	movaps xmm1, [esp + nb110_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb110_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [esp + nb110_c6]
	mulps  xmm2, [esp + nb110_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb110_Vvdwtot]
	mulps  xmm1, [esp + nb110_six]
	mulps  xmm2, [esp + nb110_twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addps  xmm3, [esp + nb110_vctot]

	movaps xmm0, [esp + nb110_dx]
	movaps xmm1, [esp + nb110_dy]
	movaps xmm2, [esp + nb110_dz]

	movaps [esp + nb110_vctot], xmm3
	movaps [esp + nb110_Vvdwtot], xmm5

	mov    edi, [ebp + nb110_faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb110_fix]
	movaps xmm4, [esp + nb110_fiy]
	movaps xmm5, [esp + nb110_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb110_fix], xmm3
	movaps [esp + nb110_fiy], xmm4
	movaps [esp + nb110_fiz], xmm5
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
	sub dword ptr [esp + nb110_innerk],  4
	jl    .nb110_finish_inner
	jmp   .nb110_unroll_loop
.nb110_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [esp + nb110_innerk],  4
	mov   edx, [esp + nb110_innerk]
	and   edx, 2
	jnz   .nb110_dopair
	jmp   .nb110_checksingle
.nb110_dopair:	
	mov esi, [ebp + nb110_charge]

    mov   ecx, [esp + nb110_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add dword ptr [esp + nb110_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 12 ;# constant 00001100 
	shufps xmm3, xmm3, 88 ;# constant 01011000 ;# xmm3(0,1) has the charges 

	mov esi, [ebp + nb110_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + nb110_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb110_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + nb110_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# constant 00001000 	
	shufps xmm6, xmm6, 13 ;# constant 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + nb110_c6], xmm4
	movaps [esp + nb110_c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + nb110_iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	mov    edi, [ebp + nb110_faction]
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb110_ix]
	movaps xmm5, [esp + nb110_iy]
	movaps xmm6, [esp + nb110_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb110_dx], xmm4
	movaps [esp + nb110_dy], xmm5
	movaps [esp + nb110_dz], xmm6
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
	movaps xmm1, [esp + nb110_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb110_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [esp + nb110_c6]
	mulps  xmm2, [esp + nb110_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb110_Vvdwtot]
	mulps  xmm1, [esp + nb110_six]
	mulps  xmm2, [esp + nb110_twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addps  xmm3, [esp + nb110_vctot]

	movaps xmm0, [esp + nb110_dx]
	movaps xmm1, [esp + nb110_dy]
	movaps xmm2, [esp + nb110_dz]

	movaps [esp + nb110_vctot], xmm3
	movaps [esp + nb110_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb110_fix]
	movaps xmm4, [esp + nb110_fiy]
	movaps xmm5, [esp + nb110_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb110_fix], xmm3
	movaps [esp + nb110_fiy], xmm4
	movaps [esp + nb110_fiz], xmm5
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

.nb110_checksingle:				
	mov   edx, [esp + nb110_innerk]
	and   edx, 1
	jnz    .nb110_dosingle
	jmp    .nb110_updateouterdata
.nb110_dosingle:			
	mov esi, [ebp + nb110_charge]
	mov edi, [ebp + nb110_pos]
	mov   ecx, [esp + nb110_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	;# xmm3(0) has the charge 	

	mov esi, [ebp + nb110_type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + nb110_vdwparam]
	shl ecx, 1
	add ecx, [esp + nb110_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# constant 11111100	
	shufps xmm6, xmm6, 253  ;# constant 11111101	
			
	movaps [esp + nb110_c6], xmm4
	movaps [esp + nb110_c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + nb110_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb110_ix]
	movaps xmm5, [esp + nb110_iy]
	movaps xmm6, [esp + nb110_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb110_dx], xmm4
	movaps [esp + nb110_dy], xmm5
	movaps [esp + nb110_dz], xmm6
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
	movaps xmm1, [esp + nb110_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb110_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [esp + nb110_c6]
	mulps  xmm2, [esp + nb110_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addss  xmm5, [esp + nb110_Vvdwtot]
	mulps  xmm1, [esp + nb110_six]
	mulps  xmm2, [esp + nb110_twelve]
	subps  xmm2, xmm1
	addps  xmm2, xmm3
	mulps  xmm4, xmm2	;# xmm4=total fscal 
	addss  xmm3, [esp + nb110_vctot]
	
	mov    edi, [ebp + nb110_faction]

	movaps xmm0, [esp + nb110_dx]
	movaps xmm1, [esp + nb110_dy]
	movaps xmm2, [esp + nb110_dz]

	movss [esp + nb110_vctot], xmm3
	movss [esp + nb110_Vvdwtot], xmm5

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb110_fix]
	movaps xmm4, [esp + nb110_fiy]
	movaps xmm5, [esp + nb110_fiz]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movaps [esp + nb110_fix], xmm3
	movaps [esp + nb110_fiy], xmm4
	movaps [esp + nb110_fiz], xmm5
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
.nb110_updateouterdata:
	mov   ecx, [esp + nb110_ii3]
	mov   edi, [ebp + nb110_faction]
	mov   esi, [ebp + nb110_fshift]
	mov   edx, [esp + nb110_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb110_fix]
	movaps xmm1, [esp + nb110_fiy]
	movaps xmm2, [esp + nb110_fiz]

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
	mov esi, [esp + nb110_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb110_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb110_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb110_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb110_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb110_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb110_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb110_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb110_n], esi
        jmp .nb110_outer
.nb110_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb110_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb110_end
        ;# non-zero, do one more workunit
        jmp   .nb110_threadloop
.nb110_end:
	emms

	mov eax, [esp + nb110_nouter]
	mov ebx, [esp + nb110_ninner]
	mov ecx, [ebp + nb110_outeriter]
	mov edx, [ebp + nb110_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb110_salign]
	add esp, eax
	add esp,  340
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret



.globl nb_kernel110nf_ia32_sse
.globl _nb_kernel110nf_ia32_sse
nb_kernel110nf_ia32_sse:	
_nb_kernel110nf_ia32_sse:	
.equiv          nb110nf_p_nri,          8
.equiv          nb110nf_iinr,           12
.equiv          nb110nf_jindex,         16
.equiv          nb110nf_jjnr,           20
.equiv          nb110nf_shift,          24
.equiv          nb110nf_shiftvec,       28
.equiv          nb110nf_fshift,         32
.equiv          nb110nf_gid,            36
.equiv          nb110nf_pos,            40
.equiv          nb110nf_faction,        44
.equiv          nb110nf_charge,         48
.equiv          nb110nf_p_facel,        52
.equiv          nb110nf_p_krf,          56
.equiv          nb110nf_p_crf,          60
.equiv          nb110nf_Vc,             64
.equiv          nb110nf_type,           68
.equiv          nb110nf_p_ntype,        72
.equiv          nb110nf_vdwparam,       76
.equiv          nb110nf_Vvdw,           80
.equiv          nb110nf_p_tabscale,     84
.equiv          nb110nf_VFtab,          88
.equiv          nb110nf_invsqrta,       92
.equiv          nb110nf_dvda,           96
.equiv          nb110nf_p_gbtabscale,   100
.equiv          nb110nf_GBtab,          104
.equiv          nb110nf_p_nthreads,     108
.equiv          nb110nf_count,          112
.equiv          nb110nf_mtx,            116
.equiv          nb110nf_outeriter,      120
.equiv          nb110nf_inneriter,      124
.equiv          nb110nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb110nf_ix,             0
.equiv          nb110nf_iy,             16
.equiv          nb110nf_iz,             32
.equiv          nb110nf_iq,             48
.equiv          nb110nf_c6,             64
.equiv          nb110nf_c12,            80
.equiv          nb110nf_vctot,          96
.equiv          nb110nf_Vvdwtot,        112
.equiv          nb110nf_half,           128
.equiv          nb110nf_three,          144
.equiv          nb110nf_is3,            160
.equiv          nb110nf_ii3,            164
.equiv          nb110nf_ntia,           168
.equiv          nb110nf_innerjjnr,      172
.equiv          nb110nf_innerk,         176
.equiv          nb110nf_n,              180
.equiv          nb110nf_nn1,            184
.equiv          nb110nf_nri,            188
.equiv          nb110nf_facel,          192
.equiv          nb110nf_ntype,          196
.equiv          nb110nf_nouter,         200
.equiv          nb110nf_ninner,         204
.equiv          nb110nf_salign,         208
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  216		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb110nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb110nf_p_nri]
	mov esi, [ebp + nb110nf_p_facel]
	mov edi, [ebp + nb110nf_p_ntype]
	mov ecx, [ecx]
	mov esi, [esi]
	mov edi, [edi]
	mov [esp + nb110nf_nri], ecx
	mov [esp + nb110nf_facel], esi
	mov [esp + nb110nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb110nf_nouter], eax
	mov [esp + nb110nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb110nf_half], eax
	movss xmm1, [esp + nb110nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb110nf_half],  xmm1
	movaps [esp + nb110nf_three],  xmm3

.nb110nf_threadloop:
        mov   esi, [ebp + nb110nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb110nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb110nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb110nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb110nf_n], eax
        mov [esp + nb110nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb110nf_outerstart
        jmp .nb110nf_end

.nb110nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb110nf_nouter]
	mov [esp + nb110nf_nouter], ebx

.nb110nf_outer:
	mov   eax, [ebp + nb110nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb110nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb110nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb110nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx +esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb110nf_charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [esp + nb110nf_facel]
	shufps xmm3, xmm3, 0

    	mov   edx, [ebp + nb110nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb110nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb110nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb110nf_pos]    ;# eax = base of pos[]  

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + nb110nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [esp + nb110nf_ix], xmm0
	movaps [esp + nb110nf_iy], xmm1
	movaps [esp + nb110nf_iz], xmm2

	mov   [esp + nb110nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb110nf_vctot], xmm4
	movaps [esp + nb110nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb110nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb110nf_pos]
	mov   eax, [ebp + nb110nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb110nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb110nf_ninner]
	mov   [esp + nb110nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb110nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb110nf_unroll_loop
	jmp   .nb110nf_finish_inner
.nb110nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb110nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	add dword ptr [esp + nb110nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb110nf_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + nb110nf_iq]
	shufps xmm3, xmm6, 0
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov esi, [ebp + nb110nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]
	mov esi, [ebp + nb110nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb110nf_ntia]
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

	movaps [esp + nb110nf_c6], xmm4
	movaps [esp + nb110nf_c12], xmm6
	
	mov esi, [ebp + nb110nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [esp + nb110nf_ix]
	movaps xmm5, [esp + nb110nf_iy]
	movaps xmm6, [esp + nb110nf_iz]

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
	movaps xmm1, [esp + nb110nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb110nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [esp + nb110nf_c6]
	mulps  xmm2, [esp + nb110nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb110nf_Vvdwtot]
	addps  xmm3, [esp + nb110nf_vctot]
	movaps [esp + nb110nf_vctot], xmm3
	movaps [esp + nb110nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb110nf_innerk],  4
	jl    .nb110nf_finish_inner
	jmp   .nb110nf_unroll_loop
.nb110nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [esp + nb110nf_innerk],  4
	mov   edx, [esp + nb110nf_innerk]
	and   edx, 2
	jnz   .nb110nf_dopair
	jmp   .nb110nf_checksingle
.nb110nf_dopair:	
	mov esi, [ebp + nb110nf_charge]

    	mov   ecx, [esp + nb110nf_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add dword ptr [esp + nb110nf_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 12 ;# constant 00001100 
	shufps xmm3, xmm3, 88 ;# constant 01011000 ;# xmm3(0,1) has the charges 

	mov esi, [ebp + nb110nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [esi + ecx*4]
	mov edx, [esi + edx*4]	
	mov esi, [ebp + nb110nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [esp + nb110nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [esi + ecx*4]
	movhps xmm6, [esi + edx*4]
	mov edi, [ebp + nb110nf_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# constant 00001000 	
	shufps xmm6, xmm6, 13 ;# constant 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [esp + nb110nf_c6], xmm4
	movaps [esp + nb110nf_c12], xmm6	
			
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + nb110nf_iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb110nf_ix]
	movaps xmm5, [esp + nb110nf_iy]
	movaps xmm6, [esp + nb110nf_iz]

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
	movaps xmm1, [esp + nb110nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb110nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [esp + nb110nf_c6]
	mulps  xmm2, [esp + nb110nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addps  xmm5, [esp + nb110nf_Vvdwtot]
	addps  xmm3, [esp + nb110nf_vctot]
	movaps [esp + nb110nf_vctot], xmm3
	movaps [esp + nb110nf_Vvdwtot], xmm5

.nb110nf_checksingle:				
	mov   edx, [esp + nb110nf_innerk]
	and   edx, 1
	jnz    .nb110nf_dosingle
	jmp    .nb110nf_updateouterdata
.nb110nf_dosingle:			
	mov esi, [ebp + nb110nf_charge]
	mov edi, [ebp + nb110nf_pos]
	mov   ecx, [esp + nb110nf_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	;# xmm3(0) has the charge 	

	mov esi, [ebp + nb110nf_type]
	mov ecx, eax
	mov ecx, [esi + ecx*4]	
	mov esi, [ebp + nb110nf_vdwparam]
	shl ecx, 1
	add ecx, [esp + nb110nf_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [esi + ecx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# constant 11111100	
	shufps xmm6, xmm6, 253  ;# constant 11111101	
			
	movaps [esp + nb110nf_c6], xmm4
	movaps [esp + nb110nf_c12], xmm6	
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + nb110nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb110nf_ix]
	movaps xmm5, [esp + nb110nf_iy]
	movaps xmm6, [esp + nb110nf_iz]

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
	movaps xmm1, [esp + nb110nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb110nf_half]
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm1, xmm4
	mulps  xmm1, xmm4
	mulps  xmm1, xmm4	;# xmm1=rinvsix 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm1, [esp + nb110nf_c6]
	mulps  xmm2, [esp + nb110nf_c12]
	movaps xmm5, xmm2
	subps  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addss  xmm5, [esp + nb110nf_Vvdwtot]
	addss  xmm3, [esp + nb110nf_vctot]
	movss [esp + nb110nf_vctot], xmm3
	movss [esp + nb110nf_Vvdwtot], xmm5

.nb110nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb110nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb110nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb110nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb110nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [esp + nb110nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb110nf_Vvdw]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb110nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb110nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb110nf_n], esi
        jmp .nb110nf_outer
.nb110nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb110nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb110nf_end
        ;# non-zero, do one more workunit
        jmp   .nb110nf_threadloop
.nb110nf_end:
	emms

	mov eax, [esp + nb110nf_nouter]
	mov ebx, [esp + nb110nf_ninner]
	mov ecx, [ebp + nb110nf_outeriter]
	mov edx, [ebp + nb110nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb110nf_salign]
	add esp, eax
	add esp,  216
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

