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



.globl nb_kernel200_ia32_sse
.globl _nb_kernel200_ia32_sse
nb_kernel200_ia32_sse:	
_nb_kernel200_ia32_sse:	
.equiv          nb200_p_nri,            8
.equiv          nb200_iinr,             12
.equiv          nb200_jindex,           16
.equiv          nb200_jjnr,             20
.equiv          nb200_shift,            24
.equiv          nb200_shiftvec,         28
.equiv          nb200_fshift,           32
.equiv          nb200_gid,              36
.equiv          nb200_pos,              40
.equiv          nb200_faction,          44
.equiv          nb200_charge,           48
.equiv          nb200_p_facel,          52
.equiv          nb200_argkrf,           56
.equiv          nb200_argcrf,           60
.equiv          nb200_Vc,               64
.equiv          nb200_type,             68
.equiv          nb200_p_ntype,          72
.equiv          nb200_vdwparam,         76
.equiv          nb200_Vvdw,             80
.equiv          nb200_p_tabscale,       84
.equiv          nb200_VFtab,            88
.equiv          nb200_invsqrta,         92
.equiv          nb200_dvda,             96
.equiv          nb200_p_gbtabscale,     100
.equiv          nb200_GBtab,            104
.equiv          nb200_p_nthreads,       108
.equiv          nb200_count,            112
.equiv          nb200_mtx,              116
.equiv          nb200_outeriter,        120
.equiv          nb200_inneriter,        124
.equiv          nb200_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb200_ix,               0
.equiv          nb200_iy,               16
.equiv          nb200_iz,               32
.equiv          nb200_iq,               48
.equiv          nb200_dx,               64
.equiv          nb200_dy,               80
.equiv          nb200_dz,               96
.equiv          nb200_vctot,            112
.equiv          nb200_fix,              128
.equiv          nb200_fiy,              144
.equiv          nb200_fiz,              160
.equiv          nb200_half,             176
.equiv          nb200_three,            192
.equiv          nb200_two,              208
.equiv          nb200_krf,              224
.equiv          nb200_crf,              240
.equiv          nb200_is3,              256
.equiv          nb200_ii3,              260
.equiv          nb200_innerjjnr,        264
.equiv          nb200_innerk,           268
.equiv          nb200_n,                272
.equiv          nb200_nn1,              276
.equiv          nb200_nri,              280
.equiv          nb200_facel,            284
.equiv          nb200_nouter,           288
.equiv          nb200_ninner,           292
.equiv          nb200_salign,           296
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  300		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb200_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb200_p_nri]
	mov esi, [ebp + nb200_p_facel]
	mov ecx, [ecx]
	mov esi, [esi]
	mov [esp + nb200_nri], ecx
	mov [esp + nb200_facel], esi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb200_nouter], eax
	mov [esp + nb200_ninner], eax


	mov esi, [ebp + nb200_argkrf]
	mov edi, [ebp + nb200_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	movaps [esp + nb200_krf], xmm5
	shufps xmm6, xmm6, 0
	movaps [esp + nb200_crf], xmm6

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb200_half], eax
	movss xmm1, [esp + nb200_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb200_half],  xmm1
	movaps [esp + nb200_two],  xmm2
	movaps [esp + nb200_three],  xmm3

.nb200_threadloop:
        mov   esi, [ebp + nb200_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb200_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb200_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb200_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb200_n], eax
        mov [esp + nb200_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb200_outerstart
        jmp .nb200_end

	;# assume we have at least one i particle - start directly 	
.nb200_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb200_nouter]
	mov [esp + nb200_nouter], ebx

.nb200_outer:
	mov   eax, [ebp + nb200_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb200_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb200_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb200_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb200_charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [esp + nb200_facel]
	shufps xmm3, xmm3, 0
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb200_pos]    ;# eax = base of pos[]  

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + nb200_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [esp + nb200_ix], xmm0
	movaps [esp + nb200_iy], xmm1
	movaps [esp + nb200_iz], xmm2

	mov   [esp + nb200_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb200_vctot], xmm4
	movaps [esp + nb200_fix], xmm4
	movaps [esp + nb200_fiy], xmm4
	movaps [esp + nb200_fiz], xmm4
	
	mov   eax, [ebp + nb200_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb200_pos]
	mov   edi, [ebp + nb200_faction]	
	mov   eax, [ebp + nb200_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb200_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb200_ninner]
	mov   [esp + nb200_ninner], ecx
	add   edx, 0
	mov   [esp + nb200_innerk], edx    ;# number of innerloop atoms 
	jge   .nb200_unroll_loop
	jmp   .nb200_finish_inner
.nb200_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb200_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	add dword ptr [esp + nb200_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb200_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + nb200_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  

	mov esi, [ebp + nb200_pos]       ;# base of pos[] 

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
	movaps xmm4, [esp + nb200_ix]
	movaps xmm5, [esp + nb200_iy]
	movaps xmm6, [esp + nb200_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb200_dx], xmm4
	movaps [esp + nb200_dy], xmm5
	movaps [esp + nb200_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 
	
	movaps xmm7, [esp + nb200_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb200_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb200_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 

	subps  xmm6, [esp + nb200_crf] ;# xmm6=rinv+ krsq-crf 

	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	mulps  xmm7, [esp + nb200_two]

	subps  xmm0, xmm7
	mulps  xmm3, xmm0	
	mulps  xmm4, xmm3	;# xmm4=total fscal 
	addps  xmm6, [esp + nb200_vctot]

	movaps xmm0, [esp + nb200_dx]
	movaps xmm1, [esp + nb200_dy]
	movaps xmm2, [esp + nb200_dz]

	movaps [esp + nb200_vctot], xmm6

	mov    edi, [ebp + nb200_faction]
	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb200_fix]
	movaps xmm4, [esp + nb200_fiy]
	movaps xmm5, [esp + nb200_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb200_fix], xmm3
	movaps [esp + nb200_fiy], xmm4
	movaps [esp + nb200_fiz], xmm5
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
	sub dword ptr [esp + nb200_innerk],  4
	jl    .nb200_finish_inner
	jmp   .nb200_unroll_loop
.nb200_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [esp + nb200_innerk],  4
	mov   edx, [esp + nb200_innerk]
	and   edx, 2
	jnz   .nb200_dopair
	jmp   .nb200_checksingle
.nb200_dopair:	
	mov esi, [ebp + nb200_charge]

    	mov   ecx, [esp + nb200_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add dword ptr [esp + nb200_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 12 ;# constant 00001100 
	shufps xmm3, xmm3, 88 ;# constant 01011000 ;# xmm3(0,1) has the charges 	

	mov edi, [ebp + nb200_pos]	
				
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + nb200_iq]

	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	mov    edi, [ebp + nb200_faction]
	;# move ix-iz to xmm4-xmm6 
	
	movaps xmm4, [esp + nb200_ix]
	movaps xmm5, [esp + nb200_iy]
	movaps xmm6, [esp + nb200_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb200_dx], xmm4
	movaps [esp + nb200_dy], xmm5
	movaps [esp + nb200_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	movaps xmm7, [esp + nb200_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb200_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb200_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 

	subps  xmm6, [esp + nb200_crf] ;# xmm6=rinv+ krsq-crf 

	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	mulps  xmm7, [esp + nb200_two]	

	subps  xmm0, xmm7
	mulps  xmm3, xmm0	

	mulps  xmm4, xmm3	;# xmm4=total fscal 
	addps  xmm6, [esp + nb200_vctot]

	movaps xmm0, [esp + nb200_dx]
	movaps xmm1, [esp + nb200_dy]
	movaps xmm2, [esp + nb200_dz]

	movaps [esp + nb200_vctot], xmm6

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb200_fix]
	movaps xmm4, [esp + nb200_fiy]
	movaps xmm5, [esp + nb200_fiz]
	addps  xmm3, xmm0
	addps  xmm4, xmm1
	addps  xmm5, xmm2
	movaps [esp + nb200_fix], xmm3
	movaps [esp + nb200_fiy], xmm4
	movaps [esp + nb200_fiz], xmm5
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

.nb200_checksingle:				
	mov   edx, [esp + nb200_innerk]
	and   edx, 1
	jnz    .nb200_dosingle
	jmp    .nb200_updateouterdata
.nb200_dosingle:			
	mov esi, [ebp + nb200_charge]
	mov edi, [ebp + nb200_pos]
	mov   ecx, [esp + nb200_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	;# xmm3(0) has the charge 		
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + nb200_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb200_ix]
	movaps xmm5, [esp + nb200_iy]
	movaps xmm6, [esp + nb200_iz]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# store dr 
	movaps [esp + nb200_dx], xmm4
	movaps [esp + nb200_dy], xmm5
	movaps [esp + nb200_dz], xmm6
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

	movaps xmm7, [esp + nb200_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb200_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb200_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 

	subps  xmm6, [esp + nb200_crf] ;# xmm6=rinv+ krsq-crf 

	mulps  xmm6, xmm3	;# xmm6=vcoul 
	mulps  xmm7, [esp + nb200_two]

	subps  xmm0, xmm7
	mulps  xmm3, xmm0
	mulps  xmm4, xmm3	;# xmm4=total fscal 
	addss  xmm6, [esp + nb200_vctot]
	
	mov    edi, [ebp + nb200_faction]

	movaps xmm0, [esp + nb200_dx]
	movaps xmm1, [esp + nb200_dy]
	movaps xmm2, [esp + nb200_dz]

	movss [esp + nb200_vctot], xmm6

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movaps xmm3, [esp + nb200_fix]
	movaps xmm4, [esp + nb200_fiy]
	movaps xmm5, [esp + nb200_fiz]
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movaps [esp + nb200_fix], xmm3
	movaps [esp + nb200_fiy], xmm4
	movaps [esp + nb200_fiz], xmm5
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
.nb200_updateouterdata:
	mov   ecx, [esp + nb200_ii3]
	mov   edi, [ebp + nb200_faction]
	mov   esi, [ebp + nb200_fshift]
	mov   edx, [esp + nb200_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [esp + nb200_fix]
	movaps xmm1, [esp + nb200_fiy]
	movaps xmm2, [esp + nb200_fiz]

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
	mov esi, [esp + nb200_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb200_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb200_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb200_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb200_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb200_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb200_n], esi
        jmp .nb200_outer
.nb200_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb200_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb200_end
        ;# non-zero, do one more workunit
        jmp   .nb200_threadloop
.nb200_end:
	emms

	mov eax, [esp + nb200_nouter]
	mov ebx, [esp + nb200_ninner]
	mov ecx, [ebp + nb200_outeriter]
	mov edx, [ebp + nb200_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb200_salign]
	add esp, eax
	add esp,  300
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret





.globl nb_kernel200nf_ia32_sse
.globl _nb_kernel200nf_ia32_sse
nb_kernel200nf_ia32_sse:	
_nb_kernel200nf_ia32_sse:	
.equiv          nb200nf_p_nri,          8
.equiv          nb200nf_iinr,           12
.equiv          nb200nf_jindex,         16
.equiv          nb200nf_jjnr,           20
.equiv          nb200nf_shift,          24
.equiv          nb200nf_shiftvec,       28
.equiv          nb200nf_fshift,         32
.equiv          nb200nf_gid,            36
.equiv          nb200nf_pos,            40
.equiv          nb200nf_faction,        44
.equiv          nb200nf_charge,         48
.equiv          nb200nf_p_facel,        52
.equiv          nb200nf_argkrf,         56
.equiv          nb200nf_argcrf,         60
.equiv          nb200nf_Vc,             64
.equiv          nb200nf_type,           68
.equiv          nb200nf_p_ntype,        72
.equiv          nb200nf_vdwparam,       76
.equiv          nb200nf_Vvdw,           80
.equiv          nb200nf_p_tabscale,     84
.equiv          nb200nf_VFtab,          88
.equiv          nb200nf_invsqrta,       92
.equiv          nb200nf_dvda,           96
.equiv          nb200nf_p_gbtabscale,   100
.equiv          nb200nf_GBtab,          104
.equiv          nb200nf_p_nthreads,     108
.equiv          nb200nf_count,          112
.equiv          nb200nf_mtx,            116
.equiv          nb200nf_outeriter,      120
.equiv          nb200nf_inneriter,      124
.equiv          nb200nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb200nf_ix,             0
.equiv          nb200nf_iy,             16
.equiv          nb200nf_iz,             32
.equiv          nb200nf_iq,             48
.equiv          nb200nf_vctot,          64
.equiv          nb200nf_half,           80
.equiv          nb200nf_three,          96
.equiv          nb200nf_krf,            112
.equiv          nb200nf_crf,            128
.equiv          nb200nf_is3,            144
.equiv          nb200nf_ii3,            148
.equiv          nb200nf_innerjjnr,      152
.equiv          nb200nf_innerk,         156
.equiv          nb200nf_n,              160
.equiv          nb200nf_nn1,            164
.equiv          nb200nf_nri,            168
.equiv          nb200nf_facel,          172
.equiv          nb200nf_nouter,         176
.equiv          nb200nf_ninner,         180
.equiv          nb200nf_salign,         184
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp,  188		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb200nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb200nf_p_nri]
	mov esi, [ebp + nb200nf_p_facel]
	mov ecx, [ecx]
	mov esi, [esi]
	mov [esp + nb200nf_nri], ecx
	mov [esp + nb200nf_facel], esi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb200nf_nouter], eax
	mov [esp + nb200nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# constant 0.5 in IEEE (hex)
	mov [esp + nb200nf_half], eax
	movss xmm1, [esp + nb200nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# constant 1.0
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# constant 2.0
	addps  xmm3, xmm2	;# constant 3.0
	movaps [esp + nb200nf_half],  xmm1
	movaps [esp + nb200nf_three],  xmm3

	mov esi, [ebp + nb200nf_argkrf]
	mov edi, [ebp + nb200nf_argcrf]
	movss xmm5, [esi]
	movss xmm6, [edi]
	shufps xmm5, xmm5, 0
	movaps [esp + nb200nf_krf], xmm5
	shufps xmm6, xmm6, 0
	movaps [esp + nb200nf_crf], xmm6

.nb200nf_threadloop:
        mov   esi, [ebp + nb200nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb200nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb200nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb200nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb200nf_n], eax
        mov [esp + nb200nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb200nf_outerstart
        jmp .nb200nf_end

.nb200nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb200nf_nouter]
	mov [esp + nb200nf_nouter], ebx

.nb200nf_outer:
	mov   eax, [ebp + nb200nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb200nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb200nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [eax + ebx*4]
	movss xmm1, [eax + ebx*4 + 4]
	movss xmm2, [eax + ebx*4 + 8] 

	mov   ecx, [ebp + nb200nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb200nf_charge]
	movss xmm3, [edx + ebx*4]	
	mulss xmm3, [esp + nb200nf_facel]
	shufps xmm3, xmm3, 0
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb200nf_pos]    ;# eax = base of pos[]  

	addss xmm0, [eax + ebx*4]
	addss xmm1, [eax + ebx*4 + 4]
	addss xmm2, [eax + ebx*4 + 8]

	movaps [esp + nb200nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [esp + nb200nf_ix], xmm0
	movaps [esp + nb200nf_iy], xmm1
	movaps [esp + nb200nf_iz], xmm2

	mov   [esp + nb200nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [esp + nb200nf_vctot], xmm4
	
	mov   eax, [ebp + nb200nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb200nf_pos]
	mov   eax, [ebp + nb200nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb200nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [esp + nb200nf_ninner]
	mov   [esp + nb200nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb200nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb200nf_unroll_loop
	jmp   .nb200nf_finish_inner
.nb200nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   edx, [esp + nb200nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	mov   ecx, [edx + 8]            
	mov   edx, [edx + 12]         ;# eax-edx=jnr1-4 
	add dword ptr [esp + nb200nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov esi, [ebp + nb200nf_charge]    ;# base of charge[] 
	
	movss xmm3, [esi + eax*4]
	movss xmm4, [esi + ecx*4]
	movss xmm6, [esi + ebx*4]
	movss xmm7, [esi + edx*4]

	movaps xmm2, [esp + nb200nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  

	mov esi, [ebp + nb200nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [esp + nb200nf_ix]
	movaps xmm5, [esp + nb200nf_iy]
	movaps xmm6, [esp + nb200nf_iz]

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
	
	movaps xmm7, [esp + nb200nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb200nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb200nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb200nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq) 
	addps  xmm6, [esp + nb200nf_vctot]
	movaps [esp + nb200nf_vctot], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb200nf_innerk],  4
	jl    .nb200nf_finish_inner
	jmp   .nb200nf_unroll_loop
.nb200nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [esp + nb200nf_innerk],  4
	mov   edx, [esp + nb200nf_innerk]
	and   edx, 2
	jnz   .nb200nf_dopair
	jmp   .nb200nf_checksingle
.nb200nf_dopair:	
	mov esi, [ebp + nb200nf_charge]

    	mov   ecx, [esp + nb200nf_innerjjnr]
	
	mov   eax, [ecx]	
	mov   ebx, [ecx + 4]              
	add dword ptr [esp + nb200nf_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [esi + eax*4]		
	movss xmm6, [esi + ebx*4]
	shufps xmm3, xmm6, 12 ;# constant 00001100 
	shufps xmm3, xmm3, 88 ;# constant 01011000 ;# xmm3(0,1) has the charges 	

	mov edi, [ebp + nb200nf_pos]	
				
	lea   eax, [eax + eax*2]
	lea   ebx, [ebx + ebx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [edi + eax*4]
	movss xmm2, [edi + eax*4 + 8]	
	movhps xmm1, [edi + ebx*4]
	movss xmm0, [edi + ebx*4 + 8]	

	mulps  xmm3, [esp + nb200nf_iq]

	xorps  xmm7,xmm7
	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	;# move ix-iz to xmm4-xmm6 
	
	movaps xmm4, [esp + nb200nf_ix]
	movaps xmm5, [esp + nb200nf_iy]
	movaps xmm6, [esp + nb200nf_iz]

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

	movaps xmm7, [esp + nb200nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb200nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb200nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb200nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 

	addps  xmm6, [esp + nb200nf_vctot]
	movaps [esp + nb200nf_vctot], xmm6

.nb200nf_checksingle:				
	mov   edx, [esp + nb200nf_innerk]
	and   edx, 1
	jnz    .nb200nf_dosingle
	jmp    .nb200nf_updateouterdata
.nb200nf_dosingle:			
	mov esi, [ebp + nb200nf_charge]
	mov edi, [ebp + nb200nf_pos]
	mov   ecx, [esp + nb200nf_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [ecx]
	movss xmm3, [esi + eax*4]	;# xmm3(0) has the charge 
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [edi + eax*4]	
	movss xmm1, [edi + eax*4 + 4]	
	movss xmm2, [edi + eax*4 + 8]	
 
	mulps  xmm3, [esp + nb200nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [esp + nb200nf_ix]
	movaps xmm5, [esp + nb200nf_iy]
	movaps xmm6, [esp + nb200nf_iz]

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

	movaps xmm7, [esp + nb200nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [esp + nb200nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [esp + nb200nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [esp + nb200nf_crf] ;# xmm6=rinv+ krsq-crf 
	mulps  xmm6, xmm3	;# xmm6=vcoul 
	addss  xmm6, [esp + nb200nf_vctot]
	movss [esp + nb200nf_vctot], xmm6

.nb200nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb200nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb200nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [esp + nb200nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   eax, [ebp + nb200nf_Vc]
	addss xmm7, [eax + edx*4] 
	;# move back to mem 
	movss [eax + edx*4], xmm7 
	
       ;# finish if last 
        mov ecx, [esp + nb200nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb200nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb200nf_n], esi
        jmp .nb200nf_outer
.nb200nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb200nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb200nf_end
        ;# non-zero, do one more workunit
        jmp   .nb200nf_threadloop
.nb200nf_end:
	emms

	mov eax, [esp + nb200nf_nouter]
	mov ebx, [esp + nb200nf_ninner]
	mov ecx, [ebp + nb200nf_outeriter]
	mov edx, [ebp + nb200nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb200nf_salign]
	add esp, eax
	add esp,  188
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret
	
