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


.globl nb_kernel010_ia32_sse2
.globl _nb_kernel010_ia32_sse2
nb_kernel010_ia32_sse2:	
_nb_kernel010_ia32_sse2:	
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
.equiv          nb010_argkrf,           56
.equiv          nb010_argcrf,           60
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
	and  eax, 0xf
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
	mov eax, 0x00000000     ;# lower half of double 2.0 IEEE (hex)
	mov ebx, 0x40000000
	mov [esp + nb010_two], eax
	mov [esp + nb010_two+4], ebx
	movsd xmm1, [esp + nb010_two]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm2, xmm1
	addpd  xmm2, xmm1       ;# 4.0
	addpd  xmm2, xmm1       ;# 6.0
	movapd xmm3, xmm2
	addpd  xmm3, xmm3       ;# 12.0
	movapd [esp + nb010_two], xmm1
	movapd [esp + nb010_six],  xmm2
	movapd [esp + nb010_twelve], xmm3

.nb010_threadloop:
        mov   esi, [ebp + nb010_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb010_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
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
	mov   eax, [ebp + nb010_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 

	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb010_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb010_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb010_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

    	mov   edx, [ebp + nb010_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb010_ntype]
    	shl   edx, 1
    	mov   [esp + nb010_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb010_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb010_ix], xmm0
	movapd [esp + nb010_iy], xmm1
	movapd [esp + nb010_iz], xmm2

	mov   [esp + nb010_ii3], ebx
	
	;# clear Vvdwtot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb010_Vvdwtot], xmm4
	movapd [esp + nb010_fix], xmm4
	movapd [esp + nb010_fiy], xmm4
	movapd [esp + nb010_fiz], xmm4
	
	mov   eax, [ebp + nb010_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb010_pos]
	mov   edi, [ebp + nb010_faction]	
	mov   eax, [ebp + nb010_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb010_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb010_ninner]
	mov   [esp + nb010_ninner], ecx
	add   edx, 0
	mov   [esp + nb010_innerk], edx    ;# number of innerloop atoms 
	
	jge   .nb010_unroll_loop
	jmp   .nb010_checksingle
.nb010_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb010_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add   dword ptr [esp + nb010_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb010_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb010_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb010_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movlpd xmm7, [esi + ebx*8]	;# c6b
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [esi + ebx*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movd  eax, mm0		
	movd  ebx, mm1	
	
	movapd [esp + nb010_c6], xmm4
	movapd [esp + nb010_c12], xmm6

	mov esi, [ebp + nb010_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	

	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]

	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb010_ix]
	movapd xmm5, [esp + nb010_iy]
	movapd xmm6, [esp + nb010_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb010_dx], xmm4 
	movapd [esp + nb010_dy], xmm5
	movapd [esp + nb010_dz], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6	;# rsq in xmm4 

	cvtpd2ps xmm6, xmm4	
	rcpps xmm6, xmm6
	cvtps2pd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [esp + nb010_two]
	movapd xmm5, xmm4
	mulpd xmm4, xmm6	;# lu*rsq 
	subpd xmm0, xmm4	;# 2-lu*rsq 
	mulpd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [esp + nb010_two]
	mulpd xmm5, xmm6	;# lu*rsq 
	subpd xmm0, xmm5	;# 2-lu*rsq 
	mulpd xmm0, xmm6	;# xmm0=rinvsq 

	movapd xmm1, xmm0
	mulpd  xmm1, xmm0
	mulpd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulpd  xmm1, [esp + nb010_c6]
	mulpd  xmm2, [esp + nb010_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [esp + nb010_Vvdwtot]
	mulpd  xmm1, [esp + nb010_six]
	mulpd  xmm2, [esp + nb010_twelve]
	subpd  xmm2, xmm1
	mulpd  xmm0, xmm2	;# xmm4=total fscal 
	movapd xmm4, xmm0
	
	movapd xmm0, [esp + nb010_dx]
	movapd xmm1, [esp + nb010_dy]
	movapd xmm2, [esp + nb010_dz]

	movapd [esp + nb010_Vvdwtot], xmm5

	mov    edi, [ebp + nb010_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb010_fix]
	movapd xmm4, [esp + nb010_fiy]
	movapd xmm5, [esp + nb010_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb010_fix], xmm3
	movapd [esp + nb010_fiy], xmm4
	movapd [esp + nb010_fiz], xmm5
	
	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	movhpd xmm3, [edi + ebx*8]
	movhpd xmm4, [edi + ebx*8 + 8]
	movhpd xmm5, [edi + ebx*8 + 16]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5
	movhpd [edi + ebx*8], xmm3
	movhpd [edi + ebx*8 + 8], xmm4
	movhpd [edi + ebx*8 + 16], xmm5
		
	;# should we do one more iteration? 
	sub   dword ptr [esp + nb010_innerk],  2
	jl    .nb010_checksingle
	jmp   .nb010_unroll_loop
.nb010_checksingle:				
	mov   edx, [esp + nb010_innerk]
	and   edx, 1
	jnz    .nb010_dosingle
	jmp    .nb010_updateouterdata
.nb010_dosingle:
	mov edi, [ebp + nb010_pos]
	mov   ecx, [esp + nb010_innerjjnr]
	mov   eax, [ecx]		

	movd  mm0, eax		;# use mmx registers as temp storage 	
	mov esi, [ebp + nb010_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb010_vdwparam]
	shl eax, 1
	mov edi, [esp + nb010_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movd  eax, mm0		
	
	movapd [esp + nb010_c6], xmm4
	movapd [esp + nb010_c12], xmm6
	
	mov esi, [ebp + nb010_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	

	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb010_ix]
	movapd xmm5, [esp + nb010_iy]
	movapd xmm6, [esp + nb010_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb010_dx], xmm4 
	movapd [esp + nb010_dy], xmm5
	movapd [esp + nb010_dz], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6	;# rsq in xmm4 

	cvtsd2ss xmm6, xmm4	
	rcpss xmm6, xmm6
	cvtss2sd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [esp + nb010_two]
	movapd xmm5, xmm4
	mulsd xmm4, xmm6	;# lu*rsq 
	subsd xmm0, xmm4	;# 2-lu*rsq 
	mulsd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [esp + nb010_two]
	mulsd xmm5, xmm6	;# lu*rsq 
	subsd xmm0, xmm5	;# 2-lu*rsq 
	mulsd xmm0, xmm6	;# xmm0=rinvsq 
	movapd xmm4, xmm0
	
	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulsd  xmm1, [esp + nb010_c6]
	mulsd  xmm2, [esp + nb010_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [esp + nb010_Vvdwtot]
	mulsd  xmm1, [esp + nb010_six]
	mulsd  xmm2, [esp + nb010_twelve]
	subsd  xmm2, xmm1
	mulsd  xmm4, xmm2	;# xmm4=total fscal 

	movapd xmm0, [esp + nb010_dx]
	movapd xmm1, [esp + nb010_dy]
	movapd xmm2, [esp + nb010_dz]

	movlpd [esp + nb010_Vvdwtot], xmm5

	mov    edi, [ebp + nb010_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movlpd xmm3, [esp + nb010_fix]
	movlpd xmm4, [esp + nb010_fiy]
	movlpd xmm5, [esp + nb010_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb010_fix], xmm3
	movlpd [esp + nb010_fiy], xmm4
	movlpd [esp + nb010_fiz], xmm5
	
	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5

.nb010_updateouterdata:
	mov   ecx, [esp + nb010_ii3]
	mov   edi, [ebp + nb010_faction]
	mov   esi, [ebp + nb010_fshift]
	mov   edx, [esp + nb010_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb010_fix]
	movapd xmm1, [esp + nb010_fiy]
	movapd xmm2, [esp + nb010_fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [edi + ecx*8]
	movsd  xmm4, [edi + ecx*8 + 8]
	movsd  xmm5, [edi + ecx*8 + 16]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [edi + ecx*8],     xmm3
	movsd  [edi + ecx*8 + 8], xmm4
	movsd  [edi + ecx*8 + 16], xmm5

	;# increment fshift force  
	movsd  xmm3, [esi + edx*8]
	movsd  xmm4, [esi + edx*8 + 8]
	movsd  xmm5, [esi + edx*8 + 16]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [esi + edx*8],     xmm3
	movsd  [esi + edx*8 + 8], xmm4
	movsd  [esi + edx*8 + 16], xmm5

	;# get n from stack
	mov esi, [esp + nb010_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb010_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb010_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 have the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb010_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 

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
	add esp, eax
	add esp, 320
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret


	


	
.globl nb_kernel010nf_ia32_sse2
.globl _nb_kernel010nf_ia32_sse2
nb_kernel010nf_ia32_sse2:	
_nb_kernel010nf_ia32_sse2:	
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
.equiv          nb010nf_argkrf,         56
.equiv          nb010nf_argcrf,         60
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
        mov eax, 0x00000000     ;# lower half of double 2.0 IEEE (hex)
        mov ebx, 0x40000000
        mov [esp + nb010nf_two], eax
        mov [esp + nb010nf_two+4], ebx
        movsd xmm1, [esp + nb010nf_two]
        shufpd xmm1, xmm1, 0    ;# splat to all elements
        movapd [esp + nb010nf_two], xmm1

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
	mov   eax, [ebp + nb010nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 

	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb010nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb010nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

    	mov   edx, [ebp + nb010nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb010nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb010nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb010nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb010nf_ix], xmm0
	movapd [esp + nb010nf_iy], xmm1
	movapd [esp + nb010nf_iz], xmm2

	mov   [esp + nb010nf_ii3], ebx
	
	;# clear Vvdwtot 
	xorpd xmm4, xmm4
	movapd [esp + nb010nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb010nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb010nf_pos]
	mov   eax, [ebp + nb010nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb010nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb010nf_ninner]
	mov   [esp + nb010nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb010nf_innerk], edx    ;# number of innerloop atoms 
	
	jge   .nb010nf_unroll_loop
	jmp   .nb010nf_checksingle
.nb010nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb010nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add   dword ptr [esp + nb010nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb010nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb010nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb010nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movlpd xmm7, [esi + ebx*8]	;# c6b
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [esi + ebx*8 + 8]	;# c6b c12b 

	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movd  eax, mm0		
	movd  ebx, mm1	
	
	movapd [esp + nb010nf_c6], xmm4
	movapd [esp + nb010nf_c12], xmm6

	mov esi, [ebp + nb010nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	

	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]

	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb010nf_ix]
	movapd xmm5, [esp + nb010nf_iy]
	movapd xmm6, [esp + nb010nf_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6	;# rsq in xmm4 

	cvtpd2ps xmm6, xmm4	
	rcpps xmm6, xmm6
	cvtps2pd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [esp + nb010nf_two]
	movapd xmm5, xmm4
	mulpd xmm4, xmm6	;# lu*rsq 
	subpd xmm0, xmm4	;# 2-lu*rsq 
	mulpd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [esp + nb010nf_two]
	mulpd xmm5, xmm6	;# lu*rsq 
	subpd xmm0, xmm5	;# 2-lu*rsq 
	mulpd xmm0, xmm6	;# xmm0=rinvsq 

	movapd xmm1, xmm0
	mulpd  xmm1, xmm0
	mulpd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulpd  xmm1, [esp + nb010nf_c6]
	mulpd  xmm2, [esp + nb010nf_c12]
	movapd xmm5, xmm2
	subpd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addpd  xmm5, [esp + nb010nf_Vvdwtot]
	movapd [esp + nb010nf_Vvdwtot], xmm5
		
	;# should we do one more iteration? 
	sub   dword ptr [esp + nb010nf_innerk],  2
	jl    .nb010nf_checksingle
	jmp   .nb010nf_unroll_loop
.nb010nf_checksingle:				
	mov   edx, [esp + nb010nf_innerk]
	and   edx, 1
	jnz    .nb010nf_dosingle
	jmp    .nb010nf_updateouterdata
.nb010nf_dosingle:
	mov edi, [ebp + nb010nf_pos]
	mov   ecx, [esp + nb010nf_innerjjnr]
	mov   eax, [ecx]		

	movd  mm0, eax		;# use mmx registers as temp storage 	
	mov esi, [ebp + nb010nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb010nf_vdwparam]
	shl eax, 1
	mov edi, [esp + nb010nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7

	movd  eax, mm0		
	
	movapd [esp + nb010nf_c6], xmm4
	movapd [esp + nb010nf_c12], xmm6
	
	mov esi, [ebp + nb010nf_pos]       ;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	

	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb010nf_ix]
	movapd xmm5, [esp + nb010nf_iy]
	movapd xmm6, [esp + nb010nf_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6	;# rsq in xmm4 

	cvtsd2ss xmm6, xmm4	
	rcpss xmm6, xmm6
	cvtss2sd xmm6, xmm6	;# lu in low xmm6 
	
	;# 1/x lookup seed in xmm6 
	movapd xmm0, [esp + nb010nf_two]
	movapd xmm5, xmm4
	mulsd xmm4, xmm6	;# lu*rsq 
	subsd xmm0, xmm4	;# 2-lu*rsq 
	mulsd xmm6, xmm0	;# (new lu) 
	
	movapd xmm0, [esp + nb010nf_two]
	mulsd xmm5, xmm6	;# lu*rsq 
	subsd xmm0, xmm5	;# 2-lu*rsq 
	mulsd xmm0, xmm6	;# xmm0=rinvsq 
	movapd xmm4, xmm0
	
	movapd xmm1, xmm0
	mulsd  xmm1, xmm0
	mulsd  xmm1, xmm0	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 

	mulsd  xmm1, [esp + nb010nf_c6]
	mulsd  xmm2, [esp + nb010nf_c12]
	movapd xmm5, xmm2
	subsd  xmm5, xmm1	;# Vvdw=Vvdw12-Vvdw6 
	addsd  xmm5, [esp + nb010nf_Vvdwtot]
	movlpd [esp + nb010nf_Vvdwtot], xmm5

.nb010nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb010nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb010nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb010nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 have the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb010nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 

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

