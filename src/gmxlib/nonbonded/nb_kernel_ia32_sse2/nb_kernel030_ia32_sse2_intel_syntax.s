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



.globl nb_kernel030_ia32_sse2
.globl _nb_kernel030_ia32_sse2
nb_kernel030_ia32_sse2:	
_nb_kernel030_ia32_sse2:	
.equiv          nb030_p_nri,            8
.equiv          nb030_iinr,             12
.equiv          nb030_jindex,           16
.equiv          nb030_jjnr,             20
.equiv          nb030_shift,            24
.equiv          nb030_shiftvec,         28
.equiv          nb030_fshift,           32
.equiv          nb030_gid,              36
.equiv          nb030_pos,              40
.equiv          nb030_faction,          44
.equiv          nb030_charge,           48
.equiv          nb030_p_facel,          52
.equiv          nb030_argkrf,           56
.equiv          nb030_argcrf,           60
.equiv          nb030_Vc,               64
.equiv          nb030_type,             68
.equiv          nb030_p_ntype,          72
.equiv          nb030_vdwparam,         76
.equiv          nb030_Vvdw,             80
.equiv          nb030_p_tabscale,       84
.equiv          nb030_VFtab,            88
.equiv          nb030_invsqrta,         92
.equiv          nb030_dvda,             96
.equiv          nb030_p_gbtabscale,     100
.equiv          nb030_GBtab,            104
.equiv          nb030_p_nthreads,       108
.equiv          nb030_count,            112
.equiv          nb030_mtx,              116
.equiv          nb030_outeriter,        120
.equiv          nb030_inneriter,        124
.equiv          nb030_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb030_ix,               0
.equiv          nb030_iy,               16
.equiv          nb030_iz,               32
.equiv          nb030_dx,               48
.equiv          nb030_dy,               64
.equiv          nb030_dz,               80
.equiv          nb030_two,              96
.equiv          nb030_tsc,              112
.equiv          nb030_c6,               128
.equiv          nb030_c12,              144
.equiv          nb030_fscal,            160
.equiv          nb030_Vvdwtot,          176
.equiv          nb030_fix,              192
.equiv          nb030_fiy,              208
.equiv          nb030_fiz,              224
.equiv          nb030_half,             240
.equiv          nb030_three,            256
.equiv          nb030_is3,              272
.equiv          nb030_ii3,              276
.equiv          nb030_ntia,             280
.equiv          nb030_innerjjnr,        284
.equiv          nb030_innerk,           288
.equiv          nb030_n,                292
.equiv          nb030_nn1,              296
.equiv          nb030_nri,              300
.equiv          nb030_ntype,            304
.equiv          nb030_nouter,           308
.equiv          nb030_ninner,           312
.equiv          nb030_salign,           316
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
	mov [esp + nb030_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb030_p_nri]
	mov edi, [ebp + nb030_p_ntype]
	mov ecx, [ecx]
	mov edi, [edi]
	mov [esp + nb030_nri], ecx
	mov [esp + nb030_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb030_nouter], eax
	mov [esp + nb030_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb030_half], eax
	mov [esp + nb030_half+4], ebx
	movsd xmm1, [esp + nb030_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb030_half], xmm1
	movapd [esp + nb030_two],  xmm2
	movapd [esp + nb030_three], xmm3
	mov eax, [ebp + nb030_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb030_tsc], xmm3

.nb030_threadloop:
        mov   esi, [ebp + nb030_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb030_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb030_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb030_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb030_n], eax
        mov [esp + nb030_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb030_outerstart
        jmp .nb030_end

.nb030_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb030_nouter]
	mov [esp + nb030_nouter], ebx

.nb030_outer:
	mov   eax, [ebp + nb030_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb030_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb030_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb030_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx + esi*4]	    ;# ebx =ii 

    	mov   edx, [ebp + nb030_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb030_ntype]
    	shl   edx, 1
    	mov   [esp + nb030_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb030_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb030_ix], xmm0
	movapd [esp + nb030_iy], xmm1
	movapd [esp + nb030_iz], xmm2

	mov   [esp + nb030_ii3], ebx
	
	;# clear tot potential and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb030_Vvdwtot], xmm4
	movapd [esp + nb030_fix], xmm4
	movapd [esp + nb030_fiy], xmm4
	movapd [esp + nb030_fiz], xmm4
	
	mov   eax, [ebp + nb030_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb030_pos]
	mov   edi, [ebp + nb030_faction]	
	mov   eax, [ebp + nb030_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb030_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb030_ninner]
	mov   [esp + nb030_ninner], ecx
	add   edx, 0
	mov   [esp + nb030_innerk], edx    ;# number of innerloop atoms 
	jge   .nb030_unroll_loop
	jmp   .nb030_checksingle
.nb030_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb030_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb030_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb030_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb030_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb030_ntia]
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

	movapd [esp + nb030_c6], xmm4
	movapd [esp + nb030_c12], xmm6
	
	mov esi, [ebp + nb030_pos]		;# base of pos[] 
	lea   eax, [eax + eax*2]	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]

	;# move nb030_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb030_ix]
	movapd xmm5, [esp + nb030_iy]
	movapd xmm6, [esp + nb030_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb030_dx], xmm4
	movapd [esp + nb030_dy], xmm5
	movapd [esp + nb030_dz], xmm6
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb030_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb030_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb030_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb030_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb030_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb030_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb030_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb030_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb030_Vvdwtot]
	movapd [esp + nb030_fscal], xmm7
	movapd [esp + nb030_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movlpd xmm3, [esi + ebx*8 + 32]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movlpd xmm3, [esi + ebx*8 + 48]	;# G2
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb030_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb030_c12]
	mulpd  xmm7, xmm4 
	mulpd  xmm5, xmm4  
	addpd  xmm7, [esp + nb030_fscal] 
	
	addpd  xmm5, [esp + nb030_Vvdwtot]
	movapd [esp + nb030_Vvdwtot], xmm5
	xorpd  xmm4, xmm4

	mulpd xmm7, [esp + nb030_tsc]
	mulpd xmm7, xmm0
	subpd  xmm4, xmm7

	movapd xmm0, [esp + nb030_dx]
	movapd xmm1, [esp + nb030_dy]
	movapd xmm2, [esp + nb030_dz]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb030_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb030_fix]
	movapd xmm4, [esp + nb030_fiy]
	movapd xmm5, [esp + nb030_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb030_fix], xmm3
	movapd [esp + nb030_fiy], xmm4
	movapd [esp + nb030_fiz], xmm5
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
	sub dword ptr [esp + nb030_innerk],  2
	jl    .nb030_checksingle
	jmp   .nb030_unroll_loop

.nb030_checksingle:				
	mov   edx, [esp + nb030_innerk]
	and   edx, 1
	jnz    .nb030_dosingle
	jmp    .nb030_updateouterdata
.nb030_dosingle:
	mov   edx, [esp + nb030_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	

	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov esi, [ebp + nb030_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb030_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb030_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		

	movapd [esp + nb030_c6], xmm4
	movapd [esp + nb030_c12], xmm6
	
	mov esi, [ebp + nb030_pos]		;# base of pos[] 
	lea   eax, [eax + eax*2]	;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move nb030_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb030_ix]
	movapd xmm5, [esp + nb030_iy]
	movapd xmm6, [esp + nb030_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb030_dx], xmm4
	movapd [esp + nb030_dy], xmm5
	movapd [esp + nb030_dz], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb030_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb030_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb030_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb030_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb030_tsc]

	movd mm0, eax
	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	shl eax, 3	

	mov  esi, [ebp + nb030_VFtab]

	;# dispersion 
	movlpd xmm4, [esi + eax*8]	;# Y1
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1 	
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	

	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1  
	unpckhpd xmm7, xmm3	;# H1  
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb030_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb030_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb030_Vvdwtot]
	movlpd [esp + nb030_fscal], xmm7
	movlpd [esp + nb030_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 

	xorpd xmm3,xmm3	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 48]	;# G1 	
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1  
	unpckhpd xmm7, xmm3	;# H1  
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb030_two]	;# two*Heps2 
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movapd xmm4, [esp + nb030_c12]
	mulsd  xmm7, xmm4 
	mulsd  xmm5, xmm4  
	addsd  xmm7, [esp + nb030_fscal] 
	
	addsd  xmm5, [esp + nb030_Vvdwtot]
	movlpd [esp + nb030_Vvdwtot], xmm5
	xorpd  xmm4, xmm4

	mulsd xmm7, [esp + nb030_tsc]
	mulsd xmm7, xmm0
	subsd  xmm4, xmm7

	movapd xmm0, [esp + nb030_dx]
	movapd xmm1, [esp + nb030_dy]
	movapd xmm2, [esp + nb030_dz]

	movd eax, mm0	

	mov    edi, [ebp + nb030_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb030_fix]
	movapd xmm4, [esp + nb030_fiy]
	movapd xmm5, [esp + nb030_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb030_fix], xmm3
	movlpd [esp + nb030_fiy], xmm4
	movlpd [esp + nb030_fiz], xmm5
	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5
.nb030_updateouterdata:
	mov   ecx, [esp + nb030_ii3]
	mov   edi, [ebp + nb030_faction]
	mov   esi, [ebp + nb030_fshift]
	mov   edx, [esp + nb030_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb030_fix]
	movapd xmm1, [esp + nb030_fiy]
	movapd xmm2, [esp + nb030_fiz]

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
	addsd xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movsd  [esi + edx*8],     xmm3
	movsd  [esi + edx*8 + 8], xmm4
	movsd  [esi + edx*8 + 16], xmm5

	;# get n from stack
	mov esi, [esp + nb030_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb030_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb030_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb030_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb030_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb030_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb030_n], esi
        jmp .nb030_outer
.nb030_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb030_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb030_end
        ;# non-zero, do one more workunit
        jmp   .nb030_threadloop
.nb030_end:
	emms

	mov eax, [esp + nb030_nouter]
	mov ebx, [esp + nb030_ninner]
	mov ecx, [ebp + nb030_outeriter]
	mov edx, [ebp + nb030_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb030_salign]
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



.globl nb_kernel030nf_ia32_sse2
.globl _nb_kernel030nf_ia32_sse2
nb_kernel030nf_ia32_sse2:	
_nb_kernel030nf_ia32_sse2:	
.equiv          nb030nf_p_nri,          8
.equiv          nb030nf_iinr,           12
.equiv          nb030nf_jindex,         16
.equiv          nb030nf_jjnr,           20
.equiv          nb030nf_shift,          24
.equiv          nb030nf_shiftvec,       28
.equiv          nb030nf_fshift,         32
.equiv          nb030nf_gid,            36
.equiv          nb030nf_pos,            40
.equiv          nb030nf_faction,        44
.equiv          nb030nf_charge,         48
.equiv          nb030nf_p_facel,        52
.equiv          nb030nf_argkrf,         56
.equiv          nb030nf_argcrf,         60
.equiv          nb030nf_Vc,             64
.equiv          nb030nf_type,           68
.equiv          nb030nf_p_ntype,        72
.equiv          nb030nf_vdwparam,       76
.equiv          nb030nf_Vvdw,           80
.equiv          nb030nf_p_tabscale,     84
.equiv          nb030nf_VFtab,          88
.equiv          nb030nf_invsqrta,       92
.equiv          nb030nf_dvda,           96
.equiv          nb030nf_p_gbtabscale,   100
.equiv          nb030nf_GBtab,          104
.equiv          nb030nf_p_nthreads,     108
.equiv          nb030nf_count,          112
.equiv          nb030nf_mtx,            116
.equiv          nb030nf_outeriter,      120
.equiv          nb030nf_inneriter,      124
.equiv          nb030nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb030nf_ix,             0
.equiv          nb030nf_iy,             16
.equiv          nb030nf_iz,             32
.equiv          nb030nf_tsc,            48
.equiv          nb030nf_c6,             64
.equiv          nb030nf_c12,            80
.equiv          nb030nf_Vvdwtot,        96
.equiv          nb030nf_half,           112
.equiv          nb030nf_three,          128
.equiv          nb030nf_is3,            144
.equiv          nb030nf_ii3,            148
.equiv          nb030nf_ntia,           152
.equiv          nb030nf_innerjjnr,      156
.equiv          nb030nf_innerk,         160
.equiv          nb030nf_n,              164
.equiv          nb030nf_nn1,            168
.equiv          nb030nf_nri,            172
.equiv          nb030nf_ntype,          176
.equiv          nb030nf_nouter,         180
.equiv          nb030nf_ninner,         184
.equiv          nb030nf_salign,         188
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
	mov [esp + nb030nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb030nf_p_nri]
	mov edi, [ebp + nb030nf_p_ntype]
	mov ecx, [ecx]
	mov edi, [edi]
	mov [esp + nb030nf_nri], ecx
	mov [esp + nb030nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb030nf_nouter], eax
	mov [esp + nb030nf_ninner], eax


	;# create constant floating-point factors on stack
	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb030nf_half], eax
	mov [esp + nb030nf_half+4], ebx
	movsd xmm1, [esp + nb030nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb030nf_half], xmm1
	movapd [esp + nb030nf_three], xmm3

	mov eax, [ebp + nb030nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb030nf_tsc], xmm3

.nb030nf_threadloop:
        mov   esi, [ebp + nb030nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb030nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb030nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb030nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb030nf_n], eax
        mov [esp + nb030nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb030nf_outerstart
        jmp .nb030nf_end

.nb030nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb030nf_nouter]
	mov [esp + nb030nf_nouter], ebx

.nb030nf_outer:
	mov   eax, [ebp + nb030nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax + esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb030nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb030nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

    	mov   edx, [ebp + nb030nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb030nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb030nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb030nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb030nf_ix], xmm0
	movapd [esp + nb030nf_iy], xmm1
	movapd [esp + nb030nf_iz], xmm2

	mov   [esp + nb030nf_ii3], ebx
	
	;# clear tot potential 
	xorpd xmm4, xmm4
	movapd [esp + nb030nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb030nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb030nf_pos]
	mov   eax, [ebp + nb030nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb030nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb030nf_ninner]
	mov   [esp + nb030nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb030nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb030nf_unroll_loop
	jmp   .nb030nf_checksingle
.nb030nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb030nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb030nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb030nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb030nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [esp + nb030nf_ntia]
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

	movapd [esp + nb030nf_c6], xmm4
	movapd [esp + nb030nf_c12], xmm6
	
	mov esi, [ebp + nb030nf_pos]		;# base of pos[] 
	lea   eax, [eax + eax*2]	;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]

	;# move nb030nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb030nf_ix]
	movapd xmm5, [esp + nb030nf_iy]
	movapd xmm6, [esp + nb030nf_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb030nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb030nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb030nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb030nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb030nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb030nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [esi + eax*8]	;# Y1 	
	movlpd xmm3, [esi + ebx*8]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movlpd xmm3, [esi + ebx*8 + 16]	;# G2
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 24]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	mulpd  xmm5, [esp + nb030nf_c6] ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [esp + nb030nf_Vvdwtot]
	movapd [esp + nb030nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1 	
	movlpd xmm3, [esi + ebx*8 + 32]	;# Y2 
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [esi + ebx*8 + 40]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movlpd xmm3, [esi + ebx*8 + 48]	;# G2
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [esi + ebx*8 + 56]	;# G2 H2 

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	mulpd  xmm5, [esp + nb030nf_c12]
	
	addpd  xmm5, [esp + nb030nf_Vvdwtot]
	movapd [esp + nb030nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [esp + nb030nf_innerk],  2
	jl    .nb030nf_checksingle
	jmp   .nb030nf_unroll_loop

.nb030nf_checksingle:				
	mov   edx, [esp + nb030nf_innerk]
	and   edx, 1
	jnz    .nb030nf_dosingle
	jmp    .nb030nf_updateouterdata
.nb030nf_dosingle:
	mov   edx, [esp + nb030nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	

	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov esi, [ebp + nb030nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb030nf_vdwparam]
	shl eax, 1	
	mov edi, [esp + nb030nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		

	movapd [esp + nb030nf_c6], xmm4
	movapd [esp + nb030nf_c12], xmm6
	
	mov esi, [ebp + nb030nf_pos]		;# base of pos[] 
	lea   eax, [eax + eax*2]	;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	;# move nb030nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb030nf_ix]
	movapd xmm5, [esp + nb030nf_iy]
	movapd xmm6, [esp + nb030nf_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	;# rsq in xmm4 

	cvtsd2ss xmm5, xmm4	
	rsqrtss xmm5, xmm5
	cvtss2sd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulsd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [esp + nb030nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb030nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb030nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb030nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb030nf_tsc]

	movd mm0, eax
	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	shl eax, 3	

	mov  esi, [ebp + nb030nf_VFtab]

	;# dispersion 
	movlpd xmm4, [esi + eax*8]	;# Y1
	movhpd xmm4, [esi + eax*8 + 8]	;# Y1 F1 
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 16]	;# G1
	movhpd xmm6, [esi + eax*8 + 24]	;# G1 H1 
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1  
	unpckhpd xmm7, xmm3	;# H1  
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [esp + nb030nf_c6];# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [esp + nb030nf_Vvdwtot]
	movlpd [esp + nb030nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [esi + eax*8 + 32]	;# Y1 
	movhpd xmm4, [esi + eax*8 + 40]	;# Y1 F1 
	xorpd xmm3,xmm3	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movlpd xmm6, [esi + eax*8 + 48]	;# G1
	movhpd xmm6, [esi + eax*8 + 56]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1  
	unpckhpd xmm7, xmm3	;# H1  
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	mulsd  xmm5, [esp + nb030nf_c12]
	
	addsd  xmm5, [esp + nb030nf_Vvdwtot]
	movlpd [esp + nb030nf_Vvdwtot], xmm5

.nb030nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb030nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb030nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb030nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb030nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb030nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb030nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb030nf_n], esi
        jmp .nb030nf_outer
.nb030nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb030nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb030nf_end
        ;# non-zero, do one more workunit
        jmp   .nb030nf_threadloop
.nb030nf_end:
	emms

	mov eax, [esp + nb030nf_nouter]
	mov ebx, [esp + nb030nf_ninner]
	mov ecx, [ebp + nb030nf_outeriter]
	mov edx, [ebp + nb030nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb030nf_salign]
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

	
