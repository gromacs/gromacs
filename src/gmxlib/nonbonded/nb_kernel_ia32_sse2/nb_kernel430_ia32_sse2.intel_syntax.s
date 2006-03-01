;#
;# $Id$
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
.equiv          .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as


.globl nb_kernel430_ia32_sse2
.globl _nb_kernel430_ia32_sse2
nb_kernel430_ia32_sse2:	
_nb_kernel430_ia32_sse2:	
.equiv          nb430_p_nri,            8
.equiv          nb430_iinr,             12
.equiv          nb430_jindex,           16
.equiv          nb430_jjnr,             20
.equiv          nb430_shift,            24
.equiv          nb430_shiftvec,         28
.equiv          nb430_fshift,           32
.equiv          nb430_gid,              36
.equiv          nb430_pos,              40
.equiv          nb430_faction,          44
.equiv          nb430_charge,           48
.equiv          nb430_p_facel,          52
.equiv          nb430_argkrf,           56
.equiv          nb430_argcrf,           60
.equiv          nb430_Vc,               64
.equiv          nb430_type,             68
.equiv          nb430_p_ntype,          72
.equiv          nb430_vdwparam,         76
.equiv          nb430_Vvdw,             80
.equiv          nb430_p_tabscale,       84
.equiv          nb430_VFtab,            88
.equiv          nb430_invsqrta,         92
.equiv          nb430_dvda,             96
.equiv          nb430_p_gbtabscale,     100
.equiv          nb430_GBtab,            104
.equiv          nb430_p_nthreads,       108
.equiv          nb430_count,            112
.equiv          nb430_mtx,              116
.equiv          nb430_outeriter,        120
.equiv          nb430_inneriter,        124
.equiv          nb430_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb430_ix,               0
.equiv          nb430_iy,               16
.equiv          nb430_iz,               32
.equiv          nb430_iq,               48
.equiv          nb430_dx,               64
.equiv          nb430_dy,               80
.equiv          nb430_dz,               96
.equiv          nb430_two,              112
.equiv          nb430_gbtsc,            128
.equiv          nb430_tsc,              144
.equiv          nb430_qq,               160
.equiv          nb430_c6,               176
.equiv          nb430_c12,              192
.equiv          nb430_fscal,            208
.equiv          nb430_vctot,            224
.equiv          nb430_Vvdwtot,          240
.equiv          nb430_fix,              256
.equiv          nb430_fiy,              272
.equiv          nb430_fiz,              288
.equiv          nb430_half,             304
.equiv          nb430_three,            320
.equiv          nb430_r,                336
.equiv          nb430_isai,             352
.equiv          nb430_isaprod,          368
.equiv          nb430_dvdasum,          384
.equiv          nb430_gbscale,          400
.equiv          nb430_ii,               416
.equiv          nb430_is3,              420
.equiv          nb430_ii3,              424
.equiv          nb430_ntia,             428
.equiv          nb430_innerjjnr,        432
.equiv          nb430_innerk,           436
.equiv          nb430_n,                440
.equiv          nb430_nn1,              444
.equiv          nb430_nri,              448
.equiv          nb430_facel,            456   ;# uses 8 bytes
.equiv          nb430_ntype,            464
.equiv          nb430_nouter,           468
.equiv          nb430_ninner,           472
.equiv          nb430_salign,           476
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 484		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb430_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb430_p_nri]
	mov esi, [ebp + nb430_p_facel]
	mov edi, [ebp + nb430_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb430_nri], ecx
	movsd [esp + nb430_facel], xmm7
	mov [esp + nb430_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb430_nouter], eax
	mov [esp + nb430_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb430_half], eax
	mov [esp + nb430_half+4], ebx
	movsd xmm1, [esp + nb430_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb430_half], xmm1
	movapd [esp + nb430_two], xmm2
	movapd [esp + nb430_three], xmm3
	mov eax, [ebp + nb430_p_tabscale]
	movsd xmm3, [eax]
	mov eax, [ebp + nb430_p_gbtabscale]
	movsd xmm4, [eax]
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb430_tsc], xmm3
	movapd [esp + nb430_gbtsc], xmm4

.nb430_threadloop:
        mov   esi, [ebp + nb430_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb430_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb430_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb430_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb430_n], eax
        mov [esp + nb430_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb430_outerstart
        jmp .nb430_end

.nb430_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb430_nouter]
	mov [esp + nb430_nouter], ebx

.nb430_outer:
	mov   eax, [ebp + nb430_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb430_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb430_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb430_iinr]       ;# ecx = pointer into iinr[]
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 
	mov   [esp + nb430_ii], ebx

	mov   edx, [ebp + nb430_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb430_facel]
	shufpd xmm3, xmm3, 0

	mov   edx, [ebp + nb430_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [edx + ebx*8]
	shufpd xmm4, xmm4, 0

    	mov   edx, [ebp + nb430_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb430_ntype]
    	shl   edx, 1
    	mov   [esp + nb430_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb430_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb430_iq], xmm3
	movapd [esp + nb430_isai], xmm4
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb430_ix], xmm0
	movapd [esp + nb430_iy], xmm1
	movapd [esp + nb430_iz], xmm2

	mov   [esp + nb430_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb430_vctot], xmm4
	movapd [esp + nb430_Vvdwtot], xmm4
	movapd [esp + nb430_dvdasum], xmm4
	movapd [esp + nb430_fix], xmm4
	movapd [esp + nb430_fiy], xmm4
	movapd [esp + nb430_fiz], xmm4
	
	mov   eax, [ebp + nb430_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb430_pos]
	mov   edi, [ebp + nb430_faction]	
	mov   eax, [ebp + nb430_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb430_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb430_ninner]
	mov   [esp + nb430_ninner], ecx
	add   edx, 0
	mov   [esp + nb430_innerk], edx    ;# number of innerloop atoms 
	jge   .nb430_unroll_loop
	jmp   .nb430_checksingle
.nb430_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb430_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [edx]
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb430_innerjjnr], 8	;# advance pointer (unrolled 2) 

	;# load isaj
	mov esi, [ebp + nb430_invsqrta]
	movlpd xmm2, [esi + eax*8]
	movhpd xmm2, [esi + ebx*8]
	mulpd  xmm2, [esp + nb430_isai]
	movapd [esp + nb430_isaprod], xmm2	
	movapd xmm1, xmm2
	mulpd xmm1, [esp + nb430_gbtsc]
	movapd [esp + nb430_gbscale], xmm1
	
	mov esi, [ebp + nb430_charge]    ;# base of charge[] 
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	mulpd xmm2, [esp + nb430_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb430_qq], xmm3	
	
	mov esi, [ebp + nb430_type]
	mov ecx, [esi + eax*4]
	mov edx, [esi + ebx*4]
	mov esi, [ebp + nb430_vdwparam]
	shl ecx, 1
	shl edx, 1
	mov edi, [esp + nb430_ntia]
	add ecx, edi
	add edx, edi

	movlpd xmm6, [esi + ecx*8]	;# c6a
	movlpd xmm7, [esi + edx*8]	;# c6b
	movhpd xmm6, [esi + ecx*8 + 8]	;# c6a c12a 
	movhpd xmm7, [esi + edx*8 + 8]	;# c6b c12b 

	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [esp + nb430_c6], xmm4
	movapd [esp + nb430_c12], xmm6
	
	mov esi, [ebp + nb430_pos]		;# base of pos[] 

	movd  mm2, eax
	movd  mm3, ebx
	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	mov    edi, [ebp + nb430_faction]
	
	;# move nb430_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb430_ix]
	movapd xmm5, [esp + nb430_iy]
	movapd xmm6, [esp + nb430_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb430_dx], xmm4
	movapd [esp + nb430_dy], xmm5
	movapd [esp + nb430_dz], xmm6
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
	movapd xmm1, [esp + nb430_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb430_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb430_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb430_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv 
	mulpd xmm4, xmm0	;# xmm4=r 
	movapd [esp + nb430_r], xmm4
	mulpd xmm4, [esp + nb430_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 

	mov  esi, [ebp + nb430_GBtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6		;# indices in eax/ebx 

	;# Coulomb 
	movapd xmm4, [esi + ecx*8]	;# Y1 F1 	
	movapd xmm3, [esi + edx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [esi + ecx*8 + 16]	;# G1 H1 	
	movapd xmm3, [esi + edx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb430_two]	;# two*Heps2 
	movapd xmm3, [esp + nb430_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
	;# get jnr from regs
	movd ecx, mm2
	movd edx, mm3
	mov esi, [ebp + nb430_dvda]
	
	;# Calculate dVda
	xorpd xmm7, xmm7
	mulpd xmm3, [esp + nb430_gbscale]
	movapd xmm6, xmm3
	mulpd  xmm6, [esp + nb430_r]
	addpd  xmm6, xmm5
	addpd  xmm5, [esp + nb430_vctot]
	movapd [esp + nb430_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subpd  xmm7, xmm6
	movapd xmm6, xmm7
	
	;# update dvdasum
	addpd  xmm7, [esp + nb430_dvdasum]
	movapd [esp + nb430_dvdasum], xmm7 

	;# update j atoms dvdaj
	movhlps xmm7, xmm6
	addsd  xmm6, [esi + ecx*8]
	addsd  xmm7, [esi + edx*8]
	movsd  [esi + ecx*8], xmm6
	movsd  [esi + edx*8], xmm7
	
	;# put scalar force on stack temporarily 
	movapd [esp + nb430_fscal], xmm3

	movapd xmm4, [esp + nb430_r]
	mulpd  xmm4, [esp + nb430_tsc]
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 3		;# idx *= 8

	mov  esi, [ebp + nb430_VFtab]

	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6		;# indices in eax/ebx 

	;# Dispersion 
	movapd xmm4, [esi + ecx*8]	;# Y1 F1 	
	movapd xmm3, [esi + edx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [esi + ecx*8 + 16]	;# G1 H1 	
	movapd xmm3, [esi + edx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb430_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb430_c6]
	mulpd  xmm7, xmm4	 ;# fijD 
	mulpd  xmm5, xmm4	 ;# Vvdw6
	mulpd  xmm7, [esp + nb430_tsc]
	addpd  xmm7, [esp + nb430_fscal] ;# add to fscal 

	;# put scalar force back on stack Update Vvdwtot directly 
	addpd  xmm5, [esp + nb430_Vvdwtot]
	movapd [esp + nb430_fscal], xmm7
	movapd [esp + nb430_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [esi + ecx*8 + 32]	;# Y1 F1 	
	movapd xmm3, [esi + edx*8 + 32]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [esi + ecx*8 + 48]	;# G1 H1 	
	movapd xmm3, [esi + edx*8 + 48]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb430_two]	;# two*Heps2 
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb430_c12]
	mulpd  xmm7, xmm4 ;# fijR 
	mulpd  xmm5, xmm4 ;# Vvdw12 
	mulpd  xmm7, [esp + nb430_tsc]
	addpd  xmm7, [esp + nb430_fscal] 
	
	addpd  xmm5, [esp + nb430_Vvdwtot]
	movapd [esp + nb430_Vvdwtot], xmm5
	xorpd  xmm4, xmm4

	mulpd xmm7, xmm0
	subpd xmm4, xmm7

	movapd xmm0, [esp + nb430_dx]
	movapd xmm1, [esp + nb430_dy]
	movapd xmm2, [esp + nb430_dz]

	mov    edi, [ebp + nb430_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb430_fix]
	movapd xmm4, [esp + nb430_fiy]
	movapd xmm5, [esp + nb430_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb430_fix], xmm3
	movapd [esp + nb430_fiy], xmm4
	movapd [esp + nb430_fiz], xmm5
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
	sub dword ptr [esp + nb430_innerk],  2
	jl    .nb430_checksingle
	jmp   .nb430_unroll_loop
.nb430_checksingle:
	mov   edx, [esp + nb430_innerk]
	and   edx, 1
	jnz    .nb430_dosingle
	jmp    .nb430_updateouterdata
.nb430_dosingle:
	mov esi, [ebp + nb430_charge]
	mov edx, [ebp + nb430_invsqrta]
	mov edi, [ebp + nb430_pos]
	mov   ecx, [esp + nb430_innerjjnr]
	mov   eax, [ecx]	

	xorpd  xmm6, xmm6
	movapd xmm7, xmm6
	movsd  xmm7, [edx + eax*8]
	movlpd xmm6, [esi + eax*8]	;# xmm6(0) has the charge
	mulsd  xmm7, [esp + nb430_isai]
	movapd [esp + nb430_isaprod], xmm7
	movapd xmm1, xmm7
	mulpd xmm1, [esp + nb430_gbtsc]
	movapd [esp + nb430_gbscale], xmm1
	
	mulsd  xmm7, [esp + nb430_iq]
	mulsd  xmm6, xmm7
	movapd [esp + nb430_qq], xmm6
	
	mov esi, [ebp + nb430_type]
	mov edx, [esi + eax*4]
	mov esi, [ebp + nb430_vdwparam]
	shl edx, 1
	mov edi, [esp + nb430_ntia]
	add edx, edi

	movlpd xmm6, [esi + edx*8]	;# c6a
	movhpd xmm6, [esi + edx*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [esp + nb430_c6], xmm4
	movapd [esp + nb430_c12], xmm6
	
	mov esi, [ebp + nb430_pos]		;# base of pos[]
	
	movd  mm2, eax
	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	mov    edi, [ebp + nb430_faction]

	;# move nb430_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb430_ix]
	movapd xmm5, [esp + nb430_iy]
	movapd xmm6, [esp + nb430_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb430_dx], xmm4
	movapd [esp + nb430_dy], xmm5
	movapd [esp + nb430_dz], xmm6
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
	movapd xmm1, [esp + nb430_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb430_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb430_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb430_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
	movsd [esp + nb430_r], xmm4
	mulsd xmm4, [esp + nb430_gbscale]
	
	cvttsd2si edx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, edx
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl edx, 2		;# idx *= 4 
	mov  esi, [ebp + nb430_GBtab]

	;# Coulomb 
	movapd xmm4, [esi + edx*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [esi + edx*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb430_two]	;# two*Heps2 
	movapd xmm3, [esp + nb430_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
	;# get jnr from regs
	movd ebx, mm2
	mov esi, [ebp + nb430_dvda]
	
	;# Calculate dVda
	xorpd xmm7, xmm7
	mulsd xmm3, [esp + nb430_gbscale]
	movsd xmm6, xmm3
	mulsd  xmm6, [esp + nb430_r]
	addsd  xmm6, xmm5
	addsd  xmm5, [esp + nb430_vctot]
	movsd [esp + nb430_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subpd xmm7, xmm6
	movsd xmm6, xmm7
	
	;# update dvdasum
	addsd  xmm7, [esp + nb430_dvdasum]
	movsd [esp + nb430_dvdasum], xmm7 

	;# update j atoms dvdaj
	addsd  xmm6, [esi + ebx*8]
	movsd  [esi + ebx*8], xmm6
	
	;# put scalar force on stack temporarily 
	movsd [esp + nb430_fscal], xmm3

	movsd xmm4, [esp + nb430_r]
	mulsd  xmm4, [esp + nb430_tsc]
	cvttsd2si edx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, edx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl edx, 3

	mov  esi, [ebp + nb430_VFtab]

	;# Dispersion 
	movapd xmm4, [esi + edx*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [esi + edx*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb430_two]	;# two*Heps2 
	movapd xmm3, [esp + nb430_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb430_c6]
	mulsd  xmm7, xmm4	 ;# fijD 
	mulsd  xmm5, xmm4	 ;# Vvdw6
	mulpd  xmm7, [esp + nb430_tsc]
	addsd  xmm7, [esp + nb430_fscal] ;# add to fscal 

	;# put scalar force back on stack Update Vvdwtot directly 
	addsd  xmm5, [esp + nb430_Vvdwtot]
	movlpd [esp + nb430_fscal], xmm7
	movlpd [esp + nb430_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [esi + edx*8 + 32]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [esi + edx*8 + 48]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb430_two]	;# two*Heps2 
	movapd xmm3, [esp + nb430_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movapd xmm4, [esp + nb430_c12]
	mulsd  xmm7, xmm4 ;# fijR 
	mulsd  xmm5, xmm4 ;# Vvdw12 
	mulpd  xmm7, [esp + nb430_tsc]
	addsd  xmm7, [esp + nb430_fscal] 
	
	addsd  xmm5, [esp + nb430_Vvdwtot]
	movlpd [esp + nb430_Vvdwtot], xmm5
	xorpd  xmm4, xmm4

	mulsd xmm7, xmm0
	subsd xmm4, xmm7

	movapd xmm0, [esp + nb430_dx]
	movapd xmm1, [esp + nb430_dy]
	movapd xmm2, [esp + nb430_dz]

	mov    edi, [ebp + nb430_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb430_fix]
	movapd xmm4, [esp + nb430_fiy]
	movapd xmm5, [esp + nb430_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb430_fix], xmm3
	movlpd [esp + nb430_fiy], xmm4
	movlpd [esp + nb430_fiz], xmm5
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
.nb430_updateouterdata:
	mov   ecx, [esp + nb430_ii3]
	mov   edi, [ebp + nb430_faction]
	mov   esi, [ebp + nb430_fshift]
	mov   edx, [esp + nb430_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb430_fix]
	movapd xmm1, [esp + nb430_fiy]
	movapd xmm2, [esp + nb430_fiz]

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
	mov esi, [esp + nb430_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb430_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb430_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb430_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb430_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb430_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate dVda and update it 
	movapd xmm7, [esp + nb430_dvdasum]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	mov edx, [esp + nb430_ii]
	mov eax, [ebp + nb430_dvda]
	addsd xmm7, [eax + edx*8]
	movsd [eax + edx*8], xmm7
	
        ;# finish if last 
        mov ecx, [esp + nb430_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb430_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb430_n], esi
        jmp .nb430_outer
.nb430_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb430_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb430_end
        ;# non-zero, do one more workunit
        jmp   .nb430_threadloop
.nb430_end:
	emms

	mov eax, [esp + nb430_nouter]
	mov ebx, [esp + nb430_ninner]
	mov ecx, [ebp + nb430_outeriter]
	mov edx, [ebp + nb430_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb430_salign]
	add esp, eax
	add esp, 484
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




	
.globl nb_kernel430nf_ia32_sse2
.globl _nb_kernel430nf_ia32_sse2
nb_kernel430nf_ia32_sse2:	
_nb_kernel430nf_ia32_sse2:	
.equiv          nb430nf_p_nri,          8
.equiv          nb430nf_iinr,           12
.equiv          nb430nf_jindex,         16
.equiv          nb430nf_jjnr,           20
.equiv          nb430nf_shift,          24
.equiv          nb430nf_shiftvec,       28
.equiv          nb430nf_fshift,         32
.equiv          nb430nf_gid,            36
.equiv          nb430nf_pos,            40
.equiv          nb430nf_faction,        44
.equiv          nb430nf_charge,         48
.equiv          nb430nf_p_facel,        52
.equiv          nb430nf_argkrf,         56
.equiv          nb430nf_argcrf,         60
.equiv          nb430nf_Vc,             64
.equiv          nb430nf_type,           68
.equiv          nb430nf_p_ntype,        72
.equiv          nb430nf_vdwparam,       76
.equiv          nb430nf_Vvdw,           80
.equiv          nb430nf_p_tabscale,     84
.equiv          nb430nf_VFtab,          88
.equiv          nb430nf_invsqrta,       92
.equiv          nb430nf_dvda,           96
.equiv          nb430nf_p_gbtabscale,   100
.equiv          nb430nf_GBtab,          104
.equiv          nb430nf_p_nthreads,     108
.equiv          nb430nf_count,          112
.equiv          nb430nf_mtx,            116
.equiv          nb430nf_outeriter,      120
.equiv          nb430nf_inneriter,      124
.equiv          nb430nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb430nf_ix,             0
.equiv          nb430nf_iy,             16
.equiv          nb430nf_iz,             32
.equiv          nb430nf_iq,             48
.equiv          nb430nf_gbtsc,          64
.equiv          nb430nf_tsc,            80
.equiv          nb430nf_qq,             96
.equiv          nb430nf_c6,             112
.equiv          nb430nf_c12,            128
.equiv          nb430nf_vctot,          144
.equiv          nb430nf_Vvdwtot,        160
.equiv          nb430nf_half,           176
.equiv          nb430nf_three,          192
.equiv          nb430nf_r,              208
.equiv          nb430nf_isai,           224
.equiv          nb430nf_isaprod,        240
.equiv          nb430nf_gbscale,        256
.equiv          nb430nf_is3,            272
.equiv          nb430nf_ii3,            276
.equiv          nb430nf_ntia,           280
.equiv          nb430nf_innerjjnr,      284
.equiv          nb430nf_innerk,         288
.equiv          nb430nf_n,              292
.equiv          nb430nf_nn1,            296
.equiv          nb430nf_nri,            300
.equiv          nb430nf_facel,          304   ;# uses 8 bytes
.equiv          nb430nf_ntype,          312
.equiv          nb430nf_nouter,         316
.equiv          nb430nf_ninner,         320
.equiv          nb430nf_salign,         324
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 328		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb430nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb430nf_p_nri]
	mov esi, [ebp + nb430nf_p_facel]
	mov edi, [ebp + nb430nf_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb430nf_nri], ecx
	movsd [esp + nb430nf_facel], xmm7
	mov [esp + nb430nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb430nf_nouter], eax
	mov [esp + nb430nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb430nf_half], eax
	mov [esp + nb430nf_half+4], ebx
	movsd xmm1, [esp + nb430nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb430nf_half], xmm1
	movapd [esp + nb430nf_three], xmm3
	mov eax, [ebp + nb430nf_p_tabscale]
	movsd xmm3, [eax]
	mov eax, [ebp + nb430nf_p_gbtabscale]
	movsd xmm4, [eax]
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [esp + nb430nf_tsc], xmm3
	movapd [esp + nb430nf_gbtsc], xmm4

.nb430nf_threadloop:
        mov   esi, [ebp + nb430nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb430nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb430nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb430nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb430nf_n], eax
        mov [esp + nb430nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb430nf_outerstart
        jmp .nb430nf_end

.nb430nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb430nf_nouter]
	mov [esp + nb430nf_nouter], ebx

.nb430nf_outer:
	mov   eax, [ebp + nb430nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb430nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb430nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb430nf_iinr]       ;# ecx = pointer into iinr[]
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb430nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb430nf_facel]
	shufpd xmm3, xmm3, 0

	mov   edx, [ebp + nb430nf_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [edx + ebx*8]
	shufpd xmm4, xmm4, 0

    	mov   edx, [ebp + nb430nf_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb430nf_ntype]
    	shl   edx, 1
    	mov   [esp + nb430nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb430nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb430nf_iq], xmm3
	movapd [esp + nb430nf_isai], xmm4	
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb430nf_ix], xmm0
	movapd [esp + nb430nf_iy], xmm1
	movapd [esp + nb430nf_iz], xmm2

	mov   [esp + nb430nf_ii3], ebx
	
	;# clear vctot
	xorpd xmm4, xmm4
	movapd [esp + nb430nf_vctot], xmm4
	movapd [esp + nb430nf_Vvdwtot], xmm4

	mov   eax, [ebp + nb430nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb430nf_pos]
	mov   edi, [ebp + nb430nf_faction]	
	mov   eax, [ebp + nb430nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb430nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb430nf_ninner]
	mov   [esp + nb430nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb430nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb430nf_unroll_loop
	jmp   .nb430nf_checksingle
.nb430nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb430nf_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [edx]
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb430nf_innerjjnr], 8	;# advance pointer (unrolled 2) 

	;# load isaj
	mov esi, [ebp + nb430nf_invsqrta]
	movlpd xmm2, [esi + eax*8]
	movhpd xmm2, [esi + ebx*8]
	mulpd  xmm2, [esp + nb430nf_isai]
	movapd [esp + nb430nf_isaprod], xmm2	
	movapd xmm1, xmm2
	mulpd xmm1, [esp + nb430nf_gbtsc]
	movapd [esp + nb430nf_gbscale], xmm1
	
	mov esi, [ebp + nb430nf_charge]    ;# base of charge[] 
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	mulpd xmm2, [esp + nb430nf_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb430nf_qq], xmm3	
	
	mov esi, [ebp + nb430nf_type]
	mov ecx, [esi + eax*4]
	mov edx, [esi + ebx*4]
	mov esi, [ebp + nb430nf_vdwparam]
	shl ecx, 1
	shl edx, 1
	mov edi, [esp + nb430nf_ntia]
	add ecx, edi
	add edx, edi

	movlpd xmm6, [esi + ecx*8]	;# c6a
	movlpd xmm7, [esi + edx*8]	;# c6b
	movhpd xmm6, [esi + ecx*8 + 8]	;# c6a c12a 
	movhpd xmm7, [esi + edx*8 + 8]	;# c6b c12b 

	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [esp + nb430nf_c6], xmm4
	movapd [esp + nb430nf_c12], xmm6
	
	mov esi, [ebp + nb430nf_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	mov    edi, [ebp + nb430nf_faction]
	
	;# move nb430nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb430nf_ix]
	movapd xmm5, [esp + nb430nf_iy]
	movapd xmm6, [esp + nb430nf_iz]

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
	movapd xmm1, [esp + nb430nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb430nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb430nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb430nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv 
	mulpd xmm4, xmm0	;# xmm4=r 
	movapd [esp + nb430nf_r], xmm4
	mulpd xmm4, [esp + nb430nf_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 

	mov  esi, [ebp + nb430nf_GBtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6		;# indices in eax/ebx 

	;# Coulomb 
	movapd xmm4, [esi + ecx*8]	;# Y1 F1 	
	movapd xmm3, [esi + edx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [esi + ecx*8 + 16]	;# G1 H1 	
	movapd xmm3, [esi + edx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb430nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	addpd  xmm5, [esp + nb430nf_vctot]
	movapd [esp + nb430nf_vctot], xmm5
	
	movapd xmm4, [esp + nb430nf_r]
	mulpd  xmm4, [esp + nb430nf_tsc]
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 3		;# idx *= 8

	mov  esi, [ebp + nb430nf_VFtab]

	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6		;# indices in eax/ebx 

	;# Dispersion 
	movapd xmm4, [esi + ecx*8]	;# Y1 F1 	
	movapd xmm3, [esi + edx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [esi + ecx*8 + 16]	;# G1 H1 	
	movapd xmm3, [esi + edx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	mulpd  xmm5, [esp + nb430nf_c6]	 ;# Vvdw6
	addpd  xmm5, [esp + nb430nf_Vvdwtot]
	movapd [esp + nb430nf_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [esi + ecx*8 + 32]	;# Y1 F1 	
	movapd xmm3, [esi + edx*8 + 32]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [esi + ecx*8 + 48]	;# G1 H1 	
	movapd xmm3, [esi + edx*8 + 48]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	mulpd  xmm5, [esp + nb430nf_c12] ;# Vvdw12 
	addpd  xmm5, [esp + nb430nf_Vvdwtot]
	movapd [esp + nb430nf_Vvdwtot], xmm5
	xorpd  xmm4, xmm4
	
	;# should we do one more iteration? 
	sub dword ptr [esp + nb430nf_innerk],  2
	jl    .nb430nf_checksingle
	jmp   .nb430nf_unroll_loop
.nb430nf_checksingle:
	mov   edx, [esp + nb430nf_innerk]
	and   edx, 1
	jnz    .nb430nf_dosingle
	jmp    .nb430nf_updateouterdata
.nb430nf_dosingle:
	mov esi, [ebp + nb430nf_charge]
	mov edx, [ebp + nb430nf_invsqrta]
	mov edi, [ebp + nb430nf_pos]
	mov   ecx, [esp + nb430nf_innerjjnr]
	mov   eax, [ecx]	

	xorpd  xmm6, xmm6
	movapd xmm7, xmm6
	movsd  xmm7, [edx + eax*8]
	movlpd xmm6, [esi + eax*8]	;# xmm6(0) has the charge
	mulsd  xmm7, [esp + nb430nf_isai]
	movapd [esp + nb430nf_isaprod], xmm7
	movapd xmm1, xmm7
	mulpd xmm1, [esp + nb430nf_gbtsc]
	movapd [esp + nb430nf_gbscale], xmm1
	
	mulsd  xmm7, [esp + nb430nf_iq]
	mulsd  xmm6, xmm7
	movapd [esp + nb430nf_qq], xmm6
	
	mov esi, [ebp + nb430nf_type]
	mov edx, [esi + eax*4]
	mov esi, [ebp + nb430nf_vdwparam]
	shl edx, 1
	mov edi, [esp + nb430nf_ntia]
	add edx, edi

	movlpd xmm6, [esi + edx*8]	;# c6a
	movhpd xmm6, [esi + edx*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [esp + nb430nf_c6], xmm4
	movapd [esp + nb430nf_c12], xmm6
	
	mov esi, [ebp + nb430nf_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]

	mov    edi, [ebp + nb430nf_faction]

	;# move nb430nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb430nf_ix]
	movapd xmm5, [esp + nb430nf_iy]
	movapd xmm6, [esp + nb430nf_iz]

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
	movapd xmm1, [esp + nb430nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb430nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb430nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb430nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
	movsd [esp + nb430nf_r], xmm4
	mulsd xmm4, [esp + nb430nf_gbscale]
	
	cvttsd2si edx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, edx
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl edx, 2		;# idx *= 4 
	mov  esi, [ebp + nb430nf_GBtab]

	;# Coulomb 
	movapd xmm4, [esi + edx*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [esi + edx*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb430nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	addsd  xmm5, [esp + nb430nf_vctot]
	movsd [esp + nb430nf_vctot], xmm5 

	movsd xmm4, [esp + nb430nf_r]
	mulsd  xmm4, [esp + nb430nf_tsc]
	cvttsd2si edx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, edx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2

	shl edx, 3

	mov  esi, [ebp + nb430nf_VFtab]

	;# Dispersion 
	movapd xmm4, [esi + edx*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [esi + edx*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [esp + nb430nf_c6]	 ;# Vvdw6
	addsd  xmm5, [esp + nb430nf_Vvdwtot]
	movlpd [esp + nb430nf_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [esi + edx*8 + 32]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [esi + edx*8 + 48]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, [esp + nb430nf_c12] ;# Vvdw12 
	addsd  xmm5, [esp + nb430nf_Vvdwtot]
	movlpd [esp + nb430nf_Vvdwtot], xmm5
.nb430nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb430nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb430nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb430nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb430nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb430nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb430nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
		
        ;# finish if last 
        mov ecx, [esp + nb430nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb430nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb430nf_n], esi
        jmp .nb430nf_outer
.nb430nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb430nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb430nf_end
        ;# non-zero, do one more workunit
        jmp   .nb430nf_threadloop
.nb430nf_end:
	emms

	mov eax, [esp + nb430nf_nouter]
	mov ebx, [esp + nb430nf_ninner]
	mov ecx, [ebp + nb430nf_outeriter]
	mov edx, [ebp + nb430nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb430nf_salign]
	add esp, eax
	add esp, 328
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret

