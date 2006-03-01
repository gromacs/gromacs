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





.globl nb_kernel410_ia32_sse2
.globl _nb_kernel410_ia32_sse2
nb_kernel410_ia32_sse2:	
_nb_kernel410_ia32_sse2:	
.equiv          nb410_p_nri,            8
.equiv          nb410_iinr,             12
.equiv          nb410_jindex,           16
.equiv          nb410_jjnr,             20
.equiv          nb410_shift,            24
.equiv          nb410_shiftvec,         28
.equiv          nb410_fshift,           32
.equiv          nb410_gid,              36
.equiv          nb410_pos,              40
.equiv          nb410_faction,          44
.equiv          nb410_charge,           48
.equiv          nb410_p_facel,          52
.equiv          nb410_argkrf,           56
.equiv          nb410_argcrf,           60
.equiv          nb410_Vc,               64
.equiv          nb410_type,             68
.equiv          nb410_p_ntype,          72
.equiv          nb410_vdwparam,         76
.equiv          nb410_Vvdw,             80
.equiv          nb410_p_tabscale,       84
.equiv          nb410_VFtab,            88
.equiv          nb410_invsqrta,         92
.equiv          nb410_dvda,             96
.equiv          nb410_p_gbtabscale,     100
.equiv          nb410_GBtab,            104
.equiv          nb410_p_nthreads,       108
.equiv          nb410_count,            112
.equiv          nb410_mtx,              116
.equiv          nb410_outeriter,        120
.equiv          nb410_inneriter,        124
.equiv          nb410_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb410_ix,               0
.equiv          nb410_iy,               16
.equiv          nb410_iz,               32
.equiv          nb410_iq,               48
.equiv          nb410_dx,               64
.equiv          nb410_dy,               80
.equiv          nb410_dz,               96
.equiv          nb410_two,              112
.equiv          nb410_six,              128
.equiv          nb410_twelve,           144
.equiv          nb410_gbtsc,            160
.equiv          nb410_qq,               176
.equiv          nb410_c6,               192
.equiv          nb410_c12,              208
.equiv          nb410_fscal,            224
.equiv          nb410_vctot,            240
.equiv          nb410_Vvdwtot,          256
.equiv          nb410_fix,              272
.equiv          nb410_fiy,              288
.equiv          nb410_fiz,              304
.equiv          nb410_half,             320
.equiv          nb410_three,            336
.equiv          nb410_r,                352
.equiv          nb410_isai,             368
.equiv          nb410_isaprod,          384
.equiv          nb410_dvdasum,          400
.equiv          nb410_gbscale,          416
.equiv          nb410_ii,               432
.equiv          nb410_is3,              436
.equiv          nb410_ii3,              440
.equiv          nb410_ntia,             444
.equiv          nb410_innerjjnr,        448
.equiv          nb410_innerk,           452
.equiv          nb410_n,                456
.equiv          nb410_nn1,              460
.equiv          nb410_nri,              464
.equiv          nb410_facel,            472   ;# uses 8 bytes
.equiv          nb410_ntype,            480
.equiv          nb410_nouter,           484
.equiv          nb410_ninner,           488
.equiv          nb410_salign,           492
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 496		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb410_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb410_p_nri]
	mov esi, [ebp + nb410_p_facel]
	mov edi, [ebp + nb410_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb410_nri], ecx
	movsd [esp + nb410_facel], xmm7
	mov [esp + nb410_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb410_nouter], eax
	mov [esp + nb410_ninner], eax


	mov eax, [ebp + nb410_p_gbtabscale]
	movsd xmm5, [eax]
	shufpd xmm5, xmm5, 0
	movapd [esp + nb410_gbtsc], xmm5
	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb410_half], eax
	mov [esp + nb410_half+4], ebx
	movsd xmm1, [esp + nb410_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# 6.0
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# 12.0
	movapd [esp + nb410_half], xmm1
	movapd [esp + nb410_two], xmm2
	movapd [esp + nb410_three], xmm3
	movapd [esp + nb410_six], xmm4
	movapd [esp + nb410_twelve], xmm5

.nb410_threadloop:
        mov   esi, [ebp + nb410_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb410_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb410_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb410_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb410_n], eax
        mov [esp + nb410_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb410_outerstart
        jmp .nb410_end

.nb410_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb410_nouter]
	mov [esp + nb410_nouter], ebx

.nb410_outer:
	mov   eax, [ebp + nb410_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb410_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb410_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb410_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 
	mov   [esp + nb410_ii], ebx

	mov   edx, [ebp + nb410_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb410_facel]
	shufpd xmm3, xmm3, 0

	mov   edx, [ebp + nb410_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [edx + ebx*8]
	shufpd xmm4, xmm4, 0

    	mov   edx, [ebp + nb410_type] 
    	mov   edx, [edx + ebx*4]
    	imul  edx, [esp + nb410_ntype]
    	shl   edx, 1
    	mov   [esp + nb410_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb410_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb410_iq], xmm3
	movapd [esp + nb410_isai], xmm4

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb410_ix], xmm0
	movapd [esp + nb410_iy], xmm1
	movapd [esp + nb410_iz], xmm2

	mov   [esp + nb410_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb410_vctot], xmm4
	movapd [esp + nb410_Vvdwtot], xmm4
	movapd [esp + nb410_dvdasum], xmm4
	movapd [esp + nb410_fix], xmm4
	movapd [esp + nb410_fiy], xmm4
	movapd [esp + nb410_fiz], xmm4
	
	mov   eax, [ebp + nb410_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb410_pos]
	mov   edi, [ebp + nb410_faction]	
	mov   eax, [ebp + nb410_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb410_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb410_ninner]
	mov   [esp + nb410_ninner], ecx
	add   edx, 0
	mov   [esp + nb410_innerk], edx    ;# number of innerloop atoms 
	jge   .nb410_unroll_loop
	jmp   .nb410_checksingle
.nb410_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb410_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb410_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	;# load isaj
	mov esi, [ebp + nb410_invsqrta]
	movlpd xmm2, [esi + eax*8]
	movhpd xmm2, [esi + ebx*8]
	mulpd  xmm2, [esp + nb410_isai]
	movapd [esp + nb410_isaprod], xmm2	
	movapd xmm1, xmm2
	mulpd xmm1, [esp + nb410_gbtsc]
	movapd [esp + nb410_gbscale], xmm1
	
	mov esi, [ebp + nb410_charge]    ;# base of charge[] 
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	mulpd xmm2, [esp + nb410_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb410_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb410_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb410_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb410_ntia]
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
	movapd [esp + nb410_c6], xmm4
	movapd [esp + nb410_c12], xmm6
	
	mov esi, [ebp + nb410_pos]       ;# base of pos[] 

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
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb410_ix]
	movapd xmm5, [esp + nb410_iy]
	movapd xmm6, [esp + nb410_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb410_dx], xmm4
	movapd [esp + nb410_dy], xmm5
	movapd [esp + nb410_dz], xmm6
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
	movapd xmm1, [esp + nb410_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb410_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb410_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb410_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	movapd [esp + nb410_r], xmm4
	mulpd xmm4, [esp + nb410_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb410_GBtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [esi + eax*8]	;# Y1 F1 	
	movapd xmm3, [esi + ebx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [esi + eax*8 + 16]	;# G1 H1 	
	movapd xmm3, [esi + ebx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb410_two]	;# two*Heps2 
	movapd xmm3, [esp + nb410_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
	;# get jnr from regs
	movd ecx, mm2
	movd edx, mm3
	mov esi, [ebp + nb410_dvda]
	
	;# Calculate dVda
	xorpd xmm7, xmm7
	mulpd xmm3, [esp + nb410_gbscale]
	movapd xmm6, xmm3
	mulpd  xmm6, [esp + nb410_r]
	addpd  xmm6, xmm5
	addpd  xmm5, [esp + nb410_vctot]
	movapd [esp + nb410_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subpd  xmm7, xmm6
	movapd xmm6, xmm7
	
	;# update dvdasum
	addpd  xmm7, [esp + nb410_dvdasum]
	movapd [esp + nb410_dvdasum], xmm7 

	;# update j atoms dvdaj
	movhlps xmm7, xmm6
	addsd  xmm6, [esi + ecx*8]
	addsd  xmm7, [esi + edx*8]
	movsd  [esi + ecx*8], xmm6
	movsd  [esi + edx*8], xmm7
	
	;# L-J 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm0	;# xmm4=rinvsq 

	movapd xmm6, xmm4
	mulpd  xmm6, xmm4
	
	mulpd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulpd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulpd  xmm6, [esp + nb410_c6]
	mulpd  xmm4, [esp + nb410_c12]
	movapd xmm7, [esp + nb410_Vvdwtot]
	addpd  xmm7, xmm4
	mulpd  xmm4, [esp + nb410_twelve]
	subpd  xmm7, xmm6
	mulpd  xmm6, [esp + nb410_six]
	movapd [esp + nb410_Vvdwtot], xmm7
	subpd  xmm4, xmm6
	mulpd  xmm4, xmm0
	subpd  xmm4, xmm3
	mulpd  xmm4, xmm0

	movapd xmm0, [esp + nb410_dx]
	movapd xmm1, [esp + nb410_dy]
	movapd xmm2, [esp + nb410_dz]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb410_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb410_fix]
	movapd xmm4, [esp + nb410_fiy]
	movapd xmm5, [esp + nb410_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb410_fix], xmm3
	movapd [esp + nb410_fiy], xmm4
	movapd [esp + nb410_fiz], xmm5
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
	sub dword ptr [esp + nb410_innerk],  2
	jl    .nb410_checksingle
	jmp   .nb410_unroll_loop
.nb410_checksingle:
	mov   edx, [esp + nb410_innerk]
	and   edx, 1
	jnz    .nb410_dosingle
	jmp    .nb410_updateouterdata
.nb410_dosingle:
	mov esi, [ebp + nb410_charge]
	mov edx, [ebp + nb410_invsqrta]
	mov edi, [ebp + nb410_pos]
	mov   ecx, [esp + nb410_innerjjnr]
	mov   eax, [ecx]
	
	xorpd  xmm6, xmm6
	movapd xmm7, xmm6
	movsd  xmm7, [edx + eax*8]
	movlpd xmm6, [esi + eax*8]	;# xmm6(0) has the charge
	mulsd  xmm7, [esp + nb410_isai]
	movapd [esp + nb410_isaprod], xmm7
	movapd xmm1, xmm7
	mulpd xmm1, [esp + nb410_gbtsc]
	movapd [esp + nb410_gbscale], xmm1
	
	mulsd  xmm7, [esp + nb410_iq]
	mulsd  xmm6, xmm7
	movapd [esp + nb410_qq], xmm6
		
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov esi, [ebp + nb410_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb410_vdwparam]
	shl eax, 1
	mov edi, [esp + nb410_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [esp + nb410_c6], xmm4
	movapd [esp + nb410_c12], xmm6
	
	mov esi, [ebp + nb410_pos]       ;# base of pos[]
	
	movd  mm2, eax
	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb410_ix]
	movapd xmm5, [esp + nb410_iy]
	movapd xmm6, [esp + nb410_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb410_dx], xmm4
	movapd [esp + nb410_dy], xmm5
	movapd [esp + nb410_dz], xmm6
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
	movapd xmm1, [esp + nb410_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb410_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb410_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb410_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	movapd [esp + nb410_r], xmm4
	mulsd xmm4, [esp + nb410_gbscale]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb410_GBtab]

	movapd xmm4, [esi + eax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [esi + eax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb410_two]	;# two*Heps2 
	movapd xmm3, [esp + nb410_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
	;# get jnr from regs
	movd ebx, mm2
	mov esi, [ebp + nb410_dvda]
	
	;# Calculate dVda
	xorpd xmm7, xmm7
	mulsd xmm3, [esp + nb410_gbscale]
	movsd xmm6, xmm3
	mulsd  xmm6, [esp + nb410_r]
	addsd  xmm6, xmm5
	addsd  xmm5, [esp + nb410_vctot]
	movsd [esp + nb410_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subpd xmm7, xmm7
	movsd xmm6, xmm7
	
	;# update dvdasum
	addsd  xmm7, [esp + nb410_dvdasum]
	movsd [esp + nb410_dvdasum], xmm7 

	;# update j atoms dvdaj
	addsd  xmm6, [esi + ebx*8]
	movsd  [esi + ebx*8], xmm6
	
	;# L-J 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm0	;# xmm4=rinvsq 


	movapd xmm6, xmm4
	mulsd  xmm6, xmm4

	mulsd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulsd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulsd  xmm6, [esp + nb410_c6]
	mulsd  xmm4, [esp + nb410_c12]
	movapd xmm7, [esp + nb410_Vvdwtot]
	addsd  xmm7, xmm4
	mulsd  xmm4, [esp + nb410_twelve]
	subsd  xmm7, xmm6
	mulsd  xmm6, [esp + nb410_six]
	movlpd [esp + nb410_Vvdwtot], xmm7
	subsd  xmm4, xmm6
	mulsd  xmm4, xmm0
	subsd  xmm4, xmm3
	mulsd  xmm4, xmm0

	movapd xmm0, [esp + nb410_dx]
	movapd xmm1, [esp + nb410_dy]
	movapd xmm2, [esp + nb410_dz]

	movd eax, mm0	

	mov    edi, [ebp + nb410_faction]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb410_fix]
	movapd xmm4, [esp + nb410_fiy]
	movapd xmm5, [esp + nb410_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb410_fix], xmm3
	movlpd [esp + nb410_fiy], xmm4
	movlpd [esp + nb410_fiz], xmm5
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
		
.nb410_updateouterdata:
	mov   ecx, [esp + nb410_ii3]
	mov   edi, [ebp + nb410_faction]
	mov   esi, [ebp + nb410_fshift]
	mov   edx, [esp + nb410_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb410_fix]
	movapd xmm1, [esp + nb410_fiy]
	movapd xmm2, [esp + nb410_fiz]

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
	mov esi, [esp + nb410_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb410_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb410_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb410_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb410_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb410_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate dVda and update it 
	movapd xmm7, [esp + nb410_dvdasum]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	mov edx, [esp + nb410_ii]
	mov eax, [ebp + nb410_dvda]
	addsd xmm7, [eax + edx*8]
	movsd [eax + edx*8], xmm7
	
        ;# finish if last 
        mov ecx, [esp + nb410_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb410_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb410_n], esi
        jmp .nb410_outer
.nb410_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb410_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb410_end
        ;# non-zero, do one more workunit
        jmp   .nb410_threadloop
.nb410_end:
	emms

	mov eax, [esp + nb410_nouter]
	mov ebx, [esp + nb410_ninner]
	mov ecx, [ebp + nb410_outeriter]
	mov edx, [ebp + nb410_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb410_salign]
	add esp, eax
	add esp, 496
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret







.globl nb_kernel410nf_ia32_sse2
.globl _nb_kernel410nf_ia32_sse2
nb_kernel410nf_ia32_sse2:	
_nb_kernel410nf_ia32_sse2:	
.equiv          nb410nf_p_nri,          8
.equiv          nb410nf_iinr,           12
.equiv          nb410nf_jindex,         16
.equiv          nb410nf_jjnr,           20
.equiv          nb410nf_shift,          24
.equiv          nb410nf_shiftvec,       28
.equiv          nb410nf_fshift,         32
.equiv          nb410nf_gid,            36
.equiv          nb410nf_pos,            40
.equiv          nb410nf_faction,        44
.equiv          nb410nf_charge,         48
.equiv          nb410nf_p_facel,        52
.equiv          nb410nf_argkrf,         56
.equiv          nb410nf_argcrf,         60
.equiv          nb410nf_Vc,             64
.equiv          nb410nf_type,           68
.equiv          nb410nf_p_ntype,        72
.equiv          nb410nf_vdwparam,       76
.equiv          nb410nf_Vvdw,           80
.equiv          nb410nf_p_tabscale,     84
.equiv          nb410nf_VFtab,          88
.equiv          nb410nf_invsqrta,       92
.equiv          nb410nf_dvda,           96
.equiv          nb410nf_p_gbtabscale,   100
.equiv          nb410nf_GBtab,          104
.equiv          nb410nf_p_nthreads,     108
.equiv          nb410nf_count,          112
.equiv          nb410nf_mtx,            116
.equiv          nb410nf_outeriter,      120
.equiv          nb410nf_inneriter,      124
.equiv          nb410nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb410nf_ix,             0
.equiv          nb410nf_iy,             16
.equiv          nb410nf_iz,             32
.equiv          nb410nf_iq,             48
.equiv          nb410nf_two,            64
.equiv          nb410nf_gbtsc,          80
.equiv          nb410nf_qq,             96
.equiv          nb410nf_c6,             112
.equiv          nb410nf_c12,            128
.equiv          nb410nf_vctot,          144
.equiv          nb410nf_Vvdwtot,        160
.equiv          nb410nf_half,           176
.equiv          nb410nf_three,          192
.equiv          nb410nf_r,              208
.equiv          nb410nf_isai,           224
.equiv          nb410nf_isaprod,        240
.equiv          nb410nf_gbscale,        256
.equiv          nb410nf_ii,             272
.equiv          nb410nf_is3,            276
.equiv          nb410nf_ii3,            280
.equiv          nb410nf_ntia,           284
.equiv          nb410nf_innerjjnr,      288
.equiv          nb410nf_innerk,         292
.equiv          nb410nf_n,              296
.equiv          nb410nf_nn1,            300
.equiv          nb410nf_nri,            304
.equiv          nb410nf_facel,          312   ;# uses 8 bytes
.equiv          nb410nf_ntype,          320
.equiv          nb410nf_nouter,         324
.equiv          nb410nf_ninner,         328
.equiv          nb410nf_salign,         332
	push ebp
	mov ebp,esp	
    	push eax
    	push ebx
    	push ecx
    	push edx
	push esi
	push edi
	sub esp, 336		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb410nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb410nf_p_nri]
	mov esi, [ebp + nb410nf_p_facel]
	mov edi, [ebp + nb410nf_p_ntype]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov edi, [edi]
	mov [esp + nb410nf_nri], ecx
	movsd [esp + nb410nf_facel], xmm7
	mov [esp + nb410nf_ntype], edi

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb410nf_nouter], eax
	mov [esp + nb410nf_ninner], eax


	mov eax, [ebp + nb410nf_p_gbtabscale]
	movsd xmm5, [eax]
	shufpd xmm5, xmm5, 0
	movapd [esp + nb410nf_gbtsc], xmm5
	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb410nf_half], eax
	mov [esp + nb410nf_half+4], ebx
	movsd xmm1, [esp + nb410nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb410nf_half], xmm1
	movapd [esp + nb410nf_two], xmm2
	movapd [esp + nb410nf_three], xmm3

.nb410nf_threadloop:
        mov   esi, [ebp + nb410nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb410nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb410nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb410nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb410nf_n], eax
        mov [esp + nb410nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
		mov esi, eax				;# copy n to esi
        jg  .nb410nf_outerstart
        jmp .nb410nf_end

.nb410nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb410nf_nouter]
	mov [esp + nb410nf_nouter], ebx

.nb410nf_outer:
	mov   eax, [ebp + nb410nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb410nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb410nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb410nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 
	mov   [esp + nb410nf_ii], ebx

	mov   edx, [ebp + nb410nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb410nf_facel]
	shufpd xmm3, xmm3, 0

	mov   edx, [ebp + nb410nf_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [edx + ebx*8]
	shufpd xmm4, xmm4, 0

   	mov   edx, [ebp + nb410nf_type] 
   	mov   edx, [edx + ebx*4]
   	imul  edx, [esp + nb410nf_ntype]
   	shl   edx, 1
    mov   [esp + nb410nf_ntia], edx
		
	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb410nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb410nf_iq], xmm3
	movapd [esp + nb410nf_isai], xmm4

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb410nf_ix], xmm0
	movapd [esp + nb410nf_iy], xmm1
	movapd [esp + nb410nf_iz], xmm2

	mov   [esp + nb410nf_ii3], ebx
	
	;# clear vctot and Vvdwtot
	xorpd xmm4, xmm4
	movapd [esp + nb410nf_vctot], xmm4
	movapd [esp + nb410nf_Vvdwtot], xmm4
	
	mov   eax, [ebp + nb410nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb410nf_pos]
	mov   edi, [ebp + nb410nf_faction]	
	mov   eax, [ebp + nb410nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb410nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb410nf_ninner]
	mov   [esp + nb410nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb410nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb410nf_unroll_loop
	jmp   .nb410nf_checksingle
.nb410nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb410nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [edx]	
	mov   ebx, [edx + 4]              
	add dword ptr [esp + nb410nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	;# load isaj
	mov esi, [ebp + nb410nf_invsqrta]
	movlpd xmm2, [esi + eax*8]
	movhpd xmm2, [esi + ebx*8]
	mulpd  xmm2, [esp + nb410nf_isai]
	movapd [esp + nb410nf_isaprod], xmm2	
	movapd xmm1, xmm2
	mulpd xmm1, [esp + nb410nf_gbtsc]
	movapd [esp + nb410nf_gbscale], xmm1
	
	mov esi, [ebp + nb410nf_charge]    ;# base of charge[] 
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	mulpd xmm2, [esp + nb410nf_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb410nf_qq], xmm3	
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov esi, [ebp + nb410nf_type]
	mov eax, [esi + eax*4]
	mov ebx, [esi + ebx*4]
	mov esi, [ebp + nb410nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [esp + nb410nf_ntia]
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
	movapd [esp + nb410nf_c6], xmm4
	movapd [esp + nb410nf_c12], xmm6
	
	mov esi, [ebp + nb410nf_pos]       ;# base of pos[] 

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
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb410nf_ix]
	movapd xmm5, [esp + nb410nf_iy]
	movapd xmm6, [esp + nb410nf_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square dr 
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
	movapd xmm1, [esp + nb410nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb410nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb410nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb410nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	movapd [esp + nb410nf_r], xmm4
	mulpd xmm4, [esp + nb410nf_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb410nf_GBtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

	movapd xmm4, [esi + eax*8]	;# Y1 F1 	
	movapd xmm3, [esi + ebx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [esi + eax*8 + 16]	;# G1 H1 	
	movapd xmm3, [esi + ebx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb410nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  

	addpd  xmm5, [esp + nb410nf_vctot]
	movapd [esp + nb410nf_vctot], xmm5 

	;# L-J 
	movapd xmm4, xmm0
	mulpd  xmm4, xmm0	;# xmm4=rinvsq 

	movapd xmm6, xmm4
	mulpd  xmm6, xmm4
	
	mulpd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulpd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulpd  xmm6, [esp + nb410nf_c6]
	mulpd  xmm4, [esp + nb410nf_c12]
	movapd xmm7, [esp + nb410nf_Vvdwtot]
	addpd  xmm7, xmm4
	subpd  xmm7, xmm6
	movapd [esp + nb410nf_Vvdwtot], xmm7

	;# should we do one more iteration? 
	sub dword ptr [esp + nb410nf_innerk],  2
	jl    .nb410nf_checksingle
	jmp   .nb410nf_unroll_loop
.nb410nf_checksingle:
	mov   edx, [esp + nb410nf_innerk]
	and   edx, 1
	jnz    .nb410nf_dosingle
	jmp    .nb410nf_updateouterdata
.nb410nf_dosingle:
	mov esi, [ebp + nb410nf_charge]
	mov edx, [ebp + nb410nf_invsqrta]
	mov edi, [ebp + nb410nf_pos]
	mov   ecx, [esp + nb410nf_innerjjnr]
	mov   eax, [ecx]
	
	xorpd  xmm6, xmm6
	movapd xmm7, xmm6
	movsd  xmm7, [edx + eax*8]
	movlpd xmm6, [esi + eax*8]	;# xmm6(0) has the charge
	mulsd  xmm7, [esp + nb410nf_isai]
	movapd [esp + nb410nf_isaprod], xmm7
	movapd xmm1, xmm7
	mulpd xmm1, [esp + nb410nf_gbtsc]
	movapd [esp + nb410nf_gbscale], xmm1
	
	mulsd  xmm7, [esp + nb410nf_iq]
	mulsd  xmm6, xmm7
	movapd [esp + nb410nf_qq], xmm6
		
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov esi, [ebp + nb410nf_type]
	mov eax, [esi + eax*4]
	mov esi, [ebp + nb410nf_vdwparam]
	shl eax, 1
	mov edi, [esp + nb410nf_ntia]
	add eax, edi

	movlpd xmm6, [esi + eax*8]	;# c6a
	movhpd xmm6, [esi + eax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movapd [esp + nb410nf_c6], xmm4
	movapd [esp + nb410nf_c12], xmm6
	
	mov esi, [ebp + nb410nf_pos]       ;# base of pos[]
	
	movd  mm2, eax
	lea   eax, [eax + eax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb410nf_ix]
	movapd xmm5, [esp + nb410nf_iy]
	movapd xmm6, [esp + nb410nf_iz]

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
	movapd xmm1, [esp + nb410nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb410nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb410nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb410nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	movapd [esp + nb410nf_r], xmm4
	mulsd xmm4, [esp + nb410nf_gbscale]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb410nf_GBtab]

	movapd xmm4, [esi + eax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [esi + eax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb410nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  

	addsd  xmm5, [esp + nb410nf_vctot]
	movsd [esp + nb410nf_vctot], xmm5 

	;# L-J 
	movapd xmm4, xmm0
	mulsd  xmm4, xmm0	;# xmm4=rinvsq 


	movapd xmm6, xmm4
	mulsd  xmm6, xmm4

	mulsd  xmm6, xmm4	;# xmm6=rinvsix 
	movapd xmm4, xmm6
	mulsd  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulsd  xmm6, [esp + nb410nf_c6]
	mulsd  xmm4, [esp + nb410nf_c12]
	movapd xmm7, [esp + nb410nf_Vvdwtot]
	addsd  xmm7, xmm4
	subsd  xmm7, xmm6
	movlpd [esp + nb410nf_Vvdwtot], xmm7

.nb410nf_updateouterdata:
	mov   ecx, [esp + nb410nf_ii3]
	mov   edx, [esp + nb410nf_is3]

	;# get n from stack
	mov esi, [esp + nb410nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb410nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb410nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb410nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [esp + nb410nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   eax, [ebp + nb410nf_Vvdw]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb410nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb410nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb410nf_n], esi
        jmp .nb410nf_outer
.nb410nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb410nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb410nf_end
        ;# non-zero, do one more workunit
        jmp   .nb410nf_threadloop
.nb410nf_end:
	emms

	mov eax, [esp + nb410nf_nouter]
	mov ebx, [esp + nb410nf_ninner]
	mov ecx, [ebp + nb410nf_outeriter]
	mov edx, [ebp + nb410nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb410nf_salign]
	add esp, eax
	add esp, 336
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret


