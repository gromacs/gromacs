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




.globl nb_kernel400_ia32_sse2
.globl _nb_kernel400_ia32_sse2
nb_kernel400_ia32_sse2:	
_nb_kernel400_ia32_sse2:	
.equiv          nb400_p_nri,            8
.equiv          nb400_iinr,             12
.equiv          nb400_jindex,           16
.equiv          nb400_jjnr,             20
.equiv          nb400_shift,            24
.equiv          nb400_shiftvec,         28
.equiv          nb400_fshift,           32
.equiv          nb400_gid,              36
.equiv          nb400_pos,              40
.equiv          nb400_faction,          44
.equiv          nb400_charge,           48
.equiv          nb400_p_facel,          52
.equiv          nb400_argkrf,           56
.equiv          nb400_argcrf,           60
.equiv          nb400_Vc,               64
.equiv          nb400_type,             68
.equiv          nb400_p_ntype,          72
.equiv          nb400_vdwparam,         76
.equiv          nb400_Vvdw,             80
.equiv          nb400_p_tabscale,       84
.equiv          nb400_VFtab,            88
.equiv          nb400_invsqrta,         92
.equiv          nb400_dvda,             96
.equiv          nb400_p_gbtabscale,     100
.equiv          nb400_GBtab,            104
.equiv          nb400_p_nthreads,       108
.equiv          nb400_count,            112
.equiv          nb400_mtx,              116
.equiv          nb400_outeriter,        120
.equiv          nb400_inneriter,        124
.equiv          nb400_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb400_ix,               0
.equiv          nb400_iy,               16
.equiv          nb400_iz,               32
.equiv          nb400_iq,               48
.equiv          nb400_dx,               64
.equiv          nb400_dy,               80
.equiv          nb400_dz,               96
.equiv          nb400_two,              112
.equiv          nb400_gbtsc,            128
.equiv          nb400_qq,               144
.equiv          nb400_r,                160
.equiv          nb400_vctot,            176
.equiv          nb400_fix,              192
.equiv          nb400_fiy,              208
.equiv          nb400_fiz,              224
.equiv          nb400_half,             240
.equiv          nb400_three,            256
.equiv          nb400_isai,             272
.equiv          nb400_isaprod,          288
.equiv          nb400_dvdasum,          304
.equiv          nb400_gbscale,          320
.equiv          nb400_is3,              336
.equiv          nb400_ii3,              340
.equiv          nb400_ii,               344
.equiv          nb400_innerjjnr,        348
.equiv          nb400_innerk,           352
.equiv          nb400_n,                356
.equiv          nb400_nn1,              360
.equiv          nb400_nri,              364
.equiv          nb400_facel,            368   ;# uses 8 bytes
.equiv          nb400_nouter,           376
.equiv          nb400_ninner,           380
.equiv          nb400_salign,           384
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 388		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb400_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb400_p_nri]
	mov esi, [ebp + nb400_p_facel]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov [esp + nb400_nri], ecx
	movsd [esp + nb400_facel], xmm7

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb400_nouter], eax
	mov [esp + nb400_ninner], eax


	mov eax, [ebp + nb400_p_gbtabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb400_gbtsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb400_half], eax
	mov [esp + nb400_half + 4], ebx
	movsd xmm1, [esp + nb400_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb400_half], xmm1
	movapd [esp + nb400_two], xmm2
	movapd [esp + nb400_three], xmm3

.nb400_threadloop:
        mov   esi, [ebp + nb400_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb400_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb400_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb400_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb400_n], eax
        mov [esp + nb400_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb400_outerstart
        jmp .nb400_end

.nb400_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb400_nouter]
	mov [esp + nb400_nouter], ebx

.nb400_outer:
	mov   eax, [ebp + nb400_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb400_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb400_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb400_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 
	mov   [esp + nb400_ii], ebx
	
	mov   edx, [ebp + nb400_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb400_facel]
	shufpd xmm3, xmm3, 0

	mov   edx, [ebp + nb400_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [edx + ebx*8]
	shufpd xmm4, xmm4, 0

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb400_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb400_iq], xmm3
	movapd [esp + nb400_isai], xmm4
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb400_ix], xmm0
	movapd [esp + nb400_iy], xmm1
	movapd [esp + nb400_iz], xmm2

	mov   [esp + nb400_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb400_vctot], xmm4
	movapd [esp + nb400_dvdasum], xmm4
	movapd [esp + nb400_fix], xmm4
	movapd [esp + nb400_fiy], xmm4
	movapd [esp + nb400_fiz], xmm4
	
	mov   eax, [ebp + nb400_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb400_pos]
	mov   edi, [ebp + nb400_faction]	
	mov   eax, [ebp + nb400_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb400_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb400_ninner]
	mov   [esp + nb400_ninner], ecx
	add   edx, 0
	mov   [esp + nb400_innerk], edx    ;# number of innerloop atoms 
	jge   .nb400_unroll_loop
	jmp   .nb400_checksingle
.nb400_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb400_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [edx]
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb400_innerjjnr], 8	;# advance pointer (unrolled 2) 

	;# load isaj
	mov esi, [ebp + nb400_invsqrta]
	movlpd xmm2, [esi + eax*8]
	movhpd xmm2, [esi + ebx*8]
	mulpd  xmm2, [esp + nb400_isai]
	movapd [esp + nb400_isaprod], xmm2	
	movapd xmm1, xmm2
	mulpd xmm1, [esp + nb400_gbtsc]
	movapd [esp + nb400_gbscale], xmm1
	
	mov esi, [ebp + nb400_charge]    ;# base of charge[] 
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	mulpd xmm2, [esp + nb400_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb400_qq], xmm3	
	
	mov esi, [ebp + nb400_pos]		;# base of pos[] 

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

	mov    edi, [ebp + nb400_faction]
	
	;# move nb400_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb400_ix]
	movapd xmm5, [esp + nb400_iy]
	movapd xmm6, [esp + nb400_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb400_dx], xmm4
	movapd [esp + nb400_dy], xmm5
	movapd [esp + nb400_dz], xmm6
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
	movapd xmm1, [esp + nb400_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb400_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb400_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb400_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	movapd [esp + nb400_r], xmm4
	mulpd xmm4, [esp + nb400_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb400_GBtab]
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
	mulpd  xmm7, [esp + nb400_two]	;# two*Heps2 
	movapd xmm3, [esp + nb400_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
	;# get jnr from regs
	movd ecx, mm2
	movd edx, mm3
	mov esi, [ebp + nb400_dvda]
	
	;# Calculate dVda
	xorpd xmm7, xmm7
	mulpd xmm3, [esp + nb400_gbscale]
	movapd xmm6, xmm3
	mulpd  xmm6, [esp + nb400_r]
	addpd  xmm6, xmm5
	addpd  xmm5, [esp + nb400_vctot]
	movapd [esp + nb400_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subpd  xmm7, xmm6
	movapd xmm6, xmm7
	
	;# update dvdasum
	addpd  xmm7, [esp + nb400_dvdasum]
	movapd [esp + nb400_dvdasum], xmm7 

	;# update j atoms dvdaj
	movhlps xmm7, xmm6
	addsd  xmm6, [esi + ecx*8]
	addsd  xmm7, [esi + edx*8]
	movsd  [esi + ecx*8], xmm6
	movsd  [esi + edx*8], xmm7
	
	xorpd  xmm4, xmm4

	mulpd xmm3, xmm0
	subpd  xmm4, xmm3

	movapd xmm0, [esp + nb400_dx]
	movapd xmm1, [esp + nb400_dy]
	movapd xmm2, [esp + nb400_dz]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb400_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb400_fix]
	movapd xmm4, [esp + nb400_fiy]
	movapd xmm5, [esp + nb400_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb400_fix], xmm3
	movapd [esp + nb400_fiy], xmm4
	movapd [esp + nb400_fiz], xmm5
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
	sub dword ptr [esp + nb400_innerk],  2
	jl    .nb400_checksingle
	jmp   .nb400_unroll_loop
.nb400_checksingle:
	mov   edx, [esp + nb400_innerk]
	and   edx, 1
	jnz    .nb400_dosingle
	jmp    .nb400_updateouterdata
.nb400_dosingle:
	mov esi, [ebp + nb400_charge]
	mov edx, [ebp + nb400_invsqrta]
	mov edi, [ebp + nb400_pos]
	mov   ecx, [esp + nb400_innerjjnr]
	mov   eax, [ecx]	
	xorpd  xmm6, xmm6
	movapd xmm7, xmm6
	movsd  xmm7, [edx + eax*8]
	movlpd xmm6, [esi + eax*8]	;# xmm6(0) has the charge
	mulsd  xmm7, [esp + nb400_isai]
	movapd [esp + nb400_isaprod], xmm7
	movapd xmm1, xmm7
	mulpd xmm1, [esp + nb400_gbtsc]
	movapd [esp + nb400_gbscale], xmm1
	
	mulsd  xmm7, [esp + nb400_iq]
	mulsd  xmm6, xmm7
	movapd [esp + nb400_qq], xmm6

	movd  mm2, eax
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [edi + eax*8]
	movlpd xmm1, [edi + eax*8 + 8]
	movlpd xmm2, [edi + eax*8 + 16]

	;# move nb400_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb400_ix]
	movapd xmm5, [esp + nb400_iy]
	movapd xmm6, [esp + nb400_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb400_dx], xmm4
	movapd [esp + nb400_dy], xmm5
	movapd [esp + nb400_dz], xmm6
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
	movapd xmm1, [esp + nb400_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb400_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb400_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb400_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	movapd [esp + nb400_r], xmm4
	mulsd xmm4, [esp + nb400_gbscale]
	
	movd mm0, eax	

	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb400_GBtab]

	;# Coulomb 
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
	;# table ready in xmm4-xmm7 

	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb400_two]	;# two*Heps2 
	movapd xmm3, [esp + nb400_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq
	;# get jnr from regs
	movd ebx, mm2
	mov esi, [ebp + nb400_dvda]
	
	;# Calculate dVda
	mulsd xmm3, [esp + nb400_gbscale]
	movsd xmm6, xmm3
	mulsd  xmm6, [esp + nb400_r]
	addsd  xmm6, xmm5
	addsd  xmm5, [esp + nb400_vctot]
	movsd [esp + nb400_vctot], xmm5 

	;# xmm6=(vcoul+fijC*r)
	subpd  xmm7, xmm6
	movsd xmm6, xmm7
	
	;# update dvdasum
	addsd  xmm7, [esp + nb400_dvdasum]
	movsd [esp + nb400_dvdasum], xmm7 

	;# update j atoms dvdaj
	addsd  xmm6, [esi + ebx*8]
	movsd  [esi + ebx*8], xmm6
	
	xorpd xmm4, xmm4
	movd eax, mm0

	mulsd xmm3, xmm0
	subsd  xmm4, xmm3
	mov    edi, [ebp + nb400_faction]

	movsd xmm0, [esp + nb400_dx]
	movsd xmm1, [esp + nb400_dy]
	movsd xmm2, [esp + nb400_dz]

	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movsd xmm3, [esp + nb400_fix]
	movsd xmm4, [esp + nb400_fiy]
	movsd xmm5, [esp + nb400_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb400_fix], xmm3
	movlpd [esp + nb400_fiy], xmm4
	movlpd [esp + nb400_fiz], xmm5
	;# update fj 
	movlpd xmm3, [edi + eax*8]
	movlpd xmm4, [edi + eax*8 + 8]
	movlpd xmm5, [edi + eax*8 + 16]
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2
	movlpd [edi + eax*8], xmm3
	movlpd [edi + eax*8 + 8], xmm4
	movlpd [edi + eax*8 + 16], xmm5

.nb400_updateouterdata:
	mov   ecx, [esp + nb400_ii3]
	mov   edi, [ebp + nb400_faction]
	mov   esi, [ebp + nb400_fshift]
	mov   edx, [esp + nb400_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb400_fix]
	movapd xmm1, [esp + nb400_fiy]
	movapd xmm2, [esp + nb400_fiz]

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
	mov esi, [esp + nb400_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb400_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb400_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb400_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
	;# accumulate dVda and update it 
	movapd xmm7, [esp + nb400_dvdasum]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	mov edx, [esp + nb400_ii]
	mov eax, [ebp + nb400_dvda]
	addsd xmm7, [eax + edx*8]
	movsd [eax + edx*8], xmm7
	
        ;# finish if last 
        mov ecx, [esp + nb400_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb400_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb400_n], esi
        jmp .nb400_outer
.nb400_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb400_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb400_end
        ;# non-zero, do one more workunit
        jmp   .nb400_threadloop
.nb400_end:
	emms

	mov eax, [esp + nb400_nouter]
	mov ebx, [esp + nb400_ninner]
	mov ecx, [ebp + nb400_outeriter]
	mov edx, [ebp + nb400_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb400_salign]
	add esp, eax
	add esp, 388
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret






.globl nb_kernel400nf_ia32_sse2
.globl _nb_kernel400nf_ia32_sse2
nb_kernel400nf_ia32_sse2:	
_nb_kernel400nf_ia32_sse2:	
.equiv          nb400nf_p_nri,          8
.equiv          nb400nf_iinr,           12
.equiv          nb400nf_jindex,         16
.equiv          nb400nf_jjnr,           20
.equiv          nb400nf_shift,          24
.equiv          nb400nf_shiftvec,       28
.equiv          nb400nf_fshift,         32
.equiv          nb400nf_gid,            36
.equiv          nb400nf_pos,            40
.equiv          nb400nf_faction,        44
.equiv          nb400nf_charge,         48
.equiv          nb400nf_p_facel,        52
.equiv          nb400nf_argkrf,         56
.equiv          nb400nf_argcrf,         60
.equiv          nb400nf_Vc,             64
.equiv          nb400nf_type,           68
.equiv          nb400nf_p_ntype,        72
.equiv          nb400nf_vdwparam,       76
.equiv          nb400nf_Vvdw,           80
.equiv          nb400nf_p_tabscale,     84
.equiv          nb400nf_VFtab,          88
.equiv          nb400nf_invsqrta,       92
.equiv          nb400nf_dvda,           96
.equiv          nb400nf_p_gbtabscale,   100
.equiv          nb400nf_GBtab,          104
.equiv          nb400nf_p_nthreads,     108
.equiv          nb400nf_count,          112
.equiv          nb400nf_mtx,            116
.equiv          nb400nf_outeriter,      120
.equiv          nb400nf_inneriter,      124
.equiv          nb400nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb400nf_ix,             0
.equiv          nb400nf_iy,             16
.equiv          nb400nf_iz,             32
.equiv          nb400nf_iq,             48
.equiv          nb400nf_gbtsc,          64
.equiv          nb400nf_qq,             80
.equiv          nb400nf_vctot,          96
.equiv          nb400nf_half,           112
.equiv          nb400nf_three,          128
.equiv          nb400nf_isai,           144
.equiv          nb400nf_isaprod,        160
.equiv          nb400nf_gbscale,        176
.equiv          nb400nf_is3,            192
.equiv          nb400nf_ii3,            196
.equiv          nb400nf_innerjjnr,      200
.equiv          nb400nf_innerk,         204
.equiv          nb400nf_n,              208
.equiv          nb400nf_nn1,            212
.equiv          nb400nf_nri,            216
.equiv          nb400nf_facel,          224   ;# uses 8 bytes
.equiv          nb400nf_nouter,         232
.equiv          nb400nf_ninner,         236
.equiv          nb400nf_salign,         240
	push ebp
	mov ebp,esp	
	push eax
	push ebx
	push ecx
	push edx
	push esi
	push edi
	sub esp, 244		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb400nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb400nf_p_nri]
	mov esi, [ebp + nb400nf_p_facel]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov [esp + nb400nf_nri], ecx
	movsd [esp + nb400nf_facel], xmm7

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb400nf_nouter], eax
	mov [esp + nb400nf_ninner], eax


	mov eax, [ebp + nb400nf_p_gbtabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb400nf_gbtsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb400nf_half], eax
	mov [esp + nb400nf_half + 4], ebx
	movsd xmm1, [esp + nb400nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb400nf_half], xmm1
	movapd [esp + nb400nf_three], xmm3

.nb400nf_threadloop:
        mov   esi, [ebp + nb400nf_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb400nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock cmpxchg [esi], ebx                 ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb400nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb400nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb400nf_n], eax
        mov [esp + nb400nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb400nf_outerstart
        jmp .nb400nf_end

.nb400nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb400nf_nouter]
	mov [esp + nb400nf_nouter], ebx

.nb400nf_outer:
	mov   eax, [ebp + nb400nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb400nf_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb400nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb400nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb400nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb400nf_facel]
	shufpd xmm3, xmm3, 0

	mov   edx, [ebp + nb400nf_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [edx + ebx*8]
	shufpd xmm4, xmm4, 0

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb400nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb400nf_iq], xmm3
	movapd [esp + nb400nf_isai], xmm4
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb400nf_ix], xmm0
	movapd [esp + nb400nf_iy], xmm1
	movapd [esp + nb400nf_iz], xmm2

	mov   [esp + nb400nf_ii3], ebx
	
	;# clear vctot
	xorpd xmm4, xmm4
	movapd [esp + nb400nf_vctot], xmm4
	
	mov   eax, [ebp + nb400nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb400nf_pos]
	mov   edi, [ebp + nb400nf_faction]	
	mov   eax, [ebp + nb400nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb400nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb400nf_ninner]
	mov   [esp + nb400nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb400nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb400nf_unroll_loop
	jmp   .nb400nf_checksingle
.nb400nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb400nf_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [edx]
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb400nf_innerjjnr], 8	;# advance pointer (unrolled 2) 

	;# load isa2
	mov esi, [ebp + nb400nf_invsqrta]
	movlpd xmm2, [esi + eax*8]
	movhpd xmm2, [esi + ebx*8]
	mulpd  xmm2, [esp + nb400nf_isai]
	movapd [esp + nb400nf_isaprod], xmm2	
	movapd xmm1, xmm2
	mulpd xmm1, [esp + nb400nf_gbtsc]
	movapd [esp + nb400nf_gbscale], xmm1
	
	mov esi, [ebp + nb400nf_charge]    ;# base of charge[] 
	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	mulpd xmm2, [esp + nb400nf_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb400nf_qq], xmm3	
	
	mov esi, [ebp + nb400nf_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	mov    edi, [ebp + nb400nf_faction]
	
	;# move nb400nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb400nf_ix]
	movapd xmm5, [esp + nb400nf_iy]
	movapd xmm6, [esp + nb400nf_iz]

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
	movapd xmm1, [esp + nb400nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb400nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb400nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb400nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb400nf_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb400nf_GBtab]
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
	movapd xmm3, [esp + nb400nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	addpd  xmm5, [esp + nb400nf_vctot]
	movapd [esp + nb400nf_vctot], xmm5  
			
	;# should we do one more iteration? 
	sub dword ptr [esp + nb400nf_innerk],  2
	jl    .nb400nf_checksingle
	jmp   .nb400nf_unroll_loop
.nb400nf_checksingle:
	mov   edx, [esp + nb400nf_innerk]
	and   edx, 1
	jnz    .nb400nf_dosingle
	jmp    .nb400nf_updateouterdata
.nb400nf_dosingle:
	mov esi, [ebp + nb400nf_charge]
	mov edx, [ebp + nb400nf_invsqrta]
	mov edi, [ebp + nb400nf_pos]
	mov   ecx, [esp + nb400nf_innerjjnr]
	mov   eax, [ecx]	
	xorpd  xmm6, xmm6
	movapd xmm7, xmm6
	movsd  xmm7, [edx + eax*8]
	movlpd xmm6, [esi + eax*8]	;# xmm6(0) has the charge
	mulsd  xmm7, [esp + nb400nf_isai]
	movapd [esp + nb400nf_isaprod], xmm7
	movapd xmm1, xmm7
	mulpd xmm1, [esp + nb400nf_gbtsc]
	movapd [esp + nb400nf_gbscale], xmm1
	
	mulsd  xmm7, [esp + nb400nf_iq]
	mulsd  xmm6, xmm7
	movapd [esp + nb400nf_qq], xmm6
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [edi + eax*8]
	movlpd xmm1, [edi + eax*8 + 8]
	movlpd xmm2, [edi + eax*8 + 16]

	;# move nb400nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb400nf_ix]
	movapd xmm5, [esp + nb400nf_iy]
	movapd xmm6, [esp + nb400nf_iz]

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
	movapd xmm1, [esp + nb400nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb400nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb400nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb400nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb400nf_gbscale]
	
	movd mm0, eax	

	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb400nf_GBtab]

	;# Coulomb 
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
	;# table ready in xmm4-xmm7 

	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb400nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	addsd  xmm5, [esp + nb400nf_vctot]
	movsd [esp + nb400nf_vctot], xmm5
	
.nb400nf_updateouterdata:
	;# get n from stack
	mov esi, [esp + nb400nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb400nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb400nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb400nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb400nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jecxz .nb400nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb400nf_n], esi
        jmp .nb400nf_outer
.nb400nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb400nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jecxz .nb400nf_end
        ;# non-zero, do one more workunit
        jmp   .nb400nf_threadloop
.nb400nf_end:
	emms

	mov eax, [esp + nb400nf_nouter]
	mov ebx, [esp + nb400nf_ninner]
	mov ecx, [ebp + nb400nf_outeriter]
	mov edx, [ebp + nb400nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb400nf_salign]
	add esp, eax
	add esp, 244
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret



