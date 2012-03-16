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



.globl nb_kernel300_ia32_sse2
.globl _nb_kernel300_ia32_sse2
nb_kernel300_ia32_sse2:	
_nb_kernel300_ia32_sse2:	
.equiv          nb300_p_nri,            8
.equiv          nb300_iinr,             12
.equiv          nb300_jindex,           16
.equiv          nb300_jjnr,             20
.equiv          nb300_shift,            24
.equiv          nb300_shiftvec,         28
.equiv          nb300_fshift,           32
.equiv          nb300_gid,              36
.equiv          nb300_pos,              40
.equiv          nb300_faction,          44
.equiv          nb300_charge,           48
.equiv          nb300_p_facel,          52
.equiv          nb300_argkrf,           56
.equiv          nb300_argcrf,           60
.equiv          nb300_Vc,               64
.equiv          nb300_type,             68
.equiv          nb300_p_ntype,          72
.equiv          nb300_vdwparam,         76
.equiv          nb300_Vvdw,             80
.equiv          nb300_p_tabscale,       84
.equiv          nb300_VFtab,            88
.equiv          nb300_invsqrta,         92
.equiv          nb300_dvda,             96
.equiv          nb300_p_gbtabscale,     100
.equiv          nb300_GBtab,            104
.equiv          nb300_p_nthreads,       108
.equiv          nb300_count,            112
.equiv          nb300_mtx,              116
.equiv          nb300_outeriter,        120
.equiv          nb300_inneriter,        124
.equiv          nb300_work,             128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb300_ix,               0
.equiv          nb300_iy,               16
.equiv          nb300_iz,               32
.equiv          nb300_iq,               48
.equiv          nb300_dx,               64
.equiv          nb300_dy,               80
.equiv          nb300_dz,               96
.equiv          nb300_two,              112
.equiv          nb300_tsc,              128
.equiv          nb300_qq,               144
.equiv          nb300_fs,               160
.equiv          nb300_vctot,            176
.equiv          nb300_fix,              192
.equiv          nb300_fiy,              208
.equiv          nb300_fiz,              224
.equiv          nb300_half,             240
.equiv          nb300_three,            256
.equiv          nb300_is3,              272
.equiv          nb300_ii3,              276
.equiv          nb300_innerjjnr,        280
.equiv          nb300_innerk,           284
.equiv          nb300_n,                288
.equiv          nb300_nn1,              292
.equiv          nb300_nri,              296
.equiv          nb300_facel,            304   ;# uses 8 bytes
.equiv          nb300_nouter,           312
.equiv          nb300_ninner,           316
.equiv          nb300_salign,           320
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 324		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb300_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb300_p_nri]
	mov esi, [ebp + nb300_p_facel]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov [esp + nb300_nri], ecx
	movsd [esp + nb300_facel], xmm7

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb300_nouter], eax
	mov [esp + nb300_ninner], eax


	mov eax, [ebp + nb300_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb300_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb300_half], eax
	mov [esp + nb300_half+4], ebx
	movsd xmm1, [esp + nb300_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb300_half], xmm1
	movapd [esp + nb300_two], xmm2
	movapd [esp + nb300_three], xmm3

.nb300_threadloop:
        mov   esi, [ebp + nb300_count]          ;# pointer to sync counter
        mov   eax, [esi]
.nb300_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb300_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb300_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb300_n], eax
        mov [esp + nb300_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb300_outerstart
        jmp .nb300_end

.nb300_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb300_nouter]
	mov [esp + nb300_nouter], ebx

.nb300_outer:
	mov   eax, [ebp + nb300_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 
	mov   [esp + nb300_is3],ebx    	;# store is3 

	mov   eax, [ebp + nb300_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb300_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb300_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb300_facel]
	shufpd xmm3, xmm3, 0

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb300_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb300_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb300_ix], xmm0
	movapd [esp + nb300_iy], xmm1
	movapd [esp + nb300_iz], xmm2

	mov   [esp + nb300_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [esp + nb300_vctot], xmm4
	movapd [esp + nb300_fix], xmm4
	movapd [esp + nb300_fiy], xmm4
	movapd [esp + nb300_fiz], xmm4
	
	mov   eax, [ebp + nb300_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb300_pos]
	mov   edi, [ebp + nb300_faction]	
	mov   eax, [ebp + nb300_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb300_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb300_ninner]
	mov   [esp + nb300_ninner], ecx
	add   edx, 0
	mov   [esp + nb300_innerk], edx    ;# number of innerloop atoms 
	jge   .nb300_unroll_loop
	jmp   .nb300_checksingle
.nb300_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb300_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [edx]
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb300_innerjjnr], 8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb300_charge]    ;# base of charge[] 

	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm2, [esp + nb300_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb300_qq], xmm3	
	
	mov esi, [ebp + nb300_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	mov    edi, [ebp + nb300_faction]
	
	;# move nb300_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb300_ix]
	movapd xmm5, [esp + nb300_iy]
	movapd xmm6, [esp + nb300_iz]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# store dr 
	movapd [esp + nb300_dx], xmm4
	movapd [esp + nb300_dy], xmm5
	movapd [esp + nb300_dz], xmm6
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
	movapd xmm1, [esp + nb300_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb300_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb300_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb300_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb300_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  esi, [ebp + nb300_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm7, [esp + nb300_two]	;# two*Heps2 
	movapd xmm3, [esp + nb300_qq]
	addpd  xmm7, xmm6
	addpd  xmm7, xmm5 ;# xmm7=FF 
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulpd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [esp + nb300_vctot]
	movapd [esp + nb300_vctot], xmm5 

	xorpd  xmm4, xmm4

	mulpd xmm3, [esp + nb300_tsc]
	mulpd xmm3, xmm0
	subpd  xmm4, xmm3

	movapd xmm0, [esp + nb300_dx]
	movapd xmm1, [esp + nb300_dy]
	movapd xmm2, [esp + nb300_dz]

	movd eax, mm0	
	movd ebx, mm1

	mov    edi, [ebp + nb300_faction]
	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb300_fix]
	movapd xmm4, [esp + nb300_fiy]
	movapd xmm5, [esp + nb300_fiz]
	addpd  xmm3, xmm0
	addpd  xmm4, xmm1
	addpd  xmm5, xmm2
	movapd [esp + nb300_fix], xmm3
	movapd [esp + nb300_fiy], xmm4
	movapd [esp + nb300_fiz], xmm5
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
	sub dword ptr [esp + nb300_innerk],  2
	jl    .nb300_checksingle
	jmp   .nb300_unroll_loop
.nb300_checksingle:
	mov   edx, [esp + nb300_innerk]
	and   edx, 1
	jnz    .nb300_dosingle
	jmp    .nb300_updateouterdata
.nb300_dosingle:
	mov esi, [ebp + nb300_charge]
	mov edi, [ebp + nb300_pos]
	mov   ecx, [esp + nb300_innerjjnr]
	mov   eax, [ecx]	
	xorpd  xmm6, xmm6
	movlpd xmm6, [esi + eax*8]	;# xmm6(0) has the charge 	
	mulsd  xmm6, [esp + nb300_iq]
	movapd [esp + nb300_qq], xmm6
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [edi + eax*8]
	movlpd xmm1, [edi + eax*8 + 8]
	movlpd xmm2, [edi + eax*8 + 16]

	;# move nb300_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb300_ix]
	movapd xmm5, [esp + nb300_iy]
	movapd xmm6, [esp + nb300_iz]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# store dr 
	movapd [esp + nb300_dx], xmm4
	movapd [esp + nb300_dy], xmm5
	movapd [esp + nb300_dz], xmm6
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
	movapd xmm1, [esp + nb300_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb300_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb300_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb300_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb300_tsc]
	
	movd mm0, eax	

	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb300_VFtab]

	;# Coulomb 
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
	;# table ready in xmm4-xmm7 

	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [esp + nb300_two]	;# two*Heps2 
	movapd xmm3, [esp + nb300_qq]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [esp + nb300_vctot]
	movsd [esp + nb300_vctot], xmm5 

	xorpd xmm4, xmm4
	movd eax, mm0

	mulpd xmm3, [esp + nb300_tsc]
	mulpd xmm3, xmm0
	subpd  xmm4, xmm3
	mov    edi, [ebp + nb300_faction]

	movapd xmm0, [esp + nb300_dx]
	movapd xmm1, [esp + nb300_dy]
	movapd xmm2, [esp + nb300_dz]

	mulpd  xmm0, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# now update f_i 
	movapd xmm3, [esp + nb300_fix]
	movapd xmm4, [esp + nb300_fiy]
	movapd xmm5, [esp + nb300_fiz]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm5, xmm2
	movlpd [esp + nb300_fix], xmm3
	movlpd [esp + nb300_fiy], xmm4
	movlpd [esp + nb300_fiz], xmm5
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

.nb300_updateouterdata:
	mov   ecx, [esp + nb300_ii3]
	mov   edi, [ebp + nb300_faction]
	mov   esi, [ebp + nb300_fshift]
	mov   edx, [esp + nb300_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [esp + nb300_fix]
	movapd xmm1, [esp + nb300_fiy]
	movapd xmm2, [esp + nb300_fiz]

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
	mov esi, [esp + nb300_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb300_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb300_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb300_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb300_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb300_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb300_n], esi
        jmp .nb300_outer
.nb300_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb300_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb300_end
        ;# non-zero, do one more workunit
        jmp   .nb300_threadloop
.nb300_end:
	emms

	mov eax, [esp + nb300_nouter]
	mov ebx, [esp + nb300_ninner]
	mov ecx, [ebp + nb300_outeriter]
	mov edx, [ebp + nb300_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb300_salign]
	add esp, eax
	add esp, 324
	pop edi
	pop esi
    	pop edx
    	pop ecx
    	pop ebx
    	pop eax
	leave
	ret




	

.globl nb_kernel300nf_ia32_sse2
.globl _nb_kernel300nf_ia32_sse2
nb_kernel300nf_ia32_sse2:	
_nb_kernel300nf_ia32_sse2:	
.equiv          nb300nf_p_nri,          8
.equiv          nb300nf_iinr,           12
.equiv          nb300nf_jindex,         16
.equiv          nb300nf_jjnr,           20
.equiv          nb300nf_shift,          24
.equiv          nb300nf_shiftvec,       28
.equiv          nb300nf_fshift,         32
.equiv          nb300nf_gid,            36
.equiv          nb300nf_pos,            40
.equiv          nb300nf_faction,        44
.equiv          nb300nf_charge,         48
.equiv          nb300nf_p_facel,        52
.equiv          nb300nf_argkrf,         56
.equiv          nb300nf_argcrf,         60
.equiv          nb300nf_Vc,             64
.equiv          nb300nf_type,           68
.equiv          nb300nf_p_ntype,        72
.equiv          nb300nf_vdwparam,       76
.equiv          nb300nf_Vvdw,           80
.equiv          nb300nf_p_tabscale,     84
.equiv          nb300nf_VFtab,          88
.equiv          nb300nf_invsqrta,       92
.equiv          nb300nf_dvda,           96
.equiv          nb300nf_p_gbtabscale,   100
.equiv          nb300nf_GBtab,          104
.equiv          nb300nf_p_nthreads,     108
.equiv          nb300nf_count,          112
.equiv          nb300nf_mtx,            116
.equiv          nb300nf_outeriter,      120
.equiv          nb300nf_inneriter,      124
.equiv          nb300nf_work,           128
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb300nf_ix,             0
.equiv          nb300nf_iy,             16
.equiv          nb300nf_iz,             32
.equiv          nb300nf_iq,             48
.equiv          nb300nf_tsc,            64
.equiv          nb300nf_qq,             80
.equiv          nb300nf_vctot,          96
.equiv          nb300nf_half,           112
.equiv          nb300nf_three,          128
.equiv          nb300nf_is3,            144
.equiv          nb300nf_ii3,            148
.equiv          nb300nf_innerjjnr,      152
.equiv          nb300nf_innerk,         156
.equiv          nb300nf_n,              160
.equiv          nb300nf_nn1,            164
.equiv          nb300nf_nri,            168
.equiv          nb300nf_facel,          176   ;# uses 8 bytes
.equiv          nb300nf_nouter,         184
.equiv          nb300nf_ninner,         188
.equiv          nb300nf_salign,         192
	push ebp
	mov ebp,esp	
    push eax
    push ebx
    push ecx
    push edx
	push esi
	push edi
	sub esp, 196		;# local stack space 
	mov  eax, esp
	and  eax, 0xf
	sub esp, eax
	mov [esp + nb300nf_salign], eax

	emms

	;# Move args passed by reference to stack
	mov ecx, [ebp + nb300nf_p_nri]
	mov esi, [ebp + nb300nf_p_facel]
	mov ecx, [ecx]
	movsd xmm7, [esi]
	mov [esp + nb300nf_nri], ecx
	movsd [esp + nb300nf_facel], xmm7

	;# zero iteration counters
	mov eax, 0
	mov [esp + nb300nf_nouter], eax
	mov [esp + nb300nf_ninner], eax


	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double 0.5 IEEE (hex)
	mov ebx, 0x3fe00000
	mov [esp + nb300nf_half], eax
	mov [esp + nb300nf_half+4], ebx
	movsd xmm1, [esp + nb300nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# 1.0
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# 2.0
	addpd  xmm3, xmm2	;# 3.0
	movapd [esp + nb300nf_half], xmm1
	movapd [esp + nb300nf_three], xmm3
	mov eax, [ebp + nb300nf_p_tabscale]
	movsd xmm3, [eax]
	shufpd xmm3, xmm3, 0
	movapd [esp + nb300nf_tsc], xmm3

.nb300nf_threadloop:
        mov   esi, [ebp + nb300nf_count]        ;# pointer to sync counter
        mov   eax, [esi]
.nb300nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [esi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb300nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [esp + nb300nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [esp + nb300nf_n], eax
        mov [esp + nb300nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb300nf_outerstart
        jmp .nb300nf_end

.nb300nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [esp + nb300nf_nouter]
	mov [esp + nb300nf_nouter], ebx

.nb300nf_outer:
	mov   eax, [ebp + nb300nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [eax+esi*4]		;# ebx=shift[n] 
	
	lea   ebx, [ebx + ebx*2]    ;# ebx=3*is 

	mov   eax, [ebp + nb300nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [eax + ebx*8]
	movsd xmm1, [eax + ebx*8 + 8]
	movsd xmm2, [eax + ebx*8 + 16] 

	mov   ecx, [ebp + nb300nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [ecx+esi*4]	    ;# ebx =ii 

	mov   edx, [ebp + nb300nf_charge]
	movsd xmm3, [edx + ebx*8]	
	mulsd xmm3, [esp + nb300nf_facel]
	shufpd xmm3, xmm3, 0

	lea   ebx, [ebx + ebx*2]	;# ebx = 3*ii=ii3 
	mov   eax, [ebp + nb300nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [eax + ebx*8]
	addsd xmm1, [eax + ebx*8 + 8]
	addsd xmm2, [eax + ebx*8 + 16]

	movapd [esp + nb300nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [esp + nb300nf_ix], xmm0
	movapd [esp + nb300nf_iy], xmm1
	movapd [esp + nb300nf_iz], xmm2

	mov   [esp + nb300nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [esp + nb300nf_vctot], xmm4
	
	mov   eax, [ebp + nb300nf_jindex]
	mov   ecx, [eax + esi*4]	     ;# jindex[n] 
	mov   edx, [eax + esi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   esi, [ebp + nb300nf_pos]
	mov   eax, [ebp + nb300nf_jjnr]
	shl   ecx, 2
	add   eax, ecx
	mov   [esp + nb300nf_innerjjnr], eax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [esp + nb300nf_ninner]
	mov   [esp + nb300nf_ninner], ecx
	add   edx, 0
	mov   [esp + nb300nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb300nf_unroll_loop
	jmp   .nb300nf_checksingle
.nb300nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   edx, [esp + nb300nf_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [edx]
	mov   ebx, [edx + 4]
	add dword ptr [esp + nb300nf_innerjjnr], 8	;# advance pointer (unrolled 2) 

	mov esi, [ebp + nb300nf_charge]    ;# base of charge[] 

	movlpd xmm3, [esi + eax*8]
	movhpd xmm3, [esi + ebx*8]

	movapd xmm2, [esp + nb300nf_iq]
	mulpd  xmm3, xmm2
	movapd [esp + nb300nf_qq], xmm3	
	
	mov esi, [ebp + nb300nf_pos]		;# base of pos[] 

	lea   eax, [eax + eax*2]     ;# replace jnr with j3 
	lea   ebx, [ebx + ebx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [esi + eax*8]
	movlpd xmm1, [esi + eax*8 + 8]
	movlpd xmm2, [esi + eax*8 + 16]
	movhpd xmm0, [esi + ebx*8]
	movhpd xmm1, [esi + ebx*8 + 8]
	movhpd xmm2, [esi + ebx*8 + 16]		

	;# move nb300nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb300nf_ix]
	movapd xmm5, [esp + nb300nf_iy]
	movapd xmm6, [esp + nb300nf_iz]

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
	movapd xmm1, [esp + nb300nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb300nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb300nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb300nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [esp + nb300nf_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb300nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb300nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul  
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [esp + nb300nf_vctot]
	movapd [esp + nb300nf_vctot], xmm5 

	;# should we do one more iteration? 
	sub dword ptr [esp + nb300nf_innerk],  2
	jl    .nb300nf_checksingle
	jmp   .nb300nf_unroll_loop
.nb300nf_checksingle:
	mov   edx, [esp + nb300nf_innerk]
	and   edx, 1
	jnz    .nb300nf_dosingle
	jmp    .nb300nf_updateouterdata
.nb300nf_dosingle:
	mov esi, [ebp + nb300nf_charge]
	mov edi, [ebp + nb300nf_pos]
	mov   ecx, [esp + nb300nf_innerjjnr]
	mov   eax, [ecx]	
	xorpd  xmm6, xmm6
	movlpd xmm6, [esi + eax*8]	;# xmm6(0) has the charge 	
	mulsd  xmm6, [esp + nb300nf_iq]
	movapd [esp + nb300nf_qq], xmm6
		
	lea   eax, [eax + eax*2]
	
	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [edi + eax*8]
	movlpd xmm1, [edi + eax*8 + 8]
	movlpd xmm2, [edi + eax*8 + 16]

	;# move nb300nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [esp + nb300nf_ix]
	movapd xmm5, [esp + nb300nf_iy]
	movapd xmm6, [esp + nb300nf_iz]

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
	movapd xmm1, [esp + nb300nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb300nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [esp + nb300nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [esp + nb300nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [esp + nb300nf_tsc]
	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	
	mov  esi, [ebp + nb300nf_VFtab]

	;# Coulomb 
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
	;# table ready in xmm4-xmm7 

	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [esp + nb300nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [esp + nb300nf_vctot]
	movsd [esp + nb300nf_vctot], xmm5 

.nb300nf_updateouterdata:
	;# get group index for i particle 
	;# get n from stack
	mov esi, [esp + nb300nf_n]
        ;# get group index for i particle 
        mov   edx, [ebp + nb300nf_gid]      	;# base of gid[]
        mov   edx, [edx + esi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [esp + nb300nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   eax, [ebp + nb300nf_Vc]
	addsd xmm7, [eax + edx*8] 
	;# move back to mem 
	movsd [eax + edx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [esp + nb300nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb300nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [esp + nb300nf_n], esi
        jmp .nb300nf_outer
.nb300nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [esp + nb300nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb300nf_end
        ;# non-zero, do one more workunit
        jmp   .nb300nf_threadloop
.nb300nf_end:
	emms

	mov eax, [esp + nb300nf_nouter]
	mov ebx, [esp + nb300nf_ninner]
	mov ecx, [ebp + nb300nf_outeriter]
	mov edx, [ebp + nb300nf_inneriter]
	mov [ecx], eax
	mov [edx], ebx

	mov eax, [esp + nb300nf_salign]
	add esp, eax
	add esp, 196
	pop edi
	pop esi
    pop edx
    pop ecx
    pop ebx
    pop eax
	leave
	ret

