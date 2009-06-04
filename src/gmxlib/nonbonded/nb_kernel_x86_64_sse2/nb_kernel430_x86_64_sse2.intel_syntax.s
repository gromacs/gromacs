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


.globl nb_kernel430_x86_64_sse2
.globl _nb_kernel430_x86_64_sse2
nb_kernel430_x86_64_sse2:	
_nb_kernel430_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb430_fshift,           16
.equiv          nb430_gid,              24
.equiv          nb430_pos,              32
.equiv          nb430_faction,          40
.equiv          nb430_charge,           48
.equiv          nb430_p_facel,          56
.equiv          nb430_argkrf,           64
.equiv          nb430_argcrf,           72
.equiv          nb430_Vc,               80
.equiv          nb430_type,             88
.equiv          nb430_p_ntype,          96
.equiv          nb430_vdwparam,         104
.equiv          nb430_Vvdw,             112
.equiv          nb430_p_tabscale,       120
.equiv          nb430_VFtab,            128
.equiv          nb430_invsqrta,         136
.equiv          nb430_dvda,             144
.equiv          nb430_p_gbtabscale,     152
.equiv          nb430_GBtab,            160
.equiv          nb430_p_nthreads,       168
.equiv          nb430_count,            176
.equiv          nb430_mtx,              184
.equiv          nb430_outeriter,        192
.equiv          nb430_inneriter,        200
.equiv          nb430_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb430_ix,               0
.equiv          nb430_iy,               16
.equiv          nb430_iz,               32
.equiv          nb430_iq,               48
.equiv          nb430_dx,               64
.equiv          nb430_dy,               80
.equiv          nb430_dz,               96
.equiv          nb430_eps,              112
.equiv          nb430_gbtsc,            128
.equiv          nb430_tsc,              144
.equiv          nb430_qq,               160
.equiv          nb430_c6,               176
.equiv          nb430_c12,              192
.equiv          nb430_epsgb,            208
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
.equiv          nb430_rinv,             416
.equiv          nb430_nri,              432
.equiv          nb430_iinr,             440
.equiv          nb430_jindex,           448
.equiv          nb430_jjnr,             456
.equiv          nb430_shift,            464
.equiv          nb430_shiftvec,         472
.equiv          nb430_facel,            480
.equiv          nb430_innerjjnr,        488
.equiv          nb430_ii,               496
.equiv          nb430_is3,              500
.equiv          nb430_ii3,              504
.equiv          nb430_ntia,             508
.equiv          nb430_innerk,           512
.equiv          nb430_n,                516
.equiv          nb430_nn1,              520
.equiv          nb430_ntype,            524
.equiv          nb430_nouter,           528
.equiv          nb430_ninner,           532

	push rbp
	mov  rbp, rsp
	push rbx

	
	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 536		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb430_nouter], eax
	mov [rsp + nb430_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb430_nri], edi
	mov [rsp + nb430_iinr], rsi
	mov [rsp + nb430_jindex], rdx
	mov [rsp + nb430_jjnr], rcx
	mov [rsp + nb430_shift], r8
	mov [rsp + nb430_shiftvec], r9
	mov rdi, [rbp + nb430_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb430_ntype], edi
	mov rsi, [rbp + nb430_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb430_facel], xmm0

	mov rax, [rbp + nb430_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb430_tsc], xmm3

	mov rbx, [rbp + nb430_p_gbtabscale]
	movsd xmm4, [rbx]
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb430_gbtsc], xmm4

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb430_half], eax
	mov [rsp + nb430_half+4], ebx
	movsd xmm1, [rsp + nb430_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb430_half], xmm1
	movapd [rsp + nb430_three], xmm3

.nb430_threadloop:
        mov   rsi, [rbp + nb430_count]          ;# pointer to sync counter
        mov   eax, [rsi]
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
        mov ecx, [rsp + nb430_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb430_n], eax
        mov [rsp + nb430_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb430_outerstart
        jmp .nb430_end

.nb430_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb430_nouter]
	mov [rsp + nb430_nouter], ebx

.nb430_outer:
	mov   rax, [rsp + nb430_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb430_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb430_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb430_iinr]       ;# rcx = pointer into iinr[]
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 
	mov   [rsp + nb430_ii], ebx

	mov   rdx, [rbp + nb430_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb430_facel]
	shufpd xmm3, xmm3, 0

	mov   rdx, [rbp + nb430_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [rdx + rbx*8]
	shufpd xmm4, xmm4, 0

    	mov   rdx, [rbp + nb430_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb430_ntype]
    	shl   edx, 1
    	mov   [rsp + nb430_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb430_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb430_iq], xmm3
	movapd [rsp + nb430_isai], xmm4
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb430_ix], xmm0
	movapd [rsp + nb430_iy], xmm1
	movapd [rsp + nb430_iz], xmm2

	mov   [rsp + nb430_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb430_vctot], xmm4
	movapd [rsp + nb430_Vvdwtot], xmm4
	movapd [rsp + nb430_dvdasum], xmm4
	movapd [rsp + nb430_fix], xmm4
	movapd [rsp + nb430_fiy], xmm4
	movapd [rsp + nb430_fiz], xmm4
	
	mov   rax, [rsp + nb430_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb430_pos]
	mov   rdi, [rbp + nb430_faction]	
	mov   rax, [rsp + nb430_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb430_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb430_ninner]
	mov   [rsp + nb430_ninner], ecx
	add   edx, 0
	mov   [rsp + nb430_innerk], edx    ;# number of innerloop atoms 
	jge   .nb430_unroll_loop
	jmp   .nb430_checksingle
.nb430_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb430_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	mov   ebx, [rdx + 4]
	add qword ptr [rsp + nb430_innerjjnr], 8	;# advance pointer (unrolled 2) 

		
	mov rsi, [rbp + nb430_pos]		;# base of pos[] 

	lea   r10, [rax + rax*2]     ;# j3 
	lea   r11, [rbx + rbx*2]	

	;# move two coordinates to xmm4-xmm6 
	movlpd xmm4, [rsi + r10*8]
	movlpd xmm5, [rsi + r10*8 + 8]
	movlpd xmm6, [rsi + r10*8 + 16]
	movhpd xmm4, [rsi + r11*8]
	movhpd xmm5, [rsi + r11*8 + 8]
	movhpd xmm6, [rsi + r11*8 + 16]		
	
	;# calc dr 
	subpd xmm4, [rsp + nb430_ix]
	subpd xmm5, [rsp + nb430_iy]
	subpd xmm6, [rsp + nb430_iz]

	;# store dr 
	movapd [rsp + nb430_dx], xmm4
	movapd [rsp + nb430_dy], xmm5
	movapd [rsp + nb430_dz], xmm6
    
	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	;# load isaj
	mov rsi, [rbp + nb430_invsqrta]
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	mulpd  xmm3, [rsp + nb430_isai]
	movapd [rsp + nb430_isaprod], xmm3
	movapd xmm6, xmm3
	mulpd xmm3, [rsp + nb430_gbtsc]
	movapd [rsp + nb430_gbscale], xmm3
	
    ;#invsqrt
	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb430_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb430_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

    mulpd  xmm6, [rsp + nb430_iq]
	mov rsi, [rbp + nb430_charge]    ;# base of charge[] 
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	mulpd  xmm3, xmm6
	movapd [rsp + nb430_qq], xmm3	

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb430_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb430_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv 
	mulpd xmm4, xmm0	;# xmm4=r 
	movapd [rsp + nb430_r], xmm4
	movapd [rsp + nb430_rinv], xmm0

	mov rsi, [rbp + nb430_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	shl r8d, 1
	shl r9d, 1
	mov edi, [rsp + nb430_ntia]
	add r8d, edi
	add r9d, edi

    movapd xmm8, xmm4 ;# r
	mulpd xmm4, [rsp + nb430_gbscale]
	mulpd xmm8, [rsp + nb430_tsc]
    
    ;# truncate and convert to integers
    cvttpd2pi mm0, xmm4  ;# gb
    cvttpd2pi mm1, xmm8  ;# lj
    
    ;# convert back to float
    cvtpi2pd  xmm6, mm0   ;# gb
    cvtpi2pd  xmm10, mm1  ;# lj
    
    ;# multiply by 4 and 8, respectively
    pslld   mm0, 2   ;# gb
    pslld   mm1, 3   ;# lj

    ;# move to integer registers
    movd    r12d, mm0       ;# gb
    movd    r14d, mm1      ;# lj
	psrlq mm0, 32
	psrlq mm1, 32
    movd    r13d, mm0      ;# gb
    movd    r15d, mm1     ;# lj
    ;# GB indices: r10-11   LJ indices: r12-r13

    ;# calculate eps
    subpd     xmm4, xmm6   ;# gb
    subpd     xmm8, xmm10  ;# lj
    movapd    [rsp + nb430_epsgb], xmm4 ;# gb eps
    movapd    [rsp + nb430_eps], xmm8 ;# lj eps
    
	mov  rsi, [rbp + nb430_GBtab]
	mov  rdi, [rbp + nb430_VFtab]

    ;# load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
    movapd xmm0,  [rsi + r12*8]        ;# Y1c F1c
    movapd xmm12, [rsi + r13*8]        ;# Y2c F2c
    movapd xmm4,  [rdi + r14*8]        ;# Y1d F1d
    movapd xmm13, [rdi + r15*8]        ;# Y2d F2d
    movapd xmm8,  [rdi + r14*8 + 32]   ;# Y1r F1r
    movapd xmm14, [rdi + r15*8 + 32]   ;# Y2r F2r
	movapd xmm1, xmm0
	movapd xmm5, xmm4
	movapd xmm9, xmm8
	unpcklpd xmm0, xmm12	;# Y1c Y2c 
	unpckhpd xmm1, xmm12	;# F1c F2c 
	unpcklpd xmm4, xmm13	;# Y1d Y2d 
	unpckhpd xmm5, xmm13	;# F1d F2d 
	unpcklpd xmm8, xmm14	;# Y1r Y2r 
	unpckhpd xmm9, xmm14	;# F1r F2r 
    
    movapd xmm2,  [rsi + r12*8 + 16]   ;# G1c H1c
    movapd xmm12, [rsi + r13*8 + 16]   ;# G2c H2c
    movapd xmm6,  [rdi + r14*8 + 16]   ;# G1d H1d
    movapd xmm13, [rdi + r15*8 + 16]   ;# G2d H2d
    movapd xmm10, [rdi + r14*8 + 48]   ;# G1r H1r
    movapd xmm14, [rdi + r15*8 + 48]   ;# G2r H2r
	movapd xmm3, xmm2
	movapd xmm7, xmm6
	movapd xmm11, xmm10
	unpcklpd xmm2, xmm12	;# G1c G2c 
	unpckhpd xmm3, xmm12	;# H1c H2c 
	unpcklpd xmm6, xmm13	;# G1d G2d 
	unpckhpd xmm7, xmm13	;# H1d H2d 
	unpcklpd xmm10, xmm14	;# G1r G2r 
	unpckhpd xmm11, xmm14	;# H1r H2r 
    ;# table data ready. Coul GB in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11
    mov rdi, [rbp + nb430_vdwparam]
    
    movapd xmm12, [rsp + nb430_epsgb]
    movapd xmm13, [rsp + nb430_eps]
    
    mulpd  xmm3, xmm12   ;# Heps
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm13
    mulpd  xmm2, xmm12     ;# Geps
    mulpd  xmm6, xmm13
    mulpd  xmm10, xmm13
    mulpd  xmm3, xmm12   ;# Heps2
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm13

    movlpd xmm14, [rdi + r8*8]
    movlpd xmm15, [rdi + r8*8 + 8]
    
    addpd  xmm1, xmm2   ;# F+Geps
    addpd  xmm5, xmm6
    addpd  xmm9, xmm10 
    addpd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addpd  xmm5, xmm7
    addpd  xmm9, xmm11 
    addpd  xmm3, xmm3    ;# 2*Heps2
    addpd  xmm7, xmm7
    addpd  xmm11, xmm11
    movhpd xmm14, [rdi + r9*8]
    movhpd xmm15, [rdi + r9*8 + 8]
    
    addpd  xmm3, xmm2    ;# 2*Heps2+Geps
    addpd  xmm7, xmm6  
    addpd  xmm11, xmm10
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, xmm12   ;# eps*Fp
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm13
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, [rsp + nb430_qq]   ;# VV*qq = vcoul
    mulpd  xmm5, xmm14   ;# vnb6
    mulpd  xmm9, xmm15   ;# vnb12
    mulpd  xmm3, [rsp + nb430_qq]    ;# FF*qq = fij
    mulpd  xmm7, xmm14   ;# fijD
    mulpd  xmm11, xmm15   ;#fijR

    addpd  xmm11, xmm7 ;# fijD+fijR
    mulpd  xmm11, [rsp + nb430_tsc] ;# (fijD+fijR)*tabscale
    
    ;# accumulate Vvdwtot
    addpd  xmm5, [rsp + nb430_Vvdwtot]
    addpd  xmm5, xmm9
    movapd [rsp + nb430_Vvdwtot], xmm5

	mov rsi, [rbp + nb430_dvda]
	
	;# Calculate dVda
	mulpd xmm3, [rsp + nb430_gbscale]   ;# fijC=qq*FF*gbscale
	movapd xmm6, xmm3 
	mulpd  xmm6, [rsp + nb430_r]
	addpd  xmm6, xmm1   ;# vcoul+fijC*r

    addpd  xmm3, xmm11  ;# fijC+fijD+fijR
    
    ;# increment vctot
	addpd  xmm1, [rsp + nb430_vctot]
    movapd [rsp + nb430_vctot], xmm1

	;# xmm6=(vcoul+fijC*r)
	xorpd  xmm7, xmm7
	subpd  xmm7, xmm6
	movapd xmm6, xmm7
	
    ;# the fj's - start by combiningg forces from memory 
    mov rdi, [rbp + nb430_faction]
	movlpd xmm0, [rdi + r10*8]
	movlpd xmm1, [rdi + r10*8 + 8]
	movlpd xmm2, [rdi + r10*8 + 16]
	movhpd xmm0, [rdi + r11*8]
	movhpd xmm1, [rdi + r11*8 + 8]
	movhpd xmm2, [rdi + r11*8 + 16]

	;# update dvdasum 
	addpd  xmm7, [rsp + nb430_dvdasum]
    movapd [rsp + nb430_dvdasum], xmm7

	;# update j atoms dvdaj
	movhlps xmm7, xmm6
	addsd  xmm6, [rsi + rax*8]
	addsd  xmm7, [rsi + rbx*8]
	movsd  [rsi + rax*8], xmm6
	movsd  [rsi + rbx*8], xmm7

	xorpd  xmm4, xmm4	
	mulpd xmm3, [rsp + nb430_rinv]
	subpd  xmm4, xmm3

    movapd  xmm9, xmm4
    movapd  xmm10, xmm4
    movapd  xmm11, xmm4
    
    mulpd  xmm9, [rsp + nb430_dx]
    mulpd  xmm10, [rsp + nb430_dy]
    mulpd  xmm11, [rsp + nb430_dz]    

	addpd xmm0, xmm9
	addpd xmm1, xmm10
	addpd xmm2, xmm11

	;# accumulate i forces
    addpd xmm9, [rsp + nb430_fix]
    addpd xmm10, [rsp + nb430_fiy]
    addpd xmm11, [rsp + nb430_fiz]

	movlpd [rdi + r10*8], xmm0
	movlpd [rdi + r10*8 + 8], xmm1
	movlpd [rdi + r10*8 + 16], xmm2

    movapd [rsp + nb430_fix], xmm9
    movapd [rsp + nb430_fiy], xmm10
    movapd [rsp + nb430_fiz], xmm11

	movhpd [rdi + r11*8], xmm0
	movhpd [rdi + r11*8 + 8], xmm1
	movhpd [rdi + r11*8 + 16], xmm2
	
    ;# should we do one more iteration? 
	sub dword ptr [rsp + nb430_innerk],  2
	jl    .nb430_checksingle
	jmp   .nb430_unroll_loop
.nb430_checksingle:
	mov   edx, [rsp + nb430_innerk]
	and   edx, 1
	jnz    .nb430_dosingle
	jmp    .nb430_updateouterdata
.nb430_dosingle:
	mov rsi, [rbp + nb430_charge]
	mov rdx, [rbp + nb430_invsqrta]
	mov rdi, [rbp + nb430_pos]
	mov   rcx, [rsp + nb430_innerjjnr]
	mov   eax, [rcx]	

	;# load isaj
	mov rsi, [rbp + nb430_invsqrta]
	movsd xmm2, [rsi + rax*8]
	mulsd  xmm2, [rsp + nb430_isai]
	movapd [rsp + nb430_isaprod], xmm2	
	movapd xmm1, xmm2
	mulsd xmm1, [rsp + nb430_gbtsc]
	movapd [rsp + nb430_gbscale], xmm1

    mulsd xmm2, [rsp + nb430_iq]
	mov rsi, [rbp + nb430_charge]    ;# base of charge[] 
	movsd xmm3, [rsi + rax*8]
	mulsd  xmm3, xmm2
	movapd [rsp + nb430_qq], xmm3	
	
	mov rsi, [rbp + nb430_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb430_vdwparam]
	shl r8d, 1
	mov edi, [rsp + nb430_ntia]
	add r8d, edi

	movsd xmm4, [rsi + r8*8]	
	movsd xmm6, [rsi + r8*8 + 8]
	movapd [rsp + nb430_c6], xmm4
	movapd [rsp + nb430_c12], xmm6
		
	mov rsi, [rbp + nb430_pos]		;# base of pos[] 

	lea   r10, [rax + rax*2]     ;# j3 

	;# move coordinate to xmm4-xmm6 
	movsd xmm4, [rsi + r10*8]
	movsd xmm5, [rsi + r10*8 + 8]
	movsd xmm6, [rsi + r10*8 + 16]

	mov    rdi, [rbp + nb430_faction]
	
	;# calc dr 
	subsd xmm4, [rsp + nb430_ix]
	subsd xmm5, [rsp + nb430_iy]
	subsd xmm6, [rsp + nb430_iz]

	;# store dr 
	movapd [rsp + nb430_dx], xmm4
	movapd [rsp + nb430_dy], xmm5
	movapd [rsp + nb430_dz], xmm6
    
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
	movapd xmm1, [rsp + nb430_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb430_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb430_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb430_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv 
	mulsd xmm4, xmm0	;# xmm4=r 
	movapd [rsp + nb430_r], xmm4
	movapd [rsp + nb430_rinv], xmm0

    movapd xmm8, xmm4 ;# r
	mulsd xmm4, [rsp + nb430_gbscale]
	mulsd xmm8, [rsp + nb430_tsc]
    
    ;# truncate and convert to integers
    cvttsd2si r12d, xmm4  ;# gb
    cvttsd2si r14d, xmm8  ;# lj
    
    ;# convert back to float
    cvtsi2sd  xmm6, r12d   ;# gb
    cvtsi2sd  xmm10, r14d  ;# lj
    
    ;# multiply by 4 and 8, respectively
    shl    r12d, 2   ;# gb
    shl    r14d, 3   ;# lj

    ;# GB indices: r10   LJ indices: r12

    ;# calculate eps
    subsd     xmm4, xmm6   ;# gb
    subsd     xmm8, xmm10  ;# lj
    movapd    [rsp + nb430_epsgb], xmm4 ;# gb eps
    movapd    [rsp + nb430_eps], xmm8 ;# lj eps
    
	mov  rsi, [rbp + nb430_GBtab]
	mov  rdi, [rbp + nb430_VFtab]

    ;# load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
    movapd xmm0,  [rsi + r12*8]        ;# Y1c F1c
    movapd xmm4,  [rdi + r14*8]        ;# Y1d F1d
    movapd xmm8,  [rdi + r14*8 + 32]   ;# Y1r F1r
	movhlps xmm1, xmm0
	movhlps xmm5, xmm4
	movhlps xmm9, xmm8
    
    movapd xmm2,  [rsi + r12*8 + 16]   ;# G1c H1c
    movapd xmm6,  [rdi + r14*8 + 16]   ;# G1d H1d
    movapd xmm10, [rdi + r14*8 + 48]   ;# G1r H1r
	movhlps xmm3, xmm2
	movhlps xmm7, xmm6
	movhlps xmm11, xmm10
    ;# table data ready. Coul GB in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11
    
    movapd xmm12, [rsp + nb430_epsgb]
    movapd xmm13, [rsp + nb430_eps]
    
    mulsd  xmm3, xmm12   ;# Heps
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm13
    mulsd  xmm2, xmm12     ;# Geps
    mulsd  xmm6, xmm13
    mulsd  xmm10, xmm13
    mulsd  xmm3, xmm12   ;# Heps2
    mulsd  xmm7, xmm13
    mulsd  xmm11, xmm13

    addsd  xmm1, xmm2   ;# F+Geps
    addsd  xmm5, xmm6
    addsd  xmm9, xmm10 
    addsd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addsd  xmm5, xmm7
    addsd  xmm9, xmm11 
    addsd  xmm3, xmm3    ;# 2*Heps2
    addsd  xmm7, xmm7
    addsd  xmm11, xmm11
    addsd  xmm3, xmm2    ;# 2*Heps2+Geps
    addsd  xmm7, xmm6  
    addsd  xmm11, xmm10
    addsd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm7, xmm5
    addsd  xmm11, xmm9
    mulsd  xmm1, xmm12   ;# eps*Fp
    mulsd  xmm5, xmm13
    mulsd  xmm9, xmm13
    addsd  xmm1, xmm0     ;# VV
    addsd  xmm5, xmm4
    addsd  xmm9, xmm8
    mulsd  xmm1, [rsp + nb430_qq]   ;# VV*qq = vcoul
    mulsd  xmm5, [rsp + nb430_c6]   ;# vnb6
    mulsd  xmm9, [rsp + nb430_c12]   ;# vnb12
    mulsd  xmm3, [rsp + nb430_qq]    ;# FF*qq = fij
    mulsd  xmm7, [rsp + nb430_c6]   ;# fijD
    mulsd  xmm11, [rsp + nb430_c12]   ;#fijR

    addsd  xmm11, xmm7 ;# fijD+fijR
    mulsd  xmm11, [rsp + nb430_tsc] ;# (fijD+fijR)*tabscale
    
    ;# accumulate Vvdwtot
    addsd  xmm5, [rsp + nb430_Vvdwtot]
    addsd  xmm5, xmm9
    movsd [rsp + nb430_Vvdwtot], xmm5

	mov rsi, [rbp + nb430_dvda]
	
	;# Calculate dVda
	mulsd xmm3, [rsp + nb430_gbscale]   ;# fijC=qq*FF*gbscale
	movapd xmm6, xmm3 
	mulsd  xmm6, [rsp + nb430_r]
	addsd  xmm6, xmm1   ;# vcoul+fijC*r

    addsd  xmm3, xmm11  ;# fijC+fijD+fijR
    
    ;# increment vctot
	addsd  xmm1, [rsp + nb430_vctot]
    movsd [rsp + nb430_vctot], xmm1

	;# xmm6=(vcoul+fijC*r)
	xorpd  xmm7, xmm7
	subsd  xmm7, xmm6
	movapd xmm6, xmm7
	
	;# update dvdasum 
	addsd  xmm7, [rsp + nb430_dvdasum]
    movsd [rsp + nb430_dvdasum], xmm7

	;# update j atoms dvdaj
	addsd  xmm6, [rsi + rax*8]
	movsd  [rsi + rax*8], xmm6

	xorpd  xmm4, xmm4	
	mulsd xmm3, [rsp + nb430_rinv]
	subsd  xmm4, xmm3

    movapd  xmm9, xmm4
    movapd  xmm10, xmm4
    movapd  xmm11, xmm4
    
    mulsd  xmm9, [rsp + nb430_dx]
    mulsd  xmm10, [rsp + nb430_dy]
    mulsd  xmm11, [rsp + nb430_dz]
    
    movapd xmm3, xmm9
    movapd xmm4, xmm10
    movapd xmm5, xmm11
    
	;# accumulate i forces
    addsd xmm9, [rsp + nb430_fix]
    addsd xmm10, [rsp + nb430_fiy]
    addsd xmm11, [rsp + nb430_fiz]
    movsd [rsp + nb430_fix], xmm9
    movsd [rsp + nb430_fiy], xmm10
    movsd [rsp + nb430_fiz], xmm11

    mov rdi, [rbp + nb430_faction]
	;# the fj's - start by accumulating forces from memory 
	addsd xmm3,   [rdi + r10*8]
	addsd xmm4,  [rdi + r10*8 + 8]
	addsd xmm5,  [rdi + r10*8 + 16]
	movsd [rdi + r10*8], xmm3
	movsd [rdi + r10*8 + 8], xmm4
	movsd [rdi + r10*8 + 16], xmm5
	
.nb430_updateouterdata:
	mov   ecx, [rsp + nb430_ii3]
	mov   rdi, [rbp + nb430_faction]
	mov   rsi, [rbp + nb430_fshift]
	mov   edx, [rsp + nb430_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb430_fix]
	movapd xmm1, [rsp + nb430_fiy]
	movapd xmm2, [rsp + nb430_fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8]
	movsd  xmm4, [rdi + rcx*8 + 8]
	movsd  xmm5, [rdi + rcx*8 + 16]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8],     xmm3
	movsd  [rdi + rcx*8 + 8], xmm4
	movsd  [rdi + rcx*8 + 16], xmm5

	;# increment fshift force  
	movsd  xmm3, [rsi + rdx*8]
	movsd  xmm4, [rsi + rdx*8 + 8]
	movsd  xmm5, [rsi + rdx*8 + 16]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rsi + rdx*8],     xmm3
	movsd  [rsi + rdx*8 + 8], xmm4
	movsd  [rsi + rdx*8 + 16], xmm5

	;# get n from stack
	mov esi, [rsp + nb430_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb430_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb430_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb430_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb430_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb430_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate dVda and update it 
	movapd xmm7, [rsp + nb430_dvdasum]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	mov edx, [rsp + nb430_ii]
	mov rax, [rbp + nb430_dvda]
	addsd xmm7, [rax + rdx*8]
	movsd [rax + rdx*8], xmm7
	
        ;# finish if last 
        mov ecx, [rsp + nb430_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb430_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb430_n], esi
        jmp .nb430_outer
.nb430_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb430_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb430_end
        ;# non-zero, do one more workunit
        jmp   .nb430_threadloop
.nb430_end:
	mov eax, [rsp + nb430_nouter]
	mov ebx, [rsp + nb430_ninner]
	mov rcx, [rbp + nb430_outeriter]
	mov rdx, [rbp + nb430_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 536
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret





	
.globl nb_kernel430nf_x86_64_sse2
.globl _nb_kernel430nf_x86_64_sse2
nb_kernel430nf_x86_64_sse2:	
_nb_kernel430nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb430nf_fshift,         16
.equiv          nb430nf_gid,            24
.equiv          nb430nf_pos,            32
.equiv          nb430nf_faction,        40
.equiv          nb430nf_charge,         48
.equiv          nb430nf_p_facel,        56
.equiv          nb430nf_argkrf,         64
.equiv          nb430nf_argcrf,         72
.equiv          nb430nf_Vc,             80
.equiv          nb430nf_type,           88
.equiv          nb430nf_p_ntype,        96
.equiv          nb430nf_vdwparam,       104
.equiv          nb430nf_Vvdw,           112
.equiv          nb430nf_p_tabscale,     120
.equiv          nb430nf_VFtab,          128
.equiv          nb430nf_invsqrta,       136
.equiv          nb430nf_dvda,           144
.equiv          nb430nf_p_gbtabscale,   152
.equiv          nb430nf_GBtab,          160
.equiv          nb430nf_p_nthreads,     168
.equiv          nb430nf_count,          176
.equiv          nb430nf_mtx,            184
.equiv          nb430nf_outeriter,      192
.equiv          nb430nf_inneriter,      200
.equiv          nb430nf_work,           208
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
.equiv          nb430nf_nri,            272
.equiv          nb430nf_iinr,           280
.equiv          nb430nf_jindex,         288
.equiv          nb430nf_jjnr,           296
.equiv          nb430nf_shift,          304
.equiv          nb430nf_shiftvec,       312
.equiv          nb430nf_facel,          320
.equiv          nb430nf_innerjjnr,      328
.equiv          nb430nf_is3,            336
.equiv          nb430nf_ii3,            340
.equiv          nb430nf_ntia,           344
.equiv          nb430nf_innerk,         348
.equiv          nb430nf_n,              352
.equiv          nb430nf_nn1,            356
.equiv          nb430nf_ntype,          360
.equiv          nb430nf_nouter,         364
.equiv          nb430nf_ninner,         368
	push rbp
	mov  rbp, rsp
	push rbx
	
	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 392		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb430nf_nouter], eax
	mov [rsp + nb430nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb430nf_nri], edi
	mov [rsp + nb430nf_iinr], rsi
	mov [rsp + nb430nf_jindex], rdx
	mov [rsp + nb430nf_jjnr], rcx
	mov [rsp + nb430nf_shift], r8
	mov [rsp + nb430nf_shiftvec], r9
	mov rdi, [rbp + nb430nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb430nf_ntype], edi
	mov rsi, [rbp + nb430nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb430nf_facel], xmm0

	mov rax, [rbp + nb430nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb430nf_tsc], xmm3

	mov rbx, [rbp + nb430nf_p_gbtabscale]
	movsd xmm4, [rbx]
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb430nf_gbtsc], xmm4

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb430nf_half], eax
	mov [rsp + nb430nf_half+4], ebx
	movsd xmm1, [rsp + nb430nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb430nf_half], xmm1
	movapd [rsp + nb430nf_three], xmm3

.nb430nf_threadloop:
        mov   rsi, [rbp + nb430nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
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
        mov ecx, [rsp + nb430nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb430nf_n], eax
        mov [rsp + nb430nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb430nf_outerstart
        jmp .nb430nf_end

.nb430nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb430nf_nouter]
	mov [rsp + nb430nf_nouter], ebx

.nb430nf_outer:
	mov   rax, [rsp + nb430nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb430nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb430nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb430nf_iinr]       ;# rcx = pointer into iinr[]
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb430nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb430nf_facel]
	shufpd xmm3, xmm3, 0

	mov   rdx, [rbp + nb430nf_invsqrta]	;# load invsqrta[ii]
	movsd xmm4, [rdx + rbx*8]
	shufpd xmm4, xmm4, 0

    	mov   rdx, [rbp + nb430nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb430nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb430nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb430nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb430nf_iq], xmm3
	movapd [rsp + nb430nf_isai], xmm4	
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb430nf_ix], xmm0
	movapd [rsp + nb430nf_iy], xmm1
	movapd [rsp + nb430nf_iz], xmm2

	mov   [rsp + nb430nf_ii3], ebx
	
	;# clear vctot
	xorpd xmm4, xmm4
	movapd [rsp + nb430nf_vctot], xmm4
	movapd [rsp + nb430nf_Vvdwtot], xmm4

	mov   rax, [rsp + nb430nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb430nf_pos]
	mov   rdi, [rbp + nb430nf_faction]	
	mov   rax, [rsp + nb430nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb430nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb430nf_ninner]
	mov   [rsp + nb430nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb430nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb430nf_unroll_loop
	jmp   .nb430nf_checksingle
.nb430nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb430nf_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	mov   ebx, [rdx + 4]
	add qword ptr [rsp + nb430nf_innerjjnr], 8	;# advance pointer (unrolled 2) 

	;# load isaj
	mov rsi, [rbp + nb430nf_invsqrta]
	movlpd xmm2, [rsi + rax*8]
	movhpd xmm2, [rsi + rbx*8]
	mulpd  xmm2, [rsp + nb430nf_isai]
	movapd [rsp + nb430nf_isaprod], xmm2	
	movapd xmm1, xmm2
	mulpd xmm1, [rsp + nb430nf_gbtsc]
	movapd [rsp + nb430nf_gbscale], xmm1
	
	mov rsi, [rbp + nb430nf_charge]    ;# base of charge[] 
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	mulpd xmm2, [rsp + nb430nf_iq]
	mulpd  xmm3, xmm2
	movapd [rsp + nb430nf_qq], xmm3	
	
	mov rsi, [rbp + nb430nf_type]
	mov ecx, [rsi + rax*4]
	mov edx, [rsi + rbx*4]
	mov rsi, [rbp + nb430nf_vdwparam]
	shl ecx, 1
	shl edx, 1
	mov edi, [rsp + nb430nf_ntia]
	add ecx, edi
	add edx, edi

	movlpd xmm6, [rsi + rcx*8]	;# c6a
	movlpd xmm7, [rsi + rdx*8]	;# c6b
	movhpd xmm6, [rsi + rcx*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rdx*8 + 8]	;# c6b c12b 

	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb430nf_c6], xmm4
	movapd [rsp + nb430nf_c12], xmm6
	
	mov rsi, [rbp + nb430nf_pos]		;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	mov    rdi, [rbp + nb430nf_faction]
	
	;# move nb430nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb430nf_ix]
	movapd xmm5, [rsp + nb430nf_iy]
	movapd xmm6, [rsp + nb430nf_iz]

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
	movapd xmm1, [rsp + nb430nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb430nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb430nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb430nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv 
	mulpd xmm4, xmm0	;# xmm4=r 
	movapd [rsp + nb430nf_r], xmm4
	mulpd xmm4, [rsp + nb430nf_gbscale]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 

	mov  rsi, [rbp + nb430nf_GBtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6		;# indices in eax/ebx 

	;# Coulomb 
	movapd xmm4, [rsi + rcx*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rdx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rcx*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rdx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb430nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
	addpd  xmm5, [rsp + nb430nf_vctot]
	movapd [rsp + nb430nf_vctot], xmm5
	
	movapd xmm4, [rsp + nb430nf_r]
	mulpd  xmm4, [rsp + nb430nf_tsc]
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 3		;# idx *= 8

	mov  rsi, [rbp + nb430nf_VFtab]

	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6		;# indices in eax/ebx 

	;# Dispersion 
	movapd xmm4, [rsi + rcx*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rdx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rcx*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rdx*8 + 16]	;# G2 H2 
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

	mulpd  xmm5, [rsp + nb430nf_c6]	 ;# Vvdw6
	addpd  xmm5, [rsp + nb430nf_Vvdwtot]
	movapd [rsp + nb430nf_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [rsi + rcx*8 + 32]	;# Y1 F1 	
	movapd xmm3, [rsi + rdx*8 + 32]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rcx*8 + 48]	;# G1 H1 	
	movapd xmm3, [rsi + rdx*8 + 48]	;# G2 H2 
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

	mulpd  xmm5, [rsp + nb430nf_c12] ;# Vvdw12 
	addpd  xmm5, [rsp + nb430nf_Vvdwtot]
	movapd [rsp + nb430nf_Vvdwtot], xmm5
	xorpd  xmm4, xmm4
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb430nf_innerk],  2
	jl    .nb430nf_checksingle
	jmp   .nb430nf_unroll_loop
.nb430nf_checksingle:
	mov   edx, [rsp + nb430nf_innerk]
	and   edx, 1
	jnz    .nb430nf_dosingle
	jmp    .nb430nf_updateouterdata
.nb430nf_dosingle:
	mov rsi, [rbp + nb430nf_charge]
	mov rdx, [rbp + nb430nf_invsqrta]
	mov rdi, [rbp + nb430nf_pos]
	mov   rcx, [rsp + nb430nf_innerjjnr]
	mov   eax, [rcx]	

	xorpd  xmm6, xmm6
	movapd xmm7, xmm6
	movsd  xmm7, [rdx + rax*8]
	movlpd xmm6, [rsi + rax*8]	;# xmm6(0) has the charge
	mulsd  xmm7, [rsp + nb430nf_isai]
	movapd [rsp + nb430nf_isaprod], xmm7
	movapd xmm1, xmm7
	mulpd xmm1, [rsp + nb430nf_gbtsc]
	movapd [rsp + nb430nf_gbscale], xmm1
	
	mulsd  xmm7, [rsp + nb430nf_iq]
	mulsd  xmm6, xmm7
	movapd [rsp + nb430nf_qq], xmm6
	
	mov rsi, [rbp + nb430nf_type]
	mov edx, [rsi + rax*4]
	mov rsi, [rbp + nb430nf_vdwparam]
	shl edx, 1
	mov edi, [rsp + nb430nf_ntia]
	add edx, edi

	movlpd xmm6, [rsi + rdx*8]	;# c6a
	movhpd xmm6, [rsi + rdx*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb430nf_c6], xmm4
	movapd [rsp + nb430nf_c12], xmm6
	
	mov rsi, [rbp + nb430nf_pos]		;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	mov    rdi, [rbp + nb430nf_faction]

	;# move nb430nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb430nf_ix]
	movapd xmm5, [rsp + nb430nf_iy]
	movapd xmm6, [rsp + nb430nf_iz]

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
	movapd xmm1, [rsp + nb430nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb430nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb430nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb430nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
	movsd [rsp + nb430nf_r], xmm4
	mulsd xmm4, [rsp + nb430nf_gbscale]
	
	cvttsd2si edx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, edx
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl edx, 2		;# idx *= 4 
	mov  rsi, [rbp + nb430nf_GBtab]

	;# Coulomb 
	movapd xmm4, [rsi + rdx*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rdx*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb430nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	addsd  xmm5, [rsp + nb430nf_vctot]
	movsd [rsp + nb430nf_vctot], xmm5 

	movsd xmm4, [rsp + nb430nf_r]
	mulsd  xmm4, [rsp + nb430nf_tsc]
	cvttsd2si edx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, edx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2

	shl edx, 3

	mov  rsi, [rbp + nb430nf_VFtab]

	;# Dispersion 
	movapd xmm4, [rsi + rdx*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rdx*8 + 16]	;# G1 H1 	
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

	mulsd  xmm5, [rsp + nb430nf_c6]	 ;# Vvdw6
	addsd  xmm5, [rsp + nb430nf_Vvdwtot]
	movlpd [rsp + nb430nf_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [rsi + rdx*8 + 32]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rdx*8 + 48]	;# G1 H1 	
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
	mulsd  xmm5, [rsp + nb430nf_c12] ;# Vvdw12 
	addsd  xmm5, [rsp + nb430nf_Vvdwtot]
	movlpd [rsp + nb430nf_Vvdwtot], xmm5
.nb430nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb430nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb430nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb430nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb430nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb430nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb430nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb430nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb430nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb430nf_n], esi
        jmp .nb430nf_outer
.nb430nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb430nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb430nf_end
        ;# non-zero, do one more workunit
        jmp   .nb430nf_threadloop
.nb430nf_end:
	mov eax, [rsp + nb430nf_nouter]
	mov ebx, [rsp + nb430nf_ninner]
	mov rcx, [rbp + nb430nf_outeriter]
	mov rdx, [rbp + nb430nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 392
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret



