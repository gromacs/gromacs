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
%macro .equiv                  2
   %1 equ %2
%endmacro
; .endif                   # End of NASM-specific block
; .intel_syntax noprefix   # Line only read by gnu as





.globl nb_kernel430_x86_64_sse
.globl _nb_kernel430_x86_64_sse
nb_kernel430_x86_64_sse:	
_nb_kernel430_x86_64_sse:	
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
	;# bottom of stack is cache-aligned for sse use 
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

	sub rsp, 552		;# local variable stack space (n*16+8)

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
	movss xmm0, [rsi]
	movss [rsp + nb430_facel], xmm0

	mov rax, [rbp + nb430_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb430_tsc], xmm3

	mov rbx, [rbp + nb430_p_gbtabscale]
	movss xmm4, [rbx]
	shufps xmm4, xmm4, 0
	movaps [rsp + nb430_gbtsc], xmm4


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb430_half], eax
	movss xmm1, [rsp + nb430_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb430_half],  xmm1
	movaps [rsp + nb430_three],  xmm3

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
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb430_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb430_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb430_iinr]       ;# rcx = pointer into iinr[]
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 
	mov   [rsp + nb430_ii], ebx

	mov   rdx, [rbp + nb430_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb430_facel]
	shufps xmm3, xmm3, 0

	mov   rdx, [rbp + nb430_invsqrta]	;# load invsqrta[ii]
	movss xmm4, [rdx + rbx*4]
	shufps xmm4, xmm4, 0

    	mov   rdx, [rbp + nb430_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb430_ntype]
    	shl   edx, 1
    	mov   [rsp + nb430_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb430_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb430_iq], xmm3
	movaps [rsp + nb430_isai], xmm4
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb430_ix], xmm0
	movaps [rsp + nb430_iy], xmm1
	movaps [rsp + nb430_iz], xmm2

	mov   [rsp + nb430_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb430_vctot], xmm4
	movaps [rsp + nb430_Vvdwtot], xmm4
	movaps [rsp + nb430_dvdasum], xmm4
	movaps [rsp + nb430_fix], xmm4
	movaps [rsp + nb430_fiy], xmm4
	movaps [rsp + nb430_fiz], xmm4
	
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
	sub   edx,  4
	add   ecx, [rsp + nb430_ninner]
	mov   [rsp + nb430_ninner], ecx
	add   edx, 0
	mov   [rsp + nb430_innerk], edx    ;# number of innerloop atoms
	
	jge   .nb430_unroll_loop
	jmp   .nb430_finish_inner
.nb430_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb430_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb430_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	;# load isaj
	mov rsi, [rbp + nb430_invsqrta]
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]
	movaps xmm2, [rsp + nb430_isai]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all isaj in xmm3 
	mulps  xmm2, xmm3
	
	movaps [rsp + nb430_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb430_gbtsc]
	movaps [rsp + nb430_gbscale], xmm1
	
	mov rsi, [rbp + nb430_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	mulps xmm2, [rsp + nb430_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm3, xmm2
	movaps [rsp + nb430_qq], xmm3	
	
    ;# vdw parameters
	mov rsi, [rbp + nb430_type]
	mov r12d, [rsi + rax*4]
	mov r13d, [rsi + rbx*4]
	mov r14d, [rsi + rcx*4]
	mov r15d, [rsi + rdx*4]
	shl r12d, 1	
	shl r13d, 1	
	shl r14d, 1	
	shl r15d, 1	
    mov edi, [rsp + nb430_ntia]
	add r12d, edi
	add r13d, edi
	add r14d, edi
	add r15d, edi

	mov rsi, [rbp + nb430_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movlps xmm7, [rsi + r14*4]
	movhps xmm3, [rsi + r13*4]
	movhps xmm7, [rsi + r15*4]

	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101

    movaps [rsp + nb430_c6], xmm0
    movaps [rsp + nb430_c12], xmm3
    
	mov rsi, [rbp + nb430_pos]       ;# base of pos[] 
		
	lea   r8, [rax + rax*2]     ;# jnr
	lea   r9, [rbx + rbx*2]	
	lea   r10, [rcx + rcx*2]    
	lea   r11, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + r8*4]
	movlps xmm5, [rsi + r10*4]
	movss xmm2, [rsi + r8*4 + 8]
	movss xmm6, [rsi + r10*4 + 8]

	movhps xmm4, [rsi + r9*4]
	movhps xmm5, [rsi + r11*4]

	movss xmm0, [rsi + r9*4 + 8]
	movss xmm1, [rsi + r11*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

	;# calc dr 
	subps xmm0, [rsp + nb430_ix]
	subps xmm1, [rsp + nb430_iy]
	subps xmm2, [rsp + nb430_iz]

	;# store dr 
	movaps [rsp + nb430_dx], xmm0
	movaps [rsp + nb430_dy], xmm1
	movaps [rsp + nb430_dz], xmm2

    movd mm0, r8  ;# store j3
    movd mm1, r9
    movd mm2, r10
    movd mm3, r11

	;# square it 
	mulps xmm0,xmm0
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
    movaps xmm4, xmm0
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb430_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb430_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r
	movaps [rsp + nb430_r], xmm4
    movaps [rsp + nb430_rinv], xmm0
    
    movaps xmm8, xmm4    ;# r
	mulps xmm4, [rsp + nb430_gbscale] ;# rgbtab
    mulps xmm8, [rsp + nb430_tsc]    ;# rtab
    
    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm4  ;# gb
    cvttps2dq xmm9, xmm8  ;# lj
    
    ;# convert back to float
    cvtdq2ps  xmm6, xmm5   ;# gb
    cvtdq2ps  xmm10, xmm9  ;# lj
    
    ;# multiply by 4 and 8, respectively
    pslld   xmm5, 2   ;# gb
    pslld   xmm9, 3   ;# lj

    ;# move to integer registers
    movhlps xmm7, xmm5     ;# gb
    movhlps xmm11, xmm9    ;# lj
    movd    r8d, xmm5       ;# gb
    movd    r12d, xmm9      ;# lj
    movd    r10d, xmm7      ;# gb
    movd    r14d, xmm11     ;# lj
    pshufd  xmm5, xmm5, 1  ;# gb
    pshufd  xmm9, xmm9, 1  ;# lj
    pshufd  xmm7, xmm7, 1  ;# gb
    pshufd  xmm11, xmm11, 1 ;# lj
    movd    r9d, xmm5       ;# gb
    movd    r13d, xmm9      ;# lj
    movd    r11d, xmm7      ;# gb
    movd    r15d, xmm11     ;# lj
    ;# GB indices: r8-r11   LJ indices: r12-r15
    
    ;# calculate eps
    subps     xmm4, xmm6   ;# gb
    subps     xmm8, xmm10  ;# lj
    movaps    [rsp + nb430_epsgb], xmm4 ;# gb eps
    movaps    [rsp + nb430_eps], xmm8 ;# lj eps
    
	mov  rsi, [rbp + nb430_GBtab]
	mov  rdi, [rbp + nb430_VFtab]

    ;# load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
   	movlps xmm1, [rsi + r8*4]         ;# Y1c F1c 
   	movlps xmm5, [rdi + r12*4]        ;# Y1d F1d 
   	movlps xmm9, [rdi + r12*4 + 16]   ;# Y1r F1r 

	movlps xmm3, [rsi + r10*4]        ;# Y3c F3c 
	movlps xmm7, [rdi + r14*4]        ;# Y3d F3d 
	movlps xmm11, [rdi + r14*4 + 16]  ;# Y3r F3r 

	movhps xmm1, [rsi + r9*4]         ;# Y1c F1c Y2c F2c
	movhps xmm5, [rdi + r13*4]        ;# Y1d F1d Y2d F2d
	movhps xmm9, [rdi + r13*4 + 16]   ;# Y1r F1r Y2r F2r

	movhps xmm3, [rsi + r11*4]        ;# Y3c F3c Y4c F4c
	movhps xmm7, [rdi + r15*4]        ;# Y3d F3d Y4d F4d
	movhps xmm11, [rdi + r15*4 + 16]  ;# Y3r F3r Y4r F4r

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000   => Y1c Y2c Y3c Y4c
	shufps xmm4, xmm7, 136  ;# 10001000   => Y1d Y2d Y3d Y4d
	shufps xmm8, xmm11, 136  ;# 10001000  => Y1r Y2r Y3r Y4r
	shufps xmm1, xmm3, 221  ;# 11011101   => F1c F2c F3c F4c
	shufps xmm5, xmm7, 221  ;# 11011101   => F1d F2d F3d F4d
	shufps xmm9, xmm11, 221  ;# 11011101  => F1r F2r F3r F4r
    
   	movlps xmm3, [rsi + r8*4 + 8]      ;# G1c H1c 
   	movlps xmm7, [rdi + r12*4 + 8]     ;# G1d H1d 
   	movlps xmm11, [rdi + r12*4 + 24]   ;# G1r H1r 

	movlps xmm12, [rsi + r10*4 + 8]    ;# G3c H3c 
	movlps xmm13, [rdi + r14*4 + 8]    ;# G3d H3d 
	movlps xmm14, [rdi + r14*4 + 24]   ;# G3r H3r 

	movhps xmm3, [rsi + r9*4 + 8]      ;# G1c H1c G2c H2c
	movhps xmm7, [rdi + r13*4 + 8]     ;# G1d H1d G2d H2d
	movhps xmm11, [rdi + r13*4 + 24]   ;# G1r H1r G2r H2r

	movhps xmm12, [rsi + r11*4 + 8]    ;# G3c H3c G4c H4c
	movhps xmm13, [rdi + r15*4 + 8]    ;# G3d H3d G4d H4d
	movhps xmm14, [rdi + r15*4 + 24]   ;# G3r H3r G4r H4r
    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000  => G1c G2c G3c G4c
	shufps xmm6, xmm13, 136  ;# 10001000  => G1d G2d G3d G4d
	shufps xmm10, xmm14, 136  ;# 10001000 => G1r G2r G3r G4r
	shufps xmm3, xmm12, 221  ;# 11011101  => H1c H2c H3c H4c
	shufps xmm7, xmm13, 221  ;# 11011101  => H1d H2d H3d H4d
	shufps xmm11, xmm14, 221  ;# 11011101 => H1r H2r H3r H4r
    ;# table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps xmm12, [rsp + nb430_epsgb]
    movaps xmm13, [rsp + nb430_eps]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13
    mulps  xmm2, xmm12     ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm13
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm13
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, [rsp + nb430_qq]   ;# VV*qq = vcoul
    mulps  xmm5, [rsp + nb430_c6]   ;# vnb6
    mulps  xmm9, [rsp + nb430_c12]   ;# vnb12
    mulps  xmm3, [rsp + nb430_qq]    ;# FF*qq = fij
    mulps  xmm7, [rsp + nb430_c6]   ;# fijD
    mulps  xmm11, [rsp + nb430_c12]   ;#fijR

    addps  xmm11, xmm7 ;# fijD+fijR
    mulps  xmm11, [rsp + nb430_tsc] ;# (fijD+fijR)*tabscale
    
    ;# accumulate Vvdwtot
    addps  xmm5, [rsp + nb430_Vvdwtot]
    addps  xmm5, xmm9
    movaps [rsp + nb430_Vvdwtot], xmm5

	mov rsi, [rbp + nb430_dvda]
	
	;# Calculate dVda
	mulps xmm3, [rsp + nb430_gbscale]   ;# fijC=qq*FF*gbscale
	movaps xmm6, xmm3 
	mulps  xmm6, [rsp + nb430_r]
	addps  xmm6, xmm1   ;# vcoul+fijC*r

    addps  xmm3, xmm11  ;# fijC+fijD+fijR
    
    ;# increment vctot
	addps  xmm1, [rsp + nb430_vctot]
    movaps [rsp + nb430_vctot], xmm1

	;# xmm6=(vcoul+fijC*r)
	xorps  xmm7, xmm7
	subps  xmm7, xmm6
	movaps xmm6, xmm7
	
	;# update dvdasum 
	addps  xmm7, [rsp + nb430_dvdasum]
    movaps [rsp + nb430_dvdasum], xmm7

	;# update j atoms dvdaj
	movhlps xmm7, xmm6
	movaps  xmm5, xmm6
	movaps  xmm4, xmm7
	shufps  xmm5, xmm5, 0x1
	shufps  xmm4, xmm4, 0x1

	;# xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
	addss  xmm6, [rsi + rax*4]
	addss  xmm5, [rsi + rbx*4]
	addss  xmm7, [rsi + rcx*4]
	addss  xmm4, [rsi + rdx*4]
	movss  [rsi + rax*4], xmm6
	movss  [rsi + rbx*4], xmm5
	movss  [rsi + rcx*4], xmm7
	movss  [rsi + rdx*4], xmm4

	xorps  xmm4, xmm4	
	mulps xmm3, [rsp + nb430_rinv]
	subps  xmm4, xmm3

    movd r8, mm0   ;# fetch j3
    movd r9, mm1
    movd r10, mm2
    movd r11, mm3

    movaps  xmm9, xmm4
    movaps  xmm10, xmm4
    movaps  xmm11, xmm4
    
    mulps  xmm9, [rsp + nb430_dx]
    mulps  xmm10, [rsp + nb430_dy]
    mulps  xmm11, [rsp + nb430_dz]
    
	;# accumulate i forces
    movaps xmm12, [rsp + nb430_fix]
    movaps xmm13, [rsp + nb430_fiy]
    movaps xmm14, [rsp + nb430_fiz]
    addps xmm12, xmm9
    addps xmm13, xmm10
    addps xmm14, xmm11
    movaps [rsp + nb430_fix], xmm12
    movaps [rsp + nb430_fiy], xmm13
    movaps [rsp + nb430_fiz], xmm14

	mov rsi, [rbp + nb430_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + r8*4] ;# x1 y1 - -
	movlps xmm1, [rsi + r10*4] ;# x3 y3 - -
	movhps xmm0, [rsi + r9*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + r11*4] ;# x3 y3 x4 y4

    movaps xmm8, xmm9
    unpcklps xmm9, xmm10 ;# x1 y1 x2 y2
    unpckhps xmm8, xmm10 ;# x3 y3 x4 y4
    
    ;# update fjx and fjy
	addps  xmm0, xmm9
	addps  xmm1, xmm8
	
	movlps [rsi + r8*4], xmm0
	movlps [rsi + r10*4], xmm1
	movhps [rsi + r9*4], xmm0
	movhps [rsi + r11*4], xmm1
    
    ;# xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd  xmm10, xmm11, 1  ;# fjz2 - - -
    movhlps xmm9,  xmm11     ;# fjz3 - - -
    pshufd  xmm8,  xmm11, 3  ;# fjz4 - - -
    
	addss  xmm11, [rsi + r8*4 + 8]
	addss  xmm10, [rsi + r9*4 + 8]
	addss  xmm9,  [rsi + r10*4 + 8]
	addss  xmm8,  [rsi + r11*4 + 8]    
	movss  [rsi + r8*4 + 8], xmm11
	movss  [rsi + r9*4 + 8], xmm10
	movss  [rsi + r10*4 + 8], xmm9
	movss  [rsi + r11*4 + 8], xmm8
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb430_innerk],  4
	jl    .nb430_finish_inner
	jmp   .nb430_unroll_loop
.nb430_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb430_innerk],  4
	mov   edx, [rsp + nb430_innerk]
	and   edx, 2
	jnz   .nb430_dopair
	jmp   .nb430_checksingle
.nb430_dopair:	
	mov   rcx, [rsp + nb430_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb430_innerjjnr],  8

	;# load isaj
	mov rsi, [rbp + nb430_invsqrta]
	movss xmm3, [rsi + rax*4]
	movss xmm6, [rsi + rbx*4]
	movaps xmm2, [rsp + nb430_isai]
    unpcklps xmm3, xmm6
	mulps  xmm2, xmm3
    movaps [rsp + nb430_isaprod], xmm2	
    
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb430_gbtsc]
	movaps [rsp + nb430_gbscale], xmm1
	
	mov rsi, [rbp + nb430_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm6, [rsi + rbx*4]
    unpcklps xmm3, xmm6
	mulps xmm2, [rsp + nb430_iq]
	mulps  xmm3, xmm2
	movaps [rsp + nb430_qq], xmm3	
	
    ;# vdw parameters
	mov rsi, [rbp + nb430_type]
	mov r12d, [rsi + rax*4]
	mov r13d, [rsi + rbx*4]
	shl r12d, 1	
	shl r13d, 1	
    mov edi, [rsp + nb430_ntia]
	add r12d, edi
	add r13d, edi

	mov rsi, [rbp + nb430_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movhps xmm3, [rsi + r13*4]

    xorps xmm7, xmm7
	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101

    movaps [rsp + nb430_c6], xmm0
    movaps [rsp + nb430_c12], xmm3
    
	mov rsi, [rbp + nb430_pos]       ;# base of pos[] 
		
	lea   r8, [rax + rax*2]     ;# j3
	lea   r9, [rbx + rbx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm0, [rsi + r8*4]	;# x1 y1 - - 
	movlps xmm1, [rsi + r9*4]	;# x2 y2 - - 

	movss xmm2, [rsi + r8*4 + 8]	;# z1 - - - 
	movss xmm7, [rsi + r9*4 + 8]	;# z2 - - - 

    unpcklps xmm0, xmm1 ;# x1 x2 y1 y2
    movhlps  xmm1, xmm0 ;# y1 y2 -  -
    unpcklps xmm2, xmm7 ;# z1 z2 -  -
    
	;# calc dr 
	subps xmm0, [rsp + nb430_ix]
	subps xmm1, [rsp + nb430_iy]
	subps xmm2, [rsp + nb430_iz]

	;# store dr 
	movaps [rsp + nb430_dx], xmm0
	movaps [rsp + nb430_dy], xmm1
	movaps [rsp + nb430_dz], xmm2

	;# square it 
	mulps xmm0,xmm0
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
    movaps xmm4, xmm0
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb430_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb430_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r
	movaps [rsp + nb430_r], xmm4
    movaps [rsp + nb430_rinv], xmm0
    
    movaps xmm8, xmm4    ;# r
	mulps xmm4, [rsp + nb430_gbscale] ;# rgbtab
    mulps xmm8, [rsp + nb430_tsc]    ;# rtab
    
    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm4  ;# gb
    cvttps2dq xmm9, xmm8  ;# lj
    
    ;# convert back to float
    cvtdq2ps  xmm6, xmm5   ;# gb
    cvtdq2ps  xmm10, xmm9  ;# lj
    
    ;# multiply by 4 and 8, respectively
    pslld   xmm5, 2   ;# gb
    pslld   xmm9, 3   ;# lj

    ;# move to integer registers
    movd    r12d, xmm5       ;# gb
    movd    r14d, xmm9      ;# lj
    pshufd  xmm5, xmm5, 1   ;# gb
    pshufd  xmm9, xmm9, 1   ;# lj
    movd    r13d, xmm5       ;# gb
    movd    r15d, xmm9      ;# lj
    ;# GB indices: r12-r13   LJ indices: r14-r15
    
    ;# calculate eps
    subps     xmm4, xmm6   ;# gb
    subps     xmm8, xmm10  ;# lj
    movaps    [rsp + nb430_epsgb], xmm4 ;# gb eps
    movaps    [rsp + nb430_eps], xmm8 ;# lj eps
    
	mov  rsi, [rbp + nb430_GBtab]
	mov  rdi, [rbp + nb430_VFtab]

    ;# load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
   	movlps xmm0, [rsi + r12*4]       ;# Y1c F1c
	movlps xmm1, [rsi + r13*4]       ;# Y2c F2c
   	movlps xmm4, [rdi + r14*4]       ;# Y1d F1d  
	movlps xmm5, [rdi + r15*4]       ;# Y2d F2d
   	movlps xmm8, [rdi + r14*4 + 16]  ;# Y1r F1r
	movlps xmm9, [rdi + r15*4 + 16]  ;# Y2r F2r

    unpcklps xmm0, xmm1
    movhlps  xmm1, xmm0
    unpcklps xmm4, xmm5
    movhlps  xmm5, xmm4
    unpcklps xmm8, xmm9
    movhlps  xmm9, xmm8
   	movlps xmm2, [rsi + r12*4 + 8]    ;# G1c H1c
	movlps xmm3, [rsi + r13*4 + 8]    ;# G2c H2c
   	movlps xmm6, [rdi + r14*4 + 8]    ;# G1d H1d  
	movlps xmm7, [rdi + r15*4 + 8]    ;# G2d H2d
   	movlps xmm10, [rdi + r14*4 + 24]  ;# G1r H1r
	movlps xmm11, [rdi + r15*4 + 24]  ;# G2r H2r
    unpcklps xmm2, xmm3
    movhlps  xmm3, xmm2
    unpcklps xmm6, xmm7
    movhlps  xmm7, xmm6
    unpcklps xmm10, xmm11
    movhlps  xmm11, xmm10
    ;# table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps xmm12, [rsp + nb430_epsgb]
    movaps xmm13, [rsp + nb430_eps]
        
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13
    mulps  xmm2, xmm12     ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm13
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13

    addps  xmm1, xmm2   ;# F+Geps
    addps  xmm5, xmm6
    addps  xmm9, xmm10 
    addps  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addps  xmm5, xmm7
    addps  xmm9, xmm11 
    addps  xmm3, xmm3    ;# 2*Heps2
    addps  xmm7, xmm7
    addps  xmm11, xmm11
    addps  xmm3, xmm2    ;# 2*Heps2+Geps
    addps  xmm7, xmm6  
    addps  xmm11, xmm10
    addps  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm7, xmm5
    addps  xmm11, xmm9
    mulps  xmm1, xmm12   ;# eps*Fp
    mulps  xmm5, xmm13
    mulps  xmm9, xmm13
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, [rsp + nb430_qq]   ;# VV*qq = vcoul
    mulps  xmm5, [rsp + nb430_c6]   ;# vnb6
    mulps  xmm9, [rsp + nb430_c12]   ;# vnb12
    mulps  xmm3, [rsp + nb430_qq]    ;# FF*qq = fij
    mulps  xmm7, [rsp + nb430_c6]   ;# fijD
    mulps  xmm11, [rsp + nb430_c12]   ;#fijR

    addps  xmm11, xmm7 ;# fijD+fijR
    mulps  xmm11, [rsp + nb430_tsc] ;# (fijD+fijR)*tabscale
    
    ;# accumulate Vvdwtot
    addps  xmm5, [rsp + nb430_Vvdwtot]
    addps  xmm5, xmm9
    movlps [rsp + nb430_Vvdwtot], xmm5

	mov rsi, [rbp + nb430_dvda]
	
	;# Calculate dVda
	mulps xmm3, [rsp + nb430_gbscale]   ;# fijC=qq*FF*gbscale
	movaps xmm6, xmm3 
	mulps  xmm6, [rsp + nb430_r]
	addps  xmm6, xmm1   ;# vcoul+fijC*r

    addps  xmm3, xmm11  ;# fijC+fijD+fijR
    
    ;# increment vctot
	addps  xmm1, [rsp + nb430_vctot]
    movlps [rsp + nb430_vctot], xmm1

	;# xmm6=(vcoul+fijC*r)
	xorps  xmm7, xmm7
	subps  xmm7, xmm6
	movaps xmm6, xmm7
	
	;# update dvdasum 
	addps  xmm7, [rsp + nb430_dvdasum]
    movlps [rsp + nb430_dvdasum], xmm7

	;# update j atoms dvdaj
	movaps  xmm5, xmm6
	shufps  xmm5, xmm5, 0x1

	;# xmm6=dvdaj1 xmm5=dvdaj2 
	addss  xmm6, [rsi + rax*4]
	addss  xmm5, [rsi + rbx*4]
	movss  [rsi + rax*4], xmm6
	movss  [rsi + rbx*4], xmm5

	xorps  xmm4, xmm4	
	mulps xmm3, [rsp + nb430_rinv]
	subps  xmm4, xmm3

    movaps  xmm9, xmm4
    movaps  xmm10, xmm4
    movaps  xmm11, xmm4
    
    mulps  xmm9, [rsp + nb430_dx]
    mulps  xmm10, [rsp + nb430_dy]
    mulps  xmm11, [rsp + nb430_dz]
    
    
	;# accumulate i forces
    movaps xmm12, [rsp + nb430_fix]
    movaps xmm13, [rsp + nb430_fiy]
    movaps xmm14, [rsp + nb430_fiz]
    addps xmm12, xmm9
    addps xmm13, xmm10
    addps xmm14, xmm11
    movlps [rsp + nb430_fix], xmm12
    movlps [rsp + nb430_fiy], xmm13
    movlps [rsp + nb430_fiz], xmm14

	mov rsi, [rbp + nb430_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + r8*4] ;# x1 y1 - -
	movhps xmm0, [rsi + r9*4] ;# x1 y1 x2 y2

    unpcklps xmm9, xmm10  ;# x1 y1 x2 y2
    addps    xmm0, xmm9

	movlps [rsi + r8*4], xmm0
	movhps [rsi + r9*4], xmm0
    
    ;# z forces
    pshufd xmm8, xmm11, 1
    addss  xmm11, [rsi + r8*4 + 8] 
    addss  xmm8,  [rsi + r9*4 + 8]
    movss  [rsi + r8*4 + 8], xmm11
    movss  [rsi + r9*4 + 8], xmm8

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
	movss xmm3, [rsi + rax*4]
	movaps xmm2, [rsp + nb430_isai]
	mulss  xmm2, xmm3
    movaps [rsp + nb430_isaprod], xmm2	
    
	movaps xmm1, xmm2
	mulss xmm1, [rsp + nb430_gbtsc]
	movaps [rsp + nb430_gbscale], xmm1
	
	mov rsi, [rbp + nb430_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	mulss xmm2, [rsp + nb430_iq]
	mulss  xmm3, xmm2
	movaps [rsp + nb430_qq], xmm3	
	
    ;# vdw parameters
	mov rsi, [rbp + nb430_type]
	mov r12d, [rsi + rax*4]
	shl r12d, 1	
    mov edi, [rsp + nb430_ntia]
	add r12d, edi

	mov rsi, [rbp + nb430_vdwparam]
	movss xmm0, [rsi + r12*4]
	movss xmm3, [rsi + r12*4 + 4]
    movaps [rsp + nb430_c6], xmm0
    movaps [rsp + nb430_c12], xmm3
    
	mov rsi, [rbp + nb430_pos]       ;# base of pos[] 
		
	lea   r8, [rax + rax*2]     ;# j3

	;# move four coordinates to xmm0-xmm2 	
    movss  xmm0, [rsi + r8*4]
    movss  xmm1, [rsi + r8*4 + 4]
    movss  xmm2, [rsi + r8*4 + 8]
    
	;# calc dr 
	subss xmm0, [rsp + nb430_ix]
	subss xmm1, [rsp + nb430_iy]
	subss xmm2, [rsp + nb430_iz]

	;# store dr 
	movaps [rsp + nb430_dx], xmm0
	movaps [rsp + nb430_dy], xmm1
	movaps [rsp + nb430_dz], xmm2

	;# square it 
	mulss xmm0,xmm0
	mulss xmm1,xmm1
	mulss xmm2,xmm2
	addss xmm0, xmm1
	addss xmm0, xmm2
    movaps xmm4, xmm0
	;# rsq in xmm4 

	rsqrtss xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulss xmm5, xmm5
	movaps xmm1, [rsp + nb430_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb430_half]
	subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 
	mulss xmm4, xmm0	;# xmm4=r
	movaps [rsp + nb430_r], xmm4
    movaps [rsp + nb430_rinv], xmm0
    
    movaps xmm8, xmm4    ;# r
	mulss xmm4, [rsp + nb430_gbscale] ;# rgbtab
    mulss xmm8, [rsp + nb430_tsc]    ;# rtab
    
    ;# truncate and convert to integers
    cvttss2si r12d, xmm4  ;# gb
    cvttss2si r14d, xmm8  ;# lj
    
    ;# convert back to float
    cvtsi2ss  xmm6, r12d   ;# gb
    cvtsi2ss  xmm10, r14d  ;# lj
    
    ;# multiply by 4 and 8, respectively
    shl   r12d, 2   ;# gb
    shl   r14d, 3   ;# lj

    ;# GB index: r12   LJ indices: r14
    
    ;# calculate eps
    subss     xmm4, xmm6   ;# gb
    subss     xmm8, xmm10  ;# lj
    movaps    [rsp + nb430_epsgb], xmm4 ;# gb eps
    movaps    [rsp + nb430_eps], xmm8 ;# lj eps
    
	mov  rsi, [rbp + nb430_GBtab]
	mov  rdi, [rbp + nb430_VFtab]

    ;# load GB table data to xmm0-xmm3, disp to xmm4-xmm7, rep. to xmm8-xmm11
    movss  xmm0,  [rsi + r12*4]
    movss  xmm1,  [rsi + r12*4 + 4]
    movss  xmm2,  [rsi + r12*4 + 8]
    movss  xmm3,  [rsi + r12*4 + 12]
    movss  xmm4,  [rdi + r14*4]
    movss  xmm5,  [rdi + r14*4 + 4]
    movss  xmm6,  [rdi + r14*4 + 8]
    movss  xmm7,  [rdi + r14*4 + 12]
    movss  xmm8,  [rdi + r14*4 + 16]
    movss  xmm9,  [rdi + r14*4 + 20]
    movss  xmm10, [rdi + r14*4 + 24]
    movss  xmm11, [rdi + r14*4 + 28]
    ;# table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps xmm12, [rsp + nb430_epsgb]
    movaps xmm13, [rsp + nb430_eps]
    
    mulss  xmm3, xmm12   ;# Heps
    mulss  xmm7, xmm13
    mulss  xmm11, xmm13
    mulss  xmm2, xmm12     ;# Geps
    mulss  xmm6, xmm13
    mulss  xmm10, xmm13
    mulss  xmm3, xmm12   ;# Heps2
    mulss  xmm7, xmm13
    mulss  xmm11, xmm13

    addss  xmm1, xmm2   ;# F+Geps
    addss  xmm5, xmm6
    addss  xmm9, xmm10 
    addss  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addss  xmm5, xmm7
    addss  xmm9, xmm11 
    addss  xmm3, xmm3    ;# 2*Heps2
    addss  xmm7, xmm7
    addss  xmm11, xmm11
    addss  xmm3, xmm2    ;# 2*Heps2+Geps
    addss  xmm7, xmm6  
    addss  xmm11, xmm10
    addss  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addss  xmm7, xmm5
    addss  xmm11, xmm9
    mulss  xmm1, xmm12   ;# eps*Fp
    mulss  xmm5, xmm13
    mulss  xmm9, xmm13
    addss  xmm1, xmm0     ;# VV
    addss  xmm5, xmm4
    addss  xmm9, xmm8
    mulss  xmm1, [rsp + nb430_qq]   ;# VV*qq = vcoul
    mulss  xmm5, [rsp + nb430_c6]   ;# vnb6
    mulss  xmm9, [rsp + nb430_c12]   ;# vnb12
    mulss  xmm3, [rsp + nb430_qq]    ;# FF*qq = fij
    mulss  xmm7, [rsp + nb430_c6]   ;# fijD
    mulss  xmm11, [rsp + nb430_c12]   ;#fijR

    addss  xmm11, xmm7 ;# fijD+fijR
    mulss  xmm11, [rsp + nb430_tsc] ;# (fijD+fijR)*tabscale
    
    ;# accumulate Vvdwtot
    addss  xmm5, [rsp + nb430_Vvdwtot]
    addss  xmm5, xmm9
    movss [rsp + nb430_Vvdwtot], xmm5

	mov rsi, [rbp + nb430_dvda]
	
	;# Calculate dVda
	mulss xmm3, [rsp + nb430_gbscale]   ;# fijC=qq*FF*gbscale
	movaps xmm6, xmm3 
	mulss  xmm6, [rsp + nb430_r]
	addss  xmm6, xmm1   ;# vcoul+fijC*r

    addss  xmm3, xmm11  ;# fijC+fijD+fijR
    
    ;# increment vctot
	addss  xmm1, [rsp + nb430_vctot]
    movss [rsp + nb430_vctot], xmm1

	;# xmm6=(vcoul+fijC*r)
	xorps  xmm7, xmm7
	subss  xmm7, xmm6
	movaps xmm6, xmm7
	
	;# update dvdasum 
	addss  xmm7, [rsp + nb430_dvdasum]
    movss [rsp + nb430_dvdasum], xmm7

	;# update j atoms dvdaj

	;# xmm6=dvdaj1
	addss  xmm6, [rsi + rax*4]
	movss  [rsi + rax*4], xmm6

	xorps  xmm4, xmm4	
	mulss xmm3, [rsp + nb430_rinv]
	subss  xmm4, xmm3

    movss  xmm9, xmm4
    movss  xmm10, xmm4
    movss  xmm11, xmm4
    
    mulss  xmm9, [rsp + nb430_dx]
    mulss  xmm10, [rsp + nb430_dy]
    mulss  xmm11, [rsp + nb430_dz]
    
	;# accumulate i forces
    movaps xmm12, [rsp + nb430_fix]
    movaps xmm13, [rsp + nb430_fiy]
    movaps xmm14, [rsp + nb430_fiz]
    addss xmm12, xmm9
    addss xmm13, xmm10
    addss xmm14, xmm11
    movss [rsp + nb430_fix], xmm12
    movss [rsp + nb430_fiy], xmm13
    movss [rsp + nb430_fiz], xmm14

	mov rsi, [rbp + nb430_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + r8*4]
    addss  xmm10, [rsi + r8*4 + 4]
    addss  xmm11, [rsi + r8*4 + 8]
    movss  [rsi + r8*4],     xmm9
    movss  [rsi + r8*4 + 4], xmm10
    movss  [rsi + r8*4 + 8], xmm11
    
.nb430_updateouterdata:
	mov   ecx, [rsp + nb430_ii3]
	mov   rdi, [rbp + nb430_faction]
	mov   rsi, [rbp + nb430_fshift]
	mov   edx, [rsp + nb430_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb430_fix]
	movaps xmm1, [rsp + nb430_fiy]
	movaps xmm2, [rsp + nb430_fiz]

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
	movss  xmm3, [rdi + rcx*4]
	movss  xmm4, [rdi + rcx*4 + 4]
	movss  xmm5, [rdi + rcx*4 + 8]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4],     xmm3
	movss  [rdi + rcx*4 + 4], xmm4
	movss  [rdi + rcx*4 + 8], xmm5

	;# increment fshift force  
	movss  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 4]
	movss  xmm5, [rsi + rdx*4 + 8]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rsi + rdx*4],     xmm3
	movss  [rsi + rdx*4 + 4], xmm4
	movss  [rsi + rdx*4 + 8], xmm5

	;# get n from stack
	mov esi, [rsp + nb430_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb430_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb430_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb430_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb430_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb430_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate dVda and update it 
	movaps xmm7, [rsp + nb430_dvdasum]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
	
	mov edx, [rsp + nb430_ii]
	mov rax, [rbp + nb430_dvda]
	addss xmm7, [rax + rdx*4]
	movss [rax + rdx*4], xmm7
	
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

	add rsp, 552
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret





.globl nb_kernel430nf_x86_64_sse
.globl _nb_kernel430nf_x86_64_sse
nb_kernel430nf_x86_64_sse:	
_nb_kernel430nf_x86_64_sse:	
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
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb430nf_isai,           208
.equiv          nb430nf_isaprod,        224
.equiv          nb430nf_gbscale,        240
.equiv          nb430nf_r,              256
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
	movss xmm0, [rsi]
	movss [rsp + nb430nf_facel], xmm0

	mov rax, [rbp + nb430nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb430nf_tsc], xmm3

	mov rbx, [rbp + nb430nf_p_gbtabscale]
	movss xmm4, [rbx]
	shufps xmm4, xmm4, 0
	movaps [rsp + nb430nf_gbtsc], xmm4

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb430nf_half], eax
	movss xmm1, [rsp + nb430nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb430nf_half],  xmm1
	movaps [rsp + nb430nf_three],  xmm3

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
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb430nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb430nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb430nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb430nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb430nf_facel]
	shufps xmm3, xmm3, 0

	mov   rdx, [rbp + nb430nf_invsqrta]	;# load invsqrta[ii]
	movss xmm4, [rdx + rbx*4]
	shufps xmm4, xmm4, 0

    	mov   rdx, [rbp + nb430nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb430nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb430nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb430nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb430nf_iq], xmm3
	movaps [rsp + nb430nf_isai], xmm4
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb430nf_ix], xmm0
	movaps [rsp + nb430nf_iy], xmm1
	movaps [rsp + nb430nf_iz], xmm2

	mov   [rsp + nb430nf_ii3], ebx
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb430nf_vctot], xmm4
	movaps [rsp + nb430nf_Vvdwtot], xmm4
	
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
	sub   edx,  4
	add   ecx, [rsp + nb430nf_ninner]
	mov   [rsp + nb430nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb430nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb430nf_unroll_loop
	jmp   .nb430nf_finish_inner
.nb430nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb430nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb430nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	;# load isa2
	mov rsi, [rbp + nb430nf_invsqrta]
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]
	movaps xmm2, [rsp + nb430nf_isai]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm2, xmm3
	
	movaps [rsp + nb430nf_isaprod], xmm2
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb430nf_gbtsc]
	movaps [rsp + nb430nf_gbscale], xmm1
	
	mov rsi, [rbp + nb430nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	mulps xmm2, [rsp + nb430nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm3, xmm2
	movaps [rsp + nb430nf_qq], xmm3	

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov rsi, [rbp + nb430nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb430nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb430nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm6, xmm7, 221  ;# 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb430nf_c6], xmm4
	movaps [rsp + nb430nf_c12], xmm6
	
	mov rsi, [rbp + nb430nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	

	;# move four coordinates to xmm0-xmm2 	

	movlps xmm4, [rsi + rax*4]
	movlps xmm5, [rsi + rcx*4]
	movss xmm2, [rsi + rax*4 + 8]
	movss xmm6, [rsi + rcx*4 + 8]

	movhps xmm4, [rsi + rbx*4]
	movhps xmm5, [rsi + rdx*4]

	movss xmm0, [rsi + rbx*4 + 8]
	movss xmm1, [rsi + rdx*4 + 8]

	shufps xmm2, xmm0, 0
	shufps xmm6, xmm1, 0
	
	movaps xmm0, xmm4
	movaps xmm1, xmm4

	shufps xmm2, xmm6, 136  ;# 10001000
	
	shufps xmm0, xmm5, 136  ;# 10001000
	shufps xmm1, xmm5, 221  ;# 11011101		

	;# move ix-iz to xmm4-xmm6 
	movaps xmm4, [rsp + nb430nf_ix]
	movaps xmm5, [rsp + nb430nf_iy]
	movaps xmm6, [rsp + nb430nf_iz]

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
	movaps xmm1, [rsp + nb430nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb430nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r
	movaps [rsp + nb430nf_r], xmm4
	mulps xmm4, [rsp + nb430nf_gbscale]

	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 2
	pslld mm7, 2

	movd mm0, eax	
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx

	mov  rsi, [rbp + nb430nf_GBtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# load coulomb table
	movaps xmm4, [rsi + rax*4]
	movaps xmm5, [rsi + rbx*4]
	movaps xmm6, [rsi + rcx*4]
	movaps xmm7, [rsi + rdx*4]
	;# transpose, using xmm3 for scratch
	movaps xmm3, xmm6
	shufps xmm3, xmm7, 0xEE 
	shufps xmm6, xmm7, 0x44
	movaps xmm7, xmm4
	shufps xmm7, xmm5, 0xEE
	shufps xmm4, xmm5, 0x44
	movaps xmm5, xmm4
	shufps xmm5, xmm6, 0xDD
	shufps xmm4, xmm6, 0x88
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x88
	shufps xmm7, xmm3, 0xDD
	;# coulomb table ready, in xmm4-xmm7  		
	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	movaps xmm3, [rsp + nb430nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	addps  xmm5, [rsp + nb430nf_vctot]
	movaps [rsp + nb430nf_vctot], xmm5

	
	movaps xmm4, [rsp + nb430nf_r]
	mulps xmm4, [rsp + nb430nf_tsc]
	
	movhlps xmm5, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm5	;# mm6/mm7 contain lu indices 
	cvtpi2ps xmm6, mm6
	cvtpi2ps xmm5, mm7
	movlhps xmm6, xmm5
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	pslld mm7, 3
	
	mov  rsi, [rbp + nb430nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movaps xmm4, [rsi + rax*4]
	movaps xmm5, [rsi + rbx*4]
	movaps xmm6, [rsi + rcx*4]
	movaps xmm7, [rsi + rdx*4]
	;# transpose, using xmm3 for scratch
	movaps xmm3, xmm6
	shufps xmm3, xmm7, 0xEE 
	shufps xmm6, xmm7, 0x44
	movaps xmm7, xmm4
	shufps xmm7, xmm5, 0xEE
	shufps xmm4, xmm5, 0x44
	movaps xmm5, xmm4
	shufps xmm5, xmm6, 0xDD
	shufps xmm4, xmm6, 0x88
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x88
	shufps xmm7, xmm3, 0xDD
	;# dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, [rsp + nb430nf_c6]	 ;# Vvdw6
	addps  xmm5, [rsp + nb430nf_Vvdwtot]
	movaps [rsp + nb430nf_Vvdwtot], xmm5

	;# repulsion 
	movaps xmm4, [rsi + rax*4 + 16]
	movaps xmm5, [rsi + rbx*4 + 16]
	movaps xmm6, [rsi + rcx*4 + 16]
	movaps xmm7, [rsi + rdx*4 + 16]
	;# transpose, using xmm3 for scratch
	movaps xmm3, xmm6
	shufps xmm3, xmm7, 0xEE 
	shufps xmm6, xmm7, 0x44
	movaps xmm7, xmm4
	shufps xmm7, xmm5, 0xEE
	shufps xmm4, xmm5, 0x44
	movaps xmm5, xmm4
	shufps xmm5, xmm6, 0xDD
	shufps xmm4, xmm6, 0x88
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 0x88
	shufps xmm7, xmm3, 0xDD
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	mulps  xmm5, [rsp + nb430nf_c12] ;# Vvdw12
	addps  xmm5, [rsp + nb430nf_Vvdwtot]
	movaps [rsp + nb430nf_Vvdwtot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb430nf_innerk],  4
	jl    .nb430nf_finish_inner
	jmp   .nb430nf_unroll_loop
.nb430nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb430nf_innerk],  4
	mov   edx, [rsp + nb430nf_innerk]
	and   edx, 2
	jnz   .nb430nf_dopair
	jmp   .nb430nf_checksingle
.nb430nf_dopair:	

	mov   rcx, [rsp + nb430nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb430nf_innerjjnr],  8	

	xorps xmm2, xmm2
	movaps xmm6, xmm2
	
	;# load isa2
	mov rsi, [rbp + nb430nf_invsqrta]
	movss xmm2, [rsi + rax*4]
	movss xmm3, [rsi + rbx*4]
	unpcklps xmm2, xmm3	;# isa2 in xmm3(0,1)
	mulps  xmm2, [rsp + nb430nf_isai]
	movaps [rsp + nb430nf_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb430nf_gbtsc]
	movaps [rsp + nb430nf_gbscale], xmm1	
	
	mov rsi, [rbp + nb430nf_charge]    ;# base of charge[] 	
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	unpcklps xmm3, xmm6 ;# 00001000 ;# xmm3(0,1) has the charges 

	mulps  xmm2, [rsp + nb430nf_iq]
	mulps  xmm3, xmm2
	movaps [rsp + nb430nf_qq], xmm3

	mov rsi, [rbp + nb430nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]	
	mov rsi, [rbp + nb430nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb430nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [rsi + rcx*4]
	movhps xmm6, [rsi + rdx*4]
	mov rdi, [rbp + nb430nf_pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# 00001000 	
	shufps xmm6, xmm6, 13 ;# 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [rsp + nb430nf_c6], xmm4
	movaps [rsp + nb430nf_c12], xmm6	
	
	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# 10001000
	
	shufps xmm0, xmm0, 136  ;# 10001000
	shufps xmm1, xmm1, 221  ;# 11011101
	
	mov    rdi, [rbp + nb430nf_faction]
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb430nf_ix]
	movaps xmm5, [rsp + nb430nf_iy]
	movaps xmm6, [rsp + nb430nf_iz]

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
	movaps xmm1, [rsp + nb430nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb430nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	movaps [rsp + nb430nf_r], xmm4
	mulps xmm4, [rsp + nb430nf_gbscale]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb430nf_GBtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	;# load coulomb table
	movaps xmm4, [rsi + rcx*4]
	movaps xmm7, [rsi + rdx*4]
	;# transpose, using xmm3 for scratch
	movaps xmm6, xmm4
	unpcklps xmm4, xmm7  	;# Y1 Y2 F1 F2 
	unpckhps xmm6, xmm7     ;# G1 G2 H1 H2
	movhlps  xmm5, xmm4    	;# F1 F2 
	movhlps  xmm7, xmm6     ;# H1 H2
	;# coulomb table ready, in xmm4-xmm7  	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	movaps xmm3, [rsp + nb430nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	addps  xmm5, [rsp + nb430nf_vctot]
	movaps [rsp + nb430nf_vctot], xmm5 

	movaps xmm4, [rsp + nb430nf_r]
	mulps xmm4, [rsp + nb430nf_tsc]
	
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3
	
	mov  rsi, [rbp + nb430nf_VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6
	
	;# dispersion 
	movaps xmm4, [rsi + rcx*4]
	movaps xmm7, [rsi + rdx*4]
	;# transpose, using xmm3 for scratch
	movaps xmm6, xmm4
	unpcklps xmm4, xmm7  	;# Y1 Y2 F1 F2 
	unpckhps xmm6, xmm7     ;# G1 G2 H1 H2
	movhlps  xmm5, xmm4    	;# F1 F2 
	movhlps  xmm7, xmm6     ;# H1 H2
	;# dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	mulps  xmm5, [rsp + nb430nf_c6]	 ;# Vvdw6 
	addps  xmm5, [rsp + nb430nf_Vvdwtot]
	movaps [rsp + nb430nf_Vvdwtot], xmm5

	;# repulsion 
	movaps xmm4, [rsi + rcx*4 + 16]
	movaps xmm7, [rsi + rdx*4 + 16]
	;# transpose, using xmm3 for scratch
	movaps xmm6, xmm4
	unpcklps xmm4, xmm7  	;# Y1 Y2 F1 F2 
	unpckhps xmm6, xmm7     ;# G1 G2 H1 H2
	movhlps  xmm5, xmm4    	;# F1 F2 
	movhlps  xmm7, xmm6     ;# H1 H2
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	mulps  xmm5, [rsp + nb430nf_c12] ;# Vvdw12 
	
	addps  xmm5, [rsp + nb430nf_Vvdwtot]
	movaps [rsp + nb430nf_Vvdwtot], xmm5
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
	xorps  xmm2, xmm2
	movaps xmm6, xmm2
	movss xmm2, [rdx + rax*4]	;# isa2
	mulss xmm2, [rsp + nb430nf_isai]
	movss [rsp + nb430nf_isaprod], xmm2	
	movss xmm1, xmm2
	mulss xmm1, [rsp + nb430nf_gbtsc]
	movss [rsp + nb430nf_gbscale], xmm1	
	
	mulss  xmm2, [rsp + nb430nf_iq]
	movss xmm6, [rsi + rax*4]	;# xmm6(0) has the charge 	
	mulss  xmm6, xmm2
	movss [rsp + nb430nf_qq], xmm6
	
	mov rsi, [rbp + nb430nf_type]
	mov ecx, eax
	mov ecx, [rsi + rcx*4]	
	mov rsi, [rbp + nb430nf_vdwparam]
	shl ecx, 1
	add ecx, [rsp + nb430nf_ntia]
	movlps xmm6, [rsi + rcx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# 11111100	
	shufps xmm6, xmm6, 253  ;# 11111101	
	
	movss [rsp + nb430nf_c6], xmm4
	movss [rsp + nb430nf_c12], xmm6	
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	 
	
	movss xmm4, [rsp + nb430nf_ix]
	movss xmm5, [rsp + nb430nf_iy]
	movss xmm6, [rsp + nb430nf_iz]

	;# calc dr 
	subss xmm4, xmm0
	subss xmm5, xmm1
	subss xmm6, xmm2

	;# square it 
	mulss xmm4,xmm4
	mulss xmm5,xmm5
	mulss xmm6,xmm6
	addss xmm4, xmm5
	addss xmm4, xmm6
	;# rsq in xmm4 

	rsqrtss xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulss xmm5, xmm5
	movss xmm1, [rsp + nb430nf_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movss xmm0, [rsp + nb430nf_half]
	subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 

	mulss xmm4, xmm0	;# xmm4=r 
	movaps [rsp + nb430nf_r], xmm4
	mulss xmm4, [rsp + nb430nf_gbscale]

	cvttss2si ebx, xmm4     ;# mm6 contain lu indices 
	cvtsi2ss xmm6, ebx
	subss xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 2

	mov  rsi, [rbp + nb430nf_GBtab]
	
	movaps xmm4, [rsi + rbx*4]	
	movhlps xmm6, xmm4
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1
	shufps xmm7, xmm7, 1
	;# table ready in xmm4-xmm7 

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	movss xmm3, [rsp + nb430nf_qq]
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
	mulss  xmm5, xmm3 ;# vcoul=qq*VV  
	addss  xmm5, [rsp + nb430nf_vctot]
	movss [rsp + nb430nf_vctot], xmm5
	
	movss xmm4, [rsp + nb430nf_r]
	mulps xmm4, [rsp + nb430nf_tsc]
	
	cvttss2si ebx, xmm4
	cvtsi2ss xmm6, ebx
	subss xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 3
	mov  rsi, [rbp + nb430nf_VFtab]
	
	;# dispersion 
	movaps xmm4, [rsi + rbx*4]	
	movhlps xmm6, xmm4
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1
	shufps xmm7, xmm7, 1
	;# table ready in xmm4-xmm7 
	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
	mulss  xmm5, [rsp + nb430nf_c6]	 ;# Vvdw6
	addss  xmm5, [rsp + nb430nf_Vvdwtot]
	movss [rsp + nb430nf_Vvdwtot], xmm5

	;# repulsion 
	movaps xmm4, [rsi + rbx*4 + 16]	
	movhlps xmm6, xmm4
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1
	shufps xmm7, xmm7, 1
	;# table ready in xmm4-xmm7 
	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	mulss  xmm5, [rsp + nb430nf_c12] ;# Vvdw12 
	
	addss  xmm5, [rsp + nb430nf_Vvdwtot]
	movss [rsp + nb430nf_Vvdwtot], xmm5

.nb430nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb430nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb430nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb430nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb430nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb430nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb430nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
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



	

