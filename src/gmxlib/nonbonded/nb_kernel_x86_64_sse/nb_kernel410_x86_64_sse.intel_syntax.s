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




.globl nb_kernel410_x86_64_sse
.globl _nb_kernel410_x86_64_sse
nb_kernel410_x86_64_sse:	
_nb_kernel410_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb410_fshift,           16
.equiv          nb410_gid,              24
.equiv          nb410_pos,              32
.equiv          nb410_faction,          40
.equiv          nb410_charge,           48
.equiv          nb410_p_facel,          56
.equiv          nb410_argkrf,           64
.equiv          nb410_argcrf,           72
.equiv          nb410_Vc,               80
.equiv          nb410_type,             88
.equiv          nb410_p_ntype,          96
.equiv          nb410_vdwparam,         104
.equiv          nb410_Vvdw,             112
.equiv          nb410_p_tabscale,       120
.equiv          nb410_VFtab,            128
.equiv          nb410_invsqrta,         136
.equiv          nb410_dvda,             144
.equiv          nb410_p_gbtabscale,     152
.equiv          nb410_GBtab,            160
.equiv          nb410_p_nthreads,       168
.equiv          nb410_count,            176
.equiv          nb410_mtx,              184
.equiv          nb410_outeriter,        192
.equiv          nb410_inneriter,        200
.equiv          nb410_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb410_nri,              432
.equiv          nb410_iinr,             440
.equiv          nb410_jindex,           448
.equiv          nb410_jjnr,             456
.equiv          nb410_shift,            464
.equiv          nb410_shiftvec,         472
.equiv          nb410_facel,            480
.equiv          nb410_innerjjnr,        488
.equiv          nb410_is3,              496
.equiv          nb410_ii3,              500
.equiv          nb410_ii,               504
.equiv          nb410_ntia,             508
.equiv          nb410_innerk,           512
.equiv          nb410_n,                516
.equiv          nb410_nn1,              520
.equiv          nb410_ntype,            524
.equiv          nb410_nouter,           528
.equiv          nb410_ninner,           532
.equiv          nb410_jnra,             536
.equiv          nb410_jnrb,             540
.equiv          nb410_jnrc,             544
.equiv          nb410_jnrd,             548

	push rbp
	mov  rbp, rsp
	push rbx

	
	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 568		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb410_nouter], eax
	mov [rsp + nb410_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb410_nri], edi
	mov [rsp + nb410_iinr], rsi
	mov [rsp + nb410_jindex], rdx
	mov [rsp + nb410_jjnr], rcx
	mov [rsp + nb410_shift], r8
	mov [rsp + nb410_shiftvec], r9
	mov rdi, [rbp + nb410_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb410_ntype], edi
	mov rsi, [rbp + nb410_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb410_facel], xmm0

	mov rbx, [rbp + nb410_p_gbtabscale]
	movss xmm4, [rbx]
	shufps xmm4, xmm4, 0
	movaps [rsp + nb410_gbtsc], xmm4


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb410_half], eax
	movss xmm1, [rsp + nb410_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps xmm4, xmm3
	addps  xmm4, xmm4	;# six
	movaps xmm5, xmm4
	addps  xmm5, xmm5	;# twelve
	movaps [rsp + nb410_half],  xmm1
	movaps [rsp + nb410_two],  xmm2
	movaps [rsp + nb410_three],  xmm3
	movaps [rsp + nb410_six],  xmm4
	movaps [rsp + nb410_twelve],  xmm5

.nb410_threadloop:
        mov   rsi, [rbp + nb410_count]          ;# pointer to sync counter
        mov   eax, [rsi]
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
        mov ecx, [rsp + nb410_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb410_n], eax
        mov [rsp + nb410_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb410_outerstart
        jmp .nb410_end

.nb410_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb410_nouter]
	mov [rsp + nb410_nouter], ebx

.nb410_outer:
	mov   rax, [rsp + nb410_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb410_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb410_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb410_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 
	mov   [rsp + nb410_ii], ebx

	mov   rdx, [rbp + nb410_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb410_facel]
	shufps xmm3, xmm3, 0

	mov   rdx, [rbp + nb410_invsqrta]	;# load invsqrta[ii]
	movss xmm4, [rdx + rbx*4]
	shufps xmm4, xmm4, 0

    	mov   rdx, [rbp + nb410_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb410_ntype]
    	shl   edx, 1
    	mov   [rsp + nb410_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb410_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb410_iq], xmm3
	movaps [rsp + nb410_isai], xmm4
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb410_ix], xmm0
	movaps [rsp + nb410_iy], xmm1
	movaps [rsp + nb410_iz], xmm2

	mov   [rsp + nb410_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm13, xmm13
	movaps xmm12, xmm13
	movaps [rsp + nb410_Vvdwtot], xmm13
	movaps [rsp + nb410_dvdasum], xmm13
	movaps xmm14, xmm13
	movaps xmm15, xmm13
	
	mov   rax, [rsp + nb410_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb410_pos]
	mov   rdi, [rbp + nb410_faction]	
	mov   rax, [rsp + nb410_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb410_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb410_ninner]
	mov   [rsp + nb410_ninner], ecx
	add   edx, 0
	mov   [rsp + nb410_innerk], edx    ;# number of innerloop atoms 
	jge   .nb410_unroll_loop
	jmp   .nb410_finish_inner
.nb410_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb410_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb410_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	;# load isaj
	mov rsi, [rbp + nb410_invsqrta]
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]
	movaps xmm2, [rsp + nb410_isai]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all isaj in xmm3 
	mulps  xmm2, xmm3
	
	movaps [rsp + nb410_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb410_gbtsc]
	movaps [rsp + nb410_gbscale], xmm1
	
	mov rsi, [rbp + nb410_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	mulps xmm2, [rsp + nb410_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm3, xmm2
	movaps [rsp + nb410_qq], xmm3	
	
    ;# vdw parameters
	mov rsi, [rbp + nb410_type]
	mov r12d, [rsi + rax*4]
	mov r13d, [rsi + rbx*4]
	mov r14d, [rsi + rcx*4]
	mov r15d, [rsi + rdx*4]
	shl r12d, 1	
	shl r13d, 1	
	shl r14d, 1	
	shl r15d, 1	
    mov edi, [rsp + nb410_ntia]
	add r12d, edi
	add r13d, edi
	add r14d, edi
	add r15d, edi

	mov rsi, [rbp + nb410_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movlps xmm7, [rsi + r14*4]
	movhps xmm3, [rsi + r13*4]
	movhps xmm7, [rsi + r15*4]

	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101

    movaps [rsp + nb410_c6], xmm0
    movaps [rsp + nb410_c12], xmm3
    
	mov rsi, [rbp + nb410_pos]       ;# base of pos[] 
	
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
	subps xmm0, [rsp + nb410_ix]
	subps xmm1, [rsp + nb410_iy]
	subps xmm2, [rsp + nb410_iz]

	;# store dr 
	movaps [rsp + nb410_dx], xmm0
	movaps [rsp + nb410_dy], xmm1
	movaps [rsp + nb410_dz], xmm2

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
	movaps xmm1, [rsp + nb410_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb410_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r
	movaps [rsp + nb410_r], xmm4
	mulps xmm4, [rsp + nb410_gbscale]

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm4
    
    ;# convert back to float
    cvtdq2ps  xmm6, xmm5
    
    ;# multiply by 4
    pslld   xmm5, 2

    ;# move to integer registers
    movhlps xmm7, xmm5
    movd    r12d, xmm5
    movd    r14d, xmm7
    pshufd  xmm5, xmm5, 1
    pshufd  xmm7, xmm7, 1
    movd    r13d, xmm5
    movd    r15d, xmm7
    
    ;# calculate eps
    subps     xmm4, xmm6
    movaps    xmm1, xmm4 ;#eps
    
	mov  rsi, [rbp + nb410_GBtab]

    movaps xmm9, xmm0 ;# rinv
    mulps  xmm9, xmm9 ;# rinvsq
    movaps xmm10, xmm9 ;# rinvsq
    mulps  xmm10, xmm10 ;# rinv4
    mulps  xmm10, xmm9 ;# rinv6
    movaps xmm11, xmm10 
    mulps  xmm11, xmm11 ;# rinv12

    ;# load table data
   	movlps xmm5, [rsi + r12*4]
	movlps xmm7, [rsi + r14*4]
	movhps xmm5, [rsi + r13*4]
	movhps xmm7, [rsi + r15*4]

    movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

    mulps  xmm10, [rsp + nb410_c6]    ;# vvdw6=c6*rinv6
	mulps  xmm11, [rsp + nb410_c12]   ;# vvdw12=c12*rinv12     

	movaps xmm9, xmm11
	subps  xmm11, xmm10	;# Vvdw=Vvdw12-Vvdw6

    ;# add potential to vvdwtot 
	addps  xmm11, [rsp + nb410_Vvdwtot]
    movaps [rsp + nb410_Vvdwtot], xmm11
    
	movlps xmm7, [rsi + r12*4 + 8]   
	movlps xmm8, [rsi + r14*4 + 8]
	movhps xmm7, [rsi + r13*4 + 8]
	movhps xmm8, [rsi + r15*4 + 8]

    movaps xmm6, xmm7
    
	shufps xmm6, xmm8, 136  ;# 10001000
	shufps xmm7, xmm8, 221  ;# 11011101
    ;# table data ready in xmm4-xmm7

    mulps  xmm7, xmm1   ;# Heps
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm1	;# Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	addps  xmm7, xmm7	;# two*Heps2 
	movaps xmm3, [rsp + nb410_qq]
	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC

    ;# LJ forces
    mulps  xmm10, [rsp + nb410_six]
    mulps  xmm9, [rsp + nb410_twelve]
    subps  xmm9, xmm10
    mulps  xmm9, xmm0 ;# (12*vnb12-6*vnb6)*rinv

	mov rsi, [rbp + nb410_dvda]
	
	;# Calculate dVda
	xorps  xmm7, xmm7
	mulps xmm3, [rsp + nb410_gbscale]
	movaps xmm6, xmm3
	mulps  xmm6, [rsp + nb410_r]
	addps  xmm6, xmm5
    
    ;# increment vctot (sum in xmm12)
	addps  xmm12, xmm5

	;# xmm6=(vcoul+fijC*r)
	subps  xmm7, xmm6
	movaps xmm6, xmm7
	
    ;# update dvdasum
    addps  xmm7, [rsp + nb410_dvdasum]
    movaps [rsp + nb410_dvdasum], xmm7 

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

    subps  xmm9, xmm3
    mulps  xmm9, xmm0 ;# fscal

    movaps  xmm10, xmm9
    movaps  xmm11, xmm9

    mulps   xmm9, [rsp + nb410_dx]
    mulps   xmm10, [rsp + nb410_dy]
    mulps   xmm11, [rsp + nb410_dz]
    
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

	mov rsi, [rbp + nb410_faction]
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
	sub dword ptr [rsp + nb410_innerk],  4
	jl    .nb410_finish_inner
	jmp   .nb410_unroll_loop
.nb410_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb410_innerk],  4
	mov   edx, [rsp + nb410_innerk]
	and   edx, 2
	jnz   .nb410_dopair
	jmp   .nb410_checksingle
.nb410_dopair:	
	mov   rcx, [rsp + nb410_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb410_innerjjnr],  8

	;# load isaj
	mov rsi, [rbp + nb410_invsqrta]
	movss xmm2, [rsi + rax*4]
	movss xmm6, [rsi + rbx*4]
    unpcklps xmm2, xmm6

	mulps  xmm2, [rsp + nb410_isai]
	
	movaps [rsp + nb410_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb410_gbtsc]
	movaps [rsp + nb410_gbscale], xmm1
	
    mulps xmm2, [rsp + nb410_iq]
	mov rsi, [rbp + nb410_charge]    ;# base of charge[] 
	movss xmm3, [rsi + rax*4]
	movss xmm6, [rsi + rbx*4]
    unpcklps xmm3, xmm6
    

	mulps xmm3, xmm2
	movaps [rsp + nb410_qq], xmm3	
	
     ;# vdw parameters
	mov rsi, [rbp + nb410_type]
	mov r12d, [rsi + rax*4]
	mov r13d, [rsi + rbx*4]
	shl r12d, 1	
	shl r13d, 1	
    mov edi, [rsp + nb410_ntia]
	add r12d, edi
	add r13d, edi

	mov rsi, [rbp + nb410_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movhps xmm3, [rsi + r13*4]

    xorps xmm7, xmm7
	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101

    movaps [rsp + nb410_c6], xmm0
    movaps [rsp + nb410_c12], xmm3

	mov rsi, [rbp + nb410_pos]       ;# base of pos[] 
	
	lea   r8, [rax + rax*2]     ;# j3
	lea   r9, [rbx + rbx*2]	

	;# move four coordinates to xmm0-xmm2 	
	movlps xmm4, [rsi + r8*4]	;# x1 y1 - - 
	movlps xmm5, [rsi + r9*4]	;# x2 y2 - - 

	movss xmm6, [rsi + r8*4 + 8]	;# z1 - - - 
	movss xmm7, [rsi + r9*4 + 8]	;# z2 - - - 

    unpcklps xmm4, xmm5 ;# x1 x2 y1 y2
    movhlps  xmm5, xmm4 ;# y1 y2 -  -
    unpcklps xmm6, xmm7 ;# z1 z2 -  -
    
	;# calc dr 
	subps xmm4, [rsp + nb410_ix]
	subps xmm5, [rsp + nb410_iy]
	subps xmm6, [rsp + nb410_iz]

	;# store dr 
	movaps [rsp + nb410_dx], xmm4
	movaps [rsp + nb410_dy], xmm5
	movaps [rsp + nb410_dz], xmm6

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
	movaps xmm1, [rsp + nb410_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb410_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r
	movaps [rsp + nb410_r], xmm4
	mulps xmm4, [rsp + nb410_gbscale]

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm4
    
    ;# convert back to float
    cvtdq2ps  xmm6, xmm5
    
    ;# multiply by 4
    pslld   xmm5, 2

    ;# move to integer registers
    movd    r12d, xmm5
    pshufd  xmm5, xmm5, 1
    movd    r13d, xmm5
    
    ;# calculate eps
    subps     xmm4, xmm6
    movaps    xmm1, xmm4 ;#eps
    
	mov  rsi, [rbp + nb410_GBtab]

    movaps xmm9, xmm0 ;# rinv
    mulps  xmm9, xmm9 ;# rinvsq
    movaps xmm10, xmm9 ;# rinvsq
    mulps  xmm10, xmm10 ;# rinv4
    mulps  xmm10, xmm9 ;# rinv6
    movaps xmm11, xmm10 
    mulps  xmm11, xmm11 ;# rinv12

    ;# load table data
   	movlps xmm4, [rsi + r12*4]  ;# Y1 F1
	movlps xmm5, [rsi + r13*4]  ;# Y2 F2
    unpcklps xmm4, xmm5         ;# Y1 Y2 F1 F2
    movhlps  xmm5, xmm4         ;# F1 F2
    
    mulps  xmm10, [rsp + nb410_c6]    ;# vvdw6=c6*rinv6
	mulps  xmm11, [rsp + nb410_c12]   ;# vvdw12=c12*rinv12     

	movaps xmm9, xmm11
	subps  xmm11, xmm10	;# Vvdw=Vvdw12-Vvdw6

    ;# add potential to vvdwtot 
	addps  xmm11, [rsp + nb410_Vvdwtot]
    movlps [rsp + nb410_Vvdwtot], xmm11
    
   	movlps xmm6, [rsi + r12*4 + 8]   ;# G1 H1
	movlps xmm7, [rsi + r13*4 + 8]   ;# G2 H2
    unpcklps xmm6, xmm7              ;# G1 G2
    movhlps  xmm7, xmm6              ;# H1 H2
    ;# table data ready in xmm4-xmm7

    mulps  xmm7, xmm1   ;# Heps
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm1	;# Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	addps  xmm7, xmm7	;# two*Heps2 
	movaps xmm3, [rsp + nb410_qq]

	addps  xmm7, xmm6
	addps  xmm7, xmm5 ;# xmm7=FF 
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	mulps  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC

    ;# LJ forces
    mulps  xmm10, [rsp + nb410_six]
    mulps  xmm9,  [rsp + nb410_twelve]
    subps  xmm9, xmm10
    mulps  xmm9, xmm0 ;# (12*vnb12-6*vnb6)*rinv

    ;# zero upper part of vcoul 
    xorps xmm2, xmm2
    movlhps xmm5, xmm2
    
	mov rsi, [rbp + nb410_dvda]
	
	;# Calculate dVda
	xorps  xmm7, xmm7
	mulps xmm3, [rsp + nb410_gbscale]
	movaps xmm6, xmm3
	mulps  xmm6, [rsp + nb410_r]
	addps  xmm6, xmm5
    
    xorps  xmm4, xmm4
    ;# increment vctot (sum in xmm12)
	addps  xmm12, xmm5

	;# xmm6=(vcoul+fijC*r)
	subps  xmm7, xmm6
	movaps xmm6, xmm7

    ;# zero upper half of dvda
    movlhps xmm7, xmm4
    
    ;# update dvdasum
    addps  xmm7, [rsp + nb410_dvdasum]
    movaps [rsp + nb410_dvdasum], xmm7 

	;# update j atoms dvdaj
	movaps  xmm5, xmm6
	shufps  xmm5, xmm5, 0x1

	;# xmm6=dvdaj1 xmm5=dvdaj2 xmm7=dvdaj3 xmm4=dvdaj4
	addss  xmm6, [rsi + rax*4]
	addss  xmm5, [rsi + rbx*4]
	movss  [rsi + rax*4], xmm6
	movss  [rsi + rbx*4], xmm5

    xorps xmm7, xmm7
    
    subps  xmm9, xmm3
    mulps  xmm9, xmm0 ;# fscal

    movaps  xmm10, xmm9
    movaps  xmm11, xmm9

    mulps   xmm9, [rsp + nb410_dx]
    mulps   xmm10, [rsp + nb410_dy]
    mulps   xmm11, [rsp + nb410_dz]

    movlhps  xmm9, xmm7
    movlhps  xmm10, xmm7
    movlhps  xmm11, xmm7
    
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

	mov rsi, [rbp + nb410_faction]
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

.nb410_checksingle:				
	mov   edx, [rsp + nb410_innerk]
	and   edx, 1
	jnz    .nb410_dosingle
	jmp    .nb410_updateouterdata
.nb410_dosingle:
	mov rsi, [rbp + nb410_charge]
	mov rdx, [rbp + nb410_invsqrta]
	mov rdi, [rbp + nb410_pos]
	mov   rcx, [rsp + nb410_innerjjnr]
	mov   eax, [rcx]	

	;# load isaj
	mov rsi, [rbp + nb410_invsqrta]
	movss xmm3, [rsi + rax*4]
	movaps xmm2, [rsp + nb410_isai]
	mulss  xmm2, xmm3
	
	movss [rsp + nb410_isaprod], xmm2	
	movaps xmm1, xmm2
	mulss xmm1, [rsp + nb410_gbtsc]
	movss [rsp + nb410_gbscale], xmm1
	
    mulss xmm2, [rsp + nb410_iq]
	mov rsi, [rbp + nb410_charge]    ;# base of charge[] 

	movss xmm3, [rsi + rax*4]
	mulss xmm3, xmm2
	movss [rsp + nb410_qq], xmm3	
	
    ;# vdw parameters
	mov rsi, [rbp + nb410_type]
	mov r12d, [rsi + rax*4]
	shl r12d, 1	
    mov edi, [rsp + nb410_ntia]
	add r12d, edi

	mov rsi, [rbp + nb410_vdwparam]
	movss xmm0, [rsi + r12*4]
	movss xmm3, [rsi + r12*4 + 4]
    movaps [rsp + nb410_c6], xmm0
    movaps [rsp + nb410_c12], xmm3

	mov rsi, [rbp + nb410_pos]       ;# base of pos[] 
	
	lea   r8, [rax + rax*2]     ;# jnr 

	;# move four coordinates to xmm0-xmm2 	
	movss xmm4, [rsi + r8*4]	
	movss xmm5, [rsi + r8*4 + 4]	
	movss xmm6, [rsi + r8*4 + 8]
    
	;# calc dr 
	subss xmm4, [rsp + nb410_ix]
	subss xmm5, [rsp + nb410_iy]
	subss xmm6, [rsp + nb410_iz]

	;# store dr 
	movaps [rsp + nb410_dx], xmm4
	movaps [rsp + nb410_dy], xmm5
	movaps [rsp + nb410_dz], xmm6

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
	movaps xmm1, [rsp + nb410_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb410_half]
	subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 
	mulss xmm4, xmm0	;# xmm4=r
	movaps [rsp + nb410_r], xmm4
	mulss xmm4, [rsp + nb410_gbscale]

    ;# truncate and convert to integers
    cvttss2si r12d, xmm4
    
    ;# convert back to float
    cvtsi2ss  xmm6, r12d
    
    ;# multiply by 4
    shl   r12d, 2

    ;# calculate eps
    subss     xmm4, xmm6
    movaps    xmm1, xmm4 ;#eps
    
	mov  rsi, [rbp + nb410_GBtab]

    movaps xmm9, xmm0 ;# rinv
    mulss  xmm9, xmm9 ;# rinvsq
    movaps xmm10, xmm9 ;# rinvsq
    mulss  xmm10, xmm10 ;# rinv4
    mulss  xmm10, xmm9 ;# rinv6
    movaps xmm11, xmm10 
    mulss  xmm11, xmm11 ;# rinv12

    ;# load table data
   	movss xmm4, [rsi + r12*4]
	movss xmm5, [rsi + r12*4 + 4]
   	movss xmm6, [rsi + r12*4 + 8]
	movss xmm7, [rsi + r12*4 + 12]
    ;# table data ready in xmm4-xmm7

    mulss  xmm10, [rsp + nb410_c6]    ;# vvdw6=c6*rinv6
	mulss  xmm11, [rsp + nb410_c12]   ;# vvdw12=c12*rinv12     

	movaps xmm9, xmm11
	subss  xmm11, xmm10	;# Vvdw=Vvdw12-Vvdw6

    ;# add potential to vvdwtot 
	addss  xmm11, [rsp + nb410_Vvdwtot]
    movss [rsp + nb410_Vvdwtot], xmm11    

    mulss  xmm7, xmm1   ;# Heps
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm1	;# Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	addss  xmm7, xmm7	;# two*Heps2 
	movss  xmm3, [rsp + nb410_qq]
	addss  xmm7, xmm6
	addss  xmm7, xmm5 ;# xmm7=FF 
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
	mulss  xmm5, xmm3 ;# vcoul=qq*VV  
	mulss  xmm3, xmm7 ;# fijC=FF*qq 
	;# at this point xmm5 contains vcoul and xmm3 fijC

    ;# LJ forces
    mulss  xmm10, [rsp + nb410_six]
    mulss  xmm9,  [rsp + nb410_twelve]
    subss  xmm9, xmm10
    mulss  xmm9, xmm0 ;# (12*vnb12-6*vnb6)*rinv

	mov rsi, [rbp + nb410_dvda]
	
	;# Calculate dVda
	xorps  xmm7, xmm7
	mulss xmm3, [rsp + nb410_gbscale]
	movaps xmm6, xmm3
	mulss  xmm6, [rsp + nb410_r]
	addss  xmm6, xmm5
    
    ;# increment vctot (sum in xmm12)
	addss  xmm12, xmm5

	;# xmm6=(vcoul+fijC*r)
	subss  xmm7, xmm6
	movaps xmm6, xmm7

    ;# update dvdasum
    addss  xmm7, [rsp + nb410_dvdasum]
    movss [rsp + nb410_dvdasum], xmm7 

	;# update j atoms dvdaj
	addss  xmm6, [rsi + rax*4]
	movss  [rsi + rax*4], xmm6

    subss  xmm9, xmm3
    mulss  xmm9, xmm0 ;# fscal

    movaps  xmm10, xmm9
    movaps  xmm11, xmm9

    mulss   xmm9, [rsp + nb410_dx]
    mulss   xmm10, [rsp + nb410_dy]
    mulss   xmm11, [rsp + nb410_dz]
    
	;# accumulate i forces
    addss xmm13, xmm9
    addss xmm14, xmm10
    addss xmm15, xmm11

	mov rsi, [rbp + nb410_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + r8*4]
    addss  xmm10, [rsi + r8*4 + 4]
    addss  xmm11, [rsi + r8*4 + 8]
    movss  [rsi + r8*4],     xmm9
    movss  [rsi + r8*4 + 4], xmm10
    movss  [rsi + r8*4 + 8], xmm11
    
.nb410_updateouterdata:
	mov   ecx, [rsp + nb410_ii3]
	mov   rdi, [rbp + nb410_faction]
	mov   rsi, [rbp + nb410_fshift]
	mov   edx, [rsp + nb410_is3]

	;# accumulate i forces in xmm13, xmm14, xmm15
	movhlps xmm0, xmm13
	movhlps xmm1, xmm14
	movhlps xmm2, xmm15
	addps  xmm0, xmm13
	addps  xmm1, xmm14
	addps  xmm2, xmm15 
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
	mov esi, [rsp + nb410_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb410_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	;# accumulate 
	movhlps xmm6, xmm12
	addps  xmm12, xmm6	;# pos 0-1 in xmm12 have the sum now 
	movaps xmm6, xmm12
	shufps xmm6, xmm6, 1
	addss  xmm12, xmm6

	;# add earlier value from mem 
	mov   rax, [rbp + nb410_Vc]
	addss xmm12, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm12
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb410_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb410_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate dVda and update it 
	movaps xmm7, [rsp + nb410_dvdasum]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
	
	mov edx, [rsp + nb410_ii]
	mov rax, [rbp + nb410_dvda]
	addss xmm7, [rax + rdx*4]
	movss [rax + rdx*4], xmm7
	
        ;# finish if last 
        mov ecx, [rsp + nb410_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb410_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb410_n], esi
        jmp .nb410_outer
.nb410_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb410_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb410_end
        ;# non-zero, do one more workunit
        jmp   .nb410_threadloop
.nb410_end:

	mov eax, [rsp + nb410_nouter]
	mov ebx, [rsp + nb410_ninner]
	mov rcx, [rbp + nb410_outeriter]
	mov rdx, [rbp + nb410_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 568
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret



.globl nb_kernel410nf_x86_64_sse
.globl _nb_kernel410nf_x86_64_sse
nb_kernel410nf_x86_64_sse:	
_nb_kernel410nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb410nf_fshift,         16
.equiv          nb410nf_gid,            24
.equiv          nb410nf_pos,            32
.equiv          nb410nf_faction,        40
.equiv          nb410nf_charge,         48
.equiv          nb410nf_p_facel,        56
.equiv          nb410nf_argkrf,         64
.equiv          nb410nf_argcrf,         72
.equiv          nb410nf_Vc,             80
.equiv          nb410nf_type,           88
.equiv          nb410nf_p_ntype,        96
.equiv          nb410nf_vdwparam,       104
.equiv          nb410nf_Vvdw,           112
.equiv          nb410nf_p_tabscale,     120
.equiv          nb410nf_VFtab,          128
.equiv          nb410nf_invsqrta,       136
.equiv          nb410nf_dvda,           144
.equiv          nb410nf_p_gbtabscale,   152
.equiv          nb410nf_GBtab,          160
.equiv          nb410nf_p_nthreads,     168
.equiv          nb410nf_count,          176
.equiv          nb410nf_mtx,            184
.equiv          nb410nf_outeriter,      192
.equiv          nb410nf_inneriter,      200
.equiv          nb410nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb410nf_ix,             0
.equiv          nb410nf_iy,             16
.equiv          nb410nf_iz,             32
.equiv          nb410nf_iq,             48
.equiv          nb410nf_gbtsc,          64
.equiv          nb410nf_qq,             80
.equiv          nb410nf_c6,             96
.equiv          nb410nf_c12,            112
.equiv          nb410nf_vctot,          128
.equiv          nb410nf_Vvdwtot,        144
.equiv          nb410nf_half,           160
.equiv          nb410nf_three,          176
.equiv          nb410nf_isai,           192
.equiv          nb410nf_isaprod,        208
.equiv          nb410nf_gbscale,        224
.equiv          nb410nf_nri,            240
.equiv          nb410nf_iinr,           248
.equiv          nb410nf_jindex,         256
.equiv          nb410nf_jjnr,           264
.equiv          nb410nf_shift,          272
.equiv          nb410nf_shiftvec,       280
.equiv          nb410nf_facel,          288
.equiv          nb410nf_innerjjnr,      296
.equiv          nb410nf_is3,            304
.equiv          nb410nf_ii3,            308
.equiv          nb410nf_ntia,           312
.equiv          nb410nf_innerk,         316
.equiv          nb410nf_n,              320
.equiv          nb410nf_nn1,            324
.equiv          nb410nf_ntype,          328
.equiv          nb410nf_nouter,         332
.equiv          nb410nf_ninner,         336

	push rbp
	mov  rbp, rsp
	push rbx

	
	emms

        push r12
        push r13
        push r14
        push r15

	sub rsp, 360		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb410nf_nouter], eax
	mov [rsp + nb410nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb410nf_nri], edi
	mov [rsp + nb410nf_iinr], rsi
	mov [rsp + nb410nf_jindex], rdx
	mov [rsp + nb410nf_jjnr], rcx
	mov [rsp + nb410nf_shift], r8
	mov [rsp + nb410nf_shiftvec], r9
	mov rdi, [rbp + nb410nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb410nf_ntype], edi
	mov rsi, [rbp + nb410nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb410nf_facel], xmm0

	mov rbx, [rbp + nb410nf_p_gbtabscale]
	movss xmm4, [rbx]
	shufps xmm4, xmm4, 0
	movaps [rsp + nb410nf_gbtsc],  xmm4


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb410nf_half], eax
	movss xmm1, [rsp + nb410nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb410nf_half],  xmm1
	movaps [rsp + nb410nf_three],  xmm3

.nb410nf_threadloop:
        mov   rsi, [rbp + nb410nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
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
        mov ecx, [rsp + nb410nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb410nf_n], eax
        mov [rsp + nb410nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb410nf_outerstart
        jmp .nb410nf_end

.nb410nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb410nf_nouter]
	mov [rsp + nb410nf_nouter], ebx

.nb410nf_outer:
	mov   rax, [rsp + nb410nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb410nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb410nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb410nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii
	
	mov   rdx, [rbp + nb410nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb410nf_facel]
	shufps xmm3, xmm3, 0

	mov   rdx, [rbp + nb410nf_invsqrta]	;# load invsqrta[ii]
	movss xmm4, [rdx + rbx*4]
	shufps xmm4, xmm4, 0

    	mov   rdx, [rbp + nb410nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb410nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb410nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb410nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb410nf_iq], xmm3
	movaps [rsp + nb410nf_isai], xmm4
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb410nf_ix], xmm0
	movaps [rsp + nb410nf_iy], xmm1
	movaps [rsp + nb410nf_iz], xmm2

	mov   [rsp + nb410nf_ii3], ebx
	
	;# clear vctot
	xorps xmm4, xmm4
	movaps [rsp + nb410nf_vctot], xmm4
	movaps [rsp + nb410nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb410nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb410nf_pos]
	mov   rdi, [rbp + nb410nf_faction]	
	mov   rax, [rsp + nb410nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb410nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb410nf_ninner]
	mov   [rsp + nb410nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb410nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb410nf_unroll_loop
	jmp   .nb410nf_finish_inner
.nb410nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb410nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb410nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	;# load isa2
	mov rsi, [rbp + nb410nf_invsqrta]
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]
	movaps xmm2, [rsp + nb410nf_isai]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm2, xmm3
	
	movaps [rsp + nb410nf_isaprod], xmm2
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb410nf_gbtsc]
	movaps [rsp + nb410nf_gbscale], xmm1
	
	mov rsi, [rbp + nb410nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	mulps xmm2, [rsp + nb410nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm3, xmm2
	movaps [rsp + nb410nf_qq], xmm3	

	movd mm0, eax
	movd mm1, ebx
	movd mm2, ecx
	movd mm3, edx
	
	mov rsi, [rbp + nb410nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb410nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb410nf_ntia]
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

	movaps [rsp + nb410nf_c6], xmm4
	movaps [rsp + nb410nf_c12], xmm6
	
	mov rsi, [rbp + nb410nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [rsp + nb410nf_ix]
	movaps xmm5, [rsp + nb410nf_iy]
	movaps xmm6, [rsp + nb410nf_iz]

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
	movaps xmm1, [rsp + nb410nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb410nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb410nf_gbscale]

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

	mov  rsi, [rbp + nb410nf_GBtab]
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
	movaps xmm3, [rsp + nb410nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# update vctot
	addps  xmm5, [rsp + nb410nf_vctot]
	movaps [rsp + nb410nf_vctot], xmm5	
	
	;# L-J 
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	;# xmm4=rinvsq 

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	mulps  xmm6, xmm4	;# xmm6=rinvsix 
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulps  xmm6, [rsp + nb410nf_c6]
	mulps  xmm4, [rsp + nb410nf_c12]
	movaps xmm7, [rsp + nb410nf_Vvdwtot]
	addps  xmm7, xmm4
	subps  xmm7, xmm6
	movaps [rsp + nb410nf_Vvdwtot], xmm7
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb410nf_innerk],  4
	jl    .nb410nf_finish_inner
	jmp   .nb410nf_unroll_loop
.nb410nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb410nf_innerk],  4
	mov   edx, [rsp + nb410nf_innerk]
	and   edx, 2
	jnz   .nb410nf_dopair
	jmp   .nb410nf_checksingle
.nb410nf_dopair:	
	mov   rcx, [rsp + nb410nf_innerjjnr]
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb410nf_innerjjnr],  8

	xorps xmm2, xmm2
	movaps xmm6, xmm2
	
	;# load isa2
	mov rsi, [rbp + nb410nf_invsqrta]
	movss xmm2, [rsi + rax*4]
	movss xmm3, [rsi + rbx*4]
	unpcklps xmm2, xmm3	;# isa2 in xmm3(0,1)
	mulps  xmm2, [rsp + nb410nf_isai]
	movaps [rsp + nb410nf_isaprod], xmm2	
	movaps xmm1, xmm2
	mulps xmm1, [rsp + nb410nf_gbtsc]
	movaps [rsp + nb410nf_gbscale], xmm1	
	
	mov rsi, [rbp + nb410nf_charge]    ;# base of charge[] 	
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	unpcklps xmm3, xmm6 ;# 00001000 ;# xmm3(0,1) has the charges 

	mulps  xmm2, [rsp + nb410nf_iq]
	mulps  xmm3, xmm2
	movaps [rsp + nb410nf_qq], xmm3

	mov rsi, [rbp + nb410nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]	
	mov rsi, [rbp + nb410nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb410nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [rsi + rcx*4]
	movhps xmm6, [rsi + rdx*4]
	mov rdi, [rbp + nb410nf_pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# 00001000 	
	shufps xmm6, xmm6, 13 ;# 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [rsp + nb410nf_c6], xmm4
	movaps [rsp + nb410nf_c12], xmm6	
	
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
	
	mov    rdi, [rbp + nb410nf_faction]
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb410nf_ix]
	movaps xmm5, [rsp + nb410nf_iy]
	movaps xmm6, [rsp + nb410nf_iz]

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
	movaps xmm1, [rsp + nb410nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb410nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb410nf_gbscale]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb410nf_GBtab]
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
	movaps xmm3, [rsp + nb410nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  

	addps  xmm5, [rsp + nb410nf_vctot]
	movaps [rsp + nb410nf_vctot], xmm5
	
	;# L-J 
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	;# xmm4=rinvsq 

	;# at this point mm5 contains vcoul and mm3 fijC 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	mulps  xmm6, xmm4	;# xmm6=rinvsix 
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulps  xmm6, [rsp + nb410nf_c6]
	mulps  xmm4, [rsp + nb410nf_c12]
	movaps xmm7, [rsp + nb410nf_Vvdwtot]
	addps  xmm7, xmm4
	subps  xmm7, xmm6
	movaps [rsp + nb410nf_Vvdwtot], xmm7
	
.nb410nf_checksingle:				
	mov   edx, [rsp + nb410nf_innerk]
	and   edx, 1
	jnz    .nb410nf_dosingle
	jmp    .nb410nf_updateouterdata
.nb410nf_dosingle:
	mov rsi, [rbp + nb410nf_charge]
	mov rdx, [rbp + nb410nf_invsqrta]
	mov rdi, [rbp + nb410nf_pos]
	mov   rcx, [rsp + nb410nf_innerjjnr]
	mov   eax, [rcx]	
	xorps  xmm2, xmm2
	movaps xmm6, xmm2
	movss xmm2, [rdx + rax*4]	;# isa2
	mulss xmm2, [rsp + nb410nf_isai]
	movss [rsp + nb410nf_isaprod], xmm2	
	movss xmm1, xmm2
	mulss xmm1, [rsp + nb410nf_gbtsc]
	movss [rsp + nb410nf_gbscale], xmm1	
	
	mulss  xmm2, [rsp + nb410nf_iq]
	movss xmm6, [rsi + rax*4]	;# xmm6(0) has the charge 	
	mulss  xmm6, xmm2
	movss [rsp + nb410nf_qq], xmm6
	
	mov rsi, [rbp + nb410nf_type]
	mov ecx, eax
	mov ecx, [rsi + rcx*4]	
	mov rsi, [rbp + nb410nf_vdwparam]
	shl ecx, 1
	add ecx, [rsp + nb410nf_ntia]
	movlps xmm6, [rsi + rcx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# 11111100	
	shufps xmm6, xmm6, 253  ;# 11111101	
	
	movaps [rsp + nb410nf_c6], xmm4
	movaps [rsp + nb410nf_c12], xmm6	
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	 
	
	movaps xmm4, [rsp + nb410nf_ix]
	movaps xmm5, [rsp + nb410nf_iy]
	movaps xmm6, [rsp + nb410nf_iz]

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
	movss xmm1, [rsp + nb410nf_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movss xmm0, [rsp + nb410nf_half]
	subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 

	mulss xmm4, xmm0	;# xmm4=r 
	mulss xmm4, [rsp + nb410nf_gbscale]

	cvttss2si ebx, xmm4     ;# mm6 contain lu indices 
	cvtsi2ss xmm6, ebx
	subss xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 

	shl ebx, 2
	mov  rsi, [rbp + nb410nf_GBtab]
	
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
	movss xmm3, [rsp + nb410nf_qq]
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
	mulss  xmm5, xmm3 ;# vcoul=qq*VV  
	addss  xmm5, [rsp + nb410nf_vctot]
	movss [rsp + nb410nf_vctot], xmm5 	
	
	;# L-J 
	movaps xmm4, xmm0
	mulss  xmm4, xmm0	;# xmm4=rinvsq 

	movaps xmm6, xmm4
	mulss  xmm6, xmm4

	mulss  xmm6, xmm4	;# xmm6=rinvsix 
	movaps xmm4, xmm6
	mulss  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulss  xmm6, [rsp + nb410nf_c6]
	mulss  xmm4, [rsp + nb410nf_c12]
	movss xmm7, [rsp + nb410nf_Vvdwtot]
	addps  xmm7, xmm4
	subps  xmm7, xmm6
	movss [rsp + nb410nf_Vvdwtot], xmm7
	
.nb410nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb410nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb410nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb410nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb410nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb410nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb410nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb410nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb410nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb410nf_n], esi
        jmp .nb410nf_outer
.nb410nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb410nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb410nf_end
        ;# non-zero, do one more workunit
        jmp   .nb410nf_threadloop
.nb410nf_end:
	
	mov eax, [rsp + nb410nf_nouter]
	mov ebx, [rsp + nb410nf_ninner]
	mov rcx, [rbp + nb410nf_outeriter]
	mov rdx, [rbp + nb410nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 360
	emms


        pop r15
        pop r14
        pop r13
        pop r12

	pop rbx
	pop	rbp
	ret



