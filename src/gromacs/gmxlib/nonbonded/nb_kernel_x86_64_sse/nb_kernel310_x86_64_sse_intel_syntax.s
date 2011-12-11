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



.globl nb_kernel310_x86_64_sse
.globl _nb_kernel310_x86_64_sse
nb_kernel310_x86_64_sse:	
_nb_kernel310_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb310_fshift,           16
.equiv          nb310_gid,              24
.equiv          nb310_pos,              32
.equiv          nb310_faction,          40
.equiv          nb310_charge,           48
.equiv          nb310_p_facel,          56
.equiv          nb310_argkrf,           64
.equiv          nb310_argcrf,           72
.equiv          nb310_Vc,               80
.equiv          nb310_type,             88
.equiv          nb310_p_ntype,          96
.equiv          nb310_vdwparam,         104
.equiv          nb310_Vvdw,             112
.equiv          nb310_p_tabscale,       120
.equiv          nb310_VFtab,            128
.equiv          nb310_invsqrta,         136
.equiv          nb310_dvda,             144
.equiv          nb310_p_gbtabscale,     152
.equiv          nb310_GBtab,            160
.equiv          nb310_p_nthreads,       168
.equiv          nb310_count,            176
.equiv          nb310_mtx,              184
.equiv          nb310_outeriter,        192
.equiv          nb310_inneriter,        200
.equiv          nb310_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb310_ix,               0
.equiv          nb310_iy,               16
.equiv          nb310_iz,               32
.equiv          nb310_iq,               48
.equiv          nb310_dx,               64
.equiv          nb310_dy,               80
.equiv          nb310_dz,               96
.equiv          nb310_two,              112
.equiv          nb310_six,              128
.equiv          nb310_twelve,           144
.equiv          nb310_tsc,              160
.equiv          nb310_qq,               176
.equiv          nb310_c6,               192
.equiv          nb310_c12,              208
.equiv          nb310_fscal,            224
.equiv          nb310_vctot,            240
.equiv          nb310_Vvdwtot,          256
.equiv          nb310_fix,              272
.equiv          nb310_fiy,              288
.equiv          nb310_fiz,              304
.equiv          nb310_half,             320
.equiv          nb310_three,            336
.equiv          nb310_nri,              352
.equiv          nb310_iinr,             360
.equiv          nb310_jindex,           368
.equiv          nb310_jjnr,             376
.equiv          nb310_shift,            384
.equiv          nb310_shiftvec,         392
.equiv          nb310_facel,            400
.equiv          nb310_innerjjnr,        408
.equiv          nb310_is3,              416
.equiv          nb310_ii3,              420
.equiv          nb310_ntia,             424
.equiv          nb310_innerk,           428
.equiv          nb310_n,                432
.equiv          nb310_nn1,              436
.equiv          nb310_ntype,            440
.equiv          nb310_nouter,           444
.equiv          nb310_ninner,           448

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
	emms
	sub rsp, 464		;# local variable stack space (n*16+8)
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb310_nouter], eax
	mov [rsp + nb310_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb310_nri], edi
	mov [rsp + nb310_iinr], rsi
	mov [rsp + nb310_jindex], rdx
	mov [rsp + nb310_jjnr], rcx
	mov [rsp + nb310_shift], r8
	mov [rsp + nb310_shiftvec], r9
	mov rdi, [rbp + nb310_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb310_ntype], edi
	mov rsi, [rbp + nb310_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb310_facel], xmm0


	mov rax, [rbp + nb310_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb310_tsc], xmm3


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb310_half], eax
	movss xmm1, [rsp + nb310_half]
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
	movaps [rsp + nb310_half],  xmm1
	movaps [rsp + nb310_two],  xmm2
	movaps [rsp + nb310_three],  xmm3
	movaps [rsp + nb310_six],  xmm4
	movaps [rsp + nb310_twelve],  xmm5

.nb310_threadloop:
        mov   rsi, [rbp + nb310_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb310_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
	    cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb310_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb310_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb310_n], eax
        mov [rsp + nb310_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb310_outerstart
        jmp .nb310_end

.nb310_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb310_nouter]
	mov [rsp + nb310_nouter], ebx

.nb310_outer:
	mov   rax, [rsp + nb310_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb310_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb310_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb310_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb310_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb310_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb310_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb310_ntype]
    	shl   edx, 1
    	mov   [rsp + nb310_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb310_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb310_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb310_ix], xmm0
	movaps [rsp + nb310_iy], xmm1
	movaps [rsp + nb310_iz], xmm2

	mov   [rsp + nb310_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm15, xmm15
	movaps [rsp + nb310_vctot], xmm15
	movaps [rsp + nb310_Vvdwtot], xmm15
	movaps xmm14, xmm15
	movaps xmm13, xmm15
	
	mov   rax, [rsp + nb310_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb310_pos]
	mov   rdi, [rbp + nb310_faction]	
	mov   rax, [rsp + nb310_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb310_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb310_ninner]
	mov   [rsp + nb310_ninner], ecx
	add   edx, 0
	mov   [rsp + nb310_innerk], edx    ;# number of innerloop atoms 
	jge   .nb310_unroll_loop
	jmp   .nb310_finish_inner
.nb310_unroll_loop:	
	;# quad-unrolled innerloop here 
	mov   rdx, [rsp + nb310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	mov   r10d, [rdx + 8]            
	mov   r11d, [rdx + 12]         ;# eax-edx=jnr1-4 
    
	add qword ptr [rsp + nb310_innerjjnr],  16 ;# advance pointer (unrolled 4) 


	lea   rax, [r8 + r8*2]     ;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	

	lea   rcx, [r10 + r10*2]     ;# replace jnr with j3 
	lea   rdx, [r11 + r11*2]	

	;# load coordinates
	mov rdi, [rbp + nb310_pos]
    
	movlps xmm1, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + rbx*4]	;# x2 y2 - - 
	movlps xmm3, [rdi + rcx*4]	;# x3 y3 - -
	movlps xmm4, [rdi + rdx*4]	;# x4 y4 - -

	movss xmm5, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + rbx*4 + 8]	;# z2 - - - 
	movss xmm7, [rdi + rcx*4 + 8]	;# z3 - - - 
	movss xmm8, [rdi + rdx*4 + 8]	;# z4 - - - 

    unpcklps xmm1, xmm3 ;# x1 x3 y1 y3
    unpcklps xmm2, xmm4 ;# x2 x4 y2 y4
    unpcklps xmm5, xmm7 ;# z1 z3 -  -
    unpcklps xmm6, xmm8 ;# z2 z4 -  -
	mov rsi, [rbp + nb310_charge]

    movaps xmm3, xmm1

    unpcklps xmm1, xmm2 ;# x1 x2 x3 x4
    unpckhps xmm3, xmm2 ;# y1 y2 y3 y4
    unpcklps xmm5, xmm6 ;# z1 z2 z3 z4

	;# calc dr  
	subps xmm1, [rsp + nb310_ix]
	subps xmm3, [rsp + nb310_iy]
	subps xmm5, [rsp + nb310_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm1
    movaps xmm10, xmm3
    movaps xmm11, xmm5
    
	movss xmm0, [rsi + r8*4]
	movss xmm2, [rsi + r10*4]
	movss xmm6, [rsi + r9*4]
	movss xmm8, [rsi + r11*4]

	;# square it 
	mulps xmm1, xmm1
	mulps xmm3, xmm3
	mulps xmm5, xmm5
	addps xmm3, xmm1
	addps xmm3, xmm5
	;# rsq in xmm3
	mov rsi, [rbp + nb310_type]
    
    unpcklps xmm0, xmm2
    unpcklps xmm6, xmm8
    
    unpcklps xmm0, xmm6
    

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm3
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb310_three]
	mulps xmm5, xmm3	;# rsq*lu*lu 	
    subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm1, [rsp + nb310_half]	
    ;# xmm1=rinv
    ;# xmm3=rsq
	mulps xmm0, [rsp + nb310_iq] 
    
    ;# vdw types
	mov r8d, [rsi + r8*4]
	mov r9d, [rsi + r9*4]
	mov r10d, [rsi + r10*4]
	mov r11d, [rsi + r11*4]

    mulps xmm3, xmm1 ;# r
    mulps xmm3, [rsp + nb310_tsc] ;# rtab
    movaps [rsp + nb310_qq], xmm0

    ;# truncate and convert to integers
    cvttps2dq xmm2, xmm3
    
	shl r8d, 1	
	shl r9d, 1	
	shl r10d, 1	
	shl r11d, 1	

    ;# convert back to float
    cvtdq2ps  xmm0, xmm2
    
    mov edi, [rsp + nb310_ntia]
	add r8d, edi
	add r9d, edi
	add r10d, edi
	add r11d, edi

    ;# multiply by 4
    pslld   xmm2, 2

    ;# move to integer registers
    movhlps xmm7, xmm2
    movd    r12d, xmm2
    movd    r14d, xmm7
    pshufd  xmm2, xmm2, 1
    pshufd  xmm7, xmm7, 1
    movd    r13d, xmm2
    movd    r15d, xmm7
    

    ;# calculate eps
    subps     xmm3, xmm0

	mov rsi, [rbp + nb310_vdwparam]
	movlps xmm7, [rsi + r8*4]
	movlps xmm8, [rsi + r10*4]
	movhps xmm7, [rsi + r9*4]
	movhps xmm8, [rsi + r11*4]

	movaps xmm12, xmm7
	shufps xmm12, xmm8, 136  ;# 10001000
	shufps xmm7, xmm8, 221  ;# 11011101

    movaps [rsp + nb310_c6], xmm12
    movaps [rsp + nb310_c12], xmm7

	mov rsi, [rbp + nb310_VFtab]
    ;# load table data
   	movlps xmm5, [rsi + r12*4]
	movlps xmm7, [rsi + r14*4]
	movhps xmm5, [rsi + r13*4]
	movhps xmm7, [rsi + r15*4]

    movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
    
    movaps xmm0, xmm1 ;# rinv
    mulps  xmm0, xmm0 ;# rinvsq
    movaps xmm2, xmm0 ;# rinvsq
    mulps  xmm2, xmm2 ;# rinv4
    mulps  xmm2, xmm0 ;# rinv6
    movaps xmm12, xmm2 
    mulps  xmm12, xmm12 ;# rinv12
    
	movlps xmm7, [rsi + r12*4 + 8]   
	movlps xmm8, [rsi + r14*4 + 8]
	movhps xmm7, [rsi + r13*4 + 8]
	movhps xmm8, [rsi + r15*4 + 8]

    movaps xmm6, xmm7

    mulps  xmm2, [rsp + nb310_c6]    ;# vvdw6=c6*rinv6
	mulps  xmm12, [rsp + nb310_c12]   ;# vvdw12=c12*rinv12     

	movaps xmm0, xmm12
	subps  xmm12, xmm2	;# Vvdw=Vvdw12-Vvdw6

    ;# add potential to vvdwtot 
	addps  xmm12, [rsp + nb310_Vvdwtot]
    movaps [rsp + nb310_Vvdwtot], xmm12
    
	shufps xmm6, xmm8, 136  ;# 10001000
	shufps xmm7, xmm8, 221  ;# 11011101
    ;# table data ready in xmm4-xmm7
    
    mulps xmm7, xmm3   ;# Heps
    mulps  xmm6, xmm3  ;# Geps
    mulps xmm7, xmm3   ;# Heps2

    addps  xmm5, xmm6   ;# F+Geps
    addps  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addps  xmm7, xmm7   ;# 2*Heps2
    addps  xmm7, xmm6   ;# 2*Heps2+Geps
    addps  xmm7, xmm5   ;# FF = Fp + 2*Heps2 + Geps
    mulps  xmm5, xmm3   ;# eps*Fp
    addps  xmm5, xmm4   ;# VV
    mulps  xmm5, [rsp + nb310_qq]   ;# VV*qq=vcoul
    mulps  xmm7, [rsp + nb310_qq]   ;# FF*qq=fijC

    ;# LJ forces
    mulps  xmm2, [rsp + nb310_six]
    mulps  xmm0, [rsp + nb310_twelve]
    subps  xmm0, xmm2
    mulps  xmm0, xmm1 ;# (12*vnb12-6*vnb6)*rinv

    ;# add potential to vctot 
	addps  xmm5, [rsp + nb310_vctot]
    movaps [rsp + nb310_vctot], xmm5

    mulps  xmm7, [rsp + nb310_tsc]
    subps  xmm0, xmm7
    
    mulps  xmm0, xmm1  ;# fscal

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm0
	mulps  xmm10, xmm0
	mulps  xmm11, xmm0

	mov rsi, [rbp + nb310_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movlps xmm1, [rsi + rcx*4] ;# x3 y3 - -
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + rdx*4] ;# x3 y3 x4 y4

	;# xmm0-xmm2 contains tx-tz (partial force) 
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

    movaps xmm8, xmm9
    unpcklps xmm9, xmm10 ;# x1 y1 x2 y2
    unpckhps xmm8, xmm10 ;# x3 y3 x4 y4
    
    ;# update fjx and fjy
	addps  xmm0, xmm9
	addps  xmm1, xmm8
	
	movlps [rsi + rax*4], xmm0
	movlps [rsi + rcx*4], xmm1
	movhps [rsi + rbx*4], xmm0
	movhps [rsi + rdx*4], xmm1
    
    ;# xmm11: fjz1 fjz2 fjz3 fjz4
    pshufd  xmm10, xmm11, 1  ;# fjz2 - - -
    movhlps xmm9,  xmm11     ;# fjz3 - - -
    pshufd  xmm8,  xmm11, 3  ;# fjz4 - - -
    
	addss  xmm11, [rsi + rax*4 + 8]
	addss  xmm10, [rsi + rbx*4 + 8]
	addss  xmm9,  [rsi + rcx*4 + 8]
	addss  xmm8,  [rsi + rdx*4 + 8]    
	movss  [rsi + rax*4 + 8], xmm11
	movss  [rsi + rbx*4 + 8], xmm10
	movss  [rsi + rcx*4 + 8], xmm9
	movss  [rsi + rdx*4 + 8], xmm8
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb310_innerk],  4
	jl    .nb310_finish_inner
	jmp   .nb310_unroll_loop
.nb310_finish_inner:
    ;# check if at least two particles remain 
    add dword ptr [rsp + nb310_innerk],  4
    mov   edx, [rsp + nb310_innerk]
    and   edx, 2
    jnz   .nb310_dopair
    jmp   .nb310_checksingle
.nb310_dopair:  
	;# twice-unrolled innerloop here 
	mov   rdx, [rsp + nb310_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
    
	add qword ptr [rsp + nb310_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb310_charge]
	movss xmm0, [rsi + rax*4]
	movss xmm2, [rsi + rbx*4]

    unpcklps xmm0, xmm2  ;# jqa jqb 
	mulps xmm0, [rsp + nb310_iq] 
    movaps [rsp + nb310_qq], xmm0

	mov rsi, [rbp + nb310_type]
    ;# vdw parameters
	mov r12d, [rsi + rax*4]
	mov r13d, [rsi + rbx*4]
	shl r12d, 1	
	shl r13d, 1	
    mov edi, [rsp + nb310_ntia]
	add r12d, edi
	add r13d, edi

	mov rsi, [rbp + nb310_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movhps xmm3, [rsi + r13*4]

    xorps  xmm7, xmm7
	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101

    movaps [rsp + nb310_c6], xmm0
    movaps [rsp + nb310_c12], xmm3

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# load coordinates
	mov rdi, [rbp + nb310_pos]
    
	movlps xmm4, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm5, [rdi + rbx*4]	;# x2 y2 - - 

	movss xmm6, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm7, [rdi + rbx*4 + 8]	;# z2 - - - 

    unpcklps xmm4, xmm5 ;# x1 x2 y1 y2
    movhlps  xmm5, xmm4 ;# y1 y2 -  -
    unpcklps xmm6, xmm7 ;# z1 z2 -  -

	;# calc dr  
	subps xmm4, [rsp + nb310_ix]
	subps xmm5, [rsp + nb310_iy]
	subps xmm6, [rsp + nb310_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm4
    movaps xmm10, xmm5
    movaps xmm11, xmm6
    
	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm4, xmm5
	addps xmm4, xmm6
	;# rsq in xmm4 

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm4
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb310_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 	
    subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm1, [rsp + nb310_half]		
    ;# xmm1=rinv
    movaps xmm3, xmm4
    ;# xmm3=rsq 

    mulps xmm3, xmm1 ;# r
    mulps xmm3, [rsp + nb310_tsc] ;# rtab
    
    ;# truncate and convert to integers
    cvttps2dq xmm2, xmm3
    
    ;# convert back to float
    cvtdq2ps  xmm0, xmm2
    
    ;# multiply by 4
    pslld   xmm2, 2

    ;# move to integer registers
    movd    r12d, xmm2
    pshufd  xmm2, xmm2, 1
    movd    r13d, xmm2
    
    ;# calculate eps
    subps     xmm3, xmm0

	mov rsi, [rbp + nb310_VFtab]
    ;# load table data
   	movlps xmm4, [rsi + r12*4]
	movlps xmm5, [rsi + r13*4]
    unpcklps xmm4, xmm5
    movhlps xmm5, xmm4
    
    movaps xmm0, xmm1 ;# rinv
    mulps  xmm0, xmm0 ;# rinvsq
    movaps xmm2, xmm0 ;# rinvsq
    mulps  xmm2, xmm2 ;# rinv4
    mulps  xmm2, xmm0 ;# rinv6
    movaps xmm12, xmm2 
    mulps  xmm12, xmm12 ;# rinv12

   	movlps xmm6, [rsi + r12*4 + 8]
	movlps xmm7, [rsi + r13*4 + 8]
    unpcklps xmm6, xmm7
    movhlps xmm7, xmm6
    ;# table data ready in xmm4-xmm7

    mulps  xmm2, [rsp + nb310_c6]    ;# vvdw6=c6*rinv6
	mulps  xmm12, [rsp + nb310_c12]   ;# vvdw12=c12*rinv12     

	movaps xmm0, xmm12
	subps  xmm12, xmm2	;# Vvdw=Vvdw12-Vvdw6

    ;# add potential to vvdwtot 
	addps  xmm12, [rsp + nb310_Vvdwtot]
    movlps [rsp + nb310_Vvdwtot], xmm12
    
    mulps xmm7, xmm3   ;# Heps
    mulps  xmm6, xmm3  ;# Geps
    mulps xmm7, xmm3   ;# Heps2

    addps  xmm5, xmm6   ;# F+Geps
    addps  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addps  xmm7, xmm7   ;# 2*Heps2
    addps  xmm7, xmm6   ;# 2*Heps2+Geps
    addps  xmm7, xmm5   ;# FF = Fp + 2*Heps2 + Geps
    mulps  xmm5, xmm3   ;# eps*Fp
    addps  xmm5, xmm4   ;# VV
    mulps  xmm5, [rsp + nb310_qq]   ;# VV*qq=vcoul
    mulps  xmm7, [rsp + nb310_qq]   ;# FF*qq=fijC

    ;# LJ forces
    mulps  xmm2, [rsp + nb310_six]
    mulps  xmm0, [rsp + nb310_twelve]
    subps  xmm0, xmm2
    mulps  xmm0, xmm1 ;# (12*vnb12-6*vnb6)*rinv

    ;# add potential to vctot 
	addps  xmm5, [rsp + nb310_vctot]
    movlps [rsp + nb310_vctot], xmm5

    xorps xmm8, xmm8
    
    mulps  xmm7, [rsp + nb310_tsc]
    subps  xmm0, xmm7
    
    mulps  xmm0, xmm1  ;# fscal

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm0
	mulps  xmm10, xmm0
	mulps  xmm11, xmm0

    movlhps xmm9, xmm8
    movlhps xmm10, xmm8
    movlhps xmm11, xmm8
    
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

	mov rsi, [rbp + nb310_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2

    unpcklps xmm9, xmm10  ;# x1 y1 x2 y2
    addps    xmm0, xmm9

	movlps [rsi + rax*4], xmm0
	movhps [rsi + rbx*4], xmm0
    
    ;# z forces
    pshufd xmm8, xmm11, 1
    addss  xmm11, [rsi + rax*4 + 8] 
    addss  xmm8,  [rsi + rbx*4 + 8]
    movss  [rsi + rax*4 + 8], xmm11
    movss  [rsi + rbx*4 + 8], xmm8

.nb310_checksingle:                             
    mov   edx, [rsp + nb310_innerk]
    and   edx, 1
    jnz    .nb310_dosingle
    jmp    .nb310_updateouterdata

.nb310_dosingle:	
    mov rcx, [rsp + nb310_innerjjnr]
	mov   eax, [rcx]	            

	mov rsi, [rbp + nb310_charge]
	movss xmm0, [rsi + rax*4]

	mulss xmm0, [rsp + nb310_iq] 
    movaps [rsp + nb310_qq], xmm0

	mov rsi, [rbp + nb310_type]
    ;# vdw parameters
	mov r12d, [rsi + rax*4]
	shl r12d, 1	
    mov edi, [rsp + nb310_ntia]
	add r12d, edi

	mov rsi, [rbp + nb310_vdwparam]
	movss xmm0, [rsi + r12*4]
	movss xmm3, [rsi + r12*4 + 4]

    movaps [rsp + nb310_c6], xmm0
    movaps [rsp + nb310_c12], xmm3

	lea   rax, [rax + rax*2]        ;# replace jnr with j3 

	mov rdi, [rbp + nb310_pos]
	movss xmm4, [rdi + rax*4]	    ;# x1 - - - 
	movss xmm5, [rdi + rax*4 + 4]    ;# y2 - - - 
	movss xmm6, [rdi + rax*4 + 8]    ;# 13 - - - 

	;# calc dr  
	subss xmm4, [rsp + nb310_ix]
	subss xmm5, [rsp + nb310_iy]
	subss xmm6, [rsp + nb310_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm4
    movaps xmm10, xmm5
    movaps xmm11, xmm6
    
	;# square it 
	mulss xmm4,xmm4
	mulss xmm5,xmm5
	mulss xmm6,xmm6
	addss xmm4, xmm5
	addss xmm4, xmm6
	;# rsq in xmm4 

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtss xmm5, xmm4
	movaps xmm2, xmm5
	mulss xmm5, xmm5
	movaps xmm1, [rsp + nb310_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 	
    subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm1, [rsp + nb310_half]		
    ;# xmm1=rinv
    movaps xmm3, xmm4
    ;# xmm3=rsq 

    mulss xmm3, xmm1 ;# r
    mulss xmm3, [rsp + nb310_tsc] ;# rtab
    
    ;# truncate and convert to integers
    cvttss2si r12d, xmm3
    
    ;# convert back to float
    cvtsi2ss  xmm0, r12d
    
    ;# multiply by 4
    shl       r12d, 2

    ;# calculate eps
    subss     xmm3, xmm0

	mov rsi, [rbp + nb310_VFtab]
    
    movaps xmm0, xmm1 ;# rinv
    mulss  xmm0, xmm0 ;# rinvsq
    movaps xmm2, xmm0 ;# rinvsq
    mulss  xmm2, xmm2 ;# rinv4
    mulss  xmm2, xmm0 ;# rinv6
    movaps xmm12, xmm2 
    mulss  xmm12, xmm12 ;# rinv12

    ;# load table data
   	movss xmm4, [rsi + r12*4]
	movss xmm5, [rsi + r12*4 + 4]
   	movss xmm6, [rsi + r12*4 + 8]
	movss xmm7, [rsi + r12*4 + 12]
    ;# table data ready in xmm4-xmm7

    mulss  xmm2, [rsp + nb310_c6]    ;# vvdw6=c6*rinv6
	mulss  xmm12, [rsp + nb310_c12]   ;# vvdw12=c12*rinv12     

	movaps xmm0, xmm12
	subss  xmm12, xmm2	;# Vvdw=Vvdw12-Vvdw6

    ;# add potential to vvdwtot 
	addss  xmm12, [rsp + nb310_Vvdwtot]
    movss [rsp + nb310_Vvdwtot], xmm12
    
    mulss xmm7, xmm3   ;# Heps
    mulss  xmm6, xmm3  ;# Geps
    mulss xmm7, xmm3   ;# Heps2

    addss  xmm5, xmm6   ;# F+Geps
    addss  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addss  xmm7, xmm7   ;# 2*Heps2
    addss  xmm7, xmm6   ;# 2*Heps2+Geps
    addss  xmm7, xmm5   ;# FF = Fp + 2*Heps2 + Geps
    mulss  xmm5, xmm3   ;# eps*Fp
    addss  xmm5, xmm4   ;# VV
    mulss  xmm5, [rsp + nb310_qq]   ;# VV*qq=vcoul
    mulss  xmm7, [rsp + nb310_qq]   ;# FF*qq=fijC

    ;# LJ forces
    mulss  xmm2, [rsp + nb310_six]
    mulss  xmm0, [rsp + nb310_twelve]
    subss  xmm0, xmm2
    mulss  xmm0, xmm1 ;# (12*vnb12-6*vnb6)*rinv

    ;# add potential to vctot 
	addss  xmm5, [rsp + nb310_vctot]
    movss [rsp + nb310_vctot], xmm5
    
    mulss  xmm7, [rsp + nb310_tsc]
    subss  xmm0, xmm7
    
    mulss  xmm0, xmm1  ;# fscal

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulss  xmm9, xmm0
	mulss  xmm10, xmm0
	mulss  xmm11, xmm0
    
	;# accumulate i forces
    addss xmm13, xmm9
    addss xmm14, xmm10
    addss xmm15, xmm11

	mov rsi, [rbp + nb310_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + rax*4]
    addss  xmm10, [rsi + rax*4 + 4]
    addss  xmm11, [rsi + rax*4 + 8]
    movss  [rsi + rax*4],     xmm9
    movss  [rsi + rax*4 + 4], xmm10
    movss  [rsi + rax*4 + 8], xmm11
    
.nb310_updateouterdata:
	mov   ecx, [rsp + nb310_ii3]
	mov   rdi, [rbp + nb310_faction]
	mov   rsi, [rbp + nb310_fshift]
	mov   edx, [rsp + nb310_is3]

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
	mov esi, [rsp + nb310_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb310_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb310_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb310_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb310_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb310_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb310_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb310_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb310_n], esi
        jmp .nb310_outer
.nb310_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb310_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb310_end
        ;# non-zero, do one more workunit
        jmp   .nb310_threadloop
.nb310_end:

	mov eax, [rsp + nb310_nouter]
	mov ebx, [rsp + nb310_ninner]
	mov rcx, [rbp + nb310_outeriter]
	mov rdx, [rbp + nb310_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 464
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret





.globl nb_kernel310nf_x86_64_sse
.globl _nb_kernel310nf_x86_64_sse
nb_kernel310nf_x86_64_sse:	
_nb_kernel310nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb310nf_fshift,         16
.equiv          nb310nf_gid,            24
.equiv          nb310nf_pos,            32
.equiv          nb310nf_faction,        40
.equiv          nb310nf_charge,         48
.equiv          nb310nf_p_facel,        56
.equiv          nb310nf_argkrf,         64
.equiv          nb310nf_argcrf,         72
.equiv          nb310nf_Vc,             80
.equiv          nb310nf_type,           88
.equiv          nb310nf_p_ntype,        96
.equiv          nb310nf_vdwparam,       104
.equiv          nb310nf_Vvdw,           112
.equiv          nb310nf_p_tabscale,     120
.equiv          nb310nf_VFtab,          128
.equiv          nb310nf_invsqrta,       136
.equiv          nb310nf_dvda,           144
.equiv          nb310nf_p_gbtabscale,   152
.equiv          nb310nf_GBtab,          160
.equiv          nb310nf_p_nthreads,     168
.equiv          nb310nf_count,          176
.equiv          nb310nf_mtx,            184
.equiv          nb310nf_outeriter,      192
.equiv          nb310nf_inneriter,      200
.equiv          nb310nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb310nf_ix,             0
.equiv          nb310nf_iy,             16
.equiv          nb310nf_iz,             32
.equiv          nb310nf_iq,             48
.equiv          nb310nf_tsc,            64
.equiv          nb310nf_qq,             80
.equiv          nb310nf_c6,             96
.equiv          nb310nf_c12,            112
.equiv          nb310nf_vctot,          128
.equiv          nb310nf_Vvdwtot,        144
.equiv          nb310nf_half,           160
.equiv          nb310nf_three,          176
.equiv          nb310nf_nri,            192
.equiv          nb310nf_iinr,           200
.equiv          nb310nf_jindex,         208
.equiv          nb310nf_jjnr,           216
.equiv          nb310nf_shift,          224
.equiv          nb310nf_shiftvec,       232
.equiv          nb310nf_facel,          240
.equiv          nb310nf_innerjjnr,      248
.equiv          nb310nf_is3,            256
.equiv          nb310nf_ii3,            260
.equiv          nb310nf_ntia,           264
.equiv          nb310nf_innerk,         268
.equiv          nb310nf_n,              272
.equiv          nb310nf_nn1,            276
.equiv          nb310nf_ntype,          280
.equiv          nb310nf_nouter,         284
.equiv          nb310nf_ninner,         288

	push rbp
	mov  rbp, rsp
    
    ;# Push integer registers on stack
	push rbx
    push rsi
    push rdi
    push r12
    push r13
    push r14
    push r15

    ;# Make room for registers xmm6-xmm15 (10 registers=160 bytes)
    sub rsp, 168
    
    ;# Save xmm registers to stack
    movaps [rsp      ], xmm6
    movaps [rsp + 16 ], xmm7
    movaps [rsp + 32 ], xmm8
    movaps [rsp + 48 ], xmm9
    movaps [rsp + 64 ], xmm10
    movaps [rsp + 80 ], xmm11
    movaps [rsp + 96 ], xmm12
    movaps [rsp + 112], xmm13
    movaps [rsp + 128], xmm14
    movaps [rsp + 144], xmm15
    
	emms
	sub rsp, 304		;# local variable stack space (n*16+8)
; .if 0    # block below only read by NASM - special calling convention on win64
%ifidn __OUTPUT_FORMAT__, win64
    ;# Adjust rbp to account for shadow space (32) & two extra args (2*8) on stack
    add rbp, 48
    ;# Adjust stack pointer for different alignment
    ;# Move around arguments to fit AMD64 convention below
    ;# AMD64 passes args in: rdi,rsi,rdx,rcx,r8,r9 + stack
    ;# win64 passes args in: rcx,rdx,r8,r9         + stack
    mov rdi, rcx
    mov rsi, rdx
    mov rdx, r8
    mov rcx, r9
    mov r8,  [rbp]
    mov r9,  [rbp + 8]
%endif
; .endif   # end NASM- and win64-specific block

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb310nf_nouter], eax
	mov [rsp + nb310nf_ninner], eax
	
	mov edi, [rdi]
	mov [rsp + nb310nf_nri], edi
	mov [rsp + nb310nf_iinr], rsi
	mov [rsp + nb310nf_jindex], rdx
	mov [rsp + nb310nf_jjnr], rcx
	mov [rsp + nb310nf_shift], r8
	mov [rsp + nb310nf_shiftvec], r9
	mov rdi, [rbp + nb310nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb310nf_ntype], edi
	mov rsi, [rbp + nb310nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb310nf_facel], xmm0

	mov rax, [rbp + nb310nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb310nf_tsc],  xmm3	

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb310nf_half], eax
	movss xmm1, [rsp + nb310nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb310nf_half],  xmm1
	movaps [rsp + nb310nf_three],  xmm3

.nb310nf_threadloop:
        mov   rsi, [rbp + nb310nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb310nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb310nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb310nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb310nf_n], eax
        mov [rsp + nb310nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb310nf_outerstart
        jmp .nb310nf_end

.nb310nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb310nf_nouter]
	mov [rsp + nb310nf_nouter], ebx

.nb310nf_outer:
	mov   rax, [rsp + nb310nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb310nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb310nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb310nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb310nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb310nf_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb310nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb310nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb310nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb310nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb310nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb310nf_ix], xmm0
	movaps [rsp + nb310nf_iy], xmm1
	movaps [rsp + nb310nf_iz], xmm2

	mov   [rsp + nb310nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb310nf_vctot], xmm4
	movaps [rsp + nb310nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb310nf_jindex]
	mov   rcx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb310nf_pos]
	mov   rax, [rsp + nb310nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb310nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb310nf_ninner]
	mov   [rsp + nb310nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb310nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb310nf_unroll_loop
	jmp   .nb310nf_finish_inner
.nb310nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb310nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb310nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb310nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm2, [rsp + nb310nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	mulps  xmm3, xmm2
	movd  mm2, ecx
	movd  mm3, edx

	movaps [rsp + nb310nf_qq], xmm3
	
	mov rsi, [rbp + nb310nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb310nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb310nf_ntia]
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

	movaps [rsp + nb310nf_c6], xmm4
	movaps [rsp + nb310nf_c12], xmm6
	
	mov rsi, [rbp + nb310nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [rsp + nb310nf_ix]
	movaps xmm5, [rsp + nb310nf_iy]
	movaps xmm6, [rsp + nb310nf_iz]

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
	movaps xmm1, [rsp + nb310nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb310nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb310nf_tsc]

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

	mov  rsi, [rbp + nb310nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# coulomb table ready, in xmm4-xmm7  	
	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	movaps xmm3, [rsp + nb310nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# L-J 
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	;# xmm4=rinvsq 

	;# at this point mm5 contains vcoul  
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addps  xmm5, [rsp + nb310nf_vctot]
	movaps xmm6, xmm4
 	mulps  xmm6, xmm4
	movaps [rsp + nb310nf_vctot], xmm5 

	mulps  xmm6, xmm4	;# xmm6=rinvsix 
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulps  xmm6, [rsp + nb310nf_c6]
	mulps  xmm4, [rsp + nb310nf_c12]
	movaps xmm7, [rsp + nb310nf_Vvdwtot]
	addps  xmm7, xmm4
	subps  xmm7, xmm6
	movaps [rsp + nb310nf_Vvdwtot], xmm7


	;# should we do one more iteration? 
	sub dword ptr [rsp + nb310nf_innerk],  4
	jl    .nb310nf_finish_inner
	jmp   .nb310nf_unroll_loop
.nb310nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb310nf_innerk],  4
	mov   edx, [rsp + nb310nf_innerk]
	and   edx, 2
	jnz   .nb310nf_dopair
	jmp   .nb310nf_checksingle
.nb310nf_dopair:	
	mov rsi, [rbp + nb310nf_charge]
    mov   rcx, [rsp + nb310nf_innerjjnr]
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb310nf_innerjjnr],  8	
	xorps xmm7, xmm7
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 0 
	shufps xmm3, xmm3, 8 ;# 00001000 ;# xmm3(0,1) has the charges 

	mulps  xmm3, [rsp + nb310nf_iq]
	movlhps xmm3, xmm7
	movaps [rsp + nb310nf_qq], xmm3

	mov rsi, [rbp + nb310nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]	
	mov rsi, [rbp + nb310nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb310nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [rsi + rcx*4]
	movhps xmm6, [rsi + rdx*4]
	mov rdi, [rbp + nb310nf_pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# 00001000 	
	shufps xmm6, xmm6, 13 ;# 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [rsp + nb310nf_c6], xmm4
	movaps [rsp + nb310nf_c12], xmm6	
	
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
	
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb310nf_ix]
	movaps xmm5, [rsp + nb310nf_iy]
	movaps xmm6, [rsp + nb310nf_iz]

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
	movaps xmm1, [rsp + nb310nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb310nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb310nf_tsc]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb310nf_VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6

	movlps xmm5, [rsi + rcx*4]
	movhps xmm5, [rsi + rdx*4] ;# got half coulomb table 
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
	
	movlps xmm7, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rdx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 136  ;# 10001000
	shufps xmm7, xmm7, 221  ;# 11011101
	;# table ready in xmm4-xmm7 

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	movaps xmm3, [rsp + nb310nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# L-J 
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	;# xmm4=rinvsq 

	;# at this point mm5 contains vcoul  
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addps  xmm5, [rsp + nb310nf_vctot]

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	movaps [rsp + nb310nf_vctot], xmm5 

	mulps  xmm6, xmm4	;# xmm6=rinvsix 
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulps  xmm6, [rsp + nb310nf_c6]
	mulps  xmm4, [rsp + nb310nf_c12]
	movaps xmm7, [rsp + nb310nf_Vvdwtot]
	addps  xmm7, xmm4
	subps  xmm7, xmm6
	movaps [rsp + nb310nf_Vvdwtot], xmm7

.nb310nf_checksingle:				
	mov   edx, [rsp + nb310nf_innerk]
	and   edx, 1
	jnz    .nb310nf_dosingle
	jmp    .nb310nf_updateouterdata
.nb310nf_dosingle:
	mov rsi, [rbp + nb310nf_charge]
	mov rdi, [rbp + nb310nf_pos]
	mov   rcx, [rsp + nb310nf_innerjjnr]
	mov   eax, [rcx]	
	xorps  xmm6, xmm6
	movss xmm6, [rsi + rax*4]	;# xmm6(0) has the charge 	
	mulps  xmm6, [rsp + nb310nf_iq]
	movaps [rsp + nb310nf_qq], xmm6

	mov rsi, [rbp + nb310nf_type]
	mov ecx, eax
	mov ecx, [rsi + rcx*4]	
	mov rsi, [rbp + nb310nf_vdwparam]
	shl ecx, 1
	add ecx, [rsp + nb310nf_ntia]
	movlps xmm6, [rsi + rcx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# 11111100	
	shufps xmm6, xmm6, 253  ;# 11111101	
	
	movaps [rsp + nb310nf_c6], xmm4
	movaps [rsp + nb310nf_c12], xmm6	
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	 
	
	movaps xmm4, [rsp + nb310nf_ix]
	movaps xmm5, [rsp + nb310nf_iy]
	movaps xmm6, [rsp + nb310nf_iz]

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
	movaps xmm1, [rsp + nb310nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb310nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb310nf_tsc]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb310nf_VFtab]
	movd ebx, mm6
	
	movlps xmm4, [rsi + rbx*4]
	movlps xmm6, [rsi + rbx*4 + 8]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1
	shufps xmm7, xmm7, 1
	;# table ready in xmm4-xmm7 

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	movaps xmm3, [rsp + nb310nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# L-J 
	movaps xmm4, xmm0
	mulps  xmm4, xmm0	;# xmm4=rinvsq 

	;# at this point mm5 contains vcoul  
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addss  xmm5, [rsp + nb310nf_vctot]

	movaps xmm6, xmm4
	mulps  xmm6, xmm4

	movss [rsp + nb310nf_vctot], xmm5 

	mulps  xmm6, xmm4	;# xmm6=rinvsix 
	movaps xmm4, xmm6
	mulps  xmm4, xmm4	;# xmm4=rinvtwelve 
	mulps  xmm6, [rsp + nb310nf_c6]
	mulps  xmm4, [rsp + nb310nf_c12]
	movss xmm7, [rsp + nb310nf_Vvdwtot]
	addps  xmm7, xmm4
	subps  xmm7, xmm6
	movss [rsp + nb310nf_Vvdwtot], xmm7

.nb310nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb310nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb310nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb310nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb310nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb310nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb310nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb310nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb310nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb310nf_n], esi
        jmp .nb310nf_outer
.nb310nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb310nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb310nf_end
        ;# non-zero, do one more workunit
        jmp   .nb310nf_threadloop
.nb310nf_end:

	mov eax, [rsp + nb310nf_nouter]
	mov ebx, [rsp + nb310nf_ninner]
	mov rcx, [rbp + nb310nf_outeriter]
	mov rdx, [rbp + nb310nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 304
	emms

    ;# Save xmm registers to stack
    movaps xmm6,  [rsp      ]
    movaps xmm7,  [rsp + 16 ]
    movaps xmm8,  [rsp + 32 ]
    movaps xmm9,  [rsp + 48 ]
    movaps xmm10, [rsp + 64 ]
    movaps xmm11, [rsp + 80 ]
    movaps xmm12, [rsp + 96 ]
    movaps xmm13, [rsp + 112]
    movaps xmm14, [rsp + 128]
    movaps xmm15, [rsp + 144]

    ;# Reset pointers after restoring xmm6-15
    add rsp, 168

    pop r15
    pop r14
    pop r13
    pop r12
    pop rdi
    pop rsi
    pop rbx
    
	pop	rbp
	ret
