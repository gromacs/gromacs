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

	

.globl nb_kernel300_x86_64_sse
.globl _nb_kernel300_x86_64_sse
nb_kernel300_x86_64_sse:	
_nb_kernel300_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb300_fshift,           16
.equiv          nb300_gid,              24
.equiv          nb300_pos,              32
.equiv          nb300_faction,          40
.equiv          nb300_charge,           48
.equiv          nb300_p_facel,          56
.equiv          nb300_argkrf,           64
.equiv          nb300_argcrf,           72
.equiv          nb300_Vc,               80
.equiv          nb300_type,             88
.equiv          nb300_p_ntype,          96
.equiv          nb300_vdwparam,         104
.equiv          nb300_Vvdw,             112
.equiv          nb300_p_tabscale,       120
.equiv          nb300_VFtab,            128
.equiv          nb300_invsqrta,         136
.equiv          nb300_dvda,             144
.equiv          nb300_p_gbtabscale,     152
.equiv          nb300_GBtab,            160
.equiv          nb300_p_nthreads,       168
.equiv          nb300_count,            176
.equiv          nb300_mtx,              184
.equiv          nb300_outeriter,        192
.equiv          nb300_inneriter,        200
.equiv          nb300_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
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
.equiv          nb300_innerjjnr,        272
.equiv          nb300_nri,              280
.equiv          nb300_iinr,             288
.equiv          nb300_jindex,           296
.equiv          nb300_jjnr,             304
.equiv          nb300_shift,            312
.equiv          nb300_shiftvec,         320
.equiv          nb300_facel,            328
.equiv          nb300_innerk,           336
.equiv          nb300_is3,              344
.equiv          nb300_ii3,              348
.equiv          nb300_n,                352
.equiv          nb300_nn1,              356
.equiv          nb300_nouter,           360
.equiv          nb300_ninner,           364

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

	sub rsp, 368		;# local variable stack space (n*16+8)
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
	mov [rsp + nb300_nouter], eax
	mov [rsp + nb300_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb300_nri], edi
	mov [rsp + nb300_iinr], rsi
	mov [rsp + nb300_jindex], rdx
	mov [rsp + nb300_jjnr], rcx
	mov [rsp + nb300_shift], r8
	mov [rsp + nb300_shiftvec], r9
	mov rsi, [rbp + nb300_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb300_facel], xmm0

	mov rax, [rbp + nb300_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb300_tsc], xmm3

	mov r8,  [rbp + nb300_pos]    
	mov r9,  [rbp + nb300_faction] 
	mov r10, [rbp + nb300_charge] 
	mov r11, [rbp + nb300_VFtab] 

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb300_half], eax
	movss xmm1, [rsp + nb300_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb300_half],  xmm1
	movaps [rsp + nb300_two],  xmm2
	movaps [rsp + nb300_three],  xmm3

.nb300_threadloop:
        mov   rsi, [rbp + nb300_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb300_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb300_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb300_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb300_n], eax
        mov [rsp + nb300_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb300_outerstart
        jmp .nb300_end

.nb300_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb300_nouter]
	mov [rsp + nb300_nouter], ebx

.nb300_outer:
	mov   rax, [rsp + nb300_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb300_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb300_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb300_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb300_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb300_facel]
	shufps xmm3, xmm3, 0

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb300_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb300_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb300_ix], xmm0
	movaps [rsp + nb300_iy], xmm1
	movaps [rsp + nb300_iz], xmm2

	mov   [rsp + nb300_ii3], ebx
	
	;# clear vctot (xmm12) and i forces (xmm13-xmm15)
	xorps xmm12, xmm12
	movaps xmm13, xmm12
	movaps xmm14, xmm12
	movaps xmm15, xmm12
	
	mov   rax, [rsp + nb300_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb300_pos]
	mov   rdi, [rbp + nb300_faction]	
	mov   rax, [rsp + nb300_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb300_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb300_ninner]
	mov   [rsp + nb300_ninner], ecx
	add   edx, 0
	mov   [rsp + nb300_innerk], edx    ;# number of innerloop atoms 
	jge   .nb300_unroll_loop
	jmp   .nb300_finish_inner
.nb300_unroll_loop:	
	;# quad-unrolled innerloop here 
	mov   rdx, [rsp + nb300_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	mov   r10d, [rdx + 8]            
	mov   r11d, [rdx + 12]         ;# eax-edx=jnr1-4 
	mov rdi, [rbp + nb300_pos]
    
	add qword ptr [rsp + nb300_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	lea   rax, [r8 + r8*2]     ;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	
	lea   rcx, [r10 + r10*2]     ;# replace jnr with j3 
	lea   rdx, [r11 + r11*2]	

	;# load coordinates
    
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

    movaps xmm3, xmm1

    unpcklps xmm1, xmm2 ;# x1 x2 x3 x4
    unpckhps xmm3, xmm2 ;# y1 y2 y3 y4
    unpcklps xmm5, xmm6 ;# z1 z2 z3 z4

	;# calc dr  
	subps xmm1, [rsp + nb300_ix]
	subps xmm3, [rsp + nb300_iy]
	subps xmm5, [rsp + nb300_iz]

	;# store dr in xmm9-xmm11
    movaps xmm9, xmm1
    movaps xmm10, xmm3
    movaps xmm11, xmm5
    
	;# square it 
	mulps xmm1, xmm1
	mulps xmm3, xmm3
	mulps xmm5, xmm5
	addps xmm3, xmm1
	addps xmm3, xmm5
	;# rsq in xmm3

    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm3
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb300_three]
	mulps xmm5, xmm3	;# rsq*lu*lu 	
    subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm1, [rsp + nb300_half]	
    ;# xmm1=rinv
    ;# xmm3=rsq
    
    mulps xmm3, xmm1 ;# r
    mulps xmm3, [rsp + nb300_tsc] ;# rtab
    
    ;# truncate and convert to integers
    cvttps2dq xmm2, xmm3
    
    ;# convert back to float
    cvtdq2ps  xmm0, xmm2
    
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

	mov rsi, [rbp + nb300_VFtab]
    ;# load table data
   	movlps xmm5, [rsi + r12*4]
	movlps xmm7, [rsi + r14*4]
	movhps xmm5, [rsi + r13*4]
	movhps xmm7, [rsi + r15*4]

    movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
    
	movlps xmm7, [rsi + r12*4 + 8]   
	movlps xmm8, [rsi + r14*4 + 8]
	movhps xmm7, [rsi + r13*4 + 8]
	movhps xmm8, [rsi + r15*4 + 8]

    movaps xmm6, xmm7
	mov rdi, [rbp + nb300_charge]
    
	shufps xmm6, xmm8, 136  ;# 10001000
	shufps xmm7, xmm8, 221  ;# 11011101
    ;# table data ready in xmm4-xmm7
    
	movss xmm0, [rdi + r8*4]
	movss xmm2, [rdi + r10*4]

    mulps xmm7, xmm3   ;# Heps
    mulps  xmm6, xmm3  ;# Geps
    mulps xmm7, xmm3   ;# Heps2

    unpcklps xmm0, xmm2
	movss xmm2, [rdi + r9*4]
	movss xmm8, [rdi + r11*4]

    addps  xmm5, xmm6   ;# F+Geps
    addps  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addps  xmm7, xmm7   ;# 2*Heps2
    unpcklps xmm2, xmm8
    unpcklps xmm0, xmm2
    addps  xmm7, xmm6   ;# 2*Heps2+Geps
    addps  xmm7, xmm5   ;# FF = Fp + 2*Heps2 + Geps
    mulps  xmm5, xmm3   ;# eps*Fp
    mulps  xmm0, [rsp + nb300_iq]
    addps  xmm5, xmm4   ;# VV
    mulps  xmm5, xmm0   ;# VV*qq=vcoul
    mulps  xmm7, xmm0   ;# FF*qq=fijC

    ;# add potential to vctot (sum in xmm12)
	addps  xmm12, xmm5

    mulps  xmm7, [rsp + nb300_tsc]
    mulps  xmm7, xmm1  
    
    xorps  xmm4, xmm4
    subps  xmm4, xmm7   ;# fscal

	mov rsi, [rbp + nb300_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movlps xmm1, [rsi + rcx*4] ;# x3 y3 - -
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + rdx*4] ;# x3 y3 x4 y4

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm4
	mulps  xmm10, xmm4
	mulps  xmm11, xmm4

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
	sub dword ptr [rsp + nb300_innerk],  4
	jl    .nb300_finish_inner
	jmp   .nb300_unroll_loop
.nb300_finish_inner:
    ;# check if at least two particles remain 
    add dword ptr [rsp + nb300_innerk],  4
    mov   edx, [rsp + nb300_innerk]
    and   edx, 2
    jnz   .nb300_dopair
    jmp   .nb300_checksingle
.nb300_dopair:  
	;# twice-unrolled innerloop here 
	mov   rdx, [rsp + nb300_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
    
	add qword ptr [rsp + nb300_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb300_charge]
	movss xmm8, [rsi + rax*4]
	movss xmm1, [rsi + rbx*4]

    unpcklps xmm8, xmm1  ;# jqa jqb - -
	mulps xmm8, [rsp + nb300_iq]    ;#qq

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# load coordinates
	mov rdi, [rbp + nb300_pos]
    
	movlps xmm4, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm5, [rdi + rbx*4]	;# x2 y2 - - 

	movss xmm6, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm7, [rdi + rbx*4 + 8]	;# z2 - - - 

    unpcklps xmm4, xmm5 ;# x1 x2 y1 y2
    movhlps  xmm5, xmm4 ;# y1 y2 -  -
    unpcklps xmm6, xmm7 ;# z1 z2 -  -

	;# calc dr  
	subps xmm4, [rsp + nb300_ix]
	subps xmm5, [rsp + nb300_iy]
	subps xmm6, [rsp + nb300_iz]

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
	movaps xmm1, [rsp + nb300_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 	
    subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm1, [rsp + nb300_half]		
    ;# xmm1=rinv
    movaps xmm3, xmm4
    ;# xmm3=rsq 

    mulps xmm3, xmm1 ;# r
    mulps xmm3, [rsp + nb300_tsc] ;# rtab
    
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

	mov rsi, [rbp + nb300_VFtab]
    ;# load table data
   	movlps xmm4, [rsi + r12*4]
	movlps xmm5, [rsi + r13*4]
    unpcklps xmm4, xmm5
    movhlps  xmm5, xmm4
    
   	movlps xmm6, [rsi + r12*4 + 8]
	movlps xmm7, [rsi + r13*4 + 8] 
    unpcklps xmm6, xmm7
    movhlps  xmm7, xmm6
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
    mulps  xmm5, xmm8   ;# VV*qq=vcoul
    mulps  xmm7, xmm8   ;# FF*qq=fijC

    xorps xmm6, xmm6
    movlhps xmm5, xmm6
    
    ;# add potential to vctot (sum in xmm12)
	addps  xmm12, xmm5

    mulps  xmm7, [rsp + nb300_tsc]
    mulps  xmm7, xmm1  
    
    xorps  xmm4, xmm4
    subps  xmm4, xmm7   ;# fscal

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulps  xmm9, xmm4
	mulps  xmm10, xmm4
	mulps  xmm11, xmm4

    movlhps xmm9, xmm6
    movlhps xmm10, xmm6
    movlhps xmm11, xmm6
    
	;# accumulate i forces
    addps xmm13, xmm9
    addps xmm14, xmm10
    addps xmm15, xmm11

	mov rsi, [rbp + nb300_faction]
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

.nb300_checksingle:                             
    mov   edx, [rsp + nb300_innerk]
    and   edx, 1
    jnz    .nb300_dosingle
    jmp    .nb300_updateouterdata

.nb300_dosingle:	
    mov rcx, [rsp + nb300_innerjjnr]
	mov   eax, [rcx]	            

	mov rsi, [rbp + nb300_charge]
	movss xmm8, [rsi + rax*4]       ;# jq
	mulss xmm8, [rsp + nb300_iq]    ;# qq

	lea   rax, [rax + rax*2]        ;# replace jnr with j3 

	mov rdi, [rbp + nb300_pos]
	movss xmm4, [rdi + rax*4]	    ;# x1 - - - 
	movss xmm5, [rdi + rax*4 + 4]    ;# y2 - - - 
	movss xmm6, [rdi + rax*4 + 8]    ;# 13 - - - 

	;# calc dr  
	subss xmm4, [rsp + nb300_ix]
	subss xmm5, [rsp + nb300_iy]
	subss xmm6, [rsp + nb300_iz]

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
	movaps xmm1, [rsp + nb300_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 	
    subss xmm1, xmm5	;# 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm1, [rsp + nb300_half]		
    ;# xmm1=rinv
    movaps xmm3, xmm4
    ;# xmm3=rsq 

    mulss xmm3, xmm1 ;# r
    mulss xmm3, [rsp + nb300_tsc] ;# rtab
    
    ;# truncate and convert to integers
    cvttss2si r12d, xmm3
    
    ;# convert back to float
    cvtsi2ss  xmm0, r12d
    
    ;# multiply by 4
    shl       r12d, 2

    ;# calculate eps
    subss     xmm3, xmm0

	mov rsi, [rbp + nb300_VFtab]
    ;# load table data
   	movss  xmm4, [rsi + r12*4]
   	movss  xmm5, [rsi + r12*4 + 4]
   	movss  xmm6, [rsi + r12*4 + 8]
   	movss  xmm7, [rsi + r12*4 + 12]
    ;# table data ready in xmm4-xmm7
    
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
    mulss  xmm5, xmm8   ;# VV*qq=vcoul
    mulss  xmm7, xmm8   ;# FF*qq=fijC

    ;# add potential to vctot (sum in xmm12)
	addss  xmm12, xmm5

    mulss  xmm7, [rsp + nb300_tsc]
    mulss  xmm7, xmm1  
    
    xorps  xmm4, xmm4
    subss  xmm4, xmm7   ;# fscal

    ;# calculate scalar force by multiplying dx/dy/dz with fscal
	mulss  xmm9, xmm4
	mulss  xmm10, xmm4
	mulss  xmm11, xmm4

	;# accumulate i forces
    addss xmm13, xmm9
    addss xmm14, xmm10
    addss xmm15, xmm11

	mov rsi, [rbp + nb300_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + rax*4]
    addss  xmm10, [rsi + rax*4 + 4]
    addss  xmm11, [rsi + rax*4 + 8]
    movss  [rsi + rax*4],     xmm9
    movss  [rsi + rax*4 + 4], xmm10
    movss  [rsi + rax*4 + 8], xmm11
    
.nb300_updateouterdata:
	mov   ecx, [rsp + nb300_ii3]
	mov   rdi, [rbp + nb300_faction]
	mov   rsi, [rbp + nb300_fshift]
	mov   edx, [rsp + nb300_is3]

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
	mov esi, [rsp + nb300_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb300_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	;# accumulate 
	movhlps xmm6, xmm12
	addps  xmm12, xmm6	;# pos 0-1 in xmm12 have the sum now 
	movaps xmm6, xmm12
	shufps xmm6, xmm6, 1
	addss  xmm12, xmm6

	;# add earlier value from mem 
	mov   rax, [rbp + nb300_Vc]
	addss xmm12, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm12
	
        ;# finish if last 
        mov ecx, [rsp + nb300_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb300_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb300_n], esi
        jmp .nb300_outer
.nb300_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb300_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb300_end
        ;# non-zero, do one more workunit
        jmp   .nb300_threadloop
.nb300_end:
	mov eax, [rsp + nb300_nouter]
	mov ebx, [rsp + nb300_ninner]
	mov rcx, [rbp + nb300_outeriter]
	mov rdx, [rbp + nb300_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 368
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




.globl nb_kernel300nf_x86_64_sse
.globl _nb_kernel300nf_x86_64_sse
nb_kernel300nf_x86_64_sse:	
_nb_kernel300nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb300nf_fshift,         16
.equiv          nb300nf_gid,            24
.equiv          nb300nf_pos,            32
.equiv          nb300nf_faction,        40
.equiv          nb300nf_charge,         48
.equiv          nb300nf_p_facel,        56
.equiv          nb300nf_argkrf,         64
.equiv          nb300nf_argcrf,         72
.equiv          nb300nf_Vc,             80
.equiv          nb300nf_type,           88
.equiv          nb300nf_p_ntype,        96
.equiv          nb300nf_vdwparam,       104
.equiv          nb300nf_Vvdw,           112
.equiv          nb300nf_p_tabscale,     120
.equiv          nb300nf_VFtab,          128
.equiv          nb300nf_invsqrta,       136
.equiv          nb300nf_dvda,           144
.equiv          nb300nf_p_gbtabscale,   152
.equiv          nb300nf_GBtab,          160
.equiv          nb300nf_p_nthreads,     168
.equiv          nb300nf_count,          176
.equiv          nb300nf_mtx,            184
.equiv          nb300nf_outeriter,      192
.equiv          nb300nf_inneriter,      200
.equiv          nb300nf_work,           208
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
.equiv          nb300nf_nri,            160
.equiv          nb300nf_iinr,           168
.equiv          nb300nf_jindex,         176
.equiv          nb300nf_jjnr,           184
.equiv          nb300nf_shift,          192
.equiv          nb300nf_shiftvec,       200
.equiv          nb300nf_facel,          208
.equiv          nb300nf_innerk,         216
.equiv          nb300nf_n,              220
.equiv          nb300nf_nn1,            224
.equiv          nb300nf_nouter,         228
.equiv          nb300nf_ninner,         232

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
	sub rsp, 240		;# local variable stack space (n*16+8)
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
	mov [rsp + nb300nf_nouter], eax
	mov [rsp + nb300nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb300nf_nri], edi
	mov [rsp + nb300nf_iinr], rsi
	mov [rsp + nb300nf_jindex], rdx
	mov [rsp + nb300nf_jjnr], rcx
	mov [rsp + nb300nf_shift], r8
	mov [rsp + nb300nf_shiftvec], r9
	mov rsi, [rbp + nb300nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb300nf_facel], xmm0

	mov rax, [rbp + nb300nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb300nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb300nf_half], eax
	movss xmm1, [rsp + nb300nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb300nf_half],  xmm1
	movaps [rsp + nb300nf_three],  xmm3

.nb300nf_threadloop:
        mov   rsi, [rbp + nb300nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb300nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb300nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb300nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb300nf_n], eax
        mov [rsp + nb300nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb300nf_outerstart
        jmp .nb300nf_end
.nb300nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb300nf_nouter]
	mov [rsp + nb300nf_nouter], ebx

.nb300nf_outer:
	mov   rax, [rsp + nb300nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb300nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb300nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb300nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb300nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb300nf_facel]
	shufps xmm3, xmm3, 0

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb300nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb300nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb300nf_ix], xmm0
	movaps [rsp + nb300nf_iy], xmm1
	movaps [rsp + nb300nf_iz], xmm2

	mov   [rsp + nb300nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb300nf_vctot], xmm4
	
	mov   rax, [rsp + nb300nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb300nf_pos]
	mov   rax, [rsp + nb300nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb300nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb300nf_ninner]
	mov   [rsp + nb300nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb300nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb300nf_unroll_loop
	jmp   .nb300nf_finish_inner
.nb300nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb300nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb300nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb300nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm2, [rsp + nb300nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	mulps  xmm3, xmm2

	movaps [rsp + nb300nf_qq], xmm3	
	
	mov rsi, [rbp + nb300nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [rsp + nb300nf_ix]
	movaps xmm5, [rsp + nb300nf_iy]
	movaps xmm6, [rsp + nb300nf_iz]

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
	movaps xmm1, [rsp + nb300nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb300nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb300nf_tsc]

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

	mov  rsi, [rbp + nb300nf_VFtab]
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
	movaps xmm3, [rsp + nb300nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  

	;# at this point xmm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addps  xmm5, [rsp + nb300nf_vctot]
	movaps [rsp + nb300nf_vctot], xmm5 

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb300nf_innerk],  4
	jl    .nb300nf_finish_inner
	jmp   .nb300nf_unroll_loop
.nb300nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb300nf_innerk],  4
	mov   edx, [rsp + nb300nf_innerk]
	and   edx, 2
	jnz   .nb300nf_dopair
	jmp   .nb300nf_checksingle
.nb300nf_dopair:	
	mov rsi, [rbp + nb300nf_charge]

    mov   rcx, [rsp + nb300nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb300nf_innerjjnr],  8	
	xorps xmm7, xmm7
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 0 
	shufps xmm3, xmm3, 8 ;# 00001000 ;# xmm3(0,1) has the charges 

	mulps  xmm3, [rsp + nb300nf_iq]
	movlhps xmm3, xmm7
	movaps [rsp + nb300nf_qq], xmm3

	mov rdi, [rbp + nb300nf_pos]	
	
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
	
	movaps xmm4, [rsp + nb300nf_ix]
	movaps xmm5, [rsp + nb300nf_iy]
	movaps xmm6, [rsp + nb300nf_iz]

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
	movaps xmm1, [rsp + nb300nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb300nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb300nf_tsc]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb300nf_VFtab]
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
	movaps xmm3, [rsp + nb300nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addps  xmm5, [rsp + nb300nf_vctot]
	movaps [rsp + nb300nf_vctot], xmm5 

.nb300nf_checksingle:				
	mov   edx, [rsp + nb300nf_innerk]
	and   edx, 1
	jnz    .nb300nf_dosingle
	jmp    .nb300nf_updateouterdata
.nb300nf_dosingle:
	mov rsi, [rbp + nb300nf_charge]
	mov rdi, [rbp + nb300nf_pos]
	mov   rcx, [rsp + nb300nf_innerjjnr]
	mov   eax, [rcx]	
	xorps  xmm6, xmm6
	movss xmm6, [rsi + rax*4]	;# xmm6(0) has the charge 	
	mulps  xmm6, [rsp + nb300nf_iq]
	movaps [rsp + nb300nf_qq], xmm6
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	 
	
	movaps xmm4, [rsp + nb300nf_ix]
	movaps xmm5, [rsp + nb300nf_iy]
	movaps xmm6, [rsp + nb300nf_iz]

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
	movaps xmm1, [rsp + nb300nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb300nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb300nf_tsc]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb300nf_VFtab]
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
	movaps xmm3, [rsp + nb300nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV 
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addss  xmm5, [rsp + nb300nf_vctot]
	movss [rsp + nb300nf_vctot], xmm5 

.nb300nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb300nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb300nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb300nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb300nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb300nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb300nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb300nf_n], esi
        jmp .nb300nf_outer
.nb300nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb300nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb300nf_end
        ;# non-zero, do one more workunit
        jmp   .nb300nf_threadloop
.nb300nf_end:

	mov eax, [rsp + nb300nf_nouter]
	mov ebx, [rsp + nb300nf_ninner]
	mov rcx, [rbp + nb300nf_outeriter]
	mov rdx, [rbp + nb300nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 240
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
