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



.globl nb_kernel230_x86_64_sse
.globl _nb_kernel230_x86_64_sse
nb_kernel230_x86_64_sse:	
_nb_kernel230_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb230_fshift,           16
.equiv          nb230_gid,              24
.equiv          nb230_pos,              32
.equiv          nb230_faction,          40
.equiv          nb230_charge,           48
.equiv          nb230_p_facel,          56
.equiv          nb230_argkrf,           64
.equiv          nb230_argcrf,           72
.equiv          nb230_Vc,               80
.equiv          nb230_type,             88
.equiv          nb230_p_ntype,          96
.equiv          nb230_vdwparam,         104
.equiv          nb230_Vvdw,             112
.equiv          nb230_p_tabscale,       120
.equiv          nb230_VFtab,            128
.equiv          nb230_invsqrta,         136
.equiv          nb230_dvda,             144
.equiv          nb230_p_gbtabscale,     152
.equiv          nb230_GBtab,            160
.equiv          nb230_p_nthreads,       168
.equiv          nb230_count,            176
.equiv          nb230_mtx,              184
.equiv          nb230_outeriter,        192
.equiv          nb230_inneriter,        200
.equiv          nb230_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb230_ix,               0
.equiv          nb230_iy,               16
.equiv          nb230_iz,               32
.equiv          nb230_iq,               48
.equiv          nb230_dx,               64
.equiv          nb230_dy,               80
.equiv          nb230_dz,               96
.equiv          nb230_c6,               112
.equiv          nb230_c12,              128
.equiv          nb230_tsc,              144
.equiv          nb230_qq,               160
.equiv          nb230_vctot,            176
.equiv          nb230_Vvdwtot,          192
.equiv          nb230_fix,              208
.equiv          nb230_fiy,              224
.equiv          nb230_fiz,              240
.equiv          nb230_half,             256
.equiv          nb230_three,            272
.equiv          nb230_two,              288
.equiv          nb230_krf,              304
.equiv          nb230_crf,              320
.equiv          nb230_nri,              336
.equiv          nb230_iinr,             344
.equiv          nb230_jindex,           352
.equiv          nb230_jjnr,             360
.equiv          nb230_shift,            368
.equiv          nb230_shiftvec,         376
.equiv          nb230_facel,            384
.equiv          nb230_innerjjnr,        392
.equiv          nb230_is3,              400
.equiv          nb230_ii3,              404
.equiv          nb230_ntia,             408
.equiv          nb230_innerk,           412
.equiv          nb230_n,                416
.equiv          nb230_nn1,              420
.equiv          nb230_ntype,            424
.equiv          nb230_nouter,           428
.equiv          nb230_ninner,           432

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
	sub rsp, 448		;# local variable stack space (n*16)

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
	mov [rsp + nb230_nouter], eax
	mov [rsp + nb230_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb230_nri], edi
	mov [rsp + nb230_iinr], rsi
	mov [rsp + nb230_jindex], rdx
	mov [rsp + nb230_jjnr], rcx
	mov [rsp + nb230_shift], r8
	mov [rsp + nb230_shiftvec], r9
	mov rdi, [rbp + nb230_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb230_ntype], edi
	mov rsi, [rbp + nb230_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb230_facel], xmm0

	mov rax, [rbp + nb230_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb230_tsc], xmm3

	mov rsi, [rbp + nb230_argkrf]
	mov rdi, [rbp + nb230_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb230_krf], xmm1
	movaps [rsp + nb230_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb230_half], eax
	movss xmm1, [rsp + nb230_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb230_half],  xmm1
	movaps [rsp + nb230_two],  xmm2
	movaps [rsp + nb230_three],  xmm3

.nb230_threadloop:
        mov   rsi, [rbp + nb230_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb230_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb230_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb230_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb230_n], eax
        mov [rsp + nb230_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb230_outerstart
        jmp .nb230_end

.nb230_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb230_nouter]
	mov [rsp + nb230_nouter], ebx

.nb230_outer:
	mov   rax, [rsp + nb230_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb230_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb230_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb230_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb230_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb230_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb230_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb230_ntype]
    	shl   edx, 1
    	mov   [rsp + nb230_ntia], edx
		
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb230_pos]    ;# eax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb230_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb230_ix], xmm0
	movaps [rsp + nb230_iy], xmm1
	movaps [rsp + nb230_iz], xmm2

	mov   [rsp + nb230_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb230_vctot], xmm4
	movaps [rsp + nb230_Vvdwtot], xmm4
	movaps [rsp + nb230_fix], xmm4
	movaps [rsp + nb230_fiy], xmm4
	movaps [rsp + nb230_fiz], xmm4
	
	mov   rax, [rsp + nb230_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb230_pos]
	mov   rdi, [rbp + nb230_faction]	
	mov   rax, [rsp + nb230_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb230_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb230_ninner]
	mov   [rsp + nb230_ninner], ecx
	add   edx, 0
	mov   [rsp + nb230_innerk], edx    ;# number of innerloop atoms 
	jge   .nb230_unroll_loop
	jmp   .nb230_finish_inner
.nb230_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb230_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r12d, [rdx]	
	mov   r13d, [rdx + 4]              
	mov   r14d, [rdx + 8]            
	mov   r15d, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb230_innerjjnr],  16 ;# advance pointer (unrolled 4) 
    
	lea   rax, [r12 + r12*2]     ;# replace jnr with j3 
	lea   rbx, [r13 + r13*2]	
	lea   rcx, [r14 + r14*2]    
	lea   rdx, [r15 + r15*2]	

	mov rdi, [rbp + nb230_pos]
	;# load coordinates
	movlps xmm1, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + rcx*4]	;# x3 y3 - - 
	movhps xmm1, [rdi + rbx*4]	;# x2 y2 - -
	movhps xmm2, [rdi + rdx*4]	;# x4 y4 - -

	movss xmm5, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + rcx*4 + 8]	;# z2 - - - 
	movss xmm7, [rdi + rbx*4 + 8]	;# z3 - - - 
	movss xmm8, [rdi + rdx*4 + 8]	;# z4 - - - 
    movlhps xmm5, xmm7 ;# jzOa  -  jzOb  -
    movlhps xmm6, xmm8 ;# jzOc  -  jzOd -

    movaps xmm4, xmm1
    unpcklps xmm1, xmm2  ;# jxa jxc jya jyc        
    unpckhps xmm4, xmm2  ;# jxb jxd jyb jyd
    movaps xmm2, xmm1
    unpcklps xmm1, xmm4 ;# x
    unpckhps xmm2, xmm4 ;# y
    shufps   xmm5, xmm6,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

	mov rsi, [rbp + nb230_charge]

	;# calc dr  
	subps xmm1, [rsp + nb230_ix]
	subps xmm2, [rsp + nb230_iy]
	subps xmm5, [rsp + nb230_iz]

	movss xmm6, [rsi + r12*4]
	movss xmm3, [rsi + r13*4]
	movss xmm4, [rsi + r14*4]
	movss xmm0, [rsi + r15*4]

	;# store dr
    movaps [rsp + nb230_dx], xmm1
    movaps [rsp + nb230_dy], xmm2
    movaps [rsp + nb230_dz], xmm5
    
	;# square it 
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	mulps xmm5,xmm5
	addps xmm1, xmm2
	addps xmm1, xmm5

	mov rsi, [rbp + nb230_type]

    unpcklps xmm6, xmm4  ;# jqa jqc - -
    unpcklps xmm3, xmm0  ;# jqb jqd - -
    unpcklps xmm6, xmm3  ;# jqa jqb jqc jqd
	mulps xmm6, [rsp + nb230_iq] 
    movaps [rsp + nb230_qq], xmm6

	;# rsq in xmm1
    movaps xmm0, [rsp + nb230_krf]
    mulps  xmm0, xmm1    ;# krsq
    
     ;# vdw parameters
	mov r12d, [rsi + r12*4]
	mov r13d, [rsi + r13*4]
	mov r14d, [rsi + r14*4]
	mov r15d, [rsi + r15*4]

   ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb230_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm2	
	mulps xmm4, [rsp + nb230_half]	
	movaps xmm2, xmm4
	mulps  xmm1, xmm4	
    
    ;# xmm2=rinv
    ;# xmm1=r

	shl r12d, 1	
	shl r13d, 1	
	shl r14d, 1	
	shl r15d, 1	
    mov edi, [rsp + nb230_ntia]

    mulps xmm1, [rsp + nb230_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm1
    
	add r12d, edi
	add r13d, edi
	add r14d, edi
	add r15d, edi

    ;# convert back to float
    cvtdq2ps  xmm4, xmm5
    
    ;# multiply by 8
    pslld   xmm5, 3

    ;# calculate eps
    subps     xmm1, xmm4

    ;# move to integer registers
    movhlps xmm6, xmm5
    movd    r8d, xmm5
    movd    r10d, xmm6
    pshufd  xmm5, xmm5, 1
    pshufd  xmm6, xmm6, 1
    movd    r9d, xmm5
    movd    r11d, xmm6

    ;# xmm1=eps
    ;# xmm2=rinv
    
	mov rsi, [rbp + nb230_VFtab]
    ;# calculate LJ table
    movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r8*4 + 16]
    
	movlps xmm7,  [rsi + r10*4]
	movlps xmm11, [rsi + r10*4 + 16]

    movaps  xmm3, xmm2              ;# rinv
    subps   xmm3, xmm0
    subps   xmm3, xmm0              ;# rinv-2*krsq
    addps   xmm0, xmm2              ;# rinv+krsq
    subps   xmm0, [rsp + nb230_crf] ;# rinv+krsq-crf
    mulps   xmm0, [rsp + nb230_qq]  ;# vcoul=qq*(rinv+krsq-crf)
    mulps   xmm3, [rsp + nb230_qq]  ;# qq*(rinv-2*krsq)

    mulps   xmm3, xmm2              ;# rinv*qq*(rinv-2*krsq)
    
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r9*4 + 16]

    addps   xmm0, [rsp + nb230_vctot]
    movaps  [rsp + nb230_vctot], xmm0
    
	movhps xmm7,  [rsi + r11*4]
	movhps xmm11, [rsi + r11*4 + 16]
    
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101

	movlps xmm7,  [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r8*4 + 24]
    
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r10*4 + 24]

	movhps xmm7,  [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r9*4 + 24]
    
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r11*4 + 24]

    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
	mov rsi, [rbp + nb230_vdwparam]
    
    mulps  xmm7, xmm1    ;# Heps
    mulps  xmm11, xmm1 
    mulps  xmm6, xmm1   ;# Geps
    mulps  xmm10, xmm1 
    mulps  xmm7, xmm1   ;# Heps2
    mulps  xmm11, xmm1 

    ;# load c6/c12
	movlps xmm13, [rsi + r12*4]
	movlps xmm14, [rsi + r14*4]
	movhps xmm13, [rsi + r13*4]
	movhps xmm14, [rsi + r15*4]

    addps  xmm5, xmm6  ;# F+Geps
    addps  xmm9, xmm10 
    addps  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addps  xmm9, xmm11 
    addps  xmm7, xmm7    ;# 2*Heps2
    addps  xmm11, xmm11
    addps  xmm7, xmm6   ;# 2*Heps2+Geps
    addps  xmm11, xmm10

	movaps xmm12, xmm13
	shufps xmm12, xmm14, 136  ;# 10001000
	shufps xmm13, xmm14, 221  ;# 11011101

    addps  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm11, xmm9
    mulps  xmm5, xmm1  ;# eps*Fp
    mulps  xmm9, xmm1
    addps  xmm5, xmm4 ;# VV
    addps  xmm9, xmm8

	mov rsi, [rbp + nb230_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movlps xmm1, [rsi + rcx*4] ;# x3 y3 - -

    mulps  xmm5, xmm12  ;# VV*c6 = vnb6
    mulps  xmm9, xmm13  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb230_Vvdwtot]
    movaps [rsp + nb230_Vvdwtot], xmm5
        
    mulps  xmm7, xmm12   ;# FF*c6 = fnb6
    mulps  xmm11, xmm13   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + rdx*4] ;# x3 y3 x4 y4

    mulps  xmm7, [rsp + nb230_tsc]
    subps  xmm3, xmm7
    mulps  xmm3, xmm2   ;# fscal

    movaps xmm9, xmm3
    movaps xmm10, xmm3
    movaps xmm11, xmm3
    
    movaps xmm12, [rsp + nb230_fix]
    movaps xmm13, [rsp + nb230_fiy]
    movaps xmm14, [rsp + nb230_fiz]
    
    mulps  xmm9,  [rsp + nb230_dx]
    mulps  xmm10, [rsp + nb230_dy]
    mulps  xmm11, [rsp + nb230_dz]

    ;# accumulate i forces
    addps xmm12, xmm9
    addps xmm13, xmm10
    addps xmm14, xmm11
    movaps [rsp + nb230_fix], xmm12
    movaps [rsp + nb230_fiy], xmm13
    movaps [rsp + nb230_fiz], xmm14
    
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
	sub dword ptr [rsp + nb230_innerk],  4
	jl    .nb230_finish_inner
	jmp   .nb230_unroll_loop
.nb230_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb230_innerk],  4
	mov   edx, [rsp + nb230_innerk]
	and   edx, 2
	jnz   .nb230_dopair
	jmp   .nb230_checksingle
.nb230_dopair:	
    mov   rcx, [rsp + nb230_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb230_innerjjnr],  8	

	mov rsi, [rbp + nb230_charge]
	movss xmm0, [rsi + rax*4]
	movss xmm2, [rsi + rbx*4]

    unpcklps xmm0, xmm2  ;# jqa jqb 
	mulps xmm0, [rsp + nb230_iq] 
    movaps [rsp + nb230_qq], xmm0

	mov rsi, [rbp + nb230_type]
    ;# vdw parameters
	mov r12d, [rsi + rax*4]
	mov r13d, [rsi + rbx*4]
	shl r12d, 1	
	shl r13d, 1	
    mov edi, [rsp + nb230_ntia]
	add r12d, edi
	add r13d, edi

	mov rsi, [rbp + nb230_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movhps xmm3, [rsi + r13*4]

    xorps  xmm7, xmm7
	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101

    movaps [rsp + nb230_c6], xmm0
    movaps [rsp + nb230_c12], xmm3
    
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# load coordinates
	mov rdi, [rbp + nb230_pos]
    
	movlps xmm1, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + rbx*4]	;# x2 y2 - - 

	movss xmm5, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + rbx*4 + 8]	;# z2 - - - 

    unpcklps xmm1, xmm2 ;# x1 x2 y1 y2
    movhlps  xmm2, xmm1 ;# y1 y2 -  -
    unpcklps xmm5, xmm6 ;# z1 z2 -  -

	;# calc dr  
	subps xmm1, [rsp + nb230_ix]
	subps xmm2, [rsp + nb230_iy]
	subps xmm5, [rsp + nb230_iz]

	;# store dr
    movaps [rsp + nb230_dx], xmm1
    movaps [rsp + nb230_dy], xmm2
    movaps [rsp + nb230_dz], xmm5
    
	;# square it 
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	mulps xmm5,xmm5
	addps xmm1, xmm2
	addps xmm1, xmm5

	;# rsq in xmm1
    
    movaps xmm0, [rsp + nb230_krf]
    mulps  xmm0, xmm1    ;# krsq
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb230_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm2	
	mulps xmm4, [rsp + nb230_half]	
	movaps xmm2, xmm4
	mulps  xmm1, xmm4	
    ;# xmm2=rinv
    ;# xmm1=r

    mulps xmm1, [rsp + nb230_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm1
    
    ;# convert back to float
    cvtdq2ps  xmm4, xmm5
    
    ;# multiply by 8
    pslld   xmm5, 3

    ;# calculate eps
    subps     xmm1, xmm4

    ;# move to integer registers
    movd    r8d, xmm5
    pshufd  xmm5, xmm5, 1
    movd    r9d, xmm5

    ;# xmm1=eps
    ;# xmm2=rinv
    
	mov rsi, [rbp + nb230_VFtab]
    ;# calculate LJ table
    movlps xmm4, [rsi + r8*4] 
   	movlps xmm5, [rsi + r9*4] 

    unpcklps xmm4, xmm5 
    movhlps  xmm5, xmm4 
    
    movlps xmm6, [rsi + r8*4 + 8] 
   	movlps xmm7, [rsi + r9*4 + 8] 
    unpcklps xmm6, xmm7 
    movhlps  xmm7, xmm6 

    movaps  xmm3, xmm2              ;# rinv
    subps   xmm3, xmm0
    subps   xmm3, xmm0              ;# rinv-2*krsq
    addps   xmm0, xmm2              ;# rinv+krsq
    subps   xmm0, [rsp + nb230_crf] ;# rinv+krsq-crf
    mulps   xmm0, [rsp + nb230_qq]  ;# vcoul=qq*(rinv+krsq-crf)
    mulps   xmm3, [rsp + nb230_qq]  ;# qq*(rinv-2*krsq)

    mulps   xmm3, xmm2              ;# rinv*qq*(rinv-2*krsq)
    
    
    movlps xmm8, [rsi + r8*4 + 16] 
   	movlps xmm9, [rsi + r9*4 + 16] 

    unpcklps xmm8, xmm9 
    movhlps  xmm9, xmm8 
    
    addps   xmm0, [rsp + nb230_vctot]
    movlps  [rsp + nb230_vctot], xmm0
    
    movlps xmm10, [rsi + r8*4 + 24] 
   	movlps xmm11, [rsi + r9*4 + 24] 

    unpcklps xmm10, xmm11
    movhlps  xmm11, xmm10
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11
    
    mulps  xmm7, xmm1    ;# Heps
    mulps  xmm11, xmm1 
    mulps  xmm6, xmm1   ;# Geps
    mulps  xmm10, xmm1 
    mulps  xmm7, xmm1   ;# Heps2
    mulps  xmm11, xmm1 
    addps  xmm5, xmm6  ;# F+Geps
    addps  xmm9, xmm10 
    addps  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addps  xmm9, xmm11 
    addps  xmm7, xmm7    ;# 2*Heps2
    addps  xmm11, xmm11
    addps  xmm7, xmm6   ;# 2*Heps2+Geps
    addps  xmm11, xmm10
    
    addps  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addps  xmm11, xmm9
    mulps  xmm5, xmm1  ;# eps*Fp
    mulps  xmm9, xmm1
    movaps xmm12, [rsp + nb230_c6]
    movaps xmm13, [rsp + nb230_c12]
    addps  xmm5, xmm4 ;# VV
    addps  xmm9, xmm8

    mulps  xmm5, xmm12  ;# VV*c6 = vnb6
    mulps  xmm9, xmm13  ;# VV*c12 = vnb12
    addps  xmm5, xmm9
    addps  xmm5, [rsp + nb230_Vvdwtot]
    movlps [rsp + nb230_Vvdwtot], xmm5
        
    mulps  xmm7, xmm12   ;# FF*c6 = fnb6
    mulps  xmm11, xmm13   ;# FF*c12  = fnb12
    addps  xmm7, xmm11
    
    mulps  xmm7, [rsp + nb230_tsc]
    subps  xmm3, xmm7
    mulps  xmm3, xmm2   ;# fscal

    movaps xmm9, xmm3
    movaps xmm10, xmm3
    movaps xmm11, xmm3
    
    xorps  xmm8, xmm8

    movaps xmm12, [rsp + nb230_fix]
    movaps xmm13, [rsp + nb230_fiy]
    movaps xmm14, [rsp + nb230_fiz]
    
    mulps  xmm9,  [rsp + nb230_dx]
    mulps  xmm10, [rsp + nb230_dy]
    mulps  xmm11, [rsp + nb230_dz]

    movlhps xmm9, xmm8
    movlhps xmm10, xmm8
    movlhps xmm11, xmm8
    
    ;# accumulate i forces
    addps xmm12, xmm9
    addps xmm13, xmm10
    addps xmm14, xmm11
    movaps [rsp + nb230_fix], xmm12
    movaps [rsp + nb230_fiy], xmm13
    movaps [rsp + nb230_fiz], xmm14
    
	mov rsi, [rbp + nb230_faction]
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
	
.nb230_checksingle:				
	mov   edx, [rsp + nb230_innerk]
	and   edx, 1
	jnz    .nb230_dosingle
	jmp    .nb230_updateouterdata
.nb230_dosingle:
	mov rdi, [rbp + nb230_pos]
	mov   rcx, [rsp + nb230_innerjjnr]
	mov   eax, [rcx]	

	mov rsi, [rbp + nb230_charge]
	movss xmm0, [rsi + rax*4]

	mulss xmm0, [rsp + nb230_iq] 
    movaps [rsp + nb230_qq], xmm0

	mov rsi, [rbp + nb230_type]
    ;# vdw parameters
	mov r12d, [rsi + rax*4]
	shl r12d, 1	
    mov edi, [rsp + nb230_ntia]
	add r12d, edi

	mov rsi, [rbp + nb230_vdwparam]
	movss xmm0, [rsi + r12*4]
	movss xmm3, [rsi + r12*4 + 4]

    movaps [rsp + nb230_c6], xmm0
    movaps [rsp + nb230_c12], xmm3
    
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# load coordinates
	mov rdi, [rbp + nb230_pos]
	movss xmm1, [rdi + rax*4]	
	movss xmm2, [rdi + rax*4 + 4]
	movss xmm5, [rdi + rax*4 + 8]

	;# calc dr  
	subss xmm1, [rsp + nb230_ix]
	subss xmm2, [rsp + nb230_iy]
	subss xmm5, [rsp + nb230_iz]

	;# store dr
    movaps [rsp + nb230_dx], xmm1
    movaps [rsp + nb230_dy], xmm2
    movaps [rsp + nb230_dz], xmm5
    
	;# square it 
	mulss xmm1,xmm1
	mulss xmm2,xmm2
	mulss xmm5,xmm5
	addss xmm1, xmm2
	addss xmm1, xmm5

	;# rsq in xmm1
    movaps xmm0, [rsp + nb230_krf]
    mulps  xmm0, xmm1    ;# krsq
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtss xmm5, xmm1
	movaps xmm2, xmm5
	mulss xmm5, xmm5
	movaps xmm4, [rsp + nb230_three]
	mulss xmm5, xmm1	;# rsq*lu*lu 	
    subss xmm4, xmm5	;# 30-rsq*lu*lu 
	mulss xmm4, xmm2	
	mulss xmm4, [rsp + nb230_half]	
	movaps xmm2, xmm4
	mulss  xmm1, xmm4	
    ;# xmm2=rinv
    ;# xmm1=r

    mulss xmm1, [rsp + nb230_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttss2si r8d, xmm1
    
    ;# convert back to float
    cvtsi2ss  xmm4, r8d
    
    ;# multiply by 8
    shl       r8d, 3

    ;# calculate eps
    subss     xmm1, xmm4

    ;# xmm1=eps
    ;# xmm2=rinv
    
	mov rsi, [rbp + nb230_VFtab]
    ;# calculate LJ table
    movss xmm4, [rsi + r8*4] 
   	movss xmm5, [rsi + r8*4 + 4] 
    movss xmm6, [rsi + r8*4 + 8] 
   	movss xmm7, [rsi + r8*4 + 12] 
    movss xmm8, [rsi + r8*4 + 16] 
   	movss xmm9, [rsi + r8*4 + 20] 
    movss xmm10, [rsi + r8*4 + 24]  
   	movss xmm11, [rsi + r8*4 + 28] 
    ;# dispersion table in xmm4-xmm7, repulsion table in xmm8-xmm11

    movaps  xmm3, xmm2              ;# rinv
    subss   xmm3, xmm0
    subss   xmm3, xmm0              ;# rinv-2*krsq
    addss   xmm0, xmm2              ;# rinv+krsq
    subss   xmm0, [rsp + nb230_crf] ;# rinv+krsq-crf
    mulss   xmm0, [rsp + nb230_qq]  ;# vcoul=qq*(rinv+krsq-crf)
    mulss   xmm3, [rsp + nb230_qq]  ;# qq*(rinv-2*krsq)

    mulss   xmm3, xmm2              ;# rinv*qq*(rinv-2*krsq)
    
    addss   xmm0, [rsp + nb230_vctot]
    movss   [rsp + nb230_vctot], xmm0
    
    ;# calculate table interaction
    mulss  xmm7, xmm1    ;# Heps
    mulss  xmm11, xmm1 
    mulss  xmm6, xmm1   ;# Geps
    mulss  xmm10, xmm1 
    mulss  xmm7, xmm1   ;# Heps2
    mulss  xmm11, xmm1 
    addss  xmm5, xmm6  ;# F+Geps
    addss  xmm9, xmm10 
    addss  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addss  xmm9, xmm11 
    addss  xmm7, xmm7    ;# 2*Heps2
    addss  xmm11, xmm11
    addss  xmm7, xmm6   ;# 2*Heps2+Geps
    addss  xmm11, xmm10
    
    addss  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addss  xmm11, xmm9
    mulss  xmm5, xmm1  ;# eps*Fp
    mulss  xmm9, xmm1
    movaps xmm12, [rsp + nb230_c6]
    movaps xmm13, [rsp + nb230_c12]
    addss  xmm5, xmm4 ;# VV
    addss  xmm9, xmm8

    mulss  xmm5, xmm12  ;# VV*c6 = vnb6
    mulss  xmm9, xmm13  ;# VV*c12 = vnb12
    addss  xmm5, xmm9
    addss  xmm5, [rsp + nb230_Vvdwtot]
    movss  [rsp + nb230_Vvdwtot], xmm5
        
    mulss  xmm7, xmm12   ;# FF*c6 = fnb6
    mulss  xmm11, xmm13   ;# FF*c12  = fnb12
    addss  xmm7, xmm11
    
    mulss  xmm7, [rsp + nb230_tsc]
    subss  xmm3, xmm7
    mulss  xmm3, xmm2   ;# fscal

    movaps xmm9, xmm3
    movaps xmm10, xmm3
    movaps xmm11, xmm3
    
    movaps xmm12, [rsp + nb230_fix]
    movaps xmm13, [rsp + nb230_fiy]
    movaps xmm14, [rsp + nb230_fiz]
    
    mulss  xmm9,  [rsp + nb230_dx]
    mulss  xmm10, [rsp + nb230_dy]
    mulss  xmm11, [rsp + nb230_dz]
    
    ;# accumulate i forces
    addss xmm12, xmm9
    addss xmm13, xmm10
    addss xmm14, xmm11
    movaps [rsp + nb230_fix], xmm12
    movaps [rsp + nb230_fiy], xmm13
    movaps [rsp + nb230_fiz], xmm14
    
	mov rsi, [rbp + nb230_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + rax*4]
    addss  xmm10, [rsi + rax*4 + 4]
    addss  xmm11, [rsi + rax*4 + 8]
    movss  [rsi + rax*4],     xmm9
    movss  [rsi + rax*4 + 4], xmm10
    movss  [rsi + rax*4 + 8], xmm11
    
.nb230_updateouterdata:
	mov   ecx, [rsp + nb230_ii3]
	mov   rdi, [rbp + nb230_faction]
	mov   rsi, [rbp + nb230_fshift]
	mov   edx, [rsp + nb230_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb230_fix]
	movaps xmm1, [rsp + nb230_fiy]
	movaps xmm2, [rsp + nb230_fiz]

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
	mov esi, [rsp + nb230_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb230_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb230_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb230_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb230_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb230_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb230_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb230_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb230_n], esi
        jmp .nb230_outer
.nb230_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb230_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb230_end
        ;# non-zero, do one more workunit
        jmp   .nb230_threadloop
.nb230_end:
	mov eax, [rsp + nb230_nouter]
	mov ebx, [rsp + nb230_ninner]
	mov rcx, [rbp + nb230_outeriter]
	mov rdx, [rbp + nb230_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 448
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







.globl nb_kernel230nf_x86_64_sse
.globl _nb_kernel230nf_x86_64_sse
nb_kernel230nf_x86_64_sse:	
_nb_kernel230nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb230nf_fshift,         16
.equiv          nb230nf_gid,            24
.equiv          nb230nf_pos,            32
.equiv          nb230nf_faction,        40
.equiv          nb230nf_charge,         48
.equiv          nb230nf_p_facel,        56
.equiv          nb230nf_argkrf,         64
.equiv          nb230nf_argcrf,         72
.equiv          nb230nf_Vc,             80
.equiv          nb230nf_type,           88
.equiv          nb230nf_p_ntype,        96
.equiv          nb230nf_vdwparam,       104
.equiv          nb230nf_Vvdw,           112
.equiv          nb230nf_p_tabscale,     120
.equiv          nb230nf_VFtab,          128
.equiv          nb230nf_invsqrta,       136
.equiv          nb230nf_dvda,           144
.equiv          nb230nf_p_gbtabscale,   152
.equiv          nb230nf_GBtab,          160
.equiv          nb230nf_p_nthreads,     168
.equiv          nb230nf_count,          176
.equiv          nb230nf_mtx,            184
.equiv          nb230nf_outeriter,      192
.equiv          nb230nf_inneriter,      200
.equiv          nb230nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb230nf_ix,             0
.equiv          nb230nf_iy,             16
.equiv          nb230nf_iz,             32
.equiv          nb230nf_iq,             48
.equiv          nb230nf_c6,             64
.equiv          nb230nf_c12,            80
.equiv          nb230nf_vctot,          96
.equiv          nb230nf_Vvdwtot,        112
.equiv          nb230nf_half,           128
.equiv          nb230nf_three,          144
.equiv          nb230nf_krf,            160
.equiv          nb230nf_crf,            176
.equiv          nb230nf_tsc,            192
.equiv          nb230nf_nri,            208
.equiv          nb230nf_iinr,           216
.equiv          nb230nf_jindex,         224
.equiv          nb230nf_jjnr,           232
.equiv          nb230nf_shift,          240
.equiv          nb230nf_shiftvec,       248
.equiv          nb230nf_facel,          256
.equiv          nb230nf_innerjjnr,      264
.equiv          nb230nf_is3,            272
.equiv          nb230nf_ii3,            280
.equiv          nb230nf_ntia,           284
.equiv          nb230nf_innerk,         288
.equiv          nb230nf_n,              292
.equiv          nb230nf_nn1,            296
.equiv          nb230nf_ntype,          300
.equiv          nb230nf_nouter,         304
.equiv          nb230nf_ninner,         308

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
	sub rsp, 328		;# local variable stack space (n*16+8)

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
	mov [rsp + nb230nf_nouter], eax
	mov [rsp + nb230nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb230nf_nri], edi
	mov [rsp + nb230nf_iinr], rsi
	mov [rsp + nb230nf_jindex], rdx
	mov [rsp + nb230nf_jjnr], rcx
	mov [rsp + nb230nf_shift], r8
	mov [rsp + nb230nf_shiftvec], r9
	mov rdi, [rbp + nb230nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb230nf_ntype], edi
	mov rsi, [rbp + nb230nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb230nf_facel], xmm0

	mov rax, [rbp + nb230nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb230nf_tsc], xmm3

	mov rsi, [rbp + nb230nf_argkrf]
	mov rdi, [rbp + nb230nf_argcrf]
	movss xmm1, [rsi]
	movss xmm2, [rdi]
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	movaps [rsp + nb230nf_krf], xmm1
	movaps [rsp + nb230nf_crf], xmm2

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb230nf_half], eax
	movss xmm1, [rsp + nb230nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb230nf_half],  xmm1
	movaps [rsp + nb230nf_three],  xmm3


.nb230nf_threadloop:
        mov   rsi, [rbp + nb230nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb230nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb230nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb230nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb230nf_n], eax
        mov [rsp + nb230nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb230nf_outerstart
        jmp .nb230nf_end

.nb230nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb230nf_nouter]
	mov [rsp + nb230nf_nouter], ebx

.nb230nf_outer:
	mov   rax, [rsp + nb230nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb230nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb230nf_shiftvec]   ;# eax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb230nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb230nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb230nf_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb230nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb230nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb230nf_ntia], edx
		
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb230nf_pos]    ;# eax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb230nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb230nf_ix], xmm0
	movaps [rsp + nb230nf_iy], xmm1
	movaps [rsp + nb230nf_iz], xmm2

	mov   [rsp + nb230nf_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb230nf_vctot], xmm4
	movaps [rsp + nb230nf_Vvdwtot], xmm4
		
	mov   rax, [rsp + nb230nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb230nf_pos]
	mov   rax, [rsp + nb230nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb230nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb230nf_ninner]
	mov   [rsp + nb230nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb230nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb230nf_unroll_loop
	jmp   .nb230nf_finish_inner
.nb230nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb230nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb230nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb230nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm2, [rsp + nb230nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# constant 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx
	
	mov rsi, [rbp + nb230nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb230nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb230nf_ntia]
	add eax, edi
	add ebx, edi
	add ecx, edi
	add edx, edi

	movlps xmm6, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm6, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4]

	movaps xmm4, xmm6
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm6, xmm7, 221  ;# constant 11011101
	
	movd  eax, mm0		
	movd  ebx, mm1
	movd  ecx, mm2
	movd  edx, mm3

	movaps [rsp + nb230nf_c6], xmm4
	movaps [rsp + nb230nf_c12], xmm6
	
	mov rsi, [rbp + nb230nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	mulps xmm3, xmm2
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

	shufps xmm2, xmm6, 136  ;# constant 10001000
	
	shufps xmm0, xmm5, 136  ;# constant 10001000
	shufps xmm1, xmm5, 221  ;# constant 11011101		

	;# move ix-iz to xmm4-xmm6 
	movaps xmm4, [rsp + nb230nf_ix]
	movaps xmm5, [rsp + nb230nf_iy]
	movaps xmm6, [rsp + nb230nf_iz]

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
	
	movaps xmm7, [rsp + nb230nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb230nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb230nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 	
	movaps xmm1, xmm0
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb230nf_crf]
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 

	addps  xmm6, [rsp + nb230nf_vctot]
	movaps [rsp + nb230nf_vctot], xmm6
	
	;# LJ table
	mulps  xmm4, xmm1  ;# r
	mulps  xmm4, [rsp + nb230nf_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
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

	mov  rsi, [rbp + nb230nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movlps xmm7, [rsi + rcx*4]
	movhps xmm5, [rsi + rbx*4]
	movhps xmm7, [rsi + rdx*4] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movlps xmm3, [rsi + rcx*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb230nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 
	
	;# Update Vvdwtot directly 
	addps  xmm5, [rsp + nb230nf_Vvdwtot]
	movaps [rsp + nb230nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movlps xmm7, [rsi + rcx*4 + 16]
	movhps xmm5, [rsi + rbx*4 + 16]
	movhps xmm7, [rsi + rdx*4 + 16] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm3, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm3, [rsi + rdx*4 + 24] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb230nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 
	
	addps  xmm5, [rsp + nb230nf_Vvdwtot]
	movaps [rsp + nb230nf_Vvdwtot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb230nf_innerk],  4
	jl    .nb230nf_finish_inner
	jmp   .nb230nf_unroll_loop
.nb230nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb230nf_innerk],  4
	mov   edx, [rsp + nb230nf_innerk]
	and   edx, 2
	jnz   .nb230nf_dopair
	jmp   .nb230nf_checksingle
.nb230nf_dopair:	
	mov rsi, [rbp + nb230nf_charge]

    mov   rcx, [rsp + nb230nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb230nf_innerjjnr],  8

	xorps xmm3, xmm3
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 12 ;# constant 00001100 
	shufps xmm3, xmm3, 88 ;# constant 01011000 ;# xmm3(0,1) has the charges 

	mov rsi, [rbp + nb230nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]	
	mov rsi, [rbp + nb230nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb230nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [rsi + rcx*4]
	movhps xmm6, [rsi + rdx*4]
	mov rdi, [rbp + nb230nf_pos]	
	xorps  xmm7,xmm7
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# constant 00001000 	
	shufps xmm6, xmm6, 13 ;# constant 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [rsp + nb230nf_c6], xmm4
	movaps [rsp + nb230nf_c12], xmm6	
			
	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	;# move coordinates to xmm0-xmm2 
	movlps xmm1, [rdi + rax*4]
	movss xmm2, [rdi + rax*4 + 8]	
	movhps xmm1, [rdi + rbx*4]
	movss xmm0, [rdi + rbx*4 + 8]	

	mulps  xmm3, [rsp + nb230nf_iq]

	movlhps xmm3, xmm7
	
	shufps xmm2, xmm0, 0
	
	movaps xmm0, xmm1

	shufps xmm2, xmm2, 136  ;# constant 10001000
	
	shufps xmm0, xmm0, 136  ;# constant 10001000
	shufps xmm1, xmm1, 221  ;# constant 11011101
			
	;# move ix-iz to xmm4-xmm6 
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb230nf_ix]
	movaps xmm5, [rsp + nb230nf_iy]
	movaps xmm6, [rsp + nb230nf_iz]

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

	movaps xmm7, [rsp + nb230nf_krf]
	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb230nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb230nf_half]
	mulps  xmm7, xmm4	;# xmm7=krsq 
	subps xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 	
	movaps xmm1, xmm0
	movaps xmm6, xmm0
	addps  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subps  xmm6, [rsp + nb230nf_crf]
	mulps  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	addps  xmm6, [rsp + nb230nf_vctot]
	movaps [rsp + nb230nf_vctot], xmm6

	;# LJ table
	mulps  xmm4, xmm1  ;# r
	mulps  xmm4, [rsp + nb230nf_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	mov  rsi, [rbp + nb230nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movhps xmm5, [rsi + rbx*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movhps xmm7, [rsi + rbx*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb230nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 
	
	;# Update Vvdwtot directly 
	addps  xmm5, [rsp + nb230nf_Vvdwtot]
	movaps [rsp + nb230nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movhps xmm5, [rsi + rbx*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb230nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 
	
	addps  xmm5, [rsp + nb230nf_Vvdwtot]
	movaps [rsp + nb230nf_Vvdwtot], xmm5

.nb230nf_checksingle:				
	mov   edx, [rsp + nb230nf_innerk]
	and   edx, 1
	jnz    .nb230nf_dosingle
	jmp    .nb230nf_updateouterdata
.nb230nf_dosingle:			
	mov rsi, [rbp + nb230nf_charge]
	mov rdi, [rbp + nb230nf_pos]
	mov   rcx, [rsp + nb230nf_innerjjnr]
	xorps xmm3, xmm3
	mov   eax, [rcx]
	movss xmm3, [rsi + rax*4]	;# xmm3(0) has the charge 	

	mov rsi, [rbp + nb230nf_type]
	mov ecx, eax
	mov ecx, [rsi + rcx*4]	
	mov rsi, [rbp + nb230nf_vdwparam]
	shl ecx, 1
	add ecx, [rsp + nb230nf_ntia]
	xorps  xmm6, xmm6
	movlps xmm6, [rsi + rcx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# constant 11111100	
	shufps xmm6, xmm6, 253  ;# constant 11111101	
			
	movaps [rsp + nb230nf_c6], xmm4
	movaps [rsp + nb230nf_c12], xmm6	
		
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	
 
	mulps  xmm3, [rsp + nb230nf_iq]
	
	xorps   xmm7, xmm7
	
	movaps xmm4, [rsp + nb230nf_ix]
	movaps xmm5, [rsp + nb230nf_iy]
	movaps xmm6, [rsp + nb230nf_iz]

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

	movss xmm7, [rsp + nb230nf_krf]
	rsqrtss xmm5, xmm4
	;# lookup seed in xmm5 
	movss xmm2, xmm5
	mulss xmm5, xmm5
	movss xmm1, [rsp + nb230nf_three]
	mulss xmm5, xmm4	;# rsq*lu*lu 			
	movss xmm0, [rsp + nb230nf_half]
	mulss  xmm7, xmm4	;# xmm7=krsq 
	subss xmm1, xmm5	;# constant 30-rsq*lu*lu 
	mulss xmm1, xmm2	
	mulss xmm0, xmm1	;# xmm0=rinv 	
	movss xmm1, xmm0
	movss xmm6, xmm0
	addss  xmm6, xmm7	;# xmm6=rinv+ krsq 
	subss  xmm6, [rsp + nb230nf_crf]
	mulss  xmm6, xmm3	;# xmm6=vcoul=qq*(rinv+ krsq-crf) 
	addss  xmm6, [rsp + nb230nf_vctot]
	movss [rsp + nb230nf_vctot], xmm6
	
	;# LJ table
	mulss  xmm4, xmm1  ;# r
	mulss  xmm4, [rsp + nb230nf_tsc] ;# rtab
	
	movaps xmm0, xmm1 ;# copy of rinv
	cvttps2pi mm6, xmm4
	cvtpi2ps xmm6, mm6
	subss xmm4, xmm6	
	movss xmm1, xmm4	;# xmm1=eps 
	movss xmm2, xmm1	
	mulss  xmm2, xmm2	;# xmm2=eps2 
	pslld mm6, 3

	movd mm0, eax	

	mov  rsi, [rbp + nb230nf_VFtab]
	movd eax, mm6
	
	;# dispersion 
	movlps xmm5, [rsi + rax*4]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101
	
	movlps xmm7, [rsi + rax*4 + 8]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# dispersion table ready, in xmm4-xmm7 	

	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 

	movss  xmm4, [rsp + nb230nf_c6]
	mulss  xmm5, xmm4	 ;# Vvdw6 
	
	;# Update Vvdwtot directly 
	addss  xmm5, [rsp + nb230nf_Vvdwtot]
	movss [rsp + nb230nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 16]
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# constant 10001000
	shufps xmm5, xmm7, 221  ;# constant 11011101

	movlps xmm7, [rsi + rax*4 + 24]
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# constant 10001000
	shufps xmm7, xmm3, 221  ;# constant 11011101
	;# table ready, in xmm4-xmm7 	
	mulss  xmm6, xmm1	;# xmm6=Geps 
	mulss  xmm7, xmm2	;# xmm7=Heps2 
	addss  xmm5, xmm6
	addss  xmm5, xmm7	;# xmm5=Fp 	
	mulss  xmm5, xmm1 ;# xmm5=eps*Fp 
	addss  xmm5, xmm4 ;# xmm5=VV 
 	
	movss  xmm4, [rsp + nb230nf_c12]
	mulss  xmm5, xmm4 ;# Vvdw12 
	
	addss  xmm5, [rsp + nb230nf_Vvdwtot]
	movss [rsp + nb230nf_Vvdwtot], xmm5


.nb230nf_updateouterdata:

	;# get n from stack
	mov esi, [rsp + nb230nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb230nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb230nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb230nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb230nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb230nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb230nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb230nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb230nf_n], esi
        jmp .nb230nf_outer
.nb230nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb230nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb230nf_end
        ;# non-zero, do one more workunit
        jmp   .nb230nf_threadloop
.nb230nf_end:

	mov eax, [rsp + nb230nf_nouter]
	mov ebx, [rsp + nb230nf_ninner]
	mov rcx, [rbp + nb230nf_outeriter]
	mov rdx, [rbp + nb230nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 328
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
