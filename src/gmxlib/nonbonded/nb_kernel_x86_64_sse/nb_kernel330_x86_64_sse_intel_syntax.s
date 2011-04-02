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
%define .long		ddq
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




.globl nb_kernel330_x86_64_sse
.globl _nb_kernel330_x86_64_sse
nb_kernel330_x86_64_sse:	
_nb_kernel330_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb330_fshift,           16
.equiv          nb330_gid,              24
.equiv          nb330_pos,              32
.equiv          nb330_faction,          40
.equiv          nb330_charge,           48
.equiv          nb330_p_facel,          56
.equiv          nb330_argkrf,           64
.equiv          nb330_argcrf,           72
.equiv          nb330_Vc,               80
.equiv          nb330_type,             88
.equiv          nb330_p_ntype,          96
.equiv          nb330_vdwparam,         104
.equiv          nb330_Vvdw,             112
.equiv          nb330_p_tabscale,       120
.equiv          nb330_VFtab,            128
.equiv          nb330_invsqrta,         136
.equiv          nb330_dvda,             144
.equiv          nb330_p_gbtabscale,     152
.equiv          nb330_GBtab,            160
.equiv          nb330_p_nthreads,       168
.equiv          nb330_count,            176
.equiv          nb330_mtx,              184
.equiv          nb330_outeriter,        192
.equiv          nb330_inneriter,        200
.equiv          nb330_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb330_ix,               0
.equiv          nb330_iy,               16
.equiv          nb330_iz,               32
.equiv          nb330_iq,               48
.equiv          nb330_dx,               64
.equiv          nb330_dy,               80
.equiv          nb330_dz,               96
.equiv          nb330_rinv,             112
.equiv          nb330_tsc,              128
.equiv          nb330_qq,               144
.equiv          nb330_c6,               160
.equiv          nb330_c12,              176
.equiv          nb330_eps,              192
.equiv          nb330_vctot,            208
.equiv          nb330_Vvdwtot,          224
.equiv          nb330_fix,              240
.equiv          nb330_fiy,              256
.equiv          nb330_fiz,              272
.equiv          nb330_half,             288
.equiv          nb330_three,            304
.equiv          nb330_nri,              320
.equiv          nb330_iinr,             328
.equiv          nb330_jindex,           336
.equiv          nb330_jjnr,             344
.equiv          nb330_shift,            352
.equiv          nb330_shiftvec,         360
.equiv          nb330_facel,            368
.equiv          nb330_innerjjnr,        376
.equiv          nb330_is3,              384
.equiv          nb330_ii3,              388
.equiv          nb330_ntia,             392
.equiv          nb330_innerk,           396
.equiv          nb330_n,                400
.equiv          nb330_nn1,              404
.equiv          nb330_ntype,            408
.equiv          nb330_nouter,           412
.equiv          nb330_ninner,           416

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
	sub rsp, 432		;# local variable stack space (n*16+8)
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
	mov [rsp + nb330_nouter], eax
	mov [rsp + nb330_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb330_nri], edi
	mov [rsp + nb330_iinr], rsi
	mov [rsp + nb330_jindex], rdx
	mov [rsp + nb330_jjnr], rcx
	mov [rsp + nb330_shift], r8
	mov [rsp + nb330_shiftvec], r9
	mov rdi, [rbp + nb330_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb330_ntype], edi
	mov rsi, [rbp + nb330_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb330_facel], xmm0

	mov rax, [rbp + nb330_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb330_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb330_half], eax
	movss xmm1, [rsp + nb330_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb330_half],  xmm1
	movaps [rsp + nb330_three],  xmm3

.nb330_threadloop:
        mov   rsi, [rbp + nb330_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb330_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb330_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb330_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb330_n], eax
        mov [rsp + nb330_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb330_outerstart
        jmp .nb330_end

.nb330_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb330_nouter]
	mov [rsp + nb330_nouter], ebx

.nb330_outer:
	mov   rax, [rsp + nb330_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb330_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb330_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb330_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb330_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb330_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb330_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb330_ntype]
    	shl   edx, 1
    	mov   [rsp + nb330_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb330_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb330_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb330_ix], xmm0
	movaps [rsp + nb330_iy], xmm1
	movaps [rsp + nb330_iz], xmm2

	mov   [rsp + nb330_ii3], ebx
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb330_vctot], xmm4
	movaps [rsp + nb330_Vvdwtot], xmm4
	movaps [rsp + nb330_fix], xmm4
	movaps [rsp + nb330_fiy], xmm4
	movaps [rsp + nb330_fiz], xmm4
	
	mov   rax, [rsp + nb330_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

        mov esi, [rsp + nb330_charge]
        mov edi, [rsp + nb330_ntia]

	mov   rsi, [rbp + nb330_pos]
	mov   rdi, [rbp + nb330_faction]	
	mov   rax, [rsp + nb330_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb330_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb330_ninner]
	mov   [rsp + nb330_ninner], ecx
	add   edx, 0
	mov   [rsp + nb330_innerk], edx    ;# number of innerloop atoms 
	jge   .nb330_unroll_loop
	jmp   .nb330_finish_inner
.nb330_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb330_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r12d, [rdx]	
	mov   r13d, [rdx + 4]              
	mov   r14d, [rdx + 8]            
	mov   r15d, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb330_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	lea   rax, [r12 + r12*2]     ;# replace jnr with j3 
	lea   rbx, [r13 + r13*2]	
	lea   rcx, [r14 + r14*2]    
	lea   rdx, [r15 + r15*2]	

	mov rdi, [rbp + nb330_pos]
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
	mov rsi, [rbp + nb330_charge]

	;# calc dr  
	subps xmm1, [rsp + nb330_ix]
	subps xmm2, [rsp + nb330_iy]
	subps xmm5, [rsp + nb330_iz]

	;# store dr
    movaps [rsp + nb330_dx], xmm1
    movaps [rsp + nb330_dy], xmm2
    movaps [rsp + nb330_dz], xmm5
    
	movss xmm0, [rsi + r12*4]
	movss xmm3, [rsi + r14*4]
	movss xmm4, [rsi + r13*4]
	movss xmm6, [rsi + r15*4]
    
	;# square it 
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	mulps xmm5,xmm5
	addps xmm1, xmm2
	addps xmm1, xmm5
	;# rsq in xmm1

    unpcklps xmm0, xmm3  ;# jqa jqc - -
    unpcklps xmm4, xmm6  ;# jqb jqd - -
    unpcklps xmm0, xmm4  ;# jqa jqb jqc jqd
	mulps xmm0, [rsp + nb330_iq] 
    movaps [rsp + nb330_qq], xmm0

	mov rsi, [rbp + nb330_type]
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm12, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb330_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm12	
	mulps xmm4, [rsp + nb330_half]	
	movaps xmm15, xmm4
	mulps  xmm1, xmm4	
    ;# xmm15=rinv
    ;# xmm1=r

    ;# vdw parameters
	mov r12d, [rsi + r12*4]
	mov r13d, [rsi + r13*4]
	mov r14d, [rsi + r14*4]
	mov r15d, [rsi + r15*4]
 
    mulps xmm1, [rsp + nb330_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm1
    
	shl r12d, 1	
	shl r13d, 1	
	shl r14d, 1	
	shl r15d, 1	
    mov edi, [rsp + nb330_ntia]

    ;# convert back to float
    cvtdq2ps  xmm4, xmm5
    
    ;# multiply by 4
    pslld   xmm5, 2

    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm6, xmm5
    pslld   xmm5, 1
    paddd   xmm5, xmm6

    ;# calculate eps
    subps     xmm1, xmm4

	add r12d, edi
	add r13d, edi
	add r14d, edi
	add r15d, edi

    ;# move to integer registers
    movhlps xmm6, xmm5
    movd    r8d, xmm5
    movd    r10d, xmm6
    pshufd  xmm5, xmm5, 1
    pshufd  xmm6, xmm6, 1
    movd    r9d, xmm5
    movd    r11d, xmm6

	mov rsi, [rbp + nb330_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movlps xmm7, [rsi + r14*4]
	movhps xmm3, [rsi + r13*4]
	movhps xmm7, [rsi + r15*4]

	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101

    movaps [rsp + nb330_c6], xmm0
    movaps [rsp + nb330_c12], xmm3

    movaps [rsp + nb330_eps], xmm1
    ;# xmm15=rinv
    
	mov rsi, [rbp + nb330_VFtab]
    ;# load Coulomb and LJ table data in parallel
   	movlps xmm1, [rsi + r8*4]        ;# Y1c F1c 
   	movlps xmm5, [rsi + r8*4 + 16]   ;# Y1d F1d 
   	movlps xmm9, [rsi + r8*4 + 32]   ;# Y1r F1r 

	movlps xmm3, [rsi + r10*4]       ;# Y3c F3c 
	movlps xmm7, [rsi + r10*4 + 16]  ;# Y3d F3d 
	movlps xmm11, [rsi + r10*4 + 32] ;# Y3r F3r 

	movhps xmm1, [rsi + r9*4]        ;# Y1c F1c Y2c F2c
	movhps xmm5, [rsi + r9*4 + 16]   ;# Y1d F1d Y2d F2d
	movhps xmm9, [rsi + r9*4 + 32]   ;# Y1r F1r Y2r F2r

	movhps xmm3, [rsi + r11*4]       ;# Y3c F3c Y4c F4c
	movhps xmm7, [rsi + r11*4 + 16]  ;# Y3d F3d Y4d F4d
	movhps xmm11, [rsi + r11*4 + 32] ;# Y3r F3r Y4r F4r

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000   => Y1c Y2c Y3c Y4c
	shufps xmm4, xmm7, 136  ;# 10001000   => Y1d Y2d Y3d Y4d
	shufps xmm8, xmm11, 136  ;# 10001000  => Y1r Y2r Y3r Y4r
	shufps xmm1, xmm3, 221  ;# 11011101   => F1c F2c F3c F4c
	shufps xmm5, xmm7, 221  ;# 11011101   => F1d F2d F3d F4d
	shufps xmm9, xmm11, 221  ;# 11011101  => F1r F2r F3r F4r
    
   	movlps xmm3, [rsi + r8*4 + 8]     ;# G1c H1c 
   	movlps xmm7, [rsi + r8*4 + 24]    ;# G1d H1d 
   	movlps xmm11, [rsi + r8*4 + 40]   ;# G1r H1r 

	movlps xmm12, [rsi + r10*4 + 8]   ;# G3c H3c 
	movlps xmm13, [rsi + r10*4 + 24]  ;# G3d H3d 
	movlps xmm14, [rsi + r10*4 + 40]  ;# G3r H3r 

	movhps xmm3, [rsi + r9*4 + 8]     ;# G1c H1c G2c H2c
	movhps xmm7, [rsi + r9*4 + 24]    ;# G1d H1d G2d H2d
	movhps xmm11, [rsi + r9*4 + 40]   ;# G1r H1r G2r H2r

	movhps xmm12, [rsi + r11*4 + 8]   ;# G3c H3c G4c H4c
	movhps xmm13, [rsi + r11*4 + 24]  ;# G3d H3d G4d H4d
	movhps xmm14, [rsi + r11*4 + 40]  ;# G3r H3r G4r H4r


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

    movaps xmm12, [rsp + nb330_eps]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm12
    mulps  xmm11, xmm12
    mulps  xmm2, xmm12     ;# Geps
    mulps  xmm6, xmm12
    mulps  xmm10, xmm12
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm12
    mulps  xmm11, xmm12

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
    mulps  xmm5, xmm12
    mulps  xmm9, xmm12
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, [rsp + nb330_qq]   ;# VV*qq = vcoul
    mulps  xmm5, [rsp + nb330_c6]   ;# vnb6
    mulps  xmm9, [rsp + nb330_c12]   ;# vnb12
    mulps  xmm3, [rsp + nb330_qq]    ;# FF*qq = fij
    mulps  xmm7, [rsp + nb330_c6]   ;# fijD
    mulps  xmm11, [rsp + nb330_c12]   ;#fijR

    ;# accumulate vctot
    addps  xmm1, [rsp + nb330_vctot]

    ;# accumulate Vvdwtot
    addps  xmm5, [rsp + nb330_Vvdwtot]
    addps  xmm5, xmm9
    xorps  xmm9, xmm9

    movaps [rsp + nb330_vctot], xmm1
    movaps [rsp + nb330_Vvdwtot], xmm5

    
	mov rsi, [rbp + nb330_faction]
	;# the fj's - start by accumulating x & y forces from memory 
	movlps xmm0, [rsi + rax*4] ;# x1 y1 - -
	movlps xmm1, [rsi + rcx*4] ;# x3 y3 - -
    addps  xmm3, xmm7
    addps  xmm3, xmm11
    mulps  xmm3, xmm15
    mulps  xmm3, [rsp + nb330_tsc]  ;# fscal
    
    subps  xmm9, xmm3
    movaps xmm10, xmm9
    movaps xmm11, xmm9
    
	movhps xmm0, [rsi + rbx*4] ;# x1 y1 x2 y2
	movhps xmm1, [rsi + rdx*4] ;# x3 y3 x4 y4

    movaps xmm12, [rsp + nb330_fix]
    movaps xmm13, [rsp + nb330_fiy]
    movaps xmm14, [rsp + nb330_fiz]
    
    mulps  xmm9,  [rsp + nb330_dx]
    mulps  xmm10, [rsp + nb330_dy]
    mulps  xmm11, [rsp + nb330_dz]

    ;# accumulate i forces
    addps xmm12, xmm9
    addps xmm13, xmm10
    addps xmm14, xmm11
    movaps [rsp + nb330_fix], xmm12
    movaps [rsp + nb330_fiy], xmm13
    movaps [rsp + nb330_fiz], xmm14
    
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
	sub dword ptr [rsp + nb330_innerk],  4
	jl    .nb330_finish_inner
	jmp   .nb330_unroll_loop
.nb330_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb330_innerk],  4
	mov   edx, [rsp + nb330_innerk]
	and   edx, 2
	jnz   .nb330_dopair
	jmp   .nb330_checksingle
.nb330_dopair:	
    mov   rcx, [rsp + nb330_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb330_innerjjnr],  8	

	mov rsi, [rbp + nb330_charge]
	movss xmm0, [rsi + rax*4]
	movss xmm2, [rsi + rbx*4]

    unpcklps xmm0, xmm2  ;# jqa jqb 
	mulps xmm0, [rsp + nb330_iq] 
    movaps [rsp + nb330_qq], xmm0

	mov rsi, [rbp + nb330_type]
    ;# vdw parameters
	mov r12d, [rsi + rax*4]
	mov r13d, [rsi + rbx*4]
	shl r12d, 1	
	shl r13d, 1	
    mov edi, [rsp + nb330_ntia]
	add r12d, edi
	add r13d, edi

	mov rsi, [rbp + nb330_vdwparam]
	movlps xmm3, [rsi + r12*4]
	movhps xmm3, [rsi + r13*4]

    xorps  xmm7, xmm7
	movaps xmm0, xmm3
	shufps xmm0, xmm7, 136  ;# 10001000
	shufps xmm3, xmm7, 221  ;# 11011101

    movaps [rsp + nb330_c6], xmm0
    movaps [rsp + nb330_c12], xmm3
    
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# load coordinates
	mov rdi, [rbp + nb330_pos]
    
	movlps xmm1, [rdi + rax*4]	;# x1 y1 - - 
	movlps xmm2, [rdi + rbx*4]	;# x2 y2 - - 

	movss xmm5, [rdi + rax*4 + 8]	;# z1 - - - 
	movss xmm6, [rdi + rbx*4 + 8]	;# z2 - - - 

    unpcklps xmm1, xmm2 ;# x1 x2 y1 y2
    movhlps  xmm2, xmm1 ;# y1 y2 -  -
    unpcklps xmm5, xmm6 ;# z1 z2 -  -

	;# calc dr  
	subps xmm1, [rsp + nb330_ix]
	subps xmm2, [rsp + nb330_iy]
	subps xmm5, [rsp + nb330_iz]

	;# store dr
    movaps [rsp + nb330_dx], xmm1
    movaps [rsp + nb330_dy], xmm2
    movaps [rsp + nb330_dz], xmm5
    
	;# square it 
	mulps xmm1,xmm1
	mulps xmm2,xmm2
	mulps xmm5,xmm5
	addps xmm1, xmm2
	addps xmm1, xmm5

	;# rsq in xmm1
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtps xmm5, xmm1
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm4, [rsp + nb330_three]
	mulps xmm5, xmm1	;# rsq*lu*lu 	
    subps xmm4, xmm5	;# 30-rsq*lu*lu 
	mulps xmm4, xmm2	
	mulps xmm4, [rsp + nb330_half]	
	movaps xmm15, xmm4
	mulps  xmm1, xmm4	
    ;# xmm15=rinv
    ;# xmm1=r

    mulps xmm1, [rsp + nb330_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttps2dq xmm5, xmm1
    
    ;# convert back to float
    cvtdq2ps  xmm4, xmm5
    
    ;# multiply by 4
    pslld   xmm5, 2

    ;# multiply by three (copy, mult. by two, add back)
    movaps  xmm6, xmm5
    pslld   xmm5, 1
    paddd   xmm5, xmm6

    ;# calculate eps
    subps     xmm1, xmm4

    ;# move to integer registers
    movd    r8d, xmm5
    pshufd  xmm5, xmm5, 1
    movd    r9d, xmm5

    movaps [rsp + nb330_eps], xmm1
    ;# xmm15=rinv
    
	mov rsi, [rbp + nb330_VFtab]
    ;# load Coulomb and LJ table data in parallel
   	movlps xmm0, [rsi + r8*4]        ;# Y1c F1c 
   	movlps xmm4, [rsi + r8*4 + 16]   ;# Y1d F1d 
   	movlps xmm8, [rsi + r8*4 + 32]   ;# Y1r F1r 

	movlps xmm1, [rsi + r9*4]        ;# Y2c F2c
	movlps xmm5, [rsi + r9*4 + 16]   ;# Y2d F2d
	movlps xmm9, [rsi + r9*4 + 32]  ;# Y2r F2r

    unpcklps xmm0, xmm1 ;# Y1c Y2c F1c F2c
    unpcklps xmm4, xmm5 ;# Y1d Y2d F1d F2d
    unpcklps xmm8, xmm9 ;# Y1r Y2r F1r F2r
    movhlps  xmm1, xmm0 ;# F1c F2c
    movhlps  xmm5, xmm4 ;# F1d F2d
    movhlps  xmm9, xmm8 ;# F1r F2r
    
   	movlps xmm2, [rsi + r8*4 + 8]    ;# G1c H1c 
   	movlps xmm6, [rsi + r8*4 + 24]   ;# G1d H1d 
   	movlps xmm10, [rsi + r8*4 + 40]  ;# G1r H1r 

	movlps xmm3, [rsi + r9*4 + 8]    ;# G2c H2c
	movlps xmm7, [rsi + r9*4 + 24]   ;# G2d H2d
	movlps xmm11, [rsi + r9*4 + 40]  ;# G2r H2r

    unpcklps xmm2, xmm3   ;# G1c G2c H1c H2c
    unpcklps xmm6, xmm7   ;# G1d G2d H1d H2d
    unpcklps xmm10, xmm11 ;# G1r G2r H1r H2r
    movhlps  xmm3, xmm2   ;# H1c H2c
    movhlps  xmm7, xmm6   ;# H1d H2d
    movhlps  xmm11, xmm10 ;# H1r H2r
    ;# table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps xmm12, [rsp + nb330_eps]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm12
    mulps  xmm11, xmm12
    mulps  xmm2, xmm12     ;# Geps
    mulps  xmm6, xmm12
    mulps  xmm10, xmm12
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm12
    mulps  xmm11, xmm12

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
    mulps  xmm5, xmm12
    mulps  xmm9, xmm12
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, [rsp + nb330_qq]   ;# VV*qq = vcoul
    mulps  xmm5, [rsp + nb330_c6]   ;# vnb6
    mulps  xmm9, [rsp + nb330_c12]   ;# vnb12
    mulps  xmm3, [rsp + nb330_qq]    ;# FF*qq = fij
    mulps  xmm7, [rsp + nb330_c6]   ;# fijD
    mulps  xmm11, [rsp + nb330_c12]   ;#fijR

    ;# accumulate vctot
    addps  xmm1, [rsp + nb330_vctot]
    movlps [rsp + nb330_vctot], xmm1

    ;# accumulate Vvdwtot
    addps  xmm5, [rsp + nb330_Vvdwtot]
    addps  xmm5, xmm9
    movlps [rsp + nb330_Vvdwtot], xmm5

    xorps  xmm9, xmm9
    
    addps  xmm3, xmm7
    addps  xmm3, xmm11
    mulps  xmm3, xmm15
    mulps  xmm3, [rsp + nb330_tsc]  ;# fscal
    
    subps  xmm9, xmm3
    movaps xmm10, xmm9
    movaps xmm11, xmm9
    
    movaps xmm12, [rsp + nb330_fix]
    movaps xmm13, [rsp + nb330_fiy]
    movaps xmm14, [rsp + nb330_fiz]
    
    mulps  xmm9,  [rsp + nb330_dx]
    mulps  xmm10, [rsp + nb330_dy]
    mulps  xmm11, [rsp + nb330_dz]

    ;# accumulate i forces
    addps xmm12, xmm9
    addps xmm13, xmm10
    addps xmm14, xmm11
    movlps [rsp + nb330_fix], xmm12
    movlps [rsp + nb330_fiy], xmm13
    movlps [rsp + nb330_fiz], xmm14
    
	mov rsi, [rbp + nb330_faction]
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
	
.nb330_checksingle:				
	mov   edx, [rsp + nb330_innerk]
	and   edx, 1
	jnz    .nb330_dosingle
	jmp    .nb330_updateouterdata
.nb330_dosingle:
	mov rdi, [rbp + nb330_pos]
	mov   rcx, [rsp + nb330_innerjjnr]
	mov   eax, [rcx]	

	mov rsi, [rbp + nb330_charge]
	movss xmm0, [rsi + rax*4]

	mulss xmm0, [rsp + nb330_iq] 
    movaps [rsp + nb330_qq], xmm0

    ;# vdw parameters
	mov rsi, [rbp + nb330_type]
	mov r12d, [rsi + rax*4]
	shl r12d, 1	
    mov edi, [rsp + nb330_ntia]
	add r12d, edi

	mov rsi, [rbp + nb330_vdwparam]
	movss xmm0, [rsi + r12*4]
	movss xmm3, [rsi + r12*4 + 4]

    movaps [rsp + nb330_c6], xmm0
    movaps [rsp + nb330_c12], xmm3
    
	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	mov rdi, [rbp + nb330_pos]
	;# load coordinates
	movss xmm1, [rdi + rax*4]	
	movss xmm2, [rdi + rax*4 + 4]
	movss xmm5, [rdi + rax*4 + 8]

	;# calc dr  
	subss xmm1, [rsp + nb330_ix]
	subss xmm2, [rsp + nb330_iy]
	subss xmm5, [rsp + nb330_iz]

	;# store dr
    movaps [rsp + nb330_dx], xmm1
    movaps [rsp + nb330_dy], xmm2
    movaps [rsp + nb330_dz], xmm5
    
	;# square it 
	mulss xmm1,xmm1
	mulss xmm2,xmm2
	mulss xmm5,xmm5
	addss xmm1, xmm2
	addss xmm1, xmm5

	;# rsq in xmm1
    
    ;# calculate rinv=1/sqrt(rsq)
	rsqrtss xmm5, xmm1
	movaps xmm2, xmm5
	mulss xmm5, xmm5
	movaps xmm4, [rsp + nb330_three]
	mulss xmm5, xmm1	;# rsq*lu*lu 	
    subss xmm4, xmm5	;# 30-rsq*lu*lu 
	mulss xmm4, xmm2	
	mulss xmm4, [rsp + nb330_half]	
	movaps xmm15, xmm4
	mulss  xmm1, xmm4	
    ;# xmm15=rinv
    ;# xmm1=r

    mulss xmm1, [rsp + nb330_tsc] ;# rtab

    ;# truncate and convert to integers
    cvttss2si r8d, xmm1
    
    ;# convert back to float
    cvtsi2ss  xmm4, r8d
    
    ;# multiply by 4
    shl       r8d, 2

    ;# calculate eps
    subss     xmm1, xmm4

    ;# mult. by 3
	lea   r8, [r8 + r8*2]  

    movaps [rsp + nb330_eps], xmm1
    ;# xmm15=rinv
    
	mov rsi, [rbp + nb330_VFtab]
    ;# load Coulomb and LJ table data in parallel
    movss  xmm0,  [rsi + r8*4]
    movss  xmm1,  [rsi + r8*4 + 4]
    movss  xmm2,  [rsi + r8*4 + 8]
    movss  xmm3,  [rsi + r8*4 + 12]
    movss  xmm4,  [rsi + r8*4 + 16]
    movss  xmm5,  [rsi + r8*4 + 20]
    movss  xmm6,  [rsi + r8*4 + 24]
    movss  xmm7,  [rsi + r8*4 + 28]
    movss  xmm8,  [rsi + r8*4 + 32]
    movss  xmm9,  [rsi + r8*4 + 36]
    movss  xmm10, [rsi + r8*4 + 40]
    movss  xmm11, [rsi + r8*4 + 44]
    ;# table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movaps xmm12, [rsp + nb330_eps]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm12
    mulps  xmm11, xmm12
    mulps  xmm2, xmm12     ;# Geps
    mulps  xmm6, xmm12
    mulps  xmm10, xmm12
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm12
    mulps  xmm11, xmm12

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
    mulps  xmm1, [rsp + nb330_eps]   ;# eps*Fp
    mulps  xmm5, [rsp + nb330_eps]
    mulps  xmm9, [rsp + nb330_eps]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, [rsp + nb330_qq]   ;# VV*qq = vcoul
    mulps  xmm5, [rsp + nb330_c6]   ;# vnb6
    mulps  xmm9, [rsp + nb330_c12]   ;# vnb12
    mulps  xmm3, [rsp + nb330_qq]    ;# FF*qq = fij
    mulps  xmm7, [rsp + nb330_c6]   ;# fijD
    mulps  xmm11, [rsp + nb330_c12]   ;#fijR

    ;# accumulate vctot
    addps  xmm1, [rsp + nb330_vctot]
    movaps [rsp + nb330_vctot], xmm1

    ;# accumulate Vvdwtot
    addps  xmm5, [rsp + nb330_Vvdwtot]
    addps  xmm5, xmm9
    movaps [rsp + nb330_Vvdwtot], xmm5

    xorps  xmm9, xmm9
    
    addps  xmm3, xmm7
    addps  xmm3, xmm11
    mulss  xmm3, xmm15
    mulps  xmm3, [rsp + nb330_tsc]  ;# fscal
    
    subps  xmm9, xmm3
    movaps xmm10, xmm9
    movaps xmm11, xmm9
    
    movaps xmm12, [rsp + nb330_fix]
    movaps xmm13, [rsp + nb330_fiy]
    movaps xmm14, [rsp + nb330_fiz]
    
    mulss  xmm9,  [rsp + nb330_dx]
    mulss  xmm10, [rsp + nb330_dy]
    mulss  xmm11, [rsp + nb330_dz]
    
    ;# accumulate i forces
    addss xmm12, xmm9
    addss xmm13, xmm10
    addss xmm14, xmm11
    movaps [rsp + nb330_fix], xmm12
    movaps [rsp + nb330_fiy], xmm13
    movaps [rsp + nb330_fiz], xmm14
    
	mov rsi, [rbp + nb330_faction]
    ;# add to j forces
    addss  xmm9,  [rsi + rax*4]
    addss  xmm10, [rsi + rax*4 + 4]
    addss  xmm11, [rsi + rax*4 + 8]
    movss  [rsi + rax*4],     xmm9
    movss  [rsi + rax*4 + 4], xmm10
    movss  [rsi + rax*4 + 8], xmm11
    
.nb330_updateouterdata:
	mov   ecx, [rsp + nb330_ii3]
	mov   rdi, [rbp + nb330_faction]
	mov   rsi, [rbp + nb330_fshift]
	mov   edx, [rsp + nb330_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb330_fix]
	movaps xmm1, [rsp + nb330_fiy]
	movaps xmm2, [rsp + nb330_fiz]

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
	mov esi, [rsp + nb330_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb330_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb330_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb330_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb330_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb330_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb330_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb330_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb330_n], esi
        jmp .nb330_outer
.nb330_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb330_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb330_end
        ;# non-zero, do one more workunit
        jmp   .nb330_threadloop
.nb330_end:

	mov eax, [rsp + nb330_nouter]
	mov ebx, [rsp + nb330_ninner]
	mov rcx, [rbp + nb330_outeriter]
	mov rdx, [rbp + nb330_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 432
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

	


.globl nb_kernel330nf_x86_64_sse
.globl _nb_kernel330nf_x86_64_sse
nb_kernel330nf_x86_64_sse:	
_nb_kernel330nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb330nf_fshift,         16
.equiv          nb330nf_gid,            24
.equiv          nb330nf_pos,            32
.equiv          nb330nf_faction,        40
.equiv          nb330nf_charge,         48
.equiv          nb330nf_p_facel,        56
.equiv          nb330nf_argkrf,         64
.equiv          nb330nf_argcrf,         72
.equiv          nb330nf_Vc,             80
.equiv          nb330nf_type,           88
.equiv          nb330nf_p_ntype,        96
.equiv          nb330nf_vdwparam,       104
.equiv          nb330nf_Vvdw,           112
.equiv          nb330nf_p_tabscale,     120
.equiv          nb330nf_VFtab,          128
.equiv          nb330nf_invsqrta,       136
.equiv          nb330nf_dvda,           144
.equiv          nb330nf_p_gbtabscale,   152
.equiv          nb330nf_GBtab,          160
.equiv          nb330nf_p_nthreads,     168
.equiv          nb330nf_count,          176
.equiv          nb330nf_mtx,            184
.equiv          nb330nf_outeriter,      192
.equiv          nb330nf_inneriter,      200
.equiv          nb330nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb330nf_ix,             0
.equiv          nb330nf_iy,             16
.equiv          nb330nf_iz,             32
.equiv          nb330nf_iq,             48
.equiv          nb330nf_tsc,            64
.equiv          nb330nf_qq,             80
.equiv          nb330nf_c6,             96
.equiv          nb330nf_c12,            112
.equiv          nb330nf_vctot,          128
.equiv          nb330nf_Vvdwtot,        144
.equiv          nb330nf_half,           160
.equiv          nb330nf_three,          176
.equiv          nb330nf_nri,            192
.equiv          nb330nf_iinr,           200
.equiv          nb330nf_jindex,         208
.equiv          nb330nf_jjnr,           216
.equiv          nb330nf_shift,          224
.equiv          nb330nf_shiftvec,       232
.equiv          nb330nf_facel,          240
.equiv          nb330nf_innerjjnr,      248
.equiv          nb330nf_is3,            256
.equiv          nb330nf_ii3,            260
.equiv          nb330nf_ntia,           264
.equiv          nb330nf_innerk,         268
.equiv          nb330nf_n,              272
.equiv          nb330nf_nn1,            276
.equiv          nb330nf_ntype,          280
.equiv          nb330nf_nouter,         284
.equiv          nb330nf_ninner,         288

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
	mov [rsp + nb330nf_nouter], eax
	mov [rsp + nb330nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb330nf_nri], edi
	mov [rsp + nb330nf_iinr], rsi
	mov [rsp + nb330nf_jindex], rdx
	mov [rsp + nb330nf_jjnr], rcx
	mov [rsp + nb330nf_shift], r8
	mov [rsp + nb330nf_shiftvec], r9
	mov rdi, [rbp + nb330nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb330nf_ntype], edi
	mov rsi, [rbp + nb330nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb330nf_facel], xmm0


	mov rax, [rbp + nb330nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb330nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb330nf_half], eax
	movss xmm1, [rsp + nb330nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb330nf_half],  xmm1
	movaps [rsp + nb330nf_three],  xmm3

.nb330nf_threadloop:
        mov   rsi, [rbp + nb330nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb330nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb330nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb330nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb330nf_n], eax
        mov [rsp + nb330nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb330nf_outerstart
        jmp .nb330nf_end

.nb330nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb330nf_nouter]
	mov [rsp + nb330nf_nouter], ebx

.nb330nf_outer:
	mov   rax, [rsp + nb330nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb330nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb330nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb330nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb330nf_charge]
	movss xmm3, [rdx + rbx*4]	
	mulss xmm3, [rsp + nb330nf_facel]
	shufps xmm3, xmm3, 0

    	mov   rdx, [rbp + nb330nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb330nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb330nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb330nf_pos]    ;# rax = base of pos[]  

	addss xmm0, [rax + rbx*4]
	addss xmm1, [rax + rbx*4 + 4]
	addss xmm2, [rax + rbx*4 + 8]

	movaps [rsp + nb330nf_iq], xmm3
	
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0

	movaps [rsp + nb330nf_ix], xmm0
	movaps [rsp + nb330nf_iy], xmm1
	movaps [rsp + nb330nf_iz], xmm2

	mov   [rsp + nb330nf_ii3], ebx
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb330nf_vctot], xmm4
	movaps [rsp + nb330nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb330nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb330nf_pos]	
	mov   rax, [rsp + nb330nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb330nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb330nf_ninner]
	mov   [rsp + nb330nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb330nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb330nf_unroll_loop
	jmp   .nb330nf_finish_inner
.nb330nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb330nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	add qword ptr [rsp + nb330nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb330nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	movaps xmm2, [rsp + nb330nf_iq]
	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	mulps  xmm3, xmm2
	movd  mm2, ecx
	movd  mm3, edx

	movaps [rsp + nb330nf_qq], xmm3
	
	mov rsi, [rbp + nb330nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]
	mov rsi, [rbp + nb330nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb330nf_ntia]
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

	movaps [rsp + nb330nf_c6], xmm4
	movaps [rsp + nb330nf_c12], xmm6
	
	mov rsi, [rbp + nb330nf_pos]       ;# base of pos[] 

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
	movaps xmm4, [rsp + nb330nf_ix]
	movaps xmm5, [rsp + nb330nf_iy]
	movaps xmm6, [rsp + nb330nf_iz]

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
	movaps xmm1, [rsp + nb330nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb330nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb330nf_tsc]

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

	mov  rsi, [rbp + nb330nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ecx, mm7
	psrlq mm7, 32
	movd ebx, mm6
	movd edx, mm7

	lea   rax, [rax + rax*2]
	lea   rbx, [rbx + rbx*2]
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]
	
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
	movaps xmm3, [rsp + nb330nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addps  xmm5, [rsp + nb330nf_vctot]
	movaps [rsp + nb330nf_vctot], xmm5 

	;# dispersion 
	movlps xmm5, [rsi + rax*4 + 16]
	movlps xmm7, [rsi + rcx*4 + 16]
	movhps xmm5, [rsi + rbx*4 + 16]
	movhps xmm7, [rsi + rdx*4 + 16] ;# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101
	
	movlps xmm7, [rsi + rax*4 + 24]
	movlps xmm3, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rbx*4 + 24]
	movhps xmm3, [rsi + rdx*4 + 24] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb330nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 
	;# put scalar force on stack 
	addps  xmm5, [rsp + nb330nf_Vvdwtot]
	movaps [rsp + nb330nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rax*4 + 32]
	movlps xmm7, [rsi + rcx*4 + 32]
	movhps xmm5, [rsi + rbx*4 + 32]
	movhps xmm7, [rsi + rdx*4 + 32] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rax*4 + 40]
	movlps xmm3, [rsi + rcx*4 + 40]
	movhps xmm7, [rsi + rbx*4 + 40]
	movhps xmm3, [rsi + rdx*4 + 40] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb330nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 
	addps  xmm5, [rsp + nb330nf_Vvdwtot]
	movaps [rsp + nb330nf_Vvdwtot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb330nf_innerk],  4
	jl    .nb330nf_finish_inner
	jmp   .nb330nf_unroll_loop
.nb330nf_finish_inner:
	;# check if at least two particles remain 
	add dword ptr [rsp + nb330nf_innerk],  4
	mov   edx, [rsp + nb330nf_innerk]
	and   edx, 2
	jnz   .nb330nf_dopair
	jmp   .nb330nf_checksingle
.nb330nf_dopair:	
	mov rsi, [rbp + nb330nf_charge]

    	mov   rcx, [rsp + nb330nf_innerjjnr]
	
	mov   eax, [rcx]	
	mov   ebx, [rcx + 4]              
	add qword ptr [rsp + nb330nf_innerjjnr],  8	
	xorps xmm7, xmm7
	movss xmm3, [rsi + rax*4]		
	movss xmm6, [rsi + rbx*4]
	shufps xmm3, xmm6, 0 
	shufps xmm3, xmm3, 8 ;# 00001000 ;# xmm3(0,1) has the charges 

	mulps  xmm3, [rsp + nb330nf_iq]
	movlhps xmm3, xmm7
	movaps [rsp + nb330nf_qq], xmm3

	mov rsi, [rbp + nb330nf_type]
	mov   ecx, eax
	mov   edx, ebx
	mov ecx, [rsi + rcx*4]
	mov edx, [rsi + rdx*4]	
	mov rsi, [rbp + nb330nf_vdwparam]
	shl ecx, 1	
	shl edx, 1	
	mov edi, [rsp + nb330nf_ntia]
	add ecx, edi
	add edx, edi
	movlps xmm6, [rsi + rcx*4]
	movhps xmm6, [rsi + rdx*4]
	mov rdi, [rbp + nb330nf_pos]	
	
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 8 ;# 00001000 	
	shufps xmm6, xmm6, 13 ;# 00001101
	movlhps xmm4, xmm7
	movlhps xmm6, xmm7
	
	movaps [rsp + nb330nf_c6], xmm4
	movaps [rsp + nb330nf_c12], xmm6	
	
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
	
	movaps xmm4, [rsp + nb330nf_ix]
	movaps xmm5, [rsp + nb330nf_iy]
	movaps xmm6, [rsp + nb330nf_iz]

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
	movaps xmm1, [rsp + nb330nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb330nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb330nf_tsc]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb330nf_VFtab]
	movd ecx, mm6
	psrlq mm6, 32
	movd edx, mm6
	lea   rcx, [rcx + rcx*2]
	lea   rdx, [rdx + rdx*2]

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
	movaps xmm3, [rsp + nb330nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV 
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addps  xmm5, [rsp + nb330nf_vctot]
	movaps [rsp + nb330nf_vctot], xmm5 

	;# dispersion 
	movlps xmm5, [rsi + rcx*4 + 16]
	movhps xmm5, [rsi + rdx*4 + 16];# got half dispersion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm4, 136  ;# 10001000
	shufps xmm5, xmm5, 221  ;# 11011101
	
	movlps xmm7, [rsi + rcx*4 + 24]
	movhps xmm7, [rsi + rdx*4 + 24] ;# other half of dispersion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 136  ;# 10001000
	shufps xmm7, xmm7, 221  ;# 11011101
	;# dispersion table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb330nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 
	;# put scalar force on stack 
	addps  xmm5, [rsp + nb330nf_Vvdwtot]
	movaps [rsp + nb330nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm5, [rsi + rcx*4 + 32]
	movhps xmm5, [rsi + rdx*4 + 32] ;# got half repulsion table 
	movaps xmm4, xmm5
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm5, xmm7, 221  ;# 11011101

	movlps xmm7, [rsi + rcx*4 + 40]
	movhps xmm7, [rsi + rdx*4 + 40] ;# other half of repulsion table 
	movaps xmm6, xmm7
	shufps xmm6, xmm3, 136  ;# 10001000
	shufps xmm7, xmm3, 221  ;# 11011101
	;# table ready, in xmm4-xmm7 	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb330nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 
	addps  xmm5, [rsp + nb330nf_Vvdwtot]
	movaps [rsp + nb330nf_Vvdwtot], xmm5
	
.nb330nf_checksingle:				
	mov   edx, [rsp + nb330nf_innerk]
	and   edx, 1
	jnz    .nb330nf_dosingle
	jmp    .nb330nf_updateouterdata
.nb330nf_dosingle:
	mov rsi, [rbp + nb330nf_charge]
	mov rdi, [rbp + nb330nf_pos]
	mov   rcx, [rsp + nb330nf_innerjjnr]
	mov   eax, [rcx]	
	xorps  xmm6, xmm6
	movss xmm6, [rsi + rax*4]	;# xmm6(0) has the charge 	
	mulps  xmm6, [rsp + nb330nf_iq]
	movaps [rsp + nb330nf_qq], xmm6

	mov rsi, [rbp + nb330nf_type]
	mov ecx, eax
	mov ecx, [rsi + rcx*4]	
	mov rsi, [rbp + nb330nf_vdwparam]
	shl ecx, 1
	add ecx, [rsp + nb330nf_ntia]
	movlps xmm6, [rsi + rcx*4]
	movaps xmm4, xmm6
	shufps xmm4, xmm4, 252  ;# 11111100	
	shufps xmm6, xmm6, 253  ;# 11111101	
	
	movaps [rsp + nb330nf_c6], xmm4
	movaps [rsp + nb330nf_c12], xmm6	
	
	lea   rax, [rax + rax*2]
	
	;# move coordinates to xmm0-xmm2 
	movss xmm0, [rdi + rax*4]	
	movss xmm1, [rdi + rax*4 + 4]	
	movss xmm2, [rdi + rax*4 + 8]	 
	
	movaps xmm4, [rsp + nb330nf_ix]
	movaps xmm5, [rsp + nb330nf_iy]
	movaps xmm6, [rsp + nb330nf_iz]

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
	movaps xmm1, [rsp + nb330nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb330nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	mulps xmm4, xmm0	;# xmm4=r 
	mulps xmm4, [rsp + nb330nf_tsc]

	cvttps2pi mm6, xmm4     ;# mm6 contain lu indices 
	cvtpi2ps xmm6, mm6
	subps xmm4, xmm6	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1	
	mulps  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 2

	mov  rsi, [rbp + nb330nf_VFtab]
	movd ebx, mm6
	
	lea  rbx, [rbx + rbx*2]
	
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
	movaps xmm3, [rsp + nb330nf_qq]
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
	mulps  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addss  xmm5, [rsp + nb330nf_vctot]
	movss [rsp + nb330nf_vctot], xmm5 

	;# dispersion 
	movlps xmm4, [rsi + rbx*4 + 16]
	movlps xmm6, [rsi + rbx*4 + 24]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1
	shufps xmm7, xmm7, 1
	;# table ready in xmm4-xmm7 
	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 

	movaps xmm4, [rsp + nb330nf_c6]
	mulps  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack  
	addss  xmm5, [rsp + nb330nf_Vvdwtot]
	movss [rsp + nb330nf_Vvdwtot], xmm5

	;# repulsion 
	movlps xmm4, [rsi + rbx*4 + 32]
	movlps xmm6, [rsi + rbx*4 + 40]
	movaps xmm5, xmm4
	movaps xmm7, xmm6
	shufps xmm5, xmm5, 1
	shufps xmm7, xmm7, 1
	;# table ready in xmm4-xmm7 
	
	mulps  xmm6, xmm1	;# xmm6=Geps 
	mulps  xmm7, xmm2	;# xmm7=Heps2 
	addps  xmm5, xmm6
	addps  xmm5, xmm7	;# xmm5=Fp 	
	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
	addps  xmm5, xmm4 ;# xmm5=VV 
 	
	movaps xmm4, [rsp + nb330nf_c12]
	mulps  xmm5, xmm4 ;# Vvdw12 
	addss  xmm5, [rsp + nb330nf_Vvdwtot]
	movss [rsp + nb330nf_Vvdwtot], xmm5
	
.nb330nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb330nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb330nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb330nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb330nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
	;# accumulate total lj energy and update it 
	movaps xmm7, [rsp + nb330nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb330nf_Vvdw]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb330nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb330nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb330nf_n], esi
        jmp .nb330nf_outer
.nb330nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb330nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb330nf_end
        ;# non-zero, do one more workunit
        jmp   .nb330nf_threadloop
.nb330nf_end:

	mov eax, [rsp + nb330nf_nouter]
	mov ebx, [rsp + nb330nf_ninner]
	mov rcx, [rbp + nb330nf_outeriter]
	mov rdx, [rbp + nb330nf_inneriter]
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
