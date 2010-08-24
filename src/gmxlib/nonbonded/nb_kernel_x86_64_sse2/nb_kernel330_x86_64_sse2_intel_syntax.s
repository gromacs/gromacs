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

.globl nb_kernel330_x86_64_sse2
.globl _nb_kernel330_x86_64_sse2
nb_kernel330_x86_64_sse2:	
_nb_kernel330_x86_64_sse2:	
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
	;# bottom of stack is cache-aligned for sse2 use 
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

	emms
	sub rsp, 432		;# local variable stack space (n*16+8)

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
	movsd xmm0, [rsi]
	movsd [rsp + nb330_facel], xmm0

	mov rax, [rbp + nb330_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb330_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb330_half], eax
	mov [rsp + nb330_half+4], ebx
	movsd xmm1, [rsp + nb330_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb330_half], xmm1
	movapd [rsp + nb330_three], xmm3

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
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb330_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb330_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb330_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb330_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb330_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb330_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb330_ntype]
    	shl   edx, 1
    	mov   [rsp + nb330_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb330_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb330_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb330_ix], xmm0
	movapd [rsp + nb330_iy], xmm1
	movapd [rsp + nb330_iz], xmm2

	mov   [rsp + nb330_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb330_vctot], xmm4
	movapd [rsp + nb330_Vvdwtot], xmm4
	movapd [rsp + nb330_fix], xmm4
	movapd [rsp + nb330_fiy], xmm4
	movapd [rsp + nb330_fiz], xmm4
	
	mov   rax, [rsp + nb330_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb330_pos]
	mov   rdi, [rbp + nb330_faction]	
	mov   rax, [rsp + nb330_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb330_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb330_ninner]
	mov   [rsp + nb330_ninner], ecx
	add   edx, 0
	mov   [rsp + nb330_innerk], edx    ;# number of innerloop atoms 
	jge   .nb330_unroll_loop
	jmp   .nb330_checksingle
.nb330_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb330_innerjjnr]   ;# pointer to jjnr[k] 
	mov   r12d, [rdx]
	mov   r13d, [rdx + 4]
	add qword ptr [rsp + nb330_innerjjnr], 8	;# advance pointer (unrolled 2) 

	
	mov rsi, [rbp + nb330_pos]		;# base of pos[] 

	lea   rax, [r12 + r12*2]     ;# replace jnr with j3 
	lea   rbx, [r13 + r13*2]	

	;# move two coordinates to xmm4-xmm6 
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
	movhpd xmm4, [rsi + rbx*8]
	movhpd xmm5, [rsi + rbx*8 + 8]
	movhpd xmm6, [rsi + rbx*8 + 16]		
	
	;# calc dr 
	subpd xmm4, [rsp + nb330_ix]
	subpd xmm5, [rsp + nb330_iy]
	subpd xmm6, [rsp + nb330_iz]

	;# store dr 
	movapd [rsp + nb330_dx], xmm4
	movapd [rsp + nb330_dy], xmm5
	movapd [rsp + nb330_dz], xmm6

	mov rsi, [rbp + nb330_charge]    ;# base of charge[] 

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 
	movlpd xmm3, [rsi + r12*8]

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 
	movhpd xmm3, [rsi + r13*8]
	mov rsi, [rbp + nb330_type]

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb330_three]
	mulpd  xmm3, [rsp + nb330_iq]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb330_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 
	movapd [rsp + nb330_qq], xmm3	

	mov r8d, [rsi + r12*4]
	mov r9d, [rsi + r13*4]

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb330_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb330_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb330_tsc]
    movapd xmm15, xmm0  ;# copy of rinv
	shl r8d, 1
	shl r9d, 1
	mov edi, [rsp + nb330_ntia]
    
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	add r8d, edi
	add r9d, edi
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd [rsp + nb330_eps], xmm4
	
	pslld mm6, 2		;# idx *= 4 

	mov  rsi, [rbp + nb330_VFtab]
	movd r10d, mm6
	psrlq mm6, 32
	movd r11d, mm6		;# indices in r10,r11
     
	lea   r10, [r10 + r10*2]	;# idx*=3 (12 total now) 
	lea   r11, [r11 + r11*2]	;# idx*=3 (12 total now) 

    ;# load Coulomb and LJ table data in parallel
    movapd xmm0,  [rsi + r10*8]        ;# Y1c F1c
    movapd xmm12, [rsi + r11*8]        ;# Y2c F2c
    movapd xmm4,  [rsi + r10*8 + 32]   ;# Y1d F1d
    movapd xmm13, [rsi + r11*8 + 32]   ;# Y2d F2d
    movapd xmm8,  [rsi + r10*8 + 64]   ;# Y1r F1r
    movapd xmm14, [rsi + r11*8 + 64]   ;# Y2r F2r
	movapd xmm1, xmm0
	movapd xmm5, xmm4
	movapd xmm9, xmm8
	unpcklpd xmm0, xmm12	;# Y1c Y2c 
	unpckhpd xmm1, xmm12	;# F1c F2c 
	unpcklpd xmm4, xmm13	;# Y1d Y2d 
	unpckhpd xmm5, xmm13	;# F1d F2d 
	unpcklpd xmm8, xmm14	;# Y1r Y2r 
	unpckhpd xmm9, xmm14	;# F1r F2r 
    
    movapd xmm2,  [rsi + r10*8 + 16]   ;# G1c H1c
    movapd xmm12, [rsi + r11*8 + 16]   ;# G2c H2c
    movapd xmm6,  [rsi + r10*8 + 48]   ;# G1d H1d
    movapd xmm13, [rsi + r11*8 + 48]   ;# G2d H2d
    movapd xmm10, [rsi + r10*8 + 80]   ;# G1r H1r
    movapd xmm14, [rsi + r11*8 + 80]   ;# G2r H2r
	movapd xmm3, xmm2
	movapd xmm7, xmm6
	movapd xmm11, xmm10
	unpcklpd xmm2, xmm12	;# G1c G2c 
	unpckhpd xmm3, xmm12	;# H1c H2c 
	unpcklpd xmm6, xmm13	;# G1d G2d 
	unpckhpd xmm7, xmm13	;# H1d H2d 
	unpcklpd xmm10, xmm14	;# G1r G2r 
	unpckhpd xmm11, xmm14	;# H1r H2r 
    ;# table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11
	mov rsi, [rbp + nb330_vdwparam]

    movapd xmm12, [rsp + nb330_eps]
    
    mulpd  xmm3, xmm12   ;# Heps
    mulpd  xmm7, xmm12
    mulpd  xmm11, xmm12
    mulpd  xmm2, xmm12     ;# Geps
    mulpd  xmm6, xmm12
    mulpd  xmm10, xmm12
    mulpd  xmm3, xmm12   ;# Heps2
    mulpd  xmm7, xmm12
    mulpd  xmm11, xmm12

    movlpd xmm13, [rsi + r8*8]   ;# c6
    movlpd xmm14, [rsi + r8*8 + 8]  ;# c12
    
    addpd  xmm1, xmm2   ;# F+Geps
    addpd  xmm5, xmm6
    addpd  xmm9, xmm10 
    addpd  xmm1, xmm3   ;# F+Geps+Heps2 = Fp
    addpd  xmm5, xmm7
    addpd  xmm9, xmm11 
    addpd  xmm3, xmm3    ;# 2*Heps2
    addpd  xmm7, xmm7
    addpd  xmm11, xmm11
    addpd  xmm3, xmm2    ;# 2*Heps2+Geps
    addpd  xmm7, xmm6  
    addpd  xmm11, xmm10
    movhpd xmm13, [rsi + r9*8]   ;# c6
    movhpd xmm14, [rsi + r9*8 + 8]  ;# c12
    
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, xmm12   ;# eps*Fp
    mulpd  xmm5, xmm12
    mulpd  xmm9, xmm12
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, [rsp + nb330_qq]   ;# VV*qq = vcoul
    mulpd  xmm5, xmm13   ;# vnb6
    mulpd  xmm9, xmm14   ;# vnb12
    mulpd  xmm3, [rsp + nb330_qq]    ;# FF*qq = fij
    mulpd  xmm7, xmm13   ;# fijD
    mulpd  xmm11, xmm14   ;#fijR

    ;# accumulate vctot
    addpd  xmm1, [rsp + nb330_vctot]
    movapd [rsp + nb330_vctot], xmm1

    ;# accumulate Vvdwtot
    addpd  xmm5, [rsp + nb330_Vvdwtot]
    addpd  xmm5, xmm9
    movapd [rsp + nb330_Vvdwtot], xmm5

    xorpd  xmm9, xmm9
    
    ;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb330_faction]
	movlpd xmm4, [rdi + rax*8]
	movlpd xmm5, [rdi + rax*8 + 8]
	movlpd xmm6, [rdi + rax*8 + 16]
	movhpd xmm4, [rdi + rbx*8]
	movhpd xmm5, [rdi + rbx*8 + 8]
	movhpd xmm6, [rdi + rbx*8 + 16]

    addpd  xmm3, xmm7
    addpd  xmm3, xmm11
    mulpd  xmm3, xmm15
    mulpd  xmm3, [rsp + nb330_tsc]  ;# fscal
    
    subpd  xmm9, xmm3
    movapd xmm10, xmm9
    movapd xmm11, xmm9
    
    movapd xmm12, [rsp + nb330_fix]
    movapd xmm13, [rsp + nb330_fiy]
    movapd xmm14, [rsp + nb330_fiz]
    
    mulpd  xmm9,  [rsp + nb330_dx]
    mulpd  xmm10, [rsp + nb330_dy]
    mulpd  xmm11, [rsp + nb330_dz]

	addpd xmm4, xmm9
	addpd xmm5, xmm10
	addpd xmm6, xmm11

    ;# accumulate i forces
    addpd xmm12, xmm9
    addpd xmm13, xmm10
    addpd xmm14, xmm11
	movlpd [rdi + rax*8], xmm4
	movlpd [rdi + rax*8 + 8], xmm5
	movlpd [rdi + rax*8 + 16], xmm6

    movapd [rsp + nb330_fix], xmm12
    movapd [rsp + nb330_fiy], xmm13
    movapd [rsp + nb330_fiz], xmm14
 

	movhpd [rdi + rbx*8], xmm4
	movhpd [rdi + rbx*8 + 8], xmm5
	movhpd [rdi + rbx*8 + 16], xmm6
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb330_innerk],  2
	jl    .nb330_checksingle
	jmp   .nb330_unroll_loop
.nb330_checksingle:
	mov   edx, [rsp + nb330_innerk]
	and   edx, 1
	jnz    .nb330_dosingle
	jmp    .nb330_updateouterdata
.nb330_dosingle:
	mov rsi, [rbp + nb330_charge]
	mov rdi, [rbp + nb330_pos]
	mov   rcx, [rsp + nb330_innerjjnr]
	mov   eax, [rcx]	

	movsd xmm3, [rsi + rax*8]
	mulsd  xmm3, [rsp + nb330_iq]
	movapd [rsp + nb330_qq], xmm3	

	mov rsi, [rbp + nb330_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb330_vdwparam]
	shl r8d, 1
	mov edi, [rsp + nb330_ntia]
	add r8d, edi

	movsd xmm4, [rsi + r8*8]
	movsd xmm6, [rsi + r8*8 + 8]
	movapd [rsp + nb330_c6], xmm4
	movapd [rsp + nb330_c12], xmm6
	
	mov rsi, [rbp + nb330_pos]		;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinate to xmm4-xmm6 
	movsd xmm4, [rsi + rax*8]
	movsd xmm5, [rsi + rax*8 + 8]
	movsd xmm6, [rsi + rax*8 + 16]

	mov    rdi, [rbp + nb330_faction]
	
	;# calc dr 
	subsd xmm4, [rsp + nb330_ix]
	subsd xmm5, [rsp + nb330_iy]
	subsd xmm6, [rsp + nb330_iz]

	;# store dr 
	movapd [rsp + nb330_dx], xmm4
	movapd [rsp + nb330_dy], xmm5
	movapd [rsp + nb330_dz], xmm6

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
	movapd xmm1, [rsp + nb330_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb330_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb330_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb330_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb330_tsc]
    movapd xmm15, xmm0  ;# copy of rinv
    
	cvttsd2si r10d, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, r10d
	subsd xmm4, xmm5
	movapd [rsp + nb330_eps], xmm4
	
	shl   r10d, 2		;# idx *= 4 

	mov  rsi, [rbp + nb330_VFtab]
     
	lea   r10, [r10 + r10*2]	;# idx*=3 (12 total now) 

    ;# load Coulomb and LJ table data in parallel    
    movapd xmm0,  [rsi + r10*8]        ;# Y1c F1c
    movapd xmm4,  [rsi + r10*8 + 32]   ;# Y1d F1d
    movapd xmm8,  [rsi + r10*8 + 64]   ;# Y1r F1r
	movhlps xmm1, xmm0
	movhlps xmm5, xmm4
	movhlps xmm9, xmm8

    movapd xmm2,  [rsi + r10*8 + 16]   ;# G1c H1c
    movapd xmm6,  [rsi + r10*8 + 48]   ;# G1d H1d
    movapd xmm10, [rsi + r10*8 + 80]   ;# G1r H1r
	movhlps xmm3, xmm2
	movhlps xmm7, xmm6
	movhlps xmm11, xmm10
    ;# table data ready. Coul in xmm0-xmm3 , disp in xmm4-xmm7 , rep. in xmm8-xmm11

    movapd xmm12, [rsp + nb330_eps]
    
    mulsd  xmm3, xmm12   ;# Heps
    mulsd  xmm7, xmm12
    mulsd  xmm11, xmm12
    mulsd  xmm2, xmm12     ;# Geps
    mulsd  xmm6, xmm12
    mulsd  xmm10, xmm12
    mulsd  xmm3, xmm12   ;# Heps2
    mulsd  xmm7, xmm12
    mulsd  xmm11, xmm12

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
    mulsd  xmm1, [rsp + nb330_eps]   ;# eps*Fp
    mulsd  xmm5, [rsp + nb330_eps]
    mulsd  xmm9, [rsp + nb330_eps]
    addsd  xmm1, xmm0     ;# VV
    addsd  xmm5, xmm4
    addsd  xmm9, xmm8
    mulsd  xmm1, [rsp + nb330_qq]   ;# VV*qq = vcoul
    mulsd  xmm5, [rsp + nb330_c6]   ;# vnb6
    mulsd  xmm9, [rsp + nb330_c12]   ;# vnb12
    mulsd  xmm3, [rsp + nb330_qq]    ;# FF*qq = fij
    mulsd  xmm7, [rsp + nb330_c6]   ;# fijD
    mulsd  xmm11, [rsp + nb330_c12]   ;#fijR

    ;# accumulate vctot
    addsd  xmm1, [rsp + nb330_vctot]
    movsd [rsp + nb330_vctot], xmm1

    ;# accumulate Vvdwtot
    addsd  xmm5, [rsp + nb330_Vvdwtot]
    addsd  xmm5, xmm9
    movsd [rsp + nb330_Vvdwtot], xmm5

    xorpd  xmm9, xmm9
    
    addsd  xmm3, xmm7
    addsd  xmm3, xmm11
    mulsd  xmm3, xmm15
    mulsd  xmm3, [rsp + nb330_tsc]  ;# fscal
    
    subsd  xmm9, xmm3
    movapd xmm10, xmm9
    movapd xmm11, xmm9
    
    movapd xmm12, [rsp + nb330_fix]
    movapd xmm13, [rsp + nb330_fiy]
    movapd xmm14, [rsp + nb330_fiz]
    
    mulsd  xmm9,  [rsp + nb330_dx]
    mulsd  xmm10, [rsp + nb330_dy]
    mulsd  xmm11, [rsp + nb330_dz]

    ;# accumulate i forces
    addsd xmm12, xmm9
    addsd xmm13, xmm10
    addsd xmm14, xmm11
    movsd [rsp + nb330_fix], xmm12
    movsd [rsp + nb330_fiy], xmm13
    movsd [rsp + nb330_fiz], xmm14
 
	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb330_faction]
	addsd xmm9,   [rdi + rax*8]
	addsd xmm10,  [rdi + rax*8 + 8]
	addsd xmm11,  [rdi + rax*8 + 16]
	movsd [rdi + rax*8], xmm9
	movsd [rdi + rax*8 + 8], xmm10
	movsd [rdi + rax*8 + 16], xmm11
	
.nb330_updateouterdata:
	mov   ecx, [rsp + nb330_ii3]
	mov   rdi, [rbp + nb330_faction]
	mov   rsi, [rbp + nb330_fshift]
	mov   edx, [rsp + nb330_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb330_fix]
	movapd xmm1, [rsp + nb330_fiy]
	movapd xmm2, [rsp + nb330_fiz]

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
	mov esi, [rsp + nb330_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb330_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb330_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb330_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb330_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb330_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
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



.globl nb_kernel330nf_x86_64_sse2
.globl _nb_kernel330nf_x86_64_sse2
nb_kernel330nf_x86_64_sse2:	
_nb_kernel330nf_x86_64_sse2:	
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

	emms
	sub rsp, 304		;# local variable stack space (n*16+8)

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
	movsd xmm0, [rsi]
	movsd [rsp + nb330nf_facel], xmm0

	mov rax, [rbp + nb330nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb330nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb330nf_half], eax
	mov [rsp + nb330nf_half+4], ebx
	movsd xmm1, [rsp + nb330nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb330nf_half], xmm1
	movapd [rsp + nb330nf_three], xmm3

.nb330nf_threadloop:
        mov   rsi, [rbp + nb330nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb330nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
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
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb330nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb330nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb330nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb330nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb330nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb330nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb330nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb330nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb330nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb330nf_ix], xmm0
	movapd [rsp + nb330nf_iy], xmm1
	movapd [rsp + nb330nf_iz], xmm2

	mov   [rsp + nb330nf_ii3], ebx
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb330nf_vctot], xmm4
	movapd [rsp + nb330nf_Vvdwtot], xmm4
	
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
	sub   edx,  2
	add   ecx, [rsp + nb330nf_ninner]
	mov   [rsp + nb330nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb330nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb330nf_unroll_loop
	jmp   .nb330nf_checksingle
.nb330nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb330nf_innerjjnr]   ;# pointer to jjnr[k] 
	mov   eax, [rdx]
	mov   ebx, [rdx + 4]
	add qword ptr [rsp + nb330nf_innerjjnr], 8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb330nf_charge]    ;# base of charge[] 

	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	movapd xmm2, [rsp + nb330nf_iq]
	mulpd  xmm3, xmm2
	movapd [rsp + nb330nf_qq], xmm3	

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb330nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb330nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [rsp + nb330nf_ntia]
	add eax, edi
	add ebx, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movlpd xmm7, [rsi + rbx*8]	;# c6b
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + rbx*8 + 8]	;# c6b c12b 

	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb330nf_c6], xmm4
	movapd [rsp + nb330nf_c12], xmm6
	
	mov rsi, [rbp + nb330nf_pos]		;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	;# move nb330nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb330nf_ix]
	movapd xmm5, [rsp + nb330nf_iy]
	movapd xmm6, [rsp + nb330nf_iz]

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
	movapd xmm1, [rsp + nb330nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb330nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb330nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb330nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb330nf_tsc]

	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 

	mov  rsi, [rbp + nb330nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 
	lea   rbx, [rbx + rbx*2]	;# idx*=3 (12 total now) 

	;# Coulomb 
	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 16]	;# G2 H2 
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# coulomb table ready, in xmm4-xmm7  		
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb330nf_qq]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV 
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addpd  xmm5, [rsp + nb330nf_vctot]
	movapd [rsp + nb330nf_vctot], xmm5 

	;# Dispersion 
	movapd xmm4, [rsi + rax*8 + 32]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8 + 32]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 48]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 48]	;# G2 H2 
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

	mulpd  xmm5, [rsp + nb330nf_c6] ;# Vvdw6 

	addpd  xmm5, [rsp + nb330nf_Vvdwtot]
	movapd [rsp + nb330nf_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [rsi + rax*8 + 64]	;# Y1 F1 	
	movapd xmm3, [rsi + rbx*8 + 64]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movapd xmm6, [rsi + rax*8 + 80]	;# G1 H1 	
	movapd xmm3, [rsi + rbx*8 + 80]	;# G2 H2 
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

	mulpd  xmm5, [rsp + nb330nf_c12] ;# Vvdw12 
	
	addpd  xmm5, [rsp + nb330nf_Vvdwtot]
	movapd [rsp + nb330nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb330nf_innerk],  2
	jl    .nb330nf_checksingle
	jmp   .nb330nf_unroll_loop
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
	xorpd  xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]	;# xmm6(0) has the charge 	
	mulpd  xmm3, [rsp + nb330nf_iq]
	movapd [rsp + nb330nf_qq], xmm3

	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov rsi, [rbp + nb330nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb330nf_vdwparam]
	shl eax, 1
	mov edi, [rsp + nb330nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb330nf_c6], xmm4
	movapd [rsp + nb330nf_c12], xmm6
	
	mov rsi, [rbp + nb330nf_pos]		;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move nb330nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb330nf_ix]
	movapd xmm5, [rsp + nb330nf_iy]
	movapd xmm6, [rsp + nb330nf_iz]

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
	movapd xmm1, [rsp + nb330nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb330nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb330nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb330nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb330nf_tsc]

	movd mm0, eax	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb330nf_VFtab]
	lea   rax, [rax + rax*2]	;# idx*=3 (12 total now) 

	;# Coulomb 
	movapd xmm4, [rsi + rax*8]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 16]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb330nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	;# at this point mm5 contains vcoul 
	;# increment vcoul - then we can get rid of mm5 
	;# update vctot 
	addsd  xmm5, [rsp + nb330nf_vctot]
	movlpd [rsp + nb330nf_vctot], xmm5 

	;# Dispersion 
	movapd xmm4, [rsi + rax*8 + 32]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 48]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb330nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [rsp + nb330nf_c6]	 ;# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb330nf_Vvdwtot]
	movlpd [rsp + nb330nf_Vvdwtot], xmm5

	;# Repulsion 
	movapd xmm4, [rsi + rax*8 + 64]	;# Y1 F1 	
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + rax*8 + 80]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# Dispersion table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	movapd xmm3, [rsp + nb330nf_qq]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [rsp + nb330nf_c12] ;# Vvdw12 
	
	addsd  xmm5, [rsp + nb330nf_Vvdwtot]
	movlpd [rsp + nb330nf_Vvdwtot], xmm5
	
.nb330nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb330nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb330nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb330nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb330nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb330nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb330nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
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
