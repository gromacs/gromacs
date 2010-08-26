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



.globl nb_kernel030_x86_64_sse2
.globl _nb_kernel030_x86_64_sse2
nb_kernel030_x86_64_sse2:	
_nb_kernel030_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb030_fshift,           16
.equiv          nb030_gid,              24
.equiv          nb030_pos,              32
.equiv          nb030_faction,          40
.equiv          nb030_charge,           48
.equiv          nb030_p_facel,          56
.equiv          nb030_argkrf,           64
.equiv          nb030_argcrf,           72
.equiv          nb030_Vc,               80
.equiv          nb030_type,             88
.equiv          nb030_p_ntype,          96
.equiv          nb030_vdwparam,         104
.equiv          nb030_Vvdw,             112
.equiv          nb030_p_tabscale,       120
.equiv          nb030_VFtab,            128
.equiv          nb030_invsqrta,         136
.equiv          nb030_dvda,             144
.equiv          nb030_p_gbtabscale,     152
.equiv          nb030_GBtab,            160
.equiv          nb030_p_nthreads,       168
.equiv          nb030_count,            176
.equiv          nb030_mtx,              184
.equiv          nb030_outeriter,        192
.equiv          nb030_inneriter,        200
.equiv          nb030_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb030_ix,               0
.equiv          nb030_iy,               16
.equiv          nb030_iz,               32
.equiv          nb030_dx,               48
.equiv          nb030_dy,               64
.equiv          nb030_dz,               80
.equiv          nb030_two,              96
.equiv          nb030_tsc,              112
.equiv          nb030_c6,               128
.equiv          nb030_c12,              144
.equiv          nb030_fscal,            160
.equiv          nb030_Vvdwtot,          176
.equiv          nb030_fix,              192
.equiv          nb030_fiy,              208
.equiv          nb030_fiz,              224
.equiv          nb030_half,             240
.equiv          nb030_three,            256
.equiv          nb030_is3,              272
.equiv          nb030_ii3,              276
.equiv          nb030_nri,              280
.equiv          nb030_iinr,             288
.equiv          nb030_jindex,           296
.equiv          nb030_jjnr,             304
.equiv          nb030_shift,            312
.equiv          nb030_shiftvec,         320
.equiv          nb030_innerjjnr,        328
.equiv          nb030_ntia,             336
.equiv          nb030_innerk,           340
.equiv          nb030_n,                344
.equiv          nb030_nn1,              348
.equiv          nb030_ntype,            352
.equiv          nb030_nouter,           356
.equiv          nb030_ninner,           360

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
	sub rsp, 368		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb030_nouter], eax
	mov [rsp + nb030_ninner], eax



	mov edi, [rdi]
	mov [rsp + nb030_nri], edi
	mov [rsp + nb030_iinr], rsi
	mov [rsp + nb030_jindex], rdx
	mov [rsp + nb030_jjnr], rcx
	mov [rsp + nb030_shift], r8
	mov [rsp + nb030_shiftvec], r9
	mov rdi, [rbp + nb030_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb030_ntype], edi

	mov rax, [rbp + nb030_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb030_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb030_half], eax
	mov [rsp + nb030_half+4], ebx
	movsd xmm1, [rsp + nb030_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb030_half], xmm1
	movapd [rsp + nb030_two],  xmm2
	movapd [rsp + nb030_three], xmm3

.nb030_threadloop:
        mov   rsi, [rbp + nb030_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb030_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb030_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb030_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb030_n], eax
        mov [rsp + nb030_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb030_outerstart
        jmp .nb030_end

.nb030_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb030_nouter]
	mov [rsp + nb030_nouter], ebx

.nb030_outer:
	mov   rax, [rsp + nb030_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb030_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb030_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb030_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

    	mov   rdx, [rbp + nb030_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb030_ntype]
    	shl   edx, 1
    	mov   [rsp + nb030_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb030_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb030_ix], xmm0
	movapd [rsp + nb030_iy], xmm1
	movapd [rsp + nb030_iz], xmm2

	mov   [rsp + nb030_ii3], ebx
	
	;# clear tot potential and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb030_Vvdwtot], xmm4
	movapd [rsp + nb030_fix], xmm4
	movapd [rsp + nb030_fiy], xmm4
	movapd [rsp + nb030_fiz], xmm4
	
	mov   rax, [rsp + nb030_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb030_pos]
	mov   rdi, [rbp + nb030_faction]	
	mov   rax, [rsp + nb030_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb030_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb030_ninner]
	mov   [rsp + nb030_ninner], ecx
	add   edx, 0
	mov   [rsp + nb030_innerk], edx    ;# number of innerloop atoms 
	jge   .nb030_unroll_loop
	jmp   .nb030_checksingle
.nb030_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb030_innerjjnr]     ;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]
	add qword ptr [rsp + nb030_innerjjnr],  8 ;# advance pointer (unrolled 2) 
	
	mov rsi, [rbp + nb030_pos]		;# base of pos[] 
	lea   rax, [r8 + r8*2]	;# replace jnr with j3 
	lea   rbx, [r9 + r9*2]	

	;# move two coordinates to xmm4-xmm6 	
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
	movhpd xmm4, [rsi + rbx*8]
	movhpd xmm5, [rsi + rbx*8 + 8]
	movhpd xmm6, [rsi + rbx*8 + 16]

	;# calc dr 
	subpd xmm4, [rsp + nb030_ix]
	subpd xmm5, [rsp + nb030_iy]
	subpd xmm6, [rsp + nb030_iz]

	;# store dr 
	movapd [rsp + nb030_dx], xmm4
	movapd [rsp + nb030_dy], xmm5
	movapd [rsp + nb030_dz], xmm6
    
	mov rsi, [rbp + nb030_type]

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm4, xmm5
	addpd xmm4, xmm6
	;# rsq in xmm4 

	mov r8d, [rsi + r8*4]
	mov r9d, [rsi + r9*4]

	cvtpd2ps xmm5, xmm4	
	rsqrtps xmm5, xmm5
	cvtps2pd xmm2, xmm5	;# lu in low xmm2 

	;# lookup seed in xmm2 
	movapd xmm5, xmm2	;# copy of lu 
	mulpd xmm2, xmm2	;# lu*lu 
	movapd xmm1, [rsp + nb030_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb030_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	shl r8d, 1	
	shl r9d, 1	

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb030_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb030_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm2, xmm0	;# xmm2=iter2 of rinv (new lu) 
	
	mov edi, [rsp + nb030_ntia]
	add r8d, edi
	add r9d, edi

    mulpd xmm4, xmm2	;# xmm4=r 
	mulpd xmm4, [rsp + nb030_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4
    ;# xmm1=eps 
    ;# xmm2=rinv
    movapd xmm3, xmm4   ;# eps
	pslld mm6, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb030_VFtab]
	movd r10d, mm6
	psrlq mm6, 32
	movd r11d, mm6

    ;# indices in r10, r11. Load dispersion and repulsion tables in parallel.
	movapd xmm4, [rsi + r10*8]          ;# Y1d F1d	
	movapd xmm12, [rsi + r11*8]         ;# Y2d F2d 
	movapd xmm8, [rsi + r10*8 + 32]     ;# Y1r F1r 	
	movapd xmm13, [rsi + r11*8 + 32]	;# Y2r F2r 
	movapd xmm5, xmm4
	movapd xmm9, xmm8
	unpcklpd xmm4, xmm12	;# Y1d Y2d 
	unpckhpd xmm5, xmm12	;# F1d F2d 
	unpcklpd xmm8, xmm13	;# Y1r Y2r 
	unpckhpd xmm9, xmm13	;# F1r F2r 

	movapd xmm6, [rsi + r10*8 + 16]     ;# G1d H1d 	
	movapd xmm12, [rsi + r11*8 + 16]	;# G2d H2d 
	movapd xmm10, [rsi + r10*8 + 48]	;# G1r H1r 	
	movapd xmm13, [rsi + r11*8 + 48]	;# G2r H2r 
	movapd xmm7, xmm6
	movapd xmm11, xmm10
	unpcklpd xmm6, xmm12	;# G1d G2d 
	unpckhpd xmm7, xmm12	;# H1d H2d 
	unpcklpd xmm10, xmm13	;# G1r G2r 
	unpckhpd xmm11, xmm13	;# H1r H2r 
	;# tables ready, in xmm4-xmm7 and xmm8-xmm11
	mov rsi, [rbp + nb030_vdwparam]

    mulpd  xmm7, xmm1    ;# Heps
    mulpd  xmm11, xmm1 
	movlpd xmm12, [rsi + r8*8]	;# c6a
	movlpd xmm0, [rsi + r9*8]	;# c6b
    mulpd  xmm6, xmm1   ;# Geps
    mulpd  xmm10, xmm1 
    mulpd  xmm7, xmm1   ;# Heps2
    mulpd  xmm11, xmm1 

	movhpd xmm12, [rsi + r8*8 + 8]	;# c6a c12a 
	movhpd xmm0, [rsi + r9*8 + 8]	;# c6b c12b 

    addpd  xmm5, xmm6  ;# F+Geps
    addpd  xmm9, xmm10 
    addpd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addpd  xmm9, xmm11 
    addpd  xmm7, xmm7    ;# 2*Heps2
    addpd  xmm11, xmm11
    addpd  xmm7, xmm6   ;# 2*Heps2+Geps
    addpd  xmm11, xmm10
    
	movapd xmm13, xmm12
	unpcklpd xmm12, xmm0
	unpckhpd xmm13, xmm0

    addpd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm11, xmm9
    mulpd  xmm5, xmm1  ;# eps*Fp
    mulpd  xmm9, xmm1
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulpd  xmm9, xmm13  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb030_Vvdwtot]
    movapd [rsp + nb030_Vvdwtot], xmm5
        
    mulpd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulpd  xmm11, xmm13   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    
    mov rdi, [rbp + nb030_faction]
	;# the fj's - start by combining forces from memory 
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]

    mulpd  xmm7, [rsp + nb030_tsc]
    mulpd  xmm7, xmm2
    xorpd  xmm9, xmm9
    
    subpd  xmm9, xmm7
    movapd xmm10, xmm9
    movapd xmm11, xmm9
    
	movhpd xmm3, [rdi + rbx*8]
	movhpd xmm4, [rdi + rbx*8 + 8]
	movhpd xmm5, [rdi + rbx*8 + 16]

    movapd xmm12, [rsp + nb030_fix]
    movapd xmm13, [rsp + nb030_fiy]
    movapd xmm14, [rsp + nb030_fiz]
    
    mulpd  xmm9,  [rsp + nb030_dx]
    mulpd  xmm10, [rsp + nb030_dy]
    mulpd  xmm11, [rsp + nb030_dz]

    ;# accumulate i forces
    addpd xmm12, xmm9
    addpd xmm13, xmm10
    addpd xmm14, xmm11
    movapd [rsp + nb030_fix], xmm12
    movapd [rsp + nb030_fiy], xmm13
    movapd [rsp + nb030_fiz], xmm14
    
	addpd xmm3, xmm9
	addpd xmm4, xmm10
	addpd xmm5, xmm11
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5
	movhpd [rdi + rbx*8], xmm3
	movhpd [rdi + rbx*8 + 8], xmm4
	movhpd [rdi + rbx*8 + 16], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb030_innerk],  2
	jl    .nb030_checksingle
	jmp   .nb030_unroll_loop

.nb030_checksingle:				
	mov   edx, [rsp + nb030_innerk]
	and   edx, 1
	jnz    .nb030_dosingle
	jmp    .nb030_updateouterdata
.nb030_dosingle:
	mov   rdx, [rsp + nb030_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb030_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb030_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb030_ntia]
	add r8d, edi

	movsd xmm4, [rsi + r8*8]	
	movsd xmm6, [rsi + r8*8 + 8]
	movapd [rsp + nb030_c6], xmm4
	movapd [rsp + nb030_c12], xmm6
	
	mov rsi, [rbp + nb030_pos]		;# base of pos[] 
	lea   rax, [rax + rax*2]	;# replace jnr with j3 

	;# move two coordinates to xmm4-xmm6 	
	movsd xmm4, [rsi + rax*8]
	movsd xmm5, [rsi + rax*8 + 8]
	movsd xmm6, [rsi + rax*8 + 16]

	;# calc dr 
	subsd xmm4, [rsp + nb030_ix]
	subsd xmm5, [rsp + nb030_iy]
	subsd xmm6, [rsp + nb030_iz]

	;# store dr 
	movapd [rsp + nb030_dx], xmm4
	movapd [rsp + nb030_dy], xmm5
	movapd [rsp + nb030_dz], xmm6
    
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
	movapd xmm1, [rsp + nb030_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb030_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb030_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb030_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm2, xmm0	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm2	;# xmm4=r 
	mulsd xmm4, [rsp + nb030_tsc]
	
	cvttsd2si r10d, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, r10d
	subsd xmm4, xmm5
	movapd xmm1, xmm4
    ;# xmm1=eps 
    ;# xmm2=rinv
    movapd xmm3, xmm4   ;# eps
	shl    r10d, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb030_VFtab]

    ;# indices in r10, r11. Load dispersion and repulsion tables in parallel.
    movsd  xmm4,  [rsi + r10*8]
    movsd  xmm5,  [rsi + r10*8 + 8]
    movsd  xmm6,  [rsi + r10*8 + 16]
    movsd  xmm7,  [rsi + r10*8 + 24]
    movsd  xmm8,  [rsi + r10*8 + 32]
    movsd  xmm9,  [rsi + r10*8 + 40]
    movsd  xmm10, [rsi + r10*8 + 48]
    movsd  xmm11, [rsi + r10*8 + 56]
	;# tables ready, in xmm4-xmm7 and xmm8-xmm11

    mulsd  xmm7, xmm1    ;# Heps
    mulsd  xmm11, xmm1 
    mulsd  xmm6, xmm1   ;# Geps
    mulsd  xmm10, xmm1 
    mulsd  xmm7, xmm1   ;# Heps2
    mulsd  xmm11, xmm1 
    addsd  xmm5, xmm6  ;# F+Geps
    addsd  xmm9, xmm10 
    addsd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addsd  xmm9, xmm11 
    addsd  xmm7, xmm7    ;# 2*Heps2
    addsd  xmm11, xmm11
    addsd  xmm7, xmm6   ;# 2*Heps2+Geps
    addsd  xmm11, xmm10
    
    addsd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addsd  xmm11, xmm9
    mulsd  xmm5, xmm1  ;# eps*Fp
    mulsd  xmm9, xmm1
    movapd xmm12, [rsp + nb030_c6]
    movapd xmm13, [rsp + nb030_c12]
    addsd  xmm5, xmm4 ;# VV
    addsd  xmm9, xmm8

    mulsd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulsd  xmm9, xmm13  ;# VV*c12 = vnb12
    addsd  xmm5, xmm9
    addsd  xmm5, [rsp + nb030_Vvdwtot]
    movsd [rsp + nb030_Vvdwtot], xmm5
        
    mulsd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulsd  xmm11, xmm13   ;# FF*c12  = fnb12
    addsd  xmm7, xmm11
    
    mulsd  xmm7, [rsp + nb030_tsc]
    mulsd  xmm7, xmm2
    xorpd  xmm9, xmm9
    
    subsd  xmm9, xmm7
    movapd xmm10, xmm9
    movapd xmm11, xmm9
    
    movapd xmm12, [rsp + nb030_fix]
    movapd xmm13, [rsp + nb030_fiy]
    movapd xmm14, [rsp + nb030_fiz]
    
    mulsd  xmm9,  [rsp + nb030_dx]
    mulsd  xmm10, [rsp + nb030_dy]
    mulsd  xmm11, [rsp + nb030_dz]

    ;# accumulate i forces
    addsd xmm12, xmm9
    addsd xmm13, xmm10
    addsd xmm14, xmm11
    movsd [rsp + nb030_fix], xmm12
    movsd [rsp + nb030_fiy], xmm13
    movsd [rsp + nb030_fiz], xmm14
    
	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb030_faction]
	addsd xmm9,  [rdi + rax*8]
	addsd xmm10, [rdi + rax*8 + 8]
	addsd xmm11, [rdi + rax*8 + 16]
	movsd [rdi + rax*8], xmm9
	movsd [rdi + rax*8 + 8], xmm10
	movsd [rdi + rax*8 + 16], xmm11
	
.nb030_updateouterdata:
	mov   ecx, [rsp + nb030_ii3]
	mov   rdi, [rbp + nb030_faction]
	mov   rsi, [rbp + nb030_fshift]
	mov   edx, [rsp + nb030_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb030_fix]
	movapd xmm1, [rsp + nb030_fiy]
	movapd xmm2, [rsp + nb030_fiz]

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
	subsd xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rsi + rdx*8],     xmm3
	movsd  [rsi + rdx*8 + 8], xmm4
	movsd  [rsi + rdx*8 + 16], xmm5

	;# get n from stack
	mov esi, [rsp + nb030_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb030_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb030_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb030_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb030_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb030_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb030_n], esi
        jmp .nb030_outer
.nb030_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb030_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb030_end
        ;# non-zero, do one more workunit
        jmp   .nb030_threadloop
.nb030_end:

	mov eax, [rsp + nb030_nouter]
	mov ebx, [rsp + nb030_ninner]
	mov rcx, [rbp + nb030_outeriter]
	mov rdx, [rbp + nb030_inneriter]
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





.globl nb_kernel030nf_x86_64_sse2
.globl _nb_kernel030nf_x86_64_sse2
nb_kernel030nf_x86_64_sse2:	
_nb_kernel030nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb030nf_fshift,         16
.equiv          nb030nf_gid,            24
.equiv          nb030nf_pos,            32
.equiv          nb030nf_faction,        40
.equiv          nb030nf_charge,         48
.equiv          nb030nf_p_facel,        56
.equiv          nb030nf_argkrf,         64
.equiv          nb030nf_argcrf,         72
.equiv          nb030nf_Vc,             80
.equiv          nb030nf_type,           88
.equiv          nb030nf_p_ntype,        96
.equiv          nb030nf_vdwparam,       104
.equiv          nb030nf_Vvdw,           112
.equiv          nb030nf_p_tabscale,     120
.equiv          nb030nf_VFtab,          128
.equiv          nb030nf_invsqrta,       136
.equiv          nb030nf_dvda,           144
.equiv          nb030nf_p_gbtabscale,   152
.equiv          nb030nf_GBtab,          160
.equiv          nb030nf_p_nthreads,     168
.equiv          nb030nf_count,          176
.equiv          nb030nf_mtx,            184
.equiv          nb030nf_outeriter,      192
.equiv          nb030nf_inneriter,      200
.equiv          nb030nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb030nf_ix,             0
.equiv          nb030nf_iy,             16
.equiv          nb030nf_iz,             32
.equiv          nb030nf_tsc,            48
.equiv          nb030nf_c6,             64
.equiv          nb030nf_c12,            80
.equiv          nb030nf_Vvdwtot,        96
.equiv          nb030nf_half,           112
.equiv          nb030nf_three,          128
.equiv          nb030nf_is3,            144
.equiv          nb030nf_ii3,            148
.equiv          nb030nf_nri,            152
.equiv          nb030nf_iinr,           160
.equiv          nb030nf_jindex,         168
.equiv          nb030nf_jjnr,           176
.equiv          nb030nf_shift,          184
.equiv          nb030nf_shiftvec,       192
.equiv          nb030nf_innerjjnr,      200
.equiv          nb030nf_ntia,           208
.equiv          nb030nf_innerk,         212
.equiv          nb030nf_n,              216
.equiv          nb030nf_nn1,            220
.equiv          nb030nf_ntype,          224
.equiv          nb030nf_nouter,         228
.equiv          nb030nf_ninner,         232

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
	sub rsp, 240		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb030nf_nouter], eax
	mov [rsp + nb030nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb030nf_nri], edi
	mov [rsp + nb030nf_iinr], rsi
	mov [rsp + nb030nf_jindex], rdx
	mov [rsp + nb030nf_jjnr], rcx
	mov [rsp + nb030nf_shift], r8
	mov [rsp + nb030nf_shiftvec], r9
	mov rdi, [rbp + nb030nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb030nf_ntype], edi

	mov rax, [rbp + nb030nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb030nf_tsc], xmm3

	;# create constant floating-point factors on stack
	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb030nf_half], eax
	mov [rsp + nb030nf_half+4], ebx
	movsd xmm1, [rsp + nb030nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb030nf_half], xmm1
	movapd [rsp + nb030nf_three], xmm3

.nb030nf_threadloop:
        mov   rsi, [rbp + nb030nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb030nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb030nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb030nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb030nf_n], eax
        mov [rsp + nb030nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb030nf_outerstart
        jmp .nb030nf_end

.nb030nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb030nf_nouter]
	mov [rsp + nb030nf_nouter], ebx

.nb030nf_outer:
	mov   rax, [rsp + nb030nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb030nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb030nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

    	mov   rdx, [rbp + nb030nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb030nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb030nf_ntia], edx
	
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb030nf_pos]    ;# rax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb030nf_ix], xmm0
	movapd [rsp + nb030nf_iy], xmm1
	movapd [rsp + nb030nf_iz], xmm2

	mov   [rsp + nb030nf_ii3], ebx
	
	;# clear tot potential 
	xorpd xmm4, xmm4
	movapd [rsp + nb030nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb030nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb030nf_pos]
	mov   rax, [rsp + nb030nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb030nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb030nf_ninner]
	mov   [rsp + nb030nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb030nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb030nf_unroll_loop
	jmp   .nb030nf_checksingle
.nb030nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb030nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]
	add qword ptr [rsp + nb030nf_innerjjnr],  8 ;# advance pointer (unrolled 2) 

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb030nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb030nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb030nf_ntia]
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

	movapd [rsp + nb030nf_c6], xmm4
	movapd [rsp + nb030nf_c12], xmm6
	
	mov rsi, [rbp + nb030nf_pos]		;# base of pos[] 
	lea   rax, [rax + rax*2]	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]

	;# move nb030nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb030nf_ix]
	movapd xmm5, [rsp + nb030nf_iy]
	movapd xmm6, [rsp + nb030nf_iz]

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
	movapd xmm1, [rsp + nb030nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb030nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb030nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb030nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb030nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	movd mm0, eax	
	movd mm1, ebx

	mov  rsi, [rbp + nb030nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
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
	;# dispersion table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 

	mulpd  xmm5, [rsp + nb030nf_c6] ;# Vvdw6 

	;# Update Vvdwtot directly 
	addpd  xmm5, [rsp + nb030nf_Vvdwtot]
	movapd [rsp + nb030nf_Vvdwtot], xmm5

	;# repulsion 
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
	
	;# table ready, in xmm4-xmm7 	
	mulpd  xmm6, xmm1	;# xmm6=Geps 
	mulpd  xmm7, xmm2	;# xmm7=Heps2 
	addpd  xmm5, xmm6
	addpd  xmm5, xmm7	;# xmm5=Fp 	
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	
	mulpd  xmm5, [rsp + nb030nf_c12]
	
	addpd  xmm5, [rsp + nb030nf_Vvdwtot]
	movapd [rsp + nb030nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb030nf_innerk],  2
	jl    .nb030nf_checksingle
	jmp   .nb030nf_unroll_loop

.nb030nf_checksingle:				
	mov   edx, [rsp + nb030nf_innerk]
	and   edx, 1
	jnz    .nb030nf_dosingle
	jmp    .nb030nf_updateouterdata
.nb030nf_dosingle:
	mov   rdx, [rsp + nb030nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	movd  mm0, eax		;# use mmx registers as temp storage 
	
	mov rsi, [rbp + nb030nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb030nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb030nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		

	movapd [rsp + nb030nf_c6], xmm4
	movapd [rsp + nb030nf_c12], xmm6
	
	mov rsi, [rbp + nb030nf_pos]		;# base of pos[] 
	lea   rax, [rax + rax*2]	;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move nb030nf_ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb030nf_ix]
	movapd xmm5, [rsp + nb030nf_iy]
	movapd xmm6, [rsp + nb030nf_iz]

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
	movapd xmm1, [rsp + nb030nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb030nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb030nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb030nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=iter2 of rinv (new lu) 
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb030nf_tsc]

	movd mm0, eax
	
	cvttsd2si eax, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, eax
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	shl eax, 3	

	mov  rsi, [rbp + nb030nf_VFtab]

	;# dispersion 
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
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	mulsd  xmm5, [rsp + nb030nf_c6];# Vvdw6 

	;# Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb030nf_Vvdwtot]
	movlpd [rsp + nb030nf_Vvdwtot], xmm5

	;# repulsion 
	movapd xmm4, [rsi + rax*8 + 32]	;# Y1 F1 
	xorpd xmm3,xmm3	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movapd xmm6, [rsi + rax*8 + 48]	;# G1 H1 	
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1  
	unpckhpd xmm7, xmm3	;# H1  
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	mulsd  xmm5, [rsp + nb030nf_c12]
	
	addsd  xmm5, [rsp + nb030nf_Vvdwtot]
	movlpd [rsp + nb030nf_Vvdwtot], xmm5

.nb030nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb030nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb030nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb030nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
	
	;# add earlier value from mem 
	mov   rax, [rbp + nb030nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb030nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb030nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb030nf_n], esi
        jmp .nb030nf_outer
.nb030nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb030nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb030nf_end
        ;# non-zero, do one more workunit
        jmp   .nb030nf_threadloop
.nb030nf_end:
	mov eax, [rsp + nb030nf_nouter]
	mov ebx, [rsp + nb030nf_ninner]
	mov rcx, [rbp + nb030nf_outeriter]
	mov rdx, [rbp + nb030nf_inneriter]
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
