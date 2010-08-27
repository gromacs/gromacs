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



.globl nb_kernel130_x86_64_sse2
.globl _nb_kernel130_x86_64_sse2
nb_kernel130_x86_64_sse2:	
_nb_kernel130_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb130_fshift,           16
.equiv          nb130_gid,              24
.equiv          nb130_pos,              32
.equiv          nb130_faction,          40
.equiv          nb130_charge,           48
.equiv          nb130_p_facel,          56
.equiv          nb130_argkrf,           64
.equiv          nb130_argcrf,           72  
.equiv          nb130_Vc,               80
.equiv          nb130_type,             88
.equiv          nb130_p_ntype,          96
.equiv          nb130_vdwparam,         104
.equiv          nb130_Vvdw,             112
.equiv          nb130_p_tabscale,       120
.equiv          nb130_VFtab,            128
.equiv          nb130_invsqrta,         136
.equiv          nb130_dvda,             144
.equiv          nb130_p_gbtabscale,     152
.equiv          nb130_GBtab,            160
.equiv          nb130_p_nthreads,       168
.equiv          nb130_count,            176
.equiv          nb130_mtx,              184
.equiv          nb130_outeriter,        192
.equiv          nb130_inneriter,        200
.equiv          nb130_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb130_ix,               0
.equiv          nb130_iy,               16
.equiv          nb130_iz,               32
.equiv          nb130_iq,               48
.equiv          nb130_dx,               64
.equiv          nb130_dy,               80
.equiv          nb130_dz,               96
.equiv          nb130_c6,               112
.equiv          nb130_c12,              128
.equiv          nb130_tsc,              144
.equiv          nb130_fstmp,            160
.equiv          nb130_vctot,            176
.equiv          nb130_Vvdwtot,          192
.equiv          nb130_fix,              208
.equiv          nb130_fiy,              224
.equiv          nb130_fiz,              240
.equiv          nb130_half,             256
.equiv          nb130_three,            272
.equiv          nb130_two,              288
.equiv          nb130_krf,              304
.equiv          nb130_crf,              320
.equiv          nb130_nri,              336
.equiv          nb130_iinr,             344
.equiv          nb130_jindex,           352
.equiv          nb130_jjnr,             360
.equiv          nb130_shift,            368
.equiv          nb130_shiftvec,         376
.equiv          nb130_facel,            384
.equiv          nb130_innerjjnr,        392
.equiv          nb130_is3,              400
.equiv          nb130_ii3,              404
.equiv          nb130_ntia,             408
.equiv          nb130_innerk,           412
.equiv          nb130_n,                416
.equiv          nb130_nn1,              420
.equiv          nb130_ntype,            424
.equiv          nb130_nouter,           428
.equiv          nb130_ninner,           432

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
	sub rsp, 448		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb130_nouter], eax
	mov [rsp + nb130_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb130_nri], edi
	mov [rsp + nb130_iinr], rsi
	mov [rsp + nb130_jindex], rdx
	mov [rsp + nb130_jjnr], rcx
	mov [rsp + nb130_shift], r8
	mov [rsp + nb130_shiftvec], r9
	mov rdi, [rbp + nb130_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb130_ntype], edi
	mov rsi, [rbp + nb130_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb130_facel], xmm0

	mov rax, [rbp + nb130_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb130_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb130_half], eax
	mov [rsp + nb130_half+4], ebx
	movsd xmm1, [rsp + nb130_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb130_half], xmm1
	movapd [rsp + nb130_two], xmm2
	movapd [rsp + nb130_three], xmm3

.nb130_threadloop:
        mov   rsi, [rbp + nb130_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb130_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb130_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb130_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb130_n], eax
        mov [rsp + nb130_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb130_outerstart
        jmp .nb130_end

.nb130_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb130_nouter]
	mov [rsp + nb130_nouter], ebx

.nb130_outer:
	mov   rax, [rsp + nb130_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb130_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb130_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb130_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb130_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb130_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb130_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb130_ntype]
    	shl   edx, 1
    	mov   [rsp + nb130_ntia], edx
		
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb130_pos]    ;# eax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb130_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb130_ix], xmm0
	movapd [rsp + nb130_iy], xmm1
	movapd [rsp + nb130_iz], xmm2

	mov   [rsp + nb130_ii3], ebx
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb130_vctot], xmm4
	movapd [rsp + nb130_Vvdwtot], xmm4
	movapd [rsp + nb130_fix], xmm4
	movapd [rsp + nb130_fiy], xmm4
	movapd [rsp + nb130_fiz], xmm4
	
	mov   rax, [rsp + nb130_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb130_pos]
	mov   rdi, [rbp + nb130_faction]	
	mov   rax, [rsp + nb130_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb130_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb130_ninner]
	mov   [rsp + nb130_ninner], ecx
	add   edx, 0
	mov   [rsp + nb130_innerk], edx    ;# number of innerloop atoms 
	jge   .nb130_unroll_loop
	jmp   .nb130_checksingle
.nb130_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb130_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add qword ptr [rsp + nb130_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb130_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	mulpd xmm3, [rsp + nb130_iq]		;# qq 
	
	mov rsi, [rbp + nb130_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb130_vdwparam]
	shl r8d, 1
	shl r9d, 1
	mov edi, [rsp + nb130_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6
	movhpd xmm6, [rsi + r9*8]	
	movlpd xmm7, [rsi + r8*8 + 8]	;# c12 
	movhpd xmm7, [rsi + r9*8 + 8]	
	movapd [rsp + nb130_c6], xmm6
	movapd [rsp + nb130_c12], xmm7
	
	mov rsi, [rbp + nb130_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm4-xmm6 	
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
	movhpd xmm4, [rsi + rbx*8]
	movhpd xmm5, [rsi + rbx*8 + 8]
	movhpd xmm6, [rsi + rbx*8 + 16]		
	
	;# calc dr 
	subpd xmm4, [rsp + nb130_ix]
	subpd xmm5, [rsp + nb130_iy]
	subpd xmm6, [rsp + nb130_iz]

	;# store dr 
	movapd [rsp + nb130_dx], xmm4
	movapd [rsp + nb130_dy], xmm5
	movapd [rsp + nb130_dz], xmm6
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
	movapd xmm1, [rsp + nb130_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb130_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb130_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb130_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm6, xmm0
	movapd xmm1, xmm4
	mulpd  xmm6, xmm3   ;# vcoul = rinv*qq
	movapd xmm3, xmm6
	mulpd  xmm3, xmm0
	
    ;# fstmp in xmm3
    
	addpd  xmm6, [rsp + nb130_vctot]
	movapd [rsp + nb130_vctot], xmm6

	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb130_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 

	pslld mm6, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb130_VFtab]
	movd r8d, mm6
	psrlq mm6, 32
	movd r9d, mm6

	;# load both disp. and rep. tables in parallel
	movlpd xmm4,  [rsi + r8*8]	    
	movlpd xmm5,  [rsi + r8*8 + 8]	
	movlpd xmm6,  [rsi + r8*8 + 16]	
	movlpd xmm7,  [rsi + r8*8 + 24]	
	movlpd xmm8,  [rsi + r8*8 + 32]	
	movlpd xmm9,  [rsi + r8*8 + 40]	
	movlpd xmm10, [rsi + r8*8 + 48]	
	movlpd xmm11, [rsi + r8*8 + 56]	
	movhpd xmm4,  [rsi + r9*8]	    
	movhpd xmm5,  [rsi + r9*8 + 8]	
	movhpd xmm6,  [rsi + r9*8 + 16]	
	movhpd xmm7,  [rsi + r9*8 + 24]	
	movhpd xmm8,  [rsi + r9*8 + 32]	
	movhpd xmm9,  [rsi + r9*8 + 40]	
	movhpd xmm10, [rsi + r9*8 + 48]	
	movhpd xmm11, [rsi + r9*8 + 56]	
	;# dispersion table ready in xmm4-xmm7, repulsion in xmm8-xmm11
    
    mulpd  xmm7, xmm1    ;# Heps
    mulpd  xmm11, xmm1 
    mulpd  xmm6, xmm1   ;# Geps
    mulpd  xmm10, xmm1 
    mulpd  xmm7, xmm1   ;# Heps2
    mulpd  xmm11, xmm1 
    addpd  xmm5, xmm6  ;# F+Geps
    addpd  xmm9, xmm10 
    addpd  xmm5, xmm7   ;# F+Geps+Heps2 = Fp
    addpd  xmm9, xmm11 
    addpd  xmm7, xmm7    ;# 2*Heps2
    addpd  xmm11, xmm11
    addpd  xmm7, xmm6   ;# 2*Heps2+Geps
    addpd  xmm11, xmm10
    
    addpd  xmm7, xmm5  ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm11, xmm9
    mulpd  xmm5, xmm1  ;# eps*Fp
    mulpd  xmm9, xmm1
    movapd xmm12, [rsp + nb130_c6]
    movapd xmm13, [rsp + nb130_c12]
    addpd  xmm5, xmm4 ;# VV
    addpd  xmm9, xmm8

    mulpd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulpd  xmm9, xmm13  ;# VV*c12 = vnb12
    addpd  xmm5, xmm9
    addpd  xmm5, [rsp + nb130_Vvdwtot]
    movapd [rsp + nb130_Vvdwtot], xmm5
        
    mulpd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulpd  xmm11, xmm13   ;# FF*c12  = fnb12
    addpd  xmm7, xmm11
    
    mulpd  xmm7, [rsp + nb130_tsc]
    subpd  xmm3, xmm7
    mulpd  xmm3, xmm0   ;# fscal

    movapd xmm9, xmm3
    movapd xmm10, xmm3
    movapd xmm11, xmm3
    
    movapd xmm12, [rsp + nb130_fix]
    movapd xmm13, [rsp + nb130_fiy]
    movapd xmm14, [rsp + nb130_fiz]
    
    mulpd  xmm9,  [rsp + nb130_dx]
    mulpd  xmm10, [rsp + nb130_dy]
    mulpd  xmm11, [rsp + nb130_dz]

    ;# accumulate i forces
    addpd xmm12, xmm9
    addpd xmm13, xmm10
    addpd xmm14, xmm11
    movapd [rsp + nb130_fix], xmm12
    movapd [rsp + nb130_fiy], xmm13
    movapd [rsp + nb130_fiz], xmm14
    
	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb130_faction]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	movhpd xmm3, [rdi + rbx*8]
	movhpd xmm4, [rdi + rbx*8 + 8]
	movhpd xmm5, [rdi + rbx*8 + 16]
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
	sub dword ptr [rsp + nb130_innerk],  2
	jl    .nb130_checksingle
	jmp   .nb130_unroll_loop

.nb130_checksingle:				
	mov   edx, [rsp + nb130_innerk]
	and   edx, 1
	jnz    .nb130_dosingle
	jmp    .nb130_updateouterdata
.nb130_dosingle:			
	mov rsi, [rbp + nb130_charge]
	mov rdi, [rbp + nb130_pos]
	mov   rcx, [rsp + nb130_innerjjnr]
	mov   eax, [rcx]

	mov rsi, [rbp + nb130_charge]    ;# base of charge[] 
	
	movsd xmm3, [rsi + rax*8]
	mulsd xmm3, [rsp + nb130_iq]		;# qq 
	
	mov rsi, [rbp + nb130_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb130_vdwparam]
	shl r8d, 1
	mov edi, [rsp + nb130_ntia]
	add r8d, edi

	movsd xmm4, [rsi + r8*8]	
	movsd xmm6, [rsi + r8*8 + 8]
	movapd [rsp + nb130_c6], xmm4
	movapd [rsp + nb130_c12], xmm6
	
	mov rsi, [rbp + nb130_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinate to xmm4-xmm6 	
	movsd xmm4, [rsi + rax*8]
	movsd xmm5, [rsi + rax*8 + 8]
	movsd xmm6, [rsi + rax*8 + 16]
	
	;# calc dr 
	subsd xmm4, [rsp + nb130_ix]
	subsd xmm5, [rsp + nb130_iy]
	subsd xmm6, [rsp + nb130_iz]

	;# store dr 
	movapd [rsp + nb130_dx], xmm4
	movapd [rsp + nb130_dy], xmm5
	movapd [rsp + nb130_dz], xmm6
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
	movapd xmm1, [rsp + nb130_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb130_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb130_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb130_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm6, xmm0
	movapd xmm1, xmm4
	mulsd  xmm6, xmm3   ;# vcoul = rinv*qq
	movapd xmm3, xmm6
	mulsd  xmm3, xmm0
	
    ;# fstmp in xmm3
    
	addsd  xmm6, [rsp + nb130_vctot]
	movsd [rsp + nb130_vctot], xmm6

	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb130_tsc]
	
	cvttsd2si r8d, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, r8d
	subsd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 

	shl     r8d, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb130_VFtab]

	;# load both disp. and rep. tables in parallel
    movsd xmm4,  [rsi + r8*8]
    movsd xmm5,  [rsi + r8*8 + 8]
    movsd xmm6,  [rsi + r8*8 + 16]
    movsd xmm7,  [rsi + r8*8 + 24]
    movsd xmm8,  [rsi + r8*8 + 32]
    movsd xmm9,  [rsi + r8*8 + 40]
    movsd xmm10, [rsi + r8*8 + 48]
    movsd xmm11, [rsi + r8*8 + 56]
	;# dispersion table ready in xmm4-xmm7, repulsion in xmm8-xmm11
    
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
    movapd xmm12, [rsp + nb130_c6]
    movapd xmm13, [rsp + nb130_c12]
    addsd  xmm5, xmm4 ;# VV
    addsd  xmm9, xmm8

    mulsd  xmm5, xmm12  ;# VV*c6 = vnb6
    mulsd  xmm9, xmm13  ;# VV*c12 = vnb12
    addsd  xmm5, xmm9
    addsd  xmm5, [rsp + nb130_Vvdwtot]
    movsd [rsp + nb130_Vvdwtot], xmm5
        
    mulsd  xmm7, xmm12   ;# FF*c6 = fnb6
    mulsd  xmm11, xmm13   ;# FF*c12  = fnb12
    addsd  xmm7, xmm11
    
    mulsd  xmm7, [rsp + nb130_tsc]
    subsd  xmm3, xmm7
    mulsd  xmm3, xmm0   ;# fscal

    movapd xmm9, xmm3
    movapd xmm10, xmm3
    movapd xmm11, xmm3
    
    movapd xmm12, [rsp + nb130_fix]
    movapd xmm13, [rsp + nb130_fiy]
    movapd xmm14, [rsp + nb130_fiz]
    
    mulsd  xmm9,  [rsp + nb130_dx]
    mulsd  xmm10, [rsp + nb130_dy]
    mulsd  xmm11, [rsp + nb130_dz]

    ;# accumulate i forces
    addsd xmm12, xmm9
    addsd xmm13, xmm10
    addsd xmm14, xmm11
    movsd [rsp + nb130_fix], xmm12
    movsd [rsp + nb130_fiy], xmm13
    movsd [rsp + nb130_fiz], xmm14
    
	;# the fj's - start by accumulating forces from memory 
    mov rdi, [rbp + nb130_faction]
	addsd xmm9,  [rdi + rax*8]
	addsd xmm10, [rdi + rax*8 + 8]
	addsd xmm11, [rdi + rax*8 + 16]
	movsd [rdi + rax*8], xmm9
	movsd [rdi + rax*8 + 8], xmm10
	movsd [rdi + rax*8 + 16], xmm11
	
.nb130_updateouterdata:
	mov   ecx, [rsp + nb130_ii3]
	mov   rdi, [rbp + nb130_faction]
	mov   rsi, [rbp + nb130_fshift]
	mov   edx, [rsp + nb130_is3]

	;# accumulate i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb130_fix]
	movapd xmm1, [rsp + nb130_fiy]
	movapd xmm2, [rsp + nb130_fiz]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addpd  xmm0, xmm3
	addpd  xmm1, xmm4
	addpd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

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
	mov esi, [rsp + nb130_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb130_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb130_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb130_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb130_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb130_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb130_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb130_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb130_n], esi
        jmp .nb130_outer
.nb130_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb130_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb130_end
        ;# non-zero, do one more workunit
        jmp   .nb130_threadloop
.nb130_end:
	mov eax, [rsp + nb130_nouter]
	mov ebx, [rsp + nb130_ninner]
	mov rcx, [rbp + nb130_outeriter]
	mov rdx, [rbp + nb130_inneriter]
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

	


.globl nb_kernel130nf_x86_64_sse2
.globl _nb_kernel130nf_x86_64_sse2
nb_kernel130nf_x86_64_sse2:	
_nb_kernel130nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb130nf_fshift,         16
.equiv          nb130nf_gid,            24
.equiv          nb130nf_pos,            32
.equiv          nb130nf_faction,        40
.equiv          nb130nf_charge,         48
.equiv          nb130nf_p_facel,        56
.equiv          nb130nf_argkrf,         64
.equiv          nb130nf_argcrf,         72
.equiv          nb130nf_Vc,             80
.equiv          nb130nf_type,           88
.equiv          nb130nf_p_ntype,        96
.equiv          nb130nf_vdwparam,       104
.equiv          nb130nf_Vvdw,           112
.equiv          nb130nf_p_tabscale,     120
.equiv          nb130nf_VFtab,          128
.equiv          nb130nf_invsqrta,       136
.equiv          nb130nf_dvda,           144
.equiv          nb130nf_p_gbtabscale,   152
.equiv          nb130nf_GBtab,          160
.equiv          nb130nf_p_nthreads,     168
.equiv          nb130nf_count,          176
.equiv          nb130nf_mtx,            184
.equiv          nb130nf_outeriter,      192
.equiv          nb130nf_inneriter,      200
.equiv          nb130nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb130nf_ix,             0
.equiv          nb130nf_iy,             16
.equiv          nb130nf_iz,             32
.equiv          nb130nf_iq,             48
.equiv          nb130nf_c6,             64
.equiv          nb130nf_c12,            80
.equiv          nb130nf_vctot,          96
.equiv          nb130nf_Vvdwtot,        112
.equiv          nb130nf_half,           128
.equiv          nb130nf_three,          144
.equiv          nb130nf_krf,            160
.equiv          nb130nf_crf,            176
.equiv          nb130nf_tsc,            192
.equiv          nb130nf_nri,            208
.equiv          nb130nf_iinr,           216
.equiv          nb130nf_jindex,         224
.equiv          nb130nf_jjnr,           232
.equiv          nb130nf_shift,          240
.equiv          nb130nf_shiftvec,       248
.equiv          nb130nf_facel,          256
.equiv          nb130nf_innerjjnr,      264
.equiv          nb130nf_is3,            272
.equiv          nb130nf_ii3,            280
.equiv          nb130nf_ntia,           284
.equiv          nb130nf_innerk,         288
.equiv          nb130nf_n,              292
.equiv          nb130nf_nn1,            296
.equiv          nb130nf_ntype,          300
.equiv          nb130nf_nouter,         304
.equiv          nb130nf_ninner,         308

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
	sub rsp, 320		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb130nf_nouter], eax
	mov [rsp + nb130nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb130nf_nri], edi
	mov [rsp + nb130nf_iinr], rsi
	mov [rsp + nb130nf_jindex], rdx
	mov [rsp + nb130nf_jjnr], rcx
	mov [rsp + nb130nf_shift], r8
	mov [rsp + nb130nf_shiftvec], r9
	mov rdi, [rbp + nb130nf_p_ntype]
	mov edi, [rdi]
	mov [rsp + nb130nf_ntype], edi
	mov rsi, [rbp + nb130nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb130nf_facel], xmm0

	mov rax, [rbp + nb130nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb130nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb130nf_half], eax
	mov [rsp + nb130nf_half+4], ebx
	movsd xmm1, [rsp + nb130nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb130nf_half], xmm1
	movapd [rsp + nb130nf_three], xmm3

.nb130nf_threadloop:
        mov   rsi, [rbp + nb130nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb130nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb130nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb130nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb130nf_n], eax
        mov [rsp + nb130nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb130nf_outerstart
        jmp .nb130nf_end

.nb130nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb130nf_nouter]
	mov [rsp + nb130nf_nouter], ebx

.nb130nf_outer:
	mov   rax, [rsp + nb130nf_shift]      ;# eax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb130nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb130nf_shiftvec]   ;# eax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb130nf_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	mov   rdx, [rbp + nb130nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	mulsd xmm3, [rsp + nb130nf_facel]
	shufpd xmm3, xmm3, 0

    	mov   rdx, [rbp + nb130nf_type] 
    	mov   edx, [rdx + rbx*4]
    	imul  edx, [rsp + nb130nf_ntype]
    	shl   edx, 1
    	mov   [rsp + nb130nf_ntia], edx
		
	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb130nf_pos]    ;# eax = base of pos[]  

	addsd xmm0, [rax + rbx*8]
	addsd xmm1, [rax + rbx*8 + 8]
	addsd xmm2, [rax + rbx*8 + 16]

	movapd [rsp + nb130nf_iq], xmm3
	
	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0

	movapd [rsp + nb130nf_ix], xmm0
	movapd [rsp + nb130nf_iy], xmm1
	movapd [rsp + nb130nf_iz], xmm2

	mov   [rsp + nb130nf_ii3], ebx
	
	;# clear vctot
	xorpd xmm4, xmm4
	movapd [rsp + nb130nf_vctot], xmm4
	movapd [rsp + nb130nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb130nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb130nf_pos]
	mov   rdi, [rbp + nb130nf_faction]	
	mov   rax, [rsp + nb130nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb130nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb130nf_ninner]
	mov   [rsp + nb130nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb130nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb130nf_unroll_loop
	jmp   .nb130nf_checksingle
.nb130nf_unroll_loop:	
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb130nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	add qword ptr [rsp + nb130nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb130nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]

	movapd xmm5, [rsp + nb130nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	
	mov rsi, [rbp + nb130nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb130nf_vdwparam]
	shl eax, 1
	shl ebx, 1
	mov edi, [rsp + nb130nf_ntia]
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
	movapd [rsp + nb130nf_c6], xmm4
	movapd [rsp + nb130nf_c12], xmm6
	
	mov rsi, [rbp + nb130nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb130nf_ix]
	movapd xmm5, [rsp + nb130nf_iy]
	movapd xmm6, [rsp + nb130nf_iz]

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
	movapd xmm1, [rsp + nb130nf_three]
	mulpd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb130nf_half]
	subpd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulpd xmm1, xmm5	
	mulpd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulpd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb130nf_three]
	mulpd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb130nf_half]
	subpd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulpd xmm2, xmm5	
	mulpd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm6, xmm0
	mulpd  xmm6, xmm3   ;# vcoul
	
	addpd  xmm6, [rsp + nb130nf_vctot]
	movapd [rsp + nb130nf_vctot], xmm6

	;# LJ table interaction. xmm0=rinv, xmm4=rsq
	
	mulpd xmm4, xmm0	;# xmm4=r 
	mulpd xmm4, [rsp + nb130nf_tsc]
	
	cvttpd2pi mm6, xmm4	;# mm6 = lu idx 
	cvtpi2pd xmm5, mm6
	subpd xmm4, xmm5
	movapd xmm1, xmm4	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 

	pslld mm6, 3		;# idx *= 8 
	
	mov  rsi, [rbp + nb130nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6

	;# dispersion 
	movlpd xmm4, [rsi + rax*8]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 8]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 8]	;# Y2 F2 
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 16]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 16]	;# G2
	movhpd xmm6, [rsi + rax*8 + 24]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 24]	;# G2 H2 
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

	movapd xmm4, [rsp + nb130nf_c6]
	mulpd  xmm5, xmm4	 ;# Vvdw6 

	;#Update Vvdwtot directly 
	addpd  xmm5, [rsp + nb130nf_Vvdwtot]
	movapd [rsp + nb130nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rax*8 + 32]	;# Y1 	
	movlpd xmm3, [rsi + rbx*8 + 32]	;# Y2 
	movhpd xmm4, [rsi + rax*8 + 40]	;# Y1 F1 	
	movhpd xmm3, [rsi + rbx*8 + 40]	;# Y2 F2 

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rax*8 + 48]	;# G1
	movlpd xmm3, [rsi + rbx*8 + 48]	;# G2
	movhpd xmm6, [rsi + rax*8 + 56]	;# G1 H1 	
	movhpd xmm3, [rsi + rbx*8 + 56]	;# G2 H2 

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
	
	movapd xmm4, [rsp + nb130nf_c12]
	mulpd  xmm5, xmm4  
	
	addpd  xmm5, [rsp + nb130nf_Vvdwtot]
	movapd [rsp + nb130nf_Vvdwtot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb130nf_innerk],  2
	jl    .nb130nf_checksingle
	jmp   .nb130nf_unroll_loop

.nb130nf_checksingle:				
	mov   edx, [rsp + nb130nf_innerk]
	and   edx, 1
	jnz    .nb130nf_dosingle
	jmp    .nb130nf_updateouterdata
.nb130nf_dosingle:			
	mov rsi, [rbp + nb130nf_charge]
	mov rdi, [rbp + nb130nf_pos]
	mov   rcx, [rsp + nb130nf_innerjjnr]
	xorpd xmm3, xmm3
	mov   eax, [rcx]

	movlpd xmm3, [rsi + rax*8]
	movapd xmm5, [rsp + nb130nf_iq]
	mulpd xmm3, xmm5		;# qq 
	
	movd  mm0, eax		;# use mmx registers as temp storage 
	mov rsi, [rbp + nb130nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb130nf_vdwparam]
	shl eax, 1
	mov edi, [rsp + nb130nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 
	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0		
	movapd [rsp + nb130nf_c6], xmm4
	movapd [rsp + nb130nf_c12], xmm6
	
	mov rsi, [rbp + nb130nf_pos]       ;# base of pos[] 

	lea rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	
	;# move ix-iz to xmm4-xmm6 
	movapd xmm4, [rsp + nb130nf_ix]
	movapd xmm5, [rsp + nb130nf_iy]
	movapd xmm6, [rsp + nb130nf_iz]

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
	movapd xmm1, [rsp + nb130nf_three]
	mulsd xmm2, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb130nf_half]
	subsd xmm1, xmm2	;# 30-rsq*lu*lu 
	mulsd xmm1, xmm5	
	mulsd xmm1, xmm0	;# xmm0=iter1 of rinv (new lu) 

	movapd xmm5, xmm1	;# copy of lu 
	mulsd xmm1, xmm1	;# lu*lu 
	movapd xmm2, [rsp + nb130nf_three]
	mulsd xmm1, xmm4	;# rsq*lu*lu 			
	movapd xmm0, [rsp + nb130nf_half]
	subsd xmm2, xmm1	;# 30-rsq*lu*lu 
	mulsd xmm2, xmm5	
	mulsd xmm0, xmm2	;# xmm0=rinv 
	movapd xmm6, xmm0
	movapd xmm1, xmm4
	mulsd  xmm6, xmm3	;# xmm6=vcoul=qq*rinv
	
	addsd  xmm6, [rsp + nb130nf_vctot]
	movsd [rsp + nb130nf_vctot], xmm6

	;# LJ table interaction. xmm0=rinv, cmm4=rsq
	
	mulsd xmm4, xmm0	;# xmm4=r 
	mulsd xmm4, [rsp + nb130nf_tsc]
	
	cvttsd2si ebx, xmm4	;# mm6 = lu idx 
	cvtsi2sd xmm5, ebx
	subsd xmm4, xmm5
	movsd xmm1, xmm4	;# xmm1=eps 
	movsd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 

	shl   ebx, 3

	mov  rsi, [rbp + nb130nf_VFtab]

	;# dispersion 
	movlpd xmm4, [rsi + rbx*8]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 8]	;# Y1 F1 	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 16]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 24]	;# G1 H1 	
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	;# dispersion table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 

	movsd xmm4, [rsp + nb130nf_c6]
	mulsd  xmm5, xmm4	 ;# Vvdw6 

	;# put scalar force on stack Update Vvdwtot directly 
	addsd  xmm5, [rsp + nb130nf_Vvdwtot]
	movsd [rsp + nb130nf_Vvdwtot], xmm5

	;# repulsion 
	movlpd xmm4, [rsi + rbx*8 + 32]	;# Y1 	
	movhpd xmm4, [rsi + rbx*8 + 40]	;# Y1 F1 	

	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 Y2 
	unpckhpd xmm5, xmm3	;# F1 F2 

	movlpd xmm6, [rsi + rbx*8 + 48]	;# G1
	movhpd xmm6, [rsi + rbx*8 + 56]	;# G1 H1 	

	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 G2 
	unpckhpd xmm7, xmm3	;# H1 H2 
	
	;# table ready, in xmm4-xmm7 	
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	
	movsd xmm4, [rsp + nb130nf_c12]
	mulsd  xmm5, xmm4  
	
	addsd  xmm5, [rsp + nb130nf_Vvdwtot]
	movsd [rsp + nb130nf_Vvdwtot], xmm5
	
.nb130nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb130nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb130nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb130nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb130nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb130nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb130nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb130nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb130nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb130nf_n], esi
        jmp .nb130nf_outer
.nb130nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb130nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb130nf_end
        ;# non-zero, do one more workunit
        jmp   .nb130nf_threadloop
.nb130nf_end:
	mov eax, [rsp + nb130nf_nouter]
	mov ebx, [rsp + nb130nf_ninner]
	mov rcx, [rbp + nb130nf_outeriter]
	mov rdx, [rbp + nb130nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 320
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
