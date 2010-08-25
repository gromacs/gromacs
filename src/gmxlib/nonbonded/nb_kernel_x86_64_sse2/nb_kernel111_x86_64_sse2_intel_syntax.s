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


.globl nb_kernel111_x86_64_sse2
.globl _nb_kernel111_x86_64_sse2
nb_kernel111_x86_64_sse2:	
_nb_kernel111_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb111_fshift,           16
.equiv          nb111_gid,              24
.equiv          nb111_pos,              32
.equiv          nb111_faction,          40
.equiv          nb111_charge,           48
.equiv          nb111_p_facel,          56
.equiv          nb111_argkrf,           64
.equiv          nb111_argcrf,           72
.equiv          nb111_Vc,               80
.equiv          nb111_type,             88
.equiv          nb111_p_ntype,          96
.equiv          nb111_vdwparam,         104
.equiv          nb111_Vvdw,             112
.equiv          nb111_p_tabscale,       120
.equiv          nb111_VFtab,            128
.equiv          nb111_invsqrta,         136
.equiv          nb111_dvda,             144
.equiv          nb111_p_gbtabscale,     152
.equiv          nb111_GBtab,            160
.equiv          nb111_p_nthreads,       168
.equiv          nb111_count,            176
.equiv          nb111_mtx,              184
.equiv          nb111_outeriter,        192
.equiv          nb111_inneriter,        200
.equiv          nb111_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb111_ixO,              0
.equiv          nb111_iyO,              16
.equiv          nb111_izO,              32
.equiv          nb111_ixH1,             48
.equiv          nb111_iyH1,             64
.equiv          nb111_izH1,             80
.equiv          nb111_ixH2,             96
.equiv          nb111_iyH2,             112
.equiv          nb111_izH2,             128
.equiv          nb111_iqO,              144
.equiv          nb111_iqH,              160
.equiv          nb111_dxO,              176
.equiv          nb111_dyO,              192
.equiv          nb111_dzO,              208
.equiv          nb111_dxH1,             224
.equiv          nb111_dyH1,             240
.equiv          nb111_dzH1,             256
.equiv          nb111_dxH2,             272
.equiv          nb111_dyH2,             288
.equiv          nb111_dzH2,             304
.equiv          nb111_qqO,              320
.equiv          nb111_qqH,              336
.equiv          nb111_c6,               352
.equiv          nb111_c12,              368
.equiv          nb111_six,              384
.equiv          nb111_twelve,           400
.equiv          nb111_vctot,            416
.equiv          nb111_Vvdwtot,          432
.equiv          nb111_fixO,             448
.equiv          nb111_fiyO,             464
.equiv          nb111_fizO,             480
.equiv          nb111_fixH1,            496
.equiv          nb111_fiyH1,            512
.equiv          nb111_fizH1,            528
.equiv          nb111_fixH2,            544
.equiv          nb111_fiyH2,            560
.equiv          nb111_fizH2,            576
.equiv          nb111_fjx,              592
.equiv          nb111_fjy,              608
.equiv          nb111_fjz,              624
.equiv          nb111_half,             640
.equiv          nb111_three,            656
.equiv          nb111_is3,              672
.equiv          nb111_ii3,              676
.equiv          nb111_nri,              692
.equiv          nb111_iinr,             700
.equiv          nb111_jindex,           708
.equiv          nb111_jjnr,             716
.equiv          nb111_shift,            724
.equiv          nb111_shiftvec,         732
.equiv          nb111_facel,            740
.equiv          nb111_innerjjnr,        748
.equiv          nb111_ntia,             756
.equiv          nb111_innerk,           760
.equiv          nb111_n,                764
.equiv          nb111_nn1,              768
.equiv          nb111_nouter,           772
.equiv          nb111_ninner,           776

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
	sub rsp, 784		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb111_nouter], eax
	mov [rsp + nb111_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb111_nri], edi
	mov [rsp + nb111_iinr], rsi
	mov [rsp + nb111_jindex], rdx
	mov [rsp + nb111_jjnr], rcx
	mov [rsp + nb111_shift], r8
	mov [rsp + nb111_shiftvec], r9
	mov rsi, [rbp + nb111_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb111_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb111_half], eax
	mov [rsp + nb111_half+4], ebx
	movsd xmm1, [rsp + nb111_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd xmm4, xmm3
	addpd  xmm4, xmm4       ;# six
	movapd xmm5, xmm4
	addpd  xmm5, xmm5       ;# twelve
	movapd [rsp + nb111_half], xmm1
	movapd [rsp + nb111_three], xmm3
	movapd [rsp + nb111_six], xmm4
	movapd [rsp + nb111_twelve], xmm5

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb111_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb111_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb111_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb111_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb111_iqO], xmm3
	movapd [rsp + nb111_iqH], xmm4
	
	mov   rdx, [rbp + nb111_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb111_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb111_ntia], ecx		
.nb111_threadloop:
        mov   rsi, [rbp + nb111_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb111_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb111_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb111_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb111_n], eax
        mov [rsp + nb111_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb111_outerstart
        jmp .nb111_end

.nb111_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb111_nouter]
	mov [rsp + nb111_nouter], ebx

.nb111_outer:
	mov   rax, [rsp + nb111_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb111_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb111_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb111_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb111_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb111_ii3], ebx

	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb111_ixO], xmm3
	movapd [rsp + nb111_iyO], xmm4
	movapd [rsp + nb111_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 24]
	addsd xmm1, [rax + rbx*8 + 32]
	addsd xmm2, [rax + rbx*8 + 40]		
	addsd xmm3, [rax + rbx*8 + 48]
	addsd xmm4, [rax + rbx*8 + 56]
	addsd xmm5, [rax + rbx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb111_ixH1], xmm0
	movapd [rsp + nb111_iyH1], xmm1
	movapd [rsp + nb111_izH1], xmm2
	movapd [rsp + nb111_ixH2], xmm3
	movapd [rsp + nb111_iyH2], xmm4
	movapd [rsp + nb111_izH2], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb111_vctot], xmm4
	movapd [rsp + nb111_Vvdwtot], xmm4
	movapd [rsp + nb111_fixO], xmm4
	movapd [rsp + nb111_fiyO], xmm4
	movapd [rsp + nb111_fizO], xmm4
	movapd [rsp + nb111_fixH1], xmm4
	movapd [rsp + nb111_fiyH1], xmm4
	movapd [rsp + nb111_fizH1], xmm4
	movapd [rsp + nb111_fixH2], xmm4
	movapd [rsp + nb111_fiyH2], xmm4
	movapd [rsp + nb111_fizH2], xmm4
	
	mov   rax, [rsp + nb111_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb111_pos]
	mov   rdi, [rbp + nb111_faction]	
	mov   rax, [rsp + nb111_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb111_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb111_ninner]
	mov   [rsp + nb111_ninner], ecx
	add   edx, 0
	mov   [rsp + nb111_innerk], edx    ;# number of innerloop atoms 
	jge   .nb111_unroll_loop
	jmp   .nb111_checksingle
.nb111_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb111_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb111_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb111_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb111_iqO]
	mulpd  xmm4, [rsp + nb111_iqH]

	movapd  [rsp + nb111_qqO], xmm3
	movapd  [rsp + nb111_qqH], xmm4
	
	mov rsi, [rbp + nb111_type]
	mov r8d, [rsi + rax*4]
	mov r9d, [rsi + rbx*4]
	mov rsi, [rbp + nb111_vdwparam]
	shl r8d, 1	
	shl r9d, 1	
	mov edi, [rsp + nb111_ntia]
	add r8d, edi
	add r9d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movlpd xmm7, [rsi + r9*8]	;# c6b
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 
	movhpd xmm7, [rsi + r9*8 + 8]	;# c6b c12b 
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb111_c6], xmm4
	movapd [rsp + nb111_c12], xmm6
	
	mov rsi, [rbp + nb111_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move j coordinates to local temp variables 
    movlpd xmm0, [rsi + rax*8] 
    movlpd xmm1, [rsi + rax*8 + 8] 
    movlpd xmm2, [rsi + rax*8 + 16] 
    movhpd xmm0, [rsi + rbx*8] 
    movhpd xmm1, [rsi + rbx*8 + 8] 
    movhpd xmm2, [rsi + rbx*8 + 16] 

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    movapd xmm3, xmm0
    movapd xmm4, xmm1
    movapd xmm5, xmm2
    movapd xmm6, xmm0
    movapd xmm7, xmm1
    movapd xmm8, xmm2
    
    subpd xmm0, [rsp + nb111_ixO]
    subpd xmm1, [rsp + nb111_iyO]
    subpd xmm2, [rsp + nb111_izO]
    subpd xmm3, [rsp + nb111_ixH1]
    subpd xmm4, [rsp + nb111_iyH1]
    subpd xmm5, [rsp + nb111_izH1]
    subpd xmm6, [rsp + nb111_ixH2]
    subpd xmm7, [rsp + nb111_iyH2]
    subpd xmm8, [rsp + nb111_izH2]
    
	movapd [rsp + nb111_dxO], xmm0
	movapd [rsp + nb111_dyO], xmm1
	movapd [rsp + nb111_dzO], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb111_dxH1], xmm3
	movapd [rsp + nb111_dyH1], xmm4
	movapd [rsp + nb111_dzH1], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb111_dxH2], xmm6
	movapd [rsp + nb111_dyH2], xmm7
	movapd [rsp + nb111_dzH2], xmm8
	mulpd  xmm6, xmm6
	mulpd  xmm7, xmm7
	mulpd  xmm8, xmm8
	addpd  xmm0, xmm1
	addpd  xmm0, xmm2
	addpd  xmm3, xmm4
	addpd  xmm3, xmm5
    addpd  xmm6, xmm7
    addpd  xmm6, xmm8

	;# start doing invsqrt for j atoms
    cvtpd2ps xmm1, xmm0
    cvtpd2ps xmm4, xmm3
    cvtpd2ps xmm7, xmm6
	rsqrtps xmm1, xmm1
	rsqrtps xmm4, xmm4
    rsqrtps xmm7, xmm7
    cvtps2pd xmm1, xmm1
    cvtps2pd xmm4, xmm4
    cvtps2pd xmm7, xmm7
	
	movapd  xmm2, xmm1
	movapd  xmm5, xmm4
    movapd  xmm8, xmm7
    
	mulpd   xmm1, xmm1 ;# lu*lu
	mulpd   xmm4, xmm4 ;# lu*lu
    mulpd   xmm7, xmm7 ;# lu*lu
		
	movapd  xmm9, [rsp + nb111_three]
	movapd  xmm10, xmm9
    movapd  xmm11, xmm9

	mulpd   xmm1, xmm0 ;# rsq*lu*lu
	mulpd   xmm4, xmm3 ;# rsq*lu*lu 
    mulpd   xmm7, xmm6 ;# rsq*lu*lu
	
	subpd   xmm9, xmm1
	subpd   xmm10, xmm4
    subpd   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm2
	mulpd   xmm10, xmm5
    mulpd   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb111_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvO
	mulpd   xmm10, xmm15 ;# first iteration for rinvH1
    mulpd   xmm11, xmm15 ;# first iteration for rinvH2

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb111_three]
	movapd  xmm4, xmm1
    movapd  xmm7, xmm1

	mulpd   xmm2, xmm0 ;# rsq*lu*lu
	mulpd   xmm5, xmm3 ;# rsq*lu*lu 
    mulpd   xmm8, xmm6 ;# rsq*lu*lu
	
	subpd   xmm1, xmm2
	subpd   xmm4, xmm5
    subpd   xmm7, xmm8 ;# 3-rsq*lu*lu

	mulpd   xmm9, xmm1
	mulpd   xmm10, xmm4
    mulpd   xmm11, xmm7 ;# lu*(3-rsq*lu*lu)

	movapd  xmm15, [rsp + nb111_half]
	mulpd   xmm9, xmm15  ;#  rinvO 
	mulpd   xmm10, xmm15 ;#   rinvH1
    mulpd   xmm11, xmm15 ;#   rinvH2
	
	;# interactions 
    movapd xmm0, xmm9
    movapd xmm1, xmm10
    movapd xmm2, xmm11
    mulpd  xmm9, xmm9    ;# rinvsq
    mulpd  xmm10, xmm10
    mulpd  xmm11, xmm11
    movapd xmm12, xmm9
    mulpd  xmm12, xmm12 ;# rinv4
    mulpd  xmm12, xmm9  ;# rinv6
    mulpd  xmm0, [rsp + nb111_qqO] 
    mulpd  xmm1, [rsp + nb111_qqH] 
    mulpd  xmm2, [rsp + nb111_qqH] 
    movapd xmm13, xmm12 ;# rinv6
    mulpd  xmm12, xmm12 ;# rinv12
	mulpd  xmm13, [rsp + nb111_c6]
	mulpd  xmm12, [rsp + nb111_c12]
    movapd xmm14, xmm12
    subpd  xmm14, xmm13
    
	addpd  xmm14, [rsp + nb111_Vvdwtot]
	mulpd  xmm13, [rsp + nb111_six]
	mulpd  xmm12, [rsp + nb111_twelve]
	movapd [rsp + nb111_Vvdwtot], xmm14
    subpd  xmm12, xmm13 ;# LJ fscal        

    addpd  xmm12, xmm0
    
    mulpd  xmm9, xmm12
    mulpd  xmm10, xmm1
    mulpd  xmm11, xmm2
    
    addpd xmm0, [rsp + nb111_vctot] 
    addpd xmm1, xmm2
    addpd xmm0, xmm1
    movapd [rsp + nb111_vctot], xmm0
    
    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb111_faction]
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

    movapd xmm7, xmm9
    movapd xmm8, xmm9
    movapd xmm13, xmm11
    movapd xmm14, xmm11
    movapd xmm15, xmm11
    movapd xmm11, xmm10
    movapd xmm12, xmm10

	mulpd xmm7, [rsp + nb111_dxO]
	mulpd xmm8, [rsp + nb111_dyO]
	mulpd xmm9, [rsp + nb111_dzO]
	mulpd xmm10, [rsp + nb111_dxH1]
	mulpd xmm11, [rsp + nb111_dyH1]
	mulpd xmm12, [rsp + nb111_dzH1]
	mulpd xmm13, [rsp + nb111_dxH2]
	mulpd xmm14, [rsp + nb111_dyH2]
	mulpd xmm15, [rsp + nb111_dzH2]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb111_fixO]
    addpd xmm8, [rsp + nb111_fiyO]
    addpd xmm9, [rsp + nb111_fizO]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb111_fixH1]
    addpd xmm11, [rsp + nb111_fiyH1]
    addpd xmm12, [rsp + nb111_fizH1]

    addpd xmm0, xmm13
    addpd xmm1, xmm14
    addpd xmm2, xmm15
    addpd xmm13, [rsp + nb111_fixH2]
    addpd xmm14, [rsp + nb111_fiyH2]
    addpd xmm15, [rsp + nb111_fizH2]

    movapd [rsp + nb111_fixO], xmm7
    movapd [rsp + nb111_fiyO], xmm8
    movapd [rsp + nb111_fizO], xmm9
    movapd [rsp + nb111_fixH1], xmm10
    movapd [rsp + nb111_fiyH1], xmm11
    movapd [rsp + nb111_fizH1], xmm12
    movapd [rsp + nb111_fixH2], xmm13
    movapd [rsp + nb111_fiyH2], xmm14
    movapd [rsp + nb111_fizH2], xmm15
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8],      xmm0
	movlpd [rdi + rax*8 + 8],  xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8],      xmm0
	movhpd [rdi + rbx*8 + 8],  xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb111_innerk],  2
	jl   .nb111_checksingle
	jmp  .nb111_unroll_loop
.nb111_checksingle:	
	mov   edx, [rsp + nb111_innerk]
	and   edx, 1
	jnz  .nb111_dosingle
	jmp  .nb111_updateouterdata
.nb111_dosingle:
	mov   rdx, [rsp + nb111_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb111_innerjjnr],  4	

	mov rsi, [rbp + nb111_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb111_iqO]
	mulpd  xmm4, [rsp + nb111_iqH]

	movapd  [rsp + nb111_qqO], xmm3
	movapd  [rsp + nb111_qqH], xmm4
	
	mov rsi, [rbp + nb111_type]
	mov r8d, [rsi + rax*4]
	mov rsi, [rbp + nb111_vdwparam]
	shl r8d, 1	
	mov edi, [rsp + nb111_ntia]
	add r8d, edi

	movlpd xmm6, [rsi + r8*8]	;# c6a
	movhpd xmm6, [rsi + r8*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movapd [rsp + nb111_c6], xmm4
	movapd [rsp + nb111_c12], xmm6
	
	mov rsi, [rbp + nb111_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 	
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6

	;# calc dr 
	subsd xmm4, [rsp + nb111_ixO]
	subsd xmm5, [rsp + nb111_iyO]
	subsd xmm6, [rsp + nb111_izO]

	;# store dr 
	movapd [rsp + nb111_dxO], xmm4
	movapd [rsp + nb111_dyO], xmm5
	movapd [rsp + nb111_dzO], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move j coords to xmm4-xmm6 
	movapd xmm4, xmm0
	movapd xmm5, xmm1
	movapd xmm6, xmm2

	;# calc dr 
	subsd xmm4, [rsp + nb111_ixH1]
	subsd xmm5, [rsp + nb111_iyH1]
	subsd xmm6, [rsp + nb111_izH1]

	;# store dr 
	movapd [rsp + nb111_dxH1], xmm4
	movapd [rsp + nb111_dyH1], xmm5
	movapd [rsp + nb111_dzH1], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move j coords to xmm3-xmm5
	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	;# calc dr 
	subsd xmm3, [rsp + nb111_ixH2]
	subsd xmm4, [rsp + nb111_iyH2]
	subsd xmm5, [rsp + nb111_izH2]

	;# store dr 
	movapd [rsp + nb111_dxH2], xmm3
	movapd [rsp + nb111_dyH2], xmm4
	movapd [rsp + nb111_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb111_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb111_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb111_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb111_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulsd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movapd xmm1, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm7, [rsp + nb111_qqO]	;# xmm7=vcoul 
	
	mulsd  xmm1, [rsp + nb111_c6]
	mulsd  xmm2, [rsp + nb111_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb111_Vvdwtot]
	mulsd  xmm1, [rsp + nb111_six]
	mulsd  xmm2, [rsp + nb111_twelve]
	subsd  xmm2, xmm1
	addsd  xmm2, xmm7	
	mulsd  xmm4, xmm2	;# total fsO in xmm4 

	addsd  xmm7, [rsp + nb111_vctot]
	
	movsd [rsp + nb111_Vvdwtot], xmm3
	movsd [rsp + nb111_vctot], xmm7

	movapd xmm0, [rsp + nb111_dxO]
	movapd xmm1, [rsp + nb111_dyO]
	movapd xmm2, [rsp + nb111_dzO]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update O forces 
	movapd xmm3, [rsp + nb111_fixO]
	movapd xmm4, [rsp + nb111_fiyO]
	movapd xmm7, [rsp + nb111_fizO]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb111_fixO], xmm3
	movsd [rsp + nb111_fiyO], xmm4
	movsd [rsp + nb111_fizO], xmm7
	;# update j forces with water O 
	movsd [rsp + nb111_fjx], xmm0
	movsd [rsp + nb111_fjy], xmm1
	movsd [rsp + nb111_fjz], xmm2

	;# H1 interactions 
	movapd  xmm4, xmm6	
	mulsd   xmm4, xmm4	;# xmm6=rinv, xmm4=rinvsq 
	mulsd  xmm6, [rsp + nb111_qqH]	;# xmm6=vcoul 
	mulsd  xmm4, xmm6		;# total fsH1 in xmm4 
	
	addsd  xmm6, [rsp + nb111_vctot]

	movapd xmm0, [rsp + nb111_dxH1]
	movapd xmm1, [rsp + nb111_dyH1]
	movapd xmm2, [rsp + nb111_dzH1]
	movsd [rsp + nb111_vctot], xmm6
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb111_fixH1]
	movapd xmm4, [rsp + nb111_fiyH1]
	movapd xmm7, [rsp + nb111_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb111_fixH1], xmm3
	movsd [rsp + nb111_fiyH1], xmm4
	movsd [rsp + nb111_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb111_fjx]
	addsd  xmm1, [rsp + nb111_fjy]
	addsd  xmm2, [rsp + nb111_fjz]
	movsd [rsp + nb111_fjx], xmm0
	movsd [rsp + nb111_fjy], xmm1
	movsd [rsp + nb111_fjz], xmm2

	;# H2 interactions 
	movapd  xmm4, xmm5	
	mulsd   xmm4, xmm4	;# xmm5=rinv, xmm4=rinvsq 
	mulsd  xmm5, [rsp + nb111_qqH]	;# xmm5=vcoul 
	mulsd  xmm4, xmm5		;# total fsH1 in xmm4 
	
	addsd  xmm5, [rsp + nb111_vctot]

	movapd xmm0, [rsp + nb111_dxH2]
	movapd xmm1, [rsp + nb111_dyH2]
	movapd xmm2, [rsp + nb111_dzH2]
	movsd [rsp + nb111_vctot], xmm5
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb111_fixH2]
	movapd xmm4, [rsp + nb111_fiyH2]
	movapd xmm7, [rsp + nb111_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movsd [rsp + nb111_fixH2], xmm3
	movsd [rsp + nb111_fiyH2], xmm4
	movsd [rsp + nb111_fizH2], xmm7

	mov rdi, [rbp + nb111_faction]
	;# update j forces 
	addsd  xmm0, [rsp + nb111_fjx]
	addsd  xmm1, [rsp + nb111_fjy]
	addsd  xmm2, [rsp + nb111_fjz]
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb111_updateouterdata:
	mov   ecx, [rsp + nb111_ii3]
	mov   rdi, [rbp + nb111_faction]
	mov   rsi, [rbp + nb111_fshift]
	mov   edx, [rsp + nb111_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb111_fixO]
	movapd xmm1, [rsp + nb111_fiyO]
	movapd xmm2, [rsp + nb111_fizO]

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

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6,xmm1 

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb111_fixH1]
	movapd xmm1, [rsp + nb111_fiyH1]
	movapd xmm2, [rsp + nb111_fizH1]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 24]
	movsd  xmm4, [rdi + rcx*8 + 32]
	movsd  xmm5, [rdi + rcx*8 + 40]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 24], xmm3
	movsd  [rdi + rcx*8 + 32], xmm4
	movsd  [rdi + rcx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb111_fixH2]
	movapd xmm1, [rsp + nb111_fiyH2]
	movapd xmm2, [rsp + nb111_fizH2]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 48]
	movsd  xmm4, [rdi + rcx*8 + 56]
	movsd  xmm5, [rdi + rcx*8 + 64]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 48], xmm3
	movsd  [rdi + rcx*8 + 56], xmm4
	movsd  [rdi + rcx*8 + 64], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addsd xmm7, xmm2
	unpcklpd xmm0, xmm1
	addpd xmm6, xmm0

	;# increment fshift force 
	movlpd xmm3, [rsi + rdx*8]
	movhpd xmm3, [rsi + rdx*8 + 8]
	movsd  xmm4, [rsi + rdx*8 + 16]
	subpd  xmm3, xmm6
	subsd  xmm4, xmm7
	movlpd [rsi + rdx*8],      xmm3
	movhpd [rsi + rdx*8 + 8],  xmm3
	movsd  [rsi + rdx*8 + 16], xmm4

	;# get n from stack
	mov esi, [rsp + nb111_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb111_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb111_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb111_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb111_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb111_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
       ;# finish if last 
        mov ecx, [rsp + nb111_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb111_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb111_n], esi
        jmp .nb111_outer
.nb111_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb111_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb111_end
        ;# non-zero, do one more workunit
        jmp   .nb111_threadloop
.nb111_end:
	mov eax, [rsp + nb111_nouter]
	mov ebx, [rsp + nb111_ninner]
	mov rcx, [rbp + nb111_outeriter]
	mov rdx, [rbp + nb111_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 784
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





.globl nb_kernel111nf_x86_64_sse2
.globl _nb_kernel111nf_x86_64_sse2
nb_kernel111nf_x86_64_sse2:	
_nb_kernel111nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb111nf_fshift,         16
.equiv          nb111nf_gid,            24
.equiv          nb111nf_pos,            32
.equiv          nb111nf_faction,        40
.equiv          nb111nf_charge,         48
.equiv          nb111nf_p_facel,        56
.equiv          nb111nf_argkrf,         64
.equiv          nb111nf_argcrf,         72
.equiv          nb111nf_Vc,             80
.equiv          nb111nf_type,           88
.equiv          nb111nf_p_ntype,        96
.equiv          nb111nf_vdwparam,       104
.equiv          nb111nf_Vvdw,           112
.equiv          nb111nf_p_tabscale,     120
.equiv          nb111nf_VFtab,          128
.equiv          nb111nf_invsqrta,       136
.equiv          nb111nf_dvda,           144
.equiv          nb111nf_p_gbtabscale,   152
.equiv          nb111nf_GBtab,          160
.equiv          nb111nf_p_nthreads,     168
.equiv          nb111nf_count,          176
.equiv          nb111nf_mtx,            184
.equiv          nb111nf_outeriter,      192
.equiv          nb111nf_inneriter,      200
.equiv          nb111nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb111nf_ixO,            0
.equiv          nb111nf_iyO,            16
.equiv          nb111nf_izO,            32
.equiv          nb111nf_ixH1,           48
.equiv          nb111nf_iyH1,           64
.equiv          nb111nf_izH1,           80
.equiv          nb111nf_ixH2,           96
.equiv          nb111nf_iyH2,           112
.equiv          nb111nf_izH2,           128
.equiv          nb111nf_iqO,            144
.equiv          nb111nf_iqH,            160
.equiv          nb111nf_qqO,            176
.equiv          nb111nf_qqH,            192
.equiv          nb111nf_c6,             208
.equiv          nb111nf_c12,            224
.equiv          nb111nf_vctot,          240
.equiv          nb111nf_Vvdwtot,        256
.equiv          nb111nf_half,           272
.equiv          nb111nf_three,          288
.equiv          nb111nf_is3,            304
.equiv          nb111nf_ii3,            308
.equiv          nb111nf_nri,            312
.equiv          nb111nf_iinr,           320
.equiv          nb111nf_jindex,         328
.equiv          nb111nf_jjnr,           336
.equiv          nb111nf_shift,          344
.equiv          nb111nf_shiftvec,       352
.equiv          nb111nf_facel,          360
.equiv          nb111nf_innerjjnr,      368
.equiv          nb111nf_ntia,           376
.equiv          nb111nf_innerk,         380
.equiv          nb111nf_n,              384
.equiv          nb111nf_nn1,            388
.equiv          nb111nf_nouter,         392
.equiv          nb111nf_ninner,         396

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
	sub rsp, 400		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb111nf_nouter], eax
	mov [rsp + nb111nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb111nf_nri], edi
	mov [rsp + nb111nf_iinr], rsi
	mov [rsp + nb111nf_jindex], rdx
	mov [rsp + nb111nf_jjnr], rcx
	mov [rsp + nb111nf_shift], r8
	mov [rsp + nb111nf_shiftvec], r9
	mov rsi, [rbp + nb111nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb111nf_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb111nf_half], eax
	mov [rsp + nb111nf_half+4], ebx
	movsd xmm1, [rsp + nb111nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb111nf_half], xmm1
	movapd [rsp + nb111nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb111nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb111nf_charge]
	movsd xmm3, [rdx + rbx*8]	
	movsd xmm4, [rdx + rbx*8 + 8]	
	mov rsi, [rbp + nb111nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb111nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb111nf_iqO], xmm3
	movapd [rsp + nb111nf_iqH], xmm4
	
	mov   rdx, [rbp + nb111nf_type]
	mov   ecx, [rdx + rbx*4]
	shl   ecx, 1
	mov rdi, [rbp + nb111nf_p_ntype]
	imul  ecx, [rdi]      ;# rcx = ntia = 2*ntype*type[ii0] 
	mov   [rsp + nb111nf_ntia], ecx		
.nb111nf_threadloop:
        mov   rsi, [rbp + nb111nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb111nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb111nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb111nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb111nf_n], eax
        mov [rsp + nb111nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb111nf_outerstart
        jmp .nb111nf_end

.nb111nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb111nf_nouter]
	mov [rsp + nb111nf_nouter], ebx

.nb111nf_outer:
	mov   rax, [rsp + nb111nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb111nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb111nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb111nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb111nf_ii3], ebx

	addsd xmm3, [rax + rbx*8]
	addsd xmm4, [rax + rbx*8 + 8]
	addsd xmm5, [rax + rbx*8 + 16]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb111nf_ixO], xmm3
	movapd [rsp + nb111nf_iyO], xmm4
	movapd [rsp + nb111nf_izO], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 24]
	addsd xmm1, [rax + rbx*8 + 32]
	addsd xmm2, [rax + rbx*8 + 40]		
	addsd xmm3, [rax + rbx*8 + 48]
	addsd xmm4, [rax + rbx*8 + 56]
	addsd xmm5, [rax + rbx*8 + 64]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb111nf_ixH1], xmm0
	movapd [rsp + nb111nf_iyH1], xmm1
	movapd [rsp + nb111nf_izH1], xmm2
	movapd [rsp + nb111nf_ixH2], xmm3
	movapd [rsp + nb111nf_iyH2], xmm4
	movapd [rsp + nb111nf_izH2], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb111nf_vctot], xmm4
	movapd [rsp + nb111nf_Vvdwtot], xmm4
	
	mov   rax, [rsp + nb111nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb111nf_pos]
	mov   rax, [rsp + nb111nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb111nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb111nf_ninner]
	mov   [rsp + nb111nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb111nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb111nf_unroll_loop
	jmp   .nb111nf_checksingle
.nb111nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb111nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb111nf_innerjjnr],  8	;# advance pointer (unrolled 2) 

	mov rsi, [rbp + nb111nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb111nf_iqO]
	mulpd  xmm4, [rsp + nb111nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx

	movapd  [rsp + nb111nf_qqO], xmm3
	movapd  [rsp + nb111nf_qqH], xmm4
	
	mov rsi, [rbp + nb111nf_type]
	mov eax, [rsi + rax*4]
	mov ebx, [rsi + rbx*4]
	mov rsi, [rbp + nb111nf_vdwparam]
	shl eax, 1	
	shl ebx, 1	
	mov edi, [rsp + nb111nf_ntia]
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
	movapd [rsp + nb111nf_c6], xmm4
	movapd [rsp + nb111nf_c12], xmm6
	
	mov rsi, [rbp + nb111nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb111nf_ixO]
	movapd xmm5, [rsp + nb111nf_iyO]
	movapd xmm6, [rsp + nb111nf_izO]

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
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb111nf_ixH1]
	movapd xmm5, [rsp + nb111nf_iyH1]
	movapd xmm6, [rsp + nb111nf_izH1]

	;# calc dr 
	subpd xmm4, xmm0
	subpd xmm5, xmm1
	subpd xmm6, xmm2

	;# square it 
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	mulpd xmm6,xmm6
	addpd xmm6, xmm5
	addpd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [rsp + nb111nf_ixH2]
	movapd xmm4, [rsp + nb111nf_iyH2]
	movapd xmm5, [rsp + nb111nf_izH2]

	;# calc dr 
	subpd xmm3, xmm0
	subpd xmm4, xmm1
	subpd xmm5, xmm2

	;# square it 
	mulpd xmm3,xmm3
	mulpd xmm4,xmm4
	mulpd xmm5,xmm5
	addpd xmm5, xmm4
	addpd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111nf_three]
	subpd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb111nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111nf_three]
	subpd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb111nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111nf_three]
	subpd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb111nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulpd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movapd xmm1, xmm4
	mulpd  xmm1, xmm4
	mulpd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulpd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulpd  xmm7, [rsp + nb111nf_qqO]	;# xmm7=vcoul 
	
	mulpd  xmm1, [rsp + nb111nf_c6]
	mulpd  xmm2, [rsp + nb111nf_c12]
	movapd xmm3, xmm2
	subpd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addpd  xmm3, [rsp + nb111nf_Vvdwtot]
	addpd  xmm7, [rsp + nb111nf_vctot]	
	movapd [rsp + nb111nf_Vvdwtot], xmm3
	movapd [rsp + nb111nf_vctot], xmm7

	;# H1 interactions 
	mulpd  xmm6, [rsp + nb111nf_qqH]	;# xmm6=vcoul 
	addpd  xmm6, [rsp + nb111nf_vctot]
	movapd [rsp + nb111nf_vctot], xmm6

	;# H2 interactions 
	mulpd  xmm5, [rsp + nb111nf_qqH]	;# xmm5=vcoul 
	addpd  xmm5, [rsp + nb111nf_vctot]
	movapd [rsp + nb111nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb111nf_innerk],  2
	jl    .nb111nf_checksingle
	jmp   .nb111nf_unroll_loop
.nb111nf_checksingle:	
	mov   edx, [rsp + nb111nf_innerk]
	and   edx, 1
	jnz   .nb111nf_dosingle
	jmp   .nb111nf_updateouterdata
.nb111nf_dosingle:
	mov   rdx, [rsp + nb111nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb111nf_innerjjnr],  4	

	mov rsi, [rbp + nb111nf_charge]    ;# base of charge[] 

	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3
	mulpd  xmm3, [rsp + nb111nf_iqO]
	mulpd  xmm4, [rsp + nb111nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 

	movapd  [rsp + nb111nf_qqO], xmm3
	movapd  [rsp + nb111nf_qqH], xmm4
	
	mov rsi, [rbp + nb111nf_type]
	mov eax, [rsi + rax*4]
	mov rsi, [rbp + nb111nf_vdwparam]
	shl eax, 1	
	mov edi, [rsp + nb111nf_ntia]
	add eax, edi

	movlpd xmm6, [rsi + rax*8]	;# c6a
	movhpd xmm6, [rsi + rax*8 + 8]	;# c6a c12a 

	xorpd xmm7, xmm7
	movapd xmm4, xmm6
	unpcklpd xmm4, xmm7
	unpckhpd xmm6, xmm7
	
	movd  eax, mm0
	movd  ebx, mm1
	movapd [rsp + nb111nf_c6], xmm4
	movapd [rsp + nb111nf_c12], xmm6
	
	mov rsi, [rbp + nb111nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 

	;# move coordinates to xmm0-xmm2 
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixO-izO to xmm4-xmm6 
	movapd xmm4, [rsp + nb111nf_ixO]
	movapd xmm5, [rsp + nb111nf_iyO]
	movapd xmm6, [rsp + nb111nf_izO]

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
	movapd xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb111nf_ixH1]
	movapd xmm5, [rsp + nb111nf_iyH1]
	movapd xmm6, [rsp + nb111nf_izH1]

	;# calc dr 
	subsd xmm4, xmm0
	subsd xmm5, xmm1
	subsd xmm6, xmm2

	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm6, xmm5
	addsd xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixH2-izH2 to xmm3-xmm5  
	movapd xmm3, [rsp + nb111nf_ixH2]
	movapd xmm4, [rsp + nb111nf_iyH2]
	movapd xmm5, [rsp + nb111nf_izH2]

	;# calc dr 
	subsd xmm3, xmm0
	subsd xmm4, xmm1
	subsd xmm5, xmm2

	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm7, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111nf_three]
	subsd xmm4, xmm7	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb111nf_half] ;# rinv 
	movapd  xmm7, xmm4	;# rinvO in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm6, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111nf_three]
	subsd xmm4, xmm6	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb111nf_half] ;# rinv 
	movapd  xmm6, xmm4	;# rinvH1 in xmm6 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb111nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb111nf_half] ;# iter1 ( new lu) 

	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm5, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb111nf_three]
	subsd xmm4, xmm5	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb111nf_half] ;# rinv 
	movapd  xmm5, xmm4	;# rinvH2 in xmm5 

	;# do O interactions 
	movapd  xmm4, xmm7	
	mulsd   xmm4, xmm4	;# xmm7=rinv, xmm4=rinvsq 
	movapd xmm1, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm1, xmm4	;# xmm1=rinvsix 
	movapd xmm2, xmm1
	mulsd  xmm2, xmm2	;# xmm2=rinvtwelve 
	mulsd  xmm7, [rsp + nb111nf_qqO]	;# xmm7=vcoul 
	
	mulsd  xmm1, [rsp + nb111nf_c6]
	mulsd  xmm2, [rsp + nb111nf_c12]
	movapd xmm3, xmm2
	subsd  xmm3, xmm1	;# Vvdw=Vvdw12-Vvdw6 		
	addsd  xmm3, [rsp + nb111nf_Vvdwtot]
	addsd  xmm7, [rsp + nb111nf_vctot]
	movsd [rsp + nb111nf_Vvdwtot], xmm3
	movsd [rsp + nb111nf_vctot], xmm7

	;# H1 interactions 
	mulsd  xmm6, [rsp + nb111nf_qqH]	;# xmm6=vcoul 
	addsd  xmm6, [rsp + nb111nf_vctot]
	movsd [rsp + nb111nf_vctot], xmm6

	;# H2 interactions 
	mulsd  xmm5, [rsp + nb111nf_qqH]	;# xmm5=vcoul 
	addsd  xmm5, [rsp + nb111nf_vctot]
	movsd [rsp + nb111nf_vctot], xmm5
	
.nb111nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb111nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb111nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb111nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb111nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
	;# accumulate total lj energy and update it 
	movapd xmm7, [rsp + nb111nf_Vvdwtot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 

	;# add earlier value from mem 
	mov   rax, [rbp + nb111nf_Vvdw]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb111nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb111nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb111nf_n], esi
        jmp .nb111nf_outer
.nb111nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb111nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb111nf_end
        ;# non-zero, do one more workunit
        jmp   .nb111nf_threadloop
.nb111nf_end:
	mov eax, [rsp + nb111nf_nouter]
	mov ebx, [rsp + nb111nf_ninner]
	mov rcx, [rbp + nb111nf_outeriter]
	mov rdx, [rbp + nb111nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 400
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
