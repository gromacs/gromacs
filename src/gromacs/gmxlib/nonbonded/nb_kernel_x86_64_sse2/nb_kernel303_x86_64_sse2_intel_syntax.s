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
	
.globl nb_kernel303_x86_64_sse2
.globl _nb_kernel303_x86_64_sse2
nb_kernel303_x86_64_sse2:	
_nb_kernel303_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb303_fshift,           16
.equiv          nb303_gid,              24
.equiv          nb303_pos,              32
.equiv          nb303_faction,          40
.equiv          nb303_charge,           48
.equiv          nb303_p_facel,          56
.equiv          nb303_argkrf,           64
.equiv          nb303_argcrf,           72
.equiv          nb303_Vc,               80
.equiv          nb303_type,             88
.equiv          nb303_p_ntype,          96
.equiv          nb303_vdwparam,         104
.equiv          nb303_Vvdw,             112
.equiv          nb303_p_tabscale,       120
.equiv          nb303_VFtab,            128
.equiv          nb303_invsqrta,         136
.equiv          nb303_dvda,             144
.equiv          nb303_p_gbtabscale,     152
.equiv          nb303_GBtab,            160
.equiv          nb303_p_nthreads,       168
.equiv          nb303_count,            176
.equiv          nb303_mtx,              184
.equiv          nb303_outeriter,        192
.equiv          nb303_inneriter,        200
.equiv          nb303_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse2 use 
.equiv          nb303_ixM,              0
.equiv          nb303_iyM,              16
.equiv          nb303_izM,              32
.equiv          nb303_ixH1,             48
.equiv          nb303_iyH1,             64
.equiv          nb303_izH1,             80
.equiv          nb303_ixH2,             96
.equiv          nb303_iyH2,             112
.equiv          nb303_izH2,             128
.equiv          nb303_iqM,              144
.equiv          nb303_iqH,              160
.equiv          nb303_dxM,              176
.equiv          nb303_dyM,              192
.equiv          nb303_dzM,              208
.equiv          nb303_dxH1,             224
.equiv          nb303_dyH1,             240
.equiv          nb303_dzH1,             256
.equiv          nb303_dxH2,             272
.equiv          nb303_dyH2,             288
.equiv          nb303_dzH2,             304
.equiv          nb303_qqM,              320
.equiv          nb303_qqH,              336
.equiv          nb303_rinvM,            352
.equiv          nb303_rinvH1,           368
.equiv          nb303_rinvH2,           384
.equiv          nb303_rM,               400
.equiv          nb303_rH1,              416
.equiv          nb303_rH2,              432
.equiv          nb303_tsc,              448
.equiv          nb303_two,              464
.equiv          nb303_vctot,            480
.equiv          nb303_fixM,             496
.equiv          nb303_fiyM,             512
.equiv          nb303_fizM,             528
.equiv          nb303_fixH1,            544
.equiv          nb303_fiyH1,            560
.equiv          nb303_fizH1,            576
.equiv          nb303_fixH2,            592
.equiv          nb303_fiyH2,            608
.equiv          nb303_fizH2,            624
.equiv          nb303_fjx,              640
.equiv          nb303_fjy,              656
.equiv          nb303_fjz,              672
.equiv          nb303_half,             688
.equiv          nb303_three,            704
.equiv          nb303_is3,              720
.equiv          nb303_ii3,              724
.equiv          nb303_nri,              728
.equiv          nb303_iinr,             736
.equiv          nb303_jindex,           744
.equiv          nb303_jjnr,             752
.equiv          nb303_shift,            760
.equiv          nb303_shiftvec,         768
.equiv          nb303_facel,            776
.equiv          nb303_innerjjnr,        784
.equiv          nb303_innerk,           792
.equiv          nb303_n,                796
.equiv          nb303_nn1,              800
.equiv          nb303_nouter,           804
.equiv          nb303_ninner,           808

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
	sub rsp, 816		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb303_nouter], eax
	mov [rsp + nb303_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb303_nri], edi
	mov [rsp + nb303_iinr], rsi
	mov [rsp + nb303_jindex], rdx
	mov [rsp + nb303_jjnr], rcx
	mov [rsp + nb303_shift], r8
	mov [rsp + nb303_shiftvec], r9
	mov rsi, [rbp + nb303_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb303_facel], xmm0

	mov rax, [rbp + nb303_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb303_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb303_half], eax
	mov [rsp + nb303_half+4], ebx
	movsd xmm1, [rsp + nb303_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb303_half], xmm1
	movapd [rsp + nb303_two], xmm2
	movapd [rsp + nb303_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb303_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb303_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	
	mov rsi, [rbp + nb303_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb303_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb303_iqH], xmm3
	movapd [rsp + nb303_iqM], xmm4
	
.nb303_threadloop:
        mov   rsi, [rbp + nb303_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb303_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb303_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb303_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb303_n], eax
        mov [rsp + nb303_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb303_outerstart
        jmp .nb303_end

.nb303_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb303_nouter]
	mov [rsp + nb303_nouter], ebx

.nb303_outer:
	mov   rax, [rsp + nb303_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb303_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb303_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb303_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb303_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb303_ii3], ebx

	addsd xmm3, [rax + rbx*8 + 24]
	addsd xmm4, [rax + rbx*8 + 32]
	addsd xmm5, [rax + rbx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb303_ixH1], xmm3
	movapd [rsp + nb303_iyH1], xmm4
	movapd [rsp + nb303_izH1], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 48]
	addsd xmm1, [rax + rbx*8 + 56]
	addsd xmm2, [rax + rbx*8 + 64]		
	addsd xmm3, [rax + rbx*8 + 72]
	addsd xmm4, [rax + rbx*8 + 80]
	addsd xmm5, [rax + rbx*8 + 88]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb303_ixH2], xmm0
	movapd [rsp + nb303_iyH2], xmm1
	movapd [rsp + nb303_izH2], xmm2
	movapd [rsp + nb303_ixM], xmm3
	movapd [rsp + nb303_iyM], xmm4
	movapd [rsp + nb303_izM], xmm5
	
	;# clear vctot and i forces 
	xorpd xmm4, xmm4
	movapd [rsp + nb303_vctot], xmm4
	movapd [rsp + nb303_fixM], xmm4
	movapd [rsp + nb303_fiyM], xmm4
	movapd [rsp + nb303_fizM], xmm4
	movapd [rsp + nb303_fixH1], xmm4
	movapd [rsp + nb303_fiyH1], xmm4
	movapd [rsp + nb303_fizH1], xmm4
	movapd [rsp + nb303_fixH2], xmm4
	movapd [rsp + nb303_fiyH2], xmm4
	movapd [rsp + nb303_fizH2], xmm4
	
	mov   rax, [rsp + nb303_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb303_pos]
	mov   rdi, [rbp + nb303_faction]	
	mov   rax, [rsp + nb303_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb303_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb303_ninner]
	mov   [rsp + nb303_ninner], ecx
	add   edx, 0
	mov   [rsp + nb303_innerk], edx    ;# number of innerloop atoms 
	jge   .nb303_unroll_loop
	jmp   .nb303_checksingle
.nb303_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb303_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb303_innerjjnr],  8	;# advance pointer (unrolled 2) 
	mov rsi, [rbp + nb303_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb303_iqM]
	mulpd  xmm4, [rsp + nb303_iqH]

	movapd  [rsp + nb303_qqM], xmm3
	movapd  [rsp + nb303_qqH], xmm4	

	mov rsi, [rbp + nb303_pos]       ;# base of pos[] 

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
    
    subpd xmm0, [rsp + nb303_ixH1]
    subpd xmm1, [rsp + nb303_iyH1]
    subpd xmm2, [rsp + nb303_izH1]
    subpd xmm3, [rsp + nb303_ixH2]
    subpd xmm4, [rsp + nb303_iyH2]
    subpd xmm5, [rsp + nb303_izH2]
    subpd xmm6, [rsp + nb303_ixM]
    subpd xmm7, [rsp + nb303_iyM]
    subpd xmm8, [rsp + nb303_izM]
    
	movapd [rsp + nb303_dxH1], xmm0
	movapd [rsp + nb303_dyH1], xmm1
	movapd [rsp + nb303_dzH1], xmm2
	mulpd  xmm0, xmm0
	mulpd  xmm1, xmm1
	mulpd  xmm2, xmm2
	movapd [rsp + nb303_dxH2], xmm3
	movapd [rsp + nb303_dyH2], xmm4
	movapd [rsp + nb303_dzH2], xmm5
	mulpd  xmm3, xmm3
	mulpd  xmm4, xmm4
	mulpd  xmm5, xmm5
	movapd [rsp + nb303_dxM], xmm6
	movapd [rsp + nb303_dyM], xmm7
	movapd [rsp + nb303_dzM], xmm8
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
		
	movapd  xmm9, [rsp + nb303_three]
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

	movapd  xmm15, [rsp + nb303_half]
	mulpd   xmm9, xmm15  ;# first iteration for rinvH1
	mulpd   xmm10, xmm15 ;# first iteration for rinvH2
    mulpd   xmm11, xmm15 ;# first iteration for rinvM

    ;# second iteration step    
	movapd  xmm2, xmm9
	movapd  xmm5, xmm10
    movapd  xmm8, xmm11
    
	mulpd   xmm2, xmm2 ;# lu*lu
	mulpd   xmm5, xmm5 ;# lu*lu
    mulpd   xmm8, xmm8 ;# lu*lu
		
	movapd  xmm1, [rsp + nb303_three]
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

	movapd  xmm15, [rsp + nb303_half]
	mulpd   xmm9, xmm15  ;#  rinvH1
	mulpd   xmm10, xmm15 ;#   rinvH2
    mulpd   xmm11, xmm15 ;#   rinvM
	
	movapd  [rsp + nb303_rinvH1], xmm9
	movapd  [rsp + nb303_rinvH2], xmm10
	movapd  [rsp + nb303_rinvM], xmm11
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movapd xmm1, [rsp + nb303_tsc]
    mulpd  xmm0, xmm9  ;# r
    mulpd  xmm3, xmm10
    mulpd  xmm6, xmm11
    mulpd  xmm0, xmm1 ;# rtab
    mulpd  xmm3, xmm1
    mulpd  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttpd2dq xmm1, xmm0
    cvttpd2dq xmm4, xmm3
    cvttpd2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2pd  xmm2, xmm1
    cvtdq2pd  xmm5, xmm4
    cvtdq2pd  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    pshufd xmm13, xmm1, 1
    pshufd xmm14, xmm4, 1
    pshufd xmm15, xmm7, 1
    movd    r8d, xmm1
    movd    r10d, xmm4
    movd    r12d, xmm7
    movd    r9d, xmm13
    movd    r11d, xmm14
    movd    r13d, xmm15
        
    mov  rsi, [rbp + nb303_VFtab]

    ;# calculate eps
    subpd     xmm0, xmm2
    subpd     xmm3, xmm5
    subpd     xmm6, xmm8

    movapd    xmm12, xmm0  ;# epsH1
    movapd    xmm13, xmm3  ;# epsH2
    movapd    xmm14, xmm6  ;# epsM
     
    ;# Load LOTS of table data
    movlpd xmm0,  [rsi + r8*8]
    movlpd xmm1,  [rsi + r8*8 + 8]
    movlpd xmm2,  [rsi + r8*8 + 16]
    movlpd xmm3,  [rsi + r8*8 + 24]
    movlpd xmm4,  [rsi + r10*8]
    movlpd xmm5,  [rsi + r10*8 + 8]
    movlpd xmm6,  [rsi + r10*8 + 16]
    movlpd xmm7,  [rsi + r10*8 + 24]
    movlpd xmm8,  [rsi + r12*8]
    movlpd xmm9,  [rsi + r12*8 + 8]
    movlpd xmm10, [rsi + r12*8 + 16]
    movlpd xmm11, [rsi + r12*8 + 24]
    movhpd xmm0,  [rsi + r9*8]
    movhpd xmm1,  [rsi + r9*8 + 8]
    movhpd xmm2,  [rsi + r9*8 + 16]
    movhpd xmm3,  [rsi + r9*8 + 24]
    movhpd xmm4,  [rsi + r11*8]
    movhpd xmm5,  [rsi + r11*8 + 8]
    movhpd xmm6,  [rsi + r11*8 + 16]
    movhpd xmm7,  [rsi + r11*8 + 24]
    movhpd xmm8,  [rsi + r13*8]
    movhpd xmm9,  [rsi + r13*8 + 8]
    movhpd xmm10, [rsi + r13*8 + 16]
    movhpd xmm11, [rsi + r13*8 + 24]
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    mulpd  xmm3, xmm12   ;# Heps
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 
    mulpd  xmm2, xmm12   ;# Geps
    mulpd  xmm6, xmm13
    mulpd  xmm10, xmm14 
    mulpd  xmm3, xmm12   ;# Heps2
    mulpd  xmm7, xmm13
    mulpd  xmm11, xmm14 

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
    addpd  xmm3, xmm1   ;# FF = Fp + 2*Heps2 + Geps
    addpd  xmm7, xmm5
    addpd  xmm11, xmm9
    mulpd  xmm1, xmm12   ;# eps*Fp
    mulpd  xmm5, xmm13
    mulpd  xmm9, xmm14
    movapd xmm12, [rsp + nb303_qqH]
    movapd xmm13, [rsp + nb303_qqM]
    addpd  xmm1, xmm0     ;# VV
    addpd  xmm5, xmm4
    addpd  xmm9, xmm8
    mulpd  xmm1, xmm12   ;# VV*qq = vcoul
    mulpd  xmm5, xmm12
    mulpd  xmm9, xmm13
    mulpd  xmm3, xmm12    ;# FF*qq = fij
    mulpd  xmm7, xmm12
    mulpd  xmm11, xmm13
    
    ;# accumulate vctot
    addpd  xmm1, [rsp + nb303_vctot]
    addpd  xmm5, xmm9
    addpd  xmm1, xmm5
    movapd [rsp + nb303_vctot], xmm1
    
    movapd xmm10, [rsp + nb303_tsc]
    mulpd  xmm3, xmm10  ;# fscal
    mulpd  xmm7, xmm10
    mulpd  xmm10, xmm11
    
    xorpd xmm4, xmm4
    xorpd xmm8, xmm8
    xorpd xmm11, xmm11
    
    subpd xmm4, xmm3
    subpd xmm8, xmm7
    subpd xmm11, xmm10
 
    mulpd xmm4, [rsp + nb303_rinvH1]
    mulpd xmm8, [rsp + nb303_rinvH2]
    mulpd xmm11, [rsp + nb303_rinvM]
    
    ;# move j forces to xmm0-xmm2
    mov rdi, [rbp + nb303_faction]
	movlpd xmm0, [rdi + rax*8]
	movlpd xmm1, [rdi + rax*8 + 8]
	movlpd xmm2, [rdi + rax*8 + 16]
	movhpd xmm0, [rdi + rbx*8]
	movhpd xmm1, [rdi + rbx*8 + 8]
	movhpd xmm2, [rdi + rbx*8 + 16]

    movapd xmm3, xmm4
    movapd xmm5, xmm4
    movapd xmm7, xmm8
    movapd xmm9, xmm8
    movapd xmm10, xmm11
    movapd xmm12, xmm11

	mulpd xmm3, [rsp + nb303_dxH1]
	mulpd xmm4, [rsp + nb303_dyH1]
	mulpd xmm5, [rsp + nb303_dzH1]
	mulpd xmm7, [rsp + nb303_dxH2]
	mulpd xmm8, [rsp + nb303_dyH2]
	mulpd xmm9, [rsp + nb303_dzH2]
	mulpd xmm10, [rsp + nb303_dxM]
	mulpd xmm11, [rsp + nb303_dyM]
	mulpd xmm12, [rsp + nb303_dzM]

    addpd xmm0, xmm3
    addpd xmm1, xmm4
    addpd xmm2, xmm5
    addpd xmm3, [rsp + nb303_fixH1]
    addpd xmm4, [rsp + nb303_fiyH1]
    addpd xmm5, [rsp + nb303_fizH1]

    addpd xmm0, xmm7
    addpd xmm1, xmm8
    addpd xmm2, xmm9
    addpd xmm7, [rsp + nb303_fixH2]
    addpd xmm8, [rsp + nb303_fiyH2]
    addpd xmm9, [rsp + nb303_fizH2]

    addpd xmm0, xmm10
    addpd xmm1, xmm11
    addpd xmm2, xmm12
    addpd xmm10, [rsp + nb303_fixM]
    addpd xmm11, [rsp + nb303_fiyM]
    addpd xmm12, [rsp + nb303_fizM]

    movapd [rsp + nb303_fixH1], xmm3
    movapd [rsp + nb303_fiyH1], xmm4
    movapd [rsp + nb303_fizH1], xmm5
    movapd [rsp + nb303_fixH2], xmm7
    movapd [rsp + nb303_fiyH2], xmm8
    movapd [rsp + nb303_fizH2], xmm9
    movapd [rsp + nb303_fixM], xmm10
    movapd [rsp + nb303_fiyM], xmm11
    movapd [rsp + nb303_fizM], xmm12
   
    ;# store back j forces from xmm0-xmm2
	movlpd [rdi + rax*8], xmm0
	movlpd [rdi + rax*8 + 8], xmm1
	movlpd [rdi + rax*8 + 16], xmm2
	movhpd [rdi + rbx*8], xmm0
	movhpd [rdi + rbx*8 + 8], xmm1
	movhpd [rdi + rbx*8 + 16], xmm2

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb303_innerk],  2
	jl    .nb303_checksingle
	jmp   .nb303_unroll_loop
.nb303_checksingle:	
	mov   edx, [rsp + nb303_innerk]
	and   edx, 1
	jnz   .nb303_dosingle
	jmp   .nb303_updateouterdata
.nb303_dosingle:
	mov   rdx, [rsp + nb303_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb303_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb303_iqM]
	mulpd  xmm4, [rsp + nb303_iqH]

	movapd  [rsp + nb303_qqM], xmm3
	movapd  [rsp + nb303_qqH], xmm4	

	mov rsi, [rbp + nb303_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	;# move coordinates to xmm4-xmm6 & xmm0-xmm2 	
	movlpd xmm4, [rsi + rax*8]
	movlpd xmm5, [rsi + rax*8 + 8]
	movlpd xmm6, [rsi + rax*8 + 16]
    movapd xmm0, xmm4
    movapd xmm1, xmm5
    movapd xmm2, xmm6
    
	;# calc dr 
	subsd xmm4, [rsp + nb303_ixM]
	subsd xmm5, [rsp + nb303_iyM]
	subsd xmm6, [rsp + nb303_izM]

	;# store dr 
	movapd [rsp + nb303_dxM], xmm4
	movapd [rsp + nb303_dyM], xmm5
	movapd [rsp + nb303_dzM], xmm6
	;# square it 
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	mulsd xmm6,xmm6
	addsd xmm4, xmm5
	addsd xmm4, xmm6
	movapd xmm7, xmm4
	;# rsqM in xmm7 

	;# move j coords to xmm4-xmm6 
	movapd xmm4, xmm0
	movapd xmm5, xmm1
	movapd xmm6, xmm2

	;# calc dr 
	subsd xmm4, [rsp + nb303_ixH1]
	subsd xmm5, [rsp + nb303_iyH1]
	subsd xmm6, [rsp + nb303_izH1]

	;# store dr 
	movapd [rsp + nb303_dxH1], xmm4
	movapd [rsp + nb303_dyH1], xmm5
	movapd [rsp + nb303_dzH1], xmm6
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
	subsd xmm3, [rsp + nb303_ixH2]
	subsd xmm4, [rsp + nb303_iyH2]
	subsd xmm5, [rsp + nb303_izH2]

	;# store dr 
	movapd [rsp + nb303_dxH2], xmm3
	movapd [rsp + nb303_dyH2], xmm4
	movapd [rsp + nb303_dzH2], xmm5
	;# square it 
	mulsd xmm3,xmm3
	mulsd xmm4,xmm4
	mulsd xmm5,xmm5
	addsd xmm5, xmm4
	addsd xmm5, xmm3
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

	;# start with rsqM - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb303_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb303_half] ;# rinv 
	movapd  [rsp + nb303_rinvM], xmm4	;# rinvM in xmm4 
	mulsd   xmm7, xmm4
	movapd  [rsp + nb303_rM], xmm7	;# r in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb303_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb303_half] ;# rinv 
	movapd [rsp + nb303_rinvH1], xmm4	;# rinvH1 
	mulsd  xmm6, xmm4
	movapd [rsp + nb303_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb303_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb303_half] ;# rinv 
	movapd [rsp + nb303_rinvH2], xmm4 ;# rinv 
	mulsd xmm5, xmm4
	movapd [rsp + nb303_rH2], xmm5 ;# r 

	;# do M interactions 
	;# rM is still in xmm7 
	mulsd   xmm7, [rsp + nb303_tsc]
	cvttsd2si r8d, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, r8d
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl r8d, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303_VFtab]

	movapd xmm4, [rsi + r8*8]	;# Y1 F1 
	xorpd xmm3, xmm3	
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1  

	movapd xmm6, [rsi + r8*8 + 16]	;# G1 H1 
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb303_two]	;# two*Heps2 
	movapd xmm3, [rsp + nb303_qqM]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul - then we can get rid of mm5 
    addsd  xmm5, [rsp + nb303_vctot]
    movlpd [rsp + nb303_vctot], xmm5 
	xorpd  xmm4, xmm4

	mulsd  xmm3, [rsp + nb303_tsc]
	mulsd  xmm3, [rsp + nb303_rinvM]	
	subsd  xmm4, xmm3

	movapd xmm0, [rsp + nb303_dxM]
	movapd xmm1, [rsp + nb303_dyM]
	movapd xmm2, [rsp + nb303_dzM]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4	;# tx in xmm0-xmm2 

	;# update M forces 
	movapd xmm3, [rsp + nb303_fixM]
	movapd xmm4, [rsp + nb303_fiyM]
	movapd xmm7, [rsp + nb303_fizM]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb303_fixM], xmm3
	movlpd [rsp + nb303_fiyM], xmm4
	movlpd [rsp + nb303_fizM], xmm7
	;# update j forces with water M 
	movlpd [rsp + nb303_fjx], xmm0
	movlpd [rsp + nb303_fjy], xmm1
	movlpd [rsp + nb303_fjz], xmm2

	;# Done with M interactions - now H1! 
	movapd xmm7, [rsp + nb303_rH1]
	mulsd xmm7, [rsp + nb303_tsc]
	cvttsd2si r8d, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, r8d
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl r8d, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303_VFtab]
	
	movapd xmm4, [rsi + r8*8]	;# Y1 F1 
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1  
	unpckhpd xmm5, xmm3	;# F1  

	movapd xmm6, [rsi + r8*8 + 16]	;# G1 H1 
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb303_two]	;# two*Heps2 
	movapd xmm3, [rsp + nb303_qqH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul 
	xorpd  xmm4, xmm4
    addsd  xmm5, [rsp + nb303_vctot]
	mulsd  xmm3, [rsp + nb303_rinvH1]
    movlpd [rsp + nb303_vctot], xmm5 
	mulsd  xmm3, [rsp + nb303_tsc]
	subsd xmm4, xmm3

	movapd xmm0, [rsp + nb303_dxH1]
	movapd xmm1, [rsp + nb303_dyH1]
	movapd xmm2, [rsp + nb303_dzH1]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H1 forces 
	movapd xmm3, [rsp + nb303_fixH1]
	movapd xmm4, [rsp + nb303_fiyH1]
	movapd xmm7, [rsp + nb303_fizH1]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb303_fixH1], xmm3
	movlpd [rsp + nb303_fiyH1], xmm4
	movlpd [rsp + nb303_fizH1], xmm7
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb303_fjx]
	addsd  xmm1, [rsp + nb303_fjy]
	addsd  xmm2, [rsp + nb303_fjz]
	movlpd [rsp + nb303_fjx], xmm0
	movlpd [rsp + nb303_fjy], xmm1
	movlpd [rsp + nb303_fjz], xmm2

	;# Done with H1, finally we do H2 interactions 
	movapd xmm7, [rsp + nb303_rH2]
	mulsd   xmm7, [rsp + nb303_tsc]
	cvttsd2si r8d, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, r8d
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl r8d, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303_VFtab]

	movapd xmm4, [rsi + r8*8]	;# Y1 F1 
	xorpd xmm3, xmm3
	movapd xmm5, xmm4
	unpcklpd xmm4, xmm3	;# Y1 
	unpckhpd xmm5, xmm3	;# F1 

	movapd xmm6, [rsi + r8*8 + 16]	;# G1 H1 
	xorpd xmm3, xmm3
	movapd xmm7, xmm6
	unpcklpd xmm6, xmm3	;# G1 
	unpckhpd xmm7, xmm3	;# H1 
	;# coulomb table ready, in xmm4-xmm7  		
	mulsd  xmm6, xmm1	;# xmm6=Geps 
	mulsd  xmm7, xmm2	;# xmm7=Heps2 
	addsd  xmm5, xmm6
	addsd  xmm5, xmm7	;# xmm5=Fp 	
	mulsd  xmm7, [rsp + nb303_two]	;# two*Heps2 
	movapd xmm3, [rsp + nb303_qqH]
	addsd  xmm7, xmm6
	addsd  xmm7, xmm5 ;# xmm7=FF 
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
	mulsd  xmm3, xmm7 ;# fijC=FF*qq 
    ;# at this point mm5 contains vcoul and xmm3 fijC 
    ;# increment vcoul 
	xorpd  xmm4, xmm4
    addsd  xmm5, [rsp + nb303_vctot]
	mulsd  xmm3, [rsp + nb303_rinvH2]
    movlpd [rsp + nb303_vctot], xmm5 
	mulsd  xmm3, [rsp + nb303_tsc]
	subsd  xmm4, xmm3

	movapd xmm0, [rsp + nb303_dxH2]
	movapd xmm1, [rsp + nb303_dyH2]
	movapd xmm2, [rsp + nb303_dzH2]
	mulsd  xmm0, xmm4
	mulsd  xmm1, xmm4
	mulsd  xmm2, xmm4

	;# update H2 forces 
	movapd xmm3, [rsp + nb303_fixH2]
	movapd xmm4, [rsp + nb303_fiyH2]
	movapd xmm7, [rsp + nb303_fizH2]
	addsd  xmm3, xmm0
	addsd  xmm4, xmm1
	addsd  xmm7, xmm2
	movlpd [rsp + nb303_fixH2], xmm3
	movlpd [rsp + nb303_fiyH2], xmm4
	movlpd [rsp + nb303_fizH2], xmm7

	mov rdi, [rbp + nb303_faction]
	;# update j forces 
	;# update j forces with water H1 
	addsd  xmm0, [rsp + nb303_fjx]
	addsd  xmm1, [rsp + nb303_fjy]
	addsd  xmm2, [rsp + nb303_fjz]

	;# the fj's - start by accumulating forces from memory 
	movlpd xmm3, [rdi + rax*8]
	movlpd xmm4, [rdi + rax*8 + 8]
	movlpd xmm5, [rdi + rax*8 + 16]
	addsd xmm3, xmm0
	addsd xmm4, xmm1
	addsd xmm5, xmm2
	movlpd [rdi + rax*8], xmm3
	movlpd [rdi + rax*8 + 8], xmm4
	movlpd [rdi + rax*8 + 16], xmm5

.nb303_updateouterdata:
	mov   ecx, [rsp + nb303_ii3]
	mov   rdi, [rbp + nb303_faction]
	mov   rsi, [rbp + nb303_fshift]
	mov   edx, [rsp + nb303_is3]

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb303_fixH1]
	movapd xmm1, [rsp + nb303_fiyH1]
	movapd xmm2, [rsp + nb303_fizH1]

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
	movsd  [rdi + rcx*8 + 24],     xmm3
	movsd  [rdi + rcx*8 + 32], xmm4
	movsd  [rdi + rcx*8 + 40], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movapd xmm6, xmm0
	movsd xmm7, xmm2
	unpcklpd xmm6, xmm1

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb303_fixH2]
	movapd xmm1, [rsp + nb303_fiyH2]
	movapd xmm2, [rsp + nb303_fizH2]

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

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movapd xmm0, [rsp + nb303_fixM]
	movapd xmm1, [rsp + nb303_fiyM]
	movapd xmm2, [rsp + nb303_fizM]

	movhlps xmm3, xmm0
	movhlps xmm4, xmm1
	movhlps xmm5, xmm2
	addsd  xmm0, xmm3
	addsd  xmm1, xmm4
	addsd  xmm2, xmm5 ;# sum is in low xmm0-xmm2 

	movapd xmm3, xmm0	
	movapd xmm4, xmm1	
	movapd xmm5, xmm2	

	;# increment i force 
	movsd  xmm3, [rdi + rcx*8 + 72]
	movsd  xmm4, [rdi + rcx*8 + 80]
	movsd  xmm5, [rdi + rcx*8 + 88]
	subsd  xmm3, xmm0
	subsd  xmm4, xmm1
	subsd  xmm5, xmm2
	movsd  [rdi + rcx*8 + 72], xmm3
	movsd  [rdi + rcx*8 + 80], xmm4
	movsd  [rdi + rcx*8 + 88], xmm5

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
	mov esi, [rsp + nb303_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb303_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb303_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb303_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb303_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb303_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb303_n], esi
        jmp .nb303_outer
.nb303_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb303_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb303_end
        ;# non-zero, do one more workunit
        jmp   .nb303_threadloop
.nb303_end:
	mov eax, [rsp + nb303_nouter]
	mov ebx, [rsp + nb303_ninner]
	mov rcx, [rbp + nb303_outeriter]
	mov rdx, [rbp + nb303_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 816
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
	


.globl nb_kernel303nf_x86_64_sse2
.globl _nb_kernel303nf_x86_64_sse2
nb_kernel303nf_x86_64_sse2:	
_nb_kernel303nf_x86_64_sse2:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb303nf_fshift,         16
.equiv          nb303nf_gid,            24
.equiv          nb303nf_pos,            32
.equiv          nb303nf_faction,        40
.equiv          nb303nf_charge,         48
.equiv          nb303nf_p_facel,        56
.equiv          nb303nf_argkrf,         64
.equiv          nb303nf_argcrf,         72
.equiv          nb303nf_Vc,             80
.equiv          nb303nf_type,           88
.equiv          nb303nf_p_ntype,        96
.equiv          nb303nf_vdwparam,       104
.equiv          nb303nf_Vvdw,           112
.equiv          nb303nf_p_tabscale,     120
.equiv          nb303nf_VFtab,          128
.equiv          nb303nf_invsqrta,       136
.equiv          nb303nf_dvda,           144
.equiv          nb303nf_p_gbtabscale,   152
.equiv          nb303nf_GBtab,          160
.equiv          nb303nf_p_nthreads,     168
.equiv          nb303nf_count,          176
.equiv          nb303nf_mtx,            184
.equiv          nb303nf_outeriter,      192
.equiv          nb303nf_inneriter,      200
.equiv          nb303nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb303nf_ixM,            0
.equiv          nb303nf_iyM,            16
.equiv          nb303nf_izM,            32
.equiv          nb303nf_ixH1,           48
.equiv          nb303nf_iyH1,           64
.equiv          nb303nf_izH1,           80
.equiv          nb303nf_ixH2,           96
.equiv          nb303nf_iyH2,           112
.equiv          nb303nf_izH2,           128
.equiv          nb303nf_iqM,            144
.equiv          nb303nf_iqH,            160
.equiv          nb303nf_qqM,            176
.equiv          nb303nf_qqH,            192
.equiv          nb303nf_rinvM,          208
.equiv          nb303nf_rinvH1,         224
.equiv          nb303nf_rinvH2,         240
.equiv          nb303nf_rM,             256
.equiv          nb303nf_rH1,            272
.equiv          nb303nf_rH2,            288
.equiv          nb303nf_tsc,            304
.equiv          nb303nf_vctot,          320
.equiv          nb303nf_half,           336
.equiv          nb303nf_three,          352
.equiv          nb303nf_is3,            368
.equiv          nb303nf_ii3,            372
.equiv          nb303nf_nri,            376
.equiv          nb303nf_iinr,           384
.equiv          nb303nf_jindex,         392
.equiv          nb303nf_jjnr,           400
.equiv          nb303nf_shift,          408
.equiv          nb303nf_shiftvec,       416
.equiv          nb303nf_facel,          424
.equiv          nb303nf_innerjjnr,      432
.equiv          nb303nf_innerk,         440
.equiv          nb303nf_n,              444
.equiv          nb303nf_nn1,            448
.equiv          nb303nf_nouter,         452
.equiv          nb303nf_ninner,         456

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
	sub rsp, 464		;# local variable stack space (n*16+8)

	;# zero 32-bit iteration counters
	mov eax, 0
	mov [rsp + nb303nf_nouter], eax
	mov [rsp + nb303nf_ninner], eax
	
	mov edi, [rdi]
	mov [rsp + nb303nf_nri], edi
	mov [rsp + nb303nf_iinr], rsi
	mov [rsp + nb303nf_jindex], rdx
	mov [rsp + nb303nf_jjnr], rcx
	mov [rsp + nb303nf_shift], r8
	mov [rsp + nb303nf_shiftvec], r9
	mov rsi, [rbp + nb303nf_p_facel]
	movsd xmm0, [rsi]
	movsd [rsp + nb303nf_facel], xmm0

	mov rax, [rbp + nb303nf_p_tabscale]
	movsd xmm3, [rax]
	shufpd xmm3, xmm3, 0
	movapd [rsp + nb303nf_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x00000000     ;# lower half of double half IEEE (hex)
	mov ebx, 0x3fe00000
	mov [rsp + nb303nf_half], eax
	mov [rsp + nb303nf_half+4], ebx
	movsd xmm1, [rsp + nb303nf_half]
	shufpd xmm1, xmm1, 0    ;# splat to all elements
	movapd xmm3, xmm1
	addpd  xmm3, xmm3       ;# one
	movapd xmm2, xmm3
	addpd  xmm2, xmm2       ;# two
	addpd  xmm3, xmm2	;# three
	movapd [rsp + nb303nf_half], xmm1
	movapd [rsp + nb303nf_three], xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb303nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb303nf_charge]
	movsd xmm3, [rdx + rbx*8 + 8]	
	movsd xmm4, [rdx + rbx*8 + 24]	
	mov rsi, [rbp + nb303nf_p_facel]
	movsd xmm0, [rsi]
	movsd xmm5, [rsp + nb303nf_facel]
	mulsd  xmm3, xmm5
	mulsd  xmm4, xmm5

	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	movapd [rsp + nb303nf_iqH], xmm3
	movapd [rsp + nb303nf_iqM], xmm4
	
.nb303nf_threadloop:
        mov   rsi, [rbp + nb303nf_count]        ;# pointer to sync counter
        mov   eax, [rsi]
.nb303nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           	;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb303nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb303nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb303nf_n], eax
        mov [rsp + nb303nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb303nf_outerstart
        jmp .nb303nf_end

.nb303nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb303nf_nouter]
	mov [rsp + nb303nf_nouter], ebx

.nb303nf_outer:
	mov   rax, [rsp + nb303nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 

	mov   rax, [rsp + nb303nf_shiftvec]   ;# rax = base of shiftvec[] 

	movsd xmm0, [rax + rbx*8]
	movsd xmm1, [rax + rbx*8 + 8]
	movsd xmm2, [rax + rbx*8 + 16] 

	mov   rcx, [rsp + nb303nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    ;# ebx =ii 

	movapd xmm3, xmm0
	movapd xmm4, xmm1
	movapd xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb303nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb303nf_ii3], ebx

	addsd xmm3, [rax + rbx*8 + 24]
	addsd xmm4, [rax + rbx*8 + 32]
	addsd xmm5, [rax + rbx*8 + 40]		
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb303nf_ixH1], xmm3
	movapd [rsp + nb303nf_iyH1], xmm4
	movapd [rsp + nb303nf_izH1], xmm5

	movsd xmm3, xmm0
	movsd xmm4, xmm1
	movsd xmm5, xmm2
	addsd xmm0, [rax + rbx*8 + 48]
	addsd xmm1, [rax + rbx*8 + 56]
	addsd xmm2, [rax + rbx*8 + 64]		
	addsd xmm3, [rax + rbx*8 + 72]
	addsd xmm4, [rax + rbx*8 + 80]
	addsd xmm5, [rax + rbx*8 + 88]		

	shufpd xmm0, xmm0, 0
	shufpd xmm1, xmm1, 0
	shufpd xmm2, xmm2, 0
	shufpd xmm3, xmm3, 0
	shufpd xmm4, xmm4, 0
	shufpd xmm5, xmm5, 0
	movapd [rsp + nb303nf_ixH2], xmm0
	movapd [rsp + nb303nf_iyH2], xmm1
	movapd [rsp + nb303nf_izH2], xmm2
	movapd [rsp + nb303nf_ixM], xmm3
	movapd [rsp + nb303nf_iyM], xmm4
	movapd [rsp + nb303nf_izM], xmm5
	
	;# clear vctot 
	xorpd xmm4, xmm4
	movapd [rsp + nb303nf_vctot], xmm4
	
	mov   rax, [rsp + nb303nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb303nf_pos]
	mov   rax, [rsp + nb303nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb303nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  2
	add   ecx, [rsp + nb303nf_ninner]
	mov   [rsp + nb303nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb303nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb303nf_unroll_loop
	jmp   .nb303nf_checksingle
.nb303nf_unroll_loop:
	;# twice unrolled innerloop here 
	mov   rdx, [rsp + nb303nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]

	add qword ptr [rsp + nb303nf_innerjjnr],  8	;# advance pointer (unrolled 2) 
	mov rsi, [rbp + nb303nf_charge]    ;# base of charge[] 
	
	movlpd xmm3, [rsi + rax*8]
	movhpd xmm3, [rsi + rbx*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb303nf_iqM]
	mulpd  xmm4, [rsp + nb303nf_iqH]

	movapd  [rsp + nb303nf_qqM], xmm3
	movapd  [rsp + nb303nf_qqH], xmm4	

	mov rsi, [rbp + nb303nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	

	;# move two coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]
	movhpd xmm0, [rsi + rbx*8]
	movhpd xmm1, [rsi + rbx*8 + 8]
	movhpd xmm2, [rsi + rbx*8 + 16]		

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [rsp + nb303nf_ixM]
	movapd xmm5, [rsp + nb303nf_iyM]
	movapd xmm6, [rsp + nb303nf_izM]

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
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb303nf_ixH1]
	movapd xmm5, [rsp + nb303nf_iyH1]
	movapd xmm6, [rsp + nb303nf_izH1]

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
	movapd xmm3, [rsp + nb303nf_ixH2]
	movapd xmm4, [rsp + nb303nf_iyH2]
	movapd xmm5, [rsp + nb303nf_izH2]

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
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

	;# start with rsqM - put seed in xmm2 
	cvtpd2ps xmm2, xmm7	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303nf_three]
	mulpd   xmm2, xmm7	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb303nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303nf_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb303nf_half] ;# rinv 
	movapd  [rsp + nb303nf_rinvM], xmm4	;# rinvM in xmm4 
	mulpd   xmm7, xmm4
	movapd  [rsp + nb303nf_rM], xmm7	;# r in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtpd2ps xmm2, xmm6	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303nf_three]
	mulpd   xmm2, xmm6	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb303nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303nf_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb303nf_half] ;# rinv 
	movapd [rsp + nb303nf_rinvH1], xmm4	;# rinvH1 
	mulpd  xmm6, xmm4
	movapd [rsp + nb303nf_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtpd2ps xmm2, xmm5	
	rsqrtps xmm2, xmm2
	cvtps2pd xmm2, xmm2

	movapd  xmm3, xmm2
	mulpd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303nf_three]
	mulpd   xmm2, xmm5	;# rsq*lu*lu 
	subpd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulpd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulpd   xmm4, [rsp + nb303nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm4
	mulpd xmm4, xmm4	;# lu*lu 
	mulpd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303nf_three]
	subpd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulpd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulpd xmm4, [rsp + nb303nf_half] ;# rinv 
	movapd [rsp + nb303nf_rinvH2], xmm4 ;# rinv 
	mulpd xmm5, xmm4
	movapd [rsp + nb303nf_rH2], xmm5 ;# r 

	;# do M interactions 
	;# rM is still in xmm7 
	mulpd xmm7, [rsp + nb303nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	movapd xmm3, [rsp + nb303nf_qqM]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    addpd  xmm5, [rsp + nb303nf_vctot]
    movapd [rsp + nb303nf_vctot], xmm5 

	;# Done with M interactions - now H1! 
	movapd xmm7, [rsp + nb303nf_rH1]
	mulpd xmm7, [rsp + nb303nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	movapd xmm3, [rsp + nb303nf_qqH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul 
    addpd  xmm5, [rsp + nb303nf_vctot]
	movapd [rsp + nb303nf_vctot], xmm5
	
	;# Done with H1, finally we do H2 interactions 
	movapd xmm7, [rsp + nb303nf_rH2]
	mulpd   xmm7, [rsp + nb303nf_tsc]
	cvttpd2pi mm6, xmm7	;# mm6 = lu idx 
	cvtpi2pd xmm6, mm6
	subpd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulpd  xmm2, xmm2	;# xmm2=eps2 
	
	pslld mm6, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303nf_VFtab]
	movd eax, mm6
	psrlq mm6, 32
	movd ebx, mm6		;# indices in eax/ebx 

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
	movapd xmm3, [rsp + nb303nf_qqH]
	mulpd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addpd  xmm5, xmm4 ;# xmm5=VV 
	mulpd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul 
    addpd  xmm5, [rsp + nb303nf_vctot]
	movapd [rsp + nb303nf_vctot], xmm5

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb303nf_innerk],  2
	jl    .nb303nf_checksingle
	jmp   .nb303nf_unroll_loop
.nb303nf_checksingle:	
	mov   edx, [rsp + nb303nf_innerk]
	and   edx, 1
	jnz   .nb303nf_dosingle
	jmp   .nb303nf_updateouterdata
.nb303nf_dosingle:
	mov   rdx, [rsp + nb303nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	

	mov rsi, [rbp + nb303nf_charge]    ;# base of charge[] 
	xorpd xmm3, xmm3
	movlpd xmm3, [rsi + rax*8]
	movapd xmm4, xmm3	     
	mulpd  xmm3, [rsp + nb303nf_iqM]
	mulpd  xmm4, [rsp + nb303nf_iqH]

	movapd  [rsp + nb303nf_qqM], xmm3
	movapd  [rsp + nb303nf_qqH], xmm4	

	mov rsi, [rbp + nb303nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	;# move coordinates to xmm0-xmm2 	
	movlpd xmm0, [rsi + rax*8]
	movlpd xmm1, [rsi + rax*8 + 8]
	movlpd xmm2, [rsi + rax*8 + 16]

	;# move ixM-izM to xmm4-xmm6 
	movapd xmm4, [rsp + nb303nf_ixM]
	movapd xmm5, [rsp + nb303nf_iyM]
	movapd xmm6, [rsp + nb303nf_izM]

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
	;# rsqM in xmm7 

	;# move ixH1-izH1 to xmm4-xmm6 
	movapd xmm4, [rsp + nb303nf_ixH1]
	movapd xmm5, [rsp + nb303nf_iyH1]
	movapd xmm6, [rsp + nb303nf_izH1]

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
	movapd xmm3, [rsp + nb303nf_ixH2]
	movapd xmm4, [rsp + nb303nf_iyH2]
	movapd xmm5, [rsp + nb303nf_izH2]

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
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqM in xmm7 

	;# start with rsqM - put seed in xmm2 
	cvtsd2ss xmm2, xmm7	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303nf_three]
	mulsd   xmm2, xmm7	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb303nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm7
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303nf_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb303nf_half] ;# rinv 
	movapd  [rsp + nb303nf_rinvM], xmm4	;# rinvM in xmm4 
	mulsd   xmm7, xmm4
	movapd  [rsp + nb303nf_rM], xmm7	;# r in xmm7 
	
	;# rsqH1 - seed in xmm2 
	cvtsd2ss xmm2, xmm6	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303nf_three]
	mulsd   xmm2, xmm6	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb303nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm6
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303nf_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb303nf_half] ;# rinv 
	movapd [rsp + nb303nf_rinvH1], xmm4	;# rinvH1 
	mulsd  xmm6, xmm4
	movapd [rsp + nb303nf_rH1], xmm6	;# rH1 
	
	;# rsqH2 - seed in xmm2 
	cvtsd2ss xmm2, xmm5	
	rsqrtss xmm2, xmm2
	cvtss2sd xmm2, xmm2

	movapd  xmm3, xmm2
	mulsd   xmm2, xmm2
	movapd  xmm4, [rsp + nb303nf_three]
	mulsd   xmm2, xmm5	;# rsq*lu*lu 
	subsd   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulsd   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulsd   xmm4, [rsp + nb303nf_half] ;# iter1 ( new lu) 

	movapd xmm2, xmm5
	movapd xmm3, xmm4
	mulsd xmm4, xmm4	;# lu*lu 
	mulsd xmm2, xmm4	;# rsq*lu*lu 
	movapd xmm4, [rsp + nb303nf_three]
	subsd xmm4, xmm2	;# 3-rsq*lu*lu 
	mulsd xmm4, xmm3	;# lu*(	3-rsq*lu*lu) 
	mulsd xmm4, [rsp + nb303nf_half] ;# rinv 
	movapd [rsp + nb303nf_rinvH2], xmm4 ;# rinv 
	mulsd xmm5, xmm4
	movapd [rsp + nb303nf_rH2], xmm5 ;# r 

	;# do M interactions 
	movd mm0, eax	
	;# rM is still in xmm7 
	mulsd   xmm7, [rsp + nb303nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303nf_VFtab]

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
	movapd xmm3, [rsp + nb303nf_qqM]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul - then we can get rid of mm5 
    addsd  xmm5, [rsp + nb303nf_vctot]
    movlpd [rsp + nb303nf_vctot], xmm5 

	;# Done with M interactions - now H1! 
	movapd xmm7, [rsp + nb303nf_rH1]
	mulsd xmm7, [rsp + nb303nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303nf_VFtab]
	
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
	movapd xmm3, [rsp + nb303nf_qqH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul 
    addsd  xmm5, [rsp + nb303nf_vctot]
    movlpd [rsp + nb303nf_vctot], xmm5 


	;# Done with H1, finally we do H2 interactions 
	movapd xmm7, [rsp + nb303nf_rH2]
	mulsd   xmm7, [rsp + nb303nf_tsc]
	cvttsd2si eax, xmm7	;# mm6 = lu idx 
	cvtsi2sd xmm6, eax
	subsd xmm7, xmm6
	movapd xmm1, xmm7	;# xmm1=eps 
	movapd xmm2, xmm1	
	mulsd  xmm2, xmm2	;# xmm2=eps2 
	
	shl eax, 2		;# idx *= 4 
	mov  rsi, [rbp + nb303nf_VFtab]

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
	movapd xmm3, [rsp + nb303nf_qqH]
	mulsd  xmm5, xmm1 ;# xmm5=eps*Fp 
	addsd  xmm5, xmm4 ;# xmm5=VV 
	mulsd  xmm5, xmm3 ;# vcoul=qq*VV  
    ;# at this point mm5 contains vcoul 
    ;# increment vcoul 
    addsd  xmm5, [rsp + nb303nf_vctot]
    movlpd [rsp + nb303nf_vctot], xmm5 

.nb303nf_updateouterdata:
	;# get group index for i particle 
	;# get n from stack
	mov esi, [rsp + nb303nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb303nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movapd xmm7, [rsp + nb303nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addsd  xmm7, xmm6	;# low xmm7 has the sum now 
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb303nf_Vc]
	addsd xmm7, [rax + rdx*8] 
	;# move back to mem 
	movsd [rax + rdx*8], xmm7 
	
        ;# finish if last 
        mov ecx, [rsp + nb303nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb303nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb303nf_n], esi
        jmp .nb303nf_outer
.nb303nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb303nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb303nf_end
        ;# non-zero, do one more workunit
        jmp   .nb303nf_threadloop
.nb303nf_end:
	mov eax, [rsp + nb303nf_nouter]
	mov ebx, [rsp + nb303nf_ninner]
	mov rcx, [rbp + nb303nf_outeriter]
	mov rdx, [rbp + nb303nf_inneriter]
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
