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


.globl nb_kernel103_x86_64_sse
.globl _nb_kernel103_x86_64_sse
nb_kernel103_x86_64_sse:	
_nb_kernel103_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb103_fshift,           16
.equiv          nb103_gid,              24
.equiv          nb103_pos,              32
.equiv          nb103_faction,          40
.equiv          nb103_charge,           48
.equiv          nb103_p_facel,          56
.equiv          nb103_argkrf,           64
.equiv          nb103_argcrf,           72
.equiv          nb103_Vc,               80
.equiv          nb103_type,             88
.equiv          nb103_p_ntype,          96
.equiv          nb103_vdwparam,         104
.equiv          nb103_Vvdw,             112
.equiv          nb103_p_tabscale,       120
.equiv          nb103_VFtab,            128
.equiv          nb103_invsqrta,         136
.equiv          nb103_dvda,             144
.equiv          nb103_p_gbtabscale,     152
.equiv          nb103_GBtab,            160
.equiv          nb103_p_nthreads,       168
.equiv          nb103_count,            176
.equiv          nb103_mtx,              184
.equiv          nb103_outeriter,        192
.equiv          nb103_inneriter,        200
.equiv          nb103_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb103_ixH1,             0
.equiv          nb103_iyH1,             16
.equiv          nb103_izH1,             32
.equiv          nb103_ixH2,             48
.equiv          nb103_iyH2,             64
.equiv          nb103_izH2,             80
.equiv          nb103_ixM,              96
.equiv          nb103_iyM,              112
.equiv          nb103_izM,              128
.equiv          nb103_iqH,              144
.equiv          nb103_iqM,              160
.equiv          nb103_dxH1,             176
.equiv          nb103_dyH1,             192
.equiv          nb103_dzH1,             208
.equiv          nb103_dxH2,             224
.equiv          nb103_dyH2,             240
.equiv          nb103_dzH2,             256
.equiv          nb103_dxM,              272
.equiv          nb103_dyM,              288
.equiv          nb103_dzM,              304
.equiv          nb103_qqH,              320
.equiv          nb103_qqM,              336
.equiv          nb103_vctot,            352
.equiv          nb103_fixH1,            368
.equiv          nb103_fiyH1,            384
.equiv          nb103_fizH1,            400
.equiv          nb103_fixH2,            416
.equiv          nb103_fiyH2,            432
.equiv          nb103_fizH2,            448
.equiv          nb103_fixM,             464
.equiv          nb103_fiyM,             480
.equiv          nb103_fizM,             496
.equiv          nb103_fjx,              512
.equiv          nb103_fjy,              528
.equiv          nb103_fjz,              544
.equiv          nb103_half,             560
.equiv          nb103_three,            576
.equiv          nb103_is3,              592
.equiv          nb103_ii3,              596
.equiv          nb103_nri,              600
.equiv          nb103_innerjjnr,        608
.equiv          nb103_iinr,             616
.equiv          nb103_jindex,           624
.equiv          nb103_jjnr,             632
.equiv          nb103_shift,            640
.equiv          nb103_shiftvec,         648
.equiv          nb103_facel,            656
.equiv          nb103_innerk,           664
.equiv          nb103_n,                668
.equiv          nb103_nn1,              672
.equiv          nb103_nouter,           676
.equiv          nb103_ninner,           680

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
	sub rsp, 688            ; # local variable stack space (n*16+8)
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
	mov [rsp + nb103_nouter], eax
	mov [rsp + nb103_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb103_nri], edi
	mov [rsp + nb103_iinr], rsi
	mov [rsp + nb103_jindex], rdx
	mov [rsp + nb103_jjnr], rcx
	mov [rsp + nb103_shift], r8
	mov [rsp + nb103_shiftvec], r9
	mov rsi, [rbp + nb103_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb103_facel], xmm0

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb103_half], eax
	movss xmm1, [rsp + nb103_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb103_half],  xmm1
	movaps [rsp + nb103_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb103_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb103_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb103_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb103_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb103_iqH], xmm3
	movaps [rsp + nb103_iqM], xmm4

.nb103_threadloop:
        mov   rsi, [rbp + nb103_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb103_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb103_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb103_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb103_n], eax
        mov [rsp + nb103_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi

        jg  .nb103_outerstart
        jmp .nb103_end
	
.nb103_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb103_nouter]
	mov [rsp + nb103_nouter], ebx

.nb103_outer:
	mov   rax, [rsp + nb103_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb103_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb103_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb103_iinr]   	;# rcx = pointer into iinr[] 
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb103_pos]	;# rax = base of pos[]  
	mov   [rsp + nb103_ii3], ebx

	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]	
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb103_ixH1], xmm3
	movaps [rsp + nb103_iyH1], xmm4
	movaps [rsp + nb103_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 24]
	addss xmm1, [rax + rbx*4 + 28]
	addss xmm2, [rax + rbx*4 + 32]		
	addss xmm3, [rax + rbx*4 + 36]
	addss xmm4, [rax + rbx*4 + 40]
	addss xmm5, [rax + rbx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb103_ixH2], xmm0
	movaps [rsp + nb103_iyH2], xmm1
	movaps [rsp + nb103_izH2], xmm2
	movaps [rsp + nb103_ixM], xmm3
	movaps [rsp + nb103_iyM], xmm4
	movaps [rsp + nb103_izM], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb103_vctot], xmm4
	movaps [rsp + nb103_fixH1], xmm4
	movaps [rsp + nb103_fiyH1], xmm4
	movaps [rsp + nb103_fizH1], xmm4
	movaps [rsp + nb103_fixH2], xmm4
	movaps [rsp + nb103_fiyH2], xmm4
	movaps [rsp + nb103_fizH2], xmm4
	movaps [rsp + nb103_fixM], xmm4
	movaps [rsp + nb103_fiyM], xmm4
	movaps [rsp + nb103_fizM], xmm4
	
	mov   rax, [rsp + nb103_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb103_pos]
	mov   rdi, [rbp + nb103_faction]	
	mov   rax, [rsp + nb103_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb103_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb103_ninner]
	mov   [rsp + nb103_ninner], ecx
	add   edx, 0
	mov   [rsp + nb103_innerk], edx	;# number of innerloop atoms 

	jge   .nb103_unroll_loop
	jmp   .nb103_odd_inner
.nb103_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb103_innerjjnr] 	;# pointer to jjnr[k] 
	mov   r8d, [rdx]	
	mov   r9d, [rdx + 4]              
	mov   r10d, [rdx + 8]            
	mov   r11d, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb103_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb103_charge]	;# base of charge[] 
	
	movss xmm12, [rsi + r8*4]
	movss xmm13, [rsi + r10*4]
	movss xmm14, [rsi + r9*4]
	movss xmm15, [rsi + r11*4]

	shufps xmm12, xmm14, 0
	shufps xmm13, xmm15, 0
	shufps xmm12, xmm13, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm13, xmm12	     ;# and in xmm4 
	mulps  xmm12, [rsp + nb103_iqH]
	mulps  xmm13, [rsp + nb103_iqM]

	movaps  [rsp + nb103_qqH], xmm12
	movaps  [rsp + nb103_qqM], xmm13

	mov rsi, [rbp + nb103_pos]   	;# base of pos[] 

	lea   rax, [r8 + r8*2]     ;# j3 
	lea   rbx, [r9 + r9*2]	
	lea   rcx, [r10 + r10*2]     
	lea   rdx, [r11 + r11*2]	

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

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    
    subps xmm0, [rsp + nb103_ixH1]
    subps xmm1, [rsp + nb103_iyH1]
    subps xmm2, [rsp + nb103_izH1]
    subps xmm3, [rsp + nb103_ixH2]
    subps xmm4, [rsp + nb103_iyH2]
    subps xmm5, [rsp + nb103_izH2]
    subps xmm6, [rsp + nb103_ixM]
    subps xmm7, [rsp + nb103_iyM]
    subps xmm8, [rsp + nb103_izM]
    
	movaps [rsp + nb103_dxH1], xmm0
	movaps [rsp + nb103_dyH1], xmm1
	movaps [rsp + nb103_dzH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb103_dxH2], xmm3
	movaps [rsp + nb103_dyH2], xmm4
	movaps [rsp + nb103_dzH2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb103_dxM], xmm6
	movaps [rsp + nb103_dyM], xmm7
	movaps [rsp + nb103_dzM], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for j atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb103_three]
	movaps  xmm10, xmm9
    movaps  xmm11, xmm9

	mulps   xmm1, xmm0 ;# rsq*lu*lu
	mulps   xmm4, xmm3 ;# rsq*lu*lu 
    mulps   xmm7, xmm6 ;# rsq*lu*lu
	
	subps   xmm9, xmm1
	subps   xmm10, xmm4
    subps   xmm11, xmm7 ;# 3-rsq*lu*lu

	mulps   xmm9, xmm2
	mulps   xmm10, xmm5
    mulps   xmm11, xmm8 ;# lu*(3-rsq*lu*lu)

	movaps  xmm0, [rsp + nb103_half]
	mulps   xmm9, xmm0  ;# rinvH1
	mulps   xmm10, xmm0 ;# rinvH2
    mulps   xmm11, xmm0 ;# rinvM
	
	;# interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb103_qqH] 
    mulps  xmm1, [rsp + nb103_qqH] 
    mulps  xmm2, [rsp + nb103_qqM] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb103_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb103_vctot], xmm0
    
	;# move j forces to local temp variables 
    movlps xmm0, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm1, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm0, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm1, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm2, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm3, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm4, [rdi + rbx*4 + 8] ;# jzb
    movss  xmm5, [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm2, xmm4 ;# jza  -  jzb  -
    movlhps xmm3, xmm5 ;# jzc  -  jzd -

    shufps xmm2, xmm3,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm0: jxa jya jxb jyb 
    ;# xmm1: jxc jyc jxd jyd
    ;# xmm2: jza jzb jzc jzd

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb103_dxH1]
	mulps xmm8, [rsp + nb103_dyH1]
	mulps xmm9, [rsp + nb103_dzH1]
	mulps xmm10, [rsp + nb103_dxH2]
	mulps xmm11, [rsp + nb103_dyH2]
	mulps xmm12, [rsp + nb103_dzH2]
	mulps xmm13, [rsp + nb103_dxM]
	mulps xmm14, [rsp + nb103_dyM]
	mulps xmm15, [rsp + nb103_dzM]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb103_fixH1]
    addps xmm8, [rsp + nb103_fiyH1]
    addps xmm9, [rsp + nb103_fizH1]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb103_fixH2]
    addps xmm11, [rsp + nb103_fiyH2]
    addps xmm12, [rsp + nb103_fizH2]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb103_fixM]
    addps xmm14, [rsp + nb103_fiyM]
    addps xmm15, [rsp + nb103_fizM]

    movaps [rsp + nb103_fixH1], xmm7
    movaps [rsp + nb103_fiyH1], xmm8
    movaps [rsp + nb103_fizH1], xmm9
    movaps [rsp + nb103_fixH2], xmm10
    movaps [rsp + nb103_fiyH2], xmm11
    movaps [rsp + nb103_fizH2], xmm12
    movaps [rsp + nb103_fixM], xmm13
    movaps [rsp + nb103_fiyM], xmm14
    movaps [rsp + nb103_fizM], xmm15
    
    ;# xmm3 = fjx , xmm4 = fjy
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fzc fzd
    
    movlps [rdi + rax*4], xmm0
    movhps [rdi + rbx*4], xmm0
    movlps [rdi + rcx*4], xmm1
    movhps [rdi + rdx*4], xmm1
    movss  [rdi + rax*4 + 8], xmm2
    movss  [rdi + rcx*4 + 8], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 8], xmm2
    movss  [rdi + rdx*4 + 8], xmm3

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb103_innerk],  4
	jl    .nb103_odd_inner
	jmp   .nb103_unroll_loop
.nb103_odd_inner:	
	add dword ptr [rsp + nb103_innerk],  4
	jnz   .nb103_odd_loop
	jmp   .nb103_updateouterdata
.nb103_odd_loop:
	mov   rdx, [rsp + nb103_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb103_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb103_iqM]
	mov rsi, [rbp + nb103_charge] 
	movhps xmm4, [rsp + nb103_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb103_qqM], xmm3	;# use dummy qq for storage 

	mov rsi, [rbp + nb103_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	movss xmm0, [rsp + nb103_ixM]
	movss xmm1, [rsp + nb103_iyM]
	movss xmm2, [rsp + nb103_izM]
	
	movlps xmm6, [rsp + nb103_ixH1]
	movlps xmm7, [rsp + nb103_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm0, xmm6
	movlps xmm6, [rsp + nb103_iyH1]
	movlps xmm7, [rsp + nb103_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm1, xmm6
	movlps xmm6, [rsp + nb103_izH1]
	movlps xmm7, [rsp + nb103_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm2, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# use dummy dx for storage
	movaps [rsp + nb103_dxM], xmm3
	movaps [rsp + nb103_dyM], xmm4
	movaps [rsp + nb103_dzM], xmm5

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb103_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb103_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000	
	movaps xmm4, xmm0
	mulps  xmm4, xmm4	;# xmm4=rinvsq 
	movaps xmm3, [rsp + nb103_qqM]

	mulps  xmm3, xmm0	;# xmm3=vcoul 
	mulps  xmm4, xmm3	;# xmm4=total fscal 
	addps  xmm3, [rsp + nb103_vctot]

	movaps xmm0, [rsp + nb103_dxM]
	movaps xmm1, [rsp + nb103_dyM]
	movaps xmm2, [rsp + nb103_dzM]

	movaps [rsp + nb103_vctot], xmm3

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4
	;# xmm0-xmm2 contains tx-tz (partial force) 
	movss  xmm3, [rsp + nb103_fixM]	
	movss  xmm4, [rsp + nb103_fiyM]	
	movss  xmm5, [rsp + nb103_fizM]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb103_fixM], xmm3	
	movss  [rsp + nb103_fiyM], xmm4	
	movss  [rsp + nb103_fizM], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# 11100110
	shufps xmm5, xmm5, 230 ;# 11100110
	addss  xmm3, [rsp + nb103_fixH1]
	addss  xmm4, [rsp + nb103_fiyH1]
	addss  xmm5, [rsp + nb103_fizH1]
	movss  [rsp + nb103_fixH1], xmm3	
	movss  [rsp + nb103_fiyH1], xmm4	
	movss  [rsp + nb103_fizH1], xmm5	;# updated the H1 force 

	mov rdi, [rbp + nb103_faction]
	shufps xmm3, xmm3, 231 ;# 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# 11100111
	shufps xmm5, xmm5, 231 ;# 11100111
	addss  xmm3, [rsp + nb103_fixH2]
	addss  xmm4, [rsp + nb103_fiyH2]
	addss  xmm5, [rsp + nb103_fizH2]
	movss  [rsp + nb103_fixH2], xmm3	
	movss  [rsp + nb103_fiyH2], xmm4	
	movss  [rsp + nb103_fizH2], xmm5	;# updated the H2 force 

	;# the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
	xorps  xmm5, xmm5
	movaps xmm3, xmm0
	movlps xmm6, [rdi + rax*4]
	movss  xmm7, [rdi + rax*4 + 8]
	unpcklps xmm3, xmm1
	movlhps  xmm3, xmm5	
	unpckhps xmm0, xmm1		
	addps    xmm0, xmm3
	movhlps  xmm3, xmm0	
	addps    xmm0, xmm3	;# x,y sum in xmm0 

	movhlps  xmm1, xmm2
	addss    xmm2, xmm1
	shufps   xmm1, xmm1, 1 
	addss    xmm2, xmm1	;# z sum in xmm2 
	addps    xmm6, xmm0
	addss    xmm7, xmm2
	
	movlps [rdi + rax*4], xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec   dword ptr [rsp + nb103_innerk]
	jz    .nb103_updateouterdata
	jmp   .nb103_odd_loop
.nb103_updateouterdata:

	mov   ecx, [rsp + nb103_ii3]
	mov   rdi, [rbp + nb103_faction]
	mov   rsi, [rbp + nb103_fshift]
	mov   edx, [rsp + nb103_is3]

	;# accumulate H1 forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb103_fixH1]
	movaps xmm1, [rsp + nb103_fiyH1]
	movaps xmm2, [rsp + nb103_fizH1]

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
	movss  xmm3, [rdi + rcx*4 + 12]
	movss  xmm4, [rdi + rcx*4 + 16]
	movss  xmm5, [rdi + rcx*4 + 20]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 12], xmm3
	movss  [rdi + rcx*4 + 16], xmm4
	movss  [rdi + rcx*4 + 20], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H2 i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb103_fixH2]
	movaps xmm1, [rsp + nb103_fiyH2]
	movaps xmm2, [rsp + nb103_fizH2]

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
	movss  xmm3, [rdi + rcx*4 + 24]
	movss  xmm4, [rdi + rcx*4 + 28]
	movss  xmm5, [rdi + rcx*4 + 32]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 24], xmm3
	movss  [rdi + rcx*4 + 28], xmm4
	movss  [rdi + rcx*4 + 32], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate M i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb103_fixM]
	movaps xmm1, [rsp + nb103_fiyM]
	movaps xmm2, [rsp + nb103_fizM]

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
	movss  xmm3, [rdi + rcx*4 + 36]
	movss  xmm4, [rdi + rcx*4 + 40]
	movss  xmm5, [rdi + rcx*4 + 44]
	subss  xmm3, xmm0
	subss  xmm4, xmm1
	subss  xmm5, xmm2
	movss  [rdi + rcx*4 + 36], xmm3
	movss  [rdi + rcx*4 + 40], xmm4
	movss  [rdi + rcx*4 + 44], xmm5

	;# accumulate force in xmm6/xmm7 for fshift 
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]

	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb103_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb103_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb103_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb103_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 	

        ;# finish if last 
        mov ecx, [rsp + nb103_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi

        jz .nb103_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb103_n], esi
        jmp .nb103_outer
.nb103_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb103_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb103_end
        ;# non-zero, do one more workunit
        jmp   .nb103_threadloop
.nb103_end:

	mov eax, [rsp + nb103_nouter]
	mov ebx, [rsp + nb103_ninner]
	mov rcx, [rbp + nb103_outeriter]
	mov rdx, [rbp + nb103_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 688
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



.globl nb_kernel103nf_x86_64_sse
.globl _nb_kernel103nf_x86_64_sse
nb_kernel103nf_x86_64_sse:	
_nb_kernel103nf_x86_64_sse:	
.equiv          nb103nf_fshift,         16
.equiv          nb103nf_gid,            24
.equiv          nb103nf_pos,            32
.equiv          nb103nf_faction,        40
.equiv          nb103nf_charge,         48
.equiv          nb103nf_p_facel,        56
.equiv          nb103nf_argkrf,         64
.equiv          nb103nf_argcrf,         72
.equiv          nb103nf_Vc,             80
.equiv          nb103nf_type,           88
.equiv          nb103nf_p_ntype,        96
.equiv          nb103nf_vdwparam,       104
.equiv          nb103nf_Vvdw,           112
.equiv          nb103nf_p_tabscale,     120
.equiv          nb103nf_VFtab,          128
.equiv          nb103nf_invsqrta,       136
.equiv          nb103nf_dvda,           144
.equiv          nb103nf_p_gbtabscale,   152
.equiv          nb103nf_GBtab,          160
.equiv          nb103nf_p_nthreads,     168
.equiv          nb103nf_count,          176
.equiv          nb103nf_mtx,            184
.equiv          nb103nf_outeriter,      192
.equiv          nb103nf_inneriter,      200
.equiv          nb103nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb103nf_ixH1,           0
.equiv          nb103nf_iyH1,           16
.equiv          nb103nf_izH1,           32
.equiv          nb103nf_ixH2,           48
.equiv          nb103nf_iyH2,           64
.equiv          nb103nf_izH2,           80
.equiv          nb103nf_ixM,            96
.equiv          nb103nf_iyM,            112
.equiv          nb103nf_izM,            128
.equiv          nb103nf_iqH,            144
.equiv          nb103nf_iqM,            160
.equiv          nb103nf_vctot,          176
.equiv          nb103nf_half,           192
.equiv          nb103nf_three,          208
.equiv          nb103nf_qqH,            224
.equiv          nb103nf_qqM,            240
.equiv          nb103nf_is3,            256
.equiv          nb103nf_ii3,            260
.equiv          nb103nf_nri,            264
.equiv          nb103nf_iinr,           272
.equiv          nb103nf_jindex,         280
.equiv          nb103nf_jjnr,           288
.equiv          nb103nf_shift,          296
.equiv          nb103nf_shiftvec,       304
.equiv          nb103nf_facel,          312
.equiv          nb103nf_innerjjnr,      320
.equiv          nb103nf_innerk,         328
.equiv          nb103nf_n,              332
.equiv          nb103nf_nn1,            336
.equiv          nb103nf_nouter,         340
.equiv          nb103nf_ninner,         344

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
	sub rsp, 352
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
	mov [rsp + nb103nf_nouter], eax
	mov [rsp + nb103nf_ninner], eax


	mov edi, [rdi]
	mov [rsp + nb103nf_nri], edi
	mov [rsp + nb103nf_iinr], rsi
	mov [rsp + nb103nf_jindex], rdx
	mov [rsp + nb103nf_jjnr], rcx
	mov [rsp + nb103nf_shift], r8
	mov [rsp + nb103nf_shiftvec], r9
	mov rsi, [rbp + nb103nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb103nf_facel], xmm0


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb103nf_half], eax
	movss xmm1, [rsp + nb103nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb103nf_half],  xmm1
	movaps [rsp + nb103nf_three],  xmm3

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb103nf_iinr]   	;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]		;# ebx =ii 

	mov   rdx, [rbp + nb103nf_charge]
	movss xmm3, [rdx + rbx*4 + 4]	
	movss xmm4, [rdx + rbx*4 + 12]	
	mov rsi, [rbp + nb103nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb103nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb103nf_iqH], xmm3
	movaps [rsp + nb103nf_iqM], xmm4
	
.nb103nf_threadloop:
        mov   rsi, [rbp + nb103nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb103nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb103nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb103nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb103nf_n], eax
        mov [rsp + nb103nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb103nf_outerstart
        jmp .nb103nf_end
	
.nb103nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb103nf_nouter]
	mov [rsp + nb103nf_nouter], ebx

.nb103nf_outer:
	mov   rax, [rsp + nb103nf_shift]  	;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]	;# rbx=3*is 
	mov   [rsp + nb103nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb103nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb103nf_iinr]   	;# rcx = pointer into iinr[] 
	mov   ebx, [rcx + rsi*4]		;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb103nf_pos]	;# rax = base of pos[]  
	mov   [rsp + nb103nf_ii3], ebx

	addss xmm3, [rax + rbx*4 + 12]
	addss xmm4, [rax + rbx*4 + 16]
	addss xmm5, [rax + rbx*4 + 20]	
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb103nf_ixH1], xmm3
	movaps [rsp + nb103nf_iyH1], xmm4
	movaps [rsp + nb103nf_izH1], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 24]
	addss xmm1, [rax + rbx*4 + 28]
	addss xmm2, [rax + rbx*4 + 32]		
	addss xmm3, [rax + rbx*4 + 36]
	addss xmm4, [rax + rbx*4 + 40]
	addss xmm5, [rax + rbx*4 + 44]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb103nf_ixH2], xmm0
	movaps [rsp + nb103nf_iyH2], xmm1
	movaps [rsp + nb103nf_izH2], xmm2
	movaps [rsp + nb103nf_ixM], xmm3
	movaps [rsp + nb103nf_iyM], xmm4
	movaps [rsp + nb103nf_izM], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb103nf_vctot], xmm4

	mov   rax, [rsp + nb103nf_jindex]
	mov   ecx, [rax + rsi*4]	 	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	 	;# jindex[n+1] 
	sub   edx, ecx           	;# number of innerloop atoms 

	mov   rsi, [rbp + nb103nf_pos]
	mov   rdi, [rbp + nb103nf_faction]	
	mov   rax, [rsp + nb103nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb103nf_innerjjnr], rax 	;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb103nf_ninner]
	mov   [rsp + nb103nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb103nf_innerk], edx	;# number of innerloop atoms 
	jge   .nb103nf_unroll_loop
	jmp   .nb103nf_odd_inner
.nb103nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb103nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]     	;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb103nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb103nf_charge]	;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0
	shufps xmm4, xmm7, 0
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	 	;# and in xmm4 
	mulps  xmm3, [rsp + nb103nf_iqH]
	mulps  xmm4, [rsp + nb103nf_iqM]

	movaps  [rsp + nb103nf_qqH], xmm3
	movaps  [rsp + nb103nf_qqM], xmm4	

	mov rsi, [rbp + nb103nf_pos]   	;# base of pos[] 

	lea   rax, [rax + rax*2] 	;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2] 	;# replace jnr with j3 
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

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [rsp + nb103nf_ixH1]
	movaps xmm5, [rsp + nb103nf_iyH1]
	movaps xmm6, [rsp + nb103nf_izH1]

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
	movaps xmm7, xmm4
	;# rsqO in xmm7 

	;# move ixH2-izH2 to xmm4-xmm6 
	movaps xmm4, [rsp + nb103nf_ixH2]
	movaps xmm5, [rsp + nb103nf_iyH2]
	movaps xmm6, [rsp + nb103nf_izH2]

	;# calc dr 
	subps xmm4, xmm0
	subps xmm5, xmm1
	subps xmm6, xmm2

	;# square it 
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	mulps xmm6,xmm6
	addps xmm6, xmm5
	addps xmm6, xmm4
	;# rsqH1 in xmm6 

	;# move ixM-izM to xmm3-xmm5  
	movaps xmm3, [rsp + nb103nf_ixM]
	movaps xmm4, [rsp + nb103nf_iyM]
	movaps xmm5, [rsp + nb103nf_izM]

	;# calc dr 
	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	;# square it 
	mulps xmm3,xmm3
	mulps xmm4,xmm4
	mulps xmm5,xmm5
	addps xmm5, xmm4
	addps xmm5, xmm3
	;# rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

	;# start with rsqH1 - seed in xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb103nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb103nf_half]
	movaps  xmm7, xmm4	;# rinvH1 in xmm7 
	;# rsqH2 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb103nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb103nf_half]
	movaps  xmm6, xmm4	;# rinvH2 in xmm6 
	;# rsqM - seed in xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb103nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb103nf_half]
	movaps  xmm5, xmm4	;# rinvM in xmm5 

	;# do H1 interactions - xmm7=rinv
	mulps  xmm7, [rsp + nb103nf_qqH]	;# xmm7=vcoul 
	addps  xmm7, [rsp + nb103nf_vctot]
	movaps [rsp + nb103nf_vctot], xmm7

	;# H2 interactions - xmm6=rinv
	mulps  xmm6, [rsp + nb103nf_qqH]	;# xmm6=vcoul 
	addps  xmm6, xmm7
	movaps [rsp + nb103nf_vctot], xmm6

	;# M interactions  - xmm5=rinv
	mulps  xmm5, [rsp + nb103nf_qqM]	;# xmm5=vcoul 
	addps  xmm5, xmm6
	movaps [rsp + nb103nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb103nf_innerk],  4
	jl    .nb103nf_odd_inner
	jmp   .nb103nf_unroll_loop
.nb103nf_odd_inner:	
	add dword ptr [rsp + nb103nf_innerk],  4
	jnz   .nb103nf_odd_loop
	jmp   .nb103nf_updateouterdata
.nb103nf_odd_loop:
	mov   rdx, [rsp + nb103nf_innerjjnr] 	;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb103nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb103nf_iqM]
	mov rsi, [rbp + nb103nf_charge] 
	movhps xmm4, [rsp + nb103nf_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb103nf_qqM], xmm3	;# use dummy qq for storage 

	mov rsi, [rbp + nb103nf_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb103nf_ixM]
	movss xmm4, [rsp + nb103nf_iyM]
	movss xmm5, [rsp + nb103nf_izM]
	
	movlps xmm6, [rsp + nb103nf_ixH1]
	movlps xmm7, [rsp + nb103nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb103nf_iyH1]
	movlps xmm7, [rsp + nb103nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb103nf_izH1]
	movlps xmm7, [rsp + nb103nf_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm5, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5

	addps  xmm4, xmm3
	addps  xmm4, xmm5
	;# rsq in xmm4 

	rsqrtps xmm5, xmm4
	;# lookup seed in xmm5 
	movaps xmm2, xmm5
	mulps xmm5, xmm5
	movaps xmm1, [rsp + nb103nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb103nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000	- xmm0=rinv
	movaps xmm3, [rsp + nb103nf_qqM]
	mulps  xmm3, xmm0	;# xmm3=vcoul 
	addps  xmm3, [rsp + nb103nf_vctot]
	movaps [rsp + nb103nf_vctot], xmm3

	dec   dword ptr [rsp + nb103nf_innerk]
	jz    .nb103nf_updateouterdata
	jmp   .nb103nf_odd_loop
.nb103nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb103nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb103nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb103nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb103nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 	
	
        ;# finish if last 
        mov ecx, [rsp + nb103nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb103nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb103nf_n], esi
        jmp .nb103nf_outer
.nb103nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb103nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb103nf_end
        ;# non-zero, do one more workunit
        jmp   .nb103nf_threadloop
.nb103nf_end:

	mov eax, [rsp + nb103nf_nouter]
	mov ebx, [rsp + nb103nf_ninner]
	mov rcx, [rbp + nb103nf_outeriter]
	mov rdx, [rbp + nb103nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 352
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
