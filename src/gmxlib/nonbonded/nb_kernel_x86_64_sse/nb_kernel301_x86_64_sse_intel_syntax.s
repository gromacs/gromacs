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



.globl nb_kernel301_x86_64_sse
.globl _nb_kernel301_x86_64_sse
nb_kernel301_x86_64_sse:	
_nb_kernel301_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb301_fshift,           16
.equiv          nb301_gid,              24
.equiv          nb301_pos,              32
.equiv          nb301_faction,          40
.equiv          nb301_charge,           48
.equiv          nb301_p_facel,          56
.equiv          nb301_argkrf,           64
.equiv          nb301_argcrf,           72
.equiv          nb301_Vc,               80
.equiv          nb301_type,             88
.equiv          nb301_p_ntype,          96
.equiv          nb301_vdwparam,         104
.equiv          nb301_Vvdw,             112
.equiv          nb301_p_tabscale,       120
.equiv          nb301_VFtab,            128
.equiv          nb301_invsqrta,         136
.equiv          nb301_dvda,             144
.equiv          nb301_p_gbtabscale,     152
.equiv          nb301_GBtab,            160
.equiv          nb301_p_nthreads,       168
.equiv          nb301_count,            176
.equiv          nb301_mtx,              184
.equiv          nb301_outeriter,        192
.equiv          nb301_inneriter,        200
.equiv          nb301_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb301_ixO,              0
.equiv          nb301_iyO,              16
.equiv          nb301_izO,              32
.equiv          nb301_ixH1,             48
.equiv          nb301_iyH1,             64
.equiv          nb301_izH1,             80
.equiv          nb301_ixH2,             96
.equiv          nb301_iyH2,             112
.equiv          nb301_izH2,             128
.equiv          nb301_iqO,              144
.equiv          nb301_iqH,              160
.equiv          nb301_dxO,              176
.equiv          nb301_dyO,              192
.equiv          nb301_dzO,              208
.equiv          nb301_dxH1,             224
.equiv          nb301_dyH1,             240
.equiv          nb301_dzH1,             256
.equiv          nb301_dxH2,             272
.equiv          nb301_dyH2,             288
.equiv          nb301_dzH2,             304
.equiv          nb301_qqO,              320
.equiv          nb301_qqH,              336
.equiv          nb301_rinvO,            352
.equiv          nb301_rinvH1,           368
.equiv          nb301_rinvH2,           384
.equiv          nb301_rO,               400
.equiv          nb301_rH1,              416
.equiv          nb301_rH2,              432
.equiv          nb301_tsc,              448
.equiv          nb301_two,              464
.equiv          nb301_vctot,            480
.equiv          nb301_fixO,             496
.equiv          nb301_fiyO,             512
.equiv          nb301_fizO,             528
.equiv          nb301_fixH1,            544
.equiv          nb301_fiyH1,            560
.equiv          nb301_fizH1,            576
.equiv          nb301_fixH2,            592
.equiv          nb301_fiyH2,            608
.equiv          nb301_fizH2,            624
.equiv          nb301_epsO,             640
.equiv          nb301_epsH1,            656
.equiv          nb301_epsH2,            672
.equiv          nb301_half,             688
.equiv          nb301_three,            704
.equiv          nb301_is3,              720
.equiv          nb301_ii3,              724
.equiv          nb301_nri,              728
.equiv          nb301_iinr,             736
.equiv          nb301_jindex,           744
.equiv          nb301_jjnr,             752
.equiv          nb301_shift,            760
.equiv          nb301_shiftvec,         768
.equiv          nb301_facel,            776
.equiv          nb301_innerjjnr,        784
.equiv          nb301_innerk,           792
.equiv          nb301_n,                796
.equiv          nb301_nn1,              800
.equiv          nb301_nouter,           804
.equiv          nb301_ninner,           808

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
	sub rsp, 816		;# local variable stack space (n*16+8)
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
	mov [rsp + nb301_nouter], eax
	mov [rsp + nb301_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb301_nri], edi
	mov [rsp + nb301_iinr], rsi
	mov [rsp + nb301_jindex], rdx
	mov [rsp + nb301_jjnr], rcx
	mov [rsp + nb301_shift], r8
	mov [rsp + nb301_shiftvec], r9
	mov rsi, [rbp + nb301_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb301_facel], xmm0

	mov rax, [rbp + nb301_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb301_tsc], xmm3

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb301_half], eax
	movss xmm1, [rsp + nb301_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb301_half],  xmm1
	movaps [rsp + nb301_two],  xmm2
	movaps [rsp + nb301_three],  xmm3
	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb301_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb301_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb301_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb301_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb301_iqO], xmm3
	movaps [rsp + nb301_iqH], xmm4
	
.nb301_threadloop:
        mov   rsi, [rbp + nb301_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb301_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb301_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb301_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb301_n], eax
        mov [rsp + nb301_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb301_outerstart
        jmp .nb301_end
	
.nb301_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb301_nouter]
	mov [rsp + nb301_nouter], ebx

.nb301_outer:
	mov   rax, [rsp + nb301_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb301_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb301_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb301_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb301_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb301_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb301_ixO], xmm3
	movaps [rsp + nb301_iyO], xmm4
	movaps [rsp + nb301_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb301_ixH1], xmm0
	movaps [rsp + nb301_iyH1], xmm1
	movaps [rsp + nb301_izH1], xmm2
	movaps [rsp + nb301_ixH2], xmm3
	movaps [rsp + nb301_iyH2], xmm4
	movaps [rsp + nb301_izH2], xmm5
	
	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb301_vctot], xmm4
	movaps [rsp + nb301_fixO], xmm4
	movaps [rsp + nb301_fiyO], xmm4
	movaps [rsp + nb301_fizO], xmm4
	movaps [rsp + nb301_fixH1], xmm4
	movaps [rsp + nb301_fiyH1], xmm4
	movaps [rsp + nb301_fizH1], xmm4
	movaps [rsp + nb301_fixH2], xmm4
	movaps [rsp + nb301_fiyH2], xmm4
	movaps [rsp + nb301_fizH2], xmm4
	
	mov   rax, [rsp + nb301_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb301_pos]
	mov   rdi, [rbp + nb301_faction]	
	mov   rax, [rsp + nb301_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb301_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb301_ninner]
	mov   [rsp + nb301_ninner], ecx
	add   edx, 0
	mov   [rsp + nb301_innerk], edx    ;# number of innerloop atoms 
	jge   .nb301_unroll_loop
	jmp   .nb301_odd_inner
.nb301_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb301_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb301_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb301_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb301_iqO]
	mulps  xmm4, [rsp + nb301_iqH]

	movaps  [rsp + nb301_qqO], xmm3
	movaps  [rsp + nb301_qqH], xmm4	

	mov rsi, [rbp + nb301_pos]       ;# base of pos[] 

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

    ;# xmm0 = jx
    ;# xmm1 = jy
    ;# xmm2 = jz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb301_ixO]
    subps xmm1, [rsp + nb301_iyO]
    subps xmm2, [rsp + nb301_izO]
    subps xmm3, [rsp + nb301_ixH1]
    subps xmm4, [rsp + nb301_iyH1]
    subps xmm5, [rsp + nb301_izH1]
    subps xmm6, [rsp + nb301_ixH2]
    subps xmm7, [rsp + nb301_iyH2]
    subps xmm8, [rsp + nb301_izH2]    

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps [rsp + nb301_dxO], xmm0
	movaps [rsp + nb301_dyO], xmm1
	movaps [rsp + nb301_dzO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb301_dxH1], xmm3
	movaps [rsp + nb301_dyH1], xmm4
	movaps [rsp + nb301_dzH1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb301_dxH2], xmm6
	movaps [rsp + nb301_dyH2], xmm7
	movaps [rsp + nb301_dzH2], xmm8
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
		
	movaps  xmm9, [rsp + nb301_three]
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

	movaps  xmm4, [rsp + nb301_half]
	mulps   xmm9, xmm4  ;# rinvO 
	mulps   xmm10, xmm4 ;# rinvH1
    mulps   xmm11, xmm4 ;# rinvH2

	movaps  [rsp + nb301_rinvO], xmm9
	movaps  [rsp + nb301_rinvH1], xmm10
	movaps  [rsp + nb301_rinvH2], xmm11
	
	;# interactions 
    ;# rsq in xmm0,xmm3,xmm6  
    ;# rinv in xmm9, xmm10, xmm11

    movaps xmm1, [rsp + nb301_tsc]
    mulps  xmm0, xmm9  ;# r
    mulps  xmm3, xmm10
    mulps  xmm6, xmm11
    mulps  xmm0, xmm1 ;# rtab
    mulps  xmm3, xmm1
    mulps  xmm6, xmm1
    
    ;# truncate and convert to integers
    cvttps2dq xmm1, xmm0
    cvttps2dq xmm4, xmm3
    cvttps2dq xmm7, xmm6        

    ;# convert back to float
    cvtdq2ps  xmm2, xmm1
    cvtdq2ps  xmm5, xmm4
    cvtdq2ps  xmm8, xmm7
    
    ;# multiply by 4
    pslld   xmm1, 2
    pslld   xmm4, 2
    pslld   xmm7, 2
    
    ;# move to integer registers
    movhlps xmm13, xmm1
    movhlps xmm14, xmm4
    movhlps xmm15, xmm7
    movd    eax, xmm1
    movd    r8d, xmm4
    movd    r12d, xmm7
    movd    ecx, xmm13
    movd    r10d, xmm14
    movd    r14d, xmm15
    pshufd  xmm1, xmm1, 1
    pshufd  xmm4, xmm4, 1
    pshufd  xmm7, xmm7, 1
    pshufd  xmm13, xmm13, 1
    pshufd  xmm14, xmm14, 1
    pshufd  xmm15, xmm15, 1
    movd    ebx, xmm1
    movd    r9d, xmm4
    movd    r13d, xmm7    
    movd    edx, xmm13
    movd    r11d, xmm14
    movd    r15d, xmm15   
        
    mov  rsi, [rbp + nb301_VFtab]

    ;# calculate eps
    subps     xmm0, xmm2
    subps     xmm3, xmm5
    subps     xmm6, xmm8

    movaps    [rsp + nb301_epsO], xmm0
    movaps    [rsp + nb301_epsH1], xmm3
    movaps    [rsp + nb301_epsH2], xmm6

    ;# Load LOTS of table data
   	movlps xmm1, [rsi + rax*4]
   	movlps xmm5, [rsi + r8*4]
   	movlps xmm9, [rsi + r12*4]

	movlps xmm3, [rsi + rcx*4]
	movlps xmm7, [rsi + r10*4]
	movlps xmm11, [rsi + r14*4]

	movhps xmm1, [rsi + rbx*4]
	movhps xmm5, [rsi + r9*4]
	movhps xmm9, [rsi + r13*4]

	movhps xmm3, [rsi + rdx*4]
	movhps xmm7, [rsi + r11*4]
	movhps xmm11, [rsi + r15*4]

    movaps xmm0, xmm1
    movaps xmm4, xmm5
    movaps xmm8, xmm9
	shufps xmm0, xmm3, 136  ;# 10001000
	shufps xmm4, xmm7, 136  ;# 10001000
	shufps xmm8, xmm11, 136  ;# 10001000
	shufps xmm1, xmm3, 221  ;# 11011101
	shufps xmm5, xmm7, 221  ;# 11011101
	shufps xmm9, xmm11, 221  ;# 11011101
    
	movlps xmm3, [rsi + rax*4 + 8]
	movlps xmm7, [rsi + r8*4 + 8]
	movlps xmm11, [rsi + r12*4 + 8]
    
	movlps xmm12, [rsi + rcx*4 + 8]
	movlps xmm13, [rsi + r10*4 + 8]
	movlps xmm14, [rsi + r14*4 + 8]

	movhps xmm3, [rsi + rbx*4 + 8]
	movhps xmm7, [rsi + r9*4 + 8]
	movhps xmm11, [rsi + r13*4 + 8]
    
	movhps xmm12, [rsi + rdx*4 + 8]
	movhps xmm13, [rsi + r11*4 + 8]
	movhps xmm14, [rsi + r15*4 + 8]

    movaps xmm2, xmm3
    movaps xmm6, xmm7
    movaps xmm10, xmm11
    
	shufps xmm2, xmm12, 136  ;# 10001000
	shufps xmm6, xmm13, 136  ;# 10001000
	shufps xmm10, xmm14, 136  ;# 10001000
	shufps xmm3, xmm12, 221  ;# 11011101
	shufps xmm7, xmm13, 221  ;# 11011101
	shufps xmm11, xmm14, 221  ;# 11011101
    ;# table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11
    
    movaps xmm12, [rsp + nb301_epsO]
    movaps xmm13, [rsp + nb301_epsH1]
    movaps xmm14, [rsp + nb301_epsH2]
    
    mulps  xmm3, xmm12   ;# Heps
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 
    mulps  xmm2, xmm12   ;# Geps
    mulps  xmm6, xmm13
    mulps  xmm10, xmm14 
    mulps  xmm3, xmm12   ;# Heps2
    mulps  xmm7, xmm13
    mulps  xmm11, xmm14 

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
    mulps  xmm5, xmm13
    mulps  xmm9, xmm14
    movaps xmm12, [rsp + nb301_qqO]
    movaps xmm13, [rsp + nb301_qqH]
    addps  xmm1, xmm0     ;# VV
    addps  xmm5, xmm4
    addps  xmm9, xmm8
    mulps  xmm1, xmm12   ;# VV*qq = vcoul
    mulps  xmm5, xmm13
    mulps  xmm9, xmm13
    mulps  xmm3, xmm12    ;# FF*qq = fij
    mulps  xmm7, xmm13
    mulps  xmm11, xmm13
    
    ;# accumulate vctot
    addps  xmm1, [rsp + nb301_vctot]
    addps  xmm5, xmm9
    addps  xmm1, xmm5
    movaps [rsp + nb301_vctot], xmm1
    
    movaps xmm10, [rsp + nb301_tsc]
    mulps  xmm3, xmm10  ;# fscal
    mulps  xmm7, xmm10
    mulps  xmm10, xmm11
    
    movd eax, mm0
    movd ebx, mm1
    movd ecx, mm2
    movd edx, mm3

	;# move j forces to local temp variables 
    movlps xmm11, [rdi + rax*4] ;# jxa jya  -   -
    movlps xmm12, [rdi + rcx*4] ;# jxc jyc  -   -
    movhps xmm11, [rdi + rbx*4] ;# jxa jya jxb jyb 
    movhps xmm12, [rdi + rdx*4] ;# jxc jyc jxd jyd 

    movss  xmm13, [rdi + rax*4 + 8] ;# jza  -  -  -
    movss  xmm14, [rdi + rcx*4 + 8] ;# jzc  -  -  -
    movss  xmm2,  [rdi + rbx*4 + 8] ;# jzb
    movss  xmm5,  [rdi + rdx*4 + 8] ;# jzd
    movlhps xmm13, xmm2 ;# jza  -  jzb  -
    movlhps xmm14, xmm5 ;# jzc  -  jzd -
    
    shufps xmm13, xmm14,  136  ;# 10001000 => jza jzb jzc jzd

    ;# xmm11: jxa jya jxb jyb 
    ;# xmm12: jxc jyc jxd jyd
    ;# xmm13: jza jzb jzc jzd

    xorps  xmm0, xmm0
    xorps  xmm4, xmm4
    xorps  xmm8, xmm8

    mulps  xmm3, [rsp + nb301_rinvO]
    mulps  xmm7, [rsp + nb301_rinvH1]
    mulps  xmm10, [rsp + nb301_rinvH2]
    
    subps  xmm0, xmm3
    subps  xmm4, xmm7
    subps  xmm8, xmm10
    
    movaps xmm1, xmm0
    movaps xmm2, xmm0
    movaps xmm3, xmm4
    movaps xmm5, xmm4
    movaps xmm6, xmm8
    movaps xmm7, xmm8

	mulps xmm0, [rsp + nb301_dxO]
	mulps xmm1, [rsp + nb301_dyO]
	mulps xmm2, [rsp + nb301_dzO]
	mulps xmm3, [rsp + nb301_dxH1]
	mulps xmm4, [rsp + nb301_dyH1]
	mulps xmm5, [rsp + nb301_dzH1]
	mulps xmm6, [rsp + nb301_dxH2]
	mulps xmm7, [rsp + nb301_dyH2]
	mulps xmm8, [rsp + nb301_dzH2]

    movaps xmm14, xmm0
    movaps xmm15, xmm1
    addps xmm13,  xmm2
    addps xmm0, [rsp + nb301_fixO]
    addps xmm1, [rsp + nb301_fiyO]
    addps xmm2, [rsp + nb301_fizO]

    addps xmm14, xmm3
    addps xmm15, xmm4
    addps xmm13, xmm5
    addps xmm3, [rsp + nb301_fixH1]
    addps xmm4, [rsp + nb301_fiyH1]
    addps xmm5, [rsp + nb301_fizH1]

    addps xmm14, xmm6
    addps xmm15, xmm7
    addps xmm13, xmm8
    addps xmm6, [rsp + nb301_fixH2]
    addps xmm7, [rsp + nb301_fiyH2]
    addps xmm8, [rsp + nb301_fizH2]

    movaps [rsp + nb301_fixO], xmm0
    movaps [rsp + nb301_fiyO], xmm1
    movaps [rsp + nb301_fizO], xmm2
    movaps [rsp + nb301_fixH1], xmm3
    movaps [rsp + nb301_fiyH1], xmm4
    movaps [rsp + nb301_fizH1], xmm5
    movaps [rsp + nb301_fixH2], xmm6
    movaps [rsp + nb301_fiyH2], xmm7
    movaps [rsp + nb301_fizH2], xmm8
    
    ;# xmm14 = fjx
    ;# xmm15 = fjy
    ;# xmm13 = fjz
    movaps xmm0, xmm14
    unpcklps xmm14, xmm15
    unpckhps xmm0,  xmm15
    
    addps  xmm11, xmm14
    addps  xmm12, xmm0
    
    movhlps  xmm14, xmm13 ;# fjzc fjzd
    
    movlps [rdi + rax*4], xmm11
    movhps [rdi + rbx*4], xmm11
    movlps [rdi + rcx*4], xmm12
    movhps [rdi + rdx*4], xmm12
    movss  [rdi + rax*4 + 8], xmm13
    movss  [rdi + rcx*4 + 8], xmm14
    shufps xmm13, xmm13, 1
    shufps xmm14, xmm14, 1
    movss  [rdi + rbx*4 + 8], xmm13
    movss  [rdi + rdx*4 + 8], xmm14
    
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb301_innerk],  4
	jl    .nb301_odd_inner
	jmp   .nb301_unroll_loop
.nb301_odd_inner:	
	add dword ptr [rsp + nb301_innerk],  4
	jnz   .nb301_odd_loop
	jmp   .nb301_updateouterdata
.nb301_odd_loop:
	mov   rdx, [rsp + nb301_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb301_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb301_iqO]
	mov rsi, [rbp + nb301_charge] 
	movhps xmm4, [rsp + nb301_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb301_qqO], xmm3	;# use oxygen qq for storage 

	mov rsi, [rbp + nb301_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rax*4 + 4]
	movss xmm5, [rsi + rax*4 + 8]
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	
	movss xmm0, [rsp + nb301_ixO]
	movss xmm1, [rsp + nb301_iyO]
	movss xmm2, [rsp + nb301_izO]
	
	movlps xmm6, [rsp + nb301_ixH1]
	movlps xmm7, [rsp + nb301_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm0, xmm6
	movlps xmm6, [rsp + nb301_iyH1]
	movlps xmm7, [rsp + nb301_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm1, xmm6
	movlps xmm6, [rsp + nb301_izH1]
	movlps xmm7, [rsp + nb301_izH2]
	unpcklps xmm6, xmm7
	movlhps xmm2, xmm6

	subps xmm3, xmm0
	subps xmm4, xmm1
	subps xmm5, xmm2
	
	movaps [rsp + nb301_dxO], xmm3
	movaps [rsp + nb301_dyO], xmm4
	movaps [rsp + nb301_dzO], xmm5

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
	movaps xmm1, [rsp + nb301_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb301_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 
	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000	
	
	mulps xmm4, xmm0	;# xmm4=r 
	movaps [rsp + nb301_rinvO], xmm0
	
	mulps xmm4, [rsp + nb301_tsc]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7    ;# mm6/mm7 contain lu indices 
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm7, mm7
    	movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
	
    	movd mm0, eax   
    	movd mm1, ecx
    	movd mm2, edx

    	mov  rsi, [rbp + nb301_VFtab]
    	movd eax, mm6
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd edx, mm7

    	movlps xmm5, [rsi + rax*4]
    	movlps xmm7, [rsi + rcx*4]
    	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# 10001000
    	shufps xmm5, xmm7, 221  ;# 11011101

    	movlps xmm7, [rsi + rax*4 + 8]
    	movlps xmm3, [rsi + rcx*4 + 8]
    	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# 10001000
    	shufps xmm7, xmm3, 221  ;# 11011101
    	;# coulomb table ready, in xmm4-xmm7      
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	mulps  xmm7, [rsp + nb301_two]       ;# two*Heps2 
    	movaps xmm0, [rsp + nb301_qqO]
    	addps  xmm7, xmm6
    	addps  xmm7, xmm5 ;# xmm7=FF 
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	mulps  xmm0, xmm7 ;# fijC=FF*qq 
    	;# at this point mm5 contains vcoul and xmm0 fijC 
    	;# increment vcoul - then we can get rid of mm5 
    	addps  xmm5, [rsp + nb301_vctot]
    	movaps [rsp + nb301_vctot], xmm5

	xorps xmm4, xmm4
	mulps  xmm0, [rsp + nb301_tsc]
	mulps  xmm0, [rsp + nb301_rinvO]	
	subps  xmm4, xmm0
	
    	movd eax, mm0   
    	movd ecx, mm1
    	movd edx, mm2	
	
	movaps xmm0, [rsp + nb301_dxO]
	movaps xmm1, [rsp + nb301_dyO]
	movaps xmm2, [rsp + nb301_dzO]

	mulps  xmm0, xmm4
	mulps  xmm1, xmm4
	mulps  xmm2, xmm4 ;# xmm0-xmm2 now contains tx-tz (partial force) 
	movss  xmm3, [rsp + nb301_fixO]	
	movss  xmm4, [rsp + nb301_fiyO]	
	movss  xmm5, [rsp + nb301_fizO]	
	addss  xmm3, xmm0
	addss  xmm4, xmm1
	addss  xmm5, xmm2
	movss  [rsp + nb301_fixO], xmm3	
	movss  [rsp + nb301_fiyO], xmm4	
	movss  [rsp + nb301_fizO], xmm5	;# updated the O force now do the H's 
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	shufps xmm3, xmm3, 230 ;# 11100110	;# shift right 
	shufps xmm4, xmm4, 230 ;# 11100110
	shufps xmm5, xmm5, 230 ;# 11100110
	addss  xmm3, [rsp + nb301_fixH1]
	addss  xmm4, [rsp + nb301_fiyH1]
	addss  xmm5, [rsp + nb301_fizH1]
	movss  [rsp + nb301_fixH1], xmm3	
	movss  [rsp + nb301_fiyH1], xmm4	
	movss  [rsp + nb301_fizH1], xmm5	;# updated the H1 force 

	mov rdi, [rbp + nb301_faction]
	shufps xmm3, xmm3, 231 ;# 11100111	;# shift right 
	shufps xmm4, xmm4, 231 ;# 11100111
	shufps xmm5, xmm5, 231 ;# 11100111
	addss  xmm3, [rsp + nb301_fixH2]
	addss  xmm4, [rsp + nb301_fiyH2]
	addss  xmm5, [rsp + nb301_fizH2]
	movss  [rsp + nb301_fixH2], xmm3	
	movss  [rsp + nb301_fiyH2], xmm4	
	movss  [rsp + nb301_fizH2], xmm5	;# updated the H2 force 

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
	addss    xmm2, xmm1    ;# z sum in xmm2 
	addps    xmm6, xmm0
	addss    xmm7, xmm2
	
	movlps [rdi + rax*4],     xmm6
	movss  [rdi + rax*4 + 8], xmm7

	dec dword ptr [rsp + nb301_innerk]
	jz    .nb301_updateouterdata
	jmp   .nb301_odd_loop
.nb301_updateouterdata:
	mov   ecx, [rsp + nb301_ii3]
	mov   rdi, [rbp + nb301_faction]
	mov   rsi, [rbp + nb301_fshift]
	mov   edx, [rsp + nb301_is3]

	;# accumulate  Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb301_fixO]
	movaps xmm1, [rsp + nb301_fiyO]
	movaps xmm2, [rsp + nb301_fizO]

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

	;# accumulate force in xmm6/xmm7 for fshift 
	movaps xmm6, xmm0
	movss xmm7, xmm2
	movlhps xmm6, xmm1
	shufps  xmm6, xmm6, 8 ;# 00001000	

	;# accumulate H1i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb301_fixH1]
	movaps xmm1, [rsp + nb301_fiyH1]
	movaps xmm2, [rsp + nb301_fizH1]

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
	addss xmm7, xmm2
	movlhps xmm0, xmm1
	shufps  xmm0, xmm0, 8 ;# 00001000	
	addps   xmm6, xmm0

	;# accumulate H2i forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb301_fixH2]
	movaps xmm1, [rsp + nb301_fiyH2]
	movaps xmm2, [rsp + nb301_fizH2]

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

	;# increment fshift force  
	movlps  xmm3, [rsi + rdx*4]
	movss  xmm4, [rsi + rdx*4 + 8]
	subps  xmm3, xmm6
	subss  xmm4, xmm7
	movlps  [rsi + rdx*4],    xmm3
	movss  [rsi + rdx*4 + 8], xmm4

	;# get n from stack
	mov esi, [rsp + nb301_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb301_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb301_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb301_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 

        ;# finish if last 
        mov ecx, [rsp + nb301_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb301_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb301_n], esi
        jmp .nb301_outer
.nb301_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb301_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb301_end
        ;# non-zero, do one more workunit
        jmp   .nb301_threadloop
.nb301_end:

	mov eax, [rsp + nb301_nouter]
	mov ebx, [rsp + nb301_ninner]
	mov rcx, [rbp + nb301_outeriter]
	mov rdx, [rbp + nb301_inneriter]
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




.globl nb_kernel301nf_x86_64_sse
.globl _nb_kernel301nf_x86_64_sse
nb_kernel301nf_x86_64_sse:	
_nb_kernel301nf_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb301nf_fshift,         16
.equiv          nb301nf_gid,            24
.equiv          nb301nf_pos,            32
.equiv          nb301nf_faction,        40
.equiv          nb301nf_charge,         48
.equiv          nb301nf_p_facel,        56
.equiv          nb301nf_argkrf,         64
.equiv          nb301nf_argcrf,         72
.equiv          nb301nf_Vc,             80
.equiv          nb301nf_type,           88
.equiv          nb301nf_p_ntype,        96
.equiv          nb301nf_vdwparam,       104
.equiv          nb301nf_Vvdw,           112
.equiv          nb301nf_p_tabscale,     120
.equiv          nb301nf_VFtab,          128
.equiv          nb301nf_invsqrta,       136
.equiv          nb301nf_dvda,           144
.equiv          nb301nf_p_gbtabscale,   152
.equiv          nb301nf_GBtab,          160
.equiv          nb301nf_p_nthreads,     168
.equiv          nb301nf_count,          176
.equiv          nb301nf_mtx,            184
.equiv          nb301nf_outeriter,      192
.equiv          nb301nf_inneriter,      200
.equiv          nb301nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 
.equiv          nb301nf_ixO,            0
.equiv          nb301nf_iyO,            16
.equiv          nb301nf_izO,            32
.equiv          nb301nf_ixH1,           48
.equiv          nb301nf_iyH1,           64
.equiv          nb301nf_izH1,           80
.equiv          nb301nf_ixH2,           96
.equiv          nb301nf_iyH2,           112
.equiv          nb301nf_izH2,           128
.equiv          nb301nf_iqO,            144
.equiv          nb301nf_iqH,            160
.equiv          nb301nf_qqO,            176
.equiv          nb301nf_qqH,            192
.equiv          nb301nf_rinvO,          208
.equiv          nb301nf_rinvH1,         224
.equiv          nb301nf_rinvH2,         240
.equiv          nb301nf_rO,             256
.equiv          nb301nf_rH1,            272
.equiv          nb301nf_rH2,            288
.equiv          nb301nf_tsc,            304
.equiv          nb301nf_vctot,          320
.equiv          nb301nf_half,           336
.equiv          nb301nf_three,          352
.equiv          nb301nf_is3,            368
.equiv          nb301nf_ii3,            372
.equiv          nb301nf_nri,            376
.equiv          nb301nf_iinr,           384
.equiv          nb301nf_jindex,         392
.equiv          nb301nf_jjnr,           400
.equiv          nb301nf_shift,          408
.equiv          nb301nf_shiftvec,       416
.equiv          nb301nf_facel,          424
.equiv          nb301nf_innerjjnr,      432
.equiv          nb301nf_innerk,         440
.equiv          nb301nf_n,              444
.equiv          nb301nf_nn1,            448
.equiv          nb301nf_nouter,         452
.equiv          nb301nf_ninner,         456

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
	mov [rsp + nb301nf_nouter], eax
	mov [rsp + nb301nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb301nf_nri], edi
	mov [rsp + nb301nf_iinr], rsi
	mov [rsp + nb301nf_jindex], rdx
	mov [rsp + nb301nf_jjnr], rcx
	mov [rsp + nb301nf_shift], r8
	mov [rsp + nb301nf_shiftvec], r9
	mov rsi, [rbp + nb301nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb301nf_facel], xmm0

	mov rax, [rbp + nb301nf_p_tabscale]
	movss xmm3, [rax]
	shufps xmm3, xmm3, 0
	movaps [rsp + nb301nf_tsc], xmm3

	
	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb301nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb301nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb301nf_p_facel]
	movss xmm0, [rsi]
	movss xmm5, [rsp + nb301nf_facel]
	mulss  xmm3, xmm5
	mulss  xmm4, xmm5

	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	movaps [rsp + nb301nf_iqO], xmm3
	movaps [rsp + nb301nf_iqH], xmm4
	
	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb301nf_half], eax
	movss xmm1, [rsp + nb301nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb301nf_half],  xmm1
	movaps [rsp + nb301nf_three],  xmm3

.nb301nf_threadloop:
        mov   rsi, [rbp + nb301nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb301nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb301nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb301nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb301nf_n], eax
        mov [rsp + nb301nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb301nf_outerstart
        jmp .nb301nf_end
	
.nb301nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb301nf_nouter]
	mov [rsp + nb301nf_nouter], ebx

.nb301nf_outer:
	mov   rax, [rsp + nb301nf_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax + rsi*4]		;# rbx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb301nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb301nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb301nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb301nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb301nf_ii3], ebx

	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb301nf_ixO], xmm3
	movaps [rsp + nb301nf_iyO], xmm4
	movaps [rsp + nb301nf_izO], xmm5

	movss xmm3, xmm0
	movss xmm4, xmm1
	movss xmm5, xmm2
	addss xmm0, [rax + rbx*4 + 12]
	addss xmm1, [rax + rbx*4 + 16]
	addss xmm2, [rax + rbx*4 + 20]		
	addss xmm3, [rax + rbx*4 + 24]
	addss xmm4, [rax + rbx*4 + 28]
	addss xmm5, [rax + rbx*4 + 32]		

	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb301nf_ixH1], xmm0
	movaps [rsp + nb301nf_iyH1], xmm1
	movaps [rsp + nb301nf_izH1], xmm2
	movaps [rsp + nb301nf_ixH2], xmm3
	movaps [rsp + nb301nf_iyH2], xmm4
	movaps [rsp + nb301nf_izH2], xmm5
	
	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb301nf_vctot], xmm4
	
	mov   rax, [rsp + nb301nf_jindex]
	mov   ecx, [rax + rsi*4]	     ;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     ;# jindex[n+1] 
	sub   edx, ecx               ;# number of innerloop atoms 

	mov   rsi, [rbp + nb301nf_pos]
	mov   rax, [rsp + nb301nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb301nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb301nf_ninner]
	mov   [rsp + nb301nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb301nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb301nf_unroll_loop
	jmp   .nb301nf_odd_inner
.nb301nf_unroll_loop:
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb301nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	mov   ebx, [rdx + 4]              
	mov   ecx, [rdx + 8]            
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 

	add qword ptr [rsp + nb301nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb301nf_charge]    ;# base of charge[] 
	
	movss xmm3, [rsi + rax*4]
	movss xmm4, [rsi + rcx*4]
	movss xmm6, [rsi + rbx*4]
	movss xmm7, [rsi + rdx*4]

	shufps xmm3, xmm6, 0 
	shufps xmm4, xmm7, 0 
	shufps xmm3, xmm4, 136  ;# 10001000 ;# all charges in xmm3  
	movaps xmm4, xmm3	     ;# and in xmm4 
	mulps  xmm3, [rsp + nb301nf_iqO]
	mulps  xmm4, [rsp + nb301nf_iqH]

	movd  mm0, eax		;# use mmx registers as temp storage 
	movd  mm1, ebx
	movd  mm2, ecx
	movd  mm3, edx

	movaps  [rsp + nb301nf_qqO], xmm3
	movaps  [rsp + nb301nf_qqH], xmm4	

	mov rsi, [rbp + nb301nf_pos]       ;# base of pos[] 

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

	;# move ixO-izO to xmm4-xmm6 
	movaps xmm4, [rsp + nb301nf_ixO]
	movaps xmm5, [rsp + nb301nf_iyO]
	movaps xmm6, [rsp + nb301nf_izO]

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

	;# move ixH1-izH1 to xmm4-xmm6 
	movaps xmm4, [rsp + nb301nf_ixH1]
	movaps xmm5, [rsp + nb301nf_iyH1]
	movaps xmm6, [rsp + nb301nf_izH1]

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

	;# move ixH2-izH2 to xmm3-xmm5  
	movaps xmm3, [rsp + nb301nf_ixH2]
	movaps xmm4, [rsp + nb301nf_iyH2]
	movaps xmm5, [rsp + nb301nf_izH2]

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
	;# rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

	;# start with rsqO - seed to xmm2 	
	rsqrtps xmm2, xmm7
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb301nf_three]
	mulps   xmm2, xmm7	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb301nf_half]
	movaps  [rsp + nb301nf_rinvO], xmm4	;# rinvO in xmm4 
	mulps   xmm7, xmm4
	movaps  [rsp + nb301nf_rO], xmm7	

	;# rsqH1 - seed in xmm2 
	rsqrtps xmm2, xmm6
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb301nf_three]
	mulps   xmm2, xmm6	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb301nf_half]
	movaps  [rsp + nb301nf_rinvH1], xmm4	;# rinvH1 in xmm4 
	mulps   xmm6, xmm4
	movaps  [rsp + nb301nf_rH1], xmm6

	;# rsqH2 - seed to xmm2 
	rsqrtps xmm2, xmm5
	movaps  xmm3, xmm2
	mulps   xmm2, xmm2
	movaps  xmm4, [rsp + nb301nf_three]
	mulps   xmm2, xmm5	;# rsq*lu*lu 
	subps   xmm4, xmm2	;# 30-rsq*lu*lu 
	mulps   xmm4, xmm3	;# lu*(3-rsq*lu*lu) 
	mulps   xmm4, [rsp + nb301nf_half]
	movaps  [rsp + nb301nf_rinvH2], xmm4	;# rinvH2 in xmm4 
	mulps   xmm5, xmm4
	movaps  [rsp + nb301nf_rH2], xmm5

	;# do O interactions 
	;# rO is still in xmm7 
	mulps   xmm7, [rsp + nb301nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3

	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
	
    	movd mm0, eax   
    	movd mm1, ebx
    	movd mm2, ecx
    	movd mm3, edx

    	mov  rsi, [rbp + nb301nf_VFtab]
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
        
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp    
    	movaps xmm0, [rsp + nb301nf_qqO]
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  

    	;# at this point mm5 contains vcoul 
    	;# increment vcoul - then we can get rid of mm5 
    	addps  xmm5, [rsp + nb301nf_vctot]
    	movaps [rsp + nb301nf_vctot], xmm5
	
	;# Done with O interactions - now H1! 
	movaps xmm7, [rsp + nb301nf_rH1]
	mulps   xmm7, [rsp + nb301nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
	
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
        
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	movaps xmm0, [rsp + nb301nf_qqH]
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	;# at this point mm5 contains vcoul 
    	;# increment vcoul 
    	addps  xmm5, [rsp + nb301nf_vctot]
    	movaps [rsp + nb301nf_vctot], xmm5 

	;# Done with H1, finally we do H2 interactions 
	movaps xmm7, [rsp + nb301nf_rH2]
	mulps   xmm7, [rsp + nb301nf_tsc]
	movhlps xmm4, xmm7
	cvttps2pi mm6, xmm7
	cvttps2pi mm7, xmm4    ;# mm6/mm7 contain lu indices 
	
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm4, mm7
    	movlhps xmm3, xmm4
	
    	subps xmm7, xmm3
	movaps xmm1, xmm7	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
	
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
        
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	movaps xmm0, [rsp + nb301nf_qqH]
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	;# at this point mm5 contains vcoul 
    	;# increment vcoul 
    	addps  xmm5, [rsp + nb301nf_vctot]
    	movaps [rsp + nb301nf_vctot], xmm5
	
	;# should we do one more iteration? 
	sub dword ptr [rsp + nb301nf_innerk],  4
	jl    .nb301nf_odd_inner
	jmp   .nb301nf_unroll_loop
.nb301nf_odd_inner:	
	add dword ptr [rsp + nb301nf_innerk],  4
	jnz   .nb301nf_odd_loop
	jmp   .nb301nf_updateouterdata
.nb301nf_odd_loop:
	mov   rdx, [rsp + nb301nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb301nf_innerjjnr],  4	

 	xorps xmm4, xmm4
	movss xmm4, [rsp + nb301nf_iqO]
	mov rsi, [rbp + nb301nf_charge] 
	movhps xmm4, [rsp + nb301nf_iqH]     
	movss xmm3, [rsi + rax*4]	;# charge in xmm3 
	shufps xmm3, xmm3, 0
	mulps xmm3, xmm4
	movaps [rsp + nb301nf_qqO], xmm3	;# use oxygen qq for storage 

	mov rsi, [rbp + nb301nf_pos]
	lea   rax, [rax + rax*2]  
	
	;# move j coords to xmm0-xmm2 
	movss xmm0, [rsi + rax*4]
	movss xmm1, [rsi + rax*4 + 4]
	movss xmm2, [rsi + rax*4 + 8]
	shufps xmm0, xmm0, 0
	shufps xmm1, xmm1, 0
	shufps xmm2, xmm2, 0
	
	movss xmm3, [rsp + nb301nf_ixO]
	movss xmm4, [rsp + nb301nf_iyO]
	movss xmm5, [rsp + nb301nf_izO]
	
	movlps xmm6, [rsp + nb301nf_ixH1]
	movlps xmm7, [rsp + nb301nf_ixH2]
	unpcklps xmm6, xmm7
	movlhps xmm3, xmm6
	movlps xmm6, [rsp + nb301nf_iyH1]
	movlps xmm7, [rsp + nb301nf_iyH2]
	unpcklps xmm6, xmm7
	movlhps xmm4, xmm6
	movlps xmm6, [rsp + nb301nf_izH1]
	movlps xmm7, [rsp + nb301nf_izH2]
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
	movaps xmm1, [rsp + nb301nf_three]
	mulps xmm5, xmm4	;# rsq*lu*lu 			
	movaps xmm0, [rsp + nb301nf_half]
	subps xmm1, xmm5	;# 30-rsq*lu*lu 
	mulps xmm1, xmm2	
	mulps xmm0, xmm1	;# xmm0=rinv 

	;# a little trick to avoid NaNs: 
	;# positions 0,2,and 3 are valid, but not 1. 
	;# If it contains NaN it doesnt help to mult by 0, 
	;# So we shuffle it and copy pos 0 to pos1! 
	shufps xmm0, xmm0, 224 ;# 11100000	
	
	mulps xmm4, xmm0	;# xmm4=r 
	movaps [rsp + nb301nf_rinvO], xmm0
	
	mulps xmm4, [rsp + nb301nf_tsc]
	movhlps xmm7, xmm4
	cvttps2pi mm6, xmm4
	cvttps2pi mm7, xmm7    ;# mm6/mm7 contain lu indices 
    	cvtpi2ps xmm3, mm6
    	cvtpi2ps xmm7, mm7
    	movlhps xmm3, xmm7

	subps   xmm4, xmm3	
	movaps xmm1, xmm4	;# xmm1=eps 
	movaps xmm2, xmm1
	mulps  xmm2, xmm2	;# xmm2=eps2 
    	pslld mm6, 2
    	pslld mm7, 2
	
    	movd mm0, eax   
    	movd mm1, ecx
    	movd mm2, edx

    	mov  rsi, [rbp + nb301nf_VFtab]
    	movd eax, mm6
    	movd ecx, mm7
    	psrlq mm7, 32
    	movd edx, mm7

    	movlps xmm5, [rsi + rax*4]
    	movlps xmm7, [rsi + rcx*4]
    	movhps xmm7, [rsi + rdx*4] ;# got half coulomb table 

    	movaps xmm4, xmm5
    	shufps xmm4, xmm7, 136  ;# 10001000
    	shufps xmm5, xmm7, 221  ;# 11011101

    	movlps xmm7, [rsi + rax*4 + 8]
    	movlps xmm3, [rsi + rcx*4 + 8]
    	movhps xmm3, [rsi + rdx*4 + 8] ;# other half of coulomb table  
    	movaps xmm6, xmm7
    	shufps xmm6, xmm3, 136  ;# 10001000
    	shufps xmm7, xmm3, 221  ;# 11011101
    	;# coulomb table ready, in xmm4-xmm7      
    	mulps  xmm6, xmm1       ;# xmm6=Geps 
    	mulps  xmm7, xmm2       ;# xmm7=Heps2 
    	addps  xmm5, xmm6
    	addps  xmm5, xmm7       ;# xmm5=Fp        
    	movaps xmm0, [rsp + nb301nf_qqO]
    	mulps  xmm5, xmm1 ;# xmm5=eps*Fp 
    	addps  xmm5, xmm4 ;# xmm5=VV 
    	mulps  xmm5, xmm0 ;# vcoul=qq*VV  
    	;# at this point mm5 contains vcoul 
    	;# increment vcoul - then we can get rid of mm5 
    	addps  xmm5, [rsp + nb301nf_vctot]
    	movaps [rsp + nb301nf_vctot], xmm5

	dec dword ptr [rsp + nb301nf_innerk]
	jz    .nb301nf_updateouterdata
	jmp   .nb301nf_odd_loop
.nb301nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb301nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb301nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb301nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		
        
	;# add earlier value from mem 
	mov   rax, [rbp + nb301nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 

        ;# finish if last 
        mov ecx, [rsp + nb301nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb301nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb301nf_n], esi
        jmp .nb301nf_outer
.nb301nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb301nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb301nf_end
        ;# non-zero, do one more workunit
        jmp   .nb301nf_threadloop
.nb301nf_end:

	mov eax, [rsp + nb301nf_nouter]
	mov ebx, [rsp + nb301nf_ninner]
	mov rcx, [rbp + nb301nf_outeriter]
	mov rdx, [rbp + nb301nf_inneriter]
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
