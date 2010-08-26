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



	
.globl nb_kernel102_x86_64_sse
.globl _nb_kernel102_x86_64_sse
nb_kernel102_x86_64_sse:	
_nb_kernel102_x86_64_sse:	
;#	Room for return address and rbp (16 bytes)
.equiv          nb102_fshift,           16
.equiv          nb102_gid,              24
.equiv          nb102_pos,              32
.equiv          nb102_faction,          40
.equiv          nb102_charge,           48
.equiv          nb102_p_facel,          56
.equiv          nb102_argkrf,           64
.equiv          nb102_argcrf,           72
.equiv          nb102_Vc,               80
.equiv          nb102_type,             88
.equiv          nb102_p_ntype,          96
.equiv          nb102_vdwparam,         104
.equiv          nb102_Vvdw,             112
.equiv          nb102_p_tabscale,       120
.equiv          nb102_VFtab,            128
.equiv          nb102_invsqrta,         136
.equiv          nb102_dvda,             144
.equiv          nb102_p_gbtabscale,     152
.equiv          nb102_GBtab,            160
.equiv          nb102_p_nthreads,       168
.equiv          nb102_count,            176
.equiv          nb102_mtx,              184
.equiv          nb102_outeriter,        192
.equiv          nb102_inneriter,        200
.equiv          nb102_work,             208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 	
.equiv          nb102_ixO,              0
.equiv          nb102_iyO,              16
.equiv          nb102_izO,              32
.equiv          nb102_ixH1,             48
.equiv          nb102_iyH1,             64
.equiv          nb102_izH1,             80
.equiv          nb102_ixH2,             96
.equiv          nb102_iyH2,             112
.equiv          nb102_izH2,             128
.equiv          nb102_jxO,              144
.equiv          nb102_jyO,              160
.equiv          nb102_jzO,              176
.equiv          nb102_jxH1,             192
.equiv          nb102_jyH1,             208
.equiv          nb102_jzH1,             224
.equiv          nb102_jxH2,             240
.equiv          nb102_jyH2,             256
.equiv          nb102_jzH2,             272
.equiv          nb102_dxOO,             288
.equiv          nb102_dyOO,             304
.equiv          nb102_dzOO,             320
.equiv          nb102_dxOH1,            336
.equiv          nb102_dyOH1,            352
.equiv          nb102_dzOH1,            368
.equiv          nb102_dxOH2,            384
.equiv          nb102_dyOH2,            400
.equiv          nb102_dzOH2,            416
.equiv          nb102_dxH1O,            432
.equiv          nb102_dyH1O,            448
.equiv          nb102_dzH1O,            464
.equiv          nb102_dxH1H1,           480
.equiv          nb102_dyH1H1,           496
.equiv          nb102_dzH1H1,           512
.equiv          nb102_dxH1H2,           528
.equiv          nb102_dyH1H2,           544
.equiv          nb102_dzH1H2,           560
.equiv          nb102_dxH2O,            576
.equiv          nb102_dyH2O,            592
.equiv          nb102_dzH2O,            608
.equiv          nb102_dxH2H1,           624
.equiv          nb102_dyH2H1,           640
.equiv          nb102_dzH2H1,           656
.equiv          nb102_dxH2H2,           672
.equiv          nb102_dyH2H2,           688
.equiv          nb102_dzH2H2,           704
.equiv          nb102_qqOO,             720
.equiv          nb102_qqOH,             736
.equiv          nb102_qqHH,             752
.equiv          nb102_vctot,            768
.equiv          nb102_fixO,             784
.equiv          nb102_fiyO,             800
.equiv          nb102_fizO,             816
.equiv          nb102_fixH1,            832
.equiv          nb102_fiyH1,            848
.equiv          nb102_fizH1,            864
.equiv          nb102_fixH2,            880
.equiv          nb102_fiyH2,            896
.equiv          nb102_fizH2,            912
.equiv          nb102_fjxO,             928
.equiv          nb102_fjyO,             944
.equiv          nb102_fjzO,             960
.equiv          nb102_fjxH1,            976
.equiv          nb102_fjyH1,            992
.equiv          nb102_fjzH1,            1008
.equiv          nb102_fjxH2,            1024
.equiv          nb102_fjyH2,            1040
.equiv          nb102_fjzH2,            1056
.equiv          nb102_half,             1072
.equiv          nb102_three,            1088
.equiv          nb102_rsqOO,            1104
.equiv          nb102_rsqOH1,           1120
.equiv          nb102_rsqOH2,           1136
.equiv          nb102_rsqH1O,           1152
.equiv          nb102_rsqH1H1,          1168
.equiv          nb102_rsqH1H2,          1184
.equiv          nb102_rsqH2O,           1200
.equiv          nb102_rsqH2H1,          1216
.equiv          nb102_rsqH2H2,          1232
.equiv          nb102_rinvOO,           1248
.equiv          nb102_rinvOH1,          1264
.equiv          nb102_rinvOH2,          1280
.equiv          nb102_rinvH1O,          1296
.equiv          nb102_rinvH1H1,         1312
.equiv          nb102_rinvH1H2,         1328
.equiv          nb102_rinvH2O,          1344
.equiv          nb102_rinvH2H1,         1360
.equiv          nb102_rinvH2H2,         1376
.equiv          nb102_is3,              1392
.equiv          nb102_ii3,              1396
.equiv          nb102_innerjjnr,        1400
.equiv          nb102_nri,              1408
.equiv          nb102_iinr,             1416
.equiv          nb102_jindex,           1424
.equiv          nb102_jjnr,             1432
.equiv          nb102_shift,            1440
.equiv          nb102_shiftvec,         1448
.equiv          nb102_facel,            1456
.equiv          nb102_innerk,           1464
.equiv          nb102_n,                1472
.equiv          nb102_nn1,              1480
.equiv          nb102_nouter,           1484
.equiv          nb102_ninner,           1488
.equiv          nb102_salign,           1492

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
	sub rsp, 1504
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
	mov [rsp + nb102_nouter], eax
	mov [rsp + nb102_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb102_nri], edi
	mov [rsp + nb102_iinr], rsi
	mov [rsp + nb102_jindex], rdx
	mov [rsp + nb102_jjnr], rcx
	mov [rsp + nb102_shift], r8
	mov [rsp + nb102_shiftvec], r9
	mov rsi, [rbp + nb102_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb102_facel], xmm0	

	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb102_iinr]       ;# ecx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb102_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb102_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb102_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb102_qqOO], xmm3
	movaps [rsp + nb102_qqOH], xmm4
	movaps [rsp + nb102_qqHH], xmm5

	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb102_half], eax
	movss xmm1, [rsp + nb102_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb102_half],  xmm1
	movaps [rsp + nb102_three],  xmm3
	mov rdi, [rbp + nb102_faction]


.nb102_threadloop:
        mov   rsi, [rbp + nb102_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb102_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb102_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb102_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb102_n], eax
        mov [rsp + nb102_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi

        jg  .nb102_outerstart
        jmp .nb102_end
	
.nb102_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb102_nouter]
	mov [rsp + nb102_nouter], ebx

.nb102_outer:
	mov   rax, [rsp + nb102_shift]      ;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]			;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb102_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb102_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb102_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx+rsi*4]	    		;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb102_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb102_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb102_ixO], xmm3
	movaps [rsp + nb102_iyO], xmm4
	movaps [rsp + nb102_izO], xmm5

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
	movaps [rsp + nb102_ixH1], xmm0
	movaps [rsp + nb102_iyH1], xmm1
	movaps [rsp + nb102_izH1], xmm2
	movaps [rsp + nb102_ixH2], xmm3
	movaps [rsp + nb102_iyH2], xmm4
	movaps [rsp + nb102_izH2], xmm5

	;# clear vctot and i forces 
	xorps xmm4, xmm4
	movaps [rsp + nb102_vctot], xmm4
	movaps [rsp + nb102_fixO], xmm4
	movaps [rsp + nb102_fiyO], xmm4
	movaps [rsp + nb102_fizO], xmm4
	movaps [rsp + nb102_fixH1], xmm4
	movaps [rsp + nb102_fiyH1], xmm4
	movaps [rsp + nb102_fizH1], xmm4
	movaps [rsp + nb102_fixH2], xmm4
	movaps [rsp + nb102_fiyH2], xmm4
	movaps [rsp + nb102_fizH2], xmm4
	
	mov   rax, [rsp + nb102_jindex]
	mov   ecx, [rax +rsi*4]	     		;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     	;# jindex[n+1] 
	sub   edx, ecx               		;# number of innerloop atoms 

	mov   rsi, [rbp + nb102_pos]
	mov   rdi, [rbp + nb102_faction]	
	mov   rax, [rsp + nb102_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb102_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb102_ninner]
	mov   [rsp + nb102_ninner], ecx
	add   edx, 0
	mov   [rsp + nb102_innerk], edx    ;# number of innerloop atoms 
	jge   .nb102_unroll_loop
	jmp   .nb102_single_check
.nb102_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb102_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb102_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb102_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j O coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm1, [rsi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm0, [rsi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm1, [rsi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm2, [rsi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm3, [rsi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm2, [rsi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm3, [rsi + rdx*4 + 8] ;# jzOc  -  jzOd -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxOa jxOc jyOa jyOc        
    unpckhps xmm4, xmm1  ;# jxOb jxOd jyOb jyOd
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps xmm2, xmm3,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm0 = Ox
    ;# xmm1 = Oy
    ;# xmm2 = Oz
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    
    subps xmm0, [rsp + nb102_ixO]
    subps xmm1, [rsp + nb102_iyO]
    subps xmm2, [rsp + nb102_izO]
    subps xmm3, [rsp + nb102_ixH1]
    subps xmm4, [rsp + nb102_iyH1]
    subps xmm5, [rsp + nb102_izH1]
    subps xmm6, [rsp + nb102_ixH2]
    subps xmm7, [rsp + nb102_iyH2]
    subps xmm8, [rsp + nb102_izH2]
    
	movaps [rsp + nb102_dxOO], xmm0
	movaps [rsp + nb102_dyOO], xmm1
	movaps [rsp + nb102_dzOO], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb102_dxH1O], xmm3
	movaps [rsp + nb102_dyH1O], xmm4
	movaps [rsp + nb102_dzH1O], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb102_dxH2O], xmm6
	movaps [rsp + nb102_dyH2O], xmm7
	movaps [rsp + nb102_dzH2O], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jO atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb102_three]
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

	movaps  xmm0, [rsp + nb102_half]
	mulps   xmm9, xmm0  ;# rinvOO 
	mulps   xmm10, xmm0 ;# rinvH1O
    mulps   xmm11, xmm0 ;# rinvH2O
	
	;# O interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb102_qqOO] 
    mulps  xmm1, [rsp + nb102_qqOH] 
    mulps  xmm2, [rsp + nb102_qqOH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb102_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb102_vctot], xmm0
    
	;# move j O forces to local temp variables 
    movlps xmm0, [rdi + rax*4] ;# jxOa jyOa  -   -
    movlps xmm1, [rdi + rcx*4] ;# jxOc jyOc  -   -
    movhps xmm0, [rdi + rbx*4] ;# jxOa jyOa jxOb jyOb 
    movhps xmm1, [rdi + rdx*4] ;# jxOc jyOc jxOd jyOd 

    movss  xmm2, [rdi + rax*4 + 8] ;# jzOa  -  -  -
    movss  xmm3, [rdi + rcx*4 + 8] ;# jzOc  -  -  -
    movhps xmm2, [rdi + rbx*4 + 8] ;# jzOa  -  jzOb  -
    movhps xmm3, [rdi + rdx*4 + 8] ;# jzOc  -  jzOd -

    shufps xmm2, xmm3,  136  ;# 10001000 => jzOa jzOb jzOc jzOd

    ;# xmm0: jxOa jyOa jxOb jyOb 
    ;# xmm1: jxOc jyOc jxOd jyOd
    ;# xmm2: jzOa jzOb jzOc jzOd

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb102_dxOO]
	mulps xmm8, [rsp + nb102_dyOO]
	mulps xmm9, [rsp + nb102_dzOO]
	mulps xmm10, [rsp + nb102_dxH1O]
	mulps xmm11, [rsp + nb102_dyH1O]
	mulps xmm12, [rsp + nb102_dzH1O]
	mulps xmm13, [rsp + nb102_dxH2O]
	mulps xmm14, [rsp + nb102_dyH2O]
	mulps xmm15, [rsp + nb102_dzH2O]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb102_fixO]
    addps xmm8, [rsp + nb102_fiyO]
    addps xmm9, [rsp + nb102_fizO]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb102_fixH1]
    addps xmm11, [rsp + nb102_fiyH1]
    addps xmm12, [rsp + nb102_fizH1]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb102_fixH2]
    addps xmm14, [rsp + nb102_fiyH2]
    addps xmm15, [rsp + nb102_fizH2]

    movaps [rsp + nb102_fixO], xmm7
    movaps [rsp + nb102_fiyO], xmm8
    movaps [rsp + nb102_fizO], xmm9
    movaps [rsp + nb102_fixH1], xmm10
    movaps [rsp + nb102_fiyH1], xmm11
    movaps [rsp + nb102_fizH1], xmm12
    movaps [rsp + nb102_fixH2], xmm13
    movaps [rsp + nb102_fiyH2], xmm14
    movaps [rsp + nb102_fizH2], xmm15
    
    ;# xmm3 = fOx , xmm4 = fOy
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fOzc fOzd
    
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

	;# move j H1 coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, [rsi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, [rsi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, [rsi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, [rsi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, [rsi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, [rsi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH1a jxH1c jyH1a jyH1c        
    unpckhps xmm4, xmm1  ;# jxH1b jxH1d jyH1b jyH1d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0 = H1x
    ;# xmm1 = H1y
    ;# xmm2 = H1z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb102_ixO]
    subps xmm1, [rsp + nb102_iyO]
    subps xmm2, [rsp + nb102_izO]
    subps xmm3, [rsp + nb102_ixH1]
    subps xmm4, [rsp + nb102_iyH1]
    subps xmm5, [rsp + nb102_izH1]
    subps xmm6, [rsp + nb102_ixH2]
    subps xmm7, [rsp + nb102_iyH2]
    subps xmm8, [rsp + nb102_izH2]
    
	movaps [rsp + nb102_dxOH1], xmm0
	movaps [rsp + nb102_dyOH1], xmm1
	movaps [rsp + nb102_dzOH1], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb102_dxH1H1], xmm3
	movaps [rsp + nb102_dyH1H1], xmm4
	movaps [rsp + nb102_dzH1H1], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb102_dxH2H1], xmm6
	movaps [rsp + nb102_dyH2H1], xmm7
	movaps [rsp + nb102_dzH2H1], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH1 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb102_three]
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

	movaps  xmm0, [rsp + nb102_half]
	mulps   xmm9, xmm0  ;# rinvOH1
	mulps   xmm10, xmm0 ;# rinvH1H1
    mulps   xmm11, xmm0 ;# rinvH2H1
	
	;# H1 interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb102_qqOH] 
    mulps  xmm1, [rsp + nb102_qqHH] 
    mulps  xmm2, [rsp + nb102_qqHH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb102_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb102_vctot], xmm0
    
	;# move j H1 forces to local temp variables 
    movlps xmm0, [rdi + rax*4 + 12] ;# jxH1a jyH1a  -   -
    movlps xmm1, [rdi + rcx*4 + 12] ;# jxH1c jyH1c  -   -
    movhps xmm0, [rdi + rbx*4 + 12] ;# jxH1a jyH1a jxH1b jyH1b 
    movhps xmm1, [rdi + rdx*4 + 12] ;# jxH1c jyH1c jxH1d jyH1d 

    movss  xmm2, [rdi + rax*4 + 20] ;# jzH1a  -  -  -
    movss  xmm3, [rdi + rcx*4 + 20] ;# jzH1c  -  -  -
    movhps xmm2, [rdi + rbx*4 + 20] ;# jzH1a  -  jzH1b  -
    movhps xmm3, [rdi + rdx*4 + 20] ;# jzH1c  -  jzH1d -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jzH1a jzH1b jzH1c jzH1d

    ;# xmm0: jxH1a jyH1a jxH1b jyH1b 
    ;# xmm1: jxH1c jyH1c jxH1d jyH1d
    ;# xmm2: jzH1a jzH1b jzH1c jzH1d

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb102_dxOH1]
	mulps xmm8, [rsp + nb102_dyOH1]
	mulps xmm9, [rsp + nb102_dzOH1]
	mulps xmm10, [rsp + nb102_dxH1H1]
	mulps xmm11, [rsp + nb102_dyH1H1]
	mulps xmm12, [rsp + nb102_dzH1H1]
	mulps xmm13, [rsp + nb102_dxH2H1]
	mulps xmm14, [rsp + nb102_dyH2H1]
	mulps xmm15, [rsp + nb102_dzH2H1]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb102_fixO]
    addps xmm8, [rsp + nb102_fiyO]
    addps xmm9, [rsp + nb102_fizO]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb102_fixH1]
    addps xmm11, [rsp + nb102_fiyH1]
    addps xmm12, [rsp + nb102_fizH1]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb102_fixH2]
    addps xmm14, [rsp + nb102_fiyH2]
    addps xmm15, [rsp + nb102_fizH2]
    
    movaps [rsp + nb102_fixO], xmm7
    movaps [rsp + nb102_fiyO], xmm8
    movaps [rsp + nb102_fizO], xmm9
    movaps [rsp + nb102_fixH1], xmm10
    movaps [rsp + nb102_fiyH1], xmm11
    movaps [rsp + nb102_fizH1], xmm12
    movaps [rsp + nb102_fixH2], xmm13
    movaps [rsp + nb102_fiyH2], xmm14
    movaps [rsp + nb102_fizH2], xmm15

    ;# xmm3 = fH1x , xmm4 = fH1y
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fH1zc fH1zd
    
    movlps [rdi + rax*4 + 12], xmm0
    movhps [rdi + rbx*4 + 12], xmm0
    movlps [rdi + rcx*4 + 12], xmm1
    movhps [rdi + rdx*4 + 12], xmm1
    movss  [rdi + rax*4 + 20], xmm2
    movss  [rdi + rcx*4 + 20], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 20], xmm2
    movss  [rdi + rdx*4 + 20], xmm3

	;# move j H2 coordinates to local temp variables 
    movlps xmm0, [rsi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, [rsi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, [rsi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, [rsi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2, [rsi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3, [rsi + rcx*4 + 32] ;# jzH2c  -  -  -
    movss  xmm5, [rsi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm6, [rsi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm2, xmm5 ;# jzH2a  -  jzH2b  -
    movlhps xmm3, xmm6 ;# jzH2c  -  jzH2d -
    
    movaps xmm4, xmm0
    unpcklps xmm0, xmm1  ;# jxH2a jxH2c jyH2a jyH2c        
    unpckhps xmm4, xmm1  ;# jxH2b jxH2d jyH2b jyH2d
    movaps xmm1, xmm0
    unpcklps xmm0, xmm4 ;# x
    unpckhps xmm1, xmm4 ;# y

    shufps   xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0 = H2x
    ;# xmm1 = H2y
    ;# xmm2 = H2z
        
    movaps xmm3, xmm0
    movaps xmm4, xmm1
    movaps xmm5, xmm2
    movaps xmm6, xmm0
    movaps xmm7, xmm1
    movaps xmm8, xmm2
    
    subps xmm0, [rsp + nb102_ixO]
    subps xmm1, [rsp + nb102_iyO]
    subps xmm2, [rsp + nb102_izO]
    subps xmm3, [rsp + nb102_ixH1]
    subps xmm4, [rsp + nb102_iyH1]
    subps xmm5, [rsp + nb102_izH1]
    subps xmm6, [rsp + nb102_ixH2]
    subps xmm7, [rsp + nb102_iyH2]
    subps xmm8, [rsp + nb102_izH2]
    
	movaps [rsp + nb102_dxOH2], xmm0
	movaps [rsp + nb102_dyOH2], xmm1
	movaps [rsp + nb102_dzOH2], xmm2
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2
	movaps [rsp + nb102_dxH1H2], xmm3
	movaps [rsp + nb102_dyH1H2], xmm4
	movaps [rsp + nb102_dzH1H2], xmm5
	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	movaps [rsp + nb102_dxH2H2], xmm6
	movaps [rsp + nb102_dyH2H2], xmm7
	movaps [rsp + nb102_dzH2H2], xmm8
	mulps  xmm6, xmm6
	mulps  xmm7, xmm7
	mulps  xmm8, xmm8
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
    addps  xmm6, xmm7
    addps  xmm6, xmm8

	;# start doing invsqrt for jH2 atoms
	rsqrtps xmm1, xmm0
	rsqrtps xmm4, xmm3
    rsqrtps xmm7, xmm6
	
	movaps  xmm2, xmm1
	movaps  xmm5, xmm4
    movaps  xmm8, xmm7
    
    
	mulps   xmm1, xmm1 ;# lu*lu
	mulps   xmm4, xmm4 ;# lu*lu
    mulps   xmm7, xmm7 ;# lu*lu
		
	movaps  xmm9, [rsp + nb102_three]
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

	movaps  xmm0, [rsp + nb102_half]
	mulps   xmm9, xmm0  ;# rinvOH2 
	mulps   xmm10, xmm0 ;# rinvH1H2
    mulps   xmm11, xmm0 ;# rinvH2H2
	
	;# H2 interactions 
    movaps xmm0, xmm9
    movaps xmm1, xmm10
    movaps xmm2, xmm11
    mulps  xmm9, xmm9
    mulps  xmm10, xmm10
    mulps  xmm11, xmm11
    mulps  xmm0, [rsp + nb102_qqOH] 
    mulps  xmm1, [rsp + nb102_qqHH] 
    mulps  xmm2, [rsp + nb102_qqHH] 
    mulps  xmm9, xmm0
    mulps  xmm10, xmm1
    mulps  xmm11, xmm2
    
    addps xmm0, [rsp + nb102_vctot] 
    addps xmm1, xmm2
    addps xmm0, xmm1
    movaps [rsp + nb102_vctot], xmm0
    
	;# move j H2 forces to local temp variables 
    movlps xmm0, [rdi + rax*4 + 24] ;# jxH2a jyH2a  -   -
    movlps xmm1, [rdi + rcx*4 + 24] ;# jxH2c jyH2c  -   -
    movhps xmm0, [rdi + rbx*4 + 24] ;# jxH2a jyH2a jxH2b jyH2b 
    movhps xmm1, [rdi + rdx*4 + 24] ;# jxH2c jyH2c jxH2d jyH2d 

    movss  xmm2,  [rdi + rax*4 + 32] ;# jzH2a  -  -  -
    movss  xmm3,  [rdi + rcx*4 + 32] ;# jzH2c  -  -  -
    movss  xmm7, [rdi + rbx*4 + 32] ;# jzH2b  -  -  -
    movss  xmm8, [rdi + rdx*4 + 32] ;# jzH2d  -  -  -
    movlhps xmm2, xmm7 ;# jzH2a  -  jzH2b  -
    movlhps xmm3, xmm8 ;# jzH2c  -  jzH2d -
    
    shufps xmm2, xmm3,  136  ;# 10001000 => jzH2a jzH2b jzH2c jzH2d

    ;# xmm0: jxH2a jyH2a jxH2b jyH2b 
    ;# xmm1: jxH2c jyH2c jxH2d jyH2d
    ;# xmm2: jzH2a jzH2b jzH2c jzH2d

    movaps xmm7, xmm9
    movaps xmm8, xmm9
    movaps xmm13, xmm11
    movaps xmm14, xmm11
    movaps xmm15, xmm11
    movaps xmm11, xmm10
    movaps xmm12, xmm10

	mulps xmm7, [rsp + nb102_dxOH2]
	mulps xmm8, [rsp + nb102_dyOH2]
	mulps xmm9, [rsp + nb102_dzOH2]
	mulps xmm10, [rsp + nb102_dxH1H2]
	mulps xmm11, [rsp + nb102_dyH1H2]
	mulps xmm12, [rsp + nb102_dzH1H2]
	mulps xmm13, [rsp + nb102_dxH2H2]
	mulps xmm14, [rsp + nb102_dyH2H2]
	mulps xmm15, [rsp + nb102_dzH2H2]

    movaps xmm3, xmm7
    movaps xmm4, xmm8
    addps xmm2, xmm9
    addps xmm7, [rsp + nb102_fixO]
    addps xmm8, [rsp + nb102_fiyO]
    addps xmm9, [rsp + nb102_fizO]

    addps xmm3, xmm10
    addps xmm4, xmm11
    addps xmm2, xmm12
    addps xmm10, [rsp + nb102_fixH1]
    addps xmm11, [rsp + nb102_fiyH1]
    addps xmm12, [rsp + nb102_fizH1]

    addps xmm3, xmm13
    addps xmm4, xmm14
    addps xmm2, xmm15
    addps xmm13, [rsp + nb102_fixH2]
    addps xmm14, [rsp + nb102_fiyH2]
    addps xmm15, [rsp + nb102_fizH2]
    
    movaps [rsp + nb102_fixO], xmm7
    movaps [rsp + nb102_fiyO], xmm8
    movaps [rsp + nb102_fizO], xmm9
    movaps [rsp + nb102_fixH1], xmm10
    movaps [rsp + nb102_fiyH1], xmm11
    movaps [rsp + nb102_fizH1], xmm12
    movaps [rsp + nb102_fixH2], xmm13
    movaps [rsp + nb102_fiyH2], xmm14
    movaps [rsp + nb102_fizH2], xmm15

    ;# xmm3 = fH2x , xmm4 = fH2y
    movaps xmm5, xmm3
    unpcklps xmm3, xmm4
    unpckhps xmm5, xmm4
    
    addps xmm0, xmm3
    addps xmm1, xmm5

    movhlps  xmm3, xmm2 ;# fH2zc fH2zd
    
    movlps [rdi + rax*4 + 24], xmm0
    movhps [rdi + rbx*4 + 24], xmm0
    movlps [rdi + rcx*4 + 24], xmm1
    movhps [rdi + rdx*4 + 24], xmm1
    movss  [rdi + rax*4 + 32], xmm2
    movss  [rdi + rcx*4 + 32], xmm3
    shufps xmm2, xmm2, 1
    shufps xmm3, xmm3, 1
    movss  [rdi + rbx*4 + 32], xmm2
    movss  [rdi + rdx*4 + 32], xmm3

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb102_innerk],  4
	jl    .nb102_single_check
	jmp   .nb102_unroll_loop
    
.nb102_single_check:
	add dword ptr [rsp + nb102_innerk],  4
	jnz   .nb102_single_loop
	jmp   .nb102_updateouterdata
.nb102_single_loop:
	mov   rdx, [rsp + nb102_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	               ;# index is 32 bits 
	add qword ptr [rsp + nb102_innerjjnr],  4	

	mov rsi, [rbp + nb102_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm0, xmm0
	xorps xmm1, xmm1
	xorps xmm2, xmm2
	
	movss xmm0, [rsi + rax*4]		;# jxO  -  -  -
	movss xmm1, [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm2, [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm5, [rsi + rax*4 + 32]		;# xmm5 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm5			;# xmm7 = jzH1 jzH2   -    -

	movlhps xmm0, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm1, xmm6, 228 ;# 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm2, xmm7, 68  ;# 01000100	;# xmm5 = jzO   0   jzH1 jzH2

	;# store all j coordinates in jO  

	movaps [rsp + nb102_jxO], xmm0
	movaps [rsp + nb102_jyO], xmm1
	movaps [rsp + nb102_jzO], xmm2
	subps  xmm0, [rsp + nb102_ixO]
	subps  xmm1, [rsp + nb102_iyO]
	subps  xmm2, [rsp + nb102_izO]
	movaps [rsp + nb102_dxOO], xmm0
	movaps [rsp + nb102_dyOO], xmm1
	movaps [rsp + nb102_dzOO], xmm2
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb102_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb102_half] ;# rinv iO - j water 

	xorps   xmm1, xmm1
	movaps  xmm0, xmm3
	xorps   xmm4, xmm4
	mulps   xmm0, xmm0	;# xmm0=rinvsq 
	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [rsp + nb102_qqOO]

	movhps  xmm4, [rsp + nb102_qqOH]

	mulps   xmm3, xmm4	;# xmm3=vcoul 
	mulps   xmm0, xmm3	;# total fscal 
	addps   xmm3, [rsp + nb102_vctot]
	movaps  [rsp + nb102_vctot], xmm3	

	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [rsp + nb102_dxOO]
	mulps   xmm1, [rsp + nb102_dyOO]
	mulps   xmm2, [rsp + nb102_dzOO]
	;# initial update for j forces 
	xorps   xmm3, xmm3
	xorps   xmm4, xmm4
	xorps   xmm5, xmm5
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb102_fjxO], xmm3
	movaps  [rsp + nb102_fjyO], xmm4
	movaps  [rsp + nb102_fjzO], xmm5
	addps   xmm0, [rsp + nb102_fixO]
	addps   xmm1, [rsp + nb102_fiyO]
	addps   xmm2, [rsp + nb102_fizO]
	movaps  [rsp + nb102_fixO], xmm0
	movaps  [rsp + nb102_fiyO], xmm1
	movaps  [rsp + nb102_fizO], xmm2
	
	;# done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
    movaps  xmm0, [rsp + nb102_jxO]
    movaps  xmm1, [rsp + nb102_jyO]
    movaps  xmm2, [rsp + nb102_jzO]
    movaps  xmm3, xmm0
    movaps  xmm4, xmm1
    movaps  xmm5, xmm2
	subps  xmm0, [rsp + nb102_ixH1]
	subps  xmm1, [rsp + nb102_iyH1]
	subps  xmm2, [rsp + nb102_izH1]	
	subps  xmm3, [rsp + nb102_ixH2] 
	subps  xmm4, [rsp + nb102_iyH2] 
	subps  xmm5, [rsp + nb102_izH2] 
        
	movaps [rsp + nb102_dxH1O], xmm0
	movaps [rsp + nb102_dyH1O], xmm1
	movaps [rsp + nb102_dzH1O], xmm2
	movaps [rsp + nb102_dxH2O], xmm3
	movaps [rsp + nb102_dyH2O], xmm4
	movaps [rsp + nb102_dzH2O], xmm5
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 

	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb102_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb102_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb102_half] ;# rinv H2 - j water  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	;# do coulomb interaction 
	movaps  xmm0, xmm3
	movss   xmm6, [rsp + nb102_qqOH]
	movaps  xmm4, xmm7
	movhps  xmm6, [rsp + nb102_qqHH]
	mulps   xmm0, xmm0	;# rinvsq 
	mulps   xmm4, xmm4	;# rinvsq 
	mulps   xmm3, xmm6	;# vcoul 
	mulps   xmm7, xmm6	;# vcoul 
	movaps  xmm2, xmm3
	addps   xmm2, xmm7	;# total vcoul 
	mulps   xmm0, xmm3	;# fscal 
	
	addps   xmm2, [rsp + nb102_vctot]
	mulps   xmm7, xmm4	;# fscal 
	movaps  [rsp + nb102_vctot], xmm2
	movaps  xmm1, xmm0
	movaps  xmm2, xmm0
	mulps   xmm0, [rsp + nb102_dxH1O]
	mulps   xmm1, [rsp + nb102_dyH1O]
	mulps   xmm2, [rsp + nb102_dzH1O]
	;# update forces H1 - j water 
	movaps  xmm3, [rsp + nb102_fjxO]
	movaps  xmm4, [rsp + nb102_fjyO]
	movaps  xmm5, [rsp + nb102_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	movaps  [rsp + nb102_fjxO], xmm3
	movaps  [rsp + nb102_fjyO], xmm4
	movaps  [rsp + nb102_fjzO], xmm5
	addps   xmm0, [rsp + nb102_fixH1]
	addps   xmm1, [rsp + nb102_fiyH1]
	addps   xmm2, [rsp + nb102_fizH1]
	movaps  [rsp + nb102_fixH1], xmm0
	movaps  [rsp + nb102_fiyH1], xmm1
	movaps  [rsp + nb102_fizH1], xmm2
	;# do forces H2 - j water 
	movaps xmm0, xmm7
	movaps xmm1, xmm7
	movaps xmm2, xmm7
	mulps   xmm0, [rsp + nb102_dxH2O]
	mulps   xmm1, [rsp + nb102_dyH2O]
	mulps   xmm2, [rsp + nb102_dzH2O]
	movaps  xmm3, [rsp + nb102_fjxO]
	movaps  xmm4, [rsp + nb102_fjyO]
	movaps  xmm5, [rsp + nb102_fjzO]
	addps   xmm3, xmm0
	addps   xmm4, xmm1
	addps   xmm5, xmm2
	mov     rsi, [rbp + nb102_faction]
	movaps  [rsp + nb102_fjxO], xmm3
	movaps  [rsp + nb102_fjyO], xmm4
	movaps  [rsp + nb102_fjzO], xmm5
	addps   xmm0, [rsp + nb102_fixH2]
	addps   xmm1, [rsp + nb102_fiyH2]
	addps   xmm2, [rsp + nb102_fizH2]
	movaps  [rsp + nb102_fixH2], xmm0
	movaps  [rsp + nb102_fiyH2], xmm1
	movaps  [rsp + nb102_fizH2], xmm2

	;# update j water forces from local variables 
	movlps  xmm0, [rsi + rax*4]
	movlps  xmm1, [rsi + rax*4 + 12]
	movhps  xmm1, [rsi + rax*4 + 24]
	movaps  xmm3, [rsp + nb102_fjxO]
	movaps  xmm4, [rsp + nb102_fjyO]
	movaps  xmm5, [rsp + nb102_fjzO]
	movaps  xmm6, xmm5
	movaps  xmm7, xmm5
	shufps  xmm6, xmm6, 2 ;# 00000010
	shufps  xmm7, xmm7, 3 ;# 00000011
	addss   xmm5, [rsi + rax*4 + 8]
	addss   xmm6, [rsi + rax*4 + 20]
	addss   xmm7, [rsi + rax*4 + 32]
	movss   [rsi + rax*4 + 8], xmm5
	movss   [rsi + rax*4 + 20], xmm6
	movss   [rsi + rax*4 + 32], xmm7
	movaps   xmm5, xmm3
	unpcklps xmm3, xmm4
	unpckhps xmm5, xmm4
	addps    xmm0, xmm3
	addps    xmm1, xmm5
	movlps  [rsi + rax*4], xmm0 
	movlps  [rsi + rax*4 + 12], xmm1 
	movhps  [rsi + rax*4 + 24], xmm1 
	
	dec   dword ptr [rsp + nb102_innerk]
	jz    .nb102_updateouterdata
	jmp   .nb102_single_loop
.nb102_updateouterdata:
	mov   ecx, [rsp + nb102_ii3]
	mov   rdi, [rbp + nb102_faction]
	mov   rsi, [rbp + nb102_fshift]
	mov   edx, [rsp + nb102_is3]

	;# accumulate Oi forces in xmm0, xmm1, xmm2 
	movaps xmm0, [rsp + nb102_fixO]
	movaps xmm1, [rsp + nb102_fiyO] 
	movaps xmm2, [rsp + nb102_fizO]

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
	movaps xmm0, [rsp + nb102_fixH1]
	movaps xmm1, [rsp + nb102_fiyH1]
	movaps xmm2, [rsp + nb102_fizH1]

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
	movaps xmm0, [rsp + nb102_fixH2]
	movaps xmm1, [rsp + nb102_fiyH2]
	movaps xmm2, [rsp + nb102_fizH2]

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
	mov esi, [rsp + nb102_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb102_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb102_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb102_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 	

        ;# finish if last 
        mov ecx, [rsp + nb102_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb102_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb102_n], esi
        jmp .nb102_outer
.nb102_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb102_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb102_end
        ;# non-zero, do one more workunit
        jmp   .nb102_threadloop
.nb102_end:

	mov eax, [rsp + nb102_nouter]
	mov ebx, [rsp + nb102_ninner]
	mov rcx, [rbp + nb102_outeriter]
	mov rdx, [rbp + nb102_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 1504
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


	
.globl nb_kernel102nf_x86_64_sse
.globl _nb_kernel102nf_x86_64_sse
nb_kernel102nf_x86_64_sse:	
_nb_kernel102nf_x86_64_sse:	
.equiv          nb102nf_fshift,         16
.equiv          nb102nf_gid,            24
.equiv          nb102nf_pos,            32
.equiv          nb102nf_faction,        40
.equiv          nb102nf_charge,         48
.equiv          nb102nf_p_facel,        56
.equiv          nb102nf_argkrf,         64
.equiv          nb102nf_argcrf,         72
.equiv          nb102nf_Vc,             80
.equiv          nb102nf_type,           88
.equiv          nb102nf_p_ntype,        96
.equiv          nb102nf_vdwparam,       104
.equiv          nb102nf_Vvdw,           112
.equiv          nb102nf_p_tabscale,     120
.equiv          nb102nf_VFtab,          128
.equiv          nb102nf_invsqrta,       136
.equiv          nb102nf_dvda,           144
.equiv          nb102nf_p_gbtabscale,   152
.equiv          nb102nf_GBtab,          160
.equiv          nb102nf_p_nthreads,     168
.equiv          nb102nf_count,          176
.equiv          nb102nf_mtx,            184
.equiv          nb102nf_outeriter,      192
.equiv          nb102nf_inneriter,      200
.equiv          nb102nf_work,           208
	;# stack offsets for local variables  
	;# bottom of stack is cache-aligned for sse use 	
.equiv          nb102nf_ixO,            0
.equiv          nb102nf_iyO,            16
.equiv          nb102nf_izO,            32
.equiv          nb102nf_ixH1,           48
.equiv          nb102nf_iyH1,           64
.equiv          nb102nf_izH1,           80
.equiv          nb102nf_ixH2,           96
.equiv          nb102nf_iyH2,           112
.equiv          nb102nf_izH2,           128
.equiv          nb102nf_jxO,            144
.equiv          nb102nf_jyO,            160
.equiv          nb102nf_jzO,            176
.equiv          nb102nf_jxH1,           192
.equiv          nb102nf_jyH1,           208
.equiv          nb102nf_jzH1,           224
.equiv          nb102nf_jxH2,           240
.equiv          nb102nf_jyH2,           256
.equiv          nb102nf_jzH2,           272
.equiv          nb102nf_qqOO,           288
.equiv          nb102nf_qqOH,           304
.equiv          nb102nf_qqHH,           320
.equiv          nb102nf_vctot,          336
.equiv          nb102nf_half,           352
.equiv          nb102nf_three,          368
.equiv          nb102nf_rsqOO,          384
.equiv          nb102nf_rsqOH1,         400
.equiv          nb102nf_rsqOH2,         416
.equiv          nb102nf_rsqH1O,         432
.equiv          nb102nf_rsqH1H1,        448
.equiv          nb102nf_rsqH1H2,        464
.equiv          nb102nf_rsqH2O,         480
.equiv          nb102nf_rsqH2H1,        496
.equiv          nb102nf_rsqH2H2,        512
.equiv          nb102nf_rinvOO,         528
.equiv          nb102nf_rinvOH1,        544
.equiv          nb102nf_rinvOH2,        560
.equiv          nb102nf_rinvH1O,        576
.equiv          nb102nf_rinvH1H1,       592
.equiv          nb102nf_rinvH1H2,       608
.equiv          nb102nf_rinvH2O,        624
.equiv          nb102nf_rinvH2H1,       640
.equiv          nb102nf_rinvH2H2,       656
.equiv          nb102nf_is3,            672
.equiv          nb102nf_ii3,            676
.equiv          nb102nf_nri,            680
.equiv          nb102nf_iinr,           688
.equiv          nb102nf_jindex,         696
.equiv          nb102nf_jjnr,           704
.equiv          nb102nf_shift,          712
.equiv          nb102nf_shiftvec,       720
.equiv          nb102nf_facel,          728
.equiv          nb102nf_innerjjnr,      736
.equiv          nb102nf_innerk,         744
.equiv          nb102nf_n,              748
.equiv          nb102nf_nn1,            752
.equiv          nb102nf_nouter,         756
.equiv          nb102nf_ninner,         760

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
	sub rsp, 768	
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
	mov [rsp + nb102nf_nouter], eax
	mov [rsp + nb102nf_ninner], eax

	mov edi, [rdi]
	mov [rsp + nb102nf_nri], edi
	mov [rsp + nb102nf_iinr], rsi
	mov [rsp + nb102nf_jindex], rdx
	mov [rsp + nb102nf_jjnr], rcx
	mov [rsp + nb102nf_shift], r8
	mov [rsp + nb102nf_shiftvec], r9
	mov rsi, [rbp + nb102nf_p_facel]
	movss xmm0, [rsi]
	movss [rsp + nb102nf_facel], xmm0


	;# assume we have at least one i particle - start directly 
	mov   rcx, [rsp + nb102nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx]	    ;# ebx =ii 

	mov   rdx, [rbp + nb102nf_charge]
	movss xmm3, [rdx + rbx*4]	
	movss xmm4, xmm3	
	movss xmm5, [rdx + rbx*4 + 4]	
	mov rsi, [rbp + nb102nf_p_facel]
	movss xmm0, [rsi]
	movss xmm6, [rsp + nb102nf_facel]
	mulss  xmm3, xmm3
	mulss  xmm4, xmm5
	mulss  xmm5, xmm5
	mulss  xmm3, xmm6
	mulss  xmm4, xmm6
	mulss  xmm5, xmm6
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb102nf_qqOO], xmm3
	movaps [rsp + nb102nf_qqOH], xmm4
	movaps [rsp + nb102nf_qqHH], xmm5


	;# create constant floating-point factors on stack
	mov eax, 0x3f000000     ;# half in IEEE (hex)
	mov [rsp + nb102nf_half], eax
	movss xmm1, [rsp + nb102nf_half]
	shufps xmm1, xmm1, 0    ;# splat to all elements
	movaps xmm2, xmm1       
	addps  xmm2, xmm2	;# one
	movaps xmm3, xmm2
	addps  xmm2, xmm2	;# two
	addps  xmm3, xmm2	;# three
	movaps [rsp + nb102nf_half],  xmm1
	movaps [rsp + nb102nf_three],  xmm3

.nb102nf_threadloop:
        mov   rsi, [rbp + nb102nf_count]          ;# pointer to sync counter
        mov   eax, [rsi]
.nb102nf_spinlock:
        mov   ebx, eax                          ;# ebx=*count=nn0
        add   ebx, 1                           ;# ebx=nn1=nn0+10
        lock
        cmpxchg [rsi], ebx                      ;# write nn1 to *counter,
                                                ;# if it hasnt changed.
                                                ;# or reread *counter to eax.
        pause                                   ;# -> better p4 performance
        jnz .nb102nf_spinlock

        ;# if(nn1>nri) nn1=nri
        mov ecx, [rsp + nb102nf_nri]
        mov edx, ecx
        sub ecx, ebx
        cmovle ebx, edx                         ;# if(nn1>nri) nn1=nri
        ;# Cleared the spinlock if we got here.
        ;# eax contains nn0, ebx contains nn1.
        mov [rsp + nb102nf_n], eax
        mov [rsp + nb102nf_nn1], ebx
        sub ebx, eax                            ;# calc number of outer lists
	mov esi, eax				;# copy n to esi
        jg  .nb102nf_outerstart
        jmp .nb102nf_end
	
.nb102nf_outerstart:
	;# ebx contains number of outer iterations
	add ebx, [rsp + nb102nf_nouter]
	mov [rsp + nb102nf_nouter], ebx

.nb102nf_outer:
	mov   rax, [rsp + nb102nf_shift]      	;# rax = pointer into shift[] 
	mov   ebx, [rax+rsi*4]			;# ebx=shift[n] 
	
	lea   rbx, [rbx + rbx*2]    ;# rbx=3*is 
	mov   [rsp + nb102nf_is3],ebx    	;# store is3 

	mov   rax, [rsp + nb102nf_shiftvec]   ;# rax = base of shiftvec[] 

	movss xmm0, [rax + rbx*4]
	movss xmm1, [rax + rbx*4 + 4]
	movss xmm2, [rax + rbx*4 + 8] 

	mov   rcx, [rsp + nb102nf_iinr]       ;# rcx = pointer into iinr[] 	
	mov   ebx, [rcx + rsi*4]	    ;# ebx =ii 

	lea   rbx, [rbx + rbx*2]	;# rbx = 3*ii=ii3 
	mov   rax, [rbp + nb102nf_pos]    ;# rax = base of pos[]  
	mov   [rsp + nb102nf_ii3], ebx	
	
	movaps xmm3, xmm0
	movaps xmm4, xmm1
	movaps xmm5, xmm2
	addss xmm3, [rax + rbx*4]
	addss xmm4, [rax + rbx*4 + 4]
	addss xmm5, [rax + rbx*4 + 8]		
	shufps xmm3, xmm3, 0
	shufps xmm4, xmm4, 0
	shufps xmm5, xmm5, 0
	movaps [rsp + nb102nf_ixO], xmm3
	movaps [rsp + nb102nf_iyO], xmm4
	movaps [rsp + nb102nf_izO], xmm5

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
	movaps [rsp + nb102nf_ixH1], xmm0
	movaps [rsp + nb102nf_iyH1], xmm1
	movaps [rsp + nb102nf_izH1], xmm2
	movaps [rsp + nb102nf_ixH2], xmm3
	movaps [rsp + nb102nf_iyH2], xmm4
	movaps [rsp + nb102nf_izH2], xmm5

	;# clear vctot 
	xorps xmm4, xmm4
	movaps [rsp + nb102nf_vctot], xmm4
	
	mov   rax, [rsp + nb102nf_jindex]
	mov   ecx, [rax + rsi*4]	     	;# jindex[n] 
	mov   edx, [rax + rsi*4 + 4]	     	;# jindex[n+1] 
	sub   edx, ecx               		;# number of innerloop atoms 

	mov   rsi, [rbp + nb102nf_pos]
	mov   rax, [rsp + nb102nf_jjnr]
	shl   ecx, 2
	add   rax, rcx
	mov   [rsp + nb102nf_innerjjnr], rax     ;# pointer to jjnr[nj0] 
	mov   ecx, edx
	sub   edx,  4
	add   ecx, [rsp + nb102nf_ninner]
	mov   [rsp + nb102nf_ninner], ecx
	add   edx, 0
	mov   [rsp + nb102nf_innerk], edx    ;# number of innerloop atoms 
	jge   .nb102nf_unroll_loop
	jmp   .nb102nf_single_check
.nb102nf_unroll_loop:	
	;# quad-unroll innerloop here 
	mov   rdx, [rsp + nb102nf_innerjjnr]     ;# pointer to jjnr[k] 

	mov   eax, [rdx]	
	mov   ebx, [rdx + 4] 
	mov   ecx, [rdx + 8]
	mov   edx, [rdx + 12]         ;# eax-edx=jnr1-4 
	
	add qword ptr [rsp + nb102nf_innerjjnr],  16 ;# advance pointer (unrolled 4) 

	mov rsi, [rbp + nb102nf_pos]       ;# base of pos[] 

	lea   rax, [rax + rax*2]     ;# replace jnr with j3 
	lea   rbx, [rbx + rbx*2]	
	lea   rcx, [rcx + rcx*2]     ;# replace jnr with j3 
	lea   rdx, [rdx + rdx*2]	
	
	;# move j coordinates to local temp variables 
	movlps xmm2, [rsi + rax*4]
	movlps xmm3, [rsi + rax*4 + 12]
	movlps xmm4, [rsi + rax*4 + 24]

	movlps xmm5, [rsi + rbx*4]
	movlps xmm6, [rsi + rbx*4 + 12]
	movlps xmm7, [rsi + rbx*4 + 24]

	movhps xmm2, [rsi + rcx*4]
	movhps xmm3, [rsi + rcx*4 + 12]
	movhps xmm4, [rsi + rcx*4 + 24]

	movhps xmm5, [rsi + rdx*4]
	movhps xmm6, [rsi + rdx*4 + 12]
	movhps xmm7, [rsi + rdx*4 + 24]

	;# current state: 	
	;# xmm2= jxOa  jyOa  jxOc  jyOc 
	;# xmm3= jxH1a jyH1a jxH1c jyH1c 
	;# xmm4= jxH2a jyH2a jxH2c jyH2c 
	;# xmm5= jxOb  jyOb  jxOd  jyOd 
	;# xmm6= jxH1b jyH1b jxH1d jyH1d 
	;# xmm7= jxH2b jyH2b jxH2d jyH2d 
	
	movaps xmm0, xmm2
	movaps xmm1, xmm3
	unpcklps xmm0, xmm5	;# xmm0= jxOa  jxOb  jyOa  jyOb 
	unpcklps xmm1, xmm6	;# xmm1= jxH1a jxH1b jyH1a jyH1b 
	unpckhps xmm2, xmm5	;# xmm2= jxOc  jxOd  jyOc  jyOd 
	unpckhps xmm3, xmm6	;# xmm3= jxH1c jxH1d jyH1c jyH1d  
	movaps xmm5, xmm4
	movaps   xmm6, xmm0
	unpcklps xmm4, xmm7	;# xmm4= jxH2a jxH2b jyH2a jyH2b 		
	unpckhps xmm5, xmm7	;# xmm5= jxH2c jxH2d jyH2c jyH2d 
	movaps   xmm7, xmm1
	movlhps  xmm0, xmm2	;# xmm0= jxOa  jxOb  jxOc  jxOd  
	movaps [rsp + nb102nf_jxO], xmm0
	movhlps  xmm2, xmm6	;# xmm2= jyOa  jyOb  jyOc  jyOd 
	movaps [rsp + nb102nf_jyO], xmm2
	movlhps  xmm1, xmm3
	movaps [rsp + nb102nf_jxH1], xmm1
	movhlps  xmm3, xmm7
	movaps   xmm6, xmm4
	movaps [rsp + nb102nf_jyH1], xmm3
	movlhps  xmm4, xmm5
	movaps [rsp + nb102nf_jxH2], xmm4
	movhlps  xmm5, xmm6
	movaps [rsp + nb102nf_jyH2], xmm5

	movss  xmm0, [rsi + rax*4 + 8]
	movss  xmm1, [rsi + rax*4 + 20]
	movss  xmm2, [rsi + rax*4 + 32]

	movss  xmm3, [rsi + rcx*4 + 8]
	movss  xmm4, [rsi + rcx*4 + 20]
	movss  xmm5, [rsi + rcx*4 + 32]

	movhps xmm0, [rsi + rbx*4 + 4]
	movhps xmm1, [rsi + rbx*4 + 16]
	movhps xmm2, [rsi + rbx*4 + 28]
	
	movhps xmm3, [rsi + rdx*4 + 4]
	movhps xmm4, [rsi + rdx*4 + 16]
	movhps xmm5, [rsi + rdx*4 + 28]
	
	shufps xmm0, xmm3, 204  ;# 11001100
	shufps xmm1, xmm4, 204  ;# 11001100
	shufps xmm2, xmm5, 204  ;# 11001100
	movaps [rsp + nb102nf_jzO],  xmm0
	movaps [rsp + nb102nf_jzH1],  xmm1
	movaps [rsp + nb102nf_jzH2],  xmm2

	movaps xmm0, [rsp + nb102nf_ixO]
	movaps xmm1, [rsp + nb102nf_iyO]
	movaps xmm2, [rsp + nb102nf_izO]
	movaps xmm3, [rsp + nb102nf_ixO]
	movaps xmm4, [rsp + nb102nf_iyO]
	movaps xmm5, [rsp + nb102nf_izO]
	subps  xmm0, [rsp + nb102nf_jxO]
	subps  xmm1, [rsp + nb102nf_jyO]
	subps  xmm2, [rsp + nb102nf_jzO]
	subps  xmm3, [rsp + nb102nf_jxH1]
	subps  xmm4, [rsp + nb102nf_jyH1]
	subps  xmm5, [rsp + nb102nf_jzH1]
	
	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb102nf_rsqOO], xmm0
	movaps [rsp + nb102nf_rsqOH1], xmm3

	movaps xmm0, [rsp + nb102nf_ixO]
	movaps xmm1, [rsp + nb102nf_iyO]
	movaps xmm2, [rsp + nb102nf_izO]
	movaps xmm3, [rsp + nb102nf_ixH1]
	movaps xmm4, [rsp + nb102nf_iyH1]
	movaps xmm5, [rsp + nb102nf_izH1]
	subps  xmm0, [rsp + nb102nf_jxH2]
	subps  xmm1, [rsp + nb102nf_jyH2]
	subps  xmm2, [rsp + nb102nf_jzH2]
	subps  xmm3, [rsp + nb102nf_jxO]
	subps  xmm4, [rsp + nb102nf_jyO]
	subps  xmm5, [rsp + nb102nf_jzO]

	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb102nf_rsqOH2], xmm0
	movaps [rsp + nb102nf_rsqH1O], xmm3

	movaps xmm0, [rsp + nb102nf_ixH1]
	movaps xmm1, [rsp + nb102nf_iyH1]
	movaps xmm2, [rsp + nb102nf_izH1]
	movaps xmm3, [rsp + nb102nf_ixH1]
	movaps xmm4, [rsp + nb102nf_iyH1]
	movaps xmm5, [rsp + nb102nf_izH1]
	subps  xmm0, [rsp + nb102nf_jxH1]
	subps  xmm1, [rsp + nb102nf_jyH1]
	subps  xmm2, [rsp + nb102nf_jzH1]
	subps  xmm3, [rsp + nb102nf_jxH2]
	subps  xmm4, [rsp + nb102nf_jyH2]
	subps  xmm5, [rsp + nb102nf_jzH2]

	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm3, xmm4
	addps  xmm3, xmm5
	movaps [rsp + nb102nf_rsqH1H1], xmm0
	movaps [rsp + nb102nf_rsqH1H2], xmm3

	movaps xmm0, [rsp + nb102nf_ixH2]
	movaps xmm1, [rsp + nb102nf_iyH2]
	movaps xmm2, [rsp + nb102nf_izH2]
	movaps xmm3, [rsp + nb102nf_ixH2]
	movaps xmm4, [rsp + nb102nf_iyH2]
	movaps xmm5, [rsp + nb102nf_izH2]
	subps  xmm0, [rsp + nb102nf_jxO]
	subps  xmm1, [rsp + nb102nf_jyO]
	subps  xmm2, [rsp + nb102nf_jzO]
	subps  xmm3, [rsp + nb102nf_jxH1]
	subps  xmm4, [rsp + nb102nf_jyH1]
	subps  xmm5, [rsp + nb102nf_jzH1]

	mulps  xmm0, xmm0
	mulps  xmm1, xmm1
	mulps  xmm2, xmm2

	mulps  xmm3, xmm3
	mulps  xmm4, xmm4
	mulps  xmm5, xmm5
	addps  xmm0, xmm1
	addps  xmm0, xmm2
	addps  xmm4, xmm3
	addps  xmm4, xmm5
	movaps [rsp + nb102nf_rsqH2O], xmm0
	movaps [rsp + nb102nf_rsqH2H1], xmm4

	movaps xmm0, [rsp + nb102nf_ixH2]
	movaps xmm1, [rsp + nb102nf_iyH2]
	movaps xmm2, [rsp + nb102nf_izH2]
	subps  xmm0, [rsp + nb102nf_jxH2]
	subps  xmm1, [rsp + nb102nf_jyH2]
	subps  xmm2, [rsp + nb102nf_jzH2]

	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2
	movaps [rsp + nb102nf_rsqH2H2], xmm0
	
	;# start doing invsqrt use rsq values in xmm0, xmm4 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb102nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb102nf_half] ;# rinvH2H2 
	mulps   xmm7, [rsp + nb102nf_half] ;# rinvH2H1 
	movaps  [rsp + nb102nf_rinvH2H2], xmm3
	movaps  [rsp + nb102nf_rinvH2H1], xmm7
	
	rsqrtps xmm1, [rsp + nb102nf_rsqOO]
	rsqrtps xmm5, [rsp + nb102nf_rsqOH1]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb102nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb102nf_rsqOO]
	mulps   xmm5, [rsp + nb102nf_rsqOH1]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb102nf_half] 
	mulps   xmm7, [rsp + nb102nf_half]
	movaps  [rsp + nb102nf_rinvOO], xmm3
	movaps  [rsp + nb102nf_rinvOH1], xmm7
	
	rsqrtps xmm1, [rsp + nb102nf_rsqOH2]
	rsqrtps xmm5, [rsp + nb102nf_rsqH1O]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb102nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb102nf_rsqOH2]
	mulps   xmm5, [rsp + nb102nf_rsqH1O]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb102nf_half] 
	mulps   xmm7, [rsp + nb102nf_half]
	movaps  [rsp + nb102nf_rinvOH2], xmm3
	movaps  [rsp + nb102nf_rinvH1O], xmm7
	
	rsqrtps xmm1, [rsp + nb102nf_rsqH1H1]
	rsqrtps xmm5, [rsp + nb102nf_rsqH1H2]
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb102nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, [rsp + nb102nf_rsqH1H1]
	mulps   xmm5, [rsp + nb102nf_rsqH1H2]
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb102nf_half] 
	mulps   xmm7, [rsp + nb102nf_half]
	movaps  [rsp + nb102nf_rinvH1H1], xmm3
	movaps  [rsp + nb102nf_rinvH1H2], xmm7
	
	rsqrtps xmm1, [rsp + nb102nf_rsqH2O]
	movaps  xmm2, xmm1
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb102nf_three]
	mulps   xmm1, [rsp + nb102nf_rsqH2O]
	subps   xmm3, xmm1
	mulps   xmm3, xmm2
	mulps   xmm3, [rsp + nb102nf_half] 
	movaps  [rsp + nb102nf_rinvH2O], xmm3

	;# sum OO pot in xmm0, OH in xmm1 HH in xmm2 
	movaps xmm0, [rsp + nb102nf_rinvOO]
	movaps xmm1, [rsp + nb102nf_rinvOH1]
	movaps xmm2, [rsp + nb102nf_rinvH1H1]
	addps  xmm1, [rsp + nb102nf_rinvOH2]
	addps  xmm2, [rsp + nb102nf_rinvH1H2]
	addps  xmm1, [rsp + nb102nf_rinvH1O]
	addps  xmm2, [rsp + nb102nf_rinvH2H1]
	addps  xmm1, [rsp + nb102nf_rinvH2O]
	addps  xmm2, [rsp + nb102nf_rinvH2H2]
	
	mulps  xmm0, [rsp + nb102nf_qqOO]
	mulps  xmm1, [rsp + nb102nf_qqOH]
	mulps  xmm2, [rsp + nb102nf_qqHH]
	addps  xmm0, [rsp + nb102nf_vctot]
	addps  xmm1, xmm2
	addps  xmm0, xmm1
	movaps 	[rsp + nb102nf_vctot], xmm0

	;# should we do one more iteration? 
	sub dword ptr [rsp + nb102nf_innerk],  4
	jl    .nb102nf_single_check
	jmp   .nb102nf_unroll_loop
.nb102nf_single_check:
	add dword ptr [rsp + nb102nf_innerk],  4
	jnz   .nb102nf_single_loop
	jmp   .nb102nf_updateouterdata
.nb102nf_single_loop:
	mov   rdx, [rsp + nb102nf_innerjjnr]     ;# pointer to jjnr[k] 
	mov   eax, [rdx]	
	add qword ptr [rsp + nb102nf_innerjjnr],  4	

	mov rsi, [rbp + nb102nf_pos]
	lea   rax, [rax + rax*2]  

	;# fetch j coordinates 
	xorps xmm3, xmm3
	xorps xmm4, xmm4
	xorps xmm5, xmm5
	
	movss xmm3, [rsi + rax*4]		;# jxO  -  -  -
	movss xmm4, [rsi + rax*4 + 4]		;# jyO  -  -  -
	movss xmm5, [rsi + rax*4 + 8]		;# jzO  -  -  -  

	movlps xmm6, [rsi + rax*4 + 12]		;# xmm6 = jxH1 jyH1   -    -
	movss  xmm7, [rsi + rax*4 + 20]		;# xmm7 = jzH1   -    -    - 
	movhps xmm6, [rsi + rax*4 + 24]		;# xmm6 = jxH1 jyH1 jxH2 jyH2
	movss  xmm2, [rsi + rax*4 + 32]		;# xmm2 = jzH2   -    -    -
	
	;# have all coords, time for some shuffling.

	shufps xmm6, xmm6, 216 ;# 11011000	;# xmm6 = jxH1 jxH2 jyH1 jyH2 
	unpcklps xmm7, xmm2			;# xmm7 = jzH1 jzH2   -    -
	movaps  xmm0, [rsp + nb102nf_ixO]     
	movaps  xmm1, [rsp + nb102nf_iyO]
	movaps  xmm2, [rsp + nb102nf_izO]	
	movlhps xmm3, xmm6			;# xmm3 = jxO   0   jxH1 jxH2 
	shufps  xmm4, xmm6, 228 ;# 11100100	;# xmm4 = jyO   0   jyH1 jyH2 
	shufps  xmm5, xmm7, 68  ;# 01000100	;# xmm5 = jzO   0   jzH1 jzH2

	;# store all j coordinates in jO  
	movaps [rsp + nb102nf_jxO], xmm3
	movaps [rsp + nb102nf_jyO], xmm4
	movaps [rsp + nb102nf_jzO], xmm5
	subps  xmm0, xmm3
	subps  xmm1, xmm4
	subps  xmm2, xmm5
	
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	addps xmm0, xmm1
	addps xmm0, xmm2	;# have rsq in xmm0 
	
	;# do invsqrt 
	rsqrtps xmm1, xmm0
	movaps  xmm2, xmm1	
	mulps   xmm1, xmm1
	movaps  xmm3, [rsp + nb102nf_three]
	mulps   xmm1, xmm0
	subps   xmm3, xmm1
	mulps   xmm3, xmm2							
	mulps   xmm3, [rsp + nb102nf_half] ;# rinv iO - j water 

	xorps   xmm1, xmm1

	xorps   xmm4, xmm4

	;# fetch charges to xmm4 (temporary) 
	movss   xmm4, [rsp + nb102nf_qqOO]

	movhps  xmm4, [rsp + nb102nf_qqOH]

	mulps   xmm3, xmm4	;# xmm3=vcoul 

	addps   xmm3, [rsp + nb102nf_vctot]
	movaps  [rsp + nb102nf_vctot], xmm3	

	;# done with i O Now do i H1 & H2 simultaneously: 
	movaps  xmm0, [rsp + nb102nf_ixH1]
	movaps  xmm1, [rsp + nb102nf_iyH1]
	movaps  xmm2, [rsp + nb102nf_izH1]	
	movaps  xmm3, [rsp + nb102nf_ixH2] 
	movaps  xmm4, [rsp + nb102nf_iyH2] 
	movaps  xmm5, [rsp + nb102nf_izH2] 
	subps   xmm0, [rsp + nb102nf_jxO]
	subps   xmm1, [rsp + nb102nf_jyO]
	subps   xmm2, [rsp + nb102nf_jzO]
	subps   xmm3, [rsp + nb102nf_jxO]
	subps   xmm4, [rsp + nb102nf_jyO]
	subps   xmm5, [rsp + nb102nf_jzO]
	mulps xmm0, xmm0
	mulps xmm1, xmm1
	mulps xmm2, xmm2
	mulps xmm3, xmm3
	mulps xmm4, xmm4
	mulps xmm5, xmm5
	addps xmm0, xmm1
	addps xmm4, xmm3
	addps xmm0, xmm2	;# have rsqH1 in xmm0 
	addps xmm4, xmm5	;# have rsqH2 in xmm4 

	;# do invsqrt 
	rsqrtps xmm1, xmm0
	rsqrtps xmm5, xmm4
	movaps  xmm2, xmm1
	movaps  xmm6, xmm5
	mulps   xmm1, xmm1
	mulps   xmm5, xmm5
	movaps  xmm3, [rsp + nb102nf_three]
	movaps  xmm7, xmm3
	mulps   xmm1, xmm0
	mulps   xmm5, xmm4
	subps   xmm3, xmm1
	subps   xmm7, xmm5
	mulps   xmm3, xmm2
	mulps   xmm7, xmm6
	mulps   xmm3, [rsp + nb102nf_half] ;# rinv H1 - j water 
	mulps   xmm7, [rsp + nb102nf_half] ;# rinv H2 - j water  

	;# assemble charges in xmm6 
	xorps   xmm6, xmm6
	;# do coulomb interaction 
	movaps  xmm0, xmm3
	movss   xmm6, [rsp + nb102nf_qqOH]
	movaps  xmm4, xmm7
	movhps  xmm6, [rsp + nb102nf_qqHH]
	mulps   xmm3, xmm6	;# vcoul 
	mulps   xmm7, xmm6	;# vcoul 
	addps   xmm3, xmm7	;# total vcoul 
	addps   xmm3, [rsp + nb102nf_vctot]
	movaps  [rsp + nb102nf_vctot], xmm3
	
	dec   dword ptr [rsp + nb102nf_innerk]
	jz    .nb102nf_updateouterdata
	jmp   .nb102nf_single_loop
.nb102nf_updateouterdata:
	;# get n from stack
	mov esi, [rsp + nb102nf_n]
        ;# get group index for i particle 
        mov   rdx, [rbp + nb102nf_gid]      	;# base of gid[]
        mov   edx, [rdx + rsi*4]		;# ggid=gid[n]

	;# accumulate total potential energy and update it 
	movaps xmm7, [rsp + nb102nf_vctot]
	;# accumulate 
	movhlps xmm6, xmm7
	addps  xmm7, xmm6	;# pos 0-1 in xmm7 have the sum now 
	movaps xmm6, xmm7
	shufps xmm6, xmm6, 1
	addss  xmm7, xmm6		

	;# add earlier value from mem 
	mov   rax, [rbp + nb102nf_Vc]
	addss xmm7, [rax + rdx*4] 
	;# move back to mem 
	movss [rax + rdx*4], xmm7 	

        ;# finish if last 
        mov ecx, [rsp + nb102nf_nn1]
	;# esi already loaded with n
	inc esi
        sub ecx, esi
        jz .nb102nf_outerend

        ;# not last, iterate outer loop once more!  
        mov [rsp + nb102nf_n], esi
        jmp .nb102nf_outer
.nb102nf_outerend:
        ;# check if more outer neighborlists remain
        mov   ecx, [rsp + nb102nf_nri]
	;# esi already loaded with n above
        sub   ecx, esi
        jz .nb102nf_end
        ;# non-zero, do one more workunit
        jmp   .nb102nf_threadloop
.nb102nf_end:
	

	mov eax, [rsp + nb102nf_nouter]
	mov ebx, [rsp + nb102nf_ninner]
	mov rcx, [rbp + nb102nf_outeriter]
	mov rdx, [rbp + nb102nf_inneriter]
	mov [rcx], eax
	mov [rdx], ebx

	add rsp, 768
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
