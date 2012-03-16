##
##
## Gromacs 4.0                         Copyright (c) 1991-2003 
## David van der Spoel, Erik Lindahl
##
## This program is free software; you can redistribute it and/or
## modify it under the terms of the GNU General Public License
## as published by the Free Software Foundation; either version 2
## of the License, or (at your option) any later version.
##
## To help us fund GROMACS development, we humbly ask that you cite
## the research papers on the package. Check out http://www.gromacs.org
## 
## And Hey:
## Gnomes, ROck Monsters And Chili Sauce
##





.globl nb_kernel103_x86_64_sse
.globl _nb_kernel103_x86_64_sse
nb_kernel103_x86_64_sse:        
_nb_kernel103_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb103_fshift, 16
.set nb103_gid, 24
.set nb103_pos, 32
.set nb103_faction, 40
.set nb103_charge, 48
.set nb103_p_facel, 56
.set nb103_argkrf, 64
.set nb103_argcrf, 72
.set nb103_Vc, 80
.set nb103_type, 88
.set nb103_p_ntype, 96
.set nb103_vdwparam, 104
.set nb103_Vvdw, 112
.set nb103_p_tabscale, 120
.set nb103_VFtab, 128
.set nb103_invsqrta, 136
.set nb103_dvda, 144
.set nb103_p_gbtabscale, 152
.set nb103_GBtab, 160
.set nb103_p_nthreads, 168
.set nb103_count, 176
.set nb103_mtx, 184
.set nb103_outeriter, 192
.set nb103_inneriter, 200
.set nb103_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb103_ixH1, 0
.set nb103_iyH1, 16
.set nb103_izH1, 32
.set nb103_ixH2, 48
.set nb103_iyH2, 64
.set nb103_izH2, 80
.set nb103_ixM, 96
.set nb103_iyM, 112
.set nb103_izM, 128
.set nb103_iqH, 144
.set nb103_iqM, 160
.set nb103_dxH1, 176
.set nb103_dyH1, 192
.set nb103_dzH1, 208
.set nb103_dxH2, 224
.set nb103_dyH2, 240
.set nb103_dzH2, 256
.set nb103_dxM, 272
.set nb103_dyM, 288
.set nb103_dzM, 304
.set nb103_qqH, 320
.set nb103_qqM, 336
.set nb103_vctot, 352
.set nb103_fixH1, 368
.set nb103_fiyH1, 384
.set nb103_fizH1, 400
.set nb103_fixH2, 416
.set nb103_fiyH2, 432
.set nb103_fizH2, 448
.set nb103_fixM, 464
.set nb103_fiyM, 480
.set nb103_fizM, 496
.set nb103_fjx, 512
.set nb103_fjy, 528
.set nb103_fjz, 544
.set nb103_half, 560
.set nb103_three, 576
.set nb103_is3, 592
.set nb103_ii3, 596
.set nb103_nri, 600
.set nb103_innerjjnr, 608
.set nb103_iinr, 616
.set nb103_jindex, 624
.set nb103_jjnr, 632
.set nb103_shift, 640
.set nb103_shiftvec, 648
.set nb103_facel, 656
.set nb103_innerk, 664
.set nb103_n, 668
.set nb103_nn1, 672
.set nb103_nouter, 676
.set nb103_ninner, 680

        push %rbp
        movq %rsp,%rbp
        push %rbx

        push %r12
        push %r13
        push %r14
        push %r15

        emms
        subq $696,%rsp          # # local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb103_nouter(%rsp)
        movl %eax,nb103_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb103_nri(%rsp)
        movq %rsi,nb103_iinr(%rsp)
        movq %rdx,nb103_jindex(%rsp)
        movq %rcx,nb103_jjnr(%rsp)
        movq %r8,nb103_shift(%rsp)
        movq %r9,nb103_shiftvec(%rsp)
        movq nb103_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb103_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb103_half(%rsp)
        movss nb103_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb103_half(%rsp)
        movaps %xmm3,nb103_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb103_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb103_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss 12(%rdx,%rbx,4),%xmm4
        movq nb103_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb103_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb103_iqH(%rsp)
        movaps %xmm4,nb103_iqM(%rsp)

_nb_kernel103_x86_64_sse.nb103_threadloop: 
        movq  nb103_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel103_x86_64_sse.nb103_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel103_x86_64_sse.nb103_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb103_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb103_n(%rsp)
        movl %ebx,nb103_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi

        jg  _nb_kernel103_x86_64_sse.nb103_outerstart
        jmp _nb_kernel103_x86_64_sse.nb103_end

_nb_kernel103_x86_64_sse.nb103_outerstart: 
        ## ebx contains number of outer iterations
        addl nb103_nouter(%rsp),%ebx
        movl %ebx,nb103_nouter(%rsp)

_nb_kernel103_x86_64_sse.nb103_outer: 
        movq  nb103_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb103_is3(%rsp)      ## store is3 

        movq  nb103_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb103_iinr(%rsp),%rcx             ## rcx = pointer into iinr[] 
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb103_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb103_ii3(%rsp)

        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb103_ixH1(%rsp)
        movaps %xmm4,nb103_iyH1(%rsp)
        movaps %xmm5,nb103_izH1(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%rax,%rbx,4),%xmm0
        addss 28(%rax,%rbx,4),%xmm1
        addss 32(%rax,%rbx,4),%xmm2
        addss 36(%rax,%rbx,4),%xmm3
        addss 40(%rax,%rbx,4),%xmm4
        addss 44(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb103_ixH2(%rsp)
        movaps %xmm1,nb103_iyH2(%rsp)
        movaps %xmm2,nb103_izH2(%rsp)
        movaps %xmm3,nb103_ixM(%rsp)
        movaps %xmm4,nb103_iyM(%rsp)
        movaps %xmm5,nb103_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb103_vctot(%rsp)
        movaps %xmm4,nb103_fixH1(%rsp)
        movaps %xmm4,nb103_fiyH1(%rsp)
        movaps %xmm4,nb103_fizH1(%rsp)
        movaps %xmm4,nb103_fixH2(%rsp)
        movaps %xmm4,nb103_fiyH2(%rsp)
        movaps %xmm4,nb103_fizH2(%rsp)
        movaps %xmm4,nb103_fixM(%rsp)
        movaps %xmm4,nb103_fiyM(%rsp)
        movaps %xmm4,nb103_fizM(%rsp)

        movq  nb103_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb103_pos(%rbp),%rsi
        movq  nb103_faction(%rbp),%rdi
        movq  nb103_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb103_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb103_ninner(%rsp),%ecx
        movl  %ecx,nb103_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb103_innerk(%rsp)   ## number of innerloop atoms 

        jge   _nb_kernel103_x86_64_sse.nb103_unroll_loop
        jmp   _nb_kernel103_x86_64_sse.nb103_odd_inner
_nb_kernel103_x86_64_sse.nb103_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb103_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%r8d
        movl  4(%rdx),%r9d
        movl  8(%rdx),%r10d
        movl  12(%rdx),%r11d            ## eax-edx=jnr1-4 

        addq $16,nb103_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb103_charge(%rbp),%rsi    ## base of charge[] 

        movss (%rsi,%r8,4),%xmm12
        movss (%rsi,%r10,4),%xmm13
        movss (%rsi,%r9,4),%xmm14
        movss (%rsi,%r11,4),%xmm15

        shufps $0,%xmm14,%xmm12
        shufps $0,%xmm15,%xmm13
        shufps $136,%xmm13,%xmm12 ## 10001000 ;# all charges in xmm3  
        movaps %xmm12,%xmm13         ## and in xmm4 
        mulps  nb103_iqH(%rsp),%xmm12
        mulps  nb103_iqM(%rsp),%xmm13

        movaps  %xmm12,nb103_qqH(%rsp)
        movaps  %xmm13,nb103_qqM(%rsp)

        movq nb103_pos(%rbp),%rsi       ## base of pos[] 

        lea  (%r8,%r8,2),%rax     ## j3 
        lea  (%r9,%r9,2),%rbx
        lea  (%r10,%r10,2),%rcx
        lea  (%r11,%r11,2),%rdx

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8


    subps nb103_ixH1(%rsp),%xmm0
    subps nb103_iyH1(%rsp),%xmm1
    subps nb103_izH1(%rsp),%xmm2
    subps nb103_ixH2(%rsp),%xmm3
    subps nb103_iyH2(%rsp),%xmm4
    subps nb103_izH2(%rsp),%xmm5
    subps nb103_ixM(%rsp),%xmm6
    subps nb103_iyM(%rsp),%xmm7
    subps nb103_izM(%rsp),%xmm8

        movaps %xmm0,nb103_dxH1(%rsp)
        movaps %xmm1,nb103_dyH1(%rsp)
        movaps %xmm2,nb103_dzH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb103_dxH2(%rsp)
        movaps %xmm4,nb103_dyH2(%rsp)
        movaps %xmm5,nb103_dzH2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb103_dxM(%rsp)
        movaps %xmm7,nb103_dyM(%rsp)
        movaps %xmm8,nb103_dzM(%rsp)
        mulps  %xmm6,%xmm6
        mulps  %xmm7,%xmm7
        mulps  %xmm8,%xmm8
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
    addps  %xmm7,%xmm6
    addps  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm3,%xmm4
    rsqrtps %xmm6,%xmm7

        movaps  %xmm1,%xmm2
        movaps  %xmm4,%xmm5
    movaps  %xmm7,%xmm8

        mulps   %xmm1,%xmm1 ## lu*lu
        mulps   %xmm4,%xmm4 ## lu*lu
    mulps   %xmm7,%xmm7 ## lu*lu

        movaps  nb103_three(%rsp),%xmm9
        movaps  %xmm9,%xmm10
    movaps  %xmm9,%xmm11

        mulps   %xmm0,%xmm1 ## rsq*lu*lu
        mulps   %xmm3,%xmm4 ## rsq*lu*lu 
    mulps   %xmm6,%xmm7 ## rsq*lu*lu

        subps   %xmm1,%xmm9
        subps   %xmm4,%xmm10
    subps   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulps   %xmm2,%xmm9
        mulps   %xmm5,%xmm10
    mulps   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movaps  nb103_half(%rsp),%xmm0
        mulps   %xmm0,%xmm9 ## rinvH1
        mulps   %xmm0,%xmm10 ## rinvH2
    mulps   %xmm0,%xmm11 ## rinvM

        ## interactions 
    movaps %xmm9,%xmm0
    movaps %xmm10,%xmm1
    movaps %xmm11,%xmm2
    mulps  %xmm9,%xmm9
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  nb103_qqH(%rsp),%xmm0
    mulps  nb103_qqH(%rsp),%xmm1
    mulps  nb103_qqM(%rsp),%xmm2
    mulps  %xmm0,%xmm9
    mulps  %xmm1,%xmm10
    mulps  %xmm2,%xmm11

    addps nb103_vctot(%rsp),%xmm0
    addps %xmm2,%xmm1
    addps %xmm1,%xmm0
    movaps %xmm0,nb103_vctot(%rsp)

        ## move j forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm0 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm1 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm0 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm1 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm2    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm3    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm4    ## jzb
    movss  8(%rdi,%rdx,4),%xmm5    ## jzd
    movlhps %xmm4,%xmm2 ## jza  -  jzb  -
    movlhps %xmm5,%xmm3 ## jzc  -  jzd -

    shufps $136,%xmm3,%xmm2 ## 10001000 => jza jzb jzc jzd

    ## xmm0: jxa jya jxb jyb 
    ## xmm1: jxc jyc jxd jyd
    ## xmm2: jza jzb jzc jzd

    movaps %xmm9,%xmm7
    movaps %xmm9,%xmm8
    movaps %xmm11,%xmm13
    movaps %xmm11,%xmm14
    movaps %xmm11,%xmm15
    movaps %xmm10,%xmm11
    movaps %xmm10,%xmm12

        mulps nb103_dxH1(%rsp),%xmm7
        mulps nb103_dyH1(%rsp),%xmm8
        mulps nb103_dzH1(%rsp),%xmm9
        mulps nb103_dxH2(%rsp),%xmm10
        mulps nb103_dyH2(%rsp),%xmm11
        mulps nb103_dzH2(%rsp),%xmm12
        mulps nb103_dxM(%rsp),%xmm13
        mulps nb103_dyM(%rsp),%xmm14
        mulps nb103_dzM(%rsp),%xmm15

    movaps %xmm7,%xmm3
    movaps %xmm8,%xmm4
    addps %xmm9,%xmm2
    addps nb103_fixH1(%rsp),%xmm7
    addps nb103_fiyH1(%rsp),%xmm8
    addps nb103_fizH1(%rsp),%xmm9

    addps %xmm10,%xmm3
    addps %xmm11,%xmm4
    addps %xmm12,%xmm2
    addps nb103_fixH2(%rsp),%xmm10
    addps nb103_fiyH2(%rsp),%xmm11
    addps nb103_fizH2(%rsp),%xmm12

    addps %xmm13,%xmm3
    addps %xmm14,%xmm4
    addps %xmm15,%xmm2
    addps nb103_fixM(%rsp),%xmm13
    addps nb103_fiyM(%rsp),%xmm14
    addps nb103_fizM(%rsp),%xmm15

    movaps %xmm7,nb103_fixH1(%rsp)
    movaps %xmm8,nb103_fiyH1(%rsp)
    movaps %xmm9,nb103_fizH1(%rsp)
    movaps %xmm10,nb103_fixH2(%rsp)
    movaps %xmm11,nb103_fiyH2(%rsp)
    movaps %xmm12,nb103_fizH2(%rsp)
    movaps %xmm13,nb103_fixM(%rsp)
    movaps %xmm14,nb103_fiyM(%rsp)
    movaps %xmm15,nb103_fizM(%rsp)

    ## xmm3 = fjx , xmm4 = fjy
    movaps %xmm3,%xmm5
    unpcklps %xmm4,%xmm3
    unpckhps %xmm4,%xmm5

    addps %xmm3,%xmm0
    addps %xmm5,%xmm1

    movhlps  %xmm2,%xmm3 ## fzc fzd

    movlps %xmm0,(%rdi,%rax,4)
    movhps %xmm0,(%rdi,%rbx,4)
    movlps %xmm1,(%rdi,%rcx,4)
    movhps %xmm1,(%rdi,%rdx,4)
    movss  %xmm2,8(%rdi,%rax,4)
    movss  %xmm3,8(%rdi,%rcx,4)
    shufps $1,%xmm2,%xmm2
    shufps $1,%xmm3,%xmm3
    movss  %xmm2,8(%rdi,%rbx,4)
    movss  %xmm3,8(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb103_innerk(%rsp)
        jl    _nb_kernel103_x86_64_sse.nb103_odd_inner
        jmp   _nb_kernel103_x86_64_sse.nb103_unroll_loop
_nb_kernel103_x86_64_sse.nb103_odd_inner: 
        addl $4,nb103_innerk(%rsp)
        jnz   _nb_kernel103_x86_64_sse.nb103_odd_loop
        jmp   _nb_kernel103_x86_64_sse.nb103_updateouterdata
_nb_kernel103_x86_64_sse.nb103_odd_loop: 
        movq  nb103_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb103_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb103_iqM(%rsp),%xmm4
        movq nb103_charge(%rbp),%rsi
        movhps nb103_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb103_qqM(%rsp)    ## use dummy qq for storage 

        movq nb103_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb103_ixM(%rsp),%xmm0
        movss nb103_iyM(%rsp),%xmm1
        movss nb103_izM(%rsp),%xmm2

        movlps nb103_ixH1(%rsp),%xmm6
        movlps nb103_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb103_iyH1(%rsp),%xmm6
        movlps nb103_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb103_izH1(%rsp),%xmm6
        movlps nb103_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use dummy dx for storage
        movaps %xmm3,nb103_dxM(%rsp)
        movaps %xmm4,nb103_dyM(%rsp)
        movaps %xmm5,nb103_dzM(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb103_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb103_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      
        movaps %xmm0,%xmm4
        mulps  %xmm4,%xmm4      ## xmm4=rinvsq 
        movaps nb103_qqM(%rsp),%xmm3

        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        mulps  %xmm3,%xmm4      ## xmm4=total fscal 
        addps  nb103_vctot(%rsp),%xmm3

        movaps nb103_dxM(%rsp),%xmm0
        movaps nb103_dyM(%rsp),%xmm1
        movaps nb103_dzM(%rsp),%xmm2

        movaps %xmm3,nb103_vctot(%rsp)

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb103_fixM(%rsp),%xmm3
        movss  nb103_fiyM(%rsp),%xmm4
        movss  nb103_fizM(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb103_fixM(%rsp)
        movss  %xmm4,nb103_fiyM(%rsp)
        movss  %xmm5,nb103_fizM(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb103_fixH1(%rsp),%xmm3
        addss  nb103_fiyH1(%rsp),%xmm4
        addss  nb103_fizH1(%rsp),%xmm5
        movss  %xmm3,nb103_fixH1(%rsp)
        movss  %xmm4,nb103_fiyH1(%rsp)
        movss  %xmm5,nb103_fizH1(%rsp)          ## updated the H1 force 

        movq nb103_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb103_fixH2(%rsp),%xmm3
        addss  nb103_fiyH2(%rsp),%xmm4
        addss  nb103_fizH2(%rsp),%xmm5
        movss  %xmm3,nb103_fixH2(%rsp)
        movss  %xmm4,nb103_fiyH2(%rsp)
        movss  %xmm5,nb103_fizH2(%rsp)          ## updated the H2 force 

        ## the fj's - start by accumulating the tx/ty/tz force in xmm0, xmm1 
        xorps  %xmm5,%xmm5
        movaps %xmm0,%xmm3
        movlps (%rdi,%rax,4),%xmm6
        movss  8(%rdi,%rax,4),%xmm7
        unpcklps %xmm1,%xmm3
        movlhps  %xmm5,%xmm3
        unpckhps %xmm1,%xmm0
        addps    %xmm3,%xmm0
        movhlps  %xmm0,%xmm3
        addps    %xmm3,%xmm0    ## x,y sum in xmm0 

        movhlps  %xmm2,%xmm1
        addss    %xmm1,%xmm2
        shufps  $1,%xmm1,%xmm1
        addss    %xmm1,%xmm2    ## z sum in xmm2 
        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl  nb103_innerk(%rsp)
        jz    _nb_kernel103_x86_64_sse.nb103_updateouterdata
        jmp   _nb_kernel103_x86_64_sse.nb103_odd_loop
_nb_kernel103_x86_64_sse.nb103_updateouterdata: 

        movl  nb103_ii3(%rsp),%ecx
        movq  nb103_faction(%rbp),%rdi
        movq  nb103_fshift(%rbp),%rsi
        movl  nb103_is3(%rsp),%edx

        ## accumulate H1 forces in xmm0, xmm1, xmm2 
        movaps nb103_fixH1(%rsp),%xmm0
        movaps nb103_fiyH1(%rsp),%xmm1
        movaps nb103_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  12(%rdi,%rcx,4),%xmm3
        movss  16(%rdi,%rcx,4),%xmm4
        movss  20(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,12(%rdi,%rcx,4)
        movss  %xmm4,16(%rdi,%rcx,4)
        movss  %xmm5,20(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## 00001000       

        ## accumulate H2 i forces in xmm0, xmm1, xmm2 
        movaps nb103_fixH2(%rsp),%xmm0
        movaps nb103_fiyH2(%rsp),%xmm1
        movaps nb103_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  24(%rdi,%rcx,4),%xmm3
        movss  28(%rdi,%rcx,4),%xmm4
        movss  32(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,24(%rdi,%rcx,4)
        movss  %xmm4,28(%rdi,%rcx,4)
        movss  %xmm5,32(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## accumulate M i forces in xmm0, xmm1, xmm2 
        movaps nb103_fixM(%rsp),%xmm0
        movaps nb103_fiyM(%rsp),%xmm1
        movaps nb103_fizM(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps  %xmm3,%xmm0
        addps  %xmm4,%xmm1
        addps  %xmm5,%xmm2 ## sum is in 1/2 in xmm0-xmm2 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        shufps $1,%xmm3,%xmm3
        shufps $1,%xmm4,%xmm4
        shufps $1,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2      ## xmm0-xmm2 has single force in pos0 

        ## increment i force 
        movss  36(%rdi,%rcx,4),%xmm3
        movss  40(%rdi,%rcx,4),%xmm4
        movss  44(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,36(%rdi,%rcx,4)
        movss  %xmm4,40(%rdi,%rcx,4)
        movss  %xmm5,44(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%rsi,%rdx,4),%xmm3
        movss  8(%rsi,%rdx,4),%xmm4

        subps  %xmm6,%xmm3
        subss  %xmm7,%xmm4
        movlps  %xmm3,(%rsi,%rdx,4)
        movss  %xmm4,8(%rsi,%rdx,4)

        ## get n from stack
        movl nb103_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb103_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb103_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb103_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb103_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx

        jz _nb_kernel103_x86_64_sse.nb103_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb103_n(%rsp)
        jmp _nb_kernel103_x86_64_sse.nb103_outer
_nb_kernel103_x86_64_sse.nb103_outerend: 
        ## check if more outer neighborlists remain
        movl  nb103_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel103_x86_64_sse.nb103_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel103_x86_64_sse.nb103_threadloop
_nb_kernel103_x86_64_sse.nb103_end: 

        movl nb103_nouter(%rsp),%eax
        movl nb103_ninner(%rsp),%ebx
        movq nb103_outeriter(%rbp),%rcx
        movq nb103_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $696,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




.globl nb_kernel103nf_x86_64_sse
.globl _nb_kernel103nf_x86_64_sse
nb_kernel103nf_x86_64_sse:      
_nb_kernel103nf_x86_64_sse:     
.set nb103nf_fshift, 16
.set nb103nf_gid, 24
.set nb103nf_pos, 32
.set nb103nf_faction, 40
.set nb103nf_charge, 48
.set nb103nf_p_facel, 56
.set nb103nf_argkrf, 64
.set nb103nf_argcrf, 72
.set nb103nf_Vc, 80
.set nb103nf_type, 88
.set nb103nf_p_ntype, 96
.set nb103nf_vdwparam, 104
.set nb103nf_Vvdw, 112
.set nb103nf_p_tabscale, 120
.set nb103nf_VFtab, 128
.set nb103nf_invsqrta, 136
.set nb103nf_dvda, 144
.set nb103nf_p_gbtabscale, 152
.set nb103nf_GBtab, 160
.set nb103nf_p_nthreads, 168
.set nb103nf_count, 176
.set nb103nf_mtx, 184
.set nb103nf_outeriter, 192
.set nb103nf_inneriter, 200
.set nb103nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb103nf_ixH1, 0
.set nb103nf_iyH1, 16
.set nb103nf_izH1, 32
.set nb103nf_ixH2, 48
.set nb103nf_iyH2, 64
.set nb103nf_izH2, 80
.set nb103nf_ixM, 96
.set nb103nf_iyM, 112
.set nb103nf_izM, 128
.set nb103nf_iqH, 144
.set nb103nf_iqM, 160
.set nb103nf_vctot, 176
.set nb103nf_half, 192
.set nb103nf_three, 208
.set nb103nf_qqH, 224
.set nb103nf_qqM, 240
.set nb103nf_is3, 256
.set nb103nf_ii3, 260
.set nb103nf_nri, 264
.set nb103nf_iinr, 272
.set nb103nf_jindex, 280
.set nb103nf_jjnr, 288
.set nb103nf_shift, 296
.set nb103nf_shiftvec, 304
.set nb103nf_facel, 312
.set nb103nf_innerjjnr, 320
.set nb103nf_innerk, 328
.set nb103nf_n, 332
.set nb103nf_nn1, 336
.set nb103nf_nouter, 340
.set nb103nf_ninner, 344

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms
        subq $360,%rsp

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb103nf_nouter(%rsp)
        movl %eax,nb103nf_ninner(%rsp)


        movl (%rdi),%edi
        movl %edi,nb103nf_nri(%rsp)
        movq %rsi,nb103nf_iinr(%rsp)
        movq %rdx,nb103nf_jindex(%rsp)
        movq %rcx,nb103nf_jjnr(%rsp)
        movq %r8,nb103nf_shift(%rsp)
        movq %r9,nb103nf_shiftvec(%rsp)
        movq nb103nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb103nf_facel(%rsp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb103nf_half(%rsp)
        movss nb103nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb103nf_half(%rsp)
        movaps %xmm3,nb103nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb103nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb103nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss 12(%rdx,%rbx,4),%xmm4
        movq nb103nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb103nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb103nf_iqH(%rsp)
        movaps %xmm4,nb103nf_iqM(%rsp)

_nb_kernel103nf_x86_64_sse.nb103nf_threadloop: 
        movq  nb103nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel103nf_x86_64_sse.nb103nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel103nf_x86_64_sse.nb103nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb103nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb103nf_n(%rsp)
        movl %ebx,nb103nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel103nf_x86_64_sse.nb103nf_outerstart
        jmp _nb_kernel103nf_x86_64_sse.nb103nf_end

_nb_kernel103nf_x86_64_sse.nb103nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb103nf_nouter(%rsp),%ebx
        movl %ebx,nb103nf_nouter(%rsp)

_nb_kernel103nf_x86_64_sse.nb103nf_outer: 
        movq  nb103nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb103nf_is3(%rsp)            ## store is3 

        movq  nb103nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb103nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[] 
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb103nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb103nf_ii3(%rsp)

        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb103nf_ixH1(%rsp)
        movaps %xmm4,nb103nf_iyH1(%rsp)
        movaps %xmm5,nb103nf_izH1(%rsp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%rax,%rbx,4),%xmm0
        addss 28(%rax,%rbx,4),%xmm1
        addss 32(%rax,%rbx,4),%xmm2
        addss 36(%rax,%rbx,4),%xmm3
        addss 40(%rax,%rbx,4),%xmm4
        addss 44(%rax,%rbx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb103nf_ixH2(%rsp)
        movaps %xmm1,nb103nf_iyH2(%rsp)
        movaps %xmm2,nb103nf_izH2(%rsp)
        movaps %xmm3,nb103nf_ixM(%rsp)
        movaps %xmm4,nb103nf_iyM(%rsp)
        movaps %xmm5,nb103nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb103nf_vctot(%rsp)

        movq  nb103nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb103nf_pos(%rbp),%rsi
        movq  nb103nf_faction(%rbp),%rdi
        movq  nb103nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb103nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb103nf_ninner(%rsp),%ecx
        movl  %ecx,nb103nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb103nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel103nf_x86_64_sse.nb103nf_unroll_loop
        jmp   _nb_kernel103nf_x86_64_sse.nb103nf_odd_inner
_nb_kernel103nf_x86_64_sse.nb103nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb103nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb103nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb103nf_charge(%rbp),%rsi  ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb103nf_iqH(%rsp),%xmm3
        mulps  nb103nf_iqM(%rsp),%xmm4

        movaps  %xmm3,nb103nf_qqH(%rsp)
        movaps  %xmm4,nb103nf_qqM(%rsp)

        movq nb103nf_pos(%rbp),%rsi     ## base of pos[] 

        lea  (%rax,%rax,2),%rax        ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx
        lea  (%rcx,%rcx,2),%rcx        ## replace jnr with j3 
        lea  (%rdx,%rdx,2),%rdx

        ## move four coordinates to xmm0-xmm2   
        movlps (%rsi,%rax,4),%xmm4
        movlps (%rsi,%rcx,4),%xmm5
        movss 8(%rsi,%rax,4),%xmm2
        movss 8(%rsi,%rcx,4),%xmm6

        movhps (%rsi,%rbx,4),%xmm4
        movhps (%rsi,%rdx,4),%xmm5

        movss 8(%rsi,%rbx,4),%xmm0
        movss 8(%rsi,%rdx,4),%xmm1

        shufps $0,%xmm0,%xmm2
        shufps $0,%xmm1,%xmm6

        movaps %xmm4,%xmm0
        movaps %xmm4,%xmm1

        shufps $136,%xmm6,%xmm2 ## 10001000

        shufps $136,%xmm5,%xmm0 ## 10001000
        shufps $221,%xmm5,%xmm1 ## 11011101             

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb103nf_ixH1(%rsp),%xmm4
        movaps nb103nf_iyH1(%rsp),%xmm5
        movaps nb103nf_izH1(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm4
        addps %xmm6,%xmm4
        movaps %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH2-izH2 to xmm4-xmm6 
        movaps nb103nf_ixH2(%rsp),%xmm4
        movaps nb103nf_iyH2(%rsp),%xmm5
        movaps nb103nf_izH2(%rsp),%xmm6

        ## calc dr 
        subps %xmm0,%xmm4
        subps %xmm1,%xmm5
        subps %xmm2,%xmm6

        ## square it 
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        mulps %xmm6,%xmm6
        addps %xmm5,%xmm6
        addps %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixM-izM to xmm3-xmm5  
        movaps nb103nf_ixM(%rsp),%xmm3
        movaps nb103nf_iyM(%rsp),%xmm4
        movaps nb103nf_izM(%rsp),%xmm5

        ## calc dr 
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## square it 
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm4,%xmm5
        addps %xmm3,%xmm5
        ## rsqM in xmm5, rsqH2 in xmm6, rsqH1 in xmm7 

        ## start with rsqH1 - seed in xmm2      
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb103nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb103nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvH1 in xmm7 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb103nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb103nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH2 in xmm6 
        ## rsqM - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb103nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb103nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvM in xmm5 

        ## do H1 interactions - xmm7=rinv
        mulps  nb103nf_qqH(%rsp),%xmm7          ## xmm7=vcoul 
        addps  nb103nf_vctot(%rsp),%xmm7
        movaps %xmm7,nb103nf_vctot(%rsp)

        ## H2 interactions - xmm6=rinv
        mulps  nb103nf_qqH(%rsp),%xmm6          ## xmm6=vcoul 
        addps  %xmm7,%xmm6
        movaps %xmm6,nb103nf_vctot(%rsp)

        ## M interactions  - xmm5=rinv
        mulps  nb103nf_qqM(%rsp),%xmm5          ## xmm5=vcoul 
        addps  %xmm6,%xmm5
        movaps %xmm5,nb103nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb103nf_innerk(%rsp)
        jl    _nb_kernel103nf_x86_64_sse.nb103nf_odd_inner
        jmp   _nb_kernel103nf_x86_64_sse.nb103nf_unroll_loop
_nb_kernel103nf_x86_64_sse.nb103nf_odd_inner: 
        addl $4,nb103nf_innerk(%rsp)
        jnz   _nb_kernel103nf_x86_64_sse.nb103nf_odd_loop
        jmp   _nb_kernel103nf_x86_64_sse.nb103nf_updateouterdata
_nb_kernel103nf_x86_64_sse.nb103nf_odd_loop: 
        movq  nb103nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb103nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb103nf_iqM(%rsp),%xmm4
        movq nb103nf_charge(%rbp),%rsi
        movhps nb103nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb103nf_qqM(%rsp)          ## use dummy qq for storage 

        movq nb103nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb103nf_ixM(%rsp),%xmm3
        movss nb103nf_iyM(%rsp),%xmm4
        movss nb103nf_izM(%rsp),%xmm5

        movlps nb103nf_ixH1(%rsp),%xmm6
        movlps nb103nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb103nf_iyH1(%rsp),%xmm6
        movlps nb103nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb103nf_izH1(%rsp),%xmm6
        movlps nb103nf_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb103nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb103nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000      - xmm0=rinv
        movaps nb103nf_qqM(%rsp),%xmm3
        mulps  %xmm0,%xmm3      ## xmm3=vcoul 
        addps  nb103nf_vctot(%rsp),%xmm3
        movaps %xmm3,nb103nf_vctot(%rsp)

        decl  nb103nf_innerk(%rsp)
        jz    _nb_kernel103nf_x86_64_sse.nb103nf_updateouterdata
        jmp   _nb_kernel103nf_x86_64_sse.nb103nf_odd_loop
_nb_kernel103nf_x86_64_sse.nb103nf_updateouterdata: 
        ## get n from stack
        movl nb103nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb103nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb103nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb103nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb103nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel103nf_x86_64_sse.nb103nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb103nf_n(%rsp)
        jmp _nb_kernel103nf_x86_64_sse.nb103nf_outer
_nb_kernel103nf_x86_64_sse.nb103nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb103nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel103nf_x86_64_sse.nb103nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel103nf_x86_64_sse.nb103nf_threadloop
_nb_kernel103nf_x86_64_sse.nb103nf_end: 

        movl nb103nf_nouter(%rsp),%eax
        movl nb103nf_ninner(%rsp),%ebx
        movq nb103nf_outeriter(%rbp),%rcx
        movq nb103nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $360,%rsp
        emms

        pop %rbx
        pop    %rbp
        ret


