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






.globl nb_kernel203_x86_64_sse
.globl _nb_kernel203_x86_64_sse
nb_kernel203_x86_64_sse:        
_nb_kernel203_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb203_fshift, 16
.set nb203_gid, 24
.set nb203_pos, 32
.set nb203_faction, 40
.set nb203_charge, 48
.set nb203_p_facel, 56
.set nb203_argkrf, 64
.set nb203_argcrf, 72
.set nb203_Vc, 80
.set nb203_type, 88
.set nb203_p_ntype, 96
.set nb203_vdwparam, 104
.set nb203_Vvdw, 112
.set nb203_p_tabscale, 120
.set nb203_VFtab, 128
.set nb203_invsqrta, 136
.set nb203_dvda, 144
.set nb203_p_gbtabscale, 152
.set nb203_GBtab, 160
.set nb203_p_nthreads, 168
.set nb203_count, 176
.set nb203_mtx, 184
.set nb203_outeriter, 192
.set nb203_inneriter, 200
.set nb203_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb203_ixH1, 0
.set nb203_iyH1, 16
.set nb203_izH1, 32
.set nb203_ixH2, 48
.set nb203_iyH2, 64
.set nb203_izH2, 80
.set nb203_ixM, 96
.set nb203_iyM, 112
.set nb203_izM, 128
.set nb203_iqH, 144
.set nb203_iqM, 160
.set nb203_dxH1, 176
.set nb203_dyH1, 192
.set nb203_dzH1, 208
.set nb203_dxH2, 224
.set nb203_dyH2, 240
.set nb203_dzH2, 256
.set nb203_dxM, 272
.set nb203_dyM, 288
.set nb203_dzM, 304
.set nb203_qqH, 320
.set nb203_qqM, 336
.set nb203_vctot, 352
.set nb203_fixH1, 384
.set nb203_fiyH1, 400
.set nb203_fizH1, 416
.set nb203_fixH2, 432
.set nb203_fiyH2, 448
.set nb203_fizH2, 464
.set nb203_fixM, 480
.set nb203_fiyM, 496
.set nb203_fizM, 512
.set nb203_fjx, 528
.set nb203_fjy, 544
.set nb203_fjz, 560
.set nb203_half, 576
.set nb203_three, 592
.set nb203_two, 608
.set nb203_krf, 624
.set nb203_crf, 640
.set nb203_krsqH1, 656
.set nb203_krsqH2, 672
.set nb203_krsqM, 688
.set nb203_is3, 704
.set nb203_ii3, 708
.set nb203_innerjjnr, 712
.set nb203_nri, 720
.set nb203_iinr, 728
.set nb203_jindex, 736
.set nb203_jjnr, 744
.set nb203_shift, 752
.set nb203_shiftvec, 760
.set nb203_facel, 768
.set nb203_innerk, 776
.set nb203_n, 780
.set nb203_nn1, 784
.set nb203_nouter, 788
.set nb203_ninner, 792

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $808,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb203_nouter(%rsp)
        movl %eax,nb203_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb203_nri(%rsp)
        movq %rsi,nb203_iinr(%rsp)
        movq %rdx,nb203_jindex(%rsp)
        movq %rcx,nb203_jjnr(%rsp)
        movq %r8,nb203_shift(%rsp)
        movq %r9,nb203_shiftvec(%rsp)
        movq nb203_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb203_facel(%rsp)


        movq nb203_argkrf(%rbp),%rsi
        movq nb203_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb203_krf(%rsp)
        movaps %xmm2,nb203_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb203_half(%rsp)
        movss nb203_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb203_half(%rsp)
        movaps %xmm2,nb203_two(%rsp)
        movaps %xmm3,nb203_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb203_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb203_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss 12(%rdx,%rbx,4),%xmm4
        movq nb203_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb203_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb203_iqH(%rsp)
        movaps %xmm4,nb203_iqM(%rsp)

_nb_kernel203_x86_64_sse.nb203_threadloop: 
        movq  nb203_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel203_x86_64_sse.nb203_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel203_x86_64_sse.nb203_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb203_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb203_n(%rsp)
        movl %ebx,nb203_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel203_x86_64_sse.nb203_outerstart
        jmp _nb_kernel203_x86_64_sse.nb203_end

_nb_kernel203_x86_64_sse.nb203_outerstart: 
        ## ebx contains number of outer iterations
        addl nb203_nouter(%rsp),%ebx
        movl %ebx,nb203_nouter(%rsp)

_nb_kernel203_x86_64_sse.nb203_outer: 
        movq  nb203_shift(%rsp),%rax            ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb203_is3(%rsp)      ## store is3 

        movq  nb203_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb203_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb203_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb203_ii3(%rsp)

        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb203_ixH1(%rsp)
        movaps %xmm4,nb203_iyH1(%rsp)
        movaps %xmm5,nb203_izH1(%rsp)

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
        movaps %xmm0,nb203_ixH2(%rsp)
        movaps %xmm1,nb203_iyH2(%rsp)
        movaps %xmm2,nb203_izH2(%rsp)
        movaps %xmm3,nb203_ixM(%rsp)
        movaps %xmm4,nb203_iyM(%rsp)
        movaps %xmm5,nb203_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb203_vctot(%rsp)
        movaps %xmm4,nb203_fixH1(%rsp)
        movaps %xmm4,nb203_fiyH1(%rsp)
        movaps %xmm4,nb203_fizH1(%rsp)
        movaps %xmm4,nb203_fixH2(%rsp)
        movaps %xmm4,nb203_fiyH2(%rsp)
        movaps %xmm4,nb203_fizH2(%rsp)
        movaps %xmm4,nb203_fixM(%rsp)
        movaps %xmm4,nb203_fiyM(%rsp)
        movaps %xmm4,nb203_fizM(%rsp)

        movq  nb203_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb203_pos(%rbp),%rsi
        movq  nb203_faction(%rbp),%rdi
        movq  nb203_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb203_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb203_ninner(%rsp),%ecx
        movl  %ecx,nb203_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb203_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel203_x86_64_sse.nb203_unroll_loop
        jmp   _nb_kernel203_x86_64_sse.nb203_odd_inner
_nb_kernel203_x86_64_sse.nb203_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb203_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb203_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb203_charge(%rbp),%rsi    ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb203_iqH(%rsp),%xmm3
        mulps  nb203_iqM(%rsp),%xmm4

        movaps  %xmm3,nb203_qqH(%rsp)
        movaps  %xmm4,nb203_qqM(%rsp)

        movq nb203_pos(%rbp),%rsi       ## base of pos[] 

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

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb203_ixH1(%rsp),%xmm0
    subps nb203_iyH1(%rsp),%xmm1
    subps nb203_izH1(%rsp),%xmm2
    subps nb203_ixH2(%rsp),%xmm3
    subps nb203_iyH2(%rsp),%xmm4
    subps nb203_izH2(%rsp),%xmm5
    subps nb203_ixM(%rsp),%xmm6
    subps nb203_iyM(%rsp),%xmm7
    subps nb203_izM(%rsp),%xmm8

        movaps %xmm0,nb203_dxH1(%rsp)
        movaps %xmm1,nb203_dyH1(%rsp)
        movaps %xmm2,nb203_dzH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb203_dxH2(%rsp)
        movaps %xmm4,nb203_dyH2(%rsp)
        movaps %xmm5,nb203_dzH2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb203_dxM(%rsp)
        movaps %xmm7,nb203_dyM(%rsp)
        movaps %xmm8,nb203_dzM(%rsp)
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

        movaps  nb203_three(%rsp),%xmm9
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

        movaps  nb203_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1 
        mulps   %xmm4,%xmm10 ## rinvH2
    mulps   %xmm4,%xmm11 ## rinvM

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb203_krf(%rsp),%xmm2
    mulps  %xmm9,%xmm9  ## rinvsq
    mulps  %xmm10,%xmm10
    mulps  %xmm11,%xmm11
    mulps  %xmm2,%xmm0 ## k*rsq
    mulps  %xmm2,%xmm3
    mulps  %xmm2,%xmm6
    movaps %xmm0,%xmm2 ## copy of k*rsq
    movaps %xmm3,%xmm5
    movaps %xmm6,%xmm8
    addps  %xmm1,%xmm2 ## rinv+krsq
    addps  %xmm4,%xmm5
    addps  %xmm7,%xmm8
    movaps nb203_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb203_qqH(%rsp),%xmm12
    movaps nb203_qqM(%rsp),%xmm13
    mulps  %xmm12,%xmm2 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm12,%xmm5 ## voul=qq*(rinv+ krsq-crf)
    mulps  %xmm13,%xmm8 ## voul=qq*(rinv+ krsq-crf)
    addps  %xmm0,%xmm0 ## 2*krsq
    addps  %xmm3,%xmm3
    addps  %xmm6,%xmm6
    subps  %xmm0,%xmm1 ## rinv-2*krsq
    subps  %xmm3,%xmm4
    subps  %xmm6,%xmm7
    mulps  %xmm12,%xmm1  ## (rinv-2*krsq)*qq
    mulps  %xmm12,%xmm4
    mulps  %xmm13,%xmm7
    addps  nb203_vctot(%rsp),%xmm2
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm2
    movaps %xmm2,nb203_vctot(%rsp)

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        ## move j forces to local temp variables 
    movlps (%rdi,%rax,4),%xmm9 ## jxa jya  -   -
    movlps (%rdi,%rcx,4),%xmm10 ## jxc jyc  -   -
    movhps (%rdi,%rbx,4),%xmm9 ## jxa jya jxb jyb 
    movhps (%rdi,%rdx,4),%xmm10 ## jxc jyc jxd jyd 

    movss  8(%rdi,%rax,4),%xmm11    ## jza  -  -  -
    movss  8(%rdi,%rcx,4),%xmm12    ## jzc  -  -  -
    movss  8(%rdi,%rbx,4),%xmm6     ## jzb
    movss  8(%rdi,%rdx,4),%xmm8     ## jzd
    movlhps %xmm6,%xmm11 ## jza  -  jzb  -
    movlhps %xmm8,%xmm12 ## jzc  -  jzd -

    shufps $136,%xmm12,%xmm11 ## 10001000 => jza jzb jzc jzd

    ## xmm9: jxa jya jxb jyb 
    ## xmm10: jxc jyc jxd jyd
    ## xmm11: jza jzb jzc jzd

    movaps %xmm1,%xmm0
    movaps %xmm1,%xmm2
    movaps %xmm4,%xmm3
    movaps %xmm4,%xmm5
    movaps %xmm7,%xmm6
    movaps %xmm7,%xmm8

        mulps nb203_dxH1(%rsp),%xmm0
        mulps nb203_dyH1(%rsp),%xmm1
        mulps nb203_dzH1(%rsp),%xmm2
        mulps nb203_dxH2(%rsp),%xmm3
        mulps nb203_dyH2(%rsp),%xmm4
        mulps nb203_dzH2(%rsp),%xmm5
        mulps nb203_dxM(%rsp),%xmm6
        mulps nb203_dyM(%rsp),%xmm7
        mulps nb203_dzM(%rsp),%xmm8

    movaps %xmm0,%xmm13
    movaps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb203_fixH1(%rsp),%xmm0
    addps nb203_fiyH1(%rsp),%xmm1
    addps nb203_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb203_fixH2(%rsp),%xmm3
    addps nb203_fiyH2(%rsp),%xmm4
    addps nb203_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb203_fixM(%rsp),%xmm6
    addps nb203_fiyM(%rsp),%xmm7
    addps nb203_fizM(%rsp),%xmm8

    movaps %xmm0,nb203_fixH1(%rsp)
    movaps %xmm1,nb203_fiyH1(%rsp)
    movaps %xmm2,nb203_fizH1(%rsp)
    movaps %xmm3,nb203_fixH2(%rsp)
    movaps %xmm4,nb203_fiyH2(%rsp)
    movaps %xmm5,nb203_fizH2(%rsp)
    movaps %xmm6,nb203_fixM(%rsp)
    movaps %xmm7,nb203_fiyM(%rsp)
    movaps %xmm8,nb203_fizM(%rsp)

    ## xmm9 = fjx
    ## xmm10 = fjy
    ## xmm11 = fjz
    movaps %xmm13,%xmm15
    unpcklps %xmm14,%xmm13
    unpckhps %xmm14,%xmm15

    addps %xmm13,%xmm9
    addps %xmm15,%xmm10

    movhlps  %xmm11,%xmm12 ## fjzc fjzd

    movlps %xmm9,(%rdi,%rax,4)
    movhps %xmm9,(%rdi,%rbx,4)
    movlps %xmm10,(%rdi,%rcx,4)
    movhps %xmm10,(%rdi,%rdx,4)
    movss  %xmm11,8(%rdi,%rax,4)
    movss  %xmm12,8(%rdi,%rcx,4)
    shufps $1,%xmm11,%xmm11
    shufps $1,%xmm12,%xmm12
    movss  %xmm11,8(%rdi,%rbx,4)
    movss  %xmm12,8(%rdi,%rdx,4)

        ## should we do one more iteration? 
        subl $4,nb203_innerk(%rsp)
        jl    _nb_kernel203_x86_64_sse.nb203_odd_inner
        jmp   _nb_kernel203_x86_64_sse.nb203_unroll_loop
_nb_kernel203_x86_64_sse.nb203_odd_inner: 
        addl $4,nb203_innerk(%rsp)
        jnz   _nb_kernel203_x86_64_sse.nb203_odd_loop
        jmp   _nb_kernel203_x86_64_sse.nb203_updateouterdata
_nb_kernel203_x86_64_sse.nb203_odd_loop: 
        movq  nb203_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb203_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb203_iqM(%rsp),%xmm4
        movq nb203_charge(%rbp),%rsi
        movhps nb203_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb203_qqM(%rsp)    ## use dummy qq for storage 

        movq nb203_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        movss nb203_ixM(%rsp),%xmm0
        movss nb203_iyM(%rsp),%xmm1
        movss nb203_izM(%rsp),%xmm2

        movlps nb203_ixH1(%rsp),%xmm6
        movlps nb203_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm0
        movlps nb203_iyH1(%rsp),%xmm6
        movlps nb203_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm1
        movlps nb203_izH1(%rsp),%xmm6
        movlps nb203_izH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm2

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use dummy dx for storage
        movaps %xmm3,nb203_dxM(%rsp)
        movaps %xmm4,nb203_dyM(%rsp)
        movaps %xmm5,nb203_dzM(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 

        movaps %xmm4,%xmm0
        mulps nb203_krf(%rsp),%xmm0
        movaps %xmm0,nb203_krsqM(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb203_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb203_half(%rsp),%xmm0
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

        movaps %xmm0,%xmm1      ## xmm1=rinv 
        movaps nb203_krsqM(%rsp),%xmm3
        addps  %xmm3,%xmm0      ## xmm0=rinv+ krsq 
        subps  nb203_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb203_two(%rsp),%xmm3
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq 
        mulps  nb203_qqM(%rsp),%xmm0    ## xmm0=vcoul 
        mulps  nb203_qqM(%rsp),%xmm1    ## xmm1=coul part of fs 


        mulps  %xmm1,%xmm4      ## xmm4=total fscal 
        addps  nb203_vctot(%rsp),%xmm0
        movaps %xmm0,nb203_vctot(%rsp)

        movaps nb203_dxM(%rsp),%xmm0
        movaps nb203_dyM(%rsp),%xmm1
        movaps nb203_dzM(%rsp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2
        ## xmm0-xmm2 contains tx-tz (partial force) 
        movss  nb203_fixM(%rsp),%xmm3
        movss  nb203_fiyM(%rsp),%xmm4
        movss  nb203_fizM(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb203_fixM(%rsp)
        movss  %xmm4,nb203_fiyM(%rsp)
        movss  %xmm5,nb203_fizM(%rsp)   ## updated the O force now do the H's 
        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $230,%xmm3,%xmm3 ## 11100110      ;# shift right 
        shufps $230,%xmm4,%xmm4 ## 11100110
        shufps $230,%xmm5,%xmm5 ## 11100110
        addss  nb203_fixH1(%rsp),%xmm3
        addss  nb203_fiyH1(%rsp),%xmm4
        addss  nb203_fizH1(%rsp),%xmm5
        movss  %xmm3,nb203_fixH1(%rsp)
        movss  %xmm4,nb203_fiyH1(%rsp)
        movss  %xmm5,nb203_fizH1(%rsp)          ## updated the H1 force 

        movq nb203_faction(%rbp),%rdi
        shufps $231,%xmm3,%xmm3 ## 11100111      ;# shift right 
        shufps $231,%xmm4,%xmm4 ## 11100111
        shufps $231,%xmm5,%xmm5 ## 11100111
        addss  nb203_fixH2(%rsp),%xmm3
        addss  nb203_fiyH2(%rsp),%xmm4
        addss  nb203_fizH2(%rsp),%xmm5
        movss  %xmm3,nb203_fixH2(%rsp)
        movss  %xmm4,nb203_fiyH2(%rsp)
        movss  %xmm5,nb203_fizH2(%rsp)          ## updated the H2 force 

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

        decl nb203_innerk(%rsp)
        jz    _nb_kernel203_x86_64_sse.nb203_updateouterdata
        jmp   _nb_kernel203_x86_64_sse.nb203_odd_loop
_nb_kernel203_x86_64_sse.nb203_updateouterdata: 
        movl  nb203_ii3(%rsp),%ecx
        movq  nb203_faction(%rbp),%rdi
        movq  nb203_fshift(%rbp),%rsi
        movl  nb203_is3(%rsp),%edx

        ## accumulate  H1 i forces in xmm0, xmm1, xmm2 
        movaps nb203_fixH1(%rsp),%xmm0
        movaps nb203_fiyH1(%rsp),%xmm1
        movaps nb203_fizH1(%rsp),%xmm2

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
        movaps nb203_fixH2(%rsp),%xmm0
        movaps nb203_fiyH2(%rsp),%xmm1
        movaps nb203_fizH2(%rsp),%xmm2

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

        ## accumulate m i forces in xmm0, xmm1, xmm2 
        movaps nb203_fixM(%rsp),%xmm0
        movaps nb203_fiyM(%rsp),%xmm1
        movaps nb203_fizM(%rsp),%xmm2

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
        movl nb203_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb203_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb203_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb203_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb203_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel203_x86_64_sse.nb203_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb203_n(%rsp)
        jmp _nb_kernel203_x86_64_sse.nb203_outer
_nb_kernel203_x86_64_sse.nb203_outerend: 
        ## check if more outer neighborlists remain
        movl  nb203_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel203_x86_64_sse.nb203_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel203_x86_64_sse.nb203_threadloop
_nb_kernel203_x86_64_sse.nb203_end: 

        movl nb203_nouter(%rsp),%eax
        movl nb203_ninner(%rsp),%ebx
        movq nb203_outeriter(%rbp),%rcx
        movq nb203_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $808,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret








.globl nb_kernel203nf_x86_64_sse
.globl _nb_kernel203nf_x86_64_sse
nb_kernel203nf_x86_64_sse:      
_nb_kernel203nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb203nf_fshift, 16
.set nb203nf_gid, 24
.set nb203nf_pos, 32
.set nb203nf_faction, 40
.set nb203nf_charge, 48
.set nb203nf_p_facel, 56
.set nb203nf_argkrf, 64
.set nb203nf_argcrf, 72
.set nb203nf_Vc, 80
.set nb203nf_type, 88
.set nb203nf_p_ntype, 96
.set nb203nf_vdwparam, 104
.set nb203nf_Vvdw, 112
.set nb203nf_p_tabscale, 120
.set nb203nf_VFtab, 128
.set nb203nf_invsqrta, 136
.set nb203nf_dvda, 144
.set nb203nf_p_gbtabscale, 152
.set nb203nf_GBtab, 160
.set nb203nf_p_nthreads, 168
.set nb203nf_count, 176
.set nb203nf_mtx, 184
.set nb203nf_outeriter, 192
.set nb203nf_inneriter, 200
.set nb203nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb203nf_ixH1, 0
.set nb203nf_iyH1, 16
.set nb203nf_izH1, 32
.set nb203nf_ixH2, 48
.set nb203nf_iyH2, 64
.set nb203nf_izH2, 80
.set nb203nf_ixM, 96
.set nb203nf_iyM, 112
.set nb203nf_izM, 128
.set nb203nf_iqH, 144
.set nb203nf_iqM, 160
.set nb203nf_qqH, 176
.set nb203nf_qqM, 192
.set nb203nf_vctot, 208
.set nb203nf_half, 224
.set nb203nf_three, 240
.set nb203nf_krf, 256
.set nb203nf_crf, 272
.set nb203nf_krsqH1, 288
.set nb203nf_krsqH2, 304
.set nb203nf_krsqM, 320
.set nb203nf_is3, 336
.set nb203nf_ii3, 340
.set nb203nf_innerjjnr, 344
.set nb203nf_nri, 352
.set nb203nf_iinr, 360
.set nb203nf_jindex, 368
.set nb203nf_jjnr, 376
.set nb203nf_shift, 384
.set nb203nf_shiftvec, 392
.set nb203nf_facel, 400
.set nb203nf_innerk, 408
.set nb203nf_n, 412
.set nb203nf_nn1, 416
.set nb203nf_nouter, 420
.set nb203nf_ninner, 424

        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $440,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb203nf_nouter(%rsp)
        movl %eax,nb203nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb203nf_nri(%rsp)
        movq %rsi,nb203nf_iinr(%rsp)
        movq %rdx,nb203nf_jindex(%rsp)
        movq %rcx,nb203nf_jjnr(%rsp)
        movq %r8,nb203nf_shift(%rsp)
        movq %r9,nb203nf_shiftvec(%rsp)
        movq nb203nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb203nf_facel(%rsp)

        movq nb203nf_argkrf(%rbp),%rsi
        movq nb203nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb203nf_krf(%rsp)
        movaps %xmm2,nb203nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb203nf_half(%rsp)
        movss nb203nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb203nf_half(%rsp)
        movaps %xmm3,nb203nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb203nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx),%ebx               ## ebx =ii 

        movq  nb203nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm3
        movss 12(%rdx,%rbx,4),%xmm4
        movq nb203nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb203nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb203nf_iqH(%rsp)
        movaps %xmm4,nb203nf_iqM(%rsp)

_nb_kernel203nf_x86_64_sse.nb203nf_threadloop: 
        movq  nb203nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel203nf_x86_64_sse.nb203nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel203nf_x86_64_sse.nb203nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb203nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb203nf_n(%rsp)
        movl %ebx,nb203nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel203nf_x86_64_sse.nb203nf_outerstart
        jmp _nb_kernel203nf_x86_64_sse.nb203nf_end

_nb_kernel203nf_x86_64_sse.nb203nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb203nf_nouter(%rsp),%ebx
        movl %ebx,nb203nf_nouter(%rsp)

_nb_kernel203nf_x86_64_sse.nb203nf_outer: 
        movq  nb203nf_shift(%rsp),%rax          ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb203nf_is3(%rsp)            ## store is3 

        movq  nb203nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb203nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb203nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb203nf_ii3(%rsp)

        addss 12(%rax,%rbx,4),%xmm3
        addss 16(%rax,%rbx,4),%xmm4
        addss 20(%rax,%rbx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb203nf_ixH1(%rsp)
        movaps %xmm4,nb203nf_iyH1(%rsp)
        movaps %xmm5,nb203nf_izH1(%rsp)

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
        movaps %xmm0,nb203nf_ixH2(%rsp)
        movaps %xmm1,nb203nf_iyH2(%rsp)
        movaps %xmm2,nb203nf_izH2(%rsp)
        movaps %xmm3,nb203nf_ixM(%rsp)
        movaps %xmm4,nb203nf_iyM(%rsp)
        movaps %xmm5,nb203nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb203nf_vctot(%rsp)

        movq  nb203nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb203nf_pos(%rbp),%rsi
        movq  nb203nf_faction(%rbp),%rdi
        movq  nb203nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb203nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb203nf_ninner(%rsp),%ecx
        movl  %ecx,nb203nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb203nf_innerk(%rsp)         ## number of innerloop atoms 
        jge   _nb_kernel203nf_x86_64_sse.nb203nf_unroll_loop
        jmp   _nb_kernel203nf_x86_64_sse.nb203nf_odd_inner
_nb_kernel203nf_x86_64_sse.nb203nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb203nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb203nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb203nf_charge(%rbp),%rsi  ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb203nf_iqH(%rsp),%xmm3
        mulps  nb203nf_iqM(%rsp),%xmm4

        movaps  %xmm3,nb203nf_qqH(%rsp)
        movaps  %xmm4,nb203nf_qqM(%rsp)

        movq nb203nf_pos(%rbp),%rsi     ## base of pos[] 

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
        movaps nb203nf_ixH1(%rsp),%xmm4
        movaps nb203nf_iyH1(%rsp),%xmm5
        movaps nb203nf_izH1(%rsp),%xmm6

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
        ## rsqH1 in xmm7 

        ## move ixH2-izH2 to xmm4-xmm6 
        movaps nb203nf_ixH2(%rsp),%xmm4
        movaps nb203nf_iyH2(%rsp),%xmm5
        movaps nb203nf_izH2(%rsp),%xmm6

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
        ## rsqH2 in xmm6 

        ## move ixM-izM to xmm3-xmm5  
        movaps nb203nf_ixM(%rsp),%xmm3
        movaps nb203nf_iyM(%rsp),%xmm4
        movaps nb203nf_izM(%rsp),%xmm5

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

        movaps %xmm5,%xmm0
        movaps %xmm6,%xmm1
        movaps %xmm7,%xmm2

        mulps  nb203nf_krf(%rsp),%xmm0
        mulps  nb203nf_krf(%rsp),%xmm1
        mulps  nb203nf_krf(%rsp),%xmm2

        movaps %xmm0,nb203nf_krsqM(%rsp)
        movaps %xmm1,nb203nf_krsqH2(%rsp)
        movaps %xmm2,nb203nf_krsqH1(%rsp)

        ## start with rsqH1 - seed in xmm2      
        rsqrtps %xmm7,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203nf_three(%rsp),%xmm4
        mulps   %xmm7,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm7     ## rinvH1 in xmm7 
        ## rsqH2 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203nf_three(%rsp),%xmm4
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm6     ## rinvH2 in xmm6 
        ## rsqM - seed in xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb203nf_three(%rsp),%xmm4
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm4     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm4     ## lu*(3-rsq*lu*lu) 
        mulps   nb203nf_half(%rsp),%xmm4
        movaps  %xmm4,%xmm5     ## rinvM in xmm5 

        ## do H1 interactions - xmm7=rinv
        addps nb203nf_krsqH1(%rsp),%xmm7
        subps nb203nf_crf(%rsp),%xmm7   ## xmm7=rinv+ krsq-crf 
        mulps nb203nf_qqH(%rsp),%xmm7
        addps nb203nf_vctot(%rsp),%xmm7

        ## H2 interactions - xmm6=rinv
        addps nb203nf_krsqH2(%rsp),%xmm6
        subps nb203nf_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps nb203nf_qqH(%rsp),%xmm6
        addps %xmm7,%xmm6

        ## M interactions - xmm5=rinv
        addps nb203nf_krsqM(%rsp),%xmm5
        subps nb203nf_crf(%rsp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps nb203nf_qqM(%rsp),%xmm5
        addps %xmm6,%xmm5
        movaps %xmm5,nb203nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb203nf_innerk(%rsp)
        jl    _nb_kernel203nf_x86_64_sse.nb203nf_odd_inner
        jmp   _nb_kernel203nf_x86_64_sse.nb203nf_unroll_loop
_nb_kernel203nf_x86_64_sse.nb203nf_odd_inner: 
        addl $4,nb203nf_innerk(%rsp)
        jnz   _nb_kernel203nf_x86_64_sse.nb203nf_odd_loop
        jmp   _nb_kernel203nf_x86_64_sse.nb203nf_updateouterdata
_nb_kernel203nf_x86_64_sse.nb203nf_odd_loop: 
        movq  nb203nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb203nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4
        movss nb203nf_iqM(%rsp),%xmm4
        movq nb203nf_charge(%rbp),%rsi
        movhps nb203nf_iqH(%rsp),%xmm4
        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb203nf_qqM(%rsp)          ## use dummy qq for storage 

        movq nb203nf_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

        movss nb203nf_ixM(%rsp),%xmm3
        movss nb203nf_iyM(%rsp),%xmm4
        movss nb203nf_izM(%rsp),%xmm5

        movlps nb203nf_ixH1(%rsp),%xmm6
        movlps nb203nf_ixH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm3
        movlps nb203nf_iyH1(%rsp),%xmm6
        movlps nb203nf_iyH2(%rsp),%xmm7
        unpcklps %xmm7,%xmm6
        movlhps %xmm6,%xmm4
        movlps nb203nf_izH1(%rsp),%xmm6
        movlps nb203nf_izH2(%rsp),%xmm7
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

        movaps %xmm4,%xmm0
        mulps nb203nf_krf(%rsp),%xmm0
        movaps %xmm0,nb203nf_krsqM(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb203nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb203nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv 
        ## a little trick to avoid NaNs: 
        ## positions 0,2,and 3 are valid, but not 1. 
        ## If it contains NaN it doesnt help to mult by 0, 
        ## So we shuffle it and copy pos 0 to pos1! 
        shufps $224,%xmm0,%xmm0 ## 11100000

        ## xmm0=rinv 
        addps  nb203nf_krsqM(%rsp),%xmm0
        subps  nb203nf_crf(%rsp),%xmm0   ## xmm0=rinv+ krsq-crf 
        mulps  nb203nf_qqM(%rsp),%xmm0          ## xmm0=vcoul 
        addps  nb203nf_vctot(%rsp),%xmm0
        movaps %xmm0,nb203nf_vctot(%rsp)

        decl nb203nf_innerk(%rsp)
        jz    _nb_kernel203nf_x86_64_sse.nb203nf_updateouterdata
        jmp   _nb_kernel203nf_x86_64_sse.nb203nf_odd_loop
_nb_kernel203nf_x86_64_sse.nb203nf_updateouterdata: 
        ## get n from stack
        movl nb203nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb203nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb203nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb203nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb203nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel203nf_x86_64_sse.nb203nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb203nf_n(%rsp)
        jmp _nb_kernel203nf_x86_64_sse.nb203nf_outer
_nb_kernel203nf_x86_64_sse.nb203nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb203nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel203nf_x86_64_sse.nb203nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel203nf_x86_64_sse.nb203nf_threadloop
_nb_kernel203nf_x86_64_sse.nb203nf_end: 

        movl nb203nf_nouter(%rsp),%eax
        movl nb203nf_ninner(%rsp),%ebx
        movq nb203nf_outeriter(%rbp),%rcx
        movq nb203nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $440,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

