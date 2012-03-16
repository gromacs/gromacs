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





.globl nb_kernel213_x86_64_sse
.globl _nb_kernel213_x86_64_sse
nb_kernel213_x86_64_sse:        
_nb_kernel213_x86_64_sse:       
##      Room for return address and rbp (16 bytes)
.set nb213_fshift, 16
.set nb213_gid, 24
.set nb213_pos, 32
.set nb213_faction, 40
.set nb213_charge, 48
.set nb213_p_facel, 56
.set nb213_argkrf, 64
.set nb213_argcrf, 72
.set nb213_Vc, 80
.set nb213_type, 88
.set nb213_p_ntype, 96
.set nb213_vdwparam, 104
.set nb213_Vvdw, 112
.set nb213_p_tabscale, 120
.set nb213_VFtab, 128
.set nb213_invsqrta, 136
.set nb213_dvda, 144
.set nb213_p_gbtabscale, 152
.set nb213_GBtab, 160
.set nb213_p_nthreads, 168
.set nb213_count, 176
.set nb213_mtx, 184
.set nb213_outeriter, 192
.set nb213_inneriter, 200
.set nb213_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb213_ixO, 0
.set nb213_iyO, 16
.set nb213_izO, 32
.set nb213_ixH1, 48
.set nb213_iyH1, 64
.set nb213_izH1, 80
.set nb213_ixH2, 96
.set nb213_iyH2, 112
.set nb213_izH2, 128
.set nb213_ixM, 144
.set nb213_iyM, 160
.set nb213_izM, 176
.set nb213_iqM, 192
.set nb213_iqH, 208
.set nb213_dxO, 224
.set nb213_dyO, 240
.set nb213_dzO, 256
.set nb213_dxH1, 272
.set nb213_dyH1, 288
.set nb213_dzH1, 304
.set nb213_dxH2, 320
.set nb213_dyH2, 336
.set nb213_dzH2, 352
.set nb213_dxM, 368
.set nb213_dyM, 384
.set nb213_dzM, 400
.set nb213_qqM, 416
.set nb213_qqH, 432
.set nb213_rinvH1, 448
.set nb213_rinvH2, 464
.set nb213_rinvM, 480
.set nb213_two, 496
.set nb213_c6, 512
.set nb213_c12, 528
.set nb213_six, 544
.set nb213_twelve, 560
.set nb213_krf, 576
.set nb213_crf, 592
.set nb213_krsqH1, 608
.set nb213_krsqH2, 624
.set nb213_krsqM, 640
.set nb213_vctot, 656
.set nb213_Vvdwtot, 672
.set nb213_fixO, 688
.set nb213_fiyO, 704
.set nb213_fizO, 720
.set nb213_fixH1, 736
.set nb213_fiyH1, 752
.set nb213_fizH1, 768
.set nb213_fixH2, 784
.set nb213_fiyH2, 800
.set nb213_fizH2, 816
.set nb213_fixM, 832
.set nb213_fiyM, 848
.set nb213_fizM, 864
.set nb213_fjx, 880
.set nb213_fjy, 896
.set nb213_fjz, 912
.set nb213_half, 928
.set nb213_three, 944
.set nb213_is3, 960
.set nb213_ii3, 964
.set nb213_nri, 968
.set nb213_iinr, 976
.set nb213_jindex, 984
.set nb213_jjnr, 992
.set nb213_shift, 1000
.set nb213_shiftvec, 1008
.set nb213_facel, 1016
.set nb213_innerjjnr, 1024
.set nb213_ntia, 1032
.set nb213_innerk, 1036
.set nb213_n, 1040
.set nb213_nn1, 1044
.set nb213_nouter, 1048
.set nb213_ninner, 1052

        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1064,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb213_nouter(%rsp)
        movl %eax,nb213_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb213_nri(%rsp)
        movq %rsi,nb213_iinr(%rsp)
        movq %rdx,nb213_jindex(%rsp)
        movq %rcx,nb213_jjnr(%rsp)
        movq %r8,nb213_shift(%rsp)
        movq %r9,nb213_shiftvec(%rsp)
        movq nb213_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb213_facel(%rsp)


        movq nb213_argkrf(%rbp),%rsi
        movq nb213_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb213_krf(%rsp)
        movaps %xmm2,nb213_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb213_half(%rsp)
        movss nb213_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm3,%xmm4
        addps  %xmm4,%xmm4      ## six
        movaps %xmm4,%xmm5
        addps  %xmm5,%xmm5      ## twelve
        movaps %xmm1,nb213_half(%rsp)
        movaps %xmm2,nb213_two(%rsp)
        movaps %xmm3,nb213_three(%rsp)
        movaps %xmm4,nb213_six(%rsp)
        movaps %xmm5,nb213_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb213_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb213_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm4
        movss 12(%rdx,%rbx,4),%xmm3
        movq nb213_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb213_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb213_iqM(%rsp)
        movaps %xmm4,nb213_iqH(%rsp)

        movq  nb213_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb213_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb213_ntia(%rsp)
_nb_kernel213_x86_64_sse.nb213_threadloop: 
        movq  nb213_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel213_x86_64_sse.nb213_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel213_x86_64_sse.nb213_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb213_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb213_n(%rsp)
        movl %ebx,nb213_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel213_x86_64_sse.nb213_outerstart
        jmp _nb_kernel213_x86_64_sse.nb213_end

_nb_kernel213_x86_64_sse.nb213_outerstart: 
        ## ebx contains number of outer iterations
        addl nb213_nouter(%rsp),%ebx
        movl %ebx,nb213_nouter(%rsp)

_nb_kernel213_x86_64_sse.nb213_outer: 
        movq  nb213_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb213_is3(%rsp)      ## store is3 

        movq  nb213_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb213_iinr(%rsp),%rcx             ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb213_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb213_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb213_ixO(%rsp)
        movaps %xmm4,nb213_iyO(%rsp)
        movaps %xmm5,nb213_izO(%rsp)
        movaps %xmm6,nb213_ixH1(%rsp)
        movaps %xmm7,nb213_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb213_izH1(%rsp)
        movaps %xmm0,nb213_ixH2(%rsp)
        movaps %xmm1,nb213_iyH2(%rsp)
        movaps %xmm2,nb213_izH2(%rsp)
        movaps %xmm3,nb213_ixM(%rsp)
        movaps %xmm4,nb213_iyM(%rsp)
        movaps %xmm5,nb213_izM(%rsp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb213_vctot(%rsp)
        movaps %xmm4,nb213_Vvdwtot(%rsp)
        movaps %xmm4,nb213_fixO(%rsp)
        movaps %xmm4,nb213_fiyO(%rsp)
        movaps %xmm4,nb213_fizO(%rsp)
        movaps %xmm4,nb213_fixH1(%rsp)
        movaps %xmm4,nb213_fiyH1(%rsp)
        movaps %xmm4,nb213_fizH1(%rsp)
        movaps %xmm4,nb213_fixH2(%rsp)
        movaps %xmm4,nb213_fiyH2(%rsp)
        movaps %xmm4,nb213_fizH2(%rsp)
        movaps %xmm4,nb213_fixM(%rsp)
        movaps %xmm4,nb213_fiyM(%rsp)
        movaps %xmm4,nb213_fizM(%rsp)

        movq  nb213_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb213_pos(%rbp),%rsi
        movq  nb213_faction(%rbp),%rdi
        movq  nb213_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb213_innerjjnr(%rsp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb213_ninner(%rsp),%ecx
        movl  %ecx,nb213_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb213_innerk(%rsp)   ## number of innerloop atoms 
        jge   _nb_kernel213_x86_64_sse.nb213_unroll_loop
        jmp   _nb_kernel213_x86_64_sse.nb213_odd_inner
_nb_kernel213_x86_64_sse.nb213_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb213_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb213_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb213_charge(%rbp),%rsi    ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb213_iqM(%rsp),%xmm3
        mulps  nb213_iqH(%rsp),%xmm4

        movaps  %xmm3,nb213_qqM(%rsp)
        movaps  %xmm4,nb213_qqH(%rsp)

        movq nb213_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movl (%rsi,%rcx,4),%r10d
        movl (%rsi,%rdx,4),%r11d
        movq nb213_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        shll %r10d
        shll %r11d
        movl nb213_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d
        addl %edi,%r10d
        addl %edi,%r11d

        movlps (%rsi,%r8,4),%xmm6
        movlps (%rsi,%r10,4),%xmm7
        movhps (%rsi,%r9,4),%xmm6
        movhps (%rsi,%r11,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movaps %xmm4,nb213_c6(%rsp)
        movaps %xmm6,nb213_c12(%rsp)

        movq nb213_pos(%rbp),%rsi       ## base of pos[] 

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

    ## O interaction
    ## copy to xmm3-xmm5
    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5

    subps nb213_ixO(%rsp),%xmm3
    subps nb213_iyO(%rsp),%xmm4
    subps nb213_izO(%rsp),%xmm5

    movaps %xmm3,%xmm13
    movaps %xmm4,%xmm14
    movaps %xmm5,%xmm15

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3

    ## calc 1/rsq
    rcpps %xmm3,%xmm5
    movaps nb213_two(%rsp),%xmm4
    mulps %xmm5,%xmm3
    subps %xmm3,%xmm4
    mulps %xmm5,%xmm4       ## xmm4=rinvsq

    movaps %xmm4,%xmm3      ## rinvsq
    mulps  %xmm4,%xmm4      ## rinv4
    mulps  %xmm3,%xmm4      ## rinv6
    movaps %xmm4,%xmm5
    mulps  %xmm5,%xmm5      ## rinv12
    mulps  nb213_c6(%rsp),%xmm4
    mulps  nb213_c12(%rsp),%xmm5
    movaps %xmm5,%xmm6
    subps  %xmm4,%xmm6 ## Vvdw=vvdw12-vvdw6
    mulps  nb213_six(%rsp),%xmm4
    mulps  nb213_twelve(%rsp),%xmm5
    subps  %xmm4,%xmm5
    mulps  %xmm5,%xmm3  ## fscal

    addps  nb213_Vvdwtot(%rsp),%xmm6
    movaps %xmm6,nb213_Vvdwtot(%rsp)

    mulps  %xmm3,%xmm13 ## fx
    mulps  %xmm3,%xmm14 ## fy
    mulps  %xmm3,%xmm15 ## fz

    ## save j force temporarily
    movaps %xmm13,nb213_fjx(%rsp)
    movaps %xmm14,nb213_fjy(%rsp)
    movaps %xmm15,nb213_fjz(%rsp)

    ## increment i O force
    addps nb213_fixO(%rsp),%xmm13
    addps nb213_fiyO(%rsp),%xmm14
    addps nb213_fizO(%rsp),%xmm15
    movaps %xmm13,nb213_fixO(%rsp)
    movaps %xmm14,nb213_fiyO(%rsp)
    movaps %xmm15,nb213_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.        
    movaps %xmm0,%xmm3
    movaps %xmm1,%xmm4
    movaps %xmm2,%xmm5
    movaps %xmm0,%xmm6
    movaps %xmm1,%xmm7
    movaps %xmm2,%xmm8

    subps nb213_ixH1(%rsp),%xmm0
    subps nb213_iyH1(%rsp),%xmm1
    subps nb213_izH1(%rsp),%xmm2
    subps nb213_ixH2(%rsp),%xmm3
    subps nb213_iyH2(%rsp),%xmm4
    subps nb213_izH2(%rsp),%xmm5
    subps nb213_ixM(%rsp),%xmm6
    subps nb213_iyM(%rsp),%xmm7
    subps nb213_izM(%rsp),%xmm8

        movaps %xmm0,nb213_dxH1(%rsp)
        movaps %xmm1,nb213_dyH1(%rsp)
        movaps %xmm2,nb213_dzH1(%rsp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb213_dxH2(%rsp)
        movaps %xmm4,nb213_dyH2(%rsp)
        movaps %xmm5,nb213_dzH2(%rsp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        movaps %xmm6,nb213_dxM(%rsp)
        movaps %xmm7,nb213_dyM(%rsp)
        movaps %xmm8,nb213_dzM(%rsp)
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

        movaps  nb213_three(%rsp),%xmm9
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

        movaps  nb213_half(%rsp),%xmm4
        mulps   %xmm4,%xmm9 ## rinvH1
        mulps   %xmm4,%xmm10 ## rinvH2
    mulps   %xmm4,%xmm11 ## rinvM

        ## interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movaps %xmm9,%xmm1 ## copy of rinv
    movaps %xmm10,%xmm4
    movaps %xmm11,%xmm7
    movaps nb213_krf(%rsp),%xmm2
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
    movaps nb213_crf(%rsp),%xmm14
    subps  %xmm14,%xmm2  ## rinv+krsq-crf
    subps  %xmm14,%xmm5
    subps  %xmm14,%xmm8
    movaps nb213_qqH(%rsp),%xmm12
    movaps nb213_qqM(%rsp),%xmm13
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
    addps  nb213_vctot(%rsp),%xmm2
    addps  %xmm8,%xmm5
    addps  %xmm5,%xmm2
    movaps %xmm2,nb213_vctot(%rsp)

    mulps  %xmm9,%xmm1  ## fscal
    mulps  %xmm10,%xmm4
    mulps  %xmm11,%xmm7

        movq  nb213_faction(%rbp),%rdi
        ## move j  forces to local temp variables 
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

        mulps nb213_dxH1(%rsp),%xmm0
        mulps nb213_dyH1(%rsp),%xmm1
        mulps nb213_dzH1(%rsp),%xmm2
        mulps nb213_dxH2(%rsp),%xmm3
        mulps nb213_dyH2(%rsp),%xmm4
        mulps nb213_dzH2(%rsp),%xmm5
        mulps nb213_dxM(%rsp),%xmm6
        mulps nb213_dyM(%rsp),%xmm7
        mulps nb213_dzM(%rsp),%xmm8

    ## fetch forces from O interaction
    movaps nb213_fjx(%rsp),%xmm13
    movaps nb213_fjy(%rsp),%xmm14
    addps  nb213_fjz(%rsp),%xmm11

    addps %xmm0,%xmm13
    addps %xmm1,%xmm14
    addps %xmm2,%xmm11
    addps nb213_fixH1(%rsp),%xmm0
    addps nb213_fiyH1(%rsp),%xmm1
    addps nb213_fizH1(%rsp),%xmm2

    addps %xmm3,%xmm13
    addps %xmm4,%xmm14
    addps %xmm5,%xmm11
    addps nb213_fixH2(%rsp),%xmm3
    addps nb213_fiyH2(%rsp),%xmm4
    addps nb213_fizH2(%rsp),%xmm5

    addps %xmm6,%xmm13
    addps %xmm7,%xmm14
    addps %xmm8,%xmm11
    addps nb213_fixM(%rsp),%xmm6
    addps nb213_fiyM(%rsp),%xmm7
    addps nb213_fizM(%rsp),%xmm8

    movaps %xmm0,nb213_fixH1(%rsp)
    movaps %xmm1,nb213_fiyH1(%rsp)
    movaps %xmm2,nb213_fizH1(%rsp)
    movaps %xmm3,nb213_fixH2(%rsp)
    movaps %xmm4,nb213_fiyH2(%rsp)
    movaps %xmm5,nb213_fizH2(%rsp)
    movaps %xmm6,nb213_fixM(%rsp)
    movaps %xmm7,nb213_fiyM(%rsp)
    movaps %xmm8,nb213_fizM(%rsp)

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
        subl $4,nb213_innerk(%rsp)
        jl    _nb_kernel213_x86_64_sse.nb213_odd_inner
        jmp   _nb_kernel213_x86_64_sse.nb213_unroll_loop
_nb_kernel213_x86_64_sse.nb213_odd_inner: 
        addl $4,nb213_innerk(%rsp)
        jnz   _nb_kernel213_x86_64_sse.nb213_odd_loop
        jmp   _nb_kernel213_x86_64_sse.nb213_updateouterdata
_nb_kernel213_x86_64_sse.nb213_odd_loop: 
        movq  nb213_innerjjnr(%rsp),%rdx        ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb213_innerjjnr(%rsp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb213_iqM(%rsp),%xmm4
        movq nb213_charge(%rbp),%rsi
        movhps nb213_iqH(%rsp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb213_qqM(%rsp)    ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movq nb213_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb213_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb213_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb213_c6(%rsp)
        movaps %xmm7,nb213_c12(%rsp)

        movq nb213_pos(%rbp),%rsi
        lea (%rax,%rax,2),%rax

        movss nb213_ixO(%rsp),%xmm0
        movss nb213_iyO(%rsp),%xmm1
        movss nb213_izO(%rsp),%xmm2
        movss nb213_ixH1(%rsp),%xmm3
        movss nb213_iyH1(%rsp),%xmm4
        movss nb213_izH1(%rsp),%xmm5
        unpcklps nb213_ixH2(%rsp),%xmm0         ## ixO ixH2 - -
        unpcklps nb213_iyH2(%rsp),%xmm1         ## iyO iyH2 - -
        unpcklps nb213_izH2(%rsp),%xmm2         ## izO izH2 - -
        unpcklps nb213_ixM(%rsp),%xmm3          ## ixH1 ixM - -
        unpcklps nb213_iyM(%rsp),%xmm4          ## iyH1 iyM - -
        unpcklps nb213_izM(%rsp),%xmm5          ## izH1 izM - -
        unpcklps %xmm3,%xmm0    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm4,%xmm1    ## same for y
        unpcklps %xmm5,%xmm2    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm3
        movss 4(%rsi,%rax,4),%xmm4
        movss 8(%rsi,%rax,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5

        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5

        ## use O distances for storage
        movaps %xmm3,nb213_dxO(%rsp)
        movaps %xmm4,nb213_dyO(%rsp)
        movaps %xmm5,nb213_dzO(%rsp)

        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5

        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        ## rsq in xmm4 
        movaps %xmm4,%xmm0
        mulps nb213_krf(%rsp),%xmm0
        movaps %xmm0,nb213_krsqM(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb213_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb213_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv


        movaps %xmm0,%xmm4
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm5
        mulps  %xmm4,%xmm4      ## xmm1=rinv, xmm4=rinvsq
        movaps nb213_krsqM(%rsp),%xmm3
        addps  %xmm3,%xmm5      ## xmm0=rinv+ krsq 
        subps  nb213_crf(%rsp),%xmm5   ## xmm0=rinv+ krsq-crf 
        mulps  nb213_two(%rsp),%xmm3
        subps  %xmm3,%xmm1      ## xmm1=rinv-2*krsq
        movaps %xmm5,%xmm7
        mulps  nb213_qqM(%rsp),%xmm7    ## xmm0=vcoul 
        mulps  nb213_qqM(%rsp),%xmm1    ## xmm1=coul part of fs 
        movaps %xmm1,%xmm6

        addps  nb213_vctot(%rsp),%xmm7
        movaps %xmm7,nb213_vctot(%rsp)

        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1
        movaps %xmm1,%xmm2
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb213_c6(%rsp),%xmm1
        mulss  nb213_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        mulss  nb213_six(%rsp),%xmm1
        mulss  nb213_twelve(%rsp),%xmm4
        subss  %xmm1,%xmm4
        addss  nb213_Vvdwtot(%rsp),%xmm3
        movss  %xmm3,nb213_Vvdwtot(%rsp)
        addps  %xmm6,%xmm4
        mulps  %xmm0,%xmm4
        mulps  %xmm0,%xmm4      ## fscal

        movaps nb213_dxO(%rsp),%xmm0
        movaps nb213_dyO(%rsp),%xmm1
        movaps nb213_dzO(%rsp),%xmm2

        mulps  %xmm4,%xmm0
        mulps  %xmm4,%xmm1
        mulps  %xmm4,%xmm2 ## xmm0-xmm2 now contains tx-tz (partial force)

        movss  nb213_fixO(%rsp),%xmm3
        movss  nb213_fiyO(%rsp),%xmm4
        movss  nb213_fizO(%rsp),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,nb213_fixO(%rsp)
        movss  %xmm4,nb213_fiyO(%rsp)
        movss  %xmm5,nb213_fizO(%rsp)   ## updated the O force now do the H's

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb213_fixH1(%rsp),%xmm3
        addss  nb213_fiyH1(%rsp),%xmm4
        addss  nb213_fizH1(%rsp),%xmm5
        movss  %xmm3,nb213_fixH1(%rsp)
        movss  %xmm4,nb213_fiyH1(%rsp)
        movss  %xmm5,nb213_fizH1(%rsp)          ## updated the H1 force 

        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb213_fixH2(%rsp),%xmm3
        addss  nb213_fiyH2(%rsp),%xmm4
        addss  nb213_fizH2(%rsp),%xmm5
        movss  %xmm3,nb213_fixH2(%rsp)
        movss  %xmm4,nb213_fiyH2(%rsp)
        movss  %xmm5,nb213_fizH2(%rsp)          ## updated the H2 force 

        movq nb213_faction(%rbp),%rdi
        shufps $0x39,%xmm3,%xmm3
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  nb213_fixM(%rsp),%xmm3
        addss  nb213_fiyM(%rsp),%xmm4
        addss  nb213_fizM(%rsp),%xmm5
        movss  %xmm3,nb213_fixM(%rsp)
        movss  %xmm4,nb213_fiyM(%rsp)
        movss  %xmm5,nb213_fizM(%rsp)   ## updated the M force 

        ## the fj's - move in from mem start by acc. tx/ty/tz in xmm0, xmm1
        movlps (%rdi,%rax,4),%xmm6
        movss  8(%rdi,%rax,4),%xmm7

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addps   %xmm0,%xmm3
        addps   %xmm1,%xmm4
        addps   %xmm2,%xmm5
        movaps  %xmm3,%xmm0
        movaps  %xmm4,%xmm1
        movaps  %xmm5,%xmm2

        shufps $0x39,%xmm3,%xmm3 ## shift right 
        shufps $0x39,%xmm4,%xmm4
        shufps $0x39,%xmm5,%xmm5
        addss  %xmm3,%xmm0
        addss  %xmm4,%xmm1
        addss  %xmm5,%xmm2
        unpcklps %xmm1,%xmm0    ## x,y sum in xmm0, z sum in xmm2

        addps    %xmm0,%xmm6
        addss    %xmm2,%xmm7

        movlps %xmm6,(%rdi,%rax,4)
        movss  %xmm7,8(%rdi,%rax,4)

        decl nb213_innerk(%rsp)
        jz    _nb_kernel213_x86_64_sse.nb213_updateouterdata
        jmp   _nb_kernel213_x86_64_sse.nb213_odd_loop
_nb_kernel213_x86_64_sse.nb213_updateouterdata: 
        movl  nb213_ii3(%rsp),%ecx
        movq  nb213_faction(%rbp),%rdi
        movq  nb213_fshift(%rbp),%rsi
        movl  nb213_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb213_fixO(%rsp),%xmm0
        movaps nb213_fiyO(%rsp),%xmm1
        movaps nb213_fizO(%rsp),%xmm2

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
        movss  (%rdi,%rcx,4),%xmm3
        movss  4(%rdi,%rcx,4),%xmm4
        movss  8(%rdi,%rcx,4),%xmm5
        subss  %xmm0,%xmm3
        subss  %xmm1,%xmm4
        subss  %xmm2,%xmm5
        movss  %xmm3,(%rdi,%rcx,4)
        movss  %xmm4,4(%rdi,%rcx,4)
        movss  %xmm5,8(%rdi,%rcx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## 00001000       

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb213_fixH1(%rsp),%xmm0
        movaps nb213_fiyH1(%rsp),%xmm1
        movaps nb213_fizH1(%rsp),%xmm2

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
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## 00001000       
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb213_fixH2(%rsp),%xmm0
        movaps nb213_fiyH2(%rsp),%xmm1
        movaps nb213_fizH2(%rsp),%xmm2

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

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movaps nb213_fixM(%rsp),%xmm0
        movaps nb213_fiyM(%rsp),%xmm1
        movaps nb213_fizM(%rsp),%xmm2

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
        movl nb213_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb213_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb213_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb213_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb213_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb213_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb213_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel213_x86_64_sse.nb213_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb213_n(%rsp)
        jmp _nb_kernel213_x86_64_sse.nb213_outer
_nb_kernel213_x86_64_sse.nb213_outerend: 
        ## check if more outer neighborlists remain
        movl  nb213_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel213_x86_64_sse.nb213_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel213_x86_64_sse.nb213_threadloop
_nb_kernel213_x86_64_sse.nb213_end: 
        movl nb213_nouter(%rsp),%eax
        movl nb213_ninner(%rsp),%ebx
        movq nb213_outeriter(%rbp),%rcx
        movq nb213_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1064,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret








.globl nb_kernel213nf_x86_64_sse
.globl _nb_kernel213nf_x86_64_sse
nb_kernel213nf_x86_64_sse:      
_nb_kernel213nf_x86_64_sse:     
##      Room for return address and rbp (16 bytes)
.set nb213nf_fshift, 16
.set nb213nf_gid, 24
.set nb213nf_pos, 32
.set nb213nf_faction, 40
.set nb213nf_charge, 48
.set nb213nf_p_facel, 56
.set nb213nf_argkrf, 64
.set nb213nf_argcrf, 72
.set nb213nf_Vc, 80
.set nb213nf_type, 88
.set nb213nf_p_ntype, 96
.set nb213nf_vdwparam, 104
.set nb213nf_Vvdw, 112
.set nb213nf_p_tabscale, 120
.set nb213nf_VFtab, 128
.set nb213nf_invsqrta, 136
.set nb213nf_dvda, 144
.set nb213nf_p_gbtabscale, 152
.set nb213nf_GBtab, 160
.set nb213nf_p_nthreads, 168
.set nb213nf_count, 176
.set nb213nf_mtx, 184
.set nb213nf_outeriter, 192
.set nb213nf_inneriter, 200
.set nb213nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb213nf_ixO, 0
.set nb213nf_iyO, 16
.set nb213nf_izO, 32
.set nb213nf_ixH1, 48
.set nb213nf_iyH1, 64
.set nb213nf_izH1, 80
.set nb213nf_ixH2, 96
.set nb213nf_iyH2, 112
.set nb213nf_izH2, 128
.set nb213nf_ixM, 144
.set nb213nf_iyM, 160
.set nb213nf_izM, 176
.set nb213nf_iqM, 192
.set nb213nf_iqH, 208
.set nb213nf_qqM, 224
.set nb213nf_qqH, 240
.set nb213nf_rinvH1, 256
.set nb213nf_rinvH2, 272
.set nb213nf_rinvM, 288
.set nb213nf_two, 304
.set nb213nf_c6, 320
.set nb213nf_c12, 336
.set nb213nf_krf, 352
.set nb213nf_crf, 368
.set nb213nf_krsqH1, 384
.set nb213nf_krsqH2, 400
.set nb213nf_krsqM, 416
.set nb213nf_vctot, 432
.set nb213nf_Vvdwtot, 448
.set nb213nf_half, 464
.set nb213nf_three, 480
.set nb213nf_nri, 496
.set nb213nf_iinr, 504
.set nb213nf_jindex, 512
.set nb213nf_jjnr, 520
.set nb213nf_shift, 528
.set nb213nf_shiftvec, 536
.set nb213nf_facel, 544
.set nb213nf_innerjjnr, 552
.set nb213nf_is3, 560
.set nb213nf_ii3, 564
.set nb213nf_ntia, 568
.set nb213nf_innerk, 572
.set nb213nf_n, 576
.set nb213nf_nn1, 580
.set nb213nf_nouter, 584
.set nb213nf_ninner, 588


        push %rbp
        movq %rsp,%rbp
        push %rbx


        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $600,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb213nf_nouter(%rsp)
        movl %eax,nb213nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb213nf_nri(%rsp)
        movq %rsi,nb213nf_iinr(%rsp)
        movq %rdx,nb213nf_jindex(%rsp)
        movq %rcx,nb213nf_jjnr(%rsp)
        movq %r8,nb213nf_shift(%rsp)
        movq %r9,nb213nf_shiftvec(%rsp)
        movq nb213nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss %xmm0,nb213nf_facel(%rsp)


        movq nb213nf_argkrf(%rbp),%rsi
        movq nb213nf_argcrf(%rbp),%rdi
        movss (%rsi),%xmm1
        movss (%rdi),%xmm2
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        movaps %xmm1,nb213nf_krf(%rsp)
        movaps %xmm2,nb213nf_crf(%rsp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## half in IEEE (hex)
        movl %eax,nb213nf_half(%rsp)
        movss nb213nf_half(%rsp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## one
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## two
        addps  %xmm2,%xmm3      ## three
        movaps %xmm1,nb213nf_half(%rsp)
        movaps %xmm2,nb213nf_two(%rsp)
        movaps %xmm3,nb213nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb213nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb213nf_charge(%rbp),%rdx
        movss 4(%rdx,%rbx,4),%xmm4
        movss 12(%rdx,%rbx,4),%xmm3
        movq nb213nf_p_facel(%rbp),%rsi
        movss (%rsi),%xmm0
        movss nb213nf_facel(%rsp),%xmm5
        mulss  %xmm5,%xmm3
        mulss  %xmm5,%xmm4

        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        movaps %xmm3,nb213nf_iqM(%rsp)
        movaps %xmm4,nb213nf_iqH(%rsp)

        movq  nb213nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb213nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb213nf_ntia(%rsp)

_nb_kernel213nf_x86_64_sse.nb213nf_threadloop: 
        movq  nb213nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel213nf_x86_64_sse.nb213nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addq  $1,%rbx                          ## rbx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel213nf_x86_64_sse.nb213nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb213nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb213nf_n(%rsp)
        movl %ebx,nb213nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel213nf_x86_64_sse.nb213nf_outerstart
        jmp _nb_kernel213nf_x86_64_sse.nb213nf_end

_nb_kernel213nf_x86_64_sse.nb213nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb213nf_nouter(%rsp),%ebx
        movl %ebx,nb213nf_nouter(%rsp)

_nb_kernel213nf_x86_64_sse.nb213nf_outer: 
        movq  nb213nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx                ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx        ## rbx=3*is 
        movl  %ebx,nb213nf_is3(%rsp)            ## store is3 

        movq  nb213nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movss (%rax,%rbx,4),%xmm0
        movss 4(%rax,%rbx,4),%xmm1
        movss 8(%rax,%rbx,4),%xmm2

        movq  nb213nf_iinr(%rsp),%rcx           ## rcx = pointer into iinr[]    
        movl  (%rcx,%rsi,4),%ebx                ## ebx =ii 

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        movaps %xmm0,%xmm6
        movaps %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb213nf_pos(%rbp),%rax    ## rax = base of pos[]  
        movl  %ebx,nb213nf_ii3(%rsp)

        addss (%rax,%rbx,4),%xmm3       ## ox
        addss 4(%rax,%rbx,4),%xmm4     ## oy
        addss 8(%rax,%rbx,4),%xmm5     ## oz
        addss 12(%rax,%rbx,4),%xmm6    ## h1x
        addss 16(%rax,%rbx,4),%xmm7    ## h1y
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm7,%xmm7
        movaps %xmm3,nb213nf_ixO(%rsp)
        movaps %xmm4,nb213nf_iyO(%rsp)
        movaps %xmm5,nb213nf_izO(%rsp)
        movaps %xmm6,nb213nf_ixH1(%rsp)
        movaps %xmm7,nb213nf_iyH1(%rsp)

        movss %xmm2,%xmm6
        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 20(%rax,%rbx,4),%xmm6    ## h1z
        addss 24(%rax,%rbx,4),%xmm0    ## h2x
        addss 28(%rax,%rbx,4),%xmm1    ## h2y
        addss 32(%rax,%rbx,4),%xmm2    ## h2z
        addss 36(%rax,%rbx,4),%xmm3    ## mx
        addss 40(%rax,%rbx,4),%xmm4    ## my
        addss 44(%rax,%rbx,4),%xmm5    ## mz

        shufps $0,%xmm6,%xmm6
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm6,nb213nf_izH1(%rsp)
        movaps %xmm0,nb213nf_ixH2(%rsp)
        movaps %xmm1,nb213nf_iyH2(%rsp)
        movaps %xmm2,nb213nf_izH2(%rsp)
        movaps %xmm3,nb213nf_ixM(%rsp)
        movaps %xmm4,nb213nf_iyM(%rsp)
        movaps %xmm5,nb213nf_izM(%rsp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb213nf_vctot(%rsp)
        movaps %xmm4,nb213nf_Vvdwtot(%rsp)

        movq  nb213nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx                ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movq  nb213nf_pos(%rbp),%rsi
        movq  nb213nf_faction(%rbp),%rdi
        movq  nb213nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb213nf_innerjjnr(%rsp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb213nf_ninner(%rsp),%ecx
        movl  %ecx,nb213nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb213nf_innerk(%rsp)         ## number of innerloop atoms 


        jge   _nb_kernel213nf_x86_64_sse.nb213nf_unroll_loop
        jmp   _nb_kernel213nf_x86_64_sse.nb213nf_odd_inner
_nb_kernel213nf_x86_64_sse.nb213nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movq  nb213nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx
        movl  8(%rdx),%ecx
        movl  12(%rdx),%edx             ## eax-edx=jnr1-4 

        addq $16,nb213nf_innerjjnr(%rsp)             ## advance pointer (unrolled 4) 

        movq nb213nf_charge(%rbp),%rsi  ## base of charge[] 

        movss (%rsi,%rax,4),%xmm3
        movss (%rsi,%rcx,4),%xmm4
        movss (%rsi,%rbx,4),%xmm6
        movss (%rsi,%rdx,4),%xmm7

        shufps $0,%xmm6,%xmm3
        shufps $0,%xmm7,%xmm4
        shufps $136,%xmm4,%xmm3 ## 10001000 ;# all charges in xmm3  
        movaps %xmm3,%xmm4              ## and in xmm4 
        mulps  nb213nf_iqM(%rsp),%xmm3
        mulps  nb213nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1
        movd  %ecx,%mm2
        movd  %edx,%mm3

        movaps  %xmm3,nb213nf_qqM(%rsp)
        movaps  %xmm4,nb213nf_qqH(%rsp)

        movq nb213nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movl (%rsi,%rcx,4),%ecx
        movl (%rsi,%rdx,4),%edx
        movq nb213nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        shll %ecx
        shll %edx
        movl nb213nf_ntia(%rsp),%edi
        addq %rdi,%rax
        addq %rdi,%rbx
        addq %rdi,%rcx
        addq %rdi,%rdx

        movlps (%rsi,%rax,4),%xmm6
        movlps (%rsi,%rcx,4),%xmm7
        movhps (%rsi,%rbx,4),%xmm6
        movhps (%rsi,%rdx,4),%xmm7

        movaps %xmm6,%xmm4
        shufps $136,%xmm7,%xmm4 ## 10001000
        shufps $221,%xmm7,%xmm6 ## 11011101

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movd  %mm2,%ecx
        movd  %mm3,%edx

        movaps %xmm4,nb213nf_c6(%rsp)
        movaps %xmm6,nb213nf_c12(%rsp)

        movq nb213nf_pos(%rbp),%rsi     ## base of pos[] 

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

        ## move ixO-izO to xmm4-xmm6 
        movaps nb213nf_ixO(%rsp),%xmm4
        movaps nb213nf_iyO(%rsp),%xmm5
        movaps nb213nf_izO(%rsp),%xmm6

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

        ## move ixH1-izH1 to xmm4-xmm6 
        movaps nb213nf_ixH1(%rsp),%xmm4
        movaps nb213nf_iyH1(%rsp),%xmm5
        movaps nb213nf_izH1(%rsp),%xmm6

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

        ## move ixH2-izH2 to xmm3-xmm5  
        movaps nb213nf_ixH2(%rsp),%xmm3
        movaps nb213nf_iyH2(%rsp),%xmm4
        movaps nb213nf_izH2(%rsp),%xmm5

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

        ## move ixM-izM to xmm2-xmm4  
        movaps nb213nf_iyM(%rsp),%xmm3
        movaps nb213nf_izM(%rsp),%xmm4
        subps  %xmm1,%xmm3
        subps  %xmm2,%xmm4
        movaps nb213nf_ixM(%rsp),%xmm2
        subps  %xmm0,%xmm2

        ## square it 
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        addps %xmm3,%xmm4
        addps %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 
        movaps %xmm4,%xmm0
        movaps %xmm5,%xmm1
        movaps %xmm6,%xmm2
        mulps  nb213nf_krf(%rsp),%xmm0
        mulps  nb213nf_krf(%rsp),%xmm1
        mulps  nb213nf_krf(%rsp),%xmm2
        movaps %xmm0,nb213nf_krsqM(%rsp)
        movaps %xmm1,nb213nf_krsqH2(%rsp)
        movaps %xmm2,nb213nf_krsqH1(%rsp)

        ## rsqH1 - seed in xmm2 
        rsqrtps %xmm6,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213nf_three(%rsp),%xmm0
        mulps   %xmm6,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213nf_half(%rsp),%xmm0
        movaps  %xmm0,nb213nf_rinvH1(%rsp)      ## rinvH1 

        ## rsqH2 - seed to xmm2 
        rsqrtps %xmm5,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213nf_three(%rsp),%xmm0
        mulps   %xmm5,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213nf_half(%rsp),%xmm0
        movaps  %xmm0,nb213nf_rinvH2(%rsp)      ## rinvH2 

        ## rsqM - seed to xmm2 
        rsqrtps %xmm4,%xmm2
        movaps  %xmm2,%xmm3
        mulps   %xmm2,%xmm2
        movaps  nb213nf_three(%rsp),%xmm0
        mulps   %xmm4,%xmm2     ## rsq*lu*lu 
        subps   %xmm2,%xmm0     ## 30-rsq*lu*lu 
        mulps   %xmm3,%xmm0     ## lu*(3-rsq*lu*lu) 
        mulps   nb213nf_half(%rsp),%xmm0
        movaps  %xmm0,nb213nf_rinvM(%rsp)

        ## Do the O LJ-only interaction directly.       
        rcpps   %xmm7,%xmm2
        movaps  nb213nf_two(%rsp),%xmm1
        mulps   %xmm2,%xmm7
        subps   %xmm7,%xmm1
        mulps   %xmm1,%xmm2 ## rinvsq 
        movaps  %xmm2,%xmm0
        mulps   %xmm2,%xmm0     ## r4
        mulps   %xmm2,%xmm0     ## r6
        movaps  %xmm0,%xmm1
        mulps   %xmm1,%xmm1     ## r12
        mulps   nb213nf_c6(%rsp),%xmm0
        mulps   nb213nf_c12(%rsp),%xmm1
        movaps  %xmm1,%xmm3
        subps   %xmm0,%xmm3     ## Vvdw12-Vvdw6
        addps   nb213nf_Vvdwtot(%rsp),%xmm3
        movaps  %xmm3,nb213nf_Vvdwtot(%rsp)

        ## do H1 interactions
        movaps  nb213nf_rinvH1(%rsp),%xmm7
        addps nb213nf_krsqH1(%rsp),%xmm7
        subps nb213nf_crf(%rsp),%xmm7   ## xmm7=rinv+ krsq-crf 
        mulps nb213nf_qqH(%rsp),%xmm7
        addps nb213nf_vctot(%rsp),%xmm7

        ## H2 interactions 
        movaps  nb213nf_rinvH2(%rsp),%xmm6
        addps nb213nf_krsqH2(%rsp),%xmm6
        subps nb213nf_crf(%rsp),%xmm6   ## xmm6=rinv+ krsq-crf 
        mulps nb213nf_qqH(%rsp),%xmm6
        addps %xmm7,%xmm6

        ## M interactions 
        movaps nb213nf_rinvM(%rsp),%xmm5
        addps nb213nf_krsqM(%rsp),%xmm5
        subps nb213nf_crf(%rsp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps nb213nf_qqM(%rsp),%xmm5
        addps %xmm6,%xmm5
        movaps %xmm5,nb213nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $4,nb213nf_innerk(%rsp)
        jl    _nb_kernel213nf_x86_64_sse.nb213nf_odd_inner
        jmp   _nb_kernel213nf_x86_64_sse.nb213nf_unroll_loop
_nb_kernel213nf_x86_64_sse.nb213nf_odd_inner: 
        addl $4,nb213nf_innerk(%rsp)
        jnz   _nb_kernel213nf_x86_64_sse.nb213nf_odd_loop
        jmp   _nb_kernel213nf_x86_64_sse.nb213nf_updateouterdata
_nb_kernel213nf_x86_64_sse.nb213nf_odd_loop: 

        movq  nb213nf_innerjjnr(%rsp),%rdx      ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb213nf_innerjjnr(%rsp)

        xorps %xmm4,%xmm4       ## clear reg.
        movss nb213nf_iqM(%rsp),%xmm4
        movq nb213nf_charge(%rbp),%rsi
        movhps nb213nf_iqH(%rsp),%xmm4    ## [qM  0  qH  qH] 
        shufps $41,%xmm4,%xmm4 ## [0 qH qH qM]

        movss (%rsi,%rax,4),%xmm3       ## charge in xmm3 
        shufps $0,%xmm3,%xmm3
        mulps %xmm4,%xmm3
        movaps %xmm3,nb213nf_qqM(%rsp)          ## use dummy qq for storage 

        xorps %xmm6,%xmm6
        movq nb213nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%ebx
        movq nb213nf_vdwparam(%rbp),%rsi
        shll %ebx
        addl nb213nf_ntia(%rsp),%ebx
        movlps (%rsi,%rbx,4),%xmm6
        movaps %xmm6,%xmm7
        shufps $252,%xmm6,%xmm6 ## 11111100
        shufps $253,%xmm7,%xmm7 ## 11111101
        movaps %xmm6,nb213nf_c6(%rsp)
        movaps %xmm7,nb213nf_c12(%rsp)

        movq nb213nf_pos(%rbp),%rsi
        lea (%rax,%rax,2),%rax

        movss nb213nf_ixO(%rsp),%xmm3
        movss nb213nf_iyO(%rsp),%xmm4
        movss nb213nf_izO(%rsp),%xmm5
        movss nb213nf_ixH1(%rsp),%xmm0
        movss nb213nf_iyH1(%rsp),%xmm1
        movss nb213nf_izH1(%rsp),%xmm2
        unpcklps nb213nf_ixH2(%rsp),%xmm3       ## ixO ixH2 - -
        unpcklps nb213nf_iyH2(%rsp),%xmm4       ## iyO iyH2 - -
        unpcklps nb213nf_izH2(%rsp),%xmm5       ## izO izH2 - -
        unpcklps nb213nf_ixM(%rsp),%xmm0        ## ixH1 ixM - -
        unpcklps nb213nf_iyM(%rsp),%xmm1        ## iyH1 iyM - -
        unpcklps nb213nf_izM(%rsp),%xmm2        ## izH1 izM - -
        unpcklps %xmm0,%xmm3    ## ixO ixH1 ixH2 ixM
        unpcklps %xmm1,%xmm4    ## same for y
        unpcklps %xmm2,%xmm5    ## same for z

        ## move j coords to xmm0-xmm2 
        movss (%rsi,%rax,4),%xmm0
        movss 4(%rsi,%rax,4),%xmm1
        movss 8(%rsi,%rax,4),%xmm2
        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2

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
        mulps nb213nf_krf(%rsp),%xmm0
        movaps %xmm0,nb213nf_krsqM(%rsp)

        rsqrtps %xmm4,%xmm5
        ## lookup seed in xmm5 
        movaps %xmm5,%xmm2
        mulps %xmm5,%xmm5
        movaps nb213nf_three(%rsp),%xmm1
        mulps %xmm4,%xmm5       ## rsq*lu*lu                    
        movaps nb213nf_half(%rsp),%xmm0
        subps %xmm5,%xmm1       ## 30-rsq*lu*lu 
        mulps %xmm2,%xmm1
        mulps %xmm1,%xmm0       ## xmm0=rinv

        movaps %xmm0,%xmm1      ## xmm1=rinv
        addps  nb213nf_krsqM(%rsp),%xmm1
        subps  nb213nf_crf(%rsp),%xmm1   ## xmm0=rinv+ krsq-crf 
        mulps  nb213nf_qqM(%rsp),%xmm1
        addps  nb213nf_vctot(%rsp),%xmm1
        movaps %xmm1,nb213nf_vctot(%rsp)

        movaps %xmm0,%xmm1
        mulps  %xmm1,%xmm1
        movaps %xmm1,%xmm2
        mulss  %xmm1,%xmm1
        mulss  %xmm2,%xmm1      ## xmm1=rinvsix
        xorps  %xmm4,%xmm4
        movss  %xmm1,%xmm4
        mulss  %xmm4,%xmm4      ## xmm4=rinvtwelve 
        mulss  nb213nf_c6(%rsp),%xmm1
        mulss  nb213nf_c12(%rsp),%xmm4
        movaps %xmm4,%xmm3
        subss  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addss  nb213nf_Vvdwtot(%rsp),%xmm3
        movss  %xmm3,nb213nf_Vvdwtot(%rsp)

        decl nb213nf_innerk(%rsp)
        jz    _nb_kernel213nf_x86_64_sse.nb213nf_updateouterdata
        jmp   _nb_kernel213nf_x86_64_sse.nb213nf_odd_loop
_nb_kernel213nf_x86_64_sse.nb213nf_updateouterdata: 
        ## get n from stack
        movl nb213nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb213nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb213nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb213nf_Vc(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## accumulate total lj energy and update it 
        movaps nb213nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movq  nb213nf_Vvdw(%rbp),%rax
        addss (%rax,%rdx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%rax,%rdx,4)

        ## finish if last 
        movl nb213nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel213nf_x86_64_sse.nb213nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb213nf_n(%rsp)
        jmp _nb_kernel213nf_x86_64_sse.nb213nf_outer
_nb_kernel213nf_x86_64_sse.nb213nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb213nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel213nf_x86_64_sse.nb213nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel213nf_x86_64_sse.nb213nf_threadloop
_nb_kernel213nf_x86_64_sse.nb213nf_end: 

        movl nb213nf_nouter(%rsp),%eax
        movl nb213nf_ninner(%rsp),%ebx
        movq nb213nf_outeriter(%rbp),%rcx
        movq nb213nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $600,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret




