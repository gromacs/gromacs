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



.globl nb_kernel113_x86_64_sse2
.globl _nb_kernel113_x86_64_sse2
nb_kernel113_x86_64_sse2:       
_nb_kernel113_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb113_fshift, 16
.set nb113_gid, 24
.set nb113_pos, 32
.set nb113_faction, 40
.set nb113_charge, 48
.set nb113_p_facel, 56
.set nb113_argkrf, 64
.set nb113_argcrf, 72
.set nb113_Vc, 80
.set nb113_type, 88
.set nb113_p_ntype, 96
.set nb113_vdwparam, 104
.set nb113_Vvdw, 112
.set nb113_p_tabscale, 120
.set nb113_VFtab, 128
.set nb113_invsqrta, 136
.set nb113_dvda, 144
.set nb113_p_gbtabscale, 152
.set nb113_GBtab, 160
.set nb113_p_nthreads, 168
.set nb113_count, 176
.set nb113_mtx, 184
.set nb113_outeriter, 192
.set nb113_inneriter, 200
.set nb113_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb113_ixO, 0
.set nb113_iyO, 16
.set nb113_izO, 32
.set nb113_ixH1, 48
.set nb113_iyH1, 64
.set nb113_izH1, 80
.set nb113_ixH2, 96
.set nb113_iyH2, 112
.set nb113_izH2, 128
.set nb113_ixM, 144
.set nb113_iyM, 160
.set nb113_izM, 176
.set nb113_iqH, 192
.set nb113_iqM, 208
.set nb113_dxO, 224
.set nb113_dyO, 240
.set nb113_dzO, 256
.set nb113_dxH1, 272
.set nb113_dyH1, 288
.set nb113_dzH1, 304
.set nb113_dxH2, 320
.set nb113_dyH2, 336
.set nb113_dzH2, 352
.set nb113_dxM, 368
.set nb113_dyM, 384
.set nb113_dzM, 400
.set nb113_qqH, 416
.set nb113_qqM, 432
.set nb113_c6, 448
.set nb113_c12, 464
.set nb113_six, 480
.set nb113_twelve, 496
.set nb113_vctot, 512
.set nb113_Vvdwtot, 528
.set nb113_fixO, 544
.set nb113_fiyO, 560
.set nb113_fizO, 576
.set nb113_fixH1, 592
.set nb113_fiyH1, 608
.set nb113_fizH1, 624
.set nb113_fixH2, 640
.set nb113_fiyH2, 656
.set nb113_fizH2, 672
.set nb113_fixM, 688
.set nb113_fiyM, 704
.set nb113_fizM, 720
.set nb113_fjx, 736
.set nb113_fjy, 752
.set nb113_fjz, 768
.set nb113_half, 784
.set nb113_three, 800
.set nb113_two, 816
.set nb113_rinvH1, 832
.set nb113_rinvH2, 848
.set nb113_rinvM, 864
.set nb113_nri, 880
.set nb113_iinr, 888
.set nb113_jindex, 896
.set nb113_jjnr, 904
.set nb113_shift, 912
.set nb113_shiftvec, 920
.set nb113_facel, 928
.set nb113_innerjjnr, 936
.set nb113_is3, 944
.set nb113_ii3, 948
.set nb113_ntia, 952
.set nb113_innerk, 956
.set nb113_n, 960
.set nb113_nn1, 964
.set nb113_nouter, 968
.set nb113_ninner, 972
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $984,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb113_nouter(%rsp)
        movl %eax,nb113_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb113_nri(%rsp)
        movq %rsi,nb113_iinr(%rsp)
        movq %rdx,nb113_jindex(%rsp)
        movq %rcx,nb113_jjnr(%rsp)
        movq %r8,nb113_shift(%rsp)
        movq %r9,nb113_shiftvec(%rsp)
        movq nb113_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb113_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx   ## upper half of half
        movl %eax,nb113_half(%rsp)
        movl %ebx,nb113_half+4(%rsp)
        movsd nb113_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## six
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## twelve
        movapd %xmm1,nb113_half(%rsp)
        movapd %xmm2,nb113_two(%rsp)
        movapd %xmm3,nb113_three(%rsp)
        movapd %xmm4,nb113_six(%rsp)
        movapd %xmm5,nb113_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb113_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb113_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb113_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb113_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb113_iqH(%rsp)
        movapd %xmm4,nb113_iqM(%rsp)

        movq  nb113_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb113_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb113_ntia(%rsp)
_nb_kernel113_x86_64_sse2.nb113_threadloop: 
        movq  nb113_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel113_x86_64_sse2.nb113_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel113_x86_64_sse2.nb113_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb113_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb113_n(%rsp)
        movl %ebx,nb113_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel113_x86_64_sse2.nb113_outerstart
        jmp _nb_kernel113_x86_64_sse2.nb113_end

_nb_kernel113_x86_64_sse2.nb113_outerstart: 
        ## ebx contains number of outer iterations
        addl nb113_nouter(%rsp),%ebx
        movl %ebx,nb113_nouter(%rsp)

_nb_kernel113_x86_64_sse2.nb113_outer: 
        movq  nb113_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb113_is3(%rsp)      ## store is3 

        movq  nb113_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb113_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb113_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb113_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3       ## ox
        addsd 8(%rax,%rbx,8),%xmm4      ## oy
        addsd 16(%rax,%rbx,8),%xmm5     ## oz   
        addsd 24(%rax,%rbx,8),%xmm6     ## h1x
        addsd 32(%rax,%rbx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb113_ixO(%rsp)
        movapd %xmm4,nb113_iyO(%rsp)
        movapd %xmm5,nb113_izO(%rsp)
        movapd %xmm6,nb113_ixH1(%rsp)
        movapd %xmm7,nb113_iyH1(%rsp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%rax,%rbx,8),%xmm6    ## h1z
        addsd 48(%rax,%rbx,8),%xmm0    ## h2x
        addsd 56(%rax,%rbx,8),%xmm1    ## h2y
        addsd 64(%rax,%rbx,8),%xmm2    ## h2z
        addsd 72(%rax,%rbx,8),%xmm3    ## mx
        addsd 80(%rax,%rbx,8),%xmm4    ## my
        addsd 88(%rax,%rbx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb113_izH1(%rsp)
        movapd %xmm0,nb113_ixH2(%rsp)
        movapd %xmm1,nb113_iyH2(%rsp)
        movapd %xmm2,nb113_izH2(%rsp)
        movapd %xmm3,nb113_ixM(%rsp)
        movapd %xmm4,nb113_iyM(%rsp)
        movapd %xmm5,nb113_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb113_vctot(%rsp)
        movapd %xmm4,nb113_Vvdwtot(%rsp)
        movapd %xmm4,nb113_fixO(%rsp)
        movapd %xmm4,nb113_fiyO(%rsp)
        movapd %xmm4,nb113_fizO(%rsp)
        movapd %xmm4,nb113_fixH1(%rsp)
        movapd %xmm4,nb113_fiyH1(%rsp)
        movapd %xmm4,nb113_fizH1(%rsp)
        movapd %xmm4,nb113_fixH2(%rsp)
        movapd %xmm4,nb113_fiyH2(%rsp)
        movapd %xmm4,nb113_fizH2(%rsp)
        movapd %xmm4,nb113_fixM(%rsp)
        movapd %xmm4,nb113_fiyM(%rsp)
        movapd %xmm4,nb113_fizM(%rsp)

        movq  nb113_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb113_pos(%rbp),%rsi
        movq  nb113_faction(%rbp),%rdi
        movq  nb113_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb113_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb113_ninner(%rsp),%ecx
        movl  %ecx,nb113_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb113_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel113_x86_64_sse2.nb113_unroll_loop
        jmp   _nb_kernel113_x86_64_sse2.nb113_checksingle
_nb_kernel113_x86_64_sse2.nb113_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb113_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb113_innerjjnr(%rsp)                   ## advance pointer (unrolled 2) 

        movq nb113_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb113_iqM(%rsp),%xmm3
        mulpd  nb113_iqH(%rsp),%xmm4

        movapd  %xmm3,nb113_qqM(%rsp)
        movapd  %xmm4,nb113_qqH(%rsp)

        movq nb113_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movl (%rsi,%rbx,4),%r9d
        movq nb113_vdwparam(%rbp),%rsi
        shll %r8d
        shll %r9d
        movl nb113_ntia(%rsp),%edi
        addl %edi,%r8d
        addl %edi,%r9d

        movlpd (%rsi,%r8,8),%xmm6       ## c6a c6b
        movhpd (%rsi,%r9,8),%xmm6       ## 
        movlpd 8(%rsi,%r8,8),%xmm7      ## c12a c12b 
        movhpd 8(%rsi,%r9,8),%xmm7      ## 

        movapd %xmm6,nb113_c6(%rsp)
        movapd %xmm7,nb113_c12(%rsp)

        movq nb113_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
    movlpd (%rsi,%rax,8),%xmm0
    movlpd 8(%rsi,%rax,8),%xmm1
    movlpd 16(%rsi,%rax,8),%xmm2
    movhpd (%rsi,%rbx,8),%xmm0
    movhpd 8(%rsi,%rbx,8),%xmm1
    movhpd 16(%rsi,%rbx,8),%xmm2

    ## xmm0 = jx
    ## xmm1 = jy
    ## xmm2 = jz

    ## O interaction
    ## copy to xmm3-xmm5
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5

    subpd nb113_ixO(%rsp),%xmm3
    subpd nb113_iyO(%rsp),%xmm4
    subpd nb113_izO(%rsp),%xmm5

    movapd %xmm3,%xmm13
    movapd %xmm4,%xmm14
    movapd %xmm5,%xmm15

        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5

        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3

    ## calc 1/rsq
    cvtpd2ps %xmm3,%xmm6
    rcpps %xmm6,%xmm6
    cvtps2pd %xmm6,%xmm6    ## lu in low xmm6 

    ## 1/x lookup seed in xmm6 
    movapd nb113_two(%rsp),%xmm4
    movapd %xmm3,%xmm5      ## rsq
    mulpd %xmm6,%xmm3       ## lu*rsq 
    subpd %xmm3,%xmm4       ## 2-lu*rsq 
    mulpd %xmm4,%xmm6       ## (new lu) 

    movapd nb113_two(%rsp),%xmm4
    mulpd %xmm6,%xmm5       ## lu*rsq 
    subpd %xmm5,%xmm4       ## 2-lu*rsq 
    mulpd %xmm6,%xmm4       ## xmm4=rinvsq 

    movapd %xmm4,%xmm3      ## rinvsq
    mulpd  %xmm4,%xmm4      ## rinv4
    mulpd  %xmm3,%xmm4      ## rinv6
    movapd %xmm4,%xmm5
    mulpd  %xmm5,%xmm5      ## rinv12
    mulpd  nb113_c6(%rsp),%xmm4
    mulpd  nb113_c12(%rsp),%xmm5
    movapd %xmm5,%xmm6
    subpd  %xmm4,%xmm6 ## Vvdw=vvdw12-vvdw6
    mulpd  nb113_six(%rsp),%xmm4
    mulpd  nb113_twelve(%rsp),%xmm5
    subpd  %xmm4,%xmm5
    mulpd  %xmm5,%xmm3  ## fscal

    addpd  nb113_Vvdwtot(%rsp),%xmm6
    movapd %xmm6,nb113_Vvdwtot(%rsp)

    mulpd  %xmm3,%xmm13 ## fx
    mulpd  %xmm3,%xmm14 ## fy
    mulpd  %xmm3,%xmm15 ## fz

    ## save j force temporarily
    movapd %xmm13,nb113_fjx(%rsp)
    movapd %xmm14,nb113_fjy(%rsp)
    movapd %xmm15,nb113_fjz(%rsp)

    ## increment i O force
    addpd nb113_fixO(%rsp),%xmm13
    addpd nb113_fiyO(%rsp),%xmm14
    addpd nb113_fizO(%rsp),%xmm15
    movapd %xmm13,nb113_fixO(%rsp)
    movapd %xmm14,nb113_fiyO(%rsp)
    movapd %xmm15,nb113_fizO(%rsp)
    ## finished O LJ interaction.


    ## do H1, H2, and M interactions in parallel.
    ## xmm0-xmm2 still contain j coordinates.        
    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb113_ixH1(%rsp),%xmm0
    subpd nb113_iyH1(%rsp),%xmm1
    subpd nb113_izH1(%rsp),%xmm2
    subpd nb113_ixH2(%rsp),%xmm3
    subpd nb113_iyH2(%rsp),%xmm4
    subpd nb113_izH2(%rsp),%xmm5
    subpd nb113_ixM(%rsp),%xmm6
    subpd nb113_iyM(%rsp),%xmm7
    subpd nb113_izM(%rsp),%xmm8

        movapd %xmm0,nb113_dxH1(%rsp)
        movapd %xmm1,nb113_dyH1(%rsp)
        movapd %xmm2,nb113_dzH1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb113_dxH2(%rsp)
        movapd %xmm4,nb113_dyH2(%rsp)
        movapd %xmm5,nb113_dzH2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb113_dxM(%rsp)
        movapd %xmm7,nb113_dyM(%rsp)
        movapd %xmm8,nb113_dzM(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for j atoms
    cvtpd2ps %xmm0,%xmm1
    cvtpd2ps %xmm3,%xmm4
    cvtpd2ps %xmm6,%xmm7
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm4,%xmm4
    rsqrtps %xmm7,%xmm7
    cvtps2pd %xmm1,%xmm1
    cvtps2pd %xmm4,%xmm4
    cvtps2pd %xmm7,%xmm7

        movapd  %xmm1,%xmm2
        movapd  %xmm4,%xmm5
    movapd  %xmm7,%xmm8

        mulpd   %xmm1,%xmm1 ## lu*lu
        mulpd   %xmm4,%xmm4 ## lu*lu
    mulpd   %xmm7,%xmm7 ## lu*lu

        movapd  nb113_three(%rsp),%xmm9
        movapd  %xmm9,%xmm10
    movapd  %xmm9,%xmm11

        mulpd   %xmm0,%xmm1 ## rsq*lu*lu
        mulpd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm7 ## rsq*lu*lu

        subpd   %xmm1,%xmm9
        subpd   %xmm4,%xmm10
    subpd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulpd   %xmm2,%xmm9
        mulpd   %xmm5,%xmm10
    mulpd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb113_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvM

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb113_three(%rsp),%xmm1
        movapd  %xmm1,%xmm4
    movapd  %xmm1,%xmm7

        mulpd   %xmm0,%xmm2 ## rsq*lu*lu
        mulpd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulpd   %xmm6,%xmm8 ## rsq*lu*lu

        subpd   %xmm2,%xmm1
        subpd   %xmm5,%xmm4
    subpd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulpd   %xmm1,%xmm9
        mulpd   %xmm4,%xmm10
    mulpd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movapd  nb113_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1
        mulpd   %xmm15,%xmm10 ##   rinvH2
    mulpd   %xmm15,%xmm11 ##   rinvM

        ## interactions 
    movapd %xmm9,%xmm0
    movapd %xmm10,%xmm1
    movapd %xmm11,%xmm2
    mulpd  %xmm9,%xmm9
    mulpd  %xmm10,%xmm10
    mulpd  %xmm11,%xmm11
    mulpd  nb113_qqH(%rsp),%xmm0
    mulpd  nb113_qqH(%rsp),%xmm1
    mulpd  nb113_qqM(%rsp),%xmm2
    mulpd  %xmm0,%xmm9
    mulpd  %xmm1,%xmm10
    mulpd  %xmm2,%xmm11

    addpd nb113_vctot(%rsp),%xmm0
    addpd %xmm2,%xmm1
    addpd %xmm1,%xmm0
    movapd %xmm0,nb113_vctot(%rsp)

    ## move j forces to xmm0-xmm2
    movq nb113_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm0
        movlpd 8(%rdi,%rax,8),%xmm1
        movlpd 16(%rdi,%rax,8),%xmm2
        movhpd (%rdi,%rbx,8),%xmm0
        movhpd 8(%rdi,%rbx,8),%xmm1
        movhpd 16(%rdi,%rbx,8),%xmm2

    movapd %xmm9,%xmm7
    movapd %xmm9,%xmm8
    movapd %xmm11,%xmm13
    movapd %xmm11,%xmm14
    movapd %xmm11,%xmm15
    movapd %xmm10,%xmm11
    movapd %xmm10,%xmm12

    ## add forces from O interaction
    addpd nb113_fjx(%rsp),%xmm0
    addpd nb113_fjy(%rsp),%xmm1
    addpd nb113_fjz(%rsp),%xmm2

        mulpd nb113_dxH1(%rsp),%xmm7
        mulpd nb113_dyH1(%rsp),%xmm8
        mulpd nb113_dzH1(%rsp),%xmm9
        mulpd nb113_dxH2(%rsp),%xmm10
        mulpd nb113_dyH2(%rsp),%xmm11
        mulpd nb113_dzH2(%rsp),%xmm12
        mulpd nb113_dxM(%rsp),%xmm13
        mulpd nb113_dyM(%rsp),%xmm14
        mulpd nb113_dzM(%rsp),%xmm15

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb113_fixH1(%rsp),%xmm7
    addpd nb113_fiyH1(%rsp),%xmm8
    addpd nb113_fizH1(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb113_fixH2(%rsp),%xmm10
    addpd nb113_fiyH2(%rsp),%xmm11
    addpd nb113_fizH2(%rsp),%xmm12

    addpd %xmm13,%xmm0
    addpd %xmm14,%xmm1
    addpd %xmm15,%xmm2
    addpd nb113_fixM(%rsp),%xmm13
    addpd nb113_fiyM(%rsp),%xmm14
    addpd nb113_fizM(%rsp),%xmm15

    movapd %xmm7,nb113_fixH1(%rsp)
    movapd %xmm8,nb113_fiyH1(%rsp)
    movapd %xmm9,nb113_fizH1(%rsp)
    movapd %xmm10,nb113_fixH2(%rsp)
    movapd %xmm11,nb113_fiyH2(%rsp)
    movapd %xmm12,nb113_fizH2(%rsp)
    movapd %xmm13,nb113_fixM(%rsp)
    movapd %xmm14,nb113_fiyM(%rsp)
    movapd %xmm15,nb113_fizM(%rsp)

    ## store back j forces from xmm0-xmm2
        movlpd %xmm0,(%rdi,%rax,8)
        movlpd %xmm1,8(%rdi,%rax,8)
        movlpd %xmm2,16(%rdi,%rax,8)
        movhpd %xmm0,(%rdi,%rbx,8)
        movhpd %xmm1,8(%rdi,%rbx,8)
        movhpd %xmm2,16(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb113_innerk(%rsp)
        jl   _nb_kernel113_x86_64_sse2.nb113_checksingle
        jmp  _nb_kernel113_x86_64_sse2.nb113_unroll_loop
_nb_kernel113_x86_64_sse2.nb113_checksingle: 
        movl  nb113_innerk(%rsp),%edx
        andl  $1,%edx
        jnz  _nb_kernel113_x86_64_sse2.nb113_dosingle
        jmp  _nb_kernel113_x86_64_sse2.nb113_updateouterdata
_nb_kernel113_x86_64_sse2.nb113_dosingle: 
        movq  nb113_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb113_innerjjnr(%rsp)

        movq nb113_charge(%rbp),%rsi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb113_iqM(%rsp),%xmm3
        mulpd  nb113_iqH(%rsp),%xmm4

        movapd  %xmm3,nb113_qqM(%rsp)
        movapd  %xmm4,nb113_qqH(%rsp)

        movq nb113_type(%rbp),%rsi
        movl (%rsi,%rax,4),%r8d
        movq nb113_vdwparam(%rbp),%rsi
        shll %r8d
        movl nb113_ntia(%rsp),%edi
        addl %edi,%r8d

        movsd (%rsi,%r8,8),%xmm6        ## c6a
    movsd 8(%rsi,%r8,8),%xmm7           ## c12a         
        movapd %xmm6,nb113_c6(%rsp)
        movapd %xmm7,nb113_c12(%rsp)

        movq nb113_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2  and xmm4-xmm6
        movlpd (%rsi,%rax,8),%xmm4
        movlpd 8(%rsi,%rax,8),%xmm5
        movlpd 16(%rsi,%rax,8),%xmm6
    movapd %xmm4,%xmm0
    movapd %xmm5,%xmm1
    movapd %xmm6,%xmm2

        ## calc dr 
        subsd nb113_ixO(%rsp),%xmm4
        subsd nb113_iyO(%rsp),%xmm5
        subsd nb113_izO(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb113_dxO(%rsp)
        movapd %xmm5,nb113_dyO(%rsp)
        movapd %xmm6,nb113_dzO(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move j coords to xmm4-xmm6 
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm5
        movapd %xmm2,%xmm6

        ## calc dr 
        subsd nb113_ixH1(%rsp),%xmm4
        subsd nb113_iyH1(%rsp),%xmm5
        subsd nb113_izH1(%rsp),%xmm6

        ## store dr 
        movapd %xmm4,nb113_dxH1(%rsp)
        movapd %xmm5,nb113_dyH1(%rsp)
        movapd %xmm6,nb113_dzH1(%rsp)
        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move j coords to xmm3-xmm5
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

        ## calc dr 
        subsd nb113_ixH2(%rsp),%xmm3
        subsd nb113_iyH2(%rsp),%xmm4
        subsd nb113_izH2(%rsp),%xmm5


        ## store dr 
        movapd %xmm3,nb113_dxH2(%rsp)
        movapd %xmm4,nb113_dyH2(%rsp)
        movapd %xmm5,nb113_dzH2(%rsp)
        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5

        ## move j coords to xmm4-xmm2
        movapd %xmm0,%xmm4
        movapd %xmm1,%xmm3
    ## xmm2 already contains z

        ## calc dr 
        subsd nb113_ixM(%rsp),%xmm4
        subsd nb113_iyM(%rsp),%xmm3
        subsd nb113_izM(%rsp),%xmm2

        ## store dr 
        movapd %xmm4,nb113_dxM(%rsp)
        movapd %xmm3,nb113_dyM(%rsp)
        movapd %xmm2,nb113_dzM(%rsp)

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113_three(%rsp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb113_three(%rsp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb113_rinvH1(%rsp)

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113_three(%rsp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb113_three(%rsp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb113_rinvH2(%rsp)

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113_three(%rsp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb113_three(%rsp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,nb113_rinvM(%rsp)

        ## do O interactions directly. xmm7=rsq
        cvtsd2ss %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movapd   nb113_two(%rsp),%xmm1
        movapd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm6
        subsd   %xmm6,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb113_c6(%rsp),%xmm1
        mulsd  nb113_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb113_Vvdwtot(%rsp),%xmm3
        mulsd  nb113_six(%rsp),%xmm1
        mulsd  nb113_twelve(%rsp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movapd %xmm2,%xmm4 ## total fsO 
        movsd %xmm3,nb113_Vvdwtot(%rsp)

        movapd nb113_dxO(%rsp),%xmm0
        movapd nb113_dyO(%rsp),%xmm1
        movapd nb113_dzO(%rsp),%xmm2
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update O forces 
        movapd nb113_fixO(%rsp),%xmm3
        movapd nb113_fiyO(%rsp),%xmm4
        movapd nb113_fizO(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb113_fixO(%rsp)
        movsd %xmm4,nb113_fiyO(%rsp)
        movsd %xmm7,nb113_fizO(%rsp)
        ## update j forces with water O 
        movsd %xmm0,nb113_fjx(%rsp)
        movsd %xmm1,nb113_fjy(%rsp)
        movsd %xmm2,nb113_fjz(%rsp)

        ## H1 interactions
        movapd  nb113_rinvH1(%rsp),%xmm6
        movapd  %xmm6,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm6=rinv, xmm4=rinvsq 
        mulsd  nb113_qqH(%rsp),%xmm6    ## xmm6=vcoul 
        mulsd  %xmm6,%xmm4              ## total fsH1 in xmm4 

        addsd  nb113_vctot(%rsp),%xmm6

        movapd nb113_dxH1(%rsp),%xmm0
        movapd nb113_dyH1(%rsp),%xmm1
        movapd nb113_dzH1(%rsp),%xmm2
        movsd %xmm6,nb113_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H1 forces 
        movapd nb113_fixH1(%rsp),%xmm3
        movapd nb113_fiyH1(%rsp),%xmm4
        movapd nb113_fizH1(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb113_fixH1(%rsp)
        movsd %xmm4,nb113_fiyH1(%rsp)
        movsd %xmm7,nb113_fizH1(%rsp)
        ## update j forces with water H1 
        addsd  nb113_fjx(%rsp),%xmm0
        addsd  nb113_fjy(%rsp),%xmm1
        addsd  nb113_fjz(%rsp),%xmm2
        movsd %xmm0,nb113_fjx(%rsp)
        movsd %xmm1,nb113_fjy(%rsp)
        movsd %xmm2,nb113_fjz(%rsp)

        ## H2 interactions 
        movapd  nb113_rinvH2(%rsp),%xmm5
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulsd  nb113_qqH(%rsp),%xmm5    ## xmm5=vcoul 
        mulsd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addsd  nb113_vctot(%rsp),%xmm5

        movapd nb113_dxH2(%rsp),%xmm0
        movapd nb113_dyH2(%rsp),%xmm1
        movapd nb113_dzH2(%rsp),%xmm2
        movsd %xmm5,nb113_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update H2 forces 
        movapd nb113_fixH2(%rsp),%xmm3
        movapd nb113_fiyH2(%rsp),%xmm4
        movapd nb113_fizH2(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb113_fixH2(%rsp)
        movsd %xmm4,nb113_fiyH2(%rsp)
        movsd %xmm7,nb113_fizH2(%rsp)
        ## update j forces with water H2 
        addsd  nb113_fjx(%rsp),%xmm0
        addsd  nb113_fjy(%rsp),%xmm1
        addsd  nb113_fjz(%rsp),%xmm2
        movsd %xmm0,nb113_fjx(%rsp)
        movsd %xmm1,nb113_fjy(%rsp)
        movsd %xmm2,nb113_fjz(%rsp)

        ## M interactions 
        movapd  nb113_rinvM(%rsp),%xmm5
        movapd  %xmm5,%xmm4
        mulsd   %xmm4,%xmm4     ## xmm5=rinv, xmm4=rinvsq 
        mulsd  nb113_qqM(%rsp),%xmm5    ## xmm5=vcoul 
        mulsd  %xmm5,%xmm4              ## total fsH1 in xmm4 

        addsd  nb113_vctot(%rsp),%xmm5

        movapd nb113_dxM(%rsp),%xmm0
        movapd nb113_dyM(%rsp),%xmm1
        movapd nb113_dzM(%rsp),%xmm2
        movsd %xmm5,nb113_vctot(%rsp)
        mulsd  %xmm4,%xmm0
        mulsd  %xmm4,%xmm1
        mulsd  %xmm4,%xmm2

        ## update M forces 
        movapd nb113_fixM(%rsp),%xmm3
        movapd nb113_fiyM(%rsp),%xmm4
        movapd nb113_fizM(%rsp),%xmm7
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm7
        movsd %xmm3,nb113_fixM(%rsp)
        movsd %xmm4,nb113_fiyM(%rsp)
        movsd %xmm7,nb113_fizM(%rsp)

        movq nb113_faction(%rbp),%rdi
        ## update j forces 
        addsd  nb113_fjx(%rsp),%xmm0
        addsd  nb113_fjy(%rsp),%xmm1
        addsd  nb113_fjz(%rsp),%xmm2
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        addsd %xmm0,%xmm3
        addsd %xmm1,%xmm4
        addsd %xmm2,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)

_nb_kernel113_x86_64_sse2.nb113_updateouterdata: 
        movl  nb113_ii3(%rsp),%ecx
        movq  nb113_faction(%rbp),%rdi
        movq  nb113_fshift(%rbp),%rsi
        movl  nb113_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb113_fixO(%rsp),%xmm0
        movapd nb113_fiyO(%rsp),%xmm1
        movapd nb113_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  (%rdi,%rcx,8),%xmm3
        movsd  8(%rdi,%rcx,8),%xmm4
        movsd  16(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,(%rdi,%rcx,8)
        movsd  %xmm4,8(%rdi,%rcx,8)
        movsd  %xmm5,16(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb113_fixH1(%rsp),%xmm0
        movapd nb113_fiyH1(%rsp),%xmm1
        movapd nb113_fizH1(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%rdi,%rcx,8),%xmm3
        movsd  32(%rdi,%rcx,8),%xmm4
        movsd  40(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,24(%rdi,%rcx,8)
        movsd  %xmm4,32(%rdi,%rcx,8)
        movsd  %xmm5,40(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb113_fixH2(%rsp),%xmm0
        movapd nb113_fiyH2(%rsp),%xmm1
        movapd nb113_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  48(%rdi,%rcx,8),%xmm3
        movsd  56(%rdi,%rcx,8),%xmm4
        movsd  64(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,48(%rdi,%rcx,8)
        movsd  %xmm4,56(%rdi,%rcx,8)
        movsd  %xmm5,64(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb113_fixM(%rsp),%xmm0
        movapd nb113_fiyM(%rsp),%xmm1
        movapd nb113_fizM(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  72(%rdi,%rcx,8),%xmm3
        movsd  80(%rdi,%rcx,8),%xmm4
        movsd  88(%rdi,%rcx,8),%xmm5
        subsd  %xmm0,%xmm3
        subsd  %xmm1,%xmm4
        subsd  %xmm2,%xmm5
        movsd  %xmm3,72(%rdi,%rcx,8)
        movsd  %xmm4,80(%rdi,%rcx,8)
        movsd  %xmm5,88(%rdi,%rcx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%rsi,%rdx,8),%xmm3
        movhpd 8(%rsi,%rdx,8),%xmm3
        movsd  16(%rsi,%rdx,8),%xmm4
        subpd  %xmm6,%xmm3
        subsd  %xmm7,%xmm4
        movlpd %xmm3,(%rsi,%rdx,8)
        movhpd %xmm3,8(%rsi,%rdx,8)
        movsd  %xmm4,16(%rsi,%rdx,8)

        ## get n from stack
        movl nb113_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb113_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb113_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb113_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb113_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb113_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb113_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel113_x86_64_sse2.nb113_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb113_n(%rsp)
        jmp _nb_kernel113_x86_64_sse2.nb113_outer
_nb_kernel113_x86_64_sse2.nb113_outerend: 
        ## check if more outer neighborlists remain
        movl  nb113_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel113_x86_64_sse2.nb113_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel113_x86_64_sse2.nb113_threadloop
_nb_kernel113_x86_64_sse2.nb113_end: 
        movl nb113_nouter(%rsp),%eax
        movl nb113_ninner(%rsp),%ebx
        movq nb113_outeriter(%rbp),%rcx
        movq nb113_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $984,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel113nf_x86_64_sse2
.globl _nb_kernel113nf_x86_64_sse2
nb_kernel113nf_x86_64_sse2:     
_nb_kernel113nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb113nf_fshift, 16
.set nb113nf_gid, 24
.set nb113nf_pos, 32
.set nb113nf_faction, 40
.set nb113nf_charge, 48
.set nb113nf_p_facel, 56
.set nb113nf_argkrf, 64
.set nb113nf_argcrf, 72
.set nb113nf_Vc, 80
.set nb113nf_type, 88
.set nb113nf_p_ntype, 96
.set nb113nf_vdwparam, 104
.set nb113nf_Vvdw, 112
.set nb113nf_p_tabscale, 120
.set nb113nf_VFtab, 128
.set nb113nf_invsqrta, 136
.set nb113nf_dvda, 144
.set nb113nf_p_gbtabscale, 152
.set nb113nf_GBtab, 160
.set nb113nf_p_nthreads, 168
.set nb113nf_count, 176
.set nb113nf_mtx, 184
.set nb113nf_outeriter, 192
.set nb113nf_inneriter, 200
.set nb113nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb113nf_ixO, 0
.set nb113nf_iyO, 16
.set nb113nf_izO, 32
.set nb113nf_ixH1, 48
.set nb113nf_iyH1, 64
.set nb113nf_izH1, 80
.set nb113nf_ixH2, 96
.set nb113nf_iyH2, 112
.set nb113nf_izH2, 128
.set nb113nf_ixM, 144
.set nb113nf_iyM, 160
.set nb113nf_izM, 176
.set nb113nf_iqH, 192
.set nb113nf_iqM, 208
.set nb113nf_qqH, 224
.set nb113nf_qqM, 240
.set nb113nf_c6, 256
.set nb113nf_c12, 272
.set nb113nf_vctot, 288
.set nb113nf_Vvdwtot, 304
.set nb113nf_half, 320
.set nb113nf_three, 336
.set nb113nf_two, 352
.set nb113nf_is3, 368
.set nb113nf_ii3, 372
.set nb113nf_nri, 376
.set nb113nf_iinr, 384
.set nb113nf_jindex, 392
.set nb113nf_jjnr, 400
.set nb113nf_shift, 408
.set nb113nf_shiftvec, 416
.set nb113nf_facel, 424
.set nb113nf_innerjjnr, 432
.set nb113nf_ntia, 440
.set nb113nf_innerk, 444
.set nb113nf_n, 448
.set nb113nf_nn1, 452
.set nb113nf_nouter, 456
.set nb113nf_ninner, 460
        push %rbp
        movq %rsp,%rbp
        push %rbx

        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $472,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb113nf_nouter(%rsp)
        movl %eax,nb113nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb113nf_nri(%rsp)
        movq %rsi,nb113nf_iinr(%rsp)
        movq %rdx,nb113nf_jindex(%rsp)
        movq %rcx,nb113nf_jjnr(%rsp)
        movq %r8,nb113nf_shift(%rsp)
        movq %r9,nb113nf_shiftvec(%rsp)
        movq nb113nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb113nf_facel(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb113nf_half(%rsp)
        movl %ebx,nb113nf_half+4(%rsp)
        movsd nb113nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb113nf_half(%rsp)
        movapd %xmm2,nb113nf_two(%rsp)
        movapd %xmm3,nb113nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb113nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb113nf_charge(%rbp),%rdx
        movsd 8(%rdx,%rbx,8),%xmm3
        movsd 24(%rdx,%rbx,8),%xmm4
        movq nb113nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb113nf_facel(%rsp),%xmm5
        mulsd  %xmm5,%xmm3
        mulsd  %xmm5,%xmm4

        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        movapd %xmm3,nb113nf_iqH(%rsp)
        movapd %xmm4,nb113nf_iqM(%rsp)

        movq  nb113nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movq nb113nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        movl  %ecx,nb113nf_ntia(%rsp)

_nb_kernel113nf_x86_64_sse2.nb113nf_threadloop: 
        movq  nb113nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel113nf_x86_64_sse2.nb113nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel113nf_x86_64_sse2.nb113nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb113nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb113nf_n(%rsp)
        movl %ebx,nb113nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel113nf_x86_64_sse2.nb113nf_outerstart
        jmp _nb_kernel113nf_x86_64_sse2.nb113nf_end

_nb_kernel113nf_x86_64_sse2.nb113nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb113nf_nouter(%rsp),%ebx
        movl %ebx,nb113nf_nouter(%rsp)

_nb_kernel113nf_x86_64_sse2.nb113nf_outer: 
        movq  nb113nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb113nf_is3(%rsp)            ## store is3 

        movq  nb113nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb113nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb113nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb113nf_ii3(%rsp)

        addsd (%rax,%rbx,8),%xmm3       ## ox
        addsd 8(%rax,%rbx,8),%xmm4      ## oy
        addsd 16(%rax,%rbx,8),%xmm5     ## oz   
        addsd 24(%rax,%rbx,8),%xmm6     ## h1x
        addsd 32(%rax,%rbx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb113nf_ixO(%rsp)
        movapd %xmm4,nb113nf_iyO(%rsp)
        movapd %xmm5,nb113nf_izO(%rsp)
        movapd %xmm6,nb113nf_ixH1(%rsp)
        movapd %xmm7,nb113nf_iyH1(%rsp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%rax,%rbx,8),%xmm6    ## h1z
        addsd 48(%rax,%rbx,8),%xmm0    ## h2x
        addsd 56(%rax,%rbx,8),%xmm1    ## h2y
        addsd 64(%rax,%rbx,8),%xmm2    ## h2z
        addsd 72(%rax,%rbx,8),%xmm3    ## mx
        addsd 80(%rax,%rbx,8),%xmm4    ## my
        addsd 88(%rax,%rbx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb113nf_izH1(%rsp)
        movapd %xmm0,nb113nf_ixH2(%rsp)
        movapd %xmm1,nb113nf_iyH2(%rsp)
        movapd %xmm2,nb113nf_izH2(%rsp)
        movapd %xmm3,nb113nf_ixM(%rsp)
        movapd %xmm4,nb113nf_iyM(%rsp)
        movapd %xmm5,nb113nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb113nf_vctot(%rsp)
        movapd %xmm4,nb113nf_Vvdwtot(%rsp)

        movq  nb113nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb113nf_pos(%rbp),%rsi
        movq  nb113nf_faction(%rbp),%rdi
        movq  nb113nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb113nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb113nf_ninner(%rsp),%ecx
        movl  %ecx,nb113nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb113nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel113nf_x86_64_sse2.nb113nf_unroll_loop
        jmp   _nb_kernel113nf_x86_64_sse2.nb113nf_checksingle
_nb_kernel113nf_x86_64_sse2.nb113nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb113nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb113nf_innerjjnr(%rsp)                 ## advance pointer (unrolled 2) 

        movq nb113nf_charge(%rbp),%rsi     ## base of charge[] 

        movlpd (%rsi,%rax,8),%xmm3
        movhpd (%rsi,%rbx,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb113nf_iqM(%rsp),%xmm3
        mulpd  nb113nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 
        movd  %ebx,%mm1

        movapd  %xmm3,nb113nf_qqM(%rsp)
        movapd  %xmm4,nb113nf_qqH(%rsp)

        movq nb113nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movl (%rsi,%rbx,4),%ebx
        movq nb113nf_vdwparam(%rbp),%rsi
        shll %eax
        shll %ebx
        movl nb113nf_ntia(%rsp),%edi
        addl %edi,%eax
        addl %edi,%ebx

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movlpd (%rsi,%rbx,8),%xmm7      ## c6b
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 
        movhpd 8(%rsi,%rbx,8),%xmm7     ## c6b c12b 

        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb113nf_c6(%rsp)
        movapd %xmm6,nb113nf_c12(%rsp)

        movq nb113nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move two coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2
        movhpd (%rsi,%rbx,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm1
        movhpd 16(%rsi,%rbx,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb113nf_ixO(%rsp),%xmm4
        movapd nb113nf_iyO(%rsp),%xmm5
        movapd nb113nf_izO(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb113nf_ixH1(%rsp),%xmm4
        movapd nb113nf_iyH1(%rsp),%xmm5
        movapd nb113nf_izH1(%rsp),%xmm6

        ## calc dr 
        subpd %xmm0,%xmm4
        subpd %xmm1,%xmm5
        subpd %xmm2,%xmm6

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm6
        addpd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb113nf_ixH2(%rsp),%xmm3
        movapd nb113nf_iyH2(%rsp),%xmm4
        movapd nb113nf_izH2(%rsp),%xmm5

        ## calc dr 
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5

        ## square it 
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        addpd %xmm4,%xmm5
        addpd %xmm3,%xmm5

        ## move ixM-izM to xmm2-xmm4  
        movapd nb113nf_iyM(%rsp),%xmm3
        movapd nb113nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb113nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtpd2ps %xmm6,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113nf_three(%rsp),%xmm1
        mulpd   %xmm6,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb113nf_three(%rsp),%xmm1
        subpd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,%xmm6     ## rinvH1

        ## rsqH2 - seed in xmm2 
        cvtpd2ps %xmm5,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113nf_three(%rsp),%xmm1
        mulpd   %xmm5,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb113nf_three(%rsp),%xmm1
        subpd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,%xmm5     ## rinvH2

        ## rsqM - seed in xmm2 
        cvtpd2ps %xmm4,%xmm2
        rsqrtps %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulpd   %xmm2,%xmm2
        movapd  nb113nf_three(%rsp),%xmm1
        mulpd   %xmm4,%xmm2     ## rsq*lu*lu 
        subpd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulpd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulpd   nb113nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulpd %xmm1,%xmm1       ## lu*lu 
        mulpd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb113nf_three(%rsp),%xmm1
        subpd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulpd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulpd nb113nf_half(%rsp),%xmm1   ## rinv 
        movapd  %xmm1,%xmm4     ## rinvM

        ## calculate coulomb potentials from rinv.
        addpd   %xmm5,%xmm6     ## rinvH1+rinvH2
        mulpd   nb113nf_qqM(%rsp),%xmm4
        mulpd   nb113nf_qqH(%rsp),%xmm6
        addpd   %xmm6,%xmm4
        addpd   nb113nf_vctot(%rsp),%xmm4
        movapd  %xmm4,nb113nf_vctot(%rsp)

        ## do O interactions - rsqO is in xmm7
        cvtpd2ps %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtps2pd %xmm2,%xmm2
        movapd   nb113nf_two(%rsp),%xmm1
        movapd   %xmm1,%xmm0
        mulpd   %xmm2,%xmm7
        subpd   %xmm7,%xmm1
        mulpd   %xmm1,%xmm2 ## iter1 
        mulpd   %xmm2,%xmm6
        subpd   %xmm6,%xmm0
        mulpd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb113nf_c6(%rsp),%xmm1
        mulpd  nb113nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb113nf_Vvdwtot(%rsp),%xmm3
        movapd %xmm3,nb113nf_Vvdwtot(%rsp)
        ## should we do one more iteration? 
        subl $2,nb113nf_innerk(%rsp)
        jl   _nb_kernel113nf_x86_64_sse2.nb113nf_checksingle
        jmp  _nb_kernel113nf_x86_64_sse2.nb113nf_unroll_loop
_nb_kernel113nf_x86_64_sse2.nb113nf_checksingle: 
        addl $2,nb113nf_innerk(%rsp)
        jnz  _nb_kernel113nf_x86_64_sse2.nb113nf_dosingle
        jmp  _nb_kernel113nf_x86_64_sse2.nb113nf_updateouterdata
_nb_kernel113nf_x86_64_sse2.nb113nf_dosingle: 
        movq  nb113nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        addq $4,nb113nf_innerjjnr(%rsp)

        movq nb113nf_charge(%rbp),%rsi     ## base of charge[] 

        xorpd %xmm3,%xmm3
        movlpd (%rsi,%rax,8),%xmm3
        movapd %xmm3,%xmm4
        mulpd  nb113nf_iqM(%rsp),%xmm3
        mulpd  nb113nf_iqH(%rsp),%xmm4

        movd  %eax,%mm0         ## use mmx registers as temp storage 

        movapd  %xmm3,nb113nf_qqM(%rsp)
        movapd  %xmm4,nb113nf_qqH(%rsp)

        movq nb113nf_type(%rbp),%rsi
        movl (%rsi,%rax,4),%eax
        movq nb113nf_vdwparam(%rbp),%rsi
        shll %eax
        movl nb113nf_ntia(%rsp),%edi
        addl %edi,%eax

        movlpd (%rsi,%rax,8),%xmm6      ## c6a
        movhpd 8(%rsi,%rax,8),%xmm6     ## c6a c12a 

        xorpd %xmm7,%xmm7
        movapd %xmm6,%xmm4
        unpcklpd %xmm7,%xmm4
        unpckhpd %xmm7,%xmm6

        movd  %mm0,%eax
        movd  %mm1,%ebx
        movapd %xmm4,nb113nf_c6(%rsp)
        movapd %xmm6,nb113nf_c12(%rsp)

        movq nb113nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move coordinates to xmm0-xmm2 
        movlpd (%rsi,%rax,8),%xmm0
        movlpd 8(%rsi,%rax,8),%xmm1
        movlpd 16(%rsi,%rax,8),%xmm2

        ## move ixO-izO to xmm4-xmm6 
        movapd nb113nf_ixO(%rsp),%xmm4
        movapd nb113nf_iyO(%rsp),%xmm5
        movapd nb113nf_izO(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4
        movapd %xmm4,%xmm7
        ## rsqO in xmm7 

        ## move ixH1-izH1 to xmm4-xmm6 
        movapd nb113nf_ixH1(%rsp),%xmm4
        movapd nb113nf_iyH1(%rsp),%xmm5
        movapd nb113nf_izH1(%rsp),%xmm6

        ## calc dr 
        subsd %xmm0,%xmm4
        subsd %xmm1,%xmm5
        subsd %xmm2,%xmm6

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm6
        addsd %xmm4,%xmm6
        ## rsqH1 in xmm6 

        ## move ixH2-izH2 to xmm3-xmm5  
        movapd nb113nf_ixH2(%rsp),%xmm3
        movapd nb113nf_iyH2(%rsp),%xmm4
        movapd nb113nf_izH2(%rsp),%xmm5

        ## calc dr 
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5

        ## square it 
        mulsd %xmm3,%xmm3
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        addsd %xmm4,%xmm5
        addsd %xmm3,%xmm5
        ## move ixM-izM to xmm2-xmm4  
        movapd nb113nf_iyM(%rsp),%xmm3
        movapd nb113nf_izM(%rsp),%xmm4
        subpd  %xmm1,%xmm3
        subpd  %xmm2,%xmm4
        movapd nb113nf_ixM(%rsp),%xmm2
        subpd  %xmm0,%xmm2

        ## square it 
        mulpd %xmm2,%xmm2
        mulpd %xmm3,%xmm3
        mulpd %xmm4,%xmm4
        addpd %xmm3,%xmm4
        addpd %xmm2,%xmm4
        ## rsqM in xmm4, rsqH2 in xmm5, rsqH1 in xmm6, rsqO in xmm7 

        ## start with rsqH1 - put seed in xmm2 
        cvtsd2ss %xmm6,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113nf_three(%rsp),%xmm1
        mulsd   %xmm6,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm6       ## rsq*lu*lu 
        movapd nb113nf_three(%rsp),%xmm1
        subsd %xmm6,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,%xmm6

        ## rsqH2 - seed in xmm2 
        cvtsd2ss %xmm5,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113nf_three(%rsp),%xmm1
        mulsd   %xmm5,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm5       ## rsq*lu*lu 
        movapd nb113nf_three(%rsp),%xmm1
        subsd %xmm5,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,%xmm5

        ## rsqM - seed in xmm2 
        cvtsd2ss %xmm4,%xmm2
        rsqrtss %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2

        movapd  %xmm2,%xmm3
        mulsd   %xmm2,%xmm2
        movapd  nb113nf_three(%rsp),%xmm1
        mulsd   %xmm4,%xmm2     ## rsq*lu*lu 
        subsd   %xmm2,%xmm1     ## 30-rsq*lu*lu 
        mulsd   %xmm3,%xmm1     ## lu*(3-rsq*lu*lu) 
        mulsd   nb113nf_half(%rsp),%xmm1   ## iter1 ( new lu) 

        movapd %xmm1,%xmm3
        mulsd %xmm1,%xmm1       ## lu*lu 
        mulsd %xmm1,%xmm4       ## rsq*lu*lu 
        movapd nb113nf_three(%rsp),%xmm1
        subsd %xmm4,%xmm1       ## 3-rsq*lu*lu 
        mulsd %xmm3,%xmm1       ## lu*( 3-rsq*lu*lu) 
        mulsd nb113nf_half(%rsp),%xmm1   ## rinv 
        movapd %xmm1,%xmm4

        ## Calculate coulomb potential
        addsd  %xmm5,%xmm6      ## rinvH1+rinvH2
        mulsd  nb113nf_qqM(%rsp),%xmm4
        mulsd  nb113nf_qqH(%rsp),%xmm6
        addsd  %xmm6,%xmm4
        addsd nb113nf_vctot(%rsp),%xmm4
        movsd %xmm4,nb113nf_vctot(%rsp)

        ## do O interactions directly. xmm7=rsq
        cvtsd2ss %xmm7,%xmm2
        movapd   %xmm7,%xmm6
        rcpps    %xmm2,%xmm2
        cvtss2sd %xmm2,%xmm2
        movapd   nb113nf_two(%rsp),%xmm1
        movapd   %xmm1,%xmm0
        mulsd   %xmm2,%xmm7
        subsd   %xmm7,%xmm1
        mulsd   %xmm1,%xmm2 ## iter1 
        mulsd   %xmm2,%xmm6
        subsd   %xmm6,%xmm0
        mulsd   %xmm2,%xmm0 ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb113nf_c6(%rsp),%xmm1
        mulsd  nb113nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb113nf_Vvdwtot(%rsp),%xmm3
        movsd %xmm3,nb113nf_Vvdwtot(%rsp)

_nb_kernel113nf_x86_64_sse2.nb113nf_updateouterdata: 
        ## get n from stack
        movl nb113nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb113nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb113nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb113nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb113nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb113nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

       ## finish if last 
        movl nb113nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel113nf_x86_64_sse2.nb113nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb113nf_n(%rsp)
        jmp _nb_kernel113nf_x86_64_sse2.nb113nf_outer
_nb_kernel113nf_x86_64_sse2.nb113nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb113nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel113nf_x86_64_sse2.nb113nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel113nf_x86_64_sse2.nb113nf_threadloop
_nb_kernel113nf_x86_64_sse2.nb113nf_end: 
        movl nb113nf_nouter(%rsp),%eax
        movl nb113nf_ninner(%rsp),%ebx
        movq nb113nf_outeriter(%rbp),%rcx
        movq nb113nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $472,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret


