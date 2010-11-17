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





.globl nb_kernel314_x86_64_sse2
.globl _nb_kernel314_x86_64_sse2
nb_kernel314_x86_64_sse2:       
_nb_kernel314_x86_64_sse2:      
##      Room for return address and rbp (16 bytes)
.set nb314_fshift, 16
.set nb314_gid, 24
.set nb314_pos, 32
.set nb314_faction, 40
.set nb314_charge, 48
.set nb314_p_facel, 56
.set nb314_argkrf, 64
.set nb314_argcrf, 72
.set nb314_Vc, 80
.set nb314_type, 88
.set nb314_p_ntype, 96
.set nb314_vdwparam, 104
.set nb314_Vvdw, 112
.set nb314_p_tabscale, 120
.set nb314_VFtab, 128
.set nb314_invsqrta, 136
.set nb314_dvda, 144
.set nb314_p_gbtabscale, 152
.set nb314_GBtab, 160
.set nb314_p_nthreads, 168
.set nb314_count, 176
.set nb314_mtx, 184
.set nb314_outeriter, 192
.set nb314_inneriter, 200
.set nb314_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb314_ixO, 0
.set nb314_iyO, 16
.set nb314_izO, 32
.set nb314_ixH1, 48
.set nb314_iyH1, 64
.set nb314_izH1, 80
.set nb314_ixH2, 96
.set nb314_iyH2, 112
.set nb314_izH2, 128
.set nb314_ixM, 144
.set nb314_iyM, 160
.set nb314_izM, 176
.set nb314_jxO, 192
.set nb314_jyO, 208
.set nb314_jzO, 224
.set nb314_jxH1, 240
.set nb314_jyH1, 256
.set nb314_jzH1, 272
.set nb314_jxH2, 288
.set nb314_jyH2, 304
.set nb314_jzH2, 320
.set nb314_jxM, 336
.set nb314_jyM, 352
.set nb314_jzM, 368
.set nb314_dxOO, 384
.set nb314_dyOO, 400
.set nb314_dzOO, 416
.set nb314_dxH1H1, 432
.set nb314_dyH1H1, 448
.set nb314_dzH1H1, 464
.set nb314_dxH1H2, 480
.set nb314_dyH1H2, 496
.set nb314_dzH1H2, 512
.set nb314_dxH1M, 528
.set nb314_dyH1M, 544
.set nb314_dzH1M, 560
.set nb314_dxH2H1, 576
.set nb314_dyH2H1, 592
.set nb314_dzH2H1, 608
.set nb314_dxH2H2, 624
.set nb314_dyH2H2, 640
.set nb314_dzH2H2, 656
.set nb314_dxH2M, 672
.set nb314_dyH2M, 688
.set nb314_dzH2M, 704
.set nb314_dxMH1, 720
.set nb314_dyMH1, 736
.set nb314_dzMH1, 752
.set nb314_dxMH2, 768
.set nb314_dyMH2, 784
.set nb314_dzMH2, 800
.set nb314_dxMM, 816
.set nb314_dyMM, 832
.set nb314_dzMM, 848
.set nb314_qqMM, 864
.set nb314_qqMH, 880
.set nb314_qqHH, 896
.set nb314_two, 912
.set nb314_tsc, 928
.set nb314_c6, 944
.set nb314_c12, 960
.set nb314_vctot, 976
.set nb314_Vvdwtot, 992
.set nb314_fixO, 1008
.set nb314_fiyO, 1024
.set nb314_fizO, 1040
.set nb314_fixH1, 1056
.set nb314_fiyH1, 1072
.set nb314_fizH1, 1088
.set nb314_fixH2, 1104
.set nb314_fiyH2, 1120
.set nb314_fizH2, 1136
.set nb314_fixM, 1152
.set nb314_fiyM, 1168
.set nb314_fizM, 1184
.set nb314_fjxO, 1200
.set nb314_fjyO, 1216
.set nb314_fjzO, 1232
.set nb314_fjxH1, 1248
.set nb314_fjyH1, 1264
.set nb314_fjzH1, 1280
.set nb314_fjxH2, 1296
.set nb314_fjyH2, 1312
.set nb314_fjzH2, 1328
.set nb314_epsH1, 1344
.set nb314_epsH2, 1360
.set nb314_epsM, 1376
.set nb314_half, 1392
.set nb314_three, 1408
.set nb314_six, 1424
.set nb314_twelve, 1440
.set nb314_rsqOO, 1456
.set nb314_rsqH1H1, 1472
.set nb314_rsqH1H2, 1488
.set nb314_rsqH1M, 1504
.set nb314_rsqH2H1, 1520
.set nb314_rsqH2H2, 1536
.set nb314_rsqH2M, 1552
.set nb314_rsqMH1, 1568
.set nb314_rsqMH2, 1584
.set nb314_rsqMM, 1600
.set nb314_rinvsqOO, 1616
.set nb314_rinvH1H1, 1632
.set nb314_rinvH1H2, 1648
.set nb314_rinvH1M, 1664
.set nb314_rinvH2H1, 1680
.set nb314_rinvH2H2, 1696
.set nb314_rinvH2M, 1712
.set nb314_rinvMH1, 1728
.set nb314_rinvMH2, 1744
.set nb314_rinvMM, 1760
.set nb314_is3, 1776
.set nb314_ii3, 1780
.set nb314_nri, 1784
.set nb314_iinr, 1792
.set nb314_jindex, 1800
.set nb314_jjnr, 1808
.set nb314_shift, 1816
.set nb314_shiftvec, 1824
.set nb314_facel, 1832
.set nb314_innerjjnr, 1840
.set nb314_innerk, 1848
.set nb314_n, 1852
.set nb314_nn1, 1856
.set nb314_nouter, 1860
.set nb314_ninner, 1864
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $1880,%rsp         ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb314_nouter(%rsp)
        movl %eax,nb314_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb314_nri(%rsp)
        movq %rsi,nb314_iinr(%rsp)
        movq %rdx,nb314_jindex(%rsp)
        movq %rcx,nb314_jjnr(%rsp)
        movq %r8,nb314_shift(%rsp)
        movq %r9,nb314_shiftvec(%rsp)
        movq nb314_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb314_facel(%rsp)

        movq nb314_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb314_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb314_half(%rsp)
        movl %ebx,nb314_half+4(%rsp)
        movsd nb314_half(%rsp),%xmm1
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
        movapd %xmm1,nb314_half(%rsp)
        movapd %xmm2,nb314_two(%rsp)
        movapd %xmm3,nb314_three(%rsp)
        movapd %xmm4,nb314_six(%rsp)
        movapd %xmm5,nb314_twelve(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb314_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb314_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5
        movq nb314_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb314_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb314_qqMM(%rsp)
        movapd %xmm4,nb314_qqMH(%rsp)
        movapd %xmm5,nb314_qqHH(%rsp)

        xorpd %xmm0,%xmm0
        movq  nb314_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb314_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb314_vdwparam(%rbp),%rax
        movlpd (%rax,%rdx,8),%xmm0
        movlpd 8(%rax,%rdx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb314_c6(%rsp)
        movapd %xmm1,nb314_c12(%rsp)

_nb_kernel314_x86_64_sse2.nb314_threadloop: 
        movq  nb314_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel314_x86_64_sse2.nb314_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel314_x86_64_sse2.nb314_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb314_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb314_n(%rsp)
        movl %ebx,nb314_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel314_x86_64_sse2.nb314_outerstart
        jmp _nb_kernel314_x86_64_sse2.nb314_end

_nb_kernel314_x86_64_sse2.nb314_outerstart: 
        ## ebx contains number of outer iterations
        addl nb314_nouter(%rsp),%ebx
        movl %ebx,nb314_nouter(%rsp)

_nb_kernel314_x86_64_sse2.nb314_outer: 
        movq  nb314_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb314_is3(%rsp)      ## store is3 

        movq  nb314_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb314_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]        
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb314_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb314_ii3(%rsp)

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
        movapd %xmm3,nb314_ixO(%rsp)
        movapd %xmm4,nb314_iyO(%rsp)
        movapd %xmm5,nb314_izO(%rsp)
        movapd %xmm6,nb314_ixH1(%rsp)
        movapd %xmm7,nb314_iyH1(%rsp)

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
        movapd %xmm6,nb314_izH1(%rsp)
        movapd %xmm0,nb314_ixH2(%rsp)
        movapd %xmm1,nb314_iyH2(%rsp)
        movapd %xmm2,nb314_izH2(%rsp)
        movapd %xmm3,nb314_ixM(%rsp)
        movapd %xmm4,nb314_iyM(%rsp)
        movapd %xmm5,nb314_izM(%rsp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb314_vctot(%rsp)
        movapd %xmm4,nb314_Vvdwtot(%rsp)
        movapd %xmm4,nb314_fixO(%rsp)
        movapd %xmm4,nb314_fiyO(%rsp)
        movapd %xmm4,nb314_fizO(%rsp)
        movapd %xmm4,nb314_fixH1(%rsp)
        movapd %xmm4,nb314_fiyH1(%rsp)
        movapd %xmm4,nb314_fizH1(%rsp)
        movapd %xmm4,nb314_fixH2(%rsp)
        movapd %xmm4,nb314_fiyH2(%rsp)
        movapd %xmm4,nb314_fizH2(%rsp)
        movapd %xmm4,nb314_fixM(%rsp)
        movapd %xmm4,nb314_fiyM(%rsp)
        movapd %xmm4,nb314_fizM(%rsp)

        movq  nb314_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb314_pos(%rbp),%rsi
        movq  nb314_faction(%rbp),%rdi
        movq  nb314_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb314_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb314_ninner(%rsp),%ecx
        movl  %ecx,nb314_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb314_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel314_x86_64_sse2.nb314_unroll_loop
        jmp   _nb_kernel314_x86_64_sse2.nb314_checksingle
_nb_kernel314_x86_64_sse2.nb314_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb314_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb314_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb314_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## load j O coordinates
    movlpd (%rsi,%rax,8),%xmm4
    movlpd 8(%rsi,%rax,8),%xmm5
    movlpd 16(%rsi,%rax,8),%xmm6
    movhpd (%rsi,%rbx,8),%xmm4
    movhpd 8(%rsi,%rbx,8),%xmm5
    movhpd 16(%rsi,%rbx,8),%xmm6

    ## xmm4 = Ox
    ## xmm5 = Oy
    ## xmm6 = Oz

    subpd nb314_ixO(%rsp),%xmm4
    subpd nb314_iyO(%rsp),%xmm5
    subpd nb314_izO(%rsp),%xmm6

    ## store dx/dy/dz
    movapd %xmm4,%xmm9
    movapd %xmm5,%xmm10
    movapd %xmm6,%xmm11

        ## square it 
        mulpd %xmm4,%xmm4
        mulpd %xmm5,%xmm5
        mulpd %xmm6,%xmm6
        addpd %xmm5,%xmm4
        addpd %xmm6,%xmm4       ## rsq in xmm4 

        cvtpd2ps %xmm4,%xmm6
        rcpps %xmm6,%xmm6
        cvtps2pd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb314_two(%rsp),%xmm0
        movapd %xmm4,%xmm5
        mulpd %xmm6,%xmm4       ## lu*rsq 
        subpd %xmm4,%xmm0       ## 2-lu*rsq 
        mulpd %xmm0,%xmm6       ## (new lu) 

        movapd nb314_two(%rsp),%xmm0
        mulpd %xmm6,%xmm5       ## lu*rsq 
        subpd %xmm5,%xmm0       ## 2-lu*rsq 
        mulpd %xmm6,%xmm0       ## xmm0=rinvsq 

        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1
        mulpd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulpd  nb314_c6(%rsp),%xmm1     ## mult by c6
        mulpd  nb314_c12(%rsp),%xmm2     ## mult by c12
        movapd %xmm2,%xmm5
        subpd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        mulpd  nb314_six(%rsp),%xmm1
        mulpd  nb314_twelve(%rsp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm2,%xmm0      ## xmm0=total fscal 

    ## increment potential
    addpd  nb314_Vvdwtot(%rsp),%xmm5
    movapd %xmm5,nb314_Vvdwtot(%rsp)

        movq   nb314_faction(%rbp),%rdi
        mulpd  %xmm0,%xmm9
        mulpd  %xmm0,%xmm10
        mulpd  %xmm0,%xmm11

    movapd nb314_fixO(%rsp),%xmm0
    movapd nb314_fiyO(%rsp),%xmm1
    movapd nb314_fizO(%rsp),%xmm2

    ## accumulate i forces
    addpd %xmm9,%xmm0
    addpd %xmm10,%xmm1
    addpd %xmm11,%xmm2
    movapd %xmm0,nb314_fixO(%rsp)
    movapd %xmm1,nb314_fiyO(%rsp)
    movapd %xmm2,nb314_fizO(%rsp)

        ## the fj's - start by accumulating forces from memory 
    movq nb314_faction(%rbp),%rdi
        movlpd (%rdi,%rax,8),%xmm3
        movlpd 8(%rdi,%rax,8),%xmm4
        movlpd 16(%rdi,%rax,8),%xmm5
        movhpd (%rdi,%rbx,8),%xmm3
        movhpd 8(%rdi,%rbx,8),%xmm4
        movhpd 16(%rdi,%rbx,8),%xmm5
        addpd %xmm9,%xmm3
        addpd %xmm10,%xmm4
        addpd %xmm11,%xmm5
        movlpd %xmm3,(%rdi,%rax,8)
        movlpd %xmm4,8(%rdi,%rax,8)
        movlpd %xmm5,16(%rdi,%rax,8)
        movhpd %xmm3,(%rdi,%rbx,8)
        movhpd %xmm4,8(%rdi,%rbx,8)
        movhpd %xmm5,16(%rdi,%rbx,8)
    ## done with OO interaction

    movq nb314_pos(%rbp),%rsi
    ## move j H1 coordinates to local temp variables 
    movlpd 24(%rsi,%rax,8),%xmm0
    movlpd 32(%rsi,%rax,8),%xmm1
    movlpd 40(%rsi,%rax,8),%xmm2
    movhpd 24(%rsi,%rbx,8),%xmm0
    movhpd 32(%rsi,%rbx,8),%xmm1
    movhpd 40(%rsi,%rbx,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb314_ixH1(%rsp),%xmm0
    subpd nb314_iyH1(%rsp),%xmm1
    subpd nb314_izH1(%rsp),%xmm2
    subpd nb314_ixH2(%rsp),%xmm3
    subpd nb314_iyH2(%rsp),%xmm4
    subpd nb314_izH2(%rsp),%xmm5
    subpd nb314_ixM(%rsp),%xmm6
    subpd nb314_iyM(%rsp),%xmm7
    subpd nb314_izM(%rsp),%xmm8

        movapd %xmm0,nb314_dxH1H1(%rsp)
        movapd %xmm1,nb314_dyH1H1(%rsp)
        movapd %xmm2,nb314_dzH1H1(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH2H1(%rsp)
        movapd %xmm4,nb314_dyH2H1(%rsp)
        movapd %xmm5,nb314_dzH2H1(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb314_dxMH1(%rsp)
        movapd %xmm7,nb314_dyMH1(%rsp)
        movapd %xmm8,nb314_dzMH1(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
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

        movapd  nb314_three(%rsp),%xmm9
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

        movapd  nb314_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1H1 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2H1
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMH1 

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb314_three(%rsp),%xmm1
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

        movapd  nb314_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H1 
        mulpd   %xmm15,%xmm10 ##   rinvH2H1
    mulpd   %xmm15,%xmm11 ##   rinvMH1

        movapd  %xmm9,nb314_rinvH1H1(%rsp)
        movapd  %xmm10,nb314_rinvH2H1(%rsp)
        movapd  %xmm11,nb314_rinvMH1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb314_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb314_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb314_epsH1(%rsp)
    movapd    %xmm3,nb314_epsH2(%rsp)
    movapd    %xmm6,nb314_epsM(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb314_epsH1(%rsp),%xmm12
    movapd nb314_epsH2(%rsp),%xmm13
    movapd nb314_epsM(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb314_qqHH(%rsp),%xmm12
    movapd nb314_qqMH(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm12,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm12,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb314_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb314_vctot(%rsp)

    movapd nb314_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb314_rinvH1H1(%rsp),%xmm4
    mulpd nb314_rinvH2H1(%rsp),%xmm8
    mulpd nb314_rinvMH1(%rsp),%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb314_faction(%rbp),%rdi
        movlpd 24(%rdi,%rax,8),%xmm0
        movlpd 32(%rdi,%rax,8),%xmm1
        movlpd 40(%rdi,%rax,8),%xmm2
        movhpd 24(%rdi,%rbx,8),%xmm0
        movhpd 32(%rdi,%rbx,8),%xmm1
        movhpd 40(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb314_dxH1H1(%rsp),%xmm3
        mulpd nb314_dyH1H1(%rsp),%xmm4
        mulpd nb314_dzH1H1(%rsp),%xmm5
        mulpd nb314_dxH2H1(%rsp),%xmm7
        mulpd nb314_dyH2H1(%rsp),%xmm8
        mulpd nb314_dzH2H1(%rsp),%xmm9
        mulpd nb314_dxMH1(%rsp),%xmm10
        mulpd nb314_dyMH1(%rsp),%xmm11
        mulpd nb314_dzMH1(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb314_fixH1(%rsp),%xmm3
    addpd nb314_fiyH1(%rsp),%xmm4
    addpd nb314_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb314_fixH2(%rsp),%xmm7
    addpd nb314_fiyH2(%rsp),%xmm8
    addpd nb314_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb314_fixM(%rsp),%xmm10
    addpd nb314_fiyM(%rsp),%xmm11
    addpd nb314_fizM(%rsp),%xmm12

    movapd %xmm3,nb314_fixH1(%rsp)
    movapd %xmm4,nb314_fiyH1(%rsp)
    movapd %xmm5,nb314_fizH1(%rsp)
    movapd %xmm7,nb314_fixH2(%rsp)
    movapd %xmm8,nb314_fiyH2(%rsp)
    movapd %xmm9,nb314_fizH2(%rsp)
    movapd %xmm10,nb314_fixM(%rsp)
    movapd %xmm11,nb314_fiyM(%rsp)
    movapd %xmm12,nb314_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movlpd %xmm0,24(%rdi,%rax,8)
        movlpd %xmm1,32(%rdi,%rax,8)
        movlpd %xmm2,40(%rdi,%rax,8)
        movhpd %xmm0,24(%rdi,%rbx,8)
        movhpd %xmm1,32(%rdi,%rbx,8)
        movhpd %xmm2,40(%rdi,%rbx,8)

        ## move j H2 coordinates to local temp variables 
    movq nb314_pos(%rbp),%rsi
    movlpd 48(%rsi,%rax,8),%xmm0
    movlpd 56(%rsi,%rax,8),%xmm1
    movlpd 64(%rsi,%rax,8),%xmm2
    movhpd 48(%rsi,%rbx,8),%xmm0
    movhpd 56(%rsi,%rbx,8),%xmm1
    movhpd 64(%rsi,%rbx,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb314_ixH1(%rsp),%xmm0
    subpd nb314_iyH1(%rsp),%xmm1
    subpd nb314_izH1(%rsp),%xmm2
    subpd nb314_ixH2(%rsp),%xmm3
    subpd nb314_iyH2(%rsp),%xmm4
    subpd nb314_izH2(%rsp),%xmm5
    subpd nb314_ixM(%rsp),%xmm6
    subpd nb314_iyM(%rsp),%xmm7
    subpd nb314_izM(%rsp),%xmm8

        movapd %xmm0,nb314_dxH1H2(%rsp)
        movapd %xmm1,nb314_dyH1H2(%rsp)
        movapd %xmm2,nb314_dzH1H2(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH2H2(%rsp)
        movapd %xmm4,nb314_dyH2H2(%rsp)
        movapd %xmm5,nb314_dzH2H2(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb314_dxMH2(%rsp)
        movapd %xmm7,nb314_dyMH2(%rsp)
        movapd %xmm8,nb314_dzMH2(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
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

        movapd  nb314_three(%rsp),%xmm9
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

        movapd  nb314_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1H2 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2H2
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMH2

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb314_three(%rsp),%xmm1
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

        movapd  nb314_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1H2
        mulpd   %xmm15,%xmm10 ##   rinvH2H2
    mulpd   %xmm15,%xmm11 ##   rinvMH2

        movapd  %xmm9,nb314_rinvH1H2(%rsp)
        movapd  %xmm10,nb314_rinvH2H2(%rsp)
        movapd  %xmm11,nb314_rinvMH2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb314_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb314_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb314_epsH1(%rsp)
    movapd    %xmm3,nb314_epsH2(%rsp)
    movapd    %xmm6,nb314_epsM(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb314_epsH1(%rsp),%xmm12
    movapd nb314_epsH2(%rsp),%xmm13
    movapd nb314_epsM(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb314_qqHH(%rsp),%xmm12
    movapd nb314_qqMH(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm12,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm12,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb314_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb314_vctot(%rsp)

    movapd nb314_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb314_rinvH1H2(%rsp),%xmm4
    mulpd nb314_rinvH2H2(%rsp),%xmm8
    mulpd nb314_rinvMH2(%rsp),%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb314_faction(%rbp),%rdi
        movlpd 48(%rdi,%rax,8),%xmm0
        movlpd 56(%rdi,%rax,8),%xmm1
        movlpd 64(%rdi,%rax,8),%xmm2
        movhpd 48(%rdi,%rbx,8),%xmm0
        movhpd 56(%rdi,%rbx,8),%xmm1
        movhpd 64(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb314_dxH1H2(%rsp),%xmm3
        mulpd nb314_dyH1H2(%rsp),%xmm4
        mulpd nb314_dzH1H2(%rsp),%xmm5
        mulpd nb314_dxH2H2(%rsp),%xmm7
        mulpd nb314_dyH2H2(%rsp),%xmm8
        mulpd nb314_dzH2H2(%rsp),%xmm9
        mulpd nb314_dxMH2(%rsp),%xmm10
        mulpd nb314_dyMH2(%rsp),%xmm11
        mulpd nb314_dzMH2(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb314_fixH1(%rsp),%xmm3
    addpd nb314_fiyH1(%rsp),%xmm4
    addpd nb314_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb314_fixH2(%rsp),%xmm7
    addpd nb314_fiyH2(%rsp),%xmm8
    addpd nb314_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb314_fixM(%rsp),%xmm10
    addpd nb314_fiyM(%rsp),%xmm11
    addpd nb314_fizM(%rsp),%xmm12

    movapd %xmm3,nb314_fixH1(%rsp)
    movapd %xmm4,nb314_fiyH1(%rsp)
    movapd %xmm5,nb314_fizH1(%rsp)
    movapd %xmm7,nb314_fixH2(%rsp)
    movapd %xmm8,nb314_fiyH2(%rsp)
    movapd %xmm9,nb314_fizH2(%rsp)
    movapd %xmm10,nb314_fixM(%rsp)
    movapd %xmm11,nb314_fiyM(%rsp)
    movapd %xmm12,nb314_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movlpd %xmm0,48(%rdi,%rax,8)
        movlpd %xmm1,56(%rdi,%rax,8)
        movlpd %xmm2,64(%rdi,%rax,8)
        movhpd %xmm0,48(%rdi,%rbx,8)
        movhpd %xmm1,56(%rdi,%rbx,8)
        movhpd %xmm2,64(%rdi,%rbx,8)

        ## move j M coordinates to local temp variables 
    movq nb314_pos(%rbp),%rsi
    movlpd 72(%rsi,%rax,8),%xmm0
    movlpd 80(%rsi,%rax,8),%xmm1
    movlpd 88(%rsi,%rax,8),%xmm2
    movhpd 72(%rsi,%rbx,8),%xmm0
    movhpd 80(%rsi,%rbx,8),%xmm1
    movhpd 88(%rsi,%rbx,8),%xmm2

    ## xmm0 = Mx
    ## xmm1 = My
    ## xmm2 = Mz

    movapd %xmm0,%xmm3
    movapd %xmm1,%xmm4
    movapd %xmm2,%xmm5
    movapd %xmm0,%xmm6
    movapd %xmm1,%xmm7
    movapd %xmm2,%xmm8

    subpd nb314_ixH1(%rsp),%xmm0
    subpd nb314_iyH1(%rsp),%xmm1
    subpd nb314_izH1(%rsp),%xmm2
    subpd nb314_ixH2(%rsp),%xmm3
    subpd nb314_iyH2(%rsp),%xmm4
    subpd nb314_izH2(%rsp),%xmm5
    subpd nb314_ixM(%rsp),%xmm6
    subpd nb314_iyM(%rsp),%xmm7
    subpd nb314_izM(%rsp),%xmm8

        movapd %xmm0,nb314_dxH1M(%rsp)
        movapd %xmm1,nb314_dyH1M(%rsp)
        movapd %xmm2,nb314_dzH1M(%rsp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH2M(%rsp)
        movapd %xmm4,nb314_dyH2M(%rsp)
        movapd %xmm5,nb314_dzH2M(%rsp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        movapd %xmm6,nb314_dxMM(%rsp)
        movapd %xmm7,nb314_dyMM(%rsp)
        movapd %xmm8,nb314_dzMM(%rsp)
        mulpd  %xmm6,%xmm6
        mulpd  %xmm7,%xmm7
        mulpd  %xmm8,%xmm8
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
    addpd  %xmm7,%xmm6
    addpd  %xmm8,%xmm6

        ## start doing invsqrt for jM atoms
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

        movapd  nb314_three(%rsp),%xmm9
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

        movapd  nb314_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ## first iteration for rinvH1M 
        mulpd   %xmm15,%xmm10 ## first iteration for rinvH2M
    mulpd   %xmm15,%xmm11 ## first iteration for rinvMM

    ## second iteration step    
        movapd  %xmm9,%xmm2
        movapd  %xmm10,%xmm5
    movapd  %xmm11,%xmm8

        mulpd   %xmm2,%xmm2 ## lu*lu
        mulpd   %xmm5,%xmm5 ## lu*lu
    mulpd   %xmm8,%xmm8 ## lu*lu

        movapd  nb314_three(%rsp),%xmm1
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

        movapd  nb314_half(%rsp),%xmm15
        mulpd   %xmm15,%xmm9 ##  rinvH1M
        mulpd   %xmm15,%xmm10 ##   rinvH2M
    mulpd   %xmm15,%xmm11 ##   rinvMM


        movapd  %xmm9,nb314_rinvH1M(%rsp)
        movapd  %xmm10,nb314_rinvH2M(%rsp)
        movapd  %xmm11,nb314_rinvMM(%rsp)

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movapd nb314_tsc(%rsp),%xmm1
    mulpd  %xmm9,%xmm0 ## r
    mulpd  %xmm10,%xmm3
    mulpd  %xmm11,%xmm6
    mulpd  %xmm1,%xmm0 ## rtab
    mulpd  %xmm1,%xmm3
    mulpd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttpd2dq %xmm0,%xmm1
    cvttpd2dq %xmm3,%xmm4
    cvttpd2dq %xmm6,%xmm7

    ## convert back to float
    cvtdq2pd  %xmm1,%xmm2
    cvtdq2pd  %xmm4,%xmm5
    cvtdq2pd  %xmm7,%xmm8

    ## multiply by 4
    pslld   $2,%xmm1
    pslld   $2,%xmm4
    pslld   $2,%xmm7

    ## move to integer registers
    pshufd $1,%xmm1,%xmm13
    pshufd $1,%xmm4,%xmm14
    pshufd $1,%xmm7,%xmm15
    movd    %xmm1,%r8d
    movd    %xmm4,%r10d
    movd    %xmm7,%r12d
    movd    %xmm13,%r9d
    movd    %xmm14,%r11d
    movd    %xmm15,%r13d

    movq nb314_VFtab(%rbp),%rsi

    ## calculate eps
    subpd     %xmm2,%xmm0
    subpd     %xmm5,%xmm3
    subpd     %xmm8,%xmm6

    movapd    %xmm0,nb314_epsH1(%rsp)
    movapd    %xmm3,nb314_epsH2(%rsp)
    movapd    %xmm6,nb314_epsM(%rsp)

    ## Load LOTS of table data
    movlpd (%rsi,%r8,8),%xmm0
    movlpd 8(%rsi,%r8,8),%xmm1
    movlpd 16(%rsi,%r8,8),%xmm2
    movlpd 24(%rsi,%r8,8),%xmm3
    movlpd (%rsi,%r10,8),%xmm4
    movlpd 8(%rsi,%r10,8),%xmm5
    movlpd 16(%rsi,%r10,8),%xmm6
    movlpd 24(%rsi,%r10,8),%xmm7
    movlpd (%rsi,%r12,8),%xmm8
    movlpd 8(%rsi,%r12,8),%xmm9
    movlpd 16(%rsi,%r12,8),%xmm10
    movlpd 24(%rsi,%r12,8),%xmm11
    movhpd (%rsi,%r9,8),%xmm0
    movhpd 8(%rsi,%r9,8),%xmm1
    movhpd 16(%rsi,%r9,8),%xmm2
    movhpd 24(%rsi,%r9,8),%xmm3
    movhpd (%rsi,%r11,8),%xmm4
    movhpd 8(%rsi,%r11,8),%xmm5
    movhpd 16(%rsi,%r11,8),%xmm6
    movhpd 24(%rsi,%r11,8),%xmm7
    movhpd (%rsi,%r13,8),%xmm8
    movhpd 8(%rsi,%r13,8),%xmm9
    movhpd 16(%rsi,%r13,8),%xmm10
    movhpd 24(%rsi,%r13,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movapd nb314_epsH1(%rsp),%xmm12
    movapd nb314_epsH2(%rsp),%xmm13
    movapd nb314_epsM(%rsp),%xmm14

    mulpd  %xmm12,%xmm3  ## Heps
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11
    mulpd  %xmm12,%xmm2  ## Geps
    mulpd  %xmm13,%xmm6
    mulpd  %xmm14,%xmm10
    mulpd  %xmm12,%xmm3  ## Heps2
    mulpd  %xmm13,%xmm7
    mulpd  %xmm14,%xmm11

    addpd  %xmm2,%xmm1  ## F+Geps
    addpd  %xmm6,%xmm5
    addpd  %xmm10,%xmm9
    addpd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addpd  %xmm7,%xmm5
    addpd  %xmm11,%xmm9
    addpd  %xmm3,%xmm3   ## 2*Heps2
    addpd  %xmm7,%xmm7
    addpd  %xmm11,%xmm11
    addpd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addpd  %xmm6,%xmm7
    addpd  %xmm10,%xmm11
    addpd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addpd  %xmm5,%xmm7
    addpd  %xmm9,%xmm11
    mulpd  %xmm12,%xmm1  ## eps*Fp
    mulpd  %xmm13,%xmm5
    mulpd  %xmm14,%xmm9
    movapd nb314_qqMH(%rsp),%xmm12
    movapd nb314_qqMM(%rsp),%xmm13
    addpd  %xmm0,%xmm1    ## VV
    addpd  %xmm4,%xmm5
    addpd  %xmm8,%xmm9
    mulpd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulpd  %xmm12,%xmm5
    mulpd  %xmm13,%xmm9
    mulpd  %xmm12,%xmm3   ## FF*qq = fij
    mulpd  %xmm12,%xmm7
    mulpd  %xmm13,%xmm11

    ## accumulate vctot
    addpd  nb314_vctot(%rsp),%xmm1
    addpd  %xmm9,%xmm5
    addpd  %xmm5,%xmm1
    movapd %xmm1,nb314_vctot(%rsp)

    movapd nb314_tsc(%rsp),%xmm10
    mulpd  %xmm10,%xmm3 ## fscal
    mulpd  %xmm10,%xmm7
    mulpd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subpd %xmm3,%xmm4
    subpd %xmm7,%xmm8
    subpd %xmm10,%xmm11

    mulpd nb314_rinvH1M(%rsp),%xmm4
    mulpd nb314_rinvH2M(%rsp),%xmm8
    mulpd nb314_rinvMM(%rsp),%xmm11

    ## move j M forces to xmm0-xmm2
    movq nb314_faction(%rbp),%rdi
        movlpd 72(%rdi,%rax,8),%xmm0
        movlpd 80(%rdi,%rax,8),%xmm1
        movlpd 88(%rdi,%rax,8),%xmm2
        movhpd 72(%rdi,%rbx,8),%xmm0
        movhpd 80(%rdi,%rbx,8),%xmm1
        movhpd 88(%rdi,%rbx,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulpd nb314_dxH1M(%rsp),%xmm3
        mulpd nb314_dyH1M(%rsp),%xmm4
        mulpd nb314_dzH1M(%rsp),%xmm5
        mulpd nb314_dxH2M(%rsp),%xmm7
        mulpd nb314_dyH2M(%rsp),%xmm8
        mulpd nb314_dzH2M(%rsp),%xmm9
        mulpd nb314_dxMM(%rsp),%xmm10
        mulpd nb314_dyMM(%rsp),%xmm11
        mulpd nb314_dzMM(%rsp),%xmm12

    addpd %xmm3,%xmm0
    addpd %xmm4,%xmm1
    addpd %xmm5,%xmm2
    addpd nb314_fixH1(%rsp),%xmm3
    addpd nb314_fiyH1(%rsp),%xmm4
    addpd nb314_fizH1(%rsp),%xmm5

    addpd %xmm7,%xmm0
    addpd %xmm8,%xmm1
    addpd %xmm9,%xmm2
    addpd nb314_fixH2(%rsp),%xmm7
    addpd nb314_fiyH2(%rsp),%xmm8
    addpd nb314_fizH2(%rsp),%xmm9

    addpd %xmm10,%xmm0
    addpd %xmm11,%xmm1
    addpd %xmm12,%xmm2
    addpd nb314_fixM(%rsp),%xmm10
    addpd nb314_fiyM(%rsp),%xmm11
    addpd nb314_fizM(%rsp),%xmm12

    movapd %xmm3,nb314_fixH1(%rsp)
    movapd %xmm4,nb314_fiyH1(%rsp)
    movapd %xmm5,nb314_fizH1(%rsp)
    movapd %xmm7,nb314_fixH2(%rsp)
    movapd %xmm8,nb314_fiyH2(%rsp)
    movapd %xmm9,nb314_fizH2(%rsp)
    movapd %xmm10,nb314_fixM(%rsp)
    movapd %xmm11,nb314_fiyM(%rsp)
    movapd %xmm12,nb314_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movlpd %xmm0,72(%rdi,%rax,8)
        movlpd %xmm1,80(%rdi,%rax,8)
        movlpd %xmm2,88(%rdi,%rax,8)
        movhpd %xmm0,72(%rdi,%rbx,8)
        movhpd %xmm1,80(%rdi,%rbx,8)
        movhpd %xmm2,88(%rdi,%rbx,8)

        ## should we do one more iteration? 
        subl $2,nb314_innerk(%rsp)
        jl    _nb_kernel314_x86_64_sse2.nb314_checksingle
        jmp   _nb_kernel314_x86_64_sse2.nb314_unroll_loop
_nb_kernel314_x86_64_sse2.nb314_checksingle: 
        movl  nb314_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel314_x86_64_sse2.nb314_dosingle
        jmp   _nb_kernel314_x86_64_sse2.nb314_updateouterdata
_nb_kernel314_x86_64_sse2.nb314_dosingle: 
        movq  nb314_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb314_pos(%rbp),%rsi
        lea  (%rax,%rax,2),%rax

        ## load j O coordinates
    movsd (%rsi,%rax,8),%xmm4
    movsd 8(%rsi,%rax,8),%xmm5
    movsd 16(%rsi,%rax,8),%xmm6

    ## xmm4 = Ox
    ## xmm5 = Oy
    ## xmm6 = Oz

    subsd nb314_ixO(%rsp),%xmm4
    subsd nb314_iyO(%rsp),%xmm5
    subsd nb314_izO(%rsp),%xmm6

    ## store dx/dy/dz
    movapd %xmm4,%xmm9
    movapd %xmm5,%xmm10
    movapd %xmm6,%xmm11

        ## square it 
        mulsd %xmm4,%xmm4
        mulsd %xmm5,%xmm5
        mulsd %xmm6,%xmm6
        addsd %xmm5,%xmm4
        addsd %xmm6,%xmm4       ## rsq in xmm4 

        cvtsd2ss %xmm4,%xmm6
        rcpss %xmm6,%xmm6
        cvtss2sd %xmm6,%xmm6    ## lu in low xmm6 

        ## 1/x lookup seed in xmm6 
        movapd nb314_two(%rsp),%xmm0
        movapd %xmm4,%xmm5
        mulsd %xmm6,%xmm4       ## lu*rsq 
        subsd %xmm4,%xmm0       ## 2-lu*rsq 
        mulsd %xmm0,%xmm6       ## (new lu) 

        movapd nb314_two(%rsp),%xmm0
        mulsd %xmm6,%xmm5       ## lu*rsq 
        subsd %xmm5,%xmm0       ## 2-lu*rsq 
        mulsd %xmm6,%xmm0       ## xmm0=rinvsq 

        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 

        mulsd  nb314_c6(%rsp),%xmm1     ## mult by c6
        mulsd  nb314_c12(%rsp),%xmm2     ## mult by c12
        movapd %xmm2,%xmm5
        subsd  %xmm1,%xmm5      ## Vvdw=Vvdw12-Vvdw6 
        mulsd  nb314_six(%rsp),%xmm1
        mulsd  nb314_twelve(%rsp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm2,%xmm0      ## xmm0=total fscal 

    ## increment potential
    addsd  nb314_Vvdwtot(%rsp),%xmm5
    movsd %xmm5,nb314_Vvdwtot(%rsp)

        movq   nb314_faction(%rbp),%rdi
        mulsd  %xmm0,%xmm9
        mulsd  %xmm0,%xmm10
        mulsd  %xmm0,%xmm11

    movapd nb314_fixO(%rsp),%xmm0
    movapd nb314_fiyO(%rsp),%xmm1
    movapd nb314_fizO(%rsp),%xmm2

    ## accumulate i forces
    addsd %xmm9,%xmm0
    addsd %xmm10,%xmm1
    addsd %xmm11,%xmm2
    movsd %xmm0,nb314_fixO(%rsp)
    movsd %xmm1,nb314_fiyO(%rsp)
    movsd %xmm2,nb314_fizO(%rsp)

        ## the fj's - start by accumulating forces from memory 
    movq nb314_faction(%rbp),%rdi
        addsd (%rdi,%rax,8),%xmm9
        addsd 8(%rdi,%rax,8),%xmm10
        addsd 16(%rdi,%rax,8),%xmm11
        movsd %xmm9,(%rdi,%rax,8)
        movsd %xmm10,8(%rdi,%rax,8)
        movsd %xmm11,16(%rdi,%rax,8)
    ## done with OO interaction

        ## move j H1 coordinates to local temp variables 
    movq nb314_pos(%rbp),%rsi
    movsd 24(%rsi,%rax,8),%xmm0
    movsd 32(%rsi,%rax,8),%xmm1
    movsd 40(%rsi,%rax,8),%xmm2

    ## xmm0 = H1x
    ## xmm1 = H1y
    ## xmm2 = H1z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb314_ixH1(%rsp),%xmm0
    subsd nb314_iyH1(%rsp),%xmm1
    subsd nb314_izH1(%rsp),%xmm2
    subsd nb314_ixH2(%rsp),%xmm3
    subsd nb314_iyH2(%rsp),%xmm4
    subsd nb314_izH2(%rsp),%xmm5
    subsd nb314_ixM(%rsp),%xmm6
    subsd nb314_iyM(%rsp),%xmm7
    subsd nb314_izM(%rsp),%xmm8

        movsd %xmm0,nb314_dxH1H1(%rsp)
        movsd %xmm1,nb314_dyH1H1(%rsp)
        movsd %xmm2,nb314_dzH1H1(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb314_dxH2H1(%rsp)
        movsd %xmm4,nb314_dyH2H1(%rsp)
        movsd %xmm5,nb314_dzH2H1(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb314_dxMH1(%rsp)
        movsd %xmm7,nb314_dyMH1(%rsp)
        movsd %xmm8,nb314_dzMH1(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH1 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb314_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb314_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1H1 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2H1
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMH1 

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb314_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb314_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H1 
        mulsd   %xmm15,%xmm10 ##   rinvH2H1
    mulsd   %xmm15,%xmm11 ##   rinvMH1

        movsd  %xmm9,nb314_rinvH1H1(%rsp)
        movsd  %xmm10,nb314_rinvH2H1(%rsp)
        movsd  %xmm11,nb314_rinvMH1(%rsp)

        ## H1 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb314_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shlq  $2,%r8
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb314_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb314_epsH1(%rsp)
    movsd    %xmm3,nb314_epsH2(%rsp)
    movsd    %xmm6,nb314_epsM(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb314_epsH1(%rsp),%xmm12
    movsd nb314_epsH2(%rsp),%xmm13
    movsd nb314_epsM(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb314_qqHH(%rsp),%xmm12
    movsd nb314_qqMH(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm12,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm12,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb314_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb314_vctot(%rsp)

    movsd nb314_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb314_rinvH1H1(%rsp),%xmm4
    mulsd nb314_rinvH2H1(%rsp),%xmm8
    mulsd nb314_rinvMH1(%rsp),%xmm11

    ## move j H1 forces to xmm0-xmm2
    movq nb314_faction(%rbp),%rdi
        movsd 24(%rdi,%rax,8),%xmm0
        movsd 32(%rdi,%rax,8),%xmm1
        movsd 40(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb314_dxH1H1(%rsp),%xmm3
        mulsd nb314_dyH1H1(%rsp),%xmm4
        mulsd nb314_dzH1H1(%rsp),%xmm5
        mulsd nb314_dxH2H1(%rsp),%xmm7
        mulsd nb314_dyH2H1(%rsp),%xmm8
        mulsd nb314_dzH2H1(%rsp),%xmm9
        mulsd nb314_dxMH1(%rsp),%xmm10
        mulsd nb314_dyMH1(%rsp),%xmm11
        mulsd nb314_dzMH1(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb314_fixH1(%rsp),%xmm3
    addsd nb314_fiyH1(%rsp),%xmm4
    addsd nb314_fizH1(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb314_fixH2(%rsp),%xmm7
    addsd nb314_fiyH2(%rsp),%xmm8
    addsd nb314_fizH2(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb314_fixM(%rsp),%xmm10
    addsd nb314_fiyM(%rsp),%xmm11
    addsd nb314_fizM(%rsp),%xmm12

    movsd %xmm3,nb314_fixH1(%rsp)
    movsd %xmm4,nb314_fiyH1(%rsp)
    movsd %xmm5,nb314_fizH1(%rsp)
    movsd %xmm7,nb314_fixH2(%rsp)
    movsd %xmm8,nb314_fiyH2(%rsp)
    movsd %xmm9,nb314_fizH2(%rsp)
    movsd %xmm10,nb314_fixM(%rsp)
    movsd %xmm11,nb314_fiyM(%rsp)
    movsd %xmm12,nb314_fizM(%rsp)

    ## store back j H1 forces from xmm0-xmm2
        movsd %xmm0,24(%rdi,%rax,8)
        movsd %xmm1,32(%rdi,%rax,8)
        movsd %xmm2,40(%rdi,%rax,8)

        ## move j H2 coordinates to local temp variables 
    movq nb314_pos(%rbp),%rsi
    movsd 48(%rsi,%rax,8),%xmm0
    movsd 56(%rsi,%rax,8),%xmm1
    movsd 64(%rsi,%rax,8),%xmm2

    ## xmm0 = H2x
    ## xmm1 = H2y
    ## xmm2 = H2z

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb314_ixH1(%rsp),%xmm0
    subsd nb314_iyH1(%rsp),%xmm1
    subsd nb314_izH1(%rsp),%xmm2
    subsd nb314_ixH2(%rsp),%xmm3
    subsd nb314_iyH2(%rsp),%xmm4
    subsd nb314_izH2(%rsp),%xmm5
    subsd nb314_ixM(%rsp),%xmm6
    subsd nb314_iyM(%rsp),%xmm7
    subsd nb314_izM(%rsp),%xmm8

        movsd %xmm0,nb314_dxH1H2(%rsp)
        movsd %xmm1,nb314_dyH1H2(%rsp)
        movsd %xmm2,nb314_dzH1H2(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb314_dxH2H2(%rsp)
        movsd %xmm4,nb314_dyH2H2(%rsp)
        movsd %xmm5,nb314_dzH2H2(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb314_dxMH2(%rsp)
        movsd %xmm7,nb314_dyMH2(%rsp)
        movsd %xmm8,nb314_dzMH2(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jH2 atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb314_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb314_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1H2 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2H2
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMH2

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb314_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb314_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1H2
        mulsd   %xmm15,%xmm10 ##   rinvH2H2
    mulsd   %xmm15,%xmm11 ##   rinvMH2

        movsd  %xmm9,nb314_rinvH1H2(%rsp)
        movsd  %xmm10,nb314_rinvH2H2(%rsp)
        movsd  %xmm11,nb314_rinvMH2(%rsp)

        ## H2 interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb314_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shlq  $2,%r8
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb314_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb314_epsH1(%rsp)
    movsd    %xmm3,nb314_epsH2(%rsp)
    movsd    %xmm6,nb314_epsM(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb314_epsH1(%rsp),%xmm12
    movsd nb314_epsH2(%rsp),%xmm13
    movsd nb314_epsM(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb314_qqHH(%rsp),%xmm12
    movsd nb314_qqMH(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm12,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm12,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb314_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb314_vctot(%rsp)

    movsd nb314_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb314_rinvH1H2(%rsp),%xmm4
    mulsd nb314_rinvH2H2(%rsp),%xmm8
    mulsd nb314_rinvMH2(%rsp),%xmm11

    ## move j H2 forces to xmm0-xmm2
    movq nb314_faction(%rbp),%rdi
        movsd 48(%rdi,%rax,8),%xmm0
        movsd 56(%rdi,%rax,8),%xmm1
        movsd 64(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb314_dxH1H2(%rsp),%xmm3
        mulsd nb314_dyH1H2(%rsp),%xmm4
        mulsd nb314_dzH1H2(%rsp),%xmm5
        mulsd nb314_dxH2H2(%rsp),%xmm7
        mulsd nb314_dyH2H2(%rsp),%xmm8
        mulsd nb314_dzH2H2(%rsp),%xmm9
        mulsd nb314_dxMH2(%rsp),%xmm10
        mulsd nb314_dyMH2(%rsp),%xmm11
        mulsd nb314_dzMH2(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb314_fixH1(%rsp),%xmm3
    addsd nb314_fiyH1(%rsp),%xmm4
    addsd nb314_fizH1(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb314_fixH2(%rsp),%xmm7
    addsd nb314_fiyH2(%rsp),%xmm8
    addsd nb314_fizH2(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb314_fixM(%rsp),%xmm10
    addsd nb314_fiyM(%rsp),%xmm11
    addsd nb314_fizM(%rsp),%xmm12

    movsd %xmm3,nb314_fixH1(%rsp)
    movsd %xmm4,nb314_fiyH1(%rsp)
    movsd %xmm5,nb314_fizH1(%rsp)
    movsd %xmm7,nb314_fixH2(%rsp)
    movsd %xmm8,nb314_fiyH2(%rsp)
    movsd %xmm9,nb314_fizH2(%rsp)
    movsd %xmm10,nb314_fixM(%rsp)
    movsd %xmm11,nb314_fiyM(%rsp)
    movsd %xmm12,nb314_fizM(%rsp)

    ## store back j H2 forces from xmm0-xmm2
        movsd %xmm0,48(%rdi,%rax,8)
        movsd %xmm1,56(%rdi,%rax,8)
        movsd %xmm2,64(%rdi,%rax,8)

        ## move j M coordinates to local temp variables 
    movq nb314_pos(%rbp),%rsi
    movsd 72(%rsi,%rax,8),%xmm0
    movsd 80(%rsi,%rax,8),%xmm1
    movsd 88(%rsi,%rax,8),%xmm2

    ## xmm0 = Mx
    ## xmm1 = My
    ## xmm2 = Mz

    movsd %xmm0,%xmm3
    movsd %xmm1,%xmm4
    movsd %xmm2,%xmm5
    movsd %xmm0,%xmm6
    movsd %xmm1,%xmm7
    movsd %xmm2,%xmm8

    subsd nb314_ixH1(%rsp),%xmm0
    subsd nb314_iyH1(%rsp),%xmm1
    subsd nb314_izH1(%rsp),%xmm2
    subsd nb314_ixH2(%rsp),%xmm3
    subsd nb314_iyH2(%rsp),%xmm4
    subsd nb314_izH2(%rsp),%xmm5
    subsd nb314_ixM(%rsp),%xmm6
    subsd nb314_iyM(%rsp),%xmm7
    subsd nb314_izM(%rsp),%xmm8

        movsd %xmm0,nb314_dxH1M(%rsp)
        movsd %xmm1,nb314_dyH1M(%rsp)
        movsd %xmm2,nb314_dzH1M(%rsp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movsd %xmm3,nb314_dxH2M(%rsp)
        movsd %xmm4,nb314_dyH2M(%rsp)
        movsd %xmm5,nb314_dzH2M(%rsp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        movsd %xmm6,nb314_dxMM(%rsp)
        movsd %xmm7,nb314_dyMM(%rsp)
        movsd %xmm8,nb314_dzMM(%rsp)
        mulsd  %xmm6,%xmm6
        mulsd  %xmm7,%xmm7
        mulsd  %xmm8,%xmm8
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
    addsd  %xmm7,%xmm6
    addsd  %xmm8,%xmm6

        ## start doing invsqrt for jM atoms
    cvtsd2ss %xmm0,%xmm1
    cvtsd2ss %xmm3,%xmm4
    cvtsd2ss %xmm6,%xmm7
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm4,%xmm4
    rsqrtss %xmm7,%xmm7
    cvtss2sd %xmm1,%xmm1
    cvtss2sd %xmm4,%xmm4
    cvtss2sd %xmm7,%xmm7

        movsd  %xmm1,%xmm2
        movsd  %xmm4,%xmm5
    movsd  %xmm7,%xmm8

        mulsd   %xmm1,%xmm1 ## lu*lu
        mulsd   %xmm4,%xmm4 ## lu*lu
    mulsd   %xmm7,%xmm7 ## lu*lu

        movsd  nb314_three(%rsp),%xmm9
        movsd  %xmm9,%xmm10
    movsd  %xmm9,%xmm11

        mulsd   %xmm0,%xmm1 ## rsq*lu*lu
        mulsd   %xmm3,%xmm4 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm7 ## rsq*lu*lu

        subsd   %xmm1,%xmm9
        subsd   %xmm4,%xmm10
    subsd   %xmm7,%xmm11 ## 3-rsq*lu*lu

        mulsd   %xmm2,%xmm9
        mulsd   %xmm5,%xmm10
    mulsd   %xmm8,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb314_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ## first iteration for rinvH1M 
        mulsd   %xmm15,%xmm10 ## first iteration for rinvH2M
    mulsd   %xmm15,%xmm11 ## first iteration for rinvMM

    ## second iteration step    
        movsd  %xmm9,%xmm2
        movsd  %xmm10,%xmm5
    movsd  %xmm11,%xmm8

        mulsd   %xmm2,%xmm2 ## lu*lu
        mulsd   %xmm5,%xmm5 ## lu*lu
    mulsd   %xmm8,%xmm8 ## lu*lu

        movsd  nb314_three(%rsp),%xmm1
        movsd  %xmm1,%xmm4
    movsd  %xmm1,%xmm7

        mulsd   %xmm0,%xmm2 ## rsq*lu*lu
        mulsd   %xmm3,%xmm5 ## rsq*lu*lu 
    mulsd   %xmm6,%xmm8 ## rsq*lu*lu

        subsd   %xmm2,%xmm1
        subsd   %xmm5,%xmm4
    subsd   %xmm8,%xmm7 ## 3-rsq*lu*lu

        mulsd   %xmm1,%xmm9
        mulsd   %xmm4,%xmm10
    mulsd   %xmm7,%xmm11 ## lu*(3-rsq*lu*lu)

        movsd  nb314_half(%rsp),%xmm15
        mulsd   %xmm15,%xmm9 ##  rinvH1M
        mulsd   %xmm15,%xmm10 ##   rinvH2M
    mulsd   %xmm15,%xmm11 ##   rinvMM

        movsd  %xmm9,nb314_rinvH1M(%rsp)
        movsd  %xmm10,nb314_rinvH2M(%rsp)
        movsd  %xmm11,nb314_rinvMM(%rsp)

        ## M interactions 
    ## rsq in xmm0,xmm3,xmm6  
    ## rinv in xmm9, xmm10, xmm11

    movsd  nb314_tsc(%rsp),%xmm1
    mulsd  %xmm9,%xmm0 ## r
    mulsd  %xmm10,%xmm3
    mulsd  %xmm11,%xmm6
    mulsd  %xmm1,%xmm0 ## rtab
    mulsd  %xmm1,%xmm3
    mulsd  %xmm1,%xmm6

    ## truncate and convert to integers
    cvttsd2si %xmm0,%r8d
    cvttsd2si %xmm3,%r10d
    cvttsd2si %xmm6,%r12d

    ## convert back to float
    cvtsi2sd  %r8d,%xmm2
    cvtsi2sd  %r10d,%xmm5
    cvtsi2sd  %r12d,%xmm8

    ## multiply by 4
    shlq  $2,%r8
    shll  $2,%r10d
    shll  $2,%r12d

    movq nb314_VFtab(%rbp),%rsi

    ## calculate eps
    subsd     %xmm2,%xmm0
    subsd     %xmm5,%xmm3
    subsd     %xmm8,%xmm6

    movsd    %xmm0,nb314_epsH1(%rsp)
    movsd    %xmm3,nb314_epsH2(%rsp)
    movsd    %xmm6,nb314_epsM(%rsp)

    ## Load LOTS of table data
    movsd (%rsi,%r8,8),%xmm0
    movsd 8(%rsi,%r8,8),%xmm1
    movsd 16(%rsi,%r8,8),%xmm2
    movsd 24(%rsi,%r8,8),%xmm3
    movsd (%rsi,%r10,8),%xmm4
    movsd 8(%rsi,%r10,8),%xmm5
    movsd 16(%rsi,%r10,8),%xmm6
    movsd 24(%rsi,%r10,8),%xmm7
    movsd (%rsi,%r12,8),%xmm8
    movsd 8(%rsi,%r12,8),%xmm9
    movsd 16(%rsi,%r12,8),%xmm10
    movsd 24(%rsi,%r12,8),%xmm11
    ## table data ready in xmm0-xmm3 , xmm4-xmm7 , and xmm8-xmm11

    movsd nb314_epsH1(%rsp),%xmm12
    movsd nb314_epsH2(%rsp),%xmm13
    movsd nb314_epsM(%rsp),%xmm14

    mulsd  %xmm12,%xmm3  ## Heps
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11
    mulsd  %xmm12,%xmm2  ## Geps
    mulsd  %xmm13,%xmm6
    mulsd  %xmm14,%xmm10
    mulsd  %xmm12,%xmm3  ## Heps2
    mulsd  %xmm13,%xmm7
    mulsd  %xmm14,%xmm11

    addsd  %xmm2,%xmm1  ## F+Geps
    addsd  %xmm6,%xmm5
    addsd  %xmm10,%xmm9
    addsd  %xmm3,%xmm1  ## F+Geps+Heps2 = Fp
    addsd  %xmm7,%xmm5
    addsd  %xmm11,%xmm9
    addsd  %xmm3,%xmm3   ## 2*Heps2
    addsd  %xmm7,%xmm7
    addsd  %xmm11,%xmm11
    addsd  %xmm2,%xmm3   ## 2*Heps2+Geps
    addsd  %xmm6,%xmm7
    addsd  %xmm10,%xmm11
    addsd  %xmm1,%xmm3  ## FF = Fp + 2*Heps2 + Geps
    addsd  %xmm5,%xmm7
    addsd  %xmm9,%xmm11
    mulsd  %xmm12,%xmm1  ## eps*Fp
    mulsd  %xmm13,%xmm5
    mulsd  %xmm14,%xmm9
    movsd nb314_qqMH(%rsp),%xmm12
    movsd nb314_qqMM(%rsp),%xmm13
    addsd  %xmm0,%xmm1    ## VV
    addsd  %xmm4,%xmm5
    addsd  %xmm8,%xmm9
    mulsd  %xmm12,%xmm1  ## VV*qq = vcoul
    mulsd  %xmm12,%xmm5
    mulsd  %xmm13,%xmm9
    mulsd  %xmm12,%xmm3   ## FF*qq = fij
    mulsd  %xmm12,%xmm7
    mulsd  %xmm13,%xmm11

    ## accumulate vctot
    addsd  nb314_vctot(%rsp),%xmm1
    addsd  %xmm9,%xmm5
    addsd  %xmm5,%xmm1
    movsd %xmm1,nb314_vctot(%rsp)

    movsd nb314_tsc(%rsp),%xmm10
    mulsd  %xmm10,%xmm3 ## fscal
    mulsd  %xmm10,%xmm7
    mulsd  %xmm11,%xmm10

    xorpd %xmm4,%xmm4
    xorpd %xmm8,%xmm8
    xorpd %xmm11,%xmm11

    subsd %xmm3,%xmm4
    subsd %xmm7,%xmm8
    subsd %xmm10,%xmm11

    mulsd nb314_rinvH1M(%rsp),%xmm4
    mulsd nb314_rinvH2M(%rsp),%xmm8
    mulsd nb314_rinvMM(%rsp),%xmm11

    ## move j M forces to xmm0-xmm2
    movq nb314_faction(%rbp),%rdi
        movsd 72(%rdi,%rax,8),%xmm0
        movsd 80(%rdi,%rax,8),%xmm1
        movsd 88(%rdi,%rax,8),%xmm2

    movapd %xmm4,%xmm3
    movapd %xmm4,%xmm5
    movapd %xmm8,%xmm7
    movapd %xmm8,%xmm9
    movapd %xmm11,%xmm10
    movapd %xmm11,%xmm12

        mulsd nb314_dxH1M(%rsp),%xmm3
        mulsd nb314_dyH1M(%rsp),%xmm4
        mulsd nb314_dzH1M(%rsp),%xmm5
        mulsd nb314_dxH2M(%rsp),%xmm7
        mulsd nb314_dyH2M(%rsp),%xmm8
        mulsd nb314_dzH2M(%rsp),%xmm9
        mulsd nb314_dxMM(%rsp),%xmm10
        mulsd nb314_dyMM(%rsp),%xmm11
        mulsd nb314_dzMM(%rsp),%xmm12

    addsd %xmm3,%xmm0
    addsd %xmm4,%xmm1
    addsd %xmm5,%xmm2
    addsd nb314_fixH1(%rsp),%xmm3
    addsd nb314_fiyH1(%rsp),%xmm4
    addsd nb314_fizH1(%rsp),%xmm5

    addsd %xmm7,%xmm0
    addsd %xmm8,%xmm1
    addsd %xmm9,%xmm2
    addsd nb314_fixH2(%rsp),%xmm7
    addsd nb314_fiyH2(%rsp),%xmm8
    addsd nb314_fizH2(%rsp),%xmm9

    addsd %xmm10,%xmm0
    addsd %xmm11,%xmm1
    addsd %xmm12,%xmm2
    addsd nb314_fixM(%rsp),%xmm10
    addsd nb314_fiyM(%rsp),%xmm11
    addsd nb314_fizM(%rsp),%xmm12

    movsd %xmm3,nb314_fixH1(%rsp)
    movsd %xmm4,nb314_fiyH1(%rsp)
    movsd %xmm5,nb314_fizH1(%rsp)
    movsd %xmm7,nb314_fixH2(%rsp)
    movsd %xmm8,nb314_fiyH2(%rsp)
    movsd %xmm9,nb314_fizH2(%rsp)
    movsd %xmm10,nb314_fixM(%rsp)
    movsd %xmm11,nb314_fiyM(%rsp)
    movsd %xmm12,nb314_fizM(%rsp)

    ## store back j M forces from xmm0-xmm2
        movsd %xmm0,72(%rdi,%rax,8)
        movsd %xmm1,80(%rdi,%rax,8)
        movsd %xmm2,88(%rdi,%rax,8)

_nb_kernel314_x86_64_sse2.nb314_updateouterdata: 
        movl  nb314_ii3(%rsp),%ecx
        movq  nb314_faction(%rbp),%rdi
        movq  nb314_fshift(%rbp),%rsi
        movl  nb314_is3(%rsp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb314_fixO(%rsp),%xmm0
        movapd nb314_fiyO(%rsp),%xmm1
        movapd nb314_fizO(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

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
        movapd nb314_fixH1(%rsp),%xmm0
        movapd nb314_fiyH1(%rsp),%xmm1
        movapd nb314_fizH1(%rsp),%xmm2

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
        movapd nb314_fixH2(%rsp),%xmm0
        movapd nb314_fiyH2(%rsp),%xmm1
        movapd nb314_fizH2(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

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
        movapd nb314_fixM(%rsp),%xmm0
        movapd nb314_fiyM(%rsp),%xmm1
        movapd nb314_fizM(%rsp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5

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
        movl nb314_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb314_gid(%rbp),%rdx              ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb314_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb314_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb314_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb314_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb314_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel314_x86_64_sse2.nb314_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb314_n(%rsp)
        jmp _nb_kernel314_x86_64_sse2.nb314_outer
_nb_kernel314_x86_64_sse2.nb314_outerend: 
        ## check if more outer neighborlists remain
        movl  nb314_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel314_x86_64_sse2.nb314_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel314_x86_64_sse2.nb314_threadloop
_nb_kernel314_x86_64_sse2.nb314_end: 
        movl nb314_nouter(%rsp),%eax
        movl nb314_ninner(%rsp),%ebx
        movq nb314_outeriter(%rbp),%rcx
        movq nb314_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $1880,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret






.globl nb_kernel314nf_x86_64_sse2
.globl _nb_kernel314nf_x86_64_sse2
nb_kernel314nf_x86_64_sse2:     
_nb_kernel314nf_x86_64_sse2:    
##      Room for return address and rbp (16 bytes)
.set nb314nf_fshift, 16
.set nb314nf_gid, 24
.set nb314nf_pos, 32
.set nb314nf_faction, 40
.set nb314nf_charge, 48
.set nb314nf_p_facel, 56
.set nb314nf_argkrf, 64
.set nb314nf_argcrf, 72
.set nb314nf_Vc, 80
.set nb314nf_type, 88
.set nb314nf_p_ntype, 96
.set nb314nf_vdwparam, 104
.set nb314nf_Vvdw, 112
.set nb314nf_p_tabscale, 120
.set nb314nf_VFtab, 128
.set nb314nf_invsqrta, 136
.set nb314nf_dvda, 144
.set nb314nf_p_gbtabscale, 152
.set nb314nf_GBtab, 160
.set nb314nf_p_nthreads, 168
.set nb314nf_count, 176
.set nb314nf_mtx, 184
.set nb314nf_outeriter, 192
.set nb314nf_inneriter, 200
.set nb314nf_work, 208
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb314nf_ixO, 0
.set nb314nf_iyO, 16
.set nb314nf_izO, 32
.set nb314nf_ixH1, 48
.set nb314nf_iyH1, 64
.set nb314nf_izH1, 80
.set nb314nf_ixH2, 96
.set nb314nf_iyH2, 112
.set nb314nf_izH2, 128
.set nb314nf_ixM, 144
.set nb314nf_iyM, 160
.set nb314nf_izM, 176
.set nb314nf_jxO, 192
.set nb314nf_jyO, 208
.set nb314nf_jzO, 224
.set nb314nf_jxH1, 240
.set nb314nf_jyH1, 256
.set nb314nf_jzH1, 272
.set nb314nf_jxH2, 288
.set nb314nf_jyH2, 304
.set nb314nf_jzH2, 320
.set nb314nf_jxM, 336
.set nb314nf_jyM, 352
.set nb314nf_jzM, 368
.set nb314nf_qqMM, 384
.set nb314nf_qqMH, 400
.set nb314nf_qqHH, 416
.set nb314nf_tsc, 432
.set nb314nf_c6, 448
.set nb314nf_c12, 464
.set nb314nf_vctot, 480
.set nb314nf_Vvdwtot, 496
.set nb314nf_half, 512
.set nb314nf_three, 528
.set nb314nf_rsqOO, 544
.set nb314nf_rsqH1H1, 560
.set nb314nf_rsqH1H2, 576
.set nb314nf_rsqH1M, 592
.set nb314nf_rsqH2H1, 608
.set nb314nf_rsqH2H2, 624
.set nb314nf_rsqH2M, 640
.set nb314nf_rsqMH1, 656
.set nb314nf_rsqMH2, 672
.set nb314nf_rsqMM, 688
.set nb314nf_rinvsqOO, 704
.set nb314nf_rinvH1H1, 720
.set nb314nf_rinvH1H2, 736
.set nb314nf_rinvH1M, 752
.set nb314nf_rinvH2H1, 768
.set nb314nf_rinvH2H2, 784
.set nb314nf_rinvH2M, 800
.set nb314nf_rinvMH1, 816
.set nb314nf_rinvMH2, 832
.set nb314nf_rinvMM, 848
.set nb314nf_is3, 864
.set nb314nf_ii3, 868
.set nb314nf_nri, 872
.set nb314nf_iinr, 880
.set nb314nf_jindex, 888
.set nb314nf_jjnr, 896
.set nb314nf_shift, 904
.set nb314nf_shiftvec, 912
.set nb314nf_facel, 920
.set nb314nf_innerjjnr, 928
.set nb314nf_innerk, 936
.set nb314nf_n, 940
.set nb314nf_nn1, 944
.set nb314nf_nouter, 948
.set nb314nf_ninner, 952
        push %rbp
        movq %rsp,%rbp
        push %rbx
        emms

        push %r12
        push %r13
        push %r14
        push %r15

        subq $968,%rsp          ## local variable stack space (n*16+8)

        ## zero 32-bit iteration counters
        movl $0,%eax
        movl %eax,nb314nf_nouter(%rsp)
        movl %eax,nb314nf_ninner(%rsp)

        movl (%rdi),%edi
        movl %edi,nb314nf_nri(%rsp)
        movq %rsi,nb314nf_iinr(%rsp)
        movq %rdx,nb314nf_jindex(%rsp)
        movq %rcx,nb314nf_jjnr(%rsp)
        movq %r8,nb314nf_shift(%rsp)
        movq %r9,nb314nf_shiftvec(%rsp)
        movq nb314nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd %xmm0,nb314nf_facel(%rsp)

        movq nb314nf_p_tabscale(%rbp),%rax
        movsd (%rax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb314nf_tsc(%rsp)

        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double half IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb314nf_half(%rsp)
        movl %ebx,nb314nf_half+4(%rsp)
        movsd nb314nf_half(%rsp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## one
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## two
        addpd  %xmm2,%xmm3      ## three
        movapd %xmm1,nb314nf_half(%rsp)
        movapd %xmm3,nb314nf_three(%rsp)

        ## assume we have at least one i particle - start directly 
        movq  nb314nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx),%ebx           ## ebx =ii 

        movq  nb314nf_charge(%rbp),%rdx
        movsd 24(%rdx,%rbx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%rdx,%rbx,8),%xmm5
        movq nb314nf_p_facel(%rbp),%rsi
        movsd (%rsi),%xmm0
        movsd nb314nf_facel(%rsp),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb314nf_qqMM(%rsp)
        movapd %xmm4,nb314nf_qqMH(%rsp)
        movapd %xmm5,nb314nf_qqHH(%rsp)

        xorpd %xmm0,%xmm0
        movq  nb314nf_type(%rbp),%rdx
        movl  (%rdx,%rbx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movq nb314nf_p_ntype(%rbp),%rdi
        imull (%rdi),%ecx     ## rcx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movq  nb314nf_vdwparam(%rbp),%rax
        movlpd (%rax,%rdx,8),%xmm0
        movlpd 8(%rax,%rdx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb314nf_c6(%rsp)
        movapd %xmm1,nb314nf_c12(%rsp)

_nb_kernel314nf_x86_64_sse2.nb314nf_threadloop: 
        movq  nb314nf_count(%rbp),%rsi            ## pointer to sync counter
        movl  (%rsi),%eax
_nb_kernel314nf_x86_64_sse2.nb314nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%rsi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel314nf_x86_64_sse2.nb314nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb314nf_nri(%rsp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb314nf_n(%rsp)
        movl %ebx,nb314nf_nn1(%rsp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel314nf_x86_64_sse2.nb314nf_outerstart
        jmp _nb_kernel314nf_x86_64_sse2.nb314nf_end

_nb_kernel314nf_x86_64_sse2.nb314nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb314nf_nouter(%rsp),%ebx
        movl %ebx,nb314nf_nouter(%rsp)

_nb_kernel314nf_x86_64_sse2.nb314nf_outer: 
        movq  nb314nf_shift(%rsp),%rax        ## rax = pointer into shift[] 
        movl  (%rax,%rsi,4),%ebx        ## rbx=shift[n] 

        lea  (%rbx,%rbx,2),%rbx    ## rbx=3*is 
        movl  %ebx,nb314nf_is3(%rsp)            ## store is3 

        movq  nb314nf_shiftvec(%rsp),%rax     ## rax = base of shiftvec[] 

        movsd (%rax,%rbx,8),%xmm0
        movsd 8(%rax,%rbx,8),%xmm1
        movsd 16(%rax,%rbx,8),%xmm2

        movq  nb314nf_iinr(%rsp),%rcx         ## rcx = pointer into iinr[]      
        movl  (%rcx,%rsi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        lea  (%rbx,%rbx,2),%rbx        ## rbx = 3*ii=ii3 
        movq  nb314nf_pos(%rbp),%rax      ## rax = base of pos[]  
        movl  %ebx,nb314nf_ii3(%rsp)

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
        movapd %xmm3,nb314nf_ixO(%rsp)
        movapd %xmm4,nb314nf_iyO(%rsp)
        movapd %xmm5,nb314nf_izO(%rsp)
        movapd %xmm6,nb314nf_ixH1(%rsp)
        movapd %xmm7,nb314nf_iyH1(%rsp)

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
        movapd %xmm6,nb314nf_izH1(%rsp)
        movapd %xmm0,nb314nf_ixH2(%rsp)
        movapd %xmm1,nb314nf_iyH2(%rsp)
        movapd %xmm2,nb314nf_izH2(%rsp)
        movapd %xmm3,nb314nf_ixM(%rsp)
        movapd %xmm4,nb314nf_iyM(%rsp)
        movapd %xmm5,nb314nf_izM(%rsp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb314nf_vctot(%rsp)
        movapd %xmm4,nb314nf_Vvdwtot(%rsp)

        movq  nb314nf_jindex(%rsp),%rax
        movl  (%rax,%rsi,4),%ecx             ## jindex[n] 
        movl  4(%rax,%rsi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movq  nb314nf_pos(%rbp),%rsi
        movq  nb314nf_faction(%rbp),%rdi
        movq  nb314nf_jjnr(%rsp),%rax
        shll  $2,%ecx
        addq  %rcx,%rax
        movq  %rax,nb314nf_innerjjnr(%rsp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb314nf_ninner(%rsp),%ecx
        movl  %ecx,nb314nf_ninner(%rsp)
        addl  $0,%edx
        movl  %edx,nb314nf_innerk(%rsp)      ## number of innerloop atoms 
        jge   _nb_kernel314nf_x86_64_sse2.nb314nf_unroll_loop
        jmp   _nb_kernel314nf_x86_64_sse2.nb314nf_checksingle
_nb_kernel314nf_x86_64_sse2.nb314nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movq  nb314nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax
        movl  4(%rdx),%ebx

        addq $8,nb314nf_innerjjnr(%rsp)            ## advance pointer (unrolled 2) 

        movq nb314nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 
        lea  (%rbx,%rbx,2),%rbx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%rsi,%rax,8),%xmm0
        movlpd (%rsi,%rbx,8),%xmm2
        movhpd 8(%rsi,%rax,8),%xmm0
        movhpd 8(%rsi,%rbx,8),%xmm2
        movlpd 16(%rsi,%rax,8),%xmm3
        movlpd 16(%rsi,%rbx,8),%xmm5
        movhpd 24(%rsi,%rax,8),%xmm3
        movhpd 24(%rsi,%rbx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb314nf_jxO(%rsp)
        movapd  %xmm1,nb314nf_jyO(%rsp)
        movapd  %xmm3,nb314nf_jzO(%rsp)
        movapd  %xmm4,nb314nf_jxH1(%rsp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%rsi,%rax,8),%xmm0
        movlpd 32(%rsi,%rbx,8),%xmm2
        movhpd 40(%rsi,%rax,8),%xmm0
        movhpd 40(%rsi,%rbx,8),%xmm2
        movlpd 48(%rsi,%rax,8),%xmm3
        movlpd 48(%rsi,%rbx,8),%xmm5
        movhpd 56(%rsi,%rax,8),%xmm3
        movhpd 56(%rsi,%rbx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb314nf_jyH1(%rsp)
        movapd  %xmm1,nb314nf_jzH1(%rsp)
        movapd  %xmm3,nb314nf_jxH2(%rsp)
        movapd  %xmm4,nb314nf_jyH2(%rsp)

        ## load h2z, mx, my, mz
        movlpd 64(%rsi,%rax,8),%xmm0
        movlpd 64(%rsi,%rbx,8),%xmm2
        movhpd 72(%rsi,%rax,8),%xmm0
        movhpd 72(%rsi,%rbx,8),%xmm2
        movlpd 80(%rsi,%rax,8),%xmm3
        movlpd 80(%rsi,%rbx,8),%xmm5
        movhpd 88(%rsi,%rax,8),%xmm3
        movhpd 88(%rsi,%rbx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb314nf_jzH2(%rsp)
        movapd  %xmm1,nb314nf_jxM(%rsp)
        movapd  %xmm3,nb314nf_jyM(%rsp)
        movapd  %xmm4,nb314nf_jzM(%rsp)

        ## start calculating pairwise distances
        movapd nb314nf_ixO(%rsp),%xmm0
        movapd nb314nf_iyO(%rsp),%xmm1
        movapd nb314nf_izO(%rsp),%xmm2
        movapd nb314nf_ixH1(%rsp),%xmm3
        movapd nb314nf_iyH1(%rsp),%xmm4
        movapd nb314nf_izH1(%rsp),%xmm5
        subpd  nb314nf_jxO(%rsp),%xmm0
        subpd  nb314nf_jyO(%rsp),%xmm1
        subpd  nb314nf_jzO(%rsp),%xmm2
        subpd  nb314nf_jxH1(%rsp),%xmm3
        subpd  nb314nf_jyH1(%rsp),%xmm4
        subpd  nb314nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb314nf_rsqOO(%rsp)
        movapd %xmm3,nb314nf_rsqH1H1(%rsp)

        movapd nb314nf_ixH1(%rsp),%xmm0
        movapd nb314nf_iyH1(%rsp),%xmm1
        movapd nb314nf_izH1(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314nf_jxH2(%rsp),%xmm0
        subpd  nb314nf_jyH2(%rsp),%xmm1
        subpd  nb314nf_jzH2(%rsp),%xmm2
        subpd  nb314nf_jxM(%rsp),%xmm3
        subpd  nb314nf_jyM(%rsp),%xmm4
        subpd  nb314nf_jzM(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb314nf_rsqH1H2(%rsp)
        movapd %xmm3,nb314nf_rsqH1M(%rsp)

        movapd nb314nf_ixH2(%rsp),%xmm0
        movapd nb314nf_iyH2(%rsp),%xmm1
        movapd nb314nf_izH2(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314nf_jxH1(%rsp),%xmm0
        subpd  nb314nf_jyH1(%rsp),%xmm1
        subpd  nb314nf_jzH1(%rsp),%xmm2
        subpd  nb314nf_jxH2(%rsp),%xmm3
        subpd  nb314nf_jyH2(%rsp),%xmm4
        subpd  nb314nf_jzH2(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb314nf_rsqH2H1(%rsp)
        movapd %xmm3,nb314nf_rsqH2H2(%rsp)

        movapd nb314nf_ixH2(%rsp),%xmm0
        movapd nb314nf_iyH2(%rsp),%xmm1
        movapd nb314nf_izH2(%rsp),%xmm2
        movapd nb314nf_ixM(%rsp),%xmm3
        movapd nb314nf_iyM(%rsp),%xmm4
        movapd nb314nf_izM(%rsp),%xmm5
        subpd  nb314nf_jxM(%rsp),%xmm0
        subpd  nb314nf_jyM(%rsp),%xmm1
        subpd  nb314nf_jzM(%rsp),%xmm2
        subpd  nb314nf_jxH1(%rsp),%xmm3
        subpd  nb314nf_jyH1(%rsp),%xmm4
        subpd  nb314nf_jzH1(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb314nf_rsqH2M(%rsp)
        movapd %xmm4,nb314nf_rsqMH1(%rsp)

        movapd nb314nf_ixM(%rsp),%xmm0
        movapd nb314nf_iyM(%rsp),%xmm1
        movapd nb314nf_izM(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314nf_jxH2(%rsp),%xmm0
        subpd  nb314nf_jyH2(%rsp),%xmm1
        subpd  nb314nf_jzH2(%rsp),%xmm2
        subpd  nb314nf_jxM(%rsp),%xmm3
        subpd  nb314nf_jyM(%rsp),%xmm4
        subpd  nb314nf_jzM(%rsp),%xmm5
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb314nf_rsqMH2(%rsp)
        movapd %xmm4,nb314nf_rsqMM(%rsp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1  ## luA
        cvtps2pd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb314nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvMH2(%rsp)
        movapd %xmm5,nb314nf_rinvMM(%rsp)

        movapd nb314nf_rsqOO(%rsp),%xmm0
        movapd nb314nf_rsqH1H1(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm3   ## iter1 of  
        mulpd   nb314nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%rsp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb314nf_rinvsqOO(%rsp)
        movapd %xmm5,nb314nf_rinvH1H1(%rsp)

        movapd nb314nf_rsqH1H2(%rsp),%xmm0
        movapd nb314nf_rsqH1M(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm3   ## iter1 
        mulpd   nb314nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvH1H2(%rsp)
        movapd %xmm5,nb314nf_rinvH1M(%rsp)

        movapd nb314nf_rsqH2H1(%rsp),%xmm0
        movapd nb314nf_rsqH2H2(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb314nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvH2H1(%rsp)
        movapd %xmm5,nb314nf_rinvH2H2(%rsp)

        movapd nb314nf_rsqMH1(%rsp),%xmm0
        movapd nb314nf_rsqH2M(%rsp),%xmm4
        cvtpd2ps %xmm0,%xmm1
        cvtpd2ps %xmm4,%xmm5
        rsqrtps %xmm1,%xmm1
        rsqrtps %xmm5,%xmm5
        cvtps2pd %xmm1,%xmm1
        cvtps2pd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        mulpd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm3   ## iter1a 
        mulpd   nb314nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvMH1(%rsp)
        movapd %xmm5,nb314nf_rinvH2M(%rsp)

        ## start with OO interaction 
        movapd nb314nf_rinvsqOO(%rsp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb314nf_c6(%rsp),%xmm1
        mulpd  nb314nf_c12(%rsp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb314nf_Vvdwtot(%rsp),%xmm3
        movapd %xmm3,nb314nf_Vvdwtot(%rsp)

        ## H1-H1 interaction 
        movapd nb314nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1

        movd %eax,%mm0
        movd %ebx,%mm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## H1-H2 interaction  
        movapd nb314nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## H1-M interaction 
        movapd nb314nf_rinvH1M(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH1M(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## H2-H1 interaction 
        movapd nb314nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## H2-H2 interaction 
        movapd nb314nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## H2-M interaction 
        movapd nb314nf_rinvH2M(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH2M(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## M-H1 interaction 
        movapd nb314nf_rinvMH1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqMH1(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## M-H2 interaction 
        movapd nb314nf_rinvMH2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqMH2(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## M-M interaction 
        movapd nb314nf_rinvMM(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqMM(%rsp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%rsp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd (%rsi,%rbx,8),%xmm3      ## Y2 F2 
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd 16(%rsi,%rbx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMM(%rsp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%rsp),%xmm5
        movapd %xmm5,nb314nf_vctot(%rsp)

        ## should we do one more iteration? 
        subl $2,nb314nf_innerk(%rsp)
        jl    _nb_kernel314nf_x86_64_sse2.nb314nf_checksingle
        jmp   _nb_kernel314nf_x86_64_sse2.nb314nf_unroll_loop
_nb_kernel314nf_x86_64_sse2.nb314nf_checksingle: 
        movl  nb314nf_innerk(%rsp),%edx
        andl  $1,%edx
        jnz   _nb_kernel314nf_x86_64_sse2.nb314nf_dosingle
        jmp   _nb_kernel314nf_x86_64_sse2.nb314nf_updateouterdata
_nb_kernel314nf_x86_64_sse2.nb314nf_dosingle: 
        movq  nb314nf_innerjjnr(%rsp),%rdx       ## pointer to jjnr[k] 
        movl  (%rdx),%eax

        movq nb314nf_pos(%rbp),%rsi        ## base of pos[] 

        lea  (%rax,%rax,2),%rax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%rsi,%rax,8),%xmm0
        movhpd 8(%rsi,%rax,8),%xmm0
        movlpd 16(%rsi,%rax,8),%xmm1
        movhpd 24(%rsi,%rax,8),%xmm1
        movlpd 32(%rsi,%rax,8),%xmm2
        movhpd 40(%rsi,%rax,8),%xmm2
        movlpd 48(%rsi,%rax,8),%xmm3
        movhpd 56(%rsi,%rax,8),%xmm3
        movlpd 64(%rsi,%rax,8),%xmm4
        movhpd 72(%rsi,%rax,8),%xmm4
        movlpd 80(%rsi,%rax,8),%xmm5
        movhpd 88(%rsi,%rax,8),%xmm5
        movsd  %xmm0,nb314nf_jxO(%rsp)
        movsd  %xmm1,nb314nf_jzO(%rsp)
        movsd  %xmm2,nb314nf_jyH1(%rsp)
        movsd  %xmm3,nb314nf_jxH2(%rsp)
        movsd  %xmm4,nb314nf_jzH2(%rsp)
        movsd  %xmm5,nb314nf_jyM(%rsp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb314nf_jyO(%rsp)
        movsd  %xmm1,nb314nf_jxH1(%rsp)
        movsd  %xmm2,nb314nf_jzH1(%rsp)
        movsd  %xmm3,nb314nf_jyH2(%rsp)
        movsd  %xmm4,nb314nf_jxM(%rsp)
        movsd  %xmm5,nb314nf_jzM(%rsp)

        ## start calculating pairwise distances
        movapd nb314nf_ixO(%rsp),%xmm0
        movapd nb314nf_iyO(%rsp),%xmm1
        movapd nb314nf_izO(%rsp),%xmm2
        movapd nb314nf_ixH1(%rsp),%xmm3
        movapd nb314nf_iyH1(%rsp),%xmm4
        movapd nb314nf_izH1(%rsp),%xmm5
        subsd  nb314nf_jxO(%rsp),%xmm0
        subsd  nb314nf_jyO(%rsp),%xmm1
        subsd  nb314nf_jzO(%rsp),%xmm2
        subsd  nb314nf_jxH1(%rsp),%xmm3
        subsd  nb314nf_jyH1(%rsp),%xmm4
        subsd  nb314nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb314nf_rsqOO(%rsp)
        movapd %xmm3,nb314nf_rsqH1H1(%rsp)

        movapd nb314nf_ixH1(%rsp),%xmm0
        movapd nb314nf_iyH1(%rsp),%xmm1
        movapd nb314nf_izH1(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314nf_jxH2(%rsp),%xmm0
        subsd  nb314nf_jyH2(%rsp),%xmm1
        subsd  nb314nf_jzH2(%rsp),%xmm2
        subsd  nb314nf_jxM(%rsp),%xmm3
        subsd  nb314nf_jyM(%rsp),%xmm4
        subsd  nb314nf_jzM(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb314nf_rsqH1H2(%rsp)
        movapd %xmm3,nb314nf_rsqH1M(%rsp)

        movapd nb314nf_ixH2(%rsp),%xmm0
        movapd nb314nf_iyH2(%rsp),%xmm1
        movapd nb314nf_izH2(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314nf_jxH1(%rsp),%xmm0
        subsd  nb314nf_jyH1(%rsp),%xmm1
        subsd  nb314nf_jzH1(%rsp),%xmm2
        subsd  nb314nf_jxH2(%rsp),%xmm3
        subsd  nb314nf_jyH2(%rsp),%xmm4
        subsd  nb314nf_jzH2(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb314nf_rsqH2H1(%rsp)
        movapd %xmm3,nb314nf_rsqH2H2(%rsp)

        movapd nb314nf_ixH2(%rsp),%xmm0
        movapd nb314nf_iyH2(%rsp),%xmm1
        movapd nb314nf_izH2(%rsp),%xmm2
        movapd nb314nf_ixM(%rsp),%xmm3
        movapd nb314nf_iyM(%rsp),%xmm4
        movapd nb314nf_izM(%rsp),%xmm5
        subsd  nb314nf_jxM(%rsp),%xmm0
        subsd  nb314nf_jyM(%rsp),%xmm1
        subsd  nb314nf_jzM(%rsp),%xmm2
        subsd  nb314nf_jxH1(%rsp),%xmm3
        subsd  nb314nf_jyH1(%rsp),%xmm4
        subsd  nb314nf_jzH1(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb314nf_rsqH2M(%rsp)
        movapd %xmm4,nb314nf_rsqMH1(%rsp)

        movapd nb314nf_ixM(%rsp),%xmm0
        movapd nb314nf_iyM(%rsp),%xmm1
        movapd nb314nf_izM(%rsp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314nf_jxH2(%rsp),%xmm0
        subsd  nb314nf_jyH2(%rsp),%xmm1
        subsd  nb314nf_jzH2(%rsp),%xmm2
        subsd  nb314nf_jxM(%rsp),%xmm3
        subsd  nb314nf_jyM(%rsp),%xmm4
        subsd  nb314nf_jzM(%rsp),%xmm5
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb314nf_rsqMH2(%rsp)
        movapd %xmm4,nb314nf_rsqMM(%rsp)

        ## Invsqrt form rsq M-H2 (rsq in xmm0) and MM (rsq in xmm4) 
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1  ## luA
        cvtss2sd %xmm5,%xmm5  ## luB

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb314nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvMH2(%rsp)
        movapd %xmm5,nb314nf_rinvMM(%rsp)

        movapd nb314nf_rsqOO(%rsp),%xmm0
        movapd nb314nf_rsqH1H1(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm3   ## iter1 of  
        mulsd   nb314nf_half(%rsp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%rsp),%xmm5   ## rinv
        mulsd   %xmm1,%xmm1
        movapd %xmm1,nb314nf_rinvsqOO(%rsp)
        movapd %xmm5,nb314nf_rinvH1H1(%rsp)

        movapd nb314nf_rsqH1H2(%rsp),%xmm0
        movapd nb314nf_rsqH1M(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm3   ## iter1 
        mulsd   nb314nf_half(%rsp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvH1H2(%rsp)
        movapd %xmm5,nb314nf_rinvH1M(%rsp)

        movapd nb314nf_rsqH2H1(%rsp),%xmm0
        movapd nb314nf_rsqH2H2(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb314nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvH2H1(%rsp)
        movapd %xmm5,nb314nf_rinvH2H2(%rsp)

        movapd nb314nf_rsqMH1(%rsp),%xmm0
        movapd nb314nf_rsqH2M(%rsp),%xmm4
        cvtsd2ss %xmm0,%xmm1
        cvtsd2ss %xmm4,%xmm5
        rsqrtss %xmm1,%xmm1
        rsqrtss %xmm5,%xmm5
        cvtss2sd %xmm1,%xmm1
        cvtss2sd %xmm5,%xmm5

        movapd  %xmm1,%xmm2     ## copy of luA 
        movapd  %xmm5,%xmm6     ## copy of luB 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        mulsd   %xmm5,%xmm5     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm3   ## iter1a 
        mulsd   nb314nf_half(%rsp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%rsp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%rsp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%rsp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvMH1(%rsp)
        movapd %xmm5,nb314nf_rinvH2M(%rsp)


        ## start with OO interaction 
        movsd nb314nf_rinvsqOO(%rsp),%xmm0   ## xmm0=rinvsq
        movsd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movsd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb314nf_c6(%rsp),%xmm1
        mulsd  nb314nf_c12(%rsp),%xmm2
        movsd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb314nf_Vvdwtot(%rsp),%xmm3
        movsd %xmm3,nb314nf_Vvdwtot(%rsp)

        ## H1-H1 interaction 
        movapd nb314nf_rinvH1H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH1H1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1

        movd %eax,%mm0

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %rax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  


        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

        ## H1-H2 interaction  
        movapd nb314nf_rinvH1H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH1H2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

        ## H1-M interaction 
        movapd nb314nf_rinvH1M(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH1M(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

        ## H2-H1 interaction 
        movapd nb314nf_rinvH2H1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH2H1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

        ## H2-H2 interaction 
        movapd nb314nf_rinvH2H2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH2H2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

        ## H2-M interaction 
        movapd nb314nf_rinvH2M(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH2M(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

        ## M-H1 interaction 
        movapd nb314nf_rinvMH1(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqMH1(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

        ## M-H2 interaction 
        movapd nb314nf_rinvMH2(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqMH2(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

        ## M-M interaction 
        movapd nb314nf_rinvMM(%rsp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqMM(%rsp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%rsp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movq nb314nf_VFtab(%rbp),%rsi

        movapd (%rsi,%rax,8),%xmm4      ## Y1 F1        
        movapd %xmm4,%xmm5
        unpcklpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5

        movapd 16(%rsi,%rax,8),%xmm6    ## G1 H1        
        movapd %xmm6,%xmm7
        unpcklpd %xmm6,%xmm6
        unpckhpd %xmm7,%xmm7
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMM(%rsp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%rsp),%xmm5
        movsd %xmm5,nb314nf_vctot(%rsp)

_nb_kernel314nf_x86_64_sse2.nb314nf_updateouterdata: 
        ## get n from stack
        movl nb314nf_n(%rsp),%esi
        ## get group index for i particle 
        movq  nb314nf_gid(%rbp),%rdx            ## base of gid[]
        movl  (%rdx,%rsi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb314nf_vctot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb314nf_Vc(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## accumulate total lj energy and update it 
        movapd nb314nf_Vvdwtot(%rsp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movq  nb314nf_Vvdw(%rbp),%rax
        addsd (%rax,%rdx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%rax,%rdx,8)

        ## finish if last 
        movl nb314nf_nn1(%rsp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel314nf_x86_64_sse2.nb314nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb314nf_n(%rsp)
        jmp _nb_kernel314nf_x86_64_sse2.nb314nf_outer
_nb_kernel314nf_x86_64_sse2.nb314nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb314nf_nri(%rsp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel314nf_x86_64_sse2.nb314nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel314nf_x86_64_sse2.nb314nf_threadloop
_nb_kernel314nf_x86_64_sse2.nb314nf_end: 
        movl nb314nf_nouter(%rsp),%eax
        movl nb314nf_ninner(%rsp),%ebx
        movq nb314nf_outeriter(%rbp),%rcx
        movq nb314nf_inneriter(%rbp),%rdx
        movl %eax,(%rcx)
        movl %ebx,(%rdx)

        addq $968,%rsp
        emms


        pop %r15
        pop %r14
        pop %r13
        pop %r12

        pop %rbx
        pop    %rbp
        ret

