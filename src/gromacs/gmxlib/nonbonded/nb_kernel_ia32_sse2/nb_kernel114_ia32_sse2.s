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


.globl nb_kernel114_ia32_sse2
.globl _nb_kernel114_ia32_sse2
nb_kernel114_ia32_sse2: 
_nb_kernel114_ia32_sse2:        
.set nb114_p_nri, 8
.set nb114_iinr, 12
.set nb114_jindex, 16
.set nb114_jjnr, 20
.set nb114_shift, 24
.set nb114_shiftvec, 28
.set nb114_fshift, 32
.set nb114_gid, 36
.set nb114_pos, 40
.set nb114_faction, 44
.set nb114_charge, 48
.set nb114_p_facel, 52
.set nb114_argkrf, 56
.set nb114_argcrf, 60
.set nb114_Vc, 64
.set nb114_type, 68
.set nb114_p_ntype, 72
.set nb114_vdwparam, 76
.set nb114_Vvdw, 80
.set nb114_p_tabscale, 84
.set nb114_VFtab, 88
.set nb114_invsqrta, 92
.set nb114_dvda, 96
.set nb114_p_gbtabscale, 100
.set nb114_GBtab, 104
.set nb114_p_nthreads, 108
.set nb114_count, 112
.set nb114_mtx, 116
.set nb114_outeriter, 120
.set nb114_inneriter, 124
.set nb114_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb114_ixO, 0
.set nb114_iyO, 16
.set nb114_izO, 32
.set nb114_ixH1, 48
.set nb114_iyH1, 64
.set nb114_izH1, 80
.set nb114_ixH2, 96
.set nb114_iyH2, 112
.set nb114_izH2, 128
.set nb114_ixM, 144
.set nb114_iyM, 160
.set nb114_izM, 176
.set nb114_jxO, 192
.set nb114_jyO, 208
.set nb114_jzO, 224
.set nb114_jxH1, 240
.set nb114_jyH1, 256
.set nb114_jzH1, 272
.set nb114_jxH2, 288
.set nb114_jyH2, 304
.set nb114_jzH2, 320
.set nb114_jxM, 336
.set nb114_jyM, 352
.set nb114_jzM, 368
.set nb114_dxOO, 384
.set nb114_dyOO, 400
.set nb114_dzOO, 416
.set nb114_dxH1H1, 432
.set nb114_dyH1H1, 448
.set nb114_dzH1H1, 464
.set nb114_dxH1H2, 480
.set nb114_dyH1H2, 496
.set nb114_dzH1H2, 512
.set nb114_dxH1M, 528
.set nb114_dyH1M, 544
.set nb114_dzH1M, 560
.set nb114_dxH2H1, 576
.set nb114_dyH2H1, 592
.set nb114_dzH2H1, 608
.set nb114_dxH2H2, 624
.set nb114_dyH2H2, 640
.set nb114_dzH2H2, 656
.set nb114_dxH2M, 672
.set nb114_dyH2M, 688
.set nb114_dzH2M, 704
.set nb114_dxMH1, 720
.set nb114_dyMH1, 736
.set nb114_dzMH1, 752
.set nb114_dxMH2, 768
.set nb114_dyMH2, 784
.set nb114_dzMH2, 800
.set nb114_dxMM, 816
.set nb114_dyMM, 832
.set nb114_dzMM, 848
.set nb114_qqMM, 864
.set nb114_qqMH, 880
.set nb114_qqHH, 896
.set nb114_two, 912
.set nb114_c6, 944
.set nb114_c12, 960
.set nb114_vctot, 976
.set nb114_Vvdwtot, 992
.set nb114_fixO, 1008
.set nb114_fiyO, 1024
.set nb114_fizO, 1040
.set nb114_fixH1, 1056
.set nb114_fiyH1, 1072
.set nb114_fizH1, 1088
.set nb114_fixH2, 1104
.set nb114_fiyH2, 1120
.set nb114_fizH2, 1136
.set nb114_fixM, 1152
.set nb114_fiyM, 1168
.set nb114_fizM, 1184
.set nb114_fjxO, 1200
.set nb114_fjyO, 1216
.set nb114_fjzO, 1232
.set nb114_fjxH1, 1248
.set nb114_fjyH1, 1264
.set nb114_fjzH1, 1280
.set nb114_fjxH2, 1296
.set nb114_fjyH2, 1312
.set nb114_fjzH2, 1328
.set nb114_fjxM, 1344
.set nb114_fjyM, 1360
.set nb114_fjzM, 1376
.set nb114_half, 1392
.set nb114_three, 1408
.set nb114_six, 1424
.set nb114_twelve, 1440
.set nb114_rsqOO, 1456
.set nb114_rsqH1H1, 1472
.set nb114_rsqH1H2, 1488
.set nb114_rsqH1M, 1504
.set nb114_rsqH2H1, 1520
.set nb114_rsqH2H2, 1536
.set nb114_rsqH2M, 1552
.set nb114_rsqMH1, 1568
.set nb114_rsqMH2, 1584
.set nb114_rsqMM, 1600
.set nb114_rinvsqOO, 1616
.set nb114_rinvH1H1, 1632
.set nb114_rinvH1H2, 1648
.set nb114_rinvH1M, 1664
.set nb114_rinvH2H1, 1680
.set nb114_rinvH2H2, 1696
.set nb114_rinvH2M, 1712
.set nb114_rinvMH1, 1728
.set nb114_rinvMH2, 1744
.set nb114_rinvMM, 1760
.set nb114_is3, 1776
.set nb114_ii3, 1780
.set nb114_innerjjnr, 1784
.set nb114_innerk, 1788
.set nb114_n, 1792
.set nb114_nn1, 1796
.set nb114_nri, 1800
.set nb114_nouter, 1804
.set nb114_ninner, 1808
.set nb114_salign, 1812
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $1816,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb114_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb114_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb114_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb114_nouter(%esp)
        movl %eax,nb114_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb114_half(%esp)
        movl %ebx,nb114_half+4(%esp)
        movsd nb114_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm3,%xmm4
        addpd  %xmm4,%xmm4      ## 6.0
        movapd %xmm4,%xmm5
        addpd  %xmm5,%xmm5      ## 12.0
        movapd %xmm1,nb114_half(%esp)
        movapd %xmm2,nb114_two(%esp)
        movapd %xmm3,nb114_three(%esp)
        movapd %xmm4,nb114_six(%esp)
        movapd %xmm5,nb114_twelve(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb114_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb114_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb114_p_facel(%ebp),%esi
        movsd (%esi),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb114_qqMM(%esp)
        movapd %xmm4,nb114_qqMH(%esp)
        movapd %xmm5,nb114_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb114_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb114_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb114_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb114_c6(%esp)
        movapd %xmm1,nb114_c12(%esp)

_nb_kernel114_ia32_sse2.nb114_threadloop: 
        movl  nb114_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel114_ia32_sse2.nb114_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel114_ia32_sse2.nb114_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb114_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb114_n(%esp)
        movl %ebx,nb114_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel114_ia32_sse2.nb114_outerstart
        jmp _nb_kernel114_ia32_sse2.nb114_end

_nb_kernel114_ia32_sse2.nb114_outerstart: 
        ## ebx contains number of outer iterations
        addl nb114_nouter(%esp),%ebx
        movl %ebx,nb114_nouter(%esp)

_nb_kernel114_ia32_sse2.nb114_outer: 
        movl  nb114_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb114_is3(%esp)      ## store is3 

        movl  nb114_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb114_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb114_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb114_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb114_ixO(%esp)
        movapd %xmm4,nb114_iyO(%esp)
        movapd %xmm5,nb114_izO(%esp)
        movapd %xmm6,nb114_ixH1(%esp)
        movapd %xmm7,nb114_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb114_izH1(%esp)
        movapd %xmm0,nb114_ixH2(%esp)
        movapd %xmm1,nb114_iyH2(%esp)
        movapd %xmm2,nb114_izH2(%esp)
        movapd %xmm3,nb114_ixM(%esp)
        movapd %xmm4,nb114_iyM(%esp)
        movapd %xmm5,nb114_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb114_vctot(%esp)
        movapd %xmm4,nb114_Vvdwtot(%esp)
        movapd %xmm4,nb114_fixO(%esp)
        movapd %xmm4,nb114_fiyO(%esp)
        movapd %xmm4,nb114_fizO(%esp)
        movapd %xmm4,nb114_fixH1(%esp)
        movapd %xmm4,nb114_fiyH1(%esp)
        movapd %xmm4,nb114_fizH1(%esp)
        movapd %xmm4,nb114_fixH2(%esp)
        movapd %xmm4,nb114_fiyH2(%esp)
        movapd %xmm4,nb114_fizH2(%esp)
        movapd %xmm4,nb114_fixM(%esp)
        movapd %xmm4,nb114_fiyM(%esp)
        movapd %xmm4,nb114_fizM(%esp)

        movl  nb114_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb114_pos(%ebp),%esi
        movl  nb114_faction(%ebp),%edi
        movl  nb114_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb114_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb114_ninner(%esp),%ecx
        movl  %ecx,nb114_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb114_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel114_ia32_sse2.nb114_unroll_loop
        jmp   _nb_kernel114_ia32_sse2.nb114_checksingle
_nb_kernel114_ia32_sse2.nb114_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb114_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb114_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb114_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movlpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%eax,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm2
        movlpd 16(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%ebx,8),%xmm5
        movhpd 24(%esi,%eax,8),%xmm3
        movhpd 24(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb114_jxO(%esp)
        movapd  %xmm1,nb114_jyO(%esp)
        movapd  %xmm3,nb114_jzO(%esp)
        movapd  %xmm4,nb114_jxH1(%esp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%esi,%eax,8),%xmm0
        movlpd 32(%esi,%ebx,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm0
        movhpd 40(%esi,%ebx,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movlpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb114_jyH1(%esp)
        movapd  %xmm1,nb114_jzH1(%esp)
        movapd  %xmm3,nb114_jxH2(%esp)
        movapd  %xmm4,nb114_jyH2(%esp)

        ## load h2z, mx, my, mz
        movlpd 64(%esi,%eax,8),%xmm0
        movlpd 64(%esi,%ebx,8),%xmm2
        movhpd 72(%esi,%eax,8),%xmm0
        movhpd 72(%esi,%ebx,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 80(%esi,%ebx,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb114_jzH2(%esp)
        movapd  %xmm1,nb114_jxM(%esp)
        movapd  %xmm3,nb114_jyM(%esp)
        movapd  %xmm4,nb114_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb114_ixO(%esp),%xmm0
        movapd nb114_iyO(%esp),%xmm1
        movapd nb114_izO(%esp),%xmm2
        movapd nb114_ixH1(%esp),%xmm3
        movapd nb114_iyH1(%esp),%xmm4
        movapd nb114_izH1(%esp),%xmm5
        subpd  nb114_jxO(%esp),%xmm0
        subpd  nb114_jyO(%esp),%xmm1
        subpd  nb114_jzO(%esp),%xmm2
        subpd  nb114_jxH1(%esp),%xmm3
        subpd  nb114_jyH1(%esp),%xmm4
        subpd  nb114_jzH1(%esp),%xmm5
        movapd %xmm0,nb114_dxOO(%esp)
        movapd %xmm1,nb114_dyOO(%esp)
        movapd %xmm2,nb114_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxH1H1(%esp)
        movapd %xmm4,nb114_dyH1H1(%esp)
        movapd %xmm5,nb114_dzH1H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb114_rsqOO(%esp)
        movapd %xmm3,nb114_rsqH1H1(%esp)

        movapd nb114_ixH1(%esp),%xmm0
        movapd nb114_iyH1(%esp),%xmm1
        movapd nb114_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb114_jxH2(%esp),%xmm0
        subpd  nb114_jyH2(%esp),%xmm1
        subpd  nb114_jzH2(%esp),%xmm2
        subpd  nb114_jxM(%esp),%xmm3
        subpd  nb114_jyM(%esp),%xmm4
        subpd  nb114_jzM(%esp),%xmm5
        movapd %xmm0,nb114_dxH1H2(%esp)
        movapd %xmm1,nb114_dyH1H2(%esp)
        movapd %xmm2,nb114_dzH1H2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxH1M(%esp)
        movapd %xmm4,nb114_dyH1M(%esp)
        movapd %xmm5,nb114_dzH1M(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb114_rsqH1H2(%esp)
        movapd %xmm3,nb114_rsqH1M(%esp)

        movapd nb114_ixH2(%esp),%xmm0
        movapd nb114_iyH2(%esp),%xmm1
        movapd nb114_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb114_jxH1(%esp),%xmm0
        subpd  nb114_jyH1(%esp),%xmm1
        subpd  nb114_jzH1(%esp),%xmm2
        subpd  nb114_jxH2(%esp),%xmm3
        subpd  nb114_jyH2(%esp),%xmm4
        subpd  nb114_jzH2(%esp),%xmm5
        movapd %xmm0,nb114_dxH2H1(%esp)
        movapd %xmm1,nb114_dyH2H1(%esp)
        movapd %xmm2,nb114_dzH2H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxH2H2(%esp)
        movapd %xmm4,nb114_dyH2H2(%esp)
        movapd %xmm5,nb114_dzH2H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb114_rsqH2H1(%esp)
        movapd %xmm3,nb114_rsqH2H2(%esp)

        movapd nb114_ixH2(%esp),%xmm0
        movapd nb114_iyH2(%esp),%xmm1
        movapd nb114_izH2(%esp),%xmm2
        movapd nb114_ixM(%esp),%xmm3
        movapd nb114_iyM(%esp),%xmm4
        movapd nb114_izM(%esp),%xmm5
        subpd  nb114_jxM(%esp),%xmm0
        subpd  nb114_jyM(%esp),%xmm1
        subpd  nb114_jzM(%esp),%xmm2
        subpd  nb114_jxH1(%esp),%xmm3
        subpd  nb114_jyH1(%esp),%xmm4
        subpd  nb114_jzH1(%esp),%xmm5
        movapd %xmm0,nb114_dxH2M(%esp)
        movapd %xmm1,nb114_dyH2M(%esp)
        movapd %xmm2,nb114_dzH2M(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxMH1(%esp)
        movapd %xmm4,nb114_dyMH1(%esp)
        movapd %xmm5,nb114_dzMH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb114_rsqH2M(%esp)
        movapd %xmm4,nb114_rsqMH1(%esp)

        movapd nb114_ixM(%esp),%xmm0
        movapd nb114_iyM(%esp),%xmm1
        movapd nb114_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb114_jxH2(%esp),%xmm0
        subpd  nb114_jyH2(%esp),%xmm1
        subpd  nb114_jzH2(%esp),%xmm2
        subpd  nb114_jxM(%esp),%xmm3
        subpd  nb114_jyM(%esp),%xmm4
        subpd  nb114_jzM(%esp),%xmm5
        movapd %xmm0,nb114_dxMH2(%esp)
        movapd %xmm1,nb114_dyMH2(%esp)
        movapd %xmm2,nb114_dzMH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxMM(%esp)
        movapd %xmm4,nb114_dyMM(%esp)
        movapd %xmm5,nb114_dzMM(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb114_rsqMH2(%esp)
        movapd %xmm4,nb114_rsqMM(%esp)

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
        movapd  nb114_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm3   ## iter1 
        mulpd   nb114_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm1   ## rinv 
        mulpd   nb114_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114_rinvMH2(%esp)
        movapd %xmm5,nb114_rinvMM(%esp)

        movapd nb114_rsqOO(%esp),%xmm0
        movapd nb114_rsqH1H1(%esp),%xmm4
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
        movapd  nb114_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb114_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm1   ## rinv 
        mulpd   nb114_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb114_rinvsqOO(%esp)
        movapd %xmm5,nb114_rinvH1H1(%esp)

        movapd nb114_rsqH1H2(%esp),%xmm0
        movapd nb114_rsqH1M(%esp),%xmm4
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
        movapd  nb114_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm3   ## iter1 
        mulpd   nb114_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm1   ## rinv 
        mulpd   nb114_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114_rinvH1H2(%esp)
        movapd %xmm5,nb114_rinvH1M(%esp)

        movapd nb114_rsqH2H1(%esp),%xmm0
        movapd nb114_rsqH2H2(%esp),%xmm4
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
        movapd  nb114_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm3   ## iter1a 
        mulpd   nb114_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm1   ## rinv 
        mulpd   nb114_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114_rinvH2H1(%esp)
        movapd %xmm5,nb114_rinvH2H2(%esp)

        movapd nb114_rsqMH1(%esp),%xmm0
        movapd nb114_rsqH2M(%esp),%xmm4
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
        movapd  nb114_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm3   ## iter1a 
        mulpd   nb114_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114_half(%esp),%xmm1   ## rinv 
        mulpd   nb114_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114_rinvMH1(%esp)
        movapd %xmm5,nb114_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb114_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb114_c6(%esp),%xmm1
        mulpd  nb114_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb114_Vvdwtot(%esp),%xmm3
        mulpd  nb114_six(%esp),%xmm1
        mulpd  nb114_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm0,%xmm2
        movapd %xmm3,nb114_Vvdwtot(%esp)

        movapd %xmm2,%xmm0
        movapd %xmm2,%xmm1 ## fscal

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb114_dxOO(%esp),%xmm0
        mulpd nb114_dyOO(%esp),%xmm1
        mulpd nb114_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixO(%esp),%xmm0
        addpd nb114_fiyO(%esp),%xmm1
        addpd nb114_fizO(%esp),%xmm2
        movapd %xmm3,nb114_fjxO(%esp)
        movapd %xmm4,nb114_fjyO(%esp)
        movapd %xmm5,nb114_fjzO(%esp)
        movapd %xmm0,nb114_fixO(%esp)
        movapd %xmm1,nb114_fiyO(%esp)
        movapd %xmm2,nb114_fizO(%esp)

        ## H1-H1 interaction 
        movapd nb114_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7              ## xmm7=rinv 
        mulpd  %xmm0,%xmm0              ## xmm0=rinvsq 
        mulpd  nb114_qqHH(%esp),%xmm7   ## xmm7=vcoul
        mulpd  %xmm7,%xmm0              ## fscal
        addpd  nb114_vctot(%esp),%xmm7   ## local vctot summation variable 
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb114_dxH1H1(%esp),%xmm0
        mulpd nb114_dyH1H1(%esp),%xmm1
        mulpd nb114_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixH1(%esp),%xmm0
        addpd nb114_fiyH1(%esp),%xmm1
        addpd nb114_fizH1(%esp),%xmm2
        movapd %xmm3,nb114_fjxH1(%esp)
        movapd %xmm4,nb114_fjyH1(%esp)
        movapd %xmm5,nb114_fjzH1(%esp)
        movapd %xmm0,nb114_fixH1(%esp)
        movapd %xmm1,nb114_fiyH1(%esp)
        movapd %xmm2,nb114_fizH1(%esp)

        ## H1-H2 interaction  
        movapd nb114_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulpd  nb114_qqHH(%esp),%xmm7   ## vcoul
        mulpd  %xmm7,%xmm0
        addpd  nb114_vctot(%esp),%xmm7   ## local vctot summation variable 
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb114_dxH1H2(%esp),%xmm0
        mulpd nb114_dyH1H2(%esp),%xmm1
        mulpd nb114_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixH1(%esp),%xmm0
        addpd nb114_fiyH1(%esp),%xmm1
        addpd nb114_fizH1(%esp),%xmm2
        movapd %xmm3,nb114_fjxH2(%esp)
        movapd %xmm4,nb114_fjyH2(%esp)
        movapd %xmm5,nb114_fjzH2(%esp)
        movapd %xmm0,nb114_fixH1(%esp)
        movapd %xmm1,nb114_fiyH1(%esp)
        movapd %xmm2,nb114_fizH1(%esp)

        ## H1-M interaction 
        movapd nb114_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulpd  nb114_qqMH(%esp),%xmm7   ## vcoul
        mulpd  %xmm7,%xmm0
        addpd  nb114_vctot(%esp),%xmm7   ## local vctot summation variable 
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb114_dxH1M(%esp),%xmm0
        mulpd nb114_dyH1M(%esp),%xmm1
        mulpd nb114_dzH1M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixH1(%esp),%xmm0
        addpd nb114_fiyH1(%esp),%xmm1
        addpd nb114_fizH1(%esp),%xmm2
        movapd %xmm3,nb114_fjxM(%esp)
        movapd %xmm4,nb114_fjyM(%esp)
        movapd %xmm5,nb114_fjzM(%esp)
        movapd %xmm0,nb114_fixH1(%esp)
        movapd %xmm1,nb114_fiyH1(%esp)
        movapd %xmm2,nb114_fizH1(%esp)

        ## H2-H1 interaction 
        movapd nb114_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulpd  nb114_qqHH(%esp),%xmm7
        mulpd  %xmm7,%xmm0
        addpd  nb114_vctot(%esp),%xmm7   ## local vctot summation variable 
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxH1(%esp),%xmm3
        movapd nb114_fjyH1(%esp),%xmm4
        movapd nb114_fjzH1(%esp),%xmm5
        mulpd nb114_dxH2H1(%esp),%xmm0
        mulpd nb114_dyH2H1(%esp),%xmm1
        mulpd nb114_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixH2(%esp),%xmm0
        addpd nb114_fiyH2(%esp),%xmm1
        addpd nb114_fizH2(%esp),%xmm2
        movapd %xmm3,nb114_fjxH1(%esp)
        movapd %xmm4,nb114_fjyH1(%esp)
        movapd %xmm5,nb114_fjzH1(%esp)
        movapd %xmm0,nb114_fixH2(%esp)
        movapd %xmm1,nb114_fiyH2(%esp)
        movapd %xmm2,nb114_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb114_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulpd  nb114_qqHH(%esp),%xmm7   ## vcoul

        mulpd  %xmm7,%xmm0
        addpd  nb114_vctot(%esp),%xmm7   ## local vctot summation variable 
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxH2(%esp),%xmm3
        movapd nb114_fjyH2(%esp),%xmm4
        movapd nb114_fjzH2(%esp),%xmm5
        mulpd nb114_dxH2H2(%esp),%xmm0
        mulpd nb114_dyH2H2(%esp),%xmm1
        mulpd nb114_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixH2(%esp),%xmm0
        addpd nb114_fiyH2(%esp),%xmm1
        addpd nb114_fizH2(%esp),%xmm2
        movapd %xmm3,nb114_fjxH2(%esp)
        movapd %xmm4,nb114_fjyH2(%esp)
        movapd %xmm5,nb114_fjzH2(%esp)
        movapd %xmm0,nb114_fixH2(%esp)
        movapd %xmm1,nb114_fiyH2(%esp)
        movapd %xmm2,nb114_fizH2(%esp)

        ## H2-M interaction 
        movapd nb114_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulpd  nb114_qqMH(%esp),%xmm7   ## vcoul

        mulpd  %xmm7,%xmm0
        addpd  nb114_vctot(%esp),%xmm7
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxM(%esp),%xmm3
        movapd nb114_fjyM(%esp),%xmm4
        movapd nb114_fjzM(%esp),%xmm5
        mulpd nb114_dxH2M(%esp),%xmm0
        mulpd nb114_dyH2M(%esp),%xmm1
        mulpd nb114_dzH2M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixH2(%esp),%xmm0
        addpd nb114_fiyH2(%esp),%xmm1
        addpd nb114_fizH2(%esp),%xmm2
        movapd %xmm3,nb114_fjxM(%esp)
        movapd %xmm4,nb114_fjyM(%esp)
        movapd %xmm5,nb114_fjzM(%esp)
        movapd %xmm0,nb114_fixH2(%esp)
        movapd %xmm1,nb114_fiyH2(%esp)
        movapd %xmm2,nb114_fizH2(%esp)

        ## M-H1 interaction 
        movapd nb114_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulpd  nb114_qqMH(%esp),%xmm7

        mulpd  %xmm7,%xmm0
        addpd  nb114_vctot(%esp),%xmm7   ## local vctot summation variable 
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxH1(%esp),%xmm3
        movapd nb114_fjyH1(%esp),%xmm4
        movapd nb114_fjzH1(%esp),%xmm5
        mulpd nb114_dxMH1(%esp),%xmm0
        mulpd nb114_dyMH1(%esp),%xmm1
        mulpd nb114_dzMH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixM(%esp),%xmm0
        addpd nb114_fiyM(%esp),%xmm1
        addpd nb114_fizM(%esp),%xmm2
        movapd %xmm3,nb114_fjxH1(%esp)
        movapd %xmm4,nb114_fjyH1(%esp)
        movapd %xmm5,nb114_fjzH1(%esp)
        movapd %xmm0,nb114_fixM(%esp)
        movapd %xmm1,nb114_fiyM(%esp)
        movapd %xmm2,nb114_fizM(%esp)

        ## M-H2 interaction 
        movapd nb114_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulpd  nb114_qqMH(%esp),%xmm7

        mulpd  %xmm7,%xmm0
        addpd  nb114_vctot(%esp),%xmm7
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxH2(%esp),%xmm3
        movapd nb114_fjyH2(%esp),%xmm4
        movapd nb114_fjzH2(%esp),%xmm5
        mulpd nb114_dxMH2(%esp),%xmm0
        mulpd nb114_dyMH2(%esp),%xmm1
        mulpd nb114_dzMH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixM(%esp),%xmm0
        addpd nb114_fiyM(%esp),%xmm1
        addpd nb114_fizM(%esp),%xmm2
        movapd %xmm3,nb114_fjxH2(%esp)
        movapd %xmm4,nb114_fjyH2(%esp)
        movapd %xmm5,nb114_fjzH2(%esp)
        movapd %xmm0,nb114_fixM(%esp)
        movapd %xmm1,nb114_fiyM(%esp)
        movapd %xmm2,nb114_fizM(%esp)

        ## M-M interaction 
        movapd nb114_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        mulpd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulpd  nb114_qqMM(%esp),%xmm7   ## vcoul

        mulpd  %xmm7,%xmm0
        addpd  nb114_vctot(%esp),%xmm7   ## local vctot summation variable 
        movapd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxM(%esp),%xmm3
        movapd nb114_fjyM(%esp),%xmm4
        movapd nb114_fjzM(%esp),%xmm5
        mulpd nb114_dxMM(%esp),%xmm0
        mulpd nb114_dyMM(%esp),%xmm1
        mulpd nb114_dzMM(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb114_fixM(%esp),%xmm0
        addpd nb114_fiyM(%esp),%xmm1
        addpd nb114_fizM(%esp),%xmm2
        movapd %xmm3,nb114_fjxM(%esp)
        movapd %xmm4,nb114_fjyM(%esp)
        movapd %xmm5,nb114_fjzM(%esp)
        movapd %xmm0,nb114_fixM(%esp)
        movapd %xmm1,nb114_fiyM(%esp)
        movapd %xmm2,nb114_fizM(%esp)

        movl nb114_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## Step1 - transpose fjxO, fjyO and fjzO, fjxH1
        movapd nb114_fjxO(%esp),%xmm0
        movapd nb114_fjyO(%esp),%xmm1
        movapd nb114_fjzO(%esp),%xmm2
        movapd nb114_fjxH1(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd (%edi,%eax,8),%xmm1
        movlpd (%edi,%ebx,8),%xmm3
        movhpd 8(%edi,%eax,8),%xmm1
        movhpd 8(%edi,%ebx,8),%xmm3
        movlpd 16(%edi,%eax,8),%xmm6
        movlpd 16(%edi,%ebx,8),%xmm7
        movhpd 24(%edi,%eax,8),%xmm6
        movhpd 24(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,(%edi,%eax,8)
        movlpd %xmm3,(%edi,%ebx,8)
        movhpd %xmm1,8(%edi,%eax,8)
        movhpd %xmm3,8(%edi,%ebx,8)
        movlpd %xmm6,16(%edi,%eax,8)
        movlpd %xmm7,16(%edi,%ebx,8)
        movhpd %xmm6,24(%edi,%eax,8)
        movhpd %xmm7,24(%edi,%ebx,8)

        ## Step2 - transpose fjyH1, fjzH1 and fjxH2, fjyH2
        movapd nb114_fjyH1(%esp),%xmm0
        movapd nb114_fjzH1(%esp),%xmm1
        movapd nb114_fjxH2(%esp),%xmm2
        movapd nb114_fjyH2(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd 32(%edi,%eax,8),%xmm1
        movlpd 32(%edi,%ebx,8),%xmm3
        movhpd 40(%edi,%eax,8),%xmm1
        movhpd 40(%edi,%ebx,8),%xmm3
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 48(%edi,%ebx,8),%xmm7
        movhpd 56(%edi,%eax,8),%xmm6
        movhpd 56(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,32(%edi,%eax,8)
        movlpd %xmm3,32(%edi,%ebx,8)
        movhpd %xmm1,40(%edi,%eax,8)
        movhpd %xmm3,40(%edi,%ebx,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,48(%edi,%ebx,8)
        movhpd %xmm6,56(%edi,%eax,8)
        movhpd %xmm7,56(%edi,%ebx,8)

        ## Step3 - transpose fjzH2, fjxM and fjyM, fjzM
        movapd nb114_fjzH2(%esp),%xmm0
        movapd nb114_fjxM(%esp),%xmm1
        movapd nb114_fjyM(%esp),%xmm2
        movapd nb114_fjzM(%esp),%xmm3
        movapd %xmm0,%xmm4
        movapd %xmm2,%xmm5
        unpcklpd %xmm1,%xmm0  ## fjOxA fjOyA
        unpckhpd %xmm1,%xmm4  ## fjOxB fjOyB
        unpcklpd %xmm3,%xmm2  ## fjOzA fjH1xA
        unpckhpd %xmm3,%xmm5  ## fjOzB fjH1xB
        movlpd 64(%edi,%eax,8),%xmm1
        movlpd 64(%edi,%ebx,8),%xmm3
        movhpd 72(%edi,%eax,8),%xmm1
        movhpd 72(%edi,%ebx,8),%xmm3
        movlpd 80(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%ebx,8),%xmm7
        movhpd 88(%edi,%eax,8),%xmm6
        movhpd 88(%edi,%ebx,8),%xmm7
        addpd  %xmm0,%xmm1
        addpd  %xmm4,%xmm3
        addpd  %xmm2,%xmm6
        addpd  %xmm5,%xmm7
        movlpd %xmm1,64(%edi,%eax,8)
        movlpd %xmm3,64(%edi,%ebx,8)
        movhpd %xmm1,72(%edi,%eax,8)
        movhpd %xmm3,72(%edi,%ebx,8)
        movlpd %xmm6,80(%edi,%eax,8)
        movlpd %xmm7,80(%edi,%ebx,8)
        movhpd %xmm6,88(%edi,%eax,8)
        movhpd %xmm7,88(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb114_innerk(%esp)
        jl    _nb_kernel114_ia32_sse2.nb114_checksingle
        jmp   _nb_kernel114_ia32_sse2.nb114_unroll_loop
_nb_kernel114_ia32_sse2.nb114_checksingle: 
        movl  nb114_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel114_ia32_sse2.nb114_dosingle
        jmp   _nb_kernel114_ia32_sse2.nb114_updateouterdata
_nb_kernel114_ia32_sse2.nb114_dosingle: 
        movl  nb114_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb114_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movhpd 8(%esi,%eax,8),%xmm0
        movlpd 16(%esi,%eax,8),%xmm1
        movhpd 24(%esi,%eax,8),%xmm1
        movlpd 32(%esi,%eax,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%eax,8),%xmm4
        movlpd 80(%esi,%eax,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm5
        movsd  %xmm0,nb114_jxO(%esp)
        movsd  %xmm1,nb114_jzO(%esp)
        movsd  %xmm2,nb114_jyH1(%esp)
        movsd  %xmm3,nb114_jxH2(%esp)
        movsd  %xmm4,nb114_jzH2(%esp)
        movsd  %xmm5,nb114_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb114_jyO(%esp)
        movsd  %xmm1,nb114_jxH1(%esp)
        movsd  %xmm2,nb114_jzH1(%esp)
        movsd  %xmm3,nb114_jyH2(%esp)
        movsd  %xmm4,nb114_jxM(%esp)
        movsd  %xmm5,nb114_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb114_ixO(%esp),%xmm0
        movapd nb114_iyO(%esp),%xmm1
        movapd nb114_izO(%esp),%xmm2
        movapd nb114_ixH1(%esp),%xmm3
        movapd nb114_iyH1(%esp),%xmm4
        movapd nb114_izH1(%esp),%xmm5
        subsd  nb114_jxO(%esp),%xmm0
        subsd  nb114_jyO(%esp),%xmm1
        subsd  nb114_jzO(%esp),%xmm2
        subsd  nb114_jxH1(%esp),%xmm3
        subsd  nb114_jyH1(%esp),%xmm4
        subsd  nb114_jzH1(%esp),%xmm5
        movapd %xmm0,nb114_dxOO(%esp)
        movapd %xmm1,nb114_dyOO(%esp)
        movapd %xmm2,nb114_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxH1H1(%esp)
        movapd %xmm4,nb114_dyH1H1(%esp)
        movapd %xmm5,nb114_dzH1H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb114_rsqOO(%esp)
        movapd %xmm3,nb114_rsqH1H1(%esp)

        movapd nb114_ixH1(%esp),%xmm0
        movapd nb114_iyH1(%esp),%xmm1
        movapd nb114_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb114_jxH2(%esp),%xmm0
        subsd  nb114_jyH2(%esp),%xmm1
        subsd  nb114_jzH2(%esp),%xmm2
        subsd  nb114_jxM(%esp),%xmm3
        subsd  nb114_jyM(%esp),%xmm4
        subsd  nb114_jzM(%esp),%xmm5
        movapd %xmm0,nb114_dxH1H2(%esp)
        movapd %xmm1,nb114_dyH1H2(%esp)
        movapd %xmm2,nb114_dzH1H2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxH1M(%esp)
        movapd %xmm4,nb114_dyH1M(%esp)
        movapd %xmm5,nb114_dzH1M(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb114_rsqH1H2(%esp)
        movapd %xmm3,nb114_rsqH1M(%esp)

        movapd nb114_ixH2(%esp),%xmm0
        movapd nb114_iyH2(%esp),%xmm1
        movapd nb114_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb114_jxH1(%esp),%xmm0
        subsd  nb114_jyH1(%esp),%xmm1
        subsd  nb114_jzH1(%esp),%xmm2
        subsd  nb114_jxH2(%esp),%xmm3
        subsd  nb114_jyH2(%esp),%xmm4
        subsd  nb114_jzH2(%esp),%xmm5
        movapd %xmm0,nb114_dxH2H1(%esp)
        movapd %xmm1,nb114_dyH2H1(%esp)
        movapd %xmm2,nb114_dzH2H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxH2H2(%esp)
        movapd %xmm4,nb114_dyH2H2(%esp)
        movapd %xmm5,nb114_dzH2H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb114_rsqH2H1(%esp)
        movapd %xmm3,nb114_rsqH2H2(%esp)

        movapd nb114_ixH2(%esp),%xmm0
        movapd nb114_iyH2(%esp),%xmm1
        movapd nb114_izH2(%esp),%xmm2
        movapd nb114_ixM(%esp),%xmm3
        movapd nb114_iyM(%esp),%xmm4
        movapd nb114_izM(%esp),%xmm5
        subsd  nb114_jxM(%esp),%xmm0
        subsd  nb114_jyM(%esp),%xmm1
        subsd  nb114_jzM(%esp),%xmm2
        subsd  nb114_jxH1(%esp),%xmm3
        subsd  nb114_jyH1(%esp),%xmm4
        subsd  nb114_jzH1(%esp),%xmm5
        movapd %xmm0,nb114_dxH2M(%esp)
        movapd %xmm1,nb114_dyH2M(%esp)
        movapd %xmm2,nb114_dzH2M(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxMH1(%esp)
        movapd %xmm4,nb114_dyMH1(%esp)
        movapd %xmm5,nb114_dzMH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb114_rsqH2M(%esp)
        movapd %xmm4,nb114_rsqMH1(%esp)

        movapd nb114_ixM(%esp),%xmm0
        movapd nb114_iyM(%esp),%xmm1
        movapd nb114_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb114_jxH2(%esp),%xmm0
        subsd  nb114_jyH2(%esp),%xmm1
        subsd  nb114_jzH2(%esp),%xmm2
        subsd  nb114_jxM(%esp),%xmm3
        subsd  nb114_jyM(%esp),%xmm4
        subsd  nb114_jzM(%esp),%xmm5
        movapd %xmm0,nb114_dxMH2(%esp)
        movapd %xmm1,nb114_dyMH2(%esp)
        movapd %xmm2,nb114_dzMH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb114_dxMM(%esp)
        movapd %xmm4,nb114_dyMM(%esp)
        movapd %xmm5,nb114_dzMM(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb114_rsqMH2(%esp)
        movapd %xmm4,nb114_rsqMM(%esp)

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
        movapd  nb114_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm3   ## iter1 
        mulsd   nb114_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm1   ## rinv 
        mulsd   nb114_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114_rinvMH2(%esp)
        movapd %xmm5,nb114_rinvMM(%esp)

        movapd nb114_rsqOO(%esp),%xmm0
        movapd nb114_rsqH1H1(%esp),%xmm4
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
        movapd  nb114_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb114_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm1   ## rinv 
        mulsd   nb114_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb114_rinvsqOO(%esp)
        movapd %xmm5,nb114_rinvH1H1(%esp)

        movapd nb114_rsqH1H2(%esp),%xmm0
        movapd nb114_rsqH1M(%esp),%xmm4
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
        movapd  nb114_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm3   ## iter1 
        mulsd   nb114_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm1   ## rinv 
        mulsd   nb114_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114_rinvH1H2(%esp)
        movapd %xmm5,nb114_rinvH1M(%esp)

        movapd nb114_rsqH2H1(%esp),%xmm0
        movapd nb114_rsqH2H2(%esp),%xmm4
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
        movapd  nb114_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm3   ## iter1a 
        mulsd   nb114_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm1   ## rinv 
        mulsd   nb114_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114_rinvH2H1(%esp)
        movapd %xmm5,nb114_rinvH2H2(%esp)

        movapd nb114_rsqMH1(%esp),%xmm0
        movapd nb114_rsqH2M(%esp),%xmm4
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
        movapd  nb114_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm3   ## iter1a 
        mulsd   nb114_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114_half(%esp),%xmm1   ## rinv 
        mulsd   nb114_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114_rinvMH1(%esp)
        movapd %xmm5,nb114_rinvH2M(%esp)

        ## start with OO interaction 
        movsd nb114_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movsd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb114_c6(%esp),%xmm1
        mulsd  nb114_c12(%esp),%xmm2
        movsd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb114_Vvdwtot(%esp),%xmm3
        mulsd  nb114_six(%esp),%xmm1
        mulsd  nb114_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movsd %xmm3,nb114_Vvdwtot(%esp)

        movapd %xmm2,%xmm0
        movapd %xmm2,%xmm1


        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb114_dxOO(%esp),%xmm0
        mulsd nb114_dyOO(%esp),%xmm1
        mulsd nb114_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixO(%esp),%xmm0
        addsd nb114_fiyO(%esp),%xmm1
        addsd nb114_fizO(%esp),%xmm2
        movsd %xmm3,nb114_fjxO(%esp)
        movsd %xmm4,nb114_fjyO(%esp)
        movsd %xmm5,nb114_fjzO(%esp)
        movsd %xmm0,nb114_fixO(%esp)
        movsd %xmm1,nb114_fiyO(%esp)
        movsd %xmm2,nb114_fizO(%esp)

        ## H1-H1 interaction 
        movsd nb114_rinvH1H1(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulsd  nb114_qqHH(%esp),%xmm7

        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7
        movsd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb114_dxH1H1(%esp),%xmm0
        mulsd nb114_dyH1H1(%esp),%xmm1
        mulsd nb114_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixH1(%esp),%xmm0
        addsd nb114_fiyH1(%esp),%xmm1
        addsd nb114_fizH1(%esp),%xmm2
        movsd %xmm3,nb114_fjxH1(%esp)
        movsd %xmm4,nb114_fjyH1(%esp)
        movsd %xmm5,nb114_fjzH1(%esp)
        movsd %xmm0,nb114_fixH1(%esp)
        movsd %xmm1,nb114_fiyH1(%esp)
        movsd %xmm2,nb114_fizH1(%esp)

        ## H1-H2 interaction  
        movsd nb114_rinvH1H2(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulsd  nb114_qqHH(%esp),%xmm7   ## vcoul

        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7
        movsd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb114_dxH1H2(%esp),%xmm0
        mulsd nb114_dyH1H2(%esp),%xmm1
        mulsd nb114_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixH1(%esp),%xmm0
        addsd nb114_fiyH1(%esp),%xmm1
        addsd nb114_fizH1(%esp),%xmm2
        movsd %xmm3,nb114_fjxH2(%esp)
        movsd %xmm4,nb114_fjyH2(%esp)
        movsd %xmm5,nb114_fjzH2(%esp)
        movsd %xmm0,nb114_fixH1(%esp)
        movsd %xmm1,nb114_fiyH1(%esp)
        movsd %xmm2,nb114_fizH1(%esp)

        ## H1-M interaction 
        movsd nb114_rinvH1M(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulsd  nb114_qqMH(%esp),%xmm7   ## vcoul

        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7
        movsd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb114_dxH1M(%esp),%xmm0
        mulsd nb114_dyH1M(%esp),%xmm1
        mulsd nb114_dzH1M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixH1(%esp),%xmm0
        addsd nb114_fiyH1(%esp),%xmm1
        addsd nb114_fizH1(%esp),%xmm2
        movsd %xmm3,nb114_fjxM(%esp)
        movsd %xmm4,nb114_fjyM(%esp)
        movsd %xmm5,nb114_fjzM(%esp)
        movsd %xmm0,nb114_fixH1(%esp)
        movsd %xmm1,nb114_fiyH1(%esp)
        movsd %xmm2,nb114_fizH1(%esp)

        ## H2-H1 interaction 
        movsd nb114_rinvH2H1(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulsd  nb114_qqHH(%esp),%xmm7   ## vcoul

        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7
        movsd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxH1(%esp),%xmm3
        movapd nb114_fjyH1(%esp),%xmm4
        movapd nb114_fjzH1(%esp),%xmm5
        mulsd nb114_dxH2H1(%esp),%xmm0
        mulsd nb114_dyH2H1(%esp),%xmm1
        mulsd nb114_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixH2(%esp),%xmm0
        addsd nb114_fiyH2(%esp),%xmm1
        addsd nb114_fizH2(%esp),%xmm2
        movsd %xmm3,nb114_fjxH1(%esp)
        movsd %xmm4,nb114_fjyH1(%esp)
        movsd %xmm5,nb114_fjzH1(%esp)
        movsd %xmm0,nb114_fixH2(%esp)
        movsd %xmm1,nb114_fiyH2(%esp)
        movsd %xmm2,nb114_fizH2(%esp)

        ## H2-H2 interaction 
        movsd nb114_rinvH2H2(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulsd  nb114_qqHH(%esp),%xmm7   ## vcoul

        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7
        movsd %xmm7,nb114_vctot(%esp)

        movsd %xmm0,%xmm1
        movsd %xmm0,%xmm2

        movsd nb114_fjxH2(%esp),%xmm3
        movsd nb114_fjyH2(%esp),%xmm4
        movsd nb114_fjzH2(%esp),%xmm5
        mulsd nb114_dxH2H2(%esp),%xmm0
        mulsd nb114_dyH2H2(%esp),%xmm1
        mulsd nb114_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixH2(%esp),%xmm0
        addsd nb114_fiyH2(%esp),%xmm1
        addsd nb114_fizH2(%esp),%xmm2
        movsd %xmm3,nb114_fjxH2(%esp)
        movsd %xmm4,nb114_fjyH2(%esp)
        movsd %xmm5,nb114_fjzH2(%esp)
        movsd %xmm0,nb114_fixH2(%esp)
        movsd %xmm1,nb114_fiyH2(%esp)
        movsd %xmm2,nb114_fizH2(%esp)

        ## H2-M interaction 
        movsd nb114_rinvH2M(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 
        mulsd  nb114_qqMH(%esp),%xmm7
        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7
        movsd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxM(%esp),%xmm3
        movapd nb114_fjyM(%esp),%xmm4
        movapd nb114_fjzM(%esp),%xmm5
        mulsd nb114_dxH2M(%esp),%xmm0
        mulsd nb114_dyH2M(%esp),%xmm1
        mulsd nb114_dzH2M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixH2(%esp),%xmm0
        addsd nb114_fiyH2(%esp),%xmm1
        addsd nb114_fizH2(%esp),%xmm2
        movsd %xmm3,nb114_fjxM(%esp)
        movsd %xmm4,nb114_fjyM(%esp)
        movsd %xmm5,nb114_fjzM(%esp)
        movsd %xmm0,nb114_fixH2(%esp)
        movsd %xmm1,nb114_fiyH2(%esp)
        movsd %xmm2,nb114_fizH2(%esp)

        ## M-H1 interaction 
        movsd nb114_rinvMH1(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb114_qqMH(%esp),%xmm7

        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7   ## local vctot summation variable 
        movsd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxH1(%esp),%xmm3
        movapd nb114_fjyH1(%esp),%xmm4
        movapd nb114_fjzH1(%esp),%xmm5
        mulsd nb114_dxMH1(%esp),%xmm0
        mulsd nb114_dyMH1(%esp),%xmm1
        mulsd nb114_dzMH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixM(%esp),%xmm0
        addsd nb114_fiyM(%esp),%xmm1
        addsd nb114_fizM(%esp),%xmm2
        movsd %xmm3,nb114_fjxH1(%esp)
        movsd %xmm4,nb114_fjyH1(%esp)
        movsd %xmm5,nb114_fjzH1(%esp)
        movsd %xmm0,nb114_fixM(%esp)
        movsd %xmm1,nb114_fiyM(%esp)
        movsd %xmm2,nb114_fizM(%esp)

        ## M-H2 interaction 
        movsd nb114_rinvMH2(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb114_qqMH(%esp),%xmm7

        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7
        movsd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxH2(%esp),%xmm3
        movapd nb114_fjyH2(%esp),%xmm4
        movapd nb114_fjzH2(%esp),%xmm5
        mulsd nb114_dxMH2(%esp),%xmm0
        mulsd nb114_dyMH2(%esp),%xmm1
        mulsd nb114_dzMH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixM(%esp),%xmm0
        addsd nb114_fiyM(%esp),%xmm1
        addsd nb114_fizM(%esp),%xmm2
        movsd %xmm3,nb114_fjxH2(%esp)
        movsd %xmm4,nb114_fjyH2(%esp)
        movsd %xmm5,nb114_fjzH2(%esp)
        movsd %xmm0,nb114_fixM(%esp)
        movsd %xmm1,nb114_fiyM(%esp)
        movsd %xmm2,nb114_fizM(%esp)

        ## M-M interaction 
        movsd nb114_rinvMM(%esp),%xmm0
        movsd %xmm0,%xmm7       ## xmm7=rinv 
        mulsd  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulsd  nb114_qqMM(%esp),%xmm7

        mulsd  %xmm7,%xmm0
        addsd  nb114_vctot(%esp),%xmm7
        movsd %xmm7,nb114_vctot(%esp)

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb114_fjxM(%esp),%xmm3
        movapd nb114_fjyM(%esp),%xmm4
        movapd nb114_fjzM(%esp),%xmm5
        mulsd nb114_dxMM(%esp),%xmm0
        mulsd nb114_dyMM(%esp),%xmm1
        mulsd nb114_dzMM(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb114_fixM(%esp),%xmm0
        addsd nb114_fiyM(%esp),%xmm1
        addsd nb114_fizM(%esp),%xmm2
        movsd %xmm3,nb114_fjxM(%esp)
        movsd %xmm4,nb114_fjyM(%esp)
        movsd %xmm5,nb114_fjzM(%esp)
        movsd %xmm0,nb114_fixM(%esp)
        movsd %xmm1,nb114_fiyM(%esp)
        movsd %xmm2,nb114_fizM(%esp)

        movl nb114_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## Step1 - merge forces
        movlpd nb114_fjxO(%esp),%xmm0
        movlpd nb114_fjzO(%esp),%xmm1
        movlpd nb114_fjyH1(%esp),%xmm2
        movlpd nb114_fjxH2(%esp),%xmm3
        movlpd nb114_fjzH2(%esp),%xmm4
        movlpd nb114_fjyM(%esp),%xmm5

        movhpd nb114_fjyO(%esp),%xmm0
        movhpd nb114_fjxH1(%esp),%xmm1
        movhpd nb114_fjzH1(%esp),%xmm2
        movhpd nb114_fjyH2(%esp),%xmm3
        movhpd nb114_fjxM(%esp),%xmm4
        movhpd nb114_fjzM(%esp),%xmm5

        movlps (%edi,%eax,8),%xmm6
        movhps 8(%edi,%eax,8),%xmm6
        movlps 16(%edi,%eax,8),%xmm7
        movhps 24(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm0
        addpd  %xmm7,%xmm1
        movlps 32(%edi,%eax,8),%xmm6
        movhps 40(%edi,%eax,8),%xmm6
        movlps 48(%edi,%eax,8),%xmm7
        movhps 56(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm2
        addpd  %xmm7,%xmm3
        movlps 64(%edi,%eax,8),%xmm6
        movhps 72(%edi,%eax,8),%xmm6
        movlps 80(%edi,%eax,8),%xmm7
        movhps 88(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm4
        addpd  %xmm7,%xmm5

        movlps %xmm0,(%edi,%eax,8)
        movhps %xmm0,8(%edi,%eax,8)
        movlps %xmm1,16(%edi,%eax,8)
        movhps %xmm1,24(%edi,%eax,8)
        movlps %xmm2,32(%edi,%eax,8)
        movhps %xmm2,40(%edi,%eax,8)
        movlps %xmm3,48(%edi,%eax,8)
        movhps %xmm3,56(%edi,%eax,8)
        movlps %xmm4,64(%edi,%eax,8)
        movhps %xmm4,72(%edi,%eax,8)
        movlps %xmm5,80(%edi,%eax,8)
        movhps %xmm5,88(%edi,%eax,8)

_nb_kernel114_ia32_sse2.nb114_updateouterdata: 
        movl  nb114_ii3(%esp),%ecx
        movl  nb114_faction(%ebp),%edi
        movl  nb114_fshift(%ebp),%esi
        movl  nb114_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb114_fixO(%esp),%xmm0
        movapd nb114_fiyO(%esp),%xmm1
        movapd nb114_fizO(%esp),%xmm2

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
        movsd  (%edi,%ecx,8),%xmm3
        movsd  8(%edi,%ecx,8),%xmm4
        movsd  16(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,(%edi,%ecx,8)
        movsd  %xmm4,8(%edi,%ecx,8)
        movsd  %xmm5,16(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        movapd %xmm0,%xmm6
        movsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm6

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movapd nb114_fixH1(%esp),%xmm0
        movapd nb114_fiyH1(%esp),%xmm1
        movapd nb114_fizH1(%esp),%xmm2

        movhlps %xmm0,%xmm3
        movhlps %xmm1,%xmm4
        movhlps %xmm2,%xmm5
        addsd  %xmm3,%xmm0
        addsd  %xmm4,%xmm1
        addsd  %xmm5,%xmm2 ## sum is in low xmm0-xmm2 

        ## increment i force 
        movsd  24(%edi,%ecx,8),%xmm3
        movsd  32(%edi,%ecx,8),%xmm4
        movsd  40(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,24(%edi,%ecx,8)
        movsd  %xmm4,32(%edi,%ecx,8)
        movsd  %xmm5,40(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movapd nb114_fixH2(%esp),%xmm0
        movapd nb114_fiyH2(%esp),%xmm1
        movapd nb114_fizH2(%esp),%xmm2

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
        movsd  48(%edi,%ecx,8),%xmm3
        movsd  56(%edi,%ecx,8),%xmm4
        movsd  64(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,48(%edi,%ecx,8)
        movsd  %xmm4,56(%edi,%ecx,8)
        movsd  %xmm5,64(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## accumulate Mi forces in xmm0, xmm1, xmm2 
        movapd nb114_fixM(%esp),%xmm0
        movapd nb114_fiyM(%esp),%xmm1
        movapd nb114_fizM(%esp),%xmm2

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
        movsd  72(%edi,%ecx,8),%xmm3
        movsd  80(%edi,%ecx,8),%xmm4
        movsd  88(%edi,%ecx,8),%xmm5
        addsd  %xmm0,%xmm3
        addsd  %xmm1,%xmm4
        addsd  %xmm2,%xmm5
        movsd  %xmm3,72(%edi,%ecx,8)
        movsd  %xmm4,80(%edi,%ecx,8)
        movsd  %xmm5,88(%edi,%ecx,8)

        ## accumulate force in xmm6/xmm7 for fshift 
        addsd %xmm2,%xmm7
        unpcklpd %xmm1,%xmm0
        addpd %xmm0,%xmm6

        ## increment fshift force 
        movlpd (%esi,%edx,8),%xmm3
        movhpd 8(%esi,%edx,8),%xmm3
        movsd  16(%esi,%edx,8),%xmm4
        addpd  %xmm6,%xmm3
        addsd  %xmm7,%xmm4
        movlpd %xmm3,(%esi,%edx,8)
        movhpd %xmm3,8(%esi,%edx,8)
        movsd  %xmm4,16(%esi,%edx,8)

        ## get n from stack
        movl nb114_n(%esp),%esi
        ## get group index for i particle 
        movl  nb114_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb114_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb114_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb114_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb114_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb114_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel114_ia32_sse2.nb114_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb114_n(%esp)
        jmp _nb_kernel114_ia32_sse2.nb114_outer
_nb_kernel114_ia32_sse2.nb114_outerend: 
        ## check if more outer neighborlists remain
        movl  nb114_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel114_ia32_sse2.nb114_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel114_ia32_sse2.nb114_threadloop
_nb_kernel114_ia32_sse2.nb114_end: 
        emms

        movl nb114_nouter(%esp),%eax
        movl nb114_ninner(%esp),%ebx
        movl nb114_outeriter(%ebp),%ecx
        movl nb114_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb114_salign(%esp),%eax
        addl %eax,%esp
        addl $1816,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret





.globl nb_kernel114nf_ia32_sse2
.globl _nb_kernel114nf_ia32_sse2
nb_kernel114nf_ia32_sse2:       
_nb_kernel114nf_ia32_sse2:      
.set nb114nf_p_nri, 8
.set nb114nf_iinr, 12
.set nb114nf_jindex, 16
.set nb114nf_jjnr, 20
.set nb114nf_shift, 24
.set nb114nf_shiftvec, 28
.set nb114nf_fshift, 32
.set nb114nf_gid, 36
.set nb114nf_pos, 40
.set nb114nf_faction, 44
.set nb114nf_charge, 48
.set nb114nf_p_facel, 52
.set nb114nf_argkrf, 56
.set nb114nf_argcrf, 60
.set nb114nf_Vc, 64
.set nb114nf_type, 68
.set nb114nf_p_ntype, 72
.set nb114nf_vdwparam, 76
.set nb114nf_Vvdw, 80
.set nb114nf_p_tabscale, 84
.set nb114nf_VFtab, 88
.set nb114nf_invsqrta, 92
.set nb114nf_dvda, 96
.set nb114nf_p_gbtabscale, 100
.set nb114nf_GBtab, 104
.set nb114nf_p_nthreads, 108
.set nb114nf_count, 112
.set nb114nf_mtx, 116
.set nb114nf_outeriter, 120
.set nb114nf_inneriter, 124
.set nb114nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb114nf_ixO, 0
.set nb114nf_iyO, 16
.set nb114nf_izO, 32
.set nb114nf_ixH1, 48
.set nb114nf_iyH1, 64
.set nb114nf_izH1, 80
.set nb114nf_ixH2, 96
.set nb114nf_iyH2, 112
.set nb114nf_izH2, 128
.set nb114nf_ixM, 144
.set nb114nf_iyM, 160
.set nb114nf_izM, 176
.set nb114nf_jxO, 192
.set nb114nf_jyO, 208
.set nb114nf_jzO, 224
.set nb114nf_jxH1, 240
.set nb114nf_jyH1, 256
.set nb114nf_jzH1, 272
.set nb114nf_jxH2, 288
.set nb114nf_jyH2, 304
.set nb114nf_jzH2, 320
.set nb114nf_jxM, 336
.set nb114nf_jyM, 352
.set nb114nf_jzM, 368
.set nb114nf_qqMM, 384
.set nb114nf_qqMH, 400
.set nb114nf_qqHH, 416
.set nb114nf_two, 432
.set nb114nf_c6, 448
.set nb114nf_c12, 464
.set nb114nf_vctot, 480
.set nb114nf_Vvdwtot, 496
.set nb114nf_half, 512
.set nb114nf_three, 528
.set nb114nf_rsqOO, 544
.set nb114nf_rsqH1H1, 560
.set nb114nf_rsqH1H2, 576
.set nb114nf_rsqH1M, 592
.set nb114nf_rsqH2H1, 608
.set nb114nf_rsqH2H2, 624
.set nb114nf_rsqH2M, 640
.set nb114nf_rsqMH1, 656
.set nb114nf_rsqMH2, 672
.set nb114nf_rsqMM, 688
.set nb114nf_rinvsqOO, 704
.set nb114nf_rinvH1H1, 720
.set nb114nf_rinvH1H2, 736
.set nb114nf_rinvH1M, 752
.set nb114nf_rinvH2H1, 768
.set nb114nf_rinvH2H2, 784
.set nb114nf_rinvH2M, 800
.set nb114nf_rinvMH1, 816
.set nb114nf_rinvMH2, 832
.set nb114nf_rinvMM, 848
.set nb114nf_is3, 864
.set nb114nf_ii3, 868
.set nb114nf_innerjjnr, 872
.set nb114nf_innerk, 876
.set nb114nf_n, 880
.set nb114nf_nn1, 884
.set nb114nf_nri, 888
.set nb114nf_nouter, 892
.set nb114nf_ninner, 896
.set nb114nf_salign, 900
        pushl %ebp
        movl %esp,%ebp
        pushl %eax
        pushl %ebx
        pushl %ecx
        pushl %edx
        pushl %esi
        pushl %edi
        subl $904,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb114nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb114nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb114nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb114nf_nouter(%esp)
        movl %eax,nb114nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb114nf_half(%esp)
        movl %ebx,nb114nf_half+4(%esp)
        movsd nb114nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb114nf_half(%esp)
        movapd %xmm2,nb114nf_two(%esp)
        movapd %xmm3,nb114nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb114nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb114nf_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb114nf_p_facel(%ebp),%esi
        movsd (%esi),%xmm6
        mulsd  %xmm3,%xmm3
        mulsd  %xmm5,%xmm4
        mulsd  %xmm5,%xmm5
        mulsd  %xmm6,%xmm3
        mulsd  %xmm6,%xmm4
        mulsd  %xmm6,%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb114nf_qqMM(%esp)
        movapd %xmm4,nb114nf_qqMH(%esp)
        movapd %xmm5,nb114nf_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb114nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb114nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb114nf_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb114nf_c6(%esp)
        movapd %xmm1,nb114nf_c12(%esp)

_nb_kernel114nf_ia32_sse2.nb114nf_threadloop: 
        movl  nb114nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel114nf_ia32_sse2.nb114nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel114nf_ia32_sse2.nb114nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb114nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb114nf_n(%esp)
        movl %ebx,nb114nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel114nf_ia32_sse2.nb114nf_outerstart
        jmp _nb_kernel114nf_ia32_sse2.nb114nf_end

_nb_kernel114nf_ia32_sse2.nb114nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb114nf_nouter(%esp),%ebx
        movl %ebx,nb114nf_nouter(%esp)

_nb_kernel114nf_ia32_sse2.nb114nf_outer: 
        movl  nb114nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb114nf_is3(%esp)            ## store is3 

        movl  nb114nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb114nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb114nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb114nf_ii3(%esp)

        addsd (%eax,%ebx,8),%xmm3       ## ox
        addsd 8(%eax,%ebx,8),%xmm4      ## oy
        addsd 16(%eax,%ebx,8),%xmm5     ## oz   
        addsd 24(%eax,%ebx,8),%xmm6     ## h1x
        addsd 32(%eax,%ebx,8),%xmm7     ## h1y
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm7,%xmm7
        movapd %xmm3,nb114nf_ixO(%esp)
        movapd %xmm4,nb114nf_iyO(%esp)
        movapd %xmm5,nb114nf_izO(%esp)
        movapd %xmm6,nb114nf_ixH1(%esp)
        movapd %xmm7,nb114nf_iyH1(%esp)

        movsd %xmm2,%xmm6
        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 40(%eax,%ebx,8),%xmm6    ## h1z
        addsd 48(%eax,%ebx,8),%xmm0    ## h2x
        addsd 56(%eax,%ebx,8),%xmm1    ## h2y
        addsd 64(%eax,%ebx,8),%xmm2    ## h2z
        addsd 72(%eax,%ebx,8),%xmm3    ## mx
        addsd 80(%eax,%ebx,8),%xmm4    ## my
        addsd 88(%eax,%ebx,8),%xmm5    ## mz

        shufpd $0,%xmm6,%xmm6
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm6,nb114nf_izH1(%esp)
        movapd %xmm0,nb114nf_ixH2(%esp)
        movapd %xmm1,nb114nf_iyH2(%esp)
        movapd %xmm2,nb114nf_izH2(%esp)
        movapd %xmm3,nb114nf_ixM(%esp)
        movapd %xmm4,nb114nf_iyM(%esp)
        movapd %xmm5,nb114nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb114nf_vctot(%esp)
        movapd %xmm4,nb114nf_Vvdwtot(%esp)

        movl  nb114nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb114nf_pos(%ebp),%esi
        movl  nb114nf_faction(%ebp),%edi
        movl  nb114nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb114nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb114nf_ninner(%esp),%ecx
        movl  %ecx,nb114nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb114nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel114nf_ia32_sse2.nb114nf_unroll_loop
        jmp   _nb_kernel114nf_ia32_sse2.nb114nf_checksingle
_nb_kernel114nf_ia32_sse2.nb114nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb114nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb114nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb114nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movlpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%eax,8),%xmm0
        movhpd 8(%esi,%ebx,8),%xmm2
        movlpd 16(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%ebx,8),%xmm5
        movhpd 24(%esi,%eax,8),%xmm3
        movhpd 24(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## ox 
        unpckhpd %xmm2,%xmm1 ## oy
        unpcklpd %xmm5,%xmm3 ## ox 
        unpckhpd %xmm5,%xmm4 ## oy
        movapd  %xmm0,nb114nf_jxO(%esp)
        movapd  %xmm1,nb114nf_jyO(%esp)
        movapd  %xmm3,nb114nf_jzO(%esp)
        movapd  %xmm4,nb114nf_jxH1(%esp)

        ## load h1y, h1z, h2x, h2y 
        movlpd 32(%esi,%eax,8),%xmm0
        movlpd 32(%esi,%ebx,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm0
        movhpd 40(%esi,%ebx,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movlpd 48(%esi,%ebx,8),%xmm5
        movhpd 56(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h1y
        unpckhpd %xmm2,%xmm1 ## h1z
        unpcklpd %xmm5,%xmm3 ## h2x
        unpckhpd %xmm5,%xmm4 ## h2y
        movapd  %xmm0,nb114nf_jyH1(%esp)
        movapd  %xmm1,nb114nf_jzH1(%esp)
        movapd  %xmm3,nb114nf_jxH2(%esp)
        movapd  %xmm4,nb114nf_jyH2(%esp)

        ## load h2z, mx, my, mz
        movlpd 64(%esi,%eax,8),%xmm0
        movlpd 64(%esi,%ebx,8),%xmm2
        movhpd 72(%esi,%eax,8),%xmm0
        movhpd 72(%esi,%ebx,8),%xmm2
        movlpd 80(%esi,%eax,8),%xmm3
        movlpd 80(%esi,%ebx,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm3
        movhpd 88(%esi,%ebx,8),%xmm5
        movapd %xmm0,%xmm1
        movapd %xmm3,%xmm4
        unpcklpd %xmm2,%xmm0 ## h2z
        unpckhpd %xmm2,%xmm1 ## mx
        unpcklpd %xmm5,%xmm3 ## my
        unpckhpd %xmm5,%xmm4 ## mz
        movapd  %xmm0,nb114nf_jzH2(%esp)
        movapd  %xmm1,nb114nf_jxM(%esp)
        movapd  %xmm3,nb114nf_jyM(%esp)
        movapd  %xmm4,nb114nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb114nf_ixO(%esp),%xmm0
        movapd nb114nf_iyO(%esp),%xmm1
        movapd nb114nf_izO(%esp),%xmm2
        movapd nb114nf_ixH1(%esp),%xmm3
        movapd nb114nf_iyH1(%esp),%xmm4
        movapd nb114nf_izH1(%esp),%xmm5
        subpd  nb114nf_jxO(%esp),%xmm0
        subpd  nb114nf_jyO(%esp),%xmm1
        subpd  nb114nf_jzO(%esp),%xmm2
        subpd  nb114nf_jxH1(%esp),%xmm3
        subpd  nb114nf_jyH1(%esp),%xmm4
        subpd  nb114nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqOO(%esp)
        movapd %xmm3,nb114nf_rsqH1H1(%esp)

        movapd nb114nf_ixH1(%esp),%xmm0
        movapd nb114nf_iyH1(%esp),%xmm1
        movapd nb114nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb114nf_jxH2(%esp),%xmm0
        subpd  nb114nf_jyH2(%esp),%xmm1
        subpd  nb114nf_jzH2(%esp),%xmm2
        subpd  nb114nf_jxM(%esp),%xmm3
        subpd  nb114nf_jyM(%esp),%xmm4
        subpd  nb114nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqH1H2(%esp)
        movapd %xmm3,nb114nf_rsqH1M(%esp)

        movapd nb114nf_ixH2(%esp),%xmm0
        movapd nb114nf_iyH2(%esp),%xmm1
        movapd nb114nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb114nf_jxH1(%esp),%xmm0
        subpd  nb114nf_jyH1(%esp),%xmm1
        subpd  nb114nf_jzH1(%esp),%xmm2
        subpd  nb114nf_jxH2(%esp),%xmm3
        subpd  nb114nf_jyH2(%esp),%xmm4
        subpd  nb114nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqH2H1(%esp)
        movapd %xmm3,nb114nf_rsqH2H2(%esp)

        movapd nb114nf_ixH2(%esp),%xmm0
        movapd nb114nf_iyH2(%esp),%xmm1
        movapd nb114nf_izH2(%esp),%xmm2
        movapd nb114nf_ixM(%esp),%xmm3
        movapd nb114nf_iyM(%esp),%xmm4
        movapd nb114nf_izM(%esp),%xmm5
        subpd  nb114nf_jxM(%esp),%xmm0
        subpd  nb114nf_jyM(%esp),%xmm1
        subpd  nb114nf_jzM(%esp),%xmm2
        subpd  nb114nf_jxH1(%esp),%xmm3
        subpd  nb114nf_jyH1(%esp),%xmm4
        subpd  nb114nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqH2M(%esp)
        movapd %xmm4,nb114nf_rsqMH1(%esp)

        movapd nb114nf_ixM(%esp),%xmm0
        movapd nb114nf_iyM(%esp),%xmm1
        movapd nb114nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb114nf_jxH2(%esp),%xmm0
        subpd  nb114nf_jyH2(%esp),%xmm1
        subpd  nb114nf_jzH2(%esp),%xmm2
        subpd  nb114nf_jxM(%esp),%xmm3
        subpd  nb114nf_jyM(%esp),%xmm4
        subpd  nb114nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqMH2(%esp)
        movapd %xmm4,nb114nf_rsqMM(%esp)

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
        movapd  nb114nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb114nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb114nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114nf_rinvMH2(%esp)
        movapd %xmm5,nb114nf_rinvMM(%esp)

        movapd nb114nf_rsqOO(%esp),%xmm0
        movapd nb114nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb114nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb114nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb114nf_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb114nf_rinvsqOO(%esp)
        movapd %xmm5,nb114nf_rinvH1H1(%esp)

        movapd nb114nf_rsqH1H2(%esp),%xmm0
        movapd nb114nf_rsqH1M(%esp),%xmm4
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
        movapd  nb114nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb114nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb114nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114nf_rinvH1H2(%esp)
        movapd %xmm5,nb114nf_rinvH1M(%esp)

        movapd nb114nf_rsqH2H1(%esp),%xmm0
        movapd nb114nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb114nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb114nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb114nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114nf_rinvH2H1(%esp)
        movapd %xmm5,nb114nf_rinvH2H2(%esp)

        movapd nb114nf_rsqMH1(%esp),%xmm0
        movapd nb114nf_rsqH2M(%esp),%xmm4
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
        movapd  nb114nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb114nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb114nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114nf_rinvMH1(%esp)
        movapd %xmm5,nb114nf_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb114nf_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb114nf_c6(%esp),%xmm1
        mulpd  nb114nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb114nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb114nf_Vvdwtot(%esp)

        ## Coulomb interactions.
        ## Sum rinv for H-H interactions in xmm0, H-M in xmm1.
        ## Use xmm2 for the M-M pair
        movapd nb114nf_rinvH1H1(%esp),%xmm0
        movapd nb114nf_rinvH1M(%esp),%xmm1
        movapd nb114nf_rinvMM(%esp),%xmm2

        addpd  nb114nf_rinvH1H2(%esp),%xmm0
        addpd  nb114nf_rinvH2M(%esp),%xmm1
        addpd  nb114nf_rinvH2H1(%esp),%xmm0
        addpd  nb114nf_rinvMH1(%esp),%xmm1
        addpd  nb114nf_rinvH2H2(%esp),%xmm0
        addpd  nb114nf_rinvMH2(%esp),%xmm1
        mulpd  nb114nf_qqHH(%esp),%xmm0
        mulpd  nb114nf_qqMH(%esp),%xmm1
        mulpd  nb114nf_qqMM(%esp),%xmm2
        addpd  %xmm1,%xmm0
        addpd  nb114nf_vctot(%esp),%xmm2
        addpd  %xmm0,%xmm2
        movapd %xmm2,nb114nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb114nf_innerk(%esp)
        jl    _nb_kernel114nf_ia32_sse2.nb114nf_checksingle
        jmp   _nb_kernel114nf_ia32_sse2.nb114nf_unroll_loop
_nb_kernel114nf_ia32_sse2.nb114nf_checksingle: 
        movl  nb114nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel114nf_ia32_sse2.nb114nf_dosingle
        jmp   _nb_kernel114nf_ia32_sse2.nb114nf_updateouterdata
_nb_kernel114nf_ia32_sse2.nb114nf_dosingle: 
        movl  nb114nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb114nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 

        ## move j coordinates to local temp variables 
        ## load ox, oy, oz, h1x
        movlpd (%esi,%eax,8),%xmm0
        movhpd 8(%esi,%eax,8),%xmm0
        movlpd 16(%esi,%eax,8),%xmm1
        movhpd 24(%esi,%eax,8),%xmm1
        movlpd 32(%esi,%eax,8),%xmm2
        movhpd 40(%esi,%eax,8),%xmm2
        movlpd 48(%esi,%eax,8),%xmm3
        movhpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 72(%esi,%eax,8),%xmm4
        movlpd 80(%esi,%eax,8),%xmm5
        movhpd 88(%esi,%eax,8),%xmm5
        movsd  %xmm0,nb114nf_jxO(%esp)
        movsd  %xmm1,nb114nf_jzO(%esp)
        movsd  %xmm2,nb114nf_jyH1(%esp)
        movsd  %xmm3,nb114nf_jxH2(%esp)
        movsd  %xmm4,nb114nf_jzH2(%esp)
        movsd  %xmm5,nb114nf_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb114nf_jyO(%esp)
        movsd  %xmm1,nb114nf_jxH1(%esp)
        movsd  %xmm2,nb114nf_jzH1(%esp)
        movsd  %xmm3,nb114nf_jyH2(%esp)
        movsd  %xmm4,nb114nf_jxM(%esp)
        movsd  %xmm5,nb114nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb114nf_ixO(%esp),%xmm0
        movapd nb114nf_iyO(%esp),%xmm1
        movapd nb114nf_izO(%esp),%xmm2
        movapd nb114nf_ixH1(%esp),%xmm3
        movapd nb114nf_iyH1(%esp),%xmm4
        movapd nb114nf_izH1(%esp),%xmm5
        subsd  nb114nf_jxO(%esp),%xmm0
        subsd  nb114nf_jyO(%esp),%xmm1
        subsd  nb114nf_jzO(%esp),%xmm2
        subsd  nb114nf_jxH1(%esp),%xmm3
        subsd  nb114nf_jyH1(%esp),%xmm4
        subsd  nb114nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqOO(%esp)
        movapd %xmm3,nb114nf_rsqH1H1(%esp)

        movapd nb114nf_ixH1(%esp),%xmm0
        movapd nb114nf_iyH1(%esp),%xmm1
        movapd nb114nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb114nf_jxH2(%esp),%xmm0
        subsd  nb114nf_jyH2(%esp),%xmm1
        subsd  nb114nf_jzH2(%esp),%xmm2
        subsd  nb114nf_jxM(%esp),%xmm3
        subsd  nb114nf_jyM(%esp),%xmm4
        subsd  nb114nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqH1H2(%esp)
        movapd %xmm3,nb114nf_rsqH1M(%esp)

        movapd nb114nf_ixH2(%esp),%xmm0
        movapd nb114nf_iyH2(%esp),%xmm1
        movapd nb114nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb114nf_jxH1(%esp),%xmm0
        subsd  nb114nf_jyH1(%esp),%xmm1
        subsd  nb114nf_jzH1(%esp),%xmm2
        subsd  nb114nf_jxH2(%esp),%xmm3
        subsd  nb114nf_jyH2(%esp),%xmm4
        subsd  nb114nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqH2H1(%esp)
        movapd %xmm3,nb114nf_rsqH2H2(%esp)

        movapd nb114nf_ixH2(%esp),%xmm0
        movapd nb114nf_iyH2(%esp),%xmm1
        movapd nb114nf_izH2(%esp),%xmm2
        movapd nb114nf_ixM(%esp),%xmm3
        movapd nb114nf_iyM(%esp),%xmm4
        movapd nb114nf_izM(%esp),%xmm5
        subsd  nb114nf_jxM(%esp),%xmm0
        subsd  nb114nf_jyM(%esp),%xmm1
        subsd  nb114nf_jzM(%esp),%xmm2
        subsd  nb114nf_jxH1(%esp),%xmm3
        subsd  nb114nf_jyH1(%esp),%xmm4
        subsd  nb114nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqH2M(%esp)
        movapd %xmm4,nb114nf_rsqMH1(%esp)

        movapd nb114nf_ixM(%esp),%xmm0
        movapd nb114nf_iyM(%esp),%xmm1
        movapd nb114nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb114nf_jxH2(%esp),%xmm0
        subsd  nb114nf_jyH2(%esp),%xmm1
        subsd  nb114nf_jzH2(%esp),%xmm2
        subsd  nb114nf_jxM(%esp),%xmm3
        subsd  nb114nf_jyM(%esp),%xmm4
        subsd  nb114nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb114nf_rsqMH2(%esp)
        movapd %xmm4,nb114nf_rsqMM(%esp)

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
        movapd  nb114nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb114nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb114nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114nf_rinvMH2(%esp)
        movapd %xmm5,nb114nf_rinvMM(%esp)

        movapd nb114nf_rsqOO(%esp),%xmm0
        movapd nb114nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb114nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb114nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb114nf_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb114nf_rinvsqOO(%esp)
        movapd %xmm5,nb114nf_rinvH1H1(%esp)

        movapd nb114nf_rsqH1H2(%esp),%xmm0
        movapd nb114nf_rsqH1M(%esp),%xmm4
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
        movapd  nb114nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb114nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb114nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114nf_rinvH1H2(%esp)
        movapd %xmm5,nb114nf_rinvH1M(%esp)

        movapd nb114nf_rsqH2H1(%esp),%xmm0
        movapd nb114nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb114nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb114nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb114nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114nf_rinvH2H1(%esp)
        movapd %xmm5,nb114nf_rinvH2H2(%esp)

        movapd nb114nf_rsqMH1(%esp),%xmm0
        movapd nb114nf_rsqH2M(%esp),%xmm4
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
        movapd  nb114nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb114nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb114nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb114nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb114nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb114nf_rinvMH1(%esp)
        movapd %xmm5,nb114nf_rinvH2M(%esp)

        ## start with OO interaction 
        movsd nb114nf_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movsd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb114nf_c6(%esp),%xmm1
        mulsd  nb114nf_c12(%esp),%xmm2
        movsd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb114nf_Vvdwtot(%esp),%xmm3
        movsd %xmm3,nb114nf_Vvdwtot(%esp)

        ## Coulomb interactions.
        ## Sum rinv for H-H interactions in xmm0, H-M in xmm1.
        ## Use xmm2 for the M-M pair
        movsd nb114nf_rinvH1H1(%esp),%xmm0
        movsd nb114nf_rinvH1M(%esp),%xmm1
        movsd nb114nf_rinvMM(%esp),%xmm2

        addsd  nb114nf_rinvH1H2(%esp),%xmm0
        addsd  nb114nf_rinvH2M(%esp),%xmm1
        addsd  nb114nf_rinvH2H1(%esp),%xmm0
        addsd  nb114nf_rinvMH1(%esp),%xmm1
        addsd  nb114nf_rinvH2H2(%esp),%xmm0
        addsd  nb114nf_rinvMH2(%esp),%xmm1
        mulsd  nb114nf_qqHH(%esp),%xmm0
        mulsd  nb114nf_qqMH(%esp),%xmm1
        mulsd  nb114nf_qqMM(%esp),%xmm2
        addsd  %xmm1,%xmm0
        addsd  nb114nf_vctot(%esp),%xmm2
        addsd  %xmm0,%xmm2
        movsd %xmm2,nb114nf_vctot(%esp)


_nb_kernel114nf_ia32_sse2.nb114nf_updateouterdata: 
        ## get n from stack
        movl nb114nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb114nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb114nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb114nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb114nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb114nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb114nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel114nf_ia32_sse2.nb114nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb114nf_n(%esp)
        jmp _nb_kernel114nf_ia32_sse2.nb114nf_outer
_nb_kernel114nf_ia32_sse2.nb114nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb114nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel114nf_ia32_sse2.nb114nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel114nf_ia32_sse2.nb114nf_threadloop
_nb_kernel114nf_ia32_sse2.nb114nf_end: 
        emms

        movl nb114nf_nouter(%esp),%eax
        movl nb114nf_ninner(%esp),%ebx
        movl nb114nf_outeriter(%ebp),%ecx
        movl nb114nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb114nf_salign(%esp),%eax
        addl %eax,%esp
        addl $904,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



