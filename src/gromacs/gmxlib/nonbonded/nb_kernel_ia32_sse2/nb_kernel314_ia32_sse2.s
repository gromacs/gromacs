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



.globl nb_kernel314_ia32_sse2
.globl _nb_kernel314_ia32_sse2
nb_kernel314_ia32_sse2: 
_nb_kernel314_ia32_sse2:        
.set nb314_p_nri, 8
.set nb314_iinr, 12
.set nb314_jindex, 16
.set nb314_jjnr, 20
.set nb314_shift, 24
.set nb314_shiftvec, 28
.set nb314_fshift, 32
.set nb314_gid, 36
.set nb314_pos, 40
.set nb314_faction, 44
.set nb314_charge, 48
.set nb314_p_facel, 52
.set nb314_argkrf, 56
.set nb314_argcrf, 60
.set nb314_Vc, 64
.set nb314_type, 68
.set nb314_p_ntype, 72
.set nb314_vdwparam, 76
.set nb314_Vvdw, 80
.set nb314_p_tabscale, 84
.set nb314_VFtab, 88
.set nb314_invsqrta, 92
.set nb314_dvda, 96
.set nb314_p_gbtabscale, 100
.set nb314_GBtab, 104
.set nb314_p_nthreads, 108
.set nb314_count, 112
.set nb314_mtx, 116
.set nb314_outeriter, 120
.set nb314_inneriter, 124
.set nb314_work, 128
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
.set nb314_fjxM, 1344
.set nb314_fjyM, 1360
.set nb314_fjzM, 1376
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
.set nb314_innerjjnr, 1784
.set nb314_innerk, 1788
.set nb314_n, 1792
.set nb314_nn1, 1796
.set nb314_nri, 1800
.set nb314_nouter, 1804
.set nb314_ninner, 1808
.set nb314_salign, 1812
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
        movl %eax,nb314_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb314_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb314_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb314_nouter(%esp)
        movl %eax,nb314_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb314_half(%esp)
        movl %ebx,nb314_half+4(%esp)
        movsd nb314_half(%esp),%xmm1
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
        movapd %xmm1,nb314_half(%esp)
        movapd %xmm2,nb314_two(%esp)
        movapd %xmm3,nb314_three(%esp)
        movapd %xmm4,nb314_six(%esp)
        movapd %xmm5,nb314_twelve(%esp)
        movl nb314_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3

        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb314_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb314_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb314_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb314_p_facel(%ebp),%esi
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
        movapd %xmm3,nb314_qqMM(%esp)
        movapd %xmm4,nb314_qqMH(%esp)
        movapd %xmm5,nb314_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb314_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb314_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb314_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movlpd 8(%eax,%edx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb314_c6(%esp)
        movapd %xmm1,nb314_c12(%esp)

_nb_kernel314_ia32_sse2.nb314_threadloop: 
        movl  nb314_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel314_ia32_sse2.nb314_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel314_ia32_sse2.nb314_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb314_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb314_n(%esp)
        movl %ebx,nb314_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel314_ia32_sse2.nb314_outerstart
        jmp _nb_kernel314_ia32_sse2.nb314_end

_nb_kernel314_ia32_sse2.nb314_outerstart: 
        ## ebx contains number of outer iterations
        addl nb314_nouter(%esp),%ebx
        movl %ebx,nb314_nouter(%esp)

_nb_kernel314_ia32_sse2.nb314_outer: 
        movl  nb314_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb314_is3(%esp)      ## store is3 

        movl  nb314_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb314_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb314_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb314_ii3(%esp)

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
        movapd %xmm3,nb314_ixO(%esp)
        movapd %xmm4,nb314_iyO(%esp)
        movapd %xmm5,nb314_izO(%esp)
        movapd %xmm6,nb314_ixH1(%esp)
        movapd %xmm7,nb314_iyH1(%esp)

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
        movapd %xmm6,nb314_izH1(%esp)
        movapd %xmm0,nb314_ixH2(%esp)
        movapd %xmm1,nb314_iyH2(%esp)
        movapd %xmm2,nb314_izH2(%esp)
        movapd %xmm3,nb314_ixM(%esp)
        movapd %xmm4,nb314_iyM(%esp)
        movapd %xmm5,nb314_izM(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb314_vctot(%esp)
        movapd %xmm4,nb314_Vvdwtot(%esp)
        movapd %xmm4,nb314_fixO(%esp)
        movapd %xmm4,nb314_fiyO(%esp)
        movapd %xmm4,nb314_fizO(%esp)
        movapd %xmm4,nb314_fixH1(%esp)
        movapd %xmm4,nb314_fiyH1(%esp)
        movapd %xmm4,nb314_fizH1(%esp)
        movapd %xmm4,nb314_fixH2(%esp)
        movapd %xmm4,nb314_fiyH2(%esp)
        movapd %xmm4,nb314_fizH2(%esp)
        movapd %xmm4,nb314_fixM(%esp)
        movapd %xmm4,nb314_fiyM(%esp)
        movapd %xmm4,nb314_fizM(%esp)

        movl  nb314_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb314_pos(%ebp),%esi
        movl  nb314_faction(%ebp),%edi
        movl  nb314_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb314_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb314_ninner(%esp),%ecx
        movl  %ecx,nb314_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb314_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel314_ia32_sse2.nb314_unroll_loop
        jmp   _nb_kernel314_ia32_sse2.nb314_checksingle
_nb_kernel314_ia32_sse2.nb314_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb314_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb314_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb314_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm0,nb314_jxO(%esp)
        movapd  %xmm1,nb314_jyO(%esp)
        movapd  %xmm3,nb314_jzO(%esp)
        movapd  %xmm4,nb314_jxH1(%esp)

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
        movapd  %xmm0,nb314_jyH1(%esp)
        movapd  %xmm1,nb314_jzH1(%esp)
        movapd  %xmm3,nb314_jxH2(%esp)
        movapd  %xmm4,nb314_jyH2(%esp)

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
        movapd  %xmm0,nb314_jzH2(%esp)
        movapd  %xmm1,nb314_jxM(%esp)
        movapd  %xmm3,nb314_jyM(%esp)
        movapd  %xmm4,nb314_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb314_ixO(%esp),%xmm0
        movapd nb314_iyO(%esp),%xmm1
        movapd nb314_izO(%esp),%xmm2
        movapd nb314_ixH1(%esp),%xmm3
        movapd nb314_iyH1(%esp),%xmm4
        movapd nb314_izH1(%esp),%xmm5
        subpd  nb314_jxO(%esp),%xmm0
        subpd  nb314_jyO(%esp),%xmm1
        subpd  nb314_jzO(%esp),%xmm2
        subpd  nb314_jxH1(%esp),%xmm3
        subpd  nb314_jyH1(%esp),%xmm4
        subpd  nb314_jzH1(%esp),%xmm5
        movapd %xmm0,nb314_dxOO(%esp)
        movapd %xmm1,nb314_dyOO(%esp)
        movapd %xmm2,nb314_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH1H1(%esp)
        movapd %xmm4,nb314_dyH1H1(%esp)
        movapd %xmm5,nb314_dzH1H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb314_rsqOO(%esp)
        movapd %xmm3,nb314_rsqH1H1(%esp)

        movapd nb314_ixH1(%esp),%xmm0
        movapd nb314_iyH1(%esp),%xmm1
        movapd nb314_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314_jxH2(%esp),%xmm0
        subpd  nb314_jyH2(%esp),%xmm1
        subpd  nb314_jzH2(%esp),%xmm2
        subpd  nb314_jxM(%esp),%xmm3
        subpd  nb314_jyM(%esp),%xmm4
        subpd  nb314_jzM(%esp),%xmm5
        movapd %xmm0,nb314_dxH1H2(%esp)
        movapd %xmm1,nb314_dyH1H2(%esp)
        movapd %xmm2,nb314_dzH1H2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH1M(%esp)
        movapd %xmm4,nb314_dyH1M(%esp)
        movapd %xmm5,nb314_dzH1M(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb314_rsqH1H2(%esp)
        movapd %xmm3,nb314_rsqH1M(%esp)

        movapd nb314_ixH2(%esp),%xmm0
        movapd nb314_iyH2(%esp),%xmm1
        movapd nb314_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314_jxH1(%esp),%xmm0
        subpd  nb314_jyH1(%esp),%xmm1
        subpd  nb314_jzH1(%esp),%xmm2
        subpd  nb314_jxH2(%esp),%xmm3
        subpd  nb314_jyH2(%esp),%xmm4
        subpd  nb314_jzH2(%esp),%xmm5
        movapd %xmm0,nb314_dxH2H1(%esp)
        movapd %xmm1,nb314_dyH2H1(%esp)
        movapd %xmm2,nb314_dzH2H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH2H2(%esp)
        movapd %xmm4,nb314_dyH2H2(%esp)
        movapd %xmm5,nb314_dzH2H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb314_rsqH2H1(%esp)
        movapd %xmm3,nb314_rsqH2H2(%esp)

        movapd nb314_ixH2(%esp),%xmm0
        movapd nb314_iyH2(%esp),%xmm1
        movapd nb314_izH2(%esp),%xmm2
        movapd nb314_ixM(%esp),%xmm3
        movapd nb314_iyM(%esp),%xmm4
        movapd nb314_izM(%esp),%xmm5
        subpd  nb314_jxM(%esp),%xmm0
        subpd  nb314_jyM(%esp),%xmm1
        subpd  nb314_jzM(%esp),%xmm2
        subpd  nb314_jxH1(%esp),%xmm3
        subpd  nb314_jyH1(%esp),%xmm4
        subpd  nb314_jzH1(%esp),%xmm5
        movapd %xmm0,nb314_dxH2M(%esp)
        movapd %xmm1,nb314_dyH2M(%esp)
        movapd %xmm2,nb314_dzH2M(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxMH1(%esp)
        movapd %xmm4,nb314_dyMH1(%esp)
        movapd %xmm5,nb314_dzMH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb314_rsqH2M(%esp)
        movapd %xmm4,nb314_rsqMH1(%esp)

        movapd nb314_ixM(%esp),%xmm0
        movapd nb314_iyM(%esp),%xmm1
        movapd nb314_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314_jxH2(%esp),%xmm0
        subpd  nb314_jyH2(%esp),%xmm1
        subpd  nb314_jzH2(%esp),%xmm2
        subpd  nb314_jxM(%esp),%xmm3
        subpd  nb314_jyM(%esp),%xmm4
        subpd  nb314_jzM(%esp),%xmm5
        movapd %xmm0,nb314_dxMH2(%esp)
        movapd %xmm1,nb314_dyMH2(%esp)
        movapd %xmm2,nb314_dzMH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxMM(%esp)
        movapd %xmm4,nb314_dyMM(%esp)
        movapd %xmm5,nb314_dzMM(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb314_rsqMH2(%esp)
        movapd %xmm4,nb314_rsqMM(%esp)

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
        movapd  nb314_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm3   ## iter1 
        mulpd   nb314_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm1   ## rinv 
        mulpd   nb314_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314_rinvMH2(%esp)
        movapd %xmm5,nb314_rinvMM(%esp)

        movapd nb314_rsqOO(%esp),%xmm0
        movapd nb314_rsqH1H1(%esp),%xmm4
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
        movapd  nb314_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb314_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm1   ## rinv 
        mulpd   nb314_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb314_rinvsqOO(%esp)
        movapd %xmm5,nb314_rinvH1H1(%esp)

        movapd nb314_rsqH1H2(%esp),%xmm0
        movapd nb314_rsqH1M(%esp),%xmm4
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
        movapd  nb314_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm3   ## iter1 
        mulpd   nb314_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm1   ## rinv 
        mulpd   nb314_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314_rinvH1H2(%esp)
        movapd %xmm5,nb314_rinvH1M(%esp)

        movapd nb314_rsqH2H1(%esp),%xmm0
        movapd nb314_rsqH2H2(%esp),%xmm4
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
        movapd  nb314_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm3   ## iter1a 
        mulpd   nb314_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm1   ## rinv 
        mulpd   nb314_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314_rinvH2H1(%esp)
        movapd %xmm5,nb314_rinvH2H2(%esp)

        movapd nb314_rsqMH1(%esp),%xmm0
        movapd nb314_rsqH2M(%esp),%xmm4
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
        movapd  nb314_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm3   ## iter1a 
        mulpd   nb314_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314_half(%esp),%xmm1   ## rinv 
        mulpd   nb314_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314_rinvMH1(%esp)
        movapd %xmm5,nb314_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb314_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb314_c6(%esp),%xmm1
        mulpd  nb314_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb314_Vvdwtot(%esp),%xmm3
        mulpd  nb314_six(%esp),%xmm1
        mulpd  nb314_twelve(%esp),%xmm2
        subpd  %xmm1,%xmm2
        mulpd  %xmm0,%xmm2
        movapd %xmm3,nb314_Vvdwtot(%esp)

        movapd %xmm2,%xmm0
        movapd %xmm2,%xmm1 ## fscal

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb314_dxOO(%esp),%xmm0
        mulpd nb314_dyOO(%esp),%xmm1
        mulpd nb314_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixO(%esp),%xmm0
        addpd nb314_fiyO(%esp),%xmm1
        addpd nb314_fizO(%esp),%xmm2
        movapd %xmm3,nb314_fjxO(%esp)
        movapd %xmm4,nb314_fjyO(%esp)
        movapd %xmm5,nb314_fjzO(%esp)
        movapd %xmm0,nb314_fixO(%esp)
        movapd %xmm1,nb314_fiyO(%esp)
        movapd %xmm2,nb314_fizO(%esp)

        ## H1-H1 interaction 
        movapd nb314_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1

        movd %eax,%mm0
        movd %ebx,%mm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb314_dxH1H1(%esp),%xmm0
        mulpd nb314_dyH1H1(%esp),%xmm1
        mulpd nb314_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixH1(%esp),%xmm0
        addpd nb314_fiyH1(%esp),%xmm1
        addpd nb314_fizH1(%esp),%xmm2
        movapd %xmm3,nb314_fjxH1(%esp)
        movapd %xmm4,nb314_fjyH1(%esp)
        movapd %xmm5,nb314_fjzH1(%esp)
        movapd %xmm0,nb314_fixH1(%esp)
        movapd %xmm1,nb314_fiyH1(%esp)
        movapd %xmm2,nb314_fizH1(%esp)

        ## H1-H2 interaction  
        movapd nb314_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 


        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb314_dxH1H2(%esp),%xmm0
        mulpd nb314_dyH1H2(%esp),%xmm1
        mulpd nb314_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixH1(%esp),%xmm0
        addpd nb314_fiyH1(%esp),%xmm1
        addpd nb314_fizH1(%esp),%xmm2
        movapd %xmm3,nb314_fjxH2(%esp)
        movapd %xmm4,nb314_fjyH2(%esp)
        movapd %xmm5,nb314_fjzH2(%esp)
        movapd %xmm0,nb314_fixH1(%esp)
        movapd %xmm1,nb314_fiyH1(%esp)
        movapd %xmm2,nb314_fizH1(%esp)

        ## H1-M interaction 
        movapd nb314_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb314_dxH1M(%esp),%xmm0
        mulpd nb314_dyH1M(%esp),%xmm1
        mulpd nb314_dzH1M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixH1(%esp),%xmm0
        addpd nb314_fiyH1(%esp),%xmm1
        addpd nb314_fizH1(%esp),%xmm2
        movapd %xmm3,nb314_fjxM(%esp)
        movapd %xmm4,nb314_fjyM(%esp)
        movapd %xmm5,nb314_fjzM(%esp)
        movapd %xmm0,nb314_fixH1(%esp)
        movapd %xmm1,nb314_fiyH1(%esp)
        movapd %xmm2,nb314_fizH1(%esp)

        ## H2-H1 interaction 
        movapd nb314_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxH1(%esp),%xmm3
        movapd nb314_fjyH1(%esp),%xmm4
        movapd nb314_fjzH1(%esp),%xmm5
        mulpd nb314_dxH2H1(%esp),%xmm0
        mulpd nb314_dyH2H1(%esp),%xmm1
        mulpd nb314_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixH2(%esp),%xmm0
        addpd nb314_fiyH2(%esp),%xmm1
        addpd nb314_fizH2(%esp),%xmm2
        movapd %xmm3,nb314_fjxH1(%esp)
        movapd %xmm4,nb314_fjyH1(%esp)
        movapd %xmm5,nb314_fjzH1(%esp)
        movapd %xmm0,nb314_fixH2(%esp)
        movapd %xmm1,nb314_fiyH2(%esp)
        movapd %xmm2,nb314_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb314_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxH2(%esp),%xmm3
        movapd nb314_fjyH2(%esp),%xmm4
        movapd nb314_fjzH2(%esp),%xmm5
        mulpd nb314_dxH2H2(%esp),%xmm0
        mulpd nb314_dyH2H2(%esp),%xmm1
        mulpd nb314_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixH2(%esp),%xmm0
        addpd nb314_fiyH2(%esp),%xmm1
        addpd nb314_fizH2(%esp),%xmm2
        movapd %xmm3,nb314_fjxH2(%esp)
        movapd %xmm4,nb314_fjyH2(%esp)
        movapd %xmm5,nb314_fjzH2(%esp)
        movapd %xmm0,nb314_fixH2(%esp)
        movapd %xmm1,nb314_fiyH2(%esp)
        movapd %xmm2,nb314_fizH2(%esp)

        ## H2-M interaction 
        movapd nb314_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxM(%esp),%xmm3
        movapd nb314_fjyM(%esp),%xmm4
        movapd nb314_fjzM(%esp),%xmm5
        mulpd nb314_dxH2M(%esp),%xmm0
        mulpd nb314_dyH2M(%esp),%xmm1
        mulpd nb314_dzH2M(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixH2(%esp),%xmm0
        addpd nb314_fiyH2(%esp),%xmm1
        addpd nb314_fizH2(%esp),%xmm2
        movapd %xmm3,nb314_fjxM(%esp)
        movapd %xmm4,nb314_fjyM(%esp)
        movapd %xmm5,nb314_fjzM(%esp)
        movapd %xmm0,nb314_fixH2(%esp)
        movapd %xmm1,nb314_fiyH2(%esp)
        movapd %xmm2,nb314_fizH2(%esp)

        ## M-H1 interaction 
        movapd nb314_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxH1(%esp),%xmm3
        movapd nb314_fjyH1(%esp),%xmm4
        movapd nb314_fjzH1(%esp),%xmm5
        mulpd nb314_dxMH1(%esp),%xmm0
        mulpd nb314_dyMH1(%esp),%xmm1
        mulpd nb314_dzMH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixM(%esp),%xmm0
        addpd nb314_fiyM(%esp),%xmm1
        addpd nb314_fizM(%esp),%xmm2
        movapd %xmm3,nb314_fjxH1(%esp)
        movapd %xmm4,nb314_fjyH1(%esp)
        movapd %xmm5,nb314_fjzH1(%esp)
        movapd %xmm0,nb314_fixM(%esp)
        movapd %xmm1,nb314_fiyM(%esp)
        movapd %xmm2,nb314_fizM(%esp)

        ## M-H2 interaction 
        movapd nb314_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxH2(%esp),%xmm3
        movapd nb314_fjyH2(%esp),%xmm4
        movapd nb314_fjzH2(%esp),%xmm5
        mulpd nb314_dxMH2(%esp),%xmm0
        mulpd nb314_dyMH2(%esp),%xmm1
        mulpd nb314_dzMH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixM(%esp),%xmm0
        addpd nb314_fiyM(%esp),%xmm1
        addpd nb314_fizM(%esp),%xmm2
        movapd %xmm3,nb314_fjxH2(%esp)
        movapd %xmm4,nb314_fjyH2(%esp)
        movapd %xmm5,nb314_fjzH2(%esp)
        movapd %xmm0,nb314_fixM(%esp)
        movapd %xmm1,nb314_fiyM(%esp)
        movapd %xmm2,nb314_fizM(%esp)

        ## M-M interaction 
        movapd nb314_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 

        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        mulpd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMM(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addpd  nb314_vctot(%esp),%xmm5
        movapd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb314_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxM(%esp),%xmm3
        movapd nb314_fjyM(%esp),%xmm4
        movapd nb314_fjzM(%esp),%xmm5
        mulpd nb314_dxMM(%esp),%xmm0
        mulpd nb314_dyMM(%esp),%xmm1
        mulpd nb314_dzMM(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb314_fixM(%esp),%xmm0
        addpd nb314_fiyM(%esp),%xmm1
        addpd nb314_fizM(%esp),%xmm2
        movapd %xmm3,nb314_fjxM(%esp)
        movapd %xmm4,nb314_fjyM(%esp)
        movapd %xmm5,nb314_fjzM(%esp)
        movapd %xmm0,nb314_fixM(%esp)
        movapd %xmm1,nb314_fiyM(%esp)
        movapd %xmm2,nb314_fizM(%esp)

        movl nb314_faction(%ebp),%edi

        movd %mm0,%eax
        movd %mm1,%ebx

        ## Did all interactions - now update j forces 
        ## Step1 - transpose fjxO, fjyO and fjzO, fjxH1
        movapd nb314_fjxO(%esp),%xmm0
        movapd nb314_fjyO(%esp),%xmm1
        movapd nb314_fjzO(%esp),%xmm2
        movapd nb314_fjxH1(%esp),%xmm3
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
        movapd nb314_fjyH1(%esp),%xmm0
        movapd nb314_fjzH1(%esp),%xmm1
        movapd nb314_fjxH2(%esp),%xmm2
        movapd nb314_fjyH2(%esp),%xmm3
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
        movapd nb314_fjzH2(%esp),%xmm0
        movapd nb314_fjxM(%esp),%xmm1
        movapd nb314_fjyM(%esp),%xmm2
        movapd nb314_fjzM(%esp),%xmm3
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
        subl $2,nb314_innerk(%esp)
        jl    _nb_kernel314_ia32_sse2.nb314_checksingle
        jmp   _nb_kernel314_ia32_sse2.nb314_unroll_loop
_nb_kernel314_ia32_sse2.nb314_checksingle: 
        movl  nb314_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel314_ia32_sse2.nb314_dosingle
        jmp   _nb_kernel314_ia32_sse2.nb314_updateouterdata
_nb_kernel314_ia32_sse2.nb314_dosingle: 
        movl  nb314_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb314_pos(%ebp),%esi        ## base of pos[] 

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
        movsd  %xmm0,nb314_jxO(%esp)
        movsd  %xmm1,nb314_jzO(%esp)
        movsd  %xmm2,nb314_jyH1(%esp)
        movsd  %xmm3,nb314_jxH2(%esp)
        movsd  %xmm4,nb314_jzH2(%esp)
        movsd  %xmm5,nb314_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb314_jyO(%esp)
        movsd  %xmm1,nb314_jxH1(%esp)
        movsd  %xmm2,nb314_jzH1(%esp)
        movsd  %xmm3,nb314_jyH2(%esp)
        movsd  %xmm4,nb314_jxM(%esp)
        movsd  %xmm5,nb314_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb314_ixO(%esp),%xmm0
        movapd nb314_iyO(%esp),%xmm1
        movapd nb314_izO(%esp),%xmm2
        movapd nb314_ixH1(%esp),%xmm3
        movapd nb314_iyH1(%esp),%xmm4
        movapd nb314_izH1(%esp),%xmm5
        subsd  nb314_jxO(%esp),%xmm0
        subsd  nb314_jyO(%esp),%xmm1
        subsd  nb314_jzO(%esp),%xmm2
        subsd  nb314_jxH1(%esp),%xmm3
        subsd  nb314_jyH1(%esp),%xmm4
        subsd  nb314_jzH1(%esp),%xmm5
        movapd %xmm0,nb314_dxOO(%esp)
        movapd %xmm1,nb314_dyOO(%esp)
        movapd %xmm2,nb314_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH1H1(%esp)
        movapd %xmm4,nb314_dyH1H1(%esp)
        movapd %xmm5,nb314_dzH1H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb314_rsqOO(%esp)
        movapd %xmm3,nb314_rsqH1H1(%esp)

        movapd nb314_ixH1(%esp),%xmm0
        movapd nb314_iyH1(%esp),%xmm1
        movapd nb314_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314_jxH2(%esp),%xmm0
        subsd  nb314_jyH2(%esp),%xmm1
        subsd  nb314_jzH2(%esp),%xmm2
        subsd  nb314_jxM(%esp),%xmm3
        subsd  nb314_jyM(%esp),%xmm4
        subsd  nb314_jzM(%esp),%xmm5
        movapd %xmm0,nb314_dxH1H2(%esp)
        movapd %xmm1,nb314_dyH1H2(%esp)
        movapd %xmm2,nb314_dzH1H2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH1M(%esp)
        movapd %xmm4,nb314_dyH1M(%esp)
        movapd %xmm5,nb314_dzH1M(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb314_rsqH1H2(%esp)
        movapd %xmm3,nb314_rsqH1M(%esp)

        movapd nb314_ixH2(%esp),%xmm0
        movapd nb314_iyH2(%esp),%xmm1
        movapd nb314_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314_jxH1(%esp),%xmm0
        subsd  nb314_jyH1(%esp),%xmm1
        subsd  nb314_jzH1(%esp),%xmm2
        subsd  nb314_jxH2(%esp),%xmm3
        subsd  nb314_jyH2(%esp),%xmm4
        subsd  nb314_jzH2(%esp),%xmm5
        movapd %xmm0,nb314_dxH2H1(%esp)
        movapd %xmm1,nb314_dyH2H1(%esp)
        movapd %xmm2,nb314_dzH2H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxH2H2(%esp)
        movapd %xmm4,nb314_dyH2H2(%esp)
        movapd %xmm5,nb314_dzH2H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb314_rsqH2H1(%esp)
        movapd %xmm3,nb314_rsqH2H2(%esp)

        movapd nb314_ixH2(%esp),%xmm0
        movapd nb314_iyH2(%esp),%xmm1
        movapd nb314_izH2(%esp),%xmm2
        movapd nb314_ixM(%esp),%xmm3
        movapd nb314_iyM(%esp),%xmm4
        movapd nb314_izM(%esp),%xmm5
        subsd  nb314_jxM(%esp),%xmm0
        subsd  nb314_jyM(%esp),%xmm1
        subsd  nb314_jzM(%esp),%xmm2
        subsd  nb314_jxH1(%esp),%xmm3
        subsd  nb314_jyH1(%esp),%xmm4
        subsd  nb314_jzH1(%esp),%xmm5
        movapd %xmm0,nb314_dxH2M(%esp)
        movapd %xmm1,nb314_dyH2M(%esp)
        movapd %xmm2,nb314_dzH2M(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxMH1(%esp)
        movapd %xmm4,nb314_dyMH1(%esp)
        movapd %xmm5,nb314_dzMH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb314_rsqH2M(%esp)
        movapd %xmm4,nb314_rsqMH1(%esp)

        movapd nb314_ixM(%esp),%xmm0
        movapd nb314_iyM(%esp),%xmm1
        movapd nb314_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314_jxH2(%esp),%xmm0
        subsd  nb314_jyH2(%esp),%xmm1
        subsd  nb314_jzH2(%esp),%xmm2
        subsd  nb314_jxM(%esp),%xmm3
        subsd  nb314_jyM(%esp),%xmm4
        subsd  nb314_jzM(%esp),%xmm5
        movapd %xmm0,nb314_dxMH2(%esp)
        movapd %xmm1,nb314_dyMH2(%esp)
        movapd %xmm2,nb314_dzMH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb314_dxMM(%esp)
        movapd %xmm4,nb314_dyMM(%esp)
        movapd %xmm5,nb314_dzMM(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb314_rsqMH2(%esp)
        movapd %xmm4,nb314_rsqMM(%esp)

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
        movapd  nb314_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm3   ## iter1 
        mulsd   nb314_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm1   ## rinv 
        mulsd   nb314_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314_rinvMH2(%esp)
        movapd %xmm5,nb314_rinvMM(%esp)

        movapd nb314_rsqOO(%esp),%xmm0
        movapd nb314_rsqH1H1(%esp),%xmm4
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
        movapd  nb314_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb314_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm1   ## rinv 
        mulsd   nb314_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb314_rinvsqOO(%esp)
        movapd %xmm5,nb314_rinvH1H1(%esp)

        movapd nb314_rsqH1H2(%esp),%xmm0
        movapd nb314_rsqH1M(%esp),%xmm4
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
        movapd  nb314_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm3   ## iter1 
        mulsd   nb314_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm1   ## rinv 
        mulsd   nb314_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314_rinvH1H2(%esp)
        movapd %xmm5,nb314_rinvH1M(%esp)

        movapd nb314_rsqH2H1(%esp),%xmm0
        movapd nb314_rsqH2H2(%esp),%xmm4
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
        movapd  nb314_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm3   ## iter1a 
        mulsd   nb314_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm1   ## rinv 
        mulsd   nb314_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314_rinvH2H1(%esp)
        movapd %xmm5,nb314_rinvH2H2(%esp)

        movapd nb314_rsqMH1(%esp),%xmm0
        movapd nb314_rsqH2M(%esp),%xmm4
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
        movapd  nb314_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm3   ## iter1a 
        mulsd   nb314_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314_half(%esp),%xmm1   ## rinv 
        mulsd   nb314_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314_rinvMH1(%esp)
        movapd %xmm5,nb314_rinvH2M(%esp)

        ## start with OO interaction 
        movsd nb314_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movsd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb314_c6(%esp),%xmm1
        mulsd  nb314_c12(%esp),%xmm2
        movsd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb314_Vvdwtot(%esp),%xmm3
        mulsd  nb314_six(%esp),%xmm1
        mulsd  nb314_twelve(%esp),%xmm2
        subsd  %xmm1,%xmm2
        mulsd  %xmm0,%xmm2
        movsd %xmm3,nb314_Vvdwtot(%esp)

        movapd %xmm2,%xmm0
        movapd %xmm2,%xmm1


        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb314_dxOO(%esp),%xmm0
        mulsd nb314_dyOO(%esp),%xmm1
        mulsd nb314_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixO(%esp),%xmm0
        addsd nb314_fiyO(%esp),%xmm1
        addsd nb314_fizO(%esp),%xmm2
        movsd %xmm3,nb314_fjxO(%esp)
        movsd %xmm4,nb314_fjyO(%esp)
        movsd %xmm5,nb314_fjzO(%esp)
        movsd %xmm0,nb314_fixO(%esp)
        movsd %xmm1,nb314_fiyO(%esp)
        movsd %xmm2,nb314_fizO(%esp)

        ## H1-H1 interaction 
        movapd nb314_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1
        movd %eax,%mm0

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb314_dxH1H1(%esp),%xmm0
        mulsd nb314_dyH1H1(%esp),%xmm1
        mulsd nb314_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixH1(%esp),%xmm0
        addsd nb314_fiyH1(%esp),%xmm1
        addsd nb314_fizH1(%esp),%xmm2
        movsd %xmm3,nb314_fjxH1(%esp)
        movsd %xmm4,nb314_fjyH1(%esp)
        movsd %xmm5,nb314_fjzH1(%esp)
        movsd %xmm0,nb314_fixH1(%esp)
        movsd %xmm1,nb314_fiyH1(%esp)
        movsd %xmm2,nb314_fizH1(%esp)

        ## H1-H2 interaction  
        movapd nb314_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb314_dxH1H2(%esp),%xmm0
        mulsd nb314_dyH1H2(%esp),%xmm1
        mulsd nb314_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixH1(%esp),%xmm0
        addsd nb314_fiyH1(%esp),%xmm1
        addsd nb314_fizH1(%esp),%xmm2
        movsd %xmm3,nb314_fjxH2(%esp)
        movsd %xmm4,nb314_fjyH2(%esp)
        movsd %xmm5,nb314_fjzH2(%esp)
        movsd %xmm0,nb314_fixH1(%esp)
        movsd %xmm1,nb314_fiyH1(%esp)
        movsd %xmm2,nb314_fizH1(%esp)

        ## H1-M interaction 
        movapd nb314_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2


        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb314_dxH1M(%esp),%xmm0
        mulsd nb314_dyH1M(%esp),%xmm1
        mulsd nb314_dzH1M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixH1(%esp),%xmm0
        addsd nb314_fiyH1(%esp),%xmm1
        addsd nb314_fizH1(%esp),%xmm2
        movsd %xmm3,nb314_fjxM(%esp)
        movsd %xmm4,nb314_fjyM(%esp)
        movsd %xmm5,nb314_fjzM(%esp)
        movsd %xmm0,nb314_fixH1(%esp)
        movsd %xmm1,nb314_fiyH1(%esp)
        movsd %xmm2,nb314_fizH1(%esp)

        ## H2-H1 interaction 
        movapd nb314_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxH1(%esp),%xmm3
        movapd nb314_fjyH1(%esp),%xmm4
        movapd nb314_fjzH1(%esp),%xmm5
        mulsd nb314_dxH2H1(%esp),%xmm0
        mulsd nb314_dyH2H1(%esp),%xmm1
        mulsd nb314_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixH2(%esp),%xmm0
        addsd nb314_fiyH2(%esp),%xmm1
        addsd nb314_fizH2(%esp),%xmm2
        movsd %xmm3,nb314_fjxH1(%esp)
        movsd %xmm4,nb314_fjyH1(%esp)
        movsd %xmm5,nb314_fjzH1(%esp)
        movsd %xmm0,nb314_fixH2(%esp)
        movsd %xmm1,nb314_fiyH2(%esp)
        movsd %xmm2,nb314_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb314_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxH2(%esp),%xmm3
        movapd nb314_fjyH2(%esp),%xmm4
        movapd nb314_fjzH2(%esp),%xmm5
        mulsd nb314_dxH2H2(%esp),%xmm0
        mulsd nb314_dyH2H2(%esp),%xmm1
        mulsd nb314_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixH2(%esp),%xmm0
        addsd nb314_fiyH2(%esp),%xmm1
        addsd nb314_fizH2(%esp),%xmm2
        movsd %xmm3,nb314_fjxH2(%esp)
        movsd %xmm4,nb314_fjyH2(%esp)
        movsd %xmm5,nb314_fjzH2(%esp)
        movsd %xmm0,nb314_fixH2(%esp)
        movsd %xmm1,nb314_fiyH2(%esp)
        movsd %xmm2,nb314_fizH2(%esp)

        ## H2-M interaction 
        movapd nb314_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxM(%esp),%xmm3
        movapd nb314_fjyM(%esp),%xmm4
        movapd nb314_fjzM(%esp),%xmm5
        mulsd nb314_dxH2M(%esp),%xmm0
        mulsd nb314_dyH2M(%esp),%xmm1
        mulsd nb314_dzH2M(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixH2(%esp),%xmm0
        addsd nb314_fiyH2(%esp),%xmm1
        addsd nb314_fizH2(%esp),%xmm2
        movsd %xmm3,nb314_fjxM(%esp)
        movsd %xmm4,nb314_fjyM(%esp)
        movsd %xmm5,nb314_fjzM(%esp)
        movsd %xmm0,nb314_fixH2(%esp)
        movsd %xmm1,nb314_fiyH2(%esp)
        movsd %xmm2,nb314_fizH2(%esp)

        ## M-H1 interaction 
        movapd nb314_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxH1(%esp),%xmm3
        movapd nb314_fjyH1(%esp),%xmm4
        movapd nb314_fjzH1(%esp),%xmm5
        mulsd nb314_dxMH1(%esp),%xmm0
        mulsd nb314_dyMH1(%esp),%xmm1
        mulsd nb314_dzMH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixM(%esp),%xmm0
        addsd nb314_fiyM(%esp),%xmm1
        addsd nb314_fizM(%esp),%xmm2
        movsd %xmm3,nb314_fjxH1(%esp)
        movsd %xmm4,nb314_fjyH1(%esp)
        movsd %xmm5,nb314_fjzH1(%esp)
        movsd %xmm0,nb314_fixM(%esp)
        movsd %xmm1,nb314_fiyM(%esp)
        movsd %xmm2,nb314_fizM(%esp)

        ## M-H2 interaction 
        movapd nb314_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxH2(%esp),%xmm3
        movapd nb314_fjyH2(%esp),%xmm4
        movapd nb314_fjzH2(%esp),%xmm5
        mulsd nb314_dxMH2(%esp),%xmm0
        mulsd nb314_dyMH2(%esp),%xmm1
        mulsd nb314_dzMH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixM(%esp),%xmm0
        addsd nb314_fiyM(%esp),%xmm1
        addsd nb314_fizM(%esp),%xmm2
        movsd %xmm3,nb314_fjxH2(%esp)
        movsd %xmm4,nb314_fjyH2(%esp)
        movsd %xmm5,nb314_fjzH2(%esp)
        movsd %xmm0,nb314_fixM(%esp)
        movsd %xmm1,nb314_fiyM(%esp)
        movsd %xmm2,nb314_fizM(%esp)

        ## M-M interaction 
        movapd nb314_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   

        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb314_two(%esp),%xmm7    ## two*Heps2 
        movapd nb314_qqMM(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
        ## at this point mm5 contains vcoul and xmm3 fijC 

        addsd  nb314_vctot(%esp),%xmm5
        movsd %xmm5,nb314_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb314_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb314_fjxM(%esp),%xmm3
        movapd nb314_fjyM(%esp),%xmm4
        movapd nb314_fjzM(%esp),%xmm5
        mulsd nb314_dxMM(%esp),%xmm0
        mulsd nb314_dyMM(%esp),%xmm1
        mulsd nb314_dzMM(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb314_fixM(%esp),%xmm0
        addsd nb314_fiyM(%esp),%xmm1
        addsd nb314_fizM(%esp),%xmm2
        movsd %xmm3,nb314_fjxM(%esp)
        movsd %xmm4,nb314_fjyM(%esp)
        movsd %xmm5,nb314_fjzM(%esp)
        movsd %xmm0,nb314_fixM(%esp)
        movsd %xmm1,nb314_fiyM(%esp)
        movsd %xmm2,nb314_fizM(%esp)

        movl nb314_faction(%ebp),%edi

        movd %mm0,%eax

        ## Did all interactions - now update j forces 
        ## Step1 - merge forces
        movlpd nb314_fjxO(%esp),%xmm0
        movlpd nb314_fjzO(%esp),%xmm1
        movlpd nb314_fjyH1(%esp),%xmm2
        movlpd nb314_fjxH2(%esp),%xmm3
        movlpd nb314_fjzH2(%esp),%xmm4
        movlpd nb314_fjyM(%esp),%xmm5

        movhpd nb314_fjyO(%esp),%xmm0
        movhpd nb314_fjxH1(%esp),%xmm1
        movhpd nb314_fjzH1(%esp),%xmm2
        movhpd nb314_fjyH2(%esp),%xmm3
        movhpd nb314_fjxM(%esp),%xmm4
        movhpd nb314_fjzM(%esp),%xmm5

        movlpd (%edi,%eax,8),%xmm6
        movhpd 8(%edi,%eax,8),%xmm6
        movlpd 16(%edi,%eax,8),%xmm7
        movhpd 24(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm0
        addpd  %xmm7,%xmm1
        movlpd 32(%edi,%eax,8),%xmm6
        movhpd 40(%edi,%eax,8),%xmm6
        movlpd 48(%edi,%eax,8),%xmm7
        movhpd 56(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm2
        addpd  %xmm7,%xmm3
        movlpd 64(%edi,%eax,8),%xmm6
        movhpd 72(%edi,%eax,8),%xmm6
        movlpd 80(%edi,%eax,8),%xmm7
        movhpd 88(%edi,%eax,8),%xmm7
        addpd  %xmm6,%xmm4
        addpd  %xmm7,%xmm5

        movlpd %xmm0,(%edi,%eax,8)
        movhpd %xmm0,8(%edi,%eax,8)
        movlpd %xmm1,16(%edi,%eax,8)
        movhpd %xmm1,24(%edi,%eax,8)
        movlpd %xmm2,32(%edi,%eax,8)
        movhpd %xmm2,40(%edi,%eax,8)
        movlpd %xmm3,48(%edi,%eax,8)
        movhpd %xmm3,56(%edi,%eax,8)
        movlpd %xmm4,64(%edi,%eax,8)
        movhpd %xmm4,72(%edi,%eax,8)
        movlpd %xmm5,80(%edi,%eax,8)
        movhpd %xmm5,88(%edi,%eax,8)

_nb_kernel314_ia32_sse2.nb314_updateouterdata: 
        movl  nb314_ii3(%esp),%ecx
        movl  nb314_faction(%ebp),%edi
        movl  nb314_fshift(%ebp),%esi
        movl  nb314_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb314_fixO(%esp),%xmm0
        movapd nb314_fiyO(%esp),%xmm1
        movapd nb314_fizO(%esp),%xmm2

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
        movapd nb314_fixH1(%esp),%xmm0
        movapd nb314_fiyH1(%esp),%xmm1
        movapd nb314_fizH1(%esp),%xmm2

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
        movapd nb314_fixH2(%esp),%xmm0
        movapd nb314_fiyH2(%esp),%xmm1
        movapd nb314_fizH2(%esp),%xmm2

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
        movapd nb314_fixM(%esp),%xmm0
        movapd nb314_fiyM(%esp),%xmm1
        movapd nb314_fizM(%esp),%xmm2

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
        movl nb314_n(%esp),%esi
        ## get group index for i particle 
        movl  nb314_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb314_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb314_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb314_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb314_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb314_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel314_ia32_sse2.nb314_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb314_n(%esp)
        jmp _nb_kernel314_ia32_sse2.nb314_outer
_nb_kernel314_ia32_sse2.nb314_outerend: 
        ## check if more outer neighborlists remain
        movl  nb314_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel314_ia32_sse2.nb314_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel314_ia32_sse2.nb314_threadloop
_nb_kernel314_ia32_sse2.nb314_end: 
        emms

        movl nb314_nouter(%esp),%eax
        movl nb314_ninner(%esp),%ebx
        movl nb314_outeriter(%ebp),%ecx
        movl nb314_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb314_salign(%esp),%eax
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



.globl nb_kernel314nf_ia32_sse2
.globl _nb_kernel314nf_ia32_sse2
nb_kernel314nf_ia32_sse2:       
_nb_kernel314nf_ia32_sse2:      
.set nb314nf_p_nri, 8
.set nb314nf_iinr, 12
.set nb314nf_jindex, 16
.set nb314nf_jjnr, 20
.set nb314nf_shift, 24
.set nb314nf_shiftvec, 28
.set nb314nf_fshift, 32
.set nb314nf_gid, 36
.set nb314nf_pos, 40
.set nb314nf_faction, 44
.set nb314nf_charge, 48
.set nb314nf_p_facel, 52
.set nb314nf_argkrf, 56
.set nb314nf_argcrf, 60
.set nb314nf_Vc, 64
.set nb314nf_type, 68
.set nb314nf_p_ntype, 72
.set nb314nf_vdwparam, 76
.set nb314nf_Vvdw, 80
.set nb314nf_p_tabscale, 84
.set nb314nf_VFtab, 88
.set nb314nf_invsqrta, 92
.set nb314nf_dvda, 96
.set nb314nf_p_gbtabscale, 100
.set nb314nf_GBtab, 104
.set nb314nf_p_nthreads, 108
.set nb314nf_count, 112
.set nb314nf_mtx, 116
.set nb314nf_outeriter, 120
.set nb314nf_inneriter, 124
.set nb314nf_work, 128
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
.set nb314nf_innerjjnr, 872
.set nb314nf_innerk, 876
.set nb314nf_n, 880
.set nb314nf_nn1, 884
.set nb314nf_nri, 888
.set nb314nf_nouter, 892
.set nb314nf_ninner, 896
.set nb314nf_salign, 900
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
        movl %eax,nb314nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb314nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb314nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb314nf_nouter(%esp)
        movl %eax,nb314nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb314nf_half(%esp)
        movl %ebx,nb314nf_half+4(%esp)
        movsd nb314nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb314nf_half(%esp)
        movapd %xmm3,nb314nf_three(%esp)
        movl nb314nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3

        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb314nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb314nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb314nf_charge(%ebp),%edx
        movsd 24(%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb314nf_p_facel(%ebp),%esi
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
        movapd %xmm3,nb314nf_qqMM(%esp)
        movapd %xmm4,nb314nf_qqMH(%esp)
        movapd %xmm5,nb314nf_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb314nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb314nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb314nf_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movlpd 8(%eax,%edx,8),%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb314nf_c6(%esp)
        movapd %xmm1,nb314nf_c12(%esp)

_nb_kernel314nf_ia32_sse2.nb314nf_threadloop: 
        movl  nb314nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel314nf_ia32_sse2.nb314nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel314nf_ia32_sse2.nb314nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb314nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb314nf_n(%esp)
        movl %ebx,nb314nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel314nf_ia32_sse2.nb314nf_outerstart
        jmp _nb_kernel314nf_ia32_sse2.nb314nf_end

_nb_kernel314nf_ia32_sse2.nb314nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb314nf_nouter(%esp),%ebx
        movl %ebx,nb314nf_nouter(%esp)

_nb_kernel314nf_ia32_sse2.nb314nf_outer: 
        movl  nb314nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb314nf_is3(%esp)            ## store is3 

        movl  nb314nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb314nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        movapd %xmm0,%xmm6
        movapd %xmm1,%xmm7

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb314nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb314nf_ii3(%esp)

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
        movapd %xmm3,nb314nf_ixO(%esp)
        movapd %xmm4,nb314nf_iyO(%esp)
        movapd %xmm5,nb314nf_izO(%esp)
        movapd %xmm6,nb314nf_ixH1(%esp)
        movapd %xmm7,nb314nf_iyH1(%esp)

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
        movapd %xmm6,nb314nf_izH1(%esp)
        movapd %xmm0,nb314nf_ixH2(%esp)
        movapd %xmm1,nb314nf_iyH2(%esp)
        movapd %xmm2,nb314nf_izH2(%esp)
        movapd %xmm3,nb314nf_ixM(%esp)
        movapd %xmm4,nb314nf_iyM(%esp)
        movapd %xmm5,nb314nf_izM(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb314nf_vctot(%esp)
        movapd %xmm4,nb314nf_Vvdwtot(%esp)

        movl  nb314nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb314nf_pos(%ebp),%esi
        movl  nb314nf_faction(%ebp),%edi
        movl  nb314nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb314nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb314nf_ninner(%esp),%ecx
        movl  %ecx,nb314nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb314nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel314nf_ia32_sse2.nb314nf_unroll_loop
        jmp   _nb_kernel314nf_ia32_sse2.nb314nf_checksingle
_nb_kernel314nf_ia32_sse2.nb314nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb314nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb314nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb314nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm0,nb314nf_jxO(%esp)
        movapd  %xmm1,nb314nf_jyO(%esp)
        movapd  %xmm3,nb314nf_jzO(%esp)
        movapd  %xmm4,nb314nf_jxH1(%esp)

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
        movapd  %xmm0,nb314nf_jyH1(%esp)
        movapd  %xmm1,nb314nf_jzH1(%esp)
        movapd  %xmm3,nb314nf_jxH2(%esp)
        movapd  %xmm4,nb314nf_jyH2(%esp)

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
        movapd  %xmm0,nb314nf_jzH2(%esp)
        movapd  %xmm1,nb314nf_jxM(%esp)
        movapd  %xmm3,nb314nf_jyM(%esp)
        movapd  %xmm4,nb314nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb314nf_ixO(%esp),%xmm0
        movapd nb314nf_iyO(%esp),%xmm1
        movapd nb314nf_izO(%esp),%xmm2
        movapd nb314nf_ixH1(%esp),%xmm3
        movapd nb314nf_iyH1(%esp),%xmm4
        movapd nb314nf_izH1(%esp),%xmm5
        subpd  nb314nf_jxO(%esp),%xmm0
        subpd  nb314nf_jyO(%esp),%xmm1
        subpd  nb314nf_jzO(%esp),%xmm2
        subpd  nb314nf_jxH1(%esp),%xmm3
        subpd  nb314nf_jyH1(%esp),%xmm4
        subpd  nb314nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqOO(%esp)
        movapd %xmm3,nb314nf_rsqH1H1(%esp)

        movapd nb314nf_ixH1(%esp),%xmm0
        movapd nb314nf_iyH1(%esp),%xmm1
        movapd nb314nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314nf_jxH2(%esp),%xmm0
        subpd  nb314nf_jyH2(%esp),%xmm1
        subpd  nb314nf_jzH2(%esp),%xmm2
        subpd  nb314nf_jxM(%esp),%xmm3
        subpd  nb314nf_jyM(%esp),%xmm4
        subpd  nb314nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqH1H2(%esp)
        movapd %xmm3,nb314nf_rsqH1M(%esp)

        movapd nb314nf_ixH2(%esp),%xmm0
        movapd nb314nf_iyH2(%esp),%xmm1
        movapd nb314nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314nf_jxH1(%esp),%xmm0
        subpd  nb314nf_jyH1(%esp),%xmm1
        subpd  nb314nf_jzH1(%esp),%xmm2
        subpd  nb314nf_jxH2(%esp),%xmm3
        subpd  nb314nf_jyH2(%esp),%xmm4
        subpd  nb314nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqH2H1(%esp)
        movapd %xmm3,nb314nf_rsqH2H2(%esp)

        movapd nb314nf_ixH2(%esp),%xmm0
        movapd nb314nf_iyH2(%esp),%xmm1
        movapd nb314nf_izH2(%esp),%xmm2
        movapd nb314nf_ixM(%esp),%xmm3
        movapd nb314nf_iyM(%esp),%xmm4
        movapd nb314nf_izM(%esp),%xmm5
        subpd  nb314nf_jxM(%esp),%xmm0
        subpd  nb314nf_jyM(%esp),%xmm1
        subpd  nb314nf_jzM(%esp),%xmm2
        subpd  nb314nf_jxH1(%esp),%xmm3
        subpd  nb314nf_jyH1(%esp),%xmm4
        subpd  nb314nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqH2M(%esp)
        movapd %xmm4,nb314nf_rsqMH1(%esp)

        movapd nb314nf_ixM(%esp),%xmm0
        movapd nb314nf_iyM(%esp),%xmm1
        movapd nb314nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subpd  nb314nf_jxH2(%esp),%xmm0
        subpd  nb314nf_jyH2(%esp),%xmm1
        subpd  nb314nf_jzH2(%esp),%xmm2
        subpd  nb314nf_jxM(%esp),%xmm3
        subpd  nb314nf_jyM(%esp),%xmm4
        subpd  nb314nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqMH2(%esp)
        movapd %xmm4,nb314nf_rsqMM(%esp)

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
        movapd  nb314nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb314nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvMH2(%esp)
        movapd %xmm5,nb314nf_rinvMM(%esp)

        movapd nb314nf_rsqOO(%esp),%xmm0
        movapd nb314nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb314nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb314nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%esp),%xmm5   ## rinv
        mulpd   %xmm1,%xmm1
        movapd %xmm1,nb314nf_rinvsqOO(%esp)
        movapd %xmm5,nb314nf_rinvH1H1(%esp)

        movapd nb314nf_rsqH1H2(%esp),%xmm0
        movapd nb314nf_rsqH1M(%esp),%xmm4
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
        movapd  nb314nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb314nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvH1H2(%esp)
        movapd %xmm5,nb314nf_rinvH1M(%esp)

        movapd nb314nf_rsqH2H1(%esp),%xmm0
        movapd nb314nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb314nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb314nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvH2H1(%esp)
        movapd %xmm5,nb314nf_rinvH2H2(%esp)

        movapd nb314nf_rsqMH1(%esp),%xmm0
        movapd nb314nf_rsqH2M(%esp),%xmm4
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
        movapd  nb314nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb314nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb314nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvMH1(%esp)
        movapd %xmm5,nb314nf_rinvH2M(%esp)

        ## start with OO interaction 
        movapd nb314nf_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movapd  %xmm0,%xmm1
        mulpd   %xmm1,%xmm1 ## rinv4
        mulpd   %xmm0,%xmm1 ##rinvsix
        movapd  %xmm1,%xmm2
        mulpd   %xmm2,%xmm2 ## rinvtwelve
        mulpd  nb314nf_c6(%esp),%xmm1
        mulpd  nb314nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subpd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addpd  nb314nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb314nf_Vvdwtot(%esp)

        ## H1-H1 interaction 
        movapd nb314nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1

        movd %eax,%mm0
        movd %ebx,%mm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## H1-H2 interaction  
        movapd nb314nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## H1-M interaction 
        movapd nb314nf_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb314nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb314nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## H2-M interaction 
        movapd nb314nf_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## M-H1 interaction 
        movapd nb314nf_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## M-H2 interaction 
        movapd nb314nf_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## M-M interaction 
        movapd nb314nf_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb314nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulpd  nb314nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi
        movd %mm6,%eax
        psrlq $32,%mm6
        movd %mm6,%ebx          ## indices in eax/ebx 

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movlpd (%esi,%ebx,8),%xmm3      ## Y2
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        movhpd 8(%esi,%ebx,8),%xmm3     ## Y2 F2 

        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 Y2 
        unpckhpd %xmm3,%xmm5    ## F1 F2 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1
        movlpd 16(%esi,%ebx,8),%xmm3    ## G2
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        movhpd 24(%esi,%ebx,8),%xmm3    ## G2 H2 
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 G2 
        unpckhpd %xmm3,%xmm7    ## H1 H2 
        ## coulomb table ready, in xmm4-xmm7            
        mulpd  %xmm1,%xmm6      ## xmm6=Geps 
        mulpd  %xmm2,%xmm7      ## xmm7=Heps2 
        addpd  %xmm6,%xmm5
        addpd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMM(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addpd  nb314nf_vctot(%esp),%xmm5
        movapd %xmm5,nb314nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb314nf_innerk(%esp)
        jl    _nb_kernel314nf_ia32_sse2.nb314nf_checksingle
        jmp   _nb_kernel314nf_ia32_sse2.nb314nf_unroll_loop
_nb_kernel314nf_ia32_sse2.nb314nf_checksingle: 
        movl  nb314nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel314nf_ia32_sse2.nb314nf_dosingle
        jmp   _nb_kernel314nf_ia32_sse2.nb314nf_updateouterdata
_nb_kernel314nf_ia32_sse2.nb314nf_dosingle: 
        movl  nb314nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb314nf_pos(%ebp),%esi        ## base of pos[] 

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
        movsd  %xmm0,nb314nf_jxO(%esp)
        movsd  %xmm1,nb314nf_jzO(%esp)
        movsd  %xmm2,nb314nf_jyH1(%esp)
        movsd  %xmm3,nb314nf_jxH2(%esp)
        movsd  %xmm4,nb314nf_jzH2(%esp)
        movsd  %xmm5,nb314nf_jyM(%esp)
        unpckhpd %xmm0,%xmm0
        unpckhpd %xmm1,%xmm1
        unpckhpd %xmm2,%xmm2
        unpckhpd %xmm3,%xmm3
        unpckhpd %xmm4,%xmm4
        unpckhpd %xmm5,%xmm5
        movsd  %xmm0,nb314nf_jyO(%esp)
        movsd  %xmm1,nb314nf_jxH1(%esp)
        movsd  %xmm2,nb314nf_jzH1(%esp)
        movsd  %xmm3,nb314nf_jyH2(%esp)
        movsd  %xmm4,nb314nf_jxM(%esp)
        movsd  %xmm5,nb314nf_jzM(%esp)

        ## start calculating pairwise distances
        movapd nb314nf_ixO(%esp),%xmm0
        movapd nb314nf_iyO(%esp),%xmm1
        movapd nb314nf_izO(%esp),%xmm2
        movapd nb314nf_ixH1(%esp),%xmm3
        movapd nb314nf_iyH1(%esp),%xmm4
        movapd nb314nf_izH1(%esp),%xmm5
        subsd  nb314nf_jxO(%esp),%xmm0
        subsd  nb314nf_jyO(%esp),%xmm1
        subsd  nb314nf_jzO(%esp),%xmm2
        subsd  nb314nf_jxH1(%esp),%xmm3
        subsd  nb314nf_jyH1(%esp),%xmm4
        subsd  nb314nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqOO(%esp)
        movapd %xmm3,nb314nf_rsqH1H1(%esp)

        movapd nb314nf_ixH1(%esp),%xmm0
        movapd nb314nf_iyH1(%esp),%xmm1
        movapd nb314nf_izH1(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314nf_jxH2(%esp),%xmm0
        subsd  nb314nf_jyH2(%esp),%xmm1
        subsd  nb314nf_jzH2(%esp),%xmm2
        subsd  nb314nf_jxM(%esp),%xmm3
        subsd  nb314nf_jyM(%esp),%xmm4
        subsd  nb314nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqH1H2(%esp)
        movapd %xmm3,nb314nf_rsqH1M(%esp)

        movapd nb314nf_ixH2(%esp),%xmm0
        movapd nb314nf_iyH2(%esp),%xmm1
        movapd nb314nf_izH2(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314nf_jxH1(%esp),%xmm0
        subsd  nb314nf_jyH1(%esp),%xmm1
        subsd  nb314nf_jzH1(%esp),%xmm2
        subsd  nb314nf_jxH2(%esp),%xmm3
        subsd  nb314nf_jyH2(%esp),%xmm4
        subsd  nb314nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqH2H1(%esp)
        movapd %xmm3,nb314nf_rsqH2H2(%esp)

        movapd nb314nf_ixH2(%esp),%xmm0
        movapd nb314nf_iyH2(%esp),%xmm1
        movapd nb314nf_izH2(%esp),%xmm2
        movapd nb314nf_ixM(%esp),%xmm3
        movapd nb314nf_iyM(%esp),%xmm4
        movapd nb314nf_izM(%esp),%xmm5
        subsd  nb314nf_jxM(%esp),%xmm0
        subsd  nb314nf_jyM(%esp),%xmm1
        subsd  nb314nf_jzM(%esp),%xmm2
        subsd  nb314nf_jxH1(%esp),%xmm3
        subsd  nb314nf_jyH1(%esp),%xmm4
        subsd  nb314nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqH2M(%esp)
        movapd %xmm4,nb314nf_rsqMH1(%esp)

        movapd nb314nf_ixM(%esp),%xmm0
        movapd nb314nf_iyM(%esp),%xmm1
        movapd nb314nf_izM(%esp),%xmm2
        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        subsd  nb314nf_jxH2(%esp),%xmm0
        subsd  nb314nf_jyH2(%esp),%xmm1
        subsd  nb314nf_jzH2(%esp),%xmm2
        subsd  nb314nf_jxM(%esp),%xmm3
        subsd  nb314nf_jyM(%esp),%xmm4
        subsd  nb314nf_jzM(%esp),%xmm5
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
        movapd %xmm0,nb314nf_rsqMH2(%esp)
        movapd %xmm4,nb314nf_rsqMM(%esp)

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
        movapd  nb314nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb314nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvMH2(%esp)
        movapd %xmm5,nb314nf_rinvMM(%esp)

        movapd nb314nf_rsqOO(%esp),%xmm0
        movapd nb314nf_rsqH1H1(%esp),%xmm4
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
        movapd  nb314nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb314nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%esp),%xmm5   ## rinv
        mulsd   %xmm1,%xmm1
        movapd %xmm1,nb314nf_rinvsqOO(%esp)
        movapd %xmm5,nb314nf_rinvH1H1(%esp)

        movapd nb314nf_rsqH1H2(%esp),%xmm0
        movapd nb314nf_rsqH1M(%esp),%xmm4
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
        movapd  nb314nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb314nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvH1H2(%esp)
        movapd %xmm5,nb314nf_rinvH1M(%esp)

        movapd nb314nf_rsqH2H1(%esp),%xmm0
        movapd nb314nf_rsqH2H2(%esp),%xmm4
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
        movapd  nb314nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb314nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvH2H1(%esp)
        movapd %xmm5,nb314nf_rinvH2H2(%esp)

        movapd nb314nf_rsqMH1(%esp),%xmm0
        movapd nb314nf_rsqH2M(%esp),%xmm4
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
        movapd  nb314nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb314nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb314nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb314nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb314nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb314nf_rinvMH1(%esp)
        movapd %xmm5,nb314nf_rinvH2M(%esp)


        ## start with OO interaction 
        movsd nb314nf_rinvsqOO(%esp),%xmm0   ## xmm0=rinvsq
        movsd  %xmm0,%xmm1
        mulsd   %xmm1,%xmm1 ## rinv4
        mulsd   %xmm0,%xmm1 ##rinvsix
        movsd  %xmm1,%xmm2
        mulsd   %xmm2,%xmm2 ## rinvtwelve
        mulsd  nb314nf_c6(%esp),%xmm1
        mulsd  nb314nf_c12(%esp),%xmm2
        movsd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## Vvdw=Vvdw12-Vvdw6            
        addsd  nb314nf_Vvdwtot(%esp),%xmm3
        movsd %xmm3,nb314nf_Vvdwtot(%esp)

        ## H1-H1 interaction 
        movapd nb314nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1

        movd %eax,%mm0

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  


        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

        ## H1-H2 interaction  
        movapd nb314nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

        ## H1-M interaction 
        movapd nb314nf_rinvH1M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH1M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb314nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb314nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

        ## H2-M interaction 
        movapd nb314nf_rinvH2M(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqH2M(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

        ## M-H1 interaction 
        movapd nb314nf_rinvMH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqMH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

        ## M-H2 interaction 
        movapd nb314nf_rinvMH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqMH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

        ## M-M interaction 
        movapd nb314nf_rinvMM(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb314nf_rsqMM(%esp),%xmm1   ## xmm1=r 
        mulsd  nb314nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb314nf_VFtab(%ebp),%esi

        movsd (%esi,%eax,8),%xmm4       ## Y1   
        movsd 8(%esi,%eax,8),%xmm5      ## F1   
        movsd 16(%esi,%eax,8),%xmm6     ## G1   
        movsd 24(%esi,%eax,8),%xmm7     ## H1   
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb314nf_qqMM(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  

        addsd  nb314nf_vctot(%esp),%xmm5
        movsd %xmm5,nb314nf_vctot(%esp)

_nb_kernel314nf_ia32_sse2.nb314nf_updateouterdata: 
        ## get n from stack
        movl nb314nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb314nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb314nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb314nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb314nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb314nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb314nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel314nf_ia32_sse2.nb314nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb314nf_n(%esp)
        jmp _nb_kernel314nf_ia32_sse2.nb314nf_outer
_nb_kernel314nf_ia32_sse2.nb314nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb314nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel314nf_ia32_sse2.nb314nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel314nf_ia32_sse2.nb314nf_threadloop
_nb_kernel314nf_ia32_sse2.nb314nf_end: 
        emms

        movl nb314nf_nouter(%esp),%eax
        movl nb314nf_ninner(%esp),%ebx
        movl nb314nf_outeriter(%ebp),%ecx
        movl nb314nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb314nf_salign(%esp),%eax
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

