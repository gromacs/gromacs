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


.globl nb_kernel212_ia32_sse2
.globl _nb_kernel212_ia32_sse2
nb_kernel212_ia32_sse2: 
_nb_kernel212_ia32_sse2:        
.set nb212_p_nri, 8
.set nb212_iinr, 12
.set nb212_jindex, 16
.set nb212_jjnr, 20
.set nb212_shift, 24
.set nb212_shiftvec, 28
.set nb212_fshift, 32
.set nb212_gid, 36
.set nb212_pos, 40
.set nb212_faction, 44
.set nb212_charge, 48
.set nb212_p_facel, 52
.set nb212_argkrf, 56
.set nb212_argcrf, 60
.set nb212_Vc, 64
.set nb212_type, 68
.set nb212_p_ntype, 72
.set nb212_vdwparam, 76
.set nb212_Vvdw, 80
.set nb212_p_tabscale, 84
.set nb212_VFtab, 88
.set nb212_invsqrta, 92
.set nb212_dvda, 96
.set nb212_p_gbtabscale, 100
.set nb212_GBtab, 104
.set nb212_p_nthreads, 108
.set nb212_count, 112
.set nb212_mtx, 116
.set nb212_outeriter, 120
.set nb212_inneriter, 124
.set nb212_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb212_ixO, 0
.set nb212_iyO, 16
.set nb212_izO, 32
.set nb212_ixH1, 48
.set nb212_iyH1, 64
.set nb212_izH1, 80
.set nb212_ixH2, 96
.set nb212_iyH2, 112
.set nb212_izH2, 128
.set nb212_jxO, 144
.set nb212_jyO, 160
.set nb212_jzO, 176
.set nb212_jxH1, 192
.set nb212_jyH1, 208
.set nb212_jzH1, 224
.set nb212_jxH2, 240
.set nb212_jyH2, 256
.set nb212_jzH2, 272
.set nb212_dxOO, 288
.set nb212_dyOO, 304
.set nb212_dzOO, 320
.set nb212_dxOH1, 336
.set nb212_dyOH1, 352
.set nb212_dzOH1, 368
.set nb212_dxOH2, 384
.set nb212_dyOH2, 400
.set nb212_dzOH2, 416
.set nb212_dxH1O, 432
.set nb212_dyH1O, 448
.set nb212_dzH1O, 464
.set nb212_dxH1H1, 480
.set nb212_dyH1H1, 496
.set nb212_dzH1H1, 512
.set nb212_dxH1H2, 528
.set nb212_dyH1H2, 544
.set nb212_dzH1H2, 560
.set nb212_dxH2O, 576
.set nb212_dyH2O, 592
.set nb212_dzH2O, 608
.set nb212_dxH2H1, 624
.set nb212_dyH2H1, 640
.set nb212_dzH2H1, 656
.set nb212_dxH2H2, 672
.set nb212_dyH2H2, 688
.set nb212_dzH2H2, 704
.set nb212_qqOO, 720
.set nb212_qqOH, 736
.set nb212_qqHH, 752
.set nb212_c6, 768
.set nb212_c12, 784
.set nb212_six, 800
.set nb212_twelve, 816
.set nb212_vctot, 832
.set nb212_Vvdwtot, 848
.set nb212_fixO, 864
.set nb212_fiyO, 880
.set nb212_fizO, 896
.set nb212_fixH1, 912
.set nb212_fiyH1, 928
.set nb212_fizH1, 944
.set nb212_fixH2, 960
.set nb212_fiyH2, 976
.set nb212_fizH2, 992
.set nb212_fjxO, 1008
.set nb212_fjyO, 1024
.set nb212_fjzO, 1040
.set nb212_fjxH1, 1056
.set nb212_fjyH1, 1072
.set nb212_fjzH1, 1088
.set nb212_fjxH2, 1104
.set nb212_fjyH2, 1120
.set nb212_fjzH2, 1136
.set nb212_half, 1152
.set nb212_three, 1168
.set nb212_rsqOO, 1184
.set nb212_rsqOH1, 1200
.set nb212_rsqOH2, 1216
.set nb212_rsqH1O, 1232
.set nb212_rsqH1H1, 1248
.set nb212_rsqH1H2, 1264
.set nb212_rsqH2O, 1280
.set nb212_rsqH2H1, 1296
.set nb212_rsqH2H2, 1312
.set nb212_rinvOO, 1328
.set nb212_rinvOH1, 1344
.set nb212_rinvOH2, 1360
.set nb212_rinvH1O, 1376
.set nb212_rinvH1H1, 1392
.set nb212_rinvH1H2, 1408
.set nb212_rinvH2O, 1424
.set nb212_rinvH2H1, 1440
.set nb212_rinvH2H2, 1456
.set nb212_two, 1472
.set nb212_krf, 1488
.set nb212_crf, 1504
.set nb212_is3, 1520
.set nb212_ii3, 1524
.set nb212_innerjjnr, 1528
.set nb212_innerk, 1532
.set nb212_n, 1536
.set nb212_nn1, 1540
.set nb212_nri, 1544
.set nb212_nouter, 1548
.set nb212_ninner, 1552
.set nb212_salign, 1556
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $1560,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb212_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb212_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb212_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb212_nouter(%esp)
        movl %eax,nb212_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb212_half(%esp)
        movl %ebx,nb212_half+4(%esp)
        movsd nb212_half(%esp),%xmm1
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
        movapd %xmm1,nb212_half(%esp)
        movapd %xmm2,nb212_two(%esp)
        movapd %xmm3,nb212_three(%esp)
        movapd %xmm4,nb212_six(%esp)
        movapd %xmm5,nb212_twelve(%esp)

        movl nb212_argkrf(%ebp),%esi
        movl nb212_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb212_krf(%esp)
        movapd %xmm6,nb212_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb212_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb212_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb212_p_facel(%ebp),%esi
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
        movapd %xmm3,nb212_qqOO(%esp)
        movapd %xmm4,nb212_qqOH(%esp)
        movapd %xmm5,nb212_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb212_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb212_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb212_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb212_c6(%esp)
        movapd %xmm1,nb212_c12(%esp)

_nb_kernel212_ia32_sse2.nb212_threadloop: 
        movl  nb212_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel212_ia32_sse2.nb212_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel212_ia32_sse2.nb212_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb212_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb212_n(%esp)
        movl %ebx,nb212_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel212_ia32_sse2.nb212_outerstart
        jmp _nb_kernel212_ia32_sse2.nb212_end

_nb_kernel212_ia32_sse2.nb212_outerstart: 
        ## ebx contains number of outer iterations
        addl nb212_nouter(%esp),%ebx
        movl %ebx,nb212_nouter(%esp)

_nb_kernel212_ia32_sse2.nb212_outer: 
        movl  nb212_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb212_is3(%esp)      ## store is3 

        movl  nb212_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movlpd (%eax,%ebx,8),%xmm0
        movlpd 8(%eax,%ebx,8),%xmm1
        movlpd 16(%eax,%ebx,8),%xmm2

        movl  nb212_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb212_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb212_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb212_ixO(%esp)
        movapd %xmm4,nb212_iyO(%esp)
        movapd %xmm5,nb212_izO(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm0
        addsd 32(%eax,%ebx,8),%xmm1
        addsd 40(%eax,%ebx,8),%xmm2
        addsd 48(%eax,%ebx,8),%xmm3
        addsd 56(%eax,%ebx,8),%xmm4
        addsd 64(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb212_ixH1(%esp)
        movapd %xmm1,nb212_iyH1(%esp)
        movapd %xmm2,nb212_izH1(%esp)
        movapd %xmm3,nb212_ixH2(%esp)
        movapd %xmm4,nb212_iyH2(%esp)
        movapd %xmm5,nb212_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb212_vctot(%esp)
        movapd %xmm4,nb212_Vvdwtot(%esp)
        movapd %xmm4,nb212_fixO(%esp)
        movapd %xmm4,nb212_fiyO(%esp)
        movapd %xmm4,nb212_fizO(%esp)
        movapd %xmm4,nb212_fixH1(%esp)
        movapd %xmm4,nb212_fiyH1(%esp)
        movapd %xmm4,nb212_fizH1(%esp)
        movapd %xmm4,nb212_fixH2(%esp)
        movapd %xmm4,nb212_fiyH2(%esp)
        movapd %xmm4,nb212_fizH2(%esp)

        movl  nb212_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb212_pos(%ebp),%esi
        movl  nb212_faction(%ebp),%edi
        movl  nb212_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb212_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb212_ninner(%esp),%ecx
        movl  %ecx,nb212_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb212_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel212_ia32_sse2.nb212_unroll_loop
        jmp   _nb_kernel212_ia32_sse2.nb212_checksingle
_nb_kernel212_ia32_sse2.nb212_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb212_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb212_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb212_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movhpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%ebx,8),%xmm3
        movhpd 16(%esi,%ebx,8),%xmm4
        movhpd 24(%esi,%ebx,8),%xmm5
        movhpd 32(%esi,%ebx,8),%xmm6
        movhpd 40(%esi,%ebx,8),%xmm7
        movapd  %xmm2,nb212_jxO(%esp)
        movapd  %xmm3,nb212_jyO(%esp)
        movapd  %xmm4,nb212_jzO(%esp)
        movapd  %xmm5,nb212_jxH1(%esp)
        movapd  %xmm6,nb212_jyH1(%esp)
        movapd  %xmm7,nb212_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb212_jxH2(%esp)
        movapd  %xmm3,nb212_jyH2(%esp)
        movapd  %xmm4,nb212_jzH2(%esp)

        movapd nb212_ixO(%esp),%xmm0
        movapd nb212_iyO(%esp),%xmm1
        movapd nb212_izO(%esp),%xmm2
        movapd nb212_ixO(%esp),%xmm3
        movapd nb212_iyO(%esp),%xmm4
        movapd nb212_izO(%esp),%xmm5
        subpd  nb212_jxO(%esp),%xmm0
        subpd  nb212_jyO(%esp),%xmm1
        subpd  nb212_jzO(%esp),%xmm2
        subpd  nb212_jxH1(%esp),%xmm3
        subpd  nb212_jyH1(%esp),%xmm4
        subpd  nb212_jzH1(%esp),%xmm5
        movapd %xmm0,nb212_dxOO(%esp)
        movapd %xmm1,nb212_dyOO(%esp)
        movapd %xmm2,nb212_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb212_dxOH1(%esp)
        movapd %xmm4,nb212_dyOH1(%esp)
        movapd %xmm5,nb212_dzOH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb212_rsqOO(%esp)
        movapd %xmm3,nb212_rsqOH1(%esp)

        movapd nb212_ixO(%esp),%xmm0
        movapd nb212_iyO(%esp),%xmm1
        movapd nb212_izO(%esp),%xmm2
        movapd nb212_ixH1(%esp),%xmm3
        movapd nb212_iyH1(%esp),%xmm4
        movapd nb212_izH1(%esp),%xmm5
        subpd  nb212_jxH2(%esp),%xmm0
        subpd  nb212_jyH2(%esp),%xmm1
        subpd  nb212_jzH2(%esp),%xmm2
        subpd  nb212_jxO(%esp),%xmm3
        subpd  nb212_jyO(%esp),%xmm4
        subpd  nb212_jzO(%esp),%xmm5
        movapd %xmm0,nb212_dxOH2(%esp)
        movapd %xmm1,nb212_dyOH2(%esp)
        movapd %xmm2,nb212_dzOH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb212_dxH1O(%esp)
        movapd %xmm4,nb212_dyH1O(%esp)
        movapd %xmm5,nb212_dzH1O(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb212_rsqOH2(%esp)
        movapd %xmm3,nb212_rsqH1O(%esp)

        movapd nb212_ixH1(%esp),%xmm0
        movapd nb212_iyH1(%esp),%xmm1
        movapd nb212_izH1(%esp),%xmm2
        movapd nb212_ixH1(%esp),%xmm3
        movapd nb212_iyH1(%esp),%xmm4
        movapd nb212_izH1(%esp),%xmm5
        subpd  nb212_jxH1(%esp),%xmm0
        subpd  nb212_jyH1(%esp),%xmm1
        subpd  nb212_jzH1(%esp),%xmm2
        subpd  nb212_jxH2(%esp),%xmm3
        subpd  nb212_jyH2(%esp),%xmm4
        subpd  nb212_jzH2(%esp),%xmm5
        movapd %xmm0,nb212_dxH1H1(%esp)
        movapd %xmm1,nb212_dyH1H1(%esp)
        movapd %xmm2,nb212_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb212_dxH1H2(%esp)
        movapd %xmm4,nb212_dyH1H2(%esp)
        movapd %xmm5,nb212_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb212_rsqH1H1(%esp)
        movapd %xmm3,nb212_rsqH1H2(%esp)

        movapd nb212_ixH2(%esp),%xmm0
        movapd nb212_iyH2(%esp),%xmm1
        movapd nb212_izH2(%esp),%xmm2
        movapd nb212_ixH2(%esp),%xmm3
        movapd nb212_iyH2(%esp),%xmm4
        movapd nb212_izH2(%esp),%xmm5
        subpd  nb212_jxO(%esp),%xmm0
        subpd  nb212_jyO(%esp),%xmm1
        subpd  nb212_jzO(%esp),%xmm2
        subpd  nb212_jxH1(%esp),%xmm3
        subpd  nb212_jyH1(%esp),%xmm4
        subpd  nb212_jzH1(%esp),%xmm5
        movapd %xmm0,nb212_dxH2O(%esp)
        movapd %xmm1,nb212_dyH2O(%esp)
        movapd %xmm2,nb212_dzH2O(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb212_dxH2H1(%esp)
        movapd %xmm4,nb212_dyH2H1(%esp)
        movapd %xmm5,nb212_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb212_rsqH2O(%esp)
        movapd %xmm4,nb212_rsqH2H1(%esp)

        movapd nb212_ixH2(%esp),%xmm0
        movapd nb212_iyH2(%esp),%xmm1
        movapd nb212_izH2(%esp),%xmm2
        subpd  nb212_jxH2(%esp),%xmm0
        subpd  nb212_jyH2(%esp),%xmm1
        subpd  nb212_jzH2(%esp),%xmm2
        movapd %xmm0,nb212_dxH2H2(%esp)
        movapd %xmm1,nb212_dyH2H2(%esp)
        movapd %xmm2,nb212_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb212_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
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
        movapd  nb212_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212_half(%esp),%xmm3   ## iter1 
        mulpd   nb212_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212_half(%esp),%xmm1   ## rinv 
        mulpd   nb212_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212_rinvH2H2(%esp)
        movapd %xmm5,nb212_rinvH2H1(%esp)

        movapd nb212_rsqOO(%esp),%xmm0
        movapd nb212_rsqOH1(%esp),%xmm4
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
        movapd  nb212_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb212_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212_half(%esp),%xmm1   ## rinv 
        mulpd   nb212_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb212_rinvOO(%esp)
        movapd %xmm5,nb212_rinvOH1(%esp)

        movapd nb212_rsqOH2(%esp),%xmm0
        movapd nb212_rsqH1O(%esp),%xmm4
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
        movapd  nb212_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212_half(%esp),%xmm3   ## iter1 
        mulpd   nb212_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212_half(%esp),%xmm1   ## rinv 
        mulpd   nb212_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212_rinvOH2(%esp)
        movapd %xmm5,nb212_rinvH1O(%esp)

        movapd nb212_rsqH1H1(%esp),%xmm0
        movapd nb212_rsqH1H2(%esp),%xmm4
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
        movapd  nb212_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212_half(%esp),%xmm3   ## iter1a 
        mulpd   nb212_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212_half(%esp),%xmm1   ## rinv 
        mulpd   nb212_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212_rinvH1H1(%esp)
        movapd %xmm5,nb212_rinvH1H2(%esp)

        movapd nb212_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb212_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb212_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb212_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb212_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb212_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb212_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7              ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0              ## xmm0=rinvsq 
        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1              ## rinvsq*rinvsq 
        mulpd  %xmm0,%xmm1              ## xmm1=rinvsix 
        mulpd  nb212_rsqOO(%esp),%xmm5          ## xmm5=krsq 
        movapd %xmm5,%xmm6              ## krsq 
        addpd  %xmm7,%xmm6              ## xmm6=rinv+ krsq 
        subpd  nb212_crf(%esp),%xmm6    ## rinv+krsq-crf 

        mulpd  nb212_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb212_two(%esp),%xmm5    ## 2*krsq 
        subpd  %xmm5,%xmm7              ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqOO(%esp),%xmm7   ## xmm7 = qq*(rinv-2*krsq) 

        movapd %xmm1,%xmm2              ## rinv6 
        mulpd  %xmm2,%xmm2              ## xmm2=rinvtwelve 
        mulpd  nb212_c6(%esp),%xmm1     ## c6*rinv6 
        mulpd  nb212_c12(%esp),%xmm2    ## c12*rinv12 
        movapd %xmm2,%xmm3              ## c12*rinv12 
        subpd  %xmm1,%xmm3              ## Vvdw12-Vvdw6 
        addpd  nb212_Vvdwtot(%esp),%xmm3
        mulpd  nb212_six(%esp),%xmm1    ## 6.0*Vvdw6 
        mulpd  nb212_twelve(%esp),%xmm2         ## 12*Vvdw12 
        movapd %xmm3,nb212_Vvdwtot(%esp)
        subpd  %xmm1,%xmm2              ## 12*Vvdw12-6*Vvdw6 
        addpd  %xmm7,%xmm2              ## 12*Vvdw12-6*Vvdw6+qq*(rinv-2*krsq) 
        addpd  nb212_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm2,%xmm0              ## (12*Vvdw12-6*Vvdw6+qq*(rinv-2*krsq))*rinvsq 

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb212_dxOO(%esp),%xmm0
        mulpd nb212_dyOO(%esp),%xmm1
        mulpd nb212_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixO(%esp),%xmm0
        addpd nb212_fiyO(%esp),%xmm1
        addpd nb212_fizO(%esp),%xmm2
        movapd %xmm3,nb212_fjxO(%esp)
        movapd %xmm4,nb212_fjyO(%esp)
        movapd %xmm5,nb212_fjzO(%esp)
        movapd %xmm0,nb212_fixO(%esp)
        movapd %xmm1,nb212_fiyO(%esp)
        movapd %xmm2,nb212_fizO(%esp)

        ## O-H1 interaction 
        movapd nb212_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212_rsqOH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd  %xmm0,%xmm0
        subpd  nb212_crf(%esp),%xmm4
        mulpd  nb212_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb212_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH1  
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb212_dxOH1(%esp),%xmm0
        mulpd nb212_dyOH1(%esp),%xmm1
        mulpd nb212_dzOH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixO(%esp),%xmm0
        addpd nb212_fiyO(%esp),%xmm1
        addpd nb212_fizO(%esp),%xmm2
        movapd %xmm3,nb212_fjxH1(%esp)
        movapd %xmm4,nb212_fjyH1(%esp)
        movapd %xmm5,nb212_fjzH1(%esp)
        movapd %xmm0,nb212_fixO(%esp)
        movapd %xmm1,nb212_fiyO(%esp)
        movapd %xmm2,nb212_fizO(%esp)

        ## O-H2 interaction  
        movapd nb212_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212_rsqOH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb212_crf(%esp),%xmm4
        mulpd  nb212_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb212_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb212_dxOH2(%esp),%xmm0
        mulpd nb212_dyOH2(%esp),%xmm1
        mulpd nb212_dzOH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixO(%esp),%xmm0
        addpd nb212_fiyO(%esp),%xmm1
        addpd nb212_fizO(%esp),%xmm2
        movapd %xmm3,nb212_fjxH2(%esp)
        movapd %xmm4,nb212_fjyH2(%esp)
        movapd %xmm5,nb212_fjzH2(%esp)
        movapd %xmm0,nb212_fixO(%esp)
        movapd %xmm1,nb212_fiyO(%esp)
        movapd %xmm2,nb212_fizO(%esp)

        ## H1-O interaction 
        movapd nb212_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212_rsqH1O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb212_crf(%esp),%xmm4
        mulpd  nb212_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb212_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxO(%esp),%xmm3
        movapd nb212_fjyO(%esp),%xmm4
        movapd nb212_fjzO(%esp),%xmm5
        mulpd nb212_dxH1O(%esp),%xmm0
        mulpd nb212_dyH1O(%esp),%xmm1
        mulpd nb212_dzH1O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixH1(%esp),%xmm0
        addpd nb212_fiyH1(%esp),%xmm1
        addpd nb212_fizH1(%esp),%xmm2
        movapd %xmm3,nb212_fjxO(%esp)
        movapd %xmm4,nb212_fjyO(%esp)
        movapd %xmm5,nb212_fjzO(%esp)
        movapd %xmm0,nb212_fixH1(%esp)
        movapd %xmm1,nb212_fiyH1(%esp)
        movapd %xmm2,nb212_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb212_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb212_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb212_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb212_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxH1(%esp),%xmm3
        movapd nb212_fjyH1(%esp),%xmm4
        movapd nb212_fjzH1(%esp),%xmm5
        mulpd nb212_dxH1H1(%esp),%xmm0
        mulpd nb212_dyH1H1(%esp),%xmm1
        mulpd nb212_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixH1(%esp),%xmm0
        addpd nb212_fiyH1(%esp),%xmm1
        addpd nb212_fizH1(%esp),%xmm2
        movapd %xmm3,nb212_fjxH1(%esp)
        movapd %xmm4,nb212_fjyH1(%esp)
        movapd %xmm5,nb212_fjzH1(%esp)
        movapd %xmm0,nb212_fixH1(%esp)
        movapd %xmm1,nb212_fiyH1(%esp)
        movapd %xmm2,nb212_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb212_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb212_crf(%esp),%xmm4
        mulpd  nb212_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb212_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxH2(%esp),%xmm3
        movapd nb212_fjyH2(%esp),%xmm4
        movapd nb212_fjzH2(%esp),%xmm5
        mulpd nb212_dxH1H2(%esp),%xmm0
        mulpd nb212_dyH1H2(%esp),%xmm1
        mulpd nb212_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixH1(%esp),%xmm0
        addpd nb212_fiyH1(%esp),%xmm1
        addpd nb212_fizH1(%esp),%xmm2
        movapd %xmm3,nb212_fjxH2(%esp)
        movapd %xmm4,nb212_fjyH2(%esp)
        movapd %xmm5,nb212_fjzH2(%esp)
        movapd %xmm0,nb212_fixH1(%esp)
        movapd %xmm1,nb212_fiyH1(%esp)
        movapd %xmm2,nb212_fizH1(%esp)

        ## H2-O interaction 
        movapd nb212_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212_rsqH2O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb212_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb212_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb212_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxO(%esp),%xmm3
        movapd nb212_fjyO(%esp),%xmm4
        movapd nb212_fjzO(%esp),%xmm5
        mulpd nb212_dxH2O(%esp),%xmm0
        mulpd nb212_dyH2O(%esp),%xmm1
        mulpd nb212_dzH2O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixH2(%esp),%xmm0
        addpd nb212_fiyH2(%esp),%xmm1
        addpd nb212_fizH2(%esp),%xmm2
        movapd %xmm3,nb212_fjxO(%esp)
        movapd %xmm4,nb212_fjyO(%esp)
        movapd %xmm5,nb212_fjzO(%esp)
        movapd %xmm0,nb212_fixH2(%esp)
        movapd %xmm1,nb212_fiyH2(%esp)
        movapd %xmm2,nb212_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb212_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb212_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb212_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb212_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxH1(%esp),%xmm3
        movapd nb212_fjyH1(%esp),%xmm4
        movapd nb212_fjzH1(%esp),%xmm5
        mulpd nb212_dxH2H1(%esp),%xmm0
        mulpd nb212_dyH2H1(%esp),%xmm1
        mulpd nb212_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixH2(%esp),%xmm0
        addpd nb212_fiyH2(%esp),%xmm1
        addpd nb212_fizH2(%esp),%xmm2
        movapd %xmm3,nb212_fjxH1(%esp)
        movapd %xmm4,nb212_fjyH1(%esp)
        movapd %xmm5,nb212_fjzH1(%esp)
        movapd %xmm0,nb212_fixH2(%esp)
        movapd %xmm1,nb212_fiyH2(%esp)
        movapd %xmm2,nb212_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb212_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb212_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb212_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb212_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb212_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd %xmm0,%xmm1
        movapd %xmm6,nb212_vctot(%esp)
        movapd %xmm0,%xmm2

        movapd nb212_fjxH2(%esp),%xmm3
        movapd nb212_fjyH2(%esp),%xmm4
        movapd nb212_fjzH2(%esp),%xmm5
        mulpd nb212_dxH2H2(%esp),%xmm0
        mulpd nb212_dyH2H2(%esp),%xmm1
        mulpd nb212_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb212_fixH2(%esp),%xmm0
        addpd nb212_fiyH2(%esp),%xmm1
        addpd nb212_fizH2(%esp),%xmm2
        movapd %xmm3,nb212_fjxH2(%esp)
        movapd %xmm4,nb212_fjyH2(%esp)
        movapd %xmm5,nb212_fjzH2(%esp)
        movapd %xmm0,nb212_fixH2(%esp)
        movapd %xmm1,nb212_fiyH2(%esp)
        movapd %xmm2,nb212_fizH2(%esp)

        movl nb212_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2
        movlpd 24(%edi,%eax,8),%xmm3
        movlpd 32(%edi,%eax,8),%xmm4
        movlpd 40(%edi,%eax,8),%xmm5
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 56(%edi,%eax,8),%xmm7
        movhpd (%edi,%ebx,8),%xmm0
        movhpd 8(%edi,%ebx,8),%xmm1
        movhpd 16(%edi,%ebx,8),%xmm2
        movhpd 24(%edi,%ebx,8),%xmm3
        movhpd 32(%edi,%ebx,8),%xmm4
        movhpd 40(%edi,%ebx,8),%xmm5
        movhpd 48(%edi,%ebx,8),%xmm6
        movhpd 56(%edi,%ebx,8),%xmm7
        addpd nb212_fjxO(%esp),%xmm0
        addpd nb212_fjyO(%esp),%xmm1
        addpd nb212_fjzO(%esp),%xmm2
        addpd nb212_fjxH1(%esp),%xmm3
        addpd nb212_fjyH1(%esp),%xmm4
        addpd nb212_fjzH1(%esp),%xmm5
        addpd nb212_fjxH2(%esp),%xmm6
        addpd nb212_fjyH2(%esp),%xmm7
        movlpd %xmm0,(%edi,%eax,8)
        movlpd %xmm1,8(%edi,%eax,8)
        movlpd %xmm2,16(%edi,%eax,8)
        movlpd %xmm3,24(%edi,%eax,8)
        movlpd %xmm4,32(%edi,%eax,8)
        movlpd %xmm5,40(%edi,%eax,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,56(%edi,%eax,8)
        movhpd %xmm0,(%edi,%ebx,8)
        movhpd %xmm1,8(%edi,%ebx,8)
        movhpd %xmm2,16(%edi,%ebx,8)
        movhpd %xmm3,24(%edi,%ebx,8)
        movhpd %xmm4,32(%edi,%ebx,8)
        movhpd %xmm5,40(%edi,%ebx,8)
        movhpd %xmm6,48(%edi,%ebx,8)
        movhpd %xmm7,56(%edi,%ebx,8)

        movlpd 64(%edi,%eax,8),%xmm0
        movhpd 64(%edi,%ebx,8),%xmm0
        addpd nb212_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)
        movhpd %xmm0,64(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb212_innerk(%esp)
        jl    _nb_kernel212_ia32_sse2.nb212_checksingle
        jmp   _nb_kernel212_ia32_sse2.nb212_unroll_loop
_nb_kernel212_ia32_sse2.nb212_checksingle: 
        movl  nb212_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel212_ia32_sse2.nb212_dosingle
        jmp   _nb_kernel212_ia32_sse2.nb212_updateouterdata
_nb_kernel212_ia32_sse2.nb212_dosingle: 
        movl  nb212_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb212_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb212_jxO(%esp)
        movapd  %xmm3,nb212_jyO(%esp)
        movapd  %xmm4,nb212_jzO(%esp)
        movapd  %xmm5,nb212_jxH1(%esp)
        movapd  %xmm6,nb212_jyH1(%esp)
        movapd  %xmm7,nb212_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb212_jxH2(%esp)
        movapd  %xmm3,nb212_jyH2(%esp)
        movapd  %xmm4,nb212_jzH2(%esp)

        movapd nb212_ixO(%esp),%xmm0
        movapd nb212_iyO(%esp),%xmm1
        movapd nb212_izO(%esp),%xmm2
        movapd nb212_ixO(%esp),%xmm3
        movapd nb212_iyO(%esp),%xmm4
        movapd nb212_izO(%esp),%xmm5
        subsd  nb212_jxO(%esp),%xmm0
        subsd  nb212_jyO(%esp),%xmm1
        subsd  nb212_jzO(%esp),%xmm2
        subsd  nb212_jxH1(%esp),%xmm3
        subsd  nb212_jyH1(%esp),%xmm4
        subsd  nb212_jzH1(%esp),%xmm5
        movapd %xmm0,nb212_dxOO(%esp)
        movapd %xmm1,nb212_dyOO(%esp)
        movapd %xmm2,nb212_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb212_dxOH1(%esp)
        movapd %xmm4,nb212_dyOH1(%esp)
        movapd %xmm5,nb212_dzOH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb212_rsqOO(%esp)
        movapd %xmm3,nb212_rsqOH1(%esp)

        movapd nb212_ixO(%esp),%xmm0
        movapd nb212_iyO(%esp),%xmm1
        movapd nb212_izO(%esp),%xmm2
        movapd nb212_ixH1(%esp),%xmm3
        movapd nb212_iyH1(%esp),%xmm4
        movapd nb212_izH1(%esp),%xmm5
        subsd  nb212_jxH2(%esp),%xmm0
        subsd  nb212_jyH2(%esp),%xmm1
        subsd  nb212_jzH2(%esp),%xmm2
        subsd  nb212_jxO(%esp),%xmm3
        subsd  nb212_jyO(%esp),%xmm4
        subsd  nb212_jzO(%esp),%xmm5
        movapd %xmm0,nb212_dxOH2(%esp)
        movapd %xmm1,nb212_dyOH2(%esp)
        movapd %xmm2,nb212_dzOH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb212_dxH1O(%esp)
        movapd %xmm4,nb212_dyH1O(%esp)
        movapd %xmm5,nb212_dzH1O(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb212_rsqOH2(%esp)
        movapd %xmm3,nb212_rsqH1O(%esp)

        movapd nb212_ixH1(%esp),%xmm0
        movapd nb212_iyH1(%esp),%xmm1
        movapd nb212_izH1(%esp),%xmm2
        movapd nb212_ixH1(%esp),%xmm3
        movapd nb212_iyH1(%esp),%xmm4
        movapd nb212_izH1(%esp),%xmm5
        subsd  nb212_jxH1(%esp),%xmm0
        subsd  nb212_jyH1(%esp),%xmm1
        subsd  nb212_jzH1(%esp),%xmm2
        subsd  nb212_jxH2(%esp),%xmm3
        subsd  nb212_jyH2(%esp),%xmm4
        subsd  nb212_jzH2(%esp),%xmm5
        movapd %xmm0,nb212_dxH1H1(%esp)
        movapd %xmm1,nb212_dyH1H1(%esp)
        movapd %xmm2,nb212_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb212_dxH1H2(%esp)
        movapd %xmm4,nb212_dyH1H2(%esp)
        movapd %xmm5,nb212_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb212_rsqH1H1(%esp)
        movapd %xmm3,nb212_rsqH1H2(%esp)

        movapd nb212_ixH2(%esp),%xmm0
        movapd nb212_iyH2(%esp),%xmm1
        movapd nb212_izH2(%esp),%xmm2
        movapd nb212_ixH2(%esp),%xmm3
        movapd nb212_iyH2(%esp),%xmm4
        movapd nb212_izH2(%esp),%xmm5
        subsd  nb212_jxO(%esp),%xmm0
        subsd  nb212_jyO(%esp),%xmm1
        subsd  nb212_jzO(%esp),%xmm2
        subsd  nb212_jxH1(%esp),%xmm3
        subsd  nb212_jyH1(%esp),%xmm4
        subsd  nb212_jzH1(%esp),%xmm5
        movapd %xmm0,nb212_dxH2O(%esp)
        movapd %xmm1,nb212_dyH2O(%esp)
        movapd %xmm2,nb212_dzH2O(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb212_dxH2H1(%esp)
        movapd %xmm4,nb212_dyH2H1(%esp)
        movapd %xmm5,nb212_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb212_rsqH2O(%esp)
        movapd %xmm4,nb212_rsqH2H1(%esp)

        movapd nb212_ixH2(%esp),%xmm0
        movapd nb212_iyH2(%esp),%xmm1
        movapd nb212_izH2(%esp),%xmm2
        subsd  nb212_jxH2(%esp),%xmm0
        subsd  nb212_jyH2(%esp),%xmm1
        subsd  nb212_jzH2(%esp),%xmm2
        movapd %xmm0,nb212_dxH2H2(%esp)
        movapd %xmm1,nb212_dyH2H2(%esp)
        movapd %xmm2,nb212_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb212_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
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
        movapd  nb212_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212_half(%esp),%xmm3   ## iter1 
        mulsd   nb212_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212_half(%esp),%xmm1   ## rinv 
        mulsd   nb212_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212_rinvH2H2(%esp)
        movapd %xmm5,nb212_rinvH2H1(%esp)

        movapd nb212_rsqOO(%esp),%xmm0
        movapd nb212_rsqOH1(%esp),%xmm4
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
        movapd  nb212_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb212_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212_half(%esp),%xmm1   ## rinv 
        mulsd   nb212_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb212_rinvOO(%esp)
        movapd %xmm5,nb212_rinvOH1(%esp)

        movapd nb212_rsqOH2(%esp),%xmm0
        movapd nb212_rsqH1O(%esp),%xmm4
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
        movapd  nb212_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212_half(%esp),%xmm3   ## iter1 
        mulsd   nb212_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212_half(%esp),%xmm1   ## rinv 
        mulsd   nb212_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212_rinvOH2(%esp)
        movapd %xmm5,nb212_rinvH1O(%esp)

        movapd nb212_rsqH1H1(%esp),%xmm0
        movapd nb212_rsqH1H2(%esp),%xmm4
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
        movapd  nb212_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212_half(%esp),%xmm3   ## iter1a 
        mulsd   nb212_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212_half(%esp),%xmm1   ## rinv 
        mulsd   nb212_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212_rinvH1H1(%esp)
        movapd %xmm5,nb212_rinvH1H2(%esp)

        movapd nb212_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb212_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb212_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb212_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb212_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb212_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb212_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulsd  nb212_rsqOO(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb212_crf(%esp),%xmm6

        mulsd  nb212_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqOO(%esp),%xmm7   ## xmm7 = coul part of fscal 

        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb212_c6(%esp),%xmm1
        mulsd  nb212_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addsd  nb212_Vvdwtot(%esp),%xmm3
        mulsd  nb212_six(%esp),%xmm1
        mulsd  nb212_twelve(%esp),%xmm2
        movlpd %xmm3,nb212_Vvdwtot(%esp)
        subsd  %xmm1,%xmm2
        addsd  %xmm7,%xmm2
        addsd  nb212_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm2,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb212_dxOO(%esp),%xmm0
        mulsd nb212_dyOO(%esp),%xmm1
        mulsd nb212_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixO(%esp),%xmm0
        addsd nb212_fiyO(%esp),%xmm1
        addsd nb212_fizO(%esp),%xmm2
        movlpd %xmm3,nb212_fjxO(%esp)
        movlpd %xmm4,nb212_fjyO(%esp)
        movlpd %xmm5,nb212_fjzO(%esp)
        movlpd %xmm0,nb212_fixO(%esp)
        movlpd %xmm1,nb212_fiyO(%esp)
        movlpd %xmm2,nb212_fizO(%esp)

        ## O-H1 interaction 
        movapd nb212_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212_rsqOH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb212_crf(%esp),%xmm4
        mulsd  nb212_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH1  
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb212_dxOH1(%esp),%xmm0
        mulsd nb212_dyOH1(%esp),%xmm1
        mulsd nb212_dzOH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixO(%esp),%xmm0
        addsd nb212_fiyO(%esp),%xmm1
        addsd nb212_fizO(%esp),%xmm2
        movlpd %xmm3,nb212_fjxH1(%esp)
        movlpd %xmm4,nb212_fjyH1(%esp)
        movlpd %xmm5,nb212_fjzH1(%esp)
        movlpd %xmm0,nb212_fixO(%esp)
        movlpd %xmm1,nb212_fiyO(%esp)
        movlpd %xmm2,nb212_fizO(%esp)

        ## O-H2 interaction  
        movapd nb212_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212_rsqOH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb212_crf(%esp),%xmm4
        mulsd  nb212_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb212_dxOH2(%esp),%xmm0
        mulsd nb212_dyOH2(%esp),%xmm1
        mulsd nb212_dzOH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixO(%esp),%xmm0
        addsd nb212_fiyO(%esp),%xmm1
        addsd nb212_fizO(%esp),%xmm2
        movlpd %xmm3,nb212_fjxH2(%esp)
        movlpd %xmm4,nb212_fjyH2(%esp)
        movlpd %xmm5,nb212_fjzH2(%esp)
        movlpd %xmm0,nb212_fixO(%esp)
        movlpd %xmm1,nb212_fiyO(%esp)
        movlpd %xmm2,nb212_fizO(%esp)

        ## H1-O interaction 
        movapd nb212_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212_rsqH1O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb212_crf(%esp),%xmm4
        mulsd  nb212_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxO(%esp),%xmm3
        movapd nb212_fjyO(%esp),%xmm4
        movapd nb212_fjzO(%esp),%xmm5
        mulsd nb212_dxH1O(%esp),%xmm0
        mulsd nb212_dyH1O(%esp),%xmm1
        mulsd nb212_dzH1O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixH1(%esp),%xmm0
        addsd nb212_fiyH1(%esp),%xmm1
        addsd nb212_fizH1(%esp),%xmm2
        movlpd %xmm3,nb212_fjxO(%esp)
        movlpd %xmm4,nb212_fjyO(%esp)
        movlpd %xmm5,nb212_fjzO(%esp)
        movlpd %xmm0,nb212_fixH1(%esp)
        movlpd %xmm1,nb212_fiyH1(%esp)
        movlpd %xmm2,nb212_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb212_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb212_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb212_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxH1(%esp),%xmm3
        movapd nb212_fjyH1(%esp),%xmm4
        movapd nb212_fjzH1(%esp),%xmm5
        mulsd nb212_dxH1H1(%esp),%xmm0
        mulsd nb212_dyH1H1(%esp),%xmm1
        mulsd nb212_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixH1(%esp),%xmm0
        addsd nb212_fiyH1(%esp),%xmm1
        addsd nb212_fizH1(%esp),%xmm2
        movlpd %xmm3,nb212_fjxH1(%esp)
        movlpd %xmm4,nb212_fjyH1(%esp)
        movlpd %xmm5,nb212_fjzH1(%esp)
        movlpd %xmm0,nb212_fixH1(%esp)
        movlpd %xmm1,nb212_fiyH1(%esp)
        movlpd %xmm2,nb212_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb212_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb212_crf(%esp),%xmm4
        mulsd  nb212_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxH2(%esp),%xmm3
        movapd nb212_fjyH2(%esp),%xmm4
        movapd nb212_fjzH2(%esp),%xmm5
        mulsd nb212_dxH1H2(%esp),%xmm0
        mulsd nb212_dyH1H2(%esp),%xmm1
        mulsd nb212_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixH1(%esp),%xmm0
        addsd nb212_fiyH1(%esp),%xmm1
        addsd nb212_fizH1(%esp),%xmm2
        movlpd %xmm3,nb212_fjxH2(%esp)
        movlpd %xmm4,nb212_fjyH2(%esp)
        movlpd %xmm5,nb212_fjzH2(%esp)
        movlpd %xmm0,nb212_fixH1(%esp)
        movlpd %xmm1,nb212_fiyH1(%esp)
        movlpd %xmm2,nb212_fizH1(%esp)

        ## H2-O interaction 
        movapd nb212_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212_rsqH2O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb212_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb212_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxO(%esp),%xmm3
        movapd nb212_fjyO(%esp),%xmm4
        movapd nb212_fjzO(%esp),%xmm5
        mulsd nb212_dxH2O(%esp),%xmm0
        mulsd nb212_dyH2O(%esp),%xmm1
        mulsd nb212_dzH2O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixH2(%esp),%xmm0
        addsd nb212_fiyH2(%esp),%xmm1
        addsd nb212_fizH2(%esp),%xmm2
        movlpd %xmm3,nb212_fjxO(%esp)
        movlpd %xmm4,nb212_fjyO(%esp)
        movlpd %xmm5,nb212_fjzO(%esp)
        movlpd %xmm0,nb212_fixH2(%esp)
        movlpd %xmm1,nb212_fiyH2(%esp)
        movlpd %xmm2,nb212_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb212_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb212_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb212_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb212_fjxH1(%esp),%xmm3
        movapd nb212_fjyH1(%esp),%xmm4
        movapd nb212_fjzH1(%esp),%xmm5
        mulsd nb212_dxH2H1(%esp),%xmm0
        mulsd nb212_dyH2H1(%esp),%xmm1
        mulsd nb212_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixH2(%esp),%xmm0
        addsd nb212_fiyH2(%esp),%xmm1
        addsd nb212_fizH2(%esp),%xmm2
        movlpd %xmm3,nb212_fjxH1(%esp)
        movlpd %xmm4,nb212_fjyH1(%esp)
        movlpd %xmm5,nb212_fjzH1(%esp)
        movlpd %xmm0,nb212_fixH2(%esp)
        movlpd %xmm1,nb212_fiyH2(%esp)
        movlpd %xmm2,nb212_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb212_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb212_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb212_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb212_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb212_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd %xmm0,%xmm1
        movlpd %xmm6,nb212_vctot(%esp)
        movapd %xmm0,%xmm2

        movapd nb212_fjxH2(%esp),%xmm3
        movapd nb212_fjyH2(%esp),%xmm4
        movapd nb212_fjzH2(%esp),%xmm5
        mulsd nb212_dxH2H2(%esp),%xmm0
        mulsd nb212_dyH2H2(%esp),%xmm1
        mulsd nb212_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb212_fixH2(%esp),%xmm0
        addsd nb212_fiyH2(%esp),%xmm1
        addsd nb212_fizH2(%esp),%xmm2
        movlpd %xmm3,nb212_fjxH2(%esp)
        movlpd %xmm4,nb212_fjyH2(%esp)
        movlpd %xmm5,nb212_fjzH2(%esp)
        movlpd %xmm0,nb212_fixH2(%esp)
        movlpd %xmm1,nb212_fiyH2(%esp)
        movlpd %xmm2,nb212_fizH2(%esp)

        movl nb212_faction(%ebp),%edi
        ## Did all interactions - now update j forces 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2
        movlpd 24(%edi,%eax,8),%xmm3
        movlpd 32(%edi,%eax,8),%xmm4
        movlpd 40(%edi,%eax,8),%xmm5
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 56(%edi,%eax,8),%xmm7
        addsd nb212_fjxO(%esp),%xmm0
        addsd nb212_fjyO(%esp),%xmm1
        addsd nb212_fjzO(%esp),%xmm2
        addsd nb212_fjxH1(%esp),%xmm3
        addsd nb212_fjyH1(%esp),%xmm4
        addsd nb212_fjzH1(%esp),%xmm5
        addsd nb212_fjxH2(%esp),%xmm6
        addsd nb212_fjyH2(%esp),%xmm7
        movlpd %xmm0,(%edi,%eax,8)
        movlpd %xmm1,8(%edi,%eax,8)
        movlpd %xmm2,16(%edi,%eax,8)
        movlpd %xmm3,24(%edi,%eax,8)
        movlpd %xmm4,32(%edi,%eax,8)
        movlpd %xmm5,40(%edi,%eax,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,56(%edi,%eax,8)

        movlpd 64(%edi,%eax,8),%xmm0
        addsd nb212_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)

_nb_kernel212_ia32_sse2.nb212_updateouterdata: 
        movl  nb212_ii3(%esp),%ecx
        movl  nb212_faction(%ebp),%edi
        movl  nb212_fshift(%ebp),%esi
        movl  nb212_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb212_fixO(%esp),%xmm0
        movapd nb212_fiyO(%esp),%xmm1
        movapd nb212_fizO(%esp),%xmm2

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
        movapd nb212_fixH1(%esp),%xmm0
        movapd nb212_fiyH1(%esp),%xmm1
        movapd nb212_fizH1(%esp),%xmm2

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
        movapd nb212_fixH2(%esp),%xmm0
        movapd nb212_fiyH2(%esp),%xmm1
        movapd nb212_fizH2(%esp),%xmm2

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
        movl nb212_n(%esp),%esi
        ## get group index for i particle 
        movl  nb212_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb212_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb212_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb212_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb212_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb212_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel212_ia32_sse2.nb212_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb212_n(%esp)
        jmp _nb_kernel212_ia32_sse2.nb212_outer
_nb_kernel212_ia32_sse2.nb212_outerend: 
        ## check if more outer neighborlists remain
        movl  nb212_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel212_ia32_sse2.nb212_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel212_ia32_sse2.nb212_threadloop
_nb_kernel212_ia32_sse2.nb212_end: 
        emms

        movl nb212_nouter(%esp),%eax
        movl nb212_ninner(%esp),%ebx
        movl nb212_outeriter(%ebp),%ecx
        movl nb212_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb212_salign(%esp),%eax
        addl %eax,%esp
        addl $1560,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret




.globl nb_kernel212nf_ia32_sse2
.globl _nb_kernel212nf_ia32_sse2
nb_kernel212nf_ia32_sse2:       
_nb_kernel212nf_ia32_sse2:      
.set nb212nf_p_nri, 8
.set nb212nf_iinr, 12
.set nb212nf_jindex, 16
.set nb212nf_jjnr, 20
.set nb212nf_shift, 24
.set nb212nf_shiftvec, 28
.set nb212nf_fshift, 32
.set nb212nf_gid, 36
.set nb212nf_pos, 40
.set nb212nf_faction, 44
.set nb212nf_charge, 48
.set nb212nf_p_facel, 52
.set nb212nf_argkrf, 56
.set nb212nf_argcrf, 60
.set nb212nf_Vc, 64
.set nb212nf_type, 68
.set nb212nf_p_ntype, 72
.set nb212nf_vdwparam, 76
.set nb212nf_Vvdw, 80
.set nb212nf_p_tabscale, 84
.set nb212nf_VFtab, 88
.set nb212nf_invsqrta, 92
.set nb212nf_dvda, 96
.set nb212nf_p_gbtabscale, 100
.set nb212nf_GBtab, 104
.set nb212nf_p_nthreads, 108
.set nb212nf_count, 112
.set nb212nf_mtx, 116
.set nb212nf_outeriter, 120
.set nb212nf_inneriter, 124
.set nb212nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb212nf_ixO, 0
.set nb212nf_iyO, 16
.set nb212nf_izO, 32
.set nb212nf_ixH1, 48
.set nb212nf_iyH1, 64
.set nb212nf_izH1, 80
.set nb212nf_ixH2, 96
.set nb212nf_iyH2, 112
.set nb212nf_izH2, 128
.set nb212nf_jxO, 144
.set nb212nf_jyO, 160
.set nb212nf_jzO, 176
.set nb212nf_jxH1, 192
.set nb212nf_jyH1, 208
.set nb212nf_jzH1, 224
.set nb212nf_jxH2, 240
.set nb212nf_jyH2, 256
.set nb212nf_jzH2, 272
.set nb212nf_qqOO, 288
.set nb212nf_qqOH, 304
.set nb212nf_qqHH, 320
.set nb212nf_c6, 336
.set nb212nf_c12, 352
.set nb212nf_vctot, 368
.set nb212nf_Vvdwtot, 384
.set nb212nf_half, 400
.set nb212nf_three, 416
.set nb212nf_rsqOO, 432
.set nb212nf_rsqOH1, 448
.set nb212nf_rsqOH2, 464
.set nb212nf_rsqH1O, 480
.set nb212nf_rsqH1H1, 496
.set nb212nf_rsqH1H2, 512
.set nb212nf_rsqH2O, 528
.set nb212nf_rsqH2H1, 544
.set nb212nf_rsqH2H2, 560
.set nb212nf_rinvOO, 576
.set nb212nf_rinvOH1, 592
.set nb212nf_rinvOH2, 608
.set nb212nf_rinvH1O, 624
.set nb212nf_rinvH1H1, 640
.set nb212nf_rinvH1H2, 656
.set nb212nf_rinvH2O, 672
.set nb212nf_rinvH2H1, 688
.set nb212nf_rinvH2H2, 704
.set nb212nf_krf, 720
.set nb212nf_crf, 736
.set nb212nf_is3, 752
.set nb212nf_ii3, 756
.set nb212nf_innerjjnr, 760
.set nb212nf_innerk, 764
.set nb212nf_n, 768
.set nb212nf_nn1, 772
.set nb212nf_nri, 776
.set nb212nf_nouter, 780
.set nb212nf_ninner, 784
.set nb212nf_salign, 788
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $792,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb212nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb212nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb212nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb212nf_nouter(%esp)
        movl %eax,nb212nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb212nf_half(%esp)
        movl %ebx,nb212nf_half+4(%esp)
        movsd nb212nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb212nf_half(%esp)
        movapd %xmm3,nb212nf_three(%esp)
        movl nb212nf_argkrf(%ebp),%esi
        movl nb212nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb212nf_krf(%esp)
        movapd %xmm6,nb212nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb212nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb212nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb212nf_p_facel(%ebp),%esi
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
        movapd %xmm3,nb212nf_qqOO(%esp)
        movapd %xmm4,nb212nf_qqOH(%esp)
        movapd %xmm5,nb212nf_qqHH(%esp)

        xorpd %xmm0,%xmm0
        movl  nb212nf_type(%ebp),%edx
        movl  (%edx,%ebx,4),%ecx
        shll  %ecx
        movl  %ecx,%edx
        movl nb212nf_p_ntype(%ebp),%edi
        imull (%edi),%ecx     ## ecx = ntia = 2*ntype*type[ii0] 
        addl  %ecx,%edx
        movl  nb212nf_vdwparam(%ebp),%eax
        movlpd (%eax,%edx,8),%xmm0
        movhpd 8(%eax,%edx,8),%xmm0
        movhlps %xmm0,%xmm1
        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        movapd %xmm0,nb212nf_c6(%esp)
        movapd %xmm1,nb212nf_c12(%esp)

_nb_kernel212nf_ia32_sse2.nb212nf_threadloop: 
        movl  nb212nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel212nf_ia32_sse2.nb212nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel212nf_ia32_sse2.nb212nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb212nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb212nf_n(%esp)
        movl %ebx,nb212nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel212nf_ia32_sse2.nb212nf_outerstart
        jmp _nb_kernel212nf_ia32_sse2.nb212nf_end

_nb_kernel212nf_ia32_sse2.nb212nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb212nf_nouter(%esp),%ebx
        movl %ebx,nb212nf_nouter(%esp)

_nb_kernel212nf_ia32_sse2.nb212nf_outer: 
        movl  nb212nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb212nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movlpd (%eax,%ebx,8),%xmm0
        movlpd 8(%eax,%ebx,8),%xmm1
        movlpd 16(%eax,%ebx,8),%xmm2

        movl  nb212nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb212nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb212nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb212nf_ixO(%esp)
        movapd %xmm4,nb212nf_iyO(%esp)
        movapd %xmm5,nb212nf_izO(%esp)

        movsd %xmm0,%xmm3
        movsd %xmm1,%xmm4
        movsd %xmm2,%xmm5
        addsd 24(%eax,%ebx,8),%xmm0
        addsd 32(%eax,%ebx,8),%xmm1
        addsd 40(%eax,%ebx,8),%xmm2
        addsd 48(%eax,%ebx,8),%xmm3
        addsd 56(%eax,%ebx,8),%xmm4
        addsd 64(%eax,%ebx,8),%xmm5

        shufpd $0,%xmm0,%xmm0
        shufpd $0,%xmm1,%xmm1
        shufpd $0,%xmm2,%xmm2
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm0,nb212nf_ixH1(%esp)
        movapd %xmm1,nb212nf_iyH1(%esp)
        movapd %xmm2,nb212nf_izH1(%esp)
        movapd %xmm3,nb212nf_ixH2(%esp)
        movapd %xmm4,nb212nf_iyH2(%esp)
        movapd %xmm5,nb212nf_izH2(%esp)

        ## clear vctot & Vvdwtot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb212nf_vctot(%esp)
        movapd %xmm4,nb212nf_Vvdwtot(%esp)

        movl  nb212nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb212nf_pos(%ebp),%esi
        movl  nb212nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb212nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb212nf_ninner(%esp),%ecx
        movl  %ecx,nb212nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb212nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel212nf_ia32_sse2.nb212nf_unroll_loop
        jmp   _nb_kernel212nf_ia32_sse2.nb212nf_checksingle
_nb_kernel212nf_ia32_sse2.nb212nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb212nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb212nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb212nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx

        ## move j coordinates to local temp variables 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movhpd (%esi,%ebx,8),%xmm2
        movhpd 8(%esi,%ebx,8),%xmm3
        movhpd 16(%esi,%ebx,8),%xmm4
        movhpd 24(%esi,%ebx,8),%xmm5
        movhpd 32(%esi,%ebx,8),%xmm6
        movhpd 40(%esi,%ebx,8),%xmm7
        movapd  %xmm2,nb212nf_jxO(%esp)
        movapd  %xmm3,nb212nf_jyO(%esp)
        movapd  %xmm4,nb212nf_jzO(%esp)
        movapd  %xmm5,nb212nf_jxH1(%esp)
        movapd  %xmm6,nb212nf_jyH1(%esp)
        movapd  %xmm7,nb212nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb212nf_jxH2(%esp)
        movapd  %xmm3,nb212nf_jyH2(%esp)
        movapd  %xmm4,nb212nf_jzH2(%esp)

        movapd nb212nf_ixO(%esp),%xmm0
        movapd nb212nf_iyO(%esp),%xmm1
        movapd nb212nf_izO(%esp),%xmm2
        movapd nb212nf_ixO(%esp),%xmm3
        movapd nb212nf_iyO(%esp),%xmm4
        movapd nb212nf_izO(%esp),%xmm5
        subpd  nb212nf_jxO(%esp),%xmm0
        subpd  nb212nf_jyO(%esp),%xmm1
        subpd  nb212nf_jzO(%esp),%xmm2
        subpd  nb212nf_jxH1(%esp),%xmm3
        subpd  nb212nf_jyH1(%esp),%xmm4
        subpd  nb212nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb212nf_rsqOO(%esp)
        movapd %xmm3,nb212nf_rsqOH1(%esp)

        movapd nb212nf_ixO(%esp),%xmm0
        movapd nb212nf_iyO(%esp),%xmm1
        movapd nb212nf_izO(%esp),%xmm2
        movapd nb212nf_ixH1(%esp),%xmm3
        movapd nb212nf_iyH1(%esp),%xmm4
        movapd nb212nf_izH1(%esp),%xmm5
        subpd  nb212nf_jxH2(%esp),%xmm0
        subpd  nb212nf_jyH2(%esp),%xmm1
        subpd  nb212nf_jzH2(%esp),%xmm2
        subpd  nb212nf_jxO(%esp),%xmm3
        subpd  nb212nf_jyO(%esp),%xmm4
        subpd  nb212nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb212nf_rsqOH2(%esp)
        movapd %xmm3,nb212nf_rsqH1O(%esp)

        movapd nb212nf_ixH1(%esp),%xmm0
        movapd nb212nf_iyH1(%esp),%xmm1
        movapd nb212nf_izH1(%esp),%xmm2
        movapd nb212nf_ixH1(%esp),%xmm3
        movapd nb212nf_iyH1(%esp),%xmm4
        movapd nb212nf_izH1(%esp),%xmm5
        subpd  nb212nf_jxH1(%esp),%xmm0
        subpd  nb212nf_jyH1(%esp),%xmm1
        subpd  nb212nf_jzH1(%esp),%xmm2
        subpd  nb212nf_jxH2(%esp),%xmm3
        subpd  nb212nf_jyH2(%esp),%xmm4
        subpd  nb212nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb212nf_rsqH1H1(%esp)
        movapd %xmm3,nb212nf_rsqH1H2(%esp)

        movapd nb212nf_ixH2(%esp),%xmm0
        movapd nb212nf_iyH2(%esp),%xmm1
        movapd nb212nf_izH2(%esp),%xmm2
        movapd nb212nf_ixH2(%esp),%xmm3
        movapd nb212nf_iyH2(%esp),%xmm4
        movapd nb212nf_izH2(%esp),%xmm5
        subpd  nb212nf_jxO(%esp),%xmm0
        subpd  nb212nf_jyO(%esp),%xmm1
        subpd  nb212nf_jzO(%esp),%xmm2
        subpd  nb212nf_jxH1(%esp),%xmm3
        subpd  nb212nf_jyH1(%esp),%xmm4
        subpd  nb212nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb212nf_rsqH2O(%esp)
        movapd %xmm4,nb212nf_rsqH2H1(%esp)

        movapd nb212nf_ixH2(%esp),%xmm0
        movapd nb212nf_iyH2(%esp),%xmm1
        movapd nb212nf_izH2(%esp),%xmm2
        subpd  nb212nf_jxH2(%esp),%xmm0
        subpd  nb212nf_jyH2(%esp),%xmm1
        subpd  nb212nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb212nf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
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
        movapd  nb212nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb212nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb212nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212nf_rinvH2H2(%esp)
        movapd %xmm5,nb212nf_rinvH2H1(%esp)

        movapd nb212nf_rsqOO(%esp),%xmm0
        movapd nb212nf_rsqOH1(%esp),%xmm4
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
        movapd  nb212nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb212nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb212nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb212nf_rinvOO(%esp)
        movapd %xmm5,nb212nf_rinvOH1(%esp)

        movapd nb212nf_rsqOH2(%esp),%xmm0
        movapd nb212nf_rsqH1O(%esp),%xmm4
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
        movapd  nb212nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb212nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb212nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212nf_rinvOH2(%esp)
        movapd %xmm5,nb212nf_rinvH1O(%esp)

        movapd nb212nf_rsqH1H1(%esp),%xmm0
        movapd nb212nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb212nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb212nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb212nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb212nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212nf_rinvH1H1(%esp)
        movapd %xmm5,nb212nf_rinvH1H2(%esp)

        movapd nb212nf_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb212nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb212nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb212nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb212nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb212nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb212nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7              ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0              ## xmm0=rinvsq 
        movapd %xmm0,%xmm1
        mulpd  %xmm0,%xmm1              ## rinvsq*rinvsq 
        mulpd  %xmm0,%xmm1              ## xmm1=rinvsix 
        mulpd  nb212nf_rsqOO(%esp),%xmm5        ## xmm5=krsq 
        movapd %xmm5,%xmm6              ## krsq 
        addpd  %xmm7,%xmm6              ## xmm6=rinv+ krsq 
        subpd  nb212nf_crf(%esp),%xmm6          ## rinv+krsq-crf 

        mulpd  nb212nf_qqOO(%esp),%xmm6         ## xmm6=voul=qq*(rinv+ krsq-crf) 

        movapd %xmm1,%xmm2              ## rinv6 
        mulpd  %xmm2,%xmm2              ## xmm2=rinvtwelve 
        mulpd  nb212nf_c6(%esp),%xmm1   ## c6*rinv6 
        mulpd  nb212nf_c12(%esp),%xmm2          ## c12*rinv12 
        movapd %xmm2,%xmm3              ## c12*rinv12 
        subpd  %xmm1,%xmm3              ## Vvdw12-Vvdw6 
        addpd  nb212nf_Vvdwtot(%esp),%xmm3
        movapd %xmm3,nb212nf_Vvdwtot(%esp)
        addpd  nb212nf_vctot(%esp),%xmm6   ## local vctot summation variable 

        ## O-H1 interaction 
        movapd nb212nf_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212nf_rsqOH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd  %xmm0,%xmm0
        subpd  nb212nf_crf(%esp),%xmm4
        mulpd  nb212nf_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## O-H2 interaction  
        movapd nb212nf_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212nf_rsqOH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb212nf_crf(%esp),%xmm4
        mulpd  nb212nf_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-O interaction 
        movapd nb212nf_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212nf_rsqH1O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb212nf_crf(%esp),%xmm4
        mulpd  nb212nf_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-H1 interaction 
        movapd nb212nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212nf_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb212nf_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb212nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-H2 interaction 
        movapd nb212nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212nf_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb212nf_crf(%esp),%xmm4
        mulpd  nb212nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-O interaction 
        movapd nb212nf_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212nf_rsqH2O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb212nf_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb212nf_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H1 interaction 
        movapd nb212nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212nf_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb212nf_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb212nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H2 interaction 
        movapd nb212nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb212nf_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb212nf_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb212nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        movapd %xmm6,nb212nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb212nf_innerk(%esp)
        jl    _nb_kernel212nf_ia32_sse2.nb212nf_checksingle
        jmp   _nb_kernel212nf_ia32_sse2.nb212nf_unroll_loop
_nb_kernel212nf_ia32_sse2.nb212nf_checksingle: 
        movl  nb212nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel212nf_ia32_sse2.nb212nf_dosingle
        jmp   _nb_kernel212nf_ia32_sse2.nb212nf_updateouterdata
_nb_kernel212nf_ia32_sse2.nb212nf_dosingle: 
        movl  nb212nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb212nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb212nf_jxO(%esp)
        movapd  %xmm3,nb212nf_jyO(%esp)
        movapd  %xmm4,nb212nf_jzO(%esp)
        movapd  %xmm5,nb212nf_jxH1(%esp)
        movapd  %xmm6,nb212nf_jyH1(%esp)
        movapd  %xmm7,nb212nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb212nf_jxH2(%esp)
        movapd  %xmm3,nb212nf_jyH2(%esp)
        movapd  %xmm4,nb212nf_jzH2(%esp)

        movapd nb212nf_ixO(%esp),%xmm0
        movapd nb212nf_iyO(%esp),%xmm1
        movapd nb212nf_izO(%esp),%xmm2
        movapd nb212nf_ixO(%esp),%xmm3
        movapd nb212nf_iyO(%esp),%xmm4
        movapd nb212nf_izO(%esp),%xmm5
        subsd  nb212nf_jxO(%esp),%xmm0
        subsd  nb212nf_jyO(%esp),%xmm1
        subsd  nb212nf_jzO(%esp),%xmm2
        subsd  nb212nf_jxH1(%esp),%xmm3
        subsd  nb212nf_jyH1(%esp),%xmm4
        subsd  nb212nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb212nf_rsqOO(%esp)
        movapd %xmm3,nb212nf_rsqOH1(%esp)

        movapd nb212nf_ixO(%esp),%xmm0
        movapd nb212nf_iyO(%esp),%xmm1
        movapd nb212nf_izO(%esp),%xmm2
        movapd nb212nf_ixH1(%esp),%xmm3
        movapd nb212nf_iyH1(%esp),%xmm4
        movapd nb212nf_izH1(%esp),%xmm5
        subsd  nb212nf_jxH2(%esp),%xmm0
        subsd  nb212nf_jyH2(%esp),%xmm1
        subsd  nb212nf_jzH2(%esp),%xmm2
        subsd  nb212nf_jxO(%esp),%xmm3
        subsd  nb212nf_jyO(%esp),%xmm4
        subsd  nb212nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb212nf_rsqOH2(%esp)
        movapd %xmm3,nb212nf_rsqH1O(%esp)

        movapd nb212nf_ixH1(%esp),%xmm0
        movapd nb212nf_iyH1(%esp),%xmm1
        movapd nb212nf_izH1(%esp),%xmm2
        movapd nb212nf_ixH1(%esp),%xmm3
        movapd nb212nf_iyH1(%esp),%xmm4
        movapd nb212nf_izH1(%esp),%xmm5
        subsd  nb212nf_jxH1(%esp),%xmm0
        subsd  nb212nf_jyH1(%esp),%xmm1
        subsd  nb212nf_jzH1(%esp),%xmm2
        subsd  nb212nf_jxH2(%esp),%xmm3
        subsd  nb212nf_jyH2(%esp),%xmm4
        subsd  nb212nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb212nf_rsqH1H1(%esp)
        movapd %xmm3,nb212nf_rsqH1H2(%esp)

        movapd nb212nf_ixH2(%esp),%xmm0
        movapd nb212nf_iyH2(%esp),%xmm1
        movapd nb212nf_izH2(%esp),%xmm2
        movapd nb212nf_ixH2(%esp),%xmm3
        movapd nb212nf_iyH2(%esp),%xmm4
        movapd nb212nf_izH2(%esp),%xmm5
        subsd  nb212nf_jxO(%esp),%xmm0
        subsd  nb212nf_jyO(%esp),%xmm1
        subsd  nb212nf_jzO(%esp),%xmm2
        subsd  nb212nf_jxH1(%esp),%xmm3
        subsd  nb212nf_jyH1(%esp),%xmm4
        subsd  nb212nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb212nf_rsqH2O(%esp)
        movapd %xmm4,nb212nf_rsqH2H1(%esp)

        movapd nb212nf_ixH2(%esp),%xmm0
        movapd nb212nf_iyH2(%esp),%xmm1
        movapd nb212nf_izH2(%esp),%xmm2
        subsd  nb212nf_jxH2(%esp),%xmm0
        subsd  nb212nf_jyH2(%esp),%xmm1
        subsd  nb212nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb212nf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
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
        movapd  nb212nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb212nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb212nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212nf_rinvH2H2(%esp)
        movapd %xmm5,nb212nf_rinvH2H1(%esp)

        movapd nb212nf_rsqOO(%esp),%xmm0
        movapd nb212nf_rsqOH1(%esp),%xmm4
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
        movapd  nb212nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb212nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb212nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb212nf_rinvOO(%esp)
        movapd %xmm5,nb212nf_rinvOH1(%esp)

        movapd nb212nf_rsqOH2(%esp),%xmm0
        movapd nb212nf_rsqH1O(%esp),%xmm4
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
        movapd  nb212nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb212nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb212nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212nf_rinvOH2(%esp)
        movapd %xmm5,nb212nf_rinvH1O(%esp)

        movapd nb212nf_rsqH1H1(%esp),%xmm0
        movapd nb212nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb212nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb212nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb212nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb212nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb212nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb212nf_rinvH1H1(%esp)
        movapd %xmm5,nb212nf_rinvH1H2(%esp)

        movapd nb212nf_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb212nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb212nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb212nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb212nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb212nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb212nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulsd  nb212nf_rsqOO(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb212nf_crf(%esp),%xmm6

        mulsd  nb212nf_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 

        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=rinvtwelve 
        mulsd  nb212nf_c6(%esp),%xmm1
        mulsd  nb212nf_c12(%esp),%xmm2
        movapd %xmm2,%xmm3
        subsd  %xmm1,%xmm3      ## xmm3=Vvdw12-Vvdw6 
        addsd  nb212nf_Vvdwtot(%esp),%xmm3
        movlpd %xmm3,nb212nf_Vvdwtot(%esp)
        addsd  nb212nf_vctot(%esp),%xmm6   ## local vctot summation variable 

        ## O-H1 interaction 
        movapd nb212nf_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212nf_rsqOH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb212nf_crf(%esp),%xmm4
        mulsd  nb212nf_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## O-H2 interaction  
        movapd nb212nf_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212nf_rsqOH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb212nf_crf(%esp),%xmm4
        mulsd  nb212nf_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-O interaction 
        movapd nb212nf_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212nf_rsqH1O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb212nf_crf(%esp),%xmm4
        mulsd  nb212nf_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-H1 interaction 
        movapd nb212nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212nf_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb212nf_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb212nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H1-H2 interaction 
        movapd nb212nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212nf_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb212nf_crf(%esp),%xmm4
        mulsd  nb212nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-O interaction 
        movapd nb212nf_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212nf_rsqH2O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb212nf_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb212nf_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H1 interaction 
        movapd nb212nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212nf_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb212nf_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb212nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 

        ## H2-H2 interaction 
        movapd nb212nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb212nf_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb212nf_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb212nf_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb212nf_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        movlpd %xmm6,nb212nf_vctot(%esp)

_nb_kernel212nf_ia32_sse2.nb212nf_updateouterdata: 
        ## get n from stack
        movl nb212nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb212nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb212nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb212nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## accumulate total lj energy and update it 
        movapd nb212nf_Vvdwtot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb212nf_Vvdw(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb212nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel212nf_ia32_sse2.nb212nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb212nf_n(%esp)
        jmp _nb_kernel212nf_ia32_sse2.nb212nf_outer
_nb_kernel212nf_ia32_sse2.nb212nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb212nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel212nf_ia32_sse2.nb212nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel212nf_ia32_sse2.nb212nf_threadloop
_nb_kernel212nf_ia32_sse2.nb212nf_end: 
        emms

        movl nb212nf_nouter(%esp),%eax
        movl nb212nf_ninner(%esp),%ebx
        movl nb212nf_outeriter(%ebp),%ecx
        movl nb212nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb212nf_salign(%esp),%eax
        addl %eax,%esp
        addl $792,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret


