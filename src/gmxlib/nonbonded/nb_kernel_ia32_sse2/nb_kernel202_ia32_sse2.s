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



.globl nb_kernel202_ia32_sse2
.globl _nb_kernel202_ia32_sse2
nb_kernel202_ia32_sse2: 
_nb_kernel202_ia32_sse2:        
.set nb202_p_nri, 8
.set nb202_iinr, 12
.set nb202_jindex, 16
.set nb202_jjnr, 20
.set nb202_shift, 24
.set nb202_shiftvec, 28
.set nb202_fshift, 32
.set nb202_gid, 36
.set nb202_pos, 40
.set nb202_faction, 44
.set nb202_charge, 48
.set nb202_p_facel, 52
.set nb202_argkrf, 56
.set nb202_argcrf, 60
.set nb202_Vc, 64
.set nb202_type, 68
.set nb202_p_ntype, 72
.set nb202_vdwparam, 76
.set nb202_Vvdw, 80
.set nb202_p_tabscale, 84
.set nb202_VFtab, 88
.set nb202_invsqrta, 92
.set nb202_dvda, 96
.set nb202_p_gbtabscale, 100
.set nb202_GBtab, 104
.set nb202_p_nthreads, 108
.set nb202_count, 112
.set nb202_mtx, 116
.set nb202_outeriter, 120
.set nb202_inneriter, 124
.set nb202_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb202_ixO, 0
.set nb202_iyO, 16
.set nb202_izO, 32
.set nb202_ixH1, 48
.set nb202_iyH1, 64
.set nb202_izH1, 80
.set nb202_ixH2, 96
.set nb202_iyH2, 112
.set nb202_izH2, 128
.set nb202_jxO, 144
.set nb202_jyO, 160
.set nb202_jzO, 176
.set nb202_jxH1, 192
.set nb202_jyH1, 208
.set nb202_jzH1, 224
.set nb202_jxH2, 240
.set nb202_jyH2, 256
.set nb202_jzH2, 272
.set nb202_dxOO, 288
.set nb202_dyOO, 304
.set nb202_dzOO, 320
.set nb202_dxOH1, 336
.set nb202_dyOH1, 352
.set nb202_dzOH1, 368
.set nb202_dxOH2, 384
.set nb202_dyOH2, 400
.set nb202_dzOH2, 416
.set nb202_dxH1O, 432
.set nb202_dyH1O, 448
.set nb202_dzH1O, 464
.set nb202_dxH1H1, 480
.set nb202_dyH1H1, 496
.set nb202_dzH1H1, 512
.set nb202_dxH1H2, 528
.set nb202_dyH1H2, 544
.set nb202_dzH1H2, 560
.set nb202_dxH2O, 576
.set nb202_dyH2O, 592
.set nb202_dzH2O, 608
.set nb202_dxH2H1, 624
.set nb202_dyH2H1, 640
.set nb202_dzH2H1, 656
.set nb202_dxH2H2, 672
.set nb202_dyH2H2, 688
.set nb202_dzH2H2, 704
.set nb202_qqOO, 720
.set nb202_qqOH, 736
.set nb202_qqHH, 752
.set nb202_vctot, 768
.set nb202_fixO, 784
.set nb202_fiyO, 800
.set nb202_fizO, 816
.set nb202_fixH1, 832
.set nb202_fiyH1, 848
.set nb202_fizH1, 864
.set nb202_fixH2, 880
.set nb202_fiyH2, 896
.set nb202_fizH2, 912
.set nb202_fjxO, 928
.set nb202_fjyO, 944
.set nb202_fjzO, 960
.set nb202_fjxH1, 976
.set nb202_fjyH1, 992
.set nb202_fjzH1, 1008
.set nb202_fjxH2, 1024
.set nb202_fjyH2, 1040
.set nb202_fjzH2, 1056
.set nb202_half, 1072
.set nb202_three, 1088
.set nb202_rsqOO, 1104
.set nb202_rsqOH1, 1120
.set nb202_rsqOH2, 1136
.set nb202_rsqH1O, 1152
.set nb202_rsqH1H1, 1168
.set nb202_rsqH1H2, 1184
.set nb202_rsqH2O, 1200
.set nb202_rsqH2H1, 1216
.set nb202_rsqH2H2, 1232
.set nb202_rinvOO, 1248
.set nb202_rinvOH1, 1264
.set nb202_rinvOH2, 1280
.set nb202_rinvH1O, 1296
.set nb202_rinvH1H1, 1312
.set nb202_rinvH1H2, 1328
.set nb202_rinvH2O, 1344
.set nb202_rinvH2H1, 1360
.set nb202_rinvH2H2, 1376
.set nb202_two, 1392
.set nb202_krf, 1408
.set nb202_crf, 1424
.set nb202_is3, 1440
.set nb202_ii3, 1444
.set nb202_innerjjnr, 1448
.set nb202_innerk, 1452
.set nb202_n, 1456
.set nb202_nn1, 1460
.set nb202_nri, 1464
.set nb202_nouter, 1468
.set nb202_ninner, 1472
.set nb202_salign, 1476
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $1480,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb202_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb202_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb202_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb202_nouter(%esp)
        movl %eax,nb202_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb202_half(%esp)
        movl %ebx,nb202_half+4(%esp)
        movsd nb202_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb202_half(%esp)
        movapd %xmm2,nb202_two(%esp)
        movapd %xmm3,nb202_three(%esp)

        movl nb202_argkrf(%ebp),%esi
        movl nb202_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb202_krf(%esp)
        movapd %xmm6,nb202_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb202_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb202_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb202_p_facel(%ebp),%esi
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
        movapd %xmm3,nb202_qqOO(%esp)
        movapd %xmm4,nb202_qqOH(%esp)
        movapd %xmm5,nb202_qqHH(%esp)

_nb_kernel202_ia32_sse2.nb202_threadloop: 
        movl  nb202_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel202_ia32_sse2.nb202_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel202_ia32_sse2.nb202_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb202_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb202_n(%esp)
        movl %ebx,nb202_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel202_ia32_sse2.nb202_outerstart
        jmp _nb_kernel202_ia32_sse2.nb202_end

_nb_kernel202_ia32_sse2.nb202_outerstart: 
        ## ebx contains number of outer iterations
        addl nb202_nouter(%esp),%ebx
        movl %ebx,nb202_nouter(%esp)

_nb_kernel202_ia32_sse2.nb202_outer: 
        movl  nb202_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb202_is3(%esp)      ## store is3 

        movl  nb202_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb202_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb202_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb202_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb202_ixO(%esp)
        movapd %xmm4,nb202_iyO(%esp)
        movapd %xmm5,nb202_izO(%esp)

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
        movapd %xmm0,nb202_ixH1(%esp)
        movapd %xmm1,nb202_iyH1(%esp)
        movapd %xmm2,nb202_izH1(%esp)
        movapd %xmm3,nb202_ixH2(%esp)
        movapd %xmm4,nb202_iyH2(%esp)
        movapd %xmm5,nb202_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb202_vctot(%esp)
        movapd %xmm4,nb202_fixO(%esp)
        movapd %xmm4,nb202_fiyO(%esp)
        movapd %xmm4,nb202_fizO(%esp)
        movapd %xmm4,nb202_fixH1(%esp)
        movapd %xmm4,nb202_fiyH1(%esp)
        movapd %xmm4,nb202_fizH1(%esp)
        movapd %xmm4,nb202_fixH2(%esp)
        movapd %xmm4,nb202_fiyH2(%esp)
        movapd %xmm4,nb202_fizH2(%esp)

        movl  nb202_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb202_pos(%ebp),%esi
        movl  nb202_faction(%ebp),%edi
        movl  nb202_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb202_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb202_ninner(%esp),%ecx
        movl  %ecx,nb202_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb202_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel202_ia32_sse2.nb202_unroll_loop
        jmp   _nb_kernel202_ia32_sse2.nb202_checksingle
_nb_kernel202_ia32_sse2.nb202_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb202_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb202_innerjjnr(%esp)                   ## advance pointer (unrolled 2) 

        movl nb202_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb202_jxO(%esp)
        movapd  %xmm3,nb202_jyO(%esp)
        movapd  %xmm4,nb202_jzO(%esp)
        movapd  %xmm5,nb202_jxH1(%esp)
        movapd  %xmm6,nb202_jyH1(%esp)
        movapd  %xmm7,nb202_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb202_jxH2(%esp)
        movapd  %xmm3,nb202_jyH2(%esp)
        movapd  %xmm4,nb202_jzH2(%esp)

        movapd nb202_ixO(%esp),%xmm0
        movapd nb202_iyO(%esp),%xmm1
        movapd nb202_izO(%esp),%xmm2
        movapd nb202_ixO(%esp),%xmm3
        movapd nb202_iyO(%esp),%xmm4
        movapd nb202_izO(%esp),%xmm5
        subpd  nb202_jxO(%esp),%xmm0
        subpd  nb202_jyO(%esp),%xmm1
        subpd  nb202_jzO(%esp),%xmm2
        subpd  nb202_jxH1(%esp),%xmm3
        subpd  nb202_jyH1(%esp),%xmm4
        subpd  nb202_jzH1(%esp),%xmm5
        movapd %xmm0,nb202_dxOO(%esp)
        movapd %xmm1,nb202_dyOO(%esp)
        movapd %xmm2,nb202_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb202_dxOH1(%esp)
        movapd %xmm4,nb202_dyOH1(%esp)
        movapd %xmm5,nb202_dzOH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb202_rsqOO(%esp)
        movapd %xmm3,nb202_rsqOH1(%esp)

        movapd nb202_ixO(%esp),%xmm0
        movapd nb202_iyO(%esp),%xmm1
        movapd nb202_izO(%esp),%xmm2
        movapd nb202_ixH1(%esp),%xmm3
        movapd nb202_iyH1(%esp),%xmm4
        movapd nb202_izH1(%esp),%xmm5
        subpd  nb202_jxH2(%esp),%xmm0
        subpd  nb202_jyH2(%esp),%xmm1
        subpd  nb202_jzH2(%esp),%xmm2
        subpd  nb202_jxO(%esp),%xmm3
        subpd  nb202_jyO(%esp),%xmm4
        subpd  nb202_jzO(%esp),%xmm5
        movapd %xmm0,nb202_dxOH2(%esp)
        movapd %xmm1,nb202_dyOH2(%esp)
        movapd %xmm2,nb202_dzOH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb202_dxH1O(%esp)
        movapd %xmm4,nb202_dyH1O(%esp)
        movapd %xmm5,nb202_dzH1O(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb202_rsqOH2(%esp)
        movapd %xmm3,nb202_rsqH1O(%esp)

        movapd nb202_ixH1(%esp),%xmm0
        movapd nb202_iyH1(%esp),%xmm1
        movapd nb202_izH1(%esp),%xmm2
        movapd nb202_ixH1(%esp),%xmm3
        movapd nb202_iyH1(%esp),%xmm4
        movapd nb202_izH1(%esp),%xmm5
        subpd  nb202_jxH1(%esp),%xmm0
        subpd  nb202_jyH1(%esp),%xmm1
        subpd  nb202_jzH1(%esp),%xmm2
        subpd  nb202_jxH2(%esp),%xmm3
        subpd  nb202_jyH2(%esp),%xmm4
        subpd  nb202_jzH2(%esp),%xmm5
        movapd %xmm0,nb202_dxH1H1(%esp)
        movapd %xmm1,nb202_dyH1H1(%esp)
        movapd %xmm2,nb202_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb202_dxH1H2(%esp)
        movapd %xmm4,nb202_dyH1H2(%esp)
        movapd %xmm5,nb202_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb202_rsqH1H1(%esp)
        movapd %xmm3,nb202_rsqH1H2(%esp)

        movapd nb202_ixH2(%esp),%xmm0
        movapd nb202_iyH2(%esp),%xmm1
        movapd nb202_izH2(%esp),%xmm2
        movapd nb202_ixH2(%esp),%xmm3
        movapd nb202_iyH2(%esp),%xmm4
        movapd nb202_izH2(%esp),%xmm5
        subpd  nb202_jxO(%esp),%xmm0
        subpd  nb202_jyO(%esp),%xmm1
        subpd  nb202_jzO(%esp),%xmm2
        subpd  nb202_jxH1(%esp),%xmm3
        subpd  nb202_jyH1(%esp),%xmm4
        subpd  nb202_jzH1(%esp),%xmm5
        movapd %xmm0,nb202_dxH2O(%esp)
        movapd %xmm1,nb202_dyH2O(%esp)
        movapd %xmm2,nb202_dzH2O(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb202_dxH2H1(%esp)
        movapd %xmm4,nb202_dyH2H1(%esp)
        movapd %xmm5,nb202_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb202_rsqH2O(%esp)
        movapd %xmm4,nb202_rsqH2H1(%esp)

        movapd nb202_ixH2(%esp),%xmm0
        movapd nb202_iyH2(%esp),%xmm1
        movapd nb202_izH2(%esp),%xmm2
        subpd  nb202_jxH2(%esp),%xmm0
        subpd  nb202_jyH2(%esp),%xmm1
        subpd  nb202_jzH2(%esp),%xmm2
        movapd %xmm0,nb202_dxH2H2(%esp)
        movapd %xmm1,nb202_dyH2H2(%esp)
        movapd %xmm2,nb202_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb202_rsqH2H2(%esp)

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
        movapd  nb202_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202_half(%esp),%xmm3   ## iter1 
        mulpd   nb202_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202_half(%esp),%xmm1   ## rinv 
        mulpd   nb202_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202_rinvH2H2(%esp)
        movapd %xmm5,nb202_rinvH2H1(%esp)

        movapd nb202_rsqOO(%esp),%xmm0
        movapd nb202_rsqOH1(%esp),%xmm4
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
        movapd  nb202_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb202_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202_half(%esp),%xmm1   ## rinv 
        mulpd   nb202_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb202_rinvOO(%esp)
        movapd %xmm5,nb202_rinvOH1(%esp)

        movapd nb202_rsqOH2(%esp),%xmm0
        movapd nb202_rsqH1O(%esp),%xmm4
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
        movapd  nb202_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202_half(%esp),%xmm3   ## iter1 
        mulpd   nb202_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202_half(%esp),%xmm1   ## rinv 
        mulpd   nb202_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202_rinvOH2(%esp)
        movapd %xmm5,nb202_rinvH1O(%esp)

        movapd nb202_rsqH1H1(%esp),%xmm0
        movapd nb202_rsqH1H2(%esp),%xmm4
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
        movapd  nb202_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202_half(%esp),%xmm3   ## iter1a 
        mulpd   nb202_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202_half(%esp),%xmm1   ## rinv 
        mulpd   nb202_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202_rinvH1H1(%esp)
        movapd %xmm5,nb202_rinvH1H2(%esp)

        movapd nb202_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb202_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb202_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb202_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb202_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb202_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb202_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        mulpd  %xmm0,%xmm0      ## rinvsq 
        mulpd  nb202_rsqOO(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addpd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subpd  nb202_crf(%esp),%xmm6

        mulpd  nb202_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqOO(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addpd  nb202_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulpd  %xmm7,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb202_dxOO(%esp),%xmm0
        mulpd nb202_dyOO(%esp),%xmm1
        mulpd nb202_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixO(%esp),%xmm0
        addpd nb202_fiyO(%esp),%xmm1
        addpd nb202_fizO(%esp),%xmm2
        movapd %xmm3,nb202_fjxO(%esp)
        movapd %xmm4,nb202_fjyO(%esp)
        movapd %xmm5,nb202_fjzO(%esp)
        movapd %xmm0,nb202_fixO(%esp)
        movapd %xmm1,nb202_fiyO(%esp)
        movapd %xmm2,nb202_fizO(%esp)

        ## O-H1 interaction 
        movapd nb202_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb202_rsqOH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd  %xmm0,%xmm0
        subpd  nb202_crf(%esp),%xmm4
        mulpd  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH1  
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb202_dxOH1(%esp),%xmm0
        mulpd nb202_dyOH1(%esp),%xmm1
        mulpd nb202_dzOH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixO(%esp),%xmm0
        addpd nb202_fiyO(%esp),%xmm1
        addpd nb202_fizO(%esp),%xmm2
        movapd %xmm3,nb202_fjxH1(%esp)
        movapd %xmm4,nb202_fjyH1(%esp)
        movapd %xmm5,nb202_fjzH1(%esp)
        movapd %xmm0,nb202_fixO(%esp)
        movapd %xmm1,nb202_fiyO(%esp)
        movapd %xmm2,nb202_fizO(%esp)

        ## O-H2 interaction  
        movapd nb202_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb202_rsqOH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb202_crf(%esp),%xmm4
        mulpd  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb202_dxOH2(%esp),%xmm0
        mulpd nb202_dyOH2(%esp),%xmm1
        mulpd nb202_dzOH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixO(%esp),%xmm0
        addpd nb202_fiyO(%esp),%xmm1
        addpd nb202_fizO(%esp),%xmm2
        movapd %xmm3,nb202_fjxH2(%esp)
        movapd %xmm4,nb202_fjyH2(%esp)
        movapd %xmm5,nb202_fjzH2(%esp)
        movapd %xmm0,nb202_fixO(%esp)
        movapd %xmm1,nb202_fiyO(%esp)
        movapd %xmm2,nb202_fizO(%esp)

        ## H1-O interaction 
        movapd nb202_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb202_rsqH1O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb202_crf(%esp),%xmm4
        mulpd  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxO(%esp),%xmm3
        movapd nb202_fjyO(%esp),%xmm4
        movapd nb202_fjzO(%esp),%xmm5
        mulpd nb202_dxH1O(%esp),%xmm0
        mulpd nb202_dyH1O(%esp),%xmm1
        mulpd nb202_dzH1O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixH1(%esp),%xmm0
        addpd nb202_fiyH1(%esp),%xmm1
        addpd nb202_fizH1(%esp),%xmm2
        movapd %xmm3,nb202_fjxO(%esp)
        movapd %xmm4,nb202_fjyO(%esp)
        movapd %xmm5,nb202_fjzO(%esp)
        movapd %xmm0,nb202_fixH1(%esp)
        movapd %xmm1,nb202_fiyH1(%esp)
        movapd %xmm2,nb202_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb202_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb202_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb202_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxH1(%esp),%xmm3
        movapd nb202_fjyH1(%esp),%xmm4
        movapd nb202_fjzH1(%esp),%xmm5
        mulpd nb202_dxH1H1(%esp),%xmm0
        mulpd nb202_dyH1H1(%esp),%xmm1
        mulpd nb202_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixH1(%esp),%xmm0
        addpd nb202_fiyH1(%esp),%xmm1
        addpd nb202_fizH1(%esp),%xmm2
        movapd %xmm3,nb202_fjxH1(%esp)
        movapd %xmm4,nb202_fjyH1(%esp)
        movapd %xmm5,nb202_fjzH1(%esp)
        movapd %xmm0,nb202_fixH1(%esp)
        movapd %xmm1,nb202_fiyH1(%esp)
        movapd %xmm2,nb202_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb202_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb202_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulpd %xmm0,%xmm0
        subpd  nb202_crf(%esp),%xmm4
        mulpd  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxH2(%esp),%xmm3
        movapd nb202_fjyH2(%esp),%xmm4
        movapd nb202_fjzH2(%esp),%xmm5
        mulpd nb202_dxH1H2(%esp),%xmm0
        mulpd nb202_dyH1H2(%esp),%xmm1
        mulpd nb202_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixH1(%esp),%xmm0
        addpd nb202_fiyH1(%esp),%xmm1
        addpd nb202_fizH1(%esp),%xmm2
        movapd %xmm3,nb202_fjxH2(%esp)
        movapd %xmm4,nb202_fjyH2(%esp)
        movapd %xmm5,nb202_fjzH2(%esp)
        movapd %xmm0,nb202_fixH1(%esp)
        movapd %xmm1,nb202_fiyH1(%esp)
        movapd %xmm2,nb202_fizH1(%esp)

        ## H2-O interaction 
        movapd nb202_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb202_rsqH2O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb202_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxO(%esp),%xmm3
        movapd nb202_fjyO(%esp),%xmm4
        movapd nb202_fjzO(%esp),%xmm5
        mulpd nb202_dxH2O(%esp),%xmm0
        mulpd nb202_dyH2O(%esp),%xmm1
        mulpd nb202_dzH2O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixH2(%esp),%xmm0
        addpd nb202_fiyH2(%esp),%xmm1
        addpd nb202_fizH2(%esp),%xmm2
        movapd %xmm3,nb202_fjxO(%esp)
        movapd %xmm4,nb202_fjyO(%esp)
        movapd %xmm5,nb202_fjzO(%esp)
        movapd %xmm0,nb202_fixH2(%esp)
        movapd %xmm1,nb202_fiyH2(%esp)
        movapd %xmm2,nb202_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb202_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb202_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb202_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxH1(%esp),%xmm3
        movapd nb202_fjyH1(%esp),%xmm4
        movapd nb202_fjzH1(%esp),%xmm5
        mulpd nb202_dxH2H1(%esp),%xmm0
        mulpd nb202_dyH2H1(%esp),%xmm1
        mulpd nb202_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixH2(%esp),%xmm0
        addpd nb202_fiyH2(%esp),%xmm1
        addpd nb202_fizH2(%esp),%xmm2
        movapd %xmm3,nb202_fjxH1(%esp)
        movapd %xmm4,nb202_fjyH1(%esp)
        movapd %xmm5,nb202_fjzH1(%esp)
        movapd %xmm0,nb202_fixH2(%esp)
        movapd %xmm1,nb202_fiyH2(%esp)
        movapd %xmm2,nb202_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb202_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulpd  nb202_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addpd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subpd  nb202_crf(%esp),%xmm4
        mulpd %xmm0,%xmm0
        mulpd  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulpd  nb202_two(%esp),%xmm5
        subpd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulpd  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addpd  %xmm4,%xmm6      ## add to local vctot 
        mulpd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd %xmm0,%xmm1
        movapd %xmm6,nb202_vctot(%esp)
        movapd %xmm0,%xmm2

        movapd nb202_fjxH2(%esp),%xmm3
        movapd nb202_fjyH2(%esp),%xmm4
        movapd nb202_fjzH2(%esp),%xmm5
        mulpd nb202_dxH2H2(%esp),%xmm0
        mulpd nb202_dyH2H2(%esp),%xmm1
        mulpd nb202_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb202_fixH2(%esp),%xmm0
        addpd nb202_fiyH2(%esp),%xmm1
        addpd nb202_fizH2(%esp),%xmm2
        movapd %xmm3,nb202_fjxH2(%esp)
        movapd %xmm4,nb202_fjyH2(%esp)
        movapd %xmm5,nb202_fjzH2(%esp)
        movapd %xmm0,nb202_fixH2(%esp)
        movapd %xmm1,nb202_fiyH2(%esp)
        movapd %xmm2,nb202_fizH2(%esp)

        movl nb202_faction(%ebp),%edi

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
        addpd nb202_fjxO(%esp),%xmm0
        addpd nb202_fjyO(%esp),%xmm1
        addpd nb202_fjzO(%esp),%xmm2
        addpd nb202_fjxH1(%esp),%xmm3
        addpd nb202_fjyH1(%esp),%xmm4
        addpd nb202_fjzH1(%esp),%xmm5
        addpd nb202_fjxH2(%esp),%xmm6
        addpd nb202_fjyH2(%esp),%xmm7
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
        addpd nb202_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)
        movhpd %xmm0,64(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb202_innerk(%esp)
        jl    _nb_kernel202_ia32_sse2.nb202_checksingle
        jmp   _nb_kernel202_ia32_sse2.nb202_unroll_loop
_nb_kernel202_ia32_sse2.nb202_checksingle: 
        movl  nb202_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel202_ia32_sse2.nb202_dosingle
        jmp   _nb_kernel202_ia32_sse2.nb202_updateouterdata
_nb_kernel202_ia32_sse2.nb202_dosingle: 
        movl  nb202_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb202_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb202_jxO(%esp)
        movapd  %xmm3,nb202_jyO(%esp)
        movapd  %xmm4,nb202_jzO(%esp)
        movapd  %xmm5,nb202_jxH1(%esp)
        movapd  %xmm6,nb202_jyH1(%esp)
        movapd  %xmm7,nb202_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb202_jxH2(%esp)
        movapd  %xmm3,nb202_jyH2(%esp)
        movapd  %xmm4,nb202_jzH2(%esp)

        movapd nb202_ixO(%esp),%xmm0
        movapd nb202_iyO(%esp),%xmm1
        movapd nb202_izO(%esp),%xmm2
        movapd nb202_ixO(%esp),%xmm3
        movapd nb202_iyO(%esp),%xmm4
        movapd nb202_izO(%esp),%xmm5
        subsd  nb202_jxO(%esp),%xmm0
        subsd  nb202_jyO(%esp),%xmm1
        subsd  nb202_jzO(%esp),%xmm2
        subsd  nb202_jxH1(%esp),%xmm3
        subsd  nb202_jyH1(%esp),%xmm4
        subsd  nb202_jzH1(%esp),%xmm5
        movapd %xmm0,nb202_dxOO(%esp)
        movapd %xmm1,nb202_dyOO(%esp)
        movapd %xmm2,nb202_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb202_dxOH1(%esp)
        movapd %xmm4,nb202_dyOH1(%esp)
        movapd %xmm5,nb202_dzOH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb202_rsqOO(%esp)
        movapd %xmm3,nb202_rsqOH1(%esp)

        movapd nb202_ixO(%esp),%xmm0
        movapd nb202_iyO(%esp),%xmm1
        movapd nb202_izO(%esp),%xmm2
        movapd nb202_ixH1(%esp),%xmm3
        movapd nb202_iyH1(%esp),%xmm4
        movapd nb202_izH1(%esp),%xmm5
        subsd  nb202_jxH2(%esp),%xmm0
        subsd  nb202_jyH2(%esp),%xmm1
        subsd  nb202_jzH2(%esp),%xmm2
        subsd  nb202_jxO(%esp),%xmm3
        subsd  nb202_jyO(%esp),%xmm4
        subsd  nb202_jzO(%esp),%xmm5
        movapd %xmm0,nb202_dxOH2(%esp)
        movapd %xmm1,nb202_dyOH2(%esp)
        movapd %xmm2,nb202_dzOH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb202_dxH1O(%esp)
        movapd %xmm4,nb202_dyH1O(%esp)
        movapd %xmm5,nb202_dzH1O(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb202_rsqOH2(%esp)
        movapd %xmm3,nb202_rsqH1O(%esp)

        movapd nb202_ixH1(%esp),%xmm0
        movapd nb202_iyH1(%esp),%xmm1
        movapd nb202_izH1(%esp),%xmm2
        movapd nb202_ixH1(%esp),%xmm3
        movapd nb202_iyH1(%esp),%xmm4
        movapd nb202_izH1(%esp),%xmm5
        subsd  nb202_jxH1(%esp),%xmm0
        subsd  nb202_jyH1(%esp),%xmm1
        subsd  nb202_jzH1(%esp),%xmm2
        subsd  nb202_jxH2(%esp),%xmm3
        subsd  nb202_jyH2(%esp),%xmm4
        subsd  nb202_jzH2(%esp),%xmm5
        movapd %xmm0,nb202_dxH1H1(%esp)
        movapd %xmm1,nb202_dyH1H1(%esp)
        movapd %xmm2,nb202_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb202_dxH1H2(%esp)
        movapd %xmm4,nb202_dyH1H2(%esp)
        movapd %xmm5,nb202_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb202_rsqH1H1(%esp)
        movapd %xmm3,nb202_rsqH1H2(%esp)

        movapd nb202_ixH2(%esp),%xmm0
        movapd nb202_iyH2(%esp),%xmm1
        movapd nb202_izH2(%esp),%xmm2
        movapd nb202_ixH2(%esp),%xmm3
        movapd nb202_iyH2(%esp),%xmm4
        movapd nb202_izH2(%esp),%xmm5
        subsd  nb202_jxO(%esp),%xmm0
        subsd  nb202_jyO(%esp),%xmm1
        subsd  nb202_jzO(%esp),%xmm2
        subsd  nb202_jxH1(%esp),%xmm3
        subsd  nb202_jyH1(%esp),%xmm4
        subsd  nb202_jzH1(%esp),%xmm5
        movapd %xmm0,nb202_dxH2O(%esp)
        movapd %xmm1,nb202_dyH2O(%esp)
        movapd %xmm2,nb202_dzH2O(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb202_dxH2H1(%esp)
        movapd %xmm4,nb202_dyH2H1(%esp)
        movapd %xmm5,nb202_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb202_rsqH2O(%esp)
        movapd %xmm4,nb202_rsqH2H1(%esp)

        movapd nb202_ixH2(%esp),%xmm0
        movapd nb202_iyH2(%esp),%xmm1
        movapd nb202_izH2(%esp),%xmm2
        subsd  nb202_jxH2(%esp),%xmm0
        subsd  nb202_jyH2(%esp),%xmm1
        subsd  nb202_jzH2(%esp),%xmm2
        movapd %xmm0,nb202_dxH2H2(%esp)
        movapd %xmm1,nb202_dyH2H2(%esp)
        movapd %xmm2,nb202_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb202_rsqH2H2(%esp)

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
        movapd  nb202_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202_half(%esp),%xmm3   ## iter1 
        mulsd   nb202_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202_half(%esp),%xmm1   ## rinv 
        mulsd   nb202_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202_rinvH2H2(%esp)
        movapd %xmm5,nb202_rinvH2H1(%esp)

        movapd nb202_rsqOO(%esp),%xmm0
        movapd nb202_rsqOH1(%esp),%xmm4
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
        movapd  nb202_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb202_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202_half(%esp),%xmm1   ## rinv 
        mulsd   nb202_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb202_rinvOO(%esp)
        movapd %xmm5,nb202_rinvOH1(%esp)

        movapd nb202_rsqOH2(%esp),%xmm0
        movapd nb202_rsqH1O(%esp),%xmm4
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
        movapd  nb202_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202_half(%esp),%xmm3   ## iter1 
        mulsd   nb202_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202_half(%esp),%xmm1   ## rinv 
        mulsd   nb202_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202_rinvOH2(%esp)
        movapd %xmm5,nb202_rinvH1O(%esp)

        movapd nb202_rsqH1H1(%esp),%xmm0
        movapd nb202_rsqH1H2(%esp),%xmm4
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
        movapd  nb202_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202_half(%esp),%xmm3   ## iter1a 
        mulsd   nb202_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202_half(%esp),%xmm1   ## rinv 
        mulsd   nb202_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202_rinvH1H1(%esp)
        movapd %xmm5,nb202_rinvH1H2(%esp)

        movapd nb202_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb202_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb202_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb202_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb202_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb202_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb202_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        mulsd  %xmm0,%xmm0
        movapd %xmm0,%xmm1
        mulsd  %xmm0,%xmm1
        mulsd  %xmm0,%xmm1      ## xmm1=rinvsix 
        mulsd  nb202_rsqOO(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm6
        addsd  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subsd  nb202_crf(%esp),%xmm6

        mulsd  nb202_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqOO(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addsd  nb202_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulsd  %xmm7,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb202_dxOO(%esp),%xmm0
        mulsd nb202_dyOO(%esp),%xmm1
        mulsd nb202_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixO(%esp),%xmm0
        addsd nb202_fiyO(%esp),%xmm1
        addsd nb202_fizO(%esp),%xmm2
        movlpd %xmm3,nb202_fjxO(%esp)
        movlpd %xmm4,nb202_fjyO(%esp)
        movlpd %xmm5,nb202_fjzO(%esp)
        movlpd %xmm0,nb202_fixO(%esp)
        movlpd %xmm1,nb202_fiyO(%esp)
        movlpd %xmm2,nb202_fizO(%esp)

        ## O-H1 interaction 
        movapd nb202_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb202_rsqOH1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb202_crf(%esp),%xmm4
        mulsd  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH1  
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb202_dxOH1(%esp),%xmm0
        mulsd nb202_dyOH1(%esp),%xmm1
        mulsd nb202_dzOH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixO(%esp),%xmm0
        addsd nb202_fiyO(%esp),%xmm1
        addsd nb202_fizO(%esp),%xmm2
        movlpd %xmm3,nb202_fjxH1(%esp)
        movlpd %xmm4,nb202_fjyH1(%esp)
        movlpd %xmm5,nb202_fjzH1(%esp)
        movlpd %xmm0,nb202_fixO(%esp)
        movlpd %xmm1,nb202_fiyO(%esp)
        movlpd %xmm2,nb202_fizO(%esp)

        ## O-H2 interaction  
        movapd nb202_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb202_rsqOH2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd  %xmm0,%xmm0
        subsd  nb202_crf(%esp),%xmm4
        mulsd  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb202_dxOH2(%esp),%xmm0
        mulsd nb202_dyOH2(%esp),%xmm1
        mulsd nb202_dzOH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixO(%esp),%xmm0
        addsd nb202_fiyO(%esp),%xmm1
        addsd nb202_fizO(%esp),%xmm2
        movlpd %xmm3,nb202_fjxH2(%esp)
        movlpd %xmm4,nb202_fjyH2(%esp)
        movlpd %xmm5,nb202_fjzH2(%esp)
        movlpd %xmm0,nb202_fixO(%esp)
        movlpd %xmm1,nb202_fiyO(%esp)
        movlpd %xmm2,nb202_fizO(%esp)

        ## H1-O interaction 
        movapd nb202_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb202_rsqH1O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=rinv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb202_crf(%esp),%xmm4
        mulsd  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxO(%esp),%xmm3
        movapd nb202_fjyO(%esp),%xmm4
        movapd nb202_fjzO(%esp),%xmm5
        mulsd nb202_dxH1O(%esp),%xmm0
        mulsd nb202_dyH1O(%esp),%xmm1
        mulsd nb202_dzH1O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixH1(%esp),%xmm0
        addsd nb202_fiyH1(%esp),%xmm1
        addsd nb202_fizH1(%esp),%xmm2
        movlpd %xmm3,nb202_fjxO(%esp)
        movlpd %xmm4,nb202_fjyO(%esp)
        movlpd %xmm5,nb202_fjzO(%esp)
        movlpd %xmm0,nb202_fixH1(%esp)
        movlpd %xmm1,nb202_fiyH1(%esp)
        movlpd %xmm2,nb202_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb202_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb202_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb202_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxH1(%esp),%xmm3
        movapd nb202_fjyH1(%esp),%xmm4
        movapd nb202_fjzH1(%esp),%xmm5
        mulsd nb202_dxH1H1(%esp),%xmm0
        mulsd nb202_dyH1H1(%esp),%xmm1
        mulsd nb202_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixH1(%esp),%xmm0
        addsd nb202_fiyH1(%esp),%xmm1
        addsd nb202_fizH1(%esp),%xmm2
        movlpd %xmm3,nb202_fjxH1(%esp)
        movlpd %xmm4,nb202_fjyH1(%esp)
        movlpd %xmm5,nb202_fjzH1(%esp)
        movlpd %xmm0,nb202_fixH1(%esp)
        movlpd %xmm1,nb202_fiyH1(%esp)
        movlpd %xmm2,nb202_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb202_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb202_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        mulsd %xmm0,%xmm0
        subsd  nb202_crf(%esp),%xmm4
        mulsd  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxH2(%esp),%xmm3
        movapd nb202_fjyH2(%esp),%xmm4
        movapd nb202_fjzH2(%esp),%xmm5
        mulsd nb202_dxH1H2(%esp),%xmm0
        mulsd nb202_dyH1H2(%esp),%xmm1
        mulsd nb202_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixH1(%esp),%xmm0
        addsd nb202_fiyH1(%esp),%xmm1
        addsd nb202_fizH1(%esp),%xmm2
        movlpd %xmm3,nb202_fjxH2(%esp)
        movlpd %xmm4,nb202_fjyH2(%esp)
        movlpd %xmm5,nb202_fjzH2(%esp)
        movlpd %xmm0,nb202_fixH1(%esp)
        movlpd %xmm1,nb202_fiyH1(%esp)
        movlpd %xmm2,nb202_fizH1(%esp)

        ## H2-O interaction 
        movapd nb202_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb202_rsqH2O(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb202_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxO(%esp),%xmm3
        movapd nb202_fjyO(%esp),%xmm4
        movapd nb202_fjzO(%esp),%xmm5
        mulsd nb202_dxH2O(%esp),%xmm0
        mulsd nb202_dyH2O(%esp),%xmm1
        mulsd nb202_dzH2O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixH2(%esp),%xmm0
        addsd nb202_fiyH2(%esp),%xmm1
        addsd nb202_fizH2(%esp),%xmm2
        movlpd %xmm3,nb202_fjxO(%esp)
        movlpd %xmm4,nb202_fjyO(%esp)
        movlpd %xmm5,nb202_fjzO(%esp)
        movlpd %xmm0,nb202_fixH2(%esp)
        movlpd %xmm1,nb202_fiyH2(%esp)
        movlpd %xmm2,nb202_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb202_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb202_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb202_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd nb202_fjxH1(%esp),%xmm3
        movapd nb202_fjyH1(%esp),%xmm4
        movapd nb202_fjzH1(%esp),%xmm5
        mulsd nb202_dxH2H1(%esp),%xmm0
        mulsd nb202_dyH2H1(%esp),%xmm1
        mulsd nb202_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixH2(%esp),%xmm0
        addsd nb202_fiyH2(%esp),%xmm1
        addsd nb202_fizH2(%esp),%xmm2
        movlpd %xmm3,nb202_fjxH1(%esp)
        movlpd %xmm4,nb202_fjyH1(%esp)
        movlpd %xmm5,nb202_fjzH1(%esp)
        movlpd %xmm0,nb202_fixH2(%esp)
        movlpd %xmm1,nb202_fiyH2(%esp)
        movlpd %xmm2,nb202_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb202_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm7      ## xmm7=rinv 
        movapd nb202_krf(%esp),%xmm5
        movapd %xmm0,%xmm1
        mulsd  nb202_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movapd %xmm5,%xmm4
        addsd  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subsd  nb202_crf(%esp),%xmm4
        mulsd %xmm0,%xmm0
        mulsd  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq) 
        mulsd  nb202_two(%esp),%xmm5
        subsd  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulsd  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addsd  %xmm4,%xmm6      ## add to local vctot 
        mulsd %xmm7,%xmm0       ## fsOH2 
        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        movapd %xmm0,%xmm1
        movlpd %xmm6,nb202_vctot(%esp)
        movapd %xmm0,%xmm2

        movapd nb202_fjxH2(%esp),%xmm3
        movapd nb202_fjyH2(%esp),%xmm4
        movapd nb202_fjzH2(%esp),%xmm5
        mulsd nb202_dxH2H2(%esp),%xmm0
        mulsd nb202_dyH2H2(%esp),%xmm1
        mulsd nb202_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb202_fixH2(%esp),%xmm0
        addsd nb202_fiyH2(%esp),%xmm1
        addsd nb202_fizH2(%esp),%xmm2
        movlpd %xmm3,nb202_fjxH2(%esp)
        movlpd %xmm4,nb202_fjyH2(%esp)
        movlpd %xmm5,nb202_fjzH2(%esp)
        movlpd %xmm0,nb202_fixH2(%esp)
        movlpd %xmm1,nb202_fiyH2(%esp)
        movlpd %xmm2,nb202_fizH2(%esp)

        movl nb202_faction(%ebp),%edi
        ## Did all interactions - now update j forces 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2
        movlpd 24(%edi,%eax,8),%xmm3
        movlpd 32(%edi,%eax,8),%xmm4
        movlpd 40(%edi,%eax,8),%xmm5
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 56(%edi,%eax,8),%xmm7
        addsd nb202_fjxO(%esp),%xmm0
        addsd nb202_fjyO(%esp),%xmm1
        addsd nb202_fjzO(%esp),%xmm2
        addsd nb202_fjxH1(%esp),%xmm3
        addsd nb202_fjyH1(%esp),%xmm4
        addsd nb202_fjzH1(%esp),%xmm5
        addsd nb202_fjxH2(%esp),%xmm6
        addsd nb202_fjyH2(%esp),%xmm7
        movlpd %xmm0,(%edi,%eax,8)
        movlpd %xmm1,8(%edi,%eax,8)
        movlpd %xmm2,16(%edi,%eax,8)
        movlpd %xmm3,24(%edi,%eax,8)
        movlpd %xmm4,32(%edi,%eax,8)
        movlpd %xmm5,40(%edi,%eax,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,56(%edi,%eax,8)

        movlpd 64(%edi,%eax,8),%xmm0
        addsd nb202_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)

_nb_kernel202_ia32_sse2.nb202_updateouterdata: 
        movl  nb202_ii3(%esp),%ecx
        movl  nb202_faction(%ebp),%edi
        movl  nb202_fshift(%ebp),%esi
        movl  nb202_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb202_fixO(%esp),%xmm0
        movapd nb202_fiyO(%esp),%xmm1
        movapd nb202_fizO(%esp),%xmm2

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
        movapd nb202_fixH1(%esp),%xmm0
        movapd nb202_fiyH1(%esp),%xmm1
        movapd nb202_fizH1(%esp),%xmm2

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
        movapd nb202_fixH2(%esp),%xmm0
        movapd nb202_fiyH2(%esp),%xmm1
        movapd nb202_fizH2(%esp),%xmm2

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
        movl nb202_n(%esp),%esi
        ## get group index for i particle 
        movl  nb202_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb202_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb202_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb202_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel202_ia32_sse2.nb202_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb202_n(%esp)
        jmp _nb_kernel202_ia32_sse2.nb202_outer
_nb_kernel202_ia32_sse2.nb202_outerend: 
        ## check if more outer neighborlists remain
        movl  nb202_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel202_ia32_sse2.nb202_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel202_ia32_sse2.nb202_threadloop
_nb_kernel202_ia32_sse2.nb202_end: 
        emms

        movl nb202_nouter(%esp),%eax
        movl nb202_ninner(%esp),%ebx
        movl nb202_outeriter(%ebp),%ecx
        movl nb202_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb202_salign(%esp),%eax
        addl %eax,%esp
        addl $1480,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret




.globl nb_kernel202nf_ia32_sse2
.globl _nb_kernel202nf_ia32_sse2
nb_kernel202nf_ia32_sse2:       
_nb_kernel202nf_ia32_sse2:      
.set nb202nf_p_nri, 8
.set nb202nf_iinr, 12
.set nb202nf_jindex, 16
.set nb202nf_jjnr, 20
.set nb202nf_shift, 24
.set nb202nf_shiftvec, 28
.set nb202nf_fshift, 32
.set nb202nf_gid, 36
.set nb202nf_pos, 40
.set nb202nf_faction, 44
.set nb202nf_charge, 48
.set nb202nf_p_facel, 52
.set nb202nf_argkrf, 56
.set nb202nf_argcrf, 60
.set nb202nf_Vc, 64
.set nb202nf_type, 68
.set nb202nf_p_ntype, 72
.set nb202nf_vdwparam, 76
.set nb202nf_Vvdw, 80
.set nb202nf_p_tabscale, 84
.set nb202nf_VFtab, 88
.set nb202nf_invsqrta, 92
.set nb202nf_dvda, 96
.set nb202nf_p_gbtabscale, 100
.set nb202nf_GBtab, 104
.set nb202nf_p_nthreads, 108
.set nb202nf_count, 112
.set nb202nf_mtx, 116
.set nb202nf_outeriter, 120
.set nb202nf_inneriter, 124
.set nb202nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb202nf_ixO, 0
.set nb202nf_iyO, 16
.set nb202nf_izO, 32
.set nb202nf_ixH1, 48
.set nb202nf_iyH1, 64
.set nb202nf_izH1, 80
.set nb202nf_ixH2, 96
.set nb202nf_iyH2, 112
.set nb202nf_izH2, 128
.set nb202nf_jxO, 144
.set nb202nf_jyO, 160
.set nb202nf_jzO, 176
.set nb202nf_jxH1, 192
.set nb202nf_jyH1, 208
.set nb202nf_jzH1, 224
.set nb202nf_jxH2, 240
.set nb202nf_jyH2, 256
.set nb202nf_jzH2, 272
.set nb202nf_qqOO, 288
.set nb202nf_qqOH, 304
.set nb202nf_qqHH, 320
.set nb202nf_vctot, 336
.set nb202nf_half, 352
.set nb202nf_three, 368
.set nb202nf_rsqOO, 384
.set nb202nf_rsqOH1, 400
.set nb202nf_rsqOH2, 416
.set nb202nf_rsqH1O, 432
.set nb202nf_rsqH1H1, 448
.set nb202nf_rsqH1H2, 464
.set nb202nf_rsqH2O, 480
.set nb202nf_rsqH2H1, 496
.set nb202nf_rsqH2H2, 512
.set nb202nf_rinvOO, 528
.set nb202nf_rinvOH1, 544
.set nb202nf_rinvOH2, 560
.set nb202nf_rinvH1O, 576
.set nb202nf_rinvH1H1, 592
.set nb202nf_rinvH1H2, 608
.set nb202nf_rinvH2O, 624
.set nb202nf_rinvH2H1, 640
.set nb202nf_rinvH2H2, 656
.set nb202nf_krf, 672
.set nb202nf_crf, 688
.set nb202nf_is3, 704
.set nb202nf_ii3, 708
.set nb202nf_innerjjnr, 712
.set nb202nf_innerk, 716
.set nb202nf_n, 720
.set nb202nf_nn1, 724
.set nb202nf_nri, 728
.set nb202nf_nouter, 732
.set nb202nf_ninner, 736
.set nb202nf_salign, 740
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $744,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb202nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb202nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb202nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb202nf_nouter(%esp)
        movl %eax,nb202nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb202nf_half(%esp)
        movl %ebx,nb202nf_half+4(%esp)
        movsd nb202nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb202nf_half(%esp)
        movapd %xmm3,nb202nf_three(%esp)

        movl nb202nf_argkrf(%ebp),%esi
        movl nb202nf_argcrf(%ebp),%edi
        movsd (%esi),%xmm5
        movsd (%edi),%xmm6
        shufpd $0,%xmm5,%xmm5
        shufpd $0,%xmm6,%xmm6
        movapd %xmm5,nb202nf_krf(%esp)
        movapd %xmm6,nb202nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb202nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb202nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb202nf_p_facel(%ebp),%esi
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
        movapd %xmm3,nb202nf_qqOO(%esp)
        movapd %xmm4,nb202nf_qqOH(%esp)
        movapd %xmm5,nb202nf_qqHH(%esp)

_nb_kernel202nf_ia32_sse2.nb202nf_threadloop: 
        movl  nb202nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel202nf_ia32_sse2.nb202nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel202nf_ia32_sse2.nb202nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb202nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb202nf_n(%esp)
        movl %ebx,nb202nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel202nf_ia32_sse2.nb202nf_outerstart
        jmp _nb_kernel202nf_ia32_sse2.nb202nf_end

_nb_kernel202nf_ia32_sse2.nb202nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb202nf_nouter(%esp),%ebx
        movl %ebx,nb202nf_nouter(%esp)

_nb_kernel202nf_ia32_sse2.nb202nf_outer: 
        movl  nb202nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 

        movl  nb202nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb202nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx    ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb202nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb202nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb202nf_ixO(%esp)
        movapd %xmm4,nb202nf_iyO(%esp)
        movapd %xmm5,nb202nf_izO(%esp)

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
        movapd %xmm0,nb202nf_ixH1(%esp)
        movapd %xmm1,nb202nf_iyH1(%esp)
        movapd %xmm2,nb202nf_izH1(%esp)
        movapd %xmm3,nb202nf_ixH2(%esp)
        movapd %xmm4,nb202nf_iyH2(%esp)
        movapd %xmm5,nb202nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb202nf_vctot(%esp)

        movl  nb202nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb202nf_pos(%ebp),%esi
        movl  nb202nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb202nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb202nf_ninner(%esp),%ecx
        movl  %ecx,nb202nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb202nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel202nf_ia32_sse2.nb202nf_unroll_loop
        jmp   _nb_kernel202nf_ia32_sse2.nb202nf_checksingle
_nb_kernel202nf_ia32_sse2.nb202nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb202nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb202nf_innerjjnr(%esp)                 ## advance pointer (unrolled 2) 

        movl nb202nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb202nf_jxO(%esp)
        movapd  %xmm3,nb202nf_jyO(%esp)
        movapd  %xmm4,nb202nf_jzO(%esp)
        movapd  %xmm5,nb202nf_jxH1(%esp)
        movapd  %xmm6,nb202nf_jyH1(%esp)
        movapd  %xmm7,nb202nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb202nf_jxH2(%esp)
        movapd  %xmm3,nb202nf_jyH2(%esp)
        movapd  %xmm4,nb202nf_jzH2(%esp)

        movapd nb202nf_ixO(%esp),%xmm0
        movapd nb202nf_iyO(%esp),%xmm1
        movapd nb202nf_izO(%esp),%xmm2
        movapd nb202nf_ixO(%esp),%xmm3
        movapd nb202nf_iyO(%esp),%xmm4
        movapd nb202nf_izO(%esp),%xmm5
        subpd  nb202nf_jxO(%esp),%xmm0
        subpd  nb202nf_jyO(%esp),%xmm1
        subpd  nb202nf_jzO(%esp),%xmm2
        subpd  nb202nf_jxH1(%esp),%xmm3
        subpd  nb202nf_jyH1(%esp),%xmm4
        subpd  nb202nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb202nf_rsqOO(%esp)
        movapd %xmm3,nb202nf_rsqOH1(%esp)

        movapd nb202nf_ixO(%esp),%xmm0
        movapd nb202nf_iyO(%esp),%xmm1
        movapd nb202nf_izO(%esp),%xmm2
        movapd nb202nf_ixH1(%esp),%xmm3
        movapd nb202nf_iyH1(%esp),%xmm4
        movapd nb202nf_izH1(%esp),%xmm5
        subpd  nb202nf_jxH2(%esp),%xmm0
        subpd  nb202nf_jyH2(%esp),%xmm1
        subpd  nb202nf_jzH2(%esp),%xmm2
        subpd  nb202nf_jxO(%esp),%xmm3
        subpd  nb202nf_jyO(%esp),%xmm4
        subpd  nb202nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb202nf_rsqOH2(%esp)
        movapd %xmm3,nb202nf_rsqH1O(%esp)

        movapd nb202nf_ixH1(%esp),%xmm0
        movapd nb202nf_iyH1(%esp),%xmm1
        movapd nb202nf_izH1(%esp),%xmm2
        movapd nb202nf_ixH1(%esp),%xmm3
        movapd nb202nf_iyH1(%esp),%xmm4
        movapd nb202nf_izH1(%esp),%xmm5
        subpd  nb202nf_jxH1(%esp),%xmm0
        subpd  nb202nf_jyH1(%esp),%xmm1
        subpd  nb202nf_jzH1(%esp),%xmm2
        subpd  nb202nf_jxH2(%esp),%xmm3
        subpd  nb202nf_jyH2(%esp),%xmm4
        subpd  nb202nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb202nf_rsqH1H1(%esp)
        movapd %xmm3,nb202nf_rsqH1H2(%esp)

        movapd nb202nf_ixH2(%esp),%xmm0
        movapd nb202nf_iyH2(%esp),%xmm1
        movapd nb202nf_izH2(%esp),%xmm2
        movapd nb202nf_ixH2(%esp),%xmm3
        movapd nb202nf_iyH2(%esp),%xmm4
        movapd nb202nf_izH2(%esp),%xmm5
        subpd  nb202nf_jxO(%esp),%xmm0
        subpd  nb202nf_jyO(%esp),%xmm1
        subpd  nb202nf_jzO(%esp),%xmm2
        subpd  nb202nf_jxH1(%esp),%xmm3
        subpd  nb202nf_jyH1(%esp),%xmm4
        subpd  nb202nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb202nf_rsqH2O(%esp)
        movapd %xmm4,nb202nf_rsqH2H1(%esp)

        movapd nb202nf_ixH2(%esp),%xmm0
        movapd nb202nf_iyH2(%esp),%xmm1
        movapd nb202nf_izH2(%esp),%xmm2
        subpd  nb202nf_jxH2(%esp),%xmm0
        subpd  nb202nf_jyH2(%esp),%xmm1
        subpd  nb202nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb202nf_rsqH2H2(%esp)

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
        movapd  nb202nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb202nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb202nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202nf_rinvH2H2(%esp)
        movapd %xmm5,nb202nf_rinvH2H1(%esp)

        movapd nb202nf_rsqOO(%esp),%xmm0
        movapd nb202nf_rsqOH1(%esp),%xmm4
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
        movapd  nb202nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb202nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb202nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb202nf_rinvOO(%esp)
        movapd %xmm5,nb202nf_rinvOH1(%esp)

        movapd nb202nf_rsqOH2(%esp),%xmm0
        movapd nb202nf_rsqH1O(%esp),%xmm4
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
        movapd  nb202nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb202nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb202nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202nf_rinvOH2(%esp)
        movapd %xmm5,nb202nf_rinvH1O(%esp)

        movapd nb202nf_rsqH1H1(%esp),%xmm0
        movapd nb202nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb202nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb202nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb202nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb202nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202nf_rinvH1H1(%esp)
        movapd %xmm5,nb202nf_rinvH1H2(%esp)

        movapd nb202nf_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb202nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb202nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb202nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb202nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb202nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb202nf_krf(%esp),%xmm6
        mulpd  nb202nf_rsqOO(%esp),%xmm6        ## xmm5=krsq 
        addpd  nb202nf_rinvOO(%esp),%xmm6       ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm6

        mulpd  nb202nf_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  nb202nf_vctot(%esp),%xmm6   ## local vctot summation variable 

        ## O-H1 interaction 
        movapd nb202nf_krf(%esp),%xmm5
        mulpd  nb202nf_rsqOH1(%esp),%xmm5       ## xmm5=krsq 
        addpd  nb202nf_rinvOH1(%esp),%xmm5      ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm5

        mulpd  nb202nf_qqOH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## O-H2 interaction 
        movapd nb202nf_krf(%esp),%xmm7
        mulpd  nb202nf_rsqOH2(%esp),%xmm7       ## xmm5=krsq 
        addpd  nb202nf_rinvOH2(%esp),%xmm7      ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm7

        mulpd  nb202nf_qqOH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H1-O interaction 
        movapd nb202nf_krf(%esp),%xmm4
        mulpd  nb202nf_rsqH1O(%esp),%xmm4       ## xmm5=krsq 
        addpd  nb202nf_rinvH1O(%esp),%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm4

        mulpd  nb202nf_qqOH(%esp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H1-H1 interaction 
        movapd nb202nf_krf(%esp),%xmm5
        mulpd  nb202nf_rsqH1H1(%esp),%xmm5      ## xmm5=krsq 
        addpd  nb202nf_rinvH1H1(%esp),%xmm5     ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm5

        mulpd  nb202nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H1-H2 interaction 
        movapd nb202nf_krf(%esp),%xmm7
        mulpd  nb202nf_rsqH1H2(%esp),%xmm7      ## xmm5=krsq 
        addpd  nb202nf_rinvH1H2(%esp),%xmm7     ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm7

        mulpd  nb202nf_qqHH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H2-O interaction 
        movapd nb202nf_krf(%esp),%xmm4
        mulpd  nb202nf_rsqH2O(%esp),%xmm4       ## xmm5=krsq 
        addpd  nb202nf_rinvH2O(%esp),%xmm4      ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm4

        mulpd  nb202nf_qqOH(%esp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H2-H1 interaction 
        movapd nb202nf_krf(%esp),%xmm5
        mulpd  nb202nf_rsqH2H1(%esp),%xmm5      ## xmm5=krsq 
        addpd  nb202nf_rinvH2H1(%esp),%xmm5     ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm5

        mulpd  nb202nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H2-H2 interaction 
        movapd nb202nf_krf(%esp),%xmm7
        mulpd  nb202nf_rsqH2H2(%esp),%xmm7      ## xmm5=krsq 
        addpd  nb202nf_rinvH2H2(%esp),%xmm7     ## xmm6=rinv+ krsq 
        subpd  nb202nf_crf(%esp),%xmm7

        mulpd  nb202nf_qqHH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addpd  %xmm7,%xmm6 ## local vctot summation variable 
        movapd %xmm6,nb202nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb202nf_innerk(%esp)
        jl    _nb_kernel202nf_ia32_sse2.nb202nf_checksingle
        jmp   _nb_kernel202nf_ia32_sse2.nb202nf_unroll_loop
_nb_kernel202nf_ia32_sse2.nb202nf_checksingle: 
        movl  nb202nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel202nf_ia32_sse2.nb202nf_dosingle
        jmp   _nb_kernel202nf_ia32_sse2.nb202nf_updateouterdata
_nb_kernel202nf_ia32_sse2.nb202nf_dosingle: 
        movl  nb202nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb202nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb202nf_jxO(%esp)
        movapd  %xmm3,nb202nf_jyO(%esp)
        movapd  %xmm4,nb202nf_jzO(%esp)
        movapd  %xmm5,nb202nf_jxH1(%esp)
        movapd  %xmm6,nb202nf_jyH1(%esp)
        movapd  %xmm7,nb202nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb202nf_jxH2(%esp)
        movapd  %xmm3,nb202nf_jyH2(%esp)
        movapd  %xmm4,nb202nf_jzH2(%esp)

        movapd nb202nf_ixO(%esp),%xmm0
        movapd nb202nf_iyO(%esp),%xmm1
        movapd nb202nf_izO(%esp),%xmm2
        movapd nb202nf_ixO(%esp),%xmm3
        movapd nb202nf_iyO(%esp),%xmm4
        movapd nb202nf_izO(%esp),%xmm5
        subsd  nb202nf_jxO(%esp),%xmm0
        subsd  nb202nf_jyO(%esp),%xmm1
        subsd  nb202nf_jzO(%esp),%xmm2
        subsd  nb202nf_jxH1(%esp),%xmm3
        subsd  nb202nf_jyH1(%esp),%xmm4
        subsd  nb202nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb202nf_rsqOO(%esp)
        movapd %xmm3,nb202nf_rsqOH1(%esp)

        movapd nb202nf_ixO(%esp),%xmm0
        movapd nb202nf_iyO(%esp),%xmm1
        movapd nb202nf_izO(%esp),%xmm2
        movapd nb202nf_ixH1(%esp),%xmm3
        movapd nb202nf_iyH1(%esp),%xmm4
        movapd nb202nf_izH1(%esp),%xmm5
        subsd  nb202nf_jxH2(%esp),%xmm0
        subsd  nb202nf_jyH2(%esp),%xmm1
        subsd  nb202nf_jzH2(%esp),%xmm2
        subsd  nb202nf_jxO(%esp),%xmm3
        subsd  nb202nf_jyO(%esp),%xmm4
        subsd  nb202nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb202nf_rsqOH2(%esp)
        movapd %xmm3,nb202nf_rsqH1O(%esp)

        movapd nb202nf_ixH1(%esp),%xmm0
        movapd nb202nf_iyH1(%esp),%xmm1
        movapd nb202nf_izH1(%esp),%xmm2
        movapd nb202nf_ixH1(%esp),%xmm3
        movapd nb202nf_iyH1(%esp),%xmm4
        movapd nb202nf_izH1(%esp),%xmm5
        subsd  nb202nf_jxH1(%esp),%xmm0
        subsd  nb202nf_jyH1(%esp),%xmm1
        subsd  nb202nf_jzH1(%esp),%xmm2
        subsd  nb202nf_jxH2(%esp),%xmm3
        subsd  nb202nf_jyH2(%esp),%xmm4
        subsd  nb202nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb202nf_rsqH1H1(%esp)
        movapd %xmm3,nb202nf_rsqH1H2(%esp)

        movapd nb202nf_ixH2(%esp),%xmm0
        movapd nb202nf_iyH2(%esp),%xmm1
        movapd nb202nf_izH2(%esp),%xmm2
        movapd nb202nf_ixH2(%esp),%xmm3
        movapd nb202nf_iyH2(%esp),%xmm4
        movapd nb202nf_izH2(%esp),%xmm5
        subsd  nb202nf_jxO(%esp),%xmm0
        subsd  nb202nf_jyO(%esp),%xmm1
        subsd  nb202nf_jzO(%esp),%xmm2
        subsd  nb202nf_jxH1(%esp),%xmm3
        subsd  nb202nf_jyH1(%esp),%xmm4
        subsd  nb202nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb202nf_rsqH2O(%esp)
        movapd %xmm4,nb202nf_rsqH2H1(%esp)

        movapd nb202nf_ixH2(%esp),%xmm0
        movapd nb202nf_iyH2(%esp),%xmm1
        movapd nb202nf_izH2(%esp),%xmm2
        subsd  nb202nf_jxH2(%esp),%xmm0
        subsd  nb202nf_jyH2(%esp),%xmm1
        subsd  nb202nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb202nf_rsqH2H2(%esp)

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
        movapd  nb202nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb202nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb202nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202nf_rinvH2H2(%esp)
        movapd %xmm5,nb202nf_rinvH2H1(%esp)

        movapd nb202nf_rsqOO(%esp),%xmm0
        movapd nb202nf_rsqOH1(%esp),%xmm4
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
        movapd  nb202nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb202nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb202nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb202nf_rinvOO(%esp)
        movapd %xmm5,nb202nf_rinvOH1(%esp)

        movapd nb202nf_rsqOH2(%esp),%xmm0
        movapd nb202nf_rsqH1O(%esp),%xmm4
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
        movapd  nb202nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb202nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb202nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202nf_rinvOH2(%esp)
        movapd %xmm5,nb202nf_rinvH1O(%esp)

        movapd nb202nf_rsqH1H1(%esp),%xmm0
        movapd nb202nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb202nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb202nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb202nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb202nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb202nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb202nf_rinvH1H1(%esp)
        movapd %xmm5,nb202nf_rinvH1H2(%esp)

        movapd nb202nf_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb202nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb202nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb202nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb202nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb202nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb202nf_krf(%esp),%xmm6
        mulsd  nb202nf_rsqOO(%esp),%xmm6        ## xmm5=krsq 
        addsd  nb202nf_rinvOO(%esp),%xmm6       ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm6

        mulsd  nb202nf_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  nb202nf_vctot(%esp),%xmm6   ## local vctot summation variable 

        ## O-H1 interaction 
        movapd nb202nf_krf(%esp),%xmm5
        mulsd  nb202nf_rsqOH1(%esp),%xmm5       ## xmm5=krsq 
        addsd  nb202nf_rinvOH1(%esp),%xmm5      ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm5

        mulsd  nb202nf_qqOH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## O-H2 interaction 
        movapd nb202nf_krf(%esp),%xmm7
        mulsd  nb202nf_rsqOH2(%esp),%xmm7       ## xmm5=krsq 
        addsd  nb202nf_rinvOH2(%esp),%xmm7      ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm7

        mulsd  nb202nf_qqOH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H1-O interaction 
        movapd nb202nf_krf(%esp),%xmm4
        mulsd  nb202nf_rsqH1O(%esp),%xmm4       ## xmm5=krsq 
        addsd  nb202nf_rinvH1O(%esp),%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm4

        mulsd  nb202nf_qqOH(%esp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H1-H1 interaction 
        movapd nb202nf_krf(%esp),%xmm5
        mulsd  nb202nf_rsqH1H1(%esp),%xmm5      ## xmm5=krsq 
        addsd  nb202nf_rinvH1H1(%esp),%xmm5     ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm5

        mulsd  nb202nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H1-H2 interaction 
        movapd nb202nf_krf(%esp),%xmm7
        mulsd  nb202nf_rsqH1H2(%esp),%xmm7      ## xmm5=krsq 
        addsd  nb202nf_rinvH1H2(%esp),%xmm7     ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm7

        mulsd  nb202nf_qqHH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 

        ## H2-O interaction 
        movapd nb202nf_krf(%esp),%xmm4
        mulsd  nb202nf_rsqH2O(%esp),%xmm4       ## xmm5=krsq 
        addsd  nb202nf_rinvH2O(%esp),%xmm4      ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm4

        mulsd  nb202nf_qqOH(%esp),%xmm4   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm4,%xmm6 ## local vctot summation variable 

        ## H2-H1 interaction 
        movapd nb202nf_krf(%esp),%xmm5
        mulsd  nb202nf_rsqH2H1(%esp),%xmm5      ## xmm5=krsq 
        addsd  nb202nf_rinvH2H1(%esp),%xmm5     ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm5

        mulsd  nb202nf_qqHH(%esp),%xmm5   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm5,%xmm6 ## local vctot summation variable 

        ## H2-H2 interaction 
        movapd nb202nf_krf(%esp),%xmm7
        mulsd  nb202nf_rsqH2H2(%esp),%xmm7      ## xmm5=krsq 
        addsd  nb202nf_rinvH2H2(%esp),%xmm7     ## xmm6=rinv+ krsq 
        subsd  nb202nf_crf(%esp),%xmm7

        mulsd  nb202nf_qqHH(%esp),%xmm7   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addsd  %xmm7,%xmm6 ## local vctot summation variable 
        movlpd %xmm6,nb202nf_vctot(%esp)

_nb_kernel202nf_ia32_sse2.nb202nf_updateouterdata: 
        ## get n from stack
        movl nb202nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb202nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb202nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb202nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb202nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel202nf_ia32_sse2.nb202nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb202nf_n(%esp)
        jmp _nb_kernel202nf_ia32_sse2.nb202nf_outer
_nb_kernel202nf_ia32_sse2.nb202nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb202nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel202nf_ia32_sse2.nb202nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel202nf_ia32_sse2.nb202nf_threadloop
_nb_kernel202nf_ia32_sse2.nb202nf_end: 
        emms

        movl nb202nf_nouter(%esp),%eax
        movl nb202nf_ninner(%esp),%ebx
        movl nb202nf_outeriter(%ebp),%ecx
        movl nb202nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb202nf_salign(%esp),%eax
        addl %eax,%esp
        addl $744,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



