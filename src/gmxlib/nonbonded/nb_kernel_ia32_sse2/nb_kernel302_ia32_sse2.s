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


.globl nb_kernel302_ia32_sse2
.globl _nb_kernel302_ia32_sse2
nb_kernel302_ia32_sse2: 
_nb_kernel302_ia32_sse2:        
.set nb302_p_nri, 8
.set nb302_iinr, 12
.set nb302_jindex, 16
.set nb302_jjnr, 20
.set nb302_shift, 24
.set nb302_shiftvec, 28
.set nb302_fshift, 32
.set nb302_gid, 36
.set nb302_pos, 40
.set nb302_faction, 44
.set nb302_charge, 48
.set nb302_p_facel, 52
.set nb302_argkrf, 56
.set nb302_argcrf, 60
.set nb302_Vc, 64
.set nb302_type, 68
.set nb302_p_ntype, 72
.set nb302_vdwparam, 76
.set nb302_Vvdw, 80
.set nb302_p_tabscale, 84
.set nb302_VFtab, 88
.set nb302_invsqrta, 92
.set nb302_dvda, 96
.set nb302_p_gbtabscale, 100
.set nb302_GBtab, 104
.set nb302_p_nthreads, 108
.set nb302_count, 112
.set nb302_mtx, 116
.set nb302_outeriter, 120
.set nb302_inneriter, 124
.set nb302_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse2 use 
.set nb302_ixO, 0
.set nb302_iyO, 16
.set nb302_izO, 32
.set nb302_ixH1, 48
.set nb302_iyH1, 64
.set nb302_izH1, 80
.set nb302_ixH2, 96
.set nb302_iyH2, 112
.set nb302_izH2, 128
.set nb302_jxO, 144
.set nb302_jyO, 160
.set nb302_jzO, 176
.set nb302_jxH1, 192
.set nb302_jyH1, 208
.set nb302_jzH1, 224
.set nb302_jxH2, 240
.set nb302_jyH2, 256
.set nb302_jzH2, 272
.set nb302_dxOO, 288
.set nb302_dyOO, 304
.set nb302_dzOO, 320
.set nb302_dxOH1, 336
.set nb302_dyOH1, 352
.set nb302_dzOH1, 368
.set nb302_dxOH2, 384
.set nb302_dyOH2, 400
.set nb302_dzOH2, 416
.set nb302_dxH1O, 432
.set nb302_dyH1O, 448
.set nb302_dzH1O, 464
.set nb302_dxH1H1, 480
.set nb302_dyH1H1, 496
.set nb302_dzH1H1, 512
.set nb302_dxH1H2, 528
.set nb302_dyH1H2, 544
.set nb302_dzH1H2, 560
.set nb302_dxH2O, 576
.set nb302_dyH2O, 592
.set nb302_dzH2O, 608
.set nb302_dxH2H1, 624
.set nb302_dyH2H1, 640
.set nb302_dzH2H1, 656
.set nb302_dxH2H2, 672
.set nb302_dyH2H2, 688
.set nb302_dzH2H2, 704
.set nb302_qqOO, 720
.set nb302_qqOH, 736
.set nb302_qqHH, 752
.set nb302_two, 768
.set nb302_tsc, 784
.set nb302_vctot, 800
.set nb302_fixO, 816
.set nb302_fiyO, 832
.set nb302_fizO, 848
.set nb302_fixH1, 864
.set nb302_fiyH1, 880
.set nb302_fizH1, 896
.set nb302_fixH2, 912
.set nb302_fiyH2, 928
.set nb302_fizH2, 944
.set nb302_fjxO, 960
.set nb302_fjyO, 976
.set nb302_fjzO, 992
.set nb302_fjxH1, 1008
.set nb302_fjyH1, 1024
.set nb302_fjzH1, 1040
.set nb302_fjxH2, 1056
.set nb302_fjyH2, 1072
.set nb302_fjzH2, 1088
.set nb302_half, 1104
.set nb302_three, 1120
.set nb302_rsqOO, 1136
.set nb302_rsqOH1, 1152
.set nb302_rsqOH2, 1168
.set nb302_rsqH1O, 1184
.set nb302_rsqH1H1, 1200
.set nb302_rsqH1H2, 1216
.set nb302_rsqH2O, 1232
.set nb302_rsqH2H1, 1248
.set nb302_rsqH2H2, 1264
.set nb302_rinvOO, 1280
.set nb302_rinvOH1, 1296
.set nb302_rinvOH2, 1312
.set nb302_rinvH1O, 1328
.set nb302_rinvH1H1, 1344
.set nb302_rinvH1H2, 1360
.set nb302_rinvH2O, 1376
.set nb302_rinvH2H1, 1392
.set nb302_rinvH2H2, 1408
.set nb302_is3, 1424
.set nb302_ii3, 1428
.set nb302_innerjjnr, 1432
.set nb302_innerk, 1436
.set nb302_n, 1440
.set nb302_nn1, 1444
.set nb302_nri, 1448
.set nb302_nouter, 1452
.set nb302_ninner, 1456
.set nb302_salign, 1460
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $1464,%esp         ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb302_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb302_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb302_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb302_nouter(%esp)
        movl %eax,nb302_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb302_half(%esp)
        movl %ebx,nb302_half+4(%esp)
        movsd nb302_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb302_half(%esp)
        movapd %xmm2,nb302_two(%esp)
        movapd %xmm3,nb302_three(%esp)

        movl nb302_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb302_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb302_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb302_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb302_p_facel(%ebp),%esi
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
        movapd %xmm3,nb302_qqOO(%esp)
        movapd %xmm4,nb302_qqOH(%esp)
        movapd %xmm5,nb302_qqHH(%esp)

_nb_kernel302_ia32_sse2.nb302_threadloop: 
        movl  nb302_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel302_ia32_sse2.nb302_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel302_ia32_sse2.nb302_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb302_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb302_n(%esp)
        movl %ebx,nb302_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel302_ia32_sse2.nb302_outerstart
        jmp _nb_kernel302_ia32_sse2.nb302_end

_nb_kernel302_ia32_sse2.nb302_outerstart: 
        ## ebx contains number of outer iterations
        addl nb302_nouter(%esp),%ebx
        movl %ebx,nb302_nouter(%esp)

_nb_kernel302_ia32_sse2.nb302_outer: 
        movl  nb302_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb302_is3(%esp)      ## store is3 

        movl  nb302_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb302_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb302_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb302_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb302_ixO(%esp)
        movapd %xmm4,nb302_iyO(%esp)
        movapd %xmm5,nb302_izO(%esp)

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
        movapd %xmm0,nb302_ixH1(%esp)
        movapd %xmm1,nb302_iyH1(%esp)
        movapd %xmm2,nb302_izH1(%esp)
        movapd %xmm3,nb302_ixH2(%esp)
        movapd %xmm4,nb302_iyH2(%esp)
        movapd %xmm5,nb302_izH2(%esp)

        ## clear vctot and i forces 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb302_vctot(%esp)
        movapd %xmm4,nb302_fixO(%esp)
        movapd %xmm4,nb302_fiyO(%esp)
        movapd %xmm4,nb302_fizO(%esp)
        movapd %xmm4,nb302_fixH1(%esp)
        movapd %xmm4,nb302_fiyH1(%esp)
        movapd %xmm4,nb302_fizH1(%esp)
        movapd %xmm4,nb302_fixH2(%esp)
        movapd %xmm4,nb302_fiyH2(%esp)
        movapd %xmm4,nb302_fizH2(%esp)

        movl  nb302_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb302_pos(%ebp),%esi
        movl  nb302_faction(%ebp),%edi
        movl  nb302_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb302_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb302_ninner(%esp),%ecx
        movl  %ecx,nb302_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb302_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel302_ia32_sse2.nb302_unroll_loop
        jmp   _nb_kernel302_ia32_sse2.nb302_checksingle
_nb_kernel302_ia32_sse2.nb302_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb302_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb302_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb302_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb302_jxO(%esp)
        movapd  %xmm3,nb302_jyO(%esp)
        movapd  %xmm4,nb302_jzO(%esp)
        movapd  %xmm5,nb302_jxH1(%esp)
        movapd  %xmm6,nb302_jyH1(%esp)
        movapd  %xmm7,nb302_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb302_jxH2(%esp)
        movapd  %xmm3,nb302_jyH2(%esp)
        movapd  %xmm4,nb302_jzH2(%esp)

        movapd nb302_ixO(%esp),%xmm0
        movapd nb302_iyO(%esp),%xmm1
        movapd nb302_izO(%esp),%xmm2
        movapd nb302_ixO(%esp),%xmm3
        movapd nb302_iyO(%esp),%xmm4
        movapd nb302_izO(%esp),%xmm5
        subpd  nb302_jxO(%esp),%xmm0
        subpd  nb302_jyO(%esp),%xmm1
        subpd  nb302_jzO(%esp),%xmm2
        subpd  nb302_jxH1(%esp),%xmm3
        subpd  nb302_jyH1(%esp),%xmm4
        subpd  nb302_jzH1(%esp),%xmm5
        movapd %xmm0,nb302_dxOO(%esp)
        movapd %xmm1,nb302_dyOO(%esp)
        movapd %xmm2,nb302_dzOO(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxOH1(%esp)
        movapd %xmm4,nb302_dyOH1(%esp)
        movapd %xmm5,nb302_dzOH1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb302_rsqOO(%esp)
        movapd %xmm3,nb302_rsqOH1(%esp)

        movapd nb302_ixO(%esp),%xmm0
        movapd nb302_iyO(%esp),%xmm1
        movapd nb302_izO(%esp),%xmm2
        movapd nb302_ixH1(%esp),%xmm3
        movapd nb302_iyH1(%esp),%xmm4
        movapd nb302_izH1(%esp),%xmm5
        subpd  nb302_jxH2(%esp),%xmm0
        subpd  nb302_jyH2(%esp),%xmm1
        subpd  nb302_jzH2(%esp),%xmm2
        subpd  nb302_jxO(%esp),%xmm3
        subpd  nb302_jyO(%esp),%xmm4
        subpd  nb302_jzO(%esp),%xmm5
        movapd %xmm0,nb302_dxOH2(%esp)
        movapd %xmm1,nb302_dyOH2(%esp)
        movapd %xmm2,nb302_dzOH2(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH1O(%esp)
        movapd %xmm4,nb302_dyH1O(%esp)
        movapd %xmm5,nb302_dzH1O(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb302_rsqOH2(%esp)
        movapd %xmm3,nb302_rsqH1O(%esp)

        movapd nb302_ixH1(%esp),%xmm0
        movapd nb302_iyH1(%esp),%xmm1
        movapd nb302_izH1(%esp),%xmm2
        movapd nb302_ixH1(%esp),%xmm3
        movapd nb302_iyH1(%esp),%xmm4
        movapd nb302_izH1(%esp),%xmm5
        subpd  nb302_jxH1(%esp),%xmm0
        subpd  nb302_jyH1(%esp),%xmm1
        subpd  nb302_jzH1(%esp),%xmm2
        subpd  nb302_jxH2(%esp),%xmm3
        subpd  nb302_jyH2(%esp),%xmm4
        subpd  nb302_jzH2(%esp),%xmm5
        movapd %xmm0,nb302_dxH1H1(%esp)
        movapd %xmm1,nb302_dyH1H1(%esp)
        movapd %xmm2,nb302_dzH1H1(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH1H2(%esp)
        movapd %xmm4,nb302_dyH1H2(%esp)
        movapd %xmm5,nb302_dzH1H2(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm4,%xmm3
        addpd  %xmm5,%xmm3
        movapd %xmm0,nb302_rsqH1H1(%esp)
        movapd %xmm3,nb302_rsqH1H2(%esp)

        movapd nb302_ixH2(%esp),%xmm0
        movapd nb302_iyH2(%esp),%xmm1
        movapd nb302_izH2(%esp),%xmm2
        movapd nb302_ixH2(%esp),%xmm3
        movapd nb302_iyH2(%esp),%xmm4
        movapd nb302_izH2(%esp),%xmm5
        subpd  nb302_jxO(%esp),%xmm0
        subpd  nb302_jyO(%esp),%xmm1
        subpd  nb302_jzO(%esp),%xmm2
        subpd  nb302_jxH1(%esp),%xmm3
        subpd  nb302_jyH1(%esp),%xmm4
        subpd  nb302_jzH1(%esp),%xmm5
        movapd %xmm0,nb302_dxH2O(%esp)
        movapd %xmm1,nb302_dyH2O(%esp)
        movapd %xmm2,nb302_dzH2O(%esp)
        mulpd  %xmm0,%xmm0
        mulpd  %xmm1,%xmm1
        mulpd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH2H1(%esp)
        movapd %xmm4,nb302_dyH2H1(%esp)
        movapd %xmm5,nb302_dzH2H1(%esp)
        mulpd  %xmm3,%xmm3
        mulpd  %xmm4,%xmm4
        mulpd  %xmm5,%xmm5
        addpd  %xmm1,%xmm0
        addpd  %xmm2,%xmm0
        addpd  %xmm3,%xmm4
        addpd  %xmm5,%xmm4
        movapd %xmm0,nb302_rsqH2O(%esp)
        movapd %xmm4,nb302_rsqH2H1(%esp)

        movapd nb302_ixH2(%esp),%xmm0
        movapd nb302_iyH2(%esp),%xmm1
        movapd nb302_izH2(%esp),%xmm2
        subpd  nb302_jxH2(%esp),%xmm0
        subpd  nb302_jyH2(%esp),%xmm1
        subpd  nb302_jzH2(%esp),%xmm2
        movapd %xmm0,nb302_dxH2H2(%esp)
        movapd %xmm1,nb302_dyH2H2(%esp)
        movapd %xmm2,nb302_dzH2H2(%esp)
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb302_rsqH2H2(%esp)

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
        movapd  nb302_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302_half(%esp),%xmm3   ## iter1 
        mulpd   nb302_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302_half(%esp),%xmm1   ## rinv 
        mulpd   nb302_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302_rinvH2H2(%esp)
        movapd %xmm5,nb302_rinvH2H1(%esp)

        movapd nb302_rsqOO(%esp),%xmm0
        movapd nb302_rsqOH1(%esp),%xmm4
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
        movapd  nb302_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb302_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302_half(%esp),%xmm1   ## rinv 
        mulpd   nb302_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb302_rinvOO(%esp)
        movapd %xmm5,nb302_rinvOH1(%esp)

        movapd nb302_rsqOH2(%esp),%xmm0
        movapd nb302_rsqH1O(%esp),%xmm4
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
        movapd  nb302_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302_half(%esp),%xmm3   ## iter1 
        mulpd   nb302_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302_half(%esp),%xmm1   ## rinv 
        mulpd   nb302_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302_rinvOH2(%esp)
        movapd %xmm5,nb302_rinvH1O(%esp)

        movapd nb302_rsqH1H1(%esp),%xmm0
        movapd nb302_rsqH1H2(%esp),%xmm4
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
        movapd  nb302_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302_half(%esp),%xmm3   ## iter1a 
        mulpd   nb302_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302_half(%esp),%xmm1   ## rinv 
        mulpd   nb302_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302_rinvH1H1(%esp)
        movapd %xmm5,nb302_rinvH1H2(%esp)

        movapd nb302_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb302_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb302_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb302_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb302_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb302_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb302_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movd %eax,%mm0
        movd %ebx,%mm1
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOO(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addpd  nb302_vctot(%esp),%xmm5
        xorpd  %xmm2,%xmm2
    movapd %xmm5,nb302_vctot(%esp)
        mulpd  nb302_tsc(%esp),%xmm3

        subpd  %xmm3,%xmm2
        mulpd  %xmm2,%xmm0      ## mult by rinv 

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb302_dxOO(%esp),%xmm0
        mulpd nb302_dyOO(%esp),%xmm1
        mulpd nb302_dzOO(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixO(%esp),%xmm0
        addpd nb302_fiyO(%esp),%xmm1
        addpd nb302_fizO(%esp),%xmm2
        movapd %xmm3,nb302_fjxO(%esp)
        movapd %xmm4,nb302_fjyO(%esp)
        movapd %xmm5,nb302_fjzO(%esp)
        movapd %xmm0,nb302_fixO(%esp)
        movapd %xmm1,nb302_fiyO(%esp)
        movapd %xmm2,nb302_fizO(%esp)

        ## O-H1 interaction 
        movapd nb302_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqOH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb302_vctot(%esp),%xmm5
    movapd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb302_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb302_dxOH1(%esp),%xmm0
        mulpd nb302_dyOH1(%esp),%xmm1
        mulpd nb302_dzOH1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixO(%esp),%xmm0
        addpd nb302_fiyO(%esp),%xmm1
        addpd nb302_fizO(%esp),%xmm2
        movapd %xmm3,nb302_fjxH1(%esp)
        movapd %xmm4,nb302_fjyH1(%esp)
        movapd %xmm5,nb302_fjzH1(%esp)
        movapd %xmm0,nb302_fixO(%esp)
        movapd %xmm1,nb302_fiyO(%esp)
        movapd %xmm2,nb302_fizO(%esp)

        ## O-H2 interaction  
        movapd nb302_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqOH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb302_vctot(%esp),%xmm5
    movapd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb302_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulpd nb302_dxOH2(%esp),%xmm0
        mulpd nb302_dyOH2(%esp),%xmm1
        mulpd nb302_dzOH2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixO(%esp),%xmm0
        addpd nb302_fiyO(%esp),%xmm1
        addpd nb302_fizO(%esp),%xmm2
        movapd %xmm3,nb302_fjxH2(%esp)
        movapd %xmm4,nb302_fjyH2(%esp)
        movapd %xmm5,nb302_fjzH2(%esp)
        movapd %xmm0,nb302_fixO(%esp)
        movapd %xmm1,nb302_fiyO(%esp)
        movapd %xmm2,nb302_fizO(%esp)

        ## H1-O interaction 
        movapd nb302_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqH1O(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb302_vctot(%esp),%xmm5
    movapd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb302_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxO(%esp),%xmm3
        movapd nb302_fjyO(%esp),%xmm4
        movapd nb302_fjzO(%esp),%xmm5
        mulpd nb302_dxH1O(%esp),%xmm0
        mulpd nb302_dyH1O(%esp),%xmm1
        mulpd nb302_dzH1O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixH1(%esp),%xmm0
        addpd nb302_fiyH1(%esp),%xmm1
        addpd nb302_fizH1(%esp),%xmm2
        movapd %xmm3,nb302_fjxO(%esp)
        movapd %xmm4,nb302_fjyO(%esp)
        movapd %xmm5,nb302_fjzO(%esp)
        movapd %xmm0,nb302_fixH1(%esp)
        movapd %xmm1,nb302_fiyH1(%esp)
        movapd %xmm2,nb302_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb302_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb302_vctot(%esp),%xmm5
    movapd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb302_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxH1(%esp),%xmm3
        movapd nb302_fjyH1(%esp),%xmm4
        movapd nb302_fjzH1(%esp),%xmm5
        mulpd nb302_dxH1H1(%esp),%xmm0
        mulpd nb302_dyH1H1(%esp),%xmm1
        mulpd nb302_dzH1H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixH1(%esp),%xmm0
        addpd nb302_fiyH1(%esp),%xmm1
        addpd nb302_fizH1(%esp),%xmm2
        movapd %xmm3,nb302_fjxH1(%esp)
        movapd %xmm4,nb302_fjyH1(%esp)
        movapd %xmm5,nb302_fjzH1(%esp)
        movapd %xmm0,nb302_fixH1(%esp)
        movapd %xmm1,nb302_fiyH1(%esp)
        movapd %xmm2,nb302_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb302_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb302_vctot(%esp),%xmm5
    movapd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb302_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxH2(%esp),%xmm3
        movapd nb302_fjyH2(%esp),%xmm4
        movapd nb302_fjzH2(%esp),%xmm5
        mulpd nb302_dxH1H2(%esp),%xmm0
        mulpd nb302_dyH1H2(%esp),%xmm1
        mulpd nb302_dzH1H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixH1(%esp),%xmm0
        addpd nb302_fiyH1(%esp),%xmm1
        addpd nb302_fizH1(%esp),%xmm2
        movapd %xmm3,nb302_fjxH2(%esp)
        movapd %xmm4,nb302_fjyH2(%esp)
        movapd %xmm5,nb302_fjzH2(%esp)
        movapd %xmm0,nb302_fixH1(%esp)
        movapd %xmm1,nb302_fiyH1(%esp)
        movapd %xmm2,nb302_fizH1(%esp)

        ## H2-O interaction 
        movapd nb302_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqH2O(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb302_vctot(%esp),%xmm5
    movapd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb302_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxO(%esp),%xmm3
        movapd nb302_fjyO(%esp),%xmm4
        movapd nb302_fjzO(%esp),%xmm5
        mulpd nb302_dxH2O(%esp),%xmm0
        mulpd nb302_dyH2O(%esp),%xmm1
        mulpd nb302_dzH2O(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixH2(%esp),%xmm0
        addpd nb302_fiyH2(%esp),%xmm1
        addpd nb302_fizH2(%esp),%xmm2
        movapd %xmm3,nb302_fjxO(%esp)
        movapd %xmm4,nb302_fjyO(%esp)
        movapd %xmm5,nb302_fjzO(%esp)
        movapd %xmm0,nb302_fixH2(%esp)
        movapd %xmm1,nb302_fiyH2(%esp)
        movapd %xmm2,nb302_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb302_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb302_vctot(%esp),%xmm5
    movapd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb302_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxH1(%esp),%xmm3
        movapd nb302_fjyH1(%esp),%xmm4
        movapd nb302_fjzH1(%esp),%xmm5
        mulpd nb302_dxH2H1(%esp),%xmm0
        mulpd nb302_dyH2H1(%esp),%xmm1
        mulpd nb302_dzH2H1(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixH2(%esp),%xmm0
        addpd nb302_fiyH2(%esp),%xmm1
        addpd nb302_fizH2(%esp),%xmm2
        movapd %xmm3,nb302_fjxH1(%esp)
        movapd %xmm4,nb302_fjyH1(%esp)
        movapd %xmm5,nb302_fjzH1(%esp)
        movapd %xmm0,nb302_fixH2(%esp)
        movapd %xmm1,nb302_fiyH2(%esp)
        movapd %xmm2,nb302_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb302_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi
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
        mulpd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqHH(%esp),%xmm3
        addpd  %xmm6,%xmm7
        addpd  %xmm5,%xmm7 ## xmm7=FF 
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulpd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addpd  nb302_vctot(%esp),%xmm5
    movapd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulpd  nb302_tsc(%esp),%xmm3
        mulpd  %xmm0,%xmm3
        subpd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxH2(%esp),%xmm3
        movapd nb302_fjyH2(%esp),%xmm4
        movapd nb302_fjzH2(%esp),%xmm5
        mulpd nb302_dxH2H2(%esp),%xmm0
        mulpd nb302_dyH2H2(%esp),%xmm1
        mulpd nb302_dzH2H2(%esp),%xmm2
        subpd %xmm0,%xmm3
        subpd %xmm1,%xmm4
        subpd %xmm2,%xmm5
        addpd nb302_fixH2(%esp),%xmm0
        addpd nb302_fiyH2(%esp),%xmm1
        addpd nb302_fizH2(%esp),%xmm2
        movapd %xmm3,nb302_fjxH2(%esp)
        movapd %xmm4,nb302_fjyH2(%esp)
        movapd %xmm5,nb302_fjzH2(%esp)
        movapd %xmm0,nb302_fixH2(%esp)
        movapd %xmm1,nb302_fiyH2(%esp)
        movapd %xmm2,nb302_fizH2(%esp)

        movl nb302_faction(%ebp),%edi

        movd %mm0,%eax
        movd %mm1,%ebx

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
        addpd nb302_fjxO(%esp),%xmm0
        addpd nb302_fjyO(%esp),%xmm1
        addpd nb302_fjzO(%esp),%xmm2
        addpd nb302_fjxH1(%esp),%xmm3
        addpd nb302_fjyH1(%esp),%xmm4
        addpd nb302_fjzH1(%esp),%xmm5
        addpd nb302_fjxH2(%esp),%xmm6
        addpd nb302_fjyH2(%esp),%xmm7
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
        addpd nb302_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)
        movhpd %xmm0,64(%edi,%ebx,8)

        ## should we do one more iteration? 
        subl $2,nb302_innerk(%esp)
        jl    _nb_kernel302_ia32_sse2.nb302_checksingle
        jmp   _nb_kernel302_ia32_sse2.nb302_unroll_loop
_nb_kernel302_ia32_sse2.nb302_checksingle: 
        movl  nb302_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel302_ia32_sse2.nb302_dosingle
        jmp   _nb_kernel302_ia32_sse2.nb302_updateouterdata
_nb_kernel302_ia32_sse2.nb302_dosingle: 
        movl  nb302_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb302_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb302_jxO(%esp)
        movapd  %xmm3,nb302_jyO(%esp)
        movapd  %xmm4,nb302_jzO(%esp)
        movapd  %xmm5,nb302_jxH1(%esp)
        movapd  %xmm6,nb302_jyH1(%esp)
        movapd  %xmm7,nb302_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb302_jxH2(%esp)
        movapd  %xmm3,nb302_jyH2(%esp)
        movapd  %xmm4,nb302_jzH2(%esp)

        movapd nb302_ixO(%esp),%xmm0
        movapd nb302_iyO(%esp),%xmm1
        movapd nb302_izO(%esp),%xmm2
        movapd nb302_ixO(%esp),%xmm3
        movapd nb302_iyO(%esp),%xmm4
        movapd nb302_izO(%esp),%xmm5
        subsd  nb302_jxO(%esp),%xmm0
        subsd  nb302_jyO(%esp),%xmm1
        subsd  nb302_jzO(%esp),%xmm2
        subsd  nb302_jxH1(%esp),%xmm3
        subsd  nb302_jyH1(%esp),%xmm4
        subsd  nb302_jzH1(%esp),%xmm5
        movapd %xmm0,nb302_dxOO(%esp)
        movapd %xmm1,nb302_dyOO(%esp)
        movapd %xmm2,nb302_dzOO(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxOH1(%esp)
        movapd %xmm4,nb302_dyOH1(%esp)
        movapd %xmm5,nb302_dzOH1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb302_rsqOO(%esp)
        movapd %xmm3,nb302_rsqOH1(%esp)

        movapd nb302_ixO(%esp),%xmm0
        movapd nb302_iyO(%esp),%xmm1
        movapd nb302_izO(%esp),%xmm2
        movapd nb302_ixH1(%esp),%xmm3
        movapd nb302_iyH1(%esp),%xmm4
        movapd nb302_izH1(%esp),%xmm5
        subsd  nb302_jxH2(%esp),%xmm0
        subsd  nb302_jyH2(%esp),%xmm1
        subsd  nb302_jzH2(%esp),%xmm2
        subsd  nb302_jxO(%esp),%xmm3
        subsd  nb302_jyO(%esp),%xmm4
        subsd  nb302_jzO(%esp),%xmm5
        movapd %xmm0,nb302_dxOH2(%esp)
        movapd %xmm1,nb302_dyOH2(%esp)
        movapd %xmm2,nb302_dzOH2(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH1O(%esp)
        movapd %xmm4,nb302_dyH1O(%esp)
        movapd %xmm5,nb302_dzH1O(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb302_rsqOH2(%esp)
        movapd %xmm3,nb302_rsqH1O(%esp)

        movapd nb302_ixH1(%esp),%xmm0
        movapd nb302_iyH1(%esp),%xmm1
        movapd nb302_izH1(%esp),%xmm2
        movapd nb302_ixH1(%esp),%xmm3
        movapd nb302_iyH1(%esp),%xmm4
        movapd nb302_izH1(%esp),%xmm5
        subsd  nb302_jxH1(%esp),%xmm0
        subsd  nb302_jyH1(%esp),%xmm1
        subsd  nb302_jzH1(%esp),%xmm2
        subsd  nb302_jxH2(%esp),%xmm3
        subsd  nb302_jyH2(%esp),%xmm4
        subsd  nb302_jzH2(%esp),%xmm5
        movapd %xmm0,nb302_dxH1H1(%esp)
        movapd %xmm1,nb302_dyH1H1(%esp)
        movapd %xmm2,nb302_dzH1H1(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH1H2(%esp)
        movapd %xmm4,nb302_dyH1H2(%esp)
        movapd %xmm5,nb302_dzH1H2(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm4,%xmm3
        addsd  %xmm5,%xmm3
        movapd %xmm0,nb302_rsqH1H1(%esp)
        movapd %xmm3,nb302_rsqH1H2(%esp)

        movapd nb302_ixH2(%esp),%xmm0
        movapd nb302_iyH2(%esp),%xmm1
        movapd nb302_izH2(%esp),%xmm2
        movapd nb302_ixH2(%esp),%xmm3
        movapd nb302_iyH2(%esp),%xmm4
        movapd nb302_izH2(%esp),%xmm5
        subsd  nb302_jxO(%esp),%xmm0
        subsd  nb302_jyO(%esp),%xmm1
        subsd  nb302_jzO(%esp),%xmm2
        subsd  nb302_jxH1(%esp),%xmm3
        subsd  nb302_jyH1(%esp),%xmm4
        subsd  nb302_jzH1(%esp),%xmm5
        movapd %xmm0,nb302_dxH2O(%esp)
        movapd %xmm1,nb302_dyH2O(%esp)
        movapd %xmm2,nb302_dzH2O(%esp)
        mulsd  %xmm0,%xmm0
        mulsd  %xmm1,%xmm1
        mulsd  %xmm2,%xmm2
        movapd %xmm3,nb302_dxH2H1(%esp)
        movapd %xmm4,nb302_dyH2H1(%esp)
        movapd %xmm5,nb302_dzH2H1(%esp)
        mulsd  %xmm3,%xmm3
        mulsd  %xmm4,%xmm4
        mulsd  %xmm5,%xmm5
        addsd  %xmm1,%xmm0
        addsd  %xmm2,%xmm0
        addsd  %xmm3,%xmm4
        addsd  %xmm5,%xmm4
        movapd %xmm0,nb302_rsqH2O(%esp)
        movapd %xmm4,nb302_rsqH2H1(%esp)

        movapd nb302_ixH2(%esp),%xmm0
        movapd nb302_iyH2(%esp),%xmm1
        movapd nb302_izH2(%esp),%xmm2
        subsd  nb302_jxH2(%esp),%xmm0
        subsd  nb302_jyH2(%esp),%xmm1
        subsd  nb302_jzH2(%esp),%xmm2
        movapd %xmm0,nb302_dxH2H2(%esp)
        movapd %xmm1,nb302_dyH2H2(%esp)
        movapd %xmm2,nb302_dzH2H2(%esp)
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb302_rsqH2H2(%esp)

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
        movapd  nb302_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302_half(%esp),%xmm3   ## iter1 
        mulsd   nb302_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302_half(%esp),%xmm1   ## rinv 
        mulsd   nb302_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302_rinvH2H2(%esp)
        movapd %xmm5,nb302_rinvH2H1(%esp)

        movapd nb302_rsqOO(%esp),%xmm0
        movapd nb302_rsqOH1(%esp),%xmm4
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
        movapd  nb302_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb302_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302_half(%esp),%xmm1   ## rinv 
        mulsd   nb302_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb302_rinvOO(%esp)
        movapd %xmm5,nb302_rinvOH1(%esp)

        movapd nb302_rsqOH2(%esp),%xmm0
        movapd nb302_rsqH1O(%esp),%xmm4
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
        movapd  nb302_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302_half(%esp),%xmm3   ## iter1 
        mulsd   nb302_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302_half(%esp),%xmm1   ## rinv 
        mulsd   nb302_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302_rinvOH2(%esp)
        movapd %xmm5,nb302_rinvH1O(%esp)

        movapd nb302_rsqH1H1(%esp),%xmm0
        movapd nb302_rsqH1H2(%esp),%xmm4
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
        movapd  nb302_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302_half(%esp),%xmm3   ## iter1a 
        mulsd   nb302_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302_half(%esp),%xmm1   ## rinv 
        mulsd   nb302_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302_rinvH1H1(%esp)
        movapd %xmm5,nb302_rinvH1H2(%esp)

        movapd nb302_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb302_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb302_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb302_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb302_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb302_rinvH2O(%esp)

        movd %eax,%mm0
        ## start with OO interaction 
        movapd nb302_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        
        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOO(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addsd  nb302_vctot(%esp),%xmm5
        xorpd  %xmm2,%xmm2
    movlpd %xmm5,nb302_vctot(%esp)
        mulsd  nb302_tsc(%esp),%xmm3

        subsd  %xmm3,%xmm2
        mulsd  %xmm2,%xmm0

        movapd %xmm0,%xmm1
        movapd %xmm0,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb302_dxOO(%esp),%xmm0
        mulsd nb302_dyOO(%esp),%xmm1
        mulsd nb302_dzOO(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixO(%esp),%xmm0
        addsd nb302_fiyO(%esp),%xmm1
        addsd nb302_fizO(%esp),%xmm2
        movlpd %xmm3,nb302_fjxO(%esp)
        movlpd %xmm4,nb302_fjyO(%esp)
        movlpd %xmm5,nb302_fjzO(%esp)
        movlpd %xmm0,nb302_fixO(%esp)
        movlpd %xmm1,nb302_fiyO(%esp)
        movlpd %xmm2,nb302_fizO(%esp)

        ## O-H1 interaction 
        movapd nb302_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqOH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb302_vctot(%esp),%xmm5
    movlpd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb302_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb302_dxOH1(%esp),%xmm0
        mulsd nb302_dyOH1(%esp),%xmm1
        mulsd nb302_dzOH1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixO(%esp),%xmm0
        addsd nb302_fiyO(%esp),%xmm1
        addsd nb302_fizO(%esp),%xmm2
        movlpd %xmm3,nb302_fjxH1(%esp)
        movlpd %xmm4,nb302_fjyH1(%esp)
        movlpd %xmm5,nb302_fjzH1(%esp)
        movlpd %xmm0,nb302_fixO(%esp)
        movlpd %xmm1,nb302_fiyO(%esp)
        movlpd %xmm2,nb302_fizO(%esp)

        ## O-H2 interaction  
        movapd nb302_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqOH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb302_vctot(%esp),%xmm5
    movlpd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb302_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        xorpd %xmm3,%xmm3
        movapd %xmm3,%xmm4
        movapd %xmm3,%xmm5
        mulsd nb302_dxOH2(%esp),%xmm0
        mulsd nb302_dyOH2(%esp),%xmm1
        mulsd nb302_dzOH2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixO(%esp),%xmm0
        addsd nb302_fiyO(%esp),%xmm1
        addsd nb302_fizO(%esp),%xmm2
        movlpd %xmm3,nb302_fjxH2(%esp)
        movlpd %xmm4,nb302_fjyH2(%esp)
        movlpd %xmm5,nb302_fjzH2(%esp)
        movlpd %xmm0,nb302_fixO(%esp)
        movlpd %xmm1,nb302_fiyO(%esp)
        movlpd %xmm2,nb302_fizO(%esp)

        ## H1-O interaction 
        movapd nb302_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqH1O(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb302_vctot(%esp),%xmm5
    movlpd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb302_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxO(%esp),%xmm3
        movapd nb302_fjyO(%esp),%xmm4
        movapd nb302_fjzO(%esp),%xmm5
        mulsd nb302_dxH1O(%esp),%xmm0
        mulsd nb302_dyH1O(%esp),%xmm1
        mulsd nb302_dzH1O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixH1(%esp),%xmm0
        addsd nb302_fiyH1(%esp),%xmm1
        addsd nb302_fizH1(%esp),%xmm2
        movlpd %xmm3,nb302_fjxO(%esp)
        movlpd %xmm4,nb302_fjyO(%esp)
        movlpd %xmm5,nb302_fjzO(%esp)
        movlpd %xmm0,nb302_fixH1(%esp)
        movlpd %xmm1,nb302_fiyH1(%esp)
        movlpd %xmm2,nb302_fizH1(%esp)

        ## H1-H1 interaction 
        movapd nb302_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb302_vctot(%esp),%xmm5
    movlpd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb302_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxH1(%esp),%xmm3
        movapd nb302_fjyH1(%esp),%xmm4
        movapd nb302_fjzH1(%esp),%xmm5
        mulsd nb302_dxH1H1(%esp),%xmm0
        mulsd nb302_dyH1H1(%esp),%xmm1
        mulsd nb302_dzH1H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixH1(%esp),%xmm0
        addsd nb302_fiyH1(%esp),%xmm1
        addsd nb302_fizH1(%esp),%xmm2
        movlpd %xmm3,nb302_fjxH1(%esp)
        movlpd %xmm4,nb302_fjyH1(%esp)
        movlpd %xmm5,nb302_fjzH1(%esp)
        movlpd %xmm0,nb302_fixH1(%esp)
        movlpd %xmm1,nb302_fiyH1(%esp)
        movlpd %xmm2,nb302_fizH1(%esp)

        ## H1-H2 interaction 
        movapd nb302_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb302_vctot(%esp),%xmm5
    movlpd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb302_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxH2(%esp),%xmm3
        movapd nb302_fjyH2(%esp),%xmm4
        movapd nb302_fjzH2(%esp),%xmm5
        mulsd nb302_dxH1H2(%esp),%xmm0
        mulsd nb302_dyH1H2(%esp),%xmm1
        mulsd nb302_dzH1H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixH1(%esp),%xmm0
        addsd nb302_fiyH1(%esp),%xmm1
        addsd nb302_fizH1(%esp),%xmm2
        movlpd %xmm3,nb302_fjxH2(%esp)
        movlpd %xmm4,nb302_fjyH2(%esp)
        movlpd %xmm5,nb302_fjzH2(%esp)
        movlpd %xmm0,nb302_fixH1(%esp)
        movlpd %xmm1,nb302_fiyH1(%esp)
        movlpd %xmm2,nb302_fizH1(%esp)

        ## H2-O interaction 
        movapd nb302_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqH2O(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqOH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb302_vctot(%esp),%xmm5
    movlpd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb302_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxO(%esp),%xmm3
        movapd nb302_fjyO(%esp),%xmm4
        movapd nb302_fjzO(%esp),%xmm5
        mulsd nb302_dxH2O(%esp),%xmm0
        mulsd nb302_dyH2O(%esp),%xmm1
        mulsd nb302_dzH2O(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixH2(%esp),%xmm0
        addsd nb302_fiyH2(%esp),%xmm1
        addsd nb302_fizH2(%esp),%xmm2
        movlpd %xmm3,nb302_fjxO(%esp)
        movlpd %xmm4,nb302_fjyO(%esp)
        movlpd %xmm5,nb302_fjzO(%esp)
        movlpd %xmm0,nb302_fixH2(%esp)
        movlpd %xmm1,nb302_fiyH2(%esp)
        movlpd %xmm2,nb302_fizH2(%esp)

        ## H2-H1 interaction 
        movapd nb302_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb302_vctot(%esp),%xmm5
    movlpd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb302_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxH1(%esp),%xmm3
        movapd nb302_fjyH1(%esp),%xmm4
        movapd nb302_fjzH1(%esp),%xmm5
        mulsd nb302_dxH2H1(%esp),%xmm0
        mulsd nb302_dyH2H1(%esp),%xmm1
        mulsd nb302_dzH2H1(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixH2(%esp),%xmm0
        addsd nb302_fiyH2(%esp),%xmm1
        addsd nb302_fizH2(%esp),%xmm2
        movlpd %xmm3,nb302_fjxH1(%esp)
        movlpd %xmm4,nb302_fjyH1(%esp)
        movlpd %xmm5,nb302_fjzH1(%esp)
        movlpd %xmm0,nb302_fixH2(%esp)
        movlpd %xmm1,nb302_fiyH2(%esp)
        movlpd %xmm2,nb302_fizH2(%esp)

        ## H2-H2 interaction 
        movapd nb302_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        mulsd  nb302_two(%esp),%xmm7    ## two*Heps2 
        movapd nb302_qqHH(%esp),%xmm3
        addsd  %xmm6,%xmm7
        addsd  %xmm5,%xmm7 ## xmm7=FF 
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
        mulsd  %xmm7,%xmm3 ## fijC=FF*qq 
    ## at this point mm5 contains vcoul and xmm3 fijC 

    addsd  nb302_vctot(%esp),%xmm5
    movlpd %xmm5,nb302_vctot(%esp)
        xorpd  %xmm1,%xmm1
        mulsd  nb302_tsc(%esp),%xmm3
        mulsd  %xmm0,%xmm3
        subsd  %xmm3,%xmm1

        movapd %xmm1,%xmm0
        movapd %xmm1,%xmm2

        movapd nb302_fjxH2(%esp),%xmm3
        movapd nb302_fjyH2(%esp),%xmm4
        movapd nb302_fjzH2(%esp),%xmm5
        mulsd nb302_dxH2H2(%esp),%xmm0
        mulsd nb302_dyH2H2(%esp),%xmm1
        mulsd nb302_dzH2H2(%esp),%xmm2
        subsd %xmm0,%xmm3
        subsd %xmm1,%xmm4
        subsd %xmm2,%xmm5
        addsd nb302_fixH2(%esp),%xmm0
        addsd nb302_fiyH2(%esp),%xmm1
        addsd nb302_fizH2(%esp),%xmm2
        movlpd %xmm3,nb302_fjxH2(%esp)
        movlpd %xmm4,nb302_fjyH2(%esp)
        movlpd %xmm5,nb302_fjzH2(%esp)
        movlpd %xmm0,nb302_fixH2(%esp)
        movlpd %xmm1,nb302_fiyH2(%esp)
        movlpd %xmm2,nb302_fizH2(%esp)

        movl nb302_faction(%ebp),%edi

        movd %mm0,%eax

        ## Did all interactions - now update j forces 
        movlpd (%edi,%eax,8),%xmm0
        movlpd 8(%edi,%eax,8),%xmm1
        movlpd 16(%edi,%eax,8),%xmm2
        movlpd 24(%edi,%eax,8),%xmm3
        movlpd 32(%edi,%eax,8),%xmm4
        movlpd 40(%edi,%eax,8),%xmm5
        movlpd 48(%edi,%eax,8),%xmm6
        movlpd 56(%edi,%eax,8),%xmm7
        addsd nb302_fjxO(%esp),%xmm0
        addsd nb302_fjyO(%esp),%xmm1
        addsd nb302_fjzO(%esp),%xmm2
        addsd nb302_fjxH1(%esp),%xmm3
        addsd nb302_fjyH1(%esp),%xmm4
        addsd nb302_fjzH1(%esp),%xmm5
        addsd nb302_fjxH2(%esp),%xmm6
        addsd nb302_fjyH2(%esp),%xmm7
        movlpd %xmm0,(%edi,%eax,8)
        movlpd %xmm1,8(%edi,%eax,8)
        movlpd %xmm2,16(%edi,%eax,8)
        movlpd %xmm3,24(%edi,%eax,8)
        movlpd %xmm4,32(%edi,%eax,8)
        movlpd %xmm5,40(%edi,%eax,8)
        movlpd %xmm6,48(%edi,%eax,8)
        movlpd %xmm7,56(%edi,%eax,8)

        movlpd 64(%edi,%eax,8),%xmm0
        addsd nb302_fjzH2(%esp),%xmm0
        movlpd %xmm0,64(%edi,%eax,8)

_nb_kernel302_ia32_sse2.nb302_updateouterdata: 
        movl  nb302_ii3(%esp),%ecx
        movl  nb302_faction(%ebp),%edi
        movl  nb302_fshift(%ebp),%esi
        movl  nb302_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movapd nb302_fixO(%esp),%xmm0
        movapd nb302_fiyO(%esp),%xmm1
        movapd nb302_fizO(%esp),%xmm2

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
        movapd nb302_fixH1(%esp),%xmm0
        movapd nb302_fiyH1(%esp),%xmm1
        movapd nb302_fizH1(%esp),%xmm2

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
        movapd nb302_fixH2(%esp),%xmm0
        movapd nb302_fiyH2(%esp),%xmm1
        movapd nb302_fizH2(%esp),%xmm2

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
        movl nb302_n(%esp),%esi
        ## get group index for i particle 
        movl  nb302_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb302_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb302_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb302_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel302_ia32_sse2.nb302_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb302_n(%esp)
        jmp _nb_kernel302_ia32_sse2.nb302_outer
_nb_kernel302_ia32_sse2.nb302_outerend: 
        ## check if more outer neighborlists remain
        movl  nb302_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel302_ia32_sse2.nb302_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel302_ia32_sse2.nb302_threadloop
_nb_kernel302_ia32_sse2.nb302_end: 
        emms

        movl nb302_nouter(%esp),%eax
        movl nb302_ninner(%esp),%ebx
        movl nb302_outeriter(%ebp),%ecx
        movl nb302_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb302_salign(%esp),%eax
        addl %eax,%esp
        addl $1464,%esp
        popl %edi
        popl %esi
    popl %edx
    popl %ecx
    popl %ebx
    popl %eax
        leave
        ret


.globl nb_kernel302nf_ia32_sse2
.globl _nb_kernel302nf_ia32_sse2
nb_kernel302nf_ia32_sse2:       
_nb_kernel302nf_ia32_sse2:      
.set nb302nf_p_nri, 8
.set nb302nf_iinr, 12
.set nb302nf_jindex, 16
.set nb302nf_jjnr, 20
.set nb302nf_shift, 24
.set nb302nf_shiftvec, 28
.set nb302nf_fshift, 32
.set nb302nf_gid, 36
.set nb302nf_pos, 40
.set nb302nf_faction, 44
.set nb302nf_charge, 48
.set nb302nf_p_facel, 52
.set nb302nf_argkrf, 56
.set nb302nf_argcrf, 60
.set nb302nf_Vc, 64
.set nb302nf_type, 68
.set nb302nf_p_ntype, 72
.set nb302nf_vdwparam, 76
.set nb302nf_Vvdw, 80
.set nb302nf_p_tabscale, 84
.set nb302nf_VFtab, 88
.set nb302nf_invsqrta, 92
.set nb302nf_dvda, 96
.set nb302nf_p_gbtabscale, 100
.set nb302nf_GBtab, 104
.set nb302nf_p_nthreads, 108
.set nb302nf_count, 112
.set nb302nf_mtx, 116
.set nb302nf_outeriter, 120
.set nb302nf_inneriter, 124
.set nb302nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb302nf_ixO, 0
.set nb302nf_iyO, 16
.set nb302nf_izO, 32
.set nb302nf_ixH1, 48
.set nb302nf_iyH1, 64
.set nb302nf_izH1, 80
.set nb302nf_ixH2, 96
.set nb302nf_iyH2, 112
.set nb302nf_izH2, 128
.set nb302nf_jxO, 144
.set nb302nf_jyO, 160
.set nb302nf_jzO, 176
.set nb302nf_jxH1, 192
.set nb302nf_jyH1, 208
.set nb302nf_jzH1, 224
.set nb302nf_jxH2, 240
.set nb302nf_jyH2, 256
.set nb302nf_jzH2, 272
.set nb302nf_qqOO, 288
.set nb302nf_qqOH, 304
.set nb302nf_qqHH, 320
.set nb302nf_tsc, 336
.set nb302nf_vctot, 352
.set nb302nf_half, 368
.set nb302nf_three, 384
.set nb302nf_rsqOO, 400
.set nb302nf_rsqOH1, 416
.set nb302nf_rsqOH2, 432
.set nb302nf_rsqH1O, 448
.set nb302nf_rsqH1H1, 464
.set nb302nf_rsqH1H2, 480
.set nb302nf_rsqH2O, 496
.set nb302nf_rsqH2H1, 512
.set nb302nf_rsqH2H2, 528
.set nb302nf_rinvOO, 544
.set nb302nf_rinvOH1, 560
.set nb302nf_rinvOH2, 576
.set nb302nf_rinvH1O, 592
.set nb302nf_rinvH1H1, 608
.set nb302nf_rinvH1H2, 624
.set nb302nf_rinvH2O, 640
.set nb302nf_rinvH2H1, 656
.set nb302nf_rinvH2H2, 672
.set nb302nf_is3, 688
.set nb302nf_ii3, 692
.set nb302nf_innerjjnr, 696
.set nb302nf_innerk, 700
.set nb302nf_n, 704
.set nb302nf_nn1, 708
.set nb302nf_nri, 712
.set nb302nf_nouter, 716
.set nb302nf_ninner, 720
.set nb302nf_salign, 724
        pushl %ebp
        movl %esp,%ebp
    pushl %eax
    pushl %ebx
    pushl %ecx
    pushl %edx
        pushl %esi
        pushl %edi
        subl $728,%esp          ## local stack space 
        movl %esp,%eax
        andl $0xf,%eax
        subl %eax,%esp
        movl %eax,nb302nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb302nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb302nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb302nf_nouter(%esp)
        movl %eax,nb302nf_ninner(%esp)


        ## create constant floating-point factors on stack
        movl $0x00000000,%eax   ## lower half of double 0.5 IEEE (hex)
        movl $0x3fe00000,%ebx
        movl %eax,nb302nf_half(%esp)
        movl %ebx,nb302nf_half+4(%esp)
        movsd nb302nf_half(%esp),%xmm1
        shufpd $0,%xmm1,%xmm1  ## splat to all elements
        movapd %xmm1,%xmm3
        addpd  %xmm3,%xmm3      ## 1.0
        movapd %xmm3,%xmm2
        addpd  %xmm2,%xmm2      ## 2.0
        addpd  %xmm2,%xmm3      ## 3.0
        movapd %xmm1,nb302nf_half(%esp)
        movapd %xmm3,nb302nf_three(%esp)

        movl nb302nf_p_tabscale(%ebp),%eax
        movsd (%eax),%xmm3
        shufpd $0,%xmm3,%xmm3
        movapd %xmm3,nb302nf_tsc(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb302nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb302nf_charge(%ebp),%edx
        movsd (%edx,%ebx,8),%xmm3
        movsd %xmm3,%xmm4
        movsd 8(%edx,%ebx,8),%xmm5
        movl nb302nf_p_facel(%ebp),%esi
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
        movapd %xmm3,nb302nf_qqOO(%esp)
        movapd %xmm4,nb302nf_qqOH(%esp)
        movapd %xmm5,nb302nf_qqHH(%esp)

_nb_kernel302nf_ia32_sse2.nb302nf_threadloop: 
        movl  nb302nf_count(%ebp),%esi          ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel302nf_ia32_sse2.nb302nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                           ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel302nf_ia32_sse2.nb302nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb302nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb302nf_n(%esp)
        movl %ebx,nb302nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel302nf_ia32_sse2.nb302nf_outerstart
        jmp _nb_kernel302nf_ia32_sse2.nb302nf_end

_nb_kernel302nf_ia32_sse2.nb302nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb302nf_nouter(%esp),%ebx
        movl %ebx,nb302nf_nouter(%esp)

_nb_kernel302nf_ia32_sse2.nb302nf_outer: 
        movl  nb302nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb302nf_is3(%esp)            ## store is3 

        movl  nb302nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movsd (%eax,%ebx,8),%xmm0
        movsd 8(%eax,%ebx,8),%xmm1
        movsd 16(%eax,%ebx,8),%xmm2

        movl  nb302nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb302nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb302nf_ii3(%esp)

        movapd %xmm0,%xmm3
        movapd %xmm1,%xmm4
        movapd %xmm2,%xmm5
        addsd (%eax,%ebx,8),%xmm3
        addsd 8(%eax,%ebx,8),%xmm4
        addsd 16(%eax,%ebx,8),%xmm5
        shufpd $0,%xmm3,%xmm3
        shufpd $0,%xmm4,%xmm4
        shufpd $0,%xmm5,%xmm5
        movapd %xmm3,nb302nf_ixO(%esp)
        movapd %xmm4,nb302nf_iyO(%esp)
        movapd %xmm5,nb302nf_izO(%esp)

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
        movapd %xmm0,nb302nf_ixH1(%esp)
        movapd %xmm1,nb302nf_iyH1(%esp)
        movapd %xmm2,nb302nf_izH1(%esp)
        movapd %xmm3,nb302nf_ixH2(%esp)
        movapd %xmm4,nb302nf_iyH2(%esp)
        movapd %xmm5,nb302nf_izH2(%esp)

        ## clear vctot 
        xorpd %xmm4,%xmm4
        movapd %xmm4,nb302nf_vctot(%esp)

        movl  nb302nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx             ## jindex[n] 
        movl  4(%eax,%esi,4),%edx            ## jindex[n+1] 
        subl  %ecx,%edx              ## number of innerloop atoms 

        movl  nb302nf_pos(%ebp),%esi
        movl  nb302nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb302nf_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $2,%edx
        addl  nb302nf_ninner(%esp),%ecx
        movl  %ecx,nb302nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb302nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel302nf_ia32_sse2.nb302nf_unroll_loop
        jmp   _nb_kernel302nf_ia32_sse2.nb302nf_checksingle
_nb_kernel302nf_ia32_sse2.nb302nf_unroll_loop: 
        ## twice unrolled innerloop here 
        movl  nb302nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        movl  4(%edx),%ebx

        addl $8,nb302nf_innerjjnr(%esp)            ## advance pointer (unrolled 2) 

        movl nb302nf_pos(%ebp),%esi        ## base of pos[] 

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
        movapd  %xmm2,nb302nf_jxO(%esp)
        movapd  %xmm3,nb302nf_jyO(%esp)
        movapd  %xmm4,nb302nf_jzO(%esp)
        movapd  %xmm5,nb302nf_jxH1(%esp)
        movapd  %xmm6,nb302nf_jyH1(%esp)
        movapd  %xmm7,nb302nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movhpd 48(%esi,%ebx,8),%xmm2
        movhpd 56(%esi,%ebx,8),%xmm3
        movhpd 64(%esi,%ebx,8),%xmm4
        movapd  %xmm2,nb302nf_jxH2(%esp)
        movapd  %xmm3,nb302nf_jyH2(%esp)
        movapd  %xmm4,nb302nf_jzH2(%esp)

        movapd nb302nf_ixO(%esp),%xmm0
        movapd nb302nf_iyO(%esp),%xmm1
        movapd nb302nf_izO(%esp),%xmm2
        movapd nb302nf_ixO(%esp),%xmm3
        movapd nb302nf_iyO(%esp),%xmm4
        movapd nb302nf_izO(%esp),%xmm5
        subpd  nb302nf_jxO(%esp),%xmm0
        subpd  nb302nf_jyO(%esp),%xmm1
        subpd  nb302nf_jzO(%esp),%xmm2
        subpd  nb302nf_jxH1(%esp),%xmm3
        subpd  nb302nf_jyH1(%esp),%xmm4
        subpd  nb302nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb302nf_rsqOO(%esp)
        movapd %xmm3,nb302nf_rsqOH1(%esp)

        movapd nb302nf_ixO(%esp),%xmm0
        movapd nb302nf_iyO(%esp),%xmm1
        movapd nb302nf_izO(%esp),%xmm2
        movapd nb302nf_ixH1(%esp),%xmm3
        movapd nb302nf_iyH1(%esp),%xmm4
        movapd nb302nf_izH1(%esp),%xmm5
        subpd  nb302nf_jxH2(%esp),%xmm0
        subpd  nb302nf_jyH2(%esp),%xmm1
        subpd  nb302nf_jzH2(%esp),%xmm2
        subpd  nb302nf_jxO(%esp),%xmm3
        subpd  nb302nf_jyO(%esp),%xmm4
        subpd  nb302nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb302nf_rsqOH2(%esp)
        movapd %xmm3,nb302nf_rsqH1O(%esp)

        movapd nb302nf_ixH1(%esp),%xmm0
        movapd nb302nf_iyH1(%esp),%xmm1
        movapd nb302nf_izH1(%esp),%xmm2
        movapd nb302nf_ixH1(%esp),%xmm3
        movapd nb302nf_iyH1(%esp),%xmm4
        movapd nb302nf_izH1(%esp),%xmm5
        subpd  nb302nf_jxH1(%esp),%xmm0
        subpd  nb302nf_jyH1(%esp),%xmm1
        subpd  nb302nf_jzH1(%esp),%xmm2
        subpd  nb302nf_jxH2(%esp),%xmm3
        subpd  nb302nf_jyH2(%esp),%xmm4
        subpd  nb302nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb302nf_rsqH1H1(%esp)
        movapd %xmm3,nb302nf_rsqH1H2(%esp)

        movapd nb302nf_ixH2(%esp),%xmm0
        movapd nb302nf_iyH2(%esp),%xmm1
        movapd nb302nf_izH2(%esp),%xmm2
        movapd nb302nf_ixH2(%esp),%xmm3
        movapd nb302nf_iyH2(%esp),%xmm4
        movapd nb302nf_izH2(%esp),%xmm5
        subpd  nb302nf_jxO(%esp),%xmm0
        subpd  nb302nf_jyO(%esp),%xmm1
        subpd  nb302nf_jzO(%esp),%xmm2
        subpd  nb302nf_jxH1(%esp),%xmm3
        subpd  nb302nf_jyH1(%esp),%xmm4
        subpd  nb302nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb302nf_rsqH2O(%esp)
        movapd %xmm4,nb302nf_rsqH2H1(%esp)

        movapd nb302nf_ixH2(%esp),%xmm0
        movapd nb302nf_iyH2(%esp),%xmm1
        movapd nb302nf_izH2(%esp),%xmm2
        subpd  nb302nf_jxH2(%esp),%xmm0
        subpd  nb302nf_jyH2(%esp),%xmm1
        subpd  nb302nf_jzH2(%esp),%xmm2
        mulpd %xmm0,%xmm0
        mulpd %xmm1,%xmm1
        mulpd %xmm2,%xmm2
        addpd %xmm1,%xmm0
        addpd %xmm2,%xmm0
        movapd %xmm0,nb302nf_rsqH2H2(%esp)

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
        movapd  nb302nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb302nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb302nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvH2H2(%esp)
        movapd %xmm5,nb302nf_rinvH2H1(%esp)

        movapd nb302nf_rsqOO(%esp),%xmm0
        movapd nb302nf_rsqOH1(%esp),%xmm4
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
        movapd  nb302nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%esp),%xmm3   ## iter1 of  
        mulpd   nb302nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb302nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb302nf_rinvOO(%esp)
        movapd %xmm5,nb302nf_rinvOH1(%esp)

        movapd nb302nf_rsqOH2(%esp),%xmm0
        movapd nb302nf_rsqH1O(%esp),%xmm4
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
        movapd  nb302nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%esp),%xmm3   ## iter1 
        mulpd   nb302nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb302nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvOH2(%esp)
        movapd %xmm5,nb302nf_rinvH1O(%esp)

        movapd nb302nf_rsqH1H1(%esp),%xmm0
        movapd nb302nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb302nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subpd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%esp),%xmm3   ## iter1a 
        mulpd   nb302nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        mulpd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulpd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subpd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulpd   nb302nf_half(%esp),%xmm1   ## rinv 
        mulpd   nb302nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvH1H1(%esp)
        movapd %xmm5,nb302nf_rinvH1H2(%esp)

        movapd nb302nf_rsqH2O(%esp),%xmm0
        cvtpd2ps %xmm0,%xmm1
        rsqrtps %xmm1,%xmm1
        cvtps2pd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulpd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb302nf_three(%esp),%xmm3
        mulpd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subpd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb302nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulpd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb302nf_three(%esp),%xmm1
        mulpd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subpd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulpd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulpd   nb302nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb302nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb302nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqOO(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## O-H1 interaction 
        movapd nb302nf_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqOH1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqOH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## O-H2 interaction  
        movapd nb302nf_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqOH2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqOH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## H1-O interaction 
        movapd nb302nf_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH1O(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1

        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqOH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## H1-H1 interaction 
        movapd nb302nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## H1-H2 interaction 
        movapd nb302nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## H2-O interaction 
        movapd nb302nf_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH2O(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqOH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb302nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb302nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulpd  nb302nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulpd  nb302nf_tsc(%esp),%xmm1
        cvttpd2pi %xmm1,%mm6    ## mm6 = lu idx 
        cvtpi2pd %mm6,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        pslld $2,%mm6           ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi
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
        movapd nb302nf_qqHH(%esp),%xmm3
        mulpd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addpd  %xmm4,%xmm5 ## xmm5=VV 
        mulpd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addpd  nb302nf_vctot(%esp),%xmm5
    movapd %xmm5,nb302nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $2,nb302nf_innerk(%esp)
        jl    _nb_kernel302nf_ia32_sse2.nb302nf_checksingle
        jmp   _nb_kernel302nf_ia32_sse2.nb302nf_unroll_loop
_nb_kernel302nf_ia32_sse2.nb302nf_checksingle: 
        movl  nb302nf_innerk(%esp),%edx
        andl  $1,%edx
        jnz   _nb_kernel302nf_ia32_sse2.nb302nf_dosingle
        jmp   _nb_kernel302nf_ia32_sse2.nb302nf_updateouterdata
_nb_kernel302nf_ia32_sse2.nb302nf_dosingle: 
        movl  nb302nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax

        movl nb302nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        movlpd (%esi,%eax,8),%xmm2
        movlpd 8(%esi,%eax,8),%xmm3
        movlpd 16(%esi,%eax,8),%xmm4
        movlpd 24(%esi,%eax,8),%xmm5
        movlpd 32(%esi,%eax,8),%xmm6
        movlpd 40(%esi,%eax,8),%xmm7
        movapd  %xmm2,nb302nf_jxO(%esp)
        movapd  %xmm3,nb302nf_jyO(%esp)
        movapd  %xmm4,nb302nf_jzO(%esp)
        movapd  %xmm5,nb302nf_jxH1(%esp)
        movapd  %xmm6,nb302nf_jyH1(%esp)
        movapd  %xmm7,nb302nf_jzH1(%esp)
        movlpd 48(%esi,%eax,8),%xmm2
        movlpd 56(%esi,%eax,8),%xmm3
        movlpd 64(%esi,%eax,8),%xmm4
        movapd  %xmm2,nb302nf_jxH2(%esp)
        movapd  %xmm3,nb302nf_jyH2(%esp)
        movapd  %xmm4,nb302nf_jzH2(%esp)

        movapd nb302nf_ixO(%esp),%xmm0
        movapd nb302nf_iyO(%esp),%xmm1
        movapd nb302nf_izO(%esp),%xmm2
        movapd nb302nf_ixO(%esp),%xmm3
        movapd nb302nf_iyO(%esp),%xmm4
        movapd nb302nf_izO(%esp),%xmm5
        subsd  nb302nf_jxO(%esp),%xmm0
        subsd  nb302nf_jyO(%esp),%xmm1
        subsd  nb302nf_jzO(%esp),%xmm2
        subsd  nb302nf_jxH1(%esp),%xmm3
        subsd  nb302nf_jyH1(%esp),%xmm4
        subsd  nb302nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb302nf_rsqOO(%esp)
        movapd %xmm3,nb302nf_rsqOH1(%esp)

        movapd nb302nf_ixO(%esp),%xmm0
        movapd nb302nf_iyO(%esp),%xmm1
        movapd nb302nf_izO(%esp),%xmm2
        movapd nb302nf_ixH1(%esp),%xmm3
        movapd nb302nf_iyH1(%esp),%xmm4
        movapd nb302nf_izH1(%esp),%xmm5
        subsd  nb302nf_jxH2(%esp),%xmm0
        subsd  nb302nf_jyH2(%esp),%xmm1
        subsd  nb302nf_jzH2(%esp),%xmm2
        subsd  nb302nf_jxO(%esp),%xmm3
        subsd  nb302nf_jyO(%esp),%xmm4
        subsd  nb302nf_jzO(%esp),%xmm5
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
        movapd %xmm0,nb302nf_rsqOH2(%esp)
        movapd %xmm3,nb302nf_rsqH1O(%esp)

        movapd nb302nf_ixH1(%esp),%xmm0
        movapd nb302nf_iyH1(%esp),%xmm1
        movapd nb302nf_izH1(%esp),%xmm2
        movapd nb302nf_ixH1(%esp),%xmm3
        movapd nb302nf_iyH1(%esp),%xmm4
        movapd nb302nf_izH1(%esp),%xmm5
        subsd  nb302nf_jxH1(%esp),%xmm0
        subsd  nb302nf_jyH1(%esp),%xmm1
        subsd  nb302nf_jzH1(%esp),%xmm2
        subsd  nb302nf_jxH2(%esp),%xmm3
        subsd  nb302nf_jyH2(%esp),%xmm4
        subsd  nb302nf_jzH2(%esp),%xmm5
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
        movapd %xmm0,nb302nf_rsqH1H1(%esp)
        movapd %xmm3,nb302nf_rsqH1H2(%esp)

        movapd nb302nf_ixH2(%esp),%xmm0
        movapd nb302nf_iyH2(%esp),%xmm1
        movapd nb302nf_izH2(%esp),%xmm2
        movapd nb302nf_ixH2(%esp),%xmm3
        movapd nb302nf_iyH2(%esp),%xmm4
        movapd nb302nf_izH2(%esp),%xmm5
        subsd  nb302nf_jxO(%esp),%xmm0
        subsd  nb302nf_jyO(%esp),%xmm1
        subsd  nb302nf_jzO(%esp),%xmm2
        subsd  nb302nf_jxH1(%esp),%xmm3
        subsd  nb302nf_jyH1(%esp),%xmm4
        subsd  nb302nf_jzH1(%esp),%xmm5
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
        movapd %xmm0,nb302nf_rsqH2O(%esp)
        movapd %xmm4,nb302nf_rsqH2H1(%esp)

        movapd nb302nf_ixH2(%esp),%xmm0
        movapd nb302nf_iyH2(%esp),%xmm1
        movapd nb302nf_izH2(%esp),%xmm2
        subsd  nb302nf_jxH2(%esp),%xmm0
        subsd  nb302nf_jyH2(%esp),%xmm1
        subsd  nb302nf_jzH2(%esp),%xmm2
        mulsd %xmm0,%xmm0
        mulsd %xmm1,%xmm1
        mulsd %xmm2,%xmm2
        addsd %xmm1,%xmm0
        addsd %xmm2,%xmm0
        movapd %xmm0,nb302nf_rsqH2H2(%esp)

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
        movapd  nb302nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb302nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb302nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvH2H2(%esp)
        movapd %xmm5,nb302nf_rinvH2H1(%esp)

        movapd nb302nf_rsqOO(%esp),%xmm0
        movapd nb302nf_rsqOH1(%esp),%xmm4
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
        movapd  nb302nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%esp),%xmm3   ## iter1 of  
        mulsd   nb302nf_half(%esp),%xmm7   ## iter1 of  

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb302nf_half(%esp),%xmm5   ## rinv
        movapd %xmm1,nb302nf_rinvOO(%esp)
        movapd %xmm5,nb302nf_rinvOH1(%esp)

        movapd nb302nf_rsqOH2(%esp),%xmm0
        movapd nb302nf_rsqH1O(%esp),%xmm4
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
        movapd  nb302nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%esp),%xmm3   ## iter1 
        mulsd   nb302nf_half(%esp),%xmm7   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb302nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvOH2(%esp)
        movapd %xmm5,nb302nf_rinvH1O(%esp)

        movapd nb302nf_rsqH1H1(%esp),%xmm0
        movapd nb302nf_rsqH1H2(%esp),%xmm4
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
        movapd  nb302nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm5     ## rsqB*luB*luB         
        movapd  %xmm3,%xmm7
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        subsd   %xmm5,%xmm7     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm7     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%esp),%xmm3   ## iter1a 
        mulsd   nb302nf_half(%esp),%xmm7   ## iter1b 

        movapd  %xmm3,%xmm2     ## copy of luA 
        movapd  %xmm7,%xmm6     ## copy of luB 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        mulsd   %xmm7,%xmm7     ## luB*luB 
        movapd  nb302nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        mulsd   %xmm4,%xmm7     ## rsqB*luB*luB         
        movapd  %xmm1,%xmm5
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        subsd   %xmm7,%xmm5     ## 3-rsqB*luB*luB 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   %xmm6,%xmm5     ## luB*(3-rsqB*luB*luB) 
        mulsd   nb302nf_half(%esp),%xmm1   ## rinv 
        mulsd   nb302nf_half(%esp),%xmm5   ## rinv 
        movapd %xmm1,nb302nf_rinvH1H1(%esp)
        movapd %xmm5,nb302nf_rinvH1H2(%esp)

        movapd nb302nf_rsqH2O(%esp),%xmm0
        cvtsd2ss %xmm0,%xmm1
        rsqrtss %xmm1,%xmm1
        cvtss2sd %xmm1,%xmm1

        movapd  %xmm1,%xmm2     ## copy of luA 
        mulsd   %xmm1,%xmm1     ## luA*luA 
        movapd  nb302nf_three(%esp),%xmm3
        mulsd   %xmm0,%xmm1     ## rsqA*luA*luA 
        subsd   %xmm1,%xmm3     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm3     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb302nf_half(%esp),%xmm3   ## iter1 

        movapd  %xmm3,%xmm2     ## copy of luA 
        mulsd   %xmm3,%xmm3     ## luA*luA 
        movapd  nb302nf_three(%esp),%xmm1
        mulsd   %xmm0,%xmm3     ## rsqA*luA*luA 
        subsd   %xmm3,%xmm1     ## 3-rsqA*luA*luA 
        mulsd   %xmm2,%xmm1     ## luA*(3-rsqA*luA*luA) 
        mulsd   nb302nf_half(%esp),%xmm1   ## rinv 
        movapd %xmm1,nb302nf_rinvH2O(%esp)

        ## start with OO interaction 
        movapd nb302nf_rinvOO(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqOO(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        
        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1  
        unpckhpd %xmm3,%xmm7    ## H1  
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOO(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 
    ## increment vcoul - then we can get rid of mm5 
    ## update vctot 
    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

        ## O-H1 interaction 
        movapd nb302nf_rinvOH1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqOH1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1  
        unpckhpd %xmm3,%xmm5    ## F1  

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul  

    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

        ## O-H2 interaction  
        movapd nb302nf_rinvOH2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqOH2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

        ## H1-O interaction 
        movapd nb302nf_rinvH1O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH1O(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1

        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

        ## H1-H1 interaction 
        movapd nb302nf_rinvH1H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH1H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

        ## H1-H2 interaction 
        movapd nb302nf_rinvH1H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH1H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

        ## H2-O interaction 
        movapd nb302nf_rinvH2O(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH2O(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqOH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

        ## H2-H1 interaction 
        movapd nb302nf_rinvH2H1(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH2H1(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subpd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulpd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

        ## H2-H2 interaction 
        movapd nb302nf_rinvH2H2(%esp),%xmm0
        movapd %xmm0,%xmm1
        mulsd  nb302nf_rsqH2H2(%esp),%xmm1   ## xmm1=r 
        mulsd  nb302nf_tsc(%esp),%xmm1
        cvttsd2si %xmm1,%eax    ## mm6 = lu idx 
        cvtsi2sd %eax,%xmm6
        subsd %xmm6,%xmm1       ## xmm1=eps 
        movapd %xmm1,%xmm2
        mulsd  %xmm2,%xmm2      ## xmm2=eps2 

        shll $2,%eax            ## idx *= 4 
        movl nb302nf_VFtab(%ebp),%esi

        movlpd (%esi,%eax,8),%xmm4      ## Y1   
        movhpd 8(%esi,%eax,8),%xmm4     ## Y1 F1        

        xorpd %xmm3,%xmm3
        movapd %xmm4,%xmm5
        unpcklpd %xmm3,%xmm4    ## Y1 
        unpckhpd %xmm3,%xmm5    ## F1 

        movlpd 16(%esi,%eax,8),%xmm6    ## G1   
        movhpd 24(%esi,%eax,8),%xmm6    ## G1 H1        

        xorpd %xmm3,%xmm3
        movapd %xmm6,%xmm7
        unpcklpd %xmm3,%xmm6    ## G1 
        unpckhpd %xmm3,%xmm7    ## H1 
        ## coulomb table ready, in xmm4-xmm7            
        mulsd  %xmm1,%xmm6      ## xmm6=Geps 
        mulsd  %xmm2,%xmm7      ## xmm7=Heps2 
        addsd  %xmm6,%xmm5
        addsd  %xmm7,%xmm5      ## xmm5=Fp      
        movapd nb302nf_qqHH(%esp),%xmm3
        mulsd  %xmm1,%xmm5 ## xmm5=eps*Fp 
        addsd  %xmm4,%xmm5 ## xmm5=VV 
        mulsd  %xmm3,%xmm5 ## vcoul=qq*VV  
    ## at this point mm5 contains vcoul 

    addsd  nb302nf_vctot(%esp),%xmm5
    movlpd %xmm5,nb302nf_vctot(%esp)

_nb_kernel302nf_ia32_sse2.nb302nf_updateouterdata: 
        ## get group index for i particle 
        ## get n from stack
        movl nb302nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb302nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movapd nb302nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addsd  %xmm6,%xmm7      ## low xmm7 has the sum now 

        ## add earlier value from mem 
        movl  nb302nf_Vc(%ebp),%eax
        addsd (%eax,%edx,8),%xmm7
        ## move back to mem 
        movsd %xmm7,(%eax,%edx,8)

        ## finish if last 
        movl nb302nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel302nf_ia32_sse2.nb302nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb302nf_n(%esp)
        jmp _nb_kernel302nf_ia32_sse2.nb302nf_outer
_nb_kernel302nf_ia32_sse2.nb302nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb302nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel302nf_ia32_sse2.nb302nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel302nf_ia32_sse2.nb302nf_threadloop
_nb_kernel302nf_ia32_sse2.nb302nf_end: 
        emms

        movl nb302nf_nouter(%esp),%eax
        movl nb302nf_ninner(%esp),%ebx
        movl nb302nf_outeriter(%ebp),%ecx
        movl nb302nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb302nf_salign(%esp),%eax
        addl %eax,%esp
        addl $728,%esp
        popl %edi
        popl %esi
        popl %edx
        popl %ecx
        popl %ebx
        popl %eax
        leave
        ret



