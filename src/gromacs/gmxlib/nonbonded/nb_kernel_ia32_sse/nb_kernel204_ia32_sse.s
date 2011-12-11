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



.globl nb_kernel204_ia32_sse
.globl _nb_kernel204_ia32_sse
nb_kernel204_ia32_sse:  
_nb_kernel204_ia32_sse: 
.set nb204_p_nri, 8
.set nb204_iinr, 12
.set nb204_jindex, 16
.set nb204_jjnr, 20
.set nb204_shift, 24
.set nb204_shiftvec, 28
.set nb204_fshift, 32
.set nb204_gid, 36
.set nb204_pos, 40
.set nb204_faction, 44
.set nb204_charge, 48
.set nb204_p_facel, 52
.set nb204_argkrf, 56
.set nb204_argcrf, 60
.set nb204_Vc, 64
.set nb204_type, 68
.set nb204_p_ntype, 72
.set nb204_vdwparam, 76
.set nb204_Vvdw, 80
.set nb204_p_tabscale, 84
.set nb204_VFtab, 88
.set nb204_invsqrta, 92
.set nb204_dvda, 96
.set nb204_p_gbtabscale, 100
.set nb204_GBtab, 104
.set nb204_p_nthreads, 108
.set nb204_count, 112
.set nb204_mtx, 116
.set nb204_outeriter, 120
.set nb204_inneriter, 124
.set nb204_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb204_ixH1, 0
.set nb204_iyH1, 16
.set nb204_izH1, 32
.set nb204_ixH2, 48
.set nb204_iyH2, 64
.set nb204_izH2, 80
.set nb204_ixM, 96
.set nb204_iyM, 112
.set nb204_izM, 128
.set nb204_jxH1, 144
.set nb204_jyH1, 160
.set nb204_jzH1, 176
.set nb204_jxH2, 192
.set nb204_jyH2, 208
.set nb204_jzH2, 224
.set nb204_jxM, 240
.set nb204_jyM, 256
.set nb204_jzM, 272
.set nb204_dxH1H1, 288
.set nb204_dyH1H1, 304
.set nb204_dzH1H1, 320
.set nb204_dxH1H2, 336
.set nb204_dyH1H2, 352
.set nb204_dzH1H2, 368
.set nb204_dxH1M, 384
.set nb204_dyH1M, 400
.set nb204_dzH1M, 416
.set nb204_dxH2H1, 432
.set nb204_dyH2H1, 448
.set nb204_dzH2H1, 464
.set nb204_dxH2H2, 480
.set nb204_dyH2H2, 496
.set nb204_dzH2H2, 512
.set nb204_dxH2M, 528
.set nb204_dyH2M, 544
.set nb204_dzH2M, 560
.set nb204_dxMH1, 576
.set nb204_dyMH1, 592
.set nb204_dzMH1, 608
.set nb204_dxMH2, 624
.set nb204_dyMH2, 640
.set nb204_dzMH2, 656
.set nb204_dxMM, 672
.set nb204_dyMM, 688
.set nb204_dzMM, 704
.set nb204_qqHH, 720
.set nb204_qqMH, 736
.set nb204_qqMM, 752
.set nb204_vctot, 768
.set nb204_fixH1, 784
.set nb204_fiyH1, 800
.set nb204_fizH1, 816
.set nb204_fixH2, 832
.set nb204_fiyH2, 848
.set nb204_fizH2, 864
.set nb204_fixM, 880
.set nb204_fiyM, 896
.set nb204_fizM, 912
.set nb204_fjxH1, 928
.set nb204_fjyH1, 944
.set nb204_fjzH1, 960
.set nb204_fjxH2, 976
.set nb204_fjyH2, 992
.set nb204_fjzH2, 1008
.set nb204_fjxM, 1024
.set nb204_fjyM, 1040
.set nb204_fjzM, 1056
.set nb204_fjzMb, 1060
.set nb204_fjzMc, 1064
.set nb204_fjzMd, 1068
.set nb204_half, 1072
.set nb204_three, 1088
.set nb204_rsqH1H1, 1104
.set nb204_rsqH1H2, 1120
.set nb204_rsqH1M, 1136
.set nb204_rsqH2H1, 1152
.set nb204_rsqH2H2, 1168
.set nb204_rsqH2M, 1184
.set nb204_rsqMH1, 1200
.set nb204_rsqMH2, 1216
.set nb204_rsqMM, 1232
.set nb204_rinvH1H1, 1248
.set nb204_rinvH1H2, 1264
.set nb204_rinvH1M, 1280
.set nb204_rinvH2H1, 1296
.set nb204_rinvH2H2, 1312
.set nb204_rinvH2M, 1328
.set nb204_rinvMH1, 1344
.set nb204_rinvMH2, 1360
.set nb204_rinvMM, 1376
.set nb204_two, 1392
.set nb204_krf, 1408
.set nb204_crf, 1424
.set nb204_is3, 1440
.set nb204_ii3, 1444
.set nb204_innerjjnr, 1448
.set nb204_innerk, 1452
.set nb204_n, 1456
.set nb204_nn1, 1460
.set nb204_nri, 1464
.set nb204_nouter, 1468
.set nb204_ninner, 1472
.set nb204_salign, 1476
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
        movl %eax,nb204_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb204_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb204_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb204_nouter(%esp)
        movl %eax,nb204_ninner(%esp)


        movl nb204_argkrf(%ebp),%esi
        movl nb204_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb204_krf(%esp)
        movaps %xmm6,nb204_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb204_half(%esp)
        movss nb204_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb204_half(%esp)
        movaps %xmm2,nb204_two(%esp)
        movaps %xmm3,nb204_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb204_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb204_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%edx,%ebx,4),%xmm5
        movl nb204_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb204_qqHH(%esp)
        movaps %xmm4,nb204_qqMH(%esp)
        movaps %xmm5,nb204_qqMM(%esp)

_nb_kernel204_ia32_sse.nb204_threadloop: 
        movl  nb204_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel204_ia32_sse.nb204_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel204_ia32_sse.nb204_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb204_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb204_n(%esp)
        movl %ebx,nb204_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel204_ia32_sse.nb204_outerstart
        jmp _nb_kernel204_ia32_sse.nb204_end

_nb_kernel204_ia32_sse.nb204_outerstart: 
        ## ebx contains number of outer iterations
        addl nb204_nouter(%esp),%ebx
        movl %ebx,nb204_nouter(%esp)

_nb_kernel204_ia32_sse.nb204_outer: 
        movl  nb204_shift(%ebp),%eax            ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb204_is3(%esp)      ## store is3 

        movl  nb204_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb204_iinr(%ebp),%ecx             ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb204_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb204_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb204_ixH1(%esp)
        movaps %xmm4,nb204_iyH1(%esp)
        movaps %xmm5,nb204_izH1(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%eax,%ebx,4),%xmm0
        addss 28(%eax,%ebx,4),%xmm1
        addss 32(%eax,%ebx,4),%xmm2
        addss 36(%eax,%ebx,4),%xmm3
        addss 40(%eax,%ebx,4),%xmm4
        addss 44(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb204_ixH2(%esp)
        movaps %xmm1,nb204_iyH2(%esp)
        movaps %xmm2,nb204_izH2(%esp)
        movaps %xmm3,nb204_ixM(%esp)
        movaps %xmm4,nb204_iyM(%esp)
        movaps %xmm5,nb204_izM(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb204_vctot(%esp)
        movaps %xmm4,nb204_fixH1(%esp)
        movaps %xmm4,nb204_fiyH1(%esp)
        movaps %xmm4,nb204_fizH1(%esp)
        movaps %xmm4,nb204_fixH2(%esp)
        movaps %xmm4,nb204_fiyH2(%esp)
        movaps %xmm4,nb204_fizH2(%esp)
        movaps %xmm4,nb204_fixM(%esp)
        movaps %xmm4,nb204_fiyM(%esp)
        movaps %xmm4,nb204_fizM(%esp)

        movl  nb204_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb204_pos(%ebp),%esi
        movl  nb204_faction(%ebp),%edi
        movl  nb204_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb204_innerjjnr(%esp)        ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb204_ninner(%esp),%ecx
        movl  %ecx,nb204_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb204_innerk(%esp)   ## number of innerloop atoms 
        jge   _nb_kernel204_ia32_sse.nb204_unroll_loop
        jmp   _nb_kernel204_ia32_sse.nb204_single_check
_nb_kernel204_ia32_sse.nb204_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb204_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb204_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb204_pos(%ebp),%esi       ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps 12(%esi,%eax,4),%xmm2
        movlps 24(%esi,%eax,4),%xmm3
        movlps 36(%esi,%eax,4),%xmm4

        movlps 12(%esi,%ebx,4),%xmm5
        movlps 24(%esi,%ebx,4),%xmm6
        movlps 36(%esi,%ebx,4),%xmm7

        movhps 12(%esi,%ecx,4),%xmm2
        movhps 24(%esi,%ecx,4),%xmm3
        movhps 36(%esi,%ecx,4),%xmm4

        movhps 12(%esi,%edx,4),%xmm5
        movhps 24(%esi,%edx,4),%xmm6
        movhps 36(%esi,%edx,4),%xmm7

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0
        unpcklps %xmm6,%xmm1
        unpckhps %xmm5,%xmm2
        unpckhps %xmm6,%xmm3
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4
        unpckhps %xmm7,%xmm5
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0
        movaps %xmm0,nb204_jxH1(%esp)
        movhlps  %xmm6,%xmm2
        movaps %xmm2,nb204_jyH1(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb204_jxH2(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb204_jyH2(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb204_jxM(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb204_jyM(%esp)

        movss  20(%esi,%eax,4),%xmm0
        movss  32(%esi,%eax,4),%xmm1
        movss  44(%esi,%eax,4),%xmm2

        movss  20(%esi,%ecx,4),%xmm3
        movss  32(%esi,%ecx,4),%xmm4
        movss  44(%esi,%ecx,4),%xmm5

        movhps 16(%esi,%ebx,4),%xmm0
        movhps 28(%esi,%ebx,4),%xmm1
        movhps 40(%esi,%ebx,4),%xmm2

        movhps 16(%esi,%edx,4),%xmm3
        movhps 28(%esi,%edx,4),%xmm4
        movhps 40(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb204_jzH1(%esp)
        movaps %xmm1,nb204_jzH2(%esp)
        movaps %xmm2,nb204_jzM(%esp)

        movaps nb204_ixH1(%esp),%xmm0
        movaps nb204_iyH1(%esp),%xmm1
        movaps nb204_izH1(%esp),%xmm2
        movaps nb204_ixH1(%esp),%xmm3
        movaps nb204_iyH1(%esp),%xmm4
        movaps nb204_izH1(%esp),%xmm5
        subps  nb204_jxH1(%esp),%xmm0
        subps  nb204_jyH1(%esp),%xmm1
        subps  nb204_jzH1(%esp),%xmm2
        subps  nb204_jxH2(%esp),%xmm3
        subps  nb204_jyH2(%esp),%xmm4
        subps  nb204_jzH2(%esp),%xmm5
        movaps %xmm0,nb204_dxH1H1(%esp)
        movaps %xmm1,nb204_dyH1H1(%esp)
        movaps %xmm2,nb204_dzH1H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb204_dxH1H2(%esp)
        movaps %xmm4,nb204_dyH1H2(%esp)
        movaps %xmm5,nb204_dzH1H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb204_rsqH1H1(%esp)
        movaps %xmm3,nb204_rsqH1H2(%esp)

        movaps nb204_ixH1(%esp),%xmm0
        movaps nb204_iyH1(%esp),%xmm1
        movaps nb204_izH1(%esp),%xmm2
        movaps nb204_ixH2(%esp),%xmm3
        movaps nb204_iyH2(%esp),%xmm4
        movaps nb204_izH2(%esp),%xmm5
        subps  nb204_jxM(%esp),%xmm0
        subps  nb204_jyM(%esp),%xmm1
        subps  nb204_jzM(%esp),%xmm2
        subps  nb204_jxH1(%esp),%xmm3
        subps  nb204_jyH1(%esp),%xmm4
        subps  nb204_jzH1(%esp),%xmm5
        movaps %xmm0,nb204_dxH1M(%esp)
        movaps %xmm1,nb204_dyH1M(%esp)
        movaps %xmm2,nb204_dzH1M(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb204_dxH2H1(%esp)
        movaps %xmm4,nb204_dyH2H1(%esp)
        movaps %xmm5,nb204_dzH2H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb204_rsqH1M(%esp)
        movaps %xmm3,nb204_rsqH2H1(%esp)

        movaps nb204_ixH2(%esp),%xmm0
        movaps nb204_iyH2(%esp),%xmm1
        movaps nb204_izH2(%esp),%xmm2
        movaps nb204_ixH2(%esp),%xmm3
        movaps nb204_iyH2(%esp),%xmm4
        movaps nb204_izH2(%esp),%xmm5
        subps  nb204_jxH2(%esp),%xmm0
        subps  nb204_jyH2(%esp),%xmm1
        subps  nb204_jzH2(%esp),%xmm2
        subps  nb204_jxM(%esp),%xmm3
        subps  nb204_jyM(%esp),%xmm4
        subps  nb204_jzM(%esp),%xmm5
        movaps %xmm0,nb204_dxH2H2(%esp)
        movaps %xmm1,nb204_dyH2H2(%esp)
        movaps %xmm2,nb204_dzH2H2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb204_dxH2M(%esp)
        movaps %xmm4,nb204_dyH2M(%esp)
        movaps %xmm5,nb204_dzH2M(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb204_rsqH2H2(%esp)
        movaps %xmm3,nb204_rsqH2M(%esp)

        movaps nb204_ixM(%esp),%xmm0
        movaps nb204_iyM(%esp),%xmm1
        movaps nb204_izM(%esp),%xmm2
        movaps nb204_ixM(%esp),%xmm3
        movaps nb204_iyM(%esp),%xmm4
        movaps nb204_izM(%esp),%xmm5
        subps  nb204_jxH1(%esp),%xmm0
        subps  nb204_jyH1(%esp),%xmm1
        subps  nb204_jzH1(%esp),%xmm2
        subps  nb204_jxH2(%esp),%xmm3
        subps  nb204_jyH2(%esp),%xmm4
        subps  nb204_jzH2(%esp),%xmm5
        movaps %xmm0,nb204_dxMH1(%esp)
        movaps %xmm1,nb204_dyMH1(%esp)
        movaps %xmm2,nb204_dzMH1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb204_dxMH2(%esp)
        movaps %xmm4,nb204_dyMH2(%esp)
        movaps %xmm5,nb204_dzMH2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb204_rsqMH1(%esp)
        movaps %xmm4,nb204_rsqMH2(%esp)

        movaps nb204_ixM(%esp),%xmm0
        movaps nb204_iyM(%esp),%xmm1
        movaps nb204_izM(%esp),%xmm2
        subps  nb204_jxM(%esp),%xmm0
        subps  nb204_jyM(%esp),%xmm1
        subps  nb204_jzM(%esp),%xmm2
        movaps %xmm0,nb204_dxMM(%esp)
        movaps %xmm1,nb204_dyMM(%esp)
        movaps %xmm2,nb204_dzMM(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb204_rsqMM(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204_half(%esp),%xmm3   ## rinvH2H2 
        mulps   nb204_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb204_rinvMM(%esp)
        movaps  %xmm7,nb204_rinvMH2(%esp)

        rsqrtps nb204_rsqH1H1(%esp),%xmm1
        rsqrtps nb204_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204_rsqH1H1(%esp),%xmm1
        mulps   nb204_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204_half(%esp),%xmm3
        mulps   nb204_half(%esp),%xmm7
        movaps  %xmm3,nb204_rinvH1H1(%esp)
        movaps  %xmm7,nb204_rinvH1H2(%esp)

        rsqrtps nb204_rsqH1M(%esp),%xmm1
        rsqrtps nb204_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204_rsqH1M(%esp),%xmm1
        mulps   nb204_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204_half(%esp),%xmm3
        mulps   nb204_half(%esp),%xmm7
        movaps  %xmm3,nb204_rinvH1M(%esp)
        movaps  %xmm7,nb204_rinvH2H1(%esp)

        rsqrtps nb204_rsqH2H2(%esp),%xmm1
        rsqrtps nb204_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204_rsqH2H2(%esp),%xmm1
        mulps   nb204_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204_half(%esp),%xmm3
        mulps   nb204_half(%esp),%xmm7
        movaps  %xmm3,nb204_rinvH2H2(%esp)
        movaps  %xmm7,nb204_rinvH2M(%esp)

        rsqrtps nb204_rsqMH1(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb204_three(%esp),%xmm3
        mulps   nb204_rsqMH1(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb204_half(%esp),%xmm3
        movaps  %xmm3,nb204_rinvMH1(%esp)

        ## start with H1-H1 interaction 
        movaps nb204_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb204_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulps  nb204_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb204_crf(%esp),%xmm6
        mulps  nb204_qqHH(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulps nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addps  nb204_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulps  %xmm7,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb204_dxH1H1(%esp),%xmm0
        mulps nb204_dyH1H1(%esp),%xmm1
        mulps nb204_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixH1(%esp),%xmm0
        addps nb204_fiyH1(%esp),%xmm1
        addps nb204_fizH1(%esp),%xmm2
        movaps %xmm3,nb204_fjxH1(%esp)
        movaps %xmm4,nb204_fjyH1(%esp)
        movaps %xmm5,nb204_fjzH1(%esp)
        movaps %xmm0,nb204_fixH1(%esp)
        movaps %xmm1,nb204_fiyH1(%esp)
        movaps %xmm2,nb204_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb204_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb204_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb204_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb204_crf(%esp),%xmm4
        mulps  %xmm0,%xmm0
        mulps  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH1  
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb204_dxH1H2(%esp),%xmm0
        mulps nb204_dyH1H2(%esp),%xmm1
        mulps nb204_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixH1(%esp),%xmm0
        addps nb204_fiyH1(%esp),%xmm1
        addps nb204_fizH1(%esp),%xmm2
        movaps %xmm3,nb204_fjxH2(%esp)
        movaps %xmm4,nb204_fjyH2(%esp)
        movaps %xmm5,nb204_fjzH2(%esp)
        movaps %xmm0,nb204_fixH1(%esp)
        movaps %xmm1,nb204_fiyH1(%esp)
        movaps %xmm2,nb204_fizH1(%esp)

        ## H1-M interaction  
        movaps nb204_rinvH1M(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=Rinv 
        movaps nb204_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb204_rsqH1M(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb204_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb204_dxH1M(%esp),%xmm0
        mulps nb204_dyH1M(%esp),%xmm1
        mulps nb204_dzH1M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixH1(%esp),%xmm0
        addps nb204_fiyH1(%esp),%xmm1
        addps nb204_fizH1(%esp),%xmm2
        movaps %xmm3,nb204_fjxM(%esp)
        movaps %xmm4,nb204_fjyM(%esp)
        movaps %xmm5,nb204_fjzM(%esp)
        movaps %xmm0,nb204_fixH1(%esp)
        movaps %xmm1,nb204_fiyH1(%esp)
        movaps %xmm2,nb204_fizH1(%esp)

        ## H2-H1 interaction 
        movaps nb204_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb204_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb204_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb204_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb204_fjxH1(%esp),%xmm3
        movaps nb204_fjyH1(%esp),%xmm4
        movaps nb204_fjzH1(%esp),%xmm5
        mulps nb204_dxH2H1(%esp),%xmm0
        mulps nb204_dyH2H1(%esp),%xmm1
        mulps nb204_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixH2(%esp),%xmm0
        addps nb204_fiyH2(%esp),%xmm1
        addps nb204_fizH2(%esp),%xmm2
        movaps %xmm3,nb204_fjxH1(%esp)
        movaps %xmm4,nb204_fjyH1(%esp)
        movaps %xmm5,nb204_fjzH1(%esp)
        movaps %xmm0,nb204_fixH2(%esp)
        movaps %xmm1,nb204_fiyH2(%esp)
        movaps %xmm2,nb204_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb204_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb204_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb204_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb204_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb204_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb204_fjxH2(%esp),%xmm3
        movaps nb204_fjyH2(%esp),%xmm4
        movaps nb204_fjzH2(%esp),%xmm5
        mulps nb204_dxH2H2(%esp),%xmm0
        mulps nb204_dyH2H2(%esp),%xmm1
        mulps nb204_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixH2(%esp),%xmm0
        addps nb204_fiyH2(%esp),%xmm1
        addps nb204_fizH2(%esp),%xmm2
        movaps %xmm3,nb204_fjxH2(%esp)
        movaps %xmm4,nb204_fjyH2(%esp)
        movaps %xmm5,nb204_fjzH2(%esp)
        movaps %xmm0,nb204_fixH2(%esp)
        movaps %xmm1,nb204_fiyH2(%esp)
        movaps %xmm2,nb204_fizH2(%esp)

        ## H2-M interaction 
        movaps nb204_rinvH2M(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb204_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb204_rsqH2M(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb204_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb204_fjxM(%esp),%xmm3
        movaps nb204_fjyM(%esp),%xmm4
        movaps nb204_fjzM(%esp),%xmm5
        mulps nb204_dxH2M(%esp),%xmm0
        mulps nb204_dyH2M(%esp),%xmm1
        mulps nb204_dzH2M(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixH2(%esp),%xmm0
        addps nb204_fiyH2(%esp),%xmm1
        addps nb204_fizH2(%esp),%xmm2
        movaps %xmm3,nb204_fjxM(%esp)
        movaps %xmm4,nb204_fjyM(%esp)
        movaps %xmm5,nb204_fjzM(%esp)
        movaps %xmm0,nb204_fixH2(%esp)
        movaps %xmm1,nb204_fiyH2(%esp)
        movaps %xmm2,nb204_fizH2(%esp)

        ## M-H1 interaction 
        movaps nb204_rinvMH1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb204_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb204_rsqMH1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb204_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb204_fjxH1(%esp),%xmm3
        movaps nb204_fjyH1(%esp),%xmm4
        movaps nb204_fjzH1(%esp),%xmm5
        mulps nb204_dxMH1(%esp),%xmm0
        mulps nb204_dyMH1(%esp),%xmm1
        mulps nb204_dzMH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixM(%esp),%xmm0
        addps nb204_fiyM(%esp),%xmm1
        addps nb204_fizM(%esp),%xmm2
        movaps %xmm3,nb204_fjxH1(%esp)
        movaps %xmm4,nb204_fjyH1(%esp)
        movaps %xmm5,nb204_fjzH1(%esp)
        movaps %xmm0,nb204_fixM(%esp)
        movaps %xmm1,nb204_fiyM(%esp)
        movaps %xmm2,nb204_fizM(%esp)

        ## M-H2 interaction 
        movaps nb204_rinvMH2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb204_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb204_rsqMH2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb204_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb204_qqMH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqMH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb204_fjxH2(%esp),%xmm3
        movaps nb204_fjyH2(%esp),%xmm4
        movaps nb204_fjzH2(%esp),%xmm5
        mulps nb204_dxMH2(%esp),%xmm0
        mulps nb204_dyMH2(%esp),%xmm1
        mulps nb204_dzMH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixM(%esp),%xmm0
        addps nb204_fiyM(%esp),%xmm1
        addps nb204_fizM(%esp),%xmm2
        movaps %xmm3,nb204_fjxH2(%esp)
        movaps %xmm4,nb204_fjyH2(%esp)
        movaps %xmm5,nb204_fjzH2(%esp)
        movaps %xmm0,nb204_fixM(%esp)
        movaps %xmm1,nb204_fiyM(%esp)
        movaps %xmm2,nb204_fizM(%esp)

        ## M-M interaction 
        movaps nb204_rinvMM(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb204_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb204_rsqMM(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb204_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb204_qqMM(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb204_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb204_qqMM(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps %xmm0,%xmm1
        movaps %xmm6,nb204_vctot(%esp)
        movaps %xmm0,%xmm2

        movaps nb204_fjxM(%esp),%xmm3
        movaps nb204_fjyM(%esp),%xmm4
        movaps nb204_fjzM(%esp),%xmm5
        mulps nb204_dxMM(%esp),%xmm0
        mulps nb204_dyMM(%esp),%xmm1
        mulps nb204_dzMM(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb204_fixM(%esp),%xmm0
        addps nb204_fiyM(%esp),%xmm1
        addps nb204_fizM(%esp),%xmm2
        movaps %xmm3,nb204_fjxM(%esp)
        movaps %xmm4,nb204_fjyM(%esp)
        movaps %xmm5,nb204_fjzM(%esp)
        movaps %xmm0,nb204_fixM(%esp)
        movaps %xmm1,nb204_fiyM(%esp)
        movaps %xmm2,nb204_fizM(%esp)

        movl nb204_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## At this stage forces are still on the stack, in positions:
        ## fjxH1, fjyH1, fjzH1, ... , fjzM.
        ## Each position is a quadruplet of forces for the four 
        ## corresponding j waters, so we need to transpose them before
        ## adding to the memory positions.
        ## 
        ## This _used_ to be a simple transpose, but the resulting high number
        ## of unaligned 128-bit load/stores might trigger a possible hardware 
        ## bug on Athlon and Opteron chips, so I have worked around it
        ## to use 64-bit load/stores instead. The performance hit should be
        ## very modest, since the 128-bit unaligned memory instructions were
        ## slow anyway. 

        ## 4 j waters with three atoms each - first do 1st Hydrogen X & Y forces for 4 j particles 
        movaps nb204_fjxH1(%esp),%xmm0   ## xmm0= fjxH1a  fjxH1b  fjxH1c  fjxH1d 
        movaps nb204_fjyH1(%esp),%xmm2   ## xmm1= fjyH1a  fjyH1b  fjyH1c  fjyH1d
        movlps 12(%edi,%eax,4),%xmm3
        movlps 12(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjxH1a  fjyH1a  fjxH1b  fjyH1b
        unpckhps %xmm2,%xmm1       ## xmm1= fjxH1c  fjyH1c  fjxH1d  fjyH1d
        movhps 12(%edi,%ebx,4),%xmm3
        movhps 12(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,12(%edi,%eax,4)
        movlps %xmm4,12(%edi,%ecx,4)
        movhps %xmm3,12(%edi,%ebx,4)
        movhps %xmm4,12(%edi,%edx,4)

        ## constant 1st Hydrogen Z & 2nd hydrogen X forces for 4 j particles 
        movaps nb204_fjzH1(%esp),%xmm0    ## xmm0= fjzH1a   fjzH1b   fjzH1c   fjzH1d 
        movaps nb204_fjxH2(%esp),%xmm2   ## xmm1= fjxH2a  fjxH2b  fjxH2c  fjxH2d
        movlps 20(%edi,%eax,4),%xmm3
        movlps 20(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjzH1a  fjxH2a  fjzH1b  fjxH2b
        unpckhps %xmm2,%xmm1       ## xmm1= fjzH1c  fjxH2c  fjzH1d  fjxH2d
        movhps 20(%edi,%ebx,4),%xmm3
        movhps 20(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,20(%edi,%eax,4)
        movlps %xmm4,20(%edi,%ecx,4)
        movhps %xmm3,20(%edi,%ebx,4)
        movhps %xmm4,20(%edi,%edx,4)

        ## constant 2nd hydrogen Y & Z forces for 4 j particles 
        movaps nb204_fjyH2(%esp),%xmm0   ## xmm0= fjyH2a  fjyH2b  fjyH2c  fjyH2d 
        movaps nb204_fjzH2(%esp),%xmm2   ## xmm1= fjzH2a  fjzH2b  fjzH2c  fjzH2d
        movlps 28(%edi,%eax,4),%xmm3
        movlps 28(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjyH2a  fjzH2a  fjyH2b  fjzH2b
        unpckhps %xmm2,%xmm1            ## xmm1= fjyH2c  fjzH2c  fjyH2d  fjzH2d
        movhps 28(%edi,%ebx,4),%xmm3
        movhps 28(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,28(%edi,%eax,4)
        movlps %xmm4,28(%edi,%ecx,4)
        movhps %xmm3,28(%edi,%ebx,4)
        movhps %xmm4,28(%edi,%edx,4)

        ## Dummy (M) X & Y forces for 4 j particles 
        movaps nb204_fjxM(%esp),%xmm0   ## xmm0= fjxMa  fjxMb  fjxMc  fjxMd 
        movaps nb204_fjyM(%esp),%xmm2   ## xmm1= fjyMa  fjyMb  fjyMc  fjyMd
        movlps 36(%edi,%eax,4),%xmm3
        movlps 36(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjxMa  fjyMa  fjxMb  fjyMb
        unpckhps %xmm2,%xmm1            ## xmm1= fjxMc  fjyMc  fjxMd  fjyMd
        movhps 36(%edi,%ebx,4),%xmm3
        movhps 36(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,36(%edi,%eax,4)
        movlps %xmm4,36(%edi,%ecx,4)
        movhps %xmm3,36(%edi,%ebx,4)
        movhps %xmm4,36(%edi,%edx,4)


        ## Dummy (M) Z forces for 4 j particles 
        ## Just load the four Z coords into one reg. each
        movss 44(%edi,%eax,4),%xmm4
        movss 44(%edi,%ebx,4),%xmm5
        movss 44(%edi,%ecx,4),%xmm6
        movss 44(%edi,%edx,4),%xmm7
        ## add what we have on the stack
        addss nb204_fjzM(%esp),%xmm4
        addss nb204_fjzMb(%esp),%xmm5
        addss nb204_fjzMc(%esp),%xmm6
        addss nb204_fjzMd(%esp),%xmm7
        ## store back
        movss %xmm4,44(%edi,%eax,4)
        movss %xmm5,44(%edi,%ebx,4)
        movss %xmm6,44(%edi,%ecx,4)
        movss %xmm7,44(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb204_innerk(%esp)
        jl    _nb_kernel204_ia32_sse.nb204_single_check
        jmp   _nb_kernel204_ia32_sse.nb204_unroll_loop
_nb_kernel204_ia32_sse.nb204_single_check: 
        addl $4,nb204_innerk(%esp)
        jnz   _nb_kernel204_ia32_sse.nb204_single_loop
        jmp   _nb_kernel204_ia32_sse.nb204_updateouterdata
_nb_kernel204_ia32_sse.nb204_single_loop: 
        movl  nb204_innerjjnr(%esp),%edx        ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb204_innerjjnr(%esp)

        movl nb204_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5
        movss 36(%esi,%eax,4),%xmm3             ## jxM  -  -  -
        movss 40(%esi,%eax,4),%xmm4             ## jyM  -  -  -
        movss 44(%esi,%eax,4),%xmm5             ## jzM  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb204_ixM(%esp),%xmm0
        movaps  nb204_iyM(%esp),%xmm1
        movaps  nb204_izM(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM
        movaps %xmm3,nb204_jxM(%esp)
        movaps %xmm4,nb204_jyM(%esp)
        movaps %xmm5,nb204_jzM(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        movaps %xmm0,nb204_dxMM(%esp)
        movaps %xmm1,nb204_dyMM(%esp)
        movaps %xmm2,nb204_dzMM(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        mulps   nb204_krf(%esp),%xmm6   ## xmm6=krsq 
        movaps  %xmm1,%xmm2
        movaps  %xmm6,%xmm7     ## xmm7=krsq 
        mulps   %xmm1,%xmm1
        movaps  nb204_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb204_half(%esp),%xmm3   ## rinv iO - j water 

        addps   %xmm3,%xmm6     ## xmm6=rinv+ krsq 
        mulps   nb204_two(%esp),%xmm7
        subps   nb204_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        subps   %xmm7,%xmm3     ## xmm3=rinv-2*krsq 
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   nb204_qqMM(%esp),%xmm4
        movhps  nb204_qqMH(%esp),%xmm4

        mulps %xmm4,%xmm6       ## vcoul  
        mulps %xmm4,%xmm3       ## coul part of fs  


        addps   nb204_vctot(%esp),%xmm6
        mulps   %xmm3,%xmm0     ## total fscal 
        movaps  %xmm6,nb204_vctot(%esp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb204_dxMM(%esp),%xmm0
        mulps   nb204_dyMM(%esp),%xmm1
        mulps   nb204_dzMM(%esp),%xmm2

        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb204_fjxM(%esp)
        movaps  %xmm4,nb204_fjyM(%esp)
        movaps  %xmm5,nb204_fjzM(%esp)
        addps   nb204_fixM(%esp),%xmm0
        addps   nb204_fiyM(%esp),%xmm1
        addps   nb204_fizM(%esp),%xmm2
        movaps  %xmm0,nb204_fixM(%esp)
        movaps  %xmm1,nb204_fiyM(%esp)
        movaps  %xmm2,nb204_fizM(%esp)


        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb204_ixH1(%esp),%xmm0
        movaps  nb204_iyH1(%esp),%xmm1
        movaps  nb204_izH1(%esp),%xmm2
        movaps  nb204_ixH2(%esp),%xmm3
        movaps  nb204_iyH2(%esp),%xmm4
        movaps  nb204_izH2(%esp),%xmm5
        subps   nb204_jxM(%esp),%xmm0
        subps   nb204_jyM(%esp),%xmm1
        subps   nb204_jzM(%esp),%xmm2
        subps   nb204_jxM(%esp),%xmm3
        subps   nb204_jyM(%esp),%xmm4
        subps   nb204_jzM(%esp),%xmm5
        movaps %xmm0,nb204_dxH1M(%esp)
        movaps %xmm1,nb204_dyH1M(%esp)
        movaps %xmm2,nb204_dzH1M(%esp)
        movaps %xmm3,nb204_dxH2M(%esp)
        movaps %xmm4,nb204_dyH2M(%esp)
        movaps %xmm5,nb204_dzH2M(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb204_half(%esp),%xmm7   ## rinv H2 - j water  

        mulps nb204_krf(%esp),%xmm0   ## krsq 
        mulps nb204_krf(%esp),%xmm4   ## krsq  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb204_qqMH(%esp),%xmm6
        movhps  nb204_qqHH(%esp),%xmm6
        movaps  %xmm0,%xmm1
        movaps  %xmm4,%xmm5
        addps   %xmm3,%xmm0     ## krsq+ rinv 
        addps   %xmm7,%xmm4     ## krsq+ rinv 
        subps   nb204_crf(%esp),%xmm0
        subps   nb204_crf(%esp),%xmm4
        mulps   nb204_two(%esp),%xmm1
        mulps   nb204_two(%esp),%xmm5
        mulps   %xmm6,%xmm0     ## vcoul 
        mulps   %xmm6,%xmm4     ## vcoul 
        addps   %xmm0,%xmm4
        addps   nb204_vctot(%esp),%xmm4
        movaps  %xmm4,nb204_vctot(%esp)
        movaps  %xmm3,%xmm0
        movaps  %xmm7,%xmm4
        mulps   %xmm3,%xmm3
        mulps   %xmm7,%xmm7
        subps   %xmm1,%xmm0
        subps   %xmm5,%xmm4
        mulps   %xmm6,%xmm0
        mulps   %xmm6,%xmm4
        mulps   %xmm3,%xmm0     ## fscal 
        mulps   %xmm4,%xmm7     ## fscal 

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb204_dxH1M(%esp),%xmm0
        mulps   nb204_dyH1M(%esp),%xmm1
        mulps   nb204_dzH1M(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb204_fjxM(%esp),%xmm3
        movaps  nb204_fjyM(%esp),%xmm4
        movaps  nb204_fjzM(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb204_fjxM(%esp)
        movaps  %xmm4,nb204_fjyM(%esp)
        movaps  %xmm5,nb204_fjzM(%esp)
        addps   nb204_fixH1(%esp),%xmm0
        addps   nb204_fiyH1(%esp),%xmm1
        addps   nb204_fizH1(%esp),%xmm2
        movaps  %xmm0,nb204_fixH1(%esp)
        movaps  %xmm1,nb204_fiyH1(%esp)
        movaps  %xmm2,nb204_fizH1(%esp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb204_dxH2M(%esp),%xmm0
        mulps   nb204_dyH2M(%esp),%xmm1
        mulps   nb204_dzH2M(%esp),%xmm2
        movaps  nb204_fjxM(%esp),%xmm3
        movaps  nb204_fjyM(%esp),%xmm4
        movaps  nb204_fjzM(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    nb204_faction(%ebp),%esi
        movaps  %xmm3,nb204_fjxM(%esp)
        movaps  %xmm4,nb204_fjyM(%esp)
        movaps  %xmm5,nb204_fjzM(%esp)
        addps   nb204_fixH2(%esp),%xmm0
        addps   nb204_fiyH2(%esp),%xmm1
        addps   nb204_fizH2(%esp),%xmm2
        movaps  %xmm0,nb204_fixH2(%esp)
        movaps  %xmm1,nb204_fiyH2(%esp)
        movaps  %xmm2,nb204_fizH2(%esp)

        ## update j water forces from local variables 
        movlps  36(%esi,%eax,4),%xmm0
        movlps  12(%esi,%eax,4),%xmm1
        movhps  24(%esi,%eax,4),%xmm1
        movaps  nb204_fjxM(%esp),%xmm3
        movaps  nb204_fjyM(%esp),%xmm4
        movaps  nb204_fjzM(%esp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## constant 00000010
        shufps $3,%xmm7,%xmm7 ## constant 00000011
        addss   44(%esi,%eax,4),%xmm5
        addss   20(%esi,%eax,4),%xmm6
        addss   32(%esi,%eax,4),%xmm7
        movss   %xmm5,44(%esi,%eax,4)
        movss   %xmm6,20(%esi,%eax,4)
        movss   %xmm7,32(%esi,%eax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,36(%esi,%eax,4)
        movlps  %xmm1,12(%esi,%eax,4)
        movhps  %xmm1,24(%esi,%eax,4)

        decl nb204_innerk(%esp)
        jz    _nb_kernel204_ia32_sse.nb204_updateouterdata
        jmp   _nb_kernel204_ia32_sse.nb204_single_loop
_nb_kernel204_ia32_sse.nb204_updateouterdata: 
        movl  nb204_ii3(%esp),%ecx
        movl  nb204_faction(%ebp),%edi
        movl  nb204_fshift(%ebp),%esi
        movl  nb204_is3(%esp),%edx

        ## accumulate  H1 i forces in xmm0, xmm1, xmm2 
        movaps nb204_fixH1(%esp),%xmm0
        movaps nb204_fiyH1(%esp),%xmm1
        movaps nb204_fizH1(%esp),%xmm2

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
        movss  12(%edi,%ecx,4),%xmm3
        movss  16(%edi,%ecx,4),%xmm4
        movss  20(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,12(%edi,%ecx,4)
        movss  %xmm4,16(%edi,%ecx,4)
        movss  %xmm5,20(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb204_fixH2(%esp),%xmm0
        movaps nb204_fiyH2(%esp),%xmm1
        movaps nb204_fizH2(%esp),%xmm2

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
        movss  24(%edi,%ecx,4),%xmm3
        movss  28(%edi,%ecx,4),%xmm4
        movss  32(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,24(%edi,%ecx,4)
        movss  %xmm4,28(%edi,%ecx,4)
        movss  %xmm5,32(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate M i forces in xmm0, xmm1, xmm2 
        movaps nb204_fixM(%esp),%xmm0
        movaps nb204_fiyM(%esp),%xmm1
        movaps nb204_fizM(%esp),%xmm2

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
        movss  36(%edi,%ecx,4),%xmm3
        movss  40(%edi,%ecx,4),%xmm4
        movss  44(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,36(%edi,%ecx,4)
        movss  %xmm4,40(%edi,%ecx,4)
        movss  %xmm5,44(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb204_n(%esp),%esi
        ## get group index for i particle 
        movl  nb204_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb204_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb204_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb204_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel204_ia32_sse.nb204_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb204_n(%esp)
        jmp _nb_kernel204_ia32_sse.nb204_outer
_nb_kernel204_ia32_sse.nb204_outerend: 
        ## check if more outer neighborlists remain
        movl  nb204_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel204_ia32_sse.nb204_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel204_ia32_sse.nb204_threadloop
_nb_kernel204_ia32_sse.nb204_end: 
        emms

        movl nb204_nouter(%esp),%eax
        movl nb204_ninner(%esp),%ebx
        movl nb204_outeriter(%ebp),%ecx
        movl nb204_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb204_salign(%esp),%eax
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




.globl nb_kernel204nf_ia32_sse
.globl _nb_kernel204nf_ia32_sse
nb_kernel204nf_ia32_sse:        
_nb_kernel204nf_ia32_sse:       
.set nb204nf_p_nri, 8
.set nb204nf_iinr, 12
.set nb204nf_jindex, 16
.set nb204nf_jjnr, 20
.set nb204nf_shift, 24
.set nb204nf_shiftvec, 28
.set nb204nf_fshift, 32
.set nb204nf_gid, 36
.set nb204nf_pos, 40
.set nb204nf_faction, 44
.set nb204nf_charge, 48
.set nb204nf_p_facel, 52
.set nb204nf_argkrf, 56
.set nb204nf_argcrf, 60
.set nb204nf_Vc, 64
.set nb204nf_type, 68
.set nb204nf_p_ntype, 72
.set nb204nf_vdwparam, 76
.set nb204nf_Vvdw, 80
.set nb204nf_p_tabscale, 84
.set nb204nf_VFtab, 88
.set nb204nf_invsqrta, 92
.set nb204nf_dvda, 96
.set nb204nf_p_gbtabscale, 100
.set nb204nf_GBtab, 104
.set nb204nf_p_nthreads, 108
.set nb204nf_count, 112
.set nb204nf_mtx, 116
.set nb204nf_outeriter, 120
.set nb204nf_inneriter, 124
.set nb204nf_work, 128
        ## stack offsets for local variables  
        ## bottom of stack is cache-aligned for sse use 
.set nb204nf_ixH1, 0
.set nb204nf_iyH1, 16
.set nb204nf_izH1, 32
.set nb204nf_ixH2, 48
.set nb204nf_iyH2, 64
.set nb204nf_izH2, 80
.set nb204nf_ixM, 96
.set nb204nf_iyM, 112
.set nb204nf_izM, 128
.set nb204nf_jxH1, 144
.set nb204nf_jyH1, 160
.set nb204nf_jzH1, 176
.set nb204nf_jxH2, 192
.set nb204nf_jyH2, 208
.set nb204nf_jzH2, 224
.set nb204nf_jxM, 240
.set nb204nf_jyM, 256
.set nb204nf_jzM, 272
.set nb204nf_qqHH, 288
.set nb204nf_qqMH, 304
.set nb204nf_qqMM, 320
.set nb204nf_vctot, 336
.set nb204nf_half, 352
.set nb204nf_three, 368
.set nb204nf_rsqH1H1, 384
.set nb204nf_rsqH1H2, 400
.set nb204nf_rsqH1M, 416
.set nb204nf_rsqH2H1, 432
.set nb204nf_rsqH2H2, 448
.set nb204nf_rsqH2M, 464
.set nb204nf_rsqMH1, 480
.set nb204nf_rsqMH2, 496
.set nb204nf_rsqMM, 512
.set nb204nf_rinvH1H1, 528
.set nb204nf_rinvH1H2, 544
.set nb204nf_rinvH1M, 560
.set nb204nf_rinvH2H1, 576
.set nb204nf_rinvH2H2, 592
.set nb204nf_rinvH2M, 608
.set nb204nf_rinvMH1, 624
.set nb204nf_rinvMH2, 640
.set nb204nf_rinvMM, 656
.set nb204nf_krf, 672
.set nb204nf_crf, 688
.set nb204nf_is3, 704
.set nb204nf_ii3, 708
.set nb204nf_innerjjnr, 712
.set nb204nf_innerk, 716
.set nb204nf_n, 720
.set nb204nf_nn1, 724
.set nb204nf_nri, 728
.set nb204nf_nouter, 732
.set nb204nf_ninner, 736
.set nb204nf_salign, 740
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
        movl %eax,nb204nf_salign(%esp)

        emms

        ## Move args passed by reference to stack
        movl nb204nf_p_nri(%ebp),%ecx
        movl (%ecx),%ecx
        movl %ecx,nb204nf_nri(%esp)

        ## zero iteration counters
        movl $0,%eax
        movl %eax,nb204nf_nouter(%esp)
        movl %eax,nb204nf_ninner(%esp)


        movl nb204nf_argkrf(%ebp),%esi
        movl nb204nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb204nf_krf(%esp)
        movaps %xmm6,nb204nf_crf(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb204nf_half(%esp)
        movss nb204nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb204nf_half(%esp)
        movaps %xmm3,nb204nf_three(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb204nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx),%ebx               ## ebx =ii 

        movl  nb204nf_charge(%ebp),%edx
        movss 4(%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 12(%edx,%ebx,4),%xmm5
        movl nb204nf_p_facel(%ebp),%esi
        movss (%esi),%xmm6
        mulss  %xmm3,%xmm3
        mulss  %xmm5,%xmm4
        mulss  %xmm5,%xmm5
        mulss  %xmm6,%xmm3
        mulss  %xmm6,%xmm4
        mulss  %xmm6,%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb204nf_qqHH(%esp)
        movaps %xmm4,nb204nf_qqMH(%esp)
        movaps %xmm5,nb204nf_qqMM(%esp)

_nb_kernel204nf_ia32_sse.nb204nf_threadloop: 
        movl  nb204nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel204nf_ia32_sse.nb204nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel204nf_ia32_sse.nb204nf_spinlock

        ## if(nn1>nri) nn1=nri
        movl nb204nf_nri(%esp),%ecx
        movl %ecx,%edx
        subl %ebx,%ecx
        cmovlel %edx,%ebx                       ## if(nn1>nri) nn1=nri
        ## Cleared the spinlock if we got here.
        ## eax contains nn0, ebx contains nn1.
        movl %eax,nb204nf_n(%esp)
        movl %ebx,nb204nf_nn1(%esp)
        subl %eax,%ebx                          ## calc number of outer lists
        movl %eax,%esi                          ## copy n to esi
        jg  _nb_kernel204nf_ia32_sse.nb204nf_outerstart
        jmp _nb_kernel204nf_ia32_sse.nb204nf_end

_nb_kernel204nf_ia32_sse.nb204nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb204nf_nouter(%esp),%ebx
        movl %ebx,nb204nf_nouter(%esp)

_nb_kernel204nf_ia32_sse.nb204nf_outer: 
        movl  nb204nf_shift(%ebp),%eax          ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx        ## ebx=3*is 
        movl  %ebx,nb204nf_is3(%esp)            ## store is3 

        movl  nb204nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb204nf_iinr(%ebp),%ecx           ## ecx = pointer into iinr[]    
        movl  (%ecx,%esi,4),%ebx                ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb204nf_pos(%ebp),%eax    ## eax = base of pos[]  
        movl  %ebx,nb204nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm3
        addss 16(%eax,%ebx,4),%xmm4
        addss 20(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb204nf_ixH1(%esp)
        movaps %xmm4,nb204nf_iyH1(%esp)
        movaps %xmm5,nb204nf_izH1(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 24(%eax,%ebx,4),%xmm0
        addss 28(%eax,%ebx,4),%xmm1
        addss 32(%eax,%ebx,4),%xmm2
        addss 36(%eax,%ebx,4),%xmm3
        addss 40(%eax,%ebx,4),%xmm4
        addss 44(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb204nf_ixH2(%esp)
        movaps %xmm1,nb204nf_iyH2(%esp)
        movaps %xmm2,nb204nf_izH2(%esp)
        movaps %xmm3,nb204nf_ixM(%esp)
        movaps %xmm4,nb204nf_iyM(%esp)
        movaps %xmm5,nb204nf_izM(%esp)

        ## clear vctot 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb204nf_vctot(%esp)

        movl  nb204nf_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                 ## number of innerloop atoms 

        movl  nb204nf_pos(%ebp),%esi
        movl  nb204nf_faction(%ebp),%edi
        movl  nb204nf_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb204nf_innerjjnr(%esp)      ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb204nf_ninner(%esp),%ecx
        movl  %ecx,nb204nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb204nf_innerk(%esp)         ## number of innerloop atoms 
        jge   _nb_kernel204nf_ia32_sse.nb204nf_unroll_loop
        jmp   _nb_kernel204nf_ia32_sse.nb204nf_single_check
_nb_kernel204nf_ia32_sse.nb204nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb204nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx             ## eax-edx=jnr1-4 

        addl $16,nb204nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb204nf_pos(%ebp),%esi     ## base of pos[] 

        leal  (%eax,%eax,2),%eax        ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx        ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps 12(%esi,%eax,4),%xmm2
        movlps 24(%esi,%eax,4),%xmm3
        movlps 36(%esi,%eax,4),%xmm4

        movlps 12(%esi,%ebx,4),%xmm5
        movlps 24(%esi,%ebx,4),%xmm6
        movlps 36(%esi,%ebx,4),%xmm7

        movhps 12(%esi,%ecx,4),%xmm2
        movhps 24(%esi,%ecx,4),%xmm3
        movhps 36(%esi,%ecx,4),%xmm4

        movhps 12(%esi,%edx,4),%xmm5
        movhps 24(%esi,%edx,4),%xmm6
        movhps 36(%esi,%edx,4),%xmm7

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0
        unpcklps %xmm6,%xmm1
        unpckhps %xmm5,%xmm2
        unpckhps %xmm6,%xmm3
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4
        unpckhps %xmm7,%xmm5
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0
        movaps %xmm0,nb204nf_jxH1(%esp)
        movhlps  %xmm6,%xmm2
        movaps %xmm2,nb204nf_jyH1(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb204nf_jxH2(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb204nf_jyH2(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb204nf_jxM(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb204nf_jyM(%esp)

        movss  20(%esi,%eax,4),%xmm0
        movss  32(%esi,%eax,4),%xmm1
        movss  44(%esi,%eax,4),%xmm2

        movss  20(%esi,%ecx,4),%xmm3
        movss  32(%esi,%ecx,4),%xmm4
        movss  44(%esi,%ecx,4),%xmm5

        movhps 16(%esi,%ebx,4),%xmm0
        movhps 28(%esi,%ebx,4),%xmm1
        movhps 40(%esi,%ebx,4),%xmm2

        movhps 16(%esi,%edx,4),%xmm3
        movhps 28(%esi,%edx,4),%xmm4
        movhps 40(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb204nf_jzH1(%esp)
        movaps %xmm1,nb204nf_jzH2(%esp)
        movaps %xmm2,nb204nf_jzM(%esp)

        movaps nb204nf_ixH1(%esp),%xmm0
        movaps nb204nf_iyH1(%esp),%xmm1
        movaps nb204nf_izH1(%esp),%xmm2
        movaps nb204nf_ixH1(%esp),%xmm3
        movaps nb204nf_iyH1(%esp),%xmm4
        movaps nb204nf_izH1(%esp),%xmm5
        subps  nb204nf_jxH1(%esp),%xmm0
        subps  nb204nf_jyH1(%esp),%xmm1
        subps  nb204nf_jzH1(%esp),%xmm2
        subps  nb204nf_jxH2(%esp),%xmm3
        subps  nb204nf_jyH2(%esp),%xmm4
        subps  nb204nf_jzH2(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb204nf_rsqH1H1(%esp)
        movaps %xmm3,nb204nf_rsqH1H2(%esp)

        movaps nb204nf_ixH1(%esp),%xmm0
        movaps nb204nf_iyH1(%esp),%xmm1
        movaps nb204nf_izH1(%esp),%xmm2
        movaps nb204nf_ixH2(%esp),%xmm3
        movaps nb204nf_iyH2(%esp),%xmm4
        movaps nb204nf_izH2(%esp),%xmm5
        subps  nb204nf_jxM(%esp),%xmm0
        subps  nb204nf_jyM(%esp),%xmm1
        subps  nb204nf_jzM(%esp),%xmm2
        subps  nb204nf_jxH1(%esp),%xmm3
        subps  nb204nf_jyH1(%esp),%xmm4
        subps  nb204nf_jzH1(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb204nf_rsqH1M(%esp)
        movaps %xmm3,nb204nf_rsqH2H1(%esp)

        movaps nb204nf_ixH2(%esp),%xmm0
        movaps nb204nf_iyH2(%esp),%xmm1
        movaps nb204nf_izH2(%esp),%xmm2
        movaps nb204nf_ixH2(%esp),%xmm3
        movaps nb204nf_iyH2(%esp),%xmm4
        movaps nb204nf_izH2(%esp),%xmm5
        subps  nb204nf_jxH2(%esp),%xmm0
        subps  nb204nf_jyH2(%esp),%xmm1
        subps  nb204nf_jzH2(%esp),%xmm2
        subps  nb204nf_jxM(%esp),%xmm3
        subps  nb204nf_jyM(%esp),%xmm4
        subps  nb204nf_jzM(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb204nf_rsqH2H2(%esp)
        movaps %xmm3,nb204nf_rsqH2M(%esp)

        movaps nb204nf_ixM(%esp),%xmm0
        movaps nb204nf_iyM(%esp),%xmm1
        movaps nb204nf_izM(%esp),%xmm2
        movaps nb204nf_ixM(%esp),%xmm3
        movaps nb204nf_iyM(%esp),%xmm4
        movaps nb204nf_izM(%esp),%xmm5
        subps  nb204nf_jxH1(%esp),%xmm0
        subps  nb204nf_jyH1(%esp),%xmm1
        subps  nb204nf_jzH1(%esp),%xmm2
        subps  nb204nf_jxH2(%esp),%xmm3
        subps  nb204nf_jyH2(%esp),%xmm4
        subps  nb204nf_jzH2(%esp),%xmm5
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb204nf_rsqMH1(%esp)
        movaps %xmm4,nb204nf_rsqMH2(%esp)

        movaps nb204nf_ixM(%esp),%xmm0
        movaps nb204nf_iyM(%esp),%xmm1
        movaps nb204nf_izM(%esp),%xmm2
        subps  nb204nf_jxM(%esp),%xmm0
        subps  nb204nf_jyM(%esp),%xmm1
        subps  nb204nf_jzM(%esp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb204nf_rsqMM(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%esp),%xmm3   ## rinvH2H2 
        mulps   nb204nf_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb204nf_rinvMM(%esp)
        movaps  %xmm7,nb204nf_rinvMH2(%esp)

        rsqrtps nb204nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb204nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204nf_rsqH1H1(%esp),%xmm1
        mulps   nb204nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%esp),%xmm3
        mulps   nb204nf_half(%esp),%xmm7
        movaps  %xmm3,nb204nf_rinvH1H1(%esp)
        movaps  %xmm7,nb204nf_rinvH1H2(%esp)

        rsqrtps nb204nf_rsqH1M(%esp),%xmm1
        rsqrtps nb204nf_rsqH2H1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204nf_rsqH1M(%esp),%xmm1
        mulps   nb204nf_rsqH2H1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%esp),%xmm3
        mulps   nb204nf_half(%esp),%xmm7
        movaps  %xmm3,nb204nf_rinvH1M(%esp)
        movaps  %xmm7,nb204nf_rinvH2H1(%esp)

        rsqrtps nb204nf_rsqH2H2(%esp),%xmm1
        rsqrtps nb204nf_rsqH2M(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb204nf_rsqH2H2(%esp),%xmm1
        mulps   nb204nf_rsqH2M(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%esp),%xmm3
        mulps   nb204nf_half(%esp),%xmm7
        movaps  %xmm3,nb204nf_rinvH2H2(%esp)
        movaps  %xmm7,nb204nf_rinvH2M(%esp)

        rsqrtps nb204nf_rsqMH1(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb204nf_three(%esp),%xmm3
        mulps   nb204nf_rsqMH1(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb204nf_half(%esp),%xmm3
        movaps  %xmm3,nb204nf_rinvMH1(%esp)

        ## Coulomb interactions 
        ## add all H-H rsq in xmm2, H-M rsq xmm4
        ## H-H rinv in xmm0, H-M in xmm1
        movaps nb204nf_rinvH1H1(%esp),%xmm0
        movaps nb204nf_rinvH1M(%esp),%xmm1
        movaps nb204nf_rsqH1H1(%esp),%xmm2
        movaps nb204nf_rsqH1M(%esp),%xmm4
        addps  nb204nf_rinvH1H2(%esp),%xmm0
        addps  nb204nf_rinvH2M(%esp),%xmm1
        addps  nb204nf_rsqH1H2(%esp),%xmm2
        addps  nb204nf_rsqH2M(%esp),%xmm4
        addps  nb204nf_rinvH2H1(%esp),%xmm0
        addps  nb204nf_rinvMH1(%esp),%xmm1
        addps  nb204nf_rsqH2H1(%esp),%xmm2
        addps  nb204nf_rsqMH1(%esp),%xmm4
        addps  nb204nf_rinvH2H2(%esp),%xmm0
        addps  nb204nf_rinvMH2(%esp),%xmm1
        addps  nb204nf_rsqH2H2(%esp),%xmm2
        addps  nb204nf_rsqMH2(%esp),%xmm4
        movaps nb204nf_krf(%esp),%xmm5
        movaps nb204nf_crf(%esp),%xmm6

        ## calc 4*crf in xmm7
        movaps %xmm6,%xmm7
        addps  %xmm7,%xmm7
        addps  %xmm7,%xmm7
        mulps  %xmm5,%xmm2 ## H-H krsq
        mulps  %xmm5,%xmm4 ## H-M krsq
        addps  %xmm2,%xmm0 ## H-H rinv+krsq
        addps  %xmm4,%xmm1 ## H-M rinv+krsq
        subps  %xmm7,%xmm0 ## H-H rinv+krsq-crf
        subps  %xmm7,%xmm1 ## H-M rinv+krsq-crf
        mulps  nb204nf_qqHH(%esp),%xmm0
        mulps  nb204nf_qqMH(%esp),%xmm1
        addps  %xmm1,%xmm0
        addps  nb204nf_vctot(%esp),%xmm0
        ## M-M interaction
        movaps nb204nf_rinvMM(%esp),%xmm4
        mulps  nb204nf_rsqMM(%esp),%xmm5   ## krsq
        addps  %xmm4,%xmm5                 ## rinv+krsq
        subps nb204nf_crf(%esp),%xmm5   ## xmm5=rinv+ krsq-crf 
        mulps nb204nf_qqMM(%esp),%xmm5
        addps %xmm0,%xmm5
        movaps %xmm5,nb204nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb204nf_innerk(%esp)
        jl    _nb_kernel204nf_ia32_sse.nb204nf_single_check
        jmp   _nb_kernel204nf_ia32_sse.nb204nf_unroll_loop
_nb_kernel204nf_ia32_sse.nb204nf_single_check: 
        addl $4,nb204nf_innerk(%esp)
        jnz   _nb_kernel204nf_ia32_sse.nb204nf_single_loop
        jmp   _nb_kernel204nf_ia32_sse.nb204nf_updateouterdata
_nb_kernel204nf_ia32_sse.nb204nf_single_loop: 
        movl  nb204nf_innerjjnr(%esp),%edx      ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb204nf_innerjjnr(%esp)

        movl nb204nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5
        movss 36(%esi,%eax,4),%xmm3             ## jxM  -  -  -
        movss 40(%esi,%eax,4),%xmm4             ## jyM  -  -  -
        movss 44(%esi,%eax,4),%xmm5             ## jzM  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb204nf_ixM(%esp),%xmm0
        movaps  nb204nf_iyM(%esp),%xmm1
        movaps  nb204nf_izM(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxM   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyM   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzM   0   jzH1 jzH2

        ## store all j coordinates in jM
        movaps %xmm3,nb204nf_jxM(%esp)
        movaps %xmm4,nb204nf_jyM(%esp)
        movaps %xmm5,nb204nf_jzM(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        mulps   nb204nf_krf(%esp),%xmm6   ## xmm6=krsq 
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb204nf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb204nf_half(%esp),%xmm3   ## rinv iO - j water 

        addps   %xmm3,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb204nf_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 

        xorps   %xmm4,%xmm4
        ## fetch charges to xmm4 (temporary) 
        movss   nb204nf_qqMM(%esp),%xmm4
        movhps  nb204nf_qqMH(%esp),%xmm4
        mulps %xmm4,%xmm6       ## vcoul  
        addps   nb204nf_vctot(%esp),%xmm6
        movaps  %xmm6,nb204nf_vctot(%esp)

        ## done with i M Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb204nf_ixH1(%esp),%xmm0
        movaps  nb204nf_iyH1(%esp),%xmm1
        movaps  nb204nf_izH1(%esp),%xmm2
        movaps  nb204nf_ixH2(%esp),%xmm3
        movaps  nb204nf_iyH2(%esp),%xmm4
        movaps  nb204nf_izH2(%esp),%xmm5
        subps   nb204nf_jxM(%esp),%xmm0
        subps   nb204nf_jyM(%esp),%xmm1
        subps   nb204nf_jzM(%esp),%xmm2
        subps   nb204nf_jxM(%esp),%xmm3
        subps   nb204nf_jyM(%esp),%xmm4
        subps   nb204nf_jzM(%esp),%xmm5
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        mulps %xmm3,%xmm3
        mulps %xmm4,%xmm4
        mulps %xmm5,%xmm5
        addps %xmm1,%xmm0
        addps %xmm3,%xmm4
        addps %xmm2,%xmm0       ## have rsqH1 in xmm0 
        addps %xmm5,%xmm4       ## have rsqH2 in xmm4 

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb204nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb204nf_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb204nf_half(%esp),%xmm7   ## rinv H2 - j water  

        mulps nb204nf_krf(%esp),%xmm0   ## krsq 
        mulps nb204nf_krf(%esp),%xmm4   ## krsq  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb204nf_qqMH(%esp),%xmm6
        movhps  nb204nf_qqHH(%esp),%xmm6
        addps   %xmm3,%xmm0     ## krsq+ rinv 
        addps   %xmm7,%xmm4     ## krsq+ rinv 
        subps   nb204nf_crf(%esp),%xmm0
        subps   nb204nf_crf(%esp),%xmm4
        mulps   %xmm6,%xmm0     ## vcoul 
        mulps   %xmm6,%xmm4     ## vcoul 
        addps   %xmm0,%xmm4
        addps   nb204nf_vctot(%esp),%xmm4
        movaps  %xmm4,nb204nf_vctot(%esp)
        decl nb204nf_innerk(%esp)
        jz    _nb_kernel204nf_ia32_sse.nb204nf_updateouterdata
        jmp   _nb_kernel204nf_ia32_sse.nb204nf_single_loop
_nb_kernel204nf_ia32_sse.nb204nf_updateouterdata: 
        ## get n from stack
        movl nb204nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb204nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb204nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb204nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb204nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel204nf_ia32_sse.nb204nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb204nf_n(%esp)
        jmp _nb_kernel204nf_ia32_sse.nb204nf_outer
_nb_kernel204nf_ia32_sse.nb204nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb204nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel204nf_ia32_sse.nb204nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel204nf_ia32_sse.nb204nf_threadloop
_nb_kernel204nf_ia32_sse.nb204nf_end: 
        emms

        movl nb204nf_nouter(%esp),%eax
        movl nb204nf_ninner(%esp),%ebx
        movl nb204nf_outeriter(%ebp),%ecx
        movl nb204nf_inneriter(%ebp),%edx
        movl %eax,(%ecx)
        movl %ebx,(%edx)

        movl nb204nf_salign(%esp),%eax
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

