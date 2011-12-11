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




.globl nb_kernel202_ia32_sse
.globl _nb_kernel202_ia32_sse
nb_kernel202_ia32_sse:  
_nb_kernel202_ia32_sse: 
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
        ## bottom of stack is cache-aligned for sse use 
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
.set nb202_fjzH2b, 1060
.set nb202_fjzH2c, 1064
.set nb202_fjzH2d, 1068
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


        movl nb202_argkrf(%ebp),%esi
        movl nb202_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb202_krf(%esp)
        movaps %xmm6,nb202_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb202_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb202_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%edx,%ebx,4),%xmm5
        movl nb202_p_facel(%ebp),%esi
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
        movaps %xmm3,nb202_qqOO(%esp)
        movaps %xmm4,nb202_qqOH(%esp)
        movaps %xmm5,nb202_qqHH(%esp)

        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb202_half(%esp)
        movss nb202_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb202_half(%esp)
        movaps %xmm2,nb202_two(%esp)
        movaps %xmm3,nb202_three(%esp)

_nb_kernel202_ia32_sse.nb202_threadloop: 
        movl  nb202_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel202_ia32_sse.nb202_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel202_ia32_sse.nb202_spinlock

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
        jg  _nb_kernel202_ia32_sse.nb202_outerstart
        jmp _nb_kernel202_ia32_sse.nb202_end

_nb_kernel202_ia32_sse.nb202_outerstart: 
        ## ebx contains number of outer iterations
        addl nb202_nouter(%esp),%ebx
        movl %ebx,nb202_nouter(%esp)

_nb_kernel202_ia32_sse.nb202_outer: 
        movl  nb202_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx        ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb202_is3(%esp)      ## store is3 

        movl  nb202_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb202_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]        
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb202_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb202_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb202_ixO(%esp)
        movaps %xmm4,nb202_iyO(%esp)
        movaps %xmm5,nb202_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb202_ixH1(%esp)
        movaps %xmm1,nb202_iyH1(%esp)
        movaps %xmm2,nb202_izH1(%esp)
        movaps %xmm3,nb202_ixH2(%esp)
        movaps %xmm4,nb202_iyH2(%esp)
        movaps %xmm5,nb202_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb202_vctot(%esp)
        movaps %xmm4,nb202_fixO(%esp)
        movaps %xmm4,nb202_fiyO(%esp)
        movaps %xmm4,nb202_fizO(%esp)
        movaps %xmm4,nb202_fixH1(%esp)
        movaps %xmm4,nb202_fiyH1(%esp)
        movaps %xmm4,nb202_fizH1(%esp)
        movaps %xmm4,nb202_fixH2(%esp)
        movaps %xmm4,nb202_fiyH2(%esp)
        movaps %xmm4,nb202_fizH2(%esp)

        movl  nb202_jindex(%ebp),%eax
        movl  (%eax,%esi,4),%ecx                ## jindex[n] 
        movl  4(%eax,%esi,4),%edx               ## jindex[n+1] 
        subl  %ecx,%edx                         ## number of innerloop atoms 

        movl  nb202_pos(%ebp),%esi
        movl  nb202_faction(%ebp),%edi
        movl  nb202_jjnr(%ebp),%eax
        shll  $2,%ecx
        addl  %ecx,%eax
        movl  %eax,nb202_innerjjnr(%esp)       ## pointer to jjnr[nj0] 
        movl  %edx,%ecx
        subl  $4,%edx
        addl  nb202_ninner(%esp),%ecx
        movl  %ecx,nb202_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb202_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel202_ia32_sse.nb202_unroll_loop
        jmp   _nb_kernel202_ia32_sse.nb202_single_check
_nb_kernel202_ia32_sse.nb202_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb202_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb202_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb202_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps (%esi,%eax,4),%xmm2
        movlps 12(%esi,%eax,4),%xmm3
        movlps 24(%esi,%eax,4),%xmm4

        movlps (%esi,%ebx,4),%xmm5
        movlps 12(%esi,%ebx,4),%xmm6
        movlps 24(%esi,%ebx,4),%xmm7

        movhps (%esi,%ecx,4),%xmm2
        movhps 12(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ecx,4),%xmm4

        movhps (%esi,%edx,4),%xmm5
        movhps 12(%esi,%edx,4),%xmm6
        movhps 24(%esi,%edx,4),%xmm7

        ## current state:       
        ## xmm2= jxOa  jyOa  jxOc  jyOc 
        ## xmm3= jxH1a jyH1a jxH1c jyH1c 
        ## xmm4= jxH2a jyH2a jxH2c jyH2c 
        ## xmm5= jxOb  jyOb  jxOd  jyOd 
        ## xmm6= jxH1b jyH1b jxH1d jyH1d 
        ## xmm7= jxH2b jyH2b jxH2d jyH2d 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxOa  jxOb  jyOa  jyOb 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH1a jxH1b jyH1a jyH1b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxOc  jxOd  jyOc  jyOd 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH1c jxH1d jyH1c jyH1d 
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxH2a jxH2b jyH2a jyH2b                
        unpckhps %xmm7,%xmm5    ## xmm5= jxH2c jxH2d jyH2c jyH2d 
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxOa  jxOb  jxOc  jxOd 
        movaps %xmm0,nb202_jxO(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb202_jyO(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb202_jxH1(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb202_jyH1(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb202_jxH2(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb202_jyH2(%esp)

        movss  8(%esi,%eax,4),%xmm0
        movss  20(%esi,%eax,4),%xmm1
        movss  32(%esi,%eax,4),%xmm2

        movss  8(%esi,%ecx,4),%xmm3
        movss  20(%esi,%ecx,4),%xmm4
        movss  32(%esi,%ecx,4),%xmm5

        movhps 4(%esi,%ebx,4),%xmm0
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 28(%esi,%ebx,4),%xmm2

        movhps 4(%esi,%edx,4),%xmm3
        movhps 16(%esi,%edx,4),%xmm4
        movhps 28(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb202_jzO(%esp)
        movaps %xmm1,nb202_jzH1(%esp)
        movaps %xmm2,nb202_jzH2(%esp)

        movaps nb202_ixO(%esp),%xmm0
        movaps nb202_iyO(%esp),%xmm1
        movaps nb202_izO(%esp),%xmm2
        movaps nb202_ixO(%esp),%xmm3
        movaps nb202_iyO(%esp),%xmm4
        movaps nb202_izO(%esp),%xmm5
        subps  nb202_jxO(%esp),%xmm0
        subps  nb202_jyO(%esp),%xmm1
        subps  nb202_jzO(%esp),%xmm2
        subps  nb202_jxH1(%esp),%xmm3
        subps  nb202_jyH1(%esp),%xmm4
        subps  nb202_jzH1(%esp),%xmm5
        movaps %xmm0,nb202_dxOO(%esp)
        movaps %xmm1,nb202_dyOO(%esp)
        movaps %xmm2,nb202_dzOO(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb202_dxOH1(%esp)
        movaps %xmm4,nb202_dyOH1(%esp)
        movaps %xmm5,nb202_dzOH1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb202_rsqOO(%esp)
        movaps %xmm3,nb202_rsqOH1(%esp)

        movaps nb202_ixO(%esp),%xmm0
        movaps nb202_iyO(%esp),%xmm1
        movaps nb202_izO(%esp),%xmm2
        movaps nb202_ixH1(%esp),%xmm3
        movaps nb202_iyH1(%esp),%xmm4
        movaps nb202_izH1(%esp),%xmm5
        subps  nb202_jxH2(%esp),%xmm0
        subps  nb202_jyH2(%esp),%xmm1
        subps  nb202_jzH2(%esp),%xmm2
        subps  nb202_jxO(%esp),%xmm3
        subps  nb202_jyO(%esp),%xmm4
        subps  nb202_jzO(%esp),%xmm5
        movaps %xmm0,nb202_dxOH2(%esp)
        movaps %xmm1,nb202_dyOH2(%esp)
        movaps %xmm2,nb202_dzOH2(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb202_dxH1O(%esp)
        movaps %xmm4,nb202_dyH1O(%esp)
        movaps %xmm5,nb202_dzH1O(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb202_rsqOH2(%esp)
        movaps %xmm3,nb202_rsqH1O(%esp)

        movaps nb202_ixH1(%esp),%xmm0
        movaps nb202_iyH1(%esp),%xmm1
        movaps nb202_izH1(%esp),%xmm2
        movaps nb202_ixH1(%esp),%xmm3
        movaps nb202_iyH1(%esp),%xmm4
        movaps nb202_izH1(%esp),%xmm5
        subps  nb202_jxH1(%esp),%xmm0
        subps  nb202_jyH1(%esp),%xmm1
        subps  nb202_jzH1(%esp),%xmm2
        subps  nb202_jxH2(%esp),%xmm3
        subps  nb202_jyH2(%esp),%xmm4
        subps  nb202_jzH2(%esp),%xmm5
        movaps %xmm0,nb202_dxH1H1(%esp)
        movaps %xmm1,nb202_dyH1H1(%esp)
        movaps %xmm2,nb202_dzH1H1(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb202_dxH1H2(%esp)
        movaps %xmm4,nb202_dyH1H2(%esp)
        movaps %xmm5,nb202_dzH1H2(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm4,%xmm3
        addps  %xmm5,%xmm3
        movaps %xmm0,nb202_rsqH1H1(%esp)
        movaps %xmm3,nb202_rsqH1H2(%esp)

        movaps nb202_ixH2(%esp),%xmm0
        movaps nb202_iyH2(%esp),%xmm1
        movaps nb202_izH2(%esp),%xmm2
        movaps nb202_ixH2(%esp),%xmm3
        movaps nb202_iyH2(%esp),%xmm4
        movaps nb202_izH2(%esp),%xmm5
        subps  nb202_jxO(%esp),%xmm0
        subps  nb202_jyO(%esp),%xmm1
        subps  nb202_jzO(%esp),%xmm2
        subps  nb202_jxH1(%esp),%xmm3
        subps  nb202_jyH1(%esp),%xmm4
        subps  nb202_jzH1(%esp),%xmm5
        movaps %xmm0,nb202_dxH2O(%esp)
        movaps %xmm1,nb202_dyH2O(%esp)
        movaps %xmm2,nb202_dzH2O(%esp)
        mulps  %xmm0,%xmm0
        mulps  %xmm1,%xmm1
        mulps  %xmm2,%xmm2
        movaps %xmm3,nb202_dxH2H1(%esp)
        movaps %xmm4,nb202_dyH2H1(%esp)
        movaps %xmm5,nb202_dzH2H1(%esp)
        mulps  %xmm3,%xmm3
        mulps  %xmm4,%xmm4
        mulps  %xmm5,%xmm5
        addps  %xmm1,%xmm0
        addps  %xmm2,%xmm0
        addps  %xmm3,%xmm4
        addps  %xmm5,%xmm4
        movaps %xmm0,nb202_rsqH2O(%esp)
        movaps %xmm4,nb202_rsqH2H1(%esp)

        movaps nb202_ixH2(%esp),%xmm0
        movaps nb202_iyH2(%esp),%xmm1
        movaps nb202_izH2(%esp),%xmm2
        subps  nb202_jxH2(%esp),%xmm0
        subps  nb202_jyH2(%esp),%xmm1
        subps  nb202_jzH2(%esp),%xmm2
        movaps %xmm0,nb202_dxH2H2(%esp)
        movaps %xmm1,nb202_dyH2H2(%esp)
        movaps %xmm2,nb202_dzH2H2(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb202_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202_half(%esp),%xmm3   ## rinvH2H2 
        mulps   nb202_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb202_rinvH2H2(%esp)
        movaps  %xmm7,nb202_rinvH2H1(%esp)

        rsqrtps nb202_rsqOO(%esp),%xmm1
        rsqrtps nb202_rsqOH1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202_rsqOO(%esp),%xmm1
        mulps   nb202_rsqOH1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202_half(%esp),%xmm3
        mulps   nb202_half(%esp),%xmm7
        movaps  %xmm3,nb202_rinvOO(%esp)
        movaps  %xmm7,nb202_rinvOH1(%esp)

        rsqrtps nb202_rsqOH2(%esp),%xmm1
        rsqrtps nb202_rsqH1O(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202_rsqOH2(%esp),%xmm1
        mulps   nb202_rsqH1O(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202_half(%esp),%xmm3
        mulps   nb202_half(%esp),%xmm7
        movaps  %xmm3,nb202_rinvOH2(%esp)
        movaps  %xmm7,nb202_rinvH1O(%esp)

        rsqrtps nb202_rsqH1H1(%esp),%xmm1
        rsqrtps nb202_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202_rsqH1H1(%esp),%xmm1
        mulps   nb202_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202_half(%esp),%xmm3
        mulps   nb202_half(%esp),%xmm7
        movaps  %xmm3,nb202_rinvH1H1(%esp)
        movaps  %xmm7,nb202_rinvH1H2(%esp)

        rsqrtps nb202_rsqH2O(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb202_three(%esp),%xmm3
        mulps   nb202_rsqH2O(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb202_half(%esp),%xmm3
        movaps  %xmm3,nb202_rinvH2O(%esp)

        ## start with OO interaction 
        movaps nb202_rinvOO(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        mulps  %xmm0,%xmm0      ## xmm0=rinvsq 

        mulps  nb202_rsqOO(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm6
        addps  %xmm7,%xmm6      ## xmm6=rinv+ krsq 
        subps  nb202_crf(%esp),%xmm6
        mulps  nb202_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        mulps nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqOO(%esp),%xmm7   ## xmm7 = coul part of fscal 

        addps  nb202_vctot(%esp),%xmm6   ## local vctot summation variable 
        mulps  %xmm7,%xmm0

        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb202_dxOO(%esp),%xmm0
        mulps nb202_dyOO(%esp),%xmm1
        mulps nb202_dzOO(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixO(%esp),%xmm0
        addps nb202_fiyO(%esp),%xmm1
        addps nb202_fizO(%esp),%xmm2
        movaps %xmm3,nb202_fjxO(%esp)
        movaps %xmm4,nb202_fjyO(%esp)
        movaps %xmm5,nb202_fjzO(%esp)
        movaps %xmm0,nb202_fixO(%esp)
        movaps %xmm1,nb202_fiyO(%esp)
        movaps %xmm2,nb202_fizO(%esp)

        ## O-H1 interaction 
        movaps nb202_rinvOH1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb202_rsqOH1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb202_crf(%esp),%xmm4
        mulps  %xmm0,%xmm0
        mulps  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH1  
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb202_dxOH1(%esp),%xmm0
        mulps nb202_dyOH1(%esp),%xmm1
        mulps nb202_dzOH1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixO(%esp),%xmm0
        addps nb202_fiyO(%esp),%xmm1
        addps nb202_fizO(%esp),%xmm2
        movaps %xmm3,nb202_fjxH1(%esp)
        movaps %xmm4,nb202_fjyH1(%esp)
        movaps %xmm5,nb202_fjzH1(%esp)
        movaps %xmm0,nb202_fixO(%esp)
        movaps %xmm1,nb202_fiyO(%esp)
        movaps %xmm2,nb202_fizO(%esp)

        ## O-H2 interaction  
        movaps nb202_rinvOH2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb202_rsqOH2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb202_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        xorps %xmm3,%xmm3
        movaps %xmm3,%xmm4
        movaps %xmm3,%xmm5
        mulps nb202_dxOH2(%esp),%xmm0
        mulps nb202_dyOH2(%esp),%xmm1
        mulps nb202_dzOH2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixO(%esp),%xmm0
        addps nb202_fiyO(%esp),%xmm1
        addps nb202_fizO(%esp),%xmm2
        movaps %xmm3,nb202_fjxH2(%esp)
        movaps %xmm4,nb202_fjyH2(%esp)
        movaps %xmm5,nb202_fjzH2(%esp)
        movaps %xmm0,nb202_fixO(%esp)
        movaps %xmm1,nb202_fiyO(%esp)
        movaps %xmm2,nb202_fizO(%esp)

        ## H1-O interaction 
        movaps nb202_rinvH1O(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb202_rsqH1O(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb202_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb202_fjxO(%esp),%xmm3
        movaps nb202_fjyO(%esp),%xmm4
        movaps nb202_fjzO(%esp),%xmm5
        mulps nb202_dxH1O(%esp),%xmm0
        mulps nb202_dyH1O(%esp),%xmm1
        mulps nb202_dzH1O(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixH1(%esp),%xmm0
        addps nb202_fiyH1(%esp),%xmm1
        addps nb202_fizH1(%esp),%xmm2
        movaps %xmm3,nb202_fjxO(%esp)
        movaps %xmm4,nb202_fjyO(%esp)
        movaps %xmm5,nb202_fjzO(%esp)
        movaps %xmm0,nb202_fixH1(%esp)
        movaps %xmm1,nb202_fiyH1(%esp)
        movaps %xmm2,nb202_fizH1(%esp)

        ## H1-H1 interaction 
        movaps nb202_rinvH1H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb202_rsqH1H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb202_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb202_fjxH1(%esp),%xmm3
        movaps nb202_fjyH1(%esp),%xmm4
        movaps nb202_fjzH1(%esp),%xmm5
        mulps nb202_dxH1H1(%esp),%xmm0
        mulps nb202_dyH1H1(%esp),%xmm1
        mulps nb202_dzH1H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixH1(%esp),%xmm0
        addps nb202_fiyH1(%esp),%xmm1
        addps nb202_fizH1(%esp),%xmm2
        movaps %xmm3,nb202_fjxH1(%esp)
        movaps %xmm4,nb202_fjyH1(%esp)
        movaps %xmm5,nb202_fjzH1(%esp)
        movaps %xmm0,nb202_fixH1(%esp)
        movaps %xmm1,nb202_fiyH1(%esp)
        movaps %xmm2,nb202_fizH1(%esp)

        ## H1-H2 interaction 
        movaps nb202_rinvH1H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb202_rsqH1H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb202_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb202_fjxH2(%esp),%xmm3
        movaps nb202_fjyH2(%esp),%xmm4
        movaps nb202_fjzH2(%esp),%xmm5
        mulps nb202_dxH1H2(%esp),%xmm0
        mulps nb202_dyH1H2(%esp),%xmm1
        mulps nb202_dzH1H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixH1(%esp),%xmm0
        addps nb202_fiyH1(%esp),%xmm1
        addps nb202_fizH1(%esp),%xmm2
        movaps %xmm3,nb202_fjxH2(%esp)
        movaps %xmm4,nb202_fjyH2(%esp)
        movaps %xmm5,nb202_fjzH2(%esp)
        movaps %xmm0,nb202_fixH1(%esp)
        movaps %xmm1,nb202_fiyH1(%esp)
        movaps %xmm2,nb202_fizH1(%esp)

        ## H2-O interaction 
        movaps nb202_rinvH2O(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb202_rsqH2O(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb202_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb202_qqOH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqOH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb202_fjxO(%esp),%xmm3
        movaps nb202_fjyO(%esp),%xmm4
        movaps nb202_fjzO(%esp),%xmm5
        mulps nb202_dxH2O(%esp),%xmm0
        mulps nb202_dyH2O(%esp),%xmm1
        mulps nb202_dzH2O(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixH2(%esp),%xmm0
        addps nb202_fiyH2(%esp),%xmm1
        addps nb202_fizH2(%esp),%xmm2
        movaps %xmm3,nb202_fjxO(%esp)
        movaps %xmm4,nb202_fjyO(%esp)
        movaps %xmm5,nb202_fjzO(%esp)
        movaps %xmm0,nb202_fixH2(%esp)
        movaps %xmm1,nb202_fiyH2(%esp)
        movaps %xmm2,nb202_fizH2(%esp)

        ## H2-H1 interaction 
        movaps nb202_rinvH2H1(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb202_rsqH2H1(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb202_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps nb202_fjxH1(%esp),%xmm3
        movaps nb202_fjyH1(%esp),%xmm4
        movaps nb202_fjzH1(%esp),%xmm5
        mulps nb202_dxH2H1(%esp),%xmm0
        mulps nb202_dyH2H1(%esp),%xmm1
        mulps nb202_dzH2H1(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixH2(%esp),%xmm0
        addps nb202_fiyH2(%esp),%xmm1
        addps nb202_fizH2(%esp),%xmm2
        movaps %xmm3,nb202_fjxH1(%esp)
        movaps %xmm4,nb202_fjyH1(%esp)
        movaps %xmm5,nb202_fjzH1(%esp)
        movaps %xmm0,nb202_fixH2(%esp)
        movaps %xmm1,nb202_fiyH2(%esp)
        movaps %xmm2,nb202_fizH2(%esp)

        ## H2-H2 interaction 
        movaps nb202_rinvH2H2(%esp),%xmm0
        movaps %xmm0,%xmm7      ## xmm7=rinv 
        movaps nb202_krf(%esp),%xmm5
        movaps %xmm0,%xmm1
        mulps  nb202_rsqH2H2(%esp),%xmm5   ## xmm5=krsq 
        movaps %xmm5,%xmm4
        addps  %xmm7,%xmm4      ## xmm4=r inv+ krsq 
        subps  nb202_crf(%esp),%xmm4
        mulps %xmm0,%xmm0
        mulps  nb202_qqHH(%esp),%xmm4   ## xmm4=voul=qq*(rinv+ krsq-crf) 
        mulps  nb202_two(%esp),%xmm5
        subps  %xmm5,%xmm7      ## xmm7=rinv-2*krsq 
        mulps  nb202_qqHH(%esp),%xmm7   ## xmm7 = coul part of fscal 
        addps  %xmm4,%xmm6      ## add to local vctot 
        mulps %xmm7,%xmm0       ## fsOH2 
        movaps %xmm0,%xmm1
        movaps %xmm0,%xmm2

        movaps %xmm0,%xmm1
        movaps %xmm6,nb202_vctot(%esp)
        movaps %xmm0,%xmm2

        movaps nb202_fjxH2(%esp),%xmm3
        movaps nb202_fjyH2(%esp),%xmm4
        movaps nb202_fjzH2(%esp),%xmm5
        mulps nb202_dxH2H2(%esp),%xmm0
        mulps nb202_dyH2H2(%esp),%xmm1
        mulps nb202_dzH2H2(%esp),%xmm2
        subps %xmm0,%xmm3
        subps %xmm1,%xmm4
        subps %xmm2,%xmm5
        addps nb202_fixH2(%esp),%xmm0
        addps nb202_fiyH2(%esp),%xmm1
        addps nb202_fizH2(%esp),%xmm2
        movaps %xmm3,nb202_fjxH2(%esp)
        movaps %xmm4,nb202_fjyH2(%esp)
        movaps %xmm5,nb202_fjzH2(%esp)
        movaps %xmm0,nb202_fixH2(%esp)
        movaps %xmm1,nb202_fiyH2(%esp)
        movaps %xmm2,nb202_fizH2(%esp)

        movl nb202_faction(%ebp),%edi

        ## Did all interactions - now update j forces 
        ## At this stage forces are still on the stack, in positions:
        ## fjxO, fjyO, fjzO, ... , fjzH2.
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


        ## 4 j waters with three atoms each - first do Oxygen X & Y forces for 4 j particles 
        movaps nb202_fjxO(%esp),%xmm0   ## xmm0= fjxOa  fjxOb  fjxOc  fjxOd 
        movaps nb202_fjyO(%esp),%xmm2   ## xmm1= fjyOa  fjyOb  fjyOc  fjyOd
        movlps (%edi,%eax,4),%xmm3
        movlps (%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjxOa  fjyOa  fjxOb  fjyOb
        unpckhps %xmm2,%xmm1       ## xmm1= fjxOc  fjyOc  fjxOd  fjyOd
        movhps (%edi,%ebx,4),%xmm3
        movhps (%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,(%edi,%eax,4)
        movlps %xmm4,(%edi,%ecx,4)
        movhps %xmm3,(%edi,%ebx,4)
        movhps %xmm4,(%edi,%edx,4)

        ## Oxygen Z & first hydrogen X forces for 4 j particles 
        movaps nb202_fjzO(%esp),%xmm0    ## xmm0= fjzOa   fjzOb   fjzOc   fjzOd 
        movaps nb202_fjxH1(%esp),%xmm2   ## xmm1= fjxH1a  fjxH1b  fjxH1c  fjxH1d
        movlps 8(%edi,%eax,4),%xmm3
        movlps 8(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0       ## xmm0= fjzOa  fjxH1a  fjzOb  fjxH1b
        unpckhps %xmm2,%xmm1       ## xmm1= fjzOc  fjxH1c  fjzOd  fjxH1d
        movhps 8(%edi,%ebx,4),%xmm3
        movhps 8(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,8(%edi,%eax,4)
        movlps %xmm4,8(%edi,%ecx,4)
        movhps %xmm3,8(%edi,%ebx,4)
        movhps %xmm4,8(%edi,%edx,4)


        ## First hydrogen Y & Z forces for 4 j particles 
        movaps nb202_fjyH1(%esp),%xmm0    ## xmm0= fjyH1a  fjyH1b  fjyH1c  fjyH1d 
        movaps nb202_fjzH1(%esp),%xmm2   ## xmm1= fjzH1a  fjzH1b  fjzH1c  fjzH1d
        movlps 16(%edi,%eax,4),%xmm3
        movlps 16(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjyH1a  fjzH1a  fjyH1b  fjzH1b
        unpckhps %xmm2,%xmm1            ## xmm1= fjyH1c  fjzH1c  fjyH1d  fjzH1d
        movhps 16(%edi,%ebx,4),%xmm3
        movhps 16(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,16(%edi,%eax,4)
        movlps %xmm4,16(%edi,%ecx,4)
        movhps %xmm3,16(%edi,%ebx,4)
        movhps %xmm4,16(%edi,%edx,4)


        ## Second hydrogen X & Y forces for 4 j particles 
        movaps nb202_fjxH2(%esp),%xmm0    ## xmm0= fjxH2a  fjxH2b  fjxH2c  fjxH2d 
        movaps nb202_fjyH2(%esp),%xmm2   ## xmm1= fjyH2a  fjyH2b  fjyH2c  fjyH2d
        movlps 24(%edi,%eax,4),%xmm3
        movlps 24(%edi,%ecx,4),%xmm4
        movaps %xmm0,%xmm1
        unpcklps %xmm2,%xmm0            ## xmm0= fjxH2a  fjyH2a  fjxH2b  fjyH2b
        unpckhps %xmm2,%xmm1            ## xmm1= fjxH2c  fjyH2c  fjxH2d  fjyH2d
        movhps 24(%edi,%ebx,4),%xmm3
        movhps 24(%edi,%edx,4),%xmm4
        addps  %xmm0,%xmm3
        addps  %xmm1,%xmm4
        movlps %xmm3,24(%edi,%eax,4)
        movlps %xmm4,24(%edi,%ecx,4)
        movhps %xmm3,24(%edi,%ebx,4)
        movhps %xmm4,24(%edi,%edx,4)


        ## Second hydrogen Z forces for 4 j particles 
        ## Just load the four Z coords into one reg. each
        movss 32(%edi,%eax,4),%xmm4
        movss 32(%edi,%ebx,4),%xmm5
        movss 32(%edi,%ecx,4),%xmm6
        movss 32(%edi,%edx,4),%xmm7
        ## add what we have on the stack
        addss nb202_fjzH2(%esp),%xmm4
        addss nb202_fjzH2b(%esp),%xmm5
        addss nb202_fjzH2c(%esp),%xmm6
        addss nb202_fjzH2d(%esp),%xmm7
        ## store back
        movss %xmm4,32(%edi,%eax,4)
        movss %xmm5,32(%edi,%ebx,4)
        movss %xmm6,32(%edi,%ecx,4)
        movss %xmm7,32(%edi,%edx,4)

        ## should we do one more iteration? 
        subl $4,nb202_innerk(%esp)
        jl    _nb_kernel202_ia32_sse.nb202_single_check
        jmp   _nb_kernel202_ia32_sse.nb202_unroll_loop
_nb_kernel202_ia32_sse.nb202_single_check: 
        addl $4,nb202_innerk(%esp)
        jnz   _nb_kernel202_ia32_sse.nb202_single_loop
        jmp   _nb_kernel202_ia32_sse.nb202_updateouterdata
_nb_kernel202_ia32_sse.nb202_single_loop: 
        movl  nb202_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb202_innerjjnr(%esp)

        movl nb202_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%esi,%eax,4),%xmm3               ## jxO  -  -  -
        movss 4(%esi,%eax,4),%xmm4              ## jyO  -  -  -
        movss 8(%esi,%eax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb202_ixO(%esp),%xmm0
        movaps  nb202_iyO(%esp),%xmm1
        movaps  nb202_izO(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb202_jxO(%esp)
        movaps %xmm4,nb202_jyO(%esp)
        movaps %xmm5,nb202_jzO(%esp)
        subps  %xmm3,%xmm0
        subps  %xmm4,%xmm1
        subps  %xmm5,%xmm2
        movaps %xmm0,nb202_dxOO(%esp)
        movaps %xmm1,nb202_dyOO(%esp)
        movaps %xmm2,nb202_dzOO(%esp)
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0       ## have rsq in xmm0 

        movaps %xmm0,%xmm6

        ## do invsqrt 
        rsqrtps %xmm0,%xmm1
        mulps   nb202_krf(%esp),%xmm6   ## xmm6=krsq 
        movaps  %xmm1,%xmm2
        movaps  %xmm6,%xmm7    ## xmm7=krsq 
        mulps   %xmm1,%xmm1
        movaps  nb202_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb202_half(%esp),%xmm3   ## rinv iO - j water 



        addps   %xmm3,%xmm6     ## xmm6=rinv+ krsq 
        mulps   nb202_two(%esp),%xmm7
        subps   nb202_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        subps   %xmm7,%xmm3     ## xmm3=rinv-2*krsq 
        xorps   %xmm4,%xmm4
        mulps   %xmm0,%xmm0     ## xmm0=rinvsq 
        ## fetch charges to xmm4 (temporary) 
        movss   nb202_qqOO(%esp),%xmm4
        movhps  nb202_qqOH(%esp),%xmm4

        mulps %xmm4,%xmm6       ## vcoul  
        mulps %xmm4,%xmm3       ## coul part of fs  


        addps   nb202_vctot(%esp),%xmm6
        mulps   %xmm3,%xmm0     ## total fscal 
        movaps  %xmm6,nb202_vctot(%esp)

        movaps  %xmm0,%xmm1
        movaps  %xmm0,%xmm2
        mulps   nb202_dxOO(%esp),%xmm0
        mulps   nb202_dyOO(%esp),%xmm1
        mulps   nb202_dzOO(%esp),%xmm2

        ## initial update for j forces 
        xorps   %xmm3,%xmm3
        xorps   %xmm4,%xmm4
        xorps   %xmm5,%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb202_fjxO(%esp)
        movaps  %xmm4,nb202_fjyO(%esp)
        movaps  %xmm5,nb202_fjzO(%esp)
        addps   nb202_fixO(%esp),%xmm0
        addps   nb202_fiyO(%esp),%xmm1
        addps   nb202_fizO(%esp),%xmm2
        movaps  %xmm0,nb202_fixO(%esp)
        movaps  %xmm1,nb202_fiyO(%esp)
        movaps  %xmm2,nb202_fizO(%esp)


        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb202_ixH1(%esp),%xmm0
        movaps  nb202_iyH1(%esp),%xmm1
        movaps  nb202_izH1(%esp),%xmm2
        movaps  nb202_ixH2(%esp),%xmm3
        movaps  nb202_iyH2(%esp),%xmm4
        movaps  nb202_izH2(%esp),%xmm5
        subps   nb202_jxO(%esp),%xmm0
        subps   nb202_jyO(%esp),%xmm1
        subps   nb202_jzO(%esp),%xmm2
        subps   nb202_jxO(%esp),%xmm3
        subps   nb202_jyO(%esp),%xmm4
        subps   nb202_jzO(%esp),%xmm5
        movaps %xmm0,nb202_dxH1O(%esp)
        movaps %xmm1,nb202_dyH1O(%esp)
        movaps %xmm2,nb202_dzH1O(%esp)
        movaps %xmm3,nb202_dxH2O(%esp)
        movaps %xmm4,nb202_dyH2O(%esp)
        movaps %xmm5,nb202_dzH2O(%esp)
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
        movaps  nb202_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb202_half(%esp),%xmm7   ## rinv H2 - j water  

        mulps nb202_krf(%esp),%xmm0   ## krsq 
        mulps nb202_krf(%esp),%xmm4   ## krsq  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb202_qqOH(%esp),%xmm6
        movhps  nb202_qqHH(%esp),%xmm6
        movaps  %xmm0,%xmm1
        movaps  %xmm4,%xmm5
        addps   %xmm3,%xmm0     ## krsq+ rinv 
        addps   %xmm7,%xmm4     ## krsq+ rinv 
        subps   nb202_crf(%esp),%xmm0
        subps   nb202_crf(%esp),%xmm4
        mulps   nb202_two(%esp),%xmm1
        mulps   nb202_two(%esp),%xmm5
        mulps   %xmm6,%xmm0     ## vcoul 
        mulps   %xmm6,%xmm4     ## vcoul 
        addps   %xmm0,%xmm4
        addps   nb202_vctot(%esp),%xmm4
        movaps  %xmm4,nb202_vctot(%esp)
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
        mulps   nb202_dxH1O(%esp),%xmm0
        mulps   nb202_dyH1O(%esp),%xmm1
        mulps   nb202_dzH1O(%esp),%xmm2
        ## update forces H1 - j water 
        movaps  nb202_fjxO(%esp),%xmm3
        movaps  nb202_fjyO(%esp),%xmm4
        movaps  nb202_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movaps  %xmm3,nb202_fjxO(%esp)
        movaps  %xmm4,nb202_fjyO(%esp)
        movaps  %xmm5,nb202_fjzO(%esp)
        addps   nb202_fixH1(%esp),%xmm0
        addps   nb202_fiyH1(%esp),%xmm1
        addps   nb202_fizH1(%esp),%xmm2
        movaps  %xmm0,nb202_fixH1(%esp)
        movaps  %xmm1,nb202_fiyH1(%esp)
        movaps  %xmm2,nb202_fizH1(%esp)
        ## do forces H2 - j water 
        movaps %xmm7,%xmm0
        movaps %xmm7,%xmm1
        movaps %xmm7,%xmm2
        mulps   nb202_dxH2O(%esp),%xmm0
        mulps   nb202_dyH2O(%esp),%xmm1
        mulps   nb202_dzH2O(%esp),%xmm2
        movaps  nb202_fjxO(%esp),%xmm3
        movaps  nb202_fjyO(%esp),%xmm4
        movaps  nb202_fjzO(%esp),%xmm5
        subps   %xmm0,%xmm3
        subps   %xmm1,%xmm4
        subps   %xmm2,%xmm5
        movl    nb202_faction(%ebp),%esi
        movaps  %xmm3,nb202_fjxO(%esp)
        movaps  %xmm4,nb202_fjyO(%esp)
        movaps  %xmm5,nb202_fjzO(%esp)
        addps   nb202_fixH2(%esp),%xmm0
        addps   nb202_fiyH2(%esp),%xmm1
        addps   nb202_fizH2(%esp),%xmm2
        movaps  %xmm0,nb202_fixH2(%esp)
        movaps  %xmm1,nb202_fiyH2(%esp)
        movaps  %xmm2,nb202_fizH2(%esp)

        ## update j water forces from local variables 
        movlps  (%esi,%eax,4),%xmm0
        movlps  12(%esi,%eax,4),%xmm1
        movhps  24(%esi,%eax,4),%xmm1
        movaps  nb202_fjxO(%esp),%xmm3
        movaps  nb202_fjyO(%esp),%xmm4
        movaps  nb202_fjzO(%esp),%xmm5
        movaps  %xmm5,%xmm6
        movaps  %xmm5,%xmm7
        shufps $2,%xmm6,%xmm6 ## constant 00000010
        shufps $3,%xmm7,%xmm7 ## constant 00000011
        addss   8(%esi,%eax,4),%xmm5
        addss   20(%esi,%eax,4),%xmm6
        addss   32(%esi,%eax,4),%xmm7
        movss   %xmm5,8(%esi,%eax,4)
        movss   %xmm6,20(%esi,%eax,4)
        movss   %xmm7,32(%esi,%eax,4)
        movaps   %xmm3,%xmm5
        unpcklps %xmm4,%xmm3
        unpckhps %xmm4,%xmm5
        addps    %xmm3,%xmm0
        addps    %xmm5,%xmm1
        movlps  %xmm0,(%esi,%eax,4)
        movlps  %xmm1,12(%esi,%eax,4)
        movhps  %xmm1,24(%esi,%eax,4)

        decl nb202_innerk(%esp)
        jz    _nb_kernel202_ia32_sse.nb202_updateouterdata
        jmp   _nb_kernel202_ia32_sse.nb202_single_loop
_nb_kernel202_ia32_sse.nb202_updateouterdata: 
        movl  nb202_ii3(%esp),%ecx
        movl  nb202_faction(%ebp),%edi
        movl  nb202_fshift(%ebp),%esi
        movl  nb202_is3(%esp),%edx

        ## accumulate  Oi forces in xmm0, xmm1, xmm2 
        movaps nb202_fixO(%esp),%xmm0
        movaps nb202_fiyO(%esp),%xmm1
        movaps nb202_fizO(%esp),%xmm2

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
        movss  (%edi,%ecx,4),%xmm3
        movss  4(%edi,%ecx,4),%xmm4
        movss  8(%edi,%ecx,4),%xmm5
        addss  %xmm0,%xmm3
        addss  %xmm1,%xmm4
        addss  %xmm2,%xmm5
        movss  %xmm3,(%edi,%ecx,4)
        movss  %xmm4,4(%edi,%ecx,4)
        movss  %xmm5,8(%edi,%ecx,4)

        ## accumulate force in xmm6/xmm7 for fshift 
        movaps %xmm0,%xmm6
        movss %xmm2,%xmm7
        movlhps %xmm1,%xmm6
        shufps $8,%xmm6,%xmm6 ## constant 00001000      

        ## accumulate H1i forces in xmm0, xmm1, xmm2 
        movaps nb202_fixH1(%esp),%xmm0
        movaps nb202_fiyH1(%esp),%xmm1
        movaps nb202_fizH1(%esp),%xmm2

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
        addss %xmm2,%xmm7
        movlhps %xmm1,%xmm0
        shufps $8,%xmm0,%xmm0 ## constant 00001000      
        addps   %xmm0,%xmm6

        ## accumulate H2i forces in xmm0, xmm1, xmm2 
        movaps nb202_fixH2(%esp),%xmm0
        movaps nb202_fiyH2(%esp),%xmm1
        movaps nb202_fizH2(%esp),%xmm2

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

        ## increment fshift force  
        movlps  (%esi,%edx,4),%xmm3
        movss  8(%esi,%edx,4),%xmm4
        addps  %xmm6,%xmm3
        addss  %xmm7,%xmm4
        movlps  %xmm3,(%esi,%edx,4)
        movss  %xmm4,8(%esi,%edx,4)

        ## get n from stack
        movl nb202_n(%esp),%esi
        ## get group index for i particle 
        movl  nb202_gid(%ebp),%edx              ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb202_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb202_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb202_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel202_ia32_sse.nb202_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb202_n(%esp)
        jmp _nb_kernel202_ia32_sse.nb202_outer
_nb_kernel202_ia32_sse.nb202_outerend: 
        ## check if more outer neighborlists remain
        movl  nb202_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel202_ia32_sse.nb202_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel202_ia32_sse.nb202_threadloop
_nb_kernel202_ia32_sse.nb202_end: 
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








.globl nb_kernel202nf_ia32_sse
.globl _nb_kernel202nf_ia32_sse
nb_kernel202nf_ia32_sse:        
_nb_kernel202nf_ia32_sse:       
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


        movl nb202nf_argkrf(%ebp),%esi
        movl nb202nf_argcrf(%ebp),%edi
        movss (%esi),%xmm5
        movss (%edi),%xmm6
        shufps $0,%xmm5,%xmm5
        shufps $0,%xmm6,%xmm6
        movaps %xmm5,nb202nf_krf(%esp)
        movaps %xmm6,nb202nf_crf(%esp)

        ## assume we have at least one i particle - start directly 
        movl  nb202nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx),%ebx           ## ebx =ii 

        movl  nb202nf_charge(%ebp),%edx
        movss (%edx,%ebx,4),%xmm3
        movss %xmm3,%xmm4
        movss 4(%edx,%ebx,4),%xmm5
        movl nb202nf_p_facel(%ebp),%esi
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
        movaps %xmm3,nb202nf_qqOO(%esp)
        movaps %xmm4,nb202nf_qqOH(%esp)
        movaps %xmm5,nb202nf_qqHH(%esp)


        ## create constant floating-point factors on stack
        movl $0x3f000000,%eax   ## constant 0.5 in IEEE (hex)
        movl %eax,nb202nf_half(%esp)
        movss nb202nf_half(%esp),%xmm1
        shufps $0,%xmm1,%xmm1  ## splat to all elements
        movaps %xmm1,%xmm2
        addps  %xmm2,%xmm2      ## constant 1.0
        movaps %xmm2,%xmm3
        addps  %xmm2,%xmm2      ## constant 2.0
        addps  %xmm2,%xmm3      ## constant 3.0
        movaps %xmm1,nb202nf_half(%esp)
        movaps %xmm3,nb202nf_three(%esp)

_nb_kernel202nf_ia32_sse.nb202nf_threadloop: 
        movl  nb202nf_count(%ebp),%esi            ## pointer to sync counter
        movl  (%esi),%eax
_nb_kernel202nf_ia32_sse.nb202nf_spinlock: 
        movl  %eax,%ebx                         ## ebx=*count=nn0
        addl  $1,%ebx                          ## ebx=nn1=nn0+10
        lock 
        cmpxchgl %ebx,(%esi)                    ## write nn1 to *counter,
                                                ## if it hasnt changed.
                                                ## or reread *counter to eax.
        pause                                   ## -> better p4 performance
        jnz _nb_kernel202nf_ia32_sse.nb202nf_spinlock

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
        jg  _nb_kernel202nf_ia32_sse.nb202nf_outerstart
        jmp _nb_kernel202nf_ia32_sse.nb202nf_end

_nb_kernel202nf_ia32_sse.nb202nf_outerstart: 
        ## ebx contains number of outer iterations
        addl nb202nf_nouter(%esp),%ebx
        movl %ebx,nb202nf_nouter(%esp)

_nb_kernel202nf_ia32_sse.nb202nf_outer: 
        movl  nb202nf_shift(%ebp),%eax        ## eax = pointer into shift[] 
        movl  (%eax,%esi,4),%ebx                ## ebx=shift[n] 

        leal  (%ebx,%ebx,2),%ebx    ## ebx=3*is 
        movl  %ebx,nb202nf_is3(%esp)            ## store is3 

        movl  nb202nf_shiftvec(%ebp),%eax     ## eax = base of shiftvec[] 

        movss (%eax,%ebx,4),%xmm0
        movss 4(%eax,%ebx,4),%xmm1
        movss 8(%eax,%ebx,4),%xmm2

        movl  nb202nf_iinr(%ebp),%ecx         ## ecx = pointer into iinr[]      
        movl  (%ecx,%esi,4),%ebx            ## ebx =ii 

        leal  (%ebx,%ebx,2),%ebx        ## ebx = 3*ii=ii3 
        movl  nb202nf_pos(%ebp),%eax      ## eax = base of pos[]  
        movl  %ebx,nb202nf_ii3(%esp)

        movaps %xmm0,%xmm3
        movaps %xmm1,%xmm4
        movaps %xmm2,%xmm5
        addss (%eax,%ebx,4),%xmm3
        addss 4(%eax,%ebx,4),%xmm4
        addss 8(%eax,%ebx,4),%xmm5
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm3,nb202nf_ixO(%esp)
        movaps %xmm4,nb202nf_iyO(%esp)
        movaps %xmm5,nb202nf_izO(%esp)

        movss %xmm0,%xmm3
        movss %xmm1,%xmm4
        movss %xmm2,%xmm5
        addss 12(%eax,%ebx,4),%xmm0
        addss 16(%eax,%ebx,4),%xmm1
        addss 20(%eax,%ebx,4),%xmm2
        addss 24(%eax,%ebx,4),%xmm3
        addss 28(%eax,%ebx,4),%xmm4
        addss 32(%eax,%ebx,4),%xmm5

        shufps $0,%xmm0,%xmm0
        shufps $0,%xmm1,%xmm1
        shufps $0,%xmm2,%xmm2
        shufps $0,%xmm3,%xmm3
        shufps $0,%xmm4,%xmm4
        shufps $0,%xmm5,%xmm5
        movaps %xmm0,nb202nf_ixH1(%esp)
        movaps %xmm1,nb202nf_iyH1(%esp)
        movaps %xmm2,nb202nf_izH1(%esp)
        movaps %xmm3,nb202nf_ixH2(%esp)
        movaps %xmm4,nb202nf_iyH2(%esp)
        movaps %xmm5,nb202nf_izH2(%esp)

        ## clear vctot and i forces 
        xorps %xmm4,%xmm4
        movaps %xmm4,nb202nf_vctot(%esp)

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
        subl  $4,%edx
        addl  nb202nf_ninner(%esp),%ecx
        movl  %ecx,nb202nf_ninner(%esp)
        addl  $0,%edx
        movl  %edx,nb202nf_innerk(%esp)      ## number of innerloop atoms 
        jge   _nb_kernel202nf_ia32_sse.nb202nf_unroll_loop
        jmp   _nb_kernel202nf_ia32_sse.nb202nf_single_check
_nb_kernel202nf_ia32_sse.nb202nf_unroll_loop: 
        ## quad-unroll innerloop here 
        movl  nb202nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 

        movl  (%edx),%eax
        movl  4(%edx),%ebx
        movl  8(%edx),%ecx
        movl  12(%edx),%edx           ## eax-edx=jnr1-4 

        addl $16,nb202nf_innerjjnr(%esp)             ## advance pointer (unrolled 4) 

        movl nb202nf_pos(%ebp),%esi        ## base of pos[] 

        leal  (%eax,%eax,2),%eax     ## replace jnr with j3 
        leal  (%ebx,%ebx,2),%ebx
        leal  (%ecx,%ecx,2),%ecx     ## replace jnr with j3 
        leal  (%edx,%edx,2),%edx

        ## move j coordinates to local temp variables 
        movlps (%esi,%eax,4),%xmm2
        movlps 12(%esi,%eax,4),%xmm3
        movlps 24(%esi,%eax,4),%xmm4

        movlps (%esi,%ebx,4),%xmm5
        movlps 12(%esi,%ebx,4),%xmm6
        movlps 24(%esi,%ebx,4),%xmm7

        movhps (%esi,%ecx,4),%xmm2
        movhps 12(%esi,%ecx,4),%xmm3
        movhps 24(%esi,%ecx,4),%xmm4

        movhps (%esi,%edx,4),%xmm5
        movhps 12(%esi,%edx,4),%xmm6
        movhps 24(%esi,%edx,4),%xmm7

        ## current state:       
        ## xmm2= jxOa  jyOa  jxOc  jyOc 
        ## xmm3= jxH1a jyH1a jxH1c jyH1c 
        ## xmm4= jxH2a jyH2a jxH2c jyH2c 
        ## xmm5= jxOb  jyOb  jxOd  jyOd 
        ## xmm6= jxH1b jyH1b jxH1d jyH1d 
        ## xmm7= jxH2b jyH2b jxH2d jyH2d 

        movaps %xmm2,%xmm0
        movaps %xmm3,%xmm1
        unpcklps %xmm5,%xmm0    ## xmm0= jxOa  jxOb  jyOa  jyOb 
        unpcklps %xmm6,%xmm1    ## xmm1= jxH1a jxH1b jyH1a jyH1b 
        unpckhps %xmm5,%xmm2    ## xmm2= jxOc  jxOd  jyOc  jyOd 
        unpckhps %xmm6,%xmm3    ## xmm3= jxH1c jxH1d jyH1c jyH1d 
        movaps %xmm4,%xmm5
        movaps   %xmm0,%xmm6
        unpcklps %xmm7,%xmm4    ## xmm4= jxH2a jxH2b jyH2a jyH2b                
        unpckhps %xmm7,%xmm5    ## xmm5= jxH2c jxH2d jyH2c jyH2d 
        movaps   %xmm1,%xmm7
        movlhps  %xmm2,%xmm0    ## xmm0= jxOa  jxOb  jxOc  jxOd 
        movaps %xmm0,nb202nf_jxO(%esp)
        movhlps  %xmm6,%xmm2    ## xmm2= jyOa  jyOb  jyOc  jyOd 
        movaps %xmm2,nb202nf_jyO(%esp)
        movlhps  %xmm3,%xmm1
        movaps %xmm1,nb202nf_jxH1(%esp)
        movhlps  %xmm7,%xmm3
        movaps   %xmm4,%xmm6
        movaps %xmm3,nb202nf_jyH1(%esp)
        movlhps  %xmm5,%xmm4
        movaps %xmm4,nb202nf_jxH2(%esp)
        movhlps  %xmm6,%xmm5
        movaps %xmm5,nb202nf_jyH2(%esp)

        movss  8(%esi,%eax,4),%xmm0
        movss  20(%esi,%eax,4),%xmm1
        movss  32(%esi,%eax,4),%xmm2

        movss  8(%esi,%ecx,4),%xmm3
        movss  20(%esi,%ecx,4),%xmm4
        movss  32(%esi,%ecx,4),%xmm5

        movhps 4(%esi,%ebx,4),%xmm0
        movhps 16(%esi,%ebx,4),%xmm1
        movhps 28(%esi,%ebx,4),%xmm2

        movhps 4(%esi,%edx,4),%xmm3
        movhps 16(%esi,%edx,4),%xmm4
        movhps 28(%esi,%edx,4),%xmm5

        shufps $204,%xmm3,%xmm0 ## constant 11001100
        shufps $204,%xmm4,%xmm1 ## constant 11001100
        shufps $204,%xmm5,%xmm2 ## constant 11001100
        movaps %xmm0,nb202nf_jzO(%esp)
        movaps %xmm1,nb202nf_jzH1(%esp)
        movaps %xmm2,nb202nf_jzH2(%esp)

        movaps nb202nf_ixO(%esp),%xmm0
        movaps nb202nf_iyO(%esp),%xmm1
        movaps nb202nf_izO(%esp),%xmm2
        movaps nb202nf_ixO(%esp),%xmm3
        movaps nb202nf_iyO(%esp),%xmm4
        movaps nb202nf_izO(%esp),%xmm5
        subps  nb202nf_jxO(%esp),%xmm0
        subps  nb202nf_jyO(%esp),%xmm1
        subps  nb202nf_jzO(%esp),%xmm2
        subps  nb202nf_jxH1(%esp),%xmm3
        subps  nb202nf_jyH1(%esp),%xmm4
        subps  nb202nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb202nf_rsqOO(%esp)
        movaps %xmm3,nb202nf_rsqOH1(%esp)

        movaps nb202nf_ixO(%esp),%xmm0
        movaps nb202nf_iyO(%esp),%xmm1
        movaps nb202nf_izO(%esp),%xmm2
        movaps nb202nf_ixH1(%esp),%xmm3
        movaps nb202nf_iyH1(%esp),%xmm4
        movaps nb202nf_izH1(%esp),%xmm5
        subps  nb202nf_jxH2(%esp),%xmm0
        subps  nb202nf_jyH2(%esp),%xmm1
        subps  nb202nf_jzH2(%esp),%xmm2
        subps  nb202nf_jxO(%esp),%xmm3
        subps  nb202nf_jyO(%esp),%xmm4
        subps  nb202nf_jzO(%esp),%xmm5
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
        movaps %xmm0,nb202nf_rsqOH2(%esp)
        movaps %xmm3,nb202nf_rsqH1O(%esp)

        movaps nb202nf_ixH1(%esp),%xmm0
        movaps nb202nf_iyH1(%esp),%xmm1
        movaps nb202nf_izH1(%esp),%xmm2
        movaps nb202nf_ixH1(%esp),%xmm3
        movaps nb202nf_iyH1(%esp),%xmm4
        movaps nb202nf_izH1(%esp),%xmm5
        subps  nb202nf_jxH1(%esp),%xmm0
        subps  nb202nf_jyH1(%esp),%xmm1
        subps  nb202nf_jzH1(%esp),%xmm2
        subps  nb202nf_jxH2(%esp),%xmm3
        subps  nb202nf_jyH2(%esp),%xmm4
        subps  nb202nf_jzH2(%esp),%xmm5
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
        movaps %xmm0,nb202nf_rsqH1H1(%esp)
        movaps %xmm3,nb202nf_rsqH1H2(%esp)

        movaps nb202nf_ixH2(%esp),%xmm0
        movaps nb202nf_iyH2(%esp),%xmm1
        movaps nb202nf_izH2(%esp),%xmm2
        movaps nb202nf_ixH2(%esp),%xmm3
        movaps nb202nf_iyH2(%esp),%xmm4
        movaps nb202nf_izH2(%esp),%xmm5
        subps  nb202nf_jxO(%esp),%xmm0
        subps  nb202nf_jyO(%esp),%xmm1
        subps  nb202nf_jzO(%esp),%xmm2
        subps  nb202nf_jxH1(%esp),%xmm3
        subps  nb202nf_jyH1(%esp),%xmm4
        subps  nb202nf_jzH1(%esp),%xmm5
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
        movaps %xmm0,nb202nf_rsqH2O(%esp)
        movaps %xmm4,nb202nf_rsqH2H1(%esp)

        movaps nb202nf_ixH2(%esp),%xmm0
        movaps nb202nf_iyH2(%esp),%xmm1
        movaps nb202nf_izH2(%esp),%xmm2
        subps  nb202nf_jxH2(%esp),%xmm0
        subps  nb202nf_jyH2(%esp),%xmm1
        subps  nb202nf_jzH2(%esp),%xmm2
        mulps %xmm0,%xmm0
        mulps %xmm1,%xmm1
        mulps %xmm2,%xmm2
        addps %xmm1,%xmm0
        addps %xmm2,%xmm0
        movaps %xmm0,nb202nf_rsqH2H2(%esp)

        ## start doing invsqrt use rsq values in xmm0, xmm4 
        rsqrtps %xmm0,%xmm1
        rsqrtps %xmm4,%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%esp),%xmm3   ## rinvH2H2 
        mulps   nb202nf_half(%esp),%xmm7   ## rinvH2H1 
        movaps  %xmm3,nb202nf_rinvH2H2(%esp)
        movaps  %xmm7,nb202nf_rinvH2H1(%esp)

        rsqrtps nb202nf_rsqOO(%esp),%xmm1
        rsqrtps nb202nf_rsqOH1(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202nf_rsqOO(%esp),%xmm1
        mulps   nb202nf_rsqOH1(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%esp),%xmm3
        mulps   nb202nf_half(%esp),%xmm7
        movaps  %xmm3,nb202nf_rinvOO(%esp)
        movaps  %xmm7,nb202nf_rinvOH1(%esp)

        rsqrtps nb202nf_rsqOH2(%esp),%xmm1
        rsqrtps nb202nf_rsqH1O(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202nf_rsqOH2(%esp),%xmm1
        mulps   nb202nf_rsqH1O(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%esp),%xmm3
        mulps   nb202nf_half(%esp),%xmm7
        movaps  %xmm3,nb202nf_rinvOH2(%esp)
        movaps  %xmm7,nb202nf_rinvH1O(%esp)

        rsqrtps nb202nf_rsqH1H1(%esp),%xmm1
        rsqrtps nb202nf_rsqH1H2(%esp),%xmm5
        movaps  %xmm1,%xmm2
        movaps  %xmm5,%xmm6
        mulps   %xmm1,%xmm1
        mulps   %xmm5,%xmm5
        movaps  nb202nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   nb202nf_rsqH1H1(%esp),%xmm1
        mulps   nb202nf_rsqH1H2(%esp),%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%esp),%xmm3
        mulps   nb202nf_half(%esp),%xmm7
        movaps  %xmm3,nb202nf_rinvH1H1(%esp)
        movaps  %xmm7,nb202nf_rinvH1H2(%esp)

        rsqrtps nb202nf_rsqH2O(%esp),%xmm1
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb202nf_three(%esp),%xmm3
        mulps   nb202nf_rsqH2O(%esp),%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb202nf_half(%esp),%xmm3
        movaps  %xmm3,nb202nf_rinvH2O(%esp)

        ## start with OO interaction 
        movaps nb202nf_krf(%esp),%xmm6
        mulps  nb202nf_rsqOO(%esp),%xmm6   ## xmm5=krsq 
        addps  nb202nf_rinvOO(%esp),%xmm6       ## xmm6=rinv+ krsq 
        subps  nb202nf_crf(%esp),%xmm6
        mulps  nb202nf_qqOO(%esp),%xmm6   ## xmm6=voul=qq*(rinv+ krsq-crf) 
        addps  nb202nf_vctot(%esp),%xmm6   ## local vctot summation variable 

        ## O-H interactions 
        movaps nb202nf_krf(%esp),%xmm0
        movaps nb202nf_krf(%esp),%xmm1
        movaps nb202nf_krf(%esp),%xmm2
        movaps nb202nf_krf(%esp),%xmm3
        mulps  nb202nf_rsqOH1(%esp),%xmm0   ## krsq 
        mulps  nb202nf_rsqOH2(%esp),%xmm1   ## krsq 
        mulps  nb202nf_rsqH1O(%esp),%xmm2   ## krsq 
        mulps  nb202nf_rsqH2O(%esp),%xmm3   ## krsq 
        addps  nb202nf_rinvOH1(%esp),%xmm0      ## rinv+ krsq 
        addps  nb202nf_rinvOH2(%esp),%xmm1      ## rinv+ krsq 
        addps  nb202nf_rinvH1O(%esp),%xmm2      ## rinv+ krsq 
        addps  nb202nf_rinvH2O(%esp),%xmm3      ## rinv+ krsq 
        subps  nb202nf_crf(%esp),%xmm0
        subps  nb202nf_crf(%esp),%xmm1
        subps  nb202nf_crf(%esp),%xmm2
        subps  nb202nf_crf(%esp),%xmm3
        mulps  nb202nf_qqOH(%esp),%xmm0   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqOH(%esp),%xmm1   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqOH(%esp),%xmm2   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqOH(%esp),%xmm3   ## voul=qq*(rinv+ krsq-crf) 
        addps %xmm0,%xmm6
        addps %xmm2,%xmm1
        addps %xmm3,%xmm6
        addps %xmm1,%xmm6

        ## H-H interactions 
        movaps nb202nf_krf(%esp),%xmm0
        movaps nb202nf_krf(%esp),%xmm1
        movaps nb202nf_krf(%esp),%xmm2
        movaps nb202nf_krf(%esp),%xmm3
        mulps  nb202nf_rsqH1H1(%esp),%xmm0   ## krsq 
        mulps  nb202nf_rsqH1H2(%esp),%xmm1   ## krsq 
        mulps  nb202nf_rsqH2H1(%esp),%xmm2   ## krsq 
        mulps  nb202nf_rsqH2H2(%esp),%xmm3   ## krsq 
        addps  nb202nf_rinvH1H1(%esp),%xmm0     ## rinv+ krsq 
        addps  nb202nf_rinvH1H2(%esp),%xmm1     ## rinv+ krsq 
        addps  nb202nf_rinvH2H1(%esp),%xmm2     ## rinv+ krsq 
        addps  nb202nf_rinvH2H2(%esp),%xmm3     ## rinv+ krsq 
        subps  nb202nf_crf(%esp),%xmm0
        subps  nb202nf_crf(%esp),%xmm1
        subps  nb202nf_crf(%esp),%xmm2
        subps  nb202nf_crf(%esp),%xmm3
        mulps  nb202nf_qqHH(%esp),%xmm0   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqHH(%esp),%xmm1   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqHH(%esp),%xmm2   ## voul=qq*(rinv+ krsq-crf) 
        mulps  nb202nf_qqHH(%esp),%xmm3   ## voul=qq*(rinv+ krsq-crf) 
        addps %xmm0,%xmm6
        addps %xmm2,%xmm1
        addps %xmm3,%xmm6
        addps %xmm1,%xmm6
        movaps %xmm6,nb202nf_vctot(%esp)

        ## should we do one more iteration? 
        subl $4,nb202nf_innerk(%esp)
        jl    _nb_kernel202nf_ia32_sse.nb202nf_single_check
        jmp   _nb_kernel202nf_ia32_sse.nb202nf_unroll_loop
_nb_kernel202nf_ia32_sse.nb202nf_single_check: 
        addl $4,nb202nf_innerk(%esp)
        jnz   _nb_kernel202nf_ia32_sse.nb202nf_single_loop
        jmp   _nb_kernel202nf_ia32_sse.nb202nf_updateouterdata
_nb_kernel202nf_ia32_sse.nb202nf_single_loop: 
        movl  nb202nf_innerjjnr(%esp),%edx       ## pointer to jjnr[k] 
        movl  (%edx),%eax
        addl $4,nb202nf_innerjjnr(%esp)

        movl nb202nf_pos(%ebp),%esi
        leal  (%eax,%eax,2),%eax

        ## fetch j coordinates 
        xorps %xmm3,%xmm3
        xorps %xmm4,%xmm4
        xorps %xmm5,%xmm5

        movss (%esi,%eax,4),%xmm3               ## jxO  -  -  -
        movss 4(%esi,%eax,4),%xmm4              ## jyO  -  -  -
        movss 8(%esi,%eax,4),%xmm5              ## jzO  -  -  -  

        movlps 12(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1   -    -
        movss  20(%esi,%eax,4),%xmm7            ## xmm7 = jzH1   -    -    - 
        movhps 24(%esi,%eax,4),%xmm6            ## xmm6 = jxH1 jyH1 jxH2 jyH2
        movss  32(%esi,%eax,4),%xmm2            ## xmm2 = jzH2   -    -    -

        ## have all coords, time for some shuffling.

        shufps $216,%xmm6,%xmm6 ## constant 11011000     ;# xmm6 = jxH1 jxH2 jyH1 jyH2 
        unpcklps %xmm2,%xmm7                    ## xmm7 = jzH1 jzH2   -    -
        movaps  nb202nf_ixO(%esp),%xmm0
        movaps  nb202nf_iyO(%esp),%xmm1
        movaps  nb202nf_izO(%esp),%xmm2
        movlhps %xmm6,%xmm3                     ## xmm3 = jxO   0   jxH1 jxH2 
        shufps $228,%xmm6,%xmm4 ## constant 11100100    ;# xmm4 = jyO   0   jyH1 jyH2 
        shufps $68,%xmm7,%xmm5 ## constant 01000100    ;# xmm5 = jzO   0   jzH1 jzH2

        ## store all j coordinates in jO  
        movaps %xmm3,nb202nf_jxO(%esp)
        movaps %xmm4,nb202nf_jyO(%esp)
        movaps %xmm5,nb202nf_jzO(%esp)
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
        mulps   nb202nf_krf(%esp),%xmm6   ## xmm6=krsq 
        movaps  %xmm1,%xmm2
        mulps   %xmm1,%xmm1
        movaps  nb202nf_three(%esp),%xmm3
        mulps   %xmm0,%xmm1
        subps   %xmm1,%xmm3
        mulps   %xmm2,%xmm3
        mulps   nb202nf_half(%esp),%xmm3   ## rinv iO - j water 

        addps   %xmm3,%xmm6     ## xmm6=rinv+ krsq 
        subps   nb202nf_crf(%esp),%xmm6   ## xmm6=rinv+ krsq-crf 

        xorps   %xmm1,%xmm1
        movaps  %xmm3,%xmm0
        subps   %xmm7,%xmm3     ## xmm3=rinv-2*krsq 
        xorps   %xmm4,%xmm4
        ## fetch charges to xmm4 (temporary) 
        movss   nb202nf_qqOO(%esp),%xmm4
        movhps  nb202nf_qqOH(%esp),%xmm4

        mulps %xmm4,%xmm6       ## vcoul  

        addps   nb202nf_vctot(%esp),%xmm6
        movaps  %xmm6,nb202nf_vctot(%esp)

        ## done with i O Now do i H1 & H2 simultaneously first get i particle coords: 
        movaps  nb202nf_ixH1(%esp),%xmm0
        movaps  nb202nf_iyH1(%esp),%xmm1
        movaps  nb202nf_izH1(%esp),%xmm2
        movaps  nb202nf_ixH2(%esp),%xmm3
        movaps  nb202nf_iyH2(%esp),%xmm4
        movaps  nb202nf_izH2(%esp),%xmm5
        subps   nb202nf_jxO(%esp),%xmm0
        subps   nb202nf_jyO(%esp),%xmm1
        subps   nb202nf_jzO(%esp),%xmm2
        subps   nb202nf_jxO(%esp),%xmm3
        subps   nb202nf_jyO(%esp),%xmm4
        subps   nb202nf_jzO(%esp),%xmm5
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
        movaps  nb202nf_three(%esp),%xmm3
        movaps  %xmm3,%xmm7
        mulps   %xmm0,%xmm1
        mulps   %xmm4,%xmm5
        subps   %xmm1,%xmm3
        subps   %xmm5,%xmm7
        mulps   %xmm2,%xmm3
        mulps   %xmm6,%xmm7
        mulps   nb202nf_half(%esp),%xmm3   ## rinv H1 - j water 
        mulps   nb202nf_half(%esp),%xmm7   ## rinv H2 - j water  

        mulps nb202nf_krf(%esp),%xmm0   ## krsq 
        mulps nb202nf_krf(%esp),%xmm4   ## krsq  

        ## assemble charges in xmm6 
        xorps   %xmm6,%xmm6
        movss   nb202nf_qqOH(%esp),%xmm6
        movhps  nb202nf_qqHH(%esp),%xmm6
        addps   %xmm3,%xmm0     ## krsq+ rinv 
        addps   %xmm7,%xmm4     ## krsq+ rinv 
        subps   nb202nf_crf(%esp),%xmm0
        subps   nb202nf_crf(%esp),%xmm4
        mulps   %xmm6,%xmm0     ## vcoul 
        mulps   %xmm6,%xmm4     ## vcoul 
        addps   %xmm0,%xmm4
        addps   nb202nf_vctot(%esp),%xmm4
        movaps  %xmm4,nb202nf_vctot(%esp)

        decl nb202nf_innerk(%esp)
        jz    _nb_kernel202nf_ia32_sse.nb202nf_updateouterdata
        jmp   _nb_kernel202nf_ia32_sse.nb202nf_single_loop
_nb_kernel202nf_ia32_sse.nb202nf_updateouterdata: 
        ## get n from stack
        movl nb202nf_n(%esp),%esi
        ## get group index for i particle 
        movl  nb202nf_gid(%ebp),%edx            ## base of gid[]
        movl  (%edx,%esi,4),%edx                ## ggid=gid[n]

        ## accumulate total potential energy and update it 
        movaps nb202nf_vctot(%esp),%xmm7
        ## accumulate 
        movhlps %xmm7,%xmm6
        addps  %xmm6,%xmm7      ## pos 0-1 in xmm7 have the sum now 
        movaps %xmm7,%xmm6
        shufps $1,%xmm6,%xmm6
        addss  %xmm6,%xmm7

        ## add earlier value from mem 
        movl  nb202nf_Vc(%ebp),%eax
        addss (%eax,%edx,4),%xmm7
        ## move back to mem 
        movss %xmm7,(%eax,%edx,4)

        ## finish if last 
        movl nb202nf_nn1(%esp),%ecx
        ## esi already loaded with n
        incl %esi
        subl %esi,%ecx
        jz _nb_kernel202nf_ia32_sse.nb202nf_outerend

        ## not last, iterate outer loop once more!  
        movl %esi,nb202nf_n(%esp)
        jmp _nb_kernel202nf_ia32_sse.nb202nf_outer
_nb_kernel202nf_ia32_sse.nb202nf_outerend: 
        ## check if more outer neighborlists remain
        movl  nb202nf_nri(%esp),%ecx
        ## esi already loaded with n above
        subl  %esi,%ecx
        jz _nb_kernel202nf_ia32_sse.nb202nf_end
        ## non-zero, do one more workunit
        jmp   _nb_kernel202nf_ia32_sse.nb202nf_threadloop
_nb_kernel202nf_ia32_sse.nb202nf_end: 
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




